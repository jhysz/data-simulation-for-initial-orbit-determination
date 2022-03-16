sat_list = {
    'envisat': {
        'norad_id': 27386,
        'cospar_id': '0200901',
        'sic_id': '6179',
        'mass': 8000.0, # kg; TODO: compute proper value
        'cross_section': 100.0, # m2; TODO: compute proper value
        'cd': 2.0, # TODO: compute proper value
        'cr': 1.0  # TODO: compute proper value
    }
}

sc_name = 'envisat'

import orekit
vm = orekit.initVM()
# Now we have populated the namespace with the orekit classes

from orekit.pyhelpers import setup_orekit_curdir
setup_orekit_curdir()

from java.util import Arrays
from orekit import JArray_double

from org.orekit.bodies import  OneAxisEllipsoid
from org.orekit.frames import  FramesFactory
from org.orekit.data import DataProvidersManager, ZipJarCrawler
from org.orekit.time import TimeScalesFactory, AbsoluteDate
from org.orekit.orbits import KeplerianOrbit
from org.orekit.utils import Constants
from org.orekit.propagation.analytical import EcksteinHechlerPropagator
from org.orekit.propagation.analytical.tle import TLEPropagator
from org.orekit.propagation.conversion import FiniteDifferencePropagatorConverter
from org.orekit.propagation.conversion import TLEPropagatorBuilder
from datetime import datetime
from org.orekit.propagation import SpacecraftState
from org.orekit.orbits import OrbitType, PositionAngle
from org.orekit.propagation.numerical import NumericalPropagator
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
from org.orekit.utils import IERSConventions
from org.orekit.frames import FramesFactory, TopocentricFrame
from org.orekit.bodies import OneAxisEllipsoid, GeodeticPoint
from org.orekit.forces.gravity import ThirdBodyAttraction
from org.orekit.forces.radiation import IsotropicRadiationSingleCoefficient
from org.orekit.forces.radiation import SolarRadiationPressure
from org.orekit.models.earth.atmosphere import NRLMSISE00
from org.orekit.models.earth.atmosphere.data import MarshallSolarActivityFutureEstimation
# from org.orekit.forces.drag.atmosphere.data import MarshallSolarActivityFutureEstimation
# from org.orekit.forces.drag.atmosphere import NRLMSISE00
from org.orekit.forces.drag import IsotropicDrag
from org.orekit.forces.drag import DragForce
from org.orekit.bodies import CelestialBodyFactory
from org.orekit.models.earth import ReferenceEllipsoid

from java.io import File
from org.orekit.data import DataProvidersManager, DirectoryCrawler
from orekit import JArray

import matplotlib.dates as mdates
from math import radians, pi
import matplotlib.pyplot as plt
from multi_target import multi_simu, show_simu
import numpy as np

# Keplerian
# a = 7445392.44286               # Semimajor axis in m (a)
# e = 0.003113099                 # Eccentricity (e)
# i = radians(5.70529451)         # Inclination (i).
# raan = radians(90.)             # Right ascension of the ascending node (Ω).
# omega = radians(179.93295544)   # Argument of pericenter (ω).
# lv = radians(92.8717018)        # True anomaly (ν).

def orekit_prop_spacebased(a, e, i, raan, omega, lv):
    # Some Constants
    ae = Constants.WGS84_EARTH_EQUATORIAL_RADIUS
    mu = Constants.WGS84_EARTH_MU
    utc = TimeScalesFactory.getUTC()

    # Keplerian
    epochDate = AbsoluteDate(2020, 1, 26, 16, 0, 00.000, utc)
    # Inertial frame where the satellite is defined
    # inertialFrame = FramesFactory.getEME2000()
    inertialFrame = FramesFactory.getGCRF()

    # Orbit construction as Keplerian
    initialOrbit = KeplerianOrbit(a, e, i, omega, raan, lv,
                                  PositionAngle.TRUE,
                                  inertialFrame, epochDate, mu)
    satellite_mass = 100.0 # kg
    initialState = SpacecraftState(initialOrbit, satellite_mass)

    minStep = 0.001
    maxstep = 1000.0
    initStep = 60.0

    positionTolerance = 1.0
    orbitType = OrbitType.CARTESIAN
    tol = NumericalPropagator.tolerances(positionTolerance, initialOrbit, orbitType)

    integrator = DormandPrince853Integrator(minStep, maxstep,
        JArray_double.cast_(tol[0]),  # Double array of doubles needs to be casted in Python
        JArray_double.cast_(tol[1]))
    integrator.setInitialStepSize(initStep)

    propagator_num = NumericalPropagator(integrator)
    propagator_num.setOrbitType(orbitType)
    propagator_num.setInitialState(initialState)

    itrf    = FramesFactory.getITRF(IERSConventions.IERS_2010, True) # International Terrestrial Reference Frame, earth fixed
    earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                             Constants.WGS84_EARTH_FLATTENING,
                             itrf)
    gravityProvider = GravityFieldFactory.getNormalizedProvider(8, 8)
    propagator_num.addForceModel(HolmesFeatherstoneAttractionModel(earth.getBodyFrame(), gravityProvider))

    ecef = itrf
    wgs84Ellipsoid = ReferenceEllipsoid.getWgs84(ecef)


    # Moon and Sun perturbations
    moon = CelestialBodyFactory.getMoon()
    sun = CelestialBodyFactory.getSun()

    moon_3dbodyattraction = ThirdBodyAttraction(moon)
    propagator_num.addForceModel(moon_3dbodyattraction)

    sun_3dbodyattraction = ThirdBodyAttraction(sun)
    propagator_num.addForceModel(sun_3dbodyattraction)

    # Solar radiation pressure
    isotropicRadiationSingleCoeff = IsotropicRadiationSingleCoefficient(sat_list[sc_name]['cross_section'],
                                                                        sat_list[sc_name]['cr'])

    solarRadiationPressure = SolarRadiationPressure(sun, wgs84Ellipsoid.getEquatorialRadius(),
                                                    isotropicRadiationSingleCoeff)
    propagator_num.addForceModel(solarRadiationPressure)

    # Atmospheric drag

    msafe = MarshallSolarActivityFutureEstimation(
        '(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\p{Digit}\p{Digit}\p{Digit}\p{Digit}F10\.(?:txt|TXT)',
        MarshallSolarActivityFutureEstimation.StrengthLevel.AVERAGE)


    orekit_filename = 'orekit-data'
    DM = DataProvidersManager.getInstance()
    datafile = File(orekit_filename)
    if not datafile.exists():
        print('Directory :', datafile.absolutePath, ' not found')

    crawler = DirectoryCrawler(datafile)
    DM.clearProviders()
    DM.addProvider(crawler)
    # DM.feed(msafe.getSupportedNames(), msafe)  # Feeding the F10.7 bulletins to Orekit's data manager

    atmosphere = NRLMSISE00(msafe, sun, wgs84Ellipsoid)
    # # from org.orekit.forces.drag.atmosphere import DTM2000
    # # atmosphere = DTM2000(msafe, sun, wgs84Ellipsoid)

    isotropicDrag = IsotropicDrag(sat_list[sc_name]['cross_section'], sat_list[sc_name]['cd'])

    dragForce = DragForce(atmosphere, isotropicDrag)
    propagator_num.addForceModel(dragForce)


    ## Create time vector
    startDate = AbsoluteDate(2020, 1, 26, 11, 0, 00.000, utc)

    # Overall duration in seconds for extrapolation
    duration = 86400 * 3
    step_time = 60 * 3

    # Time array in orekit AbsoluteDate format
    t = [startDate.shiftedBy(float(dt)) \
            for dt in np.arange(0, duration, step_time)]

    state = [propagator_num.propagate(tt) for tt in t]
    pos = [propagator_num.propagate(tt).getPVCoordinates().getPosition() for tt in t]


    longitude = radians(21.063)
    latitude  = radians(67.878)
    altitude  = 341.0
    station = GeodeticPoint(latitude, longitude, altitude)

    station_frame = TopocentricFrame(earth, station, "Esrange")

    el_tmp = []
    az_tmp = []
    for tt in t :
        position = propagator_num.propagate(tt).getPVCoordinates().getPosition()
        el = station_frame.getElevation(position,inertialFrame,tt) * 180.0 / pi
        az = station_frame.getAzimuth(position,inertialFrame,tt) * 180.0 / pi
        el_tmp.append(el)
        az_tmp.append(az)

    plt.figure()
    plt.plot(el_tmp)
    plt.show()
    plt.figure()
    plt.plot(az_tmp)
    plt.show()

    print(pos[0:9])
    plt.figure()
    plt.plot([n.x for n in pos])
    plt.show()

    pos_x = [n.x for n in pos]
    pos_y = [n.y for n in pos]
    pos_z = [n.z for n in pos]

    np.save(file='pos_x',arr=pos_x)
    np.save(file='pos_y',arr=pos_y)
    np.save(file='pos_z',arr=pos_z)

    pos_x = np.load(file='pos_x.npy')
    pos_y = np.load(file='pos_y.npy')
    pos_z = np.load(file='pos_z.npy')

    obs_target = multi_simu(pos_x,pos_y,pos_z)
    show_simu(obs_target)