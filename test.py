# from multi_target import multi_simu,show_simu
# import numpy as np
#
# pos_x = np.load(file='pos_x.npy')
# pos_y = np.load(file='pos_y.npy')
# pos_z = np.load(file='pos_z.npy')
#
# obs_target = multi_simu(pos_x,pos_y,pos_z)
# show_simu(obs_target)

from space_orekit_prop import orekit_prop_spacebased
from ground_orekit_prop import orekit_prop_groundbased
from math import radians

# a = 7445392.44286               # Semimajor axis in m (a)
# e = 0.003113099                 # Eccentricity (e)
# i = radians(5.70529451)         # Inclination (i).
# raan = radians(90.)             # Right ascension of the ascending node (Ω).
# omega = radians(179.93295544)   # Argument of pericenter (ω).
# lv = radians(92.8717018)        # True anomaly (ν).

a = 12445392.44286               # Semimajor axis in m (a)
e = 0.001113099                 # Eccentricity (e)
i = radians(25.70529451)         # Inclination (i).
raan = radians(60.)             # Right ascension of the ascending node (Ω).
omega = radians(80.93295544)   # Argument of pericenter (ω).
lv = radians(92.8717018)        # True anomaly (ν).

print('space_based:')

orekit_prop_spacebased(a,e,i,raan,omega,lv)

print('ground_based:')

orekit_prop_groundbased(a,e,i,raan,omega,lv)