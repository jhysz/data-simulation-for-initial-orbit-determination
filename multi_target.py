import numpy as np
from joblib import Parallel, delayed
from sklearn.base import BaseEstimator
# Astrodynamics libraries
from astropy import units as u
from astropy.time import Time
from poliastro.bodies import Earth, Sun
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import cowell
from poliastro.core.perturbations import J2_perturbation
# from numba import njit
from poliastro.core.propagation import func_twobody

from astropy import units as u
from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit
import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

GM = 398600.4415 * 1e9  # JGM-3, Nominal Earth mass parameter.
R = 6378.137  # WGS-84, Nominal Earth equatorial radius.

pos_x = np.load(file='pos_x.npy')
pos_y = np.load(file='pos_y.npy')
pos_z = np.load(file='pos_z.npy')

def multi_simu(pos_x,pos_y,pos_z):

    rd = random.random()
    rand_num = np.round(rd * len(pos_x))
    rand_num = np.clip(rand_num, 0, len(pos_x)-1)
    obse_pos = np.array([pos_x[int(rand_num)], pos_y[int(rand_num)], pos_z[int(rand_num)]])

    n = 1
    obs_target = []

    while True:

        delta_d = np.array([random.gauss(0, 2000), random.gauss(0, 2000), random.gauss(0, 2000)])
        obse_pos = obse_pos + delta_d

        v = np.sqrt(GM / np.linalg.norm(obse_pos))  # The speed of the circular orbit
        v_squ = v ** 2

        r1 = random.random()
        r2 = random.random()
        r3 = random.random()
        rr = r1 + r2 + r3

        vx = random.choice((-1, 1)) * (np.sqrt(v_squ * r1 / rr) + random.uniform(-1000, 1000))
        vy = random.choice((-1, 1)) * (np.sqrt(v_squ * r2 / rr) + random.uniform(-1000, 1000))
        vz = random.choice((-1, 1)) * (np.sqrt(v_squ * r3 / rr) + random.uniform(-1000, 1000))

        r_ = obse_pos * 1e-3 * u.km
        v_ = np.array([vx, vy, vz]) * 1e-3 * u.km / u.s

        orb = Orbit.from_vectors(Earth, r_, v_)

        a_ = orb.classical()[0].to(u.m).to_value()
        e_ = orb.classical()[1].to_value()

        orb.plot(interactive=True, use_3d=True)

        if a_ - a_ * e_ > 6480000:

            a_s = orb.classical()[0].to(u.m).to_value()
            e_s = orb.classical()[1].to_value()
            i_s = orb.classical()[2].to(u.deg).to_value()
            ra_s = orb.classical()[3].to(u.deg).to_value()
            om_s = orb.classical()[4].to(u.deg).to_value()
            lv_s = orb.classical()[5].to(u.deg).to_value()

            obs_target.append([a_s, e_s, i_s, ra_s, om_s, lv_s])

            n += 1

            if len(obs_target) >= 10:
                print('the number of simulation :', len(obs_target))
                break

    np.save(file = 'obs_simulation.npy',arr=obs_target)
    return obs_target



def f(t0, u_, k):
    du_kep = func_twobody(t0, u_, k)
    ax, ay, az = J2_perturbation(
        t0, u_, k, J2=Earth.J2.value, R=Earth.R.to(u.km).value
    )
    du_ad = np.array([0, 0, 0, ax, ay, az])
    return du_kep + du_ad

def prop(state):
    s = 180
    epoch = Time("2011-08-05 16:25", scale="utc")
    r = state[0:3] * u.m
    v = state[3:6] * (u.m / u.s)
    orbit = Orbit.from_vectors(Earth, r, v, epoch=epoch)
    prop_orbit = orbit.propagate(s * u.s,
                                 method=cowell, f=f)

    prop_r, prop_v = prop_orbit.rv()
    prop_state_vect = np.concatenate([prop_r.to(u.m).to_value(), prop_v.to(u.m / u.s).to_value()])
    return prop_state_vect



def show_simu(obs_target):

    length = len(obs_target)
    num_prop = 300
    sum_ss = np.zeros((length, num_prop, 6))

    for i in range(length):
        state = obs_target[i]
        a = state[0] * u.m
        ecc = state[1] * u.one
        inc = state[2] * u.deg
        raan = state[3] * u.deg
        argp = state[4] * u.deg
        nu = state[5] * u.deg
        orbit = Orbit.from_classical(Earth, a, ecc, inc, raan, argp, nu)
        orbit.rv()
        orb_r, orb_v = orbit.rv()
        orb_state = np.concatenate([orb_r.to(u.m).to_value(), orb_v.to(u.m / u.s).to_value()])

        state_series = []
        state = orb_state

        for j in range(0, num_prop):
            state = prop(state)
            state_series.append(state)

        state_series = np.array(state_series)
        sum_ss[i] = state_series

    plt.figure()
    plt.plot(pos_y,pos_z,color = 'blue', linewidth = 5,label='space-based telescope')
    for i in range(length):
        plt.plot(sum_ss[i,:,1],sum_ss[i,:,2])
    plt.xlabel('y / m')
    plt.ylabel('z / m')
    plt.savefig('yz_plane.png')
    plt.legend()
    plt.show()

    # 定义坐标轴
    fig = plt.figure()
    ax1 = plt.axes(projection='3d')
    # ax = fig.add_subplot(111,projection='3d')  #这种方法也可以画多个子图
    ax1.plot3D(pos_x,pos_y,pos_z, color='blue', linewidth=5, label='space-based telescope')
    for i in range(length):
        ax1.plot3D(sum_ss[i, :, 0], sum_ss[i, :, 1], sum_ss[i, :, 2])  # 绘制空间曲线
    plt.legend()
    plt.savefig('3d_simulation.png')
    plt.show()