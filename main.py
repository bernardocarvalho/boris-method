# -*- coding: utf-8 -*-
###############################################################################
#
#                           The Boris method
#
#                      L. Fuster & G. Bogopolsky
#
###############################################################################

# import matplotlib as mpl
# from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
# import constants as cst
from particle import Particle
# plt.rcParams["figure.figsize"] = (9, 8)
plt.rcParams['font.size'] = 14
print('proton cyclotron rotation test case')
N = 1       # Number of particles

# Parameters and fields
charge = constants.e
# charge = cst.elemCharge
mass = constants.m_p #  cst.Mp
print(f"Mp {mass:g} {constants.m_p}")
E0 = np.array((0, 0, 0))
B0 = np.array((0, 0, 1))
w0 = np.abs(charge) * np.sqrt(np.sum(B0*B0, axis=0)) / mass
dt = 0.001 / w0       # timestep
Np = 10             # Number of cyclotronic periods

Tf = Np * 2 * np.pi / w0
Nt = int(Tf // dt)  # Number of timesteps
t = np.arange(0, Nt)*dt
x, y, z = np.zeros((Nt)), np.zeros((Nt)), np.zeros((Nt))    # positions taken by the particle along
vx, vy, vz = np.zeros((Nt)), np.zeros((Nt)), np.zeros((Nt))   # velocities

part = Particle(mass, charge)
part.initPos(0, 0, 0)
part.initSpeed(200, 0, 0)

for i in range(Nt):
    part.push(dt, E0, B0)
    x[i], y[i], z[i] = part.r
    vx[i], vy[i], vz[i] = part.v

# Diagnostics
E = 0.5 * mass * (vx**2 + vy**2 + vz**2)
vx_th = 200 * np.cos(w0*t)
vx_error = np.abs(vx - vx_th)
print(vx_error.max() / vx.max())
vy_th = 200 * np.sin(w0*t + np.pi)
vy_error = np.abs(vy - vy_th)
print(f"vy_error_max / vymax {vy_error.max() / vy.max():2g}")

# Outputs
#   Speed along x
# fig = plt.figure()
fig, axs = plt.subplots(2, 2, figsize=(9, 9), tight_layout=True)

axs[0, 0].set_title('Speed along x')
axs[0, 0].plot(t, vx, label='numeric')
axs[0, 0].plot(t, vx_th, label='analytic')
axs[0, 0].plot(t, vx_error, label='erreur absolue')
axs[0, 0].set_xlabel('$t$ [s]')
axs[0, 0].set_ylabel('$v_x$ [m/s]')
axs[0, 0].legend()
axs[0, 0].set_xlim((t[0], t[-1]))
# axs[0, 0].tight_layout()
# plt.savefig('test_vx.png')

axs[0, 1].set_title('Speed along y')
axs[0, 1].plot(t, vy, label='numeric')
axs[0, 1].plot(t, vy_th, label='analytic')
axs[0, 1].plot(t, vy_error, label='erreur absolue')
axs[0, 1].set_xlabel('$t$ [s]')
axs[0, 1].set_ylabel('$v_y$ [m/s]')
axs[0, 1].legend()
axs[0, 1].set_xlim((t[0], t[-1]))
# axs[0, 1].tight_layout()
# plt.savefig('test_vy.png')
#plt.show()

# 3D trajectory
# fig2 = plt.figure()
# axs[1, 0] = fig2.add_subplot(projection='3d')
# ax = fig.gca(projection='3d')
axs[1, 0].remove()
axs[1, 0] = fig.add_subplot(2,2,3,projection='3d')
axs[0, 1].set_title('3D trajectory')
axs[1, 0].plot(x*1e6, y*1e6, z)
axs[1, 0].set_xlabel('$x$ [µm]')
axs[1, 0].set_ylabel('$y$ [µm]')
axs[1, 0].set_zlabel('$z$ [m]')
# plt.savefig('test_3d.png')

# # Energy vs. time
axs[1, 1].plot(t, E)
axs[1, 1].set_xlabel('$t$ [s]')
axs[1, 1].set_ylabel('$E_c$ [J]')
axs[1, 1].set_xlim((t[0], t[-1]))
# plt.tight_layout()
# plt.savefig('test_E.png')
plt.show()
