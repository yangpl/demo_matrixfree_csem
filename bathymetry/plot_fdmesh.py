import numpy as np
import matplotlib.pyplot as plt

x1 = np.fromfile('fx1', dtype=np.float32)
x2 = np.fromfile('fx2', dtype=np.float32)
x3 = np.fromfile('fx3', dtype=np.float32)
#x3 = -x3[::-1]

[xx, zz] = np.meshgrid(x1, x3);

topo_x, topo_z= np.loadtxt('ftopo.txt', skiprows=0, unpack=True)
xr, yr, zr, azimuth, dip, irec = np.loadtxt('receivers.txt', skiprows=1, unpack=True);
xs, ys, zs, azimuth, dip, isrc = np.loadtxt('sources.txt', skiprows=1, unpack=True);


plt.figure()
plt.plot(xx, zz, 'k', linewidth=0.5)
plt.plot(xx.T, zz.T, 'k', linewidth=0.5)
plt.plot(topo_x, topo_z, 'r-')
plt.plot(xr, zr, 'bd', markersize=3)
plt.plot(xs, zs, 'g.')

plt.xlim(x1[0], x1[-1])
plt.ylim(x3[-1], x3[0])
plt.xlabel('X (m)')
plt.ylabel('Z (m)')

plt.savefig('fdmesh.png', dpi=600)
plt.show()
