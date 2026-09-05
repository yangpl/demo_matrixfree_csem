import numpy as np
import matplotlib.pyplot as plt

x1 = np.fromfile('fx1', dtype=np.float32)
x2 = np.fromfile('fx2', dtype=np.float32)
x3 = np.fromfile('fx3', dtype=np.float32)
#x3 = -x3[::-1]

#topo_x, topo_z= np.loadtxt('ftopo.txt', skiprows=0, unpack=True)
xr, yr, zr, azimuth, dip, irec = np.loadtxt('receivers.txt', skiprows=1, unpack=True);
xs, ys, zs, azimuth, dip, isrc = np.loadtxt('sources.txt', skiprows=1, unpack=True);


plt.figure(figsize=(8,6))

plt.subplot(211)
[xx, yy] = np.meshgrid(x1, x2);
plt.plot(xx, yy, 'k', linewidth=0.5)
plt.plot(xx.T, yy.T, 'k', linewidth=0.5)
#plt.plot(topo_x, topo_z, 'r-')
plt.plot(xr, yr, 'b.', markersize=5)
plt.plot(xs, ys, 'gd', markersize=5)
plt.text(xs+16e3, ys+200, 'Inline', fontweight='bold', ha='right', va='top')
plt.text(xs+18e3, ys-800, 'Broadside', fontweight='bold', ha='right', va='top')

plt.xlim(x1[0], x1[-1])
plt.ylim(7.515e6, 7.52e6)
#plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.title('(a)', fontweight='bold')

plt.subplot(212)
[xx, zz] = np.meshgrid(x1, x3);
plt.plot(xx, zz, 'k', linewidth=0.5)
plt.plot(xx.T, zz.T, 'k', linewidth=0.5)
plt.plot(xr, zr, 'b.', markersize=5)
plt.plot(xs, zs, 'gd', markersize=5)

plt.xlim(x1[0], x1[-1])
plt.ylim(x3[-1], x3[0])
plt.xlabel('X (m)')
plt.ylabel('Z (m)')
plt.title('(b)', fontweight='bold')

plt.tight_layout()
plt.savefig('mr3d_fdmesh.png', dpi=600)
plt.show()
