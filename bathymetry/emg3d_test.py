import emg3d
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


freq = 0.25 # Frequency (Hz)
x = np.fromfile('fx1_extended', dtype=np.float32)
y = np.fromfile('fx2_extended', dtype=np.float32)
z = np.fromfile('fx3_extended', dtype=np.float32)
z = -z[::-1]

hx = np.diff(x)
hy = np.diff(y)
hz = np.diff(z)
nx = hx.size
ny = hy.size
nz = hz.size

rho11 = np.fromfile('frho11_extended', dtype=np.float32)
rho22 = np.fromfile('frho11_extended', dtype=np.float32)
rho33 = np.fromfile('frho11_extended', dtype=np.float32)

rho11_3d = rho11.reshape((nz, ny, nx))
rho22_3d = rho22.reshape((nz, ny, nx))
rho33_3d = rho33.reshape((nz, ny, nx))

# Reverse the z-axis (first dimension) to go from bottom-to-top
rho11_reversed_3d = rho11_3d[::-1, :, :]  # or np.flip(rho11_3d, axis=0)
rho22_reversed_3d = rho22_3d[::-1, :, :]  # or np.flip(rho11_3d, axis=0)
rho33_reversed_3d = rho33_3d[::-1, :, :]  # or np.flip(rho11_3d, axis=0)

# Flatten back to vector (maintains x-y-z order)
rho11_ = rho11_reversed_3d.flatten()
rho22_ = rho22_reversed_3d.flatten()
rho33_ = rho33_reversed_3d.flatten()

# We create a mesh with the input values
grid = emg3d.TensorMesh([hx, hy, hz], origin=(x[0], y[0], z[0]))
print(grid)

model = emg3d.Model(grid, rho11_, rho22_, rho33_)
grid.plot_3d_slicer(
    model.property_x, zslice=-1900,
    xlim=(-10e3, 10e3), ylim=(-10e3, 10e3), zlim=(-5e3, 0),
    pcolor_opts={'norm': LogNorm(vmin=0.1, vmax=100),'cmap':'jet'}
)

x, y, z, azimuth, dip, iTx = np.loadtxt('sources.txt', skiprows=1, unpack=True)
src = (x, y, -z, azimuth, dip)
x, y, z, azimuth, dip, iRx = np.loadtxt('receivers.txt', skiprows=1, unpack=True)
rec = (x, y, -z, azimuth, dip)

efield = emg3d.solve_source(model, src, freq, verb=4)
hfield = emg3d.get_magnetic_field(model, efield)
Ex_rec = efield.get_receiver(rec)
Hy_rec = hfield.get_receiver(rec)

f = open("emf_freq.txt",'w')
for off, re, im in zip(x, Ex_rec.real, Ex_rec.imag):
    f.write('%e \t %e \t %e\n'%(off, re, im))
f.close()

plt.figure()
plt.subplot(211)
plt.plot(x/1e3, abs(Ex_rec))

plt.yscale('log')
plt.xlabel('X (km)')
plt.ylabel('$|E_x|$ (V/m)')

plt.subplot(212)
plt.plot(x/1e3, np.angle(Ex_rec, deg=True))
plt.xlabel('X (km)')
plt.ylabel('$angle(E_x)$ (degree)')

plt.tight_layout()
plt.savefig('result_emg3d.png')
plt.show()


