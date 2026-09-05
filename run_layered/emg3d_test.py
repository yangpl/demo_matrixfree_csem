import emg3d
import empymod
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

#-------------------------------------------------
x, y, z, azimuth, dip, irec=np.loadtxt('receivers.txt', skiprows=1, unpack=True) 
iTx, iRx, ichrec, ifreq, dre, dim=np.loadtxt('emf_0001.txt', skiprows=1, unpack=True)
dat = dre - dim*1j
ref1d = empymod.dipole(ab=11,src=[0, 0, 950], #ab=11,Ex
                       rec=[x, x*0, 1000],
                       depth=[0, 1000, 2300, 2400],
		       freqtime=[1.25,0.25,0.75,1.25],
                       res=[1e8, 0.3, 1.0, 50, 2],
                       aniso=[1,1,1,1.414,1],
                       verb=1)

freq = 1.25 # Frequency (Hz)
x = np.fromfile('fx1_extended', dtype=np.float64)
y = np.fromfile('fx2_extended', dtype=np.float64)
z = np.fromfile('fx3_extended', dtype=np.float64)
z = -z[::-1]

hx = np.diff(x)
hy = np.diff(y)
hz = np.diff(z)
nx = hx.size
ny = hy.size
nz = hz.size

print(nx)

rho11 = np.fromfile('frho11_extended', dtype=np.float64)
rho22 = np.fromfile('frho22_extended', dtype=np.float64)
rho33 = np.fromfile('frho33_extended', dtype=np.float64)

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
#grid.plot_3d_slicer(
#    model.property_x, zslice=-1900,
#    xlim=(-10e3, 10e3), ylim=(-10e3, 10e3), zlim=(-5e3, 0),
#    pcolor_opts={'norm': LogNorm(vmin=0.1, vmax=100),'cmap':'jet'}
#)

x, y, z, azimuth, dip, iTx = np.loadtxt('sources.txt', skiprows=1, unpack=True)
src = (x, y, -z, azimuth, dip)
x, y, z, azimuth, dip, iRx = np.loadtxt('receivers.txt', skiprows=1, unpack=True)
rec = (x, y, -z, azimuth, dip)

efield = emg3d.solve_source(model, src, freq, tol=1e-8/freq, verb=4)
hfield = emg3d.get_magnetic_field(model, efield)
Ex_rec = efield.get_receiver(rec)
Hy_rec = hfield.get_receiver(rec)


f = open("emf_freq.txt",'w')
for off, re, im in zip(x, Ex_rec.real, Ex_rec.imag):
    f.write('%e \t %e \t %e\n'%(off, re, im))
f.close()

idx = (ifreq==1) & (ichrec==1)

plt.figure(figsize=(10,7))
plt.subplot(221)
plt.plot(x/1e3, np.abs(ref1d[0,:]), 'r', label='$E_x^{ref}$-0.25 Hz')
plt.plot(x/1e3, abs(Ex_rec), 'r--', label='$E_x^{GMG}$-0.25 Hz')
plt.plot(x/1e3, np.abs(dat[idx]), 'b--', label='$E_x^{FDFD}$-0.25 Hz')
plt.grid()
plt.yscale('log')
plt.xlabel('X (km)')
plt.ylabel('$|E_x|$ (V/m)')
plt.title('(a) Amplitude')

plt.subplot(222)
plt.plot(x/1e3, np.angle(ref1d[0,:], deg=True), 'r', label='$E_x^{ref}$-0.25 Hz')
plt.plot(x/1e3, np.angle(Ex_rec, deg=True), 'r--', label='$E_x^{GMG}-0.25$')
plt.plot(x/1e3, np.angle(dat[idx], deg=True), 'b--', label='$E_x^{FD}$-0.25 Hz')
plt.grid()
plt.xlabel('X (km)')
plt.ylabel('$angle(E_x)$ (degree)')
plt.title('(b) Phase')

plt.subplot(223)
err = np.abs(Ex_rec/ref1d[0,:])-1
plt.plot(x/1e3, err, 'r', label='$E_x^{GMG}$-0.25 Hz')
err = np.abs(dat[idx]/ref1d[0,:])-1
plt.plot(x/1e3, err, 'b', label='$E_x^{FDFD}$-0.25 Hz')
plt.grid()
plt.ylim([-0.1,0.1])
plt.xlabel('X (km)')
plt.title('(c) Amplitude error')

plt.subplot(224)
plt.plot(x/1e3, np.angle(Ex_rec/ref1d[0,:], deg=True), 'r', label='$E_x^{GMG}$-0.25 Hz')
plt.plot(x/1e3, np.angle(dat[idx]/ref1d[0,:], deg=True), 'b', label='$E_x^{FDFD}$-0.25 Hz')
plt.grid()
plt.ylim([-5,5])
plt.xlabel('X (km)')
plt.title('(d) Phase error')

plt.tight_layout()
plt.savefig('result_emg3d.png')
plt.show()




