import empymod
import numpy as np
import xarray 
import matplotlib.pyplot as plt
             
#-------------------------------------------------
iTx, iRx, ichrec, ifreq, dre, dim=np.loadtxt('emf_0001.txt', skiprows=1, unpack=True)
dat = dre + dim*1j

xs, ys, zs, azimuth, dip, isrc=np.loadtxt('sources.txt', skiprows=1, unpack=True) 
xr, yr, zr, azimuth, dip, irec=np.loadtxt('receivers.txt', skiprows=1, unpack=True) 
nrec = int(max(irec))/2
xr -= xs

data = xarray.load_dataset('./marlim_data.nc', engine='h5netcdf')
#print(data.src_x[0::2])==xr
#print(data.src_x[1::2])==xr

dil = np.abs(getattr(data, 'data_il').data[::2, :, :] + 1j*getattr(data, 'data_il').data[1::2, :, :])
dbs = np.abs(getattr(data, 'data_bs').data[::2, :, :] + 1j*getattr(data, 'data_bs').data[1::2, :, :])

plt.figure(figsize=(12,9))
plt.subplot(221)
idx = (ifreq==1) & (ichrec==1) & (iRx<=nrec)
plt.plot(xr[irec<=nrec], np.abs(dat[idx]), 'r--', markevery=15, marker='o', label='$E_x^{FD}$-0.125 Hz')
plt.plot(xr[irec<=nrec], dil[:, 5, 0], 'b-.', label='$E_x^{ref}$-1.25 Hz')

plt.axhline(2e-15, c='k')# Plot noise level.
plt.grid()
plt.yscale('log')
plt.ylabel('|$E_x$| (V/Am)')
plt.xlabel('X (m)')
plt.title('(a) Inline amplitude', fontweight='bold')

plt.subplot(222)

idx = (ifreq==1) & (ichrec==1) & (iRx>nrec)
plt.plot(xr[irec>nrec], np.abs(dat[idx]), 'r--', markevery=15, marker='o', label='$E_x^{FD}$-0.125 Hz')
plt.plot(xr[irec>nrec], dbs[:, 5, 0], 'b-.', label='$E_x^{ref}$-1.25 Hz')


plt.axhline(2e-15, c='k')# Plot noise level.
plt.grid()
plt.yscale('log')
plt.ylabel('|$E_x$| (V/Am)')
plt.xlabel('X (m)')
plt.title('(b) Broadside amplitude', fontweight='bold')
#plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0.)
plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1.02))


plt.subplot(223)
idx = (ifreq==1) & (ichrec==1) & (iRx<=nrec)
err = np.abs(dat[idx]/dil[:, 5, 0])-1
plt.plot(xr[irec<=nrec], 100*err, 'b-', markevery=15, marker='^', label='$E_x$-1.25 Hz')

plt.ylim([-15,15])
plt.grid()
plt.ylabel('|$E_x^{FD}/E_x^{ref}$|-1 (%)')
plt.xlabel('X (m)')
plt.title('(c) Inline error', fontweight='bold')


plt.subplot(224)
idx = (ifreq==1) & (ichrec==1) & (iRx>nrec)
err = np.abs(dat[idx]/dbs[:, 5, 0])-1
plt.plot(xr[irec>nrec], 100*err, 'b-', markevery=15, marker='^', label='$E_x$-1.25 Hz')

plt.ylim([-15,15])
plt.grid()
plt.ylabel('|$E_x^{FD}/E_x^{ref}$|-1 (%)')
plt.xlabel('X (m)')
plt.title('(d) Broadside error', fontweight='bold')
plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1))


plt.tight_layout()
plt.savefig('marlim_ex.png')
plt.show()




