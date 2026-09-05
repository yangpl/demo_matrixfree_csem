import emg3d
import numpy as np
import xarray as xr

#=============================================
#load model data: coordinates and resistivities
#=============================================
data = np.load('./marlim_comp.npz')  #  or 'marlim_orig.npz'
hx = data['hx']
hy = data['hy']
hz = data['hz']
x0, y0, z0 = data['x0'] #origin of the model
x = x0 + np.concatenate(([0], np.cumsum(hx)))
y = y0 + np.concatenate(([0], np.cumsum(hy)))
z = z0 + np.concatenate(([0], np.cumsum(hz)))
z = -z[::-1]
x1 = x.astype(np.float32)
x2 = y.astype(np.float32)
x3 = z.astype(np.float32)
x1.tofile("x1");
x2.tofile("x2");
x3.tofile("x3");
print(x1[0],np.sum(hx))
print(x2[0],np.sum(hy))
print(x3[0],np.sum(hz))

res_h = data['res_h']
res_v = data['res_v']
Rh = res_h.astype(np.float32)
Rv = res_v.astype(np.float32)
#reverse z axis
Rh = Rh[:,:,::-1]
Rv = Rv[:,:,::-1]
Rh = Rh.transpose((2,1,0))
Rv = Rv.transpose((2,1,0))
Rh.tofile("Rh")
Rv.tofile("Rv")


