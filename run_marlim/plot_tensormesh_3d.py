from discretize import TensorMesh
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits import mplot3d

x1 = np.fromfile('fx1', dtype=np.float32)
x2 = np.fromfile('fx2', dtype=np.float32)
x3 = np.fromfile('fx3', dtype=np.float32)
x3 = -x3[::-1]
hx = np.diff(x1)
hy = np.diff(x2)
hz = np.diff(x3)
nx = np.size(x1)-1
ny = np.size(x2)-1
nz = np.size(x3)-1


mesh = TensorMesh([hx, hy, hz], origin=(x1[0], x2[0], x3[0]))
print(mesh)

rho = np.fromfile('frho11', dtype=np.float32)
rho = rho.reshape([nx, ny, nz], order='F')
rho = rho[:,:,::-1] #reverse the order of z-axis
mesh.plot_3d_slicer(rho,  xslice=x1[int(nx/2)], yslice=x2[int(ny/2)], zslice=-2750,
                xlim=(364325, 415825), ylim=(7490049, 7546349), zlim=(-6200, 0),
                pcolor_opts={'norm': LogNorm(vmin=0.1,vmax=1500), 'cmap':'jet'})

rho = np.fromfile('frho33', dtype=np.float32)
rho = rho.reshape([nx, ny, nz], order='F')
rho = rho[:,:,::-1] #reverse the order of z-axis
mesh.plot_3d_slicer(rho,  xslice=x1[int(nx/2)], yslice=x2[int(ny/2)], zslice=-2750,
                xlim=(364325, 415825), ylim=(7490049, 7546349), zlim=(-6200, 0),
                pcolor_opts={'norm': LogNorm(vmin=0.1,vmax=1500), 'cmap':'jet'})

plt.show()
