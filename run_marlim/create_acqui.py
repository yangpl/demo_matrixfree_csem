import numpy as np
import xarray as xr

#=============================================
#load survey data: src and rec coordinates
#=============================================
# Load survey as template
ds = xr.load_dataset('./marlim_survey.nc', engine='h5netcdf')

# Use reciprocity: rec becomes src
src = np.column_stack((ds.rec_x, ds.rec_y, -ds.rec_z, ds.rec_theta, ds.rec_dip, 1))
print(src)
header_text = "x   y   z  azimulth   dip   isrc"
fmt = ["%.6f", "%.6f", "%.6f", "%.6f", "%.6f", "%d"]
np.savetxt("sources.txt", src, fmt=fmt, header=header_text, comments='')

# Use reciprocity: src becomes rec, il=inline; bs=broadside
rec_x = ds.data_il.src_x[::2]
rec_y_il = ds.data_il.src_y
rec_z_il = ds.data_il.src_z

# Ensure same coordinates
#print(np.allclose(rec_x, ds.data_bs.src_x[::2]))

rec_y_bs = ds.data_bs.src_y
rec_z_bs = ds.data_bs.src_z


f = open("receivers.txt",'w')
f.write( "x   y   z  azimulth   dip   irec\n")
ff = open("src_rec_table.txt", "w")
ff.write( "isrc   irec\n")
irec = 0
for x, z in zip(rec_x, rec_z_il):
    y = rec_y_il
    irec = irec + 1
    f.write('%e \t %e \t %e \t 0 \t 0 \t %d\n'%(x, y, -z, irec))
    ff.write('%d \t %d\n'%(1, irec))
for x, z in zip(rec_x, rec_z_bs):
    y = rec_y_bs
    irec = irec + 1
    f.write('%e \t %e \t %e \t 0 \t 0 \t %d\n'%(x, y, -z, irec))
    ff.write('%d \t %d\n'%(1, irec))
f.close()
ff.close()

# Frequency
freqs = ds.freqs.values
print(freqs)

