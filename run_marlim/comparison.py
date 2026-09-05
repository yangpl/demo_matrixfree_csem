import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

def demod_line(d, n):
    """Get inlined/broadside and demod ([real; real] to complex)."""
    return abs(getattr(d, n).data[::2, :, :] + 1j*getattr(d, n).data[1::2, :, :])

data = xr.load_dataset('marlim_data.nc', engine='h5netcdf')

# Get offsets for plotting
offs = data.src_x[::2] - data.rec_x
offs /= 1e3

# Line styles for published data.
general = {'mec': 'k', 'mew': 0.2, 'lw': 1, 'markevery': (20, 30),
           'path_effects': [pe.Stroke(linewidth=2, foreground='k'), pe.Normal()]}
styles = {
    1.25: {'marker': '^', 'ms': 5, **general},
    1.0: {'marker': 'v', 'ms': 5, **general},
    0.75: {'marker': 'D', 'ms': 4, **general},
    0.5: {'marker': 'p', 'ms': 5, **general},
    0.25: {'marker': 's', 'ms': 4, **general},
    0.125: {'marker': 'o', 'ms': 4, **general},
}

# Initiate figure.
fig, axs = plt.subplots(3, 2, figsize=(9, 10), sharex=True, sharey=True)

# Loop over Inline/Broadside.
for iii, datname in enumerate(['data_il', 'data_bs']):

    # Get absolute values of this line.
    tdat = demod_line(data, datname)

    # Loop over components Ex, Ey, Ez.
    for ii, comp in enumerate(data.components.values[:3]):

        # Get current axis.
        ax = axs[ii, iii]

        # # Loop over frequencies for our codes.
        # for i, freq in enumerate(data.freqs.values):

        #     # Loop over our codes.
        #     for ic, dat in enumerate([tegd, tcst, tptg, tspg]):
        #         ax.plot(offs[101::-1], dat[101::-1, i, ii], '.5')
        #         ax.plot(offs[102:], dat[102:, i, ii], '.5')

        # Loop over frequencies for published responses.
        for i, freq in enumerate(data.freqs.values):           
            ax.plot(offs[101::-1], tdat[101::-1, i, ii], f"C{i}", **styles[freq], label=f"{freq}")
            ax.plot(offs[102:], tdat[102:, i, ii], f"C{i}", **styles[freq])

        # Plot noise level.
        ax.axhline(2e-15, c='k')
        
        # Grid lines and yscale.
        ax.grid(axis='y', c='0.9')
        ax.set_yscale('log')
        
# Switch off spines and move ticks.
for ax in axs.ravel():
    ax.spines['top'].set_visible(False)
for ax in axs[:, 0].ravel():
    ax.spines['right'].set_visible(False)
for ax in axs[:, 1].ravel():
    ax.yaxis.set_ticks_position('right')
    ax.yaxis.set_label_position('right')
    ax.spines['left'].set_visible(False)

# Titles and labels.
axs[0, 0].set_title("Inline")
axs[0, 1].set_title("Broadside")
axs[2, 0].set_xlabel('Offset (km)')
axs[2, 1].set_xlabel('Offset (km)')
axs[0, 0].set_ylabel('$|E_x|$ (V/m)')
axs[1, 0].set_ylabel('$|E_y|$ (V/m)')
axs[2, 0].set_ylabel('$|E_z|$ (V/m)')

# Limits.
axs[0, 0].set_xlim([offs[0], offs[-1]])

# Annotate note-worthy points with numbers.
bbox = {'fontsize': 14, 'bbox': {"boxstyle" : "circle", 'ec': 'C3', 'fc': 'w'}}
axs[1, 1].annotate("1", (8, 1e-13), c='C3', **bbox)
axs[2, 1].annotate("1", (8, 1e-13), c='C3', **bbox)
axs[2, 0].annotate("2", (-8, 1e-19), c='C3', **bbox)
axs[2, 0].annotate("2", (8, 1e-19), c='C3', **bbox)
axs[2, 1].annotate("2", (-8, 1e-19), c='C3', **bbox)
axs[2, 1].annotate("2", (8, 1e-19), c='C3', **bbox)
axs[1, 0].annotate("3", (0, 1e-16), c='C3', **bbox)

# Tight layout.
fig.tight_layout()

# Add frequency legend.
axs[1, 1].legend(
    title='Frequency (Hz)', bbox_to_anchor=(0, 1.2),
    loc='upper center', borderaxespad=0., ncol=3, framealpha=1)

# Save and show.
fig.savefig(f'results-marlim-responses.png', bbox_inches='tight', dpi=300)
fig.show()
