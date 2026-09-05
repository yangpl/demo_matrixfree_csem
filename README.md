# Matrix-free CSEM demo

This repository contains a parallel 3-D controlled-source electromagnetic
(CSEM) modelling code. The default executable applies the frequency-domain
Maxwell operator without assembling the system matrix and solves it with
BiCGStab or GMRES. MPI is used for distributed vector operations and OpenMP
is enabled in the C build.

The source tree also contains matrix-storage experiments in `Acoo.c` and
`Acsr.c`; they are compiled into the executable, but the active entry point is
the matrix-free `main()` in `main.c`.

## Requirements

The default `Makefile` expects:

- MPI C compiler and launcher: `mpicc`, `mpirun`
- GCC/Fortran runtime and OpenMP
- FFTW3
- MUMPS, PORD, ParMETIS/METIS, Scotch, ScaLAPACK, LAPACK, and BLAS
- Python 3 with NumPy and Matplotlib for plotting
- `emg3d` and `empymod` for the layered comparisons
- `xarray`/`h5netcdf` for the Marlim data-generation scripts

The GNU build uses MUMPS from `$(HOME)/Install/MUMPS_5.7.3` and system
ParMETIS, Scotch, ScaLAPACK, LAPACK, BLAS, and FFTW3 paths. Adjust the
variables in `Makefile` for another installation. An Intel/MKL configuration
is provided as `Makefile_intel`:

```bash
make -f Makefile_intel
```

## Build

From the repository root:

```bash
make
```

This creates `./main`. The default GNU build uses `mpicc`, `-O3 -g -Wall`,
OpenMP, and `-DAdd_`. Clean generated objects and the executable with:

```bash
make clean
```

## Running the solver

The executable uses `key=value` command-line parameters. Run it from the
directory containing the model and acquisition files, because file names are
interpreted relative to the current working directory:

```bash
cd run_layered
mpirun -np 1 ../main \
    mode=1 airbc=1 fvm=0 preco=1 algopt=1 \
    tol=1e-8 niter=1000 \
    freqs=0.25,0.75,1.25 \
    chsrc=Ex chrec=Ex,Ey,Ez nb=10 \
    fx1=fx1 fx2=fx2 fx3=fx3 \
    frho11=frho11 frho22=frho22 frho33=frho33 \
    fsrc=sources.txt frec=receivers.txt
```

The complete commands used by the supplied examples are in:

```text
run_layered/run.sh
run_marlim/run.sh
bathymetry/run.sh
```

The direct command above uses the intended three-frequency layered run.
The current `run_layered/run.sh` contains `freqs=1.25 0.25,0.75,1.25`,
which is split into two shell arguments and does not pass the full frequency
vector to the parser. Use `freqs=0.25,0.75,1.25` (or fix the script) when
reproducing the three-frequency plots.

For a single-process run, the MPI executable can also be invoked directly:

```bash
cd run_marlim
../main mode=1 airbc=1 preco=1 algopt=1 \
    tol=1e-8 niter=1000 freqs=0.125,0.25,0.5,0.75,1.0,1.25 \
    chsrc=Ex chrec=Ex,Ey,Ez nb=10 lextend=120e3 \
    fx1=fx1 fx2=fx2 fx3=fx3 \
    frho11=frho11 frho22=frho22 frho33=frho33 \
    fsrc=sources.txt frec=receivers.txt
```

## Solver parameters

| Parameter | Meaning | Default |
| --- | --- | --- |
| `mode` | Modelling mode (`0` modelling, `1` FWI, `2` gradient-only in the data structure) | `0` |
| `preco` | Preconditioner: `0` none, `1` LU-SGS, `2` inner GMRES, `3` A-phi | `1` |
| `algopt` | Krylov solver: `1` BiCGStab, `2` GMRES | `1` |
| `niter` | Maximum iterations; GMRES uses `niter/nrestart` restart cycles | `1000` |
| `nrestart` | GMRES restart length | `10` |
| `mp` | Inner GMRES dimension when `preco=2` | `30` |
| `tol` | Base relative tolerance; solve tolerance is `tol/frequency` | `1e-8` |
| `verb` | Verbosity; defaults to `1` on rank 0 and `0` elsewhere | rank-dependent |
| `fvm` | Discretization scaling: `0` finite difference, `1` finite volume | `0` |
| `airbc` | `1` uses air-water/Earth boundary treatment for a model without air; `0` includes air | `0` |
| `nb` | Padding cells added on each side | `15` |
| `lextend` | Physical width used to construct the padding region | `40e3` when `nb>0` |
| `istretch` | `1` geometrically stretches padding cells; `0` uses constant spacing | `1` |
| `rho_air` | Air resistivity when air is included in the extension | `1e6` |
| `reciprocity` | Swap source and receiver locations | `0` |
| `param_extended` | Write the extended grid and resistivity files | `0` |

Frequencies are sorted internally. Active channels may be any comma-separated
combination of `Ex`, `Ey`, `Ez`, `Hx`, `Hy`, and `Hz`.

## Input files

### Grid and resistivity

`fx1`, `fx2`, and `fx3` are binary `float32` coordinate-node arrays. Their
lengths are `nx+1`, `ny+1`, and `nz+1`. Coordinates must be ascending.

`frho11`, `frho22`, and `frho33` are binary `float32` cell-centred
resistivity arrays. Each contains `nx*ny*nz` values, stored in C order with
the x index fastest, followed by y and z. The three files represent the
diagonal anisotropic resistivities.

The checked-in templates contain these model sizes:

| Example | `nx` | `ny` | `nz` | Cells |
| --- | ---: | ---: | ---: | ---: |
| `run_layered` | 100 | 100 | 100 | 1,000,000 |
| `run_marlim` | 180 | 120 | 120 | 2,592,000 |
| `bathymetry` | 125 | 125 | 100 | 1,562,500 |

The solver adds `nb` cells to each side. If `param_extended=1` and
`verb>0`, it writes `fx1_extended`, `fx2_extended`, `fx3_extended` and
`frho11_extended`, `frho22_extended`, `frho33_extended` as binary `float64`
arrays.

### Acquisition

`fsrc` and `frec` are text files. The first line is skipped as a header; each
following line contains:

```text
x1  x2  x3  azimuth  dip  index
```

The current driver uses one source (`nsrc=1`) and reads all receiver rows.
The checked-in templates pass `src_rec_table.txt`, but the current
matrix-free driver does not read that parameter or file.

## Outputs

The solver writes `emf_0001.txt` in the current working directory. Its columns
are:

```text
iTx  iRx  ichrec  ifreq  emf_real  emf_imag
```

`ifreq` and `ichrec` start at 1. The complex response is
`emf_real + i*emf_imag`.

The iterative solvers also write convergence histories such as
`iterate_fbicgstab.txt` or frequency-specific files, depending on the solver
and run configuration. Existing template output files may be overwritten.

## Examples

### Layered model

`run_layered` creates a one-line acquisition with
`create_acquisition_oneline.f90`, runs three frequencies, and includes
comparisons against 1-D `empymod` and 3-D `emg3d` responses:

```bash
cd run_layered
gfortran create_acquisition_oneline.f90 -o create_acquisition
./create_acquisition
sh run.sh
python3 plot_emdata.py
python3 plot_error_airbc.py
python3 plot_normalized_misfit.py
```

`make_model_1d.py` documents how the binary layered model can be regenerated.
The comparison scripts expect the corresponding Python packages and generated
solver output.

### Marlim survey

`run_marlim` contains a 180 x 120 x 120 survey model and six frequencies:

```bash
cd run_marlim
sh run.sh
python3 plot_emdata.py
python3 comparison.py
```

The checked-in `run.sh` uses `fx1`, `fx2`, `fx3`, and
`frho11`/`frho22`/`frho33`. `create_model.py` currently writes `x1`, `x2`,
`x3`, `Rh`, and `Rv`, so those names must be converted or the generator
updated before using it to replace the checked-in model files.
`create_acqui.py` generates acquisition text files from `marlim_survey.nc`.

The large source-data files `marlim_comp.npz`, `Rh`, and `Rv` are intentionally
excluded from Git because GitHub rejects files larger than 100 MB. Keep local
copies if you need to regenerate the Marlim model; the checked-in converted
`fx*` and `frho*` files are sufficient to run `run.sh`.

### Bathymetry

`bathymetry` contains a 125 x 125 x 100 model with topographic receiver
elevations:

```bash
cd bathymetry
sh run.sh
python3 plot_emdata_mare2dem.py
python3 plot_fdmesh.py
```

The directory also contains Mare2DEM comparison inputs and plotting
utilities. `create_emdata_file.f90` is an auxiliary comparison-data generator,
not a prerequisite for the CSEM solve.

The `nx`, `ny`, and `nz` arguments shown in `bathymetry/run.sh` are not read by
the current driver; dimensions are inferred from the binary grid files.

## Coordinate and boundary notes

- Coordinates are passed in metres.
- The model files contain resistivity, not conductivity.
- `nb` controls the padded computational domain; `lextend` and `istretch`
  control the nonuniform padding grid.
- With `airbc=1`, the input model is expected not to contain an air layer and
  `nb` must be nonzero.
- Run each case from its own example directory so relative paths and output
  files remain local to that case.
  
## Developer contact

Pengliang Yang

Email: ypl.2100@gmail.com
  
