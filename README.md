# CyberShake.jl

Run the wave-propagation codes of the SCEC [CyberShake](https://github.com/SCECcode/cybershake-core)
platform from Julia, on Linux and macOS, without compiling anything.

The binaries come from [`CyberShake_jll`](https://github.com/boriskaus/CyberShake_jll.jl), which is built
with [BinaryBuilder](https://binarybuilder.org) from the Yggdrasil recipe `C/CyberShake`. They are:

| executable | what it does |
|---|---|
| `pmcl3d` | AWP-ODC finite-difference solver (MPI); computes strain Green tensors (SGTs) |
| `write_head`, `compare_head`, `strip_head`, `search_head`, `generate_point_list`, `decimate_sgt` | SGT header and point-list tools |
| `reformat_awp`, `reformat_awp_mpi`, `reformat_velocity` | convert AWP-ODC input/output files to the CyberShake layout |
| `hf_synth` | high-frequency seismogram synthesis |

Not included: the DirectSynth/RupGen codes (they need libmemcached, libcfu and embedded Python) and anything on
Windows (the sources are POSIX-only; use Windows Subsystem for Linux).

## Installation

`CyberShake_jll` is not in the General registry yet, so install it from GitHub first:

```julia
using Pkg
Pkg.add(url="https://github.com/boriskaus/CyberShake_jll.jl")
Pkg.add(url="https://github.com/boriskaus/CyberShake.jl")
```

`CyberShake_jll` is only built against MPICH. If you changed your MPI setup with
`MPIPreferences`, switch back with

```julia
using MPIPreferences
MPIPreferences.use_jll_binary("MPICH_jll")   # then restart Julia
```

## Quick start: a strain Green tensor on 2 MPI ranks

```julia
using CyberShake

model = SGTModel(
    n         = (40, 40, 40),      # grid points
    dh        = 100.0,             # grid spacing (m)
    dt        = 0.005,             # time step (s)
    nt        = 400,               # number of time steps
    vp = 3000.0, vs = 1500.0, rho = 2200.0,
    source    = (20, 20),          # point force at the free surface
    component = :x,                # direction of the force
    receivers = [(10, 20, 1), (30, 20, 1)],
    procs     = (2, 1, 1),         # MPI ranks in x, y, z
)

sgt = run_sgt(model)               # (nsamples, 6, nreceivers)
size(sgt)                          # (40, 6, 2)
sgt_components                     # ("XX", "YY", "ZZ", "XY", "XZ", "YZ")
```

`run_sgt` writes the solver input files into a temporary directory (pass `dir = "my_run"` to keep them), runs
`pmcl3d` through the `mpiexec` of the JLL and reads the output back. The solver log is in `pmcl3d.log`.
`write_inputs(model, dir)` only writes the inputs, so you can inspect or edit them
(`IN3D`, `input/media`, `input/src`, `input/cord`) and run any CyberShake executable yourself:

```julia
write_inputs(model, "my_run")
run_tool(:pmcl3d, "IN3D"; nprocs = 2, dir = "my_run")
run_tool(:hf_synth; ignorestatus = true) # prints which parameters are missing
tool(:write_head)                        # a Cmd with the library paths set, for use with `run`
reformat_awp_sgt("my_run/output_sgt/sgt", 40, "sgt_cybershake")   # 40 samples per component
```

### Notes and limitations

- The medium is homogeneous; for heterogeneous models write the media file yourself
  (`write_inputs` shows the format: one `(vp, vs, rho)` Float32 triplet per grid point) or use
  `reformat_velocity`.
- Results of runs with different `procs` agree to ~1e-7 in the interior of the model, but I saw differences
  of a few percent close to the absorbing outer boundaries when the domain is split in x or y. Compare
  receivers away from the boundaries when checking a parallel run against a serial one.
- Amplitudes scale linearly with the source time function `stf`; no physical normalisation is applied.

## Tests

```julia
using Pkg; Pkg.test("CyberShake")
```

The test suite checks that the JLL is available and every executable starts, validates the model
setup and input files, runs small SGT simulations and checks that they are finite, linear in the
source, causal (first arrivals ordered by distance, with an asymmetric source position so x/y mix-ups are
caught), identical when run on 1, 2 and 4 MPI ranks, and that `reformat_awp_sgt` permutes the components as documented.
GitHub Actions (`.github/workflows/CI.yml`) runs them on Linux (x86_64), macOS (Intel) and macOS (Apple silicon)
with Julia 1.11 and the latest release.

## Licenses

This package is MIT licensed. The CyberShake codes are distributed under the BSD 3-clause license by the
Southern California Earthquake Center; please cite Graves et al. (2011),
[doi:10.1007/s00024-010-0161-6](https://doi.org/10.1007/s00024-010-0161-6) when you use them.
