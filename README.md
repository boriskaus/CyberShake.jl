# CyberShake.jl

[![CI](https://github.com/boriskaus/CyberShake.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/boriskaus/CyberShake.jl/actions/workflows/CI.yml)

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

## Credits and how to cite

**CyberShake.jl only provides a Julia interface. All the science is done by the CyberShake software of the
[Southern California Earthquake Center (SCEC)](https://scec.org)**, <https://github.com/SCECcode/cybershake-core>
(BSD 3-Clause license, Copyright (c) 2022, Southern California Earthquake Center). See the
[CyberShake wiki](https://github.com/SCECcode/cybershake-core/wiki) for its documentation and
[CyberShake study descriptions](https://strike.scec.org/scecpedia/Comparison_of_CyberShake_Studies) for the datasets.
This package uses the version tagged `study_24_8`.

The CyberShake authors ask users to cite the CyberShake paper (their words: "Cite Code As" and "Primary Reference"):

> Graves, R., Jordan, T.H., Callaghan, S. et al. CyberShake: A Physics-Based Seismic Hazard Model for
> Southern California. *Pure Appl. Geophys.* 168, 367–381 (2011; first published online 2010).
> <https://doi.org/10.1007/s00024-010-0161-6>, SCEC Contribution 1354.

suggested text for the body of your paper,

> The research described in this article used CyberShake software (Graves et al., 2011) published under the BSD-3 license.

and this acknowledgement:

> We would like to acknowledge use of the CyberShake software provided by the Southern California Earthquake Center
> (http://scec.org) which is funded by NSF Cooperative Agreement EAR-1600087 and USGS Cooperative Agreement G17AC00047.

`CyberShake.citation()` prints this text, and [`CITATION.bib`](CITATION.bib) / [`CITATION.cff`](CITATION.cff) contain
machine-readable versions.

**Authors of the CyberShake software** (from the upstream [CREDITS.md](https://github.com/SCECcode/cybershake-core/blob/main/CREDITS.md)):
Scott Callaghan, Kevin Milner, Robert Graves, Kim Olsen, Philip Maechling, Fabio Silva, Thomas H. Jordan,
Christine Goulet, Karan Vahi, Patrick Small, Matt Rynge, Ewa Deelman and Hunter Francoeur. The AWP-ODC solver `pmcl3d` goes back
to C. Marcinkovich and K.B. Olsen (2004), with major contributions from Y.F. Cui and many others, as listed in the
header of its source file `AWP-ODC-SGT/src/pmcl3d.f`. The command-line parser `libget` that the tools link
statically is Copyright 1990 Science Applications International Corporation.

## Licenses

This Julia package is MIT licensed. The CyberShake executables in `CyberShake_jll` are built unmodified from
cybershake-core (BSD 3-Clause, see above; the license text is installed with the JLL in `share/licenses/CyberShake`),
apart from a one-line patch that removes an assignment to an unallocated array in `pmcl3d.f` which crashes
with gfortran. The JLL is built with [BinaryBuilder](https://binarybuilder.org) from the
[Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil) recipe `C/CyberShake`.
