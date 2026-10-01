const sgt_components = ("XX", "YY", "ZZ", "XY", "XZ", "YZ")

"""
    gaussian_derivative(nt, dt; f0 = 3.0, t0 = 0.4, amplitude = 1.0)

Source time function with `nt` samples spaced `dt` seconds: the time derivative of a Gaussian,
peaking around frequency `f0` (Hz) and centred at time `t0` (s).
"""
function gaussian_derivative(nt::Integer, dt::Real; f0::Real = 3.0, t0::Real = 0.4, amplitude::Real = 1.0)
    return [amplitude * -2 * (π * f0)^2 * (t - t0) * exp(-(π * f0 * (t - t0))^2) for t in dt .* (0:nt-1)]
end

"""
    SGTModel(; kwargs...)

A small strain-Green-tensor simulation with the AWP-ODC finite-difference solver `pmcl3d`.

A point force acts at the free surface (`z = 1`) at grid position `source = (ix, iy)`. The solver
records the six strain-Green-tensor components at every `receivers` grid point `(ix, iy, iz)` (1-based,
`iz = 1` is the free surface), every `sgt_skip` time steps.

# Keywords
- `n = (40, 40, 40)`: number of grid points in x, y, z.
- `dh = 100.0`: grid spacing (m).
- `dt = 0.005`: time step (s). Must satisfy the CFL condition; the solver reports the stability number.
- `nt = 400`: number of time steps. Must be a multiple of `sgt_skip`.
- `vp = 3000.0`, `vs = 1500.0`, `rho = 2200.0`: homogeneous medium (m/s, m/s, kg/m³).
- `source = n[1]÷2, n[2]÷2`: surface position of the point force.
- `component = :x`: direction of the point force, `:x`, `:y` or `:z`.
- `stf = gaussian_derivative(nt, dt)`: source time function, one value per time step.
- `receivers`: vector of `(ix, iy, iz)` tuples; by default a 3×3 grid at the surface.
- `procs = (1, 1, 1)`: MPI ranks in x, y, z. Each of `n` must be divisible by the corresponding entry.
- `nd = 10`: thickness (grid points) of the absorbing Cerjan boundary layer.
- `sgt_skip = 10`: only every `sgt_skip`-th time step is recorded.
- `viscoelastic = true`: use the visco-elastic instead of the elastic scheme.
"""
struct SGTModel
    n::NTuple{3,Int}
    dh::Float64
    dt::Float64
    nt::Int
    vp::Float64
    vs::Float64
    rho::Float64
    source::NTuple{2,Int}
    component::Symbol
    stf::Vector{Float64}
    receivers::Vector{NTuple{3,Int}}
    procs::NTuple{3,Int}
    nd::Int
    sgt_skip::Int
    viscoelastic::Bool
end

function SGTModel(; n = (40, 40, 40), dh = 100.0, dt = 0.005, nt = 400,
                  vp = 3000.0, vs = 1500.0, rho = 2200.0,
                  source = (n[1] ÷ 2, n[2] ÷ 2), component = :x,
                  stf = gaussian_derivative(nt, dt),
                  receivers = [(x, y, 1) for x in n[1] ÷ 4:n[1] ÷ 4:3 * n[1] ÷ 4 for y in n[2] ÷ 4:n[2] ÷ 4:3 * n[2] ÷ 4],
                  procs = (1, 1, 1), nd = 10, sgt_skip = 10, viscoelastic = true)
    n, procs, source = NTuple{3,Int}(n), NTuple{3,Int}(procs), NTuple{2,Int}(source)
    receivers = NTuple{3,Int}[Tuple(r) for r in receivers]

    component in (:x, :y, :z) || throw(ArgumentError("component must be :x, :y or :z, got :$component"))
    length(stf) == nt || throw(ArgumentError("stf has $(length(stf)) samples, expected nt = $nt"))
    nt % sgt_skip == 0 || throw(ArgumentError("nt = $nt must be a multiple of sgt_skip = $sgt_skip"))
    isempty(receivers) && throw(ArgumentError("at least one receiver is required"))
    for i in 1:3
        n[i] % procs[i] == 0 || throw(ArgumentError("n[$i] = $(n[i]) is not divisible by procs[$i] = $(procs[i])"))
        n[i] ÷ procs[i] > nd || throw(ArgumentError("n[$i] ÷ procs[$i] must exceed the boundary layer thickness nd = $nd"))
    end
    all(1 .<= source .<= n[1:2]) || throw(ArgumentError("source $source is outside the grid $(n[1:2])"))
    all(r -> all(1 .<= r .<= n), receivers) || throw(ArgumentError("a receiver is outside the grid $n"))
    allunique(receivers) || throw(ArgumentError("receivers must be unique"))
    vs < vp || throw(ArgumentError("vs must be smaller than vp"))

    return SGTModel(n, dh, dt, nt, vp, vs, rho, source, component, Float64.(stf), receivers,
                    procs, nd, sgt_skip, viscoelastic)
end

nsamples(m::SGTModel) = m.nt ÷ m.sgt_skip

"""
    sgt_time(model)

Times (s) of the recorded SGT samples.
"""
sgt_time(m::SGTModel) = m.dt .* m.sgt_skip .* (1:nsamples(m))

function Base.show(io::IO, m::SGTModel)
    @printf(io, "SGTModel: %dx%dx%d grid, dh = %g m, %d steps of %g s, %d receivers, %d MPI ranks",
            m.n..., m.dh, m.nt, m.dt, length(m.receivers), prod(m.procs))
end

# source columns of pmcl3d's `igreen` modes 4, 5, 6
igreen(m::SGTModel) = 3 + findfirst(==(m.component), (:x, :y, :z))

"""
    write_inputs(model, dir)

Write the input files of `pmcl3d` (`IN3D`, media, source and receiver lists) and create its output
directories in `dir`.
"""
function write_inputs(m::SGTModel, dir::AbstractString)
    mkpath(dir)
    for d in ("input", "output_sfc", "output_vlm", "output_ckp", "output_dyn", "output_sgt")
        mkpath(joinpath(dir, d))
    end

    # one (vp, vs, rho) triplet per grid point
    open(joinpath(dir, "input", "media"), "w") do io
        write(io, repeat(Float32[m.vp, m.vs, m.rho], prod(m.n)))
    end

    # point-force position, then one line per time step with the six source components
    open(joinpath(dir, "input", "src"), "w") do io
        println(io, m.source[1], " ", m.source[2], " 1")
        col = igreen(m) - 3
        for v in m.stf
            vals = zeros(6)
            vals[col] = v
            println(io, join((@sprintf("%.8e", x) for x in vals), ", "))
        end
    end

    open(joinpath(dir, "input", "cord"), "w") do io
        println(io, length(m.receivers))
        for r in m.receivers
            println(io, join(r, " "))
        end
    end

    W = nsamples(m)   # all samples are written in a single chunk at the end of the run
    tmax = (m.nt - 1) * m.dt + m.dt / 2
    open(joinpath(dir, "IN3D"), "w") do io
        println(io, "$(igreen(m)) igreen")
        println(io, "$tmax TMAX")
        println(io, "$(m.dh) DH")
        println(io, "$(m.dt) DT")
        println(io, "0 NPC")
        println(io, "$(m.nd) ND")
        println(io, "0.92 ARBC")
        println(io, "0.1 PHT")
        println(io, "1 NSRC")
        println(io, "$(m.nt) NST")
        println(io, "$(m.n[1]) NX\n$(m.n[2]) NY\n$(m.n[3]) NZ")
        println(io, "$(m.procs[1]) NPX\n$(m.procs[2]) NPY\n$(m.procs[3]) NPZ")
        println(io, "0 IFAULT\n0 CHECKPOINT\n0 ISFCVLM\n0 IMD5\n1 IVELOCITY\n1 MEDIARESTART")
        println(io, "3 NVAR\n1 IOST\n1 PARTDEG\n1 IO_OPT\n1 PERF_MEAS\n0 IDYNA\n1 SOCALQ")
        println(io, "$(Int(m.viscoelastic)) NVE")
        println(io, "0.677 MU_S\n0.525 MU_D")
        println(io, "0.01 FL\n25.0 FH\n0.5 FP")
        println(io, "$(m.nt) READ_STEP\n$W WRITE_STEP\n$W WRITE_STEP2")
        # surface and volume seismogram output are switched off by the huge time skips
        println(io, "1 NBGX\n$(m.n[1]) NEDX\n10 NSKPX\n1 NBGY\n$(m.n[2]) NEDY\n10 NSKPY")
        println(io, "1 NBGZ\n$(m.n[3]) NEDZ\n10 NSKPZ\n100000 NTISKP")
        println(io, "1 NBGX2\n$(m.n[1]) NEDX2\n10 NSKPX2\n1 NBGY2\n$(m.n[2]) NEDY2\n10 NSKPY2")
        println(io, "1 NBGZ2\n$(m.n[3]) NEDZ2\n10 NSKPZ2\n1000000 NTISKP2")
        println(io, "$(m.sgt_skip) NTISKP_SGT")
        println(io, "'output_ckp/CHKP' CHKP")
        println(io, "'output_ckp/CHKJ' CHKJ")
        println(io, "'input/src' INSRC")
        println(io, "'input/media' INVEL")
        println(io, "'input/cord' INSGT")
        println(io, "'output_sfc/XX' SXRGO\n'output_sfc/YY' SYRGO\n'output_sfc/ZZ' SZRGO")
        println(io, "'output_vlm/XY' SXRGO2\n'output_vlm/YZ' SYRGO2\n'output_vlm/XZ' SZRGO2")
        println(io, "'output_sgt/sgt' SGTGRO")
        println(io, "'output_dyn/SGSN' SGSN")
    end
    return dir
end

"""
    read_sgt(model, dir) -> Array{Float32,3}

Read the SGT file written by `pmcl3d` in `dir`. The result has size
`(nsamples, 6, nreceivers)`; the second dimension holds the components `XX YY ZZ XY XZ YZ`
(see `sgt_components`).
"""
function read_sgt(m::SGTModel, dir::AbstractString)
    file = joinpath(dir, "output_sgt", "sgt")
    dims = (nsamples(m), 6, length(m.receivers))
    isfile(file) || error("pmcl3d did not write $file")
    filesize(file) == prod(dims) * sizeof(Float32) ||
        error("$file has $(filesize(file)) bytes, expected $(prod(dims) * sizeof(Float32))")
    return open(file) do io
        read!(io, Array{Float32}(undef, dims))
    end
end

"""
    run_sgt(model; dir = mktempdir()) -> Array{Float32,3}

Run `pmcl3d` for `model` on `prod(model.procs)` MPI ranks in directory `dir` and return the
strain Green tensors, an array of size `(nsamples, 6, nreceivers)` (see [`read_sgt`](@ref)).
The solver output is kept in `dir/pmcl3d.log`.
"""
function run_sgt(m::SGTModel; dir::AbstractString = mktempdir())
    check_available()
    write_inputs(m, dir)
    log = joinpath(dir, "pmcl3d.log")
    cmd = Cmd(`$(mpi_command(:pmcl3d, prod(m.procs))) IN3D`; dir = dir)
    try
        run(pipeline(cmd; stdout = log, stderr = log))
    catch err
        err isa ProcessFailedException || rethrow()
        tail = join(last(readlines(log), 20), "\n")
        error("pmcl3d failed; last lines of $log:\n$tail")
    end
    return read_sgt(m, dir)
end
