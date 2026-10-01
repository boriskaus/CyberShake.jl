"""
    tools()

Names of the CyberShake executables shipped in `CyberShake_jll`.
"""
tools() = (:pmcl3d, :write_head, :compare_head, :strip_head, :search_head, :decimate_sgt,
           :generate_point_list, :reformat_awp, :reformat_awp_mpi, :reformat_velocity, :hf_synth)

function check_tool(name::Symbol)
    name in tools() || throw(ArgumentError("unknown CyberShake tool `$name`; available: $(join(tools(), ", "))"))
    check_available()
end

"""
    tool(name::Symbol) -> Cmd

Command for the CyberShake executable `name` (see [`tools`](@ref)), with the environment of
`CyberShake_jll` set up so that its libraries are found.
"""
function tool(name::Symbol)
    check_tool(name)
    return getfield(CyberShake_jll, name)()
end

"""
    tool_path(name::Symbol) -> String

Full path of the executable `name`.
"""
function tool_path(name::Symbol)
    check_tool(name)
    return getfield(CyberShake_jll, Symbol(name, "_path"))
end

"""
    run_tool(name::Symbol, args...; nprocs = 1, dir = pwd(), ignorestatus = false)

Run the CyberShake executable `name` with command-line arguments `args` in directory `dir`,
on `nprocs` MPI ranks. Returns the finished process.
"""
function run_tool(name::Symbol, args...; nprocs::Integer = 1, dir::AbstractString = pwd(), ignorestatus::Bool = false)
    cmd = Cmd(`$(mpi_command(name, nprocs)) $(collect(string.(args)))`; dir = dir)
    ignorestatus && (cmd = Base.ignorestatus(cmd))
    return run(cmd)
end

"""
    reformat_awp_sgt(input, nsamples, output; z_source = false)

Convert an SGT file written by `pmcl3d` (components `XX YY ZZ XY XZ YZ` per receiver, `nsamples`
time samples each) to the layout that the CyberShake post-processing expects: `XX` and `YY` are
swapped, and `XZ` and `YZ` are swapped and negated. With `z_source = true` all values are
additionally multiplied by -2, as required for the vertical source.
"""
function reformat_awp_sgt(input::AbstractString, nsamples::Integer, output::AbstractString; z_source::Bool = false)
    args = Any[input, nsamples, output]
    z_source && push!(args, "-z")
    run_tool(:reformat_awp, args...)
    return output
end
