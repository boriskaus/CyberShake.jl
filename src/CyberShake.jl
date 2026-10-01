"""
    CyberShake

Run the SCEC [CyberShake](https://github.com/SCECcode/cybershake-core) wave-propagation codes from Julia.

The binaries come from `CyberShake_jll` (Linux and macOS). The main entry point is [`run_sgt`](@ref),
which sets up and runs a small strain-Green-tensor (SGT) simulation with the AWP-ODC solver `pmcl3d`.
All other CyberShake executables are available through [`tool`](@ref) and [`run_tool`](@ref).
"""
module CyberShake

using CyberShake_jll
using Printf

export SGTModel, run_sgt, write_inputs, read_sgt, sgt_components, gaussian_derivative
export reformat_awp_sgt
export tool, run_tool, tools, mpi_available, check_available

include("mpi.jl")
include("tools.jl")
include("sgt.jl")

end
