"""
    check_available()

Throw an informative error unless `CyberShake_jll` has binaries for this platform and MPI setup.
"""
function check_available()
    CyberShake_jll.is_available() && return nothing
    if Sys.iswindows()
        error("CyberShake_jll has no Windows binaries (the CyberShake sources are POSIX-only). " *
              "Use Linux, macOS, or Windows Subsystem for Linux.")
    end
    error("CyberShake_jll is not available for this platform. It is only built against MPICH_jll; run " *
          "`using MPIPreferences; MPIPreferences.use_jll_binary(\"MPICH_jll\")` in your project, " *
          "restart Julia and try again.")
end

# `mpiexec` of whichever MPI the JLL was built against
const mpiexec, MPI_LIBPATH =
    if isdefined(CyberShake_jll, :MPICH_jll)
        CyberShake_jll.MPICH_jll.mpiexec(), CyberShake_jll.MPICH_jll.LIBPATH
    elseif isdefined(CyberShake_jll, :OpenMPI_jll)
        CyberShake_jll.OpenMPI_jll.mpiexec(), CyberShake_jll.OpenMPI_jll.LIBPATH
    elseif isdefined(CyberShake_jll, :MPItrampoline_jll)
        CyberShake_jll.MPItrampoline_jll.mpiexec(), CyberShake_jll.MPItrampoline_jll.LIBPATH
    else
        nothing, Ref("")
    end

"""
    mpi_available()

Whether an `mpiexec` was found, i.e. whether [`run_sgt`](@ref) can use more than one MPI rank.
"""
mpi_available() = mpiexec !== nothing

const pathsep = Sys.iswindows() ? ';' : ':'

"""
    mpi_command(name::Symbol, nprocs)

Command that runs the CyberShake executable `name` on `nprocs` MPI ranks.
"""
function mpi_command(name::Symbol, nprocs::Integer)
    nprocs == 1 && return tool(name)
    mpi_available() || error("No MPI launcher is available for this CyberShake_jll; use nprocs = 1")
    key = CyberShake_jll.JLLWrappers.JLLWrappers.LIBPATH_env
    libpath = join((CyberShake_jll.LIBPATH[], MPI_LIBPATH[]), pathsep)
    return `$(addenv(mpiexec, key => libpath)) -n $nprocs $(tool_path(name))`
end
