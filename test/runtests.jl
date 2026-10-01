using Test
using CyberShake
import CyberShake_jll

# Euclidean norm over the six SGT components at every sample of one receiver
trace_norm(sgt, r) = vec(sqrt.(sum(abs2, sgt[:, :, r]; dims = 2)))

# index of the first sample above `frac` of the peak
first_arrival(x; frac = 0.01) = findfirst(>(frac * maximum(x)), x)

@testset "CyberShake.jl" begin

    @testset "JLL" begin
        @test CyberShake_jll.is_available()
        @test isnothing(check_available())
        @test mpi_available()
        for name in tools()
            @test isfile(CyberShake.tool_path(name))
        end
        @test_throws ArgumentError tool(:not_a_tool)
    end

    @testset "executables start" begin
        # Without arguments every tool prints a usage message or a complaint about missing
        # parameters and exits with a small status; a missing library or a crash would show
        # up as a signal or as the loader's exit code 127.
        for name in setdiff(tools(), (:pmcl3d,))
            out = IOBuffer()
            p = run(pipeline(ignorestatus(tool(name)); stdout = out, stderr = out))
            @test p.termsignal == 0
            @test p.exitcode in 0:100
            @test !isempty(String(take!(out)))
        end
    end

    @testset "SGTModel validation" begin
        @test SGTModel() isa SGTModel
        @test_throws ArgumentError SGTModel(component = :w)
        @test_throws ArgumentError SGTModel(nt = 405)                        # not a multiple of sgt_skip
        @test_throws ArgumentError SGTModel(procs = (3, 1, 1))               # 40 not divisible by 3
        @test_throws ArgumentError SGTModel(procs = (4, 1, 1))               # subdomain thinner than nd
        @test_throws ArgumentError SGTModel(receivers = [(1, 1, 99)])        # outside the grid
        @test_throws ArgumentError SGTModel(source = (0, 5))
        @test_throws ArgumentError SGTModel(vs = 4000.0)
        @test_throws ArgumentError SGTModel(stf = zeros(3))
    end

    @testset "input files" begin
        m = SGTModel(n = (30, 40, 30), nt = 100, procs = (1, 2, 1), source = (12, 28), component = :y,
                     receivers = [(5, 6, 1), (7, 8, 3)])
        dir = write_inputs(m, mktempdir())
        in3d = read(joinpath(dir, "IN3D"), String)
        @test startswith(in3d, "5 igreen")                                    # y force
        @test occursin("30 NX\n40 NY\n30 NZ", in3d)
        @test occursin("1 NPX\n2 NPY\n1 NPZ", in3d)
        @test filesize(joinpath(dir, "input", "media")) == 3 * 4 * 30 * 40 * 30
        src = readlines(joinpath(dir, "input", "src"))
        @test src[1] == "12 28 1"
        @test length(src) == 1 + 100
        @test all(l -> count(',', l) == 5, src[2:end])
        # the force sits in the second (yy) column only
        @test all(l -> split(l, ", ")[[1, 3, 4, 5, 6]] == fill("0.00000000e+00", 5), src[2:end])
        @test readlines(joinpath(dir, "input", "cord")) == ["2", "5 6 1", "7 8 3"]
    end

    @testset "SGT simulation" begin
        m1 = SGTModel()
        sgt1 = run_sgt(m1)
        nsamp = m1.nt ÷ m1.sgt_skip
        @test size(sgt1) == (nsamp, 6, 9)
        @test eltype(sgt1) == Float32
        @test all(isfinite, sgt1)
        @test maximum(abs, sgt1) > 0

        @testset "parallel runs give the same result" begin
            # Receivers straddle the subdomain interfaces but stay well away from the absorbing
            # outer boundaries, where serial and parallel runs differ slightly.
            kw = (n = (80, 80, 40), source = (30, 30),
                  receivers = [(x, y, 1) for x in (30, 39, 40, 41, 50) for y in (30, 40, 41)])
            ref = run_sgt(SGTModel(; kw...))
            @test all(isfinite, ref)
            for procs in ((2, 1, 1), (1, 2, 1), (1, 1, 2), (2, 2, 1))
                sgt = run_sgt(SGTModel(; kw..., procs))
                @test sgt ≈ ref rtol = 1e-4 atol = 1e-5 * maximum(abs, ref)
            end
        end

        @testset "linear in the source amplitude" begin
            m3 = SGTModel(stf = 3 .* gaussian_derivative(400, 0.005))
            @test run_sgt(m3) ≈ 3 .* sgt1 rtol = 1e-4 atol = 1e-6 * maximum(abs, sgt1)
        end

        @testset "waves travel outwards from the source" begin
            # asymmetric source position, so that swapping x and y would change the answer:
            # receiver A is 6 points from the source, receiver B 18 points
            m4 = SGTModel(source = (12, 28), receivers = [(18, 28, 1), (12, 10, 1)], nt = 600,
                          stf = gaussian_derivative(600, 0.005; t0 = 0.3))
            sgt4 = run_sgt(m4)
            @test all(isfinite, sgt4)
            tA = first_arrival(trace_norm(sgt4, 1))
            tB = first_arrival(trace_norm(sgt4, 2))
            @test tA < tB
            # nothing may arrive faster than the P wave
            dist = 18 * m4.dh
            @test tB * m4.dt * m4.sgt_skip > dist / m4.vp - 0.2
        end

        @testset "sgt_time" begin
            @test CyberShake.sgt_time(m1) ≈ 0.05:0.05:2.0
        end
    end

    @testset "reformat_awp_sgt" begin
        m = SGTModel()
        sgt = run_sgt(m)
        dir = mktempdir()
        infile = joinpath(dir, "in.sgt")
        open(io -> write(io, sgt), infile, "w")
        nsamp = size(sgt, 1)

        for z in (false, true)
            outfile = joinpath(dir, "out_$z.sgt")
            reformat_awp_sgt(infile, nsamp, outfile; z_source = z)
            out = open(io -> read!(io, similar(sgt)), outfile)
            f = z ? -2 : 1
            # XX <-> YY, ZZ and XY unchanged, XZ <- -YZ, YZ <- -XZ
            @test out[:, 1, :] == f .* sgt[:, 2, :]
            @test out[:, 2, :] == f .* sgt[:, 1, :]
            @test out[:, 3, :] == f .* sgt[:, 3, :]
            @test out[:, 4, :] == f .* sgt[:, 4, :]
            @test out[:, 5, :] == -f .* sgt[:, 6, :]
            @test out[:, 6, :] == -f .* sgt[:, 5, :]
        end
    end
end
