using MPI, Test
include("../testhelpers.jl")

@testset "Balance Law Solver" begin
    tests = [
        (1, "Euler/isentropic_vortex_standalone.jl")
        (1, "Euler/isentropic_vortex_standalone_aux.jl")
        (1, "Euler/isentropic_vortex_standalone_source.jl")
        (1, "Euler/isentropic_vortex_standalone_bc.jl")
        # (1, "Euler/isentropic_vortex_standalone_integral.jl")
        (1, "Euler/isentropic_vortex_standalone_IMEX.jl")
        (1, "util/grad_test.jl")
        (1, "util/grad_test_sphere.jl")
        (1, "util/integral_test.jl")
        (1, "util/integral_test_sphere.jl")
        (1, "conservation/sphere.jl")
        (1, "compressible_Navier_Stokes/mms_bc.jl")
        (1, "soundings/runtests.jl")
    ]

    runmpi(tests, @__FILE__)

    if "linux" != lowercase(get(ENV, "TRAVIS_OS_NAME", ""))
        moretests = [
            (1, "../../examples/DGmethods_old/ex_001_periodic_advection.jl")
            (1, "../../examples/DGmethods_old/ex_002_solid_body_rotation.jl")
            (1, "../../examples/DGmethods_old/ex_003_acoustic_wave.jl")
            (1, "../../examples/DGmethods_old/ex_004_nonnegative.jl")
        ]
        runmpi(moretests, @__FILE__)
    end
end
