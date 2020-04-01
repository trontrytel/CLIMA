using Test

for submodule in ["unit_tests", "ex_1_saturation_adjustment", "ex_2_Kessler"]
    println("Testing $submodule")
    include(joinpath(submodule * ".jl"))
end
