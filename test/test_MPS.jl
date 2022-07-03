#=
PrintMPS(mps)

Prints on screen the formatted MPS structure.
=#
using Test
include("../src/DiagonalMPS.jl")
include("../src/MPS.jl")

@testset "MPS" begin
mps = RandomMPS(3, 2)
PrintMPS(mps)

end
