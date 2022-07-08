using Test
using BenchmarkTools
include("../src/DiagonalMPS.jl")
include("../src/MPS.jl")
include("../src/Orthogonalization.jl")
include("../src/Classic.jl")
include("../src/Classic.jl")

@testset "Left_rounding" begin
    K=20
    N=10
    mps = RandomMPS(K, N)
    Left_rounding!(mps,K);

    for i=1:K-1
        X=TT_core(mps.X[i])
        X=reshape(X,2*size(X,1),:)
        @test X'*X ≈ I
    end
end
