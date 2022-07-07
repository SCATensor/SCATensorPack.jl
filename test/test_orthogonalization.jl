using Test
include("../src/DiagonalMPS.jl")
include("../src/MPS.jl")
include("../src/Orthogonalization.jl")

@testset "Right_QR" begin
    K=3
    N=2
    mps = RandomMPS(K, N)
    Right_orth!(mps,K);
    for i=1:K-1
        X=TT_core(mps.X[i])
    end
end
for i=1:K
print("Block indices", i, " 1: ", mps.X[i][1].BlockIndex ,"\n")
print("Block Sizes ", i, " 1: ", mps.X[i][1].BlockSizes,"\n")
print("Block indices", i, " 2: ", mps.X[i][2].BlockIndex,"\n")
print("Block Sizes ", i, " 2: ", mps.X[i][2].BlockSizes,"\n")
end
