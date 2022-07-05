# Test getindices function
using Test
using Printf
include("../src/DiagonalMPS.jl")

@testset "optimal blocks" begin
    function print_indices(N,K)
        @printf("N: %d K: %d",N,K)
        for i=1:K
            for j=1:2
                indices::Array{Tuple{Int64,Int64},1}= GetIndices(i,N,K,j)
                @printf("X[%d][%d]: ", i,j-1)
                for i = 1 : length(indices)
                    @printf("(%d, %d), ", indices[i][1], indices[i][2])
                end
                @printf("]\n")
            end
        end
    end
   #Example 1 (in https://arxiv.org/pdf/2104.13483.pdf)
   print_indices(2,5)
   print_indices(2,3)
   print_indices(3,6)
end
