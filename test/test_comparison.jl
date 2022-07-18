#=Run comprison tests between Itensors and our code for
-MPS construction
-tensor train orthogonalisation
-tensor train compression
with respect to the number of sites/ranks
=#
import ITensors
using Random
using Test
using BenchmarkTools
include("../src/DiagonalMPS.jl")
include("../src/MPS.jl")
include("../src/Classic.jl")
include("../src/Orthogonalization.jl")
include("../src/Rounding.jl")

function wave_function(mps::MPS,k)
  K = convert(UInt16, k)
  mps_block::Array{AbstractArray{Float64}, 1} = []
  for i=1:K
    X=TT_core(mps.X[i])
    push!(mps_block,X)
  end
  mps_block_contracted=contract_mps(mps_block)
  return mps_block_contracted
end

@testset "ITensors  vs SB codes TT_SVD" begin
  iter=10
  ranks=[2^i for i=2:iter]
  K=5
  N=2

  max_ranks=[]
  SB_times_rounding=[]
  ITensors_times_rounding=[]
  max_Itensors=[]
  SB_error=[];IT_error=[]
  for j=1:length(ranks)
    #if ranks[j] <= (2^K/2)
    sites = ITensors.siteinds(2,K)
    mps_IT = ITensors.randomMPS(sites,ranks[j])
    println("Ranks ITensors:",ITensors.linkdims(mps_IT))
    push!(max_Itensors,maximum(ITensors.linkdims(mps_IT)))
    ψ=copy(mps_IT)
    T1=@elapsed ITensors.truncate!(mps_IT;cutoff =1E-12)
    push!(ITensors_times_rounding,T1)
    #E=norm(ψ-mps_IT)/norm(ψ)
    #push!(IT_error,E)


    mps = RandomMPS(K, N, ITensors.linkdims(ψ))
    println("Ranks:", mps.ranks)
    #ψ=wave_function(mps,K)
    T1=@elapsed Left_rounding!(mps,K,δ=1E-12);
    push!(SB_times_rounding,T1)
    #ψ_c=wave_function(mps,K)
    #E=norm(ψ .-ψ_c)/norm(ψ_c)
    #push!(SB_error,E)
    push!(max_ranks,maximum(mps.ranks))
    #end
  end
end

@testset "ITensors  vs SB codes TT_QR" begin
  iter=10
  ranks=[2^i for i=6:iter]
  K=10
  N=2

  max_ranks=[]
  SB_times_left=[]
  ITensors_times_left=[]
  max_Itensors=[]
  SB_error=[];IT_error=[]
  for j=1:length(ranks)
   #if ranks[j] <= (2^K/2)
    sites = ITensors.siteinds(2,K)
    mps_IT = ITensors.randomMPS(sites,ranks[j])
    println("Ranks ITensors:",ITensors.linkdims(mps_IT))
    push!(max_Itensors,maximum(ITensors.linkdims(mps_IT)))
    ψ=copy(mps_IT)
    T1=@elapsed ITensors.orthogonalize!(mps_IT,K,ortho="left")
    push!(ITensors_times_left,T1)
    #E=norm(ψ-mps_IT)/norm(ψ)
    #push!(IT_error,E)


    mps = RandomMPS(K, N, ITensors.linkdims(ψ))
    println("Ranks:", mps.ranks)
    ψ=wave_function(mps,K)
    T1= @elapsed Left_orth!(mps,K);
    push!(SB_times_left,T1)
    ψ_c=wave_function(mps,K)
    E=norm(ψ .-ψ_c)/norm(ψ_c)
    #push!(SB_error,E)
    #push!(max_ranks,maximum(mps.ranks))
    #end
  end

end
