using Test
using BenchmarkTools
include("../src/DiagonalMPS.jl")
include("../src/MPS.jl")
include("../src/Orthogonalization.jl")
include("../src/Classic.jl")

@testset "Left_QR" begin
    K=20
    N=10
    mps = RandomMPS(K, N)
    Left_orth!(mps,K);

    for i=1:K-1
        X=TT_core(mps.X[i])
        X=reshape(X,2*size(X,1),:)
        @test X'*X ≈ I
    end
end

@testset "Right_QR" begin
    K=3
    N=2
    mps = RandomMPS(K, N)
    Right_orth!(mps,K);

    for i=K:-1:2
        @show i
        X=TT_core(mps.X[i])
        X=reshape(X,size(X,1),:)
        @test X*X' ≈ I
    end
end

@testset "Time Comparision: Left_orthogonlization " begin
    K=3
    N=2
    mps = RandomMPS(K, N)
    #Construct TT mps
    tt_mps::Array{AbstractArray{Float64}, 1} = []
    for i=1:K
        X=TT_core(mps.X[i])
        push!(tt_mps,X)
    end

    #Benchmarking
    @printf("Time execution of Left QR on block sparse structure:")
    @btime Left_orth!($mps,$K);
    @printf("Time execution of classical Left QR:")
    @btime TT_qr_L!($tt_mps,$K);

    #Accuracy
    mps_block::Array{AbstractArray{Float64}, 1} = []
    for i=1:K
        X=TT_core(mps.X[i])
        push!(mps_block,X)
    end
    tt_mps_contracted= contract_mps(tt_mps)
    mps_block_contracted=contract_mps(mps_block)
    E=norm(tt_mps_contracted-mps_block_contracted)
    @printf("Error is: %g",E)
end

@testset "Time Comparision: Right_orthogonlization " begin
    K=100
    N=2
    mps = RandomMPS(K, N)
    #Construct TT mps
    tt_mps::Array{AbstractArray{Float64}, 1} = []
    for i=1:K
        X=TT_core(mps.X[i])
        push!(tt_mps,X)
    end

    #Benchmarking
    @printf("Time execution of Left QR on block sparse structure:")
    @btime Right_orth!($mps,$K);
    @printf("Time execution of classical Left QR:")
    @btime TT_qr_R!($tt_mps,$K);

    #Accuracy
    mps_block::Array{AbstractArray{Float64}, 1} = []
    for i=1:K
        X=TT_core(mps.X[i])
        push!(mps_block,X)
    end
    tt_mps_contracted= contract_mps(tt_mps)
    mps_block_contracted=contract_mps(mps_block)
    E=norm(tt_mps_contracted-mps_block_contracted)
    @printf("Error is: %g",E)
end
