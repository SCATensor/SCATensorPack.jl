#""" The functions introduced here are only for comparison purposes"""
using TensorOperations

""" TT_qr_L!(mps,k)
returns left orthogonal k-1 tensor cores  where k is the number of sites and mps
is the tensor train (not considering the sparse block structure)
"""
function TT_qr_L!(mps::Array{AbstractArray{Float64}, 1},k)
    nsites=index = convert(UInt16, k)
    for i=1:nsites-1
        r_0=size(mps[i],1)
        r_1=size(mps[i],3)
        X=reshape(mps[i],2r_0,r_1)
        F=qr(X)
        Q=Matrix(F.Q)
        R=Matrix(F.R)
        mps[i]=reshape(Matrix(F.Q),r_0,2,size(Q,2))
        @tensor mps[i+1][α,s,γ] := R[α,β] * mps[i+1][β,s,γ]
    end
    nothing
end

""" TT_qr_R!(mps,k)
returns right orthogonal k-1 tensor cores  where k is the number of sites and mps
is the tensor train (not considering the sparse block structure)
"""
function TT_qr_R!(mps::Array{AbstractArray{Float64}, 1},k)
    nsites=index = convert(UInt16, k)
    for i=nsites:-1:2
        r_0=size(mps[i],1)
        r_1=size(mps[i],3)
        X=reshape(mps[i],r_0,2*r_1)
        F=lq(X)
        Q=Matrix(F.Q)
        R=Matrix(F.L)
        mps[i]=reshape(Matrix(F.Q),size(Q,1),2,r_1)
        @tensor mps[i-1][β,s,α] := mps[i-1][β,s,γ]* R[γ,α] 
    end
    nothing
end
"""
contract_mps(tensors)
returns the k-the order tensor in R^{2^k} (k is the number of sites) given
 the tensor train denoted by tensors
"""
function contract_mps(tensors::Array{AbstractArray{Float64},1})::AbstractArray{Float64}
    full_tensor = tensors[1]
    for i in 2:length(tensors)
        tensors[i]=permutedims(tensors[i],[2,1,3])
        full_tensor = ncon((full_tensor,
        tensors[i]),
        (push!(-collect(1:ndims(full_tensor)-1),ndims(full_tensor)),
        [-ndims(full_tensor)-1, ndims(full_tensor), -ndims(full_tensor) - 2]))
    end
    @assert size(full_tensor)[1] == 1
    @assert size(full_tensor)[end] == 1
    return reshape(full_tensor, size(full_tensor)[2:end-1])
end
