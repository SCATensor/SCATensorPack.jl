"""Build Right orthogonalization and left orthogonalization of MPS with QR"""
using LinearAlgebra

"""
Right_orth(mps,K)
Computes the QR factorization from left to right of given mps structure  with K sites
"""

function Right_orth(mps::MPS,K)::MPS
    nsites = convert(UInt16, K)
    for i=1:nsites-1
         Right_orth_local!(mps.X[i],mps.X[i+1])
    end
end

"""
Right_orth_local(mps1,mps2)
Computes the QR factorization of the local tensor core mps1 and updates mps2
"""

function Right_orth_local!(mps1::Tuple{DiagonalMPS,DiagonalMPS},
    mps2::Tuple{DiagonalMPS,DiagonalMPS})
    Block::Array{Matrix{Float64}}
    N1=mps1[1].NBlocks
    N2=mps1[2].NBlocks


    L1=mps1[1].BlockIndex
    L2=[(i[1]+1,i[1]+1) for i in mps1[2].BlockIndex ]
    L2=convert(Array{Tuple{UInt16,UInt16},1},L2)
    Block,idx=createBlocks_MPS(mps1,L1,L2)
    for i=1:length(Block)
        F=qr(Block[i])
        update!(mps1,F.Q,idx)
        diagonalRXMultiplication!(mps2, F.R)
    end
end

"""
createBlocks_MPS(mps1,L1,L2)

Computes the arrays of blocks in the tensor core unfolding matrix of the
tensor train mps1
mps1: tensor core with sparse block structure mps1[1] and mps1[2]
L1: List of blockindex of mps1[1]
L2: List of blockindex of mps1[2]
"""
function createBlocks_MPS(mps::Tuple{DiagonalMPS,DiagonalMPS},
    L1::Array{Tuple{UInt16, UInt16}},L2::Array{Tuple{UInt16, UInt16}})

    Block::Array{Matrix{Float64}}=[]
    idx=Array{UInt16,1}[]

    for j=1:length(L1)
        if L1[j] in L2
            id2=findall(in([L1[j]]),L2)[1]
            push!(idx,[j,id2])
            push!(Block,[mps[1].X[j];mps[2].X[id2]])
        else
            push!(idx,[j])
            push!(Block,mps[1].X[j])
        end
    end
    for j=1:length(L2)
        if L2[j]  ∉ L1
            push!(idx,[j])
            push!(Block,mps[2].X[j])
        end
    end
    return Block, idx
end
#=Elements present in 1 and not present in 2
L1_neg=findall(!in(L2),L1)
if !isempty(L1_neg)
    for i=1:length(L1_neg)
        F=qr(mps1[1].X[L1_neg[i]])
        mps1[1].X[L1_neg[i]]= Matrix(F.Q)
        diagonalRXMultiplication!(mps2, R)
    end
end
L2_neg=findall(!in(L2),L1)
if !isempty(L2_neg)
    for i=1:length(L2_neg)
        F=qr(mps1[1].X[L2_neg[i]])
        mps1[1].X[L2_neg[i]]= Matrix(F.Q)
        diagonalRXMultiplication!(mps2, R)
    end
end
L1_pos=findall(!in(L1),L2)
L2_pos=findall(!in(L2),L1)
for i=1:length(L1_pos)
    F=qr([mps1[1].X[L2_pos[i]];mps1[1].X[L1_pos[i]] ])
    Q= Matrix(F.Q)

    diagonalRXMultiplication!(mps2, R)
end=#
