"""Build Right orthogonalization and left orthogonalization of MPS with QR"""
using LinearAlgebra

"""
Right_orth(mps,K)
Computes the QR factorization from left to right of given mps structure  with K sites
"""

function Right_orth!(mps::MPS,K)
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
    L1=mps1[1].BlockIndex
    L2=[(i[1]+1,i[1]+1) for i in mps1[2].BlockIndex ]
    L2=convert(Array{Tuple{UInt16,UInt16},1},L2)
    for j=1:length(L1)
        if L1[j] in L2
            id2=findall(in([L1[j]]),L2)[1]
            Block=[GetDiagonalBlock(mps1[2],id2);GetDiagonalBlock(mps1[1], j)]
            F=qr(Block)
            Update!(mps1,Matrix(F.Q),[id2,j],concat=true)
            DiagonalRXMultiplication!(mps2, F.R,L1[j])
        else
            Block=GetDiagonalBlock(mps1[1], j)
            F=qr(Block)
            Update!(mps1,Matrix(F.Q),[j],concat=false,flag=1)
            DiagonalRXMultiplication!(mps2, Matrix(F.R),L1[j])
        end
    end
    for j=1:length(L2)
            if L2[j]  ∉ L1
                Block=GetDiagonalBlock(mps1[2], j)
                F=qr(Block)
                Update!(mps1,F.Q,[j],concat=false,flag=2)
                DiagonalRXMultiplication!(mps2, Matrix(F.R),mps1[2].BlockIndex[j])
            end
    end
end

#="""
createBlocks_MPS(mps1,[(1,1),(2,2),(3,3)],[(2,2),(3,3)])

Computes the arrays of blocks   in the tensor core unfolding matrix of the
tensor train mps1 and returns idx=[[1],[2,1],[3,2]]
mps1: tensor core with sparse block structure mps1[1] and mps1[2]
L1: List of blockindex of mps1[1]
L2: List of blockindex of mps1[2]
"""

function createBlocks_MPS(mps::Tuple{DiagonalMPS,DiagonalMPS},
    L1::Array{Tuple{UInt16, UInt16}},L2::Array{Tuple{UInt16, UInt16}})

    Block::Array{Matrix{Float64}}=[]
    idx=Array{UInt16,1}[]
    pairs=Array{Tuple{UInt16,UInt16},1}()
    for j=1:length(L1)
        if L1[j] in L2
            id2=findall(in([L1[j]]),L2)[1]
            push!(idx,[id2,j])
            push!(Block,[GetDiagonalBlock(mps[2],id2);GetDiagonalBlock(mps[1], j)])
            push!(pairs,L1[j])
        else
            push!(idx,[j])
            push!(Block, GetDiagonalBlock(mps[1], j))
            push!(pairs,L1[j])
        end
    end
    for j=1:length(L2)
        if L2[j]  ∉ L1
            push!(idx,[j])
            push!(Block, GetDiagonalBlock(mps[2], j))
            push!(pairs,mps1[2].BlockIndex[j])
        end
    end
    return Block, idx,pairs
end=#

"""
Update!(mps1,Q,idx)

updates the tensor core blocks mps1[1] and mps1[2] by Q given the indices idx of the
matrices mps1[1].X[idx] to be replaced by orthogonal matrices Q
"""
function Update!(mps::Tuple{DiagonalMPS,DiagonalMPS},Q,id;concat=true,flag=1)
    Q_ = convert(Matrix{Float64}, Q)
    id_= convert(Array{UInt64,1}, id)
    if concat
        iter=Int(mps[1].BlockSizes[id_[1]][1])
        mps[2].X[id_[1]]=@view(Q_[1:iter,:])
        mps[1].X[id_[2]]=@view(Q_[iter+1:end,:])
    else
        if flag==1
            mps[1].X[id_[1]]=Q_
        else
            mps[2].X[id_[1]]=Q_
        end
    end
end


"""
DiagonalRXMultiplication!(mps[2], R ,(1,1))

Multiplies the provided R upper triangular matrix with the blocks of the next tensor core mps[2]
that have block indices (1,1),
"""
function DiagonalRXMultiplication!(mps::Tuple{DiagonalMPS,DiagonalMPS},R,pairs)
    R_ = convert(Matrix{Float64}, R)
    pairs_= convert(Tuple{UInt64,UInt64}, pairs)
    id1=findall(in(pairs_),mps[1].BlockIndex)
    id2=findall(in(pairs_),mps[2].BlockIndex)
    for i=1:length(id1)
        mps[1].X[id1[i]]= R_ * mps[1].X[id1[i]]
    end
    for i=1:length(id2)
        mps[2].X[id2[i]]=R_ * mps[2].X[id2[i]]
    end
end
