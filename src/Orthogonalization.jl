#+"""Build Left/right orthogonalization  of MPS cores with QR"""
using LinearAlgebra
"""
Left_orth(mps,K)
Computes the QR factorization from left to right of given mps structure  with K sites
"""
function Left_orth!(mps::MPS,K)
    nsites = convert(UInt16, K)
    for i=1:nsites-1
         Left_orth_local!(mps.X[i],mps.X[i+1])
    end
    nothing
end

"""
Left_orth_local(mps1,mps2)
Computes the QR factorization of the local tensor core mps1 and updates mps2
"""

function Left_orth_local!(mps1::Tuple{DiagonalMPS,DiagonalMPS},
    mps2::Tuple{DiagonalMPS,DiagonalMPS})
    L1=mps1[1].BlockIndex
    L2=mps1[2].BlockIndex
    for j=1:length(L1)
        if getfield.(L1, 2)[j] in getfield.(L2, 2)
            id2=findall(isequal(getfield.(L1, 2)[j]), getfield.(L2, 2))[1]
            Block=[GetDiagonalBlock(mps1[2],id2);GetDiagonalBlock(mps1[1], j)]
            F=qr(Block)
            Update_left!(mps1,Matrix(F.Q),[id2,j],concat=true)
            Left_DiagonalRXMultiplication!(mps2, Matrix(F.R),getfield.(L1, 2)[j])
        else
            Block=GetDiagonalBlock(mps1[1], j)
            F=qr(Block)
            Update_left!(mps1,Matrix(F.Q),[j],concat=false,flag=1)
            Left_DiagonalRXMultiplication!(mps2, Matrix(F.R),getfield.(L1, 2)[j])
        end
    end
    for j=1:length(L2)
            if getfield.(L2, 2)[j]   ∉ getfield.(L1, 2)
                Block=GetDiagonalBlock(mps1[2], j)
                F=qr(Block)
                Update_left!(mps1,Matrix(F.Q),[j],concat=false,flag=2)
                Left_DiagonalRXMultiplication!(mps2, Matrix(F.R),getfield.(L2, 2)[j])
            end
    end
    nothing
end


"""
Right_orth(mps,K)
Computes the QR factorization from right to left of given mps structure  with K sites
"""

function Right_orth!(mps::MPS,K)
    nsites = convert(UInt16, K)
    for i=nsites:-1:2
         Right_orth_local!(mps.X[i],mps.X[i-1])
    end
    nothing
end

"""
Right_orth_local(mps1,mps2)
Computes the QR factorization of the local tensor core mps1 and updates mps2
"""

function Right_orth_local!(mps1::Tuple{DiagonalMPS,DiagonalMPS},
    mps2::Tuple{DiagonalMPS,DiagonalMPS})
    L1=mps1[1].BlockIndex
    L2=mps1[2].BlockIndex
    for j=1:length(L1)
        if getfield.(L1, 1)[j] in getfield.(L2, 1)
            id2=findall(isequal(getfield.(L1, 1)[j]), getfield.(L2, 1))[1]
            Block=[GetDiagonalBlock(mps1[2],id2) GetDiagonalBlock(mps1[1], j)]
            F=lq(Block)
            Update_right!(mps1,Matrix(F.Q),[id2,j],concat=true)
            Right_DiagonalRXMultiplication!(mps2, Matrix(F.L),getfield.(L1, 1)[j])
        else
            Block=GetDiagonalBlock(mps1[1], j)
            F=lq(Block)
            Update_right!(mps1,Matrix(F.Q),[j],concat=false,flag=1)
            Right_DiagonalRXMultiplication!(mps2, Matrix(F.L),getfield.(L1, 1)[j])
        end
    end
    for j=1:length(L2)
            if getfield.(L2, 1)[j]   ∉ getfield.(L1, 1)
                Block=GetDiagonalBlock(mps1[2], j)
                F=lq(Block)
                Update_right!(mps1,Matrix(F.Q),[j],concat=false,flag=2)
                Right_DiagonalRXMultiplication!(mps2, Matrix(F.L),getfield.(L2, 1)[j])
            end
    end
    nothing
end

"""
Update!(mps1,Q,idx)

updates the tensor core blocks mps1[1] and mps1[2] by Q given the indices idx of the
matrices mps1[1].X[idx] to be replaced by orthogonal matrices Q
"""
function Update_left!(mps::Tuple{DiagonalMPS,DiagonalMPS},Q,id;concat=true,flag=1)
    Q_ = convert(Matrix{Float64}, Q)
    id_= convert(Array{UInt64,1}, id)
    if concat
        iter=Int(mps[2].BlockSizes[id_[1]][1])
        mps[2].X[id_[1]]=@view(Q_[1:iter,:])
        mps[1].X[id_[2]]=@view(Q_[iter+1:end,:])
        mps[1].BlockSizes[id_[2]]=UInt16.(size(Q_[iter+1:end,:]))
        mps[2].BlockSizes[id_[1]]=UInt16.(size(Q_[1:iter,:]))
    else
        if flag==1
            mps[1].X[id_[1]]=Q_
            mps[1].BlockSizes[id_[1]]=UInt16.(size(Q_))
        else
            mps[2].X[id_[1]]=Q_
            mps[2].BlockSizes[id_[1]]=UInt16.(size(Q_))
        end
    end
    nothing
end

function Update_right!(mps::Tuple{DiagonalMPS,DiagonalMPS},Q,id;concat=true,flag=1)
    Q_ = convert(Matrix{Float64}, Q)
    id_= convert(Array{UInt64,1}, id)
    if concat
        iter=Int(mps[1].BlockSizes[id_[2]][2])
        mps[1].X[id_[2]]=@view(Q_[:,1:iter])
        mps[2].X[id_[1]]=@view(Q_[:,iter+1:end])
        mps[2].BlockSizes[id_[1]]=UInt16.(size(Q_[:,iter+1:end]))
        mps[1].BlockSizes[id_[2]]=UInt16.(size(Q_[:,1:iter]))
    else
        if flag==1
            mps[1].X[id_[1]]=Q_
            mps[1].BlockSizes[id_[1]]=UInt16.(size(Q_))
        else
            mps[2].X[id_[1]]=Q_
            mps[2].BlockSizes[id_[1]]=UInt16.(size(Q_))
        end
    end
    nothing
end

"""
DiagonalRXMultiplication!(mps[2], R ,(1,1))

Multiplies the provided R upper triangular matrix with the blocks of the next tensor core mps[2]
that have block indices (1,1),
"""
function Left_DiagonalRXMultiplication!(mps::Tuple{DiagonalMPS,DiagonalMPS},R,idx)
    R_ = convert(Matrix{Float64}, R)
    id1=findall(isequal(idx), getfield.(mps[1].BlockIndex, 1))
    id2=findall(isequal(idx), getfield.(mps[2].BlockIndex, 1))
    if !isempty(id1)
        for k=1:length(id1)
            mps[1].X[id1[k]]= R_ * mps[1].X[id1[k]]
            mps[1].BlockSizes[id1[k]]=size(mps[1].X[id1[k]])
        end
    end
    if !isempty(id2)
        for k=1:length(id2)
            mps[2].X[id2[k]]= R_ * mps[2].X[id2[k]]
            mps[2].BlockSizes[id2[k]]=UInt16.(size(mps[2].X[id2[k]]))
        end
    end

    nothing
end

function Right_DiagonalRXMultiplication!(mps::Tuple{DiagonalMPS,DiagonalMPS},R,idx)
    R_ = convert(Matrix{Float64}, R)
    id1=findall(isequal(idx), getfield.(mps[1].BlockIndex, 2))
    id2=findall(isequal(idx), getfield.(mps[2].BlockIndex, 2))
    if !isempty(id1)
        for k=1:length(id1)
            mps[1].X[id1[k]]=  mps[1].X[id1[k]] * R_
            mps[1].BlockSizes[id1[k]]=size(mps[1].X[id1[k]])
        end
    end
    if !isempty(id2)
        for k=1:length(id2)
            mps[2].X[id2[k]]= mps[2].X[id2[k]] * R_
            mps[2].BlockSizes[id2[k]]=UInt16.(size(mps[2].X[id2[k]]))
        end
    end

    nothing
end
