using Printf
using BlockArrays

abstract type AbstractMPS end

"""
An MPS structure represents a set of two DiagonalMPS X_{k}[1] and X_{k}[2]
structures. Block elements are stored in the matrix array X.
"""
mutable struct MPS <: AbstractMPS
    KSites::UInt16
    NParticles::UInt16
    X::Array{Tuple{DiagonalMPS, DiagonalMPS}}

    """
    mps = MPS(2, 3, [(diag00MPS, diag01MPS), (diag10MPS, diag11MPS), (diag20MPS, diag21MPS)])

    Constructs a new MPS structure with two sites and three particles. Each of the
    DiagonalMPS passed as third argument should be previously created.
    """
    function MPS(k, n, x)
        new(k, n, x)
    end
end

"""
mps = MPS(2, 3)

Constructs a new MPS structure with two sites and three particles. Each of the
DiagonalMPS contained in the MPS are set to "undef" values.
"""
function MPS(k, n)
#=     if !isinteger(k)
        error("$(typeof(k)) is not accepted as KSites attribute. It should be a positive integer.")
        return nothing
    end

    if !isinteger(n)
        error("$(typeof(n)) is not accepted as NParticles attribute. It should be a positive integer.")
        return nothing
    end
 =#
    K = convert(UInt16, k)
    N = convert(UInt16, n)
    X_ = Array{Tuple{DiagonalMPS, DiagonalMPS}}(undef, K)

    return MPS(K, N, X_)
end

"""
mps = RandomMPS(2, 3)

Constructs a new MPS structure with two sites and three particles. Each of the
DiagonalMPS contained in the MPS are generated randomly; this is, with random
sizes in the range (2 : 2^3) and elements in the sub-blocks are set to random
floating point values.
"""
function RandomMPS(k, n)
#=     if !isinteger(k)
        error("$(typeof(k)) is not accepted as KSites attribute. It should be a positive integer.")
        return nothing
    end

    if !isinteger(n)
        error("$(typeof(n)) is not accepted as NParticles attribute. It should be a positive integer.")
        return nothing
    end
 =#


    K = convert(UInt16, k)
    N = convert(UInt16, n)
    @assert n<=k
    X_ = Array{Tuple{DiagonalMPS, DiagonalMPS}}(undef, K)
    sizesRow = Array{UInt16}(undef, N+1)
    sizesCol = Array{UInt16}(undef, N+1)
    ran = 2 : 2^6

    #sizesRow[1], sizesCol[1] = 1, rand(ran)
    #sizesRow[N+1], sizesCol[N+1] = rand(ran), 1
    #for i = 2 : N
        #sizesRow[i], sizesCol[i] = rand(ran), rand(ran)
    #end
    sizesRow[1]=1
    for i = 1 : K


        indices1=GetIndices(i,N,K,1)
        indices2=GetIndices(i,N,K,2)
        N1=length(indices1)
        N2=length(indices2)
        sizes1 = Array{Tuple{UInt16, UInt16}}(undef, N1)
        sizes2 = Array{Tuple{UInt16, UInt16}}(undef, N2)

        for j = 1 :  max(N1,N2+1)
                 sizesCol[j] =  rand(ran)
        end



        l=1
        for j=indices1[1][1]:indices1[N1][1]
            if i==K
                sizes1[l] = (sizesRow[j], 1)
            else
                sizes1[l] = (sizesRow[j], sizesCol[j])
            end
            l+=1
        end

        l=2
        for j = indices2[1][1]+1 : indices2[N2][1]+1
            if i==K
                sizes2[l - 1] = (sizesRow[j - 1], 1)
            else
                sizes2[l - 1] = (sizesRow[j - 1], sizesCol[j])
            end
            l+=1
        end
        sizesRow=copy(sizesCol)




        diag1 = RandomDiagonalMPS(1, N1, sizes1,indices1)
        diag2 = RandomDiagonalMPS(2, N2, sizes2,indices2)
        X_[i] = (diag1, diag2)

    end

    return MPS(K, N, X_)
end
"""
DiagonalQR!(diagMPS, 1)
Computes the QR factorization over the first sub-block, overwriting the X[1]
values with the Q[1] orthogonal matrix values and returning the R[1] upper
triangular matrix.
"""
function MPSDiagonalQR!(diagMPS::DiagonalMPS, i)
#=     if !isinteger(i)
        error("$(typeof(i)) is not accepted as index. It should be a positive integer.")
        return nothing
    end
 =#

    index = convert(UInt16, i)

    F = qr(diagMPS.X[index])
    diagMPS.X[index] = Matrix(F.Q)

    return Matrix(F.R)
end

"""
DiagonalRXMultiplication!(diagMPS, R, 2)
Multiplies the provided R upper triangular matrix with the X[2], overwriting the X[2]
values with the result of the multiplication.
"""
function MPSDiagonalRXMultiplication!(diagMPS::DiagonalMPS, R, i)
#=     if !isinteger(i)
        error("$(typeof(i)) is not accepted as index. It should be a positive integer.")
        return nothing
    end
 =#

    R_ = convert(Matrix{Float64}, R)
    index = convert(UInt16, i)

    diagMPS.X[index] = R_ * diagMPS.X[index]

    nothing
end

"""
TT_MPS(mps::MPS,k)
returns the tensor train format of the mps for a given number of sites
"""
function TT_MPS(mps::MPS,k)
    K = convert(UInt16, k)
    tt_mps::Array{AbstractArray{Float64}, 1} = []
    for i=1:K
        push!(tt_mps,TT_core(mps.X[i]))
    end
    return tt_mps
end

"""
TT_core(mps.X[k])
returns the complete tensor core of the site k
"""
function TT_core(core::Tuple{DiagonalMPS,DiagonalMPS})
    Sizes::Tuple{Array{UInt16,1},Array{UInt16,1}}=GetSizes(core)
    r_left=Int.(sum(Sizes[1]))
    r_right=Int.(sum(Sizes[2]))
    Row_indices=sort([getfield.(core[1].BlockIndex,1);getfield.(core[2].BlockIndex,1)])
    Col_indices=sort([getfield.(core[1].BlockIndex,2); getfield.(core[2].BlockIndex,2)])
    gap_row=0;gap_col=0
    if minimum(Row_indices) != 1
        gap_row=Row_indices[1]-1
    end
    if minimum(Col_indices) != 1
        gap_col=Col_indices[1]-1
    end
    X::Array{Float64,3}=zeros(r_left,2,r_right)
    for j=1:2
        B=BlockArray(spzeros(r_left,r_right), Int64.(Sizes[1]), Int64.(Sizes[2]))
        B_sizes=core[j].BlockSizes
        for i=1:length(B_sizes)
            T=(core[j].BlockIndex[i][1]-gap_row,core[j].BlockIndex[i][2]-gap_col)
            B[Block(T)]=GetDiagonalBlock(core[j], i)
        end
        X[:,j,:]=B;
    end
    return X
end

"""
GetSizes(mps.X[k])
returns the Row and column sizes which are the ranks (\rho_k,0,\rho_k,1,...,\rho_k,N)
with k is the number of the local core and N is the number of particles
"""
function GetSizes(core::Tuple{DiagonalMPS,DiagonalMPS})::Tuple{Array{UInt16,1},Array{UInt16,1}}
    indices=[core[1].BlockIndex; core[2].BlockIndex]
    Sizes=[core[1].BlockSizes; core[2].BlockSizes]
    sortidx=sortperm(indices)
    Sizes=Sizes[sortidx]
    sort!(indices)
    id1= indexin(unique(getfield.(indices, 1)), getfield.(indices, 1))
    id2= indexin(unique(getfield.(indices, 2)), getfield.(indices, 2))
    RowSizes=getfield.(Sizes, 1)[id1]
    ColSizes=getfield.(Sizes, 2)[id2]
    return (RowSizes,ColSizes)
end

"""
PrintMPS(mps)

Prints on screen the formatted MPS structure.
"""
function PrintMPS(mps::MPS)
    @printf("KSites: [%d] - NParticles: [%d]\n", mps.KSites, mps.NParticles)

    for i = 1 : mps.KSites
        @printf("X[%d][1]: ", i)
        PrintDiagonalMPS(mps.X[i][1])
        @printf("X[%d][2]: ", i)
        PrintDiagonalMPS(mps.X[i][2])
        @printf("\n")
    end

    nothing
end
