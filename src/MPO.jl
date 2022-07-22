using Printf
using BlockArrays
using SparseArrays

"""
Abstract type for an MPO structure.
"""
abstract type AbstractMPO end

"""
An MPO structure represents a set of four DiagonalMPS W_{00}, W_{11},
W_{01} and W_{10} structures. Block elements are stored in the array X.
"""
mutable struct MPO <: AbstractMPO
    KSites::UInt16
    MParticles::UInt16
    IDs::Array{Array{Tuple{UInt8, UInt8}}}
    W::Array{Array{DiagonalMPS, DiagonalMPS, DiagonalMPS, DiagonalMPS}}

    """
    mpo = MPO(2, 2, [[diag00MPO1, diag01MPO1, diag10MPO1, diag11MPO1], [diag00MPO2, diag01MPO2, diag10MPO2, diag11MPO2]])

    Constructs a new MPO structure with 2 sites and 2 particles. Each of the
    DiagonalMPS passed as third argument should be previously created.
    """
    function MPO(k, m, w)
        K = convert(UInt16, k)
        M = convert(UInt16, m)

        @assert (K >= M) "KSites should be greater or equal to MParticles"

        IDs{Array{Tuple{UInt8, UInt8}}}(undef, K)
        for i = 1 : K
            IDs[i]{Array{Tuple{UInt8, UInt8}}}(undef, 4)
            IDs[i][1] = (0, 0)
            IDs[i][2] = (0, 1)
            IDs[i][3] = (1, 0)
            IDs[i][4] = (1, 1)
        end

        new(K, M, w)
    end
end

"""
mpo = MPO(2, 2)

Constructs a new MPO structure with 2 sites and 2 particles. Each of the
DiagonalMPS contained in the MPO are set to "undef" values.
"""
function MPO(k, m)
#=     if !isinteger(k)
        error("$(typeof(k)) is not accepted as KSites attribute. It should be a positive integer.")
        return nothing
    end

    if !isinteger(m)
        error("$(typeof(m)) is not accepted as MParticles attribute. It should be a positive integer.")
        return nothing
    end
 =#
    K = convert(UInt16, k)
    M = convert(UInt16, m)

    @assert (K >= M) "KSites should be greater or equal to MParticles"

    W_ = Array{Array{Tuple{DiagonalMPS, DiagonalMPS}}}(undef, K)

    return MPO(K, N, W_)
end

"""
mpo = RandomMPO(2, 2)

Constructs a new MPO structure with 2 sites and 2 particles. Each of the
DiagonalMPS contained in the MPO are generated randomly; this is, with random
sizes in the range (2 : 2^3) and elements in the sub-blocks are set to random
floating point values.
"""
function RandomMPO(k, m)
#=     if !isinteger(k)
        error("$(typeof(k)) is not accepted as KSites attribute. It should be a positive integer.")
        return nothing
    end

    if !isinteger(m)
        error("$(typeof(m)) is not accepted as MParticles attribute. It should be a positive integer.")
        return nothing
    end
 =#
    K = convert(UInt16, k)
    M = convert(UInt16, m)

    @assert (K >= M) "KSites should be greater or equal to MParticles"

    W_ = Array{Array{Tuple{DiagonalMPS, DiagonalMPS}}}(undef, K)
    sizesRow = Array{UInt16}(undef, M + 1)
    sizesCol = Array{UInt16}(undef, M + 1)
    ran = 2 : 2^3

    sizesRow[1] = 1
    for i = 1 : K
        indices1 = GetIndices(i, M, K, 1)
        indices2 = GetIndices(i, M, K, 2)
        N1 = length(indices1)
        N2 = length(indices2)
        sizes1 = Array{Tuple{UInt16, UInt16}}(undef, N1)
        sizes2 = Array{Tuple{UInt16, UInt16}}(undef, N2)

        for j = 1 :  max(N1, N2 + 1)
            sizesCol[j] = rand(ran)
        end

        l = 1
        for j = indices1[1][1] : indices1[N1][1]
            sizes1[l] = (i == K) ? (sizesRow[j], 1) : (sizesRow[j], sizesCol[j])
            l += 1
        end

        l = 2
        for j = indices2[1][1] + 1 : indices2[N2][1] + 1
            sizes2[l - 1] = (i == K) ? (sizesRow[j - 1], 1) : (sizesRow[j - 1], sizesCol[j])
            l += 1
        end

        sizesRow = copy(sizesCol)
        diag1 = RandomDiagonalMPO(1, N1, sizes1, indices1)
        diag2 = RandomDiagonalMPO(2, N2, sizes2, indices2)
        W_[i] = (diag1, diag2)
    end

    return MPO(K, M, W_)
end

#= """
MPS_QR!(MPS, 1)

Computes the QR factorization over the first block concatenation, overwriting the X[1]
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
end =#

"""
concatDiagonalMPO = ConcatDiagonalMPO(mpo, 1)

Concatenates all sub-blocks from W[1]/W[2]/W[3]/W[4] DiagonalMPO structures by column
in site 1. This is, all sub-blocks belonging to the same column will be concatenated if possible,
putting the sub-block on W[1] on top of the sub-block on W[2], W[2] on top of the sub-block on W[3],
W[3] on top of the sub-block on W[4]. For example, with an MPO as:

W[1][1,1]       0           0           W[4][1,1]       0           0
    0       W[1][2,2]       0               0       W[4][2,2]       0
    0           0       W[1][3,3]           0           0       W[4][3,3]

    0           0           0               0       W[3][1,2]       0
W[2][2,1]       0           0               0           0       W[3][2,3]
    0       W[2][3,2]       0               0           0           0

The resulting concatenations will take place:

W[1][1,1]   W[1][2,2]   W[1][3,3]
W[2][2,1]   W[2][3,2]       0
    0       W[3][1,2]   W[3][2,3]
W[4][1,1]   W[4][2,2]   W[4][3,3]

If there are not sub-blocks on the column, it will put a zero matrix. After all possible
concatenations have finished, it will return the array of matrices.
"""
function ConcatDiagonalMPO(mpo::MPO, k)
#=     if !isinteger(k)
        error("$(typeof(k)) is not accepted as site number. It should be a positive integer.")
        return nothing
    end
 =#

    K = convert(UInt16, k)
    maxi = max(maximum(mpo.W[K][1].BlockIndex), maximum(mpo.W[K][2].BlockIndex))[2]
    con = Array{Matrix{Float64}}(undef, maxi)

    for i = 1 : maxi
        index1 = findall(t -> t[2] == i, mps.X[K][1].BlockIndex)
        index2 = findall(t -> t[2] == i, mps.X[K][2].BlockIndex)
        con[i] = (!isempty(index1) && !isempty(index2)) ? vcat(mps.X[K][2].X[index2[1]], mps.X[K][1].X[index1[1]]) :
                 (!isempty(index1) && isempty(index2)) ? mps.X[K][1].X[index1[1]] :
                 (isempty(index1) && !isempty(index2)) ? mps.X[K][2].X[index2[1]] :
                 Matrix(undef, 0, 0);
    end

    return con
end

"""
concatAllDiagonalMPS = ConcatAllDiagonalMPS(mps)

Concatenates all sub-blocks from all X[K][1]/X[K][2] DiagonalMPS structures by column
in all K sites. This is, for every site, all sub-blocks belonging to the same column
will be concatenated if possible, putting the sub-block on X[K][2] on top of the sub-block
on X[K][1]. For example:

X_1[1][1,1]   X_1[2][1,2]
    |           |
Single      Single
block       block

            X_2[2][1,2]
            X_2[1][2,2] X_2[2][2,3]
    |           |           |
Empty       Concat of   Single
block       two blocks  block

If there are not sub-blocks on the column, it will put an empty matrix. After all possible
concatenations have finished, it will return the array of concatenated matrices.
"""
function ConcatAllDiagonalMPS(mps::MPS)
    allCon = Array{Array{Matrix{Float64}}}(undef, mps.KSites)
    for i = 1 : mps.KSites
        allCon[i] = ConcatDiagonalMPS(mps, i)
    end

    return allCon
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
