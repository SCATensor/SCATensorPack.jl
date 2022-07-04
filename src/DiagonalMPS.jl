using LinearAlgebra
using Printf

abstract type AbstractDiagonalMPS end

"""
A DiagonalMPS structure represents a set of Matrix X_{k,N}[1] or X_{k,N}[2]
into the MPS diagonal structure, each of size (BlockSizes[i][1] X BlockSizes[i][2])
and position indices (BlockIndex[i][1] X BlockIndex[i][2]). Diagonal number (1 or 2)
is set in "DiagNumber" attribute and "NBlocks" is the total number of sub-blocks
contained in the DiagonalMPS. Block elements are stored in the matrix array X.
"""
mutable struct DiagonalMPS <: AbstractDiagonalMPS
    DiagNumber::UInt8
    NBlocks::UInt16
    X::Array{Matrix{Float64}}
    BlockIndex::Array{Tuple{UInt16, UInt16}}
    BlockSizes::Array{Tuple{UInt16, UInt16}}

    """
    diagMPS = DiagonalMPS(1, 2, [[1.0 2.0; 3.0 4.0], [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]], [(1, 1), (2, 2)], [(2, 2), (3, 3)])

    Constructs a new DiagonalMPS structure with diagonal number 1 and 2 sub-blocks:
    the first one of size 2 X 2, the second one of size 3 X 3.
    Elements in the sub-blocks are set to the matrices provided as third argument
    and indices are set to the array provided as fourth argument.
    """
    function DiagonalMPS(diagN, nblocks, x, indices, sizes)
        new(diagN, nblocks, x, indices, sizes)
    end
end

"""
diagMPS = DiagonalMPS(1, 2, [(2, 2), (3, 3)])

Constructs a new DiagonalMPS structure with diagonal number 1 and 2 sub-blocks:
the first one of size 2 X 2, the second one of size 3 X 3.
Elements in the sub-blocks are set to 0 and indices in BlockIndex array are set
to (i,i) (i.e. (1, 1), (2, 2), ...).
"""
function DiagonalMPS(diagN, nblocks, sizes)
#=     if !isinteger(diagN)
        error("$(typeof(diagN)) is not accepted as DiagNumber attribute. It should be a positive integer.")
        return nothing
    end

    if !isinteger(nblocks)
        error("$(typeof(nblocks)) is not accepted as NBlocks attribute. It should be a positive integer.")
        return nothing
    end

    if !isa(Array{Tuple{UInt, UInt}}, sizes)
        error("$(typeof(sizes)) is not accepted as BlockSizes attribute. It should be an array of positive integer tuples.")
        return nothing
    end
 =#
    dNumber = convert(UInt8, diagN)
    nBlocks = convert(UInt16, nblocks)
    bSizes = convert(Array{Tuple{UInt16, UInt16}}, sizes)

    X_ = Array{Matrix{Float64}}(undef, nBlocks)
    indices = Array{Tuple{UInt16, UInt16}}(undef, nBlocks)
    for i = 1 : nBlocks
        indices[i] = (i, i)
        X_[i] = Matrix(zeros(bSizes[i][1], bSizes[i][2]))
    end

    return DiagonalMPS(dNumber, nBlocks, X_, indices, bSizes)
end

"""
diagMPS = DiagonalMPS(1, 2, [[1.0 2.0; 3.0 4.0], [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]], [(2, 2), (3, 3)])

Constructs a new DiagonalMPS structure with diagonal number 1 and 2 sub-blocks:
the first one of size 2 X 2, the second one of size 3 X 3.
Elements in the sub-blocks are set to the matrices provided as third argument
and indices in the BlockIndex array are set to (i,i) (i.e. (1, 1), (2, 2), ...).
"""
function DiagonalMPS(diagN, nblocks, x, sizes)
#=     if !isinteger(diagN)
        error("$(typeof(diagN)) is not accepted as DiagNumber attribute. It should be a positive integer.")
        return nothing
    end

    if !isinteger(nblocks)
        error("$(typeof(nblocks)) is not accepted as NBlocks attribute. It should be a positive integer.")
        return nothing
    end

    if !isa(Array{Matrix{AbstractFloat}}, x)
        error("$(typeof(x)) is not accepted as X attribute. It should be an array of floating point matrix.")
        return nothing
    end

    if !isa(Array{Tuple{UInt, UInt}}, sizes)
        error("$(typeof(sizes)) is not accepted as BlockSizes attribute. It should be an array of positive integer tuples.")
        return nothing
    end
 =#
    dNumber = convert(UInt8, diagN)
    nBlocks = convert(UInt16, nblocks)
    X_ = convert(Array{Matrix{Float64}}, x)
    bSizes = convert(Array{Tuple{UInt16, UInt16}}, sizes)

    indices = Array{Tuple{UInt16, UInt16}}(undef, nBlocks)
    for i = 1 : nBlocks
        indices[i] = (i, i)
    end

    return DiagonalMPS(dNumber, nBlocks, X_, indices, bSizes)
end

"""
diagMPS = RandomDiagonalMPS(1, 2, [(3, 4), (2, 5)])

Constructs a new DiagonalMPS structure with diagonal number 1 and 2 sub-blocks:
the first one of size 3 X 4, the second one of size 2 X 5.
Elements in the sub-blocks are set to random floating point values.
"""
function RandomDiagonalMPS(diagN, nblocks, sizes,id)
#=     if !isinteger(diagN)
        error("$(typeof(diagN)) is not accepted as DiagNumber attribute. It should be a positive integer.")
        return nothing
    end

    if !isinteger(nblocks)
        error("$(typeof(nblocks)) is not accepted as NBlocks attribute. It should be a positive integer.")
        return nothing
    end

    if !isa(Array{Tuple{UInt, UInt}}, sizes)
        error("$(typeof(sizes)) is not accepted as BlockSizes attribute. It should be an array of positive integer tuples.")
        return nothing
    end
 =#
    dNumber = convert(UInt8, diagN)
    nBlocks = convert(UInt16, nblocks)
    bSizes = convert(Array{Tuple{UInt16, UInt16}}, sizes)
    indices=convert(Array{Tuple{Int64,Int64},1},id)

    X_ = Array{Matrix{Float64}}(undef, nBlocks)
    #indices = Array{Tuple{UInt16, UInt16}}(undef, nBlocks)

    for i = 1 : nBlocks
        #indices[i] = (i, i)
        X_[i] = rand(bSizes[i][1], bSizes[i][2])
    end

    return DiagonalMPS(dNumber, nBlocks, X_, indices, bSizes)
end

"""
SetDiagonalBlock!(diagMPS, 1, (3, 3), [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0])

Sets the size of the first sub-block to 3 X 3 and its elements to the ones provided
as third argument.
"""
function SetDiagonalBlock!(diagMPS::DiagonalMPS, i, size, NewX)
#=     if !isinteger(i)
        error("$(typeof(i)) is not accepted as index. It should be a positive integer.")
        return nothing
    end

    if !isa(Tuple{UInt, UInt}, size)
        error("$(typeof(size)) is not accepted as new size. It should be a positive integer tuple.")
        return nothing
    end

    if !isa(Matrix{AbstractFloat}, NewX)
        error("$(typeof(NewX)) is not accepted as new matrix. It should be a floating point matrix.")
        return nothing
    end
 =#

    index = convert(UInt16, i)
    newSize = convert(Tuple{UInt16, UInt16}, size)
    newX = convert(Matrix{Float64}, NewX)
    diagMPS.X[index] = newX
    diagMPS.BlockSizes[index] = newSize

    nothing
end

"""
GetDiagonalBlock(diagMPS, 1)

Gets the first sub-block.
"""
function GetDiagonalBlock(diagMPS::DiagonalMPS, i)
#=     if !isinteger(i)
        error("$(typeof(i)) is not accepted as index. It should be a positive integer.")
        return nothing
    end
 =#

    index = convert(UInt16, i)

    return diagMPS.X[index]
end

"""
SetDiagonalBlockElement!(diagMPS, 1, (2, 2), 0.5)

Sets the element located at (2, 2) of the first sub-block with value 0.5.
"""
function SetDiagonalBlockElement!(diagMPS::DiagonalMPS, i, coords, value)
#=     if !isinteger(i)
        error("$(typeof(i)) is not accepted as index. It should be a positive integer.")
        return nothing
    end

    if !isa(Tuple{UInt, UInt}, coords)
        error("$(typeof(coords)) is not accepted as coordinates. It should be a positive integer tuple.")
        return nothing
    end

    if !isa(AbstractFloat, value)
        error("$(typeof(value)) is not accepted as value. It should be a float value.")
        return nothing
    end
 =#

    index = convert(UInt16, i)
    c = convert(Tuple{UInt16, UInt16}, coords)
    val = convert(Float64, value)
    diagMPS.X[index][c[1]][c[2]] = val

    nothing
end

"""
GetDiagonalBlockElement(diagMPS, 1, (2, 2))

Gets the element located at (2, 2) of the first sub-block.
"""
function GetDiagonalBlockElement(diagMPS::DiagonalMPS, i, coords)
#=     if !isinteger(i)
        error("$(typeof(i)) is not accepted as index. It should be a positive integer.")
        return nothing
    end

    if !isa(Tuple{UInt, UInt}, coords)
        error("$(typeof(coords)) is not accepted as coordinates. It should be a positive integer tuple.")
        return nothing
    end
 =#

    index = convert(UInt16, i)
    c = convert(Tuple{UInt16, UInt16}, coords)

    return diagMPS.X[index][c[1]][c[2]]
end

"""
DiagonalQR!(diagMPS, 1)

Computes the QR factorization over the first sub-block, overwriting the X[1]
values with the Q[1] orthogonal matrix values and returning the R[1] upper
triangular matrix.
"""
function DiagonalQR!(diagMPS::DiagonalMPS, i)
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
function DiagonalRXMultiplication!(diagMPS::DiagonalMPS, R, i)
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
Author: Siwar Badreddine
GetIndices(1,3,6,1)

Gives back the sub-block indices of the first tensor core X_1[0] for 3 particles and 6 sites
with respect to the number of particles
"""
function GetIndices(i,N,K,diagN)::Array{Tuple{Int64,Int64},1}


    i  = convert(UInt16, i)
    N  = convert(UInt16, N)
    K  = convert(UInt16, K)
    diagN  = convert(UInt16, diagN)
    indices = Array{Tuple{UInt16, UInt16}}(undef)
    if i < (K-N+1)
        if i> N
            indices=[(l-(diagN-1),l-(diagN-1)) for l=diagN:N+1]
        else
            indices=[(l,l) for l=1:i]
        end
    else
        leq=i-(K-N)+1-(diagN-1)
        geq=N-(diagN-1)+1
        if i<=N
            indices=[(l,l) for l=leq:(N-i+1)+1]
        else
            indices=[(l,l) for l=leq:geq]
        end
    end

    return indices
end

"""
PrintDiagonalMPS(diagMPS)

Prints on screen the formatted DiagonalMPS structure.
"""
function PrintDiagonalMPS(diagMPS::DiagonalMPS)
    @printf("DiagNumber: [%d] - NBlocks: [%d]\n", diagMPS.DiagNumber, diagMPS.NBlocks)

    @printf("BlockIndex: [")
    for i = 1 : diagMPS.NBlocks
        @printf("(%d, %d), ", diagMPS.BlockIndex[i][1], diagMPS.BlockIndex[i][2])
    end
    @printf("]\t - \t")

    @printf("BlockSizes: [")
    for i = 1 : diagMPS.NBlocks
        @printf("(%d, %d), ", diagMPS.BlockSizes[i][1], diagMPS.BlockSizes[i][2])
    end
    @printf("]\n")

    for i = 1 : diagMPS.NBlocks
        @printf("X[%d]: ", i)
        display(diagMPS.X[i])
    end

    nothing
end
