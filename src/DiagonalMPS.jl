using LinearAlgebra
using Printf

abstract type AbstractDiagonalMPS end

"""
Enum to set diagonal numbers as 1 or 2.
"""
@enum DNumber First=1 Second=2

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
    diagMPS = DiagonalMPS(1, 2, [[1.0 2.0; 3.0 4.0], [1.0 2.0 3.0; 4.0 5.0 6.0]], [(1, 1), (2, 2)], [(2, 2), (2, 3)])

    Constructs a new DiagonalMPS structure with diagonal number 1 and 2 sub-blocks:
    the first one of size 2 X 2, the second one of size 2 X 3.
    Elements in the sub-blocks are set to the matrices provided as third argument
    and indices are set to the array provided as fourth argument.
    """
    function DiagonalMPS(diagN, nblocks, x, indices, sizes)
        dNumber = convert(UInt8, diagN)
        dNumber = (dNumber == UInt8(First) || dNumber == UInt8(Second)) ? dNumber : (dNumber < UInt8(First)) ? UInt8(First) : UInt8(Second)
        new(dNumber, nblocks, x, indices, sizes)
    end
end

"""
diagMPS = DiagonalMPS(1, 2, [(2, 2), (2, 3)])

Constructs a new DiagonalMPS structure with diagonal number 1 and 2 sub-blocks:
the first one of size 2 X 2, the second one of size 2 X 3.
Elements in the sub-blocks are set to 0 and indices in BlockIndex array are set
to (i, i) if diagonal number is 1 (i.e. (1, 1), (2, 2), ...) or (i, i+1) otherwise
(i.e. (1, 2), (2, 3), ...).
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
    dNumber = (dNumber == UInt8(First) || dNumber == UInt8(Second)) ? dNumber : (dNumber < UInt8(First)) ? UInt8(First) : UInt8(Second)
    nBlocks = convert(UInt16, nblocks)
    bSizes = convert(Array{Tuple{UInt16, UInt16}}, sizes)

    X_ = Array{Matrix{Float64}}(undef, nBlocks)
    indices = Array{Tuple{UInt16, UInt16}, 1}(undef, nBlocks)
    for i = 1 : nBlocks
        indices[i] = (i, dNumber == UInt8(First) ? i : i + 1)
        X_[i] = Matrix(zeros(bSizes[i][1], bSizes[i][2]))
    end

    return DiagonalMPS(dNumber, nBlocks, X_, indices, bSizes)
end

"""
diagMPS = DiagonalMPS(1, 2, [[1.0 2.0; 3.0 4.0], [1.0 2.0 3.0; 4.0 5.0 6.0]], [(2, 2), (2, 3)])

Constructs a new DiagonalMPS structure with diagonal number 1 and 2 sub-blocks:
the first one of size 2 X 2, the second one of size 2 X 3.
Elements in the sub-blocks are set to the matrices provided as third argument
and indices in the BlockIndex array are set to (i, i) if diagonal number is 1
(i.e. (1, 1), (2, 2), ...) or (i, i+1) otherwise (i.e. (1, 2), (2, 3), ...).
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
    dNumber = (dNumber == UInt8(First) || dNumber == UInt8(Second)) ? dNumber : (dNumber < UInt8(First)) ? UInt8(First) : UInt8(Second)
    nBlocks = convert(UInt16, nblocks)
    X_ = convert(Array{Matrix{Float64}}, x)
    bSizes = convert(Array{Tuple{UInt16, UInt16}}, sizes)

    indices = Array{Tuple{UInt16, UInt16}, 1}(undef, nBlocks)
    for i = 1 : nBlocks
        indices[i] = (i, dNumber == UInt8(First) ? i : i + 1)
    end

    return DiagonalMPS(dNumber, nBlocks, X_, indices, bSizes)
end

"""
diagMPS = RandomDiagonalMPS(1, 2, [(2, 2), (2, 3)])

Constructs a new DiagonalMPS structure with diagonal number 1 and 2 sub-blocks:
the first one of size 2 X 2, the second one of size 2 X 3.
Elements in the sub-blocks are set to random floating point values and
indices in the BlockIndex array are set to (i, i) if diagonal number is 1
(i.e. (1, 1), (2, 2), ...) or (i, i+1) otherwise (i.e. (1, 2), (2, 3), ...).
"""
function RandomDiagonalMPS(diagN, nblocks, sizes)
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
    dNumber = (dNumber == UInt8(First) || dNumber == UInt8(Second)) ? dNumber : (dNumber < UInt8(First)) ? UInt8(First) : UInt8(Second)
    nBlocks = convert(UInt16, nblocks)
    bSizes = convert(Array{Tuple{UInt16, UInt16}}, sizes)

    X_ = Array{Matrix{Float64}}(undef, nBlocks)
    indices = Array{Tuple{UInt16, UInt16}, 1}(undef, nBlocks)

    for i = 1 : nBlocks
        indices[i] = (i, dNumber == UInt8(First) ? i : i + 1)
        X_[i] = rand(bSizes[i][1], bSizes[i][2])
    end

    return DiagonalMPS(dNumber, nBlocks, X_, indices, bSizes)
end

"""
diagMPS = RandomDiagonalMPS(1, 2, [(2, 2), (2, 3)], [(1, 1), (2, 2)])

Constructs a new DiagonalMPS structure with diagonal number 1 and 2 sub-blocks:
the first one of size 2 X 2, the second one of size 2 X 3.
Elements in the sub-blocks are set to random floating point values and
indices in the BlockIndex array are set to the indices specified as fourth argument.
"""
function RandomDiagonalMPS(diagN, nblocks, sizes, index)
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

    if !isa(Array{Tuple{UInt, UInt}}, index)
        error("$(typeof(index)) is not accepted as BlockIndex attribute. It should be an array of positive integer tuples.")
        return nothing
    end
 =#
    dNumber = convert(UInt8, diagN)
    dNumber = (dNumber == UInt8(First) || dNumber == UInt8(Second)) ? dNumber : (dNumber < UInt8(First)) ? UInt8(First) : UInt8(Second)
    nBlocks = convert(UInt16, nblocks)
    bSizes = convert(Array{Tuple{UInt16, UInt16}}, sizes)
    indices = convert(Array{Tuple{UInt16, UInt16}, 1}, index)

    X_ = Array{Matrix{Float64}}(undef, nBlocks)

    for i = 1 : nBlocks
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
Author: Siwar Badreddine
indices = GetIndices(1, 3, 6, 1)

Gives back the sub-block indices of the first tensor core X_1[1] for 3 particles and 6 sites
with respect to the number of particles.
"""
function GetIndices(i, N, K, diagN)::Array{Tuple{UInt16, UInt16}}
#=     if !isinteger(i)
        error("$(typeof(i)) is not accepted as index. It should be a positive integer.")
        return nothing
    end

    if !isinteger(N)
        error("$(typeof(N)) is not accepted as particle number. It should be a positive integer.")
        return nothing
    end

    if !isinteger(K)
        error("$(typeof(K)) is not accepted as site number. It should be a positive integer.")
        return nothing
    end

    if !isinteger(diagN)
        error("$(typeof(diagN)) is not accepted as diagonal number. It should be a positive integer.")
        return nothing
    end
 =#

    i  = convert(UInt16, i)
    N  = convert(UInt16, N)
    K  = convert(UInt16, K)
    diagN  = convert(UInt8, diagN)
    diagN = (diagN == UInt8(First) || diagN == UInt8(Second)) ? diagN : (diagN < UInt8(First)) ? UInt8(First) : UInt8(Second)
    #=if diagN==2
        diagN+=1
    end=#

    indices = Array{Tuple{UInt16, UInt16}}(undef)
    if i < (K - N + 1)
        if i > N
            indices = [(l - (diagN - 1), l - (diagN - 1)) for l = diagN : N+1]
        else
            indices = [(l, l) for l = 1 : i]
        end
    else
        leq = i - (K - N) + 1 - (diagN - 1)
        geq = N - (diagN - 1) + 1
        if i <= N
            indices = [(l, l) for l = leq : (N - i + 1) + 1]
        else
            indices = [(l, l) for l = leq : geq]
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
