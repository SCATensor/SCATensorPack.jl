using LinearAlgebra
using Printf

"""
Abstract type for a DiagonalMPO structure.
"""
abstract type AbstractDiagonalMPO end

"""
Enum to set diagonal numbers as 1, 2, 3 or 4.
"""
@enum DNumber First=1 Second=2 Third=3 Fourth=4

"""
A DiagonalMPO structure represents a set of Matrix W_{k}[00],  W_{k}[01], W_{k}[10]
and W_{k}[11] into the MPO diagonal structure, each of size (BlockSizes[i][1] X BlockSizes[i][2])
and position indices (BlockIndex[i][1] X BlockIndex[i][2]). Diagonal number (1, 2, 3, 4)
is set in "DiagNumber" attribute and "NBlocks" is the total number of sub-blocks
contained in the DiagonalMPO. Block elements are stored in the matrix array W.
"""
mutable struct DiagonalMPO <: AbstractDiagonalMPO
    DiagNumber::UInt8
    NBlocks::UInt16
    W::Array{Matrix{Float64}}
    BlockIndex::Array{Tuple{UInt16, UInt16}}
    BlockSizes::Array{Tuple{UInt16, UInt16}}

    """
    diagMPO = DiagonalMPO(1, 2, [[1.0 2.0; 3.0 4.0], [1.0 2.0 3.0; 4.0 5.0 6.0]], [(1, 1), (2, 2)], [(2, 2), (2, 3)])

    Constructs a new DiagonalMPO structure with diagonal number 1 and 2 sub-blocks:
    the first one of size 2 X 2, the second one of size 2 X 3.
    Elements in the sub-blocks are set to the matrices provided as third argument
    and indices are set to the array provided as fourth argument.
    """
    function DiagonalMPS(diagN, nblocks, w, indices, sizes)
        dNumber = convert(UInt8, diagN)
        dNumber = (dNumber == UInt8(First) || dNumber == UInt8(Second) || dNumber == UInt8(Third) || dNumber == UInt8(Fourth)) ? dNumber :
                  (dNumber < UInt8(First)) ? UInt8(First) : UInt8(Fourth)
        new(dNumber, nblocks, w, indices, sizes)
    end
end

"""
diagMPO = DiagonalMPO(1, 2, [(2, 2), (2, 3)])

Constructs a new DiagonalMPO structure with diagonal number 1 and 2 sub-blocks:
the first one of size 2 X 2, the second one of size 2 X 3.
Elements in the sub-blocks are set to 0 and indices in BlockIndex array are set
according to diagonal number.
"""
function DiagonalMPO(diagN, nblocks, sizes)
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
    dNumber = (dNumber == UInt8(First) || dNumber == UInt8(Second) || dNumber == UInt8(Third) || dNumber == UInt8(Fourth)) ? dNumber :
              (dNumber < UInt8(First)) ? UInt8(First) : UInt8(Fourth)
    nBlocks = convert(UInt16, nblocks)
    bSizes = convert(Array{Tuple{UInt16, UInt16}}, sizes)

    W_ = Array{Matrix{Float64}}(undef, nBlocks)
    indices = Array{Tuple{UInt16, UInt16}, 1}(undef, nBlocks)
    for i = 1 : nBlocks
        indices[i] = (dNumber == UInt8(Second) ? i + 1 : i, dNumber == UInt8(Third) ? i + 1 : i)
        W_[i] = Matrix(zeros(bSizes[i][1], bSizes[i][2]))
    end

    return DiagonalMPO(dNumber, nBlocks, W_, indices, bSizes)
end

"""
diagMPO = DiagonalMPO(1, 2, [[1.0 2.0; 3.0 4.0], [1.0 2.0 3.0; 4.0 5.0 6.0]], [(2, 2), (2, 3)])

Constructs a new DiagonalMPO structure with diagonal number 1 and 2 sub-blocks:
the first one of size 2 X 2, the second one of size 2 X 3.
Elements in the sub-blocks are set to the matrices provided as third argument
and indices in the BlockIndex array are set according to diagonal number.
"""
function DiagonalMPO(diagN, nblocks, w, sizes)
#=     if !isinteger(diagN)
        error("$(typeof(diagN)) is not accepted as DiagNumber attribute. It should be a positive integer.")
        return nothing
    end

    if !isinteger(nblocks)
        error("$(typeof(nblocks)) is not accepted as NBlocks attribute. It should be a positive integer.")
        return nothing
    end

    if !isa(Array{Matrix{AbstractFloat}}, w)
        error("$(typeof(w)) is not accepted as W attribute. It should be an array of floating point matrix.")
        return nothing
    end

    if !isa(Array{Tuple{UInt, UInt}}, sizes)
        error("$(typeof(sizes)) is not accepted as BlockSizes attribute. It should be an array of positive integer tuples.")
        return nothing
    end
 =#
    dNumber = convert(UInt8, diagN)
    dNumber = (dNumber == UInt8(First) || dNumber == UInt8(Second) || dNumber == UInt8(Third) || dNumber == UInt8(Fourth)) ? dNumber :
              (dNumber < UInt8(First)) ? UInt8(First) : UInt8(Fourth)
    nBlocks = convert(UInt16, nblocks)
    W_ = convert(Array{Matrix{Float64}}, w)
    bSizes = convert(Array{Tuple{UInt16, UInt16}}, sizes)

    indices = Array{Tuple{UInt16, UInt16}, 1}(undef, nBlocks)
    for i = 1 : nBlocks
        indices[i] = (dNumber == UInt8(Second) ? i + 1 : i, dNumber == UInt8(Third) ? i + 1 : i)
    end

    return DiagonalMPO(dNumber, nBlocks, W_, indices, bSizes)
end

"""
diagMPO = RandomDiagonalMPO(1, 2, [(2, 2), (2, 3)])

Constructs a new DiagonalMPO structure with diagonal number 1 and 2 sub-blocks:
the first one of size 2 X 2, the second one of size 2 X 3.
Elements in the sub-blocks are set to random floating point values and
indices in the BlockIndex array are set according to diagonal number.
"""
function RandomDiagonalMPO(diagN, nblocks, sizes)
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
    dNumber = (dNumber == UInt8(First) || dNumber == UInt8(Second) || dNumber == UInt8(Third) || dNumber == UInt8(Fourth)) ? dNumber :
              (dNumber < UInt8(First)) ? UInt8(First) : UInt8(Fourth)
    nBlocks = convert(UInt16, nblocks)
    bSizes = convert(Array{Tuple{UInt16, UInt16}}, sizes)

    W_ = Array{Matrix{Float64}}(undef, nBlocks)
    indices = Array{Tuple{UInt16, UInt16}, 1}(undef, nBlocks)

    for i = 1 : nBlocks
        indices[i] = (dNumber == UInt8(Second) ? i + 1 : i, dNumber == UInt8(Third) ? i + 1 : i)
        W_[i] = rand(bSizes[i][1], bSizes[i][2])
    end

    return DiagonalMPO(dNumber, nBlocks, W_, indices, bSizes)
end

"""
diagMPO = RandomDiagonalMPO(1, 2, [(2, 2), (2, 3)], [(1, 1), (2, 2)])

Constructs a new DiagonalMPO structure with diagonal number 1 and 2 sub-blocks:
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
    dNumber = (dNumber == UInt8(First) || dNumber == UInt8(Second) || dNumber == UInt8(Third) || dNumber == UInt8(Fourth)) ? dNumber :
              (dNumber < UInt8(First)) ? UInt8(First) : UInt8(Fourth)
    nBlocks = convert(UInt16, nblocks)
    bSizes = convert(Array{Tuple{UInt16, UInt16}}, sizes)
    indices = convert(Array{Tuple{UInt16, UInt16}, 1}, index)

    W_ = Array{Matrix{Float64}}(undef, nBlocks)

    for i = 1 : nBlocks
        W_[i] = rand(bSizes[i][1], bSizes[i][2])
    end

    return DiagonalMPS(dNumber, nBlocks, W_, indices, bSizes)
end

"""
SetDiagonalBlock!(diagMPO, 1, (3, 3), [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0])

Sets the size of the first sub-block to 3 X 3 and its elements to the ones provided
as third argument.
"""
function SetDiagonalBlock!(diagMPO::DiagonalMPO, i, size, NewW)
#=     if !isinteger(i)
        error("$(typeof(i)) is not accepted as index. It should be a positive integer.")
        return nothing
    end

    if !isa(Tuple{UInt, UInt}, size)
        error("$(typeof(size)) is not accepted as new size. It should be a positive integer tuple.")
        return nothing
    end

    if !isa(Matrix{AbstractFloat}, NewW)
        error("$(typeof(NewW)) is not accepted as new matrix. It should be a floating point matrix.")
        return nothing
    end
 =#

    index = convert(UInt16, i)
    newSize = convert(Tuple{UInt16, UInt16}, size)
    newW = convert(Matrix{Float64}, NewW)

    diagMPO.W[index] = newW
    diagMPO.BlockSizes[index] = newSize

    nothing
end

"""
GetDiagonalBlock(diagMPO, 1)

Gets the first sub-block.
"""
function GetDiagonalBlock(diagMPO::DiagonalMPO, i)
#=     if !isinteger(i)
        error("$(typeof(i)) is not accepted as index. It should be a positive integer.")
        return nothing
    end
 =#

    index = convert(UInt16, i)

    return diagMPO.W[index]
end

"""
SetDiagonalBlockElement!(diagMPO, 1, (2, 2), 0.5)

Sets the element located at (2, 2) of the first sub-block with value 0.5.
"""
function SetDiagonalBlockElement!(diagMPO::DiagonalMPO, i, coords, value)
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

    diagMPO.W[index][c[1]][c[2]] = val

    nothing
end

"""
GetDiagonalBlockElement(diagMPO, 1, (2, 2))

Gets the element located at (2, 2) of the first sub-block.
"""
function GetDiagonalBlockElement(diagMPO::DiagonalMPO, i, coords)
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

    return diagMPO.W[index][c[1]][c[2]]
end

#= """
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
    diagN = (diagN == UInt8(First) || diagN == UInt8(Second) || diagN == UInt8(Third) || diagN == UInt8(Fourth)) ? diagN :
            (diagN < UInt8(First)) ? UInt8(First) : UInt8(Fourth)

    indices = Array{Tuple{UInt16, UInt16}}(undef)
    if i < (K - N + 1)
        if i > N
            #indices = [(l - (diagN - 1), l - (diagN - 1)) for l = diagN : N+1]
            indices = [(l - (diagN - 1), l ) for l = diagN : N+1]
        else
            #indices = [(l, l) for l = 1 : i]
            indices = [(l, l+(diagN-1)) for l = 1 : i]
        end
    else
        leq = i - (K - N) + 1 - (diagN - 1)
        geq = N - (diagN - 1) + 1
        if i <= N
#            indices = [(l, l) for l = leq : (N - i + 1) + 1]
            indices = [(l, l+(diagN-1)) for l = leq : (N - i + 1) + 1]
        else
#            indices = [(l, l) for l = leq : geq]
            indices = [(l, l+(diagN-1)) for l = leq : geq]
        end
    end

    return indices
end =#

"""
PrintDiagonalMPO(diagMPO)

Prints on screen the formatted DiagonalMPO structure.
"""
function PrintDiagonalMPO(diagMPO::DiagonalMPO)
    @printf("DiagNumber: [%d] - NBlocks: [%d]\n", diagMPO.DiagNumber, diagMPO.NBlocks)

    @printf("BlockIndex: [")
    for i = 1 : diagMPO.NBlocks
        @printf("(%d, %d), ", diagMPO.BlockIndex[i][1], diagMPO.BlockIndex[i][2])
    end
    @printf("]\t - \t")

    @printf("BlockSizes: [")
    for i = 1 : diagMPO.NBlocks
        @printf("(%d, %d), ", diagMPO.BlockSizes[i][1], diagMPO.BlockSizes[i][2])
    end
    @printf("]\n")

    for i = 1 : diagMPO.NBlocks
        @printf("W[%d]: ", i)
        display(diagMPO.W[i])
    end

    nothing
end
