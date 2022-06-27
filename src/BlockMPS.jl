using LinearAlgebra

"""
A BlockMPS struct represents a block (or matrix) X_{k,N}[0] or X_{k,N}[1] into
the DiagonalMPS structure, with position (Row_K, Col_N) in the global matrix.
The block size is (Height X Width) and its elements are stored in matrix X.
"""
mutable struct BlockMPS
    Row_K::Integer
    Col_N::Integer
    Height::Integer
    Width::Integer
    X::Matrix{AbstractFloat}
    
    function BlockMPS(rk::Integer, cn::Integer, h::Integer, w::Integer)
        new(rk, cn, h, w, zeros(h, w))
    end
end
