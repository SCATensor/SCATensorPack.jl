"""
A BlockMPS represents a block X_{k,l}[L] into the DiagonalMPS structure,
with position (Row, Col) in the global matrix of size Height X Width.
"""
mutable struct BlockMPS
    Row::Int64
    Col::Int64
    Height::Int64
    Width::Int64
    X::Array{Float64}
    
    function BlockMPS(r::Int64, c::Int64, h::Int64, w::Int64)
        new(r, c, h, w, Array{Float64}(undef, h, w))
    end
end
