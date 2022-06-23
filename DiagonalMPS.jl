"""
A DiagonalMPS structure represents a set of BlockMPS block X_{k}[L] into the MPS structure,
with diagonal number "Diagonal" and number of BlockMPS "NBlock", each of size (BlockHeight X BlockWidth).
"""
mutable struct DiagonalMPS
    Diagonal::Int64
    NBlock::Int64
    BlockHeight::Int64
    BlockWidth::Int64
    Offset::Int64
    X::Array{BlockMPS}

    function DiagonalMPS(d::Int64, nb::Int64, bh::Int64, bw::Int64, o::Int64)
        X = Array{BlockMPS}(undef, nb)
        for i = 1:nb
            X[i] = BlockMPS(i, i + Offset, bh, bw)
        end
    end
end
