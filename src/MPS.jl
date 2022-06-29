using Printf

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
    
    X_ = Array{Tuple{DiagonalMPS, DiagonalMPS}}(undef, K)
    sizesRow = Array{UInt16}(undef, N)
    sizesCol = Array{UInt16}(undef, N)
    ran = 2 : 2^3

    for i = 1 : N
        sizesRow[i], sizesCol[i] = rand(ran), rand(ran)
    end

    for i = 1 : K
        sizes1 = Array{Tuple{UInt16, UInt16}}(undef, N)
        sizes2 = Array{Tuple{UInt16, UInt16}}(undef, N - 1)

        for j = 1 : N
            sizes1[j] = (sizesRow[j], sizesCol[j])
        end

        for j = 2 : N
            sizes2[j - 1] = (sizesRow[j - 1], sizesCol[j])
        end

        diag1 = RandomDiagonalMPS(1, N, sizes1)
        diag2 = RandomDiagonalMPS(2, N - 1, sizes2)
        X_[i] = (diag1, diag2)
    end

    return MPS(K, N, X_)
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
