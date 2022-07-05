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
    sizesRow = Array{UInt16}(undef, N+1)
    sizesCol = Array{UInt16}(undef, N+1)
    ran = 2 : 2^3

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

"""
function TT_core(core::Tuple{DiagonalMPS,DiagonalMPS})
    e0=[0;1]
    e1=[1;0]
    for i in [1,2]
        for coordinates in core[i].BlockIndex
            

        end
    end



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
