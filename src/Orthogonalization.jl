"""Build Right orthogonalization and left orthogonalization of MPS with QR"""
using LinearAlgebra

"""
Right_orth(mps,K)
Computes the QR factorization from left to right of given mps structure  with K sites
"""

function Right_orth(mps::MPS,K)::MPS
    nsites = convert(UInt16, K)
    for i=1:nsites-1
         Right_orth_local!(mps.X[i],mps.X[i+1])
    end
end

"""
Right_orth_local(mps1,mps2)
Computes the QR factorization of the local tensor core mps1 and updates mps2
"""

function Right_orth_local!(mps1::Tuple{DiagonalMPS,DiagonalMPS},
    mps2::Tuple{DiagonalMPS,DiagonalMPS})
    Block::Array{Matrix{Float64}}
    N1=mps1[1].NBlocks
    N2=mps1[2].NBlocks

    #Elements present in 1 and not present in 2
    L1=mps1[1].BlockIndex
    L2=[(i[1]+1,i[1]+1) for i in mps1[2].BlockIndex ]

    L1_neg=findall(!in(L1),L2)
    if !isempty(L1_neg)
        for i=1:length(L1_neg)
            F=qr(mps1[1].X[L1_neg[i]])
            mps1[1].X[L1_neg[i]]= Matrix(F.Q)
            diagonalRXMultiplication!(mps2, R)
        end
    end
    L2_neg=findall(!in(L2),L1)
    if !isempty(L2_neg)
        for i=1:length(L2_neg)
            F=qr(mps1[1].X[L2_neg[i]])
            mps1[1].X[L2_neg[i]]= Matrix(F.Q)
            diagonalRXMultiplication!(mps2, R)
        end
    end
    L1_pos=findall(!in(L1),L2)
    L2_pos=findall(!in(L2),L1)
    for i=1:length(L1_pos)
        F=qr([mps1[1].X[L2_pos[i]];mps1[1].X[L1_pos[i]] ])
        Q= Matrix(F.Q)
        
        diagonalRXMultiplication!(mps2, R)
    end

end
