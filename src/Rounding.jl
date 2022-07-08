#Compressing MPS cores using truncated_SVD left and right
using LinearAlgebra
"""
Left_rounding(mps,K;δ)
Computes the SVD decomposition  from left to right of given mps structure  with K sites
δ is the truncation error
"""
function Left_rounding!(mps::MPS,K;δ::Float64=1e-12)
	nsites = convert(UInt16, K)
	for i=1:nsites-1
		Left_rounding_local!(mps.X[i],mps.X[i+1],δ=δ)
	end
	nothing
end


"""
Left_orth_local(mps1,mps2)
Computes the SVD decomposition of the local tensor core mps1 and updates mps2
"""
function Left_rounding_local!(mps1::Tuple{DiagonalMPS,DiagonalMPS},
	mps2::Tuple{DiagonalMPS,DiagonalMPS};δ=1e-12)
	L1=mps1[1].BlockIndex
	L2=mps1[2].BlockIndex
	for j=1:length(L1)
		if getfield.(L1, 2)[j] in getfield.(L2, 2)
			id2=findall(isequal(getfield.(L1, 2)[j]), getfield.(L2, 2))[1]
			Block=[GetDiagonalBlock(mps1[2],id2);GetDiagonalBlock(mps1[1], j)]
			U,SV=truncated_svd(Block,δ=δ)
			Update_left!(mps1,U,[id2,j],concat=true)
			Left_DiagonalRXMultiplication!(mps2, SV,getfield.(L1, 2)[j])
		else
			Block=GetDiagonalBlock(mps1[1], j)
			U,SV=truncated_svd(Block,δ=δ)
			Update_left!(mps1,U,[j],concat=false,flag=1)
			Left_DiagonalRXMultiplication!(mps2, SV,getfield.(L1, 2)[j])
		end
	end
	for j=1:length(L2)
		if getfield.(L2, 2)[j]   ∉ getfield.(L1, 2)
			Block=GetDiagonalBlock(mps1[2], j)
			U,SV=truncated_svd(Block,δ=δ)
			Update_left!(mps1,U,[j],concat=false,flag=2)
			Left_DiagonalRXMultiplication!(mps2,SV,getfield.(L2, 2)[j])
		end
	end
	nothing
end

"""
Left_rounding(mps,K;δ)
Computes the SVD decomposition  from left to right of given mps structure  with K sites
δ is the truncation error
"""
function Left_rounding!(mps::MPS,K;δ::Float64=1e-12)
	nsites = convert(UInt16, K)
	for i=1:nsites-1
		Left_rounding_local!(mps.X[i],mps.X[i+1],δ=δ)
	end
	nothing
end


"""
Left_orth_local(mps1,mps2)
Computes the SVD decomposition of the local tensor core mps1 and updates mps2
"""
function Left_rounding_local!(mps1::Tuple{DiagonalMPS,DiagonalMPS},
	mps2::Tuple{DiagonalMPS,DiagonalMPS};δ=1e-12)
	L1=mps1[1].BlockIndex
	L2=mps1[2].BlockIndex
	for j=1:length(L1)
		if getfield.(L1, 2)[j] in getfield.(L2, 2)
			id2=findall(isequal(getfield.(L1, 2)[j]), getfield.(L2, 2))[1]
			Block=[GetDiagonalBlock(mps1[2],id2);GetDiagonalBlock(mps1[1], j)]
			U,SV=truncated_svd(Block,δ=δ)
			Update_left!(mps1,U,[id2,j],concat=true)
			Left_DiagonalRXMultiplication!(mps2, SV,getfield.(L1, 2)[j])
		else
			Block=GetDiagonalBlock(mps1[1], j)
			U,SV=truncated_svd(Block,δ=δ)
			Update_left!(mps1,U,[j],concat=false,flag=1)
			Left_DiagonalRXMultiplication!(mps2, SV,getfield.(L1, 2)[j])
		end
	end
	for j=1:length(L2)
		if getfield.(L2, 2)[j]   ∉ getfield.(L1, 2)
			Block=GetDiagonalBlock(mps1[2], j)
			U,SV=truncated_svd(Block,δ=δ)
			Update_left!(mps1,U,[j],concat=false,flag=2)
			Left_DiagonalRXMultiplication!(mps2,SV,getfield.(L2, 2)[j])
		end
	end
	nothing
end

"""
truncated_svd_residual(m::AbstractArray{Float64,2};tol::Float64=1e-12,
degenerate::Bool=true,side::Bool=true)

Perform a truncated SVD to express the matrix m as a factor of
a tall times a wide matrix.
-tol is the Tolerance : the upper bound of the ratio |A_k - A|_F with the Frobenius norm
-degenerate means that the truncations respects the degenerate subspaces
-side: means either left or right truncation
"""
function truncated_svd(m::AbstractArray{Float64,2};δ::Float64=1e-12,
	degenerate::Bool=true,side::Bool=true)
	u,s,v=svd(m)
	s_trunc=sv_trunc(s,δ,degenerate=true);
	if side
		sv_t = (diagm(s_trunc) * transpose(v)[1:length(s_trunc),1:end])
		return u[1:end,1:length(s_trunc)], sv_t
	else
		return u[:,1:length(s_trunc)]*diagm(s_trunc), transpose(v)[1:length(s_trunc),1:end]
	end
end

"""Function sv_trunc(),
Author: Mi-song Dupuy we apply truncated svd according to a tolerance
(ref:https://github.com/msdupuy/Tensor-Train-Julia/)"""
function sv_trunc(s::Array{Float64},tol;degenerate=false,degenerate_eps=1e-10)
	if tol==0.0
		return s
	else
		d = length(s)
		i=0
		weight = 0.0
		norm2 =dot(s,s)
		while (i<d) && weight<tol^2#*norm2
			weight+=s[d-i]^2
			i+=1
		end
		if degenerate  && (d-i+1)!=d
			k=i-1;
			while k>0
				if abs(s[d-k+1]-s[d-i+1])/abs(s[d-k+1])> degenerate_eps
					break
				else
					i-=1;
					k-=1;
				end
			end
		end
		return s[1:(d-i+1)]
	end
end
