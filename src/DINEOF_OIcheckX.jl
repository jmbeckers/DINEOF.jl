"""

Calculates the relative difference between DINEOF USV decomposition and OI analysis using musquare. The estimate is done using not all columns (images) but
a list jlist. The calculation has been optimized exploiting U'U=I to calculate the norm.

"""
function DINEOF_OIcheckX(X,U,S,V,missingvalues,musquare,jlist)

M=size(X)[1]
N=size(X)[2]
NE=size(S)[1]

data=zeros(Float64,M)
XD=zeros(Float64,NE)
thetas=0.0
for j in jlist
  
		
        ipresent=setdiff(1:M,missingvalues[findall(x->x==j,missingvalues[:,2]),1])
		 			
        
		
		L=U*diagm(S)/sqrt(N)
        Lp=L[ipresent,:]
        AA=cholesky(Lp'*Lp+musquare*Matrix{Float64}(I, NE, NE))
        invAAU=inv(AA.U)
        BB=L*invAAU
        data.=0.0
		data[ipresent]=X[ipresent,j]
		XD=S.*(invAAU*BB'*data/sqrt(N) - V[j,:])
		thetas=thetas+sum(XD.^2)
end
# Output relative (OI-EOF) error variance with respect to signal variance and error variance
return (thetas/(NE*size(jlist)[1]))/(sum(S.^2)/size(S)[1]),(thetas/(NE*size(jlist)[1]))

end
