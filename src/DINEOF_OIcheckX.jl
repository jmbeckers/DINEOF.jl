"""


    reldeltX,delX=DINEOF_OIcheckX(X,U,S,V,missingvalues,musquare,jlist)
	
	Calculates the relative difference between DINEOF USV decomposition and OI analysis using musquare and the original data of X. The estimate is done using not all columns (images) but a list jlist. The calculation has been optimized exploiting U'U=I to calculate the norm.
	
# Input

* `X` : the reconstructed (filtered array)

* `U,S,V` from the EOF decomposition

* `missingvalues` : list of points where the reconstruction was necessary

* `musquare` : mu^2 used for the OI interpretation

* `jlist`: an array of the columns in which the OI projection is to be compared to the EOF reconstruction

# Output

`reldeltX`: relative difference between OI and EOF reconstructed field

`deltX` : difference between OI and EOF reconstructed field (var)

"""
function DINEOF_OIcheckX(X,U,S,V,missingvalues,musquare,jlist)

M=size(X)[1]
N=size(X)[2]
NE=size(S)[1]

data=zeros(Float64,M)
XD=zeros(Float64,NE)
thetas=0.0

L=U*diagm(S)/sqrt(N)
# Loop over some columns (images) only 
for j in jlist

		# Take the present data only
        ipresent=setdiff(1:M,missingvalues[findall(x->x==j,missingvalues[:,2]),1])
		# See paper on error calculation
		
        Lp=L[ipresent,:]
        AA=cholesky(Lp'*Lp+musquare*Matrix{Float64}(I, NE, NE))
        invAAU=inv(AA.U)
        BB=L*invAAU
        data.=0.0
		data[ipresent]=X[ipresent,j]
		# Normally the OI project is L*invAUU*invAUU*L'*data and to be compared to U S V' which is U S/sqrt(N) *invAUU * BB'*data - U S V 
		XD=S.*(invAAU*BB'*data/sqrt(N) - V[j,:])
		thetas=thetas+sum(XD.^2)
end
# Output relative (OI-EOF) error variance with respect to signal variance and error variance XOI-XEOF itself
return (thetas/(NE*size(jlist)[1]))/(sum(S.^2)/size(S)[1]),(thetas/(NE*size(jlist)[1]))

end
