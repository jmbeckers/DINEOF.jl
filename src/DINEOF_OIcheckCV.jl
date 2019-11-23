"""

Calculates the relative difference between DINEOF USV cross-validator estimator and OI analysis using musquare to calculate the cross validator value. The estimate is done using not all columns (images) but
a list jlist. 

"""
function DINEOF_OIcheckCV(X,U,S,V,missingvalues,cvpoints,musquare,jlist,cvEOF)

M=size(X)[1]
N=size(X)[2]
NE=size(S)[1]

data=zeros(Float64,M)
XDcv=zeros(Float64,M)
thetas=0.0
ipoints=0
for j in jlist
  
		
		
       
		icv=cvpoints[findall(x->x==j,cvpoints[:,2]),1]
		ipresent=setdiff(1:M,[missingvalues[findall(x->x==j,missingvalues[:,2]),1]...,icv...])
		

		 			
        
		
		L=U*diagm(S)/sqrt(N)
        Lp=L[ipresent,:]
        AA=cholesky(Lp'*Lp+musquare*Matrix{Float64}(I, NE, NE))
        invAAU=inv(AA.U)
        BB=L*invAAU
        data.=0.0
		data[ipresent]=X[ipresent,j]
		XDcv=BB*(BB'*data)-X[:,j]

		thetas=thetas+sum(XDcv[icv].^2)
		ipoints=ipoints+size(icv)[1]
end
# Output squared relative change in CV value with respect to EOF CV estimator and OI CV estimator itself
return (thetas/ipoints-cvEOF)^2/(cvEOF)^2,thetas/ipoints

end
