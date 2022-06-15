"""


     reldelCV,CVOI=DINEOF_OIcheckCV(X,U,S,V,missingvalues,cvpoints,musquare,jlist,cvEOF)

Calculates the relative difference between DINEOF USV cross-validator estimator and OI analysis using musquare to calculate the cross validator value. The estimate is done using not all columns (images) but
a list jlist. 

# Input

* `X` : the reconstructed (filtered array)

* `U,S,V` from the EOF decomposition

* `missingvalues` : list of points where the reconstruction was necessary

* `cvpoints` : list of points used for cross validation

* `musquare` : mu^2 used for the OI interpretation

* `jlist`: an array of the columns in which the OI projection is to be compared to the EOF reconstruction

* `cvEOF` : value of the cross validation estimator from EOF decomposition

# Output

`reldelCV`: relative difference between OI and EOF CV values

`CVOI` : Cross validator from OI reconstruction


"""
function DINEOF_OIcheckCV(X,U,S,V,missingvalues,cvpoints,musquare,jlist,cvEOF)

M=size(X)[1]
N=size(X)[2]
NE=size(S)[1]

data=zeros(Float64,M)
XDcv=zeros(Float64,M)
thetas=0.0
ipoints=0
L=U*diagm(S)/sqrt(N)
for j in jlist
  
		
		
       
		icv=cvpoints[findall(x->x==j,cvpoints[:,2]),1]
		ipresent=setdiff(1:M,[missingvalues[findall(x->x==j,missingvalues[:,2]),1]...,icv...])
		

		 			
        
		
		
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

if ipoints==0
# Dont penalize search for best mu if you cannot calculate CV
    @warn("No cross-validation points found in the columns analyzed")
    return 0,1E36
end

return (thetas/ipoints-cvEOF)^2/(cvEOF)^2,thetas/ipoints

end
