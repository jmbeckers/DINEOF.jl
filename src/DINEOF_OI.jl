"""


     XOI=DINEOF_OI(X,U,S,musquare)

High level implementation of OI interpolation, typically used when getting new data.

# Input

* `X` : an array with one additional dimension (for several additional images)

* `U,S` from the EOF decomposition at DINEOFrun level


# Output

`XOI` : OI analysis


"""
function DINEOF_OI

M=size(X)[1]
N=size(X)[2]
NE=size(S)[1]

data=zeros(Float64,M)
XOI=zeros(Float64,size(X))

for j 1:
  
		
		
       
		
		ipresent=setdiff(1:M,[missingvalues[findall(x->x==j,missingvalues[:,2]),1]...,icv...])
  
		
		L=U*diagm(S)/sqrt(N)
        Lp=L[ipresent,:]
        AA=cholesky(Lp'*Lp+musquare*Matrix{Float64}(I, NE, NE))
        invAAU=inv(AA.U)
        BB=L*invAAU
        data.=0.0
		data[ipresent]=X[ipresent,j]
		XOI[:,j]=BB*(BB'*data)

end

return XOI

end
