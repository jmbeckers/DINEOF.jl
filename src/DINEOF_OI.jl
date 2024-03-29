"""


     XOI=DINEOF_OI(X,U,S,musquare;jlist=1:size(X)[2];)

High level implementation of OI interpolation, typically used when getting new data.

SORRY ONLY IMPLEMENTED FOR last dimension being the "time"

# Input

* `X` : an array with one additional dimension (for several additional images)

* `U,S` from the EOF decomposition at DINEOFrun level


# Output

`XOI` : OI analysis


"""


function DINEOF_OI(XR,UR,S,musquare;jlist=1:size(XR)[end],whichgroups=[ones(Int32,ndims(XR)-1)...,2])



if size(whichgroups)[1]!=ndims(XR)
        @error("Incompatible dimensions in whichgroups")
		return XR
end
if whichgroups!=[ones(Int32,ndims(XR)-1)...,2]
        @error("Sorry, assumes x,t order")
		return XR
end


X=reshape(XR,(prod(size(XR)[1:end-1]),size(XR)[end]))
U=reshape(UR,(prod(size(UR)[1:end-1]),size(UR)[end]))
M=size(X)[1]
N=size(X)[2]
NE=size(S)[1]

#X[isnan.(X)].= 0.0
#U[isnan.(U)].= 0.0

myrefm=mean(XR[.!isnan.(XR)])
#@show myrefm,size(X),size(U)

data=zeros(Float64,M)
XOI=zeros(Float64,(prod(size(UR)[1:end-1]),size(jlist)[1]))

#@show size(XOI)

L=U*diagm(S)/sqrt(N)

L[isnan.(L)].= 0.0
#@show jlist
ij=0
for j in jlist
  
		#@show j
		ij=ij+1
  
		ipresent=.!isnan.(X[:,j])
		
        Lp=L[ipresent,:]
        AA=cholesky(Lp'*Lp+musquare*Matrix{Float64}(I, NE, NE))
        invAAU=inv(AA.U)
        BB=L*invAAU
        data.=0.0
		#myrefm=mean(X[ipresent,j])
		data[ipresent]=X[ipresent,j].-myrefm
		
		XOI[:,ij]=BB*(BB'*data).+myrefm
		XOI[isnan.(U[:,1]),ij].=NaN

end
#@show mean(XOI[.!isnan.(XOI)])

return reshape(XOI,(size(XR)[1:end-1]...,size(XOI)[end]...))

end
