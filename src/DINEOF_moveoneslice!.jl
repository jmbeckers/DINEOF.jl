


function DINEOF_moveoneslice!(XF,XR,UR,S,musquare,newslice;whichgroups=[ones(Int32,ndims(XR)-1)...,2],
minimumcoverage=(0.1, 0.1),
	cvmask="Automatic",
	cvfraction=0.01,
	cvmethod="Random",
	maxbubblesize=max.(20,0.01*[size(XR)...]),
	dimensionsforcopy=[zeros(Int32,ndims(XR)-1)...,1],
	errormap=true,
	restart=[],
    keepsvdcvrestart=true,
	eofmax=min(size(S)[1]+1,size(XF)[end]),
	eofstart=max(size(S)[1]-1,1),
	dineofmaxiter=15,
	dineoftol=0.005,
	svdmeth="svd",
	svdtol=0.000001,
	tfilter="None",
	filterintensity=1.0,
	filterrepetitions=1)

@show size(S),S,"in"


if size(whichgroups)[1]!=ndims(XR)
        @error("Incompatible dimensions in whichgroups")
		return XR
end
if whichgroups!=[ones(Int32,ndims(XR)-1)...,2]
        @error("Sorry, assumes x,t order")
		return XR
end
if size(newslice)!=size(XR)[1:end-1]
         @error("Sorry, slice has not the right size")
		 @show size(newslice),size(XR)
end


#Typically one image to add at the end and drop image at start



X=reshape(XR,(prod(size(XR)[1:end-1]),size(XR)[end]))
XA=reshape(XF,(prod(size(XF)[1:end-1]),size(XF)[end]))
U=reshape(UR,(prod(size(UR)[1:end-1]),size(UR)[end]))
M=size(X)[1]
N=size(X)[2]
NE=size(S)[1]
# Move data within the array (in place)
for j=2:N
	X[:,j-1].=X[:,j]
	XA[:,j-1].=XA[:,j]
end
X[:,end].=reshape(newslice,(prod(size(XR)[1:end-1]),))

# Fill in new image using DINEOF_OI
@show size(X),size(U),size(S),musquare,N,NE,M

newsliceoi=DINEOF_OI(X[:,end:end],U,S,musquare)
XA[:,end].=newsliceoi







# DINEOFrun with restart XA and one EOF less than previous one and limit to one more EOF then previous one
XAN,offset,Uout,Sout,Vout,cvEOF,cvarray,errmap,muout=DINEOFrun(XR;
	minimumcoverage=minimumcoverage,
	cvmask=cvmask,
	cvfraction=cvfraction,
	cvmethod=cvmethod,
	maxbubblesize=maxbubblesize,
	dimensionsforcopy=dimensionsforcopy,
	errormap=errormap,
	musquare=musquare,
	restart=XF,
    keepsvdcvrestart=keepsvdcvrestart,
	eofmax=eofmax,
	eofstart=eofstart,
	dineofmaxiter=dineofmaxiter,
	dineoftol=dineoftol,
	svdmeth=svdmeth,
	svdtol=svdtol,
	tfilter=tfilter,
	filterintensity=filterintensity,
	filterrepetitions=filterrepetitions)

# Once reconstructed, weighted average between old reconstruction and new reconstrucion. Weighting with a linear ramp ? So newest image mostly by corrent reconstrucion while older images basicall untouched by new reconstruction?
@show size(XR),size(XAN)

XAN=reshape(XAN,(prod(size(XAN)[1:end-1]),size(XAN)[end]))

for j=1:N
    w=j/N
	#@show w
	XA[:,j].= w*XAN[:,j] .+ (1-w)*XA[:,j]
	#XA[:,j].=XAN[:,j]
end

@show size(S),S
# Saving is left to the calling routine but probably the  oldest image has the best estimate (most data used) 

return offset,Uout,Sout,Vout,cvEOF,cvarray,errmap,muout

end