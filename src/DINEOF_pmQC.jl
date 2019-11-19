function DINEOF_pmQC(X,XA,errmap,musquare,ws)

# poor mans quality check. Here XA is the EOF filtered field. Normally for consistency the OI-EOF analysed field should be used 
# error maps is variance of error
# TODO: add cloud and value indicators as in Aida's paper. Exploit multidimensional approach of bubble generation to make the ND generalisation
#
validpoints=.!isnan.(X)
M=prod(size(X))

    @show size(X[validpoints])

    o=abs.(X[validpoints].-XA[validpoints])./sqrt.(max.(musquare.-errmap[validpoints],1E-10*musquare))

    om=median(o)

    del=1.4826*median(abs.(o.-om))
    
    OO=zeros(Float64,size(X))
    OO[validpoints]=abs.(o.-om)./del
    
	# will contain the number of neighbors with NaN
	OO2=zeros(Float64,size(X))
	Ifirst, Ilast = CartesianIndices(X)[1], CartesianIndices(X)[M]
    I1 = oneunit(Ifirst)
	for I in CartesianIndices(X)
	    
		for J in max(Ifirst, I-I1):min(Ilast, I+I1)
			if isnan(X[J])
				OO2[I]+=1
				
			end
		end
		
		
	end
	
	
	
	
	
	OO3=zeros(Float64,size(X))
	
	boxsize=(ws*2+1)^ndims(X)
	
	
	Ifirst, Ilast = CartesianIndices(X)[1], CartesianIndices(X)[M]
    I1 = oneunit(Ifirst)
	for I in CartesianIndices(X)
	     values=zeros(Float64,boxsize)
	     n=0
		for J in max(Ifirst, I-ws*I1):min(Ilast, I+ws*I1)
			if isnan(X[J])
				
				else
				n+=1
				values[n]=X[J]
			end
		end
		values=values[1:n]
		if n>0
		vm=median(values)
		del=1.4826*median(abs.(values.-vm))+1E-36
		if !isnan(X[I])
		  OO3[I]=abs(X[I]-vm)/del
		end
		end
		
	end
	@show OO
	@show OO2
	@show OO3
	
    return OO.+OO2.+OO3

end


