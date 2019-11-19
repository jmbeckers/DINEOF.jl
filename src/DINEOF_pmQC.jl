function DINEOF_pmQC(X,XA,errmap,musquare,ws;weights=[1.0,1.0,1.0])

# poor mans quality check. Here XA is the EOF filtered field. Normally for consistency the OI-EOF analysed field should be used 
# error maps is variance of error
# added cloud and value indicators as in Aida's paper. 
#
validpoints=.!isnan.(X)
M=prod(size(X))

if errmap==[]
weights[1]=0.0
end

OO=0.
OO2=0.
OO3=0.

weights=weights/sum(weights)

    @show size(X[validpoints])

    if weights[1]>0
	@show size(validpoints)
    o=(X[validpoints].-XA[validpoints])./sqrt.(max.(musquare.-errmap[validpoints],1E-10*musquare))

    om=median(o)
	@show om

    del=1.4826*median(abs.(o.-om))+1E-36
    @show del
    OO=zeros(Float64,size(X))
    OO[validpoints]=abs.(o.-om)./del
	end
    
	if weights[2]>0
	@show "second"
	# will contain the number of neighbors with NaN
	OO2=zeros(Float64,size(X))
	Ifirst, Ilast = CartesianIndices(X)[1], CartesianIndices(X)[M]
	@show Ilast
    I1 = oneunit(Ifirst)
	for I in CartesianIndices(X)
	    if !isnan(X[I])
		for J in max(Ifirst, I-I1):min(Ilast, I+I1)
			if isnan(X[J])
			
				OO2[I]+=1
			end
				
		end
		end
		
		
	end
	end
	
	
	@show "third"
	if weights[3]>0
	OO3=zeros(Float64,size(X))
	
	boxsize=(ws*2+1)^ndims(X)
	values=zeros(Float64,boxsize)
	
	Ifirst, Ilast = CartesianIndices(X)[1], CartesianIndices(X)[M]
	@show Ilast
    I1 = oneunit(Ifirst)
	for I in CartesianIndices(X)
		if !isnan(X[I])
	     
	     n=0
		
			for J in max(Ifirst, I-ws*I1):min(Ilast, I+ws*I1)
				if isnan(X[J])
				
				else
				n+=1
				values[n]=X[J]
				end
			end
		
			#values=values[1:n]
			if n>0
				vm=median(values[1:n])
				del=1.4826*median(abs.(values[1:n] .-vm))+1E-36
		
				OO3[I]=abs(X[I]-vm)/del
			end
		end
		
	end
	end
	#@show OO
	#@show OO2
	#@show OO3
	
    return weights[1]*OO .+ weights[2]*OO2 .+ weights[3]*OO3

end


