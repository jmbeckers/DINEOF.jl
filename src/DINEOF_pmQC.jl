"""

      OutlierIndicator=DINEOF_pmQC(X,XA,errmap,musquare,ws;weights=[1.0,1.0,1.0])
	  
# Outlier detection


# Input: 
* `X`: an NDimensional array of original data with NaN at missing points

* `XA` : an NDimensional array which contains the analysed field (the EOF reconstruction) at all points (present and absent points)

* `errmap` : an NDimensional array containing the error variance estimate. If errmap=[] the outlier detection based on residuals is not exploited

* `musquare`: estimated value (effective) musquare

* `ws`: width of the window for outlier detection based on data variance. the box size is (2*ws+1)^NDIMS

* `weights` : array of weights for the three outlier indicators (residuals, proximity to NaNs, data variability). If a weight is zero the corresponding indicator is not calculated. weights are normalized internally to have a sum of 1

* `statmeth`: if meanstd mean and std are used, otherwise median and MAD (much slower)

# Output:



* `OutlierIndicator` : an ND array with the same topology as X with values indicating suspect data with high values.


"""
function DINEOF_pmQC(X,XA,errmap,musquare,ws;weights=[1.0,1.0,1.0],statmeth="meanstd")

# poor mans quality check. Here XA is the EOF filtered field. Normally for consistency the OI-EOF analysed field should be used 
# error maps is variance of error
# added cloud and value indicators as in Aida's paper. 
#

		
		M=prod(size(X))

# indicator 1 needs the error map
		if errmap==[]
			weights[1]=0.0
		end

	OO=0.
	OO2=0.
	OO3=0.

	
# weights are forced to have sum=1
	
	weights=weights/sum(weights)

# QC on difference between data and analysis
    if weights[1]>0
	    validpoints=.!isnan.(X) .& .!isnan.(XA)
		@show size(validpoints)
		o=(X[validpoints].-XA[validpoints])./sqrt.(max.(musquare.-errmap[validpoints],1E-10*musquare))
		om=median(o)
		@show om
		del=1.4826*median(abs.(o.-om))+0.1*abs(om)
		@show del
		OO=zeros(Float64,size(X))
		OO[validpoints]=min.(abs.(o.-om)./del,10.0)
	end
# QC on touching NaNs
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
	

#QC on local deviation of point from surrounding values
# Actually quite exensive in huge arrays/dimensions
# Maybe check if mean/var is much quicker ??
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
					if n>0
					if statmeth=="meanstd"
						#vm=median(values[1:n])
						vm=mean(values[1:n])
						#del=1.4826*median(abs.(values[1:n] .-vm))+1E-36
						del=std(values[1:n])
						else
						vm=median(values[1:n])
						#m=mean(values[1:n])
						del=1.4826*median(abs.(values[1:n] .-vm))+1E-36
						#del=std(values[1:n])
					end
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


