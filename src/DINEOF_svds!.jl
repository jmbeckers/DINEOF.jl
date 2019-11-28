"""


    U,S,V,cvEOF,cvarray,musquare=DINEOF_svds!(X,missingvalues,crossvalidation;
	keeprestart=true,
	ncmax=size(X)[2]-1,
	istart=1,
	dineofmaxiter=10,
	dineoftol=0.001,
	svdmeth="svd",
	svdtol=0.000001,
	filter="None",
	filterintensity=1.0,
	filterrepetitions=1
	)

# DINEOF SVD decomposition with filling in of missing points.


# Input: 
* `X`: a two-dimensional array of size MxN with N<=M already filled in with a first guess in points with missing values

* `missingvalues`: array of size Px2 collecting the indexes for the P missing values in location missingvalues[:,1] and missingvalues[:,2]. There is no distinction between missing values and topologically excludes points. All will be filled

* `crossvalidation`:  array of size Qx2 collecting the indexes for the Q points to be used for cross validation

* `keeprestart`: Boolean. If true the method keeps track of the best reconstruction up to now while doing further iterations. 

* `ncmax`: maximum number of eigenvalues calculated

* `itstart` : defines the number of modes used during the first iteration. Then one at a time is added

* `dineofmaxiter` : Maximum number of iterations during a SVD decomposition/filling in loop

* `dineoftol`: Defines the tolerance below which the iterations are stopped. Relative change of xxxxx

For the other parameters, see `DINEOF_svd`



# Output:

* `U,S,V`: the filled SVD decomposition such that in the infiltered version U*S*V' is the best approximation of X with nele singular vectors.

* `cvEOF` : the cross validation estimator of the reconstrucion (variance of misfit at crossvalidation points)

* `cvarray` : the cross validation estimator for different number of retained EOFs 

* `musquare` : estimation of mu^2

   WARNING: X is updated at the missing data points !!!!!!



"""
function DINEOF_svds!(X,
	missingvalues,
	crossvalidation;
	keeprestart=true,
	ncmax=size(X)[2]-1,
	istart=1,
	dineofmaxiter=10,
	dineoftol=0.001,
	svdmeth="svd",
	svdtol=0.000001,
	filter="None",
	filterintensity=1.0,
	filterrepetitions=1
	)
    # X array to analyze size MxN with M>N
    # missingvalues array of indexes to missing points (P,2)
    # crossvalidation: array of indexes to cross validation points (C,2)
    # keeprestart: boolean. If true keeps track of interpolation at the best value of eofs. If false, uses the lastest estimates
    # with some irrelevant eofs and performs the final reconstruction from there
    # ncmax: maximum number of eofs
    # istart: starting point into the eof decomposition
    # Changes X !!!!
    # To get the filtered matrix, use the output to calculate USV'
	
	println("svds matrix size: $(size(X))")
    
    varmatrix=var(X)
	meanmatrix=mean(X)
	
	
    squaremiss=0
    for jjj=1:size(missingvalues)[1]
        squaremiss=squaremiss+(X[missingvalues[jjj,1],missingvalues[jjj,2]]-meanmatrix)^2 
    end
    varmatrixp=(varmatrix*prod(size(X))-squaremiss)/(prod(size(X))-size(missingvalues)[1])
    varmatrixm=squaremiss/size(missingvalues)[1]
	
	
	if isnan(varmatrix)
	
	@err("Input to svds! should not include NaN anymore")
	end
	
	
	println("svds!: variance and mean of the entry matrix: $varmatrix , $meanmatrix ; intial variance at points to fill in: $varmatrixm ")
	
	
    #@show varmatrix,meanmatrix,varmatrixp,varmatrixm
    if meanmatrix^2>0.0001*varmatrix
        @warn("You should subtract a mean value from your data")
    end
    # cross validation
    cvval=1E38
    cvbest=1E38
    ibest=1
    # Make sure SV it is saved from calculations in loops
    SV=[]
    SVU=[]
    SVS=[]
    SVV=[]
    # keeps cross validation values
    cvvalues=zeros(Float64,size(crossvalidation)[1])
    
    # keeps reconstructed values at the optimum
    if keeprestart
		forfinalloop=zeros(Float64,size(missingvalues)[1])
    end
    
    #stores all estimators
    cv=zeros(Float64,ncmax)
    
    
    # Remember original values
    for jj=1:size(crossvalidation)[1]
        i=crossvalidation[jj,1]
        j=crossvalidation[jj,2]
        cvvalues[jj]=X[i,j]
    end
    
    
    # Once a local minimum is found, moreloops loops are added. If no new minimum is found the search is stopped
    moreloops=4
    iloop=istart-1
    while moreloops>0
        iloop=iloop+1
    
       
        # Loop over maximum number of iterations. If converged before a break is forced
        for iterc=1:dineofmaxiter
            #@show svdmeth,filter,filterintensity
			# SVD decomposition on a full matrix with some parameters
            SVU,SVS,SVV,ncon,nit=DINEOF_svd(X,min(iloop,ncmax),svdmeth,svdtol;filter=filter,filterintensity=filterintensity,filterrepetitions=filterrepetitions)
            
            # Check if the standard SVD decomposition worked out (basically if the matrix was not ill behaved)
            if ncon!=min(iloop,ncmax)
                @warn("Convergence problem in svds L1, will try to continue")
                @show iloop,ncon,nit
            end
        
		# Update missing points and keep track of changes amplitude in varchange
        varchange=0
        for jj=1:size(missingvalues)[1]
            i=missingvalues[jj,1]
            j=missingvalues[jj,2]
            zz=SVU[i:i,:]*diagm(SVS)*(SVV[j:j,:])' 
                varchange=varchange+(X[i,j]-zz[1])^2
            X[i,j]=zz[1]
        end
            #@show varchange/size(missingvalues)[1],varmatrix
            
            if varchange<dineoftol^2*varmatrix*size(missingvalues)[1]
                println("Convergence for $iloop eofs, relative change $(sqrt((varchange/size(missingvalues)[1])/varmatrix)) after $iterc iterations")
                break
            end
        end
        
		# Once finished, calculate CV estimator
        cvval=0
        for jj=1:size(crossvalidation)[1]
            i=crossvalidation[jj,1]
            j=crossvalidation[jj,2]
            zz=SVU[i:i,:]*diagm(SVS)*(SVV[j:j,:])' 
            X[i,j]=zz[1]
            cvval=cvval+(X[i,j]-cvvalues[jj])^2
        end
        
		# No cv points
        if size(crossvalidation)[1]==0
		# without cv points, force estimator to decrease so that the requested number of EOFs can be reached
            cvval=1E37/iloop
        end
        #@show cvval
		tutu=cvval/max(1,size(crossvalidation)[1])
		println("Eof loop $(iloop) with mean squared misfit: $(tutu) ")
        cv[iloop]=cvval
        # Avoid too flat cross validation curve near minimum and only accept "significant decrease"
        if cvval-cvbest<=-0.0001*cvbest
            cvbest=cvval
            ibest=iloop
            # If you want a clean final loop, store here the optimal data estimates at the missing points
            if keeprestart
				for jj=1:size(missingvalues)[1]
					i=missingvalues[jj,1]
					j=missingvalues[jj,2]
					# Not sure a deepcopy is needed, but to be sure ...
					forfinalloop[jj]=deepcopy(X[i,j])
				end
            end
            moreloops=4
          else
            moreloops=moreloops-1
        end
        if iloop>=ncmax
            moreloops=0
        end
      
    end
    
    # Ok, now we have the best value for number of EOF, throw in the original data and calculate final decomposition 
    
    if keeprestart
        # put back the best estimate up to now
		for jj=1:size(missingvalues)[1]
            i=missingvalues[jj,1]
            j=missingvalues[jj,2]
            X[i,j]=forfinalloop[jj]
		end
    end
    
    # put back the cross validation values
    for jj=1:size(crossvalidation)[1]
            i=crossvalidation[jj,1]
            j=crossvalidation[jj,2]
            X[i,j]=cvvalues[jj]
    end
    # make a final eof decomposition; possibly be more severe here in terms of tolerance
    
    for iterc=1:dineofmaxiter
    varchange=0
        SVU,SVS,SVV,ncon,nit=DINEOF_svd(X,min(ibest,ncmax),svdmeth,svdtol;filter=filter,filterintensity=filterintensity,filterrepetitions=filterrepetitions)
        
        
        
        if ncon!=min(ibest,ncmax)
            @warn("Convergence problem in final svds")
            @show ibest,ncon,nit
        end
        
    for jj=1:size(missingvalues)[1]
            i=missingvalues[jj,1]
            j=missingvalues[jj,2]
            zz=SVU[i:i,:]*diagm(SVS)*(SVV[j:j,:])' 
            varchange=varchange+(X[i,j]-zz[1])^2
            X[i,j]=zz[1]
    end
        if varchange<dineoftol^2*varmatrix*size(missingvalues)[1]
                #@show ibest,varchange/size(missingvalues)[1]
                break
        end
    end
    
    # Provide variance
    if size(crossvalidation)[1]>0
	#@show size(crossvalidation)[1]
	#@show cvbest
        cvbest=cvbest/size(crossvalidation)[1]
        cv=cv./size(crossvalidation)[1]
    end
    println("Cross validation value (mean squared misfit): $cvbest for $ibest EOFs")
    #@show cv[1:iloop]
	########################################################################
	# To work on
	# Now some statistics on variances and estimates of musquare (observational error covariance)
    #@show varmatrix,sum(SVS.^2)/prod(size(X)),var(X)
    
    #varmatrixf=var(X)
    #squaremiss=0
    #for jjj=1:size(missingvalues)[1]
    #    squaremiss=squaremiss+X[missingvalues[jjj,1],missingvalues[jjj,2]]^2 
    #end
    #varmatrixpp=(varmatrixf*prod(size(X))-squaremiss)/(prod(size(X))-size(missingvalues)[1])
    
    #varmatrixfp=(sum(SVS.^2)-squaremiss)/(prod(size(X))-size(missingvalues)[1])
    
    #musquare=varmatrixp-varmatrixfp
	# Adaptive estimate of mean(diag(R))  assuming U*S*V' is equivalent of an OI analysis
	# Should be better than our version of the paper
	musquare=sum(X .* X .- X.* (SVU*diagm(SVS)*SVV') )/(prod(size(X))-size(missingvalues)[1])
    println("Estimation for musquare based on DeRozier type of analysis: $musquare")
	println("Estimation of mean error variance of reconstuctions: $(musquare-cvbest) ")
    #@show musquare,musquare/varmatrix
    if musquare<0.000001*varmatrix
        @warn("Very low level of noise ? Fraction $(musquare/varmatrix) of total variance")
        musquare=0.000001*varmatrix
		@show musquare
    end
    
    
    if sum(SVS.^2) > varmatrix*prod(size(X))
        @warn("Initial Variance has been increased for filtered matrix  by factor $(sum(SVS.^2)/(varmatrix*prod(size(X))))")
		else
        println("Explained variance  $(100*sum(SVS.^2)/(varmatrix*prod(size(X)))) percent")
    end
    # maybe add musquare to output parameters ? replace cvbest by [cvbest,musquare ].... whatever since for the moment cvbest was never used
	#
    ###########################################################################
    
    return SVU,SVS,SVV,cvbest,cv[1:iloop],musquare
end