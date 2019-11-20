"""


    U,S,V,=DINEOF_svds!()

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

For the other parameters, see DINEOF_svd



# Output:

* `U,S,V`: the filled SVD decomposition such that in the infiltered version U*S*V' is the best approximation of X with nele singular vectors.

* WARNING: X is updated at the missing data points !!!!!!



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
    
    #svdmeth="svd"
    
    #svdtol=0.000001
    varmatrix=var(X)
    squaremiss=0
    for jjj=1:size(missingvalues)[1]
        squaremiss=squaremiss+X[missingvalues[jjj,1],missingvalues[jjj,2]]^2 
    end
    varmatrixp=(varmatrix*prod(size(X))-squaremiss)/(prod(size(X))-size(missingvalues)[1])
    meanmatrix=mean(X)
    @show varmatrix,meanmatrix,varmatrixp
    if meanmatrix^2>0.00001*varmatrix
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
    
    
    
    moreloops=4
    iloop=istart-1
    while moreloops>0
        iloop=iloop+1
    
        # Try to implement loose iterations to start with (higher tolerance) and then more and more severe ones
        # Another tweek could be the starting vector of the Krylov iterations that could be retrieved from previous iteration
        #
        #if iloop>istart
        #SV=svds(X;nsv=min(iloop,ncmax),v0=SV.V[:,1])[1]
        #else
        # NEED TO ADD A CONVERGENCE LOOP FOR FIXED NUMBER OF EOF. 
        
        for iterc=1:dineofmaxiter
            @show svdmeth,filter,filterintensity
			
            SVU,SVS,SVV,ncon,nit=DINEOF_svd(X,min(iloop,ncmax),svdmeth,svdtol;filter=filter,filterintensity=filterintensity,filterrepetitions=filterrepetitions)
            
            
            if ncon!=min(iloop,ncmax)
                @warn("Convergence problem in svds")
                @show iloop,ncon,nit
            end
        #end
            varchange=0
        for jj=1:size(missingvalues)[1]
            i=missingvalues[jj,1]
            j=missingvalues[jj,2]
            zz=SVU[i:i,:]*diagm(SVS)*(SVV[j:j,:])' 
                varchange=varchange+(X[i,j]-zz[1])^2
            X[i,j]=zz[1]
        end
            @show varchange/size(missingvalues)[1],varmatrix
            
            if varchange<dineoftol^2*varmatrix*size(missingvalues)[1]
                @show iloop,varchange
                break
            end
        end
        
        cvval=0
        for jj=1:size(crossvalidation)[1]
            i=crossvalidation[jj,1]
            j=crossvalidation[jj,2]
            zz=SVU[i:i,:]*diagm(SVS)*(SVV[j:j,:])' 
            X[i,j]=zz[1]
            cvval=cvval+(X[i,j]-cvvalues[jj])^2
        end
        
        if size(crossvalidation)[1]==0
            cvval=1E37
        end
        @show cvval
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
            forfinalloop[jj]=X[i,j]
            end
            end
            
            
            
            #
            
            @show iloop,cvbest
            moreloops=4
            
        else
            moreloops=moreloops-1
        end
        if iloop>=ncmax
            moreloops=0
        end
      
    end
    
    # Ok, now we have the best value, now throw in the original data and calculate final decomposition (no rolling back for the moment
    # hoping X was not to much deteriorated due to the additional vectors)
    #
    @show ibest,cvbest
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
        
        #if svdmeth=="svd"
        #    SV,ncon,nit=svds(X;nsv=min(ibest,ncmax),tol=svdtol)[1:3]
        #    SVS=SV.S
        #    SVV=SV.V
        #    SVU=SV.U
        #end    
        
        #if svdmeth=="eig"
        #    SVS,SVV,ncon,nit=eigs(X'*X;nev=min(ibest,ncmax),tol=svdtol)[1:4]
        #    SVS=sqrt.(abs.(SVS))
        #    SVU=X*SVV*diagm(1.0 ./SVS)
        #end
        
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
                @show iloop,varchange
                break
        end
    end
    
    # Provide variance
    if size(crossvalidation)[1]>0
        cvbest=cvbest/size(crossvalidation)[1]
        cv=cv./size(crossvalidation)[1]
    end
    
    @show cv[1:iloop]
    @show varmatrix,sum(SVS.^2)/prod(size(X)),var(X),prod(size(X)),size(SVS)
    
    varmatrixf=var(X)
    squaremiss=0
    for jjj=1:size(missingvalues)[1]
        squaremiss=squaremiss+X[missingvalues[jjj,1],missingvalues[jjj,2]]^2 
    end
    varmatrixpp=(varmatrixf*prod(size(X))-squaremiss)/(prod(size(X))-size(missingvalues)[1])
    
    varmatrixfp=(sum(SVS.^2)-squaremiss)/(prod(size(X))-size(missingvalues)[1])
    
    musquare=varmatrixp-varmatrixfp
	
	newmusquare=mean(X .* X .- X.* (SVU*diagm(SVS)*SVV') )*(prod(size(X))-size(missingvalues)[1])/prod(size(X))
    @show musquare,varmatrixp,varmatrixfp,varmatrixpp,newmusquare
    
    if musquare<0.000001*varmatrixp
        @warn("Strange")
        musquare=0.000001*varmatrixp
    end
    
    
    if sum(SVS.^2) > varmatrix*prod(size(X))
        @warn("Initial Variance has been increased for filtered matrix  by factor $(sum(SVS.^2)/(var(X)*prod(size(X))))")
        
    end
    
    
    
    return SVU,SVS,SVV,cvbest,cv[1:iloop]
end