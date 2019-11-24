"""


       DINEOFrun(X,whichgroups;minimumcoverage=(0.1, 0.1),cvmask="Automatic",cvfraction=0.01,cvmethod="Random",errormap=true,musquare=0,restart=[])





"""
function DINEOFrun(X,whichgroups;minimumcoverage=(0.1, 0.1),cvmask="Automatic",cvfraction=0.01,cvmethod="Random",errormap=true,musquare=0,restart=[],
    keeprestart=true,
	ncmax=size(X)[2]-1,
	istart=1,
	dineofmaxiter=10,
	dineoftol=0.001,
	svdmeth="svd",
	svdtol=0.000001,
	filter="None",
	filterintensity=1.0,
	filterrepetitions=1)
    # X ND array with NaN at missing locations
    # whichgroups: array of 1 and 2 indicating which dimensions are collapsed together. eg [1 2 1 2] regroups dimensions 
    # 1 and 3 into i direction and 2 with 4 into j direction
    
	
	#TODOTODO take out mean !!!!!
	datamean=mean(X[.!isnan.(X)])
	datavar=var(X[.!isnan.(X)])
	println("Raw data variance and mean: $datavar and $datamean")
	println("Number of missing points (including possible masks): $(sum(isnan.(X))) out of $(prod(size(X)))")
	X=X.-datamean
	
	
	
	errmap=[]
	
    
    if cvmask=="Automatic"
        
        cvmask=DINEOF_cvmask(X,cvfraction;cvmethod=cvmethod)
        
    else
        #@show size(cvmask),size(X)
        if size(cvmask)!=size(X)
            @error("Mask for cross validation and size of X do not correspond")
			@show size(X),size(cvmask)
            return
        end
    end
    
    # Do not use cross validation in already missing points
    cvmask[isnan.(X)].=false
	
	println("Number of data points before elimination of low coverage regions is $(sum(.!isnan.(X))) and cv fraction $(sum(cvmask)/sum(.!isnan.(X)))")
	
    
    if size(whichgroups)[1]!=ndims(X)
        @error("Incompatible dimensions in whichgroups")
    end
    #@show ndims(X),size(whichgroups)[1]
    #
    g1=findall(x->x==1,whichgroups)
    g2=findall(x->x==2,whichgroups)
    if size(g1)[1]<1 || size(g2)[1]<1 || size(g1)[1]+size(g2)[1] != ndims(X)
        @error("Incorrect whichgroups")
        @show whichgroups
    end
    perminput=[g1...,g2...]
    #@show perminput
    newsize=(prod(size(X)[g1]),prod(size(X)[g2]))
    #@show size(X)
    #@show size(X)[perminput]
    sizeperminput=size(X)[perminput]
    X2D=reshape(permutedims(X,perminput),newsize)
    cv2D=reshape(permutedims(cvmask,perminput),newsize)
    #@show size(X2D)
    
    # now take out lines and columns with less the specified coverage
    # 
    #
    #@show size(isnan.(X2D))
    #@show sum(cv2D)
     
    ccov=dropdims(sum(.!isnan.(X2D),dims=1)',dims=2) .- dropdims(sum(cv2D,dims=1)',dims=2)
    lcov=dropdims(sum(.!isnan.(X2D),dims=2),dims=2) .- dropdims(sum(cv2D,dims=2),dims=2)
    #@show ccov,size(ccov)
    #@show lcov,size(lcov)
    rlow=findall(x->x<newsize[2]*minimumcoverage[1],lcov)
    clow=findall(x->x<newsize[1]*minimumcoverage[2],ccov)
    
    #@show rlow,clow
    # Usefull functions: take out rows and colums with low coverage and when finished later put NaNs back
    not_in(inds::AbstractVector{Int}, n::Int)::Vector{Int} = setdiff(1:n, [ i for i in inds ])
    not_in2(inds::AbstractVector{Int}, n::Int)::Vector{Int} = setdiff(1:n,inds)
  
    X2D=X2D[not_in(rlow, size(X2D,1)), not_in(clow, size(X2D,2))]
    cv2D=cv2D[not_in(rlow, size(cv2D,1)), not_in(clow, size(cv2D,2))]
    
    
    # now transpose matrix if necessary
    transposed=false
    if size(X2D)[1]<size(X2D)[2]
        transposed=true
        X2D=permutedims(X2D,[2,1])
        cv2D=permutedims(cv2D,[2,1])
    end
    #######################################################################
	# Some more work here to have info for users
    # Ok, now we have the working matrix. Do some stats on it:
    #
    varmatrix=var(X2D[.!isnan.(X2D)])
    missingpointsforvar=deepcopy(.!isnan.(X2D))
    #@show size(missingpointsforvar),size(X2D)
    meanmatrix=mean(X2D[.!isnan.(X2D)])
    #@show varmatrix,meanmatrix
    
    NM=sum(isnan.(X2D))
    NMCV=sum(cv2D)
    #@show NM,NMCV
	println("Number of data points after elimination of low coverage regions is $(prod(size(X2D))-NM) and cv fraction $(NMCV/(prod(size(X2D))-NM))")
	
    ##############################
	
	
    # Now create index arrays pointing to the points in the matrix with NaN
    # Probably not very elegant but hey...
    missingvalues=zeros(Int,(NM,2))
    icount=0
    for j=1:size(X2D)[2]
        for i=1:size(X2D)[1]
            if isnan(X2D[i,j])
                icount=icount+1
                missingvalues[icount,1]=i
                missingvalues[icount,2]=j
                # Put a random value with average variance of present data so that we can keep an eye on how total variane
                # behaves.
                X2D[i,j]=sqrt(varmatrix)*randn()
            end
        end
    end
    
    # Deal with cross-validation ... where?
    cvpoints=zeros(Int,(NMCV,2))
    icount=0
    for j=1:size(X2D)[2]
        for i=1:size(X2D)[1]
            if  cv2D[i,j]
                icount=icount+1
                cvpoints[icount,1]=i
                cvpoints[icount,2]=j
                
            end
        end
    end
    # free cv2D
    cv2D=[]
    
    
    
    # NEED TO ADD OPTIONAL PARAMETERS ...
    U,S,V,cva,cvb,musquareestimate=DINEOF_svds!(X2D,missingvalues,cvpoints;
	keeprestart=keeprestart,
	ncmax=ncmax,
	istart=istart,
	dineofmaxiter=dineofmaxiter,
	dineoftol=dineoftol,
	svdmeth=svdmeth,
	svdtol=svdtol,
	filter=filter,
	filterintensity=filterintensity,
	filterrepetitions=filterrepetitions)
	# Decide here on musquare, error maps and QC estimators
	# Get it back from svds and finetune with DINEOF_musquare around-way above the proposed value (inflation)
	
	if musquare==0
		
		musquare,cvopt=DINEOF_musquare(X2D,U,S,V,missingvalues,cvpoints,cva;musquarer=[0.5*musquareestimate,100*musquareestimate])[1:2]
		
		
		
		
		
		
		println("Estimated musquare $musquareestimate was inflated by factor $(musquare/musquareestimate) into $musquare")
		println("This optimal value provides OI interpolation CV estimator $cvopt")
	end
	#musquare=var(X2D)
	
	if errormap
		errmap=DINEOF_errormap(U,S,V,musquare,missingvalues)
		@show mean(errmap),cvopt-musquare
	end
	
	#
    # now roll back
    #@show X2D
    #@show size(U),size(S),size(V)
    # Transpose back if that was done (also SVD...) NEED TO CHECK IF V and be used as U without problem??? (Adjoint definition)
    if transposed
	@show "Transposed"
        TT=deepcopy(U)
        U=deepcopy(V)
        V=deepcopy(TT)
        TT=[]
        X2D=permutedims(X2D,[2,1])
		if errormap
        errmap=permutedims(errmap,[2,1])
		end
    end
    
    # put back lines and columns (also SVD)
    
    # Not so easy ....
    # Lines rlow into U columns clow into V' so into lines of V
    # Check if Alex version is quicker
    # 
    
    #X2D=DINEOF_insertNaNr(X2D,rlow)
    #X2D=DINEOF_insertNaNc(X2D,clow)
    
    #@show size(U),size(V),size(S)
    XF2D=U*diagm(S)*V'
    #@show size(XF2D)
    #@show size(missingpointsforvar)
    varmatrixf=var(XF2D[missingpointsforvar])
    
    #@show varmatrix,varmatrixf,varmatrix-varmatrixf
    
    #@show "back"
	# Filled matrix NOT necessary if DINEOF_fuse is used
    #A=fill(NaN,newsize)
    #A[not_in2(rlow, end), not_in2(clow, end)] = X2D
    #X2D=A
    
    A=fill(NaN,newsize)
    A[not_in2(rlow, end), not_in2(clow, end)] = XF2D
    XF2D=A
    
    A=fill(NaN,(newsize[1],size(U)[2]))
    A[not_in2(rlow, end), not_in2(Int64[], end)] = U
    U=A
    
    A=fill(NaN,(newsize[2],size(V)[2]))
    A[not_in2(clow, end), not_in2(Int64[], end)] = V
    V=A
    
	if errormap
    A=fill(NaN,newsize)
    A[not_in2(rlow, end), not_in2(clow, end)] = errmap
    errmap=A
    end
	
    #@show "here"
    
    #U=DINEOF_insertNaNr(U,rlow)
    #V=DINEOF_insertNaNr(V,clow)
    
    #Filtered solution. To save memory could be optional, also full matrix return could be optional
    
    
    
    # reshape
	# Since X is not calculated maybe release array X2D instead?
	X2D=[]
    #X=reshape(X2D ,sizeperminput)
    # permutation back into original form
    #@show sortperm(perminput),perminput
    #X=permutedims(X,sortperm(perminput))
    # For U and V only dimensions related to their group and make array of arrays
    # To do have a better way to store the EOFs on their grid. I think there is one level of [] too much ??
    #@show size(X)[g1],size(X)[g2],size(U),size(V)
    UG=fill([],size(U)[2])
    VG=fill([],size(U)[2])
    #@show UG,size(UG)
    for jj=1:size(U)[2]
        UG[jj]= [reshape(U[:,jj],size(X)[g1])]
        VG[jj]= [reshape(V[:,jj],size(X)[g2])]
    end
    #@show size(UG[1][1])
	if errormap
		errmap=permutedims(reshape(errmap ,sizeperminput),sortperm(perminput))
	end
    
    return datamean,permutedims(reshape(XF2D ,sizeperminput),sortperm(perminput)),UG,S,VG,cva,cvb,errmap,musquare
    # Or return 

    
end
