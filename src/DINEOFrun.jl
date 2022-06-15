"""


    XA,offset,U,S,V,cvEOF,cvarray,errmap,musquare=DINEOFrun(X,whichgroups=[ ones(Int32,ndims(X)-1)...,2];
	minimumcoverage=(0.1, 0.1),
	cvmask="Automatic",
	cvfraction=0.01,
	cvmethod="Random",
	maxbubblesize=0.01*[size(X)...],
	dimensionsforcopy=[zeros(Int32,ndims(X)-1)...,1],
	errormap=true,
	musquare=0,
	restart=[],
    keepsvdcvrestart=true,
	eofmax=size(X)[2]-1,
	eofstart=1,
	dineofmaxiter=10,
	dineoftol=0.001,
	svdmeth="svd",
	svdtol=0.000001,
	filter="None",
	filterintensity=1.0,
	filterrepetitions=1)



Provides a DINEOF reconstruction of an N-dimensional array `X`. Missing points are identified with NaN values. Since the input arrays can have more than 2 dimensions, you have to specify which dimensions need to be collapsed into one. To do so the array `whichgroups` specifies to which collapsed dimensions (1 or 2) each dimension is collapsed into. The output is the filtered field. If it is to be merged with the original data, you can use DINEOF_fuse.

# Input:

* `X`: The N-dimensional array containing the data. Missing points are NaN

* `whichgroups`: array of "1" or "2" of size ndims(X). for each dimension it tells if it goes into dimensions 1 or 2 for the SVD decomposition. Default is last dimension is group 2 only

# Optional keyword inputs with their defaults:

* `minimumcoverage=(0.1, 0.1)` : The minimum coverage in each dimension for the collapsed matrix. If the coverage is below the threshold, the line or colums is taken out

* `cvmask="Automatic"` : You can provide your own cross-validation mask. In that case cvmask is a boolean array of the same size as X with "true" on points for crossvalidation. If "Automatic", DINEOF will create the mask based on the next parameters

* `cvfraction=0.01` : fraction of points to be used for cross validation (fraction is with respect to valid points)
	
* `cvmethod="Random"` : method to create the cross validation mask. "Random" "Bubbles" or "CopyMask": If "Bubbles" are created the maximum size in each direction is specified in maxbubblesize. If "CopyMask" is used, dimensionsforcopy specifies along which dimensions tha NaN pattern of X can be copied to create a mask

* `maxbubblesize=0.01*[size(X)...]` : maximum size of the bubbles in each direction 
	
*	`dimensionsforcopy=[zeros(Int32,ndims(X)-1),1]`: array of 0 and 1. For each direction is indicates if one can move a mask of NaNs from X along that direction.

*	`errormap=true` : if false, error map returned is []

*	`musquare=0` : You can provide your own estimate of musquare to be used for OI error map calculations. If 0, DINEOF will do the estimate

*	`restart=[]` : You can provide an array of the same size of X to fill in the first guess in the missing points. If not provided, the matrix is filled randomly with a variance of the present data

*    `keepsvdcvrestart=true` : goes back to the best estimate of the reconstruction during final EOF decomposition

*	`eofmax=size(X)[2]-1`: maximum number of EOFs 

*	`eofstart=1` : number of EOFs to start with in the search of the optimum. Can be larger than 1, particularly if you had a good restart matrix

*	`dineofmaxiter=10` : Maximum Number of iterations  USV=X, X=fillfrom(USV) 

*	`dineoftol=0.001` : relative change during iterations below which one stops

*	`svdmeth="svd"` : work with SVD or with eigenvalues of X'X ("eig")

*	`svdtol=0.000001` : tolerance during svd decomposition of filled matrix (svds or eig)

*	`filter="None"` : filter to be applied to dimension 2. Note that this is the dimension in the innermost svd decomposition which might be dimension 1 of the outermost call since a transpose is performed if M<N

*	`filterintensity=1.0` : filter intensity

*	`filterrepetitions=1` : filter repetitions



# Output:
* `XA` : the analysed (filtered) data filled in in places where enough data where available (see coverage parameter) including the offset: XA=offset+U S V'

* `offset` : the value that was subtracted from the original data to center them

* `U` : array of U arrays. To get access to mode 2: U[2][1]

* `S`: array of singular values. use diagm(S) if you want to work with matrices

* `V` : array of V arrays: To get mode 3: V[3][1]

* `cvEOF` : cross validation estimator (variance of misfit at cross validation points)

* `errmap` : array of the same structure as X containing the error variance estimate of the reconstruction XA

* `musquare` : the value of mu^2 used for the OI interpretation leading to error maps

"""
function DINEOFrun(X,whichgroups=[ones(Int32,ndims(X)-1)...,2];
	minimumcoverage=(0.1, 0.1),
	cvmask="Automatic",
	cvfraction=0.01,
	cvmethod="Random",
	maxbubblesize=max.(20,0.01*[size(X)...]),
	dimensionsforcopy=[zeros(Int32,ndims(X)-1)...,1],
	errormap=true,
	musquare=0,
	restart=[],
    keepsvdcvrestart=true,
	eofmax=size(X)[2]-1,
	eofstart=1,
	dineofmaxiter=15,
	dineoftol=0.005,
	svdmeth="svd",
	svdtol=0.000001,
	filter="None",
	filterintensity=1.0,
	filterrepetitions=1)
	
	
	
    # X ND array with NaN at missing locations
    # whichgroups: array of 1 and 2 indicating which dimensions are collapsed together. eg [1 2 1 2] regroups dimensions 
    # 1 and 3 into i direction and 2 with 4 into j direction
    
	
	#For restart, transform restart matrix as X and at the end copy into X but test if necessary. Maybe drop deallocate via =[]
	datamean=mean(X[.!isnan.(X)])
	datavar=var(X[.!isnan.(X)])
	println("Raw data variance and mean: $datavar and $datamean")
	println("Number of missing points (including possible masks): $(sum(isnan.(X))) out of $(prod(size(X)))")
	X=X.-datamean
	#  Also deal with restart values SUBRATCT MEAN THERE TOO
	restart2D=[]
	if restart!=[]
		restart=restart.-datamean
		if size(restart)!=size(X)
            @error("Size for restart matrix and size of X do not correspond")
			@show size(restart),size(X)
            return
        end
	end
	###################
	
	
	errmap=[]
	
    
    if cvmask=="Automatic"
        
        cvmask=DINEOF_cvmask(X,cvfraction;cvmethod=cvmethod,dimensionsforcopy=dimensionsforcopy,maxbubblesize=maxbubblesize)
        
    else
        #@show size(cvmask),size(X)
        if size(cvmask)!=size(X)
            @error("Mask for cross validation and size of X do not correspond")
			@show size(X),size(cvmask)
            return
        end
		println("Using cvmask provided")
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
	if restart!=[]
	    println("Using restart matrix")
		restart2D=reshape(permutedims(restart,perminput),newsize)
		@show mean(restart2D[.!isnan.(restart2D)])
	end
	
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
	#NEED TO CHECK
    rlow=findall(x->x<newsize[2]*minimumcoverage[1],lcov)
    clow=findall(x->x<newsize[1]*minimumcoverage[2],ccov)
    
    #@show rlow,clow
    # Usefull functions: take out rows and colums with low coverage and when finished later put NaNs back
    not_in(inds::AbstractVector{Int}, n::Int)::Vector{Int} = setdiff(1:n, [ i for i in inds ])
    not_in2(inds::AbstractVector{Int}, n::Int)::Vector{Int} = setdiff(1:n,inds)
  
    X2D=X2D[not_in(rlow, size(X2D,1)), not_in(clow, size(X2D,2))]
    cv2D=cv2D[not_in(rlow, size(cv2D,1)), not_in(clow, size(cv2D,2))]
    if restart!=[]
		restart2D=restart2D[not_in(rlow, size(restart2D,1)), not_in(clow, size(restart2D,2))]
		@show mean(restart2D[.!isnan.(restart2D)])
	end
    
    # now transpose matrix if necessary
    transposed=false
    if size(X2D)[1]<size(X2D)[2]
        transposed=true
        X2D=permutedims(X2D,[2,1])
        cv2D=permutedims(cv2D,[2,1])
		if restart!=[]
		restart2D=permutedims(restart2D,[2,1])
		end
    end
    #######################################################################
	# Some more work here to have info for users
    # Ok, now we have the working matrix. Do some stats on it:
    #
    varmatrix=var(X2D[.!isnan.(X2D)])
    
	#missingpointsforvar=deepcopy(.!isnan.(X2D))
    #@show size(missingpointsforvar),size(X2D)
    meanmatrix=mean(X2D[.!isnan.(X2D)])
    #@show varmatrix,meanmatrix
	# Maybe better subtract also meanmatrix ????
	
	X2D=X2D .- meanmatrix
	
	if restart!=[]
	    @show mean(restart2D[.!isnan.(restart2D)])
		restart2D=restart2D.-meanmatrix
		@show mean(restart2D[.!isnan.(restart2D)])
	end
    
    NM=sum(isnan.(X2D))
    NMCV=sum(cv2D)
    #@show NM,NMCV
	println("Number of data points after elimination of low coverage regions is $(prod(size(X2D))-NM) and cv fraction $(NMCV/(prod(size(X2D))-NM))")
	
    ##############################
	
	
    # Now create index arrays pointing to the points in the matrix with NaN
    # Probably not very elegant but hey...
    missingvalues=zeros(Int,(NM,2))
    icount=0
	meanmiss=0
    for j=1:size(X2D)[2]
        for i=1:size(X2D)[1]
            if isnan(X2D[i,j])
                icount=icount+1
                missingvalues[icount,1]=i
                missingvalues[icount,2]=j
                # Put a random value with average variance of present data so that we can keep an eye on how total variane
                # behaves.
				if restart==[] || isnan(restart2D[i,j])
                X2D[i,j]=sqrt(varmatrix)*randn()
				else
				X2D[i,j]=deepcopy(restart2D[i,j])
				end
				meanmiss=meanmiss+X2D[i,j]
            end
        end
    end
	meanmiss=meanmiss/icount
    restart2D=[]
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
    
	@show mean(X2D),meanmatrix,meanmiss,datamean
    
    
    # NEED TO ADD OPTIONAL PARAMETERS ...
    U,S,V,cva,cvb,musquareestimate=DINEOF_svds!(X2D,missingvalues,cvpoints;
	keeprestart=keepsvdcvrestart,
	ncmax=eofmax,
	istart=eofstart,
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
		println("Mean error variance of reconstruction: $(mean(errmap)) ")
	end
	
	#
    # now roll back
    #@show X2D
    #@show size(U),size(S),size(V)
    # Transpose back if that was done (also SVD...) NEED TO CHECK IF V and be used as U without problem??? (Adjoint definition)
    if transposed
	println("Internally final matrix was transposed for SVD analysis")
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
    #varmatrixf=var(XF2D[missingpointsforvar])
    
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
    #UG=fill([],size(U)[2])
    #VG=fill([],size(U)[2])
    #@show UG,size(UG)
    #for jj=1:size(U)[2]
    #    UG[jj]= [reshape(U[:,jj],size(X)[g1])]
    #    VG[jj]= [reshape(V[:,jj],size(X)[g2])]
    #end
    #@show size(UG[1][1])
	U=reshape(U,(size(X)[g1]...,size(U)[2]))
	V=reshape(V,(size(X)[g2]...,size(V)[2]))
	#@show size(U),size(UG[1][1])
	if errormap
		errmap=permutedims(reshape(errmap ,sizeperminput),sortperm(perminput))
	end
    
    return permutedims(reshape(XF2D ,sizeperminput),sortperm(perminput)).+datamean.+meanmatrix,datamean+meanmatrix,U,S,V,cva,cvb,errmap,musquare
    # Or return 

    
end


