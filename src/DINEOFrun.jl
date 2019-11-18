function DINEOFrun(X,whichgroups;minimumcoverage=(0.1, 0.1),cvmask="Automatic",cvfraction=0.01,cvmethod="Random")
    # X ND array with NaN at missing locations
    # whichgroups: array of 1 and 2 indicating which dimensions are collapsed together. eg [1 2 1 2] regroups dimensions 
    # 1 and 3 into i direction and 2 with 4 into j direction
    
    
    if cvmask=="Automatic"
        
        cvmask=DINEOF_cvmask(X,cvfraction;cvmethod=cvmethod)
        
    else
        @show size(cvmask),size(X)
        if size(cvmask)!=size(X)
            @warn("Mask does not correspond")
            return
        end
    end
    
    # Do not use cross validation in already missing points
    cvmask[isnan.(X)].=false
    
    if size(whichgroups)[1]!=ndims(X)
        @warn("Incompatible dimensions")
    end
    @show ndims(X),size(whichgroups)[1]
    #
    g1=findall(x->x==1,whichgroups)
    g2=findall(x->x==2,whichgroups)
    if size(g1)[1]<1 || size(g2)[1]<1 || size(g1)[1]+size(g2)[1] != ndims(X)
        @warn("Incorrect whichgroups")
        @show whichgroups
    end
    perminput=[g1...,g2...]
    @show perminput
    newsize=(prod(size(X)[g1]),prod(size(X)[g2]))
    @show size(X)
    @show size(X)[perminput]
    sizeperminput=size(X)[perminput]
    X2D=reshape(permutedims(X,perminput),newsize)
    cv2D=reshape(permutedims(cvmask,perminput),newsize)
    @show size(X2D)
    
    # now take out lines and columns with less the specified coverage
    # 
    #
    @show size(isnan.(X2D))
    @show sum(cv2D)
     
    ccov=dropdims(sum(.!isnan.(X2D),dims=1)',dims=2) .- dropdims(sum(cv2D,dims=1)',dims=2)
    lcov=dropdims(sum(.!isnan.(X2D),dims=2),dims=2) .- dropdims(sum(cv2D,dims=2),dims=2)
    #@show ccov,size(ccov)
    #@show lcov,size(lcov)
    rlow=findall(x->x<newsize[2]*minimumcoverage[1],lcov)
    clow=findall(x->x<newsize[1]*minimumcoverage[2],ccov)
    
    @show rlow,clow
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
    
    # Ok, now we have the working matrix. Do some stats on it:
    #  
    varmatrix=var(X2D[.!isnan.(X2D)])
    missingpointsforvar=deepcopy(.!isnan.(X2D))
    @show size(missingpointsforvar),size(X2D)
    meanmatrix=mean(X2D[.!isnan.(X2D)])
    @show varmatrix,meanmatrix
    
    NM=sum(isnan.(X2D))
    NMCV=sum(cv2D)
    @show NM,NMCV
    
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
    
    
    
    
    U,S,V,cva,cvb,errmap=DINEOF_svds!(X2D,missingvalues,cvpoints)
    # now roll back
    #@show X2D
    @show size(U),size(S),size(V)
    # Transpose back if that was done (also SVD...)
    if transposed
        TT=deepcopy(U)
        U=deepcopy(V)
        V=deepcopy(TT)
        TT=[]
        X2D=permutedims(X2D,[2,1])
        errmap=permutedims(errmap,[2,1])
    end
    
    # put back lines and columns (also SVD)
    
    # Not so easy ....
    # Lines rlow into U columns clow into V' so into lines of V
    # Check if Alex version is quicker
    # 
    
    #X2D=DINEOF_insertNaNr(X2D,rlow)
    #X2D=DINEOF_insertNaNc(X2D,clow)
    
    @show size(U),size(V),size(S)
    XF2D=U*diagm(S)*V'
    @show size(XF2D)
    @show size(missingpointsforvar)
    varmatrixf=var(XF2D[missingpointsforvar])
    
    @show varmatrix,varmatrixf,varmatrix-varmatrixf
    
    @show "back"
    A=fill(NaN,newsize)
    A[not_in2(rlow, end), not_in2(clow, end)] = X2D
    X2D=A
    
    A=fill(NaN,newsize)
    A[not_in2(rlow, end), not_in2(clow, end)] = XF2D
    XF2D=A
    
    A=fill(NaN,(newsize[1],size(U)[2]))
    A[not_in2(rlow, end), not_in2(Int64[], end)] = U
    U=A
    
    A=fill(NaN,(newsize[2],size(V)[2]))
    A[not_in2(clow, end), not_in2(Int64[], end)] = V
    V=A
    
    A=fill(NaN,newsize)
    A[not_in2(rlow, end), not_in2(clow, end)] = errmap
    errmap=A
    
    @show "here"
    
    #U=DINEOF_insertNaNr(U,rlow)
    #V=DINEOF_insertNaNr(V,clow)
    
    #Filtered solution. To save memory could be optional, also full matrix return could be optional
    
    
    
    # reshape
    X=reshape(X2D ,sizeperminput)
    # permutation back into original form
    @show sortperm(perminput),perminput
    X=permutedims(X,sortperm(perminput))
    # For U and V only dimensions related to their group and make array of arrays
    # To do have a better way to store the EOFs on their grid. I think there is one level of [] too much ??
    @show size(X)[g1],size(X)[g2],size(U),size(V)
    UG=fill([],size(U)[2])
    VG=fill([],size(U)[2])
    @show UG,size(UG)
    for jj=1:size(U)[2]
        UG[jj]= [reshape(U[:,jj],size(X)[g1])]
        VG[jj]= [reshape(V[:,jj],size(X)[g2])]
    end
    @show size(UG[1][1])
    
    return X,permutedims(reshape(XF2D ,sizeperminput),sortperm(perminput)),UG,S,VG,cva,cvb,permutedims(reshape(errmap ,sizeperminput),sortperm(perminput))
    # Or return 

    
end
