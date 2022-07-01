"""


    U,S,V,ncon,nit=DINEOF_svd(X,nele,svdmeth="svd",svdtol=0.00001;tfilter="None",filterintensity=1.0,filterrepetitions=1)

# Lowest level interface to SVD decomposition


# Input: 
* `X`: a two-dimensional array of size MxN with N<=M and no missing values.

* `nele`: the number of singular vectors and values asked for (<=N)

* `svdmeth`: if "svd" uses SVD decomposition, if "eig" uses eigenvalue approach on X'X. Results should be identical but CPU demand could be different

* `svdtol`: tolerance for convergence criteria during SVD decomposition

* `tfilter`: Filter to be applied to second dimension. For the moment implementation of 

   "None" : no filter applied
   
   "vpmf" : Very poor mans filter: applies a postprocessing filter to V and adapts U and S to get an orthonormal base. It is NOT the best approximation of X or the filtered X
   
   "pmf" : poor mans filter: applies a preprocessing filter on X (or X'X which mathematically is identical). USV' is the best approximation to the filtered X

For the two filters a diffusion-like filter is applied. filterintensity=1 is the strongest filter that should be used and amounts to a 0.25,0.5,0.25 filter stencil.
filterrepetitions defines how many diffusion-like steps are used (so how wide the filter stencil becomes). Presently not variable distance between images is implemented.



# Output:

* `U,S,V`: the SVD decomposition such that in the infiltered version USV' is the best approximation of X with nele singular vectors.

* `ncon` : number of converged singular values, should be equal to nele

* `nit` : number of iterations needed to find the singular values (can be used to adapt tolerance svdtol)

"""
function DINEOF_svd(X,nele,svdmeth="svd",svdtol=0.00001;tfilter="None",filterintensity=1.0,filterrepetitions=1)
    #X=deepcopy(Xin)
    
    if svdmeth=="svd"
        if tfilter=="pmf"
        #@show "poor man filter"
        n=size(X)[2]
		# Simple filter not taking into acount different distances between columns of X*FF
		# The filter is symmetric and conserves the mean
        FF=Tridiagonal(0.25*filterintensity*ones(Float64,n-1),(1.0-filterintensity*0.5)*ones(Float64,n),0.25*filterintensity*ones(Float64,n-1))
        FF[1,1]=1.0-filterintensity*0.25
        FF[end,end]=1.0-filterintensity*0.25
		# Filter X
		#b = similar(X, size(X,2)) # pre-allocate array for storing 

        #for jjj=1:filterrepetitions
		# TODO: make in place A_mul_B! multiplication to avoid allocations ??? or play with views etc ?
		# This seems to be already much faster then X=X*FF but it modifies X????
		#for i = 1:size(X, 1)
		#  b=FF*deepcopy(X[i,:])
        #  X[i,:] = b'
        #end
		
           X=X*FF^filterrepetitions
        
        end
        nel=min(nele,size(X)[2])
            SV,ncon,nit=svds(X;nsv=nel,tol=svdtol)[1:3]
              SVS=SV.S
              SVV=SV.V
             SVU=SV.U
			
			 # SVU, SVS, SVV = tsvd(X, nele)
			 # ncon=nele
			 # nit=9999
			  
    end
    if svdmeth=="eig"
        BB=X'*X
        if tfilter=="pmf"
        #@show "poor man filter"
        n=size(X)[2]
		# Same filter as for svd application but here applied to the covariance matrix
        FF=Tridiagonal(0.25*filterintensity*ones(Float64,n-1),(1.0-filterintensity*0.5)*ones(Float64,n),0.25*filterintensity*ones(Float64,n-1))
        FF[1,1]=1.0-filterintensity*0.25
        FF[end,end]=1.0-filterintensity*0.25
        #for jjj=1:filterrepetitions
		FFT=FF^filterrepetitions
        BB=FFT*BB*FFT
        #end
        end
        
        
        SVS,SVV,ncon,nit=eigs(BB;nev=nele,tol=svdtol)[1:4]
        # Should not be necessary but rounding can lead to problems?
        SVV=real.(SVV)
        SVS=sqrt.(abs.(real.(SVS)))
        SVU=X*SVV*diagm(1.0 ./SVS)
    end
    
    if tfilter=="vpmf"
        #@show "Very poor man filter"
        n=size(X)[2]
		# Same filter as for pmf but applied a posteriori to the singular vectors V
        FF=Tridiagonal(0.25*filterintensity*ones(Float64,n-1),(1.0-filterintensity*0.5)*ones(Float64,n),0.25*filterintensity*ones(Float64,n-1))
        FF[1,1]=1.0-filterintensity*0.25
        FF[end,end]=1.0-filterintensity*0.25
        #for jjj=1:filterrepetitions
        SVV=FF^filterrepetitions*SVV
        #end
        #DINEOF_GSORTHO!(SVV)
		# Now reorthonormalize the singular vectors
        SVV=Array(qr(SVV).Q)
		# Project onto U space
        SVU=X*SVV
        #DINEOF_GSORTHO!(SVU)
		# Orthonormalize those vectors
        SVU=Array(qr(SVU).Q)
        #SVS=diag(SVU'*X*SVV)
        #SVSO=deepcopy(SVS)
		# Finally take the diagonal of U'*X*V (which is not diagonal here anymore)
        for k=1:size(SVS)[1]
            SVS[k]=SVU[:,k]'*X*SVV[:,k]
			if SVS[k]<0.0
			  
			  SVS[k]=-SVS[k]
			  SVU[:,k] .= -SVU[:,k]
			end
        end
        #@show (SVSO.-SVS)./SVSO
    end
    return SVU,SVS,SVV,ncon,nit
end