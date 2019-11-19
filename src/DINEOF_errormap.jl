"""


    errormap=DINEOF_errormap(U,S,V,musquare,missingvalues)
	
# Analysis error variance based on OI interpretation of EOF base


# Input: 


* `U,S,V` decomposition 

* `musquare` : error variance on data (including inflation for correlated data).

* `missingvalues` : array of indexes of the missing values


# Output:



* `errormap` : array of error variance of reconstruction, same dimension as USV'




"""
function DINEOF_errormap(U,S,V,musquare,missingvalues)
    # Need to add optimisation of mysquare so keep location of cross-validation points in the parameter list 
    M=size(U)[1]
    N=size(V)[1]
    NE=size(U)[2]
    @show M,N,musquare
	
	
    if N>M
        @warn("Sorry normally this function is called with N<M")
        return 1
    end
    # Loop on smaller dimension
    errmap=zeros(Float64,size(U)[1],N)
	
	qcmap=zeros(Float64,size(U)[1],N)
    for j=1:N
        # points that are present in the image j
        #@show j
        ipresent=setdiff(1:M,missingvalues[findall(x->x==j,missingvalues[:,2]),1])
		L=U*diagm(S)/sqrt(N)
        #Lp=U[ipresent,:]*diagm(S)/sqrt(N)
        #@show size(Lp)
		Lp=L[ipresent,:]
        AA=cholesky(Lp'*Lp/musquare+Matrix{Float64}(I, NE, NE))
        #BB= (U*diagm(S))*inv(AA.U)/sqrt(N)
		BB=L*inv(AA.U)
        for kk=1:NE
            errmap[:,j]=errmap[:,j]+BB[:,kk].*BB[:,kk]
        end
		
		
    end
    
	
	
    
    @show "ended"
    @show size(missingvalues),size(U),size(V),size(S),N,M,NE
    
    
    
    return errmap
    
end