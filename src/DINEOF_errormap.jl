function DINEOF_errormap(U,S,V,musquare,missingvalues)
    # Need to add optimisation of mysquare so keep location of cross-validation points in the parameter list 
    M=size(U)[1]
    N=size(V)[1]
    NE=size(U)[2]
    @show M,N
    if N>M
        @warn("Sorry normally this function is called with N<M")
        return 1
    end
    # Loop on smaller dimension
    errmap=zeros(Float64,size(U)[1],N)
    for j=1:N
        # points that are present in the image j
        #@show j
        ipresent=setdiff(1:M,missingvalues[findall(x->x==j,missingvalues[:,2]),1])
        Lp=U[ipresent,:]*diagm(S)/sqrt(N)
        #@show size(Lp)
        AA=cholesky(Lp'*Lp/musquare+Matrix{Float64}(I, NE, NE))
        BB= (U*diagm(S))*inv(AA.U)/sqrt(N)
        for kk=1:NE
            errmap[:,j]=errmap[:,j]+BB[:,kk].*BB[:,kk]
        end
    end
    
    
    @show "ended"
    @show size(missingvalues),size(U),size(V),size(S),N,M,NE
    
    
    
    return errmap
    
end