"""



"""
function DINEOF_musquare(X,U,S,V,missingvalues,cvpoints;musquaresamples=[0.01,0.1,1.0,10.,100.],musquaremethod="cvpoints")
    # version where estimator is based on all present points with one taken out approach
    # another version  use crossvalidation points
    M=size(U)[1]
    N=size(V)[1]
    NE=size(U)[2]
    NP=M*N-size(missingvalues)[1]
    NCV=size(cvpoints)[1]
    @show M*N,NP,NCV
    
    OMAii=zeros(Float64,M)
    XD=zeros(Float64,M)
    thetasvalues=zeros(Float64,size(musquaresamples))
    thetascvvalues=zeros(Float64,size(musquaresamples))
    if N>M
        @warn("Sorry normally this function is called with N<M")
        return 1
    end
    ii=0
    for musquare in musquaresamples
        ii=ii+1
        #@show musquare
        thetas=0   
        thetascv=0
    for j=1:N
        OMAii=ones(Float64,M)
        XD=zeros(Float64,M)
        # points that are present in the image j
        #@show j
        ipresent=setdiff(1:M,missingvalues[findall(x->x==j,missingvalues[:,2]),1])
            icv=setdiff(1:M,cvpoints[findall(x->x==j,cvpoints[:,2]),1])
		L=U*diagm(S)/sqrt(N)
        #Lp=U[ipresent,:]*diagm(S)/sqrt(N)
        #@show size(Lp)
		Lp=L[ipresent,:]
           # @show size(L),size(Lp)
            AA=cholesky(Lp'*Lp+musquare*Matrix{Float64}(I, NE, NE))
        #BB= (U*diagm(S))*inv(AA.U)/sqrt(N)
        BB=L*inv(AA.U)
           # @show size(BB)
        for kk=1:NE
            OMAii[:]=OMAii[:]-BB[:,kk].*BB[:,kk]
        end
            if minimum(OMAii)<0 
                @show "????"
            end
           # @show size(OMAii)
            XD=BB*(BB'*X[:,j])-X[:,j]
            #@show ipresent,XD[ipresent],OMAii[ipresent]
            thetas=thetas+sum(XD[ipresent].^2 ./ (OMAii[ipresent]).^2)
            thetascv=thetascv+sum(XD[icv].^2)
            

    end
        #@show thetas
        thetasvalues[ii]=thetas/NP
         thetascvvalues[ii]=thetascv
    end
    #@show thetasvalues
    if NCV>0
        thetascvvalues=thetascvvalues/NCV
    end
    
   # @show thetasvalues,thetascvvalues
    
    bestin=findmin(thetasvalues)[2]
    bestincv=findmin(thetascvvalues)[2]
	if musquaremethod=="cvpoints"
	  bestin=bestincv
	end
    
    
	@show bestin,bestincv
    return musquaresamples[bestin],thetascvvalues[bestin],thetasvalues,thetascvvalues
end