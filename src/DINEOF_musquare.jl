"""


    musquare,thetavar,thetaoneoutvar,thetacvvar=DINEOF_musquare(X,U,S,V,missingvalues,cvpoints;musquaresamples=[0.01,0.1,1.0,10.,100.],musquaremethod="cvpoints",maxjsamples=100)

# Estimator of mu_eff^2 


# Input: 
* `X`: a two-dimensional array of size MxN with N<=M and no missing values.

* `U,S,V` decomposition 

* `missingvalues` : list of indexes of missing values in original array 

* `cvpoints` : list of indexes used for crossvalidation

* `musquaresamples` : array of values to be tested for musquare

* `musquaremethod` : either from cvpoints ("cvpoints") or "LeaveOneOut" approach. Both methods are actually calculated and only the final value of musquare selected via this parameter

* `maxjsamples` : maximum number of (random) columns to be used for the optimization


# Output:

* `musquare`: requested value for musquare

* `thetavar` : reconstruction error variance with the optimal value of musquare if EOFs used as OI base

* `thetaoneoutvar` : array of error variance of reconstructions using the leave-one-out lemma

* `thetacvvar` : array of error variance of reconstructions using the cv points
	

"""
function DINEOF_musquare(X,U,S,V,missingvalues,cvpoints;musquaresamples=[0.01,0.1,1.0,10.,100.],musquaremethod="cvpoints",maxjsamples=100)
    # version where estimator is based on all present points with one taken out approach
    # another version  use crossvalidation points
	# Possible CPU time decrease: sample columns for j only and do not all of them !!! Should be quite efficient ?????
    M=size(U)[1]
    N=size(V)[1]
    NE=size(U)[2]
    NP=M*N-size(missingvalues)[1]
    NCV=size(cvpoints)[1]
    @show M*N,NP,NCV
    
	jlist=randperm(N)[1:min(N,maxjsamples)]
	
	ispoints=0
	icvpoints=0
	
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
		ispoints=0
		icvpoints=0
		
    for j in jlist
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
            ispoints=ispoints+sum(ipresent)
			icvpoints=icvpoints+sum(icv)

    end
        #@show thetas
        thetasvalues[ii]=thetas
         thetascvvalues[ii]=thetascv
    end
    #@show thetasvalues
    if icvpoints>0
        thetascvvalues=thetascvvalues/icvpoints
    end
    if ispoints>0
	
        thetasvalues=thetasvalues/ispoints
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