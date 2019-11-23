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
function DINEOF_musquare(X,U,S,V,missingvalues,cvpoints,cvEOF;musquaresamples=[0.01,0.1,1.0,10.,100.],musquaremethod="cvpoints",maxjsamples=100)
    # version where estimator is based on all present points with one taken out approach
    # another version  use crossvalidation points
	# Possible CPU time decrease: sample columns for j only and do not all of them !!! Should be quite efficient ?????
	# TODO: add version to be close to CV estimation of error AND use diagonal decomposition
	
    M=size(U)[1]
    N=size(V)[1]
    NE=size(U)[2]
    NP=M*N-size(missingvalues)[1]
    NCV=size(cvpoints)[1]
    #@show M*N,NP,NCV
    
	jlist=randperm(N)[1:min(N,maxjsamples)]
	
	ispoints=0
	icvpoints=0
	
    OMAii=zeros(Float64,M)
    XD=zeros(Float64,M)
    thetasvalues=zeros(Float64,size(musquaresamples))
	#Esvalues=zeros(Float64,size(musquaresamples))
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
	#	Esquare=0
		ispoints=0
		icvpoints=0
		# Ok here is the ERROR TO CORRECT: reconstruction must ELIMINATE CV points 
    for j in jlist
        OMAii=ones(Float64,M)
        XD=zeros(Float64,M)
        # points that are present in the image j
        #@show j
		# Also error on points that are icv ?
		imis=[missingvalues[findall(x->x==j,missingvalues[:,2]),1]...,cvpoints[findall(x->x==j,cvpoints[:,2]),1]...]
		
        ipresent=setdiff(1:M,missingvalues[findall(x->x==j,missingvalues[:,2]),1])
		
		ipresent=setdiff(1:M,imis)
            #icv=setdiff(1:M,cvpoints[findall(x->x==j,cvpoints[:,2]),1])
			icv=cvpoints[findall(x->x==j,cvpoints[:,2]),1]
			#@show size(imis),size(icv)
		L=U*diagm(S)/sqrt(N)
        #Lp=U[ipresent,:]*diagm(S)/sqrt(N)
        #@show size(Lp)
		Lp=L[ipresent,:]
           # @show size(L),size(Lp)
            AA=cholesky(Lp'*Lp+musquare*Matrix{Float64}(I, NE, NE))
        #BB= (U*diagm(S))*inv(AA.U)/sqrt(N)
		# CHECKKKKKK with new definition of A
		invAAU=inv(AA.U)
        BB=L*invAAU
           # @show size(BB)
		# Here One could also control the error estimate and compare it to the cv error
		#for kk=1:NE
        #    errmap[:,j]=errmap[:,j]+BB[:,kk].*BB[:,kk]
        #end
		# Nope the following is the estimated error not the actual CV values
		#for kk=1:NE
		# Esquare=Esquare+sum(BB[:,kk].*BB[:,kk])  
		#end
        for kk=1:NE
            OMAii[:]=OMAii[:]-BB[:,kk].*BB[:,kk]
        end
            if minimum(OMAii[ipresent])<0 
                @show "????",minimum(OMAii),mean(OMAii)
            end
           # @show size(OMAii)
		   # Rather compare to EOF reconstructed field? but then just calculate difference everywhere ?
            #XD=BB*(BB'*X[:,j])-X[:,j]
			#@show size(V), size(V[j,:]),size(U),size(diagm(S
			# Need to put zeros in missing points !!!!!
			data=zeros(Float64,size(X[:,j]))
			data[ipresent]=X[ipresent,j]
			#filteredf=U*diagm(S)*(V[j,:])
			#XD=BB*(BB'*data)-filteredf
			
			XD=S.*(invAAU*BB'*data/sqrt(N) - V[j,:])
			#@show sum((S.*zzz).^2), sum(XD.^2)
			XDcv=BB*(BB'*data)-X[:,j]
			#@show var(XDcv),var(X[:,j]),var(data),var(data[ipresent]-filteredf[ipresent])
			
            #@show ipresent,XD[ipresent],OMAii[ipresent]
            #thetas=thetas+sum(XD[ipresent].^2 ) #./ (OMAii[ipresent]).^2)
			# FOR THETAS OPTIMIZATION POSSIBLE SINCE BOTH terms of XD areb basically BB* times something and BB'*BB is a diagonal matrix ..
			# BB=  U S inv(AA.U) /sqrt(N)
			# so calculate inv(AA.U)'*data/sqrt(N)=? zzz=inv(AA.U)*data/sqrt(N) - V[j,:] and finally sum((S.*zzz).^2)
			#
			
			thetas=thetas+sum(XD.^2)
            thetascv=thetascv+sum(XDcv[icv].^2)
            #ispoints=ispoints+size(ipresent)[1]
			ispoints=ispoints+size(XD)[1]
			icvpoints=icvpoints+size(icv)[1]

    end
        #@show thetascv,icvpoints,thetascv/icvpoints,musquare,size(cvpoints)
        thetasvalues[ii]=thetas
         thetascvvalues[ii]=thetascv
		 #Esvalues[ii]=Esquare/(N*size(jlist)[1])
    
    #@show thetasvalues
    if icvpoints>0
        thetascvvalues[ii]=(thetascv/icvpoints-cvEOF)^2/cvEOF^2
		@time mycheck=DINEOF_OIcheckCV(X,U,S,V,missingvalues,cvpoints,musquare,jlist,cvEOF)
		@show thetascvvalues[ii],mycheck
    end
    if ispoints>0
	
        thetasvalues[ii]=((thetas/ispoints)/(sum(S.^2)/size(S)[1]))
		@time mycheck=DINEOF_OIcheckX(X,U,S,V,missingvalues,cvpoints,musquare,jlist)
		@show thetasvalues[ii],mycheck
    end
    
	end
   #@show thetasvalues,thetascvvalues #,Esvalues
    
    bestin=findmin(thetasvalues)[2]
    bestincv=findmin(thetascvvalues)[2]
	if musquaremethod=="cvpoints"
	  bestin=bestincv
	end
    
    
	#@show bestin,bestincv,thetascvvalues,size(musquaresamples)[1]
	#if bestin==1
	#@show "WTD"
	#end
	
	#if bestin==size(musquaresamples)[1]
	#@show "WTD2"
	#end
	@show thetascvvalues,thetasvalues
	if bestin==1 || bestin==size(musquaresamples)[1]
	@warn("Optimal value of musquare found at edge of search region $(musquaresamples[1]) to $(musquaresamples[end])")
	end
    return musquaresamples[bestin],thetascvvalues[bestin],thetasvalues,thetascvvalues
end