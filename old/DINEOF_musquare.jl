"""


    musquare,thetavar,thetaoneoutvar,thetacvvar=DINEOF_musquare(X,U,S,V,missingvalues,cvpoints,cvEOF;musquarer=[.01,1.].*var(X),musquarew=[1.0,1.0],maxjsamples=100)

# Estimator of mu_eff^2 


# Input: 
* `X`: a two-dimensional array of size MxN with N<=M and no missing values.

* `U,S,V` decomposition 

* `missingvalues` : list of indexes of missing values in original array 

* `cvpoints` : list of indexes used for crossvalidation

* `cvEOF` : cross validator from EOF decompisition. musquare tries to find OI with CV estimator to be close to that one if musquarew[2] is not zero

* `musquarer` : range in which to search for optimal musquare

* `musquarew` : weighting of the two criteria: OI reconstruction close to EOF decomposition or OI CV value close to EOF CV value

* `maxjsamples` : maximum number of (random) columns to be used for the optimization


# Output:

* `musquare`: requested value for musquare

* `cvOI`: cross validator value when OI is used

	

"""
function DINEOF_musquare(X,U,S,V,missingvalues,cvpoints,cvEOF;musquarer=[.01,1.].*var(X),musquarew=[1.0,1.0],maxjsamples=100)
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
		@time thetascvvalues[ii],cvOI=DINEOF_OIcheckCV(X,U,S,V,missingvalues,cvpoints,musquare,jlist,cvEOF)
		#@show thetascvvalues[ii],mycheck,thetascv,tcheck,icvpoints,icheck
    end
    if ispoints>0
	
        thetasvalues[ii]=((thetas/ispoints)/(sum(S.^2)/size(S)[1]))
		@time mycheck,aerr=DINEOF_OIcheckX(X,U,S,V,missingvalues,musquare,jlist)
		#@show thetasvalues[ii],mycheck
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
	w1opt=1
	w2opt=1
	# OTher method
	#if w2opt==0
	#    @show "OI only"
		myfun1(x)=DINEOF_OIcheckX(X,U,S,V,missingvalues,x,jlist)[1]
	#end
	#tutu=Optim.optimize( myfun ,musquaresamples[1], musquaresamples[end],rel_tol=0.01,abs_tol=1.E-7,iterations=15)
	#if w1opt==0
	#    @show "CV only"
		myfun2(x)=DINEOF_OIcheckCV(X,U,S,V,missingvalues,cvpoints,x,jlist,cvEOF)[1]
	#end
		
	
	# try combination weighted by taking into accound N*M and size(cvpoints)  since estimators themself are more or less well defined depending on amount of data?
	#if w1opt>0 && w2opt>0
	#@show "Both"
	myfunb(x) = w1opt* prod(size(X))*DINEOF_OIcheckX(X,U,S,V,missingvalues,x,jlist)[1]+
	 w2opt*size(cvpoints)[1]*DINEOF_OIcheckCV(X,U,S,V,missingvalues,cvpoints,x,jlist,cvEOF)[1]
	#end
	tutu=[]
	MF1=0
	MF2=0
	musquareopt=0.
	cvOI=0
	if w2opt==0
        @show "OI only"
	    tutu=Optim.optimize( myfun1 ,musquaresamples[1], musquaresamples[end],rel_tol=0.01,abs_tol=1.E-7,iterations=15)
		musquareopt=Optim.minimizer(tutu)
		MF1=DINEOF_OIcheckX(X,U,S,V,missingvalues,musquareopt,jlist)[1]
	end
	if w2opt==0
        @show "CV only"
	    tutu=Optim.optimize( myfun2 ,musquaresamples[1], musquaresamples[end],rel_tol=0.01,abs_tol=1.E-7,iterations=15)
		musquareopt=Optim.minimizer(tutu)
		MF2,cvOI=DINEOF_OIcheckCV(X,U,S,V,missingvalues,cvpoints,musquareopt,jlist,cvEOF)
	end
	if w1opt>0 && w2opt>0
		@show "Both"
		tutu=Optim.optimize( myfunb ,musquaresamples[1], musquaresamples[end],rel_tol=0.01,abs_tol=1.E-7,iterations=15)
		musquareopt=Optim.minimizer(tutu)
		MF1=DINEOF_OIcheckX(X,U,S,V,missingvalues,musquareopt,jlist)[1]
		MF2,cvOI=DINEOF_OIcheckCV(X,U,S,V,missingvalues,cvpoints,musquareopt,jlist,cvEOF)
	end
	
	@show tutu
	
	
	MF1J=prod(size(X))*MF1
	MF2J=MF2*(size(cvpoints)[1])
	
	@show MF1,MF2,cvEOF,musquareopt,cvOI
	#@show tutu2
	@show cvEOF
	println("CV estimator from EOF $cvEOF is now $cvOI if OI is used")
	println("Optimal musquare is $musquareopt")
	println("Relative error on reconstruction $MF1, relative error on CV estimator $MF2")
	println("The two criteria to compare OI and EOF are: reconstruction $MF1J, closest CV $MF2J")
	return musquareopt,cvOI,MF1,MF2
    #return musquaresamples[bestin],thetascvvalues[bestin],thetasvalues,thetascvvalues
	
end