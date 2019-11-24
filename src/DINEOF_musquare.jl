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


	
    M=size(U)[1]
    N=size(V)[1]
    
    
	jlist=randperm(N)[1:min(N,maxjsamples)]
	
	
    if N>M
        @warn("Sorry normally this function is called with N<M")
        return 1
    end
    w1opt=musquarew[1]
	w2opt=musquarew[2]
	
		myfun1(x)=DINEOF_OIcheckX(X,U,S,V,missingvalues,x,jlist)[1]
	
		myfun2(x)=DINEOF_OIcheckCV(X,U,S,V,missingvalues,cvpoints,x,jlist,cvEOF)[1]
	
		myfunb(x) = w1opt* prod(size(X))*DINEOF_OIcheckX(X,U,S,V,missingvalues,x,jlist)[1]+
			w2opt*size(cvpoints)[1]*DINEOF_OIcheckCV(X,U,S,V,missingvalues,cvpoints,x,jlist,cvEOF)[1]

	tutu=[]
	MF1=0
	MF2=0
	musquareopt=0.
	cvOI=0
	if w2opt==0
        println("OI-EOF error only")
	    tutu=Optim.optimize( myfun1 ,musquarer[1], musquarer[end],rel_tol=0.01,abs_tol=1.E-7,iterations=15)
		musquareopt=Optim.minimizer(tutu)
		MF1=DINEOF_OIcheckX(X,U,S,V,missingvalues,musquareopt,jlist)[1]
	end
	if w1opt==0
        println("Best CV coherence only")
	    tutu=Optim.optimize( myfun2 ,musquarer[1], musquarer[end],rel_tol=0.01,abs_tol=1.E-7,iterations=15)
		musquareopt=Optim.minimizer(tutu)
		MF2,cvOI=DINEOF_OIcheckCV(X,U,S,V,missingvalues,cvpoints,musquareopt,jlist,cvEOF)
	end
	if w1opt>0 && w2opt>0
		println("Both OI-EOF error and CV coherence")
		tutu=Optim.optimize( myfunb ,musquarer[1], musquarer[end],rel_tol=0.01,abs_tol=1.E-7,iterations=15)
		musquareopt=Optim.minimizer(tutu)
		MF1=DINEOF_OIcheckX(X,U,S,V,missingvalues,musquareopt,jlist)[1]
		MF2,cvOI=DINEOF_OIcheckCV(X,U,S,V,missingvalues,cvpoints,musquareopt,jlist,cvEOF)
	end
	
	@show tutu
	
	
	MF1J=prod(size(X))*MF1
	MF2J=MF2*(size(cvpoints)[1])
	
	#@show MF1,MF2,cvEOF,musquareopt,cvOI
	#@show tutu2
	#@show cvEOF
	println("CV estimator from EOF $cvEOF is now $cvOI if OI is used")
	println("Optimal musquare is $musquareopt")
	println("Relative error on reconstruction $MF1, relative error on CV estimator $MF2")
	println("The two criteria to compare OI and EOF are: reconstruction $MF1J, closest CV $MF2J")
	return musquareopt,cvOI,MF1,MF2

	
end