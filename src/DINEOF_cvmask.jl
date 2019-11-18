"""


    cvmask=DINEOF_cvmask(X,coverage=0.01;cvmethod="Random",maxbubblesize=0.01*[size(X)...],dimensionsforcopy=[],maximumiterations=1000)

# Using the data array X, create a cross validation mask of same size. 


# Input: 
* `X`: an N-dimensional array of data where missing values are NaN

* `coverage` : fraction of points to be used for cross validation (NaN points in X are not taken into account in this coverage estimate and cvmask will never be true on a NaN point of X)

* `cvmethod`: optional keyword argument. Current options are "Random", "Bubbles" or "CopyMask". "Random" just adds single points; "Bubbles" adds ellipsoides with maximum axes length defined in maxbubblesize. "CopyMask" takes a slive of X and uses the NaNs in this slice to create a slice of cvmask in another location. The dimension along which the slice is moves/copied is defined in array dimensionsforcopy.

* `maxbubblesize`: array which for each dimension contains the maximum size of the bubble ellipsoide

* `dimensionsforcopy`: array of boolens which for each dimensions tells whether or not a mask from NaN in X can be copied.

* `maximumiterations` : saveguard for option "CopyMask" which is not ensured to reach the requested `coverage`.

# Output:

* `cvmask`: boolean array of the same size as input array X with "true" on points to be used for cross-valitation

"""
function DINEOF_cvmask(X,coverage=0.01;cvmethod="Random",maxbubblesize=0.01*[size(X)...],dimensionsforcopy=[],maximumiterations=1000)

    # Creates cross validation points for the ND-array X
    # X already has missing points where X is NaN
    # Methods to create CV point:
    # Random : just points at random positions
    # Bubbles : Ellisoides in ND with maximum size in each direction specified in maxbubblesize
    # CopyMask: copy NaN mask from a slice of X into another slice in cross validation matrix. The dimensions along which the slices are moved are
    #  put to 1 on the array dimensionsforcopy. 
    # cvmask is "true" for the points to be possibly used for cross validation
    
    M=prod(size(X))
    cvmask=fill(false,size(X))
    
    if coverage > 0.7
        @warn("Why would you try a CV with so many points?")
        coverage=0.7
    end
    
    cvpoints=Int(ceil(coverage*sum(.!isnan.(X))))
    
    cvdone=0
    #@show cvpoints
    if cvmethod=="Random"
        batchsize=Int(ceil(0.1*cvpoints))
        @show batchsize
        while cvdone<cvpoints
           
            setofpoints=mod.(rand(Int,batchsize),M).+1
            #@show setofpoints
            cvmask[setofpoints].=true
            # take out the points which have no data
            cvmask[isnan.(X)].=false
            cvdone=sum(cvmask)
            #@show cvdone
            
        end
    end
    
    if cvmethod=="Bubbles"
        #@show maxbubblesize     
        # For the moment  implemented as ellipsoides 
        #bsize=Int(ceil(sum(maxbubblesize)/size(bubblesize)[1]/2))
        #@show bsize
        Ifirst, Ilast = CartesianIndices(X)[1], CartesianIndices(X)[M]
        #@show Ifirst,Ilast
        I1 = oneunit(Ifirst)
        #@show I1
        while cvdone<cvpoints
            position=mod(rand(Int),M)+1
            #@show CartesianIndices(X)[position]
            I=CartesianIndices(X)[position]
            #@show  LinearIndices(X)[I]
            #setofpoints=[position...]
            # Now expand around this point but how?
       
            # define ellipsoide size
            bubbles=0.5*rand(Float64,ndims(X)) .* maxbubblesize
            bubbles[bubbles.<0.5].=0.5
            # Check in the cube of largest coverage
            bsize=Int(ceil(maximum(bubbles)))
            #@show bubbles
    
            for J in max(Ifirst, I-bsize*I1):min(Ilast, I+bsize*I1)
                # and check if the point false into the ellipsoide
                
                dist=0.0
                
                diffcor=I-J
                
                for kkk=1:ndims(X)
                    dist=dist+(getindex(diffcor,kkk)/bubbles[kkk])^2
                end
                if dist<1.01
                cvmask[J]=true
                end
            
        end
        
   
            
            
            
            cvmask[isnan.(X)].=false
            cvdone=sum(cvmask)
        end
        
    end
    
    if cvmethod=="CopyMask"
        # Dimension from which to copy
        groups=dimensionsforcopy

        n=ndims(X)
        #@show n
        ntimes=0
        while cvdone<cvpoints
            myf=rand(Int,n)
            #@show myf
            myt=rand(Int,n)

            for i=1:n
                myf[i]=mod(myf[i],size(X)[i])+1
                myt[i]=mod(myt[i],size(X)[i])+1
            end
            indf=[(groups[j] == 1 ? (myf[j]) : (:)) for j = 1:n]
            indt=[(groups[j] == 1 ? (myt[j]) : (:)) for j = 1:n]

            #@show indf,indt
  
            cvmask[indt...]=cvmask[indt...] .| isnan.(X[indf...])
            cvmask[isnan.(X)].=false
            cvdone=sum(cvmask)
            ntimes=ntimes+1
            if ntimes>maximumiterations
                @warn("requested cv points not reached, maybe not enought missing points to start with")
                break
            end
                
        end
    end
    
    
    return cvmask
    
end