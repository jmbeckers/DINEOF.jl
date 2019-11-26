"""

      
	  DINEOF_fuse!(X,XA,smoothiterations=4)


# In place fusion of X and XA. 

   X is modified in place such that where X is NaN it takes the value found in XA. To smooth the transition between the two fields weighting function can be used, with 1 in NaN locations of X propagated by box averaging


# Input: 
* `X`: an NDimensional array of original data with NaN at missing points

* `XA` : an NDimensional array which contains the analysed field (the EOF reconstruction) at all points (present and absent points). It might still contain NaN if XA is on a real grid and containes topologically masked points (or points with not enough data during reconstructions). Those NaNs will be imprinted on X.

* `smoothiterations` : number of times a box filter is applied on the weights (starting with 1 in NaN positions and 0 elsewhere)

# Output: 

   in place modification of X




"""
function DINEOF_fuse!(X,XA,smoothiterations=4)
# Takes the original data and where there is a NaN in that field, uses the XA, other wise X. 
# To smooth the transition, a weightinf field is calculated by diffusing 1 initially on the NaN positions during smoothiterations
# If XA has NaN they will overide the X value by contamination. But this is to be expected as XA is only calculated in regions with sufficient coverage
# Modifies X

filter=zeros(Float64,size(X))
filter[isnan.(X)].=1.0

# iterate smoothiterations times with 1 fixed on NaN positions and then repeated box filter
for k=1:smoothiterations
	filter=DIVAnd_filter3(filter,NaN,1)
	filter[isnan.(X)].=1.0
end
#@show filter

# On NaN positions of X, take X`

X[isnan.(X)]=XA[isnan.(X)]

# On other positions take weighted average. WARNING NaNs of XA will be imprinted on X

X[.!isnan.(X)]= (1.0 .- filter[.!isnan.(X)]  ) .* X[.!isnan.(X)] .+ filter[.!isnan.(X)] .*  XA[.!isnan.(X)]

return

end