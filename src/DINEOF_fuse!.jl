
function DINEOF_fuse!(X,XA,smoothiterations=4)
# Takes the original data and where there is a NaN in that field, uses the XA, other wise X. 
# To smooth the transition, a weightinf field is calculated by diffusing 1 initially on the NaN positions during smoothiterations
# If XA has NaN they will overide the X value by contamination. But this is to be expected as XA is only calculated in regions with sufficient coverage
# Modifies X

filter=zeros(Float64,size(X))

filter[isnan.(X)].=1.0

for k=1:smoothiterations
filter=DIVAnd_filter3(filter,NaN,1)
filter[isnan.(X)].=1.0
end
#@show filter

X[isnan.(X)]=XA[isnan.(X)]

X[.!isnan.(X)]= (1.0 .- filter[.!isnan.(X)]  ) .* X[.!isnan.(X)] .+ filter[.!isnan.(X)] .*  XA[.!isnan.(X)]

return

end