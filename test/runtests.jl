using DINEOF
using JLD2
using Statistics
using FileIO
using Random
using Test

#=
# the small example is created using the example file from
https://github.com/aida-alvera/DINEOF/blob/master/SmallExample/seacoos2005.avhrr

  Random.seed!(1234); XA,offset,U,S,V,cvEOF,cvarray,errmap,musquare = @time DINEOF.DINEOFrun(X; errormap = false);
  FileIO.save("testdata.jld2","X",X,"XA",XA,"cvarray",cvarray);
=#

SmallExampleFile = joinpath(dirname(@__FILE__),"SmallExample.jld2")

if !isfile(SmallExampleFile)
    download("https://dox.ulg.ac.be/index.php/s/Vr1YNtmz16mRVkP/download",SmallExampleFile)
end

ref = FileIO.load(SmallExampleFile)
X = ref["X"]

Random.seed!(1234);
XA,offset,U,S,V,cvEOF,cvarray,errmap,musquare = @time DINEOF.DINEOFrun(X; errormap = false)

difference = (XA - ref["XA"]);

@test sqrt(mean(difference[isfinite.(difference)].^2)) ≈ 0. atol=1e-7
@test cvarray ≈ ref["cvarray"]
