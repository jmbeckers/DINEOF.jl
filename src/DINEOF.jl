module DINEOF

using LinearAlgebra
using Arpack
using Statistics

include("DINEOF_svd.jl")
export DINEOF_svd

include("DINEOF_svds!.jl")
export DINEOF_svds!

include("DINEOF_errormap.jl")
export DINEOF_errormap

include("DINEOF_musquare.jl")
export DINEOF_musquare

include("DINEOF_cvmask.jl")
export DINEOF_cvmask

include("DINEOFrun.jl")
export DINEOFrun

end
