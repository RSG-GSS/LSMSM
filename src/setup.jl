#import packages
using Distributions
using StatsFuns
using Parameters
import NaNMath
using Optim
using JLD


include("params.jl")
include("grids.jl")
include("funcs.jl")
include("solvesim.jl")
##### Include non-inplace simulation function.
include("solvesim2.jl")
include("moments.jl")
include("wrapper.jl")
##### Include parallel implementation.
include("wrapper2.jl")
