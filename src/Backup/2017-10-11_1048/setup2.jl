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
include("solvesim2.jl")
include("moments.jl")
include("wrapper2.jl")