#import packages, functions, etc. - only uncomment for the first run
include("setup.jl")

#initialize the fixed parameters
fp = fparams()
#initialize the data type that will contain the moments and related information
moments = initmom(fp)
#initial parameter guess
initp0 = [0.1,0.1,0.1,0.1,0.1,0.1,1.,0.1,0.1,0.1,0.1,0.1,1.]

#calculating the data moments and the bootstrapped weight matrix for the moments
CalcDataBootMom!(fp,moments)
println("starting to solve and simulate")
#calculating the value of the objective function (the weighted sum of squares of the difference in data and simulated moments)
@time obj = objectivefunc!(initp0,fp,moments)
#this line is for optimization
#res = optimize(p -> objectivefunc!(p,fp,moments),initp0,NelderMead(),Optim.Options(show_trace = true,iterations=15))
#@time obj = objectivefunc!(initp0,fp,moments)