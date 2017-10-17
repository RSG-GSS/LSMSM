##########################################################################
##### 10/16/2017: Modified by Miura, Hirotaka <Hirotaka.Miura@ny.frb.org>                                
##### 00/00/2017: Previously modified by                                         
##### 00/00/2017 Created by Kosar, Gizem <Gizem.Kosar@ny.frb.org>
##########################################################################   
##### Description: 
##### 	- Master program
##### Modifications:
##########################################################################

#import packages, functions, etc. - only uncomment for the first run
#include("setup.jl")
##### Run setup on all workers.
@everywhere include("setup.jl");

#initialize the fixed parameters
fp = fparams();
#initialize the data type that will contain the moments and related information
moments = initmom(fp);
#initial parameter guess
initp0 = [49.38 ,    
    -27.699581,
     1.450677,
     5.66352,  
     0.67565,  
    -0.658497, 
  	15.8,      
     1.78912,
     0.529738, 
    -1.62237,
    46.997];

#calculating the data moments and the bootstrapped weight matrix for the moments
CalcDataBootMom!(fp,moments);
#println("starting to solve and simulate")
#calculating the value of the objective function (the weighted sum of squares of the difference in data and simulated moments)
@time obj = objectivefunc!(initp0,fp,moments);
println("***** Serial result: ",obj,"\n")
@time objp = objectivefuncp(initp0,fp,moments,"parfor");
println("**** Parallel result: ",objp,"\n")
@time objp = objectivefuncp(initp0,fp,moments,"pmapbatch");
println("**** Parallel result: ",objp,"\n")
println("obj==objp?"); obj==objp
#this line is for optimization
#res = optimize(p -> objectivefunc!(p,fp,moments),initp0,NelderMead(),Optim.Options(show_trace = true))
#save("./optresults.jld", "res", res, "p0", initp0, "moments", moments)