##########################################################################   
##### 04/06/2018: Modified by Hirotaka Miura.
##### 04/05/2018: Previously modified by Hirotaka Miura.
##### 04/04/2018: Created by Gizem Kosar.
##### Description: 
##### 	- Functions to solve simulation.
##### Modifications:
#####	04/05/2018: 
#####		- Begin making adjustments for parallelization. 
#####		- Not clear why, but certain variables need @everywhere
#####			macro for code to run in parallel.
#####	04/06/2018:
#####		- Updated comments.
##########################################################################

#import packages, functions, etc. - only uncomment for the first run
#parallel = "Y-bank"
#parallel = "Y-off"
@everywhere parallel = "N"
@everywhere version = 1
@everywhere gender = "f"

@everywhere include("setup.jl")

#=
if parallel == "Y-bank"
	using ClusterManagers
	myprocs = addprocs_sge(40, queue = "background.q")
	@everywhere include("setup.jl")
elseif parallel == "Y-off"
	using ClusterManagers
	#myprocs = addprocs_sge(4)
	@everywhere include("setup.jl")
else
	include("setup.jl")
end
=#

#initialize the fixed parameters
@everywhere fp = fparams(version,gender)
#initialize the data type that will contain the moments and related information
moments = initmom(fp,version,gender)
#initial parameter guess
initp0 = 0.1*ones(51)

#calculating the data moments and the bootstrapped weight matrix for the moments
#CalcDataBootMom!(fp,moments)
println("starting to solve and simulate")

#=Random number draws for wage,hours disuility,take-up disutility, measurement error in wages, 
measurement error in assets, prob of moving, childcare cost, child transition, type prob=#
sd = 2230
srand(sd)
d = Uniform()
#this is for the wage shock and the measurement error in wage
draws = norminvcdf.(rand(d,fp.nind*fp.nper,fp.nsim,2)) 
#now adding the random variables to simulate child transition
draws = cat(3,draws, rand(d,fp.nind*fp.nper,fp.nsim))

#loading the tax and transfer system information
eitcDict = load("./dicts/eitcDict.jld")["data"]
#taxDict = load("./dicts/taxDict.jld")["data"]
ctcDict = load("./dicts/ctcDict.jld")["data"]
stDict = load("./dicts/stDict.jld")["data"]
ftDict = load("./dicts/ftDict.jld")["data"]
#ftaxp = ftaxparam()
#ctcp = ctcparam()
fgr = fgrids(fp)
lEV, lEDU = initEVEDU(fp)
@everywhere lEVr, lEDUr = initEVrEDUr(fp)
sim = initsim(fp,gender)
ccp = initccp(fp,gender)


#calculating the value of the objective function (the weighted sum of squares of the difference in data and simulated moments)
obj = objectivefunc!(initp0,fp,moments, draws, lEV, lEDU, lEVr, lEDUr,eitcDict,fgr,sim,ftDict,ctcDict,stDict,ccp,version, gender)
#@time obj = objectivefunc!(initp0,fp,moments, draws, lEV, lEDU, lEVr, lEDUr, eitcDict,fgr,sim,ftDict,ctcDict,stDict,ccp,version,gender)
#@profile obj = objectivefunc!(initp0,fp,moments, draws, lEV, lEDU, eitcDict,fgr,sim,ftDict,ctcDict,stDict)

#this line is for optimization
#res = optimize(p -> objectivefunc!(p,initp0,fp,moments, draws, lEV, lEDU, eitcDict,fgr,sim,ftDict,ctcDict,stDict),initp0,NelderMead(),Optim.Options(show_trace = true,store_trace = true))

save("./sample_results.jld", "lEV", lEV, "lEDU", lEDU, "sim", sim, "moments", moments, "obj", obj)

if parallel == "Y-bank"
	rmprocs(myprocs)
end