##########################################################################   
##### 10/13/2017: Modified by Miura, Hirotaka <Hirotaka.Miura@ny.frb.org>
##### 00/00/0000: Previously modified by 
##### 10/13/2017: Created by Kosar, Gizem <Gizem.Kosar@ny.frb.org>
##########################################################################
##### Description: 
##### 	- Module for defining new types and related functions.
##### Modifications:
##########################################################################
##### Define module name.
module Mparams
##### Import packages.
using Parameters		##### For @unpack.
##### Specify names to be exported.
export
params,
fparams,
single_1,
Oh_t,
Ch_t,
OhCh_t,
sim_t,
structureM,
initEVEDU,
initsim,
initmom,
pv2p0

#composite type to hold the parameters to be estimated
type params
    αc :: Float64
    αh0 :: Float64
    αh1 :: Float64
    αh2 :: Float64
    αh3 :: Float64
    αh4 :: Float64
    σu :: Float64
    αw0 :: Float64
    αw1 :: Float64
    αw2 :: Float64
    σw :: Float64 
end

#function to initialize the parameter structe
function params()
    p0 = params(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.)  
    return p0
end

#composite type to hold the fixed parameters of the model
type fparams 
    R::Float64
    β::Float64   
    θ::Float64  
    nper::Int           # number of periods
    enda::Int            # age at the last period
    mina::Int            # age at the first period  -- In the Blundell code mina is different from as
    maxa::Int            # age at the last period        
    ngpl::Int             # number of labor supply options
    ngu::Int           # number of grids for random shocks (joint: ngu*ngu)
    ngpk::Int             # number of grid points in assets
    maxk::Float64        # largest value of k in grid
    mink::Float64            # minumum value of k allowed
    maxngpe::Int             # max number of grid points for experience
    ratio::Float64           # ratio for the asset grid
    nind::Int               #number of individuals
    nsim::Int               #number of repetitions
    totn::Int               #total number of individuals * how many periods you track them
    nmom::Int            #the upper bound of the number of moments per section of Mstructure
    nboot::Int              #number of bootstrap iterations
end


#initializing the fparams structure
function fparams()
    fp0 = fparams(0.015, 0.98,-0.5,0,0,0,0,0,0,0,0.,0.,0,0.,0,0,0,0,0)
    fp0.R = 0.015
    fp0.β = 0.98
    fp0.θ = -0.5  
    fp0.nper  = 41            # number of periods
    fp0.enda  = 60            # age at the last period
    fp0.mina  = 20            # age at the first period  -- In the Blundell code mina is different from as
    fp0.maxa  = 60            # age at the last period        
    fp0.ngpl    = 2             # number of labor supply options
    fp0.ngu     = 6            # number of grids for random shocks (joint: ngu*ngu)
    fp0.ngpk    = 6             # number of grid points in assets
    fp0.maxk    = 20000.        # largest value of k in grid
    fp0.mink    = -100.            # minumum value of k allowed
    fp0.maxngpe = 5             # max number of grid points for experience
    fp0.ratio   = 4.0           # ratio for the asset grid
    fp0.nind    = 1000          # number of women - CHECK THIS FROM THE PSID!!!
    fp0.nsim     = 5             # number of simulations for each woman
    fp0.totn    = fp0.nind*(fp0.maxa - fp0.mina+1)
    fp0.nmom    = 25            #the upper bound of the number of moments
    fp0.nboot   = 2             #number of bootstrap iterations
    return fp0
end

#composite type to hold the expected value arrays. 
#This structure is useful in case we have a distinction between married and single individuals in a model. 
type single_1
    m0 :: Array{Float64,3}
end

#=
#composite type to hold the grid variables.
type grids
    ngpe::Array{Int,1}
    ge::Array{Int,2}
    gk::Array{Float64,1}
    gεw::Array{Float64,1}
    gεu::Array{Float64,1}
    trp::Array{Float64,2}
end
=#

#part of the simulated data type structure - state variables.
type Oh_t
    id::Int                          #ID    
    insample::Array{Int,1}           #indicator for observation in data sample
    a::Array{Int,1}                  #age
    l_1::Array{Int,1}                #previous period's labor supply choice
    k::Array{Float64,1}              #accumulated assets
    e::Array{Int,1}                  #accumulated experience
    w::Array{Float64,1}              #hourly wage rate
    inc::Array{Float64,1}            #income
    V::Array{Float64,1}              #expected present value of present and future utility, conditional on today's c,l
end
#part of the simulated data type structure - choice variables.
type Ch_t
    l::Array{Int,1}                  #labor supply choice (1= not working, 2 = working)
    c::Array{Float64,1}              #consumption
end

#part of the simulated data type structure.
type OhCh_t
    Oh::Oh_t
    Ch::Ch_t 
end

#part of the simulated data type structure.
type sim_t
    id::Int                 #PSID individual woman identifier
    yob::Int                #PSID year of birth
    aa::OhCh_t               #An array of Oh and Ch structures
end


#part of the moments type structure
type structureM
    cols::Dict{Symbol,Int} #moment names, which are also used as indices
    dtamom::Array{Float64,1} #moments from the data
    simmom::Array{Float64,1} #moments from the simulated data
    wgtcov::Array{Float64} #moment weights for the estimation algorithm
    withvar::Array{Int,1} #moments with variation    
    objfunc::Float64 #objective function -- diff'*estweight*diff
end

#=
#function to calculate the grids for the experience, assets and random shocks
function grids(fp::fparams,p0::params)
    @unpack nper, maxngpe, ngpk, ngu = fp
    gr0 = grids(zeros(nper),zeros(maxngpe,nper),Float64.(zeros(ngpk)),Float64.(zeros(ngu)),Float64.(zeros(ngu)),Float64.(zeros(ngu,ngu)))
    gr0.ngpe, gr0.ge = popge(fp)
    gr0.gk = popgk(fp)
    gr0.gεw, gr0.gεu = popgerr(fp,p0)
    gr0.trp = trprob(gr0.gεw, gr0.gεu,p0,fp)
    return gr0
end
=#

#initializing the expected value function array and the derivative of expected utility array (necessary for the Euler Equation)
function initEVEDU(fp::fparams)
    @unpack nper, ngpl, ngpk, maxngpe = fp
    lEV = Array{single_1}(nper)
    lEDU = Array{single_1}(nper) 
    for i in eachindex(lEV)
        lEV[i] = single_1(0*ones(ngpl,ngpk,maxngpe))
        lEDU[i] = single_1(0*ones(ngpl,ngpk,maxngpe))
    end
    return lEV, lEDU
end

#initializing the simulated data type structure
function initsim(fp::fparams)
    @unpack nind, nsim, nper = fp
    sim = Array{sim_t}(nind,nsim)     
    for i in eachindex(sim)                
        sim[i] = (sim_t(0,0,OhCh_t(Oh_t(0,zeros(nper),zeros(nper),zeros(nper),zeros(nper),zeros(nper),zeros(nper),zeros(nper),
            zeros(nper)),Ch_t(zeros(nper),zeros(nper)))))
    end
    data = readdlm("data.txt",'\t')
    data[data.==-9] = NaN
    data_cols = [:id, :age, :actual, :ls, :ls_1, :wage, :assets, :exper]                 
    cols = Dict{Symbol,Int}()
    for (i,k) in enumerate(data_cols); cols[k] = i end    
    for i = 1:nind
        for rep = 1:nsim
            sim[i,rep].id = i
            sim[i,rep].yob = 1986
            sim[i,rep].aa.Oh.l_1[1] = 1            
            sim[i,rep].aa.Oh.k[1] = data[(i-1)*41+1,cols[:assets]]            
            for j = 1:nper
                sim[i,rep].aa.Oh.id = i
                sim[i,rep].aa.Oh.a[j] = 20+j-1
                sim[i,rep].aa.Oh.insample[j] = 1
            end
        end
    end
    return sim
end

#initiazligin the moments type data structure
function initmom(fp::fparams)
    @unpack nmom = fp
    cols = Dict{Symbol,Int}()
    columns = [:prE, :prE25, :prE30, :prE3040, :prE4059, :prE50p,
                :UU, :EE,
                :mw, :mw25, :mw30, :mw3040, :mw4050, :mw50p, :stdw,
                :mwq25e, :mwq50e, :mwq75e, :mwq90e,
                :mk, :mk30, :mk3040, :mk4050, :mk50p, :stdk]
    for (i,k) in enumerate(columns); cols[k] = i end
    dtamom = Array(Float64,nmom)
    simmom = Array(Float64,nmom)
    wgtcov = Array{Float64}()
    withvar = Array(Int,nmom)
    objfunc = 0.

    m = structureM(cols,dtamom,simmom,wgtcov,withvar,objfunc)
    return m
end

function pv2p0(p0::Array{Float64,1})
    p = params()
    p.αc  = -exp((p0[1]-10)/100)
    p.αh0 = p0[2]
    p.αh1 = p0[3]
    p.αh2 = p0[4]
    p.αh3 = p0[5]
    p.αh4 = p0[6]
    p.σu = exp((p0[7]-10)/100)
    p.αw0 = p0[8]
    p.αw1 = p0[9]
    p.αw2 = p0[10]
    p.σw = exp((p0[11]-100)/100)
    return p
end

##### End module definition.
end