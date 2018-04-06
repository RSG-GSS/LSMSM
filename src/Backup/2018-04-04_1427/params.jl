#composite type to hold the parameters to be estimated
type params
	#utility - λc
	αc1 :: Float64
	αc2 :: Float64
	αc3 :: Float64
	αc4 :: Float64
	#utility - λh
	αh0 :: Float64
	αh12 :: Float64
    αh13 :: Float64
    αh14 :: Float64    
	αh21 :: Float64
    αh22 :: Float64
	#utility - λp
	αp11 :: Float64
	αp12 :: Float64
	αp13 :: Float64
	αp14 :: Float64
	#lower bound for assets
	αb0 :: Float64
	αb1 :: Float64
	αb2 :: Float64
	αb32 :: Float64
	αb33 :: Float64
	αb34 :: Float64
	#wage
	αw0 :: Float64
	αw11 :: Float64
	αw12 :: Float64
	αw13 :: Float64
	#αw2 :: Float64
	αw32 :: Float64
	αw33 :: Float64
	αw34 :: Float64
    #αw41 :: Float64
	#αw42 :: Float64
	#αw43 :: Float64
	#αw44 :: Float64    
	αw51 :: Float64
	αw52 :: Float64
    αw53 :: Float64
    αw54 :: Float64
    αw61 :: Float64
    αw62 :: Float64
	#random shocks
	σw :: Float64
	σmw :: Float64
end


#function to initialize the parameter structe
function params()
    p0 = params(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.)        
    return p0
end

#composite type to hold the fixed parameters of the model
type fparams 
    R::Float64              #1+r
    β::Float64              #the discount factor   
    θ::Float64              #the CRRA parameter 
    nts::Int                # levels of schooling 
    nper::Int               # max number of periods (using the least educated group: 60-16+1)    
    ngpl::Int               # number of labor supply options
    ngpl1::Int
    ngpt::Int               # number of participation options
    ngpch::Int 				# number of kids options
    ngu::Int                # number of grids for random shocks     
    ngpk::Int               # number of grid points in assets
    maxngpe::Int            # max number of grid points for experience
    ngpeitc::Int 			# number of grid points for eitc
    npolicy::Int 			# number of different tax policies
    nstate::Int 			# number of state options
    #nrace::Int 				# number of race options    
	ntype::Int 				# number of unobserved types
	enda::Int               # age at the last period of working life        
    maxa::Int               # age at the last period              
    as::Array{Int,1}        # age when entering the labor market (depends on education level)
    sysY0::Array{Int,1}
    sysY1::Array{Int,1}
    mina::Int               # minimum age (first period age for the lowest educated)
    maxk::Float64           # largest value of k in grid
    mink::Float64           # minumum value of k allowed    
    ratiok::Float64         # ratio for the asset grid
    ratioe::Float64         # ratio for the eitc grid
    nind::Int               # number of individuals
    nsim::Int               # number of repetitions
    totn::Int               # total number of individuals * how many periods you track them
    nmom::Int               # the upper bound of the number of moments per section of Mstructure    
    nboot::Int              # number of bootstrap iterations
    toti::Int 				# total number of indices for the outer structure of the lEV, lEDU
    totcti::Int
    totsi::Int
    maxeitc::Float64 		# largest value of eitc allowed
    eqs::Array{Float64,1} 	#equivalence scale based on the number of kids
    wgts::Array{Float64,1}  #weights for the Gaussian quadratures
    nds::Array{Float64,1}   #nodes for the Gaussian quadrature
    hrs::Array{Float64,1}
    minc::Array{Float64,2}
end


#initializing the fparams structure
function fparams(version::Int, gender::String)
    fp0 = fparams(0.,0.,0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,[0,0,0,0],zeros(Int,4),zeros(Int,4),0,0.,0.,0.,0.,0,0,0,0,0,0,0,0,0.,
        [0.,0.,0.],zeros(Float64,5),zeros(Float64,5),zeros(Float64,5),zeros(Float64,4,3))
    fp0.R = 1.015             # 1+r
    fp0.β = 0.98              # the discount factor
    fp0.θ = -0.5              # the CRRA parameter
    fp0.nts = 4               # levels of schooling
    fp0.nper  = 45			  # max number of periods (using the least educated group: 60-16+1)
    fp0.ngpl    = 5           # number of labor supply options
    fp0.ngpl1 = 3
    fp0.ngpch 	= 3 		  # number of kids options
    fp0.ngpt    = 2           # number of participation options 
    fp0.ngu     = 5     # number of grids for random shocks (joint: ngu*ngu)    
    fp0.ngpk    = 4           # number of grid points in assets
    fp0.maxngpe = 4           # max number of grid points for experience
    fp0.ngpeitc = 4 		  # number of grid points for eitc
    fp0.npolicy = 4 		  # number of different tax policies
    fp0.nstate = 15 		  # number of state options
    #fp0.nrace = 2 			  # number of race options
    fp0.ntype = 2 			  # number of unobserved types
    fp0.enda  = 70            # age at the last period
    fp0.maxa  = 60            # age at the last period of working life
    fp0.as = [16,18,20,22]    # age when entering the labor market (depends on education level)
    fp0.sysY0 = [1985,1987,1994,2000]
    fp0.sysY1 = [1986,1993,1999,2007]
    fp0.mina  = 16            # minimum age (first period age for the lowest educated)     
    fp0.maxk    = 800.      # largest value of k in grid
    fp0.mink    = -100.       # minumum value of k allowed    
    fp0.ratiok   = 4.0         # ratio for the asset grid
    fp0.ratioe   = 1.2         # ratio for the eitc grid
    fp0.nind    = 1452        # number of women - CHECK THIS FROM THE PSID!!!
    fp0.nsim     = 5          # number of simulations for each woman
    fp0.totn    = fp0.nind*(fp0.maxa - fp0.mina+1)
    if gender == "f"
        if version == 1        
            fp0.nmom    = 153
        elseif version == 2
            fp0.nmom = 169
        else
            fp0.nmom = 178
        end
    else
        if version == 1        
            fp0.nmom    = 152
        elseif version == 2
            fp0.nmom = 168
        else
            fp0.nmom = 177
        end
    end
    fp0.nboot   = 5           #number of bootstrap iterations
    fp0.toti 	= fp0.ngpch*fp0.nper
    fp0.totcti  = fp0.nts*fp0.nstate #*(fp0.ngpl1-1)*fp0.ngpeitc + fp0.nts*fp0.nstate#*fp0.ntype
    fp0.totsi   = fp0.nsim*fp0.nind
    fp0.maxeitc = 5.
    fp0.eqs = [1.,1.3,1.6] 		#equivalence scale based on the number of kids    
    π1_2 = 1. /(π^(1/2))  
    fp0.wgts = π1_2.*[0.0199532,0.393619,0.945309,0.393619,0.0199532] 
    fp0.nds = [-2.02018,-0.958572,2.40258e-16,0.958572,2.02018]
    fp0.hrs = [0.,755.,1520.,2020.,2720.]
    fp0.minc =  0.12*[243.76015   585.52734375    1030.68603515625;
                278.58475   601.12744140625 1058.146484375;
                122.40231   657.0127437770584   870.8584483816646;
                149.3828    739.7640067152337   990.3343786042017]
    return fp0
end


#composite type to hold the expected value arrays. 
#This structure is useful in case we have a distinction between married and single individuals in a model. 
type single_1
    m0 :: Array{Float64,5}
end

type index
	ix::Array{Int,2}
	ch::Array{Int,1}
	a::Array{Int,1}
end

type ctindex
    ix::Array{Int,2}
    #ty::Array{Int,1}
    st::Array{Int,1}
    s::Array{Int,1}
    #l1::Array{Int,1}
    #cr::Array{Int,1}
    #r::Array{Int,1}    
end

type simindex
    ix::Array{Int,2}
    rep::Array{Int,1}
    ind::Array{Int,1}
end


#composite type to hold the grid variables.
#fgrids don't depend on the parameter vector
type fgrids    
    ngpe::Array{Int,1}
    ge::Array{Int64,2}    
    geitc::Array{Float64,1}
    gk::Array{Float64}    
    iEmtx::index
    icEmtx::ctindex
    isEmtx::simindex
    trkids::Array{Float64,4}
    #adjw::Array{Float64,1}
    fe1::Array{Int64,3}
    fe1lb::Array{Int64,3}
    fe1ub::Array{Int64,3}
    #pe1::Array{Int64,3}
    #pe1lb::Array{Int64,3}
    #pe1ub::Array{Int64,3}
end


type dparams
    utilpar::Array{Float64,4}
    dutilpar::Array{Float64,4}
    eepar::Array{Float64,4}	
	boundk::Array{Float64,2}	
    nodes::Array{Float64,1}    
end

function dparams(fp::fparams, p0::params, finaln::Array{Float64,1})
	@unpack β, R, θ, eqs, ngpl, ngpch, nts, ngu, ngpt, nper, mina = fp
	dp0 = dparams(zeros(Float64,nts,ngpch,ngpl,ngpt), zeros(Float64,nts,ngpch,ngpl,ngpt),zeros(Float64,nts,ngpch,ngpl,ngpt), 
		zeros(Float64,nper+10,nts),zeros(Float64,ngu))	        
    peqs = eqs.^(-θ)
    peqs2 = eqs.^(θ/(θ-1))
    pRβ = (R*β)^(1/(θ-1.))
    λc = Array{Float64}(nts,ngpch,ngpl,ngpt)
    for ip = 1:ngpt, i=1:ngpl, ch = 1:ngpch, si = 1:nts
        λh = lambdah(p0,i)
        xh = fxh(p0,si,ch)
        xp = fxp(p0,si)
        λc[si,ch,i,ip] =exp(λh*xh+xp*convert(Float64,ip-1))
        dp0.utilpar[si,ch,i,ip] = λc[si,ch,i,ip]*peqs[ch]/θ
    end
    pλc = λc.^(1/(1-θ))
    for ch = 1:ngpch
        dp0.eepar[:,ch,:,:] = pRβ*peqs2[ch]*pλc[:,ch,:,:]        
    end
    dp0.dutilpar = dp0.utilpar*θ
    for si = 1:fp.nts, ai = 1:fp.nper+10
    	dp0.boundk[ai,si] = boundk(p0,ai+mina-1,si)
    end
    dp0.nodes = finaln
    return dp0
end


#part of the simulated data type structure - state variables.
type Oh_t
    policy::Array{Int,1} 		     #which tax system
    insample::Array{Int,1}           #indicator for observation in data sample
    insampleass::Array{Int,1}        #indicator for observation in data sample
    year::Array{Int,1} 				 #year
    a::Array{Int,1}                  #age
    l_1::Array{Int,1}                #previous period's labor supply choice
    k::Array{Float64,1}              #assets
    fte::Array{Int64,1}              #accumulated ft experience    
    eitc::Array{Float64,1} 			 #eitc payments to be received that year
    w::Array{Float64,1}              #hourly wage rate
    ch::Array{Int,1} 			     #number of children
    inc::Array{Float64,1}            #income
    V::Array{Float64,1}              #expected present value of present and future utility, conditional on today's c,l
end
#part of the simulated data type structure - choice variables.
type Ch_t
    l::Array{Int,1}                  #=labor supply choice (1= not working, 2 = working for 750 hours, 
    									3 = working for 1515 hours, 4 = working for 2000 hours, 
    									5 = working for 2775 hours) =#
    hrs::Array{Int,1}                                        
    c::Array{Float64,1}              #consumption
    p::Array{Int,1} 				 #EITC take-up (1=no take-up, 2=participates)
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
    #th::Int 				#unobserved type
    race::Int               #race
    sex::Int                #1:female, 2:male
    s::Int                  #education
    st::Int                 #state of residence
    sysA0::Array{Int,1} 	#age at which first experiences each tax system
    sysA1::Array{Int,1} 	#age at which last experiences each tax system
    aa::OhCh_t               #An array of Oh and Ch structures    
end


#part of the moments type structure
type structureM
    cols::Dict{Symbol,Int} #moment names, which are also used as indices
    dtamom::Array{Float64,1} #moments from the data
    simmom::Array{Float64,1} #moments from the simulated data
    wgtcov::Array{Float64,1} #moment weights for the estimation algorithm
    withvar::Array{Int,1} #moments with variation    
    objfunc::Float64 #objective function -- diff'*estweight*diff
end


type shocks
    εw::Array{Float64,2}    #wage    
    εmw::Array{Float64,2}   #measurement error for wage
    #εma::Array{Float64,2}   #measurement error for assets
    #εth::Array{Float64,2}   #prob of type
    #εpm::Array{Float64,3}   #prob of moving
    εtrch::Array{Float64,2} #transition for kids
end

function initindex(fp::fparams)
	@unpack toti, maxa, mina, nts, ngpch, nper = fp
	iEmtx =index(zeros(Int,ngpch,nper),zeros(Int,toti),zeros(Int,toti))
	return iEmtx
end

function initctindex(fp::fparams)
    @unpack totcti, nstate, nts = fp
    icEmtx =ctindex(zeros(Int,nstate,nts),zeros(Int,totcti),zeros(Int,totcti))
    return icEmtx
end

function initsimindex(fp::fparams)
    @unpack nsim, nind, totsi = fp
    isEmtx = simindex(zeros(Int,nsim,nind),zeros(Int,totsi),zeros(Int,totsi))
    return isEmtx
end


function fgrids(fp::fparams)
    @unpack nper, ntype, nstate, totcti, maxngpe, ngpk, ngpeitc, nts, ngpch, toti,wgts, ngpl, ngpl1, nsim, nind, totsi = fp
	fgr0 =(fgrids(zeros(Int,nper),zeros(Int64,maxngpe,nper),zeros(Float64,ngpeitc),zeros(Float64,ngpk),index(zeros(Int,ngpch,nper),zeros(Int,toti),zeros(Int,toti)),
    ctindex(zeros(Int,nstate,nts),zeros(Int,totcti),zeros(Int,totcti)),
    simindex(zeros(Int,nsim,nind),zeros(Int,totsi),zeros(Int,totsi)),
    zeros(Float64,nper,nts,ngpch,ngpch),zeros(Int64,ngpl,maxngpe,nper),zeros(Int64,ngpl,maxngpe,nper),zeros(Int64,ngpl,maxngpe,nper)))
	fgr0.ngpe, fgr0.ge = popge(fp)
	fgr0.gk = popgk(fp)
	fgr0.geitc = popgeitc(fp)
	fgr0.iEmtx = pop_Emtx(fp)
    fgr0.icEmtx = pop_cEmtx(fp)
    fgr0.isEmtx = pop_sEmtx(fp)
    fgr0.trkids = inittrkids(fp)
    fgr0.fe1, fgr0.fe1lb, fgr0.fe1ub = nxtebds(ngpl, maxngpe, fgr0.ngpe, fgr0.ge, nper)   
    #lw = length(wgts)
    #i = 1
    #temp = zeros(Float64,lw,lw,lw,lw)
    #for i1 = 1:lw, i2 = 1:lw, i3 = 1:lw,i4 = 1:lw         
    #    temp[i4,i3,i2,i1] = wgts[i4]*wgts[i3]*wgts[i2]*wgts[i1]
    #end
    #fgr0.adjw = temp[:]
    return fgr0
end

#initializing the expected value function array and the derivative of expected utility array (necessary for the Euler Equation)
function initEVEDU(fp::fparams)
    @unpack toti, ngpl1, ngpk, maxngpe, ngpeitc, nstate, nts= fp
    lEV = Array{single_1}(nstate,nts)
    lEDU = Array{single_1}(nstate,nts) 
    lEV .= single_1(zeros(Float64,ngpk,maxngpe,ngpeitc,ngpl1,toti))
    lEDU .= single_1(zeros(Float64,ngpk,maxngpe,ngpeitc,ngpl1,toti))
    return lEV, lEDU
end

function initEVrEDUr(fp::fparams)
    @unpack ngpk, nts, ngpch = fp
    lEVr = zeros(Float64,ngpk,ngpch,nts)
    lEDUr = zeros(Float64,ngpk,ngpch,nts)
    return lEVr, lEDUr
end

#=function initEVEDU(fp::fparams)
    @unpack toti, totcti, nrace, ntype, ngpl1, ngpk, maxngpe, ngpeitc, nstate = fp
    lEV = Array{Float64}(ngpeitc,ngpl1,ngpk,maxngpe,toti,totcti)
    lEDU = Array{Float64}(ngpeitc,ngpl1,ngpk,maxngpe,toti,totcti)
    return lEV, lEDU
end=#

#initializing the simulated data type structure
function initsim(fp::fparams, gender::String)
    @unpack nind, nsim, nper,as, sysY0, sysY1, npolicy = fp
    sim = Array{sim_t}(nind,nsim)     
    for i in eachindex(sim)                
        sim[i] = (sim_t(0,0,0,0,0,0,zeros(Int,4),zeros(Int,4),OhCh_t(Oh_t(zeros(Int,nper),zeros(Int,nper),zeros(Int,nper),zeros(Int,nper),zeros(Int,nper),
            zeros(Int,nper),zeros(Float64,nper),zeros(Int64,nper),zeros(Float64,nper),
            zeros(Float64,nper),zeros(Int,nper),zeros(Float64,nper),zeros(Float64,nper)),
            Ch_t(zeros(Int,nper),zeros(Int,nper),zeros(Float64,nper),zeros(Int,nper)))))
    end
    data = readdlm("./data/initial_conditions_dataset_$(gender).txt",'\t') #initial conditions    
    data_cols = [:id, :rep, :yob, :race, :education, :state, :assets]        
    data_cols2 = [:id, :age, :year, :insample, :insampleass, :policy, :sysa0, :sysa1]        
    obsW = readdlm("./data/ObsWindow_$(gender).txt",'\t') #obsWindow   
    cols = Dict{Symbol,Int}()
    cols2 = Dict{Symbol,Int}()
    for (i,k) in enumerate(data_cols); cols[k] = i end
    for (i,k) in enumerate(data_cols2); cols2[k] = i end    
    k = 1
    for i = 1:nind
        for rep = 1:nsim
            sim[i,rep].id = data[k,cols[:id]]
            sim[i,rep].yob = data[k,cols[:yob]]
            sim[i,rep].race = data[k,cols[:race]]
            sim[i,rep].s = data[k,cols[:education]]
            sim[i,rep].st = data[k,cols[:state]]
            sim[i,rep].sex = 1 #female
            sim[i,rep].sysA0[1] = min(as[sim[i,rep].s],sysY0[1]-sim[i,rep].yob) #if older than starting age in 1985, starting age for first policy is 16,18,20 or 22.
            sim[i,rep].sysA1[npolicy] = max(60,sysY1[npolicy]-sim[i,rep].yob) #if younger than 60 in 07, ending age for last policy is 60
            for sys = 1:npolicy-1
                sim[i,rep].sysA1[sys] = sysY1[sys] - sim[i,rep].yob
            end
            for sys = 2:npolicy
                sim[i,rep].sysA0[sys] = sim[i,rep].sysA1[sys-1]+1
            end                       
            sim[i,rep].aa.Oh.l_1[1:1+2*(sim[i,rep].s-1)] .= 1            
            sim[i,rep].aa.Oh.k[1:1+2*(sim[i,rep].s-1)] .= data[k,cols[:assets]]
            sim[i,rep].aa.Oh.fte[1:1+2*(sim[i,rep].s-1)] .= 0.            
            sim[i,rep].aa.Oh.eitc[1:1+2*(sim[i,rep].s-1)] .= 0.
            sim[i,rep].aa.Oh.ch[1:1+2*(sim[i,rep].s-1)] .= 1
            for j = 1:nper
                #println("$i $rep $j")
                sim[i,rep].aa.Oh.policy[j] = obsW[(i-1)*nper+j,cols2[:policy]]
                sim[i,rep].aa.Oh.insample[j] = obsW[(i-1)*nper+j,cols2[:insample]]
                sim[i,rep].aa.Oh.insampleass[j] = obsW[(i-1)*nper+j,cols2[:insampleass]]
                sim[i,rep].aa.Oh.year[j] = convert(Int,obsW[(i-1)*nper+j,cols2[:year]])
                sim[i,rep].aa.Oh.a[j] = obsW[(i-1)*nper+j,cols2[:age]]                
            end
        end
    k = k+1
    end
    return sim
end

@views function simClean!(sim::Array{sim_t,2},nind::Int, nsim::Int, nper::Int)
    for rep = 1:nsim
        for i = 1:nind
            #=if εth[i,rep]<= πth1
                sim[i,rep].th = 1
            else
                sim[i,rep].th = 2
            end    =#
            sim[i,rep].aa.Oh.w[1:1+2*(sim[i,rep].s-1)] = 0.
            sim[i,rep].aa.Oh.inc[1:1+2*(sim[i,rep].s-1)] = 0.
            sim[i,rep].aa.Oh.V[1:1+2*(sim[i,rep].s-1)] = 0.
            for j = (2*sim[i,rep].s):nper
                sim[i,rep].aa.Oh.l_1[j] = 0
                sim[i,rep].aa.Oh.k[j] = 0.
                sim[i,rep].aa.Oh.fte[j] = 0.
                sim[i,rep].aa.Oh.eitc[j] = 0.
                sim[i,rep].aa.Oh.w[j] = 0.
                sim[i,rep].aa.Oh.ch[j] = 0
                sim[i,rep].aa.Oh.inc[j] = 0.
                sim[i,rep].aa.Oh.V[j] = 0.                
            end
        end
    end
    return sim    
end

#initializing the moments type data structure
function initmom(fp::fparams,version::Int, gender::String)
    @unpack nmom = fp
    cols = Dict{Symbol,Int}()
    set1 = [30,3040,4050,50]
    set2 = ["hsd", "hsg", "colld", "collg"]
    set3 = ["0k", "1k", "2k"]
    #set4 = ["DC", "IL", "IN", "KS", "MA", "MD", "ME", "MN", "NJ", "NY", "OK", "OR", "RI", "VT", "WI", "oth" ]
    set5 = ["Q25", "Q50", "Q75", "Q90"]
    columns = [[Symbol("prd1_$(j)") for j in set2]; #Moments 1-4
                [Symbol("prd1")]; #moment 5 
    			[Symbol("prd1_$(j)") for j in set3]; #Moments 6-8
                [Symbol("prd2_$(j)") for j in set2]; #Moments 9-12
                [Symbol("prd2")]; #moment 13
                [Symbol("prd2_$(j)") for j in set3]; #Moments 14-16
                [Symbol("prd3_$(j)") for j in set2]; #Moments 17-20
                [Symbol("prd3")]; #moment 21
                [Symbol("prd3_$(j)") for j in set3]; #Moments 22-24
                [Symbol("prd4_$(j)") for j in set2]; #Moments 25-28
                [Symbol("prd4")]; #moment 29 
                [Symbol("prd4_$(j)") for j in set3]; #Moments 30-32
                [:prNE, :prEN,:ρ_hrs_lhrs]; #Moments 33-35
                [:meanhrs_emp_IL_pol1]; [:prE_IL_pol1]; [:meanhrs_emp_IN_pol1]; [:prE_IN_pol1]; [:meanhrs_emp_IN_pol1]; [:prE_IN_pol1];                
                [:meanhrs_emp_MA_pol1]; [:prE_MA_pol1]; [:meanhrs_emp_MD_pol1]; [:prE_MD_pol1]; [:meanhrs_emp_MN_pol1]; [:prE_MN_pol1];
                [:meanhrs_emp_NJ_pol1]; [:prE_NJ_pol1]; [:meanhrs_emp_NY_pol1]; [:prE_NY_pol1]; [:meanhrs_emp_OR_pol1]; [:prE_OR_pol1];
                [:meanhrs_emp_WI_pol1]; [:prE_WI_pol1]; [:meanhrs_emp_IL_pol1]; [:prE_IL_pol1]; [:meanhrs_emp_IN_pol1]; [:prE_IN_pol1];
                [:meanhrs_emp_MA_pol1]; [:prE_MA_pol1]; [:meanhrs_emp_MD_pol1]; [:prE_MD_pol1]; [:meanhrs_emp_MN_pol1]]
    if gender == "f"
        columns = [columns; [:prE_MN_pol1]]
    end       
    columns = [columns;      
                [:meanhrs_emp_NJ_pol1]; [:prE_NJ_pol1]; [:meanhrs_emp_NY_pol1]; [:prE_NY_pol1]; [:meanhrs_emp_OR_pol1]; [:prE_OR_pol1];
                [:meanhrs_emp_WI_pol1]; [:prE_WI_pol1];
                [:meanhrs_emp_IL_pol2]; [:prE_IL_pol2]; [:meanhrs_emp_IN_pol2]; [:prE_IN_pol2]; [:meanhrs_emp_IN_pol2]; [:prE_IN_pol2];                
                [:meanhrs_emp_MA_pol2]; [:prE_MA_pol2]; [:meanhrs_emp_MD_pol2]; [:prE_MD_pol2]; [:meanhrs_emp_MN_pol2]; [:prE_MN_pol2];
                [:meanhrs_emp_NJ_pol2]; [:prE_NJ_pol2]; [:meanhrs_emp_NY_pol2]; [:prE_NY_pol2]; [:meanhrs_emp_OR_pol2]; [:prE_OR_pol2];
                [:meanhrs_emp_WI_pol2]; [:prE_WI_pol2]; [:meanhrs_emp_IL_pol2]; [:prE_IL_pol2]; [:meanhrs_emp_IN_pol2]; [:prE_IN_pol2];
                [:meanhrs_emp_MA_pol2]; [:prE_MA_pol2]; [:meanhrs_emp_MD_pol2]; [:prE_MD_pol2]; [:meanhrs_emp_MN_pol2]; [:prE_MN_pol2];
                [:meanhrs_emp_NJ_pol2]; [:prE_NJ_pol2]; [:meanhrs_emp_NY_pol2]; [:prE_NY_pol2]; [:meanhrs_emp_OR_pol2]; [:prE_OR_pol2];
                [:meanhrs_emp_WI_pol2]; [:prE_WI_pol2];
                [:meanhrs_emp_IL_pol3]; [:prE_IL_pol3]; [:meanhrs_emp_IN_pol3]; [:prE_IN_pol3]; [:meanhrs_emp_IN_pol3]; [:prE_IN_pol3];                
                [:meanhrs_emp_MA_pol3]; [:prE_MA_pol3]; [:meanhrs_emp_MD_pol3]; [:prE_MD_pol3]; [:meanhrs_emp_MN_pol3]; [:prE_MN_pol3];
                [:meanhrs_emp_NJ_pol3]; [:prE_NJ_pol3]; [:meanhrs_emp_NY_pol3]; [:prE_NY_pol3]; [:meanhrs_emp_OR_pol3]; [:prE_OR_pol3];
                [:meanhrs_emp_WI_pol3]; [:prE_WI_pol3]; [:meanhrs_emp_IL_pol3]; [:prE_IL_pol3]; [:meanhrs_emp_IN_pol3]; [:prE_IN_pol3];
                [:meanhrs_emp_MA_pol3]; [:prE_MA_pol3]; [:meanhrs_emp_MD_pol3]; [:prE_MD_pol3]; [:meanhrs_emp_MN_pol3]; [:prE_MN_pol3];
                [:meanhrs_emp_NJ_pol3]; [:prE_NJ_pol3]; [:meanhrs_emp_NY_pol3]; [:prE_NY_pol3]; [:meanhrs_emp_OR_pol3]; [:prE_OR_pol3];
                [:meanhrs_emp_WI_pol3]; [:prE_WI_pol3];
                [:meanhrs_emp_IL_pol4]; [:prE_IL_pol4]; [:meanhrs_emp_IN_pol4]; [:prE_IN_pol4]; [:meanhrs_emp_IN_pol4]; [:prE_IN_pol4];                
                [:meanhrs_emp_MA_pol4]; [:prE_MA_pol4]; [:meanhrs_emp_MD_pol4]; [:prE_MD_pol4]; [:meanhrs_emp_MN_pol4]; [:prE_MN_pol4];
                [:meanhrs_emp_NJ_pol4]; [:prE_NJ_pol4]; [:meanhrs_emp_NY_pol4]; [:prE_NY_pol4]; [:meanhrs_emp_OR_pol4]; [:prE_OR_pol4];
                [:meanhrs_emp_WI_pol4]; [:prE_WI_pol4]; [:meanhrs_emp_IL_pol4]; [:prE_IL_pol4]; [:meanhrs_emp_IN_pol4]; [:prE_IN_pol4];
                [:meanhrs_emp_MA_pol4]; [:prE_MA_pol4]; [:meanhrs_emp_MD_pol4]; [:prE_MD_pol4]; [:meanhrs_emp_MN_pol4]; [:prE_MN_pol4];
                [:meanhrs_emp_NJ_pol4]; [:prE_NJ_pol4]; [:meanhrs_emp_NY_pol4]; [:prE_NY_pol4]; [:meanhrs_emp_OR_pol4]; [:prE_OR_pol4];
                [:meanhrs_emp_WI_pol4]; [:prE_WI_pol4]; #moments 36-115
                [Symbol("μw_age$(j)") for j in set1]; #moments 116-119
                [:μw,:σw]; #moments 120-121
                [Symbol("μw_d$(j)") for j in 2:4]; #moments 122-124
                [:μw_d1d2, :μw_d3d4]; #moments 125-126
                [Symbol("μw_$(j)") for j in set2]; #moments 127-130
                [:ρ_w_age, :ρ_w_kids, :ρ_w_hrs, :ρ_w_fte]; #moments 131-134
                [Symbol("μw_$(i)_ft$(j)") for j in set5 for i in set2]]  #moments 135-150
    if version == 1
        columns = [columns; [Symbol("prpt_$(j)") for j in set3]] #moments 176-178
    elseif version == 2
        columns = [columns; [Symbol("pers_$(i)_d$(j)") for j in 1:4 for i in set2];  #moments 151-166
                [Symbol("prpt_$(j)") for j in set3]] #moments 176-178
    else
        columns = [columns;
                [Symbol("pers_$(i)_d$(j)") for j in 1:4 for i in set2];  #moments 151-166
                [:μnw, :σnw, :Q10nw, :Q50nw, :Q90nw] #moments 167-171
                [Symbol("μnw_d$(j)") for j in 1:4]; #moments 172-175
                [Symbol("prpt_$(j)") for j in set3]] #moments 176-178
    end
    #fp.nmom = length(columns)   
	for (i,k) in enumerate(columns); cols[k] = i end
	dtamom = readdlm("./data/moments_$(gender)_$(version).txt",'\t')[:,1]
	simmom = Array{Float64}(nmom)
	wgtcov = readdlm("./data/weights_$(gender)_$(version).txt",'\t')[:,1]
	withvar = Array{Int}(nmom)
	objfunc = 0.
	m = structureM(cols,dtamom,simmom,wgtcov,withvar,objfunc)
	return m
end


function pv2p0(p0::Array{Float64,1})
	#bound the parameters here
    p = params()
    #utility - hours
    #i'm restricting these to be <0
    i=1
    #utility - λc    
    p.αc1 = exp((p0[i]-100)/100)
    p.αc2 = exp((p0[i+1]-100)/100)
    p.αc3 = exp((p0[i+2]-100)/100)
    p.αc4 = exp((p0[i+3]-100)/100)
    #utility - λh    
    i = i+4
    p.αh0 = (p0[i]-100)/100
    p.αh12 = (p0[i+1]-100)/100
    p.αh13 = (p0[i+2]-100)/100
    p.αh14= (p0[i+3]-100)/100
    p.αh21 = (p0[i+4]-100)/100
    p.αh22 = (p0[i+5]-100)/100
    #utility - λp
    i = i+6
    p.αp11 = (p0[i]-100)/100
    p.αp12 = (p0[i+1]-100)/100
    p.αp13 = (p0[i+2]-100)/100
    p.αp14 = (p0[i+3]-100)/100
    #lower bound for assets
    i = i+4
    p.αb0 = (p0[i]-100)/100
    p.αb1 = (p0[i+1]-100)/100
    p.αb2 = (p0[i+2]-100)/100
    p.αb32 = (p0[i+3]-100)/100
    p.αb33 = (p0[i+4]-100)/100
    p.αb34 = (p0[i+5]-100)/100
    #wage
    i = i+6
    p.αw0 = (p0[i]-100)/100
    p.αw11 = (p0[i+2]-100)/100
    p.αw12 = (p0[i+3]-100)/100
    p.αw13 = (p0[i+4]-100)/100
    #p.αw2 = (p0[i+5]-100)/100
    p.αw32 = (p0[i+6-1]-100)/100
    p.αw33 = (p0[i+7-1]-100)/100
    p.αw34 = (p0[i+8-1]-100)/100
    #p.αw41 = (p0[i+9]-100)/100
    #p.αw42 = (p0[i+10]-100)/100
    #p.αw43 = (p0[i+11]-100)/100
    #p.αw44 = (p0[i+12]-100)/100    
    p.αw51 = (p0[i+9-1]-100)/100
    p.αw52 = (p0[i+10-1]-100)/100    
    p.αw53 = (p0[i+11-1]-100)/100    
    p.αw54 = (p0[i+12-1]-100)/100    
    p.αw61 = (p0[i+13-1]-100)/100    
    p.αw62 = (p0[i+14-1]-100)/100        
    #random shocks
    i = i + 14    
    p.σw = exp((p0[i+1]-100)/100)/10
    p.σmw = exp((p0[i+2]-100)/100)/10   
    return p
end

function initshocks(fp::fparams, p::params, draws::Array{Float64,3})
    @unpack nind, nsim, nper = fp
    sh = shocks(zeros(nind*nper,nsim),zeros(nind*nper,nsim),zeros(nind*nper,nsim))
    sh.εw = p.σw.*draws[:,:,1]    
    sh.εmw = p.σmw.*draws[:,:,2]    
    sh.εtrch = draws[:,:,3]
    return sh
end

type datadict
    ix :: Array{Int64,2}
    ch :: Array{Int64,1}
    policy :: Array{Int64,1}
    steitc :: Array{Float64,2}
    fedeitc:: Array{Float64,2}
end


type datadict2
    ix::Array{Int64,2}
    ex::Array{Float64,1}
    b::Array{Float64,1}
    t::Array{Float64,1}
end

type datadict3
    ix :: Array{Int64,1}
    phthr :: Array{Float64,1}
    maxcr :: Array{Float64,1}
    phrate :: Array{Float64,1}
end

type datadict4
    ix::Array{Int64,2}
    b0::Array{Float64,1}
    b1::Array{Float64,1}
    b2::Array{Float64,1}
end

type datadict5
    b0::Array{Float64,1}
    b1::Array{Float64,1}
    b2::Array{Float64,1}
end

function initccp(fp::fparams,gender::String)	
    ccp0 = zeros(Float64,2,2,3) #ngpl=2 (grouped 2&3, 4&5), nkids = 2 (leave out 0), nts = 3 (grouped s = 3 and 4)
    for si = 1:3    	
    	ccp0[:,:,si] = readdlm("./trans_matrices/ccmatrix_$(gender)_s$(si).txt")[:,2:end]    	       
    end
    return ccp0
end



