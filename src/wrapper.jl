@views function objectivefunc!(initp0::Array{Float64,1},fp::fparams,moments::structureM,draws::Array{Float64,3},
					lEV::Array{single_1,3}, lEDU::Array{single_1,3}, lEVr::Array{Float64,3}, lEDUr::Array{Float64,3},
                    eitcDict::datadict, fgr::fgrids,
                    sim::Array{sim_t,2},ftaxp::datadict5,ctcp::datadict3,staxp::datadict4,ccp::Array{Float64,3},version::Int,gender::String)
    @unpack icEmtx = fgr
    
    sd = 2230
    srand(sd)  
    p0 = pv2p0(initp0)    
    finaln = quad(fp.nds,p0)
    #fp.ngu = size(finaln,1)
    #dgr = dgrids(fp,p0)
	       
    #random shocks
    sh = initshocks(fp,p0,draws)
    dp = dparams(fp,p0,finaln)
    #initializing the simulation structure
    sim = simClean!(sim,fp.nind,fp.nsim,fp.nper)    

    SolveSimAndCalcMom!(lEDU,lEV,lEDUr,lEVr,sim,moments,fp,p0,dp,eitcDict,fgr,finaln,ftaxp,ctcp,staxp,sh,ccp,version,gender)
        
    return moments.objfunc
    #writedlm("./logs/coeff.txt",initp0) 
end


function SolveSimAndCalcMom!(lEDU::Array{single_1,3},lEV::Array{single_1,3},lEDUr::Array{Float64,3},lEVr::Array{Float64,3},
                                sim::Array{sim_t,2},moments::structureM,fp::fparams,
                                p0::params,dp::dparams,eitcDict::datadict,fgr::fgrids,finaln::Array{Float64,1},ftaxp::datadict5,
                                ctcp::datadict3,staxp::datadict4,sh::shocks,ccp::Array{Float64,3},version::Int,gender::String)

    SolveAndSim_manySys!(lEDU,lEV,lEDUr,lEVr,sim,fp,p0,dp,eitcDict,fgr,finaln,ftaxp,ctcp,staxp,sh,ccp)    
    CalcSimMom!(sim,fp, moments,version,gender)
    return nothing 
end


function simwrapper(isx:: Int, sim::Array{sim_t,2},fp::fparams,sh::shocks,dp::dparams,eitcDict::datadict,
                        fgr::fgrids, policy::Int, ftaxp::datadict5, ctcp::datadict3, staxp::datadict4,lEV::Array{single_1,3},
                        lEDU::Array{single_1,3},lEVr::Array{Float64,3},lEDUr::Array{Float64,3},ccp::Array{Float64,3}, p0::params)
    @unpack isEmtx = fgr
    @unpack mina, maxa = fp
    
    rep = isEmtx.rep[isx]
    i = isEmtx.ind[isx]
    lb,ub = ilocate(i,fp.nper)
    age0 = max(sim[i,rep].sysA0[policy],mina)            
    age1 = min(sim[i,rep].sysA1[policy],maxa)                        
    s = sim[i,rep].s
    st = sim[i,rep].st
    ri = sim[i,rep].race
    irvw = sh.εw[lb:ub,rep]
    irvmw = sh.εmw[lb:ub,rep]                
    irvch = sh.εtrch[lb:ub,rep]
    if age1 >= age0
        #println("i $(i)", "rep $(rep)", "age0 $(age0)"," age1 $(age1)"," yob $(sim[i,rep].yob)", " initial ch $(sim[i,rep].aa.Oh.ch[age0-mina+1])", " policy $(policy)")
        #println("$(sim[i,rep].aa.Oh.ch[:])")
        sim[i,rep].aa = sim_a0toa1!(sim[i,rep].aa,age0,age1,ri,s,st,fp,sh,dp,eitcDict,fgr,policy,ftaxp,ctcp,staxp,lb,ub,irvw,irvmw,irvch,
                    lEV,lEDU,lEVr,lEDUr,ccp,p0)            
    #println("sim[1].Oh.ch  4 $(sim[1,1].aa.Oh.ch[:])")              
    end
    #println("sim[1].Oh.ch  5 $(sim[1,1].aa.Oh.ch[:])")       
    return sim       
end
    



function SolveAndSim_manySys!(lEDU::Array{single_1,3},lEV::Array{single_1,3},lEDUr::Array{Float64,3},lEVr::Array{Float64,3},
                                sim::Array{sim_t,2},fp::fparams,
                                p0::params,dp::dparams,eitcDict::datadict,fgr::fgrids,finaln::Array{Float64,1},ftaxp::datadict5,
                                ctcp::datadict3,staxp::datadict4,sh::shocks,ccp::Array{Float64,3})
    

    #solving the model
    #th = 1; sti = 1; ri = 1; s = 1; ai = 43; fei = 1; pei = 1; ki = 3; l1 = 1; cri=1
    #@time get_cEmtr!(lEDU, lEV, th, sti, ri, s, ai, fei, pei, ki, l1, cri, fp, p0, gr, eitcDict, trkids, policy)  

    for policy = 1:fp.npolicy
        println("policy $(policy)")
        #println("1 policy $(policy)","-- $(sim[1,1].aa.Oh.ch[:])")
        solve!(lEDU,lEV,lEDUr,lEVr,fp,p0,dp,eitcDict,fgr,policy,finaln,ftaxp,ctcp,staxp,ccp)

        println("solution done")
        #println("2 policy $(policy)","-- $(sim[1,1].aa.Oh.ch[:])")


		##SECOND PART TO PARALLELIZE: SIMULATION
        #sim[:] = @sync @parallel (vcat) for isx = 1:fp.totsi        
        for isx = 1:fp.totsi
            simwrapper(isx,sim,fp,sh,dp,eitcDict,fgr,policy,ftaxp,ctcp,staxp,lEV,lEDU,lEVr,lEDUr,ccp,p0)
       	end                        
    end



    #@sync @parallel (vcat) for icx = 1:6
    #for icx = 1:fp.totcti
    #    #println("icx $icx")          
    #    lEDU[icx].m0, lEV[icx].m0 = solveWL(lEDU[icx].m0,lEV[icx].m0,icx,fp,p0,dp,eitcDict,fgr,policy,finaln,ftaxp,ctcp,staxp)
    #     # lEDU[icx].m0, lEV[icx].m0 = solveWL(lEDU[icx].m0,lEV[icx].m0,icx,fp,p0,dp,dgr,eitcDict,fgr,policy)
    #end

    #=#Simulation starts here
    for rep = 1:fp.nsim
        for i = 1:fp.nind
            #println("i $i", "rep $rep")
            lb,ub = ilocate(i,fp.nper)
            irvw = sh.εw[lb:ub,rep]
            irvmw = sh.εmw[lb:ub,rep]
            #irvma  = sh.εma[lb:ub,rep]
            irvch = sh.εtrch[lb:ub,rep]
            ri = sim[i,rep].race
            s = sim[i,rep].s
            st = sim[i,rep].st
            #th = sim[i,rep].th
            #icx = icEmtx.ix[th,st,ri,s]
            icx = icEmtx.ix[st,ri,s]
            for ai = 1:fp.nper    
                #println("ai $ai")                  
                (sim[i,rep].aa = getopt!(sim[i,rep].aa, ri, s, st, ai, irvw[ai] ,irvmw[ai], #irvma[ai], 
                    irvch[ai],lEDU[icx].m0, lEV[icx].m0, fp, p0, dp, eitcDict, fgr, policy, ftaxp, ctcp,staxp))
            end
            #subroutine sim_a0toa1 ends here. (you can ignore this, this comment is here to )
        end
        #can potentially have writesimstata here.
    end=#

    return nothing
end


