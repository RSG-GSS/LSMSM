##########################################################################   
##### 10/13/2017: Modified by Miura, Hirotaka <Hirotaka.Miura@ny.frb.org>
##### 00/00/0000: Previously modified by 
##### 10/13/2017: Created by Kosar, Gizem <Gizem.Kosar@ny.frb.org>
##########################################################################
##### Description: 
##### 	- Module for executing simulation.
##### Modifications:
##########################################################################
##### Define module name.
module Mwrapper
##### Import packages.
using Mparams								##### For params.
using Mgrids								##### For grids.
using Distributions: Normal	##### For Normal().
using Msolvesim							##### For get_cEmtr!.
using Mmoments							##### For CalcSimMom!.
##### Specify names to be exported.
export
objectivefunc!

function objectivefunc!(initp0::Array{Float64,1},fp::fparams,moments::structureM)
    sd = 2230
    srand(sd)  
    p0 = pv2p0(initp0)
    for i in 1:length(initp0)
        println("$(fieldnames(p0)[i])", "\t $(initp0[i])" )
    end
    gr = grids(fp,p0)    
    #initializing the EV and EDU arrays
    lEV, lEDU = initEVEDU(fp)
    #initializing the simulation structure
    sim = initsim(fp)
    #drawing the random components for the simulation -- assuming no correlation
    derrw = Normal(0,p0.σw)
    derru = Normal(0,p0.σu)
    rvw = rand(derrw, fp.totn, fp.nsim)
    rvu = rand(derru, fp.totn, fp.nsim)
    #solving the model
    for a = fp.nper:(-1):1
        println("a $(a)")                
        for ei = 1:gr.ngpe[a] 
            #println("ei $ei")                          
            for ki = 1:fp.ngpk
                #println("ki $ki")
                for l1 = 1:fp.ngpl      
                    #println("l1 $l1")  
                    #=lEV[a].m0[l1,ki,ei], lEDU[a].m0[l1,ki,ei] = get_cEmtr(trp, ngpk, nper, ngu, ngpl, a, k, ei, l1, p0, e1, e1lb, e1ub, ge, gk, gεw, gεu, lEDU[a+1].m0[l1,:,:], lEV[a+1].m0[l1,:,:])=#
                    get_cEmtr!(lEDU, lEV, a, ei, ki, l1, fp, p0, gr)                                    
                end
            end
        end
    end

    #Simulation starts here
    for rep = 1:fp.nsim
        for i = 1:fp.nind
            #println("i $i", "rep $rep")
            irvw, irvu = ilocate(i, rep, rvw, rvu, fp)
            for ai = 1:fp.nper
                a = ai-1+fp.mina
                getopt!(sim, rep, i, a, ai, irvw, irvu, lEDU, lEV, fp, p0, gr)                
            end
            #subroutine sim_a0toa1 ends here. (you can ignore this, this comment is here to )
        end
        #can potentially have writesimstata here.
    end    
    CalcSimMom!(sim,fp, moments)
    difference = moments.dtamom - moments.simmom
    obj = sum((difference[moments.withvar.==1].^2).*(1./moments.wgtcov))
    #save("./sampleresults.jld", "lEV", lEV, "lEDU", lEDU, "p0", p0, "gr", gr, "sim", sim, "moments", moments, "obj", obj, "sd", sd)
    return obj
end

##### End module definition.
end