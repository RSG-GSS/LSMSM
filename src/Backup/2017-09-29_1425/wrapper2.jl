##########################################################################
##### User: Hirotaka Miura.                                   
##### Position: Research Analytics Associate.                                            
##### Organization: Federal Reserve Bank of New York.
##########################################################################   
##### 09/29/2017: Modified.
##### 09/29/2017: Previously modified.
##### 09/29/2017: Created.
##### Description: 
##### 	- File that includes the main routine for LSMSM.
##### Modifications:
##### 		09/29/2017: 
#####			- Duplicated from wrapper.jl.
#####			- Adjust code to use getopt() instead of getopt!().
##########################################################################
function objectivefunc!(p::Array{Float64,1},fp::fparams,moments::structureM)
    sd = 2230
    srand(sd)  
    p0 = pv2p0(p)
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
        #println("a $(a)")                
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
		println("what1")
    #Simulation starts here
    for rep = 1:fp.nsim
        for i = 1:fp.nind
            #println("i $i", "rep $rep")
            irvw, irvu = ilocate(i, rep, rvw, rvu, fp)
            for ai = 1:fp.nper
                a = ai-1+fp.mina
                ##### 9/29/17 B1HXM10: Set to non-inplace function.
								sim=getopt(sim, rep, i, a, ai, irvw, irvu, lEDU, lEV, fp, p0, gr);
								#getopt!(sim, rep, i, a, ai, irvw, irvu, lEDU, lEV, fp, p0, gr)                
            end
            #subroutine sim_a0toa1 ends here. (you can ignore this, this comment is here to )
        end
        #can potentially have writesimstata here.
    end    
    println("what2")
		CalcSimMom!(sim,fp, moments)
    difference = moments.dtamom - moments.simmom
    obj = sum((difference[moments.withvar.==1].^2).*moments.wgtcov)
    #save("./sampleresults.jld", "lEV", lEV, "lEDU", lEDU, "p0", p0, "gr", gr, "sim", sim, "moments", moments, "obj", obj, "sd", sd)
    return obj
end