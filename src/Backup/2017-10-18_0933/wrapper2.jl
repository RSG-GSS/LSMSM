##########################################################################
##### User: Hirotaka Miura.                                   
##### Position: Research Analytics Associate.                                            
##### Organization: Federal Reserve Bank of New York.
##########################################################################   
##### 10/16/2017: Modified.
##### 10/12/2017: Previously modified.
##### 10/12/2017: Created.
##### Description: 
##### 	- Define wrapper for getopt2() and parallel version of main routine.
##### Modifications:
#####	10/16/2017:
#####		- Add pmap with and without batch.
##########################################################################

##### Wrapper for getopt2().
function getopt2wrapper(sim,rep,i,a,rvw,rvu,lEDU,lEV,fp,p0,gr)
	irvw, irvu = ilocate(i, rep, rvw, rvu, fp)
	for ai = 1:fp.nper
		a = ai-1+fp.mina
		##### 9/29/17 B1HXM10: Set to non-inplace function.
		sim[i,rep]=getopt2(sim[i,rep], rep, i, a, ai, irvw, irvu, lEDU, lEV, fp, p0, gr);
		#getopt!(sim, rep, i, a, ai, irvw, irvu, lEDU, lEV, fp, p0, gr)                
	end
	return(sim[i,rep]);
end

##### Parallel version of main routine.
function objectivefuncp(initp0::Array{Float64,1},fp::fparams,moments::structureM,partype::String)
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
		println("Sim start")
		##### Execute if "parfor" specified.
		if partype=="parfor"
			##### Display prompt.
			println("Executing parfor on ",nworkers()," workers");
			##### Loop through number of simulations.
			for rep = 1:fp.nsim
				##### Parallelize using @parallel.
				sim[:,rep]=@sync @parallel (vcat) for i = 1:fp.nind
				getopt2wrapper(sim,rep,i,a,rvw,rvu,lEDU,lEV,fp,p0,gr)
				end
			end   
		##### Execute if "pmap" specified.
		elseif partype=="pmap"
			##### Display prompt.
			println("Executing pmap on ",nworkers()," workers");
			##### Loop through number of simulations.
			for rep = 1:fp.nsim
				##### Generate collections of input arguments.
				inputs=[[sim,rep,i,a,rvw,rvu,lEDU,lEV,fp,p0,gr] for i = 1:fp.nind];
				##### Parallelize using pmap without batching.
				sim[:,rep]=pmap((args)->getopt2wrapper(args...),inputs)
			end			
		##### Execute if "pmapbatch" specified.
		elseif partype=="pmapbatch"
			##### Calculate batch size.
			batchsize=ceil(Int,fp.nind/nworkers());
			##### Display prompt.
			println("Executing pmap on ",nworkers(),
			" workers with batch size of ",batchsize);
			##### Loop through number of simulations.
			for rep = 1:fp.nsim
				##### Generate collections of input arguments.
				inputs=[[sim,rep,i,a,rvw,rvu,lEDU,lEV,fp,p0,gr] for i = 1:fp.nind];
				##### Parallelize using pmap batch.
				sim[:,rep]=pmap((args)->getopt2wrapper(args...),inputs,batch_size=batchsize)
			end
		end	
    println("Sim end")
		
		#=
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
		=#
		
    CalcSimMom!(sim,fp, moments)
    difference = moments.dtamom - moments.simmom
    obj = sum((difference[moments.withvar.==1].^2).*(1./moments.wgtcov))
    #save("./sampleresults.jld", "lEV", lEV, "lEDU", lEDU, "p0", p0, "gr", gr, "sim", sim, "moments", moments, "obj", obj, "sd", sd)
    return obj
end