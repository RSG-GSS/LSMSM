##########################################################################   
##### 10/13/2017: Modified by Miura, Hirotaka <Hirotaka.Miura@ny.frb.org>
##### 10/13/2017: Previously modified by Miura, Hirotaka <Hirotaka.Miura@ny.frb.org>
##### 10/13/2017: Created by Miura, Hirotaka <Hirotaka.Miura@ny.frb.org>
##########################################################################
##### Description: 
##### 	- Non-inplace version of getopt!().
##### Modifications:
#####	10/13/2017: 
#####		- Duplicated from solvesim2.jl.
#####		- Modularized version.
##########################################################################
##### Define module name.
module Msolvesim2
##### Import packages.
using Parameters		##### For @unpack.
using Mparams				##### For params.
using Mgrids				##### For grids.
using Mfuncs				##### For util.
using Msolvesim			##### For solveEE.
##### Specify names to be exported.
export
getopt2

##### Main routine for simulation.
function getopt2(sim::sim_t, rep::Int, i::Int, a::Int, ai::Int, irvw::Array{Float64,1}, irvu::Array{Float64,1}, lEDU::Array{single_1,1}, lEV::Array{single_1,1}, fp::fparams, p0::params, gr::grids)
    @unpack θ, β, R, ngpl, maxa, ngpk = fp
    @unpack αw0, αw1, αw2, αc = p0
    @unpack ge, ngpe, gk = gr

    #declarations needed for the simulation part
    y = Array{Float64}(ngpl)    
    e1val = Array{Int}(ngpl)  

    v = Array{Float64}(ngpl)
    c = Array{Float64}(ngpl)
    k1 = Array{Float64}(ngpl)
    y[1] = 0.
    #starting from here we're now in the subroutine SimValue
    if a<fp.maxa
        ix = ai+1
    end                                
    y[2] = exp(αw0 + αw1*log(sim.aa.Oh.e[ai]+1)+αw2*(convert(Float64,sim.aa.Oh.l_1[ai])-1.) + irvw[ai])
    for l = 1:ngpl                    
        e1val[l] = sim.aa.Oh.e[ai] + (l - 1)        
        #println("ai $ai","\t l $l", "\t e1val $(e1val[l])")
        if a == maxa
            c[l] = max(R*sim.aa.Oh.k[ai]+y[l],0.001)
            v[l] = util(c[l], l-1, a, irvu[ai], p0, fp)
        else
            #bounds on next period's experience
            #println("ai+1 $(ai+1)","\t gr.ngpe[ai+1] $(gr.ngpe[ai+1])")
            #println("e1val $(e1val[l])")
            #println("first gr.ge $(gr.ge)")
            lbe1, ube1 = bounds(e1val[l],ge[1:ngpe[ai+1],ai+1],ngpe[ai+1])
            #println("lbe1 $lbe1", " ube1 $ube1")
            lbe1val = ge[lbe1,ai+1]
            ube1val = ge[ube1,ai+1]
            #optimal consumption and V conditional on labor supply
            eulerArr = LinInterp1d.(e1val[l],lbe1val,ube1val,lEDU[ix].m0[l,:,lbe1],lEDU[ix].m0[l,:,ube1])
            eulerArr = (((β*R)^(1./(θ-1.)))/exp(αc*(l.-1.)/(θ-1.)))*eulerArr - (R*sim.aa.Oh.k[ai]+y[l]-gk[:]) 
            #if ai == 1
            #    println("e $(sim.aa.Oh.e[ai])"," k $(sim.aa.Oh.k[ai])")
            #    println("irvw $(irvw[ai])", "irvu $(irvu[ai])")
            #    println("eulerArr $eulerArr")
            #end
            k1[l], lbk, ubk = solveEE(ngpk, eulerArr, gk)
            #println("k1 $k1[l]", "lbk $lbk", "ubk $ubk")
            c[l] = max(R*sim.aa.Oh.k[ai]+y[l]-k1[l],0.001)
            v[l] = (util(c[l], l-1, a, irvu[ai], p0, fp) + β*calcEVinvinv(LinInterp2d(k1[l],gk[lbk],gk[ubk],e1val[l],lbe1val,ube1val, lEV[ix].m0[l,lbk,lbe1], lEV[ix].m0[l,lbk,ube1],lEV[ix].m0[l,ubk,lbe1],lEV[ix].m0[l,ubk,ube1]),fp))
            #if ai == 1
            #    println("v $(v[l])")
            #end
        end
    end 
    #subroutine SimValue ends here.
    ###CHECK THESE
    #1) whether you need lbe1+1 and ube1+1's in the v[l]
    #2) whether instead of l1 you should write l as you already did
    #4) in eulerArr, second line, whether you should use l or l-1
    sim.aa.Oh.V[ai], act = findmax(v)    #act[1]=optV, act[2] = optl
    #println("act $act","\t e1val(act) $(e1val[act])")
    sim.aa.Ch.l[ai] = act
    if a<fp.maxa
        sim.aa.Oh.l_1[ai+1] = act
    end
    sim.aa.Oh.w[ai] = y[act]
    sim.aa.Oh.inc[ai] = y[act]
    sim.aa.Ch.c[ai] = c[act]
    if a<fp.maxa
        sim.aa.Oh.k[ai+1] = k1[act]
        sim.aa.Oh.e[ai+1] = e1val[act]
    end    
    return(sim)   
end

##### End module definition.
end