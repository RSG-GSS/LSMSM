##########################################################################
##### User: Hirotaka Miura.                                   
##### Position: Research Analytics Associate.                                            
##### Organization: Federal Reserve Bank of New York.
##########################################################################   
##### 09/29/2017: Modified.
##### 09/29/2017: Previously modified.
##### 09/29/2017: Created.
##### Description: 
##### 	- File that defines solvers for LSMSM.
##### Modifications:
##### 		09/29/2017: 
#####			- Duplicated from solvesim.jl.
#####			- Convert getopt!() into a non-inplace function.
##########################################################################
#solving the euler equation (for consumption)
function solveEE(ngpk::Int, eulerArr::Array{Float64,1}, gk::Array{Float64,1}) 
    if eulerArr[1]>0.
        lbk = 1
        ubk = 1
    else
        ubk = 2
        while true
            if eulerArr[ubk]>0. || ubk == ngpk
                break
            end
            ubk += 1
        end                          
        lbk = ubk - 1
        while true
            if eulerArr[lbk]<eulerArr[ubk]
                break
            end
            if ubk==ngpk
                break
            end
            ubk = ubk+1
        end
    end
    k1 = max(0.,LinInterp1d(0., eulerArr[lbk], eulerArr[ubk], gk[lbk], gk[ubk]))    
    return k1, lbk, ubk
end




#main routine to solve for the value functions
function get_cEmtr!(lEDU::Array{single_1,1}, lEV::Array{single_1,1}, a::Int, ei::Int, ki::Int, l_1::Int, fp::fparams, p0::params, gr::grids)
    @unpack β, R, θ, ngpk, ngu, ngpl, nper, maxngpe = fp
    @unpack ngpe, ge, gk, gεu, gεw, trp = gr
    @unpack αw0, αw1, αw2, αc = p0

    y = Array{Float64}(ngpl)
    c = Array{Float64}(ngpl)
    v = Array{Float64}(ngpl)
    eulerArr = Array{Float64}(ngpk)
    e1lb = zeros(Int64,ngpl,maxngpe)
    e1ub = zeros(Int64,ngpl,maxngpe)
    e1 = zeros(Int64,ngpl,maxngpe)

    optV = Array{Float64}(ngu,ngu)
    optDU = Array{Float64}(ngu,ngu)

    l1 = convert(Float64,l_1) - 1.


    y[1] = 0.
    e = ge[ei,a]
    age = a+19
    k = gk[ki]
    if a<nper
        e1, e1lb, e1ub = nxtebds(a, ngpl, maxngpe, ngpe, ge)
    end

    for iu = 1:ngu
        for iv = 1:ngu
            y[2] = exp(αw0 + αw1*log(e+1)+αw2*l1 + gεw[iv]) #grids for the epsilon.
            for l = 1:ngpl
                #println("l $l")
                resources = R*k+y[l]
                if (a<nper)                    
                    e1val = e1[l,ei]
                    e1lbi = e1lb[l,ei]
                    e1ubi = e1ub[l,ei]
                    e1lbval = ge[e1lbi+1,a+1]
                    e1ubval = ge[e1ubi+1,a+1]
                    eulerArr = LinInterp1d.(e1val, e1lbval, e1ubval, lEDU[a+1].m0[l,:,e1lbi],lEDU[a+1].m0[l,:,e1ubi])
                    eulerArr = (((β*R)^(1./(θ-1)))/exp(αc*(l.-1.)/(θ-1)))*eulerArr - (resources-gk[:]) 
                    k1, lbk, ubk = solveEE(ngpk, eulerArr, gk)
                    c[l] = max(resources - k1,0.00001)
                    v[l] = (util(c[l], l-1, age, gεu[iu], p0, fp) + β*LinInterp2d(k1,gk[lbk],gk[ubk],e1val,e1lbval,e1ubval, lEV[a+1].m0[l,lbk,e1lbi], lEV[a+1].m0[l,lbk,e1ubi],lEV[a+1].m0[l,ubk,e1lbi],lEV[a+1].m0[l,ubk,e1ubi]))
                else
                    c[l] = max(resources, 0.00001)
                    #println("resources $resources")
                    v[l] = util(c[l], l-1, age, gεu[iu], p0, fp) 
                    #println("util $(util(c[l], l-1, age, gεu[iu], p0))")                   
                end
            end
            #println("v $v","\t c $c")
            optV[iv,iu], act = findmax(v)
            optDU[iv,iu] = dudc(c[act], act-1, p0, fp)
        end
    end
    #println("optV $optV")
    lEV[a].m0[l_1,ki,ei] = sum([trp[iu]*optV[iu] for iu in eachindex(trp)])
    lEDU[a].m0[l_1,ki,ei] = sum([trp[iu]*optDU[iu] for iu in eachindex(trp)])
    return nothing
end


#main routine for the simulation. 
function getopt(sim::Array{sim_t,2}, rep::Int, i::Int, a::Int, ai::Int, irvw::Array{Float64,1}, irvu::Array{Float64,1}, lEDU::Array{single_1,1}, lEV::Array{single_1,1}, fp::fparams, p0::params, gr::grids)
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
    y[2] = exp(αw0 + αw1*log(sim[i,rep].aa.Oh.e[ai]+1)+αw2*(convert(Float64,sim[i,rep].aa.Oh.l_1[ai])-1.) + irvw[ai])
    for l = 1:ngpl                    
        e1val[l] = sim[i,rep].aa.Oh.e[ai] + (l - 1)        
        #println("ai $ai","\t l $l", "\t e1val $e1val[l]")
        if a == maxa
            c[l] = max(R*sim[i,rep].aa.Oh.k[ai]+y[l],0.00001)
            v[l] = util(c[l], l-1, a, irvu[ai], p0, fp)
        else
            #bounds on next period's experience
            #println("ai+1 $(ai+1)","\t gr.ngpe[ai+1] $(gr.ngpe[ai+1])")
            #println("e1val $(e1val[l])")
            #println("first gr.ge $(gr.ge)")
            lbe1, ube1 = bounds(e1val[l],ge[1:ngpe[ai+1],ai+1],ngpe[ai+1])
            lbe1val = ge[lbe1,ai+1]
            ube1val = ge[ube1,ai+1]
            #optimal consumption and V conditional on labor suppy
            eulerArr = LinInterp1d.(e1val[l],lbe1val,ube1val,lEDU[ix].m0[l,:,lbe1],lEDU[ix].m0[l,:,ube1])
            eulerArr = (((β*R)^(1./(θ-1)))/exp(αc*(l.-1.)/(θ-1.)))*eulerArr - (R*sim[i,rep].aa.Oh.k[ai]+y[l]-gk[:]) 
            k1[l], lbk, ubk = solveEE(ngpk, eulerArr, gk)
            c[l] = max(R*sim[i,rep].aa.Oh.k[ai]+y[l]-k1[l],0.00001)
            v[l] = (util(c[l], l-1, a, irvu[ai], p0, fp) + β*LinInterp2d(k1[l],gk[lbk],gk[ubk],e1val[l],lbe1val,ube1val, lEV[ix].m0[l,lbk,lbe1], lEV[ix].m0[l,lbk,ube1],lEV[ix].m0[l,ubk,lbe1],lEV[ix].m0[l,ubk,ube1]))
        end
    end 
    #subroutine SimValue ends here.
    ###CHECK THESE
    #1) whether you need lbe1+1 and ube1+1's in the v[l]
    #2) whether instead of l1 you should write l as you already did
    #4) in eulerArr, second line, whether you should use l or l-1
    sim[i,rep].aa.Oh.V[ai], act = findmax(v)    #act[1]=optV, act[2] = optl
    #println("act $act","\t e1val(act) $(e1val[act])")
    sim[i,rep].aa.Ch.l[ai] = act
    if a<fp.maxa
        sim[i,rep].aa.Oh.l_1[ai+1] = act
    end
    sim[i,rep].aa.Oh.w[ai] = y[act]
    sim[i,rep].aa.Oh.inc[ai] = y[act]
    sim[i,rep].aa.Ch.c[ai] = c[act]
    if a<fp.maxa
        sim[i,rep].aa.Oh.k[ai+1] = k1[act]
        sim[i,rep].aa.Oh.e[ai+1] = e1val[act]
    end    
    ##### 9/29/17 B1HXM10: Return sim structure with calculated values.
		return(sim)    
		#return nothing   
end

##### 9/29/17 B1HXM10: 
#####		- Change getopt() to getopt2().
#####		- Change sim::Array{sim_t,2} to sim::sim_t.
#####		- Change sim[i,rep] to sim.
#main routine for the simulation. 
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
        #println("ai $ai","\t l $l", "\t e1val $e1val[l]")
        if a == maxa
            c[l] = max(R*sim.aa.Oh.k[ai]+y[l],0.00001)
            v[l] = util(c[l], l-1, a, irvu[ai], p0, fp)
        else
            #bounds on next period's experience
            #println("ai+1 $(ai+1)","\t gr.ngpe[ai+1] $(gr.ngpe[ai+1])")
            #println("e1val $(e1val[l])")
            #println("first gr.ge $(gr.ge)")
            lbe1, ube1 = bounds(e1val[l],ge[1:ngpe[ai+1],ai+1],ngpe[ai+1])
            lbe1val = ge[lbe1,ai+1]
            ube1val = ge[ube1,ai+1]
            #optimal consumption and V conditional on labor suppy
            eulerArr = LinInterp1d.(e1val[l],lbe1val,ube1val,lEDU[ix].m0[l,:,lbe1],lEDU[ix].m0[l,:,ube1])
            eulerArr = (((β*R)^(1./(θ-1)))/exp(αc*(l.-1.)/(θ-1.)))*eulerArr - (R*sim.aa.Oh.k[ai]+y[l]-gk[:]) 
            k1[l], lbk, ubk = solveEE(ngpk, eulerArr, gk)
            c[l] = max(R*sim.aa.Oh.k[ai]+y[l]-k1[l],0.00001)
            v[l] = (util(c[l], l-1, a, irvu[ai], p0, fp) + β*LinInterp2d(k1[l],gk[lbk],gk[ubk],e1val[l],lbe1val,ube1val, lEV[ix].m0[l,lbk,lbe1], lEV[ix].m0[l,lbk,ube1],lEV[ix].m0[l,ubk,lbe1],lEV[ix].m0[l,ubk,ube1]))
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
    ##### 9/29/17 B1HXM10: Return sim structure with calculated values.
		return(sim)    
		#return nothing   
end


#master routine that calls all the functions.
# function main()
#     p0 = params()       
#     fp = fparams()
#     moments = initmom(fp)
#     println("starting to calculate data moments and the bootstrapped moments")
#     CalcDataBootMom!(fp,moments)
#     println("starting to solve and simulate")
#     @time obj = objectivefunc!(p0,fp,moments)
#     return obj, moments
# end

# function objectivefunc(p0)
#     fp = fparams()
#     srand(2230)
#     gr = grids(fp,p0)
    
# main()
# #nperprintln("first done")
# #@time main() 
# obj, moments = main()


