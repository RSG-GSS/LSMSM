##########################################################################   
##### 04/05/2018: Modified by Hirotaka Miura.
##### 04/04/2018: Previously modified by Hirotaka Miura.
##### 04/04/2018: Created by Gizem Kosar.
##### Description: 
##### 	- Functions to solve simulation.
##### Modifications:
#####	04/04/2018: 
#####		- Begin making adjustments for parallelization. 
#####	04/04/2018: 
#####		- Continue development. 
##########################################################################

#solving the euler equation (for consumption)
@views function solveEE(ngpk::Int, eulerArr::Array{Float64,1}, gk::Array{Float64,1},ki::Int, boundk::Float64) 
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
    kEE = LinInterp1d(0., eulerArr[lbk], eulerArr[ubk], gk[lbk], gk[ubk])        
    if kEE >= gk[ki] + boundk
        return kEE,lbk,ubk
    else
        k1 = gk[ki]+boundk
        lbk1, ubk1 = bounds(k1,gk,ngpk)  
        return k1,lbk1,ubk1
    end
end

@views function solveEE(ngpk::Int, eulerArr::Array{Float64,1}, gk::Array{Float64,1},ki::Float64, boundk::Float64) 
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
    kEE = LinInterp1d(0., eulerArr[lbk], eulerArr[ubk], gk[lbk], gk[ubk])        
    if kEE >= ki + boundk
        return kEE,lbk,ubk
    else
        k1 = ki+boundk
        lbk1, ubk1 = bounds(k1,gk,ngpk)  
        return k1,lbk1,ubk1
    end
end

##### Hirotaka Miura: Wrapper for parallelizing call to solveWL!() in solve!().
function solveWLWrapper(
	icx::Int,
	lEDU::Array{single_1,2},
	lEV::Array{single_1,2},
	lEDUr::Array{Float64,3},
	lEVr::Array{Float64,3},
	fp::fparams,
	p0::params,
	dp::dparams,
	eitcDict::datadict,
	fgr::fgrids,
	policy::Int,
	finaln::Array{Float64,1},
	ftaxp::datadict5,
	ctcp::datadict3,
	staxp::datadict4,
	ccp::Array{Float64,3})
	
	@unpack icEmtx = fgr;
	
	s = icEmtx.s[icx];
	sti = icEmtx.st[icx];

	for l1 = 1:fp.ngpl1, cri = 1:fp.ngpeitc
		for ai = fp.nper:(-1):fp.as[s] - fp.mina+1, fei = 1:fgr.ngpe[ai]#, pei = 1:fgr.ngpe[ai]
				for ki = 1:fp.ngpk                
					#println("1 ai $(ai)")
						lEDU[sti,s].m0, lEV[sti,s].m0= solveWL!(lEDU[sti,s].m0,lEV[sti,s].m0,lEDUr,lEVr,sti,s,l1,cri,ai,fei,ki,fp,p0,dp,eitcDict,fgr,policy,finaln,ftaxp,ctcp,staxp,ccp)
					end
			end
	end	
	
	return lEDU[sti,s].m0, lEV[sti,s].m0, sti, s
	
end
	

function solve!(lEDU::Array{single_1,2},lEV::Array{single_1,2},lEDUr::Array{Float64,3},lEVr::Array{Float64,3},
                fp::fparams,p0::params,dp::dparams,
                eitcDict::datadict,fgr::fgrids,policy::Int,finaln::Array{Float64,1},ftaxp::datadict5,
                ctcp::datadict3,staxp::datadict4,ccp::Array{Float64,3})
    @unpack icEmtx = fgr    

    #solution for retirement years
    for s = 1:fp.nts
        solveR!(lEDUr,lEVr,fp,p0,dp,fgr,s,policy)
    end
		
		tic();
		##### Display prompt.
		println("Executing parfor on ",nworkers()," workers");
		
		#res1=Array{Array{Float64,5},1}(fp.totcti);
		#res2=Array{Array{Float64,5},1}(fp.totcti);
		##### Parallelize using @parallel.
		res=@sync @parallel (vcat) for i = 1:fp.totcti
			solveWLWrapper(i,lEDU,lEV,lEDUr,lEVr,fp,p0,dp,eitcDict,fgr,policy,finaln,ftaxp,ctcp,staxp,ccp)
		end 
		
		for i=1:length(res)
		
			lEDU[res[i][3],res[i][4]].m0=res[i][1];
			
			lEV[res[i][3],res[i][4]].m0=res[i][2];
			
		end

		toc();
		
		#=
    tic();
    ##FIRST PART TO PARALLELIZE: SOLUTION
    #@sync @parallel (vcat) for icx = 1:totcti
    for icx = 1:fp.totcti    
        s = icEmtx.s[icx]
        sti = icEmtx.st[icx]
        #l1 = icEmtx.l1[icx]
        #cri = icEmtx.cr[icx]
        for l1 = 1:fp.ngpl1, cri = 1:fp.ngpeitc
        	for ai = fp.nper:(-1):fp.as[s] - fp.mina+1, fei = 1:fgr.ngpe[ai]#, pei = 1:fgr.ngpe[ai]
            	for ki = 1:fp.ngpk                
            		#println("1 ai $(ai)")
                	lEDU[sti,s].m0, lEV[sti,s].m0= solveWL!(lEDU[sti,s].m0,lEV[sti,s].m0,lEDUr,lEVr,sti,s,l1,cri,ai,fei,ki,fp,p0,dp,eitcDict,fgr,policy,finaln,ftaxp,ctcp,staxp,ccp)
                end
            end
        end
    end
		toc();
		=#
		
    return nothing
end



function solveR!(lEDUr::Array{Float64,3},lEVr::Array{Float64,3},fp::fparams,p0::params,dp::dparams,fgr::fgrids, s::Int,policy::Int)
    @unpack β, R, θ, enda, mina, nper, ngpk, ngpch, maxa, minc = fp
    @unpack gk, trkids = fgr
    @unpack boundk, eepar, utilpar, dutilpar = dp    

    v = zeros(Float64,ngpch)
    du = zeros(Float64,ngpch)
    eulerArr = zeros(Float64,ngpk)      


    for a = enda:(-1):maxa+1
        ai = a-maxa
        for ch = 1:ngpch
            for ki = 1:ngpk                       
                if a == enda                
                    c = max(R*gk[ki],minc[policy,ch])
                    #println("c $(c)","minc $(minc[policy,ch])", "Rgk $(R*gk[ki])")
                    v[ch] = util(c, θ, utilpar[s,ch,1,1])
                    du[ch] = dudc(c,θ,dutilpar[s,ch,1,1])               
                else                
                    for i = 1:ngpk
                        eulerArr[i] = eepar[s,ch,1,1]*lEDUr[i,ch,s] - (R*gk[ki]-gk[i])
                    end
                    k1, lbk, ubk = solveEE(ngpk,eulerArr,gk,ki,boundk[a-mina+1,s])
                    c = max(R*gk[ki] - k1,minc[policy,ch])
                    v[ch] = (util(c, θ, utilpar[s,ch,1,1] + β*calcEVinvinv(LinInterp1d(k1,gk[lbk],gk[ubk],lEVr[lbk,ch,s],lEVr[ubk,ch,s]),fp)))
                    du[ch] = dudc(c,θ,dutilpar[s,ch,1,1])
                end                
                #integrate family composition
                for ch_1  = 1:ngpch
                    lEVr[ki,ch_1,s] = calcEVinv(trkids[ai,s,ch_1,:]'*v,fp)
                    lEDUr[ki,ch_1,s] = calcEDUinv(trkids[ai,s,ch_1,:]'*du,fp)
                end
            end
        end
    end
    return nothing
end



function solveWL!(lEDU::Array{Float64,5},lEV::Array{Float64,5},lEDUr::Array{Float64,3},lEVr::Array{Float64,3},
                            sti::Int,s::Int,l1::Int,cri::Int,
                            ai::Int,fei::Int,ki::Int,fp::fparams,p0::params,dp::dparams,eitcDict::datadict,
                            fgr::fgrids,policy::Int,finaln::Array{Float64,1},ftaxp::datadict5,ctcp::datadict3,staxp::datadict4,ccp::Array{Float64,3})
    @unpack ngpch,ngu,ngpl,ngpeitc, nper = fp
    @unpack geitc, trkids, fe1, fe1lb, fe1ub, ge, gk, iEmtx = fgr

    cEDU = zeros(Float64,fp.ngpch)
    cEV = zeros(Float64,fp.ngpch)


    fe = ge[fei,ai] #full-time experience - actual amount    
    #pe = ge[pei,ai] #full-time experience - actual amount    
    cr = geitc[cri] #the amount of eitc    
    k = gk[ki] #level of (1+r)*net worth + eitc + ctc


    y = wage(fe,[1,2,3,4,5],s,l1,p0,finaln)   
    if ai<nper        
        cr1, cr1lb, cr1ub = nxtcrbds(ngu,ngpl, ngpch, policy, sti, ngpeitc, geitc, y, eitcDict,ctcp)        
    end 
    cc = ccare(s,ccp) 
    tt = ctaxsys(y,sti,policy,ftaxp,staxp)


    fe1val = @view fe1[:,fei,ai]
    fe1lbi = @view fe1lb[:,fei,ai]
    fe1ubi = @view fe1ub[:,fei,ai]                    
    
    

    for ch = 1:ngpch
        if ai<nper
            ix = iEmtx.ix[ch,ai+1]          
            #ix_1 = iEmtx.ix[ch,ai+1]       
            #plEDU = @view lEDU[sti,s,ix].m0[:,:,:,:]
            #plEV = @view lEV[sti,s,ix].m0[:,:,:,:]
        	cEDU[ch],cEV[ch] = get_cEmtr!(lEDU[:,:,:,:,ix],lEV[:,:,:,:,ix],y,fei,fe1val,fe1lbi,fe1ubi,ki,k,
                                    l1,cri,cr,cr1,cr1lb,cr1ub,cc[:,ch],tt,fp,p0,dp,ge[:,ai],geitc,gk,s,ch,ai,policy)
        else        	        
        	cEDU[ch],cEV[ch] = get_cEmtr!(lEDUr[:,ch,s],lEVr[:,ch,s],y,ki,k,
                                    l1,cri,cr,cc[:,ch],tt,fp,p0,dp,ge[:,ai],gk,s,ch,ai,policy)
        end        
    end

    #integrate family
    trkidssub = @view trkids[ai,s,:,:]    
    uEDU, uEV = integrate_family(cEDU,cEV,ngpch, trkidssub)
    #println("cev $(cEV)","uev $(uEV)")
    #calcEinv
    for ch = 1:ngpch
        ix = iEmtx.ix[ch,ai]
        lEDU[ki,fei,cri,l1,ix], lEV[ki,fei,cri,l1,ix] = calcEinv!(uEDU[ch],uEV[ch],fp)
    end
    #println("lEV $(lEV)")
    return lEDU, lEV
end


#main routine to solve for the value functions
function get_cEmtr!(plEDU::Array{Float64,4}, plEV::Array{Float64,4},
                        y::Array{Float64,2},
                        fei::Int, fe1::SubArray{Int64,1,Array{Int64,3}},fe1lb::SubArray{Int64,1,Array{Int64,3}},
                        fe1ub::SubArray{Int64,1,Array{Int64,3}},
                        ki::Int64,k::Float64,l1::Int64,cri::Int64,cr::Float64,cr1::Array{Float64,3},cr1lb::Array{Int64,3}, cr1ub::Array{Int64,3},
                        cc::Array{Float64,1},tt::Array{Float64,2},fp::fparams, p0::params, dp::dparams,ge::Array{Int64,1},
                        geitc::Array{Float64,1},gk::Array{Float64,1},s::Int64,ch::Int64,ai::Int64,policy::Int)
    @unpack β, R, θ, ngpk, ngu, ngpl, nper, ngpch, ngpeitc, ngpt, wgts, hrs, minc = fp
    @unpack boundk, eepar, utilpar, dutilpar = dp    
    #@unpack ge, geitc, gk = fgr

    
    c = zeros(Float64,ngpl,ngpt)
    v = zeros(Float64,ngpl,ngpt)
    du = zeros(Float64,ngpl,ngpt)
    eulerArr = zeros(Float64,ngpk)      
    optV = 0.
    optDU = 0.
    tempv = 0.
    tempc = 0.
    tempdu = 0.
    
    #println("ai $(ai)")
    if ai < nper      
        for pt = 1:ngpt                           
            for l = 1:ngpl
                fe1val = fe1[l]
                #println("fe1val $(fe1val)")
                fe1lbi = fe1lb[l]
                #println("fe1lbi $(fe1lbi)")
                fe1ubi = fe1ub[l]
                #println("fe1ubi $(fe1ubi)")
                fe1lbval = ge[fe1lbi]
                fe1ubval = ge[fe1ubi]                        
                #pe1val = pe1[l]
                #println("pe1val $(pe1val)")
                #pe1lbi = pe1lb[l]
                #println("pe1lbi $(pe1lbi)")                    
                #pe1ubi = pe1ub[l]
                #println("pe1ubi $(pe1ubi)")
                #pe1lbval = ge[pe1lbi]
                #pe1ubval = ge[pe1ubi]                        
                if l == 1
                    l1 = 1
                    resources = R*k + convert(Float64,pt-1)*cr - cc[l]
                    cr1val = 0.
                    cr1lbi = 1
                    cr1ubi = 1
                    cr1lbval = 0.
                    cr1ubval = 0.
                    for i =1:ngpk                               
                        eulerArr[i] = (LinInterp1d(convert(Float64,fe1val),convert(Float64,fe1lbval),convert(Float64,fe1ubval),                            
                            plEDU[i,fe1lbi,1,l1],plEDU[i,fe1ubi,1,l1]))
                        eulerArr[i] = eepar[s,ch,l,pt]*eulerArr[i] - (resources-gk[i])
                    end
                    k1, lbk, ubk = solveEE(ngpk,eulerArr,gk,ki,boundk[ai,s])
                    c[l,pt] = max(resources - k1,minc[policy,ch])    
                    v[l,pt] = (util(c[l,pt], θ, utilpar[s,ch,l,pt]) 
                        + (β*calcEVinvinv(LinInterp2d(k1,gk[lbk],gk[ubk],
                        convert(Float64,fe1val),convert(Float64,fe1lbval),convert(Float64,fe1ubval),                            
                        plEV[lbk,fe1lbi,1,l1],plEV[lbk,fe1ubi,1,l1],
                        plEV[ubk,fe1lbi,1,l1],plEV[ubk,fe1ubi,1,l1]),fp)))                            
                    du[l,pt] = dudc(c[l,pt],θ,dutilpar[s,ch,l,pt]) 
                else
                    if l==2 || l==3
                        l1 = 2
                    else
                        l1 = 3
                    end
                    tempc = 0.
                    tempv = 0.
                    tempdu = 0.
                    for iv = 1:ngu
                        resources = R*k+y[l,iv]*hrs[l] + convert(Float64,pt-1)*cr - cc[l] + tt[l,iv]                                         
                        cr1val = cr1[iv,l,ch]
                        cr1lbi = cr1lb[iv,l,ch]
                        cr1ubi = cr1ub[iv,l,ch]
                        cr1lbval = geitc[cr1lbi]
                        cr1ubval = geitc[cr1ubi]
                        for i =1:ngpk                               
                            eulerArr[i] = (LinInterp2d(cr1val,cr1lbval,cr1ubval,
                                convert(Float64,fe1val),convert(Float64,fe1lbval),convert(Float64,fe1ubval),                                    
                                plEDU[i,fe1lbi,cr1lbi,l1],plEDU[i,fe1ubi,cr1lbi,l1],                                
                                plEDU[i,fe1lbi,cr1ubi,l1],plEDU[i,fe1ubi,cr1ubi,l1]))
                            eulerArr[i] = eepar[s,ch,l,pt]*eulerArr[i] - (resources-gk[i])
                        end
                        k1, lbk, ubk = solveEE(ngpk,eulerArr,gk,ki,boundk[ai,s])
                        civ = max(resources - k1,minc[policy,ch])
                        tempc += (wgts[iv]*civ)
                        tempv += (wgts[iv]*(util(civ, θ, utilpar[s,ch,l,pt]) 
                            + (β*calcEVinvinv(LinInterp3d(k1,gk[lbk],gk[ubk],cr1val,cr1lbval,cr1ubval,    
                            convert(Float64,fe1val),convert(Float64,fe1lbval),convert(Float64,fe1ubval),
                            plEV[lbk,fe1lbi,cr1lbi,l1],plEV[ubk,fe1lbi,cr1lbi,l1],
                            plEV[lbk,fe1lbi,cr1ubi,l1],plEV[ubk,fe1lbi,cr1ubi,l1],
                            plEV[lbk,fe1ubi,cr1lbi,l1],plEV[ubk,fe1ubi,cr1lbi,l1],
                            plEV[lbk,fe1ubi,cr1ubi,l1],plEV[ubk,fe1ubi,cr1ubi,l1]),fp))))
                        tempdu += (wgts[iv]*dudc(civ,θ,dutilpar[s,ch,l,pt])) 
                    end
                    c[l,pt] = tempc
                    v[l,pt] = tempv
                    du[l,pt] = tempdu                    
                end
            end
        end            
    end #if ai<nper
    optV, act = findmax(v)
    optl,optpt = ind2sub(v,act)
    optDU = du[act]    

    return optDU, optV
end


#main routine to solve for the value functions
function get_cEmtr!(plEDUr::Array{Float64,1}, plEVr::Array{Float64,1},
                        y::Array{Float64,2},                        
                        ki::Int64,k::Float64,l1::Int64,cri::Int64,cr::Float64,
                        cc::Array{Float64,1},tt::Array{Float64,2},fp::fparams, p0::params, dp::dparams,ge::Array{Int64,1},
                        gk::Array{Float64,1},s::Int64,ch::Int64,ai::Int64,policy::Int)
    @unpack β, R, θ, ngpk, ngu, ngpl, nper, ngpch, ngpeitc, ngpt, wgts, hrs, minc = fp
    @unpack boundk, eepar, utilpar, dutilpar = dp    
    #@unpack ge, geitc, gk = fgr

    
    c = zeros(Float64,ngpl,ngpt)
    v = zeros(Float64,ngpl,ngpt)
    du = zeros(Float64,ngpl,ngpt)
    eulerArr = zeros(Float64,ngpk)      
    optV = 0.
    optDU = 0.
    tempv = 0.
    tempc = 0.
    tempdu = 0.
    
    #println("ai $(ai)")
    for pt = 1:ngpt
        for l = 1:ngpl
            if l == 1                        
                resources = R*k + convert(Float64,pt-1)*cr - cc[l]
                for i = 1:ngpk
                    eulerArr[i] = eepar[s,ch,l,pt]*lEDUr[i] - (resources - gk[i])
                end
                k1,lbk,ubk = solveEE(ngpk,eulerArr,gk,ki,boundk[ai,s])
                c[l,pt] = max(resources-k1,minc[policy,ch])
                v[l,pt] = (util(c[l,pt], θ, utilpar[s,ch,l,pt]) + β*calcEVinvinv(LinInterp1d(k1,gk[lbk],gk[ubk],
                    lEVr[lbk],lEVr[ubk]),fp))
                du[l,pt] = dudc(c[l,pt],θ,dutilpar[s,ch,l,pt])
                #println("3 c $(c[l,pt])","v $(v[l,pt])","du $(du[l,pt])")
            else
                tempc = 0.
                tempv = 0.
                tempdu = 0.
                for iv = 1:ngu
                    resources = R*k+y[l,iv]*hrs[l] + convert(Float64,pt-1)*cr - cc[l] + tt[l,iv]    
                    for i = 1:ngpk
                        eulerArr[i] = eepar[s,ch,l,pt]*lEDUr[i] - (resources - gk[i])
                    end
                    k1,lbk,ubk = solveEE(ngpk,eulerArr,gk,ki,boundk[ai,s])
                    civ = max(resources-k1,minc[policy,ch])
                    tempc += (wgts[iv]*civ)
                    tempv += (wgts[iv]*(util(civ, θ, utilpar[s,ch,l,pt])
                        +β*calcEVinvinv(LinInterp1d(k1,gk[lbk],gk[ubk],lEVr[lbk],lEVr[ubk]),fp))) 
                    tempdu += (wgts[iv]*dudc(civ,θ,dutilpar[s,ch,l,pt])) 
                end                                     
                c[l,pt] = tempc
                v[l,pt] = tempv
                du[l,pt] = tempdu
                #println("4 c $(c[l,pt])","v $(v[l,pt])","du $(du[l,pt])")
            end #if l == 1
        end #l loop
    end #pt loop    
    optV, act = findmax(v)
    optl,optpt = ind2sub(v,act)
    optDU = du[act]    

    return optDU, optV
end


function integrate_family(cEDU::Array{Float64,1},cEV::Array{Float64,1},ngpch::Int64, trkids::SubArray{Float64,2,Array{Float64,4}})
    uEV = zeros(Float64,ngpch)
    uEDU = zeros(Float64,ngpch)
    for ch_1 = 1:ngpch                
        uEV[ch_1] = trkids[ch_1,:]'*cEV
        uEDU[ch_1] = trkids[ch_1,:]'*cEDU        
    end
    return uEDU, uEV
end

function calcEinv!(uEDU::Float64,uEV::Float64,fp::fparams)
    plEDU_1 = calcEDUinv(uEDU, fp)
    plEV_1 = calcEVinv(uEV, fp)
    return plEDU_1, plEV_1
end

                  
function sim_a0toa1!(sim::OhCh_t,age0::Int,age1::Int,ri::Int,s::Int,st::Int,fp::fparams,sh::shocks,dp::dparams,eitcDict::datadict,
                        fgr::fgrids, policy::Int, ftaxp::datadict5, ctcp::datadict3, staxp::datadict4,lb::Int,ub::Int,
                        Arrirvw::Array{Float64,1},Arrirvmw::Array{Float64,1},Arrirvch::Array{Float64,1},lEV::Array{single_1,2},
                        lEDU::Array{single_1,2},lEVr::Array{Float64,3},lEDUr::Array{Float64,3},ccp::Array{Float64,3}, p0::params)
    @unpack as, mina = fp
    @unpack iEmtx = fgr
    #println("sim[1].Oh.ch  1 $(sim.Oh.ch[:])")

    for a = max(as[s],age0):age1
        ai = a-mina+1        
        irvw = Arrirvw[ai]
        irvmw = Arrirvmw[ai]
        irvch = Arrirvch[ai]
        ch = sim.Oh.ch[ai]
        #if ch == 0
        #    println("a $(a)"," s $(s)"," policy $(policy)")
        #end
        ix = iEmtx.ix[ch,ai]
        sim = getopt!(sim, s, st, ai, irvw, irvmw, irvch, lEDU[st,s].m0[:,:,:,:,ix], lEV[st,s].m0[:,:,:,:,ix], lEDUr[:,:,s], lEVr[:,:,s],fp, p0, dp,eitcDict, fgr, policy, ftaxp, ctcp, staxp, ccp)
    end
    #println("sim[1].Oh.ch  3 $(sim.Oh.ch[:])")
    return sim
end


#main routine for the simulation. 
function getopt!(sim::OhCh_t, s::Int, st::Int, ai::Int, irvw::Float64, irvmw::Float64, #irvma::Float64, 
            irvch::Float64, lEDU::Array{Float64,4}, lEV::Array{Float64,4}, lEDUr::Array{Float64,2}, lEVr::Array{Float64,2}, 
            fp::fparams, p0::params, dp::dparams,
            eitcDict::datadict, fgr::fgrids, policy::Int, ftaxp::datadict5, ctcp::datadict3, staxp::datadict4, ccp::Array{Float64,3})
    @unpack β, R, θ, ngpk, ngu, ngpl, nper, ngpch, ngpeitc, ngpt, wgts, hrs,minc = fp
    @unpack boundk, eepar, utilpar = dp    
    @unpack trkids, ngpe, ge, geitc, gk, iEmtx, fe1, fe1lb, fe1ub = fgr

    a = ai-1+fp.mina    
    #declarations needed for the simulation part
    y = Array{Float64}(ngpl)    
    fe1val = Array{Int64}(ngpl)
    pe1val = Array{Int64}(ngpl)
    v = Array{Float64}(ngpl,ngpt)
    c = Array{Float64}(ngpl,ngpt)
    k1 = Array{Float64}(ngpl,ngpt)
    eulerArr = Array{Float64}(ngpk)
    

    y = wage(sim.Oh.fte[ai],[1,2,3,4,5],s,sim.Oh.l_1[ai],p0,irvw)
    tt = ctaxsys(y,st,policy,ftaxp,staxp)    
    if sim.Oh.ch[ai]==1
        cc = zeros(Float64,5)
    else
        cc = ccare(sim.Oh.ch[ai],s,ccp)
    end

    #starting from here we're now in the subroutine SimValue
    if ai<nper
        #ix = iEmtx.ix[sim.Oh.ch[ai],ai+1]        
        cr1, cr1lb, cr1ub = nxtcrbdsi(ngpl,sim.Oh.ch[ai],policy,st,ngpeitc,geitc,y,eitcDict,ngpt,ctcp)
    end                                
    
    for pt = 1:ngpt
        for l = 1:ngpl                                
            if (ai<nper)
                #println("sim.Oh.fte[ai] $(sim.Oh.fte[ai])"," l $l"," p0.αe $(p0.αe)")
                fe1val[l]= nxte(sim.Oh.fte[ai],l)
                #println("e1val $(e1val)"," ngpe $(ngpe[ai])")
                fe1lbi, fe1ubi = bounds(fe1val[l],ge[1:ngpe[ai+1],ai+1],ngpe[ai+1])                
                fe1lbval = ge[fe1lbi,ai+1]
                fe1ubval = ge[fe1ubi,ai+1]
                if l == 1
                	l1 = 1
                	resources = R*sim.Oh.k[ai] + convert(Float64,pt-1)*sim.Oh.eitc[ai] - cc[l]
                	cr1val = 0.
                	cr1lbi = 1
                	cr1ubi = 1
                	cr1lbval = 0.
                	for i = 1:ngpk
                		eulerArr[i] = (LinInterp1d(convert(Float64,fe1val[l]),convert(Float64,fe1lbval),convert(Float64,fe1ubval),
                        	lEDU[i,fe1lbi,1,l1],lEDU[i,fe1ubi,1,l1]))
                    	eulerArr[i] = eepar[s,sim.Oh.ch[ai],l,pt]*eulerArr[i] - (resources-gk[i])
                    end
                    k1[l,pt], lbk,ubk = solveEE(ngpk,eulerArr,gk,sim.Oh.k[ai],boundk[ai,s])
                    c[l,pt] = max(resources-k1[l,pt],minc[policy,sim.Oh.ch[ai]])
                    v[l,pt] = (util(c[l,pt], θ, utilpar[s,sim.Oh.ch[ai],l,pt]) 
                        + (β*calcEVinvinv(LinInterp2d(k1[l,pt],gk[lbk],gk[ubk],
                        convert(Float64,fe1val[l]),convert(Float64,fe1lbval),convert(Float64,fe1ubval),                            
                        lEV[lbk,fe1lbi,1,l1],lEV[lbk,fe1ubi,1,l1],
                        lEV[ubk,fe1lbi,1,l1],lEV[ubk,fe1ubi,1,l1]),fp)))     
                else                       
                	if l==2 || l==3
                        l1 = 2
                    else
                        l1 = 3
                    end
                    resources = R*sim.Oh.k[ai] + y[l]*hrs[l] + convert(Float64,pt-1)*sim.Oh.eitc[ai] - cc[l] - tt[l]
	                cr1val = cr1[l]
	                cr1lbi = cr1lb[l]
	                cr1ubi = cr1ub[l]
	                cr1lbval = geitc[cr1lbi] 
	                cr1ubval = geitc[cr1ubi]
		            for i = 1:ngpk
		                eulerArr[i] = (LinInterp2d(cr1val,cr1lbval,cr1ubval,
		                	convert(Float64,fe1val[l]),convert(Float64,fe1lbval),convert(Float64,fe1ubval),
		                    lEDU[i,fe1lbi,cr1lbi,l1],lEDU[i,fe1ubi,cr1lbi,l1],
		                    lEDU[i,fe1lbi,cr1ubi,l1],lEDU[i,fe1ubi,cr1ubi,l1]))
		                eulerArr[i] = eepar[s,sim.Oh.ch[ai],l,pt]*eulerArr[i] - (resources-gk[i])
		            end                
                	k1[l,pt], lbk,ubk = solveEE(ngpk,eulerArr,gk,sim.Oh.k[ai],boundk[ai,s])
                	c[l,pt] = max(resources-k1[l,pt],minc[policy,sim.Oh.ch[ai]])
					v[l,pt] = (util(c[l,pt], θ, utilpar[s,sim.Oh.ch[ai],l,pt]) 
                        + (β*calcEVinvinv(LinInterp3d(cr1val,cr1lbval,cr1ubval,k1[l,pt],gk[lbk],gk[ubk],                        
                        convert(Float64,fe1val[l]),convert(Float64,fe1lbval),convert(Float64,fe1ubval),    
                        lEV[lbk,fe1lbi,cr1lbi,l1],lEV[lbk,fe1lbi,cr1ubi,l1],
                        lEV[ubk,fe1lbi,cr1lbi,l1],lEV[ubk,fe1lbi,cr1ubi,l1],                        
                        lEV[lbk,fe1ubi,cr1lbi,l1],lEV[lbk,fe1ubi,cr1ubi,l1],
                        lEV[ubk,fe1ubi,cr1lbi,l1],lEV[ubk,fe1ubi,cr1ubi,l1]),fp)))
                end   
            else
            	if l ==1
            		resources = R*sim.Oh.k[ai] + convert(Float64,pt-1)*sim.Oh.eitc[ai] - cc[l]    
                    for i = 1:ngpk
                        eulerArr[i] = eepar[s,sim.Oh.ch[ai],1,pt]*lEDUr[i,sim.Oh.ch[ai]] - (resources-gk[i])
                    end
				else
					resources = R*sim.Oh.k[ai] +y[l]*hrs[l] + convert(Float64,pt-1)*sim.Oh.eitc[ai] - cc[l] - tt[l]
                    for i = 1:ngpk
                        eulerArr[i] = eepar[s,sim.Oh.ch[ai],l,pt]*lEDUr[i,sim.Oh.ch[ai]] - (resources-gk[i])
                    end                    
                end
                k1[l,pt], lbk,ubk = solveEE(ngpk,eulerArr,gk,sim.Oh.k[ai],boundk[ai,s])
                c[l,pt] = max(resources-k1[l,pt],minc[policy,sim.Oh.ch[ai]]) 
                v[l,pt] = (util(c[l,pt], θ, utilpar[s,sim.Oh.ch[ai],l,pt]) 
                        + (β*calcEVinvinv(LinInterp1d(k1[l,pt],gk[lbk],gk[ubk],
                        lEVr[lbk,sim.Oh.ch[ai]],lEVr[ubk,sim.Oh.ch[ai]]),fp)))                     
			end
		end
	end
    #subroutine SimValue ends here.    
    sim.Oh.V[ai], act = findmax(v)    #act[1]=optV, act[2] = optl
    #println("act $act","\t e1val(act) $(e1val[act])")
    sim.Ch.l[ai], sim.Ch.p[ai] = ind2sub(v,act)    
    if ai<fp.nper
        if sim.Ch.l[ai] == 1
            sim.Oh.l_1[ai+1] = 1
        elseif sim.Ch.l[ai] == 2 || sim.Ch.l[ai] == 3
            sim.Oh.l_1[ai+1] = 2
        else
            sim.Oh.l_1[ai+1] = 3
        end
        sim.Oh.eitc[ai+1] = cr1[sim.Ch.l[ai]]
        sim.Oh.k[ai+1] = k1[sim.Ch.l[ai], sim.Ch.p[ai]] #+ irvma*(σma0 + σma1*k1[sim.Ch.l[ai], sim.Ch.p[ai]]) 
        sim.Oh.fte[ai+1] = fe1val[sim.Ch.l[ai]]        
        if irvch<=trkids[ai,s,sim.Oh.ch[ai],1]
            sim.Oh.ch[ai+1] = 1
        elseif irvch>trkids[ai,s,sim.Oh.ch[ai],1] && irvch<=trkids[ai,s,sim.Oh.ch[ai],1] + trkids[ai,s,sim.Oh.ch[ai],2]
            sim.Oh.ch[ai+1] = 2
        else
            sim.Oh.ch[ai+1] = 3
        end        
        #println("ch_1 $(sim.Oh.ch[ai+1]), ai $(ai)")
    end    
    #println("ai $(ai)")   
    sim.Ch.c[ai] = c[act]
    sim.Ch.hrs[ai] = hrs[sim.Ch.l[ai]]
    sim.Oh.w[ai] = y[sim.Ch.l[ai]] + irvmw
    sim.Oh.inc[ai] = y[sim.Ch.l[ai]] + (sim.Oh.eitc[ai])*convert(Float64,sim.Ch.p[ai]-1)
    #println("sim.Oh.ch 2 $(sim.Oh.ch[:])"," ai $(ai)")
    return sim
end





##OLD STUFF
  #=#conditional on current state and child number    
    for ch = 1:ngpch
        #tt = ctaxsys(y,ch,sti,policy)
        if ai<nper
            ix = iEmtx.ix[ch,ai+1]
        end
        for icc = 1:ngu, ip = 1:ngu, iv = 1:ngu, ih = 1:ngu #grids for the utility shock for hours, the wage shock and the participation shock
            for pt = 1:ngpt   
                #println("pt $pt") 
                for l = 1:ngpl
                    #println("l $l")
                    if l == 1
                        l1 = 1
                    elseif l==2 || l==3
                        l1 = 2
                    else
                        l1 = 3
                    end
                    resources = R*k+y[l,iv] + cr - cc[l,icc,ch] + tt[l,iv,ch]
                    if (ai<nper)                    
                        fe1val = fe1[l,fei]
                        fe1lbi = fe1lb[l,fei]
                        fe1ubi = fe1ub[l,fei]
                        fe1lbval = ge[fe1lbi,ai+1]
                        fe1ubval = ge[fe1ubi,ai+1]                        
                        cr1val = cr1[iv,l,ch]*convert(Float64,pt-1)
                        if cr1val == 0.0
                            cr1lbi = 1
                            cr1ubi = 1
                        else
                            cr1lbi = cr1lb[iv,l,ch]
                            cr1ubi = cr1ub[iv,l,ch] 
                        end
                        cr1lbval = geitc[cr1lbi]
                        cr1ubval = geitc[cr1ubi]
                        for i =1:ngpk                               
                            eulerArr[i] = (LinInterp2d(cr1val,cr1lbval,cr1ubval,convert(Float64,fe1val),
                                convert(Float64,fe1lbval),convert(Float64,fe1ubval),
                                lEDU[cr1lbi,l1,i,fe1lbi,ix],lEDU[cr1lbi,l1,i,fe1ubi,ix],
                                lEDU[cr1ubi,l1,i,fe1lbi,ix],lEDU[cr1ubi,l1,i,fe1ubi,ix]))
                            eulerArr[i] = eepar[l,ch]*eulerArr[i] - (resources-gk[i])
                        end
                        k1, lbk, ubk = solveEE(ngpk,eulerArr,gk,ki,boundk[ai,s])
                        c[l,pt] = max(resources - k1,minc)                                
                        v[l,pt] = (util(c[l,pt], l, pt, p0, θ, utilpar[l,ch],λh[ai,s,th,ih],λp[s,ip]) 
                            + (β*calcEVinvinv(LinInterp3d(cr1val,cr1lbval,cr1ubval,k1,
                                gk[lbk],gk[ubk],convert(Float64,fe1val),convert(Float64,fe1lbval),
                                convert(Float64,fe1ubval),lEV[cr1lbi,l1,lbk,fe1lbi,ix],
                                lEV[cr1ubi,l1,lbk,fe1lbi,ix],lEV[cr1lbi,l1,ubk,fe1lbi,ix],
                                lEV[cr1ubi,l1,ubk,fe1lbi,ix],lEV[cr1lbi,l1,lbk,fe1ubi,ix],
                                lEV[cr1ubi,l1,lbk,fe1ubi,ix],lEV[cr1lbi,l1,ubk,fe1ubi,ix],
                                lEV[cr1ubi,l1,ubk,fe1ubi,ix]),fp)))
                    else
                        c[l,pt] = max(resources,minc)                                    
                        v[l,pt] = util(c[l,pt], l, pt, p0, θ, utilpar[l,ch],λh[ai,s,th,ih],λp[s,ip])          
                    end
                end
            end
            optV[ih,iv,ip,icc], act = findmax(v)
            optl = act%(size(v)[1])
            #optl,optpt = ind2sub(v,act)                                               
            optDU[ih,iv,ip,icc] = dudc(c[act], θ, dutilpar[optl,ch])                
        end   
        #println("optDU $(optDU)")
        cEV[ch] = π2*(adjw[adjw.>quantile(adjw,0.2,sorted=false)]'*optV[adjw.>quantile(adjw,0.2,sorted=false)])
        cEDU[ch] = π2*(adjw[adjw.>quantile(adjw,0.2,sorted=false)]'*optDU[adjw.>quantile(adjw,0.2,sorted=false)])
        #cEV[ch] = π2*([reshape(trpwu,ngu*ngu)'*reshape(optV[ip,:,:],ngu*ngu) for ip in 1:ngu]'*trpp)
        #cEDU[ch] = π2*([reshape(trpwu,ngu*ngu)'*reshape(optDU[ip,:,:],ngu*ngu) for ip in 1:ngu]'*trpp)
    end    
    #integrate family starts here
    for ch_1 = 1:ngpch        
        ix = iEmtx.ix[ch_1,ai]
        lEV[cri,l_1,ki,fei,ix] = calcEVinv(trkids[ai,s,ch_1,:]'*cEV,fp)
        lEDU[cri,l_1,ki,fei,ix] = calcEDUinv(trkids[ai,s,ch_1,:]'*cEDU,fp)
        #lEV[ix,ri,th,sti].m0[l_1,ki,pei,fei,cri] = calcEVinv(dot(trkids[ai,s,ch_1,:],cEV[:]),fp)
        #lEDU[ix,ri,th,sti].m0[l_1,ki,pei,fei,cri] = calcEDUinv(dot(trkids[ai,s,ch_1,:],cEDU[:]),fp)
    end
    return nothing
end




=#