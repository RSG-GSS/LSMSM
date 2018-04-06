#constructing the grids for the random shocks
#=function popgerr(nguj::Int64, ngu::Int64, p0::params)
    gεw = Array{Float64}(nguj)
    gεu = Array{Float64}(nguj)
    gεp = Array{Float64}(ngu)    
    @unpack σw, σh, σp = p0

    @inbounds for i = 1:nguj
        p = (i.-0.5)/nguj
        gεw[i] = norminvcdf(0,σw,p)
        gεu[i] = norminvcdf(0,σh,p)        
    end
    @inbounds for i = 1:ngu
        p = (i. -0.5)/ngu
        gεp[i] = norminvcdf(0,σp,p)        
    return gεw, gεu, gεp
end

function popgerrcc(ngu::Int64)    
    gεcc = Array{Float64}(ngu)
    @inbounds for i = 1:ngu
        p = (i.-0.5)/ngu
        gεcc[i] = norminvcdf(0,4.11,p)        
    end
    return gεcc
end=#

#constructing the grids for the full-time and part-time experience variables
function popge(fp::fparams)
    @unpack nper, maxngpe, mina, nts, as, maxa = fp
    auxge   = [0,6,16,44]
    maxe = Array{Int64}(nper)
    for a = 1:nper
        @inbounds maxe[a] = a - (mina - 15)
    end
    ngpe = zeros(Int64,nper)
    ge  = zeros(Int64,maxngpe,nper)
    for a = 1:nper 
        for ei = 2:maxngpe
            if maxe[a]<=auxge[ei]
                @inbounds ge[ei,a] = maxe[a]
                @inbounds ngpe[a] = ei
                break
            else
                @inbounds ge[ei,a] = auxge[ei]
            end
            if ngpe[a] == 0
                @inbounds ngpe[a] = maxngpe
            end                 
        end
    end
    return ngpe, ge
end    



function popgeitc(fp::fparams)
	@unpack maxeitc, ngpeitc, ratioe = fp
	geitc = Array{Float64}(ngpeitc)
	geitc[1] = 0.
    unit = (maxeitc-0)*(1.-ratioe)/(1.-ratioe^(ngpeitc - 1.))
	for ei = 2:ngpeitc		
		@inbounds geitc[ei] = geitc[ei-1]+unit*ratioe^(ei-2)
	end
	return geitc
end

#constructing the grids for the assets
function popgk(fp::fparams)
    @unpack mink,maxk,ngpk,ratiok = fp    
    ngpkf   = convert(Float64,ngpk)
    gk      = Array{Float64}(ngpk)
    gk[1]   = mink
    unit    = (maxk - mink)*(1.-ratiok)/(1.-ratiok^(ngpkf-1.))
    for ki = 2:ngpk
        @inbounds gk[ki] = gk[ki-1] + unit*ratiok^(ki-2)
    end
    return gk
end


#constructing the bounds for next year's experience based on this period's labor supply choice
function nxtebds(ngpl::Int64, maxngpe::Int64, ngpe::Array{Int64,1}, ge::Array{Int64,2}, nper::Int64)
    fe1lb = zeros(Int64,ngpl,maxngpe,nper)
    fe1ub = zeros(Int64,ngpl,maxngpe,nper)
    fe1 = zeros(Int64,ngpl,maxngpe,nper)    
    #pe1lb = zeros(Int64,ngpl,maxngpe,nper)
    #pe1ub = zeros(Int64,ngpl,maxngpe,nper)
    #pe1 = zeros(Int64,ngpl,maxngpe,nper)        
    @inbounds for ai = 1:nper-1
        a1 = ai+1
        for ei = 1:ngpe[ai] 
            for l = 1:ngpl
                fe1[l,ei,ai] = nxte(ge[ei,ai],l)
                #println("l $l","ei $ei","e1[l,ei] $(e1[l,ei])")
                for nei = 1:ngpe[a1]
                    if abs(ge[nei,a1]-fe1[l,ei,ai]) <= 0.00001
                        fe1lb[l,ei,ai] = nei
                        fe1ub[l,ei,ai] = nei
                        break
                    elseif ge[nei,a1]>fe1[l,ei,ai]
                        fe1lb[l,ei,ai] = nei-1
                        fe1ub[l,ei,ai] = nei
                        break
                    end
                end
                #=for nei = 1:ngpe[a1]
                    if abs(ge[nei,a1]-pe1[l,ei,ai]) <= 0.00001 
                        pe1lb[l,ei,ai] = nei
                        pe1ub[l,ei,ai] = nei
                        break
                    elseif ge[nei,a1]>pe1[l,ei,ai]
                        pe1lb[l,ei,ai] = nei-1
                        pe1ub[l,ei,ai] = nei
                        break
                    end
                end=#
            end
        end
    end
    return fe1, fe1lb, fe1ub#, pe1, pe1lb, pe1ub
end



@views function nxtcrbds(ngu::Int64, ngpl::Int64, ngpch::Int, policy::Int, sti::Int, ngpeitc::Int, 
                        geitc::Array{Float64,1}, y::Array{Float64,2}, eitcDict::datadict,ctcDict::datadict3)
    cr1lb = Array{Int64}(ngu,ngpl,ngpch)
    cr1ub = Array{Int64}(ngu,ngpl,ngpch)
    cr1 = Array{Float64}(ngu,ngpl,ngpch) 
    @inbounds for chi = 1:ngpch
        for l = 1:ngpl
            for iv = 1:ngu                                
            cr1[iv,l,chi] = getCredit(eitcDict, y[l,iv], sti, policy, chi) + getCtc(y[l,iv],chi,policy,ctcDict)
                for nei = 1:ngpeitc
                    if geitc[nei] - cr1[iv,l,chi] <= 0.00001 && geitc[nei] - cr1[iv,l,chi] >= -0.00001
                        cr1lb[iv,l,chi] = nei
                        cr1ub[iv,l,chi] = nei
                        break
                    elseif geitc[nei]>cr1[iv,l,chi]
                        cr1lb[iv,l,chi] = nei-1
                        cr1ub[iv,l,chi] = nei
                        break
                    end
                end
            end
        end    
    end
    return cr1, cr1lb, cr1ub
end

@views function nxtcrbdsi(ngpl::Int64, chi::Int, policy::Int, sti::Int, ngpeitc::Int, 
                        geitc::Array{Float64,1}, y::Array{Float64,1}, eitcDict::datadict, ngpt::Int,ctcDict::datadict3)
    cr1lb = Array{Int64}(ngpl)
    cr1ub = Array{Int64}(ngpl)
    cr1 = Array{Float64}(ngpl)     
    for l = 1:ngpl                                         
        cr1[l] = getCredit(eitcDict, y[l], sti, policy, chi) + getCtc(y[l],chi,policy,ctcDict)
        for nei = 1:ngpeitc
            if geitc[nei] - cr1[l] <= 0.00001 && geitc[nei] - cr1[l] >= -0.00001
                cr1lb[l] = nei
                cr1ub[l] = nei
                break
            elseif geitc[nei]>cr1[l]
                cr1lb[l] = nei-1
                cr1ub[l] = nei
                break
            end
        end
    end                
    return cr1, cr1lb, cr1ub
end



#this function determines the interval of a value of a variable belongs to given the grid strcuture
@views function bounds(x::Float64,gx::Array{Float64,1},ngpx::Int)  
    #println("gx $gx")
    lb = 0
    ub = 0    
    i = 1  
    while true
        #println("i $i","\t lb $lb", "\t ub $ub")
        #println("x $x","\t gx[i] $(gx)")
        if x-gx[i]<0.00001 && x-gx[i]>-0.00001
            #println("here1")
            lb = i
            ub = i
            break
        elseif x<gx[i] || i==ngpx
            #println("here2")
            lb = max(i-1,1)
            ub = lb+1
            break
        end
        if i>=ngpx
            break
        end
        i=i+1
    end
    #println("lb $lb", "ub $ub")
    return lb, ub
end

@views function bounds(x::Int,gx::Array{Int,1},ngpx::Int)  
    #println("gx $gx")
    lb = 0
    ub = 0    
    i = 1  
    while true
        #println("i $i","\t lb $lb", "\t ub $ub")
        #println("x $x","\t gx[i] $(gx)")
        if convert(Float64,x-gx[i])<0.00001 && convert(Float64,x-gx[i])>-0.00001 
            #println("here1")
            lb = i
            ub = i
            break
        elseif x<gx[i] || i==ngpx
            #println("here2")
            lb = max(i-1,1)
            ub = lb+1
            break
        end
        if i>=ngpx
            break
        end
        i=i+1
    end
    #println("lb $lb", "ub $ub")
    return lb, ub
end

#this function locates the bounds and thus the index of a given value of a variable (for random shocks)
function ilocate(i::Int,nper::Int)        
    lb = (i-1)*(nper)+1
    ub = i*(nper)     
    return lb,ub
end


#calculating the probabilities for the grids of the random shocks
#=@views function trprob(gεw::Array{Float64,1}, gεu::Array{Float64,1}, gεp::Array{Float64,1}, gεcc::Array{Float64,1}, p0::params, fp::fparams)
    @unpack ngu,nguj = fp
    @unpack σw, σh, σp = p0
    derrw = Normal(0,σw)
    derru = Normal(0,σh)
    derrp = Normal(0,σp)
    derrcc = Normal(0,4.11)
    cdfw = Array{Float64}(ngu)
    cdfu = Array{Float64}(ngu)
    cdfp = Array{Float64}(ngu)
    cdfwu = Array{Float64}(ngu,ngu)
    pgridwu = Array{Float64}(ngu,ngu)
    pgridp = Array{Float64}(ngu)
    for i = 1:nguj
        cdfw[i] = cdf(derrw,gεw[i])
        cdfu[i] = cdf(derru,gεu[i])        
    end
    for i = 1:ngu
        cdfp[i] = cdf(derrp,gεp[i])
        cdfcc[i] = cdf(derrcc,gεcc[i])
    end

    @inbounds for iu = 1:nguj
        for iv = 1:nguj
            @inbounds cdfwu[iv,iu] = cdfw[iv]*cdfu[iu]
        end
    end

    @inbounds for iu = 1:ngu        
        for iv = 1:ngu            
            #println("iv $iv","iu $iu")
            if iv == 1
                if iu == 1
                    #println("here1")
                    pgridwu[iv,iu] = ((cdfwu[iv,iu] + cdfwu[iv+1,iu])/2. + (cdfwu[iv,iu+1] + cdfwu[iv+1,iu+1])/2.)/2.
                else
                    if iu < ngu
                        #println("here4")
                        pgridwu[iv,iu] = ((cdfwu[iv,iu+1] - cdfwu[iv,iu-1])/2. + (cdfwu[iv+1,iu+1] - cdfwu[iv+1,iu-1])/2.)/2.
                    else
                        #println("here5")
                        pgridwu[iv,iu] = (cdfw[iv]+cdfw[iv+1])/2. - ((cdfwu[iv,iu-1] + cdfwu[iv,iu])/2. + (cdfwu[iv+1,iu-1] + cdfwu[iv+1,iu])/2.)/2.
                    end
                end
            else
                if iv <ngu
                    if iu == 1
                        #println("here2")
                        pgridwu[iv,iu] = ((cdfwu[iv+1,iu] - cdfwu[iv-1,iu])/2. + (cdfwu[iv+1,iu+1] - cdfwu[iv-1,iu+1])/2.)/2.
                    else
                        if iu < ngu
                            #println("here8")    
                            pgridwu[iv,iu] = ((cdfwu[iv+1,iu+1] - cdfwu[iv-1,iu+1])/2. - (cdfwu[iv+1,iu-1] - cdfwu[iv-1,iu-1])/2.)/2.
                        else
                            #println("here6")
                            pgridwu[iv,iu] = (cdfw[iv+1]-cdfw[iv-1])/2. - ((cdfwu[iv+1,iu-1] + cdfwu[iv+1,iu])/2. - (cdfwu[iv-1,iu-1] + cdfwu[iv-1,iu])/2.)/2.
                        end
                    end
                else
                    if iu == 1
                        #println("here3")
                         pgridwu[iv,iu] = (cdfu[iu]+cdfu[iu+1])/2. - ((cdfwu[iv-1,iu] + cdfwu[iv,iu])/2. + (cdfwu[iv-1,iu+1] + cdfwu[iv,iu+1])/2.)/2.
                    else
                        if iu < ngu
                            #println("here7")
                            pgridwu[iv,iu] = (cdfu[iu+1]-cdfu[iu-1])/2. - ((cdfwu[iv-1,iu+1] + cdfwu[iv,iu+1])/2. - (cdfwu[iv-1,iu-1] + cdfwu[iv,iu-1])/2.)/2.
                        else
                            #println("here9")
                            pgridwu[iv,iu] = (1 - (cdfw[iv-1] + cdfw[iv])/2. - (cdfu[iu-1] + cdfu[iu])/2.
                                        + ((cdfwu[iv-1,iu-1] + cdfwu[iv,iu-1])/2. + (cdfwu[iv-1,iu] + cdfwu[iv,iu])/2.)/2.)
                        end
                    end
                end
            end
        end
    end
    pgridp[1] = (cdfp[2]+cdfp[1])/2.
    @inbounds for ip=2:ngu-1
        pgridp[ip] = (cdfp[ip+1]-cdfp[ip-1])/2.
    end
    pgridp[ngu] = 1 - (cdfp[ngu-1]+cdfp[ngu])/2.
    return pgridwu, pgridp
end=#

function pop_Emtx(fp::fparams)
	@unpack toti, mina, nper, as, maxa, ngpch = fp
	iEmtx = initindex(fp)	
	ix = 0
	@inbounds for ai = 1:maxa-mina+1
		for ki = 1:ngpch				
			ix = ix +1
			#iEmtx.ix[ki,ai,si] = ix					
            iEmtx.ix[ki,ai] = ix                 
			iEmtx.ch[ix] = ki
			iEmtx.a[ix] = ai+mina-1
			#iEmtx.s[ix] = si			
		end
	end
	return iEmtx
end

function pop_sEmtx(fp::fparams)
    @unpack nind, nsim = fp
    isEmtx = initsimindex(fp)
    ix = 0
    @inbounds for i = 1:nind
        for r = 1:nsim
            ix = ix+1
            isEmtx.ix[r,i] = ix
            isEmtx.rep[ix] = r
            isEmtx.ind[ix] = i
        end
    end
    return isEmtx
end

function pop_cEmtx(fp::fparams)
    @unpack totcti, nstate, nts= fp
    icEmtx = initctindex(fp)   
    ix = 0
    @inbounds for s =1:nts, sti = 1:nstate #, l1 = 1:ngpl1
    	ix = ix + 1
    	icEmtx.ix[sti,s] = ix
    	icEmtx.st[ix] = sti
    	icEmtx.s[ix] = s
    end
    return icEmtx
end







        #=if l1 == 1
            cri = 1
            ix = ix + 1
            icEmtx.ix[sti,s] = ix
            #icEmtx.ty[ix] = th
            icEmtx.st[ix] = sti
            #icEmtx.r[ix] = ri
            icEmtx.s[ix] = s
            icEmtx.l1[ix] = l1
            icEmtx.cr[ix] = cri
        else
            for cri = 1:ngpeitc
                ix = ix + 1
                icEmtx.ix[sti,s] = ix
                #icEmtx.ty[ix] = th
                icEmtx.st[ix] = sti
                #icEmtx.r[ix] = ri
                icEmtx.s[ix] = s
                icEmtx.l1[ix] = l1
                icEmtx.cr[ix] = cri
            end
        end =#

