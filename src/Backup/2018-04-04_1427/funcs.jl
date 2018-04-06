#function for linear interpolation
@views function LinInterp1d(x::Float64, x1::Float64, x2::Float64, fn1::Float64, fn2::Float64)
    if x1 == x2 
        return fn1
    else
        dx = (x-x1)/(x2-x1)
        return fn1+dx*(fn2-fn1)
    end
end   

#function for linear interpolation in 2 dimensions
@views function LinInterp2d(x::Float64, x1::Float64, x2::Float64, y::Float64, y1::Float64, y2::Float64, 
						fn11::Float64, fn12::Float64, fn21::Float64, fn22::Float64)
    if x1==x2
        return LinInterp1d(y,y1,y2,fn11,fn22)
    else
        if y1==y2
            return LinInterp1d(x,x1,x2,fn11,fn22)
        else
            fn1 = LinInterp1d(y,y1,y2,fn11,fn12)    #interpolate in y for the 2 points in x
            fn2 = LinInterp1d(y,y1,y2,fn21,fn22)
            return LinInterp1d(x,x1,x2,fn1,fn2)     #interpolate in x
        end
    end
end

@views function LinInterp3d(x::Float64,x1::Float64, x2::Float64, y::Float64, y1::Float64, y2::Float64, z::Float64, 
                                                z1::Float64, z2::Float64,
                                                fn111::Float64, fn211::Float64, fn121::Float64, fn221::Float64,
                                                fn112::Float64, fn212::Float64, fn122::Float64, fn222::Float64)
    if x1==x2
        return LinInterp2d(y,y1,y2,z,z1,z2,fn111,fn112,fn121,fn122)
    elseif y1==y2
        return LinInterp2d(x,x1,x2,z,z1,z2,fn111,fn112,fn211,fn212)
    elseif z1==z2
        return LinInterp2d(x,x1,x2,y,y1,y2,fn111,fn121,fn211,fn221)
    else
        fn11 = LinInterp1d(z,z1,z2,fn111,fn112) #interpolate in z for the 4 points in (x,y)        
        fn21 = LinInterp1d(z,z1,z2,fn211,fn212)
        fn12 = LinInterp1d(z,z1,z2,fn121,fn122)
        fn22 = LinInterp1d(z,z1,z2,fn221,fn222)
        fn1 = LinInterp1d(y,y1,y2,fn11,fn12) #interpolate in y for the 2 points in x
        fn2 = LinInterp1d(y,y1,y2,fn21,fn22)
        return LinInterp1d(x,x1,x2,fn1,fn2)
    end
end

function LinInterp4d{T1<:Real,T2<:Real,T3<:Real,T4<:Real}(x::T1,x1::T1, x2::T1, y::T2, y1::T2, y2::T2, z::T3, z1::T3, z2::T3,
                                                        q::T4,q1::T4,q2::T4,
                                                        fn1111::Float64, fn2111::Float64, fn1211::Float64, fn2211::Float64,
                                                        fn1121::Float64, fn2121::Float64, fn1221::Float64, fn2221::Float64,
                                                        fn1112::Float64, fn2112::Float64, fn1212::Float64, fn2212::Float64,
                                                        fn1122::Float64, fn2122::Float64, fn1222::Float64, fn2222::Float64)
    if x1==x2
        return LinInterp3d(y,y1,y2,z,z1,z2,q,q1,q2,fn1111,fn1211,fn1121,fn1221,fn1112,fn1212,fn1122,fn1222)
    elseif y1==y2
        return LinInterp3d(x,x1,x2,z,z1,z2,q,q1,q2,fn1111,fn2111,fn1121,fn2121,fn1112,fn2112,fn1122,fn2122)
    elseif z1==z2
        return LinInterp3d(x,x1,x2,y,y1,y2,q,q1,q2,fn1111,fn2111,fn1211,fn2211,fn1112,fn2112,fn1212,fn2212)
    elseif q1==q2
        return LinInterp3d(x,x1,x2,y,y1,y2,z,z1,z2,fn1111,fn2111,fn1211,fn2211,fn1121,fn2121,fn1221,fn2221)
    else
        fn111 = LinInterp1d(q,q1,q2,fn1111,fn1112) #interpolate in q for the 8 points in (x,y,z)
        fn112 = LinInterp1d(q,q1,q2,fn1121,fn1122)
        fn122 = LinInterp1d(q,q1,q2,fn1221,fn1222)
        fn121 = LinInterp1d(q,q1,q2,fn1211,fn1212)
        fn211 = LinInterp1d(q,q1,q2,fn2111,fn2112)
        fn212 = LinInterp1d(q,q1,q2,fn2121,fn2122)
        fn222 = LinInterp1d(q,q1,q2,fn2221,fn2222)
        fn221 = LinInterp1d(q,q1,q2,fn2211,fn2212)
        fn11 = LinInterp1d(z,z1,z2,fn111,fn112) #interpolate in x for the 4 points in (x,y)
        fn12 = LinInterp1d(z,z1,z2,fn121,fn122)
        fn21 = LinInterp1d(z,z1,z2,fn211,fn212)
        fn22 = LinInterp1d(z,z1,z2,fn221,fn222)
        fn1 = LinInterp1d(y,y1,y2,fn11,fn12) #interpolate in y for the 2 points in x
        fn2 = LinInterp1d(y,y1,y2,fn21,fn22)
        return LinInterp1d(x,x1,x2,fn1,fn2)
    end
end

function LinInterp5d{T1<:Real,T2<:Real,T3<:Real,T4<:Real,T5<:Real}(x::T1,x1::T1, x2::T1, y::T2, y1::T2, y2::T2, z::T3, z1::T3, z2::T3,
                                                        q::T4,q1::T4,q2::T4,r::T5,r1::T5,r2::T5,
                                                        fn11111::Float64, fn21111::Float64, fn12111::Float64, fn22111::Float64,
                                                        fn11211::Float64, fn21211::Float64, fn12211::Float64, fn22211::Float64,
                                                        fn11121::Float64, fn21121::Float64, fn12121::Float64, fn22121::Float64,
                                                        fn11221::Float64, fn21221::Float64, fn12221::Float64, fn22221::Float64,
                                                        fn11112::Float64, fn21112::Float64, fn12112::Float64, fn22112::Float64,
                                                        fn11212::Float64, fn21212::Float64, fn12212::Float64, fn22212::Float64,
                                                        fn11122::Float64, fn21122::Float64, fn12122::Float64, fn22122::Float64,
                                                        fn11222::Float64, fn21222::Float64, fn12222::Float64, fn22222::Float64)
    if x1==x2
        return LinInterp4d(y,y1,y2,z,z1,z2,q,q1,q2,r,r1,r2,fn11111,fn12111,fn11211,fn12211,fn11121,fn12121,fn11221,fn12221,
                            fn11112,fn12112,fn11212,fn12212,fn11122,fn12122,fn11222,fn12222)
    elseif y1==y2
        return LinInterp4d(x,x1,x2,z,z1,z2,q,q1,q2,r,r1,r2,fn11111,fn21111,fn11211,fn21211,fn11121,fn21121,fn11221,fn21221,
                            fn11112,fn21112,fn11212,fn21212,fn11122,fn21122,fn11222,fn21222)
    elseif z1==z2
        return LinInterp4d(x,x1,x2,y,y1,y2,q,q1,q2,r,r1,r2,fn11111,fn21111,fn12111,fn22111,fn11121,fn12121,fn12121,fn22121,
                            fn11112,fn12112,fn12112,fn22112,fn11122,fn12122,fn12122,fn22122)
    elseif q1==q2
        return LinInterp4d(x,x1,x2,y,y1,y2,z,z1,z2,r,r1,r2,fn11111,fn21111,fn12111,fn22111,fn11211,fn21211,fn12211,fn22211,
                            fn11112,fn21112,fn12112,fn22112,fn11212,fn21212,fn12212,fn22212)
    elseif r1==r2
        return LinInterp4d(x,x1,x2,y,y1,y2,z,z1,z2,q,q1,q2,fn11111,fn21111,fn12111,fn22111,fn11211,fn21211,fn12211,fn22211,
                            fn11121,fn21121,fn12121,fn22121,fn11221,fn21221,fn12221,fn22221)
    else

        fn1111 = LinInterp3d(r,r1,r2,fn11111,fn11112) #interpolate in r for th 16 points in (x,y,z,q)
        fn1112 = LinInterp3d(r,r1,r2,fn11121,fn11122)
        fn1121 = LinInterp3d(r,r1,r2,fn11211,fn11212)
        fn1122 = LinInterp3d(r,r1,r2,fn11221,fn11222)
        fn1211 = LinInterp3d(r,r1,r2,fn12111,fn12112)
        fn1212 = LinInterp3d(r,r1,r2,fn12121,fn12122)
        fn1221 = LinInterp3d(r,r1,r2,fn12211,fn12212)
        fn1222 = LinInterp3d(r,r1,r2,fn12221,fn12222)
        fn2111 = LinInterp3d(r,r1,r2,fn21111,fn21112)
        fn2112 = LinInterp3d(r,r1,r2,fn21121,fn21122)
        fn2121 = LinInterp3d(r,r1,r2,fn21211,fn21212)
        fn2122 = LinInterp3d(r,r1,r2,fn21221,fn21222)
        fn2211 = LinInterp3d(r,r1,r2,fn22111,fn22112)
        fn2212 = LinInterp3d(r,r1,r2,fn22121,fn22122)
        fn2221 = LinInterp3d(r,r1,r2,fn22211,fn22212)
        fn2222 = LinInterp3d(r,r1,r2,fn22221,fn22222)
        fn111 = LinInterp2d(q,q1,q2,fn1111,fn1112) #interpolate in q for the 8 points in (x,y,z)
        fn112 = LinInterp2d(q,q1,q2,fn1121,fn1122)
        fn122 = LinInterp2d(q,q1,q2,fn1221,fn1222)
        fn121 = LinInterp2d(q,q1,q2,fn1211,fn1212)
        fn211 = LinInterp2d(q,q1,q2,fn2111,fn2112)
        fn212 = LinInterp2d(q,q1,q2,fn2121,fn2122)
        fn222 = LinInterp2d(q,q1,q2,fn2221,fn2222)
        fn221 = LinInterp2d(q,q1,q2,fn2211,fn2212)
        fn11 = LinInterp2d(z,z1,z2,fn111,fn112) #interpolate in x for the 4 points in (x,y)
        fn12 = LinInterp2d(z,z1,z2,fn121,fn122)
        fn21 = LinInterp2d(z,z1,z2,fn211,fn212)
        fn22 = LinInterp2d(z,z1,z2,fn221,fn222)
        fn1 = LinInterp2d(y,y1,y2,fn11,fn12) #interpolate in y for the 2 points in x
        fn2 = LinInterp2d(y,y1,y2,fn21,fn22)
        return LinInterp1d(x,x1,x2,fn1,fn2)
    end
end




function lambdah(p::params, l::Int64)
    @unpack αc1, αc2, αc3, αc4 = p
    if l == 1
        res = 0
    elseif l == 2
        res = αc1
    elseif l == 3
        res = αc1 + αc2
    elseif l == 4
        res = αc1 + αc2 + αc3
    else
        res = αc1 + αc2 + αc3 + αc4
    end
    return res
end


function fxh(p::params, s::Int64, ch::Int64)
    @unpack αh0, αh12, αh13, αh14, αh21, αh22 = p
    res = αh0
    if s == 2
        res += αh12
    elseif s == 3
        res += αh13
    elseif s == 4
        res += αh14
    end
    #ch = 1 means no kids
    if ch == 2 #refers to having one kid
        res += αh21
    elseif ch == 3 #refers to having two kids
        res += αh22
    end
    return res
end

function fxp(p::params, s::Int64)
    @unpack αp11, αp12, αp13, αp14 = p
    if s == 1
        return αp11
    elseif s == 2
        return αp12
    elseif s == 3
        return αp13
    elseif s == 4
        return αp14
    end
end


#utility shifter for the marginal disutility of labor
function lambdah(p::params, age::Int64, educ::Int64, ty::Int64, epsh::Float64)
    @unpack αh01, αh02, αh1, αh2, αh3, αh42, αh43, αh44 = p
    if ty == 1
        p1 = αh01
    else
        p1 = αh02
    end
    if age <= 30
        p1 += αh1
    end
    if age >= 40
        p1 += αh2
    end
    if age >= 50
        p1 += αh3
    end
    if educ == 1
        return (p1 + epsh)
    elseif educ == 2
        return (p1 + αh42 + epsh)
    elseif educ == 3
        return (p1 + αh43 + epsh)
    else
        return (p1 + αh44 + epsh)        
    end
end

function lambdap(p::params, educ::Int64,epsp::Float64)
    @unpack αp0, αp12, αp13, αp14 = p
    if educ == 1
        return (αp0 + epsp)
    elseif educ == 2
        return (αp0 + αp12 + epsp)
    elseif educ == 3
        return (αp0 + αp13 + epsp)
    else
        return (αp0 + αp14 + epsp)
    end
end




#the utility function
function util(c::Float64, θ::Float64, utilpar::Float64)    
    #pt = 1 not participating, 2 participating
    #l = 1 - 0, 2-755, 3-1520, 4-2020, 5-2720
    #println("utilpar $(utilpar)","c $(c)","θ $(θ)")
    p1 = utilpar*((c)^θ)        
    return p1
end


#the derivative of the utility function wrt consumption
function dudc(c::Float64, θ::Float64, dutilpar::Float64)    
    if isnan(c^(θ-1))
        println("c $(c)")
    end
    return ((c)^(θ-1))*dutilpar
end


#wage equation
function wage(fe::Int64, l::Int64, educ::Int64, l1::Int64, p::params, gεwiu::Float64)
    @unpack αw0, αw11, αw12, αw13, αw32, αw33, αw34, αw51, αw52, αw53, αw54, αw61, αw62 = p
    
    if l == 1
    	return 0.
    else
	    #=if th == 1
	    	αw0 = αw01
	    else
	    	αw0 = αw02
	    end=#
	    if l == 2
	        p1 = αw0
	    elseif l == 3
	        p1 = αw0 + αw11
	    elseif l == 4
	        p1 = αw0 + αw12
	    elseif l == 5
	        p1 = αw0 + αw13
	    end
	    #=if r == 1
	        p1 += αw2
	    end=#
	    if educ == 1
	        p1 += log(fe+1.0)*αw51 #log(pe+1.0)*αw41 +
	    elseif educ == 2
	        p1 += αw32 + log(fe+1.0)*αw52 #+ log(pe+1.0)*αw42
	    elseif educ == 3
	        p1 += αw33 + log(fe+1.0)*αw53 #+ log(pe+1.0)*αw43 
	    elseif educ == 4
	        p1 += αw34 + log(fe+1.0)*αw54 #+ log(pe+1.0)*αw44
	    end
	    if l1 == 2
	        p1 += αw61
	    elseif l1 == 3
	        p1 += αw62
	    end
    	return exp(p1 + gεwiu)
    end
end

function wage(fe::Int64, l::Array{Int64,1}, educ::Int64, l1::Int64, p::params, gεwiu::Array{Float64,1})
	w = zeros(Float64,size(l,1),size(gεwiu,1))
	for iu = 1:size(gεwiu,2), i = 1:size(l,2)
		w[i,iu] = wage(fe,l[i],educ,l1,p,gεwiu[iu])
	end	
	return w
end

function wage(fe::Int64, l::Array{Int64,1}, educ::Int64, l1::Int64, p::params, gεwiu::Float64)
    w = zeros(Float64,size(l,1))
    for i = 1:length(l)
        w[i] = wage(fe,l[i],educ,l1,p,gεwiu)
    end 
    return w
end


function nxte(te::Int64,l::Int)
    if l == 0  
        fe1 = te
        #pe1 = pe
    else
        if l < 4
            fe1 = te
            #pe1 = pe + 1
        else
            fe1 = te + 1
            #pe1 = pe
        end
    end
    #println("e1 $(e1)")
    return fe1
end


#approximate inverse of EdU 
@views function calcEDUinv(Edu::Float64, fp::fparams)
    @unpack θ = fp
    return Edu^(1./(θ-1.))
end

#approximate inverse of EV
@views function calcEVinv(EV::Float64, fp::fparams)
    @unpack θ = fp
    return (θ*EV)^(1./θ)    
end

#reversing the inverse of EV
@views function calcEVinvinv(EV::Float64, fp::fparams)
    @unpack θ = fp
    return (1./θ)*(EV^θ)    
end

function boundk(p0::params, age::Int, educ::Int)
    @unpack αb0, αb1, αb2, αb32, αb33, αb34 = p0     
    agef = convert(Float64,age)    
    if educ == 1
        mink = -exp(αb0 + αb1*agef/10 + αb2*(agef^2)/100)
    elseif educ == 2
        mink = -exp(αb0 + αb1*agef/10 + αb2*(agef^2)/100 + αb32)
    elseif educ == 3
        mink = -exp(αb0 + αb1*agef/10 + αb2*(agef^2)/100 + αb33)
    else
        mink = -exp(αb0 + αb1*agef/10 + αb2*(agef^2)/100 + αb34)
    end
    return mink
end

@views function getFedCredit(eitcDict::datadict, earnings::Float64, policy::Int64, numkids::Int64)
    @unpack fedeitc = eitcDict
    ix = eitcDict.ix[policy,numkids]    
    
    fedCredit = 0.
    if earnings < fedeitc[ix,2]/1000.
        fedCredit = fedeitc[ix,1]/100*earnings   
    elseif earnings >= fedeitc[ix,2]/1000. && earnings <= fedeitc[ix,5]/1000.
        fedCredit = fedeitc[ix,3]/1000.
    elseif earnings > fedeitc[ix,5]/1000. && earnings < fedeitc[ix,6]/1000.
        fedCredit = fedeitc[ix,3]/1000. - (fedeitc[ix,4]/100)*(earnings - fedeitc[ix,5]/1000.)
    end
    return max(fedCredit,0)
end

# calculate amount of state tax credit
# given the state, federal credit amount, policy, and the number of kids
@views function getStateCredit(eitcDict::datadict, state::Int64, fedCredit::Float64, policy::Int64, numkids::Int64)
    # find a matching key for the closest number of kids
    # Ex. DC in 2000 only has state credit %s for 0 children
    # so if you have 3 children, the credit % for
    # 0 children is returned
    ix = eitcDict.ix[policy,numkids]
    creditPct = eitcDict.steitc[ix,state]
    return creditPct/100*fedCredit
end



# calculate amount of tax credit given the parameters 
# earnings, state, policy, and # of children
@views function getCredit(eitcDict::datadict, earnings::Float64, state::Int64, policy::Int64, numkids::Int64)
    fedCredit = getFedCredit(eitcDict, earnings, policy, numkids)
    stateCredit = getStateCredit(eitcDict, state, fedCredit, policy, numkids)
    return fedCredit + stateCredit
end



function inittrkids(fp::fparams)
    @unpack nper, nts, ngpch, mina = fp
    trmkids = Array{Float64}(nper,nts,ngpch,ngpch)
    for ai = 1:nper
        for si = 1:nts
            age = ai + mina - 1
            if age<=25
                trmkids[ai,si,:,:] = readdlm("./trans_matrices/trkids_a25_e$(si).txt")
            elseif age>25 && age<=30
                	trmkids[ai,si,:,:] = readdlm("./trans_matrices/trkids_a2530_e$(si).txt")
            elseif age>30 && age<=40
                trmkids[ai,si,:,:] = readdlm("./trans_matrices/trkids_a3040_e$(si).txt")
            elseif age>40 && age<=50
                trmkids[ai,si,:,:] = readdlm("./trans_matrices/trkids_a4050_e$(si).txt")
            else
                trmkids[ai,si,:,:] = readdlm("./trans_matrices/trkids_a50_e$(si).txt")
            end
        end
    end
    return trmkids
end


function ccare(l::Int, ch::Int, s::Int, ccp::Array{Float64,3})	
	if s <3		
		si = s
	else
		si = 3
	end
	if l == 1
		ccarem = 0.
	elseif l == 2 || l == 3
		ccarem = ccp[1,ch,si]
	else
		ccarem = ccp[2,ch,si]
	end
	return ccarem
end

function ccare(ch::Int, s::Int, ccp::Array{Float64,3})	
	ccarem = Array{Float64}(5)
	for l = 1:5
		ccarem[l] = ccare(l,ch-1,s,ccp)
	end
	return ccarem
end

function ccare(s::Int, ccp::Array{Float64,3})	
	ccarem = Array{Float64}(5,3)
	for ch = 1:3		
		if ch == 1
			ccarem[:,1] .= 0.
		else
			for l = 1:5
				ccarem[l,ch] = ccare(l,ch-1,s,ccp)
			end
		end
	end
	return ccarem
end


@views function quad(nds::Array{Float64,1},p0)
    @unpack σw = p0
    #Σ = [σh^2 ρhw*σh*σw 0. 0.; ρhw*σh*σw σw^2 0. 0.; 0. 0. σp^2 0.; 0. 0. 0. 4.11^2]
    #ln = length(nds)
    #adjn = zeros(Float64,ln)
    #C = cholfact(2.*Σ)    
    C = σw*(2^0.5)
    adjn = C*nds    
    return adjn
end

#=function ftax(y::Float64, ch::Int, policy::Int,taxDict::datadict2)
    br = ch == 1 ? taxDict.fedtax[1:16,policy] : taxDict.fedtax[1:16,4+policy]
    tr = ch == 1 ? taxDict.fedtax[16:32,policy] : taxDict.fedtax[16:32,4+policy]
    ty = ch == 1 ? y - taxDict.fedex[policy] - taxDict.fedded[policy,1] : y - taxDict.fedex[policy]*ch - taxDict.fedded[policy,2]  
    ft = 0.0
    for i = 2:length(br)
        ft += (ty>br[i])*(tr[i]*(br[i]-br[i-1])) + (ty<=br[i])*((ty-br[i-1])*tr[i])*(ty>br[i-1])
    end
    return ft
end=#

#=function ftax(y::Float64, ch::Int64, policy::Int64, taxDict::datadict2)
    if ch == 1
        s = 1
    else
        s = 2
    end
    ix = taxDict.ix[policy,s]
    #println("ix $(ix)")
    b = taxDict.b[ix:ix+4]
    t = taxDict.t[ix:ix+4]
    ft = 0.
    for bi = 2:16
        if y <= b[bi] && y >b[bi-1]
            ft += (y-b[bi-1])*t[bi-1]
            break
        else
            ft += (b[bi] - b[bi-1])*t[bi-1]
            break
        end
    end
    return ft
end=#

function stax(y::Float64,policy::Int64,st::Int64,stDict::datadict4)
    ix = stDict.ix[policy,st]
    stTax = stDict.b0[ix] + stDict.b1[ix]*y + stDict.b2[ix]*(y^2)
    return stTax
end


function ftax(y::Float64,policy::Int64,ftDict::datadict5)    
    ft = ftDict.b0[policy] + ftDict.b1[policy]*y + ftDict.b2[policy]*(y^2)
    return ft
end


function getCtc(y::Float64, ch::Int64,  policy::Int64, ctcDict::datadict3)
    if policy < 3
        credit = 0.
    else
        ix = ctcDict.ix[policy-2]
        phthr = ctcDict.phthr[ix]
        maxcr = ctcDict.maxcr[ix]
        phrate = ctcDict.phrate[ix]
        if y < phthr
            credit = maxcr*convert(Float64,ch-1)
        else 
            credit = maxcr*convert(Float64,ch-1) - phrate*(y - phthr)
        end
    end
    return credit
end


#=function ctc(y::Float64, ch::Int, policy::Int, ftax::Float64,ctc::Array{Float64,2})    
    if policy < 3
        ctcref = 0.
        ctcnonref = 0.
    else
        credit = (y < ctc[policy-2,5] ? ctc[policy-2,1]*ch : ctc[policy-2,1]*ch - ctc[policy-2,7]*(earnings - ctc[policy-2,5]))
        ctcnonref = min(credit, ftax)
        ctcref = (y < ctc[policy-2,2] ? 0. : min(ctc[policy-2,3]*(earnings - ctc[policy-2,2]),ctc[policy-2,1]*ch))
    end
    return (ctcref, ctcnonref)
end=#

function ctaxsys(y::Float64,st::Int,policy::Int,ftaxp::datadict5, staxp::datadict4)
    ft = ftax(y, policy, ftaxp)    
    #ctcref = getCtc(y,ch,policy,ctcp)
    st = stax(y,policy,st,staxp)
    #ft = max(ft - ctcnonref,0.)
    return ft + st
end

function ctaxsys(y::Array{Float64,2},st::Int,policy::Int,ftaxp::datadict5,staxp::datadict4)
    tt = zeros(Float64,size(y,1),size(y,2))    
    for i in 1:size(y,1), j in 1:size(y,2)
        tt[i,j] = ctaxsys(y[i,j],st,policy,ftaxp,staxp)
    end
    return tt
end

function ctaxsys(y::Array{Float64,1},st::Int,policy::Int,ftaxp::datadict5,staxp::datadict4)
    tt = zeros(Float64,length(y))    
    for i in 1:length(y)
        tt[i] = ctaxsys(y[i],st,policy,ftaxp,staxp)
    end
    return tt
end

#=function getParams(ch::Int64, index::Int64,ftaxp::ftaxparam)
    ex = ftaxp.fedex[index] 
    if ch == 1        
        b1 = ftaxp.b1s[index]; b2 = ftaxp.b2s[index]; b3 = ftaxp.b3s[index]; b4 = ftaxp.b4s[index]; b5 = ftaxp.b5s[index];
        b6 = ftaxp.b6s[index]; b7 = ftaxp.b7s[index]; b8 = ftaxp.b8s[index]; b9 = ftaxp.b9s[index]; b10 = ftaxp.b10s[index];
        b11 = ftaxp.b11s[index]; b12 = ftaxp.b12s[index]; b13 = ftaxp.b13s[index]; b14 = ftaxp.b14s[index]; 
        b15 = ftaxp.b15s[index]; b16 = ftaxp.b16s[index];
        t1 = ftaxp.t1s[index]; t2 = ftaxp.t2s[index]; t3 = ftaxp.t3s[index]; t4 = ftaxp.t4s[index]; t5 = ftaxp.t5s[index];
        t6 = ftaxp.t6s[index]; t7 = ftaxp.t7s[index]; t8 = ftaxp.t8s[index]; t9 = ftaxp.t9s[index]; t10 = ftaxp.t10s[index];
        t11 = ftaxp.t11s[index]; t12 = ftaxp.t12s[index]; t13 = ftaxp.t13s[index]; t14 = ftaxp.t14s[index]; 
        t15 = ftaxp.t15s[index]; t16 = ftaxp.t16s[index];
    else        
        b1 = ftaxp.b1h[index]; b2 = ftaxp.b2h[index]; b3 = ftaxp.b3h[index]; b4 = ftaxp.b4h[index]; b5 = ftaxp.b5h[index];
        b6 = ftaxp.b6h[index]; b7 = ftaxp.b7h[index]; b8 = ftaxp.b8h[index]; b9 = ftaxp.b9h[index]; b10 = ftaxp.b10h[index];
        b11 = ftaxp.b11h[index]; b12 = ftaxp.b12h[index]; b13 = ftaxp.b13h[index]; b14 = ftaxp.b14h[index]; 
        b15 = ftaxp.b15h[index]; b16 = ftaxp.b16h[index];
        t1 = ftaxp.t1h[index]; t2 = ftaxp.t2h[index]; t3 = ftaxp.t3h[index]; t4 = ftaxp.t4h[index]; t5 = ftaxp.t5h[index];
        t6 = ftaxp.t6h[index]; t7 = ftaxp.t7h[index]; t8 = ftaxp.t8h[index]; t9 = ftaxp.t9h[index]; t10 = ftaxp.t10h[index];
        t11 = ftaxp.t11h[index]; t12 = ftaxp.t12h[index]; t13 = ftaxp.t13h[index]; t14 = ftaxp.t14h[index]; 
        t15 = ftaxp.t15h[index]; t16 = ftaxp.t16h[index];
    end        
    return (ex, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16,
        t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16)
end =#


#=@views function nxtcr(ngu::Int64, ngpl::Int64, ngpch::Int, policy::Int, sti::Int, y::Array{Float64,2}, eitcDict::datadict)    
    cr1 = Array{Float64}(ngu,ngpl,ngpch) 
    
    @inbounds for chi = 1:ngpch
        for l = 1:ngpl
            for iv = 1:ngu                                
            	cr1[iv,l,chi] = getCredit(eitcDict, y[l,iv], sti, policy, chi)                    
            end
        end    
    end
    return cr1
end =#
