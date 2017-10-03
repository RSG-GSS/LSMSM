#constructing the grids for the random shocks
function popgerr(fp::fparams, p0::params)
    @unpack ngu = fp
    gεw = Array{Float64}(ngu)
    gεu = Array{Float64}(ngu)
    @unpack σw, σu = p0

    for i = 1:ngu
        p = (i.-0.5)/ngu
        gεw[i] = norminvcdf(0,σw,p)
        gεu[i] = norminvcdf(0,σu,p)
    end
    return gεw, gεu
end

#constructing the grids for the experience variable
function popge(fp::fparams)
    @unpack nper, maxngpe, mina = fp
    auxge   = [0;2;6;16;40]
    maxe = Array{Int64}(nper)
    ngpe = Array{Int64}(nper)
    ge  = zeros(Int64,maxngpe, nper)
    for a = 1:nper
        maxe[a] = a - (mina-19)     #use as instead of mina, if the age of entry is different
    end
    @inbounds for a = 1:nper 
        @inbounds for ei = 2:maxngpe
            if maxe[a]<=auxge[ei]
                ge[ei,a] = maxe[a]
                ngpe[a] = ei
                break
            else
                ge[ei,a] = auxge[ei]
            end
            if ngpe[a] == 0
                ngpe[a] = maxngpe
            end     
        end
    end
    return ngpe, ge
end    

#constructing the grids for the assets
function popgk(fp::fparams)
    @unpack mink,maxk,ngpk,ratio = fp    
    ngpkf   = convert(Float64,ngpk)
    gk      = Array{Float64}(ngpk)
    gk[1]   = mink
    unit    = (maxk - mink)*(1.-ratio)/(1.-ratio^(ngpkf-1.))
    for ki = 2:ngpk
        gk[ki] = gk[ki] + unit*ratio^(ki-2)
    end
    return gk
end

#constructing the bounds for next year's experience based on this period's labor supply choice
function nxtebds(a::Int64, ngpl::Int64, maxngpe::Int64, ngpe::Array{Int64,1}, ge::Array{Int64,2})
    e1lb = Array{Int64}(ngpl,maxngpe)
    e1ub = Array{Int64}(ngpl,maxngpe)
    e1 = Array{Int64}(ngpl,maxngpe)
    a1 = a+1
    for ei = 1:ngpe[a] 
        #println("ei $ei")
        for l = 1:ngpl
            e1[l,ei] = ge[ei,a] + l-1
            #println("e1 $(e1[l,ei])")
            for nei = 1:ngpe[a1]
                #println("nei $nei")
                if abs(ge[nei,a1]-e1[l,ei]) <= 0.00001
                    e1lb[l,ei] = nei
                    e1ub[l,ei] = nei
                    #println("e1lb $(e1lb[l,ei])", "e1ub $(e1ub[l,ei])")
                    break
                elseif ge[nei,a1]>e1[l,ei]
                    e1lb[l,ei] = nei-1
                    e1ub[l,ei] = nei
                    #println("e1lb $(e1lb[l,ei])", "e1ub $(e1ub[l,ei])") 
                    break
                end
            end
        end
    end
    return e1, e1lb, e1ub
end

#this function determines the interval of a value of a variable belongs to given the grid strcuture
function bounds{T<:Real}(x::T,gx::Array{T,1},ngpx::Int)  
    #println("gx $gx")
    lb = 0
    ub = 0    
    i = 1  
    while true
        #println("i $i","\t lb $lb", "\t ub $ub")
        #println("x $x","\t gx[i] $(gx)")
        if abs(x-gx[i])<0.00001
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
function ilocate(i,rep,rvw,rvu,fp)    
    @unpack maxa, mina = fp
    irvw = Array{Float64}(maxa-mina+1)
    irvu = Array{Float64}(maxa-mina+1)
    lb = (i-1)*(maxa-mina+1)+1
    ub = i*(maxa-mina+1)     
    irvw = rvw[lb:ub,rep]     
    irvu = rvu[lb:ub,rep]
    return irvw, irvu
end


#calculating the probabilities for the grids of the random shocks
function trprob(gεw::Array{Float64,1}, gεu::Array{Float64,1}, p0::params, fp::fparams)
    @unpack ngu = fp
    @unpack σw, σu = p0
    derrw = Normal(0,σw)
    derru = Normal(0,σu)
    cdfw = Array{Float64}(ngu)
    cdfu = Array{Float64}(ngu)
    cdfwu = Array{Float64}(ngu,ngu)
    pgridwu = Array{Float64}(ngu,ngu)
    for i = 1:ngu
        cdfw[i] = cdf(derrw,gεw[i])
        cdfu[i] = cdf(derru,gεu[i])
    end
    for iu = 1:ngu
        for iv = 1:ngu
            cdfwu[iv,iu] = cdfw[iv]*cdfu[iu]
        end
    end

    for iu = 1:ngu        
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
    return pgridwu
end
