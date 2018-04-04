#an intermediate function to calculate the data moments
function CalcSimMom!(sim::Array{sim_t,2},fp::fparams, moments::structureM,version::Int, gender::String)
    @unpack mina, maxa, nper, nind, nsim = fp    
    #Hours and Employment Related
    f = 0    
    #moments 1-32
    for j = 2:5 #hours choice
        for s in 1:4 #schooling
            ind = [sim[i].s .== s for i in eachindex(sim)]
            x1 = count(x->(x==j), sim[i].aa.Ch.l[a] for a in 1:nper, i in eachindex(sim[ind]) if sim[i].aa.Oh.insample[a]==1) #those who choose the hours choice
            x2 = count(x->(x!=j), sim[i].aa.Ch.l[a] for a in 1:nper, i in eachindex(sim[ind]) if sim[i].aa.Oh.insample[a]==1) #those who don't choose the hours choice
            f+=1
            moments.simmom[f] = x1/(x1+x2)            
        end
        x1 = count(x->(x==j), sim[i].aa.Ch.l[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1) #those who choose the hours choice
        x2 = count(x->(x!=j), sim[i].aa.Ch.l[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1) #those who don't
        f+=1       
        moments.simmom[f] = x1/(x1+x2)
        for ch = 1:3
            x1 = 0
            x2 = 0            
            for i in eachindex(sim)
                ind = [sim[i].aa.Oh.ch[a] == ch && sim[i].aa.Oh.insample[a] ==1 for a in 1:nper]
                x1 += count(x->(x==j), sim[i].aa.Ch.l[ind]) #those who choose the hours choice
                x2 += count(x->(x!=j), sim[i].aa.Ch.l[ind]) #those who don't choose the hours choice
            end
            f +=1
            moments.simmom[f] = x1/(x1+x2)
        end
    end
    println("f0 $(f)")
    #Employment Transitions
    #moments 33-35
    x1 = 0
    x2 = 0
    for i in eachindex(sim)
        ind = [sim[i].aa.Oh.l_1[a] == 1 && sim[i].aa.Oh.insample[a] == 1 for a in 1:nper]
        x1 += count(x->(x!=1), sim[i].aa.Ch.l[ind]) #employed
        x2 += count(x->(x==1), sim[i].aa.Ch.l[ind]) #unemployed
    end
    f+=1
    if x1 == 0 && x2 == 0
        moments.simmom[f] = -99999
    else
        moments.simmom[f] = x1/(x1+x2)
    end
    x1 = 0
    x2 = 0
    for i in eachindex(sim)
        ind = [sim[i].aa.Oh.l_1[a] != 1 && sim[i].aa.Oh.insample[a] == 1 for a in 1:nper]        
        x1 += count(x->(x==1), sim[i].aa.Ch.l[ind]) #unemployed
        x2 += count(x->(x!=1), sim[i].aa.Ch.l[ind]) #employed
    end
    f+=1
    if x1 == 0 && x2 == 0
        moments.simmom[f] = -99999
    else
        moments.simmom[f] = x1/(x1+x2)
    end
    #correlation between hours and last period's hours
    a = [sim[i].aa.Ch.l[a] for a in 2:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1 && sim[i].aa.Oh.insample[a-1]==1]
    b = [sim[i].aa.Oh.l_1[a] for a in 2:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1 && sim[i].aa.Oh.insample[a-1]==1]
    f += 1
    if std(a)==0. || std(b)==0
        moments.simmom[f] = -99999
    else
        moments.simmom[f] = cor(a,b)
    end

    #println("f1 $(f)")
    #State-level moments    
    #moments 36-115
    #for men, also drop minnesota for year<=1986
    stlist = [2,3,5,6,7,8,9,11,14,15]
    if gender == "f"            
        for pl = 1:4
            for st in stlist
                ind = [sim[i].st .== st for i in eachindex(sim)]
                f+=1        
                if isempty([sim[i].aa.Ch.hrs[a] for a in 1:nper, i in eachindex(sim[ind]) if sim[i].aa.Oh.insample[a]==1 
                    && sim[i].aa.Oh.policy[a]==pl && sim[i].aa.Ch.hrs[a]>0])
                    moments.simmom[f] = -99999
                else
                    moments.simmom[f] = (mean(sim[i].aa.Ch.hrs[a] for a in 1:nper, i in eachindex(sim[ind]) if sim[i].aa.Oh.insample[a]==1 
                        && sim[i].aa.Oh.policy[a]==pl && sim[i].aa.Ch.hrs[a]>0))
                end
                x1 = count(x->(x==1), sim[i].aa.Ch.l[a] for a in 1:nper, i in eachindex(sim[ind]) if sim[i].aa.Oh.insample[a]==1 && sim[i].aa.Oh.policy[a]==pl) #those who choose the hours choice
                x2 = count(x->(x!=1), sim[i].aa.Ch.l[a] for a in 1:nper, i in eachindex(sim[ind]) if sim[i].aa.Oh.insample[a]==1 && sim[i].aa.Oh.policy[a]==pl) #those who don't choose the hours choice
                f+=1                
                moments.simmom[f] = x1/(x1+x2)
            end
        end
    else
        for pl = 1:4
            for st in stlist
                ind = [sim[i].st .== st for i in eachindex(sim)]
                f+=1        
                moments.simmom[f] = (mean[sim[i].aa.Ch.hrs[a] for a in 1:nper, i in eachindex(sim[ind]) if sim[i].aa.Oh.insample[a]==1 
                    && sim[i].aa.Oh.policy[a]==pl && sim[i].aa.Oh.hrs[a]>0])
                if st != 7 && pl!=1
                    x1 = count(x->(x==1), sim[i].aa.Ch.l[a] for a in 1:nper, i in eachindex(sim[ind]) if sim[i].aa.Oh.insample[a]==1 && sim[i].aa.Oh.policy[a]==pl) #those who choose the hours choice
                    x2 = count(x->(x!=1), sim[i].aa.Ch.l[a] for a in 1:nper, i in eachindex(sim[ind]) if sim[i].aa.Oh.insample[a]==1 && sim[i].aa.Oh.policy[a]==pl) #those who don't choose the hours choice
                    f+=1
                    moments.simmom[f] = x1/(x1+x2)
                end
            end
        end
    end

    #println("f2 $(f)")
    #WAGES
    aget = [1:30-mina+1, 30-mina+2:40-mina, 40-mina+1:50-mina, 50-mina+1:nper, 1:nper]
    aget2 = [31-mina, 41-mina-32+mina, 51-mina-41+mina, nper+1-51+mina, nper]
    for (m,n) in enumerate(aget)
        f += 1
        moments.simmom[f] = NaNMath.mean(vec(reshape([sim[i].aa.Oh.w[j] for j in n, i in eachindex(sim) if sim[i].aa.Oh.insample[j]==1],:,1)))        
    end            
    f+=1
    moments.simmom[f] = NaNMath.std(vec(reshape([sim[i].aa.Oh.w[j] for j in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[j]==1],:,1)))        
    for l = 3:5
        f+=1
        if isempty(vec(reshape([sim[i].aa.Oh.w[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1 && 
            sim[i].aa.Ch.l[a] == l],:,1)))
            moments.simmom[f] = -99999
        else
            moments.simmom[f] = (NaNMath.mean(vec(reshape([sim[i].aa.Oh.w[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1 && 
                sim[i].aa.Ch.l[a] == l],:,1))))
        end
    end
    for l_1 = 2:3
        f+=1
        if isempty(vec(reshape([sim[i].aa.Oh.w[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1 && 
            sim[i].aa.Oh.l_1[a] == l_1],:,1)))
            moments.simmom[f] = -99999
        else
            moments.simmom[f] = (NaNMath.mean(vec(reshape([sim[i].aa.Oh.w[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1 && 
                sim[i].aa.Oh.l_1[a] == l_1],:,1))))
        end
    end
    
    for s = 1:4
        ind = [sim[i].s .== s for i in eachindex(sim)]
        f += 1
        moments.simmom[f] = NaNMath.mean(vec(reshape([sim[i].aa.Oh.w[a] for a in 1:nper, i in eachindex(sim[ind]) if sim[i].aa.Oh.insample[a]==1],:,1)))
    end
    println("f3 $(f)")
    #correlations between age and wage
    a = [sim[i].aa.Oh.w[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1]
    b = [sim[i].aa.Oh.a[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1]
    f += 1
    if std(a)==0. || std(b)==0
        moments.simmom[f] = -99999
    else
        moments.simmom[f] = cor(a,b)
    end
    #correlations between wage and the number of kids
    a = [sim[i].aa.Oh.w[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1]
    b = [sim[i].aa.Oh.ch[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1]
    f += 1
    if std(a)==0. || std(b)==0
        moments.simmom[f] = -99999
    else
        moments.simmom[f] = cor(a,b)
    end
    #correlations between wage and hours
    a = [sim[i].aa.Oh.w[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1]
    b = [sim[i].aa.Ch.hrs[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1]
    f += 1
    if std(a)==0. || std(b)==0
        moments.simmom[f] = -99999
    else
        moments.simmom[f] = cor(a,b)
    end 
    #correlations between wage and full-time experience
    a = [sim[i].aa.Oh.w[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1]
    b = [sim[i].aa.Oh.fte[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1]    
    f += 1
    if std(a)==0. || std(b)==0
        moments.simmom[f] = -99999
    else
        moments.simmom[f] = cor(a,b)
    end
    #println("f4 $(f)")
    #moments with experience quantiles
    #pteq = [0 0 0 1; 1 1 1 4; 2 3 3 5; 6 5 6 7]
    fteq = [0 1 3 4; 1 7 8 9; 13 17 15 18; 26 25 23 24]
    for s = 1:4
        ind = [sim[i].s .== s for i in eachindex(sim)]
        for q = 1:4 #loop over quartiles
            f+=1
            moments.simmom[f] = (NaNMath.mean(vec(reshape([sim[i].aa.Oh.w[a] for a in 1:nper, i in eachindex(sim[ind]) if sim[i].aa.Oh.insample[a]==1 &&
                sim[i].aa.Oh.fte[a]<=fteq[q,s]],:,1))))
        end
    end
    #println("f5 $(f)")
    if version == 1
        #take-up
        for ch = 1:3
            x1 = count(x->(x==2), sim[i].aa.Ch.p[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1 && sim[i].aa.Ch.l[a]>1 && sim[i].aa.Oh.ch[a]==ch) #those who participate in the EITC
            x2 = count(x->(x!=2), sim[i].aa.Ch.p[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1 && sim[i].aa.Ch.l[a]>1 && sim[i].aa.Oh.ch[a]==ch) #those who don't participate        
            f+=1
            moments.simmom[f] = x1/(x1+x2)            
        end
    elseif version == 2    
        #Persistence - don't include these if you don't have unobserved heterogeneity
        for s = 1:4
            for l = 2:5
                m1 = Array(Float64,length(sim))
                ind = [sim[i].s .== s for i in eachindex(sim)]
                for i in eachindex(sim[ind])                
                    x1 = count(x->(x==l), sim[i].aa.Ch.l[a] for a in 1:nper if sim[i].aa.Oh.insample[a]==1)
                    x2 = count(x->(x!=l), sim[i].aa.Ch.l[a] for a in 1:nper if sim[i].aa.Oh.insample[a]==1)
                    m1[i] = x1/(x1+x2)
                end
                f += 1
                moments.simmom[f] = mean(m1)
            end
        end
        #take-up
        for ch = 1:3
            x1 = count(x->(x==2), sim[i].aa.Ch.p[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1 && sim[i].aa.Ch.l[a]>1 && sim[i].aa.Oh.ch[a]==ch) #those who participate in the EITC
            x2 = count(x->(x!=2), sim[i].aa.Ch.p[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1 && sim[i].aa.Ch.l[a]>1 && sim[i].aa.Oh.ch[a]==ch) #those who don't participate
            f+=1
            moments.simmom[f] = x1/(x1+x2)            
        end
    else #version3
        #Persistence - don't include these if you don't have unobserved heterogeneity
        for s = 1:4
            for l = 2:5
                m1 = Array(Float64,length(sim))
                ind = [sim[i].s .== s for i in eachindex(sim)]
                for i in eachindex(sim[ind])                
                    x1 += count(x->(x==l), sim[i].aa.Ch.l[a] for a in 1:nper if sim[i].aa.Oh.insample[a]==1)
                    x2 += count(x->(x!=l), sim[i].aa.Ch.l[a] for a in 1:nper if sim[i].aa.Oh.insample[a]==1)
                    m1[i] = x1/(x1+x2)
                end
                f += 1
                moments.simmom[f] = mean(m1)
            end
        end
        #Net worth - don't include these if you are not matching net worth moments
        f+=1
        moments.simmom[f] = NaNMath.mean(vec(reshape([sim[i].aa.Oh.k[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insampleass[a]==1],:,1)))
        f+=1
        moments.simmom[f] = NaNMath.std(vec(reshape([sim[i].aa.Oh.k[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insampleass[a]==1],:,1)))
        #percentiles
        for j = 2:5
            f+=1
            moments.simmom[f] = quantile([sim[i].aa.Oh.k[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insampleass[a]==1],0.1)        
            f+=1
            moments.simmom[f] = quantile([sim[i].aa.Oh.k[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insampleass[a]==1],0.5)
            f+=1
            moments.simmom[f] = quantile([sim[i].aa.Oh.k[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insampleass[a]==1],0.9)                
        end
        #take-up
        for ch = 1:3
            x1 = count(x->(x==2), sim[i].aa.Ch.p[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1 && sim[i].aa.Ch.l[a]>1 && sim[i].aa.Oh.ch[a]==ch) #those who participate in the EITC
            x2 = count(x->(x!=2), sim[i].aa.Ch.p[a] for a in 1:nper, i in eachindex(sim) if sim[i].aa.Oh.insample[a]==1 && sim[i].aa.Ch.l[a]>1 ) #those who don't participate
            f+=1
            moments.simmom[f] = x1/(x1+x2)            
        end
    end

    difference = moments.dtamom - moments.simmom
    moments.objfunc = sum((difference.^2).*(moments.wgtcov))        
    return nothing
end
