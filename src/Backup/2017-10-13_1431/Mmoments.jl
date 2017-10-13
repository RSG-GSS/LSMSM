##########################################################################   
##### 10/13/2017: Modified by Miura, Hirotaka <Hirotaka.Miura@ny.frb.org>
##### 00/00/0000: Previously modified by 
##### 10/13/2017: Created by Kosar, Gizem <Gizem.Kosar@ny.frb.org>
##########################################################################
##### Description: 
##### 	- Module for calculating and manipulating moments.
##### Modifications:
##########################################################################
##### Define module name.
module Mmoments
##### Import packages.
using Parameters		##### For @unpack.
using Mparams				##### For params.
import NaNMath
##### Specify names to be exported.
export
CalcMom,
CalcDataBootMom!,
CalcSimMom!

#an intermediate function to calculate the data moments
function CalcMom(data::Array{Float64,2}, cols::Dict{Symbol,Int})
    mom = Array(Float64,25)
    #EMPLOYMENT      
    aget = [Array(data[:,cols[:age]].<=60), Array(data[:,cols[:age]].<=25), Array(data[:,cols[:age]].<=30), Array(data[:,cols[:age]].>30) & Array(data[:,cols[:age]].<=40), Array(data[:,cols[:age]].>40) & Array(data[:,cols[:age]].<=50), Array(data[:,cols[:age]].>50) ]
    f = 1
    for (m,n) in enumerate(aget)
        x1 = count(x->(x==2), data[Array(data[:,cols[:actual]] .== 1) & n,cols[:ls]])
        x2 = count(x->(x==1), data[Array(data[:,cols[:actual]] .== 1) & n,cols[:ls]])
        mom[m] = x1/(x1+x2)
        f += 1
    end
    #TRANSITIONS 
    x1 = count(x->(x==1), data[Array(data[:,cols[:actual]] .== 1) & Array(data[:,cols[:ls_1]].==1),cols[:ls]])
    x2 = count(x->(x==2), data[Array(data[:,cols[:actual]] .== 1) & Array(data[:,cols[:ls_1]].==1),cols[:ls]])
    if x2 == 0
        mom[f] = 50
    else
        mom[f] = x1/(x1+x2)
    end
    f += 1
    x1=0
    x2=0
    x1 = count(x->(x==2), data[Array(data[:,cols[:actual]] .== 1) & Array(data[:,cols[:ls_1]].==2),cols[:ls]])
    x2 = count(x->(x==1), data[Array(data[:,cols[:actual]] .== 1) & Array(data[:,cols[:ls_1]].==2),cols[:ls]])
    if x2 == 0
        mom[f] = 50
    else
        mom[f] = x1/(x1+x2)
    end
    f +=1 
    #WAGES
    for n in aget
        mom[f] = NaNMath.mean(data[Array(data[:,cols[:actual]] .== 1) & n,cols[:wage]])
        f += 1
    end
    mom[f] = NaNMath.std(data[data[:,cols[:actual]].==1,cols[:wage]])
    f += 1
    #mean wage by experience quantiles
    #from the simulated data --> q25 = 5, q50 = 10, q75 = 16, q90 = 19
    expt = [Array(data[:,cols[:exper]].<=5), Array(data[:,cols[:exper]].<=10), Array(data[:,cols[:exper]].<=16), Array(data[:,cols[:exper]].<=19)] #will be calculated from data
    for n in expt
        mom[f] = NaNMath.mean(data[Array(data[:,cols[:actual]] .== 1) & n,cols[:wage]])
        f += 1
    end
    #ASSETS
    aget = [Array(data[:,cols[:age]].<=60), Array(data[:,cols[:age]].<=30), Array(data[:,cols[:age]].>30) & Array(data[:,cols[:age]].<=40), Array(data[:,cols[:age]].>40) & Array(data[:,cols[:age]].<=50), Array(data[:,cols[:age]].>50) ]
    for n in aget
        mom[f] = NaNMath.mean(data[Array(data[:,cols[:actual]] .== 1) & n,cols[:assets]])
        f += 1
    end
    mom[f] = NaNMath.std(data[data[:,cols[:actual]].==1,cols[:assets]])
    f += 1
    return mom
end    


#the function to calculate the data moments as well as the bootstrapped data moments 
#for the weight matrix of the objective function. 
function CalcDataBootMom!(fp::fparams, moments::structureM)
    println("calculating bootstrapped data moment weight matrix")
    @unpack mina, maxa, nper, nind, nsim, nboot, nmom = fp
    data = readdlm("data.txt",'\t')
    data[data.==-9] = NaN
    data_cols = [:id, :age, :actual, :ls, :ls_1, :wage, :assets, :exper]                 
    cols = Dict{Symbol,Int}()
    for (i,k) in enumerate(data_cols); cols[k] = i end
    moments.dtamom = CalcMom(data,cols)
    #bootstrap
    srand(1136)
    draws = rand(1:nind,nind,nboot)
    boot_moments = Array(Float64,nmom,nboot)
    temp = Array(Float64,nmom)
    for b = 1:nboot
        #randomly re-draw with replacement from the original data
        boot_data = Array(Float64,0,8)
        rand_order = draws[:,b]
        for (qi,q) in enumerate(rand_order)
            #block of women with id==rand_order[q]
            boot_data = vcat(boot_data,data[data[:,cols[:id]].==q,:])
        end
        #calculate moments for the data draw
        boot_moments[:,b] = CalcMom(boot_data,cols)        
    end
    #which elements of the moment vector has any variation?
    #some moments defined have NaN since out of age range etc.
    for i in 1:size(boot_moments,1)
        temp[i] = std(boot_moments[i,:])
    end
    moments.withvar = (temp.>0)
    boot_moments = boot_moments[temp.>0,:]
    #create bootstrapped cov matrix, off-diagonal elements set to zero.
    for b=1:nboot
        difference = boot_moments[:,b] - reshape(mean(boot_moments,2),size(boot_moments,1))
        if b == 1
            moments.wgtcov = 1/nboot.*difference.^2
        else
            moments.wgtcov += 1/nboot.*difference.^2
        end
    end
    println("total number of moments $(size(boot_moments,1))")
    return nothing
end

#the function to calculate the simulated moments.
function CalcSimMom!(sim::Array{sim_t,2},fp::fparams, moments::structureM)
    @unpack mina, maxa, nper, nind, nsim = fp
    #EMPLOYMENT  
    aget = [1:nper, 1:25-mina+1, 1:30-mina+1, 30-mina+2:40-mina+1, 40-mina+2:50-mina+1, 50-mina+2:nper]    
    f = 1
    for (m,n) in enumerate(aget)
        x1 = count(x->(x==2), sim[i].aa.Ch.l[j] for j in n, i in eachindex(sim)) #those who work
        x2 = count(x->(x==1), sim[i].aa.Ch.l[j] for j in n, i in eachindex(sim)) #those who work
        moments.simmom[m] = x1/(x1+x2)
        f += 1
    end

    #TRANSITIONS    
    x1 = 0
    x2 = 0
    for i in eachindex(sim)
        ind = sim[i].aa.Oh.l_1 .== 1
        x1 += count(x->(x==1), sim[i].aa.Ch.l[ind])
        x2 += count(x->(x==2), sim[i].aa.Ch.l[ind])
    end
    if x2 ==0
        moments.simmom[f] =50.
    else
         moments.simmom[f] = x1/(x1+x2)
    end
    f +=1
    x1 = 0 
    x2 = 0
    for i in eachindex(sim)
        ind = sim[i].aa.Oh.l_1 .== 2
        x1 += count(x->(x==2), sim[i].aa.Ch.l[ind])
        x2 += count(x->(x==1), sim[i].aa.Ch.l[ind])
    end
    if x2 ==0
        moments.simmom[f] =50.
    else
        moments.simmom[f] = x1/(x1+x2)
    end
    f += 1

    #WAGES
    aget2 = [nper, 26-mina, 31-mina, 42-mina-32+mina, 52-mina-42+mina, nper+1-52+mina]
    for (m,n) in enumerate(aget)
        moments.simmom[f] = NaNMath.mean(vec(reshape([sim[i].aa.Oh.w[j] for j in n, i in eachindex(sim)],aget2[m]*nsim*nind,1)))
        f += 1
    end
    moments.simmom[f] = NaNMath.std(vec(reshape([sim[i].aa.Oh.w[j] for j in 1:nper, i in eachindex(sim)],nper*nsim*nind,1)))
    f +=1
    #mean wage by experience quantiles
    expt3 = [5, 10, 16, 19]
    for (j,k) in enumerate(expt3)
        m1 = Array(Float64,length(sim))
        for i in eachindex(sim)
            expind1 = sim[i].aa.Oh.e.<=k
            m1[i] = NaNMath.mean(sim[i].aa.Oh.w[expind1])
        end
        moments.simmom[f] = mean(m1)
        f+=1
    end
    #ASSETS
    aget = [1:nper, 1:30-mina+1, 30-mina+2:40-mina+1, 40-mina+2:50-mina+1, 50-mina+2:nper]
    aget2 = [nper, 31-mina, 42-mina-32+mina, 52-mina-42+mina, nper+1-52+mina]
    for (m,n) in enumerate(aget)
        moments.simmom[f] = NaNMath.mean(vec(reshape([sim[i].aa.Oh.k[j] for j in n, i in eachindex(sim)],aget2[m]*nsim*nind,1)))
        f += 1
    end        
    moments.simmom[f] = NaNMath.std(vec(reshape([sim[i].aa.Oh.k[j] for j in 1:nper, i in eachindex(sim)],nper*nsim*nind,1)))
    return nothing
end

##### End module definition.
end