# Author: Laura Zhang
# Created July 2017
# This julia script will read in parameters for the
# Earned Income Tax Credit and insert them into
# a dictionary for later use

# change to your directory
#cd("/data/eitc_data_dir/juliafunc/laura")

# load Packages
using DataFrames
using HDF5, JLD

# read in state EITC csv file to DataFrame and clean
stateData = readcsv("StateEITC.csv", skipstart=21)
stateData = stateData[:,3:17]
stateData = convert.(Float64,stateData)


# read in federal EITC csv file to DataFrame and clean
fedData = readcsv("FedEITC.csv", skipstart=16)
fedData = fedData[:,3:8]
fedData = convert.(Float64,fedData)

# Create Dictionary where keys are either an array
# [policy, num kids] mapping to an array for federal 
# credit details [credit rate, min income for max credit, max credit,
# phaseout rate, begin income for phaseout, end income for phaseout]
# or an array [state, policy, num kids] mapping to 
# a float for state credit as a % of federal credit
# Ex. [1, 1] -> [34.0, 6160, 2094, 15.98, 11290, 24396]
# Ex. [9, 2, 1998]-> 10.0
# Returns object of type Dict{Array, Union{Array,Float64} }

#=

fedtax=[0   12.07547588 20.07815812 43.69321598 63.39624837 102.8023132 0.11    0.13    0.18    0.27    0.38    0.48;
0   0   0   0   29.3717655  71.156169   0   0   0   0.15    0.28    0.31;
0   0   31.3165281  75.8514724  158.2503846 343.9876627 0   0.15    0.28    0.31    0.36    0.396;
0   7.5369726   30.6641214  74.2856409  155.0241693 337.0472199 0.1 0.15    0.25    0.28    0.33    0.35;
0   12.07547588 21.93592364 53.53580061 83.09928076 152.0688257 0.11    0.13    0.18    0.27    0.38    0.48;
0   0   0   0   41.4957375  107.0229195 0   0   0   0.15    0.28    0.31;
0   0   41.9406757  108.3415982 175.421972  343.9876627 0   0.15    0.28    0.31    0.36    0.396;
0   10.7892279  41.0919876  106.1370936 171.8532999 337.0472199 0.1 0.15    0.25    0.28    0.33    0.35]



#fedtax = [0 4.26928807  6.32354802  8.18131354  12.07547588 15.80887005 20.07815812 23.99018359 27.88434593 33.83276822 43.69321598 53.53580061 63.39624837 77.15085847 102.8023132 152.0688257 0   0.11    0.12    0.14    0.15    0.16    0.18    0.2 0.23    0.26    0.3 0.34    0.38    0.42    0.48    0.5;
#    0   0   0   0   0   0   0   0   0   0   0   0   0   0   29.3717655  71.156169   0   0   0   0   0   0   0   0   0   0   0   0   0   0.15    0.28    0.31;
#    0   0   0   0   0   0   0   0   0   0   0   0   31.3165281  75.8514724  158.2503846 343.9876627 0   0   0   0   0   0   0   0   0   0   0   0.15    0.28    0.31    0.36    0.396;
#    0   0   0   0   0   0   0   0   0   0   0   7.5369726   30.6641214  74.2856409  155.0241693 337.0472199 0   0   0   0   0   0   0   0   0   0   0.1 0.15    0.25    0.28    0.33    0.35;
#    0   0   4.26928807  8.18131354  12.07547588 16.16613265 21.93592364 27.88434593 33.83276822 43.69321598 53.53580061 63.39624837 83.09928076 112.6627609 152.0688257 201.3532014 0   0   0.11    0.12    0.14    0.17    0.18    0.2 0.24    0.28    0.32    0.35    0.42    0.45    0.48    0.5;
#    0   0   0   0   0   0   0   0   0   0   0   0   0   0   41.4957375  107.0229195 0   0   0   0   0   0   0   0   0   0   0   0   0   0.15    0.28    0.31;
#    0   0   0   0   0   0   0   0   0   0   0   0   41.9406757  108.3415982 175.421972  343.9876627 0   0   0   0   0   0   0   0   0   0   0   0.15    0.28    0.31    0.36    0.396;
#    0   0   0   0   0   0   0   0   0   0   0   10.7892279  41.0919876  106.1370936 171.8532999 337.0472199 0   0   0   0   0   0   0   0   0   0   0.1 0.15    0.25    0.28    0.33    0.35]


#fedtax = fedtax'
fedexemp = [1.85776552, 3.1031595, 3.3354882, 3.3038784]
fedded = [4.26928807  4.26928807;
    4.907322    7.21665;
    5.2503055   7.7210375;
    5.16231  7.5369726]


type datadict2
    ix::Array{Int64,2}
    ex::Array{Float64,1}
    b::Array{Float64,1}
    t::Array{Float64,1}
end



function createDictTaxes(fedtax::Array{Float64,2}, fedexemp::Array{Float64,1})
    taxDict = datadict2(zeros(Int64,4,2),zeros(Float64,4),zeros(Float64,48),zeros(Float64,48))
    ix = 0
    taxDict.ex = fedexemp
    for si = 1:2, pli = 1:4
        ix = ix +1
        taxDict.ix[pli,si] = ix
        taxDict.b[ix:ix+5] = fedtax[(si-1)*4+pli,1:6]
        taxDict.t[ix:ix+5] = fedtax[(si-1)*4+pli,7:12]
    end
    return taxDict
end
        
taxDict = createDictTaxes(fedtax,fedexemp)
save("./dicts/taxDict.jld", "data", taxDict)
=#


type datadict
    ix :: Array{Int,2}
    ch :: Array{Int,1}
    policy :: Array{Int,1}
    steitc :: Array{Float64,2}
    fedeitc:: Array{Float64,2}
end

function createDictEITC(stateData::Array{Float64,2}, fedData::Array{Float64,2})
    eitcDict = datadict(zeros(Int,4,3),zeros(Int,12),zeros(Int,12),zeros(Float64,12,15),zeros(Float64,12,6))
    ix = 0 
    for pli = 1:4, chi = 1:3
        ix = ix + 1
        eitcDict.ix[pli,chi] = ix
        eitcDict.ch[ix] = chi-1
        eitcDict.policy[ix] = pli
        eitcDict.steitc[ix,:] = stateData[ix,:]
        eitcDict.fedeitc[ix,:] = fedData[ix,:]
    end
    return eitcDict
end


eitcDict = createDictEITC(stateData, fedData)
save("eitcDict.jld", "data", eitcDict)


ctc = [0.4941464    0   0   67.94513    92.65245    135.89026   0.05;
    1.032462    10.32462    0.15    56.78541    77.43465    113.57082   0.05]   

type datadict3
    ix :: Array{Int64,1}
    phthr :: Array{Float64,1}
    maxcr :: Array{Float64,1}
    phrate :: Array{Float64,1}
end

function createDictCTC(ctc::Array{Float64,2})
    ctcDict = datadict3(zeros(Int64,2),zeros(Float64,2),zeros(Float64,2),zeros(Float64,2))
    ix = 0    
    for pli = 1:2
        ix +=1 
        ctcDict.ix[pli] = ix
        ctcDict.phthr[ix] = ctc[pli,5]
        ctcDict.maxcr[ix] = ctc[pli,1]
        ctcDict.phrate[ix] = ctc[pli,7]
    end
    return ctcDict
end

ctcDict = createDictCTC(ctc)
save("ctcDict.jld", "data", ctcDict)


type datadict4
    ix :: Array{Int64,2}
    b0 :: Array{Float64,1}
    b1 :: Array{Float64,1}
    b2 :: Array{Float64,1}
end


function createDictSttax()
    sttax = datadict4(zeros(Int64,4,15),zeros(Float64,60),zeros(Float64,60),zeros(Float64,60))
    ix = 0
    data = readdlm("coeffs_sttax_polyfit.txt",'\t')
    for st = 1:15,pl = 1:4
        ix = ix + 1
        sttax.ix[pl,st] = ix
        sttax.b0[ix] = data[(st-1)*4+pl,3]
        sttax.b1[ix] = data[(st-1)*4+pl,4]
        sttax.b2[ix] = data[(st-1)*4+pl,5]
    end
    return sttax
end

stDict = createDictSttax()
save("stDict.jld","data",stDict)


type datadict5
    b0::Array{Float64,1}
    b1::Array{Float64,1}
    b2::Array{Float64,1}
end


function createDictFedtax()
    fedtax = datadict5(zeros(Float64,4),zeros(Float64,4),zeros(Float64,4))
    data = readdlm("coeffs_fedtax_polyfit.txt",'\t')
    for pl = 1:4
        fedtax.b0[pl] = data[pl,2]
        fedtax.b1[pl] = data[pl,3]
        fedtax.b2[pl] = data[pl,4]
    end
    return fedtax
end

fedDict = createDictFedtax()
save("ftDict.jld","data",fedDict)





