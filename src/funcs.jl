#function for linear interpolation
function LinInterp1d{T<:Real}(x::T, x1::T, x2::T, fn1::Float64, fn2::Float64)
    if x1 == x2 || fn1 == fn2
        return fn1
    else
        dx = (x-x1)/(x2-x1)
        return fn1+dx*(fn2-fn1)
    end
end   

#function for linear interpolation in 2 dimensions
function LinInterp2d{T1<:Real,T2<:Real}(x::T1, x1::T1, x2::T1, y::T2, y1::T2, y2::T2, fn11::Float64, fn12::Float64, fn21::Float64, fn22::Float64)
    if x1==x2
        return LinInterp1d(y,y1,y2,fn11,fn22)
    elseif y1==y2
        return LinInterp1d(x,x1,x2,fn11,fn22)
    else
        fn1 = LinInterp1d(y,y1,y2,fn11,fn12)    #interpolate in y for the 2 points in x
        fn2 = LinInterp1d(y,y1,y2,fn21,fn22)
        return LinInterp1d(x,x1,x2,fn1,fn2)     #interpolate in x
    end
end
  
#the utility function
function util(c::Float64, l::Int64, age::Int64, epsu::Float64, p::params, fp::fparams)
    @unpack αc, αh0, αh1, αh2, αh3, αh4 = p
    @unpack θ = fp
    lf = convert(Float64,l)
    exptheta = exp(αc*lf)/θ
    #println("ac*lf $(αc*lf)","\t exp(αc*lf) $(exp(αc*lf))","\t exp(αc*lf)/θ $(exp(αc*lf)/θ)", "\t exptheta $(exptheta)")
    #println("\n c $c")
    if age <= 25
        αh = αh0 + αh1
    elseif age > 25 & age <= 30
        αh = αh0 + αh2
    elseif age > 30 & age < 40
        αh = αh0
    elseif age >= 40 & age < 50
        αh = αh0 + αh3
    else
        αh = αh0 + αh4
    end
    #println("αh $(αh)","\t exptheta $(exptheta)", "\t c^θ $(c^θ)","\t epsu $epsu \n")
    return c^θ*exptheta+lf*(αh + epsu)
end

#the derivative of the utility function wrt consumption
function dudc(c::Float64, l::Int64, p::params, fp::fparams)
    @unpack αc = p
    @unpack θ = fp
    return (c^(θ-1))*exp(αc*convert(Float64,l))
end

#wage equation
function wage(e::Int64, l1::Int64, p::params, gεwiu::Float64)
    @unpack αw0, αw1, αw2 = p
    return exp(αw0 + αw1*log(convert(Float64,e)+1.0)+αw2*convert(Float64,l1) + gεwiu)
end

#approximate inverse of EdU 
function calcEDUinv(Edu::Float64, fp::fparams)
    @unpack θ = fp
    return Edu^(1./(θ-1.))
end

#approximate inverse of EV
function calcEVinv(EV::Float64, fp::fparams)
    @unpack θ = fp
    return (θ*EV)^(1./θ)    
end

#reversing the inverse of EV
function calcEVinvinv(EV::Float64, fp::fparams)
    @unpack θ = fp
    return (1./θ)*(EV^θ)    
end











                        






