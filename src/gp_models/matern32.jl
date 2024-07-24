
k_m32(x::Float64,y::Float64,A,tau)  = A*(1+sqrt(3.0)*abs(x-y)/tau)*exp(-sqrt(3.0)*abs(x-y)/tau)
dk_m32(x::Float64,y::Float64,A,tau) = -A*3.0*(x-y)*exp(-sqrt(3.0)*abs(x-y)/tau)/tau^2
ddk_m32(x::Float64,y::Float64,A,tau) =  A*3.0*exp(-sqrt(3.0)*abs(x-y)/tau)*(tau-sqrt(3.0)*abs(x-y))/tau^3


name(m::Matern32Model) = "matern32"
param_names(m::Matern32Model) = [:A_cell,:tau_cell,:A_flucs,:tau_flucs,:A_err]
lntheta0(m::Matern32Model) = log.([0.01,0.5,0.01,2.,0.01])
terms(m::Matern32Model) = ["cell","flucs","noise"]
deriv_terms(m::Matern32Model) = terms(m)[[1,2]]


function K(m::Matern32Model,term::String,x::Array{Float64,1},y::Array{Float64,1},lntheta)
    T = terms(m)
    theta = exp.(lntheta)
    if term == T[1]
        return k_se(x[2],y[2],theta[1],theta[2])
    elseif term == T[2]
        if x[3]==y[3] # if these come from the same cell
            return k_m32(x[1],y[1],theta[3],theta[4])
        else
            0. #return theta[3]*exp(-abs(x[3]-y[3])*2.) #k_m32(0.,abs(y[1]-x[1])+3.,theta[3],theta[4])
        end
    else
        return x==y ? theta[5] : 0
    end
end

function dK(m::Matern32Model,term::String,x::Array{Float64,1},y::Array{Float64,1},lntheta)
    T = terms(m)
    theta = exp.(lntheta)
    if term == T[1]
        return dk_se(x[2],y[2],theta[1],theta[2])
    elseif term == T[2]
        if x[3]==y[3] # if these come from the same cell
            return dk_m32(x[1],y[1],theta[3],theta[4])
        else
            return 0.
        end
    else
        return false
    end

end

function ddK(m::Matern32Model,term::String,x::Array{Float64,1},y::Array{Float64,1},lntheta)
    T = terms(m)
    theta = exp.(lntheta)
    if term == T[1]
        return ddk_se(x[2],y[2],theta[1],theta[2])
    elseif term == T[2]
        if x[3]==y[3] # if these come from the same cell
            return ddk_m32(x[1],y[1],theta[3],theta[4])
        else
            return 0.
        end
    else
        return false
    end
end

#---------------------------------------------------------------------------------------

name(m::Matern32NoTrendModel) = "matern32-notrend"
param_names(m::Matern32NoTrendModel) = [:A_flucs,:tau_flucs,:A_err]
lntheta0(m::Matern32NoTrendModel) = log.([0.01,2.,0.0001])
terms(m::Matern32NoTrendModel) = ["flucs","noise"]
deriv_terms(m::Matern32NoTrendModel) = terms(m)[[1]]


function K(m::Matern32NoTrendModel,term::String,x::Array{Float64,1},y::Array{Float64,1},lntheta)
    T = terms(m)
    theta = exp.(lntheta)
    if term == T[1]
        return k_m32(x[1],y[1],theta[1],theta[2])
    else
        return x==y ? theta[3] : 0
    end
end

function dK(m::Matern32NoTrendModel,term::String,x::Array{Float64,1},y::Array{Float64,1},lntheta)
    T = terms(m)
    theta = exp.(lntheta)
    if term == T[1]
        return dk_m32(x[1],y[1],theta[1],theta[2])
    else
        return false
    end

end

function ddK(m::Matern32NoTrendModel,term::String,x::Array{Float64,1},y::Array{Float64,1},lntheta)
    T = terms(m)
    theta = exp.(lntheta)
    if term == T[1]
        return ddk_m32(x[1],y[1],theta[1],theta[2])
    else
        return false
    end
end

