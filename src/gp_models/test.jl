

name(m::TestModel) = "test_model"
param_names(m::TestModel) = [:A1,:tau1,:A1,:tau2,:A_err]
lntheta0(m::TestModel) = log.([1.,4.,1.,0.2,0.3])
terms(m::TestModel) = ["k1","k2"]
deriv_terms(m::TestModel) = []


function K(m::TestModel,term::String,x::Array{Float64,1},y::Array{Float64,1},lntheta)
    T = terms(m)
    theta = exp.(lntheta)
    if term ==T[1]
        return k_se(x[2],y[2],theta[1],theta[2])
    elseif term == T[2]
        return k_ou(x[1],y[1],theta[3],theta[4])
    else
        return x==y ? theta[5] : 0
    end
end
