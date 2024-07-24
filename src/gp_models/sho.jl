function k_sho(x::Float64,y::Float64,D,Q,omega0)
    eta = sqrt(1-1/(4*Q^2))
    dt = abs(x[1]-y[1])
    if Q<1/2
        return cosh(eta*omega0*dt) + 1/(2*eta*Q)*sinh(eta*omega0*dt)
    elseif Q>1/2
        return cos(eta*omega0*dt) + 1/(2*eta*Q)*sinh(eta*omega0*dt)
    end
end


name(m::SHOModel) = "sho"
param_names(m::SHOModel) = [:A_cell,:tau_cell,:A_flucs,:tau_flucs,:A_err]
lntheta0(m::SHOModel) = log.([0.01,0.5,0.01,2.,0.0001])
terms(m::SHOModel) = ["cell","flucs","noise"]
deriv_terms(m::SHOModel) = []


function K(m::SHOModel,term::String,x::Array{Float64,1},y::Array{Float64,1},lntheta)
    T = terms(m)
    theta = exp.(lntheta)
    if term == T[1]
        return k_se(x[2],y[2],theta[1],theta[2])
    elseif term == T[2]
        return k_sho(x[1],y[1],theta[3],theta[4],theta[5])
    else
        return x==y ? theta[6] : 0
    end
end
