

name(m::OUModel) = "ou"
param_names(m::OUModel) = [:A_cell,:tau_cell,:A_flucs,:tau_flucs,:A_err]
lntheta0(m::OUModel) = log.([0.01,0.5,0.01,2.,0.0001])
terms(m::OUModel) = ["cell","flucs","noise"]
deriv_terms(m::OUModel) = [] #terms(m)[[2,3]]


function K(m::OUModel,term::String,x::Array{Float64,1},y::Array{Float64,1},lntheta)
    T = terms(m)
    theta = exp.(lntheta)
    if term ==T[1]
        return k_se(x[2],y[2],theta[1],theta[2])
    elseif term == T[2]
        return k_iou(x[1],y[1],theta[3],theta[4])
    else
        return x==y ? theta[5] : 0
    end
end

# function dK(m::OUModel,term::String,x::Array{Float64,1},y::Array{Float64,1},lntheta)
#     T = terms(m)
#     theta = exp.(lntheta)
#     if term ==T[1]
#         return false
#     elseif term == T[2]
#         return dk_se(x[2],y[2],theta[2],theta[3])
#     elseif term == T[3]
#         return dk_iou(x[1],y[1],theta[4],theta[5])
#     else
#         return false
#     end
#
# end

# function ddK(m::OUModel,term::String,x::Array{Float64,1},y::Array{Float64,1},lntheta)
#     T = terms(m)
#     theta = exp.(lntheta)
#     if term ==T[1]
#         return false
#     elseif term == T[2]
#         return ddk_se(x[2],y[2],theta[2],theta[3])
#     elseif term == T[3]
#         return ddk_iou(x[1],y[1],theta[4],theta[5])
#     else
#         return false
#     end
# end
