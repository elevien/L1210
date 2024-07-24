
########################################################################################
# Model types
########################################################################################

struct Matern32Model <: LineageAnalysisTools.GpDecomp end
struct TestModel <: LineageAnalysisTools.GpDecomp end
struct OUModel <: LineageAnalysisTools.GpDecomp end 
struct Matern32NoTrendModel <: LineageAnalysisTools.GpDecomp end
struct CellSpecificMatern32Model <: LineageAnalysisTools.GpDecomp end


########################################################################################
# some useful kernel functions
########################################################################################

k_se(x::Float64,y::Float64,A,tau)  = A*exp(-(x-y)^2.0/tau^2)
dk_se(x::Float64,y::Float64,A,tau) = -2.0*A*exp(-(x-y)^2.0/tau^2)*(x-y)/tau^2
ddk_se(x::Float64,y::Float64,A,tau) = -4.0*A*exp(-(x-y)^2.0/tau^2)*(x-y)^2.0/tau^2.0 + A*2.0*exp(-(x-y)^2.0/tau^2)/tau
k_ou(x::Float64,y::Float64,A,tau)  = A*exp(-abs(x-y)/tau)


########################################################################################
# models
########################################################################################

include("test.jl")
include("ou.jl")
include("matern32.jl")
