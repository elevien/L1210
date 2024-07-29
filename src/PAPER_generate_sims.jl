# ----------------------------------------------------------------------------------------------------------------
# SETUP
using Pkg
println(Pkg.status())
using CSV
using StatsBase
using DataFrames
using PythonPlot
using DifferentialEquations
using LinearAlgebra
using Distributions

include("./simulation_models.jl")


function solver_output_to_dataframe(sol)
    u = hcat(sol.u...)';
    df = DataFrame(u,["M","λ","Mf","position"])
    df[:,"time"] = sol.t;
    positions = unique(df.position)
    df[:,:age] = vcat([df[df.position .==p,:time] .- df[df.position .==p,:time][1] for p in positions]...)
    df[:,:age_normed] = vcat([df[df.position .==p,:age]./df[df.position .==p,:age][end] for p in positions]...);
    df[:,:age_rounded] = round.(df[:,:age_normed],digits=1);
    return df
end   

# ----------------------------------------------------------------------------------------------------------------
# PARAMS

λ0 = 0.06
CV = .07

τ = 10
σ_tot = .07 *λ0    # overall growth rate noise
Δ = 40                  # size increment added 
σ_M = sqrt(40)           # cell control noise
tmax = 9000


# initial conditions for all lineages
M0 = Δ
Mf0 = M0 + Δ
init = [M0,λ0,Mf0,0]


# ----------------------------------------------------------------------------------------------------------------
# RUN SIMULATIONS

# loop through tau values
τ_range = collect(1:2:10)
ϕ_range = collect(0.0001:0.4:0.99)
for τ in τ_range
    for ϕ in ϕ_range
        σ_NCS = sqrt(ϕ)*σ_tot
        D = σ_NCS^2/τ
        σ_DN  = sqrt((1-ϕ)*2*13/τ)*σ_tot
        θ = [Δ,σ_DN,τ,λ0,D,σ_M]
        println("running τ = "*string(τ)*" and ϕ = "*string(ϕ))
        prob,callback = build_model(θ,init,(0,tmax))
        @time sol = solve(prob,callback = callback)
        df = solver_output_to_dataframe(sol)

        cd(dirname(@__FILE__))
        folder = "./output/sims/"
        mkpath(folder)
        CSV.write(folder*"/sim_data_tau="*string(τ)*"-phi="*string(ϕ)*".csv",df)
    end
end



