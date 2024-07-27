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
σ_gr = .07 *λ0    # growth rate noise
Δ = 40                  # size increment added 
σ_M = sqrt(40)           # cell control noise
tmax = 90000


# initial conditions for all lineages
M0 = Δ
Mf0 = M0 + Δ
init = [M0,λ0,Mf0,0]


# ----------------------------------------------------------------------------------------------------------------
# RUN SIMULATIONS

# loop through tau values
τ_range = collect(1:20)
for τ in τ_range
    D = σ_gr^2/τ            #
    θ = [Δ,τ,λ0,D,σ_M]
    println("running τ = "*string(τ))
    prob,callback = build_model(θ,init,(0,tmax))
    @time sol = solve(prob,callback = callback)
    df = solver_output_to_dataframe(sol)

    cd(dirname(@__FILE__))
    folder = "./output/sims/"
    mkpath(folder)
    CSV.write(folder*"/sim_data_"*string(τ)*".csv",df)
 
end



