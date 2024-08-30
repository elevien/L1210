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
using Setfield
include("./GrowthTraceTools.jl")
using .GrowthTraceTools
cd(dirname(@__FILE__))


# ----------------------------------------------------------------------------------------------------------------
data = CSV.read("./output/data_processed.csv",DataFrames.DataFrame);
lineages = unique(data.lineage);

# based our simulated data of one lineage
data_df = data[data.lineage .== lineages[5],:]
data_df.time = data_df.time .- data_df.time[1];
times = data_df.time;

# ----------------------------------------------------------------------------------------------------------------
# RUN SIMULATIONS

# loop through tau values
τ_range = collect(1:2:10)
ϕ_range = collect(0.0001:0.3:0.99)
n_replicates = 100
θOU = GrowthTraceTools.θOU


sims_dfs = []
global lin = 0
for τ in τ_range
    for ϕ in ϕ_range

        # we change the parameters from the defaults to keep 
        θ = @set θOU.τ = τ 
        σ_NCS = sqrt(ϕ)*GrowthTraceTools.σGR # 
        θ = @set θOU.D = σ_NCS^2/τ
        θ = @set θOU.σDN  = sqrt((1-ϕ)*2*12/τ)*GrowthTraceTools.σGR
        println("running τ = "*string(τ)*" and ϕ = "*string(ϕ))


        # -----------------------------------------------------------
        # build model and run 
        prob,callback,names = GrowthTraceTools.build_model_OU(θ,GrowthTraceTools.initOU,times)
        sol = solve(prob,callback = callback);

        # -----------------------------------------------------------
        # put in dataframe 
        df = GrowthTraceTools.solver_output_to_dataframe(sol,names)
        df = df[df.position .< max(df.position...),:]
        df[:,:ϕ] = ones(length(df.time)) .*ϕ
        df[:,:τ] = ones(length(df.time)) .*τ
        
        df[:,:lineage] = ones(length(df.time)) .* lin 

        push!(sims_dfs,df)

        # increment lineage
        global lin = lin + 1
    end
end
sims_df = vcat(sims_dfs...)

#
cd(dirname(@__FILE__))
folder = "./output/sims"
mkpath(folder)
CSV.write(folder*"/sims.csv",sims_df)





