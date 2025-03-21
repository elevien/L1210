# ########################################################################################
# SETUP
# ########################################################################################
# run with julia --project=. ./pipeline/run.jl

cd(dirname(@__FILE__))
include("./imports.jl")
include("./gp_pipeline.jl")


# preprocess exerimental data
# include("./preprocess.jl")


# ########################################################################################
# FIG 4
# ########################################################################################
# # ----------------------------------------------------------------------------------------
# # run GP on raw data
# data_source = "./../output/data_processed.csv"
# data = CSV.read(data_source,DataFrames.DataFrame);
# println(unique(data.length))
# data = data[(data.length .>= 9),:];
# model = GrowthTraceTools.Matern32Model()
# gp_pipeline(data,"./../output/gp/data",model)
# alert("Finished running GP on data")

#include("./fitting.jl")
#alert("Finished fitting GP params")


# # ----------------------------------------------------------------------------------------
# # generate simulations
#include("./fig4_sims.jl")
#alert("Finished fig4 simulations")


# # ----------------------------------------------------------------------------------------
# # run GP on simulated data
# sim_source = "./../output/sims_OUfit.csv"
# sim_data = CSV.read(sim_source,DataFrames.DataFrame);
# sim_data[:,:lnM_sum] = sim_data[:,:lnM_sum] .+ rand(Normal(0,0.001)) # add experimental noise
# sim_data = sim_data[sim_data.replicate .==1,:] # only run 1 replicate through the GP
# model = GrowthTraceTools.Matern32NoTrendModel()
# gp_pipeline(sim_data,"./../output/gp/sims",model)
# alert("Finished running GP on fig 4 data")

# ########################################################################################
# FIG 5
# ########################################################################################

# ----------------------------------------------------------------------------------------
# generate simulations
#include("./fig5_sims.jl")
#alert("**** Finished fig5 simulations ****")

# ----------------------------------------------------------------------------------------
# run GP on one replicate of sims
sim_source = "./../output/sims_models.csv"
sim_data = CSV.read(sim_source,DataFrames.DataFrame);
sim_data[:,:lnM_sum] = sim_data[:,:lnM_sum] .+ rand(Normal(0,0.001)) # add experimental noise
sim_data = sim_data[sim_data.replicate .== 1,:] # only run 1 replicate through the GP
model = GrowthTraceTools.Matern32NoTrendModel()
gp_pipeline(sim_data,"./../output/gp/sims_fig5",model)
alert("**** Finished running GP on fig 5 sims ****")


