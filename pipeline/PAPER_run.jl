# runs the GP pipeline on input data (simulated or real)
cd(dirname(@__FILE__))
include("./imports.jl")
include("./gp_pipeline.jl")

# ------------------------------------------------------------------------------------
# run on raw data
# data_source = "./../output/data_processed.csv"
# data = CSV.read(data_source,DataFrames.DataFrame);
# println(unique(data.length))
# data = data[(data.length .>= 9),:];
# model = GrowthTraceTools.Matern32Model()
# gp_pipeline(data,"./../output/gp/data",model)

# ------------------------------------------------------------------------------------
# run on simulated data
sim_source = "./../output/sims_OUfit.csv"
sim_data = CSV.read(sim_source,DataFrames.DataFrame);
sim_data[:,:lnM_sum] = sim_data[:,:lnM_sum] .+ rand(Normal(0,0.001)) # add experimental noise
sim_data = sim_data[sim_data.replicate .==1,:] # only run 1 replicate through the GP
sim_data.lineage = sim_data.lineage_original
model = GrowthTraceTools.Matern32NoTrendModel()
gp_pipeline(sim_data,"./../output/gp/sims",model)



