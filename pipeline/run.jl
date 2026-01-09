# ########################################################################################
# SETUP
# ########################################################################################
# run with julia --project=. ./pipeline/run.jl

cd(dirname(@__FILE__))
include("./imports.jl")
include("./gp_pipeline.jl")
using Random
Random.seed!(1234)

# make a folder ./../output/ if it doesn't exist
ispath("./../output/") || mkdir("./../output/")



# preprocess exerimental data
include("./preprocess.jl")

# The earlier figures are generated in the notebooks. 

# ########################################################################################
# FIG 4
# ########################################################################################
# # ----------------------------------------------------------------------------------------
# run GP on raw data
data_source = "./../output/data_processed.csv"
data = CSV.read(data_source,DataFrames.DataFrame);
println(unique(data.length))
data = data[(data.length .>= 5),:];
model = GrowthTraceTools.Matern32Model()
pred_df, param_df = gp_pipeline(data,model)
alert("Finished running GP on data")

# save to output/gp/data/preds.csv and output/gp/data/params.csv
# make folder if they don't exist 
ispath("./../output/gp/") || mkdir("./../output/gp/")
ispath("./../output/gp/data/") || mkdir("./../output/gp/data/")
# save files
CSV.write("./../output/gp/data/preds.csv",pred_df)
CSV.write("./../output/gp/data/params.csv",param_df)
alert("Saved GP outputs to ./../output/gp/data/")




# ########################################################################################
# FIG 5
# ########################################################################################

# ----------------------------------------------------------------------------------------
# # fit GP params
include("./fitting.jl")
alert("Finished fitting OU params")

# ----------------------------------------------------------------------------------------
# generate simulations 
println("**** Generating fig 5 simulations ****")
include("./fig5_sims.jl")
alert("**** Finished fig5 simulations ****")

# ----------------------------------------------------------------------------------------
# run GP on one replicate of sims
println("**** Running GP on fig 5 sims ****")
sim_source = "./../output/fig5_sims.csv"
sim_data = CSV.read(sim_source,DataFrames.DataFrame);
sim_data[:,:lnM_sum] = sim_data[:,:lnM_sum] .+ rand(Normal(0,0.001)) # add experimental noise
sim_data = sim_data[sim_data.replicate .< 20,:] 
model = GrowthTraceTools.Matern32NoTrendModel()
pred_df,param_df = gp_pipeline(sim_data,model)

function map_lineage_col!(d1,d2,name)
    d2[:,name] = vcat(
    [
        d1[d1.lineage .==l,name][1]*ones(
            length(
                d2[d2.lineage .==l,:time]
                )
            ) 
    for l in unique(d2.lineage)
    ]
    ...);
end

map_lineage_col!(sim_data,pred_df,:lineage_original);
map_lineage_col!(sim_data,pred_df,:D);
map_lineage_col!(sim_data,pred_df,:τ);

out_file = "./../output/fig5_sims"
CSV.write(out_file*"_gp_preds.csv",pred_df)
CSV.write(out_file*"_gp_params.csv",param_df)
alert("**** Finished running GP on fig 5 sims ****")


