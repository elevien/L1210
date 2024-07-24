using DataFrames,CSV,LinearAlgebra,Optim,StatsBase
using SmoothingSplines
include("GPTools.jl")



################################################################################
cd(dirname(@__FILE__))

# get data
data = CSV.read(pwd()*"/../../experimental_data/processed_data/traceupdated.csv",DataFrame);
data = data[(data.length .> 20),:];
#data = data[data.position .< 10,:];
lineages = unique(data.lineage)

# setup
ops = Optim.Options(g_tol = 1e-5,iterations = 1000,store_trace = true,show_trace = false);
models = [LineageAnalysisTools.Matern32Model()]

# run
for lin in lineages

    # some new fields are needed in the data frames
    df = data[data.lineage.==lin,:]


    df = df[df.cellcycle.=="i",:]
    df[:,:y] = df[:,:lnM_sum] .- df[:,:lnM_sum][1]
    df[:,:time] = df[:,:time] .- df[:,:time][1]
    df[:,:gen_time]= vcat([ones(length(df[df.position .==p,:].time))*(df[df.position .==p,:].time[end]-df[df.position .==p,:].time[1])
        for p in unique(df.position)]...);

    # pred_df = df[1:6:length(df.time),[:time,:age,:M,:lineage,:position,:gen_time,:label,:cell]]
    pred_df = uniform_prediction_array(df,0.3)

    opt_inds = 1:6:min(3000,length(df.time))
    obs_inds = 1:1:length(df.time)

    cd(dirname(@__FILE__))
    folder = pwd()*"/output_7-23-24/lineage_$lin"
    mkpath(folder)
    @time output = gp_predict(df,model,opt_inds,obs_inds,pred_df)
    CVS.write(output[1])
end
