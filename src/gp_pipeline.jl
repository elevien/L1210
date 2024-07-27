import Pkg;
Pkg.activate("./../")

include("L1210.jl")
include("./GP.jl")
include("./LineageAnalysisTools.jl")
using .L1210
println(Pkg.status())
using DataFrames,CSV,LinearAlgebra,Optim,StatsBase,PythonPlot


smooth=false
FIG_PATH = "./figures"




 # 1: LOAD PROCESSED DATA -----------------------------------------------------------------------------------------------
 #data_src = "/Users/elevien/Dropbox (Dartmouth College)/RESEARCH/L1210_growth_rate_fluctuations/experimental_data/processed_data/traceupdated.csv"

 data = CSV.read("./output/data_processed.csv",DataFrames.DataFrame);
 data = data[(data.length .> 10),:];
 lineages = unique(data.lineage)
 data[:,:y] = data[:,:lnM_sum]


# make individual figure
k = 1
for lin in lineages
    p = FIG_PATH*"/$lin"

    mkpath(p)
    fig,ax = subplots(figsize=(5,2))
    data_lin = data[data.lineage .== lin,:]

    ax.plot(data_lin.time,data_lin.M)
    savefig(p*"/raw_trace.pdf")
end


# # 2: SIMPLE SMOOTHING  -----------------------------------------------------------------------------------------------
# smoothing

if smooth==true
    γ = 2
    function kernel(x,y)
        exp(-(x-y)^2/γ^2)
    end

    # -----------------------------------------
    # smoothing
    # this is done by cells

    function smooth_cell(time,y)

        age = time .- time[1]
        K = hcat([[kernel(time[i],time[j]) for i in 1:length(time)] for j in 1:length(time)]...);
        H = K./ sum(K,dims=1)' ;
        y_smoothed = H*y
        gr = diff(H*y)
        return (
            age_normed = age ./ age[end],
            gr = vcat(gr,gr[end]),
            y_smoothed = y_smoothed
            )
    end

    @time df_all_cells = groupby(data, [:lineage,:cell])
    @time df_all_smoothed = transform(df_all_cells,[:time,:lnM] => smooth_cell => AsTable);



    # ------------------------------------------
    # detrending
    # this is done by lineages
    γ = .1
    function kernel_trend(x,y)
        exp(-abs(x-y)^2/γ^2)
    end

    function detrend_cell(age_normed,gr)
        K_trend = hcat([[kernel_trend(age_normed[i],age_normed[j]) for i in 1:length(age_normed)] for j in 1:length(age_normed)]...)
        gr_trend = ones(length(age_normed));
        gr_dt = ones(length(age_normed));
        gr_trend  =  (K_trend*gr ./ sum(K_trend,dims=1)')[:,1];
        gr_dt = gr .- gr_trend ;
        return (gr_trend = gr_trend,gr_dt= gr_dt)
    end

    @time df_all_lineages = groupby(df_all_smoothed, :lineage)
    @time df_all_smoothed = transform(df_all_lineages,[:age_normed,:gr] => detrend_cell => AsTable);

    CSV.write("./output/data_smoothed.csv", df_all_smoothed)
end

# # 2: GAUSSIAN PROCESSES  -----------------------------------------------------------------------------------------------

# setup
ops = Optim.Options(g_tol = 1e-5,iterations = 1000,store_trace = true,show_trace = false);
model = LineageAnalysisTools.Matern32Model()

# run
for lin in lineages

    # some new fields are needed in the data frames
    df = data[data.lineage.==lin,:]


    df = df[df.cellcycle.=="i",:]
    df[:,:y] = df[:,:lnM_sum] .- df[:,:lnM_sum][1]
    df[:,:time] = df[:,:time] .- df[:,:time][1]
    df[:,:gen_time]= vcat([ones(length(df[df.position .==p,:].time))*(df[df.position .==p,:].time[end]-df[df.position .==p,:].time[1])
        for p in unique(df.position)]...);

    pred_df = uniform_prediction_array(df,0.1)

    opt_inds = 1:6:min(3000,length(df.time))
    obs_inds = 1:1:length(df.time)

    cd(dirname(@__FILE__))
    folder = "./output/output_7-24-24/lineage_$lin"
    mkpath(folder)
    @time output = gp_predict(df,model,opt_inds,obs_inds,pred_df)
    output[1][:,:age] = pred_df[:,:age]
    output[1][:,:age_normed] = pred_df[:,:age_normed]
    CSV.write(folder*"/preds.csv",output[1])
    CSV.write(folder*"/opt_params.csv",output[2])
end