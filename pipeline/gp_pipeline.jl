
function gp_pipeline(data,out_dir,model,dt_prediction = 0.1)
    # 1: LOAD PROCESSED DATA -----------------------------------------------------------------------------------------------
    #data_src = "/Users/elevien/Dropbox (Dartmouth College)/RESEARCH/L1210_growth_rate_fluctuations/experimental_data/processed_data/traceupdated.csv"
    """
    Input dataframe must have columns: 
        - lineage (unique integer identifying each lineage)
        - length (number of cells in lineage)
        - lnM_sum (log mass)
        - time (time of mass measurements)
        - position (position of cell in lineage beginning at 1)
        - age_normed (the normalized age of each cell)
    """


    lineages = unique(data.lineage)


    # # 2: GAUSSIAN PROCESSES  -----------------------------------------------------------------------------------------------

    # setup
    ops = Optim.Options(g_tol = 1e-4,iterations = 100,store_trace = true,show_trace = false);


    # run
    pred_dfs = []
    param_dfs = []
    for lin in lineages
        println("running lineage "*string(lin)*"--------------")
        try
            # some new fields are needed in the data frames
            df = data[data.lineage .== lin, :]
            df[:, :time] = df[:, :time] .- df[:, :time][1]
    
            GrowthTraceTools.get_gen_times!(df)
    
            # now we make a uniform grid specifying the points at which to evaluate the GP
            pred_df = GrowthTraceTools.uniform_prediction_array(df, dt_prediction)
    
            # indices to use for optimizing hyperparameters
            opt_inds = 1:5:min(4000, length(df.time))
    
            # indices to use for interpolation
            obs_inds = 1:1:length(df.time)
    
            # run GP prediction 
            alg = NelderMead()
            @time d, results = GrowthTraceTools.gp_run(df, model, opt_inds, obs_inds, pred_df, ops, alg)
            params = Optim.minimizer(results)
    
            # push
            d[:, :lineage] = ones(length(d.time)) .* lin
            push!(pred_dfs, d)
            println(params)
            println("")
            push!(param_dfs, DataFrame(params = params))
            
        catch e
            println(" *** Error encountered with lineage "*string(lin)*": ", e)
            println(" *** Skipping lineage")
        end
    end
    pred_df = vcat(pred_dfs...)
    param_df = vcat(param_dfs...)

    pred_df[:,:age] = vcat([d.time .- d.time[1] for d in groupby(pred_df,[:lineage,:position])]...);

    mkpath(out_dir)
    CSV.write(out_dir*"/preds.csv",pred_df)
    CSV.write(out_dir*"/params.csv",param_df)
end


# # 2: SIMPLE SMOOTHING  -----------------------------------
# OUTDATED -- NEEDS UPDATED
# smoothing

# if smooth==true
#     γ = 2
#     function kernel(x,y)
#         exp(-(x-y)^2/γ^2)
#     end

#     # -----------------------------------------
#     # smoothing
#     # this is done by cells

#     function smooth_cell(time,y)

#         age = time .- time[1]
#         K = hcat([[kernel(time[i],time[j]) for i in 1:length(time)] for j in 1:length(time)]...);
#         H = K./ sum(K,dims=1)' ;
#         y_smoothed = H*y
#         gr = diff(H*y)
#         return (
#             age_normed = age ./ age[end],
#             gr = vcat(gr,gr[end]),
#             y_smoothed = y_smoothed
#             )
#     end

#     @time df_all_cells = groupby(data, [:lineage,:cell])
#     @time df_all_smoothed = transform(df_all_cells,[:time,:lnM] => smooth_cell => AsTable);



#     # ------------------------------------------
#     # detrending
#     # this is done by lineages
#     γ = .1
#     function kernel_trend(x,y)
#         exp(-abs(x-y)^2/γ^2)
#     end

#     function detrend_cell(age_normed,gr)
#         K_trend = hcat([[kernel_trend(age_normed[i],age_normed[j]) for i in 1:length(age_normed)] for j in 1:length(age_normed)]...)
#         gr_trend = ones(length(age_normed));
#         gr_dt = ones(length(age_normed));
#         gr_trend  =  (K_trend*gr ./ sum(K_trend,dims=1)')[:,1];
#         gr_dt = gr .- gr_trend ;
#         return (gr_trend = gr_trend,gr_dt= gr_dt)
#     end

#     @time df_all_lineages = groupby(df_all_smoothed, :lineage)
#     @time df_all_smoothed = transform(df_all_lineages,[:age_normed,:gr] => detrend_cell => AsTable);

#     CSV.write("./output/data_smoothed.csv", df_all_smoothed)
# end