function gp_pipeline(data,model,dt_prediction = 0.5)
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

    println("Type of data: ", typeof(data))
    println("Is data empty? ", isempty(data))
    println("Number of rows: ", size(data, 1))
    println("Number of columns: ", size(data, 2))
    println("Columns in data: ", names(data))
    
    # Check if data is valid
    if size(data, 2) == 0
        error("Input DataFrame has no columns")
    end
    
    if !("lineage" in names(data))
        error("No 'lineage' column found in data. Available columns are: ", join(names(data), ", "))
    end

    println("First few rows of data:")
    println(first(data, 3))

    data[!, :lineage] = Int.(data[!, :lineage])  # Convert lineage to integers
    lineages = unique(data.lineage)
    println("Unique lineages: ", lineages)
    println("Type of lineage column: ", eltype(data.lineage))


    # # 2: GAUSSIAN PROCESSES  -----------------------------------------------------------------------------------------------

    # setup
    ops = Optim.Options(g_tol = 1e-4,iterations = 100,store_trace = true,show_trace = false);


    # run
    pred_dfs = []
    param_dfs = []
    for lin in lineages
        println("running lineage "*string(lin)*"--------------")

        df = data[data.lineage .== lin, :]
        df[:, :time] = df[:, :time] .- df[:, :time][1]
        df = df[df.cellcycle .=="i",:]
        df[:,:y] = df[:,:lnM_sum] .- df[:,:lnM_sum][1]
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
        
        # catch e
        #     println(" *** Error encountered with lineage "*string(lin)*": ", e)
        #     println(" *** Skipping lineage")
        # end
    end
    pred_df = DataFrame(vcat(pred_dfs...))
    param_df = vcat(param_dfs...)

    pred_df[:,:age] = vcat([d.time .- d.time[1] for d in groupby(pred_df,[:lineage,:position])]...);

    return pred_df, param_df
end

