
module LineageAnalysisTools
    using DataFrames,LinearAlgebra,StatsBase


    abstract type GpDecomp end # should this be here


    include("./GP.jl")
    include("./gp_models/gp_models.jl")
    #include("./moving_average.jl")



end

################################################################################


function gp_predict(df,model,opt_inds,obs_inds,pred_df)

    # Remove linear trend before inference
    # Without this step we run into problems with the cholesky factorization


    alpha = cov(df[:,:y],df[:,:time])/var(df[:,:time])
    beta = mean(df[:,:y]) - alpha*mean(df[:,:time])
    y_pred = alpha*df[:,:time] .+ beta
    y = df[:,:y] .- y_pred

    # make input array for optimization ----------------------------------------
    X = Array(df[:,[:time,:age_normed,:position]])

    #
    Tb = push!([min(df[df.position .== p,:time]...) for p in unique(df.position)],df[end,:time])

    X =  hcat(X,hcat([Tb for k in 1:length(X[:,1])]...)',df.position);



    # run
    ops = Optim.Options(iterations = 5,store_trace = true,show_trace = true);
    #alg = BFGS(;alphaguess = LineSearches.InitialStatic())
    #alg = ConjugateGradient()
    alg = NelderMead()
    lnth0 =  LineageAnalysisTools.lntheta0(model)
    f = x -> LineageAnalysisTools.target(model,X[opt_inds,:],y[opt_inds],x)
    @time results = Optim.optimize(f,LineageAnalysisTools.lntheta0(model),alg;autodiff = :forward)

    # predict
    Xp= Array(pred_df[:,[:time,:age_normed,:position]])
    Tb = push!([min(pred_df[pred_df.position .== p,:time]...) for p in unique(pred_df.position)],pred_df[end,:time])
    Xp = hcat(Xp,hcat([Tb for k in 1:length(Xp[:,1])]...)',pred_df.position);


    lnth_opt = results.minimizer
    @time Z = LineageAnalysisTools.predict(model,Xp,X[obs_inds,:],y[obs_inds],lnth_opt);
    Z[:,:time] = Xp[:,1];
    Z[:,:position] = pred_df[:,:position];
    Z[:,:linear] = alpha*Xp[:,1] .+ beta

    return Z,lnth_opt
end


function uniform_prediction_array(data,dt)
    # this function takes a dataframe and a time-step and produces a new
    # dataframe with samples at uniform times
    times = collect(data[1,:time]:dt:data[end,:time])
    positions = zeros(length(times))
    gen_times = zeros(length(times))
    cells = zeros(length(times))
    for k in 1:length(times)
        positions[k] = min(data[data.time .> times[k],:position]...)
        gen_times[k] =  min(data[data.time .> times[k],:gen_time]...)
        cells[k] =  min(data[data.time .> times[k],:position]...)
    end
    age = vcat([LinRange(0,gen_times[positions .== j][1],length(positions[positions .== j])) for j in unique(positions)]...);
    age_normed = age ./ gen_times
    lineages = ones(length(times))*data[1,:lineage]
    labels = fill(data[1,:label],length(times))
    return DataFrame(hcat(times,age,age_normed,positions,lineages,cells,gen_times)
    ,[:time,:age,:age_normed,:position,:lineage,:cell,:gen_time]);
end


function get_gen_times!(df)
    df[:,:gen_time]= vcat([ones(
            length(df[df.position .==p,:].time)
                )*(df[df.position .==p,:].time[end]-df[df.position .==p,:].time[1])
        for p in unique(df.position)]...);
end

