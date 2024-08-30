
########################################################################################
# Gaussian process decomposition interface
########################################################################################



# kernel inferace which takes a model and makes the matrices for posterior inference

"""  
produces the covariance matrix of the total process (summed over kernels)
"""
function K_sum(m::GpDecomp,X::Array{Float64,2},lntheta)
    T = terms(m)
    n = length(X[:,1])
    #Ks = zeros(n,n)
    Ks = hcat([[sum([K(m,t,X[i,:],X[j,:],lntheta) for t in terms(m)]) for i in 1:n] for j in 1:n]...)
    # for t in T
    #     for i in 1:n
    #         for j in i:n
    #             Ks[i,j] = Ks[i,j] + K(m,t,X[i,:],X[j,:],lntheta)
    #             Ks[j,i] = Ks[i,j]
    #         end
    #     end
    # end
    return Ks
end

function K_mat(m::GpDecomp,term::String,X::Array{Float64,2},lntheta)
    n = length(X[:,1])
    Ks = zeros(n,n)
    for i in 1:n
        for j in i:n
            Ks[i,j] = Ks[i,j] + K(m,term,X[i,:],X[j,:],lntheta)
            Ks[j,i] = Ks[i,j]
        end
    end
    return Ks
end

function K_mat(m::GpDecomp,term::String,X1::Array{Float64,2},X2::Array{Float64,2},lntheta)

    n1 = length(X1[:,1])
    n2 = length(X2[:,1])
    Ks = zeros(n1,n2)
    for i in 1:n1
        for j in 1:n2
            Ks[i,j] = Ks[i,j] + K(m,term,X1[i,:],X2[j,:],lntheta)
        end
    end
    return Ks
end

function ddK_mat(m::GpDecomp,term::String,X::Array{Float64,2},lntheta)
    n = length(X[:,1])
    Ks = zeros(n,n)
    for i in 1:n
        for j in i:n
            Ks[i,j] = Ks[i,j] + ddK(m,term,X[i,:],X[j,:],lntheta)
            Ks[j,i] = Ks[i,j]
        end
    end
    return Ks
end

function dK_mat(m::GpDecomp,term::String,X1::Array{Float64,2},X2::Array{Float64,2},lntheta)

    n1 = length(X1[:,1])
    n2 = length(X2[:,1])
    Ks = zeros(n1,n2)
    for i in 1:n1
        for j in 1:n2
            Ks[i,j] = Ks[i,j] + dK(m,term,X1[i,:],X2[j,:],lntheta)
        end
    end
    return Ks
end

# posterior inference
function target(m::GpDecomp,Xo::Array{Float64,2},y::Array{Float64,1},lntheta)
    K = K_sum(m,Xo,lntheta)
    C = cholesky(K)
    alpha = C.U \(C.U' \ y)
    return 0.5*sum(log.(diag(C.U))) + 0.5*y'*alpha
end


"""  
makes the predictions at points Xp based on observations at points Xo and labels y 
using parameters lntheta 
"""
function predict(m::GpDecomp,Xp::Array{Float64,2},Xo::Array{Float64,2},y::Array{Float64,1},lntheta)
    T = terms(m)
    K = K_sum(m,Xo,lntheta)
    C = cholesky(K)
    alpha = C.U \(transpose(C.U) \ y)

    Z = DataFrame()  # predictions

    for t in T[T .!= "noise"]
        Kp = K_mat(m,t,Xp,lntheta)
        Kop = K_mat(m,t,Xp,Xo,lntheta)
        mu = Kop*alpha
        v = diag(Kp - Kop*(C.U \(transpose(C.U) \transpose(Kop))))
        Z[:,name(m)*"_"*t]= mu
        Z[:,name(m)*"_"*t*"_var"] = v
        # derivative
        if t in deriv_terms(m)
            ddKp = ddK_mat(m,t,Xp,lntheta)
            dKop = dK_mat(m,t,Xp,Xo,lntheta)
            mu = dKop*alpha
            v = diag(ddKp - dKop*(C.U \(transpose(C.U) \transpose(dKop))))
            Z[:,name(m)*"_"*t*"_x"]= mu
            Z[:,name(m)*"_"*t*"_x_var"] = v
        end
    end
    return Z

end



################################################################################
# GP implementation 

"""
gp_predict requires columns
    - time
    - y (log mass)
    - position 
"""

function gp_predict(df,model,opt_inds,obs_inds,pred_df)


    # Remove linear trend before inference
    # Without this step we run into problems with the cholesky factorization

    slope = cov(df[:,:y],df[:,:time])/var(df[:,:time])
    inter = mean(df[:,:y]) - slope*mean(df[:,:time])
    y_pred = slope*df[:,:time] .+ inter
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
    lnth0 =  lntheta0(model)
    f = x -> target(model,X[opt_inds,:],y[opt_inds],x)
    @time results = Optim.optimize(f,lntheta0(model),alg;autodiff = :forward)

    # predict
    Xp= Array(pred_df[:,[:time,:age_normed,:position]])
    Tb = push!([min(pred_df[pred_df.position .== p,:time]...) for p in unique(pred_df.position)],pred_df[end,:time])
    Xp = hcat(Xp,hcat([Tb for k in 1:length(Xp[:,1])]...)',pred_df.position);


    lnth_opt = results.minimizer
    @time Z = predict(model,Xp,X[obs_inds,:],y[obs_inds],lnth_opt);
    Z[:,:time] = Xp[:,1];
    Z[:,:position] = pred_df[:,:position];
    Z[:,:linear] = slope*Xp[:,1] .+ inter

    return Z,lnth_opt
end

################################################################################
# helper functions

function uniform_prediction_array(data,dt)
    # this function takes a dataframe and a time-step and produces a new
    # dataframe with samples at uniform times
    times = collect(data[1,:time]:dt:data[end,:time])
    positions = zeros(length(times))
    gen_times = zeros(length(times))
    cells = zeros(length(times))
    for k in eachindex(times)
        positions[k] = max(data[data.time .<= times[k],:position]...)
        gen_times[k] =  min(data[data.time .>= times[k],:gen_time]...)
        cells[k] =  min(data[data.time .>= times[k],:position]...)
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

