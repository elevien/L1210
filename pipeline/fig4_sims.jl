

# ----------------------------------------------------------------------------------------------------------------
# get raw data as template 
data = CSV.read("./../output/data_processed.csv",DataFrames.DataFrame);

# get GP output for fitting 
pred_dir = "./../output/gp/data" # CHECK 
data_gp = CSV.read(pred_dir*"/preds.csv",DataFrame);
data = data[data.length .>=9,:]
lineages = unique(data.lineage);
lineages

# ----------------------------------------------------------------------------------------------------------------
# code used for fitting
function fit_ar(df)
    t = df[:,:time]
    dt = mean(diff(t))*1
    x = df[1:1:end,:matern32_flucs_x]
    X = hcat(ones(length(x)-1),x[1:end-1])
    y = x[2:end]
    b = X\y
    v = var(x)

    τ =  dt/(1-b[2])
    D = v/τ
    # I just use rough values for non-growth rate related parameters
    θ = (Δ = 40.0,σDN = 0.0,τ = τ,D = D,λ0 = 0.07,σM = sqrt(10.0)) 
    return θ
end


# ----------------------------------------------------------------------------------------------------------------
# 
nreps = 10
sims = []
global k = 1
for lin in lineages
    println("replicating lineage "*string(lin)*"--------------")
    # -----------------------------------------------------------
    df_gp = data_gp[data_gp.lineage .==lin,:]
    df = data[data.lineage .==lin,:]
    θ = fit_ar(df_gp)
    println("    θ = ",θ)
    init = [θ.Δ,θ.λ0,2*θ.Δ]

    for i in ProgressBar(1:nreps)
        # build model and run 
        prob,callback,names = GrowthTraceTools.build_model_OU(θ,init,df.time)
        sol = solve(prob,EM(),callback = callback);

        # -----------------------------------------------------------
        # put in dataframe 
        sim = GrowthTraceTools.solver_output_to_dataframe(sol,names)
        sim = sim[sim.position .< max(sim.position...),:]
        sim[:,:lineage_original] = ones(length(sim.time)) .* lin
        sim[:,:replicate] = ones(length(sim.time)) .* i
        sim[:,:lineage] = ones(length(sim.time)) .* k
        push!(sims,sim)
        global k = k+1
    end
end
sims = vcat(sims...);

CSV.write("./../output/sims_OUfit.csv",sims)