# generate simulations of each lineage from OU model (no division noise) based on fits to simulated lineages
cd(dirname(@__FILE__))
include("./imports.jl")



# ----------------------------------------------------------------------------------------------------------------
# get GP output
pred_dir = "./../output/gp_8-30-24" # CHECK 
data_gp = CSV.read(pred_dir*"/preds.csv",DataFrame);
lineages = unique(data_gp.lineage);
lineages

# ----------------------------------------------------------------------------------------------------------------
# code used for fitting
function fit_ar(df)
    t = df[:,:time]
    dt = mean(diff(t))
    x = df[1:1:end,:matern32_flucs_x]
    X = hcat(ones(length(x)-1),x[1:end-1])
    y = x[2:end]
    b = X\y
    v = mean((X*b .- y) .^2)

    τ =  dt/(1-b[2])
    D = v/2/dt
    # I just use rough values for non-growth rate related parameters
    θ = (Δ = 60.0,σDN = 0.0,τ = τ,D = D,λ0 = 0.07,σM = sqrt(30.0)) 
    return θ
end


# ----------------------------------------------------------------------------------------------------------------
# 
nreps = 5
sims = []
for lin in lineages
    println("replicating lineage "*string(lin)*"--------------")
    # -----------------------------------------------------------
    df = data_gp[data_gp.lineage .==lin,:]
    θ = fit_ar(df)
    init = [θ.Δ,θ.λ0,2*θ.Δ,0]

    for i in 1:nreps
        println("      rep"*string(i))
        # build model and run 
        prob,callback,names = GrowthTraceTools.build_model_OU(θ,init,df.time)
        sol = solve(prob,callback = callback);

        # -----------------------------------------------------------
        # put in dataframe 
        sim = GrowthTraceTools.solver_output_to_dataframe(sol,names)
        sim = sim[sim.position .< max(sim.position...),:]
        sim[:,:lineage] = ones(length(sim.time)) .* lin
        sim[:,:replicate] = ones(length(sim.time)) .* i
        push!(sims,sim)
    end
end
sims = vcat(sims...);

CSV.write("./../output/sims_OUfit.csv",sims)