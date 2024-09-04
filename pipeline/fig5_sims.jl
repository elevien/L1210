
# ----------------------------------------------------------------------------------------------------------------
# get raw data as template 
data = CSV.read("./../output/data_processed.csv",DataFrames.DataFrame);

# get GP output for fitting 
pred_dir = "./../output/gp/data" # CHECK 
data_gp = CSV.read(pred_dir*"/preds.csv",DataFrame);
data = data[data.length .>=9,:]
lineages = unique(data.lineage);

# ----------------------------------------------------------------------------------------------------------------
# setup 
τ_range = collect(2:0.1:15)

# ----------------------------------------------------------------------------------------------------------------
# 
nreps = 10
sims = []
global k = 1
for τ in τ_range

    println("-----------------------------------------------------")
    println("τ = "*string(τ))
    for lin in lineages
        println(" * replicating lineage "*string(lin)*"--------------")
        # -----------------------------------------------------------
        df_gp = data_gp[data_gp.lineage .==lin,:]
        df = data[data.lineage .==lin,:]
     

        for i in ProgressBar(1:nreps)

            # -----------------------------------------------------------
            # build OU model and run (all noise in diffusion coefficient)

            # setup params
            θ = GrowthTraceTools.θOU
            θ = @set θ.D = GrowthTraceTools.σGR.^2 ./ θ.τ
            θ = @set θ.σDN = 0.0 

            # run
            init = [θ.Δ,θ.λ0,2*θ.Δ]
            prob,callback,names = GrowthTraceTools.build_model_OU(θ,init,df.time)
            sol = solve(prob,EM(),callback = callback);

           
            # put in dataframe 
            sim = GrowthTraceTools.solver_output_to_dataframe(sol,names)
            sim = sim[sim.position .< max(sim.position...),:]
            sim[:,:lineage_original] = ones(length(sim.time)) .* lin
            sim[:,:replicate] = ones(length(sim.time)) .* i
            sim[:,:lineage] = ones(length(sim.time)) .* k
            sim[:,:model] = zeros(length(sim.time))
            push!(sims,sim)

            global k = k+1

            # -----------------------------------------------------------
            # build division noise model (all noise in jumps)

            # setup params
            θ = GrowthTraceTools.θOU
            θ = @set θ.D = 0.0
            θ = @set θ.σDN = GrowthTraceTools.σGR .* 10 /θ.τ
            init = [θ.Δ,θ.λ0,2*θ.Δ]

            # run 
            prob,callback,names = GrowthTraceTools.build_model_OU(θ,init,df.time)
            sol = solve(prob,EM(),callback = callback);

            # put in dataframe 
            sim = GrowthTraceTools.solver_output_to_dataframe(sol,names)
            sim = sim[sim.position .< max(sim.position...),:]
            sim[:,:lineage_original] = ones(length(sim.time)) .* lin
            sim[:,:replicate] = ones(length(sim.time)) .* i
            sim[:,:lineage] = ones(length(sim.time)) .* k
            sim[:,:model] = ones(length(sim.time))
            push!(sims,sim)
            
            global k = k+1

        end
    end
end
sims = vcat(sims...);

CSV.write("./../output/sims_models.csv",sims)