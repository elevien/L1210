
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
τ_range = [2.0,4.0,8.0,10.0]
ϕ_range = [0,0.25,0.5,0.75,1.0]

# ----------------------------------------------------------------------------------------------------------------
# 
nreps = 5
sims = []
global k = 1
for τ in τ_range
    for ϕ in ϕ_range

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
                θ = @set θ.D = GrowthTraceTools.σGR.^2 ./ θ.τ .* ϕ 
                θ = @set θ.σDN = GrowthTraceTools.σGR .* 10 /θ.τ  .* sqrt((1-ϕ))

                # run
                init = [θ.Δ,θ.λ0 + rand(Normal(0,GrowthTraceTools.σGR)),2*θ.Δ]
                prob,callback,names = GrowthTraceTools.build_model_OU(θ,init,df.time)
                sol = solve(prob,EM(),callback = callback);

            
                # put in dataframe 
                sim = GrowthTraceTools.solver_output_to_dataframe(sol,names)
                sim = sim[sim.position .< max(sim.position...),:]
                sim[:,:lineage_original] = ones(length(sim.time)) .* lin
                sim[:,:replicate] = ones(length(sim.time)) .* i
                sim[:,:lineage] = ones(length(sim.time)) .* k
                sim[:,:ϕ] = ones(length(sim.time)) .* ϕ
                sim[:,:τ] = ones(length(sim.time)) .* τ
                push!(sims,sim)

                global k = k+1
            end
        end
    end
end
sims = vcat(sims...);

CSV.write("./../output/sims_models.csv",sims)