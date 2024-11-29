# --------------------------------------------------------------------
# This file fites the each lineage to the OU and division noise models based on the variance in time averaged growth rates


# setup and get data and gp output 
data = CSV.read("./../output/data_processed.csv",DataFrames.DataFrame);

# get GP output for fitting 
pred_dir = "./../output/gp/data" # CHECK 
data_gp = CSV.read(pred_dir*"/preds.csv",DataFrame);
data = data[data.length .>=9,:]
lineages = unique(data.lineage);

# add some columns used for fitting 
data_gp[:,:y] = vcat([d.matern32_flucs .- d.matern32_flucs[1] 
for d in groupby(data_gp,[:lineage,:position])]...);
data_gp[:,:yx] = vcat([d.matern32_flucs_x .- d.matern32_flucs_x[1] 
for d in groupby(data_gp,[:lineage,:position])]...);
data_gp[:,:ager] =  [round(x) for x in data_gp.age];


# ------------------------------------------------------------------- 
# fitting 

min_age = 6 # minimum age to use when fitting 
Xopts_OU = []
Xopts_DN = []

for k in eachindex(lineages) 

    # ------------------------
    # make var get dataframe to fit 
    d = data_gp[data_gp.lineage.== lineages[k],:]
    positions  = unique(d.position)
    df = combine(groupby(d,:ager),:y => var => :Avar)
    inds = df.ager .< min_age


    # functin to optimize 
    f = x -> sum(abs.(GrowthTraceTools.Avar_theory_OU(df.ager[inds],x[1],x[2]) .- df.Avar[inds]) .^2)
    x0 = [7.52e-5,0.08]
    results = Optim.optimize(f, x0, NelderMead())
    xopt_OU = Optim.minimizer(results)
    push!(Xopts_OU,xopt_OU)

    # ------------------------
    x0 = [7.52e-5,0.08]
    f = x -> sum(abs.(GrowthTraceTools.Avar_theory_DN(df.ager[inds],x[1],x[2]) .- df.Avar[inds]) .^2)
    results = Optim.optimize(f, x0, NelderMead())
    xopt_DN = Optim.minimizer(results)
    push!(Xopts_DN,xopt_DN)
 
end

Xopts_OU = hcat(Xopts_OU...)'
Xopts_DN = hcat(Xopts_DN...)'
Xopts = hcat(Xopts_OU,Xopts_DN)
fit_df = DataFrame(Xopts,[:D,:γ_OU,:σ2,:γ_DN])
fit_df[:,:lineage] = lineages
sort!(fit_df,:lineage)


raw_cells = combine(groupby(data,[:lineage,:position]),
:M=>(x->x[1])=>:M0,:M=>(x->x[end])=>:Mf,:time=>(x->x[end]-x[1])=>:gt)


raw_cells[:,:dM] = raw_cells.Mf .- raw_cells.M0
raw_cells[:,:λ] = log.(raw_cells.Mf ./raw_cells.M0 ) ./ raw_cells.gt
size_params = combine(
   groupby(
        raw_cells,:lineage
        ),
        :dM => mean => :Δ,
        :λ => mean => :λ0,
        :dM => std => :σM);

sort!(size_params,:lineage)

fit_df = hcat(size_params,fit_df,makeunique=true)

CSV.write("./../output/fitted_params.csv",fit_df)