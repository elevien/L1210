# add the columns needed for the Gaussian Process pipeline (PAPER_run_gp.jl)

cd(dirname(@__FILE__))
include("./imports.jl")

data = CSV.read("./../smr_data/traceupdated.csv",DataFrame)
lineages = unique(data.lineage)
lengths = vcat([
    length(
        unique(data[data.lineage .== lineages[i],:position])
        ) .*
    ones(
        length(
            data[data.lineage .== lineages[i],:time]
            )
        ) 
    for i in eachindex(lineages)]...)
data[:,:length] = lengths
ages_normed = combine(groupby(data,[:lineage,:position]),:age => (x -> x ./x[end]) => :age_normed).age_normed
data[:,:age_normed] = ages_normed
data[:,:y] =data[:,:lnM_sum] .- data[:,:lnM_sum][1]

data = data[data.cellcycle .== "i",:]

CSV.write("./../output/data_processed.csv",data)
