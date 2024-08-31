# add the columns needed for the Gaussian Process pipeline (PAPER_run_gp.jl)

cd(dirname(@__FILE__))
include("./imports.jl")

data = CSV.read("./../smr_data/traceupdated.csv",DataFrame)
lineages = unique(data.lineage)
lengths = vcat([
    length(
        unique(data_gp[data_gp.lineage .== lineages[i],:position])
        ) .*
    ones(
        length(
            data_gp[data_gp.lineage .== lineages[i],:time]
            )
        ) 
    for i in eachindex(lineages)]...)
data_gp[:,:length] = lengths
ages_normed = combine(groupby(data,[:lineage,:position]),:age => (x -> x ./x[end]) => :age_normed).age_normed
data[:,:age_normed] = ages_normed

CSV.write("./../output/data_processed.csv",data)
