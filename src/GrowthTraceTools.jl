module GrowthTraceTools
    using DataFrames
    using LinearAlgebra
    using StatsBase
    using Optim
    using ParameterHandling
    using DifferentialEquations
    using Distributions
    using Base.Threads

    abstract type GpDecomp end
    abstract type GpTerm end
    #export solver_output_to_dataframe,build_model_OU,θOU,initOU


    function map_lineage_col!(d1,d2,name)
        d2[:,name] = vcat(
        [
            d1[d1.lineage .==l,name][1]*ones(
                length(
                    d2[d2.lineage .==l,:time]
                    )
                ) 
        for l in unique(d2.lineage)
        ]
        ...);
    end

    include("./GP.jl")
    include("./gp_models/gp_models.jl")
    include("./Simulations.jl")
    include("./Theory.jl")
    include("./simulation_models/OU.jl")
    include("./processing_functions.jl")
end


