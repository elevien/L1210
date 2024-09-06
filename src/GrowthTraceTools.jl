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
    #export solver_output_to_dataframe,build_model_OU,θOU,initOU

    include("./GP.jl")
    include("./gp_models/gp_models.jl")
    include("./Simulations.jl")
    include("./simulation_models/OU.jl")
    include("./processing_functions.jl")
end


