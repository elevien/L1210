# setup used for all scripts in the pipeline
import Pkg;
println(Pkg.status())

using CSV
using StatsBase
using DataFrames
using PythonPlot
using DifferentialEquations
using LinearAlgebra
using Distributions
using Setfield
using Tables
using Optim
using ProgressBars
using Alert
include("./../src/GrowthTraceTools.jl")
using .GrowthTraceTools

