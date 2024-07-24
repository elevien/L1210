using Distributions
#using DifferentialEquations
using DataFrames,CSV,LinearAlgebra,StatsBase



###############################################################################
# OU Model
###############################################################################

function ou_make_cell(init,θ,trend)

    lambda_avg,tau,D,Delta,σ_ex,σ_ex_out,dt = values(θ)

    lambda = [init.lambda]
    M = [init.M]
    times = [init.time]

    l0 = lambda_avg + rand(Normal(0,σ_ex))

    # we now use the initial mass to determine the cell-cycle progression rate
    c = 1/lambda_avg
    f = 1/log((1 + Delta/M[1])^c)
    A_tau = 0.1
    T = 1/f + rand(Normal(0,A_tau))

    # simulate cell
    while times[end]-times[1]<T
        lambda_new = lambda[end] + (l0-lambda[end])/tau*dt + sqrt(2*D*dt)*rand(Normal())
        age = (times[end]+dt-times[1])
        M_new = M[end] + (lambda[end] + trend.(age ./T))*M[end]*dt

        push!(M,M_new)
        push!(lambda,lambda_new)
        push!(times,times[end] + dt)
    end

    ages_normed = (times .- times[1]) ./ (times[end] .- times[1])
    df = DataFrame([:time => times, :M => M,:flucs => lambda,:trend => trend.(ages_normed),:lambda=>lambda .+ trend.(ages_normed)])
    return df
end

function ou_make_lineage(init,θ,trend,L)
    cells = []
    sum = 0
    for k in 1:L
        cell = ou_make_cell(init,θ,trend)
        init = (time = cell.time[end]+θ.dt,M = cell.M[end]/2,lambda = cell.flucs[end] + rand(Normal(0,θ.σ_ex_out)))
        cell[:,:position] = ones(length(cell.time))*k
        cell[:,:age] = cell.time .- cell.time[1]
        cell[:,:age_normed] = cell[:,:age] ./ (cell.time[end] .- cell.time[1])
        cell[:,:lnM_sum] = log.(cell[:,:M]) .-log.(cell[1,:M]) .+ sum
        sum = cell[end,:lnM_sum]
        push!(cells,cell)
    end
    return vcat(cells...)
end


###############################################################################
# Mechanistic model
###############################################################################

function lm_make_cell(init,θ)

    α_R,α_Q,α_C,γ,C_crit,dt = values(θ)

    Q = [init.Q]
    R = [init.R]
    C = [init.C]
    times = [init.time]

    C_crit_cell = C_crit #+ rand(Normal(0,1000))
    while C[end]<C_crit_cell
        R_new = R[end] + α_R*R[end]*dt + sqrt(α_R*R[end]*dt)*rand(Normal())
        Q_new = Q[end] + α_Q*R[end]*dt + sqrt(α_Q*R[end]*dt)*rand(Normal())
        C_new = C[end] + α_C*R[end]*dt + sqrt(α_C*R[end]*dt)*rand(Normal())

        push!(R,R_new)
        push!(Q,Q_new)
        push!(C,C_new)
        push!(times,times[end] + dt)
    end

    S = Q .+ R .+ C
    df = DataFrame([:time => times, :Q => Q,:C=>C,:R=>R, :S=>S ])
    return df
end

function lm_make_lineage(init,θ,L)
    cells = []
    sum = 0
    for k in 1:L
        cell = lm_make_cell(init,θ)
        init = (time = cell.time[end],Q = cell.Q[end]/2*(1+θ.γ) + cell.C[end]/2,R = cell.R[end]/2*(1-θ.γ),C =0.0 )
        cell[:,:position] = ones(length(cell.time))*k
        cell[:,:age] = cell.time .- cell.time[1]
        cell[:,:age_normed] = cell[:,:age] ./ (cell.time[end] .- cell.time[1])
        cell[:,:lnM_sum] = log.(cell[:,:S]) .-log.(cell[1,:S]) .+ sum
        sum = cell[end,:lnM_sum]
        push!(cells,cell)
    end
    return vcat(cells...)
end
