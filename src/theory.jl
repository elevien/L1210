
# formulas for variance of time integrated process
Avar_theory_OU = (t,D,γ) -> D ./ γ .^2  .*(t .- (1 .- exp.( -2.0 * γ .* t)) ./ (2 * γ)) 
Avar_theory_DN = (t,v,γ) -> v .* ( 1.0 .- exp.( -t .* γ) ).^2 ./ γ^2 # v is varound at division 


# Inputs:
#   γ        = gamma
#   Td       = averaging window (T_d)
#   s2       = the FIRST term in Var(λ̄) decomposition
#   varbar   = total Var(λ̄)
#
# Outputs:
#   D, σz

function solve_D_σz_from_s2_and_varbar(γ, Td, s2, varbar)
    A = Td - (1/γ) * (1 - exp(-γ*Td))                # bracket term
    if A <= 0
        error("Invalid parameters: A = Td - (1/γ)(1-exp(-γ Td)) must be positive.")
    end

    D = (γ^2 * Td^2 / (2*A)) * s2

    Δ = varbar - s2                                   # remainder attributed to σz-term
    if Δ < 0
        error("Inconsistent inputs: varbar - s2 is negative (got $(Δ)).")
    end

    denom = abs(1 - exp(-γ*Td))
    σz = (γ * Td / denom) * sqrt(Δ)

    return D, σz
end
