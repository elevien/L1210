
# formulas for variance of time integrated process
Avar_theory_OU = (t,D,γ) -> D ./ γ .^2  .*(t .- (1 .- exp.( -2.0 * γ .* t)) ./ (2 * γ)) 
Avar_theory_DN = (t,v,γ) -> v .* ( 1.0 .- exp.( -t .* γ) ).^2 ./ γ^2

# formulas for variance of time integrated process
#s_λbar_theory_OU = (D,γ)