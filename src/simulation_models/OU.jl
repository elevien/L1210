function build_model_OU(θOU,init,times)
    Δ, σDN,τ,D ,λ0, σM  = Tuple(θOU)
    p = [τ,λ0,D,σM]



    function F(du, u, p, t)
        # Model parameters.
        τ,λ0,D,σM = p
        du[1] =  u[2]*u[1]          # mass
        du[2] =  1/τ*(λ0 - u[2])    # growth rate
        du[3] = 0.0                 # division size (Mf)
        du[4] = 0                   # cell 

        return nothing
    end

    function G(du,u,p,t)
        τ,λ0,D,σM = p
        du[1] =  0.0                       # no direct noise in mass
        du[2] =  sqrt(2*D)                 # growth noise
        du[3] =  0.0
        du[4] =  0
    end

    function condition(u, t, integrator) # Event when condition(u,t,integrator) == 0
        u[3] - u[1] # divide when Mf<= M
    end

    function affect!(integrator)
        u = integrator.u
        u[1] = u[1]/2                       # perfect division 
        u[2] = u[2] + rand(Normal(0,σDN))   # cell-division perturbation of growth
        u[3] = u[1] + rand(Normal(Δ,σM))    # adder stratagy
        u[4] = u[4] + 1                     # increment cell count
        nothing
    end

    callback = ContinuousCallback(condition, affect!);
    prob = SDEProblem(F,G,init,(min(times...),max(times...)),p,saveat=times,adaptive =false,dt=10e-5)
    return prob,callback,[:M,:λ,:Mf,:position]
end

# ------------------------------------------------
# baseline parameters uesd for simuations
# parameters

σGR = sqrt(1.32e-7)     # defualt growth rate standard deviation 
τOU = 4                 # this comes from looking at autocorrelations of GP output

θOU = (
    Δ = 40.0,       
    σDN = 0.00,         # noise injected at division 
    τ = τOU,            
    D = σGR^2/τOU,      # note that if σDN >0 we will want to adjust this
    λ0 = 0.05,          # average growth rate       
    σM = sqrt(40.0)
    )

# default initial conditions
initOU = [θOU.Δ,θOU.λ0,2*θOU.Δ,0]