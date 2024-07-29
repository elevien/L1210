function build_model(θ,init,t_range)
    c_crit = θ[1]
    σM = θ[end]
    σCS = θ[2]
    p = θ[3:end]

    function F(du, u, p, t)
        # Model parameters.
        τ,λ0,D,σM = p
        du[1] =  u[2]*u[1]
        du[2] =  1/τ*(λ0 - u[2])
        du[3] = 0.0
        du[4] = 0

        return nothing
    end

    function G(du,u,p,t)
        τ,λ0,D,σM = p
        du[1] =  0.0
        du[2] =  sqrt(2*D)
        du[3] =  0.0
        du[4] =  0
    end

    function condition(u, t, integrator) # Event when condition(u,t,integrator) == 0
        u[3] - u[1]
    end

    function affect!(integrator)
        u = integrator.u
        u[1] = u[1]/2
        u[2] = u[2] + rand(Normal(0,σCS))
        u[3] = u[1] + rand(Normal(c_crit,σM))
        u[4] = u[4] + 1
        nothing
    end

    callback = ContinuousCallback(condition, affect!);
    prob = SDEProblem(F,G,init,t_range,p)
    return prob,callback
end
