
function simulate_ou(θ,init,time)
    sol = [vcat(init,[0,time[1]])] # last element is position
    age = 0. 
    for t in 2:length(time)
        M,λ,Mf,pos = sol[end]
        dt = (time[t] - time[t-1])
        Mnew = M + λ*M * dt
        age = age + dt

        # Exact OU process distribution
        λnew = θ.λ0 + (λ - θ.λ0) * exp(-dt/θ.τ) + sqrt(θ.D * θ.τ * (1 - exp(-2*dt/θ.τ))) * randn()
        if Mnew >= Mf
            Mnew = Mnew/2
            λnew =rand(Normal(λ,θ.σDN))
            Mfnew = Mnew + rand(Normal(θ.Δ,θ.σM))
            pos = pos + 1
            age = 0.
        else
            Mfnew = Mf
        end
        push!(sol,[Mnew,λnew,Mfnew,pos,time[t]])
    end
    sol = hcat(sol...)

    names = [:M,:λ,:Mf,:position,:time]
    df = DataFrame(sol',names)

    #get_positions!(df,1)
    positions = unique(df.position)
    df[:,:age] = vcat([df[df.position .==p,:time] .- df[df.position .==p,:time][1] for p in positions]...)
    df[:,:age_normed] = vcat([df[df.position .==p,:age]./df[df.position .==p,:age][end] for p in positions]...);
    df[:,:age_rounded] = round.(df[:,:age_normed],digits=1);
    df[:,:lnM] = log.(df[:,:M])


    # add summed log masses
    y =0.
    lnM_sums = []
    for p in unique(df.position)
        Z = df[df.position .==p,:].lnM
        Z = Z .- Z[1] .+ y #.+ mean(diff(Z))
        push!(lnM_sums,Z)
        y = Z[end]
    end
    lnM_sum = vcat(lnM_sums...);
    df[:,:lnM_sum] = lnM_sum
    df[:,:y] = df[:,:lnM_sum] .- df[:,:lnM_sum][1]
    df[:,:cellcycle] = fill("i",length(df[:,:time])) 

    df[:,"length"] = length(positions) .* ones(length(lnM_sum))
    return df
end

function get_positions!(df,threshold)
    p =0 
    position = [0]
    for i in 1:(length(df.time)-1)
        diff(df.M)[i] < -threshold ? p = p+1 : p=p
        push!(position,p)
        
    end
    df[:,:position] = position
end

function solver_output_to_dataframe(sol,names)
    u = hcat(sol.u...)';
    df = DataFrame(u,names)
    df[:,"time"] = sol.t;

    get_positions!(df,1)
    positions = unique(df.position)
    df[:,:age] = vcat([df[df.position .==p,:time] .- df[df.position .==p,:time][1] for p in positions]...)
    df[:,:age_normed] = vcat([df[df.position .==p,:age]./df[df.position .==p,:age][end] for p in positions]...);
    df[:,:age_rounded] = round.(df[:,:age_normed],digits=1);
    df[:,:lnM] = log.(df[:,:M])


    # add summed log masses
    y =0.
    lnM_sums = []
    for p in unique(df.position)
        Z = df[df.position .==p,:].lnM
        Z = Z .- Z[1] .+ y #.+ mean(diff(Z))
        push!(lnM_sums,Z)
        y = Z[end]
    end
    lnM_sum = vcat(lnM_sums...);
    df[:,:lnM_sum] = lnM_sum
    df[:,:y] = df[:,:lnM_sum] .- df[:,:lnM_sum][1]
    df[:,:cellcycle] = fill("i",length(df[:,:time])) 

    df[:,"length"] = length(positions) .* ones(length(lnM_sum))
    # need to add: length, lnM_sum, position
    return df
end   
