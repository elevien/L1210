
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
