function solver_output_to_dataframe(sol,names)
    u = hcat(sol.u...)';
    df = DataFrame(u,names)
    df[:,"time"] = sol.t;
    positions = unique(df.position)
    df[:,:age] = vcat([df[df.position .==p,:time] .- df[df.position .==p,:time][1] for p in positions]...)
    df[:,:age_normed] = vcat([df[df.position .==p,:age]./df[df.position .==p,:age][end] for p in positions]...);
    df[:,:age_rounded] = round.(df[:,:age_normed],digits=1);
    return df
end   