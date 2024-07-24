
########################################################################################
# Gaussian process decomposition interface
########################################################################################

abstract type GpDecomp end

# kernel infercae
function K_sum(m::GpDecomp,X::Array{Float64,2},lntheta)
    T = terms(m)
    n = length(X[:,1])
    #Ks = zeros(n,n)
    Ks = hcat([[sum([K(m,t,X[i,:],X[j,:],lntheta) for t in terms(m)]) for i in 1:n] for j in 1:n]...)
    # for t in T
    #     for i in 1:n
    #         for j in i:n
    #             Ks[i,j] = Ks[i,j] + K(m,t,X[i,:],X[j,:],lntheta)
    #             Ks[j,i] = Ks[i,j]
    #         end
    #     end
    # end
    return Ks
end

function K_mat(m::GpDecomp,term::String,X::Array{Float64,2},lntheta)
    n = length(X[:,1])
    Ks = zeros(n,n)
    for i in 1:n
        for j in i:n
            Ks[i,j] = Ks[i,j] + K(m,term,X[i,:],X[j,:],lntheta)
            Ks[j,i] = Ks[i,j]
        end
    end
    return Ks
end

function K_mat(m::GpDecomp,term::String,X1::Array{Float64,2},X2::Array{Float64,2},lntheta)

    n1 = length(X1[:,1])
    n2 = length(X2[:,1])
    Ks = zeros(n1,n2)
    for i in 1:n1
        for j in 1:n2
            Ks[i,j] = Ks[i,j] + K(m,term,X1[i,:],X2[j,:],lntheta)
        end
    end
    return Ks
end

function ddK_mat(m::GpDecomp,term::String,X::Array{Float64,2},lntheta)
    n = length(X[:,1])
    Ks = zeros(n,n)
    for i in 1:n
        for j in i:n
            Ks[i,j] = Ks[i,j] + ddK(m,term,X[i,:],X[j,:],lntheta)
            Ks[j,i] = Ks[i,j]
        end
    end
    return Ks
end

function dK_mat(m::GpDecomp,term::String,X1::Array{Float64,2},X2::Array{Float64,2},lntheta)

    n1 = length(X1[:,1])
    n2 = length(X2[:,1])
    Ks = zeros(n1,n2)
    for i in 1:n1
        for j in 1:n2
            Ks[i,j] = Ks[i,j] + dK(m,term,X1[i,:],X2[j,:],lntheta)
        end
    end
    return Ks
end

# posterior inference
function target(m::GpDecomp,Xo::Array{Float64,2},y::Array{Float64,1},lntheta)
    K = K_sum(m,Xo,lntheta)
    C = cholesky(K)
    alpha = C.U \(C.U' \ y)
    return 0.5*sum(log.(diag(C.U))) + 0.5*y'*alpha
end


#function predict_term (m::GpDecomp,Xo::Array{Float64,2},Xp::Array{Float64,2},y::Array{Float64,1},lntheta,)


function predict(m::GpDecomp,Xp::Array{Float64,2},Xo::Array{Float64,2},y::Array{Float64,1},lntheta)
    T = terms(m)
    K = K_sum(m,Xo,lntheta)
    C = cholesky(K)
    alpha = C.U \(transpose(C.U) \ y)

    Z = DataFrame()  # predictions

    for t in T[T .!= "noise"]
        Kp = K_mat(m,t,Xp,lntheta)
        Kop = K_mat(m,t,Xp,Xo,lntheta)
        mu = Kop*alpha
        v = diag(Kp - Kop*(C.U \(transpose(C.U) \transpose(Kop))))
        Z[:,name(m)*"_"*t]= mu
        Z[:,name(m)*"_"*t*"_var"] = v
        # derivative
        if t in deriv_terms(m)
            ddKp = ddK_mat(m,t,Xp,lntheta)
            dKop = dK_mat(m,t,Xp,Xo,lntheta)
            mu = dKop*alpha
            v = diag(ddKp - dKop*(C.U \(transpose(C.U) \transpose(dKop))))
            Z[:,name(m)*"_"*t*"_x"]= mu
            Z[:,name(m)*"_"*t*"_x_var"] = v
        end
    end
    return Z

end
