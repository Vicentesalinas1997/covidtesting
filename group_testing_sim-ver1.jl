using Distributions, Plots, LinearAlgebra

function get_inc()
    μ = 1.621 #μ_l= 1.504; μ_u = 1.755  # parameter incubation period lognormal distribution
    σ = 0.418 #σ_l = 0.271; σ_u = 0.542 # parameter incubation period lognormal distribution
    d_inc  = LogNormal(μ,σ) # incubation period distribution
    t_inc  = max(2,round(rand(d_inc))) # incubation period sample
    ti_inf = min(t_inc-1,round(rand(Uniform(max(t_inc-3,t_inc/3),t_inc)))) # beginning of infeccious period  -sample
    te_inf = t_inc + round(rand(Uniform(14,28))) # ending of infeccious period - sample
    t_rec  = te_inf                     # time until recovery
    return t_inc, ti_inf, te_inf, t_rec
end

function cubic_spline(x, y)
    # check x are decreasing values
    if sum(x[1:end-1] .- x[2:end] .>= 0) > 0
        error("vector x no es estrictamente creciente")
    end
    # forces second derivatives to be 0 at both ends
    n = length(x) - 1
    idx = [1 / (x[i] - x[i-1]) for i in 2:(n+1)]
    dy = [y[i] - y[i-1] for i in 2:(n+1)]
    Ak = zeros(n + 1, n + 1)
    bk = zeros(n + 1)
    # first row of A and b
    Ak[1, 1] = 2 * idx[1]
    Ak[1, 2] = idx[1]
    bk[1] = 3 * dy[1] * idx[1]^2
    # intermediate rows of A and b
    for i in 2:n
        # for first derivate equal
        Ak[i, i - 1] = idx[i - 1]
        Ak[i, i] = 2 * (idx[i - 1] + idx[i])
        Ak[i, i + 1] = idx[i]
        bk[i] = 3 * (dy[i - 1] * idx[i - 1]^2 + dy[i] * idx[i]^2)
        # for second derivate equal
        #Ak[i, i - 1] = idx[i - 1]
        #Ak[i, i] = 2 * idx[i - 1]
        #bk[i] = 3 * dy[i - 1] * idx[i - 1]^2
    end
    # last row of A and b
    Ak[n + 1, n] = idx[n]
    Ak[n + 1, n + 1] = 2 * idx[n]
    bk[n + 1] = 3 * dy[n] * idx[n]^2
    #println("Ak = ", Ak)
    # compute k
    k = Ak \ bk
    #println("k = ", k)
    # compute a y b
    a = [k[i-1] * (x[i] - x[i-1]) - (y[i] - y[i-1]) for i in 2:(n+1)]
    b = [-k[i] * (x[i] - x[i-1]) + (y[i] - y[i-1]) for i in 2:(n+1)]
    return a, b
end

function cubic_spline_interpolate(x_int, x, y, a = nothing, b = nothing)
    if a == nothing
        a, b = cubic_spline(x, y)
    end
    if isa(x_int, Number)
        if (x_int < x[1]) | (x_int > x[end])
            #error("El valor del x a interpolar esta fuera del arreglo.")
            return 0.0
        end
        # encontrar el polinomio correspondiente del 1 a n
        #i = sum(x_int .<= x)
        i = 1
        while x_int > x[i + 1]
            i = i + 1
        end
        #println("i = ", i)
        t = (x_int - x[i]) / (x[i + 1] - x[i])
        #println("t = ", t)
        y_int = (1 - t) * y[i] + t * y[i + 1] + t * (1 - t) * ( (1 - t) * a[i] + t * b[i])
        return y_int
    else
        return [cubic_spline_interpolate(x_int[j], x, y, a, b) for j in 1:length(x_int)]
    end
end

function true_positive(t, t_inf, t_sin, t_rec, p_sinto = 0.8)
    # t : dia relativo con respecto a la expocision, tb puede ser un vector
    # vectores de valores conocidos en eje x e y
    x_vec = vcat([t_inf, (t_inf + t_sin) / 2, t_sin], [t_sin + (t_rec - t_sin) * i / 28 for i in 7:7:42])
    #x_vec = [-7, -3.5,   0,    7,   14,   21,   28,   35,   42]
    y_vec = [ 0, 0.63, 0.8, 0.79, 0.64, 0.43, 0.19, 0.08, 0.00] * (p_sinto / 0.8)
    y = cubic_spline_interpolate(t, x_vec, y_vec)

    # graficar
    if false #| true
        x_plot = [x_vec[1]:x_vec[end]]
        y_plot = cubic_spline_interpolate(x_plot, x_vec, y_vec)
        plot_interpolate(x_plot, y_plot, x_vec, y_vec)
    end

    return y
end

function simulation(N, f,G,T,α,β,γ,p_false_positive,R,test)
    pcr_M = zeros(N,f) # non-adaptive test matrix
    if test == true
        for i=1:f
            pcr_M[:,i] = [zeros(G*(i-1),1) ; ones(G,1) ; zeros(G*(f-i))]
        end
        pcr_M = Int.(pcr_M)
    end

    mQua = zeros(T)
    mInf = zeros(T)

    # Loop over replications
    for r = 1:R
        tEx =zeros(N) # time of exposure
        tIn =zeros(N) # time of infecciousness
        tSy =zeros(N) # time of development of symptoms
        tRe =zeros(N) # time of recovery

        As = rand(N).< γ # asymtomatic lottery
        Su = ones(N,T) # susceptible population
        Ex = zeros(N,T) # exposed population
        In = zeros(N,T) # infectious population
        Sy = zeros(N,T) # symptomatic population
        Re = zeros(N,T) # recovered population
        Qu = zeros(N,T) # quarantined population

        # Loop over days
        for t=1:T
            su  = Su[:,t]
            ex  = Ex[:,t]
            inf = In[:,t]
            sy  = Sy[:,t]
            re  = Re[:,t]
            qu  = Qu[:,t]

            # Contagion Dynamics
            p_inf = α * sum((1 .- qu) .* inf)/(sum(1 .- qu) + 0.001) + β # probability of infection (approx)ß
            new_inf = su .* (1 .- qu) .* (rand(N) .< p_inf)                # new infection
            for i in findall(new_inf.==1)                   # cicle over newly infected
                t_inc, ti_inf, te_inf, t_rec = get_inc()
                Su[i, Int(min(T,t+1)): T] .= 0                             # no longer susceptible
                Ex[i, Int(min(T,t+1)) : Int(min(T, t+ t_inc))] .= 1        # exposed individual
                In[i, Int(min(T,t+ti_inf)) : Int(min(T,t+te_inf))] .= 1    # infeccious individual
                Sy[i, Int(min(T,t+t_inc)) : Int(min(T,t+te_inf))] .= 1-As[i] # symptomatic individual
                Re[i, Int(min(T,t+t_rec)): T] .= 1                         # recovered individual
                Qu[i, Int(min(T,t+t_inc+2)): Int(min(T,max(t+15,t+t_rec+3)))] .= 1-As[i] #you must go to quarantine if you get symptoms
                tEx[i] = t + 1       # time of exposure
                tSy[i] = t + 1 + t_inc   # time of development of symptoms
                tIn[i] = t + 1 + ti_inf # time of infecciousness
                tRe[i] = t + 1 + t_rec   # time of recovery
            end

            # Testing
            pcr_v = pcr_M[:,mod(t-1,f)+1] # testing vectror
            cand = findall(pcr_v .== 1)
            if test == "random"
                cand = ((sum(Qu[:,1:t],dims=2)==0) .> 0)
                s_cand = Random.randperm(sum(cand))
                cand = findall(cand .==1)# sample G people at random
                pcr_v = cand[s_cand]
                pcr_v = pcr_v[1:Int(min(G,end))]
            end
            for i in cand
                if (su[i] == 1) & (rand() < p_false_positive)
                    Qu[i,min(T,t+1):min(T,t+15)] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
                end
                if su[i] == 0
                    if (qu[i] == 0) & (sum(Qu[i,1:t])==0)
                        p_positive = true_positive(t-tEx[i]+1, tIn[i]-tEx[i]+1, tSy[i]-tEx[i]+1, tRe[i]-tEx[i]+1)
                        if rand() < p_positive
                            Qu[i,Int(min(T,t+1)):Int(min(T,max(t+15,tRe[i]+3)))] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
                        end
                    end
                end
            end
        end
        mQua .+= dropdims(sum(Qu,dims=1),dims=1)
        mInf .+= dropdims(sum((1 .- Qu) .* In,dims=1),dims=1)
    end

    mQua = mQua/R
    mInf = mInf/R

    return mQua, mInf
end


N = 200             # Tamaño grupo de trabajadores
f = 10              # Frecuencia testeo (cada T+1 dias)
G = Int(round(N/f)) # tamaño grupo a testear - por ahora un numero entero
T = 365             # simulation horizon
α = 0.082   #calculado           # adjusted prob of contagion amoong coworkers
β = 0.007   #calculado           # external prob of contagion
#probabilidad de contagiarse desde los pacientes
γ = 0.30 #(antes 0.25)           # prop of asymtomatic pop
p_false_positive = 0 # probability of a false positive
R = 1000       # Replications

base_mQua, base_mInf = simulation(N, f,G,T,α,β,γ,p_false_positive,R, false)
ideal_mQua, ideal_mInf = simulation(N, 1,N,T,α,β,γ,p_false_positive,R, true)
f10_mQua, f10_mInf = simulation(N, f,G,T,α,β,γ,p_false_positive,R, true)
r10_mQua, r10_mInf = simulation(N, f,G,T,α,β,γ,p_false_positive,R, "random")

f1 =plot(1:T-1,base_mQua[1:end-1],label ="base", lw=3)
plot!(1:T-1,ideal_mQua[1:end-1],label ="ideal", lw=3)
plot!(1:T-1,f10_mQua[1:end-1],label ="10 fijo", lw=3)
plot!(1:T-1,r10_mQua[1:end-1],label ="10 rand", lw=3)
title!("Cuarentena")
xlabel!("Días")
ylabel!("# personas")

f2 =plot(1:T-1,base_mInf[1:end-1], label ="base", lw=3)
plot!(1:T-1,ideal_mInf[1:end-1], label ="ideal", lw=3)
plot!(1:T-1,f10_mInf[1:end-1], label ="10 fijo", lw=3)
plot!(1:T-1,r10_mInf[1:end-1], label ="10 rand", lw=3)
title!("Infección")
xlabel!("Días")
ylabel!("# personas")

plot(f1,f2, layout = (2,1))
