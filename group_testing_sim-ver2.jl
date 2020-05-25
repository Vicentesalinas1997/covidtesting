using Distributions, Plots, LinearAlgebra, Random

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
        #error("vector x no es estrictamente creciente")
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
    x_vec = vcat([t_inf,(t_inf + t_sin) / 4,3*(t_inf + t_sin) / 4, t_sin], [t_sin + (t_rec - t_sin) * i / 28 for i in 7:7:42])
    #x_vec = [-7, -3.5,   0,    7,   14,   21,   28,   35,   42]
    y_vec = [ 0, 0.63, 0.8, 0.79, 0.64, 0.43, 0.19, 0.08, 0.00,0.00] * (p_sinto / 0.8)
    y = cubic_spline_interpolate(t, x_vec, y_vec)

    # graficar
    if false #true
        x_plot = [x_vec[1]:x_vec[end]]
        y_plot = cubic_spline_interpolate(x_plot, x_vec, y_vec)
        plot_interpolate(x_plot, y_plot, x_vec, y_vec)
    end

    return y
end



#Mejorar funcion de probabilidad de contagiar dependiento del tiempo t, el de icubación
#,infeccion y el de recuperacion. De momento se copia la funcion true positive cno un p no tan grande
function person_p(t,t_inf,t_sin,t_rec,p_sinto)
    # t : dia relativo con respecto a la expocision, tb puede ser un vector
    # vectores de valores conocidos en eje x e y
    x_vec = vcat([t_inf, (t_inf + t_sin) / 2, t_sin], [t_sin + (t_rec - t_sin) * i / 28 for i in 7:7:42])
    #x_vec = [-7, -3.5,   0,    7,   14,   21,   28,   35,   42]
    y_vec = [ 0, 0.63, 0.8, 0.79, 0.64, 0.43, 0.19, 0.08, 0.00] * (p_sinto/0.8)
    y = cubic_spline_interpolate(t, x_vec, y_vec)
    return y
end


#Hacer una funcion de probabilidad de contagio en la casa, dependiente de horas, n integrantes
# integrantes contagiados y tiempo (pues a futuro se puede detallar las etapas de los
#contagiados en la casa agregando los tiempos de estos).



# integrantes contagiados y tiempo (pues a futuro se puede detallar las etapas de los contagiados en la casa
#agregando los tiempos de estos).
function casa_p(n,n_inf,hr)
return 0.007
end




#Funcion que entrega la probabilidad de contagiarte dentro de tu grupo de trabajo
#cerrado
function group_p(G,vp_cont,hr,W)
	#G: indicatriz del grupo de trabajo
	#vp_cont: vector con las probabilidades de contagio
	#hr: horas de trabajo
	#W: Peligro de contagio de pacientes
	V=G.*vp_cont
	p=(1-prod(-V.+1))*((2/pi)*atan(2*hr)) #Es el complemento de la probabilidad
	#de no contagiarse multiplicado por un factor que pondera las horas, este desde
	#la hora toma valores del estilo ~0.9 y mas.
	if W=="none"
		return p
	end
	if W=="low"
		return 1-(1-p)*(1-0.01) #Esta probabilidades estan sin ninguna justificacion
		#podria depender de las horas de atencion a publico y medidas de seguridad
	end
	if W=="medium"
		return 1-(1-p)*(1-0.05)
	end
	if W=="high"
		return 1-(1-p)*(1-0.1)
	end
end


function simulation(N, f,G,T,β,γ,p_false_positive,R,test,S,prob,Group,Tra)


    mQua = zeros(T)
    mInf = zeros(T)
	mNTest=zeros(T)
	mNFpositive=zeros(T)
	mNQincubacion=0

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
		NTest=zeros(T)
		NFpositive=zeros(T)
		NQincubacion=0


        # Loop over days
        for t=1:T
            su  = Su[:,t]
            ex  = Ex[:,t]
            inf = In[:,t]
            sy  = Sy[:,t]
            re  = Re[:,t]
            qu  = Qu[:,t]
	    tra= Tra[:,t]

            # Contagion Dynamics
            #p_inf = α * sum((1 .- qu) .* inf)/(sum(1 .- qu) + 0.001) + β
	    # probability of infection (approx)ß
            #new_inf = su .* (1 .- qu) .* (rand(N) .< p_inf)
	    # new infection



	    new_inf=zeros(N) #Vector de los nuevos infectados
	    con=0 #Cuantos grupos hay para hacer poll testing
	    for a=1:S
		info=Group[Int((a-1)*(N+2)+1):Int((a)*(N+2))]
		g=Int.(info[1:N]) #Indicatriz del grupo
		hr=info[N+1] #Horas del grupo
		W=info[N+2] #Riesgo del grupo
		vec=zeros(N) #Vector auxiliar para las probabilidades de los infectados
		#del grupoo
		for i in findall(g.==1)
			if (inf[i]==1) & (qu[i]==0) & (re[i]==0) & (tra[i]==1)  #Condicion esta en el grupo, infectado
				#y en cuarentena
				#Recuperacion de los tiempos
				tinf=tIn[Int(i)]-tEx[Int(i)]
				tsin=tSy[i]-tEx[i]
				trec=tRe[i]-tEx[i]
				vec[i]=person_p(t,tinf,tsin,trec,prob)
				 #Se agregan las probabilidades de los enfermos
			end
		end
		PG=group_p(g,vec,hr,W) #Vector con probabilidades de contagio dentro del grupo
		r=rand(N)
		#Reemplazar el beta y 1.2 beta por una probabilidad de contagio en casa idividual.
		new_inf+=(tra).*(-qu.+1).*su.*g.*(r.<(1-(1-PG)*(1-β))) #Vector de infectados entre todo el hospital
		#new_inf+=qu.*su.*g.*(r.<(1.01*β)) #Sumar infectados en cuarentena (falsos positivos
		#que se contagian en casa,por lo que la probabilidad de casa deberia aumentar)




		#Para el poll test
		I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*(-qu.+1).*(tra).*collect(1:N)
		JJ=shuffle!(g.*I)
		III=setdiff(JJ,[0])
		if length(III)==0
			con+=0
		else
			ind=Int(round(length(III)/5))
			r=length(III)%5
			con+=(ind+(r!=0))
		end




	   end
           for i in findall(new_inf.==1)                   # cicle over newly infected
               	t_inc, ti_inf, te_inf, t_rec = get_inc()
             	Su[i, Int(min(T,t+1)): T] .= 0                             # no longer susceptible
      	        Ex[i, Int(min(T,t+1)) : Int(min(T, t+ t_inc))] .= 1        # exposed individual
               	In[i, Int(min(T,t+ti_inf)) : Int(min(T,t+te_inf))] .= 1    # infeccious individual
               	Sy[i, Int(min(T,t+t_inc)) : Int(min(T,t+te_inf))] .= (1-As[i]) # symptomatic individual
                Re[i, Int(min(T,t+t_rec)): T] .= 1                         # recovered individual
                Qu[i, Int(min(T,t+t_inc+2)): Int(min(T,max(t+15,t+t_rec+3)))] .= (1-As[i]) #you must go to quarantine if you get symptoms
               	tEx[i] = t + 1       # time of exposure
               	tSy[i] = t + 1 + t_inc   # time of development of symptoms
                tIn[i] = t + 1 + ti_inf # time of infecciousness
                tRe[i] = t + 1 + t_rec   # time of recovery
            end

            # Testing

		#Vector de funcionarios que testear
	    cand=[]

		if mod(t-1,f)==0 #Se testea cada f dias
	    		if test==true
				I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*(tra).*(-qu.+1).*collect(1:N)
				J=shuffle!(I)
				II=setdiff(J,[0])
				cand=Int.(II[1:min(G,length(II))])
	    		end
			if test=="poll"
				h=shuffle!(collect(1:con))[1:min(G,con)]
				b=sort(h)
				for a=1:S
					info=Group[Int((a-1)*(N+2)+1):Int((a)*(N+2))]
					g=Int.(info[1:N]) #Indicatriz del grupo
					I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*(-qu.+1).*collect(1:N)
					JJ=shuffle!(g.*I)
					III=setdiff(JJ,[0])
					if length(III)>=1
						ind=Int(floor(length(III)/5))
						r=length(III)%5
						lim=(ind+(r!=0))
						c=setdiff((b.<=lim).*(b.>=1).*b,[0])
						if length(c)>=1
							b.-lim
							for cc in c
								if cc-1==ind
									cand=Int.(III[5*ind+1:5*ind+r])
								else
									cand=Int.(III[5*(cc-1)+1:5*cc])
								end
								#hacer poll testing
								s=0
								pp=1
								for i in cand
									p_positive = true_positive(t-tEx[i]+1, tIn[i]-tEx[i]+1, tSy[i]-tEx[i]+1, tRe[i]-tEx[i]+1)
									s+=su[i]
									pp=pp*(1-(1-su[i])*p_positive)
								end
								PP=1-pp
								NTest[t]+=1
								if (s==length(cand)) & (rand() < (p_false_positive))
									for i in cand
										NTest[t]+=1
										if (su[i]==1) & (rand() < p_false_positive)
                    									Qu[i,min(T,t+2):min(T,t+16)] .= 1 #if someone tests (false) positive, quarantined for 2 weeks
											NFpositive[t]+=1
										end
									end
								else
									if (rand()<=PP) & (s<=(length(cand)-1))
										for i in cand
											NTest[t]+=1
                									if su[i] == 0
                    										if (qu[i] == 0) & (sum(Qu[i,1:t])==0)
                        										p_positive = true_positive(t-tEx[i]+1, tIn[i]-tEx[i]+1, tSy[i]-tEx[i]+1, tRe[i]-tEx[i]+1)
                        										if rand() < p_positive
                            											Qu[i,Int(min(T,t+1)):Int(min(T,max(t+15,tRe[i]+3)))] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
														if ex[i]==1 #Esta en incubacion
															NQincubacion+=1
														end
													end
                        									end
                    									else
												if (rand() < p_false_positive)
                    											Qu[i,min(T,t+2):min(T,t+16)] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
													NFpositive[t]+=1
												end
											end
										end
                							end
            							end
							end
						end
					end
				end
			end

	if test!="poll"
		for i in cand
			NTest[t]+=1
                	if (su[i] == 1) & (rand() < p_false_positive)
                    		Qu[i,min(T,t+2):min(T,t+16)] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
				NFpositive[t]+=1
                	end
                	if su[i] == 0
                    		if (qu[i] == 0) & (sum(Qu[i,1:t])==0)
                        		p_positive = true_positive(t-tEx[i]+1, tIn[i]-tEx[i]+1, tSy[i]-tEx[i]+1, tRe[i]-tEx[i]+1)
                        		if rand() < p_positive
                            			Qu[i,Int(min(T,t+1)):Int(min(T,max(t+15,tRe[i]+3)))] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
							if ex[i]==1 #Esta en incubacion
								NQincubacion+=1
							end
                        		end
                    		end
                	end
            	end

        end
end

end
        mQua .+= dropdims(sum((Tra).*Qu,dims=1),dims=1)
        mInf .+= dropdims(sum(Tra.*(-Qu.+1).*In,dims=1),dims=1)
		mNFpositive .+= NFpositive
		mNQincubacion+=NQincubacion
		mNTest .+= NTest
    end

    mQua = mQua/R
    mInf = mInf/R
	mNFpositive=mNFpositive/R
	mNQincubacion=mNQincubacion/R
	mNTest=mNTest/R

    return mQua, mInf, mNFpositive, mNQincubacion, mNTest
end


N = 200             # Tamaño grupo de trabajadores
f = 10              # Frecuencia testeo (cada T+1 dias)
G = Int(round(N/f)) # tamaño grupo a testear - por ahora un numero entero
T = 147             # simulation horizon
α = 0.082   #calculado           # adjusted prob of contagion amoong coworkers
β = 0.01   #calculado           # external prob of contagion
#probabilidad de contagiarse desde los pacientes
γ = 0.30 #(antes 0.25)           # prop of asymtomatic pop
p_false_positive = 0.01 # probability of a false positive
R = 500    # Replications
S=10 #Cantidad de grupos
p=0.2 #Prob peak de infeccion

    #Creacion de un grupo de ejemplo
    for a=1:S
	global Group
	if a==1
		g1=[zeros(Int((a-1)*(N/S)))' ones(Int(N/S))' zeros(Int(N-a*(N/S)))']
		g2=8
		g3="none"
		Group=[g1 g2 g3]
	else

		g1=[zeros(Int((a-1)*(N/S)))' ones(Int(N/S))' zeros(Int(N-a*(N/S)))']
		g2=8
		g3="none"
		g=[g1 g2 g3]
		Group=[Group g]
	end
    end




Trab=ones(N,T)
Trab2=zeros(N,T)
for i in 1:Int(round(T/7))
	for j in 1:N
		if (mod(j,2)==1) & (mod(i,2)==1)
			Trab2[j,(i-1)*7+1:i*7]=ones(7)
		end
		if (mod(j,2)==0) & (mod(i,2)==0)
			Trab2[j,(i-1)*7+1:i*7]=ones(7)
		end
	end
end



base_mQua, base_mInf, base_mNFp, base_mNQinc, base_T= simulation(N, f,G,T,β,γ,p_false_positive,R, false,S,p,Group,Trab2)
ideal_mQua, ideal_mInf, ideal_mNFp, ideal_mNQinc, ideal_T = simulation(N, 1,N,T,β,γ,p_false_positive,R, true,S,p,Group,Trab2)
f10_mQua, f10_mInf, f10_mNFp, f10_mNQinc, f10_T = simulation(N, f,G,T,β,γ,p_false_positive,R, true,S,p,Group,Trab2)
f102_mQua, f102_mInf,f102_mNFp,f102_mNQinc, f102_T = simulation(N, f,N,T,β,γ,p_false_positive,R, true,S,p,Group,Trab2)
f1G_mQua, f1G_mInf,f1G_mNFp,f1G_mNQinc, f1G_T = simulation(N, 1,G,T,β,γ,p_false_positive,R, true,S,p,Group,Trab2)
poll_mQua, poll_mInf, poll_mNFp, poll_mNQinc, poll_T = simulation(N, 1,G,T,β,γ,p_false_positive,R, "poll",S,p,Group,Trab2)
pollf10_mQua, pollf10_mInf, pollf10_mNFp, pollf10_mNQinc, pollf10_T = simulation(N, f,G,T,β,γ,p_false_positive,R, "poll",S,p,Group,Trab2)


f1 =plot(1:T-1,base_mQua[1:end-1],label ="base", lw=3)
plot!(1:T-1,ideal_mQua[1:end-1],label ="ideal", lw=3)
plot!(1:T-1,f10_mQua[1:end-1],label ="f=10 y G=20", lw=3)
plot!(1:T-1,f102_mQua[1:end-1],label ="f=10 y G=200 ", lw=3)
plot!(1:T-1,f1G_mQua[1:end-1],label ="f=1 y G=20 ", lw=3)
plot!(1:T-1,poll_mQua[1:end-1],label ="poll f=1 y G=20 ", lw=3)
plot!(1:T-1,pollf10_mQua[1:end-1],label ="poll f=10 y G=20 ", lw=3)
title!("Cuarentena")
xlabel!("Días")
ylabel!("# personas")

f2 =plot(1:T-1,base_mInf[1:end-1], label ="base", lw=3)
plot!(1:T-1,ideal_mInf[1:end-1], label ="ideal", lw=3)
plot!(1:T-1,f10_mInf[1:end-1], label ="f=10 G=20", lw=3)
plot!(1:T-1,f102_mInf[1:end-1], label ="f=10 G=200", lw=3)
plot!(1:T-1,f1G_mInf[1:end-1], label ="f=1 G=20", lw=3)
plot!(1:T-1,poll_mInf[1:end-1],label ="poll f=1 y G=20 ", lw=3)
plot!(1:T-1,pollf10_mInf[1:end-1],label ="poll f=10 y G=20 ", lw=3)
title!("Infección")
xlabel!("Días")
ylabel!("# personas")

f3 =plot(1:T-1,base_T[1:end-1], label ="base", lw=3)
plot!(1:T-1,ideal_T[1:end-1], label ="ideal", lw=3)
plot!(1:T-1,f10_T[1:end-1], label ="f=10 G=20", lw=3)
plot!(1:T-1,f102_T[1:end-1], label ="f=10 G=200", lw=3)
plot!(1:T-1,f1G_T[1:end-1], label ="f=1 G=20", lw=3)
plot!(1:T-1,poll_T[1:end-1],label ="poll f=1 y G=20 ", lw=3)
plot!(1:T-1,pollf10_T[1:end-1],label ="poll f=10 y G=20 ", lw=3)
title!("Test Realizados")
xlabel!("Días")
ylabel!("# Test")

f4 =plot(1:T-1,base_mNFp[1:end-1], label ="base", lw=4)
plot!(1:T-1,ideal_mNFp[1:end-1], label ="ideal", lw=4)
plot!(1:T-1,f10_mNFp[1:end-1], label ="f=10 G=20", lw=4)
plot!(1:T-1,f102_mNFp[1:end-1], label ="f=10 G=200", lw=4)
plot!(1:T-1,f1G_mNFp[1:end-1], label ="f=1 G=20", lw=4)
plot!(1:T-1,poll_mNFp[1:end-1],label ="poll f=1 y G=20 ", lw=4)
plot!(1:T-1,pollf10_mNFp[1:end-1],label ="poll f=10 y G=20 ", lw=4)
title!("Falsos Positivos")
xlabel!("Días")
ylabel!("# Test")


print("Caso: Falsos Positivos - Descubiertos en incubacion - Test ","Base:", round(sum(base_mNFp)),"-",round(sum(base_mNQinc)),"-",round(sum(base_T))," Ideal:",round(sum(ideal_mNFp)),"-",round(sum(ideal_mNQinc)),"-",round(sum(ideal_T))," Frecuencia 10 G=20:",round(sum(f10_mNFp)),"-", round(sum(f10_mNQinc)),"-",round(sum(f10_T))," Frecuencia 10 G=200:",round(sum(f102_mNFp)),"-"
, round(sum(f102_mNQinc)),"-",round(sum(f102_T))," Frecuencia 1 G=20:",round(sum(f1G_mNFp)),"-", round(sum(f1G_mNQinc)),"-",round(sum(f1G_T))," Poll Frecuencia 1 G=20:",round(sum(poll_mNFp)),"-", round(sum(poll_mNQinc)),"-",round(sum(poll_T))," Poll Frecuencia 10 G=20:",round(sum(pollf10_mNFp)),"-", round(sum(pollf10_mNQinc)),"-",round(sum(pollf10_T)))

plot(f1,f2,f3,f4, layout = (2,2))
