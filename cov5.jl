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



##################NUEVO#############################

#Recibiendo los parametros de person_p, un número de horas de turno (exposición al contagio) y un t_peak, mediante una exponencial
#estimando lambda como el que haga que en t_peak se alcanze P(tiempo para contagiarse<=t_peak)=probabilidad peak de contagio
#obtiene la probabilidad de haberse contagiado en ese turno
function exponential_p(t,t_inf,t_sin,t_rec,p_sinto,t_peak,horas_turno)
	prob=person_p(t,t_inf,t_sin,t_rec,p_sinto) #Se calcula el peak que se alcanza en t_peak
	lambda=(log(max(1-prob,0.00001)))/t_peak #Se calcula el lambda tal que P(tiempo para contagiarse<=t_peak)=probabilidad peak de contagio
	proba=1-exp(lambda*horas_turno) #Se calula P(tiempo para contagiarse<=horas_turno)
	return proba
end




##############################################################################


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
function group_p(G,M_cont,W)
	#G: indicatriz del grupo de trabajo
	#vp_cont: vector con las probabilidades de contagio
	#W: Peligro de contagio de pacientes
	V=((M_cont).*G)
	N=length(M_cont[1,:])
	p=zeros(N)
	for i in 1:N
		p[i]=(1-prod(-V[:,i].+1)) 
	end
	#de no contagiarse multiplicado por un factor que pondera las horas, este desde
	#la hora toma valores del estilo ~0.9 y mas.
	if W=="none"
		return p
	end
	if W=="low"
		return -(-p.+1).*(1-0.01).+1 #Esta probabilidades estan sin ninguna justificacion
		#podria depender de las horas de atencion a publico y medidas de seguridad
	end
	if W=="medium"
		return 1-(1-p)*(1-0.05)
	end
	if W=="high"
		return 1-(1-p)*(1-0.1)
	end
end

#Esta funcion entrega 
#mQua: Personas en cuarentena que deberian estar trabajando.
#mQua2: Personas en cuarentena en total.
#mInf: Personas infectadas trabajando (peligrosas)
#mInf2: Personas infectadas trabajando o mandadas a cuarentena (se agregan al anterior las que estan en cuarentena, pero les corresponde el turno, son las perdidas)
#mInf3: Personas infectadas en total

function simulation(N, f,G,T,β,γ,p_false_positive,R,test,prob,Group,Tra,random,cuarentena,z,t_peak,dias_atras,MF)
#N Tamaño grupo de trabajadores
#f Frecuencia testeo (puede ser vector de largo T, en caso de testear random de tamaño G o matriz de NxT con los individuos a testear en sus fechas)
# G Solo aplica para el caso random el número de test individuales o poll (por ejemplo G=4 son 4 poll testing y todos los sub test que se generen)
# T simulation horizon
# β house prob of contagion
#γ  prop of asymtomatic pop
#p_false_positive probability of a false positive
#R Replications
#test dice si es "poll", true(individual) o false (no hace).
#S Cantidad de grupos cerrados (los grupos tendran N/S integrantes)
#prob Prob peak de infeccion
#Group un vector con indicatrices de tamaño N con los integrantes de cada grupo y riesgo de atención, tamaño (N+1)*S
#Tra: Matriz (2*N,T) de turnos de trabajo
#random: Si es random o estan fijos a quienes se testeara (no los dias esos estan en f, solo los funcionarios a testear)
#cuarentena: Si se mandan a todo el grupo en cuarentena o solo al infectado (ahora tambien mixto, solo al grupo cuando uno presenta sintomas)
#z: dias de cuarntena al grupo, no se usa en caso de solo mandar al positivo
#t_peak: tiempo en que la exponencial alcanza la probabilidad peak
#dias_atras: cuantos dias hacia atras se mira con quienes has tenido contacto estrecho
#MF: Matriz de funcionarios, contiene sus interacciones con los pacientes y entre ellos.
hr=12 
MD=zeros(N,T) #Matriz trabajo dia 
MN=zeros(N,T) ##Matriz trabajo noche 
for i in 1:N
	MD[i,:]=Tra[2*i-1,:]
	MN[i,:]=Tra[2*i,:]
end
	
S=length(Group)/(N+1) #Tamaño grupo
if (cuarentena=="grupo")|(cuarentena=="mixto")
	v=vecgroup(Group)
end

g_t=ones(T)
	for tt in 1:T
		if sum((MD[:,tt]+MN[:,tt]).>=1)>=1
			g_t[tt]=sum((MD[:,tt]+MN[:,tt]).>=1)
		end
	end




#Vectores que guardan los valores de cuarentena, infectados y números de test en cada t
    mQua = zeros(T)
    mQua2 = zeros(T)
    mFQua = zeros(T)    
    mInf = zeros(T)
    mInf2 = zeros(T)
    mInf3 = zeros(T)
    mFInf=zeros(T)
    mFInf2=zeros(T)
    mSy= zeros(T)
    mFSy=zeros(T)
	mNTest=zeros(T)
	mNFpositive=zeros(T)
m=length(MF[1,:,1])
Gf=[zeros(N)' ones(m)']'

    # Loop over replications
    for r = 1:R
	print(r) #Lo uso para saber si corre 
        tEx =zeros(N+m) # time of exposure
	tIn =zeros(N+m) # time of infecciousness
        tSy =zeros(N+m) # time of development of symptoms
        tRe =zeros(N+m) # time of recovery
        As = rand(N+m).< γ # asymtomatic lottery
        Su = ones(N+m,T) # susceptible population
        Ex = zeros(N+m,T) # exposed population
        In = zeros(N+m,T) # infectious population
        Sy = zeros(N+m,T) # symptomatic population
        Re = zeros(N+m,T) # recovered population
        Qu = zeros(N+m,T) # quarantined population
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
		new_inf=zeros(N+m) #Vector de los nuevos infectados
	    	con=0 #Cuantos grupos hay para hacer poll testing
	    	con2=0 #Cuantos grupos hay para hacer poll testing no random
		MC=zeros(N+m,N+m) #Matriz de contagio
		for a=1:S
			info=Group[Int((a-1)*(N+1)+1):Int((a)*(N+1))] #Se recupera la info codificada del grupo
			g=[Int.(info[1:N])' zeros(m)' ]' #Indicatriz del grupo
			W=info[N+1] #Riesgo del grupo
			vec=zeros(N) #Vector auxiliar para las probabilidades de los infectados
			#del grupo
			for i in findall((g.*inf.*(-qu.+1).*(-re.+1)).==1)
				for j in findall((g.*su.*(-qu.+1)).==1)
					if (MD[i,t]*MD[j,t]==1) & ((1-MN[i,t])*(1-MN[j,t])==1)  #Condicion esta en el grupo, infectado
						#y en cuarentena
						#Recuperacion de los tiempos
						tinf=tIn[Int(i)]-tEx[Int(i)]
						tsin=tSy[i]-tEx[i]
						trec=tRe[i]-tEx[i]
						MC[i,j]=exponential_p(t,tinf,tsin,trec,prob,t_peak,hr)
				 		#Se agregan las probabilidades de los enfermos
					end
					if ((1-MD[i,t])*(1-MD[j,t])==1) & (MN[i,t]*MN[j,t]==1)  #Condicion esta en el grupo, infectado
						#y en cuarentena
						#Recuperacion de los tiempos
						tinf=tIn[Int(i)]-tEx[Int(i)]
						tsin=tSy[i]-tEx[i]
						trec=tRe[i]-tEx[i]
						MC[i,j]=exponential_p(t,tinf,tsin,trec,prob,t_peak,hr)
				 		#Se agregan las probabilidades de los enfermos
					end
					if (MD[i,t]*MD[j,t]==1) & (MN[i,t]*MN[j,t]==1)  #Condicion esta en el grupo, infectado
						#y en cuarentena
						#Recuperacion de los tiempos	
						tinf=tIn[Int(i)]-tEx[Int(i)]
						tsin=tSy[i]-tEx[i]
						trec=tRe[i]-tEx[i]
						MC[i,j]=exponential_p(t,tinf,tsin,trec,prob,t_peak,2*hr)
						#Se agregan las probabilidades de los enfermos
					end
				end
				for f in findall((su[N+1:N+m].*(-qu[N+1:N+m].+1)).==1)
					if (MF[2*i-1,f,t]==1) & (MF[2*i,f,t]==0)  #Condicion esta en el grupo, infectado
						#y en cuarentena
						#Recuperacion de los tiempos
						tinf=tIn[Int(i)]-tEx[Int(i)]
						tsin=tSy[i]-tEx[i]
						trec=tRe[i]-tEx[i]
						MC[i,N+f]=exponential_p(t,tinf,tsin,trec,prob,t_peak,hr)
				 		#Se agregan las probabilidades de los enfermos
					end
					if (MF[2*i-1,f,t]==0) & (MF[2*i,f,t]==1)  #Condicion esta en el grupo, infectado
						#y en cuarentena
						#Recuperacion de los tiempos
						tinf=tIn[Int(i)]-tEx[Int(i)]
						tsin=tSy[i]-tEx[i]
						trec=tRe[i]-tEx[i]
						MC[i,N+f]=exponential_p(t,tinf,tsin,trec,prob,t_peak,hr)
				 		#Se agregan las probabilidades de los enfermos
					end
					if (MF[2*i-1,f,t]==1) & (MF[2*i,f,t]==1)  #Condicion esta en el grupo, infectado
						#y en cuarentena
						#Recuperacion de los tiempos
						tinf=tIn[Int(i)]-tEx[Int(i)]
						tsin=tSy[i]-tEx[i]
						trec=tRe[i]-tEx[i]
						MC[i,N+f]=exponential_p(t,tinf,tsin,trec,prob,t_peak,2*hr)
				 		#Se agregan las probabilidades de los enfermos
					end	
				end					
			end
			for f in findall((inf[N+1:N+m].*(-qu[N+1:N+m].+1).*(-re[N+1:N+m].+1)).==1)
				for j in findall((g.*su.*(-qu.+1)).==1)
					if (MF[2*j-1,f,t]==1) & (MF[2*j,f,t]==0)  #Condicion esta en el grupo, infectado
						#y en cuarentena
						#Recuperacion de los tiempos
						tinf=tIn[Int(N+f)]-tEx[Int(N+f)]
						tsin=tSy[N+f]-tEx[N+f]
						trec=tRe[N+f]-tEx[N+f]
						MC[N+f,j]=exponential_p(t,tinf,tsin,trec,prob,t_peak,hr)
				 		#Se agregan las probabilidades de los enfermos
					end	
					if (MF[2*j-1,f,t]==0) & (MF[2*j,f,t]==1)  #Condicion esta en el grupo, infectado
						#y en cuarentena
						#Recuperacion de los tiempos
						tinf=tIn[Int(N+f)]-tEx[Int(N+f)]
						tsin=tSy[N+f]-tEx[N+f]
						trec=tRe[N+f]-tEx[N+f]
						MC[N+f,j]=exponential_p(t,tinf,tsin,trec,prob,t_peak,hr)
				 		#Se agregan las probabilidades de los enfermos
					end
					if (MF[2*j-1,f,t]==1) & (MF[2*j,f,t]==1)  #Condicion esta en el grupo, infectado
						#y en cuarentena
						#Recuperacion de los tiempos
						tinf=tIn[Int(N+f)]-tEx[Int(N+f)]
						tsin=tSy[N+f]-tEx[N+f]
						trec=tRe[N+f]-tEx[N+f]
						MC[N+f,j]=exponential_p(t,tinf,tsin,trec,prob,t_peak,2*hr)
				 		#Se agregan las probabilidades de los enfermos
					end
				end
				
			end
			PG=[group_p(g,MC[:,1:N],W)' zeros(m)']' #Vector con probabilidades de contagio dentro del grupo
			r=rand(N+m)
			#Reemplazar el beta y 1.2 beta por una probabilidad de contagio en casa idividual.
			new_inf+=(-qu.+1).*su.*g.*(r.<(-(-PG.+1).*(1-β).+1)) #Vector de infectados del grupo que no esten en cuarentena, esten o no trabajando
			new_inf+=qu.*su.*g.*(r.<(β)) #Sumar infectados en cuarentena (falsos positivos)
			#Para el poll test, predefino el numero de grupos que se pueden formar
			#se ve quienes son el publico a testear
			if test=="poll"
				I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*(-qu.+1).*(tra).*collect(1:N) #Gente que esta trabajando, no en cuarentena y nunca a presentado sintomas
				JJ=shuffle!(g.*I) #Se mezclan
				III=setdiff(JJ,[0]) #Se eliminan los 0 Se cuentan
				if length(III)==0
					con+=0
				else
					ind=Int(round(length(III)/5))
					r=length(III)%5
					con+=(ind+(r!=0))
				end	
				#Para el poll test no random, es similar. Se hacen distintos, pues en caso random la frecuencia es un vector
				if random=="no"
					I2=f[:,t].*( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*(-qu.+1).*collect(1:N)
					JJ2=(g.*I2)
					III2=setdiff(JJ2,[0])
					if length(III2)==0
						con2+=0
					else
						ind2=Int(round(length(III2)/5))
						r2=length(III2)%5
						con2+=(ind2+(r2!=0))
					end
				end
			end
	   	end
		

############################
		
		for f in findall((inf[N+1:N+m].*(-qu[N+1:N+m].+1).*(-re[N+1:N+m].+1)).==1)
			for f2 in findall((su[N+1:N+m].*(-qu[N+1:N+m].+1)).==1)
					if (MF[2*f2-1,f,t]==1) & (MF[2*f2,f,t]==0)  #Condicion esta en el grupo, infectado
						#y en cuarentena
						#Recuperacion de los tiempos
						tinf=tIn[Int(N+f)]-tEx[Int(N+f)]
						tsin=tSy[N+f]-tEx[N+f]
						trec=tRe[N+f]-tEx[N+f]
						MC[N+f,N+f2]=exponential_p(t,tinf,tsin,trec,prob,t_peak,hr)
				 		#Se agregan las probabilidades de los enfermos
					end
					if (MF[2*f2-1,f,t]==0) & (MF[2*f2,f,t]==1)  #Condicion esta en el grupo, infectado
						#y en cuarentena
						#Recuperacion de los tiempos
						tinf=tIn[Int(N+f)]-tEx[Int(N+f)]
						tsin=tSy[N+f]-tEx[N+f]
						trec=tRe[N+f]-tEx[N+f]
						MC[N+f,N+f2]=exponential_p(t,tinf,tsin,trec,prob,t_peak,hr)
				 		#Se agregan las probabilidades de los enfermos
					end
					if (MF[2*f2-1,f,t]==1) & (MF[2*f2,f,t]==1)  #Condicion esta en el grupo, infectado
						#y en cuarentena
						#Recuperacion de los tiempos
						tinf=tIn[Int(N+f)]-tEx[Int(N+f)]
						tsin=tSy[N+f]-tEx[N+f]
						trec=tRe[N+f]-tEx[N+f]
						MC[N+f,N+f2]=exponential_p(t,tinf,tsin,trec,prob,t_peak,2*hr)
				 		#Se agregan las probabilidades de los enfermos
					end	
				end	
		end
	PGf=-prod((-MC.*Gf).+1,dims=1).+1
	r=rand(N+m)
	new_inf+=(-qu.+1).*su.*Gf.*(r.<(-(-PGf'.+1).*(1-β).+1)) #Vector de infectados del grupo que no esten en cuarentena, esten o no trabajando
	new_inf+=qu.*su.*Gf.*(r.<(β)) #Sumar infectados en cuarentena (falsos positivos)		
###############################


		
           	for j in findall(new_inf.==1)                   # cicle over newly infected
			i=j[1]
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
		if (t>=3) &((cuarentena=="mixto")|(cuarentena=="grupo"))
			for i in findall(Sy[:,t].*(Qu[:,t]).*Sy[:,t-1].*(-Qu[:,t-1].+1).*Sy[:,t-2].*(-Qu[:,t-2].+1).==1)
				for t2 in findall(MD[i,Int(max(t-dias_atras,1)):t].==1)
					Aux=MD[:,t2]				
					for y in findall(Aux.*(-qu.+1).*(Group[Int((N+1)*(v[i]-1)+1):Int((N+1)*(v[i]-1)+N)]).==1)
						Qu[y,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
					end
				end
				for t3 in findall(MN[i,Int(max(t-dias_atras,1)):t].==1)
					Aux2=MN[:,t3]
					for y in findall(Aux2.*(-qu.+1).*(Group[Int((N+1)*(v[i]-1)+1):Int((N+1)*(v[i]-1)+N)]).==1)
						Qu[y,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
					end
				end
			end
		end
            # Testing
		#Vector de funcionarios que testear
	    cand=[]
		if random=="si" #Se realiza test al azar 
			if f[t]==1 #Se testea estos dias
	    			if test==true #test individual
					I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*(tra).*(-qu.+1).*collect(1:N) #Se ven los candidatos a testear
					J=shuffle!(I) #Mezclar
					II=setdiff(J,[0]) #Eliminar ceros
					cand=Int.(II[1:min(G,length(II))])  #Se toman como candidatos a los posibles y el minimo entre ellos y G
	    			end
				if test=="poll"
					h=shuffle!(collect(1:con))[1:min(G,con)] #Se mezclan los grupos y se escogen los primeros min(G ,con)
					b=sort(h) #Se ordenan los escogidos de menor a mayor
					for a=1:S #Se recorren los grupos de trabajo, para los indices previamente escogidos realizar poll testing
						info=Group[Int((a-1)*(N+1)+1):Int((a)*(N+1))] 
						g=Int.(info[1:N]) #Indicatriz del grupo
						I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*(-qu.+1).*collect(1:N) #Se crean los candidatos del grupo a testear
						JJ=shuffle!(g.*I)
						III=setdiff(JJ,[0])
						if length(III)>=1 #De exisitr algun candidato se realiza este esquema para ir testeando
							ind=Int(floor(length(III)/5))
							r=length(III)%5
							lim=(ind+(r!=0))
							c=setdiff((b.<=lim).*(b.>=1).*b,[0])
							if length(c)>=1
								b.-lim
								for cc in c
									if cc-1==ind 	 #Por ejemplo el grupo son 14 por lo que hay 2 de 5 y uno de 4
										cand=Int.(III[5*ind+1:5*ind+r])
									else	
										cand=Int.(III[5*(cc-1)+1:5*cc]) #El grupo es exacto al dividirlo en 5
									end	
									#hacer poll testing
									s=0
									pp=1
									for i in cand #Se realizan los test
										p_positive = true_positive(t-tEx[i]+1, tIn[i]-tEx[i]+1, tSy[i]-tEx[i]+1, tRe[i]-tEx[i]+1)  #Prob de ser positivo
										s+=su[i] #Cuantos sanos hay en el poll
										pp=pp*(1-(1-su[i])*p_positive)	#Probabilidad de que no sea positivo el poll
									end
									PP=1-pp #Probabilidad de que al menos uno de los del  poll sea positivo
									NTest[t]+=1 #Se realiza un poll test
									if (s==length(cand)) & (rand() < (p_false_positive))	#Condición de que es un poll falso positivo	
										for i in cand #Se realizan los sub test para encontrar al positivo
											NTest[t]+=1
											if (su[i]==1) & (rand() < p_false_positive) #Nuevamente tenemos un falso positivo
                    										Qu[i,min(T,t+1):min(T,t+15)] .= 1 #if someone tests (false) positive, quarantined for 2 week
												if cuarentena=="grupo"
													for y in findall((-qu.+1).*(Group[Int((N+1)*(v[i]-1)+1):Int((N+1)*(v[i]-1)+N)]).==1)
														Qu[y,min(T,t+1):min(T,t+1+z)] .= 1 #quarantined for 2 weeks
													end
												end
												NFpositive[t]+=1			
											end
										end
									else
										if (rand()<=PP) & (s<=(length(cand)-1)) #Es un test de verdad con algun positivo
											for i in cand #Se realizan los sub test
												NTest[t]+=1
                										if su[i] == 0
                    											if (qu[i] == 0) & (sum(Qu[i,1:t])==0)
                        											p_positive = true_positive(t-tEx[i]+1, tIn[i]-tEx[i]+1, tSy[i]-tEx[i]+1, tRe[i]-tEx[i]+1)
                        											if rand() < p_positive #Se encuentra un individuo positivo
															Qu[i,Int(min(T,t+1)):Int(min(T,max(t+15,tRe[i]+3)))] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
															if cuarentena=="grupo"
																for y in findall((-qu.+1).*(Group[Int((N+1)*(v[i]-1)+1):Int((N+1)*(v[i]-1)+N)]).==1)
																	Qu[y,min(T,t+1):min(T,t+1+z)] .= 1 #quarantined for 2 weeks
																end
															end
														end
                        										end
                    										else
													if (rand() < p_false_positive) #Se encuentra un falso positivo
                    												Qu[i,min(T,t+1):min(T,t+15)] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
														if cuarentena=="grupo"
															for y in findall((-qu.+1).*(Group[Int((N+1)*(v[i]-1)+1):Int((N+1)*(v[i]-1)+N)]).==1)
																Qu[y,min(T,t+1):min(T,t+1+z)] .= 1 #quarantined for 2 weeks
															end
														end
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
				if test!="poll"	#Caso individual	
					for i in cand	#se recorren los factibles a testear
						NTest[t]+=1 #Se hace un test
                				if (su[i] == 1) & (rand() < p_false_positive) #Falso positivo
                    					Qu[i,min(T,t+1):min(T,t+15)] .= 1 # if someone tests (false) positive, quarantined for 2 weeks		
							if cuarentena=="grupo"
								for y in findall((-qu.+1).*(Group[Int((N+1)*(v[i]-1)+1):Int((N+1)*(v[i]-1)+N)]).==1)
									Qu[y,min(T,t+1):min(T,t+1+z)] .= 1 #quarantined for 2 weeks
								end
							end
							NFpositive[t]+=1
                				end
                				if su[i] == 0 #Esta o estuvo enfermo
                    					if (qu[i] == 0) & (sum(Qu[i,1:t])==0) #nunca estuvo en cuarentena
                        					p_positive = true_positive(t-tEx[i]+1, tIn[i]-tEx[i]+1, tSy[i]-tEx[i]+1, tRe[i]-tEx[i]+1)
                        					if rand() < p_positive #Probabilidad de que el test de positivo
        	                    					Qu[i,Int(min(T,t+1)):Int(min(T,max(t+15,tRe[i]+3)))] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
									if cuarentena=="grupo"
										for y in findall((-qu.+1).*(Group[Int((N+1)*(v[i]-1)+1):Int((N+1)*(v[i]-1)+N)]).==1)
											Qu[y,min(T,t+1):min(T,t+1+z)] .= 1 #quarantined for 2 weeks
										end
									end	
                        					end
                    					end
                				end
            				end	
        			end
			end
		end
#######################
		if random=="no" #Se realiza el esquema de test de manera similar solo que los candidatos ya estan dados por defecto y solo saca los que esten en cuarentena o ya hayan tenido sintomas en el pasado
			if sum(( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*f[:,t].*(-qu.+1))>=1	
	  	  		if test==true
					J=Int.(( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*f[:,t].*(-qu.+1).*collect(1:N+m))	
					II=setdiff(J,[0])
					cand=Int.(II)
	    			end
				if test=="poll"
					h=collect(1:con2)[1:min(G,con2)]
					b=sort(h)
					for a=1:S
						info=Group[Int((a-1)*(N+1)+1):Int((a)*(N+1))]
						g=Int.(info[1:N]) #Indicatriz del grupo
						I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*(-qu.+1).*collect(1:N).*f[:,t]
						JJ=(g.*I)
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
                	    									Qu[i,min(T,t+1):min(T,t+15)] .= 1 #if someone tests (false) positive, quarantined for 2 weeks
												if cuarentena=="grupo"
													for y in findall((-qu.+1).*(Group[Int((N+1)*(v[i]-1)+1):Int((N+1)*(v[i]-1)+N)]).==1)
														Qu[y,min(T,t+1):min(T,t+1+z)] .= 1 #quarantined for 2 weeks
													end
												end	
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
															if cuarentena=="grupo"
																for y in findall((-qu.+1).*(Group[Int((N+2)*(v[i]-1)+1):Int((N+2)*(v[i]-1)+N)]).==1)
																	Qu[y,min(T,t+1):min(T,t+1+z)] .= 1 #quarantined for 2 weeks
																end
															end
														end
                        										end
                    										else
													if (rand() < p_false_positive)
                    												Qu[i,min(T,t+1):min(T,t+15)] .= 1 #if someone tests (false) positive, quarantined for 2 weeks # if someone tests (false) positive, quarantined for 2 weeks					
														if cuarentena=="grupo"
															for y in findall((-qu.+1).*(Group[Int((N+2)*(v[i]-1)+1):Int((N+2)*(v[i]-1)+N)]).==1)
																Qu[y,min(T,t+1):min(T,t+1+z)] .= 1 #quarantined for 2 weeks
															end
														end
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
				if (test!="poll")	& (length(cand)>=1)			
					for i in cand	
						NTest[t]+=1
           				     	if (inf[i] == 0) & (rand() < p_false_positive)
			                    		Qu[i,min(T,t+1):min(T,t+15)] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
							if cuarentena=="grupo"
								for y in findall((-qu.+1).*(Group[Int((N+1)*(v[i]-1)+1):Int((N+1)*(v[i]-1)+N)]).==1)
									Qu[y,min(T,t+1):min(T,t+1+z)] .= 1 #quarantined for 2 weeks
								end
							end
							NFpositive[t]+=1
                				end
                				if inf[i] == 1
                    					#if (qu[i] == 0) & (sum(Qu[i,1:t])==0)
                        					p_positive = true_positive(t-tEx[i]+1, tIn[i]-tEx[i]+1, tSy[i]-tEx[i]+1, tRe[i]-tEx[i]+1)
            			            			if rand() < p_positive
                            						Qu[i,Int(min(T,t+1)):Int(min(T,max(t+15,tRe[i]+3)))] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
									if cuarentena=="grupo"
										for y in findall((-qu.+1).*(Group[Int((N+1)*(v[i]-1)+1):Int((N+1)*(v[i]-1)+N)]).==1)
											Qu[y,min(T,t+1):min(T,t+1+z)] .= 1 #quarantined for 2 weeks
										end
									end
                	        				end
               				     		#end
                				end
            				end
		        	end	    
			end
		end
	end
	#Se guardan los valores en cada
        mQua .+= dropdims(sum(   (((MD).*Qu[1:N,:]+(MN).*Qu[1:N,:]).>=1)     ,dims=1),dims=1)
	mQua2 .+= dropdims(sum(Qu[1:N,:],dims=1),dims=1)
	mFQua .+= dropdims(sum(Qu[N+1:N+m,:],dims=1),dims=1)
        mInf .+= dropdims(sum(      (((MD).*In[1:N,:].*(-Qu[1:N,:].+1)+(MN).*In[1:N,:].*(-Qu[1:N,:].+1)).>=1)          ,dims=1),dims=1)
        mInf2 .+= dropdims(sum(    (((MD).*In[1:N,:]+(MN).*In[1:N,:]).>=1)          ,dims=1),dims=1)
	mInf3 .+= dropdims(sum(In[1:N,:],dims=1),dims=1)
	alpha=zeros(m,T)
for i in 1:m
	alpha[i,:]=(sum(MF[:,i,:],dims=1).>=1)
end
        mFInf .+= dropdims(sum(In[N+1:N+m,:].*(-Qu[N+1:N+m,:].+1).*alpha[1:m,:],dims=1),dims=1)
        mFInf2 .+= dropdims(sum(In[N+1:N+m,:],dims=1),dims=1)
	mSy .+= dropdims(sum(Sy[1:N,:],dims=1),dims=1)
	mFSy .+= dropdims(sum(Sy[N+1:N+m,:],dims=1),dims=1)
	mNFpositive .+= NFpositive	
	mNTest .+= NTest
    end
#Se calcula un promedio de las simulaciones
    mQua = mQua/R
    mQua2 = mQua2/R
    mFQua=mFQua/R
    mInf = mInf/R
    mInf2 = mInf2/R
    mInf3 = mInf3/R
    mFInf=mFInf/R
    mFInf2=mFInf2/R	
    mSy= mSy/R
    mFSy=mFSy/R
	mNFpositive=mNFpositive/R
	mNTest=mNTest/R


    return mQua,mQua2,mInf,mInf2,mInf3, mNFpositive,  mNTest, g_t, mSy,mFQua,mFInf,mFInf2,mFSy
end


N = 154             # Tamaño grupo de trabajadores
#f = 10              # Frecuencia testeo (cada T+1 dias)
#G = Int(round(N/f)) # tamaño grupo a testear - por ahora un numero entero
T = 210              # simulation horizon
β = 0.01   #prob of contagion en casa
γ = 0.3 #(antes 0.25)           # prop of asymtomatic pop
p_false_positive = 0.01 # probability of a false positive
R = 400    # Replications
p=0.01 #Prob peak de infeccion



S=Int(N/22) #Cantidad de grupos
    #Creacion de un grupo de ejemplo
    for a=1:S
	global Group
	if a==1
		g1=[zeros(Int((a-1)*(N/S)))' ones(Int(N/S))' zeros(Int(N-a*(N/S)))']
		g2="low"
		Group=[g1 g2]
	else

		g1=[zeros(Int((a-1)*(N/S)))' ones(Int(N/S))' zeros(Int(N-a*(N/S)))']
		g2="low"
		g=[g1 g2]
		Group=[Group g]
	end
    end


#Recibe un grupo y entrega un vector con los indices de los grupos a los que pertenecen
function vecgroup(group)
S=Int(length(group)/(N+1))
A=zeros(N,S)
v=zeros(N)
for b in 1:S
	for a in 1:N
		A[a,b]=group[Int.((N+1)*(b-1)+a)]
		
	end
end
for i in 1:N
v[i]=Int.(setdiff((A[i,:]).*collect(1:S),[0]))[1]


end
return v
end

#Esquemas de trabajo
#Todos los dias
Trab=ones(N,T)
#7 dias intercaldos
Trab2=zeros(N,T)
for i in 1:Int(floor(T/7))
	for j in 1:N
		if i==Int(floor(T/7)) & (T%7!=0)
			if (mod(j,2)==1) & (mod(i,2)==0)
				Trab2[j,(i)*7+1:(i)*7+T%7]=ones(T%7)
			end
			if (mod(j,2)==0) & (mod(i,2)==1)
				Trab2[j,(i)*7+1:(i)*7+T%7]=ones(T%7)
			end
		end 
		if (mod(j,2)==1) & (mod(i,2)==1)
			Trab2[j,(i-1)*7+1:i*7]=ones(7)
			
		end
		if (mod(j,2)==0) & (mod(i,2)==0)
			Trab2[j,(i-1)*7+1:i*7]=ones(7)
		end	
	end
end
#Trabajan 1-6 de a 4 grupos de 10
Trab3=zeros(N,T)
for i in 1:Int(floor(T))
	for j in 1:7 
		if (mod(i,7)==mod(j,7))
			Trab3[Int((N/7)*(j-1)+1):Int((N/7)*j),i]=ones(Int(N/7))		
		end
	end
end
#Trabajan 2-12 de a 4 grupos de 10
Trab4=zeros(N,T)
for i in 1:Int(floor(T/2))
	for j in 1:7 
		if (mod(i,7)==mod(j,7))
			Trab4[Int((N/7)*(j-1)+1):Int((N/7)*j),2(i-1)+1:2*i]=ones(Int(N/7),2)		
		end
	end
end

#######################################
#No se esta usando pasar a la otra #####################
function frecuencia(N,T,f,d,G,random)
	if f==0
		return zeros(T)
	end
	if random=="si"
		Frec=zeros(T)
		for t in 1:T
			if (mod(t-1,f)==0) & (t+d<=T)
				Frec[t+d]=1
			end
		end
	return Frec



	end
	if random=="no"
		Frec=zeros(N,T)
		ii=1
		iii=1
		for t in 1:T
			for i in 1:Int(floor(N/G))
				if (mod(t-1,f)==0) & (t+d<=T)
					if (ii==i) & (iii==1)
						Frec[G*(i-1)+1:G*i,t+d]=ones(G)
						ii+=1
						iii=0
						if ii==Int(floor(N/G))+1
							ii=1
						end
					end
				end
			end
			iii=1
		end
		return Frec
	end
end	
	
fr0=frecuencia(N,T,0,0,200,"si")
fr1200=frecuencia(N,T,1,0,200,"si")
######################
#Frecuencias de testeos no random
#Se testean 2 grupos de los 4 que trabajan 1-7 cada dia, es decir se testean todos en 14 dias
fr1=zeros(N,T)
for i in 1:Int(floor(T))
	for j in 1:7 
		if (mod(i,7)==mod(j,7))
			if (mod(i,14)<=7)& (mod(i,14)>=1)
				fr1[Int((N/7)*(j-1)+1):Int((N/7)*(j-1)+(N/14)),i]=ones(Int(N/14))	
			else
				fr1[Int((N/7)*(j-1)+(N/14)+1):Int((N/7)*(j)),i]=ones(Int(N/14))
			end
				
		end
	end
end
#Se testean 2 grupos de los 4 que trabajan 2-12 cada dia en, es decir se testean todos en 14 dias
fr2=zeros(N,T)
for i in 1:Int(floor(T))
	for j in 1:14
		if (mod(i,14)==mod(j,14))
			fr2[Int((N/14)*(j-1)+1):Int((N/14)*(j)),i]=ones(Int(N/14))
		end
	end
end

#Se testean 2 grupos de los 4 que trabajan 1-6 cada dia antes de que trabajen, es decir se testean todos en 14 dias
fr12=zeros(N,T)
for i in 2:(Int(floor(T))+1)
	for j in 1:7 
		if (mod(i,7)==mod(j,7))
			if (mod(i,14)<=7)& (mod(i,14)>=1)
				fr12[Int((N/7)*(j-1)+1):Int((N/7)*(j-1)+(N/14)),i-1]=ones(Int(N/14))	
			else
				fr12[Int((N/7)*(j-1)+(N/14)+1):Int((N/7)*(j)),i-1]=ones(Int(N/14))
			end
				
		end
	end
end


#Se testean 2 grupos de los 4 que trabajan 2-12 cada dia antes de que trabajen, es decir se testean todos en 14 dias
fr22=zeros(N,T)
for i in 3:(Int(floor(T))+2)
	for j in 1:14 
		if (mod(i,14)==mod(j,14))
			fr22[Int((N/14)*(j-1)+1):Int((N/14)*(j)),i-2]=ones(Int(N/14))
		end
	end
end
#Se testean 2 grupos de los 4 que trabajan 2-12 cada dia antes de que trabajen, es decir se testean todos en 14 dias
fr23=zeros(N,T)
for i in 2:(Int(floor(T))+1)
	for j in 1:14 
		if (mod(i,14)==mod(j,14))
			fr23[Int((N/14)*(j-1)+1):Int((N/14)*(j)),i-1]=ones(Int(N/14))
		end
	end
end


#Se testean los 4 grupos antes de trabajar 1-6, es decir en 7 dias se testean a todos
fr3=zeros(N,T)
for i in 2:(Int(floor(T))+1)
	for j in 1:7 
		if (mod(i,7)==mod(j,7))
			fr3[Int((N/7)*(j-1)+1):Int((N/7)*j),i-1]=ones(Int(N/7))		
		end
	end
end
#Se testean los 4 grupos antes de trabajar 2-12, es decir en 14 dias se testean a todos 2 veces consecutivas
fr4=zeros(N,T)
for i in 2:(Int(floor(T/2))+1)
	for j in 1:7 
		if (mod(i,7)==mod(j,7))
			fr4[Int((N/7)*(j-1)+1):Int((N/7)*j),2(i-2)+1:2*(i-1)]=ones(Int(N/7),2)		
		end
	end
end
#Es la anterior pero se testea un dia antes de trabajar y al final del primer turno
fr42=zeros(N,T)
for i in 2:(Int(floor(T))+1)
	for j in 1:7 
		if (mod((i+1)/2,7)==mod(j,7))&(mod(i,2)==1)
			fr42[Int((N/7)*(j-1)+1):Int((N/7)*j),i-1]=ones(Int(N/7),1)		
		end
		if (mod(i/2,7)==mod(j,7))&(mod(i,2)==0)
			fr42[Int((N/7)*(j-1)+1):Int((N/7)*j),i-1]=ones(Int(N/7),1)		
		end
	end
end
#Otra variante se testea el dia anterior al turno y al final del ultimo dia de trabajo
fr43=zeros(N,T)
for i in 2:(Int(floor(T))+1)
	for j in 1:7 
		if (mod((i+1)/2,7)==mod(j,7))&(mod(i,2)==1)
			fr43[Int((N/7)*(j-1)+1):Int((N/7)*j),i-1]=ones(Int(N/7),1)		
		end
		if (mod(i/2,7)==mod(j,7))&(mod(i,2)==0)&(i+1<=T)
			fr43[Int((N/7)*(j-1)+1):Int((N/7)*j),i+1]=ones(Int(N/7),1)		
		end
	end
end
#Otra variante se testea el dia anterior al turno a los 40
fr44=zeros(N,T)
for i in 2:(Int(floor(T))+1)
	for j in 1:7 
		if (mod((i+1)/2,7)==mod(j,7))&(mod(i,2)==1)
			fr44[Int((N/7)*(j-1)+1):Int((N/7)*j),i-1]=ones(Int(N/7),1)		
		end
		
	end
end
##############################################################################

t_peak=24
##########################################
#Ejemplos de uso poll
#poll_mQua,poll_mQua2, poll_mInf,poll_mInf2,poll_mInf3, poll_mNFp, poll_T = simulation(N, fr22,4,T,β,γ,p_false_positive,R, "poll",p,Group,Trab4,"no","mixto",z,t_peak)
#poll2_mQua,poll2_mQua2, poll2_mInf,poll2_mInf2,poll2_mInf3, poll2_mNFp, poll2_T = simulation(N, fr43,8,T,β,γ,p_false_positive,R, "poll",p,Group,Trab4,"no","mixto",z,t_peak)
#########################################
#Matriz turnos 
TUR=[1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0;
1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0;
0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0;
0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0;  
0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0;
0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0;
0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0;
0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0;
0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0;
0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0;
0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1;
0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1;
1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0;
1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0;
0 1 0 1	0 0 0 0	0 0 0 0	0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1	0 0 0 0	0 0 0 0	0 0;
0 1 0 1	0 0 0 0	0 0 0 0	0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1	0 0 0 0	0 0 0 0	0 0;
0 1 0 1	0 0 0 0	0 0 0 0	0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1	0 0 0 0	0 0 0 0	0 0;
0 1 0 1	0 0 0 0	0 0 0 0	0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1	0 0 0 0	0 0 0 0	0 0;
0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0;
0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0;
0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0;
0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0;
0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0;
0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1;
0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1;
1 0 1 0 0 0 0 1 0 1 0 0 0 0 1 0 1 0 0 0 0 1 0 1 0 0 0 0 1 0 1 0 0 0 0 1 0 1 0 0 0 0;
1 0 1 0 0 0 0 1 0 1 0 0 0 0 1 0 1 0 0 0 0 1 0 1 0 0 0 0 1 0 1 0 0 0 0 1 0 1 0 0 0 0;
0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 1 0 1 0;
0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 1 0 1 0;
0 0 0 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 1 0 1 0 1 0 1 0 0 0;
0 0 0 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 1 0 1 0 1 0 1 0 0 0;
0 0 0 0 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 1 0 1 0 1 0 1 0 0;
0 0 0 0 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 1 0 1 0 1 0 1 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 1 1 0 0 0 0 0 1 1 0 0 0 0 0 1 1 0 0 0 0 0 1 1 0 0 0 0 0 1 1 0 0 0 0 0 1 1 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 1 1 0 0 0 0 0 1 1 0 0 0 0 0 1 1 0 0 0 0 0 1 1 0 0 0 0 0 1 1 0 0 0 0 0 1 1;
0 1 1 0 0 1 1 0 1 1 0 0 1 1 0 1 1 0 0 1 1 0 1 1 0 0 1 1 0 1 1 0 0 1 1 0 1 1 0 0 1 1;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1]

#Turnos replicados
Tur=zeros(2*N,T)
for i in 1:Int((2*N)/length(TUR[:,1]))
	for t in 1:Int((T)/length(TUR[1,:]))
		Tur[Int((i-1)*length(TUR[:,1])+1):Int(i*length(TUR[:,1])),Int((t-1)*length(TUR[1,:])+1):Int((t)*length(TUR[1,:]))]=TUR
	end
end
#Testear antes de los turnos
fTur=zeros(N,T)  
for i in 1:N
	for t in 2:T
		if (Tur[Int(2*i-1),t]==1) | (Tur[Int(2*i),t]==1) 
			fTur[i,t-1]=1
		end
	end
end
#no testear
fzero=zeros(N,T)

########################################################################################################################################################################
#CAmbiar estos N,T,S
S=1	
N=100
T=100
m=10
Trabajo=ones(2*N,T)
Frecuencia0=zeros(N+m,T)
Frecuencia1=ones(N+m,T)
Frecuencia10=zeros(N+m,T)
for t in 1:T
	if mod(t-1,10)==0
		Frecuencia10[:,t]=ones(N+m)
	end
end
Frecuencia102=zeros(N+m,T)
for t in 1:T
	if mod(t-1,5)==0
		Frecuencia102[:,t]=ones(N+m)
	end
end
Frecuenciar=zeros(N+m,T)
for t in 1:T
	for i in 1:10
		if mod(t-1,10)==i-1
			Frecuenciar[(i-1)*Int(N/m)+1:(i)*Int(N/m),t]=ones(Int(N/m))
			Frecuenciar[N+i,t]=1	
		end
	end
end
Frecuenciar2=zeros(N+m,T)
for t in 1:T
	for i in 1:5
		if mod(t-1,5)==i-1
			Frecuenciar2[2*(i-1)*Int(N/m)+1:(i)*2*Int(N/m),t]=ones(2*Int(N/m))
			Frecuenciar2[N+2*i-1,t]=1	
			Frecuenciar2[N+2*i,t]=1	

		end
	
	end
	
end


 #Creacion de un grupo 
    for a=1:S
	global Group1
	if a==1
		g1=[zeros(Int((a-1)*(N/S)))' ones(Int(N/S))' zeros(Int(N-a*(N/S)))']
		g2="none"
		Group1=[g1 g2]
	else

		g1=[zeros(Int((a-1)*(N/S)))' ones(Int(N/S))' zeros(Int(N-a*(N/S)))']
		g2="none"
		g=[g1 g2]
		Group1=[Group1 g]
	end
    end
#Horarios e intereracciones funcionarios
MFunc=zeros(2*(N+m),m,T)
MFunc2=zeros(2*(N+m),m,T)
for i in 1:10
	MFunc2[(i-1)*20+1:i*20,i,:]=ones(20,T)
end
MFunc3=zeros(2*(N+m),m,T)
for i in 1:10
	MFunc3[(i-1)*20+1:i*20,i,:]=ones(20,T)
end
MFunc3[2*N+1:2*(N+m),:,:]=ones(2*m,m,T)

###############################################################################
z=14
#Caso cuarentena preventiva antes del trabajo en jornada de 2-12

base_mQua,base_mQua2, base_mInf,base_mInf2,base_mInf3, base_mNFp, base_T, baseg_t, base_mSy,base_mFQua,base_mFInf,base_mFInf2,base_mFSy= simulation(N, Frecuencia0,0,T,β,γ,p_false_positive,R, false,p,Group1,Trabajo,"no","solo",z,t_peak,0,MFunc3)
ideal_mQua,ideal_mQua2, ideal_mInf,ideal_mInf2,ideal_mInf3, ideal_mNFp, ideal_T, idealg_t, ideal_mSy,ideal_mFQua,ideal_mFInf,ideal_mFInf2,ideal_mFSy = simulation(N, Frecuencia1,N,T,β,γ,p_false_positive,R, true,p,Group1,Trabajo,"no","solo",z,t_peak,0,MFunc3)
#ideal2_mQua,ideal2_mQua2, ideal2_mInf,ideal2_mInf2,ideal2_mInf3, ideal2_mNFp, ideal2_T = simulation(N, fr44,N,T,β,γ,p_false_positive,R, true,p,Group1,Trab4,"no","solo",z,t_peak,MFunc2)
f1_mQua,f1_mQua2, f1_mInf,f1_mInf2,f1_mInf3, f1_mNFp, f1_T, f1g_t, f1_mSy,f1_mFQua,f1_mFInf,f1_mFInf2,f1_mFSy = simulation(N, Frecuencia10,N,T,β,γ,p_false_positive,R, true,p,Group1,Trabajo,"no","solo",z,t_peak,0,MFunc3)
fr_mQua,fr_mQua2, fr_mInf,fr_mInf2,fr_mInf3, fr_mNFp, fr_T, frg_t, fr_mSy,fr_mFQua,fr_mFInf,fr_mFInf2,fr_mFSy = simulation(N, Frecuenciar,N,T,β,γ,p_false_positive,R, true,p,Group1,Trabajo,"no","solo",z,t_peak,0,MFunc3)



#f1 =plot(1:T-1,base_mQua[1:end-1],label ="No testear", lw=3)
#plot!(1:T-1,ideal_mQua[1:end-1],label ="Simple, Todos los días", lw=3)
#plot!(1:T-1,ideal2_mQua[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
#plot!(1:T-1,f1_mQua[1:end-1],label ="Simple, Cada 10 días", lw=3)
#title!("Cuarentena Trabajando")
#xlabel!("Días")
#ylabel!("# % personas")

f12 =plot(1:T-1,base_mQua2[1:end-1],label ="No testear", lw=3)
plot!(1:T-1,ideal_mQua2[1:end-1],label ="Simple, Todos los dias", lw=3)
#plot!(1:T-1,ideal2_mQua2[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
plot!(1:T-1,f1_mQua2[1:end-1],label ="Simple, Cada 10 dias", lw=3)
plot!(1:T-1,fr_mQua2[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
title!("Personas en Cuarentena, con N=100")
xlabel!("Dias")
ylabel!("# personas")



f2 =plot(1:T-1,base_mInf[1:end-1],label ="No testear", lw=3)
plot!(1:T-1,ideal_mInf[1:end-1],label ="Simple, Todos los dias", lw=3)
#plot!(1:T-1,ideal2_mInf[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
plot!(1:T-1,f1_mInf[1:end-1],label ="Simple, Cada 10 dias", lw=3)
plot!(1:T-1,fr_mInf[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
title!("Personas infectadas no descubiertas, con N=100")
xlabel!("Días")
ylabel!("# personas")

#f22 =plot(1:T-1,base_mInf2[1:end-1],label ="No testear", lw=3)
#plot!(1:T-1,ideal_mInf2[1:end-1],label ="Simple, Todos los dias", lw=3)
#plot!(1:T-1,ideal2_mInf2[1:end-1],label ="Simple, Cada da los 4 grupos ver2", lw=3)
#plot!(1:T-1,f1_mInf2[1:end-1],label ="Simple, Cada 10 dias", lw=3)
#title!()
#xlabel!("Días")
#ylabel!("# personas")

f23 =plot(1:T-1,base_mInf3[1:end-1],label ="No testear", lw=3)
plot!(1:T-1,ideal_mInf3[1:end-1],label ="Simple, Todos los dias", lw=3)
#plot!(1:T-1,ideal2_mInf3[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
plot!(1:T-1,f1_mInf3[1:end-1],label ="Simple, Cada 10 dias", lw=3)
plot!(1:T-1,fr_mInf3[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
title!("Personas Infectadas, con N=100")
xlabel!("Días")
ylabel!("# personas")


f3 =plot(1:T-1,base_T[1:end-1],label ="No testear", lw=3)
plot!(1:T-1,ideal_T[1:end-1],label ="Simple, Todos los dias", lw=3)
#plot!(1:T-1,ideal2_T[1:end-1],label ="Simple, Cada dia los 4 grupos ver2", lw=3)
plot!(1:T-1,f1_T[1:end-1],label ="Simple, Cada 10 dias", lw=3)
plot!(1:T-1,fr_T[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
title!("Test Realizados")
xlabel!("Días")
ylabel!("# Test")


f4 =plot(1:T-1,base_mNFp[1:end-1],label ="base", lw=3)
plot!(1:T-1,ideal_mNFp[1:end-1],label ="Simple, Todos los dias", lw=3)
#plot!(1:T-1,ideal2_mNFp[1:end-1],label ="Simple, Cada dia los 4 grupos ver2", lw=3)
plot!(1:T-1,f1_mNFp[1:end-1],label ="Simple, Cada 10 dias", lw=3)
plot!(1:T-1,fr_mNFp[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
title!("Falsos Positivos")
xlabel!("Días")
ylabel!("# Test")


f5 =plot(1:T-1,base_mSy[1:end-1],label ="base", lw=3)
plot!(1:T-1,ideal_mSy[1:end-1],label ="Simple, Todos los dias", lw=3)
plot!(1:T-1,f1_mSy[1:end-1],label ="Simple, Cada 10 dias", lw=3)
plot!(1:T-1,fr_mSy[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
title!("Personas con Sintomas, con N=100")
xlabel!("Días")
ylabel!("# Personas")



f6 =plot(1:T-1,base_mFQua[1:end-1],label ="No testear", lw=3)
plot!(1:T-1,ideal_mFQua[1:end-1],label ="Simple, Todos los dias", lw=3)
#plot!(1:T-1,ideal2_mFQua[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
plot!(1:T-1,f1_mFQua[1:end-1],label ="Simple, Cada 10 dias", lw=3)
plot!(1:T-1,fr_mFQua[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
title!("Funcionarios en Cuarentena, con m=10")
xlabel!("Dias")
ylabel!("# personas")

f7 =plot(1:T-1,base_mFSy[1:end-1],label ="base", lw=3)
plot!(1:T-1,ideal_mFSy[1:end-1],label ="Simple, Todos los dias", lw=3)
plot!(1:T-1,f1_mFSy[1:end-1],label ="Simple, Cada 10 dias", lw=3)
plot!(1:T-1,fr_mFSy[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
title!("Funcionarios con Sintomas, con m=10")
xlabel!("Días")
ylabel!("# Personas")


f8 =plot(1:T-1,base_mFInf[1:end-1],label ="No testear", lw=3)
plot!(1:T-1,ideal_mFInf[1:end-1],label ="Simple, Todos los dias", lw=3)
#plot!(1:T-1,ideal2_mFInf[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
plot!(1:T-1,f1_mFInf[1:end-1],label ="Simple, Cada 10 dias", lw=3)
plot!(1:T-1,fr_mFInf[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
title!("Funcionarios infectados trabajando no descubiertos, con m=10")
xlabel!("Días")
ylabel!("# personas")


f9 =plot(1:T-1,base_mFInf2[1:end-1],label ="No testear", lw=3)
plot!(1:T-1,ideal_mFInf2[1:end-1],label ="Simple, Todos los dias", lw=3)
#plot!(1:T-1,ideal2_mFInf2[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
plot!(1:T-1,f1_mFInf2[1:end-1],label ="Simple, Cada 10 dias", lw=3)
plot!(1:T-1,fr_mFInf2[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
title!("Funcionarios infectados, con m=10")
xlabel!("Días")
ylabel!("# personas")

print("Caso: Falsos Positivos - Test "," No testear: ", round(sum(base_mNFp)),"-",round(sum(base_T))," Simple, Cada día los 4 grupos: " ,round(sum(ideal_mNFp)),"-",round(sum(ideal_T))," Simple, Cada día 2 grupos: ",round(sum(f1_mNFp)),"-",round(sum(f1_T)))

plot(f12,f5, layout = (2,1))

savefig("CuarentenaySintomas.pdf")
plot(f6,f7, layout = (2,1))

savefig("CuarentenaySintomasFuncionarios.pdf")


plot(f2,f23, layout = (2,1))
savefig("Infectados.pdf")

plot(f8,f9, layout = (2,1))
savefig("InfectadosFuncionarios.pdf")


plot(f3,f4, layout = (2,1))
savefig("Test.pdf")




