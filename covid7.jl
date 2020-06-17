using Distributions, Plots, LinearAlgebra, Random, StatsPlots

#######################################################################################
#Recibe un grupo y entrega un vector con los indices de los grupos a los que pertenecen
function vecgroup(group)
S=Int(length(group)/(N+1))
A=zeros(N,S)
v=zeros(N)
for b in 1:S
		A[:,b]=group[Int.((N+1)*(b-1)+1:N)]
end
for i in 1:N
v[i]=Int.(setdiff((A[i,:]).*collect(1:S),[0]))[1]


end
return v
end
########################################################################


########################################################################
#############Funcion que crea los tiempos###############################
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
#######################################################################




#####################################################################################################
#########################Funcion para la aproximacion de la funcion de test positivo##################
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
#################################################################################################



#################################################################################################
##################Mejora a la funcion de aproximar###############################################
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
##############################################################################################


#############################################################################################
############################Funcion acumulado################################################
####Recibe un vector con valor diario y entrega el acumulado de los dias anteriores#########
function acu(vec)
	#vec: vector
	#Recibe un vector y entrega el acumulado
	lar=length(vec)
	vec2=zeros(lar)
	vec2[1]=vec[1]
	for i in 2:lar
		vec2[i]=vec[i]+vec2[i-1]
	end
return vec2
end

#############################################################################################




#############################################################################################
##################Funcion de test positivo####################################################
function true_positive(t, t_inf, t_sin, t_rec,asin,scalar_asint)
    # t : dia relativo con respecto a la expocision, tb puede ser un vector
	#asin: Si es asintomatico
	#scalar_asint: Cuanto pondera a la curva de los asintomaticos
    # vectores de valores conocidos en eje x e y
    x_vec = vcat([t_inf,(t_inf + t_sin) / 4,3*(t_inf + t_sin) / 4, t_sin], [t_sin + (t_rec - t_sin) * i / 28 for i in 7:7:42])
    #x_vec = [-7, -3.5,   0,    7,   14,   21,   28,   35,   42]
	if asin==0
		y_vec = [ 0, 0.63, 0.8, 0.79, 0.64, 0.43, 0.19, 0.08, 0.00,0.00]
	else
		y_vec = [ 0, 0.63, 0.8, 0.79, 0.64, 0.43, 0.19, 0.08, 0.00,0.00]*scalar_asint
	end

    y = cubic_spline_interpolate(t, x_vec, y_vec)

    # graficar
    if false #true
        x_plot = [x_vec[1]:x_vec[end]]
        y_plot = cubic_spline_interpolate(x_plot, x_vec, y_vec)
        plot_interpolate(x_plot, y_plot, x_vec, y_vec)
    end

    return y
end
##############################################################################################




##################################################################################################
####################################Funcion de probabilidad de contagio uno a uno###################
#Recibe los tiempos y una probabilidad peak
function person_p(t,t_inf,t_sin,t_rec,p_sinto,asin,scalar_asint)
    # t : dia relativo con respecto a la expocision, tb puede ser un vector
    # vectores de valores conocidos en eje x e y
	#asin: Si es asintomatico
	#scalar_asint: Cuanto pondera a la curva de los asintomaticos
    x_vec = vcat([t_inf, (t_inf + t_sin) / 2, t_sin], [t_sin + (t_rec - t_sin) * i / 28 for i in 7:7:42])
    #x_vec = [-7, -3.5,   0,    7,   14,   21,   28,   35,   42]
	if asin==0
		y_vec = [ 0, 0.63, 0.8, 0.79, 0.64, 0.43, 0.19, 0.08, 0.00] * (p_sinto/0.8)
	else
		y_vec = [ 0, 0.63, 0.8, 0.79, 0.64, 0.43, 0.19, 0.08, 0.00] * (p_sinto/0.8)*scalar_asint
	end
    y = cubic_spline_interpolate(t, x_vec, y_vec)
    return y
end
###################################################################################################




##############################################################################################
##################Funcion que crea una exponencial para los tiempos de contacto con un contagiado, en base a la
##################funcion anterior y un tiempo peak de contagio.

#Recibiendo los parametros de person_p, un número de horas de turno (exposición al contagio) y un t_peak, mediante una exponencial
#estimando lambda como el que haga que en t_peak se alcanze P(tiempo para contagiarse<=t_peak)=probabilidad peak de contagio
#obtiene la probabilidad de haberse contagiado en ese turno
function exponential_p(t,t_inf,t_sin,t_rec,p_sinto,t_peak,horas_turno,asin,scalar_asint)
	prob=person_p(t,t_inf,t_sin,t_rec,p_sinto,asin,scalar_asint) #Se calcula el peak que se alcanza en t_peak
	lambda=(log(max(1-prob,0.00001)))/t_peak #Se calcula el lambda tal que P(tiempo para contagiarse<=t_peak)=probabilidad peak de contagio
	proba=1-exp(lambda*horas_turno) #Se calula P(tiempo para contagiarse<=horas_turno)
	return proba
end
##############################################################################




################################################################################
#####################Se usa para calcular las probabilidades de contagiar a cada persona de un grupo, entrega un vector con las probabilidades de
#####################cada integrante del grupo y 0 en los que no son del grupo
#Funcion que entrega la probabilidad de contagiarte dentro de tu grupo de trabajo cerrado
function group_p(G,M_cont,W)
	#G: indicatriz del grupo de trabajo
	#vp_cont: vector con las probabilidades de contagio
	#W: Peligro de contagio de paciente
	V=((M_cont).*G)
	N=length(M_cont[1,:])
	p=zeros(N)
	for i in 1:N
		p[i]=(1-prod(-V[:,i].+1))
	end
	#de no contagiarse multiplicado por un factor que pondera las horas, este desde
	#la hora toma valores del estilo ~0.9 y mas.
#Lo siguiente corresponde al riesgo del trabajo, para los abuelitos consideramos riesgo 0
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

#########################################################################################






#########################################################################################
#Esta es la funcion principal y entrega
#mQua: Personas en cuarentena que deberian estar trabajando.
#mQua2: Personas en cuarentena en total.
#mInf: Personas infectadas trabajando (peligrosas, son las no descubiertas)
#mInf2: Personas infectadas trabajando o mandadas a cuarentena (se agregan al anterior las que estan en cuarentena, pero les corresponde el turno, son las perdidas)
#mInf3: Personas infectadas en total
#mSy: Personas con sintomas.
#mFQua: Funcionarios en cuarentena en total.
#mFInf: Funcionarios infectados trabajando (peligrosas, son los no descubiertos)
#mFInf2: Funcionarios infectados en total
#mFSy: Funcionarios con sintomas.
#maxInf: Curva con el mayor número de infectados en el peak
#maxQua: Curva cuarentena con el mayor número de infectados en el peak
#maxSy: Curva sintomas con el mayor número de infectados en el peak
#Infect: Mayor número de infectados en el peak, por cada repeticion
#VInf1: Numero de abuelitos infectados en promedio, por dia
#VInf2: Numero de funcionarios infectados en promedio, por dia
function simulation(N, f,G,T,peak_A,pA_int,pA_ext,peak_B,pB_int,pB_ext,γ,p_false_positive,R,test,Group,Tra,random,cuarentena,z,t_peak,dias_atras,MF,scalar_asint)

#N Tamaño grupo de trabajadores
#f Frecuencia testeo (puede ser vector de largo T, en caso de testear random de tamaño G o matriz de NxT con los individuos a testear en sus fechas)
# G Solo aplica para el caso random el número de test individuales o poll (por ejemplo G=4 son 4 poll testing y todos los sub test que se generen)
# T simulation horizon
#
##prob Prob peak de infeccion trbajadores externos/funcionarios en asilo

##Visitas: probabilidad contagio visita
#peak_A: Probabilidad contagio entre pacientes/trabajdores en el caso hospital

#
#
#
#γ  prop of asymtomatic pop
#p_false_positive probability of a false positive
#R Replications
#test dice si es "poll", t(individual) o false (no hace).
#Group un vector con indicatrices de tamaño N con los integrantes de cada grupo y riesgo de atención, tamaño (N+1)
#Tra: Matriz (2*N,T) de turnos de trabajo
#random: Si es random o estan fijos a quienes se testeara (no los dias esos estan en f, solo los funcionarios a testear)
#cuarentena: Si se mandan a todo el grupo en cuarentena o solo al infectado (ahora tambien mixto, solo al grupo cuando uno presenta sintomas)
#z: dias de cuarntena al grupo, no se usa en caso de solo mandar al positivo
#t_peak: tiempo en que la exponencial alcanza la probabilidad peak
#dias_atras: cuantos dias hacia atras se mira con quienes has tenido contacto estrecho
#MF: Matriz de funcionarios, contiene sus interacciones con los pacientes y entre ellos.
hr=12
if length(MF)>=1
	m=length(MF[1,:,1])
else
	m=0
end
MD=zeros(N,T) #Matriz trabajo dia
MN=zeros(N,T) ##Matriz trabajo noche
for i in 1:N
	MD[i,:]=Tra[2*i-1,:]
	MN[i,:]=Tra[2*i,:]
end
#Se crean los trabajadores del tipo B de dia y de noche
if m>=1
	MDF=zeros(N,m,T) #Matriz trabajo dia
	MNF=zeros(N,m,T) ##Matriz trabajo noche
	for i in 1:N
		MNF[i,:,:]=MF[2*i-1,:,:]
		MDF[i,:,:]=MF[2*i,:,:]
	end
else
	MDF=[]
	MNF=[]
end
S=length(Group)/(N+1) #Tamaño grupo
#Cuando se utiliza la cuarentena, con este vector se recuperan los integrantes del grupo para mandarlos a cuarentena
if (cuarentena=="grupo")|(cuarentena=="mixto")
	v=vecgroup(Group)
end


#Cuenta cuantas personas distintas estan trabajando por día
g_t=ones(T)
	for tt in 1:T
		if sum((MD[:,tt]+MN[:,tt]).>=1)>=1
			g_t[tt]=sum((MD[:,tt]+MN[:,tt]).>=1)
		end
	end




#Vectores que guardan los valores de cuarentena, infectados, sintomas y números de test en cada t
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
Gf=[zeros(N)' ones(m)']'
maxInf=zeros(T)
maxQua=zeros(T)
maxSy=zeros(T)
Infect=ones(R)*1000
VInf1=zeros(T)
VInf2=zeros(T)


    # Loop over replications
    for rep = 1:R
	print(rep) #Lo uso para saber si corre
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
	NTest=zeros(T) #Vector de test diarios
	NFpositive=zeros(T) #Falsos positivos diarios
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
		for a=1:S #Se recorren los grupos
			info=Group[Int((a-1)*(N+1)+1):Int((a)*(N+1))] #Se recupera la info codificada del grupo
			g=[Int.(info[1:N])' zeros(m)' ]' #Indicatriz del grupo
			W=info[N+1] #Riesgo del grupo
			vec=zeros(N) #Vector auxiliar para las probabilidades de los infectados
			#del grupo

########################################################################################
			for i in findall((g.*inf.*(-qu.+1).*(-re.+1)).==1)#Se recorren los abuelitos enfermos, no en cuarentena y no recuperados
				for j in findall((g.*su.*(-qu.+1)).==1) #Se recorren los abuelitos que nunca se han contagiado y no en cuarentena
					if (MD[i,t]*MD[j,t]==1) & ((1-MN[i,t])*(1-MN[j,t])==1)  #Condicion esten ambos trabajando
						#Recuperacion de los tiempos
						tinf=tIn[Int(i)]-tEx[Int(i)]
						tsin=tSy[i]-tEx[i]
						trec=tRe[i]-tEx[i]
						MC[i,j]=exponential_p(t,tinf,tsin,trec,peak_A,t_peak,hr,As[i],scalar_asint)
				 		#Se agregan las probabilidades de los enfermos
					end
					if ((1-MD[i,t])*(1-MD[j,t])==1) & (MN[i,t]*MN[j,t]==1)
						#Recuperacion de los tiempos
						tinf=tIn[Int(i)]-tEx[Int(i)]
						tsin=tSy[i]-tEx[i]
						trec=tRe[i]-tEx[i]
						MC[i,j]=exponential_p(t,tinf,tsin,trec,peak_A,t_peak,hr,As[i],scalar_asint)
				 		#Se agregan las probabilidades de los enfermos
					end
					if (MD[i,t]*MD[j,t]==1) & (MN[i,t]*MN[j,t]==1)
						#Recuperacion de los tiempos
						tinf=tIn[Int(i)]-tEx[Int(i)]
						tsin=tSy[i]-tEx[i]
						trec=tRe[i]-tEx[i]
						MC[i,j]=exponential_p(t,tinf,tsin,trec,peak_A,t_peak,2*hr,As[i],scalar_asint)
						#Se agregan las probabilidades de los enfermos
					end
				end
				if m>=1
				for f in findall((su[N+1:N+m].*(-qu[N+1:N+m].+1)).==1) #Se recorren los funcionarios sanos
					if (MF[2*i-1,f,t]==1) & (MF[2*i,f,t]==0)  #Condicion esten trabajando el mismo dia
						#Recuperacion de los tiempos
						tinf=tIn[Int(i)]-tEx[Int(i)]
						tsin=tSy[i]-tEx[i]
						trec=tRe[i]-tEx[i]
						MC[i,N+f]=exponential_p(t,tinf,tsin,trec,peak_A,t_peak,hr,As[i],scalar_asint)
				 		#Se agregan las probabilidades de los enfermos
					end
					if (MF[2*i-1,f,t]==0) & (MF[2*i,f,t]==1)
						#Recuperacion de los tiempos
						tinf=tIn[Int(i)]-tEx[Int(i)]
						tsin=tSy[i]-tEx[i]
						trec=tRe[i]-tEx[i]
						MC[i,N+f]=exponential_p(t,tinf,tsin,trec,peak_A,t_peak,hr,As[i],scalar_asint)
				 		#Se agregan las probabilidades de los enfermos
					end
					if (MF[2*i-1,f,t]==1) & (MF[2*i,f,t]==1)
						#Recuperacion de los tiempos
						tinf=tIn[Int(i)]-tEx[Int(i)]
						tsin=tSy[i]-tEx[i]
						trec=tRe[i]-tEx[i]
						MC[i,N+f]=exponential_p(t,tinf,tsin,trec,peak_A,t_peak,2*hr,As[i],scalar_asint)
				 		#Se agregan las probabilidades de los enfermos
					end
				end
				end
			end
			if m>=1
			for f in findall((inf[N+1:N+m].*(-qu[N+1:N+m].+1).*(-re[N+1:N+m].+1)).==1)#Se recorren los funcionarios enfermos, no en cuarentena y no recuperados
				for j in findall((g.*su.*(-qu.+1)).==1) #Se recorren los abuelitos enfermos
					if (MF[2*j-1,f,t]==1) & (MF[2*j,f,t]==0)
						#Recuperacion de los tiempos
						tinf=tIn[Int(N+f)]-tEx[Int(N+f)]
						tsin=tSy[N+f]-tEx[N+f]
						trec=tRe[N+f]-tEx[N+f]
						MC[N+f,j]=exponential_p(t,tinf,tsin,trec,peak_B,t_peak,hr,As[N+f],scalar_asint)
				 		#Se agregan las probabilidades de los enfermos
					end
					if (MF[2*j-1,f,t]==0) & (MF[2*j,f,t]==1)
						#Recuperacion de los tiempos
						tinf=tIn[Int(N+f)]-tEx[Int(N+f)]
						tsin=tSy[N+f]-tEx[N+f]
						trec=tRe[N+f]-tEx[N+f]
						MC[N+f,j]=exponential_p(t,tinf,tsin,trec,peak_B,t_peak,hr,As[N+f],scalar_asint)
				 		#Se agregan las probabilidades de los enfermos
					end
					if (MF[2*j-1,f,t]==1) & (MF[2*j,f,t]==1)
						#Recuperacion de los tiempos
						tinf=tIn[Int(N+f)]-tEx[Int(N+f)]
						tsin=tSy[N+f]-tEx[N+f]
						trec=tRe[N+f]-tEx[N+f]
						MC[N+f,j]=exponential_p(t,tinf,tsin,trec,peak_B,t_peak,2*hr,As[N+f],scalar_asint)
				 		#Se agregan las probabilidades de los enfermos
					end
				end

			end
			end
			PG=[group_p(g+Gf,MC[:,1:N],W)' zeros(m)']' #Vector con probabilidades de contagio dentro del grupo
			r=rand(N+m)
			mm=zeros(N+m)
			mm[1:N]=(MD[:,t]+MN[:,t])/2
			new_inf+=(-qu.+1).*su.*g.*(r.<(-(-PG.+1).*(-mm.*pA_int.+1).*(-(-mm.+1).*pA_ext.+1).+1)) #Vector de infectados del grupo que no esten en cuarentena, esten o no trabajando
			new_inf+=qu.*su.*g.*(r.<(pA_ext)) #Sumar infectados en cuarentena (falsos positivos)

#################################################################################################################
			#Para el poll test, predefino el numero de grupos que se pueden formar
			#se ve quienes son el publico a testear
			if test=="poll"
				if (random=="si") & (f[t]==1)
					I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*(-qu.+1).*(tra).*collect(1:N+m) #Gente que esta trabajando, no en cuarentena y nunca a presentado sintomas
					JJ=shuffle!(g.*I) #Se mezclan
					III=setdiff(JJ,[0]) #Se eliminan los 0 Se cuentan
					if length(III)==0
						con+=0
					else
						ind=Int(round(length(III)/5))
						r=length(III)%5
						con+=(ind+(r!=0))
					end
				end
				#Para el poll test no random, es similar. Se hacen distintos, pues en caso random la frecuencia es un vector
				if (random=="no")& (sum(f[:,t],dims=1)==1)
					I2=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*f[:,t].*(-qu.+1).*collect(1:N+m)
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

####################################################################################################################
	   	end

###################################################################################################################
############################Infeccion entre funcionarios###########################################################
		if m>=1
		for f in findall((inf[N+1:N+m].*(-qu[N+1:N+m].+1).*(-re[N+1:N+m].+1)).==1) #funcionarios infectados, no en cuarentena y no recuperados
			for f2 in findall((su[N+1:N+m].*(-qu[N+1:N+m].+1)).==1) #Funcionarios sanos y no en cuarentena
					if (MF[2*f2-1,f,t]==1) & (MF[2*f2,f,t]==0)  #Que coincidan en los turnos
						#Recuperacion de los tiempos
						tinf=tIn[Int(N+f)]-tEx[Int(N+f)]
						tsin=tSy[N+f]-tEx[N+f]
						trec=tRe[N+f]-tEx[N+f]
						MC[N+f,N+f2]=exponential_p(t,tinf,tsin,trec,peak_B,t_peak,hr,As[N+f],scalar_asint)
				 		#Se agregan las probabilidades de los enfermos
					end
					if (MF[2*f2-1,f,t]==0) & (MF[2*f2,f,t]==1)
						#Recuperacion de los tiempos
						tinf=tIn[Int(N+f)]-tEx[Int(N+f)]
						tsin=tSy[N+f]-tEx[N+f]
						trec=tRe[N+f]-tEx[N+f]
						MC[N+f,N+f2]=exponential_p(t,tinf,tsin,trec,peak_B,t_peak,hr,As[N+f],scalar_asint)
				 		#Se agregan las probabilidades de los enfermos
					end
					if (MF[2*f2-1,f,t]==1) & (MF[2*f2,f,t]==1)
						#Recuperacion de los tiempos
						tinf=tIn[Int(N+f)]-tEx[Int(N+f)]
						tsin=tSy[N+f]-tEx[N+f]
						trec=tRe[N+f]-tEx[N+f]
						MC[N+f,N+f2]=exponential_p(t,tinf,tsin,trec,peak_B,t_peak,2*hr,As[N+f],scalar_asint)
				 		#Se agregan las probabilidades de los enfermos
					end
				end
		end
		end
	PGf=-prod((-MC.*Gf).+1,dims=1).+1 #Vector de probabilidades de contagio funcionarios
	r=rand(N+m)
	new_inf+=(-qu.+1).*su.*Gf.*(r.<(-(-PGf'.+1).*(1-pB_ext).+1)) #Vector de infectados de los funcionarios que no esten en cuarentena, esten o no trabajando
	new_inf+=qu.*su.*Gf.*(r.<(pB_ext)) #Sumar funcionarios infectados en cuarentena (falsos positivos)
#################################################################################################################################

##########################################################################################
#Se realiza el proceso de cambio en los estados de los nuevos infectados
           	for j in findall(new_inf.>=1)                   # cicle over newly infected
			i=j[1]                                                  #recorre los indices
			if i<=N #Infectados nuevos por dia
				VInf1[t]+=1
			else
				VInf2[t]+=1
			end
               		t_inc, ti_inf, te_inf, t_rec = get_inc()                    #los tiempos para cada nuevo infectado
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

		if (t>=3)&((cuarentena=="mixto")|(cuarentena=="grupo"))              #Si se activa se manda a todo el grupo cuando alguien tiende sintomas
			for i in findall(Sy[:,t].*(Qu[:,t]).*Sy[:,t-1].*(-Qu[:,t-1].+1).*Sy[:,t-2].*(-Qu[:,t-2].+1).==1) #Condicion de no ser pillado hasta los sintomas
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


##########################Testing###################################################
		#Vector de testear
	    	cand=[]
###################################Caso random, se a usado, poco es cuando solo se sabe que dias y a cuantos se quiere testear,
###################################con estos numeros se saca una muestra de minimo entre tamaño G y los que trabajan al azar. De momento esta descontinuado
		if random=="si" #Se realiza test al azar
			if f[t]==1 #Se testea estos dias
	    			if test=="t" #test individual
					I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*(tra).*(-qu.+1).*collect(1:N+m) #Se ven los candidatos a testear
					J=shuffle!(I) #Mezclar
					II=setdiff(J,[0]) #Eliminar ceros
					cand=Int.(II[1:min(G,length(II))])  #Se toman como candidatos a los posibles y el minimo entre ellos y G
	    			end

				if test!="poll"	#Caso individual
					for i in cand	#se recorren los factibles a testear
						NTest[t]+=1 #Se hace un test
                				if ((su[i] == 1)|((su[i]==0)&(re[i]==1) ))& (rand() < p_false_positive) #Falso positivo
                    					Qu[i,min(T,t+1):min(T,t+15)] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
										if cuarentena=="grupo"
											for y in findall((-qu.+1).*(Group[Int((N+1)*(v[i]-1)+1):Int((N+1)*(v[i]-1)+N)]).==1)
												Qu[y,min(T,t+1):min(T,t+1+z)] .= 1 #quarantined for 2 weeks
											end
										end
										NFpositive[t]+=1
                				end
                				if (su[i] == 0)&(re[i]==0) #Esta o estuvo enfermo
                        					p_positive = true_positive(t-tEx[i]+1, tIn[i]-tEx[i]+1, tSy[i]-tEx[i]+1, tRe[i]-tEx[i]+1,As[i],scalar_asint)
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
#######################Caso no random ,usa la matriz de testeo y testea a todas esas personas, exceptuando a los que estan en cuarentena o
####################### ya tuvieron sintomas, por lo que son inmunes
		if random=="no" #Se realiza el esquema de test de manera similar solo que los candidatos ya estan dados por defecto y solo saca los que esten en cuarentena o ya hayan tenido sintomas en el pasado
			if sum( f[:,t].*(-qu.+1))>=1
	  	  		if test=="t"
					testear=
					J=Int.(( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*f[:,t].*(-qu.+1).*collect(1:N+m)) #Son los que nunca han tenido sintomas, no cuarentena y toca testear
					II=setdiff(J,[0])
					cand=Int.(II)  #Se obtiene vector de a quienes testear
	    		end
				if test=="poll" #Esta medio obsoleto
					h=collect(1:con2)[1:min(G,con2)]
					b=sort(h)
					for a=1:S
						info=Group[Int((a-1)*(N+1)+1):Int((a)*(N+1))]
						g=Int.(info[1:N]) #Indicatriz del grupo
						I=(-qu.+1).*collect(1:N).*f[:,t]
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
									pp=zeros(N)
									for i in cand

										p_positive = true_positive(t-tEx[i]+1, tIn[i]-tEx[i]+1, tSy[i]-tEx[i]+1, tRe[i]-tEx[i]+1,As[i],scalar_asint)
										s+=su[i]
										pp[i]=p_positive
									end
									rp=(rand(N).<=pp)
									for rpi in findall(rp.==0)
										if (((su[rpi]==0)&(re[rpi]==1))|(su[rpi]==1))&(rand()<p_false_positive)&(rpi in cand)
											rp[rpi]=1
										end
									end
									suma=0
									for i in cand
										suma+=rp[i]
									end
									if suma==0
										NTest[t]+=1
									else
										NTest[t]+=min(length(cand)-1,1)+length(cand)
									end
									for i in findall(rp.==1)
											if (((As[i]==1)&(su[i]==0)&(re[i]==1))|(su[i]==1))&(i in cand)
                	    						Qu[i,min(T,t+1):min(T,t+15)] .= 1 #if someone tests (false) positive, quarantined for 2 weeks
												if cuarentena=="grupo"
													for y in findall((-qu.+1).*(Group[Int((N+1)*(v[i]-1)+1):Int((N+1)*(v[i]-1)+N)]).==1)
														Qu[y,min(T,t+1):min(T,t+1+z)] .= 1 #quarantined for 2 weeks
													end
												end
												NFpositive[t]+=1
											end
											if (su[i]==0)&(re[i]==0)&(i in cand)

           												Qu[i,Int(min(T,t+1)):Int(min(T,max(t+15,tRe[i]+3)))] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
														if cuarentena=="grupo"
																for y in findall((-qu.+1).*(Group[Int((N+2)*(v[i]-1)+1):Int((N+2)*(v[i]-1)+N)]).==1)
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
				if (test=="t")	& (length(cand)>=1)	#Si al menos hay 1 para testear y no poll
					for i in cand	#se recorren
						NTest[t]+=1   #Un test nuevo
           				     		if (((As[i]==1)&(su[i]==0)&(re[i]==1))|(su[i]==1))& (rand() < p_false_positive) #Caso falso positivo
			                    		Qu[i,min(T,t+1):min(T,t+15)] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
										if cuarentena=="grupo"
											for y in findall((-qu.+1).*(Group[Int((N+1)*(v[i]-1)+1):Int((N+1)*(v[i]-1)+N)]).==1)
												Qu[y,min(T,t+1):min(T,t+1+z)] .= 1 #quarantined for 2 weeks
											end
										end
										NFpositive[t]+=1
                					end
                				if (su[i]==0)&(re[i]==0)

                        					p_positive = true_positive(t-tEx[i]+1, tIn[i]-tEx[i]+1, tSy[i]-tEx[i]+1, tRe[i]-tEx[i]+1,As[i],scalar_asint)
            			            		if rand() < p_positive #Caso positivo
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
	#Se guardan los valores en cada
        mQua .+= dropdims(sum(   (((MD).*Qu[1:N,:]+(MN).*Qu[1:N,:]).>=1)     ,dims=1),dims=1)
	mQua2 .+= dropdims(sum(Qu[1:N,:],dims=1),dims=1)
        mInf .+= dropdims(sum(      (((MD).*In[1:N,:].*(-Qu[1:N,:].+1)+(MN).*In[1:N,:].*(-Qu[1:N,:].+1)).>=1)          ,dims=1),dims=1)
        mInf2 .+= dropdims(sum(    (((MD).*In[1:N,:]+(MN).*In[1:N,:]).>=1)          ,dims=1),dims=1)
	mInf3 .+= dropdims(sum(In[1:N,:],dims=1),dims=1)
	mSy .+= dropdims(sum(Sy[1:N,:],dims=1),dims=1)

if m>=1
	alpha=zeros(m,T)
for i in 1:m
	alpha[i,:]=(sum(MF[:,i,:],dims=1).>=1)
end
	mFQua .+= dropdims(sum(Qu[N+1:N+m,:],dims=1),dims=1)
        mFInf .+= dropdims(sum(In[N+1:N+m,:].*(-Qu[N+1:N+m,:].+1).*alpha[1:m,:],dims=1),dims=1)
        mFInf2 .+= dropdims(sum(In[N+1:N+m,:],dims=1),dims=1)
	mFSy .+= dropdims(sum(Sy[N+1:N+m,:],dims=1),dims=1)
end
	mNFpositive .+= NFpositive
	mNTest .+= NTest
	maxi=max(   maximum(maxInf)  ,  maximum(sum(In[1:N,:],dims=1))   )
	si=maximum(maxInf)==(maxi)
	maxInf=(sum(In[1:N,:],dims=1))'*(1-si)+maxInf*si
 	if maxInf==sum(In[1:N,:],dims=1)'
		maxQua=sum(Qu[1:N,:],dims=1)
		maxSy=sum(Sy[1:N,:],dims=1)
	end
	Infect[rep]=maximum(sum(In[1:N,:],dims=1))
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
	VInf1=VInf1/R
	VInf2=VInf2/R


    return mQua,mQua2,mInf,mInf2,mInf3, mNFpositive,  mNTest, g_t, mSy,mFQua,mFInf,mFInf2,mFSy, maxInf, maxQua, maxSy, Infect, VInf1, VInf2
end





activador=1

if activador==0
################################################################
######################Datos caso hospital#######################
################################################################

N = 154             # Tamaño grupo de trabajadores
T = 210              # simulation horizon
pB_ext = 0.01   #prob of contagion en casa
γ = 0.3 #(antes 0.25)           # prop of asymtomatic pop
p_false_positive = 0.01 # probability of a false positive
R = 400 # Replications
peak_B=0.2 #Prob peak de infeccion entre funcionarios y funcionario a anciano
m=0      #Cantidad de funcionarios
Trabajo=ones(2*N,T)    #Dias de trabajo abuelitos
pA_int=0.01
pA_ext=0
pB_int=0
peak_A=0.2             #probabilidad peak de contagio entre abulitos/ trbaajdores hospital
t_peak=12                 #Tiempo peak de contagio
z=14			#Dias de cuarentena preventiva
########################################################################################################################################################################
S=Int(N/22) #Cantidad de grupos
    #Creacion de un grupo de ejemplo
    for a=1:S
	global Group
	if a==1
		g1=[zeros(Int((a-1)*(N/S)))' ones(Int(N/S))' zeros(Int(N-a*(N/S)))']
		g2="none"
		Group=[g1 g2]
	else

		g1=[zeros(Int((a-1)*(N/S)))' ones(Int(N/S))' zeros(Int(N-a*(N/S)))']
		g2="none"
		g=[g1 g2]
		Group=[Group g]
	end
    end

#Frecuencias
Frecuencia0=zeros(N+m,T) #Ningun dia
Frecuencia1=ones(N+m,T) #Todos los dias
Frecuencia10=zeros(N+m,T) #Cada 10 dias
for t in 1:T
	if mod(t-1,10)==0
		Frecuencia10[:,t]=ones(N+m)
	end
end

Frecuencia5=zeros(N+m,T) #Cada 5 dias
for t in 1:T
	if mod(t-1,5)==0
		Frecuencia5[:,t]=ones(N+m)
	end
end
Frecuencia102=zeros(N+m,T) #Cada 10 dias, con funcionarios todos los dias


for t in 1:T
	if mod(t-1,10)==0
		Frecuencia102[:,t]=ones(N+m)
	end
end
if m>=1
Frecuencia102[N+1:N+m,:]=ones(m,T)
end

if m>=1
Frecuenciar=zeros(N+m,T) #Cada 10 dias, repartido
for t in 1:T
	for i in 1:10
		if mod(t-1,10)==i-1
			Frecuenciar[(i-1)*Int(N/10)+1:(i)*Int(N/10),t]=ones(Int(N/10))
			if m>=1
			Frecuenciar[N+i,t]=1
			Frecuenciar[N+i+10,t]=1
			end
		end
	end
end
Frecuenciar2=zeros(N+m,T) #Cada 5 dias repartidos, FALTA ARREGLAR
for t in 1:T
	for i in 1:5
		if mod(t-1,5)==i-1
			Frecuenciar2[2*(i-1)*Int(N/10)+1:(i)*2*Int(N/10),t]=ones(2*Int(N/10))
			if m>=1
			Frecuenciar2[N+2*i-1,t]=1
			Frecuenciar2[N+2*i,t]=1
			end

		end

	end

end

Frecuenciar3=zeros(N+m,T) #Cada 10 dias, repartido, con funcionarios todos los dias
for t in 1:T
	for i in 1:10
		if mod(t-1,10)==i-1
			Frecuenciar3[(i-1)*Int(N/10)+1:(i)*Int(N/10),t]=ones(Int(N/10))
		end
	end
end
Frecuenciar3[N+1:N+m,:]=ones(m,T)
end

################################################
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

##############################################
#Horarios e intereracciones funcionarios

MFunc=zeros(2*(N+m),m,T) #Funcionarios no trabajan
MFunc2=zeros(2*(N+m),m,T) #Funcionariso se relacionan con10 abuelitos dia por medio y no entre ellos
if m>=1
for t in 1:T
	for i in 1:10
		MFunc2[(i-1)*20+1:i*20,i+(1-mod(t,2))*10,t]=ones(20)
	end
end
MFunc3=zeros(2*(N+m),m,T) #Funcionariso se relacionan con10 abuelitos dia por medio y entre los que trabajan
for t in 1:T
	for i in 1:10
		MFunc3[(i-1)*20+1:i*20,i+(1-mod(t,2))*10,t]=ones(20)
	end
end
MFunc3[2*N+1:2*N+m,1:Int(m/2),:]=ones(m,Int(m/2),T)
MFunc3[2*N+m+1:2*N+2*m,1:Int(m/2),:]=ones(m,Int(m/2),T)
end


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




###############################################################################

politicas=3 #Numero de politicas
NombresP=["Simple, Antes de trabajar","Pool, Antes de trabajar","No testear y cuarentena"]

#Se realizan las 4 politicas
#base_mQua,base_mQua2, base_mInf,base_mInf2,base_mInf3, base_mNFp, base_T, baseg_t, base_mSy,base_mFQua,base_mFInf,base_mFInf2,base_mFSy,base_maxInf,base_maxQua, base_maxSy, base_Infect, base_VInf1, base_VInf2= simulation(N, fzero,0,T,β,γ,p_false_positive,R, "f",p,Group,Tur,"no","solo",z,t_peak,5,[],visita,peak_A,1)
ideal_mQua,ideal_mQua2, ideal_mInf,ideal_mInf2,ideal_mInf3, ideal_mNFp, ideal_T, idealg_t, ideal_mSy,ideal_mFQua,ideal_mFInf,ideal_mFInf2,ideal_mFSy,ideal_maxInf,ideal_maxQua, ideal_maxSy, ideal_Infect, ideal_VInf1, ideal_VInf2 = simulation(N, fTur,N,T,peak_A,pA_int,pA_ext,peak_B,pB_int,pB_ext,γ,p_false_positive,R,"t",Group,Tur,"no","solo",z,t_peak,5,[],1)
#f1_mQua,f1_mQua2, f1_mInf,f1_mInf2,f1_mInf3, f1_mNFp, f1_T, f1g_t, f1_mSy,f1_mFQua,f1_mFInf,f1_mFInf2,f1_mFSy,f1_maxInf,f1_maxQua, f1_maxSy,f1_Infect, f1_VInf1, f1_VInf2 = simulation(N, fTur,N,T,β,γ,p_false_positive,R, "t",p,Group,Tur,"no","grupo",z,t_peak,5,[],visita,peak_A,1)
f1_mQua,f1_mQua2, f1_mInf,f1_mInf2,f1_mInf3, f1_mNFp, f1_T, f1g_t, f1_mSy,f1_mFQua,f1_mFInf,f1_mFInf2,f1_mFSy,f1_maxInf,f1_maxQua, f1_maxSy,f1_Infect, f1_VInf1, f1_VInf2 = simulation(N, fTur,N,T,peak_A,pA_int,pA_ext,peak_B,pB_int,pB_ext,γ,p_false_positive,R, "poll",Group,Tur,"no","solo",z,t_peak,5,[],1)
fr_mQua,fr_mQua2, fr_mInf,fr_mInf2,fr_mInf3, fr_mNFp, fr_T, frg_t, fr_mSy,fr_mFQua,fr_mFInf,fr_mFInf2,fr_mFSy,fr_maxInf,fr_maxQua, fr_maxSy, fr_Infect, fr_VInf1, fr_VInf2 = simulation(N, fzero,N,T,peak_A,pA_int,pA_ext,peak_B,pB_int,pB_ext,γ,p_false_positive,R,"f",Group,Tur,"no","mixto",z,t_peak,5,[],1)


###############################################################################
###############Graficos########################################################
#No se usa
f1 =#plot(1:T-1,base_mQua[1:end-1],label ="No testear", lw=3)
plot(1:T-1,ideal_mQua[1:end-1],label ="Simple, Antes de trabajar", lw=3)
#plot!(1:T-1,f1_mQua[1:end-1],label ="Simple, Antes de trabajar y cuarentena grupo", lw=3)
plot!(1:T-1,f1_mQua[1:end-1],label ="Pool, Antes de trabajar", lw=3)
plot!(1:T-1,fr_mQua[1:end-1],label ="No testear y cuarentena grupo", lw=3)
title!("Cuarentena Trabajando, con N=154")
xlabel!("Días")
ylabel!("# personas")



f12 =#plot(1:T-1,base_mQua2[1:end-1],label ="No testear", lw=3)
plot(1:T-1,ideal_mQua2[1:end-1],label ="Simple, Antes de trabajar", lw=3)
#plot!(1:T-1,ideal2_mQua2[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
#plot!(1:T-1,f1_mQua2[1:end-1],label ="Simple, Antes de trabajar y cuarentena grupo", lw=3)
plot!(1:T-1,f1_mQua2[1:end-1],label ="Pool, Antes de trabajar", lw=3)
plot!(1:T-1,fr_mQua2[1:end-1],label ="No testear y cuarentena grupo", lw=3)
title!("Personas en Cuarentena, con N=154")
xlabel!("Dias")
ylabel!("# personas")



f2 =#plot(1:T-1,base_mInf[1:end-1],label ="No testear", lw=3)
plot(1:T-1,ideal_mInf[1:end-1],label ="Simple, Antes de trabajar", lw=3)
#plot!(1:T-1,ideal2_mInf[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
#plot!(1:T-1,f1_mInf[1:end-1],label ="Simple, Antes de trabajar y cuarentena grupo", lw=3)
plot!(1:T-1,f1_mInf[1:end-1],label ="Pool, Antes de trabajar", lw=3)
plot!(1:T-1,fr_mInf[1:end-1],label ="No testear y cuarentena grupo", lw=3)
title!("Personas infectadas no descubiertas, con N=154")
xlabel!("Días")
ylabel!("# personas")

#f22 =plot(1:T-1,base_mInf2[1:end-1],label ="No testear", lw=3)
#plot!(1:T-1,ideal_mInf2[1:end-1],label ="Simple, Todos los dias", lw=3)
#plot!(1:T-1,ideal2_mInf2[1:end-1],label ="Simple, Todos los dias y cuarentena grupo", lw=3)
#plot!(1:T-1,f1_mInf2[1:end-1],label ="No testear y cuarentena grupo", lw=3)
#title!()
#xlabel!("Días")
#ylabel!("# personas")

f23 =#plot(1:T-1,base_mInf3[1:end-1],label ="No testear", lw=3)
plot(1:T-1,ideal_mInf3[1:end-1],label ="Simple, Antes de trabajar", lw=3)
#plot!(1:T-1,ideal2_mInf3[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
#plot!(1:T-1,f1_mInf3[1:end-1],label ="Simple, Antes de trabajar y cuarentena grupo", lw=3)
plot!(1:T-1,f1_mInf3[1:end-1],label ="Pool, Antes de trabajar", lw=3)
plot!(1:T-1,fr_mInf3[1:end-1],label ="No testear y cuarentena grupo", lw=3)
title!("Personas Infectadas, con N=154")
xlabel!("Días")
ylabel!("# personas")


f3 =#plot(1:T-1,base_T[1:end-1],label ="No testear", lw=3)
plot(1:T-1,ideal_T[1:end-1],label ="Simple, Antes de trabajar", lw=3)
#plot!(1:T-1,ideal2_T[1:end-1],label ="Simple, Cada dia los 4 grupos ver2", lw=3)
#plot!(1:T-1,f1_T[1:end-1],label ="Simple, Antes de trabajar y cuarentena grupo", lw=3)
plot!(1:T-1,f1_T[1:end-1],label ="Pool, Antes de trabajar", lw=3)
plot!(1:T-1,fr_T[1:end-1],label ="No testear y cuarentena grupo", lw=3)
title!("Test Realizados")
xlabel!("Días")
ylabel!("# Test")


f4 =#plot(1:T-1,base_mNFp[1:end-1],label ="base", lw=3)
plot(1:T-1,ideal_mNFp[1:end-1],label ="Simple, Antes de trabajar", lw=3)
#plot!(1:T-1,ideal2_mNFp[1:end-1],label ="Simple, Cada dia los 4 grupos ver2", lw=3)
#plot!(1:T-1,f1_mNFp[1:end-1],label ="Simple, Antes de trabajar y cuarentena grupo", lw=3)
plot!(1:T-1,f1_mNFp[1:end-1],label ="Pool, Antes de trabajar", lw=3)
plot!(1:T-1,fr_mNFp[1:end-1],label ="No testear y cuarentena grupo", lw=3)
title!("Falsos Positivos")
xlabel!("Días")
ylabel!("# Test")


f5 =#plot(1:T-1,base_mSy[1:end-1],label ="base", lw=3)
plot(1:T-1,ideal_mSy[1:end-1],label ="Simple, Antes de trabajar", lw=3)
#plot!(1:T-1,f1_mSy[1:end-1],label ="Simple, Antes de trabajar y cuarentena grupo", lw=3)
plot!(1:T-1,f1_mSy[1:end-1],label ="Pool, Antes de trabajar", lw=3)
plot!(1:T-1,fr_mSy[1:end-1],label ="No testear y cuarentena grupo", lw=3)
title!("Personas con Sintomas, con N=154")
xlabel!("Días")
ylabel!("# Personas")



f6 =#plot(1:T-1,base_mFQua[1:end-1],label ="No testear", lw=3)
plot(1:T-1,ideal_mFQua[1:end-1],label ="Simple, Antes de trabajar", lw=3)
#plot!(1:T-1,ideal2_mFQua[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
#plot!(1:T-1,f1_mFQua[1:end-1],label ="Simple, Antes de trabajar y cuarentena grupo", lw=3)
plot!(1:T-1,f1_mFQua[1:end-1],label ="Pool, Antes de trabajar", lw=3)
plot!(1:T-1,fr_mFQua[1:end-1],label ="No testear y cuarentena grupo", lw=3)
title!("Funcionarios en Cuarentena, con m=20")
xlabel!("Dias")
ylabel!("# personas")

f7 =#plot(1:T-1,base_mFSy[1:end-1],label ="base", lw=3)
plot(1:T-1,ideal_mFSy[1:end-1],label ="Simple, Antes de trabajar", lw=3)
#plot!(1:T-1,f1_mFSy[1:end-1],label ="Simple, Antes de trabajar y cuarentena grupo", lw=3)
plot!(1:T-1,f1_mFSy[1:end-1],label ="Pool, Antes de trabajar", lw=3)
plot!(1:T-1,fr_mFSy[1:end-1],label ="No testear y cuarentena grupo", lw=3)
title!("Funcionarios con Sintomas, con m=20")
xlabel!("Días")
ylabel!("# Personas")


f8 =#plot(1:T-1,base_mFInf[1:end-1],label ="No testear", lw=3)
plot(1:T-1,ideal_mFInf[1:end-1],label ="Simple, Antes de trabajar", lw=3)
#plot!(1:T-1,ideal2_mFInf[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
#plot!(1:T-1,f1_mFInf[1:end-1],label ="Simple, Antes de trabajar y cuarentena grupo", lw=3)
plot!(1:T-1,f1_mFInf[1:end-1],label ="Pool, Antes de trabajar", lw=3)
plot!(1:T-1,fr_mFInf[1:end-1],label ="No testear y cuarentena grupo", lw=3)
title!("Funcionarios infectados trabajando no descubiertos, con m=20")
xlabel!("Días")
ylabel!("# personas")


f9 =#plot(1:T-1,base_mFInf2[1:end-1],label ="No testear", lw=3)
plot(1:T-1,ideal_mFInf2[1:end-1],label ="Simple, Antes de trabajar", lw=3)
#plot!(1:T-1,ideal2_mFInf2[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
#plot!(1:T-1,f1_mFInf2[1:end-1],label ="Simple, Antes de trabajar y cuarentena grupo", lw=3)
plot!(1:T-1,f1_mFInf2[1:end-1],label ="Pool, Antes de trabajar", lw=3)
plot!(1:T-1,fr_mFInf2[1:end-1],label ="No testear y cuarentena grupo", lw=3)
title!("Funcionarios infectados, con m=20")
xlabel!("Días")
ylabel!("# personas")

ff1 =#plot(1:T-1,base_maxQua[1:end-1],label ="No testear", lw=3)
plot(1:T-1,ideal_maxQua[1:end-1],label ="Simple, Antes de trabajar", lw=3)
#plot!(1:T-1,ideal2_maxQua[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
#plot!(1:T-1,f1_maxQua[1:end-1],label ="Simple, Antes de trabajar y cuarentena grupo", lw=3)
plot!(1:T-1,f1_maxQua[1:end-1],label ="Pool, Antes de trabajar", lw=3)
plot!(1:T-1,fr_maxQua[1:end-1],label ="No testear y cuarentena grupo", lw=3)
title!("Cuarentena peor caso, con N=154 y m=0")
xlabel!("Días")
ylabel!("# personas")

ff2 =#plot(1:T-1,base_maxInf[1:end-1],label ="No testear", lw=3)
plot(1:T-1,ideal_maxInf[1:end-1],label ="Simple, Antes de trabajar", lw=3)
#plot!(1:T-1,ideal2_maxInf[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
#plot!(1:T-1,f1_maxInf[1:end-1],label ="Simple, Antes de trabajar y cuarentena grupo", lw=3)
plot!(1:T-1,f1_maxInf[1:end-1],label ="Pool, Antes de trabajar", lw=3)
plot!(1:T-1,fr_maxInf[1:end-1],label ="No testear y cuarentena grupo", lw=3)
title!("Infectados peor caso, con N=154 y con m=0")
xlabel!("Días")
ylabel!("# personas")

ff3 =#plot(1:T-1,base_maxSy[1:end-1],label ="No testear", lw=3)
plot(1:T-1,ideal_maxSy[1:end-1],label ="Simple, Antes de trabajar", lw=3)
#plot!(1:T-1,ideal2_maxSy[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
#plot!(1:T-1,f1_maxSy[1:end-1],label ="Simple, Antes de trabajar y cuarentena grupo", lw=3)
plot!(1:T-1,f1_maxSy[1:end-1],label ="Pool, Antes de trabajar", lw=3)
plot!(1:T-1,fr_maxSy[1:end-1],label ="No testear y cuarentena grupo", lw=3)
title!("Sintomas peor caso, con N=154 y con m=0")
xlabel!("Días")
ylabel!("# personas")




ff4 =#plot(base_VInf1,label ="No testear", lw=3)
plot(ideal_VInf1, label="Simple Antes de trabajar", lw=3)
#plot!(1:T-1,ideal2_maxInf[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
#plot!(f1_VInf1,label ="Simple, Antes de trabajar y cuarentena grupo", lw=3)
plot!(f1_VInf1,label ="Pool, Antes de trabajar", lw=3)
plot!(fr_VInf1,label ="No testear y cuarentena grupo", lw=3)
title!("Nuevos abuelitos infectados, con N=154")
xlabel!("Días")
ylabel!("# personas")

ff5 =#plot(base_VInf2,label ="No testear", lw=3)
plot(ideal_VInf2,label ="Simple, Antes de trabajar", lw=3)
#plot!(1:T-1,ideal2_maxSy[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
#plot!(f1_VInf2,label ="Simple, Antes de trabajar y cuarentena grupo", lw=3)
plot!(f1_VInf2,label ="Pool, Antes de trabajar", lw=3)
plot!(fr_VInf2,label ="No testear y cuarentena grupo", lw=3)
title!("Nuevos funcionarios infectados, con N=154")
xlabel!("Días")
ylabel!("# personas")


ff7 =#plot(acu(base_VInf1),label ="No testear", lw=3)
plot(acu(ideal_VInf1), label="Simple Antes de trabajar", lw=3)
#plot!(1:T-1,acu(ideal2_maxInf)[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
#plot!(acu(f1_VInf1),label ="Simple, Antes de trabajar y cuarentena grupo", lw=3)
plot!(acu(f1_VInf1),label ="Pool, Antes de trabajar", lw=3)
plot!(acu(fr_VInf1),label ="No testear y cuarentena grupo", lw=3)
title!("Acumulado nuevos abuelitos infectados, con N=154")
xlabel!("Días")
ylabel!("# personas")

ff8 =#plot(acu(base_VInf2),label ="No testear", lw=3)
plot(acu(ideal_VInf2),label ="Simple, Antes de trabajar", lw=3)
#plot!(1:T-1,ideal2_maxSy[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
#plot!(acu(f1_VInf2),label ="Simple, Antes de trabajar y cuarentena grupo", lw=3)
plot!(acu(f1_VInf2),label ="Pool, Antes de trabajar", lw=3)
plot!(acu(fr_VInf2),label ="No testear y cuarentena grupo", lw=3)
title!("Acumulado nuevos funcionarios infectados, con N=154")
xlabel!("Días")
ylabel!("# personas")

###################################################################
#############Se printean valores utiles############################
#print("Caso: Gente que se enfermo - Falsos Positivos - Test "," No testear: ",acu(base_VInf1)[T] ,"-", round(sum(base_mNFp)),"-",round(sum(base_T))," Simple, Cada día: ",acu(ideal_VInf1)[T] ,"-" ,round(sum(ideal_mNFp)),"-",round(sum(ideal_T))," Simple, Cada 10 días: ",acu(f1_VInf1)[T] ,"-",round(sum(f1_mNFp)),"-",round(sum(f1_T))," Simple, Cada 10 días repartido: ",acu(fr_VInf1)[T] ,"-",round(sum(fr_mNFp)),"-",round(sum(fr_T)))
print("Caso: Gente que se enfermo - Falsos Positivos - Test "," Simple, Cada día: ",acu(ideal_VInf1)[T] ,"-" ,round(sum(ideal_mNFp)),"-",round(sum(ideal_T)),"Poll, Antes de trabajar",acu(f1_VInf1)[T] ,"-",round(sum(f1_mNFp)),"-",round(sum(f1_T))," No testear, cuarentena: ",acu(fr_VInf1)[T] ,"-",round(sum(fr_mNFp)),"-",round(sum(fr_T)))

plot(f1,f1, layout = (2,1))

savefig("CuarentenaySintomas2.pdf")


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


plot(ff1,ff2,ff3, layout = (3,1))
savefig("PeorCaso.pdf")


plot(ff4,ff5, layout = (2,1))
savefig("Diarios.pdf")


plot(ff7,ff8, layout = (2,1))
savefig("Diarios2.pdf")

##################################################################################
#########Se realizan los histogramas de los maximos infectados por repeticion#####

infect=zeros(R,3)
#infect[:,1]=sort(base_Infect)
infect[:,1]=sort(ideal_Infect)
infect[:,2]=sort(f1_Infect)
infect[:,3]=sort(fr_Infect)
print(infect)


ffff=plot(density(infect'), label=NombresP)
title!("Maximos infectados por repetición, con N=100 y con m=20")
xlabel!("# personas")
ylabel!("Repeticiones")
savefig("Maximos.pdf")



finf1=histogram(infect[:,1],label =["No testear"])
title!("Maximos infectados por repetición, con N=100 y con m=20")
xlabel!("# personas")
ylabel!("Repeticiones")
savefig("Maximos1.pdf")

finf2=histogram(infect[:,2],label =["Simple, Todos los dias"],linealpha=0.3)
title!("Maximos infectados por repetición, con N=100 y con m=20")
xlabel!("# personas")
ylabel!("Repeticiones")
savefig("Maximos2.pdf")

finf3=histogram(infect[:,3],label =["Simple, Cada 10 dias"])
title!("Maximos infectados por repetición, con N=100 y con m=20")
xlabel!("# personas")
ylabel!("Repeticiones")
savefig("Maximos3.pdf")

#finf3=histogram(infect[:,4],label =["Simple, Cada 10 dias repartido"])
#title!("Maximos infectados por repetición, con N=100 y con m=20")
#xlabel!("# personas")
#ylabel!("Repeticiones")
#savefig("Maximos4.pdf")

end

#############################################################################
#################Datos caso asilo###########################################
if activador==1
	T = 210              # simulation horizon
	pB_int = 0   #prob of contagion en casa
	pB_ext = 0.01
	γ = 0.3 #(antes 0.25)           # prop of asymtomatic pop
	p_false_positive = 0.01 # probability of a false positive
	R = 40 # Replications
	peak_B=0.2 #Prob peak de infeccion
	S=1	#Cantidad de grupos
	N=100   #Cantidad de abuelitos
	T=100    #Cantidad de tiempo
	m=20      #Cantidad de funcionarios
	Trabajo=ones(2*N,T)    #Dias de trabajo abuelitos
	pA_int=0.01*(1/30)*0.2    #Probabilidad de contagio por visitas
	pA_ext=0
	peak_A=0.01                #probabilidad peak de contagio entre abulitos
	t_peak=24                 #Tiempo peak de contagio
	z=14			#Dias de cuarentena preventiva
	########################################################################################################################################################################
	#Frecuencias

	Frecuencia0=zeros(N+m,T) #Ningun dia
	Frecuencia1=ones(N+m,T) #Todos los dias
	Frecuencia12=ones(N+m,T) #Todos los dias
	for t in 1:T
		for i in 1:Int(m/2)
			if mod(t,2)==1
				Frecuencia12[N+Int(m/2)+1:N+m,t]=zeros(Int(m/2),1)
			else
				Frecuencia12[N+1 : N+Int(m/2),t]=zeros(Int(m/2),1)
			end
		end
	end
	Frecuencia10=zeros(N+m,T) #Cada 10 dias
	for t in 1:T
		if mod(t-1,10)==0
			Frecuencia10[:,t]=ones(N+m)
		end
	end
	Frecuencia5=zeros(N+m,T) #Cada 5 dias
	for t in 1:T
		if mod(t-1,5)==0
			Frecuencia5[:,t]=ones(N+m)
		end
	end
	Frecuencia102=zeros(N+m,T) #Cada 10 dias, con funcionarios todos los dias
	for t in 1:T
		if mod(t-1,10)==0
			Frecuencia102[:,t]=ones(N+m)
		end
	end
	Frecuencia102[N+1:N+m,:]=ones(m,T)




	Frecuenciar=zeros(N+m,T) #Cada 10 dias, repartido
	for t in 1:T
		for i in 1:10
			if mod(t-1,10)==i-1
				Frecuenciar[(i-1)*Int(2*N/m)+1:(i)*Int(2*N/m),t]=ones(Int(2*N/m))
				Frecuenciar[N+i,t]=1
				Frecuenciar[N+i+10,t]=1
			end
		end
	end
	Frecuenciar2=zeros(N+m,T) #Cada 5 dias repartidos, FALTA ARREGLAR
	for t in 1:T
		for i in 1:5
			if mod(t-1,5)==i-1
				Frecuenciar2[2*(i-1)*Int(2*N/m)+1:(i)*2*Int(2*N/m),t]=ones(2*Int(2*N/m))
				Frecuenciar2[N+2*i-1,t]=1
				Frecuenciar2[N+2*i,t]=1

			end

		end

	end

	Frecuenciar3=zeros(N+m,T) #Cada 10 dias, repartido, con funcionarios todos los dias
	for t in 1:T
		for i in 1:10
			if mod(t-1,10)==i-1
				Frecuenciar3[(i-1)*Int(2*N/m)+1:(i)*Int(2*N/m),t]=ones(Int(2*N/m))
			end
		end
	end
	Frecuenciar3[N+1:N+m,:]=ones(m,T)
	################################################
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

	Frecuenciar4=zeros(N+m,T) #Cada 10 dias, repartido, con funcionarios todos los dias
	for t in 1:T
		for i in 1:10
			if mod(t-1,10)==i-1
				Frecuenciar4[(i-1)*Int(N/10)+1:(i)*Int(N/10),t]=ones(Int(N/10))
				if mod(t,2)==0
					Frecuenciar4[N+1:N+Int(m/2),t]=ones(Int(m/2),1)
				else
					Frecuenciar4[N+Int(m/2)+1:N+m,t]=ones(Int(m/2),1)
				end
			end
		end
	end
	##############################################
	#Horarios e intereracciones funcionarios

	MFunc=zeros(2*(N+m),m,T) #Funcionarios no trabajan
	MFunc2=zeros(2*(N+m),m,T) #Funcionariso se relacionan con10 abuelitos dia por medio y no entre ellos
	for t in 1:T
		for i in 1:10
			MFunc2[(i-1)*20+1:i*20,i+(1-mod(t,2))*10,t]=ones(20)
		end
	end
	MFunc3=zeros(2*(N+m),m,T) #Funcionariso se relacionan con10 abuelitos dia por medio y entre los que trabajan
	for t in 1:T
		for i in 1:10
			MFunc3[(i-1)*20+1:i*20,i+(1-mod(t,2))*10,t]=ones(20)
		end
	end
	MFunc3[2*N+1:2*N+m,1:Int(m/2),:]=ones(m,Int(m/2),T)
	MFunc3[2*N+m+1:2*N+2*m,1:Int(m/2),:]=ones(m,Int(m/2),T)
	###############################################################################
	#visita=0
	#Se realizan las 4 politicas

politicas=4 #Numero de politicas
NombresP=["No testear","Simple, Todos los dias","Comparación","Simple, Cada 10 dias"]#Nombres politicas


base_mQua,base_mQua2, base_mInf,base_mInf2,base_mInf3, base_mNFp, base_T, baseg_t, base_mSy,base_mFQua,base_mFInf,base_mFInf2,base_mFSy,base_maxInf,base_maxQua, base_maxSy, base_Infect, base_VInf1, base_VInf2= simulation(N, Frecuencia0,0,T,peak_A,pA_int,pA_ext,peak_B,pB_int,pB_ext,γ,p_false_positive,R, false,Group1,Trabajo,"no","solo",z,t_peak,0,MFunc3,0)
ideal_mQua,ideal_mQua2, ideal_mInf,ideal_mInf2,ideal_mInf3, ideal_mNFp, ideal_T, idealg_t, ideal_mSy,ideal_mFQua,ideal_mFInf,ideal_mFInf2,ideal_mFSy,ideal_maxInf,ideal_maxQua, ideal_maxSy, ideal_Infect, ideal_VInf1, ideal_VInf2 = simulation(N, Frecuencia12,N,T,peak_A,pA_int,pA_ext,peak_B,pB_int,pB_ext,γ,p_false_positive,R, "t",Group1,Trabajo,"no","solo",z,t_peak,0,MFunc3,0)
#ideal2_mQua,ideal2_mQua2, ideal2_mInf,ideal2_mInf2,ideal2_mInf3, ideal2_mNFp, ideal2_T = simulation(N, fr44,N,T,β,γ,p_false_positive,R, true,p,Group1,Trab4,"no","solo",z,t_peak,MFunc2,visita,peak_A)
f1_mQua,f1_mQua2, f1_mInf,f1_mInf2,f1_mInf3, f1_mNFp, f1_T, f1g_t, f1_mSy,f1_mFQua,f1_mFInf,f1_mFInf2,f1_mFSy,f1_maxInf,f1_maxQua, f1_maxSy,f1_Infect, f1_VInf1, f1_VInf2 = simulation(N, Frecuenciar,N,T,peak_A,pA_int,pA_ext,peak_B,pB_int,pB_ext,γ,p_false_positive,R, "t",Group1,Trabajo,"no","solo",z,t_peak,0,MFunc3,0)
fr_mQua,fr_mQua2, fr_mInf,fr_mInf2,fr_mInf3, fr_mNFp, fr_T, frg_t, fr_mSy,fr_mFQua,fr_mFInf,fr_mFInf2,fr_mFSy,fr_maxInf,fr_maxQua, fr_maxSy, fr_Infect, fr_VInf1, fr_VInf2 = simulation(N, Frecuenciar4,N,T,peak_A,pA_int,pA_ext,peak_B,pB_int,pB_ext,γ,p_false_positive,R, "t",Group1,Trabajo,"no","solo",z,t_peak,0,MFunc3,0)



	###############################################################################
	###############Graficos########################################################
	#No se usa
#f1 =plot(1:T-1,P1_mQua[1:end-1],label =NombresP[1], lw=3)
#	if politicas>=1
#		for i in 2:politicas
#		plot!(1:T-1,P,i,_mQua[1:end-1],label =NombresP[i], lw=3)
#		end
#	end
#	title!("Cuarentena Trabajando")
#	xlabel!("Días")
#	ylabel!("# personas")



	f12 =plot(1:T-1,base_mQua2[1:end-1],label ="No testear", lw=3)
	plot!(1:T-1,ideal_mQua2[1:end-1],label ="Simple, Todos los dias", lw=3)
	#plot!(1:T-1,ideal2_mQua2[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
	plot!(1:T-1,f1_mQua2[1:end-1],label ="Comparación", lw=3)
	plot!(1:T-1,fr_mQua2[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
	title!("Personas en Cuarentena, con N=100")
	xlabel!("Dias")
	ylabel!("# personas")



	f2 =plot(1:T-1,base_mInf[1:end-1],label ="No testear", lw=3)
	plot!(1:T-1,ideal_mInf[1:end-1],label ="Simple, Todos los dias", lw=3)
	#plot!(1:T-1,ideal2_mInf[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
	plot!(1:T-1,f1_mInf[1:end-1],label ="Comparación", lw=3)
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
	plot!(1:T-1,f1_mInf3[1:end-1],label ="Comparación", lw=3)
	plot!(1:T-1,fr_mInf3[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
	title!("Personas Infectadas, con N=100")
	xlabel!("Días")
	ylabel!("# personas")


	f3 =plot(1:T-1,base_T[1:end-1],label ="No testear", lw=3)
	plot!(1:T-1,ideal_T[1:end-1],label ="Simple, Todos los dias", lw=3)
	#plot!(1:T-1,ideal2_T[1:end-1],label ="Simple, Cada dia los 4 grupos ver2", lw=3)
	plot!(1:T-1,f1_T[1:end-1],label ="Comparación", lw=3)
	plot!(1:T-1,fr_T[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
	title!("Test Realizados")
	xlabel!("Días")
	ylabel!("# Test")


	f4 =plot(1:T-1,base_mNFp[1:end-1],label ="base", lw=3)
	plot!(1:T-1,ideal_mNFp[1:end-1],label ="Simple, Todos los dias", lw=3)
	#plot!(1:T-1,ideal2_mNFp[1:end-1],label ="Simple, Cada dia los 4 grupos ver2", lw=3)
	plot!(1:T-1,f1_mNFp[1:end-1],label ="Comparación", lw=3)
	plot!(1:T-1,fr_mNFp[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
	title!("Falsos Positivos")
	xlabel!("Días")
	ylabel!("# Test")


	f5 =plot(1:T-1,base_mSy[1:end-1],label ="base", lw=3)
	plot!(1:T-1,ideal_mSy[1:end-1],label ="Simple, Todos los dias", lw=3)
	plot!(1:T-1,f1_mSy[1:end-1],label ="Comparación", lw=3)
	plot!(1:T-1,fr_mSy[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
	title!("Personas con Sintomas, con N=100")
	xlabel!("Días")
	ylabel!("# Personas")



	f6 =plot(1:T-1,base_mFQua[1:end-1],label ="No testear", lw=3)
	plot!(1:T-1,ideal_mFQua[1:end-1],label ="Simple, Todos los dias", lw=3)
	#plot!(1:T-1,ideal2_mFQua[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
	plot!(1:T-1,f1_mFQua[1:end-1],label ="Comparación", lw=3)
	plot!(1:T-1,fr_mFQua[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
	title!("Funcionarios en Cuarentena, con m=20")
	xlabel!("Dias")
	ylabel!("# personas")

	f7 =plot(1:T-1,base_mFSy[1:end-1],label ="base", lw=3)
	plot!(1:T-1,ideal_mFSy[1:end-1],label ="Simple, Todos los dias", lw=3)
	plot!(1:T-1,f1_mFSy[1:end-1],label ="Comparación", lw=3)
	plot!(1:T-1,fr_mFSy[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
	title!("Funcionarios con Sintomas, con m=20")
	xlabel!("Días")
	ylabel!("# Personas")


	f8 =plot(1:T-1,base_mFInf[1:end-1],label ="No testear", lw=3)
	plot!(1:T-1,ideal_mFInf[1:end-1],label ="Simple, Todos los dias", lw=3)
	#plot!(1:T-1,ideal2_mFInf[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
	plot!(1:T-1,f1_mFInf[1:end-1],label ="Comparación", lw=3)
	plot!(1:T-1,fr_mFInf[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
	title!("Funcionarios infectados trabajando no descubiertos, con m=20")
	xlabel!("Días")
	ylabel!("# personas")


	f9 =plot(1:T-1,base_mFInf2[1:end-1],label ="No testear", lw=3)
	plot!(1:T-1,ideal_mFInf2[1:end-1],label ="Simple, Todos los dias", lw=3)
	#plot!(1:T-1,ideal2_mFInf2[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
	plot!(1:T-1,f1_mFInf2[1:end-1],label ="Comparación", lw=3)
	plot!(1:T-1,fr_mFInf2[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
	title!("Funcionarios infectados, con m=20")
	xlabel!("Días")
	ylabel!("# personas")

	ff1 =plot(1:T-1,base_maxQua[1:end-1],label ="No testear", lw=3)
	plot!(1:T-1,ideal_maxQua[1:end-1],label ="Simple, Todos los dias", lw=3)
	#plot!(1:T-1,ideal2_maxQua[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
	plot!(1:T-1,f1_maxQua[1:end-1],label ="Comparación", lw=3)
	plot!(1:T-1,fr_maxQua[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
	title!("Cuarentena peor caso, con N=100 y m=20")
	xlabel!("Días")
	ylabel!("# personas")

	ff2 =plot(1:T-1,base_maxInf[1:end-1],label ="No testear", lw=3)
	plot!(1:T-1,ideal_maxInf[1:end-1],label ="Simple, Todos los dias", lw=3)
	#plot!(1:T-1,ideal2_maxInf[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
	plot!(1:T-1,f1_maxInf[1:end-1],label ="Comparación", lw=3)
	plot!(1:T-1,fr_maxInf[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
	title!("Infectados peor caso, con N=100 y con m=20")
	xlabel!("Días")
	ylabel!("# personas")

	ff3 =plot(1:T-1,base_maxSy[1:end-1],label ="No testear", lw=3)
	plot!(1:T-1,ideal_maxSy[1:end-1],label ="Simple, Todos los dias", lw=3)
	#plot!(1:T-1,ideal2_maxSy[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
	plot!(1:T-1,f1_maxSy[1:end-1],label ="Comparación", lw=3)
	plot!(1:T-1,fr_maxSy[1:end-1],label ="Simple, Cada 10 dias repartido", lw=3)
	title!("Sintomas peor caso, con N=100 y con m=20")
	xlabel!("Días")
	ylabel!("# personas")



	ff4 =plot(base_VInf1,label ="No testear", lw=3)
	plot!(ideal_VInf1, label="Simple Todos los dias", lw=3)
	#plot!(1:T-1,ideal2_maxInf[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
	plot!(f1_VInf1,label ="Comparación", lw=3)
	plot!(fr_VInf1,label ="Simple, Cada 10 dias repartido", lw=3)
	title!("Nuevos abuelitos infectados, con N=100 y con m=20")
	xlabel!("Días")
	ylabel!("# personas")

	ff5 =plot(base_VInf2,label ="No testear", lw=3)
	plot!(ideal_VInf2,label ="Simple, Todos los dias", lw=3)
	#plot!(1:T-1,ideal2_maxSy[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
	plot!(f1_VInf2,label ="Comparación", lw=3)
	plot!(fr_VInf2,label ="Simple, Cada 10 dias repartido", lw=3)
	title!("Nuevos funcionarios infectados, con N=100 y con m=20")
	xlabel!("Días")
	ylabel!("# personas")


	ff7 =plot(acu(base_VInf1),label ="No testear", lw=3)
	plot!(acu(ideal_VInf1), label="Simple Todos los dias", lw=3)
	#plot!(1:T-1,acu(ideal2_maxInf)[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
	plot!(acu(f1_VInf1),label ="Comparación", lw=3)
	plot!(acu(fr_VInf1),label ="Simple, Cada 10 dias repartido", lw=3)
	title!("Acumulado nuevos abuelitos infectados, con N=100 y con m=20")
	xlabel!("Días")
	ylabel!("# personas")

	ff8 =plot(acu(base_VInf2),label ="No testear", lw=3)
	plot!(acu(ideal_VInf2),label ="Simple, Todos los dias", lw=3)
	#plot!(1:T-1,ideal2_maxSy[1:end-1],label ="Simple, Cada día los 4 grupos ver2", lw=3)
	plot!(acu(f1_VInf2),label ="Comparación", lw=3)
	plot!(acu(fr_VInf2),label ="Simple, Cada 10 dias repartido", lw=3)
	title!("Acumulado nuevos funcionarios infectados, con N=100 y con m=20")
	xlabel!("Días")
	ylabel!("# personas")

	###################################################################
	#############Se printean valores utiles############################
	#print("Caso: Gente que se enfermo - Falsos Positivos - Test "," No testear: ",acu(base_VInf1)[T] ,"-", round(sum(base_mNFp)),"-",round(sum(base_T))," Simple, Cada día: ",acu(ideal_VInf1)[T] ,"-" ,round(sum(ideal_mNFp)),"-",round(sum(ideal_T))," Simple, Cada 10 días repartido: ",acu(fr_VInf1)[T] ,"-",round(sum(fr_mNFp)),"-",round(sum(fr_T)))
	print("Caso: Gente que se enfermo - Falsos Positivos - Test "," No testear: ",acu(base_VInf1)[T] ,"-", round(sum(base_mNFp)),"-",round(sum(base_T))," Simple, Cada día: ",acu(ideal_VInf1)[T] ,"-" ,round(sum(ideal_mNFp)),"-",round(sum(ideal_T))," Simple, Cada 10 días repartido: ",acu(fr_VInf1)[T] ,"-",round(sum(fr_mNFp)),"-",round(sum(fr_T)), "Comparación",acu(f1_VInf1)[T] ,"-",round(sum(f1_mNFp)),"-",round(sum(f1_T)))


	print("ideal")
	ideal_50=0
	ideal_75=0
	ideal_100=0
	for t in 2:T
		if (sum(ideal_VInf1[1:t-1])<25)&(sum(ideal_VInf1[1:t])>=25)
			ideal_50=sum(ideal_T[1:t])
			print("aa",ideal_50,"-")
		end

		if (sum(ideal_VInf1[1:t-1])<50)&(sum(ideal_VInf1[1:t])>=50)
			ideal_50=sum(ideal_T[1:t])
			print("aa",ideal_50,"-")
		end

		if (sum(ideal_VInf1[1:t-1])<75)&(sum(ideal_VInf1[1:t])>=75)
			ideal_75=sum(ideal_T[1:t])
			print("aa",ideal_75,"-")
		end

		if (sum(ideal_VInf1[1:t-1])<100)&(sum(ideal_VInf1[1:t])>=100)
			ideal_100=sum(ideal_T[1:t])
			print("aa",ideal_100,"-")
		end
	end

	print("f1")
	ideal_25=0
	ideal_50=0
	ideal_75=0
	ideal_100=0
	for t in 2:T
		if (sum(f1_VInf1[1:t-1])<25)&(sum(f1_VInf1[1:t])>=25)
			ideal_50=sum(f1_T[1:t])
			print("aa",ideal_50,"-")
		end
		if (sum(f1_VInf1[1:t-1])<50)&(sum(f1_VInf1[1:t])>=50)
			ideal_50=sum(f1_T[1:t])
			print("aa",ideal_50,"-")
		end

		if (sum(f1_VInf1[1:t-1])<75)&(sum(f1_VInf1[1:t])>=75)
			ideal_75=sum(f1_T[1:t])
			print("aa",ideal_75,"-")
		end

		if (sum(f1_VInf1[1:t-1])<100)&(sum(f1_VInf1[1:t])>=100)
			ideal_100=sum(f1_T[1:t])
			print("aa",ideal_100,"-")
		end
	end

	print("fr")
	ideal_50=0
	ideal_75=0
	ideal_100=0
	for t in 2:T
		if (sum(fr_VInf1[1:t-1])<25)&(sum(fr_VInf1[1:t])>=25)
			ideal_50=sum(fr_T[1:t])
			print("aa",ideal_50,"-")
		end
		if (sum(fr_VInf1[1:t-1])<50)&(sum(fr_VInf1[1:t])>=50)
			ideal_50=sum(fr_T[1:t])
			print("aa",ideal_50,"-")
		end

		if (sum(fr_VInf1[1:t-1])<75)&(sum(fr_VInf1[1:t])>=75)
			ideal_75=sum(fr_T[1:t])
			print("aa",ideal_75,"-")
		end

		if (sum(fr_VInf1[1:t-1])<100)&(sum(fr_VInf1[1:t])>=100)
			ideal_100=sum(fr_T[1:t])
			print("aa",ideal_100,"-")
		end
	end






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


	plot(ff1,ff2,ff3, layout = (3,1))
	savefig("PeorCaso.pdf")


	plot(ff4,ff5, layout = (2,1))
	savefig("Diarios.pdf")


	plot(ff7,ff8, layout = (2,1))
	savefig("Diarios2.pdf")


	##################################################################################
	#########Se realizan los histogramas de los maximos infectados por repeticion#####
	infect=zeros(R,4)
	infect[:,1]=sort(base_Infect)
	infect[:,2]=sort(ideal_Infect)
	infect[:,3]=sort(f1_Infect)
	infect[:,4]=sort(fr_Infect)
	print(infect)

	ffff=plot(density(infect), label=NombresP)
	title!("Maximos infectados por repetición, con N=100 y con m=20")
	xlabel!("# personas")
	ylabel!("Repeticiones")
	savefig("Maximos.pdf")

	finf1=histogram(infect[:,1],label =["No testear"])
	title!("Maximos infectados por repetición, con N=100 y con m=20")
	xlabel!("# personas")
	ylabel!("Repeticiones")
	savefig("Maximos1.pdf")

	finf2=histogram(infect[:,2],label =["Simple, Todos los dias"])
	title!("Maximos infectados por repetición, con N=100 y con m=20")
	xlabel!("# personas")
	ylabel!("Repeticiones")
	savefig("Maximos2.pdf")

	finf3=histogram(infect[:,3],label =["Simple, Cada 10 dias"])
	title!("Maximos infectados por repetición, con N=100 y con m=20")
	xlabel!("# personas")
	ylabel!("Repeticiones")
	savefig("Maximos3.pdf")

	finf3=histogram(infect[:,4],label =["Simple, Cada 10 dias repartido"])
	title!("Maximos infectados por repetición, con N=100 y con m=20")
	xlabel!("# personas")
	ylabel!("Repeticiones")
	savefig("Maximos4.pdf")




end
