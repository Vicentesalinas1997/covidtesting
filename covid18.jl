using Distributions, Plots, LinearAlgebra, Random, StatsPlots,JSON, Mustache, Printf

#######################################################################################
#Recibe un grupo y entrega un vector con los indices de los grupos a los que pertenecen
function vecgroup(group)
S=Int(length(group[1,:]))
N=Int(length(group[:,1]))
A=zeros(N,S)
v=zeros(N)
for j in 1:S
	for i in findall(group[:,j].==1)
		 group[i,j]==1
			v[i]=j
	end
end
return Int.(v)
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
	lar=length(vec[1,:])
	ler=length(vec[:,1])
	vec2=zeros(ler,lar)
	vec2[:,1]=vec[:,1]
	for i in 2:lar
		vec2[:,i]=vec[:,i]+vec2[:,i-1]
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
	lambda=(log(max(1-prob,0.000001)))/t_peak #Se calcula el lambda tal que P(tiempo para contagiarse<=t_peak)=probabilidad peak de contagio
	proba=1-exp(lambda*horas_turno) #Se calula P(tiempo para contagiarse<=horas_turno)
	return proba
end
##############################################################################




################################################################################
#####################Se usa para calcular las probabilidades de contagiar a cada persona de un grupo, entrega un vector con las probabilidades de
#####################cada integrante del grupo y 0 en los que no son del grupo
#Funcion que entrega la probabilidad de contagiarte dentro de tu grupo de trabajo cerrado
function group_p(M_cont)
	#G: indicatriz del grupo de trabajo
	#vp_cont: vector con las probabilidades de contagio
	#W: Peligro de contagio de paciente
	L=length(M_cont[1,:])
	p=zeros(L)
	p=(-prod(-M_cont.+1,dims=1)).+1
	#de no contagiarse multiplicado por un factor que pondera las horas, este desde
	return p'
end

#########################################################################################

##########################################################################################
##########################Función que crea vectores de frecuencia#########################
#Recibe, offset, cada cuanto se realiza la frecuencia y un horizonte de tiempo
function frecuenciarand(offset,f,T)
	F=zeros(T)
	for i in 1:T
		if i>=offset
			if mod(i-offset,f)==0
				F[i]=1
			end
		end
	end
return F
end
############################################################################################
###########################################################################################

##############################################################################################
###########################Funcion repetir###################################################
#Funcion que repite o corta un array hasta tener un cierto largo (al lado) y Ancho (abajo)####################################
function repetir(array,A,L)
ARRAY=zeros(A,L)
l=length(array[1,:])
a=length(array[:,1])
	for ind in 1:A
		ARRAY[ind,1:min(L,l)]=array[a-mod(-ind,a),1:min(L,l)]
	end
	if L>l
		for indi in l+1:L
			ARRAY[:,indi]=ARRAY[:,l-mod(-indi,l)]
		end
	end
	return ARRAY
end
##########################################################################################
#########################################################################################


############################################################################################################
#Funcion que recibe los grupos y la matriz de relacion entre grupos y entrega las relaciones para todos los dias
function interaciones(Group,Mrel,T,v)
largo=length(Group[:,1])
Mint=zeros(largo,largo,T)
for t in 1:T
	for i in 1:largo
		for j in i+1:largo
			Mint[i,j,t]=(rand()<=Mrel[v[i],v[j]])
			Mint[j,i,t]=Mint[i,j,t]
		end
	end
end
return Mint
end

function pool(array,num)
	if length(array)==0
		return []
	else
		pools=[]
		for i in 1:length(array[1,:])
			cand=array[:,i].*collect(1:length(array[:,i]))
			cand2=setdiff(cand,[0])
			cand3=shuffle!(cand2)
			if length(cand3)>=1
				l=Int(floor(length(cand3)/num))
				r=Int(length(cand3)%num)
				if l==0
					if r>=1
						if length(pools)>=1
							pools=[pools' [Int.(cand3[1:r])]]'
						else
							pools=[Int.(cand3[1:r])]'
						end
					end
				else
					pool=[Int.(cand3[(i-1)*num+1:i*num]) for i=1:l]
					if length(pools)>=1
						pools=[pools' pool']'
					else
						pools=pool
					end
					if r>=1
						pools=[pools' [Int.(cand3[1:r])]]'
					end
				end
			end
		end
	end

	return pools
end
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
#VInf: Numero de infectados en promedio, por dia
function simulation(N,Group,T,Tra,Mrel, f,G,peak,t_peak,p_int,p_ext,γ,p_false_positive,R,Politica,tpool,random,quienes,cuarentena,z,dias_atras,scalar_asint,test_sym,distribuir)
#N Tamaño grupo de trabajadores
#Group una matriz con columnas indicatrices de tamaño N con los integrantes de cada grupo.
#Tra; horario de los trabajadores
#Mrel: Matriz de presencia de todos los trabajadores, tiene una tercer componente temporal
#f Frecuencia testeo (puede ser vector de largo T, en caso de testear random de tamaño G o matriz de NxT con los individuos a testear en sus fechas)
# G Solo aplica para el caso random porcentaje por grupo a testear
# T simulation horizon
#peak: vector de probabilidad contagio peak
#γ:  prop of asymtomatic pop
#p_false_positive probability of a false positive
#R Replications
#Politica dice si es "Pool", "Individual" o "No testear".
#tpool: Tamanaño Pool
#Tra: Matriz (2*N,T) de turnos de trabajo
#random: Si es random o estan fijos a quienes se testeara (no los dias esos estan en f, solo los funcionarios a testear)
#quienes: trabajando, no trabajando y ambos
#cuarentena: Si se mandan a todo el grupo en cuarentena o solo al infectado (ahora tambien mixto, solo al grupo cuando uno presenta sintomas)
#z: dias de cuarntena al grupo, no se usa en caso de solo mandar al positivo
#t_peak: tiempo en que la exponencial alcanza la probabilidad peak
#dias_atras: cuantos dias hacia atras se mira con quienes has tenido contacto estrecho
#test_sym: Si se quiere testear a la gente apenas le aparescan los sintomas
#distribuir: Si se distribuye o no los test (solo aplica para random)
hr=12
MD=zeros((N),T) #Matriz trabajo dia
MN=zeros((N),T) ##Matriz trabajo noche
for i in 1:N
	MD[i,:]=Tra[2*i-1,:]
	MN[i,:]=Tra[2*i,:]
end

S=length(Group[1,:]) #Tamaño grupo
#Cuando se utiliza la cuarentena, con este vector se recuperan los integrantes del grupo para mandarlos a cuarentena
v=vecgroup(Group)
#Vectores que guardan los valores de cuarentena, infectados, sintomas y números de test en cada t
    mQua = zeros(N,T)
    mQua2 = zeros(N,T)
    mFQua = zeros(N,T)
    mInf = zeros(N,T)
    mInf2 = zeros(N,T)
    mInf3 = zeros(N,T)
    mFInf=zeros(N,T)
    mFInf2=zeros(N,T)
    mSy= zeros(N,T)
    mFSy=zeros(N,T)
	mNTest=zeros(T)
	mNTest2=zeros(T)
	mNFpositive=zeros(T)
maxInf=zeros(N,T)
maxQua=zeros(N,T)
maxSy=zeros(N,T)
Infect=ones(R)*1000
VInf=zeros(N,T)



    # Loop over replications
    for rep = 1:R
		print(rep)
	Mint=interaciones(Group,Mrel,T,v) # Matriz de interacciones para cada repetición
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
	NTest=zeros(T) #Vector de test diarios
	NTest2=zeros(T) #Vector de test diarios por la modalidad sintomaticos
	NFpositive=zeros(T) #Falsos positivos diarios
	indicatriz=zeros(N,T) #Quienes testear de cada grupo caso random
	resagados=zeros(N)
	resagados2=zeros(N)
	t_inicio=1
	vector=[]
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
			new_inf=zeros(N) #Vector de los nuevos infectados
			MC=zeros(N,N) #Matriz de contagio

#reiniciar en f[t] el resagados
			if (random=="Si")
				if f[t]==1
					resagados=zeros(N)
					t_inicio=t
				end
			end
##########################################

		for a=1:S #Se recorren los grupos
			g=Group[:,a] #Se recupera la info codificada del grupo
########################################################################################
			for i in findall((g.*inf.*(-qu.+1).*(-re.+1)).==1)#Se recorren los abuelitos enfermos, no en cuarentena y no recuperados
				casos=shuffle(Mint[:,i,t].*su.*(-qu.+1).*collect(1:N))
				casos2=Int.(setdiff(casos,[0]))
				if length(casos2)>0
				for j in  casos2 #Se recorren los abuelitos que nunca se han contagiado y no en cuarentena
						#Recuperacion de los tiempos
						tinf=tIn[Int(i)]-tEx[Int(i)]
						tsin=tSy[i]-tEx[i]
						trec=tRe[i]-tEx[i]
						pond=0 #Ponderador de horas
						if (MD[i,t]*MD[j,t]==1) & ((1-MN[i,t])*(1-MN[j,t])==1)  #Condicion esten ambos trabajando
							pond=1
						end
						if ((1-MD[i,t])*(1-MD[j,t])==1) & (MN[i,t]*MN[j,t]==1)
							pond=1
						end
						if (MD[i,t]*MD[j,t]==1) & (MN[i,t]*MN[j,t]==1)
							pond=2
						end
						MC[i,j]=exponential_p(t,tinf,tsin,trec,peak[i],t_peak[i],pond*hr,As[i],scalar_asint)
				end
				end
			end
#################################################################################################################
			#Para el pool test, predefino el numero de grupos que se pueden formar
			#se ve quienes son el publico a testear
			if Politica=="Pool"
				if (random=="Si") & (f[t]==1)
					for i in findall(g.==1)
						indicatriz[i,t]=(rand()<=G[a])
					end
					if quienes=="Trabajando"
					I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).* (MD[:,t]+MN[:,t].>=1).*(-qu.+1) #Gente que esta trabajando, no en cuarentena y nunca a presentado sintomas
					end
					if quienes=="No Trabajando"
					I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).* (MD[:,t]+MN[:,t].<=1).*(-qu.+1) #Gente que esta trabajando, no en cuarentena y nunca a presentado sintomas
					end
					if quienes=="Ambos"
					I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*(-qu.+1) #Gente que esta trabajando, no en cuarentena y nunca a presentado sintomas
					end
					resagados+=g.*I.*indicatriz[:,t]
				end
			end
####################################################################################################################
	   	end
#################################################################################################################################
PG=group_p(MC) #Vector con probabilidades de contagio dentro del grupo
r=rand(N)
mm=(MD[:,t]+MN[:,t])/2
new_inf+=(-qu.+1).*su.*(r.<(-(-PG.+1).*(-mm.*p_int.+1).*(-(-mm.+1).*p_ext.+1).+1)) #Vector de infectados del grupo que no esten en cuarentena, esten o no trabajando
#new_inf+=(-qu.+1).*su.*(r.<(-(-PG.+1).*(-mm.*p_int.+1).*(-p_ext.+1).+1)) #Vector de infectados del grupo que no esten en cuarentena, esten o no trabajando
new_inf+=qu.*su.*(r.<(p_ext)) #Sumar infectados en cuarentena (falsos positivos)

##########################################################################################
#Se realiza el proceso de cambio en los estados de los nuevos infectados
           	for j in findall(new_inf.==1)                   # cicle over newly infected
				i=j[1]
					VInf[i,t]+=1
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

		if (t>=3)&((cuarentena=="Mixto")|(cuarentena=="Grupal"))              #Si se activa se manda a todo el grupo cuando alguien tiende sintomas
			for i in findall(Sy[:,t].*(Qu[:,t]).*Sy[:,t-1].*(-Qu[:,t-1].+1).*Sy[:,t-2].*(-Qu[:,t-2].+1).==1) #Condicion de no ser pillado hasta los sintomas
				Aux=MD[Int(i),Int(max(t-dias_atras,1)):t]'.*Mint[Int(i),:,Int(max(t-dias_atras,1)):t]
					for te in Int(max(t-dias_atras,1)):t
						for j in findall(Aux[:,te-Int(max(t-dias_atras,1))+1].*(-Qu[:,te].+1).*(-qu.+1).*MD[:,te].==1)
							Qu[j,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
						end
					end
				Aux2=MN[i,Int(max(t-dias_atras,1)):t]'.*Mint[i,:,Int(max(t-dias_atras,1)):t]
						for te2 in Int(max(t-dias_atras,1)):t
							for j2 in findall(Aux2[:,te2-Int(max(t-dias_atras,1))+1].*(-Qu[:,te2].+1).*(-qu.+1).*MN[:,te2].==1)
								Qu[j2,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
							end
						end
			end
		end
##########################Testing###################################################
		#Vector de testear
	    	cand=[]
###################################Caso random, se a usado, poco es cuando solo se sabe que dias y a cuantos se quiere testear,
###################################con estos numeros se saca una muestra de minimo entre tamaño G y los que trabajan al azar. De momento esta descontinuado
		if test_sym=="Si"
			for i in findall(sy.*(-qu.+1).==1)
									NTest2[t]+=1   #Un test nuevo
									p_positive = true_positive(t-tEx[i]+1, tIn[i]-tEx[i]+1, tSy[i]-tEx[i]+1, tRe[i]-tEx[i]+1,As[i],scalar_asint)
									if rand() < p_positive #Caso positivo
											Qu[i,Int(min(T,t+1)):Int(min(T,max(t+15,tRe[i]+3)))] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
											if cuarentena=="Grupal"
												Aux=MD[Int(i),Int(max(t-dias_atras,1)):t]'.*Mint[Int(i),:,Int(max(t-dias_atras,1)):t]
													for te in Int(max(t-dias_atras,1)):t
														for j in findall(Aux[:,te-Int(max(t-dias_atras,1))+1].*(-Qu[:,te].+1).*(-qu.+1).*MD[:,te].==1)
															Qu[j,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
														end
													end
												Aux2=MN[i,Int(max(t-dias_atras,1)):t]'.*Mint[i,:,Int(max(t-dias_atras,1)):t]
														for te2 in Int(max(t-dias_atras,1)):t
															for j2 in findall(Aux2[:,te2-Int(max(t-dias_atras,1))+1].*(-Qu[:,te2].+1).*(-qu.+1).*MN[:,te2].==1)
																Qu[j2,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
															end
														end
											end
									end
			end
		end
		if random=="Si" #Se realiza test al azar
	    		if Politica=="Individual" #test individual
					if f[t]==1 #Se testea estos dias
						if quienes=="Trabajando"
							I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).* (MD[:,t]+MN[:,t].>=1).*(-qu.+1).*collect(1:N) #Gente que esta trabajando, no en cuarentena y nunca a presentado sintomas
						end
						if quienes=="No Trabajando"
							I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).* (MD[:,t]+MN[:,t].<=1).*(-qu.+1).*collect(1:N) #Gente que esta trabajando, no en cuarentena y nunca a presentado sintomas
						end
						if quienes=="Ambos"
							I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*(-qu.+1).*collect(1:N) #Gente que esta trabajando, no en cuarentena y nunca a presentado sintomas
						end
						J=shuffle!(I) #Mezclar
						II=setdiff(J,[0]) #Eliminar ceros
						cand=Int.(II)  #Se toman como candidatos a los posibles y el minimo entre ellos y G
						if distribuir=="Si"
							resagados=cand
							t_inicio=t
						end
					end
					if distribuir=="Si"
						cand=[]
						if t==1
							if quienes=="Trabajando"
								I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).* (MD[:,t]+MN[:,t].>=1).*(-qu.+1).*collect(1:N) #Gente que esta trabajando, no en cuarentena y nunca a presentado sintomas
							end
							if quienes=="No Trabajando"
								I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).* (MD[:,t]+MN[:,t].<=1).*(-qu.+1).*collect(1:N) #Gente que esta trabajando, no en cuarentena y nunca a presentado sintomas
							end
							if quienes=="Ambos"
								I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*(-qu.+1).*collect(1:N) #Gente que esta trabajando, no en cuarentena y nunca a presentado sintomas
							end
							J=shuffle!(I) #Mezclar
							II=setdiff(J,[0]) #Eliminar ceros
							cand=Int.(II)  #Se toman como candidatos a los posibles y el minimo entre ellos y G
							resagados=cand
							t_inicio=t
						end
						if length(resagados)<t-t_inicio
							cand=[]
						else
							resagados2=zeros(N)
							d=findall(f.==1)
							frec=(d[2]-d[1])[2]
							for i in 1:length(resagados)
								if t-t_inicio==mod(i-1,frec)
									resagados2[i]=resagados[i]
								end
							end
							cand=Int.(setdiff(resagados2,[0]))
						end
					end
	    		end
				if (Politica=="Pool")
					if f[t]==1
						vector=pool(resagados,tpool)
					else
						vector=[]
					end
					if distribuir=="Si"
						resagados2=zeros(N)
						d=findall(f.==1)
						frec=(d[2]-d[1])
						for i in 1:N
							if t-t_inicio==mod(i-1,frec)
								resagados2[i]=resagados[i]
							end
						end
						vector=pool(resagados2,tpool)
					end
					if length(vector)>=1
						for cand in vector
							#hacer pool testing
							pp=zeros(N)
							for i in cand
								p_positive = true_positive(t-tEx[i]+1, tIn[i]-tEx[i]+1, tSy[i]-tEx[i]+1, tRe[i]-tEx[i]+1,As[i],scalar_asint)
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
									if cuarentena=="Grupal"
										Aux=MD[Int(i),Int(max(t-dias_atras,1)):t]'.*Mint[Int(i),:,Int(max(t-dias_atras,1)):t]
										for te in Int(max(t-dias_atras,1)):t
											for j in findall(Aux[:,te-Int(max(t-dias_atras,1))+1].*(-Qu[:,te].+1).*(-qu.+1).*MD[:,te].==1)
												Qu[j,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
											end
										end
										Aux2=MN[i,Int(max(t-dias_atras,1)):t]'.*Mint[i,:,Int(max(t-dias_atras,1)):t]
										for te2 in Int(max(t-dias_atras,1)):t
											for j2 in findall(Aux2[:,te2-Int(max(t-dias_atras,1))+1].*(-Qu[:,te2].+1).*(-qu.+1).*MN[:,te2].==1)
												Qu[j2,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
											end
										end
									end
									NFpositive[t]+=1
								end
								if (su[i]==0)&(re[i]==0)&(i in cand)
	           						Qu[i,Int(min(T,t+1)):Int(min(T,max(t+15,tRe[i]+3)))] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
									if cuarentena=="Grupal"
										Aux=MD[Int(i),Int(max(t-dias_atras,1)):t]'.*Mint[Int(i),:,Int(max(t-dias_atras,1)):t]
										for te in Int(max(t-dias_atras,1)):t
											for j in findall(Aux[:,te-Int(max(t-dias_atras,1))+1].*(-Qu[:,te].+1).*(-qu.+1).*MD[:,te].==1)
												Qu[j,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
											end
										end
										Aux2=MN[i,Int(max(t-dias_atras,1)):t]'.*Mint[i,:,Int(max(t-dias_atras,1)):t]
										for te2 in Int(max(t-dias_atras,1)):t
											for j2 in findall(Aux2[:,te2-Int(max(t-dias_atras,1))+1].*(-Qu[:,te2].+1).*(-qu.+1).*MN[:,te2].==1)
												Qu[j2,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
											end
										end
									end
	                       		end
							end
						end
					end
				end
					if (Politica=="Individual")	& (length(cand)>=1)	#Si al menos hay 1 para testear y no pool
						for i in cand	#se recorren
							if rand()<=G[v[i]]
							NTest[t]+=1   #Un test nuevo
	           				     		if (((As[i]==1)&(su[i]==0)&(re[i]==1))|(su[i]==1))& (rand() < p_false_positive) #Caso falso positivo
				                    		Qu[i,min(T,t+1):min(T,t+15)] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
											if cuarentena=="Grupal"
												Aux=MD[Int(i),Int(max(t-dias_atras,1)):t]'.*Mint[Int(i),:,Int(max(t-dias_atras,1)):t]
													for te in Int(max(t-dias_atras,1)):t
														for j in findall(Aux[:,te-Int(max(t-dias_atras,1))+1].*(-Qu[:,te].+1).*(-qu.+1).*MD[:,te].==1)
															Qu[j,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
														end
													end
												Aux2=MN[i,Int(max(t-dias_atras,1)):t]'.*Mint[i,:,Int(max(t-dias_atras,1)):t]
														for te2 in Int(max(t-dias_atras,1)):t
															for j2 in findall(Aux2[:,te2-Int(max(t-dias_atras,1))+1].*(-Qu[:,te2].+1).*(-qu.+1).*MN[:,te2].==1)
																Qu[j2,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
															end
														end
											end
											NFpositive[t]+=1
	                					end
	                				if (su[i]==0)&(re[i]==0)
	                        					p_positive = true_positive(t-tEx[i]+1, tIn[i]-tEx[i]+1, tSy[i]-tEx[i]+1, tRe[i]-tEx[i]+1,As[i],scalar_asint)
	            			            		if rand() < p_positive #Caso positivo
	                            						Qu[i,Int(min(T,t+1)):Int(min(T,max(t+15,tRe[i]+3)))] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
														if cuarentena=="Grupal"
															Aux=MD[Int(i),Int(max(t-dias_atras,1)):t]'.*Mint[Int(i),:,Int(max(t-dias_atras,1)):t]
																for te in Int(max(t-dias_atras,1)):t
																	for j in findall(Aux[:,te-Int(max(t-dias_atras,1))+1].*(-Qu[:,te].+1).*(-qu.+1).*MD[:,te].==1)
																		Qu[j,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
																	end
																end
															Aux2=MN[i,Int(max(t-dias_atras,1)):t]'.*Mint[i,:,Int(max(t-dias_atras,1)):t]
																	for te2 in Int(max(t-dias_atras,1)):t
																		for j2 in findall(Aux2[:,te2-Int(max(t-dias_atras,1))+1].*(-Qu[:,te2].+1).*(-qu.+1).*MN[:,te2].==1)
																			Qu[j2,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
																		end
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
		if random=="No" #Se realiza el esquema de test de manera similar solo que los candidatos ya estan dados por defecto y solo saca los que esten en cuarentena o ya hayan tenido sintomas en el pasado
			if sum( f[:,t].*(-qu.+1))>=1
	  	  		if Politica=="Individual"
					J=Int.(( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*f[:,t].*(-qu.+1).*collect(1:N)) #Son los que nunca han tenido sintomas, no cuarentena y toca testear
					II=setdiff(J,[0])
					cand=Int.(II)  #Se obtiene vector de a quienes testear
	    		end
				if Politica=="Pool" #Esta medio obsoleto
					h=(sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*f[:,t].*(-qu.+1)
					vector=pool(h,tpool)
					if length(vector)>=1
						for cand in vector
							pp=zeros(N)
							for i in cand
								p_positive = true_positive(t-tEx[i]+1, tIn[i]-tEx[i]+1, tSy[i]-tEx[i]+1, tRe[i]-tEx[i]+1,As[i],scalar_asint)
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
									if cuarentena=="Grupal"
										Aux=MD[Int(i),Int(max(t-dias_atras,1)):t]'.*Mint[Int(i),:,Int(max(t-dias_atras,1)):t]
										for te in Int(max(t-dias_atras,1)):t
											for j in findall(Aux[:,te-Int(max(t-dias_atras,1))+1].*(-Qu[:,te].+1).*(-qu.+1).*MD[:,te].==1)
												Qu[j,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
											end
										end
										Aux2=MN[i,Int(max(t-dias_atras,1)):t]'.*Mint[i,:,Int(max(t-dias_atras,1)):t]
											for te2 in Int(max(t-dias_atras,1)):t
												for j2 in findall(Aux2[:,te2-Int(max(t-dias_atras,1))+1].*(-Qu[:,te2].+1).*(-qu.+1).*MN[:,te2].==1)
													Qu[j2,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
												end
											end
									end
									NFpositive[t]+=1
								end
								if (su[i]==0)&(re[i]==0)&(i in cand)
           							Qu[i,Int(min(T,t+1)):Int(min(T,max(t+15,tRe[i]+3)))] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
									if cuarentena=="Grupal"
										Aux=MD[Int(i),Int(max(t-dias_atras,1)):t]'.*Mint[Int(i),:,Int(max(t-dias_atras,1)):t]
											for te in Int(max(t-dias_atras,1)):t
												for j in findall(Aux[:,te-Int(max(t-dias_atras,1))+1].*(-Qu[:,te].+1).*(-qu.+1).*MD[:,te].==1)
													Qu[j,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
												end
											end
											Aux2=MN[i,Int(max(t-dias_atras,1)):t]'.*Mint[i,:,Int(max(t-dias_atras,1)):t]
											for te2 in Int(max(t-dias_atras,1)):t
												for j2 in findall(Aux2[:,te2-Int(max(t-dias_atras,1))+1].*(-Qu[:,te2].+1).*(-qu.+1).*MN[:,te2].==1)
													Qu[j2,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
												end
											end
									end
                        		end
							end
						end
					end
				end
				if (Politica=="Individual")	& (length(cand)>=1)	#Si al menos hay 1 para testear y no pool
					for i in cand	#se recorren
						NTest[t]+=1   #Un test nuevo
           				     		if (((As[i]==1)&(su[i]==0)&(re[i]==1))|(su[i]==1))& (rand() < p_false_positive) #Caso falso positivo
			                    		Qu[i,min(T,t+1):min(T,t+15)] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
										if cuarentena=="Grupal"
											Aux=MD[Int(i),Int(max(t-dias_atras,1)):t]'.*Mint[Int(i),:,Int(max(t-dias_atras,1)):t]
												for te in Int(max(t-dias_atras,1)):t
													for j in findall(Aux[:,te-Int(max(t-dias_atras,1))+1].*(-Qu[:,te].+1).*(-qu.+1).*MD[:,te].==1)
														Qu[j,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
													end
												end
											Aux2=MN[i,Int(max(t-dias_atras,1)):t]'.*Mint[i,:,Int(max(t-dias_atras,1)):t]
													for te2 in Int(max(t-dias_atras,1)):t
														for j2 in findall(Aux2[:,te2-Int(max(t-dias_atras,1))+1].*(-Qu[:,te2].+1).*(-qu.+1).*MN[:,te2].==1)
															Qu[j2,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
														end
													end
										end
										NFpositive[t]+=1
                					end
                				if (su[i]==0)&(re[i]==0)

                        					p_positive = true_positive(t-tEx[i]+1, tIn[i]-tEx[i]+1, tSy[i]-tEx[i]+1, tRe[i]-tEx[i]+1,As[i],scalar_asint)
            			            		if rand() < p_positive #Caso positivo
                            						Qu[i,Int(min(T,t+1)):Int(min(T,max(t+15,tRe[i]+3)))] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
													if cuarentena=="Grupal"
														Aux=MD[Int(i),Int(max(t-dias_atras,1)):t]'.*Mint[Int(i),:,Int(max(t-dias_atras,1)):t]
															for te in Int(max(t-dias_atras,1)):t
																for j in findall(Aux[:,te-Int(max(t-dias_atras,1))+1].*(-Qu[:,te].+1).*(-qu.+1).*MD[:,te].==1)
																	Qu[j,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
																end
															end
														Aux2=MN[i,Int(max(t-dias_atras,1)):t]'.*Mint[i,:,Int(max(t-dias_atras,1)):t]
																for te2 in Int(max(t-dias_atras,1)):t
																	for j2 in findall(Aux2[:,te2-Int(max(t-dias_atras,1))+1].*(-Qu[:,te2].+1).*(-qu.+1).*MN[:,te2].==1)
																		Qu[j2,Int(min(T,t)):Int(min(T,t+z))] .= 1 #quarantined for 2 weeks
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
	#Se guardan los valores en cada
    mQua .+=((MD.*Qu+MN.*Qu).>=1)
	mQua2 .+= Qu
	mInf .+= ((MD.*In.*(-Qu.+1)+MN.*In.*(-Qu.+1)).>=1)
    mInf2 .+= ((MD.*In+MN.*In).>=1)
	mInf3 .+= In
	mSy .+= Sy
	mNFpositive .+= NFpositive
	mNTest .+= NTest
	mNTest2 .+= NTest2
	maxi=max(   maximum(sum(maxInf,dims=1))  ,  maximum(sum(In[1:N,:],dims=1))   )
	si=maximum(sum(maxInf,dims=1))==(maxi)
	maxInf=In.*(1-si)+maxInf*si
 	if maxInf==In
		maxQua=Qu
		maxSy=Sy
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
	mNTest2=mNTest2/R
	VInf=VInf/R


    return mQua,mQua2,mInf,mInf2,mInf3, mNFpositive,  mNTest,mNTest2, mSy, maxInf, maxQua, maxSy, Infect, VInf
end
################################################################################################################################################
function leer(dict12)
	NameP=dict12["NPolitica"]
	T=dict12["Tiempo"]
	S=length(dict12["Grupos"])
	Integrantes=0
	for s in 1:S
	for i in 1:length(dict12["Grupos"][s])
		Integrantes+=dict12["Grupos"][s][i]["Numero de copias"]
	end
	end
	N=Integrantes
	Group=zeros(Integrantes,S)
	p_int=zeros(Integrantes)
	p_ext=zeros(Integrantes)
	NombreGrupo=["a" for j in 1:Integrantes]
	NombreHorario=["a" for j in 1:Integrantes]
	Horario=zeros(2*Integrantes,T)
	Mrel=zeros(S,S)
	cont=0
	cont2=1
	infopar=zeros(Integrantes)
	for s in 1:S
	for i in 1:length(dict12["Grupos"][s])
		Group[cont+1:dict12["Grupos"][s][i]["Numero de copias"]+cont,s]=ones(dict12["Grupos"][s][i]["Numero de copias"])
		p_int[cont+1:dict12["Grupos"][s][i]["Numero de copias"]+cont]=ones(dict12["Grupos"][s][i]["Numero de copias"]).*dict12["Grupos"][s][i]["Probabilidad de contagio int"]
		p_ext[cont+1:dict12["Grupos"][s][i]["Numero de copias"]+cont]=ones(dict12["Grupos"][s][i]["Numero de copias"]).*dict12["Grupos"][s][i]["Probabilidad de contagio ext"]
		NombreGrupo[cont+1:dict12["Grupos"][s][i]["Numero de copias"]+cont]=[dict12["Grupos"][s][i]["Grupo"] for j in 1:dict12["Grupos"][s][i]["Numero de copias"]]
		NombreHorario[cont+1:dict12["Grupos"][s][i]["Numero de copias"]+cont]=[dict12["Grupos"][s][i]["Nhorario"] for j in 1:dict12["Grupos"][s][i]["Numero de copias"]]
		infopar[cont2]=dict12["Grupos"][s][i]["Numero de copias"]
		cont2+=1
	for t in 1:T
		Horario[2*cont+1:2*(dict12["Grupos"][s][i]["Numero de copias"]+cont),t]=repeat(dict12["Grupos"][s][i]["Horario"][t],dict12["Grupos"][s][i]["Numero de copias"])
	end
		cont+=dict12["Grupos"][s][i]["Numero de copias"]
	end
	Mrel[s,:]=dict12["Matriz grupos"][s]
	end
	Politica=dict12["Politica de testeo"]
	TestSyn=dict12["Testear sintomaticos"]
	infopar=infopar[1:cont2-1]

	if Politica!="No testear"
		Testrand=dict12["testeo random"]
		Testnorand=dict12["testeo no random"]

		if length(Testrand)>=1
			test="Si"
			Testgrupo=zeros(S)
			tpool=Testrand["tpool"]
			Testgrupo=Testrand["Porcentaje a testear por grupo"]
			Variar=Testrand["Variar"]
			if Variar==0
				quienes="Ambos"
				distribuir="No"
			end
			if Variar==1
				quienes="Ambos"
				distribuir="Si"
			end
			if Variar==2
				quienes="Trabajando"
				distribuir="No"
			end
			if Variar==3
				quienes="No Trabajando"
				distribuir="No"
			end
			Diastest=zeros(T)
			for t in 1:T
				Diastest[t]=Testrand["Dias de testeo"][t]
			end
		end
		if length(Testnorand)>=1
			test="No"
			tpool=Testnorand["tpool"]
			Diastest=zeros(Integrantes,T)
			for t in 1:T
				Diastest[:,t]=Testnorand["Dias de testeo"][t]
			end
			variar=0
			quienes=[]
			distribuir=[]
			Testgrupo=ones(length(Group[1,:]))
		end

	else
		test="No"
		Testrand=[]
		Testnorand=[]
		tpool=0
		Variar=[]
		Diastest=zeros(N,T)
		quienes=[]
		distribuir=[]
		Testgrupo=zeros(length(Group[1,:]))

	end
	Repeticiones=dict12["Repeticiones"]
	Porasint=dict12["Porcentaje asintomaticos"]
	Probfp=dict12["Probabilidad falso positivo"]
	Cuarentena=dict12["Cuarentena"]
	Cuaren=Cuarentena["Tipo"]
	if Cuaren=="solo"
		Diasatras=0
		Diascuarentena=0

	else
		Diasatras=Cuarentena["Dias atras"]
		Diascuarentena=Cuarentena["Dias cuarentena"]
	end


return N, Group,T,Horario,Mrel, Diastest,Testgrupo,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,tpool,test,quienes,Cuaren,Diascuarentena,Diasatras,TestSyn,distribuir,NombreGrupo,NombreHorario, NameP, infopar


end

####################################################################################################################
function graficar(arrays,Group,ggrafico,nombreg,NombresP,titulo,x,y,T)
	L=length(Group[:,1])
	Indices=zeros(L)
	grafico=[]
	for g in ggrafico
		Indices+=Group[:,g]
	end
	Ind2=Int.(setdiff(Indices.*collect(1:L),[0]))
	for a in 1:length(arrays)
		arr=zeros(T)
		if length(arrays[a][1])==T
			arr=arrays[a][1]
		else
			for j in Ind2
				arr+=arrays[a][1][j,:]
			end
		end
		if (titulo=="Tests Totales") | (titulo=="Tests Sintomáticos") | (titulo=="Tests Aleatorios")
			if a==1
				grafico=scatter(1:T,arr,label=NombresP[a],lw=3)
			else
				scatter!(1:T,arr,label=NombresP[a],lw=3)
			end
		else
			if a==1
				grafico=plot(1:T,arr,label=NombresP[a],lw=3)
			else
				plot!(1:T,arr,label=NombresP[a],lw=3)
			end
		end
	end
	if (titulo=="Tests Totales") | (titulo=="Tests Sintomáticos") | (titulo=="Tests Aleatorios")
		title!(titulo)
	else
		title!(titulo*" "*string(nombreg)  )
	end
	xlabel!(x)
	ylabel!(y)
	return grafico
end




#######################################################################################################################




	###############################################################################
	#Se realizan las 4 politicas
#politicas=4 #Numero de politicas
#NombresP=["No testear","Simple, Todos los dias","Comparación","Simple, Cada 10 dias"]#Nombres politicas

Caso=string(readdir("Leer")[1])
	mv("Leer/"*Caso,"Guardar/"*Caso,force=true)
dict12 = Dict()
ff=1
for json in filter(x -> endswith(x, ".json"), readdir("Guardar/"*Caso))
global N,Group,T,Horario,Mrel, Diastest,Testgrupo,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,tpool,test,quienes,Cuaren,Diascuarentena,Diasatras,TestSyn,distribuir,graficogrupo,nombregrafico,nombregrafico15,t_peak,scalar_asint,peak,NombreGrupo,NombreHorario, NameP, infopar
open("Guardar/"*Caso*"/"*string(json), "r") do f
    global dict12
    dicttxt1 = read(f,String)  # file information to string
    dict12=JSON.parse(dicttxt1)  # parse and transform data
end

if ff==1
N,Group,T,Horario,Mrel, Diastest,Testgrupo,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,tpool,test,quienes,Cuaren,Diascuarentena,Diasatras,TestSyn,distribuir,NombreGrupo,NombreHorario,NameP,infopar=leer(dict12)
t_peak=ones(N)*24                 #Tiempo peak de contagio
scalar_asint=1          #Escalar para los asintomaticos
#peak=[(ones(100)*0.01)' (ones(20)*0.2)']' #Prob peak de infeccion
peak=ones(N)*0.2 #Prob peak de infeccion
graficogrupo=[[[i for i in 1:length(Group[1,:])]']  [[i] for i in 1:length(Group[1,:])]']
nomb=["a" for i in 1:length(Group[1,:])]
for i in 1:length(Group[1,:])
	j=findfirst(Group[:,i].==1)
	nomb[i]=NombreGrupo[j]
end
nombregrafico=[["Total"]' [i for i in nomb]']'

	global ff, base_mQua,base_mQua2, base_mInf,base_mInf2,base_mInf3, base_mNFp, base_T, base_T2,base_mSy,base_maxInf,base_maxQua, base_maxSy, base_Infect, base_VInf
base_mQua,base_mQua2, base_mInf,base_mInf2,base_mInf3, base_mNFp, base_T, base_T2,base_mSy,base_maxInf,base_maxQua, base_maxSy, base_Infect, base_VInf= simulation(N,Group,T,Horario,Mrel, Diastest,Testgrupo,peak,t_peak,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,tpool,test,quienes,Cuaren,Diascuarentena,Diasatras,scalar_asint,TestSyn,distribuir)
end
if ff==2
	global N2,Group2,T2,Horario2,Mrel2, Diastest2,Testgrupo2,p_int2,p_ext2,Porasint2,Probfp2,Repeticiones2,Politica2,tpool2,test2,quienes2,Cuaren2,Diascuarentena2,Diasatras2,TestSyn2,distribuir2,NombreGrupo2,NombreHorario2,NameP2,infopar2
	N2,Group2,T2,Horario2,Mrel2, Diastest2,Testgrupo2,p_int2,p_ext2,Porasint2,Probfp2,Repeticiones2,Politica2,tpool2,test2,quienes2,Cuaren2,Diascuarentena2,Diasatras2,TestSyn2,distribuir2,NombreGrupo2,NombreHorario2,NameP2,infopar2=leer(dict12)

	global ideal_mQua,ideal_mQua2, ideal_mInf,ideal_mInf2,ideal_mInf3, ideal_mNFp, ideal_T,ideal_T2, ideal_mSy,ideal_maxInf,ideal_maxQua, ideal_maxSy, ideal_Infect, ideal_VInf
ideal_mQua,ideal_mQua2, ideal_mInf,ideal_mInf2,ideal_mInf3, ideal_mNFp, ideal_T,ideal_T2, ideal_mSy,ideal_maxInf,ideal_maxQua, ideal_maxSy, ideal_Infect, ideal_VInf= simulation(N2,Group2,T2,Horario2,Mrel2, Diastest2,Testgrupo2,peak,t_peak,p_int2,p_ext2,Porasint2,Probfp2,Repeticiones2,Politica2,tpool2,test2,quienes2,Cuaren2,Diascuarentena2,Diasatras2,scalar_asint,TestSyn2,distribuir2)
end
if ff==3
	global N3,Group3,T3,Horario3,Mrel3, Diastest3,Testgrupo3,p_int3,p_ext3,Porasint3,Probfp3,Repeticiones3,Politica3,tpool3,test3,quienes3,Cuaren3,Diascuarentena3,Diasatras3,TestSyn3,distribuir3,NombreGrupo3,NombreHorario3,NameP3,infopar3
	N3,Group3,T3,Horario3,Mrel3, Diastest3,Testgrupo3,p_int3,p_ext3,Porasint3,Probfp3,Repeticiones3,Politica3,tpool3,test3,quienes3,Cuaren3,Diascuarentena3,Diasatras3,TestSyn3,distribuir3,NombreGrupo3,NombreHorario3,NameP3,infopar3=leer(dict12)

	global f1_mQua,f1_mQua2, f1_mInf,f1_mInf2,f1_mInf3, f1_mNFp, f1_T, f1_T2, f1_mSy,f1_maxInf,f1_maxQua, f1_maxSy,f1_Infect, f1_VInf
f1_mQua,f1_mQua2, f1_mInf,f1_mInf2,f1_mInf3, f1_mNFp, f1_T, f1_T2, f1_mSy,f1_maxInf,f1_maxQua, f1_maxSy,f1_Infect, f1_VInf = simulation(N3,Group3,T3,Horario3,Mrel3, Diastest3,Testgrupo3,peak,t_peak,p_int3,p_ext3,Porasint3,Probfp3,Repeticiones3,Politica3,tpool3,test3,quienes3,Cuaren3,Diascuarentena3,Diasatras3,scalar_asint,TestSyn3,distribuir3)
end
if ff==4
	global N4,Group4,T4,Horario4,Mrel4, Diastest4,Testgrupo4,p_int4,p_ext4,Porasint4,Probfp4,Repeticiones4,Politica4,tpool4,test4,quienes4,Cuaren4,Diascuarentena4,Diasatras4,TestSyn4,distribuir4,graficogrupo4,nombregrafico4,NombreGrupo4,NombreHorario4,NameP4,infopar4
	N4,Group4,T4,Horario4,Mrel4, Diastest4,Testgrupo4,p_int4,p_ext4,Porasint4,Probfp4,Repeticiones4,Politica4,tpool4,test4,quienes4,Cuaren4,Diascuarentena4,Diasatras4,TestSyn4,distribuir4,NombreGrupo4,NombreHorario4,NameP4,infopar4=leer(dict12)


	global fr_mQua,fr_mQua2, fr_mInf,fr_mInf2,fr_mInf3, fr_mNFp, fr_T, fr_T2,fr_mSy,fr_maxInf,fr_maxQua, fr_maxSy, fr_Infect, fr_VInf
fr_mQua,fr_mQua2, fr_mInf,fr_mInf2,fr_mInf3, fr_mNFp, fr_T, fr_T2,fr_mSy,fr_maxInf,fr_maxQua, fr_maxSy, fr_Infect, fr_VInf = simulation(N4,Group4,T4,Horario4,Mrel4, Diastest4,Testgrupo4,peak,t_peak,p_int4,p_ext4,Porasint4,Probfp4,Repeticiones4,Politica4,tpool4,test4,quienes4,Cuaren4,Diascuarentena4,Diasatras4,scalar_asint,TestSyn4,distribuir4)
end
ff+=1
end
###############################################################################
###############Graficos y PDF########################################################

####################PDF##############################################################
if ff==2
	arrays=[[base_mQua2]]
	Coma=["."]
	Coma2=[[]]
	NombresP=[NameP]
	Tipos=[Politica]
	Tpool=[tpool]
	GrupoQua=[Cuaren]
	if TestSyn=="Si"
		TestSyn="Sí"
	end
	TestS=[TestSyn]
	TestP=[sum(base_T)]
	TestP2=[sum(base_T2)]
	TestTotalP=TestP.+TestP2
	InfectadosP=[sum(base_VInf)]
	Infectados2P=[sum(base_mInf3)]
	Infectados3P=[sum(base_mInf)]
	QuaP=[sum(base_mQua)]
	maxInfectadosP=[maximum(sum(base_VInf,dims=1))]
	dmaxInfectadosP=[findfirst(isequal(maximum(sum(base_VInf,dims=1))),sum(base_VInf,dims=1))[2]]
	CuarentenapeakP=[maximum(sum(base_mQua2,dims=1))]
	dCuarentenapeakP=[findfirst(isequal(maximum(sum(base_mQua2,dims=1))),sum(base_mQua2,dims=1))[2]]
	InfectadospeakP=[maximum(sum(base_mInf3,dims=1))]
	dInfectadospeakP=[findfirst(isequal(maximum(sum(base_mInf3,dims=1))),sum(base_mInf3,dims=1))[2]]
	InfectadosTpeakP=[maximum(sum(base_mInf,dims=1))]
	dInfectadosTpeakP=[findfirst(isequal(maximum(sum(base_mInf,dims=1))),sum(base_mInf,dims=1))[2]]
	if length(findall(Diastest.==1))!=0
		DiasTest=findfirst(Diastest.==1)
		auxiliar=Diastest
		auxiliar[DiasTest]=0
		if length(findall(auxiliar.==1))!=0
			FrecTest=[findfirst(auxiliar.==1)-DiasTest-1]
			DiasTest=[DiasTest]
		else
			FrecTest=["No aplica"]
		end
	else
		DiasTest=["No aplica"]
		FrecTest=["No aplica"]
	end
	A=[Dict()]
end
if ff==3
	arrays=[[base_mQua2],[ideal_mQua2]]
	Coma=[" y ", "."]
	Coma2=[" y ",[]]
	NombresP=[NameP,NameP2]
	Tipos=[Politica,Politica2]
	Tpool=[tpool,tpool2]
	GrupoQua=[Cuaren,Cuaren2]
	if TestSyn=="Si"
		TestSyn="Sí"
	end
	if TestSyn2=="Si"
		TestSyn2="Sí"
	end
	TestS=[TestSyn,TestSyn2]
	TestP=[sum(base_T),sum(ideal_T)]
	TestP2=[sum(base_T2),sum(ideal_T2)]
	TestTotalP=TestP.+TestP2
	InfectadosP=[sum(base_VInf),sum(ideal_VInf)]
	Infectados2P=[sum(base_mInf3),sum(ideal_mInf3)]
	Infectados3P=[sum(base_mInf),sum(ideal_mInf)]
	QuaP=[sum(base_mQua),sum(ideal_mQua)]
	maxInfectadosP=[maximum(sum(base_VInf,dims=1)),maximum(sum(ideal_VInf,dims=1))]
	dmaxInfectadosP=[findfirst(isequal(maximum(sum(base_VInf,dims=1))),sum(base_VInf,dims=1))[2],findfirst(isequal(maximum(sum(ideal_VInf,dims=1))),sum(ideal_VInf,dims=1))[2]]
	CuarentenapeakP=[maximum(sum(base_mQua2,dims=1)),maximum(sum(ideal_mQua2,dims=1))]
	dCuarentenapeakP=[findfirst(isequal(maximum(sum(base_mQua2,dims=1))),sum(base_mQua2,dims=1))[2],findfirst(isequal(maximum(sum(ideal_mQua2,dims=1))),sum(ideal_mQua2,dims=1))[2]]
	InfectadospeakP=[maximum(sum(base_mInf3,dims=1)),maximum(sum(ideal_mInf3,dims=1))]
	dInfectadospeakP=[findfirst(isequal(maximum(sum(base_mInf3,dims=1))),sum(base_mInf3,dims=1))[2],findfirst(isequal(maximum(sum(ideal_mInf3,dims=1))),sum(ideal_mInf3,dims=1))[2]]
	InfectadosTpeakP=[maximum(sum(base_mInf,dims=1)),maximum(sum(ideal_mInf,dims=1))]
	dInfectadosTpeakP=[findfirst(isequal(maximum(sum(base_mInf,dims=1))),sum(base_mInf,dims=1))[2],findfirst(isequal(maximum(sum(ideal_mInf,dims=1))),sum(ideal_mInf,dims=1))[2]]
	A=[Dict(),Dict()]
	if length(findall(Diastest.==1))!=0
		DiasTest1=findfirst(Diastest.==1)
		auxiliar=Diastest
		auxiliar[DiasTest1]=0
		if length(findall(auxiliar.==1))!=0
			FrecTest1=findfirst(Diastest.==1)-DiasTest1-1
		else
			FrecTest1="No aplica"
		end
	else
		DiasTest1="No aplica"
		FrecTest1="No aplica"
	end
	if length(findall(Diastest2.==1))!=0
		DiasTest2=findfirst(Diastest2.==1)
		auxiliar=Diastest2
		auxiliar[DiasTest2]=0
		if length(findall(auxiliar.==1))!=0
			FrecTest2=findfirst(Diastest2.==1)-DiasTest2-1
		else
			FrecTest2="No aplica"
		end
	else
		DiasTest2="No aplica"
		FrecTest2="No aplica"
	end
DiasTest=[DiasTest1,DiasTest2]
FrecTest=[FrecTest1,FrecTest2]
end
if ff==4
	arrays=[[base_mQua2],[ideal_mQua2],[f1_mQua2]]
	Coma=[", "," y ","."]
	Coma2=[", "," y ",[]]
	NombresP=[NameP,NameP2,NameP3]
	Tipos=[Politica,Politica2,Politica3]
	Tpool=[tpool,tpool2,tpool3]
	GrupoQua=[Cuaren,Cuaren2,Cuaren3]
	if TestSyn=="Si"
		TestSyn="Sí"
	end
	if TestSyn2=="Si"
		TestSyn2="Sí"
	end
	if TestSyn3=="Si"
		TestSyn3="Sí"
	end
	TestS=[TestSyn,TestSyn2,TestSyn3]
	TestP=[sum(base_T),sum(ideal_T),sum(f1_T)]
	TestP2=[sum(base_T2),sum(ideal_T2),sum(f1_T2)]
	TestTotalP=TestP.+TestP2
	InfectadosP=[sum(base_VInf),sum(ideal_VInf),sum(f1_VInf)]
	Infectados2P=[sum(base_mInf3),sum(ideal_mInf3),sum(f1_mInf3)]
	Infectados3P=[sum(base_mInf),sum(ideal_mInf),sum(f1_mInf)]
	QuaP=[sum(base_mQua),sum(ideal_mQua),sum(f1_mQua)]
	maxInfectadosP=[maximum(sum(base_VInf,dims=1)),maximum(sum(ideal_VInf,dims=1)),maximum(sum(f1_VInf,dims=1))]
	dmaxInfectadosP=[findfirst(isequal(maximum(sum(base_VInf,dims=1))),sum(base_VInf,dims=1))[2],findfirst(isequal(maximum(sum(ideal_VInf,dims=1))),sum(ideal_VInf,dims=1))[2],findfirst(isequal(maximum(sum(f1_VInf,dims=1))),sum(f1_VInf,dims=1))[2]]
	CuarentenapeakP=[maximum(sum(base_mQua2,dims=1)),maximum(sum(ideal_mQua2,dims=1)),maximum(sum(f1_mQua2,dims=1))]
	dCuarentenapeakP=[findfirst(isequal(maximum(sum(base_mQua2,dims=1))),sum(base_mQua2,dims=1))[2],findfirst(isequal(maximum(sum(ideal_mQua2,dims=1))),sum(ideal_mQua2,dims=1))[2],findfirst(isequal(maximum(sum(f1_mQua2,dims=1))),sum(f1_mQua2,dims=1))[2]]
	InfectadospeakP=[maximum(sum(base_mInf3,dims=1)),maximum(sum(ideal_mInf3,dims=1)),maximum(sum(f1_mInf3,dims=1))]
	dInfectadospeakP=[findfirst(isequal(maximum(sum(base_mInf3,dims=1))),sum(base_mInf3,dims=1))[2],findfirst(isequal(maximum(sum(ideal_mInf3,dims=1))),sum(ideal_mInf3,dims=1))[2],findfirst(isequal(maximum(sum(f1_mInf3,dims=1))),sum(f1_mInf3,dims=1))[2]]
	InfectadosTpeakP=[maximum(sum(base_mInf,dims=1)),maximum(sum(ideal_mInf,dims=1)),maximum(sum(f1_mInf,dims=1))]
	dInfectadosTpeakP=[findfirst(isequal(maximum(sum(base_mInf,dims=1))),sum(base_mInf,dims=1))[2],findfirst(isequal(maximum(sum(ideal_mInf,dims=1))),sum(ideal_mInf,dims=1))[2],findfirst(isequal(maximum(sum(f1_mInf,dims=1))),sum(f1_mInf,dims=1))[2]]
	A=[Dict(),Dict(),Dict()]
	if length(findall(Diastest.==1))!=0
		DiasTest1=findfirst(Diastest.==1)
		auxiliar=Diastest
		auxiliar[DiasTest1]=0
		if length(findall(auxiliar.==1))!=0
			FrecTest1=findfirst(Diastest.==1)-DiasTest1-1
		else
			FrecTest1="No aplica"
		end
	else
		DiasTest1="No aplica"
		FrecTest1="No aplica"
	end
	if length(findall(Diastest2.==1))!=0
		DiasTest2=findfirst(Diastest2.==0)
		auxiliar=Diastest2
		auxiliar[DiasTest2]=0
		if length(findall(auxiliar.==1))!=0
			FrecTest2=findfirst(Diastest2.==1)-DiasTest2-1
		else
			FrecTest2="No aplica"
		end
	else
		DiasTest2="No aplica"
		FrecTest2="No aplica"
	end
	if length(findall(Diastest3.==1))!=0
		DiasTest3=findfirst(Diastest3.==1)
		auxiliar=Diastest3
		auxiliar[DiasTest3]=0
		if length(findall(auxiliar.==1))!=0
			FrecTest3=findfirst(Diastest3.==1)-DiasTest3-1
		else
			FrecTest3="No aplica"
		end
	else
		DiasTest3="No aplica"
		FrecTest3="No aplica"

	end
DiasTest=[DiasTest1,DiasTest2,DiasTest3]
FrecTest=[FrecTest1,FrecTest2,FrecTest3]
end
if ff==5
	arrays=[[base_mQua2],[ideal_mQua2],[f1_mQua2],[fr_mQua2]]
	Coma=[", ",", "," y ","."]
	Coma2=[", ",", "," y ",[]]
	NombresP=[NameP,NameP2,NameP3,NameP4]
	Tipos=[Politica,Politica2,Politica3,Politica4]
	Tpool=[tpool,tpool2,tpool3,tpool4]
	GrupoQua=[Cuaren,Cuaren2,Cuaren3,Cuaren4]
	if TestSyn=="Si"
		TestSyn="Sí"
	end
	if TestSyn2=="Si"
		TestSyn2="Sí"
	end
	if TestSyn3=="Si"
		TestSyn3="Sí"
	end
	if TestSyn4=="Si"
		TestSyn4="Sí"
	end
	TestS=[TestSyn,TestSyn2,TestSyn3,TestSyn4]
	TestP=[sum(base_T),sum(ideal_T),sum(f1_T),sum(fr_T)]
	TestP2=[sum(base_T2),sum(ideal_T2),sum(f1_T2),sum(fr_T2)]
	TestTotalP=TestP.+TestP2
	InfectadosP=[sum(base_VInf),sum(ideal_VInf),sum(f1_VInf),sum(fr_VInf)]
	Infectados2P=[sum(base_mInf3),sum(ideal_mInf3),sum(f1_mInf3),sum(fr_mInf3)]
	Infectados3P=[sum(base_mInf),sum(ideal_mInf),sum(f1_mInf),sum(fr_mInf)]
	QuaP=[sum(base_mQua),sum(ideal_mQua),sum(f1_mQua),sum(fr_mQua)]
	maxInfectadosP=[maximum(sum(base_VInf,dims=1)),maximum(sum(ideal_VInf,dims=1)),maximum(sum(f1_VInf,dims=1)),maximum(sum(fr_VInf,dims=1))]
	dmaxInfectadosP=[findfirst(isequal(maximum(sum(base_VInf,dims=1))),sum(base_VInf,dims=1))[2],findfirst(isequal(maximum(sum(ideal_VInf,dims=1))),sum(ideal_VInf,dims=1))[2],findfirst(isequal(maximum(sum(f1_VInf,dims=1))),sum(f1_VInf,dims=1))[2],
	findfirst(isequal(maximum(sum(fr_VInf,dims=1))),sum(fr_VInf,dims=1))[2]]
	CuarentenapeakP=[maximum(sum(base_mQua2,dims=1)),maximum(sum(ideal_mQua2,dims=1)),maximum(sum(f1_mQua2,dims=1)),maximum(sum(fr_mQua2,dims=1))]
	dCuarentenapeakP=[findfirst(isequal(maximum(sum(base_mQua2,dims=1))),sum(base_mQua2,dims=1))[2],findfirst(isequal(maximum(sum(ideal_mQua2,dims=1))),sum(ideal_mQua2,dims=1))[2],findfirst(isequal(maximum(sum(f1_mQua2,dims=1))),sum(f1_mQua2,dims=1))[2],
	findfirst(isequal(maximum(sum(fr_mQua2,dims=1))),sum(fr_mQua2,dims=1))[2]]
	InfectadospeakP=[maximum(sum(base_mInf3,dims=1)),maximum(sum(ideal_mInf3,dims=1)),maximum(sum(f1_mInf3,dims=1)),maximum(sum(fr_mInf3,dims=1))]
	dInfectadospeakP=[findfirst(isequal(maximum(sum(base_mInf3,dims=1))),sum(base_mInf3,dims=1))[2],findfirst(isequal(maximum(sum(ideal_mInf3,dims=1))),sum(ideal_mInf3,dims=1))[2],findfirst(isequal(maximum(sum(f1_mInf3,dims=1))),sum(f1_mInf3,dims=1))[2],
	findfirst(isequal(maximum(sum(fr_mInf3,dims=1))),sum(fr_mInf3,dims=1))[2]]
	InfectadosTpeakP=[maximum(sum(base_mInf,dims=1)),maximum(sum(ideal_mInf,dims=1)),maximum(sum(f1_mInf,dims=1)),maximum(sum(fr_mInf,dims=1))]
	dInfectadosTpeakP=[findfirst(isequal(maximum(sum(base_mInf,dims=1))),sum(base_mInf,dims=1))[2],findfirst(isequal(maximum(sum(ideal_mInf,dims=1))),sum(ideal_mInf,dims=1))[2],findfirst(isequal(maximum(sum(f1_mInf,dims=1))),sum(f1_mInf,dims=1))[2],
	findfirst(isequal(maximum(sum(fr_mInf,dims=1))),sum(fr_mInf,dims=1))[2]]
	A=[Dict(),Dict(),Dict(),Dict()]
	if length(findall(Diastest.==1))!=0
		DiasTest1=findfirst(Diastest.==1)
		auxiliar=Diastest
		auxiliar[Int(DiasTest1)]=0
		if length(findall(auxiliar.==1))!=0
			FrecTest1=findfirst(Diastest.==1)-DiasTest1-1
		else
			FrecTest1="No aplica"
		end
	else
		DiasTest1="No aplica"
		FrecTest1="No aplica"
	end
	if length(findall(Diastest2.==1))!=0
		DiasTest2=findfirst(Diastest2.==1)
		auxiliar=Diastest2
		auxiliar[Int(DiasTest2)]=0
		if length(findall(auxiliar.==1))!=0
			FrecTest2=findfirst(Diastest2.==1)-DiasTest2-1
		else
			FrecTest2="No aplica"
		end
	else
		DiasTest2="No aplica"
		FrecTest2="No aplica"
	end
	if length(findall(Diastest3.==1))!=0
		DiasTest3=findfirst(Diastest3.==1)
		auxiliar=Diastest3
		auxiliar[Int(DiasTest3)]=0
		if length(findall(auxiliar.==1))!=0
			FrecTest3=findfirst(Diastest3.==1)-DiasTest3-1
		else
			FrecTest3="No aplica"
		end
	else
		DiasTest3="No aplica"
		FrecTest3="No aplica"
	end

	if length(findall(Diastest4.==1))!=0
		DiasTest4=findfirst(Diastest4.==1)
		auxiliar=Diastest4
		auxiliar[Int(DiasTest4)]=0
		if length(findall(auxiliar.==1))!=0
			FrecTest4=findfirst(Diastest4.==1)-DiasTest4-1
		else
			FrecTest4="No aplica"
		end
	else
		DiasTest4="No aplica"
		FrecTest4="No aplica"
	end
DiasTest=[DiasTest1,DiasTest2,DiasTest3,DiasTest4]
FrecTest=[FrecTest1,FrecTest2,FrecTest3,FrecTest4]

end

for cont in 1:length(NombresP)
A[cont]=Dict("Politica"=>NombresP[cont],"Test"=>Int.(round.(TestP[cont])),"Test2"=>Int.(round.(TestP2[cont])),"TestTotal"=>Int.(round.(TestTotalP[cont])),
"Infectados"=>Int.(round.(InfectadosP[cont])),"maxInfectados"=>Int.(round.(maxInfectadosP[cont])),"dmaxInfectados"=>dmaxInfectadosP[cont],
"PeakQua"=>Int.(round.(CuarentenapeakP[cont])),"dPeakQua"=>dCuarentenapeakP[cont],"Infectados2"=>Int.(round.(Infectados2P[cont])),"Infectados2T"=>Int.(round.(Infectados2P[cont]/T)),
"PeakInf"=>Int.(round.(InfectadospeakP[cont])),"dPeakInf"=>dInfectadospeakP[cont],"Infectados3"=>Int.(round.(Infectados3P[cont])),"Infectados3T"=>Int.(round.(Infectados3P[cont]/T)),
"PeakInfT"=>Int.(round.(InfectadosTpeakP[cont])),"dPeakInfT"=>dInfectadosTpeakP[cont],"Coma"=>Coma[cont],"Coma2"=>Coma2[cont], "Tipo"=>Tipos[cont],"Tpool"=>Tpool[cont],
"GrupoQua"=>GrupoQua[cont],"TestS"=>TestS[cont],"DiasTest"=>DiasTest[cont],"FrecTest"=>FrecTest[cont],"Qua"=>Int.(round.(QuaP[cont])),"QuaT"=>Int.(round.(QuaP[cont]/T))   )
end
AB=[Dict()]
AB[1]["Caso"]=Caso
AB[1]["NumP"]=ff-1
AB[1]["NumG"]=[string(i) for i in 2:(length(Group[1,:])+1)]
AB[1]["G"]=string(Int(floor(length(Group[1,:])/2)))
AB[1]["N"]=N
AB[1]["T"]=T
AB[1]["Rep"]=Repeticiones
AB[1]["Porc Asint"]=Porasint
AB[1]["Superc"]=prod([" c " for i in 1:length(Group[1,:])])*"|"
AB[1]["Superc2"]=prod(["c|" for i in 1:(ff-1)])
AB[1]["MInfectados"]=maximum(Int.(round.(InfectadosP)))
AB[1]["MInfectadosP"]=NombresP[findfirst(Int.(round.(InfectadosP)).==maximum(Int.(round.(InfectadosP))))]
AB[1]["mInfectados"]=minimum(Int.(round.(InfectadosP)))
AB[1]["mInfectadosP"]=NombresP[findfirst(Int.(round.(InfectadosP)).==minimum(Int.(round.(InfectadosP))))]

AB[1]["MInfectadospeak"]=maximum(Int.(round.(InfectadospeakP)))
AB[1]["MInfectadospeakP"]=NombresP[findfirst(Int.(round.(InfectadospeakP)).==maximum(Int.(round.(InfectadospeakP))))]
AB[1]["MInfectadospeakD"]=dInfectadospeakP[findfirst(Int.(round.(InfectadospeakP)).==maximum(Int.(round.(InfectadospeakP))))]

AB[1]["mInfectadospeak"]=minimum(Int.(round.(InfectadospeakP)))
AB[1]["mInfectadospeakP"]=NombresP[findfirst(Int.(round.(InfectadospeakP)).==minimum(Int.(round.(InfectadospeakP))))]
AB[1]["mInfectadospeakD"]=dInfectadospeakP[findfirst(Int.(round.(InfectadospeakP)).==minimum(Int.(round.(InfectadospeakP))))]


AB[1]["MInfectadosTpeak"]=maximum(Int.(round.(InfectadosTpeakP)))
AB[1]["MInfectadosTpeakP"]=NombresP[findfirst(Int.(round.(InfectadosTpeakP)).==maximum(Int.(round.(InfectadosTpeakP))))]
AB[1]["MInfectadosTpeakD"]=dInfectadosTpeakP[findfirst(Int.(round.(InfectadosTpeakP)).==maximum(Int.(round.(InfectadosTpeakP))))]

AB[1]["mInfectadosTpeak"]=minimum(Int.(round.(InfectadosTpeakP)))
AB[1]["mInfectadosTpeakP"]=NombresP[findfirst(Int.(round.(InfectadosTpeakP)).==minimum(Int.(round.(InfectadosTpeakP))))]
AB[1]["mInfectadosTpeakD"]=dInfectadosTpeakP[findfirst(Int.(round.(InfectadosTpeakP)).==minimum(Int.(round.(InfectadosTpeakP))))]


AB[1]["mTest"]=minimum(Int.(round.(TestP)))
AB[1]["mTestP"]=NombresP[findfirst(Int.(round.(TestP)).==minimum(Int.(round.(TestP))))]
AB[1]["MTest"]=maximum(Int.(round.(TestP)))
AB[1]["MTestP"]=NombresP[findfirst(Int.(round.(TestP)).==maximum(Int.(round.(TestP))))]


AB[1]["mTest2"]=minimum(Int.(round.(TestP2)))
AB[1]["mTestP2"]=NombresP[findfirst(Int.(round.(TestP2)).==minimum(Int.(round.(TestP2))))]
AB[1]["MTest2"]=maximum(Int.(round.(TestP2)))
AB[1]["MTestP2"]=NombresP[findfirst(Int.(round.(TestP2)).==maximum(Int.(round.(TestP2))))]



AB[1]["mTestTotal"]=minimum(Int.(round.(TestTotalP)))
AB[1]["mTestTotalP"]=NombresP[findfirst(Int.(round.(TestTotalP)).==minimum(Int.(round.(TestTotalP))))]
AB[1]["MTestTotal"]=maximum(Int.(round.(TestTotalP)))
AB[1]["MTestTotalP"]=NombresP[findfirst(Int.(round.(TestTotalP)).==maximum(Int.(round.(TestTotalP))))]




AB[1]["MCuarentenapeak"]=maximum(Int.(round.(CuarentenapeakP)))
AB[1]["MCuarentenapeakP"]=NombresP[findfirst(Int.(round.(CuarentenapeakP)).==maximum(Int.(round.(CuarentenapeakP))))]
AB[1]["MCuarentenapeakD"]=dCuarentenapeakP[findfirst(Int.(round.(CuarentenapeakP)).==maximum(Int.(round.(CuarentenapeakP))))]

AB[1]["mCuarentenapeak"]=minimum(Int.(round.(CuarentenapeakP)))
AB[1]["mCuarentenapeakP"]=NombresP[findfirst(Int.(round.(CuarentenapeakP)).==minimum(Int.(round.(CuarentenapeakP))))]
AB[1]["mCuarentenapeakD"]=dCuarentenapeakP[findfirst(Int.(round.(CuarentenapeakP)).==minimum(Int.(round.(CuarentenapeakP))))]


C=sort(collect(Set(NombreGrupo)))
if length(C)==1
	AB[1]["Nomb Grupos"]=C[1]
end
if length(C)==2
	AB[1]["Nomb Grupos"]=C[1]*" y "*C[2]
end
if length(C)>=3
		C2=C[1:end-2]
		C3=[i*", " for i in C2]
		AB[1]["Nomb Grupos"]=prod(C3)*C[end-1]*" y "*C[end]
end

A2=[Dict() for j in 1:length(Group[1,:])]
nomb=["a" for i in 1:length(Group[1,:])]
for i in 1:length(Group[1,:])
j=findfirst(Group[:,i].==1)
nomb[i]=NombreGrupo[j]
A2[i]["Nombres Grupos"]=nomb[i]
if ff==2
	A2[i]["Test Grupos"]=string.(Int.([Testgrupo[i]]*100)).*string(\).*"%"
end
if ff==3
	A2[i]["Test Grupos"]=string.(Int.([Testgrupo[i],Testgrupo2[i]]*100)).*string(\).*"%"
end
if ff==4
	A2[i]["Test Grupos"]=string.(Int.([Testgrupo[i],Testgrupo2[i],Testgrupo3[i]]*100)).*string(\).*"%"
end
if ff==5
	A2[i]["Test Grupos"]=string.(Int.([Testgrupo[i],Testgrupo2[i],Testgrupo3[i],Testgrupo4[i]]*100)).*string(\).*"%"
end
A2[i]["Matriz rel"]=string.(Int.(Mrel[i,:]*100)).*string(\).*"%"
#print(A2[i]["Matriz rel"])
#A2[i]["Matriz rel"]=Mrel[i,:]
#print(A2[i]["Matriz rel"])
A2[i]["Ind"]=i+1
end

D=sort(collect(Set(NombreHorario)))


A3=[Dict() for j in 1:N]
for i in 1:N
A3[i]["Nombres Grupos"]=NombreGrupo[i]
A3[i]["Nombres Horario"]=NombreHorario[i]
A3[i]["Prob Cont Int"]=@sprintf("%1.5f",p_int[i])
A3[i]["Prob Cont Ext"]=@sprintf("%1.5f",p_ext[i])
A3[i]["Conteo"]=i
end

A4=[Dict() for j in 1:length(infopar)]
for i in 1:length(infopar)
A4[i]["Copias"]=Int(infopar[i])
A4[i]["Nombres Grupos"]=NombreGrupo[Int(round(sum(infopar[1:i])))]
A4[i]["Nombres Horario"]=NombreHorario[Int(round(sum(infopar[1:i])))]
A4[i]["Prob Cont Int"]=@sprintf("%1.5f",p_int[Int(round(sum(infopar[1:i])))])
A4[i]["Prob Cont Ext"]=@sprintf("%1.5f",p_ext[Int(round(sum(infopar[1:i])))])
A4[i]["Conteo"]=i
end

A5=[Dict() for j in 1:length(D)]
for i in 1:length(D)
A5[i]["NHorario"]=D[i]
A5[i]["HorarioD"]=Int.(Horario[(2*(findfirst(D[i].==NombreHorario))-1),1:7])
A5[i]["HorarioN"]=Int.(Horario[2*(findfirst(D[i].==NombreHorario)),1:7])
A5[i]["Conteo"]=i
end








#############################################################
##########################Graficar###########################
	titulo="Cuarentena"
	x="Dias"
	y="Personas"
for i in 1:length(graficogrupo)
	ggrafico=graficogrupo[i]
	nombreg=nombregrafico[i]
	g=graficar(arrays,Group,ggrafico,nombreg,NombresP,titulo,x,y,T)
	plot(g)
	savefig("Guardar/"*Caso*"/"*titulo*string(i)*".pdf")
	g=[]
end

if ff==2
	arrays=[[base_mInf]]
end
if ff==3
	arrays=[[base_mInf],[ideal_mInf]]
end
if ff==4
	arrays=[[base_mInf],[ideal_mInf],[f1_mInf]]
end
if ff==5
	arrays=[[base_mInf],[ideal_mInf],[f1_mInf],[fr_mInf]]
end
titulo="Infectados Trabajando"
x="Dias"
y="Personas"
for i in 1:length(graficogrupo)
ggrafico=graficogrupo[i]
nombreg=nombregrafico[i]
g=graficar(arrays,Group,ggrafico,nombreg,NombresP,titulo,x,y,T)
plot(g)
savefig("Guardar/"*Caso*"/InfectadosTrabajando"*string(i)*".pdf")
g=[]
end


if ff==2
	arrays=[[base_mInf3]]
end
if ff==3
	arrays=[[base_mInf3],[ideal_mInf3]]
end
if ff==4
	arrays=[[base_mInf3],[ideal_mInf3],[f1_mInf3]]
end
if ff==5
	arrays=[[base_mInf3],[ideal_mInf3],[f1_mInf3],[fr_mInf3]]
end
titulo="Infectados"
x="Dias"
y="Personas"
for i in 1:length(graficogrupo)
ggrafico=graficogrupo[i]
nombreg=nombregrafico[i]
g=graficar(arrays,Group,ggrafico,nombreg,NombresP,titulo,x,y,T)
plot(g)
savefig("Guardar/"*Caso*"/"*titulo*string(i)*".pdf")
g=[]
end

if ff==2
	arrays=[[acu(base_VInf)]]
end
if ff==3
	arrays=[[acu(base_VInf)],[acu(ideal_VInf)]]
end
if ff==4
	arrays=[[acu(base_VInf)],[acu(ideal_VInf)],[acu(f1_VInf)]]
end
if ff==5
	arrays=[[acu(base_VInf)],[acu(ideal_VInf)],[acu(f1_VInf)],[acu(fr_VInf)]]
end
titulo="Infectados Acumulado"
x="Dias"
y="Personas"
for i in 1:length(graficogrupo)
ggrafico=graficogrupo[i]
nombreg=nombregrafico[i]
g=graficar(arrays,Group,ggrafico,nombreg,NombresP,titulo,x,y,T)
plot(g,legend=:bottomright)
savefig("Guardar/"*Caso*"/InfectadosAcumulado"*string(i)*".pdf")
g=[]
end





graficogrupo=[1:length(Group[1,:])]
nombregrafico=["Total"]


if ff==2
	arrays=[[base_T]]
end
if ff==3
	arrays=[[base_T],[ideal_T]]
end
if ff==4
	arrays=[[base_T],[ideal_T],[f1_T]]
end
if ff==5
	arrays=[[base_T],[ideal_T],[f1_T],[fr_T]]
end
titulo="Tests Aleatorios"
x="Dias"
y="Tests"
for i in 1:length(graficogrupo)
ggrafico=graficogrupo[i]
nombreg=nombregrafico[i]
g=graficar(arrays,Group,ggrafico,nombreg,NombresP,titulo,x,y,T)
plot(g)
savefig("Guardar/"*Caso*"/TestsAleatorios"*string(i)*".pdf")
g=[]
end


if ff==2
	arrays=[sum(base_T)]
end
if ff==3
	arrays=[sum(base_T),sum(ideal_T)]
end
if ff==4
	arrays=[sum(base_T),sum(ideal_T),sum(f1_T)]
end
if ff==5
	arrays=[sum(base_T),sum(ideal_T),sum(f1_T),sum(fr_T)]
end
titulo="Tests Aleatorios Histograma"
y="Tests"
for i in 1:length(graficogrupo)
ggrafico=graficogrupo[i]
nombreg=nombregrafico[i]
L=length(Group[:,1])
Indices=zeros(L)
for g in ggrafico
	Indices+=Group[:,g]
end
g=plot([i-1 for i in 2:ff]',arrays',seriestype=:bar,label=NombresP,xlims=(0.4,4.6),ylims=(0,maximum( [10,maximum(arrays)]  ) ))
auxi=[(i,arrays[i]/2,Int(round(arrays[i]))) for i in 1:length(arrays)]
plot!(; annotations=auxi)
title!(titulo)
ylabel!(y)
plot(g)
savefig("Guardar/"*Caso*"/TestsAleatoriosHistograma"*string(i)*".pdf")
g=[]
end



if ff==2
	arrays=[[base_T2]]
end
if ff==3
	arrays=[[base_T2],[ideal_T2]]
end
if ff==4
	arrays=[[base_T2],[ideal_T2],[f1_T2]]
end
if ff==5
	arrays=[[base_T2],[ideal_T2],[f1_T2],[fr_T2]]
end
titulo="Tests Sintomáticos"
x="Dias"
y="Tests"
for i in 1:length(graficogrupo)
ggrafico=graficogrupo[i]
nombreg=nombregrafico[i]
g=graficar(arrays,Group,ggrafico,nombreg,NombresP,titulo,x,y,T)
plot(g)
savefig("Guardar/"*Caso*"/TestsSintomaticos"*string(i)*".pdf")
g=[]
end
if ff==2
	arrays=[sum(base_T2)]
end
if ff==3
	arrays=[sum(base_T2),sum(ideal_T2)]
end
if ff==4
	arrays=[sum(base_T2),sum(ideal_T2),sum(f1_T2)]
end
if ff==5
	arrays=[sum(base_T2),sum(ideal_T2),sum(f1_T2),sum(fr_T2)]
end
titulo="Tests Sintomáticos Histograma"
y="Tests"
for i in 1:length(graficogrupo)
ggrafico=graficogrupo[i]
nombreg=nombregrafico[i]
L=length(Group[:,1])
Indices=zeros(L)
for g in ggrafico
	Indices+=Group[:,g]
end
g=plot([i-1 for i in 2:ff]',arrays',seriestype=:bar,label=NombresP,xlims=(0.4,4.6),ylims=(0,maximum( [10,maximum(arrays)]  ) ))
auxi=[(i,arrays[i]/2,Int(round(arrays[i]))) for i in 1:length(arrays)]
plot!(; annotations=auxi)
title!(titulo)
ylabel!(y)
plot(g)
savefig("Guardar/"*Caso*"/TestsSintomaticosHistograma"*string(i)*".pdf")
g=[]
end





if ff==2
	arrays=[[base_T]]
end
if ff==3
	arrays=[[base_T],[ideal_T]]
end
if ff==4
	arrays=[[base_T],[ideal_T],[f1_T]]
end
if ff==5
	arrays=[[base_T],[ideal_T],[f1_T],[fr_T]]
end
titulo="Tests Totales"
x="Dias"
y="Tests"
for i in 1:length(graficogrupo)
ggrafico=graficogrupo[i]
nombreg=nombregrafico[i]
g=graficar(arrays,Group,ggrafico,nombreg,NombresP,titulo,x,y,T)
plot(g)
savefig("Guardar/"*Caso*"/Tests"*string(i)*".pdf")
g=[]
end


if ff==2
	arrays=[sum(base_T)]
end
if ff==3
	arrays=[sum(base_T),sum(ideal_T)]
end
if ff==4
	arrays=[sum(base_T),sum(ideal_T),sum(f1_T)]
end
if ff==5
	arrays=[sum(base_T),sum(ideal_T),sum(f1_T),sum(fr_T)]
end
titulo="Tests Totales Histograma"
y="Tests"
for i in 1:length(graficogrupo)
ggrafico=graficogrupo[i]
nombreg=nombregrafico[i]
L=length(Group[:,1])
Indices=zeros(L)
for g in ggrafico
	Indices+=Group[:,g]
end
g=plot([i-1 for i in 2:ff]',arrays',seriestype=:bar,label=NombresP,xlims=(0.4,4.6),ylims=(0,maximum( [10,maximum(arrays)]  ) ))
auxi=[(i,arrays[i]/2,Int(round(arrays[i]))) for i in 1:length(arrays)]
plot!(; annotations=auxi)
title!(titulo)
ylabel!(y)
plot(g)
savefig("Guardar/"*Caso*"/TestsHistograma"*string(i)*".pdf")
g=[]
end
############################################################################



tmpl = mt"""

\documentclass[opre]{PDF/informs3}
\usepackage[utf8]{inputenc}
\usepackage[spanish]{babel}
\usepackage{multicol}
\usepackage[dvipsnames]{xcolor}
\usepackage{verbatim,wrapfig,latexsym,graphicx,psfrag,amsfonts,amsmath,amssymb}
\usepackage{hyperref,multirow,color,subfigure,comment,pdfsync,setspace,thmtools,thm-restate}
\usepackage{enumerate,enumitem,tikz,pdfpages,xspace,multirow,algorithmic,algorithm}
%\usepackage{longtable} para tablas largas
\usepackage{array}
\newcolumntype{L}[1]{>{\raggedright\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{C}[1]{>{\centering\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\newcolumntype{R}[1]{>{\raggedleft\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}


\newcommand{\E}{\mathbb{E}}
\newcommand{\Pb}{\mathbb{P}}
\newcommand{\Real}{\mathbb{R}}
\newcommand{\Natu}{\mathbb{N}}
\newcommand{\Int}{\mathbb{Z}}
\newcommand{\set}[1]{\left\{#1\right\}}
\newcommand{\eexp}[1]{e^{#1}}
\newcommand{\abs}[1]{\left|#1\right|}
\newcommand{\bra}[1]{\left(#1\right)}
\newcommand{\norm}[1]{\left\|#1\right\|}
\newcommand{\setnorm}[2]{\|#1\|_{#2}}
\newcommand{\ind}[1]{\mathbf{1}\set{#1}}
\newcommand{\sgn}[1]{\text{sign}\bra{#1}}
\newcommand{\crs}[2]{#1 \cdot #2}
\newcommand{\eqq}{:=}
\newcommand{\eps}{\varepsilon}
\newcommand{\tablehighlight}[1]{\underline{\textbf{#1}}}
\DeclareMathOperator{\inte}{int}
\DeclareMathOperator{\conv}{conv}
\DeclareMathOperator{\supp}{supp}
\DeclareMathOperator{\diag}{diag}
\DeclareMathOperator{\spn}{span}
\DeclareMathOperator{\vol}{vol}
\DeclareMathOperator{\var}{Var}

\newcommand{\Nscr}{\mathcal{N}}
\newcommand{\Tscr}{\mathcal{T}}
\newcommand{\Sscr}{\mathcal{S}}

\TheoremsNumberedThrough
\ECRepeatTheorems

\EquationsNumberedThrough


%\makeatletter
%\def\l@section#1#2{%
%\noindent\begin{minipage}[t]{\dimexpr(.333\linewidth)-2em\relax}%
%\def\numberline##1{##1: }\parfillskip0pt\relax
%#1 \mbox{}\dotfill #2\end{minipage}\linebreak[0]\hspace{2em plus 2em}}
%\makeatletter

\begin{document}

\RUNAUTHOR{Reporte ISCI}

\RUNTITLE{Cuarentena para equipos m\'edicos}

\TITLE{}

\ARTICLEAUTHORS{
\AUTHOR{}
}

\ABSTRACT{
...
}

\KEYWORDS{}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\begin{wrapfigure}{l}{4.5cm}
\includegraphics[width=4.5cm]{PDF/logo_ISCI.pdf}
\vspace{-1cm}
\end{wrapfigure}
 \hspace{5.5cm} Santiago, \today \\
\vspace{1cm}

\begin{center}
\textbf{\large{\textcolor{Blue}{Reporte Comparación Políticas de Testeo}}}
\end{center}

%\tableofcontents
%\setcounter{subsection}{1}

\section{Introducción}
El presente informe ilustra los resultados de {{#:AB}}{{NumP}}{{/:AB}} \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {políticas} \else {política}\fi de testeo mediante la ejecución de {{#:AB}}{{Rep}}{{/:AB}} simulaciones. La instancia consta de una población total de {{#:AB}}{{N}}{{/:AB}} individuos, quienes interactúan entre sí por un total de {{#:AB}}{{T}}{{/:AB}} días. Cada individuo pertenece a uno de los siguientes grupos definidos: {{#:AB}}{{Nomb Grupos}}{{/:AB}}. Los detalles de los parámetros de la instancia se encuentran en el Anexo~\ref{sec:Anexo}, mientras que las especificaciones de \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {las políticas de testeo utilizadas se encuentran} \else {la política de testeo utilizada se encuentra}\fi en la Tabla~\ref{tab:Comparacion} de la Sección~\ref{sec:Politicas}. Los resultados de las simulaciones incluyen gráficos y valores promedio para diferentes variables de interés referente a: número de contagiados, número de test utilizados, y número de individuos en cuarentena, en las secciones 3, 4 y 5, respectivamente.

\section{Políticas}\label{sec:Politicas}
Los parámetros usados para \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {cada una de las políticas de testeo se muestran} \else {la política de testeo se muestran}\fi en la Tabla~\ref{tab:Comparacion}.
\begin{table}[H]
\centering
\begin{tabular}{|l|{{#:AB}}{{Superc2}}{{/:AB}}    }
\hline
                    \textbf{Nombre Política}    {{#:A}}& {{Politica}} {{/:A}}                \\
                    \textbf{Tipo de Test aleatorio}    {{#:A}}& {{Tipo}}{{/:A}}                 \\
            \textbf{Tamaño de grupos Test aleatorio}    {{#:A}}& {{Tpool}}{{/:A}}                 \\
					 \textbf{Utiliza Test a sintomaticos}    {{#:A}}& {{TestS}}{{/:A}}                 \\
					 \textbf{Primer día de Test aleatorio}     {{#:A}}& {{DiasTest}}{{/:A}}                 \\
					 \textbf{Frecuencia de Test aleatorio}    {{#:A}}& {{FrecTest}}{{/:A}}                 \\
                   \textbf{Cuarentena grupo}     {{#:A}}& {{GrupoQua}}{{/:A}}                 \\ \hline
					    \end{tabular}

\caption{Tabla \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {con la comparación de los par\'ametros utilizados en las políticas}\else {con los parámetros de la política}\fi. }
\label{tab:Comparacion}
\end{table}
\section{Infectados}
\subsection{Infectados Activos}
Los \textit{infectados activos} son aquellos individuos que están con el virus un día en particular, considerando tanto casos que hayan o no sido detectados.\ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 { La Figura~\ref{fig:InfectadosActivos} muestra los \textit{infectados activos promedio} (de entre todas las simulaciones) para cada día con cada una de las políticas de testeo. Se puede apreciar que la política que tiene el peak más alto de infectados es la{{#:AB}}\textbf{ {{MInfectadospeakP}},} {{/:AB}} alcanzando un total de {{#:AB}}{{MInfectadospeak}}{{/:AB}} infectados en el día {{#:AB}}{{MInfectadospeakD}}{{/:AB}}. Por otro lado la polítca con el menor peak es la de{{#:AB}}\textbf{ {{mInfectadospeakP}},} {{/:AB}} con un peak de {{#:AB}}{{mInfectadospeak}}{{/:AB}} infectados activos en el día {{#:AB}}{{mInfectadospeakD}}{{/:AB}}.} \else{ La Figura~\ref{fig:InfectadosActivos} muestra los \textit{infectados activos promedio} (de entre todas las simulaciones) para cada día. Se puede apreciar que la política{{#:AB}}\textbf{ {{MInfectadospeakP}},} {{/:AB}} alcanza un total de {{#:AB}}{{MInfectadospeak}}{{/:AB}} infectados el día {{#:AB}}{{MInfectadospeakD}}{{/:AB}}.} \fi
Los valores anteriores se muestran en la Tabla~\ref{tab:ResumenInfectadosActivos}.
\begin{table}[H]
\centering
\begin{tabular}{|c|c|c|c|}
\hline
\textbf{Política}&\textbf{Día de Peak}&\textbf{Peak de infectados}&\textbf{Promedio por día}\\ \hline
{{#:A}}\textbf{ {{Politica}} } & {{dPeakInf}} & {{PeakInf}} & {{Infectados2T}}\\  {{/:A}}  \hline
\end{tabular}
\caption{Tabla resumen \textit{infectados activos} para la \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {comparación para las {{#:AB}}{{NumP}}{{/:AB}} políticas}\else {política}\fi, con un total de  {{#:AB}}{{N}}{{/:AB}} personas. }
\label{tab:ResumenInfectadosActivos}
\end{table}

\begin{figure}[H]
\centering
	\includegraphics[scale=0.45]{Guardar/{{#:AB}}{{Caso}}{{/:AB}}/Infectados1.pdf}
  \caption{Gráfico de los \textit{infectados activos} \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {comparando las   {{#:AB}}{{NumP}}{{/:AB}} políticas}\else {para la política}\fi, con un total de  {{#:AB}}{{N}}{{/:AB}} personas.  }
  \label{fig:InfectadosActivos}
\end{figure}
\subsection{Infectados Activos No Detectados Trabajando}
Los \textit{infectados activos no detectados trabajando} corresponden a aquellos individuos que además de estar con el virus activo, no han sido detectados y están en turno de trabajo. \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 { La Figura~\ref{fig:InfectadosTrabajando} muestra los \textit{infectados activos trabajando promedio} (de entre todas las simulaciones) para cada día con cada una de las políticas de testeo. Se puede apreciar que la política que tiene el peak más alto de infectados es la{{#:AB}}\textbf{ {{MInfectadosTpeakP}},} {{/:AB}} alcanzando un total de {{#:AB}}{{MInfectadosTpeak}}{{/:AB}} infectados en el día {{#:AB}}{{MInfectadosTpeakD}}{{/:AB}}. Por otro lado la política con el menor peak es la de{{#:AB}}\textbf{ {{mInfectadosTpeakP}},} {{/:AB}} con un peak de {{#:AB}}{{mInfectadosTpeak}}{{/:AB}} infectados activos en el día {{#:AB}}{{mInfectadosTpeakD}}{{/:AB}}.} \else{La Figura~\ref{fig:InfectadosTrabajando} muestra los \textit{infectados activos trabajando promedio} (de entre todas las simulaciones) para cada día. Se puede apreciar que la política{{#:AB}}\textbf{ {{MInfectadosTpeakP}},} {{/:AB}} alcanzando un total de {{#:AB}}{{MInfectadosTpeak}}{{/:AB}} infectados en el día {{#:AB}}{{MInfectadosTpeakD}}{{/:AB}}.}\fi
Los valores anteriores se muestran en la Tabla~\ref{tab:ResumenInfectadosTrabajando}.
\begin{table}[H]
\centering
\begin{tabular}{|c|c|c|c|}
\hline
\textbf{Política}&\textbf{Día de Peak}&\textbf{Peak de infectados}&\textbf{Promedio por día}\\ \hline
{{#:A}}\textbf{ {{Politica}} } & {{dPeakInfT}} & {{PeakInfT}} & {{Infectados3T}}\\ {{/:A}}  \hline
\end{tabular}
\caption{Tabla resumen \textit{infectados trabajando} para la \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {comparación para las {{#:AB}}{{NumP}}{{/:AB}} políticas}\else {política}\fi, con un total de  {{#:AB}}{{N}}{{/:AB}} personas. }
\label{tab:ResumenInfectadosTrabajando}
\end{table}
\begin{figure}[H]
\centering
    \includegraphics[scale=0.45]{Guardar/{{#:AB}}{{Caso}}{{/:AB}}/InfectadosTrabajando1.pdf}
	\caption{Gráfico de los \textit{infectados trabajando} \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {comparando las   {{#:AB}}{{NumP}}{{/:AB}} políticas}\else {para la política}\fi, con un total de  {{#:AB}}{{N}}{{/:AB}} personas.  }
	 \label{fig:InfectadosTrabajando}
\end{figure}
\subsection{Infectados Acumulado}
Los \textit{infectados acumulado} corresponde a todos los infectados, ya sea que hayan o no detectado. \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 { La Figura~\ref{fig:InfectadosAcumulados} muestra para cada una de las políticas, la cantidad de \textit{infectados acumulados promedio} (sobre las simulaciones) en cada uno de los días. Se puede apreciar que la política con más infectados acumulados al cabo de la simulación es la política{{#:AB}}\textbf{ {{MInfectadosP}},} {{/:AB}} donde se infectan {{#:AB}}{{MInfectados}}{{/:AB}}. Por otro lado, la política con menos infectados es la{{#:AB}}\textbf{ {{mInfectadosP}},} {{/:AB}} donde se infectan {{#:AB}}{{mInfectados}}{{/:AB}}.} \else{La Figura~\ref{fig:InfectadosAcumulados} muestra la cantidad de \textit{infectados acumulados promedio} (sobre las simulaciones) en cada uno de los días. Se puede apreciar que en la política{{#:AB}}\textbf{ {{MInfectadosP}},} {{/:AB}} se infectan {{#:AB}}{{MInfectados}}{{/:AB}}.} \fi
Los valores anteriores se muestran en la Tabla~\ref{tab:ResumenInfectadosAcumulados}.
\begin{table}[H]
\centering
\begin{tabular}{|c|c|c|}
\hline
\textbf{Política}&\textbf{Día de mayor contagio}&\textbf{Total}\\ \hline
{{#:A}}\textbf{ {{Politica}} } & {{dmaxInfectados}} & {{Infectados}}\\ {{/:A}}  \hline
\end{tabular}
\caption{Tabla resumen \textit{infectados acumulados} para la \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {comparación para las {{#:AB}}{{NumP}}{{/:AB}} políticas}\else {política}\fi, con un total de  {{#:AB}}{{N}}{{/:AB}} personas. }
\label{tab:ResumenInfectadosAcumulados}
\end{table}
\begin{figure}[H]
\centering
    \includegraphics[scale=0.45]{Guardar/{{#:AB}}{{Caso}}{{/:AB}}/InfectadosAcumulado1.pdf}
	\caption{Gráfico de los \textit{infectados acumulados} \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {comparando las   {{#:AB}}{{NumP}}{{/:AB}} políticas}\else {para la política}\fi, con un total de  {{#:AB}}{{N}}{{/:AB}} personas.  }
	\label{fig:InfectadosAcumulados}
\end{figure}
\section{Test}
\subsection{Test Totales}
\ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {A continuación se detalla la cantidad de tests utilizados por cada una de las políticas estudiadas.{{#:AB}}\textbf{ {{MTestTotalP}} } {{/:AB}}es la política que utiliza un mayor número tests, utilizando {{#:AB}}{{MTestTotal}}{{/:AB}} tests. Por el contrario la política{{#:AB}}\textbf{ {{mTestTotalP}} } {{/:AB}}es la que utiliza menos tests ocupando {{#:AB}}{{mTestTotal}}{{/:AB}} tests. En la Figura~\ref{fig:Tests} se puede observar número de \textit{tests utilizados} para cada día.} \else { A continuación en la Figura~\ref{fig:Tests} se detalla la cantidad de \textit{tests diarios utilizados} por la política{{#:AB}}\textbf{ {{MTestTotalP}},} {{/:AB}} ocupó {{#:AB}}{{MTestTotal}}{{/:AB}} tests.} \fi
\begin{figure}[H]
\centering
    \includegraphics[scale=0.45]{Guardar/{{#:AB}}{{Caso}}{{/:AB}}/Tests1.pdf}
	\caption{Gráfico de los \textit{test diarios utilizados} \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {al comparar las   {{#:AB}}{{NumP}}{{/:AB}} políticas}\else {por la política}\fi, con un total de  {{#:AB}}{{N}}{{/:AB}} personas.  }
	\label{fig:Tests}
\end{figure}
En la Figura~\ref{fig:Testshistograma} se observa el número de \textit{tests utilizados} por \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {las  {{#:AB}}{{NumP}}{{/:AB}} políticas}\else {por la política}\fi.
\begin{figure}[H]
\centering
    \includegraphics[scale=0.45]{Guardar/{{#:AB}}{{Caso}}{{/:AB}}/TestsHistograma1.pdf}
	\caption{Histograma de los \textit{test totales utilizados} \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {al comparar las  {{#:AB}}{{NumP}}{{/:AB}} políticas}\else {por la política}\fi, con un total de  {{#:AB}}{{N}}{{/:AB}} personas.  }
	\label{fig:Testshistograma}
\end{figure}
Para un mayor detalle de los tests se puede consultar el Anexo, la subsección~\ref{subsec:TestsAleatorios} y  la subsección~\ref{subsec:TestsSintomaticos}
\section{Cuarentena}
\subsection{Personas en cuarentena}
\ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 { La Figura~\ref{fig:Cuarentena} muestra las \textit{personas en cuarentena en promedio} (de entre todas las simulaciones) para cada día con cada una de las políticas de testeo. Se puede apreciar que la política que tiene el peak más alto es la{{#:AB}}\textbf{ {{MCuarentenapeakP}},} {{/:AB}} alcanzando un total de {{#:AB}}{{MCuarentenapeak}}{{/:AB}} personas en cuarentena el día {{#:AB}}{{MCuarentenapeakD}}{{/:AB}}. Por otro lado la política con el menor peak es la de{{#:AB}}\textbf{ {{mCuarentenapeakP}},} {{/:AB}} con un peak de {{#:AB}}{{mCuarentenapeak}}{{/:AB}} durante el día {{#:AB}}{{mCuarentenapeakD}}{{/:AB}}.} \else{La Figura~\ref{fig:Cuarentena} muestra las \textit{personas en cuarentena en promedio} (de entre todas las simulaciones) para cada día. Se puede apreciar que la política{{#:AB}}\textbf{ {{MCuarentenapeakP}} } {{/:AB}} alcanza un total de {{#:AB}}{{MCuarentenapeak}}{{/:AB}} personas en cuarentena el día {{#:AB}}{{MCuarentenapeakD}}{{/:AB}}.} \fi
\begin{table}[H]
\centering
\begin{tabular}{|c|c|c|c|}
\hline
\textbf{Política}&\textbf{Día de peak}&\textbf{Peak de cuarentena}&\textbf{Promedio por día}\\ \hline
{{#:A}} \textbf{ {{Politica}} } & {{dPeakQua}} & {{PeakQua}} & {{QuaT}}\\  {{/:A}}  \hline
\end{tabular}
\caption{Tabla resumen \textit{personas en cuarentena} para la \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {comparación para las {{#:AB}}{{NumP}}{{/:AB}} políticas}\else {política}\fi, con un total de  {{#:AB}}{{N}}{{/:AB}} personas. }
\label{tab:ResumenPersonasencuarentena}
\end{table}
En la siguiente Figura~\ref{fig:Cuarentena} se puede observar cuantas personas están cuarentena cada día.\\
\begin{figure}[H]
\centering
    \includegraphics[scale=0.45]{Guardar/{{#:AB}}{{Caso}}{{/:AB}}/Cuarentena1.pdf}
	\caption{Gráfico de las \textit{personas en cuarentena} \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {al comparar las   {{#:AB}}{{NumP}}{{/:AB}} políticas}\else {para la política}\fi, con un total de  {{#:AB}}{{N}}{{/:AB}} personas.  }
	\label{fig:Cuarentena}
\end{figure}
\newpage
\section{Anexo}\label{sec:Anexo}
\subsection{Parámetros}
\textbf{\large{Individuales}}\\
\begin{table}[H]
\centering
\begin{tabular}{|c|c|c|c|c|}
\hline
\textbf{Copias} & \textbf{Nombre Grupo} & \textbf{Nombre Horario} & \textbf{Prob Cont Int.} & \textbf{Prrob Cont Ext.}  \\ \hline
{{#:A4}} {{Copias}} & {{Nombres Grupos}} & {{Nombres Horario}} & {{Prob Cont Int}} & {{Prob Cont Ext}} \\ \hline {{/:A4}}
\end{tabular}
\caption{Tabla resumen \textit{parámetros individuales} de \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {las {{#:AB}}{{NumP}}{{/:AB}} políticas}\else {la política}\fi, con un total de  {{#:AB}}{{N}}{{/:AB}} personas. }
\label{tab:ParametrosIndividuales}
\end{table}
\textbf{\large{Grupales}}\\
\begin{table}[H]
\centering
\begin{tabular}{|c| {{#:AB}}{{Superc}}{{/:AB}} }
\hline
\textbf{Grupos}{{#:A2}}&\textbf{ {{Nombres Grupos}} }{{/:A2}}\\ \hline
{{#:A2}} \textbf{ {{Nombres Grupos}} } {{#Matriz rel}} & {{.}} {{/Matriz rel}}\\  {{/:A2}}\hline
\end{tabular}
\caption{Tabla con la \textit{matriz de relaciones} entre grupos. }
\label{tab:MatrizRelacionesGrupo}
\end{table}
\begin{table}[H]
\centering
\begin{tabular}{ |c|{{#:AB}}{{Superc2}}{{/:AB}}    }
\hline
\textbf{Grupos}  {{#:A}}&\textbf{ {{Politica}}  }{{/:A}} \\ \hline
{{#:A2}} \textbf{ {{Nombres Grupos}} } {{#Test Grupos}} &{{.}} {{/Test Grupos}} \\ {{/:A2}}\hline
\end{tabular}
\caption{Porcentaje a testear por  grupos. }
\label{tab:Porcentaje a testerar por grupo}
\end{table}
\textbf{\large{Horarios}}\\
{{#:A5}} \begin{table}[H]
\centering
\begin{tabular}{ |c|c|c|c|c|c|c|c|    }
\hline
\textbf{Turno} & \textbf{L} & \textbf{M} & \textbf{M} & \textbf{J} & \textbf{V} & \textbf{S} & \textbf{D}  \\ \hline
\textbf{Día}  {{#HorarioD}}  & $\;\;{{.}}\;\;$ {{/HorarioD}} \\ \hline
\textbf{Noche} {{#HorarioN}}  & $\;\;{{.}}\;\;$ {{/HorarioN}} \\ \hline
\end{tabular}
\caption{\textit{Horario} {{NHorario}}. }
\label{tab:Horario{{Conteo}} }
\end{table} {{/:A5}}
\subsection{Gráficos Infectados Activos por grupo}
{{#:A2}}\begin{figure}[H]
\centering
  \includegraphics[scale=0.35]{Guardar/{{#:AB}}{{Caso}}{{/:AB}}/Infectados{{Ind}}.pdf}
  \caption{\textit{Infectados activos} {{Nombres Grupos}} }
  \label{fig:InfectadosActivos{{Ind}} }
\end{figure} {{/:A2}}
\subsection{Gráficos Infectados Trabajando por grupo}
{{#:A2}} \begin{figure}[H]
\centering
  \includegraphics[scale=0.35]{Guardar/{{#:AB}}{{Caso}}{{/:AB}}/InfectadosTrabajando{{Ind}}.pdf}
  \caption{\textit{Infectados trabajando} {{Nombres Grupos}}. }
  \label{fig:InfectadosTrabajando{{Ind}} }
\end{figure} {{/:A2}}
\subsection{Gráficos Infectados Acumulados por grupo}
{{#:A2}} \begin{figure}[H]
\centering
  \includegraphics[scale=0.35]{Guardar/{{#:AB}}{{Caso}}{{/:AB}}/InfectadosAcumulado{{Ind}}.pdf}
  \caption{\textit{Infectados Acumulados} {{Nombres Grupos}}. }
  \label{fig:InfectadosAcumulados{{Ind}} }
\end{figure} {{/:A2}}
\subsection{Tests Aleatorios}\label{subsec:TestsAleatorios}
\ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {A continuación se detalla la cantidad de tests utilizados por cada una de las políticas estudiadas.{{#:AB}}\textbf{ {{MTestP}} } {{/:AB}} es la política que utiliza un mayor número tests, utilizando {{#:AB}}{{MTest}}{{/:AB}} tests. Por el contrario la política{{#:AB}}\textbf{ {{mTestP}} } {{/:AB}} es la que utiliza menos tests ocupando {{#:AB}}{{mTest}}{{/:AB}} tests. En la Figura~\ref{fig:TestsAleatorios} se puede observar número de \textit{tests utilizados} para cada día.} \else { A continuación en la Figura~\ref{fig:TestsAleatorios} se detalla la cantidad de \textit{tests diarios utilizados} por la política{{#:AB}}\textbf{ {{MTestP}},} {{/:AB}} ocupó {{#:AB}}{{MTest}}{{/:AB}} tests.} \fi
\begin{figure}[H]
\centering
    \includegraphics[scale=0.35]{Guardar/{{#:AB}}{{Caso}}{{/:AB}}/TestsAleatorios1.pdf}
	\caption{Gráfico de los \textit{test diarios utilizados} \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {al comparar las   {{#:AB}}{{NumP}}{{/:AB}} políticas}\else{por la política}\fi, con un total de  {{#:AB}}{{N}}{{/:AB}} personas.  }
	\label{fig:TestsAleatorios}
\end{figure}
En la Figura~\ref{fig:TestsAleatorioshistograma} se observa el número de tests utilizados por \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {las  {{#:AB}}{{NumP}}{{/:AB}} políticas}\else{la política}\fi.
\begin{figure}[H]
\centering
    \includegraphics[scale=0.35]{Guardar/{{#:AB}}{{Caso}}{{/:AB}}/TestsAleatoriosHistograma1.pdf}
	\caption{Histograma de los \textit{test totales utilizados} \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {al comparar las   {{#:AB}}{{NumP}}{{/:AB}} políticas}\else{por la política}\fi, con un total de  {{#:AB}}{{N}}{{/:AB}} personas.  }
	\label{fig:TestsAleatorioshistograma}
\end{figure}
\subsection{Test Sintomáticos}\label{subsec:TestsSintomaticos}
\ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {A continuación se detalla la cantidad de tests utilizados por cada una de las políticas estudiadas.{{#:AB}}\textbf{ {{MTestP2}} } {{/:AB}} es la política que utiliza un mayor número tests, utilizando {{#:AB}}{{MTest2}}{{/:AB}} tests. Por el contrario la política{{#:AB}}\textbf{ {{mTestP2}} } {{/:AB}} es la que utiliza menos tests ocupando {{#:AB}}{{mTest2}}{{/:AB}} tests. En la Figura~\ref{fig:TestsSintomaticos} se puede observar número de \textit{tests utilizados} para cada día.} \else { A continuación en la Figura~\ref{fig:TestsSintomaticos} se detalla la cantidad de \textit{tests diarios utilizados} por la política{{#:AB}}\textbf{ {{MTestP2}},} {{/:AB}} ocupó {{#:AB}}{{MTest2}}{{/:AB}} tests.} \fi
\begin{figure}[H]
\centering
    \includegraphics[scale=0.35]{Guardar/{{#:AB}}{{Caso}}{{/:AB}}/TestsSintomaticos1.pdf}
	\caption{Gráfico de los \textit{test diarios} utilizados \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {al comparar las   {{#:AB}}{{NumP}}{{/:AB}} políticas}\else{por la política}\fi, con un total de  {{#:AB}}{{N}}{{/:AB}} personas.  }
	\label{fig:TestsSintomaticos}
\end{figure}
En la Figura~\ref{fig:TestsSintomaticoshistograma} se observa el número de tests utilizados por \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {las  {{#:AB}}{{NumP}}{{/:AB}} políticas}\else{la política}\fi.
\begin{figure}[H]
\centering
    \includegraphics[scale=0.35]{Guardar/{{#:AB}}{{Caso}}{{/:AB}}/TestsSintomaticosHistograma1.pdf}
	\caption{Histograma de los \textit{test totales} utilizados \ifnum {{#:AB}}{{NumP}}{{/:AB}}>1 {al comparar las {{#:AB}}{{NumP}}{{/:AB}} políticas}\else{por la política}\fi, con un total de  {{#:AB}}{{N}}{{/:AB}} personas.  }
	\label{fig:TestsSintomaticoshistograma}
\end{figure}
\subsection{Gráficos Cuarentena por grupos}
{{#:A2}} \begin{figure}[H]
\centering
  \includegraphics[scale=0.35]{Guardar/{{#:AB}}{{Caso}}{{/:AB}}/Cuarentena{{Ind}}.pdf}
  \caption{\textit{Personas en cuarentena} {{Nombres Grupos}}. }
  \label{fig:Cuarentena{{Ind}} }
\end{figure} {{/:A2}}


\end{document}
"""

  rendered = render(tmpl,A=A,A2=A2,A4=A4,AB=AB[1],A5=A5)

  filename = string("Guardar/"*Caso*"/Reporte.tex")
  open(filename, "w") do file
    write(file, rendered)
  end
  run(`pdflatex $filename`)
 run(`pdflatex $filename`)
 mv("Reporte.pdfsync","Guardar/"*Caso*"/Reporte.pdfsync",force=true)
 mv("Reporte.out","Guardar/"*Caso*"/Reporte.out",force=true)
 mv("Reporte.pdf","Guardar/"*Caso*"/Reporte.pdf",force=true)
 mv("Reporte.aux","Guardar/"*Caso*"/Reporte.aux",force=true)
 mv("Reporte.log","Guardar/"*Caso*"/Reporte.log",force=true)
