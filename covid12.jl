using Distributions, Plots, LinearAlgebra, Random, StatsPlots,JSON

#######################################################################################
#Recibe un grupo y entrega un vector con los indices de los grupos a los que pertenecen
function vecgroup(group)
S=Int(length(group[1,:]))
A=zeros(N,S)
v=zeros(N)
for i in 1:N
	for j in 1:S
		if group[i,j]==1
			v[i]=j
		end
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
	for i in 1:L
	p[i]=(1-prod(-M_cont[:,i].+1))
	end
	#de no contagiarse multiplicado por un factor que pondera las horas, este desde
	return p
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
#VInf1: Numero de abuelitos infectados en promedio, por dia
#VInf2: Numero de funcionarios infectados en promedio, por dia
function simulation(N,Group,T,Tra,Mrel, f,G,peak,t_peak,p_int,p_ext,γ,p_false_positive,R,Politica,random,quienes,cuarentena,z,dias_atras,scalar_asint,test_sym,distribuir,frec)

#N Tamaño grupo de trabajadores
#Group una matriz con columnas indicatrices de tamaño N con los integrantes de cada grupo.
#Tra; horario de los trabajadores
#Mrel: Matriz de presencia de todos los trabajadores, tiene una tercer componente temporal
#f Frecuencia testeo (puede ser vector de largo T, en caso de testear random de tamaño G o matriz de NxT con los individuos a testear en sus fechas)
# G Solo aplica para el caso random el número de test individuales o pool (por ejemplo G=4 son 4 pool testing y todos los sub test que se generen)
# T simulation horizon
#peak: vector de probabilidad contagio peak
#γ:  prop of asymtomatic pop
#p_false_positive probability of a false positive
#R Replications
#Politica dice si es "Pool", "Individual" o "No testear".
#Tra: Matriz (2*N,T) de turnos de trabajo
#random: Si es random o estan fijos a quienes se testeara (no los dias esos estan en f, solo los funcionarios a testear)
#quienes: trabajando, no trabajando y ambos
#cuarentena: Si se mandan a todo el grupo en cuarentena o solo al infectado (ahora tambien mixto, solo al grupo cuando uno presenta sintomas)
#z: dias de cuarntena al grupo, no se usa en caso de solo mandar al positivo
#t_peak: tiempo en que la exponencial alcanza la probabilidad peak
#dias_atras: cuantos dias hacia atras se mira con quienes has tenido contacto estrecho
#test_sym: Si se quiere testear a la gente apenas le aparescan los sintomas
#distribuir: Si se distribuye o no los test (solo aplica para random)
#frec: en cuantos dias se distribuye
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
	mNFpositive=zeros(T)
maxInf=zeros(N,T)
maxQua=zeros(N,T)
maxSy=zeros(N,T)
Infect=ones(R)*1000
VInf1=zeros(N,T)
VInf2=zeros(N,T)


    # Loop over replications
    for rep = 1:R
	print(rep) #Lo uso para saber si corre
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
	NFpositive2=zeros(T) #Falsos positivos diarios
	indicatriz=zeros(N,T) #Quienes testear de cada grupo caso random
	resagados=zeros(N)
	resagados2=zeros(N)
	t_inicio=1
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
			if (random=="si")
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
				if (random=="si") & (f[t]==1)
					for i in findall(g.==1)
						indicatriz[i,t]=(rand()<=G[a,t])
					end
					if quienes=="trabajando"
					I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).* (MD[:,t]+MN[:,t].>=1).*(-qu.+1) #Gente que esta trabajando, no en cuarentena y nunca a presentado sintomas
					end
					if quienes=="no trabajando"
					I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).* (MD[:,t]+MN[:,t].<=1).*(-qu.+1) #Gente que esta trabajando, no en cuarentena y nunca a presentado sintomas
					end
					if quienes=="ambos"
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
           	for i in findall(new_inf.==1)                   # cicle over newly infected
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
		if test_sym=="si"
			for i in findall(sy.*(-qu.+1).==1)
				NTest2[t]+=1   #Un test nuevo
						if (su[i]==0)&(re[i]==0)
									p_positive = true_positive(t-tEx[i]+1, tIn[i]-tEx[i]+1, tSy[i]-tEx[i]+1, tRe[i]-tEx[i]+1,As[i],scalar_asint)
									if rand() < p_positive #Caso positivo
											Qu[i,Int(min(T,t+1)):Int(min(T,max(t+15,tRe[i]+3)))] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
											if cuarentena=="grupo"
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
		if random=="si" #Se realiza test al azar
	    		if Politica=="Individual" #test individual
					if f[t]==1 #Se testea estos dias
						if quienes=="trabajando"
							I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).* (MD[:,t]+MN[:,t].>=1).*(-qu.+1).*collect(1:N) #Gente que esta trabajando, no en cuarentena y nunca a presentado sintomas
						end
						if quienes=="no trabajando"
							I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).* (MD[:,t]+MN[:,t].<=1).*(-qu.+1).*collect(1:N) #Gente que esta trabajando, no en cuarentena y nunca a presentado sintomas
						end
						if quienes=="ambos"
							I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*(-qu.+1).*collect(1:N) #Gente que esta trabajando, no en cuarentena y nunca a presentado sintomas
						end
						J=shuffle!(I) #Mezclar
						II=setdiff(J,[0]) #Eliminar ceros
						cand=Int.(II)  #Se toman como candidatos a los posibles y el minimo entre ellos y G
						if distribuir=="si"
							resagados=cand
							t_inicio=t
						end
					end
					if distribuir=="si"
						cand=[]
						if t==1
							if quienes=="trabajando"
								I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).* (MD[:,t]+MN[:,t].>=1).*(-qu.+1).*collect(1:N) #Gente que esta trabajando, no en cuarentena y nunca a presentado sintomas
							end
							if quienes=="no trabajando"
								I=( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).* (MD[:,t]+MN[:,t].<=1).*(-qu.+1).*collect(1:N) #Gente que esta trabajando, no en cuarentena y nunca a presentado sintomas
							end
							if quienes=="ambos"
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
						vector=pool(resagados,5)
					end
					if distribuir=="si"
						resagados2=zeros(N)
						for i in 1:N
							if t-t_inicio==mod(i-1,frec)
								resagados2[i]=resagados[i]
							end
						end
						vector=pool(resagados2,5)
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
									if cuarentena=="grupo"
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
									if cuarentena=="grupo"
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
							if rand()<=G[v[i],t]
							NTest[t]+=1   #Un test nuevo
	           				     		if (((As[i]==1)&(su[i]==0)&(re[i]==1))|(su[i]==1))& (rand() < p_false_positive) #Caso falso positivo
				                    		Qu[i,min(T,t+1):min(T,t+15)] .= 1 # if someone tests (false) positive, quarantined for 2 weeks
											if cuarentena=="grupo"
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
														if cuarentena=="grupo"
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
		if random=="no" #Se realiza el esquema de test de manera similar solo que los candidatos ya estan dados por defecto y solo saca los que esten en cuarentena o ya hayan tenido sintomas en el pasado
			if sum( f[:,t].*(-qu.+1))>=1
	  	  		if Politica=="Individual"
					J=Int.(( sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*f[:,t].*(-qu.+1).*collect(1:N)) #Son los que nunca han tenido sintomas, no cuarentena y toca testear
					II=setdiff(J,[0])
					cand=Int.(II)  #Se obtiene vector de a quienes testear
	    		end
				if Politica=="Pool" #Esta medio obsoleto
					h=(sum(Sy[:,1:max((t-2),1)],dims=2).==0 ).*f[:,t].*(-qu.+1)
					vector=pool(h,5)
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
									if cuarentena=="grupo"
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
									if cuarentena=="grupo"
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
										if cuarentena=="grupo"
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
													if cuarentena=="grupo"
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
	VInf1=VInf1/R
	VInf2=VInf2/R


    return mQua,mQua2,mInf,mInf2,mInf3, mNFpositive,  mNTest, mSy, maxInf, maxQua, maxSy, Infect, VInf1, VInf2
end
################################################################################################################################################
function leer(dict12)

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
	Horario=zeros(2*Integrantes,T)
	Mrel=zeros(S,S)
	cont=0
	for s in 1:S
	for i in 1:length(dict12["Grupos"][s])
		Group[cont+1:dict12["Grupos"][s][i]["Numero de copias"]+cont,s]=ones(dict12["Grupos"][s][i]["Numero de copias"])
		p_int[cont+1:dict12["Grupos"][s][i]["Numero de copias"]+cont]=ones(dict12["Grupos"][s][i]["Numero de copias"]).*dict12["Grupos"][s][i]["Probabilidad de contagio int"]
		p_ext[cont+1:dict12["Grupos"][s][i]["Numero de copias"]+cont]=ones(dict12["Grupos"][s][i]["Numero de copias"]).*dict12["Grupos"][s][i]["Probabilidad de contagio ext"]
	for t in 1:T
		Horario[2*cont+1:2*(dict12["Grupos"][s][i]["Numero de copias"]+cont),t]=repeat(dict12["Grupos"][s][i]["Horario"][t],dict12["Grupos"][s][i]["Numero de copias"])
	end
		cont+=dict12["Grupos"][s][i]["Numero de copias"]
	end
	Mrel[s,:]=dict12["Matriz grupos"][s]
	end
	Politica=dict12["Politica de testeo"]
	TestSyn=dict12["Testear sintomaticos"]

	if Politica!="No testear"
		Testrand=dict12["testeo random"]
		Testnorand=dict12["testeo no random"]

		if length(Testrand)>=1
			test="si"
			Testgrupo=zeros(S,T)
			for t in 1:T
				Testgrupo[:,t]=Testrand["Porcentaje a testear por grupo"][t]
			end
			Variar=Testrand["Variar"]
			if Variar==0
				quienes="ambos"
				distribuir="no"
			end
			if Variar==1
				quienes="trabajando"
				distribuir="no"
			end
			if Variar==2
				quienes="no trabajando"
				distribuir="no"
			end
			if Variar==3
				quienes="ambos"
				distribuir="si"
			end

			Diastest=zeros(T)
			for t in 1:T
				Diastest[t]=Testrand["Dias de testeo"][t]
			end
		end
		if length(Testnorand)>=1
			test="no"
			Diastest=zeros(Integrantes,T)
			for t in 1:T
				Diastest[:,t]=Testnorand["Dias de testeo"][t]
			end
			quienes=[]
			distribuir=[]
			Testgrupo=[]
		end

	else
		test="no"
		Testrand=[]
		Testnorand=[]
		Testgrupo=[]
		Variar=[]
		Diastest=zeros(Integrantes,T)
		quienes=[]
		distribuir=[]

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
	graficogrupo=dict12["Grupos grafico"]
	nombregrafico=dict12["Nombres grupos"]



return N, Group,T,Horario,Mrel, Diastest,Testgrupo,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,test,quienes,Cuaren,Diascuarentena,Diasatras,TestSyn,distribuir,graficogrupo,nombregrafico


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
		if a==1
			grafico=plot(1:T,arr,label=NombresP[a],lw=3)
		else
			plot!(1:T,arr,label=NombresP[a],lw=3)
		end
	end
	title!(titulo*string(nombreg)*" N= "*string(sum(Indices)) )
	xlabel!(x)
	ylabel!(y)
	return grafico
end




#######################################################################################################################




	###############################################################################
	#Se realizan las 4 politicas
politicas=4 #Numero de politicas
NombresP=["No testear","Simple, Todos los dias","Comparación","Simple, Cada 10 dias"]#Nombres politicas
t_peak=ones(120)*24                 #Tiempo peak de contagio
scalar_asint=1          #Escalar para los asintomaticos
peak=[(ones(100)*0.01)' (ones(20)*0.2)']' #Prob peak de infeccion




dict12 = Dict()


open("parametros1.json", "r") do f

    global dict12
    dicttxt1 = read(f,String)  # file information to string
    dict12=JSON.parse(dicttxt1)  # parse and transform data
end

N,Group,T,Horario,Mrel, Diastest,Testgrupo,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,test,quienes,Cuaren,Diascuarentena,Diasatras,TestSyn,distribuir,graficogrupo,nombregrafico=leer(dict12)
base_mQua,base_mQua2, base_mInf,base_mInf2,base_mInf3, base_mNFp, base_T, base_mSy,base_maxInf,base_maxQua, base_maxSy, base_Infect, base_VInf1, base_VInf2= simulation(N,Group,T,Horario,Mrel, Diastest,Testgrupo,peak,t_peak,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,test,quienes,Cuaren,Diascuarentena,Diasatras,scalar_asint,TestSyn,distribuir,[])



open("parametros2.json", "r") do f

    global dict12
    dicttxt1 = read(f,String)  # file information to string
    dict12=JSON.parse(dicttxt1)  # parse and transform data
end

N,Group,T,Horario,Mrel, Diastest,Testgrupo,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,test,quienes,Cuaren,Diascuarentena,Diasatras,TestSyn,distribuir,graficogrupo,nombregrafico=leer(dict12)
ideal_mQua,ideal_mQua2, ideal_mInf,ideal_mInf2,ideal_mInf3, ideal_mNFp, ideal_T, ideal_mSy,ideal_maxInf,ideal_maxQua, ideal_maxSy, ideal_Infect, ideal_VInf1, ideal_VInf2= simulation(N,Group,T,Horario,Mrel, Diastest,Testgrupo,peak,t_peak,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,test,quienes,Cuaren,Diascuarentena,Diasatras,scalar_asint,TestSyn,distribuir,[])


open("parametros3.json", "r") do f
    global dict12
    dicttxt1 = read(f,String)  # file information to string
    dict12=JSON.parse(dicttxt1)  # parse and transform data
end

N,Group,T,Horario,Mrel, Diastest,Testgrupo,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,test,quienes,Cuaren,Diascuarentena,Diasatras,TestSyn,distribuir,graficogrupo,nombregrafico=leer(dict12)
f1_mQua,f1_mQua2, f1_mInf,f1_mInf2,f1_mInf3, f1_mNFp, f1_T, f1_mSy,f1_maxInf,f1_maxQua, f1_maxSy,f1_Infect, f1_VInf1, f1_VInf2 = simulation(N,Group,T,Horario,Mrel, Diastest,Testgrupo,peak,t_peak,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,test,quienes,Cuaren,Diascuarentena,Diasatras,scalar_asint,TestSyn,distribuir,[])


open("parametros4.json", "r") do f

    global dict12
    dicttxt1 = read(f,String)  # file information to string
    dict12=JSON.parse(dicttxt1)  # parse and transform data
end

N,Group,T,Horario,Mrel, Diastest,Testgrupo,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,test,quienes,Cuaren,Diascuarentena,Diasatras,TestSyn,distribuir,graficogrupo,nombregrafico=leer(dict12)
fr_mQua,fr_mQua2, fr_mInf,fr_mInf2,fr_mInf3, fr_mNFp, fr_T, fr_mSy,fr_maxInf,fr_maxQua, fr_maxSy, fr_Infect, fr_VInf1, fr_VInf2 = simulation(N,Group,T,Horario,Mrel, Diastest,Testgrupo,peak,t_peak,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,test,quienes,Cuaren,Diascuarentena,Diasatras,scalar_asint,TestSyn,distribuir,[])




###############################################################################
###############Graficos########################################################

	arrays=[[base_mQua],[ideal_mQua],[f1_mQua],[fr_mQua]]
	titulo="Cuarentena "
	x="Dias"
	y="Personas"
for i in 1:length(graficogrupo)
	ggrafico=graficogrupo[i]
	nombreg=nombregrafico[i]
	g=graficar(arrays,Group,ggrafico,nombreg,NombresP,titulo,x,y,T)
	plot(g)
	savefig(titulo*string(i)*".pdf")
	g=[]
end


arrays=[[base_mInf],[ideal_mInf],[f1_mInf],[fr_mInf]]
titulo="Infectados Trabajando "
x="Dias"
y="Personas"
for i in 1:length(graficogrupo)
ggrafico=graficogrupo[i]
nombreg=nombregrafico[i]
g=graficar(arrays,Group,ggrafico,nombreg,NombresP,titulo,x,y,T)
plot(g)
savefig(titulo*string(i)*".pdf")
g=[]
end


arrays=[[base_mInf3],[ideal_mInf3],[f1_mInf3],[fr_mInf3]]
titulo="Infectados "
x="Dias"
y="Personas"
for i in 1:length(graficogrupo)
ggrafico=graficogrupo[i]
nombreg=nombregrafico[i]
g=graficar(arrays,Group,ggrafico,nombreg,NombresP,titulo,x,y,T)
plot(g)
savefig(titulo*string(i)*".pdf")
g=[]
end

arrays=[[base_T],[ideal_T],[f1_T],[fr_T]]
titulo="Test "
x="Dias"
y="Personas"
for i in 1:length(graficogrupo)
ggrafico=graficogrupo[i]
nombreg=nombregrafico[i]
g=graficar(arrays,Group,ggrafico,nombreg,NombresP,titulo,x,y,T)
plot(g)
savefig(titulo*string(i)*".pdf")
g=[]
end
