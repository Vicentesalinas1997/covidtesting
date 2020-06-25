using JSON

############################################################

N=100   #Cantidad de abuelitos
T=100    #Cantidad de tiempo
m=20      #Cantidad de funcionarios
trab=ones(2*N,T)    #Dias de trabajo abuelitos
t_peak=ones(N+m)*24                 #Tiempo peak de contagio
z=14			#Dias de cuarentena preventiva
dias_atras=0          #Dias a mirar hacia atras para cuarentena grupo
scalar_asint=1          #Escalar para los asintomaticos
p_int = [(ones(N)*0.01*(1/30)*0.2)' zeros(m)']'   #probabilidad interna de contagio
p_ext = [zeros(N)' (ones(m).*0.01)']'     #probabilidad externa de contagio
γ = 0.3 #(antes 0.25)           # prop of asymtomatic pop
p_false_positive = 0.01 # probability of a false positive
R = 400 # Replications
peak=[(ones(N)*0.01)' (ones(m)*0.2)']' #Prob peak de infeccion
S2=Int(m/2)     #Grupos de tipo 1
S1=Int(N/S2)   #Grupos de tipo 2
S=S1+S2      #Grupos
G=ones(S,T)           #Caso no random tamaño de grupos/podria implementarse como una cantidad diaria limite de test
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
	################################################
	#Creacion de un grupo
	Group=zeros(N+m,S)
	for i in 1:S1
		Group[1:N,i]=[zeros((i-1)*Int(N/S1))' ones(Int(N/S1))'  zeros((S1-i)*Int(N/S1))'  ]'
	end
	for j in 1:S2
		Group[N+j,j]=1
		Group[N+j+Int(m/2),j]=1
	end
	##############################################
	#Horarios e intereracciones funcionarios
	Ftrab2=zeros(2*m,T) #Funcionarios no trabajan
	Ftrab2=zeros(2*m,T) #Funcionariso se relacionan con10 abuelitos dia por medio y no entre ellos
	for t in 1:T
		if mod(t,2)==1
			Ftrab2[1:Int(m/2),t]=ones(Int(m/2))
		else
			Ftrab2[Int(m/2)+1:m,t]=ones(Int(m/2))
		end
	end
	Mrel=zeros(S,S)
	Mrel2=zeros(S,S)
	for i in 1:S
		for j in 1:S
			if i==j
				Mrel2[i,j]=1
			else
				if ((i<=S1)&(j<=S1))|((i>S1)&(j>S1))
					Mrel2[i,j]=1
				else
					if i==j-S1
						Mrel2[i,j]=1
					end
				end

			end
		end
	end
	Mrel2=ones(S,S)
	###############################################################################
	#Se realizan las 4 politicas
politicas=4 #Numero de politicas
NombresP=["No testear","Simple, Todos los dias","Comparación","Simple, Cada 10 dias"]#Nombres politicas

#################################################################################################################
#print(N,m,Group,T,trab,Ftrab2,Mrel2, Frecuencia0,G,peak,t_peak,p_int,p_ext,γ,p_false_positive,R,"no","no",[],"grupo",z,dias_atras,scalar_asint,"no","no",[])




### Write data ####

###################
#"testeo random"=>Dict("Porcentaje a testear por grupo"=>ones(S,T)*0.5,"Dias de testeo"=>ones(T) ,"Variar"=>2),
#"testeo no random"=>Dict([])


# dictionary to write
dict1 = Dict("Tiempo"=>T,"Grupos"=>[[[Dict("Numero de copias"=>10,"Probabilidad de contagio int"=>(0.01*0.2*1/30),"Probabilidad de contagio ext"=>0,
"Horario"=>trab[1:2,:])] for i=1:10]' [[Dict("Numero de copias"=>1,"Probabilidad de contagio int"=>0,"Probabilidad de contagio ext"=>0.01,
"Horario"=>Ftrab2[2*j-1:2*j,:]),Dict("Numero de copias"=>1,"Probabilidad de contagio int"=>0,"Probabilidad de contagio ext"=>0.01,
"Horario"=>Ftrab2[20+2*j-1:20+2*j,:])] for j=1:10]']',"Matriz grupos"=>Mrel2,"Politica de testeo"=>"No testear","Testear sintomaticos"=>"no",
"testeo random"=>Dict([]),
"testeo no random"=>Dict([]), "Repeticiones"=>R, "Porcentaje asintomaticos"=>0.3,
"Probabilidad falso positivo"=>0.01, "Cuarentena"=>Dict("Tipo"=>"grupo","Dias atras"=>3,"Dias cuarentena"=>14) )


dict2 = Dict("Tiempo"=>T,"Grupos"=>[[[Dict("Numero de copias"=>10,"Probabilidad de contagio int"=>(0.01*0.2*(1/30)),"Probabilidad de contagio ext"=>0,
"Horario"=>trab[1:2,:])] for i=1:10]' [[Dict("Numero de copias"=>1,"Probabilidad de contagio int"=>0,"Probabilidad de contagio ext"=>0.01,
"Horario"=>Ftrab2[2*j-1:2*j,:]),Dict("Numero de copias"=>1,"Probabilidad de contagio int"=>0,"Probabilidad de contagio ext"=>0.01,
"Horario"=>Ftrab2[20+2*j-1:20+2*j,:])] for j=1:10]']',"Matriz grupos"=>Mrel2,"Politica de testeo"=>"Pool","Testear sintomaticos"=>"no",
"testeo random"=>Dict([]),
"testeo no random"=>Dict("Dias de testeo"=>Frecuencia12), "Repeticiones"=>R, "Porcentaje asintomaticos"=>0.3,
"Probabilidad falso positivo"=>0.01, "Cuarentena"=>Dict("Tipo"=>"solo","Dias atras"=>0,"Dias cuarentena"=>0) )



dict3 = Dict("Tiempo"=>T,"Grupos"=>[[[Dict("Numero de copias"=>10,"Probabilidad de contagio int"=>(0.01*0.2*(1/30)),"Probabilidad de contagio ext"=>0,
"Horario"=>trab[1:2,:])] for i=1:10]' [[Dict("Numero de copias"=>1,"Probabilidad de contagio int"=>0,"Probabilidad de contagio ext"=>0.01,
"Horario"=>Ftrab2[2*j-1:2*j,:]),Dict("Numero de copias"=>1,"Probabilidad de contagio int"=>0,"Probabilidad de contagio ext"=>0.01,
"Horario"=>Ftrab2[20+2*j-1:20+2*j,:])] for j=1:10]']',"Matriz grupos"=>Mrel2,"Politica de testeo"=>"Pool","Testear sintomaticos"=>"no",
"testeo random"=>Dict([]),
"testeo no random"=>Dict("Dias de testeo"=>Frecuenciar), "Repeticiones"=>R, "Porcentaje asintomaticos"=>0.3,
"Probabilidad falso positivo"=>0.01, "Cuarentena"=>Dict("Tipo"=>"solo","Dias atras"=>0,"Dias cuarentena"=>0) )


dict4 = Dict("Tiempo"=>T,"Grupos"=>[[[Dict("Numero de copias"=>10,"Probabilidad de contagio int"=>(0.01*0.2*(1/30)),"Probabilidad de contagio ext"=>0,
"Horario"=>trab[1:2,:])] for i=1:10]' [[Dict("Numero de copias"=>1,"Probabilidad de contagio int"=>0,"Probabilidad de contagio ext"=>0.01,
"Horario"=>Ftrab2[2*j-1:2*j,:]),Dict("Numero de copias"=>1,"Probabilidad de contagio int"=>0,"Probabilidad de contagio ext"=>0.01,
"Horario"=>Ftrab2[20+2*j-1:20+2*j,:])] for j=1:10]']',"Matriz grupos"=>Mrel2,"Politica de testeo"=>"Pool","Testear sintomaticos"=>"no",
"testeo random"=>Dict([]),
"testeo no random"=>Dict("Dias de testeo"=>Frecuenciar4), "Repeticiones"=>R, "Porcentaje asintomaticos"=>0.3,
"Probabilidad falso positivo"=>0.01, "Cuarentena"=>Dict("Tipo"=>"solo","Dias atras"=>0,"Dias cuarentena"=>0) )

# pass data as a json string (how it shall be displayed in a file)


stringdata = JSON.json(dict1)
stringdata2 = JSON.json(dict2)
stringdata3 = JSON.json(dict3)
stringdata4 = JSON.json(dict4)



# write the file with the stringdata variable information

open("parametros1.json", "w") do f
    write(f, stringdata)
end

open("parametros2.json", "w") do f
    write(f, stringdata2)
end
open("parametros3.json", "w") do f
    write(f, stringdata3)
end
open("parametros4.json", "w") do f
    write(f, stringdata4)
end



###################

### Read data #####

###################

# create variable to write the information

dict12 = Dict()
dict22 = Dict()
dict32 = Dict()
dict42 = Dict()


open("parametros1.json", "r") do f

    global dict12
    dicttxt1 = read(f,String)  # file information to string
    dict12=JSON.parse(dicttxt1)  # parse and transform data
end

open("parametros2.json", "r") do f
    global dict22
    dicttxt2 = read(f,String)  # file information to string
    dict22=JSON.parse(dicttxt2)  # parse and transform data

end
open("parametros3.json", "r") do f
    global dict32
    dicttxt3 = read(f,String)  # file information to string
    dict32=JSON.parse(dicttxt3)  # parse and transform data
end

open("parametros4.json", "r") do f
    global dict42
    dicttxt4 = read(f,String)  # file information to string
    dict42=JSON.parse(dicttxt4)  # parse and transform data
end
