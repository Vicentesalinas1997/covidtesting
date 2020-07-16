
	###############################################################################
	#Se realizan las 4 politicas
#politicas=4 #Numero de politicas
#NombresP=["No testear","Simple, Todos los dias","Comparación","Simple, Cada 10 dias"]#Nombres politicas

Caso=string(readdir("Leer")[1])
	mv("Leer/"*Caso,"Guardar/"*Caso)
dict12 = Dict()
ff=1
for json in filter(x -> endswith(x, ".json"), readdir("Guardar/"*Caso))
global N,Group,T,Horario,Mrel, Diastest,Testgrupo,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,test,quienes,Cuaren,Diascuarentena,Diasatras,TestSyn,distribuir,graficogrupo,nombregrafico
open("Guardar/"*Caso*"/"*string(json), "r") do f
    global dict12
    dicttxt1 = read(f,String)  # file information to string
    dict12=JSON.parse(dicttxt1)  # parse and transform data
end

N,Group,T,Horario,Mrel, Diastest,Testgrupo,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,tpool,test,quienes,Cuaren,Diascuarentena,Diasatras,TestSyn,distribuir=leer(dict12)
t_peak=ones(N)*24                 #Tiempo peak de contagio
scalar_asint=1          #Escalar para los asintomaticos
#peak=[(ones(100)*0.01)' (ones(20)*0.2)']' #Prob peak de infeccion
peak=ones(N)*0.2 #Prob peak de infeccion
graficogrupo=[[[i for i in 1:length(Group[1,:])]']  [[i] for i in 1:length(Group[1,:])]']
nombregrafico=[["Total"]' ["Grupo "*string(i) for i in 1:length(Group[1,:])]']'

if ff==1
	global ff, base_mQua,base_mQua2, base_mInf,base_mInf2,base_mInf3, base_mNFp, base_T, base_T2,base_mSy,base_maxInf,base_maxQua, base_maxSy, base_Infect, base_VInf
base_mQua,base_mQua2, base_mInf,base_mInf2,base_mInf3, base_mNFp, base_T, base_T2,base_mSy,base_maxInf,base_maxQua, base_maxSy, base_Infect, base_VInf= simulation(N,Group,T,Horario,Mrel, Diastest,Testgrupo,peak,t_peak,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,tpool,test,quienes,Cuaren,Diascuarentena,Diasatras,scalar_asint,TestSyn,distribuir,[])
end
if ff==2
	global ideal_mQua,ideal_mQua2, ideal_mInf,ideal_mInf2,ideal_mInf3, ideal_mNFp, ideal_T,ideal_T2, ideal_mSy,ideal_maxInf,ideal_maxQua, ideal_maxSy, ideal_Infect, ideal_VInf
ideal_mQua,ideal_mQua2, ideal_mInf,ideal_mInf2,ideal_mInf3, ideal_mNFp, ideal_T,ideal_T2, ideal_mSy,ideal_maxInf,ideal_maxQua, ideal_maxSy, ideal_Infect, ideal_VInf= simulation(N,Group,T,Horario,Mrel, Diastest,Testgrupo,peak,t_peak,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,tpool,test,quienes,Cuaren,Diascuarentena,Diasatras,scalar_asint,TestSyn,distribuir,[])
end
if ff==3
	global f1_mQua,f1_mQua2, f1_mInf,f1_mInf2,f1_mInf3, f1_mNFp, f1_T, f1_T2, f1_mSy,f1_maxInf,f1_maxQua, f1_maxSy,f1_Infect, f1_VInf
f1_mQua,f1_mQua2, f1_mInf,f1_mInf2,f1_mInf3, f1_mNFp, f1_T, f1_T2, f1_mSy,f1_maxInf,f1_maxQua, f1_maxSy,f1_Infect, f1_VInf = simulation(N,Group,T,Horario,Mrel, Diastest,Testgrupo,peak,t_peak,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,tpool,test,quienes,Cuaren,Diascuarentena,Diasatras,scalar_asint,TestSyn,distribuir,[])
end
if ff==4
	global fr_mQua,fr_mQua2, fr_mInf,fr_mInf2,fr_mInf3, fr_mNFp, fr_T, fr_T2,fr_mSy,fr_maxInf,fr_maxQua, fr_maxSy, fr_Infect, fr_VInf
fr_mQua,fr_mQua2, fr_mInf,fr_mInf2,fr_mInf3, fr_mNFp, fr_T, fr_T2,fr_mSy,fr_maxInf,fr_maxQua, fr_maxSy, fr_Infect, fr_VInf = simulation(N,Group,T,Horario,Mrel, Diastest,Testgrupo,peak,t_peak,p_int,p_ext,Porasint,Probfp,Repeticiones,Politica,tpool,test,quienes,Cuaren,Diascuarentena,Diasatras,scalar_asint,TestSyn,distribuir,[])
end
ff+=1
end
###############################################################################
###############Graficos########################################################
if ff==2
	arrays=[[base_mQua2]]
	NombresP=["a"]
	TestP=[sum(base_T)]
	InfectadosP=[sum(base_VInf)]
	CuarentenapeakP=[maximum(sum(base_mQua2,dims=1))]
	InfectadospeakP=[maximum(sum(base_mInf3,dims=1))]
	A=[Dict()]
end
if ff==3
	arrays=[[base_mQua2],[ideal_mQua2]]
	NombresP=["a","b"]
	TestP=[sum(base_T),sum(ideal_T)]
	InfectadosP=[sum(base_VInf),sum(ideal_VInf)]
	CuarentenapeakP=[maximum(sum(base_mQua2,dims=1)),maximum(sum(ideal_mQua2,dims=1))]
	InfectadospeakP=[maximum(sum(base_mInf3,dims=1)),maximum(sum(ideal_mInf3,dims=1))]
	A=[Dict(),Dict()]
end
if ff==4
	arrays=[[base_mQua2],[ideal_mQua2],[f1_mQua2]]
	NombresP=["a","b","c"]
	TestP=[sum(base_T),sum(ideal_T),sum(f1_T)]
	InfectadosP=[sum(base_VInf),sum(ideal_VInf),sum(f1_VInf)]
	CuarentenapeakP=[maximum(sum(base_mQua2,dims=1)),maximum(sum(ideal_mQua2,dims=1)),maximum(sum(f1_mQua2,dims=1))]
	InfectadospeakP=[maximum(sum(base_mInf3,dims=1)),maximum(sum(ideal_mInf3,dims=1)),maximum(sum(f1_mInf3,dims=1))]

	A=[Dict(),Dict(),Dict()]
end
if ff==5
	arrays=[[base_mQua2],[ideal_mQua2],[f1_mQua2],[fr_mQua2]]
	NombresP=["a","b","c","d"]
	TestP=[sum(base_T),sum(ideal_T),sum(f1_T),sum(fr_T)]
	InfectadosP=[sum(base_VInf),sum(ideal_VInf),sum(f1_VInf),sum(fr_VInf)]
	CuarentenapeakP=[maximum(sum(base_mQua2,dims=1)),maximum(sum(ideal_mQua2,dims=1)),maximum(sum(f1_mQua2,dims=1)),maximum(sum(fr_mQua2,dims=1))]
	InfectadospeakP=[maximum(sum(base_mInf3,dims=1)),maximum(sum(ideal_mInf3,dims=1)),maximum(sum(f1_mInf3,dims=1)),maximum(sum(fr_mInf3,dims=1))]
	A=[Dict(),Dict(),Dict(),Dict()]
end
cont=1
for json in filter(x -> endswith(x, ".json"), readdir("Guardar/"*Caso))
	global cont
NombresP[cont]=json[1:end-5]
A[cont]=Dict("Politica"=>NombresP[cont],"Test"=>Int.(round.(TestP[cont])),"Infectados"=>Int.(round.(InfectadosP[cont])),"PeakQua"=>Int.(round.(CuarentenapeakP[cont])),"PeakInf"=>Int.(round.(InfectadospeakP[cont])),"Caso"=>[])
cont+=1
end
A[1]["Caso"]=Caso
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
titulo="Infectados_Trabajando"
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
titulo="Infectados_Acumulado"
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
titulo="Test"
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
titulo="Test_rapidos"
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


############################################################################



tmpl = mt"""

\documentclass[opre]{PDF/informs3}
\usepackage[utf8]{inputenc}
\usepackage[english]{babel}

\usepackage[dvipsnames]{xcolor}
\usepackage{multirow}
\usepackage{verbatim,wrapfig,latexsym,graphicx,psfrag,amsfonts,amsmath,amssymb}
\usepackage{hyperref,multirow,color,subfigure,comment,pdfsync,setspace,thmtools,thm-restate}
\usepackage{enumerate,enumitem,tikz,pdfpages,xspace,multirow,algorithmic,algorithm}


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


\newcommand{\brown}[1]{{\color{brown} #1}}

\TheoremsNumberedThrough
\ECRepeatTheorems

\EquationsNumberedThrough

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
 \hspace{7cm} Santiago, \today \\
\vspace{1cm}

\begin{center}
\noindent \textbf{{\Reporte comparación de politicas}}\\
\end{center}
\textbf{Comparación de Políticas}\\
Comparando las politicas se obtiene que: \\
{{#:A}} {{Politica}} tiene  {{Test}} test realizados\\ {{/:A}}.

\begin{center}
    \includegraphics[scale=0.45]{Guardar/{{#:A}}{{Caso}}{{/:A}}/Test1.pdf}
\end{center}
Comparando las politicas se obtiene que:\\
{{#:A}} {{Politica}} tuvo  {{Infectados}} infectados\\ {{/:A}}.

\begin{center}
    \includegraphics[scale=0.45]{Guardar/{{#:A}}{{Caso}}{{/:A}}/Infectados_Acumulado1.pdf}
\end{center}
Comparando las politicas se obtiene que:\\
{{#:A}} {{Politica}} tiene  {{PeakQua}} personas en cuarentena el día de mayor peak\\ {{/:A}}.

\begin{center}
    \includegraphics[scale=0.45]{Guardar/{{#:A}}{{Caso}}{{/:A}}/Cuarentena1.pdf}
\end{center}
Comparando las politicas se obtiene que:\\
{{#:A}} {{Politica}} tiene  {{PeakInf}} nuevos infectados el día de mayor peak\\ {{/:A}}.

\begin{center}
    \includegraphics[scale=0.45]{Guardar/{{#:A}}{{Caso}}{{/:A}}/Infectados1.pdf}
\end{center}


\end{document}

"""

  rendered = render(tmpl,A=A)
  filename = string("Guardar/"*Caso*"/Reporte"*Caso*".tex")
  open(filename, "w") do file
    write(file, rendered)
  end
  run(`pdflatex $filename`)

 mv("Reporte"*Caso*".pdfsync","Guardar/"*Caso*"/Reporte"*Caso*".pdfsync")
 mv("Reporte"*Caso*".out","Guardar/"*Caso*"/Reporte"*Caso*".out")
 mv("Reporte"*Caso*".pdf","Guardar/"*Caso*"/Reporte"*Caso*".pdf")
 mv("Reporte"*Caso*".log","Guardar/"*Caso*"/Reporte"*Caso*".log")
 mv("Reporte"*Caso*".aux","Guardar/"*Caso*"/Reporte"*Caso*".aux")
#t = Timer((t) -> (if length(readdir("Leer"))>=1 include("covid12.jl") else println("esperando") end), 1; interval=10)
