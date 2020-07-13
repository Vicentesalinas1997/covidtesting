using Mustache
Caso="usuario0"
A = [Dict("Politica" => "nombre Politica 1", "Infectados" => string(1),"Caso"=>Caso),
Dict("Politica" => "nombre Politica 2", "Infectados" => string(2),"Caso"=>[]),
Dict("Politica" => "nombre Politica 3", "Infectados" => string(3),"Caso"=>[]),
Dict("Politica" => "nombre Politica 4", "Infectados" => 4,"Caso"=>[])
     ]

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
Comparando las politicas se obtiene que: {{#:A}} {{Politica}} tiene  y {{Infectados}} infectados{{/:A}}.\\
{{#:A}}{{Caso}}{{/:A}}
\begin{center}
    \includegraphics[scale=0.45]{Guardar/{{#:A}}{{Caso}}{{/:A}}/Infectados 1.pdf}
\end{center}

\textbf{Motivación}\\
Debido a la necesidad observada por parte del equipo en encontrar y comparar más estrategias de prevención de contagios por covid, hemos desarrollado un programa de simulación, el cual busca ser un gran aporte para la sociedad permitiendo tomar decisiones más informados comparando estrategias que se adapten mejor a los intereses y necesidades.\\
\textbf{Problema}\\

Uno de los mayores problemas que hemos detectado es la función principal que se les esta dando a los test PCR, pues es para diferenciar los casos sospechosos y para confirmar los casos positivos, se cree que esto es por los costos que tiene realizar estas pruebas y por lo tanto se ha priorizado estos usos. Pero esto trae consigo consecuencias como son una detección tardía de los casos sintomáticos, pues según un estudio de nature ("aquí habría que citar al paper de nature") el momento de mayor contagio de los infectados en unos días previos a la aparición de los síntomas. Otra consecuencia es la aplicación de medidas de cuarentenas a grupos de contacto estrecho, con el fin de reducir el número de contagios, pero para los casos de trabajo esencial estas disminuciones de personal pueden ser sumamente costosas y perjudiciales para el funcionamiento de Hospitales, recintos de salud, entre otras.\\
\textbf{Que son los Pool testing} (Dependiendo el público al que esta dirigido deberíamos ver si lo explicamos)\\

Los Pool testing consisten en pruebas PCR a grupos de individuos a diferencia de las tradicionales que se toman de manera individual, esto permite con un menor costo abarcar en la mayoría de casos un mayor número de pruebas a individuos. La validación de esta técnica de pruebas, se encuentra detallada en el ("paper del estudio de la chile"). Al momento de utilizar este tipo de test es de suma importancia un conocimiento previo de la prevalencia esperada para escoger adecuadamente el tamaño de estos grupos, pues en caso de que uno de estos test entrege un resultado positivo habrá que realizar de manera individual el test a todos los integrantes de esa muestra con el fin identificar a la o las personas infectadas. \\
\textbf{Beneficios de la aplicación  de pruebas periódicas preventivas con Pool testing }\\

Es sumamente beneficiosa para la sociedad, pues en primer lugar permite detectar casos asintomáticos, esto consideramos que es algo sumamente importante, pues los casos asintomáticos se estima que son al menos un $30\% + \int_0^{\infty}$  de la población, sobre estos porcentajes no existe un mayor acuerdo, pues es muy difícil de determinar con las estrategias de pruebas actuales, pero la aplicación de pruebas preventivas traerá consigo una mayor información al ir descubriendo estos casos. Otro beneficio es la posibilidad de una detección temprana para los casos sintomáticos, es decir, antes de que estos alcancen su peak de infecciosidad. Esto puede ser un beneficio crucial, pues a diferencia del peak de las otras enfermedades respiratorias que es posterior a los síntomas, una cuarentena posterior a los síntomas tiene un impacto mucho menor. El uso de pool testing es la principal razón de la factibilidad de estas pruebas preventivas, pues en escenarios de baja prevalencia permite una disminución considerable en el número de test utilizados, pues al contrario de los test tradicionales en los cuales se realiza un test por persona, los pool testing permiten realizar un test a 5 o mas personas en simultaneo con una efectividad   similar. Esto es clave para la aplicación de pruebas regulares, pues el realizar estas pruebas involucra un aumento en el número de test, pero junto con la aplicación de pool testing estos números son sumamente factibles. Estas pruebas preventivas abren una gama de posibles combinaciones y estrategias que permiten tomar una mejores decisiones al momento de decidir sobre las cuarentenas y el uso estratégico de los test que se dispongan, pudiendo disminuir el personal en cuarentena al evaluar los riesgos, probar diferentes combinaciones de turnos de trabajo que al combinarlas con la aplicación periódica de test, genere una mayor seguridad para los trabajadores, entre otras.\\
\textbf{Simulaciones de políticas}\\

Para realizar estas comparaciones hemos desarrollado un programa el cual con algunos supuestos permite la simulación de un escenario con una cantidad predefinidas de trabajadores, en el cual calibrando las probabilidades de contagio basándose en los riesgos a los que estos trabajadores estarán expuestos, permite simular y comparar distintas políticas de cuarentena (grupos de contacto estrecho, individuos solos, escoger la duración de las cuarentenas) y diversas decisiones a la hora de distribuir y tomar los test (aplicar test grupales, test individuales, test a las personas que presentan síntomas, tomar test todas las semanas, entre otras). Esto permite a los usuarios tomar una mejor decisión, pues no tendrán que apostar por alguna estrategia sin medir los impactos que esta puede traer, pues podrán evaluar de manera segura las ventajas y desventajas de una política sobre otra, ya sea en el número de test esperado que se utilice, en la cantidad máxima de infectados, en las personas en cuarentena, entre otras.\\
\textbf{Explicación del seudocódigo-simulación (supuestos, parámetros-instancias)
}\\

\textbf{El programa presenta los siguiente supuesto:}\\
- Se considerara que la gente infectada no puede volver a infectarse.
- Si un individuo presenta síntomas por más de dos días, este sera enviado a cuarentena sin necesidad de realizar un test.\\
- Las personas que entran en cuarentena no son reemplazadas.\\
- Resultados de test, infecciocidad, incubación, entre otros. Son independientes de la manifestación de síntomas.\\
- Los resultados del test se  conocen al final del día.\\
- Tiempo de incubación (lognormal) y recuperación (uniforme 2-4 semanas).\\
- La probabilidad de resultado positivo sigue curva publicada en Jama, ajustada a la evolución de la infección en el paciente.\\
-La infecciosidad es una versión escalada de la curva de positividad.\\
\textbf{Principales parámetros para generar una simulación, estos son valores que se pueden modificar a interés del usuario:}\\
-\textbf{Matriz de turnos de trabajo:} día y noche.\\
-\textbf{Matriz de interacciones:} La coordenada i,j contiene la probabilidad de interacción (contagio) entre un individuos el grupo i con uno del grupo j. \\
-\textbf{Grupos:} Conjunto de individuos que coinciden en la posibilidad de interacción con los grupos.\\
-\textbf{Matriz de testeo:} Días en que se testea cada individuo.\\
-\textbf{Tipo de test:} individual, grupal, no testear.\\
-\textbf{Definir previamente los test:} La simulación permite definir previamente a que individuos testear cada día o si se desea entregar un porcentaje por grupo y de manera aleatoria se escogerá este. Siempre sujeto a si la persona esta en cuarentena.\\
-\textbf{Eficiencia test:} prob. falsos positivos (scalar) y falsos negativos (curva).\\
-\textbf{Asintomáticos:} (scalar, estático) fracción de la población.\\
-\textbf{Políticas de cuarentena:} al resultado positivo, al grupo de contactos estrecho  del positivo y/o por cierta cantidad de días.\\

Algunos parámetros más que se pueden modificar son la probabilidades de contagio por factores externos, si se desea solo testear al personal trabajando, la cantidad de días de simulación ,entre otros más. Con estos datos se realiza una simulación
en la cual se ve la evolución diaria del número de personas infectados (no descubiertos, en cuarentena y en total), personas en cuarentena, personas con síntomas y número de test realizados. Mediante la siguiente secuencia temporal:
Cada día\\
En t=0, designa qué personas serán asintomáticas.\\
Se calculan las probabilidades de infección de todos, para ese día.\\
Se realizan los contagios de los trabajadores (tiempos individuales de incubación y recuperación) .\\
Se ejecuta la política de testeo. Dependiendo de los resultados, y desarrollo de síntomas, se ejecuta la política de cuarentena
el número de infectados, personas en cuarentena, personas con síntomas y test realizados.\\
Se procede al siguiente día.\\
En base a esta información los usuarios podrían saber que políticas les son factibles y entre ellas cual se adapta mejor a sus intereses. \\
\textbf{ejemplificar con 2 instancias (asilos y centros médicos)}\\
\textbf{Asilo}\\
Para simular el piloto realizado por el ISCI en los asilos "......" y ".....", se utilizaron los siguientes parámetros: \\
- \\
-\\
-\\
-\\
-\\
Con esto se obtuvieron los siguientes gráficos:\\ (Para graficar falta definir un horario para los funcionarios y una matriz de interacción)
-\\
-\\
-\\
En estos se puede observar que ......

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\end{document}

"""

  rendered = render(tmpl,(A=A))
  filename = string("PDF/Reporte.tex")
  open(filename, "w") do file
    write(file, rendered)
  end
  run(`pdflatex $filename`)
