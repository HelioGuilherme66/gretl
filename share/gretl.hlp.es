#
add
@Contrastes
(A�adir)
A�ade variables a un modelo y contrasta su significatividad.

Se a�aden las variables seleccionadas al modelo anterior y se
estima el nuevo modelo. Si se a�ade m�s de una variable, se presenta el
estad�stico F para el contraste conjunto de significaci�n de las
variables a�adidas (s�lo para el m�todo MCO) junto con su valor p. Un
valor p inferior a 0.05 indica que los coeficientes son conjuntamente
significativos al nivel del 5 por ciento.

#
adf
@Contrastes
Contraste aumentado de Dickey-Fuller

Esta orden requiere un orden de retardos entero

Calcula los estad�sticos para dos contrastes de Dickey-Fuller. En ambos
casos la hip�tesis nula es que la variable en cuesti�n presenta una ra�z
unitaria.

El primero es un contraste t basado en el modelo

(1 - L)x(t) = m + g * x(t-1) + e(t).

La hip�tesis nula es que g = 0.

En el segundo contraste (aumentado) se estima una regresi�n no
restringida (teniendo como regresores una constante, una tendencia
temporal, el primer retardo de la variable y "orden" retardos de la
primera diferencia) y una regresi�n restringida (quitando la tendencia
temporal y el primer retardo). El estad�stico de contraste es F
calculado como

[(SCRr - SCRnr)/2]/[SCRnr/(T - k)]

donde T es el tama�o muestral y k el n�mero de par�metros del modelo no
restringido. Es necesario tener en cuenta que los valores cr�ticos para
estos estad�sticos no son los usuales.


#
ar
@Estimaci�n
Estimaci�n generalizada (autorregresiva) de Cochrane-Orcutt

Calcula las estimaciones de un modelo utilizando el procedimiento
iterativo generalizado de Cochrane-Orcutt. Las iteraciones se terminan
cuando las sucesivas sumas de cuadrados residuales no var�an m�s del
0.005 por ciento o cuando han transcurrido 20 iteraciones.

En la 'lista de retardos AR' se especifica la estructura del proceso de
error. Por ejemplo, la entrada "1 3 4" corresponde a

     u(t) = rho1*u(t-1) + rho3*u(t-3) + rho4*u(t-4) + et


#
arch
@Contrastes
Realiza un contraste de ARCH (Autoregressive Conditional
Heteroskedasticity)

Esta instrucci�n requiere un orden de retardos entero.

Contrasta la existencia de un ARCH del orden especificado, en el
modelo. Si el estad�stico LM tiene un valor p inferior a 0.10, entonces
se realiza la estimaci�n del ARCH. Si, en la regresi�n auxiliar,  la
varianza estimada de cualquier observaci�n no es positiva, entonces se
usa el correspondiente residuo al cuadrado. Despu�s se estima el modelo
original utilizando m�nimos cuadrados ponderados.



#
boxplots
@Gr�ficos
(Gr�ficos de caja)
An�lisis exploratorio de los datos.

Estos gr�ficos (debidos a Tukey y Chambers) presentan la distribuci�n de
una variable. La caja central incluye el 50 por ciento de los datos que
est�n en el medio de la distribuci�n, es decir, est� delimitada por el
primer y el tercer cuartiles. Las "patillas" se extienden desde esta
caja hasta los valores m�nimo y m�ximo. Se dibuja una l�nea a lo largo de
la caja en la situaci�n de la mediana.

En el caso de cajas recortadas, el recorte muestra los l�mites de un
intervalo de confianza aproximada del 90 por ciento para la mediana.
Este se obtiene mediante el m�todo bootstrap, que puede tardar un rato
si la serie es muy larga.

Haciendo "click" con el rat�n en ventana del gr�fico de caja aparece un
men� que permite guardar el gr�fico en formato postscript encapsulado
(EPS) o como fichero postscript de p�gina completa. Bajo el sistema X
windows (linux) tambi�n se puede guardar el gr�fico en formato XPM; bajo
MS Windows se puede copiar al porta-papeles como mapa de bits.

El men� tambi�n ofrece la posibilidad de abrir una ventana de texto que
muestra un "resumen de 5 n�meros" (m�nimo, primer cuartil, mediana,
tercer cuartil, m�ximo) y, si se ha elegido la opci�n de gr�fico
"recortado", un intervalo de confianza para la mediana.

Se pueden controlar algunos detalles de los "gr�ficos de caja" de gretl
por medio de un fichero (de texto plano) de nombre .boxplotrc que se
busca sucesivamente en el directorio de trabajo actual, el directorio
home del usuario (que corresponde a la variable de entorno HOME) y el
directorio gretl del usuario (que se puede observar y cambiar en
Archivo, Preferencias, Men� general). Las opciones que pueden cambiarse
de esta manera son: la fuente a usar al producir una salida postscript
(debe ser un nombre de fuente gen�rica postscript v�lido; por defecto es
Helv�tica), el tama�o de la fuente en puntos (tambi�n para la salida
postscript; por defecto es 12), el m�nimo y el m�ximo para el rango
del eje y, la anchura y altura del gr�fico en pixels (por defecto, 560 x
448), si se deben imprimir los valores num�ricos para los cuartiles y la
mediana (por defecto, no imprimirlos), y si los 'outliers' (puntos que
quedan m�s all� de 1.5 veces el rango intercuart�lico desde la caja
central) deber�an indicarse de forma separada (por defecto, no).
Por ejemplo:

font = Times-Roman
fontsize = 16
max = 4.0
min = 0
width = 400
height = 448
numbers = %3.2f
outliers = true

En la segunda l�nea del final, el valor asociado a "numbers" es una
cadena de formato "printf" como en el lenguaje de programaci�n C;
cuando se especifica, esto controla la inscripci�n de la mediana y los
cuartiles cerca del gr�fico de caja, si no se da una opci�n "numbers",
esos valores no se escriben. En el ejemplo, los valores se escribir�n
a tres d�gitos, con dos d�gitos despu�s del punto decimal.

No es necesario establecer todas las opciones, y el orden no importa.

Las l�neas que no sigan el modelo "clave=valor" son ignoradas, como las
l�neas que comiencen con la almohadilla, #.

Despu�s de cada variable especificada en la orden 'gr�fico de
caja' puede a�adirse una expresi�n booleana entre par�ntesis, para
delimitar la muestra para la variable en cuesti�n. Debe insertarse un
espacio entre el nombre o n�mero de la variable y la expresi�n.
Supongamos que tenemos datos de 'salarios' para hombres y mujeres y
tenemos una variable ficticia 'G�NERO' que toma valor 1 para los hombres
y 0 para las mujeres. En este caso, podemos representar gr�ficos de caja
comparativos, con las siguientes l�neas en el di�logo de la orden
'gr�fico de caja':

  salarios (G�NERO=1) salarios (G�NERO=0)

#
chow
@Contrastes
Contraste de Chow de homogeneidad estructural

Esta orden necesita un n�mero de observaciones (o fechas, si los datos
tienen fecha).

Debe ejecutarse despu�s de una regresi�n MCO. Crea una variable ficticia


que es igual a 1 desde el punto de ruptura que se especifique hasta el
final de la muestra y 0 en el resto. Tambi�n crea t�rminos de
interacci�n entre esta variable ficticia y las variables independientes
originales. Se ejecuta una regresi�n aumentada incluyendo esos t�rminos
y se calcula un estad�stico F tomando la regresi�n aumentada como 'no
restringida' y la original como 'restringida'. Este estad�stico es
adecuado para contrastar la hip�tesis nula de que no hay cambio
estructural en el citado punto de ruptura.


#
coint
@Contrastes
Contraste de cointegraci�n

Con esta orden (que necesita un orden de retardos entero) se realizan
los contrastes de Dickey-Fuller de la hip�tesis nula de que cada una de
las variables seleccionadas tiene una ra�z unitaria, considerando el
orden de retardos dado. Se estima la regresi�n cointegrante y se realiza


un contraste ADF sobre los residuos de la regresi�n.

Tambi�n se ofrece el estad�stico de Durbin-Watson de la regresi�n
cointegrante.

(Hay que se�alar que para ninguno de estos estad�sticos pueden
utilizarse las tablas estad�sticas usuales)


#
compact
@Conjunto de datos
(Compactar)
Escribiendo datos a una frecuencia inferior.

Cuando se a�ade a un conjunto de datos una serie que es de una
frecuencia superior, es necesario "compactar" esa nueva serie. Por
ejemplo, una serie mensual tendr� que ser "compactada" para introducirla


en un conjunto de datos trimestral. Se ofrecen tres opciones para el
"compactado":

1. Promediado: el valor escrito en el conjunto de datos ser� la media
aritm�tica de los valores de la serie en cuesti�n. Por ejemplo, el valor


introducido para el primer trimestre de 1990 ser� la media de los
valores de enero, febrero y marzo de 1990.

2. Valores de 'final de periodo': el valor escrito en el conjunto de
datos es el �ltimo valor del periodo correspondiente en los datos
de m�s alta frecuencia. Por ejemplo, en el primer trimestre de 1990 se
introducir�a el valor de marzo de 1990.

3. Valores de 'principio de periodo': el valor escrito en el conjunto de
datos es el primer valor del periodo correspondiente en los datos de
m�s alta frecuencia. Por ejemplo, en el primer trimestre de 1990 se
introducir�a el valor de enero de 1990.



#
corc
@Estimaci�n
Modelo de Cochrane-Orcutt

Esta orden calcula las estimaciones de un modelo usando el procedimiento
iterativo de Cochrane-Orcutt. Las iteraciones acaban cuando dos valores
sucesivos de rho no difieren en m�s de 0.001 o cuando se han realizado
ya 20 iteraciones.

La regresi�n transformada final se estima para el rango de observaci�n
primobs+1 ultobs que est� actualmente en efecto.


#
dialog box for models
@Estimaci�n
(cuadro de di�logo para los modelos)
Para elegir la variable dependiente, seleccione una variable en la lista


de la izquierda y presione el bot�n "Elegir->" que se�ala a la caja
correspondiente a "variable dependiente". Si se activa el cuadro
"Establecer por defecto" la variable que se elija ser� preseleccionada
como variable dependiente la pr�xima vez que se abra el cuadro de
di�logo de modelos. Atajo: haciendo doble "click" sobre una variable
del lado izquierdo se selecciona y se establece como variable
dependiente por defecto.

Para elegir las variables independientes, selecci�nelas a la izquierda y


presione el bot�n "A�adir->" (o haga "click" con el bot�n derecho del
rat�n). Se pueden seleccionar varias variables contiguas arrastr�ndolas
con el rat�n. Se puede seleccionar un grupo de variables no
contiguas haciendo click sobre ellas mientras se pulsa la tecla Ctrl.


#
diff
@Transformaciones
(Diferencia regular)

Se calcula la primera diferencia de cada variable de la lista
dada y el resultado se guarda en una nueva variable con prefijo "d_".
As� por ejemplo, la nueva variable
d_x = x(t) - x(t-1).

#
export
@Conjunto de datos
Exporta datos de gretl a otros formatos.

Se pueden exportar datos en formato CSV (Valores separados por comas):
estos datos pueden abrirse desde hojas de c�lculo y muchos otros
programas.

Tambi�n se pueden exportar datos a los formatos nativos de GNU R y GNU
Octave. Para m�s informaci�n sobre estos programas (ambos desarrollan
an�lisis estad�sticos avanzados) por favor visite sus propias
p�ginas web, http://www.r-project.org/ y http://www.octave.org/


#
factorized plot
@Gr�ficos
(Gr�fico con factor de separaci�n)

Esta orden necesita que Vd elija tres variables y la �ltima de ellas
ha de ser una variable ficticia (con valor 1 � 0). La variable Y se
representa con respecto a la variable X, con los puntos coloreados de
forma diferente dependiendo del valor de la tercera variable.

Por ejemplo: si vd tiene datos de salarios y nivel de educaci�n para una
muestra de varias personas; y Vd tambi�n tiene una variable ficticia con
valor 1 para los hombres y 0 para las mujeres (como en el fichero que se
suministra con gretl, data7-2). Un "gr�fico factorizado" de SALARIOS
con respecto a EDUCACI�N usando la variable SEXO como factor mostrar�
los puntos para los hombres en un color y para las mujeres en otro (con
una leyenda para identificarlos).


#
genr
@Transformaciones
Genera una nueva variable
Uso:             nuevo_nombre_variable = transformaci�n

Crea nuevas variables, normalmente por medio de transformaciones de
variables ya existentes. Ver tambi�n diff, logs, lags, ldiff, multiply
y square como atajos.

Los operadores matem�ticos que se soportan son, en orden de precedencia:
^ (exponenciaci�n); *, / y % (m�dulo o resto); + y -.

Los operadores booleanos (de nuevo en orden de precedencia) son ! (NO
l�gico), & (Y l�gico), | (O l�gico), >, <, = y != (No igual). Los
operadores booleanos pueden usarse al construir variables ficticias: por
ejemplo (x>10) devuelve 1 si x(t)>10 y en caso contrario 0.

Las funciones que se soportan pertenecen a estos grupos:

- Funciones matem�ticas standard: abs, cos, exp, int (parte entera), ln
(logaritmo natural: log es un sin�nimo), sin (seno), sqrt (ra�z
cuadrada).

- Funciones estad�sticas: mean (media aritm�tica), median (mediana), var
(varianza), sd (desviaci�n t�pica o estandard), sum, cov (covarianza),
corr (coeficiente de correlaci�n, min (m�nimo), max (m�ximo).

- Funciones de series temporales: lag (retardo), lead (adelanto), diff
(primera diferencia), ldiff (log-diferencia, o primera diferencia del
logaritmo natural).

- Miscel�neas: cum (acumulaci�n), sort (ordenaci�n), uniform
(distribuci�n uniforme), normal (distribuci�n normal), missing (devuelve
1 si la variable tiene la observaci�n perdida, en caso contrario 0),
misszero (reemplaza el c�digo de observaci�n perdida por un 0), zeromiss
(operaci�n inversa de misszero).

Todas las funciones anteriores, a excepci�n de cov, corr, uniform y
normal, toman como �nico argumento o el nombre de una variable (hay
que notar que, en una orden genr,  no es posible referirse a las
variables utilizando su n�mero de ID) o una expresi�n compuesta que se
eval�a en una variable (p.ej. ln((x1+x2)/2)). cov y corr requieren dos
argumentos (dos variables) y devuelven, respectivamente, la covarianza y
el coeficiente de correlaci�n entre las dos variables mencionadas.
uniform() y normal() no tienen argumentos y devuelven, respectivamente,
series pseudoaleatorias obtenidas a partir de las distribuciones
uniforme (0-100) y normal standard (ver tambi�n la instrucci�n seed).

Hay varias variables que se definen internamente al ejecutar una
regresi�n, que pueden usarse tambi�n en transformaciones, como son:

  $ess         suma de cuadrados de los residuos
  $rsq         R-cuadrado no corregido
  $T           n�mero de observaciones utilizado por el modelo
  $df          grados de libertad
  $trsq        TR^2 (T veces el R-cuadrado, siendo T el tama�o
               muestral)
  $sigma       desviaci�n t�pica de los residuos
  $lnl         log-verosimilitud (en modelos logit y probit)
  coeff(var)   coeficiente estimado de var
  stderr(var)  desviaci�n t�pica estimada del estimador de var
  rho(i)       coeficiente autorregresivo de i�simo orden de los
               residuos
  vcv(xi,xj)   covarianza entre los coeficientes de las variables
               xi y xj

La variable interna $nobs contiene el n�mero de observaciones del
dominio muestral actual, que puede ser o puede no ser igual al $T del
�ltimo modelo.

La variable interna $pd contiene la periodicidad o frecuencia de los datos (por 
ejemplo, 4 para datos trimestrales, 12 para mensuales).

La variable interna t hace referencia a las observaciones, comenzando en
1. As�, es posible hacer "genr dum15 = (t=15)" para generar una variable
ficticia con valor 1 para la observaci�n 15 y 0 para el resto.

Ejemplos de f�rmulas v�lidas:

   y = x1^3          [x1 al cubo]
   y=ln((x1+x2)/x3)  [argumento compuesto a funci�n ln]
   z=x>y             [hace z(t) igual a 1 si x(t) > y(t) en caso
                      contrario 0]
   y=x(-2)           [x retardada 2 periodos]
   y=x(2)            [x adelantada 2 periodos]
   y = mean(x)       [media aritm�tica]
   y = diff(x)       [y(t) = x(t) - x(t-1)]
   y = ldiff(x)      [y = ln(x(t)) - ln(x(t-1))]
                      ldiff(x) es la tasa instant�nea de
                      crecimiento de x.
   y = sort(x)       [ordena x en orden creciente y la guarda en y]
   y = -sort(-x)     [ordena x en orden decreciente]
   y = int(x)        [trunca x y guarda su valor entero como y]
   y = abs(x)        [guarda los valores absolutos de x]
   y = sum(x)        [suma los valores de x excluyendo las entradas
                      de los valores perdidos -999]
   y = cum(x)        [acumula x: y(t) es la suma de x hasta t]
   aa = $ess         [aa = suma de cuadrados de los residuos de la
                      �ltima regresi�n]
   x = coeff(sqft)   [recoge el coeficiente de sqft del
                      �ltimo modelo]
   rho4 = rho(4)     [recoge el coeficiente autorregresivo de
                      4� orden del �ltimo modelo (se supone un modelo
                      ar)]
   cv=vcv(x1, x2)    [covarianza de los coeficientes de x1 y x2
                      en el �ltimo modelo]
   x=uniform()/100   [variable pseudoaleatoria uniforme, rango 0 a 1]
   x=3*normal()      [variable pseudoaleatoria normal, con media 0
                      y desviaci�n t�pica 3]

Sugerencias sobre variables ficticias:

* Supongamos que x se codifica con valores 1, 2 o 3 y Vd desea tres
variables ficticias, d1 si x=1, 0 en caso contrario, d2=1 si x=2, 0 en
caso contrario y as� sucesivamente. Para crearlas, use las f�rmulas
d1 = (x=1), d2 = (x=2), y d3 = (x=3).
* Para obtener z = m�x(x,y) genere d=x>y y despu�s
z=(x*d)+(y*(1-d))

#
graphing
@Gr�ficos
(Gr�ficos)
generando gr�ficos de varios tipos

Gretl llama a un programa aparte, gnuplot, para generar los gr�ficos.
Gnuplot es un programa gr�fico de m�ltiples caracter�sticas con
miles de opciones. Gretl le proporciona a Vd acceso directo, v�a una
interface gr�fica, a s�lo un peque�o subconjunto de esas opciones e
intenta elegir para vd los valores adecuados; tambi�n permite que vd
tome completamente el control sobre los detalles del gr�fico si as� lo
desea.

Bajo MS Windows vd puede hacer click en la esquina de arriba a la
izquierda de la ventana del gr�fico, obteniendo as� un men� porta-papeles
le permite elegir varias cosas (incluyendo copiar el gr�fico al
porta-papeles de Windows y enviarlo a la impresora).

Para tener un control completo sobre el gr�fico, siga este
procedimiento:

- Cierre la ventana del gr�fico.
- Desde el men� de sesi�n, seleccione "A�adir el �ltimo gr�fico".
- En la ventana de iconos de sesi�n, haga click con el bot�n
derecho del rat�n sobre el icono del nuevo gr�fico y seleccione o
"Editar usando GUI" o "Editar las �rdenes de gr�fico". La entrada de
" Editar usando GUI" abre un controlador gr�fico para gnuplot que le
permite refinar varios aspectos del gr�fico. La entrada de "Editar las
�rdenes de gr�fico" abre una ventana de editor que contiene el fichero
actual de instrucciones de Gnuplot para generar el gr�fico: esto le
proporciona a vd un control completo sobre los detalles del gr�fico --si
vd conoce algo sobre gnuplot. Para m�s detalles,ver
http://ricardo.ecn.wfu.edu/gnuplot.html or www.gnuplot.org.

#
hccm
@Estimaci�n
(mcch)
Matriz de covarianzas consistente ante heterocedasticidad

Esta orden ejecuta una regresi�n donde los coeficientes se estiman
por medio del procedimiento MCO standard, pero las desviaciones t�picas
de los estimadores de los coeficientes se calculan de una manera que es
robusta ante la heterocedasticidad. Concretamente, se usa el
procedimiento "jacknife" de MacKinnon-White


#
hilu
@Estimaci�n
M�todo de Hildreth-Lu

Calcula las estimaciones de un modelo utilizando el procedimiento de
b�squeda de Hildreth-Lu (refinado mediante el m�todo de CORC
[Cochrane-Orcutt]). Se representa la suma de cuadrados de los residuos
del modelo transformado con respecto a los valores de rho
desde -0.99 hasta 0.99. La regresi�n final transformada  se calcula para
el rango de observaci�n primobs+1 ultobs que est� actualmente en efecto.


#
hsk
@Estimaci�n
(correcci�n de heterocedasticidad)
Estimaciones corregidas de Heterocedasticidad

Se ejecuta una regresi�n MCO y se guardan los residuos. El logaritmo
del cuadrado de dichos residuos se constituye como variable
dependiente en una regresi�n auxiliar, en cuyo lado derecho de la
ecuaci�n est�n las variables independientes originales y sus cuadrados.
Los valores ajustados en la regresi�n auxiliar se usan entonces para
construir una serie de ponderaciones y el modelo original se reestima
utilizando m�nimos cuadrados ponderados. Este resultado final es el que
aparece en el cuadro de resultados.

La serie de ponderaciones se forma como 1/sqrt(exp(fit)), donde "fit"
representa a los valores ajustados obtenidos de la regresi�n auxiliar.

#
lad
@Estimaci�n
(Estimador de m�nima desviaci�n absoluta)
Uso:          lad vardep varindeps

Calcula una regresi�n que minimiza la suma de las desviaciones absolutas
entre los valores observados y los valores ajustados de la variable dependiente. 
Las estimaciones de los coeficientes se obtienen utilizando el algoritmo simplex 
de Barrodale-Roberts; se muestra un aviso si la soluci�n no es �nica. Las 
desviaciones t�picas se obtienen utilizando un m�todo 'bootstrap' con 500
iteraciones.

#
lags
@Transformaciones
(retardos)

Crea nuevas variables que son valores retardados de cada una de las
variables de la lista que se suministra. El n�mero de contrapartes
retardadas para cada una de las variables listadas es igual a la
periodicidad de los datos. Por ejemplo, si la periodicidad es 4 (datos
trimestrales), se crear�n cuatro t�rminos retardados; si en la lista que
se ha suministrado est� la variable "x", la orden crea x_1 = x(t-1),
x_2 = x(t-2), x_3 = x(t-3) y x_4 = x(t-4).


#
ldiff
@Transformaciones

Se obtiene la primera diferencia del logaritmo natural de cada variable
de la lista suministrada y el resultado se guarda en una nueva variable
con el prefijo "ld_".

As� por ejemplo, la nueva variable ld_x = ln[x(t)] - ln[x(t-1)].

#
logit
@Estimaci�n
Regresi�n Logit

La variable dependiente deber�a ser una variable binaria. Se utiliza
el m�todo de m�nimos cuadrados iterativos (el m�todo EM o de
expectativa-maximizaci�n) para obtener las estimaciones
m�ximo-veros�miles de los coeficientes de las variables independientes.
Como el modelo es no lineal, las pendientes dependen de los valores de
las variables independientes: las pendientes que gretl muestra se
eval�an en las medias de dichas variables. El estad�stico Chi-cuadrado
contrasta la hip�tesis nula de que todos los coeficientes, excepto la
constante, son cero.


#
logs
@Transformaciones

Se calcula el logaritmo natural de cada una de las variables de la lista


que se suministra y el resultado se guarda en una nueva variable con 
prefijo "l_". As� por ejemplo la nueva variable l_x = ln(x).


#
loop
@Programaci�n
(Bucle)
instrucciones repetidas

Uso:            loop n�mero_de_veces
                loop while condici�n
		loop for i=principio..final
Ejemplos:       loop 1000
		loop while essdiff > .00001
		loop for i=1991..2000

Esta instrucci�n (de gui�n) abre un modo especial en el cual el programa
acepta instrucciones a repetir o un n�mero de veces espec�fico, o
mientras una condici�n se satisfaga, o para los sucesivos valores
enteros de una variable �ndice (interna) i.

Dentro de un bucle s�lo se pueden utilizar 7 instrucciones:
genr, ols, print, sim, smpl, store y summary (store no puede usarse en
un bucle "while"). Con genr y ols es posible hacer bastantes cosas.
Se sale de este modo especial de introducir instrucciones de bucle
mediante la orden "endloop": en este momento se ejecutan las �rdenes
que est�n en la "pila".

Los bucles no pueden estar anidados. La instrucci�n ols muestra un
resultado especial dependiendo del tipo de bucle. Si se especifica un
"numero_de_veces" no se muestran los resultados de cada regresi�n
individual, sino que se obtiene una salida con (a) la media de
cada coeficiente estimado a lo largo de todas las repeticiones, (b) las
desviaciones t�picas de estos coeficientes estimados, (c) la media
de las desviaciones t�picas estimadas para cada coeficiente y (d) la
desviaci�n t�pica de las desviaciones t�picas estimadas. Esto s�lo
tiene sentido si hay alg�n input aleatorio en cada paso. La instrucci�n
est� dise�ada para el an�lisis de Monte Carlo.

Si se da una condici�n "while", se muestran los resultados del modelo
especificado a partir de la �ltima vuelta del bucle: esto est� dise�ado
para m�nimos cuadrados iterativos.

La instrucci�n "print" tambi�n se comporta de forma diferente en el
contexto de un bucle con "n�mero_de_veces". En concreto esta instrucci�n
muestra la media y la desviaci�n t�pica de la variable a lo largo de las
repeticiones del bucle. Esto est� dise�ado para variables que toman un
solo valor en cada iteraci�n, por ejemplo la scr (suma de cuadrados
residual $ess ) de una regresi�n. La instrucci�n "print" se comporta de
la forma usual con las otras construcciones de bucle.

La instrucci�n "store" (se usa s�lo una de ellas por bucle y s�lo en un
bucle con "n�mero_de_veces") escribe los valores de las variables
especificadas desde cada iteraci�n del bucle a un fichero. De esta
manera se guarda una grabaci�n completa de las variables. Este fichero
de datos puede despu�s ser le�do y analizado mediante gretl.

Ejemplo de c�digo de bucle (Monte Carlo):

   genr x = uniform()
   loop 100
   genr u = normal()
   genr y = (10*x) + (20*u)
   ols y const x
   genr r2 = $rsq
   print r2
   genr a = coeff(const)
   genr b = coeff(x)
   store foo.gdt a b
   endloop

#
lmtest
@Contrastes
Contraste de Multiplicador de Lagrange

Bajo este encabezamiento se encuentran varios contrastes de hip�tesis. 
Lo que tienen en com�n es que el contraste incluye la estimaci�n de una
regresi�n auxiliar, en la que la variable dependiente es el residuo de 
alguna regresi�n "original". Entre las variables del lado derecho
se incluyen las de la regresi�n original y algunas adicionales. El
estad�stico de contraste se calcula como (tama�o muestral x
R-cuadrado) de la regresi�n auxiliar: este se distribuye como una
Chi-cuadrado con grados de libertad iguales al n�mero de variables 
adicionales, bajo la hip�tesis nula de que las variables adicionales no 
tienen poder explicativo sobre el residuo. Un valor muy alto de este 
estad�stico (valor p peque�o) sugiere que esta hip�tesis nula deber�a 
ser rechazada.


#
markers
@Conjunto de datos
(marcadores)
A�ade marcadores de caja al conjunto de datos.

Esta instrucci�n necesita el nombre del fichero que contenga los
"marcadores de caja", es decir, peque�as etiquetas que
identifican a las observaciones individuales en el conjunto de datos
(por ejemplo, nombres o c�digos de pa�ses o de ciudades). Estas
etiquetas no deber�an tener m�s de 8 caracteres. El fichero deber�a
tener un marcador por l�nea y deber�a haber tantos marcadores como 
observaciones en el conjunto de datos. Si se satisfacen estas 
condiciones y se encuentra el fichero especificado, se a�adir�n los 
"marcadores de caja"; estos se podr�n ver cuando Vd elija "Mostrar
valores" en el men� "Datos" de gretl.


#
meantest
@Contrastes
(Contraste de medias)

Calcula el estad�stico t para el contraste de la hip�tesis nula de que
las medias poblacionales, de las dos variables elegidas, son iguales. 
Tambi�n muestra su valor p. La instrucci�n puede ejecutarse con o sin el 
supuesto de que las varianzas de las dos variables son iguales (aunque 
esto supondr� una diferencia en el estad�stico de contraste s�lo si hay 
un n�mero diferente de valores no-perdidos para las dos variables).


#
missing values
@Conjunto de datos
(Valores perdidos)

Establece un valor num�rico que ser� interpretado como "valor perdido" o
"no disponible", o para una serie de datos particular (bajo el
men� de "Variable") o globalmente para el conjunto de datos completo
(bajo el men� "Muestra").

Gretl tiene su propio c�digo interno para los valores perdidos, pero a
veces los datos importados pueden emplear un c�digo diferente. Por
ejemplo, si una serie determinada est� codificada de forma que el valor
-1 indica "no disponible", se puede seleccionar "Establecer c�digo de
'valor perdido'" bajo el men� de "Variable" y escribir el valor "-1"
(sin las comillas). Gretl entonces leer� los -1 como observaciones
perdidas.


#
mpols
@Estimaci�n
(MCO de precisi�n m�ltiple)

Calcula las estimaciones de m�nimos cuadrados ordinarios utilizando 
operaciones aritm�ticas con precisi�n m�ltiple. Esta orden s�lo est� 
disponible si gretl est� configurado con soporte para GMP, la biblioteca 
de precisi�n m�ltiple GNU. Hay que se�alar que la precisi�n de los 
resultados de la regresi�n puede verse limitada por (a) la precisi�n de 
los datos que se leen del fichero y (b) cualquier transformaci�n 
realizada usando la instrucci�n genr, que trabaja utilizando operaciones 
aritm�ticas ordinarias de punto flotante y doble precisi�n.

#
nulldata
@Conjunto de datos

Establece un conjunto de datos "vac�o", que contiene s�lo una constante, 
con periodicidad 1 y el n�mero de observaciones que se especifique. Esto 
puede utilizarse por motivo de simulaci�n: algunas instrucciones genr
(p.ej. genr uniform(), genr normal(), genr time) generar�n datos 
artificiales desde cero para rellenar el conjunto de datos. La 
instrucci�n "nulldata" tambi�n puede ser �til combinada con "loop".

#
ols
@Estimaci�n
M�todo de m�nimos cuadrados ordinarios

Calcula las estimaciones de m�nimos cuadrados ordinarios de los
coeficientes del modelo especificado. Muestra los valores p para los 
estad�sticos t (a dos colas) y F. Un valor p inferior a 0.01 indica 
significatividad al nivel del 1 por ciento. Tambi�n se muestran una 
serie de estad�sticos de selecci�n de modelos.

Ver "/Temas/Estimaci�n/Cuadro de di�logo" para ayuda sobre el uso del 
cuadro de di�logo.

#
omit
@Contrastes
Omite variables de un modelo y contrasta su significatividad conjunta

Las variables elegidas se sustraen del modelo anterior y se estima el 
nuevo modelo. Si se omite m�s de una variable, se mostrar� el 
estad�stico F de Wald para las variables omitidas junto con su valor p 
(s�lo para el m�todo MCO). Un valor p inferior a 0.05 indica que los 
coeficientes son conjuntamente significativos al nivel del 5 por ciento.


#
online databases
@Conjunto de datos
(Bases de datos en l�nea)
Acceso a bases de datos v�a internet

Gretl puede acceder a las bases de datos del sitio web de gretl, en
Wake Forest University (su ordenador debe de estar conectado a 
internet para que esto funcione).

Bajo el men� "Archivo, Revisar bases de datos" seleccione la entrada
"sobre servidor". Ahora deber�a aparecer una ventana mostrando las bases 
de datos gretl disponibles en Wake Forest (dependiendo del lugar en que
Vd se encuentre y de su conexi�n a internet, esto puede tardar
unos segundos). Junto al nombre de la base de datos y a una peque�a 
descripci�n aparecer� una entrada de "Estado local": esto indica si Vd 
ha instalado la base de datos de forma local (sobre el disco duro de su 
ordenador) y si es as�, si se encuentra actualizada con la versi�n
que actualmente hay en el servidor. Si Vd ha instalado localmente una
determinada base de datos y �sta se encuentra actualizada, no hay 
ninguna ventaja por acceder a ella mediante el servidor. Pero para una 
base de datos que no est� instalada y/o actualizada, Vd puede desear un 
listado de las series de datos: haga "click" sobre "Obtener listado de
series". Esto hace aparecer una nueva ventana desde la cual se pueden 
visualizar los valores de la serie de datos que se elija, representar 
esos valores o importarlos al espacio de trabajo de gretl. Estas tareas 
pueden realizarse usando el men� "Series" o por medio del men� 
contextual que aparece el hacer "click" con el bot�n derecho del rat�n 
sobre una serie dada. Tambi�n es posible buscar entre el listado una 
variable determinada (con el men� "Buscar").

Si se desea un acceso a los datos m�s r�pido, o un acceso a los datos 
"off line", es posible seleccionar la l�nea que muestra la base de datos 
que interese, en la ventana inicial de bases de datos, y presionar el 
bot�n "Instalar". Esto har� que la base de datos se descargue en formato 
comprimido, despu�s se puede descomprimir e instalar en el disco
duro. M�s adelante podremos encontrarla bajo el men� "Archivo, Revisar 
bases de datos, Nativa gretl". (Esta caracter�stica de gretl depende de
otros proyectos de software de 'fuente abierta': la biblioteca de 
compresi�n de datos zlib y el programa descargador GNU "wget", de los 
cuales gretl ha tomado prestados algunos trozos de c�digo).

#
panel
@Conjunto de datos
Establece estructura de datos de panel

Las dos opciones disponibles aqu� son "series temporales apiladas" y 
"secciones cruzadas apiladas". Si se quiere hacer uso de la
instrucci�n "MCO combinados" y sus diagn�sticos de panel asociados, 
gretl debe saber en qu� forma est�n organizados los datos. "Series 
temporales apiladas" significa que los bloques del fichero de datos son 
series temporales para cada una de las unidades de secci�n cruzada. Por 
ejemplo, las primeras 10 filas de datos podr�an representar los valores 
de ciertas variables para el pa�s A durante 10 periodos, las siguientes 
10 filas los valores para el pa�s B durante los mismos 10 periodos, y 
as� sucesivamente. "Secciones cruzadas apiladas" significa que los 
bloques del fichero de datos son secciones cruzadas para cada uno de los 
periodos. Por ejemplo, las primeras 6 filas de datos podr�an representar 
los valores de ciertas variables para los pa�ses A a F para el a�o 1970,
las siguientes 6 filas los valores para los mismos pa�ses en 1971, y as�
sucesivamente.

Si se guarda el fichero de datos despu�s de establecer este atributo, la
informaci�n se grabar� en el fichero de datos y no ser� necesario
establecerlo de nuevo la pr�xima vez que se usen estos datos.


#
pooled
@Estimaci�n
Estimaci�n de MCO combinados

Esta instrucci�n se utiliza con datos de panel. Para sacar provecho de
ella, se deber�a especificar un modelo sin ninguna variable ficticia
representando a las unidades de secci�n cruzada. La rutina presenta
estimaciones de MCO combinados directamente, las cuales tratan las
variaciones de secci�n cruzada y de series temporales de igual forma.
Este modelo puede ser o puede no ser apropiado. Bajo el men� de
"Contrastes" en la ventana del modelo se puede encontrar una entrada
"diagn�sticos de panel", en la cual se contrastan MCO combinados contra
las principales alternativas: los modelos de efectos fijos y de efectos
aleatorios.

En el modelo de efectos fijos se a�ade una variable ficticia para todas
las unidades de secci�n cruzada excepto una, permitiendo as� al t�rmino
constante de la regresi�n variar a trav�s de las unidades (individuos).
Se presenta un estad�stico F para contrastar la significatividad
conjunta de dichas variables ficticias: si el valor p de este contraste
es peque�o, esto es una indicaci�n en contra de la hip�tesis nula (de
que el modelo simple combinado es el adecuado) y en favor del modelo de
efectos fijos.

Por otro lado, el modelo de efectos aleatorios descompone la varianza
residual en dos partes, una parte espec�fica de la unidad de secci�n
cruzada o "grupo" y la otra espec�fica de la observaci�n particular.
(Este estimador s�lo puede calcularse si el panel es suficientemente
"ancho", es decir, si el n�mero de unidades de secci�n cruzada que hay
en el conjunto de datos es superior al n�mero de par�metros a estimar).
El estad�stico LM de Breusch-Pagan contrasta la hip�tesis nula (de
nuevo, de que el estimador de MCO combinados es el adecuado) contra la
alternativa de efectos aleatorios.

Es muy posible que el modelo de MCO Combinados sea rechazado contra
ambas alternativas (efectos fijos y efectos aleatorios). �C�mo se puede
entonces determinar cu�l de los dos estimadores es m�s apropiado? El
contraste de Hausman (que tambi�n se presenta, dado que se puede estimar
el modelo de efectos fijos) da una indicaci�n en este sentido. Si
que el error espec�fico de unidad --o grupo-- est� incorrelacionado con
las variables independientes, el estimador de efectos aleatorios es m�s
eficiente que el de efectos fijos; en caso contrario el estimador de
efectos aleatorios es inconsistente y entonces ser� preferible el
estimador de efectos fijos. La hip�tesis nula para el contraste de
Hausman es que el error espec�fico de grupo no est� muy correlacionado
(y por tanto, es preferible el estimador de efectos fijos). Entonces,
para este contraste, un valor p peque�o es una indicaci�n en contra del
modelo de efectos aleatorios y a favor del de efectos fijos.

Para un desarrollo riguroso de este tema ver  "An�lisis
Econom�trico" de Greene (4� edici�n), cap�tulo 14.

#
probit
@Estimaci�n
Regresi�n Probit

La variable dependiente deber�a ser una variable binaria. Se utiliza
el m�todo de m�nimos cuadrados iterativos (el m�todo EM o de
expectativa-maximizaci�n) para obtener las estimaciones
m�ximo-veros�miles de los coeficientes de las variables independientes.
Como el modelo es no lineal, las pendientes dependen de los valores de
las variables independientes: las pendientes que gretl muestra se
eval�an en las medias de dichas variables. El estad�stico Chi-cuadrado
contrasta la hip�tesis nula de que todos los coeficientes, excepto la
constante, son cero.

#
range-mean
@Gr�ficos
Gr�fico Rango-Media

Este es un gr�fico simple para ayudar a decidir si una serie temporal, y(t), 
tiene o no varianza constante. Se toma la muestra completa t=1,...,T y se divide 
en peque�as submuestras de tama�o arbitrario k [gretl elige k=sqrt(T)]. La 
primera submuestra se forma con y(1),...,y(k), la segunda con y(k+1),...,y(2k), 
y as� sucesivamente. Para cada submuestra se calcula la media muestral y el 
rango (=m�ximo-m�nimo) y se construye un gr�fico con las medias en el eje 
horizontal y los rangos en el vertical. De esta forma, cada submuestra est� 
representada por un punto en este plano. Si la varianza de la serie fuera 
constante los rangos de las submuestras no deber�a depender de sus medias; si se 
observa que los puntos se aproximan a una recta con pendiente creciente, esto 
sugiere que la varianza de la serie aumenta cuando la media aumenta; si los 
puntos se aproximan a una recta con pendiente decreciente, esto sugiere que la 
varianza est� disminuyendo cuando la media aumenta.

Adem�s del gr�fico, gretl presenta una ventana de resultados que muestra las 
medias y los rangos para cada submuestra, el coeficiente estimado para la
pendiente en una regresi�n MCO de los rangos sobre las medias y el valor p para 
el contraste de la hip�tesis nula de que esta pendiente es cero. Si el 
coeficiente de pendiente es significativo al nivel de significaci�n del 10 por 
ciento, en el gr�fico se muestra tambi�n la recta ajustada en la regresi�n de 
los rangos sobre las medias.

#
rhodiff
@Transformaciones
Uso:            rhodiff rho listavar
Ejemplo:        rhodiff .65 2 3 4

Crea las correspondientes variable rho-diferenciadas de las variables
(dadas por n�mero o nombre) de listavar y las a�ade al conjunto de
datos.  Sea la variable v1 de la lista entonces se crea
rd_v1 = v1(t) - rho*v1(t-1).

#
scatters
@Gr�ficos
M�ltiples gr�ficos cruzados por parejas

Representa un conjunto de gr�ficos cruzados de la "variable del eje Y"
elegida con respecto a las "variables del eje X" elegidas
consecutivamente. Puede ser �til echar un vistazo a estos gr�ficos al 
hacer un an�lisis exploratorio de datos. El n�mero m�ximo de gr�ficos es 
seis; cualquier variable extra en el eje X ser� ignorada.

#
seed
@Programaci�n
Pone en marcha el generador de n�meros aleatorios

Requiere como entrada un entero. Establece la semilla para el generador 
de n�meros pseudoaleatorios utilizado por las opciones "Aleatoria 
uniforme" y "Aleatoria normal" bajo el men� de "Datos, A�adir 
variables". Por defecto, la semilla se establece, utilizando la hora del 
sistema, cuando se inicia el programa. Si se desea obtener secuencias de
n�meros pseudoaleatorios repetibles es necesario establecer la semilla
de forma manual.


#
setobs
@Conjunto de datos
Establece la frecuencia de los datos y la observaci�n inicial

Utilice esta orden para forzar al programa a interpretar el conjunto de 
datos actual como "de series temporales" o "de panel" cuando los datos 
se han le�do inicialmente como series simples sin fecha. Se necesitan 
dos par�metros: una frecuencia entera y una observaci�n inicial 
(normalmente una fecha).

Ejemplos de entradas v�lidas:

  4 1990.1       Interpretar los datos como trimestrales, 
                 comenzando en 1990, trimestre 1.
  12 1978.03     Interpretar los datos como mensuales, comenzando en 
                 marzo de 1978
  20 1.01        Frecuencia de datos 20, comenzando en la observaci�n
                 1.01 (datos de panel)
  5 72/01/10     Datos diarios (5 d�as por semana), desde 10 de enero
                 de 1972
  7 02/01/10     Datos diarios (7 d�as por semana), desde 10 de enero
                 de 2002

#
sim
@Conjunto de datos
(Simulaci�n)
Introducir datos simulados en una variable

Esta instrucci�n requiere una observaci�n inicial, una observaci�n 
final, el nombre de una variable (ya existente en el conjunto de datos) 
en la cual introducir los valores y una lista de coeficientes
autorregresivos, que pueden ser constantes num�ricas o nombres 
de variables. Por ejemplo, si en el di�logo de simulaci�n se introduce
    
    1979.2 1983.1 y 0 0.9

esto rellenar� y, desde 1979.2 hasta 1983.1 con los valores

    y(t) = 0 + 0.9 y(t-1)

De forma similar

    15 25 y 10 0.8 x

generar� desde la observaci�n 15 a la 25:

    y(t) = 10 + 0.8 y(t-1) + x(t) y(t-2)

#
sampling
@Conjunto de datos
(Muestreo)
Seleccionar una submuestra del conjunto de datos actual.

Si se elige "Muestra/Definir a partir de v.ficticia..." es necesario 
suministrar el nombre de una variable ficticia (indicador) que 
deber�a tener los valores 0 o 1 en cada observaci�n. La muestra se 
restringir� a aquellas observaciones para las cuales la variable 
ficticia tome valor 1. (Haciendo "click" en la l�nea de una variable en 
la ventana principal de los datos se insertar� el nombre de esa variable 
en la caja de di�logo).

Si se elige "Muestra/Restringir a partir de criterio..." es necesario 
proporcionar una expresi�n booleana (l�gica), del mismo tipo que 
se utilizar�a para definir una variable ficticia. Por ejemplo, la 
expresi�n "sqft > 1400" seleccionar� s�lo los casos para los que la 
variable sqft tenga un valor mayor que 1400. Las condiciones pueden 
estar concatenadas utilizando los operadores l�gicos "&" (Y l�gico) y
"|" (O l�gico).

La entrada de men� "Muestra/Quitar todas las obs. con valores 
perdidos..." redefine la muestra excluyendo todas las observaciones para 
las que los valores de una o m�s variables est�n "perdidos" (dejando 
s�lo las observaciones completas, en las que hay un valor num�rico para 
todas las variables).

Al definir una muestra a partir de una variable ficticia, una expresi�n 
booleana o por el criterio de los valores perdidos hay que se�alar que
cualquier informaci�n "estructural" en el fichero de encabezamiento de 
datos (con relaci�n a la naturaleza de 'series temporales' o 'de panel' 
de los datos) se perder�. Es posible reimponer de nuevo la estructura 
mediante "Muestra/Establecer frecuencia, observaci�n inicial..."

Para s�lo reorganizar la muestra especificando una observaci�n inicial y 
una final ver "smpl".


#
smpl
@Conjunto de datos
(Establecer el rango muestral)

Restablece el rango muestral especificando una observaci�n inicial y
una observaci�n final (Muestra/Establecer rango...). Este mecanismo se
utiliza para establecer una submuestra de una serie temporal de datos.
Las observaciones inicial y final dadas deber�an estar en un formato
consistente con la frecuencia de los datos, p.ej. "1985.1" para datos
trimestrales o ""1996.03" para datos mensuales (marzo de 1996).

#
spearman
@Estad�sticos

Calcula el coeficiente de correlaci�n por rangos de Spearman para un par 
de variables especificado. No es necesario antes ordenar y hacer el
ranking de las variables, la funci�n se encarga de ello.

El ranking autom�tico es de mayor a menor (es decir, al dato mayor se le 
asigna rango 1). Si se necesita invertir este ranking, esto se puede 
hacer creando una nueva variable que sea el negativo de la original. Por 
ejemplo:

  genr altx = -x
  spearman altx y

#
square
@Transformaciones
(Cuadrados)

Genera nuevas variables que son los cuadrados de las variables de la
lista dada. Las nuevas variables se nombran con el prefijo "sq_", as�
por ejemplo, la nueva variable sq_x = x^2 (al cuadrado).

#
store
@Conjunto de datos
(Guardar los datos)

Guarda un conjunto de datos gretl. Hay dos opciones para el formato de
los datos guardados.

(1) "Formato Standard": los datos se guardan en el formato xml de gretl.
(2) "Comprimido gzip": como arriba, pero utilizando compresi�n de tipo
     gzip. Esto ahorra espacio de disco, puede ser �til para conjuntos
     de datos grandes.
     
N�tese que si Vd. desea guardar el valor de cualquier escalar generado 
en una sesi�n de gretl (en lugar de una serie de datos), Vd. deber�a usar
la instrucci�n "store" en la ventana de consola de gretl o en un gui�n
de gretl y especificar la lista de variables a guardar.

#
tsls
@Estimaci�n
(mc2e)
M�nimos cuadrados en dos etapas

Esta orden necesita la selecci�n de dos listas de variables: las 
variables independientes que aparecer�n en el modelo dado y un conjunto 
de "instrumentos". Este �ltimo comprende las variables ex�genas y/o 
predeterminadas que pueden usarse como regresores para obtener valores 
estimados de las variables end�genas del lado derecho e la ecuaci�n del 
modelo. Si alguna de las variables del lado derecho del modelo son
ex�genas, �stas deber�an aparecer en ambas listas.


#
var
@Estimaci�n
Autorregresi�n vectorial

Esta instrucci�n necesita un orden de retardos entero. Hay que 
seleccionar las variables dependiente e independientes para la primera 
ecuaci�n del sistema; estas se permutar�n para obtener las
ecuaciones restantes. NO INCLUIR NINGUNA VARIABLE RETARDADA en la lista 
de variables independientes --se a�adir�n autom�ticamente.

En general, se realizar� una regresi�n para cada variable de la lista, 
excepto la constante, la tendencia temporal y cualquier variable 
ficticia. La salida de cada ecuaci�n incluye los contrastes F para 
restricciones cero sobre todos los retardos de cada una de las 
variables y un contraste F para el m�ximo retardo.

#
vartest
@Contrastes
(Contraste de varianzas)

Calcula el estad�stico F para la hip�tesis nula de que las varianzas
poblacionales para las dos variables elegidas son iguales y muestra su
valor p.

#
wls
@Estimaci�n
(mcp)
M�todo de m�nimos cuadrados ponderados

Sea "varpond" la variable elegida en la caja de "variable de
ponderaciones". Se ejecuta una regresi�n MCO, donde la variable
dependiente es el producto de varpond por la variable dependiente
elegida y las variables independientes  tambi�n se multiplican por
varpond. Si varpond es una variable ficticia, esto es equivalente a 
eliminar todas las observaciones que tengan el valor cero para varpond.

