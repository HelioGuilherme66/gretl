#
add
@Tests
Uso:            add listavar
Ejemplos:       add 5 7 9	          add xx yy zz

Se a�aden las variables de la lista 'listavar' al modelo anterior y se
estima el nuevo modelo. Si se a�ade m�s de una variable, se presenta el
estad�stico F para las variables a�adidas (s�lo para el m�todo MCO) junto
con su valor-p. Un valor-p inferior a 0.05 indica que los coeficientes son
conjuntamente significativos al nivel del 5 por ciento.

#
addto
@Tests
Uso:            addto model_ID listavar
Ejemplo:        addto 2 5 7 9

Funciona como la instrucci�n "add", pero aqu� hay que especificar un modelo
anterior (lo cual se hace utilizando su n�mero de ID, que se presenta al
principio de los resultados del modelo) para tomarlo como base al a�adir las
variables. En el ejemplo de arriba se a�aden las variables con n�meros 5, 7
y 9 al modelo 2.

#
adf
@Tests
Uso:	        adf orden nombrevar
Ejemplo:        adf 2 x1

Calcula los estad�sticos para dos contrastes de Dickey-Fuller. En ambos casos
la hip�tesis nula es que la variable en cuesti�n presenta una ra�z unitaria.

El primero es un contraste t basado en el modelo
(1 - L)x(t) = m + g * x(t-1) + e(t).
La hip�tesis nula es que g = 0.

En el segundo contraste (aumentado) se estima una regresi�n no restringida
(teniendo como regresores una constante, una tendencia temporal, el primer
retardo de la variable y "orden" retardos de la primera diferencia) y una
regresi�n restringida (quitando la tendencia temporal y el primer retardo).
El estad�stico de contraste es F calculado como

[(SCRr - SCRnr)/2]/[SCRnr/(T - k)]

donde T es el tama�o muestral y k el n�mero de par�metros del modelo no
restringido. Es necesario tener en cuenta que los valores cr�ticos para
estos estad�sticos no son los usuales.

#
ar
@Estimation
Uso:            ar retardos ; vardep varindep      o
                ar retardos ; -o vardep varindep
Ejemplo:        ar 1 3 4 ; y 0 x1 x2 x3

Calcula las estimaciones de un modelo utilizando el procedimiento iterativo
generalizado de Cochrane-Orcutt (ver Ramanathan, secci�n 9.5). Las iteraciones
se terminan cuando las sucesivas sumas de cuadrados residuales no var�an m�s del
0.005 por ciento o cuando han transcurrido 20 iteraciones. 'retardos' es una
lista de retardos de los residuos, que acaba con un punto y coma. En el
ejemplo de arriba el t�rmino de error se especifica como

u(t) = rho1 u(t-1) + rho3 u(t-3) + rho4 u(t-4) + et

'vardep' es la variable dependiente y 'varindep' es la lista de variables
independientes que van separadas con espacios. Utilice el n�mero 0 para un
t�rmino constante. Si se usa la opci�n -o se presentar� la matriz de varianzas
y covarianzas de los coeficientes de regresi�n. Los residuos de la regresi�n
transformada se guardan bajo el nombre 'uhat' y pueden recuperarse
usando la instrucci�n 'genr'.

#
arch
@Tests
Uso:            arch retardo vardep varindep
Ejemplos:       arch 4 1 0 2 4 6 7      o  arch 4 y 0 x1 x2 x3 x4

Esta instrucci�n contrasta la posibilidad de un ARCH en el modelo, del orden
especificado en "retardo" (que debe de ser entero). Si el estad�stico de
contraste LM tiene un valor-p inferior a 0.10 tambi�n se realiza la estimaci�n
ARCH. Si la varianza predicha de cualquier observaci�n en la regresi�n auxiliar
no es positiva, se utiliza en su lugar el correspondiente 'uhat' (residuo)
al cuadrado. Luego, se estima por m�nimos cuadrados ponderados el modelo 
original.


#
chow
@Tests
Uso:           chow obs
Ejemplos:      chow 25
               chow 1988.1

Primero debe ejecutarse una regresi�n por MCO. Crea una variable ficticia
que es igual a 1 desde el punto de corte especificado en "obs" hasta el final
de la muestra y 0 en el resto. Tambi�n crea los t�rminos de interacci�n entre
esta variable ficticia y las variables independientes originales. Se ejecuta
una regresi�n aumentada incluyendo estos t�rminos y se calcula un estad�stico F,
tomando la regresi�n aumentada como 'no restringida' y la original como
'restringida'. Este estad�stico es adecuado para contrastar la hip�tesis nula de
que no hay cambio estructural en el punto de ruptura indicado.


#
coint
@Tests
Uso:	        coint orden vardep varindep
Ejemplos:       coint 2 y x
                coint 4 y x1 x2

Realiza los contrastes de Dickey-Fuller para cada una de las variables listadas.
Para cada variable, considera la hip�tesis nula de que la variable tiene una 
ra�z
unitaria y utiliza para el contraste el orden de retardos dado.
Se estima la regresi�n cointegrante y se realiza un contraste ADF sobre
los residuos de esta regresi�n. Tambi�n se proporciona el estad�stico de 
Durbin-Watson para la regresi�n cointegrante.
	Hay que se�alar que para ninguno de estos estad�sticos de contraste
se pueden aplicar las tablas estad�sticas usuales.

#
corc
@Estimation
Uso:          corc vardep varindep       o    corc -o vardep varindep
Ejemplos:       corc 1 0 2 4 6 7                corc -o 1 0 2 4 6 7
                corc y 0 x1 x2 x3               corc -o y 0 x1 x2 x3

Calcula las estimaciones de un modelo utilizando el procedimiento iterativo de
Cochrane-Orcutt (ver Ramanathan, Secci�n 9.4) siendo 'vardep' la variable
dependiente y siendo 'varindep' una lista de variables independientes separadas
por espacios y acabando con un ;. Utilice el n�mero 0 para un t�rmino constante. 




El proceso iterativo se detiene cuando valores sucesivos de rho no difieren en 
m�s
de 0.001 o cuando han transcurrido 20 iteraciones. Si se utiliza la opci�n -o, 
se muestra la matriz de covarianzas de los coeficientes de regresi�n. La 
regresi�n transformada final se calcula para el rango de observaci�n
primobs+1 ultobs que est� actualmente en efecto. Los residuos de esta regresi�n
transformada se guardan con el nombre 'uhat'.


#
corr
@Statistics
Uso:         corr
             corr listavar

'corr' muestra los coeficientes de correlaci�n para todos los pares de variables
que hay en el conjunto de datos (los valores perdidos, que se denotan por -999,
no se tienen en cuenta). 'corr listavar' muestra los coeficientes de correlaci�n
para las variables listadas.


#
corrgm
@Statistics
Uso:          corrgm nombrevar o numerovar
              corrgm nombrevar o numerovar maxretardo

Muestra los valores de la funci�n de autocorrelaci�n para la variable
especificada (ver Ramanathan, secci�n 11.7). Es, por tanto, corr[u(t), u(t-s)],
donde u(t) es la observaci�n t-�sima de la variable u y s es el n�mero
de retardos. Tambi�n se muestran las correlaciones parciales: estas tienen
descontado el efecto de los retardos que intervienen. La instrucci�n tambi�n
representa el correlograma y calcula el estad�stico Q de Box-Pierce para
contrastar la hip�tesis nula de que la serie es 'ruido blanco'. Este se
distribuye asint�ticamente como una chi-cuadrado con un n�mero de grados
de libertad igual al n�mero de retardos utilizados.

Si se suministra un n�mero entero 'maxretardo' la largura del correlograma
se limita a, como m�ximo, ese n�mero de retardos, en caso contrario la largura
se determina autom�ticamente.


#
criteria
@Utilities
Uso:          criteria scr T k        p.ej. criteria 23.45 45 8

Dados scr (suma de cuadrados de los residuos), el n�mero de observaciones (T) 
y el n�mero de coeficientes (k), calcula los estad�sticos de selecci�n de
modelos (ver Ramanathan, Secci�n 4.4). T, k y scr pueden ser valores num�ricos
o los nombres de variables definidas previamente.

#
critical
@Utilities
Uso:            critical t gl           p.ej. critical t 20
                critical X gl           p.ej. critical X 5
                critical F gln gld      p.ej. critical F 3 37
                critical d n            p.ej. critical d 24

Si el primer par�metro es t, X o F, muestra los valores cr�ticos para la 
distribuci�n t de student, chi-cuadrado o F respectivamente, para los niveles
de significaci�n m�s comunes, utilizando los grados de libertad especificados.
Si el primer par�metro es d, muestra los valores superior e inferior del 
estad�stico de Durbin-Watson al nivel de significaci�n del 5 por ciento para el 
valor de n (n�mero de observaciones) dado y para el rango de 1 a 5 variables 
explicativas.


#
cusum
@Tests
Uso:          cusum

Debe ejecutarse despu�s de una estimaci�n MCO. Desarrolla el contraste 
CUSUM de estabilidad de los par�metros. Se obtiene una serie de errores de
predicci�n (escalados) un paso hacia adelante, ejecutando para ello una serie
de regresiones: la primera regresi�n se realiza con las primeras k
observaciones (donde k es el n�mero de par�metros en el modelo original)
y se usa para generar una predicci�n de la variable dependiente para la
observaci�n k+1; la segunda utiliza las primeras k+1 observaciones y genera
una predicci�n para la observaci�n k+2 y as� sucesivamente. Se muestra y 
representa gr�ficamente la suma acumulada de los errores de predicci�n 
escalados. Se rechaza, al nivel de significaci�n del 5 por ciento, la 
hip�tesis nula de estabilidad de los par�metros si la suma acumulada 
sale fuera de las bandas del 95 por ciento de confianza.

Tambi�n se proporciona el estad�stico t de Harvey-Collier para contrastar 
la hip�tesis nula de estabilidad de los par�metros. Ver Cap�tulo 7 del libro 
de Greene 'Econometric Analysis' para m�s detalles.

#
delete
@Dataset
Uso:          delete

Elimina la �ltima (la de n�mero m�s alto) variable del conjunto de datos actual.
Util�cela con precauci�n: no se pide confirmaci�n. Puede ser �til para eliminar 
variables ficticias temporales. S�lo se puede eliminar la �ltima variable.

#
diff
@Transformations
Uso:          diff listavar

Se toma la primera diferencia de cada variable en 'listavar' y el resultado se 
guarda en una nueva variable con el prefijo "d_". As�, "diff x y" crea las 
nuevas variables d_x=x(t)-x(t-1) y d_y = y(t) - y(t-1).

#
endloop
@Programming

Termina un bucle de simulaci�n. Ver "loop".

#
eqnprint
@Printing
Uso:            eqnprint
                eqnprint -o

Debe ejecutarse despu�s de una estimaci�n por MCO. Imprime el modelo estimado,
en forma de ecuaci�n LaTeX, a un fichero cuyo nombre tiene la estructura
"equation_N.tex", donde N es el n�mero de modelos estimados hasta el momento
en la sesi�n actual. Este puede incorporarse en un documento LaTeX. Ver
tambi�n la instrucci�n 'tabprint'.

Si se utiliza la opci�n -o, el fichero LaTeX es un documento completo
(con pre�mbulo LaTeX) listo para procesar; en caso contrario debe incluirse
en un documento (que ya tenga el pre�mbulo).


#
fcast
@Prediction
Uso:            fcast primobs ultobs nuevonombrevar
                fcast nuevonombrevar
Ejemplo:        fcast ajustados

Los valores ajustados en la �ltima regresi�n ejecutada se guardan bajo 
'nuevonombrevar'. Estos valores pueden mostrarse e representarse 
gr�ficamente. Las variables del lado derecho de la ecuaci�n son las del modelo
original. No hay posibilidad de cambiar a otras variables. Si se especifican 
las observaciones inicial ('primobs') y final ('ultobs') la predicci�n se
restringe al rango especificado. Si se ha especificado un t�rmino de error 
autorregresivo (en 'hilu','corc' y 'ar'), la predicci�n es condicional un paso 
adelante e incorpora el proceso del error.


#
fcasterr
@Prediction
Uso:            fcasterr primobs ultobs
                fcasterr primobs ultobs -o

Despu�s de estimar un modelo MCO que incluya una constante y al menos una 
variable independiente (estas restricciones pueden relajarse en alg�n punto),
se puede usar esta instrucci�n para mostrar los valores ajustados sobre el 
rango de observaci�n especificado, junto con las desviaciones t�picas estimadas
de estas predicciones y los intervalos de 95 por ciento de confianza. Si se 
utiliza la opci�n -o tambi�n se mostrar�n los resultados en gr�ficos gnuplot.


#
fit
@Prediction
Uso:		fit

Esta orden (que debe seguir a una instrucci�n de estimaci�n) es una atajo para
la instrucci�n 'fcast'. Genera valores ajustados para la muestra actual, 
basados en la �ltima regresi�n y los guarda en una serie denominada "autofit".
En modelos de series temporales tambi�n muestra un gr�fico gnuplot de los 
valores
actual y estimado de la variable dependiente contra el tiempo.

#
freq
@Statistics
Uso:          freq nombrevar (o numerovar)

Muestra la distribuci�n de frecuencias de 'nombrevar' o 'numerovar';
tambi�n se proporcionan los resultados de una contraste chi-cuadrado de 
normalidad. El estad�stico para este �ltimo es:

  tama�o_muestral * [asimetr�a^2/6 + (CURTOSIS - 3.0)^2/24.0]

Bajo la hip�tesis nula de normalidad se distribuye como una chi-cuadrado
con 2 grados de libertad.

En modo interactivo se genera un gr�fico gnuplot de la distribuci�n.


#
genr
@Dataset
Uso:          genr nuevonombrevar = formula

Crea nuevas variables, normalmente por medio de transformaciones de
variables ya existentes. Ver tambi�n 'diff', 'logs', 'lags', 'ldiff',
'multiply' y 'square' como atajos.

Las operaciones aritm�ticas permitidas son, en orden de precedencia:
^(exponenciaci�n); *, / y % (m�dulo o resto); + y -.

Los operadores booleanos (de nuevo en orden de precedencia) son:
! (NO l�gico [NOT]), & (Y l�gico [AND]), | (O l�gico [OR]), >, <, =
y los s�mbolos compuestos != (no igual), >= (mayor o igual que) y
<= (menor o igual que).  Los operadores booleanos pueden usarse al
definir variables ficticias: por ejemplo (x > 10)
produce 1 si x(t) > 10, 0 en los dem�s casos.

Las funciones permitidas pertenecen a estos grupos:

-Funciones matem�ticas standard: abs, cos, exp, int (parte entera), ln
(logaritmo natural: log es un sin�nimo), sin (seno), sqrt (ra�z cuadrada).

-Funciones estad�sticas: mean (media aritm�tica), median (mediana),
var (varianza), sd(desviaci�n t�pica o standard), sum, cov (covarianza),
corr (coeficiente de correlaci�n), min (m�nimo) y max (m�ximo).

-Funciones de series temporales: lag (retardo), lead (adelanto),
diff (primera diferencia), ldiff (log-diferencia, o primera diferencia del
logaritmo natural).

-Miscel�neas: cum (acumulaci�n), sort (ordenar), uniform, normal,
misszero (reemplazar los c�digos de 'observaci�n perdida'  por ceros),
zeromiss (operaci�n inversa de misszero), pvalue (valor de probabilidad
para un estad�stico dado contra una distribuci�n especificada)
y mpow (elevar una serie a un exponente entero utilizando aritm�tica de 
precisi�n m�ltiple).

Todas las funciones anteriores, con las excepciones de cov, corr, uniform, 
normal, pvalue y mpow toman como �nico argumento o el nombre de una variable
(hay que tener en cuenta que en 'genr' no podemos referirnos a una variable por
su n�mero de ID) o una expresi�n compuesta que se eval�a en una variable (p.ej.
ln((x1+x2)/2)). 'cov' y 'corr' necesitan dos argumentos y producen
respectivamente la covarianza y el coeficiente de correlaci�n entre dos
variables especificadas. uniform() y normal(), que no tienen argumentos,
producen series pseudo-aleatorias obtenidas a partir de la distribuci�n uniforme
(0-100) y la distribuci�n normal standard respectivamente (ver tambi�n la
instrucci�n seed). La funci�n pvalue() toma los mismos argumentos que la orden
pvalue (ver m�s abajo), pero en este contexto deben situarse comas entre sus
argumentos. La funci�n mpow toma como argumentos el nombre de una serie de datos 
y un n�mero entero positivo, que es exponente al cual se desea elevar la serie.

Adem�s de los operadores y funciones mencionados hay algunos usos especiales
de 'genr':

* genr time crea una variable de tendencia temporal (1,2,3,...) denominada time.
* genr index crea una variable �ndice (1,2,3,...) denominada index.
* genr dummy crea variables ficticias hasta la periodicidad de los datos.
  P.ej. en el caso de datos trimestrales (periodicidad 4), el programa crea
  dummy_1 = 1 para el primer trimestre y 0 en los otros trimestres, dummy_2 = 1
  para el segundo trimestre y 0 en los otros trimestres, y as� sucesivamente.
* genr paneldum crea un conjunto de variables ficticias especiales para el uso
  con un panel de datos (ver el manual de gretl para m�s detalle)
* Pueden recuperarse usando 'genr' varias variables internas que se definen
  mientras se ejecuta una regresi�n. Esto se hace de la siguiente forma:

  $ess         suma de cuadrados de los residuos
  $rsq         R-cuadrado no corregido
  $T           n�mero de observaciones utilizadas en el modelo
  $df          grados de libertad
  $trsq        TR^2 (tama�o muestral por el R-cuadrado)
  $sigma       desviaci�n t�pica de los residuos
  $lnl         log-verosimilitud (en modelos logit y probit)
  coeff(var)   coeficiente estimado para var
  stderr(var)  desviaci�n t�pica estimada para var
  rho(i)       coeficiente de autorregresi�n de orden i-�simo de los residuos
  vcv(xi,xj)   covarianza entre los coeficientes de las variables xi y xj

La variable interna $nobs contiene el n�mero de observaciones en el rango
muestral actual, que puede ser o puede no ser igual al $T del �ltimo modelo.

La variable interna $pd contiene la periodicidad o frecuencia de los datos
(por ejemplo, 4 para datos trimestrales, 12 para mensuales).

La variable interna t sirve para referirse a las observaciones, comenzando
en 1. As� uno puede hacer "genr dum15 = (t=15)" para generar una variable
ficticia con valor 1 para la observaci�n 15 y 0 en el resto.

Ejemplos de instrucciones 'genr':

  genr y = x1^3          [x1 al cubo]
  genr y=ln((x1+x2)/x3)  [argumento compuesto para una funci�n ln]
  genr z=x>y             [asigna z(t) a 1 si x(t) > y(t), en otro caso a 0]
  genr y=x(-2)           [x retardada 2 periodos]
  genr y=x(2)            [x adelantada 2 periodos]
  genr y = mean(x)       [media aritm�tica]
  genr y = diff(x)       [y(t) = x(t) - x(t-1)]
  genr y = ldiff(x)      [y = ln(x(t)) - ln(x(t-1))]
                          ldiff(x) es la tasa de crecimiento instant�nea de x.
  genr y = sort(x)       [ordena x en orden creciente y lo guarda en y]
  genr y = - sort(-x)    [ordena x en orden decreciente]
  genr y = int(x)        [trunca x y guarda su valor entero como y]
  genr y = abs(x)        [guarda los valores absolutos de x]
  genr y = sum(x)        [suma los valores de x excluyendo los valores perdidos 
-999]
  genr y = cum(x)        [acumula x: y(t) es la suma de x hasta t]
  genr aa = $ess         [aa = suma de cuadrados de los residuos de la �ltima 
regresi�n]
  genr x = coeff(sqft)   [guarda en x el coeficiente de la variable sqft 
obtenido 
                          en el �ltimo modelo]
  genr rho4 = rho(4)     [guarda en rho4 el coeficiente autorregresivo de cuarto 



orden
                          obtenido del �ltimo modelo (supone un modelo ar)]
  genr cv=vcv(x1, x2)    [covarianza entre los coeficientes de x1 y x2 en el 
�ltimo modelo]
  genr x=uniform()/100   [variable pseudo-aleatoria uniforme, de rango 0 a 1]
  genr x=3*normal()      [variable pseudo-aleatoria normal, de media 0 y desv. 
t�pica 3]
  genr x=pvalue(t,20,1.4)[valor p para 1.4, bajo la distribuci�n t con 20 grados 



de libertad]

Sugerencias sobre variables ficticias:

* Supongamos que x se codifica con los valores 1, 2, o 3 y Vd desea tres 
variables ficticias d1 = 1 si x = 1, 0 en otro caso, d2 = 1 si x = 2, y as� 
sucesivamente. Para crear estas, utilice las instrucciones 
genr d1 = (x=1), genr d2 = (x=2), y genr d3 = (x=3).

* Para obtener la serie z = m�x(x,y) haga genr d=x>y y genr z=(x*d)+(y*(1-d))

#
gnuplot
@Graphs
Uso:            gnuplot yvar1 xvar [ opci�n ]
                gnuplot yvar1 yvar2 xvar [ opci�n ]
		gnuplot yvar xvar ficticia -z

En los dos primeros casos las variables yvars se representan contra xvar.
Si se proporciona la opci�n -o el gr�fico utilizar� l�neas; si se da la
opci�n -m el gr�fico utiliza impulsos (l�neas verticales); en los dem�s
casos se usan puntos.

En el tercer caso, yvar se representa contra xvar mostrando los puntos en
diferentes colores dependiendo de si el valor de 'ficticia' es 1 o 0.

Para crear un gr�fico de serie temporal, pedir "gnuplot yvars time".
Si no existe la variable "time", se generar� autom�ticamente.  Se crear�n
variables ficticias especiales para representar datos trimestrales y mensuales.

En modo interactivo, el resultado se pasa a gnuplot para que lo muestre
en pantalla.  En modo 'batch' se graba un fichero de nombre gpttmp<n>.plt,
donde <n> es un n�mero entre 1 y 99. M�s tarde, pueden generarse los gr�ficos
usando la instrucci�n de consola "gnuplot gpttmp<n>.plt".

#
graph
@Graphs
Uso:            graph var1 var2
                graph -o var1 var2
                graph var1 var2 var3

En los dos primeros ejemplos, la variable var1 (que puede ser un nombre o un
n�mero) se representa (eje y) contra var2 (eje x). La opci�n -o har� el gr�fico
con 40 filas y 60 columnas, sin ella el gr�fico ser� de 20 por 60 (salida de
pantalla). En el tercer ejemplo, las dos, var1 y var2 se representar�n (sobre
el eje y) contra var3. Esto es especialmente �til para representar los valores
observados y predichos contra el tiempo.

#
hausman
@Tests
Uso:          hausman

Este contraste s�lo est� disponible despu�s de haber estimado un modelo 
utilizando la orden "pooled" (ver tambi�n las instrucciones "panel" y 
"setobs"). Contrasta el modelo combinado simple contra las alternativas 
principales, los modelos de efectos fijos y de efectos aleatorios.

En el modelo de efectos fijos se a�ade una variable ficticia para todas las
unidades de secci�n cruzada menos una, permitiendo al t�rmino constante de la
regresi�n que var�e a trav�s de las unidades. En el modelo de efectos 
aleatorios se descompone la varianza residual en dos partes, una parte 
espec�fica de la unidad de secci�n cruzada y la otra espec�fica de la 
observaci�n particular. (Este estimador s�lo puede calcularse si el n�mero 
de unidades de secci�n cruzada en el conjunto de datos es mayor que el n�mero
de par�metros a estimar. Se presenta el estad�stico LM de Breusch-Pagan para 
contrastar la hip�tesis nula (de que el estimador MCO combinado es adecuado)
contra la alternativa de efectos aleatorios.

El modelo de MCO combinados puede ser rechazado contra ambas alternativas, 
efectos fijos y  efectos aleatorios. Si el error espec�fico de grupo 
-o de unidad- est� incorrelacionado con las variables independientes, el 
estimador de efectos aleatorios es m�s eficiente que el estimador de efectos 
fijos; en caso contrario, el estimador de efectos aleatorios es inconsistente
y ser� preferible utilizar el estimador de efectos fijos. La hip�tesis nula 
para el contraste de Hausman es que el error espec�fico de grupo no est� 
correlacionado con las variables explicativas (y por tanto, que es preferible 
el modelo de efectos aleatorios). Un valor p bajo para este contraste es 
una indicaci�n en contra del modelo de efectos aleatorios y a favor del de 
efectos fijos.


#
hccm
@Estimation
Uso:           hccm vardep varindep
            o  hccm -o vardep varindep

(Heteroskedasticity Consistent Covariance Matrix)
Presenta las estimaciones de MCO con las desviaciones t�picas de los
coeficientes obtenidas por medio de una estimaci�n de la matriz de varianzas
y covarianzas consistente ante heterocedasticidad. Utiliza para ello el
m�todo "jacknife" de MacKinnon-White.


#
help
@Utilities
help              proporciona una lista de instrucciones gretl
help nombreinst   describe la instrucci�n 'nombreinst' (p.ej. help smpl)

#
hilu
@Estimation
Uso:            hilu vardep varindep       o    hilu -o vardep varindep
Ejemplos:       hilu 1 0 2 4 6 7                hilu -o 1 0 2 4 6 7
                hilu y 0 x1 x2 x3               hilu -o y 0 x1 x2 x3

Calcula las estimaciones de un modelo utilizando el procedimiento de
b�squeda de Hildreth-Lu (se hace el ajuste fino usando el m�todo de
Cochrane-Orcutt) siendo 'vardep' la variable dependiente y 'varindep'
una lista de variables independientes separadas por espacios. Utilice el
n�mero 0 para incluir un t�rmino constante. Se representa gr�ficamente la
suma de cuadrados de los residuos del modelo transformado contra los
valores de rho desde -0.99 hasta 0.99. Si se usa la opci�n -o se mostrar�
tambi�n la matriz de varianzas y covarianzas de los coeficientes. Finalmente,
la �ltima regresi�n transformada se estima para el rango de observaci�n
primobs+1 ultobs que est� en efecto. Los residuos de esta regresi�n
transformada se guardan bajo el nombre 'uhat'.


#
hsk
@Estimation
Uso:            hsk vardep varindep
            o   hsk -o vardep varindep

Calcula estimaciones corregidas de heterocedasticidad y sus estad�sticos 
asociados. Se ajusta una regresi�n auxiliar para el logaritmo de los 
cuadrados de los residuos (utilizando los cuadrados de las variables 
independientes, pero no sus productos cruzados) y a partir de esta 
estimaci�n se obtienen los estimadores de m�nimos cuadrados ponderados del
modelo inicial. Si se usa la opci�n -o, se mostrar� tambi�n la matriz de
varianzas y covarianzas estimada de los coeficientes de la regresi�n.


#
if
@Programming
Uso:            if condici�n_boolena
                  instrucci�n1
                  instrucci�n2 ...
                endif

Las instrucciones gretl que est�n dentro del bloque "if ... endif" se 
ejecutan si y s�lo si la condici�n booleana se eval�a como cierta (no cero).
Para conocer la sintaxis de las condiciones booleanas en gretl, ver la 
instrucci�n 'genr'. Opcionalmente, la orden "endif" puede ir precedida de
una orden "else" (en una l�nea aparte para ella sola, como las instrucciones 
"if" y "endif"), seguida de un bloque de instrucciones a ejecutar si
la condici�n booleana original se eval�a como falsa (cero). Los bloques
entre "if", "else" y "endif" pueden contener tantas �rdenes como Vd
quiera y estas condiciones pueden estar anidadas. Cada orden "if" debe 
estar emparejada con una orden "endif".


#
import
@Dataset
Uso:            import archivo_csv
                import -o archivo_box

Sin la opci�n -o, importa datos desde un archivo que tenga formato de
'valores separados por comas' (CSV), como por ejemplo los que se pueden
escribir f�cilmente desde cualquier programa de hoja de c�lculo. El
archivo deber�a tener en la primera l�nea nombres de variables y, en
el resto, una matriz de datos rectangular. Las variables deber�an estar
alineadas "por observaci�n" (una columna para cada variable; cada fila
representa una observaci�n).

Con la opci�n -o, lee un fichero de datos en formato BOX1, como los que
se obtienen utilizando el servicio de extracci�n de datos del 'US Bureau
of the Census'.

#
info
@Dataset
Uso:          info

Muestra la informaci�n contenida en el archivo de cabecera correspondiente
al fichero de datos actual. Esta informaci�n debe estar limitada entre los 
caracteres "(*" y "*)", estando situados estos marcadores en l�neas separadas.

#
labels
@Dataset
Uso:          labels

Muestra las etiquetas informativas de las variables que se hayan definido
utilizando la instrucci�n 'genr'.

#
lad
@Estimaci�n
(Estimador de m�nima desviaci�n absoluta)
Uso:          lad vardep varindeps

Calcula una regresi�n que minimiza la suma de las desviaciones absolutas entre 
los valores observados y ajustados de la variable dependiente. Las estimaciones 
de los coeficientes se obtienen utilizando el algoritmo simplex de 
Barrodale-Roberts; surge un mensaje de aviso si la soluci�n no es �nica. Las 
desviaciones t�picas se obtienen por un m�todo 'bootstrap' con 500 
extracciones.

#
lags
@Transformations
Uso:          lags listavar

Crea variables nuevas que son valores retardados de cada una de
las variables que haya en 'listavar'. El n�mero de variables retardadas
que se crean es igual a la periodicidad. Por ejemplo, si la periodicidad
fuera 4 (datos trimestrales), la orden 'lags x y' crear� x_1  = x(-1),
x_2 = x(-2), x_3 = x(-3) y x_4 = x(-4); y de igual forma para y. Estas
variables deben referenciarse de forma exacta, es decir, con el car�cter de
subrayado.


#
ldiff
@Transformations
Uso:          ldiff listavar

Se calcula la primera diferencia del logaritmo natural de cada 
variable de 'listavar' y el resultado se guarda en una variable nueva
que lleva el prefijo "ld_". As� por ejemplo, "ldiff x y" crea las 
nuevas variables

               ld_x = ln[x(t)] - ln[x(t-1)]
               ld_y = ln[y(t)] - ln[y(t-1)].

#
lmtest
@Tests
Uso:            lmtest
                lmtest -o

Debe utilizarse justo despu�s de una instrucci�n 'ols'. Calcula el 
estad�stico de contraste del multiplicador de Lagrange (LM) para las 
hip�tesis alternativas de no linealidad y de heterocedasticidad 
(Contraste de White) o, si se utiliza la opci�n -o, para correlaci�n de 
orden hasta la periodicidad. Tambi�n se muestran los coeficientes
de la regresi�n auxiliar correspondiente. (Ver cap�tulos 7, 8 y 9 del
libro de Ramanathan para m�s detalles).

S�lo se usan los cuadrados de las variables independientes, y no sus
productos cruzados. No se pueden obtener los estad�sticos de
contraste LM si la generaci�n interna de los cuadrados causa
multicolinealidad exacta.


#
logit
@Estimation
Uso:          logit vardep varindeps

Regresi�n Logit: la variable dependiente deber�a ser una variable binaria.
Se obtienen las estimaciones por m�xima verosimilitud de los coeficientes de 
'varindeps' por medio de m�nimos cuadrados iterados (M�todo EM, o de 
expectativa-maximizaci�n). Como el modelo no es lineal, las pendientes dependen
de los valores de las variables independientes: las pendientes que se muestran 
se eval�an en la media de dichas variables. El estad�stico Chi-cuadrado 
contrasta la hip�tesis nula de que todos los coeficientes, excepto la 
constante, son cero.


#
logs
@Transformations
Uso:          logs listavar

Se obtiene el logaritmo natural de cada variable en 'listavar' y el
resultado se guarda en una nueva variable con prefijo "l_". As� por ejemplo,
"logs x y" crea las nuevas variables l_x = ln(x) y l_y = ln(y).


#
loop
@Programming
Uso:            loop n�mero_de_veces
                loop while condici�n
		loop for i=principio..final
Ejemplos:       loop 1000
		loop while essdiff > .00001
		loop for i=1991..2000

Esta instrucci�n (de gui�n) da acceso a un modo especial, en el cual el
programa acepta �rdenes para repetirlas o un n�mero de veces especificado,
o mientras se satisfaga una condici�n, o para valores sucesivos de la
variable �ndice i (interna). Dentro de un bucle, s�lo se pueden utilizar
seis instrucciones: genr, ols, print, smpl, store y summary (store no puede
usarse en un bucle de 'while'). Con genr y ols se pueden hacer muchas cosas.
Se sale del modo de introducci�n de �rdenes de bucle con la instrucci�n
"endloop": en este punto se ejecutar�n las �rdenes de todo el bloque. Los
bucles construidos mediante "loop" no pueden estar anidados.

La instrucci�n ols proporciona un cuadro de resultados especial, dependiendo
del tipo de bucle. Si se especifica un "n�mero_de_veces" no se muestran los
resultados de cada regresi�n particular, en su lugar se ofrece
(a) el valor medio de cada coeficiente estimado a lo largo de todas las
    iteraciones.
(b) las desviaciones t�picas de esos coeficientes estimados.
(c) el valor medio de las desviaciones t�picas estimadas de cada coeficiente.
(d) las desviaciones t�picas de las desviaciones t�picas estimadas.
Esto tiene sentido s�lo si hay alguna entrada aleatoria en cada iteraci�n.
Esta instrucci�n est� dise�ada para el an�lisis de Monte Carlo. Si se da una
condici�n "while" se muestran los resultados del modelo especificado
obtenidos en la �ltima estimaci�n del bucle: esto est� dise�ado para m�nimos
cuadrados iterativos.

La instrucci�n "print" tambi�n se comporta de modo diferente en el contexto
de un bucle con "n�mero_de_veces". Muestra la media y desviaci�n t�pica de
la variable a lo largo de las repeticiones del bucle. Est� pensado para
usarlo con variables que toman un solo valor en cada iteraci�n, por ejemplo
la scr (suma de cuadrados de los residuos) de una regresi�n. La
instrucci�n "print" funciona de la forma usual para los dem�s tipos de bucles.

La instrucci�n "store" (s�lo se puede usar una por bucle y s�lo en los bucles
de tipo "n�mero_de_veces") escribe los valores de las variables especificadas,
en cada iteraci�n del bucle, al fichero especificado. De este modo, guarda una
copia completa de las variables. Luego puede leerse mediante gretl este fichero
de datos para analizarlo.

Ejemplo de programa de bucle (Monte Carlo):

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
meantest
@Tests
Uso:            meantest x1 x2
                meantest x1 x2 -o

Calcula el estad�stico t para contrastar la hip�tesis nula de que las medias
poblacionales de las variables x1 y x2 son iguales y muestra su valor p. Sin
la opci�n -o, el estad�stico se calcula bajo el supuesto de que las varianzas
son iguales para las dos variables; con la opci�n -o se supone que las
varianzas son distintas. (La opci�n implicara diferencia s�lo si hay
diferentes n�meros de observaciones no perdidas para las dos variables)


#
mpols
@Estimation
(M�nimos cuadrados ordinarios de alta precisi�n)
Uso:            mpols vardep varindep
Ejemplos:       ols 1 0 2 4 6 7
                ols y 0 x1 x2 x3

Calcula los estimadores de m�nimos cuadrados ordinarios con 'vardep' como
variable dependiente y 'varindep' como lista de variables independientes,
utilizando para ello operaciones aritm�ticas con precisi�n m�ltiple. Las 
variables pueden especificarse por su nombre o por su n�mero; utilice el n�mero 
cero para el t�rmino constante. Esta instrucci�n s�lo est� disponible si gretl 
se configura con soporte para GMP, la biblioteca GNU de precisi�n m�ltiple.

Para estimar un ajuste polin�mico utilizando aritm�tica de precisi�n m�ltiple al 
generar las potencias necesarias de la variable independiente utilice, por 
ejemplo, la forma

mpols y 0 x ; 2 3 4

Esto hace una regresi�n de y sobre x, x al cuadrado, x al cubo y x a la cuarta 
potencia. Es decir, los n�meros a la derecha del punto y coma(que deben ser
enteros y positivos) especifican las potencias de x a utilizar. Si se especifica 
m�s de una variable independiente, la �ltima variable antes del punto y coma es 
la que ser� elevada a las potencias que se indican.

#
multiply
@Transformations
Uso:            multiply x sufijo vars
Ejemplos:       multiply invpop pc 3 4 5 6
                multiply 1000 big x1 x2 x3

Las variables de la lista "vars" (referenciadas por nombre o n�mero) se
multiplican por x, que puede ser o un valor num�rico o el nombre de una
variable previamente definida. Los productos se nombran con el prefijo
que se suministre (m�ximo tres caracteres). Si hace falta, se truncan los
nombres de las variables originales. Por ejemplo, supongamos que Vd desea
crear las versiones "per c�pita" de algunas variables y Vd tiene la variable
"pob" (poblaci�n). Las instrucciones adecuadas son

  genr invpob = 1/pob
  multiply invpob pc renta gasto

que crear�n las variables

rentapc = renta * invpob, gastopc = gasto * invpop.

#
noecho
@Programming
Uso:          noecho

Suprime el eco normal al introducir las instrucciones, cuando se ejecuta
un gui�n gretl. Ver tambi�n la orden "print" (variante de cadena literal).


#
nulldata
@Dataset
Uso:            nulldata tama�o_serie
Ejemplo:        nulldata 100

Establece un conjunto de datos vac�o, que contiene s�lo una constante, con
periodicidad 1 y el numero de observaciones especificado. Este puede usarse
al hacer simulaci�n: algunas de las �rdenes genr (p.ej. genr uniform(),
genr normal(), genr time) generar�n datos ficticios desde cero para llenar
el conjunto de datos. La orden "nulldata" tambi�n puede ser �til al utilizarla
conjuntamente con "loop".


#
ols
@Estimation
Uso:            ols vardep varindep       o      ols -o vardep varindep
Ejemplos:       ols 1 0 2 4 6 7                  ols -o 1 0 2 4 6 7
                ols y 0 x1 x2 x3                 ols -o y 0 x1 x2 x3

Calcula las estimaciones de m�nimos cuadrados ordinarios con 'vardep'
como variable dependiente y siendo 'varindep' una lista de variables 
independientes. Con la opci�n -o se mostrar� la matriz de covarianzas 
de los coeficientes de regresi�n. Las variables pueden introducirse 
mediante nombres o n�meros. Utilice el n�mero cero para incluir un 
t�rmino constante. El programa muestra tambi�n los valores p para los 
estad�sticos t (a dos colas) y F. Un valor p inferior a 0.01 indica 
significatividad al nivel del 1 por ciento y se denota mediante tres *.
Dos * indican significatividad a niveles entre el 1 y el 5 por ciento.
Tambi�n se muestran los estad�sticos de selecci�n de modelos descritos
en el libro de Ramanathan, secci�n 4.4.

 
#
omit
@Tests
Uso:            omit varlist
Ejemplos:       omit 5 7 9
		omit xx yy zz

Se ejecuta despu�s de estimar un modelo. Las variables de la lista 'varlist'
se omitir�n del modelo anterior y se estimar� el nuevo modelo. Si se omite m�s
de una variable, se muestra el estad�stico F de Wald para las variable omitidas
junto con su valor p (s�lo para el m�todo MCO). Un valor p inferior a 0.05 
implica que los coeficientes son conjuntamente significativos al nivel del 
5 por ciento.


#
omitfrom
@Tests
Uso:            omitfrom ID_de_modelo varlist
Ejemplo:        omitfrom 2 5 7 9

Funciona como la orden "omit", pero aqu� se puede especificar un modelo previo
(usando su n�mero de ID, que se muestra al principio de los resultados del 
modelo) para tomarlo como base al omitir las variables. En el ejemplo de 
arriba se omiten las variables con n�meros 5, 7 y 9 del modelo 2.


#
open
@Dataset
Uso:          open datafile

Abre un fichero de datos. Si ya hay un fichero de datos abierto, se reemplaza 
por el nuevo. El programa intentar� detectar el formato del fichero de datos
("nativo", CSV o BOX1) y lo tratar� como corresponda.


#
panel
@Dataset
Uso:            panel
                panel -s
	        panel -c

Propone que el conjunto de datos actual sea tratado como un panel (combinando
datos de secci�n cruzada y de series temporales). Sin ninguna opci�n o con 
la opci�n -s, los datos se consideran formados por series temporales apiladas
(bloques sucesivos de datos contienen series temporales para cada unidad de 
secci�n cruzada). Con la opci�n -c, los datos se leen como datos de secci�n 
cruzada apilados (bloques sucesivos contienen secciones cruzadas para cada 
periodo temporal). Ver tambi�n la instrucci�n "setobs".


#
pergm
@Statistics
Uso:            pergm varname
                pergm varname -o

Calcula y muestra (y, en modo interactivo, representa) el espectro de la 
variable especificada. Sin la opci�n -o se obtiene el periodograma muestral;
con la opci�n -o se usa una ventana de retardos de Bartlett de tama�o
2*sqrt(tama�o muestral) para estimar el espectro (ver Cap�tulo 18 del libro
"An�lisis Econom�trico" de Greene). Cuando se muestra el periodograma 
muestral, se ofrece tambi�n un contraste t de integraci�n fraccional: la
hip�tesis nula es que el orden de integraci�n de la serie es cero.


#
plot
@Graphs
Uso:            plot x1       plot x1 x2
                plot 3 7      plot -o x1 x2

Representa ( en un gr�fico de tipo texto) los valores de los datos de 
las variables especificadas, para el rango de observaci�n actualmente en 
efecto. Cada l�nea se refiere a una observaci�n y los valores se 
representan horizontalmente. Si se utiliza la opci�n -o, x1 y x2 se
representan en la misma escala, en otro caso, cada una de ellas se escala
adecuadamente. La opci�n -o s�lo deber�a usarse si las variables tienen 
aproximadamente el mismo rango de valores (p.ej. la variable dependiente 
observada y su predicci�n)


#
pooled
@Estimation
Uso:         pooled vardep varindeps

Estima un modelo mediante MCO (ver la instrucci�n "ols" para detalles sobre 
su sintaxis) y lo marca como 'modelo de panel' (o 'combinado') de manera que
la opci�n de contraste "diagn�sticos de panel" est� disponible. Para consulta 
sobre dichos diagn�sticos ver la orden "hausman".

#
print
@Printing

Imprime (=muestra en pantalla) los valores de las variables especificadas para
el rango primobs-ultobs actual, o muestra una cadena literal.

print              escribe el conjunto de datos actual en forma tabular
print 3 6          escribe los valores de las variables n�meros 3 y 6
print x y z        escribe los valores de las variables denominadas x, y y z

Si se utiliza la opci�n -o las variables se muestran en columnas, en caso
contrario, se muestran en bloques consecutivos. Si se indica la opci�n "-t", los 
datos se muestran con 10 valores significativos.

print "Cadena literal" escribe la cadena especificada.  Las comillas finales no
son necesarias, pero se necesitan las iniciales para obtener este resultado.

#
probit
@Estimation
Uso:          probit vardep varindeps

Regresi�n probit: la variable dependiente debe ser una variable binaria. Se
calculan, mediante m�nimos cuadrados iterativos (m�todo EM � de
expectativa-maximizaci�n), los estimadores de m�xima verosimilitud de los
coeficientes de 'varindeps'. Como el modelo es no lineal, las pendientes
dependen de los valores de las variables independientes: las pendientes que
se muestran se eval�an en las medias de estas variables. El estad�stico
Chi-cuadrado contrasta la hip�tesis nula de que todos los coeficientes, excepto
el t�rmino constante, son cero.


#
pvalue
@Utilities
Uso en modo interactivo:      pvalue
Uso en modo batch (o en modo consola gretl):
     Distribuci�n normal:     pvalue 1 valor_x
     Distribuci�n t:          pvalue 2 g.l. valor_x
     Chi-cuadrado:            pvalue 3 g.l. valor_x
     Distribuci�n F:          pvalue 4 gln gld valor_x
     Distribuci�n Gamma:      pvalue 5 media varianza valor_x

Calcula el �rea a la derecha del 'valor_x' en la distribuci�n especificada.
g.l. son los grados de libertad, gln son los grados de libertad del numerador,
gld son los grados de libertad del denominador.

#
quit
@Utilities

Salir de gretl. ('q' es un atajo; 'x' sale sin preguntar si se deben guardar
los resultados.)

#
reset
@Tests
Uso:          reset

Debe utilizarse inmediatamente despu�s de una instrucci�n ols. Realiza el
contraste RESET de especificaci�n de modelos (no linealidad) de Ramsey a�adiendo
a la regresi�n el cuadrado y el cubo de los valores ajustados y calculando el
estad�stico F para la hip�tesis nula de que los par�metros de los dos t�rminos
a�adidos son cero.

#
rhodiff
@Transformations
Uso:           rhodiff listarho ; listavar
Ejemplos:      rhodiff .65 ; 2 3 4
               rhodiff r1 r2 ; x1 x2 x3

Crea las transformaciones rho-diferenciadas de las variables contenidas en
'listvar'(referenciadas mediante n�mero o nombre) y las a�ade al conjunto de
datos utilizando el prefijo # para las nuevas variables. Dada la variable
v1 en listavar y r1 y r2 en listarho, se crea

v1# = v1(t) - r1*v1(t-1) - r2*v1(t-2)

Las entradas de listarho pueden darse como valores num�ricos o mediante
los nombres de variables previamente definidas.

#
rmplot
@Graphs
(Gr�fico Rango-Media)
Uso:           rmplot nombrevar

Esta instrucci�n crea un gr�fico simple para ayudar a decidir si una serie
temporal, y(t), tiene o no varianza constante. Se toma la muestra completa
t=1,...,T y se divide en peque�as submuestras de tama�o arbitrario k [gretl
elige k=sqrt(T)]. La primera submuestra se forma con y(1),...,y(k), la segunda
con y(k+1),...,y(2k), y as� sucesivamente. Para cada submuestra se calcula la
media muestral y el rango (=m�ximo-m�nimo) y se construye un gr�fico con las
medias en el eje horizontal y los rangos en el vertical. De esta forma, cada
submuestra est� representada por un punto en este plano. Si la varianza de la
serie fuera constante los rangos de las submuestras no deber�an depender de sus
medias; si se observa que los puntos se aproximan a una recta con pendiente
creciente, esto sugiere que la varianza de la serie aumenta cuando la media
aumenta; si los puntos se aproximan a una recta con pendiente decreciente, esto
sugiere que la varianza est� disminuyendo cuando la media aumenta.

Adem�s del gr�fico, gretl presenta las medias y los rangos para cada submuestra,
el coeficiente estimado para la pendiente en una regresi�n MCO de los rangos
sobre las medias y el valor p para el contraste de la hip�tesis nula de que esta
pendiente es cero. Si el coeficiente de pendiente es significativo al nivel de
significaci�n del 10 por ciento, en el gr�fico se muestra tambi�n la recta
ajustada en la regresi�n de los rangos sobre las medias.

#
run
@Programming
Uso:          run fichero

Si el fichero "fichero" contiene instrucciones gretl, esta orden (que se
invoca desde dentro de gretl) las ejecutar� de una en una. Esta es una forma
muy �til de ejecutar instrucciones 'batch' desde una sesi�n interactiva.

#
runs
@Tests
Uso:            runs nombre_var

Realiza el contraste de rachas no param�trico sobre aleatoriedad de la
variable especificada. Si Vd desea contrastar la aleatoriedad de las
desviaciones respecto a la mediana de una variable x1 que tiene una mediana
distinta de cero, lo puede hacer de la siguiente forma

genr signx1 = x1 - median(x1)
runs signx1

#
scatters
@Graphs
Uso:            scatters vary ; listavarx    o
		scatters listavary ; varx
Ejemplos:       scatters 1 ; 2 3 4 5
                scatters 1 2 3 4 5 6 ; time

Representa (mediante gnuplot) gr�ficos bivariantes (scatters) por parejas de 
variables, de vary con respecto a todas las variables de listavarx, o de 
todas las variables de listavary con respecto a varx. En el primer ejemplo de
arriba se sit�a la variable 1 en el eje y y se dibujan cuatro gr�ficos, el 
primero con la variable 2 en el eje x, el segundo con la variable 3 en el 
eje x, y as� sucesivamente. Revisar un conjunto de gr�ficos como �ste puede 
ser interesante al hacer un an�lisis de datos exploratorio. El n�mero m�ximo 
de gr�ficos es seis; cualquier variable extra en la lista ser� ignorada.


#
seed
@Programming
Uso:            seed entero

Establece la semilla para el generador de n�meros pseudo-aleatorios para
la distribuci�n uniforme y normal (ver la instrucci�n 'genr'). Por defecto
la semilla se establece cuando se inicia el programa, dependiendo de la hora
del sistema. Si se desean obtener secuencias de n�meros pseudo-aleatorios
repetibles ser� necesario establecer la semilla de forma manual.


#
setobs
@Dataset
Uso:            setobs periodicidad primobs
Ejemplos:       setobs 4 1990.1
                setobs 12 1978.03
                setobs 20 1.01

Utilice esta orden para forzar al programa a interpretar el conjunto de datos 
actual como de series temporales o de panel, cuando los datos se han le�do 
inicialmente como series simples sin fecha. La "periodicidad" debe ser un 
n�mero entero; "primobs" es una cadena que representa la fecha o 
identificaci�n de panel de la primera observaci�n. Utilice un d�gito despu�s 
del punto en "primobs" para los datos con periodicidad menor que 10, dos 
d�gitos (con un primer cero si es necesario) para periodicidad entre 10 y 99.

En caso de utilizar datos diarios se requiere una forma especial de la cadena
"primobs", concretamente la fecha ha de seguir el patr�n YY/MM/DD, por ejemplo
si los datos comienzan el 15 de febrero de 1955 ser�a "55/02/15" . Si la 
parte YY es menor que 50 se supone que el a�o pertenece al siglo XXI, en caso
contrario se supone que est� en el siglo XX. (Con datos diarios, se aceptan
las dos periodicidades, 5 y 7)


#
setmiss
@Dataset
Uso:           setmiss -1
               setmiss 100 varx

Esta orden se utiliza para hacer que el programa interprete un valor num�rico
espec�fico (el primer par�metro de la instrucci�n) como c�digo de "valor
perdido" al usar datos importados. Si este valor es el �nico par�metro, como
en el primer ejemplo de arriba, la interpretaci�n se considerar� para todas
las series del conjunto de datos. Si se usa, como segundo par�metro, una lista
de variables (por nombre o n�mero) la interpretaci�n se restringe a las
variables especificadas. As�, en el segundo ejemplo se interpreta "100" como
"valor perdido" pero s�lo para la variable denominada varx.


#
shell
@Utilities
Uso:		! [shell command]

Un "!" al comienzo de la l�nea de instrucciones gretl se interpreta como 
una salida al "shell" del usuario. As� se pueden ejecutar instrucciones
del shell desde dentro de gretl.

#
sim
@Dataset
Uso:            sim primobs ultobs y a0 a1 a2 ...

Simula valores para y para los periodos desde "primobs" hasta "ultobs". La 
variable y debe haber sido definida antes con los valores iniciales 
apropiados. "primobs" y "ultobs" deben ser consistentes con la periodicidad.
La f�rmula que se usa es:

     y(t) = a0(t) + a1(t)*y(t-1) + a2(t)*y(t-2) + ...

ai(t) pueden ser constantes o los nombres de variables definidas previamente.

Ejemplos:

  sim 1979.2 1983.1 y 0 0.9  [genera y(t) = 0.9*y(t-1)]
  sim 15 25 y 10 0.8 x       [genera y(t) = 10 + 0.8*y(t-1) + x(t)*y(t-2)]

#
smpl
@Dataset
Uso:           smpl primobs ultobs
               smpl -o var_ficticia
               smpl -o
               smpl -r <expresi�n booleana>
	       smpl full

Restablece el rango muestral. En la primera forma, "primobs" y "ultobs"
deben ser consistentes con la periodicidad de los datos. En la segunda forma,
"var_ficticia" debe ser una variable indicador con valores 0 � 1: la
muestra se restringir� a aquellas observaciones en las que el valor
indicador sea 1. La tercera forma, smpl -o, quita todas las observaciones 
para las cuales los valores de una o m�s variables est�n 'perdidos'.
La cuarta forma, usando la opci�n -r, restringe la muestra a los casos 
que satisfagan la condici�n dada. La �ltima forma, "smpl full",
recupera el rango completo de los datos.



    smpl 3 10                para datos con periodicidad 1
    smpl 1950 1990           para datos anuales con periodicidad 1
    smpl 1960.2 1982.4       para datos trimestrales
    smpl 1960.04 1985.10     para datos mensuales
    smpl 1960.2 ;            para dejar la observaci�n final sin cambiar
    smpl ; 1984.3            para dejar la observaci�n inicial sin cambiar
    smpl -o dum1             para crear una muestra basada en "dum1"
    smpl -r sqft>1400        para restringir la muestra a los casos en los que 
                             la variable sqft tenga un valor mayor que 1400


Hay que se�alar un punto especial sobre las formas "-o" y "-r" de smpl:
cualquier informaci�n "estructural" en el fichero de cabecera de datos (que
concierna a la naturaleza de series temporales o de panel de los datos) se
pierde al ejecutar esta orden. Se puede reimponer de nuevo la estructura
utilizando la instrucci�n "setobs".


#
spearman
@Statistics
Uso:            spearman x y
                spearman x y -o

Calcula el coeficiente de correlaci�n por rangos de Spearman para las dos 
variables x e y. No es necesario ordenar las variables y asignar los 
rangos manualmente; la instrucci�n ya tiene en cuenta esto. Si se 
proporciona la opci�n -o, se muestran los datos originales junto a los
ordenados.

La ordenaci�n autom�tica se hace de mayor a menor (es decir, al dato mayor se 
le asigna rango 1). Si Vd necesita invertir ese orden, puede crear una nueva 
variable cuyos valores sean los de la variable original cambiados de signo.
Por ejemplo:

  genr altx = -x
  spearman altx y

#
square
@Transformations
Uso:          square x y       o     square -o x y

Genera nuevas variables que son los cuadrados y productos cruzados de las 
variables seleccionadas (-o crea los productos cruzados). En el ejemplo de
arriba las variables creadas ser�n sq_x = x^2, sq_y = y^2 y x_y = x * y. 
Si una de las variables es una variable ficticia no se tomar� su cuadrado,
ya que se obtendr�a la misma variable.


#
store
@Dataset
Usos:            store nombre_fichero opci�n
                 store nombre_fichero opci�n lista_var
Ejemplos:        store misdatos.gdt
                 store misdatos.csv -c
                 store misdatosbin.gdt -o 2 3 4

"nombre_fichero" es el nombre del fichero en el que se guardar�n las 
variables. Si no se proporciona "lista_var" se guardar�n los valores de 
todas las variables, en caso contrario s�lo se grabar�n al fichero las
variables especificadas.

Los valores posibles de "opci�n" son:

ninguno: los datos se guardan en formato xml
  -z  : como el anterior, pero usando compresi�n de datos tipo gzip
  -o  : los datos se guardan como binarios, en doble precisi�n
  -s  : los datos se guardan como binarios, en precisi�n simple
  -c  : los datos se guardan en formato 'valores separados por comas' (CSV), 
        que pueden leerse directamente mediante cualquier programa de 
	hoja de c�lculo.
  -r  : los datos se guardan en el formato nativo de GNU R. As� se pueden
        cargar utilizando la instrucci�n de R 'source()'.
  -m  : los datos se guardan en el formato nativo de GNU Octave. La primera
  	variable citada se toma como variable dependiente y se escribe como
	vector columna; los datos restantes se escriben como una matriz,
	denominada 'X', con una variable por columna.

  -t  : los datos se guardan en el formato "tradicional" de gretl, como en el
  	programa ESL de Ramanathan, con un fichero de datos ascii y un fichero
	de "cabecera" (header .hdr).

#
summary
@Statistics
summary          muestra los estad�sticos principales para todas las variables
summary 3 7 9    muestra los estad�sticos principales para las variables 
                 n�mero 3, 7, y 9
summary x y z    muestra los estad�sticos principales para las variables 
                 x, y y z

Como resultado se ofrecen los siguientes estad�sticos: media, desviaci�n
t�pica (dt), coeficiente de variaci�n (CV= CURTOSIS), mediana, m�nimo,
m�ximo, coeficiente de asimetr�a y exceso de curtosis.


#
tabprint
@Printing
Uso:            tabprint
                tabprint -o

Debe ejecutarse despu�s de la estimaci�n de una modelo por medio de 'ols'. 
Copia el modelo estimado en forma de entorno tabular de LaTeX, a un fichero
con nombre "model_N.tex", donde N es el n�mero de modelos estimados hasta el
momento en la sesi�n actual. Esto puede incorporarse en un documento LaTeX.
Ver tambi�n la orden 'eqnprint'.

Si se proporciona la opci�n -o el fichero que se guarda es un documento 
completo LaTeX listo para ser procesado; en caso contrario debe incluirse 
dentro de un documento.


#
testuhat
@Tests
Uso:          testuhat

Debe seguir a una orden de estimaci�n de modelos. Da la distribuci�n de 
frecuencias de los residuos del modelo y un contraste Chi-cuadrado de
normalidad.

#
tsls
@Estimation
Uso:            tsls vardep listavar1 ; listavar2       [-o es opcional]
Ejemplo:        tsls y1 0 y2 y3 x1 x2 ; 0 x1 x2 x3 x4 x5 x6

Calcula las estimaciones de los par�metros de m�nimos cuadrados en dos
etapas (MC2E). "vardep" es la variable dependiente, "listavar1" es la lista de
variables independientes (incluyendo variables end�genas del lado derecho
de la ecuaci�n) en la ecuaci�n estructural para las cuales se necesitan las
estimaciones MC2E. "listavar2" es la lista combinada de variables ex�genas y
predeterminadas en todas las ecuaciones. Si "listavar2" no es al menos tan
larga como "listavar1", el modelo no est� identificado. La opci�n -o mostrar�
la matriz de covarianzas de los coeficientes. En el ejemplo de arriba, las ys
son las variables end�genas y las xs son las variables ex�genas y
predeterminadas.


#
var
@Estimation
Uso:            var orden vardep varindep
Ejemplos:       var 4 x1 const time x2 x3
                var 3 1 0 2 3 4

Organiza y estima (v�a MCO) una autorregresi�n vectorial. El primer 
argumento especifica el orden del retardo, despu�s se proporciona la 
estructura para la primera ecuaci�n, de igual forma que en la instrucci�n 
'ols'. No hay que incluir retardos entre los elementos de la lista 
"varindep" -- se a�adir�n autom�ticamente. Se ejecutar� una regresi�n 
para cada variable de la lista, excluyendo la constante, la tendencia 
temporal y las posibles variables ficticias. Los resultados de cada 
ecuaci�n incluyen los contrastes F para restricciones cero de todos
los retardos de cada variable y un contraste F para el m�ximo retardo.


#
varlist
@Dataset
Uso:          varlist

Muestra una lista de las variables definidas actualmente. "list" y "ls"
son sin�nimos.


#
vartest
@Tests
Uso:          vartest x1 x2

Calcula el estad�stico F para la hip�tesis nula de que las varianzas 
poblacionales de las variables x1 y x2 son iguales y muestra su valor p.


#
wls
@Estimation
Uso:          wls varpesos vardep varindep            [-o opcional]

Se calculan los estimadores de m�nimos cuadrados ponderados siendo
"varpesos" la variable de ponderaciones, "vardep" la variable
dependiente y "varindep" la lista de variables independientes. M�s
concretamente, se ejecuta una regresi�n MCO de varpesos*vardep
con respecto a varpesos*varindep. 

Si la variable de ponderaciones es una variable ficticia, esto es 
equivalente a eliminar todas las observaciones que tengan valor cero 
para "varpesos".

Con la opci�n -o se mostrar� la matriz de covarianzas de los 
coeficientes. Se pueden recuperar algunas variables internas 
utilizando la instrucci�n 'genr'. Para ello es necesario utilizar 
la orden 'genr' inmediatamente despu�s de esta instrucci�n. Escriba 
"help genr" para ver m�s detalles sobre esto.

  

