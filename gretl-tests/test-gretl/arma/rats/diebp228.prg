source(noecho) regrecur.src
*
* Constructed examples of recursive analysis, pp 228-232
* Note that these all depend upon random numbers, and so will not match (exactly)
* the results in the textbook, and won't match even from one run to the next. If
* you want to be able to reproduce exactly a set of results, add the instruction
*    SEEE  big integer   (like SEED 534653) at the top of the program.
*
all 200
set x = (t)
set y1 = .5*x+%ran(10.0)
spgraph(hfields=2,vfields=2,header='Figure 9.15 Recursive Analysis: Constant Parameter Model')
scatter
# x y1
linreg y1
# x
@regrecursive(cusum,cohist=cohist,sehist=sehist,sighist=sighist) resids
set upperco = cohist(1)+sehist(1)*2.0
set lowerco = cohist(1)-sehist(1)*2.0
graph(header='Recursive Estimates') 3
# cohist(1)
# upperco 10 *
# lowerco 10 *
set upperres =  2.0*sighist
set lowerres = -2.0*sighist
graph(header='Recursive Residuals') 3
# resids
# upperres
# lowerres
spgraph(done)
*
set y2 = (.02*t)*x+%ran(10.0)
spgraph(hfields=2,vfields=2,header='Figure 9.16 Recursive Analysis: Trending Parameter Model')
scatter
# x y2
linreg y2
# x
@regrecursive(cusum,cohist=cohist,sehist=sehist,sighist=sighist) resids
set upperco = cohist(1)+sehist(1)*2.0
set lowerco = cohist(1)-sehist(1)*2.0
graph(header='Recursive Estimates') 3
# cohist(1)
# upperco 10 *
# lowerco 10 *
set upperres =  2.0*sighist
set lowerres = -2.0*sighist
graph(header='Recursive Residuals') 3
# resids
# upperres
# lowerres
spgraph(done)
*
set y3 = %if(t<=100,x,3*x)+%ran(10.0)
spgraph(hfields=2,vfields=2,header='Figure 9.17 Recursive Analysis: Breaking Parameter Model')
scatter
# x y3
linreg y3
# x
@regrecursive(cusum,cohist=cohist,sehist=sehist,sighist=sighist) resids
set upperco = cohist(1)+sehist(1)*2.0
set lowerco = cohist(1)-sehist(1)*2.0
graph(header='Recursive Estimates') 3
# cohist(1)
# upperco 10 *
# lowerco 10 *
set upperres = 2.0*sighist
set lowerres = -2.0*sighist
graph(header='Recursive Residuals') 3
# resids
# upperres
# lowerres
spgraph(done)

