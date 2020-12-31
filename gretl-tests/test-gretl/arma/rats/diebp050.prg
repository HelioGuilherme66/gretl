*
* Analysis and graphics of Anscombe's data, pp 50-52
*
all 11
open data anscombe.dat
data(format=prn,org=columns)
*
linreg yone
# constant xone
linreg ytwo
# constant xtwo
linreg ythree
# constant xthree
linreg yfour
# constant xfour
*
* SPGRAPH (Special Purpose GRAPH) is used for a variety of specialty
* graphs. Here, it is used to make a 2x2 array of graphs. HFIELDS and
* VFIELDS specify the number of horizontal and vertical fields on the
* page. The SCATTER instructions then fill in the panes (by default,
* working down by columns). The SPGRAPH(DONE) instruction at the end
* means that the full page of graphs is finished and can be displayed.
*
spgraph(hfields=2,vfields=2)
scatter(hlabel='X1',vlabel='Y1')
# xone yone
scatter(hlabel='X2',vlabel='Y2')
# xtwo ytwo
scatter(hlabel='X3',vlabel='Y3')
# xthree ythree
scatter(hlabel='X4',vlabel='Y4')
# xfour yfour
spgraph(done)

