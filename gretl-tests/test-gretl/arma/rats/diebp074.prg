*
* Examples of trends
* Pages 74, 77, 79
*
* These use some rather advanced graphics techniques to produce
* the types of effects shown in the text book. If you eliminate
* the spgraph and grtext instructions, you'll just get the
* basic graphs.
*
all 100
set mtrend = 10-.25*t
set ptrend   = -50+.8*t
*
* SPGRAPH (Special Purpose GRAPH) is used for a variety of specialty
* graphs. Here, it's used to allow GRTEXT instructions to annotate
* the graph. The graph isn't done until the SPGRAPH(DONE) instruction.
*
spgraph
graph(nodates) 2
# mtrend
# ptrend
*
* GRTEXT puts text on the previous graph. The Y option indicates the y
* value (on the scale of the graph) and ENTRY the time period at which
* the text is to be placed. If you're annotating a SCATTER graph, you
* use Y and X to set the location.
*
grtext(y=15,entry=10,align=left) 'Trend=10-.25*Time'
grtext(y=-40,entry=20,align=left) 'Trend=-50+.8*Time'
spgraph(done)
*
set pptrend = 10+.3*t+.3*t**2
set pmtrend = 10+30*t-.3*t**2
set mmtrend = 10-.4*t-.4*t**2
set mptrend = 10-25*t+.3*t**2
*
* And here SPGRAPH is being used to create a 2x2 matrix of graphs. This
* is explained in more detail on diebp050.prg.
*
spgraph(hfields=2,vfields=2)
graph(nodates,header='Trend=10+.3*Time+.3*Time**2')
# pptrend
graph(nodates,header='Trend=10-.4*Time-.4*Time**2')
# mmtrend
graph(nodates,header='Trend=10+30*Time-.3*Time**2')
# pmtrend
graph(nodates,header='Trend=10-25*Time+.3*Time**2')
# mptrend
spgraph(done)
*
set epptrend = 5*exp(.02*t)
set epmtrend = 5*exp(-.02*t)
set emmtrend = -5*exp(-.02*t)
set emptrend = -5*exp(.02*t)
*
spgraph(hfields=2,vfields=2)
graph(nodates,header='Trend=5 exp(.02*Time)')
# epptrend
graph(nodates,header='Trend=-5 exp(.02*Time)')
# emptrend
graph(nodates,header='Trend=5 exp(-.02*Time)')
# epmtrend
graph(nodates,header='Trend=-5 exp(-.02*Time)')
# emmtrend
spgraph(done)


