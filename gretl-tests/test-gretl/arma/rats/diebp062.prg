*
* GNP components graphs from ppage 62-65
*
cal 1960
all 1989:1
open data gnp.dat
data(format=prn,org=columns)
*
set ga8gff = agriculture
set ga8gr  = retail
set ga8gs  = services
set ga8gm  = manufacturing
*
* With RATS, you need to use scale=none to suppress the vertical axis value
* scale. If you want dotted scale lines running horizontally across the graph
* you can use the option "extend" but that won't work if there is no scale.
*
graph(style=bargraph,key=below,scale=none,header='Figure 3.11') 4
# ga8gff
# ga8gr
# ga8gs
# ga8gm
graph(style=stacked,key=below,scale=none,header='Figure 3.12') 4
# ga8gff
# ga8gr
# ga8gs
# ga8gm
graph(style=symbols,key=below,scale=none,header='Figure 3.13') 4
# ga8gff
# ga8gr
# ga8gs
# ga8gm
graph(style=line,key=below,header='Figure 3.14') 4
# ga8gff
# ga8gr
# ga8gs
# ga8gm
*
* This works only with RATS 5.04 or later. With earlier versions, take the
* key=attached out.
*
graph(style=line,key=attached,header='Figure 3.15',$
    klabels=||'AGRICULTURE','RETAIL','SERVICES','MANUFACTURING'||) 4
# ga8gff
# ga8gr
# ga8gs
# ga8gm
*
* This creates a dummy variable for the periods of recession. The %year function
* returns the year number of an entry, so %year(t)==1960 is 1 if t is an
* entry in 1960 and 0 if it isn't.
*
set recess = %year(t)==1960.or.%year(t)==1970.or.%year(t)==1974.or.%year(t)==1981.or.%year(t)==1982
graph(style=line,shade=recess,vlabel='GDP Component',key=attached,$
 klabels=||'Agriculture','Retail','Services','Manufacturing'||,$
 header='Figure 3.16 Components of Real GDP (Millions of Dollars, Annual)') 4
# ga8gff
# ga8gr
# ga8gs
# ga8gm



