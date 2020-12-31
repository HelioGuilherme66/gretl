set matsize 200
insheet using "unbal3.csv"
tsset id year
xtabond y
xtabond y , noconstant
