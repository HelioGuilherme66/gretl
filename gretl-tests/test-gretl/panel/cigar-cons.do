insheet using cigar.csv
xtset state year
quietly tabulate year, generate(dyear)
xtabond lnc lnp lnpn lny dyear4-dyear30, lag(1)
xtabond lnc lnp lnpn lny dyear4-dyear30, lag(1) twostep rob
* xtabond2 lnc L.lnc lnp lnpn lny dyear3-dyear30, gmm(L.lnc) iv(lnp lnpn lny dyear3-dyear30) nocons noleveleq twostep

