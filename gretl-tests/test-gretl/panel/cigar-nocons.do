insheet using cigar.csv
xtset state year
quietly tabulate year, generate(dyear)
xtabond lnc lnp lnpn lny dyear3-dyear30, lag(1) nocons
xtabond lnc lnp lnpn lny dyear3-dyear30, lag(1) nocons twostep rob
xtabond2 lnc L.lnc lnp lnpn lny dyear3-dyear30, gmm(L.lnc) iv(lnp lnpn lny dyear3-dyear30) rob nocons noleveleq twostep

