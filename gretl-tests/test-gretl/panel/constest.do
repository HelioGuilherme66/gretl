insheet using cigar.csv
xtset state year
xtabond lnc lnp lnpn lny, lag(1)
xtabond2 lnc L.lnc lnp lnpn lny, gmm(L.lnc) iv(lnp lnpn lny) noleveleq

