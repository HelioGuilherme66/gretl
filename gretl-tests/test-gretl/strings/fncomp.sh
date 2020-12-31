for f in backtick \
  grabroots \
  mprintf \
  normvars \
  panelstr \
  readfile \
  readfile2 \
  strings2 \
  strings \
  quotes ; do
  script="${f}_fn.inp"
  outfile="${f}_fn.out"
  orig="./newout/${f}.out"
  gretlcli -b $script > $outfile && \
  diff $orig $outfile
done
