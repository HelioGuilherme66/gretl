#!/bin/sh

rm -rf newout
mkdir newout

for f in *.dat ; do
  INP=${f%.dat}.inp
  OUT=${f%.dat}.out
  echo "open $f" > $INP
  cat uni.cmds >> $INP
  gretlcli -b $INP > newout/$OUT
done
