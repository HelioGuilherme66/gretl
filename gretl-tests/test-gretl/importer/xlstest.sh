#!/bin/sh

rm -f xls.errout

for f in `cat xls.list` ; do
  gretlcli -b $f 2>> xls.errout
done
