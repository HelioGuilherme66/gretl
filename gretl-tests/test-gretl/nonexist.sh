#/bin/sh

if ! gretlcli -b nonexistent.inp ; then
  echo "Failed"
else
  echo "OK"
fi
