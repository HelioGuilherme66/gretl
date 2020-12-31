#!/bin/sh

wget -r -l 1 -nd -np -nc -A *.gfn,*.zip -X invalid \
  http://ricardo.ecn.wfu.edu/gretl/cgi-bin/current_fnfiles
