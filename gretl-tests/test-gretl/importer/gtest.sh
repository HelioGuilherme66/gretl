# Convert (almost) all of the xls files to gnumeric and
# diff the gnumeric import results against the xls import.
# You should inspect gtest.out after running the script.

rm -f gtest.out

for f in *.xls ; do
   base=${f%.xls}
   echo "$base"
   if [ "$base" = "pooled_data" ] || [ "$base" = "franck" ] ; then
      echo " skipping (contains formulae)"
   else 
     ssconvert $f ${base}.gnumeric
     cat ${base}.inp | sed s/xls/gnumeric/ > gtmp.inp
     gretlcli -b gtmp.inp > gtmp.out
     diff -bB newout/${base}.out gtmp.out >> gtest.out
     rm ${base}.gnumeric gtmp.inp gtmp.out
   fi
done
   
