#!/bin/bash

if [ "x$1" != "x" ] ; then
   # executing in build tree
   INPPATH=$1
   GRETLCLI="$(pwd)/../cli/gretlcli"
else
   # executing "here" using installed gretlcli
   INPPATH="."
   GRETLCLI="gretlcli"
fi

# Remove the 'fails' file if it exists
rm -f fails

# Store the current directory path
HERE=`pwd`

# Initialize exit status variable
my_status=0

# Test scripts in the 'practice_scripts' directory
echo "*** practice scripts ***"
for f in `find $INPPATH/practice_scripts -name "*.inp"` ; do
   # Print the name of the script being tested
   bname=`basename $f`
   echo -n $bname

   # Run the script with gretlcli in batch and quiet mode
   ${GRETLCLI} -b -q -e $f > /dev/null 2>&1

   # Check if the script failed
   if [ $? != 0 ] ; then
      # Print 'Failed', update status variable, and log the failed script
      echo -e " [\e[0;31mFailed\e[0m]"
      my_status=1
      echo $INPPATH/practice_scripts/$bname >> $HERE/fails
   else
      # Print 'OK' if the script succeeded
      echo -e " [\e[0;32mOK\e[0m]"
   fi
done

# Test scripts in the 'commands', 'functions', and 'fundamentals' directories
for d in commands functions fundamentals ; do
   echo "*** $d ***"
   cd $INPPATH/test_scripts/$d
   for f in `find . -name "*.inp"` ; do
      bname=`basename $f`
      echo -n $bname
      ${GRETLCLI} -b -q -e $f > /dev/null 2>&1
      if [ $? != 0 ] ; then
           echo -e " [\e[0;31mFailed\e[0m]"
           my_status=1
           echo $INPPATH/test_scripts/$d/$bname >> $HERE/fails
      else
           echo -e " [\e[0;32mOK\e[0m]"
      fi
   done
   # Return to the original directory
   cd $HERE
done

# If there were any failures, print the names of the failed scripts
if test -f fails ; then
   echo -e "\e[0;31mFailed script(s):\e[0m"
   cat fails
else
   echo "No errors were found"
fi

# Exit with status code 0 if all scripts passed, 1 if any script failed.
exit $my_status
