#/bin/bash
#
# Checks all apps in aomp/openmpapps directory. Programs return 0 for success or a number > 0 for failure.
#
realpath=`realpath $0`
thisdir=`dirname $realpath`

rm -f check-openmpapps.txt times-openmpapps.txt passing-tests.txt failing-tests.txt make-fail.txt
echo ""
echo ""
let TotFails=0

cd $thisdir

echo ""
echo "RUNNING ALL TESTS IN: $path "
echo ""

tpath=$thisdir
echo "" > check-openmpapps.txt
echo "" >> check-openmpapps.txt
echo "" > $tpath/times-openmpapps.txt
echo "" >> $tpath/times-openmpapps.txt

echo "*************************************************************************************************" >> check-openmpapps.txt
echo "*******A non-zero return code means the app failed to compile or there was a runtime issue*******" >> check-openmpapps.txt
echo "*************************************************************************************************" >> check-openmpapps.txt
execute_makefile(){
  make clean
  # Get time in seconds since epoch
  TT0=`date '+%s'`
  make
  ret=$?
  if [ $ret != 0 ]; then
    echo " Return Code for $base: $ret"  >> $1
    (( TotFails++ ))
    echo $base >> "$thisdir"/make-fail.txt
    return
  fi
  TT1=`date '+%s'`
  make run
  ret=$?
  echo " Return Code for $base: $ret"  >> $1
  if [ $ret != 0 ]; then
    (( TotFails++ ))
    echo $base >> "$thisdir"/failing-tests.txt
  else
    echo $base >> "$thisdir"/passing-tests.txt
  fi

  TT2=`date '+%s'`

  # App Build Time
  BUILDMIN=$(( ($TT1 - $TT0)/60 ))
  BUILDSEC=$(( ($TT1 - $TT0)%60 ))

  # App Run Time
  RUNMIN=$(( ($TT2 - $TT1)/60 ))
  RUNSEC=$(( ($TT2 - $TT1)%60 ))

  # Accumalate total times for all apps
  (( TOTALBUILDSEC += (TT1 - TT0) ))
  (( TOTALRUNSEC += (TT2 - TT1) ))

  echo $base >> $tpath/times-openmpapps.txt
  echo Test start `date --date=@$TT0` >> $tpath/times-openmpapps.txt
  echo Compile end `date --date=@$TT1` >> $tpath/times-openmpapps.txt
  echo Executable end `date --date=@$TT2` >> $tpath/times-openmpapps.txt
  echo App $base Build Time: $BUILDMIN min  $BUILDSEC sec >> $tpath/times-openmpapps.txt
  echo App $base Run Time: $RUNMIN min  $RUNSEC sec >> $tpath/times-openmpapps.txt
  echo >> $tpath/times-openmpapps.txt
}

# Loop over all directories and execute the Makefile, directory levels where source code is found differs. Hence, why the conditional statements are used to determine where the Makefile is found.
for directory in ./*/; do 
  pushd "$directory" && path=$(pwd) && base=$(basename $path)

  # Apps that have Makefile on first level
  if [ $base == 'hpgmg-mp4' ] || [ $base == 'lulesh-mp4' ] ; then
    execute_makefile "../check-openmpapps.txt"

    # COMD has a Makefile in a src folder, which is named src-omp, on second level
    elif [ $base == 'comd-mp4' ] ; then
      src_dir='src-omp'
      cd $src_dir
      execute_makefile "../../check-openmpapps.txt"

    # Tests that have Makefile in a folder named src on second level
    else
      src_dir='src'
      cd $src_dir
      execute_makefile "../../check-openmpapps.txt"
    fi
  popd
done

# Convert total build/run times into minute/sec format
TOTALBUILDMIN=$(( $TOTALBUILDSEC/60 ))
TOTALBUILDSEC=$(( $TOTALBUILDSEC %60 ))
TOTALRUNMIN=$(( $TOTALRUNSEC/60 ))
TOTALRUNSEC=$(( $TOTALRUNSEC%60 ))

echo "Total Build (all apps): $TOTALBUILDMIN min $TOTALBUILDSEC sec" >> $tpath/times-openmpapps.txt
echo "Total Run Time (all apps): $TOTALRUNMIN min $TOTALRUNSEC sec" >> $tpath/times-openmpapps.txt
cat times-openmpapps.txt
echo ""
echo "--------Timings--------"
appregex="App (.*)"
newlineregex="App (.*)Run Time"
file="times-openmpapps.txt"
while read -r line; do
  if [[ "$line" =~ $appregex ]]; then
    appname="${BASH_REMATCH[1]}"
    echo $appname
  fi
  if [[ "$line" =~ $newlineregex ]]; then
    echo ""
  fi
done < "$file"
echo "Total Build (all apps): $TOTALBUILDMIN min $TOTALBUILDSEC sec"
echo "Total Run Time (all apps): $TOTALRUNMIN min $TOTALRUNSEC sec"

echo "-----------------------"
cat check-openmpapps.txt
exit $TotFails
