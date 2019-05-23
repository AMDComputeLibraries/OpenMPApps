#/bin/bash
#
# Checks all apps in aomp/openmpapps directory. Programs return 0 for success or a number > 0 for failure.
#
rm -f check-openmpapps.txt times-openmpapps.txt
echo ""
echo ""
path=$(pwd)
base=$(basename $path)
echo ""
echo "RUNNING ALL TESTS IN: $path "
echo ""

tpath=$path
echo "" > check-openmpapps.txt
echo "" >> check-openmpapps.txt
echo "" > $tpath/times-openmpapps.txt
echo "" >> $tpath/times-openmpapps.txt

echo "*************************************************************************************************" >> check-openmpapps.txt
echo "*******A non-zero return code means the app failed to compile or there was a runtime issue*******" >> check-openmpapps.txt
echo "*************************************************************************************************" >> check-openmpapps.txt

execute_makefile(){
	make clean
	TT0=`date`
	make
	TT1=`date`
	make run
	TT2=`date`
	echo $base >> $tpath/times-openmpapps.txt
	echo Test start $TT0 >> $tpath/times-openmpapps.txt
	echo Compile end $TT1 >> $tpath/times-openmpapps.txt
	echo executable end $TT2 >> $tpath/times-openmpapps.txt
        echo >> $tpath/times-openmpapps.txt
}

#Loop over all directories and execute the Makefile, directory levels where source code is found differs. Hence, why the conditional statements are used to determine where the Makefile is found.
for directory in ./*/; do 
	(cd "$directory" && path=$(pwd) && base=$(basename $path)
		
		#Apps that have Makefile on first level 
		if [ $base == 'hpgmg-mp4' ] || [ $base == 'lulesh-mp4' ] ; then
			execute_makefile
			echo " Return Code for $base: $?"  >> ../check-openmpapps.txt
			make clean
		
		#COMD has a Makefile in a src folder, which is named src-omp, on second level
		elif [ $base == 'comd-mp4' ] ; then
			src_dir='src-omp'
			cd $src_dir 
			execute_makefile
			echo " Return Code for $base: $?"  >> ../../check-openmpapps.txt
			make clean
		
		#Tests that have Makefile in a folder named src on second level	
		else
			src_dir='src'
			cd $src_dir
			execute_makefile
			echo " Return Code for $base: $?"  >> ../../check-openmpapps.txt
			make clean
		fi

	)
	
done

cat times-openmpapps.txt
cat check-openmpapps.txt
