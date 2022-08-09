###########################################################
##
## copy the condor decription file and submit condor jobs 
##
###########################################################

#!/bin/bash

runname=$1
runpath=/a/data/bleeker/sbnd/runs/$runname

mkdir -p $runpath/condor     ##create multiple directories
j=0
for file in $runpath/*generator_NU_condor*.root; do
	printf -v i "%02d" $j   ## pad variable j with 0s and save the output to variable i.
	
	cp analyzer_condor.sub ./submit/condor_$i.sub
	filename=$(basename "$file" .root)   ##get only file name without extension.
	
	logname=$runpath/condor/$filename"_analyzer.log"
	outname=$runpath/condor/$filename"_analyzer.out"
	errname=$runpath/condor/$filename"_analyzer.err"
	sed -i "s|inputfile|$file|g" ./submit/condor_$i.sub   ## sed -i option edits file in place.
	sed -i "s|logname|$logname|g" ./submit/condor_$i.sub   ## use double quote in sed to use variables; use "|" as separator when you have "\" in your variable
	sed -i "s|outname|$outname|g" ./submit/condor_$i.sub
	sed -i "s|errname|$errname|g" ./submit/condor_$i.sub

	condor_submit ./submit/condor_$i.sub
	
	((j=j+1))   
done
