#!/bin/tcsh

set Objname=$1 
set theLabel=$2

foreach filetostudy (`ls ${PWD}/${Objname} | grep .root`)

     set namebase=`echo $filetostudy |awk '{split($0,a,"_"); print a[2]}'`
     if ($namebase == $Objname) then
       set datebase=`echo $filetostudy |awk '{split($0,b,"_"); print b[3]}'` 
       set theDate=`echo $datebase |awk '{split($0,c,"."); print c[1]}'`  

       echo "Studying $namebase $theDate file"

       rm -f FittedDeltaZ.tx
       root -b -q  $PWD/FitPVResiduals.C++\(\"${PWD}/${Objname}/$filetostudy=${theLabel}\"\,\"${theDate}\"\)

       rm -f numevents.out 
       root -b -q  $PWD/check3.C\(\"${PWD}/${Objname}/$filetostudy\"\) > numevents.out
       set numevents=`tail -5 numevents.out | grep events | awk '{split($0,a," "); print a[4]}'`
       set numTracks=`tail -5 numevents.out | grep tracks | awk '{split($0,b," "); print b[4]}'`
    
       set rawfits=`tail -1 FittedDeltaZ.txt`
       set deltaz=`echo $rawfits | awk '{split($0,a,"|"); print a[1]}'`
       set sigmadeltaz=`echo $rawfits | awk '{split($0,a,"|"); print a[2]}'` 
 
       setenv flag "BAD"

       if (! -d summary.txt ) then
	 touch summary.txt 
       endif
    
       if (${numevents} > 1000) then
	setenv flag "GOOD"
       endif

       echo "- File $namebase had ${numevents} events. Fit separation: $deltaz \pm $sigmadeltaz"

       echo $theDate $deltaz $sigmadeltaz ${numevents} $flag >> summary.txt
    endif
end

mv summary.txt summary_${Objname}.txt
mkdir ./Plots_${Objname}
mv *.png ./Plots_${Objname}
mv *.pdf ./Plots_${Objname}



