#!/bin/tcsh

set coordinates = ("TPB1_xCoord" "TPB1_yCoord" "TPB1_zCoord" "TPB1_phiX" "TPB1_phiY" "TPB1_phiZ" "TPB2_xCoord" "TPB2_yCoord" "TPB2_zCoord" "TPB2_phiX" "TPB2_phiY" "TPB2_phiZ" ) 

foreach coord (`seq $#coordinates`)
    echo "coordinate:  ===>"$coordinates[$coord]
    set listToAdd = ""
    set listToAddMinus = ""
    foreach file (`ls .`)
	if($file =~ *$coordinates[$coord]* ) then
	    #echo $file "to be added"
	    set scenario = `echo $file | awk '{split($0,a,"_"); print a[7]}'`
	    set scenario = "$coordinates[$coord]$scenario"
	    if($file =~ *"-"* ) then
		set listToAddMinus = "$listToAddMinus,$file=$scenario"
	    else 
		set listToAdd = "$listToAdd,$file=$scenario"
	    endif
	endif
    end	
    set substring=`echo $listToAdd | sed 's/^.//'`
    set substringMinus=`echo $listToAddMinus | sed 's/^.//'`
   
    #echo $substring
    root -b -q $PWD/FitPVResiduals.C++\(\"$substring\"\,\"\"\)
    root -b -q $PWD/FitPVResiduals.C++\(\"$substringMinus\"\,\"\"\)  
    #root -b -q FitPVResiduals.C(\"$substring\",true,true,\"\")
end	
