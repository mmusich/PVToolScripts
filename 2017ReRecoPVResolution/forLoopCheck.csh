#!/bin/tcsh
set myDirs=(Prompt Start ReReco)

touch report.txt
set COUNT=0
set COUNTREMOVED=0
set COUNTbyFolder=0

foreach i (`seq $#myDirs`)
    echo "=================================================" 
    echo "Reading from $myDirs[$i]"
    set COUNTbyFolder=0 
    foreach inputfile (`ls $PWD/$myDirs[$i]`)
	echo $inputfile
	if ($inputfile =~ *"root"*) then
	    @ COUNT+=1
	    @ COUNTbyFolder+=1
	    set myLFN = $PWD/$myDirs[$i]/$inputfile
	    echo "Analyzing: $myLFN" 	
	    root -q -b ./checkFile.C\(\"$myLFN\"\)  > numevents.out
	    set numevents=`tail -5 numevents.out | grep events | awk '{split($0,a," "); print a[4]}'`
	    set numVertices=`tail -5 numevents.out | grep vertices | awk '{split($0,b," "); print b[4]}'`
	    set vtxeffic=`tail -5 numevents.out | grep vtx | awk '{split($0,c,":"); print c[2]}'`
	    echo "ALCARECO events: $numevents"
	    echo "ALCARECO vertices: $numVertices"
	    if($numevents>0) then 
		echo "vtx eff (ALCARECO): $vtxeffic"
	    else 
		echo "0 events => going to remove the file: /store/caf/user/musich/Alignment/PVValidation/$myDirs[$i]/$inputfile "
		@ COUNTREMOVED+=1
		rm -fr $myDirs[$i]/$inputfile
	    endif
	    echo "-------------------------------------------------" 
	    #echo $inputfile $numevents $numTracks $trkeffic >> report.txt
	    printf '%-5s %-5s %-20s %-20s %-20s\n' $COUNTbyFolder $inputfile $numevents $numVertices >> report_$myDirs[$i].txt
	endif
    end
end

echo "Removed $COUNTREMOVED files" 
echo "Total checked files: $COUNT"
