#!/bin/tcsh

set myDirs=(2017_dataReRecoTest2017ReReco)

touch report.txt
set COUNT=0
set COUNTREMOVED=0
set COUNTbyFolder=0

foreach i (`seq $#myDirs`)
    echo "=================================================" 
    echo "Reading from $myDirs[$i]"
    set COUNTbyFolder=0 
    foreach inputfile (`ls $PWD/$myDirs[$i]`)
	echo "Analyzing: $myDirs[$i]/$inputfile" 	
	if ($inputfile =~ *"root"*) then
	    @ COUNT+=1
	    @ COUNTbyFolder+=1
	    set myLFN = "$PWD/$myDirs[$i]/$inputfile"
	    echo $myLFN
	    root -q -b check.C\(\"$myLFN\"\)  > numevents.out
	    #root -q -b exec.C\(\"$myLFN\"\)  > numevents.out
	    set numevents=`tail -5 numevents.out | grep events | awk '{split($0,a," "); print a[4]}'`
	    set numTracks=`tail -5 numevents.out | grep tracks | awk '{split($0,b," "); print b[4]}'`
	    set trkeffic=`tail -5 numevents.out | grep trks    | awk '{split($0,c,":"); print c[2]}'`
	    #alias MATH 'set \!:1 = `echo "scale=5;\!:5-$" | bc -l`' 
	    echo "ALCARECO events: $numevents"
	    echo "ALCARECO tracks: $numTracks"
	    if($numevents>0) then 
		echo "trk eff (ALCARECO): $trkeffic"
	    else 
		echo "0 events!!!! => going to remove the file: /store/group/alca_trackeralign/musich/test_out/$myDirs[$i]/$inputfile "
		@ COUNTREMOVED+=1
		rm -fr $PWD/$myDirs[$i]/$inputfile
	    endif
	    echo "-------------------------------------------------" 
	    #echo $inputfile $numevents $numTracks $trkeffic >> report.txt
	    printf '%-5s %-5s %-20s %-20s %-20s\n' $COUNTbyFolder $inputfile $numevents $numTracks >> report_$myDirs[$i].txt
	endif
    end
end

echo "Removed $COUNTREMOVED files" 
echo "Total checked files: $COUNT"
