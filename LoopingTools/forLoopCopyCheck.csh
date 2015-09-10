#!/bin/tcsh

#set myDirs = (pp2015B-PCLupdate pp2015B-PCLupdateAndTemplates pp2015B-Prompt pp2015B-PromptAndTemplates pp2015B-Run1Ali pp2015B-Run2Ali0TCollisions pp2015B-Run2AliCRUZET pp2015B-Run2AliCosmics pp2015B-hp1368 pp2015B-hp1370 pp2015B-hp1375 pp2015B-hp1387 pp2015B-hp1388 pp2015B-mp1799 pp2015B-hp1389 pp2015B-mp1807)

#set myDirs=(pp2015B-PCLAlignment pp2015B-hp1394 pp2015B-mp1819 pp2015B-mp1820  pp2015B-hp1388)
#set myDirs=(pp2015B-hp1398 pp2015B-mp1826 pp2015B-mp1826NoZ)
#set myDirs=(pp2015C-PCLAlignment)

set myDirs=(pp2015C-PCLLike255019 pp2015C-PCLAlignment)

touch report.txt
set COUNT=0
set COUNTREMOVED=0
set COUNTbyFolder=0

foreach i (`seq $#myDirs`)
    echo "=================================================" 
    echo "Reading from $myDirs[$i]"
    set COUNTbyFolder=0 
    mkdir ./$myDirs[$i]
    foreach inputfile (`eos ls /store/caf/user/musich/Alignment/PVValidation/$myDirs[$i]`)
	echo "Analyzing: /store/caf/user/musich/Alignment/PVValidation/$myDirs[$i]/$inputfile" 	
	if ($inputfile =~ *"root"*) then
	    @ COUNT+=1
	    @ COUNTbyFolder+=1
	    set myLFN = `cmsPfn /store/caf/user/musich/Alignment/PVValidation/$myDirs[$i]/$inputfile`
	    root -q -b ./check.C\(\"$myLFN\"\)  > numevents.out
	    set numevents=`tail -5 numevents.out | grep events | awk '{split($0,a," "); print a[4]}'`
	    set numTracks=`tail -5 numevents.out | grep tracks | awk '{split($0,b," "); print b[4]}'`
	    set trkeffic=`tail -5 numevents.out | grep trks    | awk '{split($0,c,":"); print c[2]}'`
	    #alias MATH 'set \!:1 = `echo "scale=5;\!:5-$" | bc -l`' 
	    echo "ALCARECO events: $numevents"
	    echo "ALCARECO tracks: $numTracks"
	    if($numevents>0) then 
		echo "trk eff (ALCARECO): $trkeffic"
		cmsStage /store/caf/user/musich/Alignment/PVValidation/$myDirs[$i]/$inputfile ./$myDirs[$i]
	    else 
		echo "0 events!!!! => going to remove the file:  /store/caf/user/musich/Alignment/PVValidation/$myDirs[$i]/$inputfile "
		@ COUNTREMOVED+=1
		cmsRm /store/caf/user/musich/Alignment/PVValidation/$myDirs[$i]/$inputfile
	    endif
	    echo "-------------------------------------------------" 
	    #echo $inputfile $numevents $numTracks $trkeffic >> report.txt
	    printf '%-5s %-5s %-20s %-20s %-20s\n' $COUNTbyFolder $inputfile $numevents $numTracks >> report_$myDirs[$i].txt
	endif
    end
end

echo "Removed $COUNTREMOVED files" 
echo "Total checked files: $COUNT"
