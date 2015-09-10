#!/bin/tcsh

#set myDirs =  (pp2015B-PCLupdate pp2015B-PCLupdateAndTemplates pp2015B-Prompt pp2015B-PromptAndTemplates pp2015B-Run1Ali pp2015B-Run2Ali0TCollisions pp2015B-Run2AliCRUZET pp2015B-Run2AliCosmics pp2015B-hp1368 pp2015B-hp1370 pp2015B-hp1375 pp2015B-hp1387 pp2015B-hp1388 pp2015B-mp1799 pp2015B-hp1389 pp2015B-mp1807)

#set myDirs = (pp2015B-PCLAlignment pp2015B-hp1394 pp2015B-mp1819 pp2015B-mp1820 pp2015B-hp1388)

set myDirs = (pp2015B-hp1398 pp2015B-mp1826 pp2015B-mp1826NoZ)

set toExclude = ("251493" "251496" "251497" "251498" "251499" "251500" "252126")

foreach i (`seq $#myDirs`)
    set listToAdd = ""
    foreach file (`ls $myDirs[$i]`)
	set runbase = `echo $file | awk '{split($0,a,"_"); print a[3]}' | sed "s/.root//ig"`
	echo $runbase
	#if($file =~ {*251493*,*251496*,*251497*,*251498*,*251499*,*251500*,*252126*}) then
	if($runbase>251509) then
	    echo $file "to be added"
	    set listToAdd = "$listToAdd $myDirs[$i]/$file"
	    endif
	end	
    echo $listToAdd
    #hadd $myDirs[$i]/PVValidation_$myDirs[$i].root $myDirs[$i]/*.root
    #mv $myDirs[$i]/PVValidation_$myDirs[$i].root .                    
    hadd $myDirs[$i]/PVValidation_$myDirs[$i]_postBtrip.root $listToAdd 
    mv $myDirs[$i]/PVValidation_$myDirs[$i]_postBtrip.root .                    
end

