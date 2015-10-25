#!/bin/tcsh

#set myDirs =  (pp2015B-PCLupdate pp2015B-PCLupdateAndTemplates pp2015B-Prompt pp2015B-PromptAndTemplates pp2015B-Run1Ali pp2015B-Run2Ali0TCollisions pp2015B-Run2AliCRUZET pp2015B-Run2AliCosmics pp2015B-hp1368 pp2015B-hp1370 pp2015B-hp1375 pp2015B-hp1387 pp2015B-hp1388 pp2015B-mp1799)

set myDirs = (pp2015B-hp1389)

set toExclude = ("251493" "251496" "251497" "251498" "251499" "251500" "252126")

foreach i (`seq $#myDirs`)
    set listToAdd = ""
    foreach file (`ls $myDirs[$i]`)
	if($file =~ {*251493*,*251496*,*251497*,*251498*,*251499*,*251500*,*252126*}) then
	    echo $file "to be excluded"
		#listToAdd+=" $file"
	    endif
	end	
    #echo $listToAdd
    hadd $myDirs[$i]/PVValidation_$myDirs[$i].root $myDirs[$i]/*.root
    mv $myDirs[$i]/PVValidation_$myDirs[$i].root .
end

