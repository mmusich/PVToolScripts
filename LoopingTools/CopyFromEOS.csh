#!/bin/csh

foreach subfolder (`cmsLs /store/caf/user/musich/Alignment/PVValidation/`)
    #echo $subfolder
    if ("$subfolder" =~ *"7TeVDataOldAPE"*) then
	set namebase=`echo $subfolder |awk '{split($0,b,"/"); print b[8]}'`
 	echo "opening:" $subfolder"/"$namebase".root" 
	cmsStage ${subfolder}/${namebase}.root .
    endif
end
