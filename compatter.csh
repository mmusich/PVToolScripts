#!/bin/tcsh

set myObjects = (pcl prompt startup)

foreach i (`seq $#myObjects`)
    mkdir ./$myObjects[$i]
end

set sourceDir=/eos/cms/store/caf/user/mwassmer/AlignmentValidation/STARTUP_2018/
foreach dir (`eos ls $sourceDir`)
    set date=`echo $dir | awk '{split($0,a,"_"); print a[2]}'`
    foreach file (`eos ls $sourceDir/$dir | grep root`)
	foreach i (`seq $#myObjects`)
	    if($file =~ *"$myObjects[$i]"*) then
		echo $file $date $myObjects[$i]
		eos cp  $sourceDir/$dir/$file  ./$myObjects[$i]/PVValidation_$myObjects[$i]_$date.root
	    endif		
	end
    end
end
