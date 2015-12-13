#!/bin/tcsh

set folder=$1
set IOVbegin=(252032 254286 254319 254350 255039 256350 257039 259900 260453)
set IOVend=(253974 254318 254349 254636 256349 256514 257363 260338 260531)

foreach i (`seq $#IOVbegin`)

    set listToAdd=''

    foreach file (`ls $folder`)

	set date=`echo $file | awk '{split($0,a,"_"); print a[5]}'`
	set theDate=`echo $date |awk '{split($0,c,"."); print c[1]}'`

	#echo $theDate

	if($theDate>$IOVbegin[$i] && $theDate<$IOVend[$i]) then
	   set listToAdd = "$listToAdd $folder/$file" 
	endif
	
    end

    echo "========================================================="
    echo "IOVbegin: $IOVbegin[$i], IOVend: $IOVend[$i]| $listToAdd"
    hadd PVValidation_${folder}_$IOVbegin[$i]-$IOVend[$i].root $listToAdd
end
