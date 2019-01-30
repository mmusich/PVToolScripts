#!/bin/tcsh

#################################################
#0 run: 314882-316757 
#1 run: 316758-317526 
#2 run: 317527-317660 
#3 run: 317661-317663 
#4 run: 317664-318226 
#5 run: 318227-320376 
#6 run: 320377-321830 
#7 run: 321831-322509 
#8 run: 322510-322602 
#9 run: 322603-323231 
#10 run:323232-324244 
#11 run:324245-999999 
################################################# 

set folder=$1
set IOVbegin=(314882 316758 317527 317661 317664 318227 320377 321831 322510 322603 323232 324245)
set IOVend  =(316757 317526 317660 317663 318226 320376 321830 322509 322602 323231 324244 999999)

foreach i (`seq $#IOVbegin`)

    set listToAdd=''

    foreach file (`ls $folder`)

	set date=`echo $file | awk '{split($0,a,"_"); print a[3]}'`
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
