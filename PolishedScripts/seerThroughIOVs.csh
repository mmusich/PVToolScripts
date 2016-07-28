#!/bin/tcsh

#################################################
#0 run:  275657-275657 lumi: 20.104 /pb
#1 run:  275658-275765 lumi: 70.733 /pb
#2 run:  275766-275771 lumi: 0.722 /pb
#3 run:  275772-275777 lumi: 141.197 /pb
#4 run:  275778-275827 lumi: 149.324 /pb
#5 run:  275828-275828 lumi: 0 /pb
#6 run:  275829-275836 lumi: 366.658 /pb
#7 run:  275837-275840 lumi: 85.912 /pb
#8 run:  275841-275846 lumi: 0 /pb
#9 run:  275847-275886 lumi: 213.061 /pb
#10 run: 275887-275920 lumi: 561.873 /pb
#11 run: 275921-276241 lumi: 163.333 /pb
#12 run: 276242-276242 lumi: 275.451 /pb
#13 run: 276243-276243 lumi: 80.003 /pb
#14 run: 276244-276281 lumi: 136.344 /pb
#15 run: 276282-276282 lumi: 210.215 /pb
#16 run: 276283-276314 lumi: 171.192 /pb
#17 run: 276315-276316 lumi: 34.705 /pb
#18 run: 276317-276326 lumi: 126.045 /pb
#19 run: 276327-276351 lumi: 0 /pb
#20 run: 276352-276360 lumi: 6.31 /pb
#21 run: 276361-276362 lumi: 151.659 /pb
#22 run: 276363-276383 lumi: 208.556 /pb
#23 run: 276384-276452 lumi: 495.853 /pb
#24 run: 276453-276453 lumi: 0 /pb
#25 run: 276454-276652 lumi: 1945.69 /pb
#26 run: 276653-300000 lumi: 1360.772 /pb
################################################# 

set folder=$1
set IOVbegin=(275657 275658 275766 275772 275778 275828 275829 275837 275841 275847 275887 275921 276242 276243 276244 276282 276283 276315 276317 276327 276352 276361 276363 276384 276453 276454 276653)
set IOVend  =(275657 275765 275771 275777 275827 275828 275836 275840 275846 275886 275920 276241 276242 276243 276281 276282 276314 276316 276326 276351 276360 276362 276383 276452 276453 276652 300000)

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
