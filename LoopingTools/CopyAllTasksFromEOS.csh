#!/bin/tcsh

set COUNT=0
set REM=0
set myDir=$1

if (! -d copyReport.txt ) then
    touch copyReport.txt 
endif

foreach folder (`cmsLs /store/caf/user/musich/Alignment/PVValidation/`) 
    if("$folder" =~ *"store"*) then
	echo ">>>>>>>>>>>>>>>>>>>>>>>>>>> $folder"
	set foldernamebase=`echo $folder |awk '{split($0,b,"/"); print b[8]}'`
	#echo $foldernamebase
	set COUNTinFolder=0
	set COUNTremoved=0
	foreach inputfile (`cmsLs /store/caf/user/musich/Alignment/PVValidation/$foldernamebase`)
	    #echo $inputfile
	    set filenamebase=`echo $inputfile |awk '{split($0,b,"/"); print b[9]}'`
	    if ("$filenamebase" =~ *"$foldernamebase"*) then
		set namebase=`echo $filenamebase |awk '{split($0,b,"_"); print b[2]}'` 
		if ("$namebase" != "$foldernamebase") then
		    echo "$namebase does not match with $foldernamebase!"
		    echo "removing: /store/caf/user/musich/Alignment/PVValidation/$foldernamebase/$filenamebase"
		    cmsRm /store/caf/user/musich/Alignment/PVValidation/$foldernamebase/$filenamebase
		    @ COUNTremoved+=1
		    @ REM+=1
		else
		    echo "copying: /store/caf/user/musich/Alignment/PVValidation/$foldernamebase/$filenamebase" 
		    cmsStage /store/caf/user/musich/Alignment/PVValidation/$foldernamebase/$filenamebase .
		    if (! -d evts.log ) then
			rm -f evts.log 
		    endif
		    root -b -q $PWD/check.C\(\"${PWD}/$filenamebase\"\) > evts.log
		    set word=`cat evts.log | grep Long64 | awk '{split($0,b,")"); print b[2]}'`
		    if ($word == "0") then
			rm ${PWD}/$filenamebase
			echo "0 events: removing ${PWD}/$filenamebase"
		    endif
		    @ COUNT+=1
		    @ COUNTinFolder+=1
		endif
	    endif  
	end
	echo "$foldernamebase $COUNTinFolder files copied! $COUNTremoved files removed" >> copyReport.txt
    endif
end
echo "copied $COUNT files; removed $REM files"  
