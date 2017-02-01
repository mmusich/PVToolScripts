#!/bin/tcsh

set COUNT=0
set myDir=$1

set myObjects = (PromptGT Sept2016ReReco EOYReReco EOYReRecoSingleIOVAPE)

foreach i (`seq $#myObjects`) 
    echo $myObjects[$i]
    mkdir ./$myObjects[$i]
    foreach folder (`cmsLs /store/group/alca_trackeralign/musich/test_out | grep EOY2016Final`)
    	foreach inputfile (`cmsLs /store/group/alca_trackeralign/musich/test_out/${folder} | grep root`)
    	    if($inputfile =~ *"$myObjects[$i]"*) then
    		#echo $inputfile $myObjects[$i]
    		echo "copying: /store/group/alca_trackeralign/musich/test_out/${folder}/${inputfile} $myObjects[$i]" 
    	        cmsStage /store/group/alca_trackeralign/musich/test_out/${folder}/$inputfile ./$myObjects[$i]
    		@ COUNT+=1
    	    endif
        end
    end 
end

echo "copied $COUNT files"  
