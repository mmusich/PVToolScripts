#!/bin/tcsh

set COUNT=0
set myDir=$1

set myObjects = (PromptGT Sept2016ReReco EOYReReco EOYReRecoSingleIOVAPE)

foreach folder (`cmsLs /store/group/alca_trackeralign/musich/test_out | grep EOY2016Final`)
    foreach inputfile (`cmsLs /store/group/alca_trackeralign/musich/test_out/${folder} | grep root`)
	#echo $inputfile $myObjects[$i]
	echo "removing: /store/group/alca_trackeralign/musich/test_out/${folder}/${inputfile}" 
	eos rm /store/group/alca_trackeralign/musich/test_out/${folder}/$inputfile
    	@ COUNT+=1	
    end 
    eos rm -r /store/group/alca_trackeralign/musich/test_out/${folder}
end

echo "removed $COUNT files"  
