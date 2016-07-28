#!/bin/tcsh

set COUNT=0
set myDirs = (test_2016_07_12_12_40_35_DATA_Prompt_Run2016B-v2)

foreach i (`seq $#myDirs`)
    echo "Copying files from $myDirs[$i]"
    nsmkdir /castor/cern.ch/cms/archive/user/m/musich/FromEOS/test_out/$myDirs[$i]
    foreach file (`eos ls /store/caf/user/musich/test_out/$myDirs[$i]`)
	echo "moving $myDirs[$i]/$file to castor"
	cmsStage /store/caf/user/musich/test_out/$myDirs[$i]/$file .
	xrdcp $file "root://castorcms.cern.ch//castor/cern.ch/cms/archive/user/m/musich/FromEOS/test_out/$myDirs[$i]/$file?svcClass=archive"
        @ COUNT+=1
	rm $file
    end
end

echo "Total checked files: $COUNT"
nsls -l /castor/cern.ch/cms/archive/user/m/musich/FromEOS/test_out

