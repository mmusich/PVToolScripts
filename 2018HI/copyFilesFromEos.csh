#!/bin/tcsh
set remote=/store/group/alca_trackeralign/musich/test_out/HeavyIon2018CheckHI2018
#set remote=/store/group/alca_trackeralign/musich/test_out/2017_dataReRecoTest2017ReReco/
#set remote=/store/group/alca_trackeralign/musich/test_out/2017_dataReReco_bisTest2017ReReco

set myDirs=(PromptGT hp2803)

## create the output folders
foreach i (`seq $#myDirs`)
    if (! -d ./$myDirs[$i] ) then
	mkdir ./$myDirs[$i]
    else
	echo "$myDirs[$i] already exists!"			    
    endif
end

foreach file (`eos ls $remote`)
    foreach i (`seq $#myDirs`)
	if($file =~ *"$myDirs[$i]"*) then 
	    if(-f ./$myDirs[$i]/$file) then 
		echo "$file already exists in $myDirs[$i] ---> exiting"
	    else 
		echo "copying file $file to $myDirs[$i]" 
		xrdcp root://eoscms.cern.ch//eos/cms/$remote/$file ./$myDirs[$i]
	    endif
	endif
    end
end
