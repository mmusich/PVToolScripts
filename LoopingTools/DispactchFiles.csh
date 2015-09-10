#!/bin/tcsh

foreach filetostudy (`ls ${PWD} | grep .root`)
     echo $filetostudy
     set namebase=`echo $filetostudy |awk '{split($0,a,"_"); print a[2]}'`
     echo "Studying $namebase file"
     if(! -d ${PWD}/DATA_$namebase) then
	mkdir ${PWD}/DATA_$namebase
	mv $filetostudy ${PWD}/DATA_$namebase
     else 
	mv $filetostudy ${PWD}/DATA_$namebase
     endif
end






