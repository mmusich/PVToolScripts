#!/bin/tcsh

foreach foldertostudy (`ls ${PWD} | grep DATA`)
    echo $foldertostudy
    cd /tmp/musich/$foldertostudy
    cp /tmp/musich/loopAndFit.csh .
    cp /tmp/musich/PlotPVValidation_forDevel.C .
    set namebase=`echo $foldertostudy |awk '{split($0,a,"_"); print a[2]}'`
    set label=`echo $namebase | sed "s:HIRun2011::g"`
    echo $namebase $label 
    if($label=="") then
	set label="plain"
    endif 
    ./loopAndFit.csh $namebase $label 
end
