#!/bin/tcsh

#set myDirs =  (pp2015B-PCLupdate pp2015B-PCLupdateAndTemplates pp2015B-Prompt pp2015B-PromptAndTemplates pp2015B-Run1Ali pp2015B-Run2Ali0TCollisions pp2015B-Run2AliCRUZET pp2015B-Run2AliCosmics pp2015B-hp1368 pp2015B-hp1370 pp2015B-hp1375 pp2015B-hp1387 pp2015B-hp1388 pp2015B-mp1799)

#set myDirs = (test_2015_10_16_17_31_02_DATA_pp2015B-OfflineGT   test_2015_10_17_15_43_27_DATA_pp2015B-PromptGT  test_2015_10_17_19_22_09_DATA_pp2015C-OfflineGT  test_2015_10_17_19_47_17_DATA_pp2015C-PromptGT  test_2015_10_19_08_57_24_DATA_pp2015D-PromptGT test_2015_10_19_09_06_14_DATA_pp2015D-OffineGT   test_2015_10_20_09_43_14_DATA_pp2015Dv4-PromptGT test_2015_10_20_10_06_31_DATA_pp2015Dv4-OfflineGT  test_2015_10_22_15_56_00_DATA_pp2015Dv4-testPCL_v2 test_2015_10_22_16_16_35_DATA_pp2015D-testPCL_v2)

#set myDirs = (test_2015_12_13_16_04_39_DATA_pp2015Dv3-0T_EOY test_2015_12_13_16_09_17_DATA_pp2015Dv3-0T_v16 test_2015_12_13_16_10_24_DATA_pp2015Dv4-0T_v16 test_2015_12_13_16_15_33_DATA_pp2015Dv4-0T_EOY test_2015_12_13_16_31_55_DATA_pp2015Cv3-0T_EOY test_2015_12_13_16_47_29_DATA_pp2015Cv3-0T_v16 test_2015_12_13_16_53_07_DATA_pp2015Cv2-0T_v16 test_2015_12_13_17_33_11_DATA_pp2015Cv2-0T_EOY)

set myDirs = (test_2016_05_27_11_48_20_DATA_Prompt_May2016_v2 test_2016_06_03_17_22_42_DATA_Prompt_May2016-v2 test_2016_06_14_12_28_31_DATA_Prompt_May2016_fullML)

touch report.txt
set COUNT=0
set COUNTREMOVED=0
set COUNTbyFolder=0

foreach i (`seq $#myDirs`)
    echo "=================================================" 
    echo "Reading from $myDirs[$i]"
    set COUNTbyFolder=0
    set theOutDir=`echo $myDirs[$i] | awk '{split($0,a,"_"); print a[9]}'`
    echo $theOutDir
    #mkdir ./$myDirs[$i]
    mkdir ./$theOutDir
    foreach inputfile (`eos ls /store/caf/user/musich/test_out/$myDirs[$i]`)
	echo "Analyzing: /store/caf/user/musich/test_out/$myDirs[$i]/$inputfile" 	
	if ($inputfile =~ *"root"*) then
	    @ COUNT+=1
	    @ COUNTbyFolder+=1
	    set myLFN = root://eoscms//eos/cms/store/caf/user/musich/test_out/$myDirs[$i]/$inputfile
	    root -q -b ./check3.C\(\"$myLFN\"\)  > numevents.out
	    #alias MATH 'set \!:1 = `echo "scale=5;\!:5-$" | bc -l`' 
	    set numevents=`tail -5 numevents.out | grep events | awk '{split($0,a," "); print a[4]}'`
	    set numTracks=`tail -5 numevents.out | grep tracks | awk '{split($0,b," "); print b[4]}'`
	    set trkeffic=`tail -5 numevents.out | grep trks    | awk '{split($0,c,":"); print c[2]}'`
	    echo "ALCARECO events: $numevents"
	    echo "ALCARECO tracks: $numTracks"
	    if ("$numevents" =~ *"+0"*)  then
		set rootNevents=99999
	    else
		set rootNevents=$numevents
	    endif
	    echo $rootNevents
	    if($rootNevents>0) then 
		echo "trk eff (ALCARECO): $trkeffic"
		cmsStage /store/caf/user/musich/test_out/$myDirs[$i]/$inputfile ./$theOutDir #./$myDirs[$i]
	    else 
		echo "0 events!!!! => going to remove the file:  /store/caf/user/musich/test_out/$myDirs[$i]/$inputfile "
		@ COUNTREMOVED+=1
		eos rm /store/caf/user/musich/test_out/$myDirs[$i]/$inputfile
	    endif
	    echo "-------------------------------------------------" 
	    #echo $inputfile $numevents $numTracks $trkeffic >> report.txt
	    printf '%-5s %-5s %-20s %-20s %-20s\n' $COUNTbyFolder $inputfile $numevents $numTracks >> report_$myDirs[$i].txt
	endif
    end
end

echo "Removed $COUNTREMOVED files" 
echo "Total checked files: $COUNT"
