#!/bin/tcsh

set COUNT=0
set ObjName=$1 
set theLabel=$2
set ObjName2=$3
set theLabel2=$4
set ObjName3=$5
set theLabel3=$6
set ObjName4=$7
set theLabel4=$8


setenv CMSSW_DIR ${CMSSW_BASE}/src/Alignment/OfflineValidation/test/PVValResults
#cp ${CMSSW_DIR}/PlotPVValidation_forDevel.C .
#cp ${CMSSW_DIR}/check.C .

foreach filetostudy (`ls ${PWD}/$ObjName | grep .root`)

     #set namebase=`echo $filetostudy |awk '{split($0,a,"_"); print a[2]}'`
     #set datebase=`echo $filetostudy |awk '{split($0,b,"_"); print b[]}'`
     #set words = `echo $filetostudy:q | sed 's/_/ /g'`

     # approach to ge the last item in the list
     
     set datebase=`echo $filetostudy |awk '{n=split($0,a,"_"); print a[n]}'`
     set theDate=`echo $datebase |awk '{split($0,c,"."); print c[1]}'`
     
     #echo "$datebase $theDate"

     set filetostudy2=`echo $filetostudy| sed "s/$ObjName/$ObjName2/g" `
     set filetostudy3=`echo $filetostudy| sed "s/$ObjName/$ObjName3/g" `
     set filetostudy4=`echo $filetostudy| sed "s/$ObjName/$ObjName4/g" `
     #set filetostudy2=`echo $filetostudy |sed 's/$theLabel/$theLabel2/g'`
     rm -f FittedDeltaZ.tx
	
     echo "$theLabel $theLabel2 $theLabel3 $theLabel4"
     echo "$filetostudy $filetostudy2 $filetostudy3 $filetostudy4"		
     
     #root -b -q  $PWD/FitPVResiduals_forLoop.C++\(\"${PWD}/$ObjName/$filetostudy=${theLabel}\,${PWD}/$ObjName2/$filetostudy2=${theLabel2}\"\,true\,false\,\"$theDate\"\)

     root -b -q runMe.C\(\"${PWD}/$ObjName/$filetostudy=${theLabel}\,${PWD}/$ObjName2/$filetostudy2=${theLabel2}\,${PWD}/$ObjName3/$filetostudy3=${theLabel3}\,${PWD}/$ObjName4/$filetostudy4=${theLabel4}\"\,\"$theDate\"\)

     rm -f numevents.out 
     root -b -q  $PWD/check3.C\(\"${PWD}/$ObjName/$filetostudy\"\) > numevents.out
     set rawnumevents=`tail -1 numevents.out`
     set numevents=`echo $rawnumevents | awk '{split($0,a,")"); print a[2]}'` 
    
     set rawfits=`tail -1 FittedDeltaZ.txt`
     set deltaz=`echo $rawfits | awk '{split($0,a,"|"); print a[1]}'`
     set sigmadeltaz=`echo $rawfits | awk '{split($0,a,"|"); print a[2]}'` 
 
     setenv flag "BAD"

     if (! -d summary.txt ) then
     	 touch summary.txt 
     endif
    
     if (${numevents} > 1000) then
     	setenv flag "GOOD"
     endif

     echo "- File $theDate had ${numevents} events. Fit separation: $deltaz \pm $sigmadeltaz"

     echo $theDate $deltaz $sigmadeltaz ${numevents} $flag >> summary.txt

     @ COUNT+=1
     #if ($COUNT == 1) then
     #	   break
     #endif	    
end

mv summary.txt summary_${ObjName}.txt

set TargetOutName = "Plots_${ObjName}_vs_${ObjName2}"

mkdir ./${TargetOutName}
mkdir ./${TargetOutName}/Biases/
mkdir ./${TargetOutName}/Pulls/
mkdir ./${TargetOutName}/Biases/dzPhi
mkdir ./${TargetOutName}/Biases/dxyPhi
mkdir ./${TargetOutName}/Biases/dzEta
mkdir ./${TargetOutName}/Biases/dxyEta
mkdir ./${TargetOutName}/Fit
mkdir ./${TargetOutName}/dxyVsEta
mkdir ./${TargetOutName}/dzVsEta
mkdir ./${TargetOutName}/dxyVsPhi
mkdir ./${TargetOutName}/dzVsPhi
mkdir ./${TargetOutName}/dxyVsEtaNorm
mkdir ./${TargetOutName}/dzVsEtaNorm
mkdir ./${TargetOutName}/dxyVsPhiNorm
mkdir ./${TargetOutName}/dzVsPhiNorm

mv BiasesCanvas*     ./${TargetOutName}/Biases/
mv PullCanvas*       ./${TargetOutName}/Pulls/
mv dzPhiBiasCanvas*  ./${TargetOutName}/Biases/dzPhi
mv dxyPhiBiasCanvas* ./${TargetOutName}/Biases/dxyPhi
mv dzEtaBiasCanvas*  ./${TargetOutName}/Biases/dzEta
mv dxyEtaBiasCanvas* ./${TargetOutName}/Biases/dxyEta
mv dzPhiTrendFit*    ./${TargetOutName}/Fit
mv dxyEtaTrendNorm*  ./${TargetOutName}/dxyVsEtaNorm
mv dzEtaTrendNorm*   ./${TargetOutName}/dzVsEtaNorm
mv dxyPhiTrendNorm*  ./${TargetOutName}/dxyVsPhiNorm
mv dzPhiTrendNorm*   ./${TargetOutName}/dzVsPhiNorm
mv dxyEtaTrend*      ./${TargetOutName}/dxyVsEta
mv dzEtaTrend*       ./${TargetOutName}/dzVsEta
mv dxyPhiTrend*      ./${TargetOutName}/dxyVsPhi
mv dzPhiTrend*       ./${TargetOutName}/dzVsPhi

cp index.php ./${TargetOutName}/Biases/
cp index.php ./${TargetOutName}/Pulls/
cp index.php ./${TargetOutName}/Biases/dzPhi
cp index.php ./${TargetOutName}/Biases/dxyPhi
cp index.php ./${TargetOutName}/Biases/dzEta
cp index.php ./${TargetOutName}/Biases/dxyEta
cp index.php ./${TargetOutName}/Fit
cp index.php ./${TargetOutName}/dxyVsEta
cp index.php ./${TargetOutName}/dzVsEta
cp index.php ./${TargetOutName}/dxyVsPhi
cp index.php ./${TargetOutName}/dzVsPhi
cp index.php ./${TargetOutName}/dxyVsEtaNorm
cp index.php ./${TargetOutName}/dzVsEtaNorm
cp index.php ./${TargetOutName}/dxyVsPhiNorm
cp index.php ./${TargetOutName}/dzVsPhiNorm




