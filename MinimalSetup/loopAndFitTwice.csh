#!/bin/tcsh

set COUNT=0
set ObjName=$1 
set theLabel=$2
set ObjName2=$3
set theLabel2=$4

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
     #set filetostudy2=`echo $filetostudy |sed 's/$theLabel/$theLabel2/g'`
     rm -f FittedDeltaZ.tx
	
     echo "$theLabel $theLabel2"
     echo "$filetostudy $filetostudy2"		
     
     root -b -q  $PWD/FitPVResiduals.C++\(\"${PWD}/$ObjName/$filetostudy=${theLabel}\,${PWD}/$ObjName2/$filetostudy2=${theLabel2}\"\,true\,false\,\"$theDate\"\)

     @ COUNT+=1
     #if ($COUNT == 1) then
     #	   break
     #endif	    
end

mv summary.txt summary_${ObjName}.txt

set TargetOutName = "Plots_${ObjName}_vs_${ObjName2}"

mkdir ./${TargetOutName}
mkdir ./${TargetOutName}/Biases/
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
mkdir ./${TargetOutName}/dxyVsPt
mkdir ./${TargetOutName}/dzVsPt
mkdir ./${TargetOutName}/dxyVsPtNorm
mkdir ./${TargetOutName}/dzVsPtNorm   

mv BiasesCanvas*     ./${TargetOutName}/Biases/
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
mv dxyPtTrendNorm*   ./${TargetOutName}/dxyVsPtNorm
mv dzPtTrendNorm*    ./${TargetOutName}/dzVsPtNorm
mv dxyPtTrend*       ./${TargetOutName}/dxyVsPt
mv dzPtTrend*        ./${TargetOutName}/dzVsPt

cp index.php ./${TargetOutName}/Biases/
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
cp index.php ./${TargetOutName}/dxyVsPt      
cp index.php ./${TargetOutName}/dzVsPt	      
cp index.php ./${TargetOutName}/dxyVsPtNorm  
cp index.php ./${TargetOutName}/dzVsPtNorm   
