#!/bin/tcsh

set TargetOutName = "Plots_Misalignments_v4"

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
mkdir ./${TargetOutName}/Maps
mkdir ./${TargetOutName}/Maps/Abs
mkdir ./${TargetOutName}/Maps/Norm

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
mv d*Abs*            ./${TargetOutName}/Maps/Abs
mv d*Norm*           ./${TargetOutName}/Maps/Norm

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
cp index.php ./${TargetOutName}/Maps/Abs
cp index.php ./${TargetOutName}/Maps/Norm
