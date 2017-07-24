#!/bin/csh

set targetDir=$HOME/www/PVValidation_2017/test2

if (! -d $targetDir ) then
    mkdir $targetDir
else
    echo "$targetDir already exists; I am not going to re-create it"
endif
cd $targetDir
rm *.php 
rm *.pdf 
rm *.png
wget https://raw.githubusercontent.com/mmusich/PVToolScripts/master/PolishedScripts/index.php 
cd -
cp -pr *.png *.pdf *.php $targetDir
