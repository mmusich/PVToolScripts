#!/bin/tcsh

set mytarget = $1

mkdir  /afs/cern.ch/user/m/musich/www/PVValidation_2015/run2015CValid
mkdir /afs/cern.ch/user/m/musich/www/PVValidation_2015/run2015CValid/$mytarget
mkdir /afs/cern.ch/user/m/musich/www/PVValidation_2015/run2015CValid/$mytarget/Biases
mkdir /afs/cern.ch/user/m/musich/www/PVValidation_2015/run2015CValid/$mytarget/Normalized
mkdir /afs/cern.ch/user/m/musich/www/PVValidation_2015/run2015CValid/$mytarget/Absolute

cp /afs/cern.ch/user/m/musich/www/PVValidation_2015/candidates_for25ns/V3/Biases/index.php /afs/cern.ch/user/m/musich/www/PVValidation_2015/run2015CValid/$mytarget/Biases
cp /afs/cern.ch/user/m/musich/www/PVValidation_2015/candidates_for25ns/V3/Biases/index.php /afs/cern.ch/user/m/musich/www/PVValidation_2015/run2015CValid/$mytarget/Normalized
cp /afs/cern.ch/user/m/musich/www/PVValidation_2015/candidates_for25ns/V3/Biases/index.php /afs/cern.ch/user/m/musich/www/PVValidation_2015/run2015CValid/$mytarget/Absolute

mv d*Bias*.p* /afs/cern.ch/user/m/musich/www/PVValidation_2015/run2015CValid/$mytarget/Biases
mv d*Norm*.p* /afs/cern.ch/user/m/musich/www/PVValidation_2015/run2015CValid/$mytarget/Normalized
mv d*.p*      /afs/cern.ch/user/m/musich/www/PVValidation_2015/run2015CValid/$mytarget/Absolute

