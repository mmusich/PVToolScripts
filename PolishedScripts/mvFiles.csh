#!/bin/tcsh

set mytarget = $1

mkdir /afs/cern.ch/user/m/musich/www/PVValidation_2015/candidates_for25ns/$mytarget
mkdir /afs/cern.ch/user/m/musich/www/PVValidation_2015/candidates_for25ns/$mytarget/Biases
mkdir /afs/cern.ch/user/m/musich/www/PVValidation_2015/candidates_for25ns/$mytarget/Normalized
mkdir /afs/cern.ch/user/m/musich/www/PVValidation_2015/candidates_for25ns/$mytarget/Absolute

cp /afs/cern.ch/user/m/musich/www/PVValidation_2015/candidates_for25ns/V3/Biases/index.php /afs/cern.ch/user/m/musich/www/PVValidation_2015/candidates_for25ns/$mytarget/Biases
cp /afs/cern.ch/user/m/musich/www/PVValidation_2015/candidates_for25ns/V3/Biases/index.php /afs/cern.ch/user/m/musich/www/PVValidation_2015/candidates_for25ns/$mytarget/Normalized
cp /afs/cern.ch/user/m/musich/www/PVValidation_2015/candidates_for25ns/V3/Biases/index.php /afs/cern.ch/user/m/musich/www/PVValidation_2015/candidates_for25ns/$mytarget/Absolute

mv d*Bias*.p* /afs/cern.ch/user/m/musich/www/PVValidation_2015/candidates_for25ns/$mytarget/Biases
mv d*Norm*.p* /afs/cern.ch/user/m/musich/www/PVValidation_2015/candidates_for25ns/$mytarget/Normalized
mv d*.p*      /afs/cern.ch/user/m/musich/www/PVValidation_2015/candidates_for25ns/$mytarget/Absolute

