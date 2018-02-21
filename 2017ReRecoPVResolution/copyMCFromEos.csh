#!/bin/tcsh

mkdir $PWD/MCfiles

foreach file (`eos ls /eos/cms/store/group/alca_trackeralign/musich/test_out/PVResolutions2017MC/`)
    xrdcp root://eoscms.cern.ch//eos/cms/store/group/alca_trackeralign/musich/test_out/PVResolutions2017MC/$file $PWD/MCfiles/
end
