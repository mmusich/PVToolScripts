#!/bin/tcsh
#set remote=/store/group/alca_trackeralign/musich/test_out/2017_dataReRecoTest2017ReReco/
set remote=/store/group/alca_trackeralign/musich/test_out/2017_dataReReco_bisTest2017ReReco
foreach file (`eos ls $remote`)
    xrdcp root://eoscms.cern.ch//eos/cms/$remote/$file .
end
