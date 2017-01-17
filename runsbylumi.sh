#!/bin/sh

### command to be given: 
#$ ./runsbylumi.sh 272007 284044 500

startRun=$1
endRun=$2
lumiLimit=$3
jsonPath=$4

if [ "$startRun" == "" ] || [ "$endRun" == "" ] || [ "$lumiLimit" == "" ]; then
  echo "USAGE: ./runsbylumi.sh <startRun> <endRun> <lumiBlock (pb-1)> <json_file:optional>"
  exit 1
fi

export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$PATH

echo "Evaluating the list of runs in range [$startRun-$endRun] divided by luminosity of $lumiLimit pb-1"
if [ "$jsonPath" != "" ]; then
  echo "Using JSON file: $jsonPath"
fi

# prototype search string 
#| 284044:5451 | 10/26/16 20:46:05 | 40 | 40 | 5288300.215 | 4912825.747 |

searchString="| *[0-9]\+:[0-9]\+ *|"
if [ "$jsonPath" == "" ];then 
  # brilcalc lumi -b 'STABLE BEAMS' -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-274240_13TeV_PromptReco_Collisions16_JSON.txt -u /pb --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json
  brilcalc lumi -b 'STABLE BEAMS' --begin $startRun --end $endRun  >  ./runsbylumi.log
  else
  brilcalc lumi -b 'STABLE BEAMS' -i $jsonPath --begin $startRun --end $endRun  >  ./runsbylumi.log
fi
nRuns=`cat ./runsbylumi.log | grep -e "| *[0-9]\+:[0-9]\+ *|" | wc -l`
echo "brilcalc finished with total of $nRuns runs"
echo "Selecting run ranges by luminosity"

totalLumi=0.0
currentLumi=0.0
runList=""
nRuns=0

while read line
do
  #echo $line  
  string=`grep -e "$searchString" <<< $line`
  #echo $string
  if [ "$string" == "" ]; then
    continue
  fi
  #echo $line
  run1=`echo $string |awk '{split($0,a,"|"); print a[2]}'`
  run=`echo $run1 |awk '{split($0,b,":"); print b[1]}'`
  #echo "Run: $run"
  #echo $string
  lumi=`echo $string |awk '{split($0,a,"|"); print a[7]}'`
  unit="ub"
  #echo "Lumi: $lumi"

  if [ "$unit" == "fb" ]; then
    scale=1000
  elif [ "$unit" == "pb" ]; then
    scale=1.0
  elif [ "$unit" == "nb" ]; then
    scale=0.001
  elif [ "$unit" == "ub" ]; then
    scale=0.000001
  elif [ "$unit" == "mb" ]; then
    scale=0.000000001
  elif [ "$unit" == "b" ]; then
    scale=0.000000000001
  fi
  
  totalLumi=`echo "$totalLumi + $lumi * $scale" | bc`
  currentLumi=`echo "$currentLumi + $lumi * $scale" | bc`
  compareResult=`echo "$currentLumi >= $lumiLimit" | bc`
  #echo "Lumi: $lumi /$unit Current: $currentLumi"
  if [ $compareResult -gt 0 ]; then
    echo "$run $currentLumi $totalLumi"
    currentLumi=0.0
    runList+="$run,"
    let "nRuns = nRuns + 1"
  fi

done < ./runsbylumi.log

let "nRuns = nRuns + 1"
runList=${runList%,}

echo
echo "Runs: [ $startRun,$runList,$endRun ]"
echo "Number of Run ranges: $nRuns"
echo

echo "Total Luminosity: $totalLumi /pb"
echo "From brilcalc.py:"
cat ./runsbylumi.log | grep -e "." | tail -3

#rm ~/runsbylumi.log
