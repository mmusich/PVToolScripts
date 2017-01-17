
import csv     # imports the csv module
import sys      # imports the sys module

# how to get the lumi
# export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$PATH
# brilcalc lumi -b 'STABLE BEAMS' -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt -u /pb --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -o runs.csv

IOVs=[272007,273725,274244,274335,274420,274442,274998,275068,275309,275344,275376,275783,275847,275920,276242,276283,276363,276454,276501,276525,276581,276653,276776,276811,276870,276950,277076,277096,277168,277218,278018,278240,278315,278406,278803,278820,278957,278975,279588,279667,279715,279766,279841,279931,280015,280191,280330,280385,281693,281797,282037,282731,282807,283041,283270,283308,283408,283553,283830,283877,283934,284029,284044,300000]

d = {}
d['run:fill'] = []
d['time'] = []
d['nls'] = []
d['ncms'] = []
d['delivered'] = []
d['recorded'] = []

dictReader = csv.DictReader(open('runs.csv', 'rb'), fieldnames = ['run:fill','time','nls','ncms','delivered','recorded'], delimiter = ',', quotechar = '"')

for row in dictReader:
    #print row
    if row['run:fill'].startswith("#"):
        continue
    for key in row:
        d[key].append(row[key])

#print d['nls']

out = {}
out['IOV'] = []
out['LUMI'] = []

for i in range(len(IOVs)-1):
    sum=0
    out['IOV'].append(str(IOVs[i])+"-"+str(IOVs[i+1]))
    for j,entry in enumerate(d['run:fill']):
        sep = ':'
        run = int(entry.split(sep, 1)[0])
        if( (IOVs[i] <= int(run)) and (IOVs[i+1] > int(run))):
            print run," in IOV:",IOVs[i],"-",IOVs[i+1]," LUMI: ", d['recorded'][j]," /pb"
            sum = sum + float(d['recorded'][j])
    out['LUMI'].append(sum)
    
grandSum=0

print '#################################################'
for k,element in enumerate(out['IOV']):
    grandSum=grandSum+out['LUMI'][k]
    print k,"run:",element,"lumi:",out['LUMI'][k],"/pb | grand sum: ",grandSum/1000.

print '#################################################'
print "Grand total: ",grandSum
