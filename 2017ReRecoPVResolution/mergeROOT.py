#!/usr/bin/env python

import os
from os import listdir
from os.path import isfile, join
from collections import defaultdict

onlyfiles = [f for f in listdir("./") if isfile(join("./", f))]

iovs= [294034,296641,297179,297281,298653,299277,299443,300389,301046,302131,303790,304911,305898]

#print onlyfiles                                                                                           

dictionary=defaultdict(list)
for name in onlyfiles:
    if not name.endswith('.root'):
        continue
    file_info=name.split("_")
    theRun=file_info[2].split(".")[0]
    theAlignment=file_info[1]
    theIOVcount=0
    theIOV=0
    for iov in iovs:
        theIOVcount+=1
        if(int(theRun)>int(iov)):
            theIOV=theIOVcount

    dictionary[(theIOV,theAlignment)].append(name)
    #dictionary[name]={'run':theRun,'alignment':theAlignment,'iov:':theIOV}


for key, value in dictionary.iteritems():
    #print key,value 

    command = 'hadd pvresolution_'+str(key[1])+'_IOV'+str(key[0])+'.root'
    for f in value:
        command+=(' '+f)

    print command
    os.system(command)
