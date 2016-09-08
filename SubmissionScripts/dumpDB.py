#!/usr/bin/env python

import os
import pickle
from pprint import pprint
from prettytable import PrettyTable

##############################################
def dumpDB():
##############################################
    dbName = "runInfo.pkl"
    infos = {}

    listOfAcceptedDetails=["1.run",            
                           #"1.1runevents",
                           "2.conf",           
                           "3.gt",             
                           #"4.allFromGT",      
                           #"5.alignmentDB",    
                           "6.alignmentTag",   
                           #"7.apeDB",          
                           "8.apeTag",         
                           #"9.applyBows",      
                           #"9.1bowDB",         
                           "9.2bowTag",        
                           #"9.3ptCut",         
                           #"9.4lumilist",      
                           #"9.5applyEXTRACOND",
                           #"9.6conditions",    
                           #"9.7nfiles",        
                           #"9.8srcFiles",
                           ]

    print listOfAcceptedDetails

    if os.path.exists(dbName):
        with open(dbName) as f:
            infos = pickle.load(f)
            pprint(infos)
            headers=[]
            for key in infos:
                for detail in infos[key]:
                    #print key,infos[key][detail]
                    if key == 11:
                        if(detail in listOfAcceptedDetails):
                            headers.append(detail)
                                        
            t = PrettyTable(sorted(headers))
                        
            for key in infos:
                contents=[]
                for detail in sorted(headers):
                    contents.append(infos[key][detail])
                t.add_row(contents)

            print t

def main():
    dumpDB()
    
if __name__ == "__main__":        
    main()
