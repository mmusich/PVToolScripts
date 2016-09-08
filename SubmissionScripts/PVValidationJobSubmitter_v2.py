#!/usr/bin/env python
#
# Author:   Marco MUSICH
#
# Usage:
# python PVValidationJobSubmitter_v2.py -i <input ini file> -j <name of task> -D <dataset to use> -r (if run based) -s (if submitted)
#
# exmple:
# python PVValidationJobSubmitter_v2.py -i test_Prompt_2016B_forICHEP.ini -j Prompt_Run2016B-v2 -D /HLTPhysics/Run2016B-TkAlMinBias-PromptReco-v2/ALCARECO -r
#

'''Script that submits CMS Tracker Alignment Primary Vertex Validation workflows
'''

__author__ = 'Marco Musich'
__copyright__ = 'Copyright 2015, CERN CMS'
__credits__ = ['Ernesto Migliore', 'Salvatore Di Guida']
__license__ = 'Unknown'
__maintainer__ = 'Marco Musich'
__email__ = 'marco.musich@cern.ch'
__version__ = 1

import datetime,time
import os,sys
import copy
import string, re
import ConfigParser, json
import subprocess
from optparse import OptionParser
from subprocess import Popen, PIPE
import collections
import multiprocessing
import pickle
from pprint import pprint
from prettytable import PrettyTable
#import Utilities.General.cmssw_das_client as cmssw_das_client

##############################################
def getCommandOutput(command):
##############################################
    """This function executes `command` and returns it output.
    Arguments:
    - `command`: Shell command to be invoked by this function.
    """
    child = os.popen(command)
    data = child.read()
    err = child.close()
    if err:
        print '%s failed w/ exit code %d' % (command, err)
    return data

##############################################
def getFilesForRun(blob):
##############################################
    cmd2 = ' das_client.py --limit=0 --query \'file run='+blob[0]+' dataset='+blob[1]+'\''
    q = Popen(cmd2 , shell=True, stdout=PIPE, stderr=PIPE)
    out, err = q.communicate()
    print cmd2,'\n',out.rstrip('\n')
    outputList = out.split('\n')
    outputList.pop()
    return outputList #,err

##############################################
def getNEvents(run, dataset):
##############################################
    nEvents = subprocess.check_output(["das_client", "--limit", "0", "--query", "summary run={} dataset={} | grep summary.nevents".format(run, dataset)])
    return 0 if nEvents == "[]\n" else int(nEvents)

##############################################
def isInJSON(run,jsonfile):
##############################################
    with open(jsonfile, 'rb') as myJSON:
        jsonDATA = json.load(myJSON)
        return (run in jsonDATA)
    
##############################################
def das_client(query):
##############################################
    """
    Submit `query` to DAS client and handle possible errors.
    Further treatment of the output might be necessary.

    Arguments:
    - `query`: DAS query
    """
    for _ in xrange(3):         # maximum of 3 tries
        das_data = cmssw_das_client.get_data(query, limit = 0)
        if das_data["status"] != "error": break
    if das_data["status"] == "error":
        print_msg("DAS query '{}' failed 3 times. "
                  "The last time for the the following reason:".format(query))
        print das_data["reason"]
        sys.exit(1)
    return das_data["data"]

##############################################
def find_key(collection, key):
##############################################
    """Searches for `key` in `collection` and returns first corresponding value.

    Arguments:
    - `collection`: list of dictionaries
    - `key`: key to be searched for
    """

    for item in collection:
        if item.has_key(key):
            return item[key]
    print collection
    raise KeyError(key)

##############################################
def get_events(entity, name):
##############################################
    """Retrieves the number of events from `entity` called `name`.

    Arguments:
    - `entity`: type of entity
    - `name`: name of entity
    """

    query= "{0:s}={1:s} | grep {0:s}.nevents".format(entity, name)
    data = das_client("{0:s}={1:s} | grep {0:s}.nevents".format(entity, name))
    return int(find_key(find_key(data, entity), "nevents"))

##############################################
def to_bool(value):
##############################################
    """
       Converts 'something' to boolean. Raises exception for invalid formats
           Possible True  values: 1, True, "1", "TRue", "yes", "y", "t"
           Possible False values: 0, False, None, [], {}, "", "0", "faLse", "no", "n", "f", 0.0, ...
    """
    if str(value).lower() in ("yes", "y", "true",  "t", "1"): return True
    if str(value).lower() in ("no",  "n", "false", "f", "0", "0.0", "", "none", "[]", "{}"): return False
    raise Exception('Invalid value for boolean conversion: ' + str(value))

##############################################
def updateDB(run,runInfo):
##############################################
    dbName = "runInfo.pkl"
    infos = {}
    if os.path.exists(dbName):
        with open(dbName) as f:
            infos = pickle.load(f)

    if run not in infos:
        infos[run] = runInfo

    with open(dbName, "wb") as f:
        pickle.dump(infos, f)

##############################################
def dumpDB():
##############################################
    dbName = "runInfo.pkl"
    infos = {}
    if os.path.exists(dbName):
        with open(dbName) as f:
            infos = pickle.load(f)
            #pprint(infos)
            headers=[]
            for key in infos:
                for detail in infos[key]:
                    #print key,infos[key][detail]
                    if key == 11:
                        headers.append(detail)
                       
            t = PrettyTable(headers)
                        
            for key in infos:
                contents=[]
                for detail in infos[key]:
                    contents.append(infos[key][detail])
                t.add_row(contents)

            print t
            
####################--- Classes ---############################
class BetterConfigParser(ConfigParser.ConfigParser):

    ##############################################
    def optionxform(self, optionstr):
        return optionstr

    ##############################################
    def exists( self, section, option):
         try:
             items = self.items(section) 
         except ConfigParser.NoSectionError:
             return False
         for item in items:
             if item[0] == option:
                 return True
         return False

    ##############################################
    def __updateDict( self, dictionary, section ):
        result = dictionary
        try:
            for option in self.options( section ):
                result[option] = self.get( section, option )
            if "local"+section.title() in self.sections():
                for option in self.options( "local"+section.title() ):
                    result[option] = self.get( "local"+section.title(),option )
        except ConfigParser.NoSectionError, section:
            msg = ("%s in configuration files. This section is mandatory."
                   %(str(section).replace(":", "", 1)))
            #raise AllInOneError(msg)
        return result     

    ##############################################
    def getResultingSection( self, section, defaultDict = {}, demandPars = [] ):
        result = copy.deepcopy(defaultDict)
        for option in demandPars:
            try:
                result[option] = self.get( section, option )
            except ConfigParser.NoOptionError, globalSectionError:
                globalSection = str( globalSectionError ).split( "'" )[-2]
                splittedSectionName = section.split( ":" )
                if len( splittedSectionName ) > 1:
                    localSection = ("local"+section.split( ":" )[0].title()+":"
                                    +section.split(":")[1])
                else:
                    localSection = ("local"+section.split( ":" )[0].title())
                if self.has_section( localSection ):
                    try:
                        result[option] = self.get( localSection, option )
                    except ConfigParser.NoOptionError, option:
                        msg = ("%s. This option is mandatory."
                               %(str(option).replace(":", "", 1).replace(
                                   "section",
                                   "section '"+globalSection+"' or", 1)))
                        #raise AllInOneError(msg)
                else:
                    msg = ("%s. This option is mandatory."
                           %(str(globalSectionError).replace(":", "", 1)))
                    #raise AllInOneError(msg)
        result = self.__updateDict( result, section )
        #print result
        return result

##### method to parse the input file ################################
def ConfigSectionMap(config, section):
    the_dict = {}
    options = config.options(section)
    for option in options:
        try:
            the_dict[option] = config.get(section, option)
            if the_dict[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            the_dict[option] = None
    return the_dict

###### method to create recursively directories on EOS #############
def mkdir_eos(out_path):
    print "creating",out_path
    newpath='/'
    for dir in out_path.split('/'):
        newpath=os.path.join(newpath,dir)
        # do not issue mkdir from very top of the tree
        if newpath.find('test_out') > 0:
            #getCommandOutput("eos mkdir"+newpath)
            command="/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir "+newpath
            p = subprocess.Popen(command,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (out, err) = p.communicate()
            #print out,err
            p.wait()

    # now check that the directory exists
    command2="/afs/cern.ch/project/eos/installation/cms/bin/eos.select ls "+out_path
    p = subprocess.Popen(command2,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = p.communicate()
    p.wait()
    if p.returncode !=0:
        print out

def split(sequence, size):
##########################    
# aux generator function to split lists
# based on http://sandrotosi.blogspot.com/2011/04/python-group-list-in-sub-lists-of-n.html
# about generators see also http://stackoverflow.com/questions/231767/the-python-yield-keyword-explained
##########################
    for i in xrange(0, len(sequence), size):
        yield sequence[i:i+size] 

#############
class Job:
#############

    def __init__(self, job_id, job_name, isDA, isMC, applyBOWS, applyEXTRACOND, extraconditions, runboundary, lumilist, maxevents, gt, allFromGT, alignmentDB, alignmentTAG, apeDB, apeTAG, bowDB, bowTAG, vertextype, tracktype, applyruncontrol, ptcut, CMSSW_dir ,the_dir):
###############################
        self.job_id=job_id    
        self.batch_job_id = None 
        self.job_name=job_name
        
        self.isDA              = isDA             
        self.isMC              = isMC             
        self.applyBOWS         = applyBOWS
        self.applyEXTRACOND    = applyEXTRACOND
        self.extraCondVect     = extraconditions
        self.runboundary       = runboundary         
        self.lumilist          = lumilist         
        self.maxevents         = maxevents
        self.gt                = gt
        self.allFromGT         = allFromGT
        self.alignmentDB       = alignmentDB      
        self.alignmentTAG      = alignmentTAG     
        self.apeDB             = apeDB            
        self.apeTAG            = apeTAG           
        self.bowDB             = bowDB            
        self.bowTAG            = bowTAG           
        self.vertextype        = vertextype       
        self.tracktype         = tracktype        
        self.applyruncontrol   = applyruncontrol  
        self.ptcut             = ptcut            

        self.the_dir=the_dir
        self.CMSSW_dir=CMSSW_dir

        self.output_full_name=self.getOutputBaseName()+"_"+str(self.job_id)

        self.cfg_dir=None
        self.outputCfgName=None
        
        # LSF variables        
        self.LSF_dir=None
        self.output_LSF_name=None

        self.lfn_list=list()      

        #self.OUTDIR = "/eos/cern.ch/user/m/musich/ZbbAnalysis/test01Sept" # TODO: write a setter method
        #self.OUTDIR = self.createEOSout()

    def __del__(self):
###############################
        del self.lfn_list

    def setEOSout(self,theEOSdir):    
###############################
        self.OUTDIR = theEOSdir
          
    def getOutputBaseName(self):
########################    
        return "PVValidation_"+self.job_name
        
    def createTheCfgFile(self,lfn):
###############################
        
        # write the cfg file 
        self.cfg_dir = os.path.join(self.the_dir,"cfg")
        if not os.path.exists(self.cfg_dir):
            os.makedirs(self.cfg_dir)

        self.outputCfgName=self.output_full_name+"_cfg.py"
        fout=open(os.path.join(self.cfg_dir,self.outputCfgName),'w+b')

        # decide which template according to data/mc
        if self.isMC:
            template_cfg_file = os.path.join(self.the_dir,"PVValidation_T_cfg.py")
        else:
            template_cfg_file = os.path.join(self.the_dir,"PVValidation_T_cfg.py")

        fin = open(template_cfg_file)

        for line in fin.readlines():
            if(to_bool(self.applyEXTRACOND)):
                if 'END OF EXTRA CONDITIONS' in line:
                    for element in self.extraCondVect :
                        if("Rcd" in element):
                            params = self.extraCondVect[element].split(',')
                            
                            fout.write(" \n")
                            fout.write("          process.conditionsIn"+element+"= CalibTracker.Configuration.Common.PoolDBESSource_cfi.poolDBESSource.clone( \n")
                            fout.write("               connect = cms.string('"+params[0]+"'), \n")
                            fout.write("               toGet = cms.VPSet(cms.PSet(record = cms.string('"+element+"'), \n")
                            fout.write("                                          tag = cms.string('"+params[1]+"'), \n")
                            if (len(params)>2):
                                fout.write("                                            label = cms.string('"+params[2]+"') \n")
                            fout.write("                                           ) \n")
                            fout.write("                                 ) \n")
                            fout.write("               ) \n")
                            fout.write("          process.prefer_conditionsIn"+element+" = cms.ESPrefer(\"PoolDBESSource\", \"conditionsIn"+element+"\") \n \n") 
                        
            if self.isMC:
                if line.find("ISDATEMPLATE")!=-1:
                    line=line.replace("ISDATEMPLATE",self.isDA)
                if line.find("ISMCTEMPLATE")!=-1:
                    line=line.replace("ISMCTEMPLATE",self.isMC)
                if line.find("APPLYBOWSTEMPLATE")!=-1:
                    line=line.replace("APPLYBOWSTEMPLATE",self.applyBOWS)
                if line.find("EXTRACONDTEMPLATE")!=-1:
                    line=line.replace("EXTRACONDTEMPLATE",self.applyEXTRACOND)
                if line.find("RUNBOUNDARYTEMPLATE")!=-1:
                    line=line.replace("RUNBOUNDARYTEMPLATE",self.runboundary)  
                if line.find("LUMILISTTEMPLATE")!=-1:
                    line=line.replace("LUMILISTTEMPLATE",self.lumilist)
                if line.find("MAXEVENTSTEMPLATE")!=-1:
                    line=line.replace("MAXEVENTSTEMPLATE",self.maxevents)
                if line.find("GLOBALTAGTEMPLATE")!=-1:
                    line=line.replace("GLOBALTAGTEMPLATE",self.gt)    
                if line.find("ALLFROMGTTEMPLATE")!=-1:
                    line=line.replace("ALLFROMGTTEMPLATE",self.allFromGT)    
                if line.find("ALIGNOBJTEMPLATE")!=-1:
                    line=line.replace("ALIGNOBJTEMPLATE",self.alignmentDB)
                if line.find("GEOMTAGTEMPLATE")!=-1:
                    line=line.replace("GEOMTAGTEMPLATE",self.alignmentTAG)
                if line.find("APEOBJTEMPLATE")!=-1:
                    line=line.replace("APEOBJTEMPLATE",self.apeDB)
                if line.find("ERRORTAGTEMPLATE")!=-1:
                    line=line.replace("ERRORTAGTEMPLATE",self.apeTAG)
                if line.find("BOWSOBJECTTEMPLATE")!=-1:
                    line=line.replace("BOWSOBJECTTEMPLATE",self.bowDB)  
                if line.find("BOWSTAGTEMPLATE")!=-1:
                    line=line.replace("BOWSTAGTEMPLATE",self.bowTAG)
                if line.find("VERTEXTYPETEMPLATE")!=-1:
                    line=line.replace("VERTEXTYPETEMPLATE",self.vertextype) 
                if line.find("TRACKTYPETEMPLATE")!=-1:
                    line=line.replace("TRACKTYPETEMPLATE",self.tracktype) 
                if line.find("PTCUTTEMPLATE")!=-1:
                    line=line.replace("PTCUTTEMPLATE",self.ptcut) 
                if line.find("RUNCONTROLTEMPLATE")!=-1:
                    line=line.replace("RUNCONTROLTEMPLATE",self.applyruncontrol) 
            else:                    
                if line.find("ISDATTEMPLATE")!=-1:
                    line=line.replace("ISDATEMPLATE",self.isDA)
                if line.find("ISMCTEMPLATE")!=-1:
                    line=line.replace("ISMCTEMPLATE",self.isMC)
                if line.find("APPLYBOWSTEMPLATE")!=-1:
                    line=line.replace("APPLYBOWSTEMPLATE",self.applyBOWS)
                if line.find("EXTRACONDTEMPLATE")!=-1:
                    line=line.replace("EXTRACONDTEMPLATE",self.applyEXTRACOND)
                if line.find("RUNBOUNDARYTEMPLATE")!=-1:
                    line=line.replace("RUNBOUNDARYTEMPLATE",self.runboundary)        
                if line.find("LUMILISTTEMPLATE")!=-1:
                    line=line.replace("LUMILISTTEMPLATE",self.lumilist)
                if line.find("MAXEVENTSTEMPLATE")!=-1:
                    line=line.replace("MAXEVENTSTEMPLATE",self.maxevents)
                if line.find("GLOBALTAGTEMPLATE")!=-1:
                    line=line.replace("GLOBALTAGTEMPLATE",self.gt) 
                if line.find("ALLFROMGTTEMPLATE")!=-1:
                    line=line.replace("ALLFROMGTTEMPLATE",self.allFromGT)    
                if line.find("ALIGNOBJTEMPLATE")!=-1:
                    line=line.replace("ALIGNOBJTEMPLATE",self.alignmentDB)
                if line.find("GEOMTAGTEMPLATE")!=-1:
                    line=line.replace("GEOMTAGTEMPLATE",self.alignmentTAG)
                if line.find("APEOBJTEMPLATE")!=-1:
                    line=line.replace("APEOBJTEMPLATE",self.apeDB)
                if line.find("ERRORTAGTEMPLATE")!=-1:
                    line=line.replace("ERRORTAGTEMPLATE",self.apeTAG)
                if line.find("BOWSOBJECTTEMPLATE")!=-1:
                    line=line.replace("BOWSOBJECTTEMPLATE",self.bowDB)  
                if line.find("BOWSTAGTEMPLATE")!=-1:
                    line=line.replace("BOWSTAGTEMPLATE",self.bowTAG)
                if line.find("VERTEXTYPETEMPLATE")!=-1:
                    line=line.replace("VERTEXTYPETEMPLATE",self.vertextype) 
                if line.find("TRACKTYPETEMPLATE")!=-1:
                    line=line.replace("TRACKTYPETEMPLATE",self.tracktype) 
                if line.find("PTCUTTEMPLATE")!=-1:
                    line=line.replace("PTCUTTEMPLATE",self.ptcut) 
                if line.find("RUNCONTROLTEMPLATE")!=-1:
                    line=line.replace("RUNCONTROLTEMPLATE",self.applyruncontrol) 
                        
            if line.find("FILESOURCETEMPLATE")!=-1:
                lfn_with_quotes = map(lambda x: "\'"+x+"\'",lfn)                   
                #print "["+",".join(lfn_with_quotes)+"]"
                line=line.replace("FILESOURCETEMPLATE","["+",".join(lfn_with_quotes)+"]") 
            if line.find("OUTFILETEMPLATE")!=-1:
                line=line.replace("OUTFILETEMPLATE",self.output_full_name+".root")     
            fout.write(line)    
      
        fout.close()                
                          
    def createTheLSFFile(self):
###############################

       # directory to store the LSF to be submitted
        self.LSF_dir = os.path.join(self.the_dir,"LSF")
        if not os.path.exists(self.LSF_dir):
            os.makedirs(self.LSF_dir)

        self.output_LSF_name=self.output_full_name+".lsf"
        fout=open(os.path.join(self.LSF_dir,self.output_LSF_name),'w')
    
        job_name = self.output_full_name

        log_dir = os.path.join(self.the_dir,"log")
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)

        fout.write("#!/bin/sh \n") 
        fout.write("#BSUB -L /bin/sh\n")       
        fout.write("#BSUB -J "+job_name+"\n")
        fout.write("#BSUB -o "+os.path.join(log_dir,job_name+".log")+"\n")
        fout.write("#BSUB -q cmscaf1nd \n")
        fout.write("JobName="+job_name+" \n")
        fout.write("OUT_DIR="+self.OUTDIR+" \n")
        fout.write("LXBATCH_DIR=`pwd` \n") 
        fout.write("cd "+os.path.join(self.CMSSW_dir,"src")+" \n")
        fout.write("eval `scram runtime -sh` \n")
        fout.write("cd $LXBATCH_DIR \n") 
        fout.write("cmsRun "+os.path.join(self.cfg_dir,self.outputCfgName)+" \n")
        fout.write("ls -lh . \n")
        fout.write("for RootOutputFile in $(ls *root ); do cmsStage -f ${RootOutputFile}  ${OUT_DIR}/${RootOutputFile} ; done \n")
        fout.write("for TxtOutputFile in $(ls *txt ); do cmsStage -f ${TxtOutputFile}  ${OUT_DIR}/${TxtOutputFile} ; done \n")

        fout.close()

    def getOutputFileName(self):
############################################
        return os.path.join(self.OUTDIR,self.output_full_name+".root")
        
    def submit(self):
###############################        
        print "submit job", self.job_id
        job_name = self.output_full_name
        submitcommand1 = "chmod u+x " + os.path.join(self.LSF_dir,self.output_LSF_name)
        child1  = os.system(submitcommand1)
        #submitcommand2 = "bsub < "+os.path.join(self.LSF_dir,self.output_LSF_name)
        #child2  = os.system(submitcommand2)
        self.batch_job_id = getCommandOutput("bsub < "+os.path.join(self.LSF_dir,self.output_LSF_name))

    def getBatchjobId(self):    
############################################
       return self.batch_job_id.split("<")[1].split(">")[0] 

##############################################
def main():
##############################################

    # CMSSW section
    input_CMSSW_BASE = os.environ.get('CMSSW_BASE')
    AnalysisStep_dir = os.path.join(input_CMSSW_BASE,"src/Alignment/OfflineValidation/test")
    lib_path = os.path.abspath(AnalysisStep_dir)
    sys.path.append(lib_path)

    ## N.B.: this is dediced here once and for all
    srcFiles        = []

    desc="""This is a description of %prog."""
    parser = OptionParser(description=desc,version='%prog version 0.1')
    parser.add_option('-s','--submit',    help='job submitted',    dest='submit',     action='store_true',  default=False)
    parser.add_option('-j','--jobname',   help='task name',        dest='taskname',   action='store',       default='myTask')
    parser.add_option('-D','--dataset',   help='selected dataset', dest='data',       action='store'      , default='')
    parser.add_option('-r','--doRunBased',help='selected dataset', dest='doRunBased', action='store_true' , default=False)
    parser.add_option('-i','--input',     help='set input configuration (overrides default)', dest='inputconfig',action='store',default=None)
    parser.add_option('-b','--begin',  help='starting point',        dest='start',  action='store'     ,default='1')
    parser.add_option('-e','--end',    help='ending point',          dest='end',    action='store'     ,default='999999')
   
    (opts, args) = parser.parse_args()

    now = datetime.datetime.now()
    t = now.strftime("test_%Y_%m_%d_%H_%M_%S_DATA_")
    t+=opts.taskname
    
    USER = os.environ.get('USER')
    eosdir=os.path.join("/store/caf/user",USER,"test_out",t)
    
    if opts.submit:
        mkdir_eos(eosdir)
    else:
        print "Not going to create EOS folder. -s option has not been chosen"

    #### Initialize all the variables

    jobName         = []
    isMC            = []
    isDA            = []
    doRunBased      = []
    maxevents       = []

    gt              = []
    allFromGT       = []
    applyEXTRACOND  = []
    extraCondVect   = []      
    alignmentDB     = []
    alignmentTAG    = []
    apeDB           = []
    apeTAG          = []
    applyBOWS       = []
    bowDB           = []
    bowTAG          = []
    conditions      = []
    
    vertextype      = []
    tracktype       = []

    applyruncontrol = []
    ptcut           = []
    runboundary     = []
    lumilist        = []
      
    ConfigFile = opts.inputconfig
    
    if ConfigFile is not None:

        print "********************************************************"
        print "* Parsing from input file:", ConfigFile," "
        
        #config = ConfigParser.ConfigParser()
        #config.read(ConfigFile)

        config = BetterConfigParser()
        config.read(ConfigFile)

        #print  config.sections()

        # please notice: since in principle one wants to run on several different samples simultaneously,
        # all these inputs are vectors

        doRunBased       = opts.doRunBased

        listOfValidations = config.getResultingSection("validations")
        
        for item in listOfValidations:
            if (bool(listOfValidations[item]) == True):
                
                jobName.append(ConfigSectionMap(config,"Conditions:"+item)['jobname'])
                isDA.append(ConfigSectionMap(config,"Job")['isda'])
                isMC.append(ConfigSectionMap(config,"Job")['ismc'])
                maxevents.append(ConfigSectionMap(config,"Job")['maxevents'])

                gt.append(ConfigSectionMap(config,"Conditions:"+item)['gt'])
                allFromGT.append(ConfigSectionMap(config,"Conditions:"+item)['allFromGT'])
                applyEXTRACOND.append(ConfigSectionMap(config,"Conditions:"+item)['applyextracond'])
                conditions.append(config.getResultingSection("ExtraConditions"))
                
                alignmentDB.append(ConfigSectionMap(config,"Conditions:"+item)['alignmentdb'])
                alignmentTAG.append(ConfigSectionMap(config,"Conditions:"+item)['alignmenttag'])
                apeDB.append(ConfigSectionMap(config,"Conditions:"+item)['apedb'])
                apeTAG.append(ConfigSectionMap(config,"Conditions:"+item)['apetag'])
                applyBOWS.append(ConfigSectionMap(config,"Conditions:"+item)['applybows'])
                bowDB.append(ConfigSectionMap(config,"Conditions:"+item)['bowdb'])
                bowTAG.append(ConfigSectionMap(config,"Conditions:"+item)['bowtag'])
                
                vertextype.append(ConfigSectionMap(config,"Type")['vertextype'])     
                tracktype.append(ConfigSectionMap(config,"Type")['tracktype'])
                
                applyruncontrol.append(ConfigSectionMap(config,"Selection")['applyruncontrol'])
                ptcut.append(ConfigSectionMap(config,"Selection")['ptcut'])
                runboundary.append(ConfigSectionMap(config,"Selection")['runboundary'])
                lumilist.append(ConfigSectionMap(config,"Selection")['lumilist'])
                
    else :

        print "********************************************************"
        print "* Parsing from command line                            *"
        print "********************************************************"
          
        jobName         = ['testing']
        isDA            = ['True']   
        isMC            = ['True']
        doRunBased      = opts.doRunBased
        maxevents       = ['10000']
        
        gt              = ['74X_dataRun2_Prompt_v4']       
        allFromGT       = ['False']
        applyEXTRACOND  = ['False']
        conditions      = [[('SiPixelTemplateDBObjectRcd','frontier://FrontierProd/CMS_CONDITIONS','SiPixelTemplateDBObject_38T_2015_v3_hltvalidation')]]
        alignmentDB     = ['frontier://FrontierProd/CMS_CONDITIONS']
        alignmentTAG    = ['TrackerAlignment_Prompt']  
        apeDB           = ['frontier://FrontierProd/CMS_CONDITIONS']  
        apeTAG          = ['TrackerAlignmentExtendedErr_2009_v2_express_IOVs']
        applyBOWS       = ['True']  
        bowDB           = ['frontier://FrontierProd/CMS_CONDITIONS']  
        bowTAG          = ['TrackerSurafceDeformations_v1_express']  
        
        vertextype      = ['offlinePrimaryVertices']  
        tracktype       = ['ALCARECOTkAlMinBias']  
        
        applyruncontrol = ['False']  
        ptcut           = ['3'] 
        runboundary     = ['1']  
        lumilist        = ['']  
 
    # print some of the configuration
    
    print "********************************************************"
    print "* Configuration info *"
    print "********************************************************"
    print "- submitted   : ",opts.submit
    print "- taskname    : ",opts.taskname
    print "- Jobname     : ",jobName           
    print "- use DA      : ",isDA            
    print "- is MC       : ",isMC            
    print "- is run-based: ",doRunBased
    print "- evts/job    : ",maxevents                    
    print "- GlobatTag   : ",gt      
    print "- allFromGT?  : ",allFromGT
    print "- extraCond?  : ",applyEXTRACOND
    print "- extraCond   : ",conditions                 
    print "- Align db    : ",alignmentDB     
    print "- Align tag   : ",alignmentTAG    
    print "- APE db      : ",apeDB           
    print "- APE tag     : ",apeTAG          
    print "- use bows?   : ",applyBOWS       
    print "- K&B db      : ",bowDB
    print "- K&B tag     : ",bowTAG                        
    print "- VertexColl  : ",vertextype      
    print "- TrackColl   : ",tracktype                       
    print "- RunControl? : ",applyruncontrol 
    print "- Pt>           ",ptcut           
    print "- run=          ",runboundary     
    print "- JSON        : ",lumilist  
    print "- Out Dir     : ",eosdir

    print "********************************************************"
    print "Will run on",len(jobName),"workflows"

    myRuns = []
    mylist = {}

    if(doRunBased):
        print ">>>> This is Data!"
        print ">>>> Doing run based selection"
        cmd = 'das_client.py --limit=0 --query \'run dataset='+opts.data+'\''
        p = Popen(cmd , shell=True, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        listOfRuns=out.split('\n')
        listOfRuns.pop()
        listOfRuns.sort()
        #myRuns = listOfRuns
        print "Will run on ",len(listOfRuns),"runs"
        print listOfRuns

        mytuple=[]

        for run in listOfRuns:
            if (run<=opts.start or run>opts.end):
                continue
            #print "preparing run",run

            if not isInJSON(run,lumilist[0]):
                continue
            
            # if(int(run)%100==0):
            mytuple.append((run,opts.data))

        #print mytuple

        pool = multiprocessing.Pool(processes=10)              # start 10 worker processes
        
        count = pool.map(getFilesForRun,mytuple)
        file_info = dict(zip(listOfRuns, count))

        print file_info

        for run in listOfRuns:
            if (run<=opts.start or run>opts.end):
                continue

            if not isInJSON(run,lumilist[0]):
                continue

            #if(int(run)%100==0):
            #    print "preparing run",run
            myRuns.append(run)
            #cmd2 = ' das_client.py --limit=0 --query \'file run='+run+' dataset='+opts.data+'\''
            #q = Popen(cmd2 , shell=True, stdout=PIPE, stderr=PIPE)
            #out2, err2 = q.communicate()
        
            #out2=getFilesForRun((run,opts.data))
            #print out2
            #pool.map(getFilesForRun,run,opts.data)


            #if run in file_info:
                #mylist[run] = file_info[run]
                #print run,mylist[run]
            #mylist[run] = out2.split('\n')
            #print mylist
            #mylist[run].pop()
            #print mylist
    
    od = collections.OrderedDict(sorted(file_info.items()))
    # print od
            
    # start loop on samples
    for iConf in range(len(jobName)):
        print "This is Task n.",iConf+1,"of",len(jobName)
        
        ##  print "==========>",conditions

        # for hadd script
        scripts_dir = os.path.join(AnalysisStep_dir,"scripts")
        if not os.path.exists(scripts_dir):
            os.makedirs(scripts_dir)
        hadd_script_file = os.path.join(scripts_dir,jobName[iConf]+"_"+opts.taskname+".sh")
        fout = open(hadd_script_file,'w')

        output_file_list1=list()      
        output_file_list2=list()
        output_file_list2.append("hadd ")
              
        inputFiles = []
   
        if (to_bool(isMC[iConf])):
            print ">>>> This is MC!"
            cmd = 'das_client.py --limit=0 --query \'file dataset='+opts.data+'\''
            s = Popen(cmd , shell=True, stdout=PIPE, stderr=PIPE)
            out,err = s.communicate()
            mylist2 = out.split('\n')
            mylist2.pop()
            #print mylist
           
            splitList = split(mylist2,10)
            for files in splitList:
                inputFiles.append(files)
                myRuns.append(str(1))

        else:
            #pass
            for element in od:
                #print mylist[element]
                inputFiles.append(od[element])
                #print element,od[element]
            #print mylist

        #print inputFiles

        ## declare here the list of jobs that should be waited for
        batchJobIds = []
        mergedFile = None

        for jobN,theSrcFiles in enumerate(inputFiles):
            print "JOB:",jobN,"run",myRuns[jobN],theSrcFiles
            thejobIndex=None
         
            ##### would be nice to store but it is yet another DAS query
            #eventsInRun=getNEvents(myRuns[jobN],opts.data)
            #print eventsInRun

            runInfo = {}
            runInfo["1.run"]             = myRuns[jobN]
            #runInfo["1.1runevents"]      = getNEvents(myRuns[jobN],opts.data) 
            runInfo["2.conf"]            = jobName[iConf]
            runInfo["3.gt"]              = gt[iConf]
            runInfo["4.allFromGT"]       = allFromGT[iConf] 
            runInfo["5.alignmentDB"]     = alignmentDB[iConf]
            runInfo["6.alignmentTag"]    = alignmentTAG[iConf]
            runInfo["7.apeDB"]           = apeDB[iConf]
            runInfo["8.apeTag"]          = apeTAG[iConf]
            runInfo["9.applyBows"]       = applyBOWS[iConf]
            runInfo["9.1bowDB"]          = bowDB[iConf]
            runInfo["9.2bowTag"]         = bowTAG[iConf]
            runInfo["9.3ptCut"]          = ptcut[iConf]
            runInfo["9.4lumilist"]       = lumilist[iConf]
            runInfo["9.5applyEXTRACOND"] = applyEXTRACOND[iConf]
            runInfo["9.6conditions"]     = conditions[iConf]
            runInfo["9.7nfiles"]         = len(theSrcFiles)
            runInfo["9.8srcFiles"]       = theSrcFiles

            updateDB(((iConf+1)*10)+(jobN+1),runInfo)
            #dumpDB()

            #if(to_bool(isMC[iConf]) and (not to_bool(doRunBased))):
            if(to_bool(isMC[iConf])):
                thejobIndex=jobN
            else:
                thejobIndex=myRuns[jobN]

            #print 'the configuration is:',iConf,' theJobIndex is:',thejobIndex
            #print applyBOWS[iConf],applyEXTRACOND[iConf],conditions[iConf]

            aJob = Job(thejobIndex,
                       jobName[iConf],isDA[iConf],isMC[iConf],
                       applyBOWS[iConf],applyEXTRACOND[iConf],conditions[iConf],
                       myRuns[jobN], lumilist[iConf], maxevents[iConf],
                       gt[iConf],allFromGT[iConf],
                       alignmentDB[iConf], alignmentTAG[iConf],
                       apeDB[iConf], apeTAG[iConf],
                       bowDB[iConf], bowTAG[iConf],
                       vertextype[iConf], tracktype[iConf],
                       applyruncontrol[iConf],
                       ptcut[iConf],input_CMSSW_BASE,AnalysisStep_dir)
            
            aJob.setEOSout(eosdir)
            aJob.createTheCfgFile(theSrcFiles)
            aJob.createTheLSFFile()

            output_file_list1.append("cmsStage "+aJob.getOutputFileName()+" /tmp/$USER/"+opts.taskname+" \n")
            if jobN == 0:
                mergedFile = "/tmp/$USER/"+opts.taskname+"/"+aJob.getOutputBaseName()+" "+opts.taskname+".root"
                output_file_list2.append("/tmp/$USER/"+opts.taskname+"/"+aJob.getOutputBaseName()+opts.taskname+".root ")
            output_file_list2.append("/tmp/$USER/"+opts.taskname+"/"+os.path.split(aJob.getOutputFileName())[1]+" ")    
   
            if opts.submit:
                aJob.submit()
                batchJobIds.append(aJob.getBatchjobId())
            del aJob

        if opts.submit:
            print "********************************************************"
            for theBatchJobId in batchJobIds:
                print "theBatchJobId is: ",theBatchJobId

        fout.write("#!/bin/bash \n")
        fout.write("MAIL=$USER@mail.cern.ch \n")
        fout.write("OUT_DIR="+eosdir+"\n")
        fout.write("FILE="+str(mergedFile)+"\n")
        fout.write("echo $HOST | mail -s \"Harvesting job started\" $USER@mail.cern.ch \n")
        fout.write("cd "+os.path.join(input_CMSSW_BASE,"src")+"\n")
        fout.write("eval `scram r -sh` \n")
        fout.write("mkdir -p /tmp/$USER/"+opts.taskname+" \n")
        fout.writelines(output_file_list1)
        fout.writelines(output_file_list2)
        fout.write("\n")
        fout.write("echo \"cmsStage -f $FILE $OUT_DIR\" \n")
        fout.write("cmsStage -f $FILE $OUT_DIR \n")
        fout.write("echo \"Harvesting for "+opts.taskname+" task is complete; please find output at $OUT_DIR \" | mail -s \"Harvesting for " +opts.taskname +" compled\" $MAIL \n")

        os.system("chmod u+x "+hadd_script_file)

        harvest_conditions = '"' + " && ".join(["ended(" + jobId + ")" for jobId in batchJobIds]) + '"'
        print harvest_conditions
        lastJobCommand = "bsub -o harvester"+opts.taskname+".tmp -q 1nh -w "+harvest_conditions+" "+hadd_script_file
        print lastJobCommand
        if opts.submit:
            lastJobOutput = getCommandOutput(lastJobCommand)
            print lastJobOutput

        fout.close()
        del output_file_list1
        
if __name__ == "__main__":        
    main()


   
