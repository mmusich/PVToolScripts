#!/usr/bin/env python
'''Script that plots all files in a root file
'''

__author__ = 'Marco Musich'
__copyright__ = 'Copyright 2012, CERN CMS'
__credits__ = ['Marco Musich']
__license__ = 'Unknown'
__maintainer__ = 'Marco Musich'
__email__ = 'musich@cern.ch'
__version__ = 1

import os,sys
import ROOT
import math
from xml.dom.minidom import parse
from optparse import OptionParser

#############
class Sample:
#############
    """ class to map the Sample elements in the xml file """

    def __init__(self, is_data, root_file, label, color, marker_style):
        self.is_data       = is_data
        self.the_root_file = root_file
        self.the_label     = label

        if color == 'ROOT.kRed':
            self.the_color = ROOT.kRed
        elif color == 'ROOT.kGreen':
            self.the_color = ROOT.kGreen
        elif color == 'ROOT.kAzure+10':
            self.the_color = ROOT.kAzure+10
        elif color == 'ROOT.kBlack':
            self.the_color = ROOT.kBlack
        elif color == 'ROOT.kMagenta':
            self.the_color = ROOT.kMagenta
        elif color == 'ROOT.kCyan':
            self.the_color = ROOT.kCyan
        elif color == 'ROOT.kBlue':
            self.the_color = ROOT.kBlue

        if marker_style == 'ROOT.kDot':
            self.the_marker_style = ROOT.kDot
        elif marker_style == 'ROOT.kOpenTriangleUp':
            self.the_marker_style = ROOT.kOpenTriangleUp
        elif marker_style == 'ROOT.kOpenTriangleDown':
            self.the_marker_style = ROOT.kOpenTriangleDown
        elif marker_style == 'ROOT.kOpenCircle':
            self.the_marker_style = ROOT.kOpenCircle
        elif marker_style == 'ROOT.kFullTriangleUp':
            self.the_marker_style = ROOT.kFullTriangleUp
        elif marker_style == 'ROOT.kFullTriangleDown':
            self.the_marker_style = ROOT.kFullTriangleDown
        elif marker_style == 'ROOT.kOpenSquare':
            self.the_marker_style = ROOT.kOpenSquare
        elif marker_style == 'ROOT.FullSquare':
            self.the_marker_style = ROOT.kFullSquare
        elif marker_style == 'ROOT.kFullCircle':
            self.the_marker_style = ROOT.kFullCircle

nplots = 0
######################
def set_global_var():
######################
    global USER
    global HOME
 
    USER = os.environ.get('USER')
    HOME = os.environ.get('HOME')
 
######################
def increment_count():       
######################
    global nplots
    nplots=nplots+1

######################
def setLooks(mylabels,mycolors,mytypes,mymarkers):
######################
    global thelabels
    global thecolors
    global thetypes
    global themarkers
    thelabels = []
    thecolors = []
    thetypes  = []
    themarkers = []
    for label in mylabels:
        thelabels.append(label)
    for color in mycolors:
        thecolors.append(color)
    for thetype in mytypes:
        thetypes.append(thetype)
    for themarker in mymarkers:
        themarkers.append(themarker)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Graphics beatufication:
# - general style
# - histogram style
# - legends and statbox
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

###############
def setStyle():
###############
    ROOT.gStyle.SetTitleX(0.55)
    ROOT.gStyle.SetTitleAlign(23)
    ROOT.TH1.StatOverflows(ROOT.kTRUE)
    ROOT.gStyle.SetOptTitle(1)
    ROOT.gStyle.SetOptStat("emr")

    ROOT.gStyle.SetPadTopMargin(0.07)
    ROOT.gStyle.SetPadBottomMargin(0.14)
    ROOT.gStyle.SetPadLeftMargin(0.18)
    ROOT.gStyle.SetPadRightMargin(0.05)

    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetOptTitle(0)
    #ROOT.gStyle.SetTitleFillColor(10)
    #ROOT.gStyle.SetTitleFont(42)
    #ROOT.gStyle.SetTitleTextColor(ROOT.kBlue)
    #ROOT.gStyle.SetTitleFontSize(0.06)
    #ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetStatColor(ROOT.kWhite)
    ROOT.gStyle.SetStatFont(42)
    ROOT.gStyle.SetStatFontSize(0.05)
    ROOT.gStyle.SetStatTextColor(1)
    ROOT.gStyle.SetStatFormat("6.4g")
    ROOT.gStyle.SetStatBorderSize(1)
    ROOT.gStyle.SetPadTickX(1) 
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetOptFit(1)
    ROOT.gStyle.SetNdivisions(510)

############################
def getExtrema(h1array):
############################
    the_max = 0.
    the_min =999.

    for h1 in h1array:
        this_max = h1.GetMaximum()
        this_min = h1.GetMinimum();
        if this_max>the_max:
            the_max = this_max
        if this_min<the_min:
            the_min = this_min
        
    # print "Minimum: ", the_min ", Maximum: ", the_max
    return the_min, the_max

######################################################
def makeNicePlotStyle(hist, color, marker, the_extrema, doFill=False):
######################################################

    ROOT.TH1.StatOverflows(ROOT.kTRUE)

    #hist.SetStats(ROOT.kFALSE)  
    hist.SetStats(ROOT.kTRUE)
    hist.GetXaxis().CenterTitle(ROOT.kTRUE)
    hist.GetYaxis().CenterTitle(ROOT.kTRUE)
    hist.GetXaxis().SetTitleFont(42) 
    hist.GetYaxis().SetTitleFont(42)  
    hist.GetXaxis().SetTitleSize(0.062)
    hist.GetYaxis().SetTitleSize(0.062)
    hist.GetXaxis().SetTitleOffset(0.95)
    hist.GetYaxis().SetTitleOffset(1.4)
    hist.GetXaxis().SetLabelFont(42)
    hist.GetYaxis().SetLabelFont(42)
    hist.GetXaxis().SetLabelOffset(0.009)
    hist.GetYaxis().SetLabelOffset(0.009)

    hist.GetYaxis().SetLabelSize(.055)
    hist.GetXaxis().SetLabelSize(.055)
    hist.SetMarkerSize(1.2)
    hist.SetMarkerStyle(marker)
    hist.SetLineStyle(color)
    if(doFill):
        hist.SetFillColor(color)
        hist.SetLineColor(ROOT.kBlack)
    else:
        hist.SetLineColor(color)
    #hist.SetFillColorAlpha(color, 0.35)
    #hist.SetFillStyle(3004)
    hist.GetXaxis().SetNdivisions(505)
    hist.SetLineWidth(2)
    hist.SetMarkerColor(color)
    if(the_extrema[0]>0.):
        hist.GetYaxis().SetRangeUser(0.,the_extrema[1]*1.5)
    else:
        hist.GetYaxis().SetRangeUser(the_extrema[0]*1.5,the_extrema[1]*1.5)

    fit =  hist.GetListOfFunctions().FindObject("tmp")
    if(fit):  
        fit.Delete("")
        
######################################################
def makeNiceTLegend(hists,labels):
######################################################
    #lego = ROOT.TLegend(0.55,0.55,0.93,0.75)
    lego = ROOT.TLegend(0.20,0.75,0.35,0.91)
    lego.SetFillColor(10)
    lego.SetTextSize(0.04)
    lego.SetTextFont(42)
    lego.SetFillColor(10)
    lego.SetLineColor(10)
    lego.SetShadowColor(10)
           
    for i in xrange(0,len(hists)): 
        if("MC" not in labels[i]):
            lego.AddEntry(hists[i],labels[i],"PL")
        else:
            lego.AddEntry(hists[i],labels[i],"F") 

    return lego

######################################################
def makeNicePaveTest(isRatio=False):
######################################################
    if(isRatio):
        pt = ROOT.TPaveText(0.179,0.955,0.260,0.985,"NDC")
    else:
        pt = ROOT.TPaveText(0.2,0.955,0.260,0.985,"NDC")
    pt.SetFillColor(10)
    pt.SetTextColor(1)
    pt.SetTextFont(61)
    #pt.SetTextAlign(12)
    text1 = pt.AddText("CMS") #"Preliminary 2015 - 0T collision data");
    text1.SetTextSize(0.05)
 
    extraOverCmsTextSize  = 0.76

    if(isRatio):
        pt2 = ROOT.TPaveText(0.280,0.935,0.40,0.985,"NDC")
    else:
        pt2 = ROOT.TPaveText(0.285,0.935,0.453,0.985,"NDC")
    pt2.SetFillColor(10)
    pt2.SetTextColor(1)
    pt2.SetTextFont(52)
    #pt2.SetTextAlign(11)
    text2 = pt2.AddText("Preliminary")
    text2.SetTextSize(0.05*extraOverCmsTextSize)
    
    if(isRatio):
        pt3 = ROOT.TPaveText(0.642,0.95,0.98,0.99,"NDC")
    else:
        pt3 = ROOT.TPaveText(0.642,0.95,0.94,0.98,"NDC")
    pt3.SetFillColor(10)
    pt3.SetTextColor(1)
    pt3.SetTextFont(42)
    # pt2.SetTextAlign(11)
    text3 = pt3.AddText("0T collision data 2015")
    text3.SetTextSize(0.05*extraOverCmsTextSize)

    return pt,pt2,pt3
   
######################################################
def makeNiceStats(stats,index):
######################################################
    span=0.5

    stats.SetLineColor(thecolors[index])
    stats.SetTextColor(thecolors[index])
 
    stats.SetX1NDC(0.45+index*0.14)
    stats.SetX2NDC(0.58+index*0.14)
     
    stats.SetY1NDC(0.75)
    stats.SetY2NDC(0.90)
    stats.SetOptStat(1110)
   
#############################################
def loopIntoIt(files,key,outputFolder):
#############################################
    if (key.IsA().InheritsFrom("TDirectory")):
        #print key.GetName(), "is a directory in ",key.GetPath()
        here = key.GetPath().split(":",1)[1]
        #print here
        ROOT.gSystem.MakeDirectory(outputFolder+here)
        files[0].cd(str(key.GetPath()))
        the_dir = ROOT.gDirectory
        dirList = ROOT.gDirectory.GetListOfKeys()
        for k in dirList:
            obj_array = []
            obj1 = k.ReadObj()
            if (obj1.IsA().InheritsFrom("TH1")):#and "h2_" not in obj1.GetName()):
                if (obj1.IsA().InheritsFrom("TH2")): 
                    continue

                increment_count()
                print "processed ",nplots," histograms"
                #if (nplots>100):
                #    return

                #print obj1.GetName()," is an histogram in ",key.GetPath()

                listOfAcceptedHistograms = ["Residuals"]

                searchpath = key.GetPath().split(":/",1)[1]
                if (listOfAcceptedHistograms[0] in key.GetPath()):  #obj1.GetName()):
                    continue

                for the_file in files:    
                    the_file.cd(searchpath)
                    kf = ROOT.gDirectory.GetListOfKeys().FindObject(obj1.GetName())
                    objf = kf.ReadObj()
                    obj_array.append(objf)
          
                arr = []    
                for h1 in obj_array:     
                    if (h1.GetSumOfWeights()!=0 and ("p_" not in obj1.GetName())):
                        h1.Scale(100/h1.GetSumOfWeights())
                        ##h1.Scale(1/1000.)
                        arr.append(h1)
  
                the_extrema = getExtrema(arr)
                for h1 in obj_array:
                    if thetypes[obj_array.index(h1)]:
                        doFill = False
                    else:
                    #if (obj_array.index(h1)>0):
                        doFill = True
                        #print "themarker is :",themarkers[obj_array.index(h1)]
                    makeNicePlotStyle(h1,thecolors[obj_array.index(h1)],themarkers[obj_array.index(h1)],the_extrema,doFill)
        
                compound = str(obj1.GetName())+"_"+str(key.GetName())

                #####################################
                ##
                ## Standard plot canvas
                ##
                #####################################

                c1 = ROOT.TCanvas("c1_"+compound,"c1_"+compound,800,800)
                c1.cd()
         
                if ("p_" not in obj1.GetName()):
                    for h1 in obj_array:
                        if thetypes[obj_array.index(h1)]:
                            if obj_array.index(h1)==0:
                                h1.Draw("PE1")
                            else:
                                h1.Draw("PE1sames")
                        else:
                            if obj_array.index(h1)==0:
                                h1.Draw("HIST")
                            else: 
                                h1.Draw("HISTsames")          
                else:
                    obj_array[0].Draw("PE1")
                    for h1 in obj_array:  
                        h1.Draw("PE1sames")

                c1.Draw()

                clones = []
                for h1 in obj_array:  
                    clone = h1.Clone()
                    clones.append(clone)

                    #stats = h1.FindObject("stats")
                    #makeNiceStats(stats,obj_array.index(h1))
                    #stats.Draw("same")

                myleg = makeNiceTLegend(obj_array,thelabels)
                myleg.Draw("same")
                (p1,p2,p3) = makeNicePaveTest()
                p1.Draw("same")
                p2.Draw("same")
                p3.Draw("same")

                c1.SaveAs(outputFolder+here+"/c1_"+compound+".png")
                c1.SaveAs(outputFolder+here+"/c1_"+compound+".pdf")
                #c1.SaveAs(outputFolder+here+"/c1_"+compound+".root") 
                 
                #####################################
                ##
                ## ratio plot canvas
                ##
                #####################################
                
                c1ratio = ROOT.TCanvas("c1_ratio_%s" % c1.GetName() ,"c1_ratio_%s" % c1.GetName(),800,800)
                # upper pad
                pad = ROOT.TPad("mainPad_%s" % c1.GetName(),"mainPad_%s" % c1.GetName(),0.,0.17,1.,1.)
                pad.SetFillStyle(4000)
                ## lower pad
                pad_ratio = ROOT.TPad('ratioPad_%s' % c1.GetName(), 'ratioPad_%s' % c1.GetName() ,0.,0.,1.,0.29);
                pad_ratio.SetFillStyle(4000)
                pad_ratio.SetBottomMargin(0.45)
                pad.Draw()
                pad_ratio.Draw()
    
                pad.cd()
    
                ratio_hists = []
                denom_hists = []
                for clone in clones:  
                    ratio = clone.Clone(clone.GetName()+"_%s" % clones.index(clone) )
                    if("p_" in obj1.GetName() and obj1.IsA().InheritsFrom("TProfile")):
                        ratio = ratio.ProjectionX()
                        
                        #### old logic
                        # if (clones.index(clone)!=0):
                        #     ratio_hists.append(ratio)
                        #     if ("p_" not in obj1.GetName()):
                        #         clone.Draw("HISTsames")
                        #     else:
                        #         clone.Draw("HISTsames")
                        # else:
                        #     denom_hists.append(ratio)
                        #     if ("p_" not in obj1.GetName()):
                        #         clone.Draw("PE1")
                        #     else:
                        #         clone.Draw("PE1")

                    if thetypes[clones.index(clone)]:
                        # if this is data
                       if clones.index(clone)==0:
                           denom_hists.append(ratio)
                           clone.Draw("PE1")
                       else: 
                           ratio_hists.append(ratio)
                           clone.Draw("PE1same")
                    else:                             
                        # if this is mc
                        if clones.index(clone)==0:
                            denom_hists.append(ratio)
                            clone.Draw("HIST")
                        else: 
                            ratio_hists.append(ratio)
                            clone.Draw("HISTsame")

                    #stats = clone.FindObject("stats")
                    #makeNiceStats(stats,clones.index(clone))
                    #stats.Draw("same")

                    clone.GetXaxis().SetLabelSize(0.)
                  
                myleg.Draw("same")
                (p1r,p2r,p3r) = makeNicePaveTest(True)
                p1r.Draw("same")
                p2r.Draw("same")
                p3r.Draw("same")

                pad_ratio.cd()
                
                nbins = obj1.GetNbinsX()
                flatline = ROOT.TH1F("flatline_%s" % c1.GetName(),"flatline_%s" % c1.GetName(),
                                     nbins,obj1.GetXaxis().GetBinLowEdge(1),obj1.GetXaxis().GetBinLowEdge(nbins+1))

                flatline.SetMarkerSize(0)
                makeNicePlotStyle(flatline,ROOT.kBlue,20,(0.,2.))  
                flatline.SetLineColor(ROOT.kBlue)
                flatline.SetLineStyle(9)
                flatline.SetFillStyle(0)
                flatline.SetLineWidth(2)
                flatline.SetTitle("")
              
                for i in xrange(1,nbins+1):
                    flatline.SetBinContent(i,1)
                    flatline.SetBinError(i,0)
                    
                for ratio in ratio_hists:
                    ratio.Divide(denom_hists[0])
                    ratio.SetTitle("")
                    makeNicePlotStyle(ratio,thecolors[ratio_hists.index(ratio)+1],themarkers[ratio_hists.index(ratio)+1],(0.,2.)) 
                    ratio.SetLineColor(thecolors[ratio_hists.index(ratio)+1])
                    ratio.GetXaxis().SetTitleSize(0.19)
                    ratio.GetYaxis().SetTitleSize(0.19)
                    ratio.GetXaxis().SetTitleOffset(1.0)
                    ratio.GetYaxis().SetTitleOffset(0.45)
                    ratio.GetXaxis().SetLabelFont(42)
                    ratio.GetYaxis().SetLabelFont(42)
                    ratio.GetXaxis().SetLabelOffset(0.009)
                    ratio.GetYaxis().SetLabelOffset(0.009)
                    ratio.GetYaxis().SetLabelSize(0.16)
                    ratio.GetXaxis().SetLabelSize(0.16)
                    ratio.GetYaxis().SetNdivisions(505)
                    ratio.GetYaxis().SetRangeUser(0.,1.89)
                    ratio.GetYaxis().SetTitle("ratio")
                    if(ratio_hists.index(ratio)==0):
                        ratio.Draw("P")
                    else:
                        ratio.Draw("Psames")
                    ratio.Draw("e1sames")
                    ratio.SetStats(ROOT.kFALSE)
                    
                flatline.Draw("e1same")
                flatline.SetStats(ROOT.kFALSE)

                # #stack_ratio.GetYaxis().SetNdivisions(507) # Avoids crowded labels
                c1ratio.Draw()
                c1ratio.Update()

                c1ratio.SaveAs(outputFolder+here+"/c1_ratio_"+compound+".png")
                c1ratio.SaveAs(outputFolder+here+"/c1_ratio_"+compound+".pdf")
                #c1ratio.SaveAs(outputFolder+here+"/c1_ratio"+compound+".root") 

                #####################################
                ##
                ## log scale plot canvas
                ##
                #####################################
     
                c1log = ROOT.TCanvas("c1log_"+compound,"c1log_"+compound,800,800)
                c1log.cd().SetLogy()
                
                #print "min1:",obj1.GetMinimum()," min2:",obj2.GetMinimum()
                for h1 in obj_array:  
                    #print "min:",the_extrema[1],"max:",the_extrema[0]
                    dr = the_extrema[1] - the_extrema[0]

                    #log(max+ov) = 5/4*log(max)  ---> max+ov = 10^{5/4*log(max)}

                    if( (h1.GetMinimum()>0.)):
                        h1.GetYaxis().SetRangeUser(the_extrema[0],the_extrema[1])#+ROOT.TMath.Power(10,ldr))
                        #h1.GetYaxis().UnZoom()
                    else :
                        h1.GetYaxis().SetRangeUser(0.001,ROOT.TMath.Power(10,ROOT.TMath.Log10(10*the_extrema[1])))
                        #h1.GetYaxis().UnZoom()
                   
                    if ("p_" not in obj1.GetName()):
                        if thetypes[obj_array.index(h1)]:
                            if obj_array.index(h1)==0:
                                h1.Draw("PE1")
                            else:
                                h1.Draw("PE1sames")
                        else:
                            if obj_array.index(h1)==0:
                                h1.Draw("HIST")
                            else: 
                                h1.Draw("HISTsames")          
                    else:
                        if obj_array.index(h1)==0:
                            h1.Draw("PE1")
                        else:
                            h1.Draw("PE1sames")  
                      
                    #stats = h1.FindObject("stats")
                    #makeNiceStats(stats,obj_array.index(h1))
                    #stats.Draw("same")

                c1log.Draw()
                myleg.Draw("same")
                p1.Draw("same")
                p2.Draw("same")
                p3.Draw("same")

                c1log.SaveAs(outputFolder+here+"/c1_log_"+compound+".png")
                c1log.SaveAs(outputFolder+here+"/c1_log_"+compound+".pdf")
                #c1log.SaveAs(outputFolder+here+"/c1log_"+compound+".root") 

            elif (obj1.IsA().InheritsFrom("TDirectory")):
                print "here: ",ROOT.gDirectory.GetName()
                loopIntoIt(files,obj1,outputFolder)
                
    elif (key.IsA().InheritsFrom("TH1")):
        print key.GetName()," is an histogram"
        increment_count()
        #print "processed ",nplots," histograms"
        #if (nplots>100):
        #    return
        
        #print obj1.GetName()," is an histogram in ",key.GetPath()
        searchpath = key.GetPath().split(":/",1)[1]
        
        files[0].cd(searchpath)
        k2 = ROOT.gDirectory.GetListOfKeys().FindObject(obj1.GetName())
        obj2 = k2.ReadObj()
        
        arr = []
        arr.append(obj1)
        arr.append(obj2)
        
        the_extrema = getExtrema(arr)
        makeNicePlotStyle(obj1,ROOT.kBlue,20,the_extrema)
        makeNicePlotStyle(obj2,ROOT.kRed,20,the_extrema)
    
        compound = str(obj1.GetName())+"_"+str(key.GetName())
        c1 = ROOT.TCanvas("c1_"+compound,"c1_"+compound,800,600)
        c1.cd()
        if ("p_" not in obj1.GetName()):
            obj1.Draw("HIST")
        else:
            obj1.Draw("CP")
            obj2.Draw("CPsames")
            
        c1.Draw()
        ## syntax obj1,obj2,label1,label2
        myleg = makeNiceTLegend(obj1,obj2,obj3,"53X","70X_PESS","70X_OPT")
        myleg.Draw("same")
        
        #st1 = obj1.FindObject("stats")
        #st2 = obj2.FindObject("stats")
        #makeNiceStats(st1,ROOT.kBlue,st2,ROOT.kRed,st3,ROOT.kBlack)
        #st1.Draw("same")
        #st2.Draw("same")
        
        c1.SaveAs(outputFolder+here+"/c1_"+compound+".pdf")
        #c1.SaveAs("c1_"+compound+".root") 
        
#############################################
def main():
#############################################
    
    parser = OptionParser()
    parser.add_option("-f", "--file",  
                      action="store", type="string", dest="input_xml_file",
                      help="input XML file")
    (options, args) = parser.parse_args()

    dom = parse(options.input_xml_file)

    ### parse & build "Sample(s)"
    def handleSamples(samples, sample_list):
        for sample in samples:
            sample_list.append(handleSample(sample))


    def handleSample(sample):
        is_data = sample.getAttribute('Type') == 'data'
        s = Sample(is_data, \
               handleInputRootFile(sample.getElementsByTagName("InputRootFile")[0]), \
               handleLabel(sample.getElementsByTagName("Label")[0]), \
               handleColor(sample.getElementsByTagName("Color")[0]), \
               handleMarkerStyle(sample.getElementsByTagName("MarkerStyle")[0]) \
           )
        return s

    def handleInputRootFile(input_root_file):
        return input_root_file.firstChild.nodeValue

    def handleLabel(label):
        return label.firstChild.nodeValue

    def handleColor(color):
        return color.firstChild.nodeValue

    def handleMarkerStyle(marker_style):
        return marker_style.firstChild.nodeValue

    set_global_var()
    print nplots

    setStyle()

    # do not pop-up canvases as they are drawn
    ROOT.gROOT.SetBatch(ROOT.kTRUE) 

    outputFolder =  "just_testing2"
    ROOT.gSystem.MakeDirectory(outputFolder)
    ROOT.gStyle.SetOptStat(0)

    Samples  = []
    files    = []
    mylabels = []
    mycolors = []
    mytype   = []
    mymarker = []

    handleSamples(dom.getElementsByTagName('Sample'), Samples)
    for aSample in Samples:
        files.append(ROOT.TFile(aSample.the_root_file))
        mylabels.append(aSample.the_label)
        mycolors.append(aSample.the_color)
        mytype.append(aSample.is_data)
        mymarker.append(aSample.the_marker_style)

    setLooks(mylabels,mycolors,mytype,mymarker)

    files[0].cd()
    dirList = ROOT.gDirectory.GetListOfKeys()
    for element in dirList:
        print "evaluating:",element.GetName()
        obj = element.ReadObj()
        loopIntoIt(files,obj,outputFolder)
            
    currentDir = ROOT.gSystem.pwd()
    outDirPath = os.path.join(currentDir,outputFolder)

    print "Done. See images in:",outDirPath
  
if __name__ == "__main__":
    main()
