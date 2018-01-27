#!/usr/bin/env python
from ROOT import TFile, TTree, TH1F, TH1D, TH1, TCanvas, TChain,TGraphAsymmErrors, TMath, TH2D, TLorentzVector, TF1, AddressOf
import ROOT as ROOT
import os
import sys, optparse
from array import array
import math
import AllQuantList

ROOT.gROOT.SetBatch(True)
from bbMETQuantities import *
#from PileUpWeights import PUWeight
pileup2016file = TFile('pileUPinfo2016.root')
pileup2016histo=pileup2016file.Get('hpileUPhist')
eleReweightFile = TFile('eleTrig.root')
eleEtaPt = eleReweightFile.Get('hEffEtaPt')

ROOT.gROOT.ProcessLine('.L BTagCalibrationStandalone.cpp+') 

#ROOT.gROOT.ProcessLine('.L TheaCorrection.cpp+') 

######################################
## set up running mode of the code.
######################################

#ROOT.gROOT.ProcessLine('.L PileUpWeights.h')

#print "puweight = ",PUWEIGHT(10)
usage = "usage: %prog [options] arg1 arg2"
parser = optparse.OptionParser(usage)

## data will be true if -d is passed and will be false if -m is passed
parser.add_option("-i", "--inputfile",  dest="inputfile")
parser.add_option("-o", "--outputfile", dest="outputfile")
parser.add_option("-D", "--outputdir", dest="outputdir")

parser.add_option("-a", "--analyze", action="store_true",  dest="analyze")

parser.add_option("-e", "--efficiency", action="store_true",  dest="efficiency")
parser.add_option("-F", "--farmout", action="store_true",  dest="farmout")
parser.add_option("-t", "--table", action="store_true",  dest="table")
parser.add_option("-P", "--OtherPlots", action="store_true",  dest="OtherPlots")

parser.add_option("--csv", action="store_true",  dest="CSV")
parser.add_option("--deepcsv", action="store_true",  dest="DeepCSV")
########################################################################################################################
########################## cut values########################################################################
########################################################################################################################

parser.add_option("--dbt", action="store_true",  dest="dbt")
parser.add_option( "--dbtcut", type=float,  dest="dbtcut")

parser.add_option("--theac", action="store_true",  dest="theac")

(options, args) = parser.parse_args()

if options.farmout==None:
    isfarmout = False
else:
    isfarmout = options.farmout

if options.CSV==None:
    options.CSV = False

if options.DeepCSV==None:
    options.DeepCSV = False
    
if options.CSV: print "Using CSVv2 as b-tag discriminator."    
if options.DeepCSV: print "Using DeepCSV as b-tag discriminator."

if not options.CSV and not options.DeepCSV:
    print "Please run using --csv or --deepcsv. Exiting."
    sys.exit()

#print 'options = ',[options.inputfile]
inputfilename = options.inputfile
outputdir = options.outputdir

#print inputfilename
pathlist = inputfilename.split("/")
sizeoflist = len(pathlist)
#print ('sizeoflist = ',sizeoflist)
rootfile='tmphist'
rootfile = pathlist[sizeoflist-1]
textfile = rootfile+".txt"

#outputdir='bbMETSamples/'
if outputdir!='.': os.system('mkdir -p '+outputdir)

if options.outputfile is None or options.outputfile==rootfile:
   outputfilename = "/Output_"+rootfile
else:
   outputfilename = "/"+options.outputfile

#if isfarmout:
outfilename = outputdir + outputfilename
#else:
#    outfilename = options.outputfile    

print "Input:",options.inputfile, "; Output:", outfilename

skimmedTree = TChain("outTree")


#bbMET_tree = TTree( 'bbMET_tree', 'outputTree' )
#print isfarmout



def WhichSample(filename):
    samplename = 'all'
    if filename.find('WJets')>-1:
        samplename = 'WJETS'
    elif filename.find('ZJets')>-1 or filename.find('DYJets')>-1:
        samplename = 'ZJETS'
    elif filename.find('TT')>-1:
        samplename  = 'TT'
    else:
        samplename = 'all'
#    print samplename
    return samplename

def IsoMu20isUnPrescaled(filename):    
    if filename.find('2016D')>-1 or filename.find('2016E')>-1 or filename.find('2016F')>-1 or filename.find('2016G')>-1 or filename.find('2016H')>-1:
        return False
    else:
        return True

def TheaCorrection(puppipt=200.0,  puppieta=0.0):
    puppisd_corrGEN      = TF1("puppisd_corrGEN","[0]+[1]*pow(x*[2],-[3])");
    puppisd_corrGEN.SetParameters(
        1.00626,
        -1.06161,
        0.07999,
        1.20454
        )
    puppisd_corrRECO_cen =  TF1("puppisd_corrRECO_cen","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
    puppisd_corrRECO_cen.SetParameters(
        1.05807,
        -5.91971e-05,
        2.296e-07,
        -1.98795e-10,
        6.67382e-14,
        -7.80604e-18
        )
    
    puppisd_corrRECO_for = TF1("puppisd_corrRECO_for","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
    puppisd_corrRECO_for.SetParameters(
        1.26638,
        -0.000658496,
        9.73779e-07,
        -5.93843e-10,
        1.61619e-13,
        -1.6272e-17)
    
    genCorr  = 1.
    recoCorr = 1.
    totalWeight = 1.
    
    genCorr =  puppisd_corrGEN.Eval( puppipt )
    if ( abs(puppieta)  <= 1.3 ) :
        recoCorr = puppisd_corrRECO_cen.Eval( puppipt )
    elif( abs(puppieta) > 1.3 ) :
        recoCorr = puppisd_corrRECO_for.Eval( puppipt )
        
    totalWeight = genCorr * recoCorr
    return totalWeight
    



h_t = TH1F('h_t','h_t',2,0,2)
h_t_weight = TH1F('h_t_weight','h_t_weight',2,0,2)

samplename = 'all'
if isfarmout:
    infile = open(inputfilename)
    for ifile in infile: 
        skimmedTree.Add(ifile.rstrip())
#        samplename = WhichSample(ifile.rstrip())
        ## for histograms
        f_tmp = TFile.Open(ifile.rstrip(),'READ')
        h_tmp = f_tmp.Get('h_total')
        h_tmp_weight = f_tmp.Get('h_total_mcweight')
        h_t.Add(h_tmp)
        h_t_weight.Add(h_tmp_weight)

if not isfarmout:
    skimmedTree.Add(inputfilename)
#    samplename = WhichSample(inputfilename)
    ## for histograms
    f_tmp = TFile.Open(inputfilename,'READ')
    h_tmp = f_tmp.Get('h_total')
    h_tmp_weight = f_tmp.Get('h_total_mcweight')    
    h_t.Add(h_tmp)
    h_t_weight.Add(h_tmp_weight)

debug = False 

try:
    samplepath = str(f_tmp.Get('samplepath').GetTitle())
    print "Original source file: " + samplepath
except:
#    samplepath=inputfilename
    samplepath='TT'
    print "WARNING: Looks like the input was skimmed with an older version of SkimTree. Using " + samplepath + " as sample path. Gen pT Reweighting may NOT work."
    
samplename = WhichSample(samplepath)
print "Dataset classified as: " + samplename
UnPrescaledIsoMu20 = IsoMu20isUnPrescaled(samplepath)
#print UnPrescaledIsoMu20
#print samplename
#print

def AnalyzeDataSet():
    ## Input rootfile name
    
    #rootfilename = inputfilename
    #print (rootfilename,inputfilename)
    #f = TFile(rootfilename,'READ')
    #skimmedTree = f.Get('tree/treeMaker')
    NEntries = skimmedTree.GetEntries()
    print 'NEntries = '+str(NEntries)
    npass = 0
    
    triglist=['HLT_PFMET170_','HLT_PFMET170_NoiseCleaned','HLT_PFMET170_JetIdCleaned_v','HLT_PFMET170_HBHECleaned_v','HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v','HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v','HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v','HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v','HLT_PFMET110_PFMHT110_','HLT_IsoMu24_v','HLT_IsoTkMu24_v','HLT_Ele27_WPTight_Gsf','HLT_IsoMu20','HLT_Ele27_WPLoose_Gsf']    
        
    #print [rootfilename, NEntries]
    cutStatus={'preselection':NEntries}
    cutStatusSR1={'preselection':NEntries}
    cutStatusSR2={'preselection':NEntries}
    


    cutStatus['pfmet'] =  0 
    cutStatus['njet+nBjet'] = 0
    cutStatus['lep'] = 0
    cutStatus['jet1'] = 0
    cutStatus['jet2/3'] = 0

    cutStatusSR1['njet+nBjet'] = 0
    cutStatusSR1['lep'] = 0
    cutStatusSR1['jet1'] = 0
    cutStatusSR1['jet2'] = 0
    
    cutStatusSR2['njet+nBjet'] = 0
    cutStatusSR2['lep'] = 0
    cutStatusSR2['jet1'] = 0
    cutStatusSR2['jet2'] = 0
    cutStatusSR2['jet3'] = 0
    
#    CRCutFlow['njet+nBjet']=0
#    CRCutFlow['jetcond']=0
#    CRCutFlow['nlepcond']=0
#    CRCutFlow['zJet']=0
#    CRCutFlow['zLep']=0
#    CRCutFlow['zMass']=0
#    CRCutFlow['zrecoil']=0
#    CRCutFlow['ZdPhi']=0
    
    CRcutnames=['trig','nlep','lepconds','recoil','mass','dPhicond','nbjets','jetconds']
    regionnames=['2e1b','2mu1b','2e2b','2mu2b','1e1b','1mu1b','1e2b','1mu2b','1mu1e1b','1mu1e2b']
    for CRreg in regionnames:
        exec("CR"+CRreg+"CutFlow={'preselection':NEntries}")
        for cutname in CRcutnames:
            exec("CR"+CRreg+"CutFlow['"+cutname+"']=0")
    
    
    CRs=['ZCRSR1','ZCRSR2','WCRSR1','WCRSR2','TopCRSR1','TopCRSR2']
    
    CRStatus={'total':NEntries}
    for CRname in CRs:
        CRStatus[CRname]=0
    
    
    
    #print outfilename
    allquantities = MonoHbbQuantities(outfilename)
    allquantities.defineHisto()

#    for attr, value in allquantities.__dict__.iteritems():
#       print attr, value
#       if isinstance(value, float):
#          bbMET_tree.Branch('bbMETvariables',AddressOf(allquantities,'histo'),'histo/D')
    # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    # BTag Scale Factor Initialisation
    # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    othersys = ROOT.std.vector('string')()
    othersys.push_back('down')
    othersys.push_back('up')
    ## ThinJets
    if options.CSV:
        calib1 = ROOT.BTagCalibrationStandalone('csvv2', 'CSVv2_Moriond17_B_H.csv')
    if options.DeepCSV:
        calib1 = ROOT.BTagCalibrationStandalone('csvv2', 'DeepCSV_Moriond17_B_H.csv')
    reader1 = ROOT.BTagCalibrationStandaloneReader( 0, "central", othersys)    
    reader1.load(calib1, 0,  "comb" )  
    reader1.load(calib1, 1,  "comb" )  
    reader1.load(calib1, 2,  "incl" )  
    
#    h_total = TH1F('h_total','h_total',2,0,2)
#    h_total_mcweight = TH1F('h_total_mcweight','h_total_mcweight',2,0,2)
    
    
    
    for ievent in range(NEntries):
    #for ievent in range(501):
        
        ##
        sf_resolved1 = []
        sf_resolved2 = []
        sf_resolved3 = []
        #print "event number = ",ievent
        skimmedTree.GetEntry(ievent)
        
        ## Get all relevant branches
        
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Extract branches
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        run                        = skimmedTree.__getattr__('st_runId')
        lumi                       = skimmedTree.__getattr__('st_lumiSection')
        event                      = skimmedTree.__getattr__('st_eventId')

        #if event != 4126: continue                                
        #if lumi  != 42: continue                                
        if ievent%100==0: print (ievent)
        #trigName                   = skimmedTree.__getattr__('st_hlt_trigName')
        #trigResult                 = skimmedTree.__getattr__('st_hlt_trigResult')
        #filterName                 = skimmedTree.__getattr__('st_hlt_filterName')
        #filterResult               = skimmedTree.__getattr__('st_hlt_filterResult')
         
        pfMet                      = skimmedTree.__getattr__('st_pfMetCorrPt')
        pfMetPhi                   = skimmedTree.__getattr__('st_pfMetCorrPhi')
        
        nTHINJets                  = skimmedTree.__getattr__('st_THINnJet')
        thinjetP4                  = skimmedTree.__getattr__('st_THINjetP4')
        thinJetCSV                 = skimmedTree.__getattr__('st_THINjetCISVV2')
        #passThinJetLooseID         = skimmedTree.__getattr__('st_THINjetPassIDLoose')
        #passThinJetPUID            = skimmedTree.__getattr__('st_THINisPUJetID')
        THINjetHadronFlavor        = skimmedTree.__getattr__('st_THINjetHadronFlavor')
        thinjetNhadEF              = skimmedTree.__getattr__('st_THINjetNHadEF')
        thinjetChadEF              = skimmedTree.__getattr__('st_THINjetCHadEF')
        
        nTHINdeepCSVJets           = skimmedTree.__getattr__('st_AK4deepCSVnJet')
        thindeepCSVjetP4           = skimmedTree.__getattr__('st_AK4deepCSVjetP4')
        thinJetdeepCSV             = skimmedTree.__getattr__('st_AK4deepCSVjetDeepCSV_b')
        
        nEle                       = skimmedTree.__getattr__('st_nEle')
        eleP4                      = skimmedTree.__getattr__('st_eleP4')
        eleIsPassLoose             = skimmedTree.__getattr__('st_eleIsPassLoose')
        eleIsPassMedium            = skimmedTree.__getattr__('st_eleIsPassMedium')
        eleIsPassTight             = skimmedTree.__getattr__('st_eleIsPassTight')
        
        nMu                        = skimmedTree.__getattr__('st_nMu')
        muP4                       = skimmedTree.__getattr__('st_muP4')
        isLooseMuon                = skimmedTree.__getattr__('st_isLooseMuon')
        isMediumMuon               = skimmedTree.__getattr__('st_isMediumMuon')
        isTightMuon                = skimmedTree.__getattr__('st_isTightMuon')
        muChHadIso                 = skimmedTree.__getattr__('st_muChHadIso')
        muNeHadIso                 = skimmedTree.__getattr__('st_muNeHadIso')
        muGamIso                   = skimmedTree.__getattr__('st_muGamIso')
        muPUPt                     = skimmedTree.__getattr__('st_muPUPt')
        
        nTau                       = skimmedTree.__getattr__('st_HPSTau_n')
        tauP4                      = skimmedTree.__getattr__('st_HPSTau_4Momentum')
        #isDecayModeFinding         = skimmedTree.__getattr__('st_disc_decayModeFinding')
        #passLooseTauIso            = skimmedTree.__getattr__('st_disc_byLooseIsolationMVA3oldDMwLT')
        
        isData                     = skimmedTree.__getattr__('st_isData')
        mcWeight                   = skimmedTree.__getattr__('mcweight')
        pu_nTrueInt                = int(skimmedTree.__getattr__('st_pu_nTrueInt'))
        
        nGenPar                    = skimmedTree.__getattr__('st_nGenPar')
        genParId                   = skimmedTree.__getattr__('st_genParId')
        genMomParId                = skimmedTree.__getattr__('st_genMomParId')
        genParSt                   = skimmedTree.__getattr__('st_genParSt')
        genParP4                   = skimmedTree.__getattr__('st_genParP4')
        
        WenuRecoil                 = skimmedTree.__getattr__('WenuRecoil')
        Wenumass                   = skimmedTree.__getattr__('Wenumass')
        WenuPhi                    = skimmedTree.__getattr__('WenuPhi')
        WmunuRecoil                = skimmedTree.__getattr__('WmunuRecoil')
        Wmunumass                  = skimmedTree.__getattr__('Wmunumass')
        WmunuPhi                   = skimmedTree.__getattr__('WmunuPhi')
        ZeeRecoil                  = skimmedTree.__getattr__('ZeeRecoil')
        ZeeMass                    = skimmedTree.__getattr__('ZeeMass')
        ZeePhi                     = skimmedTree.__getattr__('ZeePhi')
        ZmumuRecoil                = skimmedTree.__getattr__('ZmumuRecoil')
        ZmumuMass                  = skimmedTree.__getattr__('ZmumuMass')
        ZmumuPhi                   = skimmedTree.__getattr__('ZmumuPhi')
        TOPRecoil                  = skimmedTree.__getattr__('TOPRecoil')
        TOPPhi                     = skimmedTree.__getattr__('TOPPhi')
        
        THINdeepCSVjetHadronFlavor = skimmedTree.__getattr__('st_THINjetHadronFlavor')
        thindeepCSVjetNhadEF = skimmedTree.__getattr__('st_THINjetNHadEF')
        thindeepCSVjetChadEF = skimmedTree.__getattr__('st_THINjetCHadEF')
        
        for trig in triglist:
            exec(trig+" = skimmedTree.__getattr__('st_"+trig+"')")
        
        try:
            MET_trig=skimmedTree.__getattr__('st_MET_trig')
            SE_trig=skimmedTree.__getattr__('st_SE_trig')
        except:
            MET_trig=True
            SE_trig=True
            if ievent==0: print "No MET_trig and SE_trig info available, the SkimmedTree seems to be from an old version. Proceeding with True for both."

#        HLT_IsoMu24                = skimmedTree.__getattr__('st_HLT_IsoMu20')     #Depreciated
#        HLT_Ele27_WPLoose_Gsf      = skimmedTree.__getattr__('st_HLT_Ele27_WPLoose_Gsf')
        
        
#        trig1 = CheckFilter(trigName, trigResult, 'HLT_PFMET170_') # added from  monojet
#        trig2 = CheckFilter(trigName, trigResult, 'HLT_PFMET170_NoiseCleaned')
#        trig3 = CheckFilter(trigName, trigResult, 'HLT_PFMET170_JetIdCleaned_v')
#        trig4 = CheckFilter(trigName, trigResult, 'HLT_PFMET170_HBHECleaned_v')
#        trig5 = CheckFilter(trigName, trigResult, 'HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v')
#        trig6 = CheckFilter(trigName, trigResult, 'HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v') #added from  tt+DM all hadronic analysis
#        trig7 = CheckFilter(trigName, trigResult, 'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v')
#        trig8 = CheckFilter(trigName, trigResult, 'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v')
#        trig9 = CheckFilter(trigName, trigResult, 'HLT_PFMET110_PFMHT110_')
#        trig10 = CheckFilter(trigName, trigResult, 'HLT_IsoMu24_v') #added from tt+DM all hadronic analysis
#        trig11 = CheckFilter(trigName, trigResult, 'HLT_IsoTkMu24_v') #added from tt+DM all hadronic analysis
#        trig12 = CheckFilter(trigName, trigResult, 'HLT_Ele27_WPTight_Gsf') #added from Siew Yan slides
##        trig13 = CheckFilter(trigName, trigResult, 'HLT_IsoMu20')   #Added from AN CR
#        trig13 = CheckFilter(trigName, trigResult, 'HLT_IsoMu20')    #Added from AN CR 2015
#        trig14 = CheckFilter(trigName, trigResult, 'HLT_Ele27_WPLoose_Gsf')   #Added from AN CR           
        
#        HLT_IsoMu24=trig10
#        HLT_IsoMu20=trig13
#        HLT_Ele27_WPLoose_Gsf=trig14
        
#        print MET_trig, SE_trig
        
        
        #**************************** REMEMBER TO CHANGE DEPENDING ON THE DATASET YOU ARE USING ********************************
        #==========================================================================
        #
        
#        if not MET_trig: continue                  # For signal and mu regions with MET dataset
        if not SE_trig: continue                    # For electron regions with SE dataset
        
        
        #============================ CAUTION =====================================
        #**************************************************************************

             
        jetSR1Info           = []
        jetSR2Info           = []
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # MC Weights ----------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        mcweight = 0.0 
        if isData==1:   mcweight =  1.0
        if not isData :
            if mcWeight<0:  mcweight = -1.0
            if mcWeight>0:  mcweight =  1.0
        
        
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        #
        
        for nc in range(nTHINJets):
            for nd in range(nTHINdeepCSVJets):
                if abs(thindeepCSVjetP4[nd].Pt() - thinjetP4[nc].Pt())<0.01:
                    THINdeepCSVjetHadronFlavor[nd] = THINjetHadronFlavor[nc]
                    thindeepCSVjetNhadEF[nd] = thinjetNhadEF[nc]
                    thindeepCSVjetChadEF[nd] = thinjetChadEF[nc]
                    
        CSVMWP=0.8484
        deepCSVMWP=0.6324
        if options.CSV:
            ##number of bjets without any selection
            mybjets=[]
            for nb in range(nTHINJets):
                if thinJetCSV[nb] > CSVMWP and abs(thinjetP4[nb].Eta())<2.4:
                    mybjets.append(nb)
            nBjets=len(mybjets)
        if options.DeepCSV:
            mybjets=[]
            for nb in range(nTHINdeepCSVJets):
                if thinJetdeepCSV[nb] > deepCSVMWP and abs(thindeepCSVjetP4[nb].Eta())<2.4:
                    mybjets.append(nb)
            nBjets=len(mybjets)

        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## PFMET Selection
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------        
        pfmetstatus = ( pfMet > 200.0 )   
#        if pfmetstatus == False : continue       
        if pfmetstatus: cutStatus['pfmet'] += 1
        #
        
#        if ZeeRecoil>200. or ZmumuRecoil > 200. or WenuRecoil > 200. or WmunuRecoil > 200. or TOPRecoil > 200.:
#            CRCutFlow['recoilprenjet']+=1 #----------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        ## Leptons Info
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ##Calculate Muon Relative PF isolation:
        MuIso = [((muChHadIso[imu]+ max(0., muNeHadIso[imu] + muGamIso[imu] - 0.5*muPUPt[imu]))/muP4[imu].Pt()) for imu in range(nMu)]
        
         # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
#        print (HLT_IsoMu24,HLT_Ele27_WPLoose_Gsf)
        
        
        myEles=[]
        myEleLooseID=[]
        myEleTightID=[]
        for iele in range(nEle):
            if eleP4[iele].Pt() < 10 : continue
            if abs(eleP4[iele].Eta()) >2.5: continue
            myEles.append(eleP4[iele])
            myEleLooseID.append(eleIsPassLoose[iele])
            myEleTightID.append(eleIsPassTight[iele])
        
        myMuos = []
        myMuLooseID=[]
        myMuTightID=[]
        myMuIso=[]
        for imu in range(nMu):
            if muP4[imu].Pt()<10 : continue
            if abs(muP4[imu].Eta()) > 2.4  : continue
#            #---Fake muon cleaner----
#            isClean=True
#            for iele in myEles[:]:
#                el_mu_dR=math.sqrt(  (  iele.Eta()-muP4[imu].Eta() )**2  + (  DeltaPhi(iele.Phi(),muP4[imu].Phi()) )**2 )
#                if el_mu_dR < 0.4:
#                    isClean=False
##                    myEles.remove(iele)     #Removes correspoding electron as well
#                    break
#            if not isClean: continue
#            ##---
            myMuos.append(muP4[imu])
            myMuLooseID.append(isLooseMuon[imu])
            myMuTightID.append(isTightMuon[imu])
            myMuIso.append(MuIso[imu])
            

        myTaus=[]
        for itau in range(nTau):
            if tauP4[itau].Pt()<18. : continue
            if abs(tauP4[itau].Eta())>2.3 : continue
            #---Fake tau cleaner----
            isClean=True
            for iele in myEles[:]:
                lep_tau_dR=math.sqrt(  (  iele.Eta()-tauP4[itau].Eta() )**2  + (  DeltaPhi(iele.Phi(),tauP4[itau].Phi()) )**2 )
                if lep_tau_dR < 0.4:
                    isClean=False
#                    myEles.remove(iele)     #Removes correspoding electron as well
                    break
            for imu in myMuos[:]:
                lep_tau_dR=math.sqrt(  (  imu.Eta()-tauP4[itau].Eta() )**2  + (  DeltaPhi(imu.Phi(),tauP4[itau].Phi()) )**2 )
                if lep_tau_dR < 0.4:
                    isClean=False
#                    myMuos.remove(imu)      #Removes correspoding muon as well
                    break
            if not isClean: continue
            ##---            
            myTaus.append(tauP4[itau])
        
        nUncleanEle=nEle
        nUncleanMu=nMu
        nUncleanTau=nTau
        
        nEle=len(myEles)
        nMu=len(myMuos)
        nTau=len(myTaus)
        
#        if nEle>0 and nMu<2 :
#            print ievent
#            print nEle
#            print "e:"
#            for iele in myEles:
#                print (iele.Eta(),iele.Phi())
#            print "mu:"
##            for imu in myMuos:
##                print (imu.Eta(),imu.Phi())
#            print
            
        
# --------------------------------------------------------------------------------------------------------------------------------------------------------
        #----------------------------------------------------------------------------------------------------------------------------------------------------------------         ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Sort jets 
        if options.CSV:
            if nTHINJets==0: continue
            
            alljetPT=[jet.Pt() for jet in thinjetP4]
            jetindex=[i for i in range(len(alljetPT))]
            
                        
            sortedjets=[jet for pt,jet in sorted(zip(alljetPT,thinjetP4), reverse=True)]      # This gives a list of jets with their pTs in descending order
            sortedindex=[jetindex for pt,jetindex in sorted(zip(alljetPT,jetindex), reverse=True)]     # Indices of jets in thinjetP4 in decscending order of jetPT
            
            j1=sortedjets[0]
            if nTHINJets>1: j2=sortedjets[1]
            if nTHINJets>2: j3=sortedjets[2]
            
            ifirstjet=sortedindex[0]
            if nTHINJets>1: isecondjet=sortedindex[1]
            if nTHINJets>2: ithirdjet=sortedindex[2]
            
        if options.DeepCSV:
            if nTHINdeepCSVJets==0: continue
            
            alljetPT=[jet.Pt() for jet in thindeepCSVjetP4]
            jetindex=[i for i in range(len(alljetPT))]
            
                        
            sortedjets=[jet for pt,jet in sorted(zip(alljetPT,thindeepCSVjetP4), reverse=True)]      # This gives a list of jets with their pTs in descending order
            sortedindex=[jetindex for pt,jetindex in sorted(zip(alljetPT,jetindex), reverse=True)]     # Indices of jets in thinjetP4 in decscending order of jetPT
            
            j1=sortedjets[0]
            if nTHINdeepCSVJets>1: j2=sortedjets[1]
            if nTHINdeepCSVJets>2: j3=sortedjets[2]
            
            ifirstjet=sortedindex[0]
            if nTHINdeepCSVJets>1: isecondjet=sortedindex[1]
            if nTHINdeepCSVJets>2: ithirdjet=sortedindex[2]
        
#        print alljetPT
#        print [jet.Pt() for jet in sortedjets]
#        print sortedindex
#        print
        
        ## 
# --------------------------------------------------------------------------------------------------------------------------------------------------------

        # SR start   
        
        ###
        #********
        #Part of data is blinded
        #********
        if ievent%20==0:
            keepevent=True
        else:
            keepevent=False
        #********************              REMEMBER TO UNBLIND AT SOME POINT!
        
        writeSR1=False
        writeSR2=False
        
        
        if options.CSV:
            SR1njetcond=False
            SR2njetcond=False
            SR1jetcond=False
            SR2jetcond=False        
                   
            
            if nEle+nMu+nTau==0:
                SRlepcond=True
            else:
                SRlepcond=False
            
            ## for SR1
             # 1 or 2 jets and 1 btagged 
            
            SRtrigstatus = HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v or HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v or HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v or HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v 
            
            if (nTHINJets == 1 or nTHINJets == 2) and pfmetstatus and SRlepcond and SRtrigstatus:
                #===CSVs before any selection===
                preselquantlist=AllQuantList.getPresel()        
                for quant in preselquantlist:
                    exec("allquantities."+quant+" = None")
                    
                allquantities.presel_jet1_csv_sr1=thinJetCSV[ifirstjet]	
                if nTHINJets>1: allquantities.presel_jet2_csv_sr1=thinJetCSV[isecondjet]
                allquantities.presel_jet1_chf_sr1=thinjetChadEF[ifirstjet]
                allquantities.presel_jet1_nhf_sr1=thinjetNhadEF[ifirstjet]		
                allquantities.FillPreSel()		
                #===
             
            if (nTHINJets == 1 or nTHINJets == 2) and nBjets==1 and SRtrigstatus: 
                SR1njetcond=True
                if pfmetstatus: cutStatusSR1['njet+nBjet'] +=1    
                if pfmetstatus: cutStatus['njet+nBjet'] += 1
                             
                SR1jetcond=True
                
                if j1.Pt() < 50.0: SR1jetcond=False
                if DeltaPhi(j1.Phi(),pfMetPhi) < 0.5: SR1jetcond=False           
                if thinjetNhadEF[ifirstjet] > 0.8 : SR1jetcond=False
                if thinjetChadEF[ifirstjet]< 0.1: SR1jetcond=False
                
                if SR1jetcond and pfmetstatus:
                    cutStatus['jet1'] += 1              # Lead jet satisfies required criteria
                    cutStatusSR1['jet1'] +=1
                
                if nTHINJets>1:                
                    if j2.Pt() < 30.0: SR1jetcond=False
                    if DeltaPhi(j2.Phi(),pfMetPhi) < 0.5: SR1jetcond=False
                
                if SR1jetcond and pfmetstatus:
                    cutStatus['jet2/3'] += 1           # Jet 2 satisfies the required criteria
                    cutStatusSR1['jet2'] +=1     

                          
                if SR1jetcond:
                    jet1pt = j1.Pt()
                    jet1phi = j1.Phi()
                    jet1eta = j1.Eta()
                    
                    if nTHINJets>1:
                        jet2pt = j2.Pt()
                        jet2phi = j2.Phi()
                        jet2eta = j2.Eta() 
                        jet2csv = thinJetCSV[isecondjet]
                        min_dPhi=min(DeltaPhi(j1.Phi(),pfMetPhi),DeltaPhi(j2.Phi(),pfMetPhi)) 
                    else:
                        jet2pt = None
                        jet2phi = None
                        jet2eta = None
                        jet2csv = None    
                        min_dPhi=DeltaPhi(j1.Phi(),pfMetPhi)                         
                    
                    jet1csv = thinJetCSV[ifirstjet]              

                    jetSR1Info.append([jet1pt,jet1eta,jet1phi,jet1csv])
                    jetSR1Info.append([jet2pt,jet2eta,jet2phi,jet2csv])
                    jetSR1Info.append(min_dPhi)
                    jetSR1Info.append(pfMet)
                    jetSR1Info.append(thinjetNhadEF[ifirstjet])
                    jetSR1Info.append(thinjetChadEF[ifirstjet])
                    writeSR1=True
              
         
         ## for SR2
            # 3 jets and 2 btagged 
            
            if (nTHINJets == 2 or nTHINJets == 3) and pfmetstatus and SRlepcond and SRtrigstatus:
                #===CSVs before any selection===	
                preselquantlist=AllQuantList.getPresel()        
                for quant in preselquantlist:
                    exec("allquantities."+quant+" = None")	
                allquantities.presel_jet1_csv_sr2=thinJetCSV[ifirstjet]
                allquantities.presel_jet2_csv_sr2=thinJetCSV[isecondjet]
                if nTHINJets>2: allquantities.presel_jet3_csv_sr2=thinJetCSV[ithirdjet]
                allquantities.presel_jet1_chf_sr2=thinjetChadEF[ifirstjet]
                allquantities.presel_jet1_nhf_sr2=thinjetNhadEF[ifirstjet]
                allquantities.FillPreSel()		
                #===
            
            if (nTHINJets == 2 or nTHINJets == 3) and nBjets==2 and SRtrigstatus:
                SR2njetcond=True
                if pfmetstatus: cutStatusSR2['njet+nBjet'] +=1
                if pfmetstatus: cutStatus['njet+nBjet'] += 1
                
                SR2jetcond=True
                
                if j1.Pt() < 50.0: SR2jetcond=False
                if DeltaPhi(j1.Phi(),pfMetPhi) < 0.5: SR2jetcond=False           
                if thinjetNhadEF[ifirstjet] > 0.8 : SR2jetcond=False
                if thinjetChadEF[ifirstjet]< 0.1: SR2jetcond=False
                
                if SR2jetcond and pfmetstatus:
                    cutStatus['jet1'] += 1              # Lead jet satisfies required criteria            
                    cutStatusSR2['jet1'] += 1
                
                if j2.Pt() < 50.0: SR2jetcond=False
                if DeltaPhi(j2.Phi(),pfMetPhi) < 0.5: SR2jetcond=False
                
                if SR2jetcond and pfmetstatus:
                    cutStatusSR2['jet2'] += 1
                
                if nTHINJets>2: 
                    if j3.Pt() < 30.0: SR2jetcond=False
                    if DeltaPhi(j3.Phi(),pfMetPhi) < 0.5: SR2jetcond=False
                
                if SR2jetcond and pfmetstatus:
                    cutStatusSR2['jet3'] += 1
                    cutStatus['jet2/3'] += 1           # The jets 2 and 3 satisfy the required criteria
                
                if SR2jetcond:
                    jet1pt = j1.Pt()
                    jet1phi = j1.Phi()
                    jet1eta = j1.Eta()
                    
                    jet2pt = j2.Pt()
                    jet2phi = j2.Phi()
                    jet2eta = j2.Eta()
                    
                    if nTHINJets>2: 
                        jet3pt = j3.Pt()
                        jet3phi = j3.Phi()
                        jet3eta = j3.Eta()
                        jet3csv = thinJetCSV[ithirdjet]
                        min_dPhi=min(DeltaPhi(j1.Phi(),pfMetPhi),DeltaPhi(j2.Phi(),pfMetPhi),DeltaPhi(j3.Phi(),pfMetPhi))
                    else:
                        jet3pt = None
                        jet3phi = None
                        jet3eta = None
                        jet3csv = None
                        min_dPhi=min(DeltaPhi(j1.Phi(),pfMetPhi),DeltaPhi(j2.Phi(),pfMetPhi))
                    
                    jet1csv = thinJetCSV[ifirstjet]
                    jet2csv = thinJetCSV[isecondjet]

                    jetSR2Info.append([jet1pt,jet1eta,jet1phi,jet1csv])
                    jetSR2Info.append([jet2pt,jet2eta,jet2phi,jet2csv])
                    jetSR2Info.append([jet3pt,jet3eta,jet3phi,jet3csv])
                    jetSR2Info.append(min_dPhi)
                    jetSR2Info.append(pfMet)
                    jetSR2Info.append(thinjetNhadEF[ifirstjet])
                    jetSR2Info.append(thinjetChadEF[ifirstjet])
                    writeSR2=True
                
            if pfmetstatus and SRlepcond and SR1jetcond:
                cutStatus['lep'] += 1
                cutStatusSR1['lep'] += 1
                
            if pfmetstatus and SRlepcond and SR2jetcond:
                cutStatus['lep'] += 1
                cutStatusSR2['lep'] += 1
            
        if options.DeepCSV:
            SR1njetcond=False
            SR2njetcond=False
            SR1jetcond=False
            SR2jetcond=False        
                   
            
            if nEle+nMu+nTau==0:
                SRlepcond=True
            else:
                SRlepcond=False
            
            ## for SR1
             # 1 or 2 jets and 1 btagged 
            
            SRtrigstatus = True#HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v or HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v
            
            if (nTHINdeepCSVJets == 1 or nTHINdeepCSVJets == 2) and pfmetstatus and SRlepcond and SRtrigstatus:
                #===CSVs before any selection===
                preselquantlist=AllQuantList.getPresel()        
                for quant in preselquantlist:
                    exec("allquantities."+quant+" = None")
                    
                allquantities.presel_jet1_deepcsv_sr1=thinJetdeepCSV[ifirstjet]	
                if nTHINdeepCSVJets>1: allquantities.presel_jet2_deepcsv_sr1=thinJetdeepCSV[isecondjet]
                allquantities.presel_jet1_chf_sr1=thindeepCSVjetChadEF[ifirstjet]
                allquantities.presel_jet1_nhf_sr1=thindeepCSVjetNhadEF[ifirstjet]		
                allquantities.FillPreSel()		
                #===
             
            if (nTHINdeepCSVJets == 1 or nTHINdeepCSVJets == 2) and nBjets==1 and SRtrigstatus: 
                SR1njetcond=True
                if pfmetstatus: cutStatusSR1['njet+nBjet'] +=1    
                if pfmetstatus: cutStatus['njet+nBjet'] += 1
                             
                SR1jetcond=True
                
                if j1.Pt() < 50.0: SR1jetcond=False
                if DeltaPhi(j1.Phi(),pfMetPhi) < 0.5: SR1jetcond=False           
                if thindeepCSVjetNhadEF[ifirstjet] > 0.8 : SR1jetcond=False
                if thindeepCSVjetChadEF[ifirstjet]< 0.1: SR1jetcond=False
                
                if SR1jetcond and pfmetstatus:
                    cutStatus['jet1'] += 1              # Lead jet satisfies required criteria
                    cutStatusSR1['jet1'] +=1
                
                if nTHINdeepCSVJets>1:                
                    if j2.Pt() < 30.0: SR1jetcond=False
                    if DeltaPhi(j2.Phi(),pfMetPhi) < 0.5: SR1jetcond=False
                
                if SR1jetcond and pfmetstatus:
                    cutStatus['jet2/3'] += 1           # Jet 2 satisfies the required criteria
                    cutStatusSR1['jet2'] +=1     

                          
                if SR1jetcond:
                    jet1pt = j1.Pt()
                    jet1phi = j1.Phi()
                    jet1eta = j1.Eta()
                    
                    if nTHINdeepCSVJets>1:
                        jet2pt = j2.Pt()
                        jet2phi = j2.Phi()
                        jet2eta = j2.Eta() 
                        jet2deepcsv = thinJetdeepCSV[isecondjet]
                        min_dPhi=min(DeltaPhi(j1.Phi(),pfMetPhi),DeltaPhi(j2.Phi(),pfMetPhi)) 
                    else:
                        jet2pt = None
                        jet2phi = None
                        jet2eta = None
                        jet2deepcsv = None    
                        min_dPhi=DeltaPhi(j1.Phi(),pfMetPhi)                         
                    
                    jet1deepcsv = thinJetdeepCSV[ifirstjet]              

                    jetSR1Info.append([jet1pt,jet1eta,jet1phi,jet1deepcsv])
                    jetSR1Info.append([jet2pt,jet2eta,jet2phi,jet2deepcsv])
                    jetSR1Info.append(min_dPhi)
                    jetSR1Info.append(pfMet)
                    jetSR1Info.append(thindeepCSVjetNhadEF[ifirstjet])
                    jetSR1Info.append(thindeepCSVjetChadEF[ifirstjet])
                    writeSR1=True
         
         ## for SR2
            # 3 jets and 2 btagged 
            
            if (nTHINdeepCSVJets == 2 or nTHINdeepCSVJets == 3) and pfmetstatus and SRlepcond and SRtrigstatus:
                #===CSVs before any selection===	
                preselquantlist=AllQuantList.getPresel()        
                for quant in preselquantlist:
                    exec("allquantities."+quant+" = None")	
                allquantities.presel_jet1_deepcsv_sr2=thinJetdeepCSV[ifirstjet]
                allquantities.presel_jet2_deepcsv_sr2=thinJetdeepCSV[isecondjet]
                if nTHINdeepCSVJets>2: allquantities.presel_jet3_deepcsv_sr2=thinJetdeepCSV[ithirdjet]
                allquantities.presel_jet1_chf_sr2=thindeepCSVjetChadEF[ifirstjet]
                allquantities.presel_jet1_nhf_sr2=thindeepCSVjetNhadEF[ifirstjet]
                allquantities.FillPreSel()		
                #===
            
            if (nTHINdeepCSVJets == 2 or nTHINdeepCSVJets == 3) and nBjets==2 and SRtrigstatus:
                SR2njetcond=True
                if pfmetstatus: cutStatusSR2['njet+nBjet'] +=1
                if pfmetstatus: cutStatus['njet+nBjet'] += 1
                
                SR2jetcond=True
                
                if j1.Pt() < 50.0: SR2jetcond=False
                if DeltaPhi(j1.Phi(),pfMetPhi) < 0.5: SR2jetcond=False           
                if thindeepCSVjetNhadEF[ifirstjet] > 0.8 : SR2jetcond=False
                if thindeepCSVjetChadEF[ifirstjet]< 0.1: SR2jetcond=False
                
                if SR2jetcond and pfmetstatus:
                    cutStatus['jet1'] += 1              # Lead jet satisfies required criteria            
                    cutStatusSR2['jet1'] += 1
                
                if j2.Pt() < 50.0: SR2jetcond=False
                if DeltaPhi(j2.Phi(),pfMetPhi) < 0.5: SR2jetcond=False
                
                if SR2jetcond and pfmetstatus:
                    cutStatusSR2['jet2'] += 1
                
                if nTHINdeepCSVJets>2: 
                    if j3.Pt() < 30.0: SR2jetcond=False
                    if DeltaPhi(j3.Phi(),pfMetPhi) < 0.5: SR2jetcond=False
                
                if SR2jetcond and pfmetstatus:
                    cutStatusSR2['jet3'] += 1
                    cutStatus['jet2/3'] += 1           # The jets 2 and 3 satisfy the required criteria
                
                if SR2jetcond:
                    jet1pt = j1.Pt()
                    jet1phi = j1.Phi()
                    jet1eta = j1.Eta()
                    
                    jet2pt = j2.Pt()
                    jet2phi = j2.Phi()
                    jet2eta = j2.Eta()
                    
                    if nTHINdeepCSVJets>2: 
                        jet3pt = j3.Pt()
                        jet3phi = j3.Phi()
                        jet3eta = j3.Eta()
                        jet3deepcsv = thinJetdeepCSV[ithirdjet]
                        min_dPhi=min(DeltaPhi(j1.Phi(),pfMetPhi),DeltaPhi(j2.Phi(),pfMetPhi),DeltaPhi(j3.Phi(),pfMetPhi))
                    else:
                        jet3pt = None
                        jet3phi = None
                        jet3eta = None
                        jet3deepcsv = None
                        min_dPhi=min(DeltaPhi(j1.Phi(),pfMetPhi),DeltaPhi(j2.Phi(),pfMetPhi))
                    
                    jet1deepcsv = thinJetdeepCSV[ifirstjet]
                    jet2deepcsv = thinJetdeepCSV[isecondjet]

                    jetSR2Info.append([jet1pt,jet1eta,jet1phi,jet1deepcsv])
                    jetSR2Info.append([jet2pt,jet2eta,jet2phi,jet2deepcsv])
                    jetSR2Info.append([jet3pt,jet3eta,jet3phi,jet3deepcsv])
                    jetSR2Info.append(min_dPhi)
                    jetSR2Info.append(pfMet)
                    jetSR2Info.append(thindeepCSVjetNhadEF[ifirstjet])
                    jetSR2Info.append(thindeepCSVjetChadEF[ifirstjet])
                    writeSR2=True
                
            if pfmetstatus and SRlepcond and SR1jetcond:
                cutStatus['lep'] += 1
                cutStatusSR1['lep'] += 1
                
            if pfmetstatus and SRlepcond and SR2jetcond:
                cutStatus['lep'] += 1
                cutStatusSR2['lep'] += 1
#        CRCutFlow['njet+nBjet']+=1
        
#        if ZeeRecoil>200. or ZmumuRecoil > 200. or WenuRecoil > 200. or WmunuRecoil > 200. or TOPRecoil > 200.:
#            CRCutFlow['recoilpostnjet']+=1
        #----------------------------------------------------------------------------------------------------------------------------------------------------------------
#        if SR1jetcond==False and SR2jetcond==False:
#            continue

# --------------------------------------------------------------------------------------------------------------------------------------------------------

        #Control Region with only category cuts
        
        preselquantlist=AllQuantList.getPresel()
        
        for quant in preselquantlist:
            exec("allquantities."+quant+" = None")
            
              
        regquants=AllQuantList.getRegionQuants()
        
        for quant in regquants:
            exec("allquantities."+quant+" = None")
        

        
        ####new conds
        jetcond=True
        SR2jet2=True
        if j1.Pt() < 50.0: jetcond=False
        if options.CSV:
            if thinjetNhadEF[ifirstjet] > 0.8 : jetcond=False
            if thinjetChadEF[ifirstjet]< 0.1: jetcond=False
        if options.DeepCSV:
            if thindeepCSVjetNhadEF[ifirstjet] > 0.8 : jetcond=False
            if thindeepCSVjetChadEF[ifirstjet]< 0.1: jetcond=False
        #       
        if options.CSV:
            if nTHINJets>=2:
                if j2.Pt() < 30.0: jetcond=False            
                
                if j2.Pt() > 50.0:
                    SR2jet2=True
                else:
                    SR2jet2=False
                    
        #            if j2.Pt() > 50.0:
        #                SR2jet2=True
        #            else:
        #                SR2jet2=False
                        
                if nTHINJets>=3:            
                    if j3.Pt() < 30.0: jetcond=False
                
        if options.DeepCSV:
            if nTHINdeepCSVJets>=2:
                if j2.Pt() < 30.0: jetcond=False            
                
                if j2.Pt() > 50.0:
                    SR2jet2=True
                else:
                    SR2jet2=False
                    
    #            if j2.Pt() > 50.0:
    #                SR2jet2=True
    #            else:
    #                SR2jet2=False
                    
            if nTHINdeepCSVJets>=3:
                if j3.Pt() < 30.0: jetcond=False

#        if jetcond: CRCutFlow['jetcond']+=1
#        ### Experimental: First 1/2/3 jets alone satisfy nBjet condition: Doesn't make any positive difference
#        
#        SR1bjetcond=False
#        SR2bjetcond=False
#        
#        if nTHINJets==1 and thinJetCSV[0]>CSVMWP: SR1bjetcond = True            
#        
#        if nTHINJets>=2:
#            nbjetin2=0
#            for ijet in [ifirstjet,isecondjet]:
#                if thinJetCSV[ijet]>CSVMWP:
#                    nbjetin2 += 1
#            if nbjetin2 == 1: SR1bjetcond = True
#            if nbjetin2 == 2: SR2bjetcond = True
#        
#        if nTHINJets>=3:
#            nbjetin3=0
#            for ijet in [ifirstjet,isecondjet,ithirdjet]:
#                if thinJetCSV[ijet]>CSVMWP:
#                    nbjetin3 += 1
#            if nbjetin3 == 2: SR2bjetcond = True
#            


# -------------------------------------------
# Z CR
# -------------------------------------------
            
        #Z CR specific bools
        
        ZdPhicond=True
        
        if ZeePhi>-10.:
            if DeltaPhi(j1.Phi(),Phi_mpi_pi(math.pi+ZeePhi)) < 0.5: ZdPhicond = False      #Added +pi to ZPhi to reverse an error in SkimTree which will be fixed in next iteration.
            if options.CSV:
                if nTHINJets>=2:
                    if DeltaPhi(j2.Phi(),Phi_mpi_pi(math.pi+ZeePhi)) < 0.5: ZdPhicond=False
                if nTHINJets>=3:
                    if DeltaPhi(j3.Phi(),Phi_mpi_pi(math.pi+ZeePhi)) < 0.5: ZdPhicond=False
            if options.DeepCSV:
                if nTHINdeepCSVJets>=2:
                    if DeltaPhi(j2.Phi(),Phi_mpi_pi(math.pi+ZeePhi)) < 0.5: ZdPhicond=False
                if nTHINdeepCSVJets>=3:
                    if DeltaPhi(j3.Phi(),Phi_mpi_pi(math.pi+ZeePhi)) < 0.5: ZdPhicond=False
        if ZmumuPhi>-10.:
            if DeltaPhi(j1.Phi(),Phi_mpi_pi(math.pi+ZmumuPhi)) < 0.5: ZdPhicond = False      
            if options.CSV:
                if nTHINJets>=2:
                    if DeltaPhi(j2.Phi(),Phi_mpi_pi(math.pi+ZmumuPhi)) < 0.5: ZdPhicond=False
                if nTHINJets>=3:
                    if DeltaPhi(j3.Phi(),Phi_mpi_pi(math.pi+ZmumuPhi)) < 0.5: ZdPhicond=False
            if options.DeepCSV:
                if nTHINdeepCSVJets>=2:
                    if DeltaPhi(j2.Phi(),Phi_mpi_pi(math.pi+ZmumuPhi)) < 0.5: ZdPhicond=False
                if nTHINdeepCSVJets>=3:
                    if DeltaPhi(j3.Phi(),Phi_mpi_pi(math.pi+ZmumuPhi)) < 0.5: ZdPhicond=False
#            
#        
#        ###

#        if ZdPhicond: CRCutFlow['ZdPhi']+=1

        if (HLT_Ele27_WPLoose_Gsf or HLT_Ele27_WPTight_Gsf):
            CR2e1bCutFlow['trig']+=1
            CR2e2bCutFlow['trig']+=1
            if nEle==2 and nMu==0:
                CR2e1bCutFlow['nlep']+=1
                CR2e2bCutFlow['nlep']+=1
                if myEles[0].Pt()>myEles[1].Pt():
                    iLeadLep=0
                    iSecondLep=1
                else:
                    iLeadLep=1
                    iSecondLep=0
                if myEles[iLeadLep].Pt() > 30. and myEleTightID[iLeadLep] and myEles[iSecondLep].Pt() > 10. and myEleLooseID[iSecondLep]:
                    CR2e1bCutFlow['lepconds']+=1
                    CR2e2bCutFlow['lepconds']+=1
                    if ZeeRecoil>200.:
                        CR2e1bCutFlow['recoil']+=1
                        CR2e2bCutFlow['recoil']+=1
                        if ZeeMass>70. and ZeeMass<110.:
                            CR2e1bCutFlow['mass']+=1  
                            CR2e2bCutFlow['mass']+=1   
                            if ZdPhicond:
                                CR2e1bCutFlow['dPhicond']+=1 
                                CR2e2bCutFlow['dPhicond']+=1                                
                                if nBjets==1:
                                    CR2e1bCutFlow['nbjets']+=1
                                    if jetcond:
                                        CR2e1bCutFlow['jetconds']+=1
                                if nBjets==2:
                                    CR2e2bCutFlow['nbjets']+=1
                                    if jetcond and SR2jet2:
                                        CR2e2bCutFlow['jetconds']+=1
        
         #2e, 1 b-tagged
        
        if nEle==2 and nMu==0 and (HLT_Ele27_WPLoose_Gsf or HLT_Ele27_WPTight_Gsf) and ZeeMass>70. and ZeeMass<110. and ZeeRecoil>200. and jetcond and ZdPhicond:
#            CRCutFlow['nlepcond']+=1
            alllepPT=[lep.Pt() for lep in myEles]
            lepindex=[i for i in range(len(myEles))]            
                        
            sortedleps=[lep for pt,lep in sorted(zip(alllepPT,myEles), reverse=True)]      # This gives a list of leps with their pTs in descending order
            sortedindex=[lepind for pt,lepind in sorted(zip(alllepPT,lepindex), reverse=True)]     # Indices of leps in thinjetP4 in decscending order of jetPT
            
            iLeadLep=sortedindex[0]
            iSecondLep=sortedindex[1]            
            
            if myEles[iLeadLep].Pt() > 30. and myEleTightID[iLeadLep] and myEles[iSecondLep].Pt() > 10. and myEleLooseID[iSecondLep]:            
            
                ZpT = math.sqrt( (myEles[iLeadLep].Px()+myEles[iSecondLep].Px())*(myEles[iLeadLep].Px()+myEles[iSecondLep].Px()) + (myEles[iLeadLep].Py()+myEles[iSecondLep].Py())*(myEles[iLeadLep].Py()+myEles[iSecondLep].Py()) )                       
                            
                if nBjets==1:
                    allquantities.reg_2e1b_Zmass = ZeeMass    
                    allquantities.reg_2e1b_ZpT=ZpT
                    
                    allquantities.reg_2e1b_hadrecoil = ZeeRecoil
                    allquantities.reg_2e1b_MET = pfMet
                    
                    allquantities.reg_2e1b_lep1_pT=myEles[iLeadLep].Pt()
                    allquantities.reg_2e1b_lep2_pT=myEles[iSecondLep].Pt()
                    
                    if options.CSV:
                        allquantities.reg_2e1b_jet1_pT=j1.Pt()
                        if nTHINJets>1: allquantities.reg_2e1b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_2e1b_jet1_eta=j1.Eta()
                        if nTHINJets>1: allquantities.reg_2e1b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_2e1b_jet1_csv = thinJetCSV[ifirstjet]
                        if nTHINJets>1: allquantities.reg_2e1b_jet2_csv = thinJetCSV[isecondjet]
                        allquantities.reg_2e1b_njet = nTHINJets
                        
                    if options.DeepCSV:
                        allquantities.reg_2e1b_jet1_pT=j1.Pt()
                        if nTHINdeepCSVJets>1: allquantities.reg_2e1b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_2e1b_jet1_eta=j1.Eta()
                        if nTHINdeepCSVJets>1: allquantities.reg_2e1b_jet2_eta=j2.Eta()
                        
                        allquantities.reg_2e1b_jet1_deepcsv = thinJetdeepCSV[ifirstjet]
                        if nTHINdeepCSVJets>1: allquantities.reg_2e1b_jet2_deepcsv = thinJetdeepCSV[isecondjet]
                        allquantities.reg_2e1b_njet = nTHINdeepCSVJets
                    
                    allquantities.reg_2e1b_ntau = nTau
                    allquantities.reg_2e1b_nele = nEle
                    allquantities.reg_2e1b_nmu = nMu
                    allquantities.reg_2e1b_nUncleanTau = nUncleanTau
                    allquantities.reg_2e1b_nUncleanEle = nUncleanEle
                    allquantities.reg_2e1b_nUncleanMu = nUncleanMu
                    
#                    #--- Tau cleaning---
#                    ncleanTau=0
#                    for itau in myTaus:
#                        lep1_tau_dR=math.sqrt(  (  myEles[iLeadLep].Eta()-itau.Eta() )**2  + (  DeltaPhi(myEles[iLeadLep].Phi(),itau.Phi()) )**2 )
#                        lep2_tau_dR=math.sqrt(  (  myEles[iSecondLep].Eta()-itau.Eta() )**2  + (  DeltaPhi(myEles[iSecondLep].Phi(),itau.Phi()) )**2 )
#                        allquantities.reg_2e1b_lep1_dR_tau=lep1_tau_dR
#                        allquantities.reg_2e1b_lep2_dR_tau=lep2_tau_dR
#                        allquantities.reg_2e1b_min_lep_dR_tau=min(lep1_tau_dR,lep2_tau_dR)
##                        print ievent
##                        print myEles[iLeadLep].Eta()
##                        print itau.Eta()
##                        print myEles[iLeadLep].Phi()
##                        print itau.Phi()
##                        print myEles[iSecondLep].Eta()
##                        print myEles[iSecondLep].Phi()
##                        print lep1_tau_dR
##                        print lep2_tau_dR
##                        print
#                        if lep1_tau_dR > 0.4 and lep2_tau_dR > 0.4: ncleanTau += 1 
#                    #---  
#                    allquantities.reg_2e1b_ntaucleaned = ncleanTau
                    
            #2e, 2 b-tagged   
                if nBjets==2 and SR2jet2:         
                    allquantities.reg_2e2b_Zmass = ZeeMass   
                    allquantities.reg_2e2b_ZpT=ZpT
                    
                    allquantities.reg_2e2b_hadrecoil = ZeeRecoil
                    allquantities.reg_2e2b_MET = pfMet
                    
                    allquantities.reg_2e2b_lep1_pT=myEles[iLeadLep].Pt()
                    allquantities.reg_2e2b_lep2_pT=myEles[iSecondLep].Pt()
                    
                    if options.CSV:
                        allquantities.reg_2e2b_jet1_pT=j1.Pt()
                        if nTHINJets>1: allquantities.reg_2e2b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_2e2b_jet1_eta=j1.Eta()
                        if nTHINJets>1: allquantities.reg_2e2b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_2e2b_jet1_csv = thinJetCSV[ifirstjet]
                        if nTHINJets>1: allquantities.reg_2e2b_jet2_csv = thinJetCSV[isecondjet]
                        
                        allquantities.reg_2e2b_njet = nTHINJets
                    
                    if options.DeepCSV:
                        allquantities.reg_2e2b_jet1_pT=j1.Pt()
                        if nTHINdeepCSVJets>1: allquantities.reg_2e2b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_2e2b_jet1_eta=j1.Eta()
                        if nTHINdeepCSVJets>1: allquantities.reg_2e2b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_2e2b_jet1_deepcsv = thinJetdeepCSV[ifirstjet]
                        if nTHINdeepCSVJets>1: allquantities.reg_2e2b_jet2_deepcsv = thinJetdeepCSV[isecondjet]
                        
                        allquantities.reg_2e2b_njet = nTHINdeepCSVJets
                    
                    allquantities.reg_2e2b_ntau = nTau
                    allquantities.reg_2e2b_nele = nEle
                    allquantities.reg_2e2b_nmu = nMu
                    allquantities.reg_2e2b_nUncleanTau = nUncleanTau
                    allquantities.reg_2e2b_nUncleanEle = nUncleanEle
                    allquantities.reg_2e2b_nUncleanMu = nUncleanMu
                    
#                    #--- Tau cleaning---
#                    ncleanTau=0
#                    for itau in myTaus:
#                        lep1_tau_dR=math.sqrt(  (  myEles[iLeadLep].Eta()-itau.Eta() )**2  + (  DeltaPhi(myEles[iLeadLep].Phi(),itau.Phi()) )**2 )
#                        lep2_tau_dR=math.sqrt(  (  myEles[iSecondLep].Eta()-itau.Eta() )**2  + (  DeltaPhi(myEles[iSecondLep].Phi(),itau.Phi()) )**2 )
#                        allquantities.reg_2e2b_lep1_dR_tau=lep1_tau_dR
#                        allquantities.reg_2e2b_lep2_dR_tau=lep2_tau_dR
#                        allquantities.reg_2e2b_min_lep_dR_tau=min(lep1_tau_dR,lep2_tau_dR)
#                        if lep1_tau_dR > 0.4 and lep2_tau_dR > 0.4: ncleanTau += 1 
#                    #---  
#                    allquantities.reg_2e2b_ntaucleaned = ncleanTau


        if ((UnPrescaledIsoMu20 and HLT_IsoMu20) or HLT_IsoMu24_v or HLT_IsoTkMu24_v):
            CR2mu1bCutFlow['trig']+=1
            CR2mu2bCutFlow['trig']+=1
            if nMu==2 and nEle==0:
                CR2mu1bCutFlow['nlep']+=1
                CR2mu2bCutFlow['nlep']+=1
                if myMuos[0].Pt()>myMuos[1].Pt():
                    iLeadLep=0
                    iSecondLep=1
                else:
                    iLeadLep=1
                    iSecondLep=0
                if myMuos[iLeadLep].Pt() > 30. and myMuTightID[iLeadLep] and myMuIso[iLeadLep]<0.15 and myMuos[iSecondLep].Pt() > 10. and myMuLooseID[iSecondLep] and myMuIso[iSecondLep]<0.25:
                    CR2mu1bCutFlow['lepconds']+=1
                    CR2mu2bCutFlow['lepconds']+=1
                    if ZmumuRecoil>200.:
                        CR2mu1bCutFlow['recoil']+=1
                        CR2mu2bCutFlow['recoil']+=1
                        if ZmumuMass>70. and ZmumuMass<110.:
                            CR2mu1bCutFlow['mass']+=1  
                            CR2mu2bCutFlow['mass']+=1   
                            if ZdPhicond:
                                CR2mu1bCutFlow['dPhicond']+=1 
                                CR2mu2bCutFlow['dPhicond']+=1                                
                                if nBjets==1:
                                    CR2mu1bCutFlow['nbjets']+=1
                                    if jetcond:
                                        CR2mu1bCutFlow['jetconds']+=1
                                if nBjets==2:
                                    CR2mu2bCutFlow['nbjets']+=1
                                    if jetcond and SR2jet2:
                                        CR2mu2bCutFlow['jetconds']+=1



        #2mu, 1 b-tagged  
        if nMu==2 and nEle==0 and ((UnPrescaledIsoMu20 and HLT_IsoMu20) or HLT_IsoMu24_v or HLT_IsoTkMu24_v) and ZmumuMass>70. and ZmumuMass<110. and ZmumuRecoil>200. and jetcond and ZdPhicond:
#            CRCutFlow['nlepcond']+=1
            alllepPT=[lep.Pt() for lep in myMuos]
            lepindex=[i for i in range(len(myMuos))]            
                        
            sortedleps=[lep for pt,lep in sorted(zip(alllepPT,myMuos), reverse=True)]      # This gives a list of leps with their pTs in descending order
            sortedindex=[lepind for pt,lepind in sorted(zip(alllepPT,lepindex), reverse=True)]     # Indices of leps in thinjetP4 in decscending order of jetPT
            
            iLeadLep=sortedindex[0]
            iSecondLep=sortedindex[1]
                
            if myMuos[iLeadLep].Pt() > 30. and myMuTightID[iLeadLep] and myMuIso[iLeadLep]<0.15 and myMuos[iSecondLep].Pt() > 10. and myMuLooseID[iSecondLep] and myMuIso[iSecondLep]<0.25:  
                        
                ZpT = math.sqrt( (myMuos[iLeadLep].Px()+myMuos[iSecondLep].Px())*(myMuos[iLeadLep].Px()+myMuos[iSecondLep].Px()) + (myMuos[iLeadLep].Py()+myMuos[iSecondLep].Py())*(myMuos[iLeadLep].Py()+myMuos[iSecondLep].Py()) )                
                            
                if  nBjets==1:         
                    allquantities.reg_2mu1b_Zmass = ZmumuMass    
                    allquantities.reg_2mu1b_ZpT=ZpT
                    
                    allquantities.reg_2mu1b_hadrecoil = ZmumuRecoil
                    allquantities.reg_2mu1b_MET = pfMet
                    
                    allquantities.reg_2mu1b_lep1_pT=myMuos[iLeadLep].Pt()
                    allquantities.reg_2mu1b_lep2_pT=myMuos[iSecondLep].Pt()
                    
                    allquantities.reg_2mu1b_lep1_iso=myMuIso[iLeadLep]
                    allquantities.reg_2mu1b_lep2_iso=myMuIso[iSecondLep]
                    
                    if options.CSV:
                        allquantities.reg_2mu1b_jet1_pT=j1.Pt()
                        if nTHINJets>1: allquantities.reg_2mu1b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_2mu1b_jet1_eta=j1.Eta()
                        if nTHINJets>1: allquantities.reg_2mu1b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_2mu1b_jet1_csv = thinJetCSV[ifirstjet]
                        if nTHINJets>1: allquantities.reg_2mu1b_jet2_csv = thinJetCSV[isecondjet]
                        
                        allquantities.reg_2mu1b_njet = nTHINJets
                    
                    if options.DeepCSV:
                        allquantities.reg_2mu1b_jet1_pT=j1.Pt()
                        if nTHINdeepCSVJets>1: allquantities.reg_2mu1b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_2mu1b_jet1_eta=j1.Eta()
                        if nTHINdeepCSVJets>1: allquantities.reg_2mu1b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_2mu1b_jet1_deepcsv = thinJetdeepCSV[ifirstjet]
                        if nTHINdeepCSVJets>1: allquantities.reg_2mu1b_jet2_deepcsv = thinJetdeepCSV[isecondjet]
                        
                        allquantities.reg_2mu1b_njet = nTHINdeepCSVJets
                    
                    allquantities.reg_2mu1b_ntau = nTau
                    allquantities.reg_2mu1b_nele = nEle
                    allquantities.reg_2mu1b_nmu = nMu
                    allquantities.reg_2mu1b_nUncleanTau = nUncleanTau
                    allquantities.reg_2mu1b_nUncleanEle = nUncleanEle
                    allquantities.reg_2mu1b_nUncleanMu = nUncleanMu
                    
#                    #--- Tau cleaning---
#                    ncleanTau=0
#                    for itau in myTaus:
#                        lep1_tau_dR=math.sqrt(  (  myMuos[iLeadLep].Eta()-itau.Eta() )**2  + (  DeltaPhi(myMuos[iLeadLep].Phi(),itau.Phi()) )**2 )
#                        lep2_tau_dR=math.sqrt(  (  myMuos[iSecondLep].Eta()-itau.Eta() )**2  + (  DeltaPhi(myMuos[iSecondLep].Phi(),itau.Phi()) )**2 )
#                        allquantities.reg_2mu1b_lep1_dR_tau=lep1_tau_dR
#                        allquantities.reg_2mu1b_lep2_dR_tau=lep2_tau_dR
#                        allquantities.reg_2mu1b_min_lep_dR_tau=min(lep1_tau_dR,lep2_tau_dR)
#                        if lep1_tau_dR > 0.4 and lep2_tau_dR > 0.4: ncleanTau += 1 
#                    #---  
#                    allquantities.reg_2mu1b_ntaucleaned = ncleanTau
                    
            #2mu, 2 b-tagged        
                if  nBjets==2 and SR2jet2:
                    allquantities.reg_2mu2b_Zmass = ZmumuMass   
                    allquantities.reg_2mu2b_ZpT=ZpT
                    
                    allquantities.reg_2mu2b_hadrecoil = ZmumuRecoil
                    allquantities.reg_2mu2b_MET = pfMet
                    
                    allquantities.reg_2mu2b_lep1_pT=myMuos[iLeadLep].Pt()
                    allquantities.reg_2mu2b_lep2_pT=myMuos[iSecondLep].Pt()
                    
                    allquantities.reg_2mu2b_lep1_iso=myMuIso[iLeadLep]
                    allquantities.reg_2mu2b_lep2_iso=myMuIso[iSecondLep]
                    
                    if options.CSV:
                        allquantities.reg_2mu2b_jet1_pT=j1.Pt()
                        if nTHINJets>1: allquantities.reg_2mu2b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_2mu2b_jet1_eta=j1.Eta()
                        if nTHINJets>1: allquantities.reg_2mu2b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_2mu2b_jet1_csv = thinJetCSV[ifirstjet]
                        if nTHINJets>1: allquantities.reg_2mu2b_jet2_csv = thinJetCSV[isecondjet]
                        
                        allquantities.reg_2mu2b_njet = nTHINJets
                        
                    if options.DeepCSV:
                        allquantities.reg_2mu2b_jet1_pT=j1.Pt()
                        if nTHINdeepCSVJets>1: allquantities.reg_2mu2b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_2mu2b_jet1_eta=j1.Eta()
                        if nTHINdeepCSVJets>1: allquantities.reg_2mu2b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_2mu2b_jet1_deepcsv = thinJetdeepCSV[ifirstjet]
                        if nTHINdeepCSVJets>1: allquantities.reg_2mu2b_jet2_deepcsv = thinJetdeepCSV[isecondjet]
                        
                        allquantities.reg_2mu2b_njet = nTHINdeepCSVJets
                        
                    allquantities.reg_2mu2b_ntau = nTau
                    allquantities.reg_2mu2b_nele = nEle
                    allquantities.reg_2mu2b_nmu = nMu
                    allquantities.reg_2mu2b_nUncleanTau = nUncleanTau
                    allquantities.reg_2mu2b_nUncleanEle = nUncleanEle
                    allquantities.reg_2mu2b_nUncleanMu = nUncleanMu
                    
#                    #--- Tau cleaning---
#                    ncleanTau=0
#                    for itau in myTaus:
#                        lep1_tau_dR=math.sqrt(  (  myMuos[iLeadLep].Eta()-itau.Eta() )**2  + (  DeltaPhi(myMuos[iLeadLep].Phi(),itau.Phi()) )**2 )
#                        lep2_tau_dR=math.sqrt(  (  myMuos[iSecondLep].Eta()-itau.Eta() )**2  + (  DeltaPhi(myMuos[iSecondLep].Phi(),itau.Phi()) )**2 )
#                        allquantities.reg_2mu2b_lep1_dR_tau=lep1_tau_dR
#                        allquantities.reg_2mu2b_lep2_dR_tau=lep2_tau_dR
#                        allquantities.reg_2mu2b_min_lep_dR_tau=min(lep1_tau_dR,lep2_tau_dR)
#                        if lep1_tau_dR > 0.4 and lep2_tau_dR > 0.4: ncleanTau += 1 
#                    #---  
#                    allquantities.reg_2mu2b_ntaucleaned = ncleanTau
        
# -------------------------------------------
# W CR
# -------------------------------------------
            
        #W CR specific bools
        
        WdPhicond=True
        
        if WenuPhi>-10.:
            if DeltaPhi(j1.Phi(),Phi_mpi_pi(math.pi+WenuPhi)) < 0.5: WdPhicond = False      #Added +pi to ZPhi to reverse an error in SkimTree which will be fixed in next iteration.
            if options.CSV:
                if nTHINJets>=2:
                    if DeltaPhi(j2.Phi(),Phi_mpi_pi(math.pi+WenuPhi)) < 0.5: WdPhicond=False
                if nTHINJets>=3:
                    if DeltaPhi(j3.Phi(),Phi_mpi_pi(math.pi+WenuPhi)) < 0.5: WdPhicond=False
            if options.DeepCSV:
                if nTHINdeepCSVJets>=2:
                    if DeltaPhi(j2.Phi(),Phi_mpi_pi(math.pi+WenuPhi)) < 0.5: WdPhicond=False
                if nTHINdeepCSVJets>=3:
                    if DeltaPhi(j3.Phi(),Phi_mpi_pi(math.pi+WenuPhi)) < 0.5: WdPhicond=False
        if WmunuPhi>-10.:
            if DeltaPhi(j1.Phi(),Phi_mpi_pi(math.pi+WmunuPhi)) < 0.5: WdPhicond = False      
            if options.CSV:
                if nTHINJets>=2:
                    if DeltaPhi(j2.Phi(),Phi_mpi_pi(math.pi+WmunuPhi)) < 0.5: WdPhicond=False
                if nTHINJets>=3:
                    if DeltaPhi(j3.Phi(),Phi_mpi_pi(math.pi+WmunuPhi)) < 0.5: WdPhicond=False
            if options.DeepCSV:
                if nTHINdeepCSVJets>=2:
                    if DeltaPhi(j2.Phi(),Phi_mpi_pi(math.pi+WmunuPhi)) < 0.5: WdPhicond=False
                if nTHINdeepCSVJets>=3:
                    if DeltaPhi(j3.Phi(),Phi_mpi_pi(math.pi+WmunuPhi)) < 0.5: WdPhicond=False
                    
  
  #Cutflow                  
        if (HLT_Ele27_WPLoose_Gsf or HLT_Ele27_WPTight_Gsf):
            CR1e1bCutFlow['trig']+=1
            CR1e2bCutFlow['trig']+=1
            if nEle==1 and nMu==0:
                CR1e1bCutFlow['nlep']+=1
                CR1e2bCutFlow['nlep']+=1
                if myEles[0].Pt() > 30. and myEleTightID[0]:
                    CR1e1bCutFlow['lepconds']+=1
                    CR1e2bCutFlow['lepconds']+=1
                    if WenuRecoil>200.:
                        CR1e1bCutFlow['recoil']+=1
                        CR1e2bCutFlow['recoil']+=1
                        if Wenumass>50. and  Wenumass<160.:
                            CR1e1bCutFlow['mass']+=1  
                            CR1e2bCutFlow['mass']+=1   
                            if WdPhicond:
                                CR1e1bCutFlow['dPhicond']+=1 
                                CR1e2bCutFlow['dPhicond']+=1                                
                                if nBjets==1:
                                    CR1e1bCutFlow['nbjets']+=1
                                    if jetcond:
                                        CR1e1bCutFlow['jetconds']+=1
                                if nBjets==2:
                                    CR1e2bCutFlow['nbjets']+=1
                                    if jetcond and SR2jet2:
                                        CR1e2bCutFlow['jetconds']+=1            
        
        #1e, 1 b-tagged
        if nEle==1 and nMu==0 and (HLT_Ele27_WPLoose_Gsf or HLT_Ele27_WPTight_Gsf) and WenuRecoil>200. and jetcond and WdPhicond and Wenumass>50. and Wenumass<160.:
#            
            iLeadLep=0          
            
            if myEles[iLeadLep].Pt() > 30. and myEleTightID[iLeadLep]:             
                
                WpT = math.sqrt( ( pfMet*math.cos(pfMetPhi) + myEles[iLeadLep].Px())**2 + ( pfMet*math.sin(pfMetPhi) + myEles[iLeadLep].Py())**2)
                            
                if nBjets==1:
                    allquantities.reg_1e1b_Wmass = Wenumass    
                    allquantities.reg_1e1b_WpT=WpT
                    
                    allquantities.reg_1e1b_hadrecoil = WenuRecoil
                    allquantities.reg_1e1b_MET = pfMet
                    
                    allquantities.reg_1e1b_lep1_pT=myEles[iLeadLep].Pt()
                    
                    allquantities.reg_1e1b_jet1_pT=j1.Pt()
                    
                    if options.CSV:
                        if nTHINJets>1: allquantities.reg_1e1b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_1e1b_jet1_eta=j1.Eta()
                        if nTHINJets>1: allquantities.reg_1e1b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_1e1b_jet1_csv = thinJetCSV[ifirstjet]
                        if nTHINJets>1: allquantities.reg_1e1b_jet2_csv = thinJetCSV[isecondjet]
                        
                        allquantities.reg_1e1b_njet = nTHINJets
                    if options.DeepCSV:
                        if nTHINdeepCSVJets>1: allquantities.reg_1e1b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_1e1b_jet1_eta=j1.Eta()
                        if nTHINdeepCSVJets>1: allquantities.reg_1e1b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_1e1b_jet1_deepcsv = thinJetdeepCSV[ifirstjet]
                        if nTHINdeepCSVJets>1: allquantities.reg_1e1b_jet2_deepcsv = thinJetdeepCSV[isecondjet]
                        
                        allquantities.reg_1e1b_njet = nTHINdeepCSVJets

                    allquantities.reg_1e1b_ntau = nTau
                    allquantities.reg_1e1b_nele = nEle
                    allquantities.reg_1e1b_nmu = nMu
                    allquantities.reg_1e1b_nUncleanTau = nUncleanTau
                    allquantities.reg_1e1b_nUncleanEle = nUncleanEle
                    allquantities.reg_1e1b_nUncleanMu = nUncleanMu
                    
            #1e, 2 b-tagged   
                if nBjets==2 and SR2jet2:         
                    allquantities.reg_1e2b_Wmass = Wenumass   
                    allquantities.reg_1e2b_WpT=WpT
                    
                    allquantities.reg_1e2b_hadrecoil = WenuRecoil
                    allquantities.reg_1e2b_MET = pfMet
                    
                    allquantities.reg_1e2b_lep1_pT=myEles[iLeadLep].Pt()
                    
                    allquantities.reg_1e2b_jet1_pT=j1.Pt()
                    if options.CSV:
                        if nTHINJets>1: allquantities.reg_1e2b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_1e2b_jet1_eta=j1.Eta()
                        if nTHINJets>1: allquantities.reg_1e2b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_1e2b_jet1_csv = thinJetCSV[ifirstjet]
                        if nTHINJets>1: allquantities.reg_1e2b_jet2_csv = thinJetCSV[isecondjet]
                        
                        allquantities.reg_1e2b_njet = nTHINJets
                    if options.DeepCSV:
                        if nTHINdeepCSVJets>1: allquantities.reg_1e2b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_1e2b_jet1_eta=j1.Eta()
                        if nTHINdeepCSVJets>1: allquantities.reg_1e2b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_1e2b_jet1_deepcsv = thinJetdeepCSV[ifirstjet]
                        if nTHINdeepCSVJets>1: allquantities.reg_1e2b_jet2_deepcsv = thinJetdeepCSV[isecondjet]
                        
                        allquantities.reg_1e2b_njet = nTHINdeepCSVJets
                    allquantities.reg_1e2b_ntau = nTau
                    allquantities.reg_1e2b_nele = nEle
                    allquantities.reg_1e2b_nmu = nMu
                    allquantities.reg_1e2b_nUncleanTau = nUncleanTau
                    allquantities.reg_1e2b_nUncleanEle = nUncleanEle
                    allquantities.reg_1e2b_nUncleanMu = nUncleanMu
                    
  #Cutflow                  
        if ((UnPrescaledIsoMu20 and HLT_IsoMu20) or HLT_IsoMu24_v or HLT_IsoTkMu24_v):
            CR1mu1bCutFlow['trig']+=1
            CR1mu2bCutFlow['trig']+=1
            if nEle==0 and nMu==1:
                CR1mu1bCutFlow['nlep']+=1
                CR1mu2bCutFlow['nlep']+=1
                if myMuos[0].Pt() > 30. and myMuTightID[0]:
                    CR1mu1bCutFlow['lepconds']+=1
                    CR1mu2bCutFlow['lepconds']+=1
                    if WmunuRecoil>200.:
                        CR1mu1bCutFlow['recoil']+=1
                        CR1mu2bCutFlow['recoil']+=1
                        if Wmunumass>50. and Wmunumass<160.:
                            CR1mu1bCutFlow['mass']+=1  
                            CR1mu2bCutFlow['mass']+=1   
                            if WdPhicond:
                                CR1mu1bCutFlow['dPhicond']+=1 
                                CR1mu2bCutFlow['dPhicond']+=1                                
                                if nBjets==1:
                                    CR1mu1bCutFlow['nbjets']+=1
                                    if jetcond:
                                        CR1mu1bCutFlow['jetconds']+=1
                                if nBjets==2:
                                    CR1mu2bCutFlow['nbjets']+=1
                                    if jetcond and SR2jet2:
                                        CR1mu2bCutFlow['jetconds']+=1                  
        #1mu, 1 b-tagged  
        if nMu==1 and nEle==0 and ((UnPrescaledIsoMu20 and HLT_IsoMu20) or HLT_IsoMu24_v or HLT_IsoTkMu24_v) and WmunuRecoil>200. and jetcond and WdPhicond and Wmunumass>50. and Wmunumass<160.:
            iLeadLep=0
                
            if myMuos[iLeadLep].Pt() > 30. and myMuTightID[iLeadLep]:       # and myMuIso[iLeadLep]<0.15
                        
                WpT = math.sqrt( ( pfMet*math.cos(pfMetPhi) + myMuos[iLeadLep].Px())**2 + ( pfMet*math.sin(pfMetPhi) + myMuos[iLeadLep].Py())**2)                
                            
                if  nBjets==1:         
                    allquantities.reg_1mu1b_Wmass = Wmunumass    
                    allquantities.reg_1mu1b_WpT=WpT
                    
                    allquantities.reg_1mu1b_hadrecoil = WmunuRecoil
                    allquantities.reg_1mu1b_MET = pfMet
                    
                    allquantities.reg_1mu1b_lep1_pT=myMuos[iLeadLep].Pt()                    
                    allquantities.reg_1mu1b_lep1_iso=myMuIso[iLeadLep]
                    
                    allquantities.reg_1mu1b_jet1_pT=j1.Pt()
                    if options.CSV:
                        if nTHINJets>1: allquantities.reg_1mu1b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_1mu1b_jet1_eta=j1.Eta()
                        if nTHINJets>1: allquantities.reg_1mu1b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_1mu1b_jet1_csv = thinJetCSV[ifirstjet]
                        if nTHINJets>1: allquantities.reg_1mu1b_jet2_csv = thinJetCSV[isecondjet]
                        
                        allquantities.reg_1mu1b_njet = nTHINJets
                    if options.DeepCSV:
                        if nTHINdeepCSVJets>1: allquantities.reg_1mu1b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_1mu1b_jet1_eta=j1.Eta()
                        if nTHINdeepCSVJets>1: allquantities.reg_1mu1b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_1mu1b_jet1_deepcsv = thinJetdeepCSV[ifirstjet]
                        if nTHINdeepCSVJets>1: allquantities.reg_1mu1b_jet2_deepcsv = thinJetdeepCSV[isecondjet]
                        
                        allquantities.reg_1mu1b_njet = nTHINdeepCSVJets
                    allquantities.reg_1mu1b_ntau = nTau
                    allquantities.reg_1mu1b_nele = nEle
                    allquantities.reg_1mu1b_nmu = nMu
                    allquantities.reg_1mu1b_nUncleanTau = nUncleanTau
                    allquantities.reg_1mu1b_nUncleanEle = nUncleanEle
                    allquantities.reg_1mu1b_nUncleanMu = nUncleanMu
                    
            #1mu, 2 b-tagged        
                if  nBjets==2 and SR2jet2:
                    allquantities.reg_1mu2b_Wmass = Wmunumass   
                    allquantities.reg_1mu2b_WpT=WpT
                    
                    allquantities.reg_1mu2b_hadrecoil = WmunuRecoil
                    allquantities.reg_1mu2b_MET = pfMet
                    
                    allquantities.reg_1mu2b_lep1_pT=myMuos[iLeadLep].Pt()
                    allquantities.reg_1mu2b_lep1_iso=myMuIso[iLeadLep]
                    
                    allquantities.reg_1mu2b_jet1_pT=j1.Pt()
                    if options.CSV:
                        if nTHINJets>1: allquantities.reg_1mu2b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_1mu2b_jet1_eta=j1.Eta()
                        if nTHINJets>1: allquantities.reg_1mu2b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_1mu2b_jet1_csv = thinJetCSV[ifirstjet]
                        if nTHINJets>1: allquantities.reg_1mu2b_jet2_csv = thinJetCSV[isecondjet]
                        
                        allquantities.reg_1mu2b_njet = nTHINJets
                    if options.DeepCSV:
                        if nTHINdeepCSVJets>1: allquantities.reg_1mu2b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_1mu2b_jet1_eta=j1.Eta()
                        if nTHINdeepCSVJets>1: allquantities.reg_1mu2b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_1mu2b_jet1_deepcsv = thinJetdeepCSV[ifirstjet]
                        if nTHINJets>1: allquantities.reg_1mu2b_jet2_deepcsv = thinJetdeepCSV[isecondjet]
                        
                        allquantities.reg_1mu2b_njet = nTHINdeepCSVJets
                    allquantities.reg_1mu2b_ntau = nTau
                    allquantities.reg_1mu2b_nele = nEle
                    allquantities.reg_1mu2b_nmu = nMu
                    allquantities.reg_1mu2b_nUncleanTau = nUncleanTau
                    allquantities.reg_1mu2b_nUncleanEle = nUncleanEle
                    allquantities.reg_1mu2b_nUncleanMu = nUncleanMu
                    
# -------------------------------------------
# Top CR
# -------------------------------------------
            
        #Top CR specific bools
        
        TopdPhicond=True
        
        if TOPPhi>-10.:
            if DeltaPhi(j1.Phi(),Phi_mpi_pi(math.pi+TOPPhi)) < 0.5: TopdPhicond = False      #Added +pi to ZPhi to reverse an error in SkimTree which will be fixed in next iteration.
            if options.CSV:
                if nTHINJets>=2:
                    if DeltaPhi(j2.Phi(),Phi_mpi_pi(math.pi+TOPPhi)) < 0.5: TopdPhicond=False
                if nTHINJets>=3:
                    if DeltaPhi(j3.Phi(),Phi_mpi_pi(math.pi+TOPPhi)) < 0.5: TopdPhicond=False
            if options.DeepCSV:
                if nTHINdeepCSVJets>=2:
                    if DeltaPhi(j2.Phi(),Phi_mpi_pi(math.pi+TOPPhi)) < 0.5: TopdPhicond=False
                if nTHINdeepCSVJets>=3:
                    if DeltaPhi(j3.Phi(),Phi_mpi_pi(math.pi+TOPPhi)) < 0.5: TopdPhicond=False
                    
 #Cutflow                  
        if ((UnPrescaledIsoMu20 and HLT_IsoMu20) or HLT_IsoMu24_v or HLT_IsoTkMu24_v):
            CR1mu1e1bCutFlow['trig']+=1
            CR1mu1e2bCutFlow['trig']+=1
            if nEle==1 and nMu==1:
                CR1mu1e1bCutFlow['nlep']+=1
                CR1mu1e2bCutFlow['nlep']+=1
                if myEles[0].Pt() > 30. and myEleTightID[0] and myMuos[0].Pt() > 30. and myMuTightID[0] and myMuIso[0]<0.15:
                    CR1mu1e1bCutFlow['lepconds']+=1
                    CR1mu1e2bCutFlow['lepconds']+=1
                    if TOPRecoil>200.:
                        CR1mu1e1bCutFlow['recoil']+=1
                        CR1mu1e2bCutFlow['recoil']+=1
                        CR1mu1e1bCutFlow['mass']+=1  
                        CR1mu1e2bCutFlow['mass']+=1   
                        if TopdPhicond:
                            CR1mu1e1bCutFlow['dPhicond']+=1 
                            CR1mu1e2bCutFlow['dPhicond']+=1                                
                            if nBjets==1:
                                CR1mu1e1bCutFlow['nbjets']+=1
                                if jetcond:
                                    CR1mu1e1bCutFlow['jetconds']+=1
                            if nBjets==2:
                                CR1mu1e2bCutFlow['nbjets']+=1
                                if jetcond and SR2jet2:
                                    CR1mu1e2bCutFlow['jetconds']+=1             
                
                 
        #1mu, 1e, 1 b-tagged
        if nEle==1 and nMu==1 and ((UnPrescaledIsoMu20 and HLT_IsoMu20) or HLT_IsoMu24_v or HLT_IsoTkMu24_v) and TOPRecoil>200. and jetcond and TopdPhicond:     
        
            if myEles[0].Pt() > 30. and myEleTightID[0] and myMuos[0].Pt() > 30. and myMuTightID[0] and myMuIso[0]<0.15:
            
                if myEles[0].Pt() > myMuos[0].Pt():
                    EleLead=True
                else:
                    EleLead=False
                            
                if nBjets==1:
                    
                    allquantities.reg_1mu1e1b_hadrecoil = TOPRecoil
                    allquantities.reg_1mu1e1b_MET = pfMet
                    
                    if EleLead:
                        allquantities.reg_1mu1e1b_lep1_pT=myEles[0].Pt()
                        allquantities.reg_1mu1e1b_lep2_pT=myMuos[0].Pt()
                        allquantities.reg_1mu1e1b_lep2_iso=myMuIso[0]
                    else:
                        allquantities.reg_1mu1e1b_lep2_pT=myEles[0].Pt()
                        allquantities.reg_1mu1e1b_lep1_pT=myMuos[0].Pt()
                        allquantities.reg_1mu1e1b_lep1_iso=myMuIso[0]                    
                    
                    allquantities.reg_1mu1e1b_e_pT=myEles[0].Pt()
                    allquantities.reg_1mu1e1b_mu_pT=myMuos[0].Pt()
                    allquantities.reg_1mu1e1b_mu_iso=myMuIso[0]
                    
                    allquantities.reg_1mu1e1b_jet1_pT=j1.Pt()
                    if options.CSV:
                        if nTHINJets>1: allquantities.reg_1mu1e1b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_1mu1e1b_jet1_eta=j1.Eta()
                        if nTHINJets>1: allquantities.reg_1mu1e1b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_1mu1e1b_jet1_csv = thinJetCSV[ifirstjet]
                        if nTHINJets>1: allquantities.reg_1mu1e1b_jet2_csv = thinJetCSV[isecondjet]
                        
                        allquantities.reg_1mu1e1b_njet = nTHINJets
                    if options.DeepCSV:
                        if nTHINdeepCSVJets>1: allquantities.reg_1mu1e1b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_1mu1e1b_jet1_eta=j1.Eta()
                        if nTHINdeepCSVJets>1: allquantities.reg_1mu1e1b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_1mu1e1b_jet1_deepcsv = thinJetdeepCSV[ifirstjet]
                        if nTHINJets>1: allquantities.reg_1mu1e1b_jet2_deepcsv = thinJetdeepCSV[isecondjet]
                        
                        allquantities.reg_1mu1e1b_njet = nTHINdeepCSVJets
                    allquantities.reg_1mu1e1b_ntau = nTau
                    allquantities.reg_1mu1e1b_nele = nEle
                    allquantities.reg_1mu1e1b_nmu = nMu
                    allquantities.reg_1mu1e1b_nUncleanTau = nUncleanTau
                    allquantities.reg_1mu1e1b_nUncleanEle = nUncleanEle
                    allquantities.reg_1mu1e1b_nUncleanMu = nUncleanMu
                    
            #1mu, 1e, 2 b-tagged   
                if nBjets==2 and SR2jet2:         
                    
                    allquantities.reg_1mu1e2b_hadrecoil = TOPRecoil
                    allquantities.reg_1mu1e2b_MET = pfMet
                    
                    if EleLead:
                        allquantities.reg_1mu1e2b_lep1_pT=myEles[0].Pt()
                        allquantities.reg_1mu1e2b_lep2_pT=myMuos[0].Pt()
                        allquantities.reg_1mu1e2b_lep2_iso=myMuIso[0]
                    else:
                        allquantities.reg_1mu1e2b_lep2_pT=myEles[0].Pt()
                        allquantities.reg_1mu1e2b_lep1_pT=myMuos[0].Pt()
                        allquantities.reg_1mu1e2b_lep1_iso=myMuIso[0]                    
                    
                    allquantities.reg_1mu1e2b_e_pT=myEles[0].Pt()
                    allquantities.reg_1mu1e2b_mu_pT=myMuos[0].Pt()
                    allquantities.reg_1mu1e2b_mu_iso=myMuIso[0]
                    
                    allquantities.reg_1mu1e2b_jet1_pT=j1.Pt()
                    if options.CSV:
                        if nTHINJets>1: allquantities.reg_1mu1e2b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_1mu1e2b_jet1_eta=j1.Eta()
                        if nTHINJets>1: allquantities.reg_1mu1e2b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_1mu1e2b_jet1_csv = thinJetCSV[ifirstjet]
                        if nTHINJets>1: allquantities.reg_1mu1e2b_jet2_csv = thinJetCSV[isecondjet]
                        
                        allquantities.reg_1mu1e2b_njet = nTHINJets
                    if options.DeepCSV:
                        if nTHINdeepCSVJets>1: allquantities.reg_1mu1e2b_jet2_pT=j2.Pt()
                        
                        allquantities.reg_1mu1e2b_jet1_eta=j1.Eta()
                        if nTHINdeepCSVJets>1: allquantities.reg_1mu1e2b_jet2_eta=j2.Eta()
                    
                        allquantities.reg_1mu1e2b_jet1_deepcsv = thinJetdeepCSV[ifirstjet]
                        if nTHINJets>1: allquantities.reg_1mu1e2b_jet2_deepcsv = thinJetdeepCSV[isecondjet]
                        
                        allquantities.reg_1mu1e2b_njet = nTHINdeepCSVJets
                    allquantities.reg_1mu1e2b_ntau = nTau
                    allquantities.reg_1mu1e2b_nele = nEle
                    allquantities.reg_1mu1e2b_nmu = nMu
                    allquantities.reg_1mu1e2b_nUncleanTau = nUncleanTau
                    allquantities.reg_1mu1e2b_nUncleanEle = nUncleanEle
                    allquantities.reg_1mu1e2b_nUncleanMu = nUncleanMu
                
                
                
        
#        allquantities.FillRegionHisto()
        
        #Control Regions
        
#        if ZeeRecoil!=-1. or ZmumuRecoil!=-1. or WenuRecoil!=-1. or  WmunuRecoil!=-1. or  TOPRecoil!=-1. :
#            print nEle
#            print nMu
#            print ZeeRecoil
#            print ZmumuRecoil
#            print WenuRecoil
#            print WmunuRecoil
#            print TOPRecoil
#            print
        
        
        
#        if ZeeRecoil>200.:
#        
#            zeeCR=True
#            CRCutFlow['zrecoil']+=1
#            
#            if nEle!=2 or nMu!=0 or nTau!=0:zeeCR=False
#            if zeeCR: CRCutFlow['zLep']+=1
#            
#            if SR1jetcond==False and SR2jetcond==False: zeeCR=False
#            if zeeCR: CRCutFlow['zJet']+=1            
#            
#            if ZeeMass < 70. or ZeeMass > 110: zeeCR=False
#            if zeeCR: CRCutFlow['zMass']+=1
#            
#            
#            
#            
#        if ZmumuRecoil>200.:
#        
#            zmumuCR=True
#            CRCutFlow['zrecoil']+=1            
#            
#            if nMu!=2 or nEle!=0 or nTau!=0: zmumuCR=False
#            if zmumuCR: CRCutFlow['zLep']+=1
#            
#            if SR1jetcond==False and SR2jetcond==False: zmumuCR=False
#            if zmumuCR: CRCutFlow['zJet']+=1
#        
#            if ZmumuMass < 70. or ZmumuMass > 110: zmumuCR=False
#            if zmumuCR: CRCutFlow['zMass']+=1
#        
#        
#        
#        if SR1jetcond==False and SR2jetcond==False:   continue
        
        
        
        #=================================================================
#        #  Z control region
#        #=================================================================
#        zCR=True
#        
#        zCREle=False
#        zCRMu=False
#        
#        if nEle==2 and nMu==0 and nTau==0:
#            zCREle=True
#            LepP4=myEles
#            isLoose=myEleLooseID
#            isTight=myEleTightID
#            zmass=ZeeMass
#            hadrecoil=ZeeRecoil
#            CRCutFlow['nlepcond']+=1
#        elif nMu==2 and nEle==0 and nTau==0:
#            zCRMu=True
#            LepP4=myMuos
#            isLoose=myMuLooseID
#            isTight=myMuTightID
#            zmass=ZmumuMass
#            hadrecoil=ZmumuRecoil
#            CRCutFlow['nlepcond']+=1
#        else:
#            zCR=False
#        
#        if zCR:                                                 # Just to reduce reduntant computation
#            if LepP4[0].Pt() > LepP4[1].Pt():
#                iLeadLep=0
#                iSecondLep=1
#            else:
#                iLeadLep=1
#                iSecondLep=0
#            
#            # Leading lepton conditions:
#            if LepP4[iLeadLep].Pt() < 30.: zCR=False
##            print "isLoose: "+str(isLoose[iLeadLep])            #To see the data type in case of bugs (if any)
#            if not isTight[iLeadLep]: zCR=False
#            
#            # Sub-leading lepton conditions:
#            if LepP4[iSecondLep].Pt() < 10.: zCR=False
#            if not isLoose[iSecondLep]: zCR=False
#            
#           
#            Zmu1pT = None
#            Zmu1eta = None
#            Zmu1phi = None
#            Zele1pT = None
#            Zele1eta = None
#            Zele1phi = None
#            Zmu2pT = None
#            Zmu2eta = None
#            Zmu2phi = None
#            Zele2pT = None
#            Zele2eta = None
#            Zele2phi = None
#            Zmu1Iso  = None
#            Zmu2Iso  = None
#            if zCRMu:                                           # Special isolation requirement for Muon                
#                if myMuIso[iLeadLep] > 0.15: zCR=False             
#                if myMuIso[iSecondLep] > 0.25: zCR=False
#                Zmu1Iso = myMuIso[iLeadLep]
#                Zmu2Iso = myMuIso[iSecondLep]
#                Zmu1pT  = LepP4[iLeadLep].Pt()
#                Zmu1eta = LepP4[iLeadLep].Eta()
#                Zmu1phi = LepP4[iLeadLep].Phi()
#                Zmu2pT = LepP4[iSecondLep].Pt()
#                Zmu2eta = LepP4[iSecondLep].Eta()
#                Zmu2phi = LepP4[iSecondLep].Phi()                
#                
#            if zCREle:
#                Zele1pT  = LepP4[iLeadLep].Pt()
#                Zele1eta = LepP4[iLeadLep].Eta()
#                Zele1phi = LepP4[iLeadLep].Phi()
#                Zele2pT =  LepP4[iSecondLep].Pt()
#                Zele2eta = LepP4[iSecondLep].Eta()
#                Zele2phi = LepP4[iSecondLep].Phi()
#                
#            ZpT = math.sqrt( (LepP4[iLeadLep].Px()+LepP4[iSecondLep].Px())*(LepP4[iLeadLep].Px()+LepP4[iSecondLep].Px()) + (LepP4[iLeadLep].Py()+LepP4[iSecondLep].Py())*(LepP4[iLeadLep].Py()+LepP4[iSecondLep].Py()) )        #Identical for Mu and Ele
#            
#            # Z Mass condition:
#            if zmass <= 70. or zmass >= 110.: zCR=False             
#            
#            # Hadronic recoil:
#            if hadrecoil <= 200.: zCR=False 
#            
#        if zCR:
#            if SR1jetcond:
#                CRStatus['ZCRSR1']+=1
#            if SR2jetcond:
#                CRStatus['ZCRSR2']+=1
#                
#                
#        #=================================================================
#        #  W control region
#        #=================================================================       
#        wCR=True
#        
#        wCREle=False
#        wCRMu=False
#        
#        if nEle==1 and nMu==0 and nTau==0:
#            wCREle=True
#            LepP4=myEles
#            isTight=myEleTightID
#            wmass=Wenumass
#            hadrecoil=WenuRecoil
#        elif nMu==1 and nEle==0 and nTau==0:
#            wCRMu=True
#            LepP4=myMuos
#            isTight=myMuTightID
#            wmass=Wmunumass
#            hadrecoil=WmunuRecoil
#        else:
#            wCR=False
#            
#        if wCR:        
#            # Leading lepton conditions:
#            if LepP4[0].Pt() < 30.: wCR=False
#            if not isTight[0]: wCR=False
#            Wmu1pT = None
#            Wmu1eta = None
#            Wmu1phi = None
#            Wele1pT = None
#            Wele1eta = None
#            Wele1phi = None
#            Wmu1Iso  = None
#            if wCRMu:
#                if myMuIso[0] > 0.15: wCR=False
#                Wmu1Iso = myMuIso[0]
#                Wmu1pT  = LepP4[0].Pt()
#                Wmu1eta = LepP4[0].Eta()
#                Wmu1phi = LepP4[0].Phi()                
#            
#            if wCREle:
#                Wele1pT  = LepP4[0].Pt()
#                Wele1eta = LepP4[0].Eta()
#                Wele1phi = LepP4[0].Phi()
#                
#            WpT = math.sqrt( ( pfMet*math.cos(pfMetPhi) + LepP4[0].Px())*( pfMet*math.cos(pfMetPhi) + LepP4[0].Px()) + ( pfMet*math.sin(pfMetPhi) + LepP4[0].Py())*( pfMet*math.sin(pfMetPhi) + LepP4[0].Py()) )        #Identical for Mu and Ele
#            # W Mass condition:
#            if wmass <= 50. or wmass >= 160.: wCR=False    
#            
#            # Hadronic recoil:
#            if hadrecoil <= 200.: wCR=False 
#            
#        if wCR:
#            if SR1jetcond:
#                CRStatus['WCRSR1']+=1
#            if SR2jetcond:
#                CRStatus['WCRSR2']+=1
#                
#                
#        #=================================================================
#        #  Top control region
#        #================================================================= 
#        TopCR=False
#        
#        if nEle==1 and nMu==1 and nTau==0:
#            TopCR=True
#            
#            # Muon
#            if myMuos[0].Pt() < 30.: TopCR=False
#            if myMuIso[0] > 0.15: TopCR=False                  
#            if not myMuTightID[0]: TopCR=False                   
#            
#            # Electron
#            if myEles[0].Pt() < 30.: TopCR=False
#            if not myEleTightID[0]: TopCR=False
#            
#            # Hadronic recoil:
#            if TOPRecoil <= 200.: TopCR=False
#            
#            TOPmu1Iso = myMuIso[0]
#            
#            TOPmu1pT   = myMuos[0].Pt()
#            TOPmu1eta  = myMuos[0].Eta()
#            TOPmu1phi  = myMuos[0].Phi()
#            TOPele1pT  = myEles[0].Pt()
#            TOPele1eta = myEles[0].Eta()
#            TOPele1phi = myEles[0].Phi()
#            
#        if TopCR:
#            if SR1jetcond:
#                CRStatus['TopCRSR1']+=1
#            if SR2jetcond:
#                CRStatus['TopCRSR2']+=1
        
        zCR=False
        wCR=False
        TopCR=False
        
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## min DPhi
        ## nT<hinJets
        ## b-jet Veto
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Lepton Veto
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        nleptons_ = (len(myTaus) + len(myMuos) + len(myEles))
        
        #if not (nleptons_ >= nlepton) : continue 
        #if not (nleptons_ < nLepton) : continue  
        
        #regime = False
        #if isboosted: regime = True
        #if isresolved: regime = False
        
        
        #if regime: 
         #   mt_ = MT(fatjetP4[HIndex].Pt(), pfMet, Phi_mpi_pi(pfMetPhi-fatjetP4[HIndex].Phi()) )
        #if not regime: 
         #   mt_ = MT(HiggsInfo[0][3], pfMet, Phi_mpi_pi(pfMetPhi-HiggsInfo[0][4]) )
        
        #if mt_ < 450.0: continue 
        
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Photon Veto
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----to be added in future---------------------------------------------------------------------------------------------------------------------------------------
        
        
        if pfmetstatus and SRlepcond and (SR1jetcond or SR2jetcond): npass = npass + 1
        
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
#       ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        genpTReweighting = 1.0 
        if isData==1:   genpTReweighting  =  1.0
        if not isData :  genpTReweighting = GenWeightProducer(samplename, nGenPar, genParId, genMomParId, genParSt,genParP4)
#        print genpTReweighting
        
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Ele reweight
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        eleweight = 1.0
        for iele in range(nEle):
            elept = eleP4[iele].Pt()
            eleeta = eleP4[iele].Eta()
            xbin = eleEtaPt.GetXaxis().FindBin(eleeta)
            ybin = eleEtaPt.GetYaxis().FindBin(elept)
            eleweight *= eleEtaPt.GetBinContent(xbin,ybin)

        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Pileup weight
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
#        allpuweights = PUWeight()
#        len_puweight = len(allpuweights)
        puweight = 0.0
        if isData: puweight = 1.0 
        if not isData:
#            if pu_nTrueInt  <= len_puweight: puweight = allpuweights[pu_nTrueInt-1]
#            if pu_nTrueInt  > len_puweight : puweight = 0.0 
            if pu_nTrueInt < 100:
                puweight = pileup2016histo.GetBinContent(pu_nTrueInt)
            else:
                puweight = 1.
#            print pu_nTrueInt
            
        #print (len_puweight, pu_nTrueInt, puweight)
        
        
        allweights = puweight * mcweight * genpTReweighting * eleweight
#        print puweight
#        print mcWeight
#        print mcweight
#        print genpTReweighting
#        print allweights
#        print
                
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## BTag Scale Factor 
        if options.CSV:
            if SR1njetcond:
                ij = ifirstjet
                if nTHINJets>1: jj = isecondjet
                
                flav1 = jetflav(THINjetHadronFlavor[ij])
                if nTHINJets>1: flav2 = jetflav(THINjetHadronFlavor[jj])

    #            print ("ij, flav, pt, eta, ",ij, flav1, thinjetP4[ij].Pt(), thinjetP4[ij].Eta())
                reader1.eval_auto_bounds('central', 0, 1.2, 50.)
                sf_resolved1 = weightbtag(reader1, flav1, thinjetP4[ij].Pt(), thinjetP4[ij].Eta())
                if nTHINJets>1: sf_resolved2 = weightbtag(reader1, flav2, thinjetP4[jj].Pt(), thinjetP4[jj].Eta())
                
    #            print (sf_resolved1, sf_resolved2)
            elif SR2njetcond:
                ij = ifirstjet
                jj = isecondjet
                if nTHINJets>2: jk = ithirdjet
             
                flav1 = jetflav(THINjetHadronFlavor[ij])
                flav2 = jetflav(THINjetHadronFlavor[jj])
                if nTHINJets>2: flav3 = jetflav(THINjetHadronFlavor[jj])

    #            print ("ij, flav, pt, eta, ",ij, flav1, thinjetP4[ij].Pt(), thinjetP4[ij].Eta())
                reader1.eval_auto_bounds('central', 0, 1.2, 50.)
                sf_resolved1 = weightbtag(reader1, flav1, thinjetP4[ij].Pt(), thinjetP4[ij].Eta())
                sf_resolved2 = weightbtag(reader1, flav2, thinjetP4[jj].Pt(), thinjetP4[jj].Eta())
                if nTHINJets>2: sf_resolved3 = weightbtag(reader1, flav3, thinjetP4[jk].Pt(), thinjetP4[jk].Eta())
                
    #            print (sf_resolved1, sf_resolved2, sf_resolved3)
                
            if SR1njetcond:
                allweights = allweights * sf_resolved1[0]
                if nTHINJets>1: 
                    allweights = allweights * sf_resolved2[0]
            if SR2njetcond:
                allweights = allweights * sf_resolved1[0] * sf_resolved2[0]
                if nTHINJets>2: 
                    allweights = allweights * sf_resolved3[0]

        if options.DeepCSV:
            if SR1njetcond:
                ij = ifirstjet
                if nTHINdeepCSVJets>1: jj = isecondjet
                
                flav1 = jetflav(THINdeepCSVjetHadronFlavor[ij])
                if nTHINdeepCSVJets>1: flav2 = jetflav(THINdeepCSVjetHadronFlavor[jj])

    #            print ("ij, flav, pt, eta, ",ij, flav1, thinjetP4[ij].Pt(), thinjetP4[ij].Eta())
                reader1.eval_auto_bounds('central', 0, 1.2, 50.)
                sf_resolved1 = weightbtag(reader1, flav1, thindeepCSVjetP4[ij].Pt(), thindeepCSVjetP4[ij].Eta())
                if nTHINdeepCSVJets>1: sf_resolved2 = weightbtag(reader1, flav2, thindeepCSVjetP4[jj].Pt(), thindeepCSVjetP4[jj].Eta())
                
    #            print (sf_resolved1, sf_resolved2)
            elif SR2njetcond:
                ij = ifirstjet
                jj = isecondjet
                if nTHINdeepCSVJets>2: jk = ithirdjet
             
                flav1 = jetflav(THINdeepCSVjetHadronFlavor[ij])
                flav2 = jetflav(THINdeepCSVjetHadronFlavor[jj])
                if nTHINdeepCSVJets>2: flav3 = jetflav(THINdeepCSVjetHadronFlavor[jj])

    #            print ("ij, flav, pt, eta, ",ij, flav1, thinjetP4[ij].Pt(), thinjetP4[ij].Eta())
                reader1.eval_auto_bounds('central', 0, 1.2, 50.)
                sf_resolved1 = weightbtag(reader1, flav1, thindeepCSVjetP4[ij].Pt(), thindeepCSVjetP4[ij].Eta())
                sf_resolved2 = weightbtag(reader1, flav2, thindeepCSVjetP4[jj].Pt(), thindeepCSVjetP4[jj].Eta())
                if nTHINdeepCSVJets>2: sf_resolved3 = weightbtag(reader1, flav3, thindeepCSVjetP4[jk].Pt(), thindeepCSVjetP4[jk].Eta())
                
    #            print (sf_resolved1, sf_resolved2, sf_resolved3)
                
            if SR1njetcond:
                allweights = allweights * sf_resolved1[0]
                if nTHINdeepCSVJets>1: 
                    allweights = allweights * sf_resolved2[0]
            if SR2njetcond:
                allweights = allweights * sf_resolved1[0] * sf_resolved2[0]
                if nTHINdeepCSVJets>2: 
                    allweights = allweights * sf_resolved3[0]

        if isData: allweights = 1.0 
        
        
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        
        
        allquantities.met        = pfMet

        allquantities.N_e             = len(myEles)
        allquantities.N_mu            = len(myMuos)
        allquantities.N_tau           = len(myTaus)
        allquantities.N_Pho           = 0
        allquantities.N_b             = len(mybjets)
        allquantities.N_j             = nTHINJets

        allquantities.weight    = allweights
        allquantities.totalevents = 1
        
        allquantlist=AllQuantList.getAll()
        
        for quant in allquantlist:
            exec("allquantities."+quant+" = None")                              # Presets all quantities to None  
                 
        if SR1jetcond and pfmetstatus and SRlepcond and keepevent and writeSR1:            
           allquantities.jet1_pT_sr1     = jetSR1Info[0][0]
           allquantities.jet1_eta_sr1    = jetSR1Info[0][1]
           allquantities.jet1_phi_sr1    = jetSR1Info[0][2]
           if options.CSV:
               allquantities.jet1_csv_sr1    = jetSR1Info[0][3]
           if options.DeepCSV:
               allquantities.jet1_deepcsv_sr1    = jetSR1Info[0][3]
           allquantities.jet2_pT_sr1     = jetSR1Info[1][0]
           allquantities.jet2_eta_sr1    = jetSR1Info[1][1]
           allquantities.jet2_phi_sr1    = jetSR1Info[1][2]
           if options.CSV:
               allquantities.jet2_csv_sr1    = jetSR1Info[1][3]
           if options.DeepCSV:
               allquantities.jet2_deepcsv_sr1    = jetSR1Info[1][3]
           allquantities.min_dPhi_sr1    = jetSR1Info[2]
           allquantities.met_sr1         = jetSR1Info[3]
           allquantities.jet1_nhf_sr1    = jetSR1Info[4]
           allquantities.jet1_chf_sr1    = jetSR1Info[5]
           
           
        elif SR2jetcond and pfmetstatus and SRlepcond and keepevent and writeSR2:
           allquantities.jet1_pT_sr2     = jetSR2Info[0][0]
           allquantities.jet1_eta_sr2    = jetSR2Info[0][1]
           allquantities.jet1_phi_sr2    = jetSR2Info[0][2]
           if options.CSV:
               allquantities.jet1_csv_sr2    = jetSR2Info[0][3]
           if options.DeepCSV:
               allquantities.jet1_deepcsv_sr2    = jetSR2Info[0][3]
           
           allquantities.jet2_pT_sr2     = jetSR2Info[1][0]
           allquantities.jet2_eta_sr2    = jetSR2Info[1][1]
           allquantities.jet2_phi_sr2    = jetSR2Info[1][2]
           if options.CSV:
               allquantities.jet2_csv_sr2    = jetSR2Info[1][3]
           if options.DeepCSV:
               allquantities.jet2_deepcsv_sr2    = jetSR2Info[1][3]
           
           allquantities.jet3_pT_sr2     = jetSR2Info[2][0]
           allquantities.jet3_eta_sr2    = jetSR2Info[2][1]
           allquantities.jet3_phi_sr2    = jetSR2Info[2][2]
           if options.CSV:
               allquantities.jet3_csv_sr2    = jetSR2Info[2][3]
           if options.DeepCSV:
               allquantities.jet3_deepcsv_sr2    = jetSR2Info[2][3]
           
           allquantities.min_dPhi_sr2    = jetSR2Info[3]      
           allquantities.met_sr2         = jetSR2Info[4]   
           allquantities.jet1_nhf_sr2    = jetSR2Info[5]
           allquantities.jet1_chf_sr2    = jetSR2Info[6]  
                       

            
        #print (allquantities.regime, allquantities.met,allquantities.mass )
        allquantities.FillRegionHisto()
        allquantities.FillHisto()
        
        
        
    

    #print cutStatus    
    NEntries_Weight = h_t_weight.Integral()
    NEntries_total  = h_t.Integral()
    cutStatus['total'] = int(NEntries_total)
    cutStatusSR1['total'] = int(NEntries_total)
    cutStatusSR2['total'] = int(NEntries_total)
    
    for CRreg in regionnames:
        exec("CR"+CRreg+"CutFlow['total']=int(NEntries_total)")
    
    
    cutStatusSR1['pfmet'] = cutStatus['pfmet']
    cutStatusSR2['pfmet'] = cutStatus['pfmet']
    print "Total events =", int(NEntries_total)
    print "Preselected events=", cutStatus['preselection']
    print "Selected events =", npass
    
    # Cutflow
    cutflowTable=""
    cutflowHeader=""
    cutflowvalues=[]    
    cutflownames=['total','preselection','pfmet','njet+nBjet','jet1','jet2/3','lep']
    for cutflowname in cutflownames:   
        cutflowvalues.append(cutStatus[cutflowname])
        cutflowTable += str(cutStatus[cutflowname])+" "
        cutflowHeader += cutflowname+" "    
    
    cutflowvaluesSR1=[]
    cutflownamesSR1=['total','preselection','pfmet','njet+nBjet','jet1','jet2','lep']
    for cutflowname in cutflownamesSR1:   
        cutflowvaluesSR1.append(cutStatusSR1[cutflowname])
        
    cutflowvaluesSR2=[]
    cutflownamesSR2=['total','preselection','pfmet','njet+nBjet','jet1','jet2','jet3','lep']
    for cutflowname in cutflownamesSR2:   
        cutflowvaluesSR2.append(cutStatusSR2[cutflowname])
    
    # CR counts
    CRTable=""
    CRHeader=""
    CRvalues=[]
    for CRname in CRs:   
        CRvalues.append(CRStatus[CRname])
        CRTable += str(CRStatus[CRname])+" "
        CRHeader += CRname+" "
    
    #Cutflows for SR1 and SR2
    print "\nCutflow:"
    print cutStatus
    print
    print "SR1:"
    print cutStatusSR1
    print
    print "SR2:"
    print cutStatusSR2
    print
    print "CRs:"
    print CRStatus
#    print
#    print "CR Cutflow:"
#    print CRCutFlow
    print
    
    CRcutflowvaluesSet=[]
    CRcutnames=['total','preselection']+CRcutnames
    for CRreg in regionnames:
        CFvalues=[]
        for cutname in CRcutnames:
            exec("CFvalues.append(CR"+CRreg+"CutFlow['"+cutname+"'])") 
        CRcutflowvaluesSet.append(CFvalues)    
    
    allquantities.WriteHisto((NEntries_total,NEntries_Weight,npass,cutflowvalues,cutflownames,cutflowvaluesSR1,cutflownamesSR1,cutflowvaluesSR2,cutflownamesSR2,CRvalues,CRs,regionnames,CRcutnames,CRcutflowvaluesSet))
    
    if NEntries > 0:
        eff=round(float(npass/float(NEntries_total)),5)
    else:
        eff = "NA"
    print "efficiency =", eff  
    
    os.system("mkdir -p "+outputdir+'/efficiencyfiles/')
        
    f = open(outputdir+'/efficiencyfiles/'+textfile, 'w')
    f.write(str(eff)+"\n\n#Cutflow Table:\n"+cutflowHeader[:-1]+"\n"+cutflowTable[:-1]+"\n\n#CR Table:\n"+CRHeader[:-1]+"\n"+CRTable[:-1])
    print "ROOT file written to", outfilename
    print "Log written to "+outputdir+'/efficiencyfiles/'+textfile
    f.close()
    print "Completed."



        
def CheckFilter(filterName, filterResult,filtercompare):
    ifilter_=0
    filter1 = False
    for ifilter in filterName:
        filter1 = (ifilter.find(filtercompare) != -1)  & (bool(filterResult[ifilter_]) == True)   
        if filter1: break
        ifilter_ = ifilter_ + 1
    return filter1





######################################
######################################
######################################
def MakeTable():
    print "called MakeTable"
    files= [inputfilename]
    legend=legendTemplate
    prefix="V_met_"
    effnamelist = [prefix + ihisto  for ihisto in namelist]
    inputfile={}
    histList=[]
    for ifile_ in range(len(files)):
        print ("opening file  "+files[ifile_])
        inputfile[ifile_] = TFile( files[ifile_] )
        print "fetching histograms"
        for ihisto_ in range(len(effnamelist)):
            histo = inputfile[ifile_].Get(effnamelist[ihisto_])
            histList.append(histo)
    
    for ih in range(len(histList)):
        eff = ("%0.4f" % float(histList[ih].Integral()/histList[0].Integral()))
        toprint =  legendTemplate[ih] + " & " + str(eff) + " \\\\"
        print toprint


def DeltaR(p4_1, p4_2):
    eta1 = p4_1.Eta()
    eta2 = p4_2.Eta()
    eta = eta1 - eta2
    eta_2 = eta * eta
    
    phi1 = p4_1.Phi()
    phi2 = p4_2.Phi()
    phi = Phi_mpi_pi(phi1-phi2)
    phi_2 = phi * phi

    return math.sqrt(eta_2 + phi_2)
    
def DeltaPhi(phi1,phi2):
   phi = Phi_mpi_pi(phi1-phi2)
   
   return abs(phi)
   
    
def Phi_mpi_pi(x):
    kPI = 3.14159265358979323846
    kTWOPI = 2 * kPI
    
    while (x >= kPI): x = x - kTWOPI;
    while (x < -kPI): x = x + kTWOPI;
    return x;
    
def weightbtag(reader, flav, pt, eta):
    sf_c = reader.eval_auto_bounds('central', flav, eta, pt) 
    sf_low = reader.eval_auto_bounds('down', flav, eta, pt)
    sf_up  = reader.eval_auto_bounds('up', flav, eta, pt)
    btagsf = [sf_c, sf_low, sf_up]
    return btagsf

def jetflav(flav):
    if flav == 5: 
        flavor = 0
    elif flav == 4:
        flavor = 1
    else:
        flavor = 2
    return flavor




def GenWeightProducer(sample,nGenPar, genParId, genMomParId, genParSt,genParP4):
    pt__=0;
    #print " inside gen weight "
    k2=1.0
    #################
    # WJets
    #################
    if sample=="WJETS":
        goodLepID = []
        for ig in range(nGenPar): 
            PID    = genParId[ig]
            momPID = genMomParId[ig]
            status = genParSt[ig]
            if ( (abs(PID) != 11) & (abs(PID) != 12) &  (abs(PID) != 13) & (abs(PID) != 14) &  (abs(PID) != 15) &  (abs(PID) != 16) ): continue
            if ( ( (status != 1) & (abs(PID) != 15)) | ( (status != 2) & (abs(PID) == 15)) ): continue
            if ( (abs(momPID) != 24) & (momPID != PID) ): continue
            goodLepID.append(ig)

        if len(goodLepID) == 2 :
            l4_thisLep = genParP4[goodLepID[0]]
            l4_thatLep = genParP4[goodLepID[1]]
            l4_z = l4_thisLep + l4_thatLep
            
            pt = l4_z.Pt()
            pt__ = pt

            k2 = -0.830041 + 7.93714 *TMath.Power( pt - (-877.978) ,(-0.213831) ) ;
    
    #################        
    #ZJets
    #################
    if sample == "ZJETS":
        goodLepID = []
        for ig in range(nGenPar):
            PID    = genParId[ig]
            momPID = genMomParId[ig]
            status = genParSt[ig]


            if ( (abs(PID) != 12) &  (abs(PID) != 14) &  (abs(PID) != 16) ) : continue
            if ( status != 1 ) : continue 
            if ( (momPID != 23) & (momPID != PID) ) : continue
            goodLepID.append(ig)
        
        if len(goodLepID) == 2 :
            l4_thisLep = genParP4[goodLepID[0]]
            l4_thatLep = genParP4[goodLepID[1]]
            l4_z = l4_thisLep + l4_thatLep
            pt = l4_z.Pt()
            k2 = -0.180805 + 6.04146 *TMath.Power( pt - (-759.098) ,(-0.242556) ) ;

    #################        
    #TTBar
    #################        
    if (sample=="TT"):
        goodLepID = []
        for ig in range(nGenPar):
            PID    = genParId[ig]
            momPID = genMomParId[ig]
            status = genParSt[ig]
            if ( abs(PID) == 6) :
                goodLepID.append(ig)
        if(len(goodLepID)==2):
            l4_thisLep = genParP4[goodLepID[0]]
            l4_thatLep = genParP4[goodLepID[1]]
            pt1 = TMath.Min(400.0, l4_thisLep.Pt())
            pt2 = TMath.Min(400.0, l4_thatLep.Pt())
            
            w1 = TMath.Exp(0.156 - 0.00137*pt1);
            w2 = TMath.Exp(0.156 - 0.00137*pt2);
            k2 =  1.001*TMath.Sqrt(w1*w2);
            
    if(sample=="all"):
        k2 = 1.0
        
    return k2


def MT(Pt, met, dphi):
    return ROOT.TMath.Sqrt( 2 * Pt * met * (1.0 - ROOT.TMath.Cos(dphi)) )

if __name__ == "__main__":
    ## analyze the tree and make histograms and all the 2D plots and Efficiency plots. 
    if options.analyze:
        print "now calling analyzedataset"
        AnalyzeDataSet()
    
    
    

