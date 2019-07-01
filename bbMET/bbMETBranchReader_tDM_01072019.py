#!/usr/bin/env python
from ROOT import TFile, TTree, TH1F, TH1D, TH1, TCanvas, TChain,TGraphAsymmErrors, TMath, TH2D, TLorentzVector, TF1, AddressOf
import ROOT as ROOT
import os
import random
import sys, optparse
from array import array
import math
import AllQuantList, kfactor

ROOT.gROOT.SetBatch(True)
from bbMETQuantities import *
#pileup2016file = TFile('pileUPinfo2016.root')
pileup2016file = TFile('PU_Reweight_2016.root')
#pileup2016histo=pileup2016file.Get('hpileUPhist')
pileup2016histo=pileup2016file.Get('pileup')

#Electron Trigger reweights
eleTrigReweightFile = TFile('scalefactors/electron_Trigger_eleTrig.root')
eleTrig_hEffEtaPt = eleTrigReweightFile.Get('hEffEtaPt')

#Electron Reconstruction efficiency. Scale factors for 80X
eleRecoSFsFile = TFile('scalefactors/electron_Reco_SFs_egammaEffi_txt_EGM2D.root')
eleRecoSF_EGamma_SF2D = eleRecoSFsFile.Get('EGamma_SF2D')

#Loose electron ID SFs
eleLooseIDSFsFile = TFile('scalefactors/electron_Loose_ID_SFs_egammaEffi_txt_EGM2D.root')
eleLooseIDSF_EGamma_SF2D = eleLooseIDSFsFile.Get('EGamma_SF2D')

#Tight electron ID SFs
phoTightIDSFsFile = TFile('scalefactors/photon_Tight_ID_SFs_egammaEffi_txt_EGM2D.root')
phoTightIDSF_EGamma_SF2D = phoTightIDSFsFile.Get('EGamma_SF2D')

#Loose photon ID SFs
phoLooseIDSFsFile = TFile('scalefactors/photon_Loose_ID_SFs_egammaEffi_txt_EGM2D.root')
phoLooseIDSF_EGamma_SF2D = phoLooseIDSFsFile.Get('EGamma_SF2D')

#Tight photon ID SFs
eleTightIDSFsFile = TFile('scalefactors/electron_Tight_ID_SFs_egammaEffi_txt_EGM2D.root')
eleTightIDSF_EGamma_SF2D = eleTightIDSFsFile.Get('EGamma_SF2D')

# Veto cut-based electron ID SFs
eleVetoCutBasedIDSFsFile = TFile('scalefactors/electron_Veto_cut-based_ID_SFs_egammaEffi_txt_EGM2D.root')
eleVetoCutBasedIDSF_egammaEffi_txt_EGM2D = eleVetoCutBasedIDSFsFile.Get('EGamma_SF2D')

#Muon Trigger SFs
#BCDEF
muonTrigSFsRunBCDEFFile = TFile('scalefactors/muon_single_lepton_trigger_EfficienciesAndSF_RunBtoF.root')
muonTrigSFs_EfficienciesAndSF_RunBtoF = muonTrigSFsRunBCDEFFile.Get('IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio')
#GH
muonTrigSFsRunGHFile = TFile('scalefactors/muon_single_lepton_trigger_EfficienciesAndSF_Period4.root')
muonTrigSFs_EfficienciesAndSF_Period4 = muonTrigSFsRunBCDEFFile.Get('IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio')

#Muon ID SFs
#BCDEF
muonIDSFsBCDEFFile = TFile('scalefactors/muon_ID_SFs_EfficienciesAndSF_BCDEF.root')
muonLooseIDSFs_EfficienciesAndSF_BCDEF = muonIDSFsBCDEFFile.Get('MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio')
muonTightIDSFs_EfficienciesAndSF_BCDEF = muonIDSFsBCDEFFile.Get('MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio')
#GH
muonIDSFsGHFile = TFile('scalefactors/muon_ID_SFs_EfficienciesAndSF_GH.root')
muonLooseIDSFs_EfficienciesAndSF_GH = muonIDSFsGHFile.Get('MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio')
muonTightIDSFs_EfficienciesAndSF_GH = muonIDSFsGHFile.Get('MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio')

#Muon Iso SFs
#BCDEF
muonIsoSFsBCDEFFile = TFile('scalefactors/muon_Iso_SFs_EfficienciesAndSF_BCDEF.root')
muonLooseIsoSFs_EfficienciesAndSF_BCDEF = muonIsoSFsBCDEFFile.Get('LooseISO_LooseID_pt_eta/abseta_pt_ratio')
muonTightIsoSFs_EfficienciesAndSF_BCDEF = muonIsoSFsBCDEFFile.Get('TightISO_TightID_pt_eta/abseta_pt_ratio')
#GH
muonIsoSFsGHFile = TFile('scalefactors/muon_Iso_SFs_EfficienciesAndSF_GH.root')
muonLooseIsoSFs_EfficienciesAndSF_GH = muonIsoSFsGHFile.Get('LooseISO_LooseID_pt_eta/abseta_pt_ratio')
muonTightIsoSFs_EfficienciesAndSF_GH = muonIsoSFsGHFile.Get('TightISO_TightID_pt_eta/abseta_pt_ratio')

#Muon Tracking SFs
muonTrackingSFsFile = TFile('scalefactors/muon_Tracking_SFs_Tracking_EfficienciesAndSF_BCDEFGH.root')
muonTrackingSFs_EfficienciesAndSF_BCDEFGH = muonTrackingSFsFile.Get('ratio_eff_aeta_dr030e030_corr')


#MET Trigger reweights
metTrigEff_zmmfile = TFile('scalefactors/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root')
metTrig_firstmethod = metTrigEff_zmmfile.Get('hden_monojet_recoil_clone_passed')

metTrigEff_secondfile = TFile('scalefactors/metTriggerEfficiency_recoil_monojet_TH1F.root')
metTrig_secondmethod = metTrigEff_secondfile.Get('hden_monojet_recoil_clone_passed')

#btag_eff = TFile('btag_eff_forweight.root')
#btag_eff_light = btag_eff.Get('efficiency_light')
#btag_eff_b = btag_eff.Get('efficiency_b')
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

parser.add_option("--se", action="store_true",  dest="SE")
parser.add_option("--met", action="store_true",  dest="MET")
parser.add_option("--sp", action="store_true",  dest="SP")

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

if options.SE==None:
    options.SE==False

if options.MET==None:
    options.MET==False

if options.SP==None:
    options.SP==False


if options.SE: print "Using SingleElectron dataset."
if options.SP: print "Using SinglePhoton dataset."
if options.MET: print "Using MET dataset."

#if not options.SE and not options.MET and not options.SP:
    #print "Please run using --se or --met or --sp. Exiting."
    #sys.exit()

if options.CSV: print "Using CSVv2 as b-tag discriminator."
if options.DeepCSV: print "Using DeepCSV as b-tag discriminator."

if not options.CSV and not options.DeepCSV:
    print "Please run using --csv or --deepcsv. Exiting."
    sys.exit()

applydPhicut=False

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
    if not isfarmout:
        outputfilename = "/Output_"+rootfile
    else:
        outputfilename = "/Output_"+rootfile.split('.')[0]+".root"
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
    #print "my file :  ", filename
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

triglist=["HLT_PFMET170_BeamHaloCleaned_v","HLT_PFMET170_HBHE_BeamHaloCleaned_v","HLT_PFMET170_NotCleaned_v","HLT_PFMET170_NoiseCleaned_v","HLT_PFMET170_JetIdCleaned_v","HLT_PFMET170_HBHECleaned_v","HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v","HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v","HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v","HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v","HLT_PFMET110_PFMHT110_IDTight_v","HLT_IsoMu24_v","HLT_IsoTkMu24_v","HLT_IsoMu27_v","HLT_IsoTkMu27_v","HLT_Ele27_WPTight_Gsf_v","HLT_Ele105_CaloIdVT_GsfTrkIdT_v","HLT_Ele115_CaloIdVT_GsfTrkIdT_v","HLT_Ele32_WPTight_Gsf_v","HLT_IsoMu20_v","HLT_Ele27_eta2p1_WPTight_Gsf_v","HLT_Ele27_WPLoose_Gsf_v","HLT_Ele32_eta2p1_WPTight_Gsf_v","HLT_Photon165_HE10_v","HLT_Photon175_v"]

h_t = TH1F('h_t','h_t',2,0,2)
h_t_weight = TH1F('h_t_weight','h_t_weight',2,0,2)

samplename = 'all'
if isfarmout:
    infile = open(inputfilename)
    failcount=0
    for ifile in infile:
        try:
            print (ifile.rstrip())
            f_tmp = TFile.Open(ifile.rstrip(),'READ')
            if f_tmp.IsZombie():            # or fileIsCorr(ifile.rstrip()):
                failcount += 1
                continue
            skimmedTree.Add(ifile.rstrip())
            h_tmp = f_tmp.Get('h_total')
            h_tmp_weight = f_tmp.Get('h_total_mcweight')
            h_t.Add(h_tmp)
            h_t_weight.Add(h_tmp_weight)
        except:
            failcount += 1
    if failcount>0: print "Could not read %d files. Skipping them." %failcount

if not isfarmout:
    skimmedTree.Add(inputfilename)
    f_tmp = TFile.Open(inputfilename,'READ')
    h_tmp = f_tmp.Get('h_total')
    h_tmp_weight = f_tmp.Get('h_total_mcweight')
    h_t.Add(h_tmp)
    h_t_weight.Add(h_tmp_weight)

debug = False

try:
    #print ('samplepath',str(f_tmp).split('/')[-2].strip('_000')+'.root')
    samplepath = str(f_tmp.Get('samplepath').GetTitle())
    if not isfarmout: print "Original source file: " + samplepath
except:
#    samplepath=inputfilename
    samplepath='TT'
    print "WARNING: Looks like the input was skimmed with an older version of SkimTree. Using " + samplepath + " as sample path. Gen pT Reweighting may NOT work."

samplename = WhichSample(samplepath)
print "Dataset classified as: " + samplename
UnPrescaledIsoMu20 = IsoMu20isUnPrescaled(samplepath)
#print UnPrescaledIsoMu20
#print samplename
#btag_eff = TFile('btag_eff_forweight.root')
print('Eff_Files_v2/'+str(f_tmp).split('/')[-2].strip('_000')+'.root')
data_file = False
if '190610_09' in 'Eff_Files/'+str(f_tmp).split('/')[-2].strip('_000')+'.root': data_file=True
if 'SingleElectron' in 'Eff_Files/'+str(f_tmp).split('/')[-2].strip('_000')+'.root': data_file=True
if data_file:
    btag_eff = TFile('btag_eff_forweight.root')
else:
    btag_eff = TFile('Eff_Files_v2/'+str(f_tmp).split('/')[-2].strip('_000')+'.root')

eff_b_pass = btag_eff.Get('efficiency_beff_pass')
eff_b_fail = btag_eff.Get('efficiency_beff_fail')

eff_c_pass = btag_eff.Get('efficiency_ceff_pass')
eff_c_fail  = btag_eff.Get('efficiency_ceff_fail')

eff_light_pass = btag_eff.Get('efficiency_lighteff_pass')
eff_light_fail = btag_eff.Get('efficiency_lighteff_fail')
#btag_eff_light = btag_eff.Get('efficiency_light')
#btag_eff_b = btag_eff.Get('efficiency_b')

def AnalyzeDataSet():
    ## Input rootfile name

    #rootfilename = inputfilename
    #print (rootfilename,inputfilename)
    #f = TFile(rootfilename,'READ')
    #skimmedTree = f.Get('tree/treeMaker')
    NEntries = skimmedTree.GetEntries()
    print 'NEntries(total) =', skimmedTree.GetEntries()
    #NEntries = 1000
    print 'NEntries = '+str(NEntries)
    npass = 0
    npass_we = 0
    npass_wmu = 0


    #print [rootfilename, NEntries]
    cutStatus={'preselection':NEntries}
    cutStatusSR1={'preselection':NEntries}
    cutStatusSR2={'preselection':NEntries}

    cutflownames=['trig','MET','dPhicond','njets','nbjets','jet1','jet2/3','nlep']
    for SRname in cutflownames:
        cutStatus[SRname] = 0

    cutflownamesSR1=['trig','MET','dPhicond','njets','nbjets','jet1','jet2','nlep']
    for SRname in cutflownamesSR1:
        cutStatusSR1[SRname] = 0

    cutflownamesSR2=['trig','MET','dPhicond','njets','nbjets','jet1','jet2','jet3','nlep']
    for SRname in cutflownamesSR2:
        cutStatusSR2[SRname] = 0


    #CRcutnames=['trigger','leptons=1','njets>=3','nbjets=0','MET>250','mT<160']
    CRcutnames=['trigger','leptons=1','njets>=3','nbjets=0','mT<160']
    regionnames=['2e1b','2mu1b','2e2b','2mu2b','1e1b','1mu1b','1e2b','1mu2b','1etop1b','1mutop1b','1etop2b','1mutop2b','1mu1e1b','1mu1e2b','1gamma1b','1gamma2b','QCD1b','QCD2b']
    for CRreg in regionnames:
        exec("CR"+CRreg+"CutFlow={'preselection':NEntries}")
        for cutname in CRcutnames:
            #if ('top' in CRreg or '1e1b' in CRreg or '1e2b' in CRreg or '1mu1b' in CRreg or '1mu2b' in CRreg) and (not '1mu1e' in CRreg ) and 'njets' in cutname:
            #    cutname = 'add_Jet'
            exec("CR"+CRreg+"CutFlow['"+cutname+"']=0")


    CRs=['ZCRSR1','ZCRSR2','WCRSR1','WCRSR2','TOPCRSR1','TOPCRSR2','TopCRSR1','TopCRSR2', 'GammaCRSR1','GammaCRSR2']

    CRStatus={'total':NEntries}
    for CRname in CRs:
        CRStatus[CRname]=0

    # ---CR Summary---
    regNames=['1#mu1b','1e1b','1#mu2b','1e2b','1#mutop1b','1etop1b','1#mutop2b','1etop2b','2#mu1b','2e1b','2#mu2b','2e2b','1#mu1e1b','1#mu1e2b']
    regNamesMu=['1#mu1b','1#mu2b','2#mu1b','2#mu2b','1#mutop1b','1#mutop2b','1#mu1e1b','1#mu1e2b']
    regNamesEle=['1e1b','1e2b','1etop1b','1etop2b','2e1b','2e2b']

    CRSummary={}
    for ireg in regNames:
        CRSummary[ireg]=0.

    CRSummaryMu={}
    for ireg in regNamesMu:
        CRSummaryMu[ireg]=0.

    CRSummaryEle={}
    for ireg in regNamesEle:
        CRSummaryEle[ireg]=0.

    #print outfilename
    allquantities = bbMETQuantities(outfilename)
    allquantities.defineHisto()


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
        calib1 = ROOT.BTagCalibrationStandalone('deepcsv', 'DeepCSV_Moriond17_B_H.csv')

    reader1 = ROOT.BTagCalibrationStandaloneReader( 0, "central", othersys)
    reader1.load(calib1, 0,  "comb" )
    reader1.load(calib1, 1,  "comb" )
    reader1.load(calib1, 2,  "incl" )

#    h_total = TH1F('h_total','h_total',2,0,2)
#    h_total_mcweight = TH1F('h_total_mcweight','h_total_mcweight',2,0,2)


    e_num = 0
    e_num_tight = 0
    for ievent in range(NEntries):
    #for ievent in range(501):

        ##
        sf_resolved1 = []
        sf_resolved2 = []
        sf_resolved3 = []
        #print "event number = ",ievent
        skimmedTree.GetEntry(ievent)

        ## Get all relevant branches
#        if True:
        try:
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Extract branches
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
            run                        = skimmedTree.__getattr__('st_runId')
            lumi                       = skimmedTree.__getattr__('st_lumiSection')
            event                      = skimmedTree.__getattr__('st_eventId')
            '''
            event_ext = True
            #event_num = [10288545,15534196,15425356,17816400,25092178,22166398,14839367,1368918,3695802,17672510,4164556,21373626,17590235,24043129,22729804,24003428,20724527,3169418,25195362]
            #event_num = [18647061]
            #event_num = [2611079,18647061,1911047,23742750,9933019,10271640,5546451,5261651,20071771,3519529,18118536,12321484,22476635,10979565,22729804,24003428,25195362,22166398,1368918,17672510,21373626]
            event_num = [18647061,9933019,12321484]#,22729804,24003428,25195362,22166398,1368918,17672510,21373626]
            for ieve in event_num:
                if event==ieve:
                    event_ext = False
            if event_ext: continue
            '''

            #if event != 4126: continue
            #if lumi  != 42: continue
            #if ievent%100==0: print (ievent)
            #trigName                   = skimmedTree.__getattr__('st_hlt_trigName')
            #trigResult                 = skimmedTree.__getattr__('st_hlt_trigResult')
            #filterName                 = skimmedTree.__getattr__('st_hlt_filterName')
            #filterResult               = skimmedTree.__getattr__('st_hlt_filterResult')

            pfMet                      = skimmedTree.__getattr__('st_pfMetCorrPt')
            pfMetPhi                   = skimmedTree.__getattr__('st_pfMetCorrPhi')
            try:
                pfMetJetResUp              = skimmedTree.__getattr__('st_pfMetUncJetResUp')
                pfMetJetResDown            = skimmedTree.__getattr__('st_pfMetUncJetResDown')
                pfMetJetEnUp               = skimmedTree.__getattr__('st_pfMetUncJetEnUp')
                pfMetJetEnDown             = skimmedTree.__getattr__('st_pfMetUncJetEnDown')
            except:
                nullpfMet=-1
                pfMetJetResUp = nullpfMet
                pfMetJetResDown = nullpfMet
                pfMetJetEnUp = nullpfMet
                pfMetJetEnDown = nullpfMet

            nTHINJets                  = skimmedTree.__getattr__('st_THINnJet')
            thinjetP4                  = skimmedTree.__getattr__('st_THINjetP4')
            thinJetCSV                 = skimmedTree.__getattr__('st_THINjetCISVV2')
            #passThinJetLooseID         = skimmedTree.__getattr__('st_THINjetPassIDLoose')
            #passThinJetPUID            = skimmedTree.__getattr__('st_THINisPUJetID')
            THINjetHadronFlavor        = skimmedTree.__getattr__('st_THINjetHadronFlavor')
            thinjetNhadEF              = skimmedTree.__getattr__('st_THINjetNHadEF')
            thinjetChadEF              = skimmedTree.__getattr__('st_THINjetCHadEF')
            thinjetNPV                 = skimmedTree.__getattr__('st_THINjetNPV')
            thinjetCorrUnc             = skimmedTree.__getattr__('st_THINjetCorrUnc')

            try:
                thinjetCEmEF               = skimmedTree.__getattr__('st_THINjetCEmEF')
                thinjetPhoEF               = skimmedTree.__getattr__('st_THINjetPhoEF')
                thinjetEleEF               = skimmedTree.__getattr__('st_THINjetEleEF')
                thinjetMuoEF               = skimmedTree.__getattr__('st_THINjetMuoEF')

            except:
                nulljet=[-1. for i in range(nTHINJets)]
                thinjetCEmEF=nulljet
                thinjetPhoEF=nulljet
                thinjetEleEF=nulljet
                thinjetMuoEF=nulljet

            nTHINdeepCSVJets           = skimmedTree.__getattr__('st_AK4deepCSVnJet')
            thindeepCSVjetP4           = skimmedTree.__getattr__('st_AK4deepCSVjetP4')
            thinJetdeepCSV             = skimmedTree.__getattr__('st_AK4deepCSVjetDeepCSV_b')
            THINdeepCSVjetHadronFlavor = skimmedTree.__getattr__('st_AK4deepCSVjetHadronFlavor')
            thindeepCSVjetNhadEF       = skimmedTree.__getattr__('st_AK4deepCSVjetNHadEF')
            thindeepCSVjetChadEF       = skimmedTree.__getattr__('st_AK4deepCSVjetCHadEF')
            thindeepCSVjetNPV          = skimmedTree.__getattr__('st_AK4deepCSVjetNPV')

            try:
                thindeepCSVjetCEmEF               = skimmedTree.__getattr__('st_AK4deepCSVjetCEmEF')
                thindeepCSVjetPhoEF               = skimmedTree.__getattr__('st_AK4deepCSVjetPhoEF')
                thindeepCSVjetEleEF               = skimmedTree.__getattr__('st_AK4deepCSVjetEleEF')
                thindeepCSVjetMuoEF               = skimmedTree.__getattr__('st_AK4deepCSVjetMuoEF')
                thindeepCSVjetCorrUnc             = skimmedTree.__getattr__('st_AK4deepCSVjetCorrUnc')

            except:
                nulljet=[-1. for i in range(nTHINdeepCSVJets)]
                thindeepCSVjetCEmEF=nulljet
                thindeepCSVjetPhoEF=nulljet
                thindeepCSVjetEleEF=nulljet
                thindeepCSVjetMuoEF=nulljet
                thindeepCSVjetCorrUnc = nulljet

            nPho                       = skimmedTree.__getattr__('st_nPho')
            phoP4                      = skimmedTree.__getattr__('st_phoP4')
            phoIsPassTight             = skimmedTree.__getattr__('st_phoIsPassTight')

            nEle                       = skimmedTree.__getattr__('st_nEle')
            eleP4                      = skimmedTree.__getattr__('st_eleP4')
            eleIsPassTight             = skimmedTree.__getattr__('st_eleIsPassTight')

            nMu                        = skimmedTree.__getattr__('st_nMu')
            muP4                       = skimmedTree.__getattr__('st_muP4')
            isTightMuon                = skimmedTree.__getattr__('st_isTightMuon')
            muIso                      = skimmedTree.__getattr__('st_muIso')

            nTau                       = skimmedTree.__getattr__('st_HPSTau_n')
            nTauTightElectron          = skimmedTree.__getattr__('st_nTauTightElectron')
            nTauTightMuon              = skimmedTree.__getattr__('st_nTauTightMuon')
            nTauTightEleMu             = skimmedTree.__getattr__('st_nTauTightEleMu')
            nTauLooseEleMu             = skimmedTree.__getattr__('st_nTauLooseEleMu')

            isData                     = skimmedTree.__getattr__('st_isData')
            mcWeight                   = skimmedTree.__getattr__('mcweight')
            pu_nTrueInt                = int(skimmedTree.__getattr__('st_pu_nTrueInt'))
            pu_nPUVert                 = int(skimmedTree.__getattr__('st_pu_nPUVert'))

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
            GammaRecoil                = skimmedTree.__getattr__('GammaRecoil')
            GammaPhi                   = skimmedTree.__getattr__('GammaPhi')

            for trig in triglist:
                exec(trig+" = skimmedTree.__getattr__('st_"+trig+"')")

        except Exception as e:
#        else:
            print e
            print "Corrupt file detected! Skipping 1 event."
            continue

        ##Define region wise triggers

        if isData:
            SRtrigstatus = HLT_PFMET170_BeamHaloCleaned_v or HLT_PFMET170_HBHE_BeamHaloCleaned_v or HLT_PFMET170_NotCleaned_v or HLT_PFMET170_NoiseCleaned_v or HLT_PFMET170_JetIdCleaned_v or HLT_PFMET170_HBHECleaned_v or HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v or HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v or HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v or HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v or HLT_PFMET110_PFMHT110_IDTight_v

            MuCRtrigstatus =  HLT_IsoMu24_v or HLT_IsoTkMu24_v or HLT_IsoTkMu27_v or HLT_IsoMu27_v

            EleCRtrigstatus = HLT_Ele27_WPTight_Gsf_v or HLT_Ele105_CaloIdVT_GsfTrkIdT_v or HLT_Ele115_CaloIdVT_GsfTrkIdT_v or HLT_Ele32_WPTight_Gsf_v or HLT_Ele27_eta2p1_WPTight_Gsf_v or HLT_Ele27_WPLoose_Gsf_v or HLT_Ele32_eta2p1_WPTight_Gsf_v

            PhotonCRtrigstatus = HLT_Photon165_HE10_v or HLT_Photon175_v

        else:
            SRtrigstatus = True
            MuCRtrigstatus = True
            EleCRtrigstatus = True
            PhotonCRtrigstatus = True


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
        ## PFMET Selection
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        pfmetstatus = ( pfMet >= 250.0 )
        #----------------------------------------------------------------------------------------------------------------------------------------------------------------

         # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
#        print (HLT_IsoMu24,HLT_Ele27_WPLoose_Gsf)

        crack_reg = (abs(eleP4[iele].Eta()) > 1.4442) and (abs(eleP4[iele].Eta()) < 1.566)
        myEles=[]
        myEleTightID=[]
        myEles_tight=[]
        myEles_Wmass = []
        for iele in range(nEle):
            e_num = e_num + 1
            if abs(eleP4[iele].Eta()) >2.5 and crack_reg: continue:
                myEles.append(eleP4[iele])
            #myEleTightID.append(eleIsPassTight[iele])
            if bool(eleIsPassTight[iele])==True and eleP4[iele].Pt() > 30. and abs(eleP4[iele].Eta())<2.1:
                e_mass = MT(eleP4[iele].Pt(),pfMet, DeltaPhi(eleP4[iele].Phi(),pfMetPhi))
                myEles_tight.append(eleP4[iele])
                myEles_Wmass.append(e_mass)
                myEleTightID.append(eleIsPassTight[iele])

        myMuos = []
        myMuTightID=[]
        myMuIso=[]
        myMuos_tight=[]
        myMuos_Wmass = []
        for imu in range(nMu):
            if abs(muP4[imu].Eta()) > 2.4  : continue
            myMuos.append(muP4[imu])
            #myMuTightID.append(isTightMuon[imu])
            myMuIso.append(muIso[imu])
            if bool(isTightMuon[imu])==True and (muIso[imu] < 0.15) and muP4[imu].Pt() > 30.:
                mu_mass = MT(muP4[imu].Pt(),pfMet, DeltaPhi(muP4[imu].Phi(),pfMetPhi))
                myMuos_tight.append(muP4[imu])
                myMuos_Wmass.append(mu_mass)
                myMuTightID.append(isTightMuon[imu])
            #print 'muID: ', bool(isTightMuon[imu]),'muIso: ',muIso[imu]

        #--------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # Jet Selection
        ## Also segregate CSV or DeepCSV collection of jets in this step itself

        CSVMWP=0.8484
        deepCSVMWP=0.6324

        mybjets=[]
        myJetCSV=[]
        myJetP4=[]
        myJetHadronFlavor=[]
        myJetNhadEF=[]
        myJetChadEF=[]
        myJetCorrUnc=[]
        myJetCEmEF=[]
        myJetPhoEF=[]
        myJetEleEF=[]
        myJetMuoEF=[]
        #print '\n'
        if options.CSV: 
            for nb in range(nTHINJets):
                TightLepVeto = (thinjetNhadEF[nb]<0.90 and thinjetCEmEF[nb]<0.90 and thinjetMuoEF[nb]<0.8 and thinjetPhoEF[nb]<0.90)
                '''
                #---Fake jet cleaner, wrt electrons and muons----
                isClean=True
                for iele in range(len(myEles)):
                    if DeltaR(myEles[iele],thinjetP4[nb]) < 0.4 :
                        isClean=False
                        break
                for imu in range(len(myMuos)):
                    if DeltaR(myMuos[imu],thinjetP4[nb]) < 0.4 :
                        isClean=False
                        break
                if not isClean: continue
                #---
                '''
                #print 'thinjetNhadEF ',thinjetNhadEF[nb],'thinjetCEmEF ', thinjetCEmEF[nb],'thinjetMuoEF ', thinjetMuoEF[nb],'thinjetPhoEF ', thinjetPhoEF[nb]
                if not TightLepVeto: continue
                #print 'jet'+str(nb)+' pT', thinjetP4[nb].Pt()
                #print 'jet'+str(nb)+' Eta', thinjetP4[nb].Eta()
                #print 'jet'+str(nb)+' phi', thinjetP4[nb].Phi()
                myJetP4.append(thinjetP4[nb])
                myJetCSV.append(thinJetCSV[nb])
                myJetHadronFlavor.append(THINjetHadronFlavor[nb])
                myJetNhadEF.append(thinjetNhadEF[nb])
                myJetChadEF.append(thinjetChadEF[nb])

                myJetCorrUnc.append(thinjetCorrUnc[nb])
                myJetCEmEF.append(thinjetCEmEF[nb])
                myJetPhoEF.append(thinjetPhoEF[nb])
                myJetEleEF.append(thinjetEleEF[nb])
                myJetMuoEF.append(thinjetMuoEF[nb])
                #print 'jetHadronFlavor: ', THINjetHadronFlavor[nb],', JetCSV: ', thinJetCSV[nb], ',jetPt: ',thinjetP4[nb].Pt(),',jetEta: ',thinjetP4[nb].Eta()
                if thinJetCSV[nb] > CSVMWP and abs(thinjetP4[nb].Eta())<2.4:
                    mybjets.append(nb)

            myJetNPV=thinjetNPV
            nUncleanJets=nTHINJets
        #print 'bjet: ', len(mybjets)
        #print'\n'i
        if options.DeepCSV:
            for nb in range(nTHINdeepCSVJets):

            #---Fake jet cleaner, wrt electrons and muons----
                isClean=True

                for iele in myEles:
                    if DeltaR(iele,thindeepCSVjetP4[nb]) < 0.4:
                        isClean=False
                        break
                for imu in myMuos:
                    if DeltaR(imu,thindeepCSVjetP4[nb]) < 0.4:
                        isClean=False
                        break


                if not isClean: continue
            #---

                myJetP4.append(thindeepCSVjetP4[nb])
                myJetCSV.append(thinJetdeepCSV[nb])
                myJetHadronFlavor.append(THINdeepCSVjetHadronFlavor[nb])
                myJetNhadEF.append(thindeepCSVjetNhadEF[nb])
                myJetChadEF.append(thindeepCSVjetChadEF[nb])

                myJetCorrUnc.append(thindeepCSVjetCorrUnc[nb])

                if ievent==0: print "Jet Energy fractions are not saved in DeepCSV collection. Saving default -1 for these."
                myJetCEmEF.append(thindeepCSVjetCEmEF[nb])
                myJetPhoEF.append(thindeepCSVjetPhoEF[nb])
                myJetEleEF.append(thindeepCSVjetEleEF[nb])
                myJetMuoEF.append(thindeepCSVjetMuoEF[nb])

                if thinJetdeepCSV[nb] > deepCSVMWP and abs(thindeepCSVjetP4[nb].Eta())<2.4:
                    mybjets.append(nb)

            myJetNPV=thindeepCSVjetNPV
            nUncleanJets=nTHINdeepCSVJets

        myPhos=[]
        myPhoTightID=[]
        for ipho in range(nPho):
            if phoP4[ipho].Pt() < 175 : continue
            #---Fake Pho cleaner----
            isClean=True
            for ijet in myJetP4[:]:
                pho_jet_dR=DeltaR(ijet,phoP4[ipho])    # math.sqrt(  (  ijet.Eta()-phoP4[ipho].Eta() )**2  + (  DeltaPhi(ijet.Phi(),phoP4[ipho].Phi()) )**2 )
                if pho_jet_dR < 0.4:
                    isClean=False
                    break
            for iele in myEles[:]:
                lep_pho_dR=DeltaR(iele,phoP4[ipho])
                if lep_pho_dR < 0.4:
                    isClean=False
                    break
            if not isClean: continue
            myPhos.append(phoP4[ipho])
            myPhoTightID.append(phoIsPassTight[ipho])



        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------

        nUncleanEle=nEle
        nUncleanMu=nMu
        nUncleanTau=nTau

        nPho=len(myPhos)
        nEle=len(myEles)
        nEle_tight=len(myEles_tight)
        nMu=len(myMuos)
        nMu_tight=len(myMuos_tight)
#        nTau=len(myTaus)

        nTau=nTauLooseEleMu

#        print "Unclean,dR-based,Loose,ETight,MTight,EMTight:",nUncleanTau,nTausDRbased,nTauLooseEleMu,nTauTightElectron,nTauTightMuon, nTauTightEleMu

        nBjets=len(mybjets)
        nJets=len(myJetCSV)
 #----------------------------------------------------------------------------------------------------------------------------------------------------------------
 #----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Sort jets

        if nJets==0: continue

        alljetPT=[jet.Pt() for jet in myJetP4]
        jetindex=[i for i in range(len(alljetPT))]

        sortedjets=[jet for pt,jet in sorted(zip(alljetPT,myJetP4), reverse=True)]      # This gives a list of jets with their pTs in descending order
        sortedindex=[jetindex for pt,jetindex in sorted(zip(alljetPT,jetindex), reverse=True)]     # Indices of jets in myJetP4 in decscending order of jetPT

        j1=sortedjets[0]
        if nJets>1: j2=sortedjets[1]
        if nJets>2: j3=sortedjets[2]

        ifirstjet=sortedindex[0]
        if nJets>1: isecondjet=sortedindex[1]
        if nJets>2: ithirdjet=sortedindex[2]

        min_dPhi_jet_MET = min(   [  DeltaPhi(jt.Phi(),pfMetPhi) for jt in myJetP4]   )
#        print alljetPT
#        print [jet.Pt() for jet in sortedjets]
#        print sortedindex
#        print

        ##
# --------------------------------------------------------------------------------------------------------------------------------------------------------

        allquantlist=AllQuantList.getAll()

        for quant in allquantlist:
            exec("allquantities."+quant+" = None")
# --------------------------------------------------------------------------------------------------------------------------------------------------------


        # SR start

        ###
        #********
        #Part of data is blinded
        #********
        if isData:
            if ievent%20==0:
                keepevent=True
            else:
                keepevent=False
        else:
            keepevent=True
        #********************              REMEMBER TO UNBLIND AT SOME POINT!

        writeSR1=False
        writeSR2=False
        isZeeCR1=False
        isZeeCR2=False
        isZmumuCR1=False
        isZmumuCR2=False
        isWenuCR1=False
        isWenuCR2=False
        isWmunuCR1=False
        isWmunuCR2=False
        isTopenuCR1=False
        isTopenuCR2=False
        isTopmunuCR1=False
        isTopmunuCR2=False
        isTopCR1=False
        isTopCR2=False
        isGammaCR1=False
        isGammaCR2=False



        SR1njetcond=False
        SR2njetcond=False
        SR1jetcond=False
        SR2jetcond=False


        if nEle+nMu+nTauLooseEleMu==0:
            SRlepcond=True
        else:
            SRlepcond=False

        ## for SR1
         # 1 or 2 jets and 1 btagged

        Histos2D=AllQuantList.getHistos2D()
        for quant in Histos2D:
            exec("allquantities."+quant+" = None")

        if (nJets == 1 or nJets == 2) and pfmetstatus and SRlepcond and SRtrigstatus:
            #===CSVs before any selection===
            preselquantlist=AllQuantList.getPresel()
            for quant in preselquantlist:
                exec("allquantities."+quant+" = None")

            if options.CSV:
                allquantities.presel_jet1_csv_sr1=myJetCSV[ifirstjet]
                if nJets>1: allquantities.presel_jet2_csv_sr1=myJetCSV[isecondjet]
            if options.DeepCSV:
                allquantities.presel_jet1_deepcsv_sr1=myJetCSV[ifirstjet]
                if nJets>1: allquantities.presel_jet2_deepcsv_sr1=myJetCSV[isecondjet]

            allquantities.presel_jet1_chf_sr1=myJetChadEF[ifirstjet]
            allquantities.presel_jet1_nhf_sr1=myJetNhadEF[ifirstjet]
            allquantities.FillPreSel()



            if options.CSV:
                allquantities.csv_vs_dPhi_sr1 = [min_dPhi_jet_MET,myJetCSV[ifirstjet]]
            if options.DeepCSV:
                allquantities.deepcsv_vs_dPhi_sr1 = [min_dPhi_jet_MET,myJetCSV[ifirstjet]]
            #===
        if (nJets == 1 or nJets == 2):
            SR1jetcond=True
        if (nJets == 1 or nJets == 2) and nBjets==1:
            SR1njetcond=True

        SR1_Cut1_nJets          =   nJets == 1 or nJets == 2
        SR1_Cut2_nBjets         =   nBjets==1
        SR1_Cut3_trigstatus     =   SRtrigstatus
        SR1_Cut4_jet1           =   j1.Pt() > 50.0 and myJetNhadEF[ifirstjet] < 0.8 and myJetChadEF[ifirstjet] > 0.1
        if nJets>1:
            SR1_Cut5_jet2       =   j2.Pt() > 30.0
        else:
            SR1_Cut5_jet2       =   True
        SR1_Cut6_dPhi_jet_MET   =   min_dPhi_jet_MET > 0.5
        SR1_Cut7_nLep           =   nEle+nMu+nTauLooseEleMu == 0
        SR1_Cut8_pfMET          =   pfmetstatus




     ## for SR2
        # 3 jets and 2 btagged

        if (nJets == 2 or nJets == 3) and pfmetstatus and SRlepcond and SRtrigstatus:
            #===CSVs before any selection===
            preselquantlist=AllQuantList.getPresel()
            for quant in preselquantlist:
                exec("allquantities."+quant+" = None")

            if options.CSV:
                allquantities.presel_jet1_csv_sr2=myJetCSV[ifirstjet]
                allquantities.presel_jet2_csv_sr2=myJetCSV[isecondjet]
                if nJets>2: allquantities.presel_jet3_csv_sr2=myJetCSV[ithirdjet]

            if options.DeepCSV:
                allquantities.presel_jet1_deepcsv_sr2=myJetCSV[ifirstjet]
                allquantities.presel_jet2_deepcsv_sr2=myJetCSV[isecondjet]
                if nJets>2: allquantities.presel_jet3_deepcsv_sr2=myJetCSV[ithirdjet]

            allquantities.presel_jet1_chf_sr2=myJetChadEF[ifirstjet]
            allquantities.presel_jet1_nhf_sr2=myJetNhadEF[ifirstjet]
            allquantities.FillPreSel()


            #===
            if options.CSV:
                csv_ = min(myJetCSV[ifirstjet],myJetCSV[isecondjet])
                allquantities.csv_vs_dPhi_sr2 = [min_dPhi_jet_MET,csv_]

            if options.DeepCSV:
                deepcsv_ = min(myJetCSV[ifirstjet],myJetCSV[isecondjet])
                allquantities.deepcsv_vs_dPhi_sr2 = [min_dPhi_jet_MET,deepcsv_]

        if (nJets == 2 or nJets == 3):
            SR2jetcond=True

        if (nJets == 2 or nJets == 3) and nBjets==2:
            SR2njetcond=True

        SR2_Cut1_nJets          =   nJets == 2 or nJets == 3
        SR2_Cut2_nBjets         =   nBjets==2
        SR2_Cut3_trigstatus     =   SRtrigstatus
        SR2_Cut4_jet1           =   j1.Pt() > 50.0 and myJetNhadEF[ifirstjet] < 0.8 and myJetChadEF[ifirstjet] > 0.1
        if nJets>1:
            SR2_Cut5_jet2           =   j2.Pt() > 50.0
        else:
            SR2_Cut5_jet2           =   False
        if nJets>2:
            SR2_Cut6_jet3       =   j3.Pt() > 30.0
        else:
            SR2_Cut6_jet3       =   True
        SR2_Cut7_dPhi_jet_MET   =   min_dPhi_jet_MET > 0.5
        SR2_Cut8_nLep           =   nEle+nMu+nTauLooseEleMu == 0
        SR2_Cut9_pfMET          =   pfmetstatus






# --------------------------------------------------------------------------------------------------------------------------------------------------------

        #Control Regions

# --------------------------------------------------------------------------------------------------------------------------------------------------------

        preselquantlist=AllQuantList.getPresel()

        for quant in preselquantlist:
            exec("allquantities."+quant+" = None")


        regquants=AllQuantList.getRegionQuants()

        for quant in regquants:
            exec("allquantities."+quant+" = None")

        # Histos2D=AllQuantList.getHistos2D()
        # for quant in Histos2D:
        #     exec("allquantities."+quant+" = None")



        ####new conds
        jetcond=True
        SR2jet2=True

        #if myJetNhadEF[ifirstjet] > 0.8 : jetcond=False
        #if myJetChadEF[ifirstjet]< 0.1: jetcond=False
        for ijet in sortedjets:
            if ijet.Pt() < 30.0: jetcond=False
        #for i in sortedindex:
        #    if myJetNhadEF[i] > 0.8 : jetcond=False
        #    if myJetChadEF[i]< 0.1: jetcond=False


# -------------------------------------------
# W CR
# -------------------------------------------
        iLeadLep = 0
        #print run,' : ',lumi,' : ',event,' : nJets ', nJets,' : nBjets ',nBjets
        ''' 
        if nJets==2:
            print nJets,' : ',run, ' : ', lumi, ' : ', event, ' : ',myEles_tight[iLeadLep].Pt(),' : ',myEles_tight[iLeadLep].Eta(),' : ',myEles_tight[iLeadLep].Phi(),' : ',j1.Pt(),' : ',j1.Eta(),' : ',j1.Phi(),' : ',j2.Pt(),' : ',j2.Eta(),' : ',j2.Phi(),' : ',pfMet,' : ',myEles_Wmass[iLeadLep]
        elif nJets >=3:
            print nJets,' : ',run, ' : ', lumi, ' : ', event, ' : ',myEles_tight[iLeadLep].Pt(),' : ',myEles_tight[iLeadLep].Eta(),' : ',myEles_tight[iLeadLep].Phi(),' : ',j1.Pt(),' : ',j1.Eta(),' : ',j1.Phi(),' : ',j2.Pt(),' : ',j2.Eta(),' : ',j2.Phi(),' : ',j3.Pt(),' : ',j3.Eta(),' : ',j3.Phi(),' : ',pfMet,' : ',myEles_Wmass[iLeadLep]
        #print run, ' : ', lumi, ' : ', event, ' : ',myEles_tight[iLeadLep].Pt(),' : ',myEles_tight[iLeadLep].Eta(),' : ', myEles[iLeadLep].Phi(),' : ',j1.Pt(),' : ',j1.Eta(),' : ',j1.Phi(),' : ',j2.Pt(),' : ',j2.Eta(),' : ',j2.Phi(),' : ',j3.Pt(),' : ',j3.Eta(),' : ',j3.Phi(),' : ',pfMet,' : ',Wenumass
        if len(myEles_tight) ==2:
            print 'ele1 pT:', myEles_tight[0].Pt(),'ele2 pT:', myEles_tight[1].Pt()
        '''
        #W CR specific bools
        tDMJetCond=False
        if nJets>=3 and nBjets==0:
            tDMJetCond=True
        Wgn_pT = getGenPt(samplename, nGenPar, genParId, genMomParId, genParSt,genParP4)
        WdPhicond=True
        WpT = 1.0

        if applydPhicut:
            if WenuPhi>-10.:
                if min( [DeltaPhi(WenuPhi,myJetP4[nb].Phi()) for nb in range(nJets)] ) < 0.5: WdPhicond=False
            if WmunuPhi>-10.:
                if min( [DeltaPhi(WmunuPhi,myJetP4[nb].Phi()) for nb in range(nJets)] ) < 0.5: WdPhicond=False

        #1e, 1 b-tagged
        #if nEle==1 and nMu==0 and nTauTightElectron==0 and EleCRtrigstatus and WenuRecoil>200. and jetcond and Wenumass>50. and Wenumass<160. and pfMet > 50.:
        #if nEle==1 and nMu==0 and nTauTightElectron==0 and EleCRtrigstatus and WenuRecoil>200. and jetcond and Wenumass<160. and pfMet > 50.:
        #if nEle==1 and nMu==0 and nTauTightElectron==0 and EleCRtrigstatus and  jetcond and Wenumass <= 160. and pfmetstatus:
        #if nEle_tight==1 and nMu==0 and EleCRtrigstatus and  jetcond and Wenumass <= 160. and pfmetstatus:
        #print nEle_tight,' : ',nMu,' : ',EleCRtrigstatus,' : ',jetcond,' : ',pfmetstatus
        if nEle_tight==1 and nMu==0 and EleCRtrigstatus and  jetcond and pfmetstatus:
            iLeadLep=0

            # if myEles[iLeadLep].Pt() > 10. and myEleTightID[iLeadLep]:
            #     WpT = math.sqrt( ( pfMet*math.cos(pfMetPhi) + myEles[iLeadLep].Px())**2 + ( pfMet*math.sin(pfMetPhi) + myEles[iLeadLep].Py())**2)
            #     if  nBjets==1 and WdPhicond and (nJets-nBjets)==0:
            #         allquantities.reg_1e1b_Wmass_ptTen = Wenumass
            #         allquantities.reg_1e1b_WpT_ptTen=WpT
            #         allquantities.reg_1e1b_hadrecoil_ptTen = WenuRecoil
            #     if  nBjets==2 and SR2jet2 and WdPhicond and (nJets-nBjets)==0:
            #         allquantities.reg_1e2b_Wmass_ptTen = Wenumass
            #         allquantities.reg_1e2b_WpT_ptTen=WpT
            #         allquantities.reg_1e2b_hadrecoil_ptTen = WenuRecoil
            #
            # if myEles[iLeadLep].Pt() > 30. and myEleLooseID[iLeadLep]:
            #     WpT = math.sqrt( ( pfMet*math.cos(pfMetPhi) + myEles[iLeadLep].Px())**2 + ( pfMet*math.sin(pfMetPhi) + myEles[iLeadLep].Py())**2)
            #     if  nBjets==1 and WdPhicond and (nJets-nBjets)==0:
            #         allquantities.reg_1e1b_Wmass_looseID = Wenumass
            #         allquantities.reg_1e1b_WpT_looseID=WpT
            #         allquantities.reg_1e1b_hadrecoil_looseID = WenuRecoil
            #     if  nBjets==2 and SR2jet2 and WdPhicond and (nJets-nBjets)==0:
            #         allquantities.reg_1e2b_Wmass_looseID = Wenumass
            #         allquantities.reg_1e2b_WpT_looseID=WpT
            #         allquantities.reg_1e2b_hadrecoil_looseID = WenuRecoil
            #print myEles_Wmass[iLeadLep]
            #print myEles_tight[iLeadLep].Pt(),' : ',myEles_Wmass[iLeadLep]
            if myEles_tight[iLeadLep].Pt() > 30. and myEles_Wmass[iLeadLep] <= 160.:

                WpT = math.sqrt( ( pfMet*math.cos(pfMetPhi) + myEles[iLeadLep].Px())**2 + ( pfMet*math.sin(pfMetPhi) + myEles[iLeadLep].Py())**2)

                # if nBjets==1 and WdPhicond:
                #     allquantities.reg_1e1b_njet_n_minus_1=nJets
                #     allquantities.reg_1e1b_unclean_njet_n_minus_1=nUncleanJets
                #
                # if nBjets==1 and SR1njetcond:
                #     allquantities.reg_1e1b_min_dPhi_jet_Recoil_n_minus_1 = min( [DeltaPhi(WenuPhi,myJetP4[nb].Phi()) for nb in range(nJets)] )

                #if nBjets==1 and SR1njetcond and WdPhicond:
                #print nBjets
                if  nJets>=3 and nBjets==0:
                    npass_we = npass_we + 1
                    #print 'run', run, ': lumi', lumi, ': event', event, ': elept ',myEles[iLeadLep].Pt(),': eleeta',myEles[iLeadLep].Eta(),': elephi',myEles[iLeadLep].Phi(),': jet1pt',j1.Pt(),': jet1eta',j1.Eta(),': jet1phi',j1.Phi(),': jet2pt',j2.Pt(),': jet2eta',j2.Eta(), 
                    #print 'run', run, ' : lumi', lumi, ' : event', event, ' : elept ',myEles[iLeadLep].Pt(),' : eleeta',myEles[iLeadLep].Eta(),' : elephi',myEles[iLeadLep].Phi(),' : jet1pt',j1.Pt(),' : jet1eta',j1.Eta(),' : jet1phi',j1.Phi(),' : jet2pt',j2.Pt(),' : jet2eta',j2.Eta(),' : jet2phi',j2.Eta(),' : jet3pt',j3.Pt(),' : jet3eta',j3.Eta(),' : jet1phi',j3.Phi(),' : met',pfMet,' : mt',Wenumass
                    #print run, ' : ', lumi, ' : ', event, ' : ',myEles[iLeadLep].Pt(),' : ',myEles[iLeadLep].Eta(),' : ',myEles[iLeadLep].Phi(),' : ',j1.Pt(),' : ',j1.Eta(),' : ',j1.Phi(),' : ',j2.Pt(),' : ',j2.Eta(),' : ',j2.Phi(),' : ',j3.Pt(),' : ',j3.Eta(),' : ',j3.Phi(),' : ',pfMet,' : ',Wenumass
                    #print run,lumi,event
                    #print 'Selected event', event
                    #print '\n'
                    allquantities.reg_1e1b_Wmass = myEles_Wmass[iLeadLep]
                    allquantities.reg_1e1b_WpT=WpT

                    allquantities.reg_1e1b_hadrecoil = WenuRecoil
                    allquantities.reg_1e1b_MET = pfMet

                    allquantities.reg_1e1b_lep1_pT=myEles[iLeadLep].Pt()

                    allquantities.reg_1e1b_jet1_pT=j1.Pt()


                    if nJets>1: allquantities.reg_1e1b_jet2_pT=j2.Pt()

                    allquantities.reg_1e1b_jet1_eta=j1.Eta()
                    if nJets>1: allquantities.reg_1e1b_jet2_eta=j2.Eta()

                    allquantities.reg_1e1b_jet1_NHadEF=myJetNhadEF[ifirstjet]
                    allquantities.reg_1e1b_jet1_CHadEF=myJetChadEF[ifirstjet]
                    allquantities.reg_1e1b_jet1_CEmEF=myJetCEmEF[ifirstjet]
                    allquantities.reg_1e1b_jet1_PhoEF=myJetPhoEF[ifirstjet]
                    allquantities.reg_1e1b_jet1_EleEF=myJetEleEF[ifirstjet]
                    allquantities.reg_1e1b_jet1_MuoEF=myJetMuoEF[ifirstjet]

                    allquantities.reg_1e1b_njet = nJets

                    if options.CSV:
                        allquantities.reg_1e1b_jet1_csv = myJetCSV[ifirstjet]
                        if nJets>1: allquantities.reg_1e1b_jet2_csv = myJetCSV[isecondjet]

                        allquantities.reg_1e1b_min_dR_jet_ele_preclean = min( [DeltaR(myEles[iLeadLep],thinjetP4[nb]) for nb in range(nTHINJets)] )          #For diagnosis

                    if options.DeepCSV:
                        allquantities.reg_1e1b_jet1_deepcsv = myJetCSV[ifirstjet]
                        if nJets>1: allquantities.reg_1e1b_jet2_deepcsv = myJetCSV[isecondjet]

                    allquantities.reg_1e1b_min_dR_jet_ele_postclean = min( [DeltaR(myEles[iLeadLep],myJetP4[nb]) for nb in range(nJets)] )

                    allquantities.reg_1e1b_min_dPhi_jet_Recoil = min( [DeltaPhi(WenuPhi,myJetP4[nb].Phi()) for nb in range(nJets)] )
                    allquantities.reg_1e1b_min_dPhi_jet_MET = min( [DeltaPhi(pfMetPhi,myJetP4[nb].Phi()) for nb in range(nJets)] )

                    allquantities.reg_1e1b_ntau = nTauTightElectron
                    allquantities.reg_1e1b_nele = nEle
                    allquantities.reg_1e1b_nmu = nMu
                    allquantities.reg_1e1b_nUncleanTau = nUncleanTau
                    allquantities.reg_1e1b_nUncleanEle = nUncleanEle
                    allquantities.reg_1e1b_nUncleanMu = nUncleanMu
                    allquantities.lep_syst_1e1b_up = WenuRecoil
                    allquantities.lep_syst_1e1b_down = WenuRecoil
                    allquantities.btag_syst_1e1b_up = WenuRecoil
                    allquantities.btag_syst_1e1b_down = WenuRecoil
                    allquantities.metTrig_syst_1e1b_up =WenuRecoil
                    allquantities.metTrig_syst_1e1b_down =WenuRecoil
                    allquantities.ewkZ_syst_1e1b_up = WenuRecoil
                    allquantities.ewkZ_syst_1e1b_down = WenuRecoil
                    allquantities.ewkW_syst_1e1b_up = WenuRecoil
                    allquantities.ewkW_syst_1e1b_down = WenuRecoil
                    allquantities.ewkTop_syst_1e1b_up = WenuRecoil
                    allquantities.ewkTop_syst_1e1b_down = WenuRecoil
                    allquantities.pho_syst_1e1b_up = WenuRecoil
                    allquantities.pho_syst_1e1b_down = WenuRecoil
                    allquantities.jec_syst_1e1b_up = WenuRecoil
                    allquantities.jec_syst_1e1b_down = WenuRecoil
                    allquantities.jer_syst_1e1b_up = WenuRecoil
                    allquantities.jer_syst_1e1b_down = WenuRecoil
                    isWenuCR1 = True




        #print '\n'
        #1mu, 1 b-tagged
        #if nMu==1 and nEle==0 and nTauTightMuon==0 and MuCRtrigstatus and WmunuRecoil>200. and jetcond and Wmunumass>50. and Wmunumass<160. and pfMet > 50.:
        #if nMu==1 and nEle==0 and nTauTightMuon==0 and MuCRtrigstatus and WmunuRecoil>200. and jetcond and Wmunumass<160. and pfMet > 50.:
        #if nMu==1 and nEle==0 and  nTauTightMuon==0 and MuCRtrigstatus and jetcond and Wmunumass <= 160. and pfmetstatus:
        if nMu_tight==1 and nEle==0 and MuCRtrigstatus and jetcond and pfmetstatus:
            iLeadLep=0

            # if myMuos[iLeadLep].Pt() > 10. and myMuTightID[iLeadLep] and myMuIso[iLeadLep]<0.15:
            #     WpT = math.sqrt( ( pfMet*math.cos(pfMetPhi) + myMuos[iLeadLep].Px())**2 + ( pfMet*math.sin(pfMetPhi) + myMuos[iLeadLep].Py())**2)
            #     if  nBjets==1 and WdPhicond and (nJets-nBjets)==0:
            #         allquantities.reg_1mu1b_Wmass_ptTen = Wmunumass
            #         allquantities.reg_1mu1b_WpT_ptTen=WpT
            #         allquantities.reg_1mu1b_hadrecoil_ptTen = WmunuRecoil
            #     if  nBjets==2 and SR2jet2 and WdPhicond and (nJets-nBjets)==0:
            #         allquantities.reg_1mu2b_Wmass_ptTen = Wmunumass
            #         allquantities.reg_1mu2b_WpT_ptTen=WpT
            #         allquantities.reg_1mu2b_hadrecoil_ptTen = WmunuRecoil
            #
            # if myMuos[iLeadLep].Pt() > 30. and myMuLooseID[iLeadLep] and myMuIso[iLeadLep]<0.15:
            #     WpT = math.sqrt( ( pfMet*math.cos(pfMetPhi) + myMuos[iLeadLep].Px())**2 + ( pfMet*math.sin(pfMetPhi) + myMuos[iLeadLep].Py())**2)
            #     if  nBjets==1 and WdPhicond and (nJets-nBjets)==0:
            #         allquantities.reg_1mu1b_Wmass_looseID = Wmunumass
            #         allquantities.reg_1mu1b_WpT_looseID=WpT
            #         allquantities.reg_1mu1b_hadrecoil_looseID = WmunuRecoil
            #     if  nBjets==2 and SR2jet2 and WdPhicond and (nJets-nBjets)==0:
            #         allquantities.reg_1mu2b_Wmass_looseID = Wmunumass
            #         allquantities.reg_1mu2b_WpT_looseID=WpT
            #         allquantities.reg_1mu2b_hadrecoil_looseID = WmunuRecoil

            #if myMuos[iLeadLep].Pt() > 30. and myMuTightID[iLeadLep] and myMuIso[iLeadLep]<0.15:
            if myMuos_tight[iLeadLep].Pt() > 30. and myMuos_Wmass[iLeadLep] <= 160.:

                WpT = math.sqrt( ( pfMet*math.cos(pfMetPhi) + myMuos[iLeadLep].Px())**2 + ( pfMet*math.sin(pfMetPhi) + myMuos[iLeadLep].Py())**2)

                # if nBjets==1 and WdPhicond:
                #     allquantities.reg_1mu1b_njet_n_minus_1=nJets
                #     allquantities.reg_1mu1b_unclean_njet_n_minus_1=nUncleanJets
                #
                # if  nBjets==1 and SR1njetcond:
                #     allquantities.reg_1mu1b_min_dPhi_jet_Recoil_n_minus_1 = min( [DeltaPhi(WmunuPhi,myJetP4[nb].Phi()) for nb in range(nJets)] )

                #if  nBjets==1 and SR1njetcond and WdPhicond:
                #if  nBjets==1 and WdPhicond and (nJets-nBjets)==0:
                if  nJets>=3 and nBjets==0:
                    npass_wmu = npass_wmu + 1
                    print run,lumi,event
                    #print 'In Wmu CR : run', run, '  lumi', lumi, '  event', event
                    allquantities.reg_1mu1b_Wmass = myMuos_Wmass[iLeadLep]
                    allquantities.reg_1mu1b_WpT=WpT

                    allquantities.reg_1mu1b_hadrecoil = WmunuRecoil
                    allquantities.reg_1mu1b_MET = pfMet

                    allquantities.reg_1mu1b_lep1_pT=myMuos[iLeadLep].Pt()
                    allquantities.reg_1mu1b_lep1_iso=myMuIso[iLeadLep]

                    allquantities.reg_1mu1b_jet1_pT=j1.Pt()

                    if nJets>1: allquantities.reg_1mu1b_jet2_pT=j2.Pt()

                    allquantities.reg_1mu1b_jet1_eta=j1.Eta()
                    if nJets>1: allquantities.reg_1mu1b_jet2_eta=j2.Eta()

                    allquantities.reg_1mu1b_jet1_NHadEF=myJetNhadEF[ifirstjet]
                    allquantities.reg_1mu1b_jet1_CHadEF=myJetChadEF[ifirstjet]
                    allquantities.reg_1mu1b_jet1_CEmEF=myJetCEmEF[ifirstjet]
                    allquantities.reg_1mu1b_jet1_PhoEF=myJetPhoEF[ifirstjet]
                    allquantities.reg_1mu1b_jet1_EleEF=myJetEleEF[ifirstjet]
                    allquantities.reg_1mu1b_jet1_MuoEF=myJetMuoEF[ifirstjet]

                    allquantities.reg_1mu1b_njet = nJets

                    if options.CSV:
                        allquantities.reg_1mu1b_jet1_csv = myJetCSV[ifirstjet]
                        if nJets>1: allquantities.reg_1mu1b_jet2_csv = myJetCSV[isecondjet]

                    if options.DeepCSV:
                        allquantities.reg_1mu1b_jet1_deepcsv = myJetCSV[ifirstjet]
                        if nJets>1: allquantities.reg_1mu1b_jet2_deepcsv = myJetCSV[isecondjet]

                    allquantities.reg_1mu1b_min_dPhi_jet_Recoil = min( [DeltaPhi(WmunuPhi,myJetP4[nb].Phi()) for nb in range(nJets)] )
                    allquantities.reg_1mu1b_min_dPhi_jet_MET = min( [DeltaPhi(pfMetPhi,myJetP4[nb].Phi()) for nb in range(nJets)] )

                    allquantities.reg_1mu1b_ntau = nTauTightMuon
                    allquantities.reg_1mu1b_nele = nEle
                    allquantities.reg_1mu1b_nmu = nMu
                    allquantities.reg_1mu1b_nUncleanTau = nUncleanTau
                    allquantities.reg_1mu1b_nUncleanEle = nUncleanEle
                    allquantities.reg_1mu1b_nUncleanMu = nUncleanMu
                    allquantities.lep_syst_1mu1b_up = WmunuRecoil
                    allquantities.lep_syst_1mu1b_down = WmunuRecoil
                    allquantities.btag_syst_1mu1b_up = WmunuRecoil
                    allquantities.btag_syst_1mu1b_down = WmunuRecoil
                    allquantities.metTrig_syst_1mu1b_up = WmunuRecoil
                    allquantities.metTrig_syst_1mu1b_down = WmunuRecoil
                    allquantities.ewkZ_syst_1mu1b_up = WmunuRecoil
                    allquantities.ewkZ_syst_1mu1b_down = WmunuRecoil
                    allquantities.ewkW_syst_1mu1b_up = WmunuRecoil
                    allquantities.ewkW_syst_1mu1b_down = WmunuRecoil
                    allquantities.ewkTop_syst_1mu1b_up = WmunuRecoil
                    allquantities.ewkTop_syst_1mu1b_down = WmunuRecoil
                    allquantities.pho_syst_1mu1b_up = WmunuRecoil
                    allquantities.pho_syst_1mu1b_down = WmunuRecoil
                    allquantities.jec_syst_1mu1b_up = WmunuRecoil
                    allquantities.jec_syst_1mu1b_down = WmunuRecoil
                    allquantities.jer_syst_1mu1b_up = WmunuRecoil
                    allquantities.jer_syst_1mu1b_down = WmunuRecoil
                    isWmunuCR1 = True



# -------------------------------------------
# Top 1lep CR
# -------------------------------------------

        #W CR specific bools

        WdPhicond=True

        if applydPhicut:
            if WenuPhi>-10.:
                if min( [DeltaPhi(WenuPhi,myJetP4[nb].Phi()) for nb in range(nJets)] ) < 0.5: WdPhicond=False
            if WmunuPhi>-10.:
                if min( [DeltaPhi(WmunuPhi,myJetP4[nb].Phi()) for nb in range(nJets)] ) < 0.5: WdPhicond=False



# -------------------------------------------
# Old Top CR
# -------------------------------------------

        #Top CR specific bools

        TopdPhicond=True


# -------------------------------------------
# Gamma CR
# -------------------------------------------

        #Gamma CR specific bools

        GammaPhicond=True




        # QCD CR

        # 1b

        QCD1b_Cut1_nJets        =   SR1_Cut1_nJets
        QCD1b_Cut2_nBjets       =   SR1_Cut2_nBjets
        QCD1b_Cut3_trigstatus   =   SR1_Cut3_trigstatus
        QCD1b_Cut4_jet1         =   SR1_Cut4_jet1

        QCD1b_Cut5_jet2         =   SR1_Cut5_jet2
        QCD1b_Cut6_dPhi_jet_MET =   not SR1_Cut6_dPhi_jet_MET
        QCD1b_Cut7_nLep         =   SR1_Cut7_nLep
        QCD1b_Cut8_pfMET        =   SR1_Cut8_pfMET





        #2b

        QCD2b_Cut1_nJets        =   SR2_Cut1_nJets
        QCD2b_Cut2_nBjets       =   SR2_Cut2_nBjets
        QCD2b_Cut3_trigstatus   =   SR2_Cut3_trigstatus
        QCD2b_Cut4_jet1         =   SR2_Cut4_jet1
        QCD2b_Cut5_jet2         =   SR2_Cut5_jet2

        QCD2b_Cut6_jet3         =   SR2_Cut6_jet3
        QCD2b_Cut7_dPhi_jet_MET =   not SR2_Cut7_dPhi_jet_MET
        QCD2b_Cut8_nLep         =   SR2_Cut8_nLep
        QCD2b_Cut9_pfMET        =   SR2_Cut9_pfMET


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
#        nleptons_ = (len(myTaus) + len(myMuos) + len(myEles))

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
        if not isData :
            gn_pT = getGenPt(samplename, nGenPar, genParId, genMomParId, genParSt,genParP4)
            #print gn_pT
            if samplename=="TT":
                genpTReweighting = GenWeightProducer(samplename, nGenPar, genParId, genMomParId, genParSt,genParP4)
            else:
                genpTReweighting = kfactor.getEWKW(gn_pT)*kfactor.getQCDW(gn_pT)
        #genpTReweighting = 1.0
        #print 'genpTReweighting: ', genpTReweighting
        #----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## MET reweights
        #----------------------------------------------------------------------------------------------------------------------------------------------------------------
        metTrig_firstmethodReweight=1.0
        metTrig_secondmethodReweight=1.0
        metTrig_firstmethodReweight_up=1.0
        metTrig_firstmethodReweight_down=1.0


        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Muon reweight
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        #
        uni = random.uniform(0., 1.)
        muonTrig_SF = 1.0
        muonTrig_SF_systUP=1.0
        muonTrig_SF_systDOWN=1.0
        if nMu == 1:
            mupt = muP4[0].Pt()
            abeta = abs(muP4[0].Eta())
        if nMu == 2:
            if muP4[0].Pt() > muP4[1].Pt():
                leadmu=0
            else:
                leadmu=1
            mupt = muP4[leadmu].Pt()
            abeta = abs(muP4[leadmu].Eta())
        if nMu==1 or nMu==2:
            if uni < 0.54:
                xbin = muonTrigSFs_EfficienciesAndSF_RunBtoF.GetXaxis().FindBin(abeta)
                ybin = muonTrigSFs_EfficienciesAndSF_RunBtoF.GetYaxis().FindBin(mupt)
                muonTrig_SF *= muonTrigSFs_EfficienciesAndSF_RunBtoF.GetBinContent(xbin,ybin)
                muonTrig_SF_systUP *= muonTrigSFs_EfficienciesAndSF_RunBtoF.GetBinContent(xbin,ybin) + muonTrigSFs_EfficienciesAndSF_RunBtoF.GetBinErrorUp(xbin,ybin)
                muonTrig_SF_systDOWN *= muonTrigSFs_EfficienciesAndSF_RunBtoF.GetBinContent(xbin,ybin) - muonTrigSFs_EfficienciesAndSF_RunBtoF.GetBinErrorLow(xbin,ybin)
            elif uni > 0.54:
                xbin = muonTrigSFs_EfficienciesAndSF_Period4.GetXaxis().FindBin(abeta)
                ybin = muonTrigSFs_EfficienciesAndSF_Period4.GetYaxis().FindBin(mupt)
                muonTrig_SF *= muonTrigSFs_EfficienciesAndSF_Period4.GetBinContent(xbin,ybin)
                muonTrig_SF_systUP *= muonTrigSFs_EfficienciesAndSF_Period4.GetBinContent(xbin,ybin) + muonTrigSFs_EfficienciesAndSF_Period4.GetBinErrorUp(xbin,ybin)
                muonTrig_SF_systDOWN *= muonTrigSFs_EfficienciesAndSF_Period4.GetBinContent(xbin,ybin) - muonTrigSFs_EfficienciesAndSF_Period4.GetBinErrorLow(xbin,ybin)
        muIDSF_loose = 1.0
        muIDSF_loose_systUP=1.0
        muIDSF_loose_systDOWN=1.0
        muIDSF_tight = 1.0
        muIDSF_tight_systUP=1.0
        muIDSF_tight_systDOWN=1.0
        for imu in range(nMu):
            mupt = muP4[imu].Pt()
            abeta = abs(muP4[imu].Eta())
            if uni < 0.54:
                if mupt > 30:
                    xbin = muonTightIDSFs_EfficienciesAndSF_BCDEF.GetXaxis().FindBin(abeta)
                    ybin = muonTightIDSFs_EfficienciesAndSF_BCDEF.GetYaxis().FindBin(mupt)
                    muIDSF_tight *= muonTightIDSFs_EfficienciesAndSF_BCDEF.GetBinContent(xbin,ybin)
                    muIDSF_tight_systUP *= (muonTightIDSFs_EfficienciesAndSF_BCDEF.GetBinContent(xbin,ybin) + muonTightIDSFs_EfficienciesAndSF_BCDEF.GetBinErrorUp(xbin,ybin))
                    muIDSF_tight_systDOWN *= (muonTightIDSFs_EfficienciesAndSF_BCDEF.GetBinContent(xbin,ybin) - muonTightIDSFs_EfficienciesAndSF_BCDEF.GetBinErrorLow(xbin,ybin))
                else:
                    xbin = muonLooseIDSFs_EfficienciesAndSF_BCDEF.GetXaxis().FindBin(abeta)
                    ybin = muonLooseIDSFs_EfficienciesAndSF_BCDEF.GetYaxis().FindBin(mupt)
                    muIDSF_loose *= muonLooseIDSFs_EfficienciesAndSF_BCDEF.GetBinContent(xbin,ybin)
                    muIDSF_loose_systUP *= (muonLooseIDSFs_EfficienciesAndSF_BCDEF.GetBinContent(xbin,ybin) + muonLooseIDSFs_EfficienciesAndSF_BCDEF.GetBinErrorUp(xbin,ybin))
                    muIDSF_loose_systDOWN *= (muonLooseIDSFs_EfficienciesAndSF_BCDEF.GetBinContent(xbin,ybin) - muonLooseIDSFs_EfficienciesAndSF_BCDEF.GetBinErrorLow(xbin,ybin))
            if uni > 0.54:
                if mupt > 30:
                    xbin = muonTightIDSFs_EfficienciesAndSF_GH.GetXaxis().FindBin(abeta)
                    ybin = muonTightIDSFs_EfficienciesAndSF_GH.GetYaxis().FindBin(mupt)
                    muIDSF_tight *= muonTightIDSFs_EfficienciesAndSF_GH.GetBinContent(xbin,ybin)
                    muIDSF_tight_systUP *= (muonTightIDSFs_EfficienciesAndSF_GH.GetBinContent(xbin,ybin) - muonTightIDSFs_EfficienciesAndSF_GH.GetBinErrorUp(xbin,ybin))
                    muIDSF_tight_systDOWN *= (muonTightIDSFs_EfficienciesAndSF_GH.GetBinContent(xbin,ybin) + muonTightIDSFs_EfficienciesAndSF_GH.GetBinErrorLow(xbin,ybin))
                else:
                    xbin = muonLooseIDSFs_EfficienciesAndSF_GH.GetXaxis().FindBin(abeta)
                    ybin = muonLooseIDSFs_EfficienciesAndSF_GH.GetYaxis().FindBin(mupt)
                    muIDSF_loose *= muonLooseIDSFs_EfficienciesAndSF_GH.GetBinContent(xbin,ybin)
                    muIDSF_loose_systUP *= (muonLooseIDSFs_EfficienciesAndSF_GH.GetBinContent(xbin,ybin) + muonLooseIDSFs_EfficienciesAndSF_GH.GetBinErrorUp(xbin,ybin))
                    muIDSF_loose_systDOWN *= (muonLooseIDSFs_EfficienciesAndSF_GH.GetBinContent(xbin,ybin) - muonLooseIDSFs_EfficienciesAndSF_GH.GetBinErrorLow(xbin,ybin))

        muIsoSF_loose = 1.0
        muIsoSF_loose_systUP=1.0
        muIsoSF_loose_systDOWN=1.0
        muIsoSF_tight = 1.0
        muIsoSF_tight_systUP=1.0
        muIsoSF_tight_systDOWN=1.0
        for imu in range(nMu):
            mupt = muP4[imu].Pt()
            abeta = abs(muP4[imu].Eta())
            muiso = myMuIso[imu]
            if uni < 0.54:
                if muiso < 0.15:
                    xbin = muonTightIsoSFs_EfficienciesAndSF_BCDEF.GetXaxis().FindBin(abeta)
                    ybin = muonTightIsoSFs_EfficienciesAndSF_BCDEF.GetYaxis().FindBin(mupt)
                    muIsoSF_tight *= muonTightIsoSFs_EfficienciesAndSF_BCDEF.GetBinContent(xbin,ybin)
                    muIsoSF_tight_systUP *= (muonTightIsoSFs_EfficienciesAndSF_BCDEF.GetBinContent(xbin,ybin) + muonTightIsoSFs_EfficienciesAndSF_BCDEF.GetBinErrorUp(xbin,ybin))
                    muIsoSF_tight_systDOWN *= (muonTightIsoSFs_EfficienciesAndSF_BCDEF.GetBinContent(xbin,ybin) - muonTightIsoSFs_EfficienciesAndSF_BCDEF.GetBinErrorLow(xbin,ybin))
                elif muiso < 0.25:
                    xbin = muonLooseIsoSFs_EfficienciesAndSF_BCDEF.GetXaxis().FindBin(abeta)
                    ybin = muonLooseIsoSFs_EfficienciesAndSF_BCDEF.GetYaxis().FindBin(mupt)
                    muIsoSF_loose *= muonLooseIsoSFs_EfficienciesAndSF_BCDEF.GetBinContent(xbin,ybin)
                    muIsoSF_loose_systUP *= (muonLooseIsoSFs_EfficienciesAndSF_BCDEF.GetBinContent(xbin,ybin) + muonLooseIsoSFs_EfficienciesAndSF_BCDEF.GetBinErrorUp(xbin,ybin))
                    muIsoSF_loose_systDOWN *= (muonLooseIsoSFs_EfficienciesAndSF_BCDEF.GetBinContent(xbin,ybin) - muonLooseIsoSFs_EfficienciesAndSF_BCDEF.GetBinErrorLow(xbin,ybin))
            if uni > 0.54:
                if muiso < 0.15:
                    xbin = muonTightIsoSFs_EfficienciesAndSF_GH.GetXaxis().FindBin(abeta)
                    ybin = muonTightIsoSFs_EfficienciesAndSF_GH.GetYaxis().FindBin(mupt)
                    muIsoSF_tight *= muonTightIsoSFs_EfficienciesAndSF_GH.GetBinContent(xbin,ybin)
                    muIsoSF_tight_systUP *= (muonTightIsoSFs_EfficienciesAndSF_GH.GetBinContent(xbin,ybin) + muonTightIsoSFs_EfficienciesAndSF_GH.GetBinErrorUp(xbin,ybin))
                    muIsoSF_tight_systDOWN *= (muonTightIsoSFs_EfficienciesAndSF_GH.GetBinContent(xbin,ybin) - muonTightIsoSFs_EfficienciesAndSF_GH.GetBinErrorLow(xbin,ybin))
                elif muiso < 0.25:
                    xbin = muonLooseIsoSFs_EfficienciesAndSF_GH.GetXaxis().FindBin(abeta)
                    ybin = muonLooseIsoSFs_EfficienciesAndSF_GH.GetYaxis().FindBin(mupt)
                    muIsoSF_loose *= muonLooseIsoSFs_EfficienciesAndSF_GH.GetBinContent(xbin,ybin)
                    muIsoSF_loose_systUP *= (muonLooseIsoSFs_EfficienciesAndSF_GH.GetBinContent(xbin,ybin) + muonLooseIsoSFs_EfficienciesAndSF_GH.GetBinErrorUp(xbin,ybin))
                    muIsoSF_loose_systDOWN *= (muonLooseIsoSFs_EfficienciesAndSF_GH.GetBinContent(xbin,ybin) - muonLooseIsoSFs_EfficienciesAndSF_GH.GetBinErrorLow(xbin,ybin))

        muTracking_SF = 1.0
        muTracking_SF_systUP=1.0
        muTracking_SF_systDOWN=1.0
        for imu in range(nMu):
            abeta = abs(muP4[imu].Eta())
            muTracking_SF *= muonTrackingSFs_EfficienciesAndSF_BCDEFGH.Eval(abeta)
            ybin = muonTrackingSFs_EfficienciesAndSF_BCDEFGH.GetYaxis().FindBin(abeta)
            muTracking_SF_systUP *= (muonTrackingSFs_EfficienciesAndSF_BCDEFGH.Eval(abeta) + muonTrackingSFs_EfficienciesAndSF_BCDEFGH.GetErrorYhigh(ybin))
            muTracking_SF_systDOWN *= (muonTrackingSFs_EfficienciesAndSF_BCDEFGH.Eval(abeta) - muonTrackingSFs_EfficienciesAndSF_BCDEFGH.GetErrorYlow(ybin))


        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Electron reweight
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        eleTrig_reweight = 1.0
        eleTrig_reweight_systUP = 1.0
        eleTrig_reweight_systDOWN = 1.0
        if nEle == 1:
            elept = eleP4[0].Pt()
            eleeta = eleP4[0].Eta()
        if nEle == 2:
            if eleP4[0].Pt() > eleP4[1].Pt():
                leadele=0
            else:
                leadele=1
            elept = eleP4[leadele].Pt()
            eleeta = eleP4[leadele].Eta()
        if nEle==1 or nEle==2:
            xbin = eleTrig_hEffEtaPt.GetXaxis().FindBin(eleeta)
            ybin = eleTrig_hEffEtaPt.GetYaxis().FindBin(elept)
            eleTrig_reweight *= eleTrig_hEffEtaPt.GetBinContent(xbin,ybin)
            eleTrig_reweight_systUP *= (eleTrig_hEffEtaPt.GetBinContent(xbin,ybin) + eleTrig_hEffEtaPt.GetBinErrorUp(xbin,ybin))
            eleTrig_reweight_systDOWN *= (eleTrig_hEffEtaPt.GetBinContent(xbin,ybin) - eleTrig_hEffEtaPt.GetBinErrorLow(xbin,ybin))
#            print 'eleTrig_reweight_systUP, eleTrig_reweight, eleTrig_reweight_systDOWN', eleTrig_reweight_systUP, eleTrig_reweight, eleTrig_reweight_systDOWN


        eleRecoSF = 1.0
        eleRecoSF_systUP = 1.0
        eleRecoSF_systDOWN = 1.0
        for iele in range(nEle):
            elept = eleP4[iele].Pt()
            eleeta = eleP4[iele].Eta()
            xbin = eleRecoSF_EGamma_SF2D.GetXaxis().FindBin(eleeta)
            ybin = eleRecoSF_EGamma_SF2D.GetYaxis().FindBin(elept)
            eleRecoSF *= eleRecoSF_EGamma_SF2D.GetBinContent(xbin,ybin)
            eleRecoSF_systUP *= (eleRecoSF_EGamma_SF2D.GetBinContent(xbin,ybin) + eleRecoSF_EGamma_SF2D.GetBinErrorUp(xbin,ybin))
            eleRecoSF_systDOWN *= (eleRecoSF_EGamma_SF2D.GetBinContent(xbin,ybin) - eleRecoSF_EGamma_SF2D.GetBinErrorLow(xbin,ybin))
#            print 'eleRecoSF_systUP, eleRecoSF, eleRecoSF_systDOWN', eleRecoSF_systUP, eleRecoSF, eleRecoSF_systDOWN


        eleIDSF_loose = 1.0
        eleIDSF_loose_systUP = 1.0
        eleIDSF_loose_systDOWN = 1.0
        eleIDSF_tight = 1.0
        eleIDSF_tight_systUP = 1.0
        eleIDSF_tight_systDOWN = 1.0
        for iele in range(nEle):
            elept = eleP4[iele].Pt()
            eleeta = eleP4[iele].Eta()
            if elept > 30:
                xbin = eleTightIDSF_EGamma_SF2D.GetXaxis().FindBin(eleeta)
                ybin = eleTightIDSF_EGamma_SF2D.GetYaxis().FindBin(elept)
                eleIDSF_tight *= eleTightIDSF_EGamma_SF2D.GetBinContent(xbin,ybin)
                eleIDSF_tight_systUP *= (eleTightIDSF_EGamma_SF2D.GetBinContent(xbin,ybin) + eleTightIDSF_EGamma_SF2D.GetBinErrorUp(xbin,ybin))
                eleIDSF_tight_systDOWN *= (eleTightIDSF_EGamma_SF2D.GetBinContent(xbin,ybin) - eleTightIDSF_EGamma_SF2D.GetBinErrorLow(xbin,ybin))
            else:
                xbin = eleLooseIDSF_EGamma_SF2D.GetXaxis().FindBin(eleeta)
                ybin = eleLooseIDSF_EGamma_SF2D.GetYaxis().FindBin(elept)
                eleIDSF_loose *= eleLooseIDSF_EGamma_SF2D.GetBinContent(xbin,ybin)
                eleIDSF_loose_systUP *= (eleLooseIDSF_EGamma_SF2D.GetBinContent(xbin,ybin) + eleLooseIDSF_EGamma_SF2D.GetBinErrorUp(xbin,ybin))
                eleIDSF_loose_systDOWN *= (eleLooseIDSF_EGamma_SF2D.GetBinContent(xbin,ybin) - eleLooseIDSF_EGamma_SF2D.GetBinErrorLow(xbin,ybin))

        eleVetoCutBasedIDSF = 1.0
        for iele in range(nEle):
            elept = eleP4[iele].Pt()
            eleeta = eleP4[iele].Eta()
            xbin = eleVetoCutBasedIDSF_egammaEffi_txt_EGM2D.GetXaxis().FindBin(eleeta)
            ybin = eleVetoCutBasedIDSF_egammaEffi_txt_EGM2D.GetYaxis().FindBin(elept)
            eleVetoCutBasedIDSF *= eleVetoCutBasedIDSF_egammaEffi_txt_EGM2D.GetBinContent(xbin,ybin)

        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Photon reweight
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        phoIDSF_loose = 1.0
        phoIDSF_loose_systUP = 1.0
        phoIDSF_loose_systDOWN = 1.0
        phoIDSF_tight = 1.0
        phoIDSF_tight_systUP = 1.0
        phoIDSF_tight_systDOWN = 1.0
        for pho in range(nPho):
            phopt = phoP4[ipho].Pt()
            phoeta = phoP4[ipho].Eta()
            if phopt > 175:
                xbin = phoTightIDSF_EGamma_SF2D.GetXaxis().FindBin(phoeta)
                ybin = phoTightIDSF_EGamma_SF2D.GetYaxis().FindBin(phopt)
                phoIDSF_tight *= phoTightIDSF_EGamma_SF2D.GetBinContent(xbin,ybin)
                phoIDSF_tight_systUP *= (phoTightIDSF_EGamma_SF2D.GetBinContent(xbin,ybin) + phoTightIDSF_EGamma_SF2D.GetBinErrorUp(xbin,ybin))
                phoIDSF_tight_systDOWN *= (phoTightIDSF_EGamma_SF2D.GetBinContent(xbin,ybin) - phoTightIDSF_EGamma_SF2D.GetBinErrorLow(xbin,ybin))
            else:
                xbin = phoLooseIDSF_EGamma_SF2D.GetXaxis().FindBin(phoeta)
                ybin = phoLooseIDSF_EGamma_SF2D.GetYaxis().FindBin(phopt)
                phoIDSF_loose *= phoLooseIDSF_EGamma_SF2D.GetBinContent(xbin,ybin)
                phoIDSF_loose_systUP *= (phoLooseIDSF_EGamma_SF2D.GetBinContent(xbin,ybin) + phoLooseIDSF_EGamma_SF2D.GetBinErrorUp(xbin,ybin))
                phoIDSF_loose_systDOWN *= (phoLooseIDSF_EGamma_SF2D.GetBinContent(xbin,ybin) - phoLooseIDSF_EGamma_SF2D.GetBinErrorLow(xbin,ybin))

        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## JEC uncertainity
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        #print ('\n')
        jecUncUP = 1.0
        jecUncDOWN = 1.0
        jecUnc = 1.0
        for jet in range(nJets):
            jecUnc *= myJetCorrUnc[jet]
        jecUncUP = 1+jecUnc
        jecUncDOWN = 1-jecUnc
        #print jecUncUP,jecUncDOWN
        #print ('\n')
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Pileup weight
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
#        allpuweights = PUWeight()
#        len_puweight = len(allpuweights)
        puweight = 0.0
        if isData: puweight = 1.0
        if not isData:
#            if pu_nTrueInt  <= len_puweight: puweight = allpuweights[pu_nTrueInt-1]
#            if pu_nTrueInt  > len_puweight : puweight = 0.0
            if pu_nTrueInt < 100 and pu_nTrueInt > 0:
                puweight = pileup2016histo.GetBinContent(pu_nTrueInt)
            else:
                puweight = 1.


#            print pu_nTrueInt

        #print 'puweight:  ', puweight
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Total weight
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        if puweight == 0.0:
#            print 'Warning:: puweight is 0, setting it to 1'
            puweight = 1.0

        if genpTReweighting == 0.0:
#            print 'Warning:: genpTReweighting is 0, setting it to 1'
            genpTReweighting = 1.0
        #print 'genpTReweighting: ',genpTReweighting

        if metTrig_firstmethodReweight == 0.0:
#            print 'Warning:: metTrig_Reweight is 0, setting it to 1'
            metTrig_firstmethodReweight = 1.0

        if metTrig_firstmethodReweight_up == 0.0:
#            print 'Warning:: metTrig_Reweight is 0, setting it to 1'
            metTrig_firstmethodReweight_up = 1.0

        if metTrig_firstmethodReweight_down == 0.0:
#            print 'Warning:: metTrig_Reweight is 0, setting it to 1'
            metTrig_firstmethodReweight_down = 1.0

        muweights = muonTrig_SF * muIDSF_loose * muIDSF_tight * muIsoSF_loose * muIsoSF_tight * muTracking_SF
        if muweights == 0.0:
#            print 'Warning:: muon weight is 0, setting it to 1'
            muweights = 1.0
        #print 'muweights:  ', muweights

        #eleweights = eleTrig_reweight * eleRecoSF * eleIDSF_loose * eleIDSF_tight * eleVetoCutBasedIDSF
        eleweights = eleTrig_reweight * eleRecoSF * eleIDSF_loose * eleIDSF_tight
        if eleweights == 0.0:
#            print 'Warning:: electron weight is 0, setting it to 1'
            eleweights = 1.0
        #print 'eleweights: ',eleweights
        phoweights = phoIDSF_loose * phoIDSF_tight
        if phoweights == 0.0:
#            print 'Warning:: photon weight is 0, setting it to 1'
            phoweights = 1.0

        allweights = puweight * mcweight * genpTReweighting * eleweights * metTrig_firstmethodReweight * muweights*phoweights


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        temp_weight_withOutBtag = allweights
        ## BTag Scale Factor
        if SR1njetcond and False:
            ij = ifirstjet
            if nJets>1: jj = isecondjet

            flav1 = jetflav(myJetHadronFlavor[ij])
            if nJets>1: flav2 = jetflav(myJetHadronFlavor[jj])

#            print ("ij, flav, pt, eta, ",ij, flav1, myJetP4[ij].Pt(), myJetP4[ij].Eta())
            reader1.eval_auto_bounds('central', 0, 1.2, 50.)
            sf_resolved1 = weightbtag(reader1, flav1, myJetP4[ij].Pt(), myJetP4[ij].Eta())
            if nJets>1: sf_resolved2 = weightbtag(reader1, flav2, myJetP4[jj].Pt(), myJetP4[jj].Eta())

#            print (sf_resolved1, sf_resolved2)
        elif SR2njetcond and False:
            ij = ifirstjet
            jj = isecondjet
            if nJets>2: jk = ithirdjet

            flav1 = jetflav(myJetHadronFlavor[ij])
            flav2 = jetflav(myJetHadronFlavor[jj])
            if nJets>2: flav3 = jetflav(myJetHadronFlavor[jj])

#            print ("ij, flav, pt, eta, ",ij, flav1, myJetP4[ij].Pt(), myJetP4[ij].Eta())
            reader1.eval_auto_bounds('central', 0, 1.2, 50.)
            sf_resolved1 = weightbtag(reader1, flav1, myJetP4[ij].Pt(), myJetP4[ij].Eta())
            sf_resolved2 = weightbtag(reader1, flav2, myJetP4[jj].Pt(), myJetP4[jj].Eta())
            if nJets>2: sf_resolved3 = weightbtag(reader1, flav3, myJetP4[jk].Pt(), myJetP4[jk].Eta())


        if tDMJetCond:
            P_MC = 1.0
            P_Data = 1.0
            btagweight = 1.0
            for i in range(nJets):
                nontag_eff = getBeff(myJetP4[i],myJetHadronFlavor[i])
                P_MC *= (1-nontag_eff)
                reader1.eval_auto_bounds('central', 0, 2.4, 30.)
                SF_jet = []
                SF_jet=weightbtag(reader1, jetflav(myJetHadronFlavor[jet]), myJetP4[jet].Pt(), myJetP4[jet].Eta())
                P_Data *= (1 - SF_jet[0] * nontag_eff)
                #print 'eff: ',getBeff(myJetP4[i])
            if P_MC > 0:
                btagweight = P_Data/P_MC
            allweights = allweights * btagweight

        #print 'btag_weight:  ', allweights/temp_weight_withOutBtag
        ## for not applying btag weights
        #allweights = temp_weight_withOutBtag
        #print 'allweights without genpTReweighting: ', (allweights/genpTReweighting)
        #print 'allweights with genpTReweighting: ', (allweights)
        temp_original_weight  = allweights
        #print 'All weight step1:  ', allweights
        allweights_jec_up = temp_original_weight*jecUncUP
        allquantities.weight_jec_up = allweights_jec_up
        allweights_jec_down= temp_original_weight*jecUncDOWN
        allquantities.weight_jec_down = allweights_jec_down
        allweights_ewkW_down = temp_original_weight
        allweights_ewkW_up = temp_original_weight
        allweights_ewkZ_down = temp_original_weight
        allweights_ewkZ_up = temp_original_weight
        allweights_ewkTop_down = temp_original_weight
        allweights_ewkTop_up = temp_original_weight
        allweights_metTrig_up = temp_original_weight
        allweights_metTrig_down = temp_original_weight

        allweights_metTrig_up = (allweights/metTrig_firstmethodReweight)*metTrig_firstmethodReweight_up
        allweights_metTrig_down = (allweights/metTrig_firstmethodReweight)*metTrig_firstmethodReweight_down
        temp_weight_withBtag = allweights/(eleweights*muweights)
        temp_weight_withBtag_noPhoSF = allweights/phoweights

        if isData: allweights = 1.0
        allweights_noPU = allweights/puweight


#----------------------------------------------------------------------------------------------------------------------------------------------------------------

        if samplename=="WJETS":
            allweights_ewkW_down = temp_original_weight/genpTReweighting
            allweights_ewkW_up = temp_original_weight*genpTReweighting
            allquantities.weight_ewkW_up  = allweights_ewkW_up
            allquantities.weight_ewkW_down  =  allweights_ewkW_down

        if samplename == "ZJETS":
            allweights_ewkZ_down = temp_original_weight/genpTReweighting
            allweights_ewkZ_up = temp_original_weight*genpTReweighting
            allquantities.weight_ewkZ_up  = allweights_ewkZ_up
            allquantities.weight_ewkZ_down  =  allweights_ewkZ_down
        if samplename == "TT":
            allweights_ewkTop_down = temp_original_weight*0.5
            allweights_ewkTop_up = temp_original_weight*1.5
            allquantities.weight_ewkTop_up  = allweights_ewkTop_up
            allquantities.weight_ewkTop_down  =  allweights_ewkTop_down
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

        if EleCRtrigstatus:
            CR1e1bCutFlow['trigger']+=allweights
            if nEle_tight==1 and nMu==0 and myEles_tight[0].Pt() > 30.:# and myEleTightID[0]:
                CR1e1bCutFlow['leptons=1']+=allweights
                if nJets>=3:
                    CR1e1bCutFlow['njets>=3']+=allweights
                    if nBjets==0:
                        CR1e1bCutFlow['nbjets=0']+=allweights
                        #if pfMet >= 250.:
                            #CR1e1bCutFlow['MET>250']+=allweights
                        if Wenumass <= 160.:
                            CR1e1bCutFlow['mT<160']+=allweights

  #Cutflow
        if MuCRtrigstatus:
            CR1mu1bCutFlow['trigger']+=allweights
            if nEle==0 and nMu_tight==1 and myMuos_tight[0].Pt() > 30.:# and myMuTightID[0] and myMuIso[0]<0.15:
                CR1mu1bCutFlow['leptons=1']+=allweights
                if nJets>=3:
                    CR1mu1bCutFlow['njets>=3']+=allweights
                    if nBjets==0:
                        CR1mu1bCutFlow['nbjets=0']+=allweights
                        #if pfMet >= 250.:
                            #CR1mu1bCutFlow['MET>250']+=allweights
                        if Wmunumass <= 160.:
                            CR1mu1bCutFlow['mT<160']+=allweights


 #Cutflow

 # ---
        #CR Summary
        if isZeeCR1:
            CRSummary['2e1b']+=allweights
            CRSummaryEle['2e1b']+=allweights
        if isZeeCR2:
            CRSummary['2e2b']+=allweights
            CRSummaryEle['2e2b']+=allweights
        if isZmumuCR1:
            CRSummary['2#mu1b']+=allweights
            CRSummaryMu['2#mu1b']+=allweights
        if isZmumuCR2:
            CRSummary['2#mu2b']+=allweights
            CRSummaryMu['2#mu2b']+=allweights

        if isWenuCR1:
            CRSummary['1e1b']+=allweights
            CRSummaryEle['1e1b']+=allweights
        if isWenuCR2:
            CRSummary['1e2b']+=allweights
            CRSummaryEle['1e2b']+=allweights
        if isWmunuCR1:
            CRSummary['1#mu1b']+=allweights
            CRSummaryMu['1#mu1b']+=allweights
        if isWmunuCR2:
            CRSummary['1#mu2b']+=allweights
            CRSummaryMu['1#mu2b']+=allweights

        if isTopenuCR1:
            CRSummary['1etop1b']+=allweights
            CRSummaryEle['1etop1b']+=allweights
        if isTopenuCR2:
            CRSummary['1etop2b']+=allweights
            CRSummaryEle['1etop2b']+=allweights
        if isTopmunuCR1:
            CRSummary['1#mutop1b']+=allweights
            CRSummaryMu['1#mutop1b']+=allweights
        if isTopmunuCR2:
            CRSummary['1#mutop2b']+=allweights
            CRSummaryMu['1#mutop2b']+=allweights

        if isTopCR1:
            CRSummary['1#mu1e1b']+=allweights
            CRSummaryMu['1#mu1e1b']+=allweights
        if isTopCR2:
            CRSummary['1#mu1e2b']+=allweights
            CRSummaryMu['1#mu1e2b']+=allweights


 #--------------------------------------------------------------------------------------------------------------------------------------------------------------------
        allquantities.met             = pfMet
        allquantities.N_e             = nEle
        allquantities.N_mu            = nMu
        allquantities.N_tau           = nTauLooseEleMu
        allquantities.N_Pho           = nPho
        allquantities.N_b             = nBjets
        allquantities.N_j             = nJets
        allquantities.weight          = allweights
        #print ' All weight Step 3:  ', allweights
        allquantities.weight_NoPU     = allweights_noPU
        allquantities.weight_met_up   = allweights_metTrig_up
        allquantities.weight_met_down = allweights_metTrig_down
        allquantities.totalevents     = 1

        btag_sysnum=0
        for btag_sysnum in[1,2]:
            allweights = temp_weight_withOutBtag
            if SR1njetcond and False:
                if sf_resolved1[btag_sysnum]==0.0: sf_resolved1[btag_sysnum]=1.0
                allweights = allweights*sf_resolved1[btag_sysnum]
                if nJets>1:
                    if sf_resolved2[btag_sysnum]==0.0: sf_resolved2[btag_sysnum]=1.0
                    allweights = allweights *sf_resolved2[btag_sysnum]
            if SR2njetcond and False:
                if sf_resolved1[btag_sysnum]==0.0: sf_resolved1[btag_sysnum]=1.0
                if sf_resolved2[btag_sysnum]==0.0: sf_resolved2[btag_sysnum]=1.0
                allweights = allweights * sf_resolved1[btag_sysnum] * sf_resolved2[btag_sysnum]
                if nJets>2:
                    if sf_resolved3[btag_sysnum]==0.0: sf_resolved3[btag_sysnum]=1.0
                    allweights = allweights * sf_resolved3[btag_sysnum]
#                print 'btag central weight', temp_original_weight
                if btag_sysnum==2:
                    allquantities.weight_btag_up = allweights
#                    print 'btag up weight', allweights
                if btag_sysnum==1:
                    allquantities.weight_btag_down = allweights
#                    print 'btag down weight', allweights
        allweights = temp_weight_withBtag
        muweights_systUP = muonTrig_SF_systUP * muIDSF_loose_systUP * muIDSF_tight_systUP * muIsoSF_loose_systUP * muIsoSF_tight_systUP * muTracking_SF_systUP
        muweights_systDOWN = muonTrig_SF_systDOWN * muIDSF_loose_systDOWN * muIDSF_tight_systDOWN * muIsoSF_loose_systDOWN * muIsoSF_tight_systDOWN * muTracking_SF_systDOWN
        eleweights_systUP = eleTrig_reweight_systUP * eleRecoSF_systUP * eleIDSF_loose_systUP * eleIDSF_tight_systUP
        eleweights_systDOWN = eleTrig_reweight_systDOWN * eleRecoSF_systDOWN * eleIDSF_loose_systDOWN * eleIDSF_tight_systDOWN
        if muweights_systUP == 0.0:
            muweights_systUP = 1.0
        if muweights_systDOWN == 0.0:
            muweights_systDOWN = 1.0
        if eleweights_systUP == 0.0:
            eleweights_systUP = 1.0
        if eleweights_systDOWN == 0.0:
            eleweights_systDOWN = 1.0
#                print 'muweights_systUP, eleweights_systUP',muweights_systUP ,eleweights_systUP
#                print 'muweights_systDOWN, eleweights_systDOWN',muweights_systDOWN,eleweights_systDOWN
        allweights = allweights * muweights_systUP * eleweights_systUP
        allquantities.weight_lep_up = allweights
#        if abs(allweights-temp_weight_withBtag) > 0.001:
#            print 'lep up value', allweights
        allweights = temp_weight_withBtag
#        if abs(allweights-temp_weight_withBtag) > 0.001:
#            print 'lep central weight', allweights
        allweights = allweights * muweights_systDOWN * eleweights_systDOWN
        allquantities.weight_lep_down = allweights
#        if abs(allweights-temp_weight_withBtag) > 0.001:
#            print 'lep down value', allweights
        allweights = temp_weight_withBtag_noPhoSF
        phoweights_systUP = phoIDSF_loose_systUP * phoIDSF_tight_systUP
        allquantities.weight_pho_up = allweights*phoweights_systUP
        phoweights_systDOWN = phoIDSF_loose_systDOWN * phoIDSF_tight_systDOWN
        allquantities.weight_pho_down = allweights*phoweights_systDOWN
        nPV = myJetNPV

        allquantities.PuReweightPV = nPV
        allquantities.noPuReweightPV = nPV

        if nMu==2 and nEle==0 and MuCRtrigstatus:
            allquantities.mu_PuReweightPV = nPV
            allquantities.mu_noPuReweightPV = nPV

        if nEle==2 and nMu==0 and EleCRtrigstatus:
            allquantities.ele_PuReweightPV = nPV
            allquantities.ele_noPuReweightPV = nPV

        if nPho==1 and nEle==0 and nMu==0 and PhotonCRtrigstatus:
            allquantities.pho_PuReweightPV = nPV
            allquantities.pho_noPuReweightPV = nPV

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


    #cutStatusSR1['pfmet'] = cutStatus['pfmet']
    #cutStatusSR2['pfmet'] = cutStatus['pfmet']
    print "Total events =", int(NEntries_total)
    print "Preselected events=", cutStatus['preselection']
    print "Selected events =", npass
    print "Selected events in We =", npass_we
    print "Selected events in Wmu =", npass_wmu
    # Cutflow
    cutflowTable=""
    cutflowHeader=""
    cutflowvalues=[]
    #cutflownames=['total','preselection','pfmet','njet+nBjet','jet1','jet2/3','lep']
    for cutflowname in cutflownames:
        cutflowvalues.append(cutStatus[cutflowname])
        cutflowTable += str(cutStatus[cutflowname])+" "
        cutflowHeader += cutflowname+" "

    cutflowvaluesSR1=[]
    #cutflownamesSR1=['total','preselection','pfmet','njet+nBjet','jet1','jet2','lep']
    for cutflowname in cutflownamesSR1:
        cutflowvaluesSR1.append(cutStatusSR1[cutflowname])

    cutflowvaluesSR2=[]
    #cutflownamesSR2=['total','preselection','pfmet','njet+nBjet','jet1','jet2','jet3','lep']
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
    #print
    #print "CRs:"
    #print CRStatus
#    print
#    print "CR Cutflow:"
#    print CRCutFlow
    print

    CRcutflowvaluesSet=[]
    CRcutnames=['total','preselection']+CRcutnames
    for CRreg in regionnames:
        CFvalues=[]
        for cutname in CRcutnames:
            #if ('top' in CRreg or '1e1b' in CRreg or '1e2b' in CRreg or '1mu1b' in CRreg or '1mu2b' in CRreg) and (not '1mu1e' in CRreg ) and 'njets' in cutname:
            #    cutname = 'add_Jet'
            exec("CFvalues.append(CR"+CRreg+"CutFlow['"+cutname+"'])")
        CRcutflowvaluesSet.append(CFvalues)

    print CRSummary

    allquantities.WriteHisto((NEntries_total,NEntries_Weight,npass,cutflowvalues,cutflownames,cutflowvaluesSR1,cutflownamesSR1, cutflowvaluesSR2,cutflownamesSR2,CRvalues,CRs,regionnames,CRcutnames,CRcutflowvaluesSet, CRSummary,regNames, CRSummaryMu,regNamesMu, CRSummaryEle,regNamesEle))

    if NEntries > 0 and NEntries_total > 0:
        eff=round(float(npass/float(NEntries_total)),5)
        print "efficiency =", eff
    else:
        eff = "NA"
        print "efficiency =", eff

#    os.system("mkdir -p "+outputdir+'/efficiencyfiles/')

#    f = open(outputdir+'/efficiencyfiles/'+textfile, 'w')
#    f.write(str(eff)+"\n\n#Cutflow Table:\n"+cutflowHeader[:-1]+"\n"+cutflowTable[:-1]+"\n\n#CR Table:\n"+CRHeader[:-1]+"\n"+CRTable[:-1])
    print "ROOT file written to", outfilename
#    print "Log written to "+outputdir+'/efficiencyfiles/'+textfile
#    f.close()
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

def DeltaR(P4_1,P4_2):
    return math.sqrt(  (  P4_1.Eta()-P4_2.Eta() )**2  + (  DeltaPhi(P4_1.Phi(),P4_2.Phi()) )**2 )


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

def getGenPt(sample,nGenPar, genParId, genMomParId, genParSt,genParP4):
    pt=0.0
    #################
    # WJets
    #################
    if sample=="WJETS":
    #if True:
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

    #################
    #ZJets
    #################
    if sample == "ZJETS":
    #if True:
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

    return pt

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
            w1 = TMath.Exp(0.0615 - 0.0005*l4_thisLep.Pt());
            w2 = TMath.Exp(0.0615 - 0.0005*l4_thatLep.Pt());
            k2 =  TMath.Sqrt(w1*w2)

    if(sample=="all"):
        k2 = 1.0

    return k2

def getBeff(P4):
    if flav == 5:
        xbin = eff_b_pass.GetXaxis().FindBin(P4.Eta())
        ybin = eff_b_pass.GetYaxis().FindBin(P4.Pt())
        btag_eff = eff_b_fail.GetBinContent(xbin,ybin)
        return btag_eff
    elif flav == 4:
        xbin = eff_c_pass.GetXaxis().FindBin(P4.Eta())
        ybin = eff_c_pass.GetYaxis().FindBin(P4.Pt())
        ctag_eff = eff_c_fail.GetBinContent(xbin,ybin)
        return ctag_eff
    elif flav!=4 and flav!=5:
        xbin = eff_light_pass.GetXaxis().FindBin(P4.Eta())
        ybin = eff_light_pass.GetYaxis().FindBin(P4.Pt())
        lighttag_eff = eff_light_fail.GetBinContent(xbin,ybin)
        return lighttag_eff

def MT(Pt, met, dphi):
    return ROOT.TMath.Sqrt( 2 * Pt * met * (1.0 - ROOT.TMath.Cos(dphi)) )

if __name__ == "__main__":
    ## analyze the tree and make histograms and all the 2D plots and Efficiency plots.
    if options.analyze:
        print "now calling analyzedataset"
        AnalyzeDataSet()
