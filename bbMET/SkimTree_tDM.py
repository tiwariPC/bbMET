#!/usr/bin/env python
from ROOT import TFile, TTree, TH1F, TH1D, TH1, TCanvas, TChain,TGraphAsymmErrors, TMath, TH2D, TLorentzVector, AddressOf, gROOT, TNamed
import ROOT as ROOT
import os
import sys, optparse
from array import array
import math
import numpy as numpy_

ROOT.gROOT.LoadMacro("Loader.h+")

outfilename= 'SkimmedTree.root'
PUPPI = True
CA15  = False
usage = "usage: %prog [options] arg1 arg2"
parser = optparse.OptionParser(usage)

parser.add_option("-i", "--inputfile",  dest="inputfile")
parser.add_option("-F", "--farmout", action="store_true",  dest="farmout")

(options, args) = parser.parse_args()

if options.farmout==None:
    isfarmout = False
else:
    isfarmout = options.farmout
inputfilename = options.inputfile

skimmedTree = TChain("tree/treeMaker")

if isfarmout:
    infile = open(inputfilename)
    failcount=0
    for ifile in infile:
        try:
            f_tmp = TFile.Open(ifile.rstrip(),'READ')
            print f_tmp
            if f_tmp.IsZombie():            # or fileIsCorr(ifile.rstrip()):
                failcount += 1
                continue
            skimmedTree.Add(ifile.rstrip())
        except:
            failcount += 1
    if failcount>0: print "Could not read %d files. Skipping them." %failcount
if not isfarmout:
    skimmedTree.Add(inputfilename)
        #    samplename = WhichSample(inputfilename)

#skimmedTree = TChain("tree/treeMaker")

#skimmedTree.Add(sys.argv[1])

def arctan(x,y):
    corr=0
    if (x>0 and y>=0) or (x>0 and y<0):
        corr=0
    elif x<0 and y>=0:
        corr=math.pi
    elif x<0 and y<0:
        corr=-math.pi
    if x!=0.:
        return math.atan(y/x)+corr
    else:
        return math.pi/2+corr

def getPT(P4):
    return P4.Pt()

def AnalyzeDataSet():
    CSVMWP=0.8484
    DCSVMWP=0.6324
    NEntries = skimmedTree.GetEntries()
#    NEntries = 1000
    h_total = TH1F('h_total','h_total',2,0,2)
    h_total_mcweight = TH1F('h_total_mcweight','h_total_mcweight',2,0,2)

    triglist=["HLT_PFMET170_BeamHaloCleaned_v","HLT_PFMET170_HBHE_BeamHaloCleaned_v","HLT_PFMET170_NotCleaned_v","HLT_PFMET170_NoiseCleaned_v","HLT_PFMET170_JetIdCleaned_v","HLT_PFMET170_HBHECleaned_v","HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v","HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v","HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v","HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v","HLT_PFMET110_PFMHT110_IDTight_v","HLT_IsoMu24_v","HLT_IsoTkMu24_v","HLT_IsoMu27_v","HLT_IsoTkMu27_v","HLT_Ele27_WPTight_Gsf_v","HLT_Ele105_CaloIdVT_GsfTrkIdT_v","HLT_Ele115_CaloIdVT_GsfTrkIdT_v","HLT_Ele32_WPTight_Gsf_v","HLT_IsoMu20_v","HLT_Ele27_eta2p1_WPTight_Gsf_v","HLT_Ele27_WPLoose_Gsf_v","HLT_Ele32_eta2p1_WPTight_Gsf_v","HLT_Photon165_HE10_v","HLT_Photon175_v","HLT_Ele105_CaloIdVT_GsfTrkIdT_v"]


    outfile = TFile(outfilename,'RECREATE')

    outTree = TTree( 'outTree', 'tree branches' )
    print str(f_tmp)
    if isfarmout:
        samplepath = TNamed('samplepath', str(f_tmp).split('"')[1])
    else:
        samplepath = TNamed('samplepath', str(inputfilename))

    st_runId                  = numpy_.zeros(1, dtype=int)
    st_lumiSection            = array( 'L', [ 0 ] )
    st_eventId                = array( 'L', [ 0 ] )
    st_pfMetCorrPt            = array( 'f', [ 0. ] )
    st_pfMetCorrPhi           = array( 'f', [ 0. ] )
    st_pfMetUncJetResUp       = ROOT.std.vector('float')()
    st_pfMetUncJetResDown     = ROOT.std.vector('float')()
    st_pfMetUncJetEnUp        = ROOT.std.vector('float')()
    st_pfMetUncJetEnDown      = ROOT.std.vector('float')()
    st_isData           = array( 'b', [ 0 ] )

    for trigs in triglist:
        exec("st_"+trigs+"  = array( 'b', [ 0 ] )")

    maxn = 10

    st_THINnJet                     = array( 'L', [ 0 ] ) #ROOT.std.vector('int')()
    st_THINjetP4                    = ROOT.std.vector('TLorentzVector')()
    st_THINjetCISVV2                = ROOT.std.vector('float')()
    st_THINjetHadronFlavor          = ROOT.std.vector('int')()
    st_THINjetNHadEF                = ROOT.std.vector('float')()
    st_THINjetCHadEF                = ROOT.std.vector('float')()

    st_THINjetCEmEF                 = ROOT.std.vector('float')()
    st_THINjetPhoEF                 = ROOT.std.vector('float')()
    st_THINjetEleEF                 = ROOT.std.vector('float')()
    st_THINjetMuoEF                 = ROOT.std.vector('float')()
    st_THINjetCorrUnc               = ROOT.std.vector('float')()

    st_AK4deepCSVnJet               = array( 'L', [ 0 ] ) #ROOT.std.vector('int')()
    st_AK4deepCSVjetP4              = ROOT.std.vector('TLorentzVector')()
    st_AK4deepCSVjetDeepCSV_b       = ROOT.std.vector('float')()
    st_AK4deepCSVjetHadronFlavor    = ROOT.std.vector('int')()
    st_AK4deepCSVjetNHadEF          = ROOT.std.vector('float')()
    st_AK4deepCSVjetCHadEF          = ROOT.std.vector('float')()
    st_AK4deepCSVjetCorrUnc         = ROOT.std.vector('float')()


    st_nEle                = array( 'L', [ 0 ] ) #ROOT.std.vector('int')()
    st_eleP4               = ROOT.std.vector('TLorentzVector')()
    st_eleIsPassTight      = ROOT.std.vector('bool')()

    st_nPho                = array( 'L', [ 0 ] ) #ROOT.std.vector('int')()
    st_phoP4               = ROOT.std.vector('TLorentzVector')()
    st_phoIsPassTight      = ROOT.std.vector('bool')()

    st_nMu= array( 'L', [ 0 ] ) #ROOT.std.vector('int')()
    st_muP4                = ROOT.std.vector('TLorentzVector')()
    st_isTightMuon         = ROOT.std.vector('bool')()
    st_muIso            = ROOT.std.vector('float')()


    st_HPSTau_n= array( 'L', [ 0 ] ) #ROOT.std.vector('int')()
    st_nTauTightElectron= array( 'L', [ 0 ] )
    st_nTauTightMuon= array( 'L', [ 0 ] )
    st_nTauTightEleMu= array( 'L', [ 0 ] )
    st_nTauLooseEleMu= array( 'L', [ 0 ] )


    mcweight = array( 'f', [ 0 ] )
    st_pu_nTrueInt= array( 'f', [ 0 ] ) #ROOT.std.vector('std::vector<float>')()
    st_pu_nPUVert= array( 'f', [ 0 ] )
    st_THINjetNPV= array( 'f', [ 0 ] ) #ROOT.std.vector('std::vector<float>')()
    st_AK4deepCSVjetNPV= array( 'f', [ 0 ] )

    st_nGenPar = array( 'L', [ 0 ] )
    st_genParId = ROOT.std.vector('int')()
    st_genMomParId = ROOT.std.vector('int')()
    st_genParSt = ROOT.std.vector('int')()
    st_genParP4 = ROOT.std.vector('TLorentzVector')()

    WenuRecoil = array( 'f', [ 0. ] )
    Wenumass = array( 'f', [ 0. ] )
    WenuPhi = array( 'f', [ 0. ] )

    WmunuRecoil = array( 'f', [ 0. ] )
    Wmunumass = array( 'f', [ 0. ] )
    WmunuPhi = array( 'f', [ 0. ] )

    ZeeRecoil = array( 'f', [ 0. ] )
    ZeeMass = array( 'f', [ 0. ] )
    ZeePhi = array( 'f', [ 0. ] )

    ZmumuRecoil = array( 'f', [ 0. ] )
    ZmumuMass = array( 'f', [ 0. ] )
    ZmumuPhi = array( 'f', [ 0. ] )

    TOPRecoil = array( 'f', [ 0. ] )
    TOPPhi = array( 'f', [ 0. ] )

    GammaRecoil = array('f',[0.])
    GammaPhi = array( 'f', [ 0. ] )

    outTree.Branch( 'st_runId', st_runId , 'st_runId/L')
    outTree.Branch( 'st_lumiSection', st_lumiSection , 'st_lumiSection/L')
    outTree.Branch( 'st_eventId',  st_eventId, 'st_eventId/L')
    outTree.Branch( 'st_pfMetCorrPt', st_pfMetCorrPt , 'st_pfMetCorrPt/F')
    outTree.Branch( 'st_pfMetCorrPhi', st_pfMetCorrPhi , 'st_pfMetCorrPhi/F')
    outTree.Branch( 'st_pfMetUncJetResUp', st_pfMetUncJetResUp)
    outTree.Branch( 'st_pfMetUncJetResDown', st_pfMetUncJetResDown)
    outTree.Branch( 'st_pfMetUncJetEnUp', st_pfMetUncJetEnUp )
    outTree.Branch( 'st_pfMetUncJetEnDown', st_pfMetUncJetEnDown)
    outTree.Branch( 'st_isData', st_isData , 'st_isData/O')

    for trigs in triglist:
        exec("outTree.Branch( 'st_"+trigs+"', st_"+trigs+" , 'st_"+trigs+"/O')")

#    outTree.Branch( 'st_MET_trig', st_MET_trig , 'st_MET_trig/O')
#    outTree.Branch( 'st_SE_trig', st_SE_trig , 'st_SE_trig/O')
#    outTree.Branch( 'st_SP_trig', st_SP_trig , 'st_SP_trig/O')
#    outTree.Branch( 'st_HLT_IsoMu20', st_HLT_IsoMu20 , 'st_HLT_IsoMu20/O')
#    outTree.Branch( 'st_HLT_Ele27_WPLoose_Gsf', st_HLT_Ele27_WPLoose_Gsf , 'st_HLT_Ele27_WPLoose_Gsf/O')

    outTree.Branch( 'st_THINnJet',st_THINnJet, 'st_THINnJet/L' )
    outTree.Branch( 'st_THINjetP4',st_THINjetP4 )
    outTree.Branch( 'st_THINjetCISVV2',st_THINjetCISVV2 )
    outTree.Branch( 'st_THINjetHadronFlavor',st_THINjetHadronFlavor )
    outTree.Branch( 'st_THINjetNHadEF',st_THINjetNHadEF )
    outTree.Branch( 'st_THINjetCHadEF',st_THINjetCHadEF )

    outTree.Branch( 'st_THINjetCEmEF',st_THINjetCEmEF )
    outTree.Branch( 'st_THINjetPhoEF',st_THINjetPhoEF )
    outTree.Branch( 'st_THINjetEleEF',st_THINjetEleEF )
    outTree.Branch( 'st_THINjetMuoEF',st_THINjetMuoEF )
    outTree.Branch('st_THINjetCorrUnc', st_THINjetCorrUnc)

    outTree.Branch( 'st_AK4deepCSVnJet',st_AK4deepCSVnJet, 'st_AK4deepCSVnJet/L' )
    outTree.Branch( 'st_AK4deepCSVjetP4',st_AK4deepCSVjetP4 )
    outTree.Branch( 'st_AK4deepCSVjetDeepCSV_b',st_AK4deepCSVjetDeepCSV_b )
    outTree.Branch( 'st_AK4deepCSVjetHadronFlavor',st_AK4deepCSVjetHadronFlavor )
    outTree.Branch( 'st_AK4deepCSVjetNHadEF',st_AK4deepCSVjetNHadEF )
    outTree.Branch( 'st_AK4deepCSVjetCHadEF',st_AK4deepCSVjetCHadEF )
    outTree.Branch( 'st_AK4deepCSVjetCorrUnc', st_AK4deepCSVjetCorrUnc)

    outTree.Branch( 'st_nEle',st_nEle , 'st_nEle/L')
    outTree.Branch( 'st_eleP4',st_eleP4 )
    outTree.Branch( 'st_eleIsPassTight', st_eleIsPassTight)#, 'st_eleIsPassTight/O' )

    outTree.Branch( 'st_nPho',st_nPho , 'st_nPho/L')
    outTree.Branch( 'st_phoP4',st_phoP4 )
    outTree.Branch( 'st_phoIsPassTight', st_phoIsPassTight)#, 'st_phoIsPassTight/O' )


    outTree.Branch( 'st_nMu',st_nMu , 'st_nMu/L')
    outTree.Branch( 'st_muP4',st_muP4 )
    outTree.Branch( 'st_isTightMuon', st_isTightMuon)#, 'st_isTightMuon/O' )
    outTree.Branch( 'st_muIso', st_muIso)#, 'st_muIso/F')

    outTree.Branch( 'st_HPSTau_n', st_HPSTau_n, 'st_HPSTau_n/L')
    outTree.Branch( 'st_nTauTightElectron', st_nTauTightElectron, 'st_nTauTightElectron/L')
    outTree.Branch( 'st_nTauTightMuon', st_nTauTightMuon, 'st_nTauTightMuon/L')
    outTree.Branch( 'st_nTauTightEleMu', st_nTauTightEleMu, 'st_nTauTightEleMu/L')
    outTree.Branch( 'st_nTauLooseEleMu', st_nTauLooseEleMu, 'st_nTauLooseEleMu/L')

    outTree.Branch( 'st_pu_nTrueInt', st_pu_nTrueInt, 'st_pu_nTrueInt/F')
    outTree.Branch( 'st_pu_nPUVert', st_pu_nPUVert, 'st_pu_nPUVert/F')
    outTree.Branch( 'st_AK4deepCSVjetNPV', st_AK4deepCSVjetNPV, 'st_AK4deepCSVjetNPV/F')
    outTree.Branch( 'st_THINjetNPV', st_THINjetNPV, 'st_THINjetNPV/F')
    outTree.Branch( 'mcweight', mcweight, 'mcweight/F')
    outTree.Branch( 'st_nGenPar',st_nGenPar,'st_nGenPar/L' )  #nGenPar/I
    outTree.Branch( 'st_genParId',st_genParId )  #vector<int>
    outTree.Branch( 'st_genMomParId',st_genMomParId )
    outTree.Branch( 'st_genParSt',st_genParSt )
    outTree.Branch( 'st_genParP4', st_genParP4)

    outTree.Branch( 'WenuRecoil', WenuRecoil, 'WenuRecoil/F')
    outTree.Branch( 'Wenumass', Wenumass, 'Wenumass/F')
    outTree.Branch( 'WenuPhi', WenuPhi, 'WenuPhi/F')

    outTree.Branch( 'WmunuRecoil', WmunuRecoil, 'WmunuRecoil/F')
    outTree.Branch( 'Wmunumass', Wmunumass, 'Wmunumass/F')
    outTree.Branch( 'WmunuPhi', WmunuPhi, 'WmunuPhi/F')

    outTree.Branch( 'ZeeRecoil', ZeeRecoil, 'ZeeRecoil/F')
    outTree.Branch( 'ZeeMass', ZeeMass, 'ZeeMass/F')
    outTree.Branch( 'ZeePhi', ZeePhi, 'ZeePhi/F')

    outTree.Branch( 'ZmumuRecoil', ZmumuRecoil, 'ZmumuRecoil/F')
    outTree.Branch( 'ZmumuMass', ZmumuMass, 'ZmumuMass/F')
    outTree.Branch( 'ZmumuPhi', ZmumuPhi, 'ZmumuPhi/F')

    outTree.Branch( 'TOPRecoil', TOPRecoil, 'TOPRecoil/F')
    outTree.Branch( 'TOPPhi', TOPPhi, 'TOPPhi/F')

    outTree.Branch( 'GammaRecoil', GammaRecoil, 'GammaRecoil/F')
    outTree.Branch( 'GammaPhi', GammaPhi, 'GammaPhi/F')

    #if len(sys.argv)>2:
     #   NEntries=int(sys.argv[2])
      #  print "WARNING: Running in TEST MODE"

    for ievent in range(NEntries):

#    print "\n*****\nWARNING: *Test run* Processing 200 events only.\n*****\n"
#    for ievent in range(200):
        if ievent%100==0: print "Processed "+str(ievent)+" of "+str(NEntries)+" events."
        skimmedTree.GetEntry(ievent)
        ## Get all relevant branches
        run                        = skimmedTree.__getattr__('runId')
        lumi                       = skimmedTree.__getattr__('lumiSection')
        event                      = skimmedTree.__getattr__('eventId')
#        print "Run:"+str(run)+"; Lumi:"+str(lumi)+"; Event:"+str(event)
        trigName                   = skimmedTree.__getattr__('hlt_trigName')
        trigResult                 = skimmedTree.__getattr__('hlt_trigResult')
        filterName                 = skimmedTree.__getattr__('hlt_filterName')
        filterResult               = skimmedTree.__getattr__('hlt_filterResult')
        BadPFMuonFilter            = skimmedTree.__getattr__('hlt_filterbadPFMuon')


        pfMet                      = skimmedTree.__getattr__('pfMetCorrPt')
        pfMetPhi                   = skimmedTree.__getattr__('pfMetCorrPhi')
        pfMetJetUnc                = skimmedTree.__getattr__('pfMetCorrUnc')


        nTHINJets                  = skimmedTree.__getattr__('THINnJet')
        thinjetP4                  = skimmedTree.__getattr__('THINjetP4')
        thinJetCSV                 = skimmedTree.__getattr__('THINjetCISVV2')
        passThinJetLooseID         = skimmedTree.__getattr__('THINjetPassIDLoose')
        THINjetHadronFlavor        = skimmedTree.__getattr__('THINjetHadronFlavor')
        THINjetNPV                 = skimmedTree.__getattr__('THINjetNPV')         #int()
        thinjetNhadEF              = skimmedTree.__getattr__('THINjetNHadEF')
        thinjetChadEF              = skimmedTree.__getattr__('THINjetCHadEF')
        thinjetCEmEF               = skimmedTree.__getattr__('THINjetCEmEF')
        thinjetPhoEF               = skimmedTree.__getattr__('THINjetPhoEF')
        thinjetEleEF               = skimmedTree.__getattr__('THINjetEleEF')
        thinjetMuoEF               = skimmedTree.__getattr__('THINjetMuoEF')
        thinjetCorrUnc             = skimmedTree.__getattr__('THINjetCorrUncUp')

        nTHINdeepCSVJets           = skimmedTree.__getattr__('AK4deepCSVnJet')
        thindeepCSVjetP4           = skimmedTree.__getattr__('AK4deepCSVjetP4')
        thinJetdeepCSV             = skimmedTree.__getattr__('AK4deepCSVjetDeepCSV_b')
        THINdeepCSVjetHadronFlavor = skimmedTree.__getattr__('AK4deepCSVjetHadronFlavor')
        thindeepCSVjetNhadEF       = skimmedTree.__getattr__('AK4deepCSVjetNHadEF')
        thindeepCSVjetChadEF       = skimmedTree.__getattr__('AK4deepCSVjetCHadEF')
        THINdeepCSVjetNPV          = skimmedTree.__getattr__('AK4deepCSVjetNPV')

        try:
            thindeepCSVjetCorrUnc      = skimmedTree.__getattr__('AK4deepCSVjetCorrUncUp')
        except:
            if ievent==0: print "\n**********WARNING: Looks like the ntuple is from an older version, as DeepCSV jet correction Unc is missing. DeepCSV jet correction Unc will NOT be stored.**********\n"
            thindeepCSVjetCorrUnc = 1.


        try:
            thindeepCSVJetLooseID      = skimmedTree.__getattr__('AK4deepCSVjetPassIDLoose')
        except:
            if ievent==0: print "\n**********WARNING: Looks like the ntuple is from an older version, as DeepCSV Loose ID is missing. DeepCSV jet ID will NOT be stored.**********\n"
            thindeepCSVJetLooseID = None

        nEle                       = skimmedTree.__getattr__('nEle')
        eleP4                      = skimmedTree.__getattr__('eleP4')
        eleIsPassLoose             = skimmedTree.__getattr__('eleIsPassLoose')
        eleIsPassMedium            = skimmedTree.__getattr__('eleIsPassMedium')
        eleIsPassTight             = skimmedTree.__getattr__('eleIsPassTight')
        eleCharge                  = skimmedTree.__getattr__('eleCharge')

        nMu                        = skimmedTree.__getattr__('nMu')
        muP4                       = skimmedTree.__getattr__('muP4')
        isLooseMuon                = skimmedTree.__getattr__('isLooseMuon')
        isMediumMuon               = skimmedTree.__getattr__('isMediumMuon')
        isTightMuon                = skimmedTree.__getattr__('isTightMuon')
        muChHadIso                 = skimmedTree.__getattr__('muChHadIso')
        muNeHadIso                 = skimmedTree.__getattr__('muNeHadIso')
        muGamIso                   = skimmedTree.__getattr__('muGamIso')
        muPUPt                     = skimmedTree.__getattr__('muPUPt')
        muCharge                   = skimmedTree.__getattr__('muCharge')

        nTau                       = skimmedTree.__getattr__('HPSTau_n')
        tauP4                      = skimmedTree.__getattr__('HPSTau_4Momentum')
        isDecayModeFinding         = skimmedTree.__getattr__('disc_decayModeFinding')
        passLooseTauIso            = skimmedTree.__getattr__('disc_byLooseIsolationMVA3oldDMwLT')

        disc_againstElectronLoose  = skimmedTree.__getattr__('disc_againstElectronLooseMVA5')
        disc_againstElectronMedium = skimmedTree.__getattr__('disc_againstElectronMediumMVA5')
        disc_againstElectronTight  = skimmedTree.__getattr__('disc_againstElectronTightMVA5')
        disc_againstMuonLoose      = skimmedTree.__getattr__('disc_againstMuonLoose3')
        disc_againstMuonTight      = skimmedTree.__getattr__('disc_againstMuonTight3')

        isData                     = skimmedTree.__getattr__('isData')
        mcWeight                   = skimmedTree.__getattr__('mcWeight')
        pu_nTrueInt                = skimmedTree.__getattr__('pu_nTrueInt')         #int()
        pu_nPUVert                 = skimmedTree.__getattr__('pu_nPUVert')

        nPho                       = skimmedTree.__getattr__('nPho')
        phoP4                      = skimmedTree.__getattr__('phoP4')
        phoIsPassLoose             = skimmedTree.__getattr__('phoIsPassLoose')
        phoIsPassMedium            = skimmedTree.__getattr__('phoIsPassMedium')
        phoIsPassTight             = skimmedTree.__getattr__('phoIsPassTight')

        nGenPar                    = skimmedTree.__getattr__('nGenPar')
        genParId                   = skimmedTree.__getattr__('genParId')
        genMomParId                = skimmedTree.__getattr__('genMomParId')
        genParSt                   = skimmedTree.__getattr__('genParSt')
        genParP4                   = skimmedTree.__getattr__('genParP4')


        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # MC Weights ----------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        #
        mcweight[0] = 0.0
        if isData==1:   mcweight[0] =  1.0
        if not isData :
            if mcWeight<0:  mcweight[0] = -1.0
            if mcWeight>0:  mcweight[0] =  1.0


        h_total.Fill(1.);
        h_total_mcweight.Fill(1.,mcweight[0]);

        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Trigger selection ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        #
        trigstatus=False
        for itrig in range(len(triglist)):
            exec(triglist[itrig]+" = CheckFilter(trigName, trigResult, " + "'" + triglist[itrig] + "')")        #Runs the above commented-off code dynamically.
            exec("if "+triglist[itrig]+": trigstatus=True")                                                     #If any of the trigs is true, the event is kept.
#            exec("trig"+str(itrig+1)+"="+triglist[itrig])                                                       #Saves them as trig1, trig2, etc. #Deprecated
            exec("st_"+triglist[itrig]+"[0]="+triglist[itrig])                                                  #Adds to SkimmedTree output.

        if not isData: trigstatus=True

        if not trigstatus: continue

        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Filter selection
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------

        ###################################MET FILTER###################################
        filterstatus = False
        filter1 = False; filter2 = False;filter3 = False;filter4 = False; filter5 = False; filter6 = False; filter7 =False; filter8 = False
        ifilter_=0
        filter1 = CheckFilter(filterName, filterResult, 'Flag_HBHENoiseFilter')
        filter2 = CheckFilter(filterName, filterResult, 'Flag_globalSuperTightHalo2016Filter')
        filter3 = CheckFilter(filterName, filterResult, 'Flag_eeBadScFilter')
        filter4 = CheckFilter(filterName, filterResult, 'Flag_goodVertices')
        filter5 = CheckFilter(filterName, filterResult, 'Flag_EcalDeadCellTriggerPrimitiveFilter')
        filter6 = BadPFMuonFilter
        #filter7 = BadChCandidate
        filter8 = CheckFilter(filterName, filterResult, 'Flag_HBHENoiseIsoFilter')
        if not isData:
            filterstatus = True
        if isData:
            filterstatus = filter1 & filter2 & filter3 & filter4 & filter5 & filter6 & filter8

        if filterstatus == False: continue
        ###################################MET FILTER###################################


        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## PFMET Selection
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------

#        if samplename=="all":
        pfmetstatus = ( pfMet >= 250.0 )
        if pfmetstatus == False : continue
        print 'pfMet: ', pfMet

        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------

        thinjetpassindex=[]
        nBjets=0
        for ithinjet in range(nTHINJets):
            j1 = thinjetP4[ithinjet]
            #if (j1.Pt() > 30.0)&(abs(j1.Eta())<2.4)&(bool(passThinJetLooseID[ithinjet])==True)&(bool(passThinJetPUID[ithinjet]) == True):
            if (j1.Pt() > 30.0)&(abs(j1.Eta())<4.5)&(bool(passThinJetLooseID[ithinjet])==True):
                thinjetpassindex.append(ithinjet)
                if thinJetCSV[ithinjet] > CSVMWP and abs(j1.Eta())<2.4 : nBjets += 1

        thindCSVjetpassindex=[]
        ndBjets=0


        for jthinjet in range(nTHINdeepCSVJets):
            j1 = thindeepCSVjetP4[jthinjet]
            if (j1.Pt() > 30.0):
                if thindeepCSVJetLooseID==None:
                    deepCSVJetLooseID=True
                else:
                    deepCSVJetLooseID=bool(passThinJetLooseID[jthinjet])

            if (j1.Pt() > 30.0)&(abs(j1.Eta())<4.5) and deepCSVJetLooseID: #  &(bool(passThinJetLooseID[jthinjet])==True):
                thindCSVjetpassindex.append(jthinjet)
            if thinJetdeepCSV[jthinjet] > DCSVMWP and abs(j1.Eta())<2.4 : ndBjets += 1


        if len(thinjetpassindex) < 1 and len(thindCSVjetpassindex) < 1 : continue

        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Photon Veto
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        myPhos=[]
        for ipho in range(nPho):
            if (phoP4[ipho].Pt() > 15.) & (abs(phoP4[ipho].Eta()) <2.5) & (bool(phoIsPassLoose[ipho]) == True):
                myPhos.append(ipho)

        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Electron Veto
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        myEles=[]
        for iele in range(nEle):
            if (eleP4[iele].Pt() > 10. ) & (abs(eleP4[iele].Eta()) <2.5) & (bool(eleIsPassLoose[iele]) == True) :
                myEles.append(iele)

        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Muon Veto
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        myMuos = []
        myMuonIso = {}
        for imu in range(nMu):
            if (muP4[imu].Pt()>10.) & (abs(muP4[imu].Eta()) < 2.4) & (bool(isLooseMuon[imu]) == True):
                relPFIso = (muChHadIso[imu]+ max(0., muNeHadIso[imu] + muGamIso[imu] - 0.5*muPUPt[imu]))/muP4[imu].Pt()
                if relPFIso<0.25 :
                    myMuos.append(imu)
                    myMuonIso[imu]=relPFIso

        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        ## Tau Veto
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        myTaus=[]
        nTausDRbased=0
        myTausTightElectron=[]
        myTausTightMuon=[]
        myTausTightEleMu=[]
        myTausLooseEleMu=[]
        for itau in range(nTau):
            if (tauP4[itau].Pt()>18.) & (abs(tauP4[itau].Eta())<2.3) & (bool(isDecayModeFinding[itau]) == True) & (bool(passLooseTauIso[itau]) == True):
                myTaus.append(itau)
                if disc_againstElectronLoose!=None:
                    if disc_againstElectronTight[itau] and disc_againstMuonLoose[itau]:
                        myTausTightElectron.append(tauP4[itau])
                    if disc_againstMuonTight[itau] and disc_againstElectronLoose[itau]:
                        myTausTightMuon.append(tauP4[itau])
                    if disc_againstMuonTight[itau] and disc_againstElectronTight[itau]:
                        myTausTightEleMu.append(tauP4[itau])
                    if disc_againstMuonLoose[itau] and disc_againstElectronLoose[itau]:
                        myTausLooseEleMu.append(tauP4[itau])
                #---Fake tau cleaner----
                isClean=True
                for iele in myEles:
                    lep_tau_dR=DeltaR(eleP4[iele],tauP4[itau])    # math.sqrt(  (  iele.Eta()-tauP4[itau].Eta() )**2  + (  DeltaPhi(iele.Phi(),tauP4[itau].Phi()) )**2 )
                    if lep_tau_dR < 0.4:
                        isClean=False
                        break
                for imu in myMuos:
                    lep_tau_dR=DeltaR(muP4[imu],tauP4[itau])
                    if lep_tau_dR < 0.4:
                        isClean=False
                        break
                if isClean: nTausDRbased+=1

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------



        st_runId[0]             = long(run)
        st_lumiSection[0]       = lumi
        st_eventId[0]           = event
#        print '-----------'+str(st_runId)+", "+str(st_lumiSection)+", "+str(st_eventId)
        st_isData[0]            = isData
        st_pfMetCorrPt[0]       = pfMet
        st_pfMetCorrPhi[0]      = pfMetPhi

        st_pfMetUncJetResUp.clear()
        st_pfMetUncJetResDown.clear()

        st_pfMetUncJetEnUp.clear()
        st_pfMetUncJetEnDown.clear()

        st_THINjetP4.clear()
        st_THINjetCISVV2.clear()
        st_THINjetHadronFlavor.clear()
        st_THINjetNHadEF.clear()
        st_THINjetCHadEF.clear()

        st_THINjetCEmEF.clear()
        st_THINjetPhoEF.clear()
        st_THINjetEleEF.clear()
        st_THINjetMuoEF.clear()
        st_THINjetCorrUnc.clear()

        st_AK4deepCSVjetP4.clear()
        st_AK4deepCSVjetDeepCSV_b.clear()
        st_AK4deepCSVjetHadronFlavor.clear()
        st_AK4deepCSVjetNHadEF.clear()
        st_AK4deepCSVjetCHadEF.clear()
        st_AK4deepCSVjetCorrUnc.clear()

        st_eleP4.clear()
        st_eleIsPassTight.clear()

        st_muP4.clear()
        st_isTightMuon.clear()
        st_muIso.clear()

        st_phoP4.clear()
        st_phoIsPassTight.clear()

        st_genParId.clear()
        st_genMomParId.clear()
        st_genParSt.clear()
        st_genParP4.clear()


        st_THINnJet[0] = len(thinjetpassindex)
        for ithinjet in thinjetpassindex:
            st_THINjetP4.push_back(thinjetP4[ithinjet])
            st_THINjetCISVV2.push_back(thinJetCSV[ithinjet])
            st_THINjetHadronFlavor.push_back(THINjetHadronFlavor[ithinjet])
            st_THINjetNHadEF.push_back(thinjetNhadEF[ithinjet])
            st_THINjetCHadEF.push_back(thinjetChadEF[ithinjet])

            st_THINjetCEmEF.push_back(thinjetCEmEF[ithinjet])
            st_THINjetPhoEF.push_back(thinjetPhoEF[ithinjet])
            st_THINjetEleEF.push_back(thinjetEleEF[ithinjet])
            st_THINjetMuoEF.push_back(thinjetMuoEF[ithinjet])
            st_THINjetCorrUnc.push_back(thinjetCorrUnc[ithinjet])

        try:
            st_AK4deepCSVnJet[0] = len(thindCSVjetpassindex)
            for ithinjet in thindCSVjetpassindex:
                st_AK4deepCSVjetP4.push_back(thindeepCSVjetP4[ithinjet])
                st_AK4deepCSVjetDeepCSV_b.push_back(thinJetdeepCSV[ithinjet])
                st_AK4deepCSVjetHadronFlavor.push_back(THINdeepCSVjetHadronFlavor[ithinjet])
                st_AK4deepCSVjetNHadEF.push_back(thindeepCSVjetNhadEF[ithinjet])
                st_AK4deepCSVjetCHadEF.push_back(thindeepCSVjetChadEF[ithinjet])
                st_AK4deepCSVjetCorrUnc.push_back(thindeepCSVjetCorrUnc[ithinjet])
        except:
            pass

        st_nEle[0] = len(myEles)
        for iele in myEles:
            st_eleP4.push_back(eleP4[iele])
            st_eleIsPassTight.push_back(bool(eleIsPassTight[iele]))

        st_nMu[0] = len(myMuos)
        for imu in myMuos:
            st_muP4.push_back(muP4[imu])
            st_isTightMuon.push_back(bool(isTightMuon[imu]))
            st_muIso.push_back(myMuonIso[imu])

        st_HPSTau_n[0] = len(myTaus)
        st_nTauTightElectron[0] = len(myTausTightElectron)
        st_nTauTightMuon[0] = len(myTausTightMuon)
        st_nTauTightEleMu[0] = len(myTausTightEleMu)
        st_nTauLooseEleMu[0] = len(myTausLooseEleMu)

        st_nPho[0]=nPho
        for ipho in range(nPho):
            st_phoP4.push_back(phoP4[ipho])
            st_phoIsPassTight.push_back(bool(phoIsPassTight[ipho]))


        st_pu_nTrueInt[0] = pu_nTrueInt
        st_pu_nPUVert[0] = pu_nPUVert
        st_THINjetNPV[0] = THINjetNPV
#        try:
        st_AK4deepCSVjetNPV[0] = THINdeepCSVjetNPV
#        except:
#            pass
#        print pu_nTrueInt
#        print st_pu_nTrueInt[0]
        st_nGenPar[0] =  nGenPar
        for igp in range(nGenPar):
            st_genParId.push_back(genParId[igp])
            st_genMomParId.push_back(genMomParId[igp])
            st_genParSt.push_back(genParSt[igp])
            st_genParP4.push_back(genParP4[igp])

        st_pfMetUncJetResUp.push_back(pfMetJetUnc[0])
        st_pfMetUncJetResDown.push_back(pfMetJetUnc[1])

        st_pfMetUncJetEnUp.push_back(pfMetJetUnc[2])
        st_pfMetUncJetEnDown.push_back(pfMetJetUnc[3])

        ## Fill variables for the CRs.
        WenuRecoil[0] = -1.0
        Wenumass[0] = -1.0
        WenuPhi[0] = -10.

        WmunuRecoil[0] = -1.0
        Wmunumass[0] = -1.0
        WmunuPhi[0] = -10.

        ZeeMass[0] = -1.0
        ZeeRecoil[0] = -1.0
        ZeePhi[0] = -10.

        ZmumuMass[0] = -1.0
        ZmumuRecoil[0] = -1.0
        ZmumuPhi[0] = -10.

        TOPRecoil[0] = -1.0
        TOPPhi[0] = -10.

        GammaRecoil[0] = -1.0
        GammaPhi[0]  = -10.

# ------------------
# Z CR
# ------------------

        ## for dielectron
        # if len(myEles) == 2:
        #
        #     iele1=myEles[0]
        #     iele2=myEles[1]
        #     p4_ele1 = eleP4[iele1]
        #     p4_ele2 = eleP4[iele2]
        #     if eleCharge[iele1]*eleCharge[iele2]<0:
        #         ee_mass = ( p4_ele1 + p4_ele2 ).M()
        #         zeeRecoilPx = -( pfMet*math.cos(pfMetPhi) + p4_ele1.Px() + p4_ele2.Px())
        #         zeeRecoilPy = -( pfMet*math.sin(pfMetPhi) + p4_ele1.Py() + p4_ele2.Py())
        #         ZeeRecoilPt =  math.sqrt(zeeRecoilPx**2  +  zeeRecoilPy**2)
        #         if ee_mass > 60.0 and ee_mass < 120.0 and ZeeRecoilPt > 200.:
        #             ZeeRecoil[0] = ZeeRecoilPt
        #             ZeeMass[0] = ee_mass
        #             ZeePhi[0] = arctan(zeeRecoilPx,zeeRecoilPy)
        #
        # ## for dimu
        # if len(myMuos) == 2:
        #     imu1=myMuos[0]
        #     imu2=myMuos[1]
        #     p4_mu1 = muP4[imu1]
        #     p4_mu2 = muP4[imu2]
        #     if muCharge[imu1]*muCharge[imu2]<0:
        #         mumu_mass = ( p4_mu1 + p4_mu2 ).M()
        #         zmumuRecoilPx = -( pfMet*math.cos(pfMetPhi) + p4_mu1.Px() + p4_mu2.Px())
        #         zmumuRecoilPy = -( pfMet*math.sin(pfMetPhi) + p4_mu1.Py() + p4_mu2.Py())
        #         ZmumuRecoilPt =  math.sqrt(zmumuRecoilPx**2  +  zmumuRecoilPy**2)
        #         if mumu_mass > 60.0 and mumu_mass < 120.0 and ZmumuRecoilPt > 200.:
        #             ZmumuRecoil[0] = ZmumuRecoilPt
        #             ZmumuMass[0] = mumu_mass
        #             ZmumuPhi[0] = arctan(zmumuRecoilPx,zmumuRecoilPy)
        #
        # if len(myEles) == 2:
        #     ZRecoilstatus =(ZeeRecoil[0] > 200)
        # elif len(myMuos) == 2:
        #     ZRecoilstatus =(ZmumuRecoil[0] > 200)
        # else:
        #     ZRecoilstatus=False


# ------------------
# W CR
# ------------------

        ## for Single electron
        if len(myEles) == 1:
           ele1 = myEles[0]
           p4_ele1 = eleP4[ele1]

           e_mass = MT(p4_ele1.Pt(),pfMet, DeltaPhi(p4_ele1.Phi(),pfMetPhi)) #transverse mass defined as sqrt{2pT*MET*(1-cos(dphi)}

           WenuRecoilPx = -( pfMet*math.cos(pfMetPhi) + p4_ele1.Px())
           WenuRecoilPy = -( pfMet*math.sin(pfMetPhi) + p4_ele1.Py())
           WenuRecoilPt = math.sqrt(WenuRecoilPx**2  +  WenuRecoilPy**2)
           #if WenuRecoilPt > 250.:
           print 'WenuRecoilPt: ', WenuRecoilPt
           WenuRecoil[0] = WenuRecoilPt
           Wenumass[0] = e_mass
           WenuPhi[0] = arctan(WenuRecoilPx,WenuRecoilPy)

        ## for Single muon
        if len(myMuos) == 1:
           mu1 = myMuos[0]
           p4_mu1 = muP4[mu1]

           mu_mass = MT(p4_mu1.Pt(),pfMet, DeltaPhi(p4_mu1.Phi(),pfMetPhi)) #transverse mass defined as sqrt{2pT*MET*(1-cos(dphi)}

           WmunuRecoilPx = -( pfMet*math.cos(pfMetPhi) + p4_mu1.Px())
           WmunuRecoilPy = -( pfMet*math.sin(pfMetPhi) + p4_mu1.Py())
           WmunuRecoilPt = math.sqrt(WmunuRecoilPx**2  +  WmunuRecoilPy**2)
           #if WmunuRecoilPt > 250.:
           print 'WmunuRecoilPt: ', WmunuRecoilPt
           WmunuRecoil[0] = WmunuRecoilPt
           Wmunumass[0] = mu_mass
           WmunuPhi[0] = arctan(WmunuRecoilPx,WmunuRecoilPy)


        if len(myEles) == 1:
            WRecoilstatus =(WenuRecoil[0] > 200.)
        elif len(myMuos) == 1:
            WRecoilstatus =(WmunuRecoil[0] > 200.)
        else:
            WRecoilstatus=False
        #if pfmetstatus==False and ZRecoilstatus==False and WRecoilstatus==False: continue
        #if pfmetstatus==False : continue
        outTree.Fill()

    h_total_mcweight.Write()
    h_total.Write()
    samplepath.Write()
    outfile.Write()


def CheckFilter(filterName, filterResult,filtercompare):
    ifilter_=0
    filter1 = False
    for ifilter in filterName:
        filter1 = (ifilter.find(filtercompare) != -1)  & (bool(filterResult[ifilter_]) == True)
        if filter1: break
        ifilter_ = ifilter_ + 1
    return filter1

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

def Phi_mpi_pi(x):
    kPI = 3.14159265358979323846
    kTWOPI = 2 * kPI

    while (x >= kPI): x = x - kTWOPI;
    while (x < -kPI): x = x + kTWOPI;
    return x;

def DeltaPhi(phi1,phi2):
   phi = Phi_mpi_pi(phi1-phi2)

   return abs(phi)

def CheckFilter(filterName, filterResult,filtercompare):
    ifilter_=0
    filter1 = False
    for ifilter in filterName:
        filter1 = (ifilter.find(filtercompare) != -1)  & (bool(filterResult[ifilter_]) == True)
        if filter1: break
        ifilter_ = ifilter_ + 1
    return filter1


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
            #print "inside WJ loop pdgid", PID
            #print ("if status =",      (abs(PID) != 11),( abs(PID) != 12),(  abs(PID) != 13 ),(  abs(PID) != 14),(  abs(PID) != 15),(  abs(PID) != 16))
            #print "and of if status ", ( (abs(PID) != 11) & (abs(PID) != 12) &  (abs(PID) != 13) & (abs(PID) != 14) &  (abs(PID) != 15) &  (abs(PID) != 16) )

            if ( (abs(PID) != 11) & (abs(PID) != 12) &  (abs(PID) != 13) & (abs(PID) != 14) &  (abs(PID) != 15) &  (abs(PID) != 16) ): continue
            #print "lepton found"
            if ( ( (status != 1) & (abs(PID) != 15)) | ( (status != 2) & (abs(PID) == 15)) ): continue
            #print "tau found"
            if ( (abs(momPID) != 24) & (momPID != PID) ): continue
            #print "W found"
            #print "aftrer WJ if statement"
            goodLepID.append(ig)
        #print "length = ",len(goodLepID)
        if len(goodLepID) == 2 :
            l4_thisLep = genParP4[goodLepID[0]]
            l4_thatLep = genParP4[goodLepID[1]]
            l4_z = l4_thisLep + l4_thatLep

            pt = l4_z.Pt()
            pt__ = pt
            print " pt inside "
            k2 = -0.830041 + 7.93714 *TMath.Power( pt - (-877.978) ,(-0.213831) ) ;

    #################
    #ZJets
    #################
    if sample == "ZJETS":
        print " inside zjets "
        goodLepID = []
        for ig in range(nGenPar):
         #   print " inside loop "
            PID    = genParId[ig]
            momPID = genMomParId[ig]
            status = genParSt[ig]
          #  print " after vars "

            if ( (abs(PID) != 12) &  (abs(PID) != 14) &  (abs(PID) != 16) ) : continue
            if ( status != 1 ) : continue
            if ( (momPID != 23) & (momPID != PID) ) : continue
            goodLepID.append(ig)

        if len(goodLepID) == 2 :
            l4_thisLep = genParP4[goodLepID[0]]
            l4_thatLep = genParP4[goodLepID[1]]
            l4_z = l4_thisLep + l4_thatLep
            pt = l4_z.Pt()
            print " pt inside "
            k2 = -0.180805 + 6.04146 *TMath.Power( pt - (-759.098) ,(-0.242556) ) ;

    #################
    #TTBar
    #################
    if (sample=="TT"):
        print " inside ttbar "
        goodLepID = []
        for ig in range(nGenPar):
            print "inside TT loop "
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
    AnalyzeDataSet()
