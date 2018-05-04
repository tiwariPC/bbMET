#!/usr/bin/env python
from ROOT import TFile, TTree, TH1F, TH1D, TH1, TCanvas, TChain,TGraphAsymmErrors, TMath, TH2D, TLorentzVector, AddressOf, gROOT, TNamed
import ROOT as ROOT
import os
import sys, optparse
from array import array
import math
import numpy as numpy_

ROOT.gROOT.LoadMacro("Loader.h+")
outfilename= "ntupleHistos.root"

ntuple = TChain("tree/treeMaker")
ntuple.Add(sys.argv[1])

NEntries = ntuple.GetEntries()

if len(sys.argv)>2:
    #if sys.argv[2]=="test":
    NEntries=int(sys.argv[2])
    print "WARNING: Running in TEST MODE"

print 'NEntries = '+str(NEntries)

CSVMWP=0.8484

h_total         =TH1F('h_total','h_total',2,0,2)
h_met           =TH1F('h_met_',  'h_met_',  100,0.,1000.)
h_njet          =TH1F('h_njet_',  'h_njet_',  20,0.,20.)
h_nbjet         =TH1F('h_nbjet_',  'h_nbjet_',  6,0.,6.)
h_jet1_pT       =TH1F('h_jet1_pT_',  'h_jet1_pT_',  80,0.,800.)
h_jet2_pT       =TH1F('h_jet2_pT_',  'h_jet2_pT_',  80,0.,400.)
h_jet3_pT       =TH1F('h_jet3_pT_',  'h_jet3_pT_',  80,0.,400.)
h_jet1_Chad     =TH1F('h_jet1_Chad_',  'h_jet1_Chad_',  20,0.,1.)
h_jet1_Nhad     =TH1F('h_jet1_Nhad_',  'h_jet1_Nhad_',  20,0.,1.)
h_jet_mindPhi   =TH1F('h_jet_mindPhi_',  'h_jet_mindPhi_',  70,0.,3.5)
h_jet_eta       =TH1F('h_jet_eta_',  'h_jet_eta_',  70,-3.5,3.5)
h_bjet_eta      =TH1F('h_bjet_eta_',  'h_bjet_eta_',  70,-3.5,3.5)
h_nPho          =TH1F('h_nPho_',  'h_nPho_',  8,0.,8.)
h_nTau          =TH1F('h_nTau_',  'h_nTau_',  5,0.,5.)

def Phi_mpi_pi(x):
    kPI = 3.14159265358979323846
    kTWOPI = 2 * kPI    
    while (x >= kPI): x = x - kTWOPI
    while (x < -kPI): x = x + kTWOPI
    return x

def DeltaPhi(phi1,phi2):
    phi = Phi_mpi_pi(phi1-phi2)

    return abs(phi)

for ievent in range(NEntries):
    if ievent%100==0: print "Processed %d of %d events..." %(ievent,NEntries)
    
    ntuple.GetEntry(ievent)
    pfMet                      = ntuple.__getattr__('pfMetCorrPt')
    pfMetPhi                   = ntuple.__getattr__('pfMetCorrPhi')
    
    nTHINJets                  = ntuple.__getattr__('THINnJet')
    thinjetP4                  = ntuple.__getattr__('THINjetP4')
    thinJetCSV                 = ntuple.__getattr__('THINjetCISVV2')
    passThinJetLooseID         = ntuple.__getattr__('THINjetPassIDLoose')
    thinjetNhadEF              = ntuple.__getattr__('THINjetNHadEF')
    thinjetChadEF              = ntuple.__getattr__('THINjetCHadEF')
    
#    nTau                       = ntuple.__getattr__('HPSTau_n')
    nPho                       = ntuple.__getattr__('nPho')
    
    mcWeight                   = ntuple.__getattr__('mcWeight')
    
    
#Jet cleaning
    
    thinjetpassindex=[]
    thinjetpassP4=[]
    nBjets=0
    bjetP4=[]
    for ithinjet in range(nTHINJets):
        j1 = thinjetP4[ithinjet]
        if (j1.Pt() > 30.0)&(abs(j1.Eta())<2.4)&(bool(passThinJetLooseID[ithinjet])==True):
            thinjetpassindex.append(ithinjet)
            thinjetpassP4.append(j1)
            if thinJetCSV[ithinjet] > CSVMWP:
                nBjets += 1
                bjetP4.append(thinjetP4[ithinjet])
    nTHINJets=len(thinjetpassindex)
    
    alljetPT=[jet.Pt() for jet in thinjetpassP4]
    jetindex=[i for i in range(len(alljetPT))]    
                
    sortedjets=[jet for pt,jet in sorted(zip(alljetPT,thinjetpassP4), reverse=True)]                 # This gives a list of jets with their pTs in descending order
    sortedindex=[jetindex for pt,jetindex in sorted(zip(alljetPT,jetindex), reverse=True)]           # Indices of jets in thinjetP4 in decscending order of jetPT
    
    if nTHINJets>0: j1=sortedjets[0]
    if nTHINJets>1: j2=sortedjets[1]
    if nTHINJets>2: j3=sortedjets[2]
    
    if nTHINJets>0: ifirstjet=sortedindex[0]
#    if nTHINJets>1: isecondjet=sortedindex[1]
#    if nTHINJets>2: ithirdjet=sortedindex[2]
    
#Fill

    if nTHINJets>0: h_jet1_pT.Fill(j1.Pt(),mcWeight)
    if nTHINJets>1: h_jet2_pT.Fill(j2.Pt(),mcWeight)
    if nTHINJets>2: h_jet3_pT.Fill(j3.Pt(),mcWeight)
    
    h_met.Fill(pfMet,mcWeight)
    h_njet.Fill(nTHINJets,mcWeight)
    h_nbjet.Fill(nBjets,mcWeight)
    
    try:
        if nTHINJets>0:
            h_jet1_Chad.Fill(thinjetChadEF[thinjetpassindex[ifirstjet]],mcWeight)
            h_jet1_Nhad.Fill(thinjetNhadEF[thinjetpassindex[ifirstjet]],mcWeight)
    except:
        print len(thinjetpassP4)
        print len(thinjetpassindex)
        print ifirstjet
        sys.exit()
    
    h_nPho.Fill(nPho,mcWeight)
#    h_nTau.Fill(nTau,mcWeight)
    
    for ijet in thinjetpassP4:
        h_jet_eta.Fill(ijet.Eta(),mcWeight)
    for ibjet in bjetP4:
        h_bjet_eta.Fill(ibjet.Eta(),mcWeight)
    
    if nTHINJets > 0:
        mindphi=10
        for ithinjet in range(nTHINJets):
            dp=DeltaPhi(thinjetpassP4[ithinjet].Phi(),pfMetPhi)
            if dp < mindphi: mindphi = dp
        h_jet_mindPhi.Fill(mindphi,mcWeight)
    
#Write

f = TFile(outfilename,'RECREATE')
f.cd()
h_total.SetBinContent(1,NEntries)
h_total.Write()

h_met.Write()
h_njet.Write()
h_nbjet.Write()
h_jet1_pT.Write()
h_jet2_pT.Write()
h_jet3_pT.Write()
h_jet1_Chad.Write()
h_jet1_Nhad.Write()
h_jet_mindPhi.Write()
h_jet_eta.Write()
h_bjet_eta.Write()
#h_nTau.Write()
h_nPho.Write()

print "Histograms written to %s." %outfilename
