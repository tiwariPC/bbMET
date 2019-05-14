import os
import sys
import datetime
import sys, optparse
## Ratio is added Data/MC
## Template macro is fed to a python variable
## 1.)  is created on DateBase
## 2.) Starting Extension of your Dir..Like
## in a day you want 2 directories jsut
## change the DirPreName
## Monika Mittal Khuarana
## Raman Khurana


#if len(sys.argv) < 2 :
#    print "insufficiency inputs provided, please provide the directory with input files"
##just for argument, the input file path is explicitly provided in line no 101
#if len(sys.argv) ==2 :
#    print "plotting from directory ",sys.argv[1]
#    inputdirname = sys.argv[1]

usage = "usage: %prog [options] arg1 arg2"
parser = optparse.OptionParser(usage)

parser.add_option("-d", "--data", dest="datasetname")
parser.add_option("-s", "--sr", action="store_true", dest="plotSRs")
parser.add_option("-m", "--mu", action="store_true", dest="plotMuRegs")
parser.add_option("-e", "--ele", action="store_true", dest="plotEleRegs")
parser.add_option("-p", "--pho", action="store_true", dest="plotPhoRegs")
parser.add_option("-q", "--qcd", action="store_true", dest="plotQCDRegs")
parser.add_option("-v", "--verbose", action="store_true", dest="verbose")

(options, args) = parser.parse_args()

if options.plotSRs==None:
    makeSRplots = False
else:
    makeSRplots = options.plotSRs

if options.plotMuRegs==None:
    makeMuCRplots = False
else:
    makeMuCRplots = options.plotMuRegs

if options.plotEleRegs==None:
    makeEleCRplots = False
else:
    makeEleCRplots = options.plotEleRegs

if options.plotPhoRegs==None:
    makePhoCRplots = False
else:
    makePhoCRplots = options.plotPhoRegs

if options.plotQCDRegs==None:
    makeQCDCRplots = False
else:
    makeQCDCRplots = options.plotQCDRegs

if options.verbose==None:
    verbose = False
else:
    verbose = options.verbose

if options.datasetname.upper()=="SE":
    dtset="SE"
elif options.datasetname.upper()=="SP":
    dtset="SP"
elif options.datasetname.upper()=="SM":
    dtset="SM"
else:
    dtset="MET"

print "Using dataset "+dtset

datestr = datetime.date.today().strftime("%d%m%Y")

macro='''
#include <ctime>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "TStyle.h"
#include "TString.h"
#include "TRegexp.h"

void Plot() {
   time_t now = time(0);
   tm * ltm = localtime( & now);
   TString dirpathname;
   bool varbin = true;

   TString DirPreName = "/afs/cern.ch/work/p/ptiwari/bb+DM_analysis/2016_bbMET_WCRSplit/bbMETplot/Scripts/test/";
   dirpathname = "'''+datestr+'''"; //.Form("%d%1.2d%d",ltm->tm_mday,1 + ltm->tm_mon,1900 + ltm->tm_year);

   system("mkdir -p  " + DirPreName + dirpathname + "/bbMETROOT");
   system("mkdir -p  " + DirPreName + dirpathname + "/bbMETPdf");
   system("mkdir -p  " + DirPreName + dirpathname + "/bbMETPng");

   ofstream mout;
   mout.open(DirPreName + dirpathname + "/HISTPATH" + dirpathname + "Integral.txt", std::ios::app);
   ofstream rout;
   rout.open(DirPreName + dirpathname + "/HISTPATH" + dirpathname + "Integral.html", std::ios::app);
   ofstream tableout;
   tableout.open(DirPreName + dirpathname + "/HISTPATH" + dirpathname + "IntegralWithError.txt", std::ios::app);
   TString outputshapefilename = DirPreName + dirpathname + "/HISTPATH.root";
   TFile *fshape = new TFile(outputshapefilename, "RECREATE");

   if (VARIABLEBINS) {
      ofstream metbinsout_1;
      ofstream metbinsout_2;
      ofstream metbinsout_3;

      system("mkdir -p  " + DirPreName + "METBIN_1");
      system("mkdir -p  " + DirPreName + "METBIN_2");
      system("mkdir -p  " + DirPreName + "METBIN_3");

      metbinsout_1.open(DirPreName + "METBIN_1/HISTPATH" + dirpathname + "Integral.txt", std::ios::app);
      metbinsout_2.open(DirPreName + "METBIN_2/HISTPATH" + dirpathname + "Integral.txt", std::ios::app);
      metbinsout_3.open(DirPreName + "METBIN_3/HISTPATH" + dirpathname + "Integral.txt", std::ios::app);
   }

   gROOT -> ProcessLine(".L tdrstyle.C");
   //setTDRStyle();
   gStyle -> SetOptStat(0);
   gStyle -> SetOptTitle(0);
   gStyle -> SetFrameLineWidth(3);
   //gStyle->SetErrorX(0);
   gStyle -> SetLineWidth(1);

   //Provide luminosity of total data
   //cout << endl << "*************** WARNING: Assuming incomplete data: full luminosity is not used. ***************" << endl << endl; //Adjust the last factor in the next line according to available data.
   float lumi = 35.9 * 1000; // 35.540082604 * 1000; // 3.1054555906 * 1000 ; ////11706.;// * 0.96838; // from the twiki https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2016Analysis  //36.773
   float luminosity = 35.9; // It will print on your plots too

   std::vector < TString > filenameBkgString;
   std::vector < TString > filenameDataString;
   std::vector < TString > filenameSigString;
   //Change here Directories of the file

   TH1F *DIBOSON;
   TH1F *TT;
   TH1F *WJets;
   TH1F *DYJets;
   TH1F *ZJets;
   TH1F *STop;
   TH1F *QCD;
   //TH1F*  data_obs;
   TString filenamepath("../../../bbMET/BR_Condor_Farmout_09052019/hadd_outputs_correctedTopPt/");

   // Diboson WW 0 1 2
   filenameBkgString.push_back(filenamepath + "Output_WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8.root"); // 49.997
   filenameBkgString.push_back(filenamepath + "Output_WWTo2L2Nu_13TeV-powheg.root"); //12.178
   filenameBkgString.push_back(filenamepath + "Output_WWTo4Q_4f_13TeV_amcatnloFXFX_madspin_pythia8.root"); //51.723

   // Diboson WZ 3 4 5 6 7
   filenameBkgString.push_back(filenamepath + "Output_WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8.root"); //10.71
   filenameBkgString.push_back(filenamepath + "Output_WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8.root"); //3.0330
   filenameBkgString.push_back(filenamepath + "Output_WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root"); //5.5950
   filenameBkgString.push_back(filenamepath + "Output_WZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8.root"); //6.4880
   filenameBkgString.push_back(filenamepath + "Output_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root"); //4.4297

   // Diboson ZZ 8 9 10 11
   filenameBkgString.push_back(filenamepath + "Output_ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root"); //3.22
   filenameBkgString.push_back(filenamepath + "Output_ZZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8.root"); //4.04
   filenameBkgString.push_back(filenamepath + "Output_ZZTo4L_13TeV_powheg_pythia8.root"); //1.2120
   filenameBkgString.push_back(filenamepath + "Output_ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8.root"); //6.842

   // Diboson TTV 12 13 14 15
   filenameBkgString.push_back(filenamepath + "Output_TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root"); //0.2043
   filenameBkgString.push_back(filenamepath + "Output_TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root"); // 0.4062
   filenameBkgString.push_back(filenamepath + "Output_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root"); //0.2529
   filenameBkgString.push_back(filenamepath + "Output_TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root"); //0.5297

   //ZJets High pt DYSample 16,17,18,19,20,21,22
   filenameBkgString.push_back(filenamepath + "Output_ZJetsToNuNu_HT-100To200_13TeV-madgraph_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_ZJetsToNuNu_HT-200To400_13TeV-madgraph_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_ZJetsToNuNu_HT-400To600_13TeV-madgraph_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_ZJetsToNuNu_HT-600To800_13TeV-madgraph_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_ZJetsToNuNu_HT-800To1200_13TeV-madgraph_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_ZJetsToNuNu_HT-1200To2500_13TeV-madgraph_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph_MC25ns_LegacyMC_20170328.root");

   //DYJets High pt DYSample 23,24,25,26,27,28,29,
   filenameBkgString.push_back(filenamepath + "Output_DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_LegacyMC_20170328.root");

   // WJets in Bins  30,32,33,34,35,36,37
   filenameBkgString.push_back(filenamepath + "Output_WJetsToLNu_HT-70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_LegacyMC_20170328.root");

   // Single Top 38,39,40,41,42
   filenameBkgString.push_back(filenamepath + "Output_ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1_MC25ns_LegacyMC_2017.root");
   filenameBkgString.push_back(filenamepath + "Output_ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1_MC25ns_LegacyMC.root");
   filenameBkgString.push_back(filenamepath + "Output_ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_MC25ns_LegacyMC_20170328.root");
   filenameBkgString.push_back(filenamepath + "Output_ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_MC25ns_LegacyMC_20170328.root");

   // TTJets 43
   filenameBkgString.push_back(filenamepath + "Output_TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8.root");

   // QCD 44,45,46,47,48,49,50,51,52
   filenameBkgString.push_back(filenamepath + "Output_QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"); //dummy
   filenameBkgString.push_back(filenamepath + "Output_QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"); //dummy
   filenameBkgString.push_back(filenamepath + "Output_QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"); //dummy
   filenameBkgString.push_back(filenamepath + "Output_QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"); //dummy
   filenameBkgString.push_back(filenamepath + "Output_QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
   filenameBkgString.push_back(filenamepath + "Output_QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
   filenameBkgString.push_back(filenamepath + "Output_QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
   filenameBkgString.push_back(filenamepath + "Output_QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
   filenameBkgString.push_back(filenamepath + "Output_QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
   //

   // not used so far
   TString filenamesigpath("../../../bbMET/BR_Condor_Farmout/signal/");
   //bbMET Signal Sample 53 - 90
   filenameSigString.push_back(filenamesigpath + "Output_scalar_NLO_Mchi-50_Mphi-400.root");
   filenameSigString.push_back(filenamesigpath + "Output_scalar_NLO_Mchi-50_Mphi-350.root");
   filenameSigString.push_back(filenamesigpath + "Output_scalar_NLO_Mchi-50_Mphi-300.root");
   filenameSigString.push_back(filenamesigpath + "Output_scalar_NLO_Mchi-450_Mphi-1000.root");
   filenameSigString.push_back(filenamesigpath + "Output_scalar_NLO_Mchi-350_Mphi-750.root");
   filenameSigString.push_back(filenamesigpath + "Output_scalar_NLO_Mchi-350_Mphi-1000.root");
   filenameSigString.push_back(filenamesigpath + "Output_scalar_NLO_Mchi-1_Mphi-750.root");
   filenameSigString.push_back(filenamesigpath + "Output_scalar_NLO_Mchi-1_Mphi-500.root");
   filenameSigString.push_back(filenamesigpath + "Output_scalar_NLO_Mchi-1_Mphi-400.root");
   filenameSigString.push_back(filenamesigpath + "Output_scalar_NLO_Mchi-1_Mphi-350.root");
   filenameSigString.push_back(filenamesigpath + "Output_scalar_NLO_Mchi-1_Mphi-300.root");
   filenameSigString.push_back(filenamesigpath + "Output_scalar_NLO_Mchi-10_Mphi-50.root");
   filenameSigString.push_back(filenamesigpath + "Output_scalar_NLO_Mchi-10_Mphi-10.root");
   filenameSigString.push_back(filenamesigpath + "Output_scalar_NLO_Mchi-10_Mphi-100.root");
   filenameSigString.push_back(filenamesigpath + "Output_scalar_NLO_Mchi-100_Mphi-500.root");
   filenameSigString.push_back(filenamesigpath + "Output_scalar_NLO_Mchi-100_Mphi-400.root");
   filenameSigString.push_back(filenamesigpath + "Output_scalar_NLO_Mchi-100_Mphi-350.root");
   filenameSigString.push_back(filenamesigpath + "Output_scalar_NLO_Mchi-100_Mphi-300.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-50_Mphi-50.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-50_Mphi-500.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-50_Mphi-400.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-50_Mphi-350.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-50_Mphi-300.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-50_Mphi-200.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-450_Mphi-1000.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-1_Mphi-50.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-1_Mphi-500.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-1_Mphi-400.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-1_Mphi-350.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-1_Mphi-100.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-1_Mphi-1000.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-10_Mphi-50.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-10_Mphi-10.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-10_Mphi-100.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-100_Mphi-750.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-100_Mphi-500.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-100_Mphi-400.root");
   filenameSigString.push_back(filenamesigpath + "Output_pseudo_NLO_Mchi-100_Mphi-350.root");

   //

   TString filenamedatapath("../../../bbMET/BR_Condor_Farmout_09052019/hadd_outputs/");
   //Data File 0
   filenameDataString.push_back(filenamedatapath + "data_combined_'''+dtset+'''.root");

   //const int n_integral = (int)filenameString.size();

   TString histnameString("HISTNAME");

   TFile * fIn;
   const int nfiles_bkg = (int) filenameBkgString.size();
   float Integral[nfiles_bkg], Integral_Error[nfiles_bkg];
   //kfactor * lo crossection
   //check it once
   float Xsec[nfiles_bkg];

   Xsec[0] = 49.997; // WW
   Xsec[1] = 12.178; // WW
   Xsec[2] = 51.723; // WW

   Xsec[3] = 10.71; // WZ
   Xsec[4] = 3.0330; // WZ
   Xsec[5] = 5.5950; // WZ
   Xsec[6] = 6.4880; // WZ
   Xsec[7] = 4.4297; // WZ

   Xsec[8] = 3.22; // ZZ
   Xsec[9] = 4.04; // ZZ
   Xsec[10] = 1.2120; // ZZ
   Xsec[11] = 6.842; // ZZ

   Xsec[12] = 0.2043; // TTV
   Xsec[13] = 0.4062; // TTV
   Xsec[14] = 0.2529; // TTV
   Xsec[15] = 0.5297; // TTV

   //float Znunu=1.23;
   float Znunu = 1.0;
   Xsec[16] = Znunu * 280.35; // Znunu HT 100-200
   Xsec[17] = Znunu * 77.67; // Znunu HT 200-400
   Xsec[18] = Znunu * 10.73; // Znunu HT 400-600
   Xsec[19] = Znunu * 2.559; // Znunu HT 600-800
   Xsec[20] = Znunu * 1.1796; // Znunu HT 800-1200
   Xsec[21] = Znunu * 0.28833; // Znunu HT 1200-2500
   Xsec[22] = Znunu * 0.006945; // Znunu HT 2500-inf

   Xsec[23] = Znunu * 147.4; // DYJetsToLL HT 100-200
   Xsec[24] = Znunu * 40.99; // DYJetsToLL HT 200-400
   Xsec[25] = Znunu * 5.678; // DYJetsToLL HT 400-600
   Xsec[26] = Znunu * 1.367; // DYJetsToLL HT 600-800
   Xsec[27] = Znunu * 0.6304; // DYJetsToLL HT 800-1200
   Xsec[28] = Znunu * 0.1514; // DYJetsToLL HT 1200-2500
   Xsec[29] = Znunu * 0.003565; // DYJetsToLL HT 2500-inf

   float wjets = 1.0;
   //float wjets = 1.21;
   Xsec[30] = wjets * 1372; // WJets HT 70-100     ***not available in twiki*** https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
   Xsec[31] = wjets * 1345; // WJets HT 100-200
   Xsec[32] = wjets * 1345; // WJets HT 200-400
   Xsec[33] = wjets * 48.91; // WJets HT 400-600
   Xsec[34] = wjets * 12.05; // WJets HT 600-800
   Xsec[35] = wjets * 5.501; // WJets HT 800-1200
   Xsec[36] = wjets * 1.329; // WJets HT 1200-2500
   Xsec[37] = wjets * 0.0322; // WJets HT 2500-Inf

   Xsec[38] = 136.02; // single top t-channel_top_4f_inclusiveDecays       ***not available in twiki***
   Xsec[39] = 80.95; // single top t-channel_antitop_4f_inclusiveDecays   ***not available in twiki***
   Xsec[40] = 3.36; // single top s-channel_4f_leptonDecays
   Xsec[41] = 35.85; // single top tW_top_5f_inclusiveDecays
   Xsec[42] = 35.85; // single top tW_antitop_5f_inclusiveDecays

   Xsec[43] = 364.35; // ttbar     ***not available in twiki***

   float QSF=1.0;
   Xsec[44] = 0; //246300000;              // QCD_HT50to100
   Xsec[45] = 0; //27990000;               // QCD_HT100to200
   Xsec[46] = 0; //1712000 * QCDSF;                // QCD_HT200to300
   Xsec[47] = 0; //347700 * QCDSF;                 // QCD_HT300to500
   Xsec[48] = QSF * 32100 * QCDSF; // QCD_HT500to700
   Xsec[49] = QSF * 6831 * QCDSF; // QCD_HT700to1000
   Xsec[50] = QSF * 1207 * QCDSF; // QCD_HT1000to1500
   Xsec[51] = QSF * 119.9 * QCDSF; // QCD_HT1500to2000
   Xsec[52] = QSF * 25.24 * QCDSF; // QCD_HT2000toInf


   TH1F * h_mc[nfiles_bkg];
   float normalization[nfiles_bkg];
   TH1F * h_data;
   TH1F * h_temp;
   TH1F * hnew;
   TH1F * h_total;

   for (int i = 0; i < 53; i++) {
      if (VERBOSE) cout << "Reading file #" << to_string(i + 1) << ": " << filenameBkgString[i] << endl;
      fIn = new TFile(filenameBkgString[i], "READ");
      if (varbin && (histnameString.Index("MET") != string::npos || histnameString.Index("met") != string::npos)) {
         std::cout << "Setting variable bining for recoil" << std::endl;
         h_temp = (TH1F * ) fIn -> Get(histnameString);
         Double_t bins[16] = {250,270,290,310,330,350,370,390,410,430,450,470,490,510,530,550};
         std::cout << "start rebining histogram" << std::endl;
         std::cout << "Integral of first hist" << h_temp -> Integral() << std::endl;
         TH1F * h_ = (TH1F * ) h_temp -> Rebin(15, "hnew", bins);
         h_ -> SetBinContent(15, h_ -> GetBinContent(15) + h_ -> GetBinContent(16));
         h_ -> SetBinContent(16, 0.);
         std::cout << "Integral of new hist" << h_ -> Integral() << std::endl;
         std::cout << "done rebining histogram" << std::endl;
         h_mc[i] = (TH1F * ) h_ -> Clone();
         std::cout << "done cloning" << std::endl;
      }
      else {
         h_mc[i] = (TH1F * ) fIn -> Get(histnameString);
         h_mc[i] -> Rebin(REBIN);
         h_mc[i] -> Sumw2();
      }

      h_total = (TH1F * ) fIn -> Get("h_total");
      if (h_total -> Integral() > 0) normalization[i] = (lumi * Xsec[i]) / (h_total -> Integral() * BLINDFACTOR);
      else normalization[i] = 0;
      cout<<"h_total :" << h_total -> Integral() << std::endl;

      Integral[i] = h_mc[i] -> Integral();
      if (Integral[i] <= 0) Integral_Error[i] = 0.0;
      if (Integral[i] > 0) Integral_Error[i] = TMath::Sqrt(Integral[i]) * normalization[i];
      h_mc[i] -> Scale(normalization[i]);
   }

   fIn = new TFile(filenameDataString[0], "READ");

   if (varbin && (histnameString.Index("MET") != string::npos || histnameString.Index("met") != string::npos)) {
      std::cout << "Setting variable bining for recoil" << std::endl;
      h_temp = (TH1F * ) fIn -> Get(histnameString);
      Double_t bins[16] = {250,270,290,310,330,350,370,390,410,430,450,470,490,510,530,550};
      //Int_t binnum = sizeof(bins)/sizeof(Double_t) - 1;
      TH1F * h_ = (TH1F * ) h_temp -> Rebin(16 - 1, "hnew", bins);
      h_ -> SetBinContent(15, h_ -> GetBinContent(15) + h_ -> GetBinContent(16));
      h_ -> SetBinContent(16, 0.);
      h_data = (TH1F * ) h_ -> Clone();
    }
    else {
      h_data = (TH1F * ) fIn -> Get(histnameString);
      h_data -> Rebin(REBIN);
      h_data -> Sumw2();
   }

   data_obs = (TH1F * ) fIn -> Get(histnameString);

   DIBOSON = (TH1F * ) h_mc[0] -> Clone();
   for (int dibo1 = 1; dibo1 < 16; dibo1++) {
      DIBOSON -> Add(h_mc[dibo1]);
   }

   ZJets = (TH1F * ) h_mc[16] -> Clone();
   for (int zjets1 = 17; zjets1 < 23; zjets1++) {
      ZJets -> Add(h_mc[zjets1]);
   }

   DYJets = (TH1F * ) h_mc[23] -> Clone();
   for (int DYjets = 24; DYjets < 30; DYjets++) {
      DYJets -> Add(h_mc[DYjets]);
   }

   WJets = (TH1F * ) h_mc[30] -> Clone();
   for (int wjets1 = 31; wjets1 < 38; wjets1++) {
      WJets -> Add(h_mc[wjets1]);
   }

   STop = (TH1F * ) h_mc[38] -> Clone();
   for (int ttjets = 39; ttjets < 43; ttjets++) {
      STop -> Add(h_mc[ttjets]);
   }

   TT = (TH1F * ) h_mc[43] -> Clone();

   QCD = (TH1F * ) h_mc[44] -> Clone();
   for (int qcd = 45; qcd < 53; qcd++) {
      QCD -> Add(h_mc[qcd]);
   }

   float ZJetsCount = ZJets -> Integral();
   float DYJetsCount = DYJets -> Integral();
   float WJetsCount = WJets -> Integral();
   float STopCount = STop -> Integral();
   float TTCount = TT -> Integral();
   float VVCount = DIBOSON -> Integral();
   float QCDCount = QCD -> Integral();

   TString DYLegend, WLegend, GLegend, ZLegend, STLegend, TTLegend, VVLegend, QCDLegend;

   //if (ISCUTFLOW) {
   if (1) {
      DYLegend = "Z(ll) + jets";
      WLegend = "W(l#nu) + jets";
      ZLegend = "Z(#nu#nu) + jets";
      STLegend = "Single t";
      TTLegend = "Top";
      VVLegend = "VV, VH";
      QCDLegend = "Multijet";
   } else {
      DYLegend = "Z(ll) + jets: " + std::to_string(int(DYJetsCount));
      WLegend = "W(l#nu) + jets: " + std::to_string(int(WJetsCount));
      ZLegend = "Z(#nu#nu) + jets: " + std::to_string(int(ZJetsCount));
      STLegend = "Single t: " + std::to_string(int(STopCount));
      TTLegend = "Top: " + std::to_string(int(TTCount));
      VVLegend = "VV: " + std::to_string(int(VVCount));
      QCDLegend = "QCD Multijet: " + std::to_string(int(QCDCount));
   }

   //NORATIOPLOT=0;

   //Legend
   TLegend * legend;

   float zj_i = ZJets -> Integral();
   float dyj_i = DYJets -> Integral();
   float wj_i = WJets -> Integral();
   float tt_i = TT -> Integral();
   float st_i = STop -> Integral();
   float db_i = DIBOSON -> Integral();
   float qc_i = QCD -> Integral();

   float mcsum = zj_i + dyj_i + wj_i + tt_i + st_i + db_i + qc_i;

   legend = new TLegend(0.60, 0.70, 0.94, 0.94, NULL, "brNDC");
   legend -> SetTextSize(0.020);
   legend -> SetBorderSize(0);
   legend -> SetLineColor(1);
   legend -> SetLineStyle(1);
   legend -> SetLineWidth(1);
   legend -> SetFillColor(0);
   legend -> SetFillStyle(0);
   legend -> SetTextFont(42);
   legend -> SetNColumns(2);

   cout << "Total:" << mcsum << endl;

   float legendthres = 0.008;

   if (!NORATIOPLOT) legend -> AddEntry(h_data, "Data", "PEL");
   if (dyj_i / mcsum > legendthres) legend -> AddEntry(DYJets, DYLegend, "f");
   if (zj_i / mcsum > legendthres) legend -> AddEntry(ZJets, ZLegend, "f");
   if (wj_i / mcsum > legendthres) legend -> AddEntry(WJets, WLegend, "f");
   if (tt_i / mcsum > legendthres) legend -> AddEntry(TT, TTLegend, "f");
   if (st_i / mcsum > 0.) legend -> AddEntry(STop, STLegend, "f");
   if (db_i / mcsum > 0.) legend -> AddEntry(DIBOSON, VVLegend, "f");
   if (qc_i / mcsum > 0.) legend -> AddEntry(QCD, QCDLegend, "f");

   //============== CANVAS DECLARATION ===================
   TCanvas * c12 = new TCanvas("Hist", "Hist", 0, 0, 1000, 1000);

   //==================Stack==============================
   THStack * hs = new THStack("hs", " ");

   //Colors for Histos
   DYJets -> SetFillColor(kGreen + 2);
   DYJets -> SetLineWidth(0);

   ZJets -> SetFillColor(kAzure + 1);
   ZJets -> SetLineWidth(0);

   DIBOSON -> SetFillColor(kBlue + 2);
   DIBOSON -> SetLineWidth(0);

   TT -> SetFillColor(kOrange - 2);
   TT -> SetLineWidth(0);

   WJets -> SetFillColor(kViolet - 3);
   WJets -> SetLineWidth(0);

   STop -> SetFillColor(kOrange + 1);
   STop -> SetLineWidth(0);

   QCD -> SetFillColor(kGray + 1);
   QCD -> SetLineWidth(0);


   hs -> Add(DIBOSON, "hist");
   hs -> Add(QCD, "hist");
   hs -> Add(STop, "hist");
   hs -> Add(TT, "hist");
   hs -> Add(WJets, "hist");
   hs -> Add(ZJets, "hist");
   hs -> Add(DYJets, "hist");

   h_data -> SetMarkerColor(kBlack);
   h_data -> SetMarkerStyle(20);
   //float maxi = h_data->GetMaximum();

   TH1F * Stackhist = (TH1F * ) hs -> GetStack() -> Last();

   hasNoEvents = false;
   float maxi = Stackhist -> GetMaximum();
   if (Stackhist -> GetEntries() == 0) {
      hasNoEvents = true;
      cout << "=============================" << endl << "No events found!" << endl << "=============================" << endl;
      fstream empfile("Empty.txt", ios::app);
      empfile << "HISTNAME" << endl;
      empfile.close();
   }

   TH1F * h_err;
   h_err = (TH1F * ) h_data -> Clone("h_err");
   h_err = (TH1F * ) h_mc[0] -> Clone("h_err");
   h_err -> Sumw2();
   h_err -> Reset();
   for (int imc = 1; imc < 53; imc++) {
      h_err -> Add(h_mc[imc]);
   }
   Stackhist -> SetLineWidth(2);

   c12 -> SetLogy(ISLOG);

   // Upper canvas declaration
   TPad * c1_2 = NULL;
   if (NORATIOPLOT) {
      c1_2 = new TPad("c1_2", "newpad", 0, 0.05, 1, 1); //0.993);
      c1_2 -> SetRightMargin(0.04);
   } else {
      c1_2 = new TPad("c1_2", "newpad", 0, 0.28, 1, 1);
   }

   c1_2 -> SetBottomMargin(0.03);
   c1_2 -> SetTopMargin(0.06);
   c1_2 -> SetLogy(ISLOG);
   if (VARIABLEBINS) {
      c1_2 -> SetLogx(0);
   }
   c1_2 -> Draw();
   c1_2 -> cd();

   hs -> Draw();

   std::cout << " PREFITDIR " << std::endl;
   TH1F * h_prefit;
   TFile * fprefit;
   TString prefitfilename = "PREFITDIR/HISTPATH.root";
   if (DRAWPREFIT) {
      fprefit = new TFile(prefitfilename, "READ");
      h_prefit = (TH1F * ) fprefit -> Get("bkgSum");
      std::cout << " inside prefit loop " << h_prefit -> Integral() << std::endl;
      h_prefit -> SetLineColor(kRed);
      h_prefit -> SetLineWidth(3);
      h_prefit -> SetFillColor(0);

   }

   gStyle -> SetHistTopMargin(0.);

   TH1F * Stackhist1 = (TH1F * ) hs -> GetStack() -> Last();
   h_err -> Draw("E2 SAME");
   h_err -> Sumw2();
   h_err -> SetFillColor(kGray + 3);
   h_err -> SetLineColor(kGray + 3);
   h_err -> SetMarkerSize(0);
   h_err -> SetFillStyle(3013);

   h_data -> SetLineColor(1);

   if (!NORATIOPLOT) {
      h_data -> Draw("same p e1");
   }

   if (ISLOG == 1) hs -> SetMinimum(110);
   //if(ISLOG==1)    hs->SetMaximum(15000);
   if (ISLOG == 1) hs -> SetMaximum(maxi * 1.35);
   if (ISLOG == 0) hs -> SetMaximum(maxi * 1.35);
   //if(ISLOG==0)   hs->SetMaximum(4600);
   if (ISLOG == 0) hs -> SetMinimum(100);

   double binofwidth = h_mc[0] -> GetBinWidth(1);
   TString binwidth_;
   binwidth_.Form("%1.1f", binofwidth);

   if (!hasNoEvents) {
      hs->GetXaxis()->SetTickLength(0.07);
      hs -> GetXaxis();
      hs -> GetXaxis() -> SetNdivisions(508);
      if (NORATIOPLOT) {
         hs -> GetXaxis() -> SetTitle("XAXISLABEL");
         hs -> GetXaxis() -> SetTitleFont(42);
         hs -> GetXaxis() -> SetLabelFont(42);
         hs -> GetXaxis() -> SetLabelOffset(.01);
         hs -> GetXaxis() -> SetLabelSize(0.03);
         hs -> GetYaxis() -> SetTitle("Events");
         hs -> GetYaxis() -> SetTitleSize(0.05);
         //  hs->GetYaxis()->SetTitleOffset(1);
         hs -> GetYaxis() -> SetTitleFont(42);
         hs -> GetYaxis() -> SetLabelFont(42);
         hs -> GetYaxis() -> SetLabelSize(.03);
         hs -> GetYaxis() -> SetMoreLogLabels();
      } else {
         //hs->GetXaxis()->SetTitle("XAXISLABEL");
         hs -> GetXaxis() -> SetTitleSize(0.00);
         hs -> GetXaxis() -> SetTitleOffset(0.00);
         hs -> GetXaxis() -> SetTitleFont(42);
         hs -> GetXaxis() -> SetLabelFont(42);
         hs -> GetXaxis() -> SetLabelOffset(.01);
         hs -> GetXaxis() -> SetLabelSize(0.04);
         hs -> GetYaxis() -> SetTitle("Events");
         hs -> GetYaxis() -> SetTitleSize(0.045);
         hs -> GetYaxis() -> SetTitleOffset(1);
         hs -> GetYaxis() -> SetTitleFont(42);
         hs -> GetYaxis() -> SetLabelFont(42);
         hs -> GetYaxis() -> SetLabelSize(.03);
         hs -> GetYaxis() -> SetMoreLogLabels();
      }
      hs -> GetXaxis() -> SetRangeUser(XMIN, XMAX);
      hs -> GetXaxis() -> SetNdivisions(508);

      legend -> AddEntry(h_err, "Stat. Unc.", "f");

      //Legend
      TLegend * legendsig;

      legendsig = new TLegend(0.57, 0.5, 0.94, 0.65, NULL, "brNDC");
      legendsig -> SetTextSize(0.030);
      legendsig -> SetBorderSize(0);
      legendsig -> SetLineColor(1);
      legendsig -> SetLineStyle(1);
      legendsig -> SetLineWidth(1);
      legendsig -> SetFillColor(0);
      legendsig -> SetFillStyle(0);
      legendsig -> SetTextFont(42);

      legend -> Draw("same");
      legendsig -> Draw("same");

      //===========================Latex=================//

      TH1F * h_MC_all = new TH1F( * ((TH1F * )(hs -> GetStack() -> Last()))); // To get all MC event count

      //TString latexCMSname= "CMS #it{#bf{Preliminary}}";//"CMS";// #it{#bf{Preliminary}}";
      TString latexCMSname = "";
      TString latexPreCMSname = "#bf{CMS} #it{Preliminary}";

      TString MC_count;
      TString data_count;

      if (ISCUTFLOW) {
         MC_count = "Sel. MC=" + std::to_string(int(h_MC_all -> GetBinContent(h_MC_all -> GetNbinsX())));
         data_count = "; Sel. Data=" + std::to_string(int(h_data -> GetBinContent(h_data -> GetNbinsX())));
      } else {
         MC_count = "MC=" + std::to_string(int(h_MC_all -> Integral()));
         data_count = "; Data=" + std::to_string(int(h_data -> Integral()));
      }

      TString hasQSF = "";

      if (QCDSF == 1) {
         hasQSF = "";
      } else {
         hasQSF = "; w/QCD-SF";
      }

      //TString latexPreCMSname= ""; // "DM+bb: CMS Preliminary: "+MC_count+data_count+hasQSF;

      TString latexnamemiddle;
      latexnamemiddle.Form("%1.1f fb^{-1}", luminosity);
      TString latexnamepost = " (13 TeV)";
      //TString latexname = latexnamepre+latexnamemiddle+latexnamepost;
      TString latexname = latexnamemiddle + latexnamepost;
      TString histolabel;
      //histolabel = "bbMET";

      histolabel = "HISTOLABEL";

      //std::cout <<"HISTOLABEL"<<std::endl;

      TLatex * t2a;
      TLatex * t2b;
      TLatex * t2c;
      TLatex * t2d;

      t2c = new TLatex(0.15, 0.97, latexPreCMSname);
      t2c -> SetTextSize(0.045);

      t2a = new TLatex(0.7, 0.97, latexname);
      t2a -> SetTextSize(0.040);

      t2b = new TLatex(0.22, 0.88, latexCMSname);
      t2b -> SetTextSize(0.03);

      t2d = new TLatex(0.46, 0.9, histolabel);
      t2d -> SetTextSize(0.045);

      t2a -> SetTextAlign(12);
      t2a -> SetNDC(kTRUE);
      t2a -> SetTextFont(42);

      t2b -> SetTextAlign(12);
      t2b -> SetNDC(kTRUE);
      t2b -> SetTextFont(61);

      t2c -> SetTextAlign(12);
      t2c -> SetNDC(kTRUE);
      t2c -> SetTextFont(42);

      t2d -> SetTextAlign(12);
      t2d -> SetNDC(kTRUE);
      t2d -> SetTextFont(42);

      t2a -> Draw("same");
      t2b -> Draw("same");
      t2c -> Draw("same");
      t2d -> Draw("same");

      TH1F * ratiostaterr = (TH1F * ) h_err -> Clone("ratiostaterr");
      ratiostaterr -> Sumw2();
      ratiostaterr -> SetStats(0);
      ratiostaterr -> SetMinimum(0);
      ratiostaterr -> SetMarkerSize(0);
      ratiostaterr -> SetFillColor(kBlack);
      ratiostaterr -> SetFillStyle(3013);

      for (Int_t i = 0; i < h_err -> GetNbinsX() + 2; i++) {
         ratiostaterr -> SetBinContent(i, 1.0);

         if (h_err -> GetBinContent(i) > 1e-6) { //< not empty
            double binerror = h_err -> GetBinError(i) / h_err -> GetBinContent(i);
            ratiostaterr -> SetBinError(i, binerror);
            //     cout << "bin:" <<i << "binerror:" <<h_err->GetBinContent(i)<< "errorband:" << binerror <<std::endl;
         } else {
            ratiostaterr -> SetBinError(i, 999.);

         }

      }

      TH1F * ratiosysterr = (TH1F * ) ratiostaterr -> Clone("ratiosysterr");
      ratiosysterr -> Sumw2();
      ratiosysterr -> SetMarkerSize(0);
      //ratiosysterr->SetFillColor(kYellow-4);
      // final plot
      ratiosysterr -> SetFillColor(kGray);
      //ratiosysterr->SetFillStyle(3002);
      ratiosysterr -> SetFillStyle(1001);

      for (Int_t i = 0; i < h_err -> GetNbinsX() + 2; i++) {
         if (h_err -> GetBinContent(i) > 1e-6) { //< not empty
            double binerror2 = (pow(h_err -> GetBinError(i), 2) +
               pow(0.25 * WJets -> GetBinContent(i), 2) +
               //   pow(0.20 * WJets->GetBinContent(i), 2) +
               pow(0.25 * ZJets -> GetBinContent(i), 2) +
               pow(0.20 * DYJets -> GetBinContent(i), 2) +
               pow(0.25 * TT -> GetBinContent(i), 2) +
               //   pow(0.20 * GJets->GetBinContent(i), 2) +
               pow(0.20 * QCD -> GetBinContent(i), 2) +
               pow(0.25 * STop -> GetBinContent(i), 2) +
               pow(0.20 * DIBOSON -> GetBinContent(i), 2));
            double binerror = sqrt(binerror2);
            ratiosysterr -> SetBinError(i, binerror / h_err -> GetBinContent(i));
         }
      }

      TLegend * ratioleg = new TLegend(0.6, 0.88, 0.89, 0.98);
      //ratioleg->SetFillColor(0);
      ratioleg -> SetLineColor(0);
      ratioleg -> SetShadowColor(0);
      ratioleg -> SetTextFont(42);
      ratioleg -> SetTextSize(0.09);
      ratioleg -> SetBorderSize(1);
      ratioleg -> SetNColumns(2);
      //ratioleg->SetTextSize(0.07);
      //ratioleg->AddEntry(ratiosysterr, "stat + syst", "f");  //Pred. uncert. (stat + syst)
      ratioleg -> AddEntry(ratiostaterr, "stat", "f"); //Pred. uncert. (stat)

      TLegend * ratioleg1 = new TLegend(0.35, 0.35, 0.94, 0.94);
      ratioleg1 -> SetFillColor(0);
      ratioleg1 -> SetLineColor(0);
      ratioleg1 -> SetShadowColor(0);
      ratioleg1 -> SetTextFont(42);
      ratioleg1 -> SetTextSize(0.07);
      ratioleg1 -> SetBorderSize(1);
      ratioleg1 -> SetNColumns(2);
      //ratioleg->SetTextSize(0.07);
      //ratioleg1->AddEntry(ratiosysterr, "Pre-fit", "f");
      //ratioleg1->AddEntry(ratiostaterr, "Postfit", "f");

      //For DATA:
      TH1F * DataMC;
      TH1F * DataMCPre;

      // Lower Tpad Decalaration
      if (!NORATIOPLOT) {
         c12 -> cd();
         DataMC = (TH1F * ) h_data -> Clone();
         DataMCPre = (TH1F * ) h_data -> Clone();
         DataMC -> Divide(Stackhist);
         DataMC -> GetYaxis() -> SetTitle("Data/Pred.");
         DataMC -> GetYaxis() -> SetTitleSize(0.1);
         DataMC -> GetYaxis() -> SetTitleOffset(0.42);
         DataMC -> GetYaxis() -> SetTitleFont(42);
         DataMC -> GetYaxis() -> SetLabelSize(0.08);
         DataMC -> GetYaxis() -> CenterTitle();
         DataMC -> GetXaxis() -> SetTitle("XAXISLABEL");
         DataMC -> GetXaxis() -> SetLabelSize(0.1);
         DataMC -> GetXaxis() -> SetTitleSize(0.1);
         DataMC -> GetXaxis() -> SetTitleOffset(1);
         DataMC -> GetXaxis() -> SetTitleFont(42);
         DataMC -> GetXaxis() -> SetTickLength(0.07);
         DataMC -> GetXaxis() -> SetLabelFont(42);
         DataMC -> GetYaxis() -> SetLabelFont(42);

      }

      TPad * c1_1 = new TPad("c1_1", "newpad", 0, 0.00, 1, 0.3);
      if (!NORATIOPLOT) c1_1 -> Draw();
      c1_1 -> cd();
      c1_1 -> Range(-7.862408, -629.6193, 53.07125, 486.5489);
      c1_1 -> SetFillColor(0);
      c1_1 -> SetTicky(1);
      c1_1 -> SetLeftMargin(0.1);
      c1_1 -> SetRightMargin(0.1);
      c1_1 -> SetTopMargin(0.0); //0.0
      c1_1 -> SetBottomMargin(0.32);
      c1_1 -> SetFrameFillStyle(0);
      c1_1 -> SetFrameBorderMode(0);
      c1_1 -> SetFrameFillStyle(0);
      c1_1 -> SetFrameBorderMode(0);
      c1_1 -> SetLogy(0);

      if (!NORATIOPLOT) {
         if (VARIABLEBINS) {
            c1_1 -> SetLogx(0);
            DataMC -> GetXaxis() -> SetMoreLogLabels();
            DataMC -> GetXaxis() -> SetNoExponent();
            DataMC -> GetXaxis() -> SetNdivisions(508);
         }
         DataMC -> GetXaxis() -> SetRangeUser(XMIN, XMAX);
         DataMC -> SetMarkerSize(0.7);
         DataMC -> SetMarkerStyle(20);
         DataMC -> SetMarkerColor(1);
         DataMCPre -> SetMarkerSize(0.7);
         DataMCPre -> SetMarkerStyle(20);
         DataMCPre -> SetMarkerColor(kRed);
         DataMCPre -> SetLineColor(kRed);

         DataMC -> Draw("P e1");
         //DataMCPre->Draw("P e1 same");
         //ratiosysterr->Draw("e2 same");
         ratiostaterr -> Draw("e2 same");
         DataMC -> Draw("P e1 same");
         //DataMCPre->Draw("P e1 same");

         DataMC -> Draw("P e1 same");
         DataMC -> SetMinimum(-0.2);
         DataMC -> SetMaximum(2.1);
         DataMC -> GetXaxis() -> SetNdivisions(508);
         DataMC -> GetYaxis() -> SetNdivisions(505);
         TLine * line0 = new TLine(XMIN, 1, XMAX, 1);
         line0 -> SetLineStyle(2);
         line0 -> Draw("same");
         c1_1 -> SetGridy();
         ratioleg -> Draw("same");
      }

      if (TEXTINFILE) {

         //=======================================================================
         //Calculating the contribution of each background in particular range
         // As Data DY(ee) diboson TTjets WWJets
         TAxis * xaxis = h_mc[0] -> GetXaxis();
         Int_t binxmin = xaxis -> FindBin(XMIN);
         Int_t binxmax = xaxis -> FindBin(XMAX);

         float dyjets = h_mc[3] -> Integral() + h_mc[4] -> Integral() + h_mc[5] -> Integral() + h_mc[6] -> Integral() + h_mc[7] -> Integral() + h_mc[8] -> Integral() + h_mc[9] -> Integral() + h_mc[10] -> Integral();
         float dyjets_error = TMath::Sqrt(pow(Integral_Error[3], 2) + pow(Integral_Error[4], 2) + pow(Integral_Error[5], 2) + pow(Integral_Error[6], 2) + pow(Integral_Error[7], 2) + pow(Integral_Error[8], 2) + pow(Integral_Error[9], 2) + pow(Integral_Error[10], 2));

         float diboson_ = h_mc[0] -> Integral() + h_mc[1] -> Integral() + h_mc[2] -> Integral();
         float diboson_error = TMath::Sqrt(pow(Integral_Error[0], 2) + pow(Integral_Error[1], 2) + pow(Integral_Error[2], 2));

         float st_ = h_mc[19] -> Integral() + h_mc[20] -> Integral() + h_mc[21] -> Integral() + h_mc[22] -> Integral() + h_mc[22] -> Integral();
         float st_error = TMath::Sqrt(pow(Integral_Error[20], 2) + pow(Integral_Error[21], 2) + pow(Integral_Error[22], 2) + pow(Integral_Error[22], 2) + pow(Integral_Error[19], 2));

         float wjets = h_mc[11] -> Integral() + h_mc[12] -> Integral() + h_mc[13] -> Integral() + h_mc[14] -> Integral() + h_mc[15] -> Integral() + h_mc[16] -> Integral() + h_mc[17] -> Integral() + h_mc[18] -> Integral();
         float wjets_error = TMath::Sqrt(pow(Integral_Error[11], 2) + pow(Integral_Error[12], 2) + pow(Integral_Error[13], 2) + pow(Integral_Error[14], 2) + pow(Integral_Error[15], 2) + pow(Integral_Error[16], 2) + pow(Integral_Error[17], 2) + pow(Integral_Error[18], 2));

         float zjets = h_mc[3] -> Integral() + h_mc[4] -> Integral() + h_mc[5] -> Integral() + h_mc[6] -> Integral() + h_mc[7] -> Integral() + h_mc[8] -> Integral() + h_mc[9] -> Integral() + h_mc[10] -> Integral();
         float zjets_error = TMath::Sqrt(pow(Integral_Error[3], 2) + pow(Integral_Error[4], 2) + pow(Integral_Error[5], 2) + pow(Integral_Error[6], 2) + pow(Integral_Error[7], 2) + pow(Integral_Error[8], 2) + pow(Integral_Error[9], 2) + pow(Integral_Error[10], 2));

         mout << "HISTPATH" << " a b" << std::endl;
         mout << " DATA " << h_data -> Integral() << " 0" << std::endl;
         mout << " DIBOSON " << diboson_ << " " << diboson_error << std::endl;
         mout << " SingleT " << st_ << " " << st_error << std::endl;
         mout << " WJETS " << wjets << " " << wjets_error << std::endl;
         mout << " ZJETS " << zjets << " " << zjets_error << std::endl;
         mout << " DYJETS " << dyjets << " " << dyjets_error << std::endl;

         mout << "========= ======================== =====================" << std::endl;

         // --------------------- table output --------------------
         tableout.precision(3);
         tableout << " Z \\\\rightarrow \\\\nu \\\\nu+Jets & " << dyjets << " \\\\pm " << dyjets_error << "\\\\\\\\" << std::endl;
         tableout << " Z \\\\rightarrow ll + Jets & " << zjets << " \\\\pm " << zjets_error << "\\\\\\\\" << std::endl;
         tableout << " st  & " << st_ << " \\\\pm " << st_error << "\\\\\\\\" << std::endl;
         tableout << " W+Jets & " << wjets << " \\\\pm " << wjets_error << "\\\\\\\\" << std::endl;
         tableout << " WW/WZ/ZZ & " << diboson_ << " \\\\pm " << diboson_error << "\\\\\\\\" << std::endl;

         tableout << " DATA  & " << h_data -> Integral() << std::endl;

         float a = wjets;
         float b = st_;
         //float c = h_data->Integral() - (diboson_ + zh + dyjets);

         tableout << "a " << a << " " << " b " << b << " " << " c " << "c" << std::endl;
         tableout << a << "  " << b << "  " << diboson_ << "  " << zjets << "  " << dyjets << std::endl;
         tableout << " total_bkg " << a + b + diboson_ + zjets + dyjets << std::endl;
         tableout << " " << std::endl;
      }
      c12 -> Draw();

      if (ISLOG == 0) {
         c12 -> SaveAs(DirPreName + dirpathname + "/bbMETPdf/HISTPATH.pdf");
         c12 -> SaveAs(DirPreName + dirpathname + "/bbMETPng/HISTPATH.png");
         rout << "<hr/>" << std::endl;
         rout<<"<table class=\\"\\"> <tr><td><img src=\\""<<"DYPng/HISTPATH.png\\" height=\\"400\\" width=\\"400\\"></td>   </tr> </table>"<<std::endl;

      }

      if (ISLOG == 1) {
         c12 -> SaveAs(DirPreName + dirpathname + "/bbMETPdf/HISTPATH_log.pdf");
         c12 -> SaveAs(DirPreName + dirpathname + "/bbMETPng/HISTPATH_log.png");
         cout << "Saved." << endl;
      }

      fshape -> cd();
      //Save root files for datacards
      Stackhist -> SetNameTitle("bkgSum", "bkgSum");
      Stackhist -> Write();

      DIBOSON -> SetNameTitle("DIBOSON", "DIBOSON");
      DIBOSON -> Write();

      ZJets -> SetNameTitle("ZJets", "ZJets");
      ZJets -> Write();

      QCD -> SetNameTitle("QCD", "QCD");
      QCD -> Write();

      STop -> SetNameTitle("STop", "STop");
      STop -> Write();

      TT -> SetNameTitle("TT", "TT");
      TT -> Write();

      WJets -> SetNameTitle("WJets", "WJets");
      WJets -> Write();

      DYJets -> SetNameTitle("DYJets", "DYJets");
      DYJets -> Write();

      data_obs -> SetNameTitle("data_obs", "data_obs");
      data_obs -> Write();

      fshape -> Write();
      fshape -> Close();

      if (VARIABLEBINS) {
         system("cp " + outputshapefilename + " " + DirPreName + "METBIN_1");
         system("cp " + outputshapefilename + " " + DirPreName + "METBIN_2");
         system("cp " + outputshapefilename + " " + DirPreName + "METBIN_3");
      }
   }
}
'''
## template macro ends here

TemplateOverlapMacro = open('TemplateOverlapMacro.C','w')
TemplateOverlapMacro.write(macro)
TemplateOverlapMacro.close()

def makeplot(inputs):
    QCDSF=1
    if not 'QCD' in inputs[1]:
        if '1b' in inputs[1] or 'sr1' in inputs[1].lower():
            QCDSF=1.2508
        elif '2b' in inputs[1] or 'sr2' in inputs[1].lower():
            QCDSF=0.9458
    #print QCDSF

    TemplateOverlapMacro = open('TemplateOverlapMacro.C','r')
    NewPlot       = open('Plot.C','w')
    for line in TemplateOverlapMacro:
        line = line.replace("HISTPATH",inputs[0])
        line = line.replace("HISTNAME",inputs[1])
        line = line.replace("XAXISLABEL",inputs[2])
        line = line.replace("XMIN",inputs[3])
        line = line.replace("XMAX",inputs[4])
        line = line.replace("REBIN",inputs[5])
        line = line.replace("ISLOG",inputs[6])

        line = line.replace("QCDSF",str(QCDSF))
        line = line.replace("VERBOSE",str(int(verbose)))

        HistName=inputs[1]
        if 'h_reg_' in HistName:
            histolabel=HistName.split('_')[2]
        elif '_sr1' in HistName:
            histolabel="SR1"
        elif '_sr2' in HistName:
            histolabel="SR2"
        else:
            histolabel=""

        line = line.replace("HISTOLABEL",histolabel)


        if len(inputs) > 7 :
            line = line.replace("ISCUTFLOW", inputs[7])
        else :
            line = line.replace("ISCUTFLOW", "0")

        if len(inputs) > 8 :
            line = line.replace("BLINDFACTOR", inputs[8])
        else :
            line = line.replace("BLINDFACTOR", "1")


        if len(inputs) > 9 :
            line = line.replace("NORATIOPLOT", inputs[9])
        else :
            line = line.replace("NORATIOPLOT", "0")
        if len(inputs) >10 :
            line = line.replace("VARIABLEBINS", inputs[10])
        else:
            line = line.replace("VARIABLEBINS", "0")
        if len(inputs) > 11 :
            line = line.replace("TEXTINFILE", inputs[11])
        else :
            line = line.replace("TEXTINFILE", "0")
        if len(inputs) > 12 :
            line = line.replace(".pdf",str(inputs[12]+".pdf"))
            line = line.replace(".png",str(inputs[12]+".png"))
        if len(inputs) > 13:
            line = line.replace("PREFITDIR",inputs[13]) ## PreFitMET or PreFitMass
            line = line.replace("DRAWPREFIT","1") ## 0 or 1
        else:
            line = line.replace("DRAWPREFIT","0")
        #print line
        NewPlot.write(line)
    NewPlot.close()
    os.system('root -l -b -q  Plot.C')
    print

##########Start Adding your plots here


dirnames=['']

emplist=open("Empty.txt","w")
emplist.close()

srblindfactor='1'
srnodata='1'

for dirname in dirnames:

    regions=[]
    PUreg=[]

#    if makeMuCRplots and makeEleCRplots:
#        regions=['2e1b','2mu1b','2e2b','2mu2b','1e1b','1mu1b','1e2b','1mu2b','1mu1e1b','1mu1e2b']
    if makeMuCRplots:
        regions+=['1mu1b']#,'1mutop1b','1mutop2b','2mu1b','2mu2b','1mu2b']
        PUreg+=['mu_']
    if makeEleCRplots:
        regions+=['1e1b']#,'1e2b','1etop1b','1etop2b','2e1b','2e2b']
        PUreg+=['ele_']
    if makePhoCRplots:
        regions+=['1gamma1b','1gamma2b']
        PUreg+=['pho_']
    if makeQCDCRplots:
        regions+=['QCD1b','QCD2b']
        PUreg+=[]
#    else:
#        regions=[]
#        PUreg=[]

    makeplot([dirname+"CRSum",'h_CRSum_','','0.','10.','1','1'])
    if makeMuCRplots: makeplot([dirname+"CRSumMu",'h_CRSumMu_','','0.','6.','1','1'])
    if makeEleCRplots: makeplot([dirname+"CRSumEle",'h_CRSumEle_','','0.','4.','1','1'])

    for dt in PUreg:
        makeplot([dirname+dt+"PuReweightPV",'h_'+dt+'PuReweightPV_','nPV after PU reweighting','0.','50.','1','0'])
        makeplot([dirname+dt+"noPuReweightPV",'h_'+dt+'noPuReweightPV_','nPV before PU reweighting','0.','50.','1','0'])
#        makeplot([dirname+dt+"PuReweightnPVert",'h_'+dt+'PuReweightnPVert_','nPV after PU reweighting (nPVert): '+dt,'0.','100.','100','0'])
#        makeplot([dirname+dt+"noPuReweightnPVert",'h_'+dt+'noPuReweightnPVert_','nPV before PU reweighting (nPVert): '+dt,'0.','100.','100','0'])

# Cutflow plots:
    if makeSRplots:
        makeplot([dirname+"cutflow",'h_cutflow_','Cutflow','0.','10','1','1','1',srblindfactor,srnodata])
        makeplot([dirname+"cutflow_SR1",'h_cutflow_SR1_','SR1 Cutflow','0.','10','1','1','1',srblindfactor,srnodata])
        makeplot([dirname+"cutflow_SR2",'h_cutflow_SR2_','SR2 Cutflow','0.','10','1','1','1',srblindfactor,srnodata])

    for reg in regions:
        makeplot([dirname+"cutflow_"+reg,'h_cutflow_'+reg+'_',reg+' Cutflow','0.','13','1','1','1'])

#Linear plots:
    if makeSRplots:
        '''makeplot([dirname+"jet1_eta_sr1",'h_jet1_eta_sr1_','jet 1 #eta','-3.','3.','1','0','0',srblindfactor,srnodata])
        makeplot([dirname+"jet2_eta_sr1",'h_jet2_eta_sr1_','jet 2 #eta','-3.','3.','1','0','0',srblindfactor,srnodata])

        makeplot([dirname+"jet1_eta_sr2",'h_jet1_eta_sr2_','jet 1 #eta','-3.','3.','1','0','0',srblindfactor,srnodata])
        makeplot([dirname+"jet2_eta_sr2",'h_jet2_eta_sr2_','jet 2 #eta','-3.','3.','1','0','0',srblindfactor,srnodata])
        makeplot([dirname+"jet3_eta_sr2",'h_jet3_eta_sr2_','jet 3 #eta','-3.','3.','1','0','0',srblindfactor,srnodata])

        makeplot([dirname+"jet1_csv_sr1",'h_jet1_csv_sr1_','jet 1 csv','0.','1.','1','0','0',srblindfactor,srnodata])
        makeplot([dirname+"jet2_csv_sr1",'h_jet2_csv_sr1_','jet 2 csv','0.','1.','1','0','0',srblindfactor,srnodata])

        makeplot([dirname+"jet1_deepcsv_sr1",'h_jet1_deepcsv_sr1_','jet 1 deepcsv','0.','1.','1','0','0',srblindfactor,srnodata])
        makeplot([dirname+"jet2_deepcsv_sr1",'h_jet2_deepcsv_sr1_','jet 2 deepcsv','0.','1.','1','0','0',srblindfactor,srnodata])

        #makeplot([dirname+"presel_jet1_csv_sr1",'h_presel_jet1_csv_sr1_','jet 1 csv before selection','0.','1.','1','0','0',srblindfactor,srnodata])
        #makeplot([dirname+"presel_jet2_csv_sr1",'h_presel_jet2_csv_sr1_','jet 2 csv before selection','0.','1.','1','0','0',srblindfactor,srnodata])

        makeplot([dirname+"jet1_csv_sr2",'h_jet1_csv_sr2_','jet 1 csv','0.','1.','1','0','0',srblindfactor,srnodata])
        makeplot([dirname+"jet2_csv_sr2",'h_jet2_csv_sr2_','jet 2 csv','0.','1.','1','0','0',srblindfactor,srnodata])
        makeplot([dirname+"jet3_csv_sr2",'h_jet3_csv_sr2_','jet 3 csv','0.','1.','1','0','0',srblindfactor,srnodata])

        makeplot([dirname+"jet1_deepcsv_sr2",'h_jet1_deepcsv_sr2_','jet 1 deepcsv','0.','1.','1','0','0',srblindfactor,srnodata])
        makeplot([dirname+"jet2_deepcsv_sr2",'h_jet2_deepcsv_sr2_','jet 2 deepcsv','0.','1.','1','0','0',srblindfactor,srnodata])
        makeplot([dirname+"jet3_deepcsv_sr2",'h_jet3_deepcsv_sr2_','jet 3 deepcsv','0.','1.','1','0','0',srblindfactor,srnodata])

        #makeplot([dirname+"presel_jet1_csv_sr2",'h_presel_jet1_csv_sr2_','jet 1 csv before selection','0.','1.','1','0','0',srblindfactor,srnodata])
        #makeplot([dirname+"presel_jet2_csv_sr2",'h_presel_jet2_csv_sr2_','jet 2 csv before selection','0.','1.','1','0','0',srblindfactor,srnodata])
        #makeplot([dirname+"presel_jet3_csv_sr2",'h_presel_jet3_csv_sr2_','jet 3 csv before selection','0.','1.','1','0','0',srblindfactor,srnodata])

        #makeplot([dirname+"presel_jet1_chf_sr1",'h_presel_jet1_chf_sr1_','jet 1 CHadFrac before selection','0.','1.','1','0','0',srblindfactor,srnodata])
        #makeplot([dirname+"presel_jet1_chf_sr2",'h_presel_jet1_chf_sr2_','jet 1 CHadFrac before selection','0.','1.','1','0','0',srblindfactor,srnodata])
        #makeplot([dirname+"presel_jet1_nhf_sr1",'h_presel_jet1_nhf_sr1_','jet 1 NHadFrac before selection','0.','1.','1','0','0',srblindfactor,srnodata])
        #makeplot([dirname+"presel_jet1_nhf_sr2",'h_presel_jet1_nhf_sr2_','jet 1 NHadFrac before selection','0.','1.','1','0','0',srblindfactor,srnodata])

        makeplot([dirname+"jet1_chf_sr1",'h_jet1_chf_sr1_','jet 1 CHadFrac','0.','1.','1','0','0',srblindfactor,srnodata])
        makeplot([dirname+"jet1_chf_sr2",'h_jet1_chf_sr2_','jet 1 CHadFrac','0.','1.','1','0','0',srblindfactor,srnodata])
        makeplot([dirname+"jet1_nhf_sr1",'h_jet1_nhf_sr1_','jet 1 NHadFrac','0.','1.','1','0','0',srblindfactor,srnodata])
        makeplot([dirname+"jet1_nhf_sr2",'h_jet1_nhf_sr2_','jet 1 NHadFrac','0.','1.','1','0','0',srblindfactor,srnodata])'''

    ##for CRs
    for reg in regions:
        if reg[0]=='2': makeplot([dirname+"reg_"+reg+"_Zmass",'h_reg_'+reg+'_Zmass_','Z candidate mass (GeV)','70.','110.',reg[-2],'0'])
        if reg[0]=='1': makeplot([dirname+"reg_"+reg+"_Wmass",'h_reg_'+reg+'_Wmass_','W candidate m_{T} (GeV)','0.','160.','1','0'])
        ''' makeplot([dirname+"reg_"+reg+"_jet1_eta",'h_reg_'+reg+'_jet1_eta_','Lead Jet #eta','-3.5','3.5','1','0'])
        makeplot([dirname+"reg_"+reg+"_jet2_eta",'h_reg_'+reg+'_jet2_eta_','Second Jet #eta','-3.5','3.5','1','0'])

        makeplot([dirname+"reg_"+reg+"_jet1_deepcsv",'h_reg_'+reg+'_jet1_deepcsv_','Lead Jet deepCSV','0.','1.','1','0'])
        makeplot([dirname+"reg_"+reg+"_jet1_csv",'h_reg_'+reg+'_jet1_csv_','Lead Jet CSV','0.','1.','1','0'])
#        makeplot([dirname+"reg_"+reg+"_jet2_csv",'h_reg_'+reg+'_jet2_csv_','Second Jet CSV','0.','1.','1','0'])'''


#Log plots:

###For SR
    if makeSRplots:
        '''makeplot([dirname+"jet1_pT_sr1",'h_jet1_pT_sr1_','jet 1 p_{T} (GeV)','0.','800.','1','1','0',srblindfactor,srnodata])
        makeplot([dirname+"jet2_pT_sr1",'h_jet2_pT_sr1_','jet 2 p_{T} (GeV)','0.','400.','1','1','0',srblindfactor,srnodata])


        makeplot([dirname+"jet1_pT_sr2",'h_jet1_pT_sr2_','jet 1 p_{T} (GeV)','0.','800.','1','1','0',srblindfactor,srnodata])
        makeplot([dirname+"jet2_pT_sr2",'h_jet2_pT_sr2_','jet 2 p_{T} (GeV)','0.','400.','1','1','0',srblindfactor,srnodata])
        makeplot([dirname+"jet3_pT_sr2",'h_jet3_pT_sr2_','jet 3 p_{T} (GeV)','0.','400.','1','1','0',srblindfactor,srnodata])

        makeplot([dirname+"min_dPhi_sr1",'h_min_dPhi_sr1_','min #Delta #phi','0.','3.2','1','1','0',srblindfactor,srnodata])
        makeplot([dirname+"min_dPhi_sr2",'h_min_dPhi_sr2_','min #Delta #phi','0.','3.2','1','1','0',srblindfactor,srnodata])'''

        makeplot([dirname+"met_sr1",'h_met_sr1_','Missing Transverse Energy (GeV)','200.','1000','2','1','0',srblindfactor,srnodata])
        makeplot([dirname+"met_sr2",'h_met_sr2_','Missing Transverse Energy (GeV)','200.','1000','2','1','0',srblindfactor,srnodata])

    # Region based
    for reg in regions:
        if reg[0]=='2': makeplot([dirname+"reg_"+reg+"_ZpT",'h_reg_'+reg+'_ZpT_','Z candidate p_{T} (GeV)','0.','800.',reg[-2],'1'])
        if reg[0]=='1': makeplot([dirname+"reg_"+reg+"_WpT",'h_reg_'+reg+'_WpT_','W candidate p_{T} (GeV)','0.','800.','1','1'])
        makeplot([dirname+"reg_"+reg+"_hadrecoil",'h_reg_'+reg+'_hadrecoil_','Hadronic Recoil (GeV)','0.','1000.','1','1'])

        '''makeplot([dirname+"reg_"+reg+"_jet1_NHadEF",'h_reg_'+reg+'_jet1_NHadEF_','Lead jet neutral hadronic fraction','0.','1.','1','1'])
        makeplot([dirname+"reg_"+reg+"_jet1_CHadEF",'h_reg_'+reg+'_jet1_CHadEF_','Lead jet charged hadronic fraction','0.','1.','1','1'])
        makeplot([dirname+"reg_"+reg+"_jet1_CEmEF",'h_reg_'+reg+'_jet1_CEmEF_','Lead jet charged EM fraction','0.','1.','1','1'])
        makeplot([dirname+"reg_"+reg+"_jet1_PhoEF",'h_reg_'+reg+'_jet1_PhoEF_','Lead jet Photon fraction','0.','1.','1','1'])
        makeplot([dirname+"reg_"+reg+"_jet1_EleEF",'h_reg_'+reg+'_jet1_EleEF_','Lead jet Electron fraction','0.','1.','1','1'])
        makeplot([dirname+"reg_"+reg+"_jet1_MuoEF",'h_reg_'+reg+'_jet1_MuoEF_','Lead jet Muon fraction','0.','1.','1','1'])'''

        if not 'QCD' in reg:
            makeplot([dirname+"reg_"+reg+"_MET",'h_reg_'+reg+'_MET_','Real MET (GeV)','250.','550.','20','1'])
        else:
            makeplot([dirname+"reg_"+reg+"_MET",'h_reg_'+reg+'_MET_','Real MET (GeV)','250.','550.','20','1'])
        makeplot([dirname+"reg_"+reg+"_njet",'h_reg_'+reg+'_njet_','Number of Jets','-1','5','1','1'])

        if reg[:2]=='1e':
            makeplot([dirname+"reg_"+reg+"_njet_n_minus_1",'h_reg_'+reg+'_njet_n_minus_1_','Number of Jets (n-1 cuts plot)','-1','12','1','1'])
            '''makeplot([dirname+"reg_"+reg+"_unclean_njet_n_minus_1",'h_reg_'+reg+'_unclean_njet_n_minus_1_','Number of Jets without cleaning (n-1 cuts plot)','-1','12','1','1'])

            makeplot([dirname+"reg_"+reg+"_min_dR_jet_ele_preclean",'h_reg_'+reg+'_min_dR_jet_ele_preclean_','min dR between jets and electron before jet cleaning','0.','6.','1','1'])
            makeplot([dirname+"reg_"+reg+"_min_dR_jet_ele_postclean",'h_reg_'+reg+'_min_dR_jet_ele_postclean_','min dR between jets and electron after jet cleaning','0.','6.','1','1'])
        makeplot([dirname+"reg_"+reg+"_min_dPhi_jet_Recoil",'h_reg_'+reg+'_min_dPhi_jet_Recoil_','min d #phi between jets and recoil','0.','6.','1','1'])
        makeplot([dirname+"reg_"+reg+"_min_dPhi_jet_MET",'h_reg_'+reg+'_min_dPhi_jet_MET_','min d #phi between jets and real MET','0.','6.','1','1'])

        makeplot([dirname+"reg_"+reg+"_min_dPhi_jet_Recoil_n_minus_1",'h_reg_'+reg+'_min_dPhi_jet_Recoil_n_minus_1_','min d #phi between jets and recoil before d#phi cut','0.','6.','1','1'])'''

        makeplot([dirname+"reg_"+reg+"_ntau",'h_reg_'+reg+'_ntau_','Number of Taus','-1','4','1','1'])
        makeplot([dirname+"reg_"+reg+"_nUncleanTau",'h_reg_'+reg+'_nUncleanTau_','Number of Taus (before cleaning)','-1','6','1','1'])
#            makeplot([dirname+"reg_"+reg+"_ntaucleaned",'h_reg_'+reg+'_ntaucleaned_','Number of Taus (after Tau cleaning)','-1','4','5','1'])
        makeplot([dirname+"reg_"+reg+"_nele",'h_reg_'+reg+'_nele_','Number of Electrons','-1','5','1','1'])
        if reg[0]=='2': makeplot([dirname+"reg_"+reg+"_npho",'h_reg_'+reg+'_npho_','Number of Photons','-1','5','1','1'])
#            makeplot([dirname+"reg_"+reg+"_nUncleanEle",'h_reg_'+reg+'_nUncleanEle_','Number of Eles (before cleaning)','-1','6','7','1'])
        makeplot([dirname+"reg_"+reg+"_nmu",'h_reg_'+reg+'_nmu_','Number of Muons','-1','5','1','1'])
#            makeplot([dirname+"reg_"+reg+"_nUncleanMu",'h_reg_'+reg+'_nUncleanMu_','Number of Muons (before cleaning)','-1','6','7','1'])
        makeplot([dirname+"reg_"+reg+"_lep1_pT",'h_reg_'+reg+'_lep1_pT_','Lead Lepton p_{T} (GeV)','0.','500.','1','1'])
        makeplot([dirname+"reg_"+reg+"_lep2_pT",'h_reg_'+reg+'_lep2_pT_','Second Lepton p_{T} (GeV)','0.','250.','1','1'])
        '''if reg.startswith('1mu1e'): makeplot([dirname+"reg_"+reg+"_e_pT",'h_reg_'+reg+'_e_pT_','Electron p_{T} (GeV)','0.','500.','1','1'])         #Top
        if reg.startswith('1mu1e'): makeplot([dirname+"reg_"+reg+"_mu_pT",'h_reg_'+reg+'_mu_pT_','Muon p_{T} (GeV)','0.','500.','1','1'])           #Top
        makeplot([dirname+"reg_"+reg+"_pho_pT",'h_reg_'+reg+'_pho_pT_','Photon p_{T} (GeV)','0.','500.','1','1'])
        if reg[1]=='m': makeplot([dirname+"reg_"+reg+"_lep1_iso",'h_reg_'+reg+'_lep1_iso_','Lead Lepton isolation','0.','800.','1','1'])
        if reg[1]=='m': makeplot([dirname+"reg_"+reg+"_lep2_iso",'h_reg_'+reg+'_lep2_iso_','Second Lepton isolation','0.','800.','1','1'])
        if reg.startswith('1mu1e'): makeplot([dirname+"reg_"+reg+"_mu_iso",'h_reg_'+reg+'_mu_iso_','Muon isolation','0.','800.','1','1'])            #Top
        makeplot([dirname+"reg_"+reg+"_jet1_pT",'h_reg_'+reg+'_jet1_pT_','Lead Jet p_{T} (GeV)','0.','800.','1','1'])
        makeplot([dirname+"reg_"+reg+"_jet2_pT",'h_reg_'+reg+'_jet2_pT_','Second Jet p_{T} (GeV)','0.','400.','1','1'])
#            makeplot([dirname+"reg_"+reg+"_lep1_dR_tau",'h_reg_'+reg+'_lep1_dR_tau_','dR b/w tau and lead lepton','0.','6.','120','1'])
#            makeplot([dirname+"reg_"+reg+"_lep2_dR_tau",'h_reg_'+reg+'_lep2_dR_tau_','dR b/w tau and second lepton','0.','6.','120','1'])
#            makeplot([dirname+"reg_"+reg+"_min_lep_dR_tau",'h_reg_'+reg+'_min_lep_dR_tau_','minimum dR b/w tau and leptons','0.','6.','120','1'])'''
