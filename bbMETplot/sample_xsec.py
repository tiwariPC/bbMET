import os, sys
#//--------------------------------------------------------------------------------------
def getXsec(samplename):
    Xsec = 1.0
    #print (samplename)
    if 'DYJetsToLL_M-50_HT-70to100'   in samplename:   Xsec  = 169.9
    elif 'DYJetsToLL_M-50_HT-100to200'   in samplename:   Xsec  = 147.4
    elif 'DYJetsToLL_M-50_HT-200to400' in samplename:   Xsec  = 40.99
    elif 'DYJetsToLL_M-50_HT-400to600'   in samplename: Xsec  = 5.678
    elif 'DYJetsToLL_M-50_HT-600to800'   in samplename: Xsec  = 1.367
    elif 'DYJetsToLL_M-50_HT-800to1200'  in samplename: Xsec  = 0.6304
    elif 'DYJetsToLL_M-50_HT-1200to2500' in samplename: Xsec  = 0.1514
    elif 'DYJetsToLL_M-50_HT-2500toInf'  in samplename: Xsec  = 0.003565
    elif 'ZJetsToNuNu_HT-100To200'   in samplename:     Xsec  = 280.35
    elif 'ZJetsToNuNu_HT-200To400'   in samplename:     Xsec  = 77.67
    elif 'ZJetsToNuNu_HT-400To600'   in samplename:     Xsec  = 10.73
    elif 'ZJetsToNuNu_HT-600To800'   in samplename:     Xsec  = 2.559
    elif 'ZJetsToNuNu_HT-800To1200'  in samplename:     Xsec  = 1.1796
    elif 'ZJetsToNuNu_HT-1200To2500' in samplename:     Xsec  = 0.28833
    elif 'ZJetsToNuNu_HT-2500ToInf'  in samplename:     Xsec  = 0.006945
    elif 'WJetsToLNu_HT-70To100'    in samplename: Xsec  = 1372.0
    elif 'WJetsToLNu_HT-100To200'   in samplename: Xsec  = 1345.0
    elif 'WJetsToLNu_HT-200To400'   in samplename: Xsec  = 359.7
    elif 'WJetsToLNu_HT-400To600'   in samplename: Xsec  = 48.91
    elif 'WJetsToLNu_HT-600To800'   in samplename: Xsec  = 12.05
    elif 'WJetsToLNu_HT-800To1200'  in samplename: Xsec  = 5.501
    elif 'WJetsToLNu_HT-1200To2500' in samplename: Xsec  = 1.329
    elif 'WJetsToLNu_HT-2500ToInf'  in samplename: Xsec  = 0.03216
    elif 'GJets_HT-40To100'    in samplename: Xsec  = 20790
    elif 'GJets_HT-100To200'   in samplename: Xsec  = 9238
    elif 'GJets_HT-200To400'   in samplename: Xsec  = 2305
    elif 'GJets_HT-400To600'   in samplename: Xsec  = 274.4
    elif 'GJets_HT-600ToInf'   in samplename: Xsec  = 93.46
    elif 'QCD_HT500to700'    in samplename: Xsec  = 32100
    elif 'QCD_HT700to1000'   in samplename: Xsec  = 6831
    elif 'QCD_HT1000to1500'  in samplename: Xsec  = 1207
    elif 'QCD_HT1500to2000'  in samplename: Xsec  = 119.9
    elif 'QCD_HT2000toInf'   in samplename: Xsec  = 25.24
    elif 'TT_TuneCUETP8M2T4' in samplename: Xsec = 831.76
    elif 'TTToSemilepton'    in samplename: Xsec = 364.35
    elif 'TTTo2L2Nu'         in samplename: Xsec = 87.31
    elif 'WWTo1L1Nu2Q_13TeV' in samplename: Xsec = 49.997
    elif 'WWTo2L2Nu_13TeV'   in samplename: Xsec = 12.178
    elif 'WWTo4Q_4f_13TeV'   in samplename: Xsec = 51.723
    elif 'WZTo1L1Nu2Q_13TeV' in samplename: Xsec = 10.71
    elif 'WZTo1L3Nu_13TeV'   in samplename: Xsec = 3.0330
    elif 'WZTo2L2Q_13TeV'    in samplename: Xsec = 5.5950
    elif 'WZTo2Q2Nu_13TeV'   in samplename: Xsec = 6.4880
    elif 'WZTo3LNu'          in samplename: Xsec = 4.4297
    elif 'ZZTo2L2Q_13TeV'    in samplename: Xsec = 3.22
    elif 'ZZTo2Q2Nu_13TeV'   in samplename: Xsec = 4.04
    elif 'ZZTo4L_13TeV'      in samplename: Xsec = 1.2120
    elif 'ZZTo4Q_13TeV'      in samplename: Xsec = 6.842
    elif 'ST_s-channel_4f_leptonDecays' in samplename: Xsec = 3.36
    elif 'ST_t-channel_antitop_4f'      in samplename: Xsec = 80.95
    elif 'ST_t-channel_top_4f'          in samplename: Xsec = 136.02
    elif 'ST_tW_antitop_5f'             in samplename: Xsec = 35.85
    elif 'ST_tW_top_5f'                 in samplename: Xsec = 35.85

    return Xsec
