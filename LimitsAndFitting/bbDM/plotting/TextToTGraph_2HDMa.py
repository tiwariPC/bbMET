from ROOT import TGraph, TFile, TGraphAsymmErrors
from array import array
import os


inputstring_ = ['2HDM'] #['scalar','pseudo']


xsec_dict={}
for iline in open('xsec2HDM.txt'):
    masspoint = str(iline.rstrip().split()[0])
    xsec = str(iline.rstrip().split()[1])
    key_ = masspoint
    xsec_dict[key_] = xsec
print xsec_dict

for inputstring in inputstring_:
    filename = '../bin/limits_bbDM2016'+inputstring+'.txt'
    f = open(filename,"r") 
    med=array('f')
    mchi=array('f')
    expm2=array('f')
    expm1=array('f')
    expmed=array('f')
    expp1=array('f')
    expp2=array('f')
    obs=array('f')
    errx=array('f')
    for line in f:
        if line.rstrip().split()[0]=='400': continue
#        masspointstr = 'NLO_'+inputstring+'_mphi_'+line.rstrip().split()[0]+'_mchi_'+line.rstrip().split()[1]
        masspointstr = 'NLO_pseudo_mphi_'+line.rstrip().split()[0]+'_mchi_'+line.rstrip().split()[1]
        xsec_ = float(xsec_dict[masspointstr] )
        
        med.append(float(line.rstrip().split()[0]))
        mchi.append(float(line.rstrip().split()[1]))
        expm2.append(float(line.rstrip().split()[4])/xsec_-float(line.rstrip().split()[2])/xsec_)
        expm1.append(float(line.rstrip().split()[4])/xsec_-float(line.rstrip().split()[3])/xsec_)
        expmed.append(float(line.rstrip().split()[4])/xsec_)
        expp1.append(float(line.rstrip().split()[5])/xsec_-float(line.rstrip().split()[4])/xsec_)
        expp2.append(float(line.rstrip().split()[6])/xsec_-float(line.rstrip().split()[4])/xsec_)
        obs.append(float(line.rstrip().split()[7])/xsec_)
        errx.append(0.0)
        print masspointstr, float(line.rstrip().split()[3]), xsec_, float(line.rstrip().split()[3])/xsec_
        
    g_exp2  = TGraphAsymmErrors(int(len(med)), med, expmed, errx, errx, expm2, expp2 )   ;  g_exp2.SetName("exp2")
    g_exp1  = TGraphAsymmErrors(int(len(med)), med, expmed, errx, errx, expm1, expp1 )   ;  g_exp1.SetName("exp1")
    g_expmed = TGraphAsymmErrors(int(len(med)), med, expmed)   ;  g_expmed.SetName("expmed")
    g_obs    = TGraphAsymmErrors(int(len(med)), med, obs   )   ;  g_obs.SetName("obs")


#g_expm2  = TGraphAsymmErrors(int(len(med)), med, expm2 )   ;  g_expm2.SetName("expm2")
#g_expm1  = TGraphAsymmErrors(int(len(med)), med, expm1 )   ;  g_expm1.SetName("expm1")
#g_expmed = TGraphAsymmErrors(int(len(med)), med, expmed)   ;  g_expmed.SetName("expmed")
#g_expp1  = TGraphAsymmErrors(int(len(med)), med, expp1 )   ;  g_expp1.SetName("expp1")
#g_expp2  = TGraphAsymmErrors(int(len(med)), med, expp2 )   ;  g_expp2.SetName("expp2")
#g_obs    = TGraphAsymmErrors(int(len(med)), med, obs   )   ;  g_obs.SetName("obs")
#

    f1 = TFile(str(inputstring)+'.root','RECREATE')
    #g_expm2.Write() 
    #g_expm1.Write() 
    #g_expmed.Write()
    #g_expp1.Write() 
    #g_expp2.Write() 
    #g_obs.Write() 

    g_exp2.Write()
    g_exp1.Write()
    g_expmed.Write()
    g_obs.Write()
    
    f1.Write()
    f1.Close()



