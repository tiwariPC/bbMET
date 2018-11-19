import sys,os,array
if len(sys.argv)!=3:
    print "Usage: python CombinedRootMaker.py directorypath syst_path"
    sys.exit()

from ROOT import *
gROOT.SetBatch(True)

fold=sys.argv[1].strip('/')
fold_syst=sys.argv[2].strip('/')

def setHistStyle(h_temp2,bins,newname):
    h_temp=h_temp2.Rebin(len(bins)-1,"h_temp",array.array('d',bins))
    h_temp.SetName(newname)
    h_temp.SetTitle(newname)
    h_temp.SetLineWidth(1)
    h_temp.SetBinContent(len(bins)-1,h_temp.GetBinContent(len(bins)-1)+h_temp.GetBinContent(len(bins))) #Add overflow bin content to last bin
    h_temp.SetBinContent(len(bins),0.)
    h_temp.GetXaxis().SetRangeUser(200,1000)
    h_temp.SetMarkerColor(kBlack);
    h_temp.SetMarkerStyle(2);
    return h_temp

CRnames=['Topenu','Topmunu','Wenu','Wmunu','Zee','Zmumu','Gamma','SR','signal']

fdict={}
for CR in CRnames:
    for reg in ['1b','2b']:
        fdict[CR+"_"+reg]=TFile("DataCardRootFiles/bbDM_2016_"+CR+"_"+reg+".root","RECREATE")

f=TFile("DataCardRootFiles/AllMETHistos.root","RECREATE")
f.cd()

inCRfiles=[fold+"/"+i for i in os.listdir(fold) if i.endswith('hadrecoil.root')]
inSRfiles=[fold+"/"+i for i in ['met_sr1.root','met_sr2.root']]
inSystfiles=[fold_syst+"/"+i for i in os.listdir(fold_syst) if "syst" in i]

bins=[200.,250.,300.,400.,500.,700.,1000.,2000.]

def getCRcat(infile):

    flnamesplit=infile.split('/')[-1].split('.')[0].split('_')
    if 'sr1' in infile:
        CR="SR"
        category='1b'
    elif 'sr2' in infile:
        CR="SR"
        category='2b'
    else:
        reg=flnamesplit[1]
        category=reg[-2:]
        if reg.startswith('1mutop'):
            CR='Topmunu'
        if reg.startswith('1etop'):
            CR='Topenu'
        elif reg.startswith('1e') and not reg.startswith('1etop'):
            CR='Wenu'
        elif reg.startswith('1mu') and not reg.startswith('1mutop'):
            CR='Wmunu'
        elif reg.startswith('2e'):
            CR='Zee'
        elif reg.startswith('2mu'):
            CR='Zmumu'
        elif reg.startswith('1gamma'):
            CR='Gamma'
    return CR,category

#def getSortKey(infile):
#    CR,category=getCRcat(infile)
#    newname="bbDM_2016_"+CR+"_"+sampname+"_"+category
#    flnamesplit=infile.split('/')[-1].split('.')[0].split('_')
#    if 'syst' in infile:
#        newname += "_"+flnamesplit[2]+"_"+flnamesplit[4]
#    return newname

SRCRhistos=[]

for infile in sorted(inCRfiles+inSRfiles+inSystfiles):
    fin=TFile(infile,"READ")

    flnamesplit=infile.split('/')[-1].split('.')[0].split('_')

    samplist=['DIBOSON','ZJets','GJets','STop','TT','WJets','DYJets','QCD']

    if 'sr' in infile or 'syst' in infile:
        samplist.append('bkgSum')
    else:
        samplist.append('data_obs')
    #try:
    CR,category=getCRcat(infile)
    #except:
        #print infile
    h_tot=fin.Get('bkgSum')


    print
    if not 'syst' in infile:
    #    print infile
        print "Region: "+CR+category
        try:
            tot=h_tot.Integral()
        except Exception as e:
            print e
            print "WARNING: Skipping for "+CR+category+".******************************************************************************"
            continue
        print "Total = "+str(tot)
    else:
#        if "met" in flnamesplit[2] : continue               ####Since MET is not added yet ## TEMPORARY
        print "Region: "+CR+category+" "+flnamesplit[2]+" "+flnamesplit[4]
        try:
            tot=h_tot.Integral()
        except Exception as e:
            print e
            print "WARNING: Skipping for "+CR+category+" "+flnamesplit[2]+" "+flnamesplit[4]+".******************************************************************************"
            continue

    for samp in samplist:
        h_temp2=fin.Get(samp)

        if (samp=='bkgSum' or samp=='data_obs') and not 'syst' in infile:
            sampname='data_obs'
        elif samp=='bkgSum' and 'syst' in infile:
            sampname='bkgSum'
        else:
            sampname=samp

        sampname=sampname.replace("TT","Top")

        newname=CR+"_"+category+"_"+sampname

        shortname=sampname

        if 'syst' in infile:
            newname += "_"+flnamesplit[2]+"_"+flnamesplit[4]
            shortname += "_"+flnamesplit[2]+"_"+flnamesplit[4]

        h_temp=setHistStyle(h_temp2,bins,shortname)
        #if CR=="SR": h_temp.Scale(20)
        sel=h_temp.Integral()
        fdict[CR+"_"+category].cd()
        h_temp.Write()

        h_temp=setHistStyle(h_temp2,bins,newname)
        #if CR=="SR": h_temp.Scale(20)
        f.cd()
        h_temp.Write()

#        except:
#            sel=0.
#            print "Skipped "+infile+" "+samp
        if not 'syst' in infile:
            if tot!=0:
                frac=sel/tot
            else:
                frac=0.
            if samp!="data_obs" and samp!="bkgSum": print "    Sample = " + samp+": Count = %.2f, Fraction = %.4f"%(sel,frac)
    fin.Close()


##Signal
CrossSec_2HDMa={'Mphi_100':6.65229,'Mphi_350':0.38053,'Mphi_400':0.23712,'Mphi_500':0.17824}
def getCross_2HDMa(file):
    ma=file.split('/')[-1].split('.')[0].split('_')[4].split('-')[-1]
    XSec=CrossSec_2HDMa['Mphi_'+str(ma)]
    return XSec,ma

CrossSec_DMsimp={'Mphi_100':0.02072,'Mphi_350':0.001237,'Mphi_400':0.0005773,'Mphi_500':0.0002802}
def getCross_DMsimp(file):
    ma=file.split('/')[-1].split('.')[0].split('_')[4].split('-')[-1]
    XSec=CrossSec_DMsimp['Mphi_'+str(ma)]
    return XSec,ma

LOFiles=[fold+"/LO/"+i for i in os.listdir(fold+"/LO/") if i.endswith('.root')]
NLOFiles=[fold+"/NLO/"+i for i in os.listdir(fold+"/NLO/") if i.endswith('.root')]
ttDMFiles=[fold+"/ttDM/"+i for i in os.listdir(fold+"/ttDM/") if i.endswith('.root')]

regions=['2e1b','2mu1b','2e2b','2mu2b','1e1b','1mu1b','1e2b','1mu2b','1mutop1b','1mutop2b','1etop1b','1etop2b']

lumi=35900.

for infile in LOFiles+NLOFiles+ttDMFiles:
    fin=TFile(infile,"READ")
    h_total=fin.Get('h_total_weight')
    tot=h_total.Integral()

    Mchi=''
    Mphi=''

    for partname in infile.split('/')[-1].split('.')[0].split('_'):
        if partname.startswith('Mchi'): Mchi=partname.replace('-','_')
        if partname.startswith('Mphi'): Mphi=partname.replace('-','_')
    print "Total = "+str(tot)
    print ('Mchi ', Mchi,'  Mphi ', Mphi)
    if 'ttDM' in infile:
        samp='tt'
    elif 'NLO' in infile:
        samp='bbNLO'
    else:
        samp='bbLO'

    if "pseudo" in infile:
        samp+="_pseudo"
    else:
        samp+="_scalar"
    if (Mphi=='Mphi_50' or Mphi=='Mphi_1000' or Mphi=='Mphi_300' or Mphi=='Mphi_750'): continue
    if ( Mchi != 'Mchi_1'): continue
    print ('file ', infile)
    xsec_2HDMa,MH4=getCross_2HDMa(infile)
    xsec_2HDMa = float(xsec_2HDMa)

    xsec_DMsimp,MH4=getCross_DMsimp(infile)
    xsec_DMsimp = float(xsec_DMsimp)
    # if ('Mphi-50.' in infile or 'Mphi-1000.' in infile or 'Mphi-300.'  in infile) and ( 'Mchi-1_' not in infile) : continue
    #     xsec_2HDMa = 1.
    #     xsec_DMsimp = 1.

    Xlist=[xsec_DMsimp,xsec_2HDMa]

    # Store SR1 and SR2 in the signal_*b file
    for sr,category in [['sr1','1b'],['sr2','2b']]:
        h_temp2=fin.Get('h_met_'+sr+'_')
        for i in range(len(Xlist)):
            if i==0:
                newname=samp+"_"+category+"_"+Mchi+"_"+Mphi
                h_temp=setHistStyle(h_temp2,bins,newname)
                sel=h_temp.Integral()
                h_temp.Scale(lumi*Xlist[i]/tot)

                fdict['signal_'+category].cd()
                h_temp.Write()

                f.cd()
                h_temp.Write()
            elif i==1:
                newname=samp+"_"+category+"_"+Mchi+"_"+Mphi+"_2HDMa"
                h_temp=setHistStyle(h_temp2,bins,newname)
                sel=h_temp.Integral()
                h_temp.Scale(lumi*Xlist[i]/tot)

                fdict['signal_'+category].cd()
                h_temp.Write()

                f.cd()
                h_temp.Write()
        if samp.startswith("bb"):
            print (samp)
            #h_temp=setHistStyle(h_temp2,bins,samp[2:]+"_mphi_"+Mphi.split('-')[1]+"_mchi_"+Mchi.split('-')[1])
            h_temp=setHistStyle(h_temp2,bins,samp[2:]+"_mphi_"+Mphi.split('_')[1]+"_mchi_"+Mchi.split('_')[1])
            h_temp.Scale(lumi/tot)

            fdict['SR_'+category].cd()
            h_temp.Write()

        if samp=="bbNLO_scalar" and Mchi=="Mchi_1": #and Mphi=="Mphi_300":
            h_temp=setHistStyle(h_temp2,bins,"sig")
            h_temp.Scale(lumi/tot)

            fdict['SR_'+category].cd()
            h_temp.Write()

        print "    "+sr.upper()+": Count %.2f, Sel Eff. = %.8f"%(sel,sel/tot)

    # Store all CR

    for reg in regions:
        h_temp2=fin.Get('h_reg_'+reg+'_'+'hadrecoil_')

        if samp.startswith("bb"):
            pseudofile=fold+"/"+"reg_"+reg+"_hadrecoil.root"
            CR,category=getCRcat(pseudofile)
            h_temp=setHistStyle(h_temp2,bins,samp[2:]+"_mphi_"+Mphi.split('_')[1]+"_mchi_"+Mchi.split('_')[1])
            h_temp.Scale(lumi/tot)
            if h_temp.Integral()==0.:
                for ibin in range(h_temp.GetSize()-1):
                    h_temp.SetBinContent(ibin,6e-4)
            if h_temp.Integral()<0.:
                h_temp.Scale(6e-4/h_temp.Integral())

            fdict[CR+'_'+category].cd()
            h_temp.Write()

        #if samp.startswith("bb"):
            #h_temp=setHistStyle(h_temp2,bins,samp[2:]+"_mphi_"+Mphi.split('-')[1]+"_mchi_"+Mchi.split('-')[1])
            #h_temp.Scale(lumi/tot)

            #for CR in CRnames:
                #fdict[CR+'_'+category].cd()
                #h_temp.Write()

        #if samp=="bbNLO_scalar" and Mchi=="Mchi-1" and Mphi=="Mphi-300":
            #h_temp=setHistStyle(h_temp2,bins,"sig")
            #h_temp.Scale(lumi/tot)

            #for CR in CRnames:
                #fdict[CR+'_'+category].cd()
                #h_temp.Write()



for CR in CRnames:
    for reg in ['1b','2b']:
        fdict[CR+"_"+reg].Close()
f.Close()
