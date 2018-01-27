def getAll():
    allquantlist=[]
            
    for region in ['sr']: #,'Zmumucr','Zeecr','Wmucr','Wecr','TOPcr']:                               # Jets: Makes all combinations of region, jet number, etc.
        for jetprop in ['pT','eta','phi','csv']:
            for jetnum in [1,2]:
                allquantlist.append('jet'+str(jetnum)+"_"+jetprop+"_"+region+"1")
                
            for jetnum in [1,2,3]:
                allquantlist.append('jet'+str(jetnum)+"_"+jetprop+"_"+region+"2")                                         

    for nSR in ['1','2']:                                                   # Leptons: Makes all combinations of region, lepton number, etc.
#        for lep in ['mu','el']:
#            props = ['pT','eta','phi']
#            if lep=='mu':
#              props.append('iso')  
#              for lepprop in props:
#                 for nCR in ['1','2']:       # For ZCR, because Z has 2 mu or 2 ele
#                       allquantlist.append(lep+nCR+"_"+lepprop+"_Zmumucr"+nSR)
#                 for region in ['Wecr','Wmucr','TOPcr']:          # For ZCR, only 1 mu and/or ele
#                       allquantlist.append(lep+"1_"+lepprop+"_"+region+nSR)
#            if lep=='el':
#              for lepprop in props:
#                 for nCR in ['1','2']:       # For ZCR, because Z has 2 mu or 2 ele
#                       allquantlist.append(lep+nCR+"_"+lepprop+"_Zeecr"+nSR)
#                 for region in ['Wecr','Wmucr','TOPcr']:          # For ZCR, only 1 mu and/or ele
#                       allquantlist.append(lep+"1_"+lepprop+"_"+region+nSR)        
        
        for quantname in ['met','jet1_nhf','jet1_chf']:
            allquantlist.append(quantname+"_sr"+nSR)

    
    allquantlist+=['min_dPhi_sr1','min_dPhi_sr2'] #,'ZhadronRecoil1mumu','ZhadronRecoil1ee','Zmass1mumu','Zmass1ee','ZpT1mumu','ZpT1ee','WhadronRecoil1mu','WhadronRecoil1e','WpT1mu','WpT1e','Wmass1mu','Wmass1e','WpT1','TOPRecoil1','ZhadronRecoil2mumu','ZhadronRecoil2ee','Zmass2mumu','Zmass2ee','ZpT2mumu','ZpT2ee','WhadronRecoil2mu','WhadronRecoil2e','Wmass2mu','Wmass2e','WpT2mu','WpT2e','TOPRecoil2']
    
    return allquantlist

def getPresel():
    preselquantlist=['presel_jet1_csv_sr1','presel_jet2_csv_sr1','presel_jet1_csv_sr2','presel_jet2_csv_sr2','presel_jet3_csv_sr2','presel_jet1_deepcsv_sr1','presel_jet2_deepcsv_sr1','presel_jet1_deepcsv_sr2','presel_jet2_deepcsv_sr2','presel_jet3_deepcsv_sr2']
    preselquantlist.append('presel_jet1_chf_sr1')
    preselquantlist.append('presel_jet1_chf_sr2')
    preselquantlist.append('presel_jet1_nhf_sr1')
    preselquantlist.append('presel_jet1_nhf_sr2')
    
    return preselquantlist
    
def getRegionQuants():

    regquants=[]
    
    #Z CR
    regions=['2e1b','2mu1b','2e2b','2mu2b']
    varlist=['Zmass','ZpT','hadrecoil','MET','lep1_pT','lep2_pT','lep1_iso','lep2_iso','jet1_pT','jet2_pT','jet1_eta','jet2_eta','jet1_csv','jet2_csv','jet1_deepcsv','jet2_deepcsv','njet','ntau','nele','nmu','nUncleanEle','nUncleanMu','nUncleanTau']#,'lep1_dR_tau','lep2_dR_tau','min_lep_dR_tau','ntaucleaned']
    
    for reg in regions:
        for var in varlist:
            regquants.append("reg_"+reg+"_"+var)
            
    #W CR       
    regions=['1e1b','1mu1b','1e2b','1mu2b']
    varlist=['Wmass','WpT','hadrecoil','MET','lep1_pT','lep1_iso','jet1_pT','jet2_pT','jet1_eta','jet2_eta','jet1_csv','jet2_csv','jet1_deepcsv','jet2_deepcsv','njet','ntau','nele','nmu','nUncleanEle','nUncleanMu','nUncleanTau']
    
    for reg in regions:
        for var in varlist:
            regquants.append("reg_"+reg+"_"+var)
    
    #Top CR       
    regions=['1mu1e1b','1mu1e2b']
    varlist=['hadrecoil','MET','lep1_pT','lep2_pT','lep1_iso','lep2_iso','e_pT','mu_pT','mu_iso','jet1_pT','jet2_pT','jet1_eta','jet2_eta','jet1_csv','jet2_csv','jet1_deepcsv','jet2_deepcsv','njet','ntau','nele','nmu','nUncleanEle','nUncleanMu','nUncleanTau']
    
    for reg in regions:
        for var in varlist:
            regquants.append("reg_"+reg+"_"+var)
    
    return regquants   
    
    
    
