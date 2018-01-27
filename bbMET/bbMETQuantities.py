from ROOT import TFile, TTree, TH1F, TH1D, TH1, TCanvas, TChain,TGraphAsymmErrors, TMath, TH2D, TH2F
import ROOT as ROOT
import AllQuantList

class MonoHbbQuantities:

    def __init__(self, rootfilename):
    
        allquantlist=AllQuantList.getAll()
        preselquantlist=AllQuantList.getPresel()
        regquants=AllQuantList.getRegionQuants()
        
        for quant in allquantlist:
            exec("self."+quant+" = None")  
            exec("self.h_"+quant+" = []")  
            
        for quant in preselquantlist:
            exec("self."+quant+" = None")  
            exec("self.h_"+quant+" = []")
            
        for quant in regquants:
            exec("self."+quant+" = None")  
            exec("self.h_"+quant+" = []")
                   
        self.rootfilename = rootfilename
        #self.allquantities = allquantities
        #self.regime   =  True
        
        self.met      =  -999.0
        self.h_met         =  []
        #self.h_met_rebin         =  []
        
        #self.mass     =  -999.0
        #self.h_mass        =  []
        
#        self.csv1     =  -999.0
#        self.h_csv1        =  []
#        
#        self.csv2     =  -999.0
#        self.h_csv2        =  []
        
#        self.mt              = -999.
        #self.dPhi            = -999.
        self.N_e             = -10
        self.N_mu            = -10
        self.N_tau           = -10
        self.N_Pho           = -10
        self.N_b             = -10
        self.N_j             = -10
      
        self.h_mt              = []
        #self.h_dPhi            = []
        self.h_N_e             = []
        self.h_N_mu            = []
        self.h_N_tau           = []
        self.h_N_Pho           = []
        self.h_N_b             = []
        self.h_N_j             = []
        #self.h_mass            = []
        self.h_met_pdf         = []
        self.h_met_muR         = []
        self.h_met_muF         = []
        
        ## 2d histograms 
        self.h_met_vs_mass     = []
        
        self.weight   = 1.0 

        self.weight_pdf   = []
        self.weight_muR   = []
        self.weight_muF   = []

        self.h_total   = []
        self.h_total_weight   = []
        self.h_npass   = []
        
    def defineHisto(self):
        self.h_total.append(TH1F('h_total','h_total',2,0,2))
        self.h_total_weight.append(TH1F('h_total_weight','h_total_weight',2,0,2))
        self.h_npass.append(TH1F('h_npass','h_nass',2,0,2))
        
#        self.h_cutflow=TH1F('h_cutflow_','h_cutflow_',7, 0, 7)                          # Cutflow     

        self.h_met.append(TH1F('h_met_',  'h_met_',  1000,0.,1000.))
        

        #metbins_ = [200,350,500,1000]
        #self.h_met_rebin.append(TH1F('h_met_rebin_'+postname,  'h_met_rebin'+postname,  3, array(('d'),metbins_)))

        #self.h_mass.append(TH1F('h_mass_'+postname, 'h_mass_'+postname, 400,0.,400.))

        self.h_met_vs_mass.append(TH2F('h_met_vs_mass_', 'h_met_vs_mass_', 1000, 0., 1000., 250, 0, 250.))

#        self.h_csv1.append(TH1F('h_csv1_', 'h_csv1_', 20,0.,1.))
#        self.h_csv2.append(TH1F('h_csv2_', 'h_csv2_', 20,0.,1.))
        #self.h_mt.append(TH1F('h_mt_'+postname,'h_mt_'+postname,100,400.,1400.))
        #self.h_dPhi.append(TH1F('h_dPhi_'+postname,'h_dPhi_'+postname,70, -3.5, 3.5 ))
        self.h_N_e.append(TH1F('h_N_e_','h_N_e_',5,0,5))
        self.h_N_mu.append(TH1F('h_N_mu_','h_N_mu_',5,0,5))
        self.h_N_tau.append(TH1F('h_N_tau_','h_N_tau_',5,0,5))
        self.h_N_Pho.append(TH1F('h_N_Pho_','h_N_Pho_',5,0,5))
        self.h_N_b.append(TH1F('h_N_b_','h_N_b_',10,0,10))
        self.h_N_j.append(TH1F('h_N_j_','h_N_j_',10,0,10))
        
        allquantlist=AllQuantList.getAll()
        preselquantlist=AllQuantList.getPresel()
        regquants=AllQuantList.getRegionQuants()
        
        def getBins(quant):
            if 'eta' in quant:
                bins='30'
                low='-3'
                high='3'
            elif 'dPhi' in quant:
                bins='32'
                low='0'
                high='3.2'
            elif 'phi' in quant:
                bins='64'
                low='-3.2'
                high='3.2'
            elif 'csv' in quant:
                bins='50'
                low='0.'
                high='1.'
            elif 'iso' in quant:
                bins='50'
                low='0.'
                high='0.25'
            elif 'Zmass' in quant:
                bins='60'
                low='70.'
                high='110.'
            elif 'Wmass' in quant:
                bins='32'
                low='0.'
                high='160.'
            elif 'met' in quant:
                bins='20'
                low='0.'
                high='1000.'
            elif 'chf' in quant or 'nhf' in quant:
                bins='40'
                low='0.'
                high='1.'
            elif 'njet' in quant:
                bins='12'
                low='0'
                high='12'
            elif 'ntau' in quant or 'nele' in quant or 'nmu' in quant or 'nUnclean' in quant:
                bins='6'
                low='0'
                high='6'
            elif 'recoil' in quant:
                bins='10'
                low='0.'
                high='1000.'
            elif '_dR_' in quant:
                bins='120'
                low='0.'
                high='6.'
            elif 'lep1_pT' in quant or 'jet2_pT' in quant:                   
                bins='100'
                low='0.'
                high='1000.'
            elif 'lep2_pT' in quant:           
                bins='200'
                low='0.'
                high='1000.'
            else:                   # for pT, mass, etc.
                bins='50'
                low='0.'
                high='1000.'
            return bins,low,high
        
        for quant in allquantlist:
            bins,low,high=getBins(quant)         
            exec("self.h_"+quant+".append(TH1F('h_"+quant+"_','h_"+quant+"_',"+bins+","+low+","+high+"))")
        
        for quant in preselquantlist:
            bins,low,high=getBins(quant)         
            exec("self.h_"+quant+".append(TH1F('h_"+quant+"_','h_"+quant+"_',"+bins+","+low+","+high+"))")
        
        for quant in regquants:
            bins,low,high=getBins(quant)         
            exec("self.h_"+quant+".append(TH1F('h_"+quant+"_','h_"+quant+"_',"+bins+","+low+","+high+"))")
        
        h_met_pdf_tmp = []
        for ipdf in range(2):
            midname = str(ipdf)
            h_met_pdf_tmp.append(TH1F('h_met_pdf'+'_'+midname+'_',  'h_met_pdf',  1000,0.,1000.))
        self.h_met_pdf.append(h_met_pdf_tmp)
        h_met_muR_tmp = []
        for imuR in range(2):
            midname = str(imuR)
            h_met_muR_tmp.append(TH1F('h_met_muR'+'_'+midname+'_',  'h_met_muR',  1000,0.,1000.))
        self.h_met_muR.append(h_met_muR_tmp)
        h_met_muF_tmp = []
        for imuF in range(2):
            midname = str(imuF)
            h_met_muF_tmp.append(TH1F('h_met_muF'+'_'+midname+'_',  'h_met_muF',  1000,0.,1000.))
        self.h_met_muF.append(h_met_muF_tmp)

        print "Histograms defined"
        
        
    def FillPreSel(self):
        WF = self.weight
        
        preselquantlist=AllQuantList.getPresel()
        for quant in preselquantlist:
            exec("if self."+quant+" is not None: self.h_"+quant+"[0] .Fill(self."+quant+", WF)")
    
    def FillRegionHisto(self):
        WF = self.weight
        
        regquants=AllQuantList.getRegionQuants()
        for quant in regquants:
            exec("if self."+quant+" is not None: self.h_"+quant+"[0] .Fill(self."+quant+", WF)")
        
    def FillHisto(self):
        WF = self.weight
        #print "WF = ", WF
        self.h_met[0]        .Fill(self.met,       WF)
        
        
        for ipdf in range(2):
            self.h_met_pdf[0]        [ipdf].Fill(self.met,       1.0)

        for imuR in range(2):
            self.h_met_muR[0]        [imuR].Fill(self.met,       1.0)
            
        for imuF in range(2):
            self.h_met_muF[0]        [imuF].Fill(self.met,       1.0)
        

        #self.h_met_vs_mass[0] .Fill(self.met, self.mass, WF)

        #self.h_mass           Fill(self.mass,      WF)
        #self.h_csv1           .Fill(self.csv1,      WF)
        #self.h_csv2           .Fill(self.csv2,      WF)
        #self.h_mt             .Fill(self.mt,        WF)
        #self.h_dPhi           Fill(self.dPhi,      WF)
        self.h_N_e[0]            .Fill(self.N_e,       WF)
        self.h_N_mu[0]           .Fill(self.N_mu,      WF)
        self.h_N_tau[0]          .Fill(self.N_tau,     WF)
        self.h_N_Pho[0]          .Fill(self.N_Pho,     WF)
        self.h_N_b[0]            .Fill(self.N_b,       WF)
        self.h_N_j[0]            .Fill(self.N_j,       WF)
#        print len(self.h_jet1_pT_sr1)
#        print "HbbQuants: "+str(self.jet1_pT_sr2)

        allquantlist=AllQuantList.getAll()
        for quant in allquantlist:
            exec("if self."+quant+" is not None: self.h_"+quant+"[0] .Fill(self."+quant+", WF)")
        
       
    def WriteHisto(self, (nevts,nevts_weight,npass,cutflowvalues,cutflownames,cutflowvaluesSR1,cutflownamesSR1,cutflowvaluesSR2,cutflownamesSR2,CRvalues,CRnames,regionnames,CRcutnames,CRcutflowvaluesSet)):
        f = TFile(self.rootfilename,'RECREATE')
        print 
        f.cd()
        self.h_total[0].SetBinContent(1,nevts)
        self.h_total[0].Write()
        
        self.h_total_weight[0].SetBinContent(1,nevts_weight)
        self.h_total_weight[0].Write()
        
        self.h_npass[0].SetBinContent(1,npass)
        self.h_npass[0].Write()
        
        ncutflow=len(cutflowvalues)
        self.h_cutflow=TH1F('h_cutflow_','h_cutflow_',ncutflow, 0, ncutflow)                          # Cutflow         
        for icutflow in range(len(cutflowvalues)):
            self.h_cutflow.GetXaxis().SetBinLabel(icutflow+1,cutflownames[icutflow])
            self.h_cutflow.SetBinContent(icutflow+1,cutflowvalues[icutflow])
        self.h_cutflow.Write()
        
        ncutflowSR1=len(cutflowvaluesSR1)
        self.h_cutflowSR1=TH1F('h_cutflow_SR1_','h_cutflow_SR1_',ncutflowSR1, 0, ncutflowSR1)                          # Cutflow         
        for icutflow in range(len(cutflowvaluesSR1)):
            self.h_cutflowSR1.GetXaxis().SetBinLabel(icutflow+1,cutflownamesSR1[icutflow])
            self.h_cutflowSR1.SetBinContent(icutflow+1,cutflowvaluesSR1[icutflow])
        self.h_cutflowSR1.Write()
        
        ncutflowSR2=len(cutflowvaluesSR2)
        self.h_cutflowSR2=TH1F('h_cutflow_SR2_','h_cutflow_SR2_',ncutflowSR2, 0, ncutflowSR2)                          # Cutflow         
        for icutflow in range(len(cutflowvaluesSR2)):
            self.h_cutflowSR2.GetXaxis().SetBinLabel(icutflow+1,cutflownamesSR2[icutflow])
            self.h_cutflowSR2.SetBinContent(icutflow+1,cutflowvaluesSR2[icutflow])
        self.h_cutflowSR2.Write()
        
        for ireg in range(len(regionnames)):
            ncutflowCR=len(CRcutflowvaluesSet[ireg])
            self.h_cutflowCR=TH1F('h_cutflow_'+regionnames[ireg]+'_','h_cutflow_'+regionnames[ireg]+'_',ncutflowCR, 0, ncutflowCR)  
            for icutflow in range(len(CRcutflowvaluesSet[ireg])):
                self.h_cutflowCR.GetXaxis().SetBinLabel(icutflow+1,CRcutnames[icutflow])
                self.h_cutflowCR.SetBinContent(icutflow+1,CRcutflowvaluesSet[ireg][icutflow])
            self.h_cutflowCR.Write()
        
        nCR=len(CRvalues)
        self.h_CRs=TH1F('h_CRs_','h_CRs_',nCR, 0, nCR)                          # CR flow         
        for iCR in range(nCR):
            self.h_CRs.GetXaxis().SetBinLabel(iCR+1,CRnames[iCR])
            self.h_CRs.SetBinContent(iCR+1,CRvalues[iCR])
        self.h_CRs.Write()
        
        self.h_met[0].Write()
        #self.h_met_rebin[iregime].Write()
        for ipdf in range(2):
            self.h_met_pdf[0][ipdf].Write()
        for imuR in range(2):
            self.h_met_muR[0][imuR].Write()
        for imuF in range(2):
            self.h_met_muF[0][imuF].Write()

        #self.h_met_vs_mass.Write()

        #self.h_mass.Write()
        #self.h_csv1.Write()
        #self.h_csv2.Write()
        #self.h_mt.Write()
        #self.h_dPhi.Write()
        self.h_N_e[0].Write()
        self.h_N_mu[0].Write()
        self.h_N_tau[0].Write()
        self.h_N_Pho[0].Write()
        self.h_N_b[0].Write()
        self.h_N_j[0].Write()
        #self.h_mass.Write()
        
        allquantlist=AllQuantList.getAll()
        preselquantlist=AllQuantList.getPresel()
        regquants=AllQuantList.getRegionQuants()
        
        for quant in allquantlist:
            exec("self.h_"+quant+"[0].Write()")
            
        for quant in preselquantlist:
            exec("self.h_"+quant+"[0].Write()")
            
        for quant in regquants:
            exec("self.h_"+quant+"[0].Write()")
            
