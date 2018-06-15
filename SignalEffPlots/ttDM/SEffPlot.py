import os
from ROOT import *
import matplotlib.pyplot as plt, numpy as np
from matplotlib import text
from matplotlib.colors import LogNorm

dirpath="signal/"

chilist=[1,10,50,100,350,450]
philist=[10,50,100,200,300,350,400,500,750,1000]

scalarSEmatrixSR1=np.zeros(shape=(len(chilist),len(philist)))
pseudoscalarSEmatrixSR1=np.zeros(shape=(len(chilist),len(philist)))
scalarSEmatrixSR2=np.zeros(shape=(len(chilist),len(philist)))
pseudoscalarSEmatrixSR2=np.zeros(shape=(len(chilist),len(philist)))


def getFirstLastBin(h):
    first=h.GetBinContent(1)
    last=first
    i=2
    while h.GetBinContent(i) > 0.:        
        last=h.GetBinContent(i)
        i += 1
    return first, last

for med in ['scalar','pseudoscalar']:
    for ichi in range(len(chilist)):
        for iphi in range(len(philist)):
            filename=os.path.join(dirpath,"Output_"+med+"_Mchi-"+str(chilist[ichi])+"_Mphi-"+str(philist[iphi])+"_TuneCUETP8M1.root")
            print filename
            if os.path.isfile(filename):
                rootfile=TFile.Open(filename,'READ')
                
#                h_t = TH1F('h_t','h_t',2,0,2)
                #h_cutflow_SR1=rootfile.Get('h_cutflow_SR1_')
                #h_cutflow_SR2=rootfile.Get('h_cutflow_SR2_')               
                    
                #tot1,sel1 = getFirstLastBin(h_cutflow_SR1)
                #tot2,sel2 = getFirstLastBin(h_cutflow_SR2)
                
                #h_total_weight=rootfile.Get('h_total_weight')
                #tot1=h_total_weight.GetBinContent(1)
                #tot2=tot1
                
                sel1=rootfile.Get('h_met_sr1_').Integral()
                sel2=rootfile.Get('h_met_sr2_').Integral()
                
                h_total_weight=rootfile.Get('h_total_weight')
                tot1=h_total_weight.Integral()
                tot2=tot1
                
                print filename + ": Total: " + str(int(tot1))+ "; SR1: " + str(int(sel1))+ "; SR2: " + str(int(sel2))
                
                exec(med+"SEmatrixSR1[ichi][iphi]=sel1/tot1")
                exec(med+"SEmatrixSR2[ichi][iphi]=sel2/tot2")
                
for med in ['scalar','pseudoscalar']:
    for reg in ['1','2']:              
#    exec("for row in "+med+"SEmatrix: print row")
#    print 
        exec("grid="+med+"SEmatrixSR"+reg)
#Matplot
#        fig, ax = plt.subplots()   
#                  
#        exec("plt.imshow("+med+"SEmatrixSR"+reg+", cmap=plt.cm.rainbow, interpolation='nearest',norm=LogNorm(vmin=1E-3, vmax=1E-1))")
#        ax.invert_yaxis()
#        ax.set_xticks([i for i in range(len(philist))], minor=False)
#        ax.set_yticks([i for i in range(len(chilist))], minor=False)
#        ax.set_xticklabels(philist)
#        ax.set_yticklabels(chilist)
#        ax.set_xlabel("$\Phi$ Mass (GeV)")
#        ax.set_ylabel("$\chi$ Mass (GeV)")
#        plt.title(med+" SR"+reg)
#        plt.colorbar() 
#        
#        for (j,i),label in np.ndenumerate(grid):
#    #        print label
#            if not label==0.: ax.text(i,j,"%.5f"%label,fontsize=7,horizontalalignment='center',verticalalignment='center')
#        
#        plt.savefig(med+"_SR"+reg+".png")


#ROOT
        c=TCanvas("c","c",1500, 950)
        c.SetLogz()
        hist2=TH2F("hist2","",len(philist),0,len(philist),len(chilist),0,len(chilist))
        for i in range(len(philist)):
            for j in range(len(chilist)):
                hist2.SetBinContent(i+1,j+1,grid[j][i])
        hist2.SetStats(0)
        hist2.SetMaximum(1e-4)
        hist2.SetMinimum(1e-5)
        hist2.GetYaxis().SetTitle("#chi Mass (GeV)")
        hist2.GetYaxis().SetTitleOffset(1.2)
        hist2.GetXaxis().SetTitle("#phi Mass (GeV)")
        plttitle="tt+DM LO "+med+" SR"+reg
        hist2.SetTitle(plttitle)
        for i in range(len(philist)):
            hist2.GetXaxis().SetBinLabel(i+1,str(philist[i]))
        for j in range(len(chilist)):
            hist2.GetYaxis().SetBinLabel(j+1,str(chilist[j]))
        
        hist2.Draw("COLZ TEXT")
        c.Print(med+"_SR"+reg+".png")

#Text File        
        outfile=open(med+"_SE_SR"+reg+".txt",'w')  
        for line in grid:
            for SEval in line:
                outfile.write(str(SEval)+" ")
            outfile.write("\n")
        outfile.close()    

