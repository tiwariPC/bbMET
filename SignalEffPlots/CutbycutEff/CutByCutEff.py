from ROOT import *
import matplotlib.pyplot as plt
import os, numpy as np

dirpath="signal/"
mphi=[50.,100.,350.,400.,500.,1000.]
filewiseeff=[]
def getmphi(fname):
    return int(fname.split('.')[0].split('_')[4].split('-')[1])

rootlist=os.listdir(dirpath)

for fl in sorted(rootlist,key=getmphi):
    filename=os.path.join(dirpath,fl)
    f=TFile(filename,"READ")
    
    tot=f.Get('h_total_weight').Integral()
    h_cutflow=f.Get("h_cutflow_SR2_")
    
    cutwiseeff=[]
    for i in range(2,h_cutflow.GetSize()-1):
        cutwiseeff.append(h_cutflow.GetBinContent(i)/tot)
#    cutwiseeff.append(f.Get('h_met_sr2_').Integral()/tot)
        
    filewiseeff.append(cutwiseeff)
    f.Close()

print filewiseeff

for i,cts in enumerate(np.transpose(filewiseeff)):
    plt.plot(mphi,cts,'o-',label="Cut #%d"%(i+1))

plt.xscale('log')
plt.yscale('log')  
plt.xlim(40,3000)
plt.title("pseudoscalar, $M_{\chi}$ = 1 GeV")
ticks=[50.,100.,350.,500.,1000.]
plt.xticks(ticks,ticks)
plt.grid()
plt.legend()
#plt.show()
plt.savefig("out.png")
    
