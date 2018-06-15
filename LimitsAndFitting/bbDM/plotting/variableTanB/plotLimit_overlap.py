import os
import sys 
from ROOT import * 

gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
gROOT.SetBatch(1)

collist=[kRed,kGreen,kBlue,kYellow,kCyan,kMagenta,kOrange,kSpring,kTeal,kAzure,kViolet,kPink]

c = TCanvas("c","c",1500, 950)

if len(sys.argv)>=3:
    logy=int(sys.argv[2])
else:
    logy=1
    
c.SetLogy(logy)

c.SetGridx(1)
c.SetGridy(1)
if len(sys.argv)<2: 
    print ("Usage: python plotLimit_overlap.py tanb/sinp")
    sys.exit()
    
if sys.argv[1].lower()=="tanb":
    rootfiledir = "tanB"
    
if sys.argv[1].lower()=="sinp":
    rootfiledir = "sinp"
    
    
def getM(flname):
    return int(flname.split('.')[0].split('-')[1])

filelist=os.listdir(rootfiledir)

exp1strans=[]
exp1s=[]
exp=[]
obs=[]

leg = TLegend(.70, .55, .88, .88);
leg.SetFillColor(0);
leg.SetShadowColor(0);
leg.SetTextFont(42);
leg.SetTextSize(0.03);

for ifl,fl in enumerate(sorted(filelist,key=getM)):
    f = TFile(rootfiledir+"/"+fl,"read")
    print fl

    exp1strans.append(f.Get("exp1"))
    exp1strans[ifl].SetMarkerStyle(20)
    exp1strans[ifl].SetMarkerSize(1.1)
    exp1strans[ifl].SetLineWidth(2)
    exp1strans[ifl].SetFillColorAlpha(collist[ifl], 0.25);
    exp1strans[ifl].SetLineColor(collist[ifl]+1)
    xmin=1e-6
    
    if sys.argv[1].lower()=="tanb":
        exp1strans[ifl].GetXaxis().SetTitle("tan#beta");
        exp1strans[ifl].GetYaxis().SetRangeUser(.04,8000)
        c.SetLogx(1)
        
    elif sys.argv[1].lower()=="sinp":
        exp1strans[ifl].GetXaxis().SetTitle("sin#theta");
#        exp1strans[ifl].GetXaxis().SetRangeUser(0,1.1)
#        leg = TLegend(.50, .65, .75, .890);
        exp1strans[ifl].GetYaxis().SetRangeUser(40,1e5)
#        xmin=.089

    exp1strans[ifl].GetXaxis().SetTitleOffset(1.4)
    exp1strans[ifl].GetYaxis().SetTitle("#sigma/#sigma_{th}"); 
    exp1strans[ifl].GetYaxis().SetTitleOffset(1.2)
    exp1strans[ifl].GetYaxis().SetNdivisions(505);
    exp1strans[ifl].GetXaxis().SetNdivisions(505);
    exp1strans[ifl].GetYaxis().SetMoreLogLabels()
    exp1strans[ifl].GetXaxis().SetMoreLogLabels()
    if ifl==0:
        exp1strans[ifl].Draw("A 3")
    else:
        exp1strans[ifl].Draw("3 same")

#    exp1s.append(f.Get("exp1"))
#    exp1s[ifl].SetMarkerStyle(20)
#    exp1s[ifl].SetMarkerSize(1.1)
#    exp1s[ifl].SetLineWidth(2)
#    exp1s[ifl].SetFillColorAlpha(kGreen, 0.35);
#    exp1s[ifl].SetLineColor(kGreen)
#    exp1s[ifl].Draw("3 same")

    exp.append(f.Get("expmed"))
    exp[ifl].SetMarkerStyle(1)
    exp[ifl].SetMarkerSize(1.1)
    exp[ifl].SetLineWidth(3)
    exp[ifl].SetLineColor(collist[ifl]+1)
    exp[ifl].Draw("L same")

    obs.append(f.Get("obs"))
    obs[ifl].SetMarkerStyle(20)
    #obs.SetMarkerColor(4)
    obs[ifl].SetMarkerSize(1.1)
    #obs.SetLineColor(2)
    obs[ifl].SetLineWidth(3)
    obs[ifl].SetLineStyle(9)
    obs[ifl].SetLineColor(collist[ifl]+1)
    obs[ifl].Draw("P same")



#leg.AddEntry(obs, "CL_{S} Observed", "LP");
#leg.AddEntry(exp1s, "CL_{S}  Expected #pm 1#sigma", "LF");
#leg.AddEntry(exp1strans, " CL_{S}  Expected #pm 2#sigma", "LF");
    leg.AddEntry(exp1strans[ifl], "M_{H4} = "+str(getM(fl))+" GeV", "LF");

leg.Draw("same")

line = TLine(max(xmin,exp1strans[0].GetXaxis().GetXmin()),1,exp1strans[0].GetXaxis().GetXmax(),1)
line.SetLineColor(kRed)
line.SetLineWidth(4)
line.Draw()
  

latex =  TLatex();
latex.SetNDC();
latex.SetTextSize(0.04);
latex.SetTextAlign(31);
latex.SetTextAlign(11);

if sys.argv[1].lower()=="tanb":
    latex.DrawLatex(0.18, 0.93, "bb+DM 2HDM+a: sin#theta = .707");
    
if sys.argv[1].lower()=="sinp":
    latex.DrawLatex(0.18, 0.93, "bb+DM 2HDM+a: tan#beta = 1");


name = "limit_"+sys.argv[1].lower()
c.SaveAs(name+".png")
c.SaveAs(name+".pdf")
