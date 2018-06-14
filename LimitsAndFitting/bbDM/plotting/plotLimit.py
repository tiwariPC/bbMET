import os
import sys 
from ROOT import * 

gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
gROOT.SetBatch(1)

c = TCanvas("c","c",1500, 950)

if len(sys.argv)>=3:
    logy=int(sys.argv[2])
else:
    logy=1
    
c.SetLogy(logy)
c.SetLogx(1)
c.SetGridx(1)
c.SetGridy(1)
rootfilepath=""
if len(sys.argv)<2: 
    print ("tell me which model you need: s, ps, 2h, tanb50, tanb100, sinp50? ")
    sys.exit()
if sys.argv[1]=="s":
    rootfilename = "scalar.root"
    
if sys.argv[1]=="ps":
    rootfilename = "pseudo.root"
    
if sys.argv[1].lower()=="2h":
    rootfilename = "2HDM.root"
    
if sys.argv[1].lower()=="tanb50":
    rootfilename = "VariabletanB50.root"
    
if sys.argv[1].lower()=="tanb100":
    rootfilename = "VariabletanB100.root"

if sys.argv[1].lower()=="sinp50":
    rootfilename = "Variablesinp50.root"


f = TFile(rootfilepath + rootfilename,"read")


exp2s =  f.Get("exp2")
exp2s.SetMarkerStyle(20)
exp2s.SetMarkerSize(1.1)
exp2s.SetLineWidth(2)
exp2s.SetFillColor(kYellow);
exp2s.SetLineColor(kYellow)
xmin=1e-6
if sys.argv[1].lower()=="2h":
    exp2s.GetXaxis().SetTitle("M_{H4} [GeV]");
    exp2s.GetXaxis().SetRangeUser(40,700)
    exp2s.GetYaxis().SetRangeUser(.05,10)
    leg = TLegend(.30, .65, .55, .890);
    xmin=39
elif sys.argv[1].lower().startswith("tanb"):
    exp2s.GetXaxis().SetTitle("tan#beta");
    leg = TLegend(.50, .65, .75, .890);
    exp2s.GetYaxis().SetRangeUser(.04,500)
elif sys.argv[1].lower().startswith("sinp"):
    exp2s.GetXaxis().SetTitle("sin#theta");
    exp2s.GetXaxis().SetRangeUser(.09,1.5)
    leg = TLegend(.50, .65, .75, .890);
    exp2s.GetYaxis().SetRangeUser(.05,10)
    xmin=.089
else:
    exp2s.GetXaxis().SetTitle("M_{#phi} [GeV]");
    exp2s.GetXaxis().SetRangeUser(40,700)
    exp2s.GetYaxis().SetRangeUser(.05,10)
    leg = TLegend(.30, .65, .55, .890);

exp2s.GetXaxis().SetTitleOffset(1.4)
exp2s.GetYaxis().SetTitle("#sigma/#sigma_{th}"); 
exp2s.GetYaxis().SetTitleOffset(1.2)
exp2s.GetYaxis().SetNdivisions(505);
exp2s.GetXaxis().SetNdivisions(505);
exp2s.GetYaxis().SetMoreLogLabels()
exp2s.GetXaxis().SetMoreLogLabels()
exp2s.Draw("A 3")

exp1s =  f.Get("exp1")
exp1s.SetMarkerStyle(20)
exp1s.SetMarkerSize(1.1)
exp1s.SetLineWidth(2)
exp1s.SetFillColor(kGreen);
exp1s.SetLineColor(kGreen)
exp1s.Draw("3 same")

exp =  f.Get("expmed")
exp.SetMarkerStyle(1)
exp.SetMarkerSize(1.1)
exp.SetLineWidth(3)
exp.Draw("L same")

obs =  f.Get("obs")
obs.SetMarkerStyle(20)
#obs.SetMarkerColor(4)
obs.SetMarkerSize(1.1)
#obs.SetLineColor(2)
obs.SetLineWidth(3)
obs.SetLineStyle(9)
obs.Draw("P same")

#hr = c.DrawFrame(fr_left, fr_down, fr_right, fr_up, "");
#hr.SetXTitle("M_{Z'} [GeV]");
#hr.SetYTitle("95% CLs on #sigma(Z`#rightarrow#chi#bar{#chi}H)#timesBR(H#rightarrowb#bar{b})[pb]");
#hr.SetMinimum(0.001);
#hr.SetMaximum(1000);
#hr.Draw(same)



leg.SetFillColor(0);
leg.SetShadowColor(0);
leg.SetTextFont(42);
leg.SetTextSize(0.03);
leg.AddEntry(obs, "CL_{S} Observed", "LP");
leg.AddEntry(exp1s, "CL_{S}  Expected #pm 1#sigma", "LF");
leg.AddEntry(exp2s, " CL_{S}  Expected #pm 2#sigma", "LF");
leg.AddEntry(exp, " CL_{S}  Expected ", "L");

leg.Draw("same")

line = TLine(max(xmin,exp2s.GetXaxis().GetXmin()),1,exp2s.GetXaxis().GetXmax(),1)
line.SetLineColor(kRed)
line.SetLineWidth(4)
line.Draw()
  

latex =  TLatex();
latex.SetNDC();
latex.SetTextSize(0.04);
latex.SetTextAlign(31);
latex.SetTextAlign(11);
model_ = rootfilename.replace(".root","")
if model_ == "pseudo":
    latex.DrawLatex(0.18, 0.93, "bb+DM pseudoscalar");

if model_ == "scalar":
    latex.DrawLatex(0.18, 0.93, "bb+DM scalar");

if model_ == "2HDM":
    latex.DrawLatex(0.18, 0.93, "bb+DM 2HDM+a: tan#beta = 35, sin#theta = .707");
    
if model_.startswith("VariabletanB"):
    latex.DrawLatex(0.18, 0.93, "bb+DM 2HDM+a: M_{H4} = "+model_[12:]+" GeV, sin#theta = .707");
    
if model_.startswith("Variablesinp"):
    latex.DrawLatex(0.18, 0.93, "bb+DM 2HDM+a: M_{H4} = "+model_[12:]+" GeV, tan#beta = 35");


name = "limit_"+model_
c.SaveAs(name+".png")
c.SaveAs(name+".pdf")
