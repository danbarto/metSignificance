import ROOT
import os, pickle

ROOT.gStyle.Reset()
ROOT.gROOT.LoadMacro('tdrstyle.C')
ROOT.setTDRStyle()

res = pickle.load(file("PU_results.pkl"))

MC = ROOT.TH1F("MC","",48,1,49)
data = ROOT.TH1F("data","",48,1,49)

for i,a in enumerate(sorted(res)):
    MC.SetBinContent(i+1, res[a]['MC'])
    data.SetBinContent(i+1, res[a]['data'])



can = ROOT.TCanvas("can","",700,700)

MC.SetLineColor(ROOT.kOrange)
MC.SetLineWidth(2)
MC.SetMarkerSize(0)
MC.GetXaxis().SetTitle("Number of vertices")
MC.GetYaxis().SetTitle("<significance>")
MC.GetYaxis().SetRangeUser(0.,3.)

MC.Draw("hist")
data.Draw("p same")


leg = ROOT.TLegend(0.25,0.35,0.65,0.45)
leg.SetFillColor(ROOT.kWhite)
leg.SetShadowColor(ROOT.kWhite)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)
leg.AddEntry(MC,'simulation')
leg.AddEntry(data,'data')
leg.Draw()

latex2 = ROOT.TLatex()
latex2.SetNDC()
latex2.SetTextSize(0.04)
latex2.SetTextAlign(11)
latex2.DrawLatex(0.15,0.96,'CMS #bf{#it{Preliminary}}')
latex2.DrawLatex(0.88,0.96,'35.9fb^{-1}')


can.Print("/afs/hephy.at/user/d/dspitzbart/www/METSig/2016BH_Moriond17_PUdependence/PU_paper.png")
