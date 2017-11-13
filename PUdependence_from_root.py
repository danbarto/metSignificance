import ROOT
import os, pickle

ROOT.gStyle.Reset()
ROOT.gROOT.LoadMacro('tdrstyle.C')
ROOT.setTDRStyle()

#path = "~/www//METSig/2016BH_Moriond17_PUdependence_WtoENu/"#met_sig_to100_1nvert2_noNvert.root"
#path = "~/www//METSig/2016BH_Moriond17_singleEle/"
path = "~/www/METSig/2016BH_Moriond17_PUdependence/"

MC = ROOT.TH1F("MC","",48,1,49)
data = ROOT.TH1F("data","",48,1,49)

PUrange = range(1,50)
#PUrange = range(5,11,5)
res = []
for a in PUrange:
    if a < max(PUrange):
        loCut = PUrange[a-min(PUrange)]
        hiCut = PUrange[a-min(PUrange)+1]
    else: break
    key = (loCut,hiCut)
    res.append(key)


for i,a in enumerate(res):
    print
    print a,i+1
    print "Reading from file"
    print path+"met_sig_to100_%invert%i_noNvert.root"%(a[0],a[1])
    f = ROOT.TFile(path+"met_sig_to100_%invert%i_noNvert.root"%(a[0],a[1]))
    #f = ROOT.TFile(path+"Jack_njetEq0_PU_nVert0to5.root")
    can = f.Get("can")
    pad = can.GetPrimitive("pad1")
    dataH = pad.GetPrimitive("data")
    #dataH = pad.GetPrimitive("sig_Data")
    data_mean = dataH.GetMean()
    hists = pad.GetPrimitive("").GetHists()
    MCH = dataH.Clone()
    MCH.Reset()
    for hist in hists:
        MCH.Add(hist)
    mc_mean = MCH.GetMean()
    MC.SetBinContent(i+1, mc_mean)
    MC.SetBinError(i+1, MCH.GetMeanError())
    data.SetBinContent(i+1, data_mean)
    data.SetBinError(i+1, data.GetMeanError())
    f.Close()
    print mc_mean, data_mean
    del f,can,pad,dataH,MCH,hists



can = ROOT.TCanvas("can","",700,700)

MC.SetLineColor(ROOT.kOrange)
MC.SetLineWidth(2)
MC.SetMarkerSize(0)
MC.GetXaxis().SetTitle("Number of vertices")
MC.GetYaxis().SetTitle("<\\pazoscr{S}>")
MC.GetYaxis().SetRangeUser(0.,4)

MC.Draw("hist")
data.Draw("p same")


leg = ROOT.TLegend(0.65,0.80,0.95,0.90)
leg.SetFillColor(ROOT.kWhite)
leg.SetShadowColor(ROOT.kWhite)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)
leg.AddEntry(MC,'simulation')
leg.AddEntry(data,'data')
leg.Draw()

l = ROOT.TMathText()
l.SetTextSize(0.04);
l.SetTextFont(132);
l.SetTextAlign(12);
l.DrawMathText(0.5, 0.5, "\\pazocal{L} \\mathcal{L} \\mathscr{L} \\mathscr{S}")

latex2 = ROOT.TLatex()
latex2.SetNDC()
latex2.SetTextSize(0.04)
latex2.SetTextAlign(11)
latex2.DrawLatex(0.15,0.96,'CMS #bf{#it{Preliminary}} \\pazocal{L} \\mathcal{L} \\mathscr{L}')
latex2.DrawLatex(0.88,0.96,'35.9fb^{-1}')

can.Print("/afs/hephy.at/user/d/dspitzbart/www/METSig/2016BH_Moriond17_PUdependence/PU_paper.png")
can.Print("/afs/hephy.at/user/d/dspitzbart/www/METSig/2016BH_Moriond17_PUdependence/PU_paper.pdf")
#can.Print("/afs/hephy.at/user/d/dspitzbart/www/METSig/2016BH_Moriond17_PUdependence_WtoENu/PU_noNvert_WtoENu.png")
