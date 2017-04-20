import ROOT

from loadAllSamples_aprilRerun import *
from eventlist import *
from likelihood import *
#from helpers import *
from math import *
import os

ROOT.gStyle.Reset()
ROOT.gROOT.LoadMacro('tdrstyle.C')
ROOT.setTDRStyle()

#ROOT.gROOT.ProcessLineSync('.x nvertReweight.C+')

# Define working points etc
presel = 'Sum$(jet_pt>30&&abs(jet_eta)<2.5&&jet_passid)>=0'# && nvertices>=20'
types = ['Zmumu','top','EWK','Data']

# Define Labels
lumiStr = '(13TeV)'
nameStr = 'Simulation'
addStr  = 'Moriond17'

# Histograms
types = ['Zmumu','top','EWK','Data']

color = {
  'top':ROOT.kRed-9,
  'Zmumu':ROOT.kYellow-9,
  'EWK':ROOT.kAzure-9,
  'Data':ROOT.kBlack,
}

#allMCSamples = allMCSamples_noHCAL

for s in allMCSamples:
  s.setTargetLumi(35900)
  s.calculateWeight()

colors = [ROOT.kBlue+2, ROOT.kBlue-4, ROOT.kBlue-7, ROOT.kBlue-9, ROOT.kCyan-9, ROOT.kCyan-6, ROOT.kCyan-2,ROOT.kGreen+3,ROOT.kGreen-2,ROOT.kGreen-6,ROOT.kGreen-7, ROOT.kOrange-4, ROOT.kOrange+1, ROOT.kOrange+8, ROOT.kRed, ROOT.kRed+1]

can = ROOT.TCanvas('can','can',700,700)

points = []
points.append(ROOT.TGraph())
points[0].SetPoint(1,150.,30.)
#points[0].SetPoint(1,2.,40.)
points[0].GetXaxis().SetTitle('E_{T}^{miss} (GeV)')
points[0].GetYaxis().SetTitle('#sigma(E_{T}^{miss})')
#points[0].GetXaxis().SetTitle('#delta')
#points[0].GetYaxis().SetTitle('n_{vert}')
points[0].SetMarkerSize(0)
points[0].Draw('ap')


for s,sample in enumerate(allMCSamples):

    c = sample.chain
    
    #can = ROOT.TCanvas('can','can',700,700)
    c.Draw('>>eList',presel)
    elist = ROOT.gDirectory.Get("eList")
    number_events = elist.GetN()
    totWeight = 0.
    low = 0.
    reg1 = reg2 = reg3 = reg4 = 0
    #points = []
    #points.append(ROOT.TGraph())
    #points[0].SetPoint(1,150.,30.)
    ##points[0].SetPoint(1,2.,40.)
    #points[0].GetXaxis().SetTitle('E_{T}^{miss} (GeV)')
    #points[0].GetYaxis().SetTitle('#sigma(E_{T}^{miss})')
    ##points[0].GetXaxis().SetTitle('#delta')
    ##points[0].GetYaxis().SetTitle('n_{vert}')
    #points[0].SetMarkerSize(0)
    #points[0].Draw('ap')
    
    number_events = min(10000, number_events)
    print sample.name
    for i in range(number_events):
      c.GetEntry(elist.GetEntry(i))
      if i%1000==0: print "done with", i
      met_pt    = c.met_pt
      met_sig   = c.met_sig
      met_sumpt = c.met_sumpt
      nvertices = c.nvertices
      uncEnUnc  = (c.met_PFT1UnclusteredEnUp_sig - c.met_PFT1UnclusteredEnDown_sig)/c.met_sig
    
      points.append(ROOT.TGraph())
    
      #points[-1].SetPoint(0,uncEnUnc,nvertices)
      points[-1].SetPoint(0,met_pt,met_sig)
      points[-1].SetMarkerStyle(20)
      pointSize=.6
      points[-1].SetMarkerSize(pointSize)
      points[-1].SetMarkerColor(colors[s*2])

      #for a in range(0,16):
      #  if (met_sumpt/200.)<a:
      #    #if weight<0.: points[-1].SetMarkerColor(ROOT.kMagenta)
      #    points[-1].SetMarkerColor(colors[a])
      #    break
      points[-1].Draw('p same')
    
    
    
    #dot = []
    #
    #tex =  ROOT.TLatex()
    ##tex.SetNDC()
    #tex.SetTextSize(0.025)
    #tex.SetTextAlign(11)
    #
    #tex.DrawLatex(18.5,31,"#bf{uncl.Energy (GeV)}")
    ##tex.DrawLatex(0.6,0.87,"#Delta#Phi(W,l) > 1.")
    #
    #for t in range(16):
    #  dot.append(ROOT.TGraph())
    #  #dot[-1].SetPoint(0,450,550-40*t)
    #  dot[-1].SetPoint(0,15,30-0.9*t)
    #  #dot[-1].SetPoint(0,435,458-26.3*t)
    #  dot[-1].SetMarkerStyle(20)
    #  dot[-1].SetMarkerSize(2)#3
    #  dot[-1].SetMarkerColor(colors[t])
    #  dot[-1].Draw('p same')
    #  tex.DrawLatex(18.5,30-0.9*t-0.45,'#bf{['+str(t*200)+', '+str((t+1)*200)+')}')
    
    dot = []
    
    tex =  ROOT.TLatex()
    #tex.SetNDC()
    tex.SetTextSize(0.025)
    tex.SetTextAlign(11)
    
    tex.DrawLatex(18.5,31,"#bf{Process}")
    #tex.DrawLatex(0.6,0.87,"#Delta#Phi(W,l) > 1.")
    
    for t in range(len(allMCSamples)):
      dot.append(ROOT.TGraph())
      #dot[-1].SetPoint(0,450,550-40*t)
      dot[-1].SetPoint(0,15,30-0.9*t)
      #dot[-1].SetPoint(0,435,458-26.3*t)
      dot[-1].SetMarkerStyle(20)
      dot[-1].SetMarkerSize(2)#3
      dot[-1].SetMarkerColor(colors[t*2])
      dot[-1].Draw('p same')
      tex.DrawLatex(18.5,30-0.9*t-0.45,'#bf{'+allMCSamples[t].name+'}')



plotDir = "/afs/hephy.at/user/d/dspitzbart/www/METSig/2016BH_Moriond17_rerunApril/samples/"
if not os.path.isdir(plotDir): os.makedirs(plotDir)

can.Print(plotDir+'samples_njetGEq0.png')
can.Print(plotDir+'samples_njetGEq0.pdf')
can.Print(plotDir+'samples_njetGEq0.root')

#del can, points, dot
del can, points

#test = ROOT.TH2F("test","test",100,0,3,100,0,100)
#DY_M17.chain.Draw("nvert_true:((met_PFT1UnclusteredEnUp_sig - met_PFT1UnclusteredEnDown_sig)/met_sig)>>test",'Sum$(jet_pt>30&&abs(jet_eta)<2.5&&jet_passid)==0')
#test.Draw("colz")


