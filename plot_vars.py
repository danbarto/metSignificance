import ROOT

from loadAllSamples import *
from eventlist import *
from likelihood import *
#from helpers import *
from math import *

ROOT.gStyle.Reset()
ROOT.gROOT.LoadMacro('tdrstyle.C')
ROOT.setTDRStyle()


# Define working points etc
presel = 'Sum$(jet_pt>30&&abs(jet_eta)<2.5&&jet_passid)>=0'
types = ['Zmumu','top','EWK','Data']

ToPlot = [
{'var':'met_sig','nBins':50,'xMax':100,'name':'Significance'},
{'var':'met_pt','nBins':50,'xMax':100,'name':'E_{T}^{miss} (GeV)'},
{'var':'jet_pt[0]','nBins':50,'xMax':400,'name':'p_{T}(j_{1}) (GeV)'},
{'var':'jet_pt[1]','nBins':50,'xMax':400,'name':'p_{T}(j_{2}) (GeV)'},
{'var':'muon_pt[0]','nBins':50,'xMax':200,'name':'p_{T}(#mu_{1}) (GeV)'},
{'var':'muon_pt[1]','nBins':50,'xMax':200,'name':'p_{T}(#mu_{2}) (GeV)'},
#{'var':'mu_pt[0]','nBins':50,'xMax':200},
#{'var':'mu_pt[0]','nBins':50,'xMax':200},
]

# Histograms
types = ['Zmumu','top','EWK','Data']

color = {
  'top':ROOT.kRed-9,
  'Zmumu':ROOT.kYellow-9,
  'EWK':ROOT.kAzure-9,
  'Data':ROOT.kBlack,
}

for s in allMCSamples:
  s.setTargetLumi(7100)
  s.calculateWeight()

for p in ToPlot:
  var = p['var']
  xMax = p['xMax']
  nBins = p['nBins']
  varname = p['name']
  
  hists = []
  for s in allMCSamples:
    hists.append(ROOT.TH1F(s.name,s.name,nBins,0,xMax))
    hists[-1].SetLineColor(color[s.subGroup])
    hists[-1].SetFillColor(color[s.subGroup])
    s.chain.Draw(var+'>>'+s.name,'('+presel+')*'+str(s.weight),'goff')
  
  data_hist = ROOT.TH1F('data','data',nBins,0,xMax)
  #data.chain.Draw(var+'>>data',presel,'goff')
  #ICHEP.chain.Draw(var+'>>data',presel,'goff')
  data2016G.chain.Draw(var+'>>data',presel,'goff')
  
  
  totalH = ROOT.TH1F('total','total',nBins,0,xMax)
  
  for h in hists:
    h.SetBinContent(nBins, h.GetBinContent(nBins)+h.GetBinContent(nBins+1))
  totalH.SetBinContent(nBins, totalH.GetBinContent(nBins)+totalH.GetBinContent(nBins+1))
  data_hist.SetBinContent(nBins, data_hist.GetBinContent(nBins)+data_hist.GetBinContent(nBins+1))
  
  stack = ROOT.THStack()
  
  for h in hists:
    totalH.Add(h)
  #stack.Add(sig_hist['top'])
  #stack.Add(sig_hist['EWK'])
  #stack.Add(sig_hist['Zmumu'])
  
  scaleFactor = data_hist.Integral()/totalH.Integral()
  
  print 'Data MC SF, should be close to 1'
  print scaleFactor
  
  totalH.Scale(scaleFactor)
  for h in hists:
    h.Scale(scaleFactor)
    stack.Add(h)
    
  
  can = ROOT.TCanvas('can','can',700,700)
  
  pad1=ROOT.TPad("pad1","Main",0.,0.3,1.,1.)
  pad1.SetLeftMargin(0.15)
  pad1.SetBottomMargin(0.02)
  pad1.Draw()
  pad1.cd()
  pad1.SetLogy()
  
  stack.Draw('hist')
  stack.SetMinimum(1)
  #stack.SetMaximum(200000)
  stack.GetXaxis().SetLabelSize(0)
  stack.GetYaxis().SetTitle('Events')
  stack.GetYaxis().SetTitleSize(0.065)
  stack.GetYaxis().SetTitleOffset(0.9)
  
  data_hist.Draw('e1p same')
  
  MCerr = ROOT.TGraphAsymmErrors(totalH)
  MCerr.SetFillColor(ROOT.kGray+1)
  MCerr.SetFillStyle(3244)
  for p in range(0,nBins):
    MCerr.SetPointEXlow(p,xMax/(2*nBins))
    MCerr.SetPointEXhigh(p,xMax/(2*nBins))
  MCerr.Draw('2 same')
  
  data_hist.SetMarkerStyle(8)
  data_hist.SetMarkerSize(1)
  data_hist.SetLineWidth(0)
  data_hist.SetLineColor(ROOT.kBlack)
  data_hist.Draw('e0p same')
  
  leg = ROOT.TLegend(0.75,0.68,0.95,0.90)
  leg.SetFillColor(ROOT.kWhite)
  leg.SetShadowColor(ROOT.kWhite)
  leg.SetBorderSize(0)
  leg.SetTextSize(0.04)
  leg.AddEntry(data_hist,'Data')
  leg.AddEntry(hists[-1],'Z #rightarrow #mu#mu','f')
  leg.AddEntry(hists[0],'top','f')
  leg.AddEntry(hists[3],'EWK','f')
  leg.Draw()
  
  can.cd()
  
  pad2=ROOT.TPad("pad2","datavsMC",0.,0.,1.,.3)
  pad2.SetLeftMargin(0.15)
  pad2.SetBottomMargin(0.3)
  pad2.SetTopMargin(0.04)
  pad2.Draw()
  pad2.cd()
  
  ratio = ROOT.TH1F('ratio','',nBins,0,xMax)
  ratio.Sumw2()
  ratio = data_hist.Clone()
  ratio.Divide(totalH)
  ratio.SetLineColor(ROOT.kBlack)
  ratio.SetMarkerStyle(8)
  ratio.SetMarkerSize(1)
  ratio.SetLineWidth(0)
  ratio.GetXaxis().SetTitle('')
  ratio.SetMaximum(2.01)
  ratio.SetMinimum(0.)
  ratio.GetXaxis().SetTitle(varname)
  ratio.GetXaxis().SetTitleSize(0.12)
  ratio.GetXaxis().SetLabelSize(0.12)
  ratio.GetYaxis().SetLabelSize(0.12)
  ratio.GetYaxis().SetNdivisions(505)
  ratio.GetYaxis().SetTitle('Data/MC')
  ratio.GetYaxis().SetTitleSize(0.13)
  ratio.GetYaxis().SetTitleOffset(0.45)
  ratio.Draw('e0p')
  
  one = ROOT.TF1("one","1",0,xMax)
  one.SetLineColor(ROOT.kRed+1)
  one.SetLineWidth(2)
  one.Draw('same')
  
  ratio.Draw('e0p same')
  
  can.Print('/afs/hephy.at/user/d/dspitzbart/www/METSig/2016G_Nov16/'+var+'_to'+str(xMax)+'_njet0_30.png')
  can.Print('/afs/hephy.at/user/d/dspitzbart/www/METSig/2016G_Nov16/'+var+'_to'+str(xMax)+'_njet0_30.pdf')
  can.Print('/afs/hephy.at/user/d/dspitzbart/www/METSig/2016G_Nov16/'+var+'_to'+str(xMax)+'_njet0_30.root')
  
  obDelete = hists + [stack,ratio,data_hist,totalH,pad1,pad2,leg,one]
  for o in obDelete:
    o.Delete()
  del can