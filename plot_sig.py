import ROOT

from loadAllSamples import *
from eventlist import *
from likelihood import *
from helpers import *
from math import *

ROOT.gROOT.LoadMacro('tdrstyle.C')
ROOT.setTDRStyle()

import json

outfile = 'MC_tune.txt'

# Define working points etc
presel = 'Sum$(jet_pt>0)>1'

#samples = [WW,WZ,ZZ]#,ST_top,ST_antitop]
isData = False

# load chain to list
el_mc = eventlist( allMCSamples, presel )
el_data = eventlist( [data], presel )

with open('data/MC_tune.txt', 'r') as paraMC:
  parMC = json.load(paraMC)

with open('data/data_tune.txt', 'r') as paraData:
  parData = json.load(paraData)

el_mc.getLL(parMC)
el_data.getLL(parData)


types = ['Zmumu','top','EWK','Data']

sig_hist = {}
for t in types:
  sig_hist[t] = ROOT.TH1F('sig_'+t,t,50,0,100)

totalH = ROOT.TH1F('total','total',50,0,100)

nMC = 0.
for ev in el_mc.evlist: nMC += ev.weight

nData = 0.
for ev in el_data.evlist: nData += 1

scale = nData/nMC

for ev in el_mc.evlist + el_data.evlist:
  if ev.group is not 'Data':
    weight = ev.weight*scale
    totalH.Fill(ev.sig,weight)
  else: weight = 1
  sig_hist[ev.group].Fill(ev.sig, weight)

stack = ROOT.THStack()

sig_hist['EWK'].SetFillColor(ROOT.kRed-9)
sig_hist['top'].SetFillColor(ROOT.kYellow-9)
sig_hist['Zmumu'].SetFillColor(ROOT.kBlue-9)

stack.Add(sig_hist['top'])
stack.Add(sig_hist['EWK'])
stack.Add(sig_hist['Zmumu'])

const = sig_hist['Data'].GetBinContent(1)
const = const/ROOT.TMath.Exp(-(float(100)/50)/4.)
chi2_2 = ROOT.TF1("chi22","TMath::Exp(-x/2)*"+str(const),0,100)

can = ROOT.TCanvas('can','can',700,700)

pad1=ROOT.TPad("pad1","MyTitle",0.,0.3,1.,1.)
pad1.SetLeftMargin(0.15)
pad1.SetBottomMargin(0.02)
pad1.Draw()
pad1.cd()

pad1.SetLogy()

stack.Draw('hist')
stack.SetMinimum(0.1)
stack.SetMaximum(nData)
stack.GetXaxis().SetLabelSize(0)
stack.GetYaxis().SetTitle('Events')

chi2_2.SetLineColor(ROOT.kRed)
chi2_2.SetLineWidth(1)
chi2_2.Draw('same')

sig_hist['Data'].SetMarkerStyle(8)
sig_hist['Data'].SetMarkerSize(1)
sig_hist['Data'].SetLineWidth(0)
sig_hist['Data'].SetLineColor(ROOT.kBlack)
sig_hist['Data'].Draw('e0p same')

leg = ROOT.TLegend(0.75,0.72,0.95,0.92)
leg.SetFillColor(ROOT.kWhite)
leg.SetShadowColor(ROOT.kWhite)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)

leg.AddEntry(sig_hist['Data'],'Data')
leg.AddEntry(chi2_2,'#chi^{2} d.o.f 2','l')
leg.AddEntry(sig_hist['Zmumu'],'Z to #mu#mu','f')
leg.AddEntry(sig_hist['top'],'Top','f')
leg.AddEntry(sig_hist['EWK'],'Drell Yan','f')

leg.Draw()

can.cd()

pad2=ROOT.TPad("pad2","datavsMC",0.,0.,1.,.3)
pad2.SetLeftMargin(0.15)
pad2.SetBottomMargin(0.3)
pad2.SetTopMargin(0.04)
pad2.Draw()
pad2.cd()

ratio = ROOT.TH1F('ratio','',50,0,2)
ratio.Sumw2()
ratio = sig_hist['Data'].Clone()
ratio.Divide(totalH)
ratio.SetLineColor(ROOT.kBlack)
ratio.SetMarkerStyle(8)
ratio.SetMarkerSize(1)
ratio.SetLineWidth(0)
ratio.GetXaxis().SetTitle('')
ratio.SetMaximum(2.0)
ratio.SetMinimum(0.)
ratio.GetXaxis().SetTitle('E_{T}^{miss} Significance')
ratio.GetXaxis().SetTitleSize(0.12)
ratio.GetXaxis().SetLabelSize(0.12)
ratio.GetYaxis().SetLabelSize(0.12)
ratio.GetYaxis().SetNdivisions(108)

ratio.Draw('e0p')

one = ROOT.TF1("one","1",0,100)
one.SetLineColor(ROOT.kRed+1)
one.SetLineWidth(2)
one.Draw('same')

ratio.Draw('e0p same')

