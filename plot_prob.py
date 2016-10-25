import ROOT

from loadAllSamples import *
from eventlist import *
from likelihood import *
from helpers import *
from math import *

ROOT.gStyle.Reset()
ROOT.gROOT.LoadMacro('tdrstyle.C')
ROOT.setTDRStyle()
ROOT.TGaxis.SetMaxDigits(2)

import json

## Histograms
#types = ['Zmumu','top','EWK','Data']
#nBins = 50
#maxSig = 1
#
#nMC = 0.
#for ev in el_mc.evlist: nMC += ev.weight
#nData = len(el_data.evlist)
#scale = nData/nMC
#
#sig_hist = {}
#for t in types:
#  sig_hist[t] = ROOT.TH1F('sig_'+t,t,nBins,0,maxSig)
#
#totalH = ROOT.TH1F('total','total',nBins,0,maxSig)
#
#for ev in el_mc.evlist + el_data.evlist:
#  if ev.group is not 'Data':
#    weight = ev.weight*scale
#    totalH.Fill(ROOT.TMath.Prob(ev.sig,2),weight)
#  else: weight = 1
#  sig_hist[ev.group].Fill(ROOT.TMath.Prob(ev.sig,2), weight)

stack = ROOT.THStack()

#for h in sig_hist:
#  sig_hist[h].SetBinContent(nBins, sig_hist[h].GetBinContent(nBins)+sig_hist[h].GetBinContent(nBins+1))
#totalH.SetBinContent(nBins, totalH.GetBinContent(nBins)+totalH.GetBinContent(nBins+1))

sig_hist['EWK'].SetFillColor(ROOT.kRed-10)
sig_hist['top'].SetFillColor(ROOT.kYellow-9)
sig_hist['Zmumu'].SetFillColor(ROOT.kAzure-9)

stack.Add(sig_hist['top'])
stack.Add(sig_hist['EWK'])
stack.Add(sig_hist['Zmumu'])

can = ROOT.TCanvas('can','can',700,700)

pad1=ROOT.TPad("pad1","Main",0.,0.3,1.,1.)
pad1.SetLeftMargin(0.15)
pad1.SetBottomMargin(0.02)
pad1.Draw()
pad1.cd()
#pad1.SetLogy()

stack.Draw('hist')
stack.SetMinimum(0)
stack.SetMaximum(1.2*sig_hist['Data'].GetBinContent(1))
stack.GetXaxis().SetLabelSize(0)
stack.GetYaxis().SetTitle('Events')
stack.GetYaxis().SetTitleSize(0.065)
stack.GetYaxis().SetTitleOffset(0.9)

MCerr = ROOT.TGraphAsymmErrors(totalH)
MCerr.SetFillColor(ROOT.kGray+1)
MCerr.SetFillStyle(3244)
for p in range(0,nBins):
  MCerr.SetPointEXlow(p,maxSig/(2*nBins))
  MCerr.SetPointEXhigh(p,maxSig/(2*nBins))
MCerr.Draw('2 same')

sig_hist['Data'].SetMarkerStyle(8)
sig_hist['Data'].SetMarkerSize(1)
sig_hist['Data'].SetLineWidth(0)
sig_hist['Data'].SetLineColor(ROOT.kBlack)
sig_hist['Data'].Draw('e0p same')

leg = ROOT.TLegend(0.75,0.72,0.95,0.90)
leg.SetFillColor(ROOT.kWhite)
leg.SetShadowColor(ROOT.kWhite)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)
leg.AddEntry(sig_hist['Data'],'Data')
leg.AddEntry(sig_hist['Zmumu'],'Z #rightarrow #mu#mu','f')
leg.AddEntry(sig_hist['top'],'top','f')
leg.AddEntry(sig_hist['EWK'],'EWK','f')
leg.Draw()

can.cd()

pad2=ROOT.TPad("pad2","datavsMC",0.,0.,1.,.3)
pad2.SetLeftMargin(0.15)
pad2.SetBottomMargin(0.3)
pad2.SetTopMargin(0.04)
pad2.Draw()
pad2.cd()

ratio = ROOT.TH1F('ratio','',nBins,0,maxSig)
ratio.Sumw2()
ratio = sig_hist['Data'].Clone()
ratio.Divide(totalH)
ratio.SetLineColor(ROOT.kBlack)
ratio.SetMarkerStyle(8)
ratio.SetMarkerSize(1)
ratio.SetLineWidth(0)
ratio.GetXaxis().SetTitle('')
ratio.SetMaximum(1.59)
ratio.SetMinimum(0.41)
ratio.GetXaxis().SetTitle('#chi^{2} Probability')
ratio.GetXaxis().SetTitleSize(0.12)
ratio.GetXaxis().SetLabelSize(0.12)
ratio.GetYaxis().SetLabelSize(0.12)
ratio.GetYaxis().SetNdivisions(505)
ratio.GetYaxis().SetTitle('Data/MC')
ratio.GetYaxis().SetTitleSize(0.13)
ratio.GetYaxis().SetTitleOffset(0.45)
ratio.Draw('e0p')

one = ROOT.TF1("one","1",0,maxSig)
one.SetLineColor(ROOT.kRed+1)
one.SetLineWidth(2)
one.Draw('same')

ratio.Draw('e0p same')

