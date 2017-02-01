import ROOT

from loadAllSamples import *
from eventlist import *
from likelihood import *
#from helpers import *
from math import *

ROOT.gStyle.Reset()
ROOT.gROOT.LoadMacro('tdrstyle.C')
ROOT.setTDRStyle()
#ROOT.gStyle.SetOptStat(111111)

import json

# Define working points etc
presel = 'Sum$(jet_pt>30&&abs(jet_eta)<2.5&&jet_passid)>=0'

samplesMC   = allMCSamples
samplesData = [data]
#samplesData = [data]
isData = False


# load chain to list
el_data = eventlist( samplesData, presel, isData=True, tiny=False )
el_MC   = eventlist( samplesMC, presel, isData=False, tiny=False  )

del samplesData, samplesMC

el_data.getPileUpDist()
el_MC.getPileUpDist()

el_MC.doPileUpReweight(el_data.PUhist)

el_MC.doSmearing()
#el_data.doJetSmearing(False)

with open('data/MC_2016BCDEFG_0jet30_smear.txt', 'r') as paraMC:
  parMC = json.load(paraMC)

with open('data/data_2016BCDEFG_0jet30_smear.txt', 'r') as paraData:
  parData = json.load(paraData)

el_MC.getLL(parMC)
el_data.getLL(parData)

# Histograms
types = ['Zmumu','top','EWK','Data']
nBins = 50
maxSig = 100

nMC = 0.
for ev in el_MC.evlist: nMC += ev.weight
nData = len(el_data.evlist)
scale = nData/nMC

sig_hist = {}
for t in types:
  sig_hist[t] = ROOT.TH1F('sig_'+t,t,nBins,0,maxSig)

totalH = ROOT.TH1F('total','total',nBins,0,maxSig)

for ev in el_MC.evlist + el_data.evlist:
  if ev.group is not 'Data':
    weight = ev.weight*scale
    totalH.Fill(ev.sig,weight)
  else: weight = 1
  sig_hist[ev.group].Fill(ev.sig, weight)

stack = ROOT.THStack()

for h in sig_hist:
  sig_hist[h].SetBinContent(nBins, sig_hist[h].GetBinContent(nBins)+sig_hist[h].GetBinContent(nBins+1))
totalH.SetBinContent(nBins, totalH.GetBinContent(nBins)+totalH.GetBinContent(nBins+1))

sig_hist['EWK'].SetFillColor(ROOT.kAzure-9)
sig_hist['EWK'].SetLineColor(ROOT.kAzure-9)

sig_hist['top'].SetFillColor(ROOT.kRed-9)
sig_hist['top'].SetLineColor(ROOT.kRed-9)

sig_hist['Zmumu'].SetFillColor(ROOT.kYellow-9)
sig_hist['Zmumu'].SetLineColor(ROOT.kYellow-9)


stack.Add(sig_hist['top'])
stack.Add(sig_hist['EWK'])
stack.Add(sig_hist['Zmumu'])

const = sig_hist['Data'].GetBinContent(1)
const = const/ROOT.TMath.Exp(-(float(maxSig)/nBins)/4.)
chi2_2 = ROOT.TF1("chi22","TMath::Exp(-x/2)*"+str(const),0,maxSig)

can = ROOT.TCanvas('can','can',700,700)

pad1=ROOT.TPad("pad1","Main",0.,0.3,1.,1.)
pad1.SetLeftMargin(0.15)
pad1.SetBottomMargin(0.02)
pad1.Draw()
pad1.cd()
pad1.SetLogy()

stack.Draw('hist')
stack.SetMinimum(1)
#stack.SetMaximum(nData)
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

chi2_2.SetLineColor(ROOT.kRed)
chi2_2.SetLineWidth(1)
chi2_2.Draw('same')

sig_hist['Data'].SetMarkerStyle(8)
sig_hist['Data'].SetMarkerSize(1)
sig_hist['Data'].SetLineWidth(0)
sig_hist['Data'].SetLineColor(ROOT.kBlack)
sig_hist['Data'].Draw('e0p same')

leg = ROOT.TLegend(0.75,0.68,0.95,0.90)
leg.SetFillColor(ROOT.kWhite)
leg.SetShadowColor(ROOT.kWhite)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)
leg.AddEntry(sig_hist['Data'],'Data')
leg.AddEntry(chi2_2,'#chi^{2} d.o.f 2','l')
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
ratio.SetMaximum(2.01)
ratio.SetMinimum(0.)
ratio.GetXaxis().SetTitle('E_{T}^{miss} Significance')
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

can.Print('/afs/hephy.at/user/d/dspitzbart/www/METSig/2016BCDEFG_Nov16/significanceTune_njet0_30_smear_highStat_red3_smearedInPlot.png')
can.Print('/afs/hephy.at/user/d/dspitzbart/www/METSig/2016BCDEFG_Nov16/significanceTune_njet0_30_smear_highStat_red3_smearedInPlot.pdf')
can.Print('/afs/hephy.at/user/d/dspitzbart/www/METSig/2016BCDEFG_Nov16/significanceTune_njet0_30_smear_highStat_red3_smearedInPlot.root')

