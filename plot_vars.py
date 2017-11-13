import ROOT

from loadAllSamples_August import *
#from loadAllSamples_puppi import *
from eventlist import *
from likelihood import *
#from helpers import *
from math import *
import os
import pickle

ROOT.gStyle.Reset()
ROOT.gROOT.LoadMacro('tdrstyle.C')
ROOT.setTDRStyle()

ROOT.gROOT.ProcessLineSync('.L nvertReweight.C+')

# Define working points etc
#presel = 'Sum$(jet_pt>30&&abs(jet_eta)<2.5&&jet_passid)>=0&&(abs(jet_eta)>1.9)'
#presel = 'Sum$(jet_pt>30&&abs(jet_eta)<2.5&&jet_passid)==1&&(abs(jet_eta)<0.8)'
presel = 'Sum$(jet_pt>50&&abs(jet_eta)<2.5&&jet_passid)>=0'#&&nvertices<10'
#presel = 'Sum$(jet_pt>30&&abs(jet_eta)<2.5&&jet_passid)==0 && Sum$(ele_pt>80&&abs(ele_eta)<2.1&&ele_tight==1)==1 && Sum$(ele_pt>20&&abs(ele_eta)<2.4&&ele_loose==1)==1'#&&nvertices<10'

filename= "_njetGEq0_nvertReweight"

types = ['Zmumu','top','EWK','Data', 'QCD','gamma']

# Define Labels
lumiStr = '35.9'
nameStr = 'Preliminary'
addStr  = '2016B-H'

addUnc = False
if addUnc:
    systUnc = pickle.load(file("uncertaintyCache_noStat.pkl"))

# Define variables to plot
ToPlot = [
#{'var':'jet_sf*jet_sigmapt',   'nBins':100,'xMin':0., 'xMax':.5,'name':'#sigma_{j}'},
#{'var':'met_sig',   'nBins':50,'xMin':0, 'xMax':50,'name':'Significance'},
#{'var':'met_sig',   'nBins':50,'xMin':0, 'xMax':200,'name':'Significance'},
#{'var':'cov_xx',   'nBins':50,'xMin':0, 'xMax':4000,'name':'C_{xx}'},
#{'var':'sqrt(cov_xx)',   'nBins':30,'xMin':0, 'xMax':60,'name':'#sigma_{xx}'},
#{'var':'sqrt(abs(cov_xy))',   'nBins':30,'xMin':0, 'xMax':60,'name':'#sigma_{xy}'},
#{'var':'cov_yy',   'nBins':50,'xMin':0, 'xMax':1000,'name':'C_{yy}'},
#{'var':'met_pt',    'nBins':200,'xMin':0, 'xMax':400,'name':'E_{T}^{miss} (GeV)'},
{'var':'met_pt',    'nBins':50,'xMin':0, 'xMax':100,'name':'E_{T}^{miss} (GeV)'},
#{'var':'met_pt',    'nBins':25,'xMin':0, 'xMax':200,'name':'E_{T}^{miss} (GeV)'},
#{'var':'ele_pt[0]',    'nBins':25,'xMin':0, 'xMax':500,'name':'p_{T}^{l} (GeV)'},
#{'var':'ele_phi[0]',    'nBins':32,'xMin':-3.2, 'xMax':3.2,'name':'#phi(l)'},
#{'var':'abs(ele_eta[0])',    'nBins':24,'xMin':0, 'xMax':2.4,'name':'#eta(l)'},
#{'var':'jet_pt', 'nBins':50,'xMin':0, 'xMax':400,'name':'p_{T}(j_{i}) (GeV)'},
#{'var':'jet_pt*jet_sigmapt+jet_sf', 'nBins':50,'xMin':0, 'xMax':100,'name':'p_{T}(j_{i}) (GeV)'},
#{'var':'jet_pt*jet_sigmapt*jet_sf', 'nBins':50,'xMin':0, 'xMax':100,'name':'p_{T}(j_{i}) (GeV)'},
#{'var':'jet_phi*jet_sigmaphi', 'nBins':100,'xMin':-0.5, 'xMax':0.5,'name':'#sigma_{j}'},
#{'var':'jet_phi', 'nBins':50,'xMin':-3.2, 'xMax':3.2,'name':'#phi(j_{i})'},
{'var':'jet_pt[0]', 'nBins':50,'xMin':0, 'xMax':400,'name':'p_{T}(j_{1}) (GeV)'},
{'var':'jet_pt[1]', 'nBins':50,'xMin':0, 'xMax':400,'name':'p_{T}(j_{2}) (GeV)'},
{'var':'jet_phi[0]', 'nBins':50,'xMin':-3.2, 'xMax':3.2,'name':'#phi(j_{1})'},
{'var':'jet_phi[1]', 'nBins':50,'xMin':-3.2, 'xMax':3.2,'name':'#phi(j_{2})'},
{'var':'jet_eta[0]', 'nBins':50,'xMin':-5, 'xMax':5,'name':'#eta(j_{1})'},
{'var':'jet_eta[1]', 'nBins':50,'xMin':-5, 'xMax':5,'name':'#eta(j_{2})'},
#{'var':'muon_pt[0]','nBins':50,'xMin':0, 'xMax':200,'name':'p_{T}(#mu_{1}) (GeV)'},
#{'var':'muon_pt[1]','nBins':50,'xMin':0, 'xMax':200,'name':'p_{T}(#mu_{2}) (GeV)'},
{'var':'nvertices','nBins':100,'xMin':0, 'xMax':100,'name':'N_{vert}'},
#{'var':'nvert_true','nBins':100,'xMin':0, 'xMax':100,'name':'N_{vert}'},
#{'var':'met_sumpt','nBins':100,'xMin':0, 'xMax':3000,'name':'E_{T} unclustered (GeV)'},
#{'var':'sqrt(met_sumpt)','nBins':30,'xMin':0, 'xMax':60,'name':'sqrt(E_{T} unclustered) sqrt(GeV)'},
#{'var':'tune_sumpt',   'nBins':100,'xMin':0., 'xMax':2000.,'name':'#sigma_{j}'}
#{'var':'tune_jet_sigmapt_eta1_dataShift',   'nBins':100,'xMin':0., 'xMax':100.,'name':'#sigma_{j}'}
#{'var':'jet_sigmapt[0]*jet_pt[0]',   'nBins':100,'xMin':0., 'xMax':100.,'name':'#sigma_{j} (GeV)'}
#{'var':'tune_met_sumpt','nBins':100,'xMin':0, 'xMax':3000,'name':'E_{T} unclustered (GeV)'},
#{'var':'sqrt(2*muon_pt[0]*muon_pt[1]*(cosh(muon_eta[0]-muon_eta[1])-cos(muon_phi[0]-muon_phi[1])))','nBins':60,'xMin':60,'xMax':120,'name':'M(#mu_{1},#mu_{2} (GeV)'},
]

# Histograms
#types = ['wjets','QCD','gamma','DY','top','EWK','Data']
types = ['Zmumu','top','EWK','Data']

color = {
  'wjets':ROOT.kCyan+1,
  'QCD':ROOT.kGreen+2,
  'gamma':ROOT.kOrange,
  'top':ROOT.kRed-9,
  'Zmumu':ROOT.kYellow-9,
  'EWK':ROOT.kAzure-9,
  'Data':ROOT.kBlack,
}

#allMCSamples = allMCSamples_noHCAL
allMCSamples = allNLOSamples

for s in allMCSamples:
  s.setTargetLumi(35900)
  s.calculateWeight()

for p in ToPlot:
  var = p['var']
  xMax = p['xMax']
  xMin = p['xMin']
  nBins = p['nBins']
  varname = p['name']
  
  hists = []
  for s in allMCSamples:
    hists.append(ROOT.TH1F(s.name,s.name,nBins,xMin,xMax))
    hists[-1].SetLineColor(color[s.subGroup])
    hists[-1].SetFillColor(color[s.subGroup])
    #s.chain.Draw(var+'>>'+s.name,'('+presel+')*('+str(s.weight)+')','goff')
    #s.chain.Draw(var+'>>'+s.name,'('+presel+')*('+str(s.weight)+'*nvertReweight(nvertices))','goff')
    s.chain.Draw(var+'>>'+s.name,'('+presel+')*('+str(s.weight)+'*nvertReweight(nvertices))*mcweight','goff')

    #s.chain.Draw('(-1.20669)**2+0.611**2*met_sumpt>>'+s.name,'('+presel+')*('+str(s.weight)+'*nvertReweight(nvertices))','goff')
    #s.chain.Draw('jet_sigmapt*jet_sf*jet_pt*1.403>>'+s.name,'('+presel+')*('+str(s.weight)+'*nvertReweight(nvertices))','goff')
    #s.chain.Draw('jet_sigmapt*jet_pt*1.505>>'+s.name,'('+presel+')*('+str(s.weight)+'*nvertReweight(nvertices))','goff')


  data_hist = ROOT.TH1F('data','data',nBins,xMin,xMax)
  data.chain.Draw(var+'>>data',presel,'goff')
  #ICHEP.chain.Draw(var+'>>data',presel,'goff')
  #data.chain.Draw(var+'>>data',presel,'goff')
  #data.chain.Draw('(-0.76865)**2+(0.579)**2*met_sumpt>>data',presel,'goff')
  #data.chain.Draw('(-0.77+jet_sigmapt*jet_pt*1.481)>>data',presel,'goff')
  #data2016G.chain.Draw('-0.5066+0.5786*sqrt(met_sumpt)>>data',presel,'goff')
  
  data_mean = data_hist.GetMean()
  print 'Data mean:',data_mean
  
  totalH = ROOT.TH1F('total','total',nBins,xMin,xMax)
  stack = ROOT.THStack()
  
  for h in hists:
    totalH.Add(h)
  
  mc_mean = totalH.GetMean()
  print 'MC mean:',mc_mean
  
  print data_hist.Integral(1,nBins+1)
  print totalH.Integral(1,nBins+1)
  scaleFactor = data_hist.Integral(1,nBins+1)/totalH.Integral(1,nBins+1)
  
  for h in hists:
    h.SetBinContent(nBins, h.GetBinContent(nBins)+h.GetBinContent(nBins+1))
  totalH.SetBinContent(nBins, totalH.GetBinContent(nBins)+totalH.GetBinContent(nBins+1))
  data_hist.SetBinContent(nBins, data_hist.GetBinContent(nBins)+data_hist.GetBinContent(nBins+1))
  #stack.Add(sig_hist['top'])
  #stack.Add(sig_hist['EWK'])
  #stack.Add(sig_hist['Zmumu'])
  
  #scaleFactor = data_hist.Integral(1,nBins)/totalH.Integral(1,nBins)
  
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
  stack.SetMinimum(50)
  #stack.SetMaximum(200000)
  stack.GetXaxis().SetLabelSize(0)
  stack.GetYaxis().SetTitle('Events')
  stack.GetYaxis().SetTitleSize(0.065)
  stack.GetYaxis().SetTitleOffset(0.9)
  
  data_hist.Draw('e1p same')
  totalHUnc = totalH.Clone()
#  systUncs = {"JEC":totalH.Clone(),"JER":totalH.Clone(),"unclustered":totalH.Clone()}
  MCerr = ROOT.TGraphAsymmErrors(totalH)
  MCerr.SetFillColor(ROOT.kGray+1)
  MCerr.SetFillStyle(3244)
  for p in range(0,nBins):
    y = totalH.GetBinContent(p+1)
    yErr = totalH.GetBinError(p+1)
    if addUnc:
        yErr = y*systUnc[presel][var]['total'].GetBinError(p+1)
    totalHUnc.SetBinError(p+1,yErr)
    MCerr.SetPointEYlow(p, yErr)
    MCerr.SetPointEYhigh(p, yErr)
    MCerr.SetPointEXlow(p,xMax/(2*nBins))
    MCerr.SetPointEXhigh(p,xMax/(2*nBins))
  MCerr.Draw('2 same')
  
  data_hist.SetMarkerStyle(8)
  data_hist.SetMarkerSize(1)
  data_hist.SetLineWidth(0)
  data_hist.SetLineColor(ROOT.kBlack)
  data_hist.Draw('e0p same')
  
  #leg = ROOT.TLegend(0.75,0.64,0.95,0.92)
  leg = ROOT.TLegend(0.75,0.68,0.95,0.90)
  leg.SetFillColor(ROOT.kWhite)
  leg.SetShadowColor(ROOT.kWhite)
  leg.SetBorderSize(0)
  leg.SetTextSize(0.04)
  leg.AddEntry(data_hist,'Data')
  leg.AddEntry(hists[-1],'Z #rightarrow #mu#mu','f')
  leg.AddEntry(hists[0],'top','f')
  leg.AddEntry(hists[3],'EWK','f')
  #leg.AddEntry(hists[-1],'W+jets','f')
  #leg.AddEntry(hists[-2],'QCD multijet','f')
  #leg.AddEntry(hists[7],'#gamma+jets','f')
  #leg.AddEntry(hists[6],'Drell-Yan','f')
  #leg.AddEntry(hists[3],'EWK','f')
  #leg.AddEntry(hists[0],'top','f')
  leg.Draw()
  
  latex2 = ROOT.TLatex()
  latex2.SetNDC()
  latex2.SetTextSize(0.04)
  latex2.SetTextAlign(11)
  latex2.DrawLatex(0.15,0.96,'CMS #bf{#it{'+nameStr+'}}')
  latex2.DrawLatex(0.5,0.96,addStr)
  latex2.DrawLatex(0.88,0.96,lumiStr+'fb^{-1}')
  
  can.cd()
  
  pad2=ROOT.TPad("pad2","datavsMC",0.,0.,1.,.3)
  pad2.SetLeftMargin(0.15)
  pad2.SetBottomMargin(0.3)
  pad2.SetTopMargin(0.04)
  pad2.Draw()
  pad2.cd()
  
  ratio = ROOT.TH1F('ratio','',nBins,xMin,xMax)
  totalHUnc.Divide(totalH)
  MCerrRatio = ROOT.TGraphAsymmErrors(totalHUnc)
  
  #colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kOrange]
  #systUncs = {}
  #for i,u in enumerate(systUnc[presel][var]):
  #  print u
  #  print colors[i]
  #  print systUnc[presel][var][u].GetNbinsX()
  #  tmp = systUnc[presel][var][u].Clone()
  #  for p in range(0,nBins):
  #      tmp.SetBinContent(p+1,1)
  #      tmp.SetBinError(p+1,systUnc[presel][var][u].GetBinContent(p+1))
  #  systUncs[u] = ROOT.TGraphAsymmErrors(tmp)
  #  systUncs[u].SetFillColor(colors[i])
  #  systUncs[u].SetFillStyle(3244)
  #  del tmp
  ratio.Sumw2()
  ratio = data_hist.Clone()
  ratio.Divide(totalH)
  ratio.SetLineColor(ROOT.kBlack)
  ratio.SetMarkerStyle(8)
  ratio.SetMarkerSize(1)
  ratio.SetLineWidth(1)
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
  MCerrRatio.SetFillColor(ROOT.kGray+1)
  MCerrRatio.SetFillStyle(3244)
  for p in range(0,nBins):
    #ratio.SetBinError(p+1,sqrt(data_hist.GetBinContent(p+1))/data_hist.GetBinContent(p+1))
    MCerrRatio.SetPointEXlow(p,xMax/(2*nBins))
    MCerrRatio.SetPointEXhigh(p,xMax/(2*nBins))
  #MCerrRatio.Draw('2 same')
  ratio.Draw('e0p')
  #MCerrRatio.Draw('2 same')
  colors = [ROOT.kGray, ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kOrange]
  if addUnc:
    systUncs = {}
    for i,u in enumerate(["total","unclustered","stat","JER","JES"]):
      systUncs[u] = ROOT.TGraphAsymmErrors(systUnc[presel][var][u])
      systUncs[u].SetFillColor(colors[i])
      systUncs[u].SetFillStyle(3244)
      systUncs[u].SetLineWidth(0)
      systUncs[u].SetMarkerSize(0)
      for p in range(0,nBins):
          systUncs[u].SetPointEXlow(p,xMax/(2*nBins))
          systUncs[u].SetPointEXhigh(p,xMax/(2*nBins))
      systUncs[u].Draw('2 same')
  #MCerrRatio.Draw('2 same')

  one = ROOT.TF1("one","1",xMin,xMax)
  one.SetLineColor(ROOT.kRed+1)
  one.SetLineWidth(2)
  one.Draw('same')
  
  ratio.Draw('e0p same')
  dummy = ROOT.TH1F()
  pad1.cd()
  if addUnc:
    leg2 = ROOT.TLegend(0.55,0.64,0.75,0.90)
    leg2.SetFillColor(ROOT.kWhite)
    leg2.SetShadowColor(ROOT.kWhite)
    leg2.SetBorderSize(0)
    leg2.SetTextSize(0.04)
    leg2.AddEntry(dummy,'Uncertainties','')
    for u in ["total","unclustered","stat","JER","JES"]:
      leg2.AddEntry(systUncs[u],u)
    leg2.Draw()

  
  plotDir = "/afs/hephy.at/user/d/dspitzbart/www/METSig/2016BH_August_v1_NLO/"
  if not os.path.isdir(plotDir): os.makedirs(plotDir)
  
  can.Print(plotDir+var+'_to'+str(xMax)+filename+'.png')
  can.Print(plotDir+var+'_to'+str(xMax)+filename+'.pdf')
  can.Print(plotDir+var+'_to'+str(xMax)+filename+'.root')
  
  obDelete = hists + [stack,ratio,data_hist,totalH,pad1,pad2,leg,one]
  for o in obDelete:
    o.Delete()
  del can
