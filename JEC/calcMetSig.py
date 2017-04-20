import ROOT

from metSignificance.tools.metSig import *
from metSignificance.tools.samples import samples as sample

#ROOT.gStyle.Reset()
ROOT.gROOT.LoadMacro('tdrstyle.C')
ROOT.setTDRStyle()

mS = metSig()

basedir = '/afs/cern.ch/work/d/dspitzba/CMS/MET/CMSSW_8_0_26_patch1/src/CMGTools/ObjectStudies/cfg/'
DY = sample(name='DY', treeName='METtree', subGroup='Zmumu', rootfiles=basedir+'DY_Test_2/DYJetsToLL_M50/METtree.root',    skimreport=basedir+'DY_Test_2/DYJetsToLL_M50/SkimReport.txt')
tt = sample(name='ttbar dilep', treeName='METtree', subGroup='top', rootfiles=basedir+'TTJets_DiLepton_T1/TTJets_DiLepton/METtree.root',    skimreport=basedir+'TTJets_DiLepton_T1/TTJets_DiLepton/SkimReport.txt')

#c = ROOT.TChain("METtree")
#c.Add("/afs/hephy.at/work/d/dspitzbart/MET/test/CMSSW_8_0_25/src/CMGTools/ObjectStudies/cfg/Test/DYJetsToLL_M50_7/METtree.root")
#c.Add("/afs/cern.ch/work/d/dspitzba/CMS/MET/CMSSW_8_0_26_patch1/src/CMGTools/ObjectStudies/cfg/Test/DYJetsToLL_M50_1/METtree.root")

#for ev in range(c.GetEntries()):
#    print
#    print ev
#    c.GetEntry(ev)
#    print "{:10.2f}{:10.2f}{:15.4f}{:10.0f}".format(c.met_pt,c.met_phi,c.met_sig,len(c.jet_pt))
#    mS.getMetSig(c)
#    mS.getMetSig(c, jec='up')
#    mS.getMetSig(c, jec='down')
#
#del mS
color = {
  'top':ROOT.kRed-9,
  'Zmumu':ROOT.kYellow-9,
  'EWK':ROOT.kAzure-9,
  'Data':ROOT.kBlack,
}

nBins = 25
maxSig = 50

h_metSig        = ROOT.TH1F('h_metSig','central',nBins,0,maxSig)
h_metSig_up     = ROOT.TH1F('h_metSig_up','up',nBins,0,maxSig)
h_metSig_down   = ROOT.TH1F('h_metSig_down','down',nBins,0,maxSig)

samples = [DY, tt]

h_sig = []

for s in samples:
    h_tmp = ROOT.TH1F('h_'+s.name,s.name,nBins,0,maxSig)
    h_tmp.SetFillColor(color[s.subGroup])
    h_tmp.SetLineColor(color[s.subGroup])
    print
    print "Working on %i events for sample %s"%(s.chain.GetEntries(), s.name)
    for ev in range(s.chain.GetEntries()):
    
        if ev%10000==0: print "Done with %i"%(ev)
        s.chain.GetEntry(ev)
        mS_central  = mS.getMetSig(s.chain)
        mS_up       = mS.getMetSig(s.chain, jec='sumptup')
        mS_down     = mS.getMetSig(s.chain, jec='sumptdown')
        
        weight = s.chain.xsec * s.chain.genWeight / s.sumweight
        if mS_central>maxSig: mS_central = maxSig-0.1
        if mS_up>maxSig: mS_up = maxSig-0.1
        if mS_down>maxSig: mS_down = maxSig-0.1
        h_tmp.Fill(mS_central, weight)
        h_metSig.Fill(mS_central, weight)
        h_metSig_up.Fill(mS_up, weight)
        h_metSig_down.Fill(mS_down, weight)
#        if ev > 5000: break
        
    h_sig.append(h_tmp)
    
del mS

can = ROOT.TCanvas('can','can',700,700)

stack = ROOT.THStack()
for h in h_sig:
    stack.Add(h)

h_metSig.SetLineWidth(2)
h_metSig.SetLineStyle(1)

h_metSig_up.SetLineWidth(2)
h_metSig_up.SetLineStyle(2)

h_metSig_down.SetLineWidth(2)
h_metSig_down.SetLineStyle(2)


pad1=ROOT.TPad("pad1","Main",0.,0.3,1.,1.)
pad1.SetLeftMargin(0.15)
pad1.SetBottomMargin(0.02)
pad1.Draw()
pad1.cd()
pad1.SetLogy()

stack.Draw("hist")
h_metSig.Draw("hist same")
h_metSig_up.Draw("hist same")
h_metSig_down.Draw("hist same")

nameStr = "Simulation"
lumiStr = "(MC)"
addStr=''
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

totalH = ROOT.TH1F('totalH','',nBins,0,maxSig)
for p in range(1,nBins+1):
    totalH.SetBinContent(p,1)

totalH.SetMaximum(1.5)
totalH.SetMinimum(0.5)

totalH.GetXaxis().SetTitle('E_{T}^{miss} Significance')
totalH.GetXaxis().SetTitleSize(0.12)
totalH.GetXaxis().SetLabelSize(0.12)
totalH.GetYaxis().SetLabelSize(0.12)
totalH.GetYaxis().SetNdivisions(505)
totalH.GetYaxis().SetTitle('Uncertainty')
totalH.GetYaxis().SetTitleSize(0.13)
totalH.GetYaxis().SetTitleOffset(0.45)


totalH.Draw('hist')

MCerr = ROOT.TGraphAsymmErrors(totalH)
MCerr.SetFillColor(ROOT.kGray+1)
MCerr.SetFillStyle(3001)
for p in range(0,nBins):
    if h_metSig.GetBinContent(p+1)>0:
        err = 0.5 * ( abs(h_metSig_down.GetBinContent(p+1)-h_metSig.GetBinContent(p+1))/h_metSig.GetBinContent(p+1) + abs(h_metSig_up.GetBinContent(p+1)-h_metSig.GetBinContent(p+1))/h_metSig.GetBinContent(p+1))
    else: err = 1
    if err == 0:
        print "Error was 0"
        err = 1
    print p, err
    MCerr.SetPointEXlow(p,1)
    MCerr.SetPointEXhigh(p,1)
    MCerr.SetPointEYlow(p,err)
    MCerr.SetPointEYhigh(p,err)
MCerr.Draw('2 same')

can.Print("/afs/hephy.at/user/d/dspitzbart/www/METSig/unc_est_sumpt2.png")
can.Print("/afs/hephy.at/user/d/dspitzbart/www/METSig/unc_est_sumpt2.root")
can.Print("/afs/hephy.at/user/d/dspitzbart/www/METSig/unc_est_sumpt2.pdf")


