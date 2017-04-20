import ROOT
from loadAllSamples import *

allSamples = [DY,TTJets_DiLepton,TTJets_SingleLeptonFromTbar,TTJets_SingleLeptonFromT,T_tWch,TBar_tWch,TToLep_tch,TBarToLep_tch,TToLep_sch,WWTo2L2Nu,WZTo2L2Q,WZTo3LNu,ZZTo2L2Nu,ZZTo2L2Q,ZZTo4L,WWW,WWZ,WZZ,ZZZ]

#sigSamples = [TTJets_DiLepton]
#bkgSamples = [DY,TTJets_SingleLeptonFromTbar,TTJets_SingleLeptonFromT,T_tWch,TBar_tWch,TToLep_tch,TBarToLep_tch,TToLep_sch,WWTo2L2Nu,WZTo2L2Q,WZTo3LNu,ZZTo2L2Nu,ZZTo2L2Q,ZZTo4L,WWW,WWZ,WZZ,ZZZ]

sigSamples = [TTJets_DiLepton,TTJets_SingleLeptonFromTbar,TTJets_SingleLeptonFromT,T_tWch,TBar_tWch,TToLep_tch,TBarToLep_tch,TToLep_sch,WWTo2L2Nu,WZTo3LNu,ZZTo2L2Nu]
bkgSamples = [DY,WZTo2L2Q,ZZTo2L2Q,ZZTo4L]

h_sig_num = ROOT.TH1F("h_sig_num","",1,0,2)
h_sig_den = ROOT.TH1F("h_sig_den","",1,0,2)
h_bkg_num = ROOT.TH1F("h_bkg_num","",1,0,2)
h_bkg_den = ROOT.TH1F("h_bkg_den","",1,0,2)

#presel = "(Sum$(lep_mediumMuonId>=1&&lep_pt>20&&abs(lep_miniRelIso)<0.2&&abs(lep_pdgId)==13)==2)&&nJet40==0)"
presel = "(Sum$(lep_mediumMuonId>=1&&lep_pt>20&&abs(lep_pdgId)==13)==2&&nJet40==0)"

sig_den = 0.
for s in sigSamples:
    s.chain.Draw("(1)>>h_sig_den",presel+"*xsec*genWeight/"+str(s.sumweight))
    sig_den += h_sig_den.Integral()

bkg_den = 0.
for s in bkgSamples:
    s.chain.Draw("(1)>>h_bkg_den",presel+"*xsec*genWeight/"+str(s.sumweight))
    bkg_den += h_bkg_den.Integral()


workingPoints = range(0,20)

metSigRes = {}



print "MET Significance discriminator"
print "{:>5}{:>10}{:>10}".format("WP","SigEff","BkgEff")
for i,wp in enumerate(workingPoints):
    sig_num = 0.
    for s in sigSamples:
        s.chain.Draw("(1)>>h_sig_num",'('+presel+"&&met_sig>"+str(wp)+")*xsec*genWeight/"+str(s.sumweight))
        sig_num += h_sig_num.Integral()
    bkg_num = 0.
    for s in bkgSamples:
        s.chain.Draw("(1)>>h_bkg_num",'('+presel+"&&met_sig>"+str(wp)+")*xsec*genWeight/"+str(s.sumweight))
        bkg_num += h_bkg_num.Integral()
    metSigRes[wp] = {"sigEff":sig_num/sig_den, "bkgEff":bkg_num/bkg_den}
    print "{:5}{:11.3f}{:11.3f}".format(wp,sig_num/sig_den,bkg_num/bkg_den)


workingPoints = range(0,100,5)


metRes = {}

print
print "MET pT discriminator"
print "{:>5}{:>10}{:>10}".format("WP","SigEff","BkgEff")
for wp in workingPoints:
    sig_num = 0.
    for s in sigSamples:
        s.chain.Draw("(1)>>h_sig_num",'('+presel+"&&met_pt>"+str(wp)+")*xsec*genWeight/"+str(s.sumweight))
        sig_num += h_sig_num.Integral()
    bkg_num = 0.
    for s in bkgSamples:
        s.chain.Draw("(1)>>h_bkg_num",'('+presel+"&&met_pt>"+str(wp)+")*xsec*genWeight/"+str(s.sumweight))
        bkg_num += h_bkg_num.Integral()
    metRes[wp] = {"sigEff":sig_num/sig_den, "bkgEff":bkg_num/bkg_den}
    print "{:5}{:11.3f}{:11.3f}".format(wp,sig_num/sig_den,bkg_num/bkg_den)


can = ROOT.TCanvas('can','can',700,700)
#can.SetLogy()
workingPoints = range(0,20)
metSigRoc = ROOT.TGraph()
for i,wp in enumerate(workingPoints):
    metSigRoc.SetPoint(i,metSigRes[wp]["sigEff"],1-metSigRes[wp]["bkgEff"])
metSigRoc.SetPoint(len(workingPoints),0,1)
metSigRoc.SetLineColor(ROOT.kBlue)
metSigRoc.GetYaxis().SetTitle("Background rejection")
metSigRoc.GetXaxis().SetTitle("Signal efficiency")
metSigRoc.GetYaxis().SetRangeUser(0.5,1)
metSigRoc.SetFillColor(0)
metSigRoc.Draw()

workingPoints = range(0,100,5)
metPtRoc = ROOT.TGraph()
for i,wp in enumerate(workingPoints):
    metPtRoc.SetPoint(i,metRes[wp]["sigEff"],1-metRes[wp]["bkgEff"])
metPtRoc.SetPoint(len(workingPoints),0,1)
metPtRoc.SetLineColor(ROOT.kRed)
metPtRoc.SetFillColor(0)

metPtRoc.Draw("same")

leg = ROOT.TLegend(0.15,0.3,0.35,0.4)
leg.SetFillColor(ROOT.kWhite)
leg.SetShadowColor(ROOT.kWhite)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)
leg.AddEntry(metSigRoc,'MET Significance')
leg.AddEntry(metPtRoc,'MET')
leg.Draw()

can.Print("/afs/hephy.at/user/d/dspitzbart/www/METSig/roc_njetEq0_lep.png")
