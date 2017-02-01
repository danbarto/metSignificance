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

sig_den = 0.
for s in sigSamples:
    s.chain.Draw("(1)>>h_sig_den","(1)*xsec*genWeight/"+str(s.sumweight))
    sig_den += h_sig_den.Integral()

bkg_den = 0.
for s in bkgSamples:
    s.chain.Draw("(1)>>h_bkg_den","(1)*xsec*genWeight/"+str(s.sumweight))
    bkg_den += h_bkg_den.Integral()


workingPoints = range(0,20)

metSigRes = {}

print "MET Significance discriminator"
print "{:>5}{:>10}{:>10}".format("WP","SigEff","BkgEff")
for wp in workingPoints:
    sig_num = 0.
    for s in sigSamples:
        s.chain.Draw("(1)>>h_sig_num","(met_sig>"+str(wp)+")*xsec*genWeight/"+str(s.sumweight))
        sig_num += h_sig_num.Integral()
    bkg_num = 0.
    for s in bkgSamples:
        s.chain.Draw("(1)>>h_bkg_num","(met_sig>"+str(wp)+")*xsec*genWeight/"+str(s.sumweight))
        bkg_num += h_bkg_num.Integral()
    metSigRes["wp"] = {"sigEff":sig_num/sig_den, "bkgEff":bkg_num/bkg_den}
    print "{:5}{:10.2f}{:10.2f}".format(wp,sig_num/sig_den,bkg_num/bkg_den)

workingPoints = range(0,100,5)


metRes = {}

print
print "MET pT discriminator"
print "{:>5}{:>10}{:>10}".format("WP","SigEff","BkgEff")
for wp in workingPoints:
    sig_num = 0.
    for s in sigSamples:
        s.chain.Draw("(1)>>h_sig_num","(met_pt>"+str(wp)+")*xsec*genWeight/"+str(s.sumweight))
        sig_num += h_sig_num.Integral()
    bkg_num = 0.
    for s in bkgSamples:
        s.chain.Draw("(1)>>h_bkg_num","(met_pt>"+str(wp)+")*xsec*genWeight/"+str(s.sumweight))
        bkg_num += h_bkg_num.Integral()
    metRes["wp"] = {"sigEff":sig_num/sig_den, "bkgEff":bkg_num/bkg_den}
    print "{:5}{:10.2f}{:10.2f}".format(wp,sig_num/sig_den,bkg_num/bkg_den)

