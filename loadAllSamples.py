# load all samples in here to keep tune and plot scripts clean
from metSignificance.tools.samples import samples as sample
import math

# define xsecs (not included in samples for now)
xsec_dyjetstoll  = 1921.8*3
xsec_ttjets      = 831.76#*(3*0.108)**2
xsec_ww          = 118.7
xsec_wz          = 47.13
xsec_zz          = 16.523
xsec_wjetstolnu  = 61526.7
xsec_st_top      = 35.9
xsec_st_antitop  = 35.9

# define unskimmed number of events for weight calculation
nevts_dyjetstoll          = 49877138
nevts_ttjets              = 10259872
#nevts_ttjets              = 6058236
nevts_wjetstolnu          = 9908534
#nevts_wjetstolnu          = 28210360
nevts_ww                  = 993214
nevts_wz                  = 1000000
nevts_zz                  = 989312
nevts_st_top              = 998400
nevts_st_antitop          = 985000

basedir = '/afs/hephy.at/data/dspitzbart02/MetSig/tuples/'

# load samples

DY          = sample('DY',            xsec=xsec_dyjetstoll, nEvents=nevts_dyjetstoll, subGroup='Zmumu', rootfiles=basedir+'MC_ICHEP/DYJets/*.root')
TTJets      = sample('ttJets',        xsec=xsec_ttjets,     nEvents=nevts_ttjets,     subGroup='top',   rootfiles=basedir+'MC_ICHEP/TTJets/*.root')
WW          = sample('WW',            xsec=xsec_ww,         nEvents=nevts_ww,         subGroup='EWK',   rootfiles=basedir+'MC_ICHEP/WW/*.root')
WZ          = sample('WZ',            xsec=xsec_wz,         nEvents=nevts_wz,         subGroup='EWK',   rootfiles=basedir+'MC_ICHEP/WZ/*.root')
ZZ          = sample('ZZ',            xsec=xsec_zz,         nEvents=nevts_zz,         subGroup='EWK',   rootfiles=basedir+'MC_ICHEP/ZZ/*.root')
WJets       = sample('WJets',         xsec=xsec_wjetstolnu, nEvents=nevts_wjetstolnu, subGroup='EWK',   rootfiles=basedir+'MC_ICHEP/WJets/*.root')
ST_top      = sample('sTop_top',      xsec=xsec_st_top,     nEvents=nevts_st_top,     subGroup='top',   rootfiles=basedir+'MC_ICHEP/ST_top/*.root')
ST_antitop  = sample('sTop_antitop',  xsec=xsec_st_antitop, nEvents=nevts_st_antitop, subGroup='top',   rootfiles=basedir+'MC_ICHEP/ST_antitop/*.root')

nevts_M17_dyjetstoll          = 49144274
nevts_M17_ttjets              = 10139950
nevts_M17_wjetstolnu          = 24120319
nevts_M17_ww                  = 994012
nevts_M17_wz                  = 1000000
nevts_M17_zz                  = 990064
nevts_M17_st_top              = 6952830
nevts_M17_st_antitop          = 6933094

basedir_new = '/afs/hephy.at/data/dspitzbart02/MetSig/tune_tuples_M17/'

DY_M17          = sample('DY',            xsec=xsec_dyjetstoll, nEvents=nevts_M17_dyjetstoll, subGroup='Zmumu', rootfiles=basedir_new+'crab_MetSig_DYJets_Moriond17_v4/*.root')
TTJets_M17      = sample('ttJets',        xsec=xsec_ttjets,     nEvents=nevts_M17_ttjets,     subGroup='top',   rootfiles=basedir_new+'crab_MetSig_TTJets_Moriond17_v4/*.root')
WW_M17          = sample('WW',            xsec=xsec_ww,         nEvents=nevts_M17_ww,         subGroup='EWK',   rootfiles=basedir_new+'crab_MetSig_WW_Moriond17_v3/*.root')
WZ_M17          = sample('WZ',            xsec=xsec_wz,         nEvents=nevts_M17_wz,         subGroup='EWK',   rootfiles=basedir_new+'crab_MetSig_WZ_Moriond17_v4/*.root')
ZZ_M17          = sample('ZZ',            xsec=xsec_zz,         nEvents=nevts_M17_zz,         subGroup='EWK',   rootfiles=basedir_new+'crab_MetSig_ZZ_Moriond17_v4/*.root')
WJets_M17       = sample('WJets',         xsec=xsec_wjetstolnu, nEvents=nevts_M17_wjetstolnu, subGroup='EWK',   rootfiles=basedir_new+'crab_MetSig_WJets_Moriond17_v3/*.root')
ST_top_M17      = sample('sTop_top',      xsec=xsec_st_top,     nEvents=nevts_M17_st_top,     subGroup='top',   rootfiles=basedir_new+'crab_MetSig_ST_tW_top_Moriond17_v3/*.root')
ST_antitop_M17  = sample('sTop_antitop',  xsec=xsec_st_antitop, nEvents=nevts_M17_st_antitop, subGroup='top',   rootfiles=basedir_new+'crab_MetSig_ST_tW_antitop_Moriond17_v3/*.root')

#allMCSamples = [TTJets,ST_top,ST_antitop,WW,WZ,ZZ,WJets,DY]
#allMCSamples = [TTJets_M17, ST_top_M17, ST_antitop_M17, WW_M17, WZ_M17,ZZ_M17,WJets_M17,DY_M17]
allMCSamples = [TTJets_M17, WZ_M17, ZZ_M17, DY_M17]

for s in allMCSamples:
  s.setTargetLumi(36000)
  s.calculateWeight()

#data = sample('data', xsec=1, subGroup='Data', isData=True, rootfiles='../Ntuples/Zmumu/20160708/Data/*.root')
#data = sample('data', xsec=1, subGroup='Data', isData=True, rootfiles=basedir+'crab_MetSig_Data_2016G_promptReco/*.root')
data2016B = sample('data', xsec=1, subGroup='Data', isData=True, rootfiles=basedir_new+'DoubleMuon/crab_MetSig_DataMoriond17_2016B_v3/*.root')
data2016C = sample('data', xsec=1, subGroup='Data', isData=True, rootfiles=basedir_new+'DoubleMuon/crab_MetSig_DataMoriond17_2016C_v3/*.root')
data2016D = sample('data', xsec=1, subGroup='Data', isData=True, rootfiles=basedir_new+'DoubleMuon/crab_MetSig_DataMoriond17_2016D_v3/*.root')
data2016E = sample('data', xsec=1, subGroup='Data', isData=True, rootfiles=basedir_new+'DoubleMuon/crab_MetSig_DataMoriond17_2016E_v3/*.root')
data2016F = sample('data', xsec=1, subGroup='Data', isData=True, rootfiles=basedir_new+'DoubleMuon/crab_MetSig_DataMoriond17_2016F_v3/*.root')
data2016G = sample('data', xsec=1, subGroup='Data', isData=True, rootfiles=basedir_new+'DoubleMuon/crab_MetSig_DataMoriond17_2016G_v3/*.root')
data2016H2= sample('data', xsec=1, subGroup='Data', isData=True, rootfiles=basedir_new+'DoubleMuon/crab_MetSig_DataMoriond17_2016H-v2_v3/*.root')
data2016H3= sample('data', xsec=1, subGroup='Data', isData=True, rootfiles=basedir_new+'DoubleMuon/crab_MetSig_DataMoriond17_2016H-v3_v3/*.root')



#allICHEPSamples = [data2016B,data2016C,data2016D]
allDataSamples = [data2016B,data2016C,data2016D,data2016E,data2016F,data2016G,data2016H2,data2016H3]
#ICHEP = sample('data', xsec=1, subGroup='Data', isData=True)
data = sample('data', xsec=1, subGroup='Data', isData=True)
for s in allDataSamples:
  data.chain.Add(s.chain)

#for s in allICHEPSamples:
#  ICHEP.chain.Add(s.chain)


