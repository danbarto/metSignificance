# load all samples in here to keep tune and plot scripts clean
from samples import samples as sample
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

basedir = '/afs/hephy.at/data/dspitzbart01/MetSig/tuples/'

# load samples

DY          = sample('DY',            xsec=xsec_dyjetstoll, nEvents=nevts_dyjetstoll, subGroup='Zmumu', rootfiles=basedir+'MC_ICHEP/DYJets/*.root')
TTJets      = sample('ttJets',        xsec=xsec_ttjets,     nEvents=nevts_ttjets,     subGroup='top',   rootfiles=basedir+'MC_ICHEP/TTJets/*.root')
WW          = sample('WW',            xsec=xsec_ww,         nEvents=nevts_ww,         subGroup='EWK',   rootfiles=basedir+'MC_ICHEP/WW/*.root')
WZ          = sample('WZ',            xsec=xsec_wz,         nEvents=nevts_wz,         subGroup='EWK',   rootfiles=basedir+'MC_ICHEP/WZ/*.root')
ZZ          = sample('ZZ',            xsec=xsec_zz,         nEvents=nevts_zz,         subGroup='EWK',   rootfiles=basedir+'MC_ICHEP/ZZ/*.root')
WJets       = sample('WJets',         xsec=xsec_wjetstolnu, nEvents=nevts_wjetstolnu, subGroup='EWK',   rootfiles=basedir+'MC_ICHEP/WJets/*.root')
ST_top      = sample('sTop_top',      xsec=xsec_st_top,     nEvents=nevts_st_top,     subGroup='top',   rootfiles=basedir+'MC_ICHEP/ST_top/*.root')
ST_antitop  = sample('sTop_antitop',  xsec=xsec_st_antitop, nEvents=nevts_st_antitop, subGroup='top',   rootfiles=basedir+'MC_ICHEP/ST_antitop/*.root')


allMCSamples = [TTJets,ST_top,ST_antitop,WW,WZ,ZZ,WJets,DY]

for s in allMCSamples:
  s.setTargetLumi(12900)
  s.calculateWeight()

#data = sample('data', xsec=1, subGroup='Data', isData=True, rootfiles='../Ntuples/Zmumu/20160708/Data/*.root')
#data = sample('data', xsec=1, subGroup='Data', isData=True, rootfiles=basedir+'crab_MetSig_Data_2016G_promptReco/*.root')
data2016B = sample('data', xsec=1, subGroup='Data', isData=True, rootfiles=basedir+'crab_MetSig_Data_2016B/*.root')
data2016C = sample('data', xsec=1, subGroup='Data', isData=True, rootfiles=basedir+'crab_MetSig_Data_2016C/*.root')
data2016D = sample('data', xsec=1, subGroup='Data', isData=True, rootfiles=basedir+'crab_MetSig_Data_2016D/*.root')
data2016G = sample('data', xsec=1, subGroup='Data', isData=True, rootfiles=basedir+'crab_MetSig_Data_2016G/*.root')

allICHEPSamples = [data2016B,data2016C,data2016D]
ICHEP = sample('data', xsec=1, subGroup='Data', isData=True)
for s in allICHEPSamples:
  ICHEP.chain.Add(s.chain)


