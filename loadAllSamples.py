# load all samples in here to keep tune and plot scripts clean
from samples import samples as sample
import math

# define xsecs (not included in samples for now)
xsec_dyjetstoll  = 1921.8*3
xsec_ttjets      = 831.76*(3*0.108)**2
xsec_ww          = 118.7
xsec_wz          = 46.74
xsec_zz          = 15.99
xsec_wjetstolnu  = 20508.90*3
xsec_st_top      = 35.9
xsec_st_antitop  = 35.9

# define unskimmed number of events for weight calculation
nevts_dyjetstoll          = 49877138
nevts_ttjets              = 6058236
nevts_wjetstolnu          = 28210360
nevts_ww                  = 993214
nevts_wz                  = 1000000
nevts_zz                  = 989312
nevts_st_top              = 998400
nevts_st_antitop          = 985000



# load samples

DY          = sample('DY',            xsec=xsec_dyjetstoll, nEvents=nevts_dyjetstoll, subGroup='Zmumu', rootfiles='../Ntuples/Zmumu/20160708/DYJetsToLL/*.root')
TTJets      = sample('ttJets',        xsec=xsec_ttjets,     nEvents=nevts_ttjets,     subGroup='top',   rootfiles='../Ntuples/Zmumu/20160708/TTJets/*.root')
WW          = sample('WW',            xsec=xsec_ww,         nEvents=nevts_ww,         subGroup='EWK',   rootfiles='../Ntuples/Zmumu/20160708/WW/*.root')
WZ          = sample('WZ',            xsec=xsec_wz,         nEvents=nevts_wz,         subGroup='EWK',   rootfiles='../Ntuples/Zmumu/20160708/WZ/*.root')
ZZ          = sample('ZZ',            xsec=xsec_zz,         nEvents=nevts_zz,         subGroup='EWK',   rootfiles='../Ntuples/Zmumu/20160708/ZZ/*.root')
WJets       = sample('WJets',         xsec=xsec_wjetstolnu, nEvents=nevts_wjetstolnu, subGroup='EWK',   rootfiles='../Ntuples/Zmumu/20160708/WJetsToLNu/*.root')
ST_top      = sample('sTop_top',      xsec=xsec_st_top,     nEvents=nevts_st_top,     subGroup='top',   rootfiles='../Ntuples/Zmumu/20160708/ST_tW_top/*.root')
ST_antitop  = sample('sTop_antitop',  xsec=xsec_st_antitop, nEvents=nevts_st_antitop, subGroup='top',   rootfiles='../Ntuples/Zmumu/20160708/ST_tW_antitop/*.root')

allMCSamples = [DY,TTJets,WW,WZ,ZZ,WJets,ST_top,ST_antitop]

for s in allMCSamples:
  s.setTargetLumi(12.9)
  s.calculateWeight()

data = sample('data', xsec=1, subGroup='Data', isData=True, rootfiles='../Ntuples/Zmumu/20160708/data/*.root')

