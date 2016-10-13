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

# load samples

DY          = sample('DY', xsec=xsec_dyjetstoll, subGroup='Zmumu', rootfiles='../Ntuples/Zmumu/20160708/DYJetsToLL/*.root')
TTJets      = sample('ttJets', xsec=xsec_ttjets, subGroup='top', rootfiles='../Ntuples/Zmumu/20160708/TTJets/*.root')
WW          = sample('WW', xsec=xsec_ww, subGroup='EWK', rootfiles='../Ntuples/Zmumu/20160708/WW/*.root')
WZ          = sample('WZ', xsec=xsec_wz, subGroup='EWK', rootfiles='../Ntuples/Zmumu/20160708/WZ/*.root')
ZZ          = sample('ZZ', xsec=xsec_zz, subGroup='EWK', rootfiles='../Ntuples/Zmumu/20160708/ZZ/*.root')
WJets       = sample('WJets', xsec=xsec_wjetstolnu, subGroup='EWK', rootfiles='../Ntuples/Zmumu/20160708/WJetsToLNu/*.root')
ST_top      = sample('sTop_top', xsec=xsec_st_top, subGroup='top', rootfiles='../Ntuples/Zmumu/20160708/ST_tW_top/*.root')
ST_antitop  = sample('sTop_antitop', xsec=xsec_st_antitop, subGroup='top', rootfiles='../Ntuples/Zmumu/20160708/ST_tW_antitop/*.root')

allMCSamples = [DY,TTJets,WW,WZ,ZZ,WJets,ST_top,ST_antitop]

for s in allMCSamples:
  s.setTargetLumi(12.9)
  s.calculateWeight()

data = sample('data', xsec=1, subGroup='Data', isData=True, rootfiles='../Ntuples/Zmumu/20160708/data/*.root')

