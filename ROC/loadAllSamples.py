# load all samples in here to keep tune and plot scripts clean
from metSignificance.tools.samples import samples as sample
import math

basedir = '/afs/hephy.at/data/dspitzbart02/MetSig/METtuples_Spring16/'

# load samples

DY                          = sample('DY',                          treeName='METtree', subGroup='Zmumu', rootfiles=basedir+'DYJetsToLL_M50/METtree.root',    skimreport=basedir+'DYJetsToLL_M50/SkimReport.txt')
TTJets_DiLepton             = sample('ttJets_Dilepton',             treeName='METtree', subGroup='top',   rootfiles=basedir+'TTJets_DiLepton/METtree.root',   skimreport=basedir+'TTJets_DiLepton/SkimReport.txt')
TTJets_DiLepton_ext         = sample('ttJets_Dilepton',             treeName='METtree', subGroup='top',   rootfiles=basedir+'TTJets_DiLepton_ext/METtree.root',   skimreport=basedir+'TTJets_DiLepton_ext/SkimReport.txt')
TTJets_SingleLeptonFromTbar = sample('TTJets_SingleLeptonFromTbar', treeName='METtree', subGroup='top',   rootfiles=basedir+'TTJets_SingleLeptonFromTbar/METtree.root',   skimreport=basedir+'TTJets_SingleLeptonFromTbar/SkimReport.txt')
TTJets_SingleLeptonFromT    = sample('TTJets_SingleLeptonFromT',    treeName='METtree', subGroup='top',   rootfiles=basedir+'TTJets_SingleLeptonFromT/METtree.root',   skimreport=basedir+'TTJets_SingleLeptonFromT/SkimReport.txt')
T_tWch                      = sample('T_tWch',                      treeName='METtree', subGroup='top',   rootfiles=basedir+'T_tWch/METtree.root',   skimreport=basedir+'T_tWch/SkimReport.txt')
TBar_tWch                   = sample('TBar_tWch',                   treeName='METtree', subGroup='top',   rootfiles=basedir+'TBar_tWch/METtree.root',   skimreport=basedir+'TBar_tWch/SkimReport.txt')
TToLep_tch                  = sample('TToLep_tch',                  treeName='METtree', subGroup='top',   rootfiles=basedir+'TToLep_tch/METtree.root',   skimreport=basedir+'TToLep_tch/SkimReport.txt')
TBarToLep_tch               = sample('TBarToLep_tch',               treeName='METtree', subGroup='top',   rootfiles=basedir+'TBarToLep_tch/METtree.root',   skimreport=basedir+'TBarToLep_tch/SkimReport.txt')
TToLep_sch                  = sample('TToLep_sch',                  treeName='METtree', subGroup='top',   rootfiles=basedir+'TToLep_sch/METtree.root',   skimreport=basedir+'TToLep_sch/SkimReport.txt')
WWTo2L2Nu                   = sample('WWTo2L2Nu',                   treeName='METtree', subGroup='EWK',   rootfiles=basedir+'WWTo2L2Nu/METtree.root',   skimreport=basedir+'WWTo2L2Nu/SkimReport.txt')
WZTo2L2Q                    = sample('WZTo2L2Q',                    treeName='METtree', subGroup='EWK',   rootfiles=basedir+'WZTo2L2Q/METtree.root',   skimreport=basedir+'WZTo2L2Q/SkimReport.txt')
WZTo3LNu                    = sample('WZTo3LNu',                    treeName='METtree', subGroup='EWK',   rootfiles=basedir+'WZTo3LNu/METtree.root',   skimreport=basedir+'WZTo3LNu/SkimReport.txt')
ZZTo2L2Nu                   = sample('ZZTo2L2Nu',                   treeName='METtree', subGroup='EWK',   rootfiles=basedir+'ZZTo2L2Nu/METtree.root',   skimreport=basedir+'ZZTo2L2Nu/SkimReport.txt')
ZZTo2L2Q                    = sample('ZZTo2L2Q',                    treeName='METtree', subGroup='EWK',   rootfiles=basedir+'ZZTo2L2Q/METtree.root',   skimreport=basedir+'ZZTo2L2Q/SkimReport.txt')
ZZTo4L                      = sample('ZZTo4L',                      treeName='METtree', subGroup='EWK',   rootfiles=basedir+'ZZTo4L/METtree.root',   skimreport=basedir+'ZZTo4L/SkimReport.txt')
WWW                         = sample('WWW',                         treeName='METtree', subGroup='EWK',   rootfiles=basedir+'WWW/METtree.root',   skimreport=basedir+'WWW/SkimReport.txt')
WWZ                         = sample('WWZ',                         treeName='METtree', subGroup='EWK',   rootfiles=basedir+'WWZ/METtree.root',   skimreport=basedir+'WWZ/SkimReport.txt')
WZZ                         = sample('WZZ',                         treeName='METtree', subGroup='EWK',   rootfiles=basedir+'WZZ/METtree.root',   skimreport=basedir+'WZZ/SkimReport.txt')
ZZZ                         = sample('ZZZ',                         treeName='METtree', subGroup='EWK',   rootfiles=basedir+'ZZZ/METtree.root',   skimreport=basedir+'ZZZ/SkimReport.txt')


bkgSamples = [TTJets_DiLepton,TTJets_SingleLeptonFromTbar,TTJets_SingleLeptonFromT,T_tWch,TBar_tWch,TToLep_tch,TBarToLep_tch,TToLep_sch,WWTo2L2Nu,WZTo2L2Q,WZTo3LNu,ZZTo2L2Nu,ZZTo2L2Q,ZZTo4L,WWW,WWZ,WZZ,ZZZ]

#TTJets      = sample('ttJets',        treeName='METtree', subGroup='top',   rootfiles=basedir+'TTJets_DiLepton/*.root')
#WW          = sample('WW',            treeName='METtree', subGroup='EWK',   rootfiles=basedir+'/*.root')
#WZ          = sample('WZ',            treeName='METtree', subGroup='EWK',   rootfiles=basedir+'/*.root')
#ZZ          = sample('ZZ',            treeName='METtree', subGroup='EWK',   rootfiles=basedir+'/*.root')
#WJets       = sample('WJets',         treeName='METtree', subGroup='EWK',   rootfiles=basedir+'ets/*.root')
#ST_top      = sample('sTop_top',      treeName='METtree', subGroup='top',   rootfiles=basedir+'_top/*.root')
#ST_antitop  = sample('sTop_antitop',  treeName='METtree', subGroup='top',   rootfiles=basedir+'_antitop/*.root')




