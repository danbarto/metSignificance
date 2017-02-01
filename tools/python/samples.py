import ROOT
import subprocess

class samples:
  def __init__(self, name='default', treeName='events', isData=False, subGroup = None, xsec=1., nEvents=1, rootfiles='', skimreport=''):
    self.name     = name
    self.isData   = isData
    self.subGroup = subGroup
    self.xsec     = xsec
    self.chain    = ROOT.TChain(treeName)
    self.weight   = 1.
    self.nEvents  = nEvents
    self.nEntries = 0
    self.targetLumi = 1.
    self.chain.Add(rootfiles)
    self.skimreport = skimreport
    self.sumweight = 1
    if skimreport:
        const = 'All Events' if self.isData else 'Sum Weights'
        logfile = self.skimreport
        line = [x for x in subprocess.check_output(["cat", logfile]).split('\n') if x.count(const)]
        assert len(line)==1,"Didn't find normalization constant '%s' in  number in file %s"%(const, logfile)
        self.sumweight = float(line[0].split()[2])

#  def defineSample(self, name, treeName='events', isData=False, subGroup = None, xsec=1.):
#    self.name     = name
#    self.isData   = isData
#    self.subGroup = subGroup
#    self.xsec     = xsec
#    self.chain    = ROOT.TChain(treeName)
  def readSkimReport(self):
    const = 'All Events' if self.isData else 'Sum Weights'
    logfile = self.skimreport
    line = [x for x in subprocess.check_output(["cat", logfile]).split('\n') if x.count(const)]
    assert len(line)==1,"Didn't find normalization constant '%s' in  number in file %s"%(const, logfile)
    self.sumweight = float(line[0].split()[2])
  
  def calculateWeight(self):
    self.nEntries   = self.chain.GetEntries()
    self.weight     = self.xsec * self.targetLumi / self.nEvents

  def setTargetLumi(self,lumi):
    self.targetLumi = lumi
#    self.weight     = self.xsec * self.targetLumi / self.nEvents

  def addToChain(self,rootfile):
    self.chain.Add(rootfile)
