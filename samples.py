import ROOT

class samples:
  def __init__(self, name='default', treeName='events', isData=False, subGroup = None, xsec=1., nEvents=1, rootfiles=''):
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

#  def defineSample(self, name, treeName='events', isData=False, subGroup = None, xsec=1.):
#    self.name     = name
#    self.isData   = isData
#    self.subGroup = subGroup
#    self.xsec     = xsec
#    self.chain    = ROOT.TChain(treeName)

  def calculateWeight(self):
    self.nEntries   = self.chain.GetEntries()
    self.weight     = self.xsec * self.targetLumi / self.nEvents

  def setTargetLumi(self,lumi):
    self.targetLumi = lumi
#    self.weight     = self.xsec * self.targetLumi / self.nEvents

  def addToChain(self,rootfile):
    self.chain.Add(rootfile)
