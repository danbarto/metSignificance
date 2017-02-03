from metSig import *

mS = metSig()

c = ROOT.TChain("METtree")
c.Add("/afs/hephy.at/work/d/dspitzbart/MET/test/CMSSW_8_0_25/src/CMGTools/ObjectStudies/cfg/Test/DYJetsToLL_M50_7/METtree.root")

for ev in range(c.GetEntries()):
    print
    print ev
    c.GetEntry(ev)
    print "{:10.2f}{:10.2f}{:15.4f}{:10.0f}".format(c.met_pt,c.met_phi,c.met_sig,len(c.jet_pt))
    mS.getMetSig(c)
    mS.getMetSig(c, jec='up')
    mS.getMetSig(c, jec='down')

del mS

