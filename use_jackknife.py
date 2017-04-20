import ROOT

from metSignificance.tools.samples import samples as sample
from loadAllSamples_aprilRerun import *
from eventlist import *
from likelihood import *

from metSignificance.tools.jackknife import *

import json

testSample = sample('test',            xsec=1, subGroup='Zmumu', rootfiles='ntuple.root')
presel = "(1)"

def cartesian(pt, phi):
    return (pt*math.cos(phi), pt*math.sin(phi))

def acartesian(pt, phi):
    return (abs(pt*math.cos(phi)), abs(pt*math.sin(phi)))

def getSig(jet_pt, jet_sigmapt, jet_phi, jet_sigmaphi, jet_eta, met_pt, met_phi, met_sumpt, args):
  cov_xx       = 0
  cov_xy       = 0
  cov_yy       = 0
  i = 0
  for j in jet_pt:
    j_pt = j
    j_phi = jet_phi[i]
    j_sigmapt = jet_sigmapt[i]
    j_sigmaphi = jet_sigmaphi[i]
    index = getBin(abs(jet_eta[i]))

    cj = math.cos(j_phi)
    sj = math.sin(j_phi)
    dpt = args[index] * j_pt * j_sigmapt
    dph =               j_pt * j_sigmaphi

    dpt *= dpt
    dph *= dph

    cov_xx += dpt*cj*cj + dph*sj*sj
    cov_xy += (dpt-dph)*cj*sj
    cov_yy += dph*cj*cj + dpt*sj*sj

    i += 1

  # unclustered energy
  cov_tt = args[5]*args[5] + args[6]*args[6]*met_sumpt
  cov_xx += cov_tt
  cov_yy += cov_tt

  det = cov_xx*cov_yy - cov_xy*cov_xy

  ncov_xx =  cov_yy / det
  ncov_xy = -cov_xy / det
  ncov_yy =  cov_xx / det

  met_x = met_pt * math.cos(met_phi)
  met_y = met_pt * math.sin(met_phi)

  sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy
  return sig


def getSig2(jet_pt, jet_sigmapt, jet_phi, jet_sigmaphi, jet_eta, met_pt, met_phi, met_sumpt, args):
  cov_xx       = 0
  cov_xy       = 0
  cov_yy       = 0
  i = 0
  for j in jet_pt:
    j_pt = j
    j_phi = jet_phi[i]
    j_sigmapt = jet_sigmapt[i]
    j_sigmaphi = jet_sigmaphi[i]
    index = getBin(abs(jet_eta[i]))

    cj = math.cos(j_phi)
    sj = math.sin(j_phi)
    dpt = args[index] * j_pt * j_sigmapt
    dph =               j_pt * j_sigmaphi

    dpt *= dpt
    dph *= dph

    cov_xx += dpt*cj*cj + dph*sj*sj
    cov_xy += (dpt-dph)*cj*sj
    cov_yy += dph*cj*cj + dpt*sj*sj

    i += 1

  # unclustered energy
  cov_tt = met_sumpt
  cov_xx += cov_tt
  cov_yy += cov_tt

  det = cov_xx*cov_yy - cov_xy*cov_xy

  ncov_xx =  cov_yy / det
  ncov_xy = -cov_xy / det
  ncov_yy =  cov_xx / det

  met_x = met_pt * math.cos(met_phi)
  met_y = met_pt * math.sin(met_phi)

  sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy
  return sig

def getSig3(jet_pt, jet_sigmapt, jet_phi, jet_sigmaphi, jet_eta, met_pt, met_phi, unc_cov, args):
  cov_xx       = 0
  cov_xy       = 0
  cov_yy       = 0
  i = 0
  for j in jet_pt:
    j_pt = j
    j_phi = jet_phi[i]
    j_sigmapt = jet_sigmapt[i]
    j_sigmaphi = jet_sigmaphi[i]
    index = getBin(abs(jet_eta[i]))

    cj = math.cos(j_phi)
    sj = math.sin(j_phi)
    dpt = args[index] * j_pt * j_sigmapt
    dph =               j_pt * j_sigmaphi

    dpt *= dpt
    dph *= dph

    cov_xx += dpt*cj*cj + dph*sj*sj
    cov_xy += (dpt-dph)*cj*sj
    cov_yy += dph*cj*cj + dpt*sj*sj

    i += 1

  # unclustered energy
  cov_xy += unc_cov[0][1]
  cov_xx += unc_cov[0][0]
  cov_yy += unc_cov[1][1]

  det = cov_xx*cov_yy - cov_xy*cov_xy

  ncov_xx =  cov_yy / det
  ncov_xy = -cov_xy / det
  ncov_yy =  cov_xx / det

  met_x = met_pt * math.cos(met_phi)
  met_y = met_pt * math.sin(met_phi)

  sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy
  return sig

with open('data/MC_2016BH_njet0p.txt', 'r') as paraMC:
  parMC = json.load(paraMC)

with open('data/data_2016BH_njet0p.txt', 'r') as paraData:
  parData = json.load(paraData)


f = "{:>10}"*4+"{:>12.1f}"*15
fs= "{:>10}"*4+"{:>12}"*15
print fs.format("Ncand","Ncand3", "Ncand2","Niter","frac(0)","max pT","sumpt","nvert","sqrt(sumpt)","jackknife(1)","jackknife(Ns)","C_xx","C_yy", "MET", "METsumpt","njet","METSig", "METSigU", "METSigU2")

#s = DY_M17
#s = TTJets_M17

s = testSample
el_MC   = eventlist( [s], presel, isData=False)

for i in range(20):
    s.chain.GetEntry(i)

    cands_vec = []
    for j,c in enumerate(s.chain.cand_pt):
        if c>0:
            cands_vec.append(cartesian(s.chain.cand_pt[j], s.chain.cand_phi[j]))
    n_vert                  = s.chain.nvertices
    n_vert_true             = s.chain.nvert_true
    met                     = s.chain.met_pt
    jets30                  = [ j for j in s.chain.jet_pt if j > 30 ]
    njet = len(jets30)
    cands       = [ c for c in s.chain.cand_pt ]
    cands_min5  = [ c for c in s.chain.cand_pt if c>5 ]
    cands_min4  = [ c for c in s.chain.cand_pt if c>4 ]
    cands_min3  = [ c for c in s.chain.cand_pt if c>3 ]
    cands_min2  = [ c for c in s.chain.cand_pt if c>2 ]
    cands_min1  = [ c for c in s.chain.cand_pt if c>1 ]
    cands_min0  = [ c for c in s.chain.cand_pt if (c>0 and c<1) ]
    #print len(cands_min5), len(cands_min4), len(cands_min3), len(cands_min2)
    sig = getSig(s.chain.jet_pt, s.chain.jet_sigmapt, s.chain.jet_phi, s.chain.jet_sigmaphi, s.chain.jet_eta, s.chain.met_pt, s.chain.met_phi, s.chain.met_sumpt, parMC)    
    c_ii_1 = jackknife(cands,1)
    #c_ii_2 = jackknife(cands,2)
    Nd = min(5,int(math.sqrt(len(cands_min3)))+1)
    Nd = min(3,int(math.sqrt(len(cands_min2)))+1)
    n_iter = binomial(len(cands_min2),Nd)
    #print Nd, len(cands_min4)
    c_ii_5 = jackknife(cands_min2,Nd)
    cov = jackknifeMultiDim(cands_vec,1)
    sigU = getSig2(s.chain.jet_pt, s.chain.jet_sigmapt, s.chain.jet_phi, s.chain.jet_sigmaphi, s.chain.jet_eta, s.chain.met_pt, s.chain.met_phi, c_ii_5, parMC)
    sigU2 = getSig3(s.chain.jet_pt, s.chain.jet_sigmapt, s.chain.jet_phi, s.chain.jet_sigmaphi, s.chain.jet_eta, s.chain.met_pt, s.chain.met_phi, cov, parMC)

    
    print f.format(len(cands), len(cands_min3), len(cands_min2), n_iter, float(len(cands_min0))/len(cands), max(cands), sum(cands), n_vert, math.sqrt(sum(cands)), math.sqrt(c_ii_1), math.sqrt(c_ii_5), math.sqrt(cov[0][0]), math.sqrt(cov[1][1]), met, s.chain.met_sumpt, njet, sig, sigU,sigU2)

