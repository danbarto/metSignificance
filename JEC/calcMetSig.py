import ROOT
import math
import os

# compile the JER reader
test = ROOT.JetCorrectorParameters()
JetResolutionFile = "$CMSSW_BASE/src/JetMETCorrections/Modules/src/JetResolution.cc+"
JetResolutionFile = os.path.expandvars(JetResolutionFile)
ROOT.gROOT.ProcessLine('.L '+JetResolutionFile)

import FWCore.ParameterSet.Config as cms
from RecoMET.METProducers.METSignificanceParams_cfi import *

# Load the resolutions for data and MC, pt and phi
JERFiles = "$CMSSW_BASE/src/metSignificance/tools/data/JRDatabase/textFiles/Spring16_25nsV10_MC/"
JERFiles = os.path.expandvars(JERFiles)
res_pt  = ROOT.JME.JetResolution(JERFiles+"Spring16_25nsV6_MC_PtResolution_AK4PFchs.txt") 
res_phi = ROOT.JME.JetResolution(JERFiles+"Spring16_25nsV6_MC_PhiResolution_AK4PFchs.txt")

jetThreshold = METSignificanceParams.getParameter('jetThreshold').value()
etabins = METSignificanceParams.getParameter('jeta').value()
nEtabins = len(etabins)
jet_args = METSignificanceParams.getParameter('jpar').value()
unc_args = METSignificanceParams.getParameter('pjpar').value()

# to re-calculate the met-sig on CMG level
def calcMetSig(tree, jec='', met_pt=False, met_phi=False, sumpt=False, unclustered=''):
    if 'up' in jec.lower():
        jec = '_jecUp'
    if 'down' in jec.lower():
        jec = '_jecDown'
    
    cov_xx       = 0
    cov_xy       = 0
    cov_yy       = 0
    
    #load stuff
    jet_pt  = getattr(tree, 'jet'+jec+'_pt')
    jet_phi = getattr(tree, 'jet'+jec+'_phi')
    jet_eta = getattr(tree, 'jet'+jec+'_eta')
    
    lep_pt  = getattr(tree, 'lep_pt')
    
    if not met_pt:
        met_pt      = getattr(tree, 'met'+jec+'_pt')
    if not met_phi:
        met_phi     = getattr(tree, 'met'+jec+'_phi')
    met_sumEt   = getattr(tree, 'met_sumEt')

    rho = tree.rho
    
    nJets = len([x for x in jet_pt if x>jetThreshold])
    
    for i,j in enumerate(jet_pt):
        jet_para = ROOT.JME.JetParameters()
        j_pt    = jet_pt[i]
        j_phi   = jet_phi[i]
        j_eta   = jet_eta[i]
        j_aeta  = abs(j_eta)
        
        print j_pt, j_phi, j_eta
        jet_para.setJetEta(j_eta).setJetPt(j_pt).setRho(rho)
        j_sigmapt   = res_pt.getResolution(jet_para)
        j_sigmaphi  = res_phi.getResolution(jet_para)
        print "sigma phi", j_sigmaphi
        print "alt sigma phi", j_pt*j_sigmaphi
        
        etabin = 0
        for n in range(nEtabins):
            if j_aeta < etabins[n]: break
            etabin += 1
        print j_aeta, etabin

        cj = math.cos(j_phi)
        sj = math.sin(j_phi)
        dpt = jet_args[etabin] * j_pt * j_sigmapt
        #dph =                    j_pt * j_sigmaphi
        dph =                    j_sigmaphi

        dpt *= dpt
        dph *= dph

        cov_xx += dpt*cj*cj + dph*sj*sj
        cov_xy += (dpt-dph)*cj*sj
        cov_yy += dph*cj*cj + dpt*sj*sj
        

    # unclustered energy
    if not sumpt:
        sumpt = met_sumEt - sum([x for x in lep_pt])
    print sumpt
    print cov_xx
    print cov_yy
    cov_tt = unc_args[0]*unc_args[0] + unc_args[1]*unc_args[1]*sumpt
    cov_xx += cov_tt
    cov_yy += cov_tt

    det = cov_xx*cov_yy - cov_xy*cov_xy

    ncov_xx =  cov_yy / det
    ncov_xy = -cov_xy / det
    ncov_yy =  cov_xx / det
    
    #ncov_xy = 0
    print "normalized", ncov_xx, ncov_xy, ncov_yy
    print "test", cov_yy / math.sqrt(det)
    met_x = met_pt * math.cos(met_phi)
    met_y = met_pt * math.sin(met_phi)

    det = det
    print "Det",det
    sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy
    print sig
    return sig
