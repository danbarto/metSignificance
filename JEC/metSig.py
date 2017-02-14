import ROOT
import math
import os

import FWCore.ParameterSet.Config as cms
from RecoMET.METProducers.METSignificanceParams_cfi import *

JERFiles    = "$CMSSW_BASE/src/metSignificance/tools/data/JRDatabase/textFiles/Spring16_25nsV10_MC/"
JERFilePt   = "Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt"
JERFilePhi  = "Spring16_25nsV10_MC_PhiResolution_AK4PFchs.txt"
JERFileSF   = "Spring16_25nsV10_MC_SF_AK4PFchs.txt"

class metSig:

    def __init__(self):
        test = ROOT.JetCorrectorParameters()
        self.JetResolutionFile = "$CMSSW_BASE/src/JetMETCorrections/Modules/src/JetResolution.cc+"
        self.JetResolutionFile = os.path.expandvars(self.JetResolutionFile)
        ROOT.gROOT.ProcessLine('.L '+self.JetResolutionFile)
        
        # Load the resolutions for data and MC, pt and phi
        self.JERFiles = os.path.expandvars(JERFiles)
        self.res_pt  = ROOT.JME.JetResolution(self.JERFiles + JERFilePt)
        self.res_phi = ROOT.JME.JetResolution(self.JERFiles + JERFilePhi)
        self.jer_SF  = ROOT.JME.JetResolutionScaleFactor(self.JERFiles + JERFileSF)

        self.jetThreshold = METSignificanceParams.getParameter('jetThreshold').value()
        self.etabins = METSignificanceParams.getParameter('jeta').value()
        self.nEtabins = len(self.etabins)
        self.jet_args = METSignificanceParams.getParameter('jpar').value()
        self.unc_args = METSignificanceParams.getParameter('pjpar').value()
    
    def getJetCovariance(self, event, jec=''):
            
        if 'jecup' in jec.lower():
            jec_ = '_jecUp'
        elif 'jecdown' in jec.lower():
            jec_ = '_jecDown'
        else: jec_ = ''
        
        cov_xx       = 0
        cov_xy       = 0
        cov_yy       = 0

        #load stuff
        jet_pt  = getattr(event, 'jet'+jec_+'_pt')
        jet_phi = getattr(event, 'jet'+jec_+'_phi')
        jet_eta = getattr(event, 'jet'+jec_+'_eta')
        rho     = event.rho
        nJets   = len([x for x in jet_pt if x>self.jetThreshold])
        
        for i,j in enumerate(jet_pt):
            jet_para = ROOT.JME.JetParameters()
            j_pt    = jet_pt[i]
            j_phi   = jet_phi[i]
            j_eta   = jet_eta[i]
            j_aeta  = abs(j_eta)

            jet_para.setJetEta(j_eta).setJetPt(j_pt).setRho(rho)
            j_sigmapt   = self.res_pt.getResolution(jet_para)
            j_sigmaphi  = self.res_phi.getResolution(jet_para)
            j_SF        = self.jer_SF.getScaleFactor(jet_para)

            etabin = 0
            for n in range(self.nEtabins):
                if j_aeta < self.etabins[n]: break
                etabin += 1

            cj = math.cos(j_phi)
            sj = math.sin(j_phi)
            j_SF = 1
            dpt = self.jet_args[etabin] * j_pt * j_sigmapt * j_SF
            dph =                         j_pt * j_sigmaphi
            #dph =                    j_sigmaphi

            dpt *= dpt
            dph *= dph

            cov_xx += dpt*cj*cj + dph*sj*sj
            cov_xy += (dpt-dph)*cj*sj
            cov_yy += dph*cj*cj + dpt*sj*sj
        
        return cov_xx, cov_xy, cov_yy
        
    def getSumPt(self, event):
        cov_xx, cov_xy, cov_yy = self.getJetCovariance(event)
        sum_pt = (event.met_sig_Cxx - cov_xx - self.unc_args[0]**2)/self.unc_args[1]**2
        return sum_pt

    def getUnclusteredMetSig(self, event):
        met_pt      = getattr(event, 'met_pt')
        met_phi     = getattr(event, 'met_phi')

        sum_pt = self.getSumPt(event)
        cov_xx, cov_xy, cov_yy = self.getJetCovariance(event, jec='')
        cov_uncl = self.unc_args[0]**2 + self.unc_args[1]**2 * sum_pt
        cov_xx = cov_uncl
        cov_yy = cov_uncl
        cov_xy = 0.
        det = cov_xx * cov_yy - cov_xy*cov_xy
        ncov_xx =  cov_yy / det
        ncov_xy = -cov_xy / det
        ncov_yy =  cov_xx / det

        met_x = met_pt * math.cos(met_phi)
        met_y = met_pt * math.sin(met_phi)

        sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy
        return sig
       
    def getMetSig(self, event, jec=''):

        if 'jecup' in jec.lower():
            jec_ = '_jecUp'
        elif 'jecdown' in jec.lower():
            jec_ = '_jecDown'
        else: jec_=''
        
        met_pt      = getattr(event, 'met'+jec_+'_pt')
        met_phi     = getattr(event, 'met'+jec_+'_phi')
        
        sum_pt = self.getSumPt(event)
        if 'sumptup' in jec.lower() and sum_pt>0:
            sum_pt += math.sqrt(sum_pt)
        if 'sumptdown' in jec.lower() and sum_pt>0:
            sum_pt -= math.sqrt(sum_pt)
        cov_xx, cov_xy, cov_yy = self.getJetCovariance(event, jec=jec)
        cov_uncl = self.unc_args[0]**2 + self.unc_args[1]**2 * sum_pt
        cov_xx += cov_uncl
        cov_yy += cov_uncl

        det = cov_xx * cov_yy - cov_xy*cov_xy

        ncov_xx =  cov_yy / det
        ncov_xy = -cov_xy / det
        ncov_yy =  cov_xx / det
        
        met_x = met_pt * math.cos(met_phi)
        met_y = met_pt * math.sin(met_phi)
        
        sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy
        
        #unc_sig = self.getUnclusteredMetSig(event)
        #print "{:10}{:10.2f}{:10.2f}{:10.2f}{:15.4f}{:15.4f}{:10.1f}".format(jec,met_pt, met_phi, sum_pt, sig, unc_sig, det)
        
        return sig
