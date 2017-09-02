from __future__ import division
from heppy.framework.analyzer import Analyzer
from heppy.statistics.tree import Tree
from heppy.analyzers.ntuple import *
from heppy.particles.tlv.resonance import Resonance2 as Resonance

import math
import ROOT
from ROOT import *
import collections
from array import array

#For TMVA >>>>>>>>>>>>>>>>>>>>>
ROOT.gROOT.ProcessLine('.L /afs/cern.ch/user/r/rasmith/fcc/heppy/FCChhAnalyses/Zprime_tt/BDT_QCD.class.C+')

class TreeProducer(Analyzer):

    def beginLoop(self, setup):
        super(TreeProducer, self).beginLoop(setup)
        self.rootfile = TFile('/'.join([self.dirName,
                                        'tree.root']),
                              'recreate')
        self.tree = Tree( 'events', '')
        
	self.tree.var('weight', float)
	self.tree.var('missingET', float)
        self.tree.var('numberOfElectrons', int)
        self.tree.var('numberOfMuons', int)
	self.tree.var('numberOfFatJets', int)

	self.tree.var('Jet1_tau1', float)	
	self.tree.var('Jet1_tau2', float)
        self.tree.var('Jet1_tau3', float)
        self.tree.var('Jet2_tau1', float)
        self.tree.var('Jet2_tau2', float)
        self.tree.var('Jet2_tau3', float)
	self.tree.var('Jet1_tau32', float)
        self.tree.var('Jet1_tau31', float)
        self.tree.var('Jet1_tau21', float)
        self.tree.var('Jet2_tau32', float)
        self.tree.var('Jet2_tau31', float)
        self.tree.var('Jet2_tau21', float)

        bookParticle(self.tree, 'Jet1')
        bookParticle(self.tree, 'Jet2')
	bookParticle(self.tree, 'Jet3')
	bookParticle(self.tree, 'Jet4')

	self.tree.var('rapiditySeparation', float)
        self.tree.var('transverseMomentumAsymmetry', float)
	self.tree.var('topJetMassDifference', float)

	self.tree.var('Jet1_trimmedProngMaxPtRatio', float)
        self.tree.var('Jet1_trimmedProngMinPtRatio', float)
        self.tree.var('Jet2_trimmedProngMaxPtRatio', float)
        self.tree.var('Jet2_trimmedProngMinPtRatio', float)

	bookParticle(self.tree, 'softDroppedJet1')
	bookParticle(self.tree, 'softDroppedJet2')
	bookParticle(self.tree, 'softDroppedJet3')
	bookParticle(self.tree, 'softDroppedJet4')

	bookParticle(self.tree, 'smallJet1')
	bookParticle(self.tree, 'smallJet2')
	bookParticle(self.tree, 'smallJet3')
	bookParticle(self.tree, 'smallJet4')
        
	bookParticle(self.tree, 'Electron1')
	bookParticle(self.tree, 'Electron2')

	bookParticle(self.tree, 'Muon1')
	bookParticle(self.tree, 'Muon2')

        self.tree.var('BDTvariable_qcd', float)

	self.tree.var('zPrimeReconstructedMass', float)

    def process(self, event):
        self.tree.reset()
        jets = getattr(event, self.cfg_ana.fatjets)
	smallJets = getattr(event, self.cfg_ana.jets)
	electrons = getattr(event, self.cfg_ana.electrons)
	muons = getattr(event, self.cfg_ana.muons)

	Jet1_dR = 999
	Jet2_dR = 999
	if ( len(jets)>=2 ):
	    j1 = ROOT.TLorentzVector(); j2 = ROOT.TLorentzVector()
	    j1.SetPtEtaPhiE(jets[0].pt(), jets[0].eta(), jets[0].phi(), jets[0].e())
            j2.SetPtEtaPhiE(jets[1].pt(), jets[1].eta(), jets[1].phi(), jets[1].e())
	    if ( len(electrons)!=0 and len(muons)==0 ):
	    	e = ROOT.TLorentzVector()
	    	e.SetPtEtaPhiE(electrons[0].pt(), electrons[0].eta(), electrons[0].phi(), electrons[0].e())
	    	Jet1_dR = j1.DeltaR(e)
	    	Jet2_dR = j2.DeltaR(e)
	    if ( len(electrons)==0 and len(muons)!=0 ):
		m = ROOT.TLorentzVector()
            	m.SetPtEtaPhiE(muons[0].pt(), muons[0].eta(), muons[0].phi(), muons[0].e())
            	Jet1_dR = j1.DeltaR(m)
            	Jet2_dR = j2.DeltaR(m)
	    if ( len(electrons)!=0 and len(muons)!=0 ):
		isElectron = False; isMuon = False
		if ( electrons[0].pt() > muons[0].pt() ): isElectron = True
		else: isMuon = True
		l = ROOT.TLorentzVector()
		if isElectron: l.SetPtEtaPhiE(electrons[0].pt(), electrons[0].eta(), electrons[0].phi(), electrons[0].e())
		if isMuon: l.SetPtEtaPhiE(muons[0].pt(), muons[0].eta(), muons[0].phi(), muons[0].e())
	    	Jet1_dR = j1.DeltaR(l)
                Jet2_dR = j2.DeltaR(l)
	    
	
	if ( len(jets) >= 2 and Jet1_dR > 0.6 and Jet2_dR > 0.5 ):

	    self.tree.fill('weight' , event.weight )
            self.tree.fill('missingET', event.met.pt())
            self.tree.fill('numberOfElectrons', len(electrons))
            self.tree.fill('numberOfMuons', len(muons))
	    self.tree.fill('numberOfFatJets', len(jets))

	    self.tree.fill('rapiditySeparation', abs(jets[0].eta() - jets[1].eta()))
	    self.tree.fill('transverseMomentumAsymmetry', (jets[0].pt() - jets[1].pt())/(jets[0].pt() + jets[1].pt()))

	    self.tree.fill('Jet1_tau1' , jets[0].tau1 )
            self.tree.fill('Jet1_tau2' , jets[0].tau2 )
            self.tree.fill('Jet1_tau3' , jets[0].tau3 )
	    self.tree.fill('Jet2_tau1' , jets[1].tau1 )
            self.tree.fill('Jet2_tau2' , jets[1].tau2 )
            self.tree.fill('Jet2_tau3' , jets[1].tau3 )

	    Jet1_tau31 = -999.0
	    Jet1_tau21 = -999.0
	    Jet1_tau32 = -999.0
	    Jet2_tau31 = -999.0
            Jet2_tau21 = -999.0
            Jet2_tau32 = -999.0

            if (jets[0].tau1 != 0.0):
                Jet1_tau31 = jets[0].tau3/jets[0].tau1
                Jet1_tau21 = jets[0].tau2/jets[0].tau1 
            if (jets[0].tau2 != 0.0):
                Jet1_tau32 = jets[0].tau3/jets[0].tau2

	    if (jets[1].tau1 != 0.0):
                Jet2_tau31 = jets[1].tau3/jets[1].tau1
                Jet2_tau21 = jets[1].tau2/jets[1].tau1
	    if (jets[1].tau2 != 0.0):
                Jet2_tau32 = jets[1].tau3/jets[1].tau2

	    self.tree.fill('Jet1_tau31', Jet1_tau31)
	    self.tree.fill('Jet1_tau21', Jet1_tau21)
	    self.tree.fill('Jet1_tau32', Jet1_tau32)
	    self.tree.fill('Jet2_tau31', Jet2_tau31)
            self.tree.fill('Jet2_tau21', Jet2_tau21)
            self.tree.fill('Jet2_tau32', Jet2_tau32)

	    fillParticle(self.tree, 'Jet1', jets[0])
	    fillParticle(self.tree, 'Jet2', jets[1])

	    fillParticle(self.tree, 'softDroppedJet1', jets[0].subjetsSoftDrop[0])
	    fillParticle(self.tree, 'softDroppedJet2', jets[1].subjetsSoftDrop[0])

	    s1 = ROOT.TLorentzVector(); s2 = ROOT.TLorentzVector()
	    s1.SetPtEtaPhiE(jets[0].subjetsSoftDrop[0].pt(),
                            jets[0].subjetsSoftDrop[0].eta(),
                            jets[0].subjetsSoftDrop[0].phi(),
                            jets[0].subjetsSoftDrop[0].e())
            s2.SetPtEtaPhiE(jets[1].subjetsSoftDrop[0].pt(),
                            jets[1].subjetsSoftDrop[0].eta(),
                            jets[1].subjetsSoftDrop[0].phi(),
                            jets[1].subjetsSoftDrop[0].e())
	    self.tree.fill('topJetMassDifference', abs( s1.M() - s2.M() ))
	    
	    t1 = ROOT.TLorentzVector(); t2 = ROOT.TLorentzVector()
            t1.SetPtEtaPhiE(jets[0].pt(), jets[0].eta(), jets[0].phi(), jets[0].e())
            t2.SetPtEtaPhiE(jets[1].pt(), jets[1].eta(), jets[1].phi(), jets[1].e())
            self.tree.fill('zPrimeReconstructedMass', (t1+t2).M())            

	    Jet1_trimmedProngMaxPtRatio = -999.0
	    Jet1_trimmedProngMinPtRatio = -999.0
	    Jet2_trimmedProngMaxPtRatio = -999.0
	    Jet2_trimmedProngMinPtRatio = -999.0

	    if ( len(jets[0].subjetsTrimming) >=3 ):
                Jet1_trimmedProngMinPt = min(jets[0].subjetsTrimming[1].pt(),jets[0].subjetsTrimming[2].pt())
                Jet1_trimmedProngMaxPt = max(jets[0].subjetsTrimming[1].pt(),jets[0].subjetsTrimming[2].pt())
                Jet1_trimmedProngMaxPtRatio = Jet1_trimmedProngMaxPt/(jets[0].subjetsTrimming[0].pt())
                Jet1_trimmedProngMinPtRatio = Jet1_trimmedProngMinPt/(jets[0].subjetsTrimming[0].pt())

            if ( len(jets[1].subjetsTrimming) >=3 ):
                Jet2_trimmedProngMinPt = min(jets[1].subjetsTrimming[1].pt(),jets[1].subjetsTrimming[2].pt())
                Jet2_trimmedProngMaxPt = max(jets[1].subjetsTrimming[1].pt(),jets[1].subjetsTrimming[2].pt())
                Jet2_trimmedProngMaxPtRatio = Jet2_trimmedProngMaxPt/(jets[1].subjetsTrimming[0].pt())
                Jet2_trimmedProngMinPtRatio = Jet2_trimmedProngMinPt/(jets[1].subjetsTrimming[0].pt())

	    self.tree.fill('Jet1_trimmedProngMaxPtRatio', Jet1_trimmedProngMaxPtRatio)
	    self.tree.fill('Jet1_trimmedProngMinPtRatio', Jet1_trimmedProngMinPtRatio)
	    self.tree.fill('Jet2_trimmedProngMaxPtRatio', Jet2_trimmedProngMaxPtRatio)
	    self.tree.fill('Jet2_trimmedProngMinPtRatio', Jet2_trimmedProngMinPtRatio)

	    if ( len(jets) >= 3 ):
		fillParticle(self.tree, 'Jet3', jets[2])
		fillParticle(self.tree, 'softDroppedJet3', jets[2].subjetsSoftDrop[0])

	    if ( len(jets) >= 4 ):
                fillParticle(self.tree, 'Jet4', jets[3])
                fillParticle(self.tree, 'softDroppedJet4', jets[3].subjetsSoftDrop[0])

	    if ( len(smallJets) >= 1 ): fillParticle(self.tree, 'smallJet1', smallJets[0])
	    if ( len(smallJets) >= 2 ): fillParticle(self.tree, 'smallJet2', smallJets[1])
            if ( len(smallJets) >= 3 ): fillParticle(self.tree, 'smallJet3', smallJets[2])
            if ( len(smallJets) >= 4 ): fillParticle(self.tree, 'smallJet4', smallJets[3])    

	    if ( len(electrons) >=1 ): fillParticle(self.tree, 'Electron1', electrons[0])
	    if ( len(electrons) >=2 ): fillParticle(self.tree, 'Electron2', electrons[1])

	    if ( len(muons) >=1 ): fillParticle(self.tree, 'Muon1', muons[0])
	    if ( len(muons) >=2 ): fillParticle(self.tree, 'Muon2', muons[1])

	    ###################################
            #TMVA Stuff Starts!
	    ###################################

            ###################################
            #qcd Background >>>>>>>>>>>>>>>>>
            ###################################
	    
            varlist = [
				"Jet1_tau32",
                		"Jet2_tau32",
				"Jet1_tau21",
				"Jet2_tau21",
				"softDroppedJet1_m",
				"softDroppedJet2_m",
				"Jet1_trimmedProngMaxPtRatio",
				"Jet1_trimmedProngMinPtRatio",
				"Jet2_trimmedProngMaxPtRatio",
				"Jet2_trimmedProngMinPtRatio",
            ]            
            
	    inputs = ROOT.vector('string')()
            for v in varlist:
                inputs.push_back(v)

            mva = ROOT.ReadQCD(inputs)
	    values = ROOT.vector('double')()
	
	    values.push_back(Jet1_tau32)
	    values.push_back(Jet2_tau32)
	    values.push_back(Jet1_tau21)
	    values.push_back(Jet2_tau21)
            values.push_back(jets[0].subjetsSoftDrop[0].m())
            values.push_back(jets[1].subjetsSoftDrop[0].m())
	    values.push_back(Jet1_trimmedProngMaxPtRatio)
	    values.push_back(Jet1_trimmedProngMinPtRatio)
	    values.push_back(Jet2_trimmedProngMaxPtRatio)
	    values.push_back(Jet2_trimmedProngMinPtRatio)
   
            mva_value=mva.GetMvaValue(values)
            self.tree.fill('BDTvariable_qcd', mva_value)
            
            ###################################
            #TMVA Stuff Ends!
            ###################################
	    
            self.tree.tree.Fill()

    def write(self, setup):
        self.rootfile.Write()
        self.rootfile.Close()

