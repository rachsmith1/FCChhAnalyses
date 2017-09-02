import os
import copy
import heppy.framework.config as cfg
import sys
import logging
# next 2 lines necessary to deal with reimports from ipython
logging.shutdown()
reload(logging)
logging.basicConfig(level=logging.WARNING)
sys.path.append('/afs/cern.ch/work/h/helsens/public/FCCDicts/')
comp = cfg.Component(
    'example',
     files = ["root://eospublic.cern.ch///eos/fcc/hh/generation/DelphesEvents/fcc_v01/pp_Zprime_10TeV_ttbar/events0.root"]
)

#from heppySampleList_cms import *
from heppySampleList_fcc_v01 import *

selectedComponents = [
			pp_vv_M_5000_10000,
                        pp_vv_M_10000_15000,
                        pp_vv_M_15000_100000,
			pp_jj_M_5000_10000,
			pp_jj_M_10000_15000,
			pp_jj_M_15000_100000,
			pp_tt_M_5000_10000,
                        pp_tt_M_10000_15000,
                        pp_tt_M_15000_100000,
                        pp_Zprime_2TeV_ttbar,
			pp_Zprime_5TeV_ttbar,
                        pp_Zprime_10TeV_ttbar,
                        pp_Zprime_15TeV_ttbar,
                        pp_Zprime_20TeV_ttbar,
                        pp_Zprime_25TeV_ttbar,
                        pp_Zprime_30TeV_ttbar,
                        pp_Zprime_35TeV_ttbar,
                        pp_Zprime_40TeV_ttbar,  
		     ]

pp_vv_M_5000_10000.splitFactor = 100
pp_vv_M_10000_15000.splitFactor = 100
pp_vv_M_15000_100000.splitFactor = 100

pp_jj_M_5000_10000.splitFactor = 100
pp_jj_M_10000_15000.splitFactor = 100
pp_jj_M_15000_100000.splitFactor = 100

pp_tt_M_5000_10000.splitFactor = 100
pp_tt_M_10000_15000.splitFactor = 100
pp_tt_M_15000_100000.splitFactor = 100

pp_Zprime_2TeV_ttbar.splitFactor = 100
pp_Zprime_5TeV_ttbar.splitFactor = 100
pp_Zprime_10TeV_ttbar.splitFactor = 100
pp_Zprime_15TeV_ttbar.splitFactor = 100
pp_Zprime_20TeV_ttbar.splitFactor = 100
pp_Zprime_25TeV_ttbar.splitFactor = 100
pp_Zprime_30TeV_ttbar.splitFactor = 100
pp_Zprime_35TeV_ttbar.splitFactor = 100
pp_Zprime_40TeV_ttbar.splitFactor = 100

#selectedComponents = [comp]

from heppy.analyzers.fcc.Reader import Reader
source = cfg.Analyzer(
    Reader,

    weights = 'mcEventWeights',
    met = 'met',   
 
    electrons = 'electrons',
    muons = 'muons',
    jets = 'jets',

    fatjets = 'fatjets',
    jetsOneSubJettiness = 'jetsOneSubJettiness', 
    jetsTwoSubJettiness = 'jetsTwoSubJettiness', 
    jetsThreeSubJettiness = 'jetsThreeSubJettiness', 
    subjetsTrimmingTagged = 'subjetsTrimmingTagged', 
    subjetsTrimming = 'subjetsTrimming', 
    subjetsPruningTagged = 'subjetsPruningTagged', 
    subjetsPruning = 'subjetsPruning', 
    subjetsSoftDropTagged = 'subjetsSoftDropTagged', 
    subjetsSoftDrop = 'subjetsSoftDrop', 

)


from ROOT import gSystem
gSystem.Load("libdatamodelDict")
from EventStore import EventStore as Events

#############################
##   Reco Level Analysis   ##
#############################

#uncomment the following to go back to normal

from heppy.analyzers.Selector import Selector
# select fatjets above 500 GeV
fatjets_500 = cfg.Analyzer(
    Selector,
    'fatjets_500',
    output = 'fatjets_500',
    input_objects = 'fatjets',
    filter_func = lambda fatjet: fatjet.pt()>500.
)

# select small jets above 30 GeV
jets_30 = cfg.Analyzer(
    Selector,
    'jets_30',
    output = 'jets_30',
    input_objects = 'jets',
    filter_func = lambda jet: jet.pt()>30.
)

# select electrons above 100 GeV
electrons_100 = cfg.Analyzer(
    Selector,
    'electrons_100',
    output = 'electrons_100',
    input_objects = 'electrons',
    filter_func = lambda electron: electron.pt()>100.
)

# select muons above 100 GeV
muons_100 = cfg.Analyzer(
    Selector,
    'muons_100',
    output = 'muons_100',
    input_objects = 'muons',
    filter_func = lambda muon: muon.pt()>100.
)

# produce flat root tree containing jet substructure information
from heppy.FCChhAnalyses.Zprime_tt.TreeProducer import TreeProducer
tree = cfg.Analyzer(
    TreeProducer,
    fatjets = 'fatjets_500',
    jets = 'jets_30',
    electrons = 'electrons_100',
    muons = 'muons_100',
)


# definition of a sequence of analyzers,
# the analyzers will process each event in this order
sequence = cfg.Sequence( [
    source,
    fatjets_500,
    jets_30,
    electrons_100,
    muons_100,
    tree,
    ] )

config = cfg.Config(
    components = selectedComponents,
    sequence = sequence,
    services = [],
    events_class = Events
)

if __name__ == '__main__':
    import sys
    from heppy.framework.looper import Looper

    def next():
        loop.process(loop.iEvent+1)

    loop = Looper( 'looper', config,
                   nEvents=100,
                   nPrint=0,
                   timeReport=True)
    loop.process(6)
    print loop.event
