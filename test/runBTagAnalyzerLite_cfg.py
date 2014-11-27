
import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
import copy

###############################
####### Parameters ############
###############################

options = VarParsing ('python')

options.register('runOnData', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('outFilename', 'JetTree',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
)
options.register('reportEvery', 1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=1)"
)
options.register('wantSummary', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Print out trigger and timing summary"
)
options.register('usePFchs', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use PFchs"
)
options.register('mcGlobalTag', 'START53_V27',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "MC global tag"
)
options.register('dataGlobalTag', 'FT53_V21A_AN6',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Data global tag"
)
options.register('runSubJets', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run subjets"
)
options.register('processStdAK5Jets', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Process standard AK5 jets"
)
options.register('fatJetPtMin', 150.0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum pT for fat jets (default is 150 GeV)"
)
options.register('useExplicitJTA', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use explicit jet-track association"
)
options.register('jetAlgo', 'AntiKt',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Jet clustering algorithms (default is AntiKt)"
)
options.register('jetRadius', 0.8,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Distance parameter R for jet clustering (default is 0.8)"
)
options.register('useSVClustering', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use SV clustering"
)
options.register('useSVMomentum', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use SV momentum"
)

## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', 100)

options.parseArguments()

print "Running on data: %s"%('True' if options.runOnData else 'False')
print "Using PFchs: %s"%('True' if options.usePFchs else 'False')

## Global tag
globalTag = options.mcGlobalTag
if options.runOnData:
    globalTag = options.dataGlobalTag

## Jet energy corrections
jetCorrectionsAK5 = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
jetCorrectionsAK7 = ('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])

if not options.usePFchs:
    jetCorrectionsAK5 = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute'])
    jetCorrectionsAK7 = ('AK7PF', ['L1FastJet', 'L2Relative', 'L3Absolute'])

if options.runOnData:
    jetCorrectionsAK5[1].append('L2L3Residual')
    jetCorrectionsAK7[1].append('L2L3Residual')

## b-tag infos
bTagInfos = ['impactParameterTagInfos','secondaryVertexTagInfos','inclusiveSecondaryVertexFinderTagInfos',
    'softPFMuonsTagInfos','softPFElectronsTagInfos'
]
## b-tag discriminators
bTagDiscriminators = ['jetBProbabilityBJetTags','jetProbabilityBJetTags','trackCountingHighPurBJetTags','trackCountingHighEffBJetTags',
    #'negativeOnlyJetBProbabilityJetTags','negativeOnlyJetProbabilityJetTags','negativeTrackCountingHighEffJetTags',
    #'negativeTrackCountingHighPurJetTags','positiveOnlyJetBProbabilityJetTags','positiveOnlyJetProbabilityJetTags',
    'simpleSecondaryVertexHighEffBJetTags','simpleSecondaryVertexHighPurBJetTags',
    #'simpleSecondaryVertexNegativeHighEffBJetTags','simpleSecondaryVertexNegativeHighPurBJetTags',
    'combinedSecondaryVertexBJetTags',
    #'combinedSecondaryVertexPositiveBJetTags','combinedSecondaryVertexNegativeBJetTags',
    'softPFMuonBJetTags',
    #'positiveSoftPFMuonBJetTags','negativeSoftPFMuonBJetTags',
    'softPFElectronBJetTags',
    #'positiveSoftPFElectronBJetTags','negativeSoftPFElectronBJetTags',
    #'simpleInclusiveSecondaryVertexHighEffBJetTags','simpleInclusiveSecondaryVertexHighPurBJetTags','doubleSecondaryVertexHighEffBJetTags',
    'combinedInclusiveSecondaryVertexV2BJetTags'

]
bTagDiscriminatorsSubJets = copy.deepcopy(bTagDiscriminators)
if 'doubleSecondaryVertexHighEffBJetTags' in bTagDiscriminators:
    bTagDiscriminatorsSubJets.remove('doubleSecondaryVertexHighEffBJetTags')

## Clustering algorithm label
algoLabel = 'CA'
if options.jetAlgo == 'AntiKt':
    algoLabel = 'AK'

process = cms.Process("BTagAna")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.MessageLogger.cerr.default.limit = 10

## Input files
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        # /QCD_Pt-80to120_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v2/GEN-SIM-RECODEBUG
        #'/store/mc/Summer12_DR53X/QCD_Pt-80to120_TuneZ2star_8TeV_pythia6/GEN-SIM-RECODEBUG/PU_S10_START53_V7A-v2/0000/041C8D66-05F4-E111-B16E-003048D43656.root'
        # /QCD_Pt-300to470_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/GEN-SIM-RECODEBUG
        #'/store/mc/Summer12_DR53X/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6/GEN-SIM-RECODEBUG/PU_S10_START53_V7A-v1/0000/0A4671D2-02F4-E111-9CD8-003048C69310.root'
        # /QCD_Pt-170to300_MuEnrichedPt5_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
        #'/store/mc/Summer12_DR53X/QCD_Pt-170to300_MuEnrichedPt5_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/FEB355DA-8EE7-E111-BA8A-001EC9D83165.root'
        # /QCD_Pt-470to600_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM
        #'/store/mc/Summer12_DR53X/QCD_Pt-470to600_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v2/00000/FADB0913-1708-E211-BBB1-00261894383C.root'
        # /TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM
        '/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7C-v1/00000/FE7C71D8-DB25-E211-A93B-0025901D4C74.root'
    )
)

if options.runOnData:
    process.source.fileNames = [
        # /Jet/Run2012A-22Jan2013-v1/AOD
        #'/store/data/Run2012A/Jet/AOD/22Jan2013-v1/20000/30B21345-4172-E211-9EF3-00304867BEC0.root'
        # /BTag/Run2012A-22Jan2013-v1/AOD
        '/store/data/Run2012A/BTag/AOD/22Jan2013-v1/30000/E6959DA6-7081-E211-ABD3-002590596498.root'
    ]

if options.runOnData :
    if options.runSubJets :
        options.outFilename += '_data_subjets.root'
    else :
        options.outFilename += '_data.root'
else :
    if options.runSubJets :
        options.outFilename += '_mc_subjets.root'
    else :
        options.outFilename += '_mc.root'

## Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string(options.outFilename)
)

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(options.wantSummary) )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = globalTag + '::All'

##############################################
# Get calibrations for the CSVV2 tagger
##############################################
process.load('CondCore.DBCommon.CondDBSetup_cfi')
process.BTauMVAJetTagComputerRecord = cms.ESSource('PoolDBESSource',
    process.CondDBSetup,
    timetype = cms.string('runnumber'),
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('BTauGenericMVAJetTagComputerRcd'),
        tag = cms.string('MVAComputerContainer_53X_JetTags_v2')
    )),
    connect = cms.string('frontier://FrontierProd/CMS_COND_PAT_000'),
    BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService')
)
process.es_prefer_BTauMVAJetTagComputerRecord = cms.ESPrefer('PoolDBESSource','BTauMVAJetTagComputerRecord')


process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")

#-------------------------------------
## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.outFilename),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    # save PAT Layer 1 output; you need a '*' to
    # unpack the list of commands 'patEventContent'
    outputCommands = cms.untracked.vstring('drop *', *patEventContent)
)

#-------------------------------------
## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

postfix = "PFlow"
jetAlgo="AK5"

from PhysicsTools.PatAlgos.tools.pfTools import *
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=not options.runOnData, postfix=postfix,
          jetCorrections=jetCorrectionsAK5, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'))

## Top projections in PF2PAT
getattr(process,"pfPileUp"+postfix).checkClosestZVertex = False
getattr(process,"pfNoPileUp"+postfix).enable = options.usePFchs
getattr(process,"pfNoMuon"+postfix).enable = True
getattr(process,"pfNoElectron"+postfix).enable = True
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True

from PhysicsTools.PatAlgos.tools.coreTools import *
## Remove objects not used from the PAT sequences to speed up processing
removeSpecificPATObjects(process,names=['Electrons', 'Muons', 'Taus'],postfix=postfix)

from PhysicsTools.PatAlgos.tools.jetTools import *
## Switch the default jet collection (done in order to use the above specified b-tag infos and discriminators)
switchJetCollection(process,
    jetCollection=cms.InputTag('pfNoTau'+postfix),
    jetIdLabel='ak',
    rParam = 0.5,
    useLegacyFlavour=False,
    doJTA        = True,
    doBTagging   = True,
    btagInfo     = bTagInfos,
    btagdiscriminators = bTagDiscriminators,
    jetCorrLabel = jetCorrectionsAK5,
    doType1MET   = False,
    genJetCollection = cms.InputTag('ak5GenJetsNoNu'),
    doJetID      = False,
    postfix      = postfix
)

#-------------------------------------

#-------------------------------------
## Fat jets (Gen and Reco)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.genJetsNoNu = ca4GenJets.clone(
    jetAlgorithm = cms.string(options.jetAlgo),
    rParam = cms.double(options.jetRadius),
    src = cms.InputTag("genParticlesForJetsNoNu")
)
from RecoJets.JetProducers.ca4PFJets_cfi import ca4PFJets
process.PFJetsCHS = ca4PFJets.clone(
    jetAlgorithm = cms.string(options.jetAlgo),
    rParam = cms.double(options.jetRadius),
    src = getattr(process,"pfJets"+postfix).src,
    srcPVs = getattr(process,"pfJets"+postfix).srcPVs,
    doAreaFastjet = cms.bool(True),
    jetPtMin = cms.double(options.fatJetPtMin)
)
## Pruned fat jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.genJetsNoNuPruned = ca4GenJets.clone(
    SubJetParameters,
    jetAlgorithm = cms.string(options.jetAlgo),
    rParam = cms.double(options.jetRadius),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    usePruning = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.PFJetsCHSPruned = ak5PFJetsPruned.clone(
    jetAlgorithm = cms.string(options.jetAlgo),
    rParam = cms.double(options.jetRadius),
    src = getattr(process,"pfJets"+postfix).src,
    srcPVs = getattr(process,"pfJets"+postfix).srcPVs,
    doAreaFastjet = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(options.fatJetPtMin)
)

if options.runSubJets:
    ## PATify the above jets
    switchJetCollection(process,
        jetCollection=cms.InputTag('PFJetsCHS'),
        jetIdLabel=algoLabel,
        rParam = options.jetRadius,
        useLegacyFlavour=False,
        doJTA        = True,
        doBTagging   = True,
        btagInfo     = bTagInfos,
        btagdiscriminators = bTagDiscriminators,
        jetCorrLabel = jetCorrectionsAK7,
        doType1MET   = False,
        genJetCollection = cms.InputTag('genJetsNoNu'),
        doJetID      = False
    )
    addJetCollection(
        process,
        jetCollection=cms.InputTag('PFJetsCHSPruned'),
        algoLabel=algoLabel,
        typeLabel='PrunedPFCHS',
        getJetMCFlavour=False,
        doJTA=False,
        doBTagging=False,
        btagInfo=bTagInfos,
        btagdiscriminators=bTagDiscriminators,
        jetCorrLabel=jetCorrectionsAK7,
        doType1MET=False,
        doL1Cleaning=False,
        doL1Counters=False,
        doJetID=False,
        genJetCollection=cms.InputTag('genJetsNoNu')
    )
    addJetCollection(
        process,
        jetCollection=cms.InputTag('PFJetsCHSPruned','SubJets'),
        algoLabel=algoLabel,
        typeLabel='PrunedSubjetsPFCHS',
        rParam = options.jetRadius,
        useLegacyFlavour=False,
        doJTA=True,
        doBTagging=True,
        btagInfo=bTagInfos,
        btagdiscriminators=bTagDiscriminatorsSubJets,
        jetCorrLabel=jetCorrectionsAK5,
        doType1MET=False,
        doL1Cleaning=False,
        doL1Counters=False,
        doJetID=False,
        genJetCollection=cms.InputTag('genJetsNoNuPruned','SubJets')
    )

    ## Establish references between PAT fat jets and PAT subjets using the BoostedJetMerger
    process.selectedPatJetsPrunedPFCHSPacked = cms.EDProducer("BoostedJetMerger",
        jetSrc=cms.InputTag("selectedPatJets"+algoLabel+"PrunedPFCHS"),
        subjetSrc=cms.InputTag("selectedPatJets"+algoLabel+"PrunedSubjetsPFCHS")
    )

    ## New jet flavor still requires some cfg-level adjustments for subjets until it is better integrated into PAT
    ## Adjust the jet flavor for pruned subjets
    setattr(process,'patJetFlavourAssociation'+algoLabel+'PrunedSubjetsPFCHS', process.patJetFlavourAssociation.clone(
        groomedJets = cms.InputTag("PFJetsCHSPruned"),
        subjets = cms.InputTag("PFJetsCHSPruned", "SubJets")
    ))
    getattr(process,'patJets'+algoLabel+'PrunedSubjetsPFCHS').JetFlavourInfoSource = cms.InputTag("patJetFlavourAssociation"+algoLabel+"PrunedSubjetsPFCHS","SubJets")

    #-------------------------------------
    ## N-subjettiness
    from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

    process.Njettiness = Njettiness.clone(
        src = cms.InputTag("PFJetsCHS"),
        cone = cms.double(options.jetRadius)
    )

    process.patJets.userData.userFloats.src += ['Njettiness:tau1','Njettiness:tau2','Njettiness:tau3']

    #-------------------------------------
    ## Grooming ValueMaps
    from RecoJets.JetProducers.ca8PFJetsCHS_groomingValueMaps_cfi import ca8PFJetsCHSPrunedLinks

    process.PFJetsCHSPrunedMass = ca8PFJetsCHSPrunedLinks.clone(
        src = cms.InputTag("PFJetsCHS"),
        matched = cms.InputTag("PFJetsCHSPruned"),
        distMax = cms.double(options.jetRadius),
        value = cms.string('mass')
    )

    process.PFJetsCHSPrunedPt = ca8PFJetsCHSPrunedLinks.clone(
        src = cms.InputTag("PFJetsCHS"),
        matched = cms.InputTag("PFJetsCHSPruned"),
        distMax = cms.double(options.jetRadius),
        value = cms.string('pt')
    )

    process.patJets.userData.userFloats.src += ['PFJetsCHSPrunedMass','PFJetsCHSPrunedPt']

    #-------------------------------------
    if options.useSVClustering:
        setattr(process,'inclusiveSecondaryVertexFinderTagInfos'+algoLabel+'PrunedSubjetsPFCHS', getattr(process,'inclusiveSecondaryVertexFinderTagInfos'+algoLabel+'PrunedSubjetsPFCHS').clone(
            useSVClustering = cms.bool(True),
            useSVMomentum   = cms.bool(options.useSVMomentum), # otherwise using SV flight direction
            jetAlgorithm    = cms.string(options.jetAlgo),
            rParam          = cms.double(options.jetRadius),
            ghostRescaling  = cms.double(1e-18),
            fatJets         = cms.InputTag("PFJetsCHS"),
            groomedFatJets  = cms.InputTag("PFJetsCHSPruned")
        ))

#-------------------------------------
#from PhysicsTools.PatAlgos.tools.coreTools import * # Already imported above
## Remove objects not used from the PAT sequences to speed up processing
if options.runSubJets:
    removeAllPATObjectsBut(process, ['Jets', 'Muons'])

if options.runOnData and options.runSubJets:
    ## Remove MC matching when running over data
    removeMCMatching( process, ['All'] )

#-------------------------------------
## Add GenParticlePruner for boosted b-tagging studies
process.prunedGenParticlesBoost = cms.EDProducer('GenParticlePruner',
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
        "drop  *  ", #by default
        "keep ( status = 3 || (status>=21 && status<=29) )", #keep hard process particles
        "keep abs(pdgId) = 13 || abs(pdgId) = 15" #keep muons and taus
    )
)

#-------------------------------------

#-------------------------------------
## Produce a collection of good primary vertices
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
)

## Filter for removing scraping events
process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

## Filter for good primary vertex
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi")
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
)
#-------------------------------------

#-------------------------------------
## If using explicit jet-track association
if options.useExplicitJTA:
    from RecoJets.JetAssociationProducers.ak5JTA_cff import ak5JetTracksAssociatorExplicit
    for m in getattr(process,"patDefaultSequence"+postfix).moduleNames():
        if m.startswith('jetTracksAssociatorAtVertex'):
            print 'Switching ' + m + ' to explicit jet-track association'
            setattr( process, m, ak5JetTracksAssociatorExplicit.clone(jets = getattr(getattr(process,m),'jets')) )
    if options.runSubJets:
        for m in getattr(process,"patDefaultSequence").moduleNames():
            if m.startswith('jetTracksAssociatorAtVertex'):
                print 'Switching ' + m + ' to explicit jet-track association'
                setattr( process, m, ak5JetTracksAssociatorExplicit.clone(jets = getattr(getattr(process,m),'jets')) )

#-------------------------------------
## Add full JetFlavourInfo and TagInfos to PAT jets
for m in ['patJets'+postfix, 'patJets', 'patJets'+algoLabel+'PrunedSubjetsPFCHS']:
    if hasattr(process,m) and getattr( getattr(process,m), 'addBTagInfo' ):
        print "Switching 'addTagInfos' for " + m + " to 'True'"
        setattr( getattr(process,m), 'addTagInfos', cms.bool(True) )
    if hasattr(process,m):
        print "Switching 'addJetFlavourInfo' for " + m + " to 'True'"
        setattr( getattr(process,m), 'addJetFlavourInfo', cms.bool(True) )

#-------------------------------------
## Adapt fat jet b tagging
if options.runSubJets:
    # Set the cone size for the jet-track association to the jet radius
    process.jetTracksAssociatorAtVertex.coneSize = cms.double(options.jetRadius) # default is 0.5
    process.secondaryVertexTagInfosAOD.trackSelection.jetDeltaRMax = cms.double(options.jetRadius)   # default is 0.3
    process.secondaryVertexTagInfosAOD.vertexCuts.maxDeltaRToJetAxis = cms.double(options.jetRadius) # default is 0.5
    # Set the jet-SV dR to the jet radius
    process.inclusiveSecondaryVertexFinderTagInfosAOD.vertexCuts.maxDeltaRToJetAxis = cms.double(options.jetRadius) # default is 0.5
    process.inclusiveSecondaryVertexFinderTagInfosAOD.extSVDeltaRToJet = cms.double(options.jetRadius) # default is 0.3
    # Set the JP track dR cut to the jet radius
    process.jetProbabilityFat = process.jetProbability.clone( deltaR = cms.double(options.jetRadius) ) # default is 0.3
    process.jetProbabilityBJetTagsAOD.jetTagComputer = cms.string('jetProbabilityFat')
    # Set the JBP track dR cut to the jet radius
    process.jetBProbabilityFat = process.jetBProbability.clone( deltaR = cms.double(options.jetRadius) ) # default is 0.5
    process.jetBProbabilityBJetTagsAOD.jetTagComputer = cms.string('jetBProbabilityFat')
    # Set the CSV track dR cut to the jet radius
    process.combinedSecondaryVertexFat = process.combinedSecondaryVertex.clone()
    process.combinedSecondaryVertexFat.trackSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
    process.combinedSecondaryVertexFat.trackPseudoSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
    process.combinedSecondaryVertexBJetTagsAOD.jetTagComputer = cms.string('combinedSecondaryVertexFat')
    # Set the CSVV2 track dR cut to the jet radius
    process.combinedSecondaryVertexV2Fat = process.combinedSecondaryVertexV2.clone()
    process.combinedSecondaryVertexV2Fat.trackSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
    process.combinedSecondaryVertexV2Fat.trackPseudoSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
    process.combinedInclusiveSecondaryVertexV2BJetTagsAOD.jetTagComputer = cms.string('combinedSecondaryVertexV2Fat')

#-------------------------------------
if not options.runOnData:
    ## JP calibration for 53X MC
    process.GlobalTag.toGet = cms.VPSet(
      cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
           tag = cms.string("TrackProbabilityCalibration_2D_MC53X_v2"),
           connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
      cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
           tag = cms.string("TrackProbabilityCalibration_3D_MC53X_v2"),
           connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
    )
#-------------------------------------

#-------------------------------------
process.load("RecoBTag.BTagAnalyzerLite.BTagAnalyzerLite_cfi")
#------------------
process.btagana.produceJetTrackTree    = False ## True if you want to keep info for tracks associated to jets
process.btagana.produceJetPFLeptonTree = False ## True if you want to keep PF lepton info
process.btagana.storeMuonInfo          = False ## True if you want to keep muon info
process.btagana.storeTagVariables      = False ## True if you want to keep TagInfo TaggingVariables
process.btagana.storeCSVTagVariables   = True  ## True if you want to keep CSV TaggingVariables
process.btagana.primaryVertexColl      = cms.InputTag('goodOfflinePrimaryVertices')
process.btagana.Jets                   = cms.InputTag('selectedPatJets'+postfix)
process.btagana.muonCollectionName     = cms.InputTag('selectedPatMuons')
process.btagana.triggerTable           = cms.InputTag('TriggerResults::HLT') # Data and MC
process.btagana.svComputer             = cms.string('combinedSecondaryVertexV2')
process.btagana.svTagInfos             = cms.string('inclusiveSecondaryVertexFinder')

if options.runSubJets:
    process.btaganaSubJets = process.btagana.clone(
        storeEventInfo      = cms.bool(not options.processStdAK5Jets),
        allowJetSkipping    = cms.bool(False),
        Jets                = cms.InputTag('selectedPatJets'+algoLabel+'PrunedSubjetsPFCHS'),
        FatJets             = cms.InputTag('selectedPatJets'),
        GroomedFatJets      = cms.InputTag('selectedPatJetsPrunedPFCHSPacked'),
        runSubJets          = options.runSubJets,
        svComputer          = cms.string('combinedSecondaryVertexV2'),
        svComputerFatJets   = cms.string('combinedSecondaryVertexV2Fat'),
        svTagInfos          = cms.string('inclusiveSecondaryVertexFinder')
    )

#-------------------------------------

#-------------------------------------
## Optional MET filters:
## https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters
process.load("RecoMET.METFilters.metFilters_cff")
process.trackingFailureFilter.VertexSource = cms.InputTag('goodOfflinePrimaryVertices')
#-------------------------------------

#-------------------------------------
## Filter for HCAL laser events in prompt 2012A+B+C, snippet for "Datasets from the 2013 rereco and Multijet parked":
## https://twiki.cern.ch/twiki/bin/view/CMS/PdmVKnowFeatures#HCAL_laser_events_in_prompt_2012
process.load("EventFilter.HcalRawToDigi.hcallaserFilterFromTriggerResult_cff")
#-------------------------------------

#-------------------------------------
## Event counter
from MyAnalysis.EventCounter.eventcounter_cfi import eventCounter
process.allEvents = eventCounter.clone()
process.selectedEvents = eventCounter.clone()
#-------------------------------------

#-------------------------------------
## Define event filter sequence
process.filtSeq = cms.Sequence(
    process.noscraping
    * process.primaryVertexFilter
    * process.goodOfflinePrimaryVertices
    * process.HBHENoiseFilter
    * process.CSCTightHaloFilter
    * process.EcalDeadCellTriggerPrimitiveFilter
    * process.eeBadScFilter
    * process.ecalLaserCorrFilter
    * process.trackingFailureFilter
    * process.trkPOGFilters
)
if options.runOnData:
    process.filtSeq *= process.hcalfilter
if not options.runOnData:
    process.filtSeq = cms.Sequence( process.prunedGenParticlesBoost * process.filtSeq )

## Define jet sequences
process.genJetSeq = cms.Sequence(
    process.genJetsNoNu
    + process.genJetsNoNuPruned
)
process.jetSeq = cms.Sequence(
    process.PFJetsCHS
    + process.PFJetsCHSPruned
)

if options.runSubJets:
    process.jetSeq *= cms.Sequence(
        process.Njettiness
        + process.PFJetsCHSPrunedMass
        + process.PFJetsCHSPrunedPt
    )

if not options.runOnData:
    process.jetSeq = cms.Sequence( process.genJetSeq + process.jetSeq )


## Define combined PF2PAT + subjet sequence
process.combPF2PATSubJetSeq = cms.Sequence( getattr(process,"patPF2PATSequence"+postfix) )
if options.runSubJets:
    process.combPF2PATSubJetSeq = cms.Sequence(
        getattr(process,"patPF2PATSequence"+postfix)
        * process.jetSeq
        * getattr(process,"patDefaultSequence")
        * process.selectedPatJetsPrunedPFCHSPacked
    )

## Define analyzer sequence
process.analyzerSeq = cms.Sequence( )
if options.processStdAK5Jets:
    process.analyzerSeq += process.btagana
if options.runSubJets:
    process.analyzerSeq += process.btaganaSubJets
#-------------------------------------

#-------------------------------------
from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), postfix=postfix, sequence='patPF2PATSequence')
if options.runSubJets:
    adaptPVs(process, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), postfix='', sequence='jetSeq')
    adaptPVs(process, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), postfix='', sequence='patDefaultSequence')

#-------------------------------------
## Remove tau stuff that really shouldn't be there (probably a bug in PAT)
process.patDefaultSequencePFlow.remove(process.kt6PFJetsForRhoComputationVoronoiPFlow)
for m in getattr(process,"patDefaultSequence"+postfix).moduleNames():
    if m.startswith('hpsPFTau'):
        getattr(process,"patDefaultSequence"+postfix).remove(getattr(process,m))

if options.runSubJets:
    process.patDefaultSequence.remove(process.kt6PFJetsForRhoComputationVoronoi)
    for m in getattr(process,"patDefaultSequence").moduleNames():
        if m.startswith('hpsPFTau'):
            getattr(process,"patDefaultSequence").remove(getattr(process,m))

#-------------------------------------

process.p = cms.Path(
    process.allEvents
    * process.filtSeq
    * process.selectedEvents
    * process.combPF2PATSubJetSeq
    * process.analyzerSeq
)

# Delete predefined output module (needed for running with CRAB)
del process.out

#open('pydump.py','w').write(process.dumpPython())
