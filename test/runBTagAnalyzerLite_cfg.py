
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
options.register('runSubJets', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run subjets"
)
options.register('processStdAK4Jets', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Process standard AK4 jets"
)
options.register('fatJetPtMin', 150.0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum pT for fat jets (default is 150 GeV)"
)
options.register('useTopProjections', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use top projections"
)
options.register('miniAOD', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Running on miniAOD"
)
options.register('useExplicitJTA', False,
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
options.register('useLegacyTaggers', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use legacy taggers"
)
options.register('runIVF', False, # needs to be set to True when running over 7XY AOD where 2<=X<4
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run IVF"
)

## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', 100)

options.parseArguments()

print "Running on data: %s"%('True' if options.runOnData else 'False')
print "Running on MiniAOD: %s"%('True' if options.miniAOD else 'False')
print "Using PFchs: %s"%('True' if options.usePFchs else 'False')

## Jet energy corrections
jetCorrectionsAK4 = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
jetCorrectionsAK8 = ('AK8PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

if not options.usePFchs:
    jetCorrectionsAK4 = ('AK4PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
    jetCorrectionsAK8 = ('AK8PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

if options.runOnData:
    jetCorrectionsAK4[1].append('L2L3Residual')
    jetCorrectionsAK8[1].append('L2L3Residual')

## b-tag infos
bTagInfosLegacy = [
    'impactParameterTagInfos'
   ,'secondaryVertexTagInfos'
   ,'inclusiveSecondaryVertexFinderTagInfos'
   ,'softPFMuonsTagInfos'
   ,'softPFElectronsTagInfos'
]

bTagInfos = [
    'pfImpactParameterTagInfos'
   ,'pfSecondaryVertexTagInfos'
   ,'pfInclusiveSecondaryVertexFinderTagInfos'
   ,'softPFMuonsTagInfos'
   ,'softPFElectronsTagInfos'
]
## b-tag discriminators
bTagDiscriminatorsLegacy = [
     'jetBProbabilityBJetTags'
    ,'jetProbabilityBJetTags'
    ,'positiveOnlyJetBProbabilityBJetTags'
    ,'positiveOnlyJetProbabilityBJetTags'
    ,'negativeOnlyJetBProbabilityBJetTags'
    ,'negativeOnlyJetProbabilityBJetTags'
    ,'trackCountingHighPurBJetTags'
    ,'trackCountingHighEffBJetTags'
    ,'negativeTrackCountingHighEffBJetTags'
    ,'negativeTrackCountingHighPurBJetTags'
    ,'simpleSecondaryVertexHighEffBJetTags'
    ,'simpleSecondaryVertexHighPurBJetTags'
    ,'negativeSimpleSecondaryVertexHighEffBJetTags'
    ,'negativeSimpleSecondaryVertexHighPurBJetTags'
    ,'combinedSecondaryVertexV2BJetTags'
    ,'positiveCombinedSecondaryVertexV2BJetTags'
    ,'negativeCombinedSecondaryVertexV2BJetTags'
    ,'combinedInclusiveSecondaryVertexV2BJetTags'
    ,'positiveCombinedInclusiveSecondaryVertexV2BJetTags'
    ,'negativeCombinedInclusiveSecondaryVertexV2BJetTags'
    ,'softPFMuonBJetTags'
    ,'positiveSoftPFMuonBJetTags'
    ,'negativeSoftPFMuonBJetTags'
    ,'softPFElectronBJetTags'
    ,'positiveSoftPFElectronBJetTags'
    ,'negativeSoftPFElectronBJetTags'
]

bTagDiscriminators = [
     'pfJetBProbabilityBJetTags'
    ,'pfJetProbabilityBJetTags'
    ,'pfPositiveOnlyJetBProbabilityBJetTags'
    ,'pfPositiveOnlyJetProbabilityBJetTags'
    ,'pfNegativeOnlyJetBProbabilityBJetTags'
    ,'pfNegativeOnlyJetProbabilityBJetTags'
    ,'pfTrackCountingHighPurBJetTags'
    ,'pfTrackCountingHighEffBJetTags'
    ,'pfNegativeTrackCountingHighPurBJetTags'
    ,'pfNegativeTrackCountingHighEffBJetTags'
    ,'pfSimpleSecondaryVertexHighEffBJetTags'
    ,'pfSimpleSecondaryVertexHighPurBJetTags'
    ,'pfNegativeSimpleSecondaryVertexHighEffBJetTags'
    ,'pfNegativeSimpleSecondaryVertexHighPurBJetTags'
    ,'pfCombinedSecondaryVertexV2BJetTags'
    ,'pfPositiveCombinedSecondaryVertexV2BJetTags'
    ,'pfNegativeCombinedSecondaryVertexV2BJetTags'
    ,'pfCombinedInclusiveSecondaryVertexV2BJetTags'
    ,'pfPositiveCombinedInclusiveSecondaryVertexV2BJetTags'
    ,'pfNegativeCombinedInclusiveSecondaryVertexV2BJetTags'
    ,'softPFMuonBJetTags'
    ,'positiveSoftPFMuonBJetTags'
    ,'negativeSoftPFMuonBJetTags'
    ,'softPFElectronBJetTags'
    ,'positiveSoftPFElectronBJetTags'
    ,'negativeSoftPFElectronBJetTags'
]

## Legacy taggers not supported with MiniAOD
if options.miniAOD and options.useLegacyTaggers:
    print "WARNING: Legacy taggers not supported with MiniAOD"
    options.useLegacyTaggers = False

## If using legacy taggers
if options.useLegacyTaggers:
   bTagInfos = bTagInfosLegacy
   bTagDiscriminators = bTagDiscriminatorsLegacy

## Clustering algorithm label
algoLabel = 'CA'
if options.jetAlgo == 'AntiKt':
    algoLabel = 'AK'

## Postfix
postfix = "PFlow"
## Various collection names
genParticles = 'genParticles'
jetSource = 'pfJetsPFBRECO'+postfix
genJetCollection = 'ak4GenJetsNoNu'+postfix
pfCandidates = 'particleFlow'
pvSource = 'offlinePrimaryVertices'
svSource = 'inclusiveCandidateSecondaryVertices'
muSource = 'muons'
elSource = 'gedGsfElectrons'
patMuons = 'selectedPatMuons'
## If running on miniAOD
if options.miniAOD:
    genParticles = 'prunedGenParticles'
    jetSource = 'ak4PFJets'
    genJetCollection = 'ak4GenJetsNoNu'
    pfCandidates = 'packedPFCandidates'
    pvSource = 'offlineSlimmedPrimaryVertices'
    svSource = 'slimmedSecondaryVertices'
    muSource = 'slimmedMuons'
    elSource = 'slimmedElectrons'
    patMuons = 'slimmedMuons'

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
        #'/store/user/cvernier/boostedGen/Rad_HHto4b_M800_13TeV/AODSIM/150303_223019/0000/step4_1.root'
        # /RelValProdTTbar_13/CMSSW_7_4_0_pre7-MCRUN2_74_V7-v1/AODSIM
        '/store/relval/CMSSW_7_4_0_pre7/RelValProdTTbar_13/AODSIM/MCRUN2_74_V7-v1/00000/22E552FD-23B7-E411-B680-002618943911.root'
    )
)

if options.miniAOD:
    process.source.fileNames = [
        # /RelValTTbar_13/CMSSW_7_4_0_pre7-MCRUN2_74_V7-v1/MINIAODSIM
        '/store/relval/CMSSW_7_4_0_pre7/RelValTTbar_13/MINIAODSIM/MCRUN2_74_V7-v1/00000/B62A3865-39B7-E411-B76A-002618943880.root'
    ]

if options.runOnData:
    process.source.fileNames = [
        # /SingleMu/CMSSW_7_4_0_pre7-GR_R_74_V8A_RelVal_mu2012D-v1/RECO
        '/store/relval/CMSSW_7_4_0_pre7/SingleMu/RECO/GR_R_74_V8A_RelVal_mu2012D-v1/00000/004E151D-D8B6-E411-A889-0025905B859E.root'
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
process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary),
    allowUnscheduled = cms.untracked.bool(True)
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_' + ('data' if options.runOnData else 'mc'))

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")

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
if not options.miniAOD:
    ## PAT Configuration
    jetAlgo="AK4"

    from PhysicsTools.PatAlgos.tools.pfTools import *
    usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=not options.runOnData, postfix=postfix,
	      jetCorrections=jetCorrectionsAK4, pvCollection=cms.InputTag(pvSource))

    ## Top projections in PF2PAT
    getattr(process,"pfPileUpJME"+postfix).checkClosestZVertex = False
    getattr(process,"pfNoPileUpJME"+postfix).enable = options.usePFchs
    getattr(process,"pfNoMuonJMEPFBRECO"+postfix).enable = options.useTopProjections
    getattr(process,"pfNoElectronJMEPFBRECO"+postfix).enable = options.useTopProjections
else:
    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
    from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
    ## Select isolated collections
    process.selectedMuons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedMuons"), cut = cms.string('''abs(eta)<2.5 && pt>10. &&
       (pfIsolationR04().sumChargedHadronPt+
	max(0.,pfIsolationR04().sumNeutralHadronEt+
	pfIsolationR04().sumPhotonEt-
	0.50*pfIsolationR04().sumPUPt))/pt < 0.20 && 
	(isPFMuon && (isGlobalMuon || isTrackerMuon) )'''))
    process.selectedElectrons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedElectrons"), cut = cms.string('''abs(eta)<2.5 && pt>20. &&
	gsfTrack.isAvailable() &&
	gsfTrack.hitPattern().numberOfLostHits(\'MISSING_INNER_HITS\') < 2 &&
	(pfIsolationVariables().sumChargedHadronPt+
	max(0.,pfIsolationVariables().sumNeutralHadronEt+
	pfIsolationVariables().sumPhotonEt-
	0.5*pfIsolationVariables().sumPUPt))/pt < 0.15'''))

    ## Do projections
    process.pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))
    process.pfNoMuonCHS =  cms.EDProducer("CandPtrProjector", src = cms.InputTag("pfCHS"), veto = cms.InputTag("selectedMuons"))
    process.pfNoElectronsCHS = cms.EDProducer("CandPtrProjector", src = cms.InputTag("pfNoMuonCHS"), veto = cms.InputTag("selectedElectrons"))

    process.pfNoMuon =  cms.EDProducer("CandPtrProjector", src = cms.InputTag("packedPFCandidates"), veto = cms.InputTag("selectedMuons"))
    process.pfNoElectrons = cms.EDProducer("CandPtrProjector", src = cms.InputTag("pfNoMuon"), veto = cms.InputTag("selectedElectrons"))

    process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))
    process.ak4GenJetsNoNu = ak4GenJets.clone(src = 'packedGenParticlesForJetsNoNu')

    if options.usePFchs:
        if options.useTopProjections:
            process.ak4PFJets = ak4PFJets.clone(src = 'pfNoElectronsCHS', doAreaFastjet = True)
        else:
            process.ak4PFJets = ak4PFJets.clone(src = 'pfCHS', doAreaFastjet = True)
    else:
        if options.useTopProjections:
            process.ak4PFJets = ak4PFJets.clone(src = 'pfNoElectrons', doAreaFastjet = True)
        else:
            process.ak4PFJets = ak4PFJets.clone(src = 'packedPFCandidates', doAreaFastjet = True)

## Load standard PAT objects (here we only need PAT muons but the framework will figure out what it needs to run using the unscheduled mode)
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")

from PhysicsTools.PatAlgos.tools.jetTools import *
## Switch the default jet collection (done in order to use the above specified b-tag infos and discriminators)
switchJetCollection(
    process,
    jetSource = cms.InputTag(jetSource),
    pfCandidates = cms.InputTag(pfCandidates),
    pvSource = cms.InputTag(pvSource),
    svSource = cms.InputTag(svSource),
    muSource = cms.InputTag(muSource),
    elSource = cms.InputTag(elSource),
    btagInfos = bTagInfos,
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = jetCorrectionsAK4,
    genJetCollection = cms.InputTag(genJetCollection),
    genParticles = cms.InputTag(genParticles),
    explicitJTA = options.useExplicitJTA,
    runIVF = options.runIVF,
    postfix = postfix
)

#-------------------------------------

#-------------------------------------
## Fat jets (Gen and Reco)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.genJetsNoNu = ca4GenJets.clone(
    jetAlgorithm = cms.string(options.jetAlgo),
    rParam = cms.double(options.jetRadius),
    src = (cms.InputTag("packedGenParticlesForJetsNoNu") if options.miniAOD else cms.InputTag("genParticlesForJetsNoNu"+postfix))
)
from RecoJets.JetProducers.ca4PFJets_cfi import ca4PFJets
process.PFJetsCHS = ca4PFJets.clone(
    jetAlgorithm = cms.string(options.jetAlgo),
    rParam = cms.double(options.jetRadius),
    src = (getattr(process,"ak4PFJets").src if options.miniAOD else getattr(process,"pfJetsPFBRECO"+postfix).src),
    srcPVs = (getattr(process,"ak4PFJets").srcPVs if options.miniAOD else getattr(process,"pfJetsPFBRECO"+postfix).srcPVs),
    doAreaFastjet = cms.bool(True),
    jetPtMin = cms.double(options.fatJetPtMin)
)
## Pruned fat jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.genJetsNoNuPruned = ca4GenJets.clone(
    SubJetParameters,
    jetAlgorithm = cms.string(options.jetAlgo),
    rParam = cms.double(options.jetRadius),
    src = (cms.InputTag("packedGenParticlesForJetsNoNu") if options.miniAOD else cms.InputTag("genParticlesForJetsNoNu"+postfix)),
    usePruning = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.PFJetsCHSPruned = ak5PFJetsPruned.clone(
    jetAlgorithm = cms.string(options.jetAlgo),
    rParam = cms.double(options.jetRadius),
    src = (getattr(process,"ak4PFJets").src if options.miniAOD else getattr(process,"pfJetsPFBRECO"+postfix).src),
    srcPVs = (getattr(process,"ak4PFJets").srcPVs if options.miniAOD else getattr(process,"pfJetsPFBRECO"+postfix).srcPVs),
    doAreaFastjet = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(options.fatJetPtMin)
)

if options.runSubJets:
    ## PATify the above jets
    addJetCollection(
        process,
        labelName='PFCHS',
        jetSource=cms.InputTag('PFJetsCHS'),
        algo=algoLabel,           # needed for jet flavor clustering
        rParam=options.jetRadius, # needed for jet flavor clustering
        pfCandidates = cms.InputTag(pfCandidates),
        pvSource = cms.InputTag(pvSource),
        svSource = cms.InputTag(svSource),
        muSource = cms.InputTag(muSource),
        elSource = cms.InputTag(elSource),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        jetCorrections = jetCorrectionsAK8,
        genJetCollection = cms.InputTag('genJetsNoNu'),
        genParticles = cms.InputTag(genParticles),
        explicitJTA = options.useExplicitJTA,
        runIVF = options.runIVF,
        postfix = postfix
    )
    addJetCollection(
        process,
        labelName='PrunedPFCHS',
        jetSource=cms.InputTag('PFJetsCHSPruned'),
        algo=algoLabel,
        btagInfos=['None'],
        btagDiscriminators=['None'],
        jetCorrections=jetCorrectionsAK8,
        genJetCollection = cms.InputTag('genJetsNoNu'),
        genParticles = cms.InputTag(genParticles),
        getJetMCFlavour = False, # jet flavor disabled
        postfix = postfix
    )
    addJetCollection(
        process,
        labelName='PrunedSubjetsPFCHS',
        jetSource=cms.InputTag('PFJetsCHSPruned','SubJets'),
        algo=algoLabel,           # needed for subjet flavor clustering
        rParam=options.jetRadius, # needed for subjet flavor clustering
        pfCandidates = cms.InputTag(pfCandidates),
        pvSource = cms.InputTag(pvSource),
        svSource = cms.InputTag(svSource),
        muSource = cms.InputTag(muSource),
        elSource = cms.InputTag(elSource),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        jetCorrections = jetCorrectionsAK4,
        genJetCollection = cms.InputTag('genJetsNoNuPruned','SubJets'),
        genParticles = cms.InputTag(genParticles),
        explicitJTA = True,  # needed for subjet b tagging
        svClustering = True, # needed for subjet b tagging
        fatJets = cms.InputTag('PFJetsCHS'),              # needed for subjet flavor clustering
        groomedFatJets = cms.InputTag('PFJetsCHSPruned'), # needed for subjet flavor clustering
        runIVF = options.runIVF,
        postfix = postfix
    )

    ## Establish references between PATified fat jets and subjets using the BoostedJetMerger
    process.selectedPatJetsPrunedPFCHSPacked = cms.EDProducer("BoostedJetMerger",
        jetSrc=cms.InputTag("selectedPatJetsPrunedPFCHS"+postfix),
        subjetSrc=cms.InputTag("selectedPatJetsPrunedSubjetsPFCHS"+postfix)
    )

    #-------------------------------------
    ## N-subjettiness
    from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

    process.Njettiness = Njettiness.clone(
        src = cms.InputTag("PFJetsCHS"),
        cone = cms.double(options.jetRadius)
    )

    getattr(process,'patJetsPFCHS'+postfix).userData.userFloats.src += ['Njettiness:tau1','Njettiness:tau2','Njettiness:tau3']

    #-------------------------------------
    ## Grooming ValueMaps
    process.Pruned = cms.EDProducer("RecoJetDeltaRValueMapProducer",
        src = cms.InputTag("PFJetsCHS"),
        matched = cms.InputTag("PFJetsCHSPruned"),
        distMax = cms.double(options.jetRadius),
        values = cms.vstring('mass','pt'),
        valueLabels = cms.vstring('Mass','Pt'),
        lazyParser = cms.bool(True)
    )

    getattr(process,'patJetsPFCHS'+postfix).userData.userFloats.src += ['Pruned:Mass','Pruned:Pt']

#-------------------------------------

#-------------------------------------
if options.runOnData and options.runSubJets:
    ## Remove MC matching when running over data
    removeMCMatching( process, ['All'] )

#-------------------------------------
## Add GenParticlePruner for boosted b-tagging studies
process.prunedGenParticlesBoost = cms.EDProducer('GenParticlePruner',
    src = cms.InputTag(genParticles),
    select = cms.vstring(
        "drop  *  ", #by default
        "keep ( status = 3 || (status>=21 && status<=29) )", #keep hard process particles
        "keep abs(pdgId) = 13 || abs(pdgId) = 15" #keep muons and taus
    )
)

#-------------------------------------

#-------------------------------------
## Filter for removing scraping events
process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

## Filter for good primary vertex
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag(pvSource),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
)

#-------------------------------------
from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag(pvSource))

#-------------------------------------
## Add full JetFlavourInfo and TagInfos to PAT jets
for m in ['patJets'+postfix, 'patJetsPFCHS'+postfix, 'patJetsPrunedSubjetsPFCHS'+postfix]:
    if hasattr(process,m) and getattr( getattr(process,m), 'addBTagInfo' ):
        print "Switching 'addTagInfos' for " + m + " to 'True'"
        setattr( getattr(process,m), 'addTagInfos', cms.bool(True) )
    if hasattr(process,m):
        print "Switching 'addJetFlavourInfo' for " + m + " to 'True'"
        setattr( getattr(process,m), 'addJetFlavourInfo', cms.bool(True) )

#-------------------------------------
## Adapt fat jet b tagging
if options.runSubJets:
    if options.useLegacyTaggers:
        # Set the cone size for the jet-track association to the jet radius
        getattr(process,'jetTracksAssociatorAtVertexPFCHS'+postfix).coneSize = cms.double(options.jetRadius) # default is 0.4
        getattr(process,'secondaryVertexTagInfosPFCHS'+postfix).trackSelection.jetDeltaRMax = cms.double(options.jetRadius)   # default is 0.3
        getattr(process,'secondaryVertexTagInfosPFCHS'+postfix).vertexCuts.maxDeltaRToJetAxis = cms.double(options.jetRadius) # default is 0.4
        # Set the jet-SV dR to the jet radius
        getattr(process,'inclusiveSecondaryVertexFinderTagInfosPFCHS'+postfix).vertexCuts.maxDeltaRToJetAxis = cms.double(options.jetRadius) # default is 0.4
        getattr(process,'inclusiveSecondaryVertexFinderTagInfosPFCHS'+postfix).extSVDeltaRToJet = cms.double(options.jetRadius) # default is 0.3
        # Set the JP track dR cut to the jet radius
        process.jetProbabilityComputerFat = process.jetProbabilityComputer.clone( deltaR = cms.double(options.jetRadius) ) # default is 0.3
        getattr(process,'jetProbabilityBJetTagsPFCHS'+postfix).jetTagComputer = cms.string('jetProbabilityComputerFat')
        # Set the JBP track dR cut to the jet radius
        process.jetBProbabilityComputerFat = process.jetBProbabilityComputer.clone( deltaR = cms.double(options.jetRadius) ) # default is 0.4
        getattr(process,'jetBProbabilityBJetTagsPFCHS'+postfix).jetTagComputer = cms.string('jetBProbabilityComputerFat')
        # Set the CSV track dR cut to the jet radius
        process.combinedSecondaryVertexComputerFat = process.combinedSecondaryVertexComputer.clone()
        process.combinedSecondaryVertexComputerFat.trackSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
        process.combinedSecondaryVertexComputerFat.trackPseudoSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
        getattr(process,'combinedSecondaryVertexV2BJetTagsPFCHS'+postfix).jetTagComputer = cms.string('combinedSecondaryVertexComputerFat')
        # Set the CSVv2 track dR cut to the jet radius
        process.combinedSecondaryVertexV2ComputerFat = process.combinedSecondaryVertexV2Computer.clone()
        process.combinedSecondaryVertexV2ComputerFat.trackSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
        process.combinedSecondaryVertexV2ComputerFat.trackPseudoSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
        getattr(process,'combinedInclusiveSecondaryVertexV2BJetTagsPFCHS'+postfix).jetTagComputer = cms.string('combinedSecondaryVertexV2ComputerFat')
    else:
        # Set the cone size for the jet-track association to the jet radius
        getattr(process,'pfImpactParameterTagInfosPFCHS'+postfix).maxDeltaR = cms.double(options.jetRadius) # default is 0.4
        getattr(process,'pfSecondaryVertexTagInfosPFCHS'+postfix).trackSelection.jetDeltaRMax = cms.double(options.jetRadius)   # default is 0.3
        getattr(process,'pfSecondaryVertexTagInfosPFCHS'+postfix).vertexCuts.maxDeltaRToJetAxis = cms.double(options.jetRadius) # default is 0.4
        # Set the jet-SV dR to the jet radius
        getattr(process,'pfInclusiveSecondaryVertexFinderTagInfosPFCHS'+postfix).vertexCuts.maxDeltaRToJetAxis = cms.double(options.jetRadius) # default is 0.4
        getattr(process,'pfInclusiveSecondaryVertexFinderTagInfosPFCHS'+postfix).extSVDeltaRToJet = cms.double(options.jetRadius) # default is 0.3
        # Set the JP track dR cut to the jet radius
        process.candidateJetProbabilityComputerFat = process.candidateJetProbabilityComputer.clone( deltaR = cms.double(options.jetRadius) ) # default is 0.3
        getattr(process,'pfJetProbabilityBJetTagsPFCHS'+postfix).jetTagComputer = cms.string('candidateJetProbabilityComputerFat')
        # Set the JBP track dR cut to the jet radius
        process.candidateJetBProbabilityComputerFat = process.candidateJetBProbabilityComputer.clone( deltaR = cms.double(options.jetRadius) ) # default is 0.4
        getattr(process,'pfJetBProbabilityBJetTagsPFCHS'+postfix).jetTagComputer = cms.string('candidateJetBProbabilityComputerFat')
        # Set the CSV track dR cut to the jet radius
        process.candidateCombinedSecondaryVertexComputerFat = process.candidateCombinedSecondaryVertexComputer.clone()
        process.candidateCombinedSecondaryVertexComputerFat.trackSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
        process.candidateCombinedSecondaryVertexComputerFat.trackPseudoSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
        getattr(process,'pfCombinedSecondaryVertexV2BJetTagsPFCHS'+postfix).jetTagComputer = cms.string('candidateCombinedSecondaryVertexComputerFat')
        # Set the CSVv2 track dR cut to the jet radius
        process.candidateCombinedSecondaryVertexV2ComputerFat = process.candidateCombinedSecondaryVertexV2Computer.clone()
        process.candidateCombinedSecondaryVertexV2ComputerFat.trackSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
        process.candidateCombinedSecondaryVertexV2ComputerFat.trackPseudoSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
        getattr(process,'pfCombinedInclusiveSecondaryVertexV2BJetTagsPFCHS'+postfix).jetTagComputer = cms.string('candidateCombinedSecondaryVertexV2ComputerFat')

#-------------------------------------
from RecoBTag.BTagAnalyzerLite.bTagAnalyzerLite_cff import *
process.btagana = bTagAnalyzerLite.clone()
if options.useLegacyTaggers:
    process.btagana = bTagAnalyzerLiteLegacy.clone()
process.btagana.produceJetTrackTree    = True ## True if you want to keep info for tracks associated to jets
process.btagana.produceJetPFLeptonTree = True ## True if you want to keep PF lepton info
process.btagana.storeMuonInfo          = False ## True if you want to keep muon info
process.btagana.storeTagVariables      = False ## True if you want to keep TagInfo TaggingVariables
process.btagana.storeCSVTagVariables   = True  ## True if you want to keep CSV TaggingVariables
process.btagana.primaryVertexColl      = cms.InputTag(pvSource)
process.btagana.Jets                   = cms.InputTag('selectedPatJets'+postfix)
process.btagana.muonCollectionName     = cms.InputTag(patMuons)
process.btagana.triggerTable           = cms.InputTag('TriggerResults::HLT') # Data and MC

if options.runSubJets:
    process.btaganaSubJets = process.btagana.clone(
        storeEventInfo      = cms.bool(not options.processStdAK4Jets),
        allowJetSkipping    = cms.bool(False),
        Jets                = cms.InputTag('selectedPatJetsPrunedPFCHSPacked','SubJets'),
        FatJets             = cms.InputTag('selectedPatJetsPFCHS'+postfix),
        GroomedFatJets      = cms.InputTag('selectedPatJetsPrunedPFCHSPacked'),
        runSubJets          = options.runSubJets,
        svComputerFatJets   = cms.string('combinedSecondaryVertexV2ComputerFat' if options.useLegacyTaggers else 'candidateCombinedSecondaryVertexV2ComputerFat')
    )

#---------------------------------------

#---------------------------------------
## Event counter
from MyAnalysis.EventCounter.eventcounter_cfi import eventCounter
process.allEvents = eventCounter.clone()
process.selectedEvents = eventCounter.clone()
#---------------------------------------

#---------------------------------------
## Define event filter sequence
process.filtSeq = cms.Sequence(
    #process.noscraping
    process.primaryVertexFilter
)


## Define analyzer sequence
process.analyzerSeq = cms.Sequence( )
if options.processStdAK4Jets:
    process.analyzerSeq += process.btagana
if options.runSubJets:
    process.analyzerSeq += process.btaganaSubJets
#---------------------------------------

process.p = cms.Path(
    process.allEvents
    * process.filtSeq
    * process.selectedEvents
    * process.analyzerSeq
)

# Delete predefined output module (needed for running with CRAB)
del process.out

#open('pydump.py','w').write(process.dumpPython())
