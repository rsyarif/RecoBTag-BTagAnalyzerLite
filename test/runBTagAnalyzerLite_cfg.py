
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
options.register('mcGlobalTag', 'PHYS14_25_V1',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "MC global tag"
)
options.register('dataGlobalTag', 'GR_R_70_V2',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Data global tag"
)
options.register('runSubJets', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run subjets"
)
options.register('processStdAK4Jets', True,
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
print "Running on miniAOD: %s"%('True' if options.miniAOD else 'False')
print "Using PFchs: %s"%('True' if options.usePFchs else 'False')

## Global tag
globalTag = options.mcGlobalTag
if options.runOnData:
    globalTag = options.dataGlobalTag

## Jet energy corrections
jetCorrectionsAK4 = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
jetCorrectionsAK7 = ('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

if not options.usePFchs:
    jetCorrectionsAK4 = ('AK4PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
    jetCorrectionsAK7 = ('AK7PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

if options.runOnData:
    jetCorrectionsAK4[1].append('L2L3Residual')
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

## Postfix
postfix = "PFlow"
## Various collection names
genParticles = 'genParticles'
jetSource = 'pfJets'+postfix
genJetCollection = 'ak4GenJetsNoNu'+postfix
trackSource = 'generalTracks'
pvSource = 'offlinePrimaryVertices'
svSource = cms.InputTag('inclusiveSecondaryVertices')
muons = 'selectedPatMuons'
## If running on miniAOD
if options.miniAOD:
    genParticles = 'prunedGenParticles'
    jetSource = 'ak4PFJets'
    genJetCollection = 'ak4GenJetsNoNu'
    trackSource = 'unpackedTracksAndVertices'
    pvSource = 'unpackedTracksAndVertices'
    svSource = cms.InputTag('unpackedTracksAndVertices','secondary')
    muons = 'selectedMuons'

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
        # /RelValTTbar_13/CMSSW_7_2_0-PU50ns_PHYS14_25_V1_Phys14-v2/GEN-SIM-RECO
        '/store/relval/CMSSW_7_2_0/RelValTTbar_13/GEN-SIM-RECO/PU50ns_PHYS14_25_V1_Phys14-v2/00000/0E8D62AF-9059-E411-9066-0025905B85D0.root'
    )
)

if options.miniAOD:
    process.source.fileNames = [
        # /RelValTTbar_13/CMSSW_7_2_0-PU50ns_PHYS14_25_V1_Phys14-v2/MINIAODSIM
        '/store/relval/CMSSW_7_2_0/RelValTTbar_13/MINIAODSIM/PU50ns_PHYS14_25_V1_Phys14-v2/00000/206F2CD1-9C59-E411-B789-00261894388A.root'
    ]

if options.runOnData:
    process.source.fileNames = [
        # /JetHT/CMSSW_6_2_1-GR_R_62_V1_dvmc_RelVal_jetHT2012Cdvmc-v1/RECO
        '/store/relval/CMSSW_6_2_1/JetHT/RECO/GR_R_62_V1_dvmc_RelVal_jetHT2012Cdvmc-v1/00000/0030D71C-1E38-E311-BA7A-0025905964BE.root'
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
#process.GlobalTag.globaltag = globalTag + '::All'
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

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
if not options.miniAOD:
    ## PAT Configuration
    jetAlgo="AK4"

    from PhysicsTools.PatAlgos.tools.pfTools import *
    usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=not options.runOnData, postfix=postfix,
	      jetCorrections=jetCorrectionsAK4, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'))

    ## Top projections in PF2PAT
    getattr(process,"pfPileUpJME"+postfix).checkClosestZVertex = False
    getattr(process,"pfNoPileUpJME"+postfix).enable = options.usePFchs
    getattr(process,"pfNoMuonJME"+postfix).enable = options.useTopProjections
    getattr(process,"pfNoElectronJME"+postfix).enable = options.useTopProjections
    getattr(process,"pfNoTau"+postfix).enable = False
    getattr(process,"pfNoJet"+postfix).enable = False
else:
    ## Recreate tracks and PVs for b tagging
    process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
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
    process.pfNoElectronsCHS = cms.EDProducer("CandPtrProjector", src = cms.InputTag("pfNoMuonCHS"), veto =  cms.InputTag("selectedElectrons"))

    process.pfNoMuon =  cms.EDProducer("CandPtrProjector", src = cms.InputTag("packedPFCandidates"), veto = cms.InputTag("selectedMuons"))
    process.pfNoElectrons = cms.EDProducer("CandPtrProjector", src = cms.InputTag("pfNoMuon"), veto =  cms.InputTag("selectedElectrons"))

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

from PhysicsTools.PatAlgos.tools.jetTools import *
## Switch the default jet collection (done in order to use the above specified b-tag infos and discriminators)
switchJetCollection(
    process,
    jetSource = cms.InputTag(jetSource),
    trackSource = cms.InputTag(trackSource),
    pvSource = cms.InputTag(pvSource),
    svSource = svSource,
    btagInfos = bTagInfos,
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = jetCorrectionsAK4,
    genJetCollection = cms.InputTag(genJetCollection),
    explicitJTA=(options.useExplicitJTA and not options.miniAOD), # old-style explicit JTA does not work with miniAOD
    postfix = postfix
)
getattr(process,'patJetPartons'+postfix).particles = cms.InputTag(genParticles)
getattr(process,'patJetPartonMatch'+postfix).matched = cms.InputTag(genParticles)
process.inclusiveVertexFinder.tracks = cms.InputTag(trackSource)
process.trackVertexArbitrator.tracks = cms.InputTag(trackSource)

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
    src = (getattr(process,"ak4PFJets").src if options.miniAOD else getattr(process,"pfJets"+postfix).src),
    srcPVs = (getattr(process,"ak4PFJets").srcPVs if options.miniAOD else getattr(process,"pfJets"+postfix).srcPVs),
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
    src = (getattr(process,"ak4PFJets").src if options.miniAOD else getattr(process,"pfJets"+postfix).src),
    srcPVs = (getattr(process,"ak4PFJets").srcPVs if options.miniAOD else getattr(process,"pfJets"+postfix).srcPVs),
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
        algo=algoLabel,
        rParam=options.jetRadius,
        trackSource = cms.InputTag(trackSource),
        pvSource = cms.InputTag(pvSource),
        svSource = svSource,
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        jetCorrections = jetCorrectionsAK7,
        genJetCollection = cms.InputTag('genJetsNoNu'),
        explicitJTA=(options.useExplicitJTA and not options.miniAOD), # old-style explicit JTA does not work with miniAOD
        postfix = postfix
    )
    getattr(process,'patJetPartonMatchPFCHS'+postfix).matched = cms.InputTag(genParticles)
    addJetCollection(
        process,
        labelName='PrunedPFCHS',
        jetSource=cms.InputTag('PFJetsCHSPruned'),
        algo=algoLabel,
        btagInfos=['None'],
        btagDiscriminators=['None'],
        jetCorrections=jetCorrectionsAK7,
        genJetCollection = cms.InputTag('genJetsNoNu'),
        getJetMCFlavour = False,
        postfix = postfix
    )
    getattr(process,'patJetPartonMatchPrunedPFCHS'+postfix).matched = cms.InputTag(genParticles)
    addJetCollection(
        process,
        labelName='PrunedSubjetsPFCHS',
        jetSource=cms.InputTag('PFJetsCHSPruned','SubJets'),
        algo=algoLabel,
        rParam=options.jetRadius,
        trackSource = cms.InputTag(trackSource),
        pvSource = cms.InputTag(pvSource),
        svSource = svSource,
        btagInfos=bTagInfos,
        btagDiscriminators=bTagDiscriminatorsSubJets,
        jetCorrections=jetCorrectionsAK4,
        genJetCollection = cms.InputTag('genJetsNoNuPruned','SubJets'),
        explicitJTA=(True and not options.miniAOD), # old-style explicit JTA does not work with miniAOD
        svClustering=options.useSVClustering,
        fatJets=cms.InputTag('PFJetsCHS'),
        groomedFatJets=cms.InputTag('PFJetsCHSPruned'),
        postfix = postfix
    )
    getattr(process,'patJetPartonMatchPrunedSubjetsPFCHS'+postfix).matched = cms.InputTag(genParticles)

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

    getattr(process,'patJetsPFCHS'+postfix).userData.userFloats.src += ['PFJetsCHSPrunedMass','PFJetsCHSPrunedPt']

    #-------------------------------------
    if options.useSVClustering:
        setattr(process,'inclusiveSecondaryVertexFinderTagInfosPrunedSubjetsPFCHS'+postfix, getattr(process,'inclusiveSecondaryVertexFinderTagInfosPrunedSubjetsPFCHS'+postfix).clone(
            useSVClustering = cms.bool(True),
            useSVMomentum   = cms.bool(options.useSVMomentum), # otherwise using SV flight direction
            jetAlgorithm    = cms.string(options.jetAlgo),
            rParam          = cms.double(options.jetRadius),
            ghostRescaling  = cms.double(1e-18),
            fatJets         = cms.InputTag("PFJetsCHS"),
            groomedFatJets  = cms.InputTag("PFJetsCHSPruned")
        ))

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
## Produce a collection of good primary vertices
process.load('CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi')
process.goodOfflinePrimaryVertices.src = cms.InputTag(pvSource)

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
adaptPVs(process, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'))

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
    # Set the cone size for the jet-track association to the jet radius
    getattr(process,'jetTracksAssociatorAtVertexPFCHS'+postfix).coneSize = cms.double(options.jetRadius) # default is 0.5
    getattr(process,'secondaryVertexTagInfosPFCHS'+postfix).trackSelection.jetDeltaRMax = cms.double(options.jetRadius)   # default is 0.3
    getattr(process,'secondaryVertexTagInfosPFCHS'+postfix).vertexCuts.maxDeltaRToJetAxis = cms.double(options.jetRadius) # default is 0.5
    # Set the jet-SV dR to the jet radius
    getattr(process,'inclusiveSecondaryVertexFinderTagInfosPFCHS'+postfix).vertexCuts.maxDeltaRToJetAxis = cms.double(options.jetRadius) # default is 0.5
    getattr(process,'inclusiveSecondaryVertexFinderTagInfosPFCHS'+postfix).extSVDeltaRToJet = cms.double(options.jetRadius) # default is 0.3
    # Set the JP track dR cut to the jet radius
    process.jetProbabilityComputerFat = process.jetProbabilityComputer.clone( deltaR = cms.double(options.jetRadius) ) # default is 0.3
    getattr(process,'jetProbabilityBJetTagsPFCHS'+postfix).jetTagComputer = cms.string('jetProbabilityComputerFat')
    # Set the JBP track dR cut to the jet radius
    process.jetBProbabilityComputerFat = process.jetBProbabilityComputer.clone( deltaR = cms.double(options.jetRadius) ) # default is 0.5
    getattr(process,'jetBProbabilityBJetTagsPFCHS'+postfix).jetTagComputer = cms.string('jetBProbabilityComputerFat')
    # Set the CSV track dR cut to the jet radius
    process.combinedSecondaryVertexComputerFat = process.combinedSecondaryVertexComputer.clone()
    process.combinedSecondaryVertexComputerFat.trackSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
    process.combinedSecondaryVertexComputerFat.trackPseudoSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
    getattr(process,'combinedSecondaryVertexBJetTagsPFCHS'+postfix).jetTagComputer = cms.string('combinedSecondaryVertexComputerFat')
    # Set the CSVv2 track dR cut to the jet radius
    process.combinedSecondaryVertexV2ComputerFat = process.combinedSecondaryVertexV2Computer.clone()
    process.combinedSecondaryVertexV2ComputerFat.trackSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
    process.combinedSecondaryVertexV2ComputerFat.trackPseudoSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
    getattr(process,'combinedInclusiveSecondaryVertexV2BJetTagsPFCHS'+postfix).jetTagComputer = cms.string('combinedSecondaryVertexV2ComputerFat')

#-------------------------------------
if not options.runOnData:
    ## JP calibrations:
    process.GlobalTag.toGet = cms.VPSet(
      cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
           tag = cms.string("TrackProbabilityCalibration_2D_MC53X_v2"),
           connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
      cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
           tag = cms.string("TrackProbabilityCalibration_3D_MC53X_v2"),
           connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
    )
#---------------------------------------

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
process.btagana.muonCollectionName     = cms.InputTag(muons)
process.btagana.triggerTable           = cms.InputTag('TriggerResults::HLT') # Data and MC
process.btagana.svComputer             = cms.string('combinedSecondaryVertexV2Computer')
process.btagana.svComputerFatJets      = cms.string('combinedSecondaryVertexV2Computer')
process.btagana.svTagInfos             = cms.string('inclusiveSecondaryVertexFinder')

if options.runSubJets:
    process.btaganaSubJets = process.btagana.clone(
        storeEventInfo      = cms.bool(not options.processStdAK4Jets),
        allowJetSkipping    = cms.bool(False),
        Jets                = cms.InputTag('selectedPatJetsPrunedSubjetsPFCHS'+postfix),
        FatJets             = cms.InputTag('selectedPatJetsPFCHS'+postfix),
        GroomedFatJets      = cms.InputTag('selectedPatJetsPrunedPFCHSPacked'),
        runSubJets          = options.runSubJets,
        svComputer          = cms.string('combinedSecondaryVertexV2Computer'),
        svComputerFatJets   = cms.string('combinedSecondaryVertexV2ComputerFat'),
        svTagInfos          = cms.string('inclusiveSecondaryVertexFinder')
    )

#---------------------------------------

#---------------------------------------
## Define event filter sequence
process.filtSeq = cms.Sequence(
    #process.noscraping
    process.primaryVertexFilter
    * process.goodOfflinePrimaryVertices
)


## Define analyzer sequence
process.analyzerSeq = cms.Sequence( )
if options.processStdAK4Jets:
    process.analyzerSeq += process.btagana
if options.runSubJets:
    process.analyzerSeq += process.btaganaSubJets
#---------------------------------------

process.p = cms.Path(
    process.filtSeq
    * process.analyzerSeq
)

# Delete predefined output module (needed for running with CRAB)
del process.out

#open('pydump.py','w').write(process.dumpPython())
