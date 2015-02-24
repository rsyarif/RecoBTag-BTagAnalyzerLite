import FWCore.ParameterSet.Config as cms

btagana = cms.EDAnalyzer("BTagAnalyzerLite",
    runSubJets               = cms.bool(False),
    allowJetSkipping         = cms.bool(True),
    storeEventInfo           = cms.bool(True),
    produceJetTrackTree      = cms.bool(False),
    produceJetPFLeptonTree   = cms.bool(False),
    storeMuonInfo            = cms.bool(False),
    storeTagVariables        = cms.bool(False),
    storeCSVTagVariables     = cms.bool(True),
    microjetConesize         = cms.double(0.15), #SD parameter added by rizki
    MaxEta                   = cms.double(2.5),
    MinPt                    = cms.double(20.0),
    src                      = cms.InputTag('generator'),
    Jets                     = cms.InputTag('ak5PFJets'),
    FatJets                  = cms.InputTag('selectedPatJets'),
    GroomedFatJets           = cms.InputTag('selectedPatJetsCA8PrunedPFPacked'),
    muonCollectionName       = cms.InputTag('selectedPatMuons'),
    triggerTable             = cms.InputTag('TriggerResults'),
    prunedGenParticles       = cms.InputTag('prunedGenParticlesBoost'),
    primaryVertexColl        = cms.InputTag('offlinePrimaryVertices'),

    svComputer               = cms.string('combinedSecondaryVertex'),
    svComputerFatJets        = cms.string('combinedSecondaryVertex'),
    SDinputcard              = cms.FileInPath('ShowerDeconstruction/inputdata/input_card.dat'), #SD parameter added by rizki

    # list of taggers
    trackCHEBJetTags      = cms.string('trackCountingHighEffBJetTags'),
    trackCNegHEBJetTags   = cms.string('negativeTrackCountingHighEffJetTags'),

    trackCHPBJetTags      = cms.string('trackCountingHighPurBJetTags'),
    trackCNegHPBJetTags   = cms.string('negativeTrackCountingHighPurJetTags'),

    jetBPBJetTags          = cms.string('jetBProbabilityBJetTags'),
    jetBPNegBJetTags       = cms.string('negativeOnlyJetBProbabilityJetTags'),
    jetBPPosBJetTags       = cms.string('positiveOnlyJetBProbabilityJetTags'),

    jetPBJetTags          = cms.string('jetProbabilityBJetTags'),
    jetPNegBJetTags       = cms.string('negativeOnlyJetProbabilityJetTags'),
    jetPPosBJetTags       = cms.string('positiveOnlyJetProbabilityJetTags'),

    simpleSVHighPurBJetTags     = cms.string('simpleSecondaryVertexHighPurBJetTags'),
    simpleSVNegHighPurBJetTags  = cms.string('simpleSecondaryVertexNegativeHighPurBJetTags'),
    simpleSVHighEffBJetTags     = cms.string('simpleSecondaryVertexHighEffBJetTags'),
    simpleSVNegHighEffBJetTags  = cms.string('simpleSecondaryVertexNegativeHighEffBJetTags'),

    simpleIVFSVHighPurBJetTags = cms.string('simpleInclusiveSecondaryVertexHighPurBJetTags'),
    simpleIVFSVHighEffBJetTags = cms.string('simpleInclusiveSecondaryVertexHighEffBJetTags'),
    doubleIVFSVHighEffBJetTags = cms.string('doubleSecondaryVertexHighEffBJetTags'),

    combinedSVBJetTags    = cms.string('combinedSecondaryVertexBJetTags'),
    combinedSVNegBJetTags = cms.string('combinedSecondaryVertexNegativeBJetTags'),
    combinedSVPosBJetTags = cms.string('combinedSecondaryVertexPositiveBJetTags'),

    combinedIVFSVBJetTags      = cms.string('combinedInclusiveSecondaryVertexV2BJetTags'),
    combinedIVFSVPosBJetTags   = cms.string('combinedInclusiveSecondaryVertexPositiveBJetTags'),

    softPFMuonBJetTags        = cms.string('softPFMuonBJetTags'),
    softPFMuonNegBJetTags     = cms.string('negativeSoftPFMuonBJetTags'),
    softPFMuonPosBJetTags     = cms.string('positiveSoftPFMuonBJetTags'),

    softPFElectronBJetTags    = cms.string('softPFElectronBJetTags'),
    softPFElectronNegBJetTags = cms.string('negativeSoftPFElectronBJetTags'),
    softPFElectronPosBJetTags = cms.string('positiveSoftPFElectronBJetTags'),

    ipTagInfos               = cms.string('impactParameter'), # need to omit the 'TagInfos' part from the label
    svTagInfos               = cms.string('secondaryVertex'), # need to omit the 'TagInfos' part from the label
    softPFMuonTagInfos       = cms.string('softPFMuons'),     # need to omit the 'TagInfos' part from the label
    softPFElectronTagInfos   = cms.string('softPFElectrons'), # need to omit the 'TagInfos' part from the label

    TriggerPathNames = cms.vstring(
        "HLT_HT750_v*"
    )
)
