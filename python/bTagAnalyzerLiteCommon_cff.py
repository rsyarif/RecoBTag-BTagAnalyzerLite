import FWCore.ParameterSet.Config as cms

bTagAnalyzerLiteCommon = cms.PSet(
    runFatJets               = cms.bool(False),
    runSubJets               = cms.bool(False),
    allowJetSkipping         = cms.bool(True),
    storeEventInfo           = cms.bool(True),
    produceJetTrackTree      = cms.bool(False), ## True if you want to keep info for tracks associated to jets
    produceJetPFLeptonTree   = cms.bool(False), ## True if you want to keep PF lepton info
    storeTagVariables        = cms.bool(False), ## True if you want to keep TagInfo TaggingVariables
    storeCSVTagVariables     = cms.bool(True),  ## True if you want to keep CSV TaggingVariables
    MaxEta                   = cms.double(2.5),
    MinPt                    = cms.double(20.0),
    src                      = cms.InputTag('generator'),
    BranchNamePrefix         = cms.string(''),
    Jets                     = cms.InputTag('selectedPatJets'),
    SubJets                  = cms.VInputTag(),
    SubJetLabels             = cms.vstring(),
    muonCollectionName       = cms.InputTag('muons'),
    triggerTable             = cms.InputTag('TriggerResults'),
    prunedGenParticles       = cms.InputTag('prunedGenParticlesBoost'),
    primaryVertexColl        = cms.InputTag('offlinePrimaryVertices'),
    TriggerPathNames = cms.vstring(
        "HLT_HT750_v*"
    )
)
