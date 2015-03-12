import FWCore.ParameterSet.Config as cms

bTagAnalyzerLiteCommon = cms.PSet(
    runSubJets               = cms.bool(False),
    allowJetSkipping         = cms.bool(True),
    storeEventInfo           = cms.bool(True),
    produceJetTrackTree      = cms.bool(False), ## True if you want to keep info for tracks associated to jets
    produceJetPFLeptonTree   = cms.bool(False), ## True if you want to keep PF lepton info
    storeMuonInfo            = cms.bool(False), ## True if you want to keep muon info
    storeTagVariables        = cms.bool(False), ## True if you want to keep TagInfo TaggingVariables
    storeCSVTagVariables     = cms.bool(True),  ## True if you want to keep CSV TaggingVariables
    MaxEta                   = cms.double(2.5),
    MinPt                    = cms.double(20.0),
    src                      = cms.InputTag('generator'),
    Jets                     = cms.InputTag('selectedPatJets'),
    FatJets                  = cms.InputTag('selectedPatJets'),
    GroomedFatJets           = cms.InputTag('selectedPatJetsAK8PrunedPFPacked'),
    muonCollectionName       = cms.InputTag('selectedPatMuons'),
    triggerTable             = cms.InputTag('TriggerResults'),
    prunedGenParticles       = cms.InputTag('prunedGenParticlesBoost'),
    primaryVertexColl        = cms.InputTag('offlinePrimaryVertices'),
    TriggerPathNames = cms.vstring(
        "HLT_HT750_v*"
    )
)
