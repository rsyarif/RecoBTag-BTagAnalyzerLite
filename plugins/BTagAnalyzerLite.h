#ifndef BTagAnalyzerLite_h
#define BTagAnalyzerLite_h

// -*- C++ -*-
//
// Package:    BTagAnalyzerLite
// Class:      BTagAnalyzerLite
//
/**\class BTagAnalyzerLite BTagAnalyzerLite.cc RecoBTag/BTagAnalyzerLite/plugins/BTagAnalyzerLite.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Andrea Jeremy
//         Created:  Thu Dec 20 10:00:00 CEST 2012
//
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RegexMatch.h"
#include "FWCore/ParameterSet/interface/FileInPath.h" // added by rizki

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/SoftLeptonTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/GeometrySurface/interface/Line.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputerWrapper.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputerRecord.h"
#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"

//added by rizki - start
#include <fastjet/PseudoJet.hh>
#include <fastjet/FunctionOfPseudoJet.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "ShowerDeconstruction/SDAlgorithm/interface/Exception.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/AnalysisParameters.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/HBBModel.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/BackgroundModel.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/ISRModel.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/Deconstruct.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/Message.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/Parameters.h"
#include "ShowerDeconstruction/SDAlgorithm/interface/ParseUtils.h"
using namespace Deconstruction;
//added by rizki - end

#include "TFile.h"
#include "TTree.h"

#include <boost/regex.hpp>

#include "RecoBTag/BTagAnalyzerLite/interface/JetInfoBranches.h"
#include "RecoBTag/BTagAnalyzerLite/interface/EventInfoBranches.h"

//
// constants, enums and typedefs
//
typedef std::vector<pat::Jet> PatJetCollection;

//
// class declaration
//

struct orderByPt {
    const std::string mCorrLevel;
    orderByPt(const std::string& fCorrLevel) : mCorrLevel(fCorrLevel) {}
    bool operator ()(PatJetCollection::const_iterator const& a, PatJetCollection::const_iterator const& b) {
      if( mCorrLevel=="Uncorrected" )
        return a->correctedJet("Uncorrected").pt() > b->correctedJet("Uncorrected").pt();
      else
        return a->pt() > b->pt();
    }
};

using namespace std;
using namespace reco;

const UInt_t MAX_JETCOLLECTIONS=2;

class BTagAnalyzerLite : public edm::EDAnalyzer
{
  public:
    explicit BTagAnalyzerLite(const edm::ParameterSet&);
    ~BTagAnalyzerLite();

  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    void setTracksPV( const reco::TrackRef & trackRef, const edm::Handle<reco::VertexCollection> & pvHandle, int & iPV, float & PVweight );

    void setTracksSV( const reco::TrackRef & trackRef, const reco::SecondaryVertexTagInfo *, int & isFromSV, int & iSV, float & SVweight );

    bool NameCompatible(const std::string& pattern, const std::string& name);

    void processTrig(const edm::Handle<edm::TriggerResults>&, const std::vector<std::string>&) ;

    void processJets(const edm::Handle<PatJetCollection>&, const edm::Handle<PatJetCollection>&,
                     const edm::Event&, const edm::EventSetup&,
                     const edm::Handle<PatJetCollection>&, std::vector<int>&, const int) ;

    bool isHardProcess(const int status);

    void matchGroomedJets(const edm::Handle<PatJetCollection>& jets,
                          const edm::Handle<PatJetCollection>& matchedJets,
                          std::vector<int>& matchedIndices);

    // ----------member data ---------------------------
    std::string outputFile_;
    //std::vector< std::string > moduleLabel_;

    bool runSubJets_ ;
    bool allowJetSkipping_ ;
    bool storeEventInfo_;
    bool produceJetTrackTree_;
    bool produceJetPFLeptonTree_;
    bool storeMuonInfo_;
    bool storeTagVariables_;
    bool storeCSVTagVariables_;

    double microjetConesize_; //SD parameter added by rizki
    edm::FileInPath SDinputcard_; //added by rizki

    edm::InputTag src_;  // Generator/handronizer module label
    edm::InputTag muonCollectionName_;
    edm::InputTag prunedGenParticleCollectionName_;
    edm::InputTag triggerTable_;

    edm::InputTag JetCollectionTag_;
    edm::InputTag FatJetCollectionTag_;
    edm::InputTag GroomedFatJetCollectionTag_;

    edm::InputTag primaryVertexColl_;

    std::string jetPBJetTags_;
    std::string jetPNegBJetTags_;
    std::string jetPPosBJetTags_;

    std::string jetBPBJetTags_;
    std::string jetBPNegBJetTags_;
    std::string jetBPPosBJetTags_;

    std::string trackCHEBJetTags_;
    std::string trackCNegHEBJetTags_;

    std::string trackCHPBJetTags_;
    std::string trackCNegHPBJetTags_;

    std::string simpleSVHighEffBJetTags_;
    std::string simpleSVNegHighEffBJetTags_;
    std::string simpleSVHighPurBJetTags_;
    std::string simpleSVNegHighPurBJetTags_;

    std::string simpleIVFSVHighPurBJetTags_;
    std::string simpleIVFSVHighEffBJetTags_;
    std::string doubleIVFSVHighEffBJetTags_;

    std::string combinedSVBJetTags_;
    std::string combinedSVNegBJetTags_;
    std::string combinedSVPosBJetTags_;

    std::string combinedIVFSVBJetTags_;
    std::string combinedIVFSVPosBJetTags_;

    std::string softPFMuonBJetTags_;
    std::string softPFMuonNegBJetTags_;
    std::string softPFMuonPosBJetTags_;

    std::string softPFElectronBJetTags_;
    std::string softPFElectronNegBJetTags_;
    std::string softPFElectronPosBJetTags_;

    std::string ipTagInfos_;
    std::string svTagInfos_;
    std::string softPFMuonTagInfos_;
    std::string softPFElectronTagInfos_;
    std::string ivfTagInfos_; //added by rizki

    std::string   SVComputer_;
    std::string   SVComputerFatJets_;

    TFile*  rootFile_;
    double minJetPt_;
    double maxJetEta_;

    bool isData_;

    // trigger list
    std::vector<std::string> triggerPathNames_;

    edm::Service<TFileService> fs;

    ///////////////
    // Ntuple info

    TTree *smalltree;

    //// Event info
    EventInfoBranches EventInfo;

    //// Jet info
    JetInfoBranches JetInfo[MAX_JETCOLLECTIONS] ;

    edm::Handle<reco::VertexCollection> primaryVertex;

    const reco::Vertex *pv;
    const GenericMVAJetTagComputer *computer;

    // Generator/hadronizer type (information stored bitwise)
    unsigned int hadronizerType_;
};

#endif
