#ifndef EVENTINFOBRANCHES_H
#define EVENTINFOBRANCHES_H

#include <TTree.h>

const UInt_t nMaxPVs_= 1000;
const UInt_t nMaxPUs_= 1000;
const UInt_t nMaxTrkAll_ = 100000;
//const UInt_t nMaxSVs_= 10000;
class EventInfoBranches {

  public :

    int   nBitTrigger;
    int   BitTrigger[100];
    int   Run;
    int   Evt;
    int   LumiBlock;
    float PVz;
    float PVez;
    float GenPVz;
    float pthat;
    float mcweight;

    int   nPV;
    float PV_x[nMaxPVs_];
    float PV_y[nMaxPVs_];
    float PV_z[nMaxPVs_];
    float PV_ex[nMaxPVs_];
    float PV_ey[nMaxPVs_];
    float PV_ez[nMaxPVs_];
    float PV_chi2[nMaxPVs_];
    float PV_ndf[nMaxPVs_];
    int   PV_isgood[nMaxPVs_];
    int   PV_isfake[nMaxPVs_];

    float nPUtrue;                 // the true number of pileup interactions that have been added to the event
    int   nPU;                     // the number of pileup interactions that have been added to the event
    int   PU_bunch[nMaxPUs_];      // 0 if on time pileup, -1 or +1 if out-of-time
    float PU_z[nMaxPUs_];          // the true primary vertex position along the z axis for each added interaction
    float PU_sumpT_low[nMaxPUs_];  // the sum of the transverse momentum of the tracks originating from each interaction, where track pT > low_cut
    float PU_sumpT_high[nMaxPUs_]; // the sum of the transverse momentum of the tracks originating from each interaction, where track pT > high_cut
    int   PU_ntrks_low[nMaxPUs_];  // the number of tracks originating from each interaction, where track pT > low_cu
    int   PU_ntrks_high[nMaxPUs_]; // the number of tracks originating from each interaction, where track pT > high_cut

    int   nGenPruned;
    float GenPruned_pT[1000];
    float GenPruned_eta[1000];
    float GenPruned_phi[1000];
    float GenPruned_mass[1000];
    int   GenPruned_status[1000];
    int   GenPruned_pdgID[1000];
    int   GenPruned_mother[1000];

    int   nMuon;
    int   Muon_isGlobal[1000];
    int   Muon_isPF[1000];
    int   Muon_nTkHit[1000];
    int   Muon_nPixHit[1000];
    int   Muon_nOutHit[1000];
    int   Muon_nMuHit[1000];
    int   Muon_nMatched[1000];
    float Muon_chi2[1000];
    float Muon_chi2Tk[1000];
    float Muon_pt[1000];
    float Muon_eta[1000];
    float Muon_phi[1000];
    float Muon_vz[1000];
    float Muon_IP[1000];
    float Muon_IPsig[1000];
    float Muon_IP2D[1000];
    float Muon_IP2Dsig[1000];


	


    void RegisterTree(TTree *tree) {
      tree->Branch("nBitTrigger", &nBitTrigger,  "nBitTrigger/I");
      tree->Branch("BitTrigger" , BitTrigger  ,  "BitTrigger[nBitTrigger]/I");
      tree->Branch("Run"        , &Run        ,  "Run/I");
      tree->Branch("Evt"        , &Evt        ,  "Evt/I");
      tree->Branch("LumiBlock"  , &LumiBlock  ,  "LumiBlock/I");
      tree->Branch("pthat"      , &pthat      ,  "pthat/F");
      tree->Branch("mcweight"   , &mcweight   ,  "mcweight/F");
      tree->Branch("nPV"        , &nPV        ,  "nPV/I");
      tree->Branch("PVz"        , &PVz        ,  "PVz/F");
      tree->Branch("PVez"       , &PVez       ,  "PVez/F");
      tree->Branch("GenPVz"     , &GenPVz     ,  "GenPVz/F");

      tree->Branch("nPUtrue"      , &nPUtrue     , "nPUtrue/F");
      tree->Branch("nPU"          , &nPU         , "nPU/I"    );
      tree->Branch("PU_bunch"     , PU_bunch     , "PU_bunch[nPU]/I");
      tree->Branch("PU_z"         , PU_z         , "PU_z[nPU]/F");
      tree->Branch("PU_sumpT_low" , PU_sumpT_low , "PU_sumpT_low[nPU]/F");
      tree->Branch("PU_sumpT_high", PU_sumpT_high, "PU_sumpT_high[nPU]/F");
      tree->Branch("PU_ntrks_low" , PU_ntrks_low , "PU_ntrks_low[nPU]/I");
      tree->Branch("PU_ntrks_high", PU_ntrks_high, "PU_ntrks_high[nPU]/I");

      tree->Branch("nGenPruned",     &nGenPruned       ,"nGenPruned/I");
      tree->Branch("GenPruned_pT",     GenPruned_pT    , "GenPruned_pT[nGenPruned]/F");
      tree->Branch("GenPruned_eta",    GenPruned_eta   , "GenPruned_eta[nGenPruned]/F");
      tree->Branch("GenPruned_phi",    GenPruned_phi   , "GenPruned_phi[nGenPruned]/F");
      tree->Branch("GenPruned_mass",    GenPruned_mass   , "GenPruned_mass[nGenPruned]/F");
      tree->Branch("GenPruned_pdgID",  GenPruned_pdgID , "GenPruned_pdgID[nGenPruned]/I");
      tree->Branch("GenPruned_status", GenPruned_status, "GenPruned_status[nGenPruned]/I");
      tree->Branch("GenPruned_mother", GenPruned_mother, "GenPruned_mother[nGenPruned]/I");
    }

    void RegisterJetTrackTree(TTree *tree) {
      tree->Branch("PV_x"     , PV_x     , "PV_x[nPV]/F");
      tree->Branch("PV_y"     , PV_y     , "PV_y[nPV]/F");
      tree->Branch("PV_z"     , PV_z     , "PV_z[nPV]/F");
      tree->Branch("PV_ex"    , PV_ex    , "PV_ex[nPV]/F");
      tree->Branch("PV_ey"    , PV_ey    , "PV_ey[nPV]/F");
      tree->Branch("PV_ez"    , PV_ez    , "PV_ez[nPV]/F");
      tree->Branch("PV_chi2"  , PV_chi2  , "PV_chi2[nPV]/F");
      tree->Branch("PV_ndf"   , PV_ndf   , "PV_ndf[nPV]/F");
      tree->Branch("PV_isgood", PV_isgood, "PV_isgood[nPV]/I");
      tree->Branch("PV_isfake", PV_isfake, "PV_isfake[nPV]/I");
    }

    void RegisterMuonTree(TTree *tree) {
      tree->Branch("nMuon"        , &nMuon       , "nMuon/I");
      tree->Branch("Muon_nMuHit"  , Muon_nMuHit  , "Muon_nMuHit[nMuon]/I");
      tree->Branch("Muon_nTkHit"  , Muon_nTkHit  , "Muon_nTkHit[nMuon]/I");
      tree->Branch("Muon_nPixHit" , Muon_nPixHit , "Muon_nPixHit[nMuon]/I");
      tree->Branch("Muon_nOutHit" , Muon_nOutHit , "Muon_nOutHit[nMuon]/I");
      tree->Branch("Muon_isGlobal", Muon_isGlobal, "Muon_isGlobal[nMuon]/I");
      tree->Branch("Muon_isPF"    , Muon_isPF    , "Muon_isPF[nMuon]/I");
      tree->Branch("Muon_nMatched", Muon_nMatched, "Muon_nMatched[nMuon]/I");
      tree->Branch("Muon_chi2"    , Muon_chi2    , "Muon_chi2[nMuon]/F");
      tree->Branch("Muon_chi2Tk"  , Muon_chi2Tk  , "Muon_chi2Tk[nMuon]/F");
      tree->Branch("Muon_pt"      , Muon_pt      , "Muon_pt[nMuon]/F");
      tree->Branch("Muon_eta"     , Muon_eta     , "Muon_eta[nMuon]/F");
      tree->Branch("Muon_phi"     , Muon_phi     , "Muon_phi[nMuon]/F");
      tree->Branch("Muon_vz"      , Muon_vz      , "Muon_vz[nMuon]/F");
      tree->Branch("Muon_IP"      , Muon_IP      , "Muon_IP[nMuon]/F");
      tree->Branch("Muon_IPsig"   , Muon_IPsig   , "Muon_IPsig[nMuon]/F");
      tree->Branch("Muon_IP2D"    , Muon_IP2D    , "Muon_IP2D[nMuon]/F");
      tree->Branch("Muon_IP2Dsig" , Muon_IP2Dsig , "Muon_IP2Dsig[nMuon]/F");
    }

    //------------------------------------------------------------------------------------------------------------------

    void ReadTree(TTree *tree) {
      tree->SetBranchAddress("nBitTrigger", &nBitTrigger);
      tree->SetBranchAddress("BitTrigger" , BitTrigger  );
      tree->SetBranchAddress("Run"        , &Run        );
      tree->SetBranchAddress("Evt"        , &Evt        );
      tree->SetBranchAddress("LumiBlock"  , &LumiBlock  );
      tree->SetBranchAddress("pthat"      , &pthat      );
      tree->SetBranchAddress("mcweight"   , &mcweight   );
      tree->SetBranchAddress("nPV"        , &nPV        );
      tree->SetBranchAddress("PVz"        , &PVz        );
      tree->SetBranchAddress("PVez"       , &PVez       );
      tree->SetBranchAddress("GenPVz"     , &GenPVz     );

      tree->SetBranchAddress("nPUtrue"      , &nPUtrue     );
      tree->SetBranchAddress("nPU"          , &nPU         );
      tree->SetBranchAddress("PU_bunch"     , PU_bunch     );
      tree->SetBranchAddress("PU_z"         , PU_z         );
      tree->SetBranchAddress("PU_sumpT_low" , PU_sumpT_low );
      tree->SetBranchAddress("PU_sumpT_high", PU_sumpT_high);
      tree->SetBranchAddress("PU_ntrks_low" , PU_ntrks_low );
      tree->SetBranchAddress("PU_ntrks_high", PU_ntrks_high);

      tree->SetBranchAddress("nGenPruned",       &nGenPruned     );
      tree->SetBranchAddress("GenPruned_pT",     GenPruned_pT    );
      tree->SetBranchAddress("GenPruned_eta",    GenPruned_eta   );
      tree->SetBranchAddress("GenPruned_phi",    GenPruned_phi   );
      tree->SetBranchAddress("GenPruned_mass",    GenPruned_mass   );
      tree->SetBranchAddress("GenPruned_pdgID",  GenPruned_pdgID );
      tree->SetBranchAddress("GenPruned_status", GenPruned_status);
      tree->SetBranchAddress("GenPruned_mother", GenPruned_mother);
    }

    void ReadJetTrackTree(TTree *tree) {
      tree->SetBranchAddress("PV_x"     , PV_x     );
      tree->SetBranchAddress("PV_y"     , PV_y     );
      tree->SetBranchAddress("PV_z"     , PV_z     );
      tree->SetBranchAddress("PV_ex"    , PV_ex    );
      tree->SetBranchAddress("PV_ey"    , PV_ey    );
      tree->SetBranchAddress("PV_ez"    , PV_ez    );
      tree->SetBranchAddress("PV_chi2"  , PV_chi2  );
      tree->SetBranchAddress("PV_ndf"   , PV_ndf   );
      tree->SetBranchAddress("PV_isgood", PV_isgood);
      tree->SetBranchAddress("PV_isfake", PV_isfake);
    }

    void ReadMuonTree(TTree *tree) {
      tree->SetBranchAddress("nMuon"        , &nMuon       );
      tree->SetBranchAddress("Muon_nMuHit"  , Muon_nMuHit  );
      tree->SetBranchAddress("Muon_nTkHit"  , Muon_nTkHit  );
      tree->SetBranchAddress("Muon_nPixHit" , Muon_nPixHit );
      tree->SetBranchAddress("Muon_nOutHit" , Muon_nOutHit );
      tree->SetBranchAddress("Muon_isGlobal", Muon_isGlobal);
      tree->SetBranchAddress("Muon_isPF"    , Muon_isPF    );
      tree->SetBranchAddress("Muon_nMatched", Muon_nMatched);
      tree->SetBranchAddress("Muon_chi2"    , Muon_chi2    );
      tree->SetBranchAddress("Muon_chi2Tk"  , Muon_chi2Tk  );
      tree->SetBranchAddress("Muon_pt"      , Muon_pt      );
      tree->SetBranchAddress("Muon_eta"     , Muon_eta     );
      tree->SetBranchAddress("Muon_phi"     , Muon_phi     );
      tree->SetBranchAddress("Muon_vz"      , Muon_vz      );
      tree->SetBranchAddress("Muon_IP"      , Muon_IP      );
      tree->SetBranchAddress("Muon_IPsig"   , Muon_IPsig   );
      tree->SetBranchAddress("Muon_IP2D"    , Muon_IP2D    );
      tree->SetBranchAddress("Muon_IP2Dsig" , Muon_IP2Dsig );

/*      tree->SetBranchAddress("nSV"              ,&nSV               ) ;
      tree->SetBranchAddress("SV_x"             ,SV_x                     ) ;
      tree->SetBranchAddress("SV_y"             ,SV_y                     ) ;
      tree->SetBranchAddress("SV_z"             ,SV_z                     ) ;
      tree->SetBranchAddress("SV_ex"            ,SV_ex                  ) ;
      tree->SetBranchAddress("SV_ey"            ,SV_ey                  ) ;
      tree->SetBranchAddress("SV_ez"            ,SV_ez                  ) ;
      tree->SetBranchAddress("SV_chi2"          ,SV_chi2            ) ;
      tree->SetBranchAddress("SV_ndf"           ,SV_ndf                 ) ;
      tree->SetBranchAddress("SV_flight"        ,SV_flight          ) ;
      tree->SetBranchAddress("SV_flightErr"     ,SV_flightErr       ) ;
      tree->SetBranchAddress("SV_deltaR_jet"    ,SV_deltaR_jet      ) ;
      tree->SetBranchAddress("SV_deltaR_sum_jet",SV_deltaR_sum_jet  ) ;
      tree->SetBranchAddress("SV_deltaR_sum_dir",SV_deltaR_sum_dir  ) ;
      tree->SetBranchAddress("SV_vtx_pt"        ,SV_vtx_pt          ) ;
      tree->SetBranchAddress("SV_flight2D"      ,SV_flight2D        ) ;
      tree->SetBranchAddress("SV_flight2DErr"   ,SV_flight2DErr     ) ;
      tree->SetBranchAddress("SV_totCharge"     ,SV_totCharge       ) ;
      tree->SetBranchAddress("SV_vtxDistJetAxis",SV_vtxDistJetAxis  ) ;
      tree->SetBranchAddress("SV_nTrk"          ,SV_nTrk            ) ;
      tree->SetBranchAddress("SV_mass"          ,SV_mass            ) ;
      tree->SetBranchAddress("SV_vtx_eta"       ,SV_vtx_eta         ) ;
      tree->SetBranchAddress("SV_vtx_phi"       ,SV_vtx_phi         ) ;
*/
     
    }
};

#endif

