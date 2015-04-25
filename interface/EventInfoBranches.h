#ifndef EVENTINFOBRANCHES_H
#define EVENTINFOBRANCHES_H

#include <TTree.h>

const UInt_t nMaxPVs_= 1000;
const UInt_t nMaxPUs_= 1000;
const UInt_t nMaxTrkAll_ = 100000;

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
};

#endif

