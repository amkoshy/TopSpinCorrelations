#ifndef MINITREE_H
#define MINITREE_H

#include <TTree.h>
#include <TBranch.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include "analysisStructs.h"
#include "AnalysisConfig.h"
#include "../../common/include/classesFwd.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/AnalysisBase.h"
#include "../../common/include/analysisObjectStructs.h"

#include "../../common/include/KinematicReconstructionSolution.h"
#include "../../common/include/KinematicReconstruction.h"
#include "../../common/include/LooseKinRecoSolution.h"
#include "../../common/include/LooseKinReco.h"

#include "../../common/include/ScaleFactors.h"

#ifdef __MAKECINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

const UInt_t length_ = 1000;

// Simple structure to store info about trees which need to be created later
struct TreeSetup
{
  bool    isGen           = false;
  TString name            = "defaultSetupName";
  TString title           = "defaultSetupName";
  bool    isTTbarSample   = false;
  bool    isSystSample    = false;
  bool    isSystVariation = false;
};



struct AnalysisEntry
{
    const int* eventCounter;
    int eventCounterCorrected;
    const AnalysisBase* analysisBase;
    const AnalysisConfig* analysisConfig;
    const EventMetadata* eventMetadata;
    const RecoObjects* recoObjects;
    const CommonGenObjects* commonGenObjects;
    const TopGenObjects* topGenObjects;
    const TopPseudoObjects* topPseudoObjects;
    const ttbar::GenObjectIndices* genObjectIndices;
    const ttbar::GenLevelWeights* genLevelWeights;
    double weight;
    float pileupW;
    float pileupW_Up;
    float pileupW_Down;
    Long64_t entryNumber;
};

struct TreeVars
{
    //***********************************************************
    //******************** tree variables ***********************
    //***********************************************************
    double weight = 1.;
    // int eventCounter = 1;
    // int eventCounter;
    TBranch* b_eventCounter = NULL;

    unsigned long long eventNumber = 1;
    unsigned int lumiNumber = 1;
    unsigned int runNumber = 1;

    float leptonSF = 1.;
    float l1PrefiringWeight = 1.;
    float triggerSF = 1.;
    float btagSF = 1.;
    float pileupSF = 1.;

    int channel = 0;

    int nVtx     = -999;
    int nLeptons = -999;
    int nJets    = -999;
    int nBJets   = -999;

    float lepton1_pt  = -999.;
    float lepton1_eta = -999.;
    float lepton1_phi = -999.;
    float lepton1_m   = -999.;
    int lepton1_pdgID = -999.;

    float lepton2_pt  = -999.;
    float lepton2_eta = -999.;
    float lepton2_phi = -999.;
    float lepton2_m   = -999.;
    int lepton2_pdgID = -999.;

    float met_pt  = -999.;
    float met_phi = -999.;
    float met_significance = -999.;

    uint n_jets;
    float jets_pt[length_];
    float jets_eta[length_];
    float jets_phi[length_];
    float jets_m[length_];
    float jets_charge[length_];
    uint n_bjets;
    float bjets_pt[length_];
    float bjets_eta[length_];
    float bjets_phi[length_];
    float bjets_m[length_];
    float bjets_charge[length_];
    uint n_nonbjets;
    float nonbjets_pt[length_];
    float nonbjets_eta[length_];
    float nonbjets_phi[length_];
    float nonbjets_m[length_];
    float nonbjets_charge[length_];

    bool jets_btag[length_];
    bool jets_topMatched[length_];
    bool jets_addJetMatched[length_];
    bool jets_addJetWithNuMatched[length_];


    bool passStep3 = false;
    bool passStep4 = false;
    bool passStep5 = false;
    bool passStep6 = false;
    bool passStep7 = false;
    bool passStep8 = false;
    bool passStep8Loose = false;

    bool  hasKinRecoSolution     = false;
    float kinReco_top_pt         = -999;
    float kinReco_top_eta        = -999;
    float kinReco_top_phi        = -999;
    float kinReco_top_m          = -999;
    float kinReco_antitop_pt     = -999;
    float kinReco_antitop_eta    = -999;
    float kinReco_antitop_phi    = -999;
    float kinReco_antitop_m      = -999;
    float kinReco_lepton_pt      = -999;
    float kinReco_lepton_eta     = -999;
    float kinReco_lepton_phi     = -999;
    float kinReco_lepton_m       = -999;
    float kinReco_antilepton_pt  = -999;
    float kinReco_antilepton_eta = -999;
    float kinReco_antilepton_phi = -999;
    float kinReco_antilepton_m   = -999;
    float kinReco_b_pt           = -999;
    float kinReco_b_eta          = -999;
    float kinReco_b_phi          = -999;
    float kinReco_b_m            = -999;
    float kinReco_antib_pt       = -999;
    float kinReco_antib_eta      = -999;
    float kinReco_antib_phi      = -999;
    float kinReco_antib_m        = -999;
    float kinReco_weight        = 1.;
    uint n_kinReco_nonbjets;
    float kinReco_nonbjets_pt[length_];
    float kinReco_nonbjets_eta[length_];
    float kinReco_nonbjets_phi[length_];
    float kinReco_nonbjets_m[length_];
    float kinReco_nonbjet_pt  = -999.;
    float kinReco_nonbjet_eta = -999.;
    float kinReco_nonbjet_phi = -999.;
    float kinReco_nonbjet_m   = -999.;

    float kinReco_rho = -999;

    bool  hasLooseKinRecoSolution = false;
    float looseKinReco_ttbar_pt   = -999.;
    float looseKinReco_ttbar_eta  = -999.;
    float looseKinReco_ttbar_phi  = -999.;
    float looseKinReco_ttbar_m    = -999.;
    float looseKinReco_ww_m       = -999.;
    float looseKinReco_b1_pt      = -999.;
    float looseKinReco_b1_eta     = -999.;
    float looseKinReco_b1_phi     = -999.;
    float looseKinReco_b1_m       = -999.;
    float looseKinReco_b2_pt      = -999.;
    float looseKinReco_b2_eta     = -999.;
    float looseKinReco_b2_phi     = -999.;
    float looseKinReco_b2_m       = -999.;
    float looseKinReco_l1_pt      = -999.;
    float looseKinReco_l1_eta     = -999.;
    float looseKinReco_l1_phi     = -999.;
    float looseKinReco_l1_m       = -999.;
    float looseKinReco_l2_pt      = -999.;
    float looseKinReco_l2_eta     = -999.;
    float looseKinReco_l2_phi     = -999.;
    float looseKinReco_l2_m       = -999.;
    float looseKinReco_weight    = 1.;
    uint n_looseKinReco_nonbjets;
    float looseKinReco_nonbjets_pt[length_];
    float looseKinReco_nonbjets_eta[length_];
    float looseKinReco_nonbjets_phi[length_];
    float looseKinReco_nonbjets_m[length_];
    float looseKinReco_nonbjet_pt  = -999.;
    float looseKinReco_nonbjet_eta = -999.;
    float looseKinReco_nonbjet_phi = -999.;
    float looseKinReco_nonbjet_m   = -999.;

    float looseKinReco_rho = -999;

    uint n_gen_additional_jets;
    float gen_additional_jets_pt[length_];
    float gen_additional_jets_eta[length_];
    float gen_additional_jets_phi[length_];
    float gen_additional_jets_m[length_];

    uint n_gen_additional_jets_withNu;
    float gen_additional_jets_withNu_pt[length_];
    float gen_additional_jets_withNu_eta[length_];
    float gen_additional_jets_withNu_phi[length_];
    float gen_additional_jets_withNu_m[length_];

    float gen_additional_jet_pt  = -999.;
    float gen_additional_jet_eta = -999.;
    float gen_additional_jet_phi = -999.;
    float gen_additional_jet_m   = -999.;

    float gen_additional_jet_withNu_pt  = -999.;
    float gen_additional_jet_withNu_eta = -999.;
    float gen_additional_jet_withNu_phi = -999.;
    float gen_additional_jet_withNu_m   = -999.;

    float gen_top_pt           = -999.;
    float gen_top_eta          = -999.;
    float gen_top_phi          = -999.;
    float gen_top_m            = -999.;
    float gen_antitop_pt       = -999.;
    float gen_antitop_eta      = -999.;
    float gen_antitop_phi      = -999.;
    float gen_antitop_m        = -999.;
    float gen_lepton_pt        = -999.;
    float gen_lepton_eta       = -999.;
    float gen_lepton_phi       = -999.;
    float gen_lepton_m         = -999.;
    float gen_antilepton_pt    = -999.;
    float gen_antilepton_eta   = -999.;
    float gen_antilepton_phi   = -999.;
    float gen_antilepton_m     = -999.;
    float gen_b_pt             = -999;
    float gen_b_eta            = -999;
    float gen_b_phi            = -999;
    float gen_b_m              = -999;
    float gen_antib_pt         = -999;
    float gen_antib_eta        = -999;
    float gen_antib_phi        = -999;
    float gen_antib_m          = -999;
    float gen_neutrino_pt      = -999;
    float gen_neutrino_eta     = -999;
    float gen_neutrino_phi     = -999;
    float gen_neutrino_m       = -999;
    float gen_antineutrino_pt  = -999;
    float gen_antineutrino_eta = -999;
    float gen_antineutrino_phi = -999;
    float gen_antineutrino_m   = -999;
    uint gen_nVtx               = -999;
    uint gen_nJets              = -999;
    uint gen_nBJets             = -999;

    float gen_met_pt  = -999.;
    float gen_met_eta = -999.;
    float gen_met_phi = -999.;
    float gen_met_m   = -999.;

    float gen_rho       = -999.;
    float gen_rhoWithNu = -999.;

    float pseudo_top_pt         = -999.;
    float pseudo_top_eta        = -999.;
    float pseudo_top_phi        = -999.;
    float pseudo_top_m          = -999.;
    float pseudo_antitop_pt     = -999.;
    float pseudo_antitop_eta    = -999.;
    float pseudo_antitop_phi    = -999.;
    float pseudo_antitop_m      = -999.;
    float pseudo_lepton_pt      = -999.;
    float pseudo_lepton_eta     = -999.;
    float pseudo_lepton_phi     = -999.;
    float pseudo_lepton_m       = -999.;
    float pseudo_antilepton_pt  = -999.;
    float pseudo_antilepton_eta = -999.;
    float pseudo_antilepton_phi = -999.;
    float pseudo_antilepton_m   = -999.;
    float pseudo_b_pt           = -999.;
    float pseudo_b_eta          = -999.;
    float pseudo_b_phi          = -999.;
    float pseudo_b_m            = -999.;
    float pseudo_antib_pt       = -999.;
    float pseudo_antib_eta      = -999.;
    float pseudo_antib_phi      = -999.;
    float pseudo_antib_m        = -999.;

    float pseudo_rho = -999.;

    uint n_pseudo_additional_jets;
    float pseudo_additional_jets_pt[length_];
    float pseudo_additional_jets_eta[length_];
    float pseudo_additional_jets_phi[length_];
    float pseudo_additional_jets_m[length_];

    float pseudo_additional_jet_pt  = -999.;
    float pseudo_additional_jet_eta = -999.;
    float pseudo_additional_jet_phi = -999.;
    float pseudo_additional_jet_m   = -999.;

    uint n_gen_partonLevel_additional_jets;
    float gen_partonLevel_additional_jets_pt[length_];
    float gen_partonLevel_additional_jets_eta[length_];
    float gen_partonLevel_additional_jets_phi[length_];
    float gen_partonLevel_additional_jets_m[length_];

    float gen_partonLevel_additional_jet_pt  = -999.;
    float gen_partonLevel_additional_jet_eta = -999.;
    float gen_partonLevel_additional_jet_phi = -999.;
    float gen_partonLevel_additional_jet_m   = -999.;

    float gen_partonLevel_rho = -999;


    float var_PU_UP   = -999.;
    float var_PU_DOWN = -999.;

    // triggers
    float var_TRIG_UP   = -999.;
    float var_TRIG_DOWN = -999.;

    // leptonID
    float var_ELE_ID_UP     = -999.;
    float var_ELE_ID_DOWN   = -999.;
    float var_ELE_RECO_UP   = -999.;
    float var_ELE_RECO_DOWN = -999.;
    float var_MUON_ID_UP    = -999.;
    float var_MUON_ID_DOWN  = -999.;
    float var_MUON_ISO_UP   = -999.;
    float var_MUON_ISO_DOWN = -999.;

    float var_L1PREFIRING_UP   = -999.;
    float var_L1PREFIRING_DOWN = -999.;
    // b-tagging
    float var_BTAG_UP        = -999.;
    float var_BTAG_DOWN      = -999.;
    float var_BTAG_LJET_UP   = -999.;
    float var_BTAG_LJET_DOWN = -999.;

    float var_BTAG_MERENSCALE_UP         = -999.;
    float var_BTAG_MERENSCALE_DOWN       = -999.;
    float var_BTAG_MEFACSCALE_UP         = -999.;
    float var_BTAG_MEFACSCALE_DOWN       = -999.;
    float var_BTAG_MESCALE_UP            = -999.;
    float var_BTAG_MESCALE_DOWN          = -999.;
    float var_BTAG_BFRAG_UP              = -999.;
    float var_BTAG_BFRAG_DOWN            = -999.;
    float var_BTAG_BFRAG_PETERSON        = -999.;
    float var_BTAG_BSEMILEP_UP           = -999.;
    float var_BTAG_BSEMILEP_DOWN         = -999.;
    float var_BTAG_PDF_ALPHAS_UP         = -999.;
    float var_BTAG_PDF_ALPHAS_DOWN       = -999.;
    float var_BTAG_PSSCALE_WEIGHT_4_UP   = -999.;
    float var_BTAG_PSSCALE_WEIGHT_4_DOWN = -999.;
    float var_BTAG_PSSCALE_WEIGHT_5_UP   = -999.;
    float var_BTAG_PSSCALE_WEIGHT_5_DOWN = -999.;
    float var_BTAG_PU_UP                 = -999.;
    float var_BTAG_PU_DOWN               = -999.;
    float var_BTAG_TRIG_UP               = -999.;
    float var_BTAG_TRIG_DOWN             = -999.;
    float var_BTAG_ELE_ID_UP             = -999.;
    float var_BTAG_ELE_ID_DOWN           = -999.;
    float var_BTAG_ELE_RECO_UP           = -999.;
    float var_BTAG_ELE_RECO_DOWN         = -999.;
    float var_BTAG_MUON_ID_UP            = -999.;
    float var_BTAG_MUON_ID_DOWN          = -999.;
    float var_BTAG_MUON_ISO_UP           = -999.;
    float var_BTAG_MUON_ISO_DOWN         = -999.;
    float var_BTAG_L1PREFIRING_UP        = -999.;
    float var_BTAG_L1PREFIRING_DOWN      = -999.;

    // kin reco
    float var_KIN_UP        = -999.;
    float var_KIN_DOWN      = -999.;
    float var_LOOSEKIN_UP   = -999.;
    float var_LOOSEKIN_DOWN = -999.;

    // ME scales
    float var_MERENSCALE_UP   = -999.;
    float var_MERENSCALE_DOWN = -999.;
    float var_MEFACSCALE_UP   = -999.;
    float var_MEFACSCALE_DOWN = -999.;
    float var_MESCALE_UP      = -999.;
    float var_MESCALE_DOWN    = -999.;
    float var_BFRAG_UP        = -999.;
    float var_BFRAG_DOWN      = -999.;
    float var_BFRAG_PETERSON  = -999.;
    float var_BSEMILEP_UP     = -999.;
    float var_BSEMILEP_DOWN   = -999.;
    float var_PDF_ALPHAS_UP   = -999.;
    float var_PDF_ALPHAS_DOWN = -999.;

    uint n_PDFweights;
    float var_PDF[length_];
    uint n_PSweights;
    float var_PS[length_];

    float var_PSSCALE_WEIGHT_4_UP   = -999.;
    float var_PSSCALE_WEIGHT_4_DOWN = -999.;
    float var_PSSCALE_WEIGHT_5_UP   = -999.;
    float var_PSSCALE_WEIGHT_5_DOWN = -999.;

};

class MiniTree// : public TTree
{

    void AddCommonBranches();
    void AddRecoBranches();
    void AddKinRecoBranches();
    void AddLooseKinRecoBranches();
    void AddGenBranches();
    void AddSystematicVariationWeightBranches();

    void generatorExtraJets(const TopGenObjects& topGenObjects,std::vector<int>& genExtraJetIndices, const int bHadronIndex, const int antiBHadronIndex, const bool withNu=false);

    void ResetBranches();

  public:
    TreeVars Vars;
    MiniTree(const TreeSetup& setup_);

    void InitWriteTree();
    void Fill(const AnalysisEntry* entry);
    // pointer to ROOT TTree
    TTree* TTreePtr = NULL;
    // generator level flag (true or false for MC, false for data)
    const bool isGen_;
    // identifier of the analysis selection step
    bool isTTbarSample_;
    bool isSystSample_;
    bool isSystVar_;


    /// Set the trigger scale factors
    void SetTriggerScaleFactors_Up(const TriggerScaleFactors* const scaleFactors);
    void SetTriggerScaleFactors_Dn(const TriggerScaleFactors* const scaleFactors);

    /// Set the btag scale factors
    void SetBtagScaleFactors_Var_AllSyst(const std::map<std::tuple<Systematic::Type, Systematic::Variation, int>, BtagScaleFactors*>* mapAllSystBTagScaleFactors);
    /// Set the kinematic reconstruction
    void SetKinematicReconstruction_Up(const KinematicReconstruction* const kinematicReconstruction,
                                    const KinematicReconstructionScaleFactors* const kinematicReconstructionScaleFactors);
    /// Set the kinematic reconstruction
    void SetKinematicReconstruction_Dn(const KinematicReconstruction* const kinematicReconstruction,
                                    const KinematicReconstructionScaleFactors* const kinematicReconstructionScaleFactors);


    /// Set the loose kinematic reconstruction
    void SetLooseKinReco_Up(const LooseKinReco* const looseKinReco,
                               const LooseKinRecoScaleFactors* const looseKinRecoScaleFactors);
    /// Set the loose kinematic reconstruction
    void SetLooseKinReco_Dn(const LooseKinReco* const looseKinReco,
                               const LooseKinRecoScaleFactors* const looseKinRecoScaleFactors);

protected:
    // add instances of reweighting classes for syst variations

    /// Pointer to the kinematic reconstruction instance
    const KinematicReconstruction* kinematicReconstruction_Up_;
    const KinematicReconstruction* kinematicReconstruction_Dn_;

    /// Pointer to the kinematic reconstruction scale factors instance
    const KinematicReconstructionScaleFactors* kinematicReconstructionScaleFactors_Up_;
    const KinematicReconstructionScaleFactors* kinematicReconstructionScaleFactors_Dn_;

    /// Pointer to the loose kinematic reconstruction instance
    const LooseKinReco* looseKinReco_Up_;
    const LooseKinReco* looseKinReco_Dn_;

    /// Pointer to the loose kinematic reconstruction scale factors instance
    const LooseKinRecoScaleFactors* looseKinRecoScaleFactors_Up_;
    const LooseKinRecoScaleFactors* looseKinRecoScaleFactors_Dn_;

    /// Pointer to the pileup reweighter instance
    const PileupScaleFactors* pileupScaleFactors_Up_;
    const PileupScaleFactors* pileupScaleFactors_Dn_;

    /// Pointer to lepton scale factors instance
    const LeptonScaleFactors* leptonScaleFactors_EleUp_;
    const LeptonScaleFactors* leptonScaleFactors_EleDn_;
    const LeptonScaleFactors* leptonScaleFactors_MuonUp_;
    const LeptonScaleFactors* leptonScaleFactors_MuonDn_;

    /// Pointer to trigger scale factors instance
    const TriggerScaleFactors* triggerScaleFactors_Up_;
    const TriggerScaleFactors* triggerScaleFactors_Dn_;

    const std::map<std::tuple<Systematic::Type, Systematic::Variation, int>, BtagScaleFactors*>* mapAllSystBTagScaleFactors_;

    /// Get weight due to trigger efficiency MC-to-data scale factors
    float weightTriggerSF_Up(const int leptonXIndex, const int leptonYIndex,
                             const VLV& allLeptons, const std::vector<int>& lepPdgId)const;
    float weightTriggerSF_Dn(const int leptonXIndex, const int leptonYIndex,
                             const VLV& allLeptons, const std::vector<int>& lepPdgId)const;

     /// Get weight due to b-tagging efficiency MC-to-data scale factors
     float weightBtagSF_Var(const std::vector<int>& jetIndices,
                         const VLV& jets, const std::vector<int>& jetFlavour,
                         const std::vector<float>& btagDiscriminators, const TString& sysKey)const;

     /// Get weight due to efficiency of kinematic reconstruction
     float weightKinReco_Up()const;
     float weightKinReco_Dn()const;

     /// Get weight due to efficiency of loose kinematic reconstruction
     float weightLooseKinReco_Up()const;
     float weightLooseKinReco_Dn()const;

};

#endif // MINITREE_H
