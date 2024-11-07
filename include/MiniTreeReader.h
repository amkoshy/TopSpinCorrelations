#ifndef MINITREEREADER_H
#define MINITREEREADER_H

#include <vector>
#include <TROOT.h>
#include <TTree.h>
#include <TBranch.h>

#include <sampleHelpers.h>

using namespace std;

const UInt_t length = 1000;

class MiniTreeReader {

    public:
    MiniTreeReader(TTree& tree, const Systematic::Systematic& sys);

    //common values
    double weight;
    float leptonSF;
    float l1PrefiringWeight;
    float triggerSF;
    float btagSF;
    float pileupSF;

    // float DNN_rho;

    //reco values
    bool  passStep3;
    // bool  passStep4;
    // bool  passStep5;
    // bool  passStep6;
    // bool  passStep7;
    // bool  passStep8;
    // bool  passStep8Loose;
    float lepton1_pt;
    float lepton1_eta;
    float lepton1_phi;
    float lepton1_m;
    float lepton2_pt;
    float lepton2_eta;
    float lepton2_phi;
    float lepton2_m;
    float met_pt;
    // float met_eta;
    float met_phi;
    // float met_m;
    float met_significance;
    uint  n_jets;
    float jets_pt[length];
    float jets_eta[length];
    float jets_phi[length];
    float jets_m[length];
    bool  jets_btag[length];
    // uint n_bjets;
    // float bjets_pt[length];
    // float bjets_eta[length];
    // float bjets_phi[length];
    // float bjets_m[length];


    // kinematic Reconstruction
    bool hasKinRecoSolution;
    float kinReco_top_pt;
    float kinReco_top_eta;
    float kinReco_top_phi;
    float kinReco_top_m;
    float kinReco_antitop_pt;
    float kinReco_antitop_eta;
    float kinReco_antitop_phi;
    float kinReco_antitop_m;
    float kinReco_lepton_pt;
    float kinReco_lepton_eta;
    float kinReco_lepton_phi;
    float kinReco_lepton_m;
    float kinReco_antilepton_pt;
    float kinReco_antilepton_eta;
    float kinReco_antilepton_phi;
    float kinReco_antilepton_m;
    float kinReco_b_pt;
    float kinReco_b_eta;
    float kinReco_b_phi;
    float kinReco_b_m;
    float kinReco_antib_pt;
    float kinReco_antib_eta;
    float kinReco_antib_phi;
    float kinReco_antib_m;
    // uint n_kinReco_nonbjets;
    // float kinReco_nonbjets_pt[length];
    // float kinReco_nonbjets_eta[length];
    // float kinReco_nonbjets_phi[length];
    // float kinReco_nonbjets_m[length];
    float kinReco_nonbjet_pt;
    float kinReco_nonbjet_eta;
    float kinReco_nonbjet_phi;
    float kinReco_nonbjet_m;
    float kinReco_rho;
    // float kinReco_weight;

    // loose kinematic reconstruction
    bool hasLooseKinRecoSolution;
    float looseKinReco_ttbar_pt;
    float looseKinReco_ttbar_eta;
    float looseKinReco_ttbar_phi;
    float looseKinReco_ttbar_m;
    float looseKinReco_b1_pt;
    float looseKinReco_b1_eta;
    float looseKinReco_b1_phi;
    float looseKinReco_b1_m;
    float looseKinReco_b2_pt;
    float looseKinReco_b2_eta;
    float looseKinReco_b2_phi;
    float looseKinReco_b2_m;
    float looseKinReco_l1_pt;
    float looseKinReco_l1_eta;
    float looseKinReco_l1_phi;
    float looseKinReco_l1_m;
    float looseKinReco_l2_pt;
    float looseKinReco_l2_eta;
    float looseKinReco_l2_phi;
    float looseKinReco_l2_m;
    float looseKinReco_weight;
    uint n_looseKinReco_nonbjets;
    float looseKinReco_nonbjets_pt[length];
    float looseKinReco_nonbjets_eta[length];
    float looseKinReco_nonbjets_phi[length];
    float looseKinReco_nonbjets_m[length];
    float looseKinReco_nonbjet_pt;
    float looseKinReco_nonbjet_eta;
    float looseKinReco_nonbjet_phi;
    float looseKinReco_nonbjet_m;
    float looseKinReco_rho;

    // systematic variations
    float systVarWeight_MERENSCALE_UP;
    float systVarWeight_MERENSCALE_DOWN;
    float systVarWeight_MEFACSCALE_UP;
    float systVarWeight_MEFACSCALE_DOWN;
    float systVarWeight_MESCALE_UP;
    float systVarWeight_MESCALE_DOWN;
    float systVarWeight_BFRAG_UP;
    float systVarWeight_BFRAG_DOWN;
    float systVarWeight_BFRAG_PETERSON;
    float systVarWeight_BSEMILEP_UP;
    float systVarWeight_BSEMILEP_DOWN;

    // float systVarWeight_PS[length];
    float systVarWeight_PSSCALE_WEIGHT_4_UP;
    float systVarWeight_PSSCALE_WEIGHT_4_DOWN;
    float systVarWeight_PSSCALE_WEIGHT_5_UP;
    float systVarWeight_PSSCALE_WEIGHT_5_DOWN;

    float systVarWeight_PDF_ALPHAS_UP;
    float systVarWeight_PDF_ALPHAS_DOWN;
    float systVarWeight_PU_UP;
    float systVarWeight_PU_DOWN;
    float systVarWeight_TRIG_UP;
    float systVarWeight_TRIG_DOWN;

    float systVarWeight_ELE_ID_UP;
    float systVarWeight_ELE_ID_DOWN;
    float systVarWeight_ELE_RECO_UP;
    float systVarWeight_ELE_RECO_DOWN;
    float systVarWeight_MUON_ID_UP;
    float systVarWeight_MUON_ID_DOWN;
    float systVarWeight_MUON_ISO_UP;
    float systVarWeight_MUON_ISO_DOWN;
    float systVarWeight_L1PREFIRING_UP;
    float systVarWeight_L1PREFIRING_DOWN;

    float systVarWeight_BTAG_UP;
    float systVarWeight_BTAG_DOWN;
    float systVarWeight_BTAG_LJET_UP;
    float systVarWeight_BTAG_LJET_DOWN;

    float systVarWeight_BTAG_MERENSCALE_UP;
    float systVarWeight_BTAG_MERENSCALE_DOWN;
    float systVarWeight_BTAG_MEFACSCALE_UP;
    float systVarWeight_BTAG_MEFACSCALE_DOWN;
    float systVarWeight_BTAG_MESCALE_UP;
    float systVarWeight_BTAG_MESCALE_DOWN;
    float systVarWeight_BTAG_BFRAG_UP;
    float systVarWeight_BTAG_BFRAG_DOWN;
    float systVarWeight_BTAG_BFRAG_PETERSON;
    float systVarWeight_BTAG_BSEMILEP_UP;
    float systVarWeight_BTAG_BSEMILEP_DOWN;
    float systVarWeight_BTAG_PDF_ALPHAS_UP;
    float systVarWeight_BTAG_PDF_ALPHAS_DOWN;
    float systVarWeight_BTAG_PSSCALE_WEIGHT_4_UP;
    float systVarWeight_BTAG_PSSCALE_WEIGHT_4_DOWN;
    float systVarWeight_BTAG_PSSCALE_WEIGHT_5_UP;
    float systVarWeight_BTAG_PSSCALE_WEIGHT_5_DOWN;
    float systVarWeight_BTAG_PU_UP;
    float systVarWeight_BTAG_PU_DOWN;
    float systVarWeight_BTAG_TRIG_UP;
    float systVarWeight_BTAG_TRIG_DOWN;
    float systVarWeight_BTAG_ELE_ID_UP;
    float systVarWeight_BTAG_ELE_ID_DOWN;
    float systVarWeight_BTAG_ELE_RECO_UP;
    float systVarWeight_BTAG_ELE_RECO_DOWN;
    float systVarWeight_BTAG_MUON_ID_UP;
    float systVarWeight_BTAG_MUON_ID_DOWN;
    float systVarWeight_BTAG_MUON_ISO_UP;
    float systVarWeight_BTAG_MUON_ISO_DOWN;
    float systVarWeight_BTAG_L1PREFIRING_UP;
    float systVarWeight_BTAG_L1PREFIRING_DOWN;


    // gen information
    // uint n_gen_additional_jets;
    // TTreePtr->Branch("gen_additional_jets_p[length];
    // TTreePtr->Branch("gen_additional_jets_eta[length];
    // TTreePtr->Branch("gen_additional_jets_phi[length];
    // TTreePtr->Branch("gen_additional_jets_m[length];

    // uint n_gen_additional_jets_withNu;
    // TTreePtr->Branch("gen_additional_jets_withNu_pt[length];
    // TTreePtr->Branch("gen_additional_jets_withNu_eta[length];
    // TTreePtr->Branch("gen_additional_jets_withNu_phi[length];
    // TTreePtr->Branch("gen_additional_jets_withNu_m[length];

    float gen_additional_jet_pt;
    float gen_additional_jet_eta;
    float gen_additional_jet_phi;
    float gen_additional_jet_m;

    float gen_additional_jet_withNu_pt;
    float gen_additional_jet_withNu_eta;
    float gen_additional_jet_withNu_phi;
    float gen_additional_jet_withNu_m;

    float gen_top_pt;
    float gen_top_eta;
    float gen_top_phi;
    float gen_top_m;
    float gen_antitop_pt;
    float gen_antitop_eta;
    float gen_antitop_phi;
    float gen_antitop_m;
    float gen_lepton_pt;
    float gen_lepton_eta;
    float gen_lepton_phi;
    float gen_lepton_m;
    float gen_antilepton_pt;
    float gen_antilepton_eta;
    float gen_antilepton_phi;
    float gen_antilepton_m;
    float gen_b_pt;
    float gen_b_eta;
    float gen_b_phi;
    float gen_b_m;
    float gen_antib_pt;
    float gen_antib_eta;
    float gen_antib_phi;
    float gen_antib_m;
    float gen_neutrino_pt;
    float gen_neutrino_eta;
    float gen_neutrino_phi;
    float gen_neutrino_m;
    float gen_antineutrino_pt;
    float gen_antineutrino_eta;
    float gen_antineutrino_phi;
    float gen_antineutrino_m;

    // uint gen_nVtx;
    // uint gen_nJets;
    // uint gen_nBJets;

    float gen_met_pt;
    float gen_met_eta;
    float gen_met_phi;
    float gen_met_m;

    float gen_rho;
    float gen_rhoWithNu;

    float pseudo_top_pt;
    float pseudo_top_eta;
    float pseudo_top_phi;
    float pseudo_top_m;
    float pseudo_antitop_pt;
    float pseudo_antitop_eta;
    float pseudo_antitop_phi;
    float pseudo_antitop_m;
    float pseudo_lepton_pt;
    float pseudo_lepton_eta;
    float pseudo_lepton_phi;
    float pseudo_lepton_m;
    float pseudo_antilepton_pt;
    float pseudo_antilepton_eta;
    float pseudo_antilepton_phi;
    float pseudo_antilepton_m;
    float pseudo_b_pt;
    float pseudo_b_eta;
    float pseudo_b_phi;
    float pseudo_b_m;
    float pseudo_antib_pt;
    float pseudo_antib_eta;
    float pseudo_antib_phi;
    float pseudo_antib_m;
    // uint n_pseudo_additional_jets;
    // TTreePtr->Branch("pseudo_additional_jets_pt[length];
    // TTreePtr->Branch("pseudo_additional_jets_eta[length];
    // TTreePtr->Branch("pseudo_additional_jets_phi[length];
    // TTreePtr->Branch("pseudo_additional_jets_m[length];
    float pseudo_additional_jet_pt;
    float pseudo_additional_jet_eta;
    float pseudo_additional_jet_phi;
    float pseudo_additional_jet_m;
    float pseudo_rho;

    float gen_partonLevel_additional_jet_pt;
    float gen_partonLevel_additional_jet_eta;
    float gen_partonLevel_additional_jet_phi;
    float gen_partonLevel_additional_jet_m;
    float gen_partonLevel_rho;
};


//     // std::unique_ptr<TTreeReaderValue<unsigned long long>> eventNumber;

#endif
