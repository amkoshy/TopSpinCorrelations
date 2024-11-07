#include "MiniTree.h"
#include "analysisStructs.h"
#include "analysisObjectStructs.h"
#include "KinematicReconstructionSolution.h"
#include "LooseKinRecoSolution.h"
#include "analysisUtils.h"
#include "AnalysisBase.h"
#include "utils.h"
#include <TChain.h>
#include <TFile.h>
#include "Math/VectorUtil.h"


MiniTree::MiniTree(const TreeSetup& setup_):
isGen_(setup_.isGen),
isTTbarSample_(setup_.isTTbarSample),
isSystSample_(setup_.isSystSample),
isSystVar_(setup_.isSystVariation)
{
  TTreePtr = new TTree(setup_.name, setup_.title);
}


void MiniTree::InitWriteTree()
{
  AddCommonBranches();

  AddRecoBranches();

  AddKinRecoBranches();

  AddLooseKinRecoBranches();

  AddGenBranches();

  AddSystematicVariationWeightBranches();

}
void MiniTree::AddCommonBranches()
{
    // TTreePtr->Branch("eventCounter",      &this->Vars.eventCounter, "eventCounter/I");
    TTreePtr->Branch("eventNumber",       &this->Vars.eventNumber, "eventNumber/l");
    // TTreePtr->Branch("lumiNumber",        &this->Vars.lumiNumber, "lumiNumber/i");
    // TTreePtr->Branch("runNumber",         &this->Vars.runNumber, "runNumber/i");
    TTreePtr->Branch("weight",            &this->Vars.weight, "weight/D");
    TTreePtr->Branch("leptonSF",          &this->Vars.leptonSF, "leptonSF/F");
    TTreePtr->Branch("l1PrefiringWeight", &this->Vars.l1PrefiringWeight, "l1PrefiringWeight/F");
    TTreePtr->Branch("triggerSF",         &this->Vars.triggerSF, "triggerSF/F");
    TTreePtr->Branch("btagSF",            &this->Vars.btagSF, "btagSF/F");
    TTreePtr->Branch("pileupSF",          &this->Vars.pileupSF, "pileupSF/F");
    // TTreePtr->Branch("channel",           &this->Vars.channel, "channel/I");
}



void MiniTree::AddSystematicVariationWeightBranches(){
    TTreePtr->Branch("var_MERENSCALE_UP",             &this->Vars.var_MERENSCALE_UP,   "var_MERENSCALE_UP/F");
    TTreePtr->Branch("var_MERENSCALE_DOWN",           &this->Vars.var_MERENSCALE_DOWN, "var_MERENSCALE_DOWN/F");
    TTreePtr->Branch("var_MEFACSCALE_UP",             &this->Vars.var_MEFACSCALE_UP,   "var_MEFACSCALE_UP/F");
    TTreePtr->Branch("var_MEFACSCALE_DOWN",           &this->Vars.var_MEFACSCALE_DOWN, "var_MEFACSCALE_DOWN/F");
    TTreePtr->Branch("var_MESCALE_UP",                &this->Vars.var_MESCALE_UP,      "var_MESCALE_UP/F");
    TTreePtr->Branch("var_MESCALE_DOWN",              &this->Vars.var_MESCALE_DOWN,    "var_MESCALE_DOWN/F");
    TTreePtr->Branch("var_BFRAG_UP",                  &this->Vars.var_BFRAG_UP,        "var_BFRAG_UP/F");
    TTreePtr->Branch("var_BFRAG_DOWN",                &this->Vars.var_BFRAG_DOWN,      "var_BFRAG_DOWN/F");
    TTreePtr->Branch("var_BFRAG_PETERSON",            &this->Vars.var_BFRAG_PETERSON,  "var_BFRAG_PETERSON/F");
    TTreePtr->Branch("var_BSEMILEP_UP",               &this->Vars.var_BSEMILEP_UP,     "var_BSEMILEP_UP/F");
    TTreePtr->Branch("var_BSEMILEP_DOWN",             &this->Vars.var_BSEMILEP_DOWN,   "var_BSEMILEP_DOWN/F");
    TTreePtr->Branch("var_PDF_ALPHAS_UP",             &this->Vars.var_PDF_ALPHAS_UP,   "var_PDF_ALPHAS_UP/F");
    TTreePtr->Branch("var_PDF_ALPHAS_DOWN",           &this->Vars.var_PDF_ALPHAS_DOWN, "var_PDF_ALPHAS_DOWN/F");

    // TTreePtr->Branch("n_PSweights", &this->Vars.n_PSweights, "n_PSweights/i");
    // TTreePtr->Branch("var_PS", &this->Vars.var_PS, "var_PS[n_PSweights]/F");
    // TTreePtr->Branch("n_PDFweights", &this->Vars.n_PDFweights, "n_PDFweights/i");
    // TTreePtr->Branch("var_PDF", &this->Vars.var_PDF, "var_PDF[n_PDFweights]/F");

    TTreePtr->Branch("var_PSSCALE_WEIGHT_4_UP",       &this->Vars.var_PSSCALE_WEIGHT_4_UP,   "var_PSSCALE_WEIGHT_4_UP/F");
    TTreePtr->Branch("var_PSSCALE_WEIGHT_4_DOWN",     &this->Vars.var_PSSCALE_WEIGHT_4_DOWN, "var_PSSCALE_WEIGHT_4_DOWN/F");
    TTreePtr->Branch("var_PSSCALE_WEIGHT_5_UP",       &this->Vars.var_PSSCALE_WEIGHT_5_UP,   "var_PSSCALE_WEIGHT_5_UP/F");
    TTreePtr->Branch("var_PSSCALE_WEIGHT_5_DOWN",     &this->Vars.var_PSSCALE_WEIGHT_5_DOWN, "var_PSSCALE_WEIGHT_5_DOWN/F");

    TTreePtr->Branch("var_PU_UP",                     &this->Vars.var_PU_UP,            "var_PU_UP/F");
    TTreePtr->Branch("var_PU_DOWN",                   &this->Vars.var_PU_DOWN,          "var_PU_DOWN/F");
    TTreePtr->Branch("var_TRIG_UP",                   &this->Vars.var_TRIG_UP,          "var_TRIG_UP/F");
    TTreePtr->Branch("var_TRIG_DOWN",                 &this->Vars.var_TRIG_DOWN,        "var_TRIG_DOWN/F");
    TTreePtr->Branch("var_ELE_ID_UP",                 &this->Vars.var_ELE_ID_UP,        "var_ELE_ID_UP/F");
    TTreePtr->Branch("var_ELE_ID_DOWN",               &this->Vars.var_ELE_ID_DOWN,      "var_ELE_ID_DOWN/F");
    TTreePtr->Branch("var_ELE_RECO_UP",               &this->Vars.var_ELE_RECO_UP,      "var_ELE_RECO_UP/F");
    TTreePtr->Branch("var_ELE_RECO_DOWN",             &this->Vars.var_ELE_RECO_DOWN,    "var_ELE_RECO_DOWN/F");
    TTreePtr->Branch("var_MUON_ID_UP",                &this->Vars.var_MUON_ID_UP,       "var_MUON_ID_UP/F");
    TTreePtr->Branch("var_MUON_ID_DOWN",              &this->Vars.var_MUON_ID_DOWN,     "var_MUON_ID_DOWN/F");
    TTreePtr->Branch("var_MUON_ISO_UP",               &this->Vars.var_MUON_ISO_UP,      "var_MUON_ISO_UP/F");
    TTreePtr->Branch("var_MUON_ISO_DOWN",             &this->Vars.var_MUON_ISO_DOWN,    "var_MUON_ISO_DOWN/F");
    TTreePtr->Branch("var_L1PREFIRING_UP",            &this->Vars.var_L1PREFIRING_UP,   "var_L1PREFIRING_UP/F");
    TTreePtr->Branch("var_L1PREFIRING_DOWN",          &this->Vars.var_L1PREFIRING_DOWN, "var_L1PREFIRING_DOWN/F");
    // TTreePtr->Branch("var_KIN_UP",                    &this->Vars.var_KIN_UP, "var_KIN_UP/F");
    // TTreePtr->Branch("var_KIN_DOWN",                  &this->Vars.var_KIN_DOWN, "var_KIN_DOWN/F");
    // TTreePtr->Branch("var_LOOSEKIN_UP",               &this->Vars.var_LOOSEKIN_UP, "var_LOOSEKIN_UP/F");
    // TTreePtr->Branch("var_LOOSEKIN_DOWN",             &this->Vars.var_LOOSEKIN_DOWN, "var_LOOSEKIN_DOWN/F");
    TTreePtr->Branch("var_BTAG_UP",                   &this->Vars.var_BTAG_UP,        "var_BTAG_UP/F");
    TTreePtr->Branch("var_BTAG_DOWN",                 &this->Vars.var_BTAG_DOWN,      "var_BTAG_DOWN/F");
    TTreePtr->Branch("var_BTAG_LJET_UP",              &this->Vars.var_BTAG_LJET_UP,   "var_BTAG_LJET_UP/F");
    TTreePtr->Branch("var_BTAG_LJET_DOWN",            &this->Vars.var_BTAG_LJET_DOWN, "var_BTAG_LJET_DOWN/F");

    TTreePtr->Branch("var_BTAG_MERENSCALE_UP",         &this->Vars.var_BTAG_MERENSCALE_UP,         "var_BTAG_MERENSCALE_UP/F");
    TTreePtr->Branch("var_BTAG_MERENSCALE_DOWN",       &this->Vars.var_BTAG_MERENSCALE_DOWN,       "var_BTAG_MERENSCALE_DOWN/F");
    TTreePtr->Branch("var_BTAG_MEFACSCALE_UP",         &this->Vars.var_BTAG_MEFACSCALE_UP,         "var_BTAG_MEFACSCALE_UP/F");
    TTreePtr->Branch("var_BTAG_MEFACSCALE_DOWN",       &this->Vars.var_BTAG_MEFACSCALE_DOWN,       "var_BTAG_MEFACSCALE_DOWN/F");
    TTreePtr->Branch("var_BTAG_MESCALE_UP",            &this->Vars.var_BTAG_MESCALE_UP,            "var_BTAG_MESCALE_UP/F");
    TTreePtr->Branch("var_BTAG_MESCALE_DOWN",          &this->Vars.var_BTAG_MESCALE_DOWN,          "var_BTAG_MESCALE_DOWN/F");
    TTreePtr->Branch("var_BTAG_BFRAG_UP",              &this->Vars.var_BTAG_BFRAG_UP,              "var_BTAG_BFRAG_UP/F");
    TTreePtr->Branch("var_BTAG_BFRAG_DOWN",            &this->Vars.var_BTAG_BFRAG_DOWN,            "var_BTAG_BFRAG_DOWN/F");
    TTreePtr->Branch("var_BTAG_BFRAG_PETERSON",        &this->Vars.var_BTAG_BFRAG_PETERSON,        "var_BTAG_BFRAG_PETERSON/F");
    TTreePtr->Branch("var_BTAG_BSEMILEP_UP",           &this->Vars.var_BTAG_BSEMILEP_UP,           "var_BTAG_BSEMILEP_UP/F");
    TTreePtr->Branch("var_BTAG_BSEMILEP_DOWN",         &this->Vars.var_BTAG_BSEMILEP_DOWN,         "var_BTAG_BSEMILEP_DOWN/F");
    TTreePtr->Branch("var_BTAG_PDF_ALPHAS_UP",         &this->Vars.var_BTAG_PDF_ALPHAS_UP,         "var_BTAG_PDF_ALPHAS_UP/F");
    TTreePtr->Branch("var_BTAG_PDF_ALPHAS_DOWN",       &this->Vars.var_BTAG_PDF_ALPHAS_DOWN,       "var_BTAG_PDF_ALPHAS_DOWN/F");
    TTreePtr->Branch("var_BTAG_PSSCALE_WEIGHT_4_UP",   &this->Vars.var_BTAG_PSSCALE_WEIGHT_4_UP,   "var_BTAG_PSSCALE_WEIGHT_4_UP/F");
    TTreePtr->Branch("var_BTAG_PSSCALE_WEIGHT_4_DOWN", &this->Vars.var_BTAG_PSSCALE_WEIGHT_4_DOWN, "var_BTAG_PSSCALE_WEIGHT_4_DOWN/F");
    TTreePtr->Branch("var_BTAG_PSSCALE_WEIGHT_5_UP",   &this->Vars.var_BTAG_PSSCALE_WEIGHT_5_UP,   "var_BTAG_PSSCALE_WEIGHT_5_UP/F");
    TTreePtr->Branch("var_BTAG_PSSCALE_WEIGHT_5_DOWN", &this->Vars.var_BTAG_PSSCALE_WEIGHT_5_DOWN, "var_BTAG_PSSCALE_WEIGHT_5_DOWN/F");
    TTreePtr->Branch("var_BTAG_PU_UP",                 &this->Vars.var_BTAG_PU_UP,                 "var_BTAG_PU_UP/F");
    TTreePtr->Branch("var_BTAG_PU_DOWN",               &this->Vars.var_BTAG_PU_DOWN,               "var_BTAG_PU_DOWN/F");
    TTreePtr->Branch("var_BTAG_TRIG_UP",               &this->Vars.var_BTAG_TRIG_UP,               "var_BTAG_TRIG_UP/F");
    TTreePtr->Branch("var_BTAG_TRIG_DOWN",             &this->Vars.var_BTAG_TRIG_DOWN,             "var_BTAG_TRIG_DOWN/F");
    TTreePtr->Branch("var_BTAG_ELE_ID_UP",             &this->Vars.var_BTAG_ELE_ID_UP,             "var_BTAG_ELE_ID_UP/F");
    TTreePtr->Branch("var_BTAG_ELE_ID_DOWN",           &this->Vars.var_BTAG_ELE_ID_DOWN,           "var_BTAG_ELE_ID_DOWN/F");
    TTreePtr->Branch("var_BTAG_ELE_RECO_UP",           &this->Vars.var_BTAG_ELE_RECO_UP,           "var_BTAG_ELE_RECO_UP/F");
    TTreePtr->Branch("var_BTAG_ELE_RECO_DOWN",         &this->Vars.var_BTAG_ELE_RECO_DOWN,         "var_BTAG_ELE_RECO_DOWN/F");
    TTreePtr->Branch("var_BTAG_MUON_ID_UP",            &this->Vars.var_BTAG_MUON_ID_UP,            "var_BTAG_MUON_ID_UP/F");
    TTreePtr->Branch("var_BTAG_MUON_ID_DOWN",          &this->Vars.var_BTAG_MUON_ID_DOWN,          "var_BTAG_MUON_ID_DOWN/F");
    TTreePtr->Branch("var_BTAG_MUON_ISO_UP",           &this->Vars.var_BTAG_MUON_ISO_UP,           "var_BTAG_MUON_ISO_UP/F");
    TTreePtr->Branch("var_BTAG_MUON_ISO_DOWN",         &this->Vars.var_BTAG_MUON_ISO_DOWN,         "var_BTAG_MUON_ISO_DOWN/F");
    TTreePtr->Branch("var_BTAG_L1PREFIRING_UP",        &this->Vars.var_BTAG_L1PREFIRING_UP,        "var_BTAG_L1PREFIRING_UP/F");
    TTreePtr->Branch("var_BTAG_L1PREFIRING_DOWN",      &this->Vars.var_BTAG_L1PREFIRING_DOWN,      "var_BTAG_L1PREFIRING_DOWN/F");




}



void MiniTree::AddRecoBranches()
{
    // TTreePtr->Branch("nVtx",     &this->Vars.nVtx,     "nVtx/i");
    // TTreePtr->Branch("nLeptons", &this->Vars.nLeptons, "nLeptons/i");
    // TTreePtr->Branch("nJets",    &this->Vars.nJets,    "nJets/i");
    // TTreePtr->Branch("nBJets",   &this->Vars.nBJets,   "nBJets/i");

    TTreePtr->Branch("lepton1_pt",       &this->Vars.lepton1_pt,  "lepton1_pt/F");
    TTreePtr->Branch("lepton1_eta",      &this->Vars.lepton1_eta, "lepton1_eta/F");
    TTreePtr->Branch("lepton1_phi",      &this->Vars.lepton1_phi, "lepton1_phi/F");
    TTreePtr->Branch("lepton1_m",        &this->Vars.lepton1_m,   "lepton1_m/F");
    // TTreePtr->Branch("lepton1_pdgID", &this->Vars.lepton1_pdgID, "lepton1_pdgID/F");

    TTreePtr->Branch("lepton2_pt",       &this->Vars.lepton2_pt,  "lepton2_pt/F");
    TTreePtr->Branch("lepton2_eta",      &this->Vars.lepton2_eta, "lepton2_eta/F");
    TTreePtr->Branch("lepton2_phi",      &this->Vars.lepton2_phi, "lepton2_phi/F");
    TTreePtr->Branch("lepton2_m",        &this->Vars.lepton2_m,   "lepton2_m/F");
    // TTreePtr->Branch("lepton2_pdgID", &this->Vars.lepton2_pdgID, "lepton2_pdgID/F");

    TTreePtr->Branch("met_pt",           &this->Vars.met_pt,           "met_pt/F");
    TTreePtr->Branch("met_phi",          &this->Vars.met_phi,          "met_phi/F");
    TTreePtr->Branch("met_significance", &this->Vars.met_significance, "met_significance/F");
    TTreePtr->Branch("n_jets",           &this->Vars.n_jets,           "n_jets/i");
    TTreePtr->Branch("jets_pt",          &this->Vars.jets_pt,          "jets_pt[n_jets]/F");
    TTreePtr->Branch("jets_eta",         &this->Vars.jets_eta,         "jets_eta[n_jets]/F");
    TTreePtr->Branch("jets_phi",         &this->Vars.jets_phi,         "jets_phi[n_jets]/F");
    TTreePtr->Branch("jets_m",           &this->Vars.jets_m,           "jets_m[n_jets]/F");
    TTreePtr->Branch("jets_charge",      &this->Vars.jets_charge,      "jets_charge[n_jets]/F");
    TTreePtr->Branch("jets_btag",        &this->Vars.jets_btag,        "jets_btag[n_jets]/O");
    // TTreePtr->Branch("jets_topMatched",    &this->Vars.jets_topMatched,    "jets_topMatched[n_jets]/O");
    // TTreePtr->Branch("jets_addJetMatched", &this->Vars.jets_addJetMatched, "jets_addJetMatched[n_jets]/O");

    // TTreePtr->Branch("n_bjets",          &this->Vars.n_bjets,      "n_bjets/i");
    // TTreePtr->Branch("bjets_pt",         &this->Vars.bjets_pt,     "bjets_pt[n_bjets]/F");
    // TTreePtr->Branch("bjets_eta",        &this->Vars.bjets_eta,    "bjets_eta[n_bjets]/F");
    // TTreePtr->Branch("bjets_phi",        &this->Vars.bjets_phi,    "bjets_phi[n_bjets]/F");
    // TTreePtr->Branch("bjets_m",          &this->Vars.bjets_m,      "bjets_m[n_bjets]/F");
    // TTreePtr->Branch("bjets_charge",     &this->Vars.bjets_charge, "bjets_charge[n_bjets]/F");

    // TTreePtr->Branch("n_nonbjets",      &this->Vars.n_nonbjets,      "n_nonbjets/i");
    // TTreePtr->Branch("nonbjets_pt",     &this->Vars.nonbjets_pt,     "nonbjets_pt[n_nonbjets]/F");
    // TTreePtr->Branch("nonbjets_eta",    &this->Vars.nonbjets_eta,    "nonbjets_eta[n_nonbjets]/F");
    // TTreePtr->Branch("nonbjets_phi",    &this->Vars.nonbjets_phi,    "nonbjets_phi[n_nonbjets]/F");
    // TTreePtr->Branch("nonbjets_m",      &this->Vars.nonbjets_m,      "nonbjets_m[n_nonbjets]/F");
    // TTreePtr->Branch("nonbjets_charge", &this->Vars.nonbjets_charge, "nonbjets_charge[n_nonbjets]/F");

    TTreePtr->Branch("passStep3",         &this->Vars.passStep3,      "passStep3/O");
    TTreePtr->Branch("passStep4",         &this->Vars.passStep4,      "passStep4/O");
    TTreePtr->Branch("passStep5",         &this->Vars.passStep5,      "passStep5/O");
    TTreePtr->Branch("passStep6",         &this->Vars.passStep6,      "passStep6/O");
    TTreePtr->Branch("passStep7",         &this->Vars.passStep7,      "passStep7/O");
    TTreePtr->Branch("passStep8",         &this->Vars.passStep8,      "passStep8/O");
    TTreePtr->Branch("passStep8Loose",    &this->Vars.passStep8Loose, "passStep8Loose/O");
}

void MiniTree::AddKinRecoBranches()
{
    TTreePtr->Branch("hasKinRecoSolution",     &this->Vars.hasKinRecoSolution,  "hasKinRecoSolution/O");

    TTreePtr->Branch("kinReco_top_pt",         &this->Vars.kinReco_top_pt,      "kinReco_top_pt/F");
    TTreePtr->Branch("kinReco_top_eta",        &this->Vars.kinReco_top_eta,     "kinReco_top_eta/F");
    TTreePtr->Branch("kinReco_top_phi",        &this->Vars.kinReco_top_phi,     "kinReco_top_phi/F");
    TTreePtr->Branch("kinReco_top_m",          &this->Vars.kinReco_top_m,       "kinReco_top_m/F");
    TTreePtr->Branch("kinReco_antitop_pt",     &this->Vars.kinReco_antitop_pt,  "kinReco_antitop_pt/F");
    TTreePtr->Branch("kinReco_antitop_eta",    &this->Vars.kinReco_antitop_eta, "kinReco_antitop_eta/F");
    TTreePtr->Branch("kinReco_antitop_phi",    &this->Vars.kinReco_antitop_phi, "kinReco_antitop_phi/F");
    TTreePtr->Branch("kinReco_antitop_m",      &this->Vars.kinReco_antitop_m,   "kinReco_antitop_m/F");

    TTreePtr->Branch("kinReco_lepton_pt",      &this->Vars.kinReco_lepton_pt,      "kinReco_lepton_pt/F");
    TTreePtr->Branch("kinReco_lepton_eta",     &this->Vars.kinReco_lepton_eta,     "kinReco_lepton_eta/F");
    TTreePtr->Branch("kinReco_lepton_phi",     &this->Vars.kinReco_lepton_phi,     "kinReco_lepton_phi/F");
    TTreePtr->Branch("kinReco_lepton_m",       &this->Vars.kinReco_lepton_m,       "kinReco_lepton_m/F");
    TTreePtr->Branch("kinReco_antilepton_pt",  &this->Vars.kinReco_antilepton_pt,  "kinReco_antilepton_pt/F");
    TTreePtr->Branch("kinReco_antilepton_eta", &this->Vars.kinReco_antilepton_eta, "kinReco_antilepton_eta/F");
    TTreePtr->Branch("kinReco_antilepton_phi", &this->Vars.kinReco_antilepton_phi, "kinReco_antilepton_phi/F");
    TTreePtr->Branch("kinReco_antilepton_m",   &this->Vars.kinReco_antilepton_m,   "kinReco_antilepton_m/F");

    TTreePtr->Branch("kinReco_b_pt",           &this->Vars.kinReco_b_pt,      "kinReco_b_pt/F");
    TTreePtr->Branch("kinReco_b_eta",          &this->Vars.kinReco_b_eta,     "kinReco_b_eta/F");
    TTreePtr->Branch("kinReco_b_phi",          &this->Vars.kinReco_b_phi,     "kinReco_b_phi/F");
    TTreePtr->Branch("kinReco_b_m",            &this->Vars.kinReco_b_m,       "kinReco_b_m/F");
    TTreePtr->Branch("kinReco_antib_pt",       &this->Vars.kinReco_antib_pt,  "kinReco_antib_pt/F");
    TTreePtr->Branch("kinReco_antib_eta",      &this->Vars.kinReco_antib_eta, "kinReco_antib_eta/F");
    TTreePtr->Branch("kinReco_antib_phi",      &this->Vars.kinReco_antib_phi, "kinReco_antib_phi/F");
    TTreePtr->Branch("kinReco_antib_m",        &this->Vars.kinReco_antib_m,   "kinReco_antib_m/F");

    // TTreePtr->Branch("n_kinReco_nonbjets",     &this->Vars.n_kinReco_nonbjets,   "n_kinReco_nonbjets/i");
    // TTreePtr->Branch("kinReco_nonbjets_pt",    &this->Vars.kinReco_nonbjets_pt,  "kinReco_nonbjets_pt[n_kinReco_nonbjets]/F");
    // TTreePtr->Branch("kinReco_nonbjets_eta",   &this->Vars.kinReco_nonbjets_eta, "kinReco_nonbjets_eta[n_kinReco_nonbjets]/F");
    // TTreePtr->Branch("kinReco_nonbjets_phi",   &this->Vars.kinReco_nonbjets_phi, "kinReco_nonbjets_ph[n_kinReco_nonbjets]i/F");
    // TTreePtr->Branch("kinReco_nonbjets_m",     &this->Vars.kinReco_nonbjets_m,   "kinReco_nonbjets_m[n_kinReco_nonbjets]/F");

    TTreePtr->Branch("kinReco_nonbjet_pt",    &this->Vars.kinReco_nonbjet_pt,  "kinReco_nonbjet_pt/F");
    TTreePtr->Branch("kinReco_nonbjet_eta",   &this->Vars.kinReco_nonbjet_eta, "kinReco_nonbjet_eta/F");
    TTreePtr->Branch("kinReco_nonbjet_phi",   &this->Vars.kinReco_nonbjet_phi, "kinReco_nonbjet_phi/F");
    TTreePtr->Branch("kinReco_nonbjet_m",     &this->Vars.kinReco_nonbjet_m,   "kinReco_nonbjet_m/F");

    TTreePtr->Branch("kinReco_rho",   &this->Vars.kinReco_rho, "kinReco_rho/F");

    // TTreePtr->Branch("kinReco_weight",         &this->Vars.kinReco_weight, "kinReco_weight/F");
}

void MiniTree::AddLooseKinRecoBranches()
{
    TTreePtr->Branch("hasLooseKinRecoSolution",   &this->Vars.hasLooseKinRecoSolution, "hasLooseKinRecoSolution/O");

    TTreePtr->Branch("looseKinReco_ttbar_pt",     &this->Vars.looseKinReco_ttbar_pt,   "looseKinReco_ttbar_pt/F");
    TTreePtr->Branch("looseKinReco_ttbar_eta",    &this->Vars.looseKinReco_ttbar_eta,  "looseKinReco_ttbar_eta/F");
    TTreePtr->Branch("looseKinReco_ttbar_phi",    &this->Vars.looseKinReco_ttbar_phi,  "looseKinReco_ttbar_phi/F");
    TTreePtr->Branch("looseKinReco_ttbar_m",      &this->Vars.looseKinReco_ttbar_m,    "looseKinReco_ttbar_m/F");

    TTreePtr->Branch("looseKinReco_b1_pt",        &this->Vars.looseKinReco_b1_pt,      "looseKinReco_b1_pt/F");
    TTreePtr->Branch("looseKinReco_b1_eta",       &this->Vars.looseKinReco_b1_eta,     "looseKinReco_b1_eta/F");
    TTreePtr->Branch("looseKinReco_b1_phi",       &this->Vars.looseKinReco_b1_phi,     "looseKinReco_b1_phi/F");
    TTreePtr->Branch("looseKinReco_b1_m",         &this->Vars.looseKinReco_b1_m,       "looseKinReco_b1_m/F");

    TTreePtr->Branch("looseKinReco_b2_pt",        &this->Vars.looseKinReco_b2_pt,      "looseKinReco_b2_pt/F");
    TTreePtr->Branch("looseKinReco_b2_eta",       &this->Vars.looseKinReco_b2_eta,     "looseKinReco_b2_eta/F");
    TTreePtr->Branch("looseKinReco_b2_phi",       &this->Vars.looseKinReco_b2_phi,     "looseKinReco_b2_phi/F");
    TTreePtr->Branch("looseKinReco_b2_m",         &this->Vars.looseKinReco_b2_m,       "looseKinReco_b2_m/F");

    TTreePtr->Branch("looseKinReco_l1_pt",        &this->Vars.looseKinReco_l1_pt,      "looseKinReco_l1_pt/F");
    TTreePtr->Branch("looseKinReco_l1_eta",       &this->Vars.looseKinReco_l1_eta,     "looseKinReco_l1_eta/F");
    TTreePtr->Branch("looseKinReco_l1_phi",       &this->Vars.looseKinReco_l1_phi,     "looseKinReco_l1_phi/F");
    TTreePtr->Branch("looseKinReco_l1_m",         &this->Vars.looseKinReco_l1_m,       "looseKinReco_l1_m/F");

    TTreePtr->Branch("looseKinReco_l2_pt",        &this->Vars.looseKinReco_l2_pt,      "looseKinReco_l2_pt/F");
    TTreePtr->Branch("looseKinReco_l2_eta",       &this->Vars.looseKinReco_l2_eta,     "looseKinReco_l2_eta/F");
    TTreePtr->Branch("looseKinReco_l2_phi",       &this->Vars.looseKinReco_l2_phi,     "looseKinReco_l2_phi/F");
    TTreePtr->Branch("looseKinReco_l2_m",         &this->Vars.looseKinReco_l2_m,       "looseKinReco_l2_m/F");

    // TTreePtr->Branch("looseKinReco_weight",       &this->Vars.looseKinReco_weight, "looseKinReco_weight/F");

    // TTreePtr->Branch("n_looseKinReco_nonbjets",   &this->Vars.n_looseKinReco_nonbjets,   "n_looseKinReco_nonbjets/i");
    // TTreePtr->Branch("looseKinReco_nonbjets_pt",  &this->Vars.looseKinReco_nonbjets_pt,  "looseKinReco_nonbjets_pt[n_looseKinReco_nonbjets]/F");
    // TTreePtr->Branch("looseKinReco_nonbjets_eta", &this->Vars.looseKinReco_nonbjets_eta, "looseKinReco_nonbjets_eta[n_looseKinReco_nonbjets]/F");
    // TTreePtr->Branch("looseKinReco_nonbjets_phi", &this->Vars.looseKinReco_nonbjets_phi, "looseKinReco_nonbjets_phi[n_looseKinReco_nonbjets]/F");
    // TTreePtr->Branch("looseKinReco_nonbjets_m",   &this->Vars.looseKinReco_nonbjets_m,   "looseKinReco_nonbjets_m[n_looseKinReco_nonbjets]/F");

    TTreePtr->Branch("looseKinReco_nonbjet_pt",  &this->Vars.looseKinReco_nonbjet_pt,  "looseKinReco_nonbjet_pt/F");
    TTreePtr->Branch("looseKinReco_nonbjet_eta", &this->Vars.looseKinReco_nonbjet_eta, "looseKinReco_nonbjet_eta/F");
    TTreePtr->Branch("looseKinReco_nonbjet_phi", &this->Vars.looseKinReco_nonbjet_phi, "looseKinReco_nonbjet_phi/F");
    TTreePtr->Branch("looseKinReco_nonbjet_m",   &this->Vars.looseKinReco_nonbjet_m,   "looseKinReco_nonbjet_m/F");

    TTreePtr->Branch("looseKinReco_rho",         &this->Vars.looseKinReco_rho,         "looseKinReco_ttbar_pt/F");
}

void MiniTree::AddGenBranches()
{
    // TTreePtr->Branch("n_gen_additional_jets",   &this->Vars.n_gen_additional_jets,   "n_gen_additional_jets/i");
    // TTreePtr->Branch("gen_additional_jets_pt",  &this->Vars.gen_additional_jets_pt,  "gen_additional_jets_pt[n_looseKinReco_nonbjets/F");
    // TTreePtr->Branch("gen_additional_jets_eta", &this->Vars.gen_additional_jets_eta, "gen_additional_jets_eta[n_looseKinReco_nonbjets/F");
    // TTreePtr->Branch("gen_additional_jets_phi", &this->Vars.gen_additional_jets_phi, "gen_additional_jets_phi[n_looseKinReco_nonbjets/F");
    // TTreePtr->Branch("gen_additional_jets_m",   &this->Vars.gen_additional_jets_m,   "gen_additional_jets_m[n_looseKinReco_nonbjets/F");

    // TTreePtr->Branch("n_gen_additional_jets_withNu",   &this->Vars.n_gen_additional_jets_withNu,   "n_gen_additional_jets_withNu/i");
    // TTreePtr->Branch("gen_additional_jets_withNu_pt",  &this->Vars.gen_additional_jets_withNu_pt,  "gen_additional_jets_withNu_pt[n_gen_additional_jets_withNu/F");
    // TTreePtr->Branch("gen_additional_jets_withNu_eta", &this->Vars.gen_additional_jets_withNu_eta, "gen_additional_jets_withNu_eta[n_gen_additional_jets_withNu/F");
    // TTreePtr->Branch("gen_additional_jets_withNu_phi", &this->Vars.gen_additional_jets_withNu_phi, "gen_additional_jets_withNu_phi[n_gen_additional_jets_withNu/F");
    // TTreePtr->Branch("gen_additional_jets_withNu_m",   &this->Vars.gen_additional_jets_withNu_m,   "gen_additional_jets_withNu_m[n_gen_additional_jets_withNu/F");

    TTreePtr->Branch("gen_additional_jet_pt",         &this->Vars.gen_additional_jet_pt,         "gen_additional_jet_pt/F");
    TTreePtr->Branch("gen_additional_jet_eta",        &this->Vars.gen_additional_jet_eta,        "gen_additional_jet_eta/F");
    TTreePtr->Branch("gen_additional_jet_phi",        &this->Vars.gen_additional_jet_phi,        "gen_additional_jet_phi/F");
    TTreePtr->Branch("gen_additional_jet_m",          &this->Vars.gen_additional_jet_m,          "gen_additional_jet_m/F");

    TTreePtr->Branch("gen_additional_jet_withNu_pt",  &this->Vars.gen_additional_jet_withNu_pt,  "gen_additional_jet_withNu_pt/F");
    TTreePtr->Branch("gen_additional_jet_withNu_eta", &this->Vars.gen_additional_jet_withNu_eta, "gen_additional_jet_withNu_eta/F");
    TTreePtr->Branch("gen_additional_jet_withNu_phi", &this->Vars.gen_additional_jet_withNu_phi, "gen_additional_jet_withNu_phi/F");
    TTreePtr->Branch("gen_additional_jet_withNu_m",   &this->Vars.gen_additional_jet_withNu_m,   "gen_additional_jet_withNu_m/F");

    TTreePtr->Branch("gen_top_pt",                    &this->Vars.gen_top_pt,                    "gen_top_pt/F");
    TTreePtr->Branch("gen_top_eta",                   &this->Vars.gen_top_eta,                   "gen_top_eta/F");
    TTreePtr->Branch("gen_top_phi",                   &this->Vars.gen_top_phi,                   "gen_top_phi/F");
    TTreePtr->Branch("gen_top_m",                     &this->Vars.gen_top_m,                     "gen_top_m/F");

    TTreePtr->Branch("gen_antitop_pt",                &this->Vars.gen_antitop_pt,                "gen_antitop_pt/F");
    TTreePtr->Branch("gen_antitop_eta",               &this->Vars.gen_antitop_eta,               "gen_antitop_eta/F");
    TTreePtr->Branch("gen_antitop_phi",               &this->Vars.gen_antitop_phi,               "gen_antitop_phi/F");
    TTreePtr->Branch("gen_antitop_m",                 &this->Vars.gen_antitop_m,                 "gen_antitop_m/F");

    TTreePtr->Branch("gen_lepton_pt",                 &this->Vars.gen_lepton_pt,                 "gen_lepton_pt/F");
    TTreePtr->Branch("gen_lepton_eta",                &this->Vars.gen_lepton_eta,                "gen_lepton_eta/F");
    TTreePtr->Branch("gen_lepton_phi",                &this->Vars.gen_lepton_phi,                "gen_lepton_phi/F");
    TTreePtr->Branch("gen_lepton_m",                  &this->Vars.gen_lepton_m,                  "gen_lepton_m/F");

    TTreePtr->Branch("gen_antilepton_pt",             &this->Vars.gen_antilepton_pt,             "gen_antilepton_pt/F");
    TTreePtr->Branch("gen_antilepton_eta",            &this->Vars.gen_antilepton_eta,            "gen_antilepton_eta/F");
    TTreePtr->Branch("gen_antilepton_phi",            &this->Vars.gen_antilepton_phi,            "gen_antilepton_phi/F");
    TTreePtr->Branch("gen_antilepton_m",              &this->Vars.gen_antilepton_m,              "gen_antilepton_m/F");

    TTreePtr->Branch("gen_b_pt",                      &this->Vars.gen_b_pt,                      "gen_b_pt/F");
    TTreePtr->Branch("gen_b_eta",                     &this->Vars.gen_b_eta,                     "gen_b_eta/F");
    TTreePtr->Branch("gen_b_phi",                     &this->Vars.gen_b_phi,                     "gen_b_phi/F");
    TTreePtr->Branch("gen_b_m",                       &this->Vars.gen_b_m,                       "gen_b_m/F");

    TTreePtr->Branch("gen_antib_pt",                  &this->Vars.gen_antib_pt,                  "gen_antib_pt/F");
    TTreePtr->Branch("gen_antib_eta",                 &this->Vars.gen_antib_eta,                 "gen_antib_eta/F");
    TTreePtr->Branch("gen_antib_phi",                 &this->Vars.gen_antib_phi,                 "gen_antib_phi/F");
    TTreePtr->Branch("gen_antib_m",                   &this->Vars.gen_antib_m,                   "gen_antib_m/F");

    TTreePtr->Branch("gen_neutrino_pt",               &this->Vars.gen_neutrino_pt,               "gen_neutrino_pt/F");
    TTreePtr->Branch("gen_neutrino_eta",              &this->Vars.gen_neutrino_eta,              "gen_neutrino_eta/F");
    TTreePtr->Branch("gen_neutrino_phi",              &this->Vars.gen_neutrino_phi,              "gen_neutrino_phi/F");
    TTreePtr->Branch("gen_neutrino_m",                &this->Vars.gen_neutrino_m,                "gen_neutrino_m/F");

    TTreePtr->Branch("gen_antineutrino_pt",           &this->Vars.gen_antineutrino_pt,           "gen_antineutrino_pt/F");
    TTreePtr->Branch("gen_antineutrino_eta",          &this->Vars.gen_antineutrino_eta,          "gen_antineutrino_eta/F");
    TTreePtr->Branch("gen_antineutrino_phi",          &this->Vars.gen_antineutrino_phi,          "gen_antineutrino_phi/F");
    TTreePtr->Branch("gen_antineutrino_m",            &this->Vars.gen_antineutrino_m,            "gen_antineutrino_m/F");

    // TTreePtr->Branch("gen_nVtx",               &this->Vars.gen_nVtx,    "gen_nVtx/F");
    // TTreePtr->Branch("gen_nJets",              &this->Vars.gen_nJets,   "gen_nJets/F");
    // TTreePtr->Branch("gen_nBJets",             &this->Vars.gen_nBJets,  "gen_nBJets/F");

    TTreePtr->Branch("gen_met_pt",            &this->Vars.gen_met_pt,  "gen_met_pt/F");
    TTreePtr->Branch("gen_met_eta",           &this->Vars.gen_met_eta, "gen_met_eta/F");
    TTreePtr->Branch("gen_met_phi",           &this->Vars.gen_met_phi, "gen_met_phi/F");
    TTreePtr->Branch("gen_met_m",             &this->Vars.gen_met_m,   "gen_met_m/F");

    TTreePtr->Branch("gen_rho",               &this->Vars.gen_rho,       "gen_rho/F");
    TTreePtr->Branch("gen_rhoWithNu",         &this->Vars.gen_rhoWithNu, "gen_rhoWithNu/F");

    TTreePtr->Branch("pseudo_top_pt",         &this->Vars.pseudo_top_pt,         "pseudo_top_pt/F");
    TTreePtr->Branch("pseudo_top_eta",        &this->Vars.pseudo_top_eta,        "pseudo_top_eta/F");
    TTreePtr->Branch("pseudo_top_phi",        &this->Vars.pseudo_top_phi,        "pseudo_top_phi/F");
    TTreePtr->Branch("pseudo_top_m",          &this->Vars.pseudo_top_m,          "pseudo_top_m/F");
    TTreePtr->Branch("pseudo_antitop_pt",     &this->Vars.pseudo_antitop_pt,     "pseudo_antitop_pt/F");
    TTreePtr->Branch("pseudo_antitop_eta",    &this->Vars.pseudo_antitop_eta,    "pseudo_antitop_eta/F");
    TTreePtr->Branch("pseudo_antitop_phi",    &this->Vars.pseudo_antitop_phi,    "pseudo_antitop_phi/F");
    TTreePtr->Branch("pseudo_antitop_m",      &this->Vars.pseudo_antitop_m,      "pseudo_antitop_m/F");
    TTreePtr->Branch("pseudo_lepton_pt",      &this->Vars.pseudo_lepton_pt,      "pseudo_lepton_pt/F");
    TTreePtr->Branch("pseudo_lepton_eta",     &this->Vars.pseudo_lepton_eta,     "pseudo_lepton_eta/F");
    TTreePtr->Branch("pseudo_lepton_phi",     &this->Vars.pseudo_lepton_phi,     "pseudo_lepton_phi/F");
    TTreePtr->Branch("pseudo_lepton_m",       &this->Vars.pseudo_lepton_m,       "pseudo_lepton_m/F");
    TTreePtr->Branch("pseudo_antilepton_pt",  &this->Vars.pseudo_antilepton_pt,  "pseudo_antilepton_pt/F");
    TTreePtr->Branch("pseudo_antilepton_eta", &this->Vars.pseudo_antilepton_eta, "pseudo_antilepton_eta/F");
    TTreePtr->Branch("pseudo_antilepton_phi", &this->Vars.pseudo_antilepton_phi, "pseudo_antilepton_phi/F");
    TTreePtr->Branch("pseudo_antilepton_m",   &this->Vars.pseudo_antilepton_m,   "pseudo_antilepton_m/F");
    TTreePtr->Branch("pseudo_b_pt",           &this->Vars.pseudo_b_pt,           "pseudo_b_pt/F");
    TTreePtr->Branch("pseudo_b_eta",          &this->Vars.pseudo_b_eta,          "pseudo_b_eta/F");
    TTreePtr->Branch("pseudo_b_phi",          &this->Vars.pseudo_b_phi,          "pseudo_b_phi/F");
    TTreePtr->Branch("pseudo_b_m",            &this->Vars.pseudo_b_m,            "pseudo_b_m/F");
    TTreePtr->Branch("pseudo_antib_pt",       &this->Vars.pseudo_antib_pt,       "pseudo_antib_pt/F");
    TTreePtr->Branch("pseudo_antib_eta",      &this->Vars.pseudo_antib_eta,      "pseudo_antib_eta/F");
    TTreePtr->Branch("pseudo_antib_phi",      &this->Vars.pseudo_antib_phi,      "pseudo_antib_phi/F");
    TTreePtr->Branch("pseudo_antib_m",        &this->Vars.pseudo_antib_m,        "pseudo_antib_m/F");

    // TTreePtr->Branch("n_pseudo_additional_jets",   &this->Vars.n_pseudo_additional_jets,   "n_pseudo_additional_jets/i");
    // TTreePtr->Branch("pseudo_additional_jets_pt",  &this->Vars.pseudo_additional_jets_pt,  "pseudo_additional_jets_pt[n_pseudo_additional_jets]/F");
    // TTreePtr->Branch("pseudo_additional_jets_eta", &this->Vars.pseudo_additional_jets_eta, "pseudo_additional_jets_eta[n_pseudo_additional_jets]/F");
    // TTreePtr->Branch("pseudo_additional_jets_phi", &this->Vars.pseudo_additional_jets_phi, "pseudo_additional_jets_phi[n_pseudo_additional_jets]/F");
    // TTreePtr->Branch("pseudo_additional_jets_m",   &this->Vars.pseudo_additional_jets_m,   "pseudo_additional_jets_m[n_pseudo_additional_jets]/F");

    TTreePtr->Branch("pseudo_additional_jet_pt",           &this->Vars.pseudo_additional_jet_pt,  "pseudo_additional_jet_pt/F");
    TTreePtr->Branch("pseudo_additional_jet_eta",          &this->Vars.pseudo_additional_jet_eta, "pseudo_additional_jet_eta/F");
    TTreePtr->Branch("pseudo_additional_jet_phi",          &this->Vars.pseudo_additional_jet_phi, "pseudo_additional_jet_phi/F");
    TTreePtr->Branch("pseudo_additional_jet_m",            &this->Vars.pseudo_additional_jet_m,   "pseudo_additional_jet_m/F");

    TTreePtr->Branch("pseudo_rho",                         &this->Vars.pseudo_rho, "pseudo_rho/F");

    TTreePtr->Branch("gen_partonLevel_additional_jet_pt",  &this->Vars.gen_partonLevel_additional_jet_pt,  "gen_partonLevel_additional_jet_pt/F");
    TTreePtr->Branch("gen_partonLevel_additional_jet_eta", &this->Vars.gen_partonLevel_additional_jet_eta, "gen_partonLevel_additional_jet_eta/F");
    TTreePtr->Branch("gen_partonLevel_additional_jet_phi", &this->Vars.gen_partonLevel_additional_jet_phi, "gen_partonLevel_additional_jet_phi/F");
    TTreePtr->Branch("gen_partonLevel_additional_jet_m",   &this->Vars.gen_partonLevel_additional_jet_m,   "gen_partonLevel_additional_jet_m/F");

    TTreePtr->Branch("gen_partonLevel_rho",                &this->Vars.gen_partonLevel_rho, "gen_partonLevel_rho/F");

}

void MiniTree::Fill(const AnalysisEntry* entry)
{
    const AnalysisConfig::Selections& selections =  entry->analysisConfig->selections();

    ResetBranches();

    const VLV& allGenJets =  entry->topGenObjects->valuesSet_ ? *entry->topGenObjects->allGenJets_ : VLV();
    std::vector<int> allGenJetIndices = common::initialiseIndices(allGenJets);
    std::vector<int> genJetIndices = allGenJetIndices;
    common::selectIndices(genJetIndices, allGenJets, common::LVeta, selections.genJetEtaCut_, false);
    common::selectIndices(genJetIndices, allGenJets, common::LVeta, -selections.genJetEtaCut_);
    common::selectIndices(genJetIndices, allGenJets, common::LVpt, selections.genJetPtCut_);
    if(selections.genDeltaRLeptonJetCut_ > 0.){
        // Vector of genLeptons from which genJets need to be separated in deltaR
        VLV allGenLeptons;
        if(entry->topGenObjects->valuesSet_){
            if(entry->topGenObjects->GenLepton_) allGenLeptons.push_back(*entry->topGenObjects->GenLepton_);
            if(entry->topGenObjects->GenAntiLepton_) allGenLeptons.push_back(*entry->topGenObjects->GenAntiLepton_);
        }
        entry->analysisBase->leptonCleanedJetIndices(genJetIndices, allGenJets, allGenLeptons, selections.genDeltaRLeptonJetCut_);
    }
    common::orderIndices(genJetIndices, allGenJets, common::LVpt);

    const VLV& allGenJetsWithNu =  entry->topGenObjects->valuesSet_ ? *entry->topGenObjects->allGenJetsWithNu_ : VLV();
    std::vector<int> allGenJetWithNuIndices = common::initialiseIndices(allGenJetsWithNu);
    std::vector<int> genJetWithNuIndices = allGenJetWithNuIndices;
    common::selectIndices(genJetWithNuIndices, allGenJetsWithNu, common::LVeta, selections.genJetEtaCut_, false);
    common::selectIndices(genJetWithNuIndices, allGenJetsWithNu, common::LVeta, -selections.genJetEtaCut_);
    common::selectIndices(genJetWithNuIndices, allGenJetsWithNu, common::LVpt, selections.genJetPtCut_);
    if(selections.genDeltaRLeptonJetCut_ > 0.){
        // Vector of genLeptons from which genJets need to be separated in deltaR
        VLV allGenLeptons;
        if(entry->topGenObjects->valuesSet_){
            if(entry->topGenObjects->GenLepton_) allGenLeptons.push_back(*entry->topGenObjects->GenLepton_);
            if(entry->topGenObjects->GenAntiLepton_) allGenLeptons.push_back(*entry->topGenObjects->GenAntiLepton_);
        }
        entry->analysisBase->leptonCleanedJetIndices(genJetWithNuIndices, allGenJetsWithNu, allGenLeptons, selections.genDeltaRLeptonJetCut_);
    }
    common::orderIndices(genJetWithNuIndices, allGenJetsWithNu, common::LVpt);


    //=======================PSEUDO-TOP: BEGIN==========================//
    if(isTTbarSample_){
        if(entry->topPseudoObjects->valuesSet_){
            // Access Top signal pseudo info
            // const TopPseudoObjects& topPseudoObjects = this->getTopPseudoObjects(entry);
            const TopPseudoObjects& topPseudoObjects = *entry->topPseudoObjects;
            const VLV& allPseudoJets =  topPseudoObjects.valuesSet_ ? *topPseudoObjects.allPseudoJets_ : VLV();
            std::vector<int> allPseudoJetIndices, pseudoJetIndices;
            // Retrieve information related to pseudo-top objects
            // Pseudo jets
            allPseudoJetIndices = common::initialiseIndices(allPseudoJets);
            pseudoJetIndices = allPseudoJetIndices;
            common::selectIndices(pseudoJetIndices, allPseudoJets, common::LVeta, selections.genJetEtaCut_, false);
            common::selectIndices(pseudoJetIndices, allPseudoJets, common::LVeta, -selections.genJetEtaCut_);
            common::selectIndices(pseudoJetIndices, allPseudoJets, common::LVpt, selections.genJetPtCut_);

            if(selections.genDeltaRLeptonJetCut_ > 0.){
                // Vector of pseudoLeptons from which pseudoJets need to be separated in deltaR
                VLV allPseudoLeptons;
                if(topPseudoObjects.valuesSet_){
                    if(topPseudoObjects.PseudoLepton_) allPseudoLeptons.push_back(*topPseudoObjects.PseudoLepton_);
                    if(topPseudoObjects.PseudoAntiLepton_) allPseudoLeptons.push_back(*topPseudoObjects.PseudoAntiLepton_);
                }
                entry->analysisBase->leptonCleanedJetIndices(pseudoJetIndices, allPseudoJets, allPseudoLeptons, selections.genDeltaRLeptonJetCut_);
            }
            common::orderIndices(pseudoJetIndices, allPseudoJets, common::LVpt);

            TLorentzVector lvPseudoTop = common::LVtoTLV(*entry->topPseudoObjects->PseudoTop_);
            TLorentzVector lvPseudoTopBar = common::LVtoTLV(*entry->topPseudoObjects->PseudoAntiTop_);
            TLorentzVector lvPseudoLep = common::LVtoTLV(*entry->topPseudoObjects->PseudoLepton_);
            TLorentzVector lvPseudoLepBar = common::LVtoTLV(*entry->topPseudoObjects->PseudoAntiLepton_);
            TLorentzVector lvPseudoBot = common::LVtoTLV(*entry->topPseudoObjects->PseudoBJet_);
            TLorentzVector lvPseudoBotBar = common::LVtoTLV(*entry->topPseudoObjects->PseudoAntiBJet_);

            for(uint idx=0; idx<pseudoJetIndices.size(); ++idx){
                if(abs(allPseudoJets.at(pseudoJetIndices.at(idx)).pt() - entry->topPseudoObjects->PseudoBJet_->pt())>0.01){
                    if(abs(allPseudoJets.at(pseudoJetIndices.at(idx)).pt() - entry->topPseudoObjects->PseudoAntiBJet_->pt())>0.01){
                        Vars.pseudo_additional_jets_pt[Vars.n_pseudo_additional_jets]=(allPseudoJets.at(pseudoJetIndices.at(idx)).pt());
                        Vars.pseudo_additional_jets_eta[Vars.n_pseudo_additional_jets]=(allPseudoJets.at(pseudoJetIndices.at(idx)).eta());
                        Vars.pseudo_additional_jets_phi[Vars.n_pseudo_additional_jets]=(allPseudoJets.at(pseudoJetIndices.at(idx)).phi());
                        Vars.pseudo_additional_jets_m[Vars.n_pseudo_additional_jets]=(allPseudoJets.at(pseudoJetIndices.at(idx)).M());
                        ++Vars.n_pseudo_additional_jets;
                    }
                }
            }
            if(Vars.n_pseudo_additional_jets>0){
                Vars.pseudo_additional_jet_pt = allPseudoJets.at(pseudoJetIndices.at(0)).pt();
                Vars.pseudo_additional_jet_eta = allPseudoJets.at(pseudoJetIndices.at(0)).eta();
                Vars.pseudo_additional_jet_phi = allPseudoJets.at(pseudoJetIndices.at(0)).phi();
                Vars.pseudo_additional_jet_m = allPseudoJets.at(pseudoJetIndices.at(0)).M();
            }

            Vars.pseudo_top_pt = lvPseudoTop.Pt();
            Vars.pseudo_top_eta = lvPseudoTop.Eta();
            Vars.pseudo_top_phi = lvPseudoTop.Phi();
            Vars.pseudo_top_m = lvPseudoTop.M();
            Vars.pseudo_antitop_pt = lvPseudoTopBar.Pt();
            Vars.pseudo_antitop_eta = lvPseudoTopBar.Eta();
            Vars.pseudo_antitop_phi = lvPseudoTopBar.Phi();
            Vars.pseudo_antitop_m = lvPseudoTopBar.M();
            Vars.pseudo_lepton_pt = lvPseudoLep.Pt();
            Vars.pseudo_lepton_eta = lvPseudoLep.Eta();
            Vars.pseudo_lepton_phi = lvPseudoLep.Phi();
            Vars.pseudo_lepton_m = lvPseudoLep.M();
            Vars.pseudo_antilepton_pt = lvPseudoLepBar.Pt();
            Vars.pseudo_antilepton_eta = lvPseudoLepBar.Eta();
            Vars.pseudo_antilepton_phi = lvPseudoLepBar.Phi();
            Vars.pseudo_antilepton_m = lvPseudoLepBar.M();
            Vars.pseudo_b_pt = lvPseudoBot.Pt();
            Vars.pseudo_b_eta = lvPseudoBot.Eta();
            Vars.pseudo_b_phi = lvPseudoBot.Phi();
            Vars.pseudo_b_m = lvPseudoBot.M();
            Vars.pseudo_antib_pt = lvPseudoBotBar.Pt();
            Vars.pseudo_antib_eta = lvPseudoBotBar.Eta();
            Vars.pseudo_antib_phi = lvPseudoBotBar.Phi();
            Vars.pseudo_antib_m = lvPseudoBotBar.M();

            if(Vars.n_pseudo_additional_jets>0){
                if(Vars.pseudo_additional_jets_pt[0]>0. && lvPseudoTop.Pt()>0. && lvPseudoTopBar.Pt()>0.){
                    TLorentzVector tempRJ(0.,0.,0.,0.);
                    tempRJ.SetPtEtaPhiM(Vars.pseudo_additional_jets_pt[0], Vars.pseudo_additional_jets_eta[0], Vars.pseudo_additional_jets_phi[0], Vars.pseudo_additional_jets_m[0]);
                    Vars.pseudo_rho =  340./(lvPseudoTop+lvPseudoTopBar+tempRJ).M();
                }else{
                    Vars.pseudo_rho =  -999.;
                }
            }else{
                Vars.pseudo_rho =  -999.;
            }

        }
    }
    //=======================PSEUDO-TOP: END==========================//



    std::vector<int> genExtraJetIndices;
    std::vector<int> genExtraJetWithNuIndices;

    if(isTTbarSample_){
        if(entry->topGenObjects->valuesSet_){
            TLorentzVector lvTop    = common::LVtoTLV(*entry->topGenObjects->GenTop_);
            TLorentzVector lvTopBar = common::LVtoTLV(*entry->topGenObjects->GenAntiTop_);
            TLorentzVector lvLep    = common::LVtoTLV(*entry->topGenObjects->GenLepton_);
            TLorentzVector lvLepBar = common::LVtoTLV(*entry->topGenObjects->GenAntiLepton_);
            TLorentzVector lvNeu    = common::LVtoTLV(*entry->topGenObjects->GenNeutrino_);
            TLorentzVector lvNeuBar = common::LVtoTLV(*entry->topGenObjects->GenAntiNeutrino_);
            TLorentzVector lvBot(0.,0.,0.,0.);
            TLorentzVector lvBotBar(0.,0.,0.,0.);
            if(entry->genObjectIndices->genBjetFromTopIndex_>0 && entry->genObjectIndices->genAntiBjetFromTopIndex_>0){
                lvBot    = common::LVtoTLV(entry->topGenObjects->allGenJets_->at(entry->genObjectIndices->genBjetFromTopIndex_));
                lvBotBar = common::LVtoTLV(entry->topGenObjects->allGenJets_->at(entry->genObjectIndices->genAntiBjetFromTopIndex_));
            }

            Vars.gen_top_pt           = lvTop.Pt();
            Vars.gen_top_eta          = lvTop.Eta();
            Vars.gen_top_phi          = lvTop.Phi();
            Vars.gen_top_m            = lvTop.M();
            Vars.gen_antitop_pt       = lvTopBar.Pt();
            Vars.gen_antitop_eta      = lvTopBar.Eta();
            Vars.gen_antitop_phi      = lvTopBar.Phi();
            Vars.gen_antitop_m        = lvTopBar.M();
            Vars.gen_lepton_pt        = lvLep.Pt();
            Vars.gen_lepton_eta       = lvLep.Eta();
            Vars.gen_lepton_phi       = lvLep.Phi();
            Vars.gen_lepton_m         = lvLep.M();
            Vars.gen_antilepton_pt    = lvLepBar.Pt();
            Vars.gen_antilepton_eta   = lvLepBar.Eta();
            Vars.gen_antilepton_phi   = lvLepBar.Phi();
            Vars.gen_antilepton_m     = lvLepBar.M();
            Vars.gen_b_pt             = lvBot.Pt();
            Vars.gen_b_eta            = lvBot.Eta();
            Vars.gen_b_phi            = lvBot.Phi();
            Vars.gen_b_m              = lvBot.M();
            Vars.gen_antib_pt         = lvBotBar.Pt();
            Vars.gen_antib_eta        = lvBotBar.Eta();
            Vars.gen_antib_phi        = lvBotBar.Phi();
            Vars.gen_antib_m          = lvBotBar.M();
            Vars.gen_neutrino_pt      = lvNeu.Pt();
            Vars.gen_neutrino_eta     = lvNeu.Eta();
            Vars.gen_neutrino_phi     = lvNeu.Phi();
            Vars.gen_neutrino_m       = lvNeu.M();
            Vars.gen_antineutrino_pt  = lvNeuBar.Pt();
            Vars.gen_antineutrino_eta = lvNeuBar.Eta();
            Vars.gen_antineutrino_phi = lvNeuBar.Phi();
            Vars.gen_antineutrino_m   = lvNeuBar.M();
            // Vars.gen_nVtx = entry->analysisBase->vertMultiTrue();
            // Vars.gen_nJets = (int)(*entry->topGenObjects->allGenJets_).size(); //correct?
            // Vars.gen_nBJets =  entry->topGenObjects->genBHadJetIndex_->size(); //correct?


            TLorentzVector lvMet = TLorentzVector(0.,0.,0.,0.);
            lvMet.SetPtEtaPhiM(entry->topGenObjects->GenMet_->Pt(),0.,entry->topGenObjects->GenMet_->Phi(),0.);

            Vars.gen_met_pt  = lvMet.Pt();
            Vars.gen_met_eta = 0.;
            Vars.gen_met_phi = lvMet.Phi();
            Vars.gen_met_m   = 0.;

            //do gen_rho calculations

            // Jet matchings for ttbar system
            int genBjetWithNuFromTopIndex(-1);
            int genAntiBjetWithNuFromTopIndex(-1);
            if(entry->topGenObjects->valuesSet_){
                genBjetWithNuFromTopIndex     = entry->analysisBase->genBjetIndex(*entry->topGenObjects, 6, true);
                genAntiBjetWithNuFromTopIndex = entry->analysisBase->genBjetIndex(*entry->topGenObjects, -6, true);
            }

            if (genBjetWithNuFromTopIndex == -2) genBjetWithNuFromTopIndex = -1;// FIXME: generatorTopEvent() and generatorTTbarjetsEvent() can't handle "-2"
            if (genAntiBjetWithNuFromTopIndex == -2) genAntiBjetWithNuFromTopIndex = -1;// -1 - no b-jet found, -2 - two b-jets from one t-quark in generator


            // if(isTTbarSample_){
                // if(entry->topGenObjects->valuesSet_){
                    // Generated jets
                    // std::vector<int> genExtraJetIndices;
                    generatorExtraJets(*entry->topGenObjects,genExtraJetIndices, entry->genObjectIndices->genBjetFromTopIndex_, entry->genObjectIndices->genAntiBjetFromTopIndex_);
                    common::selectIndices(genExtraJetIndices, allGenJets, common::LVeta, selections.genJetEtaCut_, false);
                    common::selectIndices(genExtraJetIndices, allGenJets, common::LVeta, -selections.genJetEtaCut_);
                    common::selectIndices(genExtraJetIndices, allGenJets, common::LVpt, selections.lead2JetPtCut_);

                    if(selections.genDeltaRLeptonJetCut_ > 0.){
                        VLV allGenLeptons;
                        if(entry->topGenObjects->valuesSet_){
                            if(entry->topGenObjects->GenLepton_)     allGenLeptons.push_back(*entry->topGenObjects->GenLepton_);
                            if(entry->topGenObjects->GenAntiLepton_) allGenLeptons.push_back(*entry->topGenObjects->GenAntiLepton_);
                        }
                        entry->analysisBase->leptonCleanedJetIndices(genExtraJetIndices, allGenJets, allGenLeptons, selections.genDeltaRLeptonJetCut_);
                    }
                    common::orderIndices(genExtraJetIndices, allGenJets, common::LVpt);

                    // Generated jets WithNu
                    // std::vector<int> genExtraJetWithNuIndices;
                    generatorExtraJets(*entry->topGenObjects, genExtraJetWithNuIndices, genBjetWithNuFromTopIndex, genAntiBjetWithNuFromTopIndex, true);
                    common::selectIndices(genExtraJetWithNuIndices, allGenJetsWithNu, common::LVeta, selections.genJetEtaCut_, false);
                    common::selectIndices(genExtraJetWithNuIndices, allGenJetsWithNu, common::LVeta, -selections.genJetEtaCut_);
                    common::selectIndices(genExtraJetWithNuIndices, allGenJetsWithNu, common::LVpt,  selections.lead2JetPtCut_);

                    if(selections.genDeltaRLeptonJetCut_ > 0.){
                        VLV allGenLeptons;
                        if(entry->topGenObjects->valuesSet_){
                            if(entry->topGenObjects->GenLepton_)     allGenLeptons.push_back(*entry->topGenObjects->GenLepton_);
                            if(entry->topGenObjects->GenAntiLepton_) allGenLeptons.push_back(*entry->topGenObjects->GenAntiLepton_);
                        }
                        entry->analysisBase->leptonCleanedJetIndices(genExtraJetWithNuIndices, allGenJetsWithNu, allGenLeptons, selections.genDeltaRLeptonJetCut_);
                    }
                    common::orderIndices(genExtraJetWithNuIndices, allGenJetsWithNu, common::LVpt);

                    // Generated extra parton level jets
                    if(entry->topGenObjects->valuesSet_ && entry->topGenObjects->genExtraJetsPartonLevelNoTop_){
                        const VLV& allPartonLevelExtraJets =  entry->topGenObjects->valuesSet_ ? *entry->topGenObjects->genExtraJetsPartonLevelNoTop_ : VLV();
                        std::vector<int> allPartonLevelExtraJetIndices, partonLevelExtraJetIndices;
                        allPartonLevelExtraJetIndices = common::initialiseIndices(allPartonLevelExtraJets);
                        partonLevelExtraJetIndices    = allPartonLevelExtraJetIndices;
                        common::selectIndices(partonLevelExtraJetIndices, allPartonLevelExtraJets, common::LVeta, selections.genJetEtaCut_, false);
                        common::selectIndices(partonLevelExtraJetIndices, allPartonLevelExtraJets, common::LVeta, -selections.genJetEtaCut_);
                        common::selectIndices(partonLevelExtraJetIndices, allPartonLevelExtraJets, common::LVpt, selections.genJetPtCut_);


                        // if(selections.genDeltaRLeptonJetCut_ > 0.){
                        //     // Vector of pseudoLeptons from which pseudoJets need to be separated in deltaR
                        //     VLV allPseudoLeptons;
                        //     if(topPseudoObjects.valuesSet_){
                        //         if(topGenObjects.PseudoLepton_) allPseudoLeptons.push_back(*topGenObjects.PseudoLepton_);
                        //         if(topGenObjects.PseudoAntiLepton_) allPseudoLeptons.push_back(*topGenObjects.PseudoAntiLepton_);
                        //     }
                        //     entry->analysisBase->leptonCleanedJetIndices(partonLevelExtraJetIndices, allPartonLevelExtraJets, allPseudoLeptons, selections.genDeltaRLeptonJetCut_);
                        // }
                        common::orderIndices(partonLevelExtraJetIndices, allPartonLevelExtraJets, common::LVpt);
                        // Generated parton level extra jets
                        if(entry->topGenObjects->genExtraJetsPartonLevelNoTop_){
                            const auto jets_p4 = *(entry->topGenObjects->genExtraJetsPartonLevelNoTop_);
                            for(uint idx=0; idx<partonLevelExtraJetIndices.size(); ++idx){
                                  Vars.gen_partonLevel_additional_jets_pt[Vars.n_gen_partonLevel_additional_jets]=jets_p4.at(partonLevelExtraJetIndices.at(idx)).pt();
                                  Vars.gen_partonLevel_additional_jets_eta[Vars.n_gen_partonLevel_additional_jets]=jets_p4.at(partonLevelExtraJetIndices.at(idx)).eta();
                                  Vars.gen_partonLevel_additional_jets_phi[Vars.n_gen_partonLevel_additional_jets]=jets_p4.at(partonLevelExtraJetIndices.at(idx)).phi();
                                  Vars.gen_partonLevel_additional_jets_m[Vars.n_gen_partonLevel_additional_jets]=jets_p4.at(partonLevelExtraJetIndices.at(idx)).M();
                                  ++Vars.n_gen_partonLevel_additional_jets;
                            }
                            if(Vars.n_gen_partonLevel_additional_jets>0){
                                Vars.gen_partonLevel_additional_jet_pt = jets_p4.at(partonLevelExtraJetIndices.at(0)).pt();
                                Vars.gen_partonLevel_additional_jet_eta = jets_p4.at(partonLevelExtraJetIndices.at(0)).eta();
                                Vars.gen_partonLevel_additional_jet_phi = jets_p4.at(partonLevelExtraJetIndices.at(0)).phi();
                                Vars.gen_partonLevel_additional_jet_m = jets_p4.at(partonLevelExtraJetIndices.at(0)).M();
                            }
                        }
                    }

                    // Generated jets
                    if(entry->topGenObjects->allGenJets_){
                        const auto jets_p4 = *(entry->topGenObjects->allGenJets_);
                        for(uint idx=0; idx<genExtraJetIndices.size(); ++idx){
                              Vars.gen_additional_jets_pt[Vars.n_gen_additional_jets]=(jets_p4.at(genExtraJetIndices.at(idx)).pt());
                              Vars.gen_additional_jets_eta[Vars.n_gen_additional_jets]=(jets_p4.at(genExtraJetIndices.at(idx)).eta());
                              Vars.gen_additional_jets_phi[Vars.n_gen_additional_jets]=(jets_p4.at(genExtraJetIndices.at(idx)).phi());
                              Vars.gen_additional_jets_m[Vars.n_gen_additional_jets]=(jets_p4.at(genExtraJetIndices.at(idx)).M());
                              ++Vars.n_gen_additional_jets;
                        }
                        if(Vars.n_gen_additional_jets>0){
                            Vars.gen_additional_jet_pt = jets_p4.at(genExtraJetIndices.at(0)).pt();
                            Vars.gen_additional_jet_eta = jets_p4.at(genExtraJetIndices.at(0)).eta();
                            Vars.gen_additional_jet_phi = jets_p4.at(genExtraJetIndices.at(0)).phi();
                            Vars.gen_additional_jet_m = jets_p4.at(genExtraJetIndices.at(0)).M();
                        }
                    }
                    // Generated jets with nu
                    if(entry->topGenObjects->allGenJetsWithNu_){
                        const auto jets_p4 = *(entry->topGenObjects->allGenJetsWithNu_);
                        for(uint idx=0; idx<genExtraJetWithNuIndices.size(); ++idx){
                              Vars.gen_additional_jets_withNu_pt[Vars.n_gen_additional_jets_withNu]=(jets_p4.at(genExtraJetWithNuIndices.at(idx)).pt());
                              Vars.gen_additional_jets_withNu_eta[Vars.n_gen_additional_jets_withNu]=(jets_p4.at(genExtraJetWithNuIndices.at(idx)).eta());
                              Vars.gen_additional_jets_withNu_phi[Vars.n_gen_additional_jets_withNu]=(jets_p4.at(genExtraJetWithNuIndices.at(idx)).phi());
                              Vars.gen_additional_jets_withNu_m[Vars.n_gen_additional_jets_withNu]=(jets_p4.at(genExtraJetWithNuIndices.at(idx)).M());
                              ++Vars.n_gen_additional_jets_withNu;
                        }
                        if(Vars.n_gen_additional_jets_withNu>0){
                            Vars.gen_additional_jet_withNu_pt = jets_p4.at(genExtraJetWithNuIndices.at(0)).pt();
                            Vars.gen_additional_jet_withNu_eta = jets_p4.at(genExtraJetWithNuIndices.at(0)).eta();
                            Vars.gen_additional_jet_withNu_phi = jets_p4.at(genExtraJetWithNuIndices.at(0)).phi();
                            Vars.gen_additional_jet_withNu_m = jets_p4.at(genExtraJetWithNuIndices.at(0)).M();
                        }
                    }



                // }
            // }
            if(Vars.n_gen_additional_jets>0){
                      TLorentzVector tempRJ(0.,0.,0.,0.);
                      tempRJ.SetPtEtaPhiM(Vars.gen_additional_jets_pt[0],Vars.gen_additional_jets_eta[0],Vars.gen_additional_jets_phi[0],Vars.gen_additional_jets_m[0]);
                      Vars.gen_rho =  340./(common::LVtoTLV(*entry->topGenObjects->GenTop_)+common::LVtoTLV(*entry->topGenObjects->GenAntiTop_)+tempRJ).M();
            }else{
                Vars.gen_rho =  -999.;
            }
            if(Vars.n_gen_additional_jets_withNu>0){
                      TLorentzVector tempRJ(0.,0.,0.,0.);
                      tempRJ.SetPtEtaPhiM(Vars.gen_additional_jets_withNu_pt[0],Vars.gen_additional_jets_withNu_eta[0],Vars.gen_additional_jets_withNu_phi[0],Vars.gen_additional_jets_withNu_m[0]);
                      Vars.gen_rhoWithNu =  340./(common::LVtoTLV(*entry->topGenObjects->GenTop_)+common::LVtoTLV(*entry->topGenObjects->GenAntiTop_)+tempRJ).M();
            }else{
                Vars.gen_rhoWithNu =  -999.;
            }
            if(Vars.n_gen_partonLevel_additional_jets>0){
                      TLorentzVector tempRJ(0.,0.,0.,0.);
                      tempRJ.SetPtEtaPhiM(Vars.gen_partonLevel_additional_jets_pt[0],Vars.gen_partonLevel_additional_jets_eta[0],Vars.gen_partonLevel_additional_jets_phi[0],Vars.gen_partonLevel_additional_jets_m[0]);
                      Vars.gen_partonLevel_rho =  340./(common::LVtoTLV(*entry->topGenObjects->GenTop_)+common::LVtoTLV(*entry->topGenObjects->GenAntiTop_)+tempRJ).M();
            }else{
                Vars.gen_partonLevel_rho =  -999.;
            }
        }



        // }
    }

    //=======================================================
    //============= do Reco selection here again=============
    std::vector<std::vector<int> > allGenJetBhadronIndices;
    std::vector<std::vector<int> > genJetBhadronIndices;
    std::vector<int> allGenBjetIndices;
    std::vector<int> genBjetIndices;
    std::vector<std::vector<int> > allGenJetWithNuBhadronIndices;
    std::vector<std::vector<int> > genJetWithNuBhadronIndices;
    std::vector<int> allGenBjetWithNuIndices;
    std::vector<int> genBjetWithNuIndices;
    if(entry->topGenObjects->valuesSet_){
        allGenJetBhadronIndices = entry->analysisBase->matchHadronsToGenJets(allGenJetIndices, allGenJets, *entry->topGenObjects->genBHadJetIndex_);
        genJetBhadronIndices = entry->analysisBase->matchHadronsToGenJets(genJetIndices, allGenJets, *entry->topGenObjects->genBHadJetIndex_);
        allGenBjetIndices = entry->analysisBase->genBjetIndices(allGenJetBhadronIndices);
        genBjetIndices = entry->analysisBase->genBjetIndices(genJetBhadronIndices);

        allGenJetWithNuBhadronIndices = entry->analysisBase->matchHadronsToGenJets(allGenJetWithNuIndices, allGenJetsWithNu, *entry->topGenObjects->genBHadJetIndex_);
        genJetWithNuBhadronIndices = entry->analysisBase->matchHadronsToGenJets(genJetWithNuIndices, allGenJetsWithNu, *entry->topGenObjects->genBHadJetIndex_);
        allGenBjetWithNuIndices = entry->analysisBase->genBjetIndices(allGenJetWithNuBhadronIndices);
        genBjetWithNuIndices = entry->analysisBase->genBjetIndices(genJetWithNuBhadronIndices);
    }




    bool failReco=false;
    const RecoObjects& recoObjects = *entry->recoObjects;
    //===CUT===
    //selectionStep = "1a";
    // Check if event should be excluded due to MET filters
    if(entry->analysisBase->failsMetFilters(entry->entryNumber)){
        failReco=true;
    }

    //===CUT===
    //selectionStep = "1b";
    // Check if event was triggered
    if(entry->analysisBase->failsDileptonTrigger(entry->entryNumber)){
        failReco=true;

    }

    //===CUT===
    // selectionStep = "1";
    // Check if the first primary vertex is of good quality
    if(!entry->analysisBase->firstVertexIsGood(entry->entryNumber)){
        failReco=true;
    }



    // === FULL RECO OBJECT SELECTION === (can thus be used at each selection step)

    // Access objects info
    // Get allLepton indices, apply selection cuts and order them by pt (beginning with the highest value)
    const VLV& allLeptons = *recoObjects.allLeptons_;
    const std::vector<int>& lepPdgId = *recoObjects.lepPdgId_;
    const std::vector<float>& lepSCEta = *recoObjects.lepSCEta_; // FIXME: works since v037
    std::vector<int> allLeptonIndices = common::initialiseIndices(allLeptons);
    if(selections.muonIsoCut_ > -999.){
        std::vector<int> allElectronIndices = allLeptonIndices;
        std::vector<int> onlyMuonIndices = allLeptonIndices;
        std::vector<int> onlyAntiMuonIndices = allLeptonIndices;
        common::selectIndices(allElectronIndices, lepPdgId, 11, false);
        common::selectIndices(allElectronIndices, lepPdgId, -11);
        common::selectIndices(onlyMuonIndices, lepPdgId, 13);
        common::selectIndices(onlyAntiMuonIndices, lepPdgId, -13, false);
        std::vector<int> allMuonIndices = common::mergeIndices(onlyMuonIndices, onlyAntiMuonIndices);
        common::selectIndices(allMuonIndices, *recoObjects.lepIso_, selections.muonIsoCut_, false);
        allLeptonIndices.clear();
        allLeptonIndices = common::mergeIndices(allElectronIndices, allMuonIndices);
    }
    common::selectMuonIndicesID    (allLeptonIndices, recoObjects.lepID_MuonTight_   , 1, &lepPdgId);
    common::selectElectronIndicesID(allLeptonIndices, recoObjects.lepID_ElecCutBased_, 4, &lepPdgId);
    common::selectIndices(allLeptonIndices, allLeptons, common::LVeta, selections.leptonEtaCut_, false);
    common::selectIndices(allLeptonIndices, allLeptons, common::LVeta, -selections.leptonEtaCut_);
    common::selectIndices(allLeptonIndices, allLeptons, common::LVpt, selections.leptonPtCut_);
    common::orderIndices(allLeptonIndices, allLeptons, common::LVpt);

    // Get indices of leptons and antiLeptons separated by charge, and get the leading ones if they exist
    std::vector<int> leptonIndices = allLeptonIndices;
    std::vector<int> antiLeptonIndices = allLeptonIndices;
    common::selectIndices(leptonIndices, lepPdgId, 0);
    common::selectIndices(antiLeptonIndices, lepPdgId, 0, false);

    // For 8 TeV lepton pair definition
    const int numberOfLeptons = leptonIndices.size();
    const int numberOfAntiLeptons = antiLeptonIndices.size();
    const int leptonIndex = numberOfLeptons>0 ? leptonIndices.at(0) : -1;
    const int antiLeptonIndex = numberOfAntiLeptons>0 ? antiLeptonIndices.at(0) : -1;

    // For 13 TeV lepton pair definition
    const int LeptonsNumber = allLeptonIndices.size();
    const int firstLeptonIndex = LeptonsNumber>0 ? allLeptonIndices.at(0) : -1;
    const int secondLeptonIndex = LeptonsNumber>1 ? allLeptonIndices.at(1) : -1;

    // In case of an existing opposite-charge dilepton system, get indices for leading and next-to-leading lepton
    int leadingLeptonIndex(-1);
    int nLeadingLeptonIndex(-1);
    if(entry->analysisConfig->general().era_ == Era::run1_8tev && (numberOfLeptons>0 && numberOfAntiLeptons>0)){
        leadingLeptonIndex = leptonIndex;
        nLeadingLeptonIndex = antiLeptonIndex;
        common::orderIndices(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, common::LVpt);
    } else if((entry->analysisConfig->general().era_ == Era::run2_13tev_50ns || entry->analysisConfig->general().era_ == Era::run2_13tev_25ns
        || entry->analysisConfig->general().era_ == Era::run2_13tev_2015_25ns || entry->analysisConfig->general().era_ == Era::run2_13tev_2016_25ns
        || entry->analysisConfig->general().era_ == Era::run2_13tev_2016 || entry->analysisConfig->general().era_ == Era::run2_13tev_2017 || entry->analysisConfig->general().era_ == Era::run2_13tev_2018)
        && LeptonsNumber>1){
        leadingLeptonIndex = firstLeptonIndex;
        nLeadingLeptonIndex = secondLeptonIndex;
        common::orderIndices(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, common::LVpt);
    }

    const bool hasLeptonPair = entry->analysisBase->hasLeptonPair(leadingLeptonIndex, nLeadingLeptonIndex, lepPdgId);
    bool passLeadingLeptonPtCut = false;
    if(leadingLeptonIndex > -1) passLeadingLeptonPtCut = allLeptons.at(leadingLeptonIndex).pt() > entry->analysisConfig->selections().leadingLeptonPtCut_;
    const bool hasMoreThan2Leptons = LeptonsNumber > 2; // Veto based on this quantity needed in order to match fiducial phase space definition for differential xsec

    // Get two indices of the two leptons in the right order for trigger scale factor, if existing
    int leptonXIndex(leadingLeptonIndex);
    int leptonYIndex(nLeadingLeptonIndex);
    if(hasLeptonPair){
        //in ee and mumu channel leptonX must be the highest pt lepton, i.e. this is already correct
        // in emu channel leptonX must be electron
        if(std::abs(lepPdgId.at(leptonXIndex)) != std::abs(lepPdgId.at(leptonYIndex))){
            common::orderIndices(leptonYIndex, leptonXIndex, lepPdgId, true);
        }
    }

    // Get dilepton system, if existing
    const LV dummyLV(0.,0.,0.,0.);
    const LV dilepton(hasLeptonPair ? allLeptons.at(leadingLeptonIndex)+allLeptons.at(nLeadingLeptonIndex) : dummyLV);

    // Get jet indices, apply selection cuts and order them by pt (beginning with the highest value)
    const VLV& jets = *recoObjects.jets_;
    std::vector<int> jetIndices = common::initialiseIndices(jets);
    if((entry->analysisConfig->general().era_ == Era::run2_13tev_50ns || entry->analysisConfig->general().era_ == Era::run2_13tev_25ns
      || entry->analysisConfig->general().era_ == Era::run2_13tev_2015_25ns || entry->analysisConfig->general().era_ == Era::run2_13tev_2017
      || entry->analysisConfig->general().era_ == Era::run2_13tev_2016_25ns || entry->analysisConfig->general().era_ == Era::run2_13tev_2016
      || entry->analysisConfig->general().era_ == Era::run2_13tev_2018) && selections.deltaRLeptonJetCut_ > 0.){
        // Vector of leptons from which jets need to be separated in deltaR
        VLV leptonsForJetCleaning;
        for(const int index : allLeptonIndices) leptonsForJetCleaning.push_back(allLeptons.at(index));
        entry->analysisBase->leptonCleanedJetIndices(jetIndices, jets, leptonsForJetCleaning, selections.deltaRLeptonJetCut_);
    }
    std::vector<int> jetIndicesForExtraJetStudies = jetIndices; // remains the same as initial collection at 8 TeV, but includes DeltaR(l,j) cleaning at 13 TeV: before application of other cuts
    common::selectIndices(jetIndices, jets, common::LVeta, selections.jetEtaCut_, false);
    common::selectIndices(jetIndices, jets, common::LVeta, -selections.jetEtaCut_);
    common::selectIndices(jetIndices, jets, common::LVpt, selections.jetPtCut_);
    if(recoObjects.jetPFID_){

      common::selectIndices(jetIndices, *recoObjects.jetPFID_, selections.jetPFIDWorkingPoint_);
    }
    else {
      throw std::runtime_error("MiniTree::Fill -- null pointer in recoObjects.jetPFID_");
    }
    const JetPileupId::WorkingPoint jetPileupIdWorkingPoint(selections.jetPileupIdWorkingPoint_);
    if(jetPileupIdWorkingPoint != JetPileupId::none){
        const std::vector<int>& jetPileupIdFlags = *recoObjects.jetPileupIdFlag_;
        common::selectIndices(jetIndices, jetPileupIdFlags, JetPileupId::cutValue(entry->analysisConfig->selections().jetPileupIdWorkingPoint_));
    }


    common::orderIndices(jetIndices, jets, common::LVpt);
    const int numberOfJets = jetIndices.size();
    const bool has2Jets = numberOfJets > 1;

    //Anya: 2D iso
    // if(btopTag) entry->analysisBase->lepton2dIsolationIndices(allLeptonIndices, allLeptons, jets);

    // Get b-jet indices, apply selection cuts
    // and apply b-tag efficiency MC correction using random number based tag flipping (if requested correction mode is applied)
    // and order b-jets by btag discriminator (beginning with the highest value)
    const std::vector<float>& jetBtags = *recoObjects.jetBtags_;
    const std::vector<int>& jetFlavour = *entry->commonGenObjects->jetFlavour_;
    std::vector<int> bjetIndices = jetIndices;
    common::selectIndices(bjetIndices, jetBtags, entry->analysisBase->btagCutValue());
    entry->analysisBase->retagJets(bjetIndices, jetIndices, jets, jetFlavour, jetBtags);
    common::orderIndices(bjetIndices, jetBtags);
    const int numberOfBjets = bjetIndices.size();
    const bool hasBtag = numberOfBjets > 0;

    // Get MET
    // const LV& met = *recoObjects.met_;
    const LV& uncormet = *recoObjects.met_;
    LV met = uncormet;
    // if(entry->analysisConfig->selections().metAlgorithm_==Met::pf){
    //     met=common::metWithXYCorrection(met, entry->eventMetadata->runNumber_, entry->analysisBase->era_, isGen_, recoObjects.vertMulti_);
    // }
    const bool hasMetOrEmu = entry->analysisBase->channel()==Channel::emu || met.Pt() > selections.metCut_;

    const bool isZregion = dilepton.M() > 76. && dilepton.M() < 106.;

    const ttbar::RecoObjectIndices recoObjectIndices(allLeptonIndices,
                                                     leptonIndices, antiLeptonIndices,
                                                     leptonIndex, antiLeptonIndex,
                                                     leadingLeptonIndex, nLeadingLeptonIndex,
                                                     leptonXIndex, leptonYIndex,
                                                     {-99},
                                                     jetIndices, bjetIndices);

    const std::vector<int>& genJetMatchedRecoBjetIndices = entry->analysisBase->matchRecoToGenJets(jetIndices, jets, allGenBjetIndices, allGenJets);
    const std::vector<int>& genJetWithNuMatchedRecoBjetIndices = entry->analysisBase->matchRecoToGenJets(jetIndices, jets, allGenBjetWithNuIndices, allGenJetsWithNu);

    int matchedBjetFromTopIndex = entry->genObjectIndices->genBjetFromTopIndex_>=0 ? genJetMatchedRecoBjetIndices.at(entry->genObjectIndices->genBjetFromTopIndex_) : -1;
    int matchedAntiBjetFromTopIndex = entry->genObjectIndices->genAntiBjetFromTopIndex_>=0 ? genJetMatchedRecoBjetIndices.at(entry->genObjectIndices->genAntiBjetFromTopIndex_) : -1;

    // // Jet matchings for ttbar system
    // int genBjetWithNuFromTopIndex(-1);
    // int genAntiBjetWithNuFromTopIndex(-1);
    // if(entry->topGenObjects->valuesSet_){
    //     genBjetWithNuFromTopIndex = entry->analysisBase->genBjetIndex(*entry->topGenObjects, 6, true);
    //     genAntiBjetWithNuFromTopIndex = entry->analysisBase->genBjetIndex(*entry->topGenObjects, -6, true);
    // }
    //
    // if (genBjetWithNuFromTopIndex == -2) genBjetWithNuFromTopIndex = -1;// FIXME: generatorTopEvent() and generatorTTbarjetsEvent() can't handle "-2"
    // if (genAntiBjetWithNuFromTopIndex == -2) genAntiBjetWithNuFromTopIndex = -1;// -1 - no b-jet found, -2 - two b-jets from one t-quark in generator


    // uint matchedBjetFromTopIndexWithNu = genBjetWithNuFromTopIndex>=0 ? genJetWithNuMatchedRecoBjetIndices.at(genBjetWithNuFromTopIndex) : -1;
    // uint matchedAntiBjetFromTopIndexWithNu = genAntiBjetWithNuFromTopIndex>=0 ? genJetWithNuMatchedRecoBjetIndices.at(genAntiBjetWithNuFromTopIndex) : -1;

    // const ttbar::GenObjectIndices genObjectIndicesWithMatching(entry->genObjectIndices->genBjetFromTopIndex_, entry->genObjectIndices->genAntiBjetFromTopIndex_, matchedBjetFromTopIndex, matchedAntiBjetFromTopIndex, -1, -1, -1, -1,entry->genObjectIndices->genVisJetIndices_);
    // const ttbar::GenObjectIndices genObjectIndicesWithMatchingWithNu(genBjetWithNuFromTopIndex, matchedAntiBjetFromTopIndexWithNu, matchedBjetFromTopIndexWithNu, matchedAntiBjetFromTopIndexWithNu, -1, -1, -1, -1,entry->genObjectIndices->genVisJetIndices_);


    // Determine all reco level weights
    const double weightLeptonSF = entry->analysisBase->weightLeptonSF(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId, lepSCEta); // FIXME: works since v037
    const double weightLeptonSF_Ele_ID_Up = entry->analysisBase->weightLeptonSF(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId, lepSCEta, "ELE_ID_UP"); // FIXME: works since v037
    const double weightLeptonSF_Ele_ID_Dn = entry->analysisBase->weightLeptonSF(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId, lepSCEta, "ELE_ID_DOWN"); // FIXME: works since v037
    const double weightLeptonSF_Muon_ID_Up = entry->analysisBase->weightLeptonSF(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId, lepSCEta, "MUON_ID_UP"); // FIXME: works since v037
    const double weightLeptonSF_Muon_ID_Dn = entry->analysisBase->weightLeptonSF(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId, lepSCEta, "MUON_ID_DOWN"); // FIXME: works since v037
    const double weightLeptonSF_Ele_Reco_Up = entry->analysisBase->weightLeptonSF(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId, lepSCEta, "ELE_RECO_UP"); // FIXME: works since v037
    const double weightLeptonSF_Ele_Reco_Dn = entry->analysisBase->weightLeptonSF(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId, lepSCEta, "ELE_RECO_DOWN"); // FIXME: works since v037
    const double weightLeptonSF_Muon_Iso_Up = entry->analysisBase->weightLeptonSF(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId, lepSCEta, "MUON_ISO_UP"); // FIXME: works since v037
    const double weightLeptonSF_Muon_Iso_Dn = entry->analysisBase->weightLeptonSF(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId, lepSCEta, "MUON_ISO_DOWN"); // FIXME: works since v037
    double weightTriggerSF;
    double weightTriggerSF_Up;
    double weightTriggerSF_Dn;
    if(entry->analysisConfig->corrections().triggerSFReadoutMode_== 2){
      weightTriggerSF = entry->analysisBase->weightTriggerSF(leptonXIndex, leptonYIndex, allLeptons, lepPdgId);
      weightTriggerSF_Up = this->weightTriggerSF_Up(leptonXIndex, leptonYIndex, allLeptons, lepPdgId);
      weightTriggerSF_Dn = this->weightTriggerSF_Dn(leptonXIndex, leptonYIndex, allLeptons, lepPdgId);
    }else{
        weightTriggerSF = 1.0;
        weightTriggerSF_Up = 1.0;
        weightTriggerSF_Dn = 1.0;
    }
    double l1PrefiringWeight = entry->analysisBase->weightL1Prefiring(recoObjects);
    double l1PrefiringWeight_Up = entry->analysisBase->weightL1Prefiring_Up(recoObjects);
    double l1PrefiringWeight_Dn = entry->analysisBase->weightL1Prefiring_Dn(recoObjects);

    const double weightBtagSF        = entry->analysisBase->weightBtagSF(jetIndices, jets, jetFlavour, jetBtags);
    const double weightBtagSF_Up     = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "BTAG_UP");
    const double weightBtagSF_Dn     = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "BTAG_DOWN");
    const double weightBtagSF_LjetUp = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "BTAG_LJET_UP");
    const double weightBtagSF_LjetDn = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "BTAG_LJET_DOWN");

    double weightBtagSF_MERENSCALE_Up(-999.);
    double weightBtagSF_MERENSCALE_Down(-999.);
    double weightBtagSF_MEFACSCALE_Up(-999.);
    double weightBtagSF_MEFACSCALE_Down(-999.);
    double weightBtagSF_MESCALE_Up(-999.);
    double weightBtagSF_MESCALE_Down(-999.);
    double weightBtagSF_BFRAG_Up(-999.);
    double weightBtagSF_BFRAG_Down(-999.);
    double weightBtagSF_BFRAG_Peterson(-999.);
    double weightBtagSF_BSEMILEP_Up(-999.);
    double weightBtagSF_BSEMILEP_Down(-999.);
    double weightBtagSF_PDF_ALPHAS_Up(-999.);
    double weightBtagSF_PDF_ALPHAS_Down(-999.);
    double weightBtagSF_PSSCALE_WEIGHT_4_Up(-999.);
    double weightBtagSF_PSSCALE_WEIGHT_4_Down(-999.);
    double weightBtagSF_PSSCALE_WEIGHT_5_Up(-999.);
    double weightBtagSF_PSSCALE_WEIGHT_5_Down(-999.);
    double weightBtagSF_PU_Up(-999.);
    double weightBtagSF_PU_Down(-999.);
    double weightBtagSF_TRIG_Up(-999.);
    double weightBtagSF_TRIG_Down(-999.);
    double weightBtagSF_ELE_ID_Up(-999.);
    double weightBtagSF_ELE_ID_Down(-999.);
    double weightBtagSF_ELE_RECO_Up(-999.);
    double weightBtagSF_ELE_RECO_Down(-999.);
    double weightBtagSF_MUON_ID_Up(-999.);
    double weightBtagSF_MUON_ID_Down(-999.);
    double weightBtagSF_MUON_ISO_Up(-999.);
    double weightBtagSF_MUON_ISO_Down(-999.);
    double weightBtagSF_L1PREFIRING_Up(-999.);
    double weightBtagSF_L1PREFIRING_Down(-999.);

    if(!isSystSample_ && !isSystVar_){
        if(entry->analysisBase->readInMEWeights_){
            weightBtagSF_MERENSCALE_Up         = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "MERENSCALE_UP");
            weightBtagSF_MERENSCALE_Down       = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "MERENSCALE_DOWN");
            weightBtagSF_MEFACSCALE_Up         = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "MEFACSCALE_UP");
            weightBtagSF_MEFACSCALE_Down       = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "MEFACSCALE_DOWN");
            weightBtagSF_MESCALE_Up            = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "MESCALE_UP");
            weightBtagSF_MESCALE_Down          = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "MESCALE_DOWN");
        }
        if(entry->analysisBase->readInBFragWeights_){
            weightBtagSF_BFRAG_Up              = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "BFRAG_UP");
            weightBtagSF_BFRAG_Down            = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "BFRAG_DOWN");
            weightBtagSF_BFRAG_Peterson        = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "BFRAG_PETERSON");
            weightBtagSF_BSEMILEP_Up           = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "BSEMILEP_UP");
            weightBtagSF_BSEMILEP_Down         = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "BSEMILEP_DOWN");
        }
        if(entry->analysisBase->readInPDFWeights_ && entry->analysisBase->nPDFWeights_> 0){
            weightBtagSF_PDF_ALPHAS_Up         = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "PDF_ALPHAS_UP");
            weightBtagSF_PDF_ALPHAS_Down       = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "PDF_ALPHAS_DOWN");
        }
        if(entry->analysisBase->readInPSWeights_ && entry->analysisBase->nPSWeights_> 0){
            weightBtagSF_PSSCALE_WEIGHT_4_Up   = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "PSSCALE_WEIGHT_4_UP");
            weightBtagSF_PSSCALE_WEIGHT_4_Down = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "PSSCALE_WEIGHT_4_DOWN");
            weightBtagSF_PSSCALE_WEIGHT_5_Up   = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "PSSCALE_WEIGHT_5_UP");
            weightBtagSF_PSSCALE_WEIGHT_5_Down = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "PSSCALE_WEIGHT_5_DOWN");
        }
        weightBtagSF_PU_Up                 = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "PU_UP");
        weightBtagSF_PU_Down               = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "PU_DOWN");
        weightBtagSF_TRIG_Up               = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "TRIG_UP");
        weightBtagSF_TRIG_Down             = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "TRIG_DOWN");
        weightBtagSF_ELE_ID_Up             = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "ELE_ID_UP");
        weightBtagSF_ELE_ID_Down           = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "ELE_ID_DOWN");
        weightBtagSF_ELE_RECO_Up           = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "ELE_RECO_UP");
        weightBtagSF_ELE_RECO_Down         = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "ELE_RECO_DOWN");
        weightBtagSF_MUON_ID_Up            = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "MUON_ID_UP");
        weightBtagSF_MUON_ID_Down          = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "MUON_ID_DOWN");
        weightBtagSF_MUON_ISO_Up           = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "MUON_ISO_UP");
        weightBtagSF_MUON_ISO_Down         = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "MUON_ISO_DOWN");
        weightBtagSF_L1PREFIRING_Up        = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "L1PREFIRING_UP");
        weightBtagSF_L1PREFIRING_Down      = this->weightBtagSF_Var(jetIndices, jets, jetFlavour, jetBtags, "L1PREFIRING_DOWN");
    }

    const double weightKinReco = entry->analysisBase->weightKinReco();
    const double weightKinReco_Up = this->weightKinReco_Up();
    const double weightKinReco_Dn = this->weightKinReco_Dn();
    const double weightLooseKinReco = entry->analysisBase->weightLooseKinReco();


    //===CUT===
    // selectionStep = "2";
    // we need an OS lepton pair
    if (!hasLeptonPair || !passLeadingLeptonPtCut){
        failReco=true;
    }
    if (selections.vetoEventsWithMoreThan2Leptons_ == 1  && hasMoreThan2Leptons){
        failReco=true;
    }

    //===CUT===
    // selectionStep = "3";
    // with at least 20 GeV invariant mass
    if (dilepton.M() < 20.){
        failReco=true;
    }


    //=====================================
    //========= Now fill everything========

    // Vars.channel = entry->analysisBase->channelPdgIdProduct();

    Vars.passStep3=!failReco;
    // Vars.passStep4=Vars.passStep3 && ((Vars.channel== -143) || (!isZregion));
    Vars.passStep4=Vars.passStep3 && ((entry->analysisBase->channelPdgIdProduct()== -143) || (!isZregion));
    Vars.passStep5=Vars.passStep3 && Vars.passStep4 && has2Jets;
    Vars.passStep6=Vars.passStep3 && Vars.passStep4 && Vars.passStep5 && hasMetOrEmu;
    Vars.passStep7=Vars.passStep3 && Vars.passStep4 && Vars.passStep5 && Vars.passStep6 && hasBtag;

    // Vars.eventCounter = *entry->eventCounter;
    Vars.eventNumber = entry->eventMetadata->eventNumber_;
    Vars.lumiNumber = entry->eventMetadata->lumiBlock_;
    Vars.runNumber = entry->eventMetadata->runNumber_;
    Vars.weight = isGen_? entry->genLevelWeights->trueLevelWeightNoPileup_ : entry->weight;
    Vars.leptonSF = weightLeptonSF;
    Vars.l1PrefiringWeight = l1PrefiringWeight;
    Vars.triggerSF = weightTriggerSF;
    Vars.btagSF = weightBtagSF;
    Vars.pileupSF = entry->pileupW;


    if(!isSystSample_ && !isSystVar_){
        // if(isTTbarSample_){
        if(entry->analysisBase->readInPDFWeights_ && entry->analysisBase->nPDFWeights_> 0){
            if(entry->analysisBase->b_pdfWeights->GetEntryNumber() != *entry->eventCounter){
              entry->analysisBase->b_pdfWeights->GetEntry(*entry->eventCounter);
            }
            const std::vector<float>& pdfWeightsFloat = *(entry->analysisBase->pdfWeights_);
            for(unsigned int i = 0; i < pdfWeightsFloat.size(); i++){
              Vars.var_PDF[Vars.n_PDFweights]=(pdfWeightsFloat[i]);
              ++Vars.n_PDFweights;
            }
        }
        // if(isTTbarSample_){
        if(entry->analysisBase->readInPDFWeights_){
            Vars.var_PDF_ALPHAS_UP        = entry->analysisBase->weightAlphasPdf_Up(entry->entryNumber);
            Vars.var_PDF_ALPHAS_DOWN      = entry->analysisBase->weightAlphasPdf_Dn(entry->entryNumber);
            Vars.var_BTAG_PDF_ALPHAS_UP   = weightBtagSF_PDF_ALPHAS_Up;
            Vars.var_BTAG_PDF_ALPHAS_DOWN = weightBtagSF_PDF_ALPHAS_Down;
        }
        if(entry->analysisBase->readInMEWeights_){
            Vars.var_MERENSCALE_UP        = entry->analysisBase->weightMeRenScale_UP(entry->entryNumber);
            Vars.var_MERENSCALE_DOWN      = entry->analysisBase->weightMeRenScale_DN(entry->entryNumber);
            Vars.var_MEFACSCALE_UP        = entry->analysisBase->weightMeFacScale_UP(entry->entryNumber);
            Vars.var_MEFACSCALE_DOWN      = entry->analysisBase->weightMeFacScale_DN(entry->entryNumber);
            Vars.var_MESCALE_UP           = entry->analysisBase->weightMeScale_UP(entry->entryNumber);
            Vars.var_MESCALE_DOWN         = entry->analysisBase->weightMeScale_DN(entry->entryNumber);
            Vars.var_BTAG_MERENSCALE_UP   = weightBtagSF_MERENSCALE_Up;
            Vars.var_BTAG_MERENSCALE_DOWN = weightBtagSF_MERENSCALE_Down;
            Vars.var_BTAG_MEFACSCALE_UP   = weightBtagSF_MEFACSCALE_Up;
            Vars.var_BTAG_MEFACSCALE_DOWN = weightBtagSF_MEFACSCALE_Down;
            Vars.var_BTAG_MESCALE_UP      = weightBtagSF_MESCALE_Up;
            Vars.var_BTAG_MESCALE_DOWN    = weightBtagSF_MESCALE_Down;
        }
        if(entry->analysisBase->readInBFragWeights_){
            Vars.var_BFRAG_UP            = entry->analysisBase->weightBFrag_Up(entry->entryNumber);
            Vars.var_BFRAG_DOWN          = entry->analysisBase->weightBFrag_Dn(entry->entryNumber);
            Vars.var_BFRAG_PETERSON      = entry->analysisBase->weightBFrag_Peterson(entry->entryNumber);
            Vars.var_BSEMILEP_UP         = entry->analysisBase->weightBSemilep_Up(entry->entryNumber);
            Vars.var_BSEMILEP_DOWN       = entry->analysisBase->weightBSemilep_Dn(entry->entryNumber);
            Vars.var_BTAG_BFRAG_UP       = weightBtagSF_BFRAG_Up;
            Vars.var_BTAG_BFRAG_DOWN     = weightBtagSF_BFRAG_Down;
            Vars.var_BTAG_BFRAG_PETERSON = weightBtagSF_BFRAG_Peterson;
            Vars.var_BTAG_BSEMILEP_UP    = weightBtagSF_BSEMILEP_Up;
            Vars.var_BTAG_BSEMILEP_DOWN  = weightBtagSF_BSEMILEP_Down;
        }

        // const std::vector<float>& pdfWeightsFloat = *(entry->analysisBase->pdfWeights_);
        // for(unsigned int i = 0; i < pdfWeightsFloat.size(); i++){
        //   Vars.var_PDF.push_back(pdfWeightsFloat[i]);
        // }
        if(entry->analysisBase->readInPSWeights_ && entry->analysisBase->nPSWeights_> 0){
            const std::vector<float>& psWeightsFloat = *(entry->analysisBase->psWeights_);
            for(unsigned int i = 0; i < psWeightsFloat.size(); i++){
              Vars.var_PS[Vars.n_PSweights]=(psWeightsFloat[i]);
              ++Vars.n_PSweights;
            }
            if(Vars.n_PSweights>10){
                Vars.var_PSSCALE_WEIGHT_4_UP         = Vars.var_PS[6];
                Vars.var_PSSCALE_WEIGHT_4_DOWN       = Vars.var_PS[8];
                Vars.var_PSSCALE_WEIGHT_5_UP         = Vars.var_PS[7];
                Vars.var_PSSCALE_WEIGHT_5_DOWN       = Vars.var_PS[9];
                Vars.var_BTAG_PSSCALE_WEIGHT_4_UP    = weightBtagSF_PSSCALE_WEIGHT_4_Up;
                Vars.var_BTAG_PSSCALE_WEIGHT_4_DOWN  = weightBtagSF_PSSCALE_WEIGHT_4_Down;
                Vars.var_BTAG_PSSCALE_WEIGHT_5_UP    = weightBtagSF_PSSCALE_WEIGHT_5_Up;
                Vars.var_BTAG_PSSCALE_WEIGHT_5_DOWN  = weightBtagSF_PSSCALE_WEIGHT_5_Down;
            }
        }
        // }
        if(Vars.passStep3){
            Vars.var_PU_UP            = entry->pileupW_Up;
            Vars.var_PU_DOWN          = entry->pileupW_Down;
            Vars.var_TRIG_UP          = weightTriggerSF_Up;
            Vars.var_TRIG_DOWN        = weightTriggerSF_Dn;
            Vars.var_ELE_ID_UP        = weightLeptonSF_Ele_ID_Up;
            Vars.var_ELE_ID_DOWN      = weightLeptonSF_Ele_ID_Dn;
            Vars.var_ELE_RECO_UP      = weightLeptonSF_Ele_Reco_Up;
            Vars.var_ELE_RECO_DOWN    = weightLeptonSF_Ele_Reco_Dn;
            Vars.var_MUON_ID_UP       = weightLeptonSF_Muon_ID_Up;
            Vars.var_MUON_ID_DOWN     = weightLeptonSF_Muon_ID_Dn;
            Vars.var_MUON_ISO_UP      = weightLeptonSF_Muon_Iso_Up;
            Vars.var_MUON_ISO_DOWN    = weightLeptonSF_Muon_Iso_Dn;
            Vars.var_KIN_UP           = weightKinReco_Up;
            Vars.var_KIN_DOWN         = weightKinReco_Dn;
            Vars.var_LOOSEKIN_UP      = 1.;
            Vars.var_LOOSEKIN_DOWN    = 1.;
            Vars.var_BTAG_UP          = weightBtagSF_Up;
            Vars.var_BTAG_DOWN        = weightBtagSF_Dn;
            Vars.var_BTAG_LJET_UP     = weightBtagSF_LjetUp;
            Vars.var_BTAG_LJET_DOWN   = weightBtagSF_LjetDn;
            Vars.var_L1PREFIRING_UP   = l1PrefiringWeight_Up;
            Vars.var_L1PREFIRING_DOWN = l1PrefiringWeight_Dn;

            Vars.var_BTAG_PU_UP            = weightBtagSF_PU_Up;
            Vars.var_BTAG_PU_DOWN          = weightBtagSF_PU_Down;
            Vars.var_BTAG_TRIG_UP          = weightBtagSF_TRIG_Up;
            Vars.var_BTAG_TRIG_DOWN        = weightBtagSF_TRIG_Down;
            Vars.var_BTAG_ELE_ID_UP        = weightBtagSF_ELE_ID_Up;
            Vars.var_BTAG_ELE_ID_DOWN      = weightBtagSF_ELE_ID_Down;
            Vars.var_BTAG_ELE_RECO_UP      = weightBtagSF_ELE_RECO_Up;
            Vars.var_BTAG_ELE_RECO_DOWN    = weightBtagSF_ELE_RECO_Down;
            Vars.var_BTAG_MUON_ID_UP       = weightBtagSF_MUON_ID_Up;
            Vars.var_BTAG_MUON_ID_DOWN     = weightBtagSF_MUON_ID_Down;
            Vars.var_BTAG_MUON_ISO_UP      = weightBtagSF_MUON_ISO_Up;
            Vars.var_BTAG_MUON_ISO_DOWN    = weightBtagSF_MUON_ISO_Down;
            Vars.var_BTAG_UP               = weightBtagSF_Up;
            Vars.var_BTAG_DOWN             = weightBtagSF_Dn;
            Vars.var_BTAG_LJET_UP          = weightBtagSF_LjetUp;
            Vars.var_BTAG_LJET_DOWN        = weightBtagSF_LjetDn;
            Vars.var_BTAG_L1PREFIRING_UP   = weightBtagSF_L1PREFIRING_Up;
            Vars.var_BTAG_L1PREFIRING_DOWN = weightBtagSF_L1PREFIRING_Down;


        }
    }








    if(Vars.passStep3){
        Vars.nVtx = float(recoObjects.vertMultiGood_);
        // Vars.nLeptons = int(recoObjectIndices.allLeptonIndices_.size());
        // Vars.nJets = int(recoObjectIndices.jetIndices_.size());
        // Vars.nBJets = int(recoObjectIndices.bjetIndices_.size());

        // lepton kinematics
        const auto& lepton_1 = recoObjects.allLeptons_->at(recoObjectIndices.leadingLeptonIndex_);
        const auto& lepton_2 = recoObjects.allLeptons_->at(recoObjectIndices.nLeadingLeptonIndex_);

        Vars.lepton1_pt = lepton_1.Pt();
        Vars.lepton1_eta = lepton_1.Eta();
        Vars.lepton1_phi = lepton_1.Phi();
        Vars.lepton1_m = lepton_1.M();

        Vars.lepton2_pt = lepton_2.Pt();
        Vars.lepton2_eta = lepton_2.Eta();
        Vars.lepton2_phi = lepton_2.Phi();
        Vars.lepton2_m = lepton_2.M();

        if(recoObjects.lepPdgId_)
        {
          Vars.lepton1_pdgID = recoObjects.lepPdgId_->at(recoObjectIndices.leadingLeptonIndex_);
          Vars.lepton2_pdgID = recoObjects.lepPdgId_->at(recoObjectIndices.nLeadingLeptonIndex_);
        }

        // MET kinematics
        if(recoObjects.met_)
        {
          Vars.met_pt  = met.Pt();
          Vars.met_phi = met.Phi();
          Vars.met_significance = recoObjects.metsignificance_;
        }

        // jet kinematics
        if(recoObjects.jets_){
            const auto jets_p4 = *(recoObjects.jets_);
            const auto& indices_jets = recoObjectIndices.jetIndices_;
            const auto& indices_btags = recoObjectIndices.bjetIndices_;
            for(uint idx=0; idx<indices_jets.size(); ++idx){
                Vars.jets_pt[Vars.n_jets]=(jets_p4.at(indices_jets.at(idx)).pt());
                Vars.jets_eta[Vars.n_jets]=(jets_p4.at(indices_jets.at(idx)).eta());
                Vars.jets_phi[Vars.n_jets]=(jets_p4.at(indices_jets.at(idx)).phi());
                Vars.jets_m[Vars.n_jets]=(jets_p4.at(indices_jets.at(idx)).M());
                Vars.jets_charge[Vars.n_jets]=(recoObjects.jetChargeRelativePtWeighted_->at(indices_jets.at(idx)));
                if(recoObjects.jetBtags_){
                    Vars.jets_btag[Vars.n_jets]=(( std::find(indices_btags.begin(), indices_btags.end(), indices_jets.at(idx)) != indices_btags.end() ));
                    if(( std::find(indices_btags.begin(), indices_btags.end(), indices_jets.at(idx)) != indices_btags.end() )){
                        Vars.bjets_pt[Vars.n_bjets]=(jets_p4.at(indices_jets.at(idx)).pt());
                        Vars.bjets_eta[Vars.n_bjets]=(jets_p4.at(indices_jets.at(idx)).eta());
                        Vars.bjets_phi[Vars.n_bjets]=(jets_p4.at(indices_jets.at(idx)).phi());
                        Vars.bjets_m[Vars.n_bjets]=(jets_p4.at(indices_jets.at(idx)).M());
                        Vars.bjets_charge[Vars.n_bjets]=(recoObjects.jetChargeRelativePtWeighted_->at(indices_jets.at(idx)));
                        ++Vars.n_bjets;
                    }else{
                        Vars.nonbjets_pt[Vars.n_nonbjets]=(jets_p4.at(indices_jets.at(idx)).pt());
                        Vars.nonbjets_eta[Vars.n_nonbjets]=(jets_p4.at(indices_jets.at(idx)).eta());
                        Vars.nonbjets_phi[Vars.n_nonbjets]=(jets_p4.at(indices_jets.at(idx)).phi());
                        Vars.nonbjets_m[Vars.n_nonbjets]=(jets_p4.at(indices_jets.at(idx)).M());
                        Vars.nonbjets_charge[Vars.n_nonbjets]=(recoObjects.jetChargeRelativePtWeighted_->at(indices_jets.at(idx)));
                        ++Vars.n_nonbjets;
                    }
                }
                ++Vars.n_jets;
            }

            if(isTTbarSample_){
                if(entry->topGenObjects->valuesSet_){

                    // save if jets are matched to tops/ leading add jet
                    const auto& indices_jets = recoObjectIndices.jetIndices_;
                    for(uint idx=0; idx<indices_jets.size(); ++idx){
                        if(((int) indices_jets.at(idx)==matchedBjetFromTopIndex)||((int) indices_jets.at(idx)==matchedAntiBjetFromTopIndex)){
                            // Vars.jets_topMatched.push_back(true);
                            Vars.jets_topMatched[idx]=(true);
                        }else{
                            Vars.jets_topMatched[idx]=(false);
                        }
                        if(genExtraJetIndices.size()>0){
                          if((uint)indices_jets.at(idx)==(uint)genExtraJetIndices.at(0)){
                              Vars.jets_addJetMatched[idx]=(true);
                          }else{
                              Vars.jets_addJetMatched[idx]=(false);
                          }
                        }else{
                            Vars.jets_addJetMatched[idx]=(false);
                        }
                        if(genExtraJetWithNuIndices.size()>0){
                          if((uint)indices_jets.at(idx)==(uint)genExtraJetWithNuIndices.at(0)){
                              Vars.jets_addJetWithNuMatched[idx]=(true);
                          }else{
                              Vars.jets_addJetWithNuMatched[idx]=(false);
                          }
                        }else{
                            Vars.jets_addJetWithNuMatched[idx]=(false);
                        }
                    }

                }
            }
            if(Vars.n_gen_additional_jets>0){
                      TLorentzVector tempRJ(0.,0.,0.,0.);
                      tempRJ.SetPtEtaPhiM(Vars.gen_additional_jets_pt[0],Vars.gen_additional_jets_eta[0],Vars.gen_additional_jets_phi[0],Vars.gen_additional_jets_m[0]);
                      Vars.gen_rho =  340./(common::LVtoTLV(*entry->topGenObjects->GenTop_)+common::LVtoTLV(*entry->topGenObjects->GenAntiTop_)+tempRJ).M();
            }else{
                Vars.gen_rho =  -999.;
            }
            if(Vars.n_gen_additional_jets_withNu>0){
                      TLorentzVector tempRJ(0.,0.,0.,0.);
                      tempRJ.SetPtEtaPhiM(Vars.gen_additional_jets_withNu_pt[0],Vars.gen_additional_jets_withNu_eta[0],Vars.gen_additional_jets_withNu_phi[0],Vars.gen_additional_jets_withNu_m[0]);
                      Vars.gen_rhoWithNu =  340./(common::LVtoTLV(*entry->topGenObjects->GenTop_)+common::LVtoTLV(*entry->topGenObjects->GenAntiTop_)+tempRJ).M();
            }else{
                Vars.gen_rhoWithNu =  -999.;
            }
        }

        const KinematicReconstructionSolutions kinematicReconstructionSolutions = entry->analysisBase->kinematicReconstructionSolutions(leptonIndex, antiLeptonIndex, jetIndices, bjetIndices, allLeptons, jets, jetBtags, met);
        const bool hasSolution = kinematicReconstructionSolutions.numberOfSolutions();
        Vars.hasKinRecoSolution = hasSolution;

            if(Vars.hasKinRecoSolution){
                KinematicReconstructionSolution kr = kinematicReconstructionSolutions.solution();
                const LV& hypTop = kinematicReconstructionSolutions.solution().top();
                const LV& hypAntiTop = kinematicReconstructionSolutions.solution().antiTop();
                const LV& hypLepton = kinematicReconstructionSolutions.solution().lepton();
                const LV& hypAntiLepton = kinematicReconstructionSolutions.solution().antiLepton();
                const LV& hypBjet = kinematicReconstructionSolutions.solution().bjet();
                const LV& hypAntiBjet = kinematicReconstructionSolutions.solution().antiBjet();

                Vars.kinReco_top_pt = hypTop.Pt();
                Vars.kinReco_top_eta = hypTop.Eta();
                Vars.kinReco_top_phi = hypTop.Phi();
                Vars.kinReco_top_m = hypTop.M();
                Vars.kinReco_antitop_pt = hypAntiTop.Pt();
                Vars.kinReco_antitop_eta = hypAntiTop.Eta();
                Vars.kinReco_antitop_phi = hypAntiTop.Phi();
                Vars.kinReco_antitop_m = hypAntiTop.M();
                Vars.kinReco_lepton_pt = hypLepton.Pt();
                Vars.kinReco_lepton_eta = hypLepton.Eta();
                Vars.kinReco_lepton_phi = hypLepton.Phi();
                Vars.kinReco_lepton_m = hypLepton.M();
                Vars.kinReco_antilepton_pt = hypAntiLepton.Pt();
                Vars.kinReco_antilepton_eta = hypAntiLepton.Eta();
                Vars.kinReco_antilepton_phi = hypAntiLepton.Phi();
                Vars.kinReco_antilepton_m = hypAntiLepton.M();
                Vars.kinReco_b_pt = hypBjet.Pt();
                Vars.kinReco_b_eta = hypBjet.Eta();
                Vars.kinReco_b_phi = hypBjet.Phi();
                Vars.kinReco_b_m = hypBjet.M();
                Vars.kinReco_antib_pt = hypAntiBjet.Pt();
                Vars.kinReco_antib_eta = hypAntiBjet.Eta();
                Vars.kinReco_antib_phi = hypAntiBjet.Phi();
                Vars.kinReco_antib_m = hypAntiBjet.M();
                Vars.kinReco_weight = weightKinReco;

                const auto jets_p4 = *(recoObjects.jets_);
                const auto& indices_jets = recoObjectIndices.jetIndices_;
                // const auto& indices_btags = recoObjectIndices.bjetIndices_;
                for(uint idx=0; idx<indices_jets.size(); ++idx){
                    // if(jets_p4.at(indices_jets.at(idx)).pt()!=Vars.kinReco_b_pt){
                    if(abs((jets_p4.at(indices_jets.at(idx)).pt())-(Vars.kinReco_b_pt))>0.01){
                        // if(jets_p4.at(indices_jets.at(idx)).pt()!=Vars.kinReco_antib_pt){
                        if(abs((jets_p4.at(indices_jets.at(idx)).pt())-(Vars.kinReco_antib_pt))>0.01){
                            Vars.kinReco_nonbjets_pt[Vars.n_kinReco_nonbjets]=(jets_p4.at(indices_jets.at(idx)).pt());
                            Vars.kinReco_nonbjets_eta[Vars.n_kinReco_nonbjets]=(jets_p4.at(indices_jets.at(idx)).eta());
                            Vars.kinReco_nonbjets_phi[Vars.n_kinReco_nonbjets]=(jets_p4.at(indices_jets.at(idx)).phi());
                            Vars.kinReco_nonbjets_m[Vars.n_kinReco_nonbjets]=(jets_p4.at(indices_jets.at(idx)).M());
                            ++Vars.n_kinReco_nonbjets;
                        }
                    }
                }

                if(Vars.n_kinReco_nonbjets>0){
                    Vars.kinReco_nonbjet_pt  = Vars.kinReco_nonbjets_pt[0];
                    Vars.kinReco_nonbjet_eta = Vars.kinReco_nonbjets_eta[0];
                    Vars.kinReco_nonbjet_phi = Vars.kinReco_nonbjets_phi[0];
                    Vars.kinReco_nonbjet_m   = Vars.kinReco_nonbjets_m[0];
                    if(Vars.kinReco_nonbjet_pt>20.){
                        if(abs(Vars.kinReco_nonbjet_eta)<2.5){
                            TLorentzVector tempRJ(0.,0.,0.,0.);
                            tempRJ.SetPtEtaPhiM(Vars.kinReco_nonbjet_pt,Vars.kinReco_nonbjet_eta,Vars.kinReco_nonbjet_phi,Vars.kinReco_nonbjet_m);
                            Vars.kinReco_rho =  340./(common::LVtoTLV(hypTop)+common::LVtoTLV(hypAntiTop)+tempRJ).M();
                        }
                    }
                }
          }

          const LooseKinRecoSolution looseKinematicReconstructionSolution = entry->analysisBase->looseKinRecoSolutions(leptonIndex, antiLeptonIndex, jetIndices, bjetIndices, allLeptons, jets, met);
          const bool hasLooseSolution = looseKinematicReconstructionSolution.hasSolution();
          Vars.hasLooseKinRecoSolution = hasLooseSolution;
          Vars.passStep8 = Vars.passStep3 && Vars.passStep4 && Vars.passStep5 && Vars.passStep6 && Vars.passStep7 && Vars.hasKinRecoSolution;
          Vars.passStep8Loose = Vars.passStep3 && Vars.passStep4 && Vars.passStep5 && Vars.passStep6 && Vars.passStep7 && Vars.hasLooseKinRecoSolution;

          if(Vars.hasLooseKinRecoSolution){

              Vars.looseKinReco_ttbar_pt = looseKinematicReconstructionSolution.TTbar().Pt();
              Vars.looseKinReco_ttbar_eta = looseKinematicReconstructionSolution.TTbar().Eta();
              Vars.looseKinReco_ttbar_phi = looseKinematicReconstructionSolution.TTbar().Phi();
              Vars.looseKinReco_ttbar_m = looseKinematicReconstructionSolution.TTbar().M();
              Vars.looseKinReco_b1_pt = looseKinematicReconstructionSolution.getLVBBbar().first.Pt();
              Vars.looseKinReco_b1_eta = looseKinematicReconstructionSolution.getLVBBbar().first.Eta();
              Vars.looseKinReco_b1_phi = looseKinematicReconstructionSolution.getLVBBbar().first.Phi();
              Vars.looseKinReco_b1_m = looseKinematicReconstructionSolution.getLVBBbar().first.M();
              Vars.looseKinReco_b2_pt = looseKinematicReconstructionSolution.getLVBBbar().second.Pt();
              Vars.looseKinReco_b2_eta = looseKinematicReconstructionSolution.getLVBBbar().second.Eta();
              Vars.looseKinReco_b2_phi = looseKinematicReconstructionSolution.getLVBBbar().second.Phi();
              Vars.looseKinReco_b2_m = looseKinematicReconstructionSolution.getLVBBbar().second.M();
              Vars.looseKinReco_l1_pt = looseKinematicReconstructionSolution.getLVLLbar().first.Pt();
              Vars.looseKinReco_l1_eta = looseKinematicReconstructionSolution.getLVLLbar().first.Eta();
              Vars.looseKinReco_l1_phi = looseKinematicReconstructionSolution.getLVLLbar().first.Phi();
              Vars.looseKinReco_l1_m = looseKinematicReconstructionSolution.getLVLLbar().first.M();
              Vars.looseKinReco_l2_pt = looseKinematicReconstructionSolution.getLVLLbar().second.Pt();
              Vars.looseKinReco_l2_eta = looseKinematicReconstructionSolution.getLVLLbar().second.Eta();
              Vars.looseKinReco_l2_phi = looseKinematicReconstructionSolution.getLVLLbar().second.Phi();
              Vars.looseKinReco_l2_m = looseKinematicReconstructionSolution.getLVLLbar().second.M();
              Vars.looseKinReco_weight = weightLooseKinReco;

              const LV& lhypTTbar = looseKinematicReconstructionSolution.TTbar();

              const auto jets_p4 = *(recoObjects.jets_);
              const auto& indices_jets = recoObjectIndices.jetIndices_;
              // const auto& indices_btags = recoObjectIndices.bjetIndices_;
              for(uint idx=0; idx<indices_jets.size(); ++idx){
                  // if(jets_p4.at(indices_jets.at(idx)).pt()!=Vars.looseKinReco_b1_pt){
                  if(abs((jets_p4.at(indices_jets.at(idx)).pt())-(Vars.looseKinReco_b1_pt))>0.01){
                      // if(jets_p4.at(indices_jets.at(idx)).pt()!=Vars.looseKinReco_b2_pt){
                      if(abs((jets_p4.at(indices_jets.at(idx)).pt())-(Vars.looseKinReco_b2_pt))>0.01){
                          Vars.looseKinReco_nonbjets_pt[Vars.n_looseKinReco_nonbjets]=(jets_p4.at(indices_jets.at(idx)).pt());
                          Vars.looseKinReco_nonbjets_eta[Vars.n_looseKinReco_nonbjets]=(jets_p4.at(indices_jets.at(idx)).eta());
                          Vars.looseKinReco_nonbjets_phi[Vars.n_looseKinReco_nonbjets]=(jets_p4.at(indices_jets.at(idx)).phi());
                          Vars.looseKinReco_nonbjets_m[Vars.n_looseKinReco_nonbjets]=(jets_p4.at(indices_jets.at(idx)).M());
                          ++Vars.n_looseKinReco_nonbjets;
                      }
                  }
              }
              if(Vars.n_looseKinReco_nonbjets>0){
                  Vars.looseKinReco_nonbjet_pt  = Vars.looseKinReco_nonbjets_pt[0];
                  Vars.looseKinReco_nonbjet_eta = Vars.looseKinReco_nonbjets_eta[0];
                  Vars.looseKinReco_nonbjet_phi = Vars.looseKinReco_nonbjets_phi[0];
                  Vars.looseKinReco_nonbjet_m   = Vars.looseKinReco_nonbjets_m[0];
                  if(Vars.looseKinReco_nonbjet_pt>20.){
                      if(abs(Vars.looseKinReco_nonbjet_eta)<2.5){
                          TLorentzVector tempRJ(0.,0.,0.,0.);
                          tempRJ.SetPtEtaPhiM(Vars.looseKinReco_nonbjet_pt, Vars.looseKinReco_nonbjet_eta, Vars.looseKinReco_nonbjet_phi, Vars.looseKinReco_nonbjet_m);
                          Vars.looseKinReco_rho =  340./(common::LVtoTLV(lhypTTbar)+tempRJ).M();
                      }
                  }
              }
          }


    }



    TTreePtr->Fill();
}






void MiniTree::ResetBranches()
{
  // Vars.eventCounter = -999;
  Vars.eventNumber       = 0;
  Vars.lumiNumber        = 0;
  Vars.runNumber         = 0;
  Vars.weight            = 1.;
  Vars.leptonSF          = 1.;
  Vars.l1PrefiringWeight = 1.;
  Vars.triggerSF         = 1.;
  Vars.btagSF            = 1.;
  Vars.pileupSF          = 1.;
  Vars.channel           = -999.;
  Vars.nVtx              = -999.;
  Vars.nLeptons          = -999.;
  Vars.nJets             = -999.;
  Vars.nBJets            = -999.;
  Vars.lepton1_pt        = -999.;
  Vars.lepton1_eta       = -999.;
  Vars.lepton1_phi       = -999.;
  Vars.lepton1_m         = -999.;
  Vars.lepton2_pt        = -999.;
  Vars.lepton2_eta       = -999.;
  Vars.lepton2_phi       = -999.;
  Vars.lepton2_m         = -999.;
  Vars.lepton1_pdgID     = -999.;
  Vars.lepton2_pdgID     = -999.;
  Vars.met_pt            = -999.;
  Vars.met_phi           = -999.;
  Vars.met_significance  = -999.;

  Vars.n_jets = 0;
  std::fill( std::begin(Vars.jets_pt), std::end(Vars.jets_pt), 0 );
  std::fill( std::begin(Vars.jets_eta), std::end(Vars.jets_eta), 0 );
  std::fill( std::begin(Vars.jets_phi), std::end(Vars.jets_phi), 0 );
  std::fill( std::begin(Vars.jets_m), std::end(Vars.jets_m), 0 );
  std::fill( std::begin(Vars.jets_charge), std::end(Vars.jets_charge), 0 );
  std::fill( std::begin(Vars.jets_btag), std::end(Vars.jets_btag), 0 );
  std::fill( std::begin(Vars.jets_topMatched), std::end(Vars.jets_topMatched), 0 );
  std::fill( std::begin(Vars.jets_addJetMatched), std::end(Vars.jets_addJetMatched), 0 );

  Vars.n_bjets = 0;
  std::fill( std::begin(Vars.bjets_pt), std::end(Vars.bjets_pt), 0 );
  std::fill( std::begin(Vars.bjets_eta), std::end(Vars.bjets_eta), 0 );
  std::fill( std::begin(Vars.bjets_phi), std::end(Vars.bjets_phi), 0 );
  std::fill( std::begin(Vars.bjets_m), std::end(Vars.bjets_m), 0 );
  std::fill( std::begin(Vars.bjets_charge), std::end(Vars.bjets_charge), 0 );
  Vars.n_nonbjets = 0;
  std::fill( std::begin(Vars.nonbjets_pt), std::end(Vars.nonbjets_pt), 0 );
  std::fill( std::begin(Vars.nonbjets_eta), std::end(Vars.nonbjets_eta), 0 );
  std::fill( std::begin(Vars.nonbjets_phi), std::end(Vars.nonbjets_phi), 0 );
  std::fill( std::begin(Vars.nonbjets_m), std::end(Vars.nonbjets_m), 0 );
  std::fill( std::begin(Vars.nonbjets_charge), std::end(Vars.nonbjets_charge), 0 );

  Vars.passStep3      = false;
  Vars.passStep4      = false;
  Vars.passStep5      = false;
  Vars.passStep6      = false;
  Vars.passStep7      = false;
  Vars.passStep8      = false;
  Vars.passStep8Loose = false;

  Vars.hasKinRecoSolution     = false;
  Vars.kinReco_top_pt         = -999.;
  Vars.kinReco_top_eta        = -999.;
  Vars.kinReco_top_phi        = -999.;
  Vars.kinReco_top_m          = -999.;
  Vars.kinReco_antitop_pt     = -999.;
  Vars.kinReco_antitop_eta    = -999.;
  Vars.kinReco_antitop_phi    = -999.;
  Vars.kinReco_antitop_m      = -999.;
  Vars.kinReco_lepton_pt      = -999.;
  Vars.kinReco_lepton_eta     = -999.;
  Vars.kinReco_lepton_phi     = -999.;
  Vars.kinReco_lepton_m       = -999.;
  Vars.kinReco_antilepton_pt  = -999.;
  Vars.kinReco_antilepton_eta = -999.;
  Vars.kinReco_antilepton_phi = -999.;
  Vars.kinReco_antilepton_m   = -999.;
  Vars.kinReco_b_pt           = -999.;
  Vars.kinReco_b_eta          = -999.;
  Vars.kinReco_b_phi          = -999.;
  Vars.kinReco_b_m            = -999.;
  Vars.kinReco_antib_pt       = -999.;
  Vars.kinReco_antib_eta      = -999.;
  Vars.kinReco_antib_phi      = -999.;
  Vars.kinReco_antib_m        = -999.;
  Vars.kinReco_weight         = 1.;
  Vars.n_kinReco_nonbjets     = 0;
  std::fill( std::begin(Vars.kinReco_nonbjets_pt), std::end(Vars.kinReco_nonbjets_pt), 0 );
  std::fill( std::begin(Vars.kinReco_nonbjets_eta), std::end(Vars.kinReco_nonbjets_eta), 0 );
  std::fill( std::begin(Vars.kinReco_nonbjets_phi), std::end(Vars.kinReco_nonbjets_phi), 0 );
  std::fill( std::begin(Vars.kinReco_nonbjets_m), std::end(Vars.kinReco_nonbjets_m), 0 );
  Vars.kinReco_nonbjet_pt     = -999.;
  Vars.kinReco_nonbjet_eta    = -999.;
  Vars.kinReco_nonbjet_phi    = -999.;
  Vars.kinReco_nonbjet_m      = -999.;

  Vars.kinReco_rho = -999.;

  Vars.hasLooseKinRecoSolution = false;

  Vars.looseKinReco_ttbar_pt   = -999.;
  Vars.looseKinReco_ttbar_eta  = -999.;
  Vars.looseKinReco_ttbar_phi  = -999.;
  Vars.looseKinReco_ttbar_m    = -999.;
  Vars.looseKinReco_b1_pt      = -999.;
  Vars.looseKinReco_b1_eta     = -999.;
  Vars.looseKinReco_b1_phi     = -999.;
  Vars.looseKinReco_b1_m       = -999.;
  Vars.looseKinReco_b2_pt      = -999.;
  Vars.looseKinReco_b2_eta     = -999.;
  Vars.looseKinReco_b2_phi     = -999.;
  Vars.looseKinReco_b2_m       = -999.;
  Vars.looseKinReco_l1_pt      = -999.;
  Vars.looseKinReco_l1_eta     = -999.;
  Vars.looseKinReco_l1_phi     = -999.;
  Vars.looseKinReco_l1_m       = -999.;
  Vars.looseKinReco_l2_pt      = -999.;
  Vars.looseKinReco_l2_eta     = -999.;
  Vars.looseKinReco_l2_phi     = -999.;
  Vars.looseKinReco_l2_m       = -999.;
  Vars.looseKinReco_weight     = 1.;
  Vars.n_looseKinReco_nonbjets = 0;
  std::fill( std::begin(Vars.looseKinReco_nonbjets_pt), std::end(Vars.looseKinReco_nonbjets_pt), 0 );
  std::fill( std::begin(Vars.looseKinReco_nonbjets_eta), std::end(Vars.looseKinReco_nonbjets_eta), 0 );
  std::fill( std::begin(Vars.looseKinReco_nonbjets_phi), std::end(Vars.looseKinReco_nonbjets_phi), 0 );
  std::fill( std::begin(Vars.looseKinReco_nonbjets_m), std::end(Vars.looseKinReco_nonbjets_m), 0 );
  Vars.looseKinReco_nonbjet_pt  = -999.;
  Vars.looseKinReco_nonbjet_eta = -999.;
  Vars.looseKinReco_nonbjet_phi = -999.;
  Vars.looseKinReco_nonbjet_m   = -999.;

  Vars.looseKinReco_rho = -999.;

  Vars.n_gen_additional_jets = 0;
  std::fill( std::begin(Vars.gen_additional_jets_pt), std::end(Vars.gen_additional_jets_pt), 0 );
  std::fill( std::begin(Vars.gen_additional_jets_eta), std::end(Vars.gen_additional_jets_eta), 0 );
  std::fill( std::begin(Vars.gen_additional_jets_phi), std::end(Vars.gen_additional_jets_phi), 0 );
  std::fill( std::begin(Vars.gen_additional_jets_m), std::end(Vars.gen_additional_jets_m), 0 );

  Vars.n_gen_additional_jets_withNu = 0;
  std::fill( std::begin(Vars.gen_additional_jets_withNu_pt), std::end(Vars.gen_additional_jets_withNu_pt), 0 );
  std::fill( std::begin(Vars.gen_additional_jets_withNu_eta), std::end(Vars.gen_additional_jets_withNu_eta), 0 );
  std::fill( std::begin(Vars.gen_additional_jets_withNu_phi), std::end(Vars.gen_additional_jets_withNu_phi), 0 );
  std::fill( std::begin(Vars.gen_additional_jets_withNu_m), std::end(Vars.gen_additional_jets_withNu_m), 0 );

  Vars.gen_additional_jet_pt  = -999.;
  Vars.gen_additional_jet_eta = -999.;
  Vars.gen_additional_jet_phi = -999.;
  Vars.gen_additional_jet_m   = -999.;

  Vars.gen_additional_jet_withNu_pt  = -999.;
  Vars.gen_additional_jet_withNu_eta = -999.;
  Vars.gen_additional_jet_withNu_phi = -999.;
  Vars.gen_additional_jet_withNu_m   = -999.;

  Vars.gen_top_pt           = -999.;
  Vars.gen_top_eta          = -999.;
  Vars.gen_top_phi          = -999.;
  Vars.gen_top_m            = -999.;
  Vars.gen_antitop_pt       = -999.;
  Vars.gen_antitop_eta      = -999.;
  Vars.gen_antitop_phi      = -999.;
  Vars.gen_antitop_m        = -999.;
  Vars.gen_lepton_pt        = -999.;
  Vars.gen_lepton_eta       = -999.;
  Vars.gen_lepton_phi       = -999.;
  Vars.gen_lepton_m         = -999.;
  Vars.gen_antilepton_pt    = -999.;
  Vars.gen_antilepton_eta   = -999.;
  Vars.gen_antilepton_phi   = -999.;
  Vars.gen_antilepton_m     = -999.;
  Vars.gen_b_pt             = -999.;
  Vars.gen_b_eta            = -999.;
  Vars.gen_b_phi            = -999.;
  Vars.gen_b_m              = -999.;
  Vars.gen_antib_pt         = -999.;
  Vars.gen_antib_eta        = -999.;
  Vars.gen_antib_phi        = -999.;
  Vars.gen_antib_m          = -999.;
  Vars.gen_neutrino_pt      = -999.;
  Vars.gen_neutrino_eta     = -999.;
  Vars.gen_neutrino_phi     = -999.;
  Vars.gen_neutrino_m       = -999.;
  Vars.gen_antineutrino_pt  = -999.;
  Vars.gen_antineutrino_eta = -999.;
  Vars.gen_antineutrino_phi = -999.;
  Vars.gen_antineutrino_m   = -999.;
  Vars.gen_nVtx             = 0;
  Vars.gen_nJets            = 0;
  Vars.gen_nBJets           = 0;
  Vars.gen_met_pt           = -999.;
  Vars.gen_met_eta          = -999.;
  Vars.gen_met_phi          = -999.;
  Vars.gen_met_m            = -999.;
  Vars.gen_rho              = -999.;
  Vars.gen_rhoWithNu        = -999.;

  Vars.pseudo_top_pt         = -999.;
  Vars.pseudo_top_eta        = -999.;
  Vars.pseudo_top_phi        = -999.;
  Vars.pseudo_top_m          = -999.;
  Vars.pseudo_antitop_pt     = -999.;
  Vars.pseudo_antitop_eta    = -999.;
  Vars.pseudo_antitop_phi    = -999.;
  Vars.pseudo_antitop_m      = -999.;
  Vars.pseudo_lepton_pt      = -999.;
  Vars.pseudo_lepton_eta     = -999.;
  Vars.pseudo_lepton_phi     = -999.;
  Vars.pseudo_lepton_m       = -999.;
  Vars.pseudo_antilepton_pt  = -999.;
  Vars.pseudo_antilepton_eta = -999.;
  Vars.pseudo_antilepton_phi = -999.;
  Vars.pseudo_antilepton_m   = -999.;
  Vars.pseudo_b_pt           = -999.;
  Vars.pseudo_b_eta          = -999.;
  Vars.pseudo_b_phi          = -999.;
  Vars.pseudo_b_m            = -999.;
  Vars.pseudo_antib_pt       = -999.;
  Vars.pseudo_antib_eta      = -999.;
  Vars.pseudo_antib_phi      = -999.;
  Vars.pseudo_antib_m        = -999.;

  Vars.n_pseudo_additional_jets = 0;
  std::fill( std::begin(Vars.pseudo_additional_jets_pt), std::end(Vars.pseudo_additional_jets_pt), 0 );
  std::fill( std::begin(Vars.pseudo_additional_jets_eta), std::end(Vars.pseudo_additional_jets_eta), 0 );
  std::fill( std::begin(Vars.pseudo_additional_jets_phi), std::end(Vars.pseudo_additional_jets_phi), 0 );
  std::fill( std::begin(Vars.pseudo_additional_jets_m), std::end(Vars.pseudo_additional_jets_m), 0 );

  Vars.pseudo_additional_jet_pt  = -999.;
  Vars.pseudo_additional_jet_eta = -999.;
  Vars.pseudo_additional_jet_phi = -999.;
  Vars.pseudo_additional_jet_m   = -999.;

  Vars.pseudo_rho = -999.;

  Vars.n_gen_partonLevel_additional_jets = 0;
  std::fill( std::begin(Vars.gen_partonLevel_additional_jets_pt), std::end(Vars.gen_partonLevel_additional_jets_pt), 0 );
  std::fill( std::begin(Vars.gen_partonLevel_additional_jets_eta), std::end(Vars.gen_partonLevel_additional_jets_eta), 0 );
  std::fill( std::begin(Vars.gen_partonLevel_additional_jets_phi), std::end(Vars.gen_partonLevel_additional_jets_phi), 0 );
  std::fill( std::begin(Vars.gen_partonLevel_additional_jets_m), std::end(Vars.gen_partonLevel_additional_jets_m), 0 );

  Vars.gen_partonLevel_additional_jet_pt  = -999.;
  Vars.gen_partonLevel_additional_jet_eta = -999.;
  Vars.gen_partonLevel_additional_jet_phi = -999.;
  Vars.gen_partonLevel_additional_jet_m   = -999.;

  Vars.gen_partonLevel_rho = -999.;

  Vars.var_MERENSCALE_UP   = -999.;
  Vars.var_MERENSCALE_DOWN = -999.;
  Vars.var_MEFACSCALE_UP   = -999.;
  Vars.var_MEFACSCALE_DOWN = -999.;
  Vars.var_MESCALE_UP      = -999.;
  Vars.var_MESCALE_DOWN    = -999.;
  Vars.var_BFRAG_UP        = -999.;
  Vars.var_BFRAG_DOWN      = -999.;
  Vars.var_BFRAG_PETERSON  = -999.;
  Vars.var_BSEMILEP_UP     = -999.;
  Vars.var_BSEMILEP_DOWN   = -999.;
  Vars.var_PDF_ALPHAS_UP   = -999.;
  Vars.var_PDF_ALPHAS_DOWN = -999.;

  Vars.n_PSweights = 0;
  // Vars.var_PS.clear();
  std::fill( std::begin(Vars.var_PS), std::end(Vars.var_PS), 0 );
  Vars.n_PDFweights = 0;
  // Vars.var_PDF.clear();
  std::fill( std::begin(Vars.var_PDF), std::end(Vars.var_PDF), 0 );

  Vars.var_PSSCALE_WEIGHT_4_UP   = -999.;
  Vars.var_PSSCALE_WEIGHT_4_DOWN = -999.;
  Vars.var_PSSCALE_WEIGHT_5_UP   = -999.;
  Vars.var_PSSCALE_WEIGHT_5_DOWN = -999.;

  Vars.var_PU_UP            = -999.;
  Vars.var_PU_DOWN          = -999.;
  Vars.var_TRIG_UP          = -999.;
  Vars.var_TRIG_DOWN        = -999.;
  Vars.var_ELE_ID_UP        = -999.;
  Vars.var_ELE_ID_DOWN      = -999.;
  Vars.var_ELE_RECO_UP      = -999.;
  Vars.var_ELE_RECO_DOWN    = -999.;
  Vars.var_MUON_ID_UP       = -999.;
  Vars.var_MUON_ID_DOWN     = -999.;
  Vars.var_MUON_ISO_UP      = -999.;
  Vars.var_MUON_ISO_DOWN    = -999.;
  Vars.var_KIN_UP           = -999.;
  Vars.var_KIN_DOWN         = -999.;
  Vars.var_LOOSEKIN_UP      = -999.;
  Vars.var_LOOSEKIN_DOWN    = -999.;
  Vars.var_BTAG_UP          = -999.;
  Vars.var_BTAG_DOWN        = -999.;
  Vars.var_BTAG_LJET_UP     = -999.;
  Vars.var_BTAG_LJET_DOWN   = -999.;
  Vars.var_L1PREFIRING_UP   = -999.;
  Vars.var_L1PREFIRING_DOWN = -999.;

  Vars.var_BTAG_MERENSCALE_UP         = -999.;
  Vars.var_BTAG_MERENSCALE_DOWN       = -999.;
  Vars.var_BTAG_MEFACSCALE_UP         = -999.;
  Vars.var_BTAG_MEFACSCALE_DOWN       = -999.;
  Vars.var_BTAG_MESCALE_UP            = -999.;
  Vars.var_BTAG_MESCALE_DOWN          = -999.;
  Vars.var_BTAG_BFRAG_UP              = -999.;
  Vars.var_BTAG_BFRAG_DOWN            = -999.;
  Vars.var_BTAG_BFRAG_PETERSON        = -999.;
  Vars.var_BTAG_BSEMILEP_UP           = -999.;
  Vars.var_BTAG_BSEMILEP_DOWN         = -999.;
  Vars.var_BTAG_PDF_ALPHAS_UP         = -999.;
  Vars.var_BTAG_PDF_ALPHAS_DOWN       = -999.;
  Vars.var_BTAG_PSSCALE_WEIGHT_4_UP   = -999.;
  Vars.var_BTAG_PSSCALE_WEIGHT_4_DOWN = -999.;
  Vars.var_BTAG_PSSCALE_WEIGHT_5_UP   = -999.;
  Vars.var_BTAG_PSSCALE_WEIGHT_5_DOWN = -999.;
  Vars.var_BTAG_PU_UP                 = -999.;
  Vars.var_BTAG_PU_DOWN               = -999.;
  Vars.var_BTAG_TRIG_UP               = -999.;
  Vars.var_BTAG_TRIG_DOWN             = -999.;
  Vars.var_BTAG_ELE_ID_UP             = -999.;
  Vars.var_BTAG_ELE_ID_DOWN           = -999.;
  Vars.var_BTAG_ELE_RECO_UP           = -999.;
  Vars.var_BTAG_ELE_RECO_DOWN         = -999.;
  Vars.var_BTAG_MUON_ID_UP            = -999.;
  Vars.var_BTAG_MUON_ID_DOWN          = -999.;
  Vars.var_BTAG_MUON_ISO_UP           = -999.;
  Vars.var_BTAG_MUON_ISO_DOWN         = -999.;
  Vars.var_BTAG_L1PREFIRING_UP        = -999.;
  Vars.var_BTAG_L1PREFIRING_DOWN      = -999.;
}



void MiniTree::generatorExtraJets(const TopGenObjects& topGenObjects,std::vector<int>& genExtraJetIndices, const int bHadronIndex, const int antiBHadronIndex, const bool withNu)
{
    if(topGenObjects.valuesSet_)
    {
        const int BHadronIndex(bHadronIndex);
        const int AntiBHadronIndex(antiBHadronIndex);

        if(!withNu){
            for(int genJet=0; genJet<(int)(*topGenObjects.allGenJets_).size(); genJet++)
            {
                if(BHadronIndex > -1 && AntiBHadronIndex > -1 && (*topGenObjects.allGenJets_).at(BHadronIndex) != (*topGenObjects.allGenJets_).at(genJet)  && (*topGenObjects.allGenJets_).at(AntiBHadronIndex) != (*topGenObjects.allGenJets_).at(genJet))
               {
                    genExtraJetIndices.push_back(genJet);
                }
            }
        }else{
            for(int genJet=0; genJet<(int)(*topGenObjects.allGenJetsWithNu_).size(); genJet++)
            {
                if(BHadronIndex > -1 && AntiBHadronIndex > -1 && (*topGenObjects.allGenJetsWithNu_).at(BHadronIndex) != (*topGenObjects.allGenJetsWithNu_).at(genJet)  && (*topGenObjects.allGenJetsWithNu_).at(AntiBHadronIndex) != (*topGenObjects.allGenJetsWithNu_).at(genJet))
               {
                    genExtraJetIndices.push_back(genJet);
                }
            }
        }
    }
}




void MiniTree::SetTriggerScaleFactors_Up(const TriggerScaleFactors* const scaleFactors)
{
    triggerScaleFactors_Up_ = scaleFactors;
}
void MiniTree::SetTriggerScaleFactors_Dn(const TriggerScaleFactors* const scaleFactors)
{
    triggerScaleFactors_Dn_ = scaleFactors;
}
void MiniTree::SetBtagScaleFactors_Var_AllSyst(const std::map<std::tuple<Systematic::Type, Systematic::Variation, int>, BtagScaleFactors*>* mapAllSystBTagScaleFactors)
{
    mapAllSystBTagScaleFactors_ = mapAllSystBTagScaleFactors;
}

  // kinematicReconstruction as dummy to keep same structure
void MiniTree::SetKinematicReconstruction_Up(const KinematicReconstruction* const kinematicReconstruction,
                                              const KinematicReconstructionScaleFactors* const kinematicReconstructionScaleFactors)
{
    (void) kinematicReconstruction;
    kinematicReconstructionScaleFactors_Up_ = kinematicReconstructionScaleFactors;
}
void MiniTree::SetKinematicReconstruction_Dn(const KinematicReconstruction* const kinematicReconstruction,
                                              const KinematicReconstructionScaleFactors* const kinematicReconstructionScaleFactors)
{
    (void) kinematicReconstruction;
    kinematicReconstructionScaleFactors_Dn_ = kinematicReconstructionScaleFactors;
}

  // LooseKinReco as dummy to keep same structure
void MiniTree::SetLooseKinReco_Up(const LooseKinReco* const looseKinReco,
					  const LooseKinRecoScaleFactors* const looseKinRecoScaleFactors)
{
    (void) looseKinReco;
    looseKinRecoScaleFactors_Up_ = looseKinRecoScaleFactors;
}
void MiniTree::SetLooseKinReco_Dn(const LooseKinReco* const looseKinReco,
					  const LooseKinRecoScaleFactors* const looseKinRecoScaleFactors)
{
    (void) looseKinReco;
    looseKinRecoScaleFactors_Dn_ = looseKinRecoScaleFactors;
}



float MiniTree::weightTriggerSF_Up(const int leptonXIndex, const int leptonYIndex,
                                     const VLV& allLeptons, const std::vector<int>& lepPdgId)const
{
    if(!isGen_) return 1.;
    if(leptonXIndex<0 || leptonYIndex<0) return 1.;
    if(!triggerScaleFactors_Up_) return 1.;
    return triggerScaleFactors_Up_->getSF(leptonXIndex, leptonYIndex, allLeptons, lepPdgId);
}

float MiniTree::weightTriggerSF_Dn(const int leptonXIndex, const int leptonYIndex,
                                     const VLV& allLeptons, const std::vector<int>& lepPdgId)const
{
    if(!isGen_) return 1.;
    if(leptonXIndex<0 || leptonYIndex<0) return 1.;
    if(!triggerScaleFactors_Dn_) return 1.;
    return triggerScaleFactors_Dn_->getSF(leptonXIndex, leptonYIndex, allLeptons, lepPdgId);
}

float MiniTree::weightBtagSF_Var(const std::vector<int>& jetIndices,
                                  const VLV& jets, const std::vector<int>& jetFlavour,
                                  const std::vector<float>& btagDiscriminators, const TString& sysKey)const
{
    if(!isGen_) return 1.;
    if(!mapAllSystBTagScaleFactors_) return 1.;
    Systematic::Systematic temporarySystematic = Systematic::Systematic(sysKey);
    std::tuple<Systematic::Type, Systematic::Variation, int> tempTuple{temporarySystematic.type(), temporarySystematic.variation(), temporarySystematic.variationNumber()};
    return mapAllSystBTagScaleFactors_->at(tempTuple)->getSF(jetIndices, jets, jetFlavour, btagDiscriminators);
}
// float MiniTree::weightBtagSF_Up(const std::vector<int>& jetIndices,
//                                   const VLV& jets, const std::vector<int>& jetFlavour,
//                                   const std::vector<float>& btagDiscriminators)const
// {
//     if(!isGen_) return 1.;
//     if(!btagScaleFactors_Up_) return 1.;
//     return btagScaleFactors_Up_->getSF(jetIndices, jets, jetFlavour, btagDiscriminators);
// }
// float MiniTree::weightBtagSF_Dn(const std::vector<int>& jetIndices,
//                                   const VLV& jets, const std::vector<int>& jetFlavour,
//                                   const std::vector<float>& btagDiscriminators)const
// {
//     if(!isGen_) return 1.;
//     if(!btagScaleFactors_Dn_) return 1.;
//     return btagScaleFactors_Dn_->getSF(jetIndices, jets, jetFlavour, btagDiscriminators);
// }
// float MiniTree::weightBtagSF_LjetUp(const std::vector<int>& jetIndices,
//                                   const VLV& jets, const std::vector<int>& jetFlavour,
//                                   const std::vector<float>& btagDiscriminators)const
// {
//     if(!isGen_) return 1.;
//     if(!btagScaleFactors_LjetUp_) return 1.;
//     return btagScaleFactors_LjetUp_->getSF(jetIndices, jets, jetFlavour, btagDiscriminators);
// }
// float MiniTree::weightBtagSF_LjetDn(const std::vector<int>& jetIndices,
//                                   const VLV& jets, const std::vector<int>& jetFlavour,
//                                   const std::vector<float>& btagDiscriminators)const
// {
//     if(!isGen_) return 1.;
//     if(!btagScaleFactors_LjetDn_) return 1.;
//     return btagScaleFactors_LjetDn_->getSF(jetIndices, jets, jetFlavour, btagDiscriminators);
// }


float MiniTree::weightKinReco_Up()const{
    if(!isGen_) return 1.;
    if(!kinematicReconstructionScaleFactors_Up_) return 1.;
    return kinematicReconstructionScaleFactors_Up_->getSF();
}
float MiniTree::weightKinReco_Dn()const{
    if(!isGen_) return 1.;
    if(!kinematicReconstructionScaleFactors_Dn_) return 1.;
    return kinematicReconstructionScaleFactors_Dn_->getSF();
}


float MiniTree::weightLooseKinReco_Up()const{
    if(!isGen_) return 1.;
    if(!looseKinRecoScaleFactors_Up_) return 1.;
    return looseKinRecoScaleFactors_Up_->getSF();
}
float MiniTree::weightLooseKinReco_Dn()const{
    if(!isGen_) return 1.;
    if(!looseKinRecoScaleFactors_Dn_) return 1.;
    return looseKinRecoScaleFactors_Dn_->getSF();
}
