//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Feb  4 14:56:13 2017 by ROOT version 6.08/02
// from TTree ttBar_treeVariables_step0/ttBar_treeVariables
// found on file: emu_ttbarsignalplustau.root
//////////////////////////////////////////////////////////

#ifndef makeTUnfoldHisto_h
#define makeTUnfoldHisto_h
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <TUnfoldBinning.h>
#include <string>
// Header file for the classes stored in the TTree if any.

class makeTUnfoldHisto {
protected:

public :
   TFile           *fout;
   //TTree         *fChain;   //!pointer to the analyzed TTree or TChain
   //TTree         *fChain0;   //!pointer to the analyzed TTree or TChain
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   TChain          *fChain0;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxb_Pnn = 1;
   const Int_t kMaxb_Mnn = 1;
   const Int_t kMaxb_Prr = 1;
   const Int_t kMaxb_Mrr = 1;
   const Int_t kMaxb_Pkk = 1;
   const Int_t kMaxb_Mkk = 1;
   const Int_t kMaxb_Pjj = 1;
   const Int_t kMaxb_Mjj = 1;
   const Int_t kMaxb_Pqq = 1;
   const Int_t kMaxb_Mqq = 1;

   const Int_t kMaxb1n = 1;
   const Int_t kMaxb2n = 1;
   const Int_t kMaxb1r = 1;
   const Int_t kMaxb2r = 1;
   const Int_t kMaxb1k = 1;
   const Int_t kMaxb2k = 1;
   const Int_t kMaxb1j = 1;
   const Int_t kMaxb2j = 1;
   const Int_t kMaxb1q = 1;
   const Int_t kMaxb2q = 1;

   const Int_t kMaxc_kk = 1;
   const Int_t kMaxc_rr = 1;
   const Int_t kMaxc_nn = 1;
   const Int_t kMaxc_kj = 1;
   const Int_t kMaxc_rq = 1;

   const Int_t kMaxc_Prk = 1;
   const Int_t kMaxc_Mrk = 1;
   const Int_t kMaxc_Pnr = 1;
   const Int_t kMaxc_Mnr = 1;
   const Int_t kMaxc_Pnk = 1;
   const Int_t kMaxc_Mnk = 1;
   const Int_t kMaxc_Prj = 1;
   const Int_t kMaxc_Mrj = 1;

   const Int_t kMaxc_rk = 1;
   const Int_t kMaxc_kr = 1;
   const Int_t kMaxc_nr = 1;
   const Int_t kMaxc_rn = 1;
   const Int_t kMaxc_nk = 1;
   const Int_t kMaxc_kn = 1;
   const Int_t kMaxc_rj = 1;
   const Int_t kMaxc_jr = 1;

   const Int_t kMaxll_cHel = 1;
   const Int_t kMaxll_cLab = 1;
   const Int_t kMaxll_sHel_k = 1;
   const Int_t kMaxll_sHel_r = 1;

   // Declaration of leaf types
   Float_t         top_pt;
   Float_t         top_phi;
   Float_t         top_rapidity;
   // Float_t      top_arapidity;
   Float_t         tbar_pt;
   Float_t         tbar_phi;
   Float_t         tbar_rapidity;
   // Float_t      tbar_arapidity;
   Float_t         l_pt;
   Float_t         lbar_pt;
   // Float_t      e_pt;
   // Float_t      ebar_pt;
   // Float_t      mu_pt;
   // Float_t      mubar_pt;
   Float_t         b_pt;
   Float_t         bbar_pt;
   Float_t         nu_pt;
   Float_t         nubar_pt;
   Float_t         l_eta;
   Float_t         lbar_eta;
   // Float_t      e_eta;
   // Float_t      ebar_eta;
   // Float_t      mu_eta;
   // Float_t      mubar_eta;
   Float_t         b_eta;
   Float_t         bbar_eta;
   Float_t         nu_eta;
   Float_t         nubar_eta;
   Float_t         l_phi;
   Float_t         lbar_phi;
   // Float_t      e_phi;
   // Float_t      ebar_phi;
   // Float_t      mu_phi;
   // Float_t      mubar_phi;
   Float_t         b_phi;
   Float_t         bbar_phi;
   Float_t         nu_phi;
   Float_t         nubar_phi;
   Float_t         met_pt;
   Float_t         ttbar_pt;
   Float_t         ttbar_phi;
   Float_t         ttbar_rapidity;
   // Float_t      ttbar_arapidity;
   Float_t         ttbar_delta_phi;
   Float_t         ttbar_delta_eta;
   Float_t         ttbar_delta_rapidity;
   Float_t         ttbar_mass;
   Float_t         llbar_pt;
   Float_t         llbar_phi;
   Float_t         llbar_rapidity;
   // Float_t      llbar_arapidity;
   Float_t         llbar_delta_phi;
   Float_t         llbar_delta_eta;
   Float_t         llbar_delta_rapidity;
   Float_t         llbar_mass;

   Float_t         all_mass;
   Float_t         r_mass;
   Float_t         jet_multiplicity;
   Float_t         bjet_multiplicity;
   Float_t         x1;
   Float_t         x2;

   Float_t         b1n;
   Float_t         b2n;
   Float_t         b1r;
   Float_t         b2r;
   Float_t         b1k;
   Float_t         b2k;
   Float_t         b1j;
   Float_t         b2j;
   Float_t         b1q;
   Float_t         b2q;

   Float_t         c_kk;
   Float_t         c_rr;
   Float_t         c_nn;
   Float_t         c_kj;
   Float_t         c_rq;

   Float_t         c_rk;
   Float_t         c_kr;
   Float_t         c_nr;
   Float_t         c_rn;
   Float_t         c_nk;
   Float_t         c_kn;
   Float_t         c_rj;
   Float_t         c_jr;

   Float_t         ll_cHel;
   Float_t         ll_cLab;
   Float_t         ll_kNorm;
   Float_t         ll_rNorm;

   Float_t         top_scatteringangle_ttbarframe;
   
   // *************
   // Gen variables 
   // *************

   Float_t         gen_top_pt;
   Float_t         gen_top_phi;
   Float_t         gen_top_rapidity;
   // Float_t      gen_top_arapidity;
   Float_t         gen_tbar_pt;
   Float_t         gen_tbar_phi;
   Float_t         gen_tbar_rapidity;
   // Float_t      gen_tbar_arapidity;
   Float_t         gen_l_pt;
   Float_t         gen_lbar_pt;
   // Float_t      gen_e_pt;
   // Float_t      gen_ebar_pt;
   // Float_t      gen_mu_pt;
   // Float_t      gen_mubar_pt;
   Float_t         gen_b_pt;
   Float_t         gen_bbar_pt;
   Float_t         gen_nu_pt;
   Float_t         gen_nubar_pt;
   Float_t         gen_l_eta;
   Float_t         gen_lbar_eta;
   // Float_t      gen_e_eta;
   // Float_t      gen_ebar_eta;
   // Float_t      gen_mu_eta;
   // Float_t      gen_mubar_eta;
   Float_t         gen_b_eta;
   Float_t         gen_bbar_eta;
   Float_t         gen_nu_eta;
   Float_t         gen_nubar_eta;
   Float_t         gen_l_phi;
   Float_t         gen_lbar_phi;
   Float_t         gen_e_phi;
   Float_t         gen_ebar_phi;
   Float_t         gen_mu_phi;
   Float_t         gen_mubar_phi;
   Float_t         gen_b_phi;
   Float_t         gen_bbar_phi;
   Float_t         gen_nu_phi;
   Float_t         gen_nubar_phi;
   Float_t         gen_ttbar_pt;
   Float_t         gen_ttbar_phi;
   Float_t         gen_ttbar_rapidity;
   // Float_t      gen_ttbar_arapidity;
   Float_t         gen_ttbar_delta_phi;
   Float_t         gen_ttbar_delta_eta;
   Float_t         gen_ttbar_delta_rapidity;
   Float_t         gen_ttbar_mass;
   Float_t         gen_llbar_pt;
   Float_t         gen_llbar_phi;
   Float_t         gen_llbar_rapidity;
   // Float_t      gen_llbar_arapidity;
   Float_t         gen_llbar_delta_phi;
   Float_t         gen_llbar_delta_eta;
   Float_t         gen_llbar_delta_rapidity;
   Float_t         gen_llbar_mass;

   Float_t         gen_all_mass;
   Float_t         gen_r_mass;
   Float_t         gen_jet_multiplicity;
   Float_t         gen_x1;
   Float_t         gen_x2;

   Float_t         gen_b1n;
   Float_t         gen_b2n;
   Float_t         gen_b1r;
   Float_t         gen_b2r;
   Float_t         gen_b1k;
   Float_t         gen_b2k;
   Float_t         gen_b1j;
   Float_t         gen_b2j;
   Float_t         gen_b1q;
   Float_t         gen_b2q;

   Float_t         gen_c_kk;
   Float_t         gen_c_rr;
   Float_t         gen_c_nn;
   Float_t         gen_c_kj;
   Float_t         gen_c_rq;

   Float_t         gen_c_rk;
   Float_t         gen_c_kr;
   Float_t         gen_c_nr;
   Float_t         gen_c_rn;
   Float_t         gen_c_nk;
   Float_t         gen_c_kn;
   Float_t         gen_c_rj;
   Float_t         gen_c_jr;

   Float_t         gen_ll_cHel;
   Float_t         gen_ll_cLab;
   Float_t         gen_ll_kNorm;
   Float_t         gen_ll_rNorm;

   Float_t         gen_top_scatteringangle_ttbarframe;

   Int_t           entry;
   Int_t           isTopGen;
   Int_t           isKinReco;
   Float_t         trueLevelWeight;
   Float_t         eventWeight;

   // ************
   // Step 0 stuff
   // ************

   Float_t         top_pt_0;
   Float_t         top_phi_0;
   Float_t         top_rapidity_0;
   // Float_t      top_arapidity_0;
   Float_t         tbar_pt_0;
   Float_t         tbar_phi_0;
   Float_t         tbar_rapidity_0;
   // Float_t      tbar_arapidity_0;
   Float_t         l_pt_0;
   Float_t         lbar_pt_0;
   // Float_t      e_pt_0;
   // Float_t      ebar_pt_0;
   // Float_t      mu_pt_0;
   // Float_t      mubar_pt_0;
   Float_t         b_pt_0;
   Float_t         bbar_pt_0;
   Float_t         nu_pt_0;
   Float_t         nubar_pt_0;
   Float_t         l_eta_0;
   Float_t         lbar_eta_0;
   // Float_t      e_eta_0;
   // Float_t      ebar_eta_0;
   // Float_t      mu_eta_0;
   // Float_t      mubar_eta_0;
   Float_t         b_eta_0;
   Float_t         bbar_eta_0;
   Float_t         nu_eta_0;
   Float_t         nubar_eta_0;
   Float_t         l_phi_0;
   Float_t         lbar_phi_0;
   // Float_t      e_phi_0;
   // Float_t      ebar_phi_0;
   // Float_t      mu_phi_0;
   // Float_t      mubar_phi_0;
   Float_t         b_phi_0;
   Float_t         bbar_phi_0;
   Float_t         nu_phi_0;
   Float_t         nubar_phi_0;
   Float_t         met_pt_0;
   Float_t         ttbar_pt_0;
   Float_t         ttbar_phi_0;
   Float_t         ttbar_rapidity_0;
   // Float_t      ttbar_arapidity_0;
   Float_t         ttbar_delta_phi_0;
   Float_t         ttbar_delta_eta_0;
   Float_t         ttbar_delta_rapidity_0;
   Float_t         ttbar_mass_0;
   Float_t         llbar_pt_0;
   Float_t         llbar_phi_0;
   Float_t         llbar_rapidity_0;
   // Float_t      llbar_arapidity_0;
   Float_t         llbar_delta_phi_0;
   Float_t         llbar_delta_eta_0;
   Float_t         llbar_delta_rapidity_0;
   Float_t         llbar_mass_0;

   Float_t         all_mass_0;
   Float_t         r_mass_0;
   Float_t         jet_multiplicity_0;
   Float_t         bjet_multiplicity_0;
   Float_t         x1_0;
   Float_t         x2_0;

   Float_t         b1n_0;
   Float_t         b2n_0;
   Float_t         b1r_0;
   Float_t         b2r_0;
   Float_t         b1k_0;
   Float_t         b2k_0;
   Float_t         b1j_0;
   Float_t         b2j_0;
   Float_t         b1q_0;
   Float_t         b2q_0;

   Float_t         c_kk_0;
   Float_t         c_rr_0;
   Float_t         c_nn_0;
   Float_t         c_kj_0;
   Float_t         c_rq_0;

   Float_t         c_rk_0;
   Float_t         c_kr_0;
   Float_t         c_nr_0;
   Float_t         c_rn_0;
   Float_t         c_nk_0;
   Float_t         c_kn_0;
   Float_t         c_rj_0;
   Float_t         c_jr_0;

   Float_t         ll_cHel_0;
   Float_t         ll_cLab_0;
   Float_t         ll_kNorm_0;
   Float_t         ll_rNorm_0;

   Float_t         top_scatteringangle_ttbarframe_0;

   // **********
   // Gen step 0
   // **********

   Float_t         gen_top_pt_0;
   Float_t         gen_top_phi_0;
   Float_t         gen_top_rapidity_0;
   // Float_t      gen_top_arapidity_0;
   Float_t         gen_tbar_pt_0;
   Float_t         gen_tbar_phi_0;
   Float_t         gen_tbar_rapidity_0;
   // Float_t      gen_tbar_arapidity_0;
   Float_t         gen_l_pt_0;
   Float_t         gen_lbar_pt_0;
   Float_t         gen_l_pdgid_0;
   Float_t         gen_lbar_pdgid_0;
   // Float_t      gen_e_pt_0;
   // Float_t      gen_ebar_pt_0;
   // Float_t      gen_mu_pt_0;
   // Float_t      gen_mubar_pt_0;
   Float_t         gen_b_pt_0;
   Float_t         gen_bbar_pt_0;
   Float_t         gen_nu_pt_0;
   Float_t         gen_nubar_pt_0;
   Float_t         gen_l_eta_0;
   Float_t         gen_lbar_eta_0;
   // Float_t      gen_e_eta_0;
   // Float_t      gen_ebar_eta_0;
   // Float_t      gen_mu_eta_0;
   // Float_t      gen_mubar_eta_0;
   Float_t         gen_b_eta_0;
   Float_t         gen_bbar_eta_0;
   Float_t         gen_nu_eta_0;
   Float_t         gen_nubar_eta_0;
   Float_t         gen_l_phi_0;
   Float_t         gen_lbar_phi_0;
   Float_t         gen_e_phi_0;
   Float_t         gen_ebar_phi_0;
   Float_t         gen_mu_phi_0;
   Float_t         gen_mubar_phi_0;
   Float_t         gen_b_phi_0;
   Float_t         gen_bbar_phi_0;
   Float_t         gen_nu_phi_0;
   Float_t         gen_nubar_phi_0;
   Float_t         gen_ttbar_pt_0;
   Float_t         gen_ttbar_phi_0;
   Float_t         gen_ttbar_rapidity_0;
   // Float_t      gen_ttbar_arapidity_0;
   Float_t         gen_ttbar_delta_phi_0;
   Float_t         gen_ttbar_delta_eta_0;
   Float_t         gen_ttbar_delta_rapidity_0;
   Float_t         gen_ttbar_mass_0;
   Float_t         gen_llbar_pt_0;
   Float_t         gen_llbar_phi_0;
   Float_t         gen_llbar_rapidity_0;
   // Float_t      gen_llbar_arapidity_0;
   Float_t         gen_llbar_delta_phi_0;
   Float_t         gen_llbar_delta_eta_0;
   Float_t         gen_llbar_delta_rapidity_0;
   Float_t         gen_llbar_mass_0;

   Float_t         gen_all_mass_0;
   Float_t         gen_r_mass_0;
   Float_t         gen_jet_multiplicity_0;
   Float_t         gen_x1_0;
   Float_t         gen_x2_0;

   Float_t         gen_b1n_0;
   Float_t         gen_b2n_0;
   Float_t         gen_b1r_0;
   Float_t         gen_b2r_0;
   Float_t         gen_b1k_0;
   Float_t         gen_b2k_0;
   Float_t         gen_b1j_0;
   Float_t         gen_b2j_0;
   Float_t         gen_b1q_0;
   Float_t         gen_b2q_0;

   Float_t         gen_c_kk_0;
   Float_t         gen_c_rr_0;
   Float_t         gen_c_nn_0;
   Float_t         gen_c_kj_0;
   Float_t         gen_c_rq_0;

   Float_t         gen_c_rk_0;
   Float_t         gen_c_kr_0;
   Float_t         gen_c_nr_0;
   Float_t         gen_c_rn_0;
   Float_t         gen_c_nk_0;
   Float_t         gen_c_kn_0;
   Float_t         gen_c_rj_0;
   Float_t         gen_c_jr_0;

   Float_t         gen_ll_cHel_0;
   Float_t         gen_ll_cLab_0;
   Float_t         gen_ll_kNorm_0;
   Float_t         gen_ll_rNorm_0;

   Float_t         gen_top_scatteringangle_ttbarframe_0;


   Int_t           entry_0;
   Int_t           isTopGen_0;
   Int_t           isKinReco_0;
   Float_t         trueLevelWeight_0;
   Float_t         eventWeight_0;

   // ****************
   // List of branches
   // ****************

   TBranch        *b_top_pt;   //!
   TBranch        *b_top_phi;   //!
   TBranch        *b_top_rapidity;   //!
   // TBranch     *b_top_arapidity;   //!
   TBranch        *b_tbar_pt;   //!
   TBranch        *b_tbar_phi;   //!
   TBranch        *b_tbar_rapidity;   //!
   // TBranch     *b_tbar_arapidity;   //!
   TBranch        *b_l_pt;   //!
   TBranch        *b_lbar_pt;   //!
   // TBranch     *b_e_pt;   //!
   // TBranch     *b_ebar_pt;   //!
   // TBranch     *b_mu_pt;   //!
   // TBranch     *b_mubar_pt;   //!
   TBranch        *b_b_pt;   //!
   TBranch        *b_bbar_pt;   //!
   TBranch        *b_nu_pt;   //!
   TBranch        *b_nubar_pt;   //!
   TBranch        *b_l_eta;   //!
   TBranch        *b_lbar_eta;   //!
   // TBranch     *b_e_eta;   //!
   // TBranch     *b_ebar_eta;   //!
   // TBranch     *b_mu_eta;   //!
   // TBranch     *b_mubar_eta;   //!
   TBranch        *b_b_eta;   //!
   TBranch        *b_bbar_eta;   //!
   TBranch        *b_nu_eta;   //!
   TBranch        *b_nubar_eta;   //!
   TBranch        *b_l_phi;   //!
   TBranch        *b_lbar_phi;   //!
   // TBranch     *b_e_phi;   //!
   // TBranch     *b_ebar_phi;   //!
   // TBranch     *b_mu_phi;   //!
   // TBranch     *b_mubar_phi;   //!
   TBranch        *b_b_phi;   //!
   TBranch        *b_bbar_phi;   //!
   TBranch        *b_nu_phi;   //!
   TBranch        *b_nubar_phi;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_ttbar_pt;   //!
   TBranch        *b_ttbar_phi;   //!
   TBranch        *b_ttbar_rapidity;   //!
   // TBranch     *b_ttbar_arapidity;   //!
   TBranch        *b_ttbar_delta_phi;   //!
   TBranch        *b_ttbar_delta_eta;   //!
   TBranch        *b_ttbar_delta_rapidity;   //!
   TBranch        *b_ttbar_mass;   //!
   TBranch        *b_llbar_pt;   //!
   TBranch        *b_llbar_phi;   //!
   TBranch        *b_llbar_rapidity;   //!
   // TBranch     *b_llbar_arapidity;   //!
   TBranch        *b_llbar_delta_phi;   //!
   TBranch        *b_llbar_delta_eta;   //!
   TBranch        *b_llbar_delta_rapidity;   //!
   TBranch        *b_llbar_mass;   //!

   TBranch        *b_all_mass;   //!
   TBranch        *b_r_mass;   //!
   TBranch        *b_jet_multiplicity;   //!
   TBranch        *b_bjet_multiplicity;   //!
   TBranch        *b_x1;   //!
   TBranch        *b_x2;   //!

   TBranch        *b_b1n;   //!
   TBranch        *b_b2n;   //!
   TBranch        *b_b1r;   //!
   TBranch        *b_b2r;   //!
   TBranch        *b_b1k;   //!
   TBranch        *b_b2k;   //!
   TBranch        *b_b1j;   //!
   TBranch        *b_b2j;   //!
   TBranch        *b_b1q;   //!
   TBranch        *b_b2q;   //!

   TBranch        *b_c_kk;   //!
   TBranch        *b_c_rr;   //!
   TBranch        *b_c_nn;   //!
   TBranch        *b_c_kj;   //!
   TBranch        *b_c_rq;   //!

   TBranch        *b_c_rk;   //!
   TBranch        *b_c_kr;   //!
   TBranch        *b_c_nr;   //!
   TBranch        *b_c_rn;   //!
   TBranch        *b_c_nk;   //!
   TBranch        *b_c_kn;   //!
   TBranch        *b_c_rj;   //!
   TBranch        *b_c_jr;   //!

   TBranch        *b_ll_cHel;   //!
   TBranch        *b_ll_cLab;   //!
   TBranch        *b_ll_kNorm;   //!
   TBranch        *b_ll_rNorm;   //!

   TBranch        *b_top_scatteringangle_ttbarframe;

   // ********************
   // List of gen branches
   // ********************

   TBranch        *b_gen_top_pt;   //!
   TBranch        *b_gen_top_phi;   //!
   TBranch        *b_gen_top_rapidity;   //!
   // TBranch     *b_gen_top_arapidity;   //!
   TBranch        *b_gen_tbar_pt;   //!
   TBranch        *b_gen_tbar_phi;   //!
   TBranch        *b_gen_tbar_rapidity;   //!
   // TBranch     *b_gen_tbar_arapidity;   //!
   TBranch        *b_gen_l_pt;   //!
   TBranch        *b_gen_lbar_pt;   //!
   TBranch        *b_gen_l_pdgid;   //!
   TBranch        *b_gen_lbar_pdgid;   //!
   TBranch        *b_gen_e_pt;   //!
   TBranch        *b_gen_ebar_pt;   //!
   TBranch        *b_gen_mu_pt;   //!
   TBranch        *b_gen_mubar_pt;   //!
   TBranch        *b_gen_b_pt;   //!
   TBranch        *b_gen_bbar_pt;   //!
   TBranch        *b_gen_nu_pt;   //!
   TBranch        *b_gen_nubar_pt;   //!
   TBranch        *b_gen_l_eta;   //!
   TBranch        *b_gen_lbar_eta;   //!
   // TBranch     *b_gen_e_eta;   //!
   // TBranch     *b_gen_ebar_eta;   //!
   // TBranch     *b_gen_mu_eta;   //!
   // TBranch     *b_gen_mubar_eta;   //!
   TBranch        *b_gen_b_eta;   //!
   TBranch        *b_gen_bbar_eta;   //!
   TBranch        *b_gen_nu_eta;   //!
   TBranch        *b_gen_nubar_eta;   //!
   TBranch        *b_gen_l_phi;   //!
   TBranch        *b_gen_lbar_phi;   //!
   // TBranch     *b_gen_e_phi;   //!
   // TBranch     *b_gen_ebar_phi;   //!
   // TBranch     *b_gen_mu_phi;   //!
   // TBranch     *b_gen_mubar_phi;   //!
   TBranch        *b_gen_b_phi;   //!
   TBranch        *b_gen_bbar_phi;   //!
   TBranch        *b_gen_nu_phi;   //!
   TBranch        *b_gen_nubar_phi;   //!
   TBranch        *b_gen_ttbar_pt;   //!
   TBranch        *b_gen_ttbar_phi;   //!
   TBranch        *b_gen_ttbar_rapidity;   //!
   // TBranch     *b_gen_ttbar_arapidity;   //!
   TBranch        *b_gen_ttbar_delta_phi;   //!
   TBranch        *b_gen_ttbar_delta_eta;   //!
   TBranch        *b_gen_ttbar_delta_rapidity;   //!
   TBranch        *b_gen_ttbar_mass;   //!
   TBranch        *b_gen_llbar_pt;   //!
   TBranch        *b_gen_llbar_phi;   //!
   TBranch        *b_gen_llbar_rapidity;   //!
   // TBranch     *b_gen_llbar_arapidity;   //!
   TBranch        *b_gen_llbar_delta_phi;   //!
   TBranch        *b_gen_llbar_delta_eta;   //!
   TBranch        *b_gen_llbar_delta_rapidity;   //!
   TBranch        *b_gen_llbar_mass;   //!

   TBranch        *b_gen_all_mass;   //!
   TBranch        *b_gen_r_mass;   //!
   TBranch        *b_gen_jet_multiplicity;   //!
   TBranch        *b_gen_x1;   //!
   TBranch        *b_gen_x2;   //!

   TBranch        *b_gen_b1n;   //!
   TBranch        *b_gen_b2n;   //!
   TBranch        *b_gen_b1r;   //!
   TBranch        *b_gen_b2r;   //!
   TBranch        *b_gen_b1k;   //!
   TBranch        *b_gen_b2k;   //!
   TBranch        *b_gen_b1j;   //!
   TBranch        *b_gen_b2j;   //!
   TBranch        *b_gen_b1q;   //!
   TBranch        *b_gen_b2q;   //!

   TBranch        *b_gen_c_kk;   //!
   TBranch        *b_gen_c_rr;   //!
   TBranch        *b_gen_c_nn;   //!
   TBranch        *b_gen_c_kj;   //!
   TBranch        *b_gen_c_rq;   //!

   TBranch        *b_gen_c_rk;   //!
   TBranch        *b_gen_c_kr;   //!
   TBranch        *b_gen_c_nr;   //!
   TBranch        *b_gen_c_rn;   //!
   TBranch        *b_gen_c_nk;   //!
   TBranch        *b_gen_c_kn;   //!
   TBranch        *b_gen_c_rj;   //!
   TBranch        *b_gen_c_jr;   //!

   TBranch        *b_gen_ll_cHel;   //!
   TBranch        *b_gen_ll_cLab;   //!
   TBranch        *b_gen_ll_kNorm;   //!
   TBranch        *b_gen_ll_rNorm;   //!
   TBranch        *b_gen_top_scatteringangle_ttbarframe;

   TBranch        *b_entry;   //!
   TBranch        *b_isTopGen;   //!
   TBranch        *b_isKinReco;   //!
   TBranch        *b_trueLevelWeight;   //!
   TBranch        *b_eventWeight;   //!


   //makeTUnfoldHisto(TTree *tree0=0, TTree *tree8=0 );
   makeTUnfoldHisto(TChain *Chain0=0, TChain *Chain8=0 );
   virtual ~makeTUnfoldHisto();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);   
   virtual Long64_t LoadTree0(Long64_t entry);
   //virtual void   Init(TTree *tree0, TTree *tree8);
   virtual void     Init(TChain *Chain0, TChain *Chain8);
   virtual void     Loop(TH1D *wtdEvts, std::string channel, bool symmetrizeSpinVar);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void SetOutputFileName(std::string output, std::string systematics, std::string channel,  std::string era);
   void fillUnderOverFlow(TH1 *h1, TUnfoldBinning *binning, std::vector<double> v_values, double weight, bool symmetrizeSpinVar = false);
   void fillUnderOverFlow(TH2 *h2, TUnfoldBinning *xBinning, TUnfoldBinning *yBinning, std::vector<double> v_xvalues, std::vector<double> v_yvalues, double weight, bool isReconstructed = true, bool symmetrizeSpinVar = false);
   void fillUnderOverFlowDeltaWeight(TH2 *h2, TUnfoldBinning *xBinning, std::vector<double> v_xvalues, double weightGen, double weightReco, bool symmetrizeSpinVar = false);
   void BookHisto(int n_sub_bins = 1);
   void Terminate(double scalingvalue = 1);
   void Reset();
   TUnfoldBinning* CreateResidualBinning(int nbinsres, double lowerBinEdge, double upperBinEdge);
   std::pair<TUnfoldBinning*, TUnfoldBinning*> CreateXMLBinFiles(std::string varname, std::vector<std::string> axisNames, std::vector<int> gen_nbin, std::vector<int> reco_nbin, std::vector<double*> gen_bins, std::vector<double*>reco_bins, int numBinSubdivisions = 1);
   void CreateFineBinningArray(const int& input_nbin, const int& output_nbin, double input_bins[], double output_bins[]);
};

#endif

#ifdef makeTUnfoldHisto_cxx
/*makeTUnfoldHisto::makeTUnfoldHisto(TTree *tree0, TTree *tree8) : fLuminosity(0), fXsection(0), fout(0), fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree0 == 0 || tree8 == 0) {
      //should give error **	
      std::cout<<"Couldn't find any tree!!!"<<std::endl;

   }
   Init(tree0, tree8);
}*/
makeTUnfoldHisto::makeTUnfoldHisto(TChain *chain0, TChain *chain8) : fout(0), fChain(0)
{
// if parameter Chain is not specified (or zero), connect the file
// used to generate this class and read the Chain.
   if (chain0 == 0 || chain8 == 0) {
      //should give error **	
      std::cout<<"Couldn't find any Chain!!!"<<std::endl;

   }
   Init(chain0, chain8);
}

makeTUnfoldHisto::~makeTUnfoldHisto()
{
//   if (!fChain) return;
//   delete fChain->GetCurrentFile();
}

Int_t makeTUnfoldHisto::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t makeTUnfoldHisto::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}
Long64_t makeTUnfoldHisto::LoadTree0(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain0) return -5;
   Long64_t centry = fChain0->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain0->GetTreeNumber() != fCurrent) {
      fCurrent = fChain0->GetTreeNumber();
      Notify();
   }
   return centry;
}

//void makeTUnfoldHisto::Init(TTree *tree0, TTree *tree8)
void makeTUnfoldHisto::Init(TChain *Chain0, TChain *Chain8)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers to the reco tree
   //if (!tree8) return;
   if (!Chain8) return;
   //fChain = tree8;
   fChain   = Chain8;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("top_pt", &top_pt, &b_top_pt);
   fChain->SetBranchAddress("top_phi", &top_phi, &b_top_phi);
   fChain->SetBranchAddress("top_rapidity", &top_rapidity, &b_top_rapidity);
   // fChain->SetBranchAddress("top_arapidity", &top_arapidity, &b_top_arapidity);
   fChain->SetBranchAddress("tbar_pt", &tbar_pt, &b_tbar_pt);
   fChain->SetBranchAddress("tbar_phi", &tbar_phi, &b_tbar_phi);
   fChain->SetBranchAddress("tbar_rapidity", &tbar_rapidity, &b_tbar_rapidity);
   // fChain->SetBranchAddress("tbar_arapidity", &tbar_arapidity, &b_tbar_arapidity);
   fChain->SetBranchAddress("l_pt", &l_pt, &b_l_pt);
   fChain->SetBranchAddress("lbar_pt", &lbar_pt, &b_lbar_pt);
   // fChain->SetBranchAddress("e_pt", &e_pt, &b_e_pt);
   // fChain->SetBranchAddress("ebar_pt", &ebar_pt, &b_ebar_pt);
   // fChain->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
   // fChain->SetBranchAddress("mubar_pt", &mubar_pt, &b_mubar_pt);
   fChain->SetBranchAddress("b_pt", &b_pt, &b_b_pt);
   fChain->SetBranchAddress("bbar_pt", &bbar_pt, &b_bbar_pt);
   fChain->SetBranchAddress("nu_pt", &nu_pt, &b_nu_pt);
   fChain->SetBranchAddress("nubar_pt", &nubar_pt, &b_nubar_pt);
   fChain->SetBranchAddress("l_eta", &l_eta, &b_l_eta);
   fChain->SetBranchAddress("lbar_eta", &lbar_eta, &b_lbar_eta);
   // fChain->SetBranchAddress("e_eta", &e_eta, &b_e_eta);
   // fChain->SetBranchAddress("ebar_eta", &ebar_eta, &b_ebar_eta);
   // fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
   // fChain->SetBranchAddress("mubar_eta", &mubar_eta, &b_mubar_eta);
   fChain->SetBranchAddress("b_eta", &b_eta, &b_b_eta);
   fChain->SetBranchAddress("bbar_eta", &bbar_eta, &b_bbar_eta);
   fChain->SetBranchAddress("nu_eta", &nu_eta, &b_nu_eta);
   fChain->SetBranchAddress("nubar_eta", &nubar_eta, &b_nubar_eta);
   fChain->SetBranchAddress("l_phi", &l_phi, &b_l_phi);
   fChain->SetBranchAddress("lbar_phi", &lbar_phi, &b_lbar_phi);
   // fChain->SetBranchAddress("e_phi", &e_phi, &b_e_phi);
   // fChain->SetBranchAddress("ebar_phi", &ebar_phi, &b_ebar_phi);
   // fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
   // fChain->SetBranchAddress("mubar_phi", &mubar_phi, &b_mubar_phi);
   fChain->SetBranchAddress("b_phi", &b_phi, &b_b_phi);
   fChain->SetBranchAddress("bbar_phi", &bbar_phi, &b_bbar_phi);
   fChain->SetBranchAddress("nu_phi", &nu_phi, &b_nu_phi);
   fChain->SetBranchAddress("nubar_phi", &nubar_phi, &b_nubar_phi);
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("ttbar_pt", &ttbar_pt, &b_ttbar_pt);
   fChain->SetBranchAddress("ttbar_phi", &ttbar_phi, &b_ttbar_phi);
   fChain->SetBranchAddress("ttbar_rapidity", &ttbar_rapidity, &b_ttbar_rapidity);
   // fChain->SetBranchAddress("ttbar_arapidity", &ttbar_arapidity, &b_ttbar_arapidity);
   fChain->SetBranchAddress("ttbar_delta_phi", &ttbar_delta_phi, &b_ttbar_delta_phi);
   fChain->SetBranchAddress("ttbar_delta_eta", &ttbar_delta_eta, &b_ttbar_delta_eta);
   fChain->SetBranchAddress("ttbar_delta_rapidity", &ttbar_delta_rapidity, &b_ttbar_delta_rapidity);
   fChain->SetBranchAddress("ttbar_mass", &ttbar_mass, &b_ttbar_mass);
   fChain->SetBranchAddress("llbar_pt", &llbar_pt, &b_llbar_pt);
   fChain->SetBranchAddress("llbar_phi", &llbar_phi, &b_llbar_phi);
   fChain->SetBranchAddress("llbar_rapidity", &llbar_rapidity, &b_llbar_rapidity);
   // fChain->SetBranchAddress("llbar_arapidity", &llbar_arapidity, &b_llbar_arapidity);
   fChain->SetBranchAddress("llbar_delta_phi", &llbar_delta_phi, &b_llbar_delta_phi);
   fChain->SetBranchAddress("llbar_delta_eta", &llbar_delta_eta, &b_llbar_delta_eta);
   fChain->SetBranchAddress("llbar_delta_rapidity", &llbar_delta_rapidity, &b_llbar_delta_rapidity);
   fChain->SetBranchAddress("llbar_mass", &llbar_mass, &b_llbar_mass);

   fChain->SetBranchAddress("all_mass", &all_mass, &b_all_mass);
   fChain->SetBranchAddress("r_mass", &r_mass, &b_r_mass);
   fChain->SetBranchAddress("jet_multiplicity", &jet_multiplicity, &b_jet_multiplicity);
   fChain->SetBranchAddress("bjet_multiplicity", &bjet_multiplicity, &b_bjet_multiplicity);
   fChain->SetBranchAddress("x1", &x1, &b_x1);
   fChain->SetBranchAddress("x2", &x2, &b_x2);

   fChain->SetBranchAddress("b1n", &b1n, &b_b1n);
   fChain->SetBranchAddress("b2n", &b2n, &b_b2n);
   fChain->SetBranchAddress("b1r", &b1r, &b_b1r);
   fChain->SetBranchAddress("b2r", &b2r, &b_b2r);
   fChain->SetBranchAddress("b1k", &b1k, &b_b1k);
   fChain->SetBranchAddress("b2k", &b2k, &b_b2k);
   fChain->SetBranchAddress("b1j", &b1j, &b_b1j);
   fChain->SetBranchAddress("b2j", &b2j, &b_b2j);
   fChain->SetBranchAddress("b1q", &b1q, &b_b1q);
   fChain->SetBranchAddress("b2q", &b2q, &b_b2q);

   fChain->SetBranchAddress("c_kk", &c_kk, &b_c_kk);
   fChain->SetBranchAddress("c_rr", &c_rr, &b_c_rr);
   fChain->SetBranchAddress("c_nn", &c_nn, &b_c_nn);
   fChain->SetBranchAddress("c_kj", &c_kj, &b_c_kj);
   fChain->SetBranchAddress("c_rq", &c_rq, &b_c_rq);

   fChain->SetBranchAddress("c_rk", &c_rk, &b_c_rk);
   fChain->SetBranchAddress("c_kr", &c_kr, &b_c_kr);
   fChain->SetBranchAddress("c_nr", &c_nr, &b_c_nr);
   fChain->SetBranchAddress("c_rn", &c_rn, &b_c_rn);
   fChain->SetBranchAddress("c_nk", &c_nk, &b_c_nk);
   fChain->SetBranchAddress("c_kn", &c_kn, &b_c_kn);
   fChain->SetBranchAddress("c_rj", &c_rj, &b_c_rj);
   fChain->SetBranchAddress("c_jr", &c_jr, &b_c_jr);

   fChain->SetBranchAddress("ll_cHel", &ll_cHel, &b_ll_cHel);
   fChain->SetBranchAddress("ll_cLab", &ll_cLab, &b_ll_cLab);
   fChain->SetBranchAddress("ll_kNorm", &ll_kNorm, &b_ll_kNorm);
   fChain->SetBranchAddress("ll_rNorm", &ll_rNorm, &b_ll_rNorm);

   fChain->SetBranchAddress("top_scatteringangle_ttbarframe", &top_scatteringangle_ttbarframe, &b_top_scatteringangle_ttbarframe);

   // *************
   // Gen variables
   // *************

   fChain->SetBranchAddress("gen_top_pt", &gen_top_pt, &b_gen_top_pt);
   fChain->SetBranchAddress("gen_top_phi", &gen_top_phi, &b_gen_top_phi);
   fChain->SetBranchAddress("gen_top_rapidity", &gen_top_rapidity, &b_gen_top_rapidity);
   // fChain->SetBranchAddress("gen_top_arapidity", &gen_top_arapidity, &b_gen_top_arapidity);
   fChain->SetBranchAddress("gen_tbar_pt", &gen_tbar_pt, &b_gen_tbar_pt);
   fChain->SetBranchAddress("gen_tbar_phi", &gen_tbar_phi, &b_gen_tbar_phi);
   fChain->SetBranchAddress("gen_tbar_rapidity", &gen_tbar_rapidity, &b_gen_tbar_rapidity);
   // fChain->SetBranchAddress("gen_tbar_arapidity", &gen_tbar_arapidity, &b_gen_tbar_arapidity);
   fChain->SetBranchAddress("gen_l_pt", &gen_l_pt, &b_gen_l_pt);
   fChain->SetBranchAddress("gen_lbar_pt", &gen_lbar_pt, &b_gen_lbar_pt);
   // fChain->SetBranchAddress("gen_e_pt", &gen_e_pt, &b_gen_e_pt);
   // fChain->SetBranchAddress("gen_ebar_pt", &gen_ebar_pt, &b_gen_ebar_pt);
   // fChain->SetBranchAddress("gen_mu_pt", &gen_mu_pt, &b_gen_mu_pt);
   // fChain->SetBranchAddress("gen_mubar_pt", &gen_mubar_pt, &b_gen_mubar_pt);
   fChain->SetBranchAddress("gen_b_pt", &gen_b_pt, &b_gen_b_pt);
   fChain->SetBranchAddress("gen_bbar_pt", &gen_bbar_pt, &b_gen_bbar_pt);
   fChain->SetBranchAddress("gen_nu_pt", &gen_nu_pt, &b_gen_nu_pt);
   fChain->SetBranchAddress("gen_nubar_pt", &gen_nubar_pt, &b_gen_nubar_pt);
   fChain->SetBranchAddress("gen_l_eta", &gen_l_eta, &b_gen_l_eta);
   fChain->SetBranchAddress("gen_lbar_eta", &gen_lbar_eta, &b_gen_lbar_eta);
   // fChain->SetBranchAddress("gen_e_eta", &gen_e_eta, &b_gen_e_eta);
   // fChain->SetBranchAddress("gen_ebar_eta", &gen_ebar_eta, &b_gen_ebar_eta);
   // fChain->SetBranchAddress("gen_mu_eta", &gen_mu_eta, &b_gen_mu_eta);
   // fChain->SetBranchAddress("gen_mubar_eta", &gen_mubar_eta, &b_gen_mubar_eta);
   fChain->SetBranchAddress("gen_b_eta", &gen_b_eta, &b_gen_b_eta);
   fChain->SetBranchAddress("gen_bbar_eta", &gen_bbar_eta, &b_gen_bbar_eta);
   fChain->SetBranchAddress("gen_nu_eta", &gen_nu_eta, &b_gen_nu_eta);
   fChain->SetBranchAddress("gen_nubar_eta", &gen_nubar_eta, &b_gen_nubar_eta);
   fChain->SetBranchAddress("gen_l_phi", &gen_l_phi, &b_gen_l_phi);
   fChain->SetBranchAddress("gen_lbar_phi", &gen_lbar_phi, &b_gen_lbar_phi);
   // fChain->SetBranchAddress("gen_e_phi", &gen_e_phi, &b_gen_e_phi);
   // fChain->SetBranchAddress("gen_ebar_phi", &gen_ebar_phi, &b_gen_ebar_phi);
   // fChain->SetBranchAddress("gen_mu_phi", &gen_mu_phi, &b_gen_mu_phi);
   // fChain->SetBranchAddress("gen_mubar_phi", &gen_mubar_phi, &b_gen_mubar_phi);
   fChain->SetBranchAddress("gen_b_phi", &gen_b_phi, &b_gen_b_phi);
   fChain->SetBranchAddress("gen_bbar_phi", &gen_bbar_phi, &b_gen_bbar_phi);
   fChain->SetBranchAddress("gen_nu_phi", &gen_nu_phi, &b_gen_nu_phi);
   fChain->SetBranchAddress("gen_nubar_phi", &gen_nubar_phi, &b_gen_nubar_phi);
   fChain->SetBranchAddress("gen_ttbar_pt", &gen_ttbar_pt, &b_gen_ttbar_pt);
   fChain->SetBranchAddress("gen_ttbar_phi", &gen_ttbar_phi, &b_gen_ttbar_phi);
   fChain->SetBranchAddress("gen_ttbar_rapidity", &gen_ttbar_rapidity, &b_gen_ttbar_rapidity);
   // fChain->SetBranchAddress("gen_ttbar_arapidity", &gen_ttbar_arapidity, &b_gen_ttbar_arapidity);
   fChain->SetBranchAddress("gen_ttbar_delta_phi", &gen_ttbar_delta_phi, &b_gen_ttbar_delta_phi);
   fChain->SetBranchAddress("gen_ttbar_delta_eta", &gen_ttbar_delta_eta, &b_gen_ttbar_delta_eta);
   fChain->SetBranchAddress("gen_ttbar_delta_rapidity", &gen_ttbar_delta_rapidity, &b_gen_ttbar_delta_rapidity);
   fChain->SetBranchAddress("gen_ttbar_mass", &gen_ttbar_mass, &b_gen_ttbar_mass);
   fChain->SetBranchAddress("gen_llbar_pt", &gen_llbar_pt, &b_gen_llbar_pt);
   fChain->SetBranchAddress("gen_llbar_phi", &gen_llbar_phi, &b_gen_llbar_phi);
   fChain->SetBranchAddress("gen_llbar_rapidity", &gen_llbar_rapidity, &b_gen_llbar_rapidity);
   // fChain->SetBranchAddress("gen_llbar_arapidity", &gen_llbar_arapidity, &b_gen_llbar_arapidity);
   fChain->SetBranchAddress("gen_llbar_delta_phi", &gen_llbar_delta_phi, &b_gen_llbar_delta_phi);
   fChain->SetBranchAddress("gen_llbar_delta_eta", &gen_llbar_delta_eta, &b_gen_llbar_delta_eta);
   fChain->SetBranchAddress("gen_llbar_delta_rapidity", &gen_llbar_delta_rapidity, &b_gen_llbar_delta_rapidity);
   fChain->SetBranchAddress("gen_llbar_mass", &gen_llbar_mass, &b_gen_llbar_mass);

   fChain->SetBranchAddress("gen_all_mass", &gen_all_mass, &b_gen_all_mass);
   fChain->SetBranchAddress("gen_r_mass", &gen_r_mass, &b_gen_r_mass);
   fChain->SetBranchAddress("gen_jet_multiplicity", &gen_jet_multiplicity, &b_gen_jet_multiplicity);
   fChain->SetBranchAddress("gen_x1", &gen_x1, &b_gen_x1);
   fChain->SetBranchAddress("gen_x2", &gen_x2, &b_gen_x2);

   fChain->SetBranchAddress("gen_b1n", &gen_b1n, &b_gen_b1n);
   fChain->SetBranchAddress("gen_b2n", &gen_b2n, &b_gen_b2n);
   fChain->SetBranchAddress("gen_b1r", &gen_b1r, &b_gen_b1r);
   fChain->SetBranchAddress("gen_b2r", &gen_b2r, &b_gen_b2r);
   fChain->SetBranchAddress("gen_b1k", &gen_b1k, &b_gen_b1k);
   fChain->SetBranchAddress("gen_b2k", &gen_b2k, &b_gen_b2k);
   fChain->SetBranchAddress("gen_b1j", &gen_b1j, &b_gen_b1j);
   fChain->SetBranchAddress("gen_b2j", &gen_b2j, &b_gen_b2j);
   fChain->SetBranchAddress("gen_b1q", &gen_b1q, &b_gen_b1q);
   fChain->SetBranchAddress("gen_b2q", &gen_b2q, &b_gen_b2q);

   fChain->SetBranchAddress("gen_c_kk", &gen_c_kk, &b_gen_c_kk);
   fChain->SetBranchAddress("gen_c_rr", &gen_c_rr, &b_gen_c_rr);
   fChain->SetBranchAddress("gen_c_nn", &gen_c_nn, &b_gen_c_nn);
   fChain->SetBranchAddress("gen_c_kj", &gen_c_kj, &b_gen_c_kj);
   fChain->SetBranchAddress("gen_c_rq", &gen_c_rq, &b_gen_c_rq);

   fChain->SetBranchAddress("gen_c_rk", &gen_c_rk, &b_gen_c_rk);
   fChain->SetBranchAddress("gen_c_kr", &gen_c_kr, &b_gen_c_kr);
   fChain->SetBranchAddress("gen_c_nr", &gen_c_nr, &b_gen_c_nr);
   fChain->SetBranchAddress("gen_c_rn", &gen_c_rn, &b_gen_c_rn);
   fChain->SetBranchAddress("gen_c_nk", &gen_c_nk, &b_gen_c_nk);
   fChain->SetBranchAddress("gen_c_kn", &gen_c_kn, &b_gen_c_kn);
   fChain->SetBranchAddress("gen_c_rj", &gen_c_rj, &b_gen_c_rj);
   fChain->SetBranchAddress("gen_c_jr", &gen_c_jr, &b_gen_c_jr);

   fChain->SetBranchAddress("gen_ll_cHel", &gen_ll_cHel, &b_gen_ll_cHel);
   fChain->SetBranchAddress("gen_ll_cLab", &gen_ll_cLab, &b_gen_ll_cLab);
   fChain->SetBranchAddress("gen_ll_kNorm", &gen_ll_kNorm, &b_gen_ll_kNorm);
   fChain->SetBranchAddress("gen_ll_rNorm", &gen_ll_rNorm, &b_gen_ll_rNorm);

   fChain->SetBranchAddress("gen_top_scatteringangle_ttbarframe", &gen_top_scatteringangle_ttbarframe, &b_gen_top_scatteringangle_ttbarframe);

   fChain->SetBranchAddress("entry", &entry, &b_entry);
   fChain->SetBranchAddress("isTopGen", &isTopGen, &b_isTopGen);
   fChain->SetBranchAddress("isKinReco", &isKinReco, &b_isKinReco);
   fChain->SetBranchAddress("trueLevelWeight", &trueLevelWeight, &b_trueLevelWeight);
   fChain->SetBranchAddress("eventWeight", &eventWeight, &b_eventWeight);

   // Set branch addresses and branch pointers to the gen tree
   //if (!tree0) return;
   //fChain0 = tree0;
   if (!Chain0) return;
   fChain0 = Chain0;
   fCurrent = -1;
   fChain0->SetMakeClass(1);

   fChain0->SetBranchAddress("top_pt", &top_pt_0, &b_top_pt);
   fChain0->SetBranchAddress("top_phi", &top_phi_0, &b_top_phi);
   fChain0->SetBranchAddress("top_rapidity", &top_rapidity_0, &b_top_rapidity);
   // fChain0->SetBranchAddress("top_arapidity", &top_arapidity_0, &b_top_arapidity);
   fChain0->SetBranchAddress("tbar_pt", &tbar_pt_0, &b_tbar_pt);
   fChain0->SetBranchAddress("tbar_phi", &tbar_phi_0, &b_tbar_phi);
   fChain0->SetBranchAddress("tbar_rapidity", &tbar_rapidity_0, &b_tbar_rapidity);
   // fChain0->SetBranchAddress("tbar_arapidity", &tbar_arapidity_0, &b_tbar_arapidity);
   fChain0->SetBranchAddress("l_pt", &l_pt_0, &b_l_pt);
   fChain0->SetBranchAddress("lbar_pt", &lbar_pt_0, &b_lbar_pt);
   // fChain0->SetBranchAddress("e_pt", &e_pt_0, &b_e_pt);
   // fChain0->SetBranchAddress("ebar_pt", &ebar_pt_0, &b_ebar_pt);
   // fChain0->SetBranchAddress("mu_pt", &mu_pt_0, &b_mu_pt);
   // fChain0->SetBranchAddress("mubar_pt", &mubar_pt_0, &b_mubar_pt);
   fChain0->SetBranchAddress("b_pt", &b_pt_0, &b_b_pt);
   fChain0->SetBranchAddress("bbar_pt", &bbar_pt_0, &b_bbar_pt);
   fChain0->SetBranchAddress("nu_pt", &nu_pt_0, &b_nu_pt);
   fChain0->SetBranchAddress("nubar_pt", &nubar_pt_0, &b_nubar_pt);
   fChain0->SetBranchAddress("l_eta", &l_eta_0, &b_l_eta);
   fChain0->SetBranchAddress("lbar_eta", &lbar_eta_0, &b_lbar_eta);
   // fChain0->SetBranchAddress("e_eta", &e_eta_0, &b_e_eta);
   // fChain0->SetBranchAddress("ebar_eta", &ebar_eta_0, &b_ebar_eta);
   // fChain0->SetBranchAddress("mu_eta", &mu_eta_0, &b_mu_eta);
   // fChain0->SetBranchAddress("mubar_eta", &mubar_eta_0, &b_mubar_eta);
   fChain0->SetBranchAddress("b_eta", &b_eta_0, &b_b_eta);
   fChain0->SetBranchAddress("bbar_eta", &bbar_eta_0, &b_bbar_eta);
   fChain0->SetBranchAddress("nu_eta", &nu_eta_0, &b_nu_eta);
   fChain0->SetBranchAddress("nubar_eta", &nubar_eta_0, &b_nubar_eta);
   fChain0->SetBranchAddress("l_phi", &l_phi_0, &b_l_phi);
   fChain0->SetBranchAddress("lbar_phi", &lbar_phi_0, &b_lbar_phi);
   // fChain0->SetBranchAddress("e_phi", &e_phi_0, &b_e_phi);
   // fChain0->SetBranchAddress("ebar_phi", &ebar_phi_0, &b_ebar_phi);
   // fChain0->SetBranchAddress("mu_phi", &mu_phi_0, &b_mu_phi);
   // fChain0->SetBranchAddress("mubar_phi", &mubar_phi_0, &b_mubar_phi);
   fChain0->SetBranchAddress("b_phi", &b_phi_0, &b_b_phi);
   fChain0->SetBranchAddress("bbar_phi", &bbar_phi_0, &b_bbar_phi);
   fChain0->SetBranchAddress("nu_phi", &nu_phi_0, &b_nu_phi);
   fChain0->SetBranchAddress("nubar_phi", &nubar_phi_0, &b_nubar_phi);
   fChain0->SetBranchAddress("met_pt", &met_pt_0, &b_met_pt);
   fChain0->SetBranchAddress("ttbar_pt", &ttbar_pt_0, &b_ttbar_pt);
   fChain0->SetBranchAddress("ttbar_phi", &ttbar_phi_0, &b_ttbar_phi);
   fChain0->SetBranchAddress("ttbar_rapidity", &ttbar_rapidity_0, &b_ttbar_rapidity);
   // fChain0->SetBranchAddress("ttbar_arapidity", &ttbar_arapidity_0, &b_ttbar_arapidity);
   fChain0->SetBranchAddress("ttbar_delta_phi", &ttbar_delta_phi_0, &b_ttbar_delta_phi);
   fChain0->SetBranchAddress("ttbar_delta_eta", &ttbar_delta_eta_0, &b_ttbar_delta_eta);
   fChain0->SetBranchAddress("ttbar_delta_rapidity", &ttbar_delta_rapidity_0, &b_ttbar_delta_rapidity);
   fChain0->SetBranchAddress("ttbar_mass", &ttbar_mass_0, &b_ttbar_mass);
   fChain0->SetBranchAddress("llbar_pt", &llbar_pt_0, &b_llbar_pt);
   fChain0->SetBranchAddress("llbar_phi", &llbar_phi_0, &b_llbar_phi);
   fChain0->SetBranchAddress("llbar_rapidity", &llbar_rapidity_0, &b_llbar_rapidity);
   // fChain0->SetBranchAddress("llbar_arapidity", &llbar_arapidity_0, &b_llbar_arapidity);
   fChain0->SetBranchAddress("llbar_delta_phi", &llbar_delta_phi_0, &b_llbar_delta_phi);
   fChain0->SetBranchAddress("llbar_delta_eta", &llbar_delta_eta_0, &b_llbar_delta_eta);
   fChain0->SetBranchAddress("llbar_delta_rapidity", &llbar_delta_rapidity_0, &b_llbar_delta_rapidity);
   fChain0->SetBranchAddress("llbar_mass", &llbar_mass_0, &b_llbar_mass);

   fChain0->SetBranchAddress("all_mass", &all_mass_0, &b_all_mass);
   fChain0->SetBranchAddress("r_mass", &r_mass_0, &b_r_mass);
   fChain0->SetBranchAddress("jet_multiplicity", &jet_multiplicity_0, &b_jet_multiplicity);
   fChain0->SetBranchAddress("bjet_multiplicity", &bjet_multiplicity_0, &b_bjet_multiplicity);
   fChain0->SetBranchAddress("x1", &x1_0, &b_x1);
   fChain0->SetBranchAddress("x2", &x2_0, &b_x2);

   fChain0->SetBranchAddress("b1n", &b1n_0, &b_b1n);
   fChain0->SetBranchAddress("b2n", &b2n_0, &b_b2n);
   fChain0->SetBranchAddress("b1r", &b1r_0, &b_b1r);
   fChain0->SetBranchAddress("b2r", &b2r_0, &b_b2r);
   fChain0->SetBranchAddress("b1k", &b1k_0, &b_b1k);
   fChain0->SetBranchAddress("b2k", &b2k_0, &b_b2k);
   fChain0->SetBranchAddress("b1j", &b1j_0, &b_b1j);
   fChain0->SetBranchAddress("b2j", &b2j_0, &b_b2j);
   fChain0->SetBranchAddress("b1q", &b1q_0, &b_b1q);
   fChain0->SetBranchAddress("b2q", &b2q_0, &b_b2q);

   fChain0->SetBranchAddress("c_kk", &c_kk_0, &b_c_kk);
   fChain0->SetBranchAddress("c_rr", &c_rr_0, &b_c_rr);
   fChain0->SetBranchAddress("c_nn", &c_nn_0, &b_c_nn);
   fChain0->SetBranchAddress("c_kj", &c_kj_0, &b_c_kj);
   fChain0->SetBranchAddress("c_rq", &c_rq_0, &b_c_rq);

   fChain0->SetBranchAddress("c_rk", &c_rk_0, &b_c_rk);
   fChain0->SetBranchAddress("c_kr", &c_kr_0, &b_c_kr);
   fChain0->SetBranchAddress("c_nr", &c_nr_0, &b_c_nr);
   fChain0->SetBranchAddress("c_rn", &c_rn_0, &b_c_rn);
   fChain0->SetBranchAddress("c_nk", &c_nk_0, &b_c_nk);
   fChain0->SetBranchAddress("c_kn", &c_kn_0, &b_c_kn);
   fChain0->SetBranchAddress("c_rj", &c_rj_0, &b_c_rj);
   fChain0->SetBranchAddress("c_jr", &c_jr_0, &b_c_jr);

   fChain0->SetBranchAddress("ll_cHel", &ll_cHel_0, &b_ll_cHel);
   fChain0->SetBranchAddress("ll_cLab", &ll_cLab_0, &b_ll_cLab);
   fChain0->SetBranchAddress("ll_kNorm", &ll_kNorm_0, &b_ll_kNorm);
   fChain0->SetBranchAddress("ll_rNorm", &ll_rNorm_0, &b_ll_rNorm);

   fChain0->SetBranchAddress("top_scatteringangle_ttbarframe", &top_scatteringangle_ttbarframe_0, &b_top_scatteringangle_ttbarframe);

   fChain0->SetBranchAddress("gen_top_pt", &gen_top_pt_0, &b_gen_top_pt);
   fChain0->SetBranchAddress("gen_top_phi", &gen_top_phi_0, &b_gen_top_phi);
   fChain0->SetBranchAddress("gen_top_rapidity", &gen_top_rapidity_0, &b_gen_top_rapidity);
   // fChain0->SetBranchAddress("gen_top_arapidity", &gen_top_arapidity_0, &b_gen_top_arapidity);
   fChain0->SetBranchAddress("gen_tbar_pt", &gen_tbar_pt_0, &b_gen_tbar_pt);
   fChain0->SetBranchAddress("gen_tbar_phi", &gen_tbar_phi_0, &b_gen_tbar_phi);
   fChain0->SetBranchAddress("gen_tbar_rapidity", &gen_tbar_rapidity_0, &b_gen_tbar_rapidity);
   // fChain0->SetBranchAddress("gen_tbar_arapidity", &gen_tbar_arapidity_0, &b_gen_tbar_arapidity);
   fChain0->SetBranchAddress("gen_l_pt", &gen_l_pt_0, &b_gen_l_pt);
   fChain0->SetBranchAddress("gen_lbar_pt", &gen_lbar_pt_0, &b_gen_lbar_pt);
   fChain0->SetBranchAddress("gen_l_pdgid", &gen_l_pdgid_0, &b_gen_l_pdgid);
   fChain0->SetBranchAddress("gen_lbar_pdgid", &gen_lbar_pdgid_0, &b_gen_lbar_pdgid);
   // fChain0->SetBranchAddress("gen_e_pt", &gen_e_pt_0, &b_gen_e_pt);
   // fChain0->SetBranchAddress("gen_ebar_pt", &gen_ebar_pt_0, &b_gen_ebar_pt);
   // fChain0->SetBranchAddress("gen_mu_pt", &gen_mu_pt_0, &b_gen_mu_pt);
   // fChain0->SetBranchAddress("gen_mubar_pt", &gen_mubar_pt_0, &b_gen_mubar_pt);
   fChain0->SetBranchAddress("gen_b_pt", &gen_b_pt_0, &b_gen_b_pt);
   fChain0->SetBranchAddress("gen_bbar_pt", &gen_bbar_pt_0, &b_gen_bbar_pt);
   fChain0->SetBranchAddress("gen_nu_pt", &gen_nu_pt_0, &b_gen_nu_pt);
   fChain0->SetBranchAddress("gen_nubar_pt", &gen_nubar_pt_0, &b_gen_nubar_pt);
   fChain0->SetBranchAddress("gen_l_eta", &gen_l_eta_0, &b_gen_l_eta);
   fChain0->SetBranchAddress("gen_lbar_eta", &gen_lbar_eta_0, &b_gen_lbar_eta);
   // fChain0->SetBranchAddress("gen_e_eta", &gen_e_eta_0, &b_gen_e_eta);
   // fChain0->SetBranchAddress("gen_ebar_eta", &gen_ebar_eta_0, &b_gen_ebar_eta);
   // fChain0->SetBranchAddress("gen_mu_eta", &gen_mu_eta_0, &b_gen_mu_eta);
   // fChain0->SetBranchAddress("gen_mubar_eta", &gen_mubar_eta_0, &b_gen_mubar_eta);
   fChain0->SetBranchAddress("gen_b_eta", &gen_b_eta_0, &b_gen_b_eta);
   fChain0->SetBranchAddress("gen_bbar_eta", &gen_bbar_eta_0, &b_gen_bbar_eta);
   fChain0->SetBranchAddress("gen_nu_eta", &gen_nu_eta_0, &b_gen_nu_eta);
   fChain0->SetBranchAddress("gen_nubar_eta", &gen_nubar_eta_0, &b_gen_nubar_eta);
   fChain0->SetBranchAddress("gen_l_phi", &gen_l_phi_0, &b_gen_l_phi);
   fChain0->SetBranchAddress("gen_lbar_phi", &gen_lbar_phi_0, &b_gen_lbar_phi);
   // fChain0->SetBranchAddress("gen_e_phi", &gen_e_phi_0, &b_gen_e_phi);
   // fChain0->SetBranchAddress("gen_ebar_phi", &gen_ebar_phi_0, &b_gen_ebar_phi);
   // fChain0->SetBranchAddress("gen_mu_phi", &gen_mu_phi_0, &b_gen_mu_phi);
   // fChain0->SetBranchAddress("gen_mubar_phi", &gen_mubar_phi_0, &b_gen_mubar_phi);
   fChain0->SetBranchAddress("gen_b_phi", &gen_b_phi_0, &b_gen_b_phi);
   fChain0->SetBranchAddress("gen_bbar_phi", &gen_bbar_phi_0, &b_gen_bbar_phi);
   fChain0->SetBranchAddress("gen_nu_phi", &gen_nu_phi_0, &b_gen_nu_phi);
   fChain0->SetBranchAddress("gen_nubar_phi", &gen_nubar_phi_0, &b_gen_nubar_phi);
   fChain0->SetBranchAddress("gen_ttbar_pt", &gen_ttbar_pt_0, &b_gen_ttbar_pt);
   fChain0->SetBranchAddress("gen_ttbar_phi", &gen_ttbar_phi_0, &b_gen_ttbar_phi);
   fChain0->SetBranchAddress("gen_ttbar_rapidity", &gen_ttbar_rapidity_0, &b_gen_ttbar_rapidity);
   // fChain0->SetBranchAddress("gen_ttbar_arapidity", &gen_ttbar_arapidity_0, &b_gen_ttbar_arapidity);
   fChain0->SetBranchAddress("gen_ttbar_delta_phi", &gen_ttbar_delta_phi_0, &b_gen_ttbar_delta_phi);
   fChain0->SetBranchAddress("gen_ttbar_delta_eta", &gen_ttbar_delta_eta_0, &b_gen_ttbar_delta_eta);
   fChain0->SetBranchAddress("gen_ttbar_delta_rapidity", &gen_ttbar_delta_rapidity_0, &b_gen_ttbar_delta_rapidity);
   fChain0->SetBranchAddress("gen_ttbar_mass", &gen_ttbar_mass_0, &b_gen_ttbar_mass);
   fChain0->SetBranchAddress("gen_llbar_pt", &gen_llbar_pt_0, &b_gen_llbar_pt);
   fChain0->SetBranchAddress("gen_llbar_phi", &gen_llbar_phi_0, &b_gen_llbar_phi);
   fChain0->SetBranchAddress("gen_llbar_rapidity", &gen_llbar_rapidity_0, &b_gen_llbar_rapidity);
   // fChain0->SetBranchAddress("gen_llbar_arapidity", &gen_llbar_arapidity_0, &b_gen_llbar_arapidity);
   fChain0->SetBranchAddress("gen_llbar_delta_phi", &gen_llbar_delta_phi_0, &b_gen_llbar_delta_phi);
   fChain0->SetBranchAddress("gen_llbar_delta_eta", &gen_llbar_delta_eta_0, &b_gen_llbar_delta_eta);
   fChain0->SetBranchAddress("gen_llbar_delta_rapidity", &gen_llbar_delta_rapidity_0, &b_gen_llbar_delta_rapidity);
   fChain0->SetBranchAddress("gen_llbar_mass", &gen_llbar_mass_0, &b_gen_llbar_mass);

   fChain0->SetBranchAddress("gen_all_mass", &gen_all_mass_0, &b_gen_all_mass);
   fChain0->SetBranchAddress("gen_r_mass", &gen_r_mass_0, &b_gen_r_mass);
   fChain0->SetBranchAddress("gen_jet_multiplicity", &gen_jet_multiplicity_0, &b_gen_jet_multiplicity);
   fChain0->SetBranchAddress("gen_x1", &gen_x1_0, &b_gen_x1);
   fChain0->SetBranchAddress("gen_x2", &gen_x2_0, &b_gen_x2);

   fChain0->SetBranchAddress("gen_b1n", &gen_b1n_0, &b_gen_b1n);
   fChain0->SetBranchAddress("gen_b2n", &gen_b2n_0, &b_gen_b2n);
   fChain0->SetBranchAddress("gen_b1r", &gen_b1r_0, &b_gen_b1r);
   fChain0->SetBranchAddress("gen_b2r", &gen_b2r_0, &b_gen_b2r);
   fChain0->SetBranchAddress("gen_b1k", &gen_b1k_0, &b_gen_b1k);
   fChain0->SetBranchAddress("gen_b2k", &gen_b2k_0, &b_gen_b2k);
   fChain0->SetBranchAddress("gen_b1j", &gen_b1j_0, &b_gen_b1j);
   fChain0->SetBranchAddress("gen_b2j", &gen_b2j_0, &b_gen_b2j);
   fChain0->SetBranchAddress("gen_b1q", &gen_b1q_0, &b_gen_b1q);
   fChain0->SetBranchAddress("gen_b2q", &gen_b2q_0, &b_gen_b2q);

   fChain0->SetBranchAddress("gen_c_kk", &gen_c_kk_0, &b_gen_c_kk);
   fChain0->SetBranchAddress("gen_c_rr", &gen_c_rr_0, &b_gen_c_rr);
   fChain0->SetBranchAddress("gen_c_nn", &gen_c_nn_0, &b_gen_c_nn);
   fChain0->SetBranchAddress("gen_c_kj", &gen_c_kj_0, &b_gen_c_kj);
   fChain0->SetBranchAddress("gen_c_rq", &gen_c_rq_0, &b_gen_c_rq);

   fChain0->SetBranchAddress("gen_c_rk", &gen_c_rk_0, &b_gen_c_rk);
   fChain0->SetBranchAddress("gen_c_kr", &gen_c_kr_0, &b_gen_c_kr);
   fChain0->SetBranchAddress("gen_c_nr", &gen_c_nr_0, &b_gen_c_nr);
   fChain0->SetBranchAddress("gen_c_rn", &gen_c_rn_0, &b_gen_c_rn);
   fChain0->SetBranchAddress("gen_c_nk", &gen_c_nk_0, &b_gen_c_nk);
   fChain0->SetBranchAddress("gen_c_kn", &gen_c_kn_0, &b_gen_c_kn);
   fChain0->SetBranchAddress("gen_c_rj", &gen_c_rj_0, &b_gen_c_rj);
   fChain0->SetBranchAddress("gen_c_jr", &gen_c_jr_0, &b_gen_c_jr);

   fChain0->SetBranchAddress("gen_ll_cHel", &gen_ll_cHel_0, &b_gen_ll_cHel);
   fChain0->SetBranchAddress("gen_ll_cLab", &gen_ll_cLab_0, &b_gen_ll_cLab);
   fChain0->SetBranchAddress("gen_ll_kNorm", &gen_ll_kNorm_0, &b_gen_ll_kNorm);
   fChain0->SetBranchAddress("gen_ll_rNorm", &gen_ll_rNorm_0, &b_gen_ll_rNorm);

   fChain0->SetBranchAddress("gen_top_scatteringangle_ttbarframe", &gen_top_scatteringangle_ttbarframe_0, &b_gen_top_scatteringangle_ttbarframe);

   fChain0->SetBranchAddress("entry", &entry_0, &b_entry);
   fChain0->SetBranchAddress("isTopGen", &isTopGen_0, &b_isTopGen);
   fChain0->SetBranchAddress("isKinReco", &isKinReco_0, &b_isKinReco);
   fChain0->SetBranchAddress("trueLevelWeight", &trueLevelWeight_0, &b_trueLevelWeight);
   fChain0->SetBranchAddress("eventWeight", &eventWeight_0, &b_eventWeight);


   Notify();
}

Bool_t makeTUnfoldHisto::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void makeTUnfoldHisto::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t makeTUnfoldHisto::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef makeTUnfoldHisto_cxx
