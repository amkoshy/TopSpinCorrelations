#define makeTUnfoldHisto_cxx
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TXMLDocument.h>
#include <math.h>
#include <stdio.h>

#include <fstream>
#include <iostream>
#include <string>

#include "TClass.h"
#include "TKey.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TUnfold.h"
#include "TUnfoldBinning.h"
#include "TUnfoldBinningXML.h"
#include "TUnfoldDensity.h"
#include "TUnfoldSys.h"
#include "TVectorD.h"
#include "makeTUnfoldHisto.h"

/* Declare histograms to be used for unfolding */
namespace unfoldhistos {

/**
 * @brief Whether or not to create bootstrap pseudoexperiment histograms.
 * Changing doBootstrap to kFalse seems to break things
 */
const bool doBootstrap = kTRUE;

/**
 * @brief Number of pseudoexperiments to perform
 * Amandeep : Changing nPE to zero to just get the response matrices
 * Original : const int nPE = 10000;
 */
const int nPE = 1;

/**
 * @brief Number of bins in the residual histograms
 */
int nbinsres = 80;

// TODO: Creating a struct to store these binnings might be a good idea

// ********
// BINNINGS
// ********

// 1D binnings
// Kinematic variables
TUnfoldBinning *gen_top_pt_binning, *reco_top_pt_binning, *residual_top_pt_binning;
TUnfoldBinning *gen_l_pt_binning, *reco_l_pt_binning, *residual_l_pt_binning;
TUnfoldBinning *gen_lbar_pt_binning, *reco_lbar_pt_binning, *residual_lbar_pt_binning;

TUnfoldBinning *gen_ttbar_pt_binning, *reco_ttbar_pt_binning, *residual_ttbar_pt_binning;
TUnfoldBinning *gen_ttbar_mass_binning, *reco_ttbar_mass_binning, *residual_ttbar_mass_binning;
TUnfoldBinning *gen_ttbar_delta_phi_binning, *reco_ttbar_delta_phi_binning, *residual_ttbar_delta_phi_binning;
TUnfoldBinning *gen_ttbar_delta_eta_binning, *reco_ttbar_delta_eta_binning, *residual_ttbar_delta_eta_binning;
TUnfoldBinning *gen_ttbar_rapidity_binning, *reco_ttbar_rapidity_binning, *residual_ttbar_rapidity_binning;
TUnfoldBinning *gen_llbar_pt_binning, *reco_llbar_pt_binning, *residual_llbar_pt_binning;
TUnfoldBinning *gen_llbar_mass_binning, *reco_llbar_mass_binning, *residual_llbar_mass_binning;

// Spin corr variables
TUnfoldBinning *gen_b1k_binning, *reco_b1k_binning, *residual_b1k_binning;
TUnfoldBinning *gen_b2k_binning, *reco_b2k_binning, *residual_b2k_binning;
TUnfoldBinning *gen_b1j_binning, *reco_b1j_binning, *residual_b1j_binning;
TUnfoldBinning *gen_b2j_binning, *reco_b2j_binning, *residual_b2j_binning;
TUnfoldBinning *gen_b1r_binning, *reco_b1r_binning, *residual_b1r_binning;
TUnfoldBinning *gen_b2r_binning, *reco_b2r_binning, *residual_b2r_binning;
TUnfoldBinning *gen_b1q_binning, *reco_b1q_binning, *residual_b1q_binning;
TUnfoldBinning *gen_b2q_binning, *reco_b2q_binning, *residual_b2q_binning;
TUnfoldBinning *gen_b1n_binning, *reco_b1n_binning, *residual_b1n_binning;
TUnfoldBinning *gen_b2n_binning, *reco_b2n_binning, *residual_b2n_binning;

TUnfoldBinning *gen_c_kk_binning, *reco_c_kk_binning, *residual_c_kk_binning;
TUnfoldBinning *gen_c_rr_binning, *reco_c_rr_binning, *residual_c_rr_binning;
TUnfoldBinning *gen_c_nn_binning, *reco_c_nn_binning, *residual_c_nn_binning;
TUnfoldBinning *gen_c_kj_binning, *reco_c_kj_binning, *residual_c_kj_binning;
TUnfoldBinning *gen_c_rq_binning, *reco_c_rq_binning, *residual_c_rq_binning;

TUnfoldBinning *gen_c_Prk_binning, *reco_c_Prk_binning, *residual_c_Prk_binning;
TUnfoldBinning *gen_c_Mrk_binning, *reco_c_Mrk_binning, *residual_c_Mrk_binning;
TUnfoldBinning *gen_c_Pnr_binning, *reco_c_Pnr_binning, *residual_c_Pnr_binning;
TUnfoldBinning *gen_c_Mnr_binning, *reco_c_Mnr_binning, *residual_c_Mnr_binning;
TUnfoldBinning *gen_c_Pnk_binning, *reco_c_Pnk_binning, *residual_c_Pnk_binning;
TUnfoldBinning *gen_c_Mnk_binning, *reco_c_Mnk_binning, *residual_c_Mnk_binning;
TUnfoldBinning *gen_c_Prj_binning, *reco_c_Prj_binning, *residual_c_Prj_binning;
TUnfoldBinning *gen_c_Mrj_binning, *reco_c_Mrj_binning, *residual_c_Mrj_binning;

TUnfoldBinning *gen_c_han_binning, *reco_c_han_binning, *residual_c_han_binning;
TUnfoldBinning *gen_c_sca_binning, *reco_c_sca_binning, *residual_c_sca_binning;
TUnfoldBinning *gen_c_tra_binning, *reco_c_tra_binning, *residual_c_tra_binning;
TUnfoldBinning *gen_c_kjL_binning, *reco_c_kjL_binning, *residual_c_kjL_binning;
TUnfoldBinning *gen_c_rqL_binning, *reco_c_rqL_binning, *residual_c_rqL_binning;

TUnfoldBinning *gen_c_rkP_binning, *reco_c_rkP_binning, *residual_c_rkP_binning;
TUnfoldBinning *gen_c_rkM_binning, *reco_c_rkM_binning, *residual_c_rkM_binning;
TUnfoldBinning *gen_c_nrP_binning, *reco_c_nrP_binning, *residual_c_nrP_binning;
TUnfoldBinning *gen_c_nrM_binning, *reco_c_nrM_binning, *residual_c_nrM_binning;
TUnfoldBinning *gen_c_nkP_binning, *reco_c_nkP_binning, *residual_c_nkP_binning;
TUnfoldBinning *gen_c_nkM_binning, *reco_c_nkM_binning, *residual_c_nkM_binning;

TUnfoldBinning *gen_c_rk_binning, *reco_c_rk_binning, *residual_c_rk_binning;
TUnfoldBinning *gen_c_kr_binning, *reco_c_kr_binning, *residual_c_kr_binning;
TUnfoldBinning *gen_c_nr_binning, *reco_c_nr_binning, *residual_c_nr_binning;
TUnfoldBinning *gen_c_rn_binning, *reco_c_rn_binning, *residual_c_rn_binning;
TUnfoldBinning *gen_c_nk_binning, *reco_c_nk_binning, *residual_c_nk_binning;
TUnfoldBinning *gen_c_kn_binning, *reco_c_kn_binning, *residual_c_kn_binning;
TUnfoldBinning *gen_c_rj_binning, *reco_c_rj_binning, *residual_c_rj_binning;
TUnfoldBinning *gen_c_jr_binning, *reco_c_jr_binning, *residual_c_jr_binning;

TUnfoldBinning *gen_ll_cHel_binning, *reco_ll_cHel_binning, *residual_ll_cHel_binning;
TUnfoldBinning *gen_ll_cLab_binning, *reco_ll_cLab_binning, *residual_ll_cLab_binning;
TUnfoldBinning *gen_ll_kNorm_binning, *reco_ll_kNorm_binning, *residual_ll_kNorm_binning;
TUnfoldBinning *gen_ll_rNorm_binning, *reco_ll_rNorm_binning, *residual_ll_rNorm_binning;
TUnfoldBinning *gen_llbar_delta_phi_binning, *reco_llbar_delta_phi_binning, *residual_llbar_delta_phi_binning;
TUnfoldBinning *gen_llbar_delta_eta_binning, *reco_llbar_delta_eta_binning, *residual_llbar_delta_eta_binning;

// 2D binnings

// ******
// mttbar
// ******

// Amandeep : 2D binnings for spin corr variables
// TUnfoldBinning *gen_c_kk_mttbar_binning, *reco_c_kk_mttbar_binning,
// *residual_c_kk_mttbar_binning;
TUnfoldBinning *gen_b1k_exjet_binning, *reco_b1k_exjet_binning, *residual_b1k_exjet_binning;
TUnfoldBinning *gen_b2k_exjet_binning, *reco_b2k_exjet_binning, *residual_b2k_exjet_binning;
TUnfoldBinning *gen_b1j_exjet_binning, *reco_b1j_exjet_binning, *residual_b1j_exjet_binning;
TUnfoldBinning *gen_b2j_exjet_binning, *reco_b2j_exjet_binning, *residual_b2j_exjet_binning;
TUnfoldBinning *gen_b1r_exjet_binning, *reco_b1r_exjet_binning, *residual_b1r_exjet_binning;
TUnfoldBinning *gen_b2r_exjet_binning, *reco_b2r_exjet_binning, *residual_b2r_exjet_binning;
TUnfoldBinning *gen_b1q_exjet_binning, *reco_b1q_exjet_binning, *residual_b1q_exjet_binning;
TUnfoldBinning *gen_b2q_exjet_binning, *reco_b2q_exjet_binning, *residual_b2q_exjet_binning;
TUnfoldBinning *gen_b1n_exjet_binning, *reco_b1n_exjet_binning, *residual_b1n_exjet_binning;
TUnfoldBinning *gen_b2n_exjet_binning, *reco_b2n_exjet_binning, *residual_b2n_exjet_binning;

TUnfoldBinning *gen_c_kk_exjet_binning, *reco_c_kk_exjet_binning, *residual_c_kk_exjet_binning;
TUnfoldBinning *gen_c_rr_exjet_binning, *reco_c_rr_exjet_binning, *residual_c_rr_exjet_binning;
TUnfoldBinning *gen_c_nn_exjet_binning, *reco_c_nn_exjet_binning, *residual_c_nn_exjet_binning;
TUnfoldBinning *gen_c_kj_exjet_binning, *reco_c_kj_exjet_binning, *residual_c_kj_exjet_binning;
TUnfoldBinning *gen_c_rq_exjet_binning, *reco_c_rq_exjet_binning, *residual_c_rq_exjet_binning;

TUnfoldBinning *gen_c_rk_exjet_binning, *reco_c_rk_exjet_binning, *residual_c_rk_exjet_binning;
TUnfoldBinning *gen_c_kr_exjet_binning, *reco_c_kr_exjet_binning, *residual_c_kr_exjet_binning;
TUnfoldBinning *gen_c_nr_exjet_binning, *reco_c_nr_exjet_binning, *residual_c_nr_exjet_binning;
TUnfoldBinning *gen_c_rn_exjet_binning, *reco_c_rn_exjet_binning, *residual_c_rn_exjet_binning;
TUnfoldBinning *gen_c_nk_exjet_binning, *reco_c_nk_exjet_binning, *residual_c_nk_exjet_binning;
TUnfoldBinning *gen_c_kn_exjet_binning, *reco_c_kn_exjet_binning, *residual_c_kn_exjet_binning;
TUnfoldBinning *gen_c_rj_exjet_binning, *reco_c_rj_exjet_binning, *residual_c_rj_exjet_binning;
TUnfoldBinning *gen_c_jr_exjet_binning, *reco_c_jr_exjet_binning, *residual_c_jr_exjet_binning;

TUnfoldBinning *gen_c_Prk_exjet_binning, *reco_c_Prk_exjet_binning, *residual_c_Prk_exjet_binning;
TUnfoldBinning *gen_c_Mrk_exjet_binning, *reco_c_Mrk_exjet_binning, *residual_c_Mrk_exjet_binning;
TUnfoldBinning *gen_c_Pnr_exjet_binning, *reco_c_Pnr_exjet_binning, *residual_c_Pnr_exjet_binning;
TUnfoldBinning *gen_c_Mnr_exjet_binning, *reco_c_Mnr_exjet_binning, *residual_c_Mnr_exjet_binning;
TUnfoldBinning *gen_c_Pnk_exjet_binning, *reco_c_Pnk_exjet_binning, *residual_c_Pnk_exjet_binning;
TUnfoldBinning *gen_c_Mnk_exjet_binning, *reco_c_Mnk_exjet_binning, *residual_c_Mnk_exjet_binning;
TUnfoldBinning *gen_c_Prj_exjet_binning, *reco_c_Prj_exjet_binning, *residual_c_Prj_exjet_binning;
TUnfoldBinning *gen_c_Mrj_exjet_binning, *reco_c_Mrj_exjet_binning, *residual_c_Mrj_exjet_binning;

TUnfoldBinning *gen_c_han_exjet_binning, *reco_c_han_exjet_binning, *residual_c_han_exjet_binning;
TUnfoldBinning *gen_c_sca_exjet_binning, *reco_c_sca_exjet_binning, *residual_c_sca_exjet_binning;
TUnfoldBinning *gen_c_tra_exjet_binning, *reco_c_tra_exjet_binning, *residual_c_tra_exjet_binning;
TUnfoldBinning *gen_c_kjL_exjet_binning, *reco_c_kjL_exjet_binning, *residual_c_kjL_exjet_binning;
TUnfoldBinning *gen_c_rqL_exjet_binning, *reco_c_rqL_exjet_binning, *residual_c_rqL_exjet_binning;

TUnfoldBinning *gen_c_rkP_exjet_binning, *reco_c_rkP_exjet_binning, *residual_c_rkP_exjet_binning;
TUnfoldBinning *gen_c_rkM_exjet_binning, *reco_c_rkM_exjet_binning, *residual_c_rkM_exjet_binning;
TUnfoldBinning *gen_c_nrP_exjet_binning, *reco_c_nrP_exjet_binning, *residual_c_nrP_exjet_binning;
TUnfoldBinning *gen_c_nrM_exjet_binning, *reco_c_nrM_exjet_binning, *residual_c_nrM_exjet_binning;
TUnfoldBinning *gen_c_nkP_exjet_binning, *reco_c_nkP_exjet_binning, *residual_c_nkP_exjet_binning;
TUnfoldBinning *gen_c_nkM_exjet_binning, *reco_c_nkM_exjet_binning, *residual_c_nkM_exjet_binning;

TUnfoldBinning *gen_ll_cHel_exjet_binning, *reco_ll_cHel_exjet_binning, *residual_ll_cHel_exjet_binning;
TUnfoldBinning *gen_ll_cLab_exjet_binning, *reco_ll_cLab_exjet_binning, *residual_ll_cLab_exjet_binning;
TUnfoldBinning *gen_ll_kNorm_exjet_binning, *reco_ll_kNorm_exjet_binning, *residual_ll_kNorm_exjet_binning;
TUnfoldBinning *gen_ll_rNorm_exjet_binning, *reco_ll_rNorm_exjet_binning, *residual_ll_rNorm_exjet_binning;
TUnfoldBinning *gen_llbar_delta_phi_exjet_binning, *reco_llbar_delta_phi_exjet_binning,
    *residual_llbar_delta_phi_exjet_binning;
TUnfoldBinning *gen_llbar_delta_eta_exjet_binning, *reco_llbar_delta_eta_exjet_binning,
    *residual_llbar_delta_eta_exjet_binning;

// ******
// top_pt
// ******

// Amandeep : 2D binnings for spin corr variables
// TUnfoldBinning *gen_c_kk_top_pt_binning, *reco_c_kk_top_pt_binning,
// *residual_c_kk_top_pt_binning;

// TUnfoldBinning *gen_b1k_top_pt_binning, *reco_b1k_top_pt_binning,
// *residual_b1k_top_pt_binning; TUnfoldBinning *gen_b2k_top_pt_binning,
// *reco_b2k_top_pt_binning, *residual_b2k_top_pt_binning; TUnfoldBinning
// *gen_b1j_top_pt_binning, *reco_b1j_top_pt_binning,
// *residual_b1j_top_pt_binning; TUnfoldBinning *gen_b2j_top_pt_binning,
// *reco_b2j_top_pt_binning, *residual_b2j_top_pt_binning; TUnfoldBinning
// *gen_b1r_top_pt_binning, *reco_b1r_top_pt_binning,
// *residual_b1r_top_pt_binning; TUnfoldBinning *gen_b2r_top_pt_binning,
// *reco_b2r_top_pt_binning, *residual_b2r_top_pt_binning; TUnfoldBinning
// *gen_b1q_top_pt_binning, *reco_b1q_top_pt_binning,
// *residual_b1q_top_pt_binning; TUnfoldBinning *gen_b2q_top_pt_binning,
// *reco_b2q_top_pt_binning, *residual_b2q_top_pt_binning; TUnfoldBinning
// *gen_b1n_top_pt_binning, *reco_b1n_top_pt_binning,
// *residual_b1n_top_pt_binning; TUnfoldBinning *gen_b2n_top_pt_binning,
// *reco_b2n_top_pt_binning, *residual_b2n_top_pt_binning;

// TUnfoldBinning *gen_c_kk_top_pt_binning, *reco_c_kk_top_pt_binning,
// *residual_c_kk_top_pt_binning; TUnfoldBinning *gen_c_rr_top_pt_binning,
// *reco_c_rr_top_pt_binning, *residual_c_rr_top_pt_binning; TUnfoldBinning
// *gen_c_nn_top_pt_binning, *reco_c_nn_top_pt_binning,
// *residual_c_nn_top_pt_binning; TUnfoldBinning *gen_c_kj_top_pt_binning,
// *reco_c_kj_top_pt_binning, *residual_c_kj_top_pt_binning; TUnfoldBinning
// *gen_c_rq_top_pt_binning, *reco_c_rq_top_pt_binning,
// *residual_c_rq_top_pt_binning;

// TUnfoldBinning *gen_c_rk_top_pt_binning, *reco_c_rk_top_pt_binning,
// *residual_c_rk_top_pt_binning; TUnfoldBinning *gen_c_kr_top_pt_binning,
// *reco_c_kr_top_pt_binning, *residual_c_kr_top_pt_binning; TUnfoldBinning
// *gen_c_nr_top_pt_binning, *reco_c_nr_top_pt_binning,
// *residual_c_nr_top_pt_binning; TUnfoldBinning *gen_c_rn_top_pt_binning,
// *reco_c_rn_top_pt_binning, *residual_c_rn_top_pt_binning; TUnfoldBinning
// *gen_c_nk_top_pt_binning, *reco_c_nk_top_pt_binning,
// *residual_c_nk_top_pt_binning; TUnfoldBinning *gen_c_kn_top_pt_binning,
// *reco_c_kn_top_pt_binning, *residual_c_kn_top_pt_binning; TUnfoldBinning
// *gen_c_rj_top_pt_binning, *reco_c_rj_top_pt_binning,
// *residual_c_rj_top_pt_binning; TUnfoldBinning *gen_c_jr_top_pt_binning,
// *reco_c_jr_top_pt_binning, *residual_c_jr_top_pt_binning;

// TUnfoldBinning *gen_c_Prk_top_pt_binning, *reco_c_Prk_top_pt_binning,
// *residual_c_Prk_top_pt_binning; TUnfoldBinning *gen_c_Mrk_top_pt_binning,
// *reco_c_Mrk_top_pt_binning, *residual_c_Mrk_top_pt_binning; TUnfoldBinning
// *gen_c_Pnr_top_pt_binning, *reco_c_Pnr_top_pt_binning,
// *residual_c_Pnr_top_pt_binning; TUnfoldBinning *gen_c_Mnr_top_pt_binning,
// *reco_c_Mnr_top_pt_binning, *residual_c_Mnr_top_pt_binning; TUnfoldBinning
// *gen_c_Pnk_top_pt_binning, *reco_c_Pnk_top_pt_binning,
// *residual_c_Pnk_top_pt_binning; TUnfoldBinning *gen_c_Mnk_top_pt_binning,
// *reco_c_Mnk_top_pt_binning, *residual_c_Mnk_top_pt_binning; TUnfoldBinning
// *gen_c_Prj_top_pt_binning, *reco_c_Prj_top_pt_binning,
// *residual_c_Prj_top_pt_binning; TUnfoldBinning *gen_c_Mrj_top_pt_binning,
// *reco_c_Mrj_top_pt_binning, *residual_c_Mrj_top_pt_binning;

// TUnfoldBinning *gen_c_han_top_pt_binning, *reco_c_han_top_pt_binning,
// *residual_c_han_top_pt_binning; TUnfoldBinning *gen_c_sca_top_pt_binning,
// *reco_c_sca_top_pt_binning, *residual_c_sca_top_pt_binning; TUnfoldBinning
// *gen_c_tra_top_pt_binning, *reco_c_tra_top_pt_binning,
// *residual_c_tra_top_pt_binning; TUnfoldBinning *gen_c_kjL_top_pt_binning,
// *reco_c_kjL_top_pt_binning, *residual_c_kjL_top_pt_binning; TUnfoldBinning
// *gen_c_rqL_top_pt_binning, *reco_c_rqL_top_pt_binning,
// *residual_c_rqL_top_pt_binning;

// TUnfoldBinning *gen_c_rkP_top_pt_binning, *reco_c_rkP_top_pt_binning,
// *residual_c_rkP_top_pt_binning; TUnfoldBinning *gen_c_rkM_top_pt_binning,
// *reco_c_rkM_top_pt_binning, *residual_c_rkM_top_pt_binning; TUnfoldBinning
// *gen_c_nrP_top_pt_binning, *reco_c_nrP_top_pt_binning,
// *residual_c_nrP_top_pt_binning; TUnfoldBinning *gen_c_nrM_top_pt_binning,
// *reco_c_nrM_top_pt_binning, *residual_c_nrM_top_pt_binning; TUnfoldBinning
// *gen_c_nkP_top_pt_binning, *reco_c_nkP_top_pt_binning,
// *residual_c_nkP_top_pt_binning; TUnfoldBinning *gen_c_nkM_top_pt_binning,
// *reco_c_nkM_top_pt_binning, *residual_c_nkM_top_pt_binning;

// TUnfoldBinning *gen_ll_cHel_top_pt_binning, *reco_ll_cHel_top_pt_binning,
// *residual_ll_cHel_top_pt_binning; TUnfoldBinning *gen_ll_cLab_top_pt_binning,
// *reco_ll_cLab_top_pt_binning, *residual_ll_cLab_top_pt_binning;
// TUnfoldBinning *gen_ll_kNorm_top_pt_binning, *reco_ll_kNorm_top_pt_binning,
// *residual_ll_kNorm_top_pt_binning; TUnfoldBinning
// *gen_ll_rNorm_top_pt_binning, *reco_ll_rNorm_top_pt_binning,
// *residual_ll_rNorm_top_pt_binning; TUnfoldBinning
// *gen_llbar_delta_phi_top_pt_binning, *reco_llbar_delta_phi_top_pt_binning,
// *residual_llbar_delta_phi_top_pt_binning; TUnfoldBinning
// *gen_llbar_delta_eta_top_pt_binning, *reco_llbar_delta_eta_top_pt_binning,
// *residual_llbar_delta_eta_top_pt_binning;

// ********************
// top_scatteringangle_ttbarframe
// ********************

// Amandeep : 2D binnings for spin corr variables
// TUnfoldBinning *gen_c_kk_top_scatteringangle_ttbarframe_binning,
// *reco_c_kk_top_scatteringangle_ttbarframe_binning,
// *residual_c_kk_top_scatteringangle_ttbarframe_binning;

// TUnfoldBinning *gen_b1k_top_scatteringangle_ttbarframe_binning,
// *reco_b1k_top_scatteringangle_ttbarframe_binning,
// *residual_b1k_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_b2k_top_scatteringangle_ttbarframe_binning,
// *reco_b2k_top_scatteringangle_ttbarframe_binning,
// *residual_b2k_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_b1j_top_scatteringangle_ttbarframe_binning,
// *reco_b1j_top_scatteringangle_ttbarframe_binning,
// *residual_b1j_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_b2j_top_scatteringangle_ttbarframe_binning,
// *reco_b2j_top_scatteringangle_ttbarframe_binning,
// *residual_b2j_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_b1r_top_scatteringangle_ttbarframe_binning,
// *reco_b1r_top_scatteringangle_ttbarframe_binning,
// *residual_b1r_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_b2r_top_scatteringangle_ttbarframe_binning,
// *reco_b2r_top_scatteringangle_ttbarframe_binning,
// *residual_b2r_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_b1q_top_scatteringangle_ttbarframe_binning,
// *reco_b1q_top_scatteringangle_ttbarframe_binning,
// *residual_b1q_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_b2q_top_scatteringangle_ttbarframe_binning,
// *reco_b2q_top_scatteringangle_ttbarframe_binning,
// *residual_b2q_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_b1n_top_scatteringangle_ttbarframe_binning,
// *reco_b1n_top_scatteringangle_ttbarframe_binning,
// *residual_b1n_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_b2n_top_scatteringangle_ttbarframe_binning,
// *reco_b2n_top_scatteringangle_ttbarframe_binning,
// *residual_b2n_top_scatteringangle_ttbarframe_binning;

// TUnfoldBinning *gen_c_kk_top_scatteringangle_ttbarframe_binning,
// *reco_c_kk_top_scatteringangle_ttbarframe_binning,
// *residual_c_kk_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_rr_top_scatteringangle_ttbarframe_binning,
// *reco_c_rr_top_scatteringangle_ttbarframe_binning,
// *residual_c_rr_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_nn_top_scatteringangle_ttbarframe_binning,
// *reco_c_nn_top_scatteringangle_ttbarframe_binning,
// *residual_c_nn_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_kj_top_scatteringangle_ttbarframe_binning,
// *reco_c_kj_top_scatteringangle_ttbarframe_binning,
// *residual_c_kj_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_rq_top_scatteringangle_ttbarframe_binning,
// *reco_c_rq_top_scatteringangle_ttbarframe_binning,
// *residual_c_rq_top_scatteringangle_ttbarframe_binning;

// TUnfoldBinning *gen_c_Prk_top_scatteringangle_ttbarframe_binning,
// *reco_c_Prk_top_scatteringangle_ttbarframe_binning,
// *residual_c_Prk_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_Mrk_top_scatteringangle_ttbarframe_binning,
// *reco_c_Mrk_top_scatteringangle_ttbarframe_binning,
// *residual_c_Mrk_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_Pnr_top_scatteringangle_ttbarframe_binning,
// *reco_c_Pnr_top_scatteringangle_ttbarframe_binning,
// *residual_c_Pnr_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_Mnr_top_scatteringangle_ttbarframe_binning,
// *reco_c_Mnr_top_scatteringangle_ttbarframe_binning,
// *residual_c_Mnr_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_Pnk_top_scatteringangle_ttbarframe_binning,
// *reco_c_Pnk_top_scatteringangle_ttbarframe_binning,
// *residual_c_Pnk_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_Mnk_top_scatteringangle_ttbarframe_binning,
// *reco_c_Mnk_top_scatteringangle_ttbarframe_binning,
// *residual_c_Mnk_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_Prj_top_scatteringangle_ttbarframe_binning,
// *reco_c_Prj_top_scatteringangle_ttbarframe_binning,
// *residual_c_Prj_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_Mrj_top_scatteringangle_ttbarframe_binning,
// *reco_c_Mrj_top_scatteringangle_ttbarframe_binning,
// *residual_c_Mrj_top_scatteringangle_ttbarframe_binning;

// TUnfoldBinning *gen_c_han_top_scatteringangle_ttbarframe_binning,
// *reco_c_han_top_scatteringangle_ttbarframe_binning,
// *residual_c_han_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_sca_top_scatteringangle_ttbarframe_binning,
// *reco_c_sca_top_scatteringangle_ttbarframe_binning,
// *residual_c_sca_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_tra_top_scatteringangle_ttbarframe_binning,
// *reco_c_tra_top_scatteringangle_ttbarframe_binning,
// *residual_c_tra_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_kjL_top_scatteringangle_ttbarframe_binning,
// *reco_c_kjL_top_scatteringangle_ttbarframe_binning,
// *residual_c_kjL_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_rqL_top_scatteringangle_ttbarframe_binning,
// *reco_c_rqL_top_scatteringangle_ttbarframe_binning,
// *residual_c_rqL_top_scatteringangle_ttbarframe_binning;

// TUnfoldBinning *gen_c_rkP_top_scatteringangle_ttbarframe_binning,
// *reco_c_rkP_top_scatteringangle_ttbarframe_binning,
// *residual_c_rkP_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_rkM_top_scatteringangle_ttbarframe_binning,
// *reco_c_rkM_top_scatteringangle_ttbarframe_binning,
// *residual_c_rkM_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_nrP_top_scatteringangle_ttbarframe_binning,
// *reco_c_nrP_top_scatteringangle_ttbarframe_binning,
// *residual_c_nrP_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_nrM_top_scatteringangle_ttbarframe_binning,
// *reco_c_nrM_top_scatteringangle_ttbarframe_binning,
// *residual_c_nrM_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_nkP_top_scatteringangle_ttbarframe_binning,
// *reco_c_nkP_top_scatteringangle_ttbarframe_binning,
// *residual_c_nkP_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_nkM_top_scatteringangle_ttbarframe_binning,
// *reco_c_nkM_top_scatteringangle_ttbarframe_binning,
// *residual_c_nkM_top_scatteringangle_ttbarframe_binning;

// TUnfoldBinning *gen_c_rk_top_scatteringangle_ttbarframe_binning,
// *reco_c_rk_top_scatteringangle_ttbarframe_binning,
// *residual_c_rk_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_kr_top_scatteringangle_ttbarframe_binning,
// *reco_c_kr_top_scatteringangle_ttbarframe_binning,
// *residual_c_kr_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_nr_top_scatteringangle_ttbarframe_binning,
// *reco_c_nr_top_scatteringangle_ttbarframe_binning,
// *residual_c_nr_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_rn_top_scatteringangle_ttbarframe_binning,
// *reco_c_rn_top_scatteringangle_ttbarframe_binning,
// *residual_c_rn_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_nk_top_scatteringangle_ttbarframe_binning,
// *reco_c_nk_top_scatteringangle_ttbarframe_binning,
// *residual_c_nk_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_kn_top_scatteringangle_ttbarframe_binning,
// *reco_c_kn_top_scatteringangle_ttbarframe_binning,
// *residual_c_kn_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_rj_top_scatteringangle_ttbarframe_binning,
// *reco_c_rj_top_scatteringangle_ttbarframe_binning,
// *residual_c_rj_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_c_jr_top_scatteringangle_ttbarframe_binning,
// *reco_c_jr_top_scatteringangle_ttbarframe_binning,
// *residual_c_jr_top_scatteringangle_ttbarframe_binning;

// TUnfoldBinning *gen_ll_cHel_top_scatteringangle_ttbarframe_binning,
// *reco_ll_cHel_top_scatteringangle_ttbarframe_binning,
// *residual_ll_cHel_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_ll_cLab_top_scatteringangle_ttbarframe_binning,
// *reco_ll_cLab_top_scatteringangle_ttbarframe_binning,
// *residual_ll_cLab_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_ll_kNorm_top_scatteringangle_ttbarframe_binning,
// *reco_ll_kNorm_top_scatteringangle_ttbarframe_binning,
// *residual_ll_kNorm_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_ll_rNorm_top_scatteringangle_ttbarframe_binning,
// *reco_ll_rNorm_top_scatteringangle_ttbarframe_binning,
// *residual_ll_rNorm_top_scatteringangle_ttbarframe_binning; TUnfoldBinning
// *gen_llbar_delta_phi_top_scatteringangle_ttbarframe_binning,
// *reco_llbar_delta_phi_top_scatteringangle_ttbarframe_binning,
// *residual_llbar_delta_phi_top_scatteringangle_ttbarframe_binning;
// TUnfoldBinning *gen_llbar_delta_eta_top_scatteringangle_ttbarframe_binning,
// *reco_llbar_delta_eta_top_scatteringangle_ttbarframe_binning,
// *residual_llbar_delta_eta_top_scatteringangle_ttbarframe_binning;

// ***************
// ***************
// RECO HISTOGRAMS
// ***************
// ***************

// 1D Reco histograms
TH1 *hreco_top_pt, *hreco_l_pt, *hreco_lbar_pt;
TH1 *hreco_ttbar_pt, *hreco_ttbar_mass, *hreco_ttbar_delta_phi, *hreco_ttbar_delta_eta, *hreco_ttbar_rapidity,
    *hreco_llbar_pt, *hreco_llbar_mass;
TH1 *hreco_b1k, *hreco_b2k, *hreco_b1j, *hreco_b2j, *hreco_b1r, *hreco_b2r, *hreco_b1q, *hreco_b2q, *hreco_b1n,
    *hreco_b2n;
TH1 *hreco_c_kk, *hreco_c_rr, *hreco_c_nn, *hreco_c_kj, *hreco_c_rq;
TH1 *hreco_c_Prk, *hreco_c_Mrk, *hreco_c_Pnr, *hreco_c_Mnr, *hreco_c_Pnk, *hreco_c_Mnk, *hreco_c_Prj, *hreco_c_Mrj;
TH1 *hreco_c_han, *hreco_c_sca, *hreco_c_tra, *hreco_c_kjL, *hreco_c_rqL;
TH1 *hreco_c_rkP, *hreco_c_rkM, *hreco_c_nrP, *hreco_c_nrM, *hreco_c_nkP, *hreco_c_nkM;
TH1 *hreco_c_rk, *hreco_c_kr, *hreco_c_nr, *hreco_c_rn, *hreco_c_nk, *hreco_c_kn, *hreco_c_rj, *hreco_c_jr;
TH1 *hreco_ll_cHel, *hreco_ll_cLab, *hreco_ll_kNorm, *hreco_ll_rNorm, *hreco_llbar_delta_phi, *hreco_llbar_delta_eta;

// 2D Reco histograms
// Amandeep : Adding 2D reco histograms

// ******
// mttbar
// ******
TH1 *hreco_b1k_exjet, *hreco_b2k_exjet, *hreco_b1j_exjet, *hreco_b2j_exjet, *hreco_b1r_exjet, *hreco_b2r_exjet,
    *hreco_b1q_exjet, *hreco_b2q_exjet, *hreco_b1n_exjet, *hreco_b2n_exjet;
TH1 *hreco_c_kk_exjet, *hreco_c_rr_exjet, *hreco_c_nn_exjet, *hreco_c_kj_exjet, *hreco_c_rq_exjet;
TH1 *hreco_c_Prk_exjet, *hreco_c_Mrk_exjet, *hreco_c_Pnr_exjet, *hreco_c_Mnr_exjet, *hreco_c_Pnk_exjet,
    *hreco_c_Mnk_exjet, *hreco_c_Prj_exjet, *hreco_c_Mrj_exjet;
TH1 *hreco_c_han_exjet, *hreco_c_sca_exjet, *hreco_c_tra_exjet, *hreco_c_kjL_exjet, *hreco_c_rqL_exjet;
TH1 *hreco_c_rkP_exjet, *hreco_c_rkM_exjet, *hreco_c_nrP_exjet, *hreco_c_nrM_exjet, *hreco_c_nkP_exjet,
    *hreco_c_nkM_exjet;
TH1 *hreco_c_rk_exjet, *hreco_c_kr_exjet, *hreco_c_nr_exjet, *hreco_c_rn_exjet, *hreco_c_nk_exjet, *hreco_c_kn_exjet,
    *hreco_c_rj_exjet, *hreco_c_jr_exjet;
TH1 *hreco_ll_cHel_exjet, *hreco_ll_cLab_exjet, *hreco_ll_kNorm_exjet, *hreco_ll_rNorm_exjet,
    *hreco_llbar_delta_phi_exjet, *hreco_llbar_delta_eta_exjet;

// ******
// top_pt
// ******
// TH1 *hreco_b1k_top_pt  , *hreco_b2k_top_pt, *hreco_b1j_top_pt,
// *hreco_b2j_top_pt, *hreco_b1r_top_pt, *hreco_b2r_top_pt, *hreco_b1q_top_pt,
// *hreco_b2q_top_pt, *hreco_b1n_top_pt, *hreco_b2n_top_pt; TH1
// *hreco_c_kk_top_pt , *hreco_c_rr_top_pt, *hreco_c_nn_top_pt,
// *hreco_c_kj_top_pt, *hreco_c_rq_top_pt; TH1 *hreco_c_Prk_top_pt,
// *hreco_c_Mrk_top_pt, *hreco_c_Pnr_top_pt, *hreco_c_Mnr_top_pt,
// *hreco_c_Pnk_top_pt, *hreco_c_Mnk_top_pt, *hreco_c_Prj_top_pt,
// *hreco_c_Mrj_top_pt; TH1 *hreco_c_han_top_pt, * hreco_c_sca_top_pt,
// *hreco_c_tra_top_pt, *hreco_c_kjL_top_pt, *hreco_c_rqL_top_pt; TH1
// *hreco_c_rkP_top_pt, *hreco_c_rkM_top_pt, *hreco_c_nrP_top_pt,
// *hreco_c_nrM_top_pt, *hreco_c_nkP_top_pt, *hreco_c_nkM_top_pt; TH1
// *hreco_c_rk_top_pt , *hreco_c_kr_top_pt, *hreco_c_nr_top_pt,
// *hreco_c_rn_top_pt, *hreco_c_nk_top_pt, *hreco_c_kn_top_pt,
// *hreco_c_rj_top_pt, *hreco_c_jr_top_pt; TH1 *hreco_ll_cHel_top_pt,
// *hreco_ll_cLab_top_pt, *hreco_ll_kNorm_top_pt, *hreco_ll_rNorm_top_pt,
// *hreco_llbar_delta_phi_top_pt, *hreco_llbar_delta_eta_top_pt;

// // ********************
// // top_scatteringangle_ttbarframe
// // ********************
// TH1 *hreco_b1k_top_scatteringangle_ttbarframe  ,
// *hreco_b2k_top_scatteringangle_ttbarframe,
// *hreco_b1j_top_scatteringangle_ttbarframe,
// *hreco_b2j_top_scatteringangle_ttbarframe,
// *hreco_b1r_top_scatteringangle_ttbarframe,
// *hreco_b2r_top_scatteringangle_ttbarframe,
// *hreco_b1q_top_scatteringangle_ttbarframe,
// *hreco_b2q_top_scatteringangle_ttbarframe,
// *hreco_b1n_top_scatteringangle_ttbarframe,
// *hreco_b2n_top_scatteringangle_ttbarframe; TH1
// *hreco_c_kk_top_scatteringangle_ttbarframe ,
// *hreco_c_rr_top_scatteringangle_ttbarframe,
// *hreco_c_nn_top_scatteringangle_ttbarframe,
// *hreco_c_kj_top_scatteringangle_ttbarframe,
// *hreco_c_rq_top_scatteringangle_ttbarframe; TH1
// *hreco_c_Prk_top_scatteringangle_ttbarframe,
// *hreco_c_Mrk_top_scatteringangle_ttbarframe,
// *hreco_c_Pnr_top_scatteringangle_ttbarframe,
// *hreco_c_Mnr_top_scatteringangle_ttbarframe,
// *hreco_c_Pnk_top_scatteringangle_ttbarframe,
// *hreco_c_Mnk_top_scatteringangle_ttbarframe,
// *hreco_c_Prj_top_scatteringangle_ttbarframe,
// *hreco_c_Mrj_top_scatteringangle_ttbarframe; TH1
// *hreco_c_han_top_scatteringangle_ttbarframe, *
// hreco_c_sca_top_scatteringangle_ttbarframe,
// *hreco_c_tra_top_scatteringangle_ttbarframe,
// *hreco_c_kjL_top_scatteringangle_ttbarframe,
// *hreco_c_rqL_top_scatteringangle_ttbarframe; TH1
// *hreco_c_rkP_top_scatteringangle_ttbarframe,
// *hreco_c_rkM_top_scatteringangle_ttbarframe,
// *hreco_c_nrP_top_scatteringangle_ttbarframe,
// *hreco_c_nrM_top_scatteringangle_ttbarframe,
// *hreco_c_nkP_top_scatteringangle_ttbarframe,
// *hreco_c_nkM_top_scatteringangle_ttbarframe; TH1
// *hreco_c_rk_top_scatteringangle_ttbarframe ,
// *hreco_c_kr_top_scatteringangle_ttbarframe,
// *hreco_c_nr_top_scatteringangle_ttbarframe,
// *hreco_c_rn_top_scatteringangle_ttbarframe,
// *hreco_c_nk_top_scatteringangle_ttbarframe,
// *hreco_c_kn_top_scatteringangle_ttbarframe,
// *hreco_c_rj_top_scatteringangle_ttbarframe,
// *hreco_c_jr_top_scatteringangle_ttbarframe; TH1
// *hreco_ll_cHel_top_scatteringangle_ttbarframe,
// *hreco_ll_cLab_top_scatteringangle_ttbarframe,
// *hreco_ll_kNorm_top_scatteringangle_ttbarframe,
// *hreco_ll_rNorm_top_scatteringangle_ttbarframe,
// *hreco_llbar_delta_phi_top_scatteringangle_ttbarframe,
// *hreco_llbar_delta_eta_top_scatteringangle_ttbarframe;

// *****************
// *****************
// PSUEDO HISTOGRAMS
// *****************
// *****************

// 1D pseudo histograms
TH1 *hrecoBootstrap_top_pt[nPE], *hrecoBootstrap_l_pt[nPE], *hrecoBootstrap_lbar_pt[nPE];
TH1 *hrecoBootstrap_ttbar_pt[nPE], *hrecoBootstrap_ttbar_mass[nPE], *hrecoBootstrap_ttbar_delta_phi[nPE],
    *hrecoBootstrap_ttbar_delta_eta[nPE], *hrecoBootstrap_ttbar_rapidity[nPE], *hrecoBootstrap_llbar_pt[nPE],
    *hrecoBootstrap_llbar_mass[nPE];
TH1 *hrecoBootstrap_b1k[nPE], *hrecoBootstrap_b2k[nPE], *hrecoBootstrap_b1j[nPE], *hrecoBootstrap_b2j[nPE],
    *hrecoBootstrap_b1r[nPE], *hrecoBootstrap_b2r[nPE], *hrecoBootstrap_b1q[nPE], *hrecoBootstrap_b2q[nPE],
    *hrecoBootstrap_b1n[nPE], *hrecoBootstrap_b2n[nPE];
TH1 *hrecoBootstrap_c_kk[nPE], *hrecoBootstrap_c_rr[nPE], *hrecoBootstrap_c_nn[nPE], *hrecoBootstrap_c_kj[nPE],
    *hrecoBootstrap_c_rq[nPE];
TH1 *hrecoBootstrap_c_Prk[nPE], *hrecoBootstrap_c_Mrk[nPE], *hrecoBootstrap_c_Pnr[nPE], *hrecoBootstrap_c_Mnr[nPE],
    *hrecoBootstrap_c_Pnk[nPE], *hrecoBootstrap_c_Mnk[nPE], *hrecoBootstrap_c_Prj[nPE], *hrecoBootstrap_c_Mrj[nPE];
TH1 *hrecoBootstrap_c_rk[nPE], *hrecoBootstrap_c_kr[nPE], *hrecoBootstrap_c_nr[nPE], *hrecoBootstrap_c_rn[nPE],
    *hrecoBootstrap_c_nk[nPE], *hrecoBootstrap_c_kn[nPE], *hrecoBootstrap_c_rj[nPE], *hrecoBootstrap_c_jr[nPE];
TH1 *hrecoBootstrap_c_han[nPE], *hrecoBootstrap_c_sca[nPE], *hrecoBootstrap_c_tra[nPE], *hrecoBootstrap_c_kjL[nPE],
    *hrecoBootstrap_c_rqL[nPE];
TH1 *hrecoBootstrap_c_rkP[nPE], *hrecoBootstrap_c_rkM[nPE], *hrecoBootstrap_c_nrP[nPE], *hrecoBootstrap_c_nrM[nPE],
    *hrecoBootstrap_c_nkP[nPE], *hrecoBootstrap_c_nkM[nPE];
TH1 *hrecoBootstrap_ll_cHel[nPE], *hrecoBootstrap_ll_cLab[nPE], *hrecoBootstrap_ll_kNorm[nPE],
    *hrecoBootstrap_ll_rNorm[nPE], *hrecoBootstrap_llbar_delta_phi[nPE], *hrecoBootstrap_llbar_delta_eta[nPE];

// // 2D pseudo histograms
// Amanadeep : Adding 2D psuedo histograms

// ******
// mttbar
// ******
TH1 *hrecoBootstrap_b1k_exjet[nPE], *hrecoBootstrap_b2k_exjet[nPE], *hrecoBootstrap_b1j_exjet[nPE],
    *hrecoBootstrap_b2j_exjet[nPE], *hrecoBootstrap_b1r_exjet[nPE], *hrecoBootstrap_b2r_exjet[nPE],
    *hrecoBootstrap_b1q_exjet[nPE], *hrecoBootstrap_b2q_exjet[nPE], *hrecoBootstrap_b1n_exjet[nPE],
    *hrecoBootstrap_b2n_exjet[nPE];
TH1 *hrecoBootstrap_c_kk_exjet[nPE], *hrecoBootstrap_c_rr_exjet[nPE], *hrecoBootstrap_c_nn_exjet[nPE],
    *hrecoBootstrap_c_kj_exjet[nPE], *hrecoBootstrap_c_rq_exjet[nPE];
TH1 *hrecoBootstrap_c_Prk_exjet[nPE], *hrecoBootstrap_c_Mrk_exjet[nPE], *hrecoBootstrap_c_Pnr_exjet[nPE],
    *hrecoBootstrap_c_Mnr_exjet[nPE], *hrecoBootstrap_c_Pnk_exjet[nPE], *hrecoBootstrap_c_Mnk_exjet[nPE],
    *hrecoBootstrap_c_Prj_exjet[nPE], *hrecoBootstrap_c_Mrj_exjet[nPE];
TH1 *hrecoBootstrap_c_rk_exjet[nPE], *hrecoBootstrap_c_kr_exjet[nPE], *hrecoBootstrap_c_nr_exjet[nPE],
    *hrecoBootstrap_c_rn_exjet[nPE], *hrecoBootstrap_c_nk_exjet[nPE], *hrecoBootstrap_c_kn_exjet[nPE],
    *hrecoBootstrap_c_rj_exjet[nPE], *hrecoBootstrap_c_jr_exjet[nPE];
TH1 *hrecoBootstrap_c_han_exjet[nPE], *hrecoBootstrap_c_sca_exjet[nPE], *hrecoBootstrap_c_tra_exjet[nPE],
    *hrecoBootstrap_c_kjL_exjet[nPE], *hrecoBootstrap_c_rqL_exjet[nPE];
TH1 *hrecoBootstrap_c_rkP_exjet[nPE], *hrecoBootstrap_c_rkM_exjet[nPE], *hrecoBootstrap_c_nrP_exjet[nPE],
    *hrecoBootstrap_c_nrM_exjet[nPE], *hrecoBootstrap_c_nkP_exjet[nPE], *hrecoBootstrap_c_nkM_exjet[nPE];
TH1 *hrecoBootstrap_ll_cHel_exjet[nPE], *hrecoBootstrap_ll_cLab_exjet[nPE], *hrecoBootstrap_ll_kNorm_exjet[nPE],
    *hrecoBootstrap_ll_rNorm_exjet[nPE], *hrecoBootstrap_llbar_delta_phi_exjet[nPE],
    *hrecoBootstrap_llbar_delta_eta_exjet[nPE];

// ******
// top_pt
// ******
// TH1 *hrecoBootstrap_b1k_top_pt[nPE]  , *hrecoBootstrap_b2k_top_pt[nPE],
// *hrecoBootstrap_b1j_top_pt[nPE], *hrecoBootstrap_b2j_top_pt[nPE],
// *hrecoBootstrap_b1r_top_pt[nPE], *hrecoBootstrap_b2r_top_pt[nPE],
// *hrecoBootstrap_b1q_top_pt[nPE], *hrecoBootstrap_b2q_top_pt[nPE],
// *hrecoBootstrap_b1n_top_pt[nPE], *hrecoBootstrap_b2n_top_pt[nPE]; TH1
// *hrecoBootstrap_c_kk_top_pt[nPE] , *hrecoBootstrap_c_rr_top_pt[nPE],
// *hrecoBootstrap_c_nn_top_pt[nPE], *hrecoBootstrap_c_kj_top_pt[nPE],
// *hrecoBootstrap_c_rq_top_pt[nPE]; TH1 *hrecoBootstrap_c_Prk_top_pt[nPE],
// *hrecoBootstrap_c_Mrk_top_pt[nPE], *hrecoBootstrap_c_Pnr_top_pt[nPE],
// *hrecoBootstrap_c_Mnr_top_pt[nPE], *hrecoBootstrap_c_Pnk_top_pt[nPE],
// *hrecoBootstrap_c_Mnk_top_pt[nPE], *hrecoBootstrap_c_Prj_top_pt[nPE],
// *hrecoBootstrap_c_Mrj_top_pt[nPE]; TH1 *hrecoBootstrap_c_rk_top_pt[nPE] ,
// *hrecoBootstrap_c_kr_top_pt[nPE], *hrecoBootstrap_c_nr_top_pt[nPE],
// *hrecoBootstrap_c_rn_top_pt[nPE], *hrecoBootstrap_c_nk_top_pt[nPE],
// *hrecoBootstrap_c_kn_top_pt[nPE], *hrecoBootstrap_c_rj_top_pt[nPE],
// *hrecoBootstrap_c_jr_top_pt[nPE]; TH1 *hrecoBootstrap_c_han_top_pt[nPE],
// *hrecoBootstrap_c_sca_top_pt[nPE], *hrecoBootstrap_c_tra_top_pt[nPE],
// *hrecoBootstrap_c_kjL_top_pt[nPE], *hrecoBootstrap_c_rqL_top_pt[nPE]; TH1
// *hrecoBootstrap_c_rkP_top_pt[nPE], *hrecoBootstrap_c_rkM_top_pt[nPE],
// *hrecoBootstrap_c_nrP_top_pt[nPE], *hrecoBootstrap_c_nrM_top_pt[nPE],
// *hrecoBootstrap_c_nkP_top_pt[nPE], *hrecoBootstrap_c_nkM_top_pt[nPE]; TH1
// *hrecoBootstrap_ll_cHel_top_pt[nPE], *hrecoBootstrap_ll_cLab_top_pt[nPE],
// *hrecoBootstrap_ll_kNorm_top_pt[nPE], *hrecoBootstrap_ll_rNorm_top_pt[nPE],
// *hrecoBootstrap_llbar_delta_phi_top_pt[nPE],
// *hrecoBootstrap_llbar_delta_eta_top_pt[nPE];

// // ********************
// // top_scatteringangle_ttbarframe
// // ********************
// TH1 *hrecoBootstrap_b1k_top_scatteringangle_ttbarframe[nPE]  ,
// *hrecoBootstrap_b2k_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_b1j_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_b2j_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_b1r_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_b2r_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_b1q_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_b2q_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_b1n_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_b2n_top_scatteringangle_ttbarframe[nPE]; TH1
// *hrecoBootstrap_c_kk_top_scatteringangle_ttbarframe[nPE] ,
// *hrecoBootstrap_c_rr_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_nn_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_kj_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_rq_top_scatteringangle_ttbarframe[nPE]; TH1
// *hrecoBootstrap_c_Prk_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_Mrk_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_Pnr_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_Mnr_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_Pnk_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_Mnk_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_Prj_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_Mrj_top_scatteringangle_ttbarframe[nPE]; TH1
// *hrecoBootstrap_c_rk_top_scatteringangle_ttbarframe[nPE] ,
// *hrecoBootstrap_c_kr_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_nr_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_rn_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_nk_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_kn_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_rj_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_jr_top_scatteringangle_ttbarframe[nPE]; TH1
// *hrecoBootstrap_c_han_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_sca_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_tra_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_kjL_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_rqL_top_scatteringangle_ttbarframe[nPE]; TH1
// *hrecoBootstrap_c_rkP_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_rkM_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_nrP_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_nrM_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_nkP_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_c_nkM_top_scatteringangle_ttbarframe[nPE]; TH1
// *hrecoBootstrap_ll_cHel_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_ll_cLab_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_ll_kNorm_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_ll_rNorm_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_llbar_delta_phi_top_scatteringangle_ttbarframe[nPE],
// *hrecoBootstrap_llbar_delta_eta_top_scatteringangle_ttbarframe[nPE];

// **************
// **************
// GEN HISTOGRAMS
// **************
// **************

// 1D gen histograms
TH1 *hgen_top_pt, *hgen_l_pt, *hgen_lbar_pt;
TH1 *hgen_ttbar_pt, *hgen_ttbar_mass, *hgen_ttbar_delta_phi, *hgen_ttbar_delta_eta, *hgen_ttbar_rapidity,
    *hgen_llbar_pt, *hgen_llbar_mass;
TH1 *hgen_b1k, *hgen_b2k, *hgen_b1j, *hgen_b2j, *hgen_b1r, *hgen_b2r, *hgen_b1q, *hgen_b2q, *hgen_b1n, *hgen_b2n;
TH1 *hgen_c_kk, *hgen_c_rr, *hgen_c_nn, *hgen_c_kj, *hgen_c_rq;
TH1 *hgen_c_Prk, *hgen_c_Mrk, *hgen_c_Pnr, *hgen_c_Mnr, *hgen_c_Pnk, *hgen_c_Mnk, *hgen_c_Prj, *hgen_c_Mrj;
TH1 *hgen_c_rk, *hgen_c_kr, *hgen_c_nr, *hgen_c_rn, *hgen_c_nk, *hgen_c_kn, *hgen_c_rj, *hgen_c_jr;
TH1 *hgen_c_han, *hgen_c_sca, *hgen_c_tra, *hgen_c_kjL, *hgen_c_rqL;
TH1 *hgen_c_rkP, *hgen_c_rkM, *hgen_c_nrP, *hgen_c_nrM, *hgen_c_nkP, *hgen_c_nkM;
TH1 *hgen_ll_cHel, *hgen_ll_cLab, *hgen_ll_kNorm, *hgen_ll_rNorm, *hgen_llbar_delta_phi, *hgen_llbar_delta_eta;

// 2D gen histograms
// Amandeep : Adding 2D gen histograms

// ******
// mttbar
// ******
TH1 *hgen_b1k_exjet, *hgen_b2k_exjet, *hgen_b1j_exjet, *hgen_b2j_exjet, *hgen_b1r_exjet, *hgen_b2r_exjet,
    *hgen_b1q_exjet, *hgen_b2q_exjet, *hgen_b1n_exjet, *hgen_b2n_exjet;
TH1 *hgen_c_kk_exjet, *hgen_c_rr_exjet, *hgen_c_nn_exjet, *hgen_c_kj_exjet, *hgen_c_rq_exjet;
TH1 *hgen_c_Prk_exjet, *hgen_c_Mrk_exjet, *hgen_c_Pnr_exjet, *hgen_c_Mnr_exjet, *hgen_c_Pnk_exjet, *hgen_c_Mnk_exjet,
    *hgen_c_Prj_exjet, *hgen_c_Mrj_exjet;
TH1 *hgen_c_rk_exjet, *hgen_c_kr_exjet, *hgen_c_nr_exjet, *hgen_c_rn_exjet, *hgen_c_nk_exjet, *hgen_c_kn_exjet,
    *hgen_c_rj_exjet, *hgen_c_jr_exjet;
TH1 *hgen_c_han_exjet, *hgen_c_sca_exjet, *hgen_c_tra_exjet, *hgen_c_kjL_exjet, *hgen_c_rqL_exjet;
TH1 *hgen_c_rkP_exjet, *hgen_c_rkM_exjet, *hgen_c_nrP_exjet, *hgen_c_nrM_exjet, *hgen_c_nkP_exjet, *hgen_c_nkM_exjet;
TH1 *hgen_ll_cHel_exjet, *hgen_ll_cLab_exjet, *hgen_ll_kNorm_exjet, *hgen_ll_rNorm_exjet, *hgen_llbar_delta_phi_exjet,
    *hgen_llbar_delta_eta_exjet;

// ******
// top_pt
// ******
// TH1 *hgen_b1k_top_pt  , *hgen_b2k_top_pt, *hgen_b1j_top_pt, *hgen_b2j_top_pt,
// *hgen_b1r_top_pt, *hgen_b2r_top_pt, *hgen_b1q_top_pt, *hgen_b2q_top_pt,
// *hgen_b1n_top_pt, *hgen_b2n_top_pt; TH1 *hgen_c_kk_top_pt ,
// *hgen_c_rr_top_pt, *hgen_c_nn_top_pt, *hgen_c_kj_top_pt, *hgen_c_rq_top_pt;
// TH1 *hgen_c_Prk_top_pt, *hgen_c_Mrk_top_pt, *hgen_c_Pnr_top_pt,
// *hgen_c_Mnr_top_pt, *hgen_c_Pnk_top_pt, *hgen_c_Mnk_top_pt,
// *hgen_c_Prj_top_pt, *hgen_c_Mrj_top_pt; TH1 *hgen_c_rk_top_pt ,
// *hgen_c_kr_top_pt, *hgen_c_nr_top_pt, *hgen_c_rn_top_pt, *hgen_c_nk_top_pt,
// *hgen_c_kn_top_pt, *hgen_c_rj_top_pt, *hgen_c_jr_top_pt; TH1
// *hgen_c_han_top_pt, *hgen_c_sca_top_pt, *hgen_c_tra_top_pt,
// *hgen_c_kjL_top_pt, *hgen_c_rqL_top_pt; TH1 *hgen_c_rkP_top_pt,
// *hgen_c_rkM_top_pt, *hgen_c_nrP_top_pt, *hgen_c_nrM_top_pt,
// *hgen_c_nkP_top_pt, *hgen_c_nkM_top_pt; TH1 *hgen_ll_cHel_top_pt,
// *hgen_ll_cLab_top_pt, *hgen_ll_kNorm_top_pt, *hgen_ll_rNorm_top_pt,
// *hgen_llbar_delta_phi_top_pt, *hgen_llbar_delta_eta_top_pt;

// // ********************
// // top_scatteringangle_ttbarframe
// // ********************
// TH1 *hgen_b1k_top_scatteringangle_ttbarframe  ,
// *hgen_b2k_top_scatteringangle_ttbarframe,
// *hgen_b1j_top_scatteringangle_ttbarframe,
// *hgen_b2j_top_scatteringangle_ttbarframe,
// *hgen_b1r_top_scatteringangle_ttbarframe,
// *hgen_b2r_top_scatteringangle_ttbarframe,
// *hgen_b1q_top_scatteringangle_ttbarframe,
// *hgen_b2q_top_scatteringangle_ttbarframe,
// *hgen_b1n_top_scatteringangle_ttbarframe,
// *hgen_b2n_top_scatteringangle_ttbarframe; TH1
// *hgen_c_kk_top_scatteringangle_ttbarframe ,
// *hgen_c_rr_top_scatteringangle_ttbarframe,
// *hgen_c_nn_top_scatteringangle_ttbarframe,
// *hgen_c_kj_top_scatteringangle_ttbarframe,
// *hgen_c_rq_top_scatteringangle_ttbarframe; TH1
// *hgen_c_Prk_top_scatteringangle_ttbarframe,
// *hgen_c_Mrk_top_scatteringangle_ttbarframe,
// *hgen_c_Pnr_top_scatteringangle_ttbarframe,
// *hgen_c_Mnr_top_scatteringangle_ttbarframe,
// *hgen_c_Pnk_top_scatteringangle_ttbarframe,
// *hgen_c_Mnk_top_scatteringangle_ttbarframe,
// *hgen_c_Prj_top_scatteringangle_ttbarframe,
// *hgen_c_Mrj_top_scatteringangle_ttbarframe; TH1
// *hgen_c_rk_top_scatteringangle_ttbarframe ,
// *hgen_c_kr_top_scatteringangle_ttbarframe,
// *hgen_c_nr_top_scatteringangle_ttbarframe,
// *hgen_c_rn_top_scatteringangle_ttbarframe,
// *hgen_c_nk_top_scatteringangle_ttbarframe,
// *hgen_c_kn_top_scatteringangle_ttbarframe,
// *hgen_c_rj_top_scatteringangle_ttbarframe,
// *hgen_c_jr_top_scatteringangle_ttbarframe; TH1
// *hgen_c_han_top_scatteringangle_ttbarframe,
// *hgen_c_sca_top_scatteringangle_ttbarframe,
// *hgen_c_tra_top_scatteringangle_ttbarframe,
// *hgen_c_kjL_top_scatteringangle_ttbarframe,
// *hgen_c_rqL_top_scatteringangle_ttbarframe; TH1
// *hgen_c_rkP_top_scatteringangle_ttbarframe,
// *hgen_c_rkM_top_scatteringangle_ttbarframe,
// *hgen_c_nrP_top_scatteringangle_ttbarframe,
// *hgen_c_nrM_top_scatteringangle_ttbarframe,
// *hgen_c_nkP_top_scatteringangle_ttbarframe,
// *hgen_c_nkM_top_scatteringangle_ttbarframe; TH1
// *hgen_ll_cHel_top_scatteringangle_ttbarframe,
// *hgen_ll_cLab_top_scatteringangle_ttbarframe,
// *hgen_ll_kNorm_top_scatteringangle_ttbarframe,
// *hgen_ll_rNorm_top_scatteringangle_ttbarframe,
// *hgen_llbar_delta_phi_top_scatteringangle_ttbarframe,
// *hgen_llbar_delta_eta_top_scatteringangle_ttbarframe;

// *****************
// *****************
// VISGEN HISTOGRAMS
// *****************
// *****************

// 1D visible gen histograms
TH1 *hvisgen_top_pt, *hvisgen_l_pt, *hvisgen_lbar_pt;
TH1 *hvisgen_ttbar_pt, *hvisgen_ttbar_mass, *hvisgen_ttbar_delta_phi, *hvisgen_ttbar_delta_eta, *hvisgen_ttbar_rapidity,
    *hvisgen_llbar_pt, *hvisgen_llbar_mass;
TH1 *hvisgen_b1k, *hvisgen_b2k, *hvisgen_b1j, *hvisgen_b2j, *hvisgen_b1r, *hvisgen_b2r, *hvisgen_b1q, *hvisgen_b2q,
    *hvisgen_b1n, *hvisgen_b2n;
TH1 *hvisgen_c_kk, *hvisgen_c_rr, *hvisgen_c_nn, *hvisgen_c_kj, *hvisgen_c_rq;
TH1 *hvisgen_c_Prk, *hvisgen_c_Mrk, *hvisgen_c_Pnr, *hvisgen_c_Mnr, *hvisgen_c_Pnk, *hvisgen_c_Mnk, *hvisgen_c_Prj,
    *hvisgen_c_Mrj;
TH1 *hvisgen_c_rk, *hvisgen_c_kr, *hvisgen_c_nr, *hvisgen_c_rn, *hvisgen_c_nk, *hvisgen_c_kn, *hvisgen_c_rj,
    *hvisgen_c_jr;
TH1 *hvisgen_c_han, *hvisgen_c_sca, *hvisgen_c_tra, *hvisgen_c_kjL, *hvisgen_c_rqL;
TH1 *hvisgen_c_rkP, *hvisgen_c_rkM, *hvisgen_c_nrP, *hvisgen_c_nrM, *hvisgen_c_nkP, *hvisgen_c_nkM;
TH1 *hvisgen_ll_cHel, *hvisgen_ll_cLab, *hvisgen_ll_kNorm, *hvisgen_ll_rNorm, *hvisgen_llbar_delta_phi,
    *hvisgen_llbar_delta_eta;

// 2D visible gen histograms
// Amandeep : Adding 2D visible gen histograms

// ******
// mttbar
// ******
TH1 *hvisgen_b1k_exjet, *hvisgen_b2k_exjet, *hvisgen_b1j_exjet, *hvisgen_b2j_exjet, *hvisgen_b1r_exjet,
    *hvisgen_b2r_exjet, *hvisgen_b1q_exjet, *hvisgen_b2q_exjet, *hvisgen_b1n_exjet, *hvisgen_b2n_exjet;
TH1 *hvisgen_c_kk_exjet, *hvisgen_c_rr_exjet, *hvisgen_c_nn_exjet, *hvisgen_c_kj_exjet, *hvisgen_c_rq_exjet;
TH1 *hvisgen_c_Prk_exjet, *hvisgen_c_Mrk_exjet, *hvisgen_c_Pnr_exjet, *hvisgen_c_Mnr_exjet, *hvisgen_c_Pnk_exjet,
    *hvisgen_c_Mnk_exjet, *hvisgen_c_Prj_exjet, *hvisgen_c_Mrj_exjet;
TH1 *hvisgen_c_rk_exjet, *hvisgen_c_kr_exjet, *hvisgen_c_nr_exjet, *hvisgen_c_rn_exjet, *hvisgen_c_nk_exjet,
    *hvisgen_c_kn_exjet, *hvisgen_c_rj_exjet, *hvisgen_c_jr_exjet;
TH1 *hvisgen_c_han_exjet, *hvisgen_c_sca_exjet, *hvisgen_c_tra_exjet, *hvisgen_c_kjL_exjet, *hvisgen_c_rqL_exjet;
TH1 *hvisgen_c_rkP_exjet, *hvisgen_c_rkM_exjet, *hvisgen_c_nrP_exjet, *hvisgen_c_nrM_exjet, *hvisgen_c_nkP_exjet,
    *hvisgen_c_nkM_exjet;
TH1 *hvisgen_ll_cHel_exjet, *hvisgen_ll_cLab_exjet, *hvisgen_ll_kNorm_exjet, *hvisgen_ll_rNorm_exjet,
    *hvisgen_llbar_delta_phi_exjet, *hvisgen_llbar_delta_eta_exjet;

// ******
// top_pt
// ******
// TH1 *hvisgen_b1k_top_pt  , *hvisgen_b2k_top_pt, *hvisgen_b1j_top_pt,
// *hvisgen_b2j_top_pt, *hvisgen_b1r_top_pt, *hvisgen_b2r_top_pt,
// *hvisgen_b1q_top_pt, *hvisgen_b2q_top_pt, *hvisgen_b1n_top_pt,
// *hvisgen_b2n_top_pt; TH1 *hvisgen_c_kk_top_pt , *hvisgen_c_rr_top_pt,
// *hvisgen_c_nn_top_pt, *hvisgen_c_kj_top_pt, *hvisgen_c_rq_top_pt; TH1
// *hvisgen_c_Prk_top_pt, *hvisgen_c_Mrk_top_pt, *hvisgen_c_Pnr_top_pt,
// *hvisgen_c_Mnr_top_pt, *hvisgen_c_Pnk_top_pt, *hvisgen_c_Mnk_top_pt,
// *hvisgen_c_Prj_top_pt, *hvisgen_c_Mrj_top_pt; TH1 *hvisgen_c_rk_top_pt ,
// *hvisgen_c_kr_top_pt, *hvisgen_c_nr_top_pt, *hvisgen_c_rn_top_pt,
// *hvisgen_c_nk_top_pt, *hvisgen_c_kn_top_pt, *hvisgen_c_rj_top_pt,
// *hvisgen_c_jr_top_pt; TH1 *hvisgen_c_han_top_pt, *hvisgen_c_sca_top_pt,
// *hvisgen_c_tra_top_pt, *hvisgen_c_kjL_top_pt, *hvisgen_c_rqL_top_pt; TH1
// *hvisgen_c_rkP_top_pt, *hvisgen_c_rkM_top_pt, *hvisgen_c_nrP_top_pt,
// *hvisgen_c_nrM_top_pt, *hvisgen_c_nkP_top_pt, *hvisgen_c_nkM_top_pt; TH1
// *hvisgen_ll_cHel_top_pt, *hvisgen_ll_cLab_top_pt, *hvisgen_ll_kNorm_top_pt,
// *hvisgen_ll_rNorm_top_pt, *hvisgen_llbar_delta_phi_top_pt,
// *hvisgen_llbar_delta_eta_top_pt;

// // ********************
// // top_scatteringangle_ttbarframe
// // ********************
// TH1 *hvisgen_b1k_top_scatteringangle_ttbarframe  ,
// *hvisgen_b2k_top_scatteringangle_ttbarframe,
// *hvisgen_b1j_top_scatteringangle_ttbarframe,
// *hvisgen_b2j_top_scatteringangle_ttbarframe,
// *hvisgen_b1r_top_scatteringangle_ttbarframe,
// *hvisgen_b2r_top_scatteringangle_ttbarframe,
// *hvisgen_b1q_top_scatteringangle_ttbarframe,
// *hvisgen_b2q_top_scatteringangle_ttbarframe,
// *hvisgen_b1n_top_scatteringangle_ttbarframe,
// *hvisgen_b2n_top_scatteringangle_ttbarframe; TH1
// *hvisgen_c_kk_top_scatteringangle_ttbarframe ,
// *hvisgen_c_rr_top_scatteringangle_ttbarframe,
// *hvisgen_c_nn_top_scatteringangle_ttbarframe,
// *hvisgen_c_kj_top_scatteringangle_ttbarframe,
// *hvisgen_c_rq_top_scatteringangle_ttbarframe; TH1
// *hvisgen_c_Prk_top_scatteringangle_ttbarframe,
// *hvisgen_c_Mrk_top_scatteringangle_ttbarframe,
// *hvisgen_c_Pnr_top_scatteringangle_ttbarframe,
// *hvisgen_c_Mnr_top_scatteringangle_ttbarframe,
// *hvisgen_c_Pnk_top_scatteringangle_ttbarframe,
// *hvisgen_c_Mnk_top_scatteringangle_ttbarframe,
// *hvisgen_c_Prj_top_scatteringangle_ttbarframe,
// *hvisgen_c_Mrj_top_scatteringangle_ttbarframe; TH1
// *hvisgen_c_rk_top_scatteringangle_ttbarframe ,
// *hvisgen_c_kr_top_scatteringangle_ttbarframe,
// *hvisgen_c_nr_top_scatteringangle_ttbarframe,
// *hvisgen_c_rn_top_scatteringangle_ttbarframe,
// *hvisgen_c_nk_top_scatteringangle_ttbarframe,
// *hvisgen_c_kn_top_scatteringangle_ttbarframe,
// *hvisgen_c_rj_top_scatteringangle_ttbarframe,
// *hvisgen_c_jr_top_scatteringangle_ttbarframe; TH1
// *hvisgen_c_han_top_scatteringangle_ttbarframe,
// *hvisgen_c_sca_top_scatteringangle_ttbarframe,
// *hvisgen_c_tra_top_scatteringangle_ttbarframe,
// *hvisgen_c_kjL_top_scatteringangle_ttbarframe,
// *hvisgen_c_rqL_top_scatteringangle_ttbarframe; TH1
// *hvisgen_c_rkP_top_scatteringangle_ttbarframe,
// *hvisgen_c_rkM_top_scatteringangle_ttbarframe,
// *hvisgen_c_nrP_top_scatteringangle_ttbarframe,
// *hvisgen_c_nrM_top_scatteringangle_ttbarframe,
// *hvisgen_c_nkP_top_scatteringangle_ttbarframe,
// *hvisgen_c_nkM_top_scatteringangle_ttbarframe; TH1
// *hvisgen_ll_cHel_top_scatteringangle_ttbarframe,
// *hvisgen_ll_cLab_top_scatteringangle_ttbarframe,
// *hvisgen_ll_kNorm_top_scatteringangle_ttbarframe,
// *hvisgen_ll_rNorm_top_scatteringangle_ttbarframe,
// *hvisgen_llbar_delta_phi_top_scatteringangle_ttbarframe,
// *hvisgen_llbar_delta_eta_top_scatteringangle_ttbarframe;

// **********************
// **********************
// RECO VS GEN HISTOGRAMS
// **********************
// **********************

// 1D reco vs. gen histograms
TH2D *hrecoVsgen_top_pt, *hrecoVsgen_l_pt, *hrecoVsgen_lbar_pt;
TH2D *hrecoVsgen_ttbar_pt, *hrecoVsgen_ttbar_mass, *hrecoVsgen_ttbar_delta_phi, *hrecoVsgen_ttbar_delta_eta,
    *hrecoVsgen_ttbar_rapidity, *hrecoVsgen_llbar_pt, *hrecoVsgen_llbar_mass;

TH2D *hrecoVsgen_b1k, *hrecoVsgen_b2k, *hrecoVsgen_b1j, *hrecoVsgen_b2j, *hrecoVsgen_b1r, *hrecoVsgen_b2r,
    *hrecoVsgen_b1q, *hrecoVsgen_b2q, *hrecoVsgen_b1n, *hrecoVsgen_b2n;
TH2D *hrecoVsgen_c_kk, *hrecoVsgen_c_rr, *hrecoVsgen_c_nn, *hrecoVsgen_c_kj, *hrecoVsgen_c_rq;
TH2D *hrecoVsgen_c_Prk, *hrecoVsgen_c_Mrk, *hrecoVsgen_c_Pnr, *hrecoVsgen_c_Mnr, *hrecoVsgen_c_Pnk, *hrecoVsgen_c_Mnk,
    *hrecoVsgen_c_Prj, *hrecoVsgen_c_Mrj;
TH2D *hrecoVsgen_c_rk, *hrecoVsgen_c_kr, *hrecoVsgen_c_nr, *hrecoVsgen_c_rn, *hrecoVsgen_c_nk, *hrecoVsgen_c_kn,
    *hrecoVsgen_c_rj, *hrecoVsgen_c_jr;
TH2D *hrecoVsgen_c_han, *hrecoVsgen_c_sca, *hrecoVsgen_c_tra, *hrecoVsgen_c_kjL, *hrecoVsgen_c_rqL;
TH2D *hrecoVsgen_c_rkP, *hrecoVsgen_c_rkM, *hrecoVsgen_c_nrP, *hrecoVsgen_c_nrM, *hrecoVsgen_c_nkP, *hrecoVsgen_c_nkM;
TH2D *hrecoVsgen_ll_cHel, *hrecoVsgen_ll_cLab, *hrecoVsgen_ll_kNorm, *hrecoVsgen_ll_rNorm, *hrecoVsgen_llbar_delta_phi,
    *hrecoVsgen_llbar_delta_eta;

// 2D reco vs. gen histograms
// Amandeep : Adding migration matrices

// ******
// mttbar
// ******
TH2D *hrecoVsgen_b1k_exjet, *hrecoVsgen_b2k_exjet, *hrecoVsgen_b1j_exjet, *hrecoVsgen_b2j_exjet, *hrecoVsgen_b1r_exjet,
    *hrecoVsgen_b2r_exjet, *hrecoVsgen_b1q_exjet, *hrecoVsgen_b2q_exjet, *hrecoVsgen_b1n_exjet, *hrecoVsgen_b2n_exjet;
TH2D *hrecoVsgen_c_kk_exjet, *hrecoVsgen_c_rr_exjet, *hrecoVsgen_c_nn_exjet, *hrecoVsgen_c_kj_exjet,
    *hrecoVsgen_c_rq_exjet;
TH2D *hrecoVsgen_c_Prk_exjet, *hrecoVsgen_c_Mrk_exjet, *hrecoVsgen_c_Pnr_exjet, *hrecoVsgen_c_Mnr_exjet,
    *hrecoVsgen_c_Pnk_exjet, *hrecoVsgen_c_Mnk_exjet, *hrecoVsgen_c_Prj_exjet, *hrecoVsgen_c_Mrj_exjet;
TH2D *hrecoVsgen_c_rk_exjet, *hrecoVsgen_c_kr_exjet, *hrecoVsgen_c_nr_exjet, *hrecoVsgen_c_rn_exjet,
    *hrecoVsgen_c_nk_exjet, *hrecoVsgen_c_kn_exjet, *hrecoVsgen_c_rj_exjet, *hrecoVsgen_c_jr_exjet;
TH2D *hrecoVsgen_c_han_exjet, *hrecoVsgen_c_sca_exjet, *hrecoVsgen_c_tra_exjet, *hrecoVsgen_c_kjL_exjet,
    *hrecoVsgen_c_rqL_exjet;
TH2D *hrecoVsgen_c_rkP_exjet, *hrecoVsgen_c_rkM_exjet, *hrecoVsgen_c_nrP_exjet, *hrecoVsgen_c_nrM_exjet,
    *hrecoVsgen_c_nkP_exjet, *hrecoVsgen_c_nkM_exjet;
TH2D *hrecoVsgen_ll_cHel_exjet, *hrecoVsgen_ll_cLab_exjet, *hrecoVsgen_ll_kNorm_exjet, *hrecoVsgen_ll_rNorm_exjet,
    *hrecoVsgen_llbar_delta_phi_exjet, *hrecoVsgen_llbar_delta_eta_exjet;

// ******
// top_pt
// ******
// TH2D *hrecoVsgen_b1k_top_pt  , *hrecoVsgen_b2k_top_pt,
// *hrecoVsgen_b1j_top_pt, *hrecoVsgen_b2j_top_pt, *hrecoVsgen_b1r_top_pt,
// *hrecoVsgen_b2r_top_pt, *hrecoVsgen_b1q_top_pt, *hrecoVsgen_b2q_top_pt,
// *hrecoVsgen_b1n_top_pt, *hrecoVsgen_b2n_top_pt; TH2D *hrecoVsgen_c_kk_top_pt
// , *hrecoVsgen_c_rr_top_pt, *hrecoVsgen_c_nn_top_pt, *hrecoVsgen_c_kj_top_pt,
// *hrecoVsgen_c_rq_top_pt; TH2D *hrecoVsgen_c_Prk_top_pt,
// *hrecoVsgen_c_Mrk_top_pt, *hrecoVsgen_c_Pnr_top_pt, *hrecoVsgen_c_Mnr_top_pt,
// *hrecoVsgen_c_Pnk_top_pt, *hrecoVsgen_c_Mnk_top_pt, *hrecoVsgen_c_Prj_top_pt,
// *hrecoVsgen_c_Mrj_top_pt; TH2D *hrecoVsgen_c_rk_top_pt ,
// *hrecoVsgen_c_kr_top_pt, *hrecoVsgen_c_nr_top_pt, *hrecoVsgen_c_rn_top_pt,
// *hrecoVsgen_c_nk_top_pt, *hrecoVsgen_c_kn_top_pt, *hrecoVsgen_c_rj_top_pt,
// *hrecoVsgen_c_jr_top_pt; TH2D *hrecoVsgen_c_han_top_pt,
// *hrecoVsgen_c_sca_top_pt, *hrecoVsgen_c_tra_top_pt, *hrecoVsgen_c_kjL_top_pt,
// *hrecoVsgen_c_rqL_top_pt; TH2D *hrecoVsgen_c_rkP_top_pt,
// *hrecoVsgen_c_rkM_top_pt, *hrecoVsgen_c_nrP_top_pt, *hrecoVsgen_c_nrM_top_pt,
// *hrecoVsgen_c_nkP_top_pt, *hrecoVsgen_c_nkM_top_pt; TH2D
// *hrecoVsgen_ll_cHel_top_pt, *hrecoVsgen_ll_cLab_top_pt,
// *hrecoVsgen_ll_kNorm_top_pt, *hrecoVsgen_ll_rNorm_top_pt,
// *hrecoVsgen_llbar_delta_phi_top_pt, *hrecoVsgen_llbar_delta_eta_top_pt;

// // ********************
// // top_scatteringangle_ttbarframe
// // ********************
// TH2D *hrecoVsgen_b1k_top_scatteringangle_ttbarframe  ,
// *hrecoVsgen_b2k_top_scatteringangle_ttbarframe,
// *hrecoVsgen_b1j_top_scatteringangle_ttbarframe,
// *hrecoVsgen_b2j_top_scatteringangle_ttbarframe,
// *hrecoVsgen_b1r_top_scatteringangle_ttbarframe,
// *hrecoVsgen_b2r_top_scatteringangle_ttbarframe,
// *hrecoVsgen_b1q_top_scatteringangle_ttbarframe,
// *hrecoVsgen_b2q_top_scatteringangle_ttbarframe,
// *hrecoVsgen_b1n_top_scatteringangle_ttbarframe,
// *hrecoVsgen_b2n_top_scatteringangle_ttbarframe; TH2D
// *hrecoVsgen_c_kk_top_scatteringangle_ttbarframe ,
// *hrecoVsgen_c_rr_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_nn_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_kj_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_rq_top_scatteringangle_ttbarframe; TH2D
// *hrecoVsgen_c_Prk_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_Mrk_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_Pnr_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_Mnr_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_Pnk_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_Mnk_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_Prj_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_Mrj_top_scatteringangle_ttbarframe; TH2D
// *hrecoVsgen_c_rk_top_scatteringangle_ttbarframe ,
// *hrecoVsgen_c_kr_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_nr_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_rn_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_nk_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_kn_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_rj_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_jr_top_scatteringangle_ttbarframe; TH2D
// *hrecoVsgen_c_han_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_sca_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_tra_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_kjL_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_rqL_top_scatteringangle_ttbarframe; TH2D
// *hrecoVsgen_c_rkP_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_rkM_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_nrP_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_nrM_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_nkP_top_scatteringangle_ttbarframe,
// *hrecoVsgen_c_nkM_top_scatteringangle_ttbarframe; TH2D
// *hrecoVsgen_ll_cHel_top_scatteringangle_ttbarframe,
// *hrecoVsgen_ll_cLab_top_scatteringangle_ttbarframe,
// *hrecoVsgen_ll_kNorm_top_scatteringangle_ttbarframe,
// *hrecoVsgen_ll_rNorm_top_scatteringangle_ttbarframe,
// *hrecoVsgen_llbar_delta_phi_top_scatteringangle_ttbarframe,
// *hrecoVsgen_llbar_delta_eta_top_scatteringangle_ttbarframe;

// ***********
// ***********
// RESOLUTIONS
// ***********
// ***********

// 1D resolution plots
TH2D *hresolutionbins_top_pt, *hresolutionbins_l_pt, *hresolutionbins_lbar_pt;
TH2D *hresolutionbins_ttbar_pt, *hresolutionbins_ttbar_mass, *hresolutionbins_ttbar_delta_phi,
    *hresolutionbins_ttbar_delta_eta, *hresolutionbins_ttbar_rapidity, *hresolutionbins_llbar_pt,
    *hresolutionbins_llbar_mass;
TH2D *hresolutionbins_b1k, *hresolutionbins_b2k, *hresolutionbins_b1j, *hresolutionbins_b2j, *hresolutionbins_b1r,
    *hresolutionbins_b2r, *hresolutionbins_b1q, *hresolutionbins_b2q, *hresolutionbins_b1n, *hresolutionbins_b2n;
TH2D *hresolutionbins_c_kk, *hresolutionbins_c_rr, *hresolutionbins_c_nn, *hresolutionbins_c_kj, *hresolutionbins_c_rq;
TH2D *hresolutionbins_c_Prk, *hresolutionbins_c_Mrk, *hresolutionbins_c_Pnr, *hresolutionbins_c_Mnr,
    *hresolutionbins_c_Pnk, *hresolutionbins_c_Mnk, *hresolutionbins_c_Prj, *hresolutionbins_c_Mrj;
TH2D *hresolutionbins_c_rk, *hresolutionbins_c_kr, *hresolutionbins_c_nr, *hresolutionbins_c_rn, *hresolutionbins_c_nk,
    *hresolutionbins_c_kn, *hresolutionbins_c_rj, *hresolutionbins_c_jr;
TH2D *hresolutionbins_c_han, *hresolutionbins_c_sca, *hresolutionbins_c_tra, *hresolutionbins_c_kjL,
    *hresolutionbins_c_rqL;
TH2D *hresolutionbins_c_rkP, *hresolutionbins_c_rkM, *hresolutionbins_c_nrP, *hresolutionbins_c_nrM,
    *hresolutionbins_c_nkP, *hresolutionbins_c_nkM;
TH2D *hresolutionbins_ll_cHel, *hresolutionbins_ll_cLab, *hresolutionbins_ll_kNorm, *hresolutionbins_ll_rNorm,
    *hresolutionbins_llbar_delta_phi, *hresolutionbins_llbar_delta_eta;

// 2D resolution plots
// Amandeep : Adding resolution plots

// ******
// mttbar
// ******
TH2D *hresolutionbins_b1k_exjet, *hresolutionbins_b2k_exjet, *hresolutionbins_b1j_exjet, *hresolutionbins_b2j_exjet,
    *hresolutionbins_b1r_exjet, *hresolutionbins_b2r_exjet, *hresolutionbins_b1q_exjet, *hresolutionbins_b2q_exjet,
    *hresolutionbins_b1n_exjet, *hresolutionbins_b2n_exjet;
TH2D *hresolutionbins_c_kk_exjet, *hresolutionbins_c_rr_exjet, *hresolutionbins_c_nn_exjet, *hresolutionbins_c_kj_exjet,
    *hresolutionbins_c_rq_exjet;
TH2D *hresolutionbins_c_Prk_exjet, *hresolutionbins_c_Mrk_exjet, *hresolutionbins_c_Pnr_exjet,
    *hresolutionbins_c_Mnr_exjet, *hresolutionbins_c_Pnk_exjet, *hresolutionbins_c_Mnk_exjet,
    *hresolutionbins_c_Prj_exjet, *hresolutionbins_c_Mrj_exjet;
TH2D *hresolutionbins_c_rk_exjet, *hresolutionbins_c_kr_exjet, *hresolutionbins_c_nr_exjet, *hresolutionbins_c_rn_exjet,
    *hresolutionbins_c_nk_exjet, *hresolutionbins_c_kn_exjet, *hresolutionbins_c_rj_exjet, *hresolutionbins_c_jr_exjet;
TH2D *hresolutionbins_c_han_exjet, *hresolutionbins_c_sca_exjet, *hresolutionbins_c_tra_exjet,
    *hresolutionbins_c_kjL_exjet, *hresolutionbins_c_rqL_exjet;
TH2D *hresolutionbins_c_rkP_exjet, *hresolutionbins_c_rkM_exjet, *hresolutionbins_c_nrP_exjet,
    *hresolutionbins_c_nrM_exjet, *hresolutionbins_c_nkP_exjet, *hresolutionbins_c_nkM_exjet;
TH2D *hresolutionbins_ll_cHel_exjet, *hresolutionbins_ll_cLab_exjet, *hresolutionbins_ll_kNorm_exjet,
    *hresolutionbins_ll_rNorm_exjet, *hresolutionbins_llbar_delta_phi_exjet, *hresolutionbins_llbar_delta_eta_exjet;

// ******
// top_pt
// ******
// TH2D *hresolutionbins_b1k_top_pt  , *hresolutionbins_b2k_top_pt,
// *hresolutionbins_b1j_top_pt, *hresolutionbins_b2j_top_pt,
// *hresolutionbins_b1r_top_pt, *hresolutionbins_b2r_top_pt,
// *hresolutionbins_b1q_top_pt, *hresolutionbins_b2q_top_pt,
// *hresolutionbins_b1n_top_pt, *hresolutionbins_b2n_top_pt; TH2D
// *hresolutionbins_c_kk_top_pt , *hresolutionbins_c_rr_top_pt,
// *hresolutionbins_c_nn_top_pt, *hresolutionbins_c_kj_top_pt,
// *hresolutionbins_c_rq_top_pt; TH2D *hresolutionbins_c_Prk_top_pt,
// *hresolutionbins_c_Mrk_top_pt, *hresolutionbins_c_Pnr_top_pt,
// *hresolutionbins_c_Mnr_top_pt, *hresolutionbins_c_Pnk_top_pt,
// *hresolutionbins_c_Mnk_top_pt, *hresolutionbins_c_Prj_top_pt,
// *hresolutionbins_c_Mrj_top_pt; TH2D *hresolutionbins_c_rk_top_pt ,
// *hresolutionbins_c_kr_top_pt, *hresolutionbins_c_nr_top_pt,
// *hresolutionbins_c_rn_top_pt, *hresolutionbins_c_nk_top_pt,
// *hresolutionbins_c_kn_top_pt, *hresolutionbins_c_rj_top_pt,
// *hresolutionbins_c_jr_top_pt; TH2D *hresolutionbins_c_han_top_pt,
// *hresolutionbins_c_sca_top_pt, *hresolutionbins_c_tra_top_pt,
// *hresolutionbins_c_kjL_top_pt, *hresolutionbins_c_rqL_top_pt; TH2D
// *hresolutionbins_c_rkP_top_pt, *hresolutionbins_c_rkM_top_pt,
// *hresolutionbins_c_nrP_top_pt, *hresolutionbins_c_nrM_top_pt,
// *hresolutionbins_c_nkP_top_pt, *hresolutionbins_c_nkM_top_pt; TH2D
// *hresolutionbins_ll_cHel_top_pt, *hresolutionbins_ll_cLab_top_pt,
// *hresolutionbins_ll_kNorm_top_pt, *hresolutionbins_ll_rNorm_top_pt,
// *hresolutionbins_llbar_delta_phi_top_pt,
// *hresolutionbins_llbar_delta_eta_top_pt;

// // ********************
// // top_scatteringangle_ttbarframe
// // ********************
// TH2D *hresolutionbins_b1k_top_scatteringangle_ttbarframe  ,
// *hresolutionbins_b2k_top_scatteringangle_ttbarframe,
// *hresolutionbins_b1j_top_scatteringangle_ttbarframe,
// *hresolutionbins_b2j_top_scatteringangle_ttbarframe,
// *hresolutionbins_b1r_top_scatteringangle_ttbarframe,
// *hresolutionbins_b2r_top_scatteringangle_ttbarframe,
// *hresolutionbins_b1q_top_scatteringangle_ttbarframe,
// *hresolutionbins_b2q_top_scatteringangle_ttbarframe,
// *hresolutionbins_b1n_top_scatteringangle_ttbarframe,
// *hresolutionbins_b2n_top_scatteringangle_ttbarframe; TH2D
// *hresolutionbins_c_kk_top_scatteringangle_ttbarframe ,
// *hresolutionbins_c_rr_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_nn_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_kj_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_rq_top_scatteringangle_ttbarframe; TH2D
// *hresolutionbins_c_Prk_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_Mrk_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_Pnr_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_Mnr_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_Pnk_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_Mnk_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_Prj_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_Mrj_top_scatteringangle_ttbarframe; TH2D
// *hresolutionbins_c_rk_top_scatteringangle_ttbarframe ,
// *hresolutionbins_c_kr_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_nr_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_rn_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_nk_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_kn_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_rj_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_jr_top_scatteringangle_ttbarframe; TH2D
// *hresolutionbins_c_han_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_sca_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_tra_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_kjL_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_rqL_top_scatteringangle_ttbarframe; TH2D
// *hresolutionbins_c_rkP_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_rkM_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_nrP_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_nrM_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_nkP_top_scatteringangle_ttbarframe,
// *hresolutionbins_c_nkM_top_scatteringangle_ttbarframe; TH2D
// *hresolutionbins_ll_cHel_top_scatteringangle_ttbarframe,
// *hresolutionbins_ll_cLab_top_scatteringangle_ttbarframe,
// *hresolutionbins_ll_kNorm_top_scatteringangle_ttbarframe,
// *hresolutionbins_ll_rNorm_top_scatteringangle_ttbarframe,
// *hresolutionbins_llbar_delta_phi_top_scatteringangle_ttbarframe,
// *hresolutionbins_llbar_delta_eta_top_scatteringangle_ttbarframe;

}  // namespace unfoldhistos

using namespace unfoldhistos;
// Reset the values of the variables

void makeTUnfoldHisto::Reset() {
    // std::cout<<"Resetting the variables"<<std::endl;
}

void makeTUnfoldHisto::Loop(TH1D *wtdEvts, std::string channel, bool symmetrizeSpinVar) {
    if (fChain == 0 || fChain0 == 0) return;
    std::cout << "symmetrizeSpinVar: " << symmetrizeSpinVar << std::endl;
    Long64_t nentries = fChain->GetEntriesFast();
    std::cout << "nentries: " << nentries << std::endl;

    std::vector<Long64_t> recoentry;
    std::vector<double> recogenlphi;
    float trueevts(0.0);

    int prescale;
    TRandom3 *random3 = new TRandom3();
    random3->SetSeed(101);

    TString filename = fout->GetName();

    // ********************
    // ********************
    // Fill reco histograms
    // ********************
    // ********************

    // Abraham: Changed for testing. Need to revert back.
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        // for (Long64_t jentry=0; jentry<200;jentry++) {

        if (jentry % 100000 == 0) std::cout << "Processing event: " << jentry << std::endl;

        if (LoadTree(jentry) < 0) break;

        // Reset
        makeTUnfoldHisto::Reset();
        fChain->GetEntry(jentry);

        trueevts = trueevts + trueLevelWeight;

        // Process - do operations, and fill histo
        recoentry.push_back(entry);
        recogenlphi.push_back(gen_l_phi);

        // Amandeep :
        // To avoid true underflow/overflow
        if (ttbar_mass >= 2000.0) {
            ttbar_mass = 1999.0;
        }

        if (ttbar_mass <= 250.0) {
            ttbar_mass = 251.0;
        }

        if (top_pt >= 550.0) {
            top_pt = 549.0;
        }

        if (n_extraJets_iso08 >= 3.5) {
            n_extraJets_iso08 = 3.0;
        }
        // End

        // **************************
        // Fill 1D reco distributions
        // **************************

        // Kinematic
        fillUnderOverFlow(hreco_top_pt, reco_top_pt_binning, std::vector<double>{top_pt}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_l_pt, reco_l_pt_binning, std::vector<double>{l_pt}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_lbar_pt, reco_lbar_pt_binning, std::vector<double>{lbar_pt}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_ttbar_pt, reco_ttbar_pt_binning, std::vector<double>{ttbar_pt}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_ttbar_mass, reco_ttbar_mass_binning, std::vector<double>{ttbar_mass}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_ttbar_delta_phi, reco_ttbar_delta_phi_binning, std::vector<double>{ttbar_delta_phi},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_ttbar_delta_eta, reco_ttbar_delta_eta_binning, std::vector<double>{ttbar_delta_eta},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_ttbar_rapidity, reco_ttbar_rapidity_binning, std::vector<double>{ttbar_rapidity},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_llbar_pt, reco_llbar_pt_binning, std::vector<double>{llbar_pt}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_llbar_mass, reco_llbar_mass_binning, std::vector<double>{llbar_mass}, eventWeight,
                          symmetrizeSpinVar);

        // Spin corr

        // Polarizations
        fillUnderOverFlow(hreco_b1k, reco_b1k_binning, std::vector<double>{b1k}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_b2k, reco_b2k_binning, std::vector<double>{b2k}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_b1j, reco_b1j_binning, std::vector<double>{b1j}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_b2j, reco_b2j_binning, std::vector<double>{b2j}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_b1r, reco_b1r_binning, std::vector<double>{b1r}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_b2r, reco_b2r_binning, std::vector<double>{b2r}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_b1q, reco_b1q_binning, std::vector<double>{b1q}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_b2q, reco_b2q_binning, std::vector<double>{b2q}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_b1n, reco_b1n_binning, std::vector<double>{b1n}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_b2n, reco_b2n_binning, std::vector<double>{b2n}, eventWeight, symmetrizeSpinVar);

        // Diagonal elements
        fillUnderOverFlow(hreco_c_kk, reco_c_kk_binning, std::vector<double>{c_kk}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_rr, reco_c_rr_binning, std::vector<double>{c_rr}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_nn, reco_c_nn_binning, std::vector<double>{c_nn}, eventWeight, symmetrizeSpinVar);

        // Off diagonal elements
        fillUnderOverFlow(hreco_c_rk, reco_c_rk_binning, std::vector<double>{c_rk}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_kr, reco_c_kr_binning, std::vector<double>{c_kr}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_nr, reco_c_nr_binning, std::vector<double>{c_nr}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_rn, reco_c_rn_binning, std::vector<double>{c_rn}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_nk, reco_c_nk_binning, std::vector<double>{c_nk}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_kn, reco_c_kn_binning, std::vector<double>{c_kn}, eventWeight, symmetrizeSpinVar);

        fillUnderOverFlow(hreco_c_Prk, reco_c_Prk_binning, std::vector<double>{c_rk + c_kr}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_Mrk, reco_c_Mrk_binning, std::vector<double>{c_rk - c_kr}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_Pnr, reco_c_Pnr_binning, std::vector<double>{c_nr + c_rn}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_Mnr, reco_c_Mnr_binning, std::vector<double>{c_nr - c_rn}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_Pnk, reco_c_Pnk_binning, std::vector<double>{c_nk + c_kn}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_Mnk, reco_c_Mnk_binning, std::vector<double>{c_nk - c_kn}, eventWeight,
                          symmetrizeSpinVar);

        // Starred axes
        fillUnderOverFlow(hreco_c_kj, reco_c_kj_binning, std::vector<double>{c_kj}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_rq, reco_c_rq_binning, std::vector<double>{c_rq}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_rj, reco_c_rj_binning, std::vector<double>{c_rj}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_jr, reco_c_jr_binning, std::vector<double>{c_jr}, eventWeight, symmetrizeSpinVar);

        fillUnderOverFlow(hreco_c_Prj, reco_c_Prj_binning, std::vector<double>{c_rj + c_jr}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_Mrj, reco_c_Mrj_binning, std::vector<double>{c_rj - c_jr}, eventWeight,
                          symmetrizeSpinVar);

        // New variables
        fillUnderOverFlow(hreco_c_han, reco_c_han_binning, std::vector<double>{+c_kk - c_rr - c_nn}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_sca, reco_c_sca_binning, std::vector<double>{-c_kk + c_rr - c_nn}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_tra, reco_c_tra_binning, std::vector<double>{-c_kk - c_rr + c_nn}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_kjL, reco_c_kjL_binning, std::vector<double>{-c_kj - c_rr - c_nn}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_rqL, reco_c_rqL_binning, std::vector<double>{-c_kk - c_rq - c_nn}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_rkP, reco_c_rkP_binning, std::vector<double>{-c_rk - c_kr - c_nn}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_rkM, reco_c_rkM_binning, std::vector<double>{-c_rk + c_kr - c_nn}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_nrP, reco_c_nrP_binning, std::vector<double>{-c_nr - c_rn - c_kk}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_nrM, reco_c_nrM_binning, std::vector<double>{-c_nr + c_rn - c_kk}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_nkP, reco_c_nkP_binning, std::vector<double>{-c_nk - c_kn - c_rr}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_nkM, reco_c_nkM_binning, std::vector<double>{-c_nk + c_kn - c_rr}, eventWeight,
                          symmetrizeSpinVar);
        // std:: cout << "gen_n1," << n_extraJets_iso08 << "," << c_kk<<
        // std::endl; Lab frame
        fillUnderOverFlow(hreco_ll_cHel, reco_ll_cHel_binning, std::vector<double>{ll_cHel}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_ll_cLab, reco_ll_cLab_binning, std::vector<double>{ll_cLab}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_ll_kNorm, reco_ll_kNorm_binning, std::vector<double>{ll_kNorm}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_ll_rNorm, reco_ll_rNorm_binning, std::vector<double>{ll_rNorm}, eventWeight,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hreco_llbar_delta_phi, reco_llbar_delta_phi_binning, std::vector<double>{llbar_delta_phi},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_llbar_delta_eta, reco_llbar_delta_eta_binning, std::vector<double>{llbar_delta_eta},
                          eventWeight, symmetrizeSpinVar);

        // ***********************
        // Fill 2D reco histograms
        // ***********************

        // Amandeep : Adding 2D reco histograms
        // fillUnderOverFlow(hreco_c_kk_mttbar, reco_c_kk_mttbar_binning,
        // std::vector<double>{c_kk, ttbar_mass}, eventWeight,
        // symmetrizeSpinVar);

        // ******
        // mttbar
        // ******

        // Spin corr

        // Polarizations
        // std:: cout << "gen_n2," << n_extraJets_iso08 << std::endl;
        fillUnderOverFlow(hreco_b1k_exjet, reco_b1k_exjet_binning, std::vector<double>{b1k, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_b2k_exjet, reco_b2k_exjet_binning, std::vector<double>{b2k, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_b1j_exjet, reco_b1j_exjet_binning, std::vector<double>{b1j, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_b2j_exjet, reco_b2j_exjet_binning, std::vector<double>{b2j, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_b1r_exjet, reco_b1r_exjet_binning, std::vector<double>{b1r, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_b2r_exjet, reco_b2r_exjet_binning, std::vector<double>{b2r, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_b1q_exjet, reco_b1q_exjet_binning, std::vector<double>{b1q, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_b2q_exjet, reco_b2q_exjet_binning, std::vector<double>{b2q, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_b1n_exjet, reco_b1n_exjet_binning, std::vector<double>{b1n, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_b2n_exjet, reco_b2n_exjet_binning, std::vector<double>{b2n, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);

        // Diagonal elements
        fillUnderOverFlow(hreco_c_kk_exjet, reco_c_kk_exjet_binning, std::vector<double>{c_kk, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_rr_exjet, reco_c_rr_exjet_binning, std::vector<double>{c_rr, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_nn_exjet, reco_c_nn_exjet_binning, std::vector<double>{c_nn, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);

        // Off diagonal elements
        fillUnderOverFlow(hreco_c_rk_exjet, reco_c_rk_exjet_binning, std::vector<double>{c_rk, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_kr_exjet, reco_c_kr_exjet_binning, std::vector<double>{c_kr, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_nr_exjet, reco_c_nr_exjet_binning, std::vector<double>{c_nr, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_rn_exjet, reco_c_rn_exjet_binning, std::vector<double>{c_rn, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_nk_exjet, reco_c_nk_exjet_binning, std::vector<double>{c_nk, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_kn_exjet, reco_c_kn_exjet_binning, std::vector<double>{c_kn, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);

        fillUnderOverFlow(hreco_c_Prk_exjet, reco_c_Prk_exjet_binning,
                          std::vector<double>{c_rk + c_kr, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_Mrk_exjet, reco_c_Mrk_exjet_binning,
                          std::vector<double>{c_rk - c_kr, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_Pnr_exjet, reco_c_Pnr_exjet_binning,
                          std::vector<double>{c_nr + c_rn, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_Mnr_exjet, reco_c_Mnr_exjet_binning,
                          std::vector<double>{c_nr - c_rn, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_Pnk_exjet, reco_c_Pnk_exjet_binning,
                          std::vector<double>{c_nk + c_kn, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_Mnk_exjet, reco_c_Mnk_exjet_binning,
                          std::vector<double>{c_nk - c_kn, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);

        // Starred axes
        fillUnderOverFlow(hreco_c_kj_exjet, reco_c_kj_exjet_binning, std::vector<double>{c_kj, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_rq_exjet, reco_c_rq_exjet_binning, std::vector<double>{c_rq, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_rj_exjet, reco_c_rj_exjet_binning, std::vector<double>{c_rj, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_jr_exjet, reco_c_jr_exjet_binning, std::vector<double>{c_jr, n_extraJets_iso08},
                          eventWeight, symmetrizeSpinVar);

        fillUnderOverFlow(hreco_c_Prj_exjet, reco_c_Prj_exjet_binning,
                          std::vector<double>{c_rj + c_jr, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_Mrj_exjet, reco_c_Mrj_exjet_binning,
                          std::vector<double>{c_rj - c_jr, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);

        // New variables
        fillUnderOverFlow(hreco_c_han_exjet, reco_c_han_exjet_binning,
                          std::vector<double>{+c_kk - c_rr - c_nn, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_sca_exjet, reco_c_sca_exjet_binning,
                          std::vector<double>{-c_kk + c_rr - c_nn, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_tra_exjet, reco_c_tra_exjet_binning,
                          std::vector<double>{-c_kk - c_rr + c_nn, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_kjL_exjet, reco_c_kjL_exjet_binning,
                          std::vector<double>{-c_kj - c_rr - c_nn, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_rqL_exjet, reco_c_rqL_exjet_binning,
                          std::vector<double>{-c_kk - c_rq - c_nn, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_rkP_exjet, reco_c_rkP_exjet_binning,
                          std::vector<double>{-c_rk - c_kr - c_nn, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_rkM_exjet, reco_c_rkM_exjet_binning,
                          std::vector<double>{-c_rk + c_kr - c_nn, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_nrP_exjet, reco_c_nrP_exjet_binning,
                          std::vector<double>{-c_nr - c_rn - c_kk, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_nrM_exjet, reco_c_nrM_exjet_binning,
                          std::vector<double>{-c_nr + c_rn - c_kk, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_nkP_exjet, reco_c_nkP_exjet_binning,
                          std::vector<double>{-c_nk - c_kn - c_rr, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_c_nkM_exjet, reco_c_nkM_exjet_binning,
                          std::vector<double>{-c_nk + c_kn - c_rr, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);

        // Lab frame variables
        fillUnderOverFlow(hreco_ll_cHel_exjet, reco_ll_cHel_exjet_binning,
                          std::vector<double>{ll_cHel, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_ll_cLab_exjet, reco_ll_cLab_exjet_binning,
                          std::vector<double>{ll_cLab, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_ll_kNorm_exjet, reco_ll_kNorm_exjet_binning,
                          std::vector<double>{ll_kNorm, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_ll_rNorm_exjet, reco_ll_rNorm_exjet_binning,
                          std::vector<double>{ll_rNorm, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_llbar_delta_phi_exjet, reco_llbar_delta_phi_exjet_binning,
                          std::vector<double>{llbar_delta_phi, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);
        fillUnderOverFlow(hreco_llbar_delta_eta_exjet, reco_llbar_delta_eta_exjet_binning,
                          std::vector<double>{llbar_delta_eta, n_extraJets_iso08}, eventWeight, symmetrizeSpinVar);

        // ******
        // top_pt
        // ******

        // Spin corr

        // Polarizations
        // fillUnderOverFlow(hreco_b1k_top_pt, reco_b1k_top_pt_binning,
        // std::vector<double>{b1k, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_b2k_top_pt, reco_b2k_top_pt_binning,
        // std::vector<double>{b2k, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_b1j_top_pt, reco_b1j_top_pt_binning,
        // std::vector<double>{b1j, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_b2j_top_pt, reco_b2j_top_pt_binning,
        // std::vector<double>{b2j, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_b1r_top_pt, reco_b1r_top_pt_binning,
        // std::vector<double>{b1r, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_b2r_top_pt, reco_b2r_top_pt_binning,
        // std::vector<double>{b2r, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_b1q_top_pt, reco_b1q_top_pt_binning,
        // std::vector<double>{b1q, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_b2q_top_pt, reco_b2q_top_pt_binning,
        // std::vector<double>{b2q, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_b1n_top_pt, reco_b1n_top_pt_binning,
        // std::vector<double>{b1n, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_b2n_top_pt, reco_b2n_top_pt_binning,
        // std::vector<double>{b2n, top_pt}, eventWeight, symmetrizeSpinVar);

        // // Diagonal elements
        // fillUnderOverFlow(hreco_c_kk_top_pt, reco_c_kk_top_pt_binning,
        // std::vector<double>{c_kk, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_rr_top_pt, reco_c_rr_top_pt_binning,
        // std::vector<double>{c_rr, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_nn_top_pt, reco_c_nn_top_pt_binning,
        // std::vector<double>{c_nn, top_pt}, eventWeight, symmetrizeSpinVar);

        // // Off diagonal elements
        // fillUnderOverFlow(hreco_c_rk_top_pt, reco_c_rk_top_pt_binning,
        // std::vector<double>{c_rk, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_kr_top_pt, reco_c_kr_top_pt_binning,
        // std::vector<double>{c_kr, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_nr_top_pt, reco_c_nr_top_pt_binning,
        // std::vector<double>{c_nr, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_rn_top_pt, reco_c_rn_top_pt_binning,
        // std::vector<double>{c_rn, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_nk_top_pt, reco_c_nk_top_pt_binning,
        // std::vector<double>{c_nk, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_kn_top_pt, reco_c_kn_top_pt_binning,
        // std::vector<double>{c_kn, top_pt}, eventWeight, symmetrizeSpinVar);

        // fillUnderOverFlow(hreco_c_Prk_top_pt, reco_c_Prk_top_pt_binning,
        // std::vector<double>{c_rk + c_kr, top_pt}, eventWeight,
        // symmetrizeSpinVar); fillUnderOverFlow(hreco_c_Mrk_top_pt,
        // reco_c_Mrk_top_pt_binning, std::vector<double>{c_rk - c_kr, top_pt},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_Pnr_top_pt, reco_c_Pnr_top_pt_binning,
        // std::vector<double>{c_nr + c_rn, top_pt}, eventWeight,
        // symmetrizeSpinVar); fillUnderOverFlow(hreco_c_Mnr_top_pt,
        // reco_c_Mnr_top_pt_binning, std::vector<double>{c_nr - c_rn, top_pt},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_Pnk_top_pt, reco_c_Pnk_top_pt_binning,
        // std::vector<double>{c_nk + c_kn, top_pt}, eventWeight,
        // symmetrizeSpinVar); fillUnderOverFlow(hreco_c_Mnk_top_pt,
        // reco_c_Mnk_top_pt_binning, std::vector<double>{c_nk - c_kn, top_pt},
        // eventWeight, symmetrizeSpinVar);

        // // Starred axes
        // fillUnderOverFlow(hreco_c_kj_top_pt, reco_c_kj_top_pt_binning,
        // std::vector<double>{c_kj, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_rq_top_pt, reco_c_rq_top_pt_binning,
        // std::vector<double>{c_rq, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_rj_top_pt, reco_c_rj_top_pt_binning,
        // std::vector<double>{c_rj, top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_jr_top_pt, reco_c_jr_top_pt_binning,
        // std::vector<double>{c_jr, top_pt}, eventWeight, symmetrizeSpinVar);

        // fillUnderOverFlow(hreco_c_Prj_top_pt, reco_c_Prj_top_pt_binning,
        // std::vector<double>{c_rj + c_jr, top_pt}, eventWeight,
        // symmetrizeSpinVar); fillUnderOverFlow(hreco_c_Mrj_top_pt,
        // reco_c_Mrj_top_pt_binning, std::vector<double>{c_rj - c_jr, top_pt},
        // eventWeight, symmetrizeSpinVar);

        // // New variables
        // fillUnderOverFlow(hreco_c_han_top_pt, reco_c_han_top_pt_binning,
        // std::vector<double>{+c_kk - c_rr - c_nn, top_pt}, eventWeight,
        // symmetrizeSpinVar); fillUnderOverFlow(hreco_c_sca_top_pt,
        // reco_c_sca_top_pt_binning, std::vector<double>{-c_kk + c_rr - c_nn,
        // top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_tra_top_pt, reco_c_tra_top_pt_binning,
        // std::vector<double>{-c_kk - c_rr + c_nn, top_pt}, eventWeight,
        // symmetrizeSpinVar); fillUnderOverFlow(hreco_c_kjL_top_pt,
        // reco_c_kjL_top_pt_binning, std::vector<double>{-c_kj - c_rr - c_nn,
        // top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_rqL_top_pt, reco_c_rqL_top_pt_binning,
        // std::vector<double>{-c_kk - c_rq - c_nn, top_pt}, eventWeight,
        // symmetrizeSpinVar); fillUnderOverFlow(hreco_c_rkP_top_pt,
        // reco_c_rkP_top_pt_binning, std::vector<double>{-c_rk - c_kr - c_nn,
        // top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_rkM_top_pt, reco_c_rkM_top_pt_binning,
        // std::vector<double>{-c_rk + c_kr - c_nn, top_pt}, eventWeight,
        // symmetrizeSpinVar); fillUnderOverFlow(hreco_c_nrP_top_pt,
        // reco_c_nrP_top_pt_binning, std::vector<double>{-c_nr - c_rn - c_kk,
        // top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_nrM_top_pt, reco_c_nrM_top_pt_binning,
        // std::vector<double>{-c_nr + c_rn - c_kk, top_pt}, eventWeight,
        // symmetrizeSpinVar); fillUnderOverFlow(hreco_c_nkP_top_pt,
        // reco_c_nkP_top_pt_binning, std::vector<double>{-c_nk - c_kn - c_rr,
        // top_pt}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_nkM_top_pt, reco_c_nkM_top_pt_binning,
        // std::vector<double>{-c_nk + c_kn - c_rr, top_pt}, eventWeight,
        // symmetrizeSpinVar);

        // // Lab frame variables
        // fillUnderOverFlow(hreco_ll_cHel_top_pt, reco_ll_cHel_top_pt_binning,
        // std::vector<double>{ll_cHel, top_pt}, eventWeight,
        // symmetrizeSpinVar); fillUnderOverFlow(hreco_ll_cLab_top_pt,
        // reco_ll_cLab_top_pt_binning, std::vector<double>{ll_cLab, top_pt},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_ll_kNorm_top_pt,
        // reco_ll_kNorm_top_pt_binning, std::vector<double>{ll_kNorm, top_pt},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_ll_rNorm_top_pt,
        // reco_ll_rNorm_top_pt_binning, std::vector<double>{ll_rNorm, top_pt},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_llbar_delta_phi_top_pt,
        // reco_llbar_delta_phi_top_pt_binning,
        // std::vector<double>{llbar_delta_phi, top_pt}, eventWeight,
        // symmetrizeSpinVar); fillUnderOverFlow(hreco_llbar_delta_eta_top_pt,
        // reco_llbar_delta_eta_top_pt_binning,
        // std::vector<double>{llbar_delta_eta, top_pt}, eventWeight,
        // symmetrizeSpinVar);

        // ******************************
        // top_scatteringangle_ttbarframe
        // ******************************

        // Spin corr

        // Polarizations
        // fillUnderOverFlow(hreco_b1k_top_scatteringangle_ttbarframe,
        // reco_b1k_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{b1k, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_b2k_top_scatteringangle_ttbarframe,
        // reco_b2k_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{b2k, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_b1j_top_scatteringangle_ttbarframe,
        // reco_b1j_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{b1j, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_b2j_top_scatteringangle_ttbarframe,
        // reco_b2j_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{b2j, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_b1r_top_scatteringangle_ttbarframe,
        // reco_b1r_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{b1r, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_b2r_top_scatteringangle_ttbarframe,
        // reco_b2r_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{b2r, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_b1q_top_scatteringangle_ttbarframe,
        // reco_b1q_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{b1q, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_b2q_top_scatteringangle_ttbarframe,
        // reco_b2q_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{b2q, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_b1n_top_scatteringangle_ttbarframe,
        // reco_b1n_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{b1n, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_b2n_top_scatteringangle_ttbarframe,
        // reco_b2n_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{b2n, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);

        // // Diagonal elements
        // fillUnderOverFlow(hreco_c_kk_top_scatteringangle_ttbarframe,
        // reco_c_kk_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_kk, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_rr_top_scatteringangle_ttbarframe,
        // reco_c_rr_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_rr, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_nn_top_scatteringangle_ttbarframe,
        // reco_c_nn_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_nn, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);

        // Off diagonal elements
        // fillUnderOverFlow(hreco_c_rk_top_scatteringangle_ttbarframe,
        // reco_c_rk_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_rk, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_kr_top_scatteringangle_ttbarframe,
        // reco_c_kr_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_kr, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_nr_top_scatteringangle_ttbarframe,
        // reco_c_nr_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_nr, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_rn_top_scatteringangle_ttbarframe,
        // reco_c_rn_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_rn, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_nk_top_scatteringangle_ttbarframe,
        // reco_c_nk_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_nk, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_kn_top_scatteringangle_ttbarframe,
        // reco_c_kn_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_kn, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);

        // fillUnderOverFlow(hreco_c_Prk_top_scatteringangle_ttbarframe,
        // reco_c_Prk_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_rk + c_kr, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_Mrk_top_scatteringangle_ttbarframe,
        // reco_c_Mrk_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_rk - c_kr, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_Pnr_top_scatteringangle_ttbarframe,
        // reco_c_Pnr_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_nr + c_rn, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_Mnr_top_scatteringangle_ttbarframe,
        // reco_c_Mnr_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_nr - c_rn, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_Pnk_top_scatteringangle_ttbarframe,
        // reco_c_Pnk_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_nk + c_kn, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_Mnk_top_scatteringangle_ttbarframe,
        // reco_c_Mnk_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_nk - c_kn, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);

        // // Starred axes
        // fillUnderOverFlow(hreco_c_kj_top_scatteringangle_ttbarframe,
        // reco_c_kj_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_kj, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_rq_top_scatteringangle_ttbarframe,
        // reco_c_rq_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_rq, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_rj_top_scatteringangle_ttbarframe,
        // reco_c_rj_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_rj, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_jr_top_scatteringangle_ttbarframe,
        // reco_c_jr_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_jr, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);

        // fillUnderOverFlow(hreco_c_Prj_top_scatteringangle_ttbarframe,
        // reco_c_Prj_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_rj + c_jr, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_Mrj_top_scatteringangle_ttbarframe,
        // reco_c_Mrj_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{c_rj - c_jr, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);

        // New variables
        // fillUnderOverFlow(hreco_c_han_top_scatteringangle_ttbarframe,
        // reco_c_han_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{+c_kk - c_rr - c_nn,
        // top_scatteringangle_ttbarframe}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_sca_top_scatteringangle_ttbarframe,
        // reco_c_sca_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-c_kk + c_rr - c_nn,
        // top_scatteringangle_ttbarframe}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_tra_top_scatteringangle_ttbarframe,
        // reco_c_tra_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-c_kk - c_rr + c_nn,
        // top_scatteringangle_ttbarframe}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_kjL_top_scatteringangle_ttbarframe,
        // reco_c_kjL_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-c_kj - c_rr - c_nn,
        // top_scatteringangle_ttbarframe}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_rqL_top_scatteringangle_ttbarframe,
        // reco_c_rqL_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-c_kk - c_rq - c_nn,
        // top_scatteringangle_ttbarframe}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_rkP_top_scatteringangle_ttbarframe,
        // reco_c_rkP_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-c_rk - c_kr - c_nn,
        // top_scatteringangle_ttbarframe}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_rkM_top_scatteringangle_ttbarframe,
        // reco_c_rkM_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-c_rk + c_kr - c_nn,
        // top_scatteringangle_ttbarframe}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_nrP_top_scatteringangle_ttbarframe,
        // reco_c_nrP_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-c_nr - c_rn - c_kk,
        // top_scatteringangle_ttbarframe}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_nrM_top_scatteringangle_ttbarframe,
        // reco_c_nrM_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-c_nr + c_rn - c_kk,
        // top_scatteringangle_ttbarframe}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_nkP_top_scatteringangle_ttbarframe,
        // reco_c_nkP_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-c_nk - c_kn - c_rr,
        // top_scatteringangle_ttbarframe}, eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_c_nkM_top_scatteringangle_ttbarframe,
        // reco_c_nkM_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-c_nk + c_kn - c_rr,
        // top_scatteringangle_ttbarframe}, eventWeight, symmetrizeSpinVar);

        // // Lab frame variables
        // fillUnderOverFlow(hreco_ll_cHel_top_scatteringangle_ttbarframe,
        // reco_ll_cHel_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{ll_cHel, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_ll_cLab_top_scatteringangle_ttbarframe,
        // reco_ll_cLab_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{ll_cLab, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_ll_kNorm_top_scatteringangle_ttbarframe,
        // reco_ll_kNorm_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{ll_kNorm, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_ll_rNorm_top_scatteringangle_ttbarframe,
        // reco_ll_rNorm_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{ll_rNorm, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_llbar_delta_phi_top_scatteringangle_ttbarframe,
        // reco_llbar_delta_phi_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{llbar_delta_phi, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);
        // fillUnderOverFlow(hreco_llbar_delta_eta_top_scatteringangle_ttbarframe,
        // reco_llbar_delta_eta_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{llbar_delta_eta, top_scatteringangle_ttbarframe},
        // eventWeight, symmetrizeSpinVar);

        // *************************
        // Fill Bootstrap histograms
        // *************************

        // Amandeep : Do we really need to fill psuedos for the MC samples?

        // Original :
        // if(doBootstrap && (filename.Contains("ttbarsignal") ||
        // filename.Contains("ttbarbg") || filename.Contains("run"))) { //only
        // for ttbar or data Modified :
        if (doBootstrap && (filename.Contains("run"))) {  // only for data
            if (jentry == 0) {
                int seed = nentries % 10000 + int(abs(l_phi * 10000));
                std::cout << "Filling " << nPE << " bootstrap samples using seed " << seed << ":" << std::endl;
                random3->SetSeed(seed);  // so the seed is different for different samples
                                         // (which may be added together)
            }
            if (jentry % 1000 == 0) std::cout << "Processing event: " << jentry << " of " << nentries << std::endl;
            for (int iPE = 0; iPE < nPE; ++iPE) {
                Int_t temp_event_multiplicity = random3->Poisson(1);
                if (temp_event_multiplicity > 0) {
                    for (int imult = 0; imult < temp_event_multiplicity; ++imult) {
                        // Kinematic
                        fillUnderOverFlow(hrecoBootstrap_top_pt[iPE], reco_top_pt_binning, std::vector<double>{top_pt},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_l_pt[iPE], reco_l_pt_binning, std::vector<double>{l_pt},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_lbar_pt[iPE], reco_lbar_pt_binning,
                                          std::vector<double>{lbar_pt}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_ttbar_pt[iPE], reco_ttbar_pt_binning,
                                          std::vector<double>{ttbar_pt}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_ttbar_mass[iPE], reco_ttbar_mass_binning,
                                          std::vector<double>{ttbar_mass}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_ttbar_delta_phi[iPE], reco_ttbar_delta_phi_binning,
                                          std::vector<double>{ttbar_delta_phi}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_ttbar_delta_eta[iPE], reco_ttbar_delta_eta_binning,
                                          std::vector<double>{ttbar_delta_eta}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_ttbar_rapidity[iPE], reco_ttbar_rapidity_binning,
                                          std::vector<double>{ttbar_rapidity}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_llbar_pt[iPE], reco_llbar_pt_binning,
                                          std::vector<double>{llbar_pt}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_llbar_mass[iPE], reco_llbar_mass_binning,
                                          std::vector<double>{llbar_mass}, eventWeight, symmetrizeSpinVar);

                        // Polarizations
                        fillUnderOverFlow(hrecoBootstrap_b1k[iPE], reco_b1k_binning, std::vector<double>{b1k},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_b2k[iPE], reco_b2k_binning, std::vector<double>{b2k},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_b1j[iPE], reco_b1j_binning, std::vector<double>{b1j},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_b2j[iPE], reco_b2j_binning, std::vector<double>{b2j},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_b1r[iPE], reco_b1r_binning, std::vector<double>{b1r},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_b2r[iPE], reco_b2r_binning, std::vector<double>{b2r},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_b1q[iPE], reco_b1q_binning, std::vector<double>{b1q},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_b2q[iPE], reco_b2q_binning, std::vector<double>{b2q},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_b1n[iPE], reco_b1n_binning, std::vector<double>{b1n},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_b2n[iPE], reco_b2n_binning, std::vector<double>{b2n},
                                          eventWeight, symmetrizeSpinVar);

                        // Diagonal elements
                        fillUnderOverFlow(hrecoBootstrap_c_kk[iPE], reco_c_kk_binning, std::vector<double>{c_kk},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_rr[iPE], reco_c_rr_binning, std::vector<double>{c_rr},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_nn[iPE], reco_c_nn_binning, std::vector<double>{c_nn},
                                          eventWeight, symmetrizeSpinVar);

                        // Off diagonal elements
                        fillUnderOverFlow(hrecoBootstrap_c_rk[iPE], reco_c_rk_binning, std::vector<double>{c_rk},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_kr[iPE], reco_c_kr_binning, std::vector<double>{c_kr},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_nr[iPE], reco_c_nr_binning, std::vector<double>{c_nr},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_rn[iPE], reco_c_rn_binning, std::vector<double>{c_rn},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_nk[iPE], reco_c_nk_binning, std::vector<double>{c_nk},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_kn[iPE], reco_c_kn_binning, std::vector<double>{c_kn},
                                          eventWeight, symmetrizeSpinVar);

                        fillUnderOverFlow(hrecoBootstrap_c_Prk[iPE], reco_c_Prk_binning,
                                          std::vector<double>{c_rk + c_kr}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_Mrk[iPE], reco_c_Mrk_binning,
                                          std::vector<double>{c_rk - c_kr}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_Pnr[iPE], reco_c_Pnr_binning,
                                          std::vector<double>{c_nr + c_rn}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_Mnr[iPE], reco_c_Mnr_binning,
                                          std::vector<double>{c_nr - c_rn}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_Pnk[iPE], reco_c_Pnk_binning,
                                          std::vector<double>{c_nk + c_kn}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_Mnk[iPE], reco_c_Mnk_binning,
                                          std::vector<double>{c_nk - c_kn}, eventWeight, symmetrizeSpinVar);

                        // Starred axes
                        fillUnderOverFlow(hrecoBootstrap_c_kj[iPE], reco_c_kj_binning, std::vector<double>{c_kj},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_rq[iPE], reco_c_rq_binning, std::vector<double>{c_rq},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_rj[iPE], reco_c_rj_binning, std::vector<double>{c_rj},
                                          eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_jr[iPE], reco_c_jr_binning, std::vector<double>{c_jr},
                                          eventWeight, symmetrizeSpinVar);

                        fillUnderOverFlow(hrecoBootstrap_c_Prj[iPE], reco_c_Prj_binning,
                                          std::vector<double>{c_rj + c_jr}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_Mrj[iPE], reco_c_Mrj_binning,
                                          std::vector<double>{c_rj - c_jr}, eventWeight, symmetrizeSpinVar);

                        // New variables
                        fillUnderOverFlow(hrecoBootstrap_c_han[iPE], reco_c_han_binning,
                                          std::vector<double>{+c_kk - c_rr - c_nn}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_sca[iPE], reco_c_sca_binning,
                                          std::vector<double>{-c_kk + c_rr - c_nn}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_tra[iPE], reco_c_tra_binning,
                                          std::vector<double>{-c_kk - c_rr + c_nn}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_kjL[iPE], reco_c_kjL_binning,
                                          std::vector<double>{-c_kj - c_rr - c_nn}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_rqL[iPE], reco_c_rqL_binning,
                                          std::vector<double>{-c_kk - c_rq - c_nn}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_rkP[iPE], reco_c_rkP_binning,
                                          std::vector<double>{-c_rk - c_kr - c_nn}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_rkM[iPE], reco_c_rkM_binning,
                                          std::vector<double>{-c_rk + c_kr - c_nn}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_nrP[iPE], reco_c_nrP_binning,
                                          std::vector<double>{-c_nr - c_rn - c_kk}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_nrM[iPE], reco_c_nrM_binning,
                                          std::vector<double>{-c_nr + c_rn - c_kk}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_nkP[iPE], reco_c_nkP_binning,
                                          std::vector<double>{-c_nk - c_kn - c_rr}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_c_nkM[iPE], reco_c_nkM_binning,
                                          std::vector<double>{-c_nk + c_kn - c_rr}, eventWeight, symmetrizeSpinVar);

                        // Lab frame variables
                        fillUnderOverFlow(hrecoBootstrap_ll_cHel[iPE], reco_ll_cHel_binning,
                                          std::vector<double>{ll_cHel}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_ll_cLab[iPE], reco_ll_cLab_binning,
                                          std::vector<double>{ll_cLab}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_ll_kNorm[iPE], reco_ll_kNorm_binning,
                                          std::vector<double>{ll_kNorm}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_ll_rNorm[iPE], reco_ll_rNorm_binning,
                                          std::vector<double>{ll_rNorm}, eventWeight, symmetrizeSpinVar);
                        fillUnderOverFlow(hrecoBootstrap_llbar_delta_phi[iPE], reco_ll_rNorm_binning,
                                          std::vector<double>{ll_rNorm}, eventWeight, symmetrizeSpinVar);

                        // 2D pseudo histograms
                        // fillUnderOverFlow(hrecoBootstrap_c_kk_mttbar[iPE],
                        // reco_c_kk_mttbar_binning, std::vector<double>{c_kk,
                        // ttbar_mass}, eventWeight, symmetrizeSpinVar);
                        // TODO : Fill these
                    }  // imult
                }
            }  // iPE
        }
    }

    std::cout << "n_recoentry: " << recoentry.size() << std::endl;

    // Check to see which reco matches to the corresponding gen
    std::vector<Long64_t> matched_ID_reco;
    int n_matched_ID_reco_matched = 0;
    float trueevts_0(0.0);

    // Fill all the gen histograms
    Long64_t nentries0 = fChain0->GetEntriesFast();
    std::cout << "nentries0: " << nentries0 << std::endl;

    if (!(filename.Contains("ttbarsignal") || filename.Contains("ttbarbg"))) {
        std::cout << filename << ": not ttbar so skipping gen loops" << std::endl;
        nentries0 = 0;
    }

    // *******************
    // *******************
    // Fill Gen Histograms
    // *******************
    // *******************

    // Abraham: Changed for testing. Need to revert back.
    for (Long64_t jentry0 = 0; jentry0 < nentries0; jentry0++) {
        // for (Long64_t jentry0=0; jentry0<1000;jentry0++) {
        if (LoadTree0(jentry0) < 0) break;
        fChain0->GetEntry(jentry0);

        if (!isTopGen_0) continue;

        // Amandeep : Disabling now since new minitrees have this implemented
        // Code for filtering channels
        // if (channel == "ee"   && gen_l_pdgid_0*gen_lbar_pdgid_0 != -11*11)
        // continue; if (channel == "emu"  && gen_l_pdgid_0*gen_lbar_pdgid_0 !=
        // -11*13) continue; if (channel == "mumu" &&
        // gen_l_pdgid_0*gen_lbar_pdgid_0 != -13*13) continue; End

        if (jentry0 % 100000 == 0) std::cout << "Processing gen event: " << jentry0 << std::endl;

        trueevts_0 = trueevts_0 + trueLevelWeight_0;

        // Amandeep :
        // To avoid true underflow/overflow
        if (gen_ttbar_mass_0 >= 2000.0) {
            gen_ttbar_mass_0 = 1999.0;
        }

        if (gen_ttbar_mass_0 <= 250.0) {
            gen_ttbar_mass_0 = 251.0;
        }

        if (gen_top_pt_0 >= 550.0) {
            gen_top_pt_0 = 549.0;
        }
        if (gen_n_extraJets_iso08_0 >= 3.5) {
            gen_n_extraJets_iso08_0 = 3.0;
        }
        // End

        // **************************
        // Fill 1D gen distributions
        // **************************

        // Kinematic
        fillUnderOverFlow(hgen_top_pt, gen_top_pt_binning, std::vector<double>{gen_top_pt_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_l_pt, gen_l_pt_binning, std::vector<double>{gen_l_pt_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_lbar_pt, gen_lbar_pt_binning, std::vector<double>{gen_lbar_pt_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_ttbar_pt, gen_ttbar_pt_binning, std::vector<double>{gen_ttbar_pt_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_ttbar_mass, gen_ttbar_mass_binning, std::vector<double>{gen_ttbar_mass_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_ttbar_delta_phi, gen_ttbar_delta_phi_binning, std::vector<double>{gen_ttbar_delta_phi_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_ttbar_delta_eta, gen_ttbar_delta_eta_binning, std::vector<double>{gen_ttbar_delta_eta_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_ttbar_rapidity, gen_ttbar_rapidity_binning, std::vector<double>{gen_ttbar_rapidity_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_llbar_pt, gen_llbar_pt_binning, std::vector<double>{gen_llbar_pt_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_llbar_mass, gen_llbar_mass_binning, std::vector<double>{gen_llbar_mass_0},
                          trueLevelWeight_0, symmetrizeSpinVar);

        // Spin corr

        // Polarizations
        fillUnderOverFlow(hgen_b1k, gen_b1k_binning, std::vector<double>{gen_b1k_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_b2k, gen_b2k_binning, std::vector<double>{gen_b2k_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_b1j, gen_b1j_binning, std::vector<double>{gen_b1j_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_b2j, gen_b2j_binning, std::vector<double>{gen_b2j_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_b1r, gen_b1r_binning, std::vector<double>{gen_b1r_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_b2r, gen_b2r_binning, std::vector<double>{gen_b2r_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_b1q, gen_b1q_binning, std::vector<double>{gen_b1q_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_b2q, gen_b2q_binning, std::vector<double>{gen_b2q_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_b1n, gen_b1n_binning, std::vector<double>{gen_b1n_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_b2n, gen_b2n_binning, std::vector<double>{gen_b2n_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);

        // Diagonal elements
        fillUnderOverFlow(hgen_c_kk, gen_c_kk_binning, std::vector<double>{gen_c_kk_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_rr, gen_c_rr_binning, std::vector<double>{gen_c_rr_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_nn, gen_c_nn_binning, std::vector<double>{gen_c_nn_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);

        // Off diagonal elements
        fillUnderOverFlow(hgen_c_rk, gen_c_rk_binning, std::vector<double>{gen_c_rk_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_kr, gen_c_kr_binning, std::vector<double>{gen_c_kr_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_nr, gen_c_nr_binning, std::vector<double>{gen_c_nr_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_rn, gen_c_rn_binning, std::vector<double>{gen_c_rn_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_nk, gen_c_nk_binning, std::vector<double>{gen_c_nk_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_kn, gen_c_kn_binning, std::vector<double>{gen_c_kn_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);

        fillUnderOverFlow(hgen_c_Prk, gen_c_Prk_binning, std::vector<double>{gen_c_rk_0 + gen_c_kr_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_Mrk, gen_c_Mrk_binning, std::vector<double>{gen_c_rk_0 - gen_c_kr_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_Pnr, gen_c_Pnr_binning, std::vector<double>{gen_c_nr_0 + gen_c_rn_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_Mnr, gen_c_Mnr_binning, std::vector<double>{gen_c_nr_0 - gen_c_rn_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_Pnk, gen_c_Pnk_binning, std::vector<double>{gen_c_nk_0 + gen_c_kn_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_Mnk, gen_c_Mnk_binning, std::vector<double>{gen_c_nk_0 - gen_c_kn_0},
                          trueLevelWeight_0, symmetrizeSpinVar);

        // Starred axes
        fillUnderOverFlow(hgen_c_kj, gen_c_kj_binning, std::vector<double>{gen_c_kj_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_rq, gen_c_rq_binning, std::vector<double>{gen_c_rq_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_rj, gen_c_rj_binning, std::vector<double>{gen_c_rj_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_jr, gen_c_jr_binning, std::vector<double>{gen_c_jr_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);

        fillUnderOverFlow(hgen_c_Prj, gen_c_Prj_binning, std::vector<double>{gen_c_rj_0 + gen_c_jr_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_Mrj, gen_c_Mrj_binning, std::vector<double>{gen_c_rj_0 - gen_c_jr_0},
                          trueLevelWeight_0, symmetrizeSpinVar);

        // New variables
        fillUnderOverFlow(hgen_c_han, gen_c_han_binning, std::vector<double>{+gen_c_kk_0 - gen_c_rr_0 - gen_c_nn_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_sca, gen_c_sca_binning, std::vector<double>{-gen_c_kk_0 + gen_c_rr_0 - gen_c_nn_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_tra, gen_c_tra_binning, std::vector<double>{-gen_c_kk_0 - gen_c_rr_0 + gen_c_nn_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_kjL, gen_c_kjL_binning, std::vector<double>{-gen_c_kj_0 - gen_c_rr_0 - gen_c_nn_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_rqL, gen_c_rqL_binning, std::vector<double>{-gen_c_kk_0 - gen_c_rq_0 - gen_c_nn_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_rkP, gen_c_rkP_binning, std::vector<double>{-gen_c_rk_0 - gen_c_kr_0 - gen_c_nn_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_rkM, gen_c_rkM_binning, std::vector<double>{-gen_c_rk_0 + gen_c_kr_0 - gen_c_nn_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_nrP, gen_c_nrP_binning, std::vector<double>{-gen_c_nr_0 - gen_c_rn_0 - gen_c_kk_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_nrM, gen_c_nrM_binning, std::vector<double>{-gen_c_nr_0 + gen_c_rn_0 - gen_c_kk_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_nkP, gen_c_nkP_binning, std::vector<double>{-gen_c_nk_0 - gen_c_kn_0 - gen_c_rr_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_nkM, gen_c_nkM_binning, std::vector<double>{-gen_c_nk_0 + gen_c_kn_0 - gen_c_rr_0},
                          trueLevelWeight_0, symmetrizeSpinVar);

        fillUnderOverFlow(hgen_ll_cHel, gen_ll_cHel_binning, std::vector<double>{gen_ll_cHel_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_ll_cLab, gen_ll_cLab_binning, std::vector<double>{gen_ll_cLab_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_ll_kNorm, gen_ll_kNorm_binning, std::vector<double>{gen_ll_kNorm_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_ll_rNorm, gen_ll_rNorm_binning, std::vector<double>{gen_ll_rNorm_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_llbar_delta_phi, gen_llbar_delta_phi_binning, std::vector<double>{gen_llbar_delta_phi_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_llbar_delta_eta, gen_llbar_delta_eta_binning, std::vector<double>{gen_llbar_delta_eta_0},
                          trueLevelWeight_0, symmetrizeSpinVar);

        // **************************
        // Fill 2D gen distributions
        // **************************

        // Amandeep : Fill 2D gen distributions
        // fillUnderOverFlow(hgen_c_kk_mttbar, gen_c_kk_mttbar_binning,
        // std::vector<double>{gen_c_kk_0, gen_ttbar_mass_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);

        // ******
        // mttbar
        // ******

        // Polarizations
        fillUnderOverFlow(hgen_b1k_exjet, gen_b1k_exjet_binning,
                          std::vector<double>{gen_b1k_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_b2k_exjet, gen_b2k_exjet_binning,
                          std::vector<double>{gen_b2k_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_b1j_exjet, gen_b1j_exjet_binning,
                          std::vector<double>{gen_b1j_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_b2j_exjet, gen_b2j_exjet_binning,
                          std::vector<double>{gen_b2j_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_b1r_exjet, gen_b1r_exjet_binning,
                          std::vector<double>{gen_b1r_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_b2r_exjet, gen_b2r_exjet_binning,
                          std::vector<double>{gen_b2r_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_b1q_exjet, gen_b1q_exjet_binning,
                          std::vector<double>{gen_b1q_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_b2q_exjet, gen_b2q_exjet_binning,
                          std::vector<double>{gen_b2q_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_b1n_exjet, gen_b1n_exjet_binning,
                          std::vector<double>{gen_b1n_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_b2n_exjet, gen_b2n_exjet_binning,
                          std::vector<double>{gen_b2n_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);

        // Diagonal elements
        fillUnderOverFlow(hgen_c_kk_exjet, gen_c_kk_exjet_binning,
                          std::vector<double>{gen_c_kk_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_rr_exjet, gen_c_rr_exjet_binning,
                          std::vector<double>{gen_c_rr_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_nn_exjet, gen_c_nn_exjet_binning,
                          std::vector<double>{gen_c_nn_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);

        // Off diagonal elements
        fillUnderOverFlow(hgen_c_rk_exjet, gen_c_rk_exjet_binning,
                          std::vector<double>{gen_c_rk_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_kr_exjet, gen_c_kr_exjet_binning,
                          std::vector<double>{gen_c_kr_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_nr_exjet, gen_c_nr_exjet_binning,
                          std::vector<double>{gen_c_nr_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_rn_exjet, gen_c_rn_exjet_binning,
                          std::vector<double>{gen_c_rn_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_nk_exjet, gen_c_nk_exjet_binning,
                          std::vector<double>{gen_c_nk_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_kn_exjet, gen_c_kn_exjet_binning,
                          std::vector<double>{gen_c_kn_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);

        fillUnderOverFlow(hgen_c_Prk_exjet, gen_c_Prk_exjet_binning,
                          std::vector<double>{gen_c_rk_0 + gen_c_kr_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_Mrk_exjet, gen_c_Mrk_exjet_binning,
                          std::vector<double>{gen_c_rk_0 - gen_c_kr_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_Pnr_exjet, gen_c_Pnr_exjet_binning,
                          std::vector<double>{gen_c_nr_0 + gen_c_rn_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_Mnr_exjet, gen_c_Mnr_exjet_binning,
                          std::vector<double>{gen_c_nr_0 - gen_c_rn_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_Pnk_exjet, gen_c_Pnk_exjet_binning,
                          std::vector<double>{gen_c_nk_0 + gen_c_kn_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_Mnk_exjet, gen_c_Mnk_exjet_binning,
                          std::vector<double>{gen_c_nk_0 - gen_c_kn_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);

        // Starred axes
        fillUnderOverFlow(hgen_c_kj_exjet, gen_c_kj_exjet_binning,
                          std::vector<double>{gen_c_kj_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_rq_exjet, gen_c_rq_exjet_binning,
                          std::vector<double>{gen_c_rq_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_rj_exjet, gen_c_rj_exjet_binning,
                          std::vector<double>{gen_c_rj_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_jr_exjet, gen_c_jr_exjet_binning,
                          std::vector<double>{gen_c_jr_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);

        fillUnderOverFlow(hgen_c_Prj_exjet, gen_c_Prj_exjet_binning,
                          std::vector<double>{gen_c_rj_0 + gen_c_jr_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_Mrj_exjet, gen_c_Mrj_exjet_binning,
                          std::vector<double>{gen_c_rj_0 - gen_c_jr_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);

        // New variables
        fillUnderOverFlow(hgen_c_han_exjet, gen_c_han_exjet_binning,
                          std::vector<double>{+gen_c_kk_0 - gen_c_rr_0 - gen_c_nn_0, gen_n_extraJets_iso08_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_sca_exjet, gen_c_sca_exjet_binning,
                          std::vector<double>{-gen_c_kk_0 + gen_c_rr_0 - gen_c_nn_0, gen_n_extraJets_iso08_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_tra_exjet, gen_c_tra_exjet_binning,
                          std::vector<double>{-gen_c_kk_0 - gen_c_rr_0 + gen_c_nn_0, gen_n_extraJets_iso08_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_kjL_exjet, gen_c_kjL_exjet_binning,
                          std::vector<double>{-gen_c_kj_0 - gen_c_rr_0 - gen_c_nn_0, gen_n_extraJets_iso08_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_rqL_exjet, gen_c_rqL_exjet_binning,
                          std::vector<double>{-gen_c_kk_0 - gen_c_rq_0 - gen_c_nn_0, gen_n_extraJets_iso08_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_rkP_exjet, gen_c_rkP_exjet_binning,
                          std::vector<double>{-gen_c_rk_0 - gen_c_kr_0 - gen_c_nn_0, gen_n_extraJets_iso08_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_rkM_exjet, gen_c_rkM_exjet_binning,
                          std::vector<double>{-gen_c_rk_0 + gen_c_kr_0 - gen_c_nn_0, gen_n_extraJets_iso08_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_nrP_exjet, gen_c_nrP_exjet_binning,
                          std::vector<double>{-gen_c_nr_0 - gen_c_rn_0 - gen_c_kk_0, gen_n_extraJets_iso08_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_nrM_exjet, gen_c_nrM_exjet_binning,
                          std::vector<double>{-gen_c_nr_0 + gen_c_rn_0 - gen_c_kk_0, gen_n_extraJets_iso08_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_nkP_exjet, gen_c_nkP_exjet_binning,
                          std::vector<double>{-gen_c_nk_0 - gen_c_kn_0 - gen_c_rr_0, gen_n_extraJets_iso08_0},
                          trueLevelWeight_0, symmetrizeSpinVar);
        fillUnderOverFlow(hgen_c_nkM_exjet, gen_c_nkM_exjet_binning,
                          std::vector<double>{-gen_c_nk_0 + gen_c_kn_0 - gen_c_rr_0, gen_n_extraJets_iso08_0},
                          trueLevelWeight_0, symmetrizeSpinVar);

        // Lab fram variables
        fillUnderOverFlow(hgen_ll_cHel_exjet, gen_ll_cHel_exjet_binning,
                          std::vector<double>{gen_ll_cHel_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_ll_cLab_exjet, gen_ll_cLab_exjet_binning,
                          std::vector<double>{gen_ll_cLab_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_ll_kNorm_exjet, gen_ll_kNorm_exjet_binning,
                          std::vector<double>{gen_ll_kNorm_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_ll_rNorm_exjet, gen_ll_rNorm_exjet_binning,
                          std::vector<double>{gen_ll_rNorm_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_llbar_delta_phi_exjet, gen_llbar_delta_phi_exjet_binning,
                          std::vector<double>{gen_llbar_delta_phi_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);
        fillUnderOverFlow(hgen_llbar_delta_eta_exjet, gen_llbar_delta_eta_exjet_binning,
                          std::vector<double>{gen_llbar_delta_eta_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                          symmetrizeSpinVar);

        // ******
        // top_pt
        // ******

        // Polarizations
        // fillUnderOverFlow(hgen_b1k_top_pt, gen_b1k_top_pt_binning,
        // std::vector<double>{gen_b1k_0, gen_top_pt_0}, trueLevelWeight_0,
        // symmetrizeSpinVar); fillUnderOverFlow(hgen_b2k_top_pt,
        // gen_b2k_top_pt_binning, std::vector<double>{gen_b2k_0, gen_top_pt_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_b1j_top_pt, gen_b1j_top_pt_binning,
        // std::vector<double>{gen_b1j_0, gen_top_pt_0}, trueLevelWeight_0,
        // symmetrizeSpinVar); fillUnderOverFlow(hgen_b2j_top_pt,
        // gen_b2j_top_pt_binning, std::vector<double>{gen_b2j_0, gen_top_pt_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_b1r_top_pt, gen_b1r_top_pt_binning,
        // std::vector<double>{gen_b1r_0, gen_top_pt_0}, trueLevelWeight_0,
        // symmetrizeSpinVar); fillUnderOverFlow(hgen_b2r_top_pt,
        // gen_b2r_top_pt_binning, std::vector<double>{gen_b2r_0, gen_top_pt_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_b1q_top_pt, gen_b1q_top_pt_binning,
        // std::vector<double>{gen_b1q_0, gen_top_pt_0}, trueLevelWeight_0,
        // symmetrizeSpinVar); fillUnderOverFlow(hgen_b2q_top_pt,
        // gen_b2q_top_pt_binning, std::vector<double>{gen_b2q_0, gen_top_pt_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_b1n_top_pt, gen_b1n_top_pt_binning,
        // std::vector<double>{gen_b1n_0, gen_top_pt_0}, trueLevelWeight_0,
        // symmetrizeSpinVar); fillUnderOverFlow(hgen_b2n_top_pt,
        // gen_b2n_top_pt_binning, std::vector<double>{gen_b2n_0, gen_top_pt_0},
        // trueLevelWeight_0, symmetrizeSpinVar);

        // // Diagonal elements
        // fillUnderOverFlow(hgen_c_kk_top_pt, gen_c_kk_top_pt_binning,
        // std::vector<double>{gen_c_kk_0, gen_top_pt_0}, trueLevelWeight_0,
        // symmetrizeSpinVar); fillUnderOverFlow(hgen_c_rr_top_pt,
        // gen_c_rr_top_pt_binning, std::vector<double>{gen_c_rr_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_nn_top_pt, gen_c_nn_top_pt_binning,
        // std::vector<double>{gen_c_nn_0, gen_top_pt_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);

        // // Off diagonal elements
        // fillUnderOverFlow(hgen_c_rk_top_pt, gen_c_rk_top_pt_binning,
        // std::vector<double>{gen_c_rk_0, gen_top_pt_0}, trueLevelWeight_0,
        // symmetrizeSpinVar); fillUnderOverFlow(hgen_c_kr_top_pt,
        // gen_c_kr_top_pt_binning, std::vector<double>{gen_c_kr_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_nr_top_pt, gen_c_nr_top_pt_binning,
        // std::vector<double>{gen_c_nr_0, gen_top_pt_0}, trueLevelWeight_0,
        // symmetrizeSpinVar); fillUnderOverFlow(hgen_c_rn_top_pt,
        // gen_c_rn_top_pt_binning, std::vector<double>{gen_c_rn_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_nk_top_pt, gen_c_nk_top_pt_binning,
        // std::vector<double>{gen_c_nk_0, gen_top_pt_0}, trueLevelWeight_0,
        // symmetrizeSpinVar); fillUnderOverFlow(hgen_c_kn_top_pt,
        // gen_c_kn_top_pt_binning, std::vector<double>{gen_c_kn_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);

        // fillUnderOverFlow(hgen_c_Prk_top_pt, gen_c_Prk_top_pt_binning,
        // std::vector<double>{gen_c_rk_0 + gen_c_kr_0, gen_top_pt_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_Mrk_top_pt, gen_c_Mrk_top_pt_binning,
        // std::vector<double>{gen_c_rk_0 - gen_c_kr_0, gen_top_pt_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_Pnr_top_pt, gen_c_Pnr_top_pt_binning,
        // std::vector<double>{gen_c_nr_0 + gen_c_rn_0, gen_top_pt_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_Mnr_top_pt, gen_c_Mnr_top_pt_binning,
        // std::vector<double>{gen_c_nr_0 - gen_c_rn_0, gen_top_pt_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_Pnk_top_pt, gen_c_Pnk_top_pt_binning,
        // std::vector<double>{gen_c_nk_0 + gen_c_kn_0, gen_top_pt_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_Mnk_top_pt, gen_c_Mnk_top_pt_binning,
        // std::vector<double>{gen_c_nk_0 - gen_c_kn_0, gen_top_pt_0},
        // trueLevelWeight_0, symmetrizeSpinVar);

        // // Starred axes
        // fillUnderOverFlow(hgen_c_kj_top_pt, gen_c_kj_top_pt_binning,
        // std::vector<double>{gen_c_kj_0, gen_top_pt_0}, trueLevelWeight_0,
        // symmetrizeSpinVar); fillUnderOverFlow(hgen_c_rq_top_pt,
        // gen_c_rq_top_pt_binning, std::vector<double>{gen_c_rq_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_rj_top_pt, gen_c_rj_top_pt_binning,
        // std::vector<double>{gen_c_rj_0, gen_top_pt_0}, trueLevelWeight_0,
        // symmetrizeSpinVar); fillUnderOverFlow(hgen_c_jr_top_pt,
        // gen_c_jr_top_pt_binning, std::vector<double>{gen_c_jr_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);

        // fillUnderOverFlow(hgen_c_Prj_top_pt, gen_c_Prj_top_pt_binning,
        // std::vector<double>{gen_c_rj_0 + gen_c_jr_0, gen_top_pt_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_Mrj_top_pt, gen_c_Mrj_top_pt_binning,
        // std::vector<double>{gen_c_rj_0 - gen_c_jr_0, gen_top_pt_0},
        // trueLevelWeight_0, symmetrizeSpinVar);

        // // New variables
        // fillUnderOverFlow(hgen_c_han_top_pt, gen_c_han_top_pt_binning,
        // std::vector<double>{+gen_c_kk_0 - gen_c_rr_0 - gen_c_nn_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_sca_top_pt, gen_c_sca_top_pt_binning,
        // std::vector<double>{-gen_c_kk_0 + gen_c_rr_0 - gen_c_nn_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_tra_top_pt, gen_c_tra_top_pt_binning,
        // std::vector<double>{-gen_c_kk_0 - gen_c_rr_0 + gen_c_nn_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_kjL_top_pt, gen_c_kjL_top_pt_binning,
        // std::vector<double>{-gen_c_kj_0 - gen_c_rr_0 - gen_c_nn_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_rqL_top_pt, gen_c_rqL_top_pt_binning,
        // std::vector<double>{-gen_c_kk_0 - gen_c_rq_0 - gen_c_nn_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_rkP_top_pt, gen_c_rkP_top_pt_binning,
        // std::vector<double>{-gen_c_rk_0 - gen_c_kr_0 - gen_c_nn_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_rkM_top_pt, gen_c_rkM_top_pt_binning,
        // std::vector<double>{-gen_c_rk_0 + gen_c_kr_0 - gen_c_nn_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_nrP_top_pt, gen_c_nrP_top_pt_binning,
        // std::vector<double>{-gen_c_nr_0 - gen_c_rn_0 - gen_c_kk_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_nrM_top_pt, gen_c_nrM_top_pt_binning,
        // std::vector<double>{-gen_c_nr_0 + gen_c_rn_0 - gen_c_kk_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_nkP_top_pt, gen_c_nkP_top_pt_binning,
        // std::vector<double>{-gen_c_nk_0 - gen_c_kn_0 - gen_c_rr_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_nkM_top_pt, gen_c_nkM_top_pt_binning,
        // std::vector<double>{-gen_c_nk_0 + gen_c_kn_0 - gen_c_rr_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);

        // // Lab fram variables
        // fillUnderOverFlow(hgen_ll_cHel_top_pt, gen_ll_cHel_top_pt_binning,
        // std::vector<double>{gen_ll_cHel_0, gen_top_pt_0}, trueLevelWeight_0,
        // symmetrizeSpinVar); fillUnderOverFlow(hgen_ll_cLab_top_pt,
        // gen_ll_cLab_top_pt_binning, std::vector<double>{gen_ll_cLab_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_ll_kNorm_top_pt, gen_ll_kNorm_top_pt_binning,
        // std::vector<double>{gen_ll_kNorm_0, gen_top_pt_0}, trueLevelWeight_0,
        // symmetrizeSpinVar); fillUnderOverFlow(hgen_ll_rNorm_top_pt,
        // gen_ll_rNorm_top_pt_binning, std::vector<double>{gen_ll_rNorm_0,
        // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_llbar_delta_phi_top_pt,
        // gen_llbar_delta_phi_top_pt_binning,
        // std::vector<double>{gen_llbar_delta_phi_0, gen_top_pt_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_llbar_delta_eta_top_pt,
        // gen_llbar_delta_eta_top_pt_binning,
        // std::vector<double>{gen_llbar_delta_eta_0, gen_top_pt_0},
        // trueLevelWeight_0, symmetrizeSpinVar);

        // ******************************
        // top_scatteringangle_ttbarframe
        // ******************************

        // Polarizations
        // fillUnderOverFlow(hgen_b1k_top_scatteringangle_ttbarframe,
        // gen_b1k_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_b1k_0, gen_top_scatteringangle_ttbarframe_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_b2k_top_scatteringangle_ttbarframe,
        // gen_b2k_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_b2k_0, gen_top_scatteringangle_ttbarframe_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_b1j_top_scatteringangle_ttbarframe,
        // gen_b1j_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_b1j_0, gen_top_scatteringangle_ttbarframe_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_b2j_top_scatteringangle_ttbarframe,
        // gen_b2j_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_b2j_0, gen_top_scatteringangle_ttbarframe_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_b1r_top_scatteringangle_ttbarframe,
        // gen_b1r_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_b1r_0, gen_top_scatteringangle_ttbarframe_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_b2r_top_scatteringangle_ttbarframe,
        // gen_b2r_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_b2r_0, gen_top_scatteringangle_ttbarframe_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_b1q_top_scatteringangle_ttbarframe,
        // gen_b1q_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_b1q_0, gen_top_scatteringangle_ttbarframe_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_b2q_top_scatteringangle_ttbarframe,
        // gen_b2q_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_b2q_0, gen_top_scatteringangle_ttbarframe_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_b1n_top_scatteringangle_ttbarframe,
        // gen_b1n_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_b1n_0, gen_top_scatteringangle_ttbarframe_0},
        // trueLevelWeight_0, symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_b2n_top_scatteringangle_ttbarframe,
        // gen_b2n_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_b2n_0, gen_top_scatteringangle_ttbarframe_0},
        // trueLevelWeight_0, symmetrizeSpinVar);

        // // Diagonal elements
        // fillUnderOverFlow(hgen_c_kk_top_scatteringangle_ttbarframe,
        // gen_c_kk_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_kk_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_rr_top_scatteringangle_ttbarframe,
        // gen_c_rr_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_rr_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_nn_top_scatteringangle_ttbarframe,
        // gen_c_nn_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_nn_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);

        // // Off diagonal elements
        // fillUnderOverFlow(hgen_c_rk_top_scatteringangle_ttbarframe,
        // gen_c_rk_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_rk_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_kr_top_scatteringangle_ttbarframe,
        // gen_c_kr_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_kr_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_nr_top_scatteringangle_ttbarframe,
        // gen_c_nr_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_nr_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_rn_top_scatteringangle_ttbarframe,
        // gen_c_rn_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_rn_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_nk_top_scatteringangle_ttbarframe,
        // gen_c_nk_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_nk_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_kn_top_scatteringangle_ttbarframe,
        // gen_c_kn_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_kn_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);

        // fillUnderOverFlow(hgen_c_Prk_top_scatteringangle_ttbarframe,
        // gen_c_Prk_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_rk_0 + gen_c_kr_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_Mrk_top_scatteringangle_ttbarframe,
        // gen_c_Mrk_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_rk_0 - gen_c_kr_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_Pnr_top_scatteringangle_ttbarframe,
        // gen_c_Pnr_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_nr_0 + gen_c_rn_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_Mnr_top_scatteringangle_ttbarframe,
        // gen_c_Mnr_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_nr_0 - gen_c_rn_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_Pnk_top_scatteringangle_ttbarframe,
        // gen_c_Pnk_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_nk_0 + gen_c_kn_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_Mnk_top_scatteringangle_ttbarframe,
        // gen_c_Mnk_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_nk_0 - gen_c_kn_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);

        // // Starred axes
        // fillUnderOverFlow(hgen_c_kj_top_scatteringangle_ttbarframe,
        // gen_c_kj_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_kj_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_rq_top_scatteringangle_ttbarframe,
        // gen_c_rq_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_rq_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_rj_top_scatteringangle_ttbarframe,
        // gen_c_rj_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_rj_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_jr_top_scatteringangle_ttbarframe,
        // gen_c_jr_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_jr_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);

        // fillUnderOverFlow(hgen_c_Prj_top_scatteringangle_ttbarframe,
        // gen_c_Prj_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_rj_0 + gen_c_jr_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_Mrj_top_scatteringangle_ttbarframe,
        // gen_c_Mrj_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_c_rj_0 - gen_c_jr_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);

        // // New variables
        // fillUnderOverFlow(hgen_c_han_top_scatteringangle_ttbarframe,
        // gen_c_han_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{+gen_c_kk_0 - gen_c_rr_0 - gen_c_nn_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_sca_top_scatteringangle_ttbarframe,
        // gen_c_sca_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-gen_c_kk_0 + gen_c_rr_0 - gen_c_nn_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_tra_top_scatteringangle_ttbarframe,
        // gen_c_tra_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-gen_c_kk_0 - gen_c_rr_0 + gen_c_nn_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_kjL_top_scatteringangle_ttbarframe,
        // gen_c_kjL_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-gen_c_kj_0 - gen_c_rr_0 - gen_c_nn_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_rqL_top_scatteringangle_ttbarframe,
        // gen_c_rqL_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-gen_c_kk_0 - gen_c_rq_0 - gen_c_nn_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_rkP_top_scatteringangle_ttbarframe,
        // gen_c_rkP_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-gen_c_rk_0 - gen_c_kr_0 - gen_c_nn_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_rkM_top_scatteringangle_ttbarframe,
        // gen_c_rkM_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-gen_c_rk_0 + gen_c_kr_0 - gen_c_nn_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_nrP_top_scatteringangle_ttbarframe,
        // gen_c_nrP_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-gen_c_nr_0 - gen_c_rn_0 - gen_c_kk_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_nrM_top_scatteringangle_ttbarframe,
        // gen_c_nrM_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-gen_c_nr_0 + gen_c_rn_0 - gen_c_kk_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_nkP_top_scatteringangle_ttbarframe,
        // gen_c_nkP_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-gen_c_nk_0 - gen_c_kn_0 - gen_c_rr_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_c_nkM_top_scatteringangle_ttbarframe,
        // gen_c_nkM_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{-gen_c_nk_0 + gen_c_kn_0 - gen_c_rr_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);

        // // Lab fram variables
        // fillUnderOverFlow(hgen_ll_cHel_top_scatteringangle_ttbarframe,
        // gen_ll_cHel_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_ll_cHel_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_ll_cLab_top_scatteringangle_ttbarframe,
        // gen_ll_cLab_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_ll_cLab_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_ll_kNorm_top_scatteringangle_ttbarframe,
        // gen_ll_kNorm_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_ll_kNorm_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_ll_rNorm_top_scatteringangle_ttbarframe,
        // gen_ll_rNorm_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_ll_rNorm_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_llbar_delta_phi_top_scatteringangle_ttbarframe,
        // gen_llbar_delta_phi_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_llbar_delta_phi_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);
        // fillUnderOverFlow(hgen_llbar_delta_eta_top_scatteringangle_ttbarframe,
        // gen_llbar_delta_eta_top_scatteringangle_ttbarframe_binning,
        // std::vector<double>{gen_llbar_delta_eta_0,
        // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
        // symmetrizeSpinVar);

        // **********************
        // **********************
        // Fill VisGen Histograms
        // **********************
        // **********************

        // *****************************
        // Fill 1D vis-gen distributions
        // *****************************

        // Apply some phase space cuts

        if ((gen_l_pt_0 > 20 && fabs(gen_l_eta) < 2.4) && (gen_lbar_pt_0 > 20 && fabs(gen_lbar_eta) < 2.4) &&
            (gen_b_pt_0 > 30 && fabs(gen_b_eta) < 2.4) && (gen_bbar_pt_0 > 30 && fabs(gen_bbar_eta) < 2.4)) {
            // Kinematic
            fillUnderOverFlow(hvisgen_top_pt, gen_top_pt_binning, std::vector<double>{gen_top_pt_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_l_pt, gen_l_pt_binning, std::vector<double>{gen_l_pt_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_lbar_pt, gen_lbar_pt_binning, std::vector<double>{gen_lbar_pt_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_ttbar_pt, gen_ttbar_pt_binning, std::vector<double>{gen_ttbar_pt_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_ttbar_mass, gen_ttbar_mass_binning, std::vector<double>{gen_ttbar_mass_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_ttbar_delta_phi, gen_ttbar_delta_phi_binning,
                              std::vector<double>{gen_ttbar_delta_phi_0}, trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_ttbar_delta_eta, gen_ttbar_delta_eta_binning,
                              std::vector<double>{gen_ttbar_delta_eta_0}, trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_ttbar_rapidity, gen_ttbar_rapidity_binning,
                              std::vector<double>{gen_ttbar_rapidity_0}, trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_llbar_pt, gen_llbar_pt_binning, std::vector<double>{gen_llbar_pt_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_llbar_mass, gen_llbar_mass_binning, std::vector<double>{gen_llbar_mass_0},
                              trueLevelWeight_0, symmetrizeSpinVar);

            // Spin corr

            // Polarizations
            fillUnderOverFlow(hvisgen_b1k, gen_b1k_binning, std::vector<double>{gen_b1k_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_b2k, gen_b2k_binning, std::vector<double>{gen_b2k_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_b1j, gen_b1j_binning, std::vector<double>{gen_b1j_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_b2j, gen_b2j_binning, std::vector<double>{gen_b2j_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_b1r, gen_b1r_binning, std::vector<double>{gen_b1r_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_b2r, gen_b2r_binning, std::vector<double>{gen_b2r_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_b1q, gen_b1q_binning, std::vector<double>{gen_b1q_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_b2q, gen_b2q_binning, std::vector<double>{gen_b2q_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_b1n, gen_b1n_binning, std::vector<double>{gen_b1n_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_b2n, gen_b2n_binning, std::vector<double>{gen_b2n_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);

            // Diagonal elements
            fillUnderOverFlow(hvisgen_c_kk, gen_c_kk_binning, std::vector<double>{gen_c_kk_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_rr, gen_c_rr_binning, std::vector<double>{gen_c_rr_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_nn, gen_c_nn_binning, std::vector<double>{gen_c_nn_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);

            // Off diagonal elements
            fillUnderOverFlow(hvisgen_c_rk, gen_c_rk_binning, std::vector<double>{gen_c_rk_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_kr, gen_c_kr_binning, std::vector<double>{gen_c_kr_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_nr, gen_c_nr_binning, std::vector<double>{gen_c_nr_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_rn, gen_c_rn_binning, std::vector<double>{gen_c_rn_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_nk, gen_c_nk_binning, std::vector<double>{gen_c_nk_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_kn, gen_c_kn_binning, std::vector<double>{gen_c_kn_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);

            fillUnderOverFlow(hvisgen_c_Prk, gen_c_Prk_binning, std::vector<double>{gen_c_rk_0 + gen_c_kr_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_Mrk, gen_c_Mrk_binning, std::vector<double>{gen_c_rk_0 - gen_c_kr_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_Pnr, gen_c_Pnr_binning, std::vector<double>{gen_c_nr_0 + gen_c_rn_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_Mnr, gen_c_Mnr_binning, std::vector<double>{gen_c_nr_0 - gen_c_rn_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_Pnk, gen_c_Pnk_binning, std::vector<double>{gen_c_nk_0 + gen_c_kn_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_Mnk, gen_c_Mnk_binning, std::vector<double>{gen_c_nk_0 - gen_c_kn_0},
                              trueLevelWeight_0, symmetrizeSpinVar);

            // Starred axes
            fillUnderOverFlow(hvisgen_c_kj, gen_c_kj_binning, std::vector<double>{gen_c_kj_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_rq, gen_c_rq_binning, std::vector<double>{gen_c_rq_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_rj, gen_c_rj_binning, std::vector<double>{gen_c_rj_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_jr, gen_c_jr_binning, std::vector<double>{gen_c_jr_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);

            fillUnderOverFlow(hvisgen_c_Prj, gen_c_Prj_binning, std::vector<double>{gen_c_rj_0 + gen_c_jr_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_Mrj, gen_c_Mrj_binning, std::vector<double>{gen_c_rj_0 - gen_c_jr_0},
                              trueLevelWeight_0, symmetrizeSpinVar);

            // New variables
            fillUnderOverFlow(hvisgen_c_han, gen_c_han_binning,
                              std::vector<double>{+gen_c_kk_0 - gen_c_rr_0 - gen_c_nn_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_sca, gen_c_sca_binning,
                              std::vector<double>{-gen_c_kk_0 + gen_c_rr_0 - gen_c_nn_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_tra, gen_c_tra_binning,
                              std::vector<double>{-gen_c_kk_0 - gen_c_rr_0 + gen_c_nn_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_kjL, gen_c_kjL_binning,
                              std::vector<double>{-gen_c_kj_0 - gen_c_rr_0 - gen_c_nn_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_rqL, gen_c_rqL_binning,
                              std::vector<double>{-gen_c_kk_0 - gen_c_rq_0 - gen_c_nn_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_rkP, gen_c_rkP_binning,
                              std::vector<double>{-gen_c_rk_0 - gen_c_kr_0 - gen_c_nn_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_rkM, gen_c_rkM_binning,
                              std::vector<double>{-gen_c_rk_0 + gen_c_kr_0 - gen_c_nn_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_nrP, gen_c_nrP_binning,
                              std::vector<double>{-gen_c_nr_0 - gen_c_rn_0 - gen_c_kk_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_nrM, gen_c_nrM_binning,
                              std::vector<double>{-gen_c_nr_0 + gen_c_rn_0 - gen_c_kk_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_nkP, gen_c_nkP_binning,
                              std::vector<double>{-gen_c_nk_0 - gen_c_kn_0 - gen_c_rr_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_nkM, gen_c_nkM_binning,
                              std::vector<double>{-gen_c_nk_0 + gen_c_kn_0 - gen_c_rr_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);

            // Lab frame variables
            fillUnderOverFlow(hvisgen_ll_cHel, gen_ll_cHel_binning, std::vector<double>{gen_ll_cHel_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_ll_cLab, gen_ll_cLab_binning, std::vector<double>{gen_ll_cLab_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_ll_kNorm, gen_ll_kNorm_binning, std::vector<double>{gen_ll_kNorm_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_ll_rNorm, gen_ll_rNorm_binning, std::vector<double>{gen_ll_rNorm_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_llbar_delta_phi, gen_llbar_delta_phi_binning,
                              std::vector<double>{gen_llbar_delta_phi_0}, trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_llbar_delta_eta, gen_llbar_delta_eta_binning,
                              std::vector<double>{gen_llbar_delta_eta_0}, trueLevelWeight_0, symmetrizeSpinVar);

            // *****************************
            // Fill 2D vis-gen distributions
            // *****************************

            // Amandeep : Fill 2D vis-gen histograms
            // fillUnderOverFlow(hvisgen_c_kk_mttbar, gen_c_kk_mttbar_binning,
            // std::vector<double>{gen_c_kk_0, gen_ttbar_mass_0},
            // trueLevelWeight_0, symmetrizeSpinVar);

            // ******
            // mttbar
            // ******

            // Polarizations
            fillUnderOverFlow(hvisgen_b1k_exjet, gen_b1k_exjet_binning,
                              std::vector<double>{gen_b1k_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_b2k_exjet, gen_b2k_exjet_binning,
                              std::vector<double>{gen_b2k_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_b1j_exjet, gen_b1j_exjet_binning,
                              std::vector<double>{gen_b1j_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_b2j_exjet, gen_b2j_exjet_binning,
                              std::vector<double>{gen_b2j_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_b1r_exjet, gen_b1r_exjet_binning,
                              std::vector<double>{gen_b1r_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_b2r_exjet, gen_b2r_exjet_binning,
                              std::vector<double>{gen_b2r_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_b1q_exjet, gen_b1q_exjet_binning,
                              std::vector<double>{gen_b1q_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_b2q_exjet, gen_b2q_exjet_binning,
                              std::vector<double>{gen_b2q_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_b1n_exjet, gen_b1n_exjet_binning,
                              std::vector<double>{gen_b1n_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_b2n_exjet, gen_b2n_exjet_binning,
                              std::vector<double>{gen_b2n_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);

            // Diagonal elements
            fillUnderOverFlow(hvisgen_c_kk_exjet, gen_c_kk_exjet_binning,
                              std::vector<double>{gen_c_kk_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_rr_exjet, gen_c_rr_exjet_binning,
                              std::vector<double>{gen_c_rr_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_nn_exjet, gen_c_nn_exjet_binning,
                              std::vector<double>{gen_c_nn_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);

            // Off diagonal elements
            fillUnderOverFlow(hvisgen_c_rk_exjet, gen_c_rk_exjet_binning,
                              std::vector<double>{gen_c_rk_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_kr_exjet, gen_c_kr_exjet_binning,
                              std::vector<double>{gen_c_kr_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_nr_exjet, gen_c_nr_exjet_binning,
                              std::vector<double>{gen_c_nr_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_rn_exjet, gen_c_rn_exjet_binning,
                              std::vector<double>{gen_c_rn_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_nk_exjet, gen_c_nk_exjet_binning,
                              std::vector<double>{gen_c_nk_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_kn_exjet, gen_c_kn_exjet_binning,
                              std::vector<double>{gen_c_kn_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);

            fillUnderOverFlow(hvisgen_c_Prk_exjet, gen_c_Prk_exjet_binning,
                              std::vector<double>{gen_c_rk_0 + gen_c_kr_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_Mrk_exjet, gen_c_Mrk_exjet_binning,
                              std::vector<double>{gen_c_rk_0 - gen_c_kr_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_Pnr_exjet, gen_c_Pnr_exjet_binning,
                              std::vector<double>{gen_c_nr_0 + gen_c_rn_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_Mnr_exjet, gen_c_Mnr_exjet_binning,
                              std::vector<double>{gen_c_nr_0 - gen_c_rn_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_Pnk_exjet, gen_c_Pnk_exjet_binning,
                              std::vector<double>{gen_c_nk_0 + gen_c_kn_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_Mnk_exjet, gen_c_Mnk_exjet_binning,
                              std::vector<double>{gen_c_nk_0 - gen_c_kn_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);

            // Starred axes
            fillUnderOverFlow(hvisgen_c_kj_exjet, gen_c_kj_exjet_binning,
                              std::vector<double>{gen_c_kj_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_rq_exjet, gen_c_rq_exjet_binning,
                              std::vector<double>{gen_c_rq_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_rj_exjet, gen_c_rj_exjet_binning,
                              std::vector<double>{gen_c_rj_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_jr_exjet, gen_c_jr_exjet_binning,
                              std::vector<double>{gen_c_jr_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);

            fillUnderOverFlow(hvisgen_c_Prj_exjet, gen_c_Prj_exjet_binning,
                              std::vector<double>{gen_c_rj_0 + gen_c_jr_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_Mrj_exjet, gen_c_Mrj_exjet_binning,
                              std::vector<double>{gen_c_rj_0 - gen_c_jr_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);

            // New variables
            fillUnderOverFlow(hvisgen_c_han_exjet, gen_c_han_exjet_binning,
                              std::vector<double>{+gen_c_kk_0 - gen_c_rr_0 - gen_c_nn_0, gen_n_extraJets_iso08_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_sca_exjet, gen_c_sca_exjet_binning,
                              std::vector<double>{-gen_c_kk_0 + gen_c_rr_0 - gen_c_nn_0, gen_n_extraJets_iso08_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_tra_exjet, gen_c_tra_exjet_binning,
                              std::vector<double>{-gen_c_kk_0 - gen_c_rr_0 + gen_c_nn_0, gen_n_extraJets_iso08_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_kjL_exjet, gen_c_kjL_exjet_binning,
                              std::vector<double>{-gen_c_kj_0 - gen_c_rr_0 - gen_c_nn_0, gen_n_extraJets_iso08_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_rqL_exjet, gen_c_rqL_exjet_binning,
                              std::vector<double>{-gen_c_kk_0 - gen_c_rq_0 - gen_c_nn_0, gen_n_extraJets_iso08_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_rkP_exjet, gen_c_rkP_exjet_binning,
                              std::vector<double>{-gen_c_rk_0 - gen_c_kr_0 - gen_c_nn_0, gen_n_extraJets_iso08_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_rkM_exjet, gen_c_rkM_exjet_binning,
                              std::vector<double>{-gen_c_rk_0 + gen_c_kr_0 - gen_c_nn_0, gen_n_extraJets_iso08_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_nrP_exjet, gen_c_nrP_exjet_binning,
                              std::vector<double>{-gen_c_nr_0 - gen_c_rn_0 - gen_c_kk_0, gen_n_extraJets_iso08_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_nrM_exjet, gen_c_nrM_exjet_binning,
                              std::vector<double>{-gen_c_nr_0 + gen_c_rn_0 - gen_c_kk_0, gen_n_extraJets_iso08_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_nkP_exjet, gen_c_nkP_exjet_binning,
                              std::vector<double>{-gen_c_nk_0 - gen_c_kn_0 - gen_c_rr_0, gen_n_extraJets_iso08_0},
                              trueLevelWeight_0, symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_c_nkM_exjet, gen_c_nkM_exjet_binning,
                              std::vector<double>{-gen_c_nk_0 + gen_c_kn_0 - gen_c_rr_0, gen_n_extraJets_iso08_0},
                              trueLevelWeight_0, symmetrizeSpinVar);

            // Lab frame variables
            fillUnderOverFlow(hvisgen_ll_cHel_exjet, gen_ll_cHel_exjet_binning,
                              std::vector<double>{gen_ll_cHel_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_ll_cLab_exjet, gen_ll_cLab_exjet_binning,
                              std::vector<double>{gen_ll_cLab_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_ll_kNorm_exjet, gen_ll_kNorm_exjet_binning,
                              std::vector<double>{gen_ll_kNorm_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_ll_rNorm_exjet, gen_ll_rNorm_exjet_binning,
                              std::vector<double>{gen_ll_rNorm_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_llbar_delta_phi_exjet, gen_llbar_delta_phi_exjet_binning,
                              std::vector<double>{gen_llbar_delta_phi_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hvisgen_llbar_delta_eta_exjet, gen_llbar_delta_eta_exjet_binning,
                              std::vector<double>{gen_llbar_delta_eta_0, gen_n_extraJets_iso08_0}, trueLevelWeight_0,
                              symmetrizeSpinVar);

            // ******
            // top_pt
            // ******

            // // Polarizations
            // fillUnderOverFlow(hvisgen_b1k_top_pt, gen_b1k_top_pt_binning,
            // std::vector<double>{gen_b1k_0, gen_top_pt_0}, trueLevelWeight_0,
            // symmetrizeSpinVar); fillUnderOverFlow(hvisgen_b2k_top_pt,
            // gen_b2k_top_pt_binning, std::vector<double>{gen_b2k_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_b1j_top_pt, gen_b1j_top_pt_binning,
            // std::vector<double>{gen_b1j_0, gen_top_pt_0}, trueLevelWeight_0,
            // symmetrizeSpinVar); fillUnderOverFlow(hvisgen_b2j_top_pt,
            // gen_b2j_top_pt_binning, std::vector<double>{gen_b2j_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_b1r_top_pt, gen_b1r_top_pt_binning,
            // std::vector<double>{gen_b1r_0, gen_top_pt_0}, trueLevelWeight_0,
            // symmetrizeSpinVar); fillUnderOverFlow(hvisgen_b2r_top_pt,
            // gen_b2r_top_pt_binning, std::vector<double>{gen_b2r_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_b1q_top_pt, gen_b1q_top_pt_binning,
            // std::vector<double>{gen_b1q_0, gen_top_pt_0}, trueLevelWeight_0,
            // symmetrizeSpinVar); fillUnderOverFlow(hvisgen_b2q_top_pt,
            // gen_b2q_top_pt_binning, std::vector<double>{gen_b2q_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_b1n_top_pt, gen_b1n_top_pt_binning,
            // std::vector<double>{gen_b1n_0, gen_top_pt_0}, trueLevelWeight_0,
            // symmetrizeSpinVar); fillUnderOverFlow(hvisgen_b2n_top_pt,
            // gen_b2n_top_pt_binning, std::vector<double>{gen_b2n_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);

            // // Diagonal elements
            // fillUnderOverFlow(hvisgen_c_kk_top_pt, gen_c_kk_top_pt_binning,
            // std::vector<double>{gen_c_kk_0, gen_top_pt_0}, trueLevelWeight_0,
            // symmetrizeSpinVar); fillUnderOverFlow(hvisgen_c_rr_top_pt,
            // gen_c_rr_top_pt_binning, std::vector<double>{gen_c_rr_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_nn_top_pt, gen_c_nn_top_pt_binning,
            // std::vector<double>{gen_c_nn_0, gen_top_pt_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);

            // // Off diagonal elements
            // fillUnderOverFlow(hvisgen_c_rk_top_pt, gen_c_rk_top_pt_binning,
            // std::vector<double>{gen_c_rk_0, gen_top_pt_0}, trueLevelWeight_0,
            // symmetrizeSpinVar); fillUnderOverFlow(hvisgen_c_kr_top_pt,
            // gen_c_kr_top_pt_binning, std::vector<double>{gen_c_kr_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_nr_top_pt, gen_c_nr_top_pt_binning,
            // std::vector<double>{gen_c_nr_0, gen_top_pt_0}, trueLevelWeight_0,
            // symmetrizeSpinVar); fillUnderOverFlow(hvisgen_c_rn_top_pt,
            // gen_c_rn_top_pt_binning, std::vector<double>{gen_c_rn_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_nk_top_pt, gen_c_nk_top_pt_binning,
            // std::vector<double>{gen_c_nk_0, gen_top_pt_0}, trueLevelWeight_0,
            // symmetrizeSpinVar); fillUnderOverFlow(hvisgen_c_kn_top_pt,
            // gen_c_kn_top_pt_binning, std::vector<double>{gen_c_kn_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);

            // fillUnderOverFlow(hvisgen_c_Prk_top_pt, gen_c_Prk_top_pt_binning,
            // std::vector<double>{gen_c_rk_0 + gen_c_kr_0, gen_top_pt_0},
            // trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_Mrk_top_pt, gen_c_Mrk_top_pt_binning,
            // std::vector<double>{gen_c_rk_0 - gen_c_kr_0, gen_top_pt_0},
            // trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_Pnr_top_pt, gen_c_Pnr_top_pt_binning,
            // std::vector<double>{gen_c_nr_0 + gen_c_rn_0, gen_top_pt_0},
            // trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_Mnr_top_pt, gen_c_Mnr_top_pt_binning,
            // std::vector<double>{gen_c_nr_0 - gen_c_rn_0, gen_top_pt_0},
            // trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_Pnk_top_pt, gen_c_Pnk_top_pt_binning,
            // std::vector<double>{gen_c_nk_0 + gen_c_kn_0, gen_top_pt_0},
            // trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_Mnk_top_pt, gen_c_Mnk_top_pt_binning,
            // std::vector<double>{gen_c_nk_0 - gen_c_kn_0, gen_top_pt_0},
            // trueLevelWeight_0, symmetrizeSpinVar);

            // // Starred axes
            // fillUnderOverFlow(hvisgen_c_kj_top_pt, gen_c_kj_top_pt_binning,
            // std::vector<double>{gen_c_kj_0, gen_top_pt_0}, trueLevelWeight_0,
            // symmetrizeSpinVar); fillUnderOverFlow(hvisgen_c_rq_top_pt,
            // gen_c_rq_top_pt_binning, std::vector<double>{gen_c_rq_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_rj_top_pt, gen_c_rj_top_pt_binning,
            // std::vector<double>{gen_c_rj_0, gen_top_pt_0}, trueLevelWeight_0,
            // symmetrizeSpinVar); fillUnderOverFlow(hvisgen_c_jr_top_pt,
            // gen_c_jr_top_pt_binning, std::vector<double>{gen_c_jr_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);

            // fillUnderOverFlow(hvisgen_c_Prj_top_pt, gen_c_Prj_top_pt_binning,
            // std::vector<double>{gen_c_rj_0 + gen_c_jr_0, gen_top_pt_0},
            // trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_Mrj_top_pt, gen_c_Mrj_top_pt_binning,
            // std::vector<double>{gen_c_rj_0 - gen_c_jr_0, gen_top_pt_0},
            // trueLevelWeight_0, symmetrizeSpinVar);

            // // New variables
            // fillUnderOverFlow(hvisgen_c_han_top_pt, gen_c_han_top_pt_binning,
            // std::vector<double>{+gen_c_kk_0 - gen_c_rr_0 - gen_c_nn_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_sca_top_pt, gen_c_sca_top_pt_binning,
            // std::vector<double>{-gen_c_kk_0 + gen_c_rr_0 - gen_c_nn_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_tra_top_pt, gen_c_tra_top_pt_binning,
            // std::vector<double>{-gen_c_kk_0 - gen_c_rr_0 + gen_c_nn_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_kjL_top_pt, gen_c_kjL_top_pt_binning,
            // std::vector<double>{-gen_c_kj_0 - gen_c_rr_0 - gen_c_nn_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_rqL_top_pt, gen_c_rqL_top_pt_binning,
            // std::vector<double>{-gen_c_kk_0 - gen_c_rq_0 - gen_c_nn_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_rkP_top_pt, gen_c_rkP_top_pt_binning,
            // std::vector<double>{-gen_c_rk_0 - gen_c_kr_0 - gen_c_nn_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_rkM_top_pt, gen_c_rkM_top_pt_binning,
            // std::vector<double>{-gen_c_rk_0 + gen_c_kr_0 - gen_c_nn_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_nrP_top_pt, gen_c_nrP_top_pt_binning,
            // std::vector<double>{-gen_c_nr_0 - gen_c_rn_0 - gen_c_kk_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_nrM_top_pt, gen_c_nrM_top_pt_binning,
            // std::vector<double>{-gen_c_nr_0 + gen_c_rn_0 - gen_c_kk_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_nkP_top_pt, gen_c_nkP_top_pt_binning,
            // std::vector<double>{-gen_c_nk_0 - gen_c_kn_0 - gen_c_rr_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_nkM_top_pt, gen_c_nkM_top_pt_binning,
            // std::vector<double>{-gen_c_nk_0 + gen_c_kn_0 - gen_c_rr_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);

            // // Lab frame variables
            // fillUnderOverFlow(hvisgen_ll_cHel_top_pt,
            // gen_ll_cHel_top_pt_binning, std::vector<double>{gen_ll_cHel_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_ll_cLab_top_pt,
            // gen_ll_cLab_top_pt_binning, std::vector<double>{gen_ll_cLab_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_ll_kNorm_top_pt,
            // gen_ll_kNorm_top_pt_binning, std::vector<double>{gen_ll_kNorm_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_ll_rNorm_top_pt,
            // gen_ll_rNorm_top_pt_binning, std::vector<double>{gen_ll_rNorm_0,
            // gen_top_pt_0}, trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_llbar_delta_phi_top_pt,
            // gen_llbar_delta_phi_top_pt_binning,
            // std::vector<double>{gen_llbar_delta_phi_0, gen_top_pt_0},
            // trueLevelWeight_0, symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_llbar_delta_eta_top_pt,
            // gen_llbar_delta_eta_top_pt_binning,
            // std::vector<double>{gen_llbar_delta_eta_0, gen_top_pt_0},
            // trueLevelWeight_0, symmetrizeSpinVar);

            // ******************************
            // top_scatteringangle_ttbarframe
            // ******************************

            // // Polarizations
            // fillUnderOverFlow(hvisgen_b1k_top_scatteringangle_ttbarframe,
            // gen_b1k_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1k_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_b2k_top_scatteringangle_ttbarframe,
            // gen_b2k_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2k_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_b1j_top_scatteringangle_ttbarframe,
            // gen_b1j_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1j_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_b2j_top_scatteringangle_ttbarframe,
            // gen_b2j_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2j_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_b1r_top_scatteringangle_ttbarframe,
            // gen_b1r_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1r_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_b2r_top_scatteringangle_ttbarframe,
            // gen_b2r_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2r_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_b1q_top_scatteringangle_ttbarframe,
            // gen_b1q_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1q_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_b2q_top_scatteringangle_ttbarframe,
            // gen_b2q_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2q_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_b1n_top_scatteringangle_ttbarframe,
            // gen_b1n_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1n_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_b2n_top_scatteringangle_ttbarframe,
            // gen_b2n_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2n_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);

            // // Diagonal elements
            // fillUnderOverFlow(hvisgen_c_kk_top_scatteringangle_ttbarframe,
            // gen_c_kk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_kk_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_rr_top_scatteringangle_ttbarframe,
            // gen_c_rr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rr_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_nn_top_scatteringangle_ttbarframe,
            // gen_c_nn_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nn_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);

            // // Off diagonal elements
            // fillUnderOverFlow(hvisgen_c_rk_top_scatteringangle_ttbarframe,
            // gen_c_rk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rk_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_kr_top_scatteringangle_ttbarframe,
            // gen_c_kr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_kr_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_nr_top_scatteringangle_ttbarframe,
            // gen_c_nr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nr_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_rn_top_scatteringangle_ttbarframe,
            // gen_c_rn_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rn_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_nk_top_scatteringangle_ttbarframe,
            // gen_c_nk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nk_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_kn_top_scatteringangle_ttbarframe,
            // gen_c_kn_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_kn_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);

            // fillUnderOverFlow(hvisgen_c_Prk_top_scatteringangle_ttbarframe,
            // gen_c_Prk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rk_0 + gen_c_kr_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_Mrk_top_scatteringangle_ttbarframe,
            // gen_c_Mrk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rk_0 - gen_c_kr_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_Pnr_top_scatteringangle_ttbarframe,
            // gen_c_Pnr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nr_0 + gen_c_rn_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_Mnr_top_scatteringangle_ttbarframe,
            // gen_c_Mnr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nr_0 - gen_c_rn_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_Pnk_top_scatteringangle_ttbarframe,
            // gen_c_Pnk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nk_0 + gen_c_kn_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_Mnk_top_scatteringangle_ttbarframe,
            // gen_c_Mnk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nk_0 - gen_c_kn_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);

            // // Starred axes
            // fillUnderOverFlow(hvisgen_c_kj_top_scatteringangle_ttbarframe,
            // gen_c_kj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_kj_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_rq_top_scatteringangle_ttbarframe,
            // gen_c_rq_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rq_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_rj_top_scatteringangle_ttbarframe,
            // gen_c_rj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rj_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_jr_top_scatteringangle_ttbarframe,
            // gen_c_jr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_jr_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);

            // fillUnderOverFlow(hvisgen_c_Prj_top_scatteringangle_ttbarframe,
            // gen_c_Prj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rj_0 + gen_c_jr_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_Mrj_top_scatteringangle_ttbarframe,
            // gen_c_Mrj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rj_0 - gen_c_jr_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);

            // // New variables
            // fillUnderOverFlow(hvisgen_c_han_top_scatteringangle_ttbarframe,
            // gen_c_han_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{+gen_c_kk_0 - gen_c_rr_0 - gen_c_nn_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_sca_top_scatteringangle_ttbarframe,
            // gen_c_sca_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_kk_0 + gen_c_rr_0 - gen_c_nn_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_tra_top_scatteringangle_ttbarframe,
            // gen_c_tra_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_kk_0 - gen_c_rr_0 + gen_c_nn_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_kjL_top_scatteringangle_ttbarframe,
            // gen_c_kjL_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_kj_0 - gen_c_rr_0 - gen_c_nn_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_rqL_top_scatteringangle_ttbarframe,
            // gen_c_rqL_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_kk_0 - gen_c_rq_0 - gen_c_nn_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_rkP_top_scatteringangle_ttbarframe,
            // gen_c_rkP_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_rk_0 - gen_c_kr_0 - gen_c_nn_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_rkM_top_scatteringangle_ttbarframe,
            // gen_c_rkM_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_rk_0 + gen_c_kr_0 - gen_c_nn_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_nrP_top_scatteringangle_ttbarframe,
            // gen_c_nrP_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_nr_0 - gen_c_rn_0 - gen_c_kk_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_nrM_top_scatteringangle_ttbarframe,
            // gen_c_nrM_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_nr_0 + gen_c_rn_0 - gen_c_kk_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_nkP_top_scatteringangle_ttbarframe,
            // gen_c_nkP_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_nk_0 - gen_c_kn_0 - gen_c_rr_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_c_nkM_top_scatteringangle_ttbarframe,
            // gen_c_nkM_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_nk_0 + gen_c_kn_0 - gen_c_rr_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);

            // // Lab frame variables
            // fillUnderOverFlow(hvisgen_ll_cHel_top_scatteringangle_ttbarframe,
            // gen_ll_cHel_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_ll_cHel_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_ll_cLab_top_scatteringangle_ttbarframe,
            // gen_ll_cLab_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_ll_cLab_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_ll_kNorm_top_scatteringangle_ttbarframe,
            // gen_ll_kNorm_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_ll_kNorm_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_ll_rNorm_top_scatteringangle_ttbarframe,
            // gen_ll_rNorm_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_ll_rNorm_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_llbar_delta_phi_top_scatteringangle_ttbarframe,
            // gen_llbar_delta_phi_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_llbar_delta_phi_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hvisgen_llbar_delta_eta_top_scatteringangle_ttbarframe,
            // gen_llbar_delta_eta_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_llbar_delta_eta_0,
            // gen_top_scatteringangle_ttbarframe_0}, trueLevelWeight_0,
            // symmetrizeSpinVar);
        }

        Long64_t recoid(-999);

        for (Long64_t jentry = 0; jentry < recoentry.size(); jentry++) {
            if (recoentry[jentry] == entry_0 && gen_l_phi_0 == recogenlphi[jentry]) {
                recoid = jentry;
                n_matched_ID_reco_matched++;
                break;
            } else
                recoid = -999;
        }
        matched_ID_reco.push_back(recoid);
    }

    std::cout << "n_matched_ID_reco: " << matched_ID_reco.size() << std::endl;
    std::cout << "n_matched_ID_reco_matched: " << n_matched_ID_reco_matched << std::endl;

    // ***********************
    // ***********************
    // Fill Migration matrices
    // ***********************
    // ***********************

    for (Long64_t i = 0; i < matched_ID_reco.size(); i++) {
        Long64_t recoid = matched_ID_reco[i];

        if (recoid == -999) {
            // Generated but no reco match

            LoadTree0(i);
            fChain0->GetEntry(i);

            if (!isTopGen_0) continue;

            // Amandeep :
            // To avoid true underflow/overflow
            if (gen_ttbar_mass_0 >= 2000.0) {
                gen_ttbar_mass_0 = 1999.0;
            }

            if (gen_ttbar_mass_0 <= 250.0) {
                gen_ttbar_mass_0 = 251.0;
            }
            if (gen_n_extraJets_iso08_0 >= 3.5) {
                gen_n_extraJets_iso08_0 = 3.0;
            }
            // End

            // ****************************************
            // Fill Migration matrices for 1D variables
            // ****************************************

            // Kinematic
            fillUnderOverFlow(hrecoVsgen_top_pt, gen_top_pt_binning, reco_top_pt_binning,
                              std::vector<double>{gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0, false,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_l_pt, gen_l_pt_binning, reco_l_pt_binning, std::vector<double>{gen_l_pt_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_lbar_pt, gen_lbar_pt_binning, reco_lbar_pt_binning,
                              std::vector<double>{gen_lbar_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0, false,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ttbar_pt, gen_ttbar_pt_binning, reco_ttbar_pt_binning,
                              std::vector<double>{gen_ttbar_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0, false,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ttbar_mass, gen_ttbar_mass_binning, reco_ttbar_mass_binning,
                              std::vector<double>{gen_ttbar_mass_0}, std::vector<double>{0.0}, trueLevelWeight_0, false,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ttbar_delta_phi, gen_ttbar_delta_phi_binning, reco_ttbar_delta_phi_binning,
                              std::vector<double>{gen_ttbar_delta_phi_0}, std::vector<double>{0.0}, trueLevelWeight_0,
                              false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ttbar_delta_eta, gen_ttbar_delta_eta_binning, reco_ttbar_delta_eta_binning,
                              std::vector<double>{gen_ttbar_delta_eta_0}, std::vector<double>{0.0}, trueLevelWeight_0,
                              false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ttbar_rapidity, gen_ttbar_rapidity_binning, reco_ttbar_rapidity_binning,
                              std::vector<double>{gen_ttbar_rapidity_0}, std::vector<double>{0.0}, trueLevelWeight_0,
                              false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_llbar_pt, gen_llbar_pt_binning, reco_llbar_pt_binning,
                              std::vector<double>{gen_llbar_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0, false,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_llbar_mass, gen_llbar_mass_binning, reco_llbar_mass_binning,
                              std::vector<double>{gen_llbar_mass_0}, std::vector<double>{0.0}, trueLevelWeight_0, false,
                              symmetrizeSpinVar);

            // Spin corr

            // Polarizations
            fillUnderOverFlow(hrecoVsgen_b1k, gen_b1k_binning, reco_b1k_binning, std::vector<double>{gen_b1k_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2k, gen_b2k_binning, reco_b2k_binning, std::vector<double>{gen_b2k_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b1j, gen_b1j_binning, reco_b1j_binning, std::vector<double>{gen_b1j_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2j, gen_b2j_binning, reco_b2j_binning, std::vector<double>{gen_b2j_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b1r, gen_b1r_binning, reco_b1r_binning, std::vector<double>{gen_b1r_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2r, gen_b2r_binning, reco_b2r_binning, std::vector<double>{gen_b2r_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b1q, gen_b1q_binning, reco_b1q_binning, std::vector<double>{gen_b1q_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2q, gen_b2q_binning, reco_b2q_binning, std::vector<double>{gen_b2q_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b1n, gen_b1n_binning, reco_b1n_binning, std::vector<double>{gen_b1n_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2n, gen_b2n_binning, reco_b2n_binning, std::vector<double>{gen_b2n_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);

            // Diagonal elements
            fillUnderOverFlow(hrecoVsgen_c_kk, gen_c_kk_binning, reco_c_kk_binning, std::vector<double>{gen_c_kk_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rr, gen_c_rr_binning, reco_c_rr_binning, std::vector<double>{gen_c_rr_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nn, gen_c_nn_binning, reco_c_nn_binning, std::vector<double>{gen_c_nn_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);

            // Off diagonal elements
            fillUnderOverFlow(hrecoVsgen_c_rk, gen_c_rk_binning, reco_c_rk_binning, std::vector<double>{gen_c_rk_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_kr, gen_c_kr_binning, reco_c_kr_binning, std::vector<double>{gen_c_kr_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nr, gen_c_nr_binning, reco_c_nr_binning, std::vector<double>{gen_c_nr_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rn, gen_c_rn_binning, reco_c_rn_binning, std::vector<double>{gen_c_rn_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nk, gen_c_nk_binning, reco_c_nk_binning, std::vector<double>{gen_c_nk_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_kn, gen_c_kn_binning, reco_c_kn_binning, std::vector<double>{gen_c_kn_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);

            fillUnderOverFlow(hrecoVsgen_c_Prk, gen_c_Prk_binning, reco_c_Prk_binning,
                              std::vector<double>{gen_c_rk_0 + gen_c_kr_0}, std::vector<double>{0.0}, trueLevelWeight_0,
                              false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Mrk, gen_c_Mrk_binning, reco_c_Mrk_binning,
                              std::vector<double>{gen_c_rk_0 - gen_c_kr_0}, std::vector<double>{0.0}, trueLevelWeight_0,
                              false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Pnr, gen_c_Pnr_binning, reco_c_Pnr_binning,
                              std::vector<double>{gen_c_nr_0 + gen_c_rn_0}, std::vector<double>{0.0}, trueLevelWeight_0,
                              false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Mnr, gen_c_Mnr_binning, reco_c_Mnr_binning,
                              std::vector<double>{gen_c_nr_0 - gen_c_rn_0}, std::vector<double>{0.0}, trueLevelWeight_0,
                              false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Pnk, gen_c_Pnk_binning, reco_c_Pnk_binning,
                              std::vector<double>{gen_c_nk_0 + gen_c_kn_0}, std::vector<double>{0.0}, trueLevelWeight_0,
                              false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Mnk, gen_c_Mnk_binning, reco_c_Mnk_binning,
                              std::vector<double>{gen_c_nk_0 - gen_c_kn_0}, std::vector<double>{0.0}, trueLevelWeight_0,
                              false, symmetrizeSpinVar);

            // Starred axes
            fillUnderOverFlow(hrecoVsgen_c_kj, gen_c_kj_binning, reco_c_kj_binning, std::vector<double>{gen_c_kj_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rq, gen_c_rq_binning, reco_c_rq_binning, std::vector<double>{gen_c_rq_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rj, gen_c_rj_binning, reco_c_rj_binning, std::vector<double>{gen_c_rj_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_jr, gen_c_jr_binning, reco_c_jr_binning, std::vector<double>{gen_c_jr_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);

            fillUnderOverFlow(hrecoVsgen_c_Prj, gen_c_Prj_binning, reco_c_Prj_binning,
                              std::vector<double>{gen_c_rj_0 + gen_c_jr_0}, std::vector<double>{0.0}, trueLevelWeight_0,
                              false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Mrj, gen_c_Mrj_binning, reco_c_Mrj_binning,
                              std::vector<double>{gen_c_rj_0 - gen_c_jr_0}, std::vector<double>{0.0}, trueLevelWeight_0,
                              false, symmetrizeSpinVar);

            // New variables
            fillUnderOverFlow(hrecoVsgen_c_han, gen_c_han_binning, reco_c_han_binning,
                              std::vector<double>{+gen_c_kk_0 - gen_c_rr_0 - gen_c_nn_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_sca, gen_c_sca_binning, reco_c_sca_binning,
                              std::vector<double>{-gen_c_kk_0 + gen_c_rr_0 - gen_c_nn_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_tra, gen_c_tra_binning, reco_c_tra_binning,
                              std::vector<double>{-gen_c_kk_0 - gen_c_rr_0 + gen_c_nn_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_kjL, gen_c_kjL_binning, reco_c_kjL_binning,
                              std::vector<double>{-gen_c_kj_0 - gen_c_rr_0 - gen_c_nn_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rqL, gen_c_rqL_binning, reco_c_rqL_binning,
                              std::vector<double>{-gen_c_kk_0 - gen_c_rq_0 - gen_c_nn_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rkP, gen_c_rkP_binning, reco_c_rkP_binning,
                              std::vector<double>{-gen_c_rk_0 - gen_c_kr_0 - gen_c_nn_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rkM, gen_c_rkM_binning, reco_c_rkM_binning,
                              std::vector<double>{-gen_c_rk_0 + gen_c_kr_0 - gen_c_nn_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nrP, gen_c_nrP_binning, reco_c_nrP_binning,
                              std::vector<double>{-gen_c_nr_0 - gen_c_rn_0 - gen_c_kk_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nrM, gen_c_nrM_binning, reco_c_nrM_binning,
                              std::vector<double>{-gen_c_nr_0 + gen_c_rn_0 - gen_c_kk_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nkP, gen_c_nkP_binning, reco_c_nkP_binning,
                              std::vector<double>{-gen_c_nk_0 - gen_c_kn_0 - gen_c_rr_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nkM, gen_c_nkM_binning, reco_c_nkM_binning,
                              std::vector<double>{-gen_c_nk_0 + gen_c_kn_0 - gen_c_rr_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);

            // Lab frame variables
            fillUnderOverFlow(hrecoVsgen_ll_cHel, gen_ll_cHel_binning, reco_ll_cHel_binning,
                              std::vector<double>{gen_ll_cHel_0}, std::vector<double>{0.0}, trueLevelWeight_0, false,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ll_cLab, gen_ll_cLab_binning, reco_ll_cLab_binning,
                              std::vector<double>{gen_ll_cLab_0}, std::vector<double>{0.0}, trueLevelWeight_0, false,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ll_kNorm, gen_ll_kNorm_binning, reco_ll_kNorm_binning,
                              std::vector<double>{gen_ll_kNorm_0}, std::vector<double>{0.0}, trueLevelWeight_0, false,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ll_rNorm, gen_ll_rNorm_binning, reco_ll_rNorm_binning,
                              std::vector<double>{gen_ll_rNorm_0}, std::vector<double>{0.0}, trueLevelWeight_0, false,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_llbar_delta_phi, gen_llbar_delta_phi_binning, reco_llbar_delta_phi_binning,
                              std::vector<double>{gen_llbar_delta_phi_0}, std::vector<double>{0.0}, trueLevelWeight_0,
                              false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_llbar_delta_eta, gen_llbar_delta_eta_binning, reco_llbar_delta_eta_binning,
                              std::vector<double>{gen_llbar_delta_eta_0}, std::vector<double>{0.0}, trueLevelWeight_0,
                              false, symmetrizeSpinVar);

            // ****************************************
            // Fill Migration matrices for 2D variables
            // ****************************************

            // Amandeep : Fill 2D migration matrices
            // fillUnderOverFlow(hrecoVsgen_c_kk_mttbar,
            // gen_c_kk_mttbar_binning, reco_c_kk_mttbar_binning,
            // std::vector<double>{gen_c_kk_0, gen_ttbar_mass_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar);

            // ******
            // mttbar
            // ******

            // Polarizations
            fillUnderOverFlow(hrecoVsgen_b1k_exjet, gen_b1k_exjet_binning, reco_b1k_exjet_binning,
                              std::vector<double>{gen_b1k_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2k_exjet, gen_b2k_exjet_binning, reco_b2k_exjet_binning,
                              std::vector<double>{gen_b2k_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b1j_exjet, gen_b1j_exjet_binning, reco_b1j_exjet_binning,
                              std::vector<double>{gen_b1j_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2j_exjet, gen_b2j_exjet_binning, reco_b2j_exjet_binning,
                              std::vector<double>{gen_b2j_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b1r_exjet, gen_b1r_exjet_binning, reco_b1r_exjet_binning,
                              std::vector<double>{gen_b1r_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2r_exjet, gen_b2r_exjet_binning, reco_b2r_exjet_binning,
                              std::vector<double>{gen_b2r_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b1q_exjet, gen_b1q_exjet_binning, reco_b1q_exjet_binning,
                              std::vector<double>{gen_b1q_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2q_exjet, gen_b2q_exjet_binning, reco_b2q_exjet_binning,
                              std::vector<double>{gen_b2q_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b1n_exjet, gen_b1n_exjet_binning, reco_b1n_exjet_binning,
                              std::vector<double>{gen_b1n_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2n_exjet, gen_b2n_exjet_binning, reco_b2n_exjet_binning,
                              std::vector<double>{gen_b2n_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);

            // Diagonal elements
            fillUnderOverFlow(hrecoVsgen_c_kk_exjet, gen_c_kk_exjet_binning, reco_c_kk_exjet_binning,
                              std::vector<double>{gen_c_kk_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rr_exjet, gen_c_rr_exjet_binning, reco_c_rr_exjet_binning,
                              std::vector<double>{gen_c_rr_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nn_exjet, gen_c_nn_exjet_binning, reco_c_nn_exjet_binning,
                              std::vector<double>{gen_c_nn_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);

            // Off diagonal elements
            fillUnderOverFlow(hrecoVsgen_c_rk_exjet, gen_c_rk_exjet_binning, reco_c_rk_exjet_binning,
                              std::vector<double>{gen_c_rk_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_kr_exjet, gen_c_kr_exjet_binning, reco_c_kr_exjet_binning,
                              std::vector<double>{gen_c_kr_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nr_exjet, gen_c_nr_exjet_binning, reco_c_nr_exjet_binning,
                              std::vector<double>{gen_c_nr_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rn_exjet, gen_c_rn_exjet_binning, reco_c_rn_exjet_binning,
                              std::vector<double>{gen_c_rn_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nk_exjet, gen_c_nk_exjet_binning, reco_c_nk_exjet_binning,
                              std::vector<double>{gen_c_nk_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_kn_exjet, gen_c_kn_exjet_binning, reco_c_kn_exjet_binning,
                              std::vector<double>{gen_c_kn_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);

            fillUnderOverFlow(hrecoVsgen_c_Prk_exjet, gen_c_Prk_exjet_binning, reco_c_Prk_exjet_binning,
                              std::vector<double>{gen_c_rk_0 + gen_c_kr_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Mrk_exjet, gen_c_Mrk_exjet_binning, reco_c_Mrk_exjet_binning,
                              std::vector<double>{gen_c_rk_0 - gen_c_kr_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Pnr_exjet, gen_c_Pnr_exjet_binning, reco_c_Pnr_exjet_binning,
                              std::vector<double>{gen_c_nr_0 + gen_c_rn_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Mnr_exjet, gen_c_Mnr_exjet_binning, reco_c_Mnr_exjet_binning,
                              std::vector<double>{gen_c_nr_0 - gen_c_rn_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Pnk_exjet, gen_c_Pnk_exjet_binning, reco_c_Pnk_exjet_binning,
                              std::vector<double>{gen_c_nk_0 + gen_c_kn_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Mnk_exjet, gen_c_Mnk_exjet_binning, reco_c_Mnk_exjet_binning,
                              std::vector<double>{gen_c_nk_0 - gen_c_kn_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);

            // Starred axes
            fillUnderOverFlow(hrecoVsgen_c_kj_exjet, gen_c_kj_exjet_binning, reco_c_kj_exjet_binning,
                              std::vector<double>{gen_c_kj_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rq_exjet, gen_c_rq_exjet_binning, reco_c_rq_exjet_binning,
                              std::vector<double>{gen_c_rq_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rj_exjet, gen_c_rj_exjet_binning, reco_c_rj_exjet_binning,
                              std::vector<double>{gen_c_rj_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_jr_exjet, gen_c_jr_exjet_binning, reco_c_jr_exjet_binning,
                              std::vector<double>{gen_c_jr_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);

            fillUnderOverFlow(hrecoVsgen_c_Prj_exjet, gen_c_Prj_exjet_binning, reco_c_Prj_exjet_binning,
                              std::vector<double>{gen_c_rj_0 + gen_c_jr_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Mrj_exjet, gen_c_Mrj_exjet_binning, reco_c_Mrj_exjet_binning,
                              std::vector<double>{gen_c_rj_0 - gen_c_jr_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);

            // New variables
            fillUnderOverFlow(hrecoVsgen_c_han_exjet, gen_c_han_exjet_binning, reco_c_han_exjet_binning,
                              std::vector<double>{+gen_c_kk_0 - gen_c_rr_0 - gen_c_nn_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_sca_exjet, gen_c_sca_exjet_binning, reco_c_sca_exjet_binning,
                              std::vector<double>{-gen_c_kk_0 + gen_c_rr_0 - gen_c_nn_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_tra_exjet, gen_c_tra_exjet_binning, reco_c_tra_exjet_binning,
                              std::vector<double>{-gen_c_kk_0 - gen_c_rr_0 + gen_c_nn_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_kjL_exjet, gen_c_kjL_exjet_binning, reco_c_kjL_exjet_binning,
                              std::vector<double>{-gen_c_kj_0 - gen_c_rr_0 - gen_c_nn_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rqL_exjet, gen_c_rqL_exjet_binning, reco_c_rqL_exjet_binning,
                              std::vector<double>{-gen_c_kk_0 - gen_c_rq_0 - gen_c_nn_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rkP_exjet, gen_c_rkP_exjet_binning, reco_c_rkP_exjet_binning,
                              std::vector<double>{-gen_c_rk_0 - gen_c_kr_0 - gen_c_nn_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rkM_exjet, gen_c_rkM_exjet_binning, reco_c_rkM_exjet_binning,
                              std::vector<double>{-gen_c_rk_0 + gen_c_kr_0 - gen_c_nn_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nrP_exjet, gen_c_nrP_exjet_binning, reco_c_nrP_exjet_binning,
                              std::vector<double>{-gen_c_nr_0 - gen_c_rn_0 - gen_c_kk_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nrM_exjet, gen_c_nrM_exjet_binning, reco_c_nrM_exjet_binning,
                              std::vector<double>{-gen_c_nr_0 + gen_c_rn_0 - gen_c_kk_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nkP_exjet, gen_c_nkP_exjet_binning, reco_c_nkP_exjet_binning,
                              std::vector<double>{-gen_c_nk_0 - gen_c_kn_0 - gen_c_rr_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nkM_exjet, gen_c_nkM_exjet_binning, reco_c_nkM_exjet_binning,
                              std::vector<double>{-gen_c_nk_0 + gen_c_kn_0 - gen_c_rr_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);

            // Lab frame variables
            fillUnderOverFlow(hrecoVsgen_ll_cHel_exjet, gen_ll_cHel_exjet_binning, reco_ll_cHel_exjet_binning,
                              std::vector<double>{gen_ll_cHel_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ll_cLab_exjet, gen_ll_cLab_exjet_binning, reco_ll_cLab_exjet_binning,
                              std::vector<double>{gen_ll_cLab_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ll_kNorm_exjet, gen_ll_kNorm_exjet_binning, reco_ll_kNorm_exjet_binning,
                              std::vector<double>{gen_ll_kNorm_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ll_rNorm_exjet, gen_ll_rNorm_exjet_binning, reco_ll_rNorm_exjet_binning,
                              std::vector<double>{gen_ll_rNorm_0, gen_n_extraJets_iso08_0}, std::vector<double>{0.0},
                              trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_llbar_delta_phi_exjet, gen_llbar_delta_phi_exjet_binning,
                              reco_llbar_delta_phi_exjet_binning,
                              std::vector<double>{gen_llbar_delta_phi_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_llbar_delta_eta_exjet, gen_llbar_delta_eta_exjet_binning,
                              reco_llbar_delta_eta_exjet_binning,
                              std::vector<double>{gen_llbar_delta_eta_0, gen_n_extraJets_iso08_0},
                              std::vector<double>{0.0}, trueLevelWeight_0, false, symmetrizeSpinVar);

            // ******
            // top_pt
            // ******

            // // Polarizations
            // fillUnderOverFlow(hrecoVsgen_b1k_top_pt, gen_b1k_top_pt_binning,
            // reco_b1k_top_pt_binning, std::vector<double>{gen_b1k_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b2k_top_pt, gen_b2k_top_pt_binning,
            // reco_b2k_top_pt_binning, std::vector<double>{gen_b2k_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b1j_top_pt, gen_b1j_top_pt_binning,
            // reco_b1j_top_pt_binning, std::vector<double>{gen_b1j_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b2j_top_pt, gen_b2j_top_pt_binning,
            // reco_b2j_top_pt_binning, std::vector<double>{gen_b2j_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b1r_top_pt, gen_b1r_top_pt_binning,
            // reco_b1r_top_pt_binning, std::vector<double>{gen_b1r_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b2r_top_pt, gen_b2r_top_pt_binning,
            // reco_b2r_top_pt_binning, std::vector<double>{gen_b2r_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b1q_top_pt, gen_b1q_top_pt_binning,
            // reco_b1q_top_pt_binning, std::vector<double>{gen_b1q_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b2q_top_pt, gen_b2q_top_pt_binning,
            // reco_b2q_top_pt_binning, std::vector<double>{gen_b2q_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b1n_top_pt, gen_b1n_top_pt_binning,
            // reco_b1n_top_pt_binning, std::vector<double>{gen_b1n_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b2n_top_pt, gen_b2n_top_pt_binning,
            // reco_b2n_top_pt_binning, std::vector<double>{gen_b2n_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);

            // // Diagonal elements
            // fillUnderOverFlow(hrecoVsgen_c_kk_top_pt,
            // gen_c_kk_top_pt_binning, reco_c_kk_top_pt_binning,
            // std::vector<double>{gen_c_kk_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_rr_top_pt,
            // gen_c_rr_top_pt_binning, reco_c_rr_top_pt_binning,
            // std::vector<double>{gen_c_rr_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_nn_top_pt,
            // gen_c_nn_top_pt_binning, reco_c_nn_top_pt_binning,
            // std::vector<double>{gen_c_nn_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar);

            // // Off diagonal elements
            // fillUnderOverFlow(hrecoVsgen_c_rk_top_pt,
            // gen_c_rk_top_pt_binning, reco_c_rk_top_pt_binning,
            // std::vector<double>{gen_c_rk_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_kr_top_pt,
            // gen_c_kr_top_pt_binning, reco_c_kr_top_pt_binning,
            // std::vector<double>{gen_c_kr_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_nr_top_pt,
            // gen_c_nr_top_pt_binning, reco_c_nr_top_pt_binning,
            // std::vector<double>{gen_c_nr_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_rn_top_pt,
            // gen_c_rn_top_pt_binning, reco_c_rn_top_pt_binning,
            // std::vector<double>{gen_c_rn_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_nk_top_pt,
            // gen_c_nk_top_pt_binning, reco_c_nk_top_pt_binning,
            // std::vector<double>{gen_c_nk_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_kn_top_pt,
            // gen_c_kn_top_pt_binning, reco_c_kn_top_pt_binning,
            // std::vector<double>{gen_c_kn_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar);

            // fillUnderOverFlow(hrecoVsgen_c_Prk_top_pt,
            // gen_c_Prk_top_pt_binning, reco_c_Prk_top_pt_binning,
            // std::vector<double>{gen_c_rk_0 + gen_c_kr_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_Mrk_top_pt,
            // gen_c_Mrk_top_pt_binning, reco_c_Mrk_top_pt_binning,
            // std::vector<double>{gen_c_rk_0 - gen_c_kr_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_Pnr_top_pt,
            // gen_c_Pnr_top_pt_binning, reco_c_Pnr_top_pt_binning,
            // std::vector<double>{gen_c_nr_0 + gen_c_rn_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_Mnr_top_pt,
            // gen_c_Mnr_top_pt_binning, reco_c_Mnr_top_pt_binning,
            // std::vector<double>{gen_c_nr_0 - gen_c_rn_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_Pnk_top_pt,
            // gen_c_Pnk_top_pt_binning, reco_c_Pnk_top_pt_binning,
            // std::vector<double>{gen_c_nk_0 + gen_c_kn_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_Mnk_top_pt,
            // gen_c_Mnk_top_pt_binning, reco_c_Mnk_top_pt_binning,
            // std::vector<double>{gen_c_nk_0 - gen_c_kn_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar);

            // // Starred axes
            // fillUnderOverFlow(hrecoVsgen_c_kj_top_pt,
            // gen_c_kj_top_pt_binning, reco_c_kj_top_pt_binning,
            // std::vector<double>{gen_c_kj_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_rq_top_pt,
            // gen_c_rq_top_pt_binning, reco_c_rq_top_pt_binning,
            // std::vector<double>{gen_c_rq_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_rj_top_pt,
            // gen_c_rj_top_pt_binning, reco_c_rj_top_pt_binning,
            // std::vector<double>{gen_c_rj_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_jr_top_pt,
            // gen_c_jr_top_pt_binning, reco_c_jr_top_pt_binning,
            // std::vector<double>{gen_c_jr_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar);

            // fillUnderOverFlow(hrecoVsgen_c_Prj_top_pt,
            // gen_c_Prj_top_pt_binning, reco_c_Prj_top_pt_binning,
            // std::vector<double>{gen_c_rj_0 + gen_c_jr_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_Mrj_top_pt,
            // gen_c_Mrj_top_pt_binning, reco_c_Mrj_top_pt_binning,
            // std::vector<double>{gen_c_rj_0 - gen_c_jr_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar);

            // // New variables
            // fillUnderOverFlow(hrecoVsgen_c_han_top_pt,
            // gen_c_han_top_pt_binning, reco_c_han_top_pt_binning,
            // std::vector<double>{+gen_c_kk_0 - gen_c_rr_0 - gen_c_nn_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_sca_top_pt,
            // gen_c_sca_top_pt_binning, reco_c_sca_top_pt_binning,
            // std::vector<double>{-gen_c_kk_0 + gen_c_rr_0 - gen_c_nn_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_tra_top_pt,
            // gen_c_tra_top_pt_binning, reco_c_tra_top_pt_binning,
            // std::vector<double>{-gen_c_kk_0 - gen_c_rr_0 + gen_c_nn_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_kjL_top_pt,
            // gen_c_kjL_top_pt_binning, reco_c_kjL_top_pt_binning,
            // std::vector<double>{-gen_c_kj_0 - gen_c_rr_0 - gen_c_nn_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rqL_top_pt,
            // gen_c_rqL_top_pt_binning, reco_c_rqL_top_pt_binning,
            // std::vector<double>{-gen_c_kk_0 - gen_c_rq_0 - gen_c_nn_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rkP_top_pt,
            // gen_c_rkP_top_pt_binning, reco_c_rkP_top_pt_binning,
            // std::vector<double>{-gen_c_rk_0 - gen_c_kr_0 - gen_c_nn_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rkM_top_pt,
            // gen_c_rkM_top_pt_binning, reco_c_rkM_top_pt_binning,
            // std::vector<double>{-gen_c_rk_0 + gen_c_kr_0 - gen_c_nn_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nrP_top_pt,
            // gen_c_nrP_top_pt_binning, reco_c_nrP_top_pt_binning,
            // std::vector<double>{-gen_c_nr_0 - gen_c_rn_0 - gen_c_kk_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nrM_top_pt,
            // gen_c_nrM_top_pt_binning, reco_c_nrM_top_pt_binning,
            // std::vector<double>{-gen_c_nr_0 + gen_c_rn_0 - gen_c_kk_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nkP_top_pt,
            // gen_c_nkP_top_pt_binning, reco_c_nkP_top_pt_binning,
            // std::vector<double>{-gen_c_nk_0 - gen_c_kn_0 - gen_c_rr_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nkM_top_pt,
            // gen_c_nkM_top_pt_binning, reco_c_nkM_top_pt_binning,
            // std::vector<double>{-gen_c_nk_0 + gen_c_kn_0 - gen_c_rr_0,
            // gen_top_pt_0}, std::vector<double>{0.0}, trueLevelWeight_0,
            // false, symmetrizeSpinVar);

            // // Lab frame variables
            // fillUnderOverFlow(hrecoVsgen_ll_cHel_top_pt,
            // gen_ll_cHel_top_pt_binning, reco_ll_cHel_top_pt_binning,
            // std::vector<double>{gen_ll_cHel_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_ll_cLab_top_pt,
            // gen_ll_cLab_top_pt_binning, reco_ll_cLab_top_pt_binning,
            // std::vector<double>{gen_ll_cLab_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_ll_kNorm_top_pt,
            // gen_ll_kNorm_top_pt_binning, reco_ll_kNorm_top_pt_binning,
            // std::vector<double>{gen_ll_kNorm_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_ll_rNorm_top_pt,
            // gen_ll_rNorm_top_pt_binning, reco_ll_rNorm_top_pt_binning,
            // std::vector<double>{gen_ll_rNorm_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_llbar_delta_phi_top_pt,
            // gen_llbar_delta_phi_top_pt_binning,
            // reco_llbar_delta_phi_top_pt_binning,
            // std::vector<double>{gen_llbar_delta_phi_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_llbar_delta_eta_top_pt,
            // gen_llbar_delta_eta_top_pt_binning,
            // reco_llbar_delta_eta_top_pt_binning,
            // std::vector<double>{gen_llbar_delta_eta_0, gen_top_pt_0},
            // std::vector<double>{0.0}, trueLevelWeight_0, false,
            // symmetrizeSpinVar);

            // ******************************
            // top_scatteringangle_ttbarframe
            // ******************************

            // // Polarizations
            // fillUnderOverFlow(hrecoVsgen_b1k_top_scatteringangle_ttbarframe,
            // gen_b1k_top_scatteringangle_ttbarframe_binning,
            // reco_b1k_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1k_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b2k_top_scatteringangle_ttbarframe,
            // gen_b2k_top_scatteringangle_ttbarframe_binning,
            // reco_b2k_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2k_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b1j_top_scatteringangle_ttbarframe,
            // gen_b1j_top_scatteringangle_ttbarframe_binning,
            // reco_b1j_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1j_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b2j_top_scatteringangle_ttbarframe,
            // gen_b2j_top_scatteringangle_ttbarframe_binning,
            // reco_b2j_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2j_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b1r_top_scatteringangle_ttbarframe,
            // gen_b1r_top_scatteringangle_ttbarframe_binning,
            // reco_b1r_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1r_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b2r_top_scatteringangle_ttbarframe,
            // gen_b2r_top_scatteringangle_ttbarframe_binning,
            // reco_b2r_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2r_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b1q_top_scatteringangle_ttbarframe,
            // gen_b1q_top_scatteringangle_ttbarframe_binning,
            // reco_b1q_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1q_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b2q_top_scatteringangle_ttbarframe,
            // gen_b2q_top_scatteringangle_ttbarframe_binning,
            // reco_b2q_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2q_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b1n_top_scatteringangle_ttbarframe,
            // gen_b1n_top_scatteringangle_ttbarframe_binning,
            // reco_b1n_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1n_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b2n_top_scatteringangle_ttbarframe,
            // gen_b2n_top_scatteringangle_ttbarframe_binning,
            // reco_b2n_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2n_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);

            // // Diagonal elements
            // fillUnderOverFlow(hrecoVsgen_c_kk_top_scatteringangle_ttbarframe,
            // gen_c_kk_top_scatteringangle_ttbarframe_binning,
            // reco_c_kk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_kk_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rr_top_scatteringangle_ttbarframe,
            // gen_c_rr_top_scatteringangle_ttbarframe_binning,
            // reco_c_rr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rr_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nn_top_scatteringangle_ttbarframe,
            // gen_c_nn_top_scatteringangle_ttbarframe_binning,
            // reco_c_nn_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nn_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);

            // // Off diagonal elements
            // fillUnderOverFlow(hrecoVsgen_c_rk_top_scatteringangle_ttbarframe,
            // gen_c_rk_top_scatteringangle_ttbarframe_binning,
            // reco_c_rk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rk_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_kr_top_scatteringangle_ttbarframe,
            // gen_c_kr_top_scatteringangle_ttbarframe_binning,
            // reco_c_kr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_kr_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nr_top_scatteringangle_ttbarframe,
            // gen_c_nr_top_scatteringangle_ttbarframe_binning,
            // reco_c_nr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nr_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rn_top_scatteringangle_ttbarframe,
            // gen_c_rn_top_scatteringangle_ttbarframe_binning,
            // reco_c_rn_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rn_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nk_top_scatteringangle_ttbarframe,
            // gen_c_nk_top_scatteringangle_ttbarframe_binning,
            // reco_c_nk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nk_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_kn_top_scatteringangle_ttbarframe,
            // gen_c_kn_top_scatteringangle_ttbarframe_binning,
            // reco_c_kn_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_kn_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);

            // fillUnderOverFlow(hrecoVsgen_c_Prk_top_scatteringangle_ttbarframe,
            // gen_c_Prk_top_scatteringangle_ttbarframe_binning,
            // reco_c_Prk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rk_0 + gen_c_kr_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_Mrk_top_scatteringangle_ttbarframe,
            // gen_c_Mrk_top_scatteringangle_ttbarframe_binning,
            // reco_c_Mrk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rk_0 - gen_c_kr_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_Pnr_top_scatteringangle_ttbarframe,
            // gen_c_Pnr_top_scatteringangle_ttbarframe_binning,
            // reco_c_Pnr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nr_0 + gen_c_rn_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_Mnr_top_scatteringangle_ttbarframe,
            // gen_c_Mnr_top_scatteringangle_ttbarframe_binning,
            // reco_c_Mnr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nr_0 - gen_c_rn_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_Pnk_top_scatteringangle_ttbarframe,
            // gen_c_Pnk_top_scatteringangle_ttbarframe_binning,
            // reco_c_Pnk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nk_0 + gen_c_kn_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_Mnk_top_scatteringangle_ttbarframe,
            // gen_c_Mnk_top_scatteringangle_ttbarframe_binning,
            // reco_c_Mnk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nk_0 - gen_c_kn_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);

            // // Starred axes
            // fillUnderOverFlow(hrecoVsgen_c_kj_top_scatteringangle_ttbarframe,
            // gen_c_kj_top_scatteringangle_ttbarframe_binning,
            // reco_c_kj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_kj_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rq_top_scatteringangle_ttbarframe,
            // gen_c_rq_top_scatteringangle_ttbarframe_binning,
            // reco_c_rq_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rq_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rj_top_scatteringangle_ttbarframe,
            // gen_c_rj_top_scatteringangle_ttbarframe_binning,
            // reco_c_rj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rj_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_jr_top_scatteringangle_ttbarframe,
            // gen_c_jr_top_scatteringangle_ttbarframe_binning,
            // reco_c_jr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_jr_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);

            // fillUnderOverFlow(hrecoVsgen_c_Prj_top_scatteringangle_ttbarframe,
            // gen_c_Prj_top_scatteringangle_ttbarframe_binning,
            // reco_c_Prj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rj_0 + gen_c_jr_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_Mrj_top_scatteringangle_ttbarframe,
            // gen_c_Mrj_top_scatteringangle_ttbarframe_binning,
            // reco_c_Mrj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rj_0 - gen_c_jr_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);

            // // New variables
            // fillUnderOverFlow(hrecoVsgen_c_han_top_scatteringangle_ttbarframe,
            // gen_c_han_top_scatteringangle_ttbarframe_binning,
            // reco_c_han_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{+gen_c_kk_0 - gen_c_rr_0 - gen_c_nn_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_sca_top_scatteringangle_ttbarframe,
            // gen_c_sca_top_scatteringangle_ttbarframe_binning,
            // reco_c_sca_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_kk_0 + gen_c_rr_0 - gen_c_nn_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_tra_top_scatteringangle_ttbarframe,
            // gen_c_tra_top_scatteringangle_ttbarframe_binning,
            // reco_c_tra_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_kk_0 - gen_c_rr_0 + gen_c_nn_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_kjL_top_scatteringangle_ttbarframe,
            // gen_c_kjL_top_scatteringangle_ttbarframe_binning,
            // reco_c_kjL_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_kj_0 - gen_c_rr_0 - gen_c_nn_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rqL_top_scatteringangle_ttbarframe,
            // gen_c_rqL_top_scatteringangle_ttbarframe_binning,
            // reco_c_rqL_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_kk_0 - gen_c_rq_0 - gen_c_nn_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rkP_top_scatteringangle_ttbarframe,
            // gen_c_rkP_top_scatteringangle_ttbarframe_binning,
            // reco_c_rkP_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_rk_0 - gen_c_kr_0 - gen_c_nn_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rkM_top_scatteringangle_ttbarframe,
            // gen_c_rkM_top_scatteringangle_ttbarframe_binning,
            // reco_c_rkM_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_rk_0 + gen_c_kr_0 - gen_c_nn_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nrP_top_scatteringangle_ttbarframe,
            // gen_c_nrP_top_scatteringangle_ttbarframe_binning,
            // reco_c_nrP_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_nr_0 - gen_c_rn_0 - gen_c_kk_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nrM_top_scatteringangle_ttbarframe,
            // gen_c_nrM_top_scatteringangle_ttbarframe_binning,
            // reco_c_nrM_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_nr_0 + gen_c_rn_0 - gen_c_kk_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nkP_top_scatteringangle_ttbarframe,
            // gen_c_nkP_top_scatteringangle_ttbarframe_binning,
            // reco_c_nkP_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_nk_0 - gen_c_kn_0 - gen_c_rr_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nkM_top_scatteringangle_ttbarframe,
            // gen_c_nkM_top_scatteringangle_ttbarframe_binning,
            // reco_c_nkM_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_nk_0 + gen_c_kn_0 - gen_c_rr_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);

            // // Lab frame variables
            // fillUnderOverFlow(hrecoVsgen_ll_cHel_top_scatteringangle_ttbarframe,
            // gen_ll_cHel_top_scatteringangle_ttbarframe_binning,
            // reco_ll_cHel_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_ll_cHel_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_ll_cLab_top_scatteringangle_ttbarframe,
            // gen_ll_cLab_top_scatteringangle_ttbarframe_binning,
            // reco_ll_cLab_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_ll_cLab_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_ll_kNorm_top_scatteringangle_ttbarframe,
            // gen_ll_kNorm_top_scatteringangle_ttbarframe_binning,
            // reco_ll_kNorm_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_ll_kNorm_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_ll_rNorm_top_scatteringangle_ttbarframe,
            // gen_ll_rNorm_top_scatteringangle_ttbarframe_binning,
            // reco_ll_rNorm_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_ll_rNorm_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_llbar_delta_phi_top_scatteringangle_ttbarframe,
            // gen_llbar_delta_phi_top_scatteringangle_ttbarframe_binning,
            // reco_llbar_delta_phi_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_llbar_delta_phi_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_llbar_delta_eta_top_scatteringangle_ttbarframe,
            // gen_llbar_delta_eta_top_scatteringangle_ttbarframe_binning,
            // reco_llbar_delta_eta_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_llbar_delta_eta_0,
            // gen_top_scatteringangle_ttbarframe_0}, std::vector<double>{0.0},
            // trueLevelWeight_0, false, symmetrizeSpinVar);
        }
    }

    for (Long64_t i = 0; i < matched_ID_reco.size(); i++) {
        Long64_t recoid = matched_ID_reco[i];

        // If there is a matching reco to the gen
        if (recoid > -999) {
            LoadTree0(i);
            fChain0->GetEntry(i);

            if (!isTopGen_0) continue;

            LoadTree(recoid);
            fChain->GetEntry(recoid);

            // Amandeep :
            // To avoid true underflow overflow
            if (gen_ttbar_mass >= 2000.0) {
                gen_ttbar_mass = 1999.0;
            }
            if (ttbar_mass >= 2000.0) {
                ttbar_mass = 1999.0;
            }
            if (gen_ttbar_mass <= 250.0) {
                gen_ttbar_mass = 251.0;
            }
            if (ttbar_mass <= 250.0) {
                ttbar_mass = 251.0;
            }
            if (gen_top_pt >= 550.0) {
                gen_top_pt = 549.0;
            }
            if (top_pt >= 550.0) {
                top_pt = 549.0;
            }
            if (gen_n_extraJets_iso08 >= 3.5) {
                gen_n_extraJets_iso08 = 3.0;
            }
            if (n_extraJets_iso08 >= 3.5) {
                n_extraJets_iso08 = 3.0;
            }
            // End

            // ****************************************
            // Fill Migration matrices for 1D variables
            // ****************************************

            // Kinematic
            fillUnderOverFlow(hrecoVsgen_top_pt, gen_top_pt_binning, reco_top_pt_binning,
                              std::vector<double>{gen_top_pt}, std::vector<double>{top_pt}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_l_pt, gen_l_pt_binning, reco_l_pt_binning, std::vector<double>{gen_l_pt},
                              std::vector<double>{l_pt}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_lbar_pt, gen_lbar_pt_binning, reco_lbar_pt_binning,
                              std::vector<double>{gen_lbar_pt}, std::vector<double>{lbar_pt}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ttbar_pt, gen_ttbar_pt_binning, reco_ttbar_pt_binning,
                              std::vector<double>{gen_ttbar_pt}, std::vector<double>{ttbar_pt}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ttbar_mass, gen_ttbar_mass_binning, reco_ttbar_mass_binning,
                              std::vector<double>{gen_ttbar_mass}, std::vector<double>{ttbar_mass}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ttbar_delta_phi, gen_ttbar_delta_phi_binning, reco_ttbar_delta_phi_binning,
                              std::vector<double>{gen_ttbar_delta_phi}, std::vector<double>{ttbar_delta_phi},
                              eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ttbar_delta_eta, gen_ttbar_delta_eta_binning, reco_ttbar_delta_eta_binning,
                              std::vector<double>{gen_ttbar_delta_eta}, std::vector<double>{ttbar_delta_eta},
                              eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ttbar_rapidity, gen_ttbar_rapidity_binning, reco_ttbar_rapidity_binning,
                              std::vector<double>{gen_ttbar_rapidity}, std::vector<double>{ttbar_rapidity}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_llbar_pt, gen_llbar_pt_binning, reco_llbar_pt_binning,
                              std::vector<double>{gen_llbar_pt}, std::vector<double>{llbar_pt}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_llbar_mass, gen_llbar_mass_binning, reco_llbar_mass_binning,
                              std::vector<double>{gen_llbar_mass}, std::vector<double>{llbar_mass}, eventWeight, true,
                              symmetrizeSpinVar);

            // Spin corr

            // Polarizations
            fillUnderOverFlow(hrecoVsgen_b1k, gen_b1k_binning, reco_b1k_binning, std::vector<double>{gen_b1k},
                              std::vector<double>{b1k}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2k, gen_b2k_binning, reco_b2k_binning, std::vector<double>{gen_b2k},
                              std::vector<double>{b2k}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b1j, gen_b1j_binning, reco_b1j_binning, std::vector<double>{gen_b1j},
                              std::vector<double>{b1j}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2j, gen_b2j_binning, reco_b2j_binning, std::vector<double>{gen_b2j},
                              std::vector<double>{b2j}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b1r, gen_b1r_binning, reco_b1r_binning, std::vector<double>{gen_b1r},
                              std::vector<double>{b1r}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2r, gen_b2r_binning, reco_b2r_binning, std::vector<double>{gen_b2r},
                              std::vector<double>{b2r}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b1q, gen_b1q_binning, reco_b1q_binning, std::vector<double>{gen_b1q},
                              std::vector<double>{b1q}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2q, gen_b2q_binning, reco_b2q_binning, std::vector<double>{gen_b2q},
                              std::vector<double>{b2q}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b1n, gen_b1n_binning, reco_b1n_binning, std::vector<double>{gen_b1n},
                              std::vector<double>{b1n}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2n, gen_b2n_binning, reco_b2n_binning, std::vector<double>{gen_b2n},
                              std::vector<double>{b2n}, eventWeight, true, symmetrizeSpinVar);

            // Diagonal elements
            fillUnderOverFlow(hrecoVsgen_c_kk, gen_c_kk_binning, reco_c_kk_binning, std::vector<double>{gen_c_kk},
                              std::vector<double>{c_kk}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rr, gen_c_rr_binning, reco_c_rr_binning, std::vector<double>{gen_c_rr},
                              std::vector<double>{c_rr}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nn, gen_c_nn_binning, reco_c_nn_binning, std::vector<double>{gen_c_nn},
                              std::vector<double>{c_nn}, eventWeight, true, symmetrizeSpinVar);

            // Off diagonal elements
            fillUnderOverFlow(hrecoVsgen_c_rk, gen_c_rk_binning, reco_c_rk_binning, std::vector<double>{gen_c_rk},
                              std::vector<double>{c_rk}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_kr, gen_c_kr_binning, reco_c_kr_binning, std::vector<double>{gen_c_kr},
                              std::vector<double>{c_kr}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nr, gen_c_nr_binning, reco_c_nr_binning, std::vector<double>{gen_c_nr},
                              std::vector<double>{c_nr}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rn, gen_c_rn_binning, reco_c_rn_binning, std::vector<double>{gen_c_rn},
                              std::vector<double>{c_rn}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nk, gen_c_nk_binning, reco_c_nk_binning, std::vector<double>{gen_c_nk},
                              std::vector<double>{c_nk}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_kn, gen_c_kn_binning, reco_c_kn_binning, std::vector<double>{gen_c_kn},
                              std::vector<double>{c_kn}, eventWeight, true, symmetrizeSpinVar);

            fillUnderOverFlow(hrecoVsgen_c_Prk, gen_c_Prk_binning, reco_c_Prk_binning,
                              std::vector<double>{gen_c_rk + gen_c_kr}, std::vector<double>{c_rk + c_kr}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Mrk, gen_c_Mrk_binning, reco_c_Mrk_binning,
                              std::vector<double>{gen_c_rk - gen_c_kr}, std::vector<double>{c_rk - c_kr}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Pnr, gen_c_Pnr_binning, reco_c_Pnr_binning,
                              std::vector<double>{gen_c_nr + gen_c_rn}, std::vector<double>{c_nr + c_rn}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Mnr, gen_c_Mnr_binning, reco_c_Mnr_binning,
                              std::vector<double>{gen_c_nr - gen_c_rn}, std::vector<double>{c_nr - c_rn}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Pnk, gen_c_Pnk_binning, reco_c_Pnk_binning,
                              std::vector<double>{gen_c_nk + gen_c_kn}, std::vector<double>{c_nk + c_kn}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Mnk, gen_c_Mnk_binning, reco_c_Mnk_binning,
                              std::vector<double>{gen_c_nk - gen_c_kn}, std::vector<double>{c_nk - c_kn}, eventWeight,
                              true, symmetrizeSpinVar);

            // Starred axes
            fillUnderOverFlow(hrecoVsgen_c_kj, gen_c_kj_binning, reco_c_kj_binning, std::vector<double>{gen_c_kj},
                              std::vector<double>{c_kj}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rq, gen_c_rq_binning, reco_c_rq_binning, std::vector<double>{gen_c_rq},
                              std::vector<double>{c_rq}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rj, gen_c_rj_binning, reco_c_rj_binning, std::vector<double>{gen_c_rj},
                              std::vector<double>{c_rj}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_jr, gen_c_jr_binning, reco_c_jr_binning, std::vector<double>{gen_c_jr},
                              std::vector<double>{c_jr}, eventWeight, true, symmetrizeSpinVar);

            fillUnderOverFlow(hrecoVsgen_c_Prj, gen_c_Prj_binning, reco_c_Prj_binning,
                              std::vector<double>{gen_c_rj + gen_c_jr}, std::vector<double>{c_rj + c_jr}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Mrj, gen_c_Mrj_binning, reco_c_Mrj_binning,
                              std::vector<double>{gen_c_rj - gen_c_jr}, std::vector<double>{c_rj - c_jr}, eventWeight,
                              true, symmetrizeSpinVar);

            // New variables
            fillUnderOverFlow(hrecoVsgen_c_han, gen_c_han_binning, reco_c_han_binning,
                              std::vector<double>{+gen_c_kk - gen_c_rr - gen_c_nn},
                              std::vector<double>{+c_kk - c_rr - c_nn}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_sca, gen_c_sca_binning, reco_c_sca_binning,
                              std::vector<double>{-gen_c_kk + gen_c_rr - gen_c_nn},
                              std::vector<double>{-c_kk + c_rr - c_nn}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_tra, gen_c_tra_binning, reco_c_tra_binning,
                              std::vector<double>{-gen_c_kk - gen_c_rr + gen_c_nn},
                              std::vector<double>{-c_kk - c_rr + c_nn}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_kjL, gen_c_kjL_binning, reco_c_kjL_binning,
                              std::vector<double>{-gen_c_kj - gen_c_rr - gen_c_nn},
                              std::vector<double>{-c_kj - c_rr - c_nn}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rqL, gen_c_rqL_binning, reco_c_rqL_binning,
                              std::vector<double>{-gen_c_kk - gen_c_rq - gen_c_nn},
                              std::vector<double>{-c_kk - c_rq - c_nn}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rkP, gen_c_rkP_binning, reco_c_rkP_binning,
                              std::vector<double>{-gen_c_rk - gen_c_kr - gen_c_nn},
                              std::vector<double>{-c_rk - c_kr - c_nn}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rkM, gen_c_rkM_binning, reco_c_rkM_binning,
                              std::vector<double>{-gen_c_rk + gen_c_kr - gen_c_nn},
                              std::vector<double>{-c_rk + c_kr - c_nn}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nrP, gen_c_nrP_binning, reco_c_nrP_binning,
                              std::vector<double>{-gen_c_nr - gen_c_rn - gen_c_kk},
                              std::vector<double>{-c_nr - c_rn - c_kk}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nrM, gen_c_nrM_binning, reco_c_nrM_binning,
                              std::vector<double>{-gen_c_nr + gen_c_rn - gen_c_kk},
                              std::vector<double>{-c_nr + c_rn - c_kk}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nkP, gen_c_nkP_binning, reco_c_nkP_binning,
                              std::vector<double>{-gen_c_nk - gen_c_kn - gen_c_rr},
                              std::vector<double>{-c_nk - c_kn - c_rr}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nkM, gen_c_nkM_binning, reco_c_nkM_binning,
                              std::vector<double>{-gen_c_nk + gen_c_kn - gen_c_rr},
                              std::vector<double>{-c_nk + c_kn - c_rr}, eventWeight, true, symmetrizeSpinVar);

            // Lab fram variables
            fillUnderOverFlow(hrecoVsgen_ll_cHel, gen_ll_cHel_binning, reco_ll_cHel_binning,
                              std::vector<double>{gen_ll_cHel}, std::vector<double>{ll_cHel}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ll_cLab, gen_ll_cLab_binning, reco_ll_cLab_binning,
                              std::vector<double>{gen_ll_cLab}, std::vector<double>{ll_cLab}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ll_kNorm, gen_ll_kNorm_binning, reco_ll_kNorm_binning,
                              std::vector<double>{gen_ll_kNorm}, std::vector<double>{ll_kNorm}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ll_rNorm, gen_ll_rNorm_binning, reco_ll_rNorm_binning,
                              std::vector<double>{gen_ll_rNorm}, std::vector<double>{ll_rNorm}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_llbar_delta_phi, gen_llbar_delta_phi_binning, reco_llbar_delta_phi_binning,
                              std::vector<double>{gen_llbar_delta_phi}, std::vector<double>{llbar_delta_phi},
                              eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_llbar_delta_eta, gen_llbar_delta_eta_binning, reco_llbar_delta_eta_binning,
                              std::vector<double>{gen_llbar_delta_eta}, std::vector<double>{llbar_delta_eta},
                              eventWeight, true, symmetrizeSpinVar);

            // ****************************************
            // Fill Migration matrices for 2D variables
            // ****************************************

            // Amandeep : 2D migration matrix
            // fillUnderOverFlow(hrecoVsgen_c_kk_mttbar,
            // gen_c_kk_mttbar_binning, reco_c_kk_mttbar_binning,
            // std::vector<double>{gen_c_kk, gen_ttbar_mass},
            // std::vector<double>{c_kk, ttbar_mass}, eventWeight, true,
            // symmetrizeSpinVar);

            // ******
            // mttbar
            // ******

            // Polarizations
            fillUnderOverFlow(hrecoVsgen_b1k_exjet, gen_b1k_exjet_binning, reco_b1k_exjet_binning,
                              std::vector<double>{gen_b1k, gen_n_extraJets_iso08},
                              std::vector<double>{b1k, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2k_exjet, gen_b2k_exjet_binning, reco_b2k_exjet_binning,
                              std::vector<double>{gen_b2k, gen_n_extraJets_iso08},
                              std::vector<double>{b2k, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b1j_exjet, gen_b1j_exjet_binning, reco_b1j_exjet_binning,
                              std::vector<double>{gen_b1j, gen_n_extraJets_iso08},
                              std::vector<double>{b1j, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2j_exjet, gen_b2j_exjet_binning, reco_b2j_exjet_binning,
                              std::vector<double>{gen_b2j, gen_n_extraJets_iso08},
                              std::vector<double>{b2j, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b1r_exjet, gen_b1r_exjet_binning, reco_b1r_exjet_binning,
                              std::vector<double>{gen_b1r, gen_n_extraJets_iso08},
                              std::vector<double>{b1r, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2r_exjet, gen_b2r_exjet_binning, reco_b2r_exjet_binning,
                              std::vector<double>{gen_b2r, gen_n_extraJets_iso08},
                              std::vector<double>{b2r, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b1q_exjet, gen_b1q_exjet_binning, reco_b1q_exjet_binning,
                              std::vector<double>{gen_b1q, gen_n_extraJets_iso08},
                              std::vector<double>{b1q, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2q_exjet, gen_b2q_exjet_binning, reco_b2q_exjet_binning,
                              std::vector<double>{gen_b2q, gen_n_extraJets_iso08},
                              std::vector<double>{b2q, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b1n_exjet, gen_b1n_exjet_binning, reco_b1n_exjet_binning,
                              std::vector<double>{gen_b1n, gen_n_extraJets_iso08},
                              std::vector<double>{b1n, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_b2n_exjet, gen_b2n_exjet_binning, reco_b2n_exjet_binning,
                              std::vector<double>{gen_b2n, gen_n_extraJets_iso08},
                              std::vector<double>{b2n, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);

            // Diagonal elements
            fillUnderOverFlow(hrecoVsgen_c_kk_exjet, gen_c_kk_exjet_binning, reco_c_kk_exjet_binning,
                              std::vector<double>{gen_c_kk, gen_n_extraJets_iso08},
                              std::vector<double>{c_kk, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rr_exjet, gen_c_rr_exjet_binning, reco_c_rr_exjet_binning,
                              std::vector<double>{gen_c_rr, gen_n_extraJets_iso08},
                              std::vector<double>{c_rr, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nn_exjet, gen_c_nn_exjet_binning, reco_c_nn_exjet_binning,
                              std::vector<double>{gen_c_nn, gen_n_extraJets_iso08},
                              std::vector<double>{c_nn, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);

            // Off diagonal elements
            fillUnderOverFlow(hrecoVsgen_c_rk_exjet, gen_c_rk_exjet_binning, reco_c_rk_exjet_binning,
                              std::vector<double>{gen_c_rk, gen_n_extraJets_iso08},
                              std::vector<double>{c_rk, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_kr_exjet, gen_c_kr_exjet_binning, reco_c_kr_exjet_binning,
                              std::vector<double>{gen_c_kr, gen_n_extraJets_iso08},
                              std::vector<double>{c_kr, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nr_exjet, gen_c_nr_exjet_binning, reco_c_nr_exjet_binning,
                              std::vector<double>{gen_c_nr, gen_n_extraJets_iso08},
                              std::vector<double>{c_nr, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rn_exjet, gen_c_rn_exjet_binning, reco_c_rn_exjet_binning,
                              std::vector<double>{gen_c_rn, gen_n_extraJets_iso08},
                              std::vector<double>{c_rn, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nk_exjet, gen_c_nk_exjet_binning, reco_c_nk_exjet_binning,
                              std::vector<double>{gen_c_nk, gen_n_extraJets_iso08},
                              std::vector<double>{c_nk, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_kn_exjet, gen_c_kn_exjet_binning, reco_c_kn_exjet_binning,
                              std::vector<double>{gen_c_kn, gen_n_extraJets_iso08},
                              std::vector<double>{c_kn, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);

            fillUnderOverFlow(hrecoVsgen_c_Prk_exjet, gen_c_Prk_exjet_binning, reco_c_Prk_exjet_binning,
                              std::vector<double>{gen_c_rk + gen_c_kr, gen_n_extraJets_iso08},
                              std::vector<double>{c_rk + c_kr, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Mrk_exjet, gen_c_Mrk_exjet_binning, reco_c_Mrk_exjet_binning,
                              std::vector<double>{gen_c_rk - gen_c_kr, gen_n_extraJets_iso08},
                              std::vector<double>{c_rk - c_kr, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Pnr_exjet, gen_c_Pnr_exjet_binning, reco_c_Pnr_exjet_binning,
                              std::vector<double>{gen_c_nr + gen_c_rn, gen_n_extraJets_iso08},
                              std::vector<double>{c_nr + c_rn, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Mnr_exjet, gen_c_Mnr_exjet_binning, reco_c_Mnr_exjet_binning,
                              std::vector<double>{gen_c_nr - gen_c_rn, gen_n_extraJets_iso08},
                              std::vector<double>{c_nr - c_rn, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Pnk_exjet, gen_c_Pnk_exjet_binning, reco_c_Pnk_exjet_binning,
                              std::vector<double>{gen_c_nk + gen_c_kn, gen_n_extraJets_iso08},
                              std::vector<double>{c_nk + c_kn, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Mnk_exjet, gen_c_Mnk_exjet_binning, reco_c_Mnk_exjet_binning,
                              std::vector<double>{gen_c_nk - gen_c_kn, gen_n_extraJets_iso08},
                              std::vector<double>{c_nk - c_kn, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);

            // Starred axes
            fillUnderOverFlow(hrecoVsgen_c_kj_exjet, gen_c_kj_exjet_binning, reco_c_kj_exjet_binning,
                              std::vector<double>{gen_c_kj, gen_n_extraJets_iso08},
                              std::vector<double>{c_kj, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rq_exjet, gen_c_rq_exjet_binning, reco_c_rq_exjet_binning,
                              std::vector<double>{gen_c_rq, gen_n_extraJets_iso08},
                              std::vector<double>{c_rq, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rj_exjet, gen_c_rj_exjet_binning, reco_c_rj_exjet_binning,
                              std::vector<double>{gen_c_rj, gen_n_extraJets_iso08},
                              std::vector<double>{c_rj, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_jr_exjet, gen_c_jr_exjet_binning, reco_c_jr_exjet_binning,
                              std::vector<double>{gen_c_jr, gen_n_extraJets_iso08},
                              std::vector<double>{c_jr, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);

            fillUnderOverFlow(hrecoVsgen_c_Prj_exjet, gen_c_Prj_exjet_binning, reco_c_Prj_exjet_binning,
                              std::vector<double>{gen_c_rj + gen_c_jr, gen_n_extraJets_iso08},
                              std::vector<double>{c_rj + c_jr, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_Mrj_exjet, gen_c_Mrj_exjet_binning, reco_c_Mrj_exjet_binning,
                              std::vector<double>{gen_c_rj - gen_c_jr, gen_n_extraJets_iso08},
                              std::vector<double>{c_rj - c_jr, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);

            // New variables
            fillUnderOverFlow(hrecoVsgen_c_han_exjet, gen_c_han_exjet_binning, reco_c_han_exjet_binning,
                              std::vector<double>{+gen_c_kk - gen_c_rr - gen_c_nn, gen_n_extraJets_iso08},
                              std::vector<double>{c_kk - c_rr - c_nn, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_sca_exjet, gen_c_sca_exjet_binning, reco_c_sca_exjet_binning,
                              std::vector<double>{-gen_c_kk + gen_c_rr - gen_c_nn, gen_n_extraJets_iso08},
                              std::vector<double>{-c_kk + c_rr - c_nn, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_tra_exjet, gen_c_tra_exjet_binning, reco_c_tra_exjet_binning,
                              std::vector<double>{-gen_c_kk - gen_c_rr + gen_c_nn, gen_n_extraJets_iso08},
                              std::vector<double>{-c_kk - c_rr + c_nn, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_kjL_exjet, gen_c_kjL_exjet_binning, reco_c_kjL_exjet_binning,
                              std::vector<double>{-gen_c_kj - gen_c_rr - gen_c_nn, gen_n_extraJets_iso08},
                              std::vector<double>{-c_kj - c_rr - c_nn, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rqL_exjet, gen_c_rqL_exjet_binning, reco_c_rqL_exjet_binning,
                              std::vector<double>{-gen_c_kk - gen_c_rq - gen_c_nn, gen_n_extraJets_iso08},
                              std::vector<double>{-c_kk - c_rq - c_nn, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rkP_exjet, gen_c_rkP_exjet_binning, reco_c_rkP_exjet_binning,
                              std::vector<double>{-gen_c_rk - gen_c_kr - gen_c_nn, gen_n_extraJets_iso08},
                              std::vector<double>{-c_rk - c_kr - c_nn, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_rkM_exjet, gen_c_rkM_exjet_binning, reco_c_rkM_exjet_binning,
                              std::vector<double>{-gen_c_rk + gen_c_kr - gen_c_nn, gen_n_extraJets_iso08},
                              std::vector<double>{-c_rk + c_kr - c_nn, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nrP_exjet, gen_c_nrP_exjet_binning, reco_c_nrP_exjet_binning,
                              std::vector<double>{-gen_c_nr - gen_c_rn - gen_c_kk, gen_n_extraJets_iso08},
                              std::vector<double>{-c_nr - c_rn - c_kk, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nrM_exjet, gen_c_nrM_exjet_binning, reco_c_nrM_exjet_binning,
                              std::vector<double>{-gen_c_nr + gen_c_rn - gen_c_kk, gen_n_extraJets_iso08},
                              std::vector<double>{-c_nr + c_rn - c_kk, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nkP_exjet, gen_c_nkP_exjet_binning, reco_c_nkP_exjet_binning,
                              std::vector<double>{-gen_c_nk - gen_c_kn - gen_c_rr, gen_n_extraJets_iso08},
                              std::vector<double>{-c_nk - c_kn - c_rr, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_c_nkM_exjet, gen_c_nkM_exjet_binning, reco_c_nkM_exjet_binning,
                              std::vector<double>{-gen_c_nk + gen_c_kn - gen_c_rr, gen_n_extraJets_iso08},
                              std::vector<double>{-c_nk + c_kn - c_rr, n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);

            // Lab fram variables
            fillUnderOverFlow(hrecoVsgen_ll_cHel_exjet, gen_ll_cHel_exjet_binning, reco_ll_cHel_exjet_binning,
                              std::vector<double>{gen_ll_cHel, gen_n_extraJets_iso08},
                              std::vector<double>{ll_cHel, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ll_cLab_exjet, gen_ll_cLab_exjet_binning, reco_ll_cLab_exjet_binning,
                              std::vector<double>{gen_ll_cLab, gen_n_extraJets_iso08},
                              std::vector<double>{ll_cLab, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ll_kNorm_exjet, gen_ll_kNorm_exjet_binning, reco_ll_kNorm_exjet_binning,
                              std::vector<double>{gen_ll_kNorm, gen_n_extraJets_iso08},
                              std::vector<double>{ll_kNorm, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hrecoVsgen_ll_rNorm_exjet, gen_ll_rNorm_exjet_binning, reco_ll_rNorm_exjet_binning,
                              std::vector<double>{gen_ll_rNorm, gen_n_extraJets_iso08},
                              std::vector<double>{ll_rNorm, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(
                hrecoVsgen_llbar_delta_phi_exjet, gen_llbar_delta_phi_exjet_binning, reco_llbar_delta_phi_exjet_binning,
                std::vector<double>{gen_llbar_delta_phi, gen_n_extraJets_iso08},
                std::vector<double>{llbar_delta_phi, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(
                hrecoVsgen_llbar_delta_eta_exjet, gen_llbar_delta_eta_exjet_binning, reco_llbar_delta_eta_exjet_binning,
                std::vector<double>{gen_llbar_delta_eta, gen_n_extraJets_iso08},
                std::vector<double>{llbar_delta_eta, n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);

            // ******
            // top_pt
            // ******

            // // Polarizations
            // fillUnderOverFlow(hrecoVsgen_b1k_top_pt, gen_b1k_top_pt_binning,
            // reco_b1k_top_pt_binning, std::vector<double>{gen_b1k,
            // gen_top_pt}, std::vector<double>{b1k, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_b2k_top_pt,
            // gen_b2k_top_pt_binning, reco_b2k_top_pt_binning,
            // std::vector<double>{gen_b2k, gen_top_pt},
            // std::vector<double>{b2k, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_b1j_top_pt,
            // gen_b1j_top_pt_binning, reco_b1j_top_pt_binning,
            // std::vector<double>{gen_b1j, gen_top_pt},
            // std::vector<double>{b1j, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_b2j_top_pt,
            // gen_b2j_top_pt_binning, reco_b2j_top_pt_binning,
            // std::vector<double>{gen_b2j, gen_top_pt},
            // std::vector<double>{b2j, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_b1r_top_pt,
            // gen_b1r_top_pt_binning, reco_b1r_top_pt_binning,
            // std::vector<double>{gen_b1r, gen_top_pt},
            // std::vector<double>{b1r, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_b2r_top_pt,
            // gen_b2r_top_pt_binning, reco_b2r_top_pt_binning,
            // std::vector<double>{gen_b2r, gen_top_pt},
            // std::vector<double>{b2r, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_b1q_top_pt,
            // gen_b1q_top_pt_binning, reco_b1q_top_pt_binning,
            // std::vector<double>{gen_b1q, gen_top_pt},
            // std::vector<double>{b1q, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_b2q_top_pt,
            // gen_b2q_top_pt_binning, reco_b2q_top_pt_binning,
            // std::vector<double>{gen_b2q, gen_top_pt},
            // std::vector<double>{b2q, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_b1n_top_pt,
            // gen_b1n_top_pt_binning, reco_b1n_top_pt_binning,
            // std::vector<double>{gen_b1n, gen_top_pt},
            // std::vector<double>{b1n, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_b2n_top_pt,
            // gen_b2n_top_pt_binning, reco_b2n_top_pt_binning,
            // std::vector<double>{gen_b2n, gen_top_pt},
            // std::vector<double>{b2n, top_pt}, eventWeight, true,
            // symmetrizeSpinVar);

            // // Diagonal elements
            // fillUnderOverFlow(hrecoVsgen_c_kk_top_pt,
            // gen_c_kk_top_pt_binning, reco_c_kk_top_pt_binning,
            // std::vector<double>{gen_c_kk, gen_top_pt},
            // std::vector<double>{c_kk, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_rr_top_pt,
            // gen_c_rr_top_pt_binning, reco_c_rr_top_pt_binning,
            // std::vector<double>{gen_c_rr, gen_top_pt},
            // std::vector<double>{c_rr, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_nn_top_pt,
            // gen_c_nn_top_pt_binning, reco_c_nn_top_pt_binning,
            // std::vector<double>{gen_c_nn, gen_top_pt},
            // std::vector<double>{c_nn, top_pt}, eventWeight, true,
            // symmetrizeSpinVar);

            // // Off diagonal elements
            // fillUnderOverFlow(hrecoVsgen_c_rk_top_pt,
            // gen_c_rk_top_pt_binning, reco_c_rk_top_pt_binning,
            // std::vector<double>{gen_c_rk, gen_top_pt},
            // std::vector<double>{c_rk, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_kr_top_pt,
            // gen_c_kr_top_pt_binning, reco_c_kr_top_pt_binning,
            // std::vector<double>{gen_c_kr, gen_top_pt},
            // std::vector<double>{c_kr, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_nr_top_pt,
            // gen_c_nr_top_pt_binning, reco_c_nr_top_pt_binning,
            // std::vector<double>{gen_c_nr, gen_top_pt},
            // std::vector<double>{c_nr, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_rn_top_pt,
            // gen_c_rn_top_pt_binning, reco_c_rn_top_pt_binning,
            // std::vector<double>{gen_c_rn, gen_top_pt},
            // std::vector<double>{c_rn, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_nk_top_pt,
            // gen_c_nk_top_pt_binning, reco_c_nk_top_pt_binning,
            // std::vector<double>{gen_c_nk, gen_top_pt},
            // std::vector<double>{c_nk, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_kn_top_pt,
            // gen_c_kn_top_pt_binning, reco_c_kn_top_pt_binning,
            // std::vector<double>{gen_c_kn, gen_top_pt},
            // std::vector<double>{c_kn, top_pt}, eventWeight, true,
            // symmetrizeSpinVar);

            // fillUnderOverFlow(hrecoVsgen_c_Prk_top_pt,
            // gen_c_Prk_top_pt_binning, reco_c_Prk_top_pt_binning,
            // std::vector<double>{gen_c_rk + gen_c_kr, gen_top_pt},
            // std::vector<double>{c_rk + c_kr, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_Mrk_top_pt,
            // gen_c_Mrk_top_pt_binning, reco_c_Mrk_top_pt_binning,
            // std::vector<double>{gen_c_rk - gen_c_kr, gen_top_pt},
            // std::vector<double>{c_rk - c_kr, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_Pnr_top_pt,
            // gen_c_Pnr_top_pt_binning, reco_c_Pnr_top_pt_binning,
            // std::vector<double>{gen_c_nr + gen_c_rn, gen_top_pt},
            // std::vector<double>{c_nr + c_rn, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_Mnr_top_pt,
            // gen_c_Mnr_top_pt_binning, reco_c_Mnr_top_pt_binning,
            // std::vector<double>{gen_c_nr - gen_c_rn, gen_top_pt},
            // std::vector<double>{c_nr - c_rn, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_Pnk_top_pt,
            // gen_c_Pnk_top_pt_binning, reco_c_Pnk_top_pt_binning,
            // std::vector<double>{gen_c_nk + gen_c_kn, gen_top_pt},
            // std::vector<double>{c_nk + c_kn, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_Mnk_top_pt,
            // gen_c_Mnk_top_pt_binning, reco_c_Mnk_top_pt_binning,
            // std::vector<double>{gen_c_nk - gen_c_kn, gen_top_pt},
            // std::vector<double>{c_nk - c_kn, top_pt}, eventWeight, true,
            // symmetrizeSpinVar);

            // // Starred axes
            // fillUnderOverFlow(hrecoVsgen_c_kj_top_pt,
            // gen_c_kj_top_pt_binning, reco_c_kj_top_pt_binning,
            // std::vector<double>{gen_c_kj, gen_top_pt},
            // std::vector<double>{c_kj, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_rq_top_pt,
            // gen_c_rq_top_pt_binning, reco_c_rq_top_pt_binning,
            // std::vector<double>{gen_c_rq, gen_top_pt},
            // std::vector<double>{c_rq, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_rj_top_pt,
            // gen_c_rj_top_pt_binning, reco_c_rj_top_pt_binning,
            // std::vector<double>{gen_c_rj, gen_top_pt},
            // std::vector<double>{c_rj, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_jr_top_pt,
            // gen_c_jr_top_pt_binning, reco_c_jr_top_pt_binning,
            // std::vector<double>{gen_c_jr, gen_top_pt},
            // std::vector<double>{c_jr, top_pt}, eventWeight, true,
            // symmetrizeSpinVar);

            // fillUnderOverFlow(hrecoVsgen_c_Prj_top_pt,
            // gen_c_Prj_top_pt_binning, reco_c_Prj_top_pt_binning,
            // std::vector<double>{gen_c_rj + gen_c_jr, gen_top_pt},
            // std::vector<double>{c_rj + c_jr, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_c_Mrj_top_pt,
            // gen_c_Mrj_top_pt_binning, reco_c_Mrj_top_pt_binning,
            // std::vector<double>{gen_c_rj - gen_c_jr, gen_top_pt},
            // std::vector<double>{c_rj - c_jr, top_pt}, eventWeight, true,
            // symmetrizeSpinVar);

            // // New variables
            // fillUnderOverFlow(hrecoVsgen_c_han_top_pt,
            // gen_c_han_top_pt_binning, reco_c_han_top_pt_binning,
            // std::vector<double>{+gen_c_kk - gen_c_rr - gen_c_nn, gen_top_pt},
            // std::vector<double>{ c_kk - c_rr - c_nn, top_pt}, eventWeight,
            // true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_sca_top_pt,
            // gen_c_sca_top_pt_binning, reco_c_sca_top_pt_binning,
            // std::vector<double>{-gen_c_kk + gen_c_rr - gen_c_nn, gen_top_pt},
            // std::vector<double>{-c_kk + c_rr - c_nn, top_pt}, eventWeight,
            // true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_tra_top_pt,
            // gen_c_tra_top_pt_binning, reco_c_tra_top_pt_binning,
            // std::vector<double>{-gen_c_kk - gen_c_rr + gen_c_nn, gen_top_pt},
            // std::vector<double>{-c_kk - c_rr + c_nn, top_pt}, eventWeight,
            // true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_kjL_top_pt,
            // gen_c_kjL_top_pt_binning, reco_c_kjL_top_pt_binning,
            // std::vector<double>{-gen_c_kj - gen_c_rr - gen_c_nn, gen_top_pt},
            // std::vector<double>{-c_kj - c_rr - c_nn, top_pt}, eventWeight,
            // true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rqL_top_pt,
            // gen_c_rqL_top_pt_binning, reco_c_rqL_top_pt_binning,
            // std::vector<double>{-gen_c_kk - gen_c_rq - gen_c_nn, gen_top_pt},
            // std::vector<double>{-c_kk - c_rq - c_nn, top_pt}, eventWeight,
            // true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rkP_top_pt,
            // gen_c_rkP_top_pt_binning, reco_c_rkP_top_pt_binning,
            // std::vector<double>{-gen_c_rk - gen_c_kr - gen_c_nn, gen_top_pt},
            // std::vector<double>{-c_rk - c_kr - c_nn, top_pt}, eventWeight,
            // true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rkM_top_pt,
            // gen_c_rkM_top_pt_binning, reco_c_rkM_top_pt_binning,
            // std::vector<double>{-gen_c_rk + gen_c_kr - gen_c_nn, gen_top_pt},
            // std::vector<double>{-c_rk + c_kr - c_nn, top_pt}, eventWeight,
            // true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nrP_top_pt,
            // gen_c_nrP_top_pt_binning, reco_c_nrP_top_pt_binning,
            // std::vector<double>{-gen_c_nr - gen_c_rn - gen_c_kk, gen_top_pt},
            // std::vector<double>{-c_nr - c_rn - c_kk, top_pt}, eventWeight,
            // true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nrM_top_pt,
            // gen_c_nrM_top_pt_binning, reco_c_nrM_top_pt_binning,
            // std::vector<double>{-gen_c_nr + gen_c_rn - gen_c_kk, gen_top_pt},
            // std::vector<double>{-c_nr + c_rn - c_kk, top_pt}, eventWeight,
            // true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nkP_top_pt,
            // gen_c_nkP_top_pt_binning, reco_c_nkP_top_pt_binning,
            // std::vector<double>{-gen_c_nk - gen_c_kn - gen_c_rr, gen_top_pt},
            // std::vector<double>{-c_nk - c_kn - c_rr, top_pt}, eventWeight,
            // true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nkM_top_pt,
            // gen_c_nkM_top_pt_binning, reco_c_nkM_top_pt_binning,
            // std::vector<double>{-gen_c_nk + gen_c_kn - gen_c_rr, gen_top_pt},
            // std::vector<double>{-c_nk + c_kn - c_rr, top_pt}, eventWeight,
            // true, symmetrizeSpinVar);

            // // Lab fram variables
            // fillUnderOverFlow(hrecoVsgen_ll_cHel_top_pt,
            // gen_ll_cHel_top_pt_binning, reco_ll_cHel_top_pt_binning,
            // std::vector<double>{gen_ll_cHel, gen_top_pt},
            // std::vector<double>{ll_cHel, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_ll_cLab_top_pt,
            // gen_ll_cLab_top_pt_binning, reco_ll_cLab_top_pt_binning,
            // std::vector<double>{gen_ll_cLab, gen_top_pt},
            // std::vector<double>{ll_cLab, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_ll_kNorm_top_pt,
            // gen_ll_kNorm_top_pt_binning, reco_ll_kNorm_top_pt_binning,
            // std::vector<double>{gen_ll_kNorm, gen_top_pt},
            // std::vector<double>{ll_kNorm, top_pt}, eventWeight, true,
            // symmetrizeSpinVar); fillUnderOverFlow(hrecoVsgen_ll_rNorm_top_pt,
            // gen_ll_rNorm_top_pt_binning, reco_ll_rNorm_top_pt_binning,
            // std::vector<double>{gen_ll_rNorm, gen_top_pt},
            // std::vector<double>{ll_rNorm, top_pt}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_llbar_delta_phi_top_pt,
            // gen_llbar_delta_phi_top_pt_binning,
            // reco_llbar_delta_phi_top_pt_binning,
            // std::vector<double>{gen_llbar_delta_phi, gen_top_pt},
            // std::vector<double>{llbar_delta_phi, top_pt}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_llbar_delta_eta_top_pt,
            // gen_llbar_delta_eta_top_pt_binning,
            // reco_llbar_delta_eta_top_pt_binning,
            // std::vector<double>{gen_llbar_delta_eta, gen_top_pt},
            // std::vector<double>{llbar_delta_eta, top_pt}, eventWeight, true,
            // symmetrizeSpinVar);

            // ******************************
            // top_scatteringangle_ttbarframe
            // ******************************

            // // Polarizations
            // fillUnderOverFlow(hrecoVsgen_b1k_top_scatteringangle_ttbarframe,
            // gen_b1k_top_scatteringangle_ttbarframe_binning,
            // reco_b1k_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1k, gen_top_scatteringangle_ttbarframe},
            // std::vector<double>{b1k, top_scatteringangle_ttbarframe},
            // eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b2k_top_scatteringangle_ttbarframe,
            // gen_b2k_top_scatteringangle_ttbarframe_binning,
            // reco_b2k_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2k, gen_top_scatteringangle_ttbarframe},
            // std::vector<double>{b2k, top_scatteringangle_ttbarframe},
            // eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b1j_top_scatteringangle_ttbarframe,
            // gen_b1j_top_scatteringangle_ttbarframe_binning,
            // reco_b1j_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1j, gen_top_scatteringangle_ttbarframe},
            // std::vector<double>{b1j, top_scatteringangle_ttbarframe},
            // eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b2j_top_scatteringangle_ttbarframe,
            // gen_b2j_top_scatteringangle_ttbarframe_binning,
            // reco_b2j_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2j, gen_top_scatteringangle_ttbarframe},
            // std::vector<double>{b2j, top_scatteringangle_ttbarframe},
            // eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b1r_top_scatteringangle_ttbarframe,
            // gen_b1r_top_scatteringangle_ttbarframe_binning,
            // reco_b1r_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1r, gen_top_scatteringangle_ttbarframe},
            // std::vector<double>{b1r, top_scatteringangle_ttbarframe},
            // eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b2r_top_scatteringangle_ttbarframe,
            // gen_b2r_top_scatteringangle_ttbarframe_binning,
            // reco_b2r_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2r, gen_top_scatteringangle_ttbarframe},
            // std::vector<double>{b2r, top_scatteringangle_ttbarframe},
            // eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b1q_top_scatteringangle_ttbarframe,
            // gen_b1q_top_scatteringangle_ttbarframe_binning,
            // reco_b1q_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1q, gen_top_scatteringangle_ttbarframe},
            // std::vector<double>{b1q, top_scatteringangle_ttbarframe},
            // eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b2q_top_scatteringangle_ttbarframe,
            // gen_b2q_top_scatteringangle_ttbarframe_binning,
            // reco_b2q_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2q, gen_top_scatteringangle_ttbarframe},
            // std::vector<double>{b2q, top_scatteringangle_ttbarframe},
            // eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b1n_top_scatteringangle_ttbarframe,
            // gen_b1n_top_scatteringangle_ttbarframe_binning,
            // reco_b1n_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1n, gen_top_scatteringangle_ttbarframe},
            // std::vector<double>{b1n, top_scatteringangle_ttbarframe},
            // eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_b2n_top_scatteringangle_ttbarframe,
            // gen_b2n_top_scatteringangle_ttbarframe_binning,
            // reco_b2n_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2n, gen_top_scatteringangle_ttbarframe},
            // std::vector<double>{b2n, top_scatteringangle_ttbarframe},
            // eventWeight, true, symmetrizeSpinVar);

            // // Diagonal elements
            // fillUnderOverFlow(hrecoVsgen_c_kk_top_scatteringangle_ttbarframe,
            // gen_c_kk_top_scatteringangle_ttbarframe_binning,
            // reco_c_kk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_kk,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_kk,
            // top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rr_top_scatteringangle_ttbarframe,
            // gen_c_rr_top_scatteringangle_ttbarframe_binning,
            // reco_c_rr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rr,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_rr,
            // top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nn_top_scatteringangle_ttbarframe,
            // gen_c_nn_top_scatteringangle_ttbarframe_binning,
            // reco_c_nn_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nn,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_nn,
            // top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);

            // // Off diagonal elements
            // fillUnderOverFlow(hrecoVsgen_c_rk_top_scatteringangle_ttbarframe,
            // gen_c_rk_top_scatteringangle_ttbarframe_binning,
            // reco_c_rk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rk,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_rk,
            // top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_kr_top_scatteringangle_ttbarframe,
            // gen_c_kr_top_scatteringangle_ttbarframe_binning,
            // reco_c_kr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_kr,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_kr,
            // top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nr_top_scatteringangle_ttbarframe,
            // gen_c_nr_top_scatteringangle_ttbarframe_binning,
            // reco_c_nr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nr,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_nr,
            // top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rn_top_scatteringangle_ttbarframe,
            // gen_c_rn_top_scatteringangle_ttbarframe_binning,
            // reco_c_rn_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rn,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_rn,
            // top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nk_top_scatteringangle_ttbarframe,
            // gen_c_nk_top_scatteringangle_ttbarframe_binning,
            // reco_c_nk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nk,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_nk,
            // top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_kn_top_scatteringangle_ttbarframe,
            // gen_c_kn_top_scatteringangle_ttbarframe_binning,
            // reco_c_kn_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_kn,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_kn,
            // top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);

            // fillUnderOverFlow(hrecoVsgen_c_Prk_top_scatteringangle_ttbarframe,
            // gen_c_Prk_top_scatteringangle_ttbarframe_binning,
            // reco_c_Prk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rk + gen_c_kr,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_rk +
            // c_kr, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_Mrk_top_scatteringangle_ttbarframe,
            // gen_c_Mrk_top_scatteringangle_ttbarframe_binning,
            // reco_c_Mrk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rk - gen_c_kr,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_rk -
            // c_kr, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_Pnr_top_scatteringangle_ttbarframe,
            // gen_c_Pnr_top_scatteringangle_ttbarframe_binning,
            // reco_c_Pnr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nr + gen_c_rn,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_nr +
            // c_rn, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_Mnr_top_scatteringangle_ttbarframe,
            // gen_c_Mnr_top_scatteringangle_ttbarframe_binning,
            // reco_c_Mnr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nr - gen_c_rn,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_nr -
            // c_rn, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_Pnk_top_scatteringangle_ttbarframe,
            // gen_c_Pnk_top_scatteringangle_ttbarframe_binning,
            // reco_c_Pnk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nk + gen_c_kn,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_nk +
            // c_kn, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_Mnk_top_scatteringangle_ttbarframe,
            // gen_c_Mnk_top_scatteringangle_ttbarframe_binning,
            // reco_c_Mnk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nk - gen_c_kn,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_nk -
            // c_kn, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);

            // // Starred axes
            // fillUnderOverFlow(hrecoVsgen_c_kj_top_scatteringangle_ttbarframe,
            // gen_c_kj_top_scatteringangle_ttbarframe_binning,
            // reco_c_kj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_kj,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_kj,
            // top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rq_top_scatteringangle_ttbarframe,
            // gen_c_rq_top_scatteringangle_ttbarframe_binning,
            // reco_c_rq_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rq,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_rq,
            // top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rj_top_scatteringangle_ttbarframe,
            // gen_c_rj_top_scatteringangle_ttbarframe_binning,
            // reco_c_rj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rj,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_rj,
            // top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_jr_top_scatteringangle_ttbarframe,
            // gen_c_jr_top_scatteringangle_ttbarframe_binning,
            // reco_c_jr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_jr,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_jr,
            // top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);

            // fillUnderOverFlow(hrecoVsgen_c_Prj_top_scatteringangle_ttbarframe,
            // gen_c_Prj_top_scatteringangle_ttbarframe_binning,
            // reco_c_Prj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rj + gen_c_jr,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_rj +
            // c_jr, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_Mrj_top_scatteringangle_ttbarframe,
            // gen_c_Mrj_top_scatteringangle_ttbarframe_binning,
            // reco_c_Mrj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rj - gen_c_jr,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{c_rj -
            // c_jr, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);

            // // New variables
            // fillUnderOverFlow(hrecoVsgen_c_han_top_scatteringangle_ttbarframe,
            // gen_c_han_top_scatteringangle_ttbarframe_binning,
            // reco_c_han_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{+gen_c_kk - gen_c_rr - gen_c_nn,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{ c_kk -
            // c_rr - c_nn, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_sca_top_scatteringangle_ttbarframe,
            // gen_c_sca_top_scatteringangle_ttbarframe_binning,
            // reco_c_sca_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_kk + gen_c_rr - gen_c_nn,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{-c_kk +
            // c_rr - c_nn, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_tra_top_scatteringangle_ttbarframe,
            // gen_c_tra_top_scatteringangle_ttbarframe_binning,
            // reco_c_tra_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_kk - gen_c_rr + gen_c_nn,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{-c_kk -
            // c_rr + c_nn, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_kjL_top_scatteringangle_ttbarframe,
            // gen_c_kjL_top_scatteringangle_ttbarframe_binning,
            // reco_c_kjL_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_kj - gen_c_rr - gen_c_nn,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{-c_kj -
            // c_rr - c_nn, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rqL_top_scatteringangle_ttbarframe,
            // gen_c_rqL_top_scatteringangle_ttbarframe_binning,
            // reco_c_rqL_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_kk - gen_c_rq - gen_c_nn,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{-c_kk -
            // c_rq - c_nn, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rkP_top_scatteringangle_ttbarframe,
            // gen_c_rkP_top_scatteringangle_ttbarframe_binning,
            // reco_c_rkP_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_rk - gen_c_kr - gen_c_nn,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{-c_rk -
            // c_kr - c_nn, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_rkM_top_scatteringangle_ttbarframe,
            // gen_c_rkM_top_scatteringangle_ttbarframe_binning,
            // reco_c_rkM_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_rk + gen_c_kr - gen_c_nn,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{-c_rk +
            // c_kr - c_nn, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nrP_top_scatteringangle_ttbarframe,
            // gen_c_nrP_top_scatteringangle_ttbarframe_binning,
            // reco_c_nrP_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_nr - gen_c_rn - gen_c_kk,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{-c_nr -
            // c_rn - c_kk, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nrM_top_scatteringangle_ttbarframe,
            // gen_c_nrM_top_scatteringangle_ttbarframe_binning,
            // reco_c_nrM_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_nr + gen_c_rn - gen_c_kk,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{-c_nr +
            // c_rn - c_kk, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nkP_top_scatteringangle_ttbarframe,
            // gen_c_nkP_top_scatteringangle_ttbarframe_binning,
            // reco_c_nkP_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_nk - gen_c_kn - gen_c_rr,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{-c_nk -
            // c_kn - c_rr, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_c_nkM_top_scatteringangle_ttbarframe,
            // gen_c_nkM_top_scatteringangle_ttbarframe_binning,
            // reco_c_nkM_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_nk + gen_c_kn - gen_c_rr,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{-c_nk +
            // c_kn - c_rr, top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);

            // // Lab fram variables
            // fillUnderOverFlow(hrecoVsgen_ll_cHel_top_scatteringangle_ttbarframe,
            // gen_ll_cHel_top_scatteringangle_ttbarframe_binning,
            // reco_ll_cHel_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_ll_cHel,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{ll_cHel,
            // top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_ll_cLab_top_scatteringangle_ttbarframe,
            // gen_ll_cLab_top_scatteringangle_ttbarframe_binning,
            // reco_ll_cLab_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_ll_cLab,
            // gen_top_scatteringangle_ttbarframe}, std::vector<double>{ll_cLab,
            // top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_ll_kNorm_top_scatteringangle_ttbarframe,
            // gen_ll_kNorm_top_scatteringangle_ttbarframe_binning,
            // reco_ll_kNorm_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_ll_kNorm,
            // gen_top_scatteringangle_ttbarframe},
            // std::vector<double>{ll_kNorm, top_scatteringangle_ttbarframe},
            // eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_ll_rNorm_top_scatteringangle_ttbarframe,
            // gen_ll_rNorm_top_scatteringangle_ttbarframe_binning,
            // reco_ll_rNorm_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_ll_rNorm,
            // gen_top_scatteringangle_ttbarframe},
            // std::vector<double>{ll_rNorm, top_scatteringangle_ttbarframe},
            // eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_llbar_delta_phi_top_scatteringangle_ttbarframe,
            // gen_llbar_delta_phi_top_scatteringangle_ttbarframe_binning,
            // reco_llbar_delta_phi_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_llbar_delta_phi,
            // gen_top_scatteringangle_ttbarframe},
            // std::vector<double>{llbar_delta_phi,
            // top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hrecoVsgen_llbar_delta_eta_top_scatteringangle_ttbarframe,
            // gen_llbar_delta_eta_top_scatteringangle_ttbarframe_binning,
            // reco_llbar_delta_eta_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_llbar_delta_eta,
            // gen_top_scatteringangle_ttbarframe},
            // std::vector<double>{llbar_delta_eta,
            // top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);

            // ****************************************
            // Fill Migration matrices for 1D variables
            // ****************************************

            // Cases where eventWeight < trueLevelWeight implies gen events that
            // are not reconstructed Fill an underflow bin in reco where all
            // these missed events are kept (trueLevelWeight - eventWeight)

            // Kinematic
            fillUnderOverFlowDeltaWeight(hrecoVsgen_top_pt, gen_top_pt_binning, std::vector<double>{gen_top_pt},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_l_pt, gen_l_pt_binning, std::vector<double>{gen_l_pt},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_lbar_pt, gen_lbar_pt_binning, std::vector<double>{gen_lbar_pt},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_ttbar_pt, gen_ttbar_pt_binning, std::vector<double>{gen_ttbar_pt},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_ttbar_mass, gen_ttbar_mass_binning,
                                         std::vector<double>{gen_ttbar_mass}, trueLevelWeight, eventWeight,
                                         symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_ttbar_delta_phi, gen_ttbar_delta_phi_binning,
                                         std::vector<double>{gen_ttbar_delta_phi}, trueLevelWeight, eventWeight,
                                         symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_ttbar_delta_eta, gen_ttbar_delta_eta_binning,
                                         std::vector<double>{gen_ttbar_delta_eta}, trueLevelWeight, eventWeight,
                                         symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_ttbar_rapidity, gen_ttbar_rapidity_binning,
                                         std::vector<double>{gen_ttbar_rapidity}, trueLevelWeight, eventWeight,
                                         symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_llbar_pt, gen_llbar_pt_binning, std::vector<double>{gen_llbar_pt},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_llbar_mass, gen_llbar_mass_binning,
                                         std::vector<double>{gen_llbar_mass}, trueLevelWeight, eventWeight,
                                         symmetrizeSpinVar);

            // Spin corr

            // Polarizations
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b1k, gen_b1k_binning, std::vector<double>{gen_b1k}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b2k, gen_b2k_binning, std::vector<double>{gen_b2k}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b1j, gen_b1j_binning, std::vector<double>{gen_b1j}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b2j, gen_b2j_binning, std::vector<double>{gen_b2j}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b1r, gen_b1r_binning, std::vector<double>{gen_b1r}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b2r, gen_b2r_binning, std::vector<double>{gen_b2r}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b1q, gen_b1q_binning, std::vector<double>{gen_b1q}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b2q, gen_b2q_binning, std::vector<double>{gen_b2q}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b1n, gen_b1n_binning, std::vector<double>{gen_b1n}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b2n, gen_b2n_binning, std::vector<double>{gen_b2n}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);

            // Diagonal elements
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kk, gen_c_kk_binning, std::vector<double>{gen_c_kk},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rr, gen_c_rr_binning, std::vector<double>{gen_c_rr},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nn, gen_c_nn_binning, std::vector<double>{gen_c_nn},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);

            // Off diagonal elements
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rk, gen_c_rk_binning, std::vector<double>{gen_c_rk},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kr, gen_c_kr_binning, std::vector<double>{gen_c_kr},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nr, gen_c_nr_binning, std::vector<double>{gen_c_nr},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rn, gen_c_rn_binning, std::vector<double>{gen_c_rn},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nk, gen_c_nk_binning, std::vector<double>{gen_c_nk},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kn, gen_c_kn_binning, std::vector<double>{gen_c_kn},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);

            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Prk, gen_c_Prk_binning, std::vector<double>{gen_c_rk + gen_c_kr},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Mrk, gen_c_Mrk_binning, std::vector<double>{gen_c_rk - gen_c_kr},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Pnr, gen_c_Pnr_binning, std::vector<double>{gen_c_nr + gen_c_rn},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Mnr, gen_c_Mnr_binning, std::vector<double>{gen_c_nr - gen_c_rn},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Pnk, gen_c_Pnk_binning, std::vector<double>{gen_c_nk + gen_c_kn},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Mnk, gen_c_Mnk_binning, std::vector<double>{gen_c_nk - gen_c_kn},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);

            // Starred axes
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kj, gen_c_kj_binning, std::vector<double>{gen_c_kj},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rq, gen_c_rq_binning, std::vector<double>{gen_c_rq},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kj, gen_c_rj_binning, std::vector<double>{gen_c_rj},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rq, gen_c_jr_binning, std::vector<double>{gen_c_jr},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);

            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Prj, gen_c_Prj_binning, std::vector<double>{gen_c_rj + gen_c_jr},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Mrj, gen_c_Mrj_binning, std::vector<double>{gen_c_rj - gen_c_jr},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);

            // New variables
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_han, gen_c_han_binning,
                                         std::vector<double>{gen_c_kk - gen_c_rr - gen_c_nn}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_sca, gen_c_sca_binning,
                                         std::vector<double>{-gen_c_kk + gen_c_rr - gen_c_nn}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_tra, gen_c_tra_binning,
                                         std::vector<double>{-gen_c_kk - gen_c_rr + gen_c_nn}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kjL, gen_c_kjL_binning,
                                         std::vector<double>{-gen_c_kj - gen_c_rr - gen_c_nn}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rqL, gen_c_rqL_binning,
                                         std::vector<double>{-gen_c_kk - gen_c_rq - gen_c_nn}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rkP, gen_c_rkP_binning,
                                         std::vector<double>{-gen_c_rk - gen_c_kr - gen_c_nn}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rkM, gen_c_rkM_binning,
                                         std::vector<double>{-gen_c_rk + gen_c_kr - gen_c_nn}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nrP, gen_c_nrP_binning,
                                         std::vector<double>{-gen_c_nr - gen_c_rn - gen_c_kk}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nrM, gen_c_nrM_binning,
                                         std::vector<double>{-gen_c_nr + gen_c_rn - gen_c_kk}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nkP, gen_c_nkP_binning,
                                         std::vector<double>{-gen_c_nk - gen_c_kn - gen_c_rr}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nkM, gen_c_nkM_binning,
                                         std::vector<double>{-gen_c_nk + gen_c_kn - gen_c_rr}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);

            // Lab frame variables
            fillUnderOverFlowDeltaWeight(hrecoVsgen_ll_cHel, gen_ll_cHel_binning, std::vector<double>{gen_ll_cHel},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_ll_cLab, gen_ll_cLab_binning, std::vector<double>{gen_ll_cLab},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_ll_kNorm, gen_ll_kNorm_binning, std::vector<double>{gen_ll_kNorm},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_ll_rNorm, gen_ll_rNorm_binning, std::vector<double>{gen_ll_rNorm},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_llbar_delta_phi, gen_llbar_delta_phi_binning,
                                         std::vector<double>{gen_llbar_delta_phi}, trueLevelWeight, eventWeight,
                                         symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_llbar_delta_eta, gen_llbar_delta_eta_binning,
                                         std::vector<double>{gen_llbar_delta_eta}, trueLevelWeight, eventWeight,
                                         symmetrizeSpinVar);

            // ****************************************
            // Fill Migration matrices for 2D variables
            // ****************************************

            // Amandeep : Fill 2D migration matrices
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kk_mttbar,
            // gen_c_kk_mttbar_binning, std::vector<double>{gen_c_kk,
            // gen_ttbar_mass}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);

            // ******
            // mttbar
            // ******

            // Polarizations
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b1k_exjet, gen_b1k_exjet_binning,
                                         std::vector<double>{gen_b1k, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b2k_exjet, gen_b2k_exjet_binning,
                                         std::vector<double>{gen_b2k, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b1j_exjet, gen_b1j_exjet_binning,
                                         std::vector<double>{gen_b1j, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b2j_exjet, gen_b2j_exjet_binning,
                                         std::vector<double>{gen_b2j, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b1r_exjet, gen_b1r_exjet_binning,
                                         std::vector<double>{gen_b1r, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b2r_exjet, gen_b2r_exjet_binning,
                                         std::vector<double>{gen_b2r, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b1q_exjet, gen_b1q_exjet_binning,
                                         std::vector<double>{gen_b1q, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b2q_exjet, gen_b2q_exjet_binning,
                                         std::vector<double>{gen_b2q, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b1n_exjet, gen_b1n_exjet_binning,
                                         std::vector<double>{gen_b1n, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_b2n_exjet, gen_b2n_exjet_binning,
                                         std::vector<double>{gen_b2n, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);

            // Diagonal elements
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kk_exjet, gen_c_kk_exjet_binning,
                                         std::vector<double>{gen_c_kk, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rr_exjet, gen_c_rr_exjet_binning,
                                         std::vector<double>{gen_c_rr, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nn_exjet, gen_c_nn_exjet_binning,
                                         std::vector<double>{gen_c_nn, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);

            // Off diagonal elements
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rk_exjet, gen_c_rk_exjet_binning,
                                         std::vector<double>{gen_c_rk, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kr_exjet, gen_c_kr_exjet_binning,
                                         std::vector<double>{gen_c_kr, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nr_exjet, gen_c_nr_exjet_binning,
                                         std::vector<double>{gen_c_nr, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rn_exjet, gen_c_rn_exjet_binning,
                                         std::vector<double>{gen_c_rn, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nk_exjet, gen_c_nk_exjet_binning,
                                         std::vector<double>{gen_c_nk, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kn_exjet, gen_c_kn_exjet_binning,
                                         std::vector<double>{gen_c_kn, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);

            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Prk_exjet, gen_c_Prk_exjet_binning,
                                         std::vector<double>{gen_c_rk + gen_c_kr, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Mrk_exjet, gen_c_Mrk_exjet_binning,
                                         std::vector<double>{gen_c_rk - gen_c_kr, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Pnr_exjet, gen_c_Pnr_exjet_binning,
                                         std::vector<double>{gen_c_nr + gen_c_rn, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Mnr_exjet, gen_c_Mnr_exjet_binning,
                                         std::vector<double>{gen_c_nr - gen_c_rn, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Pnk_exjet, gen_c_Pnk_exjet_binning,
                                         std::vector<double>{gen_c_nk + gen_c_kn, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Mnk_exjet, gen_c_Mnk_exjet_binning,
                                         std::vector<double>{gen_c_nk - gen_c_kn, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);

            // Starred axes
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kj_exjet, gen_c_kj_exjet_binning,
                                         std::vector<double>{gen_c_kj, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rq_exjet, gen_c_rq_exjet_binning,
                                         std::vector<double>{gen_c_rq, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rj_exjet, gen_c_rj_exjet_binning,
                                         std::vector<double>{gen_c_rj, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_jr_exjet, gen_c_jr_exjet_binning,
                                         std::vector<double>{gen_c_jr, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);

            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Prj_exjet, gen_c_Prj_exjet_binning,
                                         std::vector<double>{gen_c_rj + gen_c_jr, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Mrj_exjet, gen_c_Mrj_exjet_binning,
                                         std::vector<double>{gen_c_rj - gen_c_jr, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);

            // New variables
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_han_exjet, gen_c_han_exjet_binning,
                                         std::vector<double>{gen_c_kk - gen_c_rr - gen_c_nn, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_sca_exjet, gen_c_sca_exjet_binning,
                                         std::vector<double>{-gen_c_kk + gen_c_rr - gen_c_nn, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_tra_exjet, gen_c_tra_exjet_binning,
                                         std::vector<double>{-gen_c_kk - gen_c_rr + gen_c_nn, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kjL_exjet, gen_c_kjL_exjet_binning,
                                         std::vector<double>{-gen_c_kj - gen_c_rr - gen_c_nn, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rqL_exjet, gen_c_rqL_exjet_binning,
                                         std::vector<double>{-gen_c_kk - gen_c_rq - gen_c_nn, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rkP_exjet, gen_c_rkP_exjet_binning,
                                         std::vector<double>{-gen_c_rk - gen_c_kr - gen_c_nn, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rkM_exjet, gen_c_rkM_exjet_binning,
                                         std::vector<double>{-gen_c_rk + gen_c_kr - gen_c_nn, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nrP_exjet, gen_c_nrP_exjet_binning,
                                         std::vector<double>{-gen_c_nr - gen_c_rn - gen_c_kk, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nrM_exjet, gen_c_nrM_exjet_binning,
                                         std::vector<double>{-gen_c_nr + gen_c_rn - gen_c_kk, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nkP_exjet, gen_c_nkP_exjet_binning,
                                         std::vector<double>{-gen_c_nk - gen_c_kn - gen_c_rr, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nkM_exjet, gen_c_nkM_exjet_binning,
                                         std::vector<double>{-gen_c_nk + gen_c_kn - gen_c_rr, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);

            // Lab frame variables
            fillUnderOverFlowDeltaWeight(hrecoVsgen_ll_cHel_exjet, gen_ll_cHel_exjet_binning,
                                         std::vector<double>{gen_ll_cHel, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_ll_cLab_exjet, gen_ll_cLab_exjet_binning,
                                         std::vector<double>{gen_ll_cLab, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_ll_kNorm_exjet, gen_ll_kNorm_exjet_binning,
                                         std::vector<double>{gen_ll_kNorm, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_ll_rNorm_exjet, gen_ll_rNorm_exjet_binning,
                                         std::vector<double>{gen_ll_rNorm, gen_n_extraJets_iso08}, trueLevelWeight,
                                         eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_llbar_delta_phi_exjet, gen_llbar_delta_phi_exjet_binning,
                                         std::vector<double>{gen_llbar_delta_phi, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);
            fillUnderOverFlowDeltaWeight(hrecoVsgen_llbar_delta_eta_exjet, gen_llbar_delta_eta_exjet_binning,
                                         std::vector<double>{gen_llbar_delta_eta, gen_n_extraJets_iso08},
                                         trueLevelWeight, eventWeight, symmetrizeSpinVar);

            // ******
            // top_pt
            // ******

            // // Polarizations
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b1k_top_pt,
            // gen_b1k_top_pt_binning, std::vector<double>{gen_b1k, gen_top_pt},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b2k_top_pt,
            // gen_b2k_top_pt_binning, std::vector<double>{gen_b2k, gen_top_pt},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b1j_top_pt,
            // gen_b1j_top_pt_binning, std::vector<double>{gen_b1j, gen_top_pt},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b2j_top_pt,
            // gen_b2j_top_pt_binning, std::vector<double>{gen_b2j, gen_top_pt},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b1r_top_pt,
            // gen_b1r_top_pt_binning, std::vector<double>{gen_b1r, gen_top_pt},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b2r_top_pt,
            // gen_b2r_top_pt_binning, std::vector<double>{gen_b2r, gen_top_pt},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b1q_top_pt,
            // gen_b1q_top_pt_binning, std::vector<double>{gen_b1q, gen_top_pt},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b2q_top_pt,
            // gen_b2q_top_pt_binning, std::vector<double>{gen_b2q, gen_top_pt},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b1n_top_pt,
            // gen_b1n_top_pt_binning, std::vector<double>{gen_b1n, gen_top_pt},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b2n_top_pt,
            // gen_b2n_top_pt_binning, std::vector<double>{gen_b2n, gen_top_pt},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);

            // // Diagonal elements
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kk_top_pt,
            // gen_c_kk_top_pt_binning, std::vector<double>{gen_c_kk,
            // gen_top_pt}, trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rr_top_pt,
            // gen_c_rr_top_pt_binning, std::vector<double>{gen_c_rr,
            // gen_top_pt}, trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nn_top_pt,
            // gen_c_nn_top_pt_binning, std::vector<double>{gen_c_nn,
            // gen_top_pt}, trueLevelWeight, eventWeight, symmetrizeSpinVar);

            // // Off diagonal elements
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rk_top_pt,
            // gen_c_rk_top_pt_binning, std::vector<double>{gen_c_rk,
            // gen_top_pt}, trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kr_top_pt,
            // gen_c_kr_top_pt_binning, std::vector<double>{gen_c_kr,
            // gen_top_pt}, trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nr_top_pt,
            // gen_c_nr_top_pt_binning, std::vector<double>{gen_c_nr,
            // gen_top_pt}, trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rn_top_pt,
            // gen_c_rn_top_pt_binning, std::vector<double>{gen_c_rn,
            // gen_top_pt}, trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nk_top_pt,
            // gen_c_nk_top_pt_binning, std::vector<double>{gen_c_nk,
            // gen_top_pt}, trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kn_top_pt,
            // gen_c_kn_top_pt_binning, std::vector<double>{gen_c_kn,
            // gen_top_pt}, trueLevelWeight, eventWeight, symmetrizeSpinVar);

            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Prk_top_pt,
            // gen_c_Prk_top_pt_binning, std::vector<double>{gen_c_rk +
            // gen_c_kr, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Mrk_top_pt,
            // gen_c_Mrk_top_pt_binning, std::vector<double>{gen_c_rk -
            // gen_c_kr, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Pnr_top_pt,
            // gen_c_Pnr_top_pt_binning, std::vector<double>{gen_c_nr +
            // gen_c_rn, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Mnr_top_pt,
            // gen_c_Mnr_top_pt_binning, std::vector<double>{gen_c_nr -
            // gen_c_rn, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Pnk_top_pt,
            // gen_c_Pnk_top_pt_binning, std::vector<double>{gen_c_nk +
            // gen_c_kn, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Mnk_top_pt,
            // gen_c_Mnk_top_pt_binning, std::vector<double>{gen_c_nk -
            // gen_c_kn, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);

            // // Starred axes
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kj_top_pt,
            // gen_c_kj_top_pt_binning, std::vector<double>{gen_c_kj,
            // gen_top_pt}, trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rq_top_pt,
            // gen_c_rq_top_pt_binning, std::vector<double>{gen_c_rq,
            // gen_top_pt}, trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rj_top_pt,
            // gen_c_rj_top_pt_binning, std::vector<double>{gen_c_rj,
            // gen_top_pt}, trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_jr_top_pt,
            // gen_c_jr_top_pt_binning, std::vector<double>{gen_c_jr,
            // gen_top_pt}, trueLevelWeight, eventWeight, symmetrizeSpinVar);

            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Prj_top_pt,
            // gen_c_Prj_top_pt_binning, std::vector<double>{gen_c_rj +
            // gen_c_jr, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Mrj_top_pt,
            // gen_c_Mrj_top_pt_binning, std::vector<double>{gen_c_rj -
            // gen_c_jr, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);

            // // New variables
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_han_top_pt,
            // gen_c_han_top_pt_binning, std::vector<double>{ gen_c_kk -
            // gen_c_rr - gen_c_nn, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_sca_top_pt,
            // gen_c_sca_top_pt_binning, std::vector<double>{-gen_c_kk +
            // gen_c_rr - gen_c_nn, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_tra_top_pt,
            // gen_c_tra_top_pt_binning, std::vector<double>{-gen_c_kk -
            // gen_c_rr + gen_c_nn, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kjL_top_pt,
            // gen_c_kjL_top_pt_binning, std::vector<double>{-gen_c_kj -
            // gen_c_rr - gen_c_nn, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rqL_top_pt,
            // gen_c_rqL_top_pt_binning, std::vector<double>{-gen_c_kk -
            // gen_c_rq - gen_c_nn, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rkP_top_pt,
            // gen_c_rkP_top_pt_binning, std::vector<double>{-gen_c_rk -
            // gen_c_kr - gen_c_nn, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rkM_top_pt,
            // gen_c_rkM_top_pt_binning, std::vector<double>{-gen_c_rk +
            // gen_c_kr - gen_c_nn, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nrP_top_pt,
            // gen_c_nrP_top_pt_binning, std::vector<double>{-gen_c_nr -
            // gen_c_rn - gen_c_kk, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nrM_top_pt,
            // gen_c_nrM_top_pt_binning, std::vector<double>{-gen_c_nr +
            // gen_c_rn - gen_c_kk, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nkP_top_pt,
            // gen_c_nkP_top_pt_binning, std::vector<double>{-gen_c_nk -
            // gen_c_kn - gen_c_rr, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nkM_top_pt,
            // gen_c_nkM_top_pt_binning, std::vector<double>{-gen_c_nk +
            // gen_c_kn - gen_c_rr, gen_top_pt}, trueLevelWeight, eventWeight,
            // symmetrizeSpinVar);

            // // Lab frame variables
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_ll_cHel_top_pt,
            // gen_ll_cHel_top_pt_binning, std::vector<double>{gen_ll_cHel,
            // gen_top_pt}, trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_ll_cLab_top_pt,
            // gen_ll_cLab_top_pt_binning, std::vector<double>{gen_ll_cLab,
            // gen_top_pt}, trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_ll_kNorm_top_pt,
            // gen_ll_kNorm_top_pt_binning, std::vector<double>{gen_ll_kNorm,
            // gen_top_pt}, trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_ll_rNorm_top_pt,
            // gen_ll_rNorm_top_pt_binning, std::vector<double>{gen_ll_rNorm,
            // gen_top_pt}, trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_llbar_delta_phi_top_pt,
            // gen_llbar_delta_phi_top_pt_binning,
            // std::vector<double>{gen_llbar_delta_phi, gen_top_pt},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_llbar_delta_eta_top_pt,
            // gen_llbar_delta_eta_top_pt_binning,
            // std::vector<double>{gen_llbar_delta_eta, gen_top_pt},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);

            // ******************************
            // top_scatteringangle_ttbarframe
            // ******************************

            // // Polarizations
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b1k_top_scatteringangle_ttbarframe,
            // gen_b1k_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1k, gen_top_scatteringangle_ttbarframe},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b2k_top_scatteringangle_ttbarframe,
            // gen_b2k_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2k, gen_top_scatteringangle_ttbarframe},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b1j_top_scatteringangle_ttbarframe,
            // gen_b1j_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1j, gen_top_scatteringangle_ttbarframe},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b2j_top_scatteringangle_ttbarframe,
            // gen_b2j_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2j, gen_top_scatteringangle_ttbarframe},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b1r_top_scatteringangle_ttbarframe,
            // gen_b1r_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1r, gen_top_scatteringangle_ttbarframe},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b2r_top_scatteringangle_ttbarframe,
            // gen_b2r_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2r, gen_top_scatteringangle_ttbarframe},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b1q_top_scatteringangle_ttbarframe,
            // gen_b1q_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1q, gen_top_scatteringangle_ttbarframe},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b2q_top_scatteringangle_ttbarframe,
            // gen_b2q_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2q, gen_top_scatteringangle_ttbarframe},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b1n_top_scatteringangle_ttbarframe,
            // gen_b1n_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b1n, gen_top_scatteringangle_ttbarframe},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_b2n_top_scatteringangle_ttbarframe,
            // gen_b2n_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_b2n, gen_top_scatteringangle_ttbarframe},
            // trueLevelWeight, eventWeight, symmetrizeSpinVar);

            // // Diagonal elements
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kk_top_scatteringangle_ttbarframe,
            // gen_c_kk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_kk,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rr_top_scatteringangle_ttbarframe,
            // gen_c_rr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rr,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nn_top_scatteringangle_ttbarframe,
            // gen_c_nn_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nn,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);

            // // Off diagonal elements
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rk_top_scatteringangle_ttbarframe,
            // gen_c_rk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rk,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kr_top_scatteringangle_ttbarframe,
            // gen_c_kr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_kr,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nr_top_scatteringangle_ttbarframe,
            // gen_c_nr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nr,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rn_top_scatteringangle_ttbarframe,
            // gen_c_rn_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rn,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nk_top_scatteringangle_ttbarframe,
            // gen_c_nk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nk,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kn_top_scatteringangle_ttbarframe,
            // gen_c_kn_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_kn,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);

            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Prk_top_scatteringangle_ttbarframe,
            // gen_c_Prk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rk + gen_c_kr,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Mrk_top_scatteringangle_ttbarframe,
            // gen_c_Mrk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rk - gen_c_kr,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Pnr_top_scatteringangle_ttbarframe,
            // gen_c_Pnr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nr + gen_c_rn,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Mnr_top_scatteringangle_ttbarframe,
            // gen_c_Mnr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nr - gen_c_rn,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Pnk_top_scatteringangle_ttbarframe,
            // gen_c_Pnk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nk + gen_c_kn,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Mnk_top_scatteringangle_ttbarframe,
            // gen_c_Mnk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_nk - gen_c_kn,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);

            // // Starred axes
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kj_top_scatteringangle_ttbarframe,
            // gen_c_kj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_kj,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rq_top_scatteringangle_ttbarframe,
            // gen_c_rq_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rq,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rj_top_scatteringangle_ttbarframe,
            // gen_c_rj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rj,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_jr_top_scatteringangle_ttbarframe,
            // gen_c_jr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_jr,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);

            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Prj_top_scatteringangle_ttbarframe,
            // gen_c_Prj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rj + gen_c_jr,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_Mrj_top_scatteringangle_ttbarframe,
            // gen_c_Mrj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_c_rj - gen_c_jr,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);

            // // New variables
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_han_top_scatteringangle_ttbarframe,
            // gen_c_han_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{ gen_c_kk - gen_c_rr - gen_c_nn,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_sca_top_scatteringangle_ttbarframe,
            // gen_c_sca_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_kk + gen_c_rr - gen_c_nn,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_tra_top_scatteringangle_ttbarframe,
            // gen_c_tra_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_kk - gen_c_rr + gen_c_nn,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_kjL_top_scatteringangle_ttbarframe,
            // gen_c_kjL_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_kj - gen_c_rr - gen_c_nn,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rqL_top_scatteringangle_ttbarframe,
            // gen_c_rqL_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_kk - gen_c_rq - gen_c_nn,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rkP_top_scatteringangle_ttbarframe,
            // gen_c_rkP_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_rk - gen_c_kr - gen_c_nn,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_rkM_top_scatteringangle_ttbarframe,
            // gen_c_rkM_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_rk + gen_c_kr - gen_c_nn,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nrP_top_scatteringangle_ttbarframe,
            // gen_c_nrP_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_nr - gen_c_rn - gen_c_kk,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nrM_top_scatteringangle_ttbarframe,
            // gen_c_nrM_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_nr + gen_c_rn - gen_c_kk,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nkP_top_scatteringangle_ttbarframe,
            // gen_c_nkP_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_nk - gen_c_kn - gen_c_rr,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_c_nkM_top_scatteringangle_ttbarframe,
            // gen_c_nkM_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{-gen_c_nk + gen_c_kn - gen_c_rr,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);

            // // Lab frame variables
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_ll_cHel_top_scatteringangle_ttbarframe,
            // gen_ll_cHel_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_ll_cHel,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_ll_cLab_top_scatteringangle_ttbarframe,
            // gen_ll_cLab_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_ll_cLab,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_ll_kNorm_top_scatteringangle_ttbarframe,
            // gen_ll_kNorm_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_ll_kNorm,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_ll_rNorm_top_scatteringangle_ttbarframe,
            // gen_ll_rNorm_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_ll_rNorm,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_llbar_delta_phi_top_scatteringangle_ttbarframe,
            // gen_llbar_delta_phi_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_llbar_delta_phi,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);
            // fillUnderOverFlowDeltaWeight(hrecoVsgen_llbar_delta_eta_top_scatteringangle_ttbarframe,
            // gen_llbar_delta_eta_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{gen_llbar_delta_eta,
            // gen_top_scatteringangle_ttbarframe}, trueLevelWeight,
            // eventWeight, symmetrizeSpinVar);

            // *****************************************
            // Fill Resolution matrices for 1D variables
            // *****************************************

            // Kinematic
            fillUnderOverFlow(hresolutionbins_top_pt, residual_top_pt_binning, gen_top_pt_binning,
                              std::vector<double>{top_pt - gen_top_pt}, std::vector<double>{gen_top_pt}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_l_pt, residual_l_pt_binning, gen_l_pt_binning,
                              std::vector<double>{l_pt - gen_l_pt}, std::vector<double>{gen_l_pt}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_lbar_pt, residual_lbar_pt_binning, gen_lbar_pt_binning,
                              std::vector<double>{lbar_pt - gen_lbar_pt}, std::vector<double>{gen_lbar_pt}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_ttbar_pt, residual_ttbar_pt_binning, gen_ttbar_pt_binning,
                              std::vector<double>{ttbar_pt - gen_ttbar_pt}, std::vector<double>{gen_ttbar_pt},
                              eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_ttbar_mass, residual_ttbar_mass_binning, gen_ttbar_mass_binning,
                              std::vector<double>{ttbar_mass - gen_ttbar_mass}, std::vector<double>{gen_ttbar_mass},
                              eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_ttbar_delta_phi, residual_ttbar_delta_phi_binning,
                              gen_ttbar_delta_phi_binning, std::vector<double>{ttbar_delta_phi - gen_ttbar_delta_phi},
                              std::vector<double>{gen_ttbar_delta_phi}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_ttbar_delta_eta, residual_ttbar_delta_eta_binning,
                              gen_ttbar_delta_eta_binning, std::vector<double>{ttbar_delta_eta - gen_ttbar_delta_eta},
                              std::vector<double>{gen_ttbar_delta_eta}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_ttbar_rapidity, residual_ttbar_rapidity_binning,
                              gen_ttbar_rapidity_binning, std::vector<double>{ttbar_rapidity - gen_ttbar_rapidity},
                              std::vector<double>{gen_ttbar_rapidity}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_llbar_pt, residual_llbar_pt_binning, gen_llbar_pt_binning,
                              std::vector<double>{llbar_pt - gen_llbar_pt}, std::vector<double>{gen_llbar_pt},
                              eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_llbar_mass, residual_llbar_mass_binning, gen_llbar_mass_binning,
                              std::vector<double>{llbar_mass - gen_llbar_mass}, std::vector<double>{gen_llbar_mass},
                              eventWeight, true, symmetrizeSpinVar);

            // Spin corr

            // Polarizations
            fillUnderOverFlow(hresolutionbins_b1k, residual_b1k_binning, gen_b1k_binning,
                              std::vector<double>{b1k - gen_b1k}, std::vector<double>{gen_b1k}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_b2k, residual_b2k_binning, gen_b2k_binning,
                              std::vector<double>{b2k - gen_b2k}, std::vector<double>{gen_b2k}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_b1j, residual_b1j_binning, gen_b1j_binning,
                              std::vector<double>{b1j - gen_b1j}, std::vector<double>{gen_b1j}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_b2j, residual_b2j_binning, gen_b2j_binning,
                              std::vector<double>{b2j - gen_b2j}, std::vector<double>{gen_b2j}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_b1r, residual_b1r_binning, gen_b1r_binning,
                              std::vector<double>{b1r - gen_b1r}, std::vector<double>{gen_b1r}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_b2r, residual_b2r_binning, gen_b2r_binning,
                              std::vector<double>{b2r - gen_b2r}, std::vector<double>{gen_b2r}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_b1q, residual_b1q_binning, gen_b1q_binning,
                              std::vector<double>{b1q - gen_b1q}, std::vector<double>{gen_b1q}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_b2q, residual_b2q_binning, gen_b2q_binning,
                              std::vector<double>{b2q - gen_b2q}, std::vector<double>{gen_b2q}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_b1n, residual_b1n_binning, gen_b1n_binning,
                              std::vector<double>{b1n - gen_b1n}, std::vector<double>{gen_b1n}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_b2n, residual_b2n_binning, gen_b2n_binning,
                              std::vector<double>{b2n - gen_b2n}, std::vector<double>{gen_b2n}, eventWeight, true,
                              symmetrizeSpinVar);

            // Diagonal elements
            fillUnderOverFlow(hresolutionbins_c_kk, residual_c_kk_binning, gen_c_kk_binning,
                              std::vector<double>{c_kk - gen_c_kk}, std::vector<double>{gen_c_kk}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_rr, residual_c_rr_binning, gen_c_rr_binning,
                              std::vector<double>{c_rr - gen_c_rr}, std::vector<double>{gen_c_rr}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_nn, residual_c_nn_binning, gen_c_nn_binning,
                              std::vector<double>{c_nn - gen_c_nn}, std::vector<double>{gen_c_nn}, eventWeight, true,
                              symmetrizeSpinVar);

            // Off diagonal elements
            fillUnderOverFlow(hresolutionbins_c_rk, residual_c_rk_binning, gen_c_rk_binning,
                              std::vector<double>{c_rk - gen_c_rk}, std::vector<double>{gen_c_rk}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_kr, residual_c_kr_binning, gen_c_kr_binning,
                              std::vector<double>{c_kr - gen_c_kr}, std::vector<double>{gen_c_kr}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_nr, residual_c_nr_binning, gen_c_nr_binning,
                              std::vector<double>{c_nr - gen_c_nr}, std::vector<double>{gen_c_nr}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_rn, residual_c_rn_binning, gen_c_rn_binning,
                              std::vector<double>{c_rn - gen_c_rn}, std::vector<double>{gen_c_rn}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_nk, residual_c_nk_binning, gen_c_nk_binning,
                              std::vector<double>{c_nk - gen_c_nk}, std::vector<double>{gen_c_nk}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_kn, residual_c_kn_binning, gen_c_kn_binning,
                              std::vector<double>{c_kn - gen_c_kn}, std::vector<double>{gen_c_kn}, eventWeight, true,
                              symmetrizeSpinVar);

            fillUnderOverFlow(hresolutionbins_c_Prk, residual_c_Prk_binning, gen_c_Prk_binning,
                              std::vector<double>{(c_rk + c_kr) - (gen_c_rk + gen_c_kr)},
                              std::vector<double>{gen_c_rk + gen_c_kr}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_Mrk, residual_c_Mrk_binning, gen_c_Mrk_binning,
                              std::vector<double>{(c_rk - c_kr) - (gen_c_rk - gen_c_kr)},
                              std::vector<double>{gen_c_rk - gen_c_kr}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_Pnr, residual_c_Pnr_binning, gen_c_Pnr_binning,
                              std::vector<double>{(c_nr + c_rn) - (gen_c_nr + gen_c_rn)},
                              std::vector<double>{gen_c_nr + gen_c_rn}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_Mnr, residual_c_Mnr_binning, gen_c_Mnr_binning,
                              std::vector<double>{(c_nr - c_rn) - (gen_c_nr - gen_c_rn)},
                              std::vector<double>{gen_c_nr - gen_c_rn}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_Pnk, residual_c_Pnk_binning, gen_c_Pnk_binning,
                              std::vector<double>{(c_nk + c_kn) - (gen_c_nk + gen_c_kn)},
                              std::vector<double>{gen_c_nk + gen_c_kn}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_Mnk, residual_c_Mnk_binning, gen_c_Mnk_binning,
                              std::vector<double>{(c_nk - c_kn) - (gen_c_nk - gen_c_kn)},
                              std::vector<double>{gen_c_nk - gen_c_kn}, eventWeight, true, symmetrizeSpinVar);

            // Starred variables
            fillUnderOverFlow(hresolutionbins_c_kj, residual_c_kj_binning, gen_c_kj_binning,
                              std::vector<double>{c_kj - gen_c_kj}, std::vector<double>{gen_c_kj}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_rq, residual_c_rq_binning, gen_c_rq_binning,
                              std::vector<double>{c_rq - gen_c_rq}, std::vector<double>{gen_c_rq}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_rj, residual_c_rj_binning, gen_c_rj_binning,
                              std::vector<double>{c_rj - gen_c_rj}, std::vector<double>{gen_c_rj}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_jr, residual_c_jr_binning, gen_c_jr_binning,
                              std::vector<double>{c_jr - gen_c_jr}, std::vector<double>{gen_c_jr}, eventWeight, true,
                              symmetrizeSpinVar);

            fillUnderOverFlow(hresolutionbins_c_Prj, residual_c_Prj_binning, gen_c_Prj_binning,
                              std::vector<double>{(c_rj + c_jr) - (gen_c_rj + gen_c_jr)},
                              std::vector<double>{gen_c_rj + gen_c_jr}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_Mrj, residual_c_Mrj_binning, gen_c_Mrj_binning,
                              std::vector<double>{(c_rj - c_jr) - (gen_c_rj - gen_c_jr)},
                              std::vector<double>{gen_c_rj - gen_c_jr}, eventWeight, true, symmetrizeSpinVar);

            // New variables
            fillUnderOverFlow(hresolutionbins_c_han, residual_c_han_binning, gen_c_han_binning,
                              std::vector<double>{(gen_c_kk - gen_c_rr - gen_c_nn) - (c_kk - c_rr - c_nn)},
                              std::vector<double>{gen_c_kk - gen_c_rr - gen_c_nn}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_sca, residual_c_sca_binning, gen_c_sca_binning,
                              std::vector<double>{(-gen_c_kk + gen_c_rr - gen_c_nn) - (-c_kk + c_rr - c_nn)},
                              std::vector<double>{-gen_c_kk + gen_c_rr - gen_c_nn}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_tra, residual_c_tra_binning, gen_c_tra_binning,
                              std::vector<double>{(-gen_c_kk - gen_c_rr + gen_c_nn) - (-c_kk - c_rr + c_nn)},
                              std::vector<double>{-gen_c_kk - gen_c_rr + gen_c_nn}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_kjL, residual_c_kjL_binning, gen_c_kjL_binning,
                              std::vector<double>{(-gen_c_kj - gen_c_rr - gen_c_nn) - (-c_kj - c_rr - c_nn)},
                              std::vector<double>{-gen_c_kj - gen_c_rr - gen_c_nn}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_rqL, residual_c_rqL_binning, gen_c_rqL_binning,
                              std::vector<double>{(-gen_c_kk - gen_c_rq - gen_c_nn) - (-c_kk - c_rq - c_nn)},
                              std::vector<double>{-gen_c_kk - gen_c_rq - gen_c_nn}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_rkP, residual_c_rkP_binning, gen_c_rkP_binning,
                              std::vector<double>{(-gen_c_rk - gen_c_kr - gen_c_nn) - (-c_rk - c_kr - c_nn)},
                              std::vector<double>{-gen_c_rk - gen_c_kr - gen_c_nn}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_rkM, residual_c_rkM_binning, gen_c_rkM_binning,
                              std::vector<double>{(-gen_c_rk + gen_c_kr - gen_c_nn) - (-c_rk + c_kr - c_nn)},
                              std::vector<double>{-gen_c_rk + gen_c_kr - gen_c_nn}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_nrP, residual_c_nrP_binning, gen_c_nrP_binning,
                              std::vector<double>{(-gen_c_nr - gen_c_rn - gen_c_kk) - (-c_nr - c_rn - c_kk)},
                              std::vector<double>{-gen_c_nr - gen_c_rn - gen_c_kk}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_nrM, residual_c_nrM_binning, gen_c_nrM_binning,
                              std::vector<double>{(-gen_c_nr + gen_c_rn - gen_c_kk) - (-c_nr + c_rn - c_kk)},
                              std::vector<double>{-gen_c_nr + gen_c_rn - gen_c_kk}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_nkP, residual_c_nkP_binning, gen_c_nkP_binning,
                              std::vector<double>{(-gen_c_nk - gen_c_kn - gen_c_rr) - (-c_nk - c_kn - c_rr)},
                              std::vector<double>{-gen_c_nk - gen_c_kn - gen_c_rr}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_nkM, residual_c_nkM_binning, gen_c_nkM_binning,
                              std::vector<double>{(-gen_c_nk + gen_c_kn - gen_c_rr) - (-c_nk + c_kn - c_rr)},
                              std::vector<double>{-gen_c_nk + gen_c_kn - gen_c_rr}, eventWeight, true,
                              symmetrizeSpinVar);

            // Lab frame variables
            fillUnderOverFlow(hresolutionbins_ll_cHel, residual_ll_cHel_binning, gen_ll_cHel_binning,
                              std::vector<double>{ll_cHel - gen_ll_cHel}, std::vector<double>{gen_ll_cHel}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_ll_cLab, residual_ll_cLab_binning, gen_ll_cLab_binning,
                              std::vector<double>{ll_cLab - gen_ll_cLab}, std::vector<double>{gen_ll_cLab}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_ll_kNorm, residual_ll_kNorm_binning, gen_ll_kNorm_binning,
                              std::vector<double>{ll_kNorm - gen_ll_kNorm}, std::vector<double>{gen_ll_kNorm},
                              eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_ll_rNorm, residual_ll_rNorm_binning, gen_ll_rNorm_binning,
                              std::vector<double>{ll_rNorm - gen_ll_rNorm}, std::vector<double>{gen_ll_rNorm},
                              eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_llbar_delta_phi, residual_llbar_delta_phi_binning,
                              gen_llbar_delta_phi_binning, std::vector<double>{llbar_delta_phi - gen_llbar_delta_phi},
                              std::vector<double>{gen_llbar_delta_phi}, eventWeight, true, symmetrizeSpinVar);

            // *****************************************
            // Fill Resolution matrices for 2D variables
            // *****************************************

            // Amandeep : 2D resolution matrices
            // fillUnderOverFlow(hresolutionbins_c_kk_mttbar,
            // residual_c_kk_mttbar_binning, gen_c_kk_mttbar_binning,
            // std::vector<double>{c_kk - gen_c_kk},
            // std::vector<double>{gen_c_kk, gen_ttbar_mass}, eventWeight, true,
            // symmetrizeSpinVar);

            // ******
            // mttbar
            // ******

            // Polarizations

            fillUnderOverFlow(hresolutionbins_b1k_exjet, residual_b1k_exjet_binning, gen_b1k_exjet_binning,
                              std::vector<double>{b1k - gen_b1k}, std::vector<double>{gen_b1k, gen_n_extraJets_iso08},
                              eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_b2k_exjet, residual_b2k_exjet_binning, gen_b2k_exjet_binning,
                              std::vector<double>{b2k - gen_b2k}, std::vector<double>{gen_b2k, gen_n_extraJets_iso08},
                              eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_b1j_exjet, residual_b1j_exjet_binning, gen_b1j_exjet_binning,
                              std::vector<double>{b1j - gen_b1j}, std::vector<double>{gen_b1j, gen_n_extraJets_iso08},
                              eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_b2j_exjet, residual_b2j_exjet_binning, gen_b2j_exjet_binning,
                              std::vector<double>{b2j - gen_b2j}, std::vector<double>{gen_b2j, gen_n_extraJets_iso08},
                              eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_b1r_exjet, residual_b1r_exjet_binning, gen_b1r_exjet_binning,
                              std::vector<double>{b1r - gen_b1r}, std::vector<double>{gen_b1r, gen_n_extraJets_iso08},
                              eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_b2r_exjet, residual_b2r_exjet_binning, gen_b2r_exjet_binning,
                              std::vector<double>{b2r - gen_b2r}, std::vector<double>{gen_b2r, gen_n_extraJets_iso08},
                              eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_b1q_exjet, residual_b1q_exjet_binning, gen_b1q_exjet_binning,
                              std::vector<double>{b1q - gen_b1q}, std::vector<double>{gen_b1q, gen_n_extraJets_iso08},
                              eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_b2q_exjet, residual_b2q_exjet_binning, gen_b2q_exjet_binning,
                              std::vector<double>{b2q - gen_b2q}, std::vector<double>{gen_b2q, gen_n_extraJets_iso08},
                              eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_b1n_exjet, residual_b1n_exjet_binning, gen_b1n_exjet_binning,
                              std::vector<double>{b1n - gen_b1n}, std::vector<double>{gen_b1n, gen_n_extraJets_iso08},
                              eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_b2n_exjet, residual_b2n_exjet_binning, gen_b2n_exjet_binning,
                              std::vector<double>{b2n - gen_b2n}, std::vector<double>{gen_b2n, gen_n_extraJets_iso08},
                              eventWeight, true, symmetrizeSpinVar);

            // Diagonal elements
            fillUnderOverFlow(hresolutionbins_c_kk_exjet, residual_c_kk_exjet_binning, gen_c_kk_exjet_binning,
                              std::vector<double>{c_kk - gen_c_kk},
                              std::vector<double>{gen_c_kk, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_rr_exjet, residual_c_rr_exjet_binning, gen_c_rr_exjet_binning,
                              std::vector<double>{c_rr - gen_c_rr},
                              std::vector<double>{gen_c_rr, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_nn_exjet, residual_c_nn_exjet_binning, gen_c_nn_exjet_binning,
                              std::vector<double>{c_nn - gen_c_nn},
                              std::vector<double>{gen_c_nn, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);

            // Off diagonal elements
            fillUnderOverFlow(hresolutionbins_c_rk_exjet, residual_c_rk_exjet_binning, gen_c_rk_exjet_binning,
                              std::vector<double>{c_rk - gen_c_rk},
                              std::vector<double>{gen_c_rk, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_kr_exjet, residual_c_kr_exjet_binning, gen_c_kr_exjet_binning,
                              std::vector<double>{c_kr - gen_c_kr},
                              std::vector<double>{gen_c_kr, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_nr_exjet, residual_c_nr_exjet_binning, gen_c_nr_exjet_binning,
                              std::vector<double>{c_nr - gen_c_nr},
                              std::vector<double>{gen_c_nr, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_rn_exjet, residual_c_rn_exjet_binning, gen_c_rn_exjet_binning,
                              std::vector<double>{c_rn - gen_c_rn},
                              std::vector<double>{gen_c_rn, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_nk_exjet, residual_c_nk_exjet_binning, gen_c_nk_exjet_binning,
                              std::vector<double>{c_nk - gen_c_nk},
                              std::vector<double>{gen_c_nk, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_kn_exjet, residual_c_kn_exjet_binning, gen_c_kn_exjet_binning,
                              std::vector<double>{c_kn - gen_c_kn},
                              std::vector<double>{gen_c_kn, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);

            fillUnderOverFlow(hresolutionbins_c_Prk_exjet, residual_c_Prk_exjet_binning, gen_c_Prk_exjet_binning,
                              std::vector<double>{(c_rk + c_kr) - (gen_c_rk + gen_c_kr)},
                              std::vector<double>{gen_c_rk + gen_c_kr, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_Mrk_exjet, residual_c_Mrk_exjet_binning, gen_c_Mrk_exjet_binning,
                              std::vector<double>{(c_rk - c_kr) - (gen_c_rk - gen_c_kr)},
                              std::vector<double>{gen_c_rk - gen_c_kr, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_Pnr_exjet, residual_c_Pnr_exjet_binning, gen_c_Pnr_exjet_binning,
                              std::vector<double>{(c_nr + c_rn) - (gen_c_nr + gen_c_rn)},
                              std::vector<double>{gen_c_nr + gen_c_rn, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_Mnr_exjet, residual_c_Mnr_exjet_binning, gen_c_Mnr_exjet_binning,
                              std::vector<double>{(c_nr - c_rn) - (gen_c_nr - gen_c_rn)},
                              std::vector<double>{gen_c_nr - gen_c_rn, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_Pnk_exjet, residual_c_Pnk_exjet_binning, gen_c_Pnk_exjet_binning,
                              std::vector<double>{(c_nk + c_kn) - (gen_c_nk + gen_c_kn)},
                              std::vector<double>{gen_c_nk + gen_c_kn, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_Mnk_exjet, residual_c_Mnk_exjet_binning, gen_c_Mnk_exjet_binning,
                              std::vector<double>{(c_nk - c_kn) - (gen_c_nk - gen_c_kn)},
                              std::vector<double>{gen_c_nk - gen_c_kn, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);

            // Starred axes
            fillUnderOverFlow(hresolutionbins_c_kj_exjet, residual_c_kj_exjet_binning, gen_c_kj_exjet_binning,
                              std::vector<double>{c_kj - gen_c_kj},
                              std::vector<double>{gen_c_kj, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_rq_exjet, residual_c_rq_exjet_binning, gen_c_rq_exjet_binning,
                              std::vector<double>{c_rq - gen_c_rq},
                              std::vector<double>{gen_c_rq, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_rj_exjet, residual_c_rj_exjet_binning, gen_c_rj_exjet_binning,
                              std::vector<double>{c_rj - gen_c_rj},
                              std::vector<double>{gen_c_rj, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_jr_exjet, residual_c_jr_exjet_binning, gen_c_jr_exjet_binning,
                              std::vector<double>{c_jr - gen_c_jr},
                              std::vector<double>{gen_c_jr, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);

            fillUnderOverFlow(hresolutionbins_c_Prj_exjet, residual_c_Prj_exjet_binning, gen_c_Prj_exjet_binning,
                              std::vector<double>{(c_rj + c_jr) - (gen_c_rj + gen_c_jr)},
                              std::vector<double>{gen_c_rj + gen_c_jr, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_Mrj_exjet, residual_c_Mrj_exjet_binning, gen_c_Mrj_exjet_binning,
                              std::vector<double>{(c_rj - c_jr) - (gen_c_rj - gen_c_jr)},
                              std::vector<double>{gen_c_rj - gen_c_jr, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);

            // New variables
            fillUnderOverFlow(hresolutionbins_c_han_exjet, residual_c_han_exjet_binning, gen_c_han_exjet_binning,
                              std::vector<double>{(gen_c_kk - gen_c_rr - gen_c_nn) - (c_kk - c_rr - c_nn)},
                              std::vector<double>{gen_c_kk - gen_c_rr - gen_c_nn, gen_n_extraJets_iso08}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_sca_exjet, residual_c_sca_exjet_binning, gen_c_sca_exjet_binning,
                              std::vector<double>{(-gen_c_kk + gen_c_rr - gen_c_nn) - (-c_kk + c_rr - c_nn)},
                              std::vector<double>{-gen_c_kk + gen_c_rr - gen_c_nn, gen_n_extraJets_iso08}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_tra_exjet, residual_c_tra_exjet_binning, gen_c_tra_exjet_binning,
                              std::vector<double>{(-gen_c_kk - gen_c_rr + gen_c_nn) - (-c_kk - c_rr + c_nn)},
                              std::vector<double>{-gen_c_kk - gen_c_rr + gen_c_nn, gen_n_extraJets_iso08}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_kjL_exjet, residual_c_kjL_exjet_binning, gen_c_kjL_exjet_binning,
                              std::vector<double>{(-gen_c_kj - gen_c_rr - gen_c_nn) - (-c_kj - c_rr - c_nn)},
                              std::vector<double>{-gen_c_kj - gen_c_rr - gen_c_nn, gen_n_extraJets_iso08}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_rqL_exjet, residual_c_rqL_exjet_binning, gen_c_rqL_exjet_binning,
                              std::vector<double>{(-gen_c_kk - gen_c_rq - gen_c_nn) - (-c_kk - c_rq - c_nn)},
                              std::vector<double>{-gen_c_kk - gen_c_rq - gen_c_nn, gen_n_extraJets_iso08}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_rkP_exjet, residual_c_rkP_exjet_binning, gen_c_rkP_exjet_binning,
                              std::vector<double>{(-gen_c_rk - gen_c_kr - gen_c_nn) - (-c_rk - c_kr - c_nn)},
                              std::vector<double>{-gen_c_rk - gen_c_kr - gen_c_nn, gen_n_extraJets_iso08}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_rkM_exjet, residual_c_rkM_exjet_binning, gen_c_rkM_exjet_binning,
                              std::vector<double>{(-gen_c_rk + gen_c_kr - gen_c_nn) - (-c_rk + c_kr - c_nn)},
                              std::vector<double>{-gen_c_rk + gen_c_kr - gen_c_nn, gen_n_extraJets_iso08}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_nrP_exjet, residual_c_nrP_exjet_binning, gen_c_nrP_exjet_binning,
                              std::vector<double>{(-gen_c_nr - gen_c_rn - gen_c_kk) - (-c_nr - c_rn - c_kk)},
                              std::vector<double>{-gen_c_nr - gen_c_rn - gen_c_kk, gen_n_extraJets_iso08}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_nrM_exjet, residual_c_nrM_exjet_binning, gen_c_nrM_exjet_binning,
                              std::vector<double>{(-gen_c_nr + gen_c_rn - gen_c_kk) - (-c_nr + c_rn - c_kk)},
                              std::vector<double>{-gen_c_nr + gen_c_rn - gen_c_kk, gen_n_extraJets_iso08}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_nkP_exjet, residual_c_nkP_exjet_binning, gen_c_nkP_exjet_binning,
                              std::vector<double>{(-gen_c_nk - gen_c_kn - gen_c_rr) - (-c_nk - c_kn - c_rr)},
                              std::vector<double>{-gen_c_nk - gen_c_kn - gen_c_rr, gen_n_extraJets_iso08}, eventWeight,
                              true, symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_c_nkM_exjet, residual_c_nkM_exjet_binning, gen_c_nkM_exjet_binning,
                              std::vector<double>{(-gen_c_nk + gen_c_kn - gen_c_rr) - (-c_nk + c_kn - c_rr)},
                              std::vector<double>{-gen_c_nk + gen_c_kn - gen_c_rr, gen_n_extraJets_iso08}, eventWeight,
                              true, symmetrizeSpinVar);

            // Lab frame variables
            fillUnderOverFlow(hresolutionbins_ll_cHel_exjet, residual_ll_cHel_exjet_binning, gen_ll_cHel_exjet_binning,
                              std::vector<double>{ll_cHel - gen_ll_cHel},
                              std::vector<double>{gen_ll_cHel, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_ll_cLab_exjet, residual_ll_cLab_exjet_binning, gen_ll_cLab_exjet_binning,
                              std::vector<double>{ll_cLab - gen_ll_cLab},
                              std::vector<double>{gen_ll_cLab, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_ll_kNorm_exjet, residual_ll_kNorm_exjet_binning,
                              gen_ll_kNorm_exjet_binning, std::vector<double>{ll_kNorm - gen_ll_kNorm},
                              std::vector<double>{gen_ll_kNorm, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(hresolutionbins_ll_rNorm_exjet, residual_ll_rNorm_exjet_binning,
                              gen_ll_rNorm_exjet_binning, std::vector<double>{ll_rNorm - gen_ll_rNorm},
                              std::vector<double>{gen_ll_rNorm, gen_n_extraJets_iso08}, eventWeight, true,
                              symmetrizeSpinVar);
            fillUnderOverFlow(
                hresolutionbins_llbar_delta_phi_exjet, residual_llbar_delta_phi_exjet_binning,
                gen_llbar_delta_phi_exjet_binning, std::vector<double>{llbar_delta_phi - gen_llbar_delta_phi},
                std::vector<double>{gen_llbar_delta_phi, gen_n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);
            fillUnderOverFlow(
                hresolutionbins_llbar_delta_eta_exjet, residual_llbar_delta_eta_exjet_binning,
                gen_llbar_delta_eta_exjet_binning, std::vector<double>{llbar_delta_eta - gen_llbar_delta_eta},
                std::vector<double>{gen_llbar_delta_eta, gen_n_extraJets_iso08}, eventWeight, true, symmetrizeSpinVar);

            // ******
            // top_pt
            // ******

            // // Polarizations
            // fillUnderOverFlow(hresolutionbins_b1k_top_pt,
            // residual_b1k_top_pt_binning, gen_b1k_top_pt_binning,
            // std::vector<double>{b1k - gen_b1k}, std::vector<double>{gen_b1k,
            // gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_b2k_top_pt,
            // residual_b2k_top_pt_binning, gen_b2k_top_pt_binning,
            // std::vector<double>{b2k - gen_b2k}, std::vector<double>{gen_b2k,
            // gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_b1j_top_pt,
            // residual_b1j_top_pt_binning, gen_b1j_top_pt_binning,
            // std::vector<double>{b1j - gen_b1j}, std::vector<double>{gen_b1j,
            // gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_b2j_top_pt,
            // residual_b2j_top_pt_binning, gen_b2j_top_pt_binning,
            // std::vector<double>{b2j - gen_b2j}, std::vector<double>{gen_b2j,
            // gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_b1r_top_pt,
            // residual_b1r_top_pt_binning, gen_b1r_top_pt_binning,
            // std::vector<double>{b1r - gen_b1r}, std::vector<double>{gen_b1r,
            // gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_b2r_top_pt,
            // residual_b2r_top_pt_binning, gen_b2r_top_pt_binning,
            // std::vector<double>{b2r - gen_b2r}, std::vector<double>{gen_b2r,
            // gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_b1q_top_pt,
            // residual_b1q_top_pt_binning, gen_b1q_top_pt_binning,
            // std::vector<double>{b1q - gen_b1q}, std::vector<double>{gen_b1q,
            // gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_b2q_top_pt,
            // residual_b2q_top_pt_binning, gen_b2q_top_pt_binning,
            // std::vector<double>{b2q - gen_b2q}, std::vector<double>{gen_b2q,
            // gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_b1n_top_pt,
            // residual_b1n_top_pt_binning, gen_b1n_top_pt_binning,
            // std::vector<double>{b1n - gen_b1n}, std::vector<double>{gen_b1n,
            // gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_b2n_top_pt,
            // residual_b2n_top_pt_binning, gen_b2n_top_pt_binning,
            // std::vector<double>{b2n - gen_b2n}, std::vector<double>{gen_b2n,
            // gen_top_pt}, eventWeight, true, symmetrizeSpinVar);

            // // Diagonal elements
            // fillUnderOverFlow(hresolutionbins_c_kk_top_pt,
            // residual_c_kk_top_pt_binning, gen_c_kk_top_pt_binning,
            // std::vector<double>{c_kk - gen_c_kk},
            // std::vector<double>{gen_c_kk, gen_top_pt}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_rr_top_pt,
            // residual_c_rr_top_pt_binning, gen_c_rr_top_pt_binning,
            // std::vector<double>{c_rr - gen_c_rr},
            // std::vector<double>{gen_c_rr, gen_top_pt}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_nn_top_pt,
            // residual_c_nn_top_pt_binning, gen_c_nn_top_pt_binning,
            // std::vector<double>{c_nn - gen_c_nn},
            // std::vector<double>{gen_c_nn, gen_top_pt}, eventWeight, true,
            // symmetrizeSpinVar);

            // // Off diagonal elements
            // fillUnderOverFlow(hresolutionbins_c_rk_top_pt,
            // residual_c_rk_top_pt_binning, gen_c_rk_top_pt_binning,
            // std::vector<double>{c_rk - gen_c_rk},
            // std::vector<double>{gen_c_rk, gen_top_pt}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_kr_top_pt,
            // residual_c_kr_top_pt_binning, gen_c_kr_top_pt_binning,
            // std::vector<double>{c_kr - gen_c_kr},
            // std::vector<double>{gen_c_kr, gen_top_pt}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_nr_top_pt,
            // residual_c_nr_top_pt_binning, gen_c_nr_top_pt_binning,
            // std::vector<double>{c_nr - gen_c_nr},
            // std::vector<double>{gen_c_nr, gen_top_pt}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_rn_top_pt,
            // residual_c_rn_top_pt_binning, gen_c_rn_top_pt_binning,
            // std::vector<double>{c_rn - gen_c_rn},
            // std::vector<double>{gen_c_rn, gen_top_pt}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_nk_top_pt,
            // residual_c_nk_top_pt_binning, gen_c_nk_top_pt_binning,
            // std::vector<double>{c_nk - gen_c_nk},
            // std::vector<double>{gen_c_nk, gen_top_pt}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_kn_top_pt,
            // residual_c_kn_top_pt_binning, gen_c_kn_top_pt_binning,
            // std::vector<double>{c_kn - gen_c_kn},
            // std::vector<double>{gen_c_kn, gen_top_pt}, eventWeight, true,
            // symmetrizeSpinVar);

            // fillUnderOverFlow(hresolutionbins_c_Prk_top_pt,
            // residual_c_Prk_top_pt_binning, gen_c_Prk_top_pt_binning,
            // std::vector<double>{(c_rk + c_kr) - (gen_c_rk + gen_c_kr)},
            // std::vector<double>{gen_c_rk + gen_c_kr, gen_top_pt},
            // eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_Mrk_top_pt,
            // residual_c_Mrk_top_pt_binning, gen_c_Mrk_top_pt_binning,
            // std::vector<double>{(c_rk - c_kr) - (gen_c_rk - gen_c_kr)},
            // std::vector<double>{gen_c_rk - gen_c_kr, gen_top_pt},
            // eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_Pnr_top_pt,
            // residual_c_Pnr_top_pt_binning, gen_c_Pnr_top_pt_binning,
            // std::vector<double>{(c_nr + c_rn) - (gen_c_nr + gen_c_rn)},
            // std::vector<double>{gen_c_nr + gen_c_rn, gen_top_pt},
            // eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_Mnr_top_pt,
            // residual_c_Mnr_top_pt_binning, gen_c_Mnr_top_pt_binning,
            // std::vector<double>{(c_nr - c_rn) - (gen_c_nr - gen_c_rn)},
            // std::vector<double>{gen_c_nr - gen_c_rn, gen_top_pt},
            // eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_Pnk_top_pt,
            // residual_c_Pnk_top_pt_binning, gen_c_Pnk_top_pt_binning,
            // std::vector<double>{(c_nk + c_kn) - (gen_c_nk + gen_c_kn)},
            // std::vector<double>{gen_c_nk + gen_c_kn, gen_top_pt},
            // eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_Mnk_top_pt,
            // residual_c_Mnk_top_pt_binning, gen_c_Mnk_top_pt_binning,
            // std::vector<double>{(c_nk - c_kn) - (gen_c_nk - gen_c_kn)},
            // std::vector<double>{gen_c_nk - gen_c_kn, gen_top_pt},
            // eventWeight, true, symmetrizeSpinVar);

            // // Starred axes
            // fillUnderOverFlow(hresolutionbins_c_kj_top_pt,
            // residual_c_kj_top_pt_binning, gen_c_kj_top_pt_binning,
            // std::vector<double>{c_kj - gen_c_kj},
            // std::vector<double>{gen_c_kj, gen_top_pt}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_rq_top_pt,
            // residual_c_rq_top_pt_binning, gen_c_rq_top_pt_binning,
            // std::vector<double>{c_rq - gen_c_rq},
            // std::vector<double>{gen_c_rq, gen_top_pt}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_rj_top_pt,
            // residual_c_rj_top_pt_binning, gen_c_rj_top_pt_binning,
            // std::vector<double>{c_rj - gen_c_rj},
            // std::vector<double>{gen_c_rj, gen_top_pt}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_jr_top_pt,
            // residual_c_jr_top_pt_binning, gen_c_jr_top_pt_binning,
            // std::vector<double>{c_jr - gen_c_jr},
            // std::vector<double>{gen_c_jr, gen_top_pt}, eventWeight, true,
            // symmetrizeSpinVar);

            // fillUnderOverFlow(hresolutionbins_c_Prj_top_pt,
            // residual_c_Prj_top_pt_binning, gen_c_Prj_top_pt_binning,
            // std::vector<double>{(c_rj + c_jr) - (gen_c_rj + gen_c_jr)},
            // std::vector<double>{gen_c_rj + gen_c_jr, gen_top_pt},
            // eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_Mrj_top_pt,
            // residual_c_Mrj_top_pt_binning, gen_c_Mrj_top_pt_binning,
            // std::vector<double>{(c_rj - c_jr) - (gen_c_rj - gen_c_jr)},
            // std::vector<double>{gen_c_rj - gen_c_jr, gen_top_pt},
            // eventWeight, true, symmetrizeSpinVar);

            // // New variables
            // fillUnderOverFlow(hresolutionbins_c_han_top_pt,
            // residual_c_han_top_pt_binning, gen_c_han_top_pt_binning,
            // std::vector<double>{( gen_c_kk - gen_c_rr - gen_c_nn) - ( c_kk -
            // c_rr - c_nn)}, std::vector<double>{ gen_c_kk - gen_c_rr -
            // gen_c_nn, gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_sca_top_pt,
            // residual_c_sca_top_pt_binning, gen_c_sca_top_pt_binning,
            // std::vector<double>{(-gen_c_kk + gen_c_rr - gen_c_nn) - (-c_kk +
            // c_rr - c_nn)}, std::vector<double>{-gen_c_kk + gen_c_rr -
            // gen_c_nn, gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_tra_top_pt,
            // residual_c_tra_top_pt_binning, gen_c_tra_top_pt_binning,
            // std::vector<double>{(-gen_c_kk - gen_c_rr + gen_c_nn) - (-c_kk -
            // c_rr + c_nn)}, std::vector<double>{-gen_c_kk - gen_c_rr +
            // gen_c_nn, gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_kjL_top_pt,
            // residual_c_kjL_top_pt_binning, gen_c_kjL_top_pt_binning,
            // std::vector<double>{(-gen_c_kj - gen_c_rr - gen_c_nn) - (-c_kj -
            // c_rr - c_nn)}, std::vector<double>{-gen_c_kj - gen_c_rr -
            // gen_c_nn, gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_rqL_top_pt,
            // residual_c_rqL_top_pt_binning, gen_c_rqL_top_pt_binning,
            // std::vector<double>{(-gen_c_kk - gen_c_rq - gen_c_nn) - (-c_kk -
            // c_rq - c_nn)}, std::vector<double>{-gen_c_kk - gen_c_rq -
            // gen_c_nn, gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_rkP_top_pt,
            // residual_c_rkP_top_pt_binning, gen_c_rkP_top_pt_binning,
            // std::vector<double>{(-gen_c_rk - gen_c_kr - gen_c_nn) - (-c_rk -
            // c_kr - c_nn)}, std::vector<double>{-gen_c_rk - gen_c_kr -
            // gen_c_nn, gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_rkM_top_pt,
            // residual_c_rkM_top_pt_binning, gen_c_rkM_top_pt_binning,
            // std::vector<double>{(-gen_c_rk + gen_c_kr - gen_c_nn) - (-c_rk +
            // c_kr - c_nn)}, std::vector<double>{-gen_c_rk + gen_c_kr -
            // gen_c_nn, gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_nrP_top_pt,
            // residual_c_nrP_top_pt_binning, gen_c_nrP_top_pt_binning,
            // std::vector<double>{(-gen_c_nr - gen_c_rn - gen_c_kk) - (-c_nr -
            // c_rn - c_kk)}, std::vector<double>{-gen_c_nr - gen_c_rn -
            // gen_c_kk, gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_nrM_top_pt,
            // residual_c_nrM_top_pt_binning, gen_c_nrM_top_pt_binning,
            // std::vector<double>{(-gen_c_nr + gen_c_rn - gen_c_kk) - (-c_nr +
            // c_rn - c_kk)}, std::vector<double>{-gen_c_nr + gen_c_rn -
            // gen_c_kk, gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_nkP_top_pt,
            // residual_c_nkP_top_pt_binning, gen_c_nkP_top_pt_binning,
            // std::vector<double>{(-gen_c_nk - gen_c_kn - gen_c_rr) - (-c_nk -
            // c_kn - c_rr)}, std::vector<double>{-gen_c_nk - gen_c_kn -
            // gen_c_rr, gen_top_pt}, eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_nkM_top_pt,
            // residual_c_nkM_top_pt_binning, gen_c_nkM_top_pt_binning,
            // std::vector<double>{(-gen_c_nk + gen_c_kn - gen_c_rr) - (-c_nk +
            // c_kn - c_rr)}, std::vector<double>{-gen_c_nk + gen_c_kn -
            // gen_c_rr, gen_top_pt}, eventWeight, true, symmetrizeSpinVar);

            // // Lab frame variables
            // fillUnderOverFlow(hresolutionbins_ll_cHel_top_pt,
            // residual_ll_cHel_top_pt_binning, gen_ll_cHel_top_pt_binning,
            // std::vector<double>{ll_cHel - gen_ll_cHel},
            // std::vector<double>{gen_ll_cHel, gen_top_pt}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_ll_cLab_top_pt,
            // residual_ll_cLab_top_pt_binning, gen_ll_cLab_top_pt_binning,
            // std::vector<double>{ll_cLab - gen_ll_cLab},
            // std::vector<double>{gen_ll_cLab, gen_top_pt}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_ll_kNorm_top_pt,
            // residual_ll_kNorm_top_pt_binning, gen_ll_kNorm_top_pt_binning,
            // std::vector<double>{ll_kNorm - gen_ll_kNorm},
            // std::vector<double>{gen_ll_kNorm, gen_top_pt}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_ll_rNorm_top_pt,
            // residual_ll_rNorm_top_pt_binning, gen_ll_rNorm_top_pt_binning,
            // std::vector<double>{ll_rNorm - gen_ll_rNorm},
            // std::vector<double>{gen_ll_rNorm, gen_top_pt}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_llbar_delta_phi_top_pt,
            // residual_llbar_delta_phi_top_pt_binning,
            // gen_llbar_delta_phi_top_pt_binning,
            // std::vector<double>{llbar_delta_phi - gen_llbar_delta_phi},
            // std::vector<double>{gen_llbar_delta_phi, gen_top_pt},
            // eventWeight, true, symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_llbar_delta_eta_top_pt,
            // residual_llbar_delta_eta_top_pt_binning,
            // gen_llbar_delta_eta_top_pt_binning,
            // std::vector<double>{llbar_delta_eta - gen_llbar_delta_eta},
            // std::vector<double>{gen_llbar_delta_eta, gen_top_pt},
            // eventWeight, true, symmetrizeSpinVar);

            // ******************************
            // top_scatteringangle_ttbarframe
            // ******************************

            // // Polarizations
            // fillUnderOverFlow(hresolutionbins_b1k_top_scatteringangle_ttbarframe,
            // residual_b1k_top_scatteringangle_ttbarframe_binning,
            // gen_b1k_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{b1k - gen_b1k}, std::vector<double>{gen_b1k,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_b2k_top_scatteringangle_ttbarframe,
            // residual_b2k_top_scatteringangle_ttbarframe_binning,
            // gen_b2k_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{b2k - gen_b2k}, std::vector<double>{gen_b2k,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_b1j_top_scatteringangle_ttbarframe,
            // residual_b1j_top_scatteringangle_ttbarframe_binning,
            // gen_b1j_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{b1j - gen_b1j}, std::vector<double>{gen_b1j,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_b2j_top_scatteringangle_ttbarframe,
            // residual_b2j_top_scatteringangle_ttbarframe_binning,
            // gen_b2j_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{b2j - gen_b2j}, std::vector<double>{gen_b2j,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_b1r_top_scatteringangle_ttbarframe,
            // residual_b1r_top_scatteringangle_ttbarframe_binning,
            // gen_b1r_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{b1r - gen_b1r}, std::vector<double>{gen_b1r,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_b2r_top_scatteringangle_ttbarframe,
            // residual_b2r_top_scatteringangle_ttbarframe_binning,
            // gen_b2r_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{b2r - gen_b2r}, std::vector<double>{gen_b2r,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_b1q_top_scatteringangle_ttbarframe,
            // residual_b1q_top_scatteringangle_ttbarframe_binning,
            // gen_b1q_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{b1q - gen_b1q}, std::vector<double>{gen_b1q,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_b2q_top_scatteringangle_ttbarframe,
            // residual_b2q_top_scatteringangle_ttbarframe_binning,
            // gen_b2q_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{b2q - gen_b2q}, std::vector<double>{gen_b2q,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_b1n_top_scatteringangle_ttbarframe,
            // residual_b1n_top_scatteringangle_ttbarframe_binning,
            // gen_b1n_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{b1n - gen_b1n}, std::vector<double>{gen_b1n,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_b2n_top_scatteringangle_ttbarframe,
            // residual_b2n_top_scatteringangle_ttbarframe_binning,
            // gen_b2n_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{b2n - gen_b2n}, std::vector<double>{gen_b2n,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);

            // // Diagonal elements
            // fillUnderOverFlow(hresolutionbins_c_kk_top_scatteringangle_ttbarframe,
            // residual_c_kk_top_scatteringangle_ttbarframe_binning,
            // gen_c_kk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{c_kk - gen_c_kk},
            // std::vector<double>{gen_c_kk,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_rr_top_scatteringangle_ttbarframe,
            // residual_c_rr_top_scatteringangle_ttbarframe_binning,
            // gen_c_rr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{c_rr - gen_c_rr},
            // std::vector<double>{gen_c_rr,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_nn_top_scatteringangle_ttbarframe,
            // residual_c_nn_top_scatteringangle_ttbarframe_binning,
            // gen_c_nn_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{c_nn - gen_c_nn},
            // std::vector<double>{gen_c_nn,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);

            // // Off diagonal elements
            // fillUnderOverFlow(hresolutionbins_c_rk_top_scatteringangle_ttbarframe,
            // residual_c_rk_top_scatteringangle_ttbarframe_binning,
            // gen_c_rk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{c_rk - gen_c_rk},
            // std::vector<double>{gen_c_rk,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_kr_top_scatteringangle_ttbarframe,
            // residual_c_kr_top_scatteringangle_ttbarframe_binning,
            // gen_c_kr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{c_kr - gen_c_kr},
            // std::vector<double>{gen_c_kr,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_nr_top_scatteringangle_ttbarframe,
            // residual_c_nr_top_scatteringangle_ttbarframe_binning,
            // gen_c_nr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{c_nr - gen_c_nr},
            // std::vector<double>{gen_c_nr,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_rn_top_scatteringangle_ttbarframe,
            // residual_c_rn_top_scatteringangle_ttbarframe_binning,
            // gen_c_rn_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{c_rn - gen_c_rn},
            // std::vector<double>{gen_c_rn,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_nk_top_scatteringangle_ttbarframe,
            // residual_c_nk_top_scatteringangle_ttbarframe_binning,
            // gen_c_nk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{c_nk - gen_c_nk},
            // std::vector<double>{gen_c_nk,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_kn_top_scatteringangle_ttbarframe,
            // residual_c_kn_top_scatteringangle_ttbarframe_binning,
            // gen_c_kn_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{c_kn - gen_c_kn},
            // std::vector<double>{gen_c_kn,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);

            // fillUnderOverFlow(hresolutionbins_c_Prk_top_scatteringangle_ttbarframe,
            // residual_c_Prk_top_scatteringangle_ttbarframe_binning,
            // gen_c_Prk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{(c_rk + c_kr) - (gen_c_rk + gen_c_kr)},
            // std::vector<double>{gen_c_rk + gen_c_kr,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_Mrk_top_scatteringangle_ttbarframe,
            // residual_c_Mrk_top_scatteringangle_ttbarframe_binning,
            // gen_c_Mrk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{(c_rk - c_kr) - (gen_c_rk - gen_c_kr)},
            // std::vector<double>{gen_c_rk - gen_c_kr,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_Pnr_top_scatteringangle_ttbarframe,
            // residual_c_Pnr_top_scatteringangle_ttbarframe_binning,
            // gen_c_Pnr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{(c_nr + c_rn) - (gen_c_nr + gen_c_rn)},
            // std::vector<double>{gen_c_nr + gen_c_rn,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_Mnr_top_scatteringangle_ttbarframe,
            // residual_c_Mnr_top_scatteringangle_ttbarframe_binning,
            // gen_c_Mnr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{(c_nr - c_rn) - (gen_c_nr - gen_c_rn)},
            // std::vector<double>{gen_c_nr - gen_c_rn,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_Pnk_top_scatteringangle_ttbarframe,
            // residual_c_Pnk_top_scatteringangle_ttbarframe_binning,
            // gen_c_Pnk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{(c_nk + c_kn) - (gen_c_nk + gen_c_kn)},
            // std::vector<double>{gen_c_nk + gen_c_kn,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_Mnk_top_scatteringangle_ttbarframe,
            // residual_c_Mnk_top_scatteringangle_ttbarframe_binning,
            // gen_c_Mnk_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{(c_nk - c_kn) - (gen_c_nk - gen_c_kn)},
            // std::vector<double>{gen_c_nk - gen_c_kn,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);

            // // Starred axes
            // fillUnderOverFlow(hresolutionbins_c_kj_top_scatteringangle_ttbarframe,
            // residual_c_kj_top_scatteringangle_ttbarframe_binning,
            // gen_c_kj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{c_kj - gen_c_kj},
            // std::vector<double>{gen_c_kj,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_rq_top_scatteringangle_ttbarframe,
            // residual_c_rq_top_scatteringangle_ttbarframe_binning,
            // gen_c_rq_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{c_rq - gen_c_rq},
            // std::vector<double>{gen_c_rq,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_rj_top_scatteringangle_ttbarframe,
            // residual_c_rj_top_scatteringangle_ttbarframe_binning,
            // gen_c_rj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{c_rj - gen_c_rj},
            // std::vector<double>{gen_c_rj,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_jr_top_scatteringangle_ttbarframe,
            // residual_c_jr_top_scatteringangle_ttbarframe_binning,
            // gen_c_jr_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{c_jr - gen_c_jr},
            // std::vector<double>{gen_c_jr,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);

            // fillUnderOverFlow(hresolutionbins_c_Prj_top_scatteringangle_ttbarframe,
            // residual_c_Prj_top_scatteringangle_ttbarframe_binning,
            // gen_c_Prj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{(c_rj + c_jr) - (gen_c_rj + gen_c_jr)},
            // std::vector<double>{gen_c_rj + gen_c_jr,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_Mrj_top_scatteringangle_ttbarframe,
            // residual_c_Mrj_top_scatteringangle_ttbarframe_binning,
            // gen_c_Mrj_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{(c_rj - c_jr) - (gen_c_rj - gen_c_jr)},
            // std::vector<double>{gen_c_rj - gen_c_jr,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);

            // // New variables
            // fillUnderOverFlow(hresolutionbins_c_han_top_scatteringangle_ttbarframe,
            // residual_c_han_top_scatteringangle_ttbarframe_binning,
            // gen_c_han_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{( gen_c_kk - gen_c_rr - gen_c_nn) - ( c_kk -
            // c_rr - c_nn)}, std::vector<double>{ gen_c_kk - gen_c_rr -
            // gen_c_nn, gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_sca_top_scatteringangle_ttbarframe,
            // residual_c_sca_top_scatteringangle_ttbarframe_binning,
            // gen_c_sca_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{(-gen_c_kk + gen_c_rr - gen_c_nn) - (-c_kk +
            // c_rr - c_nn)}, std::vector<double>{-gen_c_kk + gen_c_rr -
            // gen_c_nn, gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_tra_top_scatteringangle_ttbarframe,
            // residual_c_tra_top_scatteringangle_ttbarframe_binning,
            // gen_c_tra_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{(-gen_c_kk - gen_c_rr + gen_c_nn) - (-c_kk -
            // c_rr + c_nn)}, std::vector<double>{-gen_c_kk - gen_c_rr +
            // gen_c_nn, gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_kjL_top_scatteringangle_ttbarframe,
            // residual_c_kjL_top_scatteringangle_ttbarframe_binning,
            // gen_c_kjL_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{(-gen_c_kj - gen_c_rr - gen_c_nn) - (-c_kj -
            // c_rr - c_nn)}, std::vector<double>{-gen_c_kj - gen_c_rr -
            // gen_c_nn, gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_rqL_top_scatteringangle_ttbarframe,
            // residual_c_rqL_top_scatteringangle_ttbarframe_binning,
            // gen_c_rqL_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{(-gen_c_kk - gen_c_rq - gen_c_nn) - (-c_kk -
            // c_rq - c_nn)}, std::vector<double>{-gen_c_kk - gen_c_rq -
            // gen_c_nn, gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_rkP_top_scatteringangle_ttbarframe,
            // residual_c_rkP_top_scatteringangle_ttbarframe_binning,
            // gen_c_rkP_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{(-gen_c_rk - gen_c_kr - gen_c_nn) - (-c_rk -
            // c_kr - c_nn)}, std::vector<double>{-gen_c_rk - gen_c_kr -
            // gen_c_nn, gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_rkM_top_scatteringangle_ttbarframe,
            // residual_c_rkM_top_scatteringangle_ttbarframe_binning,
            // gen_c_rkM_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{(-gen_c_rk + gen_c_kr - gen_c_nn) - (-c_rk +
            // c_kr - c_nn)}, std::vector<double>{-gen_c_rk + gen_c_kr -
            // gen_c_nn, gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_nrP_top_scatteringangle_ttbarframe,
            // residual_c_nrP_top_scatteringangle_ttbarframe_binning,
            // gen_c_nrP_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{(-gen_c_nr - gen_c_rn - gen_c_kk) - (-c_nr -
            // c_rn - c_kk)}, std::vector<double>{-gen_c_nr - gen_c_rn -
            // gen_c_kk, gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_nrM_top_scatteringangle_ttbarframe,
            // residual_c_nrM_top_scatteringangle_ttbarframe_binning,
            // gen_c_nrM_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{(-gen_c_nr + gen_c_rn - gen_c_kk) - (-c_nr +
            // c_rn - c_kk)}, std::vector<double>{-gen_c_nr + gen_c_rn -
            // gen_c_kk, gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_nkP_top_scatteringangle_ttbarframe,
            // residual_c_nkP_top_scatteringangle_ttbarframe_binning,
            // gen_c_nkP_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{(-gen_c_nk - gen_c_kn - gen_c_rr) - (-c_nk -
            // c_kn - c_rr)}, std::vector<double>{-gen_c_nk - gen_c_kn -
            // gen_c_rr, gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_c_nkM_top_scatteringangle_ttbarframe,
            // residual_c_nkM_top_scatteringangle_ttbarframe_binning,
            // gen_c_nkM_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{(-gen_c_nk + gen_c_kn - gen_c_rr) - (-c_nk +
            // c_kn - c_rr)}, std::vector<double>{-gen_c_nk + gen_c_kn -
            // gen_c_rr, gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);

            // // Lab frame variables
            // fillUnderOverFlow(hresolutionbins_ll_cHel_top_scatteringangle_ttbarframe,
            // residual_ll_cHel_top_scatteringangle_ttbarframe_binning,
            // gen_ll_cHel_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{ll_cHel - gen_ll_cHel},
            // std::vector<double>{gen_ll_cHel,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_ll_cLab_top_scatteringangle_ttbarframe,
            // residual_ll_cLab_top_scatteringangle_ttbarframe_binning,
            // gen_ll_cLab_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{ll_cLab - gen_ll_cLab},
            // std::vector<double>{gen_ll_cLab,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_ll_kNorm_top_scatteringangle_ttbarframe,
            // residual_ll_kNorm_top_scatteringangle_ttbarframe_binning,
            // gen_ll_kNorm_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{ll_kNorm - gen_ll_kNorm},
            // std::vector<double>{gen_ll_kNorm,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_ll_rNorm_top_scatteringangle_ttbarframe,
            // residual_ll_rNorm_top_scatteringangle_ttbarframe_binning,
            // gen_ll_rNorm_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{ll_rNorm - gen_ll_rNorm},
            // std::vector<double>{gen_ll_rNorm,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_llbar_delta_phi_top_scatteringangle_ttbarframe,
            // residual_llbar_delta_phi_top_scatteringangle_ttbarframe_binning,
            // gen_llbar_delta_phi_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{llbar_delta_phi - gen_llbar_delta_phi},
            // std::vector<double>{gen_llbar_delta_phi,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
            // fillUnderOverFlow(hresolutionbins_llbar_delta_eta_top_scatteringangle_ttbarframe,
            // residual_llbar_delta_eta_top_scatteringangle_ttbarframe_binning,
            // gen_llbar_delta_eta_top_scatteringangle_ttbarframe_binning,
            // std::vector<double>{llbar_delta_eta - gen_llbar_delta_eta},
            // std::vector<double>{gen_llbar_delta_eta,
            // gen_top_scatteringangle_ttbarframe}, eventWeight, true,
            // symmetrizeSpinVar);
        }
    }
    wtdEvts->Write("hNrOfEvts");
}

void makeTUnfoldHisto::SetOutputFileName(std::string outputfile, std::string systematics, std::string channel,
                                         std::string era) {
    std::string dirname;
    // Amandeep : Adding support for <era>
    dirname = "UnfoldingHistos_" + era + "/" + systematics + "/" + channel;
    gSystem->mkdir(dirname.c_str(), true);
    outputfile = dirname + "/" + outputfile;
    fout = new TFile(outputfile.c_str(), "RECREATE");
    std::cout << "Output File: " << outputfile.c_str() << " successfully created" << std::endl;
}

/**
 * @brief Fills the 1D histogram correctly meaning nothing in the underflow or
 overflow
 *
 * @param h1 1D histogram to fill
 * @param binning The binning for the 1D histogram. This may be a
 multi-differential histogram that has been "unwrapped" to a 1D histogram and
 therefore requires a binning scheme
 * @param v_values The vector of values used to identify the correct bin to
 populate with the @param weight
 * @param weight The weight of the event
 * @param symmetrizeSpinVar Whether or not to symmetrize the histogram to blind
 ourselves to spin correlations
 */
void makeTUnfoldHisto::fillUnderOverFlow(TH1 *h1, TUnfoldBinning *binning, std::vector<double> v_values, double weight,
                                         bool symmetrizeSpinVar) {
    // copy to double array since TUnfoldBinning->GetGlobalBinNumber requires
    // double *
    double values[v_values.size()];
    for (int i = 0; i < v_values.size(); ++i) {
        values[i] = v_values[i];
    }

    // Symmetrize the spin correlation and polarization observables
    std::string histname;
    histname = h1->GetName();

    if (symmetrizeSpinVar == true &&
        ((histname.find("b_P") != std::string::npos) || (histname.find("b_M") != std::string::npos) ||
         (histname.find("b1") != std::string::npos) || (histname.find("b2") != std::string::npos) ||
         (histname.find("c_") != std::string::npos) || (histname.find("ll_c") != std::string::npos) ||
         (histname.find("Norm") != std::string::npos) || (histname.find("llbar_delta_phi") != std::string::npos) ||
         (histname.find("llbar_delta_eta") != std::string::npos))) {
        int binNumber = binning->GetGlobalBinNumber(values);

        int numValues = v_values.size();
        double symmetrizedValues[numValues];
        for (int i = 0; i < numValues; ++i) {
            // only symmetrize the spin correlation variable value --> first
            // value
            if (i == 0) {
                // i is the axis so c_kk, mttbar, etc.
                const TVectorD *vec_xBinEdges = binning->GetDistributionBinning(i);
                const double *xBinEdges = vec_xBinEdges->GetMatrixArray();
                symmetrizedValues[i] = xBinEdges[0] + xBinEdges[vec_xBinEdges->GetNoElements() - 1] - values[i];
            } else {  // keep other variables the same i.e. m_ttbar, top p_T etc
                symmetrizedValues[i] = values[i];
            }
        }
        int symmetricBinNumber = binning->GetGlobalBinNumber(symmetrizedValues);

        double binContent = h1->GetBinContent(binNumber);
        double binError = h1->GetBinError(binNumber);
        double symmetricBinContent = h1->GetBinContent(symmetricBinNumber);
        double symmetricBinError = h1->GetBinError(symmetricBinNumber);

        h1->SetBinContent(binNumber, binContent + 0.5 * weight);
        h1->SetBinError(binNumber, sqrt(pow(binError, 2.0) + pow(0.5 * weight, 2.0)));
        h1->SetBinContent(symmetricBinNumber, symmetricBinContent + 0.5 * weight);
        h1->SetBinError(symmetricBinNumber, sqrt(pow(symmetricBinError, 2.0) + pow(0.5 * weight, 2.0)));
    } else {
        int binNumber = binning->GetGlobalBinNumber(values);

        double binContent = h1->GetBinContent(binNumber);
        double binError = h1->GetBinError(binNumber);
        h1->SetBinContent(binNumber, binContent + weight);
        h1->SetBinError(binNumber, sqrt(pow(binError, 2.0) + pow(weight, 2.0)));
    }
}

/**
 * @brief Fills the 2D histogram correctly meaning non-reconstructed events go
 in the underflow of the y-axis
 *
 * @param h2 2D histogram that needs to be filled
 * @param xBinning The binning scheme of the 2D histogram along the x-axis. This
 may be a multi-differential histogram that has been "unwrapped" to a 2D
 histogram and therefore requires a binning scheme
 * @param yBinning The binning scheme of the 2D histogram along the y-axis. This
 may be a multi-differential histogram that has been "unwrapped" to a 2D
 histogram and therefore requires a binning scheme
 * @param v_xvalues The vector of values that are used to determine which bin
 along the x-axis needs to be filled
 * @param v_yvalues The vector of values that are used to determine which bin
 along the y-axis needs to be filled
 * @param weight The weight of the event
 * @param isReconstructed Whether this event has been reconstructed or not
 * @param symmetrizeSpinVar Whether or not to symmetrize the histogram to blind
 ourselves to spin correlations
 */
void makeTUnfoldHisto::fillUnderOverFlow(TH2 *h2, TUnfoldBinning *xBinning, TUnfoldBinning *yBinning,
                                         std::vector<double> v_xvalues, std::vector<double> v_yvalues, double weight,
                                         bool isReconstructed, bool symmetrizeSpinVar) {
    std::string histname = h2->GetName();
    // if((histname.find("b1k") != std::string::npos)){

    //     std::cout << histname << " begin2Dfill, " << h2->Integral() << "," <<
    //     v_xvalues[0] << "," <<v_xvalues[1] << "," <<weight << std::endl;
    //     std::cout << histname << " begin2Dfill, " <<
    //     h2->ProjectionY()->GetBinContent(0) << "," << v_yvalues[0] << ","
    //     <<v_yvalues[1] << "," <<h2->ProjectionY()->GetBinContent(0) << ","
    //     <<h2->GetBinContent(0,0) <<std::endl; std::cout << "," << std::endl;

    // }

    // copy to double arrays since TUnfoldBinning->GetGlobalBinNumber requires
    // double *
    double xvalues[v_xvalues.size()];
    for (int i = 0; i < v_xvalues.size(); ++i) {
        xvalues[i] = v_xvalues[i];
    }

    double yvalues[v_yvalues.size()];
    for (int i = 0; i < v_yvalues.size(); ++i) {
        yvalues[i] = v_yvalues[i];
    }

    int globalXBinNumber = xBinning->GetGlobalBinNumber(xvalues);
    int numXValues = v_xvalues.size();
    double symmetrizedXValues[numXValues];
    for (int i = 0; i < numXValues; ++i) {
        // only symmetrize the spin correlation variable value --> first value
        if (i == 0) {
            const TVectorD *vec_xBinEdges = xBinning->GetDistributionBinning(i);
            const double *xBinEdges = vec_xBinEdges->GetMatrixArray();
            symmetrizedXValues[i] = xBinEdges[0] + xBinEdges[vec_xBinEdges->GetNoElements() - 1] - xvalues[i];
        } else {
            symmetrizedXValues[i] = xvalues[i];
        }
    }
    int globalSymmetricXBinNumber = xBinning->GetGlobalBinNumber(symmetrizedXValues);

    int globalYBinNumber = yBinning->GetGlobalBinNumber(yvalues);
    int numYValues = v_yvalues.size();
    double symmetrizedYValues[numYValues];
    for (int i = 0; i < numYValues; ++i) {
        // only symmetrize the spin correlation variable value --> first value
        if (i == 0) {
            const TVectorD *vec_yBinEdges = yBinning->GetDistributionBinning(i);
            const double *yBinEdges = vec_yBinEdges->GetMatrixArray();
            symmetrizedYValues[i] = yBinEdges[0] + yBinEdges[vec_yBinEdges->GetNoElements() - 1] - yvalues[i];
        } else {
            symmetrizedYValues[i] = yvalues[i];
        }
    }
    int globalSymmetricYBinNumber = yBinning->GetGlobalBinNumber(symmetrizedYValues);

    double binContent = h2->GetBinContent(globalXBinNumber, globalYBinNumber);
    double binError = h2->GetBinError(globalXBinNumber, globalYBinNumber);
    double symmetricBinContent = h2->GetBinContent(globalSymmetricXBinNumber, globalSymmetricYBinNumber);
    double symmetricBinError = h2->GetBinError(globalSymmetricXBinNumber, globalSymmetricYBinNumber);

    // std::string histname = h2->GetName();

    // Symmetrize the spin correlation and polarization observables

    if (symmetrizeSpinVar == true &&
        ((histname.find("b_P") != std::string::npos) || (histname.find("b_M") != std::string::npos) ||
         (histname.find("b1") != std::string::npos) || (histname.find("b2") != std::string::npos) ||
         (histname.find("c_") != std::string::npos) || (histname.find("ll_c") != std::string::npos) ||
         (histname.find("Norm") != std::string::npos) || (histname.find("llbar_delta_phi") != std::string::npos))) {
        // If reconstructed fill normally

        if (isReconstructed) {
            h2->SetBinContent(globalXBinNumber, globalYBinNumber, binContent + 0.5 * weight);
            h2->SetBinError(globalXBinNumber, globalYBinNumber, sqrt(pow(binError, 2.0) + pow(0.5 * weight, 2.0)));
            h2->SetBinContent(globalSymmetricXBinNumber, globalSymmetricYBinNumber, symmetricBinContent + 0.5 * weight);
            h2->SetBinError(globalSymmetricXBinNumber, globalSymmetricYBinNumber,
                            sqrt(pow(symmetricBinError, 2.0) + pow(0.5 * weight, 2.0)));
        }

        else {
            // If not reconstructed then put in underflow of y-axis
            // Do not want the events to go to overflow so no
            // -1.0*(h2->GetYaxis()->GetBinLowEdge(0))
            globalYBinNumber = 0;
            globalSymmetricYBinNumber = 0;

            binContent = h2->GetBinContent(globalXBinNumber, globalYBinNumber);
            binError = h2->GetBinError(globalXBinNumber, globalYBinNumber);
            symmetricBinContent = h2->GetBinContent(globalSymmetricXBinNumber, globalSymmetricYBinNumber);
            symmetricBinError = h2->GetBinError(globalSymmetricXBinNumber, globalSymmetricYBinNumber);

            h2->SetBinContent(globalXBinNumber, globalYBinNumber, binContent + 0.5 * weight);
            h2->SetBinError(globalXBinNumber, globalYBinNumber, sqrt(pow(binError, 2.0) + pow(0.5 * weight, 2.0)));
            h2->SetBinContent(globalSymmetricXBinNumber, globalSymmetricYBinNumber, symmetricBinContent + 0.5 * weight);
            h2->SetBinError(globalSymmetricXBinNumber, globalSymmetricYBinNumber,
                            sqrt(pow(symmetricBinError, 2.0) + pow(0.5 * weight, 2.0)));
        }
    }

    else {
        // Do not symmetrize the spin correlation and polarization observables
        // If reconstructed fill normally
        if (isReconstructed) {
            h2->SetBinContent(globalXBinNumber, globalYBinNumber, binContent + weight);
            h2->SetBinError(globalXBinNumber, globalYBinNumber, sqrt(pow(binError, 2.0) + pow(weight, 2.0)));
        } else {
            // Fill non-reconstructed events in the Y axis underflow bin
            // Do not want the events to go to overflow so do not do:
            // -1.0*(h2->GetYaxis()->GetBinLowEdge(0))
            globalYBinNumber = 0;

            binContent = h2->GetBinContent(globalXBinNumber, globalYBinNumber);
            binError = h2->GetBinError(globalXBinNumber, globalYBinNumber);

            h2->SetBinContent(globalXBinNumber, globalYBinNumber, binContent + weight);
            h2->SetBinError(globalXBinNumber, globalYBinNumber, sqrt(pow(binError, 2.0) + pow(weight, 2.0)));
        }
    }
}

/**
 * @brief Reconstructed events have an additional weight weightReco =
 reconstructionScaleFactor * weightGen where reconstructionScaleFactor >= 1.0.
 This is to account for acceptance x efficiency effects in the reconstructed
 distribution. However, because of this when we unfold to the final GEN level
 distribution we need to correct for the fact that we corrected our
 reconstruction level distributions. Essentially we do not want to double count
 the acceptance x efficiency effects since the underflow bin was already filled
 with events that were generated but not reconstructed. To do this, we need to
 fill our reconstruction-level underflow bin with a deltaWeight = weightGen -
 weightReco >= 0. Since weights are added in quadrature, for the overall
 variance we need to manually subtract the weightReco part of the variance from
 the filled bin's variance as well.
 *
 * @param h2 The 2D histogram where the deltaWeight correction needs to be
 applied.
 * @param xBinning The binning scheme along the x-aixs (GEN) to properly handle
 multi-differential distributions that have been "flattened" to a 1D axis.
 * @param v_xvalues The vector of x values (GEN) that we are using to correct
 the weights in the y-axis underflow (RECO)
 * @param weightGen The generator-level weight for this event.
 * @param weightReco The corrected reconstruction-level weight for this event.
 * @param symmetrizeSpinVar Whether or not this distribution is currently being
 symmetrized to blind ourselves from spincorrelation.
 */

void makeTUnfoldHisto::fillUnderOverFlowDeltaWeight(TH2 *h2, TUnfoldBinning *xBinning, std::vector<double> v_xvalues,
                                                    double weightGen, double weightReco, bool symmetrizeSpinVar) {
    std::string histname;
    histname = h2->GetName();

    // if((histname.find("b1k") != std::string::npos)){

    //     std::cout << histname << " begin2DdeltaWeight, " << h2->Integral() <<
    //     "," << v_xvalues[0] << "," <<v_xvalues[1] <<  std::endl; std::cout <<
    //     histname << " begin2DdeltaWeight, " <<
    //     h2->ProjectionY()->GetBinContent(0) << "," << weightGen << ","
    //     <<weightReco << "," <<h2->ProjectionY()->GetBinContent(0) << ","
    //     <<h2->GetBinContent(0,0) <<std::endl; std::cout << "," << std::endl;

    // }
    // copy to double arrays since TUnfoldBinning->GetGlobalBinNumber requires
    // double *
    double xvalues[v_xvalues.size()];
    for (int i = 0; i < v_xvalues.size(); ++i) {
        xvalues[i] = v_xvalues[i];
    }

    int globalXBinNumber = xBinning->GetGlobalBinNumber(xvalues);
    int numXValues = v_xvalues.size();
    double symmetrizedXValues[numXValues];
    for (int i = 0; i < numXValues; ++i) {
        // only symmetrize the spin correlation variable value --> first value
        if (i == 0) {
            const TVectorD *vec_xBinEdges = xBinning->GetDistributionBinning(i);
            const double *xBinEdges = vec_xBinEdges->GetMatrixArray();
            symmetrizedXValues[i] = xBinEdges[0] + xBinEdges[vec_xBinEdges->GetNoElements() - 1] - xvalues[i];
        } else {
            symmetrizedXValues[i] = xvalues[i];
        }
    }
    int globalSymmetricXBinNumber = xBinning->GetGlobalBinNumber(symmetrizedXValues);

    // Do not want the events to go to overflow so no
    // -1.0*(h2->GetYaxis()->GetBinLowEdge(0))
    double globalYBinNumber = 0;
    double globalSymmetricYBinNumber = 0;

    double binContent = h2->GetBinContent(globalXBinNumber, globalYBinNumber);
    double binError = h2->GetBinError(globalXBinNumber, globalYBinNumber);
    double symmetricBinContent = h2->GetBinContent(globalSymmetricXBinNumber, globalSymmetricYBinNumber);
    double symmetricBinError = h2->GetBinError(globalSymmetricXBinNumber, globalSymmetricYBinNumber);

    // Symmetrize the spin correlation and polarization observables
    std::vector<int> xfactor;
    // std::string histname;
    // histname = h2->GetName();
    if (symmetrizeSpinVar == true &&
        ((histname.find("b_P") != std::string::npos) || (histname.find("b_M") != std::string::npos) ||
         (histname.find("b1") != std::string::npos) || (histname.find("b2") != std::string::npos) ||
         (histname.find("c_") != std::string::npos) || (histname.find("ll_c") != std::string::npos) ||
         (histname.find("Norm") != std::string::npos) || (histname.find("llbar_delta_phi") != std::string::npos))) {
        double new_bin_error = sqrt(pow(binError, 2.0) + pow(0.5 * weightGen, 2.0) - pow(0.5 * weightReco, 2.0));
        double new_symmetric_bin_error =
            sqrt(pow(symmetricBinError, 2.0) + pow(0.5 * weightGen, 2.0) - pow(0.5 * weightReco, 2.0));
        if (new_bin_error < 0) {
            std::cout << "WARNING in makeTUnfoldHisto::fillUnderOverFlowDeltaWeight: "
                         "negative new_bin_error^2, setting error to zero."
                      << std::endl;
            new_bin_error = 0;
        } else if (new_symmetric_bin_error < 0) {
            std::cout << "WARNING in makeTUnfoldHisto::fillUnderOverFlowDeltaWeight: "
                         "negative new_symmetric_bin_error^2, setting error to zero."
                      << std::endl;
            new_symmetric_bin_error = 0;
        }

        h2->SetBinContent(globalXBinNumber, globalYBinNumber, binContent + 0.5 * (weightGen - weightReco));
        h2->SetBinError(globalXBinNumber, globalYBinNumber, new_bin_error);
        h2->SetBinContent(globalSymmetricXBinNumber, globalSymmetricYBinNumber,
                          symmetricBinContent + 0.5 * (weightGen - weightReco));
        h2->SetBinError(globalSymmetricXBinNumber, globalSymmetricYBinNumber, new_symmetric_bin_error);
    } else {
        double new_bin_error = sqrt(pow(binError, 2.0) + pow(weightGen, 2.0) - pow(weightReco, 2.0));
        if (new_bin_error < 0) {
            std::cout << "WARNING in makeTUnfoldHisto::fillUnderOverFlowDeltaWeight: "
                         "negative new_bin_error^2, setting error to zero."
                      << std::endl;
            new_bin_error = 0;
        }
        h2->SetBinContent(globalXBinNumber, globalYBinNumber, binContent + weightGen - weightReco);
        h2->SetBinError(globalXBinNumber, globalYBinNumber, new_bin_error);
    }
}

// Booking histograms
void makeTUnfoldHisto::BookHisto(int n_sub_bins) {
    // Amandeep : Add int n_sub_bins_2D; rebin for 2D is 2, rebin for 1D is 4
    // int n_sub_bins_2D   = int(n_sub_bins * 0.5); // rebin of 2
    int n_sub_bins_2D = int(n_sub_bins);  // rebin of 4

    if (!fout) {
        std::cout << "Call setoutput filename first" << std::endl;
        exit(EXIT_FAILURE);
    }

    fout->cd();

    TH1::SetDefaultSumw2();

    // ***************
    // CUSTOM BINNINGS
    // ***************

    // Amandeep : Modifying DESY top pt binning

    // 1D top_pt binning with rebin 4
    const int gen_ntop_ptbin_final = 4;
    const int reco_ntop_ptbin_final = 8;
    double gen_top_ptbins_final[gen_ntop_ptbin_final + 1] = {0.0, 80.0, 150.0, 250.0, 550.0};
    double reco_top_ptbins_final[reco_ntop_ptbin_final + 1] = {0.0,   40.0,  80.0,  115.0, 150.0,
                                                               200.0, 250.0, 400.0, 550.0};
    const int gen_ntop_ptbin = gen_ntop_ptbin_final * n_sub_bins;
    const int reco_ntop_ptbin = reco_ntop_ptbin_final * n_sub_bins;
    double gen_top_ptbins[gen_ntop_ptbin + 1];
    double reco_top_ptbins[reco_ntop_ptbin + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_ntop_ptbin_final, gen_ntop_ptbin, gen_top_ptbins_final,
                                             gen_top_ptbins);
    makeTUnfoldHisto::CreateFineBinningArray(reco_ntop_ptbin_final, reco_ntop_ptbin, reco_top_ptbins_final,
                                             reco_top_ptbins);

    // 2D top_pt binning with rebin 2
    const int gen_ntop_ptbin_final_2D = 4;
    const int reco_ntop_ptbin_final_2D = 8;
    double gen_top_ptbins_final_2D[gen_ntop_ptbin_final_2D + 1] = {0.0, 80.0, 150.0, 250.0, 550.0};
    double reco_top_ptbins_final_2D[reco_ntop_ptbin_final_2D + 1] = {0.0,   40.0,  80.0,  115.0, 150.0,
                                                                     200.0, 250.0, 400.0, 550.0};
    const int gen_ntop_ptbin_2D = gen_ntop_ptbin_final_2D * n_sub_bins_2D;
    const int reco_ntop_ptbin_2D = reco_ntop_ptbin_final_2D * n_sub_bins_2D;
    double gen_top_ptbins_2D[gen_ntop_ptbin_2D + 1];
    double reco_top_ptbins_2D[reco_ntop_ptbin_2D + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_ntop_ptbin_final_2D, gen_ntop_ptbin_2D, gen_top_ptbins_final_2D,
                                             gen_top_ptbins_2D);
    makeTUnfoldHisto::CreateFineBinningArray(reco_ntop_ptbin_final_2D, reco_ntop_ptbin_2D, reco_top_ptbins_final_2D,
                                             reco_top_ptbins_2D);

    // Amandeep : Modifying DESY mttbar binning
    // 1D mttbar binning with rebin 4
    const int gen_nttbar_massbin_final = 4;
    const int reco_nttbar_massbin_final = 8;
    double gen_ttbar_massbins_final[gen_nttbar_massbin_final + 1] = {250.0, 450.0, 600.0, 800.0, 2000.0};
    double reco_ttbar_massbins_final[reco_nttbar_massbin_final + 1] = {250.0, 350.0, 450.0,  525.0, 600.0,
                                                                       700.0, 800.0, 1400.0, 2000.0};
    const int gen_nttbar_massbin = gen_nttbar_massbin_final * n_sub_bins;
    const int reco_nttbar_massbin = reco_nttbar_massbin_final * n_sub_bins;
    double gen_ttbar_massbins[gen_nttbar_massbin + 1];
    double reco_ttbar_massbins[reco_nttbar_massbin + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_nttbar_massbin_final, gen_nttbar_massbin, gen_ttbar_massbins_final,
                                             gen_ttbar_massbins);
    makeTUnfoldHisto::CreateFineBinningArray(reco_nttbar_massbin_final, reco_nttbar_massbin, reco_ttbar_massbins_final,
                                             reco_ttbar_massbins);

    // 2D mttbar with a rebin factor of 2 instead of 4
    const int gen_nttbar_massbin_final_2D = 4;
    const int reco_nttbar_massbin_final_2D = 8;
    double gen_ttbar_massbins_final_2D[gen_nttbar_massbin_final_2D + 1] = {250.0, 450.0, 600.0, 800.0, 2000.0};
    double reco_ttbar_massbins_final_2D[reco_nttbar_massbin_final_2D + 1] = {250.0, 350.0, 450.0,  525.0, 600.0,
                                                                             700.0, 800.0, 1400.0, 2000.0};
    const int gen_nttbar_massbin_2D = gen_nttbar_massbin_final_2D * n_sub_bins_2D;
    const int reco_nttbar_massbin_2D = reco_nttbar_massbin_final_2D * n_sub_bins_2D;
    double gen_ttbar_massbins_2D[gen_nttbar_massbin_2D + 1];
    double reco_ttbar_massbins_2D[reco_nttbar_massbin_2D + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_nttbar_massbin_final_2D, gen_nttbar_massbin_2D,
                                             gen_ttbar_massbins_final_2D, gen_ttbar_massbins_2D);
    makeTUnfoldHisto::CreateFineBinningArray(reco_nttbar_massbin_final_2D, reco_nttbar_massbin_2D,
                                             reco_ttbar_massbins_final_2D, reco_ttbar_massbins_2D);

    // Abraham: Extra-jet bining

    const int gen_nexjetbin_final = 4;
    const int reco_nexjetbin_final = 8;
    double gen_exjetbins_final[gen_nexjetbin_final + 1] = {-0.5, 0.5, 1.5, 2.5, 7.5};
    double reco_exjetbins_final[reco_nexjetbin_final + 1] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
    const int gen_nexjetbin = gen_nexjetbin_final * n_sub_bins;
    const int reco_nexjetbin = reco_nexjetbin_final * n_sub_bins;
    double gen_exjetbins[gen_nexjetbin + 1];
    double reco_exjetbins[reco_nexjetbin + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_nexjetbin_final, gen_nexjetbin, gen_exjetbins_final, gen_exjetbins);
    makeTUnfoldHisto::CreateFineBinningArray(reco_nexjetbin_final, reco_nexjetbin, reco_exjetbins_final,
                                             reco_exjetbins);

    // 2D mttbar with a rebin factor of 2 instead of 4
    const int gen_nexjetbin_final_2D = 4;
    const int reco_nexjetbin_final_2D = 4;
    double gen_exjetbins_final_2D[gen_nexjetbin_final_2D + 1] = {-0.5, 0.5, 1.5, 2.5, 3.5};
    double reco_exjetbins_final_2D[reco_nexjetbin_final_2D + 1] = {-0.5, 0.5, 1.5, 2.5, 3.5};
    const int gen_nexjetbin_2D = gen_nexjetbin_final_2D;
    const int reco_nexjetbin_2D = reco_nexjetbin_final_2D;
    double gen_exjetbins_2D[gen_nexjetbin_2D + 1] = {-0.5, 0.5, 1.5, 2.5, 3.5};
    double reco_exjetbins_2D[reco_nexjetbin_2D + 1] = {-0.5, 0.5, 1.5, 2.5, 3.5};
    // makeTUnfoldHisto::CreateFineBinningArray(gen_nexjetbin_final_2D ,
    // gen_nexjetbin_2D,  gen_exjetbins_final_2D, gen_exjetbins_2D);
    // makeTUnfoldHisto::CreateFineBinningArray(reco_nexjetbin_final_2D,
    // reco_nexjetbin_2D, reco_exjetbins_final_2D, reco_exjetbins_2D);

    // cos_theta
    // const int gen_ncsthetabin_final  = 4;
    // const int reco_ncsthetabin_final = 8;
    // double gen_csthetabins_final[gen_ncsthetabin_final+1]   = {-1.0, -0.5,
    // 0.0, 0.5, 1}; double reco_csthetabins_final[reco_ncsthetabin_final+1] =
    // {-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.50, 0.75, 1.0}; const int
    // gen_ncsthetabin = gen_ncsthetabin_final * n_sub_bins; const int
    // reco_ncsthetabin = reco_ncsthetabin_final * n_sub_bins; double
    // gen_csthetabins[gen_ncsthetabin+1]; double
    // reco_csthetabins[reco_ncsthetabin+1];
    // makeTUnfoldHisto::CreateFineBinningArray(gen_ncsthetabin_final,
    // gen_ncsthetabin, gen_csthetabins_final, gen_csthetabins);
    // makeTUnfoldHisto::CreateFineBinningArray(reco_ncsthetabin_final,
    // reco_ncsthetabin, reco_csthetabins_final, reco_csthetabins);

    // 1D delta_eta binning
    const int gen_nllbar_delta_etabin_final = 6;
    const int reco_nllbar_delta_etabin_final = 12;
    double gen_llbar_delta_etabins_final[gen_nllbar_delta_etabin_final + 1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 5.0};
    double reco_llbar_delta_etabins_final[reco_nllbar_delta_etabin_final + 1] = {
        0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.75, 5.0};
    const int gen_nllbar_delta_etabin = gen_nllbar_delta_etabin_final * n_sub_bins;
    const int reco_nllbar_delta_etabin = reco_nllbar_delta_etabin_final * n_sub_bins;
    double gen_llbar_delta_etabins[gen_nllbar_delta_etabin + 1];
    double reco_llbar_delta_etabins[reco_nllbar_delta_etabin + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_nllbar_delta_etabin_final, gen_nllbar_delta_etabin,
                                             gen_llbar_delta_etabins_final, gen_llbar_delta_etabins);
    makeTUnfoldHisto::CreateFineBinningArray(reco_nllbar_delta_etabin_final, reco_nllbar_delta_etabin,
                                             reco_llbar_delta_etabins_final, reco_llbar_delta_etabins);

    // 2D delta_eta binning
    const int gen_nllbar_delta_etabin_final_2D = 6;
    const int reco_nllbar_delta_etabin_final_2D = 12;
    double gen_llbar_delta_etabins_final_2D[gen_nllbar_delta_etabin_final_2D + 1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 5.0};
    double reco_llbar_delta_etabins_final_2D[reco_nllbar_delta_etabin_final_2D + 1] = {
        0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.75, 5.0};
    const int gen_nllbar_delta_etabin_2D = gen_nllbar_delta_etabin_final_2D * n_sub_bins_2D;
    const int reco_nllbar_delta_etabin_2D = reco_nllbar_delta_etabin_final_2D * n_sub_bins_2D;
    double gen_llbar_delta_etabins_2D[gen_nllbar_delta_etabin_2D + 1];
    double reco_llbar_delta_etabins_2D[reco_nllbar_delta_etabin_2D + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_nllbar_delta_etabin_final_2D, gen_nllbar_delta_etabin_2D,
                                             gen_llbar_delta_etabins_final_2D, gen_llbar_delta_etabins_2D);
    makeTUnfoldHisto::CreateFineBinningArray(reco_nllbar_delta_etabin_final_2D, reco_nllbar_delta_etabin_2D,
                                             reco_llbar_delta_etabins_final_2D, reco_llbar_delta_etabins_2D);

    // 1D delta_phi binning
    const int gen_nllbar_delta_phibin_final = 6;
    const int reco_nllbar_delta_phibin_final = 12;
    double gen_llbar_delta_phibins_final[gen_nllbar_delta_phibin_final + 1] = {0.,
                                                                               1. * TMath::Pi() / 6.,
                                                                               2. * TMath::Pi() / 6.,
                                                                               3. * TMath::Pi() / 6.,
                                                                               4. * TMath::Pi() / 6.,
                                                                               5. * TMath::Pi() / 6.,
                                                                               TMath::Pi()};
    double reco_llbar_delta_phibins_final[reco_nllbar_delta_phibin_final + 1] = {0.,
                                                                                 1. * TMath::Pi() / 12.,
                                                                                 2. * TMath::Pi() / 12.,
                                                                                 3. * TMath::Pi() / 12.,
                                                                                 4. * TMath::Pi() / 12.,
                                                                                 5. * TMath::Pi() / 12.,
                                                                                 6. * TMath::Pi() / 12.,
                                                                                 7. * TMath::Pi() / 12.,
                                                                                 8. * TMath::Pi() / 12.,
                                                                                 9. * TMath::Pi() / 12.,
                                                                                 10. * TMath::Pi() / 12.,
                                                                                 11. * TMath::Pi() / 12.,
                                                                                 TMath::Pi()};
    const int gen_nllbar_delta_phibin = gen_nllbar_delta_phibin_final * n_sub_bins;
    const int reco_nllbar_delta_phibin = reco_nllbar_delta_phibin_final * n_sub_bins;
    double gen_llbar_delta_phibins[gen_nllbar_delta_phibin + 1];
    double reco_llbar_delta_phibins[reco_nllbar_delta_phibin + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_nllbar_delta_phibin_final, gen_nllbar_delta_phibin,
                                             gen_llbar_delta_phibins_final, gen_llbar_delta_phibins);
    makeTUnfoldHisto::CreateFineBinningArray(reco_nllbar_delta_phibin_final, reco_nllbar_delta_phibin,
                                             reco_llbar_delta_phibins_final, reco_llbar_delta_phibins);

    // 2D delta_phi binning
    const int gen_nllbar_delta_phibin_final_2D = 6;
    const int reco_nllbar_delta_phibin_final_2D = 12;
    double gen_llbar_delta_phibins_final_2D[gen_nllbar_delta_phibin_final_2D + 1] = {0.,
                                                                                     1. * TMath::Pi() / 6.,
                                                                                     2. * TMath::Pi() / 6.,
                                                                                     3. * TMath::Pi() / 6.,
                                                                                     4. * TMath::Pi() / 6.,
                                                                                     5. * TMath::Pi() / 6.,
                                                                                     TMath::Pi()};
    double reco_llbar_delta_phibins_final_2D[reco_nllbar_delta_phibin_final_2D + 1] = {0.,
                                                                                       1. * TMath::Pi() / 12.,
                                                                                       2. * TMath::Pi() / 12.,
                                                                                       3. * TMath::Pi() / 12.,
                                                                                       4. * TMath::Pi() / 12.,
                                                                                       5. * TMath::Pi() / 12.,
                                                                                       6. * TMath::Pi() / 12.,
                                                                                       7. * TMath::Pi() / 12.,
                                                                                       8. * TMath::Pi() / 12.,
                                                                                       9. * TMath::Pi() / 12.,
                                                                                       10. * TMath::Pi() / 12.,
                                                                                       11. * TMath::Pi() / 12.,
                                                                                       TMath::Pi()};
    const int gen_nllbar_delta_phibin_2D = gen_nllbar_delta_phibin_final_2D * n_sub_bins_2D;
    const int reco_nllbar_delta_phibin_2D = reco_nllbar_delta_phibin_final_2D * n_sub_bins_2D;
    double gen_llbar_delta_phibins_2D[gen_nllbar_delta_phibin_2D + 1];
    double reco_llbar_delta_phibins_2D[reco_nllbar_delta_phibin_2D + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_nllbar_delta_phibin_final_2D, gen_nllbar_delta_phibin_2D,
                                             gen_llbar_delta_phibins_final_2D, gen_llbar_delta_phibins_2D);
    makeTUnfoldHisto::CreateFineBinningArray(reco_nllbar_delta_phibin_final_2D, reco_nllbar_delta_phibin_2D,
                                             reco_llbar_delta_phibins_final_2D, reco_llbar_delta_phibins_2D);

    // This is used for all almost all spin corr vars
    const int gen_ncbin_final = 6;
    const int reco_ncbin_final = 12;
    double gen_cbins_final[gen_ncbin_final + 1] = {-1.0, -2. / 3., -1. / 3., 0., 1. / 3., 2. / 3., 1.0};
    double reco_cbins_final[reco_ncbin_final + 1] = {-1.0,    -5. / 6., -2. / 3., -1. / 2., -1. / 3., -1. / 6., 0.,
                                                     1. / 6., 1. / 3.,  1. / 2.,  2. / 3.,  5. / 6.,  1.0};
    const int gen_ncbin = gen_ncbin_final * n_sub_bins;
    const int reco_ncbin = reco_ncbin_final * n_sub_bins;
    double gen_cbins[gen_ncbin + 1];
    double reco_cbins[reco_ncbin + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_ncbin_final, gen_ncbin, gen_cbins_final, gen_cbins);
    makeTUnfoldHisto::CreateFineBinningArray(reco_ncbin_final, reco_ncbin, reco_cbins_final, reco_cbins);

    // 2D Spin corr vars where we need rebin of 2 instead of 4
    const int gen_ncbin_final_2D = 6;
    const int reco_ncbin_final_2D = 12;
    double gen_cbins_final_2D[gen_ncbin_final_2D + 1] = {-1.0, -2. / 3., -1. / 3., 0., 1. / 3., 2. / 3., 1.0};
    double reco_cbins_final_2D[reco_ncbin_final_2D + 1] = {
        -1.0, -5. / 6., -2. / 3., -1. / 2., -1. / 3., -1. / 6., 0., 1. / 6., 1. / 3., 1. / 2., 2. / 3., 5. / 6., 1.0};
    const int gen_ncbin_2D = gen_ncbin_final_2D * n_sub_bins_2D;
    const int reco_ncbin_2D = reco_ncbin_final_2D * n_sub_bins_2D;
    double gen_cbins_2D[gen_ncbin_2D + 1];
    double reco_cbins_2D[reco_ncbin_2D + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_ncbin_final_2D, gen_ncbin_2D, gen_cbins_final_2D, gen_cbins_2D);
    makeTUnfoldHisto::CreateFineBinningArray(reco_ncbin_final_2D, reco_ncbin_2D, reco_cbins_final_2D, reco_cbins_2D);

    // Used for 1D ll cosine variables
    const int gen_ncosbin_final = 6;
    const int reco_ncosbin_final = 12;
    double gen_cosbins_final[gen_ncosbin_final + 1] = {-1.0, -2. / 3., -1. / 3., 0., 1. / 3., 2. / 3., 1.0};
    double reco_cosbins_final[reco_ncosbin_final + 1] = {-1.0,    -5. / 6., -2. / 3., -1. / 2., -1. / 3., -1. / 6., 0.,
                                                         1. / 6., 1. / 3.,  1. / 2.,  2. / 3.,  5. / 6.,  1.0};
    const int gen_ncosbin = gen_ncosbin_final * n_sub_bins;
    const int reco_ncosbin = reco_ncosbin_final * n_sub_bins;
    double gen_cosbins[gen_ncosbin + 1];
    double reco_cosbins[reco_ncosbin + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_ncosbin_final, gen_ncosbin, gen_cosbins_final, gen_cosbins);
    makeTUnfoldHisto::CreateFineBinningArray(reco_ncosbin_final, reco_ncosbin, reco_cosbins_final, reco_cosbins);

    // Use for 2D ll cosine variables
    const int gen_ncosbin_final_2D = 6;
    const int reco_ncosbin_final_2D = 12;
    double gen_cosbins_final_2D[gen_ncosbin_final_2D + 1] = {-1.0, -2. / 3., -1. / 3., 0., 1. / 3., 2. / 3., 1.0};
    double reco_cosbins_final_2D[reco_ncosbin_final_2D + 1] = {
        -1.0, -5. / 6., -2. / 3., -1. / 2., -1. / 3., -1. / 6., 0., 1. / 6., 1. / 3., 1. / 2., 2. / 3., 5. / 6., 1.0};
    const int gen_ncosbin_2D = gen_ncosbin_final_2D * n_sub_bins_2D;
    const int reco_ncosbin_2D = reco_ncosbin_final_2D * n_sub_bins_2D;
    double gen_cosbins_2D[gen_ncosbin_2D + 1];
    double reco_cosbins_2D[reco_ncosbin_2D + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_ncosbin_final_2D, gen_ncosbin_2D, gen_cosbins_final_2D,
                                             gen_cosbins_2D);
    makeTUnfoldHisto::CreateFineBinningArray(reco_ncosbin_final_2D, reco_ncosbin_2D, reco_cosbins_final_2D,
                                             reco_cosbins_2D);

    // DESY Kinematic variables
    const int gen_nl_ptbin_final = 6;
    const int reco_nl_ptbin_final = 12;
    double gen_l_ptbins_final[gen_nl_ptbin_final + 1] = {20.0, 40.0, 70.0, 120.0, 180.0, 400.0, 1200.0};
    double reco_l_ptbins_final[reco_nl_ptbin_final + 1] = {20.0,  30.0,  40.0,  55.0,  70.0,  95.0,  120.0,
                                                           150.0, 180.0, 280.0, 400.0, 600.0, 1200.0};
    const int gen_nl_ptbin = gen_nl_ptbin_final * n_sub_bins;
    const int reco_nl_ptbin = reco_nl_ptbin_final * n_sub_bins;
    double gen_l_ptbins[gen_nl_ptbin + 1];
    double reco_l_ptbins[reco_nl_ptbin + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_nl_ptbin_final, gen_nl_ptbin, gen_l_ptbins_final, gen_l_ptbins);
    makeTUnfoldHisto::CreateFineBinningArray(reco_nl_ptbin_final, reco_nl_ptbin, reco_l_ptbins_final, reco_l_ptbins);

    const int gen_nlbar_ptbin_final = 6;
    const int reco_nlbar_ptbin_final = 12;
    double gen_lbar_ptbins_final[gen_nlbar_ptbin_final + 1] = {20.0, 40.0, 70.0, 120.0, 180.0, 400.0, 1200.0};
    double reco_lbar_ptbins_final[reco_nlbar_ptbin_final + 1] = {20.0,  30.0,  40.0,  55.0,  70.0,  95.0,  120.0,
                                                                 150.0, 180.0, 280.0, 400.0, 600.0, 1200.0};
    const int gen_nlbar_ptbin = gen_nlbar_ptbin_final * n_sub_bins;
    const int reco_nlbar_ptbin = reco_nlbar_ptbin_final * n_sub_bins;
    double gen_lbar_ptbins[gen_nlbar_ptbin + 1];
    double reco_lbar_ptbins[reco_nlbar_ptbin + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_nlbar_ptbin_final, gen_nlbar_ptbin, gen_lbar_ptbins_final,
                                             gen_lbar_ptbins);
    makeTUnfoldHisto::CreateFineBinningArray(reco_nlbar_ptbin_final, reco_nlbar_ptbin, reco_lbar_ptbins_final,
                                             reco_lbar_ptbins);

    const int gen_nttbar_ptbin_final = 5;
    const int reco_nttbar_ptbin_final = 10;
    double gen_ttbar_ptbins_final[gen_nttbar_ptbin_final + 1] = {0.0, 30.0, 80.0, 170.0, 300.0, 1200.0};
    double reco_ttbar_ptbins_final[reco_nttbar_ptbin_final + 1] = {0.0,   15.0,  30.0,  55.0,  80.0,  120.0,
                                                                   170.0, 225.0, 300.0, 500.0, 1200.0};
    const int gen_nttbar_ptbin = gen_nttbar_ptbin_final * n_sub_bins;
    const int reco_nttbar_ptbin = reco_nttbar_ptbin_final * n_sub_bins;
    double gen_ttbar_ptbins[gen_nttbar_ptbin + 1];
    double reco_ttbar_ptbins[reco_nttbar_ptbin + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_nttbar_ptbin_final, gen_nttbar_ptbin, gen_ttbar_ptbins_final,
                                             gen_ttbar_ptbins);
    makeTUnfoldHisto::CreateFineBinningArray(reco_nttbar_ptbin_final, reco_nttbar_ptbin, reco_ttbar_ptbins_final,
                                             reco_ttbar_ptbins);

    const int gen_nllbar_ptbin_final = 8;
    const int reco_nllbar_ptbin_final = 16;
    double gen_llbar_ptbins_final[gen_nllbar_ptbin_final + 1] = {0.0,   10.0,  20.0,  40.0, 60.0,
                                                                 100.0, 150.0, 400.0, 457.0};
    double reco_llbar_ptbins_final[reco_nllbar_ptbin_final + 1] = {
        0.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 125.0, 150.0, 275.0, 400.0, 455.0, 457.0};
    const int gen_nllbar_ptbin = gen_nllbar_ptbin_final * n_sub_bins;
    const int reco_nllbar_ptbin = reco_nllbar_ptbin_final * n_sub_bins;
    double gen_llbar_ptbins[gen_nllbar_ptbin + 1];
    double reco_llbar_ptbins[reco_nllbar_ptbin + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_nllbar_ptbin_final, gen_nllbar_ptbin, gen_llbar_ptbins_final,
                                             gen_llbar_ptbins);
    makeTUnfoldHisto::CreateFineBinningArray(reco_nllbar_ptbin_final, reco_nllbar_ptbin, reco_llbar_ptbins_final,
                                             reco_llbar_ptbins);

    const int gen_nllbar_massbin_final = 9;
    const int reco_nllbar_massbin_final = 18;
    // double gen_llbar_massbins_final[gen_nllbar_massbin_final+1]   =
    // {20.0, 30.0, 50.0, 76.0, 106.0, 130.0, 170.0, 260.0, 400.0, 600.0};
    // double reco_llbar_massbins_final[reco_nllbar_massbin_final+1] =
    // {20.0, 25.0, 30.0, 40.0, 50.0, 63.0, 76.0, 90.0, 106.0, 118.0, 130.0,
    // 150.0, 170.0, 215.0, 260.0, 320.0, 400.0, 480.0, 600.0};
    double gen_llbar_massbins_final[gen_nllbar_massbin_final + 1] = {20.0,  30.0,  50.0,  76.0,  106.0,
                                                                     130.0, 170.0, 260.0, 400.0, 447.5};
    double reco_llbar_massbins_final[reco_nllbar_massbin_final + 1] = {20.0,  25.0,  30.0,  40.0,  50.0,  63.0,  76.0,
                                                                       90.0,  106.0, 118.0, 130.0, 150.0, 170.0, 215.0,
                                                                       260.0, 320.0, 400.0, 445.0, 447.5};
    const int gen_nllbar_massbin = gen_nllbar_massbin_final * n_sub_bins;
    const int reco_nllbar_massbin = reco_nllbar_massbin_final * n_sub_bins;
    double gen_llbar_massbins[gen_nllbar_massbin + 1];
    double reco_llbar_massbins[reco_nllbar_massbin + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_nllbar_massbin_final, gen_nllbar_massbin, gen_llbar_massbins_final,
                                             gen_llbar_massbins);
    makeTUnfoldHisto::CreateFineBinningArray(reco_nllbar_massbin_final, reco_nllbar_massbin, reco_llbar_massbins_final,
                                             reco_llbar_massbins);

    const int gen_ndphibin_final = 8;
    const int reco_ndphibin_final = 16;
    // double gen_dphibins_final[gen_ndphibin_final+1]   =
    // {0.0, 1.0, 1.6, 2.0, 2.4, 2.6, 2.8, 2.98, TMath::Pi()}; double
    // reco_dphibins_final[reco_ndphibin_final+1] = {0.0,
    // 0.5, 1.0, 1.3, 1.6, 1.8, 2.0, 2.2, 2.4, 2.5, 2.6, 2.7, 2.8, 2.89, 2.98, 3.07,
    // TMath::Pi()};
    double gen_dphibins_final[gen_ndphibin_final + 1] = {0.0, 1.0, 1.6, 2.0, 2.4, 2.6, 2.8, 3.0, 3.2};
    double reco_dphibins_final[reco_ndphibin_final + 1] = {0.0, 0.5, 1.0, 1.3, 1.6, 1.8, 2.0, 2.2, 2.4,
                                                           2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2};
    const int gen_ndphibin = gen_ndphibin_final * n_sub_bins;
    const int reco_ndphibin = reco_ndphibin_final * n_sub_bins;
    double gen_dphibins[gen_ndphibin + 1];
    double reco_dphibins[reco_ndphibin + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_ndphibin_final, gen_ndphibin, gen_dphibins_final, gen_dphibins);
    makeTUnfoldHisto::CreateFineBinningArray(reco_ndphibin_final, reco_ndphibin, reco_dphibins_final, reco_dphibins);

    // Amandeep : this one is used for ttbar deta by DESY
    const int gen_ndetabin_final = 8;
    const int reco_ndetabin_final = 16;
    double gen_detabins_final[gen_ndetabin_final + 1] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0, 5.0, 10.0};
    double reco_detabins_final[reco_ndetabin_final + 1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6,
                                                           1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0};
    const int gen_ndetabin = gen_ndetabin_final * n_sub_bins;
    const int reco_ndetabin = reco_ndetabin_final * n_sub_bins;
    double gen_detabins[gen_ndetabin + 1];
    double reco_detabins[reco_ndetabin + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_ndetabin_final, gen_ndetabin, gen_detabins_final, gen_detabins);
    makeTUnfoldHisto::CreateFineBinningArray(reco_ndetabin_final, reco_ndetabin, reco_detabins_final, reco_detabins);

    const int gen_netabin_final = 10;
    const int reco_netabin_final = 20;
    double gen_etabins_final[gen_netabin_final + 1] = {-3.2, -2.0, -1.4, -0.8, -0.4, 0.0, 0.4, 0.8, 1.4, 2.0, 3.2};
    double reco_etabins_final[reco_netabin_final + 1] = {-3.2, -2.6, -2.0, -1.7, -1.4, -1.1, -0.8,
                                                         -0.6, -0.4, -0.2, 0.0,  0.2,  0.4,  0.6,
                                                         0.8,  1.1,  1.4,  1.7,  2.0,  2.6,  3.2};
    const int gen_netabin = gen_netabin_final * n_sub_bins;
    const int reco_netabin = reco_netabin_final * n_sub_bins;
    double gen_etabins[gen_netabin + 1];
    double reco_etabins[reco_netabin + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_netabin_final, gen_netabin, gen_etabins_final, gen_etabins);
    makeTUnfoldHisto::CreateFineBinningArray(reco_netabin_final, reco_netabin, reco_etabins_final, reco_etabins);

    const int gen_nbbin_final = 6;
    const int reco_nbbin_final = 12;
    double gen_bbins_final[gen_nbbin_final + 1] = {-1.0 * 2.,    -2. / 3. * 2., -1. / 3. * 2., 0. * 2.,
                                                   1. / 3. * 2., 2. / 3. * 2.,  1.0 * 2.};
    double reco_bbins_final[reco_nbbin_final + 1] = {
        -1.0 * 2.,    -5. / 6. * 2., -2. / 3. * 2., -1. / 2. * 2., -1. / 3. * 2., -1. / 6. * 2., 0. * 2.,
        1. / 6. * 2., 1. / 3. * 2.,  1. / 2. * 2.,  2. / 3. * 2.,  5. / 6. * 2.,  1.0 * 2.};
    const int gen_nbbin = gen_nbbin_final * n_sub_bins;
    const int reco_nbbin = reco_nbbin_final * n_sub_bins;
    double gen_bbins[gen_nbbin + 1];
    double reco_bbins[reco_nbbin + 1];
    makeTUnfoldHisto::CreateFineBinningArray(gen_nbbin_final, gen_nbbin, gen_bbins_final, gen_bbins);
    makeTUnfoldHisto::CreateFineBinningArray(reco_nbbin_final, reco_nbbin, reco_bbins_final, reco_bbins);

    // const int gen_nllbar_delta_phibin_final = 12;
    // const int reco_nllbar_delta_phibin_final = 24;
    // double gen_llbar_delta_phibins_final[gen_nllbar_delta_phibin_final+1] =
    // {0., 1.*TMath::Pi()/12., 2.*TMath::Pi()/12., 3.*TMath::Pi()/12., 4.*TMath::Pi()/12.,
    // 5.*TMath::Pi()/12., 6.*TMath::Pi()/12., 7.*TMath::Pi()/12., 8.*TMath::Pi()/12.,
    // 9.*TMath::Pi()/12., 10.*TMath::Pi()/12., 11.*TMath::Pi()/12.,
    // TMath::Pi()}; double
    // reco_llbar_delta_phibins_final[reco_nllbar_delta_phibin_final+1] =
    // {0., 1.*TMath::Pi()/24., 2.*TMath::Pi()/24., 3.*TMath::Pi()/24., 4.*TMath::Pi()/24.,
    // 5.*TMath::Pi()/24., 6.*TMath::Pi()/24., 7.*TMath::Pi()/24., 8.*TMath::Pi()/24.,
    // 9.*TMath::Pi()/24., 10.*TMath::Pi()/24., 11.*TMath::Pi()/24., 12.*TMath::Pi()/24.,
    // 13.*TMath::Pi()/24., 14.*TMath::Pi()/24., 15.*TMath::Pi()/24., 16.*TMath::Pi()/24.,
    // 17.*TMath::Pi()/24., 18.*TMath::Pi()/24., 19.*TMath::Pi()/24., 20.*TMath::Pi()/24.,
    // 21.*TMath::Pi()/24., 22.*TMath::Pi()/24., 23.*TMath::Pi()/24.,
    // TMath::Pi()}; for now keep things simple by sticking to 6 bins for all
    // the spin variables

    // Make files related to binning for each of the variables
    // Each variable needs its own binning scheme even if it is the exact same
    // binning scheme as another variable This is because currently each
    // variable's binning scheme is saved in its own .xml file and read when
    // performing unfolding!
    // TODO: Perhaps change the above to allow for linking between variables to
    // the same binning scheme to save space/cleaner code

    gSystem->mkdir("binning", true);

    // ******************************
    // ******************************
    // 1D binning XML file generation
    // ******************************
    // ******************************

    // Kinematic
    std::pair<TUnfoldBinning *, TUnfoldBinning *> top_pt_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "top_pt", std::vector<std::string>{"top_pt"}, std::vector<int>{gen_ntop_ptbin},
        std::vector<int>{reco_ntop_ptbin}, std::vector<double *>{gen_top_ptbins},
        std::vector<double *>{reco_top_ptbins}, n_sub_bins);
    residual_top_pt_binning = CreateResidualBinning(nbinsres, reco_top_ptbins[0] - reco_top_ptbins[reco_ntop_ptbin],
                                                    reco_top_ptbins[reco_ntop_ptbin] - reco_top_ptbins[0]);
    gen_top_pt_binning = top_pt_binnings.first;
    reco_top_pt_binning = top_pt_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> l_pt_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "l_pt", std::vector<std::string>{"l_pt"}, std::vector<int>{gen_nl_ptbin}, std::vector<int>{reco_nl_ptbin},
        std::vector<double *>{gen_l_ptbins}, std::vector<double *>{reco_l_ptbins}, n_sub_bins);
    residual_l_pt_binning = CreateResidualBinning(nbinsres, reco_l_ptbins[0] - reco_l_ptbins[reco_nl_ptbin],
                                                  reco_l_ptbins[reco_nl_ptbin] - reco_l_ptbins[0]);
    gen_l_pt_binning = l_pt_binnings.first;
    reco_l_pt_binning = l_pt_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> lbar_pt_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "lbar_pt", std::vector<std::string>{"lbar_pt"}, std::vector<int>{gen_nlbar_ptbin},
        std::vector<int>{reco_nlbar_ptbin}, std::vector<double *>{gen_lbar_ptbins},
        std::vector<double *>{reco_lbar_ptbins}, n_sub_bins);
    residual_lbar_pt_binning = CreateResidualBinning(nbinsres, reco_lbar_ptbins[0] - reco_lbar_ptbins[reco_nlbar_ptbin],
                                                     reco_lbar_ptbins[reco_nlbar_ptbin] - reco_lbar_ptbins[0]);
    gen_lbar_pt_binning = lbar_pt_binnings.first;
    reco_lbar_pt_binning = lbar_pt_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> ttbar_pt_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "ttbar_pt", std::vector<std::string>{"ttbar_pt"}, std::vector<int>{gen_nttbar_ptbin},
        std::vector<int>{reco_nttbar_ptbin}, std::vector<double *>{gen_ttbar_ptbins},
        std::vector<double *>{reco_ttbar_ptbins}, n_sub_bins);
    residual_ttbar_pt_binning =
        CreateResidualBinning(nbinsres, reco_ttbar_ptbins[0] - reco_ttbar_ptbins[reco_nttbar_ptbin],
                              reco_ttbar_ptbins[reco_nttbar_ptbin] - reco_ttbar_ptbins[0]);
    gen_ttbar_pt_binning = ttbar_pt_binnings.first;
    reco_ttbar_pt_binning = ttbar_pt_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> ttbar_mass_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "ttbar_mass", std::vector<std::string>{"ttbar_mass"}, std::vector<int>{gen_nttbar_massbin},
        std::vector<int>{reco_nttbar_massbin}, std::vector<double *>{gen_ttbar_massbins},
        std::vector<double *>{reco_ttbar_massbins}, n_sub_bins);
    residual_ttbar_mass_binning =
        CreateResidualBinning(nbinsres, reco_ttbar_massbins[0] - reco_ttbar_massbins[reco_nttbar_massbin],
                              reco_ttbar_massbins[reco_nttbar_massbin] - reco_ttbar_massbins[0]);
    gen_ttbar_mass_binning = ttbar_mass_binnings.first;
    reco_ttbar_mass_binning = ttbar_mass_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> ttbar_delta_phi_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "ttbar_delta_phi", std::vector<std::string>{"ttbar_delta_phi"}, std::vector<int>{gen_ndphibin},
        std::vector<int>{reco_ndphibin}, std::vector<double *>{gen_dphibins}, std::vector<double *>{reco_dphibins},
        n_sub_bins);
    residual_ttbar_delta_phi_binning = CreateResidualBinning(nbinsres, reco_dphibins[0] - reco_dphibins[reco_ndphibin],
                                                             reco_dphibins[reco_ndphibin] - reco_dphibins[0]);
    gen_ttbar_delta_phi_binning = ttbar_delta_phi_binnings.first;
    reco_ttbar_delta_phi_binning = ttbar_delta_phi_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> ttbar_delta_eta_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "ttbar_delta_eta", std::vector<std::string>{"ttbar_delta_eta"}, std::vector<int>{gen_ndetabin},
        std::vector<int>{reco_ndetabin}, std::vector<double *>{gen_detabins}, std::vector<double *>{reco_detabins},
        n_sub_bins);
    residual_ttbar_delta_eta_binning = CreateResidualBinning(nbinsres, reco_detabins[0] - reco_detabins[reco_ndetabin],
                                                             reco_detabins[reco_ndetabin] - reco_detabins[0]);
    gen_ttbar_delta_eta_binning = ttbar_delta_eta_binnings.first;
    reco_ttbar_delta_eta_binning = ttbar_delta_eta_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> ttbar_rapidity_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "ttbar_rapidity", std::vector<std::string>{"ttbar_rapidity"}, std::vector<int>{gen_netabin},
        std::vector<int>{reco_netabin}, std::vector<double *>{gen_etabins}, std::vector<double *>{reco_etabins},
        n_sub_bins);
    residual_ttbar_rapidity_binning = CreateResidualBinning(nbinsres, reco_etabins[0] - reco_etabins[reco_netabin],
                                                            reco_etabins[reco_netabin] - reco_etabins[0]);
    gen_ttbar_rapidity_binning = ttbar_rapidity_binnings.first;
    reco_ttbar_rapidity_binning = ttbar_rapidity_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> llbar_pt_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "llbar_pt", std::vector<std::string>{"llbar_pt"}, std::vector<int>{gen_nllbar_ptbin},
        std::vector<int>{reco_nllbar_ptbin}, std::vector<double *>{gen_llbar_ptbins},
        std::vector<double *>{reco_llbar_ptbins}, n_sub_bins);
    residual_llbar_pt_binning =
        CreateResidualBinning(nbinsres, reco_llbar_ptbins[0] - reco_llbar_ptbins[reco_nllbar_ptbin],
                              reco_llbar_ptbins[reco_nllbar_ptbin] - reco_llbar_ptbins[0]);
    gen_llbar_pt_binning = llbar_pt_binnings.first;
    reco_llbar_pt_binning = llbar_pt_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> llbar_mass_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "llbar_mass", std::vector<std::string>{"llbar_mass"}, std::vector<int>{gen_nllbar_massbin},
        std::vector<int>{reco_nllbar_massbin}, std::vector<double *>{gen_llbar_massbins},
        std::vector<double *>{reco_llbar_massbins}, n_sub_bins);
    residual_llbar_mass_binning =
        CreateResidualBinning(nbinsres, reco_llbar_massbins[0] - reco_llbar_massbins[reco_nllbar_massbin],
                              reco_llbar_massbins[reco_nllbar_massbin] - reco_llbar_massbins[0]);
    gen_llbar_mass_binning = llbar_mass_binnings.first;
    reco_llbar_mass_binning = llbar_mass_binnings.second;

    // Spin corr
    std::pair<TUnfoldBinning *, TUnfoldBinning *> b1k_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b1k", std::vector<std::string>{"b1k"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_b1k_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_b1k_binning = b1k_binnings.first;
    reco_b1k_binning = b1k_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b2k_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b2k", std::vector<std::string>{"b2k"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_b2k_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_b2k_binning = b2k_binnings.first;
    reco_b2k_binning = b2k_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b1j_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b1j", std::vector<std::string>{"b1j"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_b1j_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_b1j_binning = b1j_binnings.first;
    reco_b1j_binning = b1j_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b2j_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b2j", std::vector<std::string>{"b2j"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_b2j_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_b2j_binning = b2j_binnings.first;
    reco_b2j_binning = b2j_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b1r_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b1r", std::vector<std::string>{"b1r"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_b1r_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_b1r_binning = b1r_binnings.first;
    reco_b1r_binning = b1r_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b2r_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b2r", std::vector<std::string>{"b2r"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_b2r_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_b2r_binning = b2r_binnings.first;
    reco_b2r_binning = b2r_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b1q_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b1q", std::vector<std::string>{"b1q"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_b1q_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_b1q_binning = b1q_binnings.first;
    reco_b1q_binning = b1q_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b2q_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b2q", std::vector<std::string>{"b2q"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_b2q_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_b2q_binning = b2q_binnings.first;
    reco_b2q_binning = b2q_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b1n_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b1n", std::vector<std::string>{"b1n"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_b1n_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_b1n_binning = b1n_binnings.first;
    reco_b1n_binning = b1n_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b2n_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b2n", std::vector<std::string>{"b2n"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_b2n_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_b2n_binning = b2n_binnings.first;
    reco_b2n_binning = b2n_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_kk_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_kk", std::vector<std::string>{"c_kk"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_kk_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_kk_binning = c_kk_binnings.first;
    reco_c_kk_binning = c_kk_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rr_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_rr", std::vector<std::string>{"c_rr"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_rr_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_rr_binning = c_rr_binnings.first;
    reco_c_rr_binning = c_rr_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nn_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_nn", std::vector<std::string>{"c_nn"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_nn_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_nn_binning = c_nn_binnings.first;
    reco_c_nn_binning = c_nn_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_kj_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_kj", std::vector<std::string>{"c_kj"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_kj_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_kj_binning = c_kj_binnings.first;
    reco_c_kj_binning = c_kj_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rq_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_rq", std::vector<std::string>{"c_rq"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_rq_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_rq_binning = c_rq_binnings.first;
    reco_c_rq_binning = c_rq_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Prk_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_Prk", std::vector<std::string>{"c_Prk"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_Prk_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_Prk_binning = c_Prk_binnings.first;
    reco_c_Prk_binning = c_Prk_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Mrk_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_Mrk", std::vector<std::string>{"c_Mrk"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_Mrk_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_Mrk_binning = c_Mrk_binnings.first;
    reco_c_Mrk_binning = c_Mrk_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Pnr_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_Pnr", std::vector<std::string>{"c_Pnr"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_Pnr_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_Pnr_binning = c_Pnr_binnings.first;
    reco_c_Pnr_binning = c_Pnr_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Mnr_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_Mnr", std::vector<std::string>{"c_Mnr"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_Mnr_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_Mnr_binning = c_Mnr_binnings.first;
    reco_c_Mnr_binning = c_Mnr_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Pnk_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_Pnk", std::vector<std::string>{"c_Pnk"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_Pnk_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_Pnk_binning = c_Pnk_binnings.first;
    reco_c_Pnk_binning = c_Pnk_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Mnk_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_Mnk", std::vector<std::string>{"c_Mnk"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_Mnk_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_Mnk_binning = c_Mnk_binnings.first;
    reco_c_Mnk_binning = c_Mnk_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rk_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_rk", std::vector<std::string>{"c_rk"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_rk_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_rk_binning = c_rk_binnings.first;
    reco_c_rk_binning = c_rk_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_kr_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_kr", std::vector<std::string>{"c_kr"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_kr_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_kr_binning = c_kr_binnings.first;
    reco_c_kr_binning = c_kr_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nr_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_nr", std::vector<std::string>{"c_nr"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_nr_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_nr_binning = c_nr_binnings.first;
    reco_c_nr_binning = c_nr_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rn_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_rn", std::vector<std::string>{"c_rn"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_rn_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_rn_binning = c_rn_binnings.first;
    reco_c_rn_binning = c_rn_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nk_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_nk", std::vector<std::string>{"c_nk"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_nk_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_nk_binning = c_nk_binnings.first;
    reco_c_nk_binning = c_nk_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_kn_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_kn", std::vector<std::string>{"c_kn"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_kn_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_kn_binning = c_kn_binnings.first;
    reco_c_kn_binning = c_kn_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_han_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_han", std::vector<std::string>{"c_han"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_han_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_han_binning = c_han_binnings.first;
    reco_c_han_binning = c_han_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_sca_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_sca", std::vector<std::string>{"c_sca"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_sca_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_sca_binning = c_sca_binnings.first;
    reco_c_sca_binning = c_sca_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_tra_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_tra", std::vector<std::string>{"c_tra"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_tra_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_tra_binning = c_tra_binnings.first;
    reco_c_tra_binning = c_tra_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_kjL_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_kjL", std::vector<std::string>{"c_kjL"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_kjL_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_kjL_binning = c_kjL_binnings.first;
    reco_c_kjL_binning = c_kjL_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rqL_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_rqL", std::vector<std::string>{"c_rqL"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_rqL_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_rqL_binning = c_rqL_binnings.first;
    reco_c_rqL_binning = c_rqL_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rkP_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_rkP", std::vector<std::string>{"c_rkP"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_rkP_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_rkP_binning = c_rkP_binnings.first;
    reco_c_rkP_binning = c_rkP_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rkM_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_rkM", std::vector<std::string>{"c_rkM"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_rkM_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_rkM_binning = c_rkM_binnings.first;
    reco_c_rkM_binning = c_rkM_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nrP_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_nrP", std::vector<std::string>{"c_nrP"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_nrP_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_nrP_binning = c_nrP_binnings.first;
    reco_c_nrP_binning = c_nrP_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nrM_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_nrM", std::vector<std::string>{"c_nrM"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_nrM_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_nrM_binning = c_nrM_binnings.first;
    reco_c_nrM_binning = c_nrM_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nkP_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_nkP", std::vector<std::string>{"c_nkP"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_nkP_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_nkP_binning = c_nkP_binnings.first;
    reco_c_nkP_binning = c_nkP_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nkM_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_nkM", std::vector<std::string>{"c_nkM"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_nkM_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_nkM_binning = c_nkM_binnings.first;
    reco_c_nkM_binning = c_nkM_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> ll_cHel_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "ll_cHel", std::vector<std::string>{"ll_cHel"}, std::vector<int>{gen_ncosbin}, std::vector<int>{reco_ncosbin},
        std::vector<double *>{gen_cosbins}, std::vector<double *>{reco_cosbins}, n_sub_bins);
    residual_ll_cHel_binning = CreateResidualBinning(nbinsres, reco_cosbins[0] - reco_cosbins[reco_ncosbin],
                                                     reco_cosbins[reco_ncosbin] - reco_cosbins[0]);
    gen_ll_cHel_binning = ll_cHel_binnings.first;
    reco_ll_cHel_binning = ll_cHel_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> ll_cLab_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "ll_cLab", std::vector<std::string>{"ll_cLab"}, std::vector<int>{gen_ncosbin}, std::vector<int>{reco_ncosbin},
        std::vector<double *>{gen_cosbins}, std::vector<double *>{reco_cosbins}, n_sub_bins);
    residual_ll_cLab_binning = CreateResidualBinning(nbinsres, reco_cosbins[0] - reco_cosbins[reco_ncosbin],
                                                     reco_cosbins[reco_ncosbin] - reco_cosbins[0]);
    gen_ll_cLab_binning = ll_cLab_binnings.first;
    reco_ll_cLab_binning = ll_cLab_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> ll_kNorm_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "ll_kNorm", std::vector<std::string>{"ll_kNorm"}, std::vector<int>{gen_ncosbin}, std::vector<int>{reco_ncosbin},
        std::vector<double *>{gen_cosbins}, std::vector<double *>{reco_cosbins}, n_sub_bins);
    residual_ll_kNorm_binning = CreateResidualBinning(nbinsres, reco_cosbins[0] - reco_cosbins[reco_ncosbin],
                                                      reco_cosbins[reco_ncosbin] - reco_cosbins[0]);
    gen_ll_kNorm_binning = ll_kNorm_binnings.first;
    reco_ll_kNorm_binning = ll_kNorm_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> ll_rNorm_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "ll_rNorm", std::vector<std::string>{"ll_rNorm"}, std::vector<int>{gen_ncosbin}, std::vector<int>{reco_ncosbin},
        std::vector<double *>{gen_cosbins}, std::vector<double *>{reco_cosbins}, n_sub_bins);
    residual_ll_rNorm_binning = CreateResidualBinning(nbinsres, reco_cosbins[0] - reco_cosbins[reco_ncosbin],
                                                      reco_cosbins[reco_ncosbin] - reco_cosbins[0]);
    gen_ll_rNorm_binning = ll_rNorm_binnings.first;
    reco_ll_rNorm_binning = ll_rNorm_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> llbar_delta_phi_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "llbar_delta_phi", std::vector<std::string>{"llbar_delta_phi"}, std::vector<int>{gen_nllbar_delta_phibin},
        std::vector<int>{reco_nllbar_delta_phibin}, std::vector<double *>{gen_llbar_delta_phibins},
        std::vector<double *>{reco_llbar_delta_phibins}, n_sub_bins);
    residual_llbar_delta_phi_binning = CreateResidualBinning(
        nbinsres, reco_llbar_delta_phibins[0] - reco_llbar_delta_phibins[reco_nllbar_delta_phibin],
        reco_llbar_delta_phibins[reco_nllbar_delta_phibin] - reco_llbar_delta_phibins[0]);
    gen_llbar_delta_phi_binning = llbar_delta_phi_binnings.first;
    reco_llbar_delta_phi_binning = llbar_delta_phi_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> llbar_delta_eta_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "llbar_delta_eta", std::vector<std::string>{"llbar_delta_eta"}, std::vector<int>{gen_nllbar_delta_etabin},
        std::vector<int>{reco_nllbar_delta_etabin}, std::vector<double *>{gen_llbar_delta_etabins},
        std::vector<double *>{reco_llbar_delta_etabins}, n_sub_bins);
    residual_llbar_delta_eta_binning = CreateResidualBinning(
        nbinsres, reco_llbar_delta_etabins[0] - reco_llbar_delta_etabins[reco_nllbar_delta_etabin],
        reco_llbar_delta_etabins[reco_nllbar_delta_etabin] - reco_llbar_delta_etabins[0]);
    gen_llbar_delta_eta_binning = llbar_delta_eta_binnings.first;
    reco_llbar_delta_eta_binning = llbar_delta_eta_binnings.second;

    // Amandeep : adding starred variable binnings here
    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rj_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_rj", std::vector<std::string>{"c_rj"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_rj_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_rj_binning = c_rj_binnings.first;
    reco_c_rj_binning = c_rj_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_jr_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_jr", std::vector<std::string>{"c_jr"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_jr_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_jr_binning = c_jr_binnings.first;
    reco_c_jr_binning = c_jr_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Prj_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_Prj", std::vector<std::string>{"c_Prj"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_Prj_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_Prj_binning = c_Prj_binnings.first;
    reco_c_Prj_binning = c_Prj_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Mrj_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_Mrj", std::vector<std::string>{"c_Mrj"}, std::vector<int>{gen_ncbin}, std::vector<int>{reco_ncbin},
        std::vector<double *>{gen_cbins}, std::vector<double *>{reco_cbins}, n_sub_bins);
    residual_c_Mrj_binning =
        CreateResidualBinning(nbinsres, reco_cbins[0] - reco_cbins[reco_ncbin], reco_cbins[reco_ncbin] - reco_cbins[0]);
    gen_c_Mrj_binning = c_Mrj_binnings.first;
    reco_c_Mrj_binning = c_Mrj_binnings.second;
    // End

    // ******************************
    // ******************************
    // 2D binning XML file generation
    // ******************************
    // ******************************

    // Amandeep : 2D XML stuff
    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_kk_mttbar_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_kk_mttbar",
    // std::vector<std::string>{"c_kk", "mttbar"}, std::vector<int>{gen_ncbin,
    // gen_nttbar_massbin}, std::vector<int>{reco_ncbin, reco_nttbar_massbin},
    // std::vector<double *>{gen_cbins, gen_ttbar_massbins}, std::vector<double
    // *>{reco_cbins, reco_ttbar_massbins}, n_sub_bins);
    // residual_c_kk_mttbar_binning = CreateResidualBinning(nbinsres,
    // reco_cbins[0]-reco_cbins[reco_ncbin],
    // reco_cbins[reco_ncbin]-reco_cbins[0]); gen_c_kk_mttbar_binning      =
    // c_kk_mttbar_binnings.first; reco_c_kk_mttbar_binning     =
    // c_kk_mttbar_binnings.second;

    // ******
    // mttbar
    // ******

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b1k_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b1k_exjet", std::vector<std::string>{"b1k", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_b1k_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                       reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_b1k_exjet_binning = b1k_exjet_binnings.first;
    reco_b1k_exjet_binning = b1k_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b2k_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b2k_exjet", std::vector<std::string>{"b2k ", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_b2k_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                       reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_b2k_exjet_binning = b2k_exjet_binnings.first;
    reco_b2k_exjet_binning = b2k_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b1j_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b1j_exjet", std::vector<std::string>{"b1j", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_b1j_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                       reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_b1j_exjet_binning = b1j_exjet_binnings.first;
    reco_b1j_exjet_binning = b1j_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b2j_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b2j_exjet", std::vector<std::string>{"b2j", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_b2j_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                       reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_b2j_exjet_binning = b2j_exjet_binnings.first;
    reco_b2j_exjet_binning = b2j_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b1r_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b1r_exjet", std::vector<std::string>{"b1r", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_b1r_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                       reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_b1r_exjet_binning = b1r_exjet_binnings.first;
    reco_b1r_exjet_binning = b1r_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b2r_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b2r_exjet", std::vector<std::string>{"b2r", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_b2r_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                       reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_b2r_exjet_binning = b2r_exjet_binnings.first;
    reco_b2r_exjet_binning = b2r_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b1q_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b1q_exjet", std::vector<std::string>{"b1q", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_b1q_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                       reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_b1q_exjet_binning = b1q_exjet_binnings.first;
    reco_b1q_exjet_binning = b1q_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b2q_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b2q_exjet", std::vector<std::string>{"b2q", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_b2q_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                       reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_b2q_exjet_binning = b2q_exjet_binnings.first;
    reco_b2q_exjet_binning = b2q_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b1n_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b1n_exjet", std::vector<std::string>{"b1n", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_b1n_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                       reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_b1n_exjet_binning = b1n_exjet_binnings.first;
    reco_b1n_exjet_binning = b1n_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> b2n_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "b2n_exjet", std::vector<std::string>{"b2n", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_b2n_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                       reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_b2n_exjet_binning = b2n_exjet_binnings.first;
    reco_b2n_exjet_binning = b2n_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_kk_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_kk_exjet", std::vector<std::string>{"c_kk", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_kk_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                        reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_kk_exjet_binning = c_kk_exjet_binnings.first;
    reco_c_kk_exjet_binning = c_kk_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rr_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_rr_exjet", std::vector<std::string>{"c_rr", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_rr_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                        reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_rr_exjet_binning = c_rr_exjet_binnings.first;
    reco_c_rr_exjet_binning = c_rr_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nn_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_nn_exjet", std::vector<std::string>{"c_nn", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_nn_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                        reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_nn_exjet_binning = c_nn_exjet_binnings.first;
    reco_c_nn_exjet_binning = c_nn_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_kj_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_kj_exjet", std::vector<std::string>{"c_kj", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_kj_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                        reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_kj_exjet_binning = c_kj_exjet_binnings.first;
    reco_c_kj_exjet_binning = c_kj_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rq_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_rq_exjet", std::vector<std::string>{"c_rq", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_rq_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                        reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_rq_exjet_binning = c_rq_exjet_binnings.first;
    reco_c_rq_exjet_binning = c_rq_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Prk_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_Prk_exjet", std::vector<std::string>{"c_Prk", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_Prk_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_Prk_exjet_binning = c_Prk_exjet_binnings.first;
    reco_c_Prk_exjet_binning = c_Prk_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Mrk_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_Mrk_exjet", std::vector<std::string>{"c_Mrk", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_Mrk_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_Mrk_exjet_binning = c_Mrk_exjet_binnings.first;
    reco_c_Mrk_exjet_binning = c_Mrk_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Pnr_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_Pnr_exjet", std::vector<std::string>{"c_Pnr", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_Pnr_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_Pnr_exjet_binning = c_Pnr_exjet_binnings.first;
    reco_c_Pnr_exjet_binning = c_Pnr_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Mnr_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_Mnr_exjet", std::vector<std::string>{"c_Mnr", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_Mnr_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_Mnr_exjet_binning = c_Mnr_exjet_binnings.first;
    reco_c_Mnr_exjet_binning = c_Mnr_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Pnk_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_Pnk_exjet", std::vector<std::string>{"c_Pnk", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_Pnk_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_Pnk_exjet_binning = c_Pnk_exjet_binnings.first;
    reco_c_Pnk_exjet_binning = c_Pnk_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Mnk_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_Mnk_exjet", std::vector<std::string>{"c_Mnk", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_Mnk_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_Mnk_exjet_binning = c_Mnk_exjet_binnings.first;
    reco_c_Mnk_exjet_binning = c_Mnk_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rk_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_rk_exjet", std::vector<std::string>{"c_rk", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_rk_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                        reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_rk_exjet_binning = c_rk_exjet_binnings.first;
    reco_c_rk_exjet_binning = c_rk_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_kr_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_kr_exjet", std::vector<std::string>{"c_kr", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_kr_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                        reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_kr_exjet_binning = c_kr_exjet_binnings.first;
    reco_c_kr_exjet_binning = c_kr_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nr_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_nr_exjet", std::vector<std::string>{"c_nr", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_nr_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                        reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_nr_exjet_binning = c_nr_exjet_binnings.first;
    reco_c_nr_exjet_binning = c_nr_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rn_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_rn_exjet", std::vector<std::string>{"c_rn", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_rn_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                        reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_rn_exjet_binning = c_rn_exjet_binnings.first;
    reco_c_rn_exjet_binning = c_rn_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nk_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_nk_exjet", std::vector<std::string>{"c_nk", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_nk_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                        reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_nk_exjet_binning = c_nk_exjet_binnings.first;
    reco_c_nk_exjet_binning = c_nk_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_kn_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_kn_exjet", std::vector<std::string>{"c_kn", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_kn_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                        reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_kn_exjet_binning = c_kn_exjet_binnings.first;
    reco_c_kn_exjet_binning = c_kn_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_han_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_han_exjet", std::vector<std::string>{"c_han", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_han_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_han_exjet_binning = c_han_exjet_binnings.first;
    reco_c_han_exjet_binning = c_han_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_sca_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_sca_exjet", std::vector<std::string>{"c_sca", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_sca_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_sca_exjet_binning = c_sca_exjet_binnings.first;
    reco_c_sca_exjet_binning = c_sca_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_tra_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_tra_exjet", std::vector<std::string>{"c_tra", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_tra_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_tra_exjet_binning = c_tra_exjet_binnings.first;
    reco_c_tra_exjet_binning = c_tra_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_kjL_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_kjL_exjet", std::vector<std::string>{"c_kjL", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_kjL_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_kjL_exjet_binning = c_kjL_exjet_binnings.first;
    reco_c_kjL_exjet_binning = c_kjL_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rqL_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_rqL_exjet", std::vector<std::string>{"c_rqL", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_rqL_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_rqL_exjet_binning = c_rqL_exjet_binnings.first;
    reco_c_rqL_exjet_binning = c_rqL_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rkP_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_rkP_exjet", std::vector<std::string>{"c_rkP", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_rkP_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_rkP_exjet_binning = c_rkP_exjet_binnings.first;
    reco_c_rkP_exjet_binning = c_rkP_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rkM_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_rkM_exjet", std::vector<std::string>{"c_rkM", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_rkM_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_rkM_exjet_binning = c_rkM_exjet_binnings.first;
    reco_c_rkM_exjet_binning = c_rkM_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nrP_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_nrP_exjet", std::vector<std::string>{"c_nrP", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_nrP_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_nrP_exjet_binning = c_nrP_exjet_binnings.first;
    reco_c_nrP_exjet_binning = c_nrP_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nrM_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_nrM_exjet", std::vector<std::string>{"c_nrM", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_nrM_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_nrM_exjet_binning = c_nrM_exjet_binnings.first;
    reco_c_nrM_exjet_binning = c_nrM_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nkP_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_nkP_exjet", std::vector<std::string>{"c_nkP", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_nkP_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_nkP_exjet_binning = c_nkP_exjet_binnings.first;
    reco_c_nkP_exjet_binning = c_nkP_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nkM_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_nkM_exjet", std::vector<std::string>{"c_nkM", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_nkM_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_nkM_exjet_binning = c_nkM_exjet_binnings.first;
    reco_c_nkM_exjet_binning = c_nkM_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> ll_cHel_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "ll_cHel_exjet", std::vector<std::string>{"ll_cHel", "exjet"},
        std::vector<int>{gen_ncosbin_2D, gen_nexjetbin_2D}, std::vector<int>{reco_ncosbin_2D, reco_nexjetbin_2D},
        std::vector<double *>{gen_cosbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cosbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_ll_cHel_exjet_binning =
        CreateResidualBinning(nbinsres, reco_cosbins_2D[0] - reco_cosbins_2D[reco_ncosbin_2D],
                              reco_cosbins_2D[reco_ncosbin_2D] - reco_cosbins_2D[0]);
    gen_ll_cHel_exjet_binning = ll_cHel_exjet_binnings.first;
    reco_ll_cHel_exjet_binning = ll_cHel_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> ll_cLab_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "ll_cLab_exjet", std::vector<std::string>{"ll_cLab", "exjet"},
        std::vector<int>{gen_ncosbin_2D, gen_nexjetbin_2D}, std::vector<int>{reco_ncosbin_2D, reco_nexjetbin_2D},
        std::vector<double *>{gen_cosbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cosbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_ll_cLab_exjet_binning =
        CreateResidualBinning(nbinsres, reco_cosbins_2D[0] - reco_cosbins_2D[reco_ncosbin_2D],
                              reco_cosbins_2D[reco_ncosbin_2D] - reco_cosbins_2D[0]);
    gen_ll_cLab_exjet_binning = ll_cLab_exjet_binnings.first;
    reco_ll_cLab_exjet_binning = ll_cLab_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> ll_kNorm_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "ll_kNorm_exjet", std::vector<std::string>{"ll_kNorm", "exjet"},
        std::vector<int>{gen_ncosbin_2D, gen_nexjetbin_2D}, std::vector<int>{reco_ncosbin_2D, reco_nexjetbin_2D},
        std::vector<double *>{gen_cosbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cosbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_ll_kNorm_exjet_binning =
        CreateResidualBinning(nbinsres, reco_cosbins_2D[0] - reco_cosbins_2D[reco_ncosbin_2D],
                              reco_cosbins_2D[reco_ncosbin_2D] - reco_cosbins_2D[0]);
    gen_ll_kNorm_exjet_binning = ll_kNorm_exjet_binnings.first;
    reco_ll_kNorm_exjet_binning = ll_kNorm_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> ll_rNorm_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "ll_rNorm_exjet", std::vector<std::string>{"ll_rNorm", "exjet"},
        std::vector<int>{gen_ncosbin_2D, gen_nexjetbin_2D}, std::vector<int>{reco_ncosbin_2D, reco_nexjetbin_2D},
        std::vector<double *>{gen_cosbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cosbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_ll_rNorm_exjet_binning =
        CreateResidualBinning(nbinsres, reco_cosbins_2D[0] - reco_cosbins_2D[reco_ncosbin_2D],
                              reco_cosbins_2D[reco_ncosbin_2D] - reco_cosbins_2D[0]);
    gen_ll_rNorm_exjet_binning = ll_rNorm_exjet_binnings.first;
    reco_ll_rNorm_exjet_binning = ll_rNorm_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> llbar_delta_phi_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "llbar_delta_phi_exjet", std::vector<std::string>{"llbar_delta_phi", "exjet"},
        std::vector<int>{gen_nllbar_delta_phibin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_nllbar_delta_phibin_2D, reco_nexjetbin_2D},
        std::vector<double *>{gen_llbar_delta_phibins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_llbar_delta_phibins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_llbar_delta_phi_exjet_binning = CreateResidualBinning(
        nbinsres, reco_llbar_delta_phibins_2D[0] - reco_llbar_delta_phibins_2D[reco_nllbar_delta_phibin_2D],
        reco_llbar_delta_phibins_2D[reco_nllbar_delta_phibin_2D] - reco_llbar_delta_phibins_2D[0]);
    gen_llbar_delta_phi_exjet_binning = llbar_delta_phi_exjet_binnings.first;
    reco_llbar_delta_phi_exjet_binning = llbar_delta_phi_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> llbar_delta_eta_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "llbar_delta_eta_exjet", std::vector<std::string>{"llbar_delta_eta", "exjet"},
        std::vector<int>{gen_nllbar_delta_etabin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_nllbar_delta_etabin_2D, reco_nexjetbin_2D},
        std::vector<double *>{gen_llbar_delta_etabins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_llbar_delta_etabins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_llbar_delta_eta_exjet_binning = CreateResidualBinning(
        nbinsres, reco_llbar_delta_etabins_2D[0] - reco_llbar_delta_etabins_2D[reco_nllbar_delta_etabin_2D],
        reco_llbar_delta_etabins_2D[reco_nllbar_delta_etabin_2D] - reco_llbar_delta_etabins_2D[0]);
    gen_llbar_delta_eta_exjet_binning = llbar_delta_eta_exjet_binnings.first;
    reco_llbar_delta_eta_exjet_binning = llbar_delta_eta_exjet_binnings.second;

    // Amandeep : adding starred variable binnings here
    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rj_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_rj_exjet", std::vector<std::string>{"c_rj", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_rj_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                        reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_rj_exjet_binning = c_rj_exjet_binnings.first;
    reco_c_rj_exjet_binning = c_rj_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_jr_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_jr_exjet", std::vector<std::string>{"c_jr", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_jr_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                        reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_jr_exjet_binning = c_jr_exjet_binnings.first;
    reco_c_jr_exjet_binning = c_jr_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Prj_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_Prj_exjet", std::vector<std::string>{"c_Prj", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_Prj_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_Prj_exjet_binning = c_Prj_exjet_binnings.first;
    reco_c_Prj_exjet_binning = c_Prj_exjet_binnings.second;

    std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Mrj_exjet_binnings = makeTUnfoldHisto::CreateXMLBinFiles(
        "c_Mrj_exjet", std::vector<std::string>{"c_Mrj", "exjet"}, std::vector<int>{gen_ncbin_2D, gen_nexjetbin_2D},
        std::vector<int>{reco_ncbin_2D, reco_nexjetbin_2D}, std::vector<double *>{gen_cbins_2D, gen_exjetbins_2D},
        std::vector<double *>{reco_cbins_2D, reco_exjetbins_2D}, n_sub_bins_2D);
    residual_c_Mrj_exjet_binning = CreateResidualBinning(nbinsres, reco_cbins_2D[0] - reco_cbins_2D[reco_ncbin_2D],
                                                         reco_cbins_2D[reco_ncbin_2D] - reco_cbins_2D[0]);
    gen_c_Mrj_exjet_binning = c_Mrj_exjet_binnings.first;
    reco_c_Mrj_exjet_binning = c_Mrj_exjet_binnings.second;

    // ******
    // top_pt
    // ******

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> b1k_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b1k_top_pt",
    // std::vector<std::string>{"b1k", "top_pt"}, std::vector<int>{gen_ncbin_2D,
    // gen_ntop_ptbin_2D}, std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D},
    // std::vector<double *>{gen_cbins_2D, gen_top_ptbins_2D},
    // std::vector<double *>{reco_cbins_2D, reco_top_ptbins_2D}, n_sub_bins_2D);
    // residual_b1k_top_pt_binning = CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_b1k_top_pt_binning =
    // b1k_top_pt_binnings.first; reco_b1k_top_pt_binning     =
    // b1k_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> b2k_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b2k_top_pt",
    // std::vector<std::string>{"b2k ", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_b2k_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_b2k_top_pt_binning =
    // b2k_top_pt_binnings.first; reco_b2k_top_pt_binning     =
    // b2k_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> b1j_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b1j_top_pt",
    // std::vector<std::string>{"b1j", "top_pt"}, std::vector<int>{gen_ncbin_2D,
    // gen_ntop_ptbin_2D}, std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D},
    // std::vector<double *>{gen_cbins_2D, gen_top_ptbins_2D},
    // std::vector<double *>{reco_cbins_2D, reco_top_ptbins_2D}, n_sub_bins_2D);
    // residual_b1j_top_pt_binning = CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_b1j_top_pt_binning =
    // b1j_top_pt_binnings.first; reco_b1j_top_pt_binning     =
    // b1j_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> b2j_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b2j_top_pt",
    // std::vector<std::string>{"b2j", "top_pt"}, std::vector<int>{gen_ncbin_2D,
    // gen_ntop_ptbin_2D}, std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D},
    // std::vector<double *>{gen_cbins_2D, gen_top_ptbins_2D},
    // std::vector<double *>{reco_cbins_2D, reco_top_ptbins_2D}, n_sub_bins_2D);
    // residual_b2j_top_pt_binning = CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_b2j_top_pt_binning =
    // b2j_top_pt_binnings.first; reco_b2j_top_pt_binning     =
    // b2j_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> b1r_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b1r_top_pt",
    // std::vector<std::string>{"b1r", "top_pt"}, std::vector<int>{gen_ncbin_2D,
    // gen_ntop_ptbin_2D}, std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D},
    // std::vector<double *>{gen_cbins_2D, gen_top_ptbins_2D},
    // std::vector<double *>{reco_cbins_2D, reco_top_ptbins_2D}, n_sub_bins_2D);
    // residual_b1r_top_pt_binning = CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_b1r_top_pt_binning =
    // b1r_top_pt_binnings.first; reco_b1r_top_pt_binning     =
    // b1r_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> b2r_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b2r_top_pt",
    // std::vector<std::string>{"b2r", "top_pt"}, std::vector<int>{gen_ncbin_2D,
    // gen_ntop_ptbin_2D}, std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D},
    // std::vector<double *>{gen_cbins_2D, gen_top_ptbins_2D},
    // std::vector<double *>{reco_cbins_2D, reco_top_ptbins_2D}, n_sub_bins_2D);
    // residual_b2r_top_pt_binning = CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_b2r_top_pt_binning =
    // b2r_top_pt_binnings.first; reco_b2r_top_pt_binning     =
    // b2r_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> b1q_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b1q_top_pt",
    // std::vector<std::string>{"b1q", "top_pt"}, std::vector<int>{gen_ncbin_2D,
    // gen_ntop_ptbin_2D}, std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D},
    // std::vector<double *>{gen_cbins_2D, gen_top_ptbins_2D},
    // std::vector<double *>{reco_cbins_2D, reco_top_ptbins_2D}, n_sub_bins_2D);
    // residual_b1q_top_pt_binning = CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_b1q_top_pt_binning  =
    // b1q_top_pt_binnings.first; reco_b1q_top_pt_binning =
    // b1q_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> b2q_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b2q_top_pt",
    // std::vector<std::string>{"b2q", "top_pt"}, std::vector<int>{gen_ncbin_2D,
    // gen_ntop_ptbin_2D}, std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D},
    // std::vector<double *>{gen_cbins_2D, gen_top_ptbins_2D},
    // std::vector<double *>{reco_cbins_2D, reco_top_ptbins_2D}, n_sub_bins_2D);
    // residual_b2q_top_pt_binning = CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_b2q_top_pt_binning  =
    // b2q_top_pt_binnings.first; reco_b2q_top_pt_binning =
    // b2q_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> b1n_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b1n_top_pt",
    // std::vector<std::string>{"b1n", "top_pt"}, std::vector<int>{gen_ncbin_2D,
    // gen_ntop_ptbin_2D}, std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D},
    // std::vector<double *>{gen_cbins_2D, gen_top_ptbins_2D},
    // std::vector<double *>{reco_cbins_2D, reco_top_ptbins_2D}, n_sub_bins_2D);
    // residual_b1n_top_pt_binning = CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_b1n_top_pt_binning  =
    // b1n_top_pt_binnings.first; reco_b1n_top_pt_binning =
    // b1n_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> b2n_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b2n_top_pt",
    // std::vector<std::string>{"b2n", "top_pt"}, std::vector<int>{gen_ncbin_2D,
    // gen_ntop_ptbin_2D}, std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D},
    // std::vector<double *>{gen_cbins_2D, gen_top_ptbins_2D},
    // std::vector<double *>{reco_cbins_2D, reco_top_ptbins_2D}, n_sub_bins_2D);
    // residual_b2n_top_pt_binning = CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_b2n_top_pt_binning  =
    // b2n_top_pt_binnings.first; reco_b2n_top_pt_binning =
    // b2n_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_kk_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_kk_top_pt",
    // std::vector<std::string>{"c_kk", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_kk_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_kk_top_pt_binning
    // = c_kk_top_pt_binnings.first; reco_c_kk_top_pt_binning =
    // c_kk_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rr_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_rr_top_pt",
    // std::vector<std::string>{"c_rr", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_rr_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_rr_top_pt_binning
    // = c_rr_top_pt_binnings.first; reco_c_rr_top_pt_binning =
    // c_rr_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nn_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_nn_top_pt",
    // std::vector<std::string>{"c_nn", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_nn_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_nn_top_pt_binning
    // = c_nn_top_pt_binnings.first; reco_c_nn_top_pt_binning =
    // c_nn_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_kj_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_kj_top_pt",
    // std::vector<std::string>{"c_kj", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_kj_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_kj_top_pt_binning
    // = c_kj_top_pt_binnings.first; reco_c_kj_top_pt_binning =
    // c_kj_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rq_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_rq_top_pt",
    // std::vector<std::string>{"c_rq", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_rq_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_rq_top_pt_binning
    // = c_rq_top_pt_binnings.first; reco_c_rq_top_pt_binning =
    // c_rq_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Prk_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_Prk_top_pt",
    // std::vector<std::string>{"c_Prk", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_Prk_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_Prk_top_pt_binning
    // = c_Prk_top_pt_binnings.first; reco_c_Prk_top_pt_binning     =
    // c_Prk_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Mrk_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_Mrk_top_pt",
    // std::vector<std::string>{"c_Mrk", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_Mrk_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_Mrk_top_pt_binning
    // = c_Mrk_top_pt_binnings.first; reco_c_Mrk_top_pt_binning     =
    // c_Mrk_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Pnr_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_Pnr_top_pt",
    // std::vector<std::string>{"c_Pnr", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_Pnr_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_Pnr_top_pt_binning
    // = c_Pnr_top_pt_binnings.first; reco_c_Pnr_top_pt_binning     =
    // c_Pnr_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Mnr_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_Mnr_top_pt",
    // std::vector<std::string>{"c_Mnr", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_Mnr_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_Mnr_top_pt_binning
    // = c_Mnr_top_pt_binnings.first; reco_c_Mnr_top_pt_binning =
    // c_Mnr_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Pnk_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_Pnk_top_pt",
    // std::vector<std::string>{"c_Pnk", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_Pnk_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_Pnk_top_pt_binning
    // = c_Pnk_top_pt_binnings.first; reco_c_Pnk_top_pt_binning =
    // c_Pnk_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Mnk_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_Mnk_top_pt",
    // std::vector<std::string>{"c_Mnk", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_Mnk_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_Mnk_top_pt_binning
    // = c_Mnk_top_pt_binnings.first; reco_c_Mnk_top_pt_binning =
    // c_Mnk_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rk_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_rk_top_pt",
    // std::vector<std::string>{"c_rk", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_rk_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_rk_top_pt_binning
    // = c_rk_top_pt_binnings.first; reco_c_rk_top_pt_binning =
    // c_rk_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_kr_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_kr_top_pt",
    // std::vector<std::string>{"c_kr", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_kr_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_kr_top_pt_binning
    // = c_kr_top_pt_binnings.first; reco_c_kr_top_pt_binning =
    // c_kr_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nr_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_nr_top_pt",
    // std::vector<std::string>{"c_nr", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_nr_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_nr_top_pt_binning
    // = c_nr_top_pt_binnings.first; reco_c_nr_top_pt_binning =
    // c_nr_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rn_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_rn_top_pt",
    // std::vector<std::string>{"c_rn", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_rn_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_rn_top_pt_binning
    // = c_rn_top_pt_binnings.first; reco_c_rn_top_pt_binning =
    // c_rn_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nk_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_nk_top_pt",
    // std::vector<std::string>{"c_nk", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_nk_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_nk_top_pt_binning
    // = c_nk_top_pt_binnings.first; reco_c_nk_top_pt_binning =
    // c_nk_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_kn_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_kn_top_pt",
    // std::vector<std::string>{"c_kn", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_kn_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_kn_top_pt_binning
    // = c_kn_top_pt_binnings.first; reco_c_kn_top_pt_binning =
    // c_kn_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_han_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_han_top_pt",
    // std::vector<std::string>{"c_han", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_han_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_han_top_pt_binning
    // = c_han_top_pt_binnings.first; reco_c_han_top_pt_binning =
    // c_han_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_sca_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_sca_top_pt",
    // std::vector<std::string>{"c_sca", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_sca_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_sca_top_pt_binning
    // = c_sca_top_pt_binnings.first; reco_c_sca_top_pt_binning =
    // c_sca_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_tra_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_tra_top_pt",
    // std::vector<std::string>{"c_tra", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_tra_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_tra_top_pt_binning
    // = c_tra_top_pt_binnings.first; reco_c_tra_top_pt_binning =
    // c_tra_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_kjL_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_kjL_top_pt",
    // std::vector<std::string>{"c_kjL", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_kjL_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_kjL_top_pt_binning
    // = c_kjL_top_pt_binnings.first; reco_c_kjL_top_pt_binning =
    // c_kjL_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rqL_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_rqL_top_pt",
    // std::vector<std::string>{"c_rqL", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_rqL_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_rqL_top_pt_binning
    // = c_rqL_top_pt_binnings.first; reco_c_rqL_top_pt_binning =
    // c_rqL_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rkP_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_rkP_top_pt",
    // std::vector<std::string>{"c_rkP", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_rkP_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_rkP_top_pt_binning
    // = c_rkP_top_pt_binnings.first; reco_c_rkP_top_pt_binning     =
    // c_rkP_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rkM_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_rkM_top_pt",
    // std::vector<std::string>{"c_rkM", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_rkM_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_rkM_top_pt_binning
    // = c_rkM_top_pt_binnings.first; reco_c_rkM_top_pt_binning     =
    // c_rkM_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nrP_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_nrP_top_pt",
    // std::vector<std::string>{"c_nrP", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_nrP_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_nrP_top_pt_binning
    // = c_nrP_top_pt_binnings.first; reco_c_nrP_top_pt_binning     =
    // c_nrP_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nrM_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_nrM_top_pt",
    // std::vector<std::string>{"c_nrM", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_nrM_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_nrM_top_pt_binning
    // = c_nrM_top_pt_binnings.first; reco_c_nrM_top_pt_binning     =
    // c_nrM_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nkP_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_nkP_top_pt",
    // std::vector<std::string>{"c_nkP", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_nkP_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_nkP_top_pt_binning
    // = c_nkP_top_pt_binnings.first; reco_c_nkP_top_pt_binning     =
    // c_nkP_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_nkM_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_nkM_top_pt",
    // std::vector<std::string>{"c_nkM", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_nkM_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_nkM_top_pt_binning
    // = c_nkM_top_pt_binnings.first; reco_c_nkM_top_pt_binning     =
    // c_nkM_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> ll_cHel_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("ll_cHel_top_pt",
    // std::vector<std::string>{"ll_cHel", "top_pt"},
    // std::vector<int>{gen_ncosbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncosbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cosbins_2D, gen_top_ptbins_2D}, std::vector<double
    // *>{reco_cosbins_2D, reco_top_ptbins_2D}, n_sub_bins_2D);
    // residual_ll_cHel_top_pt_binning = CreateResidualBinning(nbinsres,
    // reco_cosbins_2D[0]-reco_cosbins_2D[reco_ncosbin_2D],
    // reco_cosbins_2D[reco_ncosbin_2D]-reco_cosbins_2D[0]);
    // gen_ll_cHel_top_pt_binning      = ll_cHel_top_pt_binnings.first;
    // reco_ll_cHel_top_pt_binning     = ll_cHel_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> ll_cLab_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("ll_cLab_top_pt",
    // std::vector<std::string>{"ll_cLab", "top_pt"},
    // std::vector<int>{gen_ncosbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncosbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cosbins_2D, gen_top_ptbins_2D}, std::vector<double
    // *>{reco_cosbins_2D, reco_top_ptbins_2D}, n_sub_bins_2D);
    // residual_ll_cLab_top_pt_binning = CreateResidualBinning(nbinsres,
    // reco_cosbins_2D[0]-reco_cosbins_2D[reco_ncosbin_2D],
    // reco_cosbins_2D[reco_ncosbin_2D]-reco_cosbins_2D[0]);
    // gen_ll_cLab_top_pt_binning      = ll_cLab_top_pt_binnings.first;
    // reco_ll_cLab_top_pt_binning     = ll_cLab_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> ll_kNorm_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("ll_kNorm_top_pt",
    // std::vector<std::string>{"ll_kNorm", "top_pt"},
    // std::vector<int>{gen_ncosbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncosbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cosbins_2D, gen_top_ptbins_2D}, std::vector<double
    // *>{reco_cosbins_2D, reco_top_ptbins_2D}, n_sub_bins_2D);
    // residual_ll_kNorm_top_pt_binning = CreateResidualBinning(nbinsres,
    // reco_cosbins_2D[0]-reco_cosbins_2D[reco_ncosbin_2D],
    // reco_cosbins_2D[reco_ncosbin_2D]-reco_cosbins_2D[0]);
    // gen_ll_kNorm_top_pt_binning      = ll_kNorm_top_pt_binnings.first;
    // reco_ll_kNorm_top_pt_binning     = ll_kNorm_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> ll_rNorm_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("ll_rNorm_top_pt",
    // std::vector<std::string>{"ll_rNorm", "top_pt"},
    // std::vector<int>{gen_ncosbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncosbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cosbins_2D, gen_top_ptbins_2D}, std::vector<double
    // *>{reco_cosbins_2D, reco_top_ptbins_2D}, n_sub_bins_2D);
    // residual_ll_rNorm_top_pt_binning = CreateResidualBinning(nbinsres,
    // reco_cosbins_2D[0]-reco_cosbins_2D[reco_ncosbin_2D],
    // reco_cosbins_2D[reco_ncosbin_2D]-reco_cosbins_2D[0]);
    // gen_ll_rNorm_top_pt_binning      = ll_rNorm_top_pt_binnings.first;
    // reco_ll_rNorm_top_pt_binning     = ll_rNorm_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // llbar_delta_phi_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("llbar_delta_phi_top_pt",
    // std::vector<std::string>{"llbar_delta_phi", "top_pt"},
    // std::vector<int>{gen_nllbar_delta_phibin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_nllbar_delta_phibin_2D, reco_ntop_ptbin_2D},
    // std::vector<double *>{gen_llbar_delta_phibins_2D, gen_top_ptbins_2D},
    // std::vector<double *>{reco_llbar_delta_phibins_2D, reco_top_ptbins_2D},
    // n_sub_bins_2D); residual_llbar_delta_phi_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_llbar_delta_phibins_2D[0]-reco_llbar_delta_phibins_2D[reco_nllbar_delta_phibin_2D],
    // reco_llbar_delta_phibins_2D[reco_nllbar_delta_phibin_2D]-reco_llbar_delta_phibins_2D[0]);
    // gen_llbar_delta_phi_top_pt_binning      =
    // llbar_delta_phi_top_pt_binnings.first;
    // reco_llbar_delta_phi_top_pt_binning     =
    // llbar_delta_phi_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // llbar_delta_eta_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("llbar_delta_eta_top_pt",
    // std::vector<std::string>{"llbar_delta_eta", "top_pt"},
    // std::vector<int>{gen_nllbar_delta_etabin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_nllbar_delta_etabin_2D, reco_ntop_ptbin_2D},
    // std::vector<double *>{gen_llbar_delta_etabins_2D, gen_top_ptbins_2D},
    // std::vector<double *>{reco_llbar_delta_etabins_2D, reco_top_ptbins_2D},
    // n_sub_bins_2D); residual_llbar_delta_eta_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_llbar_delta_etabins_2D[0]-reco_llbar_delta_etabins_2D[reco_nllbar_delta_etabin_2D],
    // reco_llbar_delta_etabins_2D[reco_nllbar_delta_etabin_2D]-reco_llbar_delta_etabins_2D[0]);
    // gen_llbar_delta_eta_top_pt_binning      =
    // llbar_delta_eta_top_pt_binnings.first;
    // reco_llbar_delta_eta_top_pt_binning     =
    // llbar_delta_eta_top_pt_binnings.second;

    // Amandeep : adding starred variable binnings here
    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_rj_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_rj_top_pt",
    // std::vector<std::string>{"c_rj", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_rj_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_rj_top_pt_binning
    // = c_rj_top_pt_binnings.first; reco_c_rj_top_pt_binning =
    // c_rj_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_jr_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_jr_top_pt",
    // std::vector<std::string>{"c_jr", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_jr_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_jr_top_pt_binning
    // = c_jr_top_pt_binnings.first; reco_c_jr_top_pt_binning =
    // c_jr_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Prj_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_Prj_top_pt",
    // std::vector<std::string>{"c_Prj", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_Prj_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_Prj_top_pt_binning
    // = c_Prj_top_pt_binnings.first; reco_c_Prj_top_pt_binning     =
    // c_Prj_top_pt_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *> c_Mrj_top_pt_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_Mrj_top_pt",
    // std::vector<std::string>{"c_Mrj", "top_pt"},
    // std::vector<int>{gen_ncbin_2D, gen_ntop_ptbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ntop_ptbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_top_ptbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_top_ptbins_2D}, n_sub_bins_2D); residual_c_Mrj_top_pt_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]); gen_c_Mrj_top_pt_binning
    // = c_Mrj_top_pt_binnings.first; reco_c_Mrj_top_pt_binning     =
    // c_Mrj_top_pt_binnings.second;

    // ********************
    // top_scatteringangle_ttbarframe
    // ********************

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // b1k_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b1k_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"b1k", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_b1k_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_b1k_top_scatteringangle_ttbarframe_binning      =
    // b1k_top_scatteringangle_ttbarframe_binnings.first;
    // reco_b1k_top_scatteringangle_ttbarframe_binning     =
    // b1k_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // b2k_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b2k_top_pt",
    // std::vector<std::string>{"b2k ", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_b2k_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_b2k_top_scatteringangle_ttbarframe_binning      =
    // b2k_top_scatteringangle_ttbarframe_binnings.first;
    // reco_b2k_top_scatteringangle_ttbarframe_binning     =
    // b2k_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // b1j_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b1j_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"b1j", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_b1j_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_b1j_top_scatteringangle_ttbarframe_binning      =
    // b1j_top_scatteringangle_ttbarframe_binnings.first;
    // reco_b1j_top_scatteringangle_ttbarframe_binning     =
    // b1j_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // b2j_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b2j_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"b2j", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_b2j_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_b2j_top_scatteringangle_ttbarframe_binning      =
    // b2j_top_scatteringangle_ttbarframe_binnings.first;
    // reco_b2j_top_scatteringangle_ttbarframe_binning     =
    // b2j_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // b1r_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b1r_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"b1r", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_b1r_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_b1r_top_scatteringangle_ttbarframe_binning      =
    // b1r_top_scatteringangle_ttbarframe_binnings.first;
    // reco_b1r_top_scatteringangle_ttbarframe_binning     =
    // b1r_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // b2r_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b2r_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"b2r", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_b2r_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_b2r_top_scatteringangle_ttbarframe_binning      =
    // b2r_top_scatteringangle_ttbarframe_binnings.first;
    // reco_b2r_top_scatteringangle_ttbarframe_binning     =
    // b2r_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // b1q_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b1q_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"b1q", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_b1q_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_b1q_top_scatteringangle_ttbarframe_binning  =
    // b1q_top_scatteringangle_ttbarframe_binnings.first;
    // reco_b1q_top_scatteringangle_ttbarframe_binning =
    // b1q_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // b2q_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b2q_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"b2q", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_b2q_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_b2q_top_scatteringangle_ttbarframe_binning  =
    // b2q_top_scatteringangle_ttbarframe_binnings.first;
    // reco_b2q_top_scatteringangle_ttbarframe_binning =
    // b2q_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // b1n_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b1n_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"b1n", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_b1n_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_b1n_top_scatteringangle_ttbarframe_binning  =
    // b1n_top_scatteringangle_ttbarframe_binnings.first;
    // reco_b1n_top_scatteringangle_ttbarframe_binning =
    // b1n_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // b2n_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("b2n_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"b2n", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_b2n_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_b2n_top_scatteringangle_ttbarframe_binning  =
    // b2n_top_scatteringangle_ttbarframe_binnings.first;
    // reco_b2n_top_scatteringangle_ttbarframe_binning =
    // b2n_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_kk_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_kk_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_kk", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_kk_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_kk_top_scatteringangle_ttbarframe_binning  =
    // c_kk_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_kk_top_scatteringangle_ttbarframe_binning =
    // c_kk_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_rr_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_rr_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_rr", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_rr_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_rr_top_scatteringangle_ttbarframe_binning  =
    // c_rr_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_rr_top_scatteringangle_ttbarframe_binning =
    // c_rr_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_nn_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_nn_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_nn", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_nn_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_nn_top_scatteringangle_ttbarframe_binning  =
    // c_nn_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_nn_top_scatteringangle_ttbarframe_binning =
    // c_nn_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_kj_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_kj_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_kj", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_kj_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_kj_top_scatteringangle_ttbarframe_binning  =
    // c_kj_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_kj_top_scatteringangle_ttbarframe_binning =
    // c_kj_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_rq_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_rq_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_rq", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_rq_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_rq_top_scatteringangle_ttbarframe_binning  =
    // c_rq_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_rq_top_scatteringangle_ttbarframe_binning =
    // c_rq_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_Prk_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_Prk_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_Prk", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_Prk_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_Prk_top_scatteringangle_ttbarframe_binning      =
    // c_Prk_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_Prk_top_scatteringangle_ttbarframe_binning     =
    // c_Prk_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_Mrk_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_Mrk_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_Mrk", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_Mrk_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_Mrk_top_scatteringangle_ttbarframe_binning      =
    // c_Mrk_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_Mrk_top_scatteringangle_ttbarframe_binning     =
    // c_Mrk_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_Pnr_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_Pnr_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_Pnr", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_Pnr_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_Pnr_top_scatteringangle_ttbarframe_binning      =
    // c_Pnr_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_Pnr_top_scatteringangle_ttbarframe_binning     =
    // c_Pnr_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_Mnr_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_Mnr_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_Mnr", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_Mnr_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_Mnr_top_scatteringangle_ttbarframe_binning  =
    // c_Mnr_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_Mnr_top_scatteringangle_ttbarframe_binning =
    // c_Mnr_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_Pnk_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_Pnk_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_Pnk", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_Pnk_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_Pnk_top_scatteringangle_ttbarframe_binning  =
    // c_Pnk_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_Pnk_top_scatteringangle_ttbarframe_binning =
    // c_Pnk_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_Mnk_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_Mnk_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_Mnk", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_Mnk_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_Mnk_top_scatteringangle_ttbarframe_binning  =
    // c_Mnk_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_Mnk_top_scatteringangle_ttbarframe_binning =
    // c_Mnk_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_rk_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_rk_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_rk", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_rk_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_rk_top_scatteringangle_ttbarframe_binning  =
    // c_rk_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_rk_top_scatteringangle_ttbarframe_binning =
    // c_rk_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_kr_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_kr_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_kr", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_kr_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_kr_top_scatteringangle_ttbarframe_binning  =
    // c_kr_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_kr_top_scatteringangle_ttbarframe_binning =
    // c_kr_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_nr_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_nr_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_nr", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_nr_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_nr_top_scatteringangle_ttbarframe_binning  =
    // c_nr_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_nr_top_scatteringangle_ttbarframe_binning =
    // c_nr_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_rn_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_rn_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_rn", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_rn_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_rn_top_scatteringangle_ttbarframe_binning  =
    // c_rn_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_rn_top_scatteringangle_ttbarframe_binning =
    // c_rn_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_nk_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_nk_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_nk", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_nk_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_nk_top_scatteringangle_ttbarframe_binning  =
    // c_nk_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_nk_top_scatteringangle_ttbarframe_binning =
    // c_nk_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_kn_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_kn_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_kn", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_kn_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_kn_top_scatteringangle_ttbarframe_binning  =
    // c_kn_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_kn_top_scatteringangle_ttbarframe_binning =
    // c_kn_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_han_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_han_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_han", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_han_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_han_top_scatteringangle_ttbarframe_binning  =
    // c_han_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_han_top_scatteringangle_ttbarframe_binning =
    // c_han_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_sca_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_sca_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_sca", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_sca_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_sca_top_scatteringangle_ttbarframe_binning  =
    // c_sca_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_sca_top_scatteringangle_ttbarframe_binning =
    // c_sca_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_tra_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_tra_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_tra", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_tra_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_tra_top_scatteringangle_ttbarframe_binning  =
    // c_tra_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_tra_top_scatteringangle_ttbarframe_binning =
    // c_tra_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_kjL_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_kjL_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_kjL", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_kjL_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_kjL_top_scatteringangle_ttbarframe_binning  =
    // c_kjL_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_kjL_top_scatteringangle_ttbarframe_binning =
    // c_kjL_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_rqL_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_rqL_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_rqL", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_rqL_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_rqL_top_scatteringangle_ttbarframe_binning  =
    // c_rqL_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_rqL_top_scatteringangle_ttbarframe_binning =
    // c_rqL_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_rkP_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_rkP_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_rkP", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_rkP_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_rkP_top_scatteringangle_ttbarframe_binning      =
    // c_rkP_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_rkP_top_scatteringangle_ttbarframe_binning     =
    // c_rkP_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_rkM_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_rkM_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_rkM", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_rkM_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_rkM_top_scatteringangle_ttbarframe_binning      =
    // c_rkM_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_rkM_top_scatteringangle_ttbarframe_binning     =
    // c_rkM_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_nrP_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_nrP_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_nrP", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_nrP_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_nrP_top_scatteringangle_ttbarframe_binning      =
    // c_nrP_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_nrP_top_scatteringangle_ttbarframe_binning     =
    // c_nrP_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_nrM_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_nrM_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_nrM", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_nrM_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_nrM_top_scatteringangle_ttbarframe_binning      =
    // c_nrM_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_nrM_top_scatteringangle_ttbarframe_binning     =
    // c_nrM_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_nkP_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_nkP_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_nkP", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_nkP_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_nkP_top_scatteringangle_ttbarframe_binning      =
    // c_nkP_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_nkP_top_scatteringangle_ttbarframe_binning     =
    // c_nkP_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_nkM_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_nkM_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_nkM", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_nkM_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_nkM_top_scatteringangle_ttbarframe_binning      =
    // c_nkM_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_nkM_top_scatteringangle_ttbarframe_binning     =
    // c_nkM_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // ll_cHel_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("ll_cHel_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"ll_cHel", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncosbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncosbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cosbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cosbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_ll_cHel_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cosbins_2D[0]-reco_cosbins_2D[reco_ncosbin_2D],
    // reco_cosbins_2D[reco_ncosbin_2D]-reco_cosbins_2D[0]);
    // gen_ll_cHel_top_scatteringangle_ttbarframe_binning      =
    // ll_cHel_top_scatteringangle_ttbarframe_binnings.first;
    // reco_ll_cHel_top_scatteringangle_ttbarframe_binning     =
    // ll_cHel_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // ll_cLab_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("ll_cLab_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"ll_cLab", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncosbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncosbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cosbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cosbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_ll_cLab_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cosbins_2D[0]-reco_cosbins_2D[reco_ncosbin_2D],
    // reco_cosbins_2D[reco_ncosbin_2D]-reco_cosbins_2D[0]);
    // gen_ll_cLab_top_scatteringangle_ttbarframe_binning      =
    // ll_cLab_top_scatteringangle_ttbarframe_binnings.first;
    // reco_ll_cLab_top_scatteringangle_ttbarframe_binning     =
    // ll_cLab_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // ll_kNorm_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("ll_kNorm_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"ll_kNorm", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncosbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncosbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cosbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cosbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_ll_kNorm_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cosbins_2D[0]-reco_cosbins_2D[reco_ncosbin_2D],
    // reco_cosbins_2D[reco_ncosbin_2D]-reco_cosbins_2D[0]);
    // gen_ll_kNorm_top_scatteringangle_ttbarframe_binning      =
    // ll_kNorm_top_scatteringangle_ttbarframe_binnings.first;
    // reco_ll_kNorm_top_scatteringangle_ttbarframe_binning     =
    // ll_kNorm_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // ll_rNorm_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("ll_rNorm_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"ll_rNorm", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncosbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncosbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cosbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cosbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_ll_rNorm_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cosbins_2D[0]-reco_cosbins_2D[reco_ncosbin_2D],
    // reco_cosbins_2D[reco_ncosbin_2D]-reco_cosbins_2D[0]);
    // gen_ll_rNorm_top_scatteringangle_ttbarframe_binning      =
    // ll_rNorm_top_scatteringangle_ttbarframe_binnings.first;
    // reco_ll_rNorm_top_scatteringangle_ttbarframe_binning     =
    // ll_rNorm_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // llbar_delta_phi_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("llbar_delta_phi_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"llbar_delta_phi",
    // "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_nllbar_delta_phibin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_nllbar_delta_phibin_2D, reco_ncbin_2D},
    // std::vector<double *>{gen_llbar_delta_phibins_2D, gen_cbins_2D},
    // std::vector<double *>{reco_llbar_delta_phibins_2D, reco_cbins_2D},
    // n_sub_bins_2D);
    // residual_llbar_delta_phi_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_llbar_delta_phibins_2D[0]-reco_llbar_delta_phibins_2D[reco_nllbar_delta_phibin_2D],
    // reco_llbar_delta_phibins_2D[reco_nllbar_delta_phibin_2D]-reco_llbar_delta_phibins_2D[0]);
    // gen_llbar_delta_phi_top_scatteringangle_ttbarframe_binning      =
    // llbar_delta_phi_top_scatteringangle_ttbarframe_binnings.first;
    // reco_llbar_delta_phi_top_scatteringangle_ttbarframe_binning     =
    // llbar_delta_phi_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // llbar_delta_eta_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("llbar_delta_eta_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"llbar_delta_eta",
    // "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_nllbar_delta_etabin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_nllbar_delta_etabin_2D, reco_ncbin_2D},
    // std::vector<double *>{gen_llbar_delta_etabins_2D, gen_cbins_2D},
    // std::vector<double *>{reco_llbar_delta_etabins_2D, reco_cbins_2D},
    // n_sub_bins_2D);
    // residual_llbar_delta_eta_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_llbar_delta_etabins_2D[0]-reco_llbar_delta_etabins_2D[reco_nllbar_delta_etabin_2D],
    // reco_llbar_delta_etabins_2D[reco_nllbar_delta_etabin_2D]-reco_llbar_delta_etabins_2D[0]);
    // gen_llbar_delta_eta_top_scatteringangle_ttbarframe_binning      =
    // llbar_delta_eta_top_scatteringangle_ttbarframe_binnings.first;
    // reco_llbar_delta_eta_top_scatteringangle_ttbarframe_binning     =
    // llbar_delta_eta_top_scatteringangle_ttbarframe_binnings.second;

    // // Amandeep : adding starred variable binnings here
    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_rj_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_rj_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_rj", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_rj_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_rj_top_scatteringangle_ttbarframe_binning  =
    // c_rj_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_rj_top_scatteringangle_ttbarframe_binning =
    // c_rj_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_jr_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_jr_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_jr", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_jr_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_jr_top_scatteringangle_ttbarframe_binning  =
    // c_jr_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_jr_top_scatteringangle_ttbarframe_binning =
    // c_jr_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_Prj_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_Prj_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_Prj", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_Prj_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_Prj_top_scatteringangle_ttbarframe_binning      =
    // c_Prj_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_Prj_top_scatteringangle_ttbarframe_binning     =
    // c_Prj_top_scatteringangle_ttbarframe_binnings.second;

    // std::pair<TUnfoldBinning *, TUnfoldBinning *>
    // c_Mrj_top_scatteringangle_ttbarframe_binnings =
    // makeTUnfoldHisto::CreateXMLBinFiles("c_Mrj_top_scatteringangle_ttbarframe",
    // std::vector<std::string>{"c_Mrj", "top_scatteringangle_ttbarframe"},
    // std::vector<int>{gen_ncbin_2D, gen_ncbin_2D},
    // std::vector<int>{reco_ncbin_2D, reco_ncbin_2D}, std::vector<double
    // *>{gen_cbins_2D, gen_cbins_2D}, std::vector<double *>{reco_cbins_2D,
    // reco_cbins_2D}, n_sub_bins_2D);
    // residual_c_Mrj_top_scatteringangle_ttbarframe_binning =
    // CreateResidualBinning(nbinsres,
    // reco_cbins_2D[0]-reco_cbins_2D[reco_ncbin_2D],
    // reco_cbins_2D[reco_ncbin_2D]-reco_cbins_2D[0]);
    // gen_c_Mrj_top_scatteringangle_ttbarframe_binning      =
    // c_Mrj_top_scatteringangle_ttbarframe_binnings.first;
    // reco_c_Mrj_top_scatteringangle_ttbarframe_binning     =
    // c_Mrj_mttbar_binnings.second;

    // **************************
    // 1D reco histogram creation
    // **************************

    // Reconstructed info -- filled for both data and MC
    // Kinematic
    hreco_top_pt = reco_top_pt_binning->CreateHistogram("hreco_top_pt");
    hreco_l_pt = reco_l_pt_binning->CreateHistogram("hreco_l_pt");
    hreco_lbar_pt = reco_lbar_pt_binning->CreateHistogram("hreco_lbar_pt");
    hreco_ttbar_pt = reco_ttbar_pt_binning->CreateHistogram("hreco_ttbar_pt");
    hreco_ttbar_mass = reco_ttbar_mass_binning->CreateHistogram("hreco_ttbar_mass");
    hreco_ttbar_delta_phi = reco_ttbar_delta_phi_binning->CreateHistogram("hreco_ttbar_delta_phi");
    hreco_ttbar_delta_eta = reco_ttbar_delta_eta_binning->CreateHistogram("hreco_ttbar_delta_eta");
    hreco_ttbar_rapidity = reco_ttbar_rapidity_binning->CreateHistogram("hreco_ttbar_rapidity");
    hreco_llbar_pt = reco_llbar_pt_binning->CreateHistogram("hreco_llbar_pt");
    hreco_llbar_mass = reco_llbar_mass_binning->CreateHistogram("hreco_llbar_mass");

    // Spin corr
    hreco_b1k = reco_b1k_binning->CreateHistogram("hreco_b1k");
    hreco_b2k = reco_b2k_binning->CreateHistogram("hreco_b2k");
    hreco_b1j = reco_b1j_binning->CreateHistogram("hreco_b1j");
    hreco_b2j = reco_b2j_binning->CreateHistogram("hreco_b2j");
    hreco_b1r = reco_b1r_binning->CreateHistogram("hreco_b1r");
    hreco_b2r = reco_b2r_binning->CreateHistogram("hreco_b2r");
    hreco_b1q = reco_b1q_binning->CreateHistogram("hreco_b1q");
    hreco_b2q = reco_b2q_binning->CreateHistogram("hreco_b2q");
    hreco_b1n = reco_b1n_binning->CreateHistogram("hreco_b1n");
    hreco_b2n = reco_b2n_binning->CreateHistogram("hreco_b2n");

    hreco_c_kk = reco_c_kk_binning->CreateHistogram("hreco_c_kk");
    hreco_c_rr = reco_c_rr_binning->CreateHistogram("hreco_c_rr");
    hreco_c_nn = reco_c_nn_binning->CreateHistogram("hreco_c_nn");

    hreco_c_rk = reco_c_rk_binning->CreateHistogram("hreco_c_rk");
    hreco_c_kr = reco_c_kr_binning->CreateHistogram("hreco_c_kr");
    hreco_c_nr = reco_c_nr_binning->CreateHistogram("hreco_c_nr");
    hreco_c_rn = reco_c_rn_binning->CreateHistogram("hreco_c_rn");
    hreco_c_nk = reco_c_nk_binning->CreateHistogram("hreco_c_nk");
    hreco_c_kn = reco_c_kn_binning->CreateHistogram("hreco_c_kn");
    hreco_c_Prk = reco_c_Prk_binning->CreateHistogram("hreco_c_Prk");
    hreco_c_Mrk = reco_c_Mrk_binning->CreateHistogram("hreco_c_Mrk");
    hreco_c_Pnr = reco_c_Pnr_binning->CreateHistogram("hreco_c_Pnr");
    hreco_c_Mnr = reco_c_Mnr_binning->CreateHistogram("hreco_c_Mnr");
    hreco_c_Pnk = reco_c_Pnk_binning->CreateHistogram("hreco_c_Pnk");
    hreco_c_Mnk = reco_c_Mnk_binning->CreateHistogram("hreco_c_Mnk");

    hreco_c_kj = reco_c_kj_binning->CreateHistogram("hreco_c_kj");
    hreco_c_rq = reco_c_rq_binning->CreateHistogram("hreco_c_rq");
    hreco_c_rj = reco_c_rj_binning->CreateHistogram("hreco_c_rj");
    hreco_c_jr = reco_c_jr_binning->CreateHistogram("hreco_c_jr");

    hreco_c_Prj = reco_c_Prj_binning->CreateHistogram("hreco_c_Prj");
    hreco_c_Mrj = reco_c_Mrj_binning->CreateHistogram("hreco_c_Mrj");

    hreco_c_han = reco_c_han_binning->CreateHistogram("hreco_c_han");
    hreco_c_sca = reco_c_sca_binning->CreateHistogram("hreco_c_sca");
    hreco_c_tra = reco_c_tra_binning->CreateHistogram("hreco_c_tra");
    hreco_c_kjL = reco_c_kjL_binning->CreateHistogram("hreco_c_kjL");
    hreco_c_rqL = reco_c_rqL_binning->CreateHistogram("hreco_c_rqL");
    hreco_c_rkP = reco_c_rkP_binning->CreateHistogram("hreco_c_rkP");
    hreco_c_rkM = reco_c_rkM_binning->CreateHistogram("hreco_c_rkM");
    hreco_c_nrP = reco_c_nrP_binning->CreateHistogram("hreco_c_nrP");
    hreco_c_nrM = reco_c_nrM_binning->CreateHistogram("hreco_c_nrM");
    hreco_c_nkP = reco_c_nkP_binning->CreateHistogram("hreco_c_nkP");
    hreco_c_nkM = reco_c_nkM_binning->CreateHistogram("hreco_c_nkM");

    hreco_ll_cHel = reco_ll_cHel_binning->CreateHistogram("hreco_ll_cHel");
    hreco_ll_cLab = reco_ll_cLab_binning->CreateHistogram("hreco_ll_cLab");
    hreco_ll_kNorm = reco_ll_kNorm_binning->CreateHistogram("hreco_ll_kNorm");
    hreco_ll_rNorm = reco_ll_rNorm_binning->CreateHistogram("hreco_ll_rNorm");
    hreco_llbar_delta_phi = reco_llbar_delta_phi_binning->CreateHistogram("hreco_llbar_delta_phi");
    hreco_llbar_delta_eta = reco_llbar_delta_phi_binning->CreateHistogram("hreco_llbar_delta_eta");

    // **************************
    // 2D reco histogram creation
    // **************************

    // Amandeep : 2D reco histogram creation
    // hreco_c_kk_mttbar =
    // reco_c_kk_mttbar_binning->CreateHistogram("hreco_c_kk_mttbar");

    // ******
    // mttbar
    // ******

    hreco_b1k_exjet = reco_b1k_exjet_binning->CreateHistogram("hreco_b1k_exjet");
    hreco_b2k_exjet = reco_b2k_exjet_binning->CreateHistogram("hreco_b2k_exjet");
    hreco_b1j_exjet = reco_b1j_exjet_binning->CreateHistogram("hreco_b1j_exjet");
    hreco_b2j_exjet = reco_b2j_exjet_binning->CreateHistogram("hreco_b2j_exjet");
    hreco_b1r_exjet = reco_b1r_exjet_binning->CreateHistogram("hreco_b1r_exjet");
    hreco_b2r_exjet = reco_b2r_exjet_binning->CreateHistogram("hreco_b2r_exjet");
    hreco_b1q_exjet = reco_b1q_exjet_binning->CreateHistogram("hreco_b1q_exjet");
    hreco_b2q_exjet = reco_b2q_exjet_binning->CreateHistogram("hreco_b2q_exjet");
    hreco_b1n_exjet = reco_b1n_exjet_binning->CreateHistogram("hreco_b1n_exjet");
    hreco_b2n_exjet = reco_b2n_exjet_binning->CreateHistogram("hreco_b2n_exjet");

    hreco_c_kk_exjet = reco_c_kk_exjet_binning->CreateHistogram("hreco_c_kk_exjet");
    hreco_c_rr_exjet = reco_c_rr_exjet_binning->CreateHistogram("hreco_c_rr_exjet");
    hreco_c_nn_exjet = reco_c_nn_exjet_binning->CreateHistogram("hreco_c_nn_exjet");

    hreco_c_rk_exjet = reco_c_rk_exjet_binning->CreateHistogram("hreco_c_rk_exjet");
    hreco_c_kr_exjet = reco_c_kr_exjet_binning->CreateHistogram("hreco_c_kr_exjet");
    hreco_c_nr_exjet = reco_c_nr_exjet_binning->CreateHistogram("hreco_c_nr_exjet");
    hreco_c_rn_exjet = reco_c_rn_exjet_binning->CreateHistogram("hreco_c_rn_exjet");
    hreco_c_nk_exjet = reco_c_nk_exjet_binning->CreateHistogram("hreco_c_nk_exjet");
    hreco_c_kn_exjet = reco_c_kn_exjet_binning->CreateHistogram("hreco_c_kn_exjet");

    hreco_c_Prk_exjet = reco_c_Prk_exjet_binning->CreateHistogram("hreco_c_Prk_exjet");
    hreco_c_Mrk_exjet = reco_c_Mrk_exjet_binning->CreateHistogram("hreco_c_Mrk_exjet");
    hreco_c_Pnr_exjet = reco_c_Pnr_exjet_binning->CreateHistogram("hreco_c_Pnr_exjet");
    hreco_c_Mnr_exjet = reco_c_Mnr_exjet_binning->CreateHistogram("hreco_c_Mnr_exjet");
    hreco_c_Pnk_exjet = reco_c_Pnk_exjet_binning->CreateHistogram("hreco_c_Pnk_exjet");
    hreco_c_Mnk_exjet = reco_c_Mnk_exjet_binning->CreateHistogram("hreco_c_Mnk_exjet");

    hreco_c_kj_exjet = reco_c_kj_exjet_binning->CreateHistogram("hreco_c_kj_exjet");
    hreco_c_rq_exjet = reco_c_rq_exjet_binning->CreateHistogram("hreco_c_rq_exjet");
    hreco_c_rj_exjet = reco_c_rj_exjet_binning->CreateHistogram("hreco_c_rj_exjet");
    hreco_c_jr_exjet = reco_c_jr_exjet_binning->CreateHistogram("hreco_c_jr_exjet");

    hreco_c_Prj_exjet = reco_c_Prj_exjet_binning->CreateHistogram("hreco_c_Prj_exjet");
    hreco_c_Mrj_exjet = reco_c_Mrj_exjet_binning->CreateHistogram("hreco_c_Mrj_exjet");

    hreco_c_han_exjet = reco_c_han_exjet_binning->CreateHistogram("hreco_c_han_exjet");
    hreco_c_sca_exjet = reco_c_sca_exjet_binning->CreateHistogram("hreco_c_sca_exjet");
    hreco_c_tra_exjet = reco_c_tra_exjet_binning->CreateHistogram("hreco_c_tra_exjet");
    hreco_c_kjL_exjet = reco_c_kjL_exjet_binning->CreateHistogram("hreco_c_kjL_exjet");
    hreco_c_rqL_exjet = reco_c_rqL_exjet_binning->CreateHistogram("hreco_c_rqL_exjet");
    hreco_c_rkP_exjet = reco_c_rkP_exjet_binning->CreateHistogram("hreco_c_rkP_exjet");
    hreco_c_rkM_exjet = reco_c_rkM_exjet_binning->CreateHistogram("hreco_c_rkM_exjet");
    hreco_c_nrP_exjet = reco_c_nrP_exjet_binning->CreateHistogram("hreco_c_nrP_exjet");
    hreco_c_nrM_exjet = reco_c_nrM_exjet_binning->CreateHistogram("hreco_c_nrM_exjet");
    hreco_c_nkP_exjet = reco_c_nkP_exjet_binning->CreateHistogram("hreco_c_nkP_exjet");
    hreco_c_nkM_exjet = reco_c_nkM_exjet_binning->CreateHistogram("hreco_c_nkM_exjet");

    hreco_ll_cHel_exjet = reco_ll_cHel_exjet_binning->CreateHistogram("hreco_ll_cHel_exjet");
    hreco_ll_cLab_exjet = reco_ll_cLab_exjet_binning->CreateHistogram("hreco_ll_cLab_exjet");
    hreco_ll_kNorm_exjet = reco_ll_kNorm_exjet_binning->CreateHistogram("hreco_ll_kNorm_exjet");
    hreco_ll_rNorm_exjet = reco_ll_rNorm_exjet_binning->CreateHistogram("hreco_ll_rNorm_exjet");
    hreco_llbar_delta_phi_exjet = reco_llbar_delta_phi_exjet_binning->CreateHistogram("hreco_llbar_delta_phi_exjet");
    hreco_llbar_delta_eta_exjet = reco_llbar_delta_phi_exjet_binning->CreateHistogram("hreco_llbar_delta_eta_exjet");

    // ******
    // top_pt
    // ******

    // hreco_b1k_top_pt =
    // reco_b1k_top_pt_binning->CreateHistogram("hreco_b1k_top_pt");
    // hreco_b2k_top_pt =
    // reco_b2k_top_pt_binning->CreateHistogram("hreco_b2k_top_pt");
    // hreco_b1j_top_pt =
    // reco_b1j_top_pt_binning->CreateHistogram("hreco_b1j_top_pt");
    // hreco_b2j_top_pt =
    // reco_b2j_top_pt_binning->CreateHistogram("hreco_b2j_top_pt");
    // hreco_b1r_top_pt =
    // reco_b1r_top_pt_binning->CreateHistogram("hreco_b1r_top_pt");
    // hreco_b2r_top_pt =
    // reco_b2r_top_pt_binning->CreateHistogram("hreco_b2r_top_pt");
    // hreco_b1q_top_pt =
    // reco_b1q_top_pt_binning->CreateHistogram("hreco_b1q_top_pt");
    // hreco_b2q_top_pt =
    // reco_b2q_top_pt_binning->CreateHistogram("hreco_b2q_top_pt");
    // hreco_b1n_top_pt =
    // reco_b1n_top_pt_binning->CreateHistogram("hreco_b1n_top_pt");
    // hreco_b2n_top_pt =
    // reco_b2n_top_pt_binning->CreateHistogram("hreco_b2n_top_pt");

    // hreco_c_kk_top_pt =
    // reco_c_kk_top_pt_binning->CreateHistogram("hreco_c_kk_top_pt");
    // hreco_c_rr_top_pt =
    // reco_c_rr_top_pt_binning->CreateHistogram("hreco_c_rr_top_pt");
    // hreco_c_nn_top_pt =
    // reco_c_nn_top_pt_binning->CreateHistogram("hreco_c_nn_top_pt");

    // hreco_c_rk_top_pt =
    // reco_c_rk_top_pt_binning->CreateHistogram("hreco_c_rk_top_pt");
    // hreco_c_kr_top_pt =
    // reco_c_kr_top_pt_binning->CreateHistogram("hreco_c_kr_top_pt");
    // hreco_c_nr_top_pt =
    // reco_c_nr_top_pt_binning->CreateHistogram("hreco_c_nr_top_pt");
    // hreco_c_rn_top_pt =
    // reco_c_rn_top_pt_binning->CreateHistogram("hreco_c_rn_top_pt");
    // hreco_c_nk_top_pt =
    // reco_c_nk_top_pt_binning->CreateHistogram("hreco_c_nk_top_pt");
    // hreco_c_kn_top_pt =
    // reco_c_kn_top_pt_binning->CreateHistogram("hreco_c_kn_top_pt");

    // hreco_c_Prk_top_pt =
    // reco_c_Prk_top_pt_binning->CreateHistogram("hreco_c_Prk_top_pt");
    // hreco_c_Mrk_top_pt =
    // reco_c_Mrk_top_pt_binning->CreateHistogram("hreco_c_Mrk_top_pt");
    // hreco_c_Pnr_top_pt =
    // reco_c_Pnr_top_pt_binning->CreateHistogram("hreco_c_Pnr_top_pt");
    // hreco_c_Mnr_top_pt =
    // reco_c_Mnr_top_pt_binning->CreateHistogram("hreco_c_Mnr_top_pt");
    // hreco_c_Pnk_top_pt =
    // reco_c_Pnk_top_pt_binning->CreateHistogram("hreco_c_Pnk_top_pt");
    // hreco_c_Mnk_top_pt =
    // reco_c_Mnk_top_pt_binning->CreateHistogram("hreco_c_Mnk_top_pt");

    // hreco_c_kj_top_pt =
    // reco_c_kj_top_pt_binning->CreateHistogram("hreco_c_kj_top_pt");
    // hreco_c_rq_top_pt =
    // reco_c_rq_top_pt_binning->CreateHistogram("hreco_c_rq_top_pt");
    // hreco_c_rj_top_pt =
    // reco_c_rj_top_pt_binning->CreateHistogram("hreco_c_rj_top_pt");
    // hreco_c_jr_top_pt =
    // reco_c_jr_top_pt_binning->CreateHistogram("hreco_c_jr_top_pt");

    // hreco_c_Prj_top_pt =
    // reco_c_Prj_top_pt_binning->CreateHistogram("hreco_c_Prj_top_pt");
    // hreco_c_Mrj_top_pt =
    // reco_c_Mrj_top_pt_binning->CreateHistogram("hreco_c_Mrj_top_pt");

    // hreco_c_han_top_pt =
    // reco_c_han_top_pt_binning->CreateHistogram("hreco_c_han_top_pt");
    // hreco_c_sca_top_pt =
    // reco_c_sca_top_pt_binning->CreateHistogram("hreco_c_sca_top_pt");
    // hreco_c_tra_top_pt =
    // reco_c_tra_top_pt_binning->CreateHistogram("hreco_c_tra_top_pt");
    // hreco_c_kjL_top_pt =
    // reco_c_kjL_top_pt_binning->CreateHistogram("hreco_c_kjL_top_pt");
    // hreco_c_rqL_top_pt =
    // reco_c_rqL_top_pt_binning->CreateHistogram("hreco_c_rqL_top_pt");
    // hreco_c_rkP_top_pt =
    // reco_c_rkP_top_pt_binning->CreateHistogram("hreco_c_rkP_top_pt");
    // hreco_c_rkM_top_pt =
    // reco_c_rkM_top_pt_binning->CreateHistogram("hreco_c_rkM_top_pt");
    // hreco_c_nrP_top_pt =
    // reco_c_nrP_top_pt_binning->CreateHistogram("hreco_c_nrP_top_pt");
    // hreco_c_nrM_top_pt =
    // reco_c_nrM_top_pt_binning->CreateHistogram("hreco_c_nrM_top_pt");
    // hreco_c_nkP_top_pt =
    // reco_c_nkP_top_pt_binning->CreateHistogram("hreco_c_nkP_top_pt");
    // hreco_c_nkM_top_pt =
    // reco_c_nkM_top_pt_binning->CreateHistogram("hreco_c_nkM_top_pt");

    // hreco_ll_cHel_top_pt =
    // reco_ll_cHel_top_pt_binning->CreateHistogram("hreco_ll_cHel_top_pt");
    // hreco_ll_cLab_top_pt =
    // reco_ll_cLab_top_pt_binning->CreateHistogram("hreco_ll_cLab_top_pt");
    // hreco_ll_kNorm_top_pt =
    // reco_ll_kNorm_top_pt_binning->CreateHistogram("hreco_ll_kNorm_top_pt");
    // hreco_ll_rNorm_top_pt =
    // reco_ll_rNorm_top_pt_binning->CreateHistogram("hreco_ll_rNorm_top_pt");
    // hreco_llbar_delta_phi_top_pt =
    // reco_llbar_delta_phi_top_pt_binning->CreateHistogram("hreco_llbar_delta_phi_top_pt");
    // hreco_llbar_delta_eta_top_pt =
    // reco_llbar_delta_phi_top_pt_binning->CreateHistogram("hreco_llbar_delta_eta_top_pt");

    // ********************
    // top_scatteringangle_ttbarframe
    // ********************

    // hreco_b1k_top_scatteringangle_ttbarframe =
    // reco_b1k_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_b1k_top_scatteringangle_ttbarframe");
    // hreco_b2k_top_scatteringangle_ttbarframe =
    // reco_b2k_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_b2k_top_scatteringangle_ttbarframe");
    // hreco_b1j_top_scatteringangle_ttbarframe =
    // reco_b1j_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_b1j_top_scatteringangle_ttbarframe");
    // hreco_b2j_top_scatteringangle_ttbarframe =
    // reco_b2j_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_b2j_top_scatteringangle_ttbarframe");
    // hreco_b1r_top_scatteringangle_ttbarframe =
    // reco_b1r_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_b1r_top_scatteringangle_ttbarframe");
    // hreco_b2r_top_scatteringangle_ttbarframe =
    // reco_b2r_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_b2r_top_scatteringangle_ttbarframe");
    // hreco_b1q_top_scatteringangle_ttbarframe =
    // reco_b1q_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_b1q_top_scatteringangle_ttbarframe");
    // hreco_b2q_top_scatteringangle_ttbarframe =
    // reco_b2q_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_b2q_top_scatteringangle_ttbarframe");
    // hreco_b1n_top_scatteringangle_ttbarframe =
    // reco_b1n_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_b1n_top_scatteringangle_ttbarframe");
    // hreco_b2n_top_scatteringangle_ttbarframe =
    // reco_b2n_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_b2n_top_scatteringangle_ttbarframe");

    // hreco_c_kk_top_scatteringangle_ttbarframe =
    // reco_c_kk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_kk_top_scatteringangle_ttbarframe");
    // hreco_c_rr_top_scatteringangle_ttbarframe =
    // reco_c_rr_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_rr_top_scatteringangle_ttbarframe");
    // hreco_c_nn_top_scatteringangle_ttbarframe =
    // reco_c_nn_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_nn_top_scatteringangle_ttbarframe");

    // hreco_c_rk_top_scatteringangle_ttbarframe =
    // reco_c_rk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_rk_top_scatteringangle_ttbarframe");
    // hreco_c_kr_top_scatteringangle_ttbarframe =
    // reco_c_kr_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_kr_top_scatteringangle_ttbarframe");
    // hreco_c_nr_top_scatteringangle_ttbarframe =
    // reco_c_nr_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_nr_top_scatteringangle_ttbarframe");
    // hreco_c_rn_top_scatteringangle_ttbarframe =
    // reco_c_rn_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_rn_top_scatteringangle_ttbarframe");
    // hreco_c_nk_top_scatteringangle_ttbarframe =
    // reco_c_nk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_nk_top_scatteringangle_ttbarframe");
    // hreco_c_kn_top_scatteringangle_ttbarframe =
    // reco_c_kn_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_kn_top_scatteringangle_ttbarframe");

    // hreco_c_Prk_top_scatteringangle_ttbarframe =
    // reco_c_Prk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_Prk_top_scatteringangle_ttbarframe");
    // hreco_c_Mrk_top_scatteringangle_ttbarframe =
    // reco_c_Mrk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_Mrk_top_scatteringangle_ttbarframe");
    // hreco_c_Pnr_top_scatteringangle_ttbarframe =
    // reco_c_Pnr_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_Pnr_top_scatteringangle_ttbarframe");
    // hreco_c_Mnr_top_scatteringangle_ttbarframe =
    // reco_c_Mnr_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_Mnr_top_scatteringangle_ttbarframe");
    // hreco_c_Pnk_top_scatteringangle_ttbarframe =
    // reco_c_Pnk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_Pnk_top_scatteringangle_ttbarframe");
    // hreco_c_Mnk_top_scatteringangle_ttbarframe =
    // reco_c_Mnk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_Mnk_top_scatteringangle_ttbarframe");

    // hreco_c_kj_top_scatteringangle_ttbarframe =
    // reco_c_kj_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_kj_top_scatteringangle_ttbarframe");
    // hreco_c_rq_top_scatteringangle_ttbarframe =
    // reco_c_rq_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_rq_top_scatteringangle_ttbarframe");
    // hreco_c_rj_top_scatteringangle_ttbarframe =
    // reco_c_rj_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_rj_top_scatteringangle_ttbarframe");
    // hreco_c_jr_top_scatteringangle_ttbarframe =
    // reco_c_jr_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_jr_top_scatteringangle_ttbarframe");

    // hreco_c_Prj_top_scatteringangle_ttbarframe =
    // reco_c_Prj_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_Prj_top_scatteringangle_ttbarframe");
    // hreco_c_Mrj_top_scatteringangle_ttbarframe =
    // reco_c_Mrj_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_Mrj_top_scatteringangle_ttbarframe");

    // hreco_c_han_top_scatteringangle_ttbarframe =
    // reco_c_han_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_han_top_scatteringangle_ttbarframe");
    // hreco_c_sca_top_scatteringangle_ttbarframe =
    // reco_c_sca_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_sca_top_scatteringangle_ttbarframe");
    // hreco_c_tra_top_scatteringangle_ttbarframe =
    // reco_c_tra_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_tra_top_scatteringangle_ttbarframe");
    // hreco_c_kjL_top_scatteringangle_ttbarframe =
    // reco_c_kjL_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_kjL_top_scatteringangle_ttbarframe");
    // hreco_c_rqL_top_scatteringangle_ttbarframe =
    // reco_c_rqL_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_rqL_top_scatteringangle_ttbarframe");
    // hreco_c_rkP_top_scatteringangle_ttbarframe =
    // reco_c_rkP_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_rkP_top_scatteringangle_ttbarframe");
    // hreco_c_rkM_top_scatteringangle_ttbarframe =
    // reco_c_rkM_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_rkM_top_scatteringangle_ttbarframe");
    // hreco_c_nrP_top_scatteringangle_ttbarframe =
    // reco_c_nrP_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_nrP_top_scatteringangle_ttbarframe");
    // hreco_c_nrM_top_scatteringangle_ttbarframe =
    // reco_c_nrM_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_nrM_top_scatteringangle_ttbarframe");
    // hreco_c_nkP_top_scatteringangle_ttbarframe =
    // reco_c_nkP_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_nkP_top_scatteringangle_ttbarframe");
    // hreco_c_nkM_top_scatteringangle_ttbarframe =
    // reco_c_nkM_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_c_nkM_top_scatteringangle_ttbarframe");

    // hreco_ll_cHel_top_scatteringangle_ttbarframe =
    // reco_ll_cHel_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_ll_cHel_top_scatteringangle_ttbarframe");
    // hreco_ll_cLab_top_scatteringangle_ttbarframe =
    // reco_ll_cLab_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_ll_cLab_top_scatteringangle_ttbarframe");
    // hreco_ll_kNorm_top_scatteringangle_ttbarframe =
    // reco_ll_kNorm_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_ll_kNorm_top_scatteringangle_ttbarframe");
    // hreco_ll_rNorm_top_scatteringangle_ttbarframe =
    // reco_ll_rNorm_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_ll_rNorm_top_scatteringangle_ttbarframe");
    // hreco_llbar_delta_phi_top_scatteringangle_ttbarframe =
    // reco_llbar_delta_phi_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_llbar_delta_phi_top_scatteringangle_ttbarframe");
    // hreco_llbar_delta_eta_top_scatteringangle_ttbarframe =
    // reco_llbar_delta_phi_top_scatteringangle_ttbarframe_binning->CreateHistogram("hreco_llbar_delta_eta_top_scatteringangle_ttbarframe");

    TString filename = fout->GetName();

    if (doBootstrap && (filename.Contains("ttbarsignal") || filename.Contains("ttbarbg") ||
                        filename.Contains("run"))) {  // only for ttbar or data
        for (int iPE = 0; iPE < nPE; ++iPE) {
            hrecoBootstrap_top_pt[iPE] = reco_top_pt_binning->CreateHistogram(Form("hrecoBootstrap_top_pt%i", iPE));
            hrecoBootstrap_l_pt[iPE] = reco_l_pt_binning->CreateHistogram(Form("hrecoBootstrap_l_pt%i", iPE));
            hrecoBootstrap_lbar_pt[iPE] = reco_lbar_pt_binning->CreateHistogram(Form("hrecoBootstrap_lbar_pt%i", iPE));
            hrecoBootstrap_ttbar_pt[iPE] =
                reco_ttbar_pt_binning->CreateHistogram(Form("hrecoBootstrap_ttbar_pt%i", iPE));
            hrecoBootstrap_ttbar_mass[iPE] =
                reco_ttbar_mass_binning->CreateHistogram(Form("hrecoBootstrap_ttbar_mass%i", iPE));
            hrecoBootstrap_ttbar_delta_phi[iPE] =
                reco_ttbar_delta_phi_binning->CreateHistogram(Form("hrecoBootstrap_ttbar_delta_phi%i", iPE));
            hrecoBootstrap_ttbar_delta_eta[iPE] =
                reco_ttbar_delta_eta_binning->CreateHistogram(Form("hrecoBootstrap_ttbar_delta_eta%i", iPE));
            hrecoBootstrap_ttbar_rapidity[iPE] =
                reco_ttbar_rapidity_binning->CreateHistogram(Form("hrecoBootstrap_ttbar_rapidity%i", iPE));
            hrecoBootstrap_llbar_pt[iPE] =
                reco_llbar_pt_binning->CreateHistogram(Form("hrecoBootstrap_llbar_pt%i", iPE));
            hrecoBootstrap_llbar_mass[iPE] =
                reco_llbar_mass_binning->CreateHistogram(Form("hrecoBootstrap_llbar_mass%i", iPE));

            hrecoBootstrap_b1k[iPE] = reco_b1k_binning->CreateHistogram(Form("hrecoBootstrap_b1k%i", iPE));
            hrecoBootstrap_b2k[iPE] = reco_b2k_binning->CreateHistogram(Form("hrecoBootstrap_b2k%i", iPE));
            hrecoBootstrap_b1j[iPE] = reco_b1j_binning->CreateHistogram(Form("hrecoBootstrap_b1j%i", iPE));
            hrecoBootstrap_b2j[iPE] = reco_b2j_binning->CreateHistogram(Form("hrecoBootstrap_b2j%i", iPE));
            hrecoBootstrap_b1r[iPE] = reco_b1r_binning->CreateHistogram(Form("hrecoBootstrap_b1r%i", iPE));
            hrecoBootstrap_b2r[iPE] = reco_b2r_binning->CreateHistogram(Form("hrecoBootstrap_b2r%i", iPE));
            hrecoBootstrap_b1q[iPE] = reco_b1q_binning->CreateHistogram(Form("hrecoBootstrap_b1q%i", iPE));
            hrecoBootstrap_b2q[iPE] = reco_b2q_binning->CreateHistogram(Form("hrecoBootstrap_b2q%i", iPE));
            hrecoBootstrap_b1n[iPE] = reco_b1n_binning->CreateHistogram(Form("hrecoBootstrap_b1n%i", iPE));
            hrecoBootstrap_b2n[iPE] = reco_b2n_binning->CreateHistogram(Form("hrecoBootstrap_b2n%i", iPE));

            hrecoBootstrap_c_kk[iPE] = reco_c_kk_binning->CreateHistogram(Form("hrecoBootstrap_c_kk%i", iPE));
            hrecoBootstrap_c_rr[iPE] = reco_c_rr_binning->CreateHistogram(Form("hrecoBootstrap_c_rr%i", iPE));
            hrecoBootstrap_c_nn[iPE] = reco_c_nn_binning->CreateHistogram(Form("hrecoBootstrap_c_nn%i", iPE));

            hrecoBootstrap_c_rk[iPE] = reco_c_rk_binning->CreateHistogram(Form("hrecoBootstrap_c_rk%i", iPE));
            hrecoBootstrap_c_kr[iPE] = reco_c_kr_binning->CreateHistogram(Form("hrecoBootstrap_c_kr%i", iPE));
            hrecoBootstrap_c_nr[iPE] = reco_c_nr_binning->CreateHistogram(Form("hrecoBootstrap_c_nr%i", iPE));
            hrecoBootstrap_c_rn[iPE] = reco_c_rn_binning->CreateHistogram(Form("hrecoBootstrap_c_rn%i", iPE));
            hrecoBootstrap_c_nk[iPE] = reco_c_nk_binning->CreateHistogram(Form("hrecoBootstrap_c_nk%i", iPE));
            hrecoBootstrap_c_kn[iPE] = reco_c_kn_binning->CreateHistogram(Form("hrecoBootstrap_c_kn%i", iPE));

            hrecoBootstrap_c_Prk[iPE] = reco_c_Prk_binning->CreateHistogram(Form("hrecoBootstrap_c_Prk%i", iPE));
            hrecoBootstrap_c_Mrk[iPE] = reco_c_Mrk_binning->CreateHistogram(Form("hrecoBootstrap_c_Mrk%i", iPE));
            hrecoBootstrap_c_Pnr[iPE] = reco_c_Pnr_binning->CreateHistogram(Form("hrecoBootstrap_c_Pnr%i", iPE));
            hrecoBootstrap_c_Mnr[iPE] = reco_c_Mnr_binning->CreateHistogram(Form("hrecoBootstrap_c_Mnr%i", iPE));
            hrecoBootstrap_c_Pnk[iPE] = reco_c_Pnk_binning->CreateHistogram(Form("hrecoBootstrap_c_Pnk%i", iPE));
            hrecoBootstrap_c_Mnk[iPE] = reco_c_Mnk_binning->CreateHistogram(Form("hrecoBootstrap_c_Mnk%i", iPE));

            hrecoBootstrap_c_kj[iPE] = reco_c_kj_binning->CreateHistogram(Form("hrecoBootstrap_c_kj%i", iPE));
            hrecoBootstrap_c_rq[iPE] = reco_c_rq_binning->CreateHistogram(Form("hrecoBootstrap_c_rq%i", iPE));
            hrecoBootstrap_c_rj[iPE] = reco_c_rj_binning->CreateHistogram(Form("hrecoBootstrap_c_rj%i", iPE));
            hrecoBootstrap_c_jr[iPE] = reco_c_jr_binning->CreateHistogram(Form("hrecoBootstrap_c_jr%i", iPE));

            hrecoBootstrap_c_Prj[iPE] = reco_c_Prj_binning->CreateHistogram(Form("hrecoBootstrap_c_Prj%i", iPE));
            hrecoBootstrap_c_Mrj[iPE] = reco_c_Mrj_binning->CreateHistogram(Form("hrecoBootstrap_c_Mrj%i", iPE));

            hrecoBootstrap_c_han[iPE] = reco_c_han_binning->CreateHistogram(Form("hrecoBootstrap_c_han%i", iPE));
            hrecoBootstrap_c_sca[iPE] = reco_c_sca_binning->CreateHistogram(Form("hrecoBootstrap_c_sca%i", iPE));
            hrecoBootstrap_c_tra[iPE] = reco_c_tra_binning->CreateHistogram(Form("hrecoBootstrap_c_tra%i", iPE));
            hrecoBootstrap_c_kjL[iPE] = reco_c_kjL_binning->CreateHistogram(Form("hrecoBootstrap_c_kjL%i", iPE));
            hrecoBootstrap_c_rqL[iPE] = reco_c_rqL_binning->CreateHistogram(Form("hrecoBootstrap_c_rqL%i", iPE));
            hrecoBootstrap_c_rkP[iPE] = reco_c_rkP_binning->CreateHistogram(Form("hrecoBootstrap_c_rkP%i", iPE));
            hrecoBootstrap_c_rkM[iPE] = reco_c_rkM_binning->CreateHistogram(Form("hrecoBootstrap_c_rkM%i", iPE));
            hrecoBootstrap_c_nrP[iPE] = reco_c_nrP_binning->CreateHistogram(Form("hrecoBootstrap_c_nrP%i", iPE));
            hrecoBootstrap_c_nrM[iPE] = reco_c_nrM_binning->CreateHistogram(Form("hrecoBootstrap_c_nrM%i", iPE));
            hrecoBootstrap_c_nkP[iPE] = reco_c_nkP_binning->CreateHistogram(Form("hrecoBootstrap_c_nkP%i", iPE));
            hrecoBootstrap_c_nkM[iPE] = reco_c_nkM_binning->CreateHistogram(Form("hrecoBootstrap_c_nkM%i", iPE));

            hrecoBootstrap_ll_cHel[iPE] = reco_ll_cHel_binning->CreateHistogram(Form("hrecoBootstrap_ll_cHel%i", iPE));
            hrecoBootstrap_ll_cLab[iPE] = reco_ll_cLab_binning->CreateHistogram(Form("hrecoBootstrap_ll_cLab%i", iPE));
            hrecoBootstrap_ll_kNorm[iPE] =
                reco_ll_kNorm_binning->CreateHistogram(Form("hrecoBootstrap_ll_kNorm%i", iPE));
            hrecoBootstrap_ll_rNorm[iPE] =
                reco_ll_rNorm_binning->CreateHistogram(Form("hrecoBootstrap_ll_rNorm%i", iPE));
            hrecoBootstrap_llbar_delta_phi[iPE] =
                reco_llbar_delta_phi_binning->CreateHistogram(Form("hrecoBootstrap_llbar_delta_phi%i", iPE));
            hrecoBootstrap_llbar_delta_eta[iPE] =
                reco_llbar_delta_eta_binning->CreateHistogram(Form("hrecoBootstrap_llbar_delta_eta%i", iPE));

            // TODO :
            // 2D pseudo histogram creation
            // hrecoBootstrap_c_kk_mttbar[iPE] =
            // reco_c_kk_mttbar_binning->CreateHistogram(Form("hrecoBootstrap_c_kk_mttbar%i",iPE));
        }
    }

    // **************************
    // 1D gen histogram creation
    // **************************

    // Generator level info -- filled only for MC
    hgen_top_pt = gen_top_pt_binning->CreateHistogram("hgen_top_pt");
    hgen_l_pt = gen_l_pt_binning->CreateHistogram("hgen_l_pt");
    hgen_lbar_pt = gen_lbar_pt_binning->CreateHistogram("hgen_lbar_pt");
    hgen_ttbar_pt = gen_ttbar_pt_binning->CreateHistogram("hgen_ttbar_pt");
    hgen_ttbar_mass = gen_ttbar_mass_binning->CreateHistogram("hgen_ttbar_mass");
    hgen_ttbar_delta_phi = gen_ttbar_delta_phi_binning->CreateHistogram("hgen_ttbar_delta_phi");
    hgen_ttbar_delta_eta = gen_ttbar_delta_eta_binning->CreateHistogram("hgen_ttbar_delta_eta");
    hgen_ttbar_rapidity = gen_ttbar_rapidity_binning->CreateHistogram("hgen_ttbar_rapidity");
    hgen_llbar_pt = gen_llbar_pt_binning->CreateHistogram("hgen_llbar_pt");
    hgen_llbar_mass = gen_llbar_mass_binning->CreateHistogram("hgen_llbar_mass");

    hgen_b1k = gen_b1k_binning->CreateHistogram("hgen_b1k");
    hgen_b2k = gen_b2k_binning->CreateHistogram("hgen_b2k");
    hgen_b1j = gen_b1j_binning->CreateHistogram("hgen_b1j");
    hgen_b2j = gen_b2j_binning->CreateHistogram("hgen_b2j");
    hgen_b1r = gen_b1r_binning->CreateHistogram("hgen_b1r");
    hgen_b2r = gen_b2r_binning->CreateHistogram("hgen_b2r");
    hgen_b1q = gen_b1q_binning->CreateHistogram("hgen_b1q");
    hgen_b2q = gen_b2q_binning->CreateHistogram("hgen_b2q");
    hgen_b1n = gen_b1n_binning->CreateHistogram("hgen_b1n");
    hgen_b2n = gen_b2n_binning->CreateHistogram("hgen_b2n");

    hgen_c_kk = gen_c_kk_binning->CreateHistogram("hgen_c_kk");
    hgen_c_rr = gen_c_rr_binning->CreateHistogram("hgen_c_rr");
    hgen_c_nn = gen_c_nn_binning->CreateHistogram("hgen_c_nn");

    hgen_c_rk = gen_c_rk_binning->CreateHistogram("hgen_c_rk");
    hgen_c_kr = gen_c_kr_binning->CreateHistogram("hgen_c_kr");
    hgen_c_nr = gen_c_nr_binning->CreateHistogram("hgen_c_nr");
    hgen_c_rn = gen_c_rn_binning->CreateHistogram("hgen_c_rn");
    hgen_c_nk = gen_c_nk_binning->CreateHistogram("hgen_c_nk");
    hgen_c_kn = gen_c_kn_binning->CreateHistogram("hgen_c_kn");

    hgen_c_Prk = gen_c_Prk_binning->CreateHistogram("hgen_c_Prk");
    hgen_c_Mrk = gen_c_Mrk_binning->CreateHistogram("hgen_c_Mrk");
    hgen_c_Pnr = gen_c_Pnr_binning->CreateHistogram("hgen_c_Pnr");
    hgen_c_Mnr = gen_c_Mnr_binning->CreateHistogram("hgen_c_Mnr");
    hgen_c_Pnk = gen_c_Pnk_binning->CreateHistogram("hgen_c_Pnk");
    hgen_c_Mnk = gen_c_Mnk_binning->CreateHistogram("hgen_c_Mnk");

    hgen_c_kj = gen_c_kj_binning->CreateHistogram("hgen_c_kj");
    hgen_c_rq = gen_c_rq_binning->CreateHistogram("hgen_c_rq");
    hgen_c_rj = gen_c_rj_binning->CreateHistogram("hgen_c_rj");
    hgen_c_jr = gen_c_jr_binning->CreateHistogram("hgen_c_jr");

    hgen_c_Prj = gen_c_Prj_binning->CreateHistogram("hgen_c_Prj");
    hgen_c_Mrj = gen_c_Mrj_binning->CreateHistogram("hgen_c_Mrj");

    hgen_c_han = gen_c_han_binning->CreateHistogram("hgen_c_han");
    hgen_c_sca = gen_c_sca_binning->CreateHistogram("hgen_c_sca");
    hgen_c_tra = gen_c_tra_binning->CreateHistogram("hgen_c_tra");
    hgen_c_kjL = gen_c_kjL_binning->CreateHistogram("hgen_c_kjL");
    hgen_c_rqL = gen_c_rqL_binning->CreateHistogram("hgen_c_rqL");
    hgen_c_rkP = gen_c_rkP_binning->CreateHistogram("hgen_c_rkP");
    hgen_c_rkM = gen_c_rkM_binning->CreateHistogram("hgen_c_rkM");
    hgen_c_nrP = gen_c_nrP_binning->CreateHistogram("hgen_c_nrP");
    hgen_c_nrM = gen_c_nrM_binning->CreateHistogram("hgen_c_nrM");
    hgen_c_nkP = gen_c_nkP_binning->CreateHistogram("hgen_c_nkP");
    hgen_c_nkM = gen_c_nkM_binning->CreateHistogram("hgen_c_nkM");

    hgen_ll_cHel = gen_ll_cHel_binning->CreateHistogram("hgen_ll_cHel");
    hgen_ll_cLab = gen_ll_cLab_binning->CreateHistogram("hgen_ll_cLab");
    hgen_ll_kNorm = gen_ll_kNorm_binning->CreateHistogram("hgen_ll_kNorm");
    hgen_ll_rNorm = gen_ll_rNorm_binning->CreateHistogram("hgen_ll_rNorm");
    hgen_llbar_delta_phi = gen_llbar_delta_phi_binning->CreateHistogram("hgen_llbar_delta_phi");
    hgen_llbar_delta_eta = gen_llbar_delta_eta_binning->CreateHistogram("hgen_llbar_delta_eta");

    // **************************
    // 2D gen histogram creation
    // **************************

    // Amandeep : 2D gen histogram creation
    // hgen_c_kk_mttbar =
    // gen_c_kk_mttbar_binning->CreateHistogram("hgen_c_kk_mttbar");

    // ******
    // mttbar
    // ******

    hgen_b1k_exjet = gen_b1k_exjet_binning->CreateHistogram("hgen_b1k_exjet");
    hgen_b2k_exjet = gen_b2k_exjet_binning->CreateHistogram("hgen_b2k_exjet");
    hgen_b1j_exjet = gen_b1j_exjet_binning->CreateHistogram("hgen_b1j_exjet");
    hgen_b2j_exjet = gen_b2j_exjet_binning->CreateHistogram("hgen_b2j_exjet");
    hgen_b1r_exjet = gen_b1r_exjet_binning->CreateHistogram("hgen_b1r_exjet");
    hgen_b2r_exjet = gen_b2r_exjet_binning->CreateHistogram("hgen_b2r_exjet");
    hgen_b1q_exjet = gen_b1q_exjet_binning->CreateHistogram("hgen_b1q_exjet");
    hgen_b2q_exjet = gen_b2q_exjet_binning->CreateHistogram("hgen_b2q_exjet");
    hgen_b1n_exjet = gen_b1n_exjet_binning->CreateHistogram("hgen_b1n_exjet");
    hgen_b2n_exjet = gen_b2n_exjet_binning->CreateHistogram("hgen_b2n_exjet");

    hgen_c_kk_exjet = gen_c_kk_exjet_binning->CreateHistogram("hgen_c_kk_exjet");
    hgen_c_rr_exjet = gen_c_rr_exjet_binning->CreateHistogram("hgen_c_rr_exjet");
    hgen_c_nn_exjet = gen_c_nn_exjet_binning->CreateHistogram("hgen_c_nn_exjet");

    hgen_c_rk_exjet = gen_c_rk_exjet_binning->CreateHistogram("hgen_c_rk_exjet");
    hgen_c_kr_exjet = gen_c_kr_exjet_binning->CreateHistogram("hgen_c_kr_exjet");
    hgen_c_nr_exjet = gen_c_nr_exjet_binning->CreateHistogram("hgen_c_nr_exjet");
    hgen_c_rn_exjet = gen_c_rn_exjet_binning->CreateHistogram("hgen_c_rn_exjet");
    hgen_c_nk_exjet = gen_c_nk_exjet_binning->CreateHistogram("hgen_c_nk_exjet");
    hgen_c_kn_exjet = gen_c_kn_exjet_binning->CreateHistogram("hgen_c_kn_exjet");

    hgen_c_Prk_exjet = gen_c_Prk_exjet_binning->CreateHistogram("hgen_c_Prk_exjet");
    hgen_c_Mrk_exjet = gen_c_Mrk_exjet_binning->CreateHistogram("hgen_c_Mrk_exjet");
    hgen_c_Pnr_exjet = gen_c_Pnr_exjet_binning->CreateHistogram("hgen_c_Pnr_exjet");
    hgen_c_Mnr_exjet = gen_c_Mnr_exjet_binning->CreateHistogram("hgen_c_Mnr_exjet");
    hgen_c_Pnk_exjet = gen_c_Pnk_exjet_binning->CreateHistogram("hgen_c_Pnk_exjet");
    hgen_c_Mnk_exjet = gen_c_Mnk_exjet_binning->CreateHistogram("hgen_c_Mnk_exjet");

    hgen_c_kj_exjet = gen_c_kj_exjet_binning->CreateHistogram("hgen_c_kj_exjet");
    hgen_c_rq_exjet = gen_c_rq_exjet_binning->CreateHistogram("hgen_c_rq_exjet");
    hgen_c_rj_exjet = gen_c_rj_exjet_binning->CreateHistogram("hgen_c_rj_exjet");
    hgen_c_jr_exjet = gen_c_jr_exjet_binning->CreateHistogram("hgen_c_jr_exjet");

    hgen_c_Prj_exjet = gen_c_Prj_exjet_binning->CreateHistogram("hgen_c_Prj_exjet");
    hgen_c_Mrj_exjet = gen_c_Mrj_exjet_binning->CreateHistogram("hgen_c_Mrj_exjet");

    hgen_c_han_exjet = gen_c_han_exjet_binning->CreateHistogram("hgen_c_han_exjet");
    hgen_c_sca_exjet = gen_c_sca_exjet_binning->CreateHistogram("hgen_c_sca_exjet");
    hgen_c_tra_exjet = gen_c_tra_exjet_binning->CreateHistogram("hgen_c_tra_exjet");
    hgen_c_kjL_exjet = gen_c_kjL_exjet_binning->CreateHistogram("hgen_c_kjL_exjet");
    hgen_c_rqL_exjet = gen_c_rqL_exjet_binning->CreateHistogram("hgen_c_rqL_exjet");
    hgen_c_rkP_exjet = gen_c_rkP_exjet_binning->CreateHistogram("hgen_c_rkP_exjet");
    hgen_c_rkM_exjet = gen_c_rkM_exjet_binning->CreateHistogram("hgen_c_rkM_exjet");
    hgen_c_nrP_exjet = gen_c_nrP_exjet_binning->CreateHistogram("hgen_c_nrP_exjet");
    hgen_c_nrM_exjet = gen_c_nrM_exjet_binning->CreateHistogram("hgen_c_nrM_exjet");
    hgen_c_nkP_exjet = gen_c_nkP_exjet_binning->CreateHistogram("hgen_c_nkP_exjet");
    hgen_c_nkM_exjet = gen_c_nkM_exjet_binning->CreateHistogram("hgen_c_nkM_exjet");

    hgen_ll_cHel_exjet = gen_ll_cHel_exjet_binning->CreateHistogram("hgen_ll_cHel_exjet");
    hgen_ll_cLab_exjet = gen_ll_cLab_exjet_binning->CreateHistogram("hgen_ll_cLab_exjet");
    hgen_ll_kNorm_exjet = gen_ll_kNorm_exjet_binning->CreateHistogram("hgen_ll_kNorm_exjet");
    hgen_ll_rNorm_exjet = gen_ll_rNorm_exjet_binning->CreateHistogram("hgen_ll_rNorm_exjet");
    hgen_llbar_delta_phi_exjet = gen_llbar_delta_phi_exjet_binning->CreateHistogram("hgen_llbar_delta_phi_exjet");
    hgen_llbar_delta_eta_exjet = gen_llbar_delta_eta_exjet_binning->CreateHistogram("hgen_llbar_delta_eta_exjet");

    // ******
    // top_pt
    // ******

    // hgen_b1k_top_pt =
    // gen_b1k_top_pt_binning->CreateHistogram("hgen_b1k_top_pt");
    // hgen_b2k_top_pt =
    // gen_b2k_top_pt_binning->CreateHistogram("hgen_b2k_top_pt");
    // hgen_b1j_top_pt =
    // gen_b1j_top_pt_binning->CreateHistogram("hgen_b1j_top_pt");
    // hgen_b2j_top_pt =
    // gen_b2j_top_pt_binning->CreateHistogram("hgen_b2j_top_pt");
    // hgen_b1r_top_pt =
    // gen_b1r_top_pt_binning->CreateHistogram("hgen_b1r_top_pt");
    // hgen_b2r_top_pt =
    // gen_b2r_top_pt_binning->CreateHistogram("hgen_b2r_top_pt");
    // hgen_b1q_top_pt =
    // gen_b1q_top_pt_binning->CreateHistogram("hgen_b1q_top_pt");
    // hgen_b2q_top_pt =
    // gen_b2q_top_pt_binning->CreateHistogram("hgen_b2q_top_pt");
    // hgen_b1n_top_pt =
    // gen_b1n_top_pt_binning->CreateHistogram("hgen_b1n_top_pt");
    // hgen_b2n_top_pt =
    // gen_b2n_top_pt_binning->CreateHistogram("hgen_b2n_top_pt");

    // hgen_c_kk_top_pt =
    // gen_c_kk_top_pt_binning->CreateHistogram("hgen_c_kk_top_pt");
    // hgen_c_rr_top_pt =
    // gen_c_rr_top_pt_binning->CreateHistogram("hgen_c_rr_top_pt");
    // hgen_c_nn_top_pt =
    // gen_c_nn_top_pt_binning->CreateHistogram("hgen_c_nn_top_pt");

    // hgen_c_rk_top_pt =
    // gen_c_rk_top_pt_binning->CreateHistogram("hgen_c_rk_top_pt");
    // hgen_c_kr_top_pt =
    // gen_c_kr_top_pt_binning->CreateHistogram("hgen_c_kr_top_pt");
    // hgen_c_nr_top_pt =
    // gen_c_nr_top_pt_binning->CreateHistogram("hgen_c_nr_top_pt");
    // hgen_c_rn_top_pt =
    // gen_c_rn_top_pt_binning->CreateHistogram("hgen_c_rn_top_pt");
    // hgen_c_nk_top_pt =
    // gen_c_nk_top_pt_binning->CreateHistogram("hgen_c_nk_top_pt");
    // hgen_c_kn_top_pt =
    // gen_c_kn_top_pt_binning->CreateHistogram("hgen_c_kn_top_pt");

    // hgen_c_Prk_top_pt =
    // gen_c_Prk_top_pt_binning->CreateHistogram("hgen_c_Prk_top_pt");
    // hgen_c_Mrk_top_pt =
    // gen_c_Mrk_top_pt_binning->CreateHistogram("hgen_c_Mrk_top_pt");
    // hgen_c_Pnr_top_pt =
    // gen_c_Pnr_top_pt_binning->CreateHistogram("hgen_c_Pnr_top_pt");
    // hgen_c_Mnr_top_pt =
    // gen_c_Mnr_top_pt_binning->CreateHistogram("hgen_c_Mnr_top_pt");
    // hgen_c_Pnk_top_pt =
    // gen_c_Pnk_top_pt_binning->CreateHistogram("hgen_c_Pnk_top_pt");
    // hgen_c_Mnk_top_pt =
    // gen_c_Mnk_top_pt_binning->CreateHistogram("hgen_c_Mnk_top_pt");

    // hgen_c_kj_top_pt =
    // gen_c_kj_top_pt_binning->CreateHistogram("hgen_c_kj_top_pt");
    // hgen_c_rq_top_pt =
    // gen_c_rq_top_pt_binning->CreateHistogram("hgen_c_rq_top_pt");
    // hgen_c_rj_top_pt =
    // gen_c_rj_top_pt_binning->CreateHistogram("hgen_c_rj_top_pt");
    // hgen_c_jr_top_pt =
    // gen_c_jr_top_pt_binning->CreateHistogram("hgen_c_jr_top_pt");

    // hgen_c_Prj_top_pt =
    // gen_c_Prj_top_pt_binning->CreateHistogram("hgen_c_Prj_top_pt");
    // hgen_c_Mrj_top_pt =
    // gen_c_Mrj_top_pt_binning->CreateHistogram("hgen_c_Mrj_top_pt");

    // hgen_c_han_top_pt =
    // gen_c_han_top_pt_binning->CreateHistogram("hgen_c_han_top_pt");
    // hgen_c_sca_top_pt =
    // gen_c_sca_top_pt_binning->CreateHistogram("hgen_c_sca_top_pt");
    // hgen_c_tra_top_pt =
    // gen_c_tra_top_pt_binning->CreateHistogram("hgen_c_tra_top_pt");
    // hgen_c_kjL_top_pt =
    // gen_c_kjL_top_pt_binning->CreateHistogram("hgen_c_kjL_top_pt");
    // hgen_c_rqL_top_pt =
    // gen_c_rqL_top_pt_binning->CreateHistogram("hgen_c_rqL_top_pt");
    // hgen_c_rkP_top_pt =
    // gen_c_rkP_top_pt_binning->CreateHistogram("hgen_c_rkP_top_pt");
    // hgen_c_rkM_top_pt =
    // gen_c_rkM_top_pt_binning->CreateHistogram("hgen_c_rkM_top_pt");
    // hgen_c_nrP_top_pt =
    // gen_c_nrP_top_pt_binning->CreateHistogram("hgen_c_nrP_top_pt");
    // hgen_c_nrM_top_pt =
    // gen_c_nrM_top_pt_binning->CreateHistogram("hgen_c_nrM_top_pt");
    // hgen_c_nkP_top_pt =
    // gen_c_nkP_top_pt_binning->CreateHistogram("hgen_c_nkP_top_pt");
    // hgen_c_nkM_top_pt =
    // gen_c_nkM_top_pt_binning->CreateHistogram("hgen_c_nkM_top_pt");

    // hgen_ll_cHel_top_pt  =
    // gen_ll_cHel_top_pt_binning->CreateHistogram("hgen_ll_cHel_top_pt");
    // hgen_ll_cLab_top_pt  =
    // gen_ll_cLab_top_pt_binning->CreateHistogram("hgen_ll_cLab_top_pt");
    // hgen_ll_kNorm_top_pt =
    // gen_ll_kNorm_top_pt_binning->CreateHistogram("hgen_ll_kNorm_top_pt");
    // hgen_ll_rNorm_top_pt =
    // gen_ll_rNorm_top_pt_binning->CreateHistogram("hgen_ll_rNorm_top_pt");
    // hgen_llbar_delta_phi_top_pt =
    // gen_llbar_delta_phi_top_pt_binning->CreateHistogram("hgen_llbar_delta_phi_top_pt");
    // hgen_llbar_delta_eta_top_pt =
    // gen_llbar_delta_eta_top_pt_binning->CreateHistogram("hgen_llbar_delta_eta_top_pt");

    // ********************
    // top_scatteringangle_ttbarframe
    // ********************

    // hgen_b1k_top_scatteringangle_ttbarframe =
    // gen_b1k_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_b1k_top_scatteringangle_ttbarframe");
    // hgen_b2k_top_scatteringangle_ttbarframe =
    // gen_b2k_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_b2k_top_scatteringangle_ttbarframe");
    // hgen_b1j_top_scatteringangle_ttbarframe =
    // gen_b1j_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_b1j_top_scatteringangle_ttbarframe");
    // hgen_b2j_top_scatteringangle_ttbarframe =
    // gen_b2j_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_b2j_top_scatteringangle_ttbarframe");
    // hgen_b1r_top_scatteringangle_ttbarframe =
    // gen_b1r_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_b1r_top_scatteringangle_ttbarframe");
    // hgen_b2r_top_scatteringangle_ttbarframe =
    // gen_b2r_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_b2r_top_scatteringangle_ttbarframe");
    // hgen_b1q_top_scatteringangle_ttbarframe =
    // gen_b1q_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_b1q_top_scatteringangle_ttbarframe");
    // hgen_b2q_top_scatteringangle_ttbarframe =
    // gen_b2q_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_b2q_top_scatteringangle_ttbarframe");
    // hgen_b1n_top_scatteringangle_ttbarframe =
    // gen_b1n_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_b1n_top_scatteringangle_ttbarframe");
    // hgen_b2n_top_scatteringangle_ttbarframe =
    // gen_b2n_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_b2n_top_scatteringangle_ttbarframe");

    // hgen_c_kk_top_scatteringangle_ttbarframe =
    // gen_c_kk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_kk_top_scatteringangle_ttbarframe");
    // hgen_c_rr_top_scatteringangle_ttbarframe =
    // gen_c_rr_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_rr_top_scatteringangle_ttbarframe");
    // hgen_c_nn_top_scatteringangle_ttbarframe =
    // gen_c_nn_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_nn_top_scatteringangle_ttbarframe");

    // hgen_c_rk_top_scatteringangle_ttbarframe =
    // gen_c_rk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_rk_top_scatteringangle_ttbarframe");
    // hgen_c_kr_top_scatteringangle_ttbarframe =
    // gen_c_kr_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_kr_top_scatteringangle_ttbarframe");
    // hgen_c_nr_top_scatteringangle_ttbarframe =
    // gen_c_nr_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_nr_top_scatteringangle_ttbarframe");
    // hgen_c_rn_top_scatteringangle_ttbarframe =
    // gen_c_rn_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_rn_top_scatteringangle_ttbarframe");
    // hgen_c_nk_top_scatteringangle_ttbarframe =
    // gen_c_nk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_nk_top_scatteringangle_ttbarframe");
    // hgen_c_kn_top_scatteringangle_ttbarframe =
    // gen_c_kn_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_kn_top_scatteringangle_ttbarframe");

    // hgen_c_Prk_top_scatteringangle_ttbarframe =
    // gen_c_Prk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_Prk_top_scatteringangle_ttbarframe");
    // hgen_c_Mrk_top_scatteringangle_ttbarframe =
    // gen_c_Mrk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_Mrk_top_scatteringangle_ttbarframe");
    // hgen_c_Pnr_top_scatteringangle_ttbarframe =
    // gen_c_Pnr_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_Pnr_top_scatteringangle_ttbarframe");
    // hgen_c_Mnr_top_scatteringangle_ttbarframe =
    // gen_c_Mnr_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_Mnr_top_scatteringangle_ttbarframe");
    // hgen_c_Pnk_top_scatteringangle_ttbarframe =
    // gen_c_Pnk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_Pnk_top_scatteringangle_ttbarframe");
    // hgen_c_Mnk_top_scatteringangle_ttbarframe =
    // gen_c_Mnk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_Mnk_top_scatteringangle_ttbarframe");

    // hgen_c_kj_top_scatteringangle_ttbarframe =
    // gen_c_kj_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_kj_top_scatteringangle_ttbarframe");
    // hgen_c_rq_top_scatteringangle_ttbarframe =
    // gen_c_rq_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_rq_top_scatteringangle_ttbarframe");
    // hgen_c_rj_top_scatteringangle_ttbarframe =
    // gen_c_rj_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_rj_top_scatteringangle_ttbarframe");
    // hgen_c_jr_top_scatteringangle_ttbarframe =
    // gen_c_jr_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_jr_top_scatteringangle_ttbarframe");

    // hgen_c_Prj_top_scatteringangle_ttbarframe =
    // gen_c_Prj_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_Prj_top_scatteringangle_ttbarframe");
    // hgen_c_Mrj_top_scatteringangle_ttbarframe =
    // gen_c_Mrj_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_Mrj_top_scatteringangle_ttbarframe");

    // hgen_c_han_top_scatteringangle_ttbarframe =
    // gen_c_han_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_han_top_scatteringangle_ttbarframe");
    // hgen_c_sca_top_scatteringangle_ttbarframe =
    // gen_c_sca_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_sca_top_scatteringangle_ttbarframe");
    // hgen_c_tra_top_scatteringangle_ttbarframe =
    // gen_c_tra_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_tra_top_scatteringangle_ttbarframe");
    // hgen_c_kjL_top_scatteringangle_ttbarframe =
    // gen_c_kjL_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_kjL_top_scatteringangle_ttbarframe");
    // hgen_c_rqL_top_scatteringangle_ttbarframe =
    // gen_c_rqL_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_rqL_top_scatteringangle_ttbarframe");
    // hgen_c_rkP_top_scatteringangle_ttbarframe =
    // gen_c_rkP_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_rkP_top_scatteringangle_ttbarframe");
    // hgen_c_rkM_top_scatteringangle_ttbarframe =
    // gen_c_rkM_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_rkM_top_scatteringangle_ttbarframe");
    // hgen_c_nrP_top_scatteringangle_ttbarframe =
    // gen_c_nrP_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_nrP_top_scatteringangle_ttbarframe");
    // hgen_c_nrM_top_scatteringangle_ttbarframe =
    // gen_c_nrM_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_nrM_top_scatteringangle_ttbarframe");
    // hgen_c_nkP_top_scatteringangle_ttbarframe =
    // gen_c_nkP_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_nkP_top_scatteringangle_ttbarframe");
    // hgen_c_nkM_top_scatteringangle_ttbarframe =
    // gen_c_nkM_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_c_nkM_top_scatteringangle_ttbarframe");

    // hgen_ll_cHel_top_scatteringangle_ttbarframe  =
    // gen_ll_cHel_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_ll_cHel_top_scatteringangle_ttbarframe");
    // hgen_ll_cLab_top_scatteringangle_ttbarframe  =
    // gen_ll_cLab_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_ll_cLab_top_scatteringangle_ttbarframe");
    // hgen_ll_kNorm_top_scatteringangle_ttbarframe =
    // gen_ll_kNorm_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_ll_kNorm_top_scatteringangle_ttbarframe");
    // hgen_ll_rNorm_top_scatteringangle_ttbarframe =
    // gen_ll_rNorm_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_ll_rNorm_top_scatteringangle_ttbarframe");
    // hgen_llbar_delta_phi_top_scatteringangle_ttbarframe =
    // gen_llbar_delta_phi_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_llbar_delta_phi_top_scatteringangle_ttbarframe");
    // hgen_llbar_delta_eta_top_scatteringangle_ttbarframe =
    // gen_llbar_delta_eta_top_scatteringangle_ttbarframe_binning->CreateHistogram("hgen_llbar_delta_eta_top_scatteringangle_ttbarframe");

    // ****************************
    // 1D visgen histogram creation
    // ****************************

    // Visible generated info

    // Kinematic
    hvisgen_top_pt = gen_top_pt_binning->CreateHistogram("hvisgen_top_pt");
    hvisgen_l_pt = gen_l_pt_binning->CreateHistogram("hvisgen_l_pt");
    hvisgen_lbar_pt = gen_lbar_pt_binning->CreateHistogram("hvisgen_lbar_pt");
    hvisgen_ttbar_pt = gen_ttbar_pt_binning->CreateHistogram("hvisgen_ttbar_pt");
    hvisgen_ttbar_mass = gen_ttbar_mass_binning->CreateHistogram("hvisgen_ttbar_mass");
    hvisgen_ttbar_delta_phi = gen_ttbar_delta_phi_binning->CreateHistogram("hvisgen_ttbar_delta_phi");
    hvisgen_ttbar_delta_eta = gen_ttbar_delta_eta_binning->CreateHistogram("hvisgen_ttbar_delta_eta");
    hvisgen_ttbar_rapidity = gen_ttbar_rapidity_binning->CreateHistogram("hvisgen_ttbar_rapidity");
    hvisgen_llbar_pt = gen_llbar_pt_binning->CreateHistogram("hvisgen_llbar_pt");
    hvisgen_llbar_mass = gen_llbar_mass_binning->CreateHistogram("hvisgen_llbar_mass");

    // Spin corr
    hvisgen_b1k = gen_b1k_binning->CreateHistogram("hvisgen_b1k");
    hvisgen_b2k = gen_b2k_binning->CreateHistogram("hvisgen_b2k");
    hvisgen_b1j = gen_b1j_binning->CreateHistogram("hvisgen_b1j");
    hvisgen_b2j = gen_b2j_binning->CreateHistogram("hvisgen_b2j");
    hvisgen_b1r = gen_b1r_binning->CreateHistogram("hvisgen_b1r");
    hvisgen_b2r = gen_b2r_binning->CreateHistogram("hvisgen_b2r");
    hvisgen_b1q = gen_b1q_binning->CreateHistogram("hvisgen_b1q");
    hvisgen_b2q = gen_b2q_binning->CreateHistogram("hvisgen_b2q");
    hvisgen_b1n = gen_b1n_binning->CreateHistogram("hvisgen_b1n");
    hvisgen_b2n = gen_b2n_binning->CreateHistogram("hvisgen_b2n");

    hvisgen_c_kk = gen_c_kk_binning->CreateHistogram("hvisgen_c_kk");
    hvisgen_c_rr = gen_c_rr_binning->CreateHistogram("hvisgen_c_rr");
    hvisgen_c_nn = gen_c_nn_binning->CreateHistogram("hvisgen_c_nn");

    hvisgen_c_rk = gen_c_rk_binning->CreateHistogram("hvisgen_c_rk");
    hvisgen_c_kr = gen_c_kr_binning->CreateHistogram("hvisgen_c_kr");
    hvisgen_c_nr = gen_c_nr_binning->CreateHistogram("hvisgen_c_nr");
    hvisgen_c_rn = gen_c_rn_binning->CreateHistogram("hvisgen_c_rn");
    hvisgen_c_nk = gen_c_nk_binning->CreateHistogram("hvisgen_c_nk");
    hvisgen_c_kn = gen_c_kn_binning->CreateHistogram("hvisgen_c_kn");

    hvisgen_c_Prk = gen_c_Prk_binning->CreateHistogram("hvisgen_c_Prk");
    hvisgen_c_Mrk = gen_c_Mrk_binning->CreateHistogram("hvisgen_c_Mrk");
    hvisgen_c_Pnr = gen_c_Pnr_binning->CreateHistogram("hvisgen_c_Pnr");
    hvisgen_c_Mnr = gen_c_Mnr_binning->CreateHistogram("hvisgen_c_Mnr");
    hvisgen_c_Pnk = gen_c_Pnk_binning->CreateHistogram("hvisgen_c_Pnk");
    hvisgen_c_Mnk = gen_c_Mnk_binning->CreateHistogram("hvisgen_c_Mnk");

    hvisgen_c_kj = gen_c_kj_binning->CreateHistogram("hvisgen_c_kj");
    hvisgen_c_rq = gen_c_rq_binning->CreateHistogram("hvisgen_c_rq");
    hvisgen_c_rj = gen_c_rj_binning->CreateHistogram("hvisgen_c_rj");
    hvisgen_c_jr = gen_c_jr_binning->CreateHistogram("hvisgen_c_jr");

    hvisgen_c_Prj = gen_c_Prj_binning->CreateHistogram("hvisgen_c_Prj");
    hvisgen_c_Mrj = gen_c_Mrj_binning->CreateHistogram("hvisgen_c_Mrj");

    hvisgen_c_han = gen_c_han_binning->CreateHistogram("hvisgen_c_han");
    hvisgen_c_sca = gen_c_sca_binning->CreateHistogram("hvisgen_c_sca");
    hvisgen_c_tra = gen_c_tra_binning->CreateHistogram("hvisgen_c_tra");
    hvisgen_c_kjL = gen_c_kjL_binning->CreateHistogram("hvisgen_c_kjL");
    hvisgen_c_rqL = gen_c_rqL_binning->CreateHistogram("hvisgen_c_rqL");
    hvisgen_c_rkP = gen_c_rkP_binning->CreateHistogram("hvisgen_c_rkP");
    hvisgen_c_rkM = gen_c_rkM_binning->CreateHistogram("hvisgen_c_rkM");
    hvisgen_c_nrP = gen_c_nrP_binning->CreateHistogram("hvisgen_c_nrP");
    hvisgen_c_nrM = gen_c_nrM_binning->CreateHistogram("hvisgen_c_nrM");
    hvisgen_c_nkP = gen_c_nkP_binning->CreateHistogram("hvisgen_c_nkP");
    hvisgen_c_nkM = gen_c_nkM_binning->CreateHistogram("hvisgen_c_nkM");

    hvisgen_ll_cHel = gen_ll_cHel_binning->CreateHistogram("hvisgen_ll_cHel");
    hvisgen_ll_cLab = gen_ll_cLab_binning->CreateHistogram("hvisgen_ll_cLab");
    hvisgen_ll_kNorm = gen_ll_kNorm_binning->CreateHistogram("hvisgen_ll_kNorm");
    hvisgen_ll_rNorm = gen_ll_rNorm_binning->CreateHistogram("hvisgen_ll_rNorm");
    hvisgen_llbar_delta_phi = gen_llbar_delta_phi_binning->CreateHistogram("hvisgen_llbar_delta_phi");
    hvisgen_llbar_delta_eta = gen_llbar_delta_eta_binning->CreateHistogram("hvisgen_llbar_delta_eta");

    // ****************************
    // 2D visgen histogram creation
    // ****************************

    // Amandeep : 2D visible gen histogram creation
    // hvisgen_c_kk_mttbar =
    // gen_c_kk_mttbar_binning->CreateHistogram("hvisgen_c_kk_mttbar");

    // ******
    // mttbar
    // ******

    hvisgen_b1k_exjet = gen_b1k_exjet_binning->CreateHistogram("hvisgen_b1k_exjet");
    hvisgen_b2k_exjet = gen_b2k_exjet_binning->CreateHistogram("hvisgen_b2k_exjet");
    hvisgen_b1j_exjet = gen_b1j_exjet_binning->CreateHistogram("hvisgen_b1j_exjet");
    hvisgen_b2j_exjet = gen_b2j_exjet_binning->CreateHistogram("hvisgen_b2j_exjet");
    hvisgen_b1r_exjet = gen_b1r_exjet_binning->CreateHistogram("hvisgen_b1r_exjet");
    hvisgen_b2r_exjet = gen_b2r_exjet_binning->CreateHistogram("hvisgen_b2r_exjet");
    hvisgen_b1q_exjet = gen_b1q_exjet_binning->CreateHistogram("hvisgen_b1q_exjet");
    hvisgen_b2q_exjet = gen_b2q_exjet_binning->CreateHistogram("hvisgen_b2q_exjet");
    hvisgen_b1n_exjet = gen_b1n_exjet_binning->CreateHistogram("hvisgen_b1n_exjet");
    hvisgen_b2n_exjet = gen_b2n_exjet_binning->CreateHistogram("hvisgen_b2n_exjet");

    hvisgen_c_kk_exjet = gen_c_kk_exjet_binning->CreateHistogram("hvisgen_c_kk_exjet");
    hvisgen_c_rr_exjet = gen_c_rr_exjet_binning->CreateHistogram("hvisgen_c_rr_exjet");
    hvisgen_c_nn_exjet = gen_c_nn_exjet_binning->CreateHistogram("hvisgen_c_nn_exjet");

    hvisgen_c_rk_exjet = gen_c_rk_exjet_binning->CreateHistogram("hvisgen_c_rk_exjet");
    hvisgen_c_kr_exjet = gen_c_kr_exjet_binning->CreateHistogram("hvisgen_c_kr_exjet");
    hvisgen_c_nr_exjet = gen_c_nr_exjet_binning->CreateHistogram("hvisgen_c_nr_exjet");
    hvisgen_c_rn_exjet = gen_c_rn_exjet_binning->CreateHistogram("hvisgen_c_rn_exjet");
    hvisgen_c_nk_exjet = gen_c_nk_exjet_binning->CreateHistogram("hvisgen_c_nk_exjet");
    hvisgen_c_kn_exjet = gen_c_kn_exjet_binning->CreateHistogram("hvisgen_c_kn_exjet");

    hvisgen_c_Prk_exjet = gen_c_Prk_exjet_binning->CreateHistogram("hvisgen_c_Prk_exjet");
    hvisgen_c_Mrk_exjet = gen_c_Mrk_exjet_binning->CreateHistogram("hvisgen_c_Mrk_exjet");
    hvisgen_c_Pnr_exjet = gen_c_Pnr_exjet_binning->CreateHistogram("hvisgen_c_Pnr_exjet");
    hvisgen_c_Mnr_exjet = gen_c_Mnr_exjet_binning->CreateHistogram("hvisgen_c_Mnr_exjet");
    hvisgen_c_Pnk_exjet = gen_c_Pnk_exjet_binning->CreateHistogram("hvisgen_c_Pnk_exjet");
    hvisgen_c_Mnk_exjet = gen_c_Mnk_exjet_binning->CreateHistogram("hvisgen_c_Mnk_exjet");

    hvisgen_c_kj_exjet = gen_c_kj_exjet_binning->CreateHistogram("hvisgen_c_kj_exjet");
    hvisgen_c_rq_exjet = gen_c_rq_exjet_binning->CreateHistogram("hvisgen_c_rq_exjet");
    hvisgen_c_rj_exjet = gen_c_rj_exjet_binning->CreateHistogram("hvisgen_c_rj_exjet");
    hvisgen_c_jr_exjet = gen_c_jr_exjet_binning->CreateHistogram("hvisgen_c_jr_exjet");

    hvisgen_c_Prj_exjet = gen_c_Prj_exjet_binning->CreateHistogram("hvisgen_c_Prj_exjet");
    hvisgen_c_Mrj_exjet = gen_c_Mrj_exjet_binning->CreateHistogram("hvisgen_c_Mrj_exjet");

    hvisgen_c_han_exjet = gen_c_han_exjet_binning->CreateHistogram("hvisgen_c_han_exjet");
    hvisgen_c_sca_exjet = gen_c_sca_exjet_binning->CreateHistogram("hvisgen_c_sca_exjet");
    hvisgen_c_tra_exjet = gen_c_tra_exjet_binning->CreateHistogram("hvisgen_c_tra_exjet");
    hvisgen_c_kjL_exjet = gen_c_kjL_exjet_binning->CreateHistogram("hvisgen_c_kjL_exjet");
    hvisgen_c_rqL_exjet = gen_c_rqL_exjet_binning->CreateHistogram("hvisgen_c_rqL_exjet");
    hvisgen_c_rkP_exjet = gen_c_rkP_exjet_binning->CreateHistogram("hvisgen_c_rkP_exjet");
    hvisgen_c_rkM_exjet = gen_c_rkM_exjet_binning->CreateHistogram("hvisgen_c_rkM_exjet");
    hvisgen_c_nrP_exjet = gen_c_nrP_exjet_binning->CreateHistogram("hvisgen_c_nrP_exjet");
    hvisgen_c_nrM_exjet = gen_c_nrM_exjet_binning->CreateHistogram("hvisgen_c_nrM_exjet");
    hvisgen_c_nkP_exjet = gen_c_nkP_exjet_binning->CreateHistogram("hvisgen_c_nkP_exjet");
    hvisgen_c_nkM_exjet = gen_c_nkM_exjet_binning->CreateHistogram("hvisgen_c_nkM_exjet");

    hvisgen_ll_cHel_exjet = gen_ll_cHel_exjet_binning->CreateHistogram("hvisgen_ll_cHel_exjet");
    hvisgen_ll_cLab_exjet = gen_ll_cLab_exjet_binning->CreateHistogram("hvisgen_ll_cLab_exjet");
    hvisgen_ll_kNorm_exjet = gen_ll_kNorm_exjet_binning->CreateHistogram("hvisgen_ll_kNorm_exjet");
    hvisgen_ll_rNorm_exjet = gen_ll_rNorm_exjet_binning->CreateHistogram("hvisgen_ll_rNorm_exjet");
    hvisgen_llbar_delta_phi_exjet = gen_llbar_delta_phi_exjet_binning->CreateHistogram("hvisgen_llbar_delta_phi_exjet");
    hvisgen_llbar_delta_eta_exjet = gen_llbar_delta_eta_exjet_binning->CreateHistogram("hvisgen_llbar_delta_eta_exjet");

    // ******
    // top_pt
    // ******

    // hvisgen_b1k_top_pt =
    // gen_b1k_top_pt_binning->CreateHistogram("hvisgen_b1k_top_pt");
    // hvisgen_b2k_top_pt =
    // gen_b2k_top_pt_binning->CreateHistogram("hvisgen_b2k_top_pt");
    // hvisgen_b1j_top_pt =
    // gen_b1j_top_pt_binning->CreateHistogram("hvisgen_b1j_top_pt");
    // hvisgen_b2j_top_pt =
    // gen_b2j_top_pt_binning->CreateHistogram("hvisgen_b2j_top_pt");
    // hvisgen_b1r_top_pt =
    // gen_b1r_top_pt_binning->CreateHistogram("hvisgen_b1r_top_pt");
    // hvisgen_b2r_top_pt =
    // gen_b2r_top_pt_binning->CreateHistogram("hvisgen_b2r_top_pt");
    // hvisgen_b1q_top_pt =
    // gen_b1q_top_pt_binning->CreateHistogram("hvisgen_b1q_top_pt");
    // hvisgen_b2q_top_pt =
    // gen_b2q_top_pt_binning->CreateHistogram("hvisgen_b2q_top_pt");
    // hvisgen_b1n_top_pt =
    // gen_b1n_top_pt_binning->CreateHistogram("hvisgen_b1n_top_pt");
    // hvisgen_b2n_top_pt =
    // gen_b2n_top_pt_binning->CreateHistogram("hvisgen_b2n_top_pt");

    // hvisgen_c_kk_top_pt =
    // gen_c_kk_top_pt_binning->CreateHistogram("hvisgen_c_kk_top_pt");
    // hvisgen_c_rr_top_pt =
    // gen_c_rr_top_pt_binning->CreateHistogram("hvisgen_c_rr_top_pt");
    // hvisgen_c_nn_top_pt =
    // gen_c_nn_top_pt_binning->CreateHistogram("hvisgen_c_nn_top_pt");

    // hvisgen_c_rk_top_pt =
    // gen_c_rk_top_pt_binning->CreateHistogram("hvisgen_c_rk_top_pt");
    // hvisgen_c_kr_top_pt =
    // gen_c_kr_top_pt_binning->CreateHistogram("hvisgen_c_kr_top_pt");
    // hvisgen_c_nr_top_pt =
    // gen_c_nr_top_pt_binning->CreateHistogram("hvisgen_c_nr_top_pt");
    // hvisgen_c_rn_top_pt =
    // gen_c_rn_top_pt_binning->CreateHistogram("hvisgen_c_rn_top_pt");
    // hvisgen_c_nk_top_pt =
    // gen_c_nk_top_pt_binning->CreateHistogram("hvisgen_c_nk_top_pt");
    // hvisgen_c_kn_top_pt =
    // gen_c_kn_top_pt_binning->CreateHistogram("hvisgen_c_kn_top_pt");

    // hvisgen_c_Prk_top_pt =
    // gen_c_Prk_top_pt_binning->CreateHistogram("hvisgen_c_Prk_top_pt");
    // hvisgen_c_Mrk_top_pt =
    // gen_c_Mrk_top_pt_binning->CreateHistogram("hvisgen_c_Mrk_top_pt");
    // hvisgen_c_Pnr_top_pt =
    // gen_c_Pnr_top_pt_binning->CreateHistogram("hvisgen_c_Pnr_top_pt");
    // hvisgen_c_Mnr_top_pt =
    // gen_c_Mnr_top_pt_binning->CreateHistogram("hvisgen_c_Mnr_top_pt");
    // hvisgen_c_Pnk_top_pt =
    // gen_c_Pnk_top_pt_binning->CreateHistogram("hvisgen_c_Pnk_top_pt");
    // hvisgen_c_Mnk_top_pt =
    // gen_c_Mnk_top_pt_binning->CreateHistogram("hvisgen_c_Mnk_top_pt");

    // hvisgen_c_kj_top_pt =
    // gen_c_kj_top_pt_binning->CreateHistogram("hvisgen_c_kj_top_pt");
    // hvisgen_c_rq_top_pt =
    // gen_c_rq_top_pt_binning->CreateHistogram("hvisgen_c_rq_top_pt");
    // hvisgen_c_rj_top_pt =
    // gen_c_rj_top_pt_binning->CreateHistogram("hvisgen_c_rj_top_pt");
    // hvisgen_c_jr_top_pt =
    // gen_c_jr_top_pt_binning->CreateHistogram("hvisgen_c_jr_top_pt");

    // hvisgen_c_Prj_top_pt =
    // gen_c_Prj_top_pt_binning->CreateHistogram("hvisgen_c_Prj_top_pt");
    // hvisgen_c_Mrj_top_pt =
    // gen_c_Mrj_top_pt_binning->CreateHistogram("hvisgen_c_Mrj_top_pt");

    // hvisgen_c_han_top_pt =
    // gen_c_han_top_pt_binning->CreateHistogram("hvisgen_c_han_top_pt");
    // hvisgen_c_sca_top_pt =
    // gen_c_sca_top_pt_binning->CreateHistogram("hvisgen_c_sca_top_pt");
    // hvisgen_c_tra_top_pt =
    // gen_c_tra_top_pt_binning->CreateHistogram("hvisgen_c_tra_top_pt");
    // hvisgen_c_kjL_top_pt =
    // gen_c_kjL_top_pt_binning->CreateHistogram("hvisgen_c_kjL_top_pt");
    // hvisgen_c_rqL_top_pt =
    // gen_c_rqL_top_pt_binning->CreateHistogram("hvisgen_c_rqL_top_pt");
    // hvisgen_c_rkP_top_pt =
    // gen_c_rkP_top_pt_binning->CreateHistogram("hvisgen_c_rkP_top_pt");
    // hvisgen_c_rkM_top_pt =
    // gen_c_rkM_top_pt_binning->CreateHistogram("hvisgen_c_rkM_top_pt");
    // hvisgen_c_nrP_top_pt =
    // gen_c_nrP_top_pt_binning->CreateHistogram("hvisgen_c_nrP_top_pt");
    // hvisgen_c_nrM_top_pt =
    // gen_c_nrM_top_pt_binning->CreateHistogram("hvisgen_c_nrM_top_pt");
    // hvisgen_c_nkP_top_pt =
    // gen_c_nkP_top_pt_binning->CreateHistogram("hvisgen_c_nkP_top_pt");
    // hvisgen_c_nkM_top_pt =
    // gen_c_nkM_top_pt_binning->CreateHistogram("hvisgen_c_nkM_top_pt");

    // hvisgen_ll_cHel_top_pt =
    // gen_ll_cHel_top_pt_binning->CreateHistogram("hvisgen_ll_cHel_top_pt");
    // hvisgen_ll_cLab_top_pt =
    // gen_ll_cLab_top_pt_binning->CreateHistogram("hvisgen_ll_cLab_top_pt");
    // hvisgen_ll_kNorm_top_pt =
    // gen_ll_kNorm_top_pt_binning->CreateHistogram("hvisgen_ll_kNorm_top_pt");
    // hvisgen_ll_rNorm_top_pt =
    // gen_ll_rNorm_top_pt_binning->CreateHistogram("hvisgen_ll_rNorm_top_pt");
    // hvisgen_llbar_delta_phi_top_pt =
    // gen_llbar_delta_phi_top_pt_binning->CreateHistogram("hvisgen_llbar_delta_phi_top_pt");
    // hvisgen_llbar_delta_eta_top_pt =
    // gen_llbar_delta_eta_top_pt_binning->CreateHistogram("hvisgen_llbar_delta_eta_top_pt");

    // ********************
    // top_scatteringangle_ttbarframe
    // ********************

    // hvisgen_b1k_top_scatteringangle_ttbarframe =
    // gen_b1k_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_b1k_top_scatteringangle_ttbarframe");
    // hvisgen_b2k_top_scatteringangle_ttbarframe =
    // gen_b2k_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_b2k_top_scatteringangle_ttbarframe");
    // hvisgen_b1j_top_scatteringangle_ttbarframe =
    // gen_b1j_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_b1j_top_scatteringangle_ttbarframe");
    // hvisgen_b2j_top_scatteringangle_ttbarframe =
    // gen_b2j_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_b2j_top_scatteringangle_ttbarframe");
    // hvisgen_b1r_top_scatteringangle_ttbarframe =
    // gen_b1r_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_b1r_top_scatteringangle_ttbarframe");
    // hvisgen_b2r_top_scatteringangle_ttbarframe =
    // gen_b2r_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_b2r_top_scatteringangle_ttbarframe");
    // hvisgen_b1q_top_scatteringangle_ttbarframe =
    // gen_b1q_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_b1q_top_scatteringangle_ttbarframe");
    // hvisgen_b2q_top_scatteringangle_ttbarframe =
    // gen_b2q_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_b2q_top_scatteringangle_ttbarframe");
    // hvisgen_b1n_top_scatteringangle_ttbarframe =
    // gen_b1n_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_b1n_top_scatteringangle_ttbarframe");
    // hvisgen_b2n_top_scatteringangle_ttbarframe =
    // gen_b2n_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_b2n_top_scatteringangle_ttbarframe");

    // hvisgen_c_kk_top_scatteringangle_ttbarframe =
    // gen_c_kk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_kk_top_scatteringangle_ttbarframe");
    // hvisgen_c_rr_top_scatteringangle_ttbarframe =
    // gen_c_rr_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_rr_top_scatteringangle_ttbarframe");
    // hvisgen_c_nn_top_scatteringangle_ttbarframe =
    // gen_c_nn_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_nn_top_scatteringangle_ttbarframe");

    // hvisgen_c_rk_top_scatteringangle_ttbarframe =
    // gen_c_rk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_rk_top_scatteringangle_ttbarframe");
    // hvisgen_c_kr_top_scatteringangle_ttbarframe =
    // gen_c_kr_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_kr_top_scatteringangle_ttbarframe");
    // hvisgen_c_nr_top_scatteringangle_ttbarframe =
    // gen_c_nr_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_nr_top_scatteringangle_ttbarframe");
    // hvisgen_c_rn_top_scatteringangle_ttbarframe =
    // gen_c_rn_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_rn_top_scatteringangle_ttbarframe");
    // hvisgen_c_nk_top_scatteringangle_ttbarframe =
    // gen_c_nk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_nk_top_scatteringangle_ttbarframe");
    // hvisgen_c_kn_top_scatteringangle_ttbarframe =
    // gen_c_kn_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_kn_top_scatteringangle_ttbarframe");

    // hvisgen_c_Prk_top_scatteringangle_ttbarframe =
    // gen_c_Prk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_Prk_top_scatteringangle_ttbarframe");
    // hvisgen_c_Mrk_top_scatteringangle_ttbarframe =
    // gen_c_Mrk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_Mrk_top_scatteringangle_ttbarframe");
    // hvisgen_c_Pnr_top_scatteringangle_ttbarframe =
    // gen_c_Pnr_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_Pnr_top_scatteringangle_ttbarframe");
    // hvisgen_c_Mnr_top_scatteringangle_ttbarframe =
    // gen_c_Mnr_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_Mnr_top_scatteringangle_ttbarframe");
    // hvisgen_c_Pnk_top_scatteringangle_ttbarframe =
    // gen_c_Pnk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_Pnk_top_scatteringangle_ttbarframe");
    // hvisgen_c_Mnk_top_scatteringangle_ttbarframe =
    // gen_c_Mnk_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_Mnk_top_scatteringangle_ttbarframe");

    // hvisgen_c_kj_top_scatteringangle_ttbarframe =
    // gen_c_kj_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_kj_top_scatteringangle_ttbarframe");
    // hvisgen_c_rq_top_scatteringangle_ttbarframe =
    // gen_c_rq_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_rq_top_scatteringangle_ttbarframe");
    // hvisgen_c_rj_top_scatteringangle_ttbarframe =
    // gen_c_rj_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_rj_top_scatteringangle_ttbarframe");
    // hvisgen_c_jr_top_scatteringangle_ttbarframe =
    // gen_c_jr_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_jr_top_scatteringangle_ttbarframe");

    // hvisgen_c_Prj_top_scatteringangle_ttbarframe =
    // gen_c_Prj_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_Prj_top_scatteringangle_ttbarframe");
    // hvisgen_c_Mrj_top_scatteringangle_ttbarframe =
    // gen_c_Mrj_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_Mrj_top_scatteringangle_ttbarframe");

    // hvisgen_c_han_top_scatteringangle_ttbarframe =
    // gen_c_han_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_han_top_scatteringangle_ttbarframe");
    // hvisgen_c_sca_top_scatteringangle_ttbarframe =
    // gen_c_sca_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_sca_top_scatteringangle_ttbarframe");
    // hvisgen_c_tra_top_scatteringangle_ttbarframe =
    // gen_c_tra_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_tra_top_scatteringangle_ttbarframe");
    // hvisgen_c_kjL_top_scatteringangle_ttbarframe =
    // gen_c_kjL_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_kjL_top_scatteringangle_ttbarframe");
    // hvisgen_c_rqL_top_scatteringangle_ttbarframe =
    // gen_c_rqL_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_rqL_top_scatteringangle_ttbarframe");
    // hvisgen_c_rkP_top_scatteringangle_ttbarframe =
    // gen_c_rkP_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_rkP_top_scatteringangle_ttbarframe");
    // hvisgen_c_rkM_top_scatteringangle_ttbarframe =
    // gen_c_rkM_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_rkM_top_scatteringangle_ttbarframe");
    // hvisgen_c_nrP_top_scatteringangle_ttbarframe =
    // gen_c_nrP_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_nrP_top_scatteringangle_ttbarframe");
    // hvisgen_c_nrM_top_scatteringangle_ttbarframe =
    // gen_c_nrM_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_nrM_top_scatteringangle_ttbarframe");
    // hvisgen_c_nkP_top_scatteringangle_ttbarframe =
    // gen_c_nkP_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_nkP_top_scatteringangle_ttbarframe");
    // hvisgen_c_nkM_top_scatteringangle_ttbarframe =
    // gen_c_nkM_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_c_nkM_top_scatteringangle_ttbarframe");

    // hvisgen_ll_cHel_top_scatteringangle_ttbarframe =
    // gen_ll_cHel_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_ll_cHel_top_scatteringangle_ttbarframe");
    // hvisgen_ll_cLab_top_scatteringangle_ttbarframe =
    // gen_ll_cLab_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_ll_cLab_top_scatteringangle_ttbarframe");
    // hvisgen_ll_kNorm_top_scatteringangle_ttbarframe =
    // gen_ll_kNorm_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_ll_kNorm_top_scatteringangle_ttbarframe");
    // hvisgen_ll_rNorm_top_scatteringangle_ttbarframe =
    // gen_ll_rNorm_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_ll_rNorm_top_scatteringangle_ttbarframe");
    // hvisgen_llbar_delta_phi_top_scatteringangle_ttbarframe =
    // gen_llbar_delta_phi_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_llbar_delta_phi_top_scatteringangle_ttbarframe");
    // hvisgen_llbar_delta_eta_top_scatteringangle_ttbarframe =
    // gen_llbar_delta_eta_top_scatteringangle_ttbarframe_binning->CreateHistogram("hvisgen_llbar_delta_eta_top_scatteringangle_ttbarframe");

    // ***************************************
    // ***************************************
    // Response matrices -- filled only for MC
    // ***************************************
    // ***************************************

    // ********************************
    // 1D var migration matrix creation
    // ********************************

    // Kinematic
    hrecoVsgen_top_pt =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_top_pt_binning, reco_top_pt_binning, "hrecoVsgen_top_pt");
    hrecoVsgen_l_pt =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_l_pt_binning, reco_l_pt_binning, "hrecoVsgen_l_pt");
    hrecoVsgen_lbar_pt =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_lbar_pt_binning, reco_lbar_pt_binning, "hrecoVsgen_lbar_pt");
    hrecoVsgen_ttbar_pt =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_ttbar_pt_binning, reco_ttbar_pt_binning, "hrecoVsgen_ttbar_pt");
    hrecoVsgen_ttbar_mass = TUnfoldBinning::CreateHistogramOfMigrations(gen_ttbar_mass_binning, reco_ttbar_mass_binning,
                                                                        "hrecoVsgen_ttbar_mass");
    hrecoVsgen_ttbar_delta_phi = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_ttbar_delta_phi_binning, reco_ttbar_delta_phi_binning, "hrecoVsgen_ttbar_delta_phi");
    hrecoVsgen_ttbar_delta_eta = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_ttbar_delta_eta_binning, reco_ttbar_delta_eta_binning, "hrecoVsgen_ttbar_delta_eta");
    hrecoVsgen_ttbar_rapidity = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_ttbar_rapidity_binning, reco_ttbar_rapidity_binning, "hrecoVsgen_ttbar_rapidity");
    hrecoVsgen_llbar_pt =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_llbar_pt_binning, reco_llbar_pt_binning, "hrecoVsgen_llbar_pt");
    hrecoVsgen_llbar_mass = TUnfoldBinning::CreateHistogramOfMigrations(gen_llbar_mass_binning, reco_llbar_mass_binning,
                                                                        "hrecoVsgen_llbar_mass");

    // Spin corr
    hrecoVsgen_b1k = TUnfoldBinning::CreateHistogramOfMigrations(gen_b1k_binning, reco_b1k_binning, "hrecoVsgen_b1k");
    hrecoVsgen_b2k = TUnfoldBinning::CreateHistogramOfMigrations(gen_b2k_binning, reco_b2k_binning, "hrecoVsgen_b2k");
    hrecoVsgen_b1j = TUnfoldBinning::CreateHistogramOfMigrations(gen_b1j_binning, reco_b1j_binning, "hrecoVsgen_b1j");
    hrecoVsgen_b2j = TUnfoldBinning::CreateHistogramOfMigrations(gen_b2j_binning, reco_b2j_binning, "hrecoVsgen_b2j");
    hrecoVsgen_b1r = TUnfoldBinning::CreateHistogramOfMigrations(gen_b1r_binning, reco_b1r_binning, "hrecoVsgen_b1r");
    hrecoVsgen_b2r = TUnfoldBinning::CreateHistogramOfMigrations(gen_b2r_binning, reco_b2r_binning, "hrecoVsgen_b2r");
    hrecoVsgen_b1q = TUnfoldBinning::CreateHistogramOfMigrations(gen_b1q_binning, reco_b1q_binning, "hrecoVsgen_b1q");
    hrecoVsgen_b2q = TUnfoldBinning::CreateHistogramOfMigrations(gen_b2q_binning, reco_b2q_binning, "hrecoVsgen_b2q");
    hrecoVsgen_b1n = TUnfoldBinning::CreateHistogramOfMigrations(gen_b1n_binning, reco_b1n_binning, "hrecoVsgen_b1n");
    hrecoVsgen_b2n = TUnfoldBinning::CreateHistogramOfMigrations(gen_b2n_binning, reco_b2n_binning, "hrecoVsgen_b2n");

    hrecoVsgen_c_kk =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kk_binning, reco_c_kk_binning, "hrecoVsgen_c_kk");
    hrecoVsgen_c_rr =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rr_binning, reco_c_rr_binning, "hrecoVsgen_c_rr");
    hrecoVsgen_c_nn =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nn_binning, reco_c_nn_binning, "hrecoVsgen_c_nn");

    hrecoVsgen_c_rk =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rk_binning, reco_c_rk_binning, "hrecoVsgen_c_rk");
    hrecoVsgen_c_kr =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kr_binning, reco_c_kr_binning, "hrecoVsgen_c_kr");
    hrecoVsgen_c_nr =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nr_binning, reco_c_nr_binning, "hrecoVsgen_c_nr");
    hrecoVsgen_c_rn =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rn_binning, reco_c_rn_binning, "hrecoVsgen_c_rn");
    hrecoVsgen_c_nk =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nk_binning, reco_c_nk_binning, "hrecoVsgen_c_nk");
    hrecoVsgen_c_kn =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kn_binning, reco_c_kn_binning, "hrecoVsgen_c_kn");

    hrecoVsgen_c_Prk =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Prk_binning, reco_c_Prk_binning, "hrecoVsgen_c_Prk");
    hrecoVsgen_c_Mrk =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Mrk_binning, reco_c_Mrk_binning, "hrecoVsgen_c_Mrk");
    hrecoVsgen_c_Pnr =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Pnr_binning, reco_c_Pnr_binning, "hrecoVsgen_c_Pnr");
    hrecoVsgen_c_Mnr =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Mnr_binning, reco_c_Mnr_binning, "hrecoVsgen_c_Mnr");
    hrecoVsgen_c_Pnk =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Pnk_binning, reco_c_Pnk_binning, "hrecoVsgen_c_Pnk");
    hrecoVsgen_c_Mnk =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Mnk_binning, reco_c_Mnk_binning, "hrecoVsgen_c_Mnk");

    hrecoVsgen_c_kj =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kj_binning, reco_c_kj_binning, "hrecoVsgen_c_kj");
    hrecoVsgen_c_rq =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rq_binning, reco_c_rq_binning, "hrecoVsgen_c_rq");
    hrecoVsgen_c_rj =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rj_binning, reco_c_rj_binning, "hrecoVsgen_c_rj");
    hrecoVsgen_c_jr =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_jr_binning, reco_c_jr_binning, "hrecoVsgen_c_jr");

    hrecoVsgen_c_Prj =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Prj_binning, reco_c_Prj_binning, "hrecoVsgen_c_Prj");
    hrecoVsgen_c_Mrj =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Mrj_binning, reco_c_Mrj_binning, "hrecoVsgen_c_Mrj");

    hrecoVsgen_c_han =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_han_binning, reco_c_han_binning, "hrecoVsgen_c_han");
    hrecoVsgen_c_sca =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_sca_binning, reco_c_sca_binning, "hrecoVsgen_c_sca");
    hrecoVsgen_c_tra =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_tra_binning, reco_c_tra_binning, "hrecoVsgen_c_tra");
    hrecoVsgen_c_kjL =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kjL_binning, reco_c_kjL_binning, "hrecoVsgen_c_kjL");
    hrecoVsgen_c_rqL =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rqL_binning, reco_c_rqL_binning, "hrecoVsgen_c_rqL");
    hrecoVsgen_c_rkP =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rkP_binning, reco_c_rkP_binning, "hrecoVsgen_c_rkP");
    hrecoVsgen_c_rkM =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rkM_binning, reco_c_rkM_binning, "hrecoVsgen_c_rkM");
    hrecoVsgen_c_nrP =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nrP_binning, reco_c_nrP_binning, "hrecoVsgen_c_nrP");
    hrecoVsgen_c_nrM =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nrM_binning, reco_c_nrM_binning, "hrecoVsgen_c_nrM");
    hrecoVsgen_c_nkP =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nkP_binning, reco_c_nkP_binning, "hrecoVsgen_c_nkP");
    hrecoVsgen_c_nkM =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nkM_binning, reco_c_nkM_binning, "hrecoVsgen_c_nkM");

    hrecoVsgen_ll_cHel =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_ll_cHel_binning, reco_ll_cHel_binning, "hrecoVsgen_ll_cHel");
    hrecoVsgen_ll_cLab =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_ll_cLab_binning, reco_ll_cLab_binning, "hrecoVsgen_ll_cLab");
    hrecoVsgen_ll_kNorm =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_ll_kNorm_binning, reco_ll_kNorm_binning, "hrecoVsgen_ll_kNorm");
    hrecoVsgen_ll_rNorm =
        TUnfoldBinning::CreateHistogramOfMigrations(gen_ll_rNorm_binning, reco_ll_rNorm_binning, "hrecoVsgen_ll_rNorm");
    hrecoVsgen_llbar_delta_phi = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_llbar_delta_phi_binning, reco_llbar_delta_phi_binning, "hrecoVsgen_llbar_delta_phi");
    hrecoVsgen_llbar_delta_eta = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_llbar_delta_eta_binning, reco_llbar_delta_eta_binning, "hrecoVsgen_llbar_delta_eta");

    // ********************************
    // 2D var migration matrix creation
    // ********************************
    // Amandeep : 2D var migration matrix creation
    // hrecoVsgen_c_kk_mttbar =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kk_mttbar_binning,
    // reco_c_kk_mttbar_binning, "hrecoVsgen_c_kk_mttbar");

    // ******
    // mttbar
    // ******

    hrecoVsgen_b1k_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_b1k_exjet_binning, reco_b1k_exjet_binning,
                                                                       "hrecoVsgen_b1k_exjet");
    hrecoVsgen_b2k_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_b2k_exjet_binning, reco_b2k_exjet_binning,
                                                                       "hrecoVsgen_b2k_exjet");
    hrecoVsgen_b1j_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_b1j_exjet_binning, reco_b1j_exjet_binning,
                                                                       "hrecoVsgen_b1j_exjet");
    hrecoVsgen_b2j_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_b2j_exjet_binning, reco_b2j_exjet_binning,
                                                                       "hrecoVsgen_b2j_exjet");
    hrecoVsgen_b1r_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_b1r_exjet_binning, reco_b1r_exjet_binning,
                                                                       "hrecoVsgen_b1r_exjet");
    hrecoVsgen_b2r_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_b2r_exjet_binning, reco_b2r_exjet_binning,
                                                                       "hrecoVsgen_b2r_exjet");
    hrecoVsgen_b1q_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_b1q_exjet_binning, reco_b1q_exjet_binning,
                                                                       "hrecoVsgen_b1q_exjet");
    hrecoVsgen_b2q_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_b2q_exjet_binning, reco_b2q_exjet_binning,
                                                                       "hrecoVsgen_b2q_exjet");
    hrecoVsgen_b1n_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_b1n_exjet_binning, reco_b1n_exjet_binning,
                                                                       "hrecoVsgen_b1n_exjet");
    hrecoVsgen_b2n_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_b2n_exjet_binning, reco_b2n_exjet_binning,
                                                                       "hrecoVsgen_b2n_exjet");

    hrecoVsgen_c_kk_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kk_exjet_binning, reco_c_kk_exjet_binning,
                                                                        "hrecoVsgen_c_kk_exjet");
    hrecoVsgen_c_rr_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rr_exjet_binning, reco_c_rr_exjet_binning,
                                                                        "hrecoVsgen_c_rr_exjet");
    hrecoVsgen_c_nn_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nn_exjet_binning, reco_c_nn_exjet_binning,
                                                                        "hrecoVsgen_c_nn_exjet");

    hrecoVsgen_c_rk_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rk_exjet_binning, reco_c_rk_exjet_binning,
                                                                        "hrecoVsgen_c_rk_exjet");
    hrecoVsgen_c_kr_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kr_exjet_binning, reco_c_kr_exjet_binning,
                                                                        "hrecoVsgen_c_kr_exjet");
    hrecoVsgen_c_nr_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nr_exjet_binning, reco_c_nr_exjet_binning,
                                                                        "hrecoVsgen_c_nr_exjet");
    hrecoVsgen_c_rn_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rn_exjet_binning, reco_c_rn_exjet_binning,
                                                                        "hrecoVsgen_c_rn_exjet");
    hrecoVsgen_c_nk_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nk_exjet_binning, reco_c_nk_exjet_binning,
                                                                        "hrecoVsgen_c_nk_exjet");
    hrecoVsgen_c_kn_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kn_exjet_binning, reco_c_kn_exjet_binning,
                                                                        "hrecoVsgen_c_kn_exjet");

    hrecoVsgen_c_Prk_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_Prk_exjet_binning, reco_c_Prk_exjet_binning, "hrecoVsgen_c_Prk_exjet");
    hrecoVsgen_c_Mrk_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_Mrk_exjet_binning, reco_c_Mrk_exjet_binning, "hrecoVsgen_c_Mrk_exjet");
    hrecoVsgen_c_Pnr_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_Pnr_exjet_binning, reco_c_Pnr_exjet_binning, "hrecoVsgen_c_Pnr_exjet");
    hrecoVsgen_c_Mnr_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_Mnr_exjet_binning, reco_c_Mnr_exjet_binning, "hrecoVsgen_c_Mnr_exjet");
    hrecoVsgen_c_Pnk_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_Pnk_exjet_binning, reco_c_Pnk_exjet_binning, "hrecoVsgen_c_Pnk_exjet");
    hrecoVsgen_c_Mnk_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_Mnk_exjet_binning, reco_c_Mnk_exjet_binning, "hrecoVsgen_c_Mnk_exjet");

    hrecoVsgen_c_kj_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kj_exjet_binning, reco_c_kj_exjet_binning,
                                                                        "hrecoVsgen_c_kj_exjet");
    hrecoVsgen_c_rq_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rq_exjet_binning, reco_c_rq_exjet_binning,
                                                                        "hrecoVsgen_c_rq_exjet");
    hrecoVsgen_c_rj_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rj_exjet_binning, reco_c_rj_exjet_binning,
                                                                        "hrecoVsgen_c_rj_exjet");
    hrecoVsgen_c_jr_exjet = TUnfoldBinning::CreateHistogramOfMigrations(gen_c_jr_exjet_binning, reco_c_jr_exjet_binning,
                                                                        "hrecoVsgen_c_jr_exjet");

    hrecoVsgen_c_Prj_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_Prj_exjet_binning, reco_c_Prj_exjet_binning, "hrecoVsgen_c_Prj_exjet");
    hrecoVsgen_c_Mrj_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_Mrj_exjet_binning, reco_c_Mrj_exjet_binning, "hrecoVsgen_c_Mrj_exjet");

    hrecoVsgen_c_han_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_han_exjet_binning, reco_c_han_exjet_binning, "hrecoVsgen_c_han_exjet");
    hrecoVsgen_c_sca_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_sca_exjet_binning, reco_c_sca_exjet_binning, "hrecoVsgen_c_sca_exjet");
    hrecoVsgen_c_tra_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_tra_exjet_binning, reco_c_tra_exjet_binning, "hrecoVsgen_c_tra_exjet");
    hrecoVsgen_c_kjL_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_kjL_exjet_binning, reco_c_kjL_exjet_binning, "hrecoVsgen_c_kjL_exjet");
    hrecoVsgen_c_rqL_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_rqL_exjet_binning, reco_c_rqL_exjet_binning, "hrecoVsgen_c_rqL_exjet");
    hrecoVsgen_c_rkP_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_rkP_exjet_binning, reco_c_rkP_exjet_binning, "hrecoVsgen_c_rkP_exjet");
    hrecoVsgen_c_rkM_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_rkM_exjet_binning, reco_c_rkM_exjet_binning, "hrecoVsgen_c_rkM_exjet");
    hrecoVsgen_c_nrP_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_nrP_exjet_binning, reco_c_nrP_exjet_binning, "hrecoVsgen_c_nrP_exjet");
    hrecoVsgen_c_nrM_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_nrM_exjet_binning, reco_c_nrM_exjet_binning, "hrecoVsgen_c_nrM_exjet");
    hrecoVsgen_c_nkP_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_nkP_exjet_binning, reco_c_nkP_exjet_binning, "hrecoVsgen_c_nkP_exjet");
    hrecoVsgen_c_nkM_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_c_nkM_exjet_binning, reco_c_nkM_exjet_binning, "hrecoVsgen_c_nkM_exjet");

    hrecoVsgen_ll_cHel_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_ll_cHel_exjet_binning, reco_ll_cHel_exjet_binning, "hrecoVsgen_ll_cHel_exjet");
    hrecoVsgen_ll_cLab_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_ll_cLab_exjet_binning, reco_ll_cLab_exjet_binning, "hrecoVsgen_ll_cLab_exjet");
    hrecoVsgen_ll_kNorm_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_ll_kNorm_exjet_binning, reco_ll_kNorm_exjet_binning, "hrecoVsgen_ll_kNorm_exjet");
    hrecoVsgen_ll_rNorm_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_ll_rNorm_exjet_binning, reco_ll_rNorm_exjet_binning, "hrecoVsgen_ll_rNorm_exjet");
    hrecoVsgen_llbar_delta_phi_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_llbar_delta_phi_exjet_binning, reco_llbar_delta_phi_exjet_binning, "hrecoVsgen_llbar_delta_phi_exjet");
    hrecoVsgen_llbar_delta_eta_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        gen_llbar_delta_eta_exjet_binning, reco_llbar_delta_eta_exjet_binning, "hrecoVsgen_llbar_delta_eta_exjet");

    // ******
    // top_pt
    // ******

    // hrecoVsgen_b1k_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_b1k_top_pt_binning,
    // reco_b1k_top_pt_binning, "hrecoVsgen_b1k_top_pt"); hrecoVsgen_b2k_top_pt
    // = TUnfoldBinning::CreateHistogramOfMigrations(gen_b2k_top_pt_binning,
    // reco_b2k_top_pt_binning, "hrecoVsgen_b2k_top_pt"); hrecoVsgen_b1j_top_pt
    // = TUnfoldBinning::CreateHistogramOfMigrations(gen_b1j_top_pt_binning,
    // reco_b1j_top_pt_binning, "hrecoVsgen_b1j_top_pt"); hrecoVsgen_b2j_top_pt
    // = TUnfoldBinning::CreateHistogramOfMigrations(gen_b2j_top_pt_binning,
    // reco_b2j_top_pt_binning, "hrecoVsgen_b2j_top_pt"); hrecoVsgen_b1r_top_pt
    // = TUnfoldBinning::CreateHistogramOfMigrations(gen_b1r_top_pt_binning,
    // reco_b1r_top_pt_binning, "hrecoVsgen_b1r_top_pt"); hrecoVsgen_b2r_top_pt
    // = TUnfoldBinning::CreateHistogramOfMigrations(gen_b2r_top_pt_binning,
    // reco_b2r_top_pt_binning, "hrecoVsgen_b2r_top_pt"); hrecoVsgen_b1q_top_pt
    // = TUnfoldBinning::CreateHistogramOfMigrations(gen_b1q_top_pt_binning,
    // reco_b1q_top_pt_binning, "hrecoVsgen_b1q_top_pt"); hrecoVsgen_b2q_top_pt
    // = TUnfoldBinning::CreateHistogramOfMigrations(gen_b2q_top_pt_binning,
    // reco_b2q_top_pt_binning, "hrecoVsgen_b2q_top_pt"); hrecoVsgen_b1n_top_pt
    // = TUnfoldBinning::CreateHistogramOfMigrations(gen_b1n_top_pt_binning,
    // reco_b1n_top_pt_binning, "hrecoVsgen_b1n_top_pt"); hrecoVsgen_b2n_top_pt
    // = TUnfoldBinning::CreateHistogramOfMigrations(gen_b2n_top_pt_binning,
    // reco_b2n_top_pt_binning, "hrecoVsgen_b2n_top_pt");

    // hrecoVsgen_c_kk_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kk_top_pt_binning,
    // reco_c_kk_top_pt_binning, "hrecoVsgen_c_kk_top_pt");
    // hrecoVsgen_c_rr_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rr_top_pt_binning,
    // reco_c_rr_top_pt_binning, "hrecoVsgen_c_rr_top_pt");
    // hrecoVsgen_c_nn_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nn_top_pt_binning,
    // reco_c_nn_top_pt_binning, "hrecoVsgen_c_nn_top_pt");

    // hrecoVsgen_c_rk_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rk_top_pt_binning,
    // reco_c_rk_top_pt_binning, "hrecoVsgen_c_rk_top_pt");
    // hrecoVsgen_c_kr_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kr_top_pt_binning,
    // reco_c_kr_top_pt_binning, "hrecoVsgen_c_kr_top_pt");
    // hrecoVsgen_c_nr_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nr_top_pt_binning,
    // reco_c_nr_top_pt_binning, "hrecoVsgen_c_nr_top_pt");
    // hrecoVsgen_c_rn_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rn_top_pt_binning,
    // reco_c_rn_top_pt_binning, "hrecoVsgen_c_rn_top_pt");
    // hrecoVsgen_c_nk_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nk_top_pt_binning,
    // reco_c_nk_top_pt_binning, "hrecoVsgen_c_nk_top_pt");
    // hrecoVsgen_c_kn_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kn_top_pt_binning,
    // reco_c_kn_top_pt_binning, "hrecoVsgen_c_kn_top_pt");

    // hrecoVsgen_c_Prk_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Prk_top_pt_binning,
    // reco_c_Prk_top_pt_binning, "hrecoVsgen_c_Prk_top_pt");
    // hrecoVsgen_c_Mrk_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Mrk_top_pt_binning,
    // reco_c_Mrk_top_pt_binning, "hrecoVsgen_c_Mrk_top_pt");
    // hrecoVsgen_c_Pnr_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Pnr_top_pt_binning,
    // reco_c_Pnr_top_pt_binning, "hrecoVsgen_c_Pnr_top_pt");
    // hrecoVsgen_c_Mnr_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Mnr_top_pt_binning,
    // reco_c_Mnr_top_pt_binning, "hrecoVsgen_c_Mnr_top_pt");
    // hrecoVsgen_c_Pnk_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Pnk_top_pt_binning,
    // reco_c_Pnk_top_pt_binning, "hrecoVsgen_c_Pnk_top_pt");
    // hrecoVsgen_c_Mnk_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Mnk_top_pt_binning,
    // reco_c_Mnk_top_pt_binning, "hrecoVsgen_c_Mnk_top_pt");

    // hrecoVsgen_c_kj_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kj_top_pt_binning,
    // reco_c_kj_top_pt_binning, "hrecoVsgen_c_kj_top_pt");
    // hrecoVsgen_c_rq_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rq_top_pt_binning,
    // reco_c_rq_top_pt_binning, "hrecoVsgen_c_rq_top_pt");
    // hrecoVsgen_c_rj_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rj_top_pt_binning,
    // reco_c_rj_top_pt_binning, "hrecoVsgen_c_rj_top_pt");
    // hrecoVsgen_c_jr_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_jr_top_pt_binning,
    // reco_c_jr_top_pt_binning, "hrecoVsgen_c_jr_top_pt");

    // hrecoVsgen_c_Prj_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Prj_top_pt_binning,
    // reco_c_Prj_top_pt_binning, "hrecoVsgen_c_Prj_top_pt");
    // hrecoVsgen_c_Mrj_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Mrj_top_pt_binning,
    // reco_c_Mrj_top_pt_binning, "hrecoVsgen_c_Mrj_top_pt");

    // hrecoVsgen_c_han_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_han_top_pt_binning,
    // reco_c_han_top_pt_binning, "hrecoVsgen_c_han_top_pt");
    // hrecoVsgen_c_sca_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_sca_top_pt_binning,
    // reco_c_sca_top_pt_binning, "hrecoVsgen_c_sca_top_pt");
    // hrecoVsgen_c_tra_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_tra_top_pt_binning,
    // reco_c_tra_top_pt_binning, "hrecoVsgen_c_tra_top_pt");
    // hrecoVsgen_c_kjL_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kjL_top_pt_binning,
    // reco_c_kjL_top_pt_binning, "hrecoVsgen_c_kjL_top_pt");
    // hrecoVsgen_c_rqL_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rqL_top_pt_binning,
    // reco_c_rqL_top_pt_binning, "hrecoVsgen_c_rqL_top_pt");
    // hrecoVsgen_c_rkP_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rkP_top_pt_binning,
    // reco_c_rkP_top_pt_binning, "hrecoVsgen_c_rkP_top_pt");
    // hrecoVsgen_c_rkM_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rkM_top_pt_binning,
    // reco_c_rkM_top_pt_binning, "hrecoVsgen_c_rkM_top_pt");
    // hrecoVsgen_c_nrP_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nrP_top_pt_binning,
    // reco_c_nrP_top_pt_binning, "hrecoVsgen_c_nrP_top_pt");
    // hrecoVsgen_c_nrM_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nrM_top_pt_binning,
    // reco_c_nrM_top_pt_binning, "hrecoVsgen_c_nrM_top_pt");
    // hrecoVsgen_c_nkP_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nkP_top_pt_binning,
    // reco_c_nkP_top_pt_binning, "hrecoVsgen_c_nkP_top_pt");
    // hrecoVsgen_c_nkM_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nkM_top_pt_binning,
    // reco_c_nkM_top_pt_binning, "hrecoVsgen_c_nkM_top_pt");

    // hrecoVsgen_ll_cHel_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_ll_cHel_top_pt_binning,
    // reco_ll_cHel_top_pt_binning, "hrecoVsgen_ll_cHel_top_pt");
    // hrecoVsgen_ll_cLab_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_ll_cLab_top_pt_binning,
    // reco_ll_cLab_top_pt_binning, "hrecoVsgen_ll_cLab_top_pt");
    // hrecoVsgen_ll_kNorm_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_ll_kNorm_top_pt_binning,
    // reco_ll_kNorm_top_pt_binning, "hrecoVsgen_ll_kNorm_top_pt");
    // hrecoVsgen_ll_rNorm_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_ll_rNorm_top_pt_binning,
    // reco_ll_rNorm_top_pt_binning, "hrecoVsgen_ll_rNorm_top_pt");
    // hrecoVsgen_llbar_delta_phi_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_llbar_delta_phi_top_pt_binning,
    // reco_llbar_delta_phi_top_pt_binning,
    // "hrecoVsgen_llbar_delta_phi_top_pt"); hrecoVsgen_llbar_delta_eta_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_llbar_delta_eta_top_pt_binning,
    // reco_llbar_delta_eta_top_pt_binning,
    // "hrecoVsgen_llbar_delta_eta_top_pt");

    // ********************
    // top_scatteringangle_ttbarframe
    // ********************

    // hrecoVsgen_b1k_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_b1k_top_scatteringangle_ttbarframe_binning,
    // reco_b1k_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_b1k_top_scatteringangle_ttbarframe");
    // hrecoVsgen_b2k_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_b2k_top_scatteringangle_ttbarframe_binning,
    // reco_b2k_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_b2k_top_scatteringangle_ttbarframe");
    // hrecoVsgen_b1j_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_b1j_top_scatteringangle_ttbarframe_binning,
    // reco_b1j_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_b1j_top_scatteringangle_ttbarframe");
    // hrecoVsgen_b2j_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_b2j_top_scatteringangle_ttbarframe_binning,
    // reco_b2j_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_b2j_top_scatteringangle_ttbarframe");
    // hrecoVsgen_b1r_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_b1r_top_scatteringangle_ttbarframe_binning,
    // reco_b1r_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_b1r_top_scatteringangle_ttbarframe");
    // hrecoVsgen_b2r_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_b2r_top_scatteringangle_ttbarframe_binning,
    // reco_b2r_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_b2r_top_scatteringangle_ttbarframe");
    // hrecoVsgen_b1q_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_b1q_top_scatteringangle_ttbarframe_binning,
    // reco_b1q_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_b1q_top_scatteringangle_ttbarframe");
    // hrecoVsgen_b2q_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_b2q_top_scatteringangle_ttbarframe_binning,
    // reco_b2q_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_b2q_top_scatteringangle_ttbarframe");
    // hrecoVsgen_b1n_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_b1n_top_scatteringangle_ttbarframe_binning,
    // reco_b1n_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_b1n_top_scatteringangle_ttbarframe");
    // hrecoVsgen_b2n_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_b2n_top_scatteringangle_ttbarframe_binning,
    // reco_b2n_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_b2n_top_scatteringangle_ttbarframe");

    // hrecoVsgen_c_kk_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kk_top_scatteringangle_ttbarframe_binning,
    // reco_c_kk_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_kk_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_rr_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rr_top_scatteringangle_ttbarframe_binning,
    // reco_c_rr_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_rr_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_nn_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nn_top_scatteringangle_ttbarframe_binning,
    // reco_c_nn_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_nn_top_scatteringangle_ttbarframe");

    // hrecoVsgen_c_rk_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rk_top_scatteringangle_ttbarframe_binning,
    // reco_c_rk_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_rk_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_kr_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kr_top_scatteringangle_ttbarframe_binning,
    // reco_c_kr_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_kr_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_nr_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nr_top_scatteringangle_ttbarframe_binning,
    // reco_c_nr_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_nr_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_rn_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rn_top_scatteringangle_ttbarframe_binning,
    // reco_c_rn_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_rn_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_nk_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nk_top_scatteringangle_ttbarframe_binning,
    // reco_c_nk_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_nk_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_kn_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kn_top_scatteringangle_ttbarframe_binning,
    // reco_c_kn_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_kn_top_scatteringangle_ttbarframe");

    // hrecoVsgen_c_Prk_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Prk_top_scatteringangle_ttbarframe_binning,
    // reco_c_Prk_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_Prk_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_Mrk_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Mrk_top_scatteringangle_ttbarframe_binning,
    // reco_c_Mrk_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_Mrk_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_Pnr_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Pnr_top_scatteringangle_ttbarframe_binning,
    // reco_c_Pnr_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_Pnr_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_Mnr_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Mnr_top_scatteringangle_ttbarframe_binning,
    // reco_c_Mnr_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_Mnr_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_Pnk_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Pnk_top_scatteringangle_ttbarframe_binning,
    // reco_c_Pnk_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_Pnk_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_Mnk_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Mnk_top_scatteringangle_ttbarframe_binning,
    // reco_c_Mnk_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_Mnk_top_scatteringangle_ttbarframe");

    // hrecoVsgen_c_kj_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kj_top_scatteringangle_ttbarframe_binning,
    // reco_c_kj_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_kj_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_rq_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rq_top_scatteringangle_ttbarframe_binning,
    // reco_c_rq_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_rq_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_rj_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rj_top_scatteringangle_ttbarframe_binning,
    // reco_c_rj_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_rj_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_jr_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_jr_top_scatteringangle_ttbarframe_binning,
    // reco_c_jr_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_jr_top_scatteringangle_ttbarframe");

    // hrecoVsgen_c_Prj_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Prj_top_scatteringangle_ttbarframe_binning,
    // reco_c_Prj_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_Prj_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_Mrj_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_Mrj_top_scatteringangle_ttbarframe_binning,
    // reco_c_Mrj_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_Mrj_top_scatteringangle_ttbarframe");

    // hrecoVsgen_c_han_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_han_top_scatteringangle_ttbarframe_binning,
    // reco_c_han_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_han_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_sca_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_sca_top_scatteringangle_ttbarframe_binning,
    // reco_c_sca_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_sca_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_tra_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_tra_top_scatteringangle_ttbarframe_binning,
    // reco_c_tra_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_tra_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_kjL_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_kjL_top_scatteringangle_ttbarframe_binning,
    // reco_c_kjL_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_kjL_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_rqL_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rqL_top_scatteringangle_ttbarframe_binning,
    // reco_c_rqL_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_rqL_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_rkP_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rkP_top_scatteringangle_ttbarframe_binning,
    // reco_c_rkP_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_rkP_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_rkM_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_rkM_top_scatteringangle_ttbarframe_binning,
    // reco_c_rkM_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_rkM_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_nrP_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nrP_top_scatteringangle_ttbarframe_binning,
    // reco_c_nrP_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_nrP_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_nrM_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nrM_top_scatteringangle_ttbarframe_binning,
    // reco_c_nrM_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_nrM_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_nkP_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nkP_top_scatteringangle_ttbarframe_binning,
    // reco_c_nkP_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_nkP_top_scatteringangle_ttbarframe");
    // hrecoVsgen_c_nkM_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_c_nkM_top_scatteringangle_ttbarframe_binning,
    // reco_c_nkM_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_c_nkM_top_scatteringangle_ttbarframe");

    // hrecoVsgen_ll_cHel_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_ll_cHel_top_scatteringangle_ttbarframe_binning,
    // reco_ll_cHel_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_ll_cHel_top_scatteringangle_ttbarframe");
    // hrecoVsgen_ll_cLab_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_ll_cLab_top_scatteringangle_ttbarframe_binning,
    // reco_ll_cLab_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_ll_cLab_top_scatteringangle_ttbarframe");
    // hrecoVsgen_ll_kNorm_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_ll_kNorm_top_scatteringangle_ttbarframe_binning,
    // reco_ll_kNorm_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_ll_kNorm_top_scatteringangle_ttbarframe");
    // hrecoVsgen_ll_rNorm_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_ll_rNorm_top_scatteringangle_ttbarframe_binning,
    // reco_ll_rNorm_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_ll_rNorm_top_scatteringangle_ttbarframe");
    // hrecoVsgen_llbar_delta_phi_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_llbar_delta_phi_top_scatteringangle_ttbarframe_binning,
    // reco_llbar_delta_phi_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_llbar_delta_phi_top_scatteringangle_ttbarframe");
    // hrecoVsgen_llbar_delta_eta_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(gen_llbar_delta_eta_top_scatteringangle_ttbarframe_binning,
    // reco_llbar_delta_eta_top_scatteringangle_ttbarframe_binning,
    // "hrecoVsgen_llbar_delta_eta_top_scatteringangle_ttbarframe");

    // ******************
    // 1D var resolutions
    // ******************
    // Resolution plots (reco-gen) for guiding choice of binning

    // Kinematic
    hresolutionbins_top_pt = TUnfoldBinning::CreateHistogramOfMigrations(residual_top_pt_binning, gen_top_pt_binning,
                                                                         "hresolutionbins_top_pt");
    hresolutionbins_l_pt =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_l_pt_binning, gen_l_pt_binning, "hresolutionbins_l_pt");
    hresolutionbins_lbar_pt = TUnfoldBinning::CreateHistogramOfMigrations(residual_lbar_pt_binning, gen_lbar_pt_binning,
                                                                          "hresolutionbins_lbar_pt");
    hresolutionbins_ttbar_pt = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_ttbar_pt_binning, gen_ttbar_pt_binning, "hresolutionbins_ttbar_pt");
    hresolutionbins_ttbar_mass = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_ttbar_mass_binning, gen_ttbar_mass_binning, "hresolutionbins_ttbar_mass");
    hresolutionbins_ttbar_delta_phi = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_ttbar_delta_phi_binning, gen_ttbar_delta_phi_binning, "hresolutionbins_ttbar_delta_phi");
    hresolutionbins_ttbar_delta_eta = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_ttbar_delta_eta_binning, gen_ttbar_delta_eta_binning, "hresolutionbins_ttbar_delta_eta");
    hresolutionbins_ttbar_rapidity = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_ttbar_rapidity_binning, gen_ttbar_rapidity_binning, "hresolutionbins_ttbar_rapidity");
    hresolutionbins_llbar_pt = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_llbar_pt_binning, gen_llbar_pt_binning, "hresolutionbins_llbar_pt");
    hresolutionbins_llbar_mass = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_llbar_mass_binning, gen_llbar_mass_binning, "hresolutionbins_llbar_mass");

    // Spin corr
    hresolutionbins_b1k =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_b1k_binning, gen_b1k_binning, "hresolutionbins_b1k");
    hresolutionbins_b2k =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_b2k_binning, gen_b2k_binning, "hresolutionbins_b2k");
    hresolutionbins_b1j =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_b1j_binning, gen_b1j_binning, "hresolutionbins_b1j");
    hresolutionbins_b2j =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_b2j_binning, gen_b2j_binning, "hresolutionbins_b2j");
    hresolutionbins_b1r =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_b1r_binning, gen_b1r_binning, "hresolutionbins_b1r");
    hresolutionbins_b2r =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_b2r_binning, gen_b2r_binning, "hresolutionbins_b2r");
    hresolutionbins_b1q =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_b1q_binning, gen_b1q_binning, "hresolutionbins_b1q");
    hresolutionbins_b2q =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_b2q_binning, gen_b2q_binning, "hresolutionbins_b2q");
    hresolutionbins_b1n =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_b1n_binning, gen_b1n_binning, "hresolutionbins_b1n");
    hresolutionbins_b2n =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_b2n_binning, gen_b2n_binning, "hresolutionbins_b2n");

    hresolutionbins_c_kk =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_kk_binning, gen_c_kk_binning, "hresolutionbins_c_kk");
    hresolutionbins_c_rr =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rr_binning, gen_c_rr_binning, "hresolutionbins_c_rr");
    hresolutionbins_c_nn =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nn_binning, gen_c_nn_binning, "hresolutionbins_c_nn");

    hresolutionbins_c_rk =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rk_binning, gen_c_rk_binning, "hresolutionbins_c_rk");
    hresolutionbins_c_kr =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_kr_binning, gen_c_kr_binning, "hresolutionbins_c_kr");
    hresolutionbins_c_nr =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nr_binning, gen_c_nr_binning, "hresolutionbins_c_nr");
    hresolutionbins_c_rn =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rn_binning, gen_c_rn_binning, "hresolutionbins_c_rn");
    hresolutionbins_c_nk =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nk_binning, gen_c_nk_binning, "hresolutionbins_c_nk");
    hresolutionbins_c_kn =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_kn_binning, gen_c_kn_binning, "hresolutionbins_c_kn");

    hresolutionbins_c_Prk =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Prk_binning, gen_c_Prk_binning, "hresolutionbins_c_Prk");
    hresolutionbins_c_Mrk =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Mrk_binning, gen_c_Mrk_binning, "hresolutionbins_c_Mrk");
    hresolutionbins_c_Pnr =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Pnr_binning, gen_c_Pnr_binning, "hresolutionbins_c_Pnr");
    hresolutionbins_c_Mnr =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Mnr_binning, gen_c_Mnr_binning, "hresolutionbins_c_Mnr");
    hresolutionbins_c_Pnk =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Pnk_binning, gen_c_Pnk_binning, "hresolutionbins_c_Pnk");
    hresolutionbins_c_Mnk =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Mnk_binning, gen_c_Mnk_binning, "hresolutionbins_c_Mnk");

    hresolutionbins_c_kj =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_kj_binning, gen_c_kj_binning, "hresolutionbins_c_kj");
    hresolutionbins_c_rq =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rq_binning, gen_c_rq_binning, "hresolutionbins_c_rq");
    hresolutionbins_c_rj =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rj_binning, gen_c_rj_binning, "hresolutionbins_c_rj");
    hresolutionbins_c_jr =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_jr_binning, gen_c_jr_binning, "hresolutionbins_c_jr");

    hresolutionbins_c_Prj =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Prj_binning, gen_c_Prj_binning, "hresolutionbins_c_Prj");
    hresolutionbins_c_Mrj =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Mrj_binning, gen_c_Mrj_binning, "hresolutionbins_c_Mrj");

    hresolutionbins_c_han =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_han_binning, gen_c_han_binning, "hresolutionbins_c_han");
    hresolutionbins_c_sca =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_sca_binning, gen_c_sca_binning, "hresolutionbins_c_sca");
    hresolutionbins_c_tra =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_tra_binning, gen_c_tra_binning, "hresolutionbins_c_tra");
    hresolutionbins_c_kjL =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_kjL_binning, gen_c_kjL_binning, "hresolutionbins_c_kjL");
    hresolutionbins_c_rqL =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rqL_binning, gen_c_rqL_binning, "hresolutionbins_c_rqL");
    hresolutionbins_c_rkP =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rkP_binning, gen_c_rkP_binning, "hresolutionbins_c_rkP");
    hresolutionbins_c_rkM =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rkM_binning, gen_c_rkM_binning, "hresolutionbins_c_rkM");
    hresolutionbins_c_nrP =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nrP_binning, gen_c_nrP_binning, "hresolutionbins_c_nrP");
    hresolutionbins_c_nrM =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nrM_binning, gen_c_nrM_binning, "hresolutionbins_c_nrM");
    hresolutionbins_c_nkP =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nkP_binning, gen_c_nkP_binning, "hresolutionbins_c_nkP");
    hresolutionbins_c_nkM =
        TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nkM_binning, gen_c_nkM_binning, "hresolutionbins_c_nkM");

    hresolutionbins_ll_cHel = TUnfoldBinning::CreateHistogramOfMigrations(residual_ll_cHel_binning, gen_ll_cHel_binning,
                                                                          "hresolutionbins_ll_cHel");
    hresolutionbins_ll_cLab = TUnfoldBinning::CreateHistogramOfMigrations(residual_ll_cLab_binning, gen_ll_cLab_binning,
                                                                          "hresolutionbins_ll_cLab");
    hresolutionbins_ll_kNorm = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_ll_kNorm_binning, gen_ll_kNorm_binning, "hresolutionbins_ll_kNorm");
    hresolutionbins_ll_rNorm = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_ll_rNorm_binning, gen_ll_rNorm_binning, "hresolutionbins_ll_rNorm");
    hresolutionbins_llbar_delta_phi = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_llbar_delta_phi_binning, gen_llbar_delta_phi_binning, "hresolutionbins_llbar_delta_phi");
    hresolutionbins_llbar_delta_eta = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_llbar_delta_eta_binning, gen_llbar_delta_eta_binning, "hresolutionbins_llbar_delta_eta");

    // ******************
    // 2D var resolutions
    // ******************
    // Amandeep : 2D resolution histograms
    // hresolutionbins_c_kk_mttbar =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_kk_mttbar_binning,
    // gen_c_kk_mttbar_binning, "hresolutionbins_ll_cHel");

    // ******
    // mttbar
    // ******

    hresolutionbins_b1k_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_b1k_exjet_binning, gen_b1k_exjet_binning, "hresolutionbins_b1k");
    hresolutionbins_b2k_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_b2k_exjet_binning, gen_b2k_exjet_binning, "hresolutionbins_b2k");
    hresolutionbins_b1j_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_b1j_exjet_binning, gen_b1j_exjet_binning, "hresolutionbins_b1j");
    hresolutionbins_b2j_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_b2j_exjet_binning, gen_b2j_exjet_binning, "hresolutionbins_b2j");
    hresolutionbins_b1r_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_b1r_exjet_binning, gen_b1r_exjet_binning, "hresolutionbins_b1r");
    hresolutionbins_b2r_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_b2r_exjet_binning, gen_b2r_exjet_binning, "hresolutionbins_b2r");
    hresolutionbins_b1q_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_b1q_exjet_binning, gen_b1q_exjet_binning, "hresolutionbins_b1q");
    hresolutionbins_b2q_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_b2q_exjet_binning, gen_b2q_exjet_binning, "hresolutionbins_b2q");
    hresolutionbins_b1n_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_b1n_exjet_binning, gen_b1n_exjet_binning, "hresolutionbins_b1n");
    hresolutionbins_b2n_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_b2n_exjet_binning, gen_b2n_exjet_binning, "hresolutionbins_b2n");

    hresolutionbins_c_kk_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_kk_exjet_binning, gen_c_kk_exjet_binning, "hresolutionbins_c_kk");
    hresolutionbins_c_rr_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_rr_exjet_binning, gen_c_rr_exjet_binning, "hresolutionbins_c_rr");
    hresolutionbins_c_nn_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_nn_exjet_binning, gen_c_nn_exjet_binning, "hresolutionbins_c_nn");

    hresolutionbins_c_rk_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_rk_exjet_binning, gen_c_rk_exjet_binning, "hresolutionbins_c_rk");
    hresolutionbins_c_kr_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_kr_exjet_binning, gen_c_kr_exjet_binning, "hresolutionbins_c_kr");
    hresolutionbins_c_nr_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_nr_exjet_binning, gen_c_nr_exjet_binning, "hresolutionbins_c_nr");
    hresolutionbins_c_rn_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_rn_exjet_binning, gen_c_rn_exjet_binning, "hresolutionbins_c_rn");
    hresolutionbins_c_nk_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_nk_exjet_binning, gen_c_nk_exjet_binning, "hresolutionbins_c_nk");
    hresolutionbins_c_kn_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_kn_exjet_binning, gen_c_kn_exjet_binning, "hresolutionbins_c_kn");

    hresolutionbins_c_Prk_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_Prk_exjet_binning, gen_c_Prk_exjet_binning, "hresolutionbins_c_Prk");
    hresolutionbins_c_Mrk_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_Mrk_exjet_binning, gen_c_Mrk_exjet_binning, "hresolutionbins_c_Mrk");
    hresolutionbins_c_Pnr_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_Pnr_exjet_binning, gen_c_Pnr_exjet_binning, "hresolutionbins_c_Pnr");
    hresolutionbins_c_Mnr_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_Mnr_exjet_binning, gen_c_Mnr_exjet_binning, "hresolutionbins_c_Mnr");
    hresolutionbins_c_Pnk_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_Pnk_exjet_binning, gen_c_Pnk_exjet_binning, "hresolutionbins_c_Pnk");
    hresolutionbins_c_Mnk_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_Mnk_exjet_binning, gen_c_Mnk_exjet_binning, "hresolutionbins_c_Mnk");

    hresolutionbins_c_kj_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_kj_exjet_binning, gen_c_kj_exjet_binning, "hresolutionbins_c_kj");
    hresolutionbins_c_rq_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_rq_exjet_binning, gen_c_rq_exjet_binning, "hresolutionbins_c_rq");
    hresolutionbins_c_rj_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_rj_exjet_binning, gen_c_rj_exjet_binning, "hresolutionbins_c_rj");
    hresolutionbins_c_jr_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_jr_exjet_binning, gen_c_jr_exjet_binning, "hresolutionbins_c_jr");

    hresolutionbins_c_Prj_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_Prj_exjet_binning, gen_c_Prj_exjet_binning, "hresolutionbins_c_Prj");
    hresolutionbins_c_Mrj_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_Mrj_exjet_binning, gen_c_Mrj_exjet_binning, "hresolutionbins_c_Mrj");

    hresolutionbins_c_han_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_han_exjet_binning, gen_c_han_exjet_binning, "hresolutionbins_c_han");
    hresolutionbins_c_sca_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_sca_exjet_binning, gen_c_sca_exjet_binning, "hresolutionbins_c_sca");
    hresolutionbins_c_tra_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_tra_exjet_binning, gen_c_tra_exjet_binning, "hresolutionbins_c_tra");
    hresolutionbins_c_kjL_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_kjL_exjet_binning, gen_c_kjL_exjet_binning, "hresolutionbins_c_kjL");
    hresolutionbins_c_rqL_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_rqL_exjet_binning, gen_c_rqL_exjet_binning, "hresolutionbins_c_rqL");
    hresolutionbins_c_rkP_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_rkP_exjet_binning, gen_c_rkP_exjet_binning, "hresolutionbins_c_rkP");
    hresolutionbins_c_rkM_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_rkM_exjet_binning, gen_c_rkM_exjet_binning, "hresolutionbins_c_rkM");
    hresolutionbins_c_nrP_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_nrP_exjet_binning, gen_c_nrP_exjet_binning, "hresolutionbins_c_nrP");
    hresolutionbins_c_nrM_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_nrM_exjet_binning, gen_c_nrM_exjet_binning, "hresolutionbins_c_nrM");
    hresolutionbins_c_nkP_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_nkP_exjet_binning, gen_c_nkP_exjet_binning, "hresolutionbins_c_nkP");
    hresolutionbins_c_nkM_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_c_nkM_exjet_binning, gen_c_nkM_exjet_binning, "hresolutionbins_c_nkM");

    hresolutionbins_ll_cHel_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_ll_cHel_exjet_binning, gen_ll_cHel_exjet_binning, "hresolutionbins_ll_cHel");
    hresolutionbins_ll_cLab_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_ll_cLab_exjet_binning, gen_ll_cLab_exjet_binning, "hresolutionbins_ll_cLab");
    hresolutionbins_ll_kNorm_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_ll_kNorm_exjet_binning, gen_ll_kNorm_exjet_binning, "hresolutionbins_ll_kNorm");
    hresolutionbins_ll_rNorm_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_ll_rNorm_exjet_binning, gen_ll_rNorm_exjet_binning, "hresolutionbins_ll_rNorm");
    hresolutionbins_llbar_delta_phi_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_llbar_delta_phi_exjet_binning, gen_llbar_delta_phi_exjet_binning, "hresolutionbins_llbar_delta_phi");
    hresolutionbins_llbar_delta_eta_exjet = TUnfoldBinning::CreateHistogramOfMigrations(
        residual_llbar_delta_eta_exjet_binning, gen_llbar_delta_eta_exjet_binning, "hresolutionbins_llbar_delta_eta");

    // ******
    // top_pt
    // ******

    // hresolutionbins_b1k_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b1k_top_pt_binning,
    // gen_b1k_top_pt_binning, "hresolutionbins_b1k");
    // hresolutionbins_b2k_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b2k_top_pt_binning,
    // gen_b2k_top_pt_binning, "hresolutionbins_b2k");
    // hresolutionbins_b1j_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b1j_top_pt_binning,
    // gen_b1j_top_pt_binning, "hresolutionbins_b1j");
    // hresolutionbins_b2j_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b2j_top_pt_binning,
    // gen_b2j_top_pt_binning, "hresolutionbins_b2j");
    // hresolutionbins_b1r_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b1r_top_pt_binning,
    // gen_b1r_top_pt_binning, "hresolutionbins_b1r");
    // hresolutionbins_b2r_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b2r_top_pt_binning,
    // gen_b2r_top_pt_binning, "hresolutionbins_b2r");
    // hresolutionbins_b1q_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b1q_top_pt_binning,
    // gen_b1q_top_pt_binning, "hresolutionbins_b1q");
    // hresolutionbins_b2q_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b2q_top_pt_binning,
    // gen_b2q_top_pt_binning, "hresolutionbins_b2q");
    // hresolutionbins_b1n_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b1n_top_pt_binning,
    // gen_b1n_top_pt_binning, "hresolutionbins_b1n");
    // hresolutionbins_b2n_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b2n_top_pt_binning,
    // gen_b2n_top_pt_binning, "hresolutionbins_b2n");

    // hresolutionbins_c_kk_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_kk_top_pt_binning,
    // gen_c_kk_top_pt_binning, "hresolutionbins_c_kk");
    // hresolutionbins_c_rr_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rr_top_pt_binning,
    // gen_c_rr_top_pt_binning, "hresolutionbins_c_rr");
    // hresolutionbins_c_nn_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nn_top_pt_binning,
    // gen_c_nn_top_pt_binning, "hresolutionbins_c_nn");

    // hresolutionbins_c_rk_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rk_top_pt_binning,
    // gen_c_rk_top_pt_binning, "hresolutionbins_c_rk");
    // hresolutionbins_c_kr_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_kr_top_pt_binning,
    // gen_c_kr_top_pt_binning, "hresolutionbins_c_kr");
    // hresolutionbins_c_nr_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nr_top_pt_binning,
    // gen_c_nr_top_pt_binning, "hresolutionbins_c_nr");
    // hresolutionbins_c_rn_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rn_top_pt_binning,
    // gen_c_rn_top_pt_binning, "hresolutionbins_c_rn");
    // hresolutionbins_c_nk_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nk_top_pt_binning,
    // gen_c_nk_top_pt_binning, "hresolutionbins_c_nk");
    // hresolutionbins_c_kn_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_kn_top_pt_binning,
    // gen_c_kn_top_pt_binning, "hresolutionbins_c_kn");

    // hresolutionbins_c_Prk_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Prk_top_pt_binning,
    // gen_c_Prk_top_pt_binning, "hresolutionbins_c_Prk");
    // hresolutionbins_c_Mrk_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Mrk_top_pt_binning,
    // gen_c_Mrk_top_pt_binning, "hresolutionbins_c_Mrk");
    // hresolutionbins_c_Pnr_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Pnr_top_pt_binning,
    // gen_c_Pnr_top_pt_binning, "hresolutionbins_c_Pnr");
    // hresolutionbins_c_Mnr_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Mnr_top_pt_binning,
    // gen_c_Mnr_top_pt_binning, "hresolutionbins_c_Mnr");
    // hresolutionbins_c_Pnk_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Pnk_top_pt_binning,
    // gen_c_Pnk_top_pt_binning, "hresolutionbins_c_Pnk");
    // hresolutionbins_c_Mnk_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Mnk_top_pt_binning,
    // gen_c_Mnk_top_pt_binning, "hresolutionbins_c_Mnk");

    // hresolutionbins_c_kj_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_kj_top_pt_binning,
    // gen_c_kj_top_pt_binning, "hresolutionbins_c_kj");
    // hresolutionbins_c_rq_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rq_top_pt_binning,
    // gen_c_rq_top_pt_binning, "hresolutionbins_c_rq");
    // hresolutionbins_c_rj_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rj_top_pt_binning,
    // gen_c_rj_top_pt_binning, "hresolutionbins_c_rj");
    // hresolutionbins_c_jr_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_jr_top_pt_binning,
    // gen_c_jr_top_pt_binning, "hresolutionbins_c_jr");

    // hresolutionbins_c_Prj_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Prj_top_pt_binning,
    // gen_c_Prj_top_pt_binning, "hresolutionbins_c_Prj");
    // hresolutionbins_c_Mrj_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Mrj_top_pt_binning,
    // gen_c_Mrj_top_pt_binning, "hresolutionbins_c_Mrj");

    // hresolutionbins_c_han_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_han_top_pt_binning,
    // gen_c_han_top_pt_binning, "hresolutionbins_c_han");
    // hresolutionbins_c_sca_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_sca_top_pt_binning,
    // gen_c_sca_top_pt_binning, "hresolutionbins_c_sca");
    // hresolutionbins_c_tra_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_tra_top_pt_binning,
    // gen_c_tra_top_pt_binning, "hresolutionbins_c_tra");
    // hresolutionbins_c_kjL_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_kjL_top_pt_binning,
    // gen_c_kjL_top_pt_binning, "hresolutionbins_c_kjL");
    // hresolutionbins_c_rqL_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rqL_top_pt_binning,
    // gen_c_rqL_top_pt_binning, "hresolutionbins_c_rqL");
    // hresolutionbins_c_rkP_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rkP_top_pt_binning,
    // gen_c_rkP_top_pt_binning, "hresolutionbins_c_rkP");
    // hresolutionbins_c_rkM_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rkM_top_pt_binning,
    // gen_c_rkM_top_pt_binning, "hresolutionbins_c_rkM");
    // hresolutionbins_c_nrP_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nrP_top_pt_binning,
    // gen_c_nrP_top_pt_binning, "hresolutionbins_c_nrP");
    // hresolutionbins_c_nrM_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nrM_top_pt_binning,
    // gen_c_nrM_top_pt_binning, "hresolutionbins_c_nrM");
    // hresolutionbins_c_nkP_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nkP_top_pt_binning,
    // gen_c_nkP_top_pt_binning, "hresolutionbins_c_nkP");
    // hresolutionbins_c_nkM_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nkM_top_pt_binning,
    // gen_c_nkM_top_pt_binning, "hresolutionbins_c_nkM");

    // hresolutionbins_ll_cHel_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_ll_cHel_top_pt_binning,
    // gen_ll_cHel_top_pt_binning, "hresolutionbins_ll_cHel");
    // hresolutionbins_ll_cLab_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_ll_cLab_top_pt_binning,
    // gen_ll_cLab_top_pt_binning, "hresolutionbins_ll_cLab");
    // hresolutionbins_ll_kNorm_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_ll_kNorm_top_pt_binning,
    // gen_ll_kNorm_top_pt_binning, "hresolutionbins_ll_kNorm");
    // hresolutionbins_ll_rNorm_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_ll_rNorm_top_pt_binning,
    // gen_ll_rNorm_top_pt_binning, "hresolutionbins_ll_rNorm");
    // hresolutionbins_llbar_delta_phi_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_llbar_delta_phi_top_pt_binning,
    // gen_llbar_delta_phi_top_pt_binning, "hresolutionbins_llbar_delta_phi");
    // hresolutionbins_llbar_delta_eta_top_pt =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_llbar_delta_eta_top_pt_binning,
    // gen_llbar_delta_eta_top_pt_binning, "hresolutionbins_llbar_delta_eta");

    // ********************
    // top_scatteringangle_ttbarframe
    // ********************

    // hresolutionbins_b1k_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b1k_top_scatteringangle_ttbarframe_binning,
    // gen_b1k_top_scatteringangle_ttbarframe_binning, "hresolutionbins_b1k");
    // hresolutionbins_b2k_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b2k_top_scatteringangle_ttbarframe_binning,
    // gen_b2k_top_scatteringangle_ttbarframe_binning, "hresolutionbins_b2k");
    // hresolutionbins_b1j_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b1j_top_scatteringangle_ttbarframe_binning,
    // gen_b1j_top_scatteringangle_ttbarframe_binning, "hresolutionbins_b1j");
    // hresolutionbins_b2j_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b2j_top_scatteringangle_ttbarframe_binning,
    // gen_b2j_top_scatteringangle_ttbarframe_binning, "hresolutionbins_b2j");
    // hresolutionbins_b1r_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b1r_top_scatteringangle_ttbarframe_binning,
    // gen_b1r_top_scatteringangle_ttbarframe_binning, "hresolutionbins_b1r");
    // hresolutionbins_b2r_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b2r_top_scatteringangle_ttbarframe_binning,
    // gen_b2r_top_scatteringangle_ttbarframe_binning, "hresolutionbins_b2r");
    // hresolutionbins_b1q_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b1q_top_scatteringangle_ttbarframe_binning,
    // gen_b1q_top_scatteringangle_ttbarframe_binning, "hresolutionbins_b1q");
    // hresolutionbins_b2q_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b2q_top_scatteringangle_ttbarframe_binning,
    // gen_b2q_top_scatteringangle_ttbarframe_binning, "hresolutionbins_b2q");
    // hresolutionbins_b1n_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b1n_top_scatteringangle_ttbarframe_binning,
    // gen_b1n_top_scatteringangle_ttbarframe_binning, "hresolutionbins_b1n");
    // hresolutionbins_b2n_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_b2n_top_scatteringangle_ttbarframe_binning,
    // gen_b2n_top_scatteringangle_ttbarframe_binning, "hresolutionbins_b2n");

    // hresolutionbins_c_kk_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_kk_top_scatteringangle_ttbarframe_binning,
    // gen_c_kk_top_scatteringangle_ttbarframe_binning, "hresolutionbins_c_kk");
    // hresolutionbins_c_rr_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rr_top_scatteringangle_ttbarframe_binning,
    // gen_c_rr_top_scatteringangle_ttbarframe_binning, "hresolutionbins_c_rr");
    // hresolutionbins_c_nn_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nn_top_scatteringangle_ttbarframe_binning,
    // gen_c_nn_top_scatteringangle_ttbarframe_binning, "hresolutionbins_c_nn");

    // hresolutionbins_c_rk_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rk_top_scatteringangle_ttbarframe_binning,
    // gen_c_rk_top_scatteringangle_ttbarframe_binning, "hresolutionbins_c_rk");
    // hresolutionbins_c_kr_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_kr_top_scatteringangle_ttbarframe_binning,
    // gen_c_kr_top_scatteringangle_ttbarframe_binning, "hresolutionbins_c_kr");
    // hresolutionbins_c_nr_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nr_top_scatteringangle_ttbarframe_binning,
    // gen_c_nr_top_scatteringangle_ttbarframe_binning, "hresolutionbins_c_nr");
    // hresolutionbins_c_rn_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rn_top_scatteringangle_ttbarframe_binning,
    // gen_c_rn_top_scatteringangle_ttbarframe_binning, "hresolutionbins_c_rn");
    // hresolutionbins_c_nk_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nk_top_scatteringangle_ttbarframe_binning,
    // gen_c_nk_top_scatteringangle_ttbarframe_binning, "hresolutionbins_c_nk");
    // hresolutionbins_c_kn_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_kn_top_scatteringangle_ttbarframe_binning,
    // gen_c_kn_top_scatteringangle_ttbarframe_binning, "hresolutionbins_c_kn");

    // hresolutionbins_c_Prk_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Prk_top_scatteringangle_ttbarframe_binning,
    // gen_c_Prk_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_Prk");
    // hresolutionbins_c_Mrk_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Mrk_top_scatteringangle_ttbarframe_binning,
    // gen_c_Mrk_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_Mrk");
    // hresolutionbins_c_Pnr_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Pnr_top_scatteringangle_ttbarframe_binning,
    // gen_c_Pnr_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_Pnr");
    // hresolutionbins_c_Mnr_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Mnr_top_scatteringangle_ttbarframe_binning,
    // gen_c_Mnr_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_Mnr");
    // hresolutionbins_c_Pnk_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Pnk_top_scatteringangle_ttbarframe_binning,
    // gen_c_Pnk_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_Pnk");
    // hresolutionbins_c_Mnk_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Mnk_top_scatteringangle_ttbarframe_binning,
    // gen_c_Mnk_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_Mnk");

    // hresolutionbins_c_kj_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_kj_top_scatteringangle_ttbarframe_binning,
    // gen_c_kj_top_scatteringangle_ttbarframe_binning, "hresolutionbins_c_kj");
    // hresolutionbins_c_rq_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rq_top_scatteringangle_ttbarframe_binning,
    // gen_c_rq_top_scatteringangle_ttbarframe_binning, "hresolutionbins_c_rq");
    // hresolutionbins_c_rj_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rj_top_scatteringangle_ttbarframe_binning,
    // gen_c_rj_top_scatteringangle_ttbarframe_binning, "hresolutionbins_c_rj");
    // hresolutionbins_c_jr_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_jr_top_scatteringangle_ttbarframe_binning,
    // gen_c_jr_top_scatteringangle_ttbarframe_binning, "hresolutionbins_c_jr");

    // hresolutionbins_c_Prj_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Prj_top_scatteringangle_ttbarframe_binning,
    // gen_c_Prj_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_Prj");
    // hresolutionbins_c_Mrj_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_Mrj_top_scatteringangle_ttbarframe_binning,
    // gen_c_Mrj_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_Mrj");

    // hresolutionbins_c_han_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_han_top_scatteringangle_ttbarframe_binning,
    // gen_c_han_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_han");
    // hresolutionbins_c_sca_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_sca_top_scatteringangle_ttbarframe_binning,
    // gen_c_sca_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_sca");
    // hresolutionbins_c_tra_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_tra_top_scatteringangle_ttbarframe_binning,
    // gen_c_tra_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_tra");
    // hresolutionbins_c_kjL_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_kjL_top_scatteringangle_ttbarframe_binning,
    // gen_c_kjL_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_kjL");
    // hresolutionbins_c_rqL_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rqL_top_scatteringangle_ttbarframe_binning,
    // gen_c_rqL_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_rqL");
    // hresolutionbins_c_rkP_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rkP_top_scatteringangle_ttbarframe_binning,
    // gen_c_rkP_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_rkP");
    // hresolutionbins_c_rkM_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_rkM_top_scatteringangle_ttbarframe_binning,
    // gen_c_rkM_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_rkM");
    // hresolutionbins_c_nrP_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nrP_top_scatteringangle_ttbarframe_binning,
    // gen_c_nrP_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_nrP");
    // hresolutionbins_c_nrM_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nrM_top_scatteringangle_ttbarframe_binning,
    // gen_c_nrM_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_nrM");
    // hresolutionbins_c_nkP_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nkP_top_scatteringangle_ttbarframe_binning,
    // gen_c_nkP_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_nkP");
    // hresolutionbins_c_nkM_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_c_nkM_top_scatteringangle_ttbarframe_binning,
    // gen_c_nkM_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_c_nkM");

    // hresolutionbins_ll_cHel_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_ll_cHel_top_scatteringangle_ttbarframe_binning,
    // gen_ll_cHel_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_ll_cHel");
    // hresolutionbins_ll_cLab_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_ll_cLab_top_scatteringangle_ttbarframe_binning,
    // gen_ll_cLab_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_ll_cLab");
    // hresolutionbins_ll_kNorm_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_ll_kNorm_top_scatteringangle_ttbarframe_binning,
    // gen_ll_kNorm_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_ll_kNorm");
    // hresolutionbins_ll_rNorm_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_ll_rNorm_top_scatteringangle_ttbarframe_binning,
    // gen_ll_rNorm_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_ll_rNorm");
    // hresolutionbins_llbar_delta_phi_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_llbar_delta_phi_top_scatteringangle_ttbarframe_binning,
    // gen_llbar_delta_phi_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_llbar_delta_phi");
    // hresolutionbins_llbar_delta_eta_top_scatteringangle_ttbarframe =
    // TUnfoldBinning::CreateHistogramOfMigrations(residual_llbar_delta_eta_top_scatteringangle_ttbarframe_binning,
    // gen_llbar_delta_eta_top_scatteringangle_ttbarframe_binning,
    // "hresolutionbins_llbar_delta_eta");

    std::cout << "Histograms are booked !!" << std::endl;
}

/**
 * @brief Creates the TUnfoldBinning object for the residual, reco - gen,
 * histograms.
 *
 * @param numResidualBins Number of bins to bin the residuals in.
 * @param lowerBinEdge The minimum value of the variable you are computing the
 * residual in.
 * @param upperBinEdge The maximum value of the variable you are computing the
 * residual in.
 * @return TUnfoldBinning* The binning for the residual axis of the histogram.
 */
TUnfoldBinning *makeTUnfoldHisto::CreateResidualBinning(int numResidualBins, double lowerBinEdge, double upperBinEdge) {
    TUnfoldBinning *resolutionBinning = new TUnfoldBinning("resolution");
    TUnfoldBinning *resolutionDistribution = resolutionBinning->AddBinning("residual");

    double residualBinWidth = (upperBinEdge - lowerBinEdge) / numResidualBins;
    double residual_bins[numResidualBins + 1] = {};
    for (int i = 0; i < numResidualBins + 1; ++i) {
        residual_bins[i] = i * residualBinWidth + lowerBinEdge;
    }

    resolutionDistribution->AddAxis("(reco-gen)", numResidualBins, residual_bins, false, false);
    return resolutionDistribution;
}

/**
 * @brief Create the .xml binning files to be read later on for unfolding
 multi-differential distributions.
 *
 * @param varname The name of the variable that is being unfolded.
 * @param axisNames The vector of axis names that you will be adding to this
 distribution to be unfolded. Typically these are the different differential
 variables you will be unfolding in.
 * @param gen_nbin Vector of number of bins on each axis at the generator (true)
 level.
 * @param reco_nbin Vector of number of bins on each axis at the reconstructed
 (detector) level.
 * @param gen_bins The vector containing the binning for each axis at the
 generator level. Each element of the vector contains an array of doubles
 detailing the edges for each bin for that axis.
 * @param reco_bins The vector containing the binning for each axis at the
 reconstruction level. Each element of the vector contains an array of doubles
 detailing the edges for each bin for that axis.
 * @param numBinSubdivisions How many bins each bin was subdivided into. This
 helps with numerically stability of unfolding outputs and empirically has been
 shown to increase our resolution on measured spin correlation coefficients.
 * @return std::pair<TUnfoldBinning *, TUnfoldBinning *> The pair of binnings
 for this variable's distributions. First element of the pair is the generator
 distribution's binning. Second element of the pair is the detected
 distribution's binning.
 */
std::pair<TUnfoldBinning *, TUnfoldBinning *> makeTUnfoldHisto::CreateXMLBinFiles(
    std::string varname, std::vector<std::string> axisNames, std::vector<int> gen_nbin, std::vector<int> reco_nbin,
    std::vector<double *> gen_bins, std::vector<double *> reco_bins, int numBinSubdivisions) {
    // TODO: Add bin factor function generation here. See for example:
    // TF2 *userFunc=new TF2("userfunc","1./x+0.2*y^2",ptBinsCoarse[0],
    // ptBinsCoarse[NBIN_PT_COARSE],
    // etaBinsCoarse[0],etaBinsCoarse[NBIN_ETA_COARSE]);
    // signalBinning->SetBinFactorFunction(1.0,userFunc);

    // TODO: After above is implemented XML binning files should be created
    // independently of filling of histograms,
    //       in case bin factor functions need to be updated

    std::ofstream xmlOut("binning/" + varname + "_binning.xml");
    TUnfoldBinning *generatorBinning = new TUnfoldBinning("generator");
    TUnfoldBinning *genDistribution = generatorBinning->AddBinning("ttbargen");
    for (int i = 0; i < axisNames.size(); ++i) {
        genDistribution->AddAxis(axisNames[i].c_str(), gen_nbin[i], gen_bins[i], false, false);
    }

    TUnfoldBinning *detectorBinning = new TUnfoldBinning("detector");
    TUnfoldBinning *detectorDistribution = detectorBinning->AddBinning("ttbarreco");
    for (int i = 0; i < axisNames.size(); ++i) {
        detectorDistribution->AddAxis(axisNames[i].c_str(), reco_nbin[i], reco_bins[i], false, false);
    }

    TUnfoldBinningXML::ExportXML(*detectorBinning, xmlOut, kTRUE, kFALSE);
    TUnfoldBinningXML::ExportXML(*generatorBinning, xmlOut, kFALSE, kTRUE);
    TUnfoldBinningXML::WriteDTD("binning/tunfoldbinning.dtd");
    xmlOut.close();

    std::ofstream numBinSubdivisionsOutput("binning/" + varname + "_rebin.txt");
    numBinSubdivisionsOutput << numBinSubdivisions << std::endl;
    numBinSubdivisionsOutput.close();
    return std::make_pair(genDistribution, detectorDistribution);
}

void makeTUnfoldHisto::CreateFineBinningArray(const int &input_nbin, const int &output_nbin, double input_bins[],
                                              double output_bins[]) {
    int n_sub_bins = output_nbin / input_nbin;
    for (int i = 0; i < input_nbin; ++i) {
        for (int j = 0; j < n_sub_bins; ++j) {
            output_bins[i * n_sub_bins + j] =
                input_bins[i] + (input_bins[i + 1] - input_bins[i]) * double(j) / double(n_sub_bins);
        }
    }
    output_bins[input_nbin * n_sub_bins] = input_bins[input_nbin];
}

// Terminates the process
// Saves the new histograms

void makeTUnfoldHisto::Terminate(double scalingval) {
    fout->cd();
    fout->Write();
    // Maybe scale histograms here
    TIter next(fout->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next())) {
        TClass *c1 = gROOT->GetClass(key->GetClassName());
        if (!c1->InheritsFrom("TH1")) continue;
        TH1 *h1 = (TH1 *)key->ReadObj();
        std::string histoname = h1->GetName();
        if (histoname.find("reco") != std::string::npos) {
            h1->Scale(scalingval);
        } else
            continue;
    }
    fout->Close();
    fout = 0;
    std::cout << "+++++++++++++++++++++++++" << std::endl;
}
