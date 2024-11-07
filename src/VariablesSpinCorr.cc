#include <iostream>
#include <cstdlib>

#include <Math/VectorUtil.h>
#include "TLorentzVector.h"
#include <TVector3.h>

#include "VariablesSpinCorr.h"
#include "analysisStructs.h"
#include "../../common/include/classes.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/KinematicReconstructionSolution.h"
#include "VariablesPhiTT.h"
//#include "AnalyzerBaseClass.h" //FIXME: rename to AnalyzerBase.h
#include "../../common/include/AnalysisBase.h"
#include "../../common/include/mt2_bisect.h"

#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <TH1.h>
#include <TFile.h>
// ----------------------------------------- Methods for VariablesSpinCorr -----------------------------------------------------------



VariablesSpinCorr::VariablesSpinCorr():
VariablesBase(),
top_pt_(VariableFloat(name_top_pt_)),
top_phi_(VariableFloat(name_top_phi_)),
top_rapidity_(VariableFloat(name_top_rapidity_)),
//top_arapidity_(VariableFloat(name_top_arapidity_)),
top_eta_(VariableFloat(name_top_eta_)),
top_mass_(VariableFloat(name_top_mass_)),

tbar_pt_(VariableFloat(name_tbar_pt_)),
tbar_phi_(VariableFloat(name_tbar_phi_)),
tbar_rapidity_(VariableFloat(name_tbar_rapidity_)),
//tbar_arapidity_(VariableFloat(name_tbar_arapidity_)),
tbar_eta_(VariableFloat(name_tbar_eta_)),
tbar_mass_(VariableFloat(name_tbar_mass_)),

l_pt_(VariableFloat(name_l_pt_)),
l_eta_(VariableFloat(name_l_eta_)),
l_phi_(VariableFloat(name_l_phi_)),
l_mass_(VariableFloat(name_l_mass_)),
l_pdgid_(VariableFloat(name_l_pdgid_)),

lbar_pt_(VariableFloat(name_lbar_pt_)),
lbar_eta_(VariableFloat(name_lbar_eta_)),
lbar_phi_(VariableFloat(name_lbar_phi_)),
lbar_mass_(VariableFloat(name_lbar_mass_)),
lbar_pdgid_(VariableFloat(name_lbar_pdgid_)),

//e_pt_(VariableFloat(name_e_pt_)),
//e_eta_(VariableFloat(name_e_eta_)),
//e_phi_(VariableFloat(name_e_phi_)),
//e_mass_(VariableFloat(name_e_mass_)),

//ebar_pt_(VariableFloat(name_ebar_pt_)),
//ebar_eta_(VariableFloat(name_ebar_eta_)),
//ebar_phi_(VariableFloat(name_ebar_phi_)),
//ebar_mass_(VariableFloat(name_ebar_mass_)),

//mu_pt_(VariableFloat(name_mu_pt_)),
//mu_eta_(VariableFloat(name_mu_eta_)),
//mu_phi_(VariableFloat(name_mu_phi_)),
//mu_mass_(VariableFloat(name_mu_mass_)),

//mubar_pt_(VariableFloat(name_mubar_pt_)),
//mubar_eta_(VariableFloat(name_mubar_eta_)),
//mubar_phi_(VariableFloat(name_mubar_phi_)),
//mubar_mass_(VariableFloat(name_mubar_mass_)),

b_pt_(VariableFloat(name_b_pt_)),
b_eta_(VariableFloat(name_b_eta_)),
b_phi_(VariableFloat(name_b_phi_)),
b_mass_(VariableFloat(name_b_mass_)),

bbar_pt_(VariableFloat(name_bbar_pt_)),
bbar_eta_(VariableFloat(name_bbar_eta_)),
bbar_phi_(VariableFloat(name_bbar_phi_)),
bbar_mass_(VariableFloat(name_bbar_mass_)),

nu_pt_(VariableFloat(name_nu_pt_)),
nu_eta_(VariableFloat(name_nu_eta_)),
nu_phi_(VariableFloat(name_nu_phi_)),
nu_mass_(VariableFloat(name_nu_mass_)),

nubar_pt_(VariableFloat(name_nubar_pt_)),
nubar_eta_(VariableFloat(name_nubar_eta_)),
nubar_phi_(VariableFloat(name_nubar_phi_)),
nubar_mass_(VariableFloat(name_nubar_mass_)),

met_pt_(VariableFloat(name_met_pt_)),
met_phi_(VariableFloat(name_met_phi_)),
met_mass_(VariableFloat(name_met_mass_)),

ttbar_pt_(VariableFloat(name_ttbar_pt_)),
ttbar_phi_(VariableFloat(name_ttbar_phi_)),
ttbar_rapidity_(VariableFloat(name_ttbar_rapidity_)),
//ttbar_arapidity_(VariableFloat(name_ttbar_arapidity_)),
ttbar_eta_(VariableFloat(name_ttbar_eta_)),
ttbar_delta_phi_(VariableFloat(name_ttbar_delta_phi_)),
ttbar_delta_eta_(VariableFloat(name_ttbar_delta_eta_)),
ttbar_delta_rapidity_(VariableFloat(name_ttbar_delta_rapidity_)),
ttbar_mass_(VariableFloat(name_ttbar_mass_)),

llbar_pt_(VariableFloat(name_llbar_pt_)),
llbar_phi_(VariableFloat(name_llbar_phi_)),
llbar_rapidity_(VariableFloat(name_llbar_rapidity_)),
//llbar_arapidity_(VariableFloat(name_llbar_arapidity_)),
llbar_delta_phi_(VariableFloat(name_llbar_delta_phi_)),
llbar_delta_eta_(VariableFloat(name_llbar_delta_eta_)),
llbar_delta_rapidity_(VariableFloat(name_llbar_delta_rapidity_)),
llbar_mass_(VariableFloat(name_llbar_mass_)),

//bbbar_pt_(VariableFloat(name_bbbar_pt_)),
//bbbar_phi_(VariableFloat(name_bbbar_phi_)),
//bbbar_rapidity_(VariableFloat(name_bbbar_rapidity_)),
//bbbar_arapidity_(VariableFloat(name_bbbar_arapidity_)),
//bbbar_delta_phi_(VariableFloat(name_bbbar_delta_phi_)),
//bbbar_delta_eta_(VariableFloat(name_bbbar_delta_eta_)),
//bbbar_delta_rapidity_(VariableFloat(name_bbbar_delta_rapidity_)),
//bbbar_mass_(VariableFloat(name_bbbar_mass_)),

//nunubar_pt_(VariableFloat(name_nunubar_pt_)),
//nunubar_phi_(VariableFloat(name_nunubar_phi_)),
//nunubar_rapidity_(VariableFloat(name_nunubar_rapidity_)),
//nunubar_arapidity_(VariableFloat(name_nunubar_arapidity_)),
//nunubar_delta_phi_(VariableFloat(name_nunubar_delta_phi_)),
//nunubar_delta_eta_(VariableFloat(name_nunubar_delta_eta_)),
//nunubar_delta_rapidity_(VariableFloat(name_nunubar_delta_rapidity_)),
//nunubar_mass_(VariableFloat(name_nunubar_mass_)),

MT2_(VariableFloat(name_MT2_)),

all_mass_(VariableFloat(name_all_mass_)),
r_mass_(VariableFloat(name_r_mass_)),
jet_multiplicity_(VariableFloat(name_jet_multiplicity_)),
bjet_multiplicity_(VariableFloat(name_bjet_multiplicity_)),
x1_(VariableFloat(name_x1_)),
x2_(VariableFloat(name_x2_)),

//ExtraJet variables
n_extraJets_iso08_(VariableFloat(name_n_extraJets_iso08_)),
//ttbar frame variables
top_scatteringangle_ttbarframe_(VariableFloat(name_top_scatteringangle_ttbarframe_)),
//top_pt_ttbarframe_(VariableFloat(name_top_pt_ttbarframe_)),
//top_phi_ttbarframe_(VariableFloat(name_top_phi_ttbarframe_)),
//top_rapidity_ttbarframe_(VariableFloat(name_top_rapidity_ttbarframe_)),
//top_arapidity_ttbarframe_(VariableFloat(name_top_arapidity_ttbarframe_)),
//tbar_pt_ttbarframe_(VariableFloat(name_tbar_pt_ttbarframe_)),
//tbar_phi_ttbarframe_(VariableFloat(name_tbar_phi_ttbarframe_)),
//tbar_rapidity_ttbarframe_(VariableFloat(name_tbar_rapidity_ttbarframe_)),
//tbar_arapidity_ttbarframe_(VariableFloat(name_tbar_arapidity_ttbarframe_)),

//top and antitop frame variables
//l_pt_antitopframe_(VariableFloat(name_l_pt_antitopframe_)),
//l_eta_antitopframe_(VariableFloat(name_l_eta_antitopframe_)),
//l_phi_antitopframe_(VariableFloat(name_l_phi_antitopframe_)),
//lbar_pt_topframe_(VariableFloat(name_lbar_pt_topframe_)),
//lbar_eta_topframe_(VariableFloat(name_lbar_eta_topframe_)),
//lbar_phi_topframe_(VariableFloat(name_lbar_phi_topframe_)),

//Spin correlation variables
b1k_(VariableFloat(name_b1k_)),
b2k_(VariableFloat(name_b2k_)),
b1j_(VariableFloat(name_b1j_)),
b2j_(VariableFloat(name_b2j_)),
b1r_(VariableFloat(name_b1r_)),
b2r_(VariableFloat(name_b2r_)),
b1q_(VariableFloat(name_b1q_)),
b2q_(VariableFloat(name_b2q_)),
b1n_(VariableFloat(name_b1n_)),
b2n_(VariableFloat(name_b2n_)),

//b_Pkk_(VariableFloat(name_b_Pkk_)),
//b_Mkk_(VariableFloat(name_b_Mkk_)),
//b_Pjj_(VariableFloat(name_b_Pjj_)),
//b_Mjj_(VariableFloat(name_b_Mjj_)),
//b_Prr_(VariableFloat(name_b_Prr_)),
//b_Mrr_(VariableFloat(name_b_Mrr_)),
//b_Pqq_(VariableFloat(name_b_Pqq_)),
//b_Mqq_(VariableFloat(name_b_Mqq_)),
//b_Pnn_(VariableFloat(name_b_Pnn_)),
//b_Mnn_(VariableFloat(name_b_Mnn_)),

c_kk_(VariableFloat(name_c_kk_)),
c_rr_(VariableFloat(name_c_rr_)),
c_nn_(VariableFloat(name_c_nn_)),
c_kj_(VariableFloat(name_c_kj_)),
c_rq_(VariableFloat(name_c_rq_)),

c_rk_(VariableFloat(name_c_rk_)),
c_kr_(VariableFloat(name_c_kr_)),
c_nr_(VariableFloat(name_c_nr_)),
c_rn_(VariableFloat(name_c_rn_)),
c_nk_(VariableFloat(name_c_nk_)),
c_kn_(VariableFloat(name_c_kn_)),

c_rj_(VariableFloat(name_c_rj_)),
c_jr_(VariableFloat(name_c_jr_)),
//c_qk_(VariableFloat(name_c_qk_)),
//c_kq_(VariableFloat(name_c_kq_)),

//c_qj_(VariableFloat(name_c_qj_)),
//c_jq_(VariableFloat(name_c_jq_)),
//c_nq_(VariableFloat(name_c_nq_)),
//c_qn_(VariableFloat(name_c_qn_)),
//c_nj_(VariableFloat(name_c_nj_)),
//c_jn_(VariableFloat(name_c_jn_)),

//c_Prk_(VariableFloat(name_c_Prk_)),
//c_Mrk_(VariableFloat(name_c_Mrk_)),
//c_Pnr_(VariableFloat(name_c_Pnr_)),
//c_Mnr_(VariableFloat(name_c_Mnr_)),
//c_Pnk_(VariableFloat(name_c_Pnk_)),
//c_Mnk_(VariableFloat(name_c_Mnk_)),

//c_Pqj_(VariableFloat(name_c_Pqj_)),
//c_Mqj_(VariableFloat(name_c_Mqj_)),
//c_Pnq_(VariableFloat(name_c_Pnq_)),
//c_Mnq_(VariableFloat(name_c_Mnq_)),
//c_Pnj_(VariableFloat(name_c_Pnj_)),
//c_Mnj_(VariableFloat(name_c_Mnj_)),

//c_han_(VariableFloat(name_c_han_)),
//c_sca_(VariableFloat(name_c_sca_)),
//c_tra_(VariableFloat(name_c_tra_)),
//c_lij_(VariableFloat(name_c_lij_)),
//c_liq_(VariableFloat(name_c_liq_)),

//c_rkP_(VariableFloat(name_c_rkP_)),
//c_rkM_(VariableFloat(name_c_rkM_)),
//c_nrP_(VariableFloat(name_c_nrP_)),
//c_nrM_(VariableFloat(name_c_nrM_)),
//c_nkP_(VariableFloat(name_c_nkP_)),
//c_nkM_(VariableFloat(name_c_nkM_)),

//c_qjP_(VariableFloat(name_c_qjP_)),
//c_qjM_(VariableFloat(name_c_qjM_)),
//c_nqP_(VariableFloat(name_c_nqP_)),
//c_nqM_(VariableFloat(name_c_nqM_)),
//c_njP_(VariableFloat(name_c_njP_)),
//c_njM_(VariableFloat(name_c_njM_)),

ll_cHel_(VariableFloat(name_ll_cHel_)),
ll_cLab_(VariableFloat(name_ll_cLab_)),
ll_kNorm_(VariableFloat(name_ll_kNorm_)),
ll_rNorm_(VariableFloat(name_ll_rNorm_)),

gen_top_pt_(VariableFloat(name_gen_top_pt_)),
gen_top_phi_(VariableFloat(name_gen_top_phi_)),
gen_top_rapidity_(VariableFloat(name_gen_top_rapidity_)),
//gen_top_arapidity_(VariableFloat(name_gen_top_arapidity_)),
gen_top_eta_(VariableFloat(name_gen_top_eta_)),
gen_top_mass_(VariableFloat(name_gen_top_mass_)),

gen_tbar_pt_(VariableFloat(name_gen_tbar_pt_)),
gen_tbar_phi_(VariableFloat(name_gen_tbar_phi_)),
gen_tbar_rapidity_(VariableFloat(name_gen_tbar_rapidity_)),
//gen_tbar_arapidity_(VariableFloat(name_gen_tbar_arapidity_)),
gen_tbar_eta_(VariableFloat(name_gen_tbar_eta_)),
gen_tbar_mass_(VariableFloat(name_gen_tbar_mass_)),

gen_l_pt_(VariableFloat(name_gen_l_pt_)),
gen_l_eta_(VariableFloat(name_gen_l_eta_)),
gen_l_phi_(VariableFloat(name_gen_l_phi_)),
gen_l_mass_(VariableFloat(name_gen_l_mass_)),
gen_l_pdgid_(VariableFloat(name_gen_l_pdgid_)),

gen_lbar_pt_(VariableFloat(name_gen_lbar_pt_)),
gen_lbar_eta_(VariableFloat(name_gen_lbar_eta_)),
gen_lbar_phi_(VariableFloat(name_gen_lbar_phi_)),
gen_lbar_mass_(VariableFloat(name_gen_lbar_mass_)),
gen_lbar_pdgid_(VariableFloat(name_gen_lbar_pdgid_)),

//gen_e_pt_(VariableFloat(name_gen_e_pt_)),
//gen_e_eta_(VariableFloat(name_gen_e_eta_)),
//gen_e_phi_(VariableFloat(name_gen_e_phi_)),
//gen_e_mass_(VariableFloat(name_gen_e_mass_)),

//gen_ebar_pt_(VariableFloat(name_gen_ebar_pt_)),
//gen_ebar_eta_(VariableFloat(name_gen_ebar_eta_)),
//gen_ebar_phi_(VariableFloat(name_gen_ebar_phi_)),
//gen_ebar_mass_(VariableFloat(name_gen_ebar_mass_)),

//gen_mu_pt_(VariableFloat(name_gen_mu_pt_)),
//gen_mu_eta_(VariableFloat(name_gen_mu_eta_)),
//gen_mu_phi_(VariableFloat(name_gen_mu_phi_)),
//gen_mu_mass_(VariableFloat(name_gen_mu_mass_)),

//gen_mubar_pt_(VariableFloat(name_gen_mubar_pt_)),
//gen_mubar_eta_(VariableFloat(name_gen_mubar_eta_)),
//gen_mubar_phi_(VariableFloat(name_gen_mubar_phi_)),
//gen_mubar_mass_(VariableFloat(name_gen_mubar_mass_)),

gen_b_pt_(VariableFloat(name_gen_b_pt_)),
gen_b_eta_(VariableFloat(name_gen_b_eta_)),
gen_b_phi_(VariableFloat(name_gen_b_phi_)),
gen_b_mass_(VariableFloat(name_gen_b_mass_)),

gen_bbar_pt_(VariableFloat(name_gen_bbar_pt_)),
gen_bbar_eta_(VariableFloat(name_gen_bbar_eta_)),
gen_bbar_phi_(VariableFloat(name_gen_bbar_phi_)),
gen_bbar_mass_(VariableFloat(name_gen_bbar_mass_)),

gen_nu_pt_(VariableFloat(name_gen_nu_pt_)),
gen_nu_eta_(VariableFloat(name_gen_nu_eta_)),
gen_nu_phi_(VariableFloat(name_gen_nu_phi_)),
gen_nu_mass_(VariableFloat(name_gen_nu_mass_)),

gen_nubar_pt_(VariableFloat(name_gen_nubar_pt_)),
gen_nubar_eta_(VariableFloat(name_gen_nubar_eta_)),
gen_nubar_phi_(VariableFloat(name_gen_nubar_phi_)),
gen_nubar_mass_(VariableFloat(name_gen_nubar_mass_)),

gen_ttbar_pt_(VariableFloat(name_gen_ttbar_pt_)),
gen_ttbar_phi_(VariableFloat(name_gen_ttbar_phi_)),
gen_ttbar_rapidity_(VariableFloat(name_gen_ttbar_rapidity_)),
//gen_ttbar_arapidity_(VariableFloat(name_gen_ttbar_arapidity_)),
gen_ttbar_eta_(VariableFloat(name_gen_ttbar_eta_)),
gen_ttbar_delta_phi_(VariableFloat(name_gen_ttbar_delta_phi_)),
gen_ttbar_delta_eta_(VariableFloat(name_gen_ttbar_delta_eta_)),
gen_ttbar_delta_rapidity_(VariableFloat(name_gen_ttbar_delta_rapidity_)),
gen_ttbar_mass_(VariableFloat(name_gen_ttbar_mass_)),

gen_llbar_pt_(VariableFloat(name_gen_llbar_pt_)),
gen_llbar_phi_(VariableFloat(name_gen_llbar_phi_)),
gen_llbar_rapidity_(VariableFloat(name_gen_llbar_rapidity_)),
//gen_llbar_arapidity_(VariableFloat(name_gen_llbar_arapidity_)),
gen_llbar_delta_phi_(VariableFloat(name_gen_llbar_delta_phi_)),
gen_llbar_delta_eta_(VariableFloat(name_gen_llbar_delta_eta_)),
gen_llbar_delta_rapidity_(VariableFloat(name_gen_llbar_delta_rapidity_)),
gen_llbar_mass_(VariableFloat(name_gen_llbar_mass_)),

//gen_bbbar_pt_(VariableFloat(name_gen_bbbar_pt_)),
//gen_bbbar_phi_(VariableFloat(name_gen_bbbar_phi_)),
//gen_bbbar_rapidity_(VariableFloat(name_gen_bbbar_rapidity_)),
//gen_bbbar_arapidity_(VariableFloat(name_gen_bbbar_arapidity_)),
//gen_bbbar_delta_phi_(VariableFloat(name_gen_bbbar_delta_phi_)),
//gen_bbbar_delta_eta_(VariableFloat(name_gen_bbbar_delta_eta_)),
//gen_bbbar_delta_rapidity_(VariableFloat(name_gen_bbbar_delta_rapidity_)),
//gen_bbbar_mass_(VariableFloat(name_gen_bbbar_mass_)),

//gen_nunubar_pt_(VariableFloat(name_gen_nunubar_pt_)),
//gen_nunubar_phi_(VariableFloat(name_gen_nunubar_phi_)),
//gen_nunubar_rapidity_(VariableFloat(name_gen_nunubar_rapidity_)),
//gen_nunubar_arapidity_(VariableFloat(name_gen_nunubar_arapidity_)),
//gen_nunubar_delta_phi_(VariableFloat(name_gen_nunubar_delta_phi_)),
//gen_nunubar_delta_eta_(VariableFloat(name_gen_nunubar_delta_eta_)),
//gen_nunubar_delta_rapidity_(VariableFloat(name_gen_nunubar_delta_rapidity_)),
//gen_nunubar_mass_(VariableFloat(name_gen_nunubar_mass_)),

gen_MT2_(VariableFloat(name_gen_MT2_)),

gen_all_mass_(VariableFloat(name_gen_all_mass_)),
gen_r_mass_(VariableFloat(name_gen_r_mass_)),
gen_jet_multiplicity_(VariableFloat(name_gen_jet_multiplicity_)),
gen_x1_(VariableFloat(name_gen_x1_)),
gen_x2_(VariableFloat(name_gen_x2_)),
gen_production_mode_(VariableFloat(name_gen_production_mode_)),

//ExtraJet variables
gen_n_extraJets_iso08_(VariableFloat(name_gen_n_extraJets_iso08_)),

//ttbar frame variables
gen_top_scatteringangle_ttbarframe_(VariableFloat(name_gen_top_scatteringangle_ttbarframe_)),
//gen_top_pt_ttbarframe_(VariableFloat(name_gen_top_pt_ttbarframe_)),
//gen_top_phi_ttbarframe_(VariableFloat(name_gen_top_phi_ttbarframe_)),
//gen_top_rapidity_ttbarframe_(VariableFloat(name_gen_top_rapidity_ttbarframe_)),
//gen_top_arapidity_ttbarframe_(VariableFloat(name_gen_top_arapidity_ttbarframe_)),
//gen_tbar_pt_ttbarframe_(VariableFloat(name_gen_tbar_pt_ttbarframe_)),
//gen_tbar_phi_ttbarframe_(VariableFloat(name_gen_tbar_phi_ttbarframe_)),
//gen_tbar_rapidity_ttbarframe_(VariableFloat(name_gen_tbar_rapidity_ttbarframe_)),
//gen_tbar_arapidity_ttbarframe_(VariableFloat(name_gen_tbar_arapidity_ttbarframe_)),

//top and antitop frame variables
//gen_l_pt_antitopframe_(VariableFloat(name_gen_l_pt_antitopframe_)),
//gen_l_eta_antitopframe_(VariableFloat(name_gen_l_eta_antitopframe_)),
//gen_l_phi_antitopframe_(VariableFloat(name_gen_l_phi_antitopframe_)),
//gen_l_mass_(VariableFloat(name_gen_l_mass_)),
//gen_lbar_pt_topframe_(VariableFloat(name_gen_lbar_pt_topframe_)),
//gen_lbar_eta_topframe_(VariableFloat(name_gen_lbar_eta_topframe_)),
//gen_lbar_phi_topframe_(VariableFloat(name_gen_lbar_phi_topframe_)),
//gen_lbar_mass_(VariableFloat(name_gen_lbar_mass_)),

//Spin correlation variables
gen_b1k_(VariableFloat(name_gen_b1k_)),
gen_b2k_(VariableFloat(name_gen_b2k_)),
gen_b1j_(VariableFloat(name_gen_b1j_)),
gen_b2j_(VariableFloat(name_gen_b2j_)),
gen_b1r_(VariableFloat(name_gen_b1r_)),
gen_b2r_(VariableFloat(name_gen_b2r_)),
gen_b1q_(VariableFloat(name_gen_b1q_)),
gen_b2q_(VariableFloat(name_gen_b2q_)),
gen_b1n_(VariableFloat(name_gen_b1n_)),
gen_b2n_(VariableFloat(name_gen_b2n_)),

//gen_b_Pkk_(VariableFloat(name_gen_b_Pkk_)),
//gen_b_Mkk_(VariableFloat(name_gen_b_Mkk_)),
//gen_b_Pjj_(VariableFloat(name_gen_b_Pjj_)),
//gen_b_Mjj_(VariableFloat(name_gen_b_Mjj_)),
//gen_b_Prr_(VariableFloat(name_gen_b_Prr_)),
//gen_b_Mrr_(VariableFloat(name_gen_b_Mrr_)),
//gen_b_Pqq_(VariableFloat(name_gen_b_Pqq_)),
//gen_b_Mqq_(VariableFloat(name_gen_b_Mqq_)),
//gen_b_Pnn_(VariableFloat(name_gen_b_Pnn_)),
//gen_b_Mnn_(VariableFloat(name_gen_b_Mnn_)),

gen_c_kk_(VariableFloat(name_gen_c_kk_)),
gen_c_rr_(VariableFloat(name_gen_c_rr_)),
gen_c_nn_(VariableFloat(name_gen_c_nn_)),
gen_c_kj_(VariableFloat(name_gen_c_kj_)),
gen_c_rq_(VariableFloat(name_gen_c_rq_)),

gen_c_rk_(VariableFloat(name_gen_c_rk_)),
gen_c_kr_(VariableFloat(name_gen_c_kr_)),
gen_c_nr_(VariableFloat(name_gen_c_nr_)),
gen_c_rn_(VariableFloat(name_gen_c_rn_)),
gen_c_nk_(VariableFloat(name_gen_c_nk_)),
gen_c_kn_(VariableFloat(name_gen_c_kn_)),

gen_c_rj_(VariableFloat(name_gen_c_rj_)),
gen_c_jr_(VariableFloat(name_gen_c_jr_)),
//gen_c_qk_(VariableFloat(name_gen_c_qk_)),
//gen_c_kq_(VariableFloat(name_gen_c_kq_)),

//gen_c_qj_(VariableFloat(name_gen_c_qj_)),
//gen_c_jq_(VariableFloat(name_gen_c_jq_)),
//gen_c_nq_(VariableFloat(name_gen_c_nq_)),
//gen_c_qn_(VariableFloat(name_gen_c_qn_)),
//gen_c_nj_(VariableFloat(name_gen_c_nj_)),
//gen_c_jn_(VariableFloat(name_gen_c_jn_)),

//gen_c_Prk_(VariableFloat(name_gen_c_Prk_)),
//gen_c_Mrk_(VariableFloat(name_gen_c_Mrk_)),
//gen_c_Pnr_(VariableFloat(name_gen_c_Pnr_)),
//gen_c_Mnr_(VariableFloat(name_gen_c_Mnr_)),
//gen_c_Pnk_(VariableFloat(name_gen_c_Pnk_)),
//gen_c_Mnk_(VariableFloat(name_gen_c_Mnk_)),

//gen_c_Pqj_(VariableFloat(name_gen_c_Pqj_)),
//gen_c_Mqj_(VariableFloat(name_gen_c_Mqj_)),
//gen_c_Pnq_(VariableFloat(name_gen_c_Pnq_)),
//gen_c_Mnq_(VariableFloat(name_gen_c_Mnq_)),
//gen_c_Pnj_(VariableFloat(name_gen_c_Pnj_)),
//gen_c_Mnj_(VariableFloat(name_gen_c_Mnj_)),

//gen_c_han_(VariableFloat(name_gen_c_han_)),
//gen_c_sca_(VariableFloat(name_gen_c_sca_)),
//gen_c_tra_(VariableFloat(name_gen_c_tra_)),
//gen_c_lij_(VariableFloat(name_gen_c_lij_)),
//gen_c_liq_(VariableFloat(name_gen_c_liq_)),

//gen_c_rkP_(VariableFloat(name_gen_c_rkP_)),
//gen_c_rkM_(VariableFloat(name_gen_c_rkM_)),
//gen_c_nrP_(VariableFloat(name_gen_c_nrP_)),
//gen_c_nrM_(VariableFloat(name_gen_c_nrM_)),
//gen_c_nkP_(VariableFloat(name_gen_c_nkP_)),
//gen_c_nkM_(VariableFloat(name_gen_c_nkM_)),

//gen_c_qjP_(VariableFloat(name_gen_c_qjP_)),
//gen_c_qjM_(VariableFloat(name_gen_c_qjM_)),
//gen_c_nqP_(VariableFloat(name_gen_c_nqP_)),
//gen_c_nqM_(VariableFloat(name_gen_c_nqM_)),
//gen_c_njP_(VariableFloat(name_gen_c_njP_)),
//gen_c_njM_(VariableFloat(name_gen_c_njM_)),

gen_ll_cHel_(VariableFloat(name_gen_ll_cHel_)),
gen_ll_cLab_(VariableFloat(name_gen_ll_cLab_)),
gen_ll_kNorm_(VariableFloat(name_gen_ll_kNorm_)),
gen_ll_rNorm_(VariableFloat(name_gen_ll_rNorm_)),

entry_(VariableInt(name_entry_)),
isTopGen_(VariableInt(name_isTopGen_)),
TopDecayMode_(VariableInt(name_TopDecayMode_)),
isKinReco_(VariableInt(name_isKinReco_)),
eventWeight_(VariableFloat(name_eventWeight_)),
trueLevelWeight_(VariableFloat(name_trueLevelWeight_))
{}



VariablesSpinCorr::VariablesSpinCorr(const EventMetadata& eventMetadata,
                               const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                                          const TopGenObjects& topGenObjects,
                                                          const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                                                          const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                                                          const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                                                          const double& weight):
VariablesBase(weight),
top_pt_(VariableFloat(name_top_pt_)),
top_phi_(VariableFloat(name_top_phi_)),
top_rapidity_(VariableFloat(name_top_rapidity_)),
//top_arapidity_(VariableFloat(name_top_arapidity_)),
top_eta_(VariableFloat(name_top_eta_)),
top_mass_(VariableFloat(name_top_mass_)),

tbar_pt_(VariableFloat(name_tbar_pt_)),
tbar_phi_(VariableFloat(name_tbar_phi_)),
tbar_rapidity_(VariableFloat(name_tbar_rapidity_)),
//tbar_arapidity_(VariableFloat(name_tbar_arapidity_)),
tbar_eta_(VariableFloat(name_tbar_eta_)),
tbar_mass_(VariableFloat(name_tbar_mass_)),

l_pt_(VariableFloat(name_l_pt_)),
l_eta_(VariableFloat(name_l_eta_)),
l_phi_(VariableFloat(name_l_phi_)),
l_mass_(VariableFloat(name_l_mass_)),
l_pdgid_(VariableFloat(name_l_pdgid_)),

lbar_pt_(VariableFloat(name_lbar_pt_)),
lbar_eta_(VariableFloat(name_lbar_eta_)),
lbar_phi_(VariableFloat(name_lbar_phi_)),
lbar_mass_(VariableFloat(name_lbar_mass_)),
lbar_pdgid_(VariableFloat(name_lbar_pdgid_)),

//e_pt_(VariableFloat(name_e_pt_)),
//e_eta_(VariableFloat(name_e_eta_)),
//e_phi_(VariableFloat(name_e_phi_)),
//e_mass_(VariableFloat(name_e_mass_)),

//ebar_pt_(VariableFloat(name_ebar_pt_)),
//ebar_eta_(VariableFloat(name_ebar_eta_)),
//ebar_phi_(VariableFloat(name_ebar_phi_)),
//ebar_mass_(VariableFloat(name_ebar_mass_)),

//mu_pt_(VariableFloat(name_mu_pt_)),
//mu_eta_(VariableFloat(name_mu_eta_)),
//mu_phi_(VariableFloat(name_mu_phi_)),
//mu_mass_(VariableFloat(name_mu_mass_)),

//mubar_pt_(VariableFloat(name_mubar_pt_)),
//mubar_eta_(VariableFloat(name_mubar_eta_)),
//mubar_phi_(VariableFloat(name_mubar_phi_)),
//mubar_mass_(VariableFloat(name_mubar_mass_)),

b_pt_(VariableFloat(name_b_pt_)),
b_eta_(VariableFloat(name_b_eta_)),
b_phi_(VariableFloat(name_b_phi_)),
b_mass_(VariableFloat(name_b_mass_)),

bbar_pt_(VariableFloat(name_bbar_pt_)),
bbar_eta_(VariableFloat(name_bbar_eta_)),
bbar_phi_(VariableFloat(name_bbar_phi_)),
bbar_mass_(VariableFloat(name_bbar_mass_)),

nu_pt_(VariableFloat(name_nu_pt_)),
nu_eta_(VariableFloat(name_nu_eta_)),
nu_phi_(VariableFloat(name_nu_phi_)),
nu_mass_(VariableFloat(name_nu_mass_)),

nubar_pt_(VariableFloat(name_nubar_pt_)),
nubar_eta_(VariableFloat(name_nubar_eta_)),
nubar_phi_(VariableFloat(name_nubar_phi_)),
nubar_mass_(VariableFloat(name_nubar_mass_)),

met_pt_(VariableFloat(name_met_pt_)),
met_phi_(VariableFloat(name_met_phi_)),
met_mass_(VariableFloat(name_met_mass_)),

ttbar_pt_(VariableFloat(name_ttbar_pt_)),
ttbar_phi_(VariableFloat(name_ttbar_phi_)),
ttbar_rapidity_(VariableFloat(name_ttbar_rapidity_)),
//ttbar_arapidity_(VariableFloat(name_ttbar_arapidity_)),
ttbar_eta_(VariableFloat(name_ttbar_eta_)),
ttbar_delta_phi_(VariableFloat(name_ttbar_delta_phi_)),
ttbar_delta_eta_(VariableFloat(name_ttbar_delta_eta_)),
ttbar_delta_rapidity_(VariableFloat(name_ttbar_delta_rapidity_)),
ttbar_mass_(VariableFloat(name_ttbar_mass_)),

llbar_pt_(VariableFloat(name_llbar_pt_)),
llbar_phi_(VariableFloat(name_llbar_phi_)),
llbar_rapidity_(VariableFloat(name_llbar_rapidity_)),
//llbar_arapidity_(VariableFloat(name_llbar_arapidity_)),
llbar_delta_phi_(VariableFloat(name_llbar_delta_phi_)),
llbar_delta_eta_(VariableFloat(name_llbar_delta_eta_)),
llbar_delta_rapidity_(VariableFloat(name_llbar_delta_rapidity_)),
llbar_mass_(VariableFloat(name_llbar_mass_)),

//bbbar_pt_(VariableFloat(name_bbbar_pt_)),
//bbbar_phi_(VariableFloat(name_bbbar_phi_)),
//bbbar_rapidity_(VariableFloat(name_bbbar_rapidity_)),
//bbbar_arapidity_(VariableFloat(name_bbbar_arapidity_)),
//bbbar_delta_phi_(VariableFloat(name_bbbar_delta_phi_)),
//bbbar_delta_eta_(VariableFloat(name_bbbar_delta_eta_)),
//bbbar_delta_rapidity_(VariableFloat(name_bbbar_delta_rapidity_)),
//bbbar_mass_(VariableFloat(name_llbar_mass_)),

//nunubar_pt_(VariableFloat(name_nunubar_pt_)),
//nunubar_phi_(VariableFloat(name_nunubar_phi_)),
//nunubar_rapidity_(VariableFloat(name_nunubar_rapidity_)),
//nunubar_arapidity_(VariableFloat(name_nunubar_arapidity_)),
//nunubar_delta_phi_(VariableFloat(name_nunubar_delta_phi_)),
//nunubar_delta_eta_(VariableFloat(name_nunubar_delta_eta_)),
//nunubar_delta_rapidity_(VariableFloat(name_nunubar_delta_rapidity_)),
//nunubar_mass_(VariableFloat(name_nunubar_mass_)),

MT2_(VariableFloat(name_MT2_)),

all_mass_(VariableFloat(name_all_mass_)),
r_mass_(VariableFloat(name_r_mass_)),
jet_multiplicity_(VariableFloat(name_jet_multiplicity_)),
bjet_multiplicity_(VariableFloat(name_bjet_multiplicity_)),
x1_(VariableFloat(name_x1_)),
x2_(VariableFloat(name_x2_)),

//ExtraJet variables
n_extraJets_iso08_(VariableFloat(name_n_extraJets_iso08_)),

//ttbar frame variables
top_scatteringangle_ttbarframe_(VariableFloat(name_top_scatteringangle_ttbarframe_)),
//top_pt_ttbarframe_(VariableFloat(name_top_pt_ttbarframe_)),
//top_phi_ttbarframe_(VariableFloat(name_top_phi_ttbarframe_)),
//top_rapidity_ttbarframe_(VariableFloat(name_top_rapidity_ttbarframe_)),
//top_arapidity_ttbarframe_(VariableFloat(name_top_arapidity_ttbarframe_)),
//tbar_pt_ttbarframe_(VariableFloat(name_tbar_pt_ttbarframe_)),
//tbar_phi_ttbarframe_(VariableFloat(name_tbar_phi_ttbarframe_)),
//tbar_rapidity_ttbarframe_(VariableFloat(name_tbar_rapidity_ttbarframe_)),
//tbar_arapidity_ttbarframe_(VariableFloat(name_tbar_arapidity_ttbarframe_)),

// top and antitop frame variables
//l_pt_antitopframe_(VariableFloat(name_l_pt_antitopframe_)),
//l_eta_antitopframe_(VariableFloat(name_l_eta_antitopframe_)),
//l_phi_antitopframe_(VariableFloat(name_l_phi_antitopframe_)),
//lbar_pt_topframe_(VariableFloat(name_lbar_pt_topframe_)),
//lbar_eta_topframe_(VariableFloat(name_lbar_eta_topframe_)),
//lbar_phi_topframe_(VariableFloat(name_lbar_phi_topframe_)),

//Spin correlation variables
b1k_(VariableFloat(name_b1k_)),
b2k_(VariableFloat(name_b2k_)),
b1j_(VariableFloat(name_b1j_)),
b2j_(VariableFloat(name_b2j_)),
b1r_(VariableFloat(name_b1r_)),
b2r_(VariableFloat(name_b2r_)),
b1q_(VariableFloat(name_b1q_)),
b2q_(VariableFloat(name_b2q_)),
b1n_(VariableFloat(name_b1n_)),
b2n_(VariableFloat(name_b2n_)),

//b_Pkk_(VariableFloat(name_b_Pkk_)),
//b_Mkk_(VariableFloat(name_b_Mkk_)),
//b_Pjj_(VariableFloat(name_b_Pjj_)),
//b_Mjj_(VariableFloat(name_b_Mjj_)),
//b_Prr_(VariableFloat(name_b_Prr_)),
//b_Mrr_(VariableFloat(name_b_Mrr_)),
//b_Pqq_(VariableFloat(name_b_Pqq_)),
//b_Mqq_(VariableFloat(name_b_Mqq_)),
//b_Pnn_(VariableFloat(name_b_Pnn_)),
//b_Mnn_(VariableFloat(name_b_Mnn_)),

c_kk_(VariableFloat(name_c_kk_)),
c_rr_(VariableFloat(name_c_rr_)),
c_nn_(VariableFloat(name_c_nn_)),
c_kj_(VariableFloat(name_c_kj_)),
c_rq_(VariableFloat(name_c_rq_)),

c_rk_(VariableFloat(name_c_rk_)),
c_kr_(VariableFloat(name_c_kr_)),
c_nr_(VariableFloat(name_c_nr_)),
c_rn_(VariableFloat(name_c_rn_)),
c_nk_(VariableFloat(name_c_nk_)),
c_kn_(VariableFloat(name_c_kn_)),

c_rj_(VariableFloat(name_c_rj_)),
c_jr_(VariableFloat(name_c_jr_)),
//c_qk_(VariableFloat(name_c_qk_)),
//c_kq_(VariableFloat(name_c_kq_)),

//c_qj_(VariableFloat(name_c_qj_)),
//c_jq_(VariableFloat(name_c_jq_)),
//c_nq_(VariableFloat(name_c_nq_)),
//c_qn_(VariableFloat(name_c_qn_)),
//c_nj_(VariableFloat(name_c_nj_)),
//c_jn_(VariableFloat(name_c_jn_)),

//c_Prk_(VariableFloat(name_c_Prk_)),
//c_Mrk_(VariableFloat(name_c_Mrk_)),
//c_Pnr_(VariableFloat(name_c_Pnr_)),
//c_Mnr_(VariableFloat(name_c_Mnr_)),
//c_Pnk_(VariableFloat(name_c_Pnk_)),
//c_Mnk_(VariableFloat(name_c_Mnk_)),

//c_Pqj_(VariableFloat(name_c_Pqj_)),
//c_Mqj_(VariableFloat(name_c_Mqj_)),
//c_Pnq_(VariableFloat(name_c_Pnq_)),
//c_Mnq_(VariableFloat(name_c_Mnq_)),
//c_Pnj_(VariableFloat(name_c_Pnj_)),
//c_Mnj_(VariableFloat(name_c_Mnj_)),

//c_han_(VariableFloat(name_c_han_)),
//c_sca_(VariableFloat(name_c_sca_)),
//c_tra_(VariableFloat(name_c_tra_)),
//c_lij_(VariableFloat(name_c_lij_)),
//c_liq_(VariableFloat(name_c_liq_)),

//c_rkP_(VariableFloat(name_c_rkP_)),
//c_rkM_(VariableFloat(name_c_rkM_)),
//c_nrP_(VariableFloat(name_c_nrP_)),
//c_nrM_(VariableFloat(name_c_nrM_)),
//c_nkP_(VariableFloat(name_c_nkP_)),
//c_nkM_(VariableFloat(name_c_nkM_)),

//c_qjP_(VariableFloat(name_c_qjP_)),
//c_qjM_(VariableFloat(name_c_qjM_)),
//c_nqP_(VariableFloat(name_c_nqP_)),
//c_nqM_(VariableFloat(name_c_nqM_)),
//c_njP_(VariableFloat(name_c_njP_)),
//c_njM_(VariableFloat(name_c_njM_)),

ll_cHel_(VariableFloat(name_ll_cHel_)),
ll_cLab_(VariableFloat(name_ll_cLab_)),
ll_kNorm_(VariableFloat(name_ll_kNorm_)),
ll_rNorm_(VariableFloat(name_ll_rNorm_)),

gen_top_pt_(VariableFloat(name_gen_top_pt_)),
gen_top_phi_(VariableFloat(name_gen_top_phi_)),
gen_top_rapidity_(VariableFloat(name_gen_top_rapidity_)),
//gen_top_arapidity_(VariableFloat(name_gen_top_arapidity_)),
gen_top_eta_(VariableFloat(name_gen_top_eta_)),
gen_top_mass_(VariableFloat(name_gen_top_mass_)),

gen_tbar_pt_(VariableFloat(name_gen_tbar_pt_)),
gen_tbar_phi_(VariableFloat(name_gen_tbar_phi_)),
gen_tbar_rapidity_(VariableFloat(name_gen_tbar_rapidity_)),
//gen_tbar_arapidity_(VariableFloat(name_gen_tbar_arapidity_)),
gen_tbar_eta_(VariableFloat(name_gen_tbar_eta_)),
gen_tbar_mass_(VariableFloat(name_gen_tbar_mass_)),

gen_l_pt_(VariableFloat(name_gen_l_pt_)),
gen_l_eta_(VariableFloat(name_gen_l_eta_)),
gen_l_phi_(VariableFloat(name_gen_l_phi_)),
gen_l_mass_(VariableFloat(name_gen_l_mass_)),
gen_l_pdgid_(VariableFloat(name_gen_l_pdgid_)),

gen_lbar_pt_(VariableFloat(name_gen_lbar_pt_)),
gen_lbar_eta_(VariableFloat(name_gen_lbar_eta_)),
gen_lbar_phi_(VariableFloat(name_gen_lbar_phi_)),
gen_lbar_mass_(VariableFloat(name_gen_lbar_mass_)),
gen_lbar_pdgid_(VariableFloat(name_gen_lbar_pdgid_)),

//gen_e_pt_(VariableFloat(name_gen_e_pt_)),
//gen_e_eta_(VariableFloat(name_gen_e_eta_)),
//gen_e_phi_(VariableFloat(name_gen_e_phi_)),
//gen_e_mass_(VariableFloat(name_gen_e_mass_)),

//gen_ebar_pt_(VariableFloat(name_gen_ebar_pt_)),
//gen_ebar_eta_(VariableFloat(name_gen_ebar_eta_)),
//gen_ebar_phi_(VariableFloat(name_gen_ebar_phi_)),
//gen_ebar_mass_(VariableFloat(name_gen_ebar_mass_)),

//gen_mu_pt_(VariableFloat(name_gen_mu_pt_)),
//gen_mu_eta_(VariableFloat(name_gen_mu_eta_)),
//gen_mu_phi_(VariableFloat(name_gen_mu_phi_)),
//gen_mu_mass_(VariableFloat(name_gen_mu_mass_)),

//gen_mubar_pt_(VariableFloat(name_gen_mubar_pt_)),
//gen_mubar_eta_(VariableFloat(name_gen_mubar_eta_)),
//gen_mubar_phi_(VariableFloat(name_gen_mubar_phi_)),
//gen_mubar_mass_(VariableFloat(name_gen_mubar_mass_)),

gen_b_pt_(VariableFloat(name_gen_b_pt_)),
gen_b_eta_(VariableFloat(name_gen_b_eta_)),
gen_b_phi_(VariableFloat(name_gen_b_phi_)),
gen_b_mass_(VariableFloat(name_gen_b_mass_)),

gen_bbar_pt_(VariableFloat(name_gen_bbar_pt_)),
gen_bbar_eta_(VariableFloat(name_gen_bbar_eta_)),
gen_bbar_phi_(VariableFloat(name_gen_bbar_phi_)),
gen_bbar_mass_(VariableFloat(name_gen_bbar_mass_)),

gen_nu_pt_(VariableFloat(name_gen_nu_pt_)),
gen_nu_eta_(VariableFloat(name_gen_nu_eta_)),
gen_nu_phi_(VariableFloat(name_gen_nu_phi_)),
gen_nu_mass_(VariableFloat(name_gen_nu_mass_)),

gen_nubar_pt_(VariableFloat(name_gen_nubar_pt_)),
gen_nubar_eta_(VariableFloat(name_gen_nubar_eta_)),
gen_nubar_phi_(VariableFloat(name_gen_nubar_phi_)),
gen_nubar_mass_(VariableFloat(name_gen_nubar_mass_)),

gen_ttbar_pt_(VariableFloat(name_gen_ttbar_pt_)),
gen_ttbar_phi_(VariableFloat(name_gen_ttbar_phi_)),
gen_ttbar_rapidity_(VariableFloat(name_gen_ttbar_rapidity_)),
//gen_ttbar_arapidity_(VariableFloat(name_gen_ttbar_arapidity_)),
gen_ttbar_eta_(VariableFloat(name_gen_ttbar_eta_)),
gen_ttbar_delta_phi_(VariableFloat(name_gen_ttbar_delta_phi_)),
gen_ttbar_delta_eta_(VariableFloat(name_gen_ttbar_delta_eta_)),
gen_ttbar_delta_rapidity_(VariableFloat(name_gen_ttbar_delta_rapidity_)),
gen_ttbar_mass_(VariableFloat(name_gen_ttbar_mass_)),

gen_llbar_pt_(VariableFloat(name_gen_llbar_pt_)),
gen_llbar_phi_(VariableFloat(name_gen_llbar_phi_)),
gen_llbar_rapidity_(VariableFloat(name_gen_llbar_rapidity_)),
//gen_llbar_arapidity_(VariableFloat(name_gen_llbar_arapidity_)),
gen_llbar_delta_phi_(VariableFloat(name_gen_llbar_delta_phi_)),
gen_llbar_delta_eta_(VariableFloat(name_gen_llbar_delta_eta_)),
gen_llbar_delta_rapidity_(VariableFloat(name_gen_llbar_delta_rapidity_)),
gen_llbar_mass_(VariableFloat(name_gen_llbar_mass_)),

//gen_bbbar_pt_(VariableFloat(name_gen_bbbar_pt_)),
//gen_bbbar_phi_(VariableFloat(name_gen_bbbar_phi_)),
//gen_bbbar_rapidity_(VariableFloat(name_gen_bbbar_rapidity_)),
//gen_bbbar_arapidity_(VariableFloat(name_gen_bbbar_arapidity_)),
//gen_bbbar_delta_phi_(VariableFloat(name_gen_bbbar_delta_phi_)),
//gen_bbbar_delta_eta_(VariableFloat(name_gen_bbbar_delta_eta_)),
//gen_bbbar_delta_rapidity_(VariableFloat(name_gen_bbbar_delta_rapidity_)),
//gen_bbbar_mass_(VariableFloat(name_gen_bbbar_mass_)),

//gen_nunubar_pt_(VariableFloat(name_gen_nunubar_pt_)),
//gen_nunubar_phi_(VariableFloat(name_gen_nunubar_phi_)),
//gen_nunubar_rapidity_(VariableFloat(name_gen_nunubar_rapidity_)),
//gen_nunubar_arapidity_(VariableFloat(name_gen_nunubar_arapidity_)),
//gen_nunubar_delta_phi_(VariableFloat(name_gen_nunubar_delta_phi_)),
//gen_nunubar_delta_eta_(VariableFloat(name_gen_nunubar_delta_eta_)),
//gen_nunubar_delta_rapidity_(VariableFloat(name_gen_nunubar_delta_rapidity_)),
//gen_nunubar_mass_(VariableFloat(name_gen_nunubar_mass_)),

gen_MT2_(VariableFloat(name_gen_MT2_)),

gen_all_mass_(VariableFloat(name_gen_all_mass_)),
gen_r_mass_(VariableFloat(name_gen_r_mass_)),
gen_jet_multiplicity_(VariableFloat(name_gen_jet_multiplicity_)),
gen_x1_(VariableFloat(name_gen_x1_)),
gen_x2_(VariableFloat(name_gen_x2_)),
gen_production_mode_(VariableFloat(name_gen_production_mode_)),

//ExtraJet variables
gen_n_extraJets_iso08_(VariableFloat(name_gen_n_extraJets_iso08_)),

//ttbar frame variables
gen_top_scatteringangle_ttbarframe_(VariableFloat(name_gen_top_scatteringangle_ttbarframe_)),
//gen_top_pt_ttbarframe_(VariableFloat(name_gen_top_pt_ttbarframe_)),
//gen_top_phi_ttbarframe_(VariableFloat(name_gen_top_phi_ttbarframe_)),
//gen_top_rapidity_ttbarframe_(VariableFloat(name_gen_top_rapidity_ttbarframe_)),
//gen_top_arapidity_ttbarframe_(VariableFloat(name_gen_top_arapidity_ttbarframe_)),
//gen_tbar_pt_ttbarframe_(VariableFloat(name_gen_tbar_pt_ttbarframe_)),
//gen_tbar_phi_ttbarframe_(VariableFloat(name_gen_tbar_phi_ttbarframe_)),
//gen_tbar_rapidity_ttbarframe_(VariableFloat(name_gen_tbar_rapidity_ttbarframe_)),
//gen_tbar_arapidity_ttbarframe_(VariableFloat(name_gen_tbar_arapidity_ttbarframe_)),

// top and antitop frame variables
//gen_l_pt_antitopframe_(VariableFloat(name_gen_l_pt_antitopframe_)),
//gen_l_eta_antitopframe_(VariableFloat(name_gen_l_eta_antitopframe_)),
//gen_l_phi_antitopframe_(VariableFloat(name_gen_l_phi_antitopframe_)),
//gen_lbar_pt_topframe_(VariableFloat(name_gen_lbar_pt_topframe_)),
//gen_lbar_eta_topframe_(VariableFloat(name_gen_lbar_eta_topframe_)),
//gen_lbar_phi_topframe_(VariableFloat(name_gen_lbar_phi_topframe_)),

//Spin correlation variables
gen_b1k_(VariableFloat(name_gen_b1k_)),
gen_b2k_(VariableFloat(name_gen_b2k_)),
gen_b1j_(VariableFloat(name_gen_b1j_)),
gen_b2j_(VariableFloat(name_gen_b2j_)),
gen_b1r_(VariableFloat(name_gen_b1r_)),
gen_b2r_(VariableFloat(name_gen_b2r_)),
gen_b1q_(VariableFloat(name_gen_b1q_)),
gen_b2q_(VariableFloat(name_gen_b2q_)),
gen_b1n_(VariableFloat(name_gen_b1n_)),
gen_b2n_(VariableFloat(name_gen_b2n_)),

//gen_b_Pkk_(VariableFloat(name_gen_b_Pkk_)),
//gen_b_Mkk_(VariableFloat(name_gen_b_Mkk_)),
//gen_b_Pjj_(VariableFloat(name_gen_b_Pjj_)),
//gen_b_Mjj_(VariableFloat(name_gen_b_Mjj_)),
//gen_b_Prr_(VariableFloat(name_gen_b_Prr_)),
//gen_b_Mrr_(VariableFloat(name_gen_b_Mrr_)),
//gen_b_Pqq_(VariableFloat(name_gen_b_Pqq_)),
//gen_b_Mqq_(VariableFloat(name_gen_b_Mqq_)),
//gen_b_Pnn_(VariableFloat(name_gen_b_Pnn_)),
//gen_b_Mnn_(VariableFloat(name_gen_b_Mnn_)),

gen_c_kk_(VariableFloat(name_gen_c_kk_)),
gen_c_rr_(VariableFloat(name_gen_c_rr_)),
gen_c_nn_(VariableFloat(name_gen_c_nn_)),
gen_c_kj_(VariableFloat(name_gen_c_kj_)),
gen_c_rq_(VariableFloat(name_gen_c_rq_)),

gen_c_rk_(VariableFloat(name_gen_c_rk_)),
gen_c_kr_(VariableFloat(name_gen_c_kr_)),
gen_c_nr_(VariableFloat(name_gen_c_nr_)),
gen_c_rn_(VariableFloat(name_gen_c_rn_)),
gen_c_nk_(VariableFloat(name_gen_c_nk_)),
gen_c_kn_(VariableFloat(name_gen_c_kn_)),

gen_c_rj_(VariableFloat(name_gen_c_rj_)),
gen_c_jr_(VariableFloat(name_gen_c_jr_)),
//gen_c_qk_(VariableFloat(name_gen_c_qk_)),
//gen_c_kq_(VariableFloat(name_gen_c_kq_)),

//gen_c_qj_(VariableFloat(name_gen_c_qj_)),
//gen_c_jq_(VariableFloat(name_gen_c_jq_)),
//gen_c_nq_(VariableFloat(name_gen_c_nq_)),
//gen_c_qn_(VariableFloat(name_gen_c_qn_)),
//gen_c_nj_(VariableFloat(name_gen_c_nj_)),
//gen_c_jn_(VariableFloat(name_gen_c_jn_)),

//gen_c_Prk_(VariableFloat(name_gen_c_Prk_)),
//gen_c_Mrk_(VariableFloat(name_gen_c_Mrk_)),
//gen_c_Pnr_(VariableFloat(name_gen_c_Pnr_)),
//gen_c_Mnr_(VariableFloat(name_gen_c_Mnr_)),
//gen_c_Pnk_(VariableFloat(name_gen_c_Pnk_)),
//gen_c_Mnk_(VariableFloat(name_gen_c_Mnk_)),

//gen_c_Pqj_(VariableFloat(name_gen_c_Pqj_)),
//gen_c_Mqj_(VariableFloat(name_gen_c_Mqj_)),
//gen_c_Pnq_(VariableFloat(name_gen_c_Pnq_)),
//gen_c_Mnq_(VariableFloat(name_gen_c_Mnq_)),
//gen_c_Pnj_(VariableFloat(name_gen_c_Pnj_)),
//gen_c_Mnj_(VariableFloat(name_gen_c_Mnj_)),

//gen_c_han_(VariableFloat(name_gen_c_han_)),
//gen_c_sca_(VariableFloat(name_gen_c_sca_)),
//gen_c_tra_(VariableFloat(name_gen_c_tra_)),
//gen_c_lij_(VariableFloat(name_gen_c_lij_)),
//gen_c_liq_(VariableFloat(name_gen_c_liq_)),

//gen_c_rkP_(VariableFloat(name_gen_c_rkP_)),
//gen_c_rkM_(VariableFloat(name_gen_c_rkM_)),
//gen_c_nrP_(VariableFloat(name_gen_c_nrP_)),
//gen_c_nrM_(VariableFloat(name_gen_c_nrM_)),
//gen_c_nkP_(VariableFloat(name_gen_c_nkP_)),
//gen_c_nkM_(VariableFloat(name_gen_c_nkM_)),

//gen_c_qjP_(VariableFloat(name_gen_c_qjP_)),
//gen_c_qjM_(VariableFloat(name_gen_c_qjM_)),
//gen_c_nqP_(VariableFloat(name_gen_c_nqP_)),
//gen_c_nqM_(VariableFloat(name_gen_c_nqM_)),
//gen_c_njP_(VariableFloat(name_gen_c_njP_)),
//gen_c_njM_(VariableFloat(name_gen_c_njM_)),

gen_ll_cHel_(VariableFloat(name_gen_ll_cHel_)),
gen_ll_cLab_(VariableFloat(name_gen_ll_cLab_)),
gen_ll_kNorm_(VariableFloat(name_gen_ll_kNorm_)),
gen_ll_rNorm_(VariableFloat(name_gen_ll_rNorm_)),

entry_(VariableInt(name_entry_)),
isTopGen_(VariableInt(name_isTopGen_)),
TopDecayMode_(VariableInt(name_TopDecayMode_)),
isKinReco_(VariableInt(name_isKinReco_)),
eventWeight_(VariableFloat(name_eventWeight_)),
trueLevelWeight_(VariableFloat(name_trueLevelWeight_))
{
    isKinReco_.value_ = 0;
    isTopGen_.value_ = 0;
    TopDecayMode_.value_ = -1;
    entry_.value_ = -999;

    // Use utilities without namespaces
    using ROOT::Math::VectorUtil::DeltaPhi;
    using ROOT::Math::VectorUtil::DeltaR;
    using ROOT::Math::VectorUtil::Angle;
    using namespace common;
    
   //proton Energy [GeV]  
   double protonE = 6500;
    
  if(kinematicReconstructionSolutions.numberOfSolutions()){
   
   isKinReco_.value_ = 1;
   eventWeight_.value_ = recoLevelWeights.weight_;

   TLorentzVector hyptop(common::LVtoTLV(kinematicReconstructionSolutions.solution().top()));
   TLorentzVector hypantitop(common::LVtoTLV(kinematicReconstructionSolutions.solution().antiTop()));
   TLorentzVector hypttbar(hyptop+hypantitop);
   TLorentzVector hyplep(common::LVtoTLV(kinematicReconstructionSolutions.solution().lepton()));
   TLorentzVector hypantilep(common::LVtoTLV(kinematicReconstructionSolutions.solution().antiLepton()));
   TLorentzVector hypllbar(hyplep+hypantilep);
   TLorentzVector hypBjet(common::LVtoTLV(kinematicReconstructionSolutions.solution().bjet()));
   TLorentzVector hypAntiBjet(common::LVtoTLV(kinematicReconstructionSolutions.solution().antiBjet()));
   //   TLorentzVector hypbbbar(hypBjet+hypAntiBjet);
   TLorentzVector hypNeutrino(common::LVtoTLV(kinematicReconstructionSolutions.solution().neutrino()));
   TLorentzVector hypAntiNeutrino(common::LVtoTLV(kinematicReconstructionSolutions.solution().antiNeutrino()));
   TLorentzVector hypnunubar(hypNeutrino+hypAntiNeutrino);

   TVector3 TTBarFrameBoost(-1. * hypttbar.BoostVector());
   TVector3 TopFrameBoost(-1. * hyptop.BoostVector());
   TVector3 AntiTopFrameBoost(-1. * hypantitop.BoostVector());

   TLorentzVector hyptop_ttbarframe = hyptop;
   hyptop_ttbarframe.Boost(TTBarFrameBoost);
   TLorentzVector hypantitop_ttbarframe = hypantitop;
   hypantitop_ttbarframe.Boost(TTBarFrameBoost);
   TLorentzVector hyplep_antitopframe = hyplep;
   hyplep_antitopframe.Boost(AntiTopFrameBoost);
   TLorentzVector hypantilep_topframe = hypantilep;
   hypantilep_topframe.Boost(TopFrameBoost);

   //Evaluate discr based on the vars                                                                                                                                                      
   VariablesPhiTT *vars = VariablesPhiTT::fillVariables(eventMetadata, recoObjects, commonGenObjects, topGenObjects, kinematicReconstructionSolutions, recoObjectIndices, genObjectIndices, genLevelWeights, recoLevelWeights, weight);
   b1k_.value_ = vars->sol_b1k.value_;
   b2k_.value_ = vars->sol_b2k.value_;
   b1j_.value_ = vars->sol_b1j.value_;
   b2j_.value_ = vars->sol_b2j.value_;
   b1r_.value_ = vars->sol_b1r.value_;
   b2r_.value_ = vars->sol_b2r.value_;
   b1q_.value_ = vars->sol_b1q.value_;
   b2q_.value_ = vars->sol_b2q.value_;
   b1n_.value_ = vars->sol_b1n.value_;
   b2n_.value_ = vars->sol_b2n.value_;

   //b_Pkk_.value_ = vars->sol_bP_kk.value_;
   //b_Mkk_.value_ = vars->sol_bM_kk.value_;
   //b_Pjj_.value_ = vars->sol_bP_jj.value_;
   //b_Mjj_.value_ = vars->sol_bM_jj.value_;
   //b_Prr_.value_ = vars->sol_bP_rr.value_;
   //b_Mrr_.value_ = vars->sol_bM_rr.value_;
   //b_Pqq_.value_ = vars->sol_bP_qq.value_;
   //b_Mqq_.value_ = vars->sol_bM_qq.value_;
   //b_Pnn_.value_ = vars->sol_bP_nn.value_;
   //b_Mnn_.value_ = vars->sol_bM_nn.value_;

   c_kk_.value_ = vars->sol_ckk.value_;
   c_rr_.value_ = vars->sol_crr.value_;
   c_nn_.value_ = vars->sol_cnn.value_;
   c_kj_.value_ = vars->sol_ckj.value_;
   c_rq_.value_ = vars->sol_crq.value_;

   c_rk_.value_ = vars->sol_crk.value_;
   c_kr_.value_ = vars->sol_ckr.value_;
   c_nr_.value_ = vars->sol_cnr.value_;
   c_rn_.value_ = vars->sol_crn.value_;
   c_nk_.value_ = vars->sol_cnk.value_;
   c_kn_.value_ = vars->sol_ckn.value_;

   c_rj_.value_ = vars->sol_crj.value_;
   c_jr_.value_ = vars->sol_cjr.value_;
   //c_qk_.value_ = vars->sol_cqk.value_;
   //c_kq_.value_ = vars->sol_ckq.value_;

   //c_qj_.value_ = vars->sol_cqj.value_;
   //c_jq_.value_ = vars->sol_cjq.value_;
   //c_nq_.value_ = vars->sol_cnq.value_;
   //c_qn_.value_ = vars->sol_cqn.value_;
   //c_nj_.value_ = vars->sol_cnj.value_;
   //c_jn_.value_ = vars->sol_cjn.value_;

   //c_Prk_.value_ = vars->sol_cP_rk.value_;
   //c_Mrk_.value_ = vars->sol_cM_rk.value_;
   //c_Pnr_.value_ = vars->sol_cP_nr.value_;
   //c_Mnr_.value_ = vars->sol_cM_nr.value_;
   //c_Pnk_.value_ = vars->sol_cP_nk.value_;
   //c_Mnk_.value_ = vars->sol_cM_nk.value_;

   //c_Pqj_.value_ = vars->sol_cqj.value_ + vars->sol_cjq.value_;
   //c_Mqj_.value_ = vars->sol_cqj.value_ - vars->sol_cjq.value_;
   //c_Pnq_.value_ = vars->sol_cnq.value_ + vars->sol_cqn.value_;
   //c_Mnq_.value_ = vars->sol_cnq.value_ - vars->sol_cqn.value_;
   //c_Pnj_.value_ = vars->sol_cnj.value_ + vars->sol_cjn.value_;
   //c_Mnj_.value_ = vars->sol_cnj.value_ + vars->sol_cjn.value_;

   //c_han_.value_ = +vars->sol_ckk.value_ - vars->sol_crr.value_ - vars->sol_cnn.value_;
   //c_sca_.value_ = -vars->sol_ckk.value_ + vars->sol_crr.value_ - vars->sol_cnn.value_;
   //c_tra_.value_ = -vars->sol_ckk.value_ - vars->sol_crr.value_ + vars->sol_cnn.value_;
   //c_lij_.value_ = -vars->sol_ckj.value_ - vars->sol_crr.value_ - vars->sol_cnn.value_;
   //c_liq_.value_ = -vars->sol_ckk.value_ - vars->sol_crq.value_ - vars->sol_cnn.value_;

   //c_rkP_.value_ = -vars->sol_crk.value_ - vars->sol_crk.value_ - vars->sol_cnn.value_;
   //c_rkM_.value_ = -vars->sol_crk.value_ + vars->sol_crk.value_ - vars->sol_cnn.value_;
   //c_nrP_.value_ = -vars->sol_cnr.value_ - vars->sol_crn.value_ - vars->sol_ckk.value_;
   //c_nrM_.value_ = -vars->sol_cnr.value_ + vars->sol_crn.value_ - vars->sol_ckk.value_;
   //c_nkP_.value_ = -vars->sol_cnk.value_ - vars->sol_cnk.value_ - vars->sol_crr.value_;
   //c_nkM_.value_ = -vars->sol_cnk.value_ + vars->sol_cnk.value_ - vars->sol_crr.value_;

   //c_qjP_.value_ = -vars->sol_cqj.value_ - vars->sol_cqj.value_ - vars->sol_cnn.value_;
   //c_qjM_.value_ = -vars->sol_cqj.value_ + vars->sol_cqj.value_ - vars->sol_cnn.value_;
   //c_nqP_.value_ = -vars->sol_cnq.value_ - vars->sol_cqn.value_ - vars->sol_ckj.value_;
   //c_nqM_.value_ = -vars->sol_cnq.value_ + vars->sol_cqn.value_ - vars->sol_ckj.value_;
   //c_njP_.value_ = -vars->sol_cnj.value_ - vars->sol_cnj.value_ - vars->sol_crq.value_;
   //c_njM_.value_ = -vars->sol_cnj.value_ + vars->sol_cnj.value_ - vars->sol_crq.value_;

   ll_cHel_.value_ = vars->sol_cHel.value_;
   ll_cLab_.value_ = vars->sol_cLab.value_;
   ll_kNorm_.value_ = vars->sol_kNorm.value_;
   ll_rNorm_.value_ = vars->sol_rNorm.value_;

   delete vars;

   // Get MET
   LV& met = *recoObjects.met_;
   met_pt_.value_ = met.Pt();
   met_phi_.value_ = met.Phi();

   double hypMT2 = VariablesSpinCorr::getMT2Variable(hyplep, hypantilep, common::LVtoTLV(met));

   int nRecoJets=recoObjectIndices.jetIndices_.size();
   jet_multiplicity_.value_ = nRecoJets;
   
   int numberOfBjets = recoObjectIndices.bjetIndices_.size();
   bjet_multiplicity_.value_= recoObjectIndices.bjetIndices_.size();

   n_extraJets_iso08_.value_ = recoObjectIndices.extraJetIndicesIso08_.size();
   
   top_pt_.value_ = hyptop.Pt();
   top_phi_.value_ = hyptop.Phi();
   top_rapidity_.value_ = hyptop.Rapidity();
   //   top_arapidity_.value_ = fabs(hyptop.Rapidity());
   top_eta_.value_ = hyptop.Eta();
   top_mass_.value_ = hyptop.M();

   tbar_pt_.value_ = hypantitop.Pt();
   tbar_phi_.value_ = hypantitop.Phi();
   tbar_rapidity_.value_ = hypantitop.Rapidity();
   //   tbar_arapidity_.value_ = fabs(hypantitop.Rapidity());
   tbar_eta_.value_ = hypantitop.Eta();
   tbar_mass_.value_ = hypantitop.M();

   ttbar_pt_.value_ = hypttbar.Pt();
   ttbar_phi_.value_ = hypttbar.Phi();
   ttbar_rapidity_.value_ = hypttbar.Rapidity();
   //   ttbar_arapidity_.value_ = fabs(hypttbar.Rapidity());
   ttbar_eta_.value_ = hypttbar.Eta();
   ttbar_mass_.value_ = hypttbar.M();

   TVector3 ProtonBeamDirection(0., 0., 1.);
   top_scatteringangle_ttbarframe_.value_ = hyptop_ttbarframe.Vect().Unit().Dot(ProtonBeamDirection);

   //   top_pt_ttbarframe_.value_ = hyptop_ttbarframe.Pt();
   //   top_phi_ttbarframe_.value_ = hyptop_ttbarframe.Phi();
   //   top_rapidity_ttbarframe_.value_ = hyptop_ttbarframe.Rapidity();

   //   tbar_pt_ttbarframe_.value_ = hypantitop_ttbarframe.Pt();
   //   tbar_phi_ttbarframe_.value_ = hypantitop_ttbarframe.Phi();
   //   tbar_rapidity_ttbarframe_.value_ = hypantitop_ttbarframe.Rapidity();

   //   l_pt_antitopframe_.value_ = hyplep_antitopframe.Pt();
   //   l_eta_antitopframe_.value_ = hyplep_antitopframe.Eta();
   //   l_phi_antitopframe_.value_ = hyplep_antitopframe.Phi();
   //   lbar_pt_topframe_.value_ = hypantilep_topframe.Pt();
   //   lbar_eta_topframe_.value_ = hypantilep_topframe.Eta();
   //   lbar_phi_topframe_.value_ = hypantilep_topframe.Phi();


   LV allRecoObjects = kinematicReconstructionSolutions.solution().antiNeutrino() + kinematicReconstructionSolutions.solution().neutrino();
   for(int i=0;i<(int)recoObjects.allLeptons_->size();i++){allRecoObjects = allRecoObjects + (recoObjects.allLeptons_->at(i));
     //std::cout << "ttbar: " << (recoObjects.allLeptons_->at(i)).Pt() << std::endl;                                                                                                        
   }
   for(int i=0;i<(int)recoObjects.jets_->size();i++)allRecoObjects = allRecoObjects + (recoObjects.jets_->at(i));

   all_mass_.value_ = allRecoObjects.M();

   r_mass_.value_ = ttbar_mass_.value_/all_mass_.value_;

   ttbar_delta_eta_.value_ = fabs(hyptop.Eta()-hypantitop.Eta());
   ttbar_delta_rapidity_.value_ = fabs(hyptop.Rapidity()-hypantitop.Rapidity());
   ttbar_delta_phi_.value_ = fabs(hyptop.DeltaPhi(hypantitop));


   double restPzJetsSum=0; // rest means jets not from top or anti-top
   double restEJetsSum=0; 
   for(int i=0;i<nRecoJets;i++)
   {
       if(i != kinematicReconstructionSolutions.solution().bjetIndex() && i != kinematicReconstructionSolutions.solution().bjetIndex())
       {
           restEJetsSum += (*recoObjects.jets_).at(recoObjectIndices.jetIndices_.at(i)).E();
           restPzJetsSum += (*recoObjects.jets_).at(recoObjectIndices.jetIndices_.at(i)).Pz();
       }
   }
   
   //double x1 = (hyptop.E()+hypantitop.E()+restEJetsSum+hyptop.Pz()+hypantitop.Pz()+restPzJetsSum)/(2*protonE);
   //double x2 = (hyptop.E()+hypantitop.E()+restEJetsSum-hyptop.Pz()-hypantitop.Pz()-restPzJetsSum)/(2*protonE);
   double x1 = (hyptop.E()+hypantitop.E()+hyptop.Pz()+hypantitop.Pz())/(2*protonE);
   double x2 = (hyptop.E()+hypantitop.E()-hyptop.Pz()-hypantitop.Pz())/(2*protonE);
   x1_.value_ = x1;
   x2_.value_ = x2;
  
   l_pt_.value_ = hyplep.Pt();
   l_eta_.value_ = hyplep.Eta();
   l_phi_.value_ = hyplep.Phi();
   l_mass_.value_ = hyplep.M();
   if(recoObjects.lepPdgId_->at(0) < 0) l_pdgid_.value_ = recoObjects.lepPdgId_->at(0);
   else l_pdgid_.value_ = recoObjects.lepPdgId_->at(1);

   lbar_pt_.value_ = hypantilep.Pt();
   lbar_eta_.value_ = hypantilep.Eta();
   lbar_phi_.value_ = hypantilep.Phi();
   lbar_mass_.value_ = hypantilep.M();
   if(recoObjects.lepPdgId_->at(0) > 0) lbar_pdgid_.value_ = recoObjects.lepPdgId_->at(0);
   else lbar_pdgid_.value_ = recoObjects.lepPdgId_->at(1);

   llbar_mass_.value_ = hypllbar.M();
   llbar_pt_.value_=hypllbar.Pt();
   llbar_phi_.value_=hypllbar.Phi();
   llbar_rapidity_.value_=hypllbar.Rapidity();
   //   llbar_arapidity_.value_=fabs(hypllbar.Rapidity());

   llbar_delta_eta_.value_ = fabs(hyplep.Eta()-hypantilep.Eta());
   llbar_delta_rapidity_.value_ = fabs(hyplep.Rapidity()-hypantilep.Rapidity());
   llbar_delta_phi_.value_ = fabs(hyplep.DeltaPhi(hypantilep));

   b_pt_.value_ = hypBjet.Pt();
   bbar_pt_.value_ = hypAntiBjet.Pt();
   b_eta_.value_ = hypBjet.Eta();
   bbar_eta_.value_ = hypAntiBjet.Eta();
   b_phi_.value_ = hypBjet.Phi();
   bbar_phi_.value_ = hypAntiBjet.Phi();
   b_mass_.value_ = hypBjet.M();
   bbar_mass_.value_ = hypAntiBjet.M();

   //   bbbar_mass_.value_ = hypbbbar.M();
   //   bbbar_pt_.value_= hypbbbar.Pt();
   //   bbbar_phi_.value_= hypbbbar.Phi();
   //   bbbar_rapidity_.value_=hypbbbar.Rapidity();
   //   bbbar_arapidity_.value_=fabs(hypbbbar.Rapidity());

   //   bbbar_delta_eta_.value_ = fabs(hypBjet.Eta()-hypAntiBjet.Eta());
   //   bbbar_delta_rapidity_.value_ = fabs(hypBjet.Rapidity()-hypAntiBjet.Rapidity());
   //   bbbar_delta_phi_.value_ = fabs(hypBjet.DeltaPhi(hypAntiBjet));

   nu_pt_.value_ = hypNeutrino.Pt();
   nubar_pt_.value_ = hypAntiNeutrino.Pt();
   nu_eta_.value_ = hypNeutrino.Eta();
   nubar_eta_.value_ = hypAntiNeutrino.Eta();
   nu_phi_.value_ = hypNeutrino.Phi();
   nubar_phi_.value_ = hypAntiNeutrino.Phi();

   //   nunubar_mass_.value_ = hypnunubar.M();
   //   nunubar_pt_.value_= hypnunubar.Pt();
   //   nunubar_phi_.value_= hypnunubar.Phi();
   //   nunubar_rapidity_.value_=hypnunubar.Rapidity();
   //   nunubar_arapidity_.value_=fabs(hypnunubar.Rapidity());

   //   nunubar_delta_eta_.value_ = fabs(hypNeutrino.Eta()-hypAntiNeutrino.Eta());
   //   nunubar_delta_rapidity_.value_ = fabs(hypNeutrino.Rapidity()-hypAntiNeutrino.Rapidity());
   //   nunubar_delta_phi_.value_ = fabs(hypNeutrino.DeltaPhi(hypAntiNeutrino));

   MT2_.value_ = hypMT2; 

   }
   
  TopDecayMode_.value_ = topGenObjects.decayMode_;

   if(topGenObjects.valuesSet_){
       
      entry_.value_ = (Long64_t)eventMetadata.eventNumber_;
      isTopGen_.value_ = 1;
      TopDecayMode_.value_ = topGenObjects.decayMode_;
      trueLevelWeight_.value_ = genLevelWeights.trueLevelWeight_;
      
      gen_jet_multiplicity_.value_ = genObjectIndices.genVisJetIndices_.size();
      gen_n_extraJets_iso08_.value_ = genObjectIndices.genExtraJetIndicesIso08_.size();
    
      TLorentzVector gentop(common::LVtoTLV((*topGenObjects.GenTop_)));
      TLorentzVector genantitop(common::LVtoTLV((*topGenObjects.GenAntiTop_)));
      TLorentzVector genttbar(gentop+genantitop);

      TLorentzVector genlep(common::LVtoTLV((*topGenObjects.GenLepton_)));
      TLorentzVector genantilep(common::LVtoTLV((*topGenObjects.GenAntiLepton_)));
      TLorentzVector genllbar(genlep+genantilep);

      TLorentzVector genb(common::LVtoTLV((*topGenObjects.GenB_)));
      TLorentzVector genantib(common::LVtoTLV((*topGenObjects.GenAntiB_)));
      //      TLorentzVector genbbbar(genb+genantib);

      TLorentzVector gennu(common::LVtoTLV((*topGenObjects.GenNeutrino_)));
      TLorentzVector genantinu(common::LVtoTLV((*topGenObjects.GenAntiNeutrino_)));
      TLorentzVector gennunubar(gennu+genantinu);

      TLorentzVector genmet;
      genmet.SetPxPyPzE(gennunubar.Px(), gennunubar.Py(), 0.0, gennunubar.Pt());
      double genMT2 = VariablesSpinCorr::getMT2Variable(genlep, genantilep, genmet);


      TVector3 GenTTBarFrameBoost(-1. * genttbar.BoostVector());
      TVector3 GenTopFrameBoost(-1. * gentop.BoostVector());
      TVector3 GenAntiTopFrameBoost(-1. * genantitop.BoostVector());

      TLorentzVector gentop_ttbarframe = gentop;
      gentop_ttbarframe.Boost(GenTTBarFrameBoost);
      TLorentzVector genantitop_ttbarframe = genantitop;
      genantitop_ttbarframe.Boost(GenTTBarFrameBoost);
      TLorentzVector genlep_antitopframe = genlep;
      genlep_antitopframe.Boost(GenAntiTopFrameBoost);
      TLorentzVector genantilep_topframe = genantilep;
      genantilep_topframe.Boost(GenTopFrameBoost);


      TVector3 ProtonBeamDirection(0., 0., 1.);
      gen_top_scatteringangle_ttbarframe_.value_ = gentop_ttbarframe.Vect().Unit().Dot(ProtonBeamDirection);
      

      //      gen_top_pt_ttbarframe_.value_ = gentop_ttbarframe.Pt();
      //      gen_top_phi_ttbarframe_.value_ = gentop_ttbarframe.Phi();
      //      gen_top_rapidity_ttbarframe_.value_ = gentop_ttbarframe.Rapidity();

      //      gen_tbar_pt_ttbarframe_.value_ = genantitop_ttbarframe.Pt();
      //      gen_tbar_phi_ttbarframe_.value_ = genantitop_ttbarframe.Phi();
      //      gen_tbar_rapidity_ttbarframe_.value_ = genantitop_ttbarframe.Rapidity();
      
      //      gen_l_pt_antitopframe_.value_ = genlep_antitopframe.Pt();
      //      gen_l_eta_antitopframe_.value_ = genlep_antitopframe.Eta();
      //      gen_l_phi_antitopframe_.value_ = genlep_antitopframe.Phi();
      
      //      gen_lbar_pt_topframe_.value_ = genantilep_topframe.Pt();
      //      gen_lbar_eta_topframe_.value_ = genantilep_topframe.Eta();
      //      gen_lbar_phi_topframe_.value_ = genantilep_topframe.Phi();

      gen_top_pt_.value_ = gentop.Pt();
      gen_top_phi_.value_ = gentop.Phi();
      gen_top_rapidity_.value_ = gentop.Rapidity();
      //      gen_top_arapidity_.value_ = fabs(gentop.Rapidity());
      gen_top_eta_.value_ = gentop.Eta();
      gen_top_mass_.value_ = gentop.M();

      gen_tbar_pt_.value_ = genantitop.Pt();
      gen_tbar_phi_.value_ = genantitop.Phi();
      gen_tbar_rapidity_.value_ = genantitop.Rapidity();
      //      gen_tbar_arapidity_.value_ = fabs(genantitop.Rapidity());
      gen_tbar_eta_.value_ = genantitop.Eta();
      gen_tbar_mass_.value_ = genantitop.M();

      gen_ttbar_pt_.value_ = genttbar.Pt();
      gen_ttbar_phi_.value_ = genttbar.Phi();
      gen_ttbar_rapidity_.value_ = genttbar.Rapidity();
      //      gen_ttbar_arapidity_.value_ = fabs(genttbar.Rapidity());
      gen_ttbar_eta_.value_ = genttbar.Eta();
      gen_ttbar_mass_.value_ = genttbar.M();

      LV allGenObjects = *(topGenObjects.GenAntiLepton_) + *(topGenObjects.GenLepton_) + *(topGenObjects.GenMet_);
      for(int i=0;i<(int)topGenObjects.allGenJets_->size();i++)allGenObjects = allGenObjects + topGenObjects.allGenJets_->at(i);
      gen_all_mass_.value_ = allGenObjects.M();

      gen_r_mass_.value_ = gen_ttbar_mass_.value_/gen_all_mass_.value_;

      gen_ttbar_delta_eta_.value_ =fabs(gentop.Eta()-genantitop.Eta());
      gen_ttbar_delta_rapidity_.value_ =fabs(gentop.Rapidity()-genantitop.Rapidity());
      gen_ttbar_delta_phi_.value_ =fabs(gentop.DeltaPhi(genantitop));

      double gen_x1 = (gentop.E()+genantitop.E()+gentop.Pz()+genantitop.Pz())/(2*protonE);
      double gen_x2 = (gentop.E()+genantitop.E()-gentop.Pz()-genantitop.Pz())/(2*protonE);
      gen_x1_.value_ = gen_x1;
      gen_x2_.value_ = gen_x2;

      //leptons                                                                                                                                                                               
      gen_l_pt_.value_ = genlep.Pt();
      gen_lbar_pt_.value_ = genantilep.Pt();
      gen_l_eta_.value_ = genlep.Eta();
      gen_l_mass_.value_ = genlep.M();
      gen_l_pdgid_.value_ = topGenObjects.GenLeptonPdgId_;

      gen_lbar_eta_.value_ = genantilep.Eta();
      gen_l_phi_.value_ = genlep.Phi();
      gen_lbar_phi_.value_ = genantilep.Phi();
      gen_lbar_mass_.value_ = genantilep.M();
      gen_lbar_pdgid_.value_ = topGenObjects.GenAntiLeptonPdgId_;

      gen_llbar_mass_.value_=genllbar.M();
      gen_llbar_pt_.value_=genllbar.Pt();
      gen_llbar_phi_.value_=genllbar.Phi();
      gen_llbar_rapidity_.value_=genllbar.Rapidity();
      //      gen_llbar_arapidity_.value_=fabs(genllbar.Rapidity());

      gen_llbar_delta_eta_.value_ =fabs(genlep.Eta()-genantilep.Eta());
      gen_llbar_delta_rapidity_.value_ =fabs(genlep.Rapidity()-genantilep.Rapidity());
      gen_llbar_delta_phi_.value_ =fabs(genlep.DeltaPhi(genantilep));

      gen_b_pt_.value_ = genb.Pt();
      gen_bbar_pt_.value_ = genantib.Pt();
      gen_b_eta_.value_ = genb.Eta();
      gen_bbar_eta_.value_ = genantib.Eta();
      gen_b_phi_.value_ = genb.Phi();
      gen_bbar_phi_.value_ = genantib.Phi();
      gen_b_mass_.value_ = genb.M();
      gen_bbar_mass_.value_ = genantib.M();

      //      gen_bbbar_mass_.value_=genbbbar.M();
      //      gen_bbbar_pt_.value_=genbbbar.Pt();
      //      gen_bbbar_phi_.value_=genbbbar.Phi();
      //      gen_bbbar_rapidity_.value_=genbbbar.Rapidity();
      //      gen_bbbar_arapidity_.value_=fabs(genbbbar.Rapidity());

      //      gen_bbbar_delta_eta_.value_ =fabs(genb.Eta()-genantib.Eta());
      //      gen_bbbar_delta_rapidity_.value_ =fabs(genb.Rapidity()-genantib.Rapidity());
      //      gen_bbbar_delta_phi_.value_ =fabs(genb.DeltaPhi(genantib));

      //Neutrino                                                                                                                                                                              
      gen_nu_pt_.value_ = gennu.Pt();
      gen_nubar_pt_.value_ = genantinu.Pt();
      gen_nu_phi_.value_ = gennu.Phi();
      gen_nubar_phi_.value_ = genantinu.Phi();
      gen_nu_eta_.value_ = gennu.Eta();
      gen_nubar_eta_.value_ = genantinu.Eta();
      gen_nu_mass_.value_ = gennu.M();
      gen_nubar_mass_.value_ = genantinu.M();

      //      gen_nunubar_mass_.value_=gennunubar.M();
      //      gen_nunubar_pt_.value_=gennunubar.Pt();
      //      gen_nunubar_phi_.value_=gennunubar.Phi();
      //      gen_nunubar_rapidity_.value_=gennunubar.Rapidity();
      //      gen_nunubar_arapidity_.value_=fabs(gennunubar.Rapidity());

      //      gen_nunubar_delta_eta_.value_ =fabs(gennu.Eta()-genantinu.Eta());
      //      gen_nunubar_delta_rapidity_.value_ =fabs(gennu.Rapidity()-genantinu.Rapidity());
      //      gen_nunubar_delta_phi_.value_ =fabs(gennu.DeltaPhi(genantinu));

      gen_MT2_.value_ = genMT2; 

      //Evaluate discr based on the genvars                                                                                                                                                   
      VariablesPhiTT *genvars = VariablesPhiTT::fillVariables(eventMetadata, recoObjects, commonGenObjects, topGenObjects, kinematicReconstructionSolutions, recoObjectIndices, genObjectIndices, genLevelWeights, recoLevelWeights, weight);
      gen_b1k_.value_ = genvars->gen_b1k.value_;
      gen_b2k_.value_ = genvars->gen_b2k.value_;
      gen_b1j_.value_ = genvars->gen_b1j.value_;
      gen_b2j_.value_ = genvars->gen_b2j.value_;
      gen_b1r_.value_ = genvars->gen_b1r.value_;
      gen_b2r_.value_ = genvars->gen_b2r.value_;
      gen_b1q_.value_ = genvars->gen_b1q.value_;
      gen_b2q_.value_ = genvars->gen_b2q.value_;
      gen_b1n_.value_ = genvars->gen_b1n.value_;
      gen_b2n_.value_ = genvars->gen_b2n.value_;

      //gen_b_Pkk_.value_ = genvars->gen_bP_kk.value_;
      //gen_b_Mkk_.value_ = genvars->gen_bM_kk.value_;
      //gen_b_Pjj_.value_ = genvars->gen_bP_jj.value_;
      //gen_b_Mjj_.value_ = genvars->gen_bM_jj.value_;
      //gen_b_Prr_.value_ = genvars->gen_bP_rr.value_;
      //gen_b_Mrr_.value_ = genvars->gen_bM_rr.value_;
      //gen_b_Pqq_.value_ = genvars->gen_bP_qq.value_;
      //gen_b_Mqq_.value_ = genvars->gen_bM_qq.value_;
      //gen_b_Pnn_.value_ = genvars->gen_bP_nn.value_;
      //gen_b_Mnn_.value_ = genvars->gen_bM_nn.value_;

      gen_c_kk_.value_ = genvars->gen_ckk.value_;
      gen_c_rr_.value_ = genvars->gen_crr.value_;
      gen_c_nn_.value_ = genvars->gen_cnn.value_;
      gen_c_kj_.value_ = genvars->gen_ckj.value_;
      gen_c_rq_.value_ = genvars->gen_crq.value_;

      gen_c_rk_.value_ = genvars->gen_crk.value_;
      gen_c_kr_.value_ = genvars->gen_ckr.value_;
      gen_c_nr_.value_ = genvars->gen_cnr.value_;
      gen_c_rn_.value_ = genvars->gen_crn.value_;
      gen_c_nk_.value_ = genvars->gen_cnk.value_;
      gen_c_kn_.value_ = genvars->gen_ckn.value_;

      gen_c_rj_.value_ = genvars->gen_crj.value_;
      gen_c_jr_.value_ = genvars->gen_cjr.value_;
      //gen_c_qk_.value_ = genvars->gen_cqk.value_;
      //gen_c_kq_.value_ = genvars->gen_ckq.value_;

      //gen_c_qj_.value_ = genvars->gen_cqj.value_;
      //gen_c_jq_.value_ = genvars->gen_cjq.value_;
      //gen_c_nq_.value_ = genvars->gen_cnq.value_;
      //gen_c_qn_.value_ = genvars->gen_cqn.value_;
      //gen_c_nj_.value_ = genvars->gen_cnj.value_;
      //gen_c_jn_.value_ = genvars->gen_cjn.value_;

      //gen_c_Prk_.value_ = genvars->gen_cP_rk.value_;
      //gen_c_Mrk_.value_ = genvars->gen_cM_rk.value_;
      //gen_c_Pnr_.value_ = genvars->gen_cP_nr.value_;
      //gen_c_Mnr_.value_ = genvars->gen_cM_nr.value_;
      //gen_c_Pnk_.value_ = genvars->gen_cP_nk.value_;
      //gen_c_Mnk_.value_ = genvars->gen_cM_nk.value_;

      //gen_c_Pqj_.value_ = genvars->gen_cqj.value_ + genvars->gen_cjq.value_;
      //gen_c_Mqj_.value_ = genvars->gen_cqj.value_ - genvars->gen_cjq.value_;
      //gen_c_Pnq_.value_ = genvars->gen_cnq.value_ + genvars->gen_cqn.value_;
      //gen_c_Mnq_.value_ = genvars->gen_cnq.value_ - genvars->gen_cqn.value_;
      //gen_c_Pnj_.value_ = genvars->gen_cnj.value_ + genvars->gen_cjn.value_;
      //gen_c_Mnj_.value_ = genvars->gen_cnj.value_ + genvars->gen_cjn.value_;

      //gen_c_han_.value_ = +genvars->gen_ckk.value_ - genvars->gen_crr.value_ - genvars->gen_cnn.value_;
      //gen_c_sca_.value_ = -genvars->gen_ckk.value_ + genvars->gen_crr.value_ - genvars->gen_cnn.value_;
      //gen_c_tra_.value_ = -genvars->gen_ckk.value_ - genvars->gen_crr.value_ + genvars->gen_cnn.value_;
      //gen_c_lij_.value_ = -genvars->gen_ckj.value_ - genvars->gen_crr.value_ - genvars->gen_cnn.value_;
      //gen_c_liq_.value_ = -genvars->gen_ckk.value_ - genvars->gen_crq.value_ - genvars->gen_cnn.value_;

      //gen_c_rkP_.value_ = -genvars->gen_crk.value_ - genvars->gen_crk.value_ - genvars->gen_cnn.value_;
      //gen_c_rkM_.value_ = -genvars->gen_crk.value_ + genvars->gen_crk.value_ - genvars->gen_cnn.value_;
      //gen_c_nrP_.value_ = -genvars->gen_cnr.value_ - genvars->gen_crn.value_ - genvars->gen_ckk.value_;
      //gen_c_nrM_.value_ = -genvars->gen_cnr.value_ + genvars->gen_crn.value_ - genvars->gen_ckk.value_;
      //gen_c_nkP_.value_ = -genvars->gen_cnk.value_ - genvars->gen_cnk.value_ - genvars->gen_crr.value_;
      //gen_c_nkM_.value_ = -genvars->gen_cnk.value_ + genvars->gen_cnk.value_ - genvars->gen_crr.value_;

      //gen_c_qjP_.value_ = -genvars->gen_cqj.value_ - genvars->gen_cqj.value_ - genvars->gen_cnn.value_;
      //gen_c_qjM_.value_ = -genvars->gen_cqj.value_ + genvars->gen_cqj.value_ - genvars->gen_cnn.value_;
      //gen_c_nqP_.value_ = -genvars->gen_cnq.value_ - genvars->gen_cqn.value_ - genvars->gen_ckj.value_;
      //gen_c_nqM_.value_ = -genvars->gen_cnq.value_ + genvars->gen_cqn.value_ - genvars->gen_ckj.value_;
      //gen_c_njP_.value_ = -genvars->gen_cnj.value_ - genvars->gen_cnj.value_ - genvars->gen_crq.value_;
      //gen_c_njM_.value_ = -genvars->gen_cnj.value_ + genvars->gen_cnj.value_ - genvars->gen_crq.value_;

      gen_ll_cHel_.value_ = genvars->gen_cHel.value_;
      gen_ll_cLab_.value_ = genvars->gen_cLab.value_;
      gen_ll_kNorm_.value_ = genvars->gen_kNorm.value_;
      gen_ll_rNorm_.value_ = genvars->gen_rNorm.value_;

      delete genvars;
    
   }
    
}


std::vector<VariablesBase*> VariablesSpinCorr::fillVariables(const EventMetadata& eventMetadata, 
                                                          const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                                          const TopGenObjects& topGenObjects,
                                                          const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                                                          const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                                                          const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                                                          const double& weight)
{
    std::vector<VariablesBase*> result;

        VariablesBase* variables = new VariablesSpinCorr(eventMetadata, recoObjects, commonGenObjects, topGenObjects, kinematicReconstructionSolutions, recoObjectIndices, genObjectIndices, genLevelWeights, recoLevelWeights,weight);
        result.push_back(variables);
    
    return result;
}

double VariablesSpinCorr::getMT2Variable(TLorentzVector lep, TLorentzVector antilep, TLorentzVector MET){
    Double_t l_1[3], l_2[3], emiss[3];

    l_1[0]=0.0;
    l_1[1]=lep.Px();
    l_1[2]=lep.Py();

    l_2[0]=0.0;
    l_2[1]=antilep.Px();
    l_2[2]=antilep.Py();

    emiss[0]=0.0;
    emiss[1]=MET.Px();
    emiss[2]=MET.Py();

    mt2_bisect::mt2 mt2_event;
    mt2_event.set_momenta(l_1, l_2, emiss);
    mt2_event.set_mn(0.0);
    return mt2_event.get_mt2();

}
