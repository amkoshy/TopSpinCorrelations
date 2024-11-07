#include <iostream>
#include <cstdlib>

#include <TTree.h>
#include <TString.h>
#include <Rtypes.h>

#include "TreeHandlerSpinCorr.h"
#include "VariablesBase.h"
#include "VariablesSpinCorr.h"
#include "analysisStructs.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/KinematicReconstructionSolution.h"





TreeHandlerSpinCorr::TreeHandlerSpinCorr(const char* inputDir,
                                             const std::vector<TString>& selectionStepsNoCategories):
TreeHandlerBase("ttBar_", inputDir, selectionStepsNoCategories, new VariablesSpinCorr())
{
    std::cout<<"--- Beginning setting up MVA tree handler for top jets assignment\n";
    std::cout<<"=== Finishing setting up MVA tree handler for top jets assignment\n\n";

}



void TreeHandlerSpinCorr::fillVariables(const EventMetadata& eventMetadata,
                                     const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                          const TopGenObjects& topGenObjects,
                                          const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                                          const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                                          const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                                          const double& weight, const TString&,
                                          std::vector<VariablesBase*>& variables)const
{
    // Loop over all jet combinations and get MVA input variables
    const std::vector<VariablesBase*> v_variablesSpinCorr = 
            VariablesSpinCorr::fillVariables(eventMetadata, recoObjects, commonGenObjects,topGenObjects,kinematicReconstructionSolutions, recoObjectIndices,  genObjectIndices,genLevelWeights,recoLevelWeights,weight);
    
    // Fill the MVA variables
    variables.insert(variables.end(), v_variablesSpinCorr.begin(), v_variablesSpinCorr.end());
}



void TreeHandlerSpinCorr::bookBranches(TTree* tree, VariablesBase* const variables_)const
{
    VariablesSpinCorr* const variablesSpinCorr = dynamic_cast<VariablesSpinCorr*>(variables_);
    this->createBranch(tree, variablesSpinCorr->top_pt_);
    this->createBranch(tree, variablesSpinCorr->top_phi_);
    this->createBranch(tree, variablesSpinCorr->top_rapidity_);
    //    this->createBranch(tree, variablesSpinCorr->top_arapidity_);
    this->createBranch(tree, variablesSpinCorr->top_eta_);
    this->createBranch(tree, variablesSpinCorr->top_mass_);

    this->createBranch(tree, variablesSpinCorr->tbar_pt_);
    this->createBranch(tree, variablesSpinCorr->tbar_phi_);
    this->createBranch(tree, variablesSpinCorr->tbar_rapidity_);
    //    this->createBranch(tree, variablesSpinCorr->tbar_arapidity_);
    this->createBranch(tree, variablesSpinCorr->tbar_eta_);
    this->createBranch(tree, variablesSpinCorr->tbar_mass_);

    this->createBranch(tree, variablesSpinCorr->l_pt_);
    this->createBranch(tree, variablesSpinCorr->l_eta_);
    this->createBranch(tree, variablesSpinCorr->l_phi_);
    this->createBranch(tree, variablesSpinCorr->l_mass_);
    this->createBranch(tree, variablesSpinCorr->l_pdgid_);

    this->createBranch(tree, variablesSpinCorr->lbar_pt_);
    this->createBranch(tree, variablesSpinCorr->lbar_eta_);
    this->createBranch(tree, variablesSpinCorr->lbar_phi_);
    this->createBranch(tree, variablesSpinCorr->lbar_mass_);
    this->createBranch(tree, variablesSpinCorr->lbar_pdgid_);

    //    this->createBranch(tree, variablesSpinCorr->e_pt_);
    //    this->createBranch(tree, variablesSpinCorr->e_eta_);
    //    this->createBranch(tree, variablesSpinCorr->e_phi_);
    //    this->createBranch(tree, variablesSpinCorr->e_mass_);

    //    this->createBranch(tree, variablesSpinCorr->ebar_pt_);
    //    this->createBranch(tree, variablesSpinCorr->ebar_eta_);
    //    this->createBranch(tree, variablesSpinCorr->ebar_phi_);
    //    this->createBranch(tree, variablesSpinCorr->ebar_mass_);

    //    this->createBranch(tree, variablesSpinCorr->mu_pt_);
    //    this->createBranch(tree, variablesSpinCorr->mu_eta_);
    //    this->createBranch(tree, variablesSpinCorr->mu_phi_);
    //    this->createBranch(tree, variablesSpinCorr->mu_mass_);

    //    this->createBranch(tree, variablesSpinCorr->mubar_pt_);
    //    this->createBranch(tree, variablesSpinCorr->mubar_eta_);
    //    this->createBranch(tree, variablesSpinCorr->mubar_phi_);
    //    this->createBranch(tree, variablesSpinCorr->mubar_mass_);

    this->createBranch(tree, variablesSpinCorr->b_pt_);
    this->createBranch(tree, variablesSpinCorr->b_eta_);
    this->createBranch(tree, variablesSpinCorr->b_phi_);
    this->createBranch(tree, variablesSpinCorr->b_mass_);

    this->createBranch(tree, variablesSpinCorr->bbar_pt_);
    this->createBranch(tree, variablesSpinCorr->bbar_eta_);
    this->createBranch(tree, variablesSpinCorr->bbar_phi_);
    this->createBranch(tree, variablesSpinCorr->bbar_mass_);

    this->createBranch(tree, variablesSpinCorr->nu_pt_);
    this->createBranch(tree, variablesSpinCorr->nu_eta_);
    this->createBranch(tree, variablesSpinCorr->nu_phi_);
    this->createBranch(tree, variablesSpinCorr->nu_mass_);

    this->createBranch(tree, variablesSpinCorr->nubar_pt_);
    this->createBranch(tree, variablesSpinCorr->nubar_eta_);
    this->createBranch(tree, variablesSpinCorr->nubar_phi_);
    this->createBranch(tree, variablesSpinCorr->nubar_mass_);

    this->createBranch(tree, variablesSpinCorr->met_pt_);
    this->createBranch(tree, variablesSpinCorr->met_phi_);
    this->createBranch(tree, variablesSpinCorr->met_mass_);

    this->createBranch(tree, variablesSpinCorr->ttbar_pt_);
    this->createBranch(tree, variablesSpinCorr->ttbar_phi_);
    this->createBranch(tree, variablesSpinCorr->ttbar_rapidity_);
    //    this->createBranch(tree, variablesSpinCorr->ttbar_arapidity_);
    this->createBranch(tree, variablesSpinCorr->ttbar_eta_);
    this->createBranch(tree, variablesSpinCorr->ttbar_delta_phi_);
    this->createBranch(tree, variablesSpinCorr->ttbar_delta_eta_);
    this->createBranch(tree, variablesSpinCorr->ttbar_delta_rapidity_);
    this->createBranch(tree, variablesSpinCorr->ttbar_mass_);

    this->createBranch(tree, variablesSpinCorr->llbar_pt_);
    this->createBranch(tree, variablesSpinCorr->llbar_phi_);
    this->createBranch(tree, variablesSpinCorr->llbar_rapidity_);
    //    this->createBranch(tree, variablesSpinCorr->llbar_arapidity_);
    this->createBranch(tree, variablesSpinCorr->llbar_delta_phi_);
    this->createBranch(tree, variablesSpinCorr->llbar_delta_eta_);
    this->createBranch(tree, variablesSpinCorr->llbar_delta_rapidity_);
    this->createBranch(tree, variablesSpinCorr->llbar_mass_);

    //    this->createBranch(tree, variablesSpinCorr->bbbar_pt_);
    //    this->createBranch(tree, variablesSpinCorr->bbbar_phi_);
    //    this->createBranch(tree, variablesSpinCorr->bbbar_rapidity_);
    //    this->createBranch(tree, variablesSpinCorr->bbbar_arapidity_);
    //    this->createBranch(tree, variablesSpinCorr->bbbar_delta_phi_);
    //    this->createBranch(tree, variablesSpinCorr->bbbar_delta_eta_);
    //    this->createBranch(tree, variablesSpinCorr->bbbar_delta_rapidity_);
    //    this->createBranch(tree, variablesSpinCorr->bbbar_mass_);

    //    this->createBranch(tree, variablesSpinCorr->nunubar_pt_);
    //    this->createBranch(tree, variablesSpinCorr->nunubar_phi_);
    //    this->createBranch(tree, variablesSpinCorr->nunubar_rapidity_);
    //    this->createBranch(tree, variablesSpinCorr->nunubar_arapidity_);
    //    this->createBranch(tree, variablesSpinCorr->nunubar_delta_phi_);
    //    this->createBranch(tree, variablesSpinCorr->nunubar_delta_eta_);
    //    this->createBranch(tree, variablesSpinCorr->nunubar_delta_rapidity_);
    //    this->createBranch(tree, variablesSpinCorr->nunubar_mass_);

    this->createBranch(tree, variablesSpinCorr->MT2_);

    this->createBranch(tree, variablesSpinCorr->all_mass_);
    this->createBranch(tree, variablesSpinCorr->r_mass_);
    this->createBranch(tree, variablesSpinCorr->jet_multiplicity_);
    this->createBranch(tree, variablesSpinCorr->bjet_multiplicity_);
    this->createBranch(tree, variablesSpinCorr->x1_);
    this->createBranch(tree, variablesSpinCorr->x2_);

    //ExtraJet variables
    this->createBranch(tree, variablesSpinCorr->n_extraJets_iso08_);

    //ttbar frame variables
    this->createBranch(tree, variablesSpinCorr->top_scatteringangle_ttbarframe_);

    //    this->createBranch(tree, variablesSpinCorr->top_pt_ttbarframe_);
    //    this->createBranch(tree, variablesSpinCorr->top_phi_ttbarframe_);
    //    this->createBranch(tree, variablesSpinCorr->top_rapidity_ttbarframe_);
    //    this->createBranch(tree, variablesSpinCorr->top_arapidity_ttbarframe_);
    //    this->createBranch(tree, variablesSpinCorr->tbar_pt_ttbarframe_);
    //    this->createBranch(tree, variablesSpinCorr->tbar_phi_ttbarframe_);
    //    this->createBranch(tree, variablesSpinCorr->tbar_rapidity_ttbarframe_);
    //    this->createBranch(tree, variablesSpinCorr->tbar_arapidity_ttbarframe_);

    //top and antitop frame variables
    //    this->createBranch(tree, variablesSpinCorr->l_pt_antitopframe_);
    //    this->createBranch(tree, variablesSpinCorr->l_eta_antitopframe_);
    //    this->createBranch(tree, variablesSpinCorr->l_phi_antitopframe_);
    //    this->createBranch(tree, variablesSpinCorr->lbar_pt_topframe_);
    //    this->createBranch(tree, variablesSpinCorr->lbar_eta_topframe_);
    //    this->createBranch(tree, variablesSpinCorr->lbar_phi_topframe_);

    //Spin correlation variables
    this->createBranch(tree, variablesSpinCorr->b1k_);
    this->createBranch(tree, variablesSpinCorr->b2k_);
    this->createBranch(tree, variablesSpinCorr->b1j_);
    this->createBranch(tree, variablesSpinCorr->b2j_);
    this->createBranch(tree, variablesSpinCorr->b1r_);
    this->createBranch(tree, variablesSpinCorr->b2r_);
    this->createBranch(tree, variablesSpinCorr->b1q_);
    this->createBranch(tree, variablesSpinCorr->b2q_);
    this->createBranch(tree, variablesSpinCorr->b1n_);
    this->createBranch(tree, variablesSpinCorr->b2n_);

    //this->createBranch(tree, variablesSpinCorr->b_Pkk_);
    //this->createBranch(tree, variablesSpinCorr->b_Mkk_);
    //this->createBranch(tree, variablesSpinCorr->b_Pjj_);
    //this->createBranch(tree, variablesSpinCorr->b_Mjj_);
    //this->createBranch(tree, variablesSpinCorr->b_Prr_);
    //this->createBranch(tree, variablesSpinCorr->b_Mrr_);
    //this->createBranch(tree, variablesSpinCorr->b_Pqq_);
    //this->createBranch(tree, variablesSpinCorr->b_Mqq_);
    //this->createBranch(tree, variablesSpinCorr->b_Pnn_);
    //this->createBranch(tree, variablesSpinCorr->b_Mnn_);

    this->createBranch(tree, variablesSpinCorr->c_kk_);
    this->createBranch(tree, variablesSpinCorr->c_rr_);
    this->createBranch(tree, variablesSpinCorr->c_nn_);
    this->createBranch(tree, variablesSpinCorr->c_kj_);
    this->createBranch(tree, variablesSpinCorr->c_rq_);
 
    this->createBranch(tree, variablesSpinCorr->c_rk_);
    this->createBranch(tree, variablesSpinCorr->c_kr_);
    this->createBranch(tree, variablesSpinCorr->c_nr_);
    this->createBranch(tree, variablesSpinCorr->c_rn_);
    this->createBranch(tree, variablesSpinCorr->c_nk_);
    this->createBranch(tree, variablesSpinCorr->c_kn_);

    this->createBranch(tree, variablesSpinCorr->c_rj_);
    this->createBranch(tree, variablesSpinCorr->c_jr_);
    //this->createBranch(tree, variablesSpinCorr->c_qk_);
    //this->createBranch(tree, variablesSpinCorr->c_kq_);

    //this->createBranch(tree, variablesSpinCorr->c_qj_);
    //this->createBranch(tree, variablesSpinCorr->c_jq_);
    //this->createBranch(tree, variablesSpinCorr->c_nq_);
    //this->createBranch(tree, variablesSpinCorr->c_qn_);
    //this->createBranch(tree, variablesSpinCorr->c_nj_);
    //this->createBranch(tree, variablesSpinCorr->c_jn_);

    //this->createBranch(tree, variablesSpinCorr->c_han_);
    //this->createBranch(tree, variablesSpinCorr->c_sca_);
    //this->createBranch(tree, variablesSpinCorr->c_tra_);
    //this->createBranch(tree, variablesSpinCorr->c_lij_);
    //this->createBranch(tree, variablesSpinCorr->c_liq_);

    //this->createBranch(tree, variablesSpinCorr->c_rkP_);
    //this->createBranch(tree, variablesSpinCorr->c_rkM_);
    //this->createBranch(tree, variablesSpinCorr->c_nrP_);
    //this->createBranch(tree, variablesSpinCorr->c_nrM_);
    //this->createBranch(tree, variablesSpinCorr->c_nkP_);
    //this->createBranch(tree, variablesSpinCorr->c_nkM_);

    //this->createBranch(tree, variablesSpinCorr->c_qjP_);
    //this->createBranch(tree, variablesSpinCorr->c_qjM_);
    //this->createBranch(tree, variablesSpinCorr->c_nqP_);
    //this->createBranch(tree, variablesSpinCorr->c_nqM_);
    //this->createBranch(tree, variablesSpinCorr->c_njP_);
    //this->createBranch(tree, variablesSpinCorr->c_njM_);

    this->createBranch(tree, variablesSpinCorr->ll_cHel_);
    this->createBranch(tree, variablesSpinCorr->ll_cLab_);
    this->createBranch(tree, variablesSpinCorr->ll_kNorm_);
    this->createBranch(tree, variablesSpinCorr->ll_rNorm_);

    this->createBranch(tree, variablesSpinCorr->gen_top_pt_);
    this->createBranch(tree, variablesSpinCorr->gen_top_phi_);
    this->createBranch(tree, variablesSpinCorr->gen_top_rapidity_);
    //    this->createBranch(tree, variablesSpinCorr->gen_top_arapidity_);
    this->createBranch(tree, variablesSpinCorr->gen_top_eta_);
    this->createBranch(tree, variablesSpinCorr->gen_top_mass_);

    this->createBranch(tree, variablesSpinCorr->gen_tbar_pt_);
    this->createBranch(tree, variablesSpinCorr->gen_tbar_phi_);
    this->createBranch(tree, variablesSpinCorr->gen_tbar_rapidity_);
    //    this->createBranch(tree, variablesSpinCorr->gen_tbar_arapidity_);
    this->createBranch(tree, variablesSpinCorr->gen_tbar_eta_);
    this->createBranch(tree, variablesSpinCorr->gen_tbar_mass_);

    this->createBranch(tree, variablesSpinCorr->gen_l_pt_);
    this->createBranch(tree, variablesSpinCorr->gen_l_eta_); 
    this->createBranch(tree, variablesSpinCorr->gen_l_phi_);
    this->createBranch(tree, variablesSpinCorr->gen_l_mass_);
    this->createBranch(tree, variablesSpinCorr->gen_l_pdgid_);
 
    this->createBranch(tree, variablesSpinCorr->gen_lbar_pt_);
    this->createBranch(tree, variablesSpinCorr->gen_lbar_eta_);
    this->createBranch(tree, variablesSpinCorr->gen_lbar_phi_);
    this->createBranch(tree, variablesSpinCorr->gen_lbar_mass_);
    this->createBranch(tree, variablesSpinCorr->gen_lbar_pdgid_);

    //    this->createBranch(tree, variablesSpinCorr->gen_e_pt_);
    //    this->createBranch(tree, variablesSpinCorr->gen_e_eta_);
    //    this->createBranch(tree, variablesSpinCorr->gen_e_phi_);
    //    this->createBranch(tree, variablesSpinCorr->gen_e_mass_);

    //    this->createBranch(tree, variablesSpinCorr->gen_ebar_pt_);
    //    this->createBranch(tree, variablesSpinCorr->gen_ebar_eta_);
    //    this->createBranch(tree, variablesSpinCorr->gen_ebar_phi_);
    //    this->createBranch(tree, variablesSpinCorr->gen_ebar_mass_);

    //    this->createBranch(tree, variablesSpinCorr->gen_mu_pt_);
    //    this->createBranch(tree, variablesSpinCorr->gen_mu_eta_);
    //    this->createBranch(tree, variablesSpinCorr->gen_mu_phi_);
    //    this->createBranch(tree, variablesSpinCorr->gen_mu_mass_);

    //    this->createBranch(tree, variablesSpinCorr->gen_mubar_pt_);
    //    this->createBranch(tree, variablesSpinCorr->gen_mubar_eta_);
    //    this->createBranch(tree, variablesSpinCorr->gen_mubar_phi_);
    //    this->createBranch(tree, variablesSpinCorr->gen_mubar_mass_);

    this->createBranch(tree, variablesSpinCorr->gen_b_pt_);
    this->createBranch(tree, variablesSpinCorr->gen_b_eta_);
    this->createBranch(tree, variablesSpinCorr->gen_b_phi_);
    this->createBranch(tree, variablesSpinCorr->gen_b_mass_);

    this->createBranch(tree, variablesSpinCorr->gen_bbar_pt_);
    this->createBranch(tree, variablesSpinCorr->gen_bbar_eta_);
    this->createBranch(tree, variablesSpinCorr->gen_bbar_phi_);
    this->createBranch(tree, variablesSpinCorr->gen_bbar_mass_);

    this->createBranch(tree, variablesSpinCorr->gen_nu_pt_);
    this->createBranch(tree, variablesSpinCorr->gen_nu_eta_);
    this->createBranch(tree, variablesSpinCorr->gen_nu_phi_);
    this->createBranch(tree, variablesSpinCorr->gen_nu_mass_);

    this->createBranch(tree, variablesSpinCorr->gen_nubar_pt_);
    this->createBranch(tree, variablesSpinCorr->gen_nubar_eta_);
    this->createBranch(tree, variablesSpinCorr->gen_nubar_phi_);
    this->createBranch(tree, variablesSpinCorr->gen_nubar_mass_);

    this->createBranch(tree, variablesSpinCorr->gen_ttbar_pt_);
    this->createBranch(tree, variablesSpinCorr->gen_ttbar_phi_);
    this->createBranch(tree, variablesSpinCorr->gen_ttbar_rapidity_);
    //    this->createBranch(tree, variablesSpinCorr->gen_ttbar_arapidity_);
    this->createBranch(tree, variablesSpinCorr->gen_ttbar_eta_);
    this->createBranch(tree, variablesSpinCorr->gen_ttbar_delta_phi_);
    this->createBranch(tree, variablesSpinCorr->gen_ttbar_delta_eta_);
    this->createBranch(tree, variablesSpinCorr->gen_ttbar_delta_rapidity_);
    this->createBranch(tree, variablesSpinCorr->gen_ttbar_mass_);

    this->createBranch(tree, variablesSpinCorr->gen_llbar_pt_);
    this->createBranch(tree, variablesSpinCorr->gen_llbar_phi_);
    this->createBranch(tree, variablesSpinCorr->gen_llbar_rapidity_);
    //    this->createBranch(tree, variablesSpinCorr->gen_llbar_arapidity_);
    this->createBranch(tree, variablesSpinCorr->gen_llbar_delta_phi_);
    this->createBranch(tree, variablesSpinCorr->gen_llbar_delta_eta_);
    this->createBranch(tree, variablesSpinCorr->gen_llbar_delta_rapidity_);
    this->createBranch(tree, variablesSpinCorr->gen_llbar_mass_);

    //    this->createBranch(tree, variablesSpinCorr->gen_bbbar_pt_);
    //    this->createBranch(tree, variablesSpinCorr->gen_bbbar_phi_);
    //    this->createBranch(tree, variablesSpinCorr->gen_bbbar_rapidity_);
    //    this->createBranch(tree, variablesSpinCorr->gen_bbbar_arapidity_);
    //    this->createBranch(tree, variablesSpinCorr->gen_bbbar_delta_phi_);
    //    this->createBranch(tree, variablesSpinCorr->gen_bbbar_delta_eta_);
    //    this->createBranch(tree, variablesSpinCorr->gen_bbbar_delta_rapidity_);
    //    this->createBranch(tree, variablesSpinCorr->gen_bbbar_mass_);

    //    this->createBranch(tree, variablesSpinCorr->gen_nunubar_pt_);
    //    this->createBranch(tree, variablesSpinCorr->gen_nunubar_phi_);
    //    this->createBranch(tree, variablesSpinCorr->gen_nunubar_rapidity_);
    //    this->createBranch(tree, variablesSpinCorr->gen_nunubar_arapidity_);
    //    this->createBranch(tree, variablesSpinCorr->gen_nunubar_delta_phi_);
    //    this->createBranch(tree, variablesSpinCorr->gen_nunubar_delta_eta_);
    //    this->createBranch(tree, variablesSpinCorr->gen_nunubar_delta_rapidity_);
    //    this->createBranch(tree, variablesSpinCorr->gen_nunubar_mass_);

    this->createBranch(tree, variablesSpinCorr->gen_MT2_);

    this->createBranch(tree, variablesSpinCorr->gen_all_mass_);
    this->createBranch(tree, variablesSpinCorr->gen_r_mass_);
    this->createBranch(tree, variablesSpinCorr->gen_jet_multiplicity_);
    this->createBranch(tree, variablesSpinCorr->gen_x1_);
    this->createBranch(tree, variablesSpinCorr->gen_x2_);
    this->createBranch(tree, variablesSpinCorr->gen_production_mode_);

    //extraJet variables
    this->createBranch(tree, variablesSpinCorr->gen_n_extraJets_iso08_);

    //ttbar frame variables
    this->createBranch(tree, variablesSpinCorr->gen_top_scatteringangle_ttbarframe_);

    //    this->createBranch(tree, variablesSpinCorr->gen_top_pt_ttbarframe_);
    //    this->createBranch(tree, variablesSpinCorr->gen_top_phi_ttbarframe_);
    //    this->createBranch(tree, variablesSpinCorr->gen_top_rapidity_ttbarframe_);
    //    this->createBranch(tree, variablesSpinCorr->gen_top_arapidity_ttbarframe_);
    //    this->createBranch(tree, variablesSpinCorr->gen_tbar_pt_ttbarframe_);
    //    this->createBranch(tree, variablesSpinCorr->gen_tbar_phi_ttbarframe_);
    //    this->createBranch(tree, variablesSpinCorr->gen_tbar_rapidity_ttbarframe_);
    //    this->createBranch(tree, variablesSpinCorr->gen_tbar_arapidity_ttbarframe_);

    //top and antitop frame variables
    //    this->createBranch(tree, variablesSpinCorr->gen_l_pt_antitopframe_);
    //    this->createBranch(tree, variablesSpinCorr->gen_l_eta_antitopframe_);
    //    this->createBranch(tree, variablesSpinCorr->gen_l_phi_antitopframe_);
    //    this->createBranch(tree, variablesSpinCorr->gen_lbar_pt_topframe_);
    //    this->createBranch(tree, variablesSpinCorr->gen_lbar_eta_topframe_);
    //    this->createBranch(tree, variablesSpinCorr->gen_lbar_phi_topframe_);

    //Spin correlation variables                                                                                                                                                          $
    this->createBranch(tree, variablesSpinCorr->gen_b1k_);
    this->createBranch(tree, variablesSpinCorr->gen_b2k_);
    this->createBranch(tree, variablesSpinCorr->gen_b1j_);
    this->createBranch(tree, variablesSpinCorr->gen_b2j_);
    this->createBranch(tree, variablesSpinCorr->gen_b1r_);
    this->createBranch(tree, variablesSpinCorr->gen_b2r_);
    this->createBranch(tree, variablesSpinCorr->gen_b1q_);
    this->createBranch(tree, variablesSpinCorr->gen_b2q_);
    this->createBranch(tree, variablesSpinCorr->gen_b1n_);
    this->createBranch(tree, variablesSpinCorr->gen_b2n_);

    //this->createBranch(tree, variablesSpinCorr->gen_b_Pkk_);
    //this->createBranch(tree, variablesSpinCorr->gen_b_Mkk_);
    //this->createBranch(tree, variablesSpinCorr->gen_b_Pjj_);
    //this->createBranch(tree, variablesSpinCorr->gen_b_Mjj_);
    //this->createBranch(tree, variablesSpinCorr->gen_b_Prr_);
    //this->createBranch(tree, variablesSpinCorr->gen_b_Mrr_);
    //this->createBranch(tree, variablesSpinCorr->gen_b_Pqq_);
    //this->createBranch(tree, variablesSpinCorr->gen_b_Mqq_);
    //this->createBranch(tree, variablesSpinCorr->gen_b_Pnn_);
    //this->createBranch(tree, variablesSpinCorr->gen_b_Mnn_);

    this->createBranch(tree, variablesSpinCorr->gen_c_kk_);
    this->createBranch(tree, variablesSpinCorr->gen_c_rr_);
    this->createBranch(tree, variablesSpinCorr->gen_c_nn_);
    this->createBranch(tree, variablesSpinCorr->gen_c_kj_);
    this->createBranch(tree, variablesSpinCorr->gen_c_rq_);

    this->createBranch(tree, variablesSpinCorr->gen_c_rk_);
    this->createBranch(tree, variablesSpinCorr->gen_c_kr_);
    this->createBranch(tree, variablesSpinCorr->gen_c_nr_);
    this->createBranch(tree, variablesSpinCorr->gen_c_rn_);
    this->createBranch(tree, variablesSpinCorr->gen_c_nk_);
    this->createBranch(tree, variablesSpinCorr->gen_c_kn_);

    this->createBranch(tree, variablesSpinCorr->gen_c_rj_);
    this->createBranch(tree, variablesSpinCorr->gen_c_jr_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_qk_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_kq_);

    //this->createBranch(tree, variablesSpinCorr->gen_c_qj_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_jq_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_nq_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_qn_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_nj_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_jn_);

    //this->createBranch(tree, variablesSpinCorr->gen_c_Prk_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_Mrk_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_Pnr_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_Mnr_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_Pnk_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_Mnk_);

    //this->createBranch(tree, variablesSpinCorr->gen_c_Pqj_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_Mqj_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_Pnq_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_Mnq_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_Pnj_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_Mnj_);

    //this->createBranch(tree, variablesSpinCorr->gen_c_han_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_sca_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_tra_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_lij_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_liq_);

    //this->createBranch(tree, variablesSpinCorr->gen_c_rkP_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_rkM_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_nrP_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_nrM_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_nkP_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_nkM_);

    //this->createBranch(tree, variablesSpinCorr->gen_c_qjP_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_qjM_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_nqP_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_nqM_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_njP_);
    //this->createBranch(tree, variablesSpinCorr->gen_c_njM_);

    this->createBranch(tree, variablesSpinCorr->gen_ll_cHel_);
    this->createBranch(tree, variablesSpinCorr->gen_ll_cLab_);
    this->createBranch(tree, variablesSpinCorr->gen_ll_kNorm_);
    this->createBranch(tree, variablesSpinCorr->gen_ll_rNorm_);

    this->createBranch(tree, variablesSpinCorr->entry_);
    this->createBranch(tree, variablesSpinCorr->isTopGen_);
    this->createBranch(tree, variablesSpinCorr->TopDecayMode_);
    this->createBranch(tree, variablesSpinCorr->isKinReco_);
    this->createBranch(tree, variablesSpinCorr->eventWeight_);
    this->createBranch(tree, variablesSpinCorr->trueLevelWeight_);
}



void TreeHandlerSpinCorr::fillBranches(TTree* tree, const std::vector<VariablesBase*>& v_variables)
{
    
    for(const VariablesBase* variablesTmp : v_variables){
        const VariablesSpinCorr* variablesSpinCorrTmp = dynamic_cast<const VariablesSpinCorr*>(variablesTmp);
        if(!variablesSpinCorrTmp){
            std::cerr<<"ERROR in TreeHandlerSpinCorr::fillBranches()! variables are of wrong type, cannot typecast\n"
                     <<"...break\n"<<std::endl;
            exit(395);
        }
        
        *(dynamic_cast<VariablesSpinCorr*>(variables_)) = *variablesSpinCorrTmp;
        
        tree->Fill();
    }
}



void TreeHandlerSpinCorr::importBranches(TTree* tree, std::vector<VariablesBase*>& v_variables)const
{
    // Set up variables struct
    VariablesSpinCorr variablesSpinCorr;


    // Set branch addresses
    this->importBranch(tree, variablesSpinCorr.top_pt_);
    this->importBranch(tree, variablesSpinCorr.top_phi_);
    this->importBranch(tree, variablesSpinCorr.top_rapidity_);
    //    this->importBranch(tree, variablesSpinCorr.top_arapidity_);
    this->importBranch(tree, variablesSpinCorr.top_eta_);
    this->importBranch(tree, variablesSpinCorr.top_mass_);

    this->importBranch(tree, variablesSpinCorr.tbar_pt_);
    this->importBranch(tree, variablesSpinCorr.tbar_phi_);
    this->importBranch(tree, variablesSpinCorr.tbar_rapidity_);
    //    this->importBranch(tree, variablesSpinCorr.tbar_arapidity_);
    this->importBranch(tree, variablesSpinCorr.tbar_eta_);
    this->importBranch(tree, variablesSpinCorr.tbar_mass_);

    this->importBranch(tree, variablesSpinCorr.l_pt_);
    this->importBranch(tree, variablesSpinCorr.l_eta_);
    this->importBranch(tree, variablesSpinCorr.l_phi_);
    this->importBranch(tree, variablesSpinCorr.l_mass_);
    this->importBranch(tree, variablesSpinCorr.l_pdgid_);

    this->importBranch(tree, variablesSpinCorr.lbar_pt_);
    this->importBranch(tree, variablesSpinCorr.lbar_eta_);
    this->importBranch(tree, variablesSpinCorr.lbar_phi_);
    this->importBranch(tree, variablesSpinCorr.lbar_mass_);
    this->importBranch(tree, variablesSpinCorr.lbar_pdgid_);

    //    this->importBranch(tree, variablesSpinCorr.e_pt_);
    //    this->importBranch(tree, variablesSpinCorr.e_eta_);
    //    this->importBranch(tree, variablesSpinCorr.e_phi_);
    //    this->importBranch(tree, variablesSpinCorr.e_mass_);

    //    this->importBranch(tree, variablesSpinCorr.ebar_pt_);
    //    this->importBranch(tree, variablesSpinCorr.ebar_eta_);
    //    this->importBranch(tree, variablesSpinCorr.ebar_phi_);
    //    this->importBranch(tree, variablesSpinCorr.ebar_mass_);

    //    this->importBranch(tree, variablesSpinCorr.mu_pt_);
    //    this->importBranch(tree, variablesSpinCorr.mu_eta_);
    //    this->importBranch(tree, variablesSpinCorr.mu_phi_);
    //    this->importBranch(tree, variablesSpinCorr.mu_mass_);

    //    this->importBranch(tree, variablesSpinCorr.mubar_pt_);
    //    this->importBranch(tree, variablesSpinCorr.mubar_eta_);    
    //    this->importBranch(tree, variablesSpinCorr.mubar_phi_);
    //    this->importBranch(tree, variablesSpinCorr.mubar_mass_);

    this->importBranch(tree, variablesSpinCorr.b_pt_);
    this->importBranch(tree, variablesSpinCorr.b_eta_);
    this->importBranch(tree, variablesSpinCorr.b_phi_);
    this->importBranch(tree, variablesSpinCorr.b_mass_);

    this->importBranch(tree, variablesSpinCorr.bbar_pt_);
    this->importBranch(tree, variablesSpinCorr.bbar_eta_);
    this->importBranch(tree, variablesSpinCorr.bbar_phi_);
    this->importBranch(tree, variablesSpinCorr.bbar_mass_);

    this->importBranch(tree, variablesSpinCorr.nu_pt_);
    this->importBranch(tree, variablesSpinCorr.nu_eta_);
    this->importBranch(tree, variablesSpinCorr.nu_phi_);
    this->importBranch(tree, variablesSpinCorr.nu_mass_);

    this->importBranch(tree, variablesSpinCorr.nubar_pt_);
    this->importBranch(tree, variablesSpinCorr.nubar_eta_);
    this->importBranch(tree, variablesSpinCorr.nubar_phi_);
    this->importBranch(tree, variablesSpinCorr.nubar_mass_);

    this->importBranch(tree, variablesSpinCorr.met_pt_);
    this->importBranch(tree, variablesSpinCorr.met_phi_);
    this->importBranch(tree, variablesSpinCorr.met_mass_);

    this->importBranch(tree, variablesSpinCorr.ttbar_pt_);
    this->importBranch(tree, variablesSpinCorr.ttbar_phi_);
    this->importBranch(tree, variablesSpinCorr.ttbar_rapidity_);
    //    this->importBranch(tree, variablesSpinCorr.ttbar_arapidity_);
    this->importBranch(tree, variablesSpinCorr.ttbar_eta_);
    this->importBranch(tree, variablesSpinCorr.ttbar_delta_phi_);
    this->importBranch(tree, variablesSpinCorr.ttbar_delta_eta_);
    this->importBranch(tree, variablesSpinCorr.ttbar_delta_rapidity_);
    this->importBranch(tree, variablesSpinCorr.ttbar_mass_);

    this->importBranch(tree, variablesSpinCorr.llbar_pt_);
    this->importBranch(tree, variablesSpinCorr.llbar_phi_);
    this->importBranch(tree, variablesSpinCorr.llbar_rapidity_);
    //    this->importBranch(tree, variablesSpinCorr.llbar_arapidity_);
    this->importBranch(tree, variablesSpinCorr.llbar_delta_phi_);
    this->importBranch(tree, variablesSpinCorr.llbar_delta_eta_);
    this->importBranch(tree, variablesSpinCorr.llbar_delta_rapidity_);
    this->importBranch(tree, variablesSpinCorr.llbar_mass_);

    //    this->importBranch(tree, variablesSpinCorr.bbbar_pt_);
    //    this->importBranch(tree, variablesSpinCorr.bbbar_phi_);
    //    this->importBranch(tree, variablesSpinCorr.bbbar_rapidity_);
    //    this->importBranch(tree, variablesSpinCorr.bbbar_arapidity_);
    //    this->importBranch(tree, variablesSpinCorr.bbbar_delta_phi_);
    //    this->importBranch(tree, variablesSpinCorr.bbbar_delta_eta_);
    //    this->importBranch(tree, variablesSpinCorr.bbbar_delta_rapidity_);
    //    this->importBranch(tree, variablesSpinCorr.bbbar_mass_);

    //    this->importBranch(tree, variablesSpinCorr.nunubar_pt_);
    //    this->importBranch(tree, variablesSpinCorr.nunubar_phi_);
    //    this->importBranch(tree, variablesSpinCorr.nunubar_rapidity_);
    //    this->importBranch(tree, variablesSpinCorr.nunubar_arapidity_);
    //    this->importBranch(tree, variablesSpinCorr.nunubar_delta_phi_);
    //    this->importBranch(tree, variablesSpinCorr.nunubar_delta_eta_);
    //    this->importBranch(tree, variablesSpinCorr.nunubar_delta_rapidity_);
    //    this->importBranch(tree, variablesSpinCorr.nunubar_mass_);

    this->importBranch(tree, variablesSpinCorr.MT2_);

    this->importBranch(tree, variablesSpinCorr.all_mass_);
    this->importBranch(tree, variablesSpinCorr.r_mass_);
    this->importBranch(tree, variablesSpinCorr.jet_multiplicity_);
    this->importBranch(tree, variablesSpinCorr.bjet_multiplicity_);
    this->importBranch(tree, variablesSpinCorr.x1_);
    this->importBranch(tree, variablesSpinCorr.x2_);

    //extraJet variables
    this->importBranch(tree, variablesSpinCorr.n_extraJets_iso08_);

    //ttbar frame variables
    this->importBranch(tree, variablesSpinCorr.top_scatteringangle_ttbarframe_);

    //    this->importBranch(tree, variablesSpinCorr.top_pt_ttbarframe_);
    //    this->importBranch(tree, variablesSpinCorr.top_phi_ttbarframe_);
    //    this->importBranch(tree, variablesSpinCorr.top_rapidity_ttbarframe_);
    //    this->importBranch(tree, variablesSpinCorr.top_arapidity_ttbarframe_);
    //    this->importBranch(tree, variablesSpinCorr.tbar_pt_ttbarframe_);
    //    this->importBranch(tree, variablesSpinCorr.tbar_phi_ttbarframe_);
    //    this->importBranch(tree, variablesSpinCorr.tbar_rapidity_ttbarframe_);
    //    this->importBranch(tree, variablesSpinCorr.tbar_arapidity_ttbarframe_);

    //top and antitop frame variables
    //    this->importBranch(tree, variablesSpinCorr.l_pt_antitopframe_);
    //    this->importBranch(tree, variablesSpinCorr.l_eta_antitopframe_);
    //    this->importBranch(tree, variablesSpinCorr.l_phi_antitopframe_);
    //    this->importBranch(tree, variablesSpinCorr.lbar_pt_topframe_);
    //    this->importBranch(tree, variablesSpinCorr.lbar_eta_topframe_);
    //    this->importBranch(tree, variablesSpinCorr.lbar_phi_topframe_);

    //Spin correlation variables
    this->importBranch(tree, variablesSpinCorr.b1k_);
    this->importBranch(tree, variablesSpinCorr.b2k_);
    this->importBranch(tree, variablesSpinCorr.b1j_);
    this->importBranch(tree, variablesSpinCorr.b2j_);
    this->importBranch(tree, variablesSpinCorr.b1r_);
    this->importBranch(tree, variablesSpinCorr.b2r_);
    this->importBranch(tree, variablesSpinCorr.b1q_);
    this->importBranch(tree, variablesSpinCorr.b2q_);
    this->importBranch(tree, variablesSpinCorr.b1n_);
    this->importBranch(tree, variablesSpinCorr.b2n_);
 
    //this->importBranch(tree, variablesSpinCorr.b_Pkk_);
    //this->importBranch(tree, variablesSpinCorr.b_Mkk_);
    //this->importBranch(tree, variablesSpinCorr.b_Pjj_);
    //this->importBranch(tree, variablesSpinCorr.b_Mjj_);
    //this->importBranch(tree, variablesSpinCorr.b_Prr_);
    //this->importBranch(tree, variablesSpinCorr.b_Mrr_);
    //this->importBranch(tree, variablesSpinCorr.b_Pqq_);
    //this->importBranch(tree, variablesSpinCorr.b_Mqq_);
    //this->importBranch(tree, variablesSpinCorr.b_Pnn_);
    //this->importBranch(tree, variablesSpinCorr.b_Mnn_);

    this->importBranch(tree, variablesSpinCorr.c_kk_);
    this->importBranch(tree, variablesSpinCorr.c_rr_);
    this->importBranch(tree, variablesSpinCorr.c_nn_);
    this->importBranch(tree, variablesSpinCorr.c_kj_);
    this->importBranch(tree, variablesSpinCorr.c_rq_);

    this->importBranch(tree, variablesSpinCorr.c_rk_);
    this->importBranch(tree, variablesSpinCorr.c_kr_);
    this->importBranch(tree, variablesSpinCorr.c_nr_);
    this->importBranch(tree, variablesSpinCorr.c_rn_);
    this->importBranch(tree, variablesSpinCorr.c_nk_);
    this->importBranch(tree, variablesSpinCorr.c_kn_);

    this->importBranch(tree, variablesSpinCorr.c_rj_);
    this->importBranch(tree, variablesSpinCorr.c_jr_);
    //this->importBranch(tree, variablesSpinCorr.c_qk_);
    //this->importBranch(tree, variablesSpinCorr.c_kq_);

    //this->importBranch(tree, variablesSpinCorr.c_qj_);
    //this->importBranch(tree, variablesSpinCorr.c_jq_);
    //this->importBranch(tree, variablesSpinCorr.c_nq_);
    //this->importBranch(tree, variablesSpinCorr.c_qn_);
    //this->importBranch(tree, variablesSpinCorr.c_nj_);
    //this->importBranch(tree, variablesSpinCorr.c_jn_);

    //this->importBranch(tree, variablesSpinCorr.c_Prk_);
    //this->importBranch(tree, variablesSpinCorr.c_Mrk_);
    //this->importBranch(tree, variablesSpinCorr.c_Pnr_);
    //this->importBranch(tree, variablesSpinCorr.c_Mnr_);
    //this->importBranch(tree, variablesSpinCorr.c_Pnk_);
    //this->importBranch(tree, variablesSpinCorr.c_Mnk_);

    //this->importBranch(tree, variablesSpinCorr.c_Pqj_);
    //this->importBranch(tree, variablesSpinCorr.c_Mqj_);
    //this->importBranch(tree, variablesSpinCorr.c_Pnq_);
    //this->importBranch(tree, variablesSpinCorr.c_Mnq_);
    //this->importBranch(tree, variablesSpinCorr.c_Pnj_);
    //this->importBranch(tree, variablesSpinCorr.c_Mnj_);

    //this->importBranch(tree, variablesSpinCorr.c_han_);
    //this->importBranch(tree, variablesSpinCorr.c_sca_);
    //this->importBranch(tree, variablesSpinCorr.c_tra_);
    //this->importBranch(tree, variablesSpinCorr.c_lij_);
    //this->importBranch(tree, variablesSpinCorr.c_liq_);

    //this->importBranch(tree, variablesSpinCorr.c_rkP_);
    //this->importBranch(tree, variablesSpinCorr.c_rkM_);
    //this->importBranch(tree, variablesSpinCorr.c_nrP_);
    //this->importBranch(tree, variablesSpinCorr.c_nrM_);
    //this->importBranch(tree, variablesSpinCorr.c_nkP_);
    //this->importBranch(tree, variablesSpinCorr.c_nkM_);

    //this->importBranch(tree, variablesSpinCorr.c_qjP_);
    //this->importBranch(tree, variablesSpinCorr.c_qjM_);
    //this->importBranch(tree, variablesSpinCorr.c_nqP_);
    //this->importBranch(tree, variablesSpinCorr.c_nqM_);
    //this->importBranch(tree, variablesSpinCorr.c_njP_);
    //this->importBranch(tree, variablesSpinCorr.c_njM_);

    this->importBranch(tree, variablesSpinCorr.ll_cHel_);
    this->importBranch(tree, variablesSpinCorr.ll_cLab_);
    this->importBranch(tree, variablesSpinCorr.ll_kNorm_);
    this->importBranch(tree, variablesSpinCorr.ll_rNorm_);

    this->importBranch(tree, variablesSpinCorr.gen_top_pt_);
    this->importBranch(tree, variablesSpinCorr.gen_top_phi_);
    this->importBranch(tree, variablesSpinCorr.gen_top_rapidity_);
    //    this->importBranch(tree, variablesSpinCorr.gen_top_arapidity_);
    this->importBranch(tree, variablesSpinCorr.gen_top_eta_);
    this->importBranch(tree, variablesSpinCorr.gen_top_mass_);

    this->importBranch(tree, variablesSpinCorr.gen_tbar_pt_);
    this->importBranch(tree, variablesSpinCorr.gen_tbar_phi_);
    this->importBranch(tree, variablesSpinCorr.gen_tbar_rapidity_);
    //    this->importBranch(tree, variablesSpinCorr.gen_tbar_arapidity_);
    this->importBranch(tree, variablesSpinCorr.gen_tbar_eta_);
    this->importBranch(tree, variablesSpinCorr.gen_tbar_mass_);

    this->importBranch(tree, variablesSpinCorr.gen_l_pt_);
    this->importBranch(tree, variablesSpinCorr.gen_l_eta_);
    this->importBranch(tree, variablesSpinCorr.gen_l_phi_);
    this->importBranch(tree, variablesSpinCorr.gen_l_mass_);
    this->importBranch(tree, variablesSpinCorr.gen_l_pdgid_);

    this->importBranch(tree, variablesSpinCorr.gen_lbar_pt_);
    this->importBranch(tree, variablesSpinCorr.gen_lbar_eta_);
    this->importBranch(tree, variablesSpinCorr.gen_lbar_phi_);
    this->importBranch(tree, variablesSpinCorr.gen_lbar_mass_);
    this->importBranch(tree, variablesSpinCorr.gen_lbar_pdgid_);

    //    this->importBranch(tree, variablesSpinCorr.gen_e_pt_);
    //    this->importBranch(tree, variablesSpinCorr.gen_e_eta_);
    //    this->importBranch(tree, variablesSpinCorr.gen_e_phi_);
    //    this->importBranch(tree, variablesSpinCorr.gen_e_mass_);

    //    this->importBranch(tree, variablesSpinCorr.gen_ebar_pt_);
    //    this->importBranch(tree, variablesSpinCorr.gen_ebar_eta_);
    //    this->importBranch(tree, variablesSpinCorr.gen_ebar_phi_);
    //    this->importBranch(tree, variablesSpinCorr.gen_ebar_mass_);

    //    this->importBranch(tree, variablesSpinCorr.gen_mu_pt_);
    //    this->importBranch(tree, variablesSpinCorr.gen_mu_eta_);
    //    this->importBranch(tree, variablesSpinCorr.gen_mu_phi_);
    //    this->importBranch(tree, variablesSpinCorr.gen_mu_mass_);

    //    this->importBranch(tree, variablesSpinCorr.gen_mubar_pt_);
    //    this->importBranch(tree, variablesSpinCorr.gen_mubar_eta_);
    //    this->importBranch(tree, variablesSpinCorr.gen_mubar_phi_);
    //    this->importBranch(tree, variablesSpinCorr.gen_mubar_mass_);

    this->importBranch(tree, variablesSpinCorr.gen_b_pt_);
    this->importBranch(tree, variablesSpinCorr.gen_b_eta_);
    this->importBranch(tree, variablesSpinCorr.gen_b_phi_);
    this->importBranch(tree, variablesSpinCorr.gen_b_mass_);

    this->importBranch(tree, variablesSpinCorr.gen_bbar_pt_);
    this->importBranch(tree, variablesSpinCorr.gen_bbar_eta_);
    this->importBranch(tree, variablesSpinCorr.gen_bbar_phi_);
    this->importBranch(tree, variablesSpinCorr.gen_bbar_mass_);

    this->importBranch(tree, variablesSpinCorr.gen_nu_pt_);
    this->importBranch(tree, variablesSpinCorr.gen_nu_eta_);
    this->importBranch(tree, variablesSpinCorr.gen_nu_phi_);
    this->importBranch(tree, variablesSpinCorr.gen_nu_mass_);

    this->importBranch(tree, variablesSpinCorr.gen_nubar_pt_);
    this->importBranch(tree, variablesSpinCorr.gen_nubar_eta_);
    this->importBranch(tree, variablesSpinCorr.gen_nubar_phi_);
    this->importBranch(tree, variablesSpinCorr.gen_nubar_mass_);

    this->importBranch(tree, variablesSpinCorr.gen_ttbar_pt_);
    this->importBranch(tree, variablesSpinCorr.gen_ttbar_phi_);
    this->importBranch(tree, variablesSpinCorr.gen_ttbar_rapidity_);
    //    this->importBranch(tree, variablesSpinCorr.gen_ttbar_arapidity_);
    this->importBranch(tree, variablesSpinCorr.gen_ttbar_eta_);
    this->importBranch(tree, variablesSpinCorr.gen_ttbar_delta_phi_);
    this->importBranch(tree, variablesSpinCorr.gen_ttbar_delta_eta_);
    this->importBranch(tree, variablesSpinCorr.gen_ttbar_delta_rapidity_);
    this->importBranch(tree, variablesSpinCorr.gen_ttbar_mass_);

    this->importBranch(tree, variablesSpinCorr.gen_llbar_pt_);
    this->importBranch(tree, variablesSpinCorr.gen_llbar_phi_);
    this->importBranch(tree, variablesSpinCorr.gen_llbar_rapidity_);
    //    this->importBranch(tree, variablesSpinCorr.gen_llbar_arapidity_);
    this->importBranch(tree, variablesSpinCorr.gen_llbar_delta_phi_);
    this->importBranch(tree, variablesSpinCorr.gen_llbar_delta_eta_);
    this->importBranch(tree, variablesSpinCorr.gen_llbar_delta_rapidity_);
    this->importBranch(tree, variablesSpinCorr.gen_llbar_mass_);

    //    this->importBranch(tree, variablesSpinCorr.gen_bbbar_pt_);
    //    this->importBranch(tree, variablesSpinCorr.gen_bbbar_phi_);
    //    this->importBranch(tree, variablesSpinCorr.gen_bbbar_rapidity_);
    //    this->importBranch(tree, variablesSpinCorr.gen_bbbar_arapidity_);
    //    this->importBranch(tree, variablesSpinCorr.gen_bbbar_delta_phi_);
    //    this->importBranch(tree, variablesSpinCorr.gen_bbbar_delta_eta_);
    //    this->importBranch(tree, variablesSpinCorr.gen_bbbar_delta_rapidity_);
    //    this->importBranch(tree, variablesSpinCorr.gen_bbbar_mass_);

    //    this->importBranch(tree, variablesSpinCorr.gen_nunubar_pt_);
    //    this->importBranch(tree, variablesSpinCorr.gen_nunubar_phi_);
    //    this->importBranch(tree, variablesSpinCorr.gen_nunubar_rapidity_);
    //    this->importBranch(tree, variablesSpinCorr.gen_nunubar_arapidity_);
    //    this->importBranch(tree, variablesSpinCorr.gen_nunubar_delta_phi_);
    //    this->importBranch(tree, variablesSpinCorr.gen_nunubar_delta_eta_);
    //    this->importBranch(tree, variablesSpinCorr.gen_nunubar_delta_rapidity_);
    //    this->importBranch(tree, variablesSpinCorr.gen_nunubar_mass_);

    this->importBranch(tree, variablesSpinCorr.gen_MT2_);

    this->importBranch(tree, variablesSpinCorr.gen_all_mass_);
    this->importBranch(tree, variablesSpinCorr.gen_r_mass_);
    this->importBranch(tree, variablesSpinCorr.gen_jet_multiplicity_);
    this->importBranch(tree, variablesSpinCorr.gen_x1_);
    this->importBranch(tree, variablesSpinCorr.gen_x2_);

    this->importBranch(tree, variablesSpinCorr.gen_production_mode_);

    this->importBranch(tree, variablesSpinCorr.gen_n_extraJets_iso08_);
    
    //ttbar frame variables
    this->importBranch(tree, variablesSpinCorr.gen_top_scatteringangle_ttbarframe_);

    //    this->importBranch(tree, variablesSpinCorr.gen_top_pt_ttbarframe_);
    //    this->importBranch(tree, variablesSpinCorr.gen_top_phi_ttbarframe_);
    //    this->importBranch(tree, variablesSpinCorr.gen_top_rapidity_ttbarframe_);
    //    this->importBranch(tree, variablesSpinCorr.gen_top_arapidity_ttbarframe_);
    //    this->importBranch(tree, variablesSpinCorr.gen_tbar_pt_ttbarframe_);
    //    this->importBranch(tree, variablesSpinCorr.gen_tbar_phi_ttbarframe_);
    //    this->importBranch(tree, variablesSpinCorr.gen_tbar_rapidity_ttbarframe_);
    //    this->importBranch(tree, variablesSpinCorr.gen_tbar_arapidity_ttbarframe_);

    //top and antitop frame variables
    //    this->importBranch(tree, variablesSpinCorr.gen_l_pt_antitopframe_);
    //    this->importBranch(tree, variablesSpinCorr.gen_l_eta_antitopframe_);
    //    this->importBranch(tree, variablesSpinCorr.gen_l_phi_antitopframe_);
    //    this->importBranch(tree, variablesSpinCorr.gen_lbar_pt_topframe_);
    //    this->importBranch(tree, variablesSpinCorr.gen_lbar_eta_topframe_);
    //    this->importBranch(tree, variablesSpinCorr.gen_lbar_phi_topframe_);

    //Spin correlation variables                                                                                                                                                          $
    this->importBranch(tree, variablesSpinCorr.gen_b1k_);
    this->importBranch(tree, variablesSpinCorr.gen_b2k_);
    this->importBranch(tree, variablesSpinCorr.gen_b1j_);
    this->importBranch(tree, variablesSpinCorr.gen_b2j_);
    this->importBranch(tree, variablesSpinCorr.gen_b1r_);
    this->importBranch(tree, variablesSpinCorr.gen_b2r_);
    this->importBranch(tree, variablesSpinCorr.gen_b1q_);
    this->importBranch(tree, variablesSpinCorr.gen_b2q_);
    this->importBranch(tree, variablesSpinCorr.gen_b1n_);
    this->importBranch(tree, variablesSpinCorr.gen_b2n_);

    //this->importBranch(tree, variablesSpinCorr.gen_b_Pkk_);
    //this->importBranch(tree, variablesSpinCorr.gen_b_Mkk_);
    //this->importBranch(tree, variablesSpinCorr.gen_b_Pjj_);
    //this->importBranch(tree, variablesSpinCorr.gen_b_Mjj_);
    //this->importBranch(tree, variablesSpinCorr.gen_b_Prr_);
    //this->importBranch(tree, variablesSpinCorr.gen_b_Mrr_);
    //this->importBranch(tree, variablesSpinCorr.gen_b_Pqq_);
    //this->importBranch(tree, variablesSpinCorr.gen_b_Mqq_);
    //this->importBranch(tree, variablesSpinCorr.gen_b_Pnn_);
    //this->importBranch(tree, variablesSpinCorr.gen_b_Mnn_);

    this->importBranch(tree, variablesSpinCorr.gen_c_kk_);
    this->importBranch(tree, variablesSpinCorr.gen_c_rr_);
    this->importBranch(tree, variablesSpinCorr.gen_c_nn_);
    this->importBranch(tree, variablesSpinCorr.gen_c_kj_);
    this->importBranch(tree, variablesSpinCorr.gen_c_rq_);

    this->importBranch(tree, variablesSpinCorr.gen_c_rk_);
    this->importBranch(tree, variablesSpinCorr.gen_c_kr_);
    this->importBranch(tree, variablesSpinCorr.gen_c_nr_);
    this->importBranch(tree, variablesSpinCorr.gen_c_rn_);
    this->importBranch(tree, variablesSpinCorr.gen_c_nk_);
    this->importBranch(tree, variablesSpinCorr.gen_c_kn_);

    this->importBranch(tree, variablesSpinCorr.gen_c_rj_);
    this->importBranch(tree, variablesSpinCorr.gen_c_jr_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_qk_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_kq_);

    //this->importBranch(tree, variablesSpinCorr.gen_c_qj_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_jq_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_nq_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_qn_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_nj_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_jn_);

    //this->importBranch(tree, variablesSpinCorr.gen_c_Prk_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_Mrk_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_Pnr_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_Mnr_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_Pnk_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_Mnk_);

    //this->importBranch(tree, variablesSpinCorr.gen_c_Pqj_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_Mqj_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_Pnq_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_Mnq_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_Pnj_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_Mnj_);

    //this->importBranch(tree, variablesSpinCorr.gen_c_han_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_sca_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_tra_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_lij_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_liq_);

    //this->importBranch(tree, variablesSpinCorr.gen_c_rkP_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_rkM_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_nrP_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_nrM_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_nkP_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_nkM_);

    //this->importBranch(tree, variablesSpinCorr.gen_c_qjP_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_qjM_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_nqP_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_nqM_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_njP_);
    //this->importBranch(tree, variablesSpinCorr.gen_c_njM_);

    this->importBranch(tree, variablesSpinCorr.gen_ll_cHel_);
    this->importBranch(tree, variablesSpinCorr.gen_ll_cLab_);
    this->importBranch(tree, variablesSpinCorr.gen_ll_kNorm_);
    this->importBranch(tree, variablesSpinCorr.gen_ll_rNorm_);

    this->importBranch(tree, variablesSpinCorr.entry_);
    this->importBranch(tree, variablesSpinCorr.isTopGen_);
    this->importBranch(tree, variablesSpinCorr.TopDecayMode_);
    this->importBranch(tree, variablesSpinCorr.isKinReco_);
    this->importBranch(tree, variablesSpinCorr.eventWeight_);
    this->importBranch(tree, variablesSpinCorr.trueLevelWeight_);



    // Loop over all tree entries and fill vector of structs
    for(Long64_t iEntry = 0; iEntry < tree->GetEntries(); ++iEntry){
        tree->GetEntry(iEntry);
        VariablesSpinCorr* const variablesSpinCorrPtr = new VariablesSpinCorr();
        *variablesSpinCorrPtr = variablesSpinCorr;
        v_variables.push_back(variablesSpinCorrPtr);
    }
}
