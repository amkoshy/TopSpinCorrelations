#ifndef VariablesSpinCorr_h
#define VariablesSpinCorr_h

#include <map>
#include <vector>
#include <string>

#include "VariablesBase.h"
#include "analysisStructsFwd.h"
#include "../../common/include/classesFwd.h"
#include "TLorentzVector.h"

class EventMetadata;
class RecoObjects;
class TopGenObjects;
class CommonGenObjects;
class KinematicReconstructionSolutions;
namespace ttbar{
    class RecoObjectIndices;
    class GenObjectIndices;
    class GenLevelWeights;
    class RecoLevelWeights;
}









/// Class holding the input variables for one entry for MVA for top system jet identification,
/// i.e. one entry of each quantity per selected jet combination
class VariablesSpinCorr : public VariablesBase{
    
public:
    
    /// Empty constructor
    VariablesSpinCorr();
    
    /// Constructor setting up input variables from physics objects
    VariablesSpinCorr(const EventMetadata& eventMetadata, 
                   const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                                          const TopGenObjects& topGenObjects,
                                                          const KinematicReconstructionSolutions& kinRecoObjects,
                                                          const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                                                          const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                                                          const double& weight);
    
    /// Destructor
    ~VariablesSpinCorr(){}
    
    /// Fill the input structs for all jet combinations of one event
    static std::vector<VariablesBase*> fillVariables(const EventMetadata& eventMetadata,
                                                     const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                                          const TopGenObjects& topGenObjects,
                                                          const KinematicReconstructionSolutions& kinRecoObjects,
                                                          const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                                                          const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                                                          const double& weight);
    
    double getMT2Variable(TLorentzVector lep, TLorentzVector antilep, TLorentzVector MET);

    VariableFloat top_pt_;
    VariableFloat top_phi_;
    VariableFloat top_rapidity_;
    //    VariableFloat top_arapidity_;
    VariableFloat top_eta_;
    VariableFloat top_mass_;

    VariableFloat tbar_pt_;
    VariableFloat tbar_phi_;
    VariableFloat tbar_rapidity_;
    //    VariableFloat tbar_arapidity_;
    VariableFloat tbar_eta_;
    VariableFloat tbar_mass_;

    VariableFloat l_pt_;
    VariableFloat l_eta_;
    VariableFloat l_phi_;
    VariableFloat l_mass_;
    VariableFloat l_pdgid_;

    VariableFloat lbar_pt_;
    VariableFloat lbar_eta_;
    VariableFloat lbar_phi_;
    VariableFloat lbar_mass_;
    VariableFloat lbar_pdgid_;

    //    VariableFloat e_pt_;
    //    VariableFloat e_eta_;
    //    VariableFloat e_phi_;
    //    VariableFloat e_mass_;

    //    VariableFloat ebar_pt_;
    //    VariableFloat ebar_eta_;
    //    VariableFloat ebar_phi_;
    //    VariableFloat ebar_mass_;

    //    VariableFloat mu_pt_;
    //    VariableFloat mu_eta_;
    //    VariableFloat mu_phi_;
    //    VariableFloat mu_mass_;

    //    VariableFloat mubar_pt_;
    //    VariableFloat mubar_eta_;
    //    VariableFloat mubar_phi_;
    //    VariableFloat mubar_mass_;

    VariableFloat b_pt_;
    VariableFloat b_eta_;
    VariableFloat b_phi_;
    VariableFloat b_mass_;

    VariableFloat bbar_pt_;
    VariableFloat bbar_eta_;
    VariableFloat bbar_phi_;
    VariableFloat bbar_mass_;

    VariableFloat nu_pt_;
    VariableFloat nu_eta_;
    VariableFloat nu_phi_;
    VariableFloat nu_mass_;

    VariableFloat nubar_pt_;
    VariableFloat nubar_eta_;
    VariableFloat nubar_phi_;
    VariableFloat nubar_mass_;

    VariableFloat met_pt_;
    VariableFloat met_phi_;
    VariableFloat met_mass_;

    VariableFloat ttbar_pt_;
    VariableFloat ttbar_phi_;
    VariableFloat ttbar_rapidity_;
    //    VariableFloat ttbar_arapidity_;
    VariableFloat ttbar_eta_;
    VariableFloat ttbar_delta_phi_;
    VariableFloat ttbar_delta_eta_;
    VariableFloat ttbar_delta_rapidity_;
    VariableFloat ttbar_mass_;

    VariableFloat llbar_pt_;
    VariableFloat llbar_phi_;
    VariableFloat llbar_rapidity_;
    //    VariableFloat llbar_arapidity_;
    VariableFloat llbar_delta_phi_;
    VariableFloat llbar_delta_eta_;
    VariableFloat llbar_delta_rapidity_;
    VariableFloat llbar_mass_;

    //    VariableFloat bbbar_pt_;
    //    VariableFloat bbbar_phi_;
    //    VariableFloat bbbar_rapidity_; 
    //    VariableFloat bbbar_arapidity_;
    //    VariableFloat bbbar_delta_phi_;
    //    VariableFloat bbbar_delta_eta_;
    //    VariableFloat bbbar_delta_rapidity_;
    //    VariableFloat bbbar_mass_;

    //    VariableFloat nunubar_pt_;
    //    VariableFloat nunubar_phi_;
    //    VariableFloat nunubar_rapidity_; 
    //    VariableFloat nunubar_arapidity_;
    //    VariableFloat nunubar_delta_phi_;
    //    VariableFloat nunubar_delta_eta_;
    //    VariableFloat nunubar_delta_rapidity_;
    //    VariableFloat nunubar_mass_;

    VariableFloat MT2_;

    VariableFloat all_mass_;
    VariableFloat r_mass_;
    VariableFloat jet_multiplicity_;
    VariableFloat bjet_multiplicity_;
    VariableFloat x1_;
    VariableFloat x2_;

    //extraJet variables
    VariableFloat n_extraJets_iso08_;
    //ttbar frame variables
    VariableFloat top_scatteringangle_ttbarframe_;

    //    VariableFloat top_pt_ttbarframe_;
    //    VariableFloat top_phi_ttbarframe_;
    //    VariableFloat top_rapidity_ttbarframe_;
    //    VariableFloat top_arapidity_ttbarframe_;
    //    VariableFloat tbar_pt_ttbarframe_;
    //    VariableFloat tbar_phi_ttbarframe_;
    //    VariableFloat tbar_rapidity_ttbarframe_;
    //    VariableFloat tbar_arapidity_ttbarframe_;

    //top and antitop frame variables
    //    VariableFloat l_pt_antitopframe_;
    //    VariableFloat l_eta_antitopframe_;
    //    VariableFloat l_phi_antitopframe_;
    //    VariableFloat lbar_pt_topframe_;
    //    VariableFloat lbar_eta_topframe_;
    //    VariableFloat lbar_phi_topframe_;

    //Spin correlation variables
    VariableFloat b1k_;
    VariableFloat b2k_;
    VariableFloat b1j_;
    VariableFloat b2j_;
    VariableFloat b1r_;
    VariableFloat b2r_;
    VariableFloat b1q_;
    VariableFloat b2q_;
    VariableFloat b1n_;
    VariableFloat b2n_;

    //VariableFloat b_Pkk_;
    //VariableFloat b_Mkk_;
    //VariableFloat b_Pjj_;
    //VariableFloat b_Mjj_;
    //VariableFloat b_Prr_;
    //VariableFloat b_Mrr_;
    //VariableFloat b_Pqq_;
    //VariableFloat b_Mqq_;
    //VariableFloat b_Pnn_;
    //VariableFloat b_Mnn_;

    VariableFloat c_kk_;
    VariableFloat c_rr_;
    VariableFloat c_nn_;
    VariableFloat c_kj_;
    VariableFloat c_rq_;

    VariableFloat c_rk_;
    VariableFloat c_kr_;
    VariableFloat c_nr_;
    VariableFloat c_rn_;
    VariableFloat c_nk_;
    VariableFloat c_kn_;

    VariableFloat c_rj_;
    VariableFloat c_jr_;
    //VariableFloat c_qk_;
    //VariableFloat c_kq_;

    //VariableFloat c_qj_;
    //VariableFloat c_jq_;
    //VariableFloat c_nq_;
    //VariableFloat c_qn_;
    //VariableFloat c_nj_;
    //VariableFloat c_jn_;

    //VariableFloat c_Prk_;
    //VariableFloat c_Mrk_;
    //VariableFloat c_Pnr_;
    //VariableFloat c_Mnr_;
    //VariableFloat c_Pnk_;
    //VariableFloat c_Mnk_;

    //VariableFloat c_Pqj_;
    //VariableFloat c_Mqj_;
    //VariableFloat c_Pnq_;
    //VariableFloat c_Mnq_;
    //VariableFloat c_Pnj_;
    //VariableFloat c_Mnj_;

    //VariableFloat c_han_;
    //VariableFloat c_sca_;
    //VariableFloat c_tra_;
    //VariableFloat c_lij_;
    //VariableFloat c_liq_;

    //VariableFloat c_rkP_;
    //VariableFloat c_rkM_;
    //VariableFloat c_nrP_;
    //VariableFloat c_nrM_;
    //VariableFloat c_nkP_;
    //VariableFloat c_nkM_;

    //VariableFloat c_qjP_;
    //VariableFloat c_qjM_;
    //VariableFloat c_nqP_;
    //VariableFloat c_nqM_;
    //VariableFloat c_njP_;
    //VariableFloat c_njM_;

    VariableFloat ll_cHel_;
    VariableFloat ll_cLab_;
    VariableFloat ll_kNorm_;
    VariableFloat ll_rNorm_;

    VariableFloat gen_top_pt_;
    VariableFloat gen_top_phi_;
    VariableFloat gen_top_rapidity_;
    //    VariableFloat gen_top_arapidity_;
    VariableFloat gen_top_eta_;
    VariableFloat gen_top_mass_;

    VariableFloat gen_tbar_pt_;
    VariableFloat gen_tbar_phi_;
    VariableFloat gen_tbar_rapidity_;
    //    VariableFloat gen_tbar_arapidity_;
    VariableFloat gen_tbar_eta_;
    VariableFloat gen_tbar_mass_;

    VariableFloat gen_l_pt_;
    VariableFloat gen_lbar_pt_;
    //    VariableFloat gen_e_pt_;
    //    VariableFloat gen_ebar_pt_;
    //    VariableFloat gen_mu_pt_;
    //    VariableFloat gen_mubar_pt_;
    VariableFloat gen_b_pt_;
    VariableFloat gen_bbar_pt_;
    VariableFloat gen_nu_pt_;
    VariableFloat gen_nubar_pt_;
    VariableFloat gen_l_eta_;
    VariableFloat gen_lbar_eta_;
    //    VariableFloat gen_e_eta_;
    //    VariableFloat gen_ebar_eta_;
    //    VariableFloat gen_mu_eta_;
    //    VariableFloat gen_mubar_eta_;
    VariableFloat gen_b_eta_;
    VariableFloat gen_bbar_eta_;
    VariableFloat gen_nu_eta_;
    VariableFloat gen_nubar_eta_;
    VariableFloat gen_l_phi_;
    VariableFloat gen_lbar_phi_;
    //    VariableFloat gen_e_phi_;
    //    VariableFloat gen_ebar_phi_;
    //    VariableFloat gen_mu_phi_;
    //    VariableFloat gen_mubar_phi_;
    VariableFloat gen_b_phi_;
    VariableFloat gen_bbar_phi_;
    VariableFloat gen_nu_phi_;
    VariableFloat gen_nubar_phi_;
    VariableFloat gen_l_mass_;
    VariableFloat gen_lbar_mass_;
    VariableFloat gen_l_pdgid_;
    VariableFloat gen_lbar_pdgid_;
    //    VariableFloat gen_e_mass_;
    //    VariableFloat gen_ebar_mass_;
    //    VariableFloat gen_mu_mass_;
    //    VariableFloat gen_mubar_mass_;
    VariableFloat gen_b_mass_;
    VariableFloat gen_bbar_mass_;
    VariableFloat gen_nu_mass_;
    VariableFloat gen_nubar_mass_;

    VariableFloat gen_ttbar_pt_;
    VariableFloat gen_ttbar_phi_;
    VariableFloat gen_ttbar_rapidity_;
    //    VariableFloat gen_ttbar_arapidity_;
    VariableFloat gen_ttbar_eta_;
    VariableFloat gen_ttbar_delta_phi_;
    VariableFloat gen_ttbar_delta_eta_;
    VariableFloat gen_ttbar_delta_rapidity_;
    VariableFloat gen_ttbar_mass_;

    VariableFloat gen_llbar_pt_;
    VariableFloat gen_llbar_phi_;
    VariableFloat gen_llbar_rapidity_;
    //    VariableFloat gen_llbar_arapidity_;
    VariableFloat gen_llbar_delta_phi_;
    VariableFloat gen_llbar_delta_eta_;
    VariableFloat gen_llbar_delta_rapidity_;
    VariableFloat gen_llbar_mass_;

    //    VariableFloat gen_bbbar_pt_;
    //    VariableFloat gen_bbbar_phi_;
    //    VariableFloat gen_bbbar_rapidity_;
    //    VariableFloat gen_bbbar_arapidity_;
    //    VariableFloat gen_bbbar_delta_phi_;
    //    VariableFloat gen_bbbar_delta_eta_;
    //    VariableFloat gen_bbbar_delta_rapidity_;
    //    VariableFloat gen_bbbar_mass_;

    //    VariableFloat gen_nunubar_pt_;
    //    VariableFloat gen_nunubar_phi_;
    //    VariableFloat gen_nunubar_rapidity_; 
    //    VariableFloat gen_nunubar_arapidity_;
    //    VariableFloat gen_nunubar_delta_phi_;
    //    VariableFloat gen_nunubar_delta_eta_;
    //    VariableFloat gen_nunubar_delta_rapidity_;
    //    VariableFloat gen_nunubar_mass_;

    VariableFloat gen_MT2_;

    VariableFloat gen_all_mass_;
    VariableFloat gen_r_mass_;
    VariableFloat gen_jet_multiplicity_;
    VariableFloat gen_x1_;
    VariableFloat gen_x2_;
    VariableFloat gen_production_mode_;

    //extraJet variables
    VariableFloat gen_n_extraJets_iso08_;

    //ttbar frame variables
    VariableFloat gen_top_scatteringangle_ttbarframe_;

    //    VariableFloat gen_top_pt_ttbarframe_;
    //    VariableFloat gen_top_phi_ttbarframe_;
    //    VariableFloat gen_top_rapidity_ttbarframe_;
    //    VariableFloat gen_top_arapidity_ttbarframe_;
    //    VariableFloat gen_tbar_pt_ttbarframe_;
    //    VariableFloat gen_tbar_phi_ttbarframe_;
    //    VariableFloat gen_tbar_rapidity_ttbarframe_;
    //    VariableFloat gen_tbar_arapidity_ttbarframe_;

    //top and antitop frame variables
    //    VariableFloat gen_l_pt_antitopframe_;
    //    VariableFloat gen_l_eta_antitopframe_;
    //    VariableFloat gen_l_phi_antitopframe_;
    //    VariableFloat gen_lbar_pt_topframe_;
    //    VariableFloat gen_lbar_eta_topframe_;
    //    VariableFloat gen_lbar_phi_topframe_;

    //Spin correlation variables
    VariableFloat gen_b1k_;
    VariableFloat gen_b2k_;
    VariableFloat gen_b1j_;
    VariableFloat gen_b2j_;
    VariableFloat gen_b1r_;
    VariableFloat gen_b2r_;
    VariableFloat gen_b1q_;
    VariableFloat gen_b2q_;
    VariableFloat gen_b1n_;
    VariableFloat gen_b2n_;

    //VariableFloat gen_b_Pkk_;
    //VariableFloat gen_b_Mkk_;
    //VariableFloat gen_b_Pjj_;
    //VariableFloat gen_b_Mjj_;
    //VariableFloat gen_b_Prr_;
    //VariableFloat gen_b_Mrr_;
    //VariableFloat gen_b_Pqq_;
    //VariableFloat gen_b_Mqq_;
    //VariableFloat gen_b_Pnn_;
    //VariableFloat gen_b_Mnn_;

    VariableFloat gen_c_kk_;
    VariableFloat gen_c_rr_;
    VariableFloat gen_c_nn_;
    VariableFloat gen_c_kj_;
    VariableFloat gen_c_rq_;

    VariableFloat gen_c_rk_;
    VariableFloat gen_c_kr_;
    VariableFloat gen_c_nr_;
    VariableFloat gen_c_rn_;
    VariableFloat gen_c_nk_;
    VariableFloat gen_c_kn_;

    VariableFloat gen_c_rj_;
    VariableFloat gen_c_jr_;
    //VariableFloat gen_c_qk_;
    //VariableFloat gen_c_kq_;

    //VariableFloat gen_c_qj_;
    //VariableFloat gen_c_jq_;
    //VariableFloat gen_c_nq_;
    //VariableFloat gen_c_qn_;
    //VariableFloat gen_c_nj_;
    //VariableFloat gen_c_jn_;

    //VariableFloat gen_c_Prk_;
    //VariableFloat gen_c_Mrk_;
    //VariableFloat gen_c_Pnr_;
    //VariableFloat gen_c_Mnr_;
    //VariableFloat gen_c_Pnk_;
    //VariableFloat gen_c_Mnk_;

    //VariableFloat gen_c_Pqj_;
    //VariableFloat gen_c_Mqj_;
    //VariableFloat gen_c_Pnq_;
    //VariableFloat gen_c_Mnq_;
    //VariableFloat gen_c_Pnj_;
    //VariableFloat gen_c_Mnj_;

    //VariableFloat gen_c_han_;
    //VariableFloat gen_c_sca_;
    //VariableFloat gen_c_tra_;
    //VariableFloat gen_c_lij_;
    //VariableFloat gen_c_liq_;

    //VariableFloat gen_c_rkP_;
    //VariableFloat gen_c_rkM_;
    //VariableFloat gen_c_nrP_;
    //VariableFloat gen_c_nrM_;
    //VariableFloat gen_c_nkP_;
    //VariableFloat gen_c_nkM_;

    //VariableFloat gen_c_qjP_;
    //VariableFloat gen_c_qjM_;
    //VariableFloat gen_c_nqP_;
    //VariableFloat gen_c_nqM_;
    //VariableFloat gen_c_njP_;
    //VariableFloat gen_c_njM_;

    VariableFloat gen_ll_cHel_;
    VariableFloat gen_ll_cLab_;
    VariableFloat gen_ll_kNorm_;
    VariableFloat gen_ll_rNorm_;

    VariableInt entry_;
    VariableInt isTopGen_;
    VariableInt isKinReco_;
    VariableInt TopDecayMode_;
    VariableFloat eventWeight_;
    VariableFloat trueLevelWeight_;   
    
private:
    
    // The names associated to the variables
    
    static constexpr const char* name_top_pt_ = "top_pt";
    static constexpr const char* name_top_phi_ = "top_phi";
    static constexpr const char* name_top_rapidity_ = "top_rapidity";
    //    static constexpr const char* name_top_arapidity_ = "top_arapidity";
    static constexpr const char* name_top_eta_ = "top_eta";
    static constexpr const char* name_top_mass_ = "top_mass";


    static constexpr const char* name_tbar_pt_ = "tbar_pt";
    static constexpr const char* name_tbar_phi_ = "tbar_phi";
    static constexpr const char* name_tbar_rapidity_ = "tbar_rapidity";
    //    static constexpr const char* name_tbar_arapidity_ = "tbar_arapidity";
    static constexpr const char* name_tbar_eta_ = "tbar_eta";
    static constexpr const char* name_tbar_mass_ = "tbar_mass";


    static constexpr const char* name_l_pt_ = "l_pt";
    static constexpr const char* name_l_eta_ = "l_eta";
    static constexpr const char* name_l_phi_ = "l_phi";
    static constexpr const char* name_l_mass_ = "l_mass";
    static constexpr const char* name_l_pdgid_ = "l_pdgid";

    static constexpr const char* name_lbar_pt_ = "lbar_pt";
    static constexpr const char* name_lbar_eta_ = "lbar_eta";
    static constexpr const char* name_lbar_phi_ = "lbar_phi";
    static constexpr const char* name_lbar_mass_ = "lbar_mass";
    static constexpr const char* name_lbar_pdgid_ = "lbar_pdgid";


    //    static constexpr const char* name_e_pt_ = "e_pt";
    //    static constexpr const char* name_e_eta_ = "e_eta";
    //    static constexpr const char* name_e_phi_ = "e_phi";
    //    static constexpr const char* name_e_mass_ = "e_mass";

    //    static constexpr const char* name_ebar_pt_ = "ebar_pt";
    //    static constexpr const char* name_ebar_eta_ = "ebar_eta";
    //    static constexpr const char* name_ebar_phi_ = "ebar_phi";

    //    static constexpr const char* name_mu_pt_ = "mu_pt";
    //    static constexpr const char* name_mu_eta_ = "mu_eta";
    //    static constexpr const char* name_mu_phi_ = "mu_phi";
    //    static constexpr const char* name_mu_mass_ = "mu_mass";

    //    static constexpr const char* name_mubar_pt_ = "mubar_pt";
    //    static constexpr const char* name_mubar_eta_ = "mubar_eta";
    //    static constexpr const char* name_mubar_phi_ = "mubar_phi";
    //    static constexpr const char* name_mubar_mass_ = "mubar_mass";

    static constexpr const char* name_b_pt_ = "b_pt";
    static constexpr const char* name_b_eta_ = "b_eta";
    static constexpr const char* name_b_phi_ = "b_phi";
    static constexpr const char* name_b_mass_ = "b_mass";

    static constexpr const char* name_bbar_pt_ = "bbar_pt";
    static constexpr const char* name_bbar_eta_ = "bbar_eta";
    static constexpr const char* name_bbar_phi_ = "bbar_phi";
    static constexpr const char* name_bbar_mass_ = "bbar_mass";

    static constexpr const char* name_nu_pt_ = "nu_pt";
    static constexpr const char* name_nu_eta_ = "nu_eta";
    static constexpr const char* name_nu_phi_ = "nu_phi";    
    static constexpr const char* name_nu_mass_ = "nu_mass";

    static constexpr const char* name_nubar_pt_ = "nubar_pt";
    static constexpr const char* name_nubar_eta_ = "nubar_eta";
    static constexpr const char* name_nubar_phi_ = "nubar_phi";
    static constexpr const char* name_nubar_mass_ = "nubar_mass";

    static constexpr const char* name_met_pt_ = "met_pt";
    static constexpr const char* name_met_phi_ = "met_phi";
    static constexpr const char* name_met_mass_ = "met_mass";

    static constexpr const char* name_ttbar_pt_ = "ttbar_pt";
    static constexpr const char* name_ttbar_phi_ = "ttbar_phi";
    static constexpr const char* name_ttbar_rapidity_ = "ttbar_rapidity";
    //    static constexpr const char* name_ttbar_arapidity_ = "ttbar_arapidity";
    static constexpr const char* name_ttbar_delta_phi_ = "ttbar_delta_phi";
    static constexpr const char* name_ttbar_delta_eta_ = "ttbar_delta_eta";
    static constexpr const char* name_ttbar_delta_rapidity_ = "ttbar_delta_rapidity";
    static constexpr const char* name_ttbar_mass_ = "ttbar_mass";
    static constexpr const char* name_ttbar_eta_ = "ttbar_eta";

    static constexpr const char* name_llbar_pt_ = "llbar_pt";
    static constexpr const char* name_llbar_phi_ = "llbar_phi";
    static constexpr const char* name_llbar_rapidity_ = "llbar_rapidity";
    //    static constexpr const char* name_llbar_arapidity_ = "llbar_arapidity";
    static constexpr const char* name_llbar_delta_phi_ = "llbar_delta_phi";
    static constexpr const char* name_llbar_delta_eta_ = "llbar_delta_eta";
    static constexpr const char* name_llbar_delta_rapidity_ = "llbar_delta_rapidity";
    static constexpr const char* name_llbar_mass_ = "llbar_mass";

    //    static constexpr const char* name_bbbar_pt_ = "bbbar_pt";
    //    static constexpr const char* name_bbbar_phi_ = "bbbar_phi";
    //    static constexpr const char* name_bbbar_rapidity_ = "bbbar_rapidity";
    //    static constexpr const char* name_bbbar_arapidity_ = "bbbar_arapidity";
    //    static constexpr const char* name_bbbar_delta_phi_ = "bbbar_delta_phi";
    //    static constexpr const char* name_bbbar_delta_eta_ = "bbbar_delta_eta";
    //    static constexpr const char* name_bbbar_delta_rapidity_ = "bbbar_delta_rapidity";
    //    static constexpr const char* name_bbbar_mass_ = "bbbar_mass";

    //    static constexpr const char* name_nunubar_pt_ = "nunubar_pt";
    //    static constexpr const char* name_nunubar_phi_ = "nunubar_phi";
    //    static constexpr const char* name_nunubar_rapidity_ = "nunubar_rapidity";
    //    static constexpr const char* name_nunubar_arapidity_ = "nunubar_arapidity";
    //    static constexpr const char* name_nunubar_delta_phi_ = "nunubar_delta_phi";
    //    static constexpr const char* name_nunubar_delta_eta_ = "nunubar_delta_eta";
    //    static constexpr const char* name_nunubar_delta_rapidity_ = "nunubar_delta_rapidity";
    //    static constexpr const char* name_nunubar_mass_ = "nunubar_mass";

    static constexpr const char* name_MT2_ = "MT2";

    static constexpr const char* name_all_mass_ = "all_mass";
    static constexpr const char* name_r_mass_ = "r_mass";
    static constexpr const char* name_jet_multiplicity_ = "jet_multiplicity";
    static constexpr const char* name_bjet_multiplicity_ = "bjet_multiplicity";
    static constexpr const char* name_x1_ = "x1";
    static constexpr const char* name_x2_ = "x2";

    //extraJet variables
    static constexpr const char* name_n_extraJets_iso08_ = "n_extraJets_iso08";

    //ttbar frame variables
    static constexpr const char* name_top_scatteringangle_ttbarframe_ = "top_scatteringangle_ttbarframe";

    //    static constexpr const char* name_top_pt_ttbarframe_ = "top_pt_ttbarframe";
    //    static constexpr const char* name_top_phi_ttbarframe_ = "top_phi_ttbarframe";
    //    static constexpr const char* name_top_rapidity_ttbarframe_ = "top_rapidity_ttbarframe";
    //    static constexpr const char* name_top_arapidity_ttbarframe_ = "top_arapidity_ttbarframe";
    //    static constexpr const char* name_tbar_pt_ttbarframe_ = "tbar_pt_ttbarframe";
    //    static constexpr const char* name_tbar_phi_ttbarframe_ = "tbar_phi_ttbarframe";
    //    static constexpr const char* name_tbar_rapidity_ttbarframe_ = "tbar_rapidity_ttbarframe";
    //    static constexpr const char* name_tbar_arapidity_ttbarframe_ = "tbar_arapidity_ttbarframe";

    //top and antitop frame variables
    //    static constexpr const char* name_l_pt_antitopframe_ = "l_pt_antitopframe";
    //    static constexpr const char* name_l_eta_antitopframe_ = "l_eta_antitopframe";
    //    static constexpr const char* name_l_phi_antitopframe_ = "l_phi_antitopframe";
    //    static constexpr const char* name_lbar_pt_topframe_ = "lbar_pt_topframe";
    //    static constexpr const char* name_lbar_eta_topframe_ = "lbar_eta_topframe";
    //    static constexpr const char* name_lbar_phi_topframe_ = "lbar_phi_topframe";


    //Spin correlation variables
    static constexpr const char* name_b1k_ = "b1k";
    static constexpr const char* name_b2k_ = "b2k";
    static constexpr const char* name_b1j_ = "b1j";
    static constexpr const char* name_b2j_ = "b2j";
    static constexpr const char* name_b1r_ = "b1r";
    static constexpr const char* name_b2r_ = "b2r";
    static constexpr const char* name_b1q_ = "b1q";
    static constexpr const char* name_b2q_ = "b2q";
    static constexpr const char* name_b1n_ = "b1n";
    static constexpr const char* name_b2n_ = "b2n";

    //static constexpr const char* name_b_Pkk_ = "b_Pkk";
    //static constexpr const char* name_b_Mkk_ = "b_Mkk";
    //static constexpr const char* name_b_Pjj_ = "b_Pjj";
    //static constexpr const char* name_b_Mjj_ = "b_Mjj";
    //static constexpr const char* name_b_Prr_ = "b_Prr";
    //static constexpr const char* name_b_Mrr_ = "b_Mrr";
    //static constexpr const char* name_b_Pqq_ = "b_Pqq";
    //static constexpr const char* name_b_Mqq_ = "b_Mqq";
    //static constexpr const char* name_b_Pnn_ = "b_Pnn";
    //static constexpr const char* name_b_Mnn_ = "b_Mnn";

    static constexpr const char* name_c_kk_ = "c_kk";
    static constexpr const char* name_c_rr_ = "c_rr";
    static constexpr const char* name_c_nn_ = "c_nn";
    static constexpr const char* name_c_kj_ = "c_kj";
    static constexpr const char* name_c_rq_ = "c_rq";

    static constexpr const char* name_c_rk_ = "c_rk";
    static constexpr const char* name_c_kr_ = "c_kr";
    static constexpr const char* name_c_nr_ = "c_nr";
    static constexpr const char* name_c_rn_ = "c_rn";
    static constexpr const char* name_c_nk_ = "c_nk";
    static constexpr const char* name_c_kn_ = "c_kn";

    static constexpr const char* name_c_rj_ = "c_rj";
    static constexpr const char* name_c_jr_ = "c_jr";
    //static constexpr const char* name_c_qk_ = "c_qk";
    //static constexpr const char* name_c_kq_ = "c_kq";

    //static constexpr const char* name_c_qj_ = "c_qj";
    //static constexpr const char* name_c_jq_ = "c_jq";
    //static constexpr const char* name_c_nq_ = "c_nq";
    //static constexpr const char* name_c_qn_ = "c_qn";
    //static constexpr const char* name_c_nj_ = "c_nj";
    //static constexpr const char* name_c_jn_ = "c_jn";

    //static constexpr const char* name_c_Prk_ = "c_Prk";
    //static constexpr const char* name_c_Mrk_ = "c_Mrk";
    //static constexpr const char* name_c_Pnr_ = "c_Pnr";
    //static constexpr const char* name_c_Mnr_ = "c_Mnr";
    //static constexpr const char* name_c_Pnk_ = "c_Pnk";
    //static constexpr const char* name_c_Mnk_ = "c_Mnk";

    //static constexpr const char* name_c_Pqj_ = "c_Pqj";
    //static constexpr const char* name_c_Mqj_ = "c_Mqj";
    //static constexpr const char* name_c_Pnq_ = "c_Pnq";
    //static constexpr const char* name_c_Mnq_ = "c_Mnq";
    //static constexpr const char* name_c_Pnj_ = "c_Pnj";
    //static constexpr const char* name_c_Mnj_ = "c_Mnj";

    //static constexpr const char* name_c_han_ = "c_han";
    //static constexpr const char* name_c_sca_ = "c_sca";
    //static constexpr const char* name_c_tra_ = "c_tra";
    //static constexpr const char* name_c_lij_ = "c_lij";
    //static constexpr const char* name_c_liq_ = "c_liq";

    //static constexpr const char* name_c_rkP_ = "c_rkP";
    //static constexpr const char* name_c_rkM_ = "c_rkM";
    //static constexpr const char* name_c_nrP_ = "c_nrP";
    //static constexpr const char* name_c_nrM_ = "c_nrM";
    //static constexpr const char* name_c_nkP_ = "c_nkP";
    //static constexpr const char* name_c_nkM_ = "c_nkM";

    //static constexpr const char* name_c_qjP_ = "c_qjP";
    //static constexpr const char* name_c_qjM_ = "c_qjM";
    //static constexpr const char* name_c_nqP_ = "c_nqP";
    //static constexpr const char* name_c_nqM_ = "c_nqM";
    //static constexpr const char* name_c_njP_ = "c_njP";
    //static constexpr const char* name_c_njM_ = "c_njM";

    static constexpr const char* name_ll_cHel_ = "ll_cHel";
    static constexpr const char* name_ll_cLab_ = "ll_cLab";
    static constexpr const char* name_ll_kNorm_ = "ll_kNorm";
    static constexpr const char* name_ll_rNorm_ = "ll_rNorm";

    static constexpr const char* name_gen_top_pt_ = "gen_top_pt";
    static constexpr const char* name_gen_top_phi_ = "gen_top_phi";
    static constexpr const char* name_gen_top_rapidity_ = "gen_top_rapidity";
    //    static constexpr const char* name_gen_top_arapidity_ = "gen_top_arapidity";
    static constexpr const char* name_gen_top_eta_ = "gen_top_eta";
    static constexpr const char* name_gen_top_mass_ = "gen_top_mass";

    static constexpr const char* name_gen_tbar_pt_ = "gen_tbar_pt";
    static constexpr const char* name_gen_tbar_phi_ = "gen_tbar_phi";
    static constexpr const char* name_gen_tbar_rapidity_ = "gen_tbar_rapidity";
    //    static constexpr const char* name_gen_tbar_arapidity_ = "gen_tbar_arapidity";
    static constexpr const char* name_gen_tbar_eta_ = "gen_tbar_eta";
    static constexpr const char* name_gen_tbar_mass_ = "gen_tbar_mass";

    static constexpr const char* name_gen_l_pt_ = "gen_l_pt";
    static constexpr const char* name_gen_l_eta_ = "gen_l_eta";
    static constexpr const char* name_gen_l_phi_ = "gen_l_phi";
    static constexpr const char* name_gen_l_mass_ = "gen_l_mass";
    static constexpr const char* name_gen_l_pdgid_ = "gen_l_pdgid";

    static constexpr const char* name_gen_lbar_pt_ = "gen_lbar_pt";
    static constexpr const char* name_gen_lbar_eta_ = "gen_lbar_eta";
    static constexpr const char* name_gen_lbar_phi_ = "gen_lbar_phi";
    static constexpr const char* name_gen_lbar_mass_ = "gen_lbar_mass";
    static constexpr const char* name_gen_lbar_pdgid_ = "gen_lbar_pdgid";

    //    static constexpr const char* name_gen_e_pt_ = "gen_e_pt";
    //    static constexpr const char* name_gen_e_eta_ = "gen_e_eta";
    //    static constexpr const char* name_gen_e_phi_ = "gen_e_phi";
    //    static constexpr const char* name_gen_e_mass_ = "gen_e_mass";

    //    static constexpr const char* name_gen_ebar_pt_ = "gen_ebar_pt";
    //    static constexpr const char* name_gen_ebar_eta_ = "gen_ebar_eta";
    //    static constexpr const char* name_gen_ebar_phi_ = "gen_ebar_phi";
    //    static constexpr const char* name_gen_ebar_mass_ = "gen_ebar_mass";

    //    static constexpr const char* name_gen_mu_pt_ = "gen_mu_pt";
    //    static constexpr const char* name_gen_mu_eta_ = "gen_mu_eta";
    //    static constexpr const char* name_gen_mu_phi_ = "gen_mu_phi";
    //    static constexpr const char* name_gen_mu_mass_ = "gen_mu_mass";

    //    static constexpr const char* name_gen_mubar_pt_ = "gen_mubar_pt";
    //    static constexpr const char* name_gen_mubar_eta_ = "gen_mubar_eta";
    //    static constexpr const char* name_gen_mubar_phi_ = "gen_mubar_phi";
    //    static constexpr const char* name_gen_mubar_mass_ = "gen_mubar_mass";

    static constexpr const char* name_gen_b_pt_ = "gen_b_pt";
    static constexpr const char* name_gen_b_eta_ = "gen_b_eta";
    static constexpr const char* name_gen_b_phi_ = "gen_b_phi";
    static constexpr const char* name_gen_b_mass_ = "gen_b_mass";

    static constexpr const char* name_gen_bbar_pt_ = "gen_bbar_pt";
    static constexpr const char* name_gen_bbar_eta_ = "gen_bbar_eta";
    static constexpr const char* name_gen_bbar_phi_ = "gen_bbar_phi";
    static constexpr const char* name_gen_bbar_mass_ = "gen_bbar_mass";

    static constexpr const char* name_gen_nu_pt_ = "gen_nu_pt";
    static constexpr const char* name_gen_nu_eta_ = "gen_nu_eta";
    static constexpr const char* name_gen_nu_phi_ = "gen_nu_phi";
    static constexpr const char* name_gen_nu_mass_ = "gen_nu_mass";

    static constexpr const char* name_gen_nubar_pt_ = "gen_nubar_pt";
    static constexpr const char* name_gen_nubar_eta_ = "gen_nubar_eta";
    static constexpr const char* name_gen_nubar_phi_ = "gen_nubar_phi";
    static constexpr const char* name_gen_nubar_mass_ = "gen_nubar_mass";

    static constexpr const char* name_gen_ttbar_pt_ = "gen_ttbar_pt";
    static constexpr const char* name_gen_ttbar_phi_ = "gen_ttbar_phi";
    static constexpr const char* name_gen_ttbar_rapidity_ = "gen_ttbar_rapidity";
    //    static constexpr const char* name_gen_ttbar_arapidity_ = "gen_ttbar_arapidity";
    static constexpr const char* name_gen_ttbar_eta_ = "gen_ttbar_eta";
    static constexpr const char* name_gen_ttbar_delta_phi_ = "gen_ttbar_delta_phi";
    static constexpr const char* name_gen_ttbar_delta_eta_ = "gen_ttbar_delta_eta";
    static constexpr const char* name_gen_ttbar_delta_rapidity_ = "gen_ttbar_delta_rapidity";
    static constexpr const char* name_gen_ttbar_mass_ = "gen_ttbar_mass";

    static constexpr const char* name_gen_llbar_pt_ = "gen_llbar_pt";
    static constexpr const char* name_gen_llbar_phi_ = "gen_llbar_phi";
    static constexpr const char* name_gen_llbar_rapidity_ = "gen_llbar_rapidity";
    //    static constexpr const char* name_gen_llbar_arapidity_ = "gen_llbar_arapidity";
    static constexpr const char* name_gen_llbar_delta_phi_ = "gen_llbar_delta_phi";
    static constexpr const char* name_gen_llbar_delta_eta_ = "gen_llbar_delta_eta";
    static constexpr const char* name_gen_llbar_delta_rapidity_ = "gen_llbar_delta_rapidity";
    static constexpr const char* name_gen_llbar_mass_ = "gen_llbar_mass";

    //    static constexpr const char* name_gen_bbbar_pt_ = "gen_bbbar_pt";
    //    static constexpr const char* name_gen_bbbar_phi_ = "gen_bbbar_phi";
    //    static constexpr const char* name_gen_bbbar_rapidity_ = "gen_bbbar_rapidity";
    //    static constexpr const char* name_gen_bbbar_arapidity_ = "gen_bbbar_arapidity";
    //    static constexpr const char* name_gen_bbbar_delta_phi_ = "gen_bbbar_delta_phi";
    //    static constexpr const char* name_gen_bbbar_delta_eta_ = "gen_bbbar_delta_eta";
    //    static constexpr const char* name_gen_bbbar_delta_rapidity_ = "gen_bbbar_delta_rapidity";
    //    static constexpr const char* name_gen_bbbar_mass_ = "gen_bbbar_mass";

    //    static constexpr const char* name_gen_nunubar_pt_ = "gen_nunubar_pt";
    //    static constexpr const char* name_gen_nunubar_phi_ = "gen_nunubar_phi";
    //    static constexpr const char* name_gen_nunubar_rapidity_ = "gen_nunubar_rapidity";
    //    static constexpr const char* name_gen_nunubar_arapidity_ = "gen_nunubar_arapidity";
    //    static constexpr const char* name_gen_nunubar_delta_phi_ = "gen_nunubar_delta_phi";
    //    static constexpr const char* name_gen_nunubar_delta_eta_ = "gen_nunubar_delta_eta";
    //    static constexpr const char* name_gen_nunubar_delta_rapidity_ = "gen_nunubar_delta_rapidity";
    //    static constexpr const char* name_gen_nunubar_mass_ = "gen_nunubar_mass";

    static constexpr const char* name_gen_MT2_ = "gen_MT2";

    static constexpr const char* name_gen_all_mass_ = "gen_all_mass";
    static constexpr const char* name_gen_r_mass_ = "gen_r_mass";
    static constexpr const char* name_gen_jet_multiplicity_ = "gen_jet_multiplicity";
    static constexpr const char* name_gen_x1_ = "gen_x1";
    static constexpr const char* name_gen_x2_ = "gen_x2";
    static constexpr const char* name_gen_production_mode_ = "gen_production_mode";

    //extraJet variables
    static constexpr const char* name_gen_n_extraJets_iso08_ = "gen_n_extraJets_iso08";

    //ttbar frame variables
    static constexpr const char* name_gen_top_scatteringangle_ttbarframe_ = "gen_top_scatteringangle_ttbarframe";

    //    static constexpr const char* name_gen_top_pt_ttbarframe_ = "gen_top_pt_ttbarframe";
    //    static constexpr const char* name_gen_top_phi_ttbarframe_ = "gen_top_phi_ttbarframe";
    //    static constexpr const char* name_gen_top_rapidity_ttbarframe_ = "gen_top_rapidity_ttbarframe";
    //    static constexpr const char* name_gen_top_arapidity_ttbarframe_ = "gen_top_arapidity_ttbarframe";
    //    static constexpr const char* name_gen_tbar_pt_ttbarframe_ = "gen_tbar_pt_ttbarframe";
    //    static constexpr const char* name_gen_tbar_phi_ttbarframe_ = "gen_tbar_phi_ttbarframe";
    //    static constexpr const char* name_gen_tbar_rapidity_ttbarframe_ = "gen_tbar_rapidity_ttbarframe";
    //    static constexpr const char* name_gen_tbar_arapidity_ttbarframe_ = "gen_tbar_arapidity_ttbarframe";

    //top and antitop frame variables
    //    static constexpr const char* name_gen_l_pt_antitopframe_ = "gen_l_pt_antitopframe";
    //    static constexpr const char* name_gen_l_eta_antitopframe_ = "gen_l_eta_antitopframe";
    //    static constexpr const char* name_gen_l_phi_antitopframe_ = "gen_l_phi_antitopframe";
    //    static constexpr const char* name_gen_lbar_pt_topframe_ = "gen_lbar_pt_topframe";
    //    static constexpr const char* name_gen_lbar_eta_topframe_ = "gen_lbar_eta_topframe";
    //    static constexpr const char* name_gen_lbar_phi_topframe_ = "gen_lbar_phi_topframe";

    //Spin correlation variables
    static constexpr const char* name_gen_b1k_ = "gen_b1k";
    static constexpr const char* name_gen_b2k_ = "gen_b2k";
    static constexpr const char* name_gen_b1j_ = "gen_b1j";
    static constexpr const char* name_gen_b2j_ = "gen_b2j";
    static constexpr const char* name_gen_b1r_ = "gen_b1r";
    static constexpr const char* name_gen_b2r_ = "gen_b2r";
    static constexpr const char* name_gen_b1q_ = "gen_b1q";
    static constexpr const char* name_gen_b2q_ = "gen_b2q";
    static constexpr const char* name_gen_b1n_ = "gen_b1n";
    static constexpr const char* name_gen_b2n_ = "gen_b2n";

    //static constexpr const char* name_gen_b_Pkk_ = "gen_b_Pkk";
    //static constexpr const char* name_gen_b_Mkk_ = "gen_b_Mkk";
    //static constexpr const char* name_gen_b_Pjj_ = "gen_b_Pjj";
    //static constexpr const char* name_gen_b_Mjj_ = "gen_b_Mjj";
    //static constexpr const char* name_gen_b_Prr_ = "gen_b_Prr";
    //static constexpr const char* name_gen_b_Mrr_ = "gen_b_Mrr";
    //static constexpr const char* name_gen_b_Pqq_ = "gen_b_Pqq";
    //static constexpr const char* name_gen_b_Mqq_ = "gen_b_Mqq";
    //static constexpr const char* name_gen_b_Pnn_ = "gen_b_Pnn";
    //static constexpr const char* name_gen_b_Mnn_ = "gen_b_Mnn";

    static constexpr const char* name_gen_c_kk_ = "gen_c_kk";
    static constexpr const char* name_gen_c_rr_ = "gen_c_rr";
    static constexpr const char* name_gen_c_nn_ = "gen_c_nn";
    static constexpr const char* name_gen_c_kj_ = "gen_c_kj";
    static constexpr const char* name_gen_c_rq_ = "gen_c_rq";

    static constexpr const char* name_gen_c_rk_ = "gen_c_rk";
    static constexpr const char* name_gen_c_kr_ = "gen_c_kr";
    static constexpr const char* name_gen_c_nr_ = "gen_c_nr";
    static constexpr const char* name_gen_c_rn_ = "gen_c_rn";
    static constexpr const char* name_gen_c_nk_ = "gen_c_nk";
    static constexpr const char* name_gen_c_kn_ = "gen_c_kn";

    static constexpr const char* name_gen_c_rj_ = "gen_c_rj";
    static constexpr const char* name_gen_c_jr_ = "gen_c_jr";
    //static constexpr const char* name_gen_c_qk_ = "gen_c_qk";
    //static constexpr const char* name_gen_c_kq_ = "gen_c_kq";

    //static constexpr const char* name_gen_c_qj_ = "gen_c_qj";
    //static constexpr const char* name_gen_c_jq_ = "gen_c_jq";
    //static constexpr const char* name_gen_c_nq_ = "gen_c_nq";
    //static constexpr const char* name_gen_c_qn_ = "gen_c_qn";
    //static constexpr const char* name_gen_c_nj_ = "gen_c_nj";
    //static constexpr const char* name_gen_c_jn_ = "gen_c_jn";

    //static constexpr const char* name_gen_c_Prk_ = "gen_c_Prk";
    //static constexpr const char* name_gen_c_Mrk_ = "gen_c_Mrk";
    //static constexpr const char* name_gen_c_Pnr_ = "gen_c_Pnr";
    //static constexpr const char* name_gen_c_Mnr_ = "gen_c_Mnr";
    //static constexpr const char* name_gen_c_Pnk_ = "gen_c_Pnk";
    //static constexpr const char* name_gen_c_Mnk_ = "gen_c_Mnk";

    //static constexpr const char* name_gen_c_Pqj_ = "gen_c_Pqj";
    //static constexpr const char* name_gen_c_Mqj_ = "gen_c_Mqj";
    //static constexpr const char* name_gen_c_Pnq_ = "gen_c_Pnq";
    //static constexpr const char* name_gen_c_Mnq_ = "gen_c_Mnq";
    //static constexpr const char* name_gen_c_Pnj_ = "gen_c_Pnj";
    //static constexpr const char* name_gen_c_Mnj_ = "gen_c_Mnj";

    //static constexpr const char* name_gen_c_han_ = "gen_c_han";
    //static constexpr const char* name_gen_c_sca_ = "gen_c_sca";
    //static constexpr const char* name_gen_c_tra_ = "gen_c_tra";
    //static constexpr const char* name_gen_c_lij_ = "gen_c_lij";
    //static constexpr const char* name_gen_c_liq_ = "gen_c_liq";

    //static constexpr const char* name_gen_c_rkP_ = "gen_c_rkP";
    //static constexpr const char* name_gen_c_rkM_ = "gen_c_rkM";
    //static constexpr const char* name_gen_c_nrP_ = "gen_c_nrP";
    //static constexpr const char* name_gen_c_nrM_ = "gen_c_nrM";
    //static constexpr const char* name_gen_c_nkP_ = "gen_c_nkP";
    //static constexpr const char* name_gen_c_nkM_ = "gen_c_nkM";

    //static constexpr const char* name_gen_c_qjP_ = "gen_c_qjP";
    //static constexpr const char* name_gen_c_qjM_ = "gen_c_qjM";
    //static constexpr const char* name_gen_c_nqP_ = "gen_c_nqP";
    //static constexpr const char* name_gen_c_nqM_ = "gen_c_nqM";
    //static constexpr const char* name_gen_c_njP_ = "gen_c_njP";
    //static constexpr const char* name_gen_c_njM_ = "gen_c_njM";

    static constexpr const char* name_gen_ll_cHel_ = "gen_ll_cHel";
    static constexpr const char* name_gen_ll_cLab_ = "gen_ll_cLab";
    static constexpr const char* name_gen_ll_kNorm_ = "gen_ll_kNorm";
    static constexpr const char* name_gen_ll_rNorm_ = "gen_ll_rNorm";


    static constexpr const char* name_entry_ = "entry";
    static constexpr const char* name_isTopGen_ = "isTopGen";
    static constexpr const char* name_TopDecayMode_ = "TopDecayMode";
    static constexpr const char* name_isKinReco_ = "isKinReco";
    static constexpr const char* name_eventWeight_ = "eventWeight";
    static constexpr const char* name_trueLevelWeight_ = "trueLevelWeight";
    
};






#endif
