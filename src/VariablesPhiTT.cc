#include <iostream>
#include <cstdlib>
#include <cmath>

#include "VariablesPhiTT.h"
#include "analysisStructs.h"
#include "../../common/include/AnalysisBase.h"
#include "../../common/include/classes.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/KinematicReconstructionSolution.h"


// ----------------------------------------- Methods for VariablesPhiTT -----------------------------------------------------------
VariablesPhiTT::VariablesPhiTT():
VariablesBase(),

sol_pLep_pt( VariableFloat("sol_pLep_pt") ),
sol_pLep_eta( VariableFloat("sol_pLep_eta") ),
sol_pLep_phi( VariableFloat("sol_pLep_phi") ),
sol_aLep_pt( VariableFloat("sol_aLep_pt") ),
sol_aLep_eta( VariableFloat("sol_aLep_eta") ),
sol_aLep_phi( VariableFloat("sol_aLep_phi") ),

sol_pNeu_pt( VariableFloat("sol_pNeu_pt") ),
sol_pNeu_eta( VariableFloat("sol_pNeu_eta") ),
sol_pNeu_phi( VariableFloat("sol_pNeu_phi") ),
sol_aNeu_pt( VariableFloat("sol_aNeu_pt") ),
sol_aNeu_eta( VariableFloat("sol_aNeu_eta") ),
sol_aNeu_phi( VariableFloat("sol_aNeu_phi") ),

sol_pBot_pt( VariableFloat("sol_pBot_pt") ),
sol_pBot_eta( VariableFloat("sol_pBot_eta") ),
sol_pBot_phi( VariableFloat("sol_pBot_phi") ),
sol_aBot_pt( VariableFloat("sol_aBot_pt") ),
sol_aBot_eta( VariableFloat("sol_aBot_eta") ),
sol_aBot_phi( VariableFloat("sol_aBot_phi") ),

sol_pTop_pt( VariableFloat("sol_pTop_pt") ),
sol_pTop_rap( VariableFloat("sol_pTop_rap") ),
sol_pTop_phi( VariableFloat("sol_pTop_phi") ),
sol_aTop_pt( VariableFloat("sol_aTop_pt") ),
sol_aTop_rap( VariableFloat("sol_aTop_rap") ),
sol_aTop_phi( VariableFloat("sol_aTop_phi") ),

sol_TT_m( VariableFloat("sol_TT_m") ),
sol_TT_pt( VariableFloat("sol_TT_pt") ),
sol_TT_rap( VariableFloat("sol_TT_rap") ),
sol_TT_phi( VariableFloat("sol_TT_phi") ),

sol_LL_dEta( VariableFloat("sol_LL_dEta") ),
sol_LL_dPhi( VariableFloat("sol_LL_dPhi") ),
sol_LL_dR( VariableFloat("sol_LL_dR") ),
sol_cLab( VariableFloat("sol_cLab") ),

sol_b1k( VariableFloat("sol_b1k") ),
sol_b1j( VariableFloat("sol_b1j") ),
sol_b1r( VariableFloat("sol_b1r") ),
sol_b1q( VariableFloat("sol_b1q") ),
sol_b1n( VariableFloat("sol_b1n") ),

sol_b2k( VariableFloat("sol_b2k") ),
sol_b2j( VariableFloat("sol_b2j") ),
sol_b2r( VariableFloat("sol_b2r") ),
sol_b2q( VariableFloat("sol_b2q") ),
sol_b2n( VariableFloat("sol_b2n") ),

//sol_bP_kk( VariableFloat("sol_bP_kk") ),
//sol_bM_kk( VariableFloat("sol_bM_kk") ),
//sol_bP_jj( VariableFloat("sol_bP_jj") ),
//sol_bM_jj( VariableFloat("sol_bM_jj") ),
//sol_bP_rr( VariableFloat("sol_bP_rr") ),
//sol_bM_rr( VariableFloat("sol_bM_rr") ),
//sol_bP_qq( VariableFloat("sol_bP_qq") ),
//sol_bM_qq( VariableFloat("sol_bM_qq") ),
//sol_bP_nn( VariableFloat("sol_bP_nn") ),
//sol_bM_nn( VariableFloat("sol_bM_nn") ),

sol_ckk( VariableFloat("sol_ckk") ),
sol_crr( VariableFloat("sol_crr") ),
sol_cnn( VariableFloat("sol_cnn") ),
sol_ckj( VariableFloat("sol_ckj") ),
sol_crq( VariableFloat("sol_crq") ),

sol_crk( VariableFloat("sol_crk") ),
sol_ckr( VariableFloat("sol_ckr") ),
sol_cnr( VariableFloat("sol_cnr") ),
sol_crn( VariableFloat("sol_crn") ),
sol_cnk( VariableFloat("sol_cnk") ),
sol_ckn( VariableFloat("sol_ckn") ),

sol_cqk( VariableFloat("sol_cqk") ),
sol_ckq( VariableFloat("sol_ckq") ),
sol_crj( VariableFloat("sol_crj") ),
sol_cjr( VariableFloat("sol_cjr") ),

sol_cqj( VariableFloat("sol_cqj") ),
sol_cjq( VariableFloat("sol_cjq") ),
sol_cnq( VariableFloat("sol_cnq") ),
sol_cqn( VariableFloat("sol_cqn") ),
sol_cnj( VariableFloat("sol_cnj") ),
sol_cjn( VariableFloat("sol_cjn") ),

//sol_cP_rk( VariableFloat("sol_cP_rk") ),
//sol_cM_rk( VariableFloat("sol_cM_rk") ),
//sol_cP_nr( VariableFloat("sol_cP_nr") ),
//sol_cM_nr( VariableFloat("sol_cM_nr") ),
//sol_cP_nk( VariableFloat("sol_cP_nk") ),
//sol_cM_nk( VariableFloat("sol_cM_nk") ),

sol_ll_dEta( VariableFloat("sol_ll_dEta") ),
sol_ll_dPhi( VariableFloat("sol_ll_dPhi") ),
sol_ll_dR( VariableFloat("sol_ll_dR") ),
sol_cHel( VariableFloat("sol_cHel") ),

sol_kNorm( VariableFloat("sol_kNorm") ),
sol_rNorm( VariableFloat("sol_rNorm") ),
sol_nNorm( VariableFloat("sol_nNorm") ),
sol_jNorm( VariableFloat("sol_jNorm") ),
sol_qNorm( VariableFloat("sol_qNorm") ),

sol_phi0( VariableFloat("sol_phi0") ),
sol_phi1( VariableFloat("sol_phi1") ),

sol_TT_sph( VariableFloat("sol_TT_sph") ),
sol_TT_apl( VariableFloat("sol_TT_apl") ),
sol_TT_cir( VariableFloat("sol_TT_cir") ),
sol_TT_iso( VariableFloat("sol_TT_iso") ),
sol_TT_C( VariableFloat("sol_TT_C") ),
sol_TT_D( VariableFloat("sol_TT_D") ),
sol_TT_spt( VariableFloat("sol_TT_spt") ),
sol_TT_H0( VariableFloat("sol_TT_H0") ),
sol_TT_H1( VariableFloat("sol_TT_H1") ),
sol_TT_H2( VariableFloat("sol_TT_H2") ),
sol_TT_H3( VariableFloat("sol_TT_H3") ),
sol_TT_H4( VariableFloat("sol_TT_H4") ),
sol_TT_R1( VariableFloat("sol_TT_R1") ),
sol_TT_R2( VariableFloat("sol_TT_R2") ),
sol_TT_R3( VariableFloat("sol_TT_R3") ),
sol_TT_R4( VariableFloat("sol_TT_R4") ),

gen_pLep_pt( VariableFloat("gen_pLep_pt") ),
gen_pLep_eta( VariableFloat("gen_pLep_eta") ),
gen_pLep_phi( VariableFloat("gen_pLep_phi") ),
gen_aLep_pt( VariableFloat("gen_aLep_pt") ),
gen_aLep_eta( VariableFloat("gen_aLep_eta") ),
gen_aLep_phi( VariableFloat("gen_aLep_phi") ),

gen_pNeu_pt( VariableFloat("gen_pNeu_pt") ),
gen_pNeu_eta( VariableFloat("gen_pNeu_eta") ),
gen_pNeu_phi( VariableFloat("gen_pNeu_phi") ),
gen_aNeu_pt( VariableFloat("gen_aNeu_pt") ),
gen_aNeu_eta( VariableFloat("gen_aNeu_eta") ),
gen_aNeu_phi( VariableFloat("gen_aNeu_phi") ),

gen_pBot_pt( VariableFloat("gen_pBot_pt") ),
gen_pBot_eta( VariableFloat("gen_pBot_eta") ),
gen_pBot_phi( VariableFloat("gen_pBot_phi") ),
gen_aBot_pt( VariableFloat("gen_aBot_pt") ),
gen_aBot_eta( VariableFloat("gen_aBot_eta") ),
gen_aBot_phi( VariableFloat("gen_aBot_phi") ),

gen_pTop_pt( VariableFloat("gen_pTop_pt") ),
gen_pTop_rap( VariableFloat("gen_pTop_rap") ),
gen_pTop_phi( VariableFloat("gen_pTop_phi") ),
gen_aTop_pt( VariableFloat("gen_aTop_pt") ),
gen_aTop_rap( VariableFloat("gen_aTop_rap") ),
gen_aTop_phi( VariableFloat("gen_aTop_phi") ),

gen_TT_m( VariableFloat("gen_TT_m") ),
gen_TT_pt( VariableFloat("gen_TT_pt") ),
gen_TT_rap( VariableFloat("gen_TT_rap") ),
gen_TT_phi( VariableFloat("gen_TT_phi") ),

gen_LL_dEta( VariableFloat("gen_LL_dEta") ),
gen_LL_dPhi( VariableFloat("gen_LL_dPhi") ),
gen_LL_dR( VariableFloat("gen_LL_dR") ),
gen_cLab( VariableFloat("gen_cLab") ),

gen_b1k( VariableFloat("gen_b1k") ),
gen_b1j( VariableFloat("gen_b1j") ),
gen_b1r( VariableFloat("gen_b1r") ),
gen_b1q( VariableFloat("gen_b1q") ),
gen_b1n( VariableFloat("gen_b1n") ),

gen_b2k( VariableFloat("gen_b2k") ),
gen_b2j( VariableFloat("gen_b2j") ),
gen_b2r( VariableFloat("gen_b2r") ),
gen_b2q( VariableFloat("gen_b2q") ),
gen_b2n( VariableFloat("gen_b2n") ),

//gen_bP_kk( VariableFloat("gen_bP_kk") ),
//gen_bM_kk( VariableFloat("gen_bM_kk") ),

//gen_bP_jj( VariableFloat("gen_bP_jj") ),
//gen_bM_jj( VariableFloat("gen_bM_jj") ),

//gen_bP_rr( VariableFloat("gen_bP_rr") ),
//gen_bM_rr( VariableFloat("gen_bM_rr") ),

//gen_bP_qq( VariableFloat("gen_bP_qq") ),
//gen_bM_qq( VariableFloat("gen_bM_qq") ),

//gen_bP_nn( VariableFloat("gen_bP_nn") ),
//gen_bM_nn( VariableFloat("gen_bM_nn") ),

gen_ckk( VariableFloat("gen_ckk") ),
gen_crr( VariableFloat("gen_crr") ),
gen_cnn( VariableFloat("gen_cnn") ),
gen_ckj( VariableFloat("gen_ckj") ),
gen_crq( VariableFloat("gen_crq") ),

gen_crk( VariableFloat("gen_crk") ),
gen_ckr( VariableFloat("gen_ckr") ),
gen_cnr( VariableFloat("gen_cnr") ),
gen_crn( VariableFloat("gen_crn") ),
gen_cnk( VariableFloat("gen_cnk") ),
gen_ckn( VariableFloat("gen_ckn") ),

gen_crj( VariableFloat("gen_crj") ),
gen_cjr( VariableFloat("gen_cjr") ),
gen_cqk( VariableFloat("gen_cqk") ),
gen_ckq( VariableFloat("gen_ckq") ),

gen_cqj( VariableFloat("gen_cqj") ),
gen_cjq( VariableFloat("gen_cjq") ),
gen_cnq( VariableFloat("gen_cnq") ),
gen_cqn( VariableFloat("gen_cqn") ),
gen_cnj( VariableFloat("gen_cnj") ),
gen_cjn( VariableFloat("gen_cjn") ),

//gen_cP_rk( VariableFloat("gen_cP_rk") ),
//gen_cM_rk( VariableFloat("gen_cM_rk") ),

//gen_cP_nr( VariableFloat("gen_cP_nr") ),
//gen_cM_nr( VariableFloat("gen_cM_nr") ),

//gen_cP_nk( VariableFloat("gen_cP_nk") ),
//gen_cM_nk( VariableFloat("gen_cM_nk") ),

gen_ll_dEta( VariableFloat("gen_ll_dEta") ),
gen_ll_dPhi( VariableFloat("gen_ll_dPhi") ),
gen_ll_dR( VariableFloat("gen_ll_dR") ),
gen_cHel( VariableFloat("gen_cHel") ),

gen_kNorm( VariableFloat("gen_kNorm") ),
gen_rNorm( VariableFloat("gen_rNorm") ),
gen_nNorm( VariableFloat("gen_nNorm") ),
gen_jNorm( VariableFloat("gen_jNorm") ),
gen_qNorm( VariableFloat("gen_qNorm") ),


gen_phi0( VariableFloat("gen_phi0") ),
gen_phi1( VariableFloat("gen_phi1") ),

gen_TT_sph( VariableFloat("gen_TT_sph") ),
gen_TT_apl( VariableFloat("gen_TT_apl") ),
gen_TT_cir( VariableFloat("gen_TT_cir") ),
gen_TT_iso( VariableFloat("gen_TT_iso") ),
gen_TT_C( VariableFloat("gen_TT_C") ),
gen_TT_D( VariableFloat("gen_TT_D") ),
gen_TT_spt( VariableFloat("gen_TT_spt") ),
gen_TT_H0( VariableFloat("gen_TT_H0") ),
gen_TT_H1( VariableFloat("gen_TT_H1") ),
gen_TT_H2( VariableFloat("gen_TT_H2") ),
gen_TT_H3( VariableFloat("gen_TT_H3") ),
gen_TT_H4( VariableFloat("gen_TT_H4") ),
gen_TT_R1( VariableFloat("gen_TT_R1") ),
gen_TT_R2( VariableFloat("gen_TT_R2") ),
gen_TT_R3( VariableFloat("gen_TT_R3") ),
gen_TT_R4( VariableFloat("gen_TT_R4") ),

pBotMatchSolGen( VariableInt("pBotMatchSolGen") ),
aBotMatchSolGen( VariableInt("aBotMatchSolGen") ),

nRun( VariableInt("nRun") ),
nLumi( VariableInt("nLumi") ),
nEvt( VariableLong64("nEvt") ),
isKinReco( VariableInt("isKinReco") ),
isTopGen( VariableInt("isTopGen") ),
topProdMode( VariableInt("topProdMode") ),
recoWeight( VariableFloat("recoWeight") ),
genWeight( VariableFloat("genWeight") )
{}



VariablesPhiTT::VariablesPhiTT(const EventMetadata& eventMetadata,
                               const RecoObjects&,
                               const CommonGenObjects&, 
                               const TopGenObjects& topGenObjects,
                               const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                               const ttbar::RecoObjectIndices&, const ttbar::GenObjectIndices&, 
                               const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                               const double& weight):
VariablesBase(weight),

sol_pLep_pt( VariableFloat("sol_pLep_pt") ),
sol_pLep_eta( VariableFloat("sol_pLep_eta") ),
sol_pLep_phi( VariableFloat("sol_pLep_phi") ),
sol_aLep_pt( VariableFloat("sol_aLep_pt") ),
sol_aLep_eta( VariableFloat("sol_aLep_eta") ),
sol_aLep_phi( VariableFloat("sol_aLep_phi") ),

sol_pNeu_pt( VariableFloat("sol_pNeu_pt") ),
sol_pNeu_eta( VariableFloat("sol_pNeu_eta") ),
sol_pNeu_phi( VariableFloat("sol_pNeu_phi") ),
sol_aNeu_pt( VariableFloat("sol_aNeu_pt") ),
sol_aNeu_eta( VariableFloat("sol_aNeu_eta") ),
sol_aNeu_phi( VariableFloat("sol_aNeu_phi") ),

sol_pBot_pt( VariableFloat("sol_pBot_pt") ),
sol_pBot_eta( VariableFloat("sol_pBot_eta") ),
sol_pBot_phi( VariableFloat("sol_pBot_phi") ),
sol_aBot_pt( VariableFloat("sol_aBot_pt") ),
sol_aBot_eta( VariableFloat("sol_aBot_eta") ),
sol_aBot_phi( VariableFloat("sol_aBot_phi") ),

sol_pTop_pt( VariableFloat("sol_pTop_pt") ),
sol_pTop_rap( VariableFloat("sol_pTop_rap") ),
sol_pTop_phi( VariableFloat("sol_pTop_phi") ),
sol_aTop_pt( VariableFloat("sol_aTop_pt") ),
sol_aTop_rap( VariableFloat("sol_aTop_rap") ),
sol_aTop_phi( VariableFloat("sol_aTop_phi") ),

sol_TT_m( VariableFloat("sol_TT_m") ),
sol_TT_pt( VariableFloat("sol_TT_pt") ),
sol_TT_rap( VariableFloat("sol_TT_rap") ),
sol_TT_phi( VariableFloat("sol_TT_phi") ),

sol_LL_dEta( VariableFloat("sol_LL_dEta") ),
sol_LL_dPhi( VariableFloat("sol_LL_dPhi") ),
sol_LL_dR( VariableFloat("sol_LL_dR") ),
sol_cLab( VariableFloat("sol_cLab") ),

sol_b1k( VariableFloat("sol_b1k") ),
sol_b1j( VariableFloat("sol_b1j") ),
sol_b1r( VariableFloat("sol_b1r") ),
sol_b1q( VariableFloat("sol_b1q") ),
sol_b1n( VariableFloat("sol_b1n") ),

sol_b2k( VariableFloat("sol_b2k") ),
sol_b2j( VariableFloat("sol_b2j") ),
sol_b2r( VariableFloat("sol_b2r") ),
sol_b2q( VariableFloat("sol_b2q") ),
sol_b2n( VariableFloat("sol_b2n") ),

//sol_bP_kk( VariableFloat("sol_bP_kk") ),
//sol_bM_kk( VariableFloat("sol_bM_kk") ),

//sol_bP_jj( VariableFloat("sol_bP_jj") ),
//sol_bM_jj( VariableFloat("sol_bM_jj") ),
//sol_bP_rr( VariableFloat("sol_bP_rr") ),
//sol_bM_rr( VariableFloat("sol_bM_rr") ),
//sol_bP_qq( VariableFloat("sol_bP_qq") ),
//sol_bM_qq( VariableFloat("sol_bM_qq") ),
//sol_bP_nn( VariableFloat("sol_bP_nn") ),
//sol_bM_nn( VariableFloat("sol_bM_nn") ),

sol_ckk( VariableFloat("sol_ckk") ),
sol_crr( VariableFloat("sol_crr") ),
sol_cnn( VariableFloat("sol_cnn") ),
sol_ckj( VariableFloat("sol_ckj") ),
sol_crq( VariableFloat("sol_crq") ),

sol_crk( VariableFloat("sol_crk") ),
sol_ckr( VariableFloat("sol_ckr") ),
sol_cnr( VariableFloat("sol_cnr") ),
sol_crn( VariableFloat("sol_crn") ),
sol_cnk( VariableFloat("sol_cnk") ),
sol_ckn( VariableFloat("sol_ckn") ),

sol_crj( VariableFloat("sol_crj") ),
sol_cjr( VariableFloat("sol_cjr") ),
sol_cqk( VariableFloat("sol_cqk") ),
sol_ckq( VariableFloat("sol_ckq") ),

sol_cqj( VariableFloat("sol_cqj") ),
sol_cjq( VariableFloat("sol_cjq") ),
sol_cnq( VariableFloat("sol_cnq") ),
sol_cqn( VariableFloat("sol_cqn") ),
sol_cnj( VariableFloat("sol_cnj") ),
sol_cjn( VariableFloat("sol_cjn") ),

//sol_cP_rk( VariableFloat("sol_cP_rk") ),
//sol_cM_rk( VariableFloat("sol_cM_rk") ),
//sol_cP_nr( VariableFloat("sol_cP_nr") ),
//sol_cM_nr( VariableFloat("sol_cM_nr") ),
//sol_cP_nk( VariableFloat("sol_cP_nk") ),
//sol_cM_nk( VariableFloat("sol_cM_nk") ),

sol_ll_dEta( VariableFloat("sol_ll_dEta") ),
sol_ll_dPhi( VariableFloat("sol_ll_dPhi") ),
sol_ll_dR( VariableFloat("sol_ll_dR") ),
sol_cHel( VariableFloat("sol_cHel") ),

sol_kNorm( VariableFloat("sol_kNorm") ),
sol_rNorm( VariableFloat("sol_rNorm") ),
sol_nNorm( VariableFloat("sol_nNorm") ),
sol_jNorm( VariableFloat("sol_jNorm") ),
sol_qNorm( VariableFloat("sol_qNorm") ),


sol_phi0( VariableFloat("sol_phi0") ),
sol_phi1( VariableFloat("sol_phi1") ),

sol_TT_sph( VariableFloat("sol_TT_sph") ),
sol_TT_apl( VariableFloat("sol_TT_apl") ),
sol_TT_cir( VariableFloat("sol_TT_cir") ),
sol_TT_iso( VariableFloat("sol_TT_iso") ),
sol_TT_C( VariableFloat("sol_TT_C") ),
sol_TT_D( VariableFloat("sol_TT_D") ),
sol_TT_spt( VariableFloat("sol_TT_spt") ),
sol_TT_H0( VariableFloat("sol_TT_H0") ),
sol_TT_H1( VariableFloat("sol_TT_H1") ),
sol_TT_H2( VariableFloat("sol_TT_H2") ),
sol_TT_H3( VariableFloat("sol_TT_H3") ),
sol_TT_H4( VariableFloat("sol_TT_H4") ),
sol_TT_R1( VariableFloat("sol_TT_R1") ),
sol_TT_R2( VariableFloat("sol_TT_R2") ),
sol_TT_R3( VariableFloat("sol_TT_R3") ),
sol_TT_R4( VariableFloat("sol_TT_R4") ),

gen_pLep_pt( VariableFloat("gen_pLep_pt") ),
gen_pLep_eta( VariableFloat("gen_pLep_eta") ),
gen_pLep_phi( VariableFloat("gen_pLep_phi") ),
gen_aLep_pt( VariableFloat("gen_aLep_pt") ),
gen_aLep_eta( VariableFloat("gen_aLep_eta") ),
gen_aLep_phi( VariableFloat("gen_aLep_phi") ),

gen_pNeu_pt( VariableFloat("gen_pNeu_pt") ),
gen_pNeu_eta( VariableFloat("gen_pNeu_eta") ),
gen_pNeu_phi( VariableFloat("gen_pNeu_phi") ),
gen_aNeu_pt( VariableFloat("gen_aNeu_pt") ),
gen_aNeu_eta( VariableFloat("gen_aNeu_eta") ),
gen_aNeu_phi( VariableFloat("gen_aNeu_phi") ),

gen_pBot_pt( VariableFloat("gen_pBot_pt") ),
gen_pBot_eta( VariableFloat("gen_pBot_eta") ),
gen_pBot_phi( VariableFloat("gen_pBot_phi") ),
gen_aBot_pt( VariableFloat("gen_aBot_pt") ),
gen_aBot_eta( VariableFloat("gen_aBot_eta") ),
gen_aBot_phi( VariableFloat("gen_aBot_phi") ),

gen_pTop_pt( VariableFloat("gen_pTop_pt") ),
gen_pTop_rap( VariableFloat("gen_pTop_rap") ),
gen_pTop_phi( VariableFloat("gen_pTop_phi") ),
gen_aTop_pt( VariableFloat("gen_aTop_pt") ),
gen_aTop_rap( VariableFloat("gen_aTop_rap") ),
gen_aTop_phi( VariableFloat("gen_aTop_phi") ),

gen_TT_m( VariableFloat("gen_TT_m") ),
gen_TT_pt( VariableFloat("gen_TT_pt") ),
gen_TT_rap( VariableFloat("gen_TT_rap") ),
gen_TT_phi( VariableFloat("gen_TT_phi") ),

gen_LL_dEta( VariableFloat("gen_LL_dEta") ),
gen_LL_dPhi( VariableFloat("gen_LL_dPhi") ),
gen_LL_dR( VariableFloat("gen_LL_dR") ),
gen_cLab( VariableFloat("gen_cLab") ),

gen_b1k( VariableFloat("gen_b1k") ),
gen_b1j( VariableFloat("gen_b1j") ),
gen_b1r( VariableFloat("gen_b1r") ),
gen_b1q( VariableFloat("gen_b1q") ),
gen_b1n( VariableFloat("gen_b1n") ),

gen_b2k( VariableFloat("gen_b2k") ),
gen_b2j( VariableFloat("gen_b2j") ),
gen_b2r( VariableFloat("gen_b2r") ),
gen_b2q( VariableFloat("gen_b2q") ),
gen_b2n( VariableFloat("gen_b2n") ),

//gen_bP_kk( VariableFloat("gen_bP_kk") ),
//gen_bM_kk( VariableFloat("gen_bM_kk") ),

//gen_bP_jj( VariableFloat("gen_bP_jj") ),
//gen_bM_jj( VariableFloat("gen_bM_jj") ),

//gen_bP_rr( VariableFloat("gen_bP_rr") ),
//gen_bM_rr( VariableFloat("gen_bM_rr") ),

//gen_bP_qq( VariableFloat("gen_bP_qq") ),
//gen_bM_qq( VariableFloat("gen_bM_qq") ),

//gen_bP_nn( VariableFloat("gen_bP_nn") ),
//gen_bM_nn( VariableFloat("gen_bM_nn") ),

gen_ckk( VariableFloat("gen_ckk") ),
gen_crr( VariableFloat("gen_crr") ),
gen_cnn( VariableFloat("gen_cnn") ),
gen_ckj( VariableFloat("gen_ckj") ),
gen_crq( VariableFloat("gen_crq") ),

gen_crk( VariableFloat("gen_crk") ),
gen_ckr( VariableFloat("gen_ckr") ),
gen_cnr( VariableFloat("gen_cnr") ),
gen_crn( VariableFloat("gen_crn") ),
gen_cnk( VariableFloat("gen_cnk") ),
gen_ckn( VariableFloat("gen_ckn") ),

gen_crj( VariableFloat("gen_crj") ),
gen_cjr( VariableFloat("gen_cjr") ),
gen_cqk( VariableFloat("gen_cqk") ),
gen_ckq( VariableFloat("gen_ckq") ),

gen_cqj( VariableFloat("gen_cqj") ),
gen_cjq( VariableFloat("gen_cjq") ),
gen_cnq( VariableFloat("gen_cnq") ),
gen_cqn( VariableFloat("gen_cqn") ),
gen_cnj( VariableFloat("gen_cnj") ),
gen_cjn( VariableFloat("gen_cjn") ),

//gen_cP_rk( VariableFloat("gen_cP_rk") ),
//gen_cM_rk( VariableFloat("gen_cM_rk") ),
//gen_cP_nr( VariableFloat("gen_cP_nr") ),
//gen_cM_nr( VariableFloat("gen_cM_nr") ),
//gen_cP_nk( VariableFloat("gen_cP_nk") ),
//gen_cM_nk( VariableFloat("gen_cM_nk") ),

gen_ll_dEta( VariableFloat("gen_ll_dEta") ),
gen_ll_dPhi( VariableFloat("gen_ll_dPhi") ),
gen_ll_dR( VariableFloat("gen_ll_dR") ),
gen_cHel( VariableFloat("gen_cHel") ),

gen_kNorm( VariableFloat("gen_kNorm") ),
gen_rNorm( VariableFloat("gen_rNorm") ),
gen_nNorm( VariableFloat("gen_nNorm") ),
gen_jNorm( VariableFloat("gen_jNorm") ),
gen_qNorm( VariableFloat("gen_qNorm") ),

gen_phi0( VariableFloat("gen_phi0") ),
gen_phi1( VariableFloat("gen_phi1") ),

gen_TT_sph( VariableFloat("gen_TT_sph") ),
gen_TT_apl( VariableFloat("gen_TT_apl") ),
gen_TT_cir( VariableFloat("gen_TT_cir") ),
gen_TT_iso( VariableFloat("gen_TT_iso") ),
gen_TT_C( VariableFloat("gen_TT_C") ),
gen_TT_D( VariableFloat("gen_TT_D") ),
gen_TT_spt( VariableFloat("gen_TT_spt") ),
gen_TT_H0( VariableFloat("gen_TT_H0") ),
gen_TT_H1( VariableFloat("gen_TT_H1") ),
gen_TT_H2( VariableFloat("gen_TT_H2") ),
gen_TT_H3( VariableFloat("gen_TT_H3") ),
gen_TT_H4( VariableFloat("gen_TT_H4") ),
gen_TT_R1( VariableFloat("gen_TT_R1") ),
gen_TT_R2( VariableFloat("gen_TT_R2") ),
gen_TT_R3( VariableFloat("gen_TT_R3") ),
gen_TT_R4( VariableFloat("gen_TT_R4") ),

pBotMatchSolGen( VariableInt("pBotMatchSolGen") ),
aBotMatchSolGen( VariableInt("aBotMatchSolGen") ),

nRun( VariableInt("nRun") ),
nLumi( VariableInt("nLumi") ),
nEvt( VariableLong64("nEvt") ),
isKinReco( VariableInt("isKinReco") ),
isTopGen( VariableInt("isTopGen") ),
topProdMode( VariableInt("topProdMode") ),
recoWeight( VariableFloat("recoWeight") ),
genWeight( VariableFloat("genWeight") )
{
    
  isKinReco.value_ = 0;
  isTopGen.value_ = 0;
  topProdMode.value_ = -1;
  pBotMatchSolGen.value_ = -1;
  aBotMatchSolGen.value_ = -1;

  nRun.value_ = eventMetadata.runNumber_;
  nLumi.value_ = eventMetadata.lumiBlock_;
  nEvt.value_ = eventMetadata.eventNumber_;
  recoWeight.value_ = recoLevelWeights.weight_;

  TLorentzVector sol_pLep, sol_aLep;
  TLorentzVector sol_pNeu, sol_aNeu;
  TLorentzVector sol_pBot, sol_aBot;
  TLorentzVector sol_pTop, sol_aTop;
  TLorentzVector sol_TT;

  TLorentzVector gen_pLep, gen_aLep;
  TLorentzVector gen_pNeu, gen_aNeu;
  TLorentzVector gen_pBot, gen_aBot;
  TLorentzVector gen_pTop, gen_aTop;
  TLorentzVector gen_TT;

  if ( kinematicReconstructionSolutions.numberOfSolutions() ) {
   
    isKinReco.value_ = 1;

    sol_pLep = common::LVtoTLV(kinematicReconstructionSolutions.solution().lepton());
    sol_aLep = common::LVtoTLV(kinematicReconstructionSolutions.solution().antiLepton());

    sol_pNeu = common::LVtoTLV(kinematicReconstructionSolutions.solution().neutrino());
    sol_aNeu = common::LVtoTLV(kinematicReconstructionSolutions.solution().antiNeutrino());

    sol_pBot = common::LVtoTLV(kinematicReconstructionSolutions.solution().bjet());
    sol_aBot = common::LVtoTLV(kinematicReconstructionSolutions.solution().antiBjet());

    sol_pTop = common::LVtoTLV(kinematicReconstructionSolutions.solution().top());
    sol_aTop = common::LVtoTLV(kinematicReconstructionSolutions.solution().antiTop());

    sol_TT = sol_pTop + sol_aTop;


    sol_pLep_pt.value_ = ifOk( sol_pLep.Pt() );
    sol_pLep_eta.value_ = ifOk( sol_pLep.Eta() );
    sol_pLep_phi.value_ = ifOk( sol_pLep.Phi() );
    sol_aLep_pt.value_ = ifOk( sol_aLep.Pt() );
    sol_aLep_eta.value_ = ifOk( sol_aLep.Eta() );
    sol_aLep_phi.value_ = ifOk( sol_aLep.Phi() );

    sol_pNeu_pt.value_ = ifOk( sol_pNeu.Pt() );
    sol_pNeu_eta.value_ = ifOk( sol_pNeu.Eta() );
    sol_pNeu_phi.value_ = ifOk( sol_pNeu.Phi() );
    sol_aNeu_pt.value_ = ifOk( sol_aNeu.Pt() );
    sol_aNeu_eta.value_ = ifOk( sol_aNeu.Eta() );
    sol_aNeu_phi.value_ = ifOk( sol_aNeu.Phi() );

    sol_pBot_pt.value_ = ifOk( sol_pBot.Pt() );
    sol_pBot_eta.value_ = ifOk( sol_pBot.Eta() );
    sol_pBot_phi.value_ = ifOk( sol_pBot.Phi() );
    sol_aBot_pt.value_ = ifOk( sol_aBot.Pt() );
    sol_aBot_eta.value_ = ifOk( sol_aBot.Eta() );
    sol_aBot_phi.value_ = ifOk( sol_aBot.Phi() );

    sol_pTop_pt.value_ = ifOk( sol_pTop.Pt() );
    sol_pTop_rap.value_ = ifOk( sol_pTop.Rapidity() );
    sol_pTop_phi.value_ = ifOk( sol_pTop.Phi() );
    sol_aTop_pt.value_ = ifOk( sol_aTop.Pt() );
    sol_aTop_rap.value_ = ifOk( sol_aTop.Rapidity() );
    sol_aTop_phi.value_ = ifOk( sol_aTop.Phi() );

    sol_TT_m.value_ = ifOk( sol_TT.M() );
    sol_TT_pt.value_ = ifOk( sol_TT.Pt() );
    sol_TT_rap.value_ = ifOk( sol_TT.Rapidity() );
    sol_TT_phi.value_ = ifOk( sol_TT.Phi() );

    // SpinCorr variables without boosting
    sol_LL_dEta.value_ = ifOk( sol_aLep.Eta() - sol_pLep.Eta() );
    sol_LL_dPhi.value_ = ifOk( sol_aLep.DeltaPhi( sol_pLep ) );
    sol_LL_dR.value_ = ifOk( sol_aLep.DeltaR( sol_pLep ) );
    sol_cLab.value_ = ifOk( sol_aLep.Vect().Unit().Dot( sol_pLep.Vect().Unit() ) );

    if(sol_TT_m.value_ == -999999  
  ||  sol_pTop_pt.value_   == -999999  ||  sol_pTop_rap.value_   == -999999  ||  sol_pTop_phi.value_   == -999999
  ||  sol_aTop_pt.value_   == -999999  ||  sol_aTop_rap.value_   == -999999  ||  sol_aTop_phi.value_   == -999999
  ||  sol_TT_pt.value_   == -999999  ||  sol_TT_rap.value_   == -999999  ||   sol_TT_phi.value_  == -999999 
  ||  sol_aBot_pt.value_   == -999999  ||  sol_aBot_eta.value_   == -999999  ||   sol_aBot_phi.value_  == -999999 
  ||  sol_pBot_pt.value_   == -999999  ||  sol_pBot_eta.value_   == -999999  ||   sol_pBot_phi.value_  == -999999 
  ||  sol_pNeu_pt.value_   == -999999  ||  sol_pNeu_eta.value_   == -999999  ||   sol_pNeu_phi.value_  == -999999
  ||  sol_aNeu_pt.value_   == -999999  ||  sol_aNeu_eta.value_   == -999999  ||   sol_aNeu_phi.value_  == -999999 
  ||  sol_aLep_pt.value_   == -999999  ||  sol_aLep_eta.value_   == -999999  ||   sol_aLep_phi.value_  == -999999 
  ||  sol_pLep_pt.value_   == -999999  ||  sol_pLep_eta.value_   == -999999  ||   sol_pLep_phi.value_  == -999999 ) {

      kinematicReconstructionSolutions.solution().print();

    }

    compSpinCorr(sol_pLep, sol_aLep, sol_pTop, sol_aTop,
                 sol_b1k, sol_b1j, sol_b1r, sol_b1q, sol_b1n,
                 sol_b2k, sol_b2j, sol_b2r, sol_b2q, sol_b2n,
                 //sol_bP_kk, sol_bM_kk, sol_bP_jj, sol_bM_jj,
                 //sol_bP_rr, sol_bM_rr, sol_bP_qq, sol_bM_qq,
                 //sol_bP_nn, sol_bM_nn,
                 sol_ckk, sol_crr, sol_cnn,
		 sol_ckj, sol_crq,
                 sol_crk, sol_ckr, sol_cnr, sol_crn, sol_cnk, sol_ckn,
		 sol_crj, sol_cjr, sol_cqk, sol_ckq,
                 sol_cqj, sol_cjq, sol_cnq, sol_cqn, sol_cnj, sol_cjn,
                 //sol_cP_rk, sol_cM_rk, sol_cP_nr, sol_cM_nr, sol_cP_nk, sol_cM_nk,
                 sol_ll_dEta, sol_ll_dPhi, sol_ll_dR, sol_cHel,
                 sol_kNorm, sol_rNorm, sol_nNorm, sol_jNorm, sol_qNorm,
                 sol_phi0, sol_phi1);

    compEvtShp(sol_pLep, sol_aLep, sol_pNeu, sol_aNeu, sol_pBot, sol_aBot,
               sol_TT_sph, sol_TT_apl, sol_TT_cir, sol_TT_iso, sol_TT_C, sol_TT_D, sol_TT_spt,
               sol_TT_H0, sol_TT_H1, sol_TT_H2, sol_TT_H3, sol_TT_H4,
               sol_TT_R1, sol_TT_R2, sol_TT_R3, sol_TT_R4);
  }
   
  if ( topGenObjects.valuesSet_ ) {

    isTopGen.value_ = 1;
    topProdMode.value_ = topGenObjects.productionMode_; // 0 gg 1 qq 2 qg
    genWeight.value_ = ifOk( genLevelWeights.trueLevelWeight_ );

    gen_pLep = common::LVtoTLV(*topGenObjects.GenLepton_);
    gen_aLep = common::LVtoTLV(*topGenObjects.GenAntiLepton_);

    gen_pNeu = common::LVtoTLV(*topGenObjects.GenNeutrino_);
    gen_aNeu = common::LVtoTLV(*topGenObjects.GenAntiNeutrino_);

    // -1 - no b-jet found, -2 - two b-jets from one top in generator
    // use genJets if found, otherwise use quarks
    if (topGenObjects.GenTopBHadIndex_ > -1 and topGenObjects.GenAntiTopBHadIndex_ > -1) {
      gen_pBot = common::LVtoTLV( (*topGenObjects.allGenJets_).at(topGenObjects.GenTopBHadIndex_) );
      gen_aBot = common::LVtoTLV( (*topGenObjects.allGenJets_).at(topGenObjects.GenAntiTopBHadIndex_) );
    }
    else {
      gen_pBot = common::LVtoTLV(*topGenObjects.GenB_);
      gen_aBot = common::LVtoTLV(*topGenObjects.GenAntiB_);
    }

    gen_pTop = common::LVtoTLV(*topGenObjects.GenTop_);
    gen_aTop = common::LVtoTLV(*topGenObjects.GenAntiTop_);

    gen_TT = gen_pTop + gen_aTop;

    gen_pLep_pt.value_ = ifOk( gen_pLep.Pt() );
    gen_pLep_eta.value_ = ifOk( gen_pLep.Eta() );
    gen_pLep_phi.value_ = ifOk( gen_pLep.Phi() );
    gen_aLep_pt.value_ = ifOk( gen_aLep.Pt() );
    gen_aLep_eta.value_ = ifOk( gen_aLep.Eta() );
    gen_aLep_phi.value_ = ifOk( gen_aLep.Phi() );

    gen_pNeu_pt.value_ = ifOk( gen_pNeu.Pt() );
    gen_pNeu_eta.value_ = ifOk( gen_pNeu.Eta() );
    gen_pNeu_phi.value_ = ifOk( gen_pNeu.Phi() );
    gen_aNeu_pt.value_ = ifOk( gen_aNeu.Pt() );
    gen_aNeu_eta.value_ = ifOk( gen_aNeu.Eta() );
    gen_aNeu_phi.value_ = ifOk( gen_aNeu.Phi() );

    gen_pBot_pt.value_ = ifOk( gen_pBot.Pt() );
    gen_pBot_eta.value_ = ifOk( gen_pBot.Eta() );
    gen_pBot_phi.value_ = ifOk( gen_pBot.Phi() );
    gen_aBot_pt.value_ = ifOk( gen_aBot.Pt() );
    gen_aBot_eta.value_ = ifOk( gen_aBot.Eta() );
    gen_aBot_phi.value_ = ifOk( gen_aBot.Phi() );

    gen_pTop_pt.value_ = ifOk( gen_pTop.Pt() );
    gen_pTop_rap.value_ = ifOk( gen_pTop.Rapidity() );
    gen_pTop_phi.value_ = ifOk( gen_pTop.Phi() );
    gen_aTop_pt.value_ = ifOk( gen_aTop.Pt() );
    gen_aTop_rap.value_ = ifOk( gen_aTop.Rapidity() );
    gen_aTop_phi.value_ = ifOk( gen_aTop.Phi() );

    gen_TT_m.value_ = ifOk( gen_TT.M() );
    gen_TT_pt.value_ = ifOk( gen_TT.Pt() );
    gen_TT_rap.value_ = ifOk( gen_TT.Rapidity() );
    gen_TT_phi.value_ = ifOk( gen_TT.Phi() );

    // SpinCorr variables without boosting
    gen_LL_dEta.value_ = ifOk( gen_aLep.Eta() - gen_pLep.Eta() );
    gen_LL_dPhi.value_ = ifOk( gen_aLep.DeltaPhi( gen_pLep ) );
    gen_LL_dR.value_ = ifOk( gen_aLep.DeltaR( gen_pLep ) );
    gen_cLab.value_ = ifOk( gen_aLep.Vect().Unit().Dot( gen_pLep.Vect().Unit() ) );

    compSpinCorr(gen_pLep, gen_aLep, gen_pTop, gen_aTop,
                 gen_b1k, gen_b1j, gen_b1r, gen_b1q, gen_b1n,
                 gen_b2k, gen_b2j, gen_b2r, gen_b2q, gen_b2n,
                 //gen_bP_kk, gen_bM_kk, gen_bP_jj, gen_bM_jj,
                 //gen_bP_rr, gen_bM_rr, gen_bP_qq, gen_bM_qq,
                 //gen_bP_nn, gen_bM_nn,
                 gen_ckk, gen_crr, gen_cnn,
                 gen_ckj, gen_crq,
                 gen_crk, gen_ckr, gen_cnr, gen_crn, gen_cnk, gen_ckn,
		 gen_crj, gen_cjr, gen_cqk, gen_ckq,
                 gen_cqj, gen_cjq, gen_cnq, gen_cqn, gen_cnj, gen_cjn,
                 //gen_cP_rk, gen_cM_rk, gen_cP_nr, gen_cM_nr, gen_cP_nk, gen_cM_nk,
                 gen_ll_dEta, gen_ll_dPhi, gen_ll_dR, gen_cHel,
                 gen_kNorm, gen_rNorm, gen_nNorm, gen_jNorm, gen_qNorm,
                 gen_phi0, gen_phi1);

    compEvtShp(gen_pLep, gen_aLep, gen_pNeu, gen_aNeu, gen_pBot, gen_aBot,
               gen_TT_sph, gen_TT_apl, gen_TT_cir, gen_TT_iso, gen_TT_C, gen_TT_D, gen_TT_spt,
               gen_TT_H0, gen_TT_H1, gen_TT_H2, gen_TT_H3, gen_TT_H4,
               gen_TT_R1, gen_TT_R2, gen_TT_R3, gen_TT_R4);
   }

  if (isKinReco.value_ + isTopGen.value_ == 2) {
    const double dRMax = 0.3; // past matching 0.1 too tight?
    pBotMatchSolGen.value_ = (sol_pBot.DeltaR(gen_pBot) < dRMax) ? 1 : 0;
    aBotMatchSolGen.value_ = (sol_aBot.DeltaR(gen_aBot) < dRMax) ? 1 : 0;
  }
}



double VariablesPhiTT::ifOk(const double& vCheck) const
{
  double vOk = std::isfinite(vCheck) ? vCheck : -999999.;

  if(vOk == -999999) std::cout << "  Error due to non finite problem "  <<  std::endl;

  return vOk;
}



void VariablesPhiTT::compSpinCorr(const TLorentzVector& p4_pLep, const TLorentzVector& p4_aLep,
                                  const TLorentzVector& p4_pTop, const TLorentzVector& p4_aTop,
                                  VariableFloat& b1k, VariableFloat& b1j, VariableFloat& b1r, VariableFloat& b1q, VariableFloat& b1n,
                                  VariableFloat& b2k, VariableFloat& b2j, VariableFloat& b2r, VariableFloat& b2q, VariableFloat& b2n,
                                  //VariableFloat& bP_kk, VariableFloat& bM_kk,
				  //VariableFloat& bP_jj, VariableFloat& bM_jj,
                                  //VariableFloat& bP_rr, VariableFloat& bM_rr,
                                  //VariableFloat& bP_qq, VariableFloat& bM_qq,
                                  //VariableFloat& bP_nn, VariableFloat& bM_nn,
                                  VariableFloat& ckk, VariableFloat& crr, VariableFloat& cnn,
                                  VariableFloat& ckj, VariableFloat& crq,
                                  VariableFloat& c_rk, VariableFloat& c_kr,
                                  VariableFloat& c_nr, VariableFloat& c_rn,
                                  VariableFloat& c_nk, VariableFloat& c_kn,

                                  VariableFloat& c_rj, VariableFloat& c_jr,
                                  VariableFloat& c_qk, VariableFloat& c_kq,

                                  VariableFloat& c_qj, VariableFloat& c_jq,
                                  VariableFloat& c_nq, VariableFloat& c_qn,
                                  VariableFloat& c_nj, VariableFloat& c_jn,
                                  //VariableFloat& cP_rk, VariableFloat& cM_rk,
                                  //VariableFloat& cP_nr, VariableFloat& cM_nr,
                                  //VariableFloat& cP_nk, VariableFloat& cM_nk,
                                  VariableFloat& ll_dEta, VariableFloat& ll_dPhi, VariableFloat& ll_dR, VariableFloat& cHel,
                                  VariableFloat& kNorm, VariableFloat& rNorm, VariableFloat& nNorm, VariableFloat& jNorm, VariableFloat& qNorm,
                                  VariableFloat& phi0, VariableFloat& phi1)
{
  // Note: convention: capital for lab, small for hel -> LL_dPhi = dPhiLab, ll_dPhi = dPhiHel
  // The various Bernreuther bases
  TVector3 kBase, jBase, rBase, qBase, nBase;

  // Some proxy quantities not saved in
  double c_pTP, s_pTP, sY, sD;
  double ck_aLep, ck_pLep, cr_aLep, cr_pLep, cn_aLep, cn_pLep;
  double cj_aLep, cj_pLep, cq_aLep, cq_pLep;
  double crk, ckr, cnr, crn, cnk, ckn;
  double crj, cjr, cqk, ckq;
  double cqj, cjq, cnq, cqn, cnj, cjn;

  // Beam unit vector
  TVector3 p3_pPro(0., 0., 1.);

  // The bases definition: Bernreuther spinMatrix 1508.05271
  TLorentzVector p4_TT( p4_pTop + p4_aTop );
  TVector3 b4_TT( -1. * p4_TT.BoostVector() );

  TLorentzVector b4_pTop( p4_pTop );
  b4_pTop.Boost( b4_TT );

  TLorentzVector b4_aTop( p4_aTop );
  b4_aTop.Boost( b4_TT );

  TLorentzVector b4_aLep( p4_aLep );
  b4_aLep.Boost( b4_TT );
  b4_aLep.Boost( -1. * b4_pTop.BoostVector() );

  TLorentzVector b4_pLep( p4_pLep );
  b4_pLep.Boost( b4_TT );
  b4_pLep.Boost( -1. * b4_aTop.BoostVector() );

  // Calculating the top-beam angle for pTop only
  c_pTP = b4_pTop.Vect().Unit().Dot(p3_pPro);
  s_pTP = std::sqrt(1. - (c_pTP * c_pTP));

  // The signs needed to account for Bose symmetry
  sY = (c_pTP >= 0.) ? 1. : -1.;
  sD = ( std::abs(p4_pTop.Rapidity()) >= std::abs(p4_aTop.Rapidity()) ) ? 1. : -1.;

  // Define the base vectors a
  // j and q base are the k* and r* respectively
  // b is always -a
  kBase = b4_pTop.Vect().Unit();
  jBase = sD * kBase;
  rBase = ( (sY / s_pTP) * (p3_pPro - (c_pTP * kBase)) ).Unit();
  qBase = sD * rBase;
  nBase = ( (sY / s_pTP) * p3_pPro.Cross(kBase) ).Unit();

  // Find the relevant angles in these bases
  ck_aLep = b4_aLep.Vect().Unit().Dot( kBase );
  ck_pLep = b4_pLep.Vect().Unit().Dot( -1. * kBase );

  cj_aLep = b4_aLep.Vect().Unit().Dot( jBase );
  cj_pLep = b4_pLep.Vect().Unit().Dot( -1. * jBase );

  cr_aLep = b4_aLep.Vect().Unit().Dot( rBase );
  cr_pLep = b4_pLep.Vect().Unit().Dot( -1. * rBase );

  cq_aLep = b4_aLep.Vect().Unit().Dot( qBase );
  cq_pLep = b4_pLep.Vect().Unit().Dot( -1. * qBase );

  cn_aLep = b4_aLep.Vect().Unit().Dot( nBase );
  cn_pLep = b4_pLep.Vect().Unit().Dot( -1. * nBase );


  // Fill the raw angles into VarFloats
  b1k.value_ = ifOk( ck_aLep );
  b2k.value_ = ifOk( ck_pLep );

  b1j.value_ = ifOk( cj_aLep );
  b2j.value_ = ifOk( cj_pLep );

  b1r.value_ = ifOk( cr_aLep );
  b2r.value_ = ifOk( cr_pLep );

  b1q.value_ = ifOk( cq_aLep );
  b2q.value_ = ifOk( cq_pLep );

  b1n.value_ = ifOk( cn_aLep );
  b2n.value_ = ifOk( cn_pLep );

  // Now we can squeeze it all out based on table 5 page 16
  // The B1 ~ c_aLep, B2 ~ c_pLep sums
  //bP_kk.value_ = ifOk( ck_aLep + ck_pLep );
  //bM_kk.value_ = ifOk( ck_aLep - ck_pLep );

  //bP_jj.value_ = ifOk( cj_aLep + cj_pLep );
  //bM_jj.value_ = ifOk( cj_aLep - cj_pLep );

  //bP_rr.value_ = ifOk( cr_aLep + cr_pLep );
  //bM_rr.value_ = ifOk( cr_aLep - cr_pLep );

  //bP_qq.value_ = ifOk( cq_aLep + cq_pLep );
  //bM_qq.value_ = ifOk( cq_aLep - cq_pLep );

  //bP_nn.value_ = ifOk( cn_aLep + cn_pLep );
  //bM_nn.value_ = ifOk( cn_aLep - cn_pLep );

  // spinCorr coeff Cab = -9<cab>
  ckk.value_ = ifOk( ck_aLep * ck_pLep );
  crr.value_ = ifOk( cr_aLep * cr_pLep );
  cnn.value_ = ifOk( cn_aLep * cn_pLep );
  ckj.value_ = ifOk( ck_aLep * cj_pLep );
  crq.value_ = ifOk( cr_aLep * cq_pLep );

  crk = cr_aLep * ck_pLep;
  ckr = ck_aLep * cr_pLep;

  cnr = cn_aLep * cr_pLep;
  crn = cr_aLep * cn_pLep;

  cnk = cn_aLep * ck_pLep;
  ckn = ck_aLep * cn_pLep;

  crj = cr_aLep * cj_pLep;
  cjr = cj_aLep * cr_pLep;

  cqk = cq_aLep * ck_pLep;
  ckq = ck_aLep * cq_pLep;

  cqj = cq_aLep * cj_pLep;
  cjq = cj_aLep * cq_pLep;

  cnq = cn_aLep * cq_pLep;
  cqn = cq_aLep * cn_pLep;
  
  cnj = cn_aLep * cj_pLep;
  cjn = cj_aLep * cn_pLep;

  //  if(ckk.value_ == -999999) p4_pLep.Print();
  //  if(ckk.value_ == -999999) p4_aLep.Print();
  //  if(ckk.value_ == -999999) p4_pTop.Print();
  //  if(ckk.value_ == -999999) p4_aTop.Print();

  //+++++++++++++++++++++++++++++++//
  //ajeeta -- adding the cos(theta_i)*cos(theta_j) also not just their sums and differences
  c_rk.value_ = ifOk(crk);
  c_kr.value_ = ifOk(ckr);
  c_nr.value_ = ifOk(cnr);
  c_rn.value_ = ifOk(crn);
  c_nk.value_ = ifOk(cnk);
  c_kn.value_ = ifOk(ckn);

  c_rj.value_ = ifOk(crj);
  c_jr.value_ = ifOk(cjr);
  c_qk.value_ = ifOk(cqk);
  c_kq.value_ = ifOk(ckq);

  c_qj.value_ = ifOk(cqj);
  c_jq.value_ = ifOk(cjq);
  c_nq.value_ = ifOk(cnq);
  c_qn.value_ = ifOk(cqn);
  c_nj.value_ = ifOk(cnj);
  c_jn.value_ = ifOk(cjn);
  //+++++++++++++++++++++++++++++++//

  //cP_rk.value_ = ifOk( crk + ckr );
  //cM_rk.value_ = ifOk( crk - ckr );

  //cP_nr.value_ = ifOk( cnr + crn );
  //cM_nr.value_ = ifOk( cnr - crn );

  //cP_nk.value_ = ifOk( cnk + ckn );
  //cM_nk.value_ = ifOk( cnk - ckn );

  // Find also the opening angles of the lepton
  ll_dEta.value_ = ifOk( b4_aLep.Eta() - b4_pLep.Eta() );

  ll_dPhi.value_ = ifOk( b4_aLep.DeltaPhi( b4_pLep ) );
  ll_dR.value_ = ifOk( b4_aLep.DeltaR( b4_pLep ) );
  cHel.value_ = ifOk( b4_aLep.Vect().Unit().Dot( b4_pLep.Vect().Unit() ) );

  // These are the O_CP1 and O_CP2 as in page 18 (why not O_CP3 with n base too)
  TVector3 llNorm = b4_aLep.Vect().Unit().Cross( b4_pLep.Vect().Unit() );
  kNorm.value_ =  ifOk( llNorm.Dot( kBase ) );
  rNorm.value_ =  ifOk( llNorm.Dot( rBase ) );
  nNorm.value_ =  ifOk( llNorm.Dot( nBase ) );
  jNorm.value_ =  ifOk( llNorm.Dot( jBase ) );
  qNorm.value_ =  ifOk( llNorm.Dot( qBase ) );


  // Weirder angles as in 1702.06063; phi0 = phi* and phi1 = phi*_CP
  TVector3 t3_aLep = ( b4_aLep.Vect().Unit() - (b4_aLep.Vect().Unit().Dot( kBase ) * kBase) ).Unit();
  TVector3 t3_pLep = ( b4_pLep.Vect().Unit() - (b4_pLep.Vect().Unit().Dot( kBase ) * kBase) ).Unit();

  phi0.value_ = ifOk( std::acos(t3_aLep.Dot(t3_pLep)) );
  phi1.value_ = ifOk( (kNorm.value_ < 0.) ? (2. * M_PI) - phi0.value_ : phi0.value_ );
}



void VariablesPhiTT::compEvtShp(const TLorentzVector& p4_pLep, const TLorentzVector& p4_aLep,
                                const TLorentzVector& p4_pNeu, const TLorentzVector& p4_aNeu,
                                const TLorentzVector& p4_pBot, const TLorentzVector& p4_aBot,
                                VariableFloat& TT_sph, VariableFloat& TT_apl, VariableFloat& TT_cir, VariableFloat& TT_iso,
                                VariableFloat& TT_C, VariableFloat& TT_D, VariableFloat& TT_spt,
                                VariableFloat& TT_H0, VariableFloat& TT_H1, VariableFloat& TT_H2,
                                VariableFloat& TT_H3, VariableFloat& TT_H4,
                                VariableFloat& TT_R1, VariableFloat& TT_R2, VariableFloat& TT_R3, VariableFloat& TT_R4)
{
  // Calculate tt event shape considering W as jet
  // W labeling follows the charged lepton
  TLorentzVector p4_pW( p4_pLep + p4_aNeu );
  TLorentzVector p4_aW( p4_aLep + p4_pNeu );

  // Not sure if ordering matters, but pTop then aTop
  std::vector<TLorentzVector> vDaus;
  vDaus.push_back(p4_aW);
  vDaus.push_back(p4_pBot);
  vDaus.push_back(p4_pW);
  vDaus.push_back(p4_aBot);
  const EventShapeVariables TT_evtShp(vDaus);

  TT_sph.value_ = ifOk( TT_evtShp.sphericity() );
  TT_apl.value_ = ifOk( TT_evtShp.aplanarity() );
  TT_cir.value_ = ifOk( TT_evtShp.circularity() );
  TT_iso.value_ = ifOk( TT_evtShp.isotropy() );
  TT_C.value_ = ifOk( TT_evtShp.C() );
  TT_D.value_ = ifOk( TT_evtShp.D() );
  TT_spt.value_ = ifOk( TT_evtShp.transSphericity() );

  TT_H0.value_ = ifOk( TT_evtShp.H(0) );
  TT_H1.value_ = ifOk( TT_evtShp.H(1) );
  TT_H2.value_ = ifOk( TT_evtShp.H(2) );
  TT_H3.value_ = ifOk( TT_evtShp.H(3) );
  TT_H4.value_ = ifOk( TT_evtShp.H(4) );

  TT_R1.value_ = ifOk( TT_evtShp.R(1) );
  TT_R2.value_ = ifOk( TT_evtShp.R(2) );
  TT_R3.value_ = ifOk( TT_evtShp.R(3) );
  TT_R4.value_ = ifOk( TT_evtShp.R(4) );
}



VariablesPhiTT* VariablesPhiTT::fillVariables(const EventMetadata& eventMetadata, 
                                              const RecoObjects& recoObjects,
                                              const CommonGenObjects& commonGenObjects,
                                              const TopGenObjects& topGenObjects,
                                              const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                                              const ttbar::RecoObjectIndices& recoObjectIndices,
                                              const ttbar::GenObjectIndices& genObjectIndices,
                                              const ttbar::GenLevelWeights& genLevelWeights,
                                              const ttbar::RecoLevelWeights& recoLevelWeights,
                                              const double& weight)
{
    return new VariablesPhiTT(eventMetadata, recoObjects, commonGenObjects, topGenObjects, kinematicReconstructionSolutions, recoObjectIndices, genObjectIndices, genLevelWeights, recoLevelWeights, weight);
}



// ---------------------------------- Class VariablesPhiTT::EventShapeVariables -------------------------------------------

VariablesPhiTT::EventShapeVariables::EventShapeVariables(const std::vector<TLorentzVector>& inputVectors):
inputVectors_(this->transformInputVectors(inputVectors))
{}



std::vector<ROOT::Math::XYZVector> VariablesPhiTT::EventShapeVariables::transformInputVectors(const std::vector<TLorentzVector>& jets)const
{
    std::vector<ROOT::Math::XYZVector> result;
    for(size_t i = 0; i < jets.size(); ++i){
        const TLorentzVector& jetP4(jets.at(i));
        const ROOT::Math::XYZVector jetP3 = ROOT::Math::XYZVector(jetP4.Px(), jetP4.Py(), jetP4.Pz());
        result.push_back(jetP3);
    }
    return result;
}



double VariablesPhiTT::EventShapeVariables::isotropy(unsigned int numberOfSteps)const
{
    double eIn(-1.);
    double eOut(-1.);
    
    const double deltaPhi = 2*TMath::Pi()/numberOfSteps;
    double phi(0.);
    for(unsigned int i = 0; i < numberOfSteps; ++i){
        phi += deltaPhi;
        
        // Sum over inner product of unit vectors and momenta
        double sum(0.);
        for(unsigned int j = 0; j < inputVectors_.size(); ++j)
            sum += TMath::Abs(TMath::Cos(phi)*inputVectors_.at(j).x()+TMath::Sin(phi)*inputVectors_.at(j).y());
        
        if(eOut<0. || sum<eOut) eOut = sum;
        if(eIn<0. || sum>eIn) eIn = sum;
    }
    
    return (eIn-eOut)/eIn;
}



double VariablesPhiTT::EventShapeVariables::circularity(unsigned int numberOfSteps)const
{
    double circularity(-1.);
    
    double area(0.);
    for(unsigned int i = 0; i < inputVectors_.size(); ++i)
        area += TMath::Sqrt(inputVectors_.at(i).x()*inputVectors_.at(i).x()+inputVectors_.at(i).y()*inputVectors_.at(i).y());
    
    const double deltaPhi = 2*TMath::Pi()/numberOfSteps;
    double phi(0.);
    for(unsigned int i = 0; i < numberOfSteps; ++i){
        phi += deltaPhi;
        
        double sum(0.);
        for(unsigned int j = 0; j < inputVectors_.size(); ++j)
            sum += TMath::Abs(TMath::Cos(phi)*inputVectors_.at(j).x() + TMath::Sin(phi)*inputVectors_.at(j).y());
        
        const double tmpCircularity = TMath::Pi()/2.*sum/area;
        if(circularity<0. || tmpCircularity<circularity) circularity = tmpCircularity;
    }
    
    return circularity;
}



double VariablesPhiTT::EventShapeVariables::sphericity(const double& r)const
{
    const TVectorD eigenValues = this->compEigenValues(r);
    return 1.5*(eigenValues(1) + eigenValues(2));
}



double VariablesPhiTT::EventShapeVariables::aplanarity(const double& r)const
{
    const TVectorD eigenValues = this->compEigenValues(r);
    return 1.5*eigenValues(2);
}



double VariablesPhiTT::EventShapeVariables::C(const double& r)const
{
    const TVectorD eigenValues = this->compEigenValues(r);
    return 3.*(eigenValues(0)*eigenValues(1) + eigenValues(0)*eigenValues(2) + eigenValues(1)*eigenValues(2));
}



double VariablesPhiTT::EventShapeVariables::D(const double& r)const
{
    const TVectorD eigenValues = this->compEigenValues(r);
    return 27.*eigenValues(0)*eigenValues(1)*eigenValues(2);
}



double VariablesPhiTT::EventShapeVariables::transSphericity(const double& r)const
{
    const TVectorD eigenValues = this->compEigenValues(r);
    return 2.*eigenValues(1)/(eigenValues(0) + eigenValues(1));
}



double VariablesPhiTT::EventShapeVariables::H(const int order)const
{
    const FoxWolframMoments fwm(inputVectors_, order);
    return fwm.H(order);
}



double VariablesPhiTT::EventShapeVariables::R(const int order)const
{
    const FoxWolframMoments fwm(inputVectors_, order);
    return fwm.R(order);
}



TVectorD VariablesPhiTT::EventShapeVariables::compEigenValues(const double& r)const
{
    TVectorD eigenValues(3);
    const TMatrixDSym tensor = this->compMomentumTensor(r);
    if(tensor.IsSymmetric() && tensor.NonZeros()!=0) tensor.EigenVectors(eigenValues);
    return eigenValues;
}



TMatrixDSym VariablesPhiTT::EventShapeVariables::compMomentumTensor(const double& r)const
{
    TMatrixDSym momentumTensor(3);
    momentumTensor.Zero();
    if(inputVectors_.size() < 2) return momentumTensor;
    
    // Fill momentumTensor from inputVectors
    double norm(1.);
    for(int i = 0; i < (int)inputVectors_.size(); ++i){
        const ROOT::Math::XYZVector& inputVector(inputVectors_.at(i));
        const double p2 = inputVector.Dot(inputVector);
        const double pR = r==2. ? p2 : TMath::Power(p2, 0.5*r);
        norm += pR;
        const double pRminus2 = r==2. ? 1. : TMath::Power(p2, 0.5*r - 1.);
        momentumTensor(0,0) += pRminus2*inputVector.x()*inputVector.x();
        momentumTensor(0,1) += pRminus2*inputVector.x()*inputVector.y();
        momentumTensor(0,2) += pRminus2*inputVector.x()*inputVector.z();
        momentumTensor(1,0) += pRminus2*inputVector.y()*inputVector.x();
        momentumTensor(1,1) += pRminus2*inputVector.y()*inputVector.y();
        momentumTensor(1,2) += pRminus2*inputVector.y()*inputVector.z();
        momentumTensor(2,0) += pRminus2*inputVector.z()*inputVector.x();
        momentumTensor(2,1) += pRminus2*inputVector.z()*inputVector.y();
        momentumTensor(2,2) += pRminus2*inputVector.z()*inputVector.z();
    }
    
    return (1./norm)*momentumTensor;
}


// ---------------------------------- Class VariablesPhiTT::FoxWolframMoments -------------------------------------------

VariablesPhiTT::FoxWolframMoments::FoxWolframMoments(const std::vector<ROOT::Math::XYZVector>& inputVectors,
                                                                           const int maxorder):
nMoment_(maxorder+1),
fwArray_(nMoment_),
sumArray_(nMoment_)
{
    for(int i = 0; i < nMoment_; ++i){ 
        fwArray_(i) = 0.;
        sumArray_(i) = 0.;
    }
    
    this->compute(inputVectors);
}



void VariablesPhiTT::FoxWolframMoments::compute(const std::vector<ROOT::Math::XYZVector>& inputVectors)
{
    if(inputVectors.size() == 0) return;
    
    // Loop over the all particle candidates
    double s(0.);
    for(size_t i = 0; i<inputVectors.size(); ++i){
        // Candidate particle's 3-momentum
        const TVector3 p1(inputVectors.at(i).X(), inputVectors.at(i).Y(), inputVectors.at(i).Z());
        const double pmag1 = p1.Mag();

        // Loop over other particle's candidates, starting at the next one in the list
        for(size_t j = i; j < inputVectors.size(); ++j){
            // Candidate particle's 3-momentum
            const TVector3 p2(inputVectors.at(j).X(),inputVectors.at(j).Y(),inputVectors.at(j).Z());
            
            // Cosine of the angle between the two candidate particles
            const double cosPhi = TMath::Cos(p1.Angle(p2));
            
            // Contribution of this pair of track
            // (note the factor 2 : the pair enters the sum twice)
            const double pmag2 = p2.Mag();
            for(int l = 0; l < nMoment_; ++l) sumArray_(l) += 2 * pmag1 * pmag2 * this->legendre(l, 0, cosPhi);
        }
        
        // Contribution for this moment
        for(int l = 0; l < nMoment_; ++l) sumArray_(l) += pmag1*pmag1*this->legendre(l, 0, 1.);
        
        // Total energy
        s += p1.Mag();
    }

    if(s <= 0.) return;

    // Normalize Fox Wolfram Moments
    for(int i = 0; i < nMoment_; ++i) fwArray_(i) = sumArray_(i)/TMath::Power(s, 2);
}



double VariablesPhiTT::FoxWolframMoments::R(const int order)const
{
    if(this->H(0)>0. && order<nMoment_) return (this->H(order)/this->H(0));
    return 0.;
}



double VariablesPhiTT::FoxWolframMoments::legendre(const int l, const int m, const double& x)const
{
    assert(m >= 0.);
    assert(m <= l);
    assert(std::abs(x) <= 1.);
    
    double pmm(1.);
    if(m > 0){
        const double somx2 = TMath::Sqrt((1. - x) * (1. + x));
        double fact = 1.;
        for(int i = 0; i < m; ++i){
            pmm *= -fact * somx2;
            fact += 2.;
        }
    }
    if(l == m) return pmm;

    double pmmp1 = x * (2 * m + 1) * pmm;
    if(l == m + 1) return pmmp1;
    
    for(int ll = m+2; ll <= l; ++ll){
        const double pll = (x*(2 * ll - 1)*pmmp1 - (ll + m - 1)*pmm)/(ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pmmp1;
}


