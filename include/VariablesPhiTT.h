#ifndef VariablesPhiTT_h
#define VariablesPhiTT_h

#include <map>
#include <vector>
#include <string>

#include "../../common/include/classesFwd.h"

#include "VariablesBase.h"
#include "analysisStructsFwd.h"

#include <TVector.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <Math/Vector3D.h>
#include <TLorentzVector.h>
#include <Math/VectorUtil.h>

#include "TMatrixDSym.h"

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


/// Class holding the variables for A/H analysis
class VariablesPhiTT : public VariablesBase{
    
public:
    
    /// Empty constructor
    VariablesPhiTT();
    
    /// Constructor setting up input variables from physics objects
    VariablesPhiTT(const EventMetadata& eventMetadata, 
                   const RecoObjects& recoObjects,
                   const CommonGenObjects& commonGenObjects,
                   const TopGenObjects& topGenObjects,
                   const KinematicReconstructionSolutions& kinRecoObjects,
                   const ttbar::RecoObjectIndices& recoObjectIndices,
                   const ttbar::GenObjectIndices& genObjectIndices,
                   const ttbar::GenLevelWeights& genLevelWeights,
                   const ttbar::RecoLevelWeights& recoLevelWeights,
                   const double& weight);
    
    /// Destructor
    ~VariablesPhiTT(){}
    
    /// Fill the input structs for all jet combinations of one event
    static VariablesPhiTT* fillVariables(const EventMetadata& eventMetadata,
                                         const RecoObjects& recoObjects,
                                         const CommonGenObjects& commonGenObjects,
                                         const TopGenObjects& topGenObjects,
                                         const KinematicReconstructionSolutions& kinRecoObjects,
                                         const ttbar::RecoObjectIndices& recoObjectIndices,
                                         const ttbar::GenObjectIndices& genObjectIndices,
                                         const ttbar::GenLevelWeights& genLevelWeights,
                                         const ttbar::RecoLevelWeights& recoLevelWeights,
                                         const double& weight);

    // Solution variables - p particle a anti
    VariableFloat sol_pLep_pt, sol_pLep_eta, sol_pLep_phi;
    VariableFloat sol_aLep_pt, sol_aLep_eta, sol_aLep_phi;

    VariableFloat sol_pNeu_pt, sol_pNeu_eta, sol_pNeu_phi;
    VariableFloat sol_aNeu_pt, sol_aNeu_eta, sol_aNeu_phi;

    VariableFloat sol_pBot_pt, sol_pBot_eta, sol_pBot_phi;
    VariableFloat sol_aBot_pt, sol_aBot_eta, sol_aBot_phi;

    VariableFloat sol_pTop_pt, sol_pTop_rap, sol_pTop_phi;
    VariableFloat sol_aTop_pt, sol_aTop_rap, sol_aTop_phi;

    VariableFloat sol_TT_m, sol_TT_pt, sol_TT_rap, sol_TT_phi;

    VariableFloat sol_LL_dEta, sol_LL_dPhi, sol_LL_dR, sol_cLab;

    // Spin correlation stuff
    VariableFloat sol_b1k, sol_b1j, sol_b1r, sol_b1q, sol_b1n;
    VariableFloat sol_b2k, sol_b2j, sol_b2r, sol_b2q, sol_b2n;

    //VariableFloat sol_bP_kk, sol_bM_kk, sol_bP_jj, sol_bM_jj, sol_bP_rr;
    //VariableFloat sol_bM_rr, sol_bP_qq, sol_bM_qq, sol_bP_nn, sol_bM_nn;

    VariableFloat sol_ckk, sol_crr, sol_cnn;
    VariableFloat sol_ckj, sol_crq;
    //ajeeta
    VariableFloat sol_crk, sol_ckr, sol_cnr, sol_crn, sol_cnk, sol_ckn;
    VariableFloat sol_cqj, sol_cjq, sol_cnq, sol_cqn, sol_cnj, sol_cjn;
    VariableFloat sol_cqk, sol_ckq, sol_crj, sol_cjr;

    //VariableFloat sol_cP_rk, sol_cM_rk;
    //VariableFloat sol_cP_nr, sol_cM_nr;
    //VariableFloat sol_cP_nk, sol_cM_nk;

    VariableFloat sol_ll_dEta, sol_ll_dPhi, sol_ll_dR, sol_cHel;
    VariableFloat sol_kNorm, sol_rNorm, sol_nNorm, sol_jNorm, sol_qNorm;
    VariableFloat sol_phi0, sol_phi1;

    // Event shape stuff with all WbWb - the W's each considered as a single "jet"
    VariableFloat sol_TT_sph, sol_TT_apl, sol_TT_cir, sol_TT_iso;
    VariableFloat sol_TT_C, sol_TT_D, sol_TT_spt;
    VariableFloat sol_TT_H0, sol_TT_H1, sol_TT_H2, sol_TT_H3, sol_TT_H4;
    VariableFloat sol_TT_R1, sol_TT_R2, sol_TT_R3, sol_TT_R4;

    // Gen variables
    VariableFloat gen_pLep_pt, gen_pLep_eta, gen_pLep_phi;
    VariableFloat gen_aLep_pt, gen_aLep_eta, gen_aLep_phi;

    VariableFloat gen_pNeu_pt, gen_pNeu_eta, gen_pNeu_phi;
    VariableFloat gen_aNeu_pt, gen_aNeu_eta, gen_aNeu_phi;

    VariableFloat gen_pBot_pt, gen_pBot_eta, gen_pBot_phi;
    VariableFloat gen_aBot_pt, gen_aBot_eta, gen_aBot_phi;

    VariableFloat gen_pTop_pt, gen_pTop_rap, gen_pTop_phi;
    VariableFloat gen_aTop_pt, gen_aTop_rap, gen_aTop_phi;

    VariableFloat gen_TT_m, gen_TT_pt, gen_TT_rap, gen_TT_phi;

    VariableFloat gen_LL_dEta, gen_LL_dPhi, gen_LL_dR, gen_cLab;

    // Spin correlation stuff
    VariableFloat gen_b1k, gen_b1j, gen_b1r, gen_b1q, gen_b1n;
    VariableFloat gen_b2k, gen_b2j, gen_b2r, gen_b2q, gen_b2n;

    //VariableFloat gen_bP_kk, gen_bM_kk, gen_bP_jj, gen_bM_jj, gen_bP_rr;
    //VariableFloat gen_bM_rr, gen_bP_qq, gen_bM_qq, gen_bP_nn, gen_bM_nn;

    VariableFloat gen_ckk, gen_crr, gen_cnn;
    VariableFloat gen_ckj, gen_crq;

    //ajeeta
    VariableFloat gen_crk, gen_ckr, gen_cnr, gen_crn, gen_cnk, gen_ckn;
    VariableFloat gen_cqj, gen_cjq, gen_cnq, gen_cqn, gen_cnj, gen_cjn;
    VariableFloat gen_cqk, gen_ckq, gen_crj, gen_cjr;

    //VariableFloat gen_cP_rk, gen_cM_rk;
    //VariableFloat gen_cP_nr, gen_cM_nr;
    //VariableFloat gen_cP_nk, gen_cM_nk;

    VariableFloat gen_ll_dEta, gen_ll_dPhi, gen_ll_dR, gen_cHel;
    VariableFloat gen_kNorm, gen_rNorm, gen_nNorm, gen_jNorm, gen_qNorm;
    VariableFloat gen_phi0, gen_phi1;

    // Event shape stuff with all WbWb - the W's each considered as a single "jet"
    VariableFloat gen_TT_sph, gen_TT_apl, gen_TT_cir, gen_TT_iso;
    VariableFloat gen_TT_C, gen_TT_D, gen_TT_spt;
    VariableFloat gen_TT_H0, gen_TT_H1, gen_TT_H2, gen_TT_H3, gen_TT_H4;
    VariableFloat gen_TT_R1, gen_TT_R2, gen_TT_R3, gen_TT_R4;

    // b-jet matching between sol-gen
    VariableInt pBotMatchSolGen, aBotMatchSolGen;

    // General stuff
    VariableInt nRun, nLumi; 
    VariableLong64 nEvt;
    VariableInt isKinReco;
    VariableInt isTopGen, topProdMode;
    VariableFloat recoWeight;
    VariableFloat genWeight;
    
private:

    class EventShapeVariables;
    class FoxWolframMoments;

    /// Exactly what it says 
    double ifOk(const double& vCheck) const;

    /// Compute and fill the spin correlation variables
    void compSpinCorr(const TLorentzVector& p4_pLep, const TLorentzVector& p4_aLep,
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

                      VariableFloat& phi0, VariableFloat& phi1);

    /// Compute and fill "tt event shapes"
    void compEvtShp(const TLorentzVector& p4_pLep, const TLorentzVector& p4_aLep,
                    const TLorentzVector& p4_pNeu, const TLorentzVector& p4_aNeu,
                    const TLorentzVector& p4_pBot, const TLorentzVector& p4_aBot,
                    VariableFloat& TT_sph, VariableFloat& TT_apl, VariableFloat& TT_cir, VariableFloat& TT_iso,
                    VariableFloat& TT_C, VariableFloat& TT_D, VariableFloat& TT_spt,
                    VariableFloat& TT_H0, VariableFloat& TT_H1, VariableFloat& TT_H2,
                    VariableFloat& TT_H3, VariableFloat& TT_H4,
                    VariableFloat& TT_R1, VariableFloat& TT_R2, VariableFloat& TT_R3, VariableFloat& TT_R4);

};

/// Implementation for event shape variables copied off the DLBDTMvaVariablesEventClassification from ttH framework
class VariablesPhiTT::EventShapeVariables{

public:

    /// Constructor from XYZ coordinates
    explicit EventShapeVariables(const std::vector<TLorentzVector>& inputVectors);

    /// Default destructor  
    ~EventShapeVariables(){};

    /// The return value is 1 for spherical events and 0 for events linear in r-phi
    /// Needs the number of steps to determine how fine the granularity of the algorithm in phi should be
    double isotropy(const unsigned int numberOfSteps = 1000) const;

    /// The return value is 1 for spherical and 0 linear events in r-phi
    double circularity(const unsigned int numberOfSteps = 1000) const;

    /// 1.5*(v1+v2) where 0<=v1<=v2<=v3 are the eigenvalues of the momemtum tensor 
    /// sum{p_j[a]*p_j[b]}/sum{p_j**2} normalized to 1. Return values are 1 for spherical, 3/4 for 
    /// plane and 0 for linear events
    double sphericity(const double& r = 2.) const;

    /// 1.5*v1 where 0<=v1<=v2<=v3 are the eigenvalues of the momemtum tensor 
    /// sum{p_j[a]*p_j[b]}/sum{p_j**2} normalized to 1. Return values are 0.5 for spherical and 0 
    /// for plane and linear events
    double aplanarity(const double& r = 2.) const;

    /// 3.*(v1*v2+v1*v3+v2*v3) where 0<=v1<=v2<=v3 are the eigenvalues of the momemtum tensor 
    /// sum{p_j[a]*p_j[b]}/sum{p_j**2} normalized to 1. Return value is between 0 and 1 
    /// and measures the 3-jet structure of the event (C vanishes for a "perfect" 2-jet event)
    double C(const double& r = 2.) const;

    /// 27.*(v1*v2*v3) where 0<=v1<=v2<=v3 are the eigenvalues of the momemtum tensor 
    /// sum{p_j[a]*p_j[b]}/sum{p_j**2} normalized to 1. Return value is between 0 and 1 
    /// and measures the 4-jet structure of the event (D vanishes for a planar event)
    double D(const double& r = 2.) const;

    /// 2.*v2/(v1+v2) where 0<=v1<=v2<=v3 are the eigenvalues of the momemtum tensor. Return value is between 0 and 1
    /// and measures "isotropic" structore for a value of 1 and "pencile-like" limit for value of 0
    double transSphericity(const double& r = 2.) const;

    /// Calculate the Fox-Wolfram moment for a given order
    double H(const int order) const;

    /// Ratio between Fox-Wolfram moment with given order to the 0-th Fox Wolfram moment
    double R(const int order) const;
    
private:

    /// Method used to convert vector of TLorentzVector to XYZVector type
    std::vector<ROOT::Math::XYZVector> transformInputVectors(const std::vector<TLorentzVector>& jets) const;

    /// Helper function to fill the 3 dimensional vector of eigen-values;
    /// The largest (smallest) eigen-value is stored at index position 0 (2)
    TVectorD compEigenValues(const double& r = 2.) const;

    /// Helper function to fill the 3 dimensional momentum tensor from the inputVectors
    TMatrixDSym compMomentumTensor(const double& r = 2.) const;

    /// Caching of input vectors
    const std::vector<ROOT::Math::XYZVector> inputVectors_;
};


class VariablesPhiTT::FoxWolframMoments{

public:

    FoxWolframMoments(const std::vector<ROOT::Math::XYZVector>& inputVectors, const int maxorder = 4);
    ~FoxWolframMoments(){};

    double H(const int order) const {return fwArray_(order);}
    double R(const int order) const;

private:

    void compute(const std::vector<ROOT::Math::XYZVector>& inputVectors);
    double legendre(const int l, const int m, const double& x) const;

    const int nMoment_;
    TVectorD fwArray_;
    TVectorD sumArray_;
};

#endif







