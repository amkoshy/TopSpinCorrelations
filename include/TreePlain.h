#ifndef TREEPLAIN_H
#define TREEPLAIN_H

#include <TTree.h>
#include <TBranch.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <memory>
#include <map>
#include "boost/variant.hpp"
//#include "reflection.h"
//#include "classesFwd.h"
#include "AnalysisConfig.h"
#include "analysisStructs.h"
#include "../../common/include/classesFwd.h"
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


// class TTree;
// class TBranch;
// class TChain;
// class TopAnalysisPlain;
// class EventMetadata;
// class RecoObjects;
//class CommonGenObjects;
// class TopGenObjects;
// class KinematicReconstructionSolutions;

// Simple structure to store info about trees which need to be created later
struct TreePlainSetup
{
  bool isGen=false;
  TString name="defaultSetupName";
  TString title="defaultSetupName";

  bool fillRecoTree=false;
  bool fillSimpleRecoTree=false;
  bool fillGenTree=false;

  bool isTTbarSample=false;
  bool isSystSample=false;
  bool isSystVariation=false;

};



// namespace ttbar
// {
//   class RecoObjectIndices;
//   class GenObjectIndices;
//   class GenLevelWeights;
//   //class RecoLevelWeights;
// }

struct AnalysisBaseEntry
{
    const AnalysisConfig* analysisConfig;
    const int* eventCounter;
    const AnalysisBase* analysisBase;
    // const TopAnalysisPlain* analysisBase;
    const EventMetadata* eventMetadata;
    const RecoObjects* recoObjects;
    const CommonGenObjects* commonGenObjects;
    const TopGenObjects* topGenObjects;
    const TopPseudoObjects* topPseudoObjects;
    const KinematicReconstructionSolutions* kinematicReconstructionSolutions;
    const LooseKinRecoSolution* looseKinRecoSolution;
    // const RecoObjectIndices* recoObjectIndices;
    // const GenObjectIndices* genObjectIndices;
    // const GenLevelWeights* genLevelWeights;
    const ttbar::RecoObjectIndices* recoObjectIndices;
    const ttbar::GenObjectIndices* genObjectIndices;
    const ttbar::GenLevelWeights* genLevelWeights;
    // const ttbar::RecoLevelWeights*;
    const double* weight;
    const double* weightLooseKinReco;
    float lepSF[11];
    float kinRecoSF[3];
    // float trigSF[5];
    float trigSF[3];
    // float btagSF[13];
    float nominalBtagSF;
    float pileupW[3];
    float l1prefiringW[3];
    Long64_t entryNumber;

    const std::vector<int>* jetIndices;
    const VLV* jets;
    const std::vector<int>* jetFlavour;
    std::vector<float>* btagDiscriminators;
};

struct TreePlainVars
{
    //***********************************************************
    //******************** tree variables ***********************
    //***********************************************************
    int eventCounter = 0;
    TBranch* b_eventCounter = NULL;

    float weight = 1.;

    bool inFiducialPhaseSpace = 0;
    TBranch* b_inFiducialPhaseSpace = NULL;

    // variables
    float ptt = -999.;
    float ptat = -999.;
    float pttLead = -999.;
    float pttNLead = -999.;
    float pttTTRestFrame = -999.;
    float ptatTTRestFrame = -999.;
    float yt = -999.;
    float yat = -999.;
    float ytLead = -999.;
    float ytNLead = -999.;
    float mtt = -999.;
    float ytt = -999.;
    float pttt = -999.;
    float dphitt = -999.;
    float dytt = -999.;
    float detatt = -999.;
    float recoPartonMomFraction = -999.;
    float recoAntipartonMomFraction = -999.;

    //particle-level cross sections
    //reco-level
    float ptl = -999.;
    float ptal = -999.;
    float ptlLead = -999.;
    float ptlNLead = -999.;
    float etal = -999.;
    float etaal = -999.;
    float etalLead = -999.;
    float etalNLead = -999.;
    float detall = -999.;
    float dphill = -999.;
    float yll = -999.;
    float ptll = -999.;
    float mll = -999.;
    float ptllbb = -999.;
    float yllbb = -999.;
    float mllbb = -999.;
    float ptllbbmet = -999.;
    float yllbbmet = -999.;
    float mllbbmet = -999.;
    float ptbLead = -999.;
    float ptbNLead = -999.;
    float etabLead = -999.;
    float etabNLead = -999.;
    float ybb = -999.;
    float ptbb = -999.;
    float mbb = -999.;
    float mww = -999.;
    float minmlbdef1 = -999.;
    float minmlbdef2 = -999.;

    //gen-level
    float visptt = -999.;
    float visptat = -999.;
    float vispttLead = -999.;
    float vispttNLead = -999.;
    float vispttTTRestFrame = -999.;
    float visptatTTRestFrame = -999.;
    float visyt = -999.;
    float visyat = -999.;
    float visytLead = -999.;
    float visytNLead = -999.;
    float vismtt = -999.;
    float visytt = -999.;
    float vispttt = -999.;
    float visdphitt = -999.;
    float visdytt = -999.;
    float visdetatt = -999.;
    float visrecoPartonMomFraction = -999.;
    float visrecoAntipartonMomFraction = -999.;

    float visptl = -999.;
    float visptal = -999.;
    float visptlLead = -999.;
    float visptlNLead = -999.;
    float visetal = -999.;
    float visetaal = -999.;
    float visetalLead = -999.;
    float visetalNLead = -999.;
    float visdetall = -999.;
    float visdphill = -999.;
    float visyll = -999.;
    float visptll = -999.;
    float vismll = -999.;
    float visptllbb = -999.;
    float visyllbb = -999.;
    float vismllbb = -999.;
    float visyllbbmet = -999.;
    float visptllbbmet = -999.;
    float vismllbbmet = -999.;
    float visptbLead = -999.;
    float visptbNLead = -999.;
    float visetabLead = -999.;
    float visetabNLead = -999.;
    float visybb = -999.;
    float visptbb = -999.;
    float vismbb = -999.;
    float vismww = -999.;
    float visminmlb = -999.;

    // int njet;

    int njdefIso04Pt30 = -999;
    int njdefIso04Pt40 = -999;
    int njdefIso04Pt50 = -999;
    int njdefIso04Pt75 = -999;
    int njdefIso04Pt100 = -999;
    int njdefIso04Pt150 = -999;

    int njinetattbarIso04Pt30 = -999;
    int njinetattbarIso04Pt40 = -999;
    int njinetattbarIso04Pt50 = -999;
    int njinetattbarIso04Pt75 = -999;
    int njinetattbarIso04Pt100 = -999;
    int njinetattbarIso04Pt150 = -999;

    float detatej1_detattbar = -999.;

    int nvtx = -999;
    int nbj = -999;
    int naj = -999;

    float met = -999.;
    float ptj1 = -999.;
    float etaj1 = -999.;
    float ptj2 = -999.;
    float etaj2 = -999.;

    //extra jets (1 leading, 2 trailing)
    float ptej1 = -999.;
    float etaej1 = -999.;
    float ptej2 = -999.;
    float etaej2 = -999.;


    float mttmeas = -999.;
    float yttmeas = -999.;
    float ptttmeas = -999.;
    float ttbaremeas = -999.;
    float ttbaremeasnomet = -999.;
    float ttbarzmeas = -999.;
    float ttbare = -999.;
    float ttbarz = -999.;
    float trz = -999.;
    float nnbare = -999.;
    float nnbarz = -999.;
    float llbare = -999.;
    float llbarz = -999.;



    float matchedPt = -999.;
    float matchedR = -999.;
    float matchedPtFullKR[2];
    float matchedRFullKR[2];

    // float mttKR2 = -999.;
    // float yttKR2 = -999.;
    // float mwwKR2 = -999.;
    // float ptttKR2 = -999.;
    // float minmlbKR2 = -999.;
    // float mttKR6 = -999.;
    // float yttKR6 = -999.;
    // float ptttKR6 = -999.;
    // float mwwKR6 = -999.;
    // float minmlbKR6 = -999.;
    float mttKR9 = -999.;
    float yttKR9 = -999.;
    float ptttKR9 = -999.;
    float mwwKR9 = -999.;
    float minmlbKR9 = -999.;

    //ttbar+j1 vars
    float mttplusej1 = -999.;
    float ptttplusej1 = -999.;
    float r_mttplusej1_mtt = -999.;
    float r_ptej1_mtt = -999.;
    float mean_ptt_pttbar_ptej1 = -999.;

    // reco and gen level
    TLorentzVector lvTop = TLorentzVector(0.,0.,0.,0.);
    TLorentzVector lvTopBar = TLorentzVector(0.,0.,0.,0.);
    TLorentzVector lvLeadTop = TLorentzVector(0.,0.,0.,0.);
    TLorentzVector lvNLeadTop = TLorentzVector(0.,0.,0.,0.);
    TLorentzVector lvTopTTRestFrame = TLorentzVector(0.,0.,0.,0.);
    TLorentzVector lvTopBarTTRestFrame = TLorentzVector(0.,0.,0.,0.);
    TLorentzVector lvLep = TLorentzVector(0.,0.,0.,0.);
    TLorentzVector lvLepBar = TLorentzVector(0.,0.,0.,0.);
    TLorentzVector lvBot = TLorentzVector(0.,0.,0.,0.);
    TLorentzVector lvBotBar = TLorentzVector(0.,0.,0.,0.);
    TLorentzVector lvTTbarPlusLeadingExtraJet = TLorentzVector(0.,0.,0.,0.);
    TLorentzVector lvMET = TLorentzVector(0.,0.,0.,0.);

    // reco level
    TLorentzVector lvLeadLep = TLorentzVector(0.,0.,0.,0.);
    TLorentzVector lvNLeadLep = TLorentzVector(0.,0.,0.,0.);
    TLorentzVector lvLeadBot = TLorentzVector(0.,0.,0.,0.);
    TLorentzVector lvNLeadBot = TLorentzVector(0.,0.,0.,0.);

    // Leading and Trailing extra jets
    std::pair<TLorentzVector, TLorentzVector> extraJets_j1_j2 = std::pair<TLorentzVector, TLorentzVector> (TLorentzVector(0.,0.,0.,0.),TLorentzVector(0.,0.,0.,0.));
    //unsigned char indBot;
    //unsigned char indBotBar;
    //
    // std::vector<TLorentzVector> vlvJets;
    // std::vector<float> vdJetsBDiscr;
    //
    TVector2 v2MET = TVector2(0.,0.);

    // reco level only
    // std::vector<TLorentzVector> vlvBJets;
    // int nPrimVtx;

    // pileup
    // float wPileup;
    float var_PU_U = -999.;
    float var_PU_D = -999.;

    // triggers
    // float wTrig;
    float var_TRIG_U = -999.;
    float var_TRIG_D = -999.;
    // float var_TRIG_ETA_U = -999.;
    // float var_TRIG_ETA_D = -999.;

   // leptonID
   // float wLept;
   float var_LEPT_U = -999.;
   float var_LEPT_D = -999.;

   float var_ELE_ID_U = -999.;
   float var_ELE_ID_D = -999.;
   float var_ELE_RECO_U = -999.;
   float var_ELE_RECO_D = -999.;

   float var_MUON_ID_U = -999.;
   float var_MUON_ID_D = -999.;
   float var_MUON_ISO_U = -999.;
   float var_MUON_ISO_D = -999.;

    // b-tagging
    // float wBTag;
    float var_BTAG_U = -999.;
    float var_BTAG_D = -999.;
    float var_BTAG_LJET_U = -999.;
    float var_BTAG_LJET_D = -999.;
    // float var_BTAG_PT_U = -999.;
    // float var_BTAG_PT_D = -999.;
    // float var_BTAG_ETA_U = -999.;
    // float var_BTAG_ETA_D = -999.;
    // float var_BTAG_LJET_PT_U = -999.;
    // float var_BTAG_LJET_PT_D = -999.;
    // float var_BTAG_LJET_ETA_U = -999.;
    // float var_BTAG_LJET_ETA_D = -999.;

    // kin reco
    // float wKinReco;
    float var_KIN_U = -999.;
    float var_KIN_D = -999.;

    // JES, JER (optional)
    // std::vector<float> vwJEC;
    // std::vector<float> vwJECUnc;
    // std::vector<float> vwJER;
    // std::vector<float> vwJERUnc;
    //std::vector<float> vwJERUncU;
    //std::vector<float> vwJERUncD;

    // generator level only
    // int TTBarDecay;
    //float weightGen;

    // ME scales
    float var_MERENSCALE_U = -999.;
    float var_MERENSCALE_D = -999.;
    float var_MEFACSCALE_U = -999.;
    float var_MEFACSCALE_D = -999.;
    float var_MESCALE_U = -999.;
    float var_MESCALE_D = -999.;
    float var_BFRAG_U = -999.;
    float var_BFRAG_C = -999.;
    float var_BFRAG_D = -999.;
    float var_BFRAG_P = -999.;
    float var_BSEMILEP_U = -999.;
    float var_BSEMILEP_D = -999.;
    float var_PDF_ALPHAS_U = -999.;
    float var_PDF_ALPHAS_D = -999.;
    
    //PS Weights
    float var_PSSCALE_WEIGHT_4_U = -999.;
    float var_PSSCALE_WEIGHT_4_D = -999.;
    float var_PSSCALE_WEIGHT_5_U = -999.;
    float var_PSSCALE_WEIGHT_5_D = -999.;

    float var_L1PREFIRING_U = -999.;
    float var_L1PREFIRING_D = -999.;

    // conditional content
    // std::vector<float>  weightPDF;
    // float var_PDF[160];
    std::vector<float> var_PDF= std::vector<float>(160, 1.);

};

class TreePlain// : public TTree
{

    template <class T> void AddBranchVerbosely(const char* name, T* ptr);
    // void FillReco(const AnalysisBaseEntry* entry);
    // void FillGen(const AnalysisBaseEntry* entry, const float deltaRLepJet = -1.0);
    int RecoTTBarDecayFromPdgId(const int pdgid1, const int pdgid2);
    void CalculateTtbarXsecVars(TreePlainVars &outputVars, const TLorentzVector &t, const TLorentzVector &tbar);

    int SelectExtraGenJets(const VLV* jets, const std::vector<TLorentzVector> jetsOmit,
                                      float& extraJetPtMin, float& extraJetEtaAbsMax, float& extraJetDeltaRMin);
    // int SelectExtraGenJets(const VLV* jets, const std::vector<TLorentzVector> jetsOmit,
    //                                   float& extraJetPtMin, float& extraJetEtaAbsMax, float& extraJetDeltaRMin, int& extraJetNLead);
    // int TreePlain::SelectExtraGenJets(const VLV* jets, const std::vector<TLorentzVector> jetsOmit,
    //                                   float& extraJetPtMin, float& extraJetEtaAbsMax, float& extraJetDeltaRMin, int& extraJetNLead,
    //                                   onst TLorentzVector& ttbar, float* outMttj, float* outYttj, float* outRhoj, float* maxptjPtr);
    /**
     * @brief Select_N_extraJets_and_N_j_inside_ttbar: Select at the same time number of extra jets and the number of extra jets inside the ttbar (in rapidity)
     * @param njet: number of extra jets (reference to the variable to fill)
     * @param njet_inside: number of extra jets inside the ttbar (in rapidity) (reference to the variable to fill)
     * @param jets: VLV of jets objects
     * @param jetsOmit: std::vector<TLorentzVector> of jets to exclude from the selection (ttbar b jet, etc)
     * @param extraJetPtMin: min Pt required for the jets
     * @param extraJetEtaAbsMax: max Eta allowed for the jets
     * @param extraJetDeltaRMin: require min isolation from jets to exclude
     */
    void Select_N_extraJets_and_N_j_inside_ttbar(int &njet, int &njet_inside, const VLV* jets, const std::vector<TLorentzVector> jetsOmit, float eta_t, float eta_at,
                                                 float& extraJetPtMin, float& extraJetEtaAbsMax, float& extraJetDeltaRMin);
    /**
     * @brief SelectLeadingAndTrailingExtraJet
     * @param jets: VLV of jets objects
     * @param jetsOmit: std::vector<TLorentzVector> of jets to exclude from the selection (ttbar b jet, etc)
     * @param extraJetPtMin: min Pt required for the jets
     * @param extraJetEtaAbsMax: max Eta allowed for the jets
     * @param extraJetDeltaRMin: require min isolation from jets to exclude
     * @return std::pair<TLorentzVector, TLorentzVector> Leading(.first) and Trailing(.second) jets
     */
    //std::pair<TLorentzVector, TLorentzVector> SelectLeadingAndTrailingExtraJet(const VLV* jets, const std::vector<TLorentzVector> jetsOmit, float& extraJetPtMin, float& extraJetEtaAbsMax, float& extraJetDeltaRMin);
    void Compute_extraJets_observables(int &njet, int &njet_inside, std::pair<TLorentzVector, TLorentzVector> & leadingAndTrailingExtraJet,
                                       const VLV* jets, const std::vector<TLorentzVector> jetsOmit, float eta_t, float eta_at,
                                       float& extraJetPtMin, float& extraJetEtaAbsMax, float& extraJetDeltaRMin);
    
    std::pair<bool, std::pair<float, float> > IsMatched(const TLorentzVector& lv1, const TLorentzVector& lv2, const float maxDeltaPt = -1.0, const float maxDeltaR = -1.0);

    std::pair<float, float> CalculateMatchingPtR(const AnalysisBaseEntry* entry, const std::pair<TLorentzVector, TLorentzVector>& bbbar);
    std::pair<float, float> CalculateMatchingPtRFullKR(const AnalysisBaseEntry* entry, const TLorentzVector& b, const TLorentzVector& bbar);
    float GetBTagWeightForSyst(const AnalysisBaseEntry* entry, const TString& systName);
    
    void AddBasicBranches();
    void AddParticleLevelBranches();
    void AddMeasureBranches();
    void AddGenVarBranches();
    void AddPUVarBranches();
    void AddRecoVarBranches();
    void AddLooseKinRecoBranches();
    void AddJetBranches();
    void AddPseudoBranches();

    void FillBasicBranches(const AnalysisBaseEntry* entry);
    void FillMeasureBranches(const AnalysisBaseEntry* entry);
    void FillGenVarBranches(const AnalysisBaseEntry* entry);
    void FillPUVarBranches(const AnalysisBaseEntry* entry);
    void FillRecoVarBranches(const AnalysisBaseEntry* entry);
    void FillLooseKinRecoBranches(const AnalysisBaseEntry* entry);
    void FillJetBranches(const AnalysisBaseEntry* entry);
    void FillPseudoBranches(const AnalysisBaseEntry* entry);
    
  public:
    TreePlainVars Vars;

    // TreePlain(const TString& name, const TString& title, const bool flagGen);
    TreePlain(const TreePlainSetup& setup_);

// FIXME PORT==============================================

    short* cutVarShort = NULL;
    template<class T> void* SetBranchAddressVerbosely(TTree* tree, const char* name, T& ptr, TBranch** br);

    // tree to read
    TreePlain(const TString& name, const bool flagGen, const TString var = "", const bool DoParticle = false);
    ~TreePlain();
    void Add(const TString& name, const double w = 1.0, const bool doWeightVar = false, const int sampleType = -1);
    void InitReading();
    void InitVar(const TString& nameVar, const TString& nameCut, std::map<const TreePlain*, boost::variant<std::shared_ptr<int>, std::shared_ptr<float>, std::vector<std::shared_ptr<float>> >> & zTrees);

    TChain* TChainPtr = NULL;

    std::vector<bool> zDoWeightVar;
    std::vector<double> zWeight;
    std::vector<int> zSampleType;
    void AddWeight(const double w);
    void DoWeightVar(const bool doWeightVar);
    float* weightVar = NULL;
    TString weightVarName = "";
    std::map<const TString, boost::variant<std::shared_ptr<int>, std::shared_ptr<float>, std::vector<std::shared_ptr<float>> >> initVars;
    bool addFiducialInfoBranch = false;

    double Weight(bool flagApplyCut) const;
    double GetWeightOfTree(const int n) const {return zWeight[n];}
    int SampleType(int n = -1) const;
// PORT end ===============================================

    void InitWriteTree();

    // void Fill(const AnalysisBaseEntry* entry, const float deltaRLepJet = -1.0);
    void Fill(const AnalysisBaseEntry* entry);

    // pointer to ROOT TTree
    TTree* TTreePtr = NULL;

    // generator level flag (true or false for MC, false for data)
    const bool isGen_;

    // identifier of the analysis selection step
    int StepId = -1;

    std::vector<float> vExtraJetPtMin = { 30.0, 40.0, 50.0, 75.0, 100.0, 150.0 };
    // int extraJetNLead = 2;
    float extraJetDeltaRMin = 0.4;
    float extraJetEtaAbsMax = 2.4;

    /// Set the btag scale factors
    void SetBtagScaleFactors_Var_AllSyst(const std::map<std::tuple<Systematic::Type, Systematic::Variation, int>, BtagScaleFactors*>* mapAllSystBTagScaleFactors);

    // content flags
    //const bool FlagEnableJets;
    // bool FlagEnablePDFWeights = false;
    // bool FlagEnableSystExpWeights = false;
    // bool FlagEnableSystModWWeights = false;
    // bool FlagStoreJetCorr = false;
    // bool FlagStoreExtraVars = false;

    // new booleans
    bool fillRecoTree_;
    bool fillSimpleRecoTree_;
    bool fillGenTree_;

    bool isTTbarSample_;
    bool isSystSample_;
    bool isSystVariation_;

    bool fillBasicVariables_;
    bool fillParticleLevelVariables_;
    bool fillMeasureVariables_;
    bool fillLooseKinRecoVariables_;
    bool fillPUVarWeights_;
    bool fillJetVariables_;
    bool fillPseudoVariables_;
    bool fillGenVariationWeights_;
    bool fillRecoVariationWeights_;

  protected:
    const std::map<std::tuple<Systematic::Type, Systematic::Variation, int>, BtagScaleFactors*>* mapAllSystBTagScaleFactors_;

    /// Get weight due to b-tagging efficiency MC-to-data scale factors
     float weightBtagSF_Var(const std::vector<int>& jetIndices,
                         const VLV& jets, const std::vector<int>& jetFlavour,
                         const std::vector<float>& btagDiscriminators, const TString& sysKey)const;


};

#endif // TREEPLAIN_H
