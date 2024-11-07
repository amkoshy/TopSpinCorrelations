#ifndef TopAnalysisPlain_h
#define TopAnalysisPlain_h

#include "TreePlain.h"
#include "TopAnalysis.h"
// #include "../../common/include/analysisObjectStructs.h"
// #include "analysisStructs.h"
#include "TFile.h"

// class analysisObjectStructs;

// Class for plain tree production (2D analysis).
// Overrides a few methods of the base class to alter treatment in some places and to store more variables.

class TopAnalysisPlain : public TopAnalysis
{
  using TopAnalysis::TopAnalysis;

  public:
    virtual void Begin(TTree*);

    /// Set the lepton scale factors
    // void SetLeptonScaleFactors_Up(const LeptonScaleFactors* const scaleFactors);
    // void SetLeptonScaleFactors_Dn(const LeptonScaleFactors* const scaleFactors);

    /// Set the trigger scale factors
    void SetTriggerScaleFactors_Up(const TriggerScaleFactors* const scaleFactors);
    void SetTriggerScaleFactors_Dn(const TriggerScaleFactors* const scaleFactors);
    // void SetTriggerScaleFactors_EtaUp(const TriggerScaleFactors* const scaleFactors);
    // void SetTriggerScaleFactors_EtaDn(const TriggerScaleFactors* const scaleFactors);

    /// Set the btag scale factors
    void SetBtagScaleFactors_Var_AllSyst(const std::map<std::tuple<Systematic::Type, Systematic::Variation, int>, BtagScaleFactors*>* mapAllSystBTagScaleFactors);

    /// Set the kinematic reconstruction
    void SetKinematicReconstruction_Up(const KinematicReconstruction* const kinematicReconstruction,
                                      const KinematicReconstructionScaleFactors* const kinematicReconstructionScaleFactors);
    /// Set the kinematic reconstruction
    void SetKinematicReconstruction_Dn(const KinematicReconstruction* const kinematicReconstruction,
                                      const KinematicReconstructionScaleFactors* const kinematicReconstructionScaleFactors);


    /// Set the loose kinematic reconstruction
    void SetLooseKinReco_Up(const LooseKinReco* const looseKinReco,
			                     const LooseKinRecoScaleFactors* const looseKinRecoScaleFactors);
    /// Set the loose kinematic reconstruction
    void SetLooseKinReco_Dn(const LooseKinReco* const looseKinReco,
			                     const LooseKinRecoScaleFactors* const looseKinRecoScaleFactors);


    // does not work: need to be called from virtual method
    //float weightLeptonSF(const int leadingLeptonIndex, const int nLeadingLeptonIndex,
    //                                    const VLV& allLeptons, const std::vector<int>& lepPdgId) const;
    virtual bool failsTopGeneratorSelection(const Long64_t& entry) const;
    virtual void fillAll(const std::string& selectionStep,
                 const EventMetadata& eventMetadata,
                 const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                 const TopGenObjects& topGenObjects, const TopPseudoObjects& topPseudoObjects,
                 const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                 const LooseKinRecoSolution& looseKinRecoSolution,
                 const ttbar::GenObjectIndices& genObjectIndices, const ttbar::RecoObjectIndices& recoObjectIndices,
                 const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                 const double& defaultWeight,const double& defaultWeightLooseKinReco)const;
    virtual void Terminate();

    void SetPlainTreeOutputDirectory(const TString& plainTreeOutputDirectory);

    virtual void SetRunViaTau(const bool& runViaTau);

    void SetRecoTree(const bool& makeRecoTree);
    void SetSimpleRecoTree(const bool& makeSimpleRecoTree);
    void SetGenTree(const bool& makeGenTree);

    void SetTTBarSample(const bool& isTTBarMC);
    void SetSystSample(const bool& isSystMC);
    void SetSystVarSample(const bool& isSystVar);

    //std::vector<TreePlainSetup> VecTreePlainSetup;

  protected:

    // add instances of reweighting classes for syst variations

    /// Pointer to the kinematic reconstruction instance
    const KinematicReconstruction* kinematicReconstruction_Up_;
    const KinematicReconstruction* kinematicReconstruction_Dn_;

    /// Pointer to the kinematic reconstruction scale factors instance
    const KinematicReconstructionScaleFactors* kinematicReconstructionScaleFactors_Up_;
    const KinematicReconstructionScaleFactors* kinematicReconstructionScaleFactors_Dn_;

    /// Pointer to the loose kinematic reconstruction instance
    const LooseKinReco* looseKinReco_Up_;
    const LooseKinReco* looseKinReco_Dn_;

    /// Pointer to the loose kinematic reconstruction scale factors instance
    const LooseKinRecoScaleFactors* looseKinRecoScaleFactors_Up_;
    const LooseKinRecoScaleFactors* looseKinRecoScaleFactors_Dn_;

    /// Pointer to the pileup reweighter instance
    const PileupScaleFactors* pileupScaleFactors_Up_;
    const PileupScaleFactors* pileupScaleFactors_Dn_;

    /// Pointer to lepton scale factors instance
    const LeptonScaleFactors* leptonScaleFactors_Up_;
    const LeptonScaleFactors* leptonScaleFactors_Dn_;

    /// Pointer to trigger scale factors instance
    const TriggerScaleFactors* triggerScaleFactors_Up_;
    const TriggerScaleFactors* triggerScaleFactors_Dn_;
    const TriggerScaleFactors* triggerScaleFactors_EtaUp_;
    const TriggerScaleFactors* triggerScaleFactors_EtaDn_;

    /// Pointer to btag scale factors instance
    const std::map<std::tuple<Systematic::Type, Systematic::Variation, int>, BtagScaleFactors*>* mapAllSystBTagScaleFactors_;

    /// Pointer to jet energy resolution scale factors instance
    const JetEnergyScaleScaleFactors*    jetEnergyScaleScaleFactors_Up_; // JECs for standard jets (AK4)
    const JetEnergyScaleScaleFactors*    jetEnergyScaleScaleFactors_Dn_; // JECs for standard jets (AK4)
    const JetEnergyScaleScaleFactors*    jetEnergyScaleScaleFactors_L1only_Up_; // L1-only JECs for standard jets (AK4), needed in L1L2L3-L1 scheme for MET correction
    const JetEnergyScaleScaleFactors*    jetEnergyScaleScaleFactors_L1only_Dn_; // L1-only JECs for standard jets (AK4), needed in L1L2L3-L1 scheme for MET correction
    const JetEnergyScaleScaleFactors* fatjetEnergyScaleScaleFactors_Up_; // JECs for FatJets
    const JetEnergyScaleScaleFactors* fatjetEnergyScaleScaleFactors_Dn_; // JECs for FatJets
    const JetEnergyScaleScaleFactors* subjetEnergyScaleScaleFactors_Up_; // JECs for subjets of FatJets
    const JetEnergyScaleScaleFactors* subjetEnergyScaleScaleFactors_Dn_; // JECs for subjets of FatJets

    /// Get weight due to lepton efficiency MC-to-data scale factors
    float weightLeptonSF_Up(const int leadingLeptonIndex, const int nLeadingLeptonIndex,
                          const VLV& allLeptons, const std::vector<int>& lepPdgId)const;

    /// Get weight due to lepton efficiency MC-to-data scale factors. For electrons the supercluster eta is used instead of candidate eta.
    float weightLeptonSF_Up(const int leadingLeptonIndex, const int nLeadingLeptonIndex,
                          const VLV& allLeptons, const std::vector<int>& lepPdgId, const std::vector<float>& lepSCEta)const;
    /// Get weight due to lepton efficiency MC-to-data scale factors
    float weightLeptonSF_Dn(const int leadingLeptonIndex, const int nLeadingLeptonIndex,
                          const VLV& allLeptons, const std::vector<int>& lepPdgId)const;

    /// Get weight due to lepton efficiency MC-to-data scale factors. For electrons the supercluster eta is used instead of candidate eta.
    float weightLeptonSF_Dn(const int leadingLeptonIndex, const int nLeadingLeptonIndex,
                          const VLV& allLeptons, const std::vector<int>& lepPdgId, const std::vector<float>& lepSCEta)const;

    /// Get weight due to trigger efficiency MC-to-data scale factors
    float weightTriggerSF_Up(const int leptonXIndex, const int leptonYIndex,
                             const VLV& allLeptons, const std::vector<int>& lepPdgId)const;
    float weightTriggerSF_Dn(const int leptonXIndex, const int leptonYIndex,
                             const VLV& allLeptons, const std::vector<int>& lepPdgId)const;
    float weightTriggerSF_EtaUp(const int leptonXIndex, const int leptonYIndex,
                             const VLV& allLeptons, const std::vector<int>& lepPdgId)const;
    float weightTriggerSF_EtaDn(const int leptonXIndex, const int leptonYIndex,
                             const VLV& allLeptons, const std::vector<int>& lepPdgId)const;

     /// Get weight due to efficiency of kinematic reconstruction
     float weightKinReco_Up()const;
     float weightKinReco_Dn()const;

     /// Get weight due to efficiency of loose kinematic reconstruction
     float weightLooseKinReco_Up()const;
     float weightLooseKinReco_Dn()const;


  private:
    std::vector<TreePlain*> vecTreePlain_;
    TFile* fileOut_ = NULL;
    TString dirOut_ = "plainTree";

    bool fillRecoTree_;
    bool fillSimpleRecoTree_;
    bool fillGenTree_;

    bool isTTbarSample_;
    bool isSystSample_;
    bool isSystVariation_;

};

#endif
