#ifndef MiniTreeTopAnalysis_h
#define MiniTreeTopAnalysis_h

#include "MiniTree.h"
#include "TopAnalysis.h"
#include "AnalysisConfig.h"
#include "TFile.h"


class MiniTreeTopAnalysis : public TopAnalysis
{
  using TopAnalysis::TopAnalysis;

  public:

    /// Set the lepton scale factors
    void SetLeptonScaleFactors_EleUp(const LeptonScaleFactors* const scaleFactors);
    void SetLeptonScaleFactors_EleDn(const LeptonScaleFactors* const scaleFactors);
    void SetLeptonScaleFactors_MuonUp(const LeptonScaleFactors* const scaleFactors);
    void SetLeptonScaleFactors_MuonDn(const LeptonScaleFactors* const scaleFactors);

    /// Set the trigger scale factors
    void SetTriggerScaleFactors_Up(const TriggerScaleFactors* const scaleFactors);
    void SetTriggerScaleFactors_Dn(const TriggerScaleFactors* const scaleFactors);

    /// Set the btag scale factors for all systematic variations
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

   void SetRunViaTau(const bool runViaTau);
   void SetTTBarSample(const bool& isTTBarMC);
   void SetSystSample(const bool& isSystMC);
   void SetSystVarSample(const bool& isSystVar);


    virtual void Begin(TTree*);

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
    void SetTreeOutputDirectory(const TString& treeOutputDirectory);


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

        /// Pointer to lepton scale factors instance
        const LeptonScaleFactors* leptonScaleFactors_EleUp_;
        const LeptonScaleFactors* leptonScaleFactors_EleDn_;
        const LeptonScaleFactors* leptonScaleFactors_MuonUp_;
        const LeptonScaleFactors* leptonScaleFactors_MuonDn_;

        /// Pointer to trigger scale factors instance
        const TriggerScaleFactors* triggerScaleFactors_Up_;
        const TriggerScaleFactors* triggerScaleFactors_Dn_;

        const std::map<std::tuple<Systematic::Type, Systematic::Variation, int>, BtagScaleFactors*>* mapAllSystBTagScaleFactors_;

        /// Get weight due to trigger efficiency MC-to-data scale factors
        float weightTriggerSF_Up(const int leptonXIndex, const int leptonYIndex,
                                 const VLV& allLeptons, const std::vector<int>& lepPdgId)const;
        float weightTriggerSF_Dn(const int leptonXIndex, const int leptonYIndex,
                                 const VLV& allLeptons, const std::vector<int>& lepPdgId)const;


    private:
        std::vector<MiniTree*> vecMiniTree_;
        TFile* fileOut_ = NULL;
        TString dirOut_ = "miniTree";

        bool isTTbarSample_;
        bool isSystSample_;
        bool isSystVariation_;

};

#endif
