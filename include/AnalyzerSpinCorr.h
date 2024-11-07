#ifndef AnalyzerSpinCorr_h
#define AnalyzerSpinCorr_h

#include <map>
#include <vector>
#include "../common/include/AnalysisBase.h"
#include "../../common/include/classesFwd.h"
#include "analysisStructsFwd.h"

class TString;
class TH1;

#include "AnalyzerBaseClass.h"

class EventMetadata;
class RecoObjects;
class CommonGenObjects;
class TopGenObjects;
class KinematicReconstructionSolutions;
class analysisConfig;
namespace ttbar{
    class RecoLevelWeights;
    class GenLevelWeights;
    class GenObjectIndices;
    class RecoObjectIndices;
}



/// Class for basic histograms that are filled simultaneously for any step
class AnalyzerSpinCorr : public AnalyzerBaseClass{


public:

    /// Constructor
    AnalyzerSpinCorr(const AnalysisConfig& analysisConfig,const double& leptonPtCut, const double& leptonEtaCut, const double& jetPtCut, const double& jetEtaCut, const double& genDeltaRLeptonJetCut, const std::vector<TString>& selectionStepsNoCategories);

    /// Destructor
    ~AnalyzerSpinCorr(){}


private:

    const double leptonPtCut_;
    const double leptonEtaCut_;
    const double jetPtCut_;
    const double jetEtaCut_;
    const double genDeltaRLeptonJetCut_;

    /// Book all histograms for given selection step
    virtual void bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram);

    /// Fill all histograms for given selection step
    virtual void fillHistos(const EventMetadata& eventMetadata,
                            const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                            const TopGenObjects& topGenObjects,
                            const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
			    const LooseKinRecoSolution& looseKinematicReconstructionSolution,
                            const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                            const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                            const double& weight, const double& weightLooseKinReco, const TString& step,
                            std::map<TString, TH1*>& m_histogram);

    virtual double getMT2Variable(const LV& Lepton, const LV& antiLepton, LV MET);
    float getJetHT(const std::vector<int>& jetIndices, const VLV& jets)const;

    const AnalysisConfig& analysisConfig_;
};




#endif
