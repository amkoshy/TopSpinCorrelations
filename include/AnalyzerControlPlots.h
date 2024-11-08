#ifndef AnalyzerControlPlots_h
#define AnalyzerControlPlots_h

#include <map>
#include <vector>

class TString;
class TH1;

#include "AnalyzerBaseClass.h"

class EventMetadata;
class RecoObjects;
class CommonGenObjects;
class TopGenObjects;
class KinematicReconstructionSolutions;
class LooseKinRecoSolution;
namespace ttbar{
    class RecoLevelWeights;
    class GenLevelWeights;
    class GenObjectIndices;
    class RecoObjectIndices;
}



/// Class for basic histograms that are filled simultaneously for any step
class AnalyzerControlPlots : public AnalyzerBaseClass{

public:
    
    /// Constructor
    AnalyzerControlPlots(const std::vector<TString>& selectionStepsNoCategories);

    /// Destructor
    ~AnalyzerControlPlots(){}



private:

    /// Book all histograms for given selection step
    virtual void bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram);

    /// Fill all histograms for given selection step
    virtual void fillHistos(const EventMetadata& eventMetadata, const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                            const TopGenObjects& topGenObjects,
                            const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
			    const LooseKinRecoSolution& looseKinematicReconstructionSolution, 
                            const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                            const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                            const double& weight, const double& weightLooseKinReco, const TString& step,
                            std::map<TString, TH1*>& m_histogram);
};




#endif
