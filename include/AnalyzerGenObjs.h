#ifndef AnalyzerGenObjs_h
#define AnalyzerGenObjs_h

#include "AnalyzerBaseClass.h"
#include "analysisHelpers.h"
#include "analysisStructsFwd.h"

#include "../../common/include/classes.h"
#include "../../common/include/AnalysisBase.h"
#include "../../common/include/classesFwd.h"

class AnalyzerGenObjs : public AnalyzerBaseClass
{
 public:
  explicit AnalyzerGenObjs(
    const std::vector<TString>& selectionStepsNoCategories
  );

  virtual ~AnalyzerGenObjs(){}

 protected:

  virtual void bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram);

  virtual void fillHistos(
      const EventMetadata& eventMetadata, const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                              const TopGenObjects& topGenObjects,
                              const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                              const LooseKinRecoSolution& looseKinematicReconstructionSolution,
                              const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                              const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                              const double& weight, const double& weightLooseKinReco, const TString& step,
                              std::map<TString, TH1*>& m_histogram);

  virtual void endJob(const TString&, const Channel::Channel&, const Systematic::Systematic&) {}
};

#endif
