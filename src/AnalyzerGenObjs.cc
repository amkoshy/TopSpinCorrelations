#define ANALYZERGENOBJS_CC
#include   "AnalyzerGenObjs.h"

#include <iostream>
#include <utility>
#include <string>

#include <TH1D.h>
#include <TH2D.h>
#include <Math/VectorUtil.h>
#include <TLorentzVector.h>

#include "analysisStructs.h"
#include "AnalysisConfig.h"

#include "../../common/include/analysisUtils.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/classes.h"
#include "../../common/include/KinematicReconstructionSolution.h"
#include "../../common/include/AnalysisBase.h"

AnalyzerGenObjs::AnalyzerGenObjs
(
  const std::vector<TString>& selectionStepsNoCategories
)
 : AnalyzerBaseClass("gen_", selectionStepsNoCategories)
{
  std::cout << "--- Beginning setting up \"AnalyzerGenObjs\"" << "\n";
  std::cout << "=== Finishing setting up \"AnalyzerGenObjs\"" << "\n\n";
}

void AnalyzerGenObjs::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
  TString name;

  name ="GenTop_pt";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step, "GenTop_pt", 180, 0, 900 ));
  name ="GenTop_eta";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenTop_eta", 60, -6, 6));
  name ="GenAntiTop_pt";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenAntiTop_pt", 180, 0, 900));
  name ="GenAntiTop_eta";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenAntiTop_eta", 60, -6, 6));
  name ="GenLepton_pt";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenLepton_pt", 180, 0, 900));
  name ="GenLepton_eta";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenLepton_eta", 60, -6, 6));
  name ="GenAntiLepton_pt";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenAntiLepton_pt", 180, 0, 900));
  name ="GenAntiLepton_eta";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenAntiLepton_eta", 60, -6, 6));
  name ="GenNeutrino_pt";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenNeutrino_pt", 180, 0, 900));
  name ="GenNeutrino_eta";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenNeutrino_eta", 60, -6, 6));
  name ="GenAntiNeutrino_pt";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenAntiNeutrino_pt", 180, 0, 900));
  name ="GenAntiNeutrino_eta";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenAntiNeutrino_eta", 60, -6, 6));
  name ="GenB_pt";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenB_pt", 180, 0, 900));
  name ="GenB_eta";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenB_eta", 60, -6, 6));
  name ="GenAntiB_pt";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenAntiB_pt", 180, 0, 900));
  name ="GenAntiB_eta";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenAntiB_eta", 60, -6, 6));
  name ="GenMet_pt";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenMet_pt", 180, 0, 900));
  name ="GenMet_eta";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenMet_eta", 60, -6, 6));
  name ="GenJet_1_pt";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenJet_1_pt", 180, 0, 900));
  name ="GenJet_1_eta";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenJet_1_eta", 60, -6, 6));
  name ="GenJet_2_pt";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenJet_2_pt", 180, 0, 900));
  name ="GenJet_2_eta";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenJet_2_eta", 60, -6, 6));
  name ="GenJet_3_pt";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenJet_3_pt", 180, 0, 900));
  name ="GenJet_3_eta";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenJet_3_eta", 60, -6, 6));
  name ="GenJet_4_pt";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenJet_4_pt", 180, 0, 900));
  name ="GenJet_4_eta";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenJet_4_eta", 60, -6, 6));
  name ="GenJet_5_pt";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenJet_5_pt", 180, 0, 900));
  name ="GenJet_5_eta";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenJet_5_eta", 60, -6, 6));
  name ="GenJet_6_pt";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenJet_6_pt", 180, 0, 900));
  name ="GenJet_6_eta";
  m_histogram[name] = this->store(new TH1D ( prefix_+name+step,  "GenJet_6_eta", 60, -6, 6));

  return;
}

void AnalyzerGenObjs::fillHistos(
    const EventMetadata& eventMetadata, const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                            const TopGenObjects& topGenObjects,
                            const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                            const LooseKinRecoSolution& looseKinematicReconstructionSolution,
                            const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                            const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                            const double& weight, const double& weightLooseKinReco, const TString& step,
                            std::map<TString, TH1*>& m_histogram)
{

  //silence unused parameter warnings
  (void) eventMetadata;
  (void) recoObjects;
  (void) commonGenObjects;
  (void) kinematicReconstructionSolutions;
  (void) looseKinematicReconstructionSolution;
  (void) recoObjectIndices;
  (void) genObjectIndices;
  (void) recoLevelWeights;
  (void) weight;
  (void) weightLooseKinReco;
  (void) step;

  const float weightToUse = genLevelWeights.trueLevelWeight_;

  if(topGenObjects.GenTop_)
  {
    m_histogram["GenTop_pt"] ->Fill(topGenObjects.GenTop_->pt() , weightToUse);
    m_histogram["GenTop_eta"]->Fill(topGenObjects.GenTop_->eta(), weightToUse);
  }

  if(topGenObjects.GenAntiTop_)
  {
    m_histogram["GenAntiTop_pt"]->Fill(topGenObjects.GenAntiTop_->pt() , weightToUse);
    m_histogram["GenAntiTop_eta"]->Fill(topGenObjects.GenAntiTop_->eta(), weightToUse);
  }

  if(topGenObjects.GenLepton_)
  {
    m_histogram["GenLepton_pt"]->Fill(topGenObjects.GenLepton_->pt() , weightToUse);
    m_histogram["GenLepton_eta"]->Fill(topGenObjects.GenLepton_->eta(), weightToUse);
  }

  if(topGenObjects.GenAntiLepton_)
  {
    m_histogram["GenAntiLepton_pt"]->Fill(topGenObjects.GenAntiLepton_->pt() , weightToUse);
    m_histogram["GenAntiLepton_eta"]->Fill(topGenObjects.GenAntiLepton_->eta(), weightToUse);
  }

  if(topGenObjects.GenNeutrino_)
  {
    m_histogram["GenNeutrino_pt"]->Fill(topGenObjects.GenNeutrino_->pt() , weightToUse);
    m_histogram["GenNeutrino_eta"]->Fill(topGenObjects.GenNeutrino_->eta(), weightToUse);
  }

  if(topGenObjects.GenAntiNeutrino_)
  {
    m_histogram["GenAntiNeutrino_pt"]->Fill(topGenObjects.GenAntiNeutrino_->pt() , weightToUse);
    m_histogram["GenAntiNeutrino_eta"]->Fill(topGenObjects.GenAntiNeutrino_->eta(), weightToUse);
  }

  if(topGenObjects.GenB_)
  {
    m_histogram["GenB_pt"]->Fill(topGenObjects.GenB_->pt() , weightToUse);
    m_histogram["GenB_eta"]->Fill(topGenObjects.GenB_->eta(), weightToUse);
  }

  if(topGenObjects.GenAntiB_)
  {
    m_histogram["GenAntiB_pt"]->Fill(topGenObjects.GenAntiB_->pt() , weightToUse);
    m_histogram["GenAntiB_eta"]->Fill(topGenObjects.GenAntiB_->eta(), weightToUse);
  }

  if(topGenObjects.GenMet_)
  {
    m_histogram["GenMet_pt"]->Fill(topGenObjects.GenMet_->pt() , weightToUse);
    m_histogram["GenMet_eta"]->Fill(topGenObjects.GenMet_->eta(), weightToUse);
  }

  if(topGenObjects.allGenJets_)
  {
    const auto& allGenJets = *(topGenObjects.allGenJets_);

    for(uint idx=0; idx<allGenJets.size(); ++idx)
    {
      if(idx > 5){ break; }

      m_histogram["GenJet_"+std::to_string(1+idx)+"_pt"]->Fill(allGenJets.at(idx).pt() , weightToUse);
      m_histogram["GenJet_"+std::to_string(1+idx)+"_eta"]->Fill(allGenJets.at(idx).eta(), weightToUse);
    }
  }

  return;
}
