#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <cmath>
#include <limits>
#include <iomanip>
#include <vector>

#include <TH1.h>
#include <TH2.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TMath.h>
#include <Math/VectorUtil.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TStyle.h>
#include <TIterator.h>
#include <TObject.h>

#include "TopAnalysis.h"
#include "AnalyzerSpinCorr.h"
#include "analysisStructs.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"
#include "../../common/include/KinematicReconstructionSolution.h"
#include "../../common/include/mt2_bisect.h"
#include "VariablesPhiTT.h"
#include "AnalysisConfig.h"

#include "../../common/include/ScaleFactors.h"
#include "../../common/include/AnalysisBase.h"

float AnalyzerSpinCorr::getJetHT(const std::vector<int>& jetIndices, const VLV& jets)const
{
  float result = 0.;
  for(const int index : jetIndices){
    const float pt = jets.at(index).pt();
    result += pt;
  }
  return result;
}


AnalyzerSpinCorr::AnalyzerSpinCorr(const AnalysisConfig& analysisConfig,const double& leptonPtCut, const double& leptonEtaCut, const double& jetPtCut, const double& jetEtaCut, const double& genDeltaRLeptonJetCut, const std::vector<TString>& selectionStepsNoCategories):
AnalyzerBaseClass("SpinCorr_", selectionStepsNoCategories),
analysisConfig_(analysisConfig),
leptonPtCut_(leptonPtCut),
leptonEtaCut_(leptonEtaCut),
jetPtCut_(jetPtCut),
jetEtaCut_(jetEtaCut),
genDeltaRLeptonJetCut_(genDeltaRLeptonJetCut)
{
    std::cout<<"--- Beginning setting up spin correlation histograms\n";
    std::cout<<"=== Finishing setting up spin correlation histograms\n\n";
}

void AnalyzerSpinCorr::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    
  TString name;

  // Kinematic Distributions Begin
  

  // Number of vertices 

  // name = "HypvertMulti";
  // m_histogram[name] = this->store(new TH1D ( name+step, "Primary Vertex Multiplicity", 50, 0, 50 ));

  // Leptons

  // name = "HypLeptonpT";
  // m_histogram[name] = store(new TH1D ( name+step, "Lepton Hypothesis pT", 1200,0,1200 ));
  // name = "HypAntiLeptonpT";
  // m_histogram[name] = store(new TH1D ( name+step, "AntiLepton Hypothesis pT", 1200,0,1200 ));
  // name = "HypLeptonEta";
  // m_histogram[name] = store(new TH1D ( name+step, "Lepton Eta", 200,-5,5 ));
  // name = "HypAntiLeptonEta";
  // m_histogram[name] = store(new TH1D ( name+step, "AntiLepton Hypothesis Eta", 200,-5,5 ));

  // name = "HypLeptonpT_AntiTopFrame";
  // m_histogram[name] = store(new TH1D ( name+step, "Lepton pT AntiTopFrame", 80, 0, 400 ));
  // name = "HypLeptonEta_AntiTopFrame";
  // m_histogram[name] = store(new TH1D ( name+step, "Lepton Eta AntiTopFrame", 100, -5, 5 ));
  // name = "HypAntiLeptonpT_TopFrame";
  // m_histogram[name] = store(new TH1D ( name+step, "AntiLepton pT TopFrame", 80, 0, 400 ));
  // name = "HypAntiLeptonEta_TopFrame";
  // m_histogram[name] = store(new TH1D ( name+step, "AntiLepton Eta TopFrame", 100, -5, 5 ));

  // name = "HypLLBarMass";
  // m_histogram[name] = store(new TH1D ( name+step, "Mass of LLbar System", 500, 0, 1000 ));
  // name = "HypLLBarpT";
  // m_histogram[name] = store(new TH1D ( name+step, "pT of LLbar System", 200, 0, 1000 ));
  name = "HypLLBarDPhi";
  m_histogram[name] = store(new TH1D ( name+step, "#Delta#phi(Lep, AntiLep)",3200, 0., 3.2 ));
  name = "HypLLBarDEta";
  m_histogram[name] = store( new TH1D( name+step, "#Delta #Eta (Lep, AntiLep)", 400, 0, 5));

  //Declaring some essential histograms in the extra-jet phase space:begin

  if(step.Contains("extrajets")){
  name = "HypLLBarMass";
  m_histogram[name] = store(new TH1D ( name+step, "Mass of LLbar System", 500, 0, 1000 ));


  // Jets 

  name = "HypJetMulti";
  m_histogram[name] = store(new TH1D ( name+step, "Jet Multiplicity", 10, -0.5, 9.5 ));
  name = "HypJetHT";
  m_histogram[name] = store(new TH1D ( name+step, "Jet HT", 800, 0, 800));
  name = "HypJetpT";
  m_histogram[name] = store(new TH1D ( name+step, "Jet pT", 80, 0, 400 ));
  name = "HypJetEta";
  m_histogram[name] = store(new TH1D ( name+step, "Jet Eta", 200, -5, 5 )); 
  name = "HypJetRapidity";
  m_histogram[name] = store(new TH1D ( name+step, "Jet Rapidity", 100, -5, 5 ));

  // Extra Jets 

  name = "HypExtraJetMulti";
  m_histogram[name] = store(new TH1D ( name+step, "Extra Jet Multiplicity", 10, -0.5, 9.5 ));
  name = "HypExtraJetpT";
  m_histogram[name] = store(new TH1D ( name+step, "Extra Jet pT", 80, 0, 400 ));
  name = "HypExtraJetHT";
  m_histogram[name] = store(new TH1D ( name+step, "Extra Jet HT", 800, 0, 800));
  name = "HypExtraJetEta";
  m_histogram[name] = store(new TH1D ( name+step, "Extra Jet Eta", 200, -5, 5 ));
  name = "HypExtraJetRapidity";
  m_histogram[name] = store(new TH1D ( name+step, "Extra Jet Rapidity", 100, -5, 5 ));

  // B Jets

  name = "HypBjetMulti";
  m_histogram[name] = store(new TH1D ( name+step, "B-Jet Multiplicity", 10, -0.5, 9.5 ));
  name = "HypBjetpT";
  m_histogram[name] = store(new TH1D ( name+step, "B-Jet pT", 1200, 0, 1200 ));
  name = "HypBjetHT";
  m_histogram[name] = store(new TH1D ( name+step, "B-Jet HT", 800, 0, 800));
  name = "HypBjetEta";
  m_histogram[name] = store(new TH1D ( name+step, "B-Jet Eta", 200, -5, 5 ));
  name = "HypBjetRapidity";
  m_histogram[name] = store(new TH1D ( name+step, "B-Jet Rapidity", 100, -5, 5 ));
  }
  //Declaring some essential histograms in the extra-jet phase space:end

  // Jets 

  // name = "HypJetMulti";
  // m_histogram[name] = store(new TH1D ( name+step, "Jet Multiplicity", 10, -0.5, 9.5 ));
  // name = "HypJetHT";
  // m_histogram[name] = store(new TH1D ( name+step, "Jet HT", 800, 0, 800));
  // name = "HypJetpT";
  // m_histogram[name] = store(new TH1D ( name+step, "Jet pT", 80, 0, 400 ));
  //  // name = "HypJetEta";

  // Extra Jets 

  //  name = "HypExtraJetMulti";
  //  m_histogram[name] = store(new TH1D ( name+step, "Extra Jet Multiplicity", 10, -0.5, 9.5 ));
  //  name = "HypExtraJetpT";
  //  m_histogram[name] = store(new TH1D ( name+step, "Extra Jet pT", 80, 0, 400 ));
  //  name = "HypExtraJetHT";
  //  m_histogram[name] = store(new TH1D ( name+step, "Extra Jet HT", 800, 0, 800));

  // B Jets

  // name = "HypBjetMulti";
  // m_histogram[name] = store(new TH1D ( name+step, "B-Jet Multiplicity", 10, -0.5, 9.5 ));
  // name = "HypBJetpT";
  // m_histogram[name] = store(new TH1D ( name+step, "B Hypothesis pT", 1200, 0, 1200 ));
  // name = "HypBJetHT";
  // m_histogram[name] = store(new TH1D ( name+step, "B-Jet HT", 800, 0, 800));
  // name = "HypBJetEta";
  // m_histogram[name] = store(new TH1D ( name+step, "B Hypothesis Eta", 200, -5, 5 ));
  // name = "HypBJetRapidity";
  // m_histogram[name] = store(new TH1D ( name+step, "B Hypothesis Eta", 100, -5, 5 ));
  // name = "HypBJetCSVdiscriminator";
  // m_histogram[name] = store(new TH1D ( name+step, "B Hypothesis CSV discriminator", 100, 0., 1. ));
  // name = "HypAntiBJetpT";
  // m_histogram[name] = store(new TH1D ( name+step, "AntiB Hypothesis pT", 1200, 0, 1200 ));
  // name = "HypAntiBJetEta";
  // m_histogram[name] = store(new TH1D ( name+step, "AntiB Hypothesis Eta", 200, -5, 5 ));
  // name = "HypAntiBJetRapidity";
  // m_histogram[name] = store(new TH1D ( name+step, "AntiB Hypothesis Eta", 100, -5, 5 ));
  // name = "HypAntiBJetCSVdiscriminator";
  // m_histogram[name] = store(new TH1D ( name+step, "AntiB Hypothesis CSV discriminator", 100, 0., 1. ));

  // name = "HypBBBarpT";
  // m_histogram[name] = store(new TH1D ( name+step, "p_{T} of bbbarpT", 400, 0, 800));
  // name = "HypBBBarMass";
  // m_histogram[name] = store(new TH1D( name+step, "Mass of bbbarMass", 400, 0, 800));
  // name = "HypBBBarDPhi";
  // m_histogram[name] = store(new TH1D ( name+step, "#Delta#phi(bbbar)",320, 0., 3.2 ));

  // name = "HypLeptonBjetMass";
  // m_histogram[name] = store(new TH1D ( name+step, "Mass(Lep, AntiBJet)", 500, 0, 1000 ));
  // name = "HypAntiLeptonBjetMass";
  // m_histogram[name] = store(new TH1D ( name+step, "Mass(AntiLep, BJet)", 500, 0, 1000 ));

  // MET and Neutrinos

  // name = "HypMET";
  // m_histogram[name] = store(new TH1D( name+step,"MET ", 500, 0, 500));

  // name = "HypNeutrinopT";
  // m_histogram[name] = store(new TH1D( name+step, "hyp nu pT", 80, 0, 400));
  // name = "HypAntiNeutrinopT";
  // m_histogram[name] = store(new TH1D( name+step, "hyp nubar pT", 80, 0, 400));
  // name = "HypNeutrinoEta";
  // m_histogram[name] = store(new TH1D ( name+step, "Neutrino Eta", 200,-5,5 ));
  // name = "HypAntiNeutrinoEta";
  // m_histogram[name] = store(new TH1D ( name+step, "AntiNeutrino Hypothesis Eta", 200,-5,5 ));

  // Top and Anti-Top

  // name = "HypToppT";
  // m_histogram[name] = store(new TH1D ( name+step, "Top pT",1200,0,1200 ));
  // name = "HypTopEta";
  // m_histogram[name] = store(new TH1D ( name+step, "Top #eta", 100, -5, 5 ));
  // name = "HypTopMass";
  // m_histogram[name] = store(new TH1D ( name+step, "Top Mass", 80, 0, 400 ));
  // name = "HypTopRapidity";
  // m_histogram[name] = store(new TH1D ( name+step, "Top Rapidity", 400, -5, 5 ));
  // name = "HypTopRapidityAbs";
  // m_histogram[name] = store(new TH1D ( name+step, "Absolute Top Rapidity", 200, 0, 5 ));
  // name = "HypToppT_TTBarFrame";
  // m_histogram[name] = store(new TH1D ( name+step, "Top pT TTBar Frame",1200,0,1200 ));
  // name = "HypTopEta_TTBarFrame";
  // m_histogram[name] = store(new TH1D ( name+step, "Top #eta TTBar Frame", 100, -5, 5 ));
  // name = "HypTopRapidity_TTBarFrame";
  // m_histogram[name] = store(new TH1D ( name+step, "Top Rapidity TTBar Frame", 400, -5, 5 ));

  // name = "HypAntiToppT";
  // m_histogram[name] = store(new TH1D ( name+step, "AntiTop pT",1200,0,1200 ));
  // name = "HypAntiTopEta";
  // m_histogram[name] = store(new TH1D ( name+step, "AntiTop #eta", 100, -5, 5 ));
  // name = "HypAntiTopMass";
  // m_histogram[name] = store(new TH1D ( name+step, "AntiTop Mass", 80, 0, 400 ));
  // name = "HypAntiTopRapidity";
  // m_histogram[name] = store(new TH1D ( name+step, "Top Rapidity", 400, -5, 5 ));
  // name = "HypAntiTopRapidityAbs";
  // m_histogram[name] = store(new TH1D ( name+step, "Absolute Top Rapidity", 200, 0, 5 ));
  // name = "HypAntiToppT_TTBarFrame";
  // m_histogram[name] = store(new TH1D ( name+step, "AntiTop pT TTBar Frame",1200,0,1200 ));
  // name = "HypAntiTopEta_TTBarFrame";
  // m_histogram[name] = store(new TH1D ( name+step, "AntiTop #eta TTBar Frame", 100, -5, 5 ));
  // name = "HypAntiTopRapidity_TTBarFrame";
  // m_histogram[name] = store(new TH1D ( name+step, "Top Rapidity TTBar Frame", 400, -5, 5 ));

  // name = "HypTTBarRapidity";
  // m_histogram[name] = store(new TH1D ( name+step, "Rapidity of TTbar System", 400, -5, 5 ));
  // name = "HypTTBarRapidityAbs";
  // m_histogram[name] = store(new TH1D ( name+step, "Absolute Rapidity of TTbar System", 200, 0, 5 ));
  // name = "HypTTBarpT";
  // m_histogram[name] = store(new TH1D ( name+step, "pT of TTbar System", 1200, 0, 1200 ));
  name = "HypTTBarMass";
  m_histogram[name] = store(new TH1D ( name+step, "Mass of TTbar System", 3000, 0, 3000 ));
  // name = "HypTTBarDPhi";
  // m_histogram[name] = store(new TH1D( name+step, "#Delta#Phi ofTTBar", 3500, 0, 3.5));
  // name = "HypTTBarDeltaRapidity";
  // m_histogram[name] = store(new TH1D( name+step, "|y^{t}| - |y^{#bar{t}}|", 400, 0, 5));

  name = "HypScatteringAngle_LabFrame";
  m_histogram[name] = store(new TH1D ( name+step, "ScatteringAngle LabFrame",200, -1.0, 1.0 ));
  name = "HypScatteringAngle_TTBarFrame";
  m_histogram[name] = store(new TH1D ( name+step, "ScatteringAngle TTBarFrame",200, -1.0, 1.0 ));

  // name = "HypMT2";
  // m_histogram[name] = this->store(new TH1D (name+step, "M_{T2}",200, 0.0, 200.0 ));

  // Kinematic Distributions End

  //UnSymmetrized Begin
  
  name = "HypAntiLeptonBk";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(k)",200, -1.0, 1.0 ));
  name = "HypLeptonBk";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{2}(k)",200, -1.0, 1.0 ));
  name = "HypAntiLeptonBj";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(k*)",200, -1.0, 1.0 ));
  name = "HypLeptonBj";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{2}(k*)",200, -1.0, 1.0 ));
  name = "HypAntiLeptonBr";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(r)",200, -1.0, 1.0 ));
  name = "HypLeptonBr";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{2}(r)",200, -1.0, 1.0 ));
  name = "HypAntiLeptonBq";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(r*)",200, -1.0, 1.0 ));
  name = "HypLeptonBq";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{2}(r*)",200, -1.0, 1.0 ));
  name = "HypAntiLeptonBn";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(n)",200, -1.0, 1.0 ));
  name = "HypLeptonBn";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{2}(n)",200, -1.0, 1.0 ));

  //name = "HypLLBarBPnn";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(n) + B_{2}(n)",400, -2.0, 2.0 ));
  //name = "HypLLBarBMnn";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(n) - B_{2}(n)",400, -2.0, 2.0 ));
  //name = "HypLLBarBPrr";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(r) + B_{2}(r)",400, -2.0, 2.0 ));
  //name = "HypLLBarBMrr";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(r) - B_{2}(r)",400, -2.0, 2.0 ));
  //name = "HypLLBarBPkk";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(k) + B_{2}(k)",400, -2.0, 2.0 ));
  //name = "HypLLBarBMkk";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(k) - B_{2}(k)",400, -2.0, 2.0 ));
  //name = "HypLLBarBPjj";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(k*) + B_{2}(k*)",400, -2.0, 2.0 ));
  //name = "HypLLBarBMjj";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(k*) - B_{2}(k*)",400, -2.0, 2.0 ));
  //name = "HypLLBarBPqq";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(r*) + B_{2}(r*)",400, -2.0, 2.0 ));
  //name = "HypLLBarBMqq";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(r*) - B_{2}(r*)",400, -2.0, 2.0 ));

  name = "HypLLBarCkk";
  m_histogram[name] = this->store(new TH1D (name+step, "C(k,k)",200, -1.0, 1.0 ));
  name = "HypLLBarCrr";
  m_histogram[name] = this->store(new TH1D (name+step, "C(r,r)",200, -1.0, 1.0 ));
  name = "HypLLBarCnn";
  m_histogram[name] = this->store(new TH1D (name+step, "C(n,n)",200, -1.0, 1.0 ));
  name = "HypLLBarCkj";
  m_histogram[name] = this->store(new TH1D (name+step, "C(k,k*)",200, -1.0, 1.0 ));
  name = "HypLLBarCrq";
  m_histogram[name] = this->store(new TH1D (name+step, "C(r,r*)",200, -1.0, 1.0 ));

  name = "HypLLBarCrk";
  m_histogram[name] = this->store(new TH1D (name+step, "C(r,k)",200, -1.0, 1.0 ));
  name = "HypLLBarCkr";
  m_histogram[name] = this->store(new TH1D (name+step, "C(k,r)",200, -1.0, 1.0 ));
  name = "HypLLBarCnr";
  m_histogram[name] = this->store(new TH1D (name+step, "C(n,r)",200, -1.0, 1.0 ));
  name = "HypLLBarCrn";
  m_histogram[name] = this->store(new TH1D (name+step, "C(r,n)",200, -1.0, 1.0 ));
  name = "HypLLBarCnk";
  m_histogram[name] = this->store(new TH1D (name+step, "C(n,k)",200, -1.0, 1.0 ));
  name = "HypLLBarCkn";
  m_histogram[name] = this->store(new TH1D (name+step, "C(k,n)",200, -1.0, 1.0 ));

  name = "HypLLBarCrj";
  m_histogram[name] = this->store(new TH1D (name+step, "C(r,k*)",200, -1.0, 1.0 ));
  name = "HypLLBarCjr";
  m_histogram[name] = this->store(new TH1D (name+step, "C(k*,r)",200, -1.0, 1.0 ));
  //name = "HypLLBarCqk";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(r*,k)",200, -1.0, 1.0 ));
  //name = "HypLLBarCkq";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(k,r*)",200, -1.0, 1.0 ));

  //name = "HypLLBarCqj";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(r*,k*)",200, -1.0, 1.0 ));
  //name = "HypLLBarCjq";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(k*,r*)",200, -1.0, 1.0 ));
  //name = "HypLLBarCnq";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(n,r*)",200, -1.0, 1.0 ));
  //name = "HypLLBarCqn";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(r*,n)",200, -1.0, 1.0 ));
  //name = "HypLLBarCnj";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(n,k*)",200, -1.0, 1.0 ));
  //name = "HypLLBarCjn";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(k*,n)",200, -1.0, 1.0 ));

  name = "HypLLBarCPrk";
  m_histogram[name] = this->store(new TH1D (name+step, "C(r,k) + C(k,r)",400, -2.0, 2.0 ));
  name = "HypLLBarCMrk";
  m_histogram[name] = this->store(new TH1D (name+step, "C(r,k) - C(k,r)",400, -2.0, 2.0 ));
  name = "HypLLBarCPnr";
  m_histogram[name] = this->store(new TH1D (name+step, "C(n,r) + C(r,n)",400, -2.0, 2.0 ));
  name = "HypLLBarCMnr";
  m_histogram[name] = this->store(new TH1D (name+step, "C(n,r) - C(r,n)",400, -2.0, 2.0 ));
  name = "HypLLBarCPnk";
  m_histogram[name] = this->store(new TH1D (name+step, "C(n,k) + C(k,n)",400, -2.0, 2.0 ));
  name = "HypLLBarCMnk";
  m_histogram[name] = this->store(new TH1D (name+step, "C(n,k) - C(k,n)",400, -2.0, 2.0 ));

  name = "HypLLBarCPrj";
  m_histogram[name] = this->store(new TH1D (name+step, "C(r,k*) + C(k*,r)",400, -2.0, 2.0 ));
  name = "HypLLBarCMrj";
  m_histogram[name] = this->store(new TH1D (name+step, "C(r,k*) - C(k*,r)",400, -2.0, 2.0 ));
  //name = "HypLLBarCPqj";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(r*,k*) + C(k*,r*)",400, -2.0, 2.0 ));
  //name = "HypLLBarCMqj";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(r*,k*) - C(k*,r*)",400, -2.0, 2.0 ));
  //name = "HypLLBarCPnq";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(n,r*) + C(r*,n)",400, -2.0, 2.0 ));
  //name = "HypLLBarCMnq";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(n,r*) - C(r*,n)",400, -2.0, 2.0 ));
  //name = "HypLLBarCPnj";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(n,k*) + C(k*,n)",400, -2.0, 2.0 ));
  //name = "HypLLBarCMnj";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(n,k*) - C(k*,n)",400, -2.0, 2.0 ));

  name = "HypLLBarChan";
  m_histogram[name] = this->store(new TH1D (name+step, "+C(k,k) - C(r,r) - C(n,n)",400, -2.0, 2.0 ));
  name = "HypLLBarCsca";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(k,k) + C(r,r) - C(n,n)",400, -2.0, 2.0 ));
  name = "HypLLBarCtra";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(k,k) - C(r,r) + C(n,n)",400, -2.0, 2.0 ));
  name = "HypLLBarCkjL";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(k,k*) - C(r,r) - C(n,n)",400, -2.0, 2.0 ));
  name = "HypLLBarCrqL";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(k,k) - C(r,r*) - C(n,n)",400, -2.0, 2.0 ));

  name = "HypLLBarCrkP";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(r,k) - C(k,r) - C(n,n)",400, -3.0, 3.0 ));
  name = "HypLLBarCrkM";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(r,k) + C(k,r) - C(n,n)",400, -3.0, 3.0 ));
  name = "HypLLBarCnrP";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(n,r) - C(r,n) - C(k,k)",400, -3.0, 3.0 ));
  name = "HypLLBarCnrM";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(n,r) + C(r,n) - C(k,k)",400, -3.0, 3.0 ));
  name = "HypLLBarCnkP";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(n,k) - C(k,n) - C(r,r)",400, -3.0, 3.0 ));
  name = "HypLLBarCnkM";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(n,k) + C(k,n) - C(r,r)",400, -3.0, 3.0 ));

  //name = "HypLLBarCqjP";
  //m_histogram[name] = this->store(new TH1D (name+step, "-C(q,j) - C(j,q) - C(n,n)",400, -3.0, 3.0 ));
  //name = "HypLLBarCqjM";
  //m_histogram[name] = this->store(new TH1D (name+step, "-C(r,j) + C(k,q) - C(n,n)",400, -3.0, 3.0 ));
  //name = "HypLLBarCnqP";
  //m_histogram[name] = this->store(new TH1D (name+step, "-C(n,q) - C(q,n) - C(j,j)",400, -3.0, 3.0 ));
  //name = "HypLLBarCnqM";
  //m_histogram[name] = this->store(new TH1D (name+step, "-C(n,q) + C(q,n) - C(j,j)",400, -3.0, 3.0 ));
  //name = "HypLLBarCnjP";
  //m_histogram[name] = this->store(new TH1D (name+step, "-C(n,j) - C(j,n) - C(q,q)",400, -3.0, 3.0 ));
  //name = "HypLLBarCnjM";
  //m_histogram[name] = this->store(new TH1D (name+step, "-C(n,j) + C(j,n) - C(q,q)",400, -3.0, 3.0 ));

  name = "HypLLBarcHel";
  m_histogram[name] = this->store(new TH1D (name+step, "cos(#theta_{LL})(Hel)",200, -1.0, 1.0 ));
  name = "HypLLBarcLab";
  m_histogram[name] = this->store(new TH1D (name+step, "cos(#theta_{LL})(Lab)",200, -1.0, 1.0 ));

  name = "HypLLBarkNorm";
  m_histogram[name] = this->store(new TH1D (name+step, "sin(#theta_{LL})(k)",200, -1.0, 1.0 ));
  name = "HypLLBarrNorm";
  m_histogram[name] = this->store(new TH1D (name+step, "sin(#theta_{LL})(r)",200, -1.0, 1.0 ));
  name = "HypLLBarnNorm";
  m_histogram[name] = this->store(new TH1D (name+step, "sin(#theta_{LL})(n)",200, -1.0, 1.0 ));
  name = "HypLLBarjNorm";
  m_histogram[name] = this->store(new TH1D (name+step, "sin(#theta_{LL})(k*)",200, -1.0, 1.0 ));
  name = "HypLLBarqNorm";
  m_histogram[name] = this->store(new TH1D (name+step, "sin(#theta_{LL})(r*)",200, -1.0, 1.0 ));

  
  //UnSymmetrized End


  //Symmetrized Begin

  name = "HypSymAntiLeptonBk";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(k)",200, -1.0, 1.0 ));
  name = "HypSymLeptonBk";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{2}(k)",200, -1.0, 1.0 ));
  name = "HypSymAntiLeptonBj";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(k*)",200, -1.0, 1.0 ));
  name = "HypSymLeptonBj";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{2}(k*)",200, -1.0, 1.0 ));
  name = "HypSymAntiLeptonBr";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(r)",200, -1.0, 1.0 ));
  name = "HypSymLeptonBr";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{2}(r)",200, -1.0, 1.0 ));
  name = "HypSymAntiLeptonBq";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(r*)",200, -1.0, 1.0 ));
  name = "HypSymLeptonBq";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{2}(r*)",200, -1.0, 1.0 ));
  name = "HypSymAntiLeptonBn";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(n)",200, -1.0, 1.0 ));
  name = "HypSymLeptonBn";
  m_histogram[name] = this->store(new TH1D (name+step, "B_{2}(n)",200, -1.0, 1.0 ));

  //name = "HypSymLLBarBPnn";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(n) + B_{2}(n)",400, -2.0, 2.0 ));
  //name = "HypSymLLBarBMnn";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(n) - B_{2}(n)",400, -2.0, 2.0 ));
  //name = "HypSymLLBarBPrr";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(r) + B_{2}(r)",400, -2.0, 2.0 ));
  //name = "HypSymLLBarBMrr";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(r) - B_{2}(r)",400, -2.0, 2.0 ));
  //name = "HypSymLLBarBPkk";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(k) + B_{2}(k)",400, -2.0, 2.0 ));
  //name = "HypSymLLBarBMkk";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(k) - B_{2}(k)",400, -2.0, 2.0 ));
  //name = "HypSymLLBarBPjj";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(k*) + B_{2}(k*)",400, -2.0, 2.0 ));
  //name = "HypSymLLBarBMjj";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(k*) - B_{2}(k*)",400, -2.0, 2.0 ));
  //name = "HypSymLLBarBPqq";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(r*) + B_{2}(r*)",400, -2.0, 2.0 ));
  //name = "HypSymLLBarBMqq";
  //m_histogram[name] = this->store(new TH1D (name+step, "B_{1}(r*) - B_{2}(r*)",400, -2.0, 2.0 ));

  name = "HypSymLLBarCkk";
  m_histogram[name] = this->store(new TH1D (name+step, "C(k,k)",200, -1.0, 1.0 ));
  name = "HypSymLLBarCrr";
  m_histogram[name] = this->store(new TH1D (name+step, "C(r,r)",200, -1.0, 1.0 ));
  name = "HypSymLLBarCnn";
  m_histogram[name] = this->store(new TH1D (name+step, "C(n,n)",200, -1.0, 1.0 ));
  name = "HypSymLLBarCkj";
  m_histogram[name] = this->store(new TH1D (name+step, "C(k,k*)",200, -1.0, 1.0 ));
  name = "HypSymLLBarCrq";
  m_histogram[name] = this->store(new TH1D (name+step, "C(r,r*)",200, -1.0, 1.0 ));

  name = "HypSymLLBarCrk";
  m_histogram[name] = this->store(new TH1D (name+step, "C(r,k)",200, -1.0, 1.0 ));
  name = "HypSymLLBarCkr";
  m_histogram[name] = this->store(new TH1D (name+step, "C(k,r)",200, -1.0, 1.0 ));
  name = "HypSymLLBarCnr";
  m_histogram[name] = this->store(new TH1D (name+step, "C(n,r)",200, -1.0, 1.0 ));
  name = "HypSymLLBarCrn";
  m_histogram[name] = this->store(new TH1D (name+step, "C(r,n)",200, -1.0, 1.0 ));
  name = "HypSymLLBarCnk";
  m_histogram[name] = this->store(new TH1D (name+step, "C(n,k)",200, -1.0, 1.0 ));
  name = "HypSymLLBarCkn";
  m_histogram[name] = this->store(new TH1D (name+step, "C(k,n)",200, -1.0, 1.0 ));

  name = "HypSymLLBarCrj";
  m_histogram[name] = this->store(new TH1D (name+step, "C(r,k*)",200, -1.0, 1.0 ));
  name = "HypSymLLBarCjr";
  m_histogram[name] = this->store(new TH1D (name+step, "C(k*,r)",200, -1.0, 1.0 ));
  //name = "HypSymLLBarCqk";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(r*,k)",200, -1.0, 1.0 ));
  //name = "HypSymLLBarCkq";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(k,r*)",200, -1.0, 1.0 ));

  //name = "HypSymLLBarCqj";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(r*,k*)",200, -1.0, 1.0 ));
  //name = "HypSymLLBarCjq";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(k*,r*)",200, -1.0, 1.0 ));
  //name = "HypSymLLBarCnq";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(n,r*)",200, -1.0, 1.0 ));
  //name = "HypSymLLBarCqn";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(r*,n)",200, -1.0, 1.0 ));
  //name = "HypSymLLBarCnj";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(n,k*)",200, -1.0, 1.0 ));
  //name = "HypSymLLBarCjn";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(k*,n)",200, -1.0, 1.0 ));

  name = "HypSymLLBarCPrk";
  m_histogram[name] = this->store(new TH1D (name+step, "C(r,k) + C(k,r)",400, -2.0, 2.0 ));
  name = "HypSymLLBarCMrk";
  m_histogram[name] = this->store(new TH1D (name+step, "C(r,k) - C(k,r)",400, -2.0, 2.0 ));
  name = "HypSymLLBarCPnr";
  m_histogram[name] = this->store(new TH1D (name+step, "C(n,r) + C(r,n)",400, -2.0, 2.0 ));
  name = "HypSymLLBarCMnr";
  m_histogram[name] = this->store(new TH1D (name+step, "C(n,r) - C(r,n)",400, -2.0, 2.0 ));
  name = "HypSymLLBarCPnk";
  m_histogram[name] = this->store(new TH1D (name+step, "C(n,k) + C(k,n)",400, -2.0, 2.0 ));
  name = "HypSymLLBarCMnk";
  m_histogram[name] = this->store(new TH1D (name+step, "C(n,k) - C(k,n)",400, -2.0, 2.0 ));

  name = "HypSymLLBarCPrj";
  m_histogram[name] = this->store(new TH1D (name+step, "C(r,k*) + C(k*,r)",400, -2.0, 2.0 ));
  name = "HypSymLLBarCMrj";
  m_histogram[name] = this->store(new TH1D (name+step, "C(r,k*) - C(k*,r)",400, -2.0, 2.0 ));  //name = "HypSymLLBarCPqj";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(r*,k*) + C(k*,r*)",400, -2.0, 2.0 ));
  //name = "HypSymLLBarCMqj";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(r*,k*) - C(k*,r*)",400, -2.0, 2.0 ));
  //name = "HypSymLLBarCPnq";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(n,r*) + C(r*,n)",400, -2.0, 2.0 ));
  //name = "HypSymLLBarCMnq";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(n,r*) - C(r*,n)",400, -2.0, 2.0 ));
  //name = "HypSymLLBarCPnj";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(n,k*) + C(k*,n)",400, -2.0, 2.0 ));
  //name = "HypSymLLBarCMnj";
  //m_histogram[name] = this->store(new TH1D (name+step, "C(n,k*) - C(k*,n)",400, -2.0, 2.0 ));

  name = "HypSymLLBarChan";
  m_histogram[name] = this->store(new TH1D (name+step, "+C(k,k) - C(r,r) - C(n,n)",400, -2.0, 2.0 ));
  name = "HypSymLLBarCsca";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(k,k) + C(r,r) - C(n,n)",400, -2.0, 2.0 ));
  name = "HypSymLLBarCtra";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(k,k) - C(r,r) + C(n,n)",400, -2.0, 2.0 ));
  name = "HypSymLLBarCkjL";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(k,k*) - C(r,r) - C(n,n)",400, -2.0, 2.0 ));
  name = "HypSymLLBarCrqL";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(k,k) - C(r,r*) - C(n,n)",400, -2.0, 2.0 ));

  name = "HypSymLLBarCrkP";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(r,k) - C(k,r) - C(n,n)",400, -3.0, 3.0 ));
  name = "HypSymLLBarCrkM";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(r,k) + C(k,r) - C(n,n)",400, -3.0, 3.0 ));
  name = "HypSymLLBarCnrP";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(n,r) - C(r,n) - C(k,k)",400, -3.0, 3.0 ));
  name = "HypSymLLBarCnrM";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(n,r) + C(r,n) - C(k,k)",400, -3.0, 3.0 ));
  name = "HypSymLLBarCnkP";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(n,k) - C(k,n) - C(r,r)",400, -3.0, 3.0 ));
  name = "HypSymLLBarCnkM";
  m_histogram[name] = this->store(new TH1D (name+step, "-C(n,k) + C(k,n) - C(r,r)",400, -3.0, 3.0 ));

  //name = "HypSymLLBarCqjP";
  //m_histogram[name] = this->store(new TH1D (name+step, "-C(q,j) - C(j,q) - C(n,n)",400, -3.0, 3.0 ));
  //name = "HypSymLLBarCqjM";
  //m_histogram[name] = this->store(new TH1D (name+step, "-C(r,j) + C(k,q) - C(n,n)",400, -3.0, 3.0 ));
  //name = "HypSymLLBarCnqP";
  //m_histogram[name] = this->store(new TH1D (name+step, "-C(n,q) - C(q,n) - C(j,j)",400, -3.0, 3.0 ));
  //name = "HypSymLLBarCnqM";
  //m_histogram[name] = this->store(new TH1D (name+step, "-C(n,q) + C(q,n) - C(j,j)",400, -3.0, 3.0 ));
  //name = "HypSymLLBarCnjP";
  //m_histogram[name] = this->store(new TH1D (name+step, "-C(n,j) - C(j,n) - C(q,q)",400, -3.0, 3.0 ));
  //name = "HypSymLLBarCnjM";
  //m_histogram[name] = this->store(new TH1D (name+step, "-C(n,j) + C(j,n) - C(q,q)",400, -3.0, 3.0 ));

  name = "HypSymLLBarcHel";
  m_histogram[name] = this->store(new TH1D (name+step, "cos(#theta_{LL})(Hel)",200, -1.0, 1.0 ));
  name = "HypSymLLBarcLab";
  m_histogram[name] = this->store(new TH1D (name+step, "cos(#theta_{LL})(Lab)",200, -1.0, 1.0 ));

  name = "HypSymLLBarkNorm";
  m_histogram[name] = this->store(new TH1D (name+step, "sin(#theta_{LL})(k)",200, -1.0, 1.0 ));
  name = "HypSymLLBarrNorm";
  m_histogram[name] = this->store(new TH1D (name+step, "sin(#theta_{LL})(r)",200, -1.0, 1.0 ));
  name = "HypSymLLBarnNorm";
  m_histogram[name] = this->store(new TH1D (name+step, "sin(#theta_{LL})(n)",200, -1.0, 1.0 ));
  name = "HypSymLLBarjNorm";
  m_histogram[name] = this->store(new TH1D (name+step, "sin(#theta_{LL})(k*)",200, -1.0, 1.0 ));
  name = "HypSymLLBarqNorm";
  m_histogram[name] = this->store(new TH1D (name+step, "sin(#theta_{LL})(r*)",200, -1.0, 1.0 ));
  
  //Symmetrized End


}



void AnalyzerSpinCorr::fillHistos(const EventMetadata& eventMetadata,
                                      const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                      const TopGenObjects& topGenObjects,
                                      const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
				      const LooseKinRecoSolution&,
                                      const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                                      const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
				      const double& weight, const double&, const TString& step ,
                                      std::map< TString, TH1* >& m_histogram)
{ 

    /// Reference to the analysis config
    //const AnalysisConfig& analysisConfig_;

    // Access object selections from config
    const AnalysisConfig::Selections& selections = analysisConfig_.selections();

    // Use utilities without namespaces
    using ROOT::Math::VectorUtil::DeltaPhi;
    using ROOT::Math::VectorUtil::DeltaR;
    using ROOT::Math::VectorUtil::Angle;
    using namespace common;

    //spin correlation  analysis
    const KinematicReconstructionSolution solution = kinematicReconstructionSolutions.solution();
    const LV& hypLepton = solution.lepton();
    const LV& hypAntiLepton = solution.antiLepton();
    LV hypllbar(hypLepton + hypAntiLepton);
    const LV& hypTop = solution.top();
    const LV& hypAntiTop = solution.antiTop();
    LV hypttbar(hypTop+hypAntiTop);
    const LV& hypBjet = solution.bjet();
    const LV& hypAntiBjet = solution.antiBjet();
    LV hypBBbar(hypBjet+hypAntiBjet);
    const LV& hypNeutrino = solution.neutrino();
    const LV& hypAntiNeutrino = solution.antiNeutrino();
    LV NuNubar(hypNeutrino+hypAntiNeutrino);
    
    const LV& recomet = *recoObjects.met_;

    double hypMT2 = AnalyzerSpinCorr::getMT2Variable(hypLepton, hypAntiLepton, recomet);

    //Ajeeta: Spin correlation variables from JHEP 12(2015)026
    VariablesPhiTT *vars = VariablesPhiTT::fillVariables(eventMetadata, recoObjects, commonGenObjects, topGenObjects, kinematicReconstructionSolutions, recoObjectIndices, genObjectIndices, genLevelWeights, recoLevelWeights, weight);

    TLorentzVector T_hyptop(common::LVtoTLV(solution.top()));
    TLorentzVector T_hypantitop(common::LVtoTLV(solution.antiTop()));
    TLorentzVector T_hypttbar(T_hyptop+T_hypantitop);
    TLorentzVector T_hyplep(common::LVtoTLV(solution.lepton()));
    TLorentzVector T_hypantilep(common::LVtoTLV(solution.antiLepton()));
    
    TVector3 TTBarFrameBoost(-1. * T_hypttbar.BoostVector());
    TVector3 TopFrameBoost(-1. * T_hyptop.BoostVector());
    TVector3 AntiTopFrameBoost(-1. * T_hypantitop.BoostVector());

    TLorentzVector hyptop_ttbarframe = T_hyptop;
    hyptop_ttbarframe.Boost(TTBarFrameBoost);
    TLorentzVector hypantitop_ttbarframe = T_hypantitop;
    hypantitop_ttbarframe.Boost(TTBarFrameBoost);
    TLorentzVector hyplep_antitopframe = T_hyplep;
    hyplep_antitopframe.Boost(AntiTopFrameBoost);
    TLorentzVector hypantilep_topframe = T_hypantilep;
    hypantilep_topframe.Boost(TopFrameBoost);

    TVector3 ProtonBeamDirection(0., 0., 1.);

    //    const int n_ExtraJets =  recoObjectIndices.n_ExtraJets_;
    const std::vector<int> JetIndices = recoObjectIndices.jetIndices_;
    const std::vector<int> BJetIndices = recoObjectIndices.bjetIndices_;
    //    const std::vector<int> ExtraJetIndices = recoObjectIndices.ExtraJetIndices_;
    //    const std::vector<int> ExtraJetIndicesIso04 = recoObjectIndices.ExtraJetIndicesIso04_;
    const std::vector<int> ExtraJetIndicesIso08 = recoObjectIndices.extraJetIndicesIso08_;
      // std::cout<< step << " Number of Jets: " << JetIndices.size() << "\n";
      //  std::cout<< step << " Number of Extra Jets: " << ExtraJetIndicesIso08.size() << "\n";


    //    if(ExtraJetIndicesIso08.size() > 0) {

    //    }
    

    // Kinematic Begin


    // m_histogram["HypvertMulti"]->Fill(recoObjects.vertMulti_, weight);

    // m_histogram["HypLeptonpT"]->Fill( hypLepton.pt(), weight);
    // m_histogram["HypAntiLeptonpT"]->Fill( hypAntiLepton.pt(), weight);
    // m_histogram["HypLeptonEta"]->Fill( hypLepton.eta(), weight);
    // m_histogram["HypAntiLeptonEta"]->Fill( hypAntiLepton.eta(), weight);

    // m_histogram["HypLeptonpT_AntiTopFrame"]->Fill( hyplep_antitopframe.Pt(), weight);
    // m_histogram["HypLeptonEta_AntiTopFrame"]->Fill( hyplep_antitopframe.Eta(), weight);
    // m_histogram["HypAntiLeptonpT_TopFrame"]->Fill( hypantilep_topframe.Pt(), weight);
    // m_histogram["HypAntiLeptonEta_TopFrame"]->Fill( hypantilep_topframe.Eta(), weight);
    if(step.Contains("extrajets")){
      m_histogram["HypLLBarMass"]->Fill( hypllbar.mass(), weight);
      m_histogram["HypJetMulti"]->Fill( JetIndices.size(), weight);
      m_histogram["HypJetHT"]->Fill( AnalyzerSpinCorr::getJetHT(JetIndices, (*recoObjects.jets_)), weight);
      for(int index: JetIndices) {
          const LV& jet = recoObjects.jets_->at(index);
          m_histogram["HypJetpT"]->Fill( jet.Pt(), weight);
          m_histogram["HypJetEta"]->Fill( jet.Eta(), weight);
          m_histogram["HypJetRapidity"]->Fill( jet.Rapidity(), weight);
    }

      m_histogram["HypExtraJetMulti"]->Fill( ExtraJetIndicesIso08.size(), weight);
      m_histogram["HypExtraJetHT"]->Fill( AnalyzerSpinCorr::getJetHT(ExtraJetIndicesIso08, (*recoObjects.jets_)), weight);
      for(int index: ExtraJetIndicesIso08) {
          const LV& jet = recoObjects.jets_->at(index);
          m_histogram["HypExtraJetpT"]->Fill( jet.Pt(), weight);
          m_histogram["HypExtraJetEta"]->Fill( jet.Eta(), weight);
          m_histogram["HypExtraJetRapidity"]->Fill( jet.Rapidity(), weight);
    }

      m_histogram["HypBjetMulti"]->Fill( BJetIndices.size(), weight);
      m_histogram["HypBjetHT"]->Fill( AnalyzerSpinCorr::getJetHT(BJetIndices, (*recoObjects.jets_)), weight);
      for(int index: BJetIndices) {
          const LV& jet = recoObjects.jets_->at(index);
          m_histogram["HypBjetpT"]->Fill( jet.Pt(), weight);
          m_histogram["HypBjetEta"]->Fill( jet.Eta(), weight);
          m_histogram["HypBjetRapidity"]->Fill( jet.Rapidity(), weight);
    }
    }
    // m_histogram["HypLLBarMass"]->Fill( hypllbar.mass(), weight);
    // m_histogram["HypLLBarpT"]->Fill( hypllbar.pt(), weight);
    m_histogram["HypLLBarDPhi"]->Fill(std::fabs(ROOT::Math::VectorUtil::DeltaPhi(hypLepton,hypAntiLepton)), weight);
    m_histogram["HypLLBarDEta"]->Fill(std::abs(hypLepton.Eta() - hypAntiLepton.Eta()) , weight);

    // m_histogram["HypJetMulti"]->Fill( JetIndices.size(), weight);
    // m_histogram["HypJetHT"]->Fill( AnalyzerSpinCorr::getJetHT(JetIndices, (*recoObjects.jets_)), weight);
    // for(int index: JetIndices) {
        // const LV& jet = recoObjects.jets_->at(index);
        // m_histogram["HypJetpT"]->Fill( jet.Pt(), weight);
	//	m_histogram["HypJetEta"]->Fill( jet.Eta(), weight);
    // }
    /*
    m_histogram["HypExtraJetMulti"]->Fill( ExtraJetIndicesIso08.size(), weight);
    m_histogram["HypExtraJetHT"]->Fill( AnalyzerSpinCorr::getJetHT(ExtraJetIndicesIso08, (*recoObjects.jets_)), weight);
    for(int index: ExtraJetIndicesIso08) {
	const LV& extrajet = recoObjects.jets_->at(index);
	m_histogram["HypExtraJetpT"]->Fill( extrajet.Pt(), weight);
	//	m_histogram["HypExtraJetEta"]->Fill( extrajet.Eta(), weight);
//	std::cout << "index: " << index << " Pt: " << extrajet.Pt() << "\n";
	}
    */
    // m_histogram["HypBjetMulti"]->Fill( BJetIndices.size(), weight);
    // m_histogram["HypBJetHT"]->Fill( AnalyzerSpinCorr::getJetHT(BJetIndices, (*recoObjects.jets_)), weight);
    // m_histogram["HypBJetpT"]->Fill( hypBjet.pt(), weight);
    // m_histogram["HypBJetEta"]->Fill( hypBjet.eta(), weight);
    // m_histogram["HypBJetRapidity"]->Fill( hypBjet.Rapidity(), weight);
    // m_histogram["HypBJetCSVdiscriminator"]->Fill( (recoObjects.jetBtags_)->at(solution.bjetIndex()), weight);
    // m_histogram["HypAntiBJetpT"]->Fill( hypAntiBjet.pt(), weight);
    // m_histogram["HypAntiBJetEta"]->Fill( hypAntiBjet.eta(), weight);
    // m_histogram["HypAntiBJetRapidity"]->Fill( hypAntiBjet.Rapidity(), weight);
    // m_histogram["HypAntiBJetCSVdiscriminator"]->Fill( (recoObjects.jetBtags_)->at(solution.antiBjetIndex()), weight);

    // m_histogram["HypBBBarpT"]->Fill( hypBBbar.pt(), weight);
    // m_histogram["HypBBBarMass"]->Fill( hypBBbar.mass(), weight);
    // m_histogram["HypBBBarDPhi"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(hypBjet, hypAntiBjet), weight);

    // m_histogram["HypLeptonBjetMass"]->Fill( (hypAntiBjet+hypLepton).mass(), weight);
    // m_histogram["HypAntiLeptonBjetMass"]->Fill( (hypBjet+hypAntiLepton).mass(), weight);

    // m_histogram["HypMET"]->Fill( recomet.pt(), weight);

    // m_histogram["HypNeutrinopT"]->Fill( hypNeutrino.pt(), weight);
    // m_histogram["HypAntiNeutrinopT"]->Fill( hypAntiNeutrino.pt(), weight);
    // m_histogram["HypNeutrinoEta"]->Fill( hypNeutrino.eta(), weight);
    // m_histogram["HypAntiNeutrinoEta"]->Fill( hypAntiNeutrino.eta(), weight);

    // m_histogram["HypToppT"]->Fill( hypTop.pt(), weight);
    // m_histogram["HypTopEta"]->Fill( hypTop.eta(), weight);
    // m_histogram["HypTopMass"]->Fill( hypTop.mass(), weight);
    // m_histogram["HypTopRapidity"]->Fill( hypTop.Rapidity(), weight);
    // m_histogram["HypTopRapidityAbs"]->Fill( std::abs(hypTop.Rapidity()), weight);
    // m_histogram["HypToppT_TTBarFrame"]->Fill( hyptop_ttbarframe.Pt(), weight);
    // m_histogram["HypTopEta_TTBarFrame"]->Fill( hyptop_ttbarframe.Eta(), weight);
    // m_histogram["HypTopRapidity_TTBarFrame"]->Fill( hyptop_ttbarframe.Rapidity(), weight);
    // m_histogram["HypAntiToppT"]->Fill( hypAntiTop.pt(), weight);
    // m_histogram["HypAntiTopEta"]->Fill( hypAntiTop.eta(), weight);
    // m_histogram["HypAntiTopMass"]->Fill( hypAntiTop.mass(), weight);
    // m_histogram["HypAntiTopRapidity"]->Fill( hypAntiTop.Rapidity(), weight);
    // m_histogram["HypAntiTopRapidityAbs"]->Fill( std::abs(hypAntiTop.Rapidity()), weight);
    // m_histogram["HypAntiToppT_TTBarFrame"]->Fill( hypantitop_ttbarframe.Pt(), weight);
    // m_histogram["HypAntiTopEta_TTBarFrame"]->Fill( hypantitop_ttbarframe.Eta(), weight);
    // m_histogram["HypAntiTopRapidity_TTBarFrame"]->Fill( hypantitop_ttbarframe.Rapidity(), weight);

    // m_histogram["HypTTBarRapidity"]->Fill( hypttbar.Rapidity(), weight);
    // m_histogram["HypTTBarRapidityAbs"]->Fill( std::abs(hypttbar.Rapidity()), weight);
    // m_histogram["HypTTBarpT"]->Fill( hypttbar.pt(), weight);
    m_histogram["HypTTBarMass"]->Fill( hypttbar.mass(), weight);
    // m_histogram["HypTTBarDPhi"]->Fill( ROOT::Math::VectorUtil::DeltaPhi(hypTop, hypAntiTop), weight);
    // m_histogram["HypTTBarDeltaRapidity"]->Fill( std::abs(hypTop.Rapidity() - hypAntiTop.Rapidity()), weight);

    m_histogram["HypScatteringAngle_LabFrame"]->Fill( T_hyptop.Vect().Unit().Dot(ProtonBeamDirection), weight);
    m_histogram["HypScatteringAngle_TTBarFrame"]->Fill( hyptop_ttbarframe.Vect().Unit().Dot(ProtonBeamDirection), weight);

    // m_histogram["HypMT2"]->Fill(hypMT2, weight);

    //m_histogram[""]->Fill( , weight);

    // Kinematic End

    //UnSymmetrized Begin
    
    m_histogram["HypAntiLeptonBk"]->Fill(vars->sol_b1k.value_,weight);
    m_histogram["HypLeptonBk"]->Fill(vars->sol_b2k.value_,weight);
    m_histogram["HypAntiLeptonBj"]->Fill(vars->sol_b1j.value_,weight);
    m_histogram["HypLeptonBj"]->Fill(vars->sol_b2j.value_,weight);
    m_histogram["HypAntiLeptonBr"]->Fill(vars->sol_b1r.value_,weight);
    m_histogram["HypLeptonBr"]->Fill(vars->sol_b2r.value_,weight);
    m_histogram["HypAntiLeptonBq"]->Fill(vars->sol_b1q.value_,weight);
    m_histogram["HypLeptonBq"]->Fill(vars->sol_b2q.value_,weight);
    m_histogram["HypAntiLeptonBn"]->Fill(vars->sol_b1n.value_,weight);
    m_histogram["HypLeptonBn"]->Fill(vars->sol_b2n.value_,weight);

    //m_histogram["HypLLBarBPnn"]->Fill(vars->sol_bP_nn.value_,weight);
    //m_histogram["HypLLBarBMnn"]->Fill(vars->sol_bM_nn.value_,weight);
    //m_histogram["HypLLBarBPrr"]->Fill(vars->sol_bP_rr.value_,weight);
    //m_histogram["HypLLBarBMrr"]->Fill(vars->sol_bM_rr.value_,weight);
    //m_histogram["HypLLBarBPkk"]->Fill(vars->sol_bP_kk.value_,weight);
    //m_histogram["HypLLBarBMkk"]->Fill(vars->sol_bM_kk.value_,weight);
    //m_histogram["HypLLBarBPjj"]->Fill(vars->sol_bP_jj.value_,weight);
    //m_histogram["HypLLBarBMjj"]->Fill(vars->sol_bM_jj.value_,weight);
    //m_histogram["HypLLBarBPqq"]->Fill(vars->sol_bP_qq.value_,weight);
    //m_histogram["HypLLBarBMqq"]->Fill(vars->sol_bM_qq.value_,weight);

    m_histogram["HypLLBarCkk"]->Fill(vars->sol_ckk.value_,weight);
    m_histogram["HypLLBarCrr"]->Fill(vars->sol_crr.value_,weight);
    m_histogram["HypLLBarCnn"]->Fill(vars->sol_cnn.value_,weight);
    m_histogram["HypLLBarCkj"]->Fill(vars->sol_ckj.value_,weight);
    m_histogram["HypLLBarCrq"]->Fill(vars->sol_crq.value_,weight);

    m_histogram["HypLLBarCrk"]->Fill(vars->sol_crk.value_,weight);
    m_histogram["HypLLBarCkr"]->Fill(vars->sol_ckr.value_,weight);
    m_histogram["HypLLBarCnr"]->Fill(vars->sol_cnr.value_,weight);
    m_histogram["HypLLBarCrn"]->Fill(vars->sol_crn.value_,weight);
    m_histogram["HypLLBarCnk"]->Fill(vars->sol_cnk.value_,weight);
    m_histogram["HypLLBarCkn"]->Fill(vars->sol_ckn.value_,weight);

    m_histogram["HypLLBarCrj"]->Fill(vars->sol_crj.value_,weight);
    m_histogram["HypLLBarCjr"]->Fill(vars->sol_cjr.value_,weight);
    //m_histogram["HypLLBarCqk"]->Fill(vars->sol_cqk.value_,weight);
    //m_histogram["HypLLBarCkq"]->Fill(vars->sol_ckq.value_,weight);

    //m_histogram["HypLLBarCqj"]->Fill(vars->sol_cqj.value_,weight);
    //m_histogram["HypLLBarCjq"]->Fill(vars->sol_cjq.value_,weight);
    //m_histogram["HypLLBarCnq"]->Fill(vars->sol_cnq.value_,weight);
    //m_histogram["HypLLBarCqn"]->Fill(vars->sol_cqn.value_,weight);
    //m_histogram["HypLLBarCnj"]->Fill(vars->sol_cnj.value_,weight);
    //m_histogram["HypLLBarCjn"]->Fill(vars->sol_cjn.value_,weight);

    m_histogram["HypLLBarCPrk"]->Fill(vars->sol_crk.value_ + vars->sol_ckr.value_,weight);
    m_histogram["HypLLBarCMrk"]->Fill(vars->sol_crk.value_ - vars->sol_ckr.value_,weight);
    m_histogram["HypLLBarCPnr"]->Fill(vars->sol_cnr.value_ + vars->sol_crn.value_,weight);
    m_histogram["HypLLBarCMnr"]->Fill(vars->sol_cnr.value_ - vars->sol_crn.value_,weight);
    m_histogram["HypLLBarCPnk"]->Fill(vars->sol_cnk.value_ + vars->sol_ckn.value_,weight);
    m_histogram["HypLLBarCMnk"]->Fill(vars->sol_cnk.value_ - vars->sol_ckn.value_,weight);

    m_histogram["HypLLBarCPrj"]->Fill(vars->sol_crj.value_ + vars->sol_cjr.value_,weight);
    m_histogram["HypLLBarCMrj"]->Fill(vars->sol_crj.value_ - vars->sol_cjr.value_,weight);
    //m_histogram["HypLLBarCPqj"]->Fill(vars->sol_cqj.value_ + vars->sol_cjq.value_,weight);
    //m_histogram["HypLLBarCMqj"]->Fill(vars->sol_cqj.value_ - vars->sol_cjq.value_,weight);
    //m_histogram["HypLLBarCPnq"]->Fill(vars->sol_cnq.value_ + vars->sol_cqn.value_,weight);
    //m_histogram["HypLLBarCMnq"]->Fill(vars->sol_cnq.value_ - vars->sol_cqn.value_,weight);
    //m_histogram["HypLLBarCPnj"]->Fill(vars->sol_cnj.value_ + vars->sol_cjn.value_,weight);
    //m_histogram["HypLLBarCMnj"]->Fill(vars->sol_cnj.value_ + vars->sol_cjn.value_,weight);

    m_histogram["HypLLBarcHel"]->Fill(vars->sol_cHel.value_,weight);
    m_histogram["HypLLBarcLab"]->Fill(vars->sol_cLab.value_,weight);

    m_histogram["HypLLBarkNorm"]->Fill(vars->sol_kNorm.value_,weight);
    m_histogram["HypLLBarrNorm"]->Fill(vars->sol_rNorm.value_,weight);
    m_histogram["HypLLBarnNorm"]->Fill(vars->sol_nNorm.value_,weight);
    m_histogram["HypLLBarjNorm"]->Fill(vars->sol_jNorm.value_,weight);
    m_histogram["HypLLBarqNorm"]->Fill(vars->sol_qNorm.value_,weight);

    m_histogram["HypLLBarChan"]->Fill(+vars->sol_ckk.value_ - vars->sol_crr.value_ - vars->sol_cnn.value_,weight);
    m_histogram["HypLLBarCsca"]->Fill(-vars->sol_ckk.value_ + vars->sol_crr.value_ - vars->sol_cnn.value_,weight);
    m_histogram["HypLLBarCtra"]->Fill(-vars->sol_ckk.value_ - vars->sol_crr.value_ + vars->sol_cnn.value_,weight);

    m_histogram["HypLLBarCkjL"]->Fill(-vars->sol_ckj.value_ - vars->sol_crr.value_ - vars->sol_cnn.value_,weight);
    m_histogram["HypLLBarCrqL"]->Fill(-vars->sol_ckk.value_ - vars->sol_crq.value_ - vars->sol_cnn.value_,weight);

    m_histogram["HypLLBarCrkP"]->Fill(-vars->sol_crk.value_ - vars->sol_ckr.value_ - vars->sol_cnn.value_,weight);
    m_histogram["HypLLBarCrkM"]->Fill(-vars->sol_crk.value_ + vars->sol_ckr.value_ - vars->sol_cnn.value_,weight);
    m_histogram["HypLLBarCnrP"]->Fill(-vars->sol_cnr.value_ - vars->sol_crn.value_ - vars->sol_ckk.value_,weight);
    m_histogram["HypLLBarCnrM"]->Fill(-vars->sol_cnr.value_ + vars->sol_crn.value_ - vars->sol_ckk.value_,weight);
    m_histogram["HypLLBarCnkP"]->Fill(-vars->sol_cnk.value_ - vars->sol_ckn.value_ - vars->sol_crr.value_,weight);
    m_histogram["HypLLBarCnkM"]->Fill(-vars->sol_cnk.value_ + vars->sol_ckn.value_ - vars->sol_crr.value_,weight);

    //m_histogram["HypLLBarCqjP"]->Fill(-vars->sol_cqj.value_ - vars->sol_cjq.value_ - vars->sol_cnn.value_,weight);
    //m_histogram["HypLLBarCqjM"]->Fill(-vars->sol_cqj.value_ + vars->sol_cjq.value_ - vars->sol_cnn.value_,weight);
    //m_histogram["HypLLBarCnqP"]->Fill(-vars->sol_cnq.value_ - vars->sol_cqn.value_ - vars->sol_ckj.value_,weight);
    //m_histogram["HypLLBarCnqM"]->Fill(-vars->sol_cnq.value_ + vars->sol_cqn.value_ - vars->sol_ckj.value_,weight);
    //m_histogram["HypLLBarCnjP"]->Fill(-vars->sol_cnj.value_ - vars->sol_cjn.value_ - vars->sol_crq.value_,weight);
    //m_histogram["HypLLBarCnjM"]->Fill(-vars->sol_cnj.value_ + vars->sol_cjn.value_ - vars->sol_crq.value_,weight);

    //UnSymmetrized End

    //Symmetrized Begin 

    //To symmetrize for blinding - fill twice for negative and positive values with 0.5 weight
    double symweight=weight*0.5;

    m_histogram["HypSymAntiLeptonBk"]->Fill(vars->sol_b1k.value_,symweight);
    m_histogram["HypSymLeptonBk"]->Fill(vars->sol_b2k.value_,symweight);
    m_histogram["HypSymAntiLeptonBj"]->Fill(vars->sol_b1j.value_,symweight);
    m_histogram["HypSymLeptonBj"]->Fill(vars->sol_b2j.value_,symweight);
    m_histogram["HypSymAntiLeptonBr"]->Fill(vars->sol_b1r.value_,symweight);
    m_histogram["HypSymLeptonBr"]->Fill(vars->sol_b2r.value_,symweight);
    m_histogram["HypSymAntiLeptonBq"]->Fill(vars->sol_b1q.value_,symweight);
    m_histogram["HypSymLeptonBq"]->Fill(vars->sol_b2q.value_,symweight);
    m_histogram["HypSymAntiLeptonBn"]->Fill(vars->sol_b1n.value_,symweight);
    m_histogram["HypSymLeptonBn"]->Fill(vars->sol_b2n.value_,symweight);

    //m_histogram["HypSymLLBarBPnn"]->Fill(vars->sol_bP_nn.value_,symweight);
    //m_histogram["HypSymLLBarBMnn"]->Fill(vars->sol_bM_nn.value_,symweight);
    //m_histogram["HypSymLLBarBPrr"]->Fill(vars->sol_bP_rr.value_,symweight);
    //m_histogram["HypSymLLBarBMrr"]->Fill(vars->sol_bM_rr.value_,symweight);
    //m_histogram["HypSymLLBarBPkk"]->Fill(vars->sol_bP_kk.value_,symweight);
    //m_histogram["HypSymLLBarBMkk"]->Fill(vars->sol_bM_kk.value_,symweight);
    //m_histogram["HypSymLLBarBPjj"]->Fill(vars->sol_bP_jj.value_,symweight);
    //m_histogram["HypSymLLBarBMjj"]->Fill(vars->sol_bM_jj.value_,symweight);
    //m_histogram["HypSymLLBarBPqq"]->Fill(vars->sol_bP_qq.value_,symweight);
    //m_histogram["HypSymLLBarBMqq"]->Fill(vars->sol_bM_qq.value_,symweight);

    m_histogram["HypSymLLBarCkk"]->Fill(vars->sol_ckk.value_,symweight);
    m_histogram["HypSymLLBarCrr"]->Fill(vars->sol_crr.value_,symweight);
    m_histogram["HypSymLLBarCnn"]->Fill(vars->sol_cnn.value_,symweight);
    m_histogram["HypSymLLBarCkj"]->Fill(vars->sol_ckj.value_,symweight);
    m_histogram["HypSymLLBarCrq"]->Fill(vars->sol_crq.value_,symweight);

    m_histogram["HypSymLLBarCrk"]->Fill(vars->sol_crk.value_,symweight);
    m_histogram["HypSymLLBarCkr"]->Fill(vars->sol_ckr.value_,symweight);
    m_histogram["HypSymLLBarCnr"]->Fill(vars->sol_cnr.value_,symweight);
    m_histogram["HypSymLLBarCrn"]->Fill(vars->sol_crn.value_,symweight);
    m_histogram["HypSymLLBarCnk"]->Fill(vars->sol_cnk.value_,symweight);
    m_histogram["HypSymLLBarCkn"]->Fill(vars->sol_ckn.value_,symweight);  

    m_histogram["HypSymLLBarCrj"]->Fill(vars->sol_crj.value_,symweight);
    m_histogram["HypSymLLBarCjr"]->Fill(vars->sol_cjr.value_,symweight);
    //m_histogram["HypSymLLBarCqk"]->Fill(vars->sol_cqk.value_,symweight);
    //m_histogram["HypSymLLBarCkq"]->Fill(vars->sol_ckq.value_,symweight);

    //m_histogram["HypSymLLBarCqj"]->Fill(vars->sol_cqj.value_,symweight);
    //m_histogram["HypSymLLBarCjq"]->Fill(vars->sol_cjq.value_,symweight);
    //m_histogram["HypSymLLBarCnq"]->Fill(vars->sol_cnq.value_,symweight);
    //m_histogram["HypSymLLBarCqn"]->Fill(vars->sol_cqn.value_,symweight);
    //m_histogram["HypSymLLBarCnj"]->Fill(vars->sol_cnj.value_,symweight);
    //m_histogram["HypSymLLBarCjn"]->Fill(vars->sol_cjn.value_,symweight);

    m_histogram["HypSymLLBarCPrk"]->Fill(vars->sol_crk.value_ + vars->sol_ckr.value_,symweight);
    m_histogram["HypSymLLBarCMrk"]->Fill(vars->sol_crk.value_ - vars->sol_ckr.value_,symweight);
    m_histogram["HypSymLLBarCPnr"]->Fill(vars->sol_cnr.value_ + vars->sol_crn.value_,symweight);
    m_histogram["HypSymLLBarCMnr"]->Fill(vars->sol_cnr.value_ - vars->sol_crn.value_,symweight);
    m_histogram["HypSymLLBarCPnk"]->Fill(vars->sol_cnk.value_ + vars->sol_ckn.value_,symweight);
    m_histogram["HypSymLLBarCMnk"]->Fill(vars->sol_cnk.value_ - vars->sol_ckn.value_,symweight);

    m_histogram["HypSymLLBarCPrj"]->Fill(vars->sol_crj.value_ + vars->sol_cjr.value_,symweight);
    m_histogram["HypSymLLBarCMrj"]->Fill(vars->sol_crj.value_ - vars->sol_cjr.value_,symweight);
    //m_histogram["HypSymLLBarCPqj"]->Fill(vars->sol_cqj.value_ + vars->sol_cjq.value_,symweight);
    //m_histogram["HypSymLLBarCMqj"]->Fill(vars->sol_cqj.value_ - vars->sol_cjq.value_,symweight);
    //m_histogram["HypSymLLBarCPnq"]->Fill(vars->sol_cnq.value_ + vars->sol_cqn.value_,symweight);
    //m_histogram["HypSymLLBarCMnq"]->Fill(vars->sol_cnq.value_ - vars->sol_cqn.value_,symweight);
    //m_histogram["HypSymLLBarCPnj"]->Fill(vars->sol_cnj.value_ + vars->sol_cjn.value_,symweight);
    //m_histogram["HypSymLLBarCMnj"]->Fill(vars->sol_cnj.value_ + vars->sol_cjn.value_,symweight);

    m_histogram["HypSymLLBarcHel"]->Fill(vars->sol_cHel.value_,symweight);
    m_histogram["HypSymLLBarcLab"]->Fill(vars->sol_cLab.value_,symweight);

    m_histogram["HypSymLLBarkNorm"]->Fill(vars->sol_kNorm.value_,symweight);
    m_histogram["HypSymLLBarrNorm"]->Fill(vars->sol_rNorm.value_,symweight);
    m_histogram["HypSymLLBarnNorm"]->Fill(vars->sol_nNorm.value_,symweight);
    m_histogram["HypSymLLBarjNorm"]->Fill(vars->sol_jNorm.value_,symweight);
    m_histogram["HypSymLLBarqNorm"]->Fill(vars->sol_qNorm.value_,symweight);

    m_histogram["HypSymLLBarChan"]->Fill(+vars->sol_ckk.value_ - vars->sol_crr.value_ - vars->sol_cnn.value_,symweight);
    m_histogram["HypSymLLBarCsca"]->Fill(-vars->sol_ckk.value_ + vars->sol_crr.value_ - vars->sol_cnn.value_,symweight);
    m_histogram["HypSymLLBarCtra"]->Fill(-vars->sol_ckk.value_ - vars->sol_crr.value_ + vars->sol_cnn.value_,symweight);

    m_histogram["HypSymLLBarCkjL"]->Fill(-vars->sol_ckj.value_ - vars->sol_crr.value_ - vars->sol_cnn.value_,symweight);
    m_histogram["HypSymLLBarCrqL"]->Fill(-vars->sol_ckk.value_ - vars->sol_crq.value_ - vars->sol_cnn.value_,symweight);

    m_histogram["HypSymLLBarCrkP"]->Fill(-vars->sol_crk.value_ - vars->sol_ckr.value_ - vars->sol_cnn.value_,symweight);
    m_histogram["HypSymLLBarCrkM"]->Fill(-vars->sol_crk.value_ + vars->sol_ckr.value_ - vars->sol_cnn.value_,symweight);
    m_histogram["HypSymLLBarCnrP"]->Fill(-vars->sol_cnr.value_ - vars->sol_crn.value_ - vars->sol_ckk.value_,symweight);
    m_histogram["HypSymLLBarCnrM"]->Fill(-vars->sol_cnr.value_ + vars->sol_crn.value_ - vars->sol_ckk.value_,symweight);
    m_histogram["HypSymLLBarCnkP"]->Fill(-vars->sol_cnk.value_ - vars->sol_ckn.value_ - vars->sol_crr.value_,symweight);
    m_histogram["HypSymLLBarCnkM"]->Fill(-vars->sol_cnk.value_ + vars->sol_ckn.value_ - vars->sol_crr.value_,symweight);

    //m_histogram["HypLLBarCqjP"]->Fill(-vars->sol_cqj.value_ - vars->sol_cjq.value_ - vars->sol_cnn.value_,symweight);
    //m_histogram["HypLLBarCqjM"]->Fill(-vars->sol_cqj.value_ + vars->sol_cjq.value_ - vars->sol_cnn.value_,symweight);
    //m_histogram["HypLLBarCnqP"]->Fill(-vars->sol_cnq.value_ - vars->sol_cqn.value_ - vars->sol_ckj.value_,symweight);
    //m_histogram["HypLLBarCnqM"]->Fill(-vars->sol_cnq.value_ + vars->sol_cqn.value_ - vars->sol_ckj.value_,symweight);
    //m_histogram["HypLLBarCnjP"]->Fill(-vars->sol_cnj.value_ - vars->sol_cjn.value_ - vars->sol_crq.value_,symweight);
    //m_histogram["HypLLBarCnjM"]->Fill(-vars->sol_cnj.value_ + vars->sol_cjn.value_ - vars->sol_crq.value_,symweight);

    ////

    m_histogram["HypSymAntiLeptonBk"]->Fill(-1.0*(vars->sol_b1k.value_),symweight);
    m_histogram["HypSymLeptonBk"]->Fill(-1.0*(vars->sol_b2k.value_),symweight);
    m_histogram["HypSymAntiLeptonBj"]->Fill(-1.0*(vars->sol_b1j.value_),symweight);
    m_histogram["HypSymLeptonBj"]->Fill(-1.0*(vars->sol_b2j.value_),symweight);
    m_histogram["HypSymAntiLeptonBr"]->Fill(-1.0*(vars->sol_b1r.value_),symweight);
    m_histogram["HypSymLeptonBr"]->Fill(-1.0*(vars->sol_b2r.value_),symweight);
    m_histogram["HypSymAntiLeptonBq"]->Fill(-1.0*(vars->sol_b1q.value_),symweight);
    m_histogram["HypSymLeptonBq"]->Fill(-1.0*(vars->sol_b2q.value_),symweight);
    m_histogram["HypSymAntiLeptonBn"]->Fill(-1.0*(vars->sol_b1n.value_),symweight);
    m_histogram["HypSymLeptonBn"]->Fill(-1.0*(vars->sol_b2n.value_),symweight);

    //m_histogram["HypSymLLBarBPnn"]->Fill(-1.0*(vars->sol_bP_nn.value_),symweight);
    //m_histogram["HypSymLLBarBMnn"]->Fill(-1.0*(vars->sol_bM_nn.value_),symweight);
    //m_histogram["HypSymLLBarBPrr"]->Fill(-1.0*(vars->sol_bP_rr.value_),symweight);
    //m_histogram["HypSymLLBarBMrr"]->Fill(-1.0*(vars->sol_bM_rr.value_),symweight);
    //m_histogram["HypSymLLBarBPkk"]->Fill(-1.0*(vars->sol_bP_kk.value_),symweight);
    //m_histogram["HypSymLLBarBMkk"]->Fill(-1.0*(vars->sol_bM_kk.value_),symweight);
    //m_histogram["HypSymLLBarBPjj"]->Fill(-1.0*(vars->sol_bP_jj.value_),symweight);
    //m_histogram["HypSymLLBarBMjj"]->Fill(-1.0*(vars->sol_bM_jj.value_),symweight);
    //m_histogram["HypSymLLBarBPqq"]->Fill(-1.0*(vars->sol_bP_qq.value_),symweight);
    //m_histogram["HypSymLLBarBMqq"]->Fill(-1.0*(vars->sol_bM_qq.value_),symweight);

    m_histogram["HypSymLLBarCkk"]->Fill(-1.0*(vars->sol_ckk.value_),symweight);
    m_histogram["HypSymLLBarCrr"]->Fill(-1.0*(vars->sol_crr.value_),symweight);
    m_histogram["HypSymLLBarCnn"]->Fill(-1.0*(vars->sol_cnn.value_),symweight);
    m_histogram["HypSymLLBarCkj"]->Fill(-1.0*(vars->sol_ckj.value_),symweight);
    m_histogram["HypSymLLBarCrq"]->Fill(-1.0*(vars->sol_crq.value_),symweight);

    m_histogram["HypSymLLBarCrk"]->Fill(-1.0*(vars->sol_crk.value_),symweight);
    m_histogram["HypSymLLBarCkr"]->Fill(-1.0*(vars->sol_ckr.value_),symweight);
    m_histogram["HypSymLLBarCnr"]->Fill(-1.0*(vars->sol_cnr.value_),symweight);
    m_histogram["HypSymLLBarCrn"]->Fill(-1.0*(vars->sol_crn.value_),symweight);
    m_histogram["HypSymLLBarCnk"]->Fill(-1.0*(vars->sol_cnk.value_),symweight);
    m_histogram["HypSymLLBarCkn"]->Fill(-1.0*(vars->sol_ckn.value_),symweight);  

    m_histogram["HypSymLLBarCrj"]->Fill(-1.0*(vars->sol_crj.value_),symweight);
    m_histogram["HypSymLLBarCjr"]->Fill(-1.0*(vars->sol_cjr.value_),symweight);
    //m_histogram["HypSymLLBarCqk"]->Fill(-1.0*(vars->sol_cqk.value_),symweight);
    //m_histogram["HypSymLLBarCkq"]->Fill(-1.0*(vars->sol_ckq.value_),symweight);

    //m_histogram["HypSymLLBarCqj"]->Fill(-1.0*(vars->sol_cqj.value_),symweight);
    //m_histogram["HypSymLLBarCjq"]->Fill(-1.0*(vars->sol_cjq.value_),symweight);
    //m_histogram["HypSymLLBarCnq"]->Fill(-1.0*(vars->sol_cnq.value_),symweight);
    //m_histogram["HypSymLLBarCqn"]->Fill(-1.0*(vars->sol_cqn.value_),symweight);
    //m_histogram["HypSymLLBarCnj"]->Fill(-1.0*(vars->sol_cnj.value_),symweight);
    //m_histogram["HypSymLLBarCjn"]->Fill(-1.0*(vars->sol_cjn.value_),symweight);

    m_histogram["HypSymLLBarCPrk"]->Fill(-1.0*(vars->sol_crk.value_ + vars->sol_ckr.value_),symweight);
    m_histogram["HypSymLLBarCMrk"]->Fill(-1.0*(vars->sol_crk.value_ - vars->sol_ckr.value_),symweight);
    m_histogram["HypSymLLBarCPnr"]->Fill(-1.0*(vars->sol_cnr.value_ + vars->sol_crn.value_),symweight);
    m_histogram["HypSymLLBarCMnr"]->Fill(-1.0*(vars->sol_cnr.value_ - vars->sol_crn.value_),symweight);
    m_histogram["HypSymLLBarCPnk"]->Fill(-1.0*(vars->sol_cnk.value_ + vars->sol_ckn.value_),symweight);
    m_histogram["HypSymLLBarCMnk"]->Fill(-1.0*(vars->sol_cnk.value_ - vars->sol_ckn.value_),symweight);

    m_histogram["HypSymLLBarCPrj"]->Fill(-1.0*(vars->sol_crj.value_ + vars->sol_cjr.value_),symweight);
    m_histogram["HypSymLLBarCMrj"]->Fill(-1.0*(vars->sol_crj.value_ - vars->sol_cjr.value_),symweight);
    //m_histogram["HypSymLLBarCPqj"]->Fill(-1.0*(vars->sol_cqj.value_ + vars->sol_cjq.value_),symweight);
    //m_histogram["HypSymLLBarCMqj"]->Fill(-1.0*(vars->sol_cqj.value_ - vars->sol_cjq.value_),symweight);
    //m_histogram["HypSymLLBarCPnq"]->Fill(-1.0*(vars->sol_cnq.value_ + vars->sol_cqn.value_),symweight);
    //m_histogram["HypSymLLBarCMnq"]->Fill(-1.0*(vars->sol_cnq.value_ - vars->sol_cqn.value_),symweight);
    //m_histogram["HypSymLLBarCPnj"]->Fill(-1.0*(vars->sol_cnj.value_ + vars->sol_cjn.value_),symweight);
    //m_histogram["HypSymLLBarCMnj"]->Fill(-1.0*(vars->sol_cnj.value_ + vars->sol_cjn.value_),symweight);

    m_histogram["HypSymLLBarcHel"]->Fill(-1.0*(vars->sol_cHel.value_),symweight);
    m_histogram["HypSymLLBarcLab"]->Fill(-1.0*(vars->sol_cLab.value_),symweight);

    m_histogram["HypSymLLBarkNorm"]->Fill(-1.0*(vars->sol_kNorm.value_),symweight);
    m_histogram["HypSymLLBarrNorm"]->Fill(-1.0*(vars->sol_rNorm.value_),symweight);
    m_histogram["HypSymLLBarnNorm"]->Fill(-1.0*(vars->sol_nNorm.value_),symweight);
    m_histogram["HypSymLLBarjNorm"]->Fill(-1.0*(vars->sol_jNorm.value_),symweight);
    m_histogram["HypSymLLBarqNorm"]->Fill(-1.0*(vars->sol_qNorm.value_),symweight);

    m_histogram["HypSymLLBarChan"]->Fill(-1.0*(+vars->sol_ckk.value_ - vars->sol_crr.value_ - vars->sol_cnn.value_),symweight);
    m_histogram["HypSymLLBarCsca"]->Fill(-1.0*(-vars->sol_ckk.value_ + vars->sol_crr.value_ - vars->sol_cnn.value_),symweight);
    m_histogram["HypSymLLBarCtra"]->Fill(-1.0*(-vars->sol_ckk.value_ - vars->sol_crr.value_ + vars->sol_cnn.value_),symweight);

    m_histogram["HypSymLLBarCkjL"]->Fill(-1.0*(-vars->sol_ckj.value_ - vars->sol_crr.value_ - vars->sol_cnn.value_),symweight);
    m_histogram["HypSymLLBarCrqL"]->Fill(-1.0*(-vars->sol_ckk.value_ - vars->sol_crq.value_ - vars->sol_cnn.value_),symweight);

    m_histogram["HypSymLLBarCrkP"]->Fill(-1.0*(-vars->sol_crk.value_ - vars->sol_ckr.value_ - vars->sol_cnn.value_),symweight);
    m_histogram["HypSymLLBarCrkM"]->Fill(-1.0*(-vars->sol_crk.value_ + vars->sol_ckr.value_ - vars->sol_cnn.value_),symweight);
    m_histogram["HypSymLLBarCnrP"]->Fill(-1.0*(-vars->sol_cnr.value_ - vars->sol_crn.value_ - vars->sol_ckk.value_),symweight);
    m_histogram["HypSymLLBarCnrM"]->Fill(-1.0*(-vars->sol_cnr.value_ + vars->sol_crn.value_ - vars->sol_ckk.value_),symweight);
    m_histogram["HypSymLLBarCnkP"]->Fill(-1.0*(-vars->sol_cnk.value_ - vars->sol_ckn.value_ - vars->sol_crr.value_),symweight);
    m_histogram["HypSymLLBarCnkM"]->Fill(-1.0*(-vars->sol_cnk.value_ + vars->sol_ckn.value_ - vars->sol_crr.value_),symweight);

    //m_histogram["HypLLBarCqjP"]->Fill(-1.0*(-vars->sol_cqj.value_ - vars->sol_cjq.value_ - vars->sol_cnn.value_),symweight);
    //m_histogram["HypLLBarCqjM"]->Fill(-1.0*(-vars->sol_cqj.value_ + vars->sol_cjq.value_ - vars->sol_cnn.value_),symweight);
    //m_histogram["HypLLBarCnqP"]->Fill(-1.0*(-vars->sol_cnq.value_ - vars->sol_cqn.value_ - vars->sol_ckj.value_),symweight);
    //m_histogram["HypLLBarCnqM"]->Fill(-1.0*(-vars->sol_cnq.value_ + vars->sol_cqn.value_ - vars->sol_ckj.value_),symweight);
    //m_histogram["HypLLBarCnjP"]->Fill(-1.0*(-vars->sol_cnj.value_ - vars->sol_cjn.value_ - vars->sol_crq.value_),symweight);
    //m_histogram["HypLLBarCnjM"]->Fill(-1.0*(-vars->sol_cnj.value_ + vars->sol_cjn.value_ - vars->sol_crq.value_),symweight);

    //Symmetrized End
}

double AnalyzerSpinCorr::getMT2Variable(const LV& lep, const LV& antilep, LV MET){
    Double_t l_1[3], l_2[3], emiss[3];

    l_1[0]=0.0;
    l_1[1]=lep.Px();
    l_1[2]=lep.Py();

    l_2[0]=0.0;
    l_2[1]=antilep.Px();
    l_2[2]=antilep.Py();

    emiss[0]=0.0;
    emiss[1]=MET.Px();
    emiss[2]=MET.Py();

    mt2_bisect::mt2 mt2_event;
    mt2_event.set_momenta(l_1, l_2, emiss);
    mt2_event.set_mn(0.0);
    return mt2_event.get_mt2();

}
