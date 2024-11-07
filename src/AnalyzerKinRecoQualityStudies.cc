#include <map>
#include <vector>
#include <iostream>

#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <Math/VectorUtil.h>

#include "AnalyzerKinRecoQualityStudies.h"
#include "analysisStructs.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"
#include "../../common/include/KinematicReconstructionSolution.h"
#include "../../common/include/KinematicReconstruction_LSroutines.h"
#include "../../common/include/LooseKinReco.h"
#include "../../common/include/LooseKinRecoSolution.h"




AnalyzerKinRecoQualityStudies::AnalyzerKinRecoQualityStudies(const Era::Era era,
                                 const std::vector<TString>& selectionStepsNoCategories):
AnalyzerBaseClass("kinRecoQualityStudies_", selectionStepsNoCategories),
era_(era)              
{
    std::cout<<"--- Beginning setting up basic loose kin reco histograms\n";
    std::cout<<"=== Finishing setting up basic loose kin reco histograms\n\n";
}



void AnalyzerKinRecoQualityStudies::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    
    TString name;

	//Loose Kin Reco histograms
	name = "LooseKinRecoTTbarM";
	m_histogram[name] = store(new TH1D(prefix_+name+step, "ttbar mass (loose kin reco)", 3000, 0, 3000 ));
	name = "LooseKinRecoTTbarRapidity";
	m_histogram[name] = store(new TH1D(prefix_+name+step, "ttbar rapidity (loose kin reco)", 400, -5, 5 ));
	name = "LooseKinRecoTTbarPt";
	m_histogram[name] = store(new TH1D(prefix_+name+step, "ttbar Pt (loose kin reco)", 1200, 0, 1200 ));
        
        //Other Loose KinReco output i.e WW mass
        name="LooseWWmass";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "WW mass", 3000, 0, 3000 ));
        name="LooseMinMlb";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "min Mlb mass", 30, 0, 300 ));
        
        //Loose kin reco: RMS ttbar plots
        name = "LooseGenVsDeltaTrueReco_TTbarpT";
        m_histogram[name] = store(new TH2D (prefix_+name+step , "RMS vs Gen", 1000, 0, 1000, 2000, -1000, 1000 ));
        name = "LooseGenVsDeltaTrueReco_TTbarRapidity";
        m_histogram[name] = store(new TH2D ( prefix_+name+step, "RMS vs Gen", 400, -5, 5, 400, -5, 5 ));
        name = "LooseGenVsDeltaTrueReco_TTBarMass";
        m_histogram[name] = store(new TH2D ( prefix_+name+step, "RMS vs Gen", 3000, 0, 3000, 6000, -3000, 3000 ));

	//Loose kin reco: migration matricies
	name = "LooseRecoVsGenTTBarMass";
	m_histogram[name] = store(new TH2D ( prefix_+name+step , "Reco vs Gen", 2000, 0, 2000, 2000, 0, 2000 ));
	name = "LooseRecoVsGenTTBarpT";
	m_histogram[name] = store(new TH2D ( prefix_+name+step , "Reco vs Gen", 2000, 0, 2000, 2000, 0, 2000 ));
	name = "LooseRecoVsGenTTBarRapidity";
        m_histogram[name] = store(new TH2D ( prefix_+name+step , "Reco vs Gen", 200, -2.5, 2.5, 200, -2.5, 2.5 ));

	//Loose kin reco: gen histograms
	name = "LooseGenTTBarMass";
	m_histogram[name] = store(new TH1D ( prefix_+name+step , "Gen", 2000, 0, 2000 ));
	name = "LooseGenTTBarpT";
	m_histogram[name] = store(new TH1D ( prefix_+name+step , "Gen", 2000, 0, 2000 ));
	name = "LooseGenTTBarRapidity";
	m_histogram[name] = store(new TH1D ( prefix_+name+step , "Gen", 200, -2.5, 2.5 ));

	//Loose ttbar mass at different cuts on LooseMinMlb
	name = "LooseKinRecoTTbarM_LooseMinMlb40";
	m_histogram[name] = store(new TH1D(prefix_+name+step, "ttbar mass (loose kin reco), cut at LooseMinMlb<40", 3000, 0, 3000 ));
	name = "LooseKinRecoTTbarM_LooseMinMlb80";
	m_histogram[name] = store(new TH1D(prefix_+name+step, "ttbar mass (loose kin reco), cut at LooseMinMlb<80", 3000, 0, 3000 ));
	name = "LooseKinRecoTTbarM_LooseMinMlb100";
	m_histogram[name] = store(new TH1D(prefix_+name+step, "ttbar mass (loose kin reco), cut at LooseMinMlb<100", 3000, 0, 3000 ));
	name = "LooseKinRecoTTbarM_LooseMinMlb120";
	m_histogram[name] = store(new TH1D(prefix_+name+step, "ttbar mass (loose kin reco), cut at LooseMinMlb<120", 3000, 0, 3000 ));
	name = "LooseKinRecoTTbarM_LooseMinMlb180";
	m_histogram[name] = store(new TH1D(prefix_+name+step, "ttbar mass (loose kin reco), cut at LooseMinMlb<180", 3000, 0, 3000 ));

	//Full kin reco: RMS top plots
	name = "GenVsDeltaTrueReco_ToppT";
	m_histogram[name] = store(new TH2D ( prefix_+name+step, "RMS vs Gen", 1000, 0, 1000, 2000, -1000, 1000 ));
	name = "GenVsDeltaTrueReco_TopRapidity";
	m_histogram[name] = store(new TH2D ( prefix_+name+step, "RMS vs Gen", 400, -5, 5, 400, -5, 5 ));

	//Full kin reco: RMS ttbar plots
	name = "GenVsDeltaTrueReco_TTBarMass";
	m_histogram[name] = store(new TH2D ( prefix_+name+step, "RMS vs Gen", 3000, 0, 3000, 6000, -3000, 3000 ));
	name = "GenVsDeltaTrueReco_TTBarpT";
	m_histogram[name] = store(new TH2D ( prefix_+name+step, "RMS vs Gen", 2000, 0, 2000, 4000, -2000, 2000 ));
	name = "GenVsDeltaTrueReco_TTBarRapidity";
	m_histogram[name] = store(new TH2D ( prefix_+name+step, "RMS vs Gen", 400, -5, 5, 400, -5, 5 ));
	name = "GenVsDeltaTrueReco_TTBarDeltaRapidity";
	m_histogram[name] = store(new TH2D ( prefix_+name+step, "RMS vs Gen", 400, -5, 5, 400, -5, 5 ));

	//Full kin reco: migration matricies
	name = "RecoVsGenTTBarMass";
	m_histogram[name] = store(new TH2D ( prefix_+name+step, "Reco vs Gen", 2000, 0, 2000, 2000, 0, 2000 ));
	name = "RecoVsGenTTBarpT";
	m_histogram[name] = store(new TH2D ( prefix_+name+step, "Reco vs Gen", 2000, 0, 2000, 2000, 0, 2000 ));
	name = "RecoVsGenTTBarRapidity";
	m_histogram[name] = store(new TH2D ( prefix_+name+step, "Reco vs Gen", 200, -2.5, 2.5, 200, -2.5, 2.5 ));
	name = "RecoVsGenToppT";
	m_histogram[name] = store(new TH2D ( prefix_+name+step, "Reco vs Gen", 1000, 0, 1000, 1000, 0, 1000 ));
	name = "RecoVsGenTopRapidity";
	m_histogram[name] = store(new TH2D ( prefix_+name+step, "Reco vs Gen", 200, -2.5, 2.5, 200, -2.5, 2.5 ));

	//Full kin reco: gen histograms
	name = "GenTTBarMass";
	m_histogram[name] = store(new TH1D ( prefix_+name+step, "Gen", 2000, 0, 2000 ));
	name = "GenTTBarpT";
	m_histogram[name] = store(new TH1D ( prefix_+name+step, "Gen", 2000, 0, 2000 ));
	name = "GenTTBarRapidity";
	m_histogram[name] = store(new TH1D ( prefix_+name+step, "Gen", 200, -2.5, 2.5 ));
	name = "GenToppT";
	m_histogram[name] = store(new TH1D ( prefix_+name+step, "Gen", 1000, 0, 1000 ));
	name = "GenTopRapidity";
	m_histogram[name] = store(new TH1D ( prefix_+name+step, "Gen", 200, -2.5, 2.5 ));

	//Full ttbar mass at different cuts on FullMinMlb
	name = "FullKinRecoTTbarM_FullMinMlb40";
	m_histogram[name] = store(new TH1D(prefix_+name+step, "ttbar mass (loose kin reco), cut at FullMinMlb<40", 3000, 0, 3000 ));
	name = "FullKinRecoTTbarM_FullMinMlb80";
	m_histogram[name] = store(new TH1D(prefix_+name+step, "ttbar mass (loose kin reco), cut at FullMinMlb<80", 3000, 0, 3000 ));
	name = "FullKinRecoTTbarM_FullMinMlb100";
	m_histogram[name] = store(new TH1D(prefix_+name+step, "ttbar mass (loose kin reco), cut at FullMinMlb<100", 3000, 0, 3000 ));
	name = "FullKinRecoTTbarM_FullMinMlb120";
	m_histogram[name] = store(new TH1D(prefix_+name+step, "ttbar mass (loose kin reco), cut at FullMinMlb<120", 3000, 0, 3000 ));
	name = "FullKinRecoTTbarM_FullMinMlb180";
	m_histogram[name] = store(new TH1D(prefix_+name+step, "ttbar mass (loose kin reco), cut at FullMinMlb<180", 3000, 0, 3000 ));

        
}


void AnalyzerKinRecoQualityStudies::fillHistos(const EventMetadata&,
                                 const RecoObjects&, const CommonGenObjects&,
				 const TopGenObjects& topGenObjects,
				 const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
				 const LooseKinRecoSolution& looseKinematicReconstructionSolution,
				 const ttbar::RecoObjectIndices&, const ttbar::GenObjectIndices&,
				 const ttbar::GenLevelWeights&, const ttbar::RecoLevelWeights&,
				 const double& weight, const double& weightLooseKinReco, const TString&,
				 std::map< TString, TH1* >& m_histogram)
{

  if(looseKinematicReconstructionSolution.hasSolution()){
    //Loose kin reco solutions                                                                                                                                                                              
    const LV& looseKinRecoTTbar = looseKinematicReconstructionSolution.TTbar();
    const float looseKinRecoTTbarM = looseKinRecoTTbar.M();
    const float looseKinRecoTTbarRapidity = looseKinRecoTTbar.Rapidity();
    const float looseKinRecoTTbarPt = looseKinRecoTTbar.Pt();
    
    const float LooseWWmass = looseKinematicReconstructionSolution.getWWmass();
    const float LooseMinMlb = looseKinematicReconstructionSolution.getMinMlb();

    //Loose Kin Reco histograms
    m_histogram["LooseKinRecoTTbarM"]->Fill(looseKinRecoTTbarM, weightLooseKinReco);
    m_histogram["LooseKinRecoTTbarRapidity"]->Fill(looseKinRecoTTbarRapidity, weightLooseKinReco);
    m_histogram["LooseKinRecoTTbarPt"]->Fill(looseKinRecoTTbarPt, weightLooseKinReco);

    //Other Loose KinReco output i.e WW mass
    m_histogram["LooseWWmass"]->Fill(LooseWWmass, weightLooseKinReco);
    m_histogram["LooseMinMlb"]->Fill(LooseMinMlb, weightLooseKinReco);

    //Loose ttbar mass at different cuts on LooseMinMlb
    if (LooseMinMlb < 40){
      m_histogram["LooseKinRecoTTbarM_LooseMinMlb40"]->Fill(looseKinRecoTTbarM, weightLooseKinReco);
    }
    if (LooseMinMlb < 80){
      m_histogram["LooseKinRecoTTbarM_LooseMinMlb80"]->Fill(looseKinRecoTTbarM, weightLooseKinReco);
    }
    if (LooseMinMlb < 100){
      m_histogram["LooseKinRecoTTbarM_LooseMinMlb100"]->Fill(looseKinRecoTTbarM, weightLooseKinReco);
    }
    if (LooseMinMlb < 120){
      m_histogram["LooseKinRecoTTbarM_LooseMinMlb120"]->Fill(looseKinRecoTTbarM, weightLooseKinReco);
    }
    if (LooseMinMlb < 180){
      m_histogram["LooseKinRecoTTbarM_LooseMinMlb180"]->Fill(looseKinRecoTTbarM, weightLooseKinReco);
    }

    if(topGenObjects.valuesSet_){

      //Loose kin reco: RMS ttbar plots
      m_histogram["LooseGenVsDeltaTrueReco_TTBarMass"]->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).M(),((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).M()-looseKinRecoTTbar.M());
      m_histogram["LooseGenVsDeltaTrueReco_TTbarpT"]->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Pt(),((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Pt()-looseKinRecoTTbar.Pt());
      m_histogram["LooseGenVsDeltaTrueReco_TTbarRapidity"]->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Rapidity(),((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Rapidity()-looseKinRecoTTbar.Rapidity());

      //Loose kin reco: migration matricies
      m_histogram["LooseRecoVsGenTTBarMass"]->Fill(looseKinRecoTTbar.M(), ((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).M());
      m_histogram["LooseRecoVsGenTTBarpT"]->Fill(looseKinRecoTTbar.Pt(), ((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Pt());
      m_histogram["LooseRecoVsGenTTBarRapidity"]->Fill(looseKinRecoTTbar.Rapidity(), ((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Rapidity());

      //Loose kin reco: gen histograms
      m_histogram["LooseGenTTBarMass"]->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).M());
      m_histogram["LooseGenTTBarpT"]->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Pt());
      m_histogram["LooseGenTTBarRapidity"]->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Rapidity());
	
    }

  }


  if(kinematicReconstructionSolutions.numberOfSolutions()){

    const KinematicReconstructionSolution solution = kinematicReconstructionSolutions.solution();
    const LV& hypTop = solution.top();
    const LV& hypAntiTop = solution.antiTop();
    //const LV& hypTtbar = solution.ttbar();
    const LV& hypLepton = solution.lepton();
    const LV& hypAntiLepton = solution.antiLepton();
    const LV& hypBjet = solution.bjet();
    const LV& hypAntiBjet = solution.antiBjet();
    //const LV& hypNeutrino = solution.neutrino();
    //const LV& hypAntiNeutrino = solution.antiNeutrino();


    const auto lvLLbar = std::pair<LV, LV>(hypLepton, hypAntiLepton);
    const auto lvBBbar = std::pair<LV, LV>(hypBjet, hypAntiBjet);
    std::vector<float> mlbKR(4,0.);
    mlbKR[0] = (lvLLbar.first + lvBBbar.first).M();
    mlbKR[1] = (lvLLbar.first + lvBBbar.second).M();
    mlbKR[2] = (lvLLbar.second + lvBBbar.first).M();
    mlbKR[3] = (lvLLbar.second + lvBBbar.second).M();
    float mlbmax1 = std::max(mlbKR[0], mlbKR[3]);
    float mlbmax2 = std::max(mlbKR[1], mlbKR[2]);
    const float FullMinMlb = std::min(mlbmax1, mlbmax2);

    //Full ttbar mass at different cuts on FullMinMlb
    if (FullMinMlb < 40){
      m_histogram["FullKinRecoTTbarM_FullMinMlb40"]->Fill((hypTop+hypAntiTop).M(), weight);
    }
    if (FullMinMlb < 80){
      m_histogram["FullKinRecoTTbarM_FullMinMlb80"]->Fill((hypTop+hypAntiTop).M(), weight);
    }
    if (FullMinMlb < 100){
      m_histogram["FullKinRecoTTbarM_FullMinMlb100"]->Fill((hypTop+hypAntiTop).M(), weight);
    }
    if (FullMinMlb < 120){
      m_histogram["FullKinRecoTTbarM_FullMinMlb120"]->Fill((hypTop+hypAntiTop).M(), weight);
    }
    if (FullMinMlb < 180){
      m_histogram["FullKinRecoTTbarM_FullMinMlb180"]->Fill((hypTop+hypAntiTop).M(), weight);
    }
    
    if(topGenObjects.valuesSet_){

      //Full kin reco: RMS top plots 
      m_histogram["GenVsDeltaTrueReco_ToppT"]->Fill((*topGenObjects.GenTop_).Pt(),(*topGenObjects.GenTop_).Pt()-hypTop.Pt());
      m_histogram["GenVsDeltaTrueReco_ToppT"]->Fill((*topGenObjects.GenAntiTop_).Pt(),(*topGenObjects.GenAntiTop_).Pt()-hypAntiTop.Pt());
      m_histogram["GenVsDeltaTrueReco_TopRapidity"]->Fill((*topGenObjects.GenTop_).Rapidity(),(*topGenObjects.GenTop_).Rapidity()-hypTop.Rapidity());
      m_histogram["GenVsDeltaTrueReco_TopRapidity"]->Fill((*topGenObjects.GenAntiTop_).Rapidity(),(*topGenObjects.GenAntiTop_).Rapidity()-hypAntiTop.Rapidity());
    
      //Full kin reco: RMS ttbar plots
      m_histogram["GenVsDeltaTrueReco_TTBarMass"]->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).M(),((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).M()-(hypTop+hypAntiTop).M());
      m_histogram["GenVsDeltaTrueReco_TTBarpT"]->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Pt(),((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Pt()-(hypTop+hypAntiTop).Pt());
      m_histogram["GenVsDeltaTrueReco_TTBarRapidity"]->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Rapidity(),((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Rapidity()-(hypTop+hypAntiTop).Rapidity());
      m_histogram["GenVsDeltaTrueReco_TTBarDeltaRapidity"]->Fill(std::abs((*topGenObjects.GenTop_).Rapidity()) - std::abs((*topGenObjects.GenAntiTop_).Rapidity()), (std::abs((*topGenObjects.GenTop_).Rapidity()) - std::abs((*topGenObjects.GenAntiTop_).Rapidity())) - (std::abs(hypTop.Rapidity()) - std::abs(hypAntiTop.Rapidity())));
      
      //Full kin reco: migration matricies 
      m_histogram["RecoVsGenTTBarMass"]->Fill((hypTop+hypAntiTop).M(), ((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).M());
      m_histogram["RecoVsGenTTBarpT"]->Fill((hypTop+hypAntiTop).Pt(), ((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Pt());
      m_histogram["RecoVsGenTTBarRapidity"]->Fill((hypTop+hypAntiTop).Rapidity(), ((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Rapidity());
      m_histogram["RecoVsGenToppT"]->Fill(hypTop.Pt(), (*topGenObjects.GenTop_).Pt());
      m_histogram["RecoVsGenToppT"]->Fill(hypAntiTop.Pt(), (*topGenObjects.GenAntiTop_).Pt());
      m_histogram["RecoVsGenTopRapidity"]->Fill(hypTop.Rapidity(), (*topGenObjects.GenTop_).Rapidity());
      m_histogram["RecoVsGenTopRapidity"]->Fill(hypAntiTop.Rapidity(), (*topGenObjects.GenAntiTop_).Rapidity());

      //Full kin reco: gen histograms
      m_histogram["GenTTBarMass"]->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).M());
      m_histogram["GenTTBarpT"]->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Pt());
      m_histogram["GenTTBarRapidity"]->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Rapidity());
      m_histogram["GenToppT"]->Fill((*topGenObjects.GenTop_).Pt());
      m_histogram["GenToppT"]->Fill((*topGenObjects.GenAntiTop_).Pt());
      m_histogram["GenTopRapidity"]->Fill((*topGenObjects.GenTop_).Rapidity());
      m_histogram["GenTopRapidity"]->Fill((*topGenObjects.GenAntiTop_).Rapidity());
      
    }
    
  }
    
    
}








