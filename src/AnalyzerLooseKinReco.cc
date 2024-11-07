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

#include "AnalyzerLooseKinReco.h"
#include "analysisStructs.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"
#include "../../common/include/LooseKinReco.h"
#include "../../common/include/LooseKinRecoSolution.h"




AnalyzerLooseKinReco::AnalyzerLooseKinReco(const Era::Era era,
                                 const std::vector<TString>& selectionStepsNoCategories):
AnalyzerBaseClass("looseKinReco_", selectionStepsNoCategories),
era_(era)              
{
    std::cout<<"--- Beginning setting up basic loose kin reco histograms\n";
    std::cout<<"=== Finishing setting up basic loose kin reco histograms\n\n";
}



void AnalyzerLooseKinReco::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;

	//Loose Kin Reco histograms
	name = "LooseKinRecoTTbarM";
	m_histogram[name] = store(new TH1D(prefix_+name+step, "ttbar mass (loose kin reco)", 500, 0, 1500 ));
	name = "LooseKinRecoTTbarRapidity";
	m_histogram[name] = store(new TH1D(prefix_+name+step, "ttbar rapidity (loose kin reco)", 400, -5, 5 ));
	name = "LooseKinRecoTTbarPt";
	m_histogram[name] = store(new TH1D(prefix_+name+step, "ttbar Pt (loose kin reco)", 800, 0, 800 ));
        
        //Other Loose KinReco output i.e WW mass
        name="LooseWWmass";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "WW mass", 500, 0, 1500 ));
        name="LooseMinMlb";
        m_histogram[name] = store(new TH1D ( prefix_+name+step, "min Mlb mass", 30, 0, 300 ));
        
        //RMS vs TTbarPt, TTbarRapidity, TTbarMass
        name = "LooseRMSvsGenTTbarpT";
        m_histogram[name] = store(new TH2D (prefix_+name+step , "RMS vs Gen", 1000, 0, 1000, 2000, -1000, 1000 ));
        name = "LooseRMSvsGenTTbarRapidity";
        m_histogram[name] = store(new TH2D ( prefix_+name+step, "RMS vs Gen", 400, -5, 5, 400, -5, 5 ));
        name = "LooseRMSvsGenTTBarMass";
        m_histogram[name] = store(new TH2D ( prefix_+name+step, "RMS vs Gen", 3000, 0, 3000, 6000, -3000, 3000 ));
        
}



void AnalyzerLooseKinReco::fillHistos(const EventMetadata&,
                                 const RecoObjects&, const CommonGenObjects&,
				 const TopGenObjects& topGenObjects,
				 const KinematicReconstructionSolutions&, 
				 const LooseKinRecoSolution& looseKinematicReconstructionSolution,
				 const ttbar::RecoObjectIndices&, const ttbar::GenObjectIndices&,
				 const ttbar::GenLevelWeights&, const ttbar::RecoLevelWeights&,
				 const double&, const double& weightLooseKinReco, const TString&,
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

    m_histogram["LooseWWmass"]->Fill(LooseWWmass, weightLooseKinReco);
    m_histogram["LooseMinMlb"]->Fill(LooseMinMlb, weightLooseKinReco);

    m_histogram["LooseKinRecoTTbarM"]->Fill(looseKinRecoTTbarM, weightLooseKinReco);
    m_histogram["LooseKinRecoTTbarRapidity"]->Fill(looseKinRecoTTbarRapidity, weightLooseKinReco);
    m_histogram["LooseKinRecoTTbarPt"]->Fill(looseKinRecoTTbarPt, weightLooseKinReco);

    if(topGenObjects.valuesSet_){

      m_histogram["LooseRMSvsGenTTBarMass"]->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).M(),((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).M()-looseKinRecoTTbar.M());
      m_histogram["LooseRMSvsGenTTbarpT"]->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Pt(),((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Pt()-looseKinRecoTTbar.Pt());
      m_histogram["LooseRMSvsGenTTbarRapidity"]->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Rapidity(),((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Rapidity()-looseKinRecoTTbar.Rapidity());

    }

  } 
    
    
    
}








