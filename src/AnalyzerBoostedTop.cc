#include <map>
#include <vector>
#include <iostream>

#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <Math/VectorUtil.h>

#include "AnalyzerBoostedTop.h"
#include "analysisStructs.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"
#include "../../common/include/KinematicReconstructionSolution.h"
#include "../../common/include/LooseKinRecoSolution.h"






AnalyzerBoostedTop::AnalyzerBoostedTop(const std::vector<TString>& selectionStepsNoCategories):
AnalyzerBaseClass("BoostedTop_", selectionStepsNoCategories)
{
    std::cout<<"--- Beginning setting up basic histograms\n";
    std::cout<<"=== Finishing setting up basic histograms\n\n";
}



void AnalyzerBoostedTop::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
        //boosted top analysis
        std::vector<double> pTtopBins1={400, 500, 600, 700, 800, 1200};
        std::vector<double> pTtopBins2={400, 600, 800, 1200};
        name = "pTt_true_vs_reco_bin1";
        m_histogram[name] = this->store(new TH2D (prefix_+name+step,"pTt true vs reco ; pTt(reco), [GeV];pTt(true), [GeV]",5,pTtopBins1.data(),5,pTtopBins1.data()));
        name = "pTt_true_vs_reco_bin2";
        m_histogram[name] = this->store(new TH2D (prefix_+name+step,"pTt true vs reco ; pTt(reco), [GeV];pTt(true), [GeV]",3,pTtopBins2.data(),3,pTtopBins2.data()));
        
        //mlblbmet
        name = "gen_mlblbmet";
        m_histogram[name]    = this->store(new TH1D ( prefix_+name+step, "true mlblbmet;mlblbmet, GeV;Nentries", 1000, 0, 1000 ));
        name = "gen_mlblbmet_mtt";
        m_histogram[name]    = this->store(new TH1D ( prefix_+name+step, "true mlblbmet / Mtt;mlblbmet / Mtt;Nentries", 1000, 0, 1 ));
        
        name = "gen_mlblbmet_vs_pTt";
        m_histogram[name] = this->store(new TH2D (prefix_+name+step,"true mlblbmet / Mtt  vs pTt ; pTt(true), GeV; mlblbmet / Mtt (true)",1200,0,1200,1000,0,1));
        name = "gen_mlblbmet_vs_mtt";
        m_histogram[name] = this->store(new TH2D (prefix_+name+step,"true mlblbmet / Mtt  vs mtt ; mtt(true), GeV; mlblbmet / Mtt (true)",2000,0,2000,1000,0,1));
        
}



void AnalyzerBoostedTop::fillHistos(const EventMetadata& ,
                                      const RecoObjects& , const CommonGenObjects& ,
                                      const TopGenObjects& topGenObjects,
                                      const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
				      const LooseKinRecoSolution&,
                                      const ttbar::RecoObjectIndices& , const ttbar::GenObjectIndices&,
                                      const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& ,
				      const double& , const double&, const TString&,
                                      std::map< TString, TH1* >& m_histogram)
{
    
    //boosted top analysis
    
    if(topGenObjects.valuesSet_){
    
        TLorentzVector gentop = common::LVtoTLV((*topGenObjects.GenTop_));
        TLorentzVector gentopbar = common::LVtoTLV((*topGenObjects.GenAntiTop_));
        TLorentzVector genttbar = gentop + gentopbar;
        TLorentzVector gen_b = common::LVtoTLV((*topGenObjects.GenB_));
        TLorentzVector gen_bbar = common::LVtoTLV((*topGenObjects.GenAntiB_));
        //TLorentzVector gen_Wp = common::LVtoTLV((*topGenObjects.GenWPlus_));
        //TLorentzVector gen_Wm = common::LVtoTLV((*topGenObjects.GenWMinus_));
        TLorentzVector gen_al = common::LVtoTLV((*topGenObjects.GenAntiLepton_));
        TLorentzVector gen_l = common::LVtoTLV((*topGenObjects.GenLepton_));
        TLorentzVector gen_neutrino = common::LVtoTLV((*topGenObjects.GenNeutrino_));
        TLorentzVector gen_antineutrino = common::LVtoTLV((*topGenObjects.GenAntiNeutrino_));
        TLorentzVector gen_met = TLorentzVector((gen_neutrino+gen_antineutrino).Px(), (gen_neutrino+gen_antineutrino).Py(), 0., 0.);
        TLorentzVector gen_mlblbmet = gen_b+gen_bbar+gen_al+gen_l+gen_met;
        
        if(gen_b.Pt()>30 && gen_bbar.Pt()>30 && gen_al.Pt()>20 && gen_l.Pt()>20 && fabs(gen_b.Eta())<2.4&& fabs(gen_bbar.Eta())<2.4&& fabs(gen_al.Eta())<2.4&& fabs(gen_l.Eta())<2.4){
            m_histogram["gen_mlblbmet_mtt"]->Fill(gen_mlblbmet.M()/genttbar.M(),genLevelWeights.trueLevelWeight_);
            m_histogram["gen_mlblbmet"]->Fill(gen_mlblbmet.M(),genLevelWeights.trueLevelWeight_);
            
            ((TH2D*)m_histogram["gen_mlblbmet_vs_pTt"])->Fill(gentop.Pt(),gen_mlblbmet.M()/genttbar.M(),genLevelWeights.trueLevelWeight_);
            ((TH2D*)m_histogram["gen_mlblbmet_vs_mtt"])->Fill(genttbar.M(),gen_mlblbmet.M()/genttbar.M(),genLevelWeights.trueLevelWeight_);
            
        }
        
        
        if( kinematicReconstructionSolutions.numberOfSolutions()){
            
            TLorentzVector recotop = common::LVtoTLV(kinematicReconstructionSolutions.solution().top());
            TLorentzVector recotopbar = common::LVtoTLV(kinematicReconstructionSolutions.solution().antiTop());
            TLorentzVector recottbar = recotop + recotopbar;
            
            ((TH2D*)m_histogram["pTt_true_vs_reco_bin1"])->Fill(recotop.Pt(),gentop.Pt());
            ((TH2D*)m_histogram["pTt_true_vs_reco_bin2"])->Fill(recotop.Pt(),gentop.Pt());
            
        }
    }

}
