#include <map>
#include <vector>
#include <iostream>

#include <TH1.h>
#include <TH1D.h>
#include <TString.h>
#include <Math/VectorUtil.h>

#include "AnalyzerControlPlots.h"
#include "analysisStructs.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"
#include "../../common/include/KinematicReconstructionSolution.h"
#include "../../common/include/LooseKinRecoSolution.h"






AnalyzerControlPlots::AnalyzerControlPlots(const std::vector<TString>& selectionStepsNoCategories):
AnalyzerBaseClass("basic_", selectionStepsNoCategories)
{
    std::cout<<"--- Beginning setting up basic histograms\n";
    std::cout<<"=== Finishing setting up basic histograms\n\n";
}



void AnalyzerControlPlots::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;

    name ="vertMulti";
    m_histogram[name] = this->store(new TH1D ( prefix_+name+step, "Primary Vertex Multiplicity", 50, 0, 50 ));
    name ="vertMulti_noPU";
    m_histogram[name] = this->store(new TH1D ( prefix_+name+step, "Primary Vertex Multiplicity (no Pileup)", 50, 0, 50 ));

    // Leptons
    name = "pfiso_mu";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Moun isolation;Muon isolation;Events",2000,0,2));
    name = "pfiso_e";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Electron isolation;Electron isolation;Events",2000,0,2));
    name = "lepton_multiplicity";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Lepton multiplicity;N leptons;Events",21,-0.5,20.5));
    name = "lepton_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Lepton p_{t};p_{t}^{l} [GeV];Leptons",80,0,400));
    name = "lepton_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Lepton #eta;#eta^{l};Leptons",20,-2.4,2.4));
    name = "lepton_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Lepton #phi;#phi^{l};Leptons",50,-3.2,3.2));
    //name = "lepton_dxy";
    //m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Lepton d_{xy};d_{xy}^{l} [cm];Events",60,-0.25,0.25));
    //name = "lepton_dz";
    //m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Lepton d_{z};d_{z} [cm];Events",60,-1,1));
    name = "lepton2lead_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"p_{t} of 2 leading leptons;p_{t}^{l} [GeV];Leptons", 80, 0, 400));
    name = "lepton2lead_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "#eta of 2 leading leptons;#eta^{l};Leptons",20,-2.4,2.4));
    name = "lepton2lead_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "#phi of 2 leading leptons;#phi^{l};Leptons",50,-3.2,3.2));

    // Leading lepton and antilepton
    name = "lepton1st_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st Lepton p_{t};p_{t}^{l_{1}} [GeV];Leptons",50,0,250));
    name = "lepton1st_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st Lepton #eta;#eta^{l_{1}};Leptons",50,-2.6,2.6));
    name = "lepton1st_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st Lepton #phi;#phi^{l_{1}};Leptons",50,-3.2,3.2));
    name = "lepton2nd_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd Lepton p_{t};p_{t}^{l_{2}} [GeV];Leptons",50,0,250));
    name = "lepton2nd_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd Lepton #eta;#eta^{l_{2}};Leptons",50,-2.6,2.6));
    name = "lepton2nd_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd Lepton #phi;#phi^{l_{2}};Leptons",50,-3.2,3.2));
    // Leading lepton and antilepton
    name = "lepton1st_px";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st Lepton p_{x};p_{x}^{l_{1}} [GeV];Leptons",100,-250.,250));
    name = "lepton1st_py";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st Lepton p_{y};p_{y}^{l_{1}} [GeV];Leptons",100,-250.,250));
    name = "lepton1st_pz";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st Lepton p_{z};p_{z}^{l_{1}} [GeV];Leptons",100,-250.,250));
    name = "lepton1st_m";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st Lepton m;m^{l_{1}} [GeV];Leptons",50,0,250));
    name = "lepton1st_E";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st Lepton E;E^{l_{1}} [GeV];Leptons",50,0,250));
    name = "lepton2nd_px";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd Lepton p_{x};p_{x}^{l_{2}} [GeV];Leptons",100,-250.,250));
    name = "lepton2nd_py";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd Lepton p_{y};p_{y}^{l_{2}} [GeV];Leptons",100,-250.,250));
    name = "lepton2nd_pz";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd Lepton p_{z};p_{z}^{l_{2}} [GeV];Leptons",100,-250.,250));
    name = "lepton2nd_m";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd Lepton m;m^{l_{2}} [GeV];Leptons",50,0,250));
    name = "lepton2nd_E";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd Lepton E;E^{l_{2}} [GeV];Leptons",50,0,250));

    // Dilepton
    name = "dilepton_mass";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Dilepton mass;m^{l^{+}l^{-}} [GeV];Events",100,0,400));
    name = "dilepton_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Dilepton p_{t};p_{t}^{l^{+}l^{-}} [GeV];Events",50,0,300));
    name = "dilepton_rapidity";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Dilepton rapidity;y^{l^{+}l^{-}};Events",50,-2.6,2.6));
    name = "dilepton_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Dilepton #phi;#phi^{l^{+}l^{-}};Events",50,-3.2,3.2));
    name = "dilepton_deltaEta";
    //m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Dilepton #Delta#eta;#eta^{l^{+}}-#eta^{l^{-}};Events",50,-5,5));
    name = "dilepton_deltaPhi";
    //m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Dilepton #Delta#phi;#phi^{l^{+}}-#phi^{l^{-}};Events",50,-3.2,3.2));
    name = "dilepton_deltaDxy";
    //m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Dilepton {#Delta}d_{xy};d_{xy}^{l^{+}}-d_{xy}^{l^{-}} [cm];Events",40,-0.3,0.3));
    //name = "dilepton_deltaDz";
    //m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Dilepton {#Delta}d_{z};d_{z}^{l^{+}}-d_{z}^{l^{-}} [cm];Events",40,-1,1));
    name = "dilepton_alpha";
    //m_histogram[name] = this->store(new TH1D(prefix_+name+step,"Dilepton opening angle #alpha^{l^{+}l^{-}};#alpha^{l^{+}l^{-}};Events",40,0,3.2));

    // Jets
    name = "jet_multiplicity";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Jet Multiplicity;N jets;Events",21,-0.5,20.5));
    name = "jet_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Jet p_{t};p_{t}^{jet} [GeV];Jets",80,0,400));
    name = "jet_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Jet #eta;#eta^{jet};Jets",20,-2.4,2.4));
    name = "jet_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Jet #phi;#phi^{jet};Jets",50,-3.2,3.2));
    name = "jet_btagDiscriminator";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "b-tag Discriminator d;d;Jets",100,-0.1,1.1));
    name = "jet_chargeGlobalPtWeighted";
    //m_histogram[name] = this->store(new TH1D(prefix_+name+step, "jetChargeGlobalPtWeighted c_{glob}^{jet}; c_{glob}^{jet};# jets", 110, -1.1, 1.1));
    name = "jet_chargeRelativePtWeighted";
    //m_histogram[name] = this->store(new TH1D(prefix_+name+step, "jetChargeRelativePtWeighted c_{rel}^{jet}; c_{rel}^{jet};# jets", 110, -1.1, 1.1));
    name = "jet2lead_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"p_{t} of 2 leading jets;p_{t}^{jet} [GeV];Jets", 80, 0, 400));
    name = "jet2lead_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "#eta of 2 leading jets;#eta^{jet};Jets",20,-2.4,2.4));


    // individual jet information
    name = "jet1_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st jet p_{t};p_{t}^{jet_{1}} [GeV];jets",250,0,250));
    name = "jet1_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st jet #eta;#eta^{jet_{1}};jets",50,-2.6,2.6));
    name = "jet1_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st jet #phi;#phi^{jet_{1}};jets",50,-3.2,3.2));
    name = "jet1_px";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st jet p_{x};p_{x}^{jet_{1}} [GeV];jets",500,-250,250));
    name = "jet1_py";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st jet p_{y};p_{y}^{jet_{1}} [GeV];jets",500,-250,250));
    name = "jet1_pz";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st jet p_{z};p_{z}^{jet_{1}} [GeV];jets",500,-250,250));
    name = "jet1_m";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st jet m;m^{jet_{1}} [GeV];jets",250,0,250));
    name = "jet1_E";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "1-st jet E;E^{jet_{1}} [GeV];jets",250,0,250));

    name = "jet2_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd jet p_{t};p_{t}^{jet_{2}} [GeV];jets",250,0,250));
    name = "jet2_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd jet #eta;#eta^{jet_{2}};jets",50,-2.6,2.6));
    name = "jet2_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd jet #phi;#phi^{jet_{2}};jets",50,-3.2,3.2));
    name = "jet2_px";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd jet p_{x};p_{x}^{jet_{2}} [GeV];jets",500,-250,250));
    name = "jet2_py";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd jet p_{y};p_{y}^{jet_{2}} [GeV];jets",500,-250,250));
    name = "jet2_pz";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd jet p_{z};p_{z}^{jet_{2}} [GeV];jets",500,-250,250));
    name = "jet2_m";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd jet m;m^{jet_{2}} [GeV];jets",250,0,250));
    name = "jet2_E";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "2-nd jet E;E^{jet_{2}} [GeV];jets",250,0,250));

    name = "jet3_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "3-rd jet p_{t};p_{t}^{jet_{3}} [GeV];jets",250,0,250));
    name = "jet3_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "3-rd jet #eta;#eta^{jet_{3}};jets",50,-2.6,2.6));
    name = "jet3_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "3-rd jet #phi;#phi^{jet_{3}};jets",50,-3.2,3.2));
    name = "jet3_px";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "3-rd jet p_{x};p_{x}^{jet_{3}} [GeV];jets",500,-250,250));
    name = "jet3_py";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "3-rd jet p_{y};p_{y}^{jet_{3}} [GeV];jets",500,-250,250));
    name = "jet3_pz";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "3-rd jet p_{z};p_{z}^{jet_{3}} [GeV];jets",500,-250,250));
    name = "jet3_m";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "3-rd jet m;m^{jet_{3}} [GeV];jets",250,0,250));
    name = "jet3_E";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "3-rd jet E;E^{jet_{3}} [GeV];jets",250,0,250));
    /*
    name = "jet4_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "4-th jet p_{t};p_{t}^{jet_{4}} [GeV];jets",250,0,250));
    name = "jet4_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "4-th jet #eta;#eta^{jet_{4}};jets",50,-2.6,2.6));
    name = "jet4_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "4-th jet #phi;#phi^{jet_{4}};jets",50,-3.2,3.2));
    name = "jet4_px";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "4-th jet p_{x};p_{x}^{jet_{4}} [GeV];jets",500,-250,250));
    name = "jet4_py";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "4-th jet p_{y};p_{y}^{jet_{4}} [GeV];jets",500,-250,250));
    name = "jet4_pz";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "4-th jet p_{z};p_{z}^{jet_{4}} [GeV];jets",500,-250,250));
    name = "jet4_m";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "4-th jet m;m^{jet_{4}} [GeV];jets",250,0,250));
    name = "jet4_E";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "4-th jet E;E^{jet_{4}} [GeV];jets",250,0,250));
    
    name = "jet5_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "5-th jet p_{t};p_{t}^{jet_{5}} [GeV];jets",250,0,250));
    name = "jet5_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "5-th jet #eta;#eta^{jet_{5}};jets",50,-2.6,2.6));
    name = "jet5_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "5-th jet #phi;#phi^{jet_{5}};jets",50,-3.2,3.2));
    name = "jet5_px";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "5-th jet p_{x};p_{x}^{jet_{5}} [GeV];jets",500,-250,250));
    name = "jet5_py";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "5-th jet p_{y};p_{y}^{jet_{5}} [GeV];jets",500,-250,250));
    name = "jet5_pz";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "5-th jet p_{z};p_{z}^{jet_{5}} [GeV];jets",500,-250,250));
    name = "jet5_m";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "5-th jet m;m^{jet_{5}} [GeV];jets",250,0,250));
    name = "jet5_E";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "5-th jet E;E^{jet_{5}} [GeV];jets",250,0,250));
    
    name = "jet6_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "6-th jet p_{t};p_{t}^{jet_{6}} [GeV];jets",250,0,250));
    name = "jet6_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "6-th jet #eta;#eta^{jet_{6}};jets",50,-2.6,2.6));
    name = "jet6_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "6-th jet #phi;#phi^{jet_{6}};jets",50,-3.2,3.2));
    name = "jet6_px";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "6-th jet p_{x};p_{x}^{jet_{6}} [GeV];jets",500,-250,250));
    name = "jet6_py";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "6-th jet p_{y};p_{y}^{jet_{6}} [GeV];jets",500,-250,250));
    name = "jet6_pz";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "6-th jet p_{z};p_{z}^{jet_{6}} [GeV];jets",500,-250,250));
    name = "jet6_m";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "6-th jet m;m^{jet_{6}} [GeV];jets",50,0,250));
    name = "jet6_E";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "6-th jet E;E^{jet_{6}} [GeV];jets",250,0,250));
    */


    // Bjets
    name = "bjet_multiplicity";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-Jet Multiplicity;N b-jets;Events",21,-0.5,20.5));
    name = "bjet_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-Jet p_{t};p_{t}^{b-jet} [GeV];B-Jets",50,0,300));
    name = "bjet_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-Jet #eta;#eta^{b-jet};B-Jets",50,-2.6,2.6));
    name = "bjet_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-Jet #phi;#phi^{b-jet};B-Jets",50,-3.2,3.2));
    // name = "bjet_chargeRelativePtWeighted";
    // m_histogram[name] = this->store(new TH1D(prefix_+name+step, "B-JetChargeRelativePtWeighted c_{rel}^{jet}; c_{rel}^{jet};# B-Jets", 110, -1.1, 1.1));

    // Met
    name = "met_et";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Met E_{t};E_{t}^{met};Events",500,0,500));
    name = "met_phi";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Met #phi;#phi^{met};Events",100,-3.2,3.2));
    name = "met_eta";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Met #eta;#eta^{met};Events",100,-3.2,3.2));
    name = "met_px";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Met p_{x};p_{x}^{met};Events",1000,-500,500));
    name = "met_py";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Met p_{y};p_{y}^{met};Events",1000,-500,500));
    name = "met_pt";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Met p_{T};p_{T}^{met};Events",500,0.,500.));
    name = "met_m";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Met m;#m^{met};Events",500,0.,500.));
    name = "met_E";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Met E;#E^{met};Events",500,0.,500.));

    name = "met_res";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step, "Met Resolution;MET^{true}-MET^{reco};Events",200,-100,100));

}



void AnalyzerControlPlots::fillHistos(const EventMetadata&,
                                      const RecoObjects& recoObjects, const CommonGenObjects& ,
                                      const TopGenObjects& topGenObjects,
                                      const KinematicReconstructionSolutions&,
				      const LooseKinRecoSolution&,
                                      const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices&,
                                      const ttbar::GenLevelWeights&, const ttbar::RecoLevelWeights& recoLevelWeights,
                                      const double& weight, const double&, const TString&,
                                      std::map< TString, TH1* >& m_histogram)
{

    double weightNoPU = recoLevelWeights.weightNoPileup_;
    m_histogram["vertMulti"]->Fill(recoObjects.vertMulti_, weight);
    m_histogram["vertMulti_noPU"]->Fill(recoObjects.vertMulti_, weightNoPU);

    // Leptons

    for(const int index : recoObjectIndices.allLeptonIndices_){
        if(std::abs(recoObjects.lepPdgId_->at(index))==11)
            m_histogram["pfiso_e"]->Fill(recoObjects.lepIso_->at(index), weight);
        if(std::abs(recoObjects.lepPdgId_->at(index))==13)
            m_histogram["pfiso_mu"]->Fill(recoObjects.lepIso_->at(index), weight);
    }

    m_histogram["lepton_multiplicity"]->Fill(recoObjectIndices.allLeptonIndices_.size(), weight);
    for(const int index : recoObjectIndices.leptonIndices_){
        m_histogram["lepton_pt"]->Fill(recoObjects.allLeptons_->at(index).Pt(), weight);
        m_histogram["lepton_eta"]->Fill(recoObjects.allLeptons_->at(index).Eta(), weight);
        m_histogram["lepton_phi"]->Fill(recoObjects.allLeptons_->at(index).Phi(), weight);
        //m_histogram["lepton_dxy"]->Fill(recoObjects.lepDxyVertex0_->at(index), weight);
    }
    for(const int index : recoObjectIndices.antiLeptonIndices_){
        m_histogram["lepton_pt"]->Fill(recoObjects.allLeptons_->at(index).Pt(), weight);
        m_histogram["lepton_eta"]->Fill(recoObjects.allLeptons_->at(index).Eta(), weight);
        m_histogram["lepton_phi"]->Fill(recoObjects.allLeptons_->at(index).Phi(), weight);
        //m_histogram["lepton_dxy"]->Fill(recoObjects.lepDxyVertex0_->at(index), weight);
    }

    if(recoObjectIndices.allLeptonIndices_.size()>1){
        m_histogram["lepton2lead_eta"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.allLeptonIndices_.at(0)).Eta(), weight);
        m_histogram["lepton2lead_pt"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.allLeptonIndices_.at(0)).Pt(), weight);
        m_histogram["lepton2lead_phi"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.allLeptonIndices_.at(0)).Phi(), weight);
        m_histogram["lepton2lead_eta"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.allLeptonIndices_.at(1)).Eta(), weight);
        m_histogram["lepton2lead_pt"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.allLeptonIndices_.at(1)).Pt(), weight);
        m_histogram["lepton2lead_phi"]->Fill((*recoObjects.allLeptons_).at(recoObjectIndices.allLeptonIndices_.at(1)).Phi(), weight);
    }

    const int leptonIndex = recoObjectIndices.leptonIndices_.size()>0 ? recoObjectIndices.leptonIndices_.at(0) : -1;
    const int antiLeptonIndex = recoObjectIndices.antiLeptonIndices_.size()>0 ? recoObjectIndices.antiLeptonIndices_.at(0) : -1;
    const bool hasLeptonPair = (leptonIndex!=-1 && antiLeptonIndex!=-1);

    //if(hasLeptonPair){
    //    m_histogram["lepton_pt"]->Fill(recoObjects.allLeptons_->at(leptonIndex).Pt(), weight);
    //    m_histogram["lepton_pt"]->Fill(recoObjects.allLeptons_->at(antiLeptonIndex).Pt(), weight);
    //    m_histogram["lepton_eta"]->Fill(recoObjects.allLeptons_->at(leptonIndex).Eta(), weight);
    //    m_histogram["lepton_eta"]->Fill(recoObjects.allLeptons_->at(antiLeptonIndex).Eta(), weight);
    //}

    // Leading lepton and antilepton
    int leadingLeptonIndex(leptonIndex);
    int nLeadingLeptonIndex(antiLeptonIndex);
    if(hasLeptonPair){
        common::orderIndices(leadingLeptonIndex, nLeadingLeptonIndex, *recoObjects.allLeptons_, common::LVpt);

        m_histogram["lepton1st_pt"]->Fill(recoObjects.allLeptons_->at(leadingLeptonIndex).Pt(), weight);
        m_histogram["lepton1st_eta"]->Fill(recoObjects.allLeptons_->at(leadingLeptonIndex).Eta(), weight);
        m_histogram["lepton1st_phi"]->Fill(recoObjects.allLeptons_->at(leadingLeptonIndex).Phi(), weight);
        m_histogram["lepton1st_px"]->Fill(recoObjects.allLeptons_->at(leadingLeptonIndex).Px(), weight);
        m_histogram["lepton1st_py"]->Fill(recoObjects.allLeptons_->at(leadingLeptonIndex).Py(), weight);
        m_histogram["lepton1st_pz"]->Fill(recoObjects.allLeptons_->at(leadingLeptonIndex).Pz(), weight);
        m_histogram["lepton1st_m"]->Fill(recoObjects.allLeptons_->at(leadingLeptonIndex).M(), weight);
        m_histogram["lepton1st_E"]->Fill(recoObjects.allLeptons_->at(leadingLeptonIndex).E(), weight);

        m_histogram["lepton2nd_pt"]->Fill(recoObjects.allLeptons_->at(nLeadingLeptonIndex).Pt(), weight);
        m_histogram["lepton2nd_eta"]->Fill(recoObjects.allLeptons_->at(nLeadingLeptonIndex).Eta(), weight);
        m_histogram["lepton2nd_phi"]->Fill(recoObjects.allLeptons_->at(nLeadingLeptonIndex).Phi(), weight);
        m_histogram["lepton2nd_px"]->Fill(recoObjects.allLeptons_->at(nLeadingLeptonIndex).Px(), weight);
        m_histogram["lepton2nd_py"]->Fill(recoObjects.allLeptons_->at(nLeadingLeptonIndex).Py(), weight);
        m_histogram["lepton2nd_pz"]->Fill(recoObjects.allLeptons_->at(nLeadingLeptonIndex).Pz(), weight);
        m_histogram["lepton2nd_m"]->Fill(recoObjects.allLeptons_->at(nLeadingLeptonIndex).M(), weight);
        m_histogram["lepton2nd_E"]->Fill(recoObjects.allLeptons_->at(nLeadingLeptonIndex).E(), weight);
    }

    // Dilepton
    if(hasLeptonPair){
        LV dilepton(0.,0.,0.,0.);
        dilepton = recoObjects.allLeptons_->at(leadingLeptonIndex) + recoObjects.allLeptons_->at(nLeadingLeptonIndex);

        m_histogram["dilepton_mass"]->Fill(dilepton.M(), weight);
        m_histogram["dilepton_pt"]->Fill(dilepton.Pt(), weight);
        m_histogram["dilepton_rapidity"]->Fill(dilepton.Rapidity(), weight);
        m_histogram["dilepton_phi"]->Fill(dilepton.Phi(), weight);

        //m_histogram["dilepton_deltaEta"]->Fill(recoObjects.allLeptons_->at(leptonIndex).Eta() - recoObjects.allLeptons_->at(antiLeptonIndex).Eta(), weight);
        //m_histogram["dilepton_deltaPhi"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(recoObjects.allLeptons_->at(antiLeptonIndex), recoObjects.allLeptons_->at(leptonIndex)), weight);
        //m_histogram["dilepton_deltaDxy"]->Fill(recoObjects.lepDxyVertex0_->at(leptonIndex) - recoObjects.lepDxyVertex0_->at(antiLeptonIndex), weight);

        //m_histogram["dilepton_alpha"]->Fill(ROOT::Math::VectorUtil::Angle(recoObjects.allLeptons_->at(antiLeptonIndex), recoObjects.allLeptons_->at(leptonIndex)), weight);
    }


    // Jets
    m_histogram["jet_multiplicity"]->Fill(recoObjectIndices.jetIndices_.size(), weight);
    for(const int index : recoObjectIndices.jetIndices_){
        m_histogram["jet_pt"]->Fill(recoObjects.jets_->at(index).Pt(), weight);
        m_histogram["jet_eta"]->Fill(recoObjects.jets_->at(index).Eta(), weight);
        m_histogram["jet_phi"]->Fill(recoObjects.jets_->at(index).Phi(), weight);
        double btagDiscriminator = recoObjects.jetBtags_->at(index);
        if(btagDiscriminator < -0.1) btagDiscriminator = -0.05;
        m_histogram["jet_btagDiscriminator"]->Fill(btagDiscriminator, weight);
        //m_histogram["jet_chargeGlobalPtWeighted"]->Fill(recoObjects.jetChargeGlobalPtWeighted_->at(index), weight);
        //m_histogram["jet_chargeRelativePtWeighted"]->Fill(recoObjects.jetChargeRelativePtWeighted_->at(index), weight);
    }

    if(recoObjectIndices.jetIndices_.size()>1){
        m_histogram["jet2lead_eta"]->Fill((*recoObjects.jets_).at(recoObjectIndices.jetIndices_.at(0)).Eta(), weight);
        m_histogram["jet2lead_pt"]->Fill((*recoObjects.jets_).at(recoObjectIndices.jetIndices_.at(0)).Pt(), weight);
        m_histogram["jet2lead_eta"]->Fill((*recoObjects.jets_).at(recoObjectIndices.jetIndices_.at(1)).Eta(), weight);
        m_histogram["jet2lead_pt"]->Fill((*recoObjects.jets_).at(recoObjectIndices.jetIndices_.at(1)).Pt(), weight);
    }
    for(uint i_jet=0; i_jet< recoObjectIndices.jetIndices_.size(); ++i_jet){
        if(i_jet>2) { continue; }
        const std::string tag = "jet"+std::to_string(1+i_jet);
        const auto& idx = recoObjectIndices.jetIndices_.at(i_jet);

        m_histogram[tag+"_eta"]->Fill((*recoObjects.jets_).at(idx).Eta(), weight);
        m_histogram[tag+"_pt"]->Fill((*recoObjects.jets_).at(idx).Pt(), weight);
        m_histogram[tag+"_phi"]->Fill((*recoObjects.jets_).at(idx).Phi(), weight);
        m_histogram[tag+"_px"]->Fill((*recoObjects.jets_).at(idx).Px(), weight);
        m_histogram[tag+"_py"]->Fill((*recoObjects.jets_).at(idx).Py(), weight);
        m_histogram[tag+"_pz"]->Fill((*recoObjects.jets_).at(idx).Pz(), weight);
        m_histogram[tag+"_m"]->Fill((*recoObjects.jets_).at(idx).M(), weight);
        m_histogram[tag+"_E"]->Fill((*recoObjects.jets_).at(idx).E(), weight);
    }


    // Bjets
    m_histogram["bjet_multiplicity"]->Fill(recoObjectIndices.bjetIndices_.size(), weight);
//     for(const int index : recoObjectIndices.bjetIndices_){
//         m_histogram["bjet_pt"]->Fill(recoObjects.jets_->at(index).Pt(), weight);
//         m_histogram["bjet_eta"]->Fill(recoObjects.jets_->at(index).Eta(), weight);
//         m_histogram["bjet_phi"]->Fill(recoObjects.jets_->at(index).Phi(), weight);
//         m_histogram["bjet_chargeRelativePtWeighted"]->Fill(recoObjects.jetChargeRelativePtWeighted_->at(index), weight);
//     }



    // Met
    const LV& met = *recoObjects.met_;

    m_histogram["met_et"]->Fill(met.Pt(), weight);
    m_histogram["met_phi"]->Fill(met.Phi(), weight);
    m_histogram["met_eta"]->Fill(met.Eta(), weight);
    m_histogram["met_px"]->Fill(met.Px(), weight);
    m_histogram["met_py"]->Fill(met.Py(), weight);
    m_histogram["met_pt"]->Fill(met.pt(), weight);
    m_histogram["met_m"]->Fill(met.M(), weight);
    m_histogram["met_E"]->Fill(met.E(), weight);
    if(topGenObjects.valuesSet_){
        m_histogram["met_res"]->Fill(topGenObjects.GenMet_->Pt() - met.Pt(),weight);
    }

}
