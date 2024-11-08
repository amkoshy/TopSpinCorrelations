#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <cmath>
#include <limits>
#include <iomanip>

#include <TTree.h>
#include <TString.h>
#include <TMath.h>
#include <TSystem.h>
#include <Math/VectorUtil.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TIterator.h>
#include <TObject.h>

#include "TopAnalysis.h"
#include "AnalysisConfig.h"
#include "HistoListReader.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"
#include "../../common/include/ScaleFactors.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/KinematicReconstructionSolution.h"
#include "../../common/include/LooseKinRecoSolution.h"
#include "analysisStructs.h"
#include "AnalyzerBaseClass.h" //FIXME: rename to AnalyzerBase.
#include "TreeHandlerBase.h"
#include "utils.h"



///top production xsec in pb
constexpr double TOPXSEC = 234.;

///do we want to run the sync excercise?
constexpr bool RUNSYNC = false;

/// rho_s (mass measurement)
constexpr double rho0 = 340.0;

/// Flag of the boosted top analysis
constexpr bool btopTag = false;





TopAnalysis::TopAnalysis(const AnalysisConfig& analysisConfig, TTree*):
AnalysisBase(analysisConfig.general().era_,
             analysisConfig.selections().btagAlgorithm_,
             analysisConfig.selections().btagWorkingPoint_,
             analysisConfig.selections().metAlgorithm_,
             analysisConfig.sampleComposition().addSingleLeptonTriggerInfo_),
analysisConfig_(analysisConfig),
kinRecoOnTheFly_(false),
doClosureTest_(false),
closureFunction_(nullptr),
closureMaxEvents_(0),
runViaTau_(false),
binnedControlPlots_(0)
{
  applyMuonEnergyCorrections_     = analysisConfig_.corrections().applyMuonEnergyCorrections_;
  applyElectronEnergyCorrections_ = analysisConfig_.corrections().applyElectronEnergyCorrections_;
  applyBJetEnergyRegression_ = analysisConfig_.corrections().applyBJetEnergyRegression_;
  applyTopSlope2DGenLevelReweighting_ = analysisConfig_.corrections().applyTopSlope2DGenLevelReweighting_;
  applyTopPtReweighting_ = analysisConfig_.corrections().topPtReweighting_;
  propagateCorrectionsToMET_ = analysisConfig_.corrections().propagateCorrectionsToMET_;
  correctAllTTbarBR_ = false;
}



void TopAnalysis::Begin(TTree*)
{
    // Defaults from AnalysisBase
    AnalysisBase::Begin(0);

    // Set up selection steps of tree handlers
    for(TreeHandlerBase* treeHandler : v_treeHandler_){

      //        if(treeHandler) treeHandler->book();
      //ajeeta -- create TFile before TTrees
      if(treeHandler){
	treeHandler->setOutputFile(this->outputFilename(), this->channel(), this->systematic());
	treeHandler->book();
      }

    }
}



void TopAnalysis::Terminate()
{
    // Produce b-tag efficiencies if required for given correction mode
    this->produceBtagEfficiencies();

    // Do everything needed for TTree
    for(TreeHandlerBase* treeHandler : v_treeHandler_){
        if(treeHandler){
            // Produce and write tree
	  //            treeHandler->writeTrees(this->outputFilename(), this->channel(), this->systematic());
            //treeHandler->writeTrees(fOutput);

	  //ajeeta
	  TH1* h_weightedEvents = getWeightedEvents();
	  treeHandler->writeTrees(h_weightedEvents);

            // Cleanup
            treeHandler->clear();
        }
    }

    // Defaults from AnalysisBase
    AnalysisBase::Terminate();
}



void TopAnalysis::SlaveBegin(TTree*)
{
    // Defaults from AnalysisBase
    AnalysisBase::SlaveBegin(0);

    // Pile-Up Reweighting (Step 3 and Step 8)
    h_vertMulti_step3 = store(new TH1D ( "vertMulti_step3", "Primary Vertex Multiplicity", 50, 0, 50 ));
    h_vertMulti_noPU_step3 = store(new TH1D ( "vertMulti_noPU_step3", "Primary Vertex Multiplicity (no Pileup)", 50, 0, 50 ));
    h_vertMulti_step8 = store(new TH1D ( "vertMulti_step8", "Primary Vertex Multiplicity", 50, 0, 50 ));
    h_vertMulti_noPU_step8 = store(new TH1D ( "vertMulti_noPU_step8", "Primary Vertex Multiplicity (no Pileup)", 50, 0, 50 ));

    // Full Dilepton System / pre-Zpeak cut (Step 3)
    h_diLepMassFull = store(new TH1D ( "DIMFull", "DiLepton Mass (Full Range)", 100, 0, 300 ));

    // Dilepton System / pre-jetMult cut (Step 4)
    h_jetMulti_diLep = store(new TH1D ( "HypjetMulti_diLep", "Jet Multiplicity (after dilepton)", 10, -0.5, 9.5 ));
    h_jetMulti_Zwindow = store(new TH1D ( "HypjetMulti_Zwindow", "Jet Multiplicity (in Z-Window after dilepton)", 10, -0.5, 9.5 ));
    h_MET_diLep = store(new TH1D ( "MET_diLep", "Missing Transverse Energy", 500, 0, 500 ));
    h_MuonpT_diLep = store(new TH1D ( "MuonpT_diLep", "Muon pT", 80, 0, 400 ));
    h_MuonEta_diLep = store(new TH1D ( "MuonEta_diLep", "Muon Eta", 100, -5, 5 ));
    h_ElectronpT_diLep = store(new TH1D ( "ElectronpT_diLep", "Electron pT", 80, 0, 400 ));
    h_ElectronEta_diLep = store(new TH1D ( "ElectronEta_diLep", "Electron Eta", 100, -5, 5 ));
    h_LeptonpT_diLep = store(new TH1D ( "LeptonpT_diLep", "Lepton pT (after dilepton cut)", 80, 0, 400 ));
    h_LeptonEta_diLep = store(new TH1D ( "LeptonEta_diLep", "Lepton Eta (after dilepton cut)", 100, -5, 5 ));
    h_AntiLeptonpT_diLep = store(new TH1D ( "AntiLeptonpT_diLep", "Lepton pT (after dilepton cut)", 80, 0, 400 ));
    h_AntiLeptonEta_diLep = store(new TH1D ( "AntiLeptonEta_diLep", "Lepton Eta (after dilepton cut)", 100, -5, 5 ));

    // pre-MET cut (Step 5)
    h_MET_preMETcut = store(new TH1D ( "MET_preMETcut", "Missing Transverse Energy", 500, 0, 500 ));
    h_MET_Zwindow = store(new TH1D ( "MET_Zwindow", "Missing Transverse Energy (in Z-Window after jet multiplicity cut)", 500, 0, 500 ));

    // pre-btag cut (Step 6)
    h_LeptonpT_postMETcut = store(new TH1D ( "LeptonpT_postMETcut", "Lepton pT (post MET cut)", 80, 0, 400 ));
    h_LeptonEta_postMETcut = store(new TH1D ( "LeptonEta_postMETcut", "Lepton Eta (post MET cut)", 100, -5, 5 ));
    h_AntiLeptonpT_postMETcut = store(new TH1D ( "AntiLeptonpT_postMETcut", "AntiLepton pT (post MET cut)", 80, 0, 400 ));
    h_AntiLeptonEta_postMETcut = store(new TH1D ( "AntiLeptonEta_postMETcut", "AntiLepton Eta (post MET cut)", 100, -5, 5 ));
    h_MuonpT_postMETcut = store(new TH1D ( "MuonpT_postMETcut", "Muon pT (post MET cut)", 80, 0, 400 ));
    h_MuonEta_postMETcut = store(new TH1D ( "MuonEta_postMETcut", "Muon Eta (post MET cut)", 100, -5, 5 ));
    h_ElectronpT_postMETcut = store(new TH1D ( "ElectronpT_postMETcut", "Electron pT (post MET cut)", 80, 0, 400 ));
    h_ElectronEta_postMETcut = store(new TH1D ( "ElectronEta_postMETcut", "Electron Eta (post MET cut)", 100, -5, 5 ));

    h_jetMulti_noBTag = store(new TH1D ( "HypjetMulti_noBTag", "Jet Multiplicity (after dilepton before btag)", 10, -0.5, 9.5 ));
    h_BjetMulti_noBTag = store(new TH1D ( "HypBjetMulti_noBTag", "B-Jet Multiplicity before B-Tag requirement", 10, -0.5, 9.5 ));
    h_BjetMulti_Zwindow = store(new TH1D ( "HypBjetMulti_Zwindow", "B-Jet Multiplicity before B-Tag requirement (in Z-Window after MET cut)", 10, -0.5, 9.5 ));

    // pre-KinReco (Step 7)
    h_LeptonPt_BeforeKinReco = store(new TH1D ( "LeptonpT_BeforeKinReco", "Lepton pT (before kin reco)", 80, 0, 400 ));
    h_LeptonEta_BeforeKinReco = store(new TH1D ( "LeptonEta_BeforeKinReco", "Lepton #eta (before kin reco)", 80, -2.5, 2.5 ));
    h_MuonpT_BeforeKinReco = store(new TH1D ( "MuonpT_BeforeKinReco", "Muon pT (post MET cut)", 80, 0, 400 ));
    h_MuonEta_BeforeKinReco = store(new TH1D ( "MuonEta_BeforeKinReco", "Muon Eta (post MET cut)", 100, -5, 5 ));
    h_ElectronpT_BeforeKinReco = store(new TH1D ( "ElectronpT_BeforeKinReco", "Electron pT (post MET cut)", 80, 0, 400 ));
    h_ElectronEta_BeforeKinReco = store(new TH1D ( "ElectronEta_BeforeKinReco", "Electron Eta (post MET cut)", 100, -5, 5 ));
    h_jetMulti = store(new TH1D ( "HypjetMulti_step7", "Jet Multiplicity", 10, -0.5, 9.5 ));
    h_jetpT = store(new TH1D ( "jetpT", "jet pT", 80, 0, 400 ));
    h_jetHT = store(new TH1D ( "jetHT", "jet HT", 80, 0, 1000 ));
    h_MET_step7 = store(new TH1D ( "MET_step7", "Missing Transverse Energy (before kin reco)", 500, 0, 500 ));
    h_BjetMulti = store(new TH1D ( "HypBjetMulti_step7", "B-Jet Multiplicity", 10, -0.5, 9.5 ));
    h_BjetpT_BeforeKinReco= store(new TH1D ( "BjetpT_BeforeKinReco", "B pT after b-tag", 80, 0, 400 ));
    h_BjetEta_BeforeKinReco = store(new TH1D ( "BjetEta_BeforeKinReco", "b-jet eta (before kin reco)", 80, -2.5, 2.5 ));

    // post-KinReco (Step 8)
    h_diLepMassFull_fullSel = store(new TH1D ( "DIMFull_fullSel", "DiLepton Mass (Full Range)", 100, 0, 300 ));

    h_MET_step8 = store(new TH1D ( "MET_step8", "Missing Transverse Energy (after kin reco)", 500, 0, 500 ));

    h_leptonPt_AfterKinReco = store(new TH1D ( "LeptonpT_AfterKinReco", "Lepton pT (after kin reco)", 80, 0, 400 ));
    h_leptonEta_AfterKinReco = store(new TH1D ( "LeptonEta_AfterKinReco", "Lepton #eta (after kin reco)", 80, -2.5, 2.5 ));
    h_bjeteta_AfterKinReco = store(new TH1D ( "BjetEta_AfterKinReco", "b-jet eta (after kin reco)", 80, -2.5, 2.5 ));



    h_LeptonpT = store(new TH1D ( "LeptonpT", "Lepton pT", 80, 0, 400 ));
    h_LeptonEta = store(new TH1D ( "LeptonEta", "Lepton Eta", 100, -5, 5 ));
    h_AntiLeptonpT = store(new TH1D ( "AntiLeptonpT", "AntiLepton pT", 80, 0, 400 ));
    h_AntiLeptonEta = store(new TH1D ( "AntiLeptonEta", "AntiLepton Eta", 100, -5, 5 ));

    h_HypLLBarDPhi = store(new TH1D ( "HypLLBarDPhi", "#Delta#phi(Lep, AntiLep) (HYP)",3200, 0., 3.2 ));
    h_HypLLBarMass = store(new TH1D ( "HypLLBarMass", "Mass of LLbar System (HYP)", 500, 0, 1000 ));
    h_HypLLBarpT = store(new TH1D ( "HypLLBarpT", "pT of LLbar System (HYP)", 200, 0, 1000 ));

    h_HypLeptonantiBjetMass = store(new TH1D ( "HypLeptonBjetMass", "Mass(Lep, AntiBJet) (HYP)", 500, 0, 1000 ));
    h_HypAntiLeptonBjetMass = store(new TH1D ( "HypAntiLeptonBjetMass", "Mass(AntiLep, BJet) (HYP)", 500, 0, 1000 ));

    h_HypJetMult = store(new TH1D ( "HypJetMult", "Jet Multiplicity (HYP)", 26, -0.5, 25.5 ));
    h_jetMulti_kinReco = store(new TH1D ( "HypjetMulti_step8", "Jet Multiplicity", 10, -0.5, 9.5 ));
    h_bjetMulti_kinReco = store(new TH1D ( "HypBjetMulti", "BJet Multiplicity", 10, -0.5, 9.5 ));
    h_exjetMulti_kinReco = store(new TH1D ( "HypExjetMulti", "ExJet Multiplicity", 10, -0.5, 9.5 ));

    h_HypBJetpT = store(new TH1D ( "HypBJetpT", "B Hypothesis pT", 1200, 0, 1200 ));
    h_HypBJetEta = store(new TH1D ( "HypBJetEta", "B Hypothesis Eta", 200, -5, 5 ));
    h_HypBJetRapidity = store(new TH1D ( "HypBJetRapidity", "B Hypothesis Eta", 100, -5, 5 ));
    h_HypBJetCSVdiscriminator = store(new TH1D ( "HypBJetCSVdiscriminator", "B Hypothesis CSV discriminator", 100, 0., 1. ));

    h_HypAntiBJetpT = store(new TH1D ( "HypAntiBJetpT", "AntiB Hypothesis pT", 1200, 0, 1200 ));
    h_HypAntiBJetEta = store(new TH1D ( "HypAntiBJetEta", "AntiB Hypothesis Eta", 200, -5, 5 ));
    h_HypAntiBJetRapidity = store(new TH1D ( "HypAntiBJetRapidity", "AntiB Hypothesis Eta", 100, -5, 5 ));
    h_HypAntiBJetCSVdiscriminator = store(new TH1D ( "HypAntiBJetCSVdiscriminator", "AntiB Hypothesis CSV discriminator", 100, 0., 1. ));

    h_HypToppT = store(new TH1D ( "HypToppT", "Top pT",1200,0,1200 ));
    h_HypTopEta = store(new TH1D ( "HypTopEta", "Top #eta", 100, -5, 5 ));
    h_HypTopMass = store(new TH1D ( "HypTopMass", "Top Mass", 80, 0, 400 ));
    h_HypTopRapidity = store(new TH1D ( "HypTopRapidity", "Top Rapidity", 400, -5, 5 ));
    h_HypTopRapidityAbs = store(new TH1D ( "HypTopRapidityAbs", "Absolute Top Rapidity", 200, 0, 5 ));

    h_HypAntiToppT = store(new TH1D ( "HypAntiToppT", "AntiTop pT",1200,0,1200 ));
    h_HypAntiTopEta = store(new TH1D ( "HypAntiTopEta", "AntiTop #eta", 100, -5, 5 ));
    h_HypAntiTopMass = store(new TH1D ( "HypAntiTopMass", "AntiTop Mass", 80, 0, 400 ));
    h_HypAntiTopRapidity = store(new TH1D ( "HypAntiTopRapidity", "Top Rapidity", 400, -5, 5 ));
    h_HypAntiTopRapidityAbs = store(new TH1D ( "HypAntiTopRapidityAbs", "Absolute Top Rapidity", 200, 0, 5 ));

    h_HypTTBarRapidity = store(new TH1D ( "HypTTBarRapidity", "Rapidity of TTbar System (HYP)", 400, -5, 5 ));
    h_HypTTBarRapidity_Loose = store(new TH1D ( "HypTTBarRapidity_Loose", "Rapidity of TTbar System (HYP)", 400, -5, 5 ));
    h_HypTTBarRapidityAbs = store(new TH1D ( "HypTTBarRapidityAbs", "Absolute Rapidity of TTbar System (HYP)", 200, 0, 5 ));
    h_HypTTBarRapidityAbs_Loose = store(new TH1D ( "HypTTBarRapidityAbs_Loose", "Absolute Rapidity of TTbar System (HYP)", 200, 0, 5 ));
    h_HypTTBarpT = store(new TH1D ( "HypTTBarpT", "pT of TTbar System (HYP)", 1200, 0, 1200 ));
    h_HypTTBarpT_Loose = store(new TH1D ( "HypTTBarpT_Loose", "pT of TTbar System (HYP)", 1200, 0, 1200 ));
    h_HypTTBarMass = store(new TH1D ( "HypTTBarMass", "Mass of TTbar System (HYP)", 3000, 0, 3000 ));
    h_HypTTBarMass_Loose = store(new TH1D ( "HypTTBarMass_Loose", "Mass of TTbar System (HYP)", 3000, 0, 3000 ));

    //h_jetMultiAll = store(new TH1D ( "HypjetMultiAll", "Jet Multiplicity (AllJets)", 10, -0.5, 9.5 ));
    h_jetMultiXSec = store(new TH1D ( "HypjetMultiXSec", "Jet Multiplicity (for cross-section)", 10, -0.5, 9.5 ));
    h_jetMultiXSec_Loose = store(new TH1D ( "HypjetMultiXSec_Loose", "Jet Multiplicity (for cross-section) looseKR", 10, -0.5, 9.5 ));
    h_jetMultiNoPU = store(new TH1D ( "HypjetMultiNoPU", "Jet Multiplicity (No Pileup or lumi weight)", 10, -0.5, 9.5 ));
    //h_jetMultiVisTop = store(new TH1D ( "HypjetMultiVisTop", "Jet Multiplicity for Visible Top (No Pileup or lumi Weight)", 10, -0.5, 9.5 ));

    h_exjetMulti = store(new TH1D ( "HypexjetMulti", "Extra Jet Multiplicity", 10, -0.5, 9.5 ));
    //    h_exjetMulti_1 = store(new TH1D ( "HypExjetMulti_1", "ExJet Multiplicity", 10, -0.5, 9.5 ));
    //    h_exjetMulti_2 = store(new TH1D ( "HypExjetMulti_2", "ExJet Multiplicity", 10, -0.5, 9.5 ));
    
    h_VisGenTTBarMass = store(new TH1D ( "VisGenTTBarMass", "Mass of TTbar System(VisGEN)", 3000, 0, 3000 ));
    h_VisGenTTBarRapidity = store(new TH1D ( "VisGenTTBarRapidity", "Rapidity of TTbar System(VisGEN)", 400, -5, 5 ));
    h_VisGenTTBarRapidityAbs = store(new TH1D ( "VisGenTTBarRapidityAbs", "Absolute Rapidity of TTbar System(VisGEN)", 200, 0, 5 ));
    h_VisGenTTBarpT = store(new TH1D ( "VisGenTTBarpT", "pT of TTbar System(VisGEN)", 1200, 0, 1200 ));
    h_VisGenTopRapidity = store(new TH1D ( "VisGenTopRapidity", "Rapidity of Top(VisGEN)", 400, -5, 5 ));
    h_VisGenAntiTopRapidity = store(new TH1D ( "VisGenAntiTopRapidity", "Rapidity of AntiTop(VisGEN)", 400, -5, 5 ));
    h_VisGenTopRapidityAbs = store(new TH1D ( "VisGenTopRapidityAbs", "Absolute Rapidity of Top(VisGEN)", 200, 0, 5 ));
    h_VisGenAntiTopRapidityAbs = store(new TH1D ( "VisGenAntiTopRapidityAbs", "Absolute Rapidity of AntiTop(VisGEN)", 200, 0, 5 ));

    h_VisGenLLBarpT = store(new TH1D ( "VisGenLLBarpT", "pT of LLbar System(VisGEN)", 200, 0, 1000 ));
    h_VisGenLLBarMass = store(new TH1D ( "VisGenLLBarMass", "Mass of LLbar System(VisGEN)", 500, 0, 1000 ));

    h_RecoTTBarMass = store(new TH1D ( "RecoTTBarMass","Mass of TTbar System (HYP)",3000,0,3000 ));
    h_RecoTTBarRapidity = store(new TH1D ( "RecoTTBarRapidity","Rapidity of TTbar System (HYP)",400,-5,5 ));
    h_RecoTTBarRapidityAbs = store(new TH1D ( "RecoTTBarRapidityAbs","Absolute Rapidity of TTbar System (HYP)",200,0,5 ));
    h_RecoTTBarpT = store(new TH1D ( "RecoTTBarpT","pT of TTbar System (HYP)",1200,0,1200 ));
    h_RecoToppT = store(new TH1D ( "RecoToppT","pT of Top (HYP)",1200,0,1200 ));
    h_RecoAntiToppT = store(new TH1D ( "RecoAntiToppT","pT of AntiTop (HYP)",1200,0,1200 ));
    h_RecoTopRapidity = store(new TH1D ( "RecoTopRapidity","Rapidity of Top (HYP)",400,-5,5 ));
    h_RecoAntiTopRapidity = store(new TH1D ( "RecoAntiTopRapidity","Rapidity of AntiTop (HYP)",400,-5,5 ));
    h_RecoTopRapidityAbs = store(new TH1D ( "RecoTopRapidityAbs","Absolute Rapidity of Top (HYP)",200,0,5 ));
    h_RecoAntiTopRapidityAbs = store(new TH1D ( "RecoAntiTopRapidityAbs","Absolute Rapidity of AntiTop (HYP)",200,0,5 ));

    h_RecoBJetpT = store(new TH1D ( "RecoBJetpT","pT of BJet (HYP)",1200,0,1200 ));
    h_RecoAntiBJetpT = store(new TH1D ( "RecoAntiBJetpT","pT of AntiBJet (HYP)",1200,0,1200 ));
    h_RecoBJetRapidity = store(new TH1D ( "RecoBJetRapidity","Rapidity of BJet (HYP)",100,-5,5 ));
    h_RecoAntiBJetRapidity = store(new TH1D ( "RecoAntiBJetRapidity","Rapidity of AntiBJet (HYP)",100,-5,5 ));
    h_RecoBJetEta = store(new TH1D ( "RecoBJetEta","#eta of BJet (HYP)",200,-5,5 ));
    h_RecoAntiBJetEta = store(new TH1D ( "RecoAntiBJetEta","#eta of AntiBJet (HYP)",200,-5,5 ));

    h_RecoLLBarMass = store(new TH1D ( "RecoLLBarMass","Mass of LLbar System (HYP)",500,0,1000 ));
    h_RecoLLBarpT = store(new TH1D ( "RecoLLBarpT","pT of LLbar System (HYP)",200,0,1000 ));

    h_VisGenAll = store(new TH1D ( "VisGenAll", "All Visible Generated particles (IM)", 40, 0, 400 ));
    h_VisGenAll_noweight = store(new TH1D ( "VisGenAll_noweight", "All Visible Generated particles (IM)", 40, 0, 400 ));
    h_GenAll = store(new TH1D ( "GenAll", "AllGenerated particles (IM)", 40, 0, 400 ));         h_GenAll->Sumw2();
    h_GenAll_noweight = store(new TH1D ( "GenAll_noweight", "AllGenerated particles (IM)", 40, 0, 400 ));         h_GenAll_noweight->Sumw2();
    h_GenAll_RecoCuts = store(new TH1D ( "GenAll_RecoCuts", "AllGenerated particles (IM)", 40, 0, 400 ));         h_GenAll_RecoCuts->Sumw2();
    h_GenAll_RecoCuts_noweight= store(new TH1D ( "GenAll_RecoCuts_noweight", "AllGenerated particles (IM)", 40, 0, 400 ));         h_GenAll_RecoCuts_noweight->Sumw2();

    h_RecoToppTTTRestFrame = store(new TH1D ( "RecoToppTTTRestFrame", "Top pT in TTBar Rest Frame",1000,0,1000 ));
    h_HypToppTTTRestFrame = store(new TH1D ( "HypToppTTTRestFrame", "Top pT in TTBar Rest Frame (HYP)",1000,0,1000 ));
    h_VisGenToppTTTRestFrame = store(new TH1D ( "VisGenToppTTTRestFrame", "Top pT in TTBar Rest Frame (VisGEN)",1000,0,1000 ));
    h_GenRecoToppTTTRestFrame = store(new TH2D ( "GenRecoToppTTTRestFrame", "Gen/Reco (Top pT in TTBar Rest Frame)",1000,0,1000,1000,0,1000));

    h_RecoAntiToppTTTRestFrame = store(new TH1D ( "RecoAntiToppTTTRestFrame", "AntiTop pT in TTBar Rest Frame",1000,0,1000 ));
    h_HypAntiToppTTTRestFrame = store(new TH1D ( "HypAntiToppTTTRestFrame", "AntiTop pT in TTBar Rest Frame (HYP)",1000,0,1000 ));
    h_VisGenAntiToppTTTRestFrame = store(new TH1D ( "VisGenAntiToppTTTRestFrame", "AntiTop pT in TTBar Rest Frame (VisGEN)",1000,0,1000 ));
    h_GenRecoAntiToppTTTRestFrame = store(new TH2D ( "GenRecoAntiToppTTTRestFrame", "Gen/Reco (AntiTop pT in TTBar Rest Frame)",1000,0,1000,1000,0,1000));

    h_GenRecoLeptonpT = store(new TH2D ( "GenRecoLeptonpT", "Gen/Reco Matching", 1200,0,1200, 1200,0,1200 ));
    h_GenRecoAntiLeptonpT = store(new TH2D ( "GenRecoAntiLeptonpT", "Gen/Reco Matching", 1200,0,1200, 1200,0,1200 ));
    h_HypLeptonpT = store(new TH1D ( "HypLeptonpT", "Lepton Hypothesis pT", 1200,0,1200 ));
    h_HypAntiLeptonpT = store(new TH1D ( "HypAntiLeptonpT", "AntiLepton Hypothesis pT", 1200,0,1200 ));
    h_VisGenLeptonpT = store(new TH1D ( "VisGenLeptonpT", "Lepton VisGenothesis pT", 1200,0,1200 ));
    h_VisGenAntiLeptonpT = store(new TH1D ( "VisGenAntiLeptonpT", "AntiLepton VisGenothesis pT", 1200,0,1200 ));
    h_RecoLeptonpT = store(new TH1D ( "RecoLeptonpT","pT of Lepton (HYP)",1200,0,1200 ));
    h_RecoAntiLeptonpT = store(new TH1D ( "RecoAntiLeptonpT","pT of AntiLepton (HYP)",1200,0,1200 ));

    h_GenRecoLeptonEta = store(new TH2D ( "GenRecoLeptonEta", "Gen/Reco Matching", 200,-5,5, 200,-5,5 ));
    h_GenRecoAntiLeptonEta = store(new TH2D ( "GenRecoAntiLeptonEta", "Gen/Reco Matching", 200,-5,5, 200,-5,5 ));
    h_HypLeptonEta = store(new TH1D ( "HypLeptonEta", "Lepton Eta", 200,-5,5 ));
    h_HypAntiLeptonEta = store(new TH1D ( "HypAntiLeptonEta", "AntiLepton Hypothesis Eta", 200,-5,5 ));
    h_VisGenLeptonEta = store(new TH1D ( "VisGenLeptonEta", "Lepton Eta", 200,-5,5 ));
    h_VisGenAntiLeptonEta = store(new TH1D ( "VisGenAntiLeptonEta", "AntiLepton VisGenothesis Eta", 200,-5,5 ));
    h_RecoLeptonEta = store(new TH1D ( "RecoLeptonEta","Eta of Lepton (HYP)",200,-5,5 ));
    h_RecoAntiLeptonEta = store(new TH1D ( "RecoAntiLeptonEta","Eta of AntiLepton (HYP)",200,-5,5 ));

    h_VisGenToppT = store(new TH1D ( "VisGenToppT", "Top pT (VisGen)", 1200,0,1200 ));
    h_VisGenTopEta = store(new TH1D ( "VisGenTopEta", "Top Eta (VisGen)", 100, -5, 5 ));

    h_VisGenAntiToppT = store(new TH1D ( "VisGenAntiToppT", "AntiTop pT (VisGen)", 1200,0,1200 ));
    h_VisGenAntiTopEta = store(new TH1D ( "VisGenAntiTopEta", "AntiTop pT (VisGen)", 100, -5, 5 ));

    h_VisGenBJetpT = store(new TH1D ( "VisGenBJetpT", "B VisGenothesis pT", 1200, 0, 1200 ));
    h_VisGenBJetEta = store(new TH1D ( "VisGenBJetEta", "B VisGenothesis Eta", 200, -5, 5 ));
    h_VisGenBJetRapidity = store(new TH1D ( "VisGenBJetRapidity", "B VisGenothesis Rapidity", 100, -5, 5 ));

    h_VisGenAntiBJetpT = store(new TH1D ( "VisGenAntiBJetpT", "AntiB VisGenothesis pT", 1200, 0, 1200 ));
    h_VisGenAntiBJetEta = store(new TH1D ( "VisGenAntiBJetEta", "AntiB VisGenothesis Eta", 200, -5, 5 ));
    h_VisGenAntiBJetRapidity = store(new TH1D ( "VisGenAntiBJetRapidity", "AntiB VisGenothesis Rapidity", 100, -5, 5 ));

    h_GenRecoBJetpT = store(new TH2D ( "GenRecoBJetpT", "Gen/Reco Matching", 1200, 0, 1200, 1200, 0, 1200 ));
    h_GenRecoBJetEta = store(new TH2D ( "GenRecoBJetEta", "Gen/Reco Matching", 200, -5, 5, 200, -5, 5 ));
    h_GenRecoBJetRapidity = store(new TH2D ( "GenRecoBJetRapidity", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 ));

    h_GenRecoAntiBJetpT = store(new TH2D ( "GenRecoAntiBJetpT", "Gen/Reco Matching", 1200, 0, 1200, 1200, 0, 1200 ));
    h_GenRecoAntiBJetEta = store(new TH2D ( "GenRecoAntiBJetEta", "Gen/Reco Matching", 200, -5, 5, 200, -5, 5 ));
    h_GenRecoAntiBJetRapidity = store(new TH2D ( "GenRecoAntiBJetRapidity", "Gen/Reco Matching", 100, -5, 5, 100, -5, 5 ));

    h_GenRecoTopRapidity = store(new TH2D ( "GenRecoTopRapidity", "Gen/Reco Matching", 400, -5, 5, 400, -5, 5 ));
    h_GenRecoAntiTopRapidity = store(new TH2D ( "GenRecoAntiTopRapidity", "Gen/Reco Matching", 400, -5, 5, 400, -5, 5 ));
    h_GenRecoTopRapidityAbs = store(new TH2D ( "GenRecoTopRapidityAbs", "Gen/Reco Matching", 200, 0, 5, 200, 0, 5 ));
    h_GenRecoAntiTopRapidityAbs = store(new TH2D ( "GenRecoAntiTopRapidityAbs", "Gen/Reco Matching", 200, 0, 5, 200, 0, 5 ));
    h_GenRecoToppT = store(new TH2D ( "GenRecoToppT", "Gen/Reco Matching", 1200,0,1200, 1200,0,1200 ));
    h_GenRecoAntiToppT = store(new TH2D ( "GenRecoAntiToppT", "Gen/Reco Matching", 1200,0,1200, 1200,0,1200 ));

    h_GenRecoMet = store(new TH2D("GenRecoMet", "Missing ET in the event", 500, 0, 500, 500, 0, 500));
    h_VisGenMet = store(new TH1D("VisGenMet", "MET (VisGEN)", 500, 0, 500));
    h_RecoMet = store(new TH1D("RecoMet", "Reconstructed MET", 500, 0, 500));
    h_HypMet = store(new TH1D("HypMet","MET ", 500, 0, 500));

    h_GenRecoHT = store(new TH2D("GenRecoHT", "HT in the event", 800, 0, 800, 800, 0, 800));
    h_VisGenHT = store(new TH1D("VisGenHT", "HT (VisGEN)", 800, 0, 800));
    h_RecoHT = store(new TH1D("RecoHT", "Reconstructed HT", 800, 0, 800));
    h_HypHT = store(new TH1D("HypHT", "HT", 800, 0, 800));

    h_GenRecoNeutrinopT = store(new TH2D("GenRecoNeutrinopT", "Gen/Reco nu pt", 80, 0, 400, 80, 0, 400));
    h_VisGenNeutrinopT = store(new TH1D("VisGenNeutrinopT", "Nu pT (VisGEN)", 80, 0, 400));
    h_RecoNeutrinopT = store(new TH1D("RecoNeutrinopT", "reco nu pT", 80, 0, 400));
    h_HypNeutrinopT = store(new TH1D("HypNeutrinopT", "hyp nu pT", 80, 0, 400));

    h_GenRecoAntiNeutrinopT = store(new TH2D("GenRecoAntiNeutrinopT", "Gen/Reco nubar pt", 80, 0, 400, 80, 0, 400));
    h_VisGenAntiNeutrinopT = store(new TH1D("VisGenAntiNeutrinopT", "Nubar pT (VisGEN)", 80, 0, 400));
    h_RecoAntiNeutrinopT = store(new TH1D("RecoAntiNeutrinopT", "reco nubar pT", 80, 0, 400));
    h_HypAntiNeutrinopT = store(new TH1D("HypAntiNeutrinopT", "hyp nubar pT", 80, 0, 400));

    h_GenRecoTTBarDeltaPhi = store( new TH2D("GenRecoTTBarDeltaPhi", "Gen/Reco #Delta#Phi (ttbar)", 3500, 0, 3.5, 3500, 0, 3.5));
    h_RecoTTBarDeltaPhi = store( new TH1D("RecoTTBarDeltaPhi", "#Delta#Phi ofTTBar (RECO)", 3500, 0, 3.5));
    h_HypTTBarDeltaPhi = store( new TH1D("HypTTBarDeltaPhi", "#Delta#Phi ofTTBar (HYP)", 3500, 0, 3.5));
    h_VisGenTTBarDeltaPhi = store( new TH1D("VisGenTTBarDeltaPhi", "#Delta#Phi ofTTBar (VisGen)", 3500, 0, 3.5));

    h_GenRecoTTBarDeltaRapidity = store( new TH2D("GenRecoTTBarDeltaRapidity", "Gen/Reco |y^{t}| - |y^{#bar{t}}|", 400, -5, 5, 400, -5, 5));
    h_RecoTTBarDeltaRapidity = store( new TH1D("RecoTTBarDeltaRapidity", "|y^{t}| - |y^{#bar{t}}| (RECO)", 400, -5, 5));
    h_HypTTBarDeltaRapidity = store( new TH1D("HypTTBarDeltaRapidity", "|y^{t}| - |y^{#bar{t}}| (HYP)", 400, -5, 5));
    h_VisGenTTBarDeltaRapidity = store( new TH1D("VisGenTTBarDeltaRapidity", "|y^{t}| - |y^{#bar{t}}| (VisGen)", 400, -5, 5));

    h_GenRecoBBBarpT = store( new TH2D("GenRecoBBBarpT", "Gen/Reco p_{T} (bbbar)", 400, 0, 800, 400, 0, 800));
    h_RecoBBBarpT = store( new TH1D("RecoBBBarpT", "p_{T} of bbbar (RECO)", 400, 0, 800));
    h_HypBBBarpT = store( new TH1D("HypBBBarpT", "p_{T} of bbbarpT (HYP)", 400, 0, 800));
    h_VisGenBBBarpT = store( new TH1D("VisGenBBBarpT", "p_{T} of bbbarpT (VisGen)", 400, 0, 800));

    h_GenRecoBBBarMass = store( new TH2D("GenRecoBBBarMass", "Gen/Reco Mass (bbbar)", 400, 0, 800, 400, 0, 800));
    h_RecoBBBarMass = store( new TH1D("RecoBBBarMass", "Mass of bbbar (RECO)", 400, 0, 800));
    h_HypBBBarMass = store( new TH1D("HypBBBarMass", "Mass of bbbarMass (HYP)", 400, 0, 800));
    h_VisGenBBBarMass = store( new TH1D("VisGenBBBarMass", "Mass of bbbarMass (VisGen)", 400, 0, 800));

    h_GenRecoBBBarDPhi = store(new TH2D ( "GenRecoBBBarDPhi", "Gen/Reco Matching (bbbar)", 320, 0., 3.2, 320, 0., 3.2 ));
    h_HypBBBarDPhi = store(new TH1D ( "HypBBBarDPhi", "#Delta#phi(bbbar) (HYP)",320, 0., 3.2 ));
    h_VisGenBBBarDPhi = store(new TH1D ( "VisGenBBBarDPhi", "#Delta #Phi (bbbar) (VisGEN)", 320, 0., 3.2 ));
    h_RecoBBBarDPhi = store(new TH1D ( "RecoBBBarDPhi", "#Delta #Phi (bbbar) (Reco)", 320, 0., 3.2 ));

    h_GenRecoTTBarRapidity = store(new TH2D ( "GenRecoTTBarRapidity", "Rapidity of TTbar System (HYP)", 400, -5, 5, 400, -5, 5 ));
    h_GenRecoTTBarRapidityAbs = store(new TH2D ( "GenRecoTTBarRapidityAbs", "Absolute Rapidity of TTbar System (HYP)", 200, 0, 5, 200, 0, 5 ));
    h_GenRecoTTBarpT = store(new TH2D ( "GenRecoTTBarpT", "pT of TTbar System (HYP)", 1200, 0, 1200, 1200, 0, 1200 ));
    h_GenRecoTTBarMass = store(new TH2D ( "GenRecoTTBarMass", "Mass of TTbar System (HYP)", 3000, 0, 3000, 3000, 0, 3000 ));
    h_GenRecoLLBarMass = store(new TH2D ( "GenRecoLLBarMass", "Mass of LLbar System (HYP)", 500, 0, 1000, 500, 0, 1000 ));
    h_GenRecoLLBarpT = store(new TH2D ( "GenRecoLLBarpT", "pT of LLbar System (HYP)", 200, 0, 1000, 200, 0, 1000 ));

    h_NJetMatching = store(new TH1D ( "NJetMatching", "NJet Gen/Reco Matching", 5, 0, 5 ));

    h_GenRecoLLBarDPhi = store(new TH2D ( "GenRecoLLBarDPhi", "Gen/Reco Matching", 3200, 0., 3.2, 3200, 0., 3.2 ));
    h_GenRecoLeptonantiBjetMass = store(new TH2D ( "GenRecoLeptonBjetMass", "Gen/Reco Matching", 500, 0, 1000, 500, 0, 1000 ));
    h_GenRecoAntiLeptonBjetMass = store(new TH2D ( "GenRecoAntiLeptonBjetMass", "Gen/Reco Matching", 500, 0, 1000, 500, 0, 1000 ));
    h_GenRecoJetMult = store(new TH2D ( "GenRecoJetMult", "Gen/Reco Matching", 26, -0.5, 25.5, 26, -0.5, 25.5 ));

    h_VisGenLLBarDPhi = store(new TH1D ( "VisGenLLBarDPhi", "#Delta #Phi (Lep, AntiLep) (VisGEN)", 3200, 0., 3.2 ));
    h_VisGenLeptonantiBjetMass = store(new TH1D ( "VisGenLeptonBjetMass", "M(Lep, AntiBJet) (VisGEN)", 500, 0, 1000 ));
    h_VisGenAntiLeptonBjetMass = store(new TH1D ( "VisGenAntiLeptonBjetMass", "M(AntiLep, BJet) (VisGEN)", 500, 0, 1000 ));
    h_VisGenJetMult = store(new TH1D ( "VisGenJetMult", "Jet Multiplicty (VisGEN)", 26, -0.5, 25.5 ));

    h_RecoLLBarDPhi = store(new TH1D ( "RecoLLBarDPhi", "#Delta #Phi (Lep, AntiLep) (Reco)", 3200, 0., 3.2 ));
    h_RecoLeptonantiBjetMass = store(new TH1D ( "RecoLeptonBjetMass", "M(Lep, AntiBJet) (Reco)", 500, 0, 1000 ));
    h_RecoAntiLeptonBjetMass = store(new TH1D ( "RecoAntiLeptonBjetMass", "M(AntiLep, BJet) (Reco)", 500, 0, 1000 ));
    h_RecoJetMult = store(new TH1D ( "RecoJetMult", "Jet Multiplicty (Reco)", 26, -0.5, 25.5 ));

    h_GenRecoLLBarDEta = store( new TH2D("GenRecoLLBarDEta", "Gen/Reco #Delta #Eta (Lep, AntiLep)", 400, -5, 5, 400, -5, 5));
    h_RecoLLBarDEta = store( new TH1D("RecoLLBarDEta", "#Delta #Eta (Lep, AntiLep) (Reco)", 400, -5, 5));
    h_HypLLBarDEta = store( new TH1D("HypLLBarDEta", "#Delta #Eta (Lep, AntiLep) (Hyp)", 400, -5, 5));
    h_VisGenLLBarDEta = store( new TH1D("VisGenLLBarDEta", "#Delta #Eta (Lep, AntiLep) (VisGen)", 400, -5, 5));

    h_GenRecoLLBarCosPhiInTopRestFrame = store(new TH2D ( "GenRecoLLBarCosPhiIntopRestFrame", "Gen/Reco Matching", 1000, -1., 1., 1000, -1., 1. ));
    h_HypLLBarCosPhiInTopRestFrame = store(new TH1D ( "HypLLBarCosPhiIntopRestFrame", "Cos #phi(Lep, AntiLep) in RF (HYP)",1000, -1., 1. ));
    h_VisGenLLBarCosPhiInTopRestFrame = store(new TH1D ( "VisGenLLBarCosPhiIntopRestFrame", "Cos #Phi (Lep, AntiLep) (VisGEN)", 1000, -1., 1. ));
    h_RecoLLBarCosPhiInTopRestFrame = store(new TH1D ( "RecoLLBarCosPhiIntopRestFrame", "Cos #Phi (Lep, AntiLep) (Reco)", 1000, -1., 1. ));

    h_HypToppTLead = store(new TH1D ( "HypToppTLead","Leading pT Top pT",1200,0,1200 ));
    h_RecoToppTLead = store(new TH1D ( "RecoToppTLead","Leading pT Top pT",1200,0,1200 ));
    h_VisGenToppTLead = store(new TH1D ( "VisGenToppTLead","Leading pT Top pT",1200,0,1200 ));
    h_GenRecoToppTLead = store(new TH2D ( "GenRecoToppTLead", "Gen/Reco Matching", 1200,0,1200,1200,0,1200 ));

    h_HypToppTNLead = store(new TH1D ( "HypToppTNLead","NLeading pT Top pT",1200,0,1200 ));
    h_RecoToppTNLead = store(new TH1D ( "RecoToppTNLead","NLeading pT Top pT",1200,0,1200 ));
    h_VisGenToppTNLead = store(new TH1D ( "VisGenToppTNLead","NLeading pT Top pT",1200,0,1200 ));
    h_GenRecoToppTNLead = store(new TH2D ( "GenRecoToppTNLead", "Gen/Reco Matching", 1200,0,1200,1200,0,1200 ));

    h_HypTopRapidityLead = store(new TH1D ( "HypTopRapidityLead","Leading pT Top Rapidity",400,-5,5 ));
    h_RecoTopRapidityLead = store(new TH1D ( "RecoTopRapidityLead","Leading pT Top Rapidity",400,-5,5 ));
    h_VisGenTopRapidityLead = store(new TH1D ( "VisGenTopRapidityLead","Leading pT Top Rapidity",400,-5,5 ));
    h_GenRecoTopRapidityLead = store(new TH2D ( "GenRecoTopRapidityLead", "Gen/Reco Matching", 400,-5,5,400,-5,5 ));

    h_HypTopRapidityNLead = store(new TH1D ( "HypTopRapidityNLead","NLeading pT Top Rapidity",400,-5,5 ));
    h_RecoTopRapidityNLead = store(new TH1D ( "RecoTopRapidityNLead","NLeading pT Top Rapidity",400,-5,5 ));
    h_VisGenTopRapidityNLead = store(new TH1D ( "VisGenTopRapidityNLead","NLeading pT Top Rapidity",400,-5,5 ));
    h_GenRecoTopRapidityNLead = store(new TH2D ( "GenRecoTopRapidityNLead", "Gen/Reco Matching", 400,-5,5,400,-5,5 ));

    h_HypTopMassLead = store(new TH1D ( "HypTopMassLead","Leading pT Top Mass",80,0,400 ));
    h_RecoTopMassLead = store(new TH1D ( "RecoTopMassLead","Leading pT Top Mass",80,0,400 ));
    h_VisGenTopMassLead = store(new TH1D ( "VisGenTopMassLead","Leading pT Top Mass",80,0,400 ));
    h_GenRecoTopMassLead = store(new TH2D ( "GenRecoTopMassLead", "Gen/Reco Matching", 80,0,400,80,0,400 ));

    h_HypTopMassNLead = store(new TH1D ( "HypTopMassNLead","NLeading pT Top Mass",80,0,400 ));
    h_RecoTopMassNLead = store(new TH1D ( "RecoTopMassNLead","NLeading pT Top Mass",80,0,400 ));
    h_VisGenTopMassNLead = store(new TH1D ( "VisGenTopMassNLead","NLeading pT Top Mass",80,0,400 ));
    h_GenRecoTopMassNLead = store(new TH2D ( "GenRecoTopMassNLead", "Gen/Reco Matching", 80,0,400,80,0,400 ));

    h_HypLeptonpTLead = store(new TH1D ( "HypLeptonpTLead","Leading pT Lepton pT",1200,0,1200 ));
    h_RecoLeptonpTLead = store(new TH1D ( "RecoLeptonpTLead","Leading pT Lepton pT",1200,0,1200 ));
    h_VisGenLeptonpTLead = store(new TH1D ( "VisGenLeptonpTLead","Leading pT Lepton pT",1200,0,1200 ));
    h_GenRecoLeptonpTLead = store(new TH2D ( "GenRecoLeptonpTLead", "Gen/Reco Matching", 1200,0,1200,1200,0,1200 ));

    h_HypLeptonpTNLead = store(new TH1D ( "HypLeptonpTNLead","NLeading pT Lepton pT",1200,0,1200 ));
    h_RecoLeptonpTNLead = store(new TH1D ( "RecoLeptonpTNLead","NLeading pT Lepton pT",1200,0,1200 ));
    h_VisGenLeptonpTNLead = store(new TH1D ( "VisGenLeptonpTNLead","NLeading pT Lepton pT",1200,0,1200 ));
    h_GenRecoLeptonpTNLead = store(new TH2D ( "GenRecoLeptonpTNLead", "Gen/Reco Matching", 1200,0,1200,1200,0,1200 ));

    h_HypLeptonEtaLead = store(new TH1D ( "HypLeptonEtaLead","Leading pT Lepton Eta",200,-5,5 ));
    h_RecoLeptonEtaLead = store(new TH1D ( "RecoLeptonEtaLead","Leading pT Lepton Eta",200,-5,5 ));
    h_VisGenLeptonEtaLead = store(new TH1D ( "VisGenLeptonEtaLead","Leading pT Lepton Eta",200,-5,5 ));
    h_GenRecoLeptonEtaLead = store(new TH2D ( "GenRecoLeptonEtaLead", "Gen/Reco Matching", 200,-5,5,200,-5,5 ));

    h_HypLeptonEtaNLead = store(new TH1D ( "HypLeptonEtaNLead","NLeading pT Lepton Eta",200,-5,5 ));
    h_RecoLeptonEtaNLead = store(new TH1D ( "RecoLeptonEtaNLead","NLeading pT Lepton Eta",200,-5,5 ));
    h_VisGenLeptonEtaNLead = store(new TH1D ( "VisGenLeptonEtaNLead","NLeading pT Lepton Eta",200,-5,5 ));
    h_GenRecoLeptonEtaNLead = store(new TH2D ( "GenRecoLeptonEtaNLead", "Gen/Reco Matching", 200,-5,5,200,-5,5 ));

    h_HypBJetpTLead = store(new TH1D ( "HypBJetpTLead","Leading pT BJet pT",1200,0,1200 ));
    h_RecoBJetpTLead = store(new TH1D ( "RecoBJetpTLead","Leading pT BJet pT",1200,0,1200 ));
    h_VisGenBJetpTLead = store(new TH1D ( "VisGenBJetpTLead","Leading pT BJet pT",1200,0,1200 ));
    h_GenRecoBJetpTLead = store(new TH2D ( "GenRecoBJetpTLead", "Gen/Reco Matching", 1200,0,1200,1200,0,1200 ));

    h_HypBJetpTNLead = store(new TH1D ( "HypBJetpTNLead","NLeading pT BJet pT",1200,0,1200 ));
    h_RecoBJetpTNLead = store(new TH1D ( "RecoBJetpTNLead","NLeading pT BJet pT",1200,0,1200 ));
    h_VisGenBJetpTNLead = store(new TH1D ( "VisGenBJetpTNLead","NLeading pT BJet pT",1200,0,1200 ));
    h_GenRecoBJetpTNLead = store(new TH2D ( "GenRecoBJetpTNLead", "Gen/Reco Matching", 1200,0,1200,1200,0,1200 ));

    h_HypBJetEtaLead = store(new TH1D ( "HypBJetEtaLead","Leading pT BJet Eta",200,-5,5 ));
    h_RecoBJetEtaLead = store(new TH1D ( "RecoBJetEtaLead","Leading pT BJet Eta",200,-5,5 ));
    h_VisGenBJetEtaLead = store(new TH1D ( "VisGenBJetEtaLead","Leading pT BJet Eta",200,-5,5 ));
    h_GenRecoBJetEtaLead = store(new TH2D ( "GenRecoBJetEtaLead", "Gen/Reco Matching", 200,-5,5,200,-5,5 ));

    h_HypBJetEtaNLead = store(new TH1D ( "HypBJetEtaNLead","NLeading pT BJet Eta",200,-5,5 ));
    h_RecoBJetEtaNLead = store(new TH1D ( "RecoBJetEtaNLead","NLeading pT BJet Eta",200,-5,5 ));
    h_VisGenBJetEtaNLead = store(new TH1D ( "VisGenBJetEtaNLead","NLeading pT BJet Eta",200,-5,5 ));
    h_GenRecoBJetEtaNLead = store(new TH2D ( "GenRecoBJetEtaNLead", "Gen/Reco Matching", 200,-5,5,200,-5,5 ));
    //h_RecoExtraJetpT   = store(new TH1D("RecoExtraJetpT","pT of additional Jet (HYP)",400,0,400));
    h_RecoExtraJetEta  = store(new TH1D("RecoExtraJetEta","#eta of additional Jet (HYP)",100,-5,5));
    //h_RecoExtraJetpT2  = store(new TH1D("RecoExtraJetpT2","pT of additional Jet (HYP)",400,0,400));
    h_RecoExtraJetEta2 = store(new TH1D("RecoExtraJetEta2","#eta of additional Jet (HYP)",100,-5,5));
    //h_RecoExtraJetpT3  = store(new TH1D("RecoExtraJetpT3","pT of additional Jet (HYP)",400,0,400));
    h_RecoExtraJetEta3 = store(new TH1D("RecoExtraJetEta3","#eta of additional Jet (HYP)",100,-5,5));
    //h_RecoExtraJetpT4  = store(new TH1D("RecoExtraJetpT4","pT of additional Jet (HYP)",400,0,400));
    h_RecoExtraJetEta4 = store(new TH1D("RecoExtraJetEta4","#eta of additional Jet (HYP)",100,-5,5));

    //h_HypExtraJetpT   = store(new TH1D("HypExtraJetpT","pT of additional Jet",400,0,400));
    h_HypExtraJetEta  = store(new TH1D("HypExtraJetEta","#eta of additional Jet",100,-5,5));
    //h_HypExtraJetpT2  = store(new TH1D("HypExtraJetpT2","pT of additional Jet",400,0,400));
    h_HypExtraJetEta2 = store(new TH1D("HypExtraJetEta2","#eta of additional Jet",100,-5,5));
    //h_HypExtraJetpT3  = store(new TH1D("HypExtraJetpT3","pT of additional Jet",400,0,400));
    h_HypExtraJetEta3 = store(new TH1D("HypExtraJetEta3","#eta of additional Jet",100,-5,5));
    //h_HypExtraJetpT4  = store(new TH1D("HypExtraJetpT4","pT of additional Jet",400,0,400));
    h_HypExtraJetEta4 = store(new TH1D("HypExtraJetEta4","#eta of additional Jet",100,-5,5));

    h_VisGenExtraJetpT   = store(new TH1D("VisGenExtraJetpT","pT of gen additional Jet",400,0,400));
    h_VisGenExtraJetEta  = store(new TH1D("VisGenExtraJetEta","#eta of gen additional Jet",100,-5,5));
    h_VisGenExtraJetpT2  = store(new TH1D("VisGenExtraJetpT2","pT of gen additional Jet",400,0,400));
    h_VisGenExtraJetEta2 = store(new TH1D("VisGenExtraJetEta2","#eta of gen additional Jet",100,-5,5));
    h_VisGenExtraJetpT3  = store(new TH1D("VisGenExtraJetpT3","pT of gen additional Jet",400,0,400));
    h_VisGenExtraJetEta3 = store(new TH1D("VisGenExtraJetEta3","#eta of gen additional Jet",100,-5,5));
    h_VisGenExtraJetpT4  = store(new TH1D("VisGenExtraJetpT4","pT of gen additional Jet",400,0,400));
    h_VisGenExtraJetEta4 = store(new TH1D("VisGenExtraJetEta4","#eta of gen additional Jet",100,-5,5));
    h_GenRecoExtraJetpT   = store(new TH2D("GenRecoExtraJetpT","Gen/Reco pT of additional Jet",400,0,400,400,0,400));
    h_GenRecoExtraJetEta  = store(new TH2D("GenRecoExtraJetEta","Gen/Reco #eta of additional Jet",100,-5,5,100,-5,5));
    h_GenRecoExtraJetpT2  = store(new TH2D("GenRecoExtraJetpT2","Gen/Reco pT of additional Jet",400,0,400,400,0,400));
    h_GenRecoExtraJetEta2 = store(new TH2D("GenRecoExtraJetEta2","Gen/Reco #eta of additional Jet",100,-5,5,100,-5,5));
    h_GenRecoExtraJetpT3  = store(new TH2D("GenRecoExtraJetpT3","Gen/Reco pT of additional Jet",400,0,400,400,0,400));
    h_GenRecoExtraJetEta3 = store(new TH2D("GenRecoExtraJetEta3","Gen/Reco #eta of additional Jet",100,-5,5,100,-5,5));
    h_GenRecoExtraJetpT4  = store(new TH2D("GenRecoExtraJetpT4","Gen/Reco pT of additional Jet",400,0,400,400,0,400));
    h_GenRecoExtraJetEta4 = store(new TH2D("GenRecoExtraJetEta4","Gen/Reco #eta of additional Jet",100,-5,5,100,-5,5));

    h_HypExtraJetAbsEta  = store(new TH1D("HypExtraJetAbsEta","#eta of additional Jet",100,0.,5));
    h_HypExtraJetAbsEta2  = store(new TH1D("HypExtraJetAbsEta2","#eta of additional Jet",100,0.,5));
    h_HypExtraJetAbsEta3  = store(new TH1D("HypExtraJetAbsEta3","#eta of additional Jet",100,0.,5));
    h_HypExtraJetAbsEta4  = store(new TH1D("HypExtraJetAbsEta4","#eta of additional Jet",100,0.,5));
    h_VisGenExtraJetAbsEta  = store(new TH1D("VisGenExtraJetAbsEta","#eta of additional Jet",100,0.,5));
    h_VisGenExtraJetAbsEta2  = store(new TH1D("VisGenExtraJetAbsEta2","#eta of additional Jet",100,0.,5));
    h_VisGenExtraJetAbsEta3  = store(new TH1D("VisGenExtraJetAbsEta3","#eta of additional Jet",100,0.,5));
    h_VisGenExtraJetAbsEta4  = store(new TH1D("VisGenExtraJetAbsEta4","#eta of additional Jet",100,0.,5));
    h_GenRecoExtraJetAbsEta  = store(new TH2D("GenRecoExtraJetAbsEta","Gen/Reco #eta of additional Jet",100,0.,5,100,0.,5));
    h_GenRecoExtraJetAbsEta2  = store(new TH2D("GenRecoExtraJetAbsEta2","Gen/Reco #eta of additional Jet",100,0.,5,100,0.,5));
    h_GenRecoExtraJetAbsEta3  = store(new TH2D("GenRecoExtraJetAbsEta3","Gen/Reco #eta of additional Jet",100,0.,5,100,0.,5));
    h_GenRecoExtraJetAbsEta4  = store(new TH2D("GenRecoExtraJetAbsEta4","Gen/Reco #eta of additional Jet",100,0.,5,100,0.,5));

    h_RecoExtraJetHT   = store(new TH1D("RecoExtraJetHT","HT of additional Jet (Reco)",160,0,800));
    h_HypExtraJetHT   = store(new TH1D("HypExtraJetHT","HT of additional Jet (HYP)",160,0,800));
    h_VisGenExtraJetHT   = store(new TH1D("VisGenExtraJetHT","HT of additional Jet (HYP)",160,0,800));
    h_GenRecoExtraJetHT   = store(new TH2D("GenRecoExtraJetHT","Gen/Reco HT of additional Jet",160,0,800,160,0,800));

    h_GenRecoJetMultpt30 = store(new TH2D("GenRecoJetMultpt30", "Gen/Reco Matching",8,1.5,9.5,8,1.5,9.5));
    h_RecoJetMultpt30    = store(new TH1D("RecoJetMultpt30", "Jet Multiplicity (HYP)",8,1.5,9.5));
    h_HypJetMultpt30     = store(new TH1D("HypJetMultpt30", "Jet Multiplicity (HYP)",8,1.5,9.5));
    h_VisGenJetMultpt30  = store(new TH1D("VisGenJetMultpt30", "Jet Multiplicty (VisGEN)",8,1.5,9.5));

    h_GenRecoJetMultpt40 = store(new TH2D("GenRecoJetMultpt40", "Gen/Reco Matching",10,-0.5,9.5,10,-0.5,9.5));
    h_RecoJetMultpt40    = store(new TH1D("RecoJetMultpt40", "Jet Multiplicity (HYP)",10,-0.5,9.5));
    h_HypJetMultpt40     = store(new TH1D("HypJetMultpt40", "Jet Multiplicity (HYP)",10,-0.5,9.5));
    h_VisGenJetMultpt40  = store(new TH1D("VisGenJetMultpt40", "Jet Multiplicty (VisGEN)",10,-0.5,9.5));

    h_GenRecoJetMultpt60 = store(new TH2D("GenRecoJetMultpt60", "Gen/Reco Matching",10,-0.5,9.5,10,-0.5,9.5));
    h_RecoJetMultpt60    = store(new TH1D("RecoJetMultpt60", "Jet Multiplicity (HYP)",10,-0.5,9.5));
    h_HypJetMultpt60     = store(new TH1D("HypJetMultpt60", "Jet Multiplicity (HYP)",10,-0.5,9.5));
    h_VisGenJetMultpt60  = store(new TH1D("VisGenJetMultpt60", "Jet Multiplicty (VisGEN)",10,-0.5,9.5));

    h_GenRecoJetMultpt100 = store(new TH2D("GenRecoJetMultpt100", "Gen/Reco Matching",10,-0.5,9.5,10,-0.5,9.5));
    h_RecoJetMultpt100    = store(new TH1D("RecoJetMultpt100", "Jet Multiplicity (HYP)",10,-0.5,9.5));
    h_HypJetMultpt100     = store(new TH1D("HypJetMultpt100", "Jet Multiplicity (HYP)",10,-0.5,9.5));
    h_VisGenJetMultpt100  = store(new TH1D("VisGenJetMultpt100", "Jet Multiplicty (VisGEN)",10,-0.5,9.5));

    h_GenRecoDeltaRExtraJet12 = store(new TH2D("GenRecoDeltaRExtraJet12","Gen/Reco #Delta R of additional Jets",100,0.,10,100,0,10));
    h_RecoDeltaRExtraJet12 = store(new TH1D("RecoDeltaRExtraJet12","Reco #Delta R of additional Jets",100,0,10));
    h_HypDeltaRExtraJet12 = store(new TH1D("HypDeltaRExtraJet12","Hyp #DeltaR of additional Jets",100,0,10));
    h_VisGenDeltaRExtraJet12 = store(new TH1D("VisGenDeltaRExtraJet12","#Delta R of additional Jets",100,0,10));

    h_GenRecoMassExtraJet12 = store(new TH2D("GenRecoMassExtraJet12","Gen/Reco #Delta R of additional Jets",120,0,600,120,0,600));
    h_RecoMassExtraJet12 = store(new TH1D("RecoMassExtraJet12","Reco #Delta R of additional Jets",120,0,600));
    h_HypMassExtraJet12 = store(new TH1D("HypMassExtraJet12","Hyp #DeltaR of additional Jets",120,0,600));
    h_VisGenMassExtraJet12 = store(new TH1D("VisGenMassExtraJet12","#Delta R of additional Jets",120,0,600));

    h_RecoDeltaPhiExtraJet12 = store(new TH1D("RecoDeltaPhiExtraJet12","Reco #Delta#Phi of additional Jets",100,-TMath::Pi(), TMath::Pi()));
    h_HypDeltaPhiExtraJet12 = store(new TH1D("HypDeltaPhiExtraJet12","Hyp #Delta#Phi of additional Jets",100,-TMath::Pi(), TMath::Pi()));
    h_GenRecoDeltaPhiExtraJet12 = store(new TH2D("GenRecoDeltaPhiExtraJet12","Gen/Reco #Delta #phi of additional Jets",100,-TMath::Pi(), TMath::Pi(),100,0,10));
    h_VisGenDeltaPhiExtraJet12 = store(new TH1D("VisGenDeltaPhiExtraJet12","#Delta #phi of additional Jets",100,-TMath::Pi(), TMath::Pi()));

    h_RecoPhiExtraJet12 = store(new TH1D("RecoPhiExtraJet12","Reco #sum#Phi of additional Jets",100,-TMath::Pi(), TMath::Pi()));
    h_HypPhiExtraJet12 = store(new TH1D("HypPhiExtraJet12","Hyp #sum#Phi of additional Jets",100,-TMath::Pi(), TMath::Pi()));
    h_GenRecoPhiExtraJet12 = store(new TH2D("GenRecoPhiExtraJet12","Gen/Reco #Sum #phi of additional Jets",100,-TMath::Pi(), TMath::Pi(),100,-TMath::Pi(), TMath::Pi()));
    h_VisGenPhiExtraJet12 = store(new TH1D("VisGenPhiExtraJet12","#Sum #phi of additional Jets",100, -TMath::Pi(), TMath::Pi()));

    h_VisGenTTBar1stJetMass = store(new TH1D("VisGenTTBar1stJetMass","TTBar1stJetMass (VisGEN)",400,0,1)); h_VisGenTTBar1stJetMass->Sumw2();
    h_GenRecoTTBar1stJetMass = store(new TH2D("GenRecoTTBar1stJetMass","TTBar1stJetMass Gen/Reco",400,0,1,400,0,1)); h_GenRecoTTBar1stJetMass->Sumw2();
    h_RecoTTBar1stJetMass = store(new TH1D("RecoTTBar1stJetMass", "TTBar1stJetMass (Reco)",400,0,1)); h_RecoTTBar1stJetMass->Sumw2();
    h_HypTTBar1stJetMass = store(new TH1D("HypTTBar1stJetMass","TTBar1stJetMass (HYP)",400,0,1)); h_HypTTBar1stJetMass->Sumw2();

    h_VisGenTTBar0Mass = store(new TH1D("VisGenTTBar0Mass","TTBar0Mass (VisGEN)",200,0,2));
    h_GenRecoTTBar0Mass = store(new TH2D("GenRecoTTBar0Mass","TTBar0Mass Gen/Reco",200,0,2,200,0,2));
    h_RecoTTBar0Mass = store(new TH1D("RecoTTBar0Mass", "TTBar0Mass (Reco)",200,0,2));
    h_HypTTBar0Mass = store(new TH1D("HypTTBar0Mass","TTBar0Mass (HYP)",200,0,2));

    h_HypTopPartonFraction = store(new TH1D("HypTopPartonFraction","Parton Momentum Fraction (HYP)",100,0,1));
    h_VisGenTopPartonFraction = store(new TH1D("VisGenTopPartonFraction","Parton Momentum Fraction (VisGEN)",100,0,1));
    h_RecoTopPartonFraction = store(new TH1D("RecoTopPartonFraction","Parton Momentum Fraction (reco)",100,0,1));
    h_GenRecoTopPartonFraction = store(new TH2D("GenRecoTopPartonFraction","Parton Momentum Fraction (Gen/Reco)",100,0,1, 100, 0,1));
    h_HypAntiTopPartonFraction = store(new TH1D("HypAntiTopPartonFraction","AntiParton Momentum Fraction (HYP)",100,0,1));
    h_VisGenAntiTopPartonFraction = store(new TH1D("VisGenAntiTopPartonFraction","AntiParton Momentum Fraction (VisGEN)",100,0,1));
    h_RecoAntiTopPartonFraction = store(new TH1D("RecoAntiTopPartonFraction","AntiParton Momentum Fraction (reco)",100,0,1));
    h_GenRecoAntiTopPartonFraction = store(new TH2D("GenRecoAntiTopPartonFraction","AntiParton Momentum Fraction (Gen/Reco)",100,0,1, 100, 0,1));

    h_HypExtraJetpT   = store(new TH1D("HypExtraJetpT","pT of additional Jet",80,0,400));
    h_HypExtraJetpT2  = store(new TH1D("HypExtraJetpT2","pT of additional Jet",80,0,400));
    h_HypExtraJetpT3  = store(new TH1D("HypExtraJetpT3","pT of additional Jet",80,0,400));
    h_HypExtraJetpT4  = store(new TH1D("HypExtraJetpT4","pT of additional Jet",80,0,400));
    h_RecoExtraJetpT   = store(new TH1D("RecoExtraJetpT","pT of additional Jet (HYP)",80,0,400));
    h_RecoExtraJetpT2  = store(new TH1D("RecoExtraJetpT2","pT of additional Jet (HYP)",80,0,400));
    h_RecoExtraJetpT3  = store(new TH1D("RecoExtraJetpT3","pT of additional Jet (HYP)",80,0,400));
    h_RecoExtraJetpT4  = store(new TH1D("RecoExtraJetpT4","pT of additional Jet (HYP)",80,0,400));

    h_RecoAbsDeltaPhiExtraJet12 = store(new TH1D("RecoAbsDeltaPhiExtraJet12","Reco #Delta#Phi of additional Jets",10,0, TMath::Pi()));
    h_HypAbsDeltaPhiExtraJet12 = store(new TH1D("HypAbsDeltaPhiExtraJet12","Hyp #Delta#Phi of additional Jets",10,0, TMath::Pi()));
    h_RecoAbsDeltaPhiExtraJet12_eta1 = store(new TH1D("RecoAbsDeltaPhiExtraJet12_eta1","Reco #Delta#Phi of additional Jets",10,0, TMath::Pi()));
    h_HypAbsDeltaPhiExtraJet12_eta1 = store(new TH1D("HypAbsDeltaPhiExtraJet12_eta1","Hyp #Delta#Phi of additional Jets",10,0, TMath::Pi()));
    h_RecoAbsDeltaPhiExtraJet12_eta2 = store(new TH1D("RecoAbsDeltaPhiExtraJet12_eta2","Reco #Delta#Phi of additional Jets",10,0, TMath::Pi()));
    h_HypAbsDeltaPhiExtraJet12_eta2 = store(new TH1D("HypAbsDeltaPhiExtraJet12_eta2","Hyp #Delta#Phi of additional Jets",10,0, TMath::Pi()));
    h_RecoAbsDeltaPhiExtraJet12_eta3 = store(new TH1D("RecoAbsDeltaPhiExtraJet12_eta3","Reco #Delta#Phi of additional Jets",10,0, TMath::Pi()));
    h_HypAbsDeltaPhiExtraJet12_eta3 = store(new TH1D("HypAbsDeltaPhiExtraJet12_eta3","Hyp #Delta#Phi of additional Jets",10,0, TMath::Pi()));
    h_RecoAbsDeltaPhiExtraJet12_eta4 = store(new TH1D("RecoAbsDeltaPhiExtraJet12_eta4","Reco #Delta#Phi of additional Jets",10,0, TMath::Pi()));
    h_HypAbsDeltaPhiExtraJet12_eta4 = store(new TH1D("HypAbsDeltaPhiExtraJet12_eta4","Hyp #Delta#Phi of additional Jets",10,0, TMath::Pi()));
    h_RecoAbsDeltaEtaExtraJet12 = store(new TH1D("RecoAbsDeltaEtaExtraJet12","Reco #Delta#eta of additional Jets",5,0, 5));
    h_HypAbsDeltaEtaExtraJet12 = store(new TH1D("HypAbsDeltaEtaExtraJet12","Hyp #Delta#eta of additional Jets",5,0, 5));


  //=======================PSEUDO-TOP: BEGIN==========================//
    //Simple protection from non-supported pseudo-top steering parameter
    if (analysisConfig_.sampleComposition().pseudoTopMode_ < 0 || analysisConfig_.sampleComposition().pseudoTopMode_ > 3) {std::cout<< "Requested pseudo-top mode is not supported! Exit..." <<std::endl; exit(322);}

    // Configure pseudo-top histograms
    if (analysisConfig_.sampleComposition().pseudoTopMode_ > 0) {
        h_VisPseudoAll = store(new TH1D ( "VisPseudoAll", "All Visible Pseudo particles (IM)", 40, 0, 400 ));
        h_VisPseudoAll_noweight = store(new TH1D ( "VisPseudoAll_noweight", "All Visible Pseudo particles (IM)", 40, 0, 400 ));

        //Input distributions at particle level with pseudo-top definition

        h_VisPseudoToppT = store(new TH1D ( "VisPseudoToppT", "Top pT (VisPseudo)", 1200,0,1200 ));
        h_VisPseudoAntiToppT = store(new TH1D ( "VisPseudoAntiToppT", "AntiTop pT (VisPseudo)", 1200,0,1200 ));
        h_VisPseudoToppTLead = store(new TH1D ( "VisPseudoToppTLead","Leading pT Top pT",1200,0,1200 ));
        h_VisPseudoToppTNLead = store(new TH1D ( "VisPseudoToppTNLead","NLeading pT Top pT",1200,0,1200 ));
        h_VisPseudoToppTTTRestFrame = store(new TH1D ( "VisPseudoToppTTTRestFrame", "Top pT in TTBar Rest Frame (VisPSEUDO)",1000,0,1000 ));
        h_VisPseudoAntiToppTTTRestFrame = store(new TH1D ( "VisPseudoAntiToppTTTRestFrame", "AntiTop pT in TTBar Rest Frame (VisPSEUDO)",1000,0,1000 ));

        h_VisPseudoTopRapidity = store(new TH1D ( "VisPseudoTopRapidity", "Rapidity of Top(VisPSEUDO)", 400, -5, 5 ));
        h_VisPseudoAntiTopRapidity = store(new TH1D ( "VisPseudoAntiTopRapidity", "Rapidity of AntiTop(VisPSEUDO)", 400, -5, 5 ));
        h_VisPseudoTopRapidityAbs = store(new TH1D ( "VisPseudoTopRapidityAbs", "Absolute Rapidity of Top(VisPSEUDO)", 200, 0, 5 ));
        h_VisPseudoAntiTopRapidityAbs = store(new TH1D ( "VisPseudoAntiTopRapidityAbs", "Absolute Rapidity of AntiTop(VisPSEUDO)", 200, 0, 5 ));
        h_VisPseudoTopRapidityLead = store(new TH1D ( "VisPseudoTopRapidityLead","Leading pT Top Rapidity",400,-5,5 ));
        h_VisPseudoTopRapidityNLead = store(new TH1D ( "VisPseudoTopRapidityNLead","NLeading pT Top Rapidity",400,-5,5 ));

        h_VisPseudoTTBarMass = store(new TH1D ( "VisPseudoTTBarMass", "Mass of TTbar System(VisPSEUDO)", 3000, 0, 3000 ));
        h_VisPseudoTTBarRapidity = store(new TH1D ( "VisPseudoTTBarRapidity", "Rapidity of TTbar System(VisPSEUDO)", 400, -5, 5 ));
        h_VisPseudoTTBarRapidityAbs = store(new TH1D ( "VisPseudoTTBarRapidityAbs", "Absolute Rapidity of TTbar System(VisPSEUDO)", 200, 0, 5 ));
        h_VisPseudoTTBarpT = store(new TH1D ( "VisPseudoTTBarpT", "pT of TTbar System(VisPSEUDO)", 1200, 0, 1200 ));
        h_VisPseudoTTBarDeltaPhi = store( new TH1D("VisPseudoTTBarDeltaPhi", "#Delta#Phi ofTTBar (VisPseudo)", 3500, 0, 3.5));
        h_VisPseudoTTBarDeltaRapidity = store( new TH1D("VisPseudoTTBarDeltaRapidity", "|y^{t}| - |y^{#bar{t}}| (VisPseudo)", 400, -5, 5));

        h_VisPseudoJetMultpt30  = store(new TH1D("VisPseudoJetMultpt30", "Jet Multiplicty (VisPSEUDO)",8,1.5,9.5));
        h_VisPseudoJetMultpt40  = store(new TH1D("VisPseudoJetMultpt40", "Jet Multiplicty (VisPSEUDO)",10,-0.5,9.5));
        h_VisPseudoJetMultpt60  = store(new TH1D("VisPseudoJetMultpt60", "Jet Multiplicty (VisPSEUDO)",10,-0.5,9.5));
        h_VisPseudoJetMultpt100  = store(new TH1D("VisPseudoJetMultpt100", "Jet Multiplicty (VisPSEUDO)",10,-0.5,9.5));

        h_VisPseudoLeptonpT = store(new TH1D ( "VisPseudoLeptonpT", "Lepton VisPseudothesis pT", 1200,0,1200 ));
        h_VisPseudoAntiLeptonpT = store(new TH1D ( "VisPseudoAntiLeptonpT", "AntiLepton VisPseudothesis pT", 1200,0,1200 ));
        h_VisPseudoLeptonpTLead = store(new TH1D ( "VisPseudoLeptonpTLead","Leading pT Lepton pT", 1200,0,1200 ));
        h_VisPseudoLeptonpTNLead = store(new TH1D ( "VisPseudoLeptonpTNLead","NLeading pT Lepton pT", 1200,0,1200 ));
        h_VisPseudoLeptonEta = store(new TH1D ( "VisPseudoLeptonEta", "Lepton Eta", 200,-5,5 ));
        h_VisPseudoAntiLeptonEta = store(new TH1D ( "VisPseudoAntiLeptonEta", "AntiLepton VisPseudothesis Eta", 200,-5,5 ));
        h_VisPseudoLeptonEtaLead = store(new TH1D ( "VisPseudoLeptonEtaLead","Leading pT Lepton Eta",200,-5,5 ));
        h_VisPseudoLeptonEtaNLead = store(new TH1D ( "VisPseudoLeptonEtaNLead","NLeading pT Lepton Eta",200,-5,5 ));

        h_VisPseudoBJetpT = store(new TH1D ( "VisPseudoBJetpT", "B VisPseudothesis pT", 1200, 0, 1200 ));
        h_VisPseudoBJetEta = store(new TH1D ( "VisPseudoBJetEta", "B VisPseudothesis Eta", 200, -5, 5 ));
        h_VisPseudoBJetRapidity = store(new TH1D ( "VisPseudoBJetRapidity", "B VisPseudothesis Rapidity", 100, -5, 5 ));
        h_VisPseudoAntiBJetpT = store(new TH1D ( "VisPseudoAntiBJetpT", "AntiB VisPseudothesis pT", 1200, 0, 1200 ));
        h_VisPseudoAntiBJetEta = store(new TH1D ( "VisPseudoAntiBJetEta", "AntiB VisPseudothesis Eta", 200, -5, 5 ));
        h_VisPseudoAntiBJetRapidity = store(new TH1D ( "VisPseudoAntiBJetRapidity", "AntiB VisPseudothesis Rapidity", 100, -5, 5 ));
        h_VisPseudoBJetpTLead = store(new TH1D ( "VisPseudoBJetpTLead","Leading pT BJet pT",1200,0,1200 ));
        h_VisPseudoBJetEtaLead = store(new TH1D ( "VisPseudoBJetEtaLead","Leading pT BJet Eta",200,-5,5 ));
        h_VisPseudoBJetpTNLead = store(new TH1D ( "VisPseudoBJetpTNLead","NLeading pT BJet pT",1200,0,1200 ));
        h_VisPseudoBJetEtaNLead = store(new TH1D ( "VisPseudoBJetEtaNLead","NLeading pT BJet Eta",200,-5,5 ));

        h_VisPseudoNuNuBarpT = store(new TH1D ( "VisPseudoNuNuBarpT", "sum pT of two neutrinos (VisPSEUDO)", 200, 0, 1000 ));

        h_VisPseudoLLBarpT = store(new TH1D ( "VisPseudoLLBarpT", "pT of LLbar System(VisPSEUDO)", 200, 0, 1000 ));
        h_VisPseudoLLBarMass = store(new TH1D ( "VisPseudoLLBarMass", "Mass of LLbar System(VisPSEUDO)", 500, 0, 1000 ));
        h_VisPseudoLLBarDPhi = store(new TH1D ( "VisPseudoLLBarDPhi", "#Delta #Phi (Lep, AntiLep) (VisPSEUDO)", 3200, 0., 3.2 ));
        h_VisPseudoLLBarDEta = store(new TH1D ( "VisPseudoLLBarDEta", "#Delta #Eta (Lep, AntiLep) (VisPSEUDO)", 400, -5, 5));
        h_VisPseudoLLBarCosPhiInTopRestFrame = store(new TH1D ( "VisPseudoLLBarCosPhiIntopRestFrame", "Cos #Phi (Lep, AntiLep) (VisPSEUDO)", 1000, -1., 1. ));

        h_VisPseudoBBBarpT = store( new TH1D("VisPseudoBBBarpT", "p_{T} of bbbarpT (VisPseudo)", 400, 0, 800));
        h_VisPseudoBBBarMass = store( new TH1D("VisPseudoBBBarMass", "Mass of bbbarMass (VisPseudo)", 400, 0, 800));
        h_VisPseudoBBBarDPhi = store(new TH1D ( "VisPseudoBBBarDPhi", "#Delta #Phi (bbbar) (VisPSEUDO)", 320, 0., 3.2 ));

        h_VisPseudoLeptonantiBjetMass = store(new TH1D ( "VisPseudoLeptonBjetMass", "M(Lep, AntiBJet) (VisPSEUDO)", 500, 0, 1000 ));
        h_VisPseudoAntiLeptonBjetMass = store(new TH1D ( "VisPseudoAntiLeptonBjetMass", "M(AntiLep, BJet) (VisPSEUDO)", 500, 0, 1000 ));

        //Pseudo-Reco response matrices

        h_PseudoRecoToppT = store(new TH2D ( "PseudoRecoToppT", "Pseudo/Reco Matching", 1200,0,1200, 1200,0,1200 ));
        h_PseudoRecoAntiToppT = store(new TH2D ( "PseudoRecoAntiToppT", "Pseudo/Reco Matching", 1200,0,1200, 1200,0,1200 ));
        h_PseudoRecoToppTLead = store(new TH2D ( "PseudoRecoToppTLead", "Pseudo/Reco Matching", 1200,0,1200,1200,0,1200 ));
        h_PseudoRecoToppTNLead = store(new TH2D ( "PseudoRecoToppTNLead", "Pseudo/Reco Matching", 1200,0,1200,1200,0,1200 ));
        h_PseudoRecoToppTTTRestFrame = store(new TH2D ( "PseudoRecoToppTTTRestFrame", "Pseudo/Reco (Top pT in TTBar Rest Frame)",1000,0,1000,1000,0,1000));
        h_PseudoRecoAntiToppTTTRestFrame = store(new TH2D ( "PseudoRecoAntiToppTTTRestFrame", "Pseudo/Reco (AntiTop pT in TTBar Rest Frame)",1000,0,1000,1000,0,1000));

        h_PseudoRecoTopRapidity = store(new TH2D ( "PseudoRecoTopRapidity", "Pseudo/Reco Matching", 400, -5, 5, 400, -5, 5 ));
        h_PseudoRecoAntiTopRapidity = store(new TH2D ( "PseudoRecoAntiTopRapidity", "Pseudo/Reco Matching", 400, -5, 5, 400, -5, 5 ));
        h_PseudoRecoTopRapidityAbs = store(new TH2D ( "PseudoRecoTopRapidityAbs", "Pseudo/Reco Matching", 200, 0, 5, 200, 0, 5 ));
        h_PseudoRecoAntiTopRapidityAbs = store(new TH2D ( "PseudoRecoAntiTopRapidityAbs", "Pseudo/Reco Matching", 200, 0, 5, 200, 0, 5 ));
        h_PseudoRecoTopRapidityLead = store(new TH2D ( "PseudoRecoTopRapidityLead", "Pseudo/Reco Matching", 400,-5,5,400,-5,5 ));
        h_PseudoRecoTopRapidityNLead = store(new TH2D ( "PseudoRecoTopRapidityNLead", "Pseudo/Reco Matching", 400,-5,5,400,-5,5 ));

        h_PseudoRecoTTBarRapidity = store(new TH2D ( "PseudoRecoTTBarRapidity", "Rapidity of TTbar System (HYP)", 400, -5, 5, 400, -5, 5 ));
        h_PseudoRecoTTBarRapidityAbs = store(new TH2D ( "PseudoRecoTTBarRapidityAbs", "Absolute Rapidity of TTbar System (HYP)", 200, 0, 5, 200, 0, 5 ));
        h_PseudoRecoTTBarpT = store(new TH2D ( "PseudoRecoTTBarpT", "pT of TTbar System (HYP)", 1200, 0, 1200, 1200, 0, 1200 ));
        h_PseudoRecoTTBarDeltaPhi = store( new TH2D("PseudoRecoTTBarDeltaPhi", "Pseudo/Reco #Delta#Phi (ttbar)", 3500, 0, 3.5, 3500, 0, 3.5));
        h_PseudoRecoTTBarDeltaRapidity = store( new TH2D("PseudoRecoTTBarDeltaRapidity", "Pseudo/Reco |y^{t}| - |y^{#bar{t}}|", 400, -5, 5, 400, -5, 5));
        h_PseudoRecoTTBarMass = store(new TH2D ( "PseudoRecoTTBarMass", "Mass of TTbar System (HYP)", 3000, 0, 3000, 3000, 0, 3000 ));

        h_PseudoRecoJetMultpt30 = store(new TH2D("PseudoRecoJetMultpt30", "Pseudo/Reco Matching",8,1.5,9.5,8,1.5,9.5));
        h_PseudoRecoJetMultpt40 = store(new TH2D("PseudoRecoJetMultpt40", "Pseudo/Reco Matching",10,-0.5,9.5,10,-0.5,9.5));
        h_PseudoRecoJetMultpt60 = store(new TH2D("PseudoRecoJetMultpt60", "Pseudo/Reco Matching",10,-0.5,9.5,10,-0.5,9.5));
        h_PseudoRecoJetMultpt100 = store(new TH2D("PseudoRecoJetMultpt100", "Pseudo/Reco Matching",10,-0.5,9.5,10,-0.5,9.5));

        h_PseudoRecoLeptonpT = store(new TH2D ( "PseudoRecoLeptonpT", "Pseudo/Reco Matching", 1200,0,1200, 1200,0,1200 ));
        h_PseudoRecoAntiLeptonpT = store(new TH2D ( "PseudoRecoAntiLeptonpT", "Pseudo/Reco Matching", 1200,0,1200, 1200,0,1200 ));
        h_PseudoRecoLeptonpTLead = store(new TH2D ( "PseudoRecoLeptonpTLead", "Pseudo/Reco Matching", 1200,0,1200,1200,0,1200 ));
        h_PseudoRecoLeptonpTNLead = store(new TH2D ( "PseudoRecoLeptonpTNLead", "Pseudo/Reco Matching", 1200,0,1200,1200,0,1200 ));
        h_PseudoRecoLeptonEta = store(new TH2D ( "PseudoRecoLeptonEta", "Pseudo/Reco Matching", 200,-5,5, 200,-5,5 ));
        h_PseudoRecoAntiLeptonEta = store(new TH2D ( "PseudoRecoAntiLeptonEta", "Pseudo/Reco Matching", 200,-5,5, 200,-5,5 ));
        h_PseudoRecoLeptonEtaLead = store(new TH2D ( "PseudoRecoLeptonEtaLead", "Pseudo/Reco Matching", 200,-5,5,200,-5,5 ));
        h_PseudoRecoLeptonEtaNLead = store(new TH2D ( "PseudoRecoLeptonEtaNLead", "Pseudo/Reco Matching", 200,-5,5,200,-5,5 ));

        h_PseudoRecoBJetpT = store(new TH2D ( "PseudoRecoBJetpT", "Pseudo/Reco Matching", 1200, 0, 1200, 1200, 0, 1200 ));
        h_PseudoRecoBJetEta = store(new TH2D ( "PseudoRecoBJetEta", "Pseudo/Reco Matching", 200, -5, 5, 200, -5, 5 ));
        h_PseudoRecoBJetRapidity = store(new TH2D ( "PseudoRecoBJetRapidity", "Pseudo/Reco Matching", 100, -5, 5, 100, -5, 5 ));
        h_PseudoRecoAntiBJetpT = store(new TH2D ( "PseudoRecoAntiBJetpT", "Pseudo/Reco Matching", 1200, 0, 1200, 1200, 0, 1200 ));
        h_PseudoRecoAntiBJetEta = store(new TH2D ( "PseudoRecoAntiBJetEta", "Pseudo/Reco Matching", 200, -5, 5, 200, -5, 5 ));
        h_PseudoRecoAntiBJetRapidity = store(new TH2D ( "PseudoRecoAntiBJetRapidity", "Pseudo/Reco Matching", 100, -5, 5, 100, -5, 5 ));
        h_PseudoRecoBJetpTLead = store(new TH2D ( "PseudoRecoBJetpTLead", "Pseudo/Reco Matching", 1200,0,1200,1200,0,1200 ));
        h_PseudoRecoBJetEtaLead = store(new TH2D ( "PseudoRecoBJetEtaLead", "Pseudo/Reco Matching", 200,-5,5,200,-5,5 ));
        h_PseudoRecoBJetpTNLead = store(new TH2D ( "PseudoRecoBJetpTNLead", "Pseudo/Reco Matching", 1200,0,1200,1200,0,1200 ));
        h_PseudoRecoBJetEtaNLead = store(new TH2D ( "PseudoRecoBJetEtaNLead", "Pseudo/Reco Matching", 200,-5,5,200,-5,5 ));

        h_PseudoRecoLLBarMass = store(new TH2D ( "PseudoRecoLLBarMass", "Mass of LLbar System (HYP)", 500, 0, 1000, 500, 0, 1000 ));
        h_PseudoRecoLLBarpT = store(new TH2D ( "PseudoRecoLLBarpT", "pT of LLbar System (HYP)", 200, 0, 1000, 200, 0, 1000 ));
        h_PseudoRecoLLBarDPhi = store(new TH2D ( "PseudoRecoLLBarDPhi", "Pseudo/Reco Matching", 3200, 0., 3.2, 3200, 0., 3.2 ));
        h_PseudoRecoLLBarDEta = store( new TH2D( "PseudoRecoLLBarDEta", "Pseudo/Reco #Delta #Eta (Lep, AntiLep)", 400, -5, 5, 400, -5, 5));
        h_PseudoRecoLLBarCosPhiInTopRestFrame = store(new TH2D ( "PseudoRecoLLBarCosPhiIntopRestFrame", "Pseudo/Reco Matching", 1000, -1., 1., 1000, -1., 1. ));

        h_PseudoRecoBBBarpT = store( new TH2D("PseudoRecoBBBarpT", "Pseudo/Reco p_{T} (bbbar)", 400, 0, 800, 400, 0, 800));
        h_PseudoRecoBBBarMass = store( new TH2D("PseudoRecoBBBarMass", "Pseudo/Reco Mass (bbbar)", 400, 0, 800, 400, 0, 800));
        h_PseudoRecoBBBarDPhi = store(new TH2D ( "PseudoRecoBBBarDPhi", "Pseudo/Reco Matching (bbbar)", 320, 0., 3.2, 320, 0., 3.2 ));

        h_PseudoRecoLeptonantiBjetMass = store(new TH2D ( "PseudoRecoLeptonBjetMass", "Pseudo/Reco Matching", 500, 0, 1000, 500, 0, 1000 ));
        h_PseudoRecoAntiLeptonBjetMass = store(new TH2D ( "PseudoRecoAntiLeptonBjetMass", "Pseudo/Reco Matching", 500, 0, 1000, 500, 0, 1000 ));

        //Histograms for the estimation of out-of-space background at reco level

        h_HypToppT_OutOfSpace = store(new TH1D ( "HypToppT_OutOfSpace", "Out Of Space Bkgd", 1200,0,1200 ));
        h_HypAntiToppT_OutOfSpace = store(new TH1D ( "HypAntiToppT_OutOfSpace", "Out Of Space Bkgd", 1200,0,1200 ));
        h_HypToppTLead_OutOfSpace = store(new TH1D ( "HypToppTLead_OutOfSpace", "Out Of Space Bkgd", 1200,0,1200 ));
        h_HypToppTNLead_OutOfSpace = store(new TH1D ( "HypToppTNLead_OutOfSpace", "Out Of Space Bkgd", 1200,0,1200 ));
        h_HypToppTTTRestFrame_OutOfSpace = store(new TH1D ( "HypToppTTTRestFrame_OutOfSpace", "Out Of Space Bkgd", 1000,0,1000));
        h_HypAntiToppTTTRestFrame_OutOfSpace = store(new TH1D ( "HypAntiToppTTTRestFrame_OutOfSpace", "Out Of Space Bkgd", 1000,0,1000));

        h_HypTopRapidity_OutOfSpace = store(new TH1D ( "HypTopRapidity_OutOfSpace", "Out Of Space Bkgd", 400, -5, 5 ));
        h_HypAntiTopRapidity_OutOfSpace = store(new TH1D ( "HypAntiTopRapidity_OutOfSpace", "Out Of Space Bkgd", 400, -5, 5 ));
        h_HypTopRapidityAbs_OutOfSpace = store(new TH1D ( "HypTopRapidityAbs_OutOfSpace", "Out Of Space Bkgd", 200, 0, 5 ));
        h_HypAntiTopRapidityAbs_OutOfSpace = store(new TH1D ( "HypAntiTopRapidityAbs_OutOfSpace", "Out Of Space Bkgd", 200, 0, 5 ));
        h_HypTopRapidityLead_OutOfSpace = store(new TH1D ( "HypTopRapidityLead_OutOfSpace", "Out Of Space Bkgd", 400,-5,5 ));
        h_HypTopRapidityNLead_OutOfSpace = store(new TH1D ( "HypTopRapidityNLead_OutOfSpace", "Out Of Space Bkgd", 400,-5,5 ));

        h_HypTTBarRapidity_OutOfSpace = store(new TH1D ( "HypTTBarRapidity_OutOfSpace", "Out Of Space Bkgd", 400, -5, 5 ));
        h_HypTTBarRapidityAbs_OutOfSpace = store(new TH1D ( "HypTTBarRapidityAbs_OutOfSpace", "Out Of Space Bkgd", 200, 0, 5 ));
        h_HypTTBarpT_OutOfSpace = store(new TH1D ( "HypTTBarpT_OutOfSpace", "Out Of Space Bkgd", 1200, 0, 1200 ));
        h_HypTTBarDeltaPhi_OutOfSpace = store( new TH1D("HypTTBarDeltaPhi_OutOfSpace", "Out Of Space Bkgd", 3500, 0, 3.5));
        h_HypTTBarDeltaRapidity_OutOfSpace = store( new TH1D("HypTTBarDeltaRapidity_OutOfSpace", "Out Of Space Bkgd", 400, -5, 5));
        h_HypTTBarMass_OutOfSpace = store(new TH1D ( "HypTTBarMass_OutOfSpace", "Out Of Space Bkgd", 3000, 0, 3000 ));

        h_HypJetMultpt30_OutOfSpace = store(new TH1D("HypJetMultpt30_OutOfSpace", "Out Of Space Bkgd",8,1.5,9.5));
        h_HypJetMultpt40_OutOfSpace = store(new TH1D("HypJetMultpt40_OutOfSpace", "Out Of Space Bkgd",10,-0.5,9.5));
        h_HypJetMultpt60_OutOfSpace = store(new TH1D("HypJetMultpt60_OutOfSpace", "Out Of Space Bkgd",10,-0.5,9.5));
        h_HypJetMultpt100_OutOfSpace = store(new TH1D("HypJetMultpt100_OutOfSpace", "Out Of Space Bkgd",10,-0.5,9.5));

        h_HypLeptonpT_OutOfSpace = store(new TH1D ( "HypLeptonpT_OutOfSpace", "Out Of Space Bkgd", 1200,0,1200 ));
        h_HypAntiLeptonpT_OutOfSpace = store(new TH1D ( "HypAntiLeptonpT_OutOfSpace", "Out Of Space Bkgd", 1200,0,1200 ));
        h_HypLeptonpTLead_OutOfSpace = store(new TH1D ( "HypLeptonpTLead_OutOfSpace", "Out Of Space Bkgd", 1200,0,1200 ));
        h_HypLeptonpTNLead_OutOfSpace = store(new TH1D ( "HypLeptonpTNLead_OutOfSpace", "Out Of Space Bkgd", 1200,0,1200 ));
        h_HypLeptonEta_OutOfSpace = store(new TH1D ( "HypLeptonEta_OutOfSpace", "Out Of Space Bkgd", 200,-5,5 ));
        h_HypAntiLeptonEta_OutOfSpace = store(new TH1D ( "HypAntiLeptonEta_OutOfSpace", "Out Of Space Bkgd", 200,-5,5 ));
        h_HypLeptonEtaLead_OutOfSpace = store(new TH1D ( "HypLeptonEtaLead_OutOfSpace", "Out Of Space Bkgd", 200,-5,5 ));
        h_HypLeptonEtaNLead_OutOfSpace = store(new TH1D ( "HypLeptonEtaNLead_OutOfSpace", "Out Of Space Bkgd", 200,-5,5 ));

        h_HypBJetpT_OutOfSpace = store(new TH1D ( "HypBJetpT_OutOfSpace", "Out Of Space Bkgd", 1200, 0, 1200 ));
        h_HypBJetEta_OutOfSpace = store(new TH1D ( "HypBJetEta_OutOfSpace", "Out Of Space Bkgd", 200, -5, 5 ));
        h_HypBJetRapidity_OutOfSpace = store(new TH1D ( "HypBJetRapidity_OutOfSpace", "Out Of Space Bkgd", 100, -5, 5 ));
        h_HypAntiBJetpT_OutOfSpace = store(new TH1D ( "HypAntiBJetpT_OutOfSpace", "Out Of Space Bkgd", 1200, 0, 1200 ));
        h_HypAntiBJetEta_OutOfSpace = store(new TH1D ( "HypAntiBJetEta_OutOfSpace", "Out Of Space Bkgd", 200, -5, 5 ));
        h_HypAntiBJetRapidity_OutOfSpace = store(new TH1D ( "HypAntiBJetRapidity_OutOfSpace", "Out Of Space Bkgd", 100, -5, 5 ));
        h_HypBJetpTLead_OutOfSpace = store(new TH1D ( "HypBJetpTLead_OutOfSpace", "Out Of Space Bkgd", 1200,0,1200 ));
        h_HypBJetEtaLead_OutOfSpace = store(new TH1D ( "HypBJetEtaLead_OutOfSpace", "Out Of Space Bkgd", 200,-5,5 ));
        h_HypBJetpTNLead_OutOfSpace = store(new TH1D ( "HypBJetpTNLead_OutOfSpace", "Out Of Space Bkgd", 1200,0,1200 ));
        h_HypBJetEtaNLead_OutOfSpace = store(new TH1D ( "HypBJetEtaNLead_OutOfSpace", "Out Of Space Bkgd", 200,-5,5 ));

        h_HypLLBarMass_OutOfSpace = store(new TH1D ( "HypLLBarMass_OutOfSpace", "Out Of Space Bkgd", 500, 0, 1000 ));
        h_HypLLBarpT_OutOfSpace = store(new TH1D ( "HypLLBarpT_OutOfSpace", "Out Of Space Bkgd", 200, 0, 1000 ));
        h_HypLLBarDPhi_OutOfSpace = store(new TH1D ( "HypLLBarDPhi_OutOfSpace", "Out Of Space Bkgd", 3200, 0., 3.2 ));
        h_HypLLBarDEta_OutOfSpace = store(new TH1D ( "HypLLBarDEta_OutOfSpace", "Out Of Space Bkgd", 400, -5, 5 ));
        h_HypLLBarCosPhiInTopRestFrame_OutOfSpace = store(new TH1D ( "HypLLBarCosPhiIntopRestFrame_OutOfSpace", "Out Of Space Bkgd", 1000, -1., 1. ));

        h_HypBBBarpT_OutOfSpace = store( new TH1D("HypBBBarpT_OutOfSpace", "Out Of Space Bkgd", 400, 0, 800));
        h_HypBBBarMass_OutOfSpace = store( new TH1D("HypBBBarMass_OutOfSpace", "Out Of Space Bkgd", 400, 0, 800));
        h_HypBBBarDPhi_OutOfSpace = store(new TH1D ( "HypBBBarDPhi_OutOfSpace", "Out Of Space Bkgd", 320, 0., 3.2 ));

        h_HypLeptonantiBjetMass_OutOfSpace = store(new TH1D ( "HypLeptonBjetMass_OutOfSpace", "Out Of Space Bkgd", 500, 0, 1000 ));
        h_HypAntiLeptonBjetMass_OutOfSpace = store(new TH1D ( "HypAntiLeptonBjetMass_OutOfSpace", "Out Of Space Bkgd", 500, 0, 1000 ));

        //Gen-pseudo matrices
        h_GenPseudoToppT = store(new TH2D ( "GenPseudoToppT", "Gen/Pseudo Matching", 1200,0,1200, 1200,0,1200 ));
        h_GenPseudoAntiToppT = store(new TH2D ( "GenPseudoAntiToppT", "Gen/Pseudo Matching", 1200,0,1200, 1200,0,1200 ));
        h_GenPseudoTopRapidity = store(new TH2D ( "GenPseudoTopRapidity", "Gen/Pseudo Matching", 400, -5, 5, 400, -5, 5 ));
        h_GenPseudoAntiTopRapidity = store(new TH2D ( "GenPseudoAntiTopRapidity", "Gen/Pseudo Matching", 400, -5, 5, 400, -5, 5 ));
        h_GenPseudoTTBarRapidity = store(new TH2D ( "GenPseudoTTBarRapidity", "Rapidity of TTbar System (HYP)", 400, -5, 5, 400, -5, 5 ));
        h_GenPseudoTTBarpT = store(new TH2D ( "GenPseudoTTBarpT", "pT of TTbar System (HYP)", 1200, 0, 1200, 1200, 0, 1200 ));
        h_GenPseudoTTBarMass = store(new TH2D ( "GenPseudoTTBarMass", "Mass of TTbar System (HYP)", 3000, 0, 3000, 3000, 0, 3000 ));

        h_GenPseudoLeptonpT = store(new TH2D ( "GenPseudoLeptonpT", "Gen/Pseudo Matching", 1200,0,1200, 1200,0,1200 ));
        h_GenPseudoAntiLeptonpT = store(new TH2D ( "GenPseudoAntiLeptonpT", "Gen/Pseudo Matching", 1200,0,1200, 1200,0,1200 ));
        h_GenPseudoLeptonEta = store(new TH2D ( "GenPseudoLeptonEta", "Gen/Pseudo Matching", 200,-5,5, 200,-5,5 ));
        h_GenPseudoAntiLeptonEta = store(new TH2D ( "GenPseudoAntiLeptonEta", "Gen/Pseudo Matching", 200,-5,5, 200,-5,5 ));

        h_GenPseudoBJetpT = store(new TH2D ( "GenPseudoBJetpT", "Gen/Pseudo Matching", 1200, 0, 1200, 1200, 0, 1200 ));
        h_GenPseudoAntiBJetpT = store(new TH2D ( "GenPseudoAntiBJetpT", "Gen/Pseudo Matching", 1200, 0, 1200, 1200, 0, 1200 ));
        h_GenPseudoBJetEta = store(new TH2D ( "GenPseudoBJetEta", "Gen/Pseudo Matching", 200, -5, 5, 200, -5, 5 ));
        h_GenPseudoAntiBJetEta = store(new TH2D ( "GenPseudoAntiBJetEta", "Gen/Pseudo Matching", 200, -5, 5, 200, -5, 5 ));

        h_GenPseudoJetMultpt30 = store(new TH2D("GenPseudoJetMultpt30", "Gen/Pseudo Matching",8,1.5,9.5,8,1.5,9.5));

        //Resolution
        h_GenPseudoToppT_reso = store(new TH1D ( "GenPseudoToppT_reso", "Gen/Pseudo Matching Reso", 2000, -2, 2 ));
        h_GenPseudoAntiToppT_reso = store(new TH1D ( "GenPseudoAntiToppT_reso", "Gen/Pseudo Matching Reso", 2000, -2, 2 ));
        h_GenPseudoTopRapidity_reso = store(new TH1D ( "GenPseudoTopRapidity_reso", "Gen/Pseudo Matching Reso", 2000, -2, 2 ));
        h_GenPseudoAntiTopRapidity_reso = store(new TH1D ( "GenPseudoAntiTopRapidity_reso", "Gen/Pseudo Matching Reso", 2000, -2, 2 ));
        h_GenPseudoTTBarRapidity_reso = store(new TH1D ( "GenPseudoTTBarRapidity_reso", "Rapidity of TTbar System (HYP) Reso", 2000, -2, 2 ));
        h_GenPseudoTTBarpT_reso = store(new TH1D ( "GenPseudoTTBarpT_reso", "pT of TTbar System (HYP) Reso", 2000, -2, 2 ));
        h_GenPseudoTTBarMass_reso = store(new TH1D ( "GenPseudoTTBarMass_reso", "Mass of TTbar System (HYP) Reso", 2000, -2, 2 ));

        h_GenPseudoLeptonpT_reso = store(new TH1D ( "GenPseudoLeptonpT_reso", "Gen/Pseudo Matching Reso", 1000, -1, 1 ));
        h_GenPseudoAntiLeptonpT_reso = store(new TH1D ( "GenPseudoAntiLeptonpT_reso", "Gen/Pseudo Matching Reso", 1000, -1, 1 ));
        h_GenPseudoLeptonEta_reso = store(new TH1D ( "GenPseudoLeptonEta_reso", "Gen/Pseudo Matching Reso", 1000, -1, 1 ));
        h_GenPseudoAntiLeptonEta_reso = store(new TH1D ( "GenPseudoAntiLeptonEta_reso", "Gen/Pseudo Matching Reso", 1000, -1, 1 ));

        h_GenPseudoBJetpT_reso = store(new TH1D ( "GenPseudoBJetpT_reso", "Gen/Pseudo Matching Reso", 1000, -1, 1 ));
        h_GenPseudoAntiBJetpT_reso = store(new TH1D ( "GenPseudoAntiBJetpT_reso", "Gen/Pseudo Matching Reso", 1000, -1, 1 ));
        h_GenPseudoBJetEta_reso = store(new TH1D ( "GenPseudoBJetEta_reso", "Gen/Pseudo Matching Reso", 1000, -1, 1 ));
        h_GenPseudoAntiBJetEta_reso = store(new TH1D ( "GenPseudoAntiBJetEta_reso", "Gen/Pseudo Matching Reso", 1000, -1, 1 ));

        h_GenPseudoJetMultpt30_reso = store(new TH1D("GenPseudoJetMultpt30_reso", "Gen/Pseudo Matching Reso",40,-0.2, 0.2 ));

    }
    //=======================PSEUDO-TOP: END==========================//

    //    double xbin[20]={30.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,320.,360.,400.};
    double xbin[21]={20.,30.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,320.,360.,400.};
    int nbin = 20;
    h_RecoJetMultQ0   = store(new TH1D("RecoJetMultQ0","Gap fraction Q0",nbin,xbin));
    h_HypJetMultQ0  = store(new TH1D("HypJetMultQ0",         "Gap fraction Q0 (HYP)", nbin,xbin));
    h_VisGenJetMultQ0 = store(new TH1D("VisGenJetMultQ0",         "Gap fraction Q0 (VisGEN)", nbin,xbin));
    h_GenRecoJetMultQ0 = store(new TH2D("GenRecoJetMultQ0","Gap fraction Q0", nbin,xbin,nbin,xbin));

    h_RecoJetExtra2Q0   = store(new TH1D("RecoJet2Q0","Gap fraction 2nd jet",nbin,xbin));
    h_HypJetExtra2Q0  = store(new TH1D("HypJet2Q0",         "Gap fraction 2nd jet (HYP)", nbin,xbin));
    h_VisGenJetExtra2Q0 = store(new TH1D("VisGenJet2Q0",         "Gap fraction 2nd jet (VisGEN)", nbin,xbin));
    h_GenRecoJetExtra2Q0 = store(new TH2D("GenRecoJet2Q0","Gap fraction 2nd jet", nbin,xbin,nbin,xbin));

    h_RecoJetMultQsum   = store(new TH1D("RecoJetMultQsum","Gap fraction Qsum",nbin,xbin));
    h_HypJetMultQsum  = store(new TH1D("HypJetMultQsum",         "Gap fraction Qsum (HYP)", nbin,xbin));
    h_VisGenJetMultQsum = store(new TH1D("VisGenJetMultQsum",         "Gap fraction Qsum (VisGEN)", nbin,xbin));
    h_GenRecoJetMultQsum = store(new TH2D("GenRecoJetMultQsum","Gap fraction Qsum", nbin,xbin,nbin,xbin));

    h_RecoJetMultTotal  = store(new TH1D("RecoJetMultTotal",         "Jet Multiplicity (HYP)", nbin,xbin));
    h_HypJetMultTotal  = store(new TH1D("HypJetMultTotal",         "Jet Multiplicity (HYP)", nbin,xbin));
    h_VisGenJetMultTotal = store(new TH1D("VisGenJetMultTotal",         "Gap fraction Qsum (VisGEN)", nbin,xbin));
    h_GenRecoJetMultTotal = store(new TH2D("GenRecoJetMultTotal","Gap fraction Qsum", nbin,xbin,nbin,xbin));

    h_ClosureTotalWeight = store(new TH1D("ClosureTotalWeight", "Total Weights from closure test",1,0,2));
    h_PDFTotalWeight = store(new TH1D("PDFTotalWeight", "PDF Weights",1,0,2));
    h_PSTotalWeight = store(new TH1D("PSTotalWeight", "PS Weights",1,0,2));

    h_ToppTWeight = store(new TH1D("ToppTWeight", "ToppT Reweights", 200, 0.8, 1.2));

    h_PUSF = store(new TH1D("PUSF", "PU SF per event", 200, 0.5, 1.5));
    h_TrigSF = store(new TH1D("TrigSF", "Trigger SF per event", 200, 0.5, 1.5));
    h_PrefiringWeight = store(new TH1D("PrefiringWeight", "Prefiring Weight per event", 200, 0.5, 1.5));
    h_LepSF = store(new TH1D("LepSF", "Lep. Id and Isol. SF per event", 200, 0.75, 1.25));
    h_BTagSF = store(new TH1D("BTagSF", "BTagging SF per event", 200 , 0.95, 1.15 ));
    h_BTagSF->Sumw2();
    h_KinRecoSF = store(new TH1D("KinRecoSF", "Kinematic Reco. SF per event", 200, 0.5, 1.5));
    h_weights = store(new TH1D("weights", "weights per event", 300, -200, 1200));
    h_EventWeight = store(new TH1D("EventWeight", "Event SF", 600, 0, 3));


    // Map for binned control plots
    binnedControlPlots_ = new std::map<std::string, std::pair<TH1*, std::vector<std::map<std::string, TH1*> > > >;

    CreateBinnedControlPlots(h_HypToppT, h_LeptonpT);
    CreateBinnedControlPlots(h_HypToppT, h_LeptonEta);
    CreateBinnedControlPlots(h_HypToppT, h_MET_preMETcut);
    CreateBinnedControlPlots(h_HypToppT, h_diLepMassFull);

    CreateBinnedControlPlots(h_HypTopRapidity, h_LeptonpT);
    CreateBinnedControlPlots(h_HypTopRapidity, h_LeptonEta);
    CreateBinnedControlPlots(h_HypTopRapidity, h_MET_preMETcut);
    CreateBinnedControlPlots(h_HypTopRapidity, h_diLepMassFull);

    // Histograms for resolution studies

    //h_RMSvsGenToppT = store(new TH2D ( "RMSvsGenToppT", "RMS vs Gen", 1000, 0, 1000, 2000, -1000, 1000 ));
    //h_RMSvsGenTopRapidity = store(new TH2D ( "RMSvsGenTopRapidity", "RMS vs Gen", 400, -5, 5, 400, -5, 5 ));
    //h_RMSvsGenToppTLead = store(new TH2D ( "RMSvsGenToppTLead", "RMS vs Gen", 1000, 0, 1000, 2000, -1000, 1000 ));
    //h_RMSvsGenTopRapidityLead = store(new TH2D ( "RMSvsGenTopRapidityLead", "RMS vs Gen", 400, -5, 5, 400, -5, 5 ));
    //h_RMSvsGenToppTNLead = store(new TH2D ( "RMSvsGenToppTNLead", "RMS vs Gen", 1000, 0, 1000, 2000, -1000, 1000 ));
    //h_RMSvsGenTopRapidityNLead = store(new TH2D ( "RMSvsGenTopRapidityNLead", "RMS vs Gen", 400, -5, 5, 400, -5, 5 ));
    //
    //h_RMSvsGenTTBarMass = store(new TH2D ( "RMSvsGenTTBarMass", "RMS vs Gen", 3000, 0, 3000, 6000, -3000, 3000 ));
    //h_RMSvsGenTTBarpT = store(new TH2D ( "RMSvsGenTTBarpT", "RMS vs Gen", 1000, 0, 1000, 2000, -1000, 1000 ));
    //h_RMSvsGenTTBarRapidity = store(new TH2D ( "RMSvsGenTTBarRapidity", "RMS vs Gen", 400, -5, 5, 400, -5, 5 ));
    //h_RMSvsGenTTBarDPhi = store(new TH2D ( "RMSvsGenTTBarDPhi", "RMS vs Gen", 350, 0, 3.5, 700, -3.5, 3.5 ));
    //h_RMSvsGenTTBarDeltaRapidity = store(new TH2D ( "RMSvsGenTTBarDeltaRapidity", "RMS vs Gen", 400, -5, 5, 400, -5, 5 ));

    //2d cs

    h_HypTTBarRapidityvsTTBarpT = store(new TH2D ("HypTTBarRapidityvsTTBarpT","TTBarRapidity vs TTBarpT;p_{T}^{t#bar{t}} [GeV];y(t#bar{t})",400,0,400,100,-2.4,2.4));
    h_VisGenTTBarRapidityvsTTBarpT = store(new TH2D ("VisGenTTBarRapidityvsTTBarpT","TTBarRapidity vs TTBarpT;p_{T}^{t#bar{t}} [GeV];y(t#bar{t})",400,0,400,100,-2.4,2.4));
    // ...

    // Book histograms of all analyzers
    this->bookAll();
}



void TopAnalysis::SlaveTerminate()
{

    this->clearAll();

    for (auto it = binnedControlPlots_->begin(); it != binnedControlPlots_->end(); ++it) {
        delete (*it).second.first;
    }
    delete binnedControlPlots_;

    // Defaults from AnalysisBase
    AnalysisBase::SlaveTerminate();
}



Bool_t TopAnalysis::Process ( Long64_t entry )
{ 
    // Defaults from AnalysisBase
    if(!AnalysisBase::Process(entry)) return kFALSE;


    // Use utilities without namespaces
    using ROOT::Math::VectorUtil::DeltaPhi;
    using ROOT::Math::VectorUtil::DeltaR;
    using ROOT::Math::VectorUtil::Angle;
    using namespace common;


    // Entry for object structs are not yet read, so reset
    this->resetObjectStructEntry();


    // Define the selection steps as strings
    std::string selectionStep("");

    //===CUT===
    // select events on generator level and access true level weights

    // Separate DY dilepton decays in lepton flavours
    if(this->failsDrellYanGeneratorSelection(this->zDecayModes(entry))) return kTRUE;

    // Separate dileptonic ttbar decays by channel at step0 and separate via tau decays from prompt decays
   
    if(this->failsTopGeneratorSelection(entry)) return kTRUE;

    // Count events for closure test here, where no more taus are available
    if (doClosureTest_) {
        static int closureTestEventCounter = 0;
        if (++closureTestEventCounter > closureMaxEvents_) return kTRUE;
    }

    // Determine all true level weights
    const double weightBRCorrection = this->brWDecayCorrection(entry);
    const double pdfWeight = this->weightPdf(entry);
    const double psWeight = this->weightPS(entry);
    const double weightAlphasPdf = this->weightAlphasPdf(entry);
    const double weightGenerator = this->weightGenerator(entry);
    const double weightMEScale  = this->weightMeScale(entry);
    const double weightMEFacScale  = this->weightMeFacScale(entry);
    const double weightMERenScale  = this->weightMeRenScale(entry);
    const double weightBFrag  = this->weightBFrag(entry);
    const double weightBSemilep  = this->weightBSemilep(entry);
    double weightTopPt = this->weightTopPtReweighting(entry);
    //    double weightTopPt = 1.;
    //    if (this->applyTopSlope2DGenLevelReweighting_) weightTopPt = this->weightTopSlope2DGenLevelReweighting(entry);
    //    else if (this->applyTopPtReweighting_) weightTopPt = this->weightTopPtReweighting(entry);
    const double weightPU = this->weightPileup(entry);
    const double trueLevelWeightNoPileupNoClosure = weightGenerator*weightBRCorrection*pdfWeight*psWeight*weightAlphasPdf*weightMEScale*weightMEFacScale*weightMERenScale*weightBFrag*weightBSemilep*weightTopPt;
    const double trueLevelWeightNoPileup = doClosureTest_ ? this->calculateClosureTestWeight(entry) : trueLevelWeightNoPileupNoClosure;
    //    const double trueLevelWeight = trueLevelWeightNoPileup*weightPU;
    // ajeeta
    const double weightPreScale = this->weightPreScale();
    const double trueLevelWeight = trueLevelWeightNoPileup*weightPU*weightPreScale;
    
    // for ttbar inclusive decays samples we need to include the weightBRCorrection in the renormalisationWeightSum, so check for it in a simple way and set the renormalisationWeights
    //bool isTTbarInclusiveSample(!outputfilename_.Contains("fromDilepton") && samplename_=="ttbarsignalplustau");
    //if(isTTbarInclusiveSample){
        this->renormalisationWeights(trueLevelWeight, weightGenerator*weightPU*weightBRCorrection, entry);
	//}else{
        //this->renormalisationWeights(trueLevelWeight, weightGenerator*weightPU, entry);
	//}

    const ttbar::GenLevelWeights genLevelWeights(weightBRCorrection, weightPU, weightGenerator,
						 pdfWeight, psWeight, weightAlphasPdf,
						 weightMEScale, weightMEFacScale, weightMERenScale,
						 weightBFrag, weightBSemilep, weightTopPt,
                                                 trueLevelWeightNoPileup, trueLevelWeight);

    h_PDFTotalWeight->Fill(1, pdfWeight);

    h_PSTotalWeight->Fill(1, psWeight);

    h_ToppTWeight->Fill(1, weightTopPt);

    const EventMetadata eventMetadataDummy;
    const CommonGenObjects commonGenObjectsDummy;

    const KinematicReconstructionSolutions kinematicReconstructionSolutionsDummy;
    const LooseKinRecoSolution looseKinematicReconstructionSolutionsDummy;

    // Access MC general generator info
    const CommonGenObjects& commonGenObjects = this->getCommonGenObjects(entry);

    // Access Top signal generator info
    const TopGenObjects& topGenObjects = this->getTopGenObjects(entry);

    // Access object selections from config
    const AnalysisConfig::Selections& selections = analysisConfig_.selections();

    // Generated jets
    const VLV& allGenJets =  topGenObjects.valuesSet_ ? *topGenObjects.allGenJets_ : VLV();
    std::vector<int> allGenJetIndices = initialiseIndices(allGenJets);
    std::vector<int> genJetIndices = allGenJetIndices;
    selectIndices(genJetIndices, allGenJets, LVeta, selections.genJetEtaCut_, false);
    selectIndices(genJetIndices, allGenJets, LVeta, -selections.genJetEtaCut_);
    selectIndices(genJetIndices, allGenJets, LVpt, selections.genJetPtCut_);

    if(selections.genDeltaRLeptonJetCut_ > 0.){
        // Vector of genLeptons from which genJets need to be separated in deltaR
        VLV allGenLeptons;
        VLV topBQuarksForJetCleaning_gen;
        if(topGenObjects.valuesSet_){
            if(topGenObjects.GenLepton_) allGenLeptons.push_back(*topGenObjects.GenLepton_);
            if(topGenObjects.GenAntiLepton_) allGenLeptons.push_back(*topGenObjects.GenAntiLepton_);
        }
        this->leptonCleanedJetIndices(genJetIndices, allGenJets, allGenLeptons, selections.genDeltaRLeptonJetCut_);
    }
    orderIndices(genJetIndices, allGenJets, LVpt);

    // Indices for Jet pT related studies
    std::vector<int> genJet40Indices = genJetIndices;
    selectIndices(genJet40Indices, allGenJets, LVpt, 40.0);
    std::vector<int> genJet60Indices = genJet40Indices;
    selectIndices(genJet60Indices, allGenJets, LVpt, 60.0);
    std::vector<int> genJet100Indices = genJet60Indices;
    selectIndices(genJet100Indices, allGenJets, LVpt, 100.0);

    // Match for all genJets all B hadrons
    std::vector<std::vector<int> > allGenJetBhadronIndices;
    std::vector<std::vector<int> > genJetBhadronIndices;
    std::vector<int> allGenBjetIndices;
    std::vector<int> genBjetIndices;
    if(topGenObjects.valuesSet_){
        allGenJetBhadronIndices = this->matchHadronsToGenJets(allGenJetIndices, allGenJets, *topGenObjects.genBHadJetIndex_);
        genJetBhadronIndices = this->matchHadronsToGenJets(genJetIndices, allGenJets, *topGenObjects.genBHadJetIndex_);
        allGenBjetIndices = this->genBjetIndices(allGenJetBhadronIndices);
        genBjetIndices = this->genBjetIndices(genJetBhadronIndices);
    }

    // Jet matchings for ttbar system
    int genBjetFromTopIndex(-1);
    int genAntiBjetFromTopIndex(-1);
    if(topGenObjects.valuesSet_){
        genBjetFromTopIndex = this->genBjetIndex(topGenObjects, 6);
        genAntiBjetFromTopIndex = this->genBjetIndex(topGenObjects, -6);
    }

    if (genBjetFromTopIndex == -2) genBjetFromTopIndex = -1;// FIXME: generatorTopEvent() and generatorTTbarjetsEvent() can't handle "-2"
    if (genAntiBjetFromTopIndex == -2) genAntiBjetFromTopIndex = -1;// -1 - no b-jet found, -2 - two b-jets from one t-quark in generator

    //const int numberOfGenBjets = genBjetIndices.size();
    //int leadingGenBjetIndex = numberOfGenBjets>0 ? genBjetIndices.at(0) : -1;
    //int nLeadingGenBjetIndex = numberOfGenBjets>1 ? genBjetIndices.at(1) : -1;

    int leadingGenBjetFromTopIndex = genBjetFromTopIndex;
    int nLeadingGenBjetFromTopIndex = genAntiBjetFromTopIndex;

    if(leadingGenBjetFromTopIndex <= -1 && nLeadingGenBjetFromTopIndex > -1) {// can be -1 or -2
        leadingGenBjetFromTopIndex = nLeadingGenBjetFromTopIndex;
        nLeadingGenBjetFromTopIndex = -1;
    } else if (leadingGenBjetFromTopIndex > -1 && nLeadingGenBjetFromTopIndex > -1) {// can be -1 or -2
        orderIndices(leadingGenBjetFromTopIndex, nLeadingGenBjetFromTopIndex, allGenJets, LVpt);
    }

    // Get indices of B and anti-B hadrons steming from ttbar system
    int BHadronIndex = genBjetFromTopIndex;//FIXME: switch everywhere to proper namings
    int AntiBHadronIndex = genAntiBjetFromTopIndex;//FIXME: here too

    // Access ttbar dilepton generator event
    LV LeadGenTop, NLeadGenTop;
    LV LeadGenLepton, NLeadGenLepton;
    LV LeadGenBJet, NLeadGenBJet;
    double genHT = -1;
    this->generatorTopEvent(LeadGenTop, NLeadGenTop,
                            LeadGenLepton, NLeadGenLepton,
                            LeadGenBJet, NLeadGenBJet,
                            genHT,
                            BHadronIndex, AntiBHadronIndex,
                            trueLevelWeightNoPileup, trueLevelWeight, topGenObjects);

    double jetHTGen = 0.;

    std::vector<int> genVisJetIndices;
    this->generatorVisJets(topGenObjects,genVisJetIndices);

    std::vector<int> genExtraJetIndices;
    this->generatorExtraJets(topGenObjects,genExtraJetIndices, BHadronIndex,AntiBHadronIndex);
    selectIndices(genExtraJetIndices, allGenJets, LVeta, selections.genJetEtaCut_, false);
    selectIndices(genExtraJetIndices, allGenJets, LVeta, -selections.genJetEtaCut_);
    selectIndices(genExtraJetIndices, allGenJets, LVpt, 40);

    std::vector<int> selectedExtraJetIndicesIso08_gen;

    if(selections.genDeltaRLeptonJetCut_ > 0.){
        VLV allGenLeptons;
        VLV topBQuarksForJetCleaning_gen;
        if(topGenObjects.valuesSet_){
            if(topGenObjects.GenLepton_) allGenLeptons.push_back(*topGenObjects.GenLepton_);
            if(topGenObjects.GenAntiLepton_) allGenLeptons.push_back(*topGenObjects.GenAntiLepton_);
	    	    if(BHadronIndex>-1 && topGenObjects.allGenJets_) topBQuarksForJetCleaning_gen.push_back((*topGenObjects.allGenJets_).at(BHadronIndex));
	    	    if(AntiBHadronIndex>-1 && topGenObjects.allGenJets_) topBQuarksForJetCleaning_gen.push_back((*topGenObjects.allGenJets_).at(AntiBHadronIndex));         
        }
        this->leptonCleanedJetIndices(genExtraJetIndices, allGenJets, allGenLeptons, selections.genDeltaRLeptonJetCut_);
	std::vector<int> genExtraJetIndices2 = genExtraJetIndices;

	this->leptonCleanedJetIndices(genExtraJetIndices2, allGenJets, topBQuarksForJetCleaning_gen, 0.8);
	selectedExtraJetIndicesIso08_gen = genExtraJetIndices2;
    }

    int n_gen_ExtraJets = genExtraJetIndices.size(); 

    orderIndices(genExtraJetIndices, allGenJets, LVpt);

    this->generatorTTbarjetsEvent(jetHTGen,
                                  BHadronIndex, AntiBHadronIndex,
                                  trueLevelWeight,
                                  genJetIndices, genJet40Indices, genJet60Indices,
                                  genJet100Indices, genExtraJetIndices,topGenObjects);
    




    //=======================PSEUDO-TOP: BEGIN==========================//

    // Access Top signal pseudo info
    const TopPseudoObjects& topPseudoObjects = this->getTopPseudoObjects(entry);
    LV LeadPseudoTop, NLeadPseudoTop;
    LV LeadPseudoLepton, NLeadPseudoLepton;
    LV LeadPseudoBJet, NLeadPseudoBJet;
    const VLV& allPseudoJets =  topPseudoObjects.valuesSet_ ? *topPseudoObjects.allPseudoJets_ : VLV();
    std::vector<int> allPseudoJetIndices, pseudoJetIndices, pseudoJet40Indices, pseudoJet60Indices, pseudoJet100Indices;

    // Retrieve information related to pseudo-top objects
    if (analysisConfig_.sampleComposition().pseudoTopMode_ > 0) {

        // Pseudo jets
        allPseudoJetIndices = initialiseIndices(allPseudoJets);
        pseudoJetIndices = allPseudoJetIndices;
        selectIndices(pseudoJetIndices, allPseudoJets, LVeta, selections.genJetEtaCut_, false);
        selectIndices(pseudoJetIndices, allPseudoJets, LVeta, -selections.genJetEtaCut_);
        selectIndices(pseudoJetIndices, allPseudoJets, LVpt, selections.genJetPtCut_);

        if(selections.genDeltaRLeptonJetCut_ > 0.){
            // Vector of pseudoLeptons from which pseudoJets need to be separated in deltaR
            VLV allPseudoLeptons;
            if(topPseudoObjects.valuesSet_){
                if(topPseudoObjects.PseudoLepton_) allPseudoLeptons.push_back(*topPseudoObjects.PseudoLepton_);
                if(topPseudoObjects.PseudoAntiLepton_) allPseudoLeptons.push_back(*topPseudoObjects.PseudoAntiLepton_);
            }
            this->leptonCleanedJetIndices(pseudoJetIndices, allPseudoJets, allPseudoLeptons, selections.genDeltaRLeptonJetCut_);
        }
        orderIndices(pseudoJetIndices, allPseudoJets, LVpt);

        // Indices for Jet pT related studies
        pseudoJet40Indices = pseudoJetIndices;
        selectIndices(pseudoJet40Indices, allPseudoJets, LVpt, 40.0);
        pseudoJet60Indices = pseudoJet40Indices;
        selectIndices(pseudoJet60Indices, allPseudoJets, LVpt, 60.0);
        pseudoJet100Indices = pseudoJet60Indices;
        selectIndices(pseudoJet100Indices, allPseudoJets, LVpt, 100.0);

        // Access ttbar dilepton pseudotop event
        this->pseudoTopEvent(LeadPseudoTop, NLeadPseudoTop,
                                LeadPseudoLepton, NLeadPseudoLepton,
                                LeadPseudoBJet, NLeadPseudoBJet,
                                trueLevelWeightNoPileup, trueLevelWeight,
                                pseudoJetIndices, pseudoJet40Indices, pseudoJet60Indices, pseudoJet100Indices,
                                topPseudoObjects);
    }
    //=======================PSEUDO-TOP: END==========================//

    selectionStep = "0";
    //    std::vector<int> gen_extrajetCounts = {genExtraJetIndices.size(), selectedExtraJetIndicesIso04_gen.size(), selectedExtraJetIndicesIso08_gen.size()};
    

    const ttbar::GenObjectIndices genObjectIndices(BHadronIndex, AntiBHadronIndex, -1, -1, -1, -1, -1, -1, 
						   selectedExtraJetIndicesIso08_gen,
						   genVisJetIndices);    
    
    const ttbar::RecoObjectIndices recoObjectIndicesDummy({0},{0},{0},0,0,0,0,0,0,
							  {0},{0},{0});

    ttbar::RecoLevelWeights recoLevelWeightsDummy(0,0,0,0,0,0);

    const RecoObjects recoObjectsDummy;

    // Access event meta data
    const EventMetadata eventMetadata = this->getEventMetadata(entry);
   
    this->fillAll(selectionStep,
                  eventMetadata,
                  recoObjectsDummy, commonGenObjects,
                  topGenObjects, topPseudoObjects,
                  kinematicReconstructionSolutionsDummy,
                  looseKinematicReconstructionSolutionsDummy,
                  genObjectIndices, recoObjectIndicesDummy,
                  genLevelWeights, recoLevelWeightsDummy,
                  1.,1.);
    


      const RecoObjects& recoObjects = this->getRecoObjects(entry);

      selectionStep = "lhAllGenLevel";
      this->fillAll(selectionStep,
                    eventMetadata,
                    recoObjects, commonGenObjects,
                    topGenObjects, topPseudoObjects,
                    kinematicReconstructionSolutionsDummy,
                    looseKinematicReconstructionSolutionsDummy,
                    genObjectIndices, recoObjectIndicesDummy,
                    genLevelWeights, recoLevelWeightsDummy,
                    1.,1.);

    //===CUT===
    //selectionStep = "1a";
    // Check if event should be excluded due to MET filters
    if(this->failsMetFilters(entry)) return kTRUE;

    //===CUT===
    //selectionStep = "1b";
    // Check if event was triggered
    if(this->failsDileptonTrigger(entry)) return kTRUE;

    //===CUT===
    selectionStep = "1";
    // Check if the first primary vertex is of good quality
    if(!this->firstVertexIsGood(entry)) return kTRUE;

    // === FULL RECO OBJECT SELECTION === (can thus be used at each selection step)

    // Access objects info
    // const RecoObjects& recoObjects = this->getRecoObjects(entry);

    // Get allLepton indices, apply selection cuts and order them by pt (beginning with the highest value)
    const VLV& allLeptons = *recoObjects.allLeptons_;
    const std::vector<int>& lepPdgId = *recoObjects.lepPdgId_;
    const std::vector<float>& lepSCEta = *recoObjects.lepSCEta_; // FIXME: works since v037
    std::vector<int> allLeptonIndices = initialiseIndices(allLeptons);
    if(selections.muonIsoCut_ > -999.){
        std::vector<int> allElectronIndices = allLeptonIndices;
        std::vector<int> onlyMuonIndices = allLeptonIndices;
        std::vector<int> onlyAntiMuonIndices = allLeptonIndices;
        selectIndices(allElectronIndices, lepPdgId, 11, false);
        selectIndices(allElectronIndices, lepPdgId, -11);
        selectIndices(onlyMuonIndices, lepPdgId, 13);
        selectIndices(onlyAntiMuonIndices, lepPdgId, -13, false);
        std::vector<int> allMuonIndices = mergeIndices(onlyMuonIndices, onlyAntiMuonIndices);
        selectIndices(allMuonIndices, *recoObjects.lepIso_, selections.muonIsoCut_, false);
        allLeptonIndices.clear();
        allLeptonIndices = mergeIndices(allElectronIndices, allMuonIndices);
    }
    selectMuonIndicesID    (allLeptonIndices, recoObjects.lepID_MuonTight_   , 1, &lepPdgId);
    selectElectronIndicesID(allLeptonIndices, recoObjects.lepID_ElecCutBased_, 4, &lepPdgId);
    selectIndices(allLeptonIndices, allLeptons, LVeta, selections.leptonEtaCut_, false);
    selectIndices(allLeptonIndices, allLeptons, LVeta, -selections.leptonEtaCut_);
    selectIndices(allLeptonIndices, allLeptons, LVpt, selections.leptonPtCut_);
    orderIndices(allLeptonIndices, allLeptons, LVpt);

    // Get indices of leptons and antiLeptons separated by charge, and get the leading ones if they exist
    std::vector<int> leptonIndices = allLeptonIndices;
    std::vector<int> antiLeptonIndices = allLeptonIndices;
    selectIndices(leptonIndices, lepPdgId, 0);
    selectIndices(antiLeptonIndices, lepPdgId, 0, false);

    // For 8 TeV lepton pair definition
    const int numberOfLeptons = leptonIndices.size();
    const int numberOfAntiLeptons = antiLeptonIndices.size();
    const int leptonIndex = numberOfLeptons>0 ? leptonIndices.at(0) : -1;
    const int antiLeptonIndex = numberOfAntiLeptons>0 ? antiLeptonIndices.at(0) : -1;

    // For 13 TeV lepton pair definition
    const int LeptonsNumber = allLeptonIndices.size();
    const int firstLeptonIndex = LeptonsNumber>0 ? allLeptonIndices.at(0) : -1;
    const int secondLeptonIndex = LeptonsNumber>1 ? allLeptonIndices.at(1) : -1;

    // In case of an existing opposite-charge dilepton system, get indices for leading and next-to-leading lepton
    int leadingLeptonIndex(-1);
    int nLeadingLeptonIndex(-1);
    if(analysisConfig_.general().era_ == Era::run1_8tev && (numberOfLeptons>0 && numberOfAntiLeptons>0)){
        leadingLeptonIndex = leptonIndex;
        nLeadingLeptonIndex = antiLeptonIndex;
        orderIndices(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, LVpt);
    } else if((analysisConfig_.general().era_ == Era::run2_13tev_50ns || analysisConfig_.general().era_ == Era::run2_13tev_25ns
        || analysisConfig_.general().era_ == Era::run2_13tev_2015_25ns || analysisConfig_.general().era_ == Era::run2_13tev_2016_25ns
        || analysisConfig_.general().era_ == Era::run2_13tev_2016 || analysisConfig_.general().era_ == Era::run2_13tev_2017 || analysisConfig_.general().era_ == Era::run2_13tev_2018
        || analysisConfig_.general().era_ == Era::run2_UL_13tev_2016preVFP || analysisConfig_.general().era_ == Era::run2_UL_13tev_2016postVFP || analysisConfig_.general().era_ == Era::run2_UL_13tev_2017 || analysisConfig_.general().era_ == Era::run2_UL_13tev_2018)
        && LeptonsNumber>1){
        leadingLeptonIndex = firstLeptonIndex;
        nLeadingLeptonIndex = secondLeptonIndex;
        orderIndices(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, LVpt);
    }

    const bool hasLeptonPair = this->hasLeptonPair(leadingLeptonIndex, nLeadingLeptonIndex, lepPdgId);
    bool passLeadingLeptonPtCut = false;
    if(leadingLeptonIndex > -1) passLeadingLeptonPtCut = allLeptons.at(leadingLeptonIndex).pt() > analysisConfig_.selections().leadingLeptonPtCut_;
    const bool hasMoreThan2Leptons = LeptonsNumber > 2; // Veto based on this quantity needed in order to match fiducial phase space definition for differential xsec

    // Get two indices of the two leptons in the right order for trigger scale factor, if existing
    int leptonXIndex(leadingLeptonIndex);
    int leptonYIndex(nLeadingLeptonIndex);
    if(hasLeptonPair){
        //in ee and mumu channel leptonX must be the highest pt lepton, i.e. this is already correct
        // in emu channel leptonX must be electron
        if(std::abs(lepPdgId.at(leptonXIndex)) != std::abs(lepPdgId.at(leptonYIndex))){
            orderIndices(leptonYIndex, leptonXIndex, lepPdgId, true);
        }
    }

    // Get dilepton system, if existing
    const LV dummyLV(0.,0.,0.,0.);
    const LV dilepton(hasLeptonPair ? allLeptons.at(leadingLeptonIndex)+allLeptons.at(nLeadingLeptonIndex) : dummyLV);

    // Get jet indices, apply selection cuts and order them by pt (beginning with the highest value)
    const VLV& jets = *recoObjects.jets_;
    std::vector<int> jetIndices = initialiseIndices(jets);
    
    if((analysisConfig_.general().era_ == Era::run2_13tev_50ns || analysisConfig_.general().era_ == Era::run2_13tev_25ns
      || analysisConfig_.general().era_ == Era::run2_13tev_2015_25ns || analysisConfig_.general().era_ == Era::run2_13tev_2017
      || analysisConfig_.general().era_ == Era::run2_13tev_2016_25ns || analysisConfig_.general().era_ == Era::run2_13tev_2016
      || analysisConfig_.general().era_ == Era::run2_13tev_2018
	  || analysisConfig_.general().era_ == Era::run2_UL_13tev_2017 || analysisConfig_.general().era_ == Era::run2_UL_13tev_2016preVFP || analysisConfig_.general().era_ == Era::run2_UL_13tev_2016postVFP
      || analysisConfig_.general().era_ == Era::run2_UL_13tev_2018) && selections.deltaRLeptonJetCut_ > 0.){
        // Vector of leptons from which jets need to be separated in deltaR
        VLV leptonsForJetCleaning;
        for(const int index : allLeptonIndices) leptonsForJetCleaning.push_back(allLeptons.at(index));
        this->leptonCleanedJetIndices(jetIndices, jets, leptonsForJetCleaning, selections.deltaRLeptonJetCut_);
    }
    std::vector<int> jetIndicesForExtraJetStudies = jetIndices; // remains the same as initial collection at 8 TeV, but includes DeltaR(l,j) cleaning at 13 TeV: before application of other cuts

    std::vector<int> ExtraJetIndicesIso08;
    
    selectIndices(jetIndices, jets, LVeta, selections.jetEtaCut_, false);
    selectIndices(jetIndices, jets, LVeta, -selections.jetEtaCut_);
    selectIndices(jetIndices, jets, LVpt, selections.jetPtCut_);
    if(recoObjects.jetPFID_){

      selectIndices(jetIndices, *recoObjects.jetPFID_, selections.jetPFIDWorkingPoint_);
    }
    else {
      throw std::runtime_error("TopAnalysis::Process -- null pointer in recoObjects.jetPFID_");
    }
    const JetPileupId::WorkingPoint jetPileupIdWorkingPoint(selections.jetPileupIdWorkingPoint_);
    if(jetPileupIdWorkingPoint != JetPileupId::none){
        const std::vector<int>& jetPileupIdFlags = *recoObjects.jetPileupIdFlag_;
        std::vector<bool> jetFlagsPassPUID_ = utils::flagsForPassingCutInPtRange(jetPileupIdFlags, JetPileupId::cutValue(analysisConfig_.selections().jetPileupIdWorkingPoint_), jets, selections.jetPileupIdMinPt_, selections.jetPileupIdMaxPt_);
        utils::selectIndicesFromFlags(jetIndices, jetFlagsPassPUID_);

    }

    // apply jet veto maps

    std::vector<bool> jetFlagsPassVetoMap;

    const float jetVetoMapCutValue(analysisConfig_.selections().jetVetoMapCutValue_);
    if(jetVetoMapCutValue < 999.){
      jetFlagsPassVetoMap = utils::flagsForPassingJetVetoMap(jetVetoMap_, jets, analysisConfig_.selections().jetVetoMapCutValue_);
      utils::selectIndicesFromFlags(jetIndices, jetFlagsPassVetoMap);
    }


    orderIndices(jetIndices, jets, LVpt);
    const int numberOfJets = jetIndices.size();
    const bool has2Jets = numberOfJets > 1;

    //Anya: 2D iso
    if(btopTag) this->lepton2dIsolationIndices(allLeptonIndices, allLeptons, jets);

    // Get b-jet indices, apply selection cuts
    // and apply b-tag efficiency MC correction using random number based tag flipping (if requested correction mode is applied)
    // and order b-jets by btag discriminator (beginning with the highest value)
    const std::vector<float>& jetBtags = *recoObjects.jetBtags_;
    const std::vector<int>& jetFlavour = *commonGenObjects.jetFlavour_;
    std::vector<int> bjetIndices = jetIndices;
    selectIndices(bjetIndices, jetBtags, this->btagCutValue());
    this->retagJets(bjetIndices, jetIndices, jets, jetFlavour, jetBtags);
    orderIndices(bjetIndices, jetBtags);
    const int numberOfBjets = bjetIndices.size();
    const bool hasBtag = numberOfBjets > 0;

    // Get MET
    const LV& met = *recoObjects.met_;
    const bool hasMetOrEmu = this->channel()==Channel::emu || met.Pt() > selections.metCut_;
    

    const ttbar::RecoObjectIndices recoObjectIndices(allLeptonIndices,
                                                     leptonIndices, antiLeptonIndices,
                                                     leptonIndex, antiLeptonIndex,
                                                     leadingLeptonIndex, nLeadingLeptonIndex,
                                                     leptonXIndex, leptonYIndex,
						     ExtraJetIndicesIso08, 
						     jetIndices, bjetIndices);
    

    const std::vector<int>& genJetMatchedRecoBjetIndices = this->matchRecoToGenJets(jetIndices, jets, allGenBjetIndices, allGenJets);

    int matchedBjetFromTopIndex = BHadronIndex>=0 ? genJetMatchedRecoBjetIndices.at(BHadronIndex) : -1;
    int matchedAntiBjetFromTopIndex = AntiBHadronIndex>=0 ? genJetMatchedRecoBjetIndices.at(AntiBHadronIndex) : -1;

    const ttbar::GenObjectIndices genObjectIndicesWithMatching(BHadronIndex, AntiBHadronIndex, matchedBjetFromTopIndex, matchedAntiBjetFromTopIndex, -1, -1, -1, -1,
							       selectedExtraJetIndicesIso08_gen, 
							       genVisJetIndices);    


    // Determine all reco level weights
    const double weightLeptonSF = this->weightLeptonSF(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId, lepSCEta, this->systematic().name()); // FIXME: works since v037
    double weightTriggerSF;
    if(analysisConfig_.corrections().triggerSFReadoutMode_ == 1){
      weightTriggerSF = this->weightTriggerSF(leptonXIndex, leptonYIndex, allLeptons);
    }else if(analysisConfig_.corrections().triggerSFReadoutMode_== 2){
      weightTriggerSF = this->weightTriggerSF(leptonXIndex, leptonYIndex, allLeptons, lepPdgId);
    }else if(analysisConfig_.corrections().triggerSFReadoutMode_== 0){
    	weightTriggerSF = 1.0;
    }
    double l1PrefiringWeight = this->weightL1Prefiring(recoObjects);
    h_PrefiringWeight->Fill(l1PrefiringWeight, 1.);

    const double weightNoPileup = trueLevelWeightNoPileup*weightTriggerSF*weightLeptonSF*l1PrefiringWeight;
    const double weightBtagSF = this->weightBtagSF(jetIndices, jets, jetFlavour, jetBtags);
    const double weightKinReco = this->weightKinReco();
    const double weightLooseKinReco = this->weightLooseKinReco();

    // The weight to be used for filling the histograms
    //    double weight = weightNoPileup*weightPU;
    //    double weight = weightNoPileup*weightPU;

    double weight = weightNoPileup*weightPU*weightPreScale; //ajeeta

    ttbar::RecoLevelWeights recoLevelWeights(weightLeptonSF, weightTriggerSF, weightBtagSF,
                                             weightNoPileup, weightPU, weight);
   

    this->fillAll(selectionStep,
                  eventMetadata,
                  recoObjects, commonGenObjects,
                  topGenObjects, topPseudoObjects,
                  kinematicReconstructionSolutionsDummy,
                  looseKinematicReconstructionSolutionsDummy,
                  genObjectIndicesWithMatching, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight/weightLeptonSF, weight/weightLeptonSF);
   

    h_PUSF->Fill(weightPU, 1);


    //===CUT===
    selectionStep = "2";
    // we need an OS lepton pair
    if (!hasLeptonPair || !passLeadingLeptonPtCut) return kTRUE;
    if (selections.vetoEventsWithMoreThan2Leptons_ == 1  && hasMoreThan2Leptons) return kTRUE;
    this->fillAll(selectionStep,
                  eventMetadata,
                  recoObjects, commonGenObjects,
                  topGenObjects, topPseudoObjects,
                  kinematicReconstructionSolutionsDummy,
                  looseKinematicReconstructionSolutionsDummy,
                  genObjectIndicesWithMatching, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight, weight);

    //===CUT===
    selectionStep = "3";
    // with at least 20 GeV invariant mass
    if (dilepton.M() < 20.) return kTRUE;

    h_TrigSF->Fill(weightTriggerSF, 1.);
    h_LepSF->Fill(weightLeptonSF, 1.);

    h_vertMulti_step3->Fill(recoObjects.vertMulti_, weight);
    h_vertMulti_noPU_step3->Fill(recoObjects.vertMulti_, weightNoPileup);

    h_diLepMassFull->Fill(dilepton.M(), weight);

    if(topGenObjects.valuesSet_ && hasLeptonPair && has2Jets)
    {// Set of histograms needed to estimate the efficiency and acceptance requested by the TopXSection conveners
        h_GenAll_RecoCuts_noweight->Fill((*topGenObjects.GenTop_).M(), trueLevelWeightNoPileup);
        h_GenAll_RecoCuts->Fill((*topGenObjects.GenTop_).M(), weight);
    }

    // Access kinematic reconstruction info
    const KinematicReconstructionSolutions kinematicReconstructionSolutions = !this->makeBtagEfficiencies() ? this->kinematicReconstructionSolutions(leptonIndex, antiLeptonIndex, jetIndices, bjetIndices, allLeptons, jets, jetBtags, met) : kinematicReconstructionSolutionsDummy;
    const bool hasSolution = kinematicReconstructionSolutions.numberOfSolutions();

    // Access loose kinematic reconstruction info
    const LooseKinRecoSolution looseKinematicReconstructionSolution = !this->makeBtagEfficiencies() ? this->looseKinRecoSolutions(leptonIndex, antiLeptonIndex, jetIndices, bjetIndices, allLeptons, jets, met) : looseKinematicReconstructionSolutionsDummy;
    const bool hasLooseSolution = looseKinematicReconstructionSolution.hasSolution();
    double weightLooseKinRecoTot = weight*weightLooseKinReco;

    this->fillAll(selectionStep,
                  eventMetadata,
                  recoObjects, commonGenObjects,
                  topGenObjects, topPseudoObjects,
                  kinematicReconstructionSolutions,
                  looseKinematicReconstructionSolutionsDummy,
                  genObjectIndicesWithMatching, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight,weight);



    //****************************************
    // Handle inverted Z cut
    const bool isZregion = dilepton.M() > 76. && dilepton.M() < 106.;
    if ( isZregion ) {
        double fullWeights = weight;
        selectionStep = "4zWindow";
        this->fillAll(selectionStep,
                  eventMetadata,
                  recoObjects, commonGenObjects,
		  topGenObjects, topPseudoObjects,
                  kinematicReconstructionSolutions,
                  looseKinematicReconstructionSolutionsDummy,
                  genObjectIndicesWithMatching, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  fullWeights,fullWeights);
	h_jetMulti_Zwindow->Fill(numberOfJets, weight);

        if ( has2Jets ) {
            selectionStep = "5zWindow";
            this->fillAll(selectionStep,
                  eventMetadata,
                  recoObjects, commonGenObjects,
		  topGenObjects, topPseudoObjects,
                  kinematicReconstructionSolutions,
                  looseKinematicReconstructionSolutionsDummy,
                  genObjectIndicesWithMatching, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  fullWeights,fullWeights);
	    h_MET_Zwindow->Fill(met.Pt(), weight);

            if ( hasMetOrEmu ) {
                selectionStep = "6zWindow";
                this->fillAll(selectionStep,
                  eventMetadata,
                  recoObjects, commonGenObjects,
		  topGenObjects, topPseudoObjects,
                  kinematicReconstructionSolutions,
                  looseKinematicReconstructionSolutionsDummy,
                  genObjectIndicesWithMatching, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  fullWeights,fullWeights);
		h_BjetMulti_Zwindow->Fill(numberOfBjets, weight);

                if ( hasBtag ) {
                    selectionStep = "7zWindow";
                    fullWeights *= weightBtagSF;
                    this->fillAll(selectionStep,
                        eventMetadata,
                        recoObjects, commonGenObjects,
			topGenObjects, topPseudoObjects,
                        kinematicReconstructionSolutions,
                        looseKinematicReconstructionSolutionsDummy,
                        genObjectIndicesWithMatching, recoObjectIndices,
                        genLevelWeights, recoLevelWeights,
                        fullWeights,fullWeights);

                    if ( hasSolution ) {
                        fullWeights *= weightKinReco;
                        selectionStep = "8zWindow";
                        this->fillAll(selectionStep,
                            eventMetadata,
                            recoObjects, commonGenObjects,
			    topGenObjects, topPseudoObjects,
                            kinematicReconstructionSolutions,
                            looseKinematicReconstructionSolutionsDummy,
                            genObjectIndicesWithMatching, recoObjectIndices,
                            genLevelWeights, recoLevelWeights,
                            fullWeights,fullWeights);
                    }
                }
            }
        }
    }

    //=== CUT ===
    selectionStep = "4";
    //Exclude the Z window
    if (this->channel() != Channel::emu && isZregion) return kTRUE;

    this->fillAll(selectionStep,
                  eventMetadata,
                  recoObjects, commonGenObjects,
                  topGenObjects, topPseudoObjects,
                  kinematicReconstructionSolutions,
                  looseKinematicReconstructionSolutionsDummy,
                  genObjectIndicesWithMatching, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight,weight);

    h_jetMulti_diLep->Fill(numberOfJets, weight);

    h_LeptonpT_diLep->Fill((*recoObjects.allLeptons_).at(leptonIndex).Pt(), weight);
    h_AntiLeptonpT_diLep->Fill((*recoObjects.allLeptons_).at(antiLeptonIndex).Pt(), weight);
    h_LeptonEta_diLep->Fill((*recoObjects.allLeptons_).at(leptonIndex).Eta(), weight);
    h_AntiLeptonEta_diLep->Fill((*recoObjects.allLeptons_).at(antiLeptonIndex).Eta(), weight);

    h_MET_diLep->Fill(met.Pt(), weight);
    //loop over both leptons
    for (const int index : {leadingLeptonIndex, nLeadingLeptonIndex}) {
        if ( std::abs((*recoObjects.lepPdgId_).at(index)) == 11 ) {
            h_ElectronpT_diLep->Fill((*recoObjects.allLeptons_).at(index).Pt(), weight);
            h_ElectronEta_diLep->Fill((*recoObjects.allLeptons_).at(index).Eta(), weight);
        }
        else if ( std::abs((*recoObjects.lepPdgId_).at(index)) == 13 ) {
            h_MuonpT_diLep->Fill((*recoObjects.allLeptons_).at(index).Pt(), weight);
            h_MuonEta_diLep->Fill((*recoObjects.allLeptons_).at(index).Eta(), weight);
        }
    }

    //=== CUT ===
    selectionStep = "5";
    //Require at least two jets > 30 GeV (check for > 30 needed because we might have 20 GeV jets in our NTuple)
    if(!has2Jets) return kTRUE;

    this->fillAll(selectionStep,
                  eventMetadata,
                  recoObjects, commonGenObjects,
                  topGenObjects, topPseudoObjects,
                  kinematicReconstructionSolutions,
                  looseKinematicReconstructionSolutionsDummy,
                  genObjectIndicesWithMatching, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight,weight);

    h_MET_preMETcut->Fill(met.Pt(), weight);

    //=== CUT ===
    selectionStep = "6";
    //Require MET > 30 GeV in non-emu channels
    if (!hasMetOrEmu) return kTRUE;

    this->fillAll(selectionStep,
                  eventMetadata,
                  recoObjects, commonGenObjects,
                  topGenObjects, topPseudoObjects,
                  kinematicReconstructionSolutions,
                  looseKinematicReconstructionSolutionsDummy,
                  genObjectIndicesWithMatching, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight,weight);

    h_LeptonpT_postMETcut->Fill((*recoObjects.allLeptons_).at(leptonIndex).Pt(), weight);
    h_AntiLeptonpT_postMETcut->Fill((*recoObjects.allLeptons_).at(antiLeptonIndex).Pt(), weight);
    h_LeptonEta_postMETcut->Fill((*recoObjects.allLeptons_).at(leptonIndex).Eta(), weight);
    h_AntiLeptonEta_postMETcut->Fill((*recoObjects.allLeptons_).at(antiLeptonIndex).Eta(), weight);

    //loop over both leptons
    for(const int index : {leadingLeptonIndex, nLeadingLeptonIndex}){
        if ( std::abs((*recoObjects.lepPdgId_).at(index)) == 11 ) {
            h_ElectronpT_postMETcut->Fill((*recoObjects.allLeptons_).at(index).Pt(), weight);
            h_ElectronEta_postMETcut->Fill((*recoObjects.allLeptons_).at(index).Eta(), weight);
        }
        else if ( std::abs((*recoObjects.lepPdgId_).at(index)) == 13 ) {
            h_MuonpT_postMETcut->Fill((*recoObjects.allLeptons_).at(index).Pt(), weight);
            h_MuonEta_postMETcut->Fill((*recoObjects.allLeptons_).at(index).Eta(), weight);
        }
    }

    h_jetMulti_noBTag->Fill(numberOfJets, weight);
    h_BjetMulti_noBTag->Fill(numberOfBjets, weight);

    // Fill b-tagging efficiencies if required for given correction mode, and in case do not process further steps
    this->fillBtagEfficiencyHistos(jetIndices, jetBtags, jets, jetFlavour, weight);
    if(this->makeBtagEfficiencies()) return kTRUE;

    //=== CUT ===
    selectionStep = "7";
    //Require at least one b tagged jet
    if (!hasBtag) return kTRUE;

    weight *= weightBtagSF;

    // // Access loose kinematic reconstruction info
    // obsolete? moved above
    // const LooseKinRecoSolution looseKinematicReconstructionSolution = !this->makeBtagEfficiencies() ? this->looseKinRecoSolutions(leptonIndex, antiLeptonIndex, jetIndices, bjetIndices, allLeptons, jets, met) : looseKinematicReconstructionSolutionsDummy;

    this->fillAll(selectionStep,
                  eventMetadata,
                  recoObjects, commonGenObjects,
                  topGenObjects, topPseudoObjects,
                  kinematicReconstructionSolutions,
                  looseKinematicReconstructionSolutionsDummy,
                  genObjectIndicesWithMatching, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight,weight);

    h_BTagSF->Fill(weightBtagSF);

    if (RUNSYNC) {
        static int fullSelectionCounter = 0;
        if (fullSelectionCounter == 0)
        {
            std::cout << "Selected#\tRun\tEvent\tlep+\tlep-\tMll\tNJets\tjet0\tjet1\tNTags\tGenJet1\tGenJet2\tMet\tGenMet\tt/tbar_decay\n"
            << std::setprecision(2) << std::fixed;
        std::cout << "Event#" << ++fullSelectionCounter << ":\t" << eventMetadata.runNumber_ << "\t" << eventMetadata.eventNumber_ << "\t" << (*recoObjects.allLeptons_).at(antiLeptonIndex) << "\t"
            << (*recoObjects.allLeptons_).at(leptonIndex) << "\t"
            << dilepton.M() << "\t" << numberOfJets << "\t"
            << (*recoObjects.jets_).at(jetIndices.at(0)) << "\t" << (*recoObjects.jets_).at(jetIndices.at(1)) << "\t" << numberOfBjets << "\t"
            << (*commonGenObjects.associatedGenJet_).at(jetIndices.at(0)) << "\t" << (*commonGenObjects.associatedGenJet_).at(jetIndices.at(1)) << "\t"
            << met.Pt() << "\t" << (*topGenObjects.GenMet_).Pt() << "\t"
            << topDecayModeString()
            << "\n";
        }
    }

    h_BjetMulti->Fill(numberOfBjets, weight);
    h_jetMulti->Fill(numberOfJets, weight);

    double jetHT = getJetHT(jetIndices, (*recoObjects.jets_));
    h_jetHT->Fill(jetHT, weight);

    for ( size_t i = 0; i < 2; ++i ) {
        const int index = jetIndices.at(i);
        h_jetpT->Fill((*recoObjects.jets_).at(index).Pt(), weight);
    }

    h_LeptonPt_BeforeKinReco->Fill((*recoObjects.allLeptons_).at(leptonIndex).Pt(), weight);
    h_LeptonPt_BeforeKinReco->Fill((*recoObjects.allLeptons_).at(antiLeptonIndex).Pt(), weight);
    h_LeptonEta_BeforeKinReco->Fill((*recoObjects.allLeptons_).at(leptonIndex).Eta(), weight);
    h_LeptonEta_BeforeKinReco->Fill((*recoObjects.allLeptons_).at(antiLeptonIndex).Eta(), weight);

    //loop over both leptons
    for(const int index : {leadingLeptonIndex, nLeadingLeptonIndex}){
        if ( std::abs((*recoObjects.lepPdgId_).at(index)) == 11 ) {
            h_ElectronpT_BeforeKinReco->Fill((*recoObjects.allLeptons_).at(index).Pt(), weight);
            h_ElectronEta_BeforeKinReco->Fill((*recoObjects.allLeptons_).at(index).Eta(), weight);
        }
        else if ( std::abs((*recoObjects.lepPdgId_).at(index)) == 13 ) {
            h_MuonpT_BeforeKinReco->Fill((*recoObjects.allLeptons_).at(index).Pt(), weight);
            h_MuonEta_BeforeKinReco->Fill((*recoObjects.allLeptons_).at(index).Eta(), weight);
        }
    }

    h_MET_step7->Fill(met.Pt(), weight);

    for (const int index : bjetIndices) {
        h_BjetEta_BeforeKinReco->Fill((*recoObjects.jets_).at(index).Eta(), weight);
        h_BjetpT_BeforeKinReco->Fill((*recoObjects.jets_).at(index).Pt(), weight);
    }

    // obsolete? moved above
    // const bool hasLooseSolution = looseKinematicReconstructionSolution.hasSolution();
    // double weightLooseKinRecoTot = weight*weightLooseKinReco;
    if(hasLooseSolution){
      selectionStep = "7L";
      this->fillAll(selectionStep,
                    eventMetadata,
                    recoObjects, commonGenObjects,
                    topGenObjects, topPseudoObjects,
                    kinematicReconstructionSolutions,
                    looseKinematicReconstructionSolution,
                    genObjectIndicesWithMatching, recoObjectIndices,
                    genLevelWeights, recoLevelWeights,
                    weight,weightLooseKinRecoTot);

        h_jetMultiXSec_Loose->Fill(numberOfJets, weightLooseKinRecoTot);
        LV hypttbar_Loose(looseKinematicReconstructionSolution.TTbar());
        h_HypTTBarMass_Loose->Fill(hypttbar_Loose.M(), weightLooseKinRecoTot);
        h_HypTTBarRapidity_Loose->Fill(hypttbar_Loose.Rapidity(), weightLooseKinRecoTot);
        h_HypTTBarRapidityAbs_Loose->Fill(std::abs(hypttbar_Loose.Rapidity()), weightLooseKinRecoTot);
        h_HypTTBarpT_Loose->Fill(hypttbar_Loose.Pt(), weightLooseKinRecoTot);

    }



    //...
    //=== CUT ===
    selectionStep = "8";


    //Require at least one solution for the kinematic event reconstruction
    if (!hasSolution) return kTRUE;

    //// obsolete? moved above
    // double weightLooseKinRecoTot = weight*weightLooseKinReco;

    weight *= weightKinReco;


    // FIXME Jenya:
    // Use variables accessed here everywhere in the following
    const KinematicReconstructionSolution solution = kinematicReconstructionSolutions.solution();
    const LV& hypTop = solution.top();
    const LV& hypAntiTop = solution.antiTop();
    //    const LV& hypTtbar = solution.ttbar();
    const LV& hypLepton = solution.lepton();
    const LV& hypAntiLepton = solution.antiLepton();
    const LV& hypBjet = solution.bjet();
    const LV& hypAntiBjet = solution.antiBjet();
    const LV& hypNeutrino = solution.neutrino();
    const LV& hypAntiNeutrino = solution.antiNeutrino();

    // Calculate the Top Scattering Angle
    TLorentzVector T_hyptop(common::LVtoTLV(kinematicReconstructionSolutions.solution().top()));
    TLorentzVector T_hypantitop(common::LVtoTLV(kinematicReconstructionSolutions.solution().antiTop()));
    TLorentzVector T_hypttbar(T_hyptop+T_hypantitop);

    TVector3 TTBarFrameBoost(-1. * T_hypttbar.BoostVector());
    TLorentzVector T_hyptop_ttbarframe = T_hyptop;
    T_hyptop_ttbarframe.Boost(TTBarFrameBoost);

    TVector3 ProtonBeamDirection(0., 0., 1.);
    double top_scatteringangle_ttbarframe = T_hyptop_ttbarframe.Vect().Unit().Dot(ProtonBeamDirection);
    // End Calculating Top Scattering Angle

    //Anya: MblblMET cut
    if( btopTag && ((met+hypBjet+hypAntiBjet+hypAntiLepton+hypLepton).M()/(hypTop+hypAntiTop).M() < 0.5)) return kTRUE;


    
    // Extra jet code moved to here: starts
    int extrarecojet[4] = {0,0,0,0};
    int jetnumReco = -1;
    double jetHTreco = 0;
    int RecoJets = 0, RecoJets_cut40 = 0, RecoJets_cut60 = 0, RecoJets_cut100 = 0;

    double cbin[20]={25.,35.,45.,55.,65.,75.,85.,95.,110.,130.,150.,170.,190.,210.,230.,250.,270.,300.,340.,380.};

    for(const int k : jetIndicesForExtraJetStudies){
      //std::cout <<k << "  TA " <<(*recoObjects.jets_).at(k).Eta()<< std::endl;
      //if(std::fabs(jets->at(k).Eta())>1.5 || std::fabs(jets->at(k).Eta())<0.8) continue;
      if(std::fabs((*recoObjects.jets_).at(k).Eta()) > 2.4) continue;//carmen eta cuts
      //        if((*recoObjects.jets_).at(k).Pt()< selections.jetPtCut_) {continue;}
      // changed to fill Mult***
      if ((*recoObjects.jets_).at(k).Pt()> selections.lead2JetPtCut_)
        {
	  RecoJets++;
	  if(std::fabs(hypAntiBjet.Pt() - (*recoObjects.jets_).at(k).Pt())>0.1 && std::fabs(hypBjet.Pt() - (*recoObjects.jets_).at(k).Pt())>0.1 && jetnumReco<3) {
	    jetHTreco+=(*recoObjects.jets_).at(k).Pt();
	    jetnumReco++;
	    extrarecojet[jetnumReco]= k;
	    //	    ExtraJetIndices.push_back(k);
	  }
        }
      if((*recoObjects.jets_).at(k).Pt()>40.) RecoJets_cut40++;
      if((*recoObjects.jets_).at(k).Pt()>60.) RecoJets_cut60++;
      if((*recoObjects.jets_).at(k).Pt()>100.) RecoJets_cut100++;
    }

    //    n_ExtraJets = jetnumReco+1;


    int first= -1, second=-1, third=-1, fourth=-1;
    double ptjet = 0;
    for(int ord = 0; ord <= jetnumReco; ord++)
      {
        if((*recoObjects.jets_).at(extrarecojet[ord]).Pt()> ptjet) {
	  first = ord; ptjet=(*recoObjects.jets_).at(extrarecojet[ord]).Pt();
        }
      }
    ptjet = 0;
    for(int ord = 0; ord <= jetnumReco && jetnumReco>0; ord++)
      {
	if((*recoObjects.jets_).at(extrarecojet[ord]).Pt()> ptjet && (*recoObjects.jets_).at(extrarecojet[ord]).Pt()< (*recoObjects.jets_).at(extrarecojet[first]).Pt()) {
	  second = ord; ptjet=(*recoObjects.jets_).at(extrarecojet[ord]).Pt();
        }
      }
    ptjet = 0;
    for(int ord = 0; ord <= jetnumReco && jetnumReco>1; ord++)
      {
        if((*recoObjects.jets_).at(extrarecojet[ord]).Pt()> ptjet && (*recoObjects.jets_).at(extrarecojet[ord]).Pt()< (*recoObjects.jets_).at(extrarecojet[second]).Pt() ) {
	  third = ord; ptjet=(*recoObjects.jets_).at(extrarecojet[ord]).Pt();
        }
      }
    ptjet = 0;
    for(int ord = 0; ord <= jetnumReco && jetnumReco>2; ord++)
      {
        if((*recoObjects.jets_).at(extrarecojet[ord]).Pt()> ptjet && (*recoObjects.jets_).at(extrarecojet[ord]).Pt()< (*recoObjects.jets_).at(extrarecojet[third]).Pt() ) {
	  fourth = ord; ptjet=(*recoObjects.jets_).at(extrarecojet[ord]).Pt();
        }
      }

       VLV leptonsForExtraJetCleaning;
       VLV topBQuarksForJetCleaning;
       leptonsForExtraJetCleaning.push_back(solution.lepton());
       leptonsForExtraJetCleaning.push_back(solution.antiLepton());
       topBQuarksForJetCleaning.push_back(solution.bjet());
       topBQuarksForJetCleaning.push_back(solution.antiBjet());

    selectIndices(jetIndicesForExtraJetStudies, jets, LVeta, selections.jetEtaCut_, false);
    selectIndices(jetIndicesForExtraJetStudies, jets, LVeta, -selections.jetEtaCut_);
    selectIndices(jetIndicesForExtraJetStudies, jets, LVpt, 40);
    if(recoObjects.jetPFID_){

      selectIndices(jetIndicesForExtraJetStudies, *recoObjects.jetPFID_, selections.jetPFIDWorkingPoint_);
    }
    else {
      throw std::runtime_error("TopAnalysis::Process -- null pointer in recoObjects.jetPFID_");
    }


    this->leptonCleanedJetIndices(jetIndicesForExtraJetStudies, jets, leptonsForExtraJetCleaning, selections.deltaRLeptonJetCut_);
    this->leptonCleanedJetIndices(jetIndicesForExtraJetStudies, jets, topBQuarksForJetCleaning, 0.8);
    ExtraJetIndicesIso08 = jetIndicesForExtraJetStudies;
    h_exjetMulti->Fill(ExtraJetIndicesIso08.size(), weight);
    //    extrajetCounts = {n_ExtraJets,ExtraJetIndicesIso04.size(), ExtraJetIndicesIso08.size() };
    // Extra jet code moved to here: ends
 

    this->fillAll(selectionStep,
                  eventMetadata,
                  recoObjects, commonGenObjects,
                  topGenObjects, topPseudoObjects,
                  kinematicReconstructionSolutions,
                  looseKinematicReconstructionSolution,
                  genObjectIndicesWithMatching, recoObjectIndices,
                  genLevelWeights, recoLevelWeights,
                  weight,weightLooseKinRecoTot);
    
    h_leptonPt_AfterKinReco->Fill((*recoObjects.allLeptons_).at(leptonIndex).Pt(), weight);
    h_leptonPt_AfterKinReco->Fill((*recoObjects.allLeptons_).at(antiLeptonIndex).Pt(), weight);
    h_leptonEta_AfterKinReco->Fill((*recoObjects.allLeptons_).at(leptonIndex).Eta(), weight);
    h_leptonEta_AfterKinReco->Fill((*recoObjects.allLeptons_).at(antiLeptonIndex).Eta(), weight);
    h_MET_step8->Fill(met.Pt(), weight);
    for (const int index : bjetIndices)
        h_bjeteta_AfterKinReco->Fill((*recoObjects.jets_).at(index).Eta(), weight);

    h_KinRecoSF->Fill(weightKinReco, 1);
    h_weights->Fill(weight, 1);
    h_EventWeight->Fill(weight, 1);

    h_vertMulti_step8->Fill(recoObjects.vertMulti_, weight);
    h_vertMulti_noPU_step8->Fill(recoObjects.vertMulti_, weightNoPileup);

    h_jetMultiXSec->Fill(numberOfJets, weight);
    h_jetMultiNoPU->Fill(numberOfJets, weightNoPileup );
    h_diLepMassFull_fullSel->Fill(dilepton.M(), weight);

    // Find 1st (and 2nd) leading pT particles: Top, Lepton, BJetIndex
    LV LeadHypTop, NLeadHypTop;
    LV LeadHypLepton, NLeadHypLepton;
    LV LeadHypBJet, NLeadHypBJet;
    orderLV(LeadHypTop, NLeadHypTop, hypTop, hypAntiTop, LVpt);
    orderLV(LeadHypLepton, NLeadHypLepton, hypLepton, hypAntiLepton, LVpt);
    orderLV(LeadHypBJet, NLeadHypBJet, hypBjet, hypAntiBjet, LVpt);

    //create ll, bb and tt system
    LV hypllbar(hypLepton + hypAntiLepton);
    LV hypbbbar(hypBjet + hypAntiBjet);
    LV hypttbar(hypTop+hypAntiTop);

    // create top/antitop quark in the ttbar rest frame
    ROOT::Math::Boost CoMBoostHypTtbar (hypttbar.BoostToCM());
    LV top = hypTop;
    LV antitop = hypAntiTop;
    top = CoMBoostHypTtbar(top);
    antitop = CoMBoostHypTtbar(antitop);

    //create lepton/antilepton in respective top rest frames
    LV lepton_rftbar = hypLepton;
    LV antilepton_rft =  hypAntiLepton;
    this->boostLeptonsToRespectiveTopQuarkRestFrames(lepton_rftbar, antilepton_rft, hypTop, hypAntiTop);

    // Histograms for resolution studies
    //if(topGenObjects.valuesSet_){
    //    h_RMSvsGenToppT->Fill((*topGenObjects.GenTop_).Pt(),(*topGenObjects.GenTop_).Pt()-hypTop.Pt());
    //    h_RMSvsGenToppT->Fill((*topGenObjects.GenAntiTop_).Pt(),(*topGenObjects.GenAntiTop_).Pt()-hypAntiTop.Pt());
    //    h_RMSvsGenTopRapidity->Fill((*topGenObjects.GenTop_).Rapidity(),(*topGenObjects.GenTop_).Rapidity()-hypTop.Rapidity());
    //    h_RMSvsGenTopRapidity->Fill((*topGenObjects.GenAntiTop_).Rapidity(),(*topGenObjects.GenAntiTop_).Rapidity()-hypAntiTop.Rapidity());
    //    h_RMSvsGenToppTLead->Fill(LeadGenTop.Pt(),LeadGenTop.Pt()-LeadHypTop.Pt());
    //    h_RMSvsGenToppTNLead->Fill(NLeadGenTop.Pt(),NLeadGenTop.Pt()-NLeadHypTop.Pt());
    //    h_RMSvsGenTopRapidityLead->Fill(LeadGenTop.Rapidity(),LeadGenTop.Rapidity()-LeadHypTop.Rapidity());
    //    h_RMSvsGenTopRapidityNLead->Fill(NLeadGenTop.Rapidity(),NLeadGenTop.Rapidity()-NLeadHypTop.Rapidity());
    //
    //
    //    h_RMSvsGenTTBarMass->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).M(),((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).M()-(hypTop+hypAntiTop).M());
    //    h_RMSvsGenTTBarpT->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Pt(),((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Pt()-(hypTop+hypAntiTop).Pt());
    //    h_RMSvsGenTTBarRapidity->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Rapidity(),((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Rapidity()-(hypTop+hypAntiTop).Rapidity());
    //    h_RMSvsGenTTBarDPhi->Fill(std::abs(DeltaPhi((*topGenObjects.GenTop_), (*topGenObjects.GenAntiTop_))),std::abs(DeltaPhi((*topGenObjects.GenTop_), (*topGenObjects.GenAntiTop_))) - std::abs(DeltaPhi(hypTop, hypAntiTop)));
    //    h_RMSvsGenTTBarDeltaRapidity->Fill(std::abs((*topGenObjects.GenTop_).Rapidity()) - std::abs((*topGenObjects.GenAntiTop_).Rapidity()), (std::abs((*topGenObjects.GenTop_).Rapidity()) -
    //                                       std::abs((*topGenObjects.GenAntiTop_).Rapidity())) - (std::abs(hypTop.Rapidity()) - std::abs(hypAntiTop.Rapidity())));
    //}

    //First fill the reco histograms (which have no scaling factors applied)
    const double recoWeight = trueLevelWeight;
    h_RecoTTBarDeltaRapidity->Fill(std::abs(hypTop.Rapidity()) - std::abs(hypAntiTop.Rapidity()), recoWeight);
    h_RecoTTBarDeltaPhi->Fill(std::abs(DeltaPhi(hypTop, hypAntiTop)), recoWeight);
    h_RecoTTBarMass->Fill(hypttbar.M(), recoWeight);
    h_RecoTTBarRapidity->Fill(hypttbar.Rapidity(), recoWeight);
    h_RecoTTBarRapidityAbs->Fill(std::abs(hypttbar.Rapidity()), recoWeight);
    h_RecoTTBarpT->Fill(hypttbar.Pt(), recoWeight);
    h_RecoToppT->Fill(hypTop.Pt(), recoWeight);
    h_RecoAntiToppT->Fill(hypAntiTop.Pt(), recoWeight);
    h_RecoTopRapidity->Fill(hypTop.Rapidity(), recoWeight);
    h_RecoAntiTopRapidity->Fill(hypAntiTop.Rapidity(), recoWeight);
    h_RecoTopRapidityAbs->Fill(std::abs(hypTop.Rapidity()), recoWeight);
    h_RecoAntiTopRapidityAbs->Fill(std::abs(hypAntiTop.Rapidity()), recoWeight);

    h_RecoToppTTTRestFrame->Fill(top.Pt(), recoWeight);
    h_RecoAntiToppTTTRestFrame->Fill(antitop.Pt(), recoWeight);

    h_RecoLLBarMass->Fill(hypllbar.M(), recoWeight);
    h_RecoLLBarpT->Fill(hypllbar.Pt(), recoWeight);
    h_RecoLeptonpT->Fill(hypLepton.Pt(), recoWeight);
    h_RecoAntiLeptonpT->Fill(hypAntiLepton.Pt(), recoWeight);
    h_RecoLeptonEta->Fill(hypLepton.Eta(), recoWeight);
    h_RecoAntiLeptonEta->Fill(hypAntiLepton.Eta(), recoWeight);

    h_RecoMet->Fill(met.Pt(), recoWeight);
    h_RecoHT->Fill(jetHT, recoWeight);

    h_RecoNeutrinopT->Fill(hypNeutrino.Pt(), recoWeight);
    h_RecoAntiNeutrinopT->Fill(hypAntiNeutrino.Pt(), recoWeight);

    h_RecoBBBarMass->Fill(hypbbbar.M(), recoWeight);
    h_RecoBBBarpT->Fill(hypbbbar.Pt(), recoWeight);
    h_RecoBBBarDPhi->Fill(std::fabs ( DeltaPhi ( hypBjet, hypAntiBjet ) ), recoWeight);

    h_RecoBJetpT->Fill(hypBjet.Pt(), recoWeight);
    h_RecoAntiBJetpT->Fill(hypAntiBjet.Pt(), recoWeight);
    h_RecoBJetRapidity->Fill(hypBjet.Rapidity(), recoWeight);
    h_RecoAntiBJetRapidity->Fill(hypAntiBjet.Rapidity(), recoWeight);
    h_RecoBJetEta->Fill(hypBjet.Eta(), recoWeight);
    h_RecoAntiBJetEta->Fill(hypAntiBjet.Eta(), recoWeight);

    h_RecoLLBarDPhi->Fill(std::fabs ( DeltaPhi ( hypLepton, hypAntiLepton ) ), recoWeight);
    h_RecoLLBarDEta->Fill(std::abs(hypLepton.Eta()) - std::abs(hypAntiLepton.Eta()), recoWeight);
    h_RecoLLBarCosPhiInTopRestFrame->Fill(std::cos ( Angle ( lepton_rftbar, antilepton_rft ) ), recoWeight);
    h_RecoLeptonantiBjetMass->Fill(( hypLepton+hypAntiBjet ).M(), recoWeight);
    h_RecoAntiLeptonBjetMass->Fill(( hypAntiLepton+hypBjet ).M(), recoWeight);

    h_RecoToppTLead->Fill(LeadHypTop.Pt(), recoWeight);
    h_RecoToppTNLead->Fill(NLeadHypTop.Pt(), recoWeight);
    h_RecoTopRapidityLead->Fill(LeadHypTop.Rapidity(), recoWeight);
    h_RecoTopRapidityNLead->Fill(NLeadHypTop.Rapidity(), recoWeight);
    h_RecoTopMassLead->Fill(LeadHypTop.M(), recoWeight);
    h_RecoTopMassNLead->Fill(NLeadHypTop.M(), recoWeight);

    h_RecoLeptonpTLead->Fill(LeadHypLepton.Pt(), recoWeight);
    h_RecoLeptonpTNLead->Fill(NLeadHypLepton.Pt(), recoWeight);
    h_RecoLeptonEtaLead->Fill(LeadHypLepton.Eta(), recoWeight);
    h_RecoLeptonEtaNLead->Fill(NLeadHypLepton.Eta(), recoWeight);

    h_RecoBJetpTLead->Fill(LeadHypBJet.Pt(), recoWeight);
    h_RecoBJetpTNLead->Fill(NLeadHypBJet.Pt(), recoWeight);
    h_RecoBJetEtaLead->Fill(LeadHypBJet.Eta(), recoWeight);
    h_RecoBJetEtaNLead->Fill(NLeadHypBJet.Eta(), recoWeight);

    //now go to the plots
    h_HypTTBarMass->Fill(hypttbar.M(), weight);
    h_HypTTBarRapidity->Fill(hypttbar.Rapidity(), weight);
    h_HypTTBarRapidityAbs->Fill(std::abs(hypttbar.Rapidity()), weight);
    h_HypTTBarpT->Fill(hypttbar.Pt(), weight);
    h_HypTTBarDeltaPhi->Fill(std::abs(DeltaPhi(hypTop, hypAntiTop)), weight);
    h_HypTTBarDeltaRapidity->Fill(std::abs(hypTop.Rapidity()) - std::abs(hypAntiTop.Rapidity()), weight);

    h_HypBBBarMass->Fill(hypbbbar.M(), weight);
    h_HypBBBarpT->Fill(hypbbbar.Pt(), weight);

    h_HypLLBarMass->Fill(hypllbar.M(), weight);
    h_HypLLBarpT->Fill(hypllbar.Pt(), weight);

    h_HypToppTTTRestFrame->Fill(top.Pt(), weight);
    h_HypAntiToppTTTRestFrame->Fill(antitop.Pt(), weight);

    h_HypMet->Fill(met.Pt(), weight);
    h_HypHT->Fill(jetHT, weight);

    h_HypTopMass->Fill(hypTop.M(), weight);
    h_HypAntiTopMass->Fill(hypAntiTop.M(), weight);
    h_HypToppT->Fill(hypTop.Pt(), weight);
    h_HypAntiToppT->Fill(hypAntiTop.Pt(), weight);

    h_HypLeptonpT->Fill(hypLepton.Pt(), weight);
    h_HypAntiLeptonpT->Fill(hypAntiLepton.Pt(), weight);

    h_HypBJetpT->Fill(hypBjet.Pt(), weight);
    h_HypAntiBJetpT->Fill(hypAntiBjet.Pt(), weight);
    h_HypBJetRapidity->Fill(hypBjet.Rapidity(), weight);
    h_HypAntiBJetRapidity->Fill(hypAntiBjet.Rapidity(), weight);
    h_HypBBBarDPhi->Fill(std::fabs ( DeltaPhi ( hypBjet, hypAntiBjet ) ), weight);

    h_HypBJetCSVdiscriminator->Fill(jetBtags.at(solution.bjetIndex()), weight);
    h_HypAntiBJetCSVdiscriminator->Fill(jetBtags.at(solution.antiBjetIndex()), weight);

    h_HypTopRapidity->Fill(hypTop.Rapidity(), weight);
    h_HypAntiTopRapidity->Fill(hypAntiTop.Rapidity(), weight);
    h_HypTopRapidityAbs->Fill(std::abs(hypTop.Rapidity()), weight);
    h_HypAntiTopRapidityAbs->Fill(std::abs(hypAntiTop.Rapidity()), weight);

    h_HypNeutrinopT->Fill(hypNeutrino.Pt(), weight);
    h_HypAntiNeutrinopT->Fill(hypAntiNeutrino.Pt(), weight);

    h_HypTopEta->Fill(hypTop.Eta(), weight);
    h_HypAntiTopEta->Fill(hypAntiTop.Eta(), weight);
    h_HypBJetEta->Fill(hypBjet.Eta(), weight);
    h_HypAntiBJetEta->Fill(hypAntiBjet.Eta(), weight);
    h_HypLeptonEta->Fill(hypLepton.Eta(), weight);

    h_HypAntiLeptonEta->Fill(hypAntiLepton.Eta(), weight);

    h_HypLLBarDPhi->Fill(std::fabs ( DeltaPhi ( hypLepton, hypAntiLepton ) ), weight);
    h_HypLLBarDEta->Fill(std::abs(hypLepton.Eta()) - std::abs(hypAntiLepton.Eta()), weight);
    h_HypLLBarCosPhiInTopRestFrame->Fill(std::cos ( Angle ( lepton_rftbar, antilepton_rft ) ), weight);

    h_HypLeptonantiBjetMass->Fill(( hypLepton + hypAntiBjet ).M(), weight);
    h_HypAntiLeptonBjetMass->Fill(( hypAntiLepton + hypBjet ).M(), weight);

    h_HypToppTLead->Fill(LeadHypTop.Pt(), weight);
    h_HypToppTNLead->Fill(NLeadHypTop.Pt(), weight);
    h_HypTopRapidityLead->Fill(LeadHypTop.Rapidity(), weight);
    h_HypTopRapidityNLead->Fill(NLeadHypTop.Rapidity(), weight);
    h_HypTopMassLead->Fill(LeadHypTop.M(), weight);
    h_HypTopMassNLead->Fill(NLeadHypTop.M(), weight);

    h_HypLeptonpTLead->Fill(LeadHypLepton.Pt(), weight);
    h_HypLeptonpTNLead->Fill(NLeadHypLepton.Pt(), weight);
    h_HypLeptonEtaLead->Fill(LeadHypLepton.Eta(), weight);
    h_HypLeptonEtaNLead->Fill(NLeadHypLepton.Eta(), weight);

    h_HypBJetpTLead->Fill(LeadHypBJet.Pt(), weight);
    h_HypBJetpTNLead->Fill(NLeadHypBJet.Pt(), weight);
    h_HypBJetEtaLead->Fill(LeadHypBJet.Eta(), weight);
    h_HypBJetEtaNLead->Fill(NLeadHypBJet.Eta(), weight);

    // Begin TTBar Mass 2D Distributions

    if (0 < hypttbar.M() && hypttbar.M() <= 450) {
      selectionStep = "8_0to450_ttbarmass";
      this->fillAll(selectionStep,
		    eventMetadata,
		    recoObjects, commonGenObjects,
		    topGenObjects, topPseudoObjects,
		    kinematicReconstructionSolutions,
		    looseKinematicReconstructionSolution,
		    genObjectIndicesWithMatching, recoObjectIndices,
		    genLevelWeights, recoLevelWeights,
		    weight,weightLooseKinRecoTot);
    }
    else if (450 < hypttbar.M() && hypttbar.M() <= 600) {
      selectionStep = "8_450to600_ttbarmass";
      this->fillAll(selectionStep,
		    eventMetadata,
		    recoObjects, commonGenObjects,
		    topGenObjects, topPseudoObjects,
		    kinematicReconstructionSolutions,
		    looseKinematicReconstructionSolution,
		    genObjectIndicesWithMatching, recoObjectIndices,
		    genLevelWeights, recoLevelWeights,
		    weight,weightLooseKinRecoTot);
    }
    else if (600 < hypttbar.M() && hypttbar.M() <= 800) {
      selectionStep = "8_600to800_ttbarmass";
      this->fillAll(selectionStep,
		    eventMetadata,
		    recoObjects, commonGenObjects,
		    topGenObjects, topPseudoObjects,
		    kinematicReconstructionSolutions,
		    looseKinematicReconstructionSolution,
		    genObjectIndicesWithMatching, recoObjectIndices,
		    genLevelWeights, recoLevelWeights,
		    weight,weightLooseKinRecoTot);
    }
    else if (800 < hypttbar.M()) {
      selectionStep = "8_800toinf_ttbarmass";
      this->fillAll(selectionStep,
		    eventMetadata,
		    recoObjects, commonGenObjects,
		    topGenObjects, topPseudoObjects,
		    kinematicReconstructionSolutions,
		    looseKinematicReconstructionSolution,
		    genObjectIndicesWithMatching, recoObjectIndices,
		    genLevelWeights, recoLevelWeights,
		    weight,weightLooseKinRecoTot);
    }

    // End TTBar Mass 2D Distributions


    // Begin Top Scattering Angle 2D Distributions

    if (-1.0 <= top_scatteringangle_ttbarframe && top_scatteringangle_ttbarframe <= -0.5) {
      selectionStep = "8_m1tomhalf_topscatteringangle";
      this->fillAll(selectionStep,
		    eventMetadata,
		    recoObjects, commonGenObjects,
		    topGenObjects, topPseudoObjects,
		    kinematicReconstructionSolutions,
		    looseKinematicReconstructionSolution,
		    genObjectIndicesWithMatching, recoObjectIndices,
		    genLevelWeights, recoLevelWeights,
		    weight,weightLooseKinRecoTot);
    }
    else if (-0.5 < top_scatteringangle_ttbarframe && top_scatteringangle_ttbarframe <= 0.0) {
      selectionStep = "8_mhalfto0_topscatteringangle";
      this->fillAll(selectionStep,
		    eventMetadata,
		    recoObjects, commonGenObjects,
		    topGenObjects, topPseudoObjects,
		    kinematicReconstructionSolutions,
		    looseKinematicReconstructionSolution,
		    genObjectIndicesWithMatching, recoObjectIndices,
		    genLevelWeights, recoLevelWeights,
		    weight,weightLooseKinRecoTot);
    }
    else if (0.0 < top_scatteringangle_ttbarframe && top_scatteringangle_ttbarframe <= 0.5) {
      selectionStep = "8_0tophalf_topscatteringangle";
      this->fillAll(selectionStep,
		    eventMetadata,
		    recoObjects, commonGenObjects,
		    topGenObjects, topPseudoObjects,
		    kinematicReconstructionSolutions,
		    looseKinematicReconstructionSolution,
		    genObjectIndicesWithMatching, recoObjectIndices,
		    genLevelWeights, recoLevelWeights,
		    weight,weightLooseKinRecoTot);
    }
    if (0.5 < top_scatteringangle_ttbarframe && top_scatteringangle_ttbarframe <= 1.0) {
      selectionStep = "8_phalftop1_topscatteringangle";
      this->fillAll(selectionStep,
		    eventMetadata,
		    recoObjects, commonGenObjects,
		    topGenObjects, topPseudoObjects,
		    kinematicReconstructionSolutions,
		    looseKinematicReconstructionSolution,
		    genObjectIndicesWithMatching, recoObjectIndices,
		    genLevelWeights, recoLevelWeights,
		    weight,weightLooseKinRecoTot);
    }
    // End Top Scattering Angle 2D Distributions

    // Begin Extra Jets 2D Distributions
    if (ExtraJetIndicesIso08.size() == 0) {
        selectionStep = "8_0_extrajets";
        this->fillAll(selectionStep,
                eventMetadata,
                recoObjects, commonGenObjects,
                topGenObjects, topPseudoObjects,
                kinematicReconstructionSolutions,
                looseKinematicReconstructionSolution,
                genObjectIndicesWithMatching, recoObjectIndices,
                genLevelWeights, recoLevelWeights,
                weight,weightLooseKinRecoTot);
    }
    else if (ExtraJetIndicesIso08.size() == 1) {
        selectionStep = "8_1_extrajets";
        this->fillAll(selectionStep,
                eventMetadata,
                recoObjects, commonGenObjects,
                topGenObjects, topPseudoObjects,
                kinematicReconstructionSolutions,
                looseKinematicReconstructionSolution,
                genObjectIndicesWithMatching, recoObjectIndices,
                genLevelWeights, recoLevelWeights,
                weight,weightLooseKinRecoTot);
    }
    else if (ExtraJetIndicesIso08.size() == 2) {
        selectionStep = "8_2_extrajets";
        this->fillAll(selectionStep,
                eventMetadata,
                recoObjects, commonGenObjects,
                topGenObjects, topPseudoObjects,
                kinematicReconstructionSolutions,
                looseKinematicReconstructionSolution,
                genObjectIndicesWithMatching, recoObjectIndices,
                genLevelWeights, recoLevelWeights,
                weight,weightLooseKinRecoTot);
    }
    else if (ExtraJetIndicesIso08.size() >=3 ) {
        selectionStep = "8_3ormore_extrajets";
        this->fillAll(selectionStep,
                eventMetadata,
                recoObjects, commonGenObjects,
                topGenObjects, topPseudoObjects,
                kinematicReconstructionSolutions,
                looseKinematicReconstructionSolution,
                genObjectIndicesWithMatching, recoObjectIndices,
                genLevelWeights, recoLevelWeights,
                weight,weightLooseKinRecoTot);
    }


    // End Extra Jets 2D Distributions

    //Ievgen

    h_HypTTBarRapidityvsTTBarpT->Fill(hypttbar.Pt(),hypttbar.Rapidity(),weight);

    // ...
    // Parton momentum fraction
    double RecoPartonMomFraction = (hypTop.energy() - hypTop.Pz() + hypAntiTop.energy() - hypAntiTop.Pz()) / (2 * 4000);
    double RecoAntipartonMomFraction = (hypTop.energy() + hypTop.Pz() + hypAntiTop.energy() + hypAntiTop.Pz()) / (2 * 4000);
    h_HypTopPartonFraction->Fill(RecoPartonMomFraction, weight);
    h_HypAntiTopPartonFraction->Fill(RecoAntipartonMomFraction, weight);
    h_RecoTopPartonFraction->Fill(RecoPartonMomFraction, recoWeight);
    h_RecoAntiTopPartonFraction->Fill(RecoAntipartonMomFraction, recoWeight);

    //New plots from Carmen: Begin
    h_RecoTTBar0Mass->Fill(rho0/(hypttbar).M(), recoWeight);
    h_HypTTBar0Mass->Fill(rho0/(hypttbar).M(), weight);
    h_jetMulti_kinReco->Fill(numberOfJets,weight);
    //h_bjetMulti_kinReco->Fill(numberOfBJets,weight);
    h_exjetMulti_kinReco->Fill(jetnumReco+1,weight);
    //    h_exjetMulti_0->Fill(jetnumReco+1,weight);
    //    h_exjetMulti_1->Fill(ExtraJetIndicesIso04.size(),weight);
    //    h_exjetMulti_2->Fill(ExtraJetIndicesIso08.size(),weight);

    for(int q0 = 0; q0<20;q0++){
      h_RecoJetMultTotal->Fill(cbin[q0],recoWeight);
      h_HypJetMultTotal->Fill(cbin[q0],weight);
      if((first >-1 && (*recoObjects.jets_).at(extrarecojet[first]).Pt()<= cbin[q0] )|| jetnumReco <0) {h_RecoJetMultQ0->Fill(cbin[q0],recoWeight);h_HypJetMultQ0->Fill(cbin[q0],weight);}
      if((second >-1 && (*recoObjects.jets_).at(extrarecojet[second]).Pt()<= cbin[q0]) || jetnumReco <1 ) {h_RecoJetExtra2Q0->Fill(cbin[q0],recoWeight);h_HypJetExtra2Q0->Fill(cbin[q0],weight);}
      if(jetHTreco<=cbin[q0]) {h_RecoJetMultQsum->Fill(cbin[q0],recoWeight); h_HypJetMultQsum->Fill(cbin[q0],weight);}
    }
    //     //New plots from Carmen: End

    h_RecoJetMultpt30->Fill(RecoJets,recoWeight);
    h_RecoJetMultpt40->Fill(RecoJets_cut40,recoWeight);
    h_RecoJetMultpt60->Fill(RecoJets_cut60,recoWeight);
    h_RecoJetMultpt100->Fill(RecoJets_cut100,recoWeight);
    h_HypJetMultpt30->Fill(RecoJets,weight);
    h_HypJetMultpt40->Fill(RecoJets_cut40,weight);
    h_HypJetMultpt60->Fill(RecoJets_cut60,weight);
    h_HypJetMultpt100->Fill(RecoJets_cut100,weight);
    

    //make sure you have called CreateBinnedControlPlots in the SlaveBegin first
    for (const auto& i : { hypTop, hypAntiTop } ) {
        FillBinnedControlPlot(h_HypToppT, i.Pt(), h_LeptonpT, (*recoObjects.allLeptons_).at(leptonIndex).Pt(), weight);
        FillBinnedControlPlot(h_HypToppT, i.Pt(), h_LeptonpT, (*recoObjects.allLeptons_).at(antiLeptonIndex).Pt(), weight);
        FillBinnedControlPlot(h_HypToppT, i.Pt(), h_diLepMassFull, dilepton.M(), weight);
        FillBinnedControlPlot(h_HypToppT, i.Pt(), h_LeptonEta, (*recoObjects.allLeptons_).at(leptonIndex).Eta(), weight);
        FillBinnedControlPlot(h_HypToppT, i.Pt(), h_LeptonEta, (*recoObjects.allLeptons_).at(antiLeptonIndex).Eta(), weight);
        FillBinnedControlPlot(h_HypToppT, i.Pt(), h_MET_preMETcut, met.Pt(), weight);

        FillBinnedControlPlot(h_HypTopRapidity, i.Rapidity(), h_LeptonpT, (*recoObjects.allLeptons_).at(leptonIndex).Pt(), weight);
        FillBinnedControlPlot(h_HypTopRapidity, i.Rapidity(), h_LeptonpT, (*recoObjects.allLeptons_).at(antiLeptonIndex).Pt(), weight);
        FillBinnedControlPlot(h_HypTopRapidity, i.Rapidity(), h_diLepMassFull, dilepton.M(), weight);
        FillBinnedControlPlot(h_HypTopRapidity, i.Rapidity(), h_LeptonEta, (*recoObjects.allLeptons_).at(leptonIndex).Eta(), weight);
        FillBinnedControlPlot(h_HypTopRapidity, i.Rapidity(), h_LeptonEta, (*recoObjects.allLeptons_).at(antiLeptonIndex).Eta(), weight);
        FillBinnedControlPlot(h_HypTopRapidity, i.Rapidity(), h_MET_preMETcut, met.Pt(), weight);
    }

    //=== CUT ===
    //Following histograms only filled for the signal sample
    if (!topGenObjects.valuesSet_) return kTRUE;


    // top quark properties
    h_GenRecoToppT->Fill(hypTop.Pt(), (*topGenObjects.GenTop_).Pt(), weight );
    h_GenRecoAntiToppT->Fill(hypAntiTop.Pt(), (*topGenObjects.GenAntiTop_).Pt(), weight );
    h_GenRecoToppTLead->Fill(LeadHypTop.Pt(), LeadGenTop.Pt(), weight);
    h_GenRecoTopRapidity->Fill(hypTop.Rapidity(), (*topGenObjects.GenTop_).Rapidity(), weight );
    h_GenRecoAntiTopRapidity->Fill(hypAntiTop.Rapidity(), (*topGenObjects.GenAntiTop_).Rapidity(), weight );
    h_GenRecoTopRapidityAbs->Fill(std::abs(hypTop.Rapidity()), std::abs((*topGenObjects.GenTop_).Rapidity()), weight );
    h_GenRecoAntiTopRapidityAbs->Fill(std::abs(hypAntiTop.Rapidity()), std::abs((*topGenObjects.GenAntiTop_).Rapidity()), weight );
    h_GenRecoToppTNLead->Fill(NLeadHypTop.Pt(), NLeadGenTop.Pt(), weight);
    h_GenRecoTopRapidityLead->Fill(LeadHypTop.Rapidity(), LeadGenTop.Rapidity(), weight);
    h_GenRecoTopRapidityNLead->Fill(NLeadHypTop.Rapidity(), NLeadGenTop.Rapidity(), weight);
    h_GenRecoTopMassLead->Fill(LeadHypTop.M(), LeadGenTop.M(), weight);
    h_GenRecoTopMassNLead->Fill(NLeadHypTop.M(), NLeadGenTop.M(), weight);

    // ttbar properties
    LV genttbar((*topGenObjects.GenTop_) + (*topGenObjects.GenAntiTop_));
    h_GenRecoTTBarMass->Fill(hypttbar.M(), genttbar.M(), weight );
    h_GenRecoTTBarpT->Fill(hypttbar.Pt(), genttbar.Pt(), weight );
    h_GenRecoTTBarRapidity->Fill(hypttbar.Rapidity(), genttbar.Rapidity(), weight );
    h_GenRecoTTBarRapidityAbs->Fill(std::abs(hypttbar.Rapidity()), std::abs(genttbar.Rapidity()), weight );
    h_GenRecoTTBarDeltaPhi->Fill(std::abs(DeltaPhi(hypTop, hypAntiTop)),
                                 std::abs(DeltaPhi((*topGenObjects.GenTop_), (*topGenObjects.GenAntiTop_))), weight);
    h_GenRecoTTBarDeltaRapidity->Fill(std::abs(hypTop.Rapidity()) - std::abs(hypAntiTop.Rapidity()),
                                     std::abs((*topGenObjects.GenTop_).Rapidity()) - std::abs((*topGenObjects.GenAntiTop_).Rapidity()),
                                      weight);


    // create top/antitop quark in the ttbar rest frame
    ROOT::Math::Boost CoMBoostGenTtbar (genttbar.BoostToCM());
    LV gentop ((*topGenObjects.GenTop_));
    LV genantitop ((*topGenObjects.GenAntiTop_));
    gentop = CoMBoostGenTtbar(gentop);
    genantitop = CoMBoostGenTtbar(genantitop);
    h_GenRecoToppTTTRestFrame->Fill(top.Pt(), gentop.Pt(), weight);
    h_GenRecoAntiToppTTTRestFrame->Fill(antitop.Pt(), genantitop.Pt(), weight);

    //create lepton/antilepton in respective top rest frames
    LV genlepton_rftbar = (*topGenObjects.GenLepton_);
    LV genantilepton_rft =  (*topGenObjects.GenAntiLepton_);
    this->boostLeptonsToRespectiveTopQuarkRestFrames(genlepton_rftbar, genantilepton_rft, (*topGenObjects.GenTop_), (*topGenObjects.GenAntiTop_));

    // Parton momentum fraction as defined by Olaf
    double GenPartonMomFraction = ((*topGenObjects.GenTop_).energy() - (*topGenObjects.GenTop_).Pz() + (*topGenObjects.GenAntiTop_).energy() - (*topGenObjects.GenAntiTop_).Pz()) / (2 * 4000);
    double GenAntipartonMomFraction = ((*topGenObjects.GenTop_).energy() + (*topGenObjects.GenTop_).Pz() + (*topGenObjects.GenAntiTop_).energy() + (*topGenObjects.GenAntiTop_).Pz()) / (2 * 4000);
    h_GenRecoTopPartonFraction->Fill(RecoPartonMomFraction,GenPartonMomFraction, weight);
    h_GenRecoAntiTopPartonFraction->Fill(RecoAntipartonMomFraction, GenAntipartonMomFraction, weight);


    // fill object in visible phase-space (relying on generator collections)

    if( (*topGenObjects.GenLepton_).Pt()> selections.leptonPtCut_ && std::fabs((*topGenObjects.GenLepton_).Eta()) < selections.leptonEtaCut_ &&
       (*topGenObjects.GenAntiLepton_).Pt()> selections.leptonPtCut_ && std::fabs((*topGenObjects.GenAntiLepton_).Eta()) < selections.leptonEtaCut_ &&
       BHadronIndex >= 0 && AntiBHadronIndex >= 0 && (*topGenObjects.allGenJets_).at(BHadronIndex).pt() > selections.jetPtCut_ && std::fabs ( (*topGenObjects.allGenJets_).at(BHadronIndex).eta() ) < selections.jetEtaCut_ && (*topGenObjects.allGenJets_).at(AntiBHadronIndex).pt() > selections.jetPtCut_ && std::fabs ( (*topGenObjects.allGenJets_).at(AntiBHadronIndex).Eta() ) < selections.jetEtaCut_ && std::fabs(DeltaR((*topGenObjects.GenLepton_), (*topGenObjects.allGenJets_).at(BHadronIndex))) > selections.genDeltaRLeptonJetCut_  && std::fabs(DeltaR((*topGenObjects.GenLepton_), (*topGenObjects.allGenJets_).at(AntiBHadronIndex))) > selections.genDeltaRLeptonJetCut_  && std::fabs(DeltaR((*topGenObjects.GenAntiLepton_), (*topGenObjects.allGenJets_).at(BHadronIndex))) > selections.genDeltaRLeptonJetCut_  && std::fabs(DeltaR((*topGenObjects.GenAntiLepton_), (*topGenObjects.allGenJets_).at(AntiBHadronIndex))) > selections.genDeltaRLeptonJetCut_)
            {
            // extra objects: met, ht, ...
            h_GenRecoMet->Fill(met.Pt(), (*topGenObjects.GenMet_).Pt(), weight);
            h_GenRecoHT->Fill(jetHT, genHT, weight);
            h_GenRecoNeutrinopT->Fill(hypNeutrino.Pt(), (*topGenObjects.GenNeutrino_).Pt(), weight);
            h_GenRecoAntiNeutrinopT->Fill(hypAntiNeutrino.Pt(), (*topGenObjects.GenAntiNeutrino_).Pt(), weight);

            // lepton distributions
            h_GenRecoLeptonpT->Fill(hypLepton.Pt(), (*topGenObjects.GenLepton_).Pt(), weight );
            h_GenRecoAntiLeptonpT->Fill(hypAntiLepton.Pt(), (*topGenObjects.GenAntiLepton_).Pt(), weight );
            h_GenRecoLeptonEta->Fill(hypLepton.Eta(), (*topGenObjects.GenLepton_).Eta(), weight );
            h_GenRecoAntiLeptonEta->Fill(hypAntiLepton.Eta(), (*topGenObjects.GenAntiLepton_).Eta(), weight );
            h_GenRecoLeptonpTLead->Fill(LeadHypLepton.Pt(), LeadGenLepton.Pt(), weight);
            h_GenRecoLeptonpTNLead->Fill(NLeadHypLepton.Pt(), NLeadGenLepton.Pt(), weight);
            h_GenRecoLeptonEtaLead->Fill(LeadHypLepton.Eta(), LeadGenLepton.Eta(), weight);
            h_GenRecoLeptonEtaNLead->Fill(NLeadHypLepton.Eta(), NLeadGenLepton.Eta(), weight);

            // lepton-pair distributions
            LV genllbar((*topGenObjects.GenLepton_) + (*topGenObjects.GenAntiLepton_));
            h_GenRecoLLBarMass->Fill(hypllbar.M(), genllbar.M(), weight );
            h_GenRecoLLBarpT->Fill(hypllbar.Pt(), genllbar.Pt(), weight );
            h_GenRecoLLBarDPhi->Fill(
                std::fabs( DeltaPhi( hypLepton, hypAntiLepton ) ),
                std::fabs( DeltaPhi( (*topGenObjects.GenLepton_), (*topGenObjects.GenAntiLepton_) ) ),
                weight );
            h_GenRecoLLBarDEta->Fill(
                std::abs(hypLepton.Eta()) - std::abs(hypAntiLepton.Eta()),
                std::abs((*topGenObjects.GenLepton_).Eta()) - std::abs((*topGenObjects.GenAntiLepton_).Eta()),
                weight );
            h_GenRecoLLBarCosPhiInTopRestFrame->Fill(
                std::cos( Angle( lepton_rftbar, antilepton_rft ) ),
                std::cos( Angle( genlepton_rftbar, genantilepton_rft ) ),
                weight );


            // letpon-b-jet mass
            h_GenRecoLeptonantiBjetMass->Fill(( hypLepton + hypAntiBjet ).M(), ( (*topGenObjects.GenLepton_)+(*topGenObjects.allGenJets_).at(AntiBHadronIndex) ).M(), weight );
            h_GenRecoAntiLeptonBjetMass->Fill(( hypAntiLepton+hypBjet ).M(), ( (*topGenObjects.GenAntiLepton_)+(*topGenObjects.allGenJets_).at(BHadronIndex) ).M(), weight );

            // b-jet-pair distributions
            LV genbbbar ((*topGenObjects.allGenJets_).at(BHadronIndex) + (*topGenObjects.allGenJets_).at(AntiBHadronIndex));
            h_GenRecoBBBarpT->Fill(hypbbbar.Pt(), genbbbar.Pt(), weight);
            h_GenRecoBBBarMass->Fill(hypbbbar.M(), genbbbar.M(), weight);
            h_GenRecoBBBarDPhi->Fill(
                std::fabs( DeltaPhi( hypBjet, hypAntiBjet ) ),
                std::fabs( DeltaPhi( (*topGenObjects.allGenJets_).at(BHadronIndex), (*topGenObjects.allGenJets_).at(AntiBHadronIndex) ) ),
                weight );

            // bjet distributions
            h_GenRecoBJetpT->Fill(hypBjet.Pt(), (*topGenObjects.allGenJets_).at(BHadronIndex).Pt(), weight );
            h_GenRecoAntiBJetpT->Fill(hypAntiBjet.Pt(), (*topGenObjects.allGenJets_).at(AntiBHadronIndex).Pt(), weight );
            h_GenRecoBJetEta->Fill(hypBjet.Eta(), (*topGenObjects.allGenJets_).at(BHadronIndex).Eta(), weight );
            h_GenRecoAntiBJetEta->Fill(hypAntiBjet.Eta(), (*topGenObjects.allGenJets_).at(AntiBHadronIndex).Eta(), weight );
            h_GenRecoBJetRapidity->Fill(hypBjet.Rapidity(), (*topGenObjects.allGenJets_).at(BHadronIndex).Rapidity(), weight );
            h_GenRecoAntiBJetRapidity->Fill(hypAntiBjet.Rapidity(), (*topGenObjects.allGenJets_).at(AntiBHadronIndex).Rapidity(), weight );
            h_GenRecoBJetpTLead->Fill(LeadHypBJet.Pt(), LeadGenBJet.Pt(), weight);
            h_GenRecoBJetpTNLead->Fill(NLeadHypBJet.Pt(), NLeadGenBJet.Pt(), weight);
            h_GenRecoBJetEtaLead->Fill(LeadHypBJet.Eta(), LeadGenBJet.Eta(), weight);
            h_GenRecoBJetEtaNLead->Fill(NLeadHypBJet.Eta(), NLeadGenBJet.Eta(), weight);

            //rho/masstt
            h_GenRecoTTBar0Mass->Fill(rho0/hypttbar.M(),rho0/genttbar.M(),weight);
            // jet multiplicities
            h_GenRecoJetMult->Fill(numberOfJets, (*topGenObjects.allGenJets_).size(), weight );
            h_GenRecoJetMultpt30->Fill(RecoJets,genJetIndices.size(),weight);
            h_GenRecoJetMultpt40->Fill(RecoJets_cut40,genJet40Indices.size(),weight);
            h_GenRecoJetMultpt60->Fill(RecoJets_cut60,genJet60Indices.size(),weight);
            h_GenRecoJetMultpt100->Fill(RecoJets_cut100,genJet100Indices.size(),weight);


            double ptaddgen = -1000.;
            double etaaddgen = -1000.;
            double DeltaRgen = -1000.;
            double massgen = -1000.;

            if(jetnumReco>2 && fourth != -1){
                if(genExtraJetIndices.size()>3) {ptaddgen = (*topGenObjects.allGenJets_).at(genExtraJetIndices.at(3)).Pt(); etaaddgen = (*topGenObjects.allGenJets_).at(genExtraJetIndices.at(3)).Eta();}

                h_GenRecoExtraJetpT4->Fill((*recoObjects.jets_).at(extrarecojet[fourth]).Pt(),ptaddgen,weight);
                h_GenRecoExtraJetEta4->Fill((*recoObjects.jets_).at(extrarecojet[fourth]).Eta(),etaaddgen,weight);
                if(std::fabs(etaaddgen) <2.4){
                    h_GenRecoExtraJetAbsEta4->Fill(std::fabs((*recoObjects.jets_).at(extrarecojet[fourth]).Eta()),std::fabs(etaaddgen),weight);
                 } else {h_GenRecoExtraJetAbsEta4->Fill(std::fabs((*recoObjects.jets_).at(extrarecojet[fourth]).Eta()),-1000.,weight);}
            }
            if (jetnumReco>1 && third != -1){
                if(genExtraJetIndices.size()>2) {ptaddgen = (*topGenObjects.allGenJets_).at(genExtraJetIndices.at(2)).Pt(); etaaddgen = (*topGenObjects.allGenJets_).at(genExtraJetIndices.at(2)).Eta();}
                h_GenRecoExtraJetpT3->Fill((*recoObjects.jets_).at(extrarecojet[third]).Pt(),ptaddgen,weight);
                h_GenRecoExtraJetEta3->Fill((*recoObjects.jets_).at(extrarecojet[third]).Eta(),etaaddgen,weight);
                if(std::fabs(etaaddgen) <2.4){
                    h_GenRecoExtraJetAbsEta3->Fill(std::fabs((*recoObjects.jets_).at(extrarecojet[third]).Eta()),std::fabs(etaaddgen),weight);
                } else {h_GenRecoExtraJetAbsEta3->Fill(std::fabs((*recoObjects.jets_).at(extrarecojet[third]).Eta()),-1000.,weight);}
            }
            if (jetnumReco>0 && second != -1){
                if(genExtraJetIndices.size()>1) {
                ptaddgen = (*topGenObjects.allGenJets_).at(genExtraJetIndices.at(1)).Pt(); etaaddgen = (*topGenObjects.allGenJets_).at(genExtraJetIndices.at(1)).Eta();
                DeltaRgen = std::fabs(DeltaR((*topGenObjects.allGenJets_).at(genExtraJetIndices.at(1)),(*topGenObjects.allGenJets_).at(genExtraJetIndices.at(0))));
                massgen = ((*topGenObjects.allGenJets_).at(genExtraJetIndices.at(1))+(*topGenObjects.allGenJets_).at(genExtraJetIndices.at(0))).M();
                }
                h_GenRecoExtraJetpT2->Fill((*recoObjects.jets_).at(extrarecojet[second]).Pt(),ptaddgen,weight);
                h_GenRecoExtraJetEta2->Fill((*recoObjects.jets_).at(extrarecojet[second]).Eta(),etaaddgen,weight);
                if(std::fabs(etaaddgen) <2.4){
                    h_GenRecoExtraJetAbsEta2->Fill(std::fabs((*recoObjects.jets_).at(extrarecojet[second]).Eta()),std::fabs(etaaddgen),weight);
                } else {h_GenRecoExtraJetAbsEta2->Fill(std::fabs((*recoObjects.jets_).at(extrarecojet[second]).Eta()),-1000.,weight);}

                h_GenRecoDeltaRExtraJet12->Fill(std::fabs(DeltaR((*recoObjects.jets_).at(extrarecojet[second]),(*recoObjects.jets_).at(extrarecojet[first]))),DeltaRgen,weight);
                h_GenRecoMassExtraJet12->Fill(((*recoObjects.jets_).at(extrarecojet[second])+(*recoObjects.jets_).at(extrarecojet[first])).M(),massgen, weight);
                if((*recoObjects.jets_).at(extrarecojet[first]).Pt()< 60. && (*recoObjects.jets_).at(extrarecojet[first]).Eta() >0.0 &&
                (*recoObjects.jets_).at(extrarecojet[second]).Pt()< 60. && (*recoObjects.jets_).at(extrarecojet[second]).Eta() < 0.0 &&
                TMath::Abs(DeltaR((*recoObjects.jets_).at(extrarecojet[second]),(*recoObjects.jets_).at(extrarecojet[first]))) > 0.6 )
                {
                    if(genExtraJetIndices.size()>1) {
                    h_GenRecoDeltaPhiExtraJet12->Fill(DeltaPhi((*recoObjects.jets_).at(extrarecojet[first]),(*recoObjects.jets_).at(extrarecojet[second])),DeltaPhi((*topGenObjects.allGenJets_).at(genExtraJetIndices.at(1)),(*topGenObjects.allGenJets_).at(genExtraJetIndices.at(0))),weight);
                    h_GenRecoPhiExtraJet12->Fill((*recoObjects.jets_).at(extrarecojet[first]).Phi()+(*recoObjects.jets_).at(extrarecojet[second]).Phi(),(*topGenObjects.allGenJets_).at(genExtraJetIndices.at(0)).Phi()+(*topGenObjects.allGenJets_).at(genExtraJetIndices.at(1)).Phi(),weight);
                    } else {
                    h_GenRecoDeltaPhiExtraJet12->Fill(DeltaPhi((*recoObjects.jets_).at(extrarecojet[first]),(*recoObjects.jets_).at(extrarecojet[second])),-1000,weight);
                    h_GenRecoPhiExtraJet12->Fill((*recoObjects.jets_).at(extrarecojet[first]).Phi()+(*recoObjects.jets_).at(extrarecojet[second]).Phi(),-1000,weight);
                    }
                }
            }
            if (jetnumReco >-1 && first != -1){
                 if(genExtraJetIndices.size()>0) {
                    ptaddgen = (*topGenObjects.allGenJets_).at(genExtraJetIndices.at(0)).Pt(); etaaddgen = (*topGenObjects.allGenJets_).at(genExtraJetIndices.at(0)).Eta();
                    h_GenRecoTTBar1stJetMass->Fill(rho0/(hypttbar+(*recoObjects.jets_).at(extrarecojet[first])).M(),rho0/(genttbar+(*topGenObjects.allGenJets_).at(genExtraJetIndices.at(0))).M(),weight);
                    h_GenRecoExtraJetHT->Fill(jetHTreco,jetHTGen,weight);
               } else {
                    h_GenRecoExtraJetHT->Fill(jetHTreco,-1000.,weight);
                    h_GenRecoTTBar1stJetMass->Fill(rho0/(hypttbar+(*recoObjects.jets_).at(extrarecojet[first])).M(),-1000.,weight);
               }
                h_GenRecoExtraJetpT->Fill((*recoObjects.jets_).at(extrarecojet[first]).Pt(),ptaddgen,weight);
                h_GenRecoExtraJetEta->Fill((*recoObjects.jets_).at(extrarecojet[first]).Eta(),etaaddgen,weight);
                if(std::fabs(etaaddgen) <2.4){
                    h_GenRecoExtraJetAbsEta->Fill(std::fabs((*recoObjects.jets_).at(extrarecojet[first]).Eta()),std::fabs(etaaddgen),weight);
                } else {h_GenRecoExtraJetAbsEta->Fill(std::fabs((*recoObjects.jets_).at(extrarecojet[first]).Eta()),-1000.,weight);}

    //                h_GenRecoTTBar1stJetMass->Fill(rho0/(hypttbar+(*recoObjects.jets_).at(extrarecojet[first])).M(),rho0/(genttbar+(*topGenObjects.allGenJets_).at(extragenjet[0])).M(),weight);
            }

    }
    else{// fill underflow/overflow for reco objects not in vis. phase space
        // extra objects: met, ht, ...
        // check h_GenRecoExtraJetHT->Fill(jetHTreco,-1000.,weight);
        h_GenRecoMet->Fill(met.Pt(), -1000, weight);
        h_GenRecoHT->Fill(jetHT, -1000, weight);
        h_GenRecoNeutrinopT->Fill(hypNeutrino.Pt(), -1000, weight);
        h_GenRecoAntiNeutrinopT->Fill(hypAntiNeutrino.Pt(), -1000, weight);

        // lepton distributions
        h_GenRecoLeptonpT->Fill(hypLepton.Pt(), -1000, weight );
        h_GenRecoAntiLeptonpT->Fill(hypAntiLepton.Pt(), -1000, weight );
        h_GenRecoLeptonEta->Fill(hypLepton.Eta(), -1000, weight );
        h_GenRecoAntiLeptonEta->Fill(hypAntiLepton.Eta(), -1000, weight );
        h_GenRecoLeptonpTLead->Fill(LeadHypLepton.Pt(), -1000, weight);
        h_GenRecoLeptonpTNLead->Fill(NLeadHypLepton.Pt(), -1000, weight);
        h_GenRecoLeptonEtaLead->Fill(LeadHypLepton.Eta(), -1000, weight);
        h_GenRecoLeptonEtaNLead->Fill(NLeadHypLepton.Eta(), -1000, weight);

        // lepton-pair distributions
        h_GenRecoLLBarMass->Fill(hypllbar.M(), -1000, weight );
        h_GenRecoLLBarpT->Fill(hypllbar.Pt(), -1000, weight );
        h_GenRecoLLBarDPhi->Fill(
            std::fabs( DeltaPhi( hypLepton, hypAntiLepton ) ), -1000, weight );
        h_GenRecoLLBarDEta->Fill(
            std::abs(hypLepton.Eta()) - std::abs(hypAntiLepton.Eta()), -1000, weight );
        h_GenRecoLLBarCosPhiInTopRestFrame->Fill(
                std::cos( Angle( lepton_rftbar, antilepton_rft ) ), -1000, weight );

        // letpon-b-jet mass
        h_GenRecoLeptonantiBjetMass->Fill(( hypLepton + hypAntiBjet ).M(), -1000, weight );
        h_GenRecoAntiLeptonBjetMass->Fill(( hypAntiLepton+hypBjet ).M(), -1000, weight );

        // b-jet-pair distributions
        h_GenRecoBBBarpT->Fill(hypbbbar.Pt(), -1000, weight);
        h_GenRecoBBBarMass->Fill(hypbbbar.M(), -1000, weight);
        h_GenRecoBBBarDPhi->Fill(
            std::fabs( DeltaPhi( hypBjet, hypAntiBjet ) ), -1000, weight );

        // bjet distributions
        h_GenRecoBJetpT->Fill(hypBjet.Pt(), -1000, weight );
        h_GenRecoAntiBJetpT->Fill(hypAntiBjet.Pt(), -1000, weight );
        h_GenRecoBJetEta->Fill(hypBjet.Eta(), -1000, weight );
        h_GenRecoAntiBJetEta->Fill(hypAntiBjet.Eta(), -1000, weight );
        h_GenRecoBJetRapidity->Fill(hypBjet.Rapidity(), -1000, weight );
        h_GenRecoAntiBJetRapidity->Fill(hypAntiBjet.Rapidity(), -1000, weight );
        h_GenRecoBJetpTLead->Fill(LeadHypBJet.Pt(), -1000, weight);
        h_GenRecoBJetpTNLead->Fill(NLeadHypBJet.Pt(), -1000, weight);
        h_GenRecoBJetEtaLead->Fill(LeadHypBJet.Eta(), -1000, weight);
        h_GenRecoBJetEtaNLead->Fill(NLeadHypBJet.Eta(), -1000, weight);

        //rho/masstt
        h_GenRecoTTBar0Mass->Fill(rho0/hypttbar.M(), -1000, weight);

        // jet multiplicities
        h_GenRecoJetMult->Fill(numberOfJets, -1000, weight );
        h_GenRecoJetMultpt30->Fill(RecoJets,-1000,weight);
        h_GenRecoJetMultpt40->Fill(RecoJets_cut40,-1000,weight);
        h_GenRecoJetMultpt60->Fill(RecoJets_cut60,-1000,weight);
        h_GenRecoJetMultpt100->Fill(RecoJets_cut100,-1000,weight);
        if(jetnumReco>2 && fourth != -1){
            h_GenRecoExtraJetpT4->Fill((*recoObjects.jets_).at(extrarecojet[fourth]).Pt(),-1000.,weight);
            h_GenRecoExtraJetEta4->Fill((*recoObjects.jets_).at(extrarecojet[fourth]).Eta(),-1000.,weight);
            h_GenRecoExtraJetAbsEta4->Fill(std::fabs((*recoObjects.jets_).at(extrarecojet[fourth]).Eta()),-1000.,weight);
        }
        if(jetnumReco>1 && third != -1){
            h_GenRecoExtraJetpT3->Fill((*recoObjects.jets_).at(extrarecojet[third]).Pt(),-1000.,weight);
            h_GenRecoExtraJetEta3->Fill((*recoObjects.jets_).at(extrarecojet[third]).Eta(),-1000.,weight);
            h_GenRecoExtraJetAbsEta3->Fill(std::fabs((*recoObjects.jets_).at(extrarecojet[third]).Eta()),-1000.,weight);
        }
        if (jetnumReco>0 && second != -1){
            h_GenRecoExtraJetpT2->Fill((*recoObjects.jets_).at(extrarecojet[second]).Pt(),-1000.,weight);
            h_GenRecoExtraJetEta2->Fill((*recoObjects.jets_).at(extrarecojet[second]).Eta(),-1000.,weight);
            h_GenRecoExtraJetAbsEta2->Fill(std::fabs((*recoObjects.jets_).at(extrarecojet[second]).Eta()),-1000.,weight);
            h_GenRecoDeltaRExtraJet12->Fill(std::fabs(DeltaR((*recoObjects.jets_).at(extrarecojet[second]),(*recoObjects.jets_).at(extrarecojet[first]))),-1000.,weight);
            h_GenRecoMassExtraJet12->Fill(((*recoObjects.jets_).at(extrarecojet[second])+(*recoObjects.jets_).at(extrarecojet[first])).M(),-1000.,weight);
            h_GenRecoDeltaPhiExtraJet12->Fill(DeltaPhi((*recoObjects.jets_).at(extrarecojet[first]),(*recoObjects.jets_).at(extrarecojet[second])),-1000,weight);
            h_GenRecoPhiExtraJet12->Fill((*recoObjects.jets_).at(extrarecojet[first]).Phi()+(*recoObjects.jets_).at(extrarecojet[second]).Phi(),-1000,weight);
        }
        if (jetnumReco >-1 && first != -1){
            h_GenRecoExtraJetpT->Fill((*recoObjects.jets_).at(extrarecojet[first]).Pt(),-1000.,weight);
            h_GenRecoExtraJetEta->Fill((*recoObjects.jets_).at(extrarecojet[first]).Eta(),-1000.,weight);
            h_GenRecoExtraJetAbsEta->Fill(std::fabs((*recoObjects.jets_).at(extrarecojet[first]).Eta()),-1000.,weight);
            h_GenRecoExtraJetHT->Fill(jetHTreco,-1000.,weight);
            h_GenRecoTTBar1stJetMass->Fill(rho0/(hypttbar+(*recoObjects.jets_).at(extrarecojet[first])).M(),-1000.,weight);
        }
    }

    if (analysisConfig_.sampleComposition().pseudoTopMode_ == 0) return kTRUE;

    //=======================PSEUDO-TOP: BEGIN==========================//
    //Following histograms only filled for the signal sample
    if (!topPseudoObjects.valuesSet_) return kTRUE;

    // fill object in pseudo visible phase-space

    LV pseudodilepton((*topPseudoObjects.PseudoLepton_) + (*topPseudoObjects.PseudoAntiLepton_));

    if( (*topPseudoObjects.PseudoLepton_).Pt()> selections.leptonPtCut_ && std::fabs((*topPseudoObjects.PseudoLepton_).Eta()) < selections.leptonEtaCut_ &&
        (*topPseudoObjects.PseudoAntiLepton_).Pt()> selections.leptonPtCut_ && std::fabs((*topPseudoObjects.PseudoAntiLepton_).Eta()) < selections.leptonEtaCut_ &&
        (*topPseudoObjects.PseudoBJet_).pt() > selections.jetPtCut_ && std::fabs ( (*topPseudoObjects.PseudoBJet_).eta() ) < selections.jetEtaCut_ &&
        (*topPseudoObjects.PseudoAntiBJet_).pt() > selections.jetPtCut_ && std::fabs ( (*topPseudoObjects.PseudoAntiBJet_).Eta() ) < selections.jetEtaCut_ &&
        std::fabs(DeltaR((*topPseudoObjects.PseudoLepton_), (*topPseudoObjects.PseudoBJet_))) > selections.genDeltaRLeptonJetCut_  &&
        std::fabs(DeltaR((*topPseudoObjects.PseudoLepton_), (*topPseudoObjects.PseudoAntiBJet_))) > selections.genDeltaRLeptonJetCut_  &&
        std::fabs(DeltaR((*topPseudoObjects.PseudoAntiLepton_), (*topPseudoObjects.PseudoBJet_))) > selections.genDeltaRLeptonJetCut_  &&
        std::fabs(DeltaR((*topPseudoObjects.PseudoAntiLepton_), (*topPseudoObjects.PseudoAntiBJet_))) > selections.genDeltaRLeptonJetCut_ &&
        pseudodilepton.M() > 20.0)
   {
           // top quark properties
           h_PseudoRecoToppT->Fill(hypTop.Pt(), (*topPseudoObjects.PseudoTop_).Pt(), weight );
           h_PseudoRecoAntiToppT->Fill(hypAntiTop.Pt(), (*topPseudoObjects.PseudoAntiTop_).Pt(), weight );
           h_PseudoRecoTopRapidity->Fill(hypTop.Rapidity(), (*topPseudoObjects.PseudoTop_).Rapidity(), weight );
           h_PseudoRecoAntiTopRapidity->Fill(hypAntiTop.Rapidity(), (*topPseudoObjects.PseudoAntiTop_).Rapidity(), weight );
           h_PseudoRecoTopRapidityAbs->Fill(std::abs(hypTop.Rapidity()), std::abs((*topPseudoObjects.PseudoTop_).Rapidity()), weight );
           h_PseudoRecoAntiTopRapidityAbs->Fill(std::abs(hypAntiTop.Rapidity()), std::abs((*topPseudoObjects.PseudoAntiTop_).Rapidity()), weight );

           h_PseudoRecoToppTLead->Fill(LeadHypTop.Pt(), LeadPseudoTop.Pt(), weight);
           h_PseudoRecoToppTNLead->Fill(NLeadHypTop.Pt(), NLeadPseudoTop.Pt(), weight);
           h_PseudoRecoTopRapidityLead->Fill(LeadHypTop.Rapidity(), LeadPseudoTop.Rapidity(), weight);
           h_PseudoRecoTopRapidityNLead->Fill(NLeadHypTop.Rapidity(), NLeadPseudoTop.Rapidity(), weight);

           // ttbar properties
           LV pseudottbar((*topPseudoObjects.PseudoTop_) + (*topPseudoObjects.PseudoAntiTop_));
           h_PseudoRecoTTBarMass->Fill(hypttbar.M(), pseudottbar.M(), weight );
           h_PseudoRecoTTBarpT->Fill(hypttbar.Pt(), pseudottbar.Pt(), weight );
           h_PseudoRecoTTBarRapidity->Fill(hypttbar.Rapidity(), pseudottbar.Rapidity(), weight );
           h_PseudoRecoTTBarRapidityAbs->Fill(std::abs(hypttbar.Rapidity()), std::abs(pseudottbar.Rapidity()), weight );
           h_PseudoRecoTTBarDeltaPhi->Fill(std::abs(DeltaPhi(hypTop, hypAntiTop)), std::abs(DeltaPhi((*topPseudoObjects.PseudoTop_), (*topPseudoObjects.PseudoAntiTop_))), weight);
           h_PseudoRecoTTBarDeltaRapidity->Fill(std::abs(hypTop.Rapidity()) - std::abs(hypAntiTop.Rapidity()), std::abs((*topPseudoObjects.PseudoTop_).Rapidity()) - std::abs((*topPseudoObjects.PseudoAntiTop_).Rapidity()), weight);

           // create top/antitop quark in the ttbar rest frame
           ROOT::Math::Boost CoMBoostPseudoTtbar (pseudottbar.BoostToCM());
           LV pseudotop ((*topPseudoObjects.PseudoTop_));
           LV pseudoantitop ((*topPseudoObjects.PseudoAntiTop_));
           pseudotop = CoMBoostPseudoTtbar(pseudotop);
           pseudoantitop = CoMBoostPseudoTtbar(pseudoantitop);
           h_PseudoRecoToppTTTRestFrame->Fill(top.Pt(), pseudotop.Pt(), weight);
           h_PseudoRecoAntiToppTTTRestFrame->Fill(antitop.Pt(), pseudoantitop.Pt(), weight);

           //create lepton/antilepton in respective top rest frames
           LV pseudolepton_rftbar = (*topPseudoObjects.PseudoLepton_);
           LV pseudoantilepton_rft =  (*topPseudoObjects.PseudoAntiLepton_);
           this->boostLeptonsToRespectiveTopQuarkRestFrames(pseudolepton_rftbar, pseudoantilepton_rft, (*topPseudoObjects.PseudoTop_), (*topPseudoObjects.PseudoAntiTop_));

           // jet multiplicities
           h_PseudoRecoJetMultpt30->Fill(RecoJets,pseudoJetIndices.size(),weight);
           h_PseudoRecoJetMultpt40->Fill(RecoJets_cut40,pseudoJet40Indices.size(),weight);
           h_PseudoRecoJetMultpt60->Fill(RecoJets_cut60,pseudoJet60Indices.size(),weight);
           h_PseudoRecoJetMultpt100->Fill(RecoJets_cut100,pseudoJet100Indices.size(),weight);

           // lepton distributions
           h_PseudoRecoLeptonpT->Fill(hypLepton.Pt(), (*topPseudoObjects.PseudoLepton_).Pt(), weight );
           h_PseudoRecoAntiLeptonpT->Fill(hypAntiLepton.Pt(), (*topPseudoObjects.PseudoAntiLepton_).Pt(), weight );
           h_PseudoRecoLeptonEta->Fill(hypLepton.Eta(), (*topPseudoObjects.PseudoLepton_).Eta(), weight );
           h_PseudoRecoAntiLeptonEta->Fill(hypAntiLepton.Eta(), (*topPseudoObjects.PseudoAntiLepton_).Eta(), weight );
           h_PseudoRecoLeptonpTLead->Fill(LeadHypLepton.Pt(), LeadPseudoLepton.Pt(), weight);
           h_PseudoRecoLeptonpTNLead->Fill(NLeadHypLepton.Pt(), NLeadPseudoLepton.Pt(), weight);
           h_PseudoRecoLeptonEtaLead->Fill(LeadHypLepton.Eta(), LeadPseudoLepton.Eta(), weight);
           h_PseudoRecoLeptonEtaNLead->Fill(NLeadHypLepton.Eta(), NLeadPseudoLepton.Eta(), weight);

           // lepton-pair distributions
           LV pseudollbar((*topPseudoObjects.PseudoLepton_) + (*topPseudoObjects.PseudoAntiLepton_));
           h_PseudoRecoLLBarMass->Fill(hypllbar.M(), pseudollbar.M(), weight );
           h_PseudoRecoLLBarpT->Fill(hypllbar.Pt(), pseudollbar.Pt(), weight );
           h_PseudoRecoLLBarDPhi->Fill(
               std::fabs( DeltaPhi( hypLepton, hypAntiLepton ) ),
               std::fabs( DeltaPhi( (*topPseudoObjects.PseudoLepton_), (*topPseudoObjects.PseudoAntiLepton_) ) ),
               weight );
           h_PseudoRecoLLBarDEta->Fill(
               std::abs(hypLepton.Eta()) - std::abs(hypAntiLepton.Eta()),
               std::abs((*topPseudoObjects.PseudoLepton_).Eta()) - std::abs((*topPseudoObjects.PseudoAntiLepton_).Eta()),
               weight );
           h_PseudoRecoLLBarCosPhiInTopRestFrame->Fill(
               std::cos( Angle( lepton_rftbar, antilepton_rft ) ),
               std::cos( Angle( pseudolepton_rftbar, pseudoantilepton_rft ) ),
               weight );

           // b-jet-pair distributions
           LV pseudobbbar ((*topPseudoObjects.PseudoBJet_) + (*topPseudoObjects.PseudoAntiBJet_));
           h_PseudoRecoBBBarpT->Fill(hypbbbar.Pt(), pseudobbbar.Pt(), weight);
           h_PseudoRecoBBBarMass->Fill(hypbbbar.M(), pseudobbbar.M(), weight);
           h_PseudoRecoBBBarDPhi->Fill(
               std::fabs( DeltaPhi( hypBjet, hypAntiBjet ) ),
               std::fabs( DeltaPhi( (*topPseudoObjects.PseudoBJet_), (*topPseudoObjects.PseudoAntiBJet_) ) ),
               weight );

           // bjet distributions
           h_PseudoRecoBJetpT->Fill(hypBjet.Pt(), (*topPseudoObjects.PseudoBJet_).Pt(), weight );
           h_PseudoRecoAntiBJetpT->Fill(hypAntiBjet.Pt(), (*topPseudoObjects.PseudoAntiBJet_).Pt(), weight );
           h_PseudoRecoBJetEta->Fill(hypBjet.Eta(), (*topPseudoObjects.PseudoBJet_).Eta(), weight );
           h_PseudoRecoAntiBJetEta->Fill(hypAntiBjet.Eta(), (*topPseudoObjects.PseudoAntiBJet_).Eta(), weight );
           h_PseudoRecoBJetRapidity->Fill(hypBjet.Rapidity(), (*topPseudoObjects.PseudoBJet_).Rapidity(), weight );
           h_PseudoRecoAntiBJetRapidity->Fill(hypAntiBjet.Rapidity(), (*topPseudoObjects.PseudoAntiBJet_).Rapidity(), weight );
           h_PseudoRecoBJetpTLead->Fill(LeadHypBJet.Pt(), LeadPseudoBJet.Pt(), weight);
           h_PseudoRecoBJetpTNLead->Fill(NLeadHypBJet.Pt(), NLeadPseudoBJet.Pt(), weight);
           h_PseudoRecoBJetEtaLead->Fill(LeadHypBJet.Eta(), LeadPseudoBJet.Eta(), weight);
           h_PseudoRecoBJetEtaNLead->Fill(NLeadHypBJet.Eta(), NLeadPseudoBJet.Eta(), weight);

           // letpon-b-jet mass
           h_PseudoRecoLeptonantiBjetMass->Fill(( hypLepton + hypAntiBjet ).M(), ( (*topPseudoObjects.PseudoLepton_)+ (*topPseudoObjects.PseudoAntiBJet_) ).M(), weight );
           h_PseudoRecoAntiLeptonBjetMass->Fill(( hypAntiLepton+hypBjet ).M(), ( (*topPseudoObjects.PseudoAntiLepton_)+ (*topPseudoObjects.PseudoBJet_) ).M(), weight );

           //Gen-pseudo matrices
           h_GenPseudoToppT->Fill((*topPseudoObjects.PseudoTop_).Pt(), (*topGenObjects.GenTop_).Pt(), weight );
           h_GenPseudoAntiToppT->Fill((*topPseudoObjects.PseudoAntiTop_).Pt(), (*topGenObjects.GenAntiTop_).Pt(), weight );
           h_GenPseudoTopRapidity->Fill((*topPseudoObjects.PseudoTop_).Rapidity(), (*topGenObjects.GenTop_).Rapidity(), weight );
           h_GenPseudoAntiTopRapidity->Fill((*topPseudoObjects.PseudoAntiTop_).Rapidity(), (*topGenObjects.GenAntiTop_).Rapidity(), weight );
           h_GenPseudoTTBarRapidity->Fill(pseudottbar.Rapidity(), genttbar.Rapidity(), weight );
           h_GenPseudoTTBarpT->Fill(pseudottbar.Pt(), genttbar.Pt(), weight );
           h_GenPseudoTTBarMass->Fill(pseudottbar.M(), genttbar.M(), weight );

           h_GenPseudoLeptonpT->Fill((*topPseudoObjects.PseudoLepton_).Pt(), (*topGenObjects.GenLepton_).Pt(), weight );
           h_GenPseudoAntiLeptonpT->Fill((*topPseudoObjects.PseudoAntiLepton_).Pt(), (*topGenObjects.GenAntiLepton_).Pt(), weight );
           h_GenPseudoLeptonEta->Fill((*topPseudoObjects.PseudoLepton_).Eta(), (*topGenObjects.GenLepton_).Eta(), weight );
           h_GenPseudoAntiLeptonEta->Fill((*topPseudoObjects.PseudoAntiLepton_).Eta(), (*topGenObjects.GenAntiLepton_).Eta(), weight );

           h_GenPseudoJetMultpt30->Fill(pseudoJetIndices.size(),genJetIndices.size(),weight);

           //Resolution
           if((*topGenObjects.GenTop_).Pt()!=0.0) h_GenPseudoToppT_reso->Fill(((*topGenObjects.GenTop_).Pt() - (*topPseudoObjects.PseudoTop_).Pt())/(*topGenObjects.GenTop_).Pt(), weight );
           if((*topGenObjects.GenAntiTop_).Pt()!=0.0) h_GenPseudoAntiToppT_reso->Fill(((*topGenObjects.GenAntiTop_).Pt() - (*topPseudoObjects.PseudoAntiTop_).Pt())/(*topGenObjects.GenAntiTop_).Pt(), weight );
           if((*topGenObjects.GenTop_).Rapidity()!=0.0) h_GenPseudoTopRapidity_reso->Fill(((*topGenObjects.GenTop_).Rapidity() - (*topPseudoObjects.PseudoTop_).Rapidity())/std::fabs((*topGenObjects.GenTop_).Rapidity()), weight );
           if((*topGenObjects.GenAntiTop_).Rapidity()!=0.0) h_GenPseudoAntiTopRapidity_reso->Fill(((*topGenObjects.GenAntiTop_).Rapidity() - (*topPseudoObjects.PseudoAntiTop_).Rapidity())/
                    std::fabs((*topGenObjects.GenAntiTop_).Rapidity()), weight );
           if(std::fabs(genttbar.Rapidity())!=0.0) h_GenPseudoTTBarRapidity_reso->Fill((genttbar.Rapidity() - pseudottbar.Rapidity())/std::fabs(genttbar.Rapidity()), weight );
           if(genttbar.Pt()!=0.0) h_GenPseudoTTBarpT_reso->Fill((genttbar.Pt() - pseudottbar.Pt())/genttbar.Pt(), weight );
           if(genttbar.M()!=0.0) h_GenPseudoTTBarMass_reso->Fill((genttbar.M() - pseudottbar.M())/genttbar.M(), weight );

           if((*topGenObjects.GenLepton_).Pt()!=0.0) h_GenPseudoLeptonpT_reso->Fill(((*topGenObjects.GenLepton_).Pt() - (*topPseudoObjects.PseudoLepton_).Pt())/(*topGenObjects.GenLepton_).Pt(), weight );
           if((*topGenObjects.GenAntiLepton_).Pt()!=0.0) h_GenPseudoAntiLeptonpT_reso->Fill(((*topGenObjects.GenAntiLepton_).Pt() - (*topPseudoObjects.PseudoAntiLepton_).Pt())/(*topGenObjects.GenAntiLepton_).Pt(), weight );
           if(std::fabs((*topGenObjects.GenLepton_).Eta())!=0.0) h_GenPseudoLeptonEta_reso->Fill(((*topGenObjects.GenLepton_).Eta() - (*topPseudoObjects.PseudoLepton_).Eta())/std::fabs((*topGenObjects.GenLepton_).Eta()), weight );
           if(std::fabs((*topGenObjects.GenAntiLepton_).Eta())!=0.0) h_GenPseudoAntiLeptonEta_reso->Fill(((*topGenObjects.GenAntiLepton_).Eta() - (*topPseudoObjects.PseudoAntiLepton_).Eta())/
                    std::fabs((*topGenObjects.GenAntiLepton_).Eta()), weight );

           if(genJetIndices.size()!=0.0) h_GenPseudoJetMultpt30_reso->Fill((genJetIndices.size() - pseudoJetIndices.size())/genJetIndices.size(),weight);

           if(BHadronIndex >= 0 && AntiBHadronIndex >= 0 && BHadronIndex!=AntiBHadronIndex) {
               h_GenPseudoBJetpT->Fill((*topPseudoObjects.PseudoBJet_).Pt(), (*topGenObjects.allGenJets_).at(BHadronIndex).Pt(), weight );
               h_GenPseudoAntiBJetpT->Fill((*topPseudoObjects.PseudoAntiBJet_).Pt(), (*topGenObjects.allGenJets_).at(AntiBHadronIndex).Pt(), weight );
               h_GenPseudoBJetEta->Fill((*topPseudoObjects.PseudoBJet_).Eta(), (*topGenObjects.allGenJets_).at(BHadronIndex).Eta(), weight );
               h_GenPseudoAntiBJetEta->Fill((*topPseudoObjects.PseudoAntiBJet_).Eta(), (*topGenObjects.allGenJets_).at(AntiBHadronIndex).Eta(), weight );

               if((*topGenObjects.allGenJets_).at(BHadronIndex).Pt()!=0.0) h_GenPseudoBJetpT_reso->Fill(((*topGenObjects.allGenJets_).at(BHadronIndex).Pt() -
                     (*topPseudoObjects.PseudoBJet_).Pt())/(*topGenObjects.allGenJets_).at(BHadronIndex).Pt() , weight );
               if((*topGenObjects.allGenJets_).at(AntiBHadronIndex).Pt()!=0.0) h_GenPseudoAntiBJetpT_reso->Fill(((*topGenObjects.allGenJets_).at(AntiBHadronIndex).Pt() -
                     (*topPseudoObjects.PseudoAntiBJet_).Pt())/(*topGenObjects.allGenJets_).at(AntiBHadronIndex).Pt(), weight );
               if(std::fabs((*topGenObjects.allGenJets_).at(BHadronIndex).Eta())!=0.0) h_GenPseudoBJetEta_reso->Fill(((*topGenObjects.allGenJets_).at(BHadronIndex).Eta() -
                     (*topPseudoObjects.PseudoBJet_).Eta())/std::fabs((*topGenObjects.allGenJets_).at(BHadronIndex).Eta()), weight );
               if(std::fabs((*topGenObjects.allGenJets_).at(AntiBHadronIndex).Eta())!=0.0) h_GenPseudoAntiBJetEta_reso->Fill(((*topGenObjects.allGenJets_).at(AntiBHadronIndex).Eta() -
                     (*topPseudoObjects.PseudoAntiBJet_).Eta())/std::fabs((*topGenObjects.allGenJets_).at(AntiBHadronIndex).Eta()), weight );
           }

   } else if (analysisConfig_.sampleComposition().pseudoTopMode_ == 2 &&
        (*topPseudoObjects.PseudoLepton_).Pt()> 18. && std::fabs((*topPseudoObjects.PseudoLepton_).Eta()) < 2.5 &&
        (*topPseudoObjects.PseudoAntiLepton_).Pt()> 18. && std::fabs((*topPseudoObjects.PseudoAntiLepton_).Eta()) < 2.5 &&
        (*topPseudoObjects.PseudoBJet_).pt() > 27. && std::fabs ( (*topPseudoObjects.PseudoBJet_).eta() ) < 2.8 &&
        (*topPseudoObjects.PseudoAntiBJet_).pt() > 27. && std::fabs ( (*topPseudoObjects.PseudoAntiBJet_).Eta() ) < 2.8 &&
        std::fabs(DeltaR((*topPseudoObjects.PseudoLepton_), (*topPseudoObjects.PseudoBJet_))) > selections.genDeltaRLeptonJetCut_  &&
        std::fabs(DeltaR((*topPseudoObjects.PseudoLepton_), (*topPseudoObjects.PseudoAntiBJet_))) > selections.genDeltaRLeptonJetCut_  &&
        std::fabs(DeltaR((*topPseudoObjects.PseudoAntiLepton_), (*topPseudoObjects.PseudoBJet_))) > selections.genDeltaRLeptonJetCut_  &&
        std::fabs(DeltaR((*topPseudoObjects.PseudoAntiLepton_), (*topPseudoObjects.PseudoAntiBJet_))) > selections.genDeltaRLeptonJetCut_ &&
        pseudodilepton.M() > 18.0 ) {// fill underflow for reco objects not in vis. phase space
            // top quark properties
            h_PseudoRecoToppT->Fill(hypTop.Pt(), -1000, weight );
            h_PseudoRecoAntiToppT->Fill(hypAntiTop.Pt(), -1000, weight );
            h_PseudoRecoTopRapidity->Fill(hypTop.Rapidity(), -1000, weight );
            h_PseudoRecoAntiTopRapidity->Fill(hypAntiTop.Rapidity(), -1000, weight );
            h_PseudoRecoTopRapidityAbs->Fill(std::abs(hypTop.Rapidity()), -1000, weight );
            h_PseudoRecoAntiTopRapidityAbs->Fill(std::abs(hypAntiTop.Rapidity()), -1000, weight );

            h_PseudoRecoToppTLead->Fill(LeadHypTop.Pt(), -1000, weight);
            h_PseudoRecoToppTNLead->Fill(NLeadHypTop.Pt(), -1000, weight);
            h_PseudoRecoTopRapidityLead->Fill(LeadHypTop.Rapidity(), -1000, weight);
            h_PseudoRecoTopRapidityNLead->Fill(NLeadHypTop.Rapidity(), -1000, weight);

            h_PseudoRecoToppTTTRestFrame->Fill(top.Pt(), -1000, weight);
            h_PseudoRecoAntiToppTTTRestFrame->Fill(antitop.Pt(), -1000, weight);

            // ttbar properties
            h_PseudoRecoTTBarMass->Fill(hypttbar.M(), -1000, weight );
            h_PseudoRecoTTBarpT->Fill(hypttbar.Pt(), -1000, weight );
            h_PseudoRecoTTBarRapidity->Fill(hypttbar.Rapidity(), -1000, weight );
            h_PseudoRecoTTBarRapidityAbs->Fill(std::abs(hypttbar.Rapidity()), -1000, weight );
            h_PseudoRecoTTBarDeltaPhi->Fill(std::abs(DeltaPhi(hypTop, hypAntiTop)), -1000, weight);
            h_PseudoRecoTTBarDeltaRapidity->Fill(std::abs(hypTop.Rapidity()) - std::abs(hypAntiTop.Rapidity()), -1000, weight);

            // jet multiplicities
            h_PseudoRecoJetMultpt30->Fill(RecoJets, -1000, weight);
            h_PseudoRecoJetMultpt40->Fill(RecoJets_cut40, -1000, weight);
            h_PseudoRecoJetMultpt60->Fill(RecoJets_cut60, -1000, weight);
            h_PseudoRecoJetMultpt100->Fill(RecoJets_cut100, -1000, weight);

            // lepton distributions
            h_PseudoRecoLeptonpT->Fill(hypLepton.Pt(), -1000, weight );
            h_PseudoRecoAntiLeptonpT->Fill(hypAntiLepton.Pt(), -1000, weight );
            h_PseudoRecoLeptonEta->Fill(hypLepton.Eta(), -1000, weight );
            h_PseudoRecoAntiLeptonEta->Fill(hypAntiLepton.Eta(), -1000, weight );
            h_PseudoRecoLeptonpTLead->Fill(LeadHypLepton.Pt(), -1000, weight);
            h_PseudoRecoLeptonpTNLead->Fill(NLeadHypLepton.Pt(), -1000, weight);
            h_PseudoRecoLeptonEtaLead->Fill(LeadHypLepton.Eta(), -1000, weight);
            h_PseudoRecoLeptonEtaNLead->Fill(NLeadHypLepton.Eta(), -1000, weight);

            // lepton-pair distributions
            h_PseudoRecoLLBarMass->Fill(hypllbar.M(), -1000, weight );
            h_PseudoRecoLLBarpT->Fill(hypllbar.Pt(), -1000, weight );
            h_PseudoRecoLLBarDPhi->Fill( std::fabs( DeltaPhi( hypLepton, hypAntiLepton ) ), -1000, weight );
            h_PseudoRecoLLBarDEta->Fill( std::abs(hypLepton.Eta()) - std::abs(hypAntiLepton.Eta()), -1000, weight );
            h_PseudoRecoLLBarCosPhiInTopRestFrame->Fill( std::cos( Angle( lepton_rftbar, antilepton_rft ) ), -1000, weight );

            // b-jet-pair distributions
            h_PseudoRecoBBBarpT->Fill(hypbbbar.Pt(), -1000, weight);
            h_PseudoRecoBBBarMass->Fill(hypbbbar.M(), -1000, weight);
            h_PseudoRecoBBBarDPhi->Fill( std::fabs( DeltaPhi( hypBjet, hypAntiBjet ) ), -1000, weight );

            // bjet distributions
            h_PseudoRecoBJetpT->Fill(hypBjet.Pt(), -1000, weight );
            h_PseudoRecoAntiBJetpT->Fill(hypAntiBjet.Pt(), -1000, weight );
            h_PseudoRecoBJetEta->Fill(hypBjet.Eta(), -1000, weight );
            h_PseudoRecoAntiBJetEta->Fill(hypAntiBjet.Eta(), -1000, weight );
            h_PseudoRecoBJetRapidity->Fill(hypBjet.Rapidity(), -1000, weight );
            h_PseudoRecoAntiBJetRapidity->Fill(hypAntiBjet.Rapidity(), -1000, weight );
            h_PseudoRecoBJetpTLead->Fill(LeadHypBJet.Pt(), -1000, weight);
            h_PseudoRecoBJetpTNLead->Fill(NLeadHypBJet.Pt(), -1000, weight);
            h_PseudoRecoBJetEtaLead->Fill(LeadHypBJet.Eta(), -1000, weight);
            h_PseudoRecoBJetEtaNLead->Fill(NLeadHypBJet.Eta(), -1000, weight);

            // letpon-b-jet mass
            h_PseudoRecoLeptonantiBjetMass->Fill(( hypLepton + hypAntiBjet ).M(), -1000, weight );
            h_PseudoRecoAntiLeptonBjetMass->Fill(( hypAntiLepton+hypBjet ).M(), -1000, weight );

    } else if (analysisConfig_.sampleComposition().pseudoTopMode_ == 3) {// fill underflow for reco objects not in vis. phase space
            // top quark properties
            h_HypToppT_OutOfSpace->Fill(hypTop.Pt(), weight );
            h_HypAntiToppT_OutOfSpace->Fill(hypAntiTop.Pt(), weight );
            h_HypTopRapidity_OutOfSpace->Fill(hypTop.Rapidity(), weight );
            h_HypAntiTopRapidity_OutOfSpace->Fill(hypAntiTop.Rapidity(), weight );
            h_HypTopRapidityAbs_OutOfSpace->Fill(std::abs(hypTop.Rapidity()), weight );
            h_HypAntiTopRapidityAbs_OutOfSpace->Fill(std::abs(hypAntiTop.Rapidity()), weight );

            h_HypToppTLead_OutOfSpace->Fill(LeadHypTop.Pt(), weight);
            h_HypToppTNLead_OutOfSpace->Fill(NLeadHypTop.Pt(), weight);
            h_HypTopRapidityLead_OutOfSpace->Fill(LeadHypTop.Rapidity(), weight);
            h_HypTopRapidityNLead_OutOfSpace->Fill(NLeadHypTop.Rapidity(), weight);

            h_HypToppTTTRestFrame_OutOfSpace->Fill(top.Pt(), weight);
            h_HypAntiToppTTTRestFrame_OutOfSpace->Fill(antitop.Pt(), weight);

            // ttbar properties
            h_HypTTBarMass_OutOfSpace->Fill(hypttbar.M(), weight );
            h_HypTTBarpT_OutOfSpace->Fill(hypttbar.Pt(), weight );
            h_HypTTBarRapidity_OutOfSpace->Fill(hypttbar.Rapidity(), weight );
            h_HypTTBarRapidityAbs_OutOfSpace->Fill(std::abs(hypttbar.Rapidity()), weight );
            h_HypTTBarDeltaPhi_OutOfSpace->Fill(std::abs(DeltaPhi(hypTop, hypAntiTop)), weight);
            h_HypTTBarDeltaRapidity_OutOfSpace->Fill(std::abs(hypTop.Rapidity()) - std::abs(hypAntiTop.Rapidity()), weight);

            // jet multiplicities
            h_HypJetMultpt30_OutOfSpace->Fill(RecoJets, weight);
            h_HypJetMultpt40_OutOfSpace->Fill(RecoJets_cut40, weight);
            h_HypJetMultpt60_OutOfSpace->Fill(RecoJets_cut60, weight);
            h_HypJetMultpt100_OutOfSpace->Fill(RecoJets_cut100, weight);

            // lepton distributions
            h_HypLeptonpT_OutOfSpace->Fill(hypLepton.Pt(), weight );
            h_HypAntiLeptonpT_OutOfSpace->Fill(hypAntiLepton.Pt(), weight );
            h_HypLeptonEta_OutOfSpace->Fill(hypLepton.Eta(), weight );
            h_HypAntiLeptonEta_OutOfSpace->Fill(hypAntiLepton.Eta(), weight );
            h_HypLeptonpTLead_OutOfSpace->Fill(LeadHypLepton.Pt(), weight);
            h_HypLeptonpTNLead_OutOfSpace->Fill(NLeadHypLepton.Pt(), weight);
            h_HypLeptonEtaLead_OutOfSpace->Fill(LeadHypLepton.Eta(), weight);
            h_HypLeptonEtaNLead_OutOfSpace->Fill(NLeadHypLepton.Eta(), weight);

            // lepton-pair distributions
            h_HypLLBarMass_OutOfSpace->Fill(hypllbar.M(), weight );
            h_HypLLBarpT_OutOfSpace->Fill(hypllbar.Pt(), weight );
            h_HypLLBarDPhi_OutOfSpace->Fill( std::fabs( DeltaPhi( hypLepton, hypAntiLepton ) ), weight );
            h_HypLLBarDEta_OutOfSpace->Fill( std::abs(hypLepton.Eta()) - std::abs(hypAntiLepton.Eta()) , weight );
            h_HypLLBarCosPhiInTopRestFrame_OutOfSpace->Fill( std::cos( Angle( lepton_rftbar, antilepton_rft ) ), weight );

            // b-jet-pair distributions
            h_HypBBBarpT_OutOfSpace->Fill(hypbbbar.Pt(), weight);
            h_HypBBBarMass_OutOfSpace->Fill(hypbbbar.M(), weight);
            h_HypBBBarDPhi_OutOfSpace->Fill( std::fabs( DeltaPhi( hypBjet, hypAntiBjet ) ), weight );

            // bjet distributions
            h_HypBJetpT_OutOfSpace->Fill(hypBjet.Pt(), weight );
            h_HypAntiBJetpT_OutOfSpace->Fill(hypAntiBjet.Pt(), weight );
            h_HypBJetEta_OutOfSpace->Fill(hypBjet.Eta(), weight );
            h_HypAntiBJetEta_OutOfSpace->Fill(hypAntiBjet.Eta(), weight );
            h_HypBJetRapidity_OutOfSpace->Fill(hypBjet.Rapidity(), weight );
            h_HypAntiBJetRapidity_OutOfSpace->Fill(hypAntiBjet.Rapidity(), weight );
            h_HypBJetpTLead_OutOfSpace->Fill(LeadHypBJet.Pt(), weight);
            h_HypBJetpTNLead_OutOfSpace->Fill(NLeadHypBJet.Pt(), weight);
            h_HypBJetEtaLead_OutOfSpace->Fill(LeadHypBJet.Eta(), weight);
            h_HypBJetEtaNLead_OutOfSpace->Fill(NLeadHypBJet.Eta(), weight);

            // letpon-b-jet mass
            h_HypLeptonantiBjetMass_OutOfSpace->Fill(( hypLepton + hypAntiBjet ).M(), weight );
            h_HypAntiLeptonBjetMass_OutOfSpace->Fill(( hypAntiLepton+hypBjet ).M(), weight );

    }

    //=======================PSEUDO-TOP: END==========================//



    // ...

    return kTRUE;
}



const TopGenObjects& TopAnalysis::getTopGenObjects(const Long64_t& entry)const
{
    if((!isTopSignal_)||(isTtbarZSample_)) return *topGenObjects_;
    if(topGenObjects_->valuesSet_) return *topGenObjects_;

    this->GetTopSignalBranchesEntry(entry);
    return *topGenObjects_;
}


void TopAnalysis::SetRunViaTau(const bool runViaTau)
{
    runViaTau_ = runViaTau;

    //    std::cout << "   runViaTau_: " << runViaTau_ << "\n";

    // If running signalviatau, set Top Signal to true
    if(runViaTau) this->SetTopSignal(true);

    // If running bgviatau, set Top Signal to false
    //    if(runViaTau) this->SetTopSignal(false);
}



void TopAnalysis::SetClosureTest(TString closure, double slope)
{
    if (closure == "") {
        doClosureTest_ = false;
    } else {
        doClosureTest_ = true;
        if (closure == "pttop") {
            closureFunction_ = [&,slope](Long64_t entry) -> double {
                const TopGenObjects& topGenObjects = this->getTopGenObjects(entry);
                return std::max((1.+(topGenObjects.GenTop_->Pt()-100.)*slope)
                               *(1.+(topGenObjects.GenAntiTop_->Pt()-100.)*slope) , 0.1);
            };
        } else if (closure == "ytop") {
            closureFunction_ = [&,slope](Long64_t entry) -> double {
                const TopGenObjects& topGenObjects = this->getTopGenObjects(entry);
                return std::max((1.+(std::fabs(topGenObjects.GenTop_->Rapidity())-1.)*slope)
                               *(1.+(std::fabs(topGenObjects.GenAntiTop_->Rapidity()-1.))*slope) , 0.1);
            };
        } else if (closure == "nominal") {
            closureFunction_ = [](Long64_t) -> double {return 1.;};
        } else {
            std::cerr << "invalid closure test function\n";
            exit(1);
        }

        TString outputFilename = this->outputFilename();
        if (closure != "nominal") {
            outputFilename.ReplaceAll(".root", TString::Format("_fakerun_%s%.3f.root", closure.Data(), slope));
        } else {
            outputFilename.ReplaceAll(".root", TString::Format("_fakerun_%s.root", closure.Data()));
        }
        this->SetOutputfilename(outputFilename);
        std::cout << "<<< Closure test. Writing to: " << outputFilename << "\n";

        //BRANCHING FRACTION
        double br = 0.;
        const int channelPdgIdProduct = this->channelPdgIdProduct();
        if (channelPdgIdProduct == -11*11) br = 0.01147;
        else if (channelPdgIdProduct == -13*13) br = 0.01130;
        else if (channelPdgIdProduct == -11*13) br = 0.02277;
        else {
            std::cerr << "closure test channel invalid\n"; exit(1);
        }
        closureMaxEvents_ = TOPXSEC * analysisConfig_.general().luminosity_ * br;
        TString samplename = this->samplename();
        samplename.Append("_fakedata");
        this->SetSamplename(samplename);
    }
}



double TopAnalysis::calculateClosureTestWeight(const Long64_t& entry)
{
    if(!doClosureTest_ || !this->isTopSignal()) return 1.;
    const double weight = closureFunction_(entry);
    h_ClosureTotalWeight->Fill(1, weight);
    return weight;
}



bool TopAnalysis::failsTopGeneratorSelection(const Long64_t& entry)const
{
    if(!this->isTtbarPlusTauSample()) return false;
    // topDecayMode contains the decay of the top (*10) + the decay of the antitop (plus 100 or 200 for non-b decays of tops)
    // 1=hadron, 2=e, 3=mu, 4=tau->hadron, 5=tau->e, 6=tau->mu
    // i.e. 23 == top decays to e, tbar decays to mu
    const int topDecayMode = this->topDecayMode(entry) % 100;
    const bool isViaTau = topDecayMode > 40 || (topDecayMode % 10 > 4);
    bool isCorrectChannel(false);
    switch(this->channelPdgIdProduct()){
        case -11*13:
            isCorrectChannel = topDecayMode == 23 || topDecayMode == 32 //emu prompt
                            || topDecayMode == 53 || topDecayMode == 35 //e via tau, mu prompt
                            || topDecayMode == 26 || topDecayMode == 62 //e prompt, mu via tau
                            || topDecayMode == 56 || topDecayMode == 65; //both via tau
            break;
        case -11*11:
            isCorrectChannel = topDecayMode == 22  //ee prompt
                            || topDecayMode == 52 || topDecayMode == 25 //e prompt, e via tau
                            || topDecayMode == 55; //both via tau
            break;
        case -13*13:
            isCorrectChannel = topDecayMode == 33  //mumu prompt
                            || topDecayMode == 36 || topDecayMode == 63 //mu prompt, mu via tau
                            || topDecayMode == 66; //both via tau
            break;
        default:
            std::cerr<<"Invalid channel in failsTopGeneratorSelection()! Product = "<<this->channelPdgIdProduct()
                     <<"\n...break\n"<<std::endl;
            exit(213);
    };

    /*
    if ( topDecayMode == 23 || topDecayMode == 32 )    std::cout << "  emu prompt  " << "\n";
    if ( topDecayMode == 53 || topDecayMode == 35 )    std::cout << "  e via tau, mu prompt " << "\n";
    if ( topDecayMode == 26 || topDecayMode == 62 )    std::cout << "  e prompt, mu via tau " << "\n";
    if ( topDecayMode == 56 || topDecayMode == 65 )    std::cout << "  both emu via tau " << "\n";
    if ( topDecayMode == 22 )    std::cout << "  both ee prompt " << "\n";
    if ( topDecayMode == 52 || topDecayMode == 25 )    std::cout << "  e prompt, e via tau  " << "\n";
    if ( topDecayMode == 55 )    std::cout << "  both ee via tau " << "\n";
    if ( topDecayMode == 33 )    std::cout << "  both mumu prompt  " << "\n";
    if ( topDecayMode == 36 || topDecayMode == 63 )    std::cout << " mu prompt, mu via tau " << "\n";
    if ( topDecayMode == 66 )    std::cout << "  both mumu via tau " << "\n";
    */

    if(!isCorrectChannel) {
      //std::cout << "   fails TopGeneratorSelection because top decay mode is the wrong channel" << "\n";
      return true;
    }
    else if(runViaTau_ != isViaTau) {
      //std::cout << "   fails TopGeneratorSelection because decay is via tau and the mode is not via tau" << "\n";
      return true;
    }
    else {
      //std::cout << "   does not fail TopGeneratorSelection" << "\n";
      return false;
    }
}



void TopAnalysis::generatorTopEvent(LV& leadGenTop, LV& nLeadGenTop,
                                    LV& leadGenLepton, LV& nLeadGenLepton,
                                    LV& leadGenBJet, LV& nLeadGenBJet,
                                    double& genHT,
                                    const int bHadronIndex, const int antiBHadronIndex,
                                    const double trueLevelWeightNoPileup, const double trueLevelWeight, const TopGenObjects& topGenObjects)
{
    // Use utilities without namespaces
    using namespace common;
    using ROOT::Math::VectorUtil::DeltaPhi;
    using ROOT::Math::VectorUtil::Angle;

    LV& LeadGenTop(leadGenTop);
    LV& NLeadGenTop(nLeadGenTop);
    LV& LeadGenLepton(leadGenLepton);
    LV& NLeadGenLepton(nLeadGenLepton);
    LV& LeadGenBJet(leadGenBJet);
    LV& NLeadGenBJet(nLeadGenBJet);

    genHT = -1.;

    if(!topGenObjects.valuesSet_) return;

    // Access object selections from config
    const AnalysisConfig::Selections& selections = analysisConfig_.selections();

    const int BHadronIndex(bHadronIndex);
    const int AntiBHadronIndex(antiBHadronIndex);

    h_GenAll->Fill((*topGenObjects.GenTop_).M(), trueLevelWeight);
    h_GenAll_noweight->Fill((*topGenObjects.GenTop_).M(), trueLevelWeightNoPileup);

    //Begin: Select & Fill histograms with Leading pT and 2nd Leading pT: Lepton and BJet
    orderLV(LeadGenLepton, NLeadGenLepton, (*topGenObjects.GenLepton_), (*topGenObjects.GenAntiLepton_), LVpt);

    //create lepton/antilepton in respective top rest frames
    LV genlepton_rftbar = (*topGenObjects.GenLepton_);
    LV genantilepton_rft =  (*topGenObjects.GenAntiLepton_);
    this->boostLeptonsToRespectiveTopQuarkRestFrames(genlepton_rftbar, genantilepton_rft, (*topGenObjects.GenTop_), (*topGenObjects.GenAntiTop_));

    if (BHadronIndex != -1 && AntiBHadronIndex != -1 && BHadronIndex!=AntiBHadronIndex) {
        orderLV(LeadGenBJet, NLeadGenBJet, (*topGenObjects.allGenJets_).at(BHadronIndex), (*topGenObjects.allGenJets_).at(AntiBHadronIndex), LVpt);
    }

    if ( (*topGenObjects.GenLepton_).pt() > selections.leptonPtCut_ && (*topGenObjects.GenAntiLepton_).pt() > selections.leptonPtCut_ &&
        std::fabs( (*topGenObjects.GenLepton_).eta() ) < selections.leptonEtaCut_ && std::fabs ( (*topGenObjects.GenAntiLepton_).eta() ) < selections.leptonEtaCut_) {

        if ( BHadronIndex != -1 && AntiBHadronIndex != -1 && BHadronIndex!=AntiBHadronIndex ) {
            if ( (*topGenObjects.allGenJets_).at(BHadronIndex).pt() > selections.jetPtCut_ && std::fabs ( (*topGenObjects.allGenJets_).at(BHadronIndex).eta() ) < selections.jetEtaCut_ &&
                 (*topGenObjects.allGenJets_).at(AntiBHadronIndex).pt() > selections.jetPtCut_ && std::fabs ( (*topGenObjects.allGenJets_).at(AntiBHadronIndex).Eta() ) < selections.jetEtaCut_ )
            {

                h_VisGenAll->Fill((*topGenObjects.GenTop_).M(), trueLevelWeight);
                h_VisGenAll_noweight->Fill((*topGenObjects.GenTop_).M(), trueLevelWeightNoPileup);

                h_VisGenLLBarpT->Fill(( (*topGenObjects.GenLepton_) + (*topGenObjects.GenAntiLepton_) ).Pt(), trueLevelWeight );
                h_VisGenLLBarMass->Fill(( (*topGenObjects.GenLepton_) + (*topGenObjects.GenAntiLepton_) ).M(), trueLevelWeight );

                h_VisGenLeptonpT->Fill((*topGenObjects.GenLepton_).pt(), trueLevelWeight );
                h_VisGenAntiLeptonpT->Fill((*topGenObjects.GenAntiLepton_).Pt(), trueLevelWeight );

                h_VisGenLeptonEta->Fill((*topGenObjects.GenLepton_).Eta(), trueLevelWeight );
                h_VisGenAntiLeptonEta->Fill((*topGenObjects.GenAntiLepton_).Eta(), trueLevelWeight );

                h_VisGenBJetEta->Fill((*topGenObjects.allGenJets_).at(BHadronIndex).Eta(), trueLevelWeight );
                h_VisGenAntiBJetEta->Fill((*topGenObjects.allGenJets_).at(AntiBHadronIndex).Eta(), trueLevelWeight );
                h_VisGenBJetRapidity->Fill((*topGenObjects.allGenJets_).at(BHadronIndex).Rapidity(), trueLevelWeight );
                h_VisGenAntiBJetRapidity->Fill((*topGenObjects.allGenJets_).at(AntiBHadronIndex).Rapidity(), trueLevelWeight );
                h_VisGenBJetpT->Fill((*topGenObjects.allGenJets_).at(BHadronIndex).Pt(), trueLevelWeight );
                h_VisGenAntiBJetpT->Fill((*topGenObjects.allGenJets_).at(AntiBHadronIndex).Pt(), trueLevelWeight );
                h_VisGenMet->Fill((*topGenObjects.GenMet_).Pt(), trueLevelWeight);

                //for HT, count only >= selections.jetPtCut_
                std::vector<int> genJetIndices = initialiseIndices(*topGenObjects.allGenJets_);
                selectIndices(genJetIndices, *topGenObjects.allGenJets_, LVpt, selections.jetPtCut_);
                genHT = getJetHT(genJetIndices, *topGenObjects.allGenJets_);
                h_VisGenHT->Fill(genHT, trueLevelWeight);

                h_VisGenLLBarDPhi->Fill(std::fabs( DeltaPhi((*topGenObjects.GenLepton_), (*topGenObjects.GenAntiLepton_))), trueLevelWeight );
                h_VisGenLLBarDEta->Fill(std::abs((*topGenObjects.GenLepton_).Eta()) - std::abs((*topGenObjects.GenAntiLepton_).Eta()), trueLevelWeight );
                h_VisGenLLBarCosPhiInTopRestFrame->Fill(std::cos( Angle(genlepton_rftbar, genantilepton_rft)), trueLevelWeight );
                h_VisGenLeptonantiBjetMass->Fill(( (*topGenObjects.GenLepton_) + (*topGenObjects.allGenJets_).at(AntiBHadronIndex) ).M(), trueLevelWeight );
                h_VisGenAntiLeptonBjetMass->Fill(( (*topGenObjects.GenAntiLepton_) + (*topGenObjects.allGenJets_).at(BHadronIndex) ).M(), trueLevelWeight );
                h_VisGenJetMult->Fill((*topGenObjects.allGenJets_).size(), trueLevelWeight );

                h_VisGenBBBarpT->Fill(((*topGenObjects.allGenJets_).at(BHadronIndex) + (*topGenObjects.allGenJets_).at(AntiBHadronIndex)).Pt(), trueLevelWeight );
                h_VisGenBBBarMass->Fill(((*topGenObjects.allGenJets_).at(BHadronIndex) + (*topGenObjects.allGenJets_).at(AntiBHadronIndex)).M(), trueLevelWeight );
                h_VisGenBBBarDPhi->Fill(std::fabs( DeltaPhi((*topGenObjects.allGenJets_).at(BHadronIndex), (*topGenObjects.allGenJets_).at(AntiBHadronIndex))), trueLevelWeight );

                //Begin: Select & Fill histograms with Leading pT and 2nd Leading pT: Lepton and BJet
                h_VisGenLeptonpTLead->Fill(LeadGenLepton.Pt(), trueLevelWeight);
                h_VisGenLeptonpTNLead->Fill(NLeadGenLepton.Pt(), trueLevelWeight);
                h_VisGenLeptonEtaLead->Fill(LeadGenLepton.Eta(), trueLevelWeight);
                h_VisGenLeptonEtaNLead->Fill(NLeadGenLepton.Eta(), trueLevelWeight);

                h_VisGenBJetpTLead->Fill(LeadGenBJet.Pt(), trueLevelWeight);
                h_VisGenBJetpTNLead->Fill(NLeadGenBJet.Pt(), trueLevelWeight);
                h_VisGenBJetEtaLead->Fill(LeadGenBJet.Eta(), trueLevelWeight);
                h_VisGenBJetEtaNLead->Fill(NLeadGenBJet.Eta(), trueLevelWeight);
                //End: Select & Fill histograms with Leading pT and 2nd Leading pT: Lepton and BJet
            }
        }
    }

    LV genttbar((*topGenObjects.GenTop_) + (*topGenObjects.GenAntiTop_));
    h_VisGenTTBarMass->Fill(genttbar.M(), trueLevelWeight );
    h_VisGenTTBarRapidity->Fill(genttbar.Rapidity(), trueLevelWeight );
    h_VisGenTTBarRapidityAbs->Fill(std::abs(genttbar.Rapidity()), trueLevelWeight );
    h_VisGenTTBarpT->Fill(genttbar.Pt(), trueLevelWeight );

    h_VisGenTTBarRapidityvsTTBarpT->Fill(((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Pt(),((*topGenObjects.GenTop_)+(*topGenObjects.GenAntiTop_)).Rapidity(),trueLevelWeight);

    h_VisGenToppT->Fill((*topGenObjects.GenTop_).Pt(), trueLevelWeight );
    h_VisGenAntiToppT->Fill((*topGenObjects.GenAntiTop_).Pt(), trueLevelWeight );
    h_VisGenTopRapidity->Fill((*topGenObjects.GenTop_).Rapidity(), trueLevelWeight );
    h_VisGenAntiTopRapidity->Fill((*topGenObjects.GenAntiTop_).Rapidity(), trueLevelWeight );
    h_VisGenTopRapidityAbs->Fill(std::abs((*topGenObjects.GenTop_).Rapidity()), trueLevelWeight );
    h_VisGenAntiTopRapidityAbs->Fill(std::abs((*topGenObjects.GenAntiTop_).Rapidity()), trueLevelWeight );
    h_VisGenTopEta->Fill((*topGenObjects.GenTop_).Eta(), trueLevelWeight );
    h_VisGenAntiTopEta->Fill((*topGenObjects.GenAntiTop_).Eta(), trueLevelWeight );
    h_VisGenTTBarDeltaPhi->Fill(std::abs(DeltaPhi((*topGenObjects.GenTop_), (*topGenObjects.GenAntiTop_))), trueLevelWeight);
    h_VisGenTTBarDeltaRapidity->Fill(std::abs((*topGenObjects.GenTop_).Rapidity()) - std::abs((*topGenObjects.GenAntiTop_).Rapidity()), trueLevelWeight);

    h_VisGenNeutrinopT->Fill((*topGenObjects.GenNeutrino_).Pt(), trueLevelWeight);
    h_VisGenAntiNeutrinopT->Fill((*topGenObjects.GenAntiNeutrino_).Pt(), trueLevelWeight);

    /// Parton momentum fraction
    double partonMomFraction = ((*topGenObjects.GenTop_).energy() - (*topGenObjects.GenTop_).Pz() + (*topGenObjects.GenAntiTop_).energy() - (*topGenObjects.GenAntiTop_).Pz()) / (2 * 4000);
    double antipartonMomFraction = ((*topGenObjects.GenTop_).energy() + (*topGenObjects.GenTop_).Pz() + (*topGenObjects.GenAntiTop_).energy() + (*topGenObjects.GenAntiTop_).Pz()) / (2 * 4000);
    h_VisGenTopPartonFraction->Fill(partonMomFraction, trueLevelWeight);
    h_VisGenAntiTopPartonFraction->Fill(antipartonMomFraction, trueLevelWeight);

    //Begin: Fill histograms with Leading pT and 2nd Leading pT: Top
    orderLV(LeadGenTop, NLeadGenTop, (*topGenObjects.GenTop_), (*topGenObjects.GenAntiTop_), LVpt);
    h_VisGenToppTLead->Fill(LeadGenTop.Pt(), trueLevelWeight);
    h_VisGenToppTNLead->Fill(NLeadGenTop.Pt(), trueLevelWeight);
    h_VisGenTopRapidityLead->Fill(LeadGenTop.Rapidity(), trueLevelWeight);
    h_VisGenTopRapidityNLead->Fill(NLeadGenTop.Rapidity(), trueLevelWeight);
    h_VisGenTopMassLead->Fill(LeadGenTop.M(), trueLevelWeight);
    h_VisGenTopMassNLead->Fill(NLeadGenTop.M(), trueLevelWeight);
    //End: Fill histograms with Leading pT and 2nd Leading pT: Top

    // create top/antitop quark in the ttbar rest frame
    ROOT::Math::Boost CoMBoostGenTtbar (genttbar.BoostToCM());
    LV gentop ((*topGenObjects.GenTop_));
    LV genantitop ((*topGenObjects.GenAntiTop_));
    gentop = CoMBoostGenTtbar(gentop);
    genantitop = CoMBoostGenTtbar(genantitop);

    h_VisGenToppTTTRestFrame->Fill(gentop.Pt(), trueLevelWeight);
    h_VisGenAntiToppTTTRestFrame->Fill(genantitop.Pt(), trueLevelWeight);

}

void TopAnalysis::generatorTTbarjetsEvent(double& jetHTGen,
                                          const int bHadronIndex, const int antiBHadronIndex,
                                          const double trueLevelWeight,
                                          std::vector<int>& genJetIndices,std::vector<int>& genJet40Indices, std::vector<int>& genJet60Indices,
                                          std::vector<int>& genJet100Indices, std::vector<int>& genExtraJetIndices,
                                          const TopGenObjects& topGenObjects)
{
    // Use utilities without namespaces
    using namespace common;
    using ROOT::Math::VectorUtil::DeltaPhi;
    using ROOT::Math::VectorUtil::DeltaR;
    using ROOT::Math::VectorUtil::Angle;

    double cbin[20]={25.,35.,45.,55.,65.,75.,85.,95.,110.,130.,150.,170.,190.,210.,230.,250.,270.,300.,340.,380.};
    jetHTGen = 0.;

    if(!topGenObjects.valuesSet_) return;

    // Access object selections from config
    const AnalysisConfig::Selections& selections = analysisConfig_.selections();

    const int BHadronIndex(bHadronIndex);
    const int AntiBHadronIndex(antiBHadronIndex);

    if ( (*topGenObjects.GenLepton_).pt() > selections.leptonPtCut_ && (*topGenObjects.GenAntiLepton_).pt() > selections.leptonPtCut_ &&
        std::fabs( (*topGenObjects.GenLepton_).eta() ) < selections.leptonEtaCut_ && std::fabs ( (*topGenObjects.GenAntiLepton_).eta() ) < selections.leptonEtaCut_ ) {
        if ( BHadronIndex != -1 && AntiBHadronIndex != -1 && BHadronIndex!=AntiBHadronIndex) {
            if ( (*topGenObjects.allGenJets_).at(BHadronIndex).Pt() > selections.jetPtCut_ && std::fabs ( (*topGenObjects.allGenJets_).at(BHadronIndex).Eta() ) < selections.jetEtaCut_ &&
                 (*topGenObjects.allGenJets_).at(AntiBHadronIndex).Pt() > selections.jetPtCut_ && std::fabs ( (*topGenObjects.allGenJets_).at(AntiBHadronIndex).Eta() ) < selections.jetEtaCut_ &&
                 std::fabs(DeltaR((*topGenObjects.GenLepton_), (*topGenObjects.allGenJets_).at(BHadronIndex))) > selections.genDeltaRLeptonJetCut_  &&
                 std::fabs(DeltaR((*topGenObjects.GenLepton_), (*topGenObjects.allGenJets_).at(AntiBHadronIndex))) > selections.genDeltaRLeptonJetCut_  &&
                 std::fabs(DeltaR((*topGenObjects.GenAntiLepton_), (*topGenObjects.allGenJets_).at(BHadronIndex))) > selections.genDeltaRLeptonJetCut_  &&
                 std::fabs(DeltaR((*topGenObjects.GenAntiLepton_), (*topGenObjects.allGenJets_).at(AntiBHadronIndex))) > selections.genDeltaRLeptonJetCut_)
            {
                for(int genJet=0; genJet < (int) genExtraJetIndices.size();genJet++)
                {
                      jetHTGen+=(*topGenObjects.allGenJets_).at(genExtraJetIndices.at(genJet)).Pt();
                }
                h_VisGenJetMultpt30->Fill(genJetIndices.size(),trueLevelWeight);
                h_VisGenJetMultpt40->Fill(genJet40Indices.size(),trueLevelWeight);
                h_VisGenJetMultpt60->Fill(genJet60Indices.size(),trueLevelWeight);
                h_VisGenJetMultpt100->Fill(genJet100Indices.size(),trueLevelWeight);

                LV genttbar((*topGenObjects.GenTop_) + (*topGenObjects.GenAntiTop_));
                h_VisGenTTBar0Mass->Fill(rho0/genttbar.M(),trueLevelWeight);
                //std::cout << "AMK TA 4046 " << genExtraJetIndices.size() << std::endl;

                if(genExtraJetIndices.size()>3){
                    h_VisGenExtraJetpT4->Fill((*topGenObjects.allGenJets_).at(genExtraJetIndices.at(3)).Pt(),trueLevelWeight);
                    h_VisGenExtraJetEta4->Fill((*topGenObjects.allGenJets_).at(genExtraJetIndices.at(3)).Eta(),trueLevelWeight);
                    h_VisGenExtraJetAbsEta4->Fill(std::fabs((*topGenObjects.allGenJets_).at(genExtraJetIndices.at(3)).Eta()),trueLevelWeight);
                } if(genExtraJetIndices.size()>2){
                    h_VisGenExtraJetpT3->Fill((*topGenObjects.allGenJets_).at(genExtraJetIndices.at(2)).Pt(),trueLevelWeight);
                    h_VisGenExtraJetEta3->Fill((*topGenObjects.allGenJets_).at(genExtraJetIndices.at(2)).Eta(),trueLevelWeight);
                    h_VisGenExtraJetAbsEta3->Fill(std::fabs((*topGenObjects.allGenJets_).at(genExtraJetIndices.at(2)).Eta()),trueLevelWeight);
                } if(genExtraJetIndices.size()>1){
                    h_VisGenExtraJetpT2->Fill((*topGenObjects.allGenJets_).at(genExtraJetIndices.at(1)).Pt(),trueLevelWeight);
                    h_VisGenExtraJetEta2->Fill((*topGenObjects.allGenJets_).at(genExtraJetIndices.at(1)).Eta(),trueLevelWeight);
                    h_VisGenExtraJetAbsEta2->Fill(std::fabs((*topGenObjects.allGenJets_).at(genExtraJetIndices.at(1)).Eta()),trueLevelWeight);
                    h_VisGenMassExtraJet12->Fill(((*topGenObjects.allGenJets_).at(genExtraJetIndices.at(0))+(*topGenObjects.allGenJets_).at(genExtraJetIndices.at(1))).M(),trueLevelWeight);
                    h_VisGenDeltaRExtraJet12->Fill(std::fabs(DeltaR((*topGenObjects.allGenJets_).at(genExtraJetIndices.at(0)),(*topGenObjects.allGenJets_).at(genExtraJetIndices.at(1)))),trueLevelWeight);
                } if(genExtraJetIndices.size()>0){
                    h_VisGenExtraJetHT->Fill(jetHTGen,trueLevelWeight);
                    h_VisGenExtraJetpT->Fill((*topGenObjects.allGenJets_).at(genExtraJetIndices.at(0)).Pt(),trueLevelWeight);
                    h_VisGenExtraJetEta->Fill((*topGenObjects.allGenJets_).at(genExtraJetIndices.at(0)).Eta(),trueLevelWeight);
                    h_VisGenExtraJetAbsEta->Fill(std::fabs((*topGenObjects.allGenJets_).at(genExtraJetIndices.at(0)).Eta()),trueLevelWeight);
                    LV genttbar((*topGenObjects.GenTop_) + (*topGenObjects.GenAntiTop_));
                    h_VisGenTTBar1stJetMass->Fill(rho0/(genttbar+(*topGenObjects.allGenJets_).at(genExtraJetIndices.at(0))).M(),trueLevelWeight);
                }

                for(int q0 = 0; q0<20; q0++){
                    h_VisGenJetMultTotal->Fill(cbin[q0],trueLevelWeight);
                    if((genExtraJetIndices.size() > 0 && (*topGenObjects.allGenJets_).at(genExtraJetIndices.at(0)).Pt()<=cbin[q0]) || genExtraJetIndices.size()<1 ) h_VisGenJetMultQ0->Fill(cbin[q0],trueLevelWeight);
                    if((genExtraJetIndices.size() > 1 && (*topGenObjects.allGenJets_).at(genExtraJetIndices.at(1)).Pt()<=cbin[q0]) || genExtraJetIndices.size()<2 ) h_VisGenJetExtra2Q0->Fill(cbin[q0],trueLevelWeight);
                    if(jetHTGen<=cbin[q0]) h_VisGenJetMultQsum->Fill(cbin[q0],trueLevelWeight);
                }
            }
        }
    }
}

void TopAnalysis::pseudoTopEvent(LV& leadPseudoTop, LV& nLeadPseudoTop,
                                    LV& leadPseudoLepton, LV& nLeadPseudoLepton,
                                    LV& leadPseudoBJet, LV& nLeadPseudoBJet,
                                    const double trueLevelWeightNoPileup, const double trueLevelWeight,
                                    std::vector<int>& pseudoJetIndices, std::vector<int>& pseudoJet40Indices, std::vector<int>& pseudoJet60Indices, std::vector<int>& pseudoJet100Indices,
                                    const TopPseudoObjects& topPseudoObjects)
{
    // Use utilities without namespaces
    using namespace common;
    using ROOT::Math::VectorUtil::DeltaPhi;
    using ROOT::Math::VectorUtil::DeltaR;
    using ROOT::Math::VectorUtil::Angle;

    LV& LeadPseudoTop(leadPseudoTop);
    LV& NLeadPseudoTop(nLeadPseudoTop);
    LV& LeadPseudoLepton(leadPseudoLepton);
    LV& NLeadPseudoLepton(nLeadPseudoLepton);
    LV& LeadPseudoBJet(leadPseudoBJet);
    LV& NLeadPseudoBJet(nLeadPseudoBJet);

    if(!topPseudoObjects.valuesSet_) return;

    // Access object selections from config
    const AnalysisConfig::Selections& selections = analysisConfig_.selections();

    // Define leading and trailing top
    orderLV(LeadPseudoTop, NLeadPseudoTop, (*topPseudoObjects.PseudoTop_), (*topPseudoObjects.PseudoAntiTop_), LVpt);

    //Begin: Select & Fill histograms with Leading pT and 2nd Leading pT: Lepton and BJet
    orderLV(LeadPseudoLepton, NLeadPseudoLepton, (*topPseudoObjects.PseudoLepton_), (*topPseudoObjects.PseudoAntiLepton_), LVpt);
    orderLV(LeadPseudoBJet, NLeadPseudoBJet, (*topPseudoObjects.PseudoBJet_), (*topPseudoObjects.PseudoAntiBJet_), LVpt);

    LV pseudodilepton((*topPseudoObjects.PseudoLepton_) + (*topPseudoObjects.PseudoAntiLepton_));

    if ( (*topPseudoObjects.PseudoLepton_).pt() > selections.leptonPtCut_ && (*topPseudoObjects.PseudoAntiLepton_).pt() > selections.leptonPtCut_ &&
          std::fabs( (*topPseudoObjects.PseudoLepton_).eta() ) < selections.leptonEtaCut_ && std::fabs ( (*topPseudoObjects.PseudoAntiLepton_).eta() ) < selections.leptonEtaCut_ &&
          (*topPseudoObjects.PseudoBJet_).pt() > selections.jetPtCut_ && std::fabs ( (*topPseudoObjects.PseudoBJet_).eta() ) < selections.jetEtaCut_ &&
          (*topPseudoObjects.PseudoAntiBJet_).pt() > selections.jetPtCut_ && std::fabs ( (*topPseudoObjects.PseudoAntiBJet_).Eta() ) < selections.jetEtaCut_ &&
          std::fabs(DeltaR((*topPseudoObjects.PseudoLepton_), (*topPseudoObjects.PseudoBJet_))) > selections.genDeltaRLeptonJetCut_  &&
          std::fabs(DeltaR((*topPseudoObjects.PseudoLepton_), (*topPseudoObjects.PseudoAntiBJet_))) > selections.genDeltaRLeptonJetCut_  &&
          std::fabs(DeltaR((*topPseudoObjects.PseudoAntiLepton_), (*topPseudoObjects.PseudoBJet_))) > selections.genDeltaRLeptonJetCut_  &&
          std::fabs(DeltaR((*topPseudoObjects.PseudoAntiLepton_), (*topPseudoObjects.PseudoAntiBJet_))) > selections.genDeltaRLeptonJetCut_ &&
          pseudodilepton.M() > 20.0 )

       {
                h_VisPseudoAll->Fill((*topPseudoObjects.PseudoTop_).M(), trueLevelWeight);
                h_VisPseudoAll_noweight->Fill((*topPseudoObjects.PseudoTop_).M(), trueLevelWeightNoPileup);

                // Particle jet multiplicities
                h_VisPseudoJetMultpt30->Fill(pseudoJetIndices.size(), trueLevelWeight);
                h_VisPseudoJetMultpt40->Fill(pseudoJet40Indices.size(), trueLevelWeight);
                h_VisPseudoJetMultpt60->Fill(pseudoJet60Indices.size(), trueLevelWeight);
                h_VisPseudoJetMultpt100->Fill(pseudoJet100Indices.size(), trueLevelWeight);

                //create lepton/antilepton in respective top rest frames
                LV pseudolepton_rftbar = (*topPseudoObjects.PseudoLepton_);
                LV pseudoantilepton_rft =  (*topPseudoObjects.PseudoAntiLepton_);
                this->boostLeptonsToRespectiveTopQuarkRestFrames(pseudolepton_rftbar, pseudoantilepton_rft, (*topPseudoObjects.PseudoTop_), (*topPseudoObjects.PseudoAntiTop_));

                // Particle objects of the decay products
                h_VisPseudoLLBarpT->Fill(( (*topPseudoObjects.PseudoLepton_) + (*topPseudoObjects.PseudoAntiLepton_) ).Pt(), trueLevelWeight );
                h_VisPseudoLLBarMass->Fill(( (*topPseudoObjects.PseudoLepton_) + (*topPseudoObjects.PseudoAntiLepton_) ).M(), trueLevelWeight );
                h_VisPseudoLLBarDPhi->Fill(std::fabs( DeltaPhi((*topPseudoObjects.PseudoLepton_), (*topPseudoObjects.PseudoAntiLepton_))), trueLevelWeight );
                h_VisPseudoLLBarDEta->Fill(std::abs((*topPseudoObjects.PseudoLepton_).Eta()) - std::abs((*topPseudoObjects.PseudoAntiLepton_).Eta()), trueLevelWeight );
                h_VisPseudoLLBarCosPhiInTopRestFrame->Fill(std::cos( Angle(pseudolepton_rftbar, pseudoantilepton_rft)), trueLevelWeight );

                h_VisPseudoLeptonpT->Fill((*topPseudoObjects.PseudoLepton_).pt(), trueLevelWeight );
                h_VisPseudoAntiLeptonpT->Fill((*topPseudoObjects.PseudoAntiLepton_).Pt(), trueLevelWeight );

                h_VisPseudoLeptonEta->Fill((*topPseudoObjects.PseudoLepton_).Eta(), trueLevelWeight );
                h_VisPseudoAntiLeptonEta->Fill((*topPseudoObjects.PseudoAntiLepton_).Eta(), trueLevelWeight );

                h_VisPseudoBJetEta->Fill((*topPseudoObjects.PseudoBJet_).Eta(), trueLevelWeight );
                h_VisPseudoAntiBJetEta->Fill((*topPseudoObjects.PseudoAntiBJet_).Eta(), trueLevelWeight );
                h_VisPseudoBJetRapidity->Fill((*topPseudoObjects.PseudoBJet_).Rapidity(), trueLevelWeight );
                h_VisPseudoAntiBJetRapidity->Fill((*topPseudoObjects.PseudoAntiBJet_).Rapidity(), trueLevelWeight );
                h_VisPseudoBJetpT->Fill((*topPseudoObjects.PseudoBJet_).Pt(), trueLevelWeight );
                h_VisPseudoAntiBJetpT->Fill((*topPseudoObjects.PseudoAntiBJet_).Pt(), trueLevelWeight );

                h_VisPseudoBBBarpT->Fill(((*topPseudoObjects.PseudoBJet_) + (*topPseudoObjects.PseudoAntiBJet_)).Pt(), trueLevelWeight );
                h_VisPseudoBBBarMass->Fill(((*topPseudoObjects.PseudoBJet_) + (*topPseudoObjects.PseudoAntiBJet_)).M(), trueLevelWeight );
                h_VisPseudoBBBarDPhi->Fill(std::fabs( DeltaPhi((*topPseudoObjects.PseudoBJet_), (*topPseudoObjects.PseudoAntiBJet_))), trueLevelWeight );

                h_VisPseudoNuNuBarpT->Fill((*topPseudoObjects.PseudoNeutrino_).Pt() + (*topPseudoObjects.PseudoAntiNeutrino_).Pt(), trueLevelWeight );

                h_VisPseudoLeptonpTLead->Fill(LeadPseudoLepton.Pt(), trueLevelWeight);
                h_VisPseudoLeptonpTNLead->Fill(NLeadPseudoLepton.Pt(), trueLevelWeight);
                h_VisPseudoLeptonEtaLead->Fill(LeadPseudoLepton.Eta(), trueLevelWeight);
                h_VisPseudoLeptonEtaNLead->Fill(NLeadPseudoLepton.Eta(), trueLevelWeight);

                h_VisPseudoBJetpTLead->Fill(LeadPseudoBJet.Pt(), trueLevelWeight);
                h_VisPseudoBJetpTNLead->Fill(NLeadPseudoBJet.Pt(), trueLevelWeight);
                h_VisPseudoBJetEtaLead->Fill(LeadPseudoBJet.Eta(), trueLevelWeight);
                h_VisPseudoBJetEtaNLead->Fill(NLeadPseudoBJet.Eta(), trueLevelWeight);

                h_VisPseudoLeptonantiBjetMass->Fill(( (*topPseudoObjects.PseudoLepton_) + (*topPseudoObjects.PseudoAntiBJet_) ).M(), trueLevelWeight );
                h_VisPseudoAntiLeptonBjetMass->Fill(( (*topPseudoObjects.PseudoAntiLepton_) + (*topPseudoObjects.PseudoBJet_) ).M(), trueLevelWeight );

                // Pseudo top objects
                LV pseudottbar((*topPseudoObjects.PseudoTop_) + (*topPseudoObjects.PseudoAntiTop_));
                h_VisPseudoTTBarMass->Fill(pseudottbar.M(), trueLevelWeight );
                h_VisPseudoTTBarRapidity->Fill(pseudottbar.Rapidity(), trueLevelWeight );
                h_VisPseudoTTBarRapidityAbs->Fill(std::abs(pseudottbar.Rapidity()), trueLevelWeight );
                h_VisPseudoTTBarpT->Fill(pseudottbar.Pt(), trueLevelWeight );
                h_VisPseudoTTBarDeltaPhi->Fill(std::abs(DeltaPhi((*topPseudoObjects.PseudoTop_), (*topPseudoObjects.PseudoAntiTop_))), trueLevelWeight);
                h_VisPseudoTTBarDeltaRapidity->Fill(std::abs((*topPseudoObjects.PseudoTop_).Rapidity()) - std::abs((*topPseudoObjects.PseudoAntiTop_).Rapidity()), trueLevelWeight);

                h_VisPseudoToppT->Fill((*topPseudoObjects.PseudoTop_).Pt(), trueLevelWeight );
                h_VisPseudoAntiToppT->Fill((*topPseudoObjects.PseudoAntiTop_).Pt(), trueLevelWeight );
                h_VisPseudoTopRapidity->Fill((*topPseudoObjects.PseudoTop_).Rapidity(), trueLevelWeight );
                h_VisPseudoAntiTopRapidity->Fill((*topPseudoObjects.PseudoAntiTop_).Rapidity(), trueLevelWeight );
                h_VisPseudoTopRapidityAbs->Fill(std::abs((*topPseudoObjects.PseudoTop_).Rapidity()), trueLevelWeight );
                h_VisPseudoAntiTopRapidityAbs->Fill(std::abs((*topPseudoObjects.PseudoAntiTop_).Rapidity()), trueLevelWeight );

                h_VisPseudoToppTLead->Fill(LeadPseudoTop.Pt(), trueLevelWeight);
                h_VisPseudoToppTNLead->Fill(NLeadPseudoTop.Pt(), trueLevelWeight);
                h_VisPseudoTopRapidityLead->Fill(LeadPseudoTop.Rapidity(), trueLevelWeight);
                h_VisPseudoTopRapidityNLead->Fill(NLeadPseudoTop.Rapidity(), trueLevelWeight);

                // create top/antitop quark in the ttbar rest frame
                ROOT::Math::Boost CoMBoostPseudoTtbar (pseudottbar.BoostToCM());
                LV pseudotop ((*topPseudoObjects.PseudoTop_));
                LV pseudoantitop ((*topPseudoObjects.PseudoAntiTop_));
                pseudotop = CoMBoostPseudoTtbar(pseudotop);
                pseudoantitop = CoMBoostPseudoTtbar(pseudoantitop);

                h_VisPseudoToppTTTRestFrame->Fill(pseudotop.Pt(), trueLevelWeight);
                h_VisPseudoAntiToppTTTRestFrame->Fill(pseudoantitop.Pt(), trueLevelWeight);

       }  else if (analysisConfig_.sampleComposition().pseudoTopMode_ == 2 &&
        (*topPseudoObjects.PseudoLepton_).Pt()> 18. && std::fabs((*topPseudoObjects.PseudoLepton_).Eta()) < 2.5 &&
        (*topPseudoObjects.PseudoAntiLepton_).Pt()> 18. && std::fabs((*topPseudoObjects.PseudoAntiLepton_).Eta()) < 2.5 &&
        (*topPseudoObjects.PseudoBJet_).pt() > 27. && std::fabs ( (*topPseudoObjects.PseudoBJet_).eta() ) < 2.8 &&
        (*topPseudoObjects.PseudoAntiBJet_).pt() > 27. && std::fabs ( (*topPseudoObjects.PseudoAntiBJet_).Eta() ) < 2.8 &&
        std::fabs(DeltaR((*topPseudoObjects.PseudoLepton_), (*topPseudoObjects.PseudoBJet_))) > selections.genDeltaRLeptonJetCut_  &&
        std::fabs(DeltaR((*topPseudoObjects.PseudoLepton_), (*topPseudoObjects.PseudoAntiBJet_))) > selections.genDeltaRLeptonJetCut_  &&
        std::fabs(DeltaR((*topPseudoObjects.PseudoAntiLepton_), (*topPseudoObjects.PseudoBJet_))) > selections.genDeltaRLeptonJetCut_  &&
        std::fabs(DeltaR((*topPseudoObjects.PseudoAntiLepton_), (*topPseudoObjects.PseudoAntiBJet_))) > selections.genDeltaRLeptonJetCut_ &&
        pseudodilepton.M() > 18.0 ) {// fill underflow for reco objects not in vis. phase space

                // Particle jet multiplicities
                h_VisPseudoJetMultpt30->Fill(-1000.,trueLevelWeight);
                h_VisPseudoJetMultpt40->Fill(-1000.,trueLevelWeight);
                h_VisPseudoJetMultpt60->Fill(-1000.,trueLevelWeight);
                h_VisPseudoJetMultpt100->Fill(-1000.,trueLevelWeight);

                // Particle objects of the decay products
                h_VisPseudoLLBarpT->Fill(-1000., trueLevelWeight );
                h_VisPseudoLLBarMass->Fill(-1000., trueLevelWeight );
                h_VisPseudoLLBarDPhi->Fill(-1000., trueLevelWeight );
                h_VisPseudoLLBarDEta->Fill(-1000., trueLevelWeight );
                h_VisPseudoLLBarCosPhiInTopRestFrame->Fill(-1000., trueLevelWeight );

                h_VisPseudoLeptonpT->Fill(-1000., trueLevelWeight );
                h_VisPseudoAntiLeptonpT->Fill(-1000., trueLevelWeight );

                h_VisPseudoLeptonEta->Fill(-1000., trueLevelWeight );
                h_VisPseudoAntiLeptonEta->Fill(-1000., trueLevelWeight );

                h_VisPseudoBJetEta->Fill(-1000., trueLevelWeight );
                h_VisPseudoAntiBJetEta->Fill(-1000., trueLevelWeight );
                h_VisPseudoBJetRapidity->Fill(-1000., trueLevelWeight );
                h_VisPseudoAntiBJetRapidity->Fill(-1000., trueLevelWeight );
                h_VisPseudoBJetpT->Fill(-1000., trueLevelWeight );
                h_VisPseudoAntiBJetpT->Fill(-1000., trueLevelWeight );

                h_VisPseudoBBBarpT->Fill(-1000., trueLevelWeight );
                h_VisPseudoBBBarMass->Fill(-1000., trueLevelWeight );
                h_VisPseudoBBBarDPhi->Fill(-1000., trueLevelWeight );

                h_VisPseudoNuNuBarpT->Fill(-1000., trueLevelWeight );

                h_VisPseudoLeptonpTLead->Fill(-1000., trueLevelWeight);
                h_VisPseudoLeptonpTNLead->Fill(-1000., trueLevelWeight);
                h_VisPseudoLeptonEtaLead->Fill(-1000., trueLevelWeight);
                h_VisPseudoLeptonEtaNLead->Fill(-1000., trueLevelWeight);

                h_VisPseudoBJetpTLead->Fill(-1000., trueLevelWeight);
                h_VisPseudoBJetpTNLead->Fill(-1000., trueLevelWeight);
                h_VisPseudoBJetEtaLead->Fill(-1000., trueLevelWeight);
                h_VisPseudoBJetEtaNLead->Fill(-1000., trueLevelWeight);

                h_VisPseudoLeptonantiBjetMass->Fill(-1000., trueLevelWeight );
                h_VisPseudoAntiLeptonBjetMass->Fill(-1000., trueLevelWeight );

                // Pseudo top objects
                h_VisPseudoTTBarMass->Fill(-1000., trueLevelWeight );
                h_VisPseudoTTBarRapidity->Fill(-1000., trueLevelWeight );
                h_VisPseudoTTBarRapidityAbs->Fill(-1000., trueLevelWeight );
                h_VisPseudoTTBarpT->Fill(-1000., trueLevelWeight );
                h_VisPseudoTTBarDeltaPhi->Fill(-1000., trueLevelWeight);
                h_VisPseudoTTBarDeltaRapidity->Fill(-1000., trueLevelWeight);

                h_VisPseudoToppT->Fill(-1000., trueLevelWeight );
                h_VisPseudoAntiToppT->Fill(-1000., trueLevelWeight );
                h_VisPseudoTopRapidity->Fill(-1000., trueLevelWeight );
                h_VisPseudoAntiTopRapidity->Fill(-1000., trueLevelWeight );
                h_VisPseudoTopRapidityAbs->Fill(-1000., trueLevelWeight );
                h_VisPseudoAntiTopRapidityAbs->Fill(-1000., trueLevelWeight );

                h_VisPseudoToppTLead->Fill(-1000., trueLevelWeight);
                h_VisPseudoToppTNLead->Fill(-1000., trueLevelWeight);
                h_VisPseudoTopRapidityLead->Fill(-1000., trueLevelWeight);
                h_VisPseudoTopRapidityNLead->Fill(-1000., trueLevelWeight);

                h_VisPseudoToppTTTRestFrame->Fill(-1000., trueLevelWeight);
                h_VisPseudoAntiToppTTTRestFrame->Fill(-1000., trueLevelWeight);

        }
}

void TopAnalysis::generatorVisJets(const TopGenObjects& topGenObjects,std::vector<int>& genVisJetIndices)
{
    if(topGenObjects.valuesSet_)
    {
        for(int genJet=0; genJet<(int)(*topGenObjects.allGenJets_).size(); genJet++)
        {
            if(std::fabs((*topGenObjects.allGenJets_).at(genJet).Eta() ) < analysisConfig_.selections().jetEtaCut_ && (*topGenObjects.allGenJets_).at(genJet).Pt() > analysisConfig_.selections().jetPtCut_)
            {
                genVisJetIndices.push_back(genJet);
            }
        }
    }
}

void TopAnalysis::generatorExtraJets(const TopGenObjects& topGenObjects,std::vector<int>& genExtraJetIndices, const int bHadronIndex, const int antiBHadronIndex)
{
    if(topGenObjects.valuesSet_)
    {
        const int BHadronIndex(bHadronIndex);
        const int AntiBHadronIndex(antiBHadronIndex);

        for(int genJet=0; genJet<(int)(*topGenObjects.allGenJets_).size(); genJet++)
        {
            if(BHadronIndex > -1 && AntiBHadronIndex > -1 && (*topGenObjects.allGenJets_).at(BHadronIndex) != (*topGenObjects.allGenJets_).at(genJet)  && (*topGenObjects.allGenJets_).at(AntiBHadronIndex) != (*topGenObjects.allGenJets_).at(genJet))
           {
                genExtraJetIndices.push_back(genJet);
            }
        }
    }
}

void TopAnalysis::lepton2dIsolationIndices(std::vector<int>& v_index, const VLV& v_lvLep, const VLV& v_lvJet)const
{
    std::vector<int> result;

    for(const int index : v_index){
        bool reject(false);
        for(const LV& jet : v_lvJet){
            const TLorentzVector jetTLV = common::LVtoTLV(jet);
            const TLorentzVector lepTLV = common::LVtoTLV(v_lvLep.at(index));
            const double deltaR = jetTLV.DeltaR(lepTLV);
            const double ptRel = (lepTLV.Vect().Cross(jetTLV.Vect()).Mag())/(jetTLV.Vect().Mag());
            if(deltaR < 0.5 && ptRel < 15.) reject = true;
        }
        if(!reject) result.push_back(index);
    }

    v_index.clear();
    v_index = result;
}

void TopAnalysis::boostLeptonsToRespectiveTopQuarkRestFrames(LV& lv_lepton, LV& lv_antilepton, const LV& lv_top, const LV& lv_antitop)
{
    const TVector3 boostTop = -(common::LVtoTLV(lv_top)).BoostVector();
    const TVector3 boostAntiTop =  -(common::LVtoTLV(lv_antitop)).BoostVector();
    TLorentzVector tlv_lepton = common::LVtoTLV(lv_lepton);
    TLorentzVector tlv_antilepton = common::LVtoTLV(lv_antilepton);
    tlv_lepton.Boost(boostAntiTop);
    tlv_antilepton.Boost(boostTop);
    lv_lepton = common::TLVtoLV(tlv_lepton);
    lv_antilepton = common::TLVtoLV(tlv_antilepton);
}


void TopAnalysis::CreateBinnedControlPlots(TH1* h_differential, TH1* h_control, const bool fromHistoList)
{
    auto &pair = (*binnedControlPlots_)[h_differential->GetName()];
    if(fromHistoList){
        if(analysisConfig_.general().era_ == Era::run1_8tev){
            HistoListReader histoList("data/binnings/obsolete/8tev/HistoList");
            if(histoList.IsZombie()) { std::cout << "Need a HistoList to create binned control plots!\n"; exit(273); }
            pair.first = histoList.getPlotProperties(h_differential->GetName()).getClonedHistogram();
        } else if(analysisConfig_.general().era_ == Era::run2_13tev_50ns || analysisConfig_.general().era_ == Era::run2_13tev_25ns
                 || analysisConfig_.general().era_ == Era::run2_13tev_2015_25ns || analysisConfig_.general().era_ == Era::run2_13tev_2016_25ns
                 || analysisConfig_.general().era_ == Era::run2_13tev_2016 || analysisConfig_.general().era_ == Era::run2_13tev_2017
                 || analysisConfig_.general().era_ == Era::run2_13tev_2018
		  || analysisConfig_.general().era_ == Era::run2_UL_13tev_2016preVFP || analysisConfig_.general().era_ == Era::run2_UL_13tev_2016postVFP || analysisConfig_.general().era_ == Era::run2_UL_13tev_2017
                 || analysisConfig_.general().era_ == Era::run2_UL_13tev_2018) {
                    HistoListReader histoList("data/binnings/obsolete/13tev/HistoList_13TeV");
                    if(histoList.IsZombie()) { std::cout << "Need a HistoList_13TeV to create binned control plots!\n"; exit(273); }
                    pair.first = histoList.getPlotProperties(h_differential->GetName()).getClonedHistogram();

        }
        else {
            std::cerr<<"ERROR in TopAnalysis::CreateBinnedControlPlots() ! Era is not supported! "<<"\n...break\n"<<std::endl;
            exit(323);
        }
    }
    else{
        bool old = TH1::AddDirectoryStatus();
        TH1::AddDirectory(false);
        TH1* clone = static_cast<TH1*>(h_differential->Clone());
        TH1::AddDirectory(old);
        pair.first = clone;
    }
    std::string name = "bcp_";
    name.append(h_differential->GetName()).append("_bin_");
    //create maps if we are called for the first time with a certain h_differential
    if (pair.second.size() == 0) {
        for (int i = 0; i <= pair.first->GetNbinsX() + 1; ++i)
            pair.second.push_back(std::map<std::string, TH1*>());
    }
    //now really create the histograms
    for (int i = 0; i <= pair.first->GetNbinsX() + 1; ++i) {
        std::string binning =
            i == 0 ? "underflow" :
            i == pair.first->GetNbinsX() + 1 ? "overflow" :
            common::d2s(pair.first->GetBinLowEdge(i)) + " to " + common::d2s(pair.first->GetBinLowEdge(i+1));
        binning = std::string(" (") + h_differential->GetName() + " " + binning + ")";
        std::string n = name + std::to_string(i) + "_" + h_control->GetName();
        pair.second[i][h_control->GetName()] = store(
            new TH1D(n.c_str(), (std::string(h_control->GetName())+binning).c_str(),
                     h_control->GetNbinsX(), h_control->GetBinLowEdge(1),
                     h_control->GetBinLowEdge(h_control->GetNbinsX()+1)));
    }
}



void TopAnalysis::FillBinnedControlPlot(TH1* h_differential, double binvalue,
                                        TH1* h_control, double value, double weight)
{
    auto pair = (*binnedControlPlots_)[h_differential->GetName()];
    auto bin = pair.first->FindBin(binvalue);
    auto m = pair.second.at(bin);
    TH1* h = m[h_control->GetName()];
    if (!h) {
        std::cerr << "Error: please CreateBinnedControlPlots for " << h_differential->GetName()
            << " and " << h_control->GetName() << std::endl;
        exit(911);
    }
    h->Fill(value, weight);
}



void TopAnalysis::SetAllAnalyzers(std::vector<AnalyzerBaseClass*> v_analyzer)
{
    v_analyzer_ = v_analyzer;
}



void TopAnalysis::SetAllTreeHandlers(std::vector<TreeHandlerBase*> v_treeHandler)
{
    v_treeHandler_ = v_treeHandler;
}



void TopAnalysis::fillAll(const std::string& selectionStep,
                          const EventMetadata& eventMetadata,
                          const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                          const TopGenObjects& topGenObjects, const TopPseudoObjects& topPseudoObjects,
                          const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                          const LooseKinRecoSolution& looseKinRecoSolution,
                          const ttbar::GenObjectIndices& genObjectIndices, const ttbar::RecoObjectIndices& recoObjectIndices,
                          const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                          const double& defaultWeight,const double& defaultWeightLooseKinReco)const
{
    // In case b-tag efficiencies are produced, analysis output is not
    if(this->makeBtagEfficiencies()) return;

    for(AnalyzerBaseClass* analyzer : v_analyzer_){
        if(analyzer) analyzer->fill(eventMetadata,
                                    recoObjects, commonGenObjects,
				    topGenObjects,
				    kinematicReconstructionSolutions,
				    looseKinRecoSolution,
				    recoObjectIndices, genObjectIndices,
				    genLevelWeights, recoLevelWeights,
				    defaultWeight, defaultWeightLooseKinReco, selectionStep);
    }

    for(TreeHandlerBase* treeHandler : v_treeHandler_){
        if(treeHandler) treeHandler->fill(eventMetadata,
                                          recoObjects, commonGenObjects,
					  topGenObjects,
					  kinematicReconstructionSolutions,
					  recoObjectIndices, genObjectIndices,
					  genLevelWeights, recoLevelWeights,
					  defaultWeight, selectionStep);
    }

}



void TopAnalysis::bookAll()
{
    for(AnalyzerBaseClass* analyzer : v_analyzer_){
        if(analyzer) analyzer->book(fOutput);
    }
}



void TopAnalysis::clearAll()
{
    for(AnalyzerBaseClass* analyzer : v_analyzer_){
        if(analyzer) analyzer->clear();
    }
}
