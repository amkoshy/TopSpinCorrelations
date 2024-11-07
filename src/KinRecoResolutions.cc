#include <TH1F.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TChain.h>
#include <TCanvas.h>
#include <iostream>
#include "TLorentzVector.h"
#include "TMath.h"
#include <vector>
#include <string>
#include "../../common/include/classesFwd.h"
#include "../../common/include/classes.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/CommandLineParameters.h"

using namespace std;

float dphi(float genphi, float recophi){
  return 2*TMath::ASin(TMath::Sin((genphi-recophi)/2));
}

// Reco

float top_pt = 0;
TBranch* b_top_pt = 0;
float top_phi = 0;
TBranch* b_top_phi = 0;
float top_rapidity = 0;
TBranch* b_top_rapidity = 0;
float top_eta = 0;
TBranch* b_top_eta = 0;
float top_scatteringangle_ttbarframe = 0;
TBranch* b_top_scatteringangle_ttbarframe = 0;

float tbar_pt = 0;
TBranch* b_tbar_pt = 0;
float tbar_phi = 0;
TBranch* b_tbar_phi = 0;
float tbar_rapidity = 0;
TBranch* b_tbar_rapidity = 0;
float tbar_eta = 0;
TBranch* b_tbar_eta = 0;

float ttbar_pt = 0;
TBranch* b_ttbar_pt = 0;
float ttbar_phi = 0;
TBranch* b_ttbar_phi = 0;
float ttbar_rapidity = 0;
TBranch* b_ttbar_rapidity = 0;
float ttbar_eta = 0;
TBranch* b_ttbar_eta = 0;
float ttbar_mass = 0;
TBranch* b_ttbar_mass = 0;

// Gen

float gen_top_pt = 0;
TBranch* b_gen_top_pt = 0;
float gen_top_phi = 0;
TBranch* b_gen_top_phi = 0;
float gen_top_rapidity = 0;
TBranch* b_gen_top_rapidity = 0;
float gen_top_eta = 0;
TBranch* b_gen_top_eta = 0;
float gen_top_scatteringangle_ttbarframe = 0;
TBranch* b_gen_top_scatteringangle_ttbarframe = 0;

float gen_tbar_pt = 0;
TBranch* b_gen_tbar_pt = 0;
float gen_tbar_phi = 0;
TBranch* b_gen_tbar_phi = 0;
float gen_tbar_rapidity = 0;
TBranch* b_gen_tbar_rapidity = 0;
float gen_tbar_eta = 0;
TBranch* b_gen_tbar_eta = 0;

float gen_ttbar_pt = 0;
TBranch* b_gen_ttbar_pt = 0;
float gen_ttbar_phi = 0;
TBranch* b_gen_ttbar_phi = 0;
float gen_ttbar_rapidity = 0;
TBranch* b_gen_ttbar_rapidity = 0;
float gen_ttbar_eta = 0;
TBranch* b_gen_ttbar_eta = 0;
float gen_ttbar_mass = 0;
TBranch* b_gen_ttbar_mass = 0;


float weight = 0;
TBranch* b_weight = 0;


/// Main Function

void KinRecoResolutions(std::string era_, std::string prompt_or_viatau_, std::string filenumber_){

  // To run in root, example:
  // root -l -q src/KinRecoResolutions.cc++'("2017", "prompt", "0")'

  system("mkdir -p KinRecoResolutions");
  
  //TChain *step0chain = new TChain("ttBar_treeVariables_step0");
  TChain *step8chain = new TChain("ttBar_treeVariables_step8");
  TString filename;

  if (era_ == "2016preVFP" && prompt_or_viatau_ == "prompt") {
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016preVFP/Nominal/ee/ee_ttbarsignalplustau_fromDilepton_2016ULpreVFP_"+filenumber_+".root").c_str());
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016preVFP/Nominal/emu/emu_ttbarsignalplustau_fromDilepton_2016ULpreVFP_"+filenumber_+".root").c_str());
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016preVFP/Nominal/mumu/mumu_ttbarsignalplustau_fromDilepton_2016ULpreVFP_"+filenumber_+".root").c_str());

      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016preVFP/Nominal/ee/ee_ttbarsignalplustau_fromDilepton_2016ULpreVFP_"+filenumber_+".root").c_str());
      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016preVFP/Nominal/emu/emu_ttbarsignalplustau_fromDilepton_2016ULpreVFP_"+filenumber_+".root").c_str());
      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016preVFP/Nominal/mumu/mumu_ttbarsignalplustau_fromDilepton_2016ULpreVFP_"+filenumber_+".root").c_str());

      filename = "/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016preVFP/Nominal/mumu/mumu_ttbarsignalplustau_fromDilepton_2016ULpreVFP_"+filenumber_+".root";
  }
  else if (era_ == "2016preVFP" && prompt_or_viatau_ == "viatau") {
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016preVFP/Nominal/ee/ee_ttbarsignalviatau_fromDilepton_2016ULpreVFP_"+filenumber_+".root").c_str());
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016preVFP/Nominal/emu/emu_ttbarsignalviatau_fromDilepton_2016ULpreVFP_"+filenumber_+".root").c_str());
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016preVFP/Nominal/mumu/mumu_ttbarsignalviatau_fromDilepton_2016ULpreVFP_"+filenumber_+".root").c_str());

      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016preVFP/Nominal/ee/ee_ttbarsignalviatau_fromDilepton_2016ULpreVFP_"+filenumber_+".root").c_str());
      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016preVFP/Nominal/emu/emu_ttbarsignalviatau_fromDilepton_2016ULpreVFP_"+filenumber_+".root").c_str());
      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016preVFP/Nominal/mumu/mumu_ttbarsignalviatau_fromDilepton_2016ULpreVFP_"+filenumber_+".root").c_str());

      filename = "/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016preVFP/Nominal/mumu/mumu_ttbarsignalviatau_fromDilepton_2016ULpreVFP_"+filenumber_+".root";
  }
  else if (era_ == "2016postVFP" && prompt_or_viatau_ == "prompt") {
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016postVFP/Nominal/ee/ee_ttbarsignalplustau_fromDilepton_2016ULpostVFP_"+filenumber_+".root").c_str());
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016postVFP/Nominal/emu/emu_ttbarsignalplustau_fromDilepton_2016ULpostVFP_"+filenumber_+".root").c_str());
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016postVFP/Nominal/mumu/mumu_ttbarsignalplustau_fromDilepton_2016ULpostVFP_"+filenumber_+".root").c_str());

      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016postVFP/Nominal/ee/ee_ttbarsignalplustau_fromDilepton_2016ULpostVFP_"+filenumber_+".root").c_str());
      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016postVFP/Nominal/emu/emu_ttbarsignalplustau_fromDilepton_2016ULpostVFP_"+filenumber_+".root").c_str());
      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016postVFP/Nominal/mumu/mumu_ttbarsignalplustau_fromDilepton_2016ULpostVFP_"+filenumber_+".root").c_str());

      filename = "/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016postVFP/Nominal/mumu/mumu_ttbarsignalplustau_fromDilepton_2016ULpostVFP_"+filenumber_+".root";
  }
  else if (era_ == "2016postVFP" && prompt_or_viatau_ == "viatau") {
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016postVFP/Nominal/ee/ee_ttbarsignalviatau_fromDilepton_2016ULpostVFP_"+filenumber_+".root").c_str());
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016postVFP/Nominal/emu/emu_ttbarsignalviatau_fromDilepton_2016ULpostVFP_"+filenumber_+".root").c_str());
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016postVFP/Nominal/mumu/mumu_ttbarsignalviatau_fromDilepton_2016ULpostVFP_"+filenumber_+".root").c_str());

      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016postVFP/Nominal/ee/ee_ttbarsignalviatau_fromDilepton_2016ULpostVFP_"+filenumber_+".root").c_str());
      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016postVFP/Nominal/emu/emu_ttbarsignalviatau_fromDilepton_2016ULpostVFP_"+filenumber_+".root").c_str());
      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016postVFP/Nominal/mumu/mumu_ttbarsignalviatau_fromDilepton_2016ULpostVFP_"+filenumber_+".root").c_str());

      filename = "/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016postVFP/Nominal/mumu/mumu_ttbarsignalviatau_fromDilepton_2016ULpostVFP_"+filenumber_+".root";
  }
  else if (era_ == "2017" && prompt_or_viatau_ == "prompt") {
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2017/Nominal/ee/ee_ttbarsignalplustau_fromDilepton_2017UL_"+filenumber_+".root").c_str());
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2017/Nominal/emu/emu_ttbarsignalplustau_fromDilepton_2017UL_"+filenumber_+".root").c_str());
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2017/Nominal/mumu/mumu_ttbarsignalplustau_fromDilepton_2017UL_"+filenumber_+".root").c_str());

      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2017/Nominal/ee/ee_ttbarsignalplustau_fromDilepton_2017UL_"+filenumber_+".root").c_str());
      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2017/Nominal/emu/emu_ttbarsignalplustau_fromDilepton_2017UL_"+filenumber_+".root").c_str());
      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2017/Nominal/mumu/mumu_ttbarsignalplustau_fromDilepton_2017UL_"+filenumber_+".root").c_str());

      filename = "/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2017/Nominal/mumu/mumu_ttbarsignalplustau_fromDilepton_2017UL_"+filenumber_+".root";
  }
  else if (era_ == "2017" && prompt_or_viatau_ == "viatau") {
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2017/Nominal/ee/ee_ttbarsignalviatau_fromDilepton_2017UL_"+filenumber_+".root").c_str());
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2017/Nominal/emu/emu_ttbarsignalviatau_fromDilepton_2017UL_"+filenumber_+".root").c_str());
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2017/Nominal/mumu/mumu_ttbarsignalviatau_fromDilepton_2017UL_"+filenumber_+".root").c_str());

      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2017/Nominal/ee/ee_ttbarsignalviatau_fromDilepton_2017UL_"+filenumber_+".root").c_str());
      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2017/Nominal/emu/emu_ttbarsignalviatau_fromDilepton_2017UL_"+filenumber_+".root").c_str());
      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2017/Nominal/mumu/mumu_ttbarsignalviatau_fromDilepton_2017UL_"+filenumber_+".root").c_str());

      filename = "/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2017/Nominal/mumu/mumu_ttbarsignalviatau_fromDilepton_2017UL_"+filenumber_+".root";
  }
  else if (era_ == "2018" && prompt_or_viatau_ == "prompt") {
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2018/Nominal/ee/ee_ttbarsignalplustau_fromDilepton_2018UL_"+filenumber_+".root").c_str());
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2018/Nominal/emu/emu_ttbarsignalplustau_fromDilepton_2018UL_"+filenumber_+".root").c_str());
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2018/Nominal/mumu/mumu_ttbarsignalplustau_fromDilepton_2018UL_"+filenumber_+".root").c_str());

      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2018/Nominal/ee/ee_ttbarsignalplustau_fromDilepton_2018UL_"+filenumber_+".root").c_str());
      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2018/Nominal/emu/emu_ttbarsignalplustau_fromDilepton_2018UL_"+filenumber_+".root").c_str());
      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2018/Nominal/mumu/mumu_ttbarsignalplustau_fromDilepton_2018UL_"+filenumber_+".root").c_str());

      filename = "/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2018/Nominal/mumu/mumu_ttbarsignalplustau_fromDilepton_2018UL_"+filenumber_+".root";
  }
  else if (era_ == "2018" && prompt_or_viatau_ == "viatau") {
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2018/Nominal/ee/ee_ttbarsignalviatau_fromDilepton_2018UL_"+filenumber_+".root").c_str());
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2018/Nominal/emu/emu_ttbarsignalviatau_fromDilepton_2018UL_"+filenumber_+".root").c_str());
      //step0chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2018/Nominal/mumu/mumu_ttbarsignalviatau_fromDilepton_2018UL_"+filenumber_+".root").c_str());

      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2018/Nominal/ee/ee_ttbarsignalviatau_fromDilepton_2018UL_"+filenumber_+".root").c_str());
      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2018/Nominal/emu/emu_ttbarsignalviatau_fromDilepton_2018UL_"+filenumber_+".root").c_str());
      step8chain->AddFile(("/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2018/Nominal/mumu/mumu_ttbarsignalviatau_fromDilepton_2018UL_"+filenumber_+".root").c_str());

      filename = "/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2018/Nominal/mumu/mumu_ttbarsignalviatau_fromDilepton_2018UL_"+filenumber_+".root";
  }


  TFile* file = TFile::Open(filename);
  TH1D *weightedEvents = (TH1D*) file->Get("weightedEvents");


  step8chain->SetBranchAddress("top_pt", & top_pt, & b_top_pt);
  step8chain->SetBranchAddress("top_phi", & top_phi, & b_top_phi);
  step8chain->SetBranchAddress("top_rapidity", & top_rapidity, & b_top_rapidity);
  step8chain->SetBranchAddress("top_eta", & top_eta, & b_top_eta);
  step8chain->SetBranchAddress("top_scatteringangle_ttbarframe", & top_scatteringangle_ttbarframe, & b_top_scatteringangle_ttbarframe);

  step8chain->SetBranchAddress("tbar_pt", & tbar_pt, & b_tbar_pt);
  step8chain->SetBranchAddress("tbar_phi", & tbar_phi, & b_tbar_phi);
  step8chain->SetBranchAddress("tbar_rapidity", & tbar_rapidity, & b_tbar_rapidity);
  step8chain->SetBranchAddress("tbar_eta", & tbar_eta, & b_tbar_eta);

  step8chain->SetBranchAddress("ttbar_pt", & ttbar_pt, & b_ttbar_pt);
  step8chain->SetBranchAddress("ttbar_phi", & ttbar_phi, & b_ttbar_phi);
  step8chain->SetBranchAddress("ttbar_rapidity", & ttbar_rapidity, & b_ttbar_rapidity);
  step8chain->SetBranchAddress("ttbar_eta", & ttbar_eta, & b_ttbar_eta);
  step8chain->SetBranchAddress("ttbar_mass", & ttbar_mass, & b_ttbar_mass);

  step8chain->SetBranchAddress("gen_top_pt", & gen_top_pt, & b_gen_top_pt);
  step8chain->SetBranchAddress("gen_top_phi", & gen_top_phi, & b_gen_top_phi);
  step8chain->SetBranchAddress("gen_top_rapidity", & gen_top_rapidity, & b_gen_top_rapidity);
  step8chain->SetBranchAddress("gen_top_eta", & gen_top_eta, & b_gen_top_eta);
  step8chain->SetBranchAddress("gen_top_scatteringangle_ttbarframe", & gen_top_scatteringangle_ttbarframe, & b_gen_top_scatteringangle_ttbarframe);

  step8chain->SetBranchAddress("gen_tbar_pt", & gen_tbar_pt, & b_gen_tbar_pt);
  step8chain->SetBranchAddress("gen_tbar_phi", & gen_tbar_phi, & b_gen_tbar_phi);
  step8chain->SetBranchAddress("gen_tbar_rapidity", & gen_tbar_rapidity, & b_gen_tbar_rapidity);
  step8chain->SetBranchAddress("gen_tbar_eta", & gen_tbar_eta, & b_gen_tbar_eta);

  step8chain->SetBranchAddress("gen_ttbar_pt", & gen_ttbar_pt, & b_gen_ttbar_pt);
  step8chain->SetBranchAddress("gen_ttbar_phi", & gen_ttbar_phi, & b_gen_ttbar_phi);
  step8chain->SetBranchAddress("gen_ttbar_rapidity", & gen_ttbar_rapidity, & b_gen_ttbar_rapidity);
  step8chain->SetBranchAddress("gen_ttbar_eta", & gen_ttbar_eta, & b_gen_ttbar_eta);
  step8chain->SetBranchAddress("gen_ttbar_mass", & gen_ttbar_mass, & b_gen_ttbar_mass);

  step8chain->SetBranchAddress("eventWeight", & weight, & b_weight);


  TH2D *h_top_pT_genreco = new TH2D("top_pT_genreco", "p_{T}(t)^{gen} vs. p_{T}(t)^{rec} [GeV]", 40, 0, 800, 40, 0, 800);
  TH1D *h_top_pT_residual = new TH1D("top_pT_residual", "p_{T}(t)^{gen} - p_{T}(t)^{rec} [GeV]", 40, -400, 400);
  TH2D *h_top_pT_multiresidual = new TH2D("top_pT_multiresidual", "p_{T}(t)^{gen} vs. p_{T}(t)^{gen} - p_{T}(t)^{rec} [GeV]", 40, 0, 800, 40, -400, 400);

  TH2D *h_top_phi_genreco = new TH2D("top_phi_genreco", "#phi(t)^{gen} vs. #phi(t)^{rec}", 36, -3.1416, 3.1416, 36, -3.1416, 3.1416);
  TH1D *h_top_phi_residual = new TH1D("top_phi_residual", "#phi(t)^{gen} - #phi(t)^{rec}", 36, -1.578, 1.578);
  TH2D *h_top_phi_multiresidual = new TH2D("top_phi_multiresidual", "#phi(t)^{gen} vs. #phi(t)^{gen} - #phi(t)^{rec}", 36, -3.1416, 3.1416, 36, -1.578, 1.578);

  TH2D *h_top_rapidity_genreco = new TH2D("top_rapidity_genreco", "y(t)^{gen} vs. y(t)^{rec}", 40, -3.0, 3.0, 40, -3.0, 3.0);
  TH1D *h_top_rapidity_residual = new TH1D("top_rapidity_residual", "y(t)^{gen} - y(t)^{rec}", 40, -3.0, 3.0);
  TH2D *h_top_rapidity_multiresidual = new TH2D("top_rapidity_multiresidual", "y(t)^{gen} vs. y(t)^{gen} - y(t)^{rec}", 40, -3.0, 3.0, 40, -3.0, 3.0);

  TH2D *h_tbar_pT_genreco = new TH2D("tbar_pT_genreco", "p_{T}(#bar{t})^{gen} vs. p_{T}(#bar{t})^{rec} [GeV]", 40, 0, 800, 40, 0, 800);
  TH1D *h_tbar_pT_residual = new TH1D("tbar_pT_residual", "p_{T}(#bar{t})^{gen} - p_{T}(#bar{t})^{rec} [GeV]", 40, -400, 400);
  TH2D *h_tbar_pT_multiresidual = new TH2D("tbar_pT_multiresidual", "p_{T}(#bar{t})^{gen} vs. p_{T}(#bar{t})^{gen} - p_{T}(#bar{t})^{rec} [GeV]", 40, 0, 800, 40, -400, 400);

  TH2D *h_tbar_phi_genreco = new TH2D("tbar_phi_genreco", "#phi(#bar{t})^{gen} vs. #phi(#bar{t})^{rec}", 36, -3.1416, 3.1416, 36, -3.1416, 3.1416);
  TH1D *h_tbar_phi_residual = new TH1D("tbar_phi_residual", "#phi(#bar{t})^{gen} - #phi(#bar{t})^{rec}", 36, -1.578, 1.578);
  TH2D *h_tbar_phi_multiresidual = new TH2D("tbar_phi_multiresidual", "#phi(#bar{t})^{gen} vs. #phi(#bar{t})^{gen} - #phi(#bar{t})^{rec}", 36, 0, 3.1416, 36, -1.578, 1.578);

  TH2D *h_tbar_rapidity_genreco = new TH2D("tbar_rapidity_genreco", "y(#bar{t})^{gen} vs. y(#bar{t})^{rec}", 40, -3.0, 3.0, 40, -3.0, 3.0);
  TH1D *h_tbar_rapidity_residual = new TH1D("tbar_rapidity_residual", "y(#bar{t})^{gen} - y(#bar{t})^{rec}", 40, -3.0, 3.0);
  TH2D *h_tbar_rapidity_multiresidual = new TH2D("tbar_rapidity_multiresidual", "y(#bar{t})^{gen} vs. y(#bar{t})^{gen} - y(#bar{t})^{rec}", 40, -3.0, 3.0, 40, -3.0, 3.0);

  TH2D *h_top_scatteringangle_ttbarframe_genreco = new TH2D("top_scatteringangle_ttbarframe_genreco", "cos#Theta(t)^{gen} vs. cos#Theta(t)^{rec}", 40, -1.0, 1.0, 40, -1.0, 1.0);
  TH1D *h_top_scatteringangle_ttbarframe_residual = new TH1D("top_scatteringangle_ttbarframe_residual", "cos#Theta(t)^{gen} - cos#Theta(t)^{rec}", 40, -2.0, 2.0);
  TH2D *h_top_scatteringangle_ttbarframe_multiresidual = new TH2D("top_scatteringangle_ttbarframe_multiresidual", "cos#Theta(t)^{gen} vs. cos#Theta(t)^{gen} - cos#Theta(t)^{rec}", 40, -1.0, 1.0, 40, -2.0, 2.0);


  TH2D *h_ttbar_pT_genreco = new TH2D("ttbar_pT_genreco", "p_{T}(t#bar{t})^{gen} vs. p_{T}(t#bar{t})^{rec} [GeV]", 50, 0, 500, 50, 0, 500);
  TH1D *h_ttbar_pT_residual = new TH1D("ttbar_pT_residual", "p_{T}(t#bar{t})^{gen} - p_{T}(t#bar{t})^{rec} [GeV]", 50, -500, 500);
  TH2D *h_ttbar_pT_multiresidual = new TH2D("ttbar_pT_multiresidual", "p_{T}(t#bar{t})^{gen} vs. p_{T}(t#bar{t})^{gen} - p_{T}(t#bar{t})^{rec} [GeV]", 50, 0, 500, 50, -250, 250);

  TH2D *h_ttbar_phi_genreco = new TH2D("ttbar_phi_genreco", "#phi(t#bar{t})^{gen} vs. #phi(t#bar{t})^{rec}", 36, -3.1416, 3.1416, 36, -3.1416, 3.1416);
  TH1D *h_ttbar_phi_residual = new TH1D("ttbar_phi_residual", "#phi(t#bar{t})^{gen} - #phi(t#bar{t})^{rec}", 36, -1.578, 1.578);
  TH2D *h_ttbar_phi_multiresidual = new TH2D("ttbar_phi_multiresidual", "#phi(t#bar{t})^{gen} vs. #phi(t#bar{t})^{gen} - #phi(t#bar{t})^{rec}", 36, -3.1416, 3.1416, 36, -1.578, 1.578);

  TH2D *h_ttbar_rapidity_genreco = new TH2D("ttbar_rapidity_genreco", "y(t#bar{t})^{gen} vs. y(t#bar{t})^{rec}", 40, -3.0, 3.0, 40, -3.0, 3.0);
  TH1D *h_ttbar_rapidity_residual = new TH1D("ttbar_rapidity_residual", "y(t#bar{t})^{gen} - y(t#bar{t})^{rec}", 40, -3.0, 3.0);
  TH2D *h_ttbar_rapidity_multiresidual = new TH2D("ttbar_rapidity_multiresidual", "y(t#bar{t})^{gen} vs. y(t#bar{t})^{gen} - y(t#bar{t})^{rec}", 40, -3.0, 3.0, 40, -3.0, 3.0);

  TH2D *h_ttbar_mass_genreco = new TH2D("ttbar_mass_genreco","M(t#bar{t})^{gen} vs. M(t#bar{t})^{rec} [GeV]", 46, 200.0, 2000.0, 46, 200.0, 2000.0);
  TH1D *h_ttbar_mass_residual = new TH1D("ttbar_mass_residual", "M(t#bar{t})^{gen} - M(t#bar{t})^{rec} [GeV]", 50, -1000.0, 1000.0);
  TH2D *h_ttbar_mass_multiresidual = new TH2D("ttbar_mass_multiresidual","M(t#bar{t})^gen} vs. M(t#bar{t})^{gen} - M(t#bar{t})^{rec} [GeV]", 46, 200.0, 2000.0, 50, -1000.0, 1000.0);


  Long64_t nEvents = step8chain->GetEntries();
  //  nEvents = 10000;

  //Event loop
  for (Long64_t i = 0; i < nEvents; i++) {

    if (i % (nEvents/10) == 0) cout << "=== Event " << i/(nEvents/10) * 10 << "%" << endl;
    //  cout << "=== Event " << i  << endl;

    step8chain->GetEntry(i);

    h_top_pT_genreco->Fill(gen_top_pt, top_pt, weight);
    h_top_pT_residual->Fill(gen_top_pt - top_pt, weight);
    h_top_pT_multiresidual->Fill(gen_top_pt, gen_top_pt - top_pt, weight);

    h_top_phi_genreco->Fill(gen_top_phi, top_phi, weight);
    h_top_phi_residual->Fill(gen_top_phi - top_phi, weight);
    h_top_phi_multiresidual->Fill(gen_top_phi, gen_top_phi - top_phi, weight);

    h_top_rapidity_genreco->Fill(gen_top_rapidity, top_rapidity, weight);
    h_top_rapidity_residual->Fill(gen_top_rapidity - top_rapidity, weight);
    h_top_rapidity_multiresidual->Fill(gen_top_rapidity, gen_top_rapidity - top_rapidity, weight);

    h_tbar_pT_genreco->Fill(gen_tbar_pt, tbar_pt, weight);
    h_tbar_pT_residual->Fill(gen_tbar_pt - tbar_pt, weight);
    h_tbar_pT_multiresidual->Fill(gen_tbar_pt, gen_tbar_pt - tbar_pt, weight);

    h_tbar_phi_genreco->Fill(gen_tbar_phi, tbar_phi, weight);
    h_tbar_phi_residual->Fill(gen_tbar_phi - tbar_phi, weight);
    h_tbar_phi_multiresidual->Fill(gen_tbar_phi, gen_tbar_phi - tbar_phi, weight);

    h_tbar_rapidity_genreco->Fill(gen_tbar_rapidity, tbar_rapidity, weight);
    h_tbar_rapidity_residual->Fill(gen_tbar_rapidity - tbar_rapidity, weight);
    h_tbar_rapidity_multiresidual->Fill(gen_tbar_rapidity, gen_tbar_rapidity - tbar_rapidity, weight);

    h_top_scatteringangle_ttbarframe_genreco->Fill(gen_top_scatteringangle_ttbarframe, top_scatteringangle_ttbarframe, weight);
    h_top_scatteringangle_ttbarframe_residual->Fill(gen_top_scatteringangle_ttbarframe - top_scatteringangle_ttbarframe, weight);
    h_top_scatteringangle_ttbarframe_multiresidual->Fill(gen_top_scatteringangle_ttbarframe, gen_top_scatteringangle_ttbarframe - top_scatteringangle_ttbarframe, weight);

    h_ttbar_pT_genreco->Fill(gen_ttbar_pt, ttbar_pt, weight);
    h_ttbar_pT_residual->Fill(gen_ttbar_pt - ttbar_pt, weight);
    h_ttbar_pT_multiresidual->Fill(gen_ttbar_pt, gen_ttbar_pt - ttbar_pt, weight);

    h_ttbar_phi_genreco->Fill(gen_ttbar_phi, ttbar_phi, weight);
    h_ttbar_phi_residual->Fill(gen_ttbar_phi - ttbar_phi, weight);
    h_ttbar_phi_multiresidual->Fill(gen_ttbar_phi, gen_ttbar_phi - ttbar_phi, weight);

    h_ttbar_rapidity_genreco->Fill(gen_ttbar_rapidity, ttbar_rapidity, weight);
    h_ttbar_rapidity_residual->Fill(gen_ttbar_rapidity - ttbar_rapidity, weight);
    h_ttbar_rapidity_multiresidual->Fill(gen_ttbar_rapidity, gen_ttbar_rapidity - ttbar_rapidity, weight);

    h_ttbar_mass_genreco->Fill(gen_ttbar_mass, ttbar_mass, weight);
    h_ttbar_mass_residual->Fill(gen_ttbar_mass - ttbar_mass, weight);
    h_ttbar_mass_multiresidual->Fill(gen_ttbar_mass, gen_ttbar_mass - ttbar_mass, weight);

  }


  TFile *output_hists = new TFile(("KinRecoResolutions/output_hists_"+prompt_or_viatau_+"_"+era_+"_"+filenumber_+".root").c_str(),"RECREATE");

  h_top_pT_genreco->Write();
  h_top_pT_residual->Write();
  h_top_pT_multiresidual->Write();

  h_top_phi_genreco->Write();
  h_top_phi_residual->Write();
  h_top_phi_multiresidual->Write();

  h_top_rapidity_genreco->Write();
  h_top_rapidity_residual->Write();
  h_top_rapidity_multiresidual->Write();

  h_tbar_pT_genreco->Write();
  h_tbar_pT_residual->Write();
  h_tbar_pT_multiresidual->Write();

  h_tbar_phi_genreco->Write();
  h_tbar_phi_residual->Write();
  h_tbar_phi_multiresidual->Write();

  h_tbar_rapidity_genreco->Write();
  h_tbar_rapidity_residual->Write();
  h_tbar_rapidity_multiresidual->Write();

  h_top_scatteringangle_ttbarframe_genreco->Write();
  h_top_scatteringangle_ttbarframe_residual->Write();
  h_top_scatteringangle_ttbarframe_multiresidual->Write();

  h_ttbar_pT_genreco->Write();
  h_ttbar_pT_residual->Write();
  h_ttbar_pT_multiresidual->Write();

  h_ttbar_phi_genreco->Write();
  h_ttbar_phi_residual->Write();
  h_ttbar_phi_multiresidual->Write();

  h_ttbar_rapidity_genreco->Write();
  h_ttbar_rapidity_residual->Write();
  h_ttbar_rapidity_multiresidual->Write();

  h_ttbar_mass_genreco->Write();
  h_ttbar_mass_residual->Write();
  h_ttbar_mass_multiresidual->Write();


  weightedEvents->Write();


  output_hists->Close();

}


int main(int argc, char** argv) {

  CLParameter<std::string> opt_y("y", "Era", false, 1, 1);
  CLParameter<std::string> opt_t("t", "Prompt or via Tau", false, 1, 1);
  CLParameter<std::string> opt_n("n", "File Number", false, 1, 1);
  CLAnalyser::interpretGlobal(argc, argv);


  KinRecoResolutions(opt_y[0], opt_t[0], opt_n[0]);

}
