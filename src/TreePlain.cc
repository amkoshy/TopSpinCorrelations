#include "TreePlain.h"
#include "analysisStructs.h"
#include "analysisObjectStructs.h"
#include "KinematicReconstructionSolution.h"
#include "analysisUtils.h"
#include "AnalysisBase.h"
#include "utils.h"
#include <TChain.h>
#include <TFile.h>
#include "Math/VectorUtil.h"
#include <memory>
#include <map>
#include "boost/variant.hpp"
//#include "cassert"


// TreePlain::TreePlain(const TString& name, const TString& title, const bool flagGen):
TreePlain::TreePlain(const TreePlainSetup& setup_):
  // TODO move tree creation to InitWriteTree()
  //TTree(name.Data(), title.Data()),
  isGen_(setup_.isGen),
  fillRecoTree_(setup_.fillRecoTree),
  fillSimpleRecoTree_(setup_.fillSimpleRecoTree),
  fillGenTree_(setup_.fillGenTree),
  isTTbarSample_(setup_.isTTbarSample),
  isSystSample_(setup_.isSystSample),
  isSystVariation_(setup_.isSystVariation),
  fillBasicVariables_(false),
  fillMeasureVariables_(false),
  fillLooseKinRecoVariables_(false),
  fillPUVarWeights_(false),
  fillJetVariables_(false),
  fillPseudoVariables_(false),
  fillGenVariationWeights_(false),
  fillRecoVariationWeights_(false)
{
  TTreePtr = new TTree(setup_.name, setup_.title);
  fillBasicVariables_         = (fillRecoTree_ || fillGenTree_);
  fillParticleLevelVariables_ = fillRecoTree_;
  fillMeasureVariables_       = (fillGenTree_) && (isGen_) && (isTTbarSample_);
  fillGenVariationWeights_    = (isTTbarSample_ && !isSystSample_ && !isSystVariation_) && (isGen_);
  fillPUVarWeights_           = (!isSystSample_ && !isSystVariation_) && (isGen_);
  fillRecoVariationWeights_   = (!isSystSample_ && !isSystVariation_) && (fillRecoTree_ || fillSimpleRecoTree_) && (isGen_);
  fillLooseKinRecoVariables_  = (fillSimpleRecoTree_);
  fillJetVariables_           = true;
  fillPseudoVariables_        = isTTbarSample_;
}



// FIXME PORT=====================================================
TreePlain::TreePlain(const TString& name, const bool flagGen, const TString var, const bool DoParticle) :
  isGen_(flagGen)
{
  TChainPtr = new TChain(name);
  TTreePtr = TChainPtr;
  weightVarName = var;
  addFiducialInfoBranch = DoParticle;

}

TreePlain::~TreePlain()
{

  initVars.clear();
  //printf("~TreePlain()\n");
  //delete TTreePtr;
  if(TChainPtr)
    delete TChainPtr;
}

// TODO use argument?
void TreePlain::InitReading()
{
  //int ret = TChainPtr->LoadTree(0);
  //printf("ret = %d\n", ret);
  TChainPtr->LoadTree(0);
  //TChainPtr->SetMakeClass(1);
  TChainPtr->SetBranchStatus("*", 0);

  TBranch** b_pptr = NULL;

  //TChainPtr->SetBranchAddress("eventCounter", &Vars.eventCounter);
  SetBranchAddressVerbosely(TChainPtr, "eventCounter", Vars.eventCounter, &Vars.b_eventCounter);
  SetBranchAddressVerbosely(TChainPtr, "weight", Vars.weight, b_pptr);
  if(addFiducialInfoBranch)
    SetBranchAddressVerbosely(TChainPtr, "inFiducialPhaseSpace", Vars.inFiducialPhaseSpace, &Vars.b_inFiducialPhaseSpace);

  if(weightVarName == "" || weightVarName == "Nominal")
    return;


  if(weightVarName == "BTAG_UP")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_BTAG_U", Vars.var_BTAG_U, b_pptr);
    weightVar = &Vars.var_BTAG_U;
  }
  else if(weightVarName == "BTAG_DOWN")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_BTAG_D", Vars.var_BTAG_D, b_pptr);
    weightVar = &Vars.var_BTAG_D;
  }
  else if(weightVarName == "BTAG_LJET_UP")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_BTAG_LJET_U", Vars.var_BTAG_LJET_U, b_pptr);
    weightVar = &Vars.var_BTAG_LJET_U;
  }
  else if(weightVarName == "BTAG_LJET_DOWN")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_BTAG_LJET_D", Vars.var_BTAG_LJET_D, b_pptr);
    weightVar = &Vars.var_BTAG_LJET_D;
  }
  // else if(weightVarName == "BTAG_PT_UP")
  // {
  //   SetBranchAddressVerbosely(TChainPtr, "var_BTAG_PT_U", Vars.var_BTAG_PT_U, b_pptr);
  //   weightVar = &Vars.var_BTAG_PT_U;
  // }
  // else if(weightVarName == "BTAG_PT_DOWN")
  // {
  //   SetBranchAddressVerbosely(TChainPtr, "var_BTAG_PT_D", Vars.var_BTAG_PT_D, b_pptr);
  //   weightVar = &Vars.var_BTAG_PT_D;
  // }
  // else if(weightVarName == "BTAG_ETA_UP")
  // {
  //   SetBranchAddressVerbosely(TChainPtr, "var_BTAG_ETA_U", Vars.var_BTAG_ETA_U, b_pptr);
  //   weightVar = &Vars.var_BTAG_ETA_U;
  // }
  // else if(weightVarName == "BTAG_ETA_DOWN")
  // {
  //   SetBranchAddressVerbosely(TChainPtr, "var_BTAG_ETA_D", Vars.var_BTAG_ETA_D, b_pptr);
  //   weightVar = &Vars.var_BTAG_ETA_D;
  // }
  // else if(weightVarName == "BTAG_LJET_PT_UP")
  // {
  //   SetBranchAddressVerbosely(TChainPtr, "var_BTAG_LJET_PT_U", Vars.var_BTAG_LJET_PT_U, b_pptr);
  //   weightVar = &Vars.var_BTAG_LJET_PT_U;
  // }
  // else if(weightVarName == "BTAG_LJET_PT_DOWN")
  // {
  //   SetBranchAddressVerbosely(TChainPtr, "var_BTAG_LJET_PT_D", Vars.var_BTAG_LJET_PT_D, b_pptr);
  //   weightVar = &Vars.var_BTAG_LJET_PT_D;
  // }
  // else if(weightVarName == "BTAG_LJET_ETA_UP")
  // {
  //   SetBranchAddressVerbosely(TChainPtr, "var_BTAG_LJET_ETA_U", Vars.var_BTAG_LJET_ETA_U, b_pptr);
  //   weightVar = &Vars.var_BTAG_LJET_ETA_U;
  // }
  // else if(weightVarName == "BTAG_LJET_ETA_DOWN")
  // {
  //   SetBranchAddressVerbosely(TChainPtr, "var_BTAG_LJET_ETA_D", Vars.var_BTAG_LJET_ETA_D, b_pptr);
  //   weightVar = &Vars.var_BTAG_LJET_ETA_D;
  // }

  else if(weightVarName == "KIN_UP")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_KIN_U", Vars.var_KIN_U, b_pptr);
    weightVar = &Vars.var_KIN_U;
  }
  else if(weightVarName == "KIN_DOWN")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_KIN_D", Vars.var_KIN_D, b_pptr);
    weightVar = &Vars.var_KIN_D;
  }
  // else if(weightVarName == "LEPT_UP")
  // {
  //   SetBranchAddressVerbosely(TChainPtr, "var_LEPT_U", Vars.var_LEPT_U, b_pptr);
  //   weightVar = &Vars.var_LEPT_U;
  // }
  // else if(weightVarName == "LEPT_DOWN")
  // {
  //   SetBranchAddressVerbosely(TChainPtr, "var_LEPT_D", Vars.var_LEPT_D, b_pptr);
  //   weightVar = &Vars.var_LEPT_D;
  // }
  else if(weightVarName == "ELE_ID_UP"){
      SetBranchAddressVerbosely(TChainPtr, "var_ELE_ID_U", Vars.var_ELE_ID_U, b_pptr);
      weightVar = &Vars.var_ELE_ID_U;
  }
  else if(weightVarName == "ELE_ID_DOWN"){
      SetBranchAddressVerbosely(TChainPtr, "var_ELE_ID_D", Vars.var_ELE_ID_D, b_pptr);
      weightVar = &Vars.var_ELE_ID_D;
  }
  else if(weightVarName == "ELE_RECO_UP"){
        SetBranchAddressVerbosely(TChainPtr, "var_ELE_RECO_U", Vars.var_ELE_RECO_U, b_pptr);
        weightVar = &Vars.var_ELE_RECO_U;
  }
  else if(weightVarName == "ELE_RECO_DOWN"){
        SetBranchAddressVerbosely(TChainPtr, "var_ELE_RECO_D", Vars.var_ELE_RECO_D, b_pptr);
        weightVar = &Vars.var_ELE_RECO_D;
  }
  else if(weightVarName == "MUON_ID_UP"){
          SetBranchAddressVerbosely(TChainPtr, "var_MUON_ID_U", Vars.var_MUON_ID_U, b_pptr);
          weightVar = &Vars.var_MUON_ID_U;
  }
  else if(weightVarName == "MUON_ID_DOWN"){
          SetBranchAddressVerbosely(TChainPtr, "var_MUON_ID_D", Vars.var_MUON_ID_D, b_pptr);
          weightVar = &Vars.var_MUON_ID_D;
  }
  else if(weightVarName == "MUON_ISO_UP"){
            SetBranchAddressVerbosely(TChainPtr, "var_MUON_ISO_U", Vars.var_MUON_ISO_U, b_pptr);
            weightVar = &Vars.var_MUON_ISO_U;
  }
  else if(weightVarName == "MUON_ISO_DOWN"){
            SetBranchAddressVerbosely(TChainPtr, "var_MUON_ISO_D", Vars.var_MUON_ISO_D, b_pptr);
            weightVar = &Vars.var_MUON_ISO_D;
  }
  else if(weightVarName == "PU_UP")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_PU_U", Vars.var_PU_U, b_pptr);
    weightVar = &Vars.var_PU_U;
  }
  else if(weightVarName == "PU_DOWN")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_PU_D", Vars.var_PU_D, b_pptr);
    weightVar = &Vars.var_PU_D;
  }
  else if(weightVarName == "TRIG_UP")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_TRIG_U", Vars.var_TRIG_U, b_pptr);
    weightVar = &Vars.var_TRIG_U;
  }
  else if(weightVarName == "TRIG_DOWN")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_TRIG_D", Vars.var_TRIG_D, b_pptr);
    weightVar = &Vars.var_TRIG_D;
  }
  // else if(weightVarName == "TRIG_ETA_UP")
  // {
  //   SetBranchAddressVerbosely(TChainPtr, "var_TRIG_ETA_U", Vars.var_TRIG_ETA_U, b_pptr);
  //   weightVar = &Vars.var_TRIG_ETA_U;
  // }
  // else if(weightVarName == "TRIG_ETA_DOWN")
  // {
  //   SetBranchAddressVerbosely(TChainPtr, "var_TRIG_ETA_D", Vars.var_TRIG_ETA_D, b_pptr);
  //   weightVar = &Vars.var_TRIG_ETA_D;
  // }
  else if(weightVarName == "MEFACSCALE_UP")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_MEFACSCALE_U", Vars.var_MEFACSCALE_U, b_pptr);
    weightVar = &Vars.var_MEFACSCALE_U;
  }
  else if(weightVarName == "MEFACSCALE_DOWN")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_MEFACSCALE_D", Vars.var_MEFACSCALE_D, b_pptr);
    weightVar = &Vars.var_MEFACSCALE_D;
  }
  else if(weightVarName == "MERENSCALE_UP")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_MERENSCALE_U", Vars.var_MERENSCALE_U, b_pptr);
    weightVar = &Vars.var_MERENSCALE_U;
  }
  else if(weightVarName == "MERENSCALE_DOWN")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_MERENSCALE_D", Vars.var_MERENSCALE_D, b_pptr);
    weightVar = &Vars.var_MERENSCALE_D;
  }
  else if(weightVarName == "MESCALE_UP")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_MESCALE_U", Vars.var_MESCALE_U, b_pptr);
    weightVar = &Vars.var_MESCALE_U;
  }
  else if(weightVarName == "MESCALE_DOWN")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_MESCALE_D", Vars.var_MESCALE_D, b_pptr);
    weightVar = &Vars.var_MESCALE_D;
  }
  else if(weightVarName == "BFRAG_UP")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_BFRAG_U", Vars.var_BFRAG_U, b_pptr);
    weightVar = &Vars.var_BFRAG_U;
  }
  else if(weightVarName == "BFRAG_DOWN")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_BFRAG_D", Vars.var_BFRAG_D, b_pptr);
    weightVar = &Vars.var_BFRAG_D;
  }
  else if(weightVarName == "BFRAG_PETERSON")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_BFRAG_P", Vars.var_BFRAG_P, b_pptr);
    weightVar = &Vars.var_BFRAG_P;
  }
  // else if(weightVarName == "BFRAG_CENTRAL")
  // {
  //   SetBranchAddressVerbosely(TChainPtr, "var_BFRAG_C", Vars.var_BFRAG_C, b_pptr);
  //   weightVar = &Vars.var_BFRAG_C;
  // }
  else if(weightVarName == "BSEMILEP_UP")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_BSEMILEP_U", Vars.var_BSEMILEP_U, b_pptr);
    weightVar = &Vars.var_BSEMILEP_U;
  }
  else if(weightVarName == "BSEMILEP_DOWN")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_BSEMILEP_D", Vars.var_BSEMILEP_D, b_pptr);
    weightVar = &Vars.var_BSEMILEP_D;
  }
  else if(weightVarName == "PDF_ALPHAS_UP")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_PDF_ALPHAS_U", Vars.var_PDF_ALPHAS_U, b_pptr);
    weightVar = &Vars.var_PDF_ALPHAS_U;
  }
  else if(weightVarName == "PDF_ALPHAS_DOWN")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_PDF_ALPHAS_D", Vars.var_PDF_ALPHAS_D, b_pptr);
    weightVar = &Vars.var_PDF_ALPHAS_D;
  }
  else if(weightVarName == "PSSCALE_WEIGHT_4_UP"){
      SetBranchAddressVerbosely(TChainPtr, "var_PSSCALE_WEIGHT_4_U", Vars.var_PSSCALE_WEIGHT_4_U, b_pptr);
      weightVar = &Vars.var_PSSCALE_WEIGHT_4_U;
  }
  else if(weightVarName == "PSSCALE_WEIGHT_4_DOWN"){
      SetBranchAddressVerbosely(TChainPtr, "var_PSSCALE_WEIGHT_4_D", Vars.var_PSSCALE_WEIGHT_4_D, b_pptr);
      weightVar = &Vars.var_PSSCALE_WEIGHT_4_D;
  }
  else if(weightVarName == "PSSCALE_WEIGHT_5_UP"){
        SetBranchAddressVerbosely(TChainPtr, "var_PSSCALE_WEIGHT_5_U", Vars.var_PSSCALE_WEIGHT_5_U, b_pptr);
        weightVar = &Vars.var_PSSCALE_WEIGHT_5_U;
  }
  else if(weightVarName == "PSSCALE_WEIGHT_5_DOWN"){
        SetBranchAddressVerbosely(TChainPtr, "var_PSSCALE_WEIGHT_5_D", Vars.var_PSSCALE_WEIGHT_5_D, b_pptr);
        weightVar = &Vars.var_PSSCALE_WEIGHT_5_D;
  }
  else if(weightVarName == "L1PREFIRING_UP")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_L1PREFIRING_U", Vars.var_L1PREFIRING_U, b_pptr);
    weightVar = &Vars.var_L1PREFIRING_U;
  }
  else if(weightVarName == "L1PREFIRING_DOWN")
  {
    SetBranchAddressVerbosely(TChainPtr, "var_L1PREFIRING_D", Vars.var_L1PREFIRING_D, b_pptr);
    weightVar = &Vars.var_L1PREFIRING_D;
  }
  /*else if(TString(weightVarName, 4) == "PDF_")    //FIXME required
  {
    int n = TString(weightVarName.Data() + 4).Atoi();
    SetBranchAddressVerbosely(TChainPtr, TString::Format("var_PDF_%d", n), *(Vars.var_PDF + n), b_pptr);
    weightVar = Vars.var_PDF + n;
    // (here shift Vars.var_PDF + n is not necessary)
  }*/
  else
    throw std::runtime_error(TString::Format("Error in TreePlain::InitReading(): unknown weightVarName %s", weightVarName.Data()).Data());
}

void TreePlain::Add(const TString& name, const double w, const bool doWeightVar, const int sampleType)
{
  assert(TChainPtr);
  //printf("adding tree with name = %s\n", name.Data());
  printf("adding tree with name = %s weight = %f\n", name.Data(), w);
  int ntrees = TChainPtr->GetNtrees();
  int ret = TChainPtr->Add(name, 0);
  (void)ret;
  //printf("TreePlain::Add(): name = %s  w = %e  sampleTyp = %d\n", name.Data(), w, sampleType);
  //printf("ret = %d\n", ret);
  // if exactly one tree added, add corresponding weight
  if(TChainPtr->GetNtrees() == ntrees + 1)
  {
    AddWeight(w);
    DoWeightVar(doWeightVar);
    zSampleType.push_back(sampleType);
  }
  else if(TChainPtr->GetNtrees() == ntrees)
    ;// empty tree was not addeed
  else
    throw std::runtime_error(TString::Format("Error in TreePlain::Add(): wrong number of trees [%d] after adding %s", TChainPtr->GetNtrees() - ntrees, name.Data()).Data());
}


void TreePlain::InitVar(const TString& nameVar, const TString& nameCut, std::map<const TreePlain*, boost::variant<std::shared_ptr<int>, std::shared_ptr<float>, std::vector<std::shared_ptr<float>> >> & zTrees)
{
  TBranch** b_pptr = NULL;

  TObjString* readIsSignalFromFile = dynamic_cast<TObjString*>(TChainPtr->GetCurrentFile()->Get("isSignal"));
  const bool isSignalSample = readIsSignalFromFile->GetString() == "1";

  //TODO: fix
  //this->cutVarFloat = NULL;
  //this->cutVarShort = NULL;
  //if(nameCut != "")
  //{
  //  if(nameCut.BeginsWith("nj"))
  //  {
  //    //printf("enabling cut %s\n", nameCut.Data());
  //    short* tmp = new short;
  //    short* tmpReturned = (short*) SetBranchAddressVerbosely(TChainPtr, nameCut.Data(), *tmp, b_pptr);
  //    if(tmp != tmpReturned)
  //    {
  //      // branch existed, another pointer returned
  //      //printf("branch existed\n");
  //      this->cutVarShort = tmpReturned;
  //      delete tmp;
  //    }
  //    else
  //    {
  //      //printf("branch not existed\n");
  //      this->cutVarShort = tmp;
  //    }
  //  }
  //  else
  //    //SetBranchAddressVerbosely(TChainPtr, nameCut, this->cutVarFloat, b_pptr);
  //    throw std::logic_error(TString::Format("Error: cut in plainTree on variable %s not implemented", nameCut.Data()));
  //}

  // TODO fix signed, enabled jets
  if(nameVar == "ptt" || nameVar == "visptt"
     || nameVar == "ptat" || nameVar == "visptat"
     || nameVar == "pttLead" || nameVar == "vispttLead"
     || nameVar == "pttNLead" || nameVar == "vispttNLead"
     || nameVar == "pttTTRestFrame" || nameVar == "vispttTTRestFrame"
     || nameVar == "ptatTTRestFrame" || nameVar == "visptatTTRestFrame"
     || nameVar == "mtt" || nameVar == "vismtt"
     || nameVar == "pttt" || nameVar == "vispttt"
     || nameVar == "ytt" || nameVar == "ytts" || nameVar == "visytt" || nameVar == "visytts"
     || nameVar == "yt" || nameVar == "yts" || nameVar == "visyt" || nameVar == "visyts"
     || nameVar == "yat" || nameVar == "yats" || nameVar == "visyat" || nameVar == "visyats"
     || nameVar == "ytLead" || nameVar == "visytLead"
     || nameVar == "ytNLead" || nameVar == "visytNLead"
     || nameVar == "detatt" || nameVar == "visdetatt"
     || nameVar == "dytt" || nameVar == "visdytt"
     || nameVar == "dphitt" || nameVar == "visdphitt"
     || nameVar == "recoPartonMomFraction" || nameVar == "visrecoPartonMomFraction"
     || nameVar == "recoAntipartonMomFraction" || nameVar == "visrecoAntipartonMomFraction"
     || nameVar == "visetal"
     || nameVar == "visetaal"
     || nameVar == "visetalLead"
     || nameVar == "visetalNLead"
     || nameVar == "visdetall"
     || nameVar == "visdphill"
     || nameVar == "visptl"
     || nameVar == "visptal"
     || nameVar == "visptlLead"
     || nameVar == "visptlNLead"
     || nameVar == "visyll"
     || nameVar == "visptll"
     || nameVar == "vismll"
     || nameVar == "visetabLead"
     || nameVar == "visetabNLead"
     || nameVar == "visptbLead"
     || nameVar == "visptbNLead"
     || nameVar == "visybb"
     || nameVar == "visptbb"
     || nameVar == "vismbb"
     || nameVar == "visyllbb"
     || nameVar == "visptllbb"
     || nameVar == "vismllbb"
     || nameVar == "visyllbbmet"
     || nameVar == "visptllbbmet"
     || nameVar == "vismllbbmet"
     || nameVar == "visminmlb"
     || nameVar == "mttplusej1"
     || nameVar == "ptttplusej1"
     || nameVar == "r_mttplusej1_mtt"
     || nameVar == "r_ptej1_mtt"
     || nameVar == "mean_ptt_pttbar_ptej1"
     || nameVar == "detatej1_detattbar"){

    TString name = nameVar;
    if(nameVar.BeginsWith("vis")){
      if(!(isGen_ && isSignalSample))
	name.ReplaceAll("vis", "");
    }

    if(name == "yts")
      name.ReplaceAll("yts", "yt");
    if(name == "visyts")
      name.ReplaceAll("visyts", "visyt");

    if(name == "yats")
      name.ReplaceAll("yats", "yat");
    if(name == "visyats")
      name.ReplaceAll("visyats", "visyat");

    if(name == "ytts")
      name.ReplaceAll("ytts", "ytt");
    if(name == "visytts")
      name.ReplaceAll("visytts", "visytt");

    if(name == "minmlb")
      name.ReplaceAll("minmlb", "minmlbdef1");


    if(initVars.find(name) != initVars.end())
      zTrees.insert({this, initVars[name]});
    else {
      auto varPtr = std::make_shared<float>(0.);
      SetBranchAddressVerbosely(TChainPtr, name, *varPtr, b_pptr);
      initVars.insert({name, varPtr});
      zTrees.insert({this, varPtr});

    }

  }
  else if(nameVar.BeginsWith("yttKR") || nameVar.BeginsWith("yttsKR") || nameVar.BeginsWith("visyttKR") || nameVar.BeginsWith("visyttsKR")) {

    TString name = nameVar;
    bool DoParticle = false;
    if(nameVar.BeginsWith("vis")){
      DoParticle = true;
      if(!(isGen_ && isSignalSample))
	name.ReplaceAll("vis", "");
    }

    if(DoParticle && isSignalSample)
      name = isGen_ ? "visytt" : name;
    else
      name = isGen_ ? "ytt" : name;
    if(name.BeginsWith("yttsKR"))
      name.ReplaceAll("yttsKR", "yttKR");
    if(name.BeginsWith("visyttsKR"))
      name.ReplaceAll("visyttsKR", "visyttKR");

    if(initVars.find(name) != initVars.end())
      zTrees.insert({this, initVars[name]});
    else {
      auto varPtr = std::make_shared<float>(0.);
      SetBranchAddressVerbosely(TChainPtr, name, *varPtr, b_pptr);
      initVars.insert({name, varPtr});
      zTrees.insert({this, varPtr});

    }

  }
  else if(nameVar.BeginsWith("mttKR") || nameVar.BeginsWith("vismttKR")) {

    TString name = nameVar;
    bool DoParticle = false;
    if(nameVar.BeginsWith("vis")){
      DoParticle = true;
      if(!(isGen_ && isSignalSample))
	name.ReplaceAll("vis", "");
    }

    if(DoParticle && isSignalSample)
      name = isGen_ ? "vismtt" : name;
    else
      name = isGen_ ? "mtt" : name;


    if(initVars.find(name) != initVars.end())
      zTrees.insert({this, initVars[name]});
    else {
      auto varPtr = std::make_shared<float>(0.);
      SetBranchAddressVerbosely(TChainPtr, name, *varPtr, b_pptr);
      initVars.insert({name, varPtr});
      zTrees.insert({this, varPtr});

    }

  }
  else if(nameVar.BeginsWith("ptttKR") || nameVar.BeginsWith("visptttKR")) {

    TString name = nameVar;
    bool DoParticle = false;
    if(nameVar.BeginsWith("vis")){
      DoParticle = true;
      if(!(isGen_ && isSignalSample))
	name.ReplaceAll("vis", "");
    }

    if(DoParticle && isSignalSample)
      name = isGen_ ? "vispttt" : name;
    else
      name = isGen_ ? "pttt" : name;

    if(initVars.find(name) != initVars.end())
      zTrees.insert({this, initVars[name]});
    else {
      auto varPtr = std::make_shared<float>(0.);
      SetBranchAddressVerbosely(TChainPtr, name, *varPtr, b_pptr);
      initVars.insert({name, varPtr});
      zTrees.insert({this, varPtr});

    }

  }
  else if(nameVar == "visratioPtbLeadt"
	  || nameVar == "visratioPtbNLeadt"
	  || nameVar == "visratioPtlt"
	  || nameVar == "ratioPttMtt" || nameVar == "visratioPttMtt"
	  || nameVar == "ratioPtttMtt" || nameVar == "visratioPtttMtt")
  {

    TString name1;
    TString name2;
    if(nameVar == "visratioPtbLeadt") {
      name1 = "visptbLead";
      name2 = "visptt";
    }
    else if(nameVar == "visratioPtbNLeadt"){
      name1 = "visptbNLead";
      name2 = "visptt";
    }
    else if(nameVar == "visratioPtlt"){
      name1 = "visptl";
      name2 = "visptt";
    }
    else if(nameVar == "ratioPttMtt"){
      name1 = "ptt";
      name2 = "mtt";
    }
    else if(nameVar == "visratioPttMtt"){
      name1 = "visptt";
      name2 = "vismtt";
    }
    else if(nameVar == "ratioPtttMtt"){
      name1 = "pttt";
      name2 = "mtt";
    }
    else if(nameVar == "visratioPtttMtt"){
      name1 = "vispttt";
      name2 = "vismtt";
    }

    if(nameVar.BeginsWith("vis")){
      if(!(isGen_ && isSignalSample)){
	name1.ReplaceAll("vis", "");
	name2.ReplaceAll("vis", "");
      }
    }

    bool found_name1 = false;
    bool found_name2 = false;

    if(initVars.find(name1) != initVars.end())
      found_name1 = true;

    if(initVars.find(name2) != initVars.end())
      found_name2 = true;

    if(found_name1 && !found_name2){

	auto varPtr2 = std::make_shared<float>(0.);
	SetBranchAddressVerbosely(TChainPtr, name2, *varPtr2, b_pptr);
	std::vector<std::shared_ptr<float>> vvarPtr;
	vvarPtr.push_back(boost::get<std::shared_ptr<float>>(initVars[name1]));
	vvarPtr.push_back(varPtr2);
	initVars.insert({name2, varPtr2});
	zTrees.insert({this, vvarPtr});

    }
    else if(!found_name1 && found_name2)
      {
	auto varPtr1 = std::make_shared<float>(0.);
	SetBranchAddressVerbosely(TChainPtr, name1, *varPtr1, b_pptr);
	std::vector<std::shared_ptr<float>> vvarPtr;
	vvarPtr.push_back(varPtr1);
	vvarPtr.push_back(boost::get<std::shared_ptr<float>>(initVars[name2]));
        initVars.insert({name1, varPtr1});
	zTrees.insert({this, vvarPtr});
      }
    else if(found_name1 && found_name2)
      {
	std::vector<std::shared_ptr<float>> vvarPtr;
	vvarPtr.push_back(boost::get<std::shared_ptr<float>>(initVars[name1]));
        vvarPtr.push_back(boost::get<std::shared_ptr<float>>(initVars[name2]));
	zTrees.insert({this, vvarPtr});
      }
    else {

      auto varPtr1 = std::make_shared<float>(0.);
      auto varPtr2 = std::make_shared<float>(0.);
      SetBranchAddressVerbosely(TChainPtr, name1, *varPtr1, b_pptr);
      SetBranchAddressVerbosely(TChainPtr, name2, *varPtr2, b_pptr);
      std::vector<std::shared_ptr<float>> vvarPtr;
      vvarPtr.push_back(varPtr1);
      vvarPtr.push_back(varPtr2);
      initVars.insert({name1, varPtr1});
      initVars.insert({name2, varPtr2});
      zTrees.insert({this, vvarPtr});
    }

  }
  else if(nameVar.BeginsWith("nj")
	  || nameVar == "nvtx"
	  || nameVar == "nbj"
	  || nameVar == "naj")
  {
    //SetBranchAddressVerbosely(TChainPtr, name, this->Vars.dummyShort, b_pptr);
    //return &this->Vars.dummyShort;
    // TODO prevent memory leak

    if(initVars.find(nameVar) != initVars.end())
      zTrees.insert({this, initVars[nameVar]});
    else {
      auto varPtr = std::make_shared<int>(0.);
      SetBranchAddressVerbosely(TChainPtr, nameVar, *varPtr, b_pptr);
      initVars.insert({nameVar, varPtr});
      zTrees.insert({this, varPtr});

    }
  }
  else if(nameVar.BeginsWith("rhoj")
	  || nameVar.BeginsWith("mttj")
	  || nameVar.BeginsWith("yttj")
	  || nameVar == "met"
	  || nameVar == "nbj"
	  || nameVar == "ptj1"
	  || nameVar == "etaj1"
	  || nameVar == "ptj2"
	  || nameVar == "etaj2"
	  || nameVar == "ptej1"
	  || nameVar == "ptej2"
	  || nameVar == "etaej1"
	  || nameVar == "etaej2")
  {
    //printf("name = %s\n", name.Data());
    //SetBranchAddressVerbosely(TChainPtr, name, this->Vars.dummyFloat, b_pptr);
    //return &this->Vars.dummyFloat;
    // TODO prevent memory leak

    if(initVars.find(nameVar) != initVars.end())
      zTrees.insert({this, initVars[nameVar]});
    else {
      auto varPtr = std::make_shared<float>(0.5);
      SetBranchAddressVerbosely(TChainPtr, nameVar, *varPtr, b_pptr);
      initVars.insert({nameVar, varPtr});
      zTrees.insert({this, varPtr});

    }
  }
  else
    throw std::logic_error(TString::Format("Error in TreePlain::GetVar(): unknown variable %s\n", nameVar.Data()).Data());
}

double TreePlain::Weight(bool flagApplyCut) const
{
  // TODO somehow it is slower with saved weight
  //static int entryRead = -1;
  static double weight = -1.0;
  //if(entryRead == TTreePtr->GetReadEntry())
  //  return weight;
  weight = Vars.weight;
  weight *= zWeight[TChainPtr->GetTreeNumber()];
  if(weightVar && zDoWeightVar[TChainPtr->GetTreeNumber()])
    weight *= (*weightVar);
  // TODO investigate source of true level weight = 0
  if(weight != weight)
  {
    //printf("dupa\n");
    weight = 0.0;
  }
  // apply optional cut
  if(flagApplyCut && this->cutVarShort && !(*this->cutVarShort))
  {
    //printf("(*this->cutVarShort) = %f\n", (int)(*this->cutVarShort));
    weight = 0.0;
  }
  //if(this->cutVarFloat && !(*this->cutVarFloat))
  //  weight = 0.0;

  //entryRead = TTreePtr->GetReadEntry();
  //weight = 1.0;
  //printf("weight: %e\n", weight);
  return weight;
}

int TreePlain::SampleType(int n) const
{
  if(n < 0)
    n = TChainPtr->GetTreeNumber();
  return zSampleType[n];
}

void TreePlain::AddWeight(const double w)
{
  assert(TChainPtr);
  if(TChainPtr->GetNtrees() != (int)(zWeight.size() + 1))
  {
    printf("Error in TreePlain::AddWeight(): TChainPtr->GetNtrees() = %d, zWeight.size() = %ld\n", TChainPtr->GetNtrees(), zWeight.size());
    throw;
    return;
  }
  //printf("AddWeight(): pushing %f\n", w);
  zWeight.push_back(w);
}

void TreePlain::DoWeightVar(const bool doWeightVar)
{
  assert(TChainPtr);
  if(TChainPtr->GetNtrees() != (int)(zDoWeightVar.size() + 1))
  {
    printf("Error in TreePlain::DoWeightVar(): TChainPtr->GetNtrees() = %d, zDoWeightVar.size() = %ld\n", TChainPtr->GetNtrees(), zDoWeightVar.size());
    throw;
    return;
  }
  //printf("DoWeightVar(): pushing %d\n", doWeightVar);
  zDoWeightVar.push_back(doWeightVar);
}

template<class T>
void* TreePlain::SetBranchAddressVerbosely(TTree* tree, const char* name, T& ptr, TBranch** br)
{
  // do not anything if branch is already active (prevent seeting address of branch to different pointers)
  int statusBranch = tree->GetBranchStatus(name);
  if(statusBranch)
  {
    //printf("statusBranch = %d\n", statusBranch);
    assert(statusBranch == 1);
    TBranch* br = tree->FindBranch(name);
    return br->GetAddress();
  }
  //printf("Accessing branch %20s ... ", name);
  tree->SetBranchStatus(name, 1);
  auto ret = tree->SetBranchAddress(name, &ptr, br);
  if(ret == 0)
    ;//printf("success\n");
  else
  {
    printf("Accessing branch %20s ... ", name);printf("failed\n");
    printf("%s\n", tree->GetCurrentFile()->GetName());
    throw std::logic_error(TString::Format("Error while accessing branch %s: returned %d\n", name, ret).Data());
  }
  return &ptr;
}

// PORT end ==============================


void TreePlain::InitWriteTree()
{

  // common for everything
  this->AddBranchVerbosely("eventCounter", &this->Vars.eventCounter);
  this->AddBranchVerbosely("weight", &this->Vars.weight);

  if(fillBasicVariables_){
    AddBasicBranches();
  }
  if(fillParticleLevelVariables_){
    AddParticleLevelBranches();
  }
  if(fillMeasureVariables_){
    AddMeasureBranches();
  }
  if(fillGenVariationWeights_){
    AddGenVarBranches();
  }
  if(fillPUVarWeights_){
    AddPUVarBranches();
  }
  if(fillRecoVariationWeights_){
    AddRecoVarBranches();
  }
  if(fillLooseKinRecoVariables_){
    AddLooseKinRecoBranches();
  }
  if(fillJetVariables_){
    AddJetBranches();
  }
  if(fillPseudoVariables_){
    AddPseudoBranches();
  }
}



void TreePlain::AddBasicBranches()
{
  this->AddBranchVerbosely("ptt", &this->Vars.ptt);
  this->AddBranchVerbosely("ptat", &this->Vars.ptat);
  this->AddBranchVerbosely("pttLead", &this->Vars.pttLead);
  this->AddBranchVerbosely("pttNLead", &this->Vars.pttNLead);
  this->AddBranchVerbosely("pttTTRestFrame", &this->Vars.pttTTRestFrame);
  this->AddBranchVerbosely("ptatTTRestFrame", &this->Vars.ptatTTRestFrame);
  this->AddBranchVerbosely("yt", &this->Vars.yt);
  this->AddBranchVerbosely("yat", &this->Vars.yat);
  this->AddBranchVerbosely("ytLead", &this->Vars.ytLead);
  this->AddBranchVerbosely("ytNLead", &this->Vars.ytNLead);
  this->AddBranchVerbosely("mtt", &this->Vars.mtt);
  this->AddBranchVerbosely("ytt", &this->Vars.ytt);
  this->AddBranchVerbosely("pttt", &this->Vars.pttt);
  this->AddBranchVerbosely("dphitt", &this->Vars.dphitt);
  this->AddBranchVerbosely("detatt", &this->Vars.detatt);
  this->AddBranchVerbosely("dytt", &this->Vars.dytt);
  this->AddBranchVerbosely("recoPartonMomFraction", &this->Vars.recoPartonMomFraction);
  this->AddBranchVerbosely("recoAntipartonMomFraction", &this->Vars.recoAntipartonMomFraction);
  this->AddBranchVerbosely("nvtx", &this->Vars.nvtx);
  this->AddBranchVerbosely("nbj", &this->Vars.nbj);
  this->AddBranchVerbosely("naj", &this->Vars.naj);
  this->AddBranchVerbosely("met", &this->Vars.met);
  this->AddBranchVerbosely("ptj1", &this->Vars.ptj1);
  this->AddBranchVerbosely("etaj1", &this->Vars.etaj1);
  this->AddBranchVerbosely("ptj2", &this->Vars.ptj2);
  this->AddBranchVerbosely("etaj2", &this->Vars.etaj2);
}

void TreePlain::AddParticleLevelBranches()
{
    this->AddBranchVerbosely("ptl", &this->Vars.ptl);
    this->AddBranchVerbosely("ptal", &this->Vars.ptal);
    this->AddBranchVerbosely("ptlLead", &this->Vars.ptlLead);
    this->AddBranchVerbosely("ptlNLead", &this->Vars.ptlNLead);
    this->AddBranchVerbosely("etal", &this->Vars.etal);
    this->AddBranchVerbosely("etaal", &this->Vars.etaal);
    this->AddBranchVerbosely("etalLead", &this->Vars.etalLead);
    this->AddBranchVerbosely("etalNLead", &this->Vars.etalNLead);
    this->AddBranchVerbosely("detall", &this->Vars.detall);
    this->AddBranchVerbosely("dphill", &this->Vars.dphill);
    this->AddBranchVerbosely("ptll", &this->Vars.ptll);
    this->AddBranchVerbosely("yll", &this->Vars.yll);
    this->AddBranchVerbosely("mll", &this->Vars.mll);
    this->AddBranchVerbosely("ptllbb", &this->Vars.ptllbb);
    this->AddBranchVerbosely("yllbb", &this->Vars.yllbb);
    this->AddBranchVerbosely("mllbb", &this->Vars.mllbb);
    this->AddBranchVerbosely("ptllbbmet", &this->Vars.ptllbbmet);
    this->AddBranchVerbosely("yllbbmet", &this->Vars.yllbbmet);
    this->AddBranchVerbosely("mllbbmet", &this->Vars.mllbbmet);
    this->AddBranchVerbosely("ptbLead", &this->Vars.ptbLead);
    this->AddBranchVerbosely("ptbNLead", &this->Vars.ptbNLead);
    this->AddBranchVerbosely("etabLead", &this->Vars.etabLead);
    this->AddBranchVerbosely("etabNLead", &this->Vars.etabNLead);
    this->AddBranchVerbosely("ybb", &this->Vars.ybb);
    this->AddBranchVerbosely("ptbb", &this->Vars.ptbb);
    this->AddBranchVerbosely("mbb", &this->Vars.mbb);
    this->AddBranchVerbosely("matchedPtFullKR", &this->Vars.matchedPtFullKR[0]);
    this->AddBranchVerbosely("matchedRFullKR", &this->Vars.matchedRFullKR[0]);
    this->AddBranchVerbosely("matchedPtFullKR2", &this->Vars.matchedPtFullKR[1]);
    this->AddBranchVerbosely("matchedRFullKR2", &this->Vars.matchedRFullKR[1]);
    this->AddBranchVerbosely("mww", &this->Vars.mww);
    this->AddBranchVerbosely("minmlbdef1", &this->Vars.minmlbdef1);
    this->AddBranchVerbosely("minmlbdef2", &this->Vars.minmlbdef2);
}

void TreePlain::AddMeasureBranches()
{
  this->AddBranchVerbosely("mttmeas", &this->Vars.mttmeas);
  this->AddBranchVerbosely("yttmeas", &this->Vars.yttmeas);
  this->AddBranchVerbosely("ptttmeas", &this->Vars.ptttmeas);
  this->AddBranchVerbosely("ttbaremeas", &this->Vars.ttbaremeas);
  this->AddBranchVerbosely("ttbarzmeas", &this->Vars.ttbarzmeas);
  this->AddBranchVerbosely("ttbare", &this->Vars.ttbare);
  this->AddBranchVerbosely("ttbarz", &this->Vars.ttbarz);
  this->AddBranchVerbosely("trz", &this->Vars.trz);
  this->AddBranchVerbosely("nnbare", &this->Vars.nnbare);
  this->AddBranchVerbosely("nnbarz", &this->Vars.nnbarz);
  this->AddBranchVerbosely("ttbaremeasnomet", &this->Vars.ttbaremeasnomet);
  this->AddBranchVerbosely("llbare", &this->Vars.llbare);
  this->AddBranchVerbosely("llbarz", &this->Vars.llbarz);

}

void TreePlain::AddGenVarBranches()
{
  this->AddBranchVerbosely("var_MERENSCALE_U", &this->Vars.var_MERENSCALE_U);
  this->AddBranchVerbosely("var_MERENSCALE_D", &this->Vars.var_MERENSCALE_D);
  this->AddBranchVerbosely("var_MEFACSCALE_U", &this->Vars.var_MEFACSCALE_U);
  this->AddBranchVerbosely("var_MEFACSCALE_D", &this->Vars.var_MEFACSCALE_D);
  this->AddBranchVerbosely("var_MESCALE_U", &this->Vars.var_MESCALE_U);
  this->AddBranchVerbosely("var_MESCALE_D", &this->Vars.var_MESCALE_D);
  this->AddBranchVerbosely("var_BFRAG_U", &this->Vars.var_BFRAG_U);
  // this->AddBranchVerbosely("var_BFRAG_C", &this->Vars.var_BFRAG_C);
  this->AddBranchVerbosely("var_BFRAG_D", &this->Vars.var_BFRAG_D);
  this->AddBranchVerbosely("var_BFRAG_P", &this->Vars.var_BFRAG_P);
  this->AddBranchVerbosely("var_BSEMILEP_U", &this->Vars.var_BSEMILEP_U);
  this->AddBranchVerbosely("var_BSEMILEP_D", &this->Vars.var_BSEMILEP_D);
  this->AddBranchVerbosely("var_PDF_ALPHAS_U", &this->Vars.var_PDF_ALPHAS_U);
  this->AddBranchVerbosely("var_PDF_ALPHAS_D", &this->Vars.var_PDF_ALPHAS_D);
  this->AddBranchVerbosely("var_PSSCALE_WEIGHT_4_U", &this->Vars.var_PSSCALE_WEIGHT_4_U);
  this->AddBranchVerbosely("var_PSSCALE_WEIGHT_4_D", &this->Vars.var_PSSCALE_WEIGHT_4_D);
  this->AddBranchVerbosely("var_PSSCALE_WEIGHT_5_U", &this->Vars.var_PSSCALE_WEIGHT_5_U);
  this->AddBranchVerbosely("var_PSSCALE_WEIGHT_5_D", &this->Vars.var_PSSCALE_WEIGHT_5_D);

  for(int w = 1; w <= 159; w++){
    this->AddBranchVerbosely(TString::Format("var_PDF_%d", w), &this->Vars.var_PDF[w]);
  }
}

void TreePlain::AddPUVarBranches()
{
  this->AddBranchVerbosely("var_PU_U", &this->Vars.var_PU_U);
  this->AddBranchVerbosely("var_PU_D", &this->Vars.var_PU_D);
}
void TreePlain::AddRecoVarBranches()
{
  this->AddBranchVerbosely("var_TRIG_U", &this->Vars.var_TRIG_U);
  this->AddBranchVerbosely("var_TRIG_D", &this->Vars.var_TRIG_D);
  // this->AddBranchVerbosely("var_TRIG_ETA_U", &this->Vars.var_TRIG_ETA_U);
  // this->AddBranchVerbosely("var_TRIG_ETA_D", &this->Vars.var_TRIG_ETA_D);
  // this->AddBranchVerbosely("var_LEPT_U", &this->Vars.var_LEPT_U);
  // this->AddBranchVerbosely("var_LEPT_D", &this->Vars.var_LEPT_D);
  this->AddBranchVerbosely("var_ELE_ID_U", &this->Vars.var_ELE_ID_U);
  this->AddBranchVerbosely("var_ELE_ID_D", &this->Vars.var_ELE_ID_D);
  this->AddBranchVerbosely("var_ELE_RECO_U", &this->Vars.var_ELE_RECO_U);
  this->AddBranchVerbosely("var_ELE_RECO_D", &this->Vars.var_ELE_RECO_D);
  this->AddBranchVerbosely("var_MUON_ID_U", &this->Vars.var_MUON_ID_U);
  this->AddBranchVerbosely("var_MUON_ID_D", &this->Vars.var_MUON_ID_D);
  this->AddBranchVerbosely("var_MUON_ISO_U", &this->Vars.var_MUON_ISO_U);
  this->AddBranchVerbosely("var_MUON_ISO_D", &this->Vars.var_MUON_ISO_D);
  this->AddBranchVerbosely("var_KIN_U", &this->Vars.var_KIN_U);
  this->AddBranchVerbosely("var_KIN_D", &this->Vars.var_KIN_D);
  this->AddBranchVerbosely("var_BTAG_U", &this->Vars.var_BTAG_U);
  this->AddBranchVerbosely("var_BTAG_D", &this->Vars.var_BTAG_D);
  this->AddBranchVerbosely("var_BTAG_LJET_U", &this->Vars.var_BTAG_LJET_U);
  this->AddBranchVerbosely("var_BTAG_LJET_D", &this->Vars.var_BTAG_LJET_D);
  // this->AddBranchVerbosely("var_BTAG_PT_U", &this->Vars.var_BTAG_PT_U);
  // this->AddBranchVerbosely("var_BTAG_PT_D", &this->Vars.var_BTAG_PT_D);
  // this->AddBranchVerbosely("var_BTAG_ETA_U", &this->Vars.var_BTAG_ETA_U);
  // this->AddBranchVerbosely("var_BTAG_ETA_D", &this->Vars.var_BTAG_ETA_D);
  // this->AddBranchVerbosely("var_BTAG_LJET_PT_U", &this->Vars.var_BTAG_LJET_PT_U);
  // this->AddBranchVerbosely("var_BTAG_LJET_PT_D", &this->Vars.var_BTAG_LJET_PT_D);
  // this->AddBranchVerbosely("var_BTAG_LJET_ETA_U", &this->Vars.var_BTAG_LJET_ETA_U);
  // this->AddBranchVerbosely("var_BTAG_LJET_ETA_D", &this->Vars.var_BTAG_LJET_ETA_D);
  this->AddBranchVerbosely("var_L1PREFIRING_U", &this->Vars.var_L1PREFIRING_U);
  this->AddBranchVerbosely("var_L1PREFIRING_D", &this->Vars.var_L1PREFIRING_D);
}
void TreePlain::AddLooseKinRecoBranches()
{
  this->AddBranchVerbosely("matchedPt", &this->Vars.matchedPt);
  this->AddBranchVerbosely("matchedR", &this->Vars.matchedR);
  // this->AddBranchVerbosely("mttKR2", &this->Vars.mttKR2);
  // this->AddBranchVerbosely("yttKR2", &this->Vars.yttKR2);
  // this->AddBranchVerbosely("mwwKR2", &this->Vars.mwwKR2);
  // this->AddBranchVerbosely("ptttKR2", &this->Vars.ptttKR2);
  // this->AddBranchVerbosely("minmlbKR2", &this->Vars.minmlbKR2);
  // this->AddBranchVerbosely("mttKR6", &this->Vars.mttKR6);
  // this->AddBranchVerbosely("yttKR6", &this->Vars.yttKR6);
  // this->AddBranchVerbosely("ptttKR6", &this->Vars.ptttKR6);
  // this->AddBranchVerbosely("mwwKR6", &this->Vars.mwwKR6);
  // this->AddBranchVerbosely("minmlbKR6", &this->Vars.minmlbKR6);
  this->AddBranchVerbosely("mttKR9", &this->Vars.mttKR9);
  this->AddBranchVerbosely("yttKR9", &this->Vars.yttKR9);
  this->AddBranchVerbosely("ptttKR9", &this->Vars.ptttKR9);
  this->AddBranchVerbosely("mwwKR9", &this->Vars.mwwKR9);
  this->AddBranchVerbosely("minmlbKR9", &this->Vars.minmlbKR9);
}
void TreePlain::AddJetBranches()
{
  this->AddBranchVerbosely("njdefIso04Pt30", &this->Vars.njdefIso04Pt30);
  this->AddBranchVerbosely("njdefIso04Pt40", &this->Vars.njdefIso04Pt40);
  this->AddBranchVerbosely("njdefIso04Pt50", &this->Vars.njdefIso04Pt50);
  this->AddBranchVerbosely("njdefIso04Pt75", &this->Vars.njdefIso04Pt75);
  this->AddBranchVerbosely("njdefIso04Pt100", &this->Vars.njdefIso04Pt100);
  this->AddBranchVerbosely("njdefIso04Pt150", &this->Vars.njdefIso04Pt150);

  //leading and trailing extra jets
  this->AddBranchVerbosely("ptej1", &this->Vars.ptej1);
  this->AddBranchVerbosely("ptej2", &this->Vars.ptej2);
  this->AddBranchVerbosely("etaej1", &this->Vars.etaej1);
  this->AddBranchVerbosely("etaej2", &this->Vars.etaej2);

    //ttbar + ej1
    this->AddBranchVerbosely("mttplusej1", &this->Vars.mttplusej1);
    this->AddBranchVerbosely("ptttplusej1", &this->Vars.ptttplusej1);
    this->AddBranchVerbosely("r_mttplusej1_mtt", &this->Vars.r_mttplusej1_mtt);
    this->AddBranchVerbosely("r_ptej1_mtt", &this->Vars.r_ptej1_mtt);


    // Not possible for looseKR
    if(!fillLooseKinRecoVariables_){
      //multiplicity of jets inside ttbar (in pseudo-rapidity)
      this->AddBranchVerbosely("njinetattbarIso04Pt30", &this->Vars.njinetattbarIso04Pt30);
      this->AddBranchVerbosely("njinetattbarIso04Pt40", &this->Vars.njinetattbarIso04Pt40);
      this->AddBranchVerbosely("njinetattbarIso04Pt50", &this->Vars.njinetattbarIso04Pt50);
      this->AddBranchVerbosely("njinetattbarIso04Pt75", &this->Vars.njinetattbarIso04Pt75);
      this->AddBranchVerbosely("njinetattbarIso04Pt100", &this->Vars.njinetattbarIso04Pt100);
      this->AddBranchVerbosely("njinetattbarIso04Pt150", &this->Vars.njinetattbarIso04Pt150);

      //ttbar + ej1
      this->AddBranchVerbosely("mean_ptt_pttbar_ptej1", &this->Vars.mean_ptt_pttbar_ptej1);

      //see if jets within the ttbar system
      this->AddBranchVerbosely("detatej1_detattbar", &this->Vars.detatej1_detattbar);
      
    }
}

void TreePlain::AddPseudoBranches()
{
  if(fillRecoTree_ || fillGenTree_ || fillSimpleRecoTree_)
    this->AddBranchVerbosely("inFiducialPhaseSpace", &this->Vars.inFiducialPhaseSpace);

  if(fillGenTree_){
    this->AddBranchVerbosely("visptt", &this->Vars.visptt);
    this->AddBranchVerbosely("visptat", &this->Vars.visptat);
    this->AddBranchVerbosely("vispttLead", &this->Vars.vispttLead);
    this->AddBranchVerbosely("vispttNLead", &this->Vars.vispttNLead);
    this->AddBranchVerbosely("vispttTTRestFrame", &this->Vars.vispttTTRestFrame);
    this->AddBranchVerbosely("visptatTTRestFrame", &this->Vars.visptatTTRestFrame);
    this->AddBranchVerbosely("visyt", &this->Vars.visyt);
    this->AddBranchVerbosely("visyat", &this->Vars.visyat);
    this->AddBranchVerbosely("visytLead", &this->Vars.visytLead);
    this->AddBranchVerbosely("visytNLead", &this->Vars.visytNLead);
    this->AddBranchVerbosely("vismtt", &this->Vars.vismtt);
    this->AddBranchVerbosely("visytt", &this->Vars.visytt);
    this->AddBranchVerbosely("vispttt", &this->Vars.vispttt);
    this->AddBranchVerbosely("visdphitt", &this->Vars.visdphitt);
    this->AddBranchVerbosely("visdetatt", &this->Vars.visdetatt);
    this->AddBranchVerbosely("visdytt", &this->Vars.visdytt);
    this->AddBranchVerbosely("visrecoPartonMomFraction", &this->Vars.visrecoPartonMomFraction);
    this->AddBranchVerbosely("visrecoAntipartonMomFraction", &this->Vars.visrecoAntipartonMomFraction);
    this->AddBranchVerbosely("visptl", &this->Vars.visptl);
    this->AddBranchVerbosely("visptal", &this->Vars.visptal);
    this->AddBranchVerbosely("visptlLead", &this->Vars.visptlLead);
    this->AddBranchVerbosely("visptlNLead", &this->Vars.visptlNLead);
    this->AddBranchVerbosely("visetal", &this->Vars.visetal);
    this->AddBranchVerbosely("visetaal", &this->Vars.visetaal);
    this->AddBranchVerbosely("visetalLead", &this->Vars.visetalLead);
    this->AddBranchVerbosely("visetalNLead", &this->Vars.visetalNLead);
    this->AddBranchVerbosely("visdetall", &this->Vars.visdetall);
    this->AddBranchVerbosely("visdphill", &this->Vars.visdphill);
    this->AddBranchVerbosely("visyll", &this->Vars.visyll);
    this->AddBranchVerbosely("visptll", &this->Vars.visptll);
    this->AddBranchVerbosely("vismll", &this->Vars.vismll);
    this->AddBranchVerbosely("visptllbb", &this->Vars.visptllbb);
    this->AddBranchVerbosely("visyllbb", &this->Vars.visyllbb);
    this->AddBranchVerbosely("vismllbb", &this->Vars.vismllbb);
    this->AddBranchVerbosely("visptllbbmet", &this->Vars.visptllbbmet);
    this->AddBranchVerbosely("visyllbbmet", &this->Vars.visyllbbmet);
    this->AddBranchVerbosely("vismllbbmet", &this->Vars.vismllbbmet);
    this->AddBranchVerbosely("visptbLead", &this->Vars.visptbLead);
    this->AddBranchVerbosely("visptbNLead", &this->Vars.visptbNLead);
    this->AddBranchVerbosely("visetabLead", &this->Vars.visetabLead);
    this->AddBranchVerbosely("visetabNLead", &this->Vars.visetabNLead);
    this->AddBranchVerbosely("visybb", &this->Vars.visybb);
    this->AddBranchVerbosely("visptbb", &this->Vars.visptbb);
    this->AddBranchVerbosely("vismbb", &this->Vars.vismbb);
    this->AddBranchVerbosely("vismww", &this->Vars.vismww);
    this->AddBranchVerbosely("visminmlb", &this->Vars.visminmlb);
  }
}

// void TreePlain::Fill(const AnalysisBaseEntry* entry, const float deltaRLepJet)
void TreePlain::Fill(const AnalysisBaseEntry* entry)
{
  // std::cout<<fillLooseKinRecoVariables_<<std::endl;
  // common for everything
  Vars.eventCounter = *entry->eventCounter;
  // weight needs to be adjusted for simpleKinRecoTree
  Vars.weight = fillGenTree_ ? entry->genLevelWeights->trueLevelWeight_ : *entry->weight;

  // if(IsGen_)
  // {
  //   Vars.weight = entry->genLevelWeights->trueLevelWeight_;
  //   this->FillGen(entry, deltaRLepJet);
  // }
  // else
  // {
  //   Vars.weight = *entry->weight;
  //   this->FillReco(entry);
  // }
  if(fillBasicVariables_){
    FillBasicBranches(entry);
  }
  if(fillMeasureVariables_){
    if(entry->topGenObjects->valuesSet_) FillMeasureBranches(entry);
  }
  if(fillGenVariationWeights_){
    FillGenVarBranches(entry);
  }
  if(fillPUVarWeights_){
    FillPUVarBranches(entry);
  }
  if(fillRecoVariationWeights_){
    FillRecoVarBranches(entry);
  }
  if(fillLooseKinRecoVariables_){
    FillLooseKinRecoBranches(entry);
  }
  if(fillJetVariables_){
    FillJetBranches(entry);
  }
  if(fillPseudoVariables_){
    FillPseudoBranches(entry);
  }
  TTreePtr->Fill();
}





void TreePlain::FillBasicBranches(const AnalysisBaseEntry* entry)
{
  //leading and trailing top
  LV LeadTop, NLeadTop;
  //leading and trailing lepton
  LV LeadLep, NLeadLep;
  //leading and trailing b-jet
  LV LeadBot, NLeadBot;
  // create top/antitop quark in the ttbar rest frame
  LV ttbarSystem;
  LV topRestFrameTTbar, antiTopRestFrameTTbar;
  if(fillRecoTree_){
    Vars.lvTop = common::LVtoTLV(entry->kinematicReconstructionSolutions->solution().top());
    Vars.lvTopBar = common::LVtoTLV(entry->kinematicReconstructionSolutions->solution().antiTop());
    Vars.lvLep = common::LVtoTLV(entry->kinematicReconstructionSolutions->solution().lepton());
    Vars.lvLepBar = common::LVtoTLV(entry->kinematicReconstructionSolutions->solution().antiLepton());
    Vars.lvBot = common::LVtoTLV(entry->kinematicReconstructionSolutions->solution().bjet());
    Vars.lvBotBar = common::LVtoTLV(entry->kinematicReconstructionSolutions->solution().antiBjet());
    Vars.lvMET = common::LVtoTLV(*entry->recoObjects->met_);

    //leading and trailing top
    common::orderLV(LeadTop, NLeadTop, entry->kinematicReconstructionSolutions->solution().top(), entry->kinematicReconstructionSolutions->solution().antiTop(), common::LVpt);
    Vars.lvLeadTop = common::LVtoTLV(LeadTop);
    Vars.lvNLeadTop = common::LVtoTLV(NLeadTop);

    //leading and trailing Lepton
    common::orderLV(LeadLep, NLeadLep, entry->kinematicReconstructionSolutions->solution().lepton(), entry->kinematicReconstructionSolutions->solution().antiLepton(), common::LVpt);
    Vars.lvLeadLep = common::LVtoTLV(LeadLep);
    Vars.lvNLeadLep = common::LVtoTLV(NLeadLep);

    //leading and trailing b-jet
    common::orderLV(LeadBot, NLeadBot, entry->kinematicReconstructionSolutions->solution().bjet(), entry->kinematicReconstructionSolutions->solution().antiBjet(), common::LVpt);
    Vars.lvLeadBot = common::LVtoTLV(LeadBot);
    Vars.lvNLeadBot = common::LVtoTLV(NLeadBot);

    // create top/antitop quark in the ttbar rest frame
    ttbarSystem = entry->kinematicReconstructionSolutions->solution().top() + entry->kinematicReconstructionSolutions->solution().antiTop();
    ROOT::Math::Boost CoMBoostHypTtbar (ttbarSystem.BoostToCM());
    topRestFrameTTbar = CoMBoostHypTtbar(entry->kinematicReconstructionSolutions->solution().top());
    antiTopRestFrameTTbar = CoMBoostHypTtbar(entry->kinematicReconstructionSolutions->solution().antiTop());
    Vars.lvTopTTRestFrame = common::LVtoTLV(topRestFrameTTbar);
    Vars.lvTopBarTTRestFrame = common::LVtoTLV(antiTopRestFrameTTbar);

    Vars.nvtx = entry->recoObjects->vertMultiGood_; //use vertMulti_ or vertMultiGood_?
    Vars.nbj = entry->recoObjectIndices->bjetIndices_.size();//? correct
    Vars.naj = entry->recoObjectIndices->jetIndices_.size();//? correct

    //quantities for particle-level only

    if(isGen_ && isTTbarSample_ && entry->topGenObjects->valuesSet_){
      std::pair<float, float> mat = CalculateMatchingPtR(entry, std::pair<TLorentzVector, TLorentzVector> (Vars.lvBot,Vars.lvBotBar));
      Vars.matchedPtFullKR[0] = mat.first;
      Vars.matchedRFullKR[0] = mat.second;
      // below taking into account correctlepton assignment (i.e. b assigned to bbar and vice versa does not work)
      std::pair<double, double> matFullKR = CalculateMatchingPtRFullKR(entry, Vars.lvBot, Vars.lvBotBar);
      Vars.matchedPtFullKR[1] = matFullKR.first;
      Vars.matchedRFullKR[1] = matFullKR.second;
    }

    TLorentzVector dilepton = Vars.lvLep + Vars.lvLepBar;
    TLorentzVector bbar = Vars.lvBot + Vars.lvBotBar;
    Vars.ptl = Vars.lvLep.Pt();
    Vars.ptal = Vars.lvLepBar.Pt();
    Vars.ptlLead = Vars.lvLeadLep.Pt();
    Vars.ptlNLead = Vars.lvNLeadLep.Pt();
    Vars.etal = Vars.lvLep.PseudoRapidity();
    Vars.etaal = Vars.lvLepBar.PseudoRapidity();
    Vars.etalLead = Vars.lvLeadLep.PseudoRapidity();
    Vars.etalNLead = Vars.lvNLeadLep.PseudoRapidity();
    Vars.detall = TMath::Abs(Vars.lvLep.PseudoRapidity()) - TMath::Abs(Vars.lvLepBar.PseudoRapidity());
    Vars.dphill = TMath::Abs(Vars.lvLep.DeltaPhi(Vars.lvLepBar));
    Vars.ptll = dilepton.Pt();
    Vars.yll = dilepton.Rapidity();
    Vars.mll = dilepton.M();
    Vars.ptllbb = (dilepton+bbar).Pt();
    Vars.yllbb = (dilepton+bbar).Rapidity();
    Vars.mllbb = (dilepton+bbar).M();
    Vars.ptllbbmet = (dilepton+bbar+Vars.lvMET).Pt();
    Vars.yllbbmet = (dilepton+bbar+Vars.lvMET).Rapidity();
    Vars.mllbbmet = (dilepton+bbar+Vars.lvMET).M();
    Vars.ptbLead = Vars.lvLeadBot.Pt();
    Vars.ptbNLead = Vars.lvNLeadBot.Pt();
    Vars.etabLead = Vars.lvLeadBot.PseudoRapidity();
    Vars.etabNLead = Vars.lvNLeadBot.PseudoRapidity();
    Vars.ptbb = bbar.Pt();
    Vars.ybb = bbar.Rapidity();
    Vars.mbb = bbar.M();

    //vars from kinReco
    Vars.mww = ((Vars.lvTop + Vars.lvTopBar) - Vars.lvBot - Vars.lvBotBar).M();
    Vars.minmlbdef1 = std::min((Vars.lvLep + Vars.lvBotBar).M(), (Vars.lvLepBar + Vars.lvBot).M());
    //if using same definition for mlb as for loose kin reco
    std::vector<float> mlb(4,0.);
    mlb[0] = (Vars.lvLep + Vars.lvBot).M();
    mlb[1] = (Vars.lvLep +  Vars.lvBotBar).M();
    mlb[2] = (Vars.lvLepBar + Vars.lvBot).M();
    mlb[3] = (Vars.lvLepBar + Vars.lvBotBar).M();
    float mlbmax1 = std::max(mlb[0], mlb[3]);
    float mlbmax2 = std::max(mlb[1], mlb[2]);
    Vars.minmlbdef2 = std::min(mlbmax1, mlbmax2);
  }
  if(fillGenTree_){
    if(isTTbarSample_){
      if(entry->topGenObjects->valuesSet_){
        Vars.lvTop = common::LVtoTLV(*entry->topGenObjects->GenTop_);
        Vars.lvTopBar = common::LVtoTLV(*entry->topGenObjects->GenAntiTop_);
        Vars.lvLep = common::LVtoTLV(*entry->topGenObjects->GenLepton_);
        Vars.lvLepBar = common::LVtoTLV(*entry->topGenObjects->GenAntiLepton_);
        Vars.lvBot = common::LVtoTLV(*entry->topGenObjects->GenB_);
        Vars.lvBotBar = common::LVtoTLV(*entry->topGenObjects->GenAntiB_);
        Vars.lvMET = common::LVtoTLV(*entry->topGenObjects->GenMet_);

	//leading and trailing top
	common::orderLV(LeadTop, NLeadTop, *entry->topGenObjects->GenTop_, *entry->topGenObjects->GenAntiTop_, common::LVpt);
	Vars.lvLeadTop = common::LVtoTLV(LeadTop);
	Vars.lvNLeadTop = common::LVtoTLV(NLeadTop);

	// create top/antitop quark in the ttbar rest frame
	ttbarSystem = *entry->topGenObjects->GenTop_ + *entry->topGenObjects->GenAntiTop_;
	ROOT::Math::Boost CoMBoostGenTtbar (ttbarSystem.BoostToCM());
	topRestFrameTTbar = CoMBoostGenTtbar(*entry->topGenObjects->GenTop_);
	antiTopRestFrameTTbar = CoMBoostGenTtbar(*entry->topGenObjects->GenAntiTop_);
	Vars.lvTopTTRestFrame = common::LVtoTLV(topRestFrameTTbar);
	Vars.lvTopBarTTRestFrame = common::LVtoTLV(antiTopRestFrameTTbar);

        Vars.nbj = entry->topGenObjects->genBHadJetIndex_->size(); //correct?
        Vars.naj = (int)(*entry->topGenObjects->allGenJets_).size(); //correct?
      }
      Vars.nvtx = entry->analysisBase->vertMultiTrue();
    }

  }
	  TLorentzVector ttbar = Vars.lvTop + Vars.lvTopBar;
	  Vars.ptt = Vars.lvTop.Pt();
	  Vars.ptat = Vars.lvTopBar.Pt();
	  Vars.pttLead = Vars.lvLeadTop.Pt();
	  Vars.pttNLead = Vars.lvNLeadTop.Pt();
	  Vars.pttTTRestFrame = Vars.lvTopTTRestFrame.Pt();
	  Vars.ptatTTRestFrame = Vars.lvTopBarTTRestFrame.Pt();
	  Vars.yt = Vars.lvTop.Rapidity();
	  Vars.yat = Vars.lvTopBar.Rapidity();
	  Vars.ytLead = Vars.lvLeadTop.Rapidity();
	  Vars.ytNLead = Vars.lvNLeadTop.Rapidity();
	  Vars.mtt = ttbar.M();
	  Vars.ytt = ttbar.Rapidity();
	  Vars.pttt = ttbar.Pt();
	  Vars.dphitt = TMath::Abs(Vars.lvTop.DeltaPhi(Vars.lvTopBar));
	  Vars.dytt = TMath::Abs(Vars.lvTop.Rapidity()) - TMath::Abs(Vars.lvTopBar.Rapidity());
	  Vars.detatt = TMath::Abs(Vars.lvTop.PseudoRapidity() - Vars.lvTopBar.PseudoRapidity());
	  Vars.met = Vars.lvMET.Pt();
	  Vars.recoPartonMomFraction = (Vars.lvTop.Energy() - Vars.lvTop.Pz() + Vars.lvTopBar.Energy() - Vars.lvTopBar.Pz()) / (2 * 6500);
	  Vars.recoAntipartonMomFraction = (Vars.lvTop.Energy() + Vars.lvTop.Pz() + Vars.lvTopBar.Energy() + Vars.lvTopBar.Pz()) / (2 * 6500);

	  //leading and sub-leading bjets
	  if(Vars.lvBot.Pt() > Vars.lvBotBar.Pt())
		{
		Vars.ptj1 = Vars.lvBot.Pt();
		Vars.etaj1 = Vars.lvBot.Eta();
		Vars.ptj2 = Vars.lvBotBar.Pt();
		Vars.etaj2 = Vars.lvBotBar.Eta();
	        }
	  else
	       {
	       Vars.ptj1 = Vars.lvBotBar.Pt();
	       Vars.etaj1 = Vars.lvBotBar.Eta();
	       Vars.ptj2 = Vars.lvBot.Pt();
	       Vars.etaj2 = Vars.lvBot.Eta();
	       }
  // }
}


void TreePlain::FillMeasureBranches(const AnalysisBaseEntry* entry)
{
	TLorentzVector ttbar = Vars.lvTop + Vars.lvTopBar;
	Vars.ttbare = ttbar.E();
	Vars.ttbarz = ttbar.Pz();

	Vars.v2MET = TVector2(entry->topGenObjects->GenMet_->Px(), entry->topGenObjects->GenMet_->Py());
	TLorentzVector lvGenMET(Vars.v2MET.Px(), Vars.v2MET.Py(), 0.0, Vars.v2MET.Mod());
	TLorentzVector nu = common::LVtoTLV(*entry->topGenObjects->GenNeutrino_);
	TLorentzVector nubar = common::LVtoTLV(*entry->topGenObjects->GenAntiNeutrino_);
	const TLorentzVector nnbar = nu + nubar;
	Vars.nnbare = nnbar.E();
	Vars.nnbarz = nnbar.Pz();

	TLorentzVector nunubarxy(nnbar.X(), nnbar.Y(), 0.0, TMath::Sqrt(nnbar.X() * nnbar.X() + nnbar.Y() * nnbar.Y()));

	TLorentzVector ttbarMeas = Vars.lvLep + Vars.lvLepBar + Vars.lvBot + Vars.lvBotBar + nunubarxy;
	Vars.mttmeas = ttbarMeas.M();
	Vars.yttmeas = ttbarMeas.Rapidity();
	Vars.ptttmeas = ttbarMeas.Pt();
	Vars.ttbaremeas = ttbarMeas.E();
	Vars.ttbarzmeas = ttbarMeas.Pz();

	TLorentzVector ttbarMeasNoMET = ttbarMeas - lvGenMET;
	Vars.ttbaremeasnomet = ttbarMeasNoMET.E();

	const TLorentzVector llbar = ttbar - Vars.lvBot - Vars.lvBotBar - nnbar;
	Vars.llbare = llbar.E();
	Vars.llbarz = llbar.Pz();

}


void TreePlain::FillGenVarBranches(const AnalysisBaseEntry* entry)
{
  // Long64_t entryNumber = entry->analysisBase->b_eventNumber->GetEntryNumber();
  Vars.var_MERENSCALE_U = entry->analysisBase->weightMeRenScale_UP(entry->entryNumber);
  Vars.var_MERENSCALE_D = entry->analysisBase->weightMeRenScale_DN(entry->entryNumber);
  Vars.var_MEFACSCALE_U = entry->analysisBase->weightMeFacScale_UP(entry->entryNumber);
  Vars.var_MEFACSCALE_D = entry->analysisBase->weightMeFacScale_DN(entry->entryNumber);
  Vars.var_MESCALE_U = entry->analysisBase->weightMeScale_UP(entry->entryNumber);
  Vars.var_MESCALE_D = entry->analysisBase->weightMeScale_DN(entry->entryNumber);
  Vars.var_BFRAG_U = entry->analysisBase->weightBFrag_Up(entry->entryNumber);
  // Vars.var_BFRAG_C = entry->analysisBase->weightBFrag_Central(entry->entryNumber);
  Vars.var_BFRAG_D = entry->analysisBase->weightBFrag_Dn(entry->entryNumber);
  Vars.var_BFRAG_P = entry->analysisBase->weightBFrag_Peterson(entry->entryNumber);
  Vars.var_BSEMILEP_U = entry->analysisBase->weightBSemilep_Up(entry->entryNumber);
  Vars.var_BSEMILEP_D = entry->analysisBase->weightBSemilep_Dn(entry->entryNumber);
  Vars.var_PDF_ALPHAS_U = entry->analysisBase->weightAlphasPdf_Up(entry->entryNumber);
  Vars.var_PDF_ALPHAS_D = entry->analysisBase->weightAlphasPdf_Dn(entry->entryNumber);

  if(entry->analysisBase->b_pdfWeights->GetEntryNumber() != *entry->eventCounter){
    entry->analysisBase->b_pdfWeights->GetEntry(*entry->eventCounter);
  }
  const std::vector<float>& pdfWeightsFloat = *(entry->analysisBase->pdfWeights_);
  for(unsigned int i = 0; i < pdfWeightsFloat.size(); i++){
    Vars.var_PDF[i]=pdfWeightsFloat[i];
  }
  
  const std::vector<float>& psWeightsFloat = *(entry->analysisBase->psWeights_);
  /*
  for(unsigned int i = 0; i < psWeightsFloat.size(); i++){
    Vars.var_PS.push_back(psWeightsFloat[i]);
  }*/

  Vars.var_PSSCALE_WEIGHT_4_U = psWeightsFloat.at(6)/psWeightsFloat.at(1);
  Vars.var_PSSCALE_WEIGHT_4_D = psWeightsFloat.at(8)/psWeightsFloat.at(1);
  Vars.var_PSSCALE_WEIGHT_5_U = psWeightsFloat.at(7)/psWeightsFloat.at(1);
  Vars.var_PSSCALE_WEIGHT_5_D = psWeightsFloat.at(9)/psWeightsFloat.at(1);

  // Use correct BtagEff for each syst
  if (fillRecoVariationWeights_) {
    Vars.var_MERENSCALE_U  *= GetBTagWeightForSyst(entry, "MERENSCALE_UP");
    Vars.var_MERENSCALE_D  *= GetBTagWeightForSyst(entry, "MERENSCALE_DOWN");
    Vars.var_MEFACSCALE_U  *= GetBTagWeightForSyst(entry, "MEFACSCALE_UP");
    Vars.var_MEFACSCALE_D  *= GetBTagWeightForSyst(entry, "MEFACSCALE_DOWN");
    Vars.var_MESCALE_U  *= GetBTagWeightForSyst(entry, "MESCALE_UP");
    Vars.var_MESCALE_D  *= GetBTagWeightForSyst(entry, "MESCALE_DOWN");
    Vars.var_BFRAG_U  *= GetBTagWeightForSyst(entry, "BFRAG_UP");
    Vars.var_BFRAG_D  *= GetBTagWeightForSyst(entry, "BFRAG_DOWN");
    Vars.var_BFRAG_P  *= GetBTagWeightForSyst(entry, "BFRAG_PETERSON");
    // Vars.var_BFRAG_C  *= GetBTagWeightForSyst(entry, "BFRAG_CENTRAL");
    Vars.var_BSEMILEP_U  *= GetBTagWeightForSyst(entry, "BSEMILEP_UP");
    Vars.var_BSEMILEP_D  *= GetBTagWeightForSyst(entry, "BSEMILEP_DOWN");
    Vars.var_PDF_ALPHAS_U  *= GetBTagWeightForSyst(entry, "PDF_ALPHAS_UP");
    Vars.var_PDF_ALPHAS_D  *= GetBTagWeightForSyst(entry, "PDF_ALPHAS_DOWN");
    Vars.var_PSSCALE_WEIGHT_4_U *= GetBTagWeightForSyst(entry, "PSSCALE_WEIGHT_4_UP");
    Vars.var_PSSCALE_WEIGHT_4_D *= GetBTagWeightForSyst(entry, "PSSCALE_WEIGHT_4_DOWN");
    Vars.var_PSSCALE_WEIGHT_5_U *= GetBTagWeightForSyst(entry, "PSSCALE_WEIGHT_5_UP");
    Vars.var_PSSCALE_WEIGHT_5_D *= GetBTagWeightForSyst(entry, "PSSCALE_WEIGHT_5_DOWN");
  }
}

void TreePlain::FillPUVarBranches(const AnalysisBaseEntry* entry)
{
  Vars.var_PU_U = entry->pileupW[0] / entry->pileupW[1];
  Vars.var_PU_D = entry->pileupW[2] / entry->pileupW[1];
  
  // Use correct BtagEff for each syst
  if (fillRecoVariationWeights_) {
    Vars.var_PU_U *= GetBTagWeightForSyst(entry, "PU_UP");
    Vars.var_PU_D *= GetBTagWeightForSyst(entry, "PU_DOWN");
  }
}

void TreePlain::FillRecoVarBranches(const AnalysisBaseEntry* entry)
{
  Vars.var_TRIG_U = entry->trigSF[0]/ entry->trigSF[1];
  Vars.var_TRIG_D = entry->trigSF[2]/ entry->trigSF[1];
  // Vars.var_TRIG_ETA_U = entry->trigSF[3]/ entry->trigSF[1];
  // Vars.var_TRIG_ETA_D = entry->trigSF[4]/ entry->trigSF[1];

  Vars.var_LEPT_U = entry->lepSF[0]/ entry->lepSF[1];
  Vars.var_LEPT_D = entry->lepSF[2]/ entry->lepSF[1];
    
  Vars.var_ELE_ID_U = entry->lepSF[3]/ entry->lepSF[1];
  Vars.var_ELE_ID_D = entry->lepSF[4]/ entry->lepSF[1];
  Vars.var_ELE_RECO_U = entry->lepSF[5]/ entry->lepSF[1];
  Vars.var_ELE_RECO_D = entry->lepSF[6]/ entry->lepSF[1];

  Vars.var_MUON_ID_U = entry->lepSF[7]/ entry->lepSF[1];
  Vars.var_MUON_ID_D = entry->lepSF[8]/ entry->lepSF[1];
  Vars.var_MUON_ISO_U = entry->lepSF[9]/ entry->lepSF[1];
  Vars.var_MUON_ISO_D = entry->lepSF[10]/ entry->lepSF[1];

  Vars.var_KIN_U = entry->kinRecoSF[0]/ entry->kinRecoSF[1];
  Vars.var_KIN_D = entry->kinRecoSF[2]/ entry->kinRecoSF[1];

  Vars.var_L1PREFIRING_U = entry->l1prefiringW[0]/entry->l1prefiringW[1];
  Vars.var_L1PREFIRING_D = entry->l1prefiringW[2]/entry->l1prefiringW[1];

  Vars.var_BTAG_U  = GetBTagWeightForSyst(entry, "BTAG_UP");
  Vars.var_BTAG_D = GetBTagWeightForSyst(entry, "BTAG_DOWN");
  Vars.var_BTAG_LJET_U = GetBTagWeightForSyst(entry, "BTAG_LJET_UP");
  Vars.var_BTAG_LJET_D = GetBTagWeightForSyst(entry, "BTAG_LJET_DOWN");

  // Use correct BtagEff for each syst
  Vars.var_TRIG_U *= GetBTagWeightForSyst(entry, "TRIG_UP");
  Vars.var_TRIG_D *= GetBTagWeightForSyst(entry, "TRIG_DOWN");
  // Vars.var_LEPT_U *= GetBTagWeightForSyst(entry, "LEPT_UP");
  // Vars.var_LEPT_D *= GetBTagWeightForSyst(entry, "LEPT_DOWN");
  Vars.var_ELE_ID_U *= GetBTagWeightForSyst(entry, "ELE_ID_UP");
  Vars.var_ELE_ID_D *= GetBTagWeightForSyst(entry, "ELE_ID_DOWN");
  Vars.var_ELE_RECO_U *= GetBTagWeightForSyst(entry, "ELE_RECO_UP");
  Vars.var_ELE_RECO_D *= GetBTagWeightForSyst(entry, "ELE_RECO_DOWN");
  Vars.var_MUON_ID_U *= GetBTagWeightForSyst(entry, "MUON_ID_UP");
  Vars.var_MUON_ID_D *= GetBTagWeightForSyst(entry, "MUON_ID_DOWN");
  Vars.var_MUON_ISO_U *= GetBTagWeightForSyst(entry, "MUON_ISO_UP");
  Vars.var_MUON_ISO_D *= GetBTagWeightForSyst(entry, "MUON_ISO_DOWN");
  Vars.var_L1PREFIRING_U *= GetBTagWeightForSyst(entry, "L1PREFIRING_UP");
  Vars.var_L1PREFIRING_D *= GetBTagWeightForSyst(entry, "L1PREFIRING_DOWN");
  Vars.var_KIN_U *= GetBTagWeightForSyst(entry, "KIN_UP");
  Vars.var_KIN_D *= GetBTagWeightForSyst(entry, "KIN_DOWN");

}
void TreePlain::FillLooseKinRecoBranches(const AnalysisBaseEntry* entry)
{
  if(isGen_ && isTTbarSample_ && entry->topGenObjects->valuesSet_){
    std::pair<float, float> mat = CalculateMatchingPtR(entry, entry->looseKinRecoSolution->getLVBBbar());
    Vars.matchedPt = mat.first;
    Vars.matchedR = mat.second;
  }
  Vars.mttKR9 = entry->looseKinRecoSolution->TTbar().M();
  Vars.yttKR9 = entry->looseKinRecoSolution->TTbar().Rapidity();
  Vars.ptttKR9 = entry->looseKinRecoSolution->TTbar().Pt();
  Vars.mwwKR9 = entry->looseKinRecoSolution->getWWmass();
  Vars.minmlbKR9 = entry->looseKinRecoSolution->getMinMlb();
}
void TreePlain::FillJetBranches(const AnalysisBaseEntry* entry)
{
  //float minPtCut = vExtraJetPtMin.at(0);
  std::pair<TLorentzVector, TLorentzVector> leadingAndTrailingExtraJetIso04Pt30, leadingAndTrailingExtraJetIso04Pt40, leadingAndTrailingExtraJetIso04Pt50,
                                            leadingAndTrailingExtraJetIso04Pt75, leadingAndTrailingExtraJetIso04Pt100, leadingAndTrailingExtraJetIso04Pt150;

  std::vector<TLorentzVector> jetsOmit = {};
  const VLV* inputJets = {};

  bool compute_extrajets_vars = false;

  if(fillRecoTree_){
    compute_extrajets_vars = true;
    inputJets = entry->recoObjects->jets_;
    jetsOmit = { Vars.lvLep, Vars.lvLepBar, Vars.lvBot, Vars.lvBotBar };
  }
  if(fillSimpleRecoTree_){
    compute_extrajets_vars = true;

    const std::pair<TLorentzVector, TLorentzVector>& llbar = entry->looseKinRecoSolution->getLVLLbar();
    const std::pair<TLorentzVector, TLorentzVector>& bbbar = entry->looseKinRecoSolution->getLVBBbar();

    inputJets = entry->recoObjects->jets_;
    jetsOmit = { llbar.first, llbar.second, bbbar.first, bbbar.second };
  }
  if(fillGenTree_){
    inputJets = entry->topGenObjects->allGenJets_;
    
    if(isTTbarSample_){
      jetsOmit = { Vars.lvLep, Vars.lvLepBar, Vars.lvBot, Vars.lvBotBar };
      if(entry->topGenObjects->valuesSet_) compute_extrajets_vars = true;
    }
    else{
      // std::vector<TLorentzVector> jetsOmitEmpty = {  };
      jetsOmit = {  };  
      if(entry->topGenObjects->valuesSet_) compute_extrajets_vars = true;
    }
  }


  if (compute_extrajets_vars){
    TreePlain::Compute_extraJets_observables(Vars.njdefIso04Pt30,Vars.njinetattbarIso04Pt30,leadingAndTrailingExtraJetIso04Pt30,
                                             inputJets, jetsOmit,Vars.lvTop.Eta(),Vars.lvTopBar.Eta(),
                                             vExtraJetPtMin.at(0), extraJetEtaAbsMax, extraJetDeltaRMin);
    TreePlain::Compute_extraJets_observables(Vars.njdefIso04Pt40,Vars.njinetattbarIso04Pt40,leadingAndTrailingExtraJetIso04Pt40,
                                            inputJets, jetsOmit,Vars.lvTop.Eta(),Vars.lvTopBar.Eta(),
                                            vExtraJetPtMin.at(1), extraJetEtaAbsMax, extraJetDeltaRMin);
    TreePlain::Compute_extraJets_observables(Vars.njdefIso04Pt50,Vars.njinetattbarIso04Pt50,leadingAndTrailingExtraJetIso04Pt50,
                                            inputJets, jetsOmit,Vars.lvTop.Eta(),Vars.lvTopBar.Eta(),
                                            vExtraJetPtMin.at(2), extraJetEtaAbsMax, extraJetDeltaRMin);
    TreePlain::Compute_extraJets_observables(Vars.njdefIso04Pt75,Vars.njinetattbarIso04Pt75,leadingAndTrailingExtraJetIso04Pt75,
                                            inputJets, jetsOmit,Vars.lvTop.Eta(),Vars.lvTopBar.Eta(),
                                            vExtraJetPtMin.at(3), extraJetEtaAbsMax, extraJetDeltaRMin);
    TreePlain::Compute_extraJets_observables(Vars.njdefIso04Pt100,Vars.njinetattbarIso04Pt100,leadingAndTrailingExtraJetIso04Pt100,
                                            inputJets, jetsOmit,Vars.lvTop.Eta(),Vars.lvTopBar.Eta(),
                                            vExtraJetPtMin.at(4), extraJetEtaAbsMax, extraJetDeltaRMin);
    TreePlain::Compute_extraJets_observables(Vars.njdefIso04Pt150,Vars.njinetattbarIso04Pt150,leadingAndTrailingExtraJetIso04Pt150,
                                            inputJets, jetsOmit,Vars.lvTop.Eta(),Vars.lvTopBar.Eta(),
                                            vExtraJetPtMin.at(5), extraJetEtaAbsMax, extraJetDeltaRMin);

    // nj_default_cut and Vars.extraJets_j1_j2 must come from the same jet cuts
    Vars.extraJets_j1_j2 = leadingAndTrailingExtraJetIso04Pt30;
    int nj_default_cut = Vars.njdefIso04Pt30;

    Vars.ptej1 = Vars.extraJets_j1_j2.first.Pt();
    Vars.etaej1 = Vars.extraJets_j1_j2.first.Eta();
    Vars.ptej2 = Vars.extraJets_j1_j2.second.Pt();
    Vars.etaej2 = Vars.extraJets_j1_j2.second.Eta();

    TLorentzVector ttbar_solution = TLorentzVector(0.,0.,0.,0.);

    if (fillSimpleRecoTree_) ttbar_solution = common::LVtoTLV(entry->looseKinRecoSolution->TTbar());
    else ttbar_solution = Vars.lvTop + Vars.lvTopBar;

    //no extra jet case
    if (nj_default_cut < 1) {
        Vars.ptej1 = Vars.etaej1 = Vars.ptej2 = Vars.etaej2 = -999.;
        //std::cout<< "not filling extra jets variables" << std::endl;
        //Vars.lvTTbarPlusLeadingExtraJet = ttbar_solution;
    }
    //at least 1 extra jet case
    else {
        if (nj_default_cut < 2) Vars.ptej2 = Vars.etaej2 = -999.;

        Vars.lvTTbarPlusLeadingExtraJet = ttbar_solution + Vars.extraJets_j1_j2.first;
        Vars.mttplusej1 = Vars.lvTTbarPlusLeadingExtraJet.M();
        Vars.ptttplusej1 = Vars.lvTTbarPlusLeadingExtraJet.Pt();
        Vars.r_mttplusej1_mtt = Vars.lvTTbarPlusLeadingExtraJet.M()/ttbar_solution.M();
        Vars.r_ptej1_mtt = Vars.ptej1/ttbar_solution.M();
        if (!fillSimpleRecoTree_) {
          Vars.mean_ptt_pttbar_ptej1 = (Vars.lvTop.Pt() + Vars.lvTopBar.Pt() + Vars.ptej1)/3.;
          Vars.detatej1_detattbar = (Vars.lvTop.Eta() - Vars.extraJets_j1_j2.first.Eta())/(Vars.lvTop.Eta() - Vars.lvTopBar.Eta());
        }
    }

    // std::cout << Vars.ptej1 << " " << Vars.etaej1 << " " << Vars.ptej2 << " " << Vars.etaej2 <<std::endl;
  }

}

void TreePlain::FillPseudoBranches(const AnalysisBaseEntry* entry)
{
  //pseudo particles
  TLorentzVector lvPseudoTop = TLorentzVector(0.,0.,0.,0.), lvPseudoTopBar = TLorentzVector(0.,0.,0.,0.);
  TLorentzVector lvPseudoLep = TLorentzVector(0.,0.,0.,0.), lvPseudoLepBar = TLorentzVector(0.,0.,0.,0.);
  TLorentzVector lvPseudoBot = TLorentzVector(0.,0.,0.,0.), lvPseudoBotBar = TLorentzVector(0.,0.,0.,0.);
  TLorentzVector lvPseudoNu = TLorentzVector(0.,0.,0.,0.), lvPseudoNuBar = TLorentzVector(0.,0.,0.,0.);
  //leading and trailing top, lepton, b-jet
  LV LeadPseudoTop = LV(0.,0.,0.,0.), NLeadPseudoTop = LV(0.,0.,0.,0.);
  TLorentzVector lvLeadPseudoTop = TLorentzVector(0.,0.,0.,0.), lvNLeadPseudoTop = TLorentzVector(0.,0.,0.,0.);
  LV LeadPseudoLep = LV(0.,0.,0.,0.), NLeadPseudoLep = LV(0.,0.,0.,0.);
  TLorentzVector lvLeadPseudoLep = TLorentzVector(0.,0.,0.,0.), lvNLeadPseudoLep = TLorentzVector(0.,0.,0.,0.);
  LV LeadPseudoBot = LV(0.,0.,0.,0.), NLeadPseudoBot = LV(0.,0.,0.,0.);
  TLorentzVector lvLeadPseudoBot = TLorentzVector(0.,0.,0.,0.), lvNLeadPseudoBot = TLorentzVector(0.,0.,0.,0.);
  // create top/antitop quark in the ttbar rest frame
  LV PseudoTTbarSystem = LV(0.,0.,0.,0.);
  LV PseudoTopTTRestFrame = LV(0.,0.,0.,0.), PseudoTopBarTTRestFrame = LV(0.,0.,0.,0.);
  TLorentzVector lvPseudoTopTTRestFrame = TLorentzVector(0.,0.,0.,0.), lvPseudoTopBarTTRestFrame = TLorentzVector(0.,0.,0.,0.);

  if(entry->topPseudoObjects->valuesSet_){
    lvPseudoTop = common::LVtoTLV(*entry->topPseudoObjects->PseudoTop_);
    lvPseudoTopBar = common::LVtoTLV(*entry->topPseudoObjects->PseudoAntiTop_);
    lvPseudoLep = common::LVtoTLV(*entry->topPseudoObjects->PseudoLepton_);
    lvPseudoLepBar = common::LVtoTLV(*entry->topPseudoObjects->PseudoAntiLepton_);
    lvPseudoBot = common::LVtoTLV(*entry->topPseudoObjects->PseudoBJet_);
    lvPseudoBotBar = common::LVtoTLV(*entry->topPseudoObjects->PseudoAntiBJet_);
    lvPseudoNu = common::LVtoTLV(*entry->topPseudoObjects->PseudoNeutrino_);
    lvPseudoNuBar = common::LVtoTLV(*entry->topPseudoObjects->PseudoAntiNeutrino_);

    //leading and trailing top, lepton and b-jets
    common::orderLV(LeadPseudoTop, NLeadPseudoTop, *entry->topPseudoObjects->PseudoTop_, *entry->topPseudoObjects->PseudoAntiTop_, common::LVpt);
    lvLeadPseudoTop = common::LVtoTLV(LeadPseudoTop);
    lvNLeadPseudoTop = common::LVtoTLV(NLeadPseudoTop);

    common::orderLV(LeadPseudoLep, NLeadPseudoLep, *entry->topPseudoObjects->PseudoLepton_, *entry->topPseudoObjects->PseudoAntiLepton_, common::LVpt);
    lvLeadPseudoLep = common::LVtoTLV(LeadPseudoLep);
    lvNLeadPseudoLep = common::LVtoTLV(NLeadPseudoLep);

    common::orderLV(LeadPseudoBot, NLeadPseudoBot, *entry->topPseudoObjects->PseudoBJet_, *entry->topPseudoObjects->PseudoAntiBJet_, common::LVpt);
    lvLeadPseudoBot = common::LVtoTLV(LeadPseudoBot);
    lvNLeadPseudoBot = common::LVtoTLV(NLeadPseudoBot);

    // create top/antitop quark in the ttbar rest frame
    PseudoTTbarSystem = *entry->topPseudoObjects->PseudoTop_ + *entry->topPseudoObjects->PseudoAntiTop_;
    ROOT::Math::Boost CoMBoostGenTtbar (PseudoTTbarSystem.BoostToCM());
    PseudoTopTTRestFrame = CoMBoostGenTtbar(*entry->topPseudoObjects->PseudoTop_);
    PseudoTopBarTTRestFrame = CoMBoostGenTtbar(*entry->topPseudoObjects->PseudoAntiTop_);
    lvPseudoTopTTRestFrame = common::LVtoTLV(PseudoTopTTRestFrame);
    lvPseudoTopBarTTRestFrame = common::LVtoTLV(PseudoTopBarTTRestFrame);

    // Access object selections from config
    const AnalysisConfig::Selections& selections =  entry->analysisConfig->selections();

    // Use utilities without namespaces

    using namespace common;
    using ROOT::Math::VectorUtil::DeltaPhi;
    using ROOT::Math::VectorUtil::DeltaR;
    using ROOT::Math::VectorUtil::Angle;

    TLorentzVector lvPseudoTTbarSystem = lvPseudoTop + lvPseudoTopBar;
    TLorentzVector lvPseudodilepton = lvPseudoLep + lvPseudoLepBar;
    TLorentzVector lvPseudobbar = lvPseudoBot + lvPseudoBotBar;
    TLorentzVector lvPseudoNuNubar = lvPseudoNu + lvPseudoNuBar;

    TVector2 v2MET = TVector2(lvPseudoNuNubar.Px(), lvPseudoNuNubar.Py());
    TLorentzVector lvPseudoMET(v2MET.Px(), v2MET.Py(), 0.0, v2MET.Mod());

    if(fillGenTree_) {
      Vars.visptt = lvPseudoTop.Pt();
      Vars.visptat = lvPseudoTopBar.Pt();
      Vars.vispttLead = lvLeadPseudoTop.Pt();
      Vars.vispttNLead = lvNLeadPseudoTop.Pt();
      Vars.vispttTTRestFrame = lvPseudoTopTTRestFrame.Pt();
      Vars.visptatTTRestFrame = lvPseudoTopBarTTRestFrame.Pt();
      Vars.visyt = lvPseudoTop.Rapidity();
      Vars.visyat = lvPseudoTopBar.Rapidity();
      Vars.visytLead = lvLeadPseudoTop.Rapidity();
      Vars.visytNLead = lvNLeadPseudoTop.Rapidity();
      Vars.vismtt = lvPseudoTTbarSystem.M();
      Vars.visytt = lvPseudoTTbarSystem.Rapidity();
      Vars.vispttt = lvPseudoTTbarSystem.Pt();
      Vars.visdphitt = TMath::Abs(lvPseudoTop.DeltaPhi(lvPseudoTopBar));
      Vars.visdytt = TMath::Abs(lvPseudoTop.Rapidity()) - TMath::Abs(lvPseudoTopBar.Rapidity());
      Vars.visdetatt = TMath::Abs(lvPseudoTop.PseudoRapidity() - lvPseudoTopBar.PseudoRapidity());
      Vars.visrecoPartonMomFraction = (lvPseudoTop.Energy() - lvPseudoTop.Pz() + lvPseudoTopBar.Energy() - lvPseudoTopBar.Pz()) / (2 * 6500);
      Vars.visrecoAntipartonMomFraction = (lvPseudoTop.Energy() + lvPseudoTop.Pz() + lvPseudoTopBar.Energy() + lvPseudoTopBar.Pz()) / (2 * 6500);

      Vars.visptl = lvPseudoLep.Pt();
      Vars.visptal = lvPseudoLepBar.Pt();
      Vars.visptlLead = lvLeadPseudoLep.Pt();
      Vars.visptlNLead = lvNLeadPseudoLep.Pt();
      Vars.visetal = lvPseudoLep.PseudoRapidity();
      Vars.visetaal = lvPseudoLepBar.PseudoRapidity();
      Vars.visetalLead = lvLeadPseudoLep.PseudoRapidity();
      Vars.visetalNLead = lvNLeadPseudoLep.PseudoRapidity();
      Vars.visdetall = TMath::Abs(lvPseudoLep.PseudoRapidity()) - TMath::Abs(lvPseudoLepBar.PseudoRapidity());
      Vars.visdphill = TMath::Abs(lvPseudoLep.DeltaPhi(lvPseudoLepBar));;
      Vars.visptll = lvPseudodilepton.Pt();
      Vars.visyll = lvPseudodilepton.Rapidity();
      Vars.vismll = lvPseudodilepton.M();
      Vars.visptllbb = (lvPseudodilepton+lvPseudobbar).Pt();
      Vars.visyllbb = (lvPseudodilepton+lvPseudobbar).Rapidity();
      Vars.vismllbb = (lvPseudodilepton+lvPseudobbar).M();
      Vars.visptllbbmet = (lvPseudodilepton+lvPseudobbar+lvPseudoMET).Pt();
      Vars.visyllbbmet = (lvPseudodilepton+lvPseudobbar+lvPseudoMET).Rapidity();
      Vars.vismllbbmet = (lvPseudodilepton+lvPseudobbar+lvPseudoMET).M();
      Vars.visptbLead = lvLeadPseudoBot.Pt();
      Vars.visptbNLead = lvNLeadPseudoBot.Pt();
      Vars.visetabLead = lvLeadPseudoBot.PseudoRapidity();
      Vars.visetabNLead = lvNLeadPseudoBot.PseudoRapidity();
      Vars.visptbb = lvPseudobbar.Pt();
      Vars.visybb = lvPseudobbar.Rapidity();
      Vars.vismbb = lvPseudobbar.M();
      Vars.vismww = (lvPseudoTTbarSystem - lvPseudoBot - lvPseudoBotBar).M();
      Vars.visminmlb = std::min((lvPseudoLep + lvPseudoBotBar).M(), (lvPseudoLepBar + lvPseudoBot).M());

    }

    if ( lvPseudoLep.Pt() > selections.leptonPtCut_ && lvPseudoLepBar.Pt() > selections.leptonPtCut_ &&
	 std::fabs( lvPseudoLep.PseudoRapidity() ) < selections.leptonEtaCut_ && std::fabs ( lvPseudoLepBar.PseudoRapidity() ) < selections.leptonEtaCut_ &&
	 lvPseudoBot.Pt() > selections.jetPtCut_ && std::fabs ( lvPseudoBot.PseudoRapidity() ) < selections.jetEtaCut_ &&
	 lvPseudoBotBar.Pt() > selections.jetPtCut_ && std::fabs ( lvPseudoBotBar.PseudoRapidity() ) < selections.jetEtaCut_ &&
	 std::fabs(DeltaR(lvPseudoLep, lvPseudoBot)) > selections.genDeltaRLeptonJetCut_  &&
	 std::fabs(DeltaR(lvPseudoLep, lvPseudoBotBar)) > selections.genDeltaRLeptonJetCut_  &&
	 std::fabs(DeltaR(lvPseudoLepBar, lvPseudoBot)) > selections.genDeltaRLeptonJetCut_  &&
	 std::fabs(DeltaR(lvPseudoLepBar, lvPseudoBotBar)) > selections.genDeltaRLeptonJetCut_ &&
	 lvPseudodilepton.M() > 20.0 ) {


      Vars.inFiducialPhaseSpace = true;

      }
    else {

      Vars.inFiducialPhaseSpace = false;
    }
  }
}

int TreePlain::RecoTTBarDecayFromPdgId(const int pdgid1, const int pdgid2)
{
  int decay = -1;

  if(TMath::Abs(pdgid1) == 11)
    decay = 2;
  else if(TMath::Abs(pdgid1) == 13)
    decay = 3;
  else
    decay = -1;

  decay *= 10;

  if(TMath::Abs(pdgid2) == 11)
    decay += 2;
  else if(TMath::Abs(pdgid2) == 13)
    decay += 3;
  else
    decay += -1;

  return decay;
}

template<class T>
void TreePlain::AddBranchVerbosely(const char* name, T* ptr)
{
  //printf("Adding branch %20s ... ", name);
  //auto ret = this->Branch(name, ptr);
  auto ret = TTreePtr->Branch(name, ptr);
  if(ret)
    ;//printf("success\n");
  else
  {
    printf("Adding branch %20s ... ", name);printf("failed\n");
    throw std::logic_error(TString::Format("Error while adding branch %s\n", name).Data());
  }
}

std::pair<bool, std::pair<float, float> > TreePlain::IsMatched(const TLorentzVector& lv1, const TLorentzVector& lv2, const float maxDeltaPt, const float maxDeltaR)
{
  float deltaPt = TMath::Abs(lv1.Pt() - lv2.Pt());
  float deltaR = ROOT::Math::VectorUtil::DeltaR(lv1, lv2);
  bool matched = (maxDeltaPt < 0.0 || (deltaPt < maxDeltaPt)) && (maxDeltaR < 0.0 || (deltaR < maxDeltaR));
  return std::pair<bool, std::pair<float, float> >(matched, std::pair<float, float>(deltaPt, deltaR));
}

// std::pair<float, float> TreePlain::CalculateMatchingPtR(const ZTreeNtupleVars& inVars, const std::pair<TLorentzVector, TLorentzVector>& bbbar)
std::pair<float, float> TreePlain::CalculateMatchingPtR(const AnalysisBaseEntry* entry, const std::pair<TLorentzVector, TLorentzVector>& bbbar)
{
  const auto& mat11 = IsMatched(common::LVtoTLV(*entry->topGenObjects->GenB_), bbbar.first);
  const auto& mat12 = IsMatched(common::LVtoTLV(*entry->topGenObjects->GenAntiB_), bbbar.first);
  const auto& mat21 = IsMatched(common::LVtoTLV(*entry->topGenObjects->GenB_), bbbar.second);
  const auto& mat22 = IsMatched(common::LVtoTLV(*entry->topGenObjects->GenAntiB_), bbbar.second);
  const float matMinPt1 = std::min(mat11.second.first / entry->topGenObjects->GenB_->Pt(), mat12.second.first / entry->topGenObjects->GenAntiB_->Pt());
  const float matMinPt2 = std::min(mat21.second.first / entry->topGenObjects->GenB_->Pt(), mat22.second.first / entry->topGenObjects->GenAntiB_->Pt());
  const float matMinPt = std::max(matMinPt1, matMinPt2);
  const float matMinR1 = std::min(mat11.second.second, mat12.second.second);
  const float matMinR2 = std::min(mat21.second.second, mat22.second.second);
  const float matMinR = std::max(matMinR1, matMinR2);
  return std::pair<float, float>(matMinPt, matMinR);
}

std::pair<float, float> TreePlain::CalculateMatchingPtRFullKR(const AnalysisBaseEntry* entry, const TLorentzVector& b, const TLorentzVector& bbar)
{
  const auto& mat1 = IsMatched(common::LVtoTLV(*entry->topGenObjects->GenB_), b);
  const auto& mat2 = IsMatched(common::LVtoTLV(*entry->topGenObjects->GenAntiB_), bbar);
  const float matMinPt = std::min(mat1.second.first / entry->topGenObjects->GenB_->Pt(), mat2.second.first / entry->topGenObjects->GenAntiB_->Pt());
  const float matMinR = std::min(mat1.second.second, mat2.second.second);
  return std::pair<float, float>(matMinPt, matMinR);
}



int TreePlain::SelectExtraGenJets(const VLV* jets, const std::vector<TLorentzVector> jetsOmit,
                                  float& extraJetPtMin, float& extraJetEtaAbsMax, float& extraJetDeltaRMin)
// int TreePlain::SelectExtraGenJets(const VLV* jets, const std::vector<TLorentzVector> jetsOmit,
//                                   float& extraJetPtMin, float& extraJetEtaAbsMax, float& extraJetDeltaRMin, int& extraJetNLead)
// int TreePlain::SelectExtraGenJets(const VLV* jets, const std::vector<TLorentzVector> jetsOmit,
//                                   float& extraJetPtMin, float& extraJetEtaAbsMax, float& extraJetDeltaRMin, int& extraJetNLead,
//                                   const TLorentzVector& ttbar, float* outMttj, float* outYttj, float* outRhoj, float* maxptjPtr)
{
  //printf("config: pt = %f iso = %f lead = %d\n", config.extraJetPtMin, config.extraJetDeltaRMin, config.extraJetNLead);
  int njet = 0;
  // TLorentzVector sumLeadJets;
  // int njetLead = 0;
  for(unsigned int i = 0; i < jets->size(); i++)
  {
    // make sure jets are sorted in pT
    assert(i == 0 || jets->at(i).Pt() <= jets->at(i - 1).Pt());

    const auto& jet = jets->at(i);
    if(jet.Pt() > extraJetPtMin && TMath::Abs(jet.Eta()) < extraJetEtaAbsMax)
    {
      // remove jets close to leptons and b/bbar
      if(extraJetDeltaRMin > 0)
      {
        bool flagBreak = false;
        for(unsigned int jj = 0; jj < jetsOmit.size(); jj++)
        {
          if(ROOT::Math::VectorUtil::DeltaR(jet, jetsOmit[jj]) < extraJetDeltaRMin)
          {
            flagBreak = true;
            break;
          }
        }
        if(flagBreak)
          continue;
      }
      // if(extraJetNLead < 0 || njet < extraJetNLead)
      // {
      //   // sumLeadJets += common::LVtoTLV(jet);
      //   // njetLead++;
      // }
      njet++;
    }
  }
  //sumLeadJets.Print();
  // TLorentzVector ttj = ttbar + sumLeadJets;
  // float yttj = ttj.Rapidity();
  // float mttj = ttj.M();
  //
  // if(outYttj)
  //   *outYttj = yttj;
  // if(outMttj)
  //   *outMttj = mttj;
  // if(outRhoj)
  // {
  //   float rhoj = Rhoj(mttj, njetLead, extraJetPtMin);
  //   *outRhoj = rhoj;
  // }
  // if(maxptjPtr)
  //   *maxptjPtr = (jets->size() == 0) ? 0.0 : jets->at(0).Pt();

  return njet;
}

//Select at the same time number off extra jets and the number of extra jets inside the ttbar (in rapidity)
void TreePlain::Select_N_extraJets_and_N_j_inside_ttbar(int &njet, int &njet_inside, const VLV* jets, const std::vector<TLorentzVector> jetsOmit, float eta_t, float eta_at,
                                  float& extraJetPtMin, float& extraJetEtaAbsMax, float& extraJetDeltaRMin)
{
  //printf("config: pt = %f iso = %f lead = %d\n", config.extraJetPtMin, config.extraJetDeltaRMin, config.extraJetNLead);
  njet = njet_inside = 0;
  // TLorentzVector sumLeadJets;
  // int njetLead = 0;
  for(unsigned int i = 0; i < jets->size(); i++)
  {
    // make sure jets are sorted in pT
    assert(i == 0 || jets->at(i).Pt() <= jets->at(i - 1).Pt());

    const auto& jet = jets->at(i);
    if(jet.Pt() > extraJetPtMin && TMath::Abs(jet.Eta()) < extraJetEtaAbsMax)
    {
      // remove jets close to leptons and b/bbar
      if(extraJetDeltaRMin > 0)
      {
        bool flagBreak = false;
        for(unsigned int jj = 0; jj < jetsOmit.size(); jj++)
        {
          if(ROOT::Math::VectorUtil::DeltaR(jet, jetsOmit[jj]) < extraJetDeltaRMin)
          {
            flagBreak = true;
            break;
          }
        }
        if(flagBreak)
          continue;
      }
      njet++;
      if ((jet.Eta() - eta_t)*(jet.Eta() - eta_at) < 0) njet_inside++;
    }
  }
  if (njet == 0) njet_inside=-10;
}

void TreePlain::Compute_extraJets_observables(int &njet, int &njet_inside, std::pair<TLorentzVector, TLorentzVector> & leadingAndTrailingExtraJet, 
                                              const VLV* jets, const std::vector<TLorentzVector> jetsOmit, float eta_t, float eta_at,
                                              float& extraJetPtMin, float& extraJetEtaAbsMax, float& extraJetDeltaRMin)
{
  //printf("config: pt = %f iso = %f lead = %d\n", config.extraJetPtMin, config.extraJetDeltaRMin, config.extraJetNLead);
  njet = njet_inside = 0;
  TLorentzVector leadingJet;
  TLorentzVector trailingJet;
  bool leading_jet_found = false;
  bool trailing_jet_found = false;

  // TLorentzVector sumLeadJets;
  // int njetLead = 0;
  for(unsigned int i = 0; i < jets->size(); i++)
  {
    // make sure jets are sorted in pT
    assert(i == 0 || jets->at(i).Pt() <= jets->at(i - 1).Pt());

    const auto& jet = jets->at(i);
    if(jet.Pt() > extraJetPtMin && TMath::Abs(jet.Eta()) < extraJetEtaAbsMax){
      // remove jets close to leptons and b/bbar
      if(extraJetDeltaRMin > 0)
      {
        bool flagBreak = false;
        for(unsigned int jj = 0; jj < jetsOmit.size(); jj++)
        {
          if(ROOT::Math::VectorUtil::DeltaR(jet, jetsOmit[jj]) < extraJetDeltaRMin)
          {
            flagBreak = true;
            break;
          }
        }
        if(flagBreak)
          continue;
      }
      njet++;
      if ((jet.Eta() - eta_t)*(jet.Eta() - eta_at) < 0) njet_inside++;

      if (!leading_jet_found) {
                leadingJet = common::LVtoTLV(jet);
                leading_jet_found = true;
                // std::cout << "Found leading jet:" << leadingJet.Pt() << std::endl;
      }
      else if (!trailing_jet_found){
                trailingJet = common::LVtoTLV(jet);
                trailing_jet_found = true;
                // std::cout << "Found sud-leading jet:" << trailingJet.Pt() << std::endl;
      }
    }
  }
  if (njet < 1) {
    njet_inside=-10;
  }

  if (!leading_jet_found){
      leadingJet.SetPtEtaPhiE(-999., -999., -999., -999.);
  }
  if (!trailing_jet_found){
      trailingJet.SetPtEtaPhiE(-999., -999., -999., -999.);
  }

  leadingAndTrailingExtraJet = std::pair<TLorentzVector, TLorentzVector> (leadingJet,trailingJet);

}

float TreePlain::GetBTagWeightForSyst(const AnalysisBaseEntry* entry, const TString& systName){
    float nominalBtagSf = entry->nominalBtagSF;
    bool printDebug = false;

    float correctedWeight = weightBtagSF_Var(*entry->jetIndices, *entry->jets, *entry->jetFlavour, *entry->btagDiscriminators, systName)/nominalBtagSf;

    if (printDebug) std::cout << "GetBTagWeightForSyst returned: " << correctedWeight << " for syst: " << systName << std::endl 
                              << "        nominalBtagSf: " << nominalBtagSf << ";  systBtagSf: " << weightBtagSF_Var(*entry->jetIndices, *entry->jets, *entry->jetFlavour, *entry->btagDiscriminators, systName)
                              << std::endl;
    return correctedWeight;
}


void TreePlain::SetBtagScaleFactors_Var_AllSyst(const std::map<std::tuple<Systematic::Type, Systematic::Variation, int>, BtagScaleFactors*>* mapAllSystBTagScaleFactors)
{
    mapAllSystBTagScaleFactors_ = mapAllSystBTagScaleFactors;
}

float TreePlain::weightBtagSF_Var(const std::vector<int>& jetIndices,
                                  const VLV& jets, const std::vector<int>& jetFlavour,
                                  const std::vector<float>& btagDiscriminators, const TString& sysKey)const
{
    if(!isGen_) return 1.;
    if(!mapAllSystBTagScaleFactors_) return 1.;
    Systematic::Systematic temporarySystematic = Systematic::Systematic(sysKey);
    std::tuple<Systematic::Type, Systematic::Variation, int> tempTuple{temporarySystematic.type(), temporarySystematic.variation(), temporarySystematic.variationNumber()};
    return mapAllSystBTagScaleFactors_->at(tempTuple)->getSF(jetIndices, jets, jetFlavour, btagDiscriminators);
}
/*
std::pair<TLorentzVector, TLorentzVector> TreePlain::SelectLeadingAndTrailingExtraJet(const VLV* jets, const std::vector<TLorentzVector> jetsOmit,
                                                                                        float& extraJetPtMin, float& extraJetEtaAbsMax, float& extraJetDeltaRMin){
    TLorentzVector leadingJet;
    TLorentzVector trailingJet;

    bool leading_jet_found = false;
    bool trailing_jet_found = false;

    // std::cout << "Selecting leading jet" << std::endl;

    // make sure jets are sorted in pT
    for(unsigned int i = 0; i < jets->size(); i++){
        
        if (leading_jet_found && trailing_jet_found) break; 

        assert(i == 0 || jets->at(i).Pt() <= jets->at(i - 1).Pt());

        const auto& jet = jets->at(i);
        // std::cout << jet.Pt() << std::endl;
        if(jet.Pt() > extraJetPtMin && TMath::Abs(jet.Eta()) < extraJetEtaAbsMax){
            // remove jets close to leptons and b/bbar
            if(extraJetDeltaRMin > 0){
                bool flagBreak = false;
                for(unsigned int jj = 0; jj < jetsOmit.size(); jj++){
                  if(ROOT::Math::VectorUtil::DeltaR(jet, jetsOmit[jj]) < extraJetDeltaRMin){
                      flagBreak = true;
                      break;
                  }
                }
                if(flagBreak) continue;
            }
            if (!leading_jet_found) {
                leadingJet = common::LVtoTLV(jet);
                leading_jet_found = true;
                // std::cout << "Found leading jet:" << leadingJet.Pt() << std::endl;
            }
            else if (!trailing_jet_found){
                trailingJet = common::LVtoTLV(jet);
                trailing_jet_found = true;
                // std::cout << "Found sud-leading jet:" << trailingJet.Pt() << std::endl;
            }
            else if (leading_jet_found && trailing_jet_found) break;
        }
    }

    if (!leading_jet_found){
      leadingJet.SetPtEtaPhiE(-999., -999., -999., -999.);
    }
    if (!trailing_jet_found){
      trailingJet.SetPtEtaPhiE(-999., -999., -999., -999.);
    }
    //std::cout << "Leading Jet: " << leadingJet.Pt() << "; SubLeading Jet: " << trailingJet.Pt() << "; Number of Jets: " << jets->size() << std::endl;

    return std::pair<TLorentzVector, TLorentzVector> (leadingJet,trailingJet);
}*/

// void TreePlain::FillReco(const AnalysisBaseEntry* entry)
// {
//   // no access to gen objects
// #define topGenObjects NOACCESS
// #define genObjectIndices NOACCESS
// #define genLevelWeights NOACCESS
//
//   // ttbar
//   Vars.lvTop = common::LVtoTLV(entry->kinematicReconstructionSolutions->solution().top());
//   Vars.lvTopBar = common::LVtoTLV(entry->kinematicReconstructionSolutions->solution().antiTop());
//   Vars.lvLep = common::LVtoTLV(entry->kinematicReconstructionSolutions->solution().lepton());
//   Vars.lvLepBar = common::LVtoTLV(entry->kinematicReconstructionSolutions->solution().antiLepton());
//   Vars.lvBot = common::LVtoTLV(entry->kinematicReconstructionSolutions->solution().bjet());
//   Vars.lvBotBar = common::LVtoTLV(entry->kinematicReconstructionSolutions->solution().antiBjet());
//   //Vars.indBot = entry->kinematicReconstructionSolutions->solution().bjetIndex();
//   //Vars.indBotBar = entry->kinematicReconstructionSolutions->solution().antiBjetIndex();
//   //Vars.idLep1 = entry->recoObjects->lepPdgId_->at(0);
//   //Vars.idLep2 = entry->recoObjects->lepPdgId_->at(1);
//   // Vars.TTBarDecay = this->RecoTTBarDecayFromPdgId(entry->recoObjects->lepPdgId_->at(0), entry->recoObjects->lepPdgId_->at(1));
//   // TODO sanity check of ttbar decay kinematics ?
//
//   // jets
//   // Vars.vlvJets.clear();
//   // const int njets = entry->recoObjectIndices->jetIndices_.size();
//   // Vars.vlvJets.resize(njets);
//   // for(int i = 0; i < njets; i++)
//   // {
//   //   const int index = entry->recoObjectIndices->jetIndices_[i];
//   //   Vars.vlvJets[i] = common::LVtoTLV(entry->recoObjects->jets_->at(index));
//   // }
//
//   // b jets
//   // Vars.vlvBJets.clear();
//   // const int nbjets = entry->recoObjectIndices->bjetIndices_.size();
//   // Vars.vlvBJets.resize(nbjets);
//   // for(int i = 0; i < nbjets; i++)
//   // {
//   //   const int index = entry->recoObjectIndices->bjetIndices_[i];
//   //   Vars.vlvBJets[i] = common::LVtoTLV(entry->recoObjects->jets_->at(index));
//   // }
//
//   // MET
//   // Vars.v2MET = TVector2(entry->recoObjects->met_->Px(), entry->recoObjects->met_->Py());
//   //assert(TMath::Abs(entry->recoObjects->met_->Pz()) < 1e-4);
//   //assert(TMath::Abs(entry->recoObjects->met_->E()) < 1e-4);
//   //assert("assert in here?");
//
//   // primary vertices
//   // Vars.nPrimVtx = entry->recoObjects->vertMulti_;
//
//   //printf("rec entry: %d  FlagEnablePDFWeights = %d\n", *entry->eventCounter, FlagEnablePDFWeights);
//
//   // cross section variables
//   TLorentzVector t = Vars.lvTop;
//   TLorentzVector tbar = Vars.lvTopBar;
//   this->CalculateTtbarXsecVars(Vars, t, tbar);
//   //outVarsPtr->njet = jets.size() - 2;
//   // Vars.njet = njets;
//
// #undef topGenObjects
// #undef genObjectIndices
// #undef genLevelWeights
// }

// void TreePlain::FillGen(const AnalysisBaseEntry* entry, const float deltaRLepJet)
// {
//   // no access to reco objects
// #define kinematicReconstructionSolutions NOACCESS
// #define recoObjectIndices NOACCESS
// #define recoObjects NOACCESS
//
//   Vars.weight = entry->genLevelWeights->trueLevelWeight_;
//   // Vars.TTBarDecay = entry->topGenObjects->decayMode_;
//
//   // ttbar
//   Vars.lvTop = common::LVtoTLV(*entry->topGenObjects->GenTop_);
//   Vars.lvTopBar = common::LVtoTLV(*entry->topGenObjects->GenAntiTop_);
//   Vars.lvLep = common::LVtoTLV(*entry->topGenObjects->GenLepton_);
//   Vars.lvLepBar = common::LVtoTLV(*entry->topGenObjects->GenAntiLepton_);
//   Vars.lvBot = common::LVtoTLV(*entry->topGenObjects->GenB_);
//   Vars.lvBotBar = common::LVtoTLV(*entry->topGenObjects->GenAntiB_);
//   //Vars.indBot = entry->genObjectIndices->genBjetFromTopIndex_;
//   //Vars.indBotBar = entry->genObjectIndices->genAntiBjetFromTopIndex_;
//   // Vars.TTBarDecay = entry->topGenObjects->decayMode_;
//
//   // TODO sanity check of ttbar decay kinematics ?
//
//   // jets
//   // Vars.vlvJets.clear();
//   // const int njets = entry->genObjectIndices->genVisJetIndices_.size();
//   // Vars.vlvJets.resize(njets);
//   // TODO how to access jet flavour?
//   /*if(this->FlagStoreExtraVars)
//   {
//     Vars.vdJetsBDiscr.clear();
//     Vars.vdJetsBDiscr.resize(njets);
//   }*/
//   // int nj = 0; // number of extra jets separated from leptons
//   // for(int i = 0; i < njets; i++)
//   // {
//   //   const int index = entry->genObjectIndices->genVisJetIndices_[i];
//   //   Vars.vlvJets[i] = common::LVtoTLV(entry->topGenObjects->allGenJets_->at(index));
//   //   // TODO how to access jet flavour?
//   //   //if(this->FlagStoreExtraVars)
//   //   //  Vars.vdJetsBDiscr[i] = ((entry->Vars.jetHadronFlavour)[i] == 5) ? 1.0 : 0.0;
//   //
//   //   // remove jets close to leptons
//   //   if(deltaRLepJet > 0)
//   //   {
//   //     if(ROOT::Math::VectorUtil::DeltaR(Vars.vlvJets[i], Vars.lvLep) < deltaRLepJet)
//   //       continue;
//   //     if(ROOT::Math::VectorUtil::DeltaR(Vars.vlvJets[i], Vars.lvLepBar) < deltaRLepJet)
//   //       continue;
//   //   }
//   //   nj++;
//   // }
//
//   // b jets
//   /*Vars.vlvBJets.clear();
//   const int nbjets = entry->genObjectIndices->genBjetFromTopIndex_.size();
//   Vars.vlvBJets.resize(nbjets);
//   for(int index = 0; index < nbjets; index++)
//     Vars.vlvBJets[index] = common::LVtoTLV(entry->topGenObjects->allGenJets_->at(index));*/
//
//   // MET
//   // Vars.v2MET = TVector2(entry->topGenObjects->GenMet_->Px(), entry->topGenObjects->GenMet_->Py());
//   //assert(TMath::Abs(entry->topGenObjects->GenMet_->Pz()) < 1e-4);
//   //assert(TMath::Abs(entry->topGenObjects->GenMet_->E()) < 1e-4);
//
//   // primary vertices
//   //nPrimVtx = entry->recoObjects->vertMulti_;
//
//   //printf("gen entry: %d  FlagEnablePDFWeights = %d\n", *entry->eventCounter, FlagEnablePDFWeights);
//   if(FlagEnablePDFWeights){
//     // get branch entry if necessary
//     if(entry->analysisBase->b_pdfWeights->GetEntryNumber() != *entry->eventCounter)
//       entry->analysisBase->b_pdfWeights->GetEntry(*entry->eventCounter);
//
//     const std::vector<double>& pdfWeightsFloat = *(entry->analysisBase->pdfWeights_);
//     //printf("writing PDFs size = %ld\n", pdfWeightsFloat.size());
//     // std::vector<float>& pdfWeightsFloat = this->Vars.weightPDF;
//     // pdfWeightsFloat.clear();
//     // pdfWeightsFloat.resize(pdfWeightsFloat.size());
//     for(unsigned int i = 0; i < pdfWeightsFloat.size(); i++){
//       // pdfWeightsFloat[i] = pdfWeightsFloat[i];
//       this->Vars.var_PDF[i]=pdfWeightsFloat[i];
//       // outputVars.var_PDF[w] = (*inVars.pdfWeights)[w] / ((w == 101) ? 1.0 : norm);
//     }
//   }
//
//   // cross section variables
//   TLorentzVector t = Vars.lvTop;
//   TLorentzVector tbar = Vars.lvTopBar;
//   this->CalculateTtbarXsecVars(Vars, t, tbar);
//   //Vars.njet = njets;
//   // Vars.njet = nj;
//
// #undef kinematicReconstructionSolutions
// #undef recoObjectIndices
// #undef recoObjects
// }

/*
void TreePlain::CalculateTtbarXsecVars(TreePlainVars& outputVars, const TLorentzVector& t, const TLorentzVector& tbar)
{
  TLorentzVector ttbar = t + tbar;
  outputVars.ptt = t.Pt();
  outputVars.yt = t.Rapidity();
  outputVars.mtt = ttbar.M();
  outputVars.ytt = ttbar.Rapidity();
  outputVars.pttt = ttbar.Pt();
  outputVars.dphitt = TMath::Abs(t.DeltaPhi(tbar));
  outputVars.detatt = TMath::Abs(t.PseudoRapidity() - tbar.PseudoRapidity());
}
*/
