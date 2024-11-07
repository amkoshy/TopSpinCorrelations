#include "TopAnalysisPlain.h"
#include "AnalyzerBaseClass.h"
#include "TreeHandlerBase.h"
// #include "TreePlain.h"
#include "ScaleFactors.h"
// #include "analysisObjectStructs.h"
#include "AnalysisConfig.h"
// #include <TFile.h>
#include <TTree.h>
#include "../../common/include/KinematicReconstructionSolution.h"
#include "../../common/include/KinematicReconstruction.h"
#include "../../common/include/LooseKinRecoSolution.h"
#include "../../common/include/LooseKinReco.h"

#include <TH1.h>

// void TopAnalysisPlain::SetLeptonScaleFactors_Up(const LeptonScaleFactors* const scaleFactors)
// {
//     leptonScaleFactors_Up_ = scaleFactors;
// }
// void TopAnalysisPlain::SetLeptonScaleFactors_Dn(const LeptonScaleFactors* const scaleFactors)
// {
//     leptonScaleFactors_Dn_ = scaleFactors;
// }
void TopAnalysisPlain::SetTriggerScaleFactors_Up(const TriggerScaleFactors* const scaleFactors)
{
    triggerScaleFactors_Up_ = scaleFactors;
}
void TopAnalysisPlain::SetTriggerScaleFactors_Dn(const TriggerScaleFactors* const scaleFactors)
{
    triggerScaleFactors_Dn_ = scaleFactors;
}
// void TopAnalysisPlain::SetTriggerScaleFactors_EtaUp(const TriggerScaleFactors* const scaleFactors)
// {
//     triggerScaleFactors_EtaUp_ = scaleFactors;
// }
// void TopAnalysisPlain::SetTriggerScaleFactors_EtaDn(const TriggerScaleFactors* const scaleFactors)
// {
//     triggerScaleFactors_EtaDn_ = scaleFactors;
// }

void TopAnalysisPlain::SetBtagScaleFactors_Var_AllSyst(const std::map<std::tuple<Systematic::Type, Systematic::Variation, int>, BtagScaleFactors*>* mapAllSystBTagScaleFactors){
    mapAllSystBTagScaleFactors_ = mapAllSystBTagScaleFactors;
}

  // kinematicReconstruction as dummy to keep same structure
void TopAnalysisPlain::SetKinematicReconstruction_Up(const KinematicReconstruction* const kinematicReconstruction,
                                              const KinematicReconstructionScaleFactors* const kinematicReconstructionScaleFactors)
{
    (void) kinematicReconstruction;
    kinematicReconstructionScaleFactors_Up_ = kinematicReconstructionScaleFactors;
}
void TopAnalysisPlain::SetKinematicReconstruction_Dn(const KinematicReconstruction* const kinematicReconstruction,
                                              const KinematicReconstructionScaleFactors* const kinematicReconstructionScaleFactors)
{
    (void) kinematicReconstruction;
    kinematicReconstructionScaleFactors_Dn_ = kinematicReconstructionScaleFactors;
}

  // LooseKinReco as dummy to keep same structure
void TopAnalysisPlain::SetLooseKinReco_Up(const LooseKinReco* const looseKinReco,
					  const LooseKinRecoScaleFactors* const looseKinRecoScaleFactors)
{
    (void) looseKinReco;
    looseKinRecoScaleFactors_Up_ = looseKinRecoScaleFactors;
}
void TopAnalysisPlain::SetLooseKinReco_Dn(const LooseKinReco* const looseKinReco,
					  const LooseKinRecoScaleFactors* const looseKinRecoScaleFactors)
{
    (void) looseKinReco;
    looseKinRecoScaleFactors_Dn_ = looseKinRecoScaleFactors;
}

// float TopAnalysisPlain::weightLeptonSF_Up(const int leadingLeptonIndex, const int nLeadingLeptonIndex,
//                                     const VLV& allLeptons, const std::vector<int>& lepPdgId)const
// {
//     if(!isMC_) return 1.;
//     if(leadingLeptonIndex<0 || nLeadingLeptonIndex<0) return 1.;
//     if(!leptonScaleFactors_Up_) return 1.;
//     return leptonScaleFactors_Up_->getSFDilepton(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId);
// }
//
//
// float TopAnalysisPlain::weightLeptonSF_Up(const int leadingLeptonIndex, const int nLeadingLeptonIndex,
//                                     const VLV& allLeptons, const std::vector<int>& lepPdgId, const std::vector<float>& lepSCEta)const
// {
//     if(!isMC_) return 1.;
//     if(leadingLeptonIndex<0 || nLeadingLeptonIndex<0) return 1.;
//     if(!leptonScaleFactors_Up_) return 1.;
//     return leptonScaleFactors_Up_->getSFDilepton(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId, lepSCEta);
// }
// float TopAnalysisPlain::weightLeptonSF_Dn(const int leadingLeptonIndex, const int nLeadingLeptonIndex,
//                                     const VLV& allLeptons, const std::vector<int>& lepPdgId)const
// {
//     if(!isMC_) return 1.;
//     if(leadingLeptonIndex<0 || nLeadingLeptonIndex<0) return 1.;
//     if(!leptonScaleFactors_Dn_) return 1.;
//     return leptonScaleFactors_Dn_->getSFDilepton(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId);
// }
//
//
// float TopAnalysisPlain::weightLeptonSF_Dn(const int leadingLeptonIndex, const int nLeadingLeptonIndex,
//                                     const VLV& allLeptons, const std::vector<int>& lepPdgId, const std::vector<float>& lepSCEta)const
// {
//     if(!isMC_) return 1.;
//     if(leadingLeptonIndex<0 || nLeadingLeptonIndex<0) return 1.;
//     if(!leptonScaleFactors_Dn_) return 1.;
//     return leptonScaleFactors_Dn_->getSFDilepton(leadingLeptonIndex, nLeadingLeptonIndex, allLeptons, lepPdgId, lepSCEta);
// }
//
float TopAnalysisPlain::weightTriggerSF_Up(const int leptonXIndex, const int leptonYIndex,
                                     const VLV& allLeptons, const std::vector<int>& lepPdgId)const
{
    if(!isMC_) return 1.;
    if(leptonXIndex<0 || leptonYIndex<0) return 1.;
    if(!triggerScaleFactors_Up_) return 1.;
    return triggerScaleFactors_Up_->getSF(leptonXIndex, leptonYIndex, allLeptons, lepPdgId);
}

float TopAnalysisPlain::weightTriggerSF_Dn(const int leptonXIndex, const int leptonYIndex,
                                     const VLV& allLeptons, const std::vector<int>& lepPdgId)const
{
    if(!isMC_) return 1.;
    if(leptonXIndex<0 || leptonYIndex<0) return 1.;
    if(!triggerScaleFactors_Dn_) return 1.;
    return triggerScaleFactors_Dn_->getSF(leptonXIndex, leptonYIndex, allLeptons, lepPdgId);
}
// float TopAnalysisPlain::weightTriggerSF_EtaUp(const int leptonXIndex, const int leptonYIndex,
//                                      const VLV& allLeptons, const std::vector<int>& lepPdgId)const
// {
//     if(!isMC_) return 1.;
//     if(leptonXIndex<0 || leptonYIndex<0) return 1.;
//     if(!triggerScaleFactors_EtaUp_) return 1.;
//     return triggerScaleFactors_EtaUp_->getSF(leptonXIndex, leptonYIndex, allLeptons, lepPdgId);
// }
// float TopAnalysisPlain::weightTriggerSF_EtaDn(const int leptonXIndex, const int leptonYIndex,
//                                      const VLV& allLeptons, const std::vector<int>& lepPdgId)const
// {
//     if(!isMC_) return 1.;
//     if(leptonXIndex<0 || leptonYIndex<0) return 1.;
//     if(!triggerScaleFactors_EtaDn_) return 1.;
//     return triggerScaleFactors_EtaDn_->getSF(leptonXIndex, leptonYIndex, allLeptons, lepPdgId);
// }


// float TopAnalysisPlain::weightBtagSF_PtUp(const std::vector<int>& jetIndices,
//                                   const VLV& jets, const std::vector<int>& jetFlavour,
//                                   const std::vector<float>& btagDiscriminators)const{
//     if(!isMC_) return 1.;
//     if(!btagScaleFactors_) return 1.;
//     return btagScaleFactors_PtUp_->getSF(jetIndices, jets, jetFlavour, btagDiscriminators);
// }
// float TopAnalysisPlain::weightBtagSF_PtDn(const std::vector<int>& jetIndices,
//                                   const VLV& jets, const std::vector<int>& jetFlavour,
//                                   const std::vector<float>& btagDiscriminators)const{
//     if(!isMC_) return 1.;
//     if(!btagScaleFactors_) return 1.;
//     return btagScaleFactors_PtDn_->getSF(jetIndices, jets, jetFlavour, btagDiscriminators);
// }
//
// float TopAnalysisPlain::weightBtagSF_EtaUp(const std::vector<int>& jetIndices,
//                                   const VLV& jets, const std::vector<int>& jetFlavour,
//                                   const std::vector<float>& btagDiscriminators)const{
//     if(!isMC_) return 1.;
//     if(!btagScaleFactors_) return 1.;
//     return btagScaleFactors_EtaUp_->getSF(jetIndices, jets, jetFlavour, btagDiscriminators);
// }
// float TopAnalysisPlain::weightBtagSF_EtaDn(const std::vector<int>& jetIndices,
//                                   const VLV& jets, const std::vector<int>& jetFlavour,
//                                   const std::vector<float>& btagDiscriminators)const{
//     if(!isMC_) return 1.;
//     if(!btagScaleFactors_) return 1.;
//     return btagScaleFactors_EtaDn_->getSF(jetIndices, jets, jetFlavour, btagDiscriminators);
// }
//
// float TopAnalysisPlain::weightBtagSF_LjetPtUp(const std::vector<int>& jetIndices,
//                                   const VLV& jets, const std::vector<int>& jetFlavour,
//                                   const std::vector<float>& btagDiscriminators)const{
//     if(!isMC_) return 1.;
//     if(!btagScaleFactors_) return 1.;
//     return btagScaleFactors_LjetPtUp_->getSF(jetIndices, jets, jetFlavour, btagDiscriminators);
// }
// float TopAnalysisPlain::weightBtagSF_LjetPtDn(const std::vector<int>& jetIndices,
//                                   const VLV& jets, const std::vector<int>& jetFlavour,
//                                   const std::vector<float>& btagDiscriminators)const{
//     if(!isMC_) return 1.;
//     if(!btagScaleFactors_) return 1.;
//     return btagScaleFactors_LjetPtDn_->getSF(jetIndices, jets, jetFlavour, btagDiscriminators);
// }
//
// float TopAnalysisPlain::weightBtagSF_LjetEtaUp(const std::vector<int>& jetIndices,
//                                   const VLV& jets, const std::vector<int>& jetFlavour,
//                                   const std::vector<float>& btagDiscriminators)const{
//     if(!isMC_) return 1.;
//     if(!btagScaleFactors_) return 1.;
//     return btagScaleFactors_LjetEtaUp_->getSF(jetIndices, jets, jetFlavour, btagDiscriminators);
// }
// float TopAnalysisPlain::weightBtagSF_LjetEtaDn(const std::vector<int>& jetIndices,
//                                   const VLV& jets, const std::vector<int>& jetFlavour,
//                                   const std::vector<float>& btagDiscriminators)const{
//     if(!isMC_) return 1.;
//     if(!btagScaleFactors_) return 1.;
//     return btagScaleFactors_LjetEtaDn_->getSF(jetIndices, jets, jetFlavour, btagDiscriminators);
// }



float TopAnalysisPlain::weightKinReco_Up()const{
    if(!isMC_) return 1.;
    if(!kinematicReconstructionScaleFactors_Up_) return 1.;
    return kinematicReconstructionScaleFactors_Up_->getSF();
}
float TopAnalysisPlain::weightKinReco_Dn()const{
    if(!isMC_) return 1.;
    if(!kinematicReconstructionScaleFactors_Dn_) return 1.;
    return kinematicReconstructionScaleFactors_Dn_->getSF();
}


float TopAnalysisPlain::weightLooseKinReco_Up()const{
    if(!isMC_) return 1.;
    if(!looseKinRecoScaleFactors_Up_) return 1.;
    return looseKinRecoScaleFactors_Up_->getSF();
}
float TopAnalysisPlain::weightLooseKinReco_Dn()const{
    if(!isMC_) return 1.;
    if(!looseKinRecoScaleFactors_Dn_) return 1.;
    return looseKinRecoScaleFactors_Dn_->getSF();
}

void TopAnalysisPlain::SetPlainTreeOutputDirectory(const TString& plainTreeOutputDirectory)
{
   dirOut_ = plainTreeOutputDirectory;
}

// reimplemented to store true level info for ttbar other decays
void TopAnalysisPlain::SetRunViaTau(const bool& runViaTau)
{
  runViaTau_ = runViaTau;
  //if(runViaTau_)
  //this->SetTopSignal(true);
}

void TopAnalysisPlain::SetRecoTree(const bool& makeRecoTree)
{
  fillRecoTree_=makeRecoTree;
}
void TopAnalysisPlain::SetSimpleRecoTree(const bool& makeSimpleRecoTree)
{
  fillSimpleRecoTree_=makeSimpleRecoTree;
}
void TopAnalysisPlain::SetGenTree(const bool& makeGenTree)
{
  fillGenTree_=makeGenTree;
}

void TopAnalysisPlain::SetTTBarSample(const bool& isTTBarMC)
{
  isTTbarSample_=isTTBarMC;
}
void TopAnalysisPlain::SetSystSample(const bool& isSystMC)
{
  isSystSample_=isSystMC;
}
void TopAnalysisPlain::SetSystVarSample(const bool& isSystVar)
{
  isSystVariation_=isSystVar;
}


// reimplemented to treat ttbar other decays as signal
bool TopAnalysisPlain::failsTopGeneratorSelection(const Long64_t& entry)const
{
  int e = entry;

    if(!this->isTtbarPlusTauSample()) return false;
    // topDecayMode contains the decay of the top (*10) + the decay of the antitop (plus 100 or 200 for non-b decays of tops)
    // 1=hadron, 2=e, 3=mu, 4=tau->hadron, 5=tau->e, 6=tau->mu
    // i.e. 23 == top decays to e, tbar decays to mu
    const int topDecayMode = this->topDecayMode(e) % 100;
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

    const bool isBackgroundInSignalSample = !isCorrectChannel || isViaTau;
    if(runViaTau_ != isBackgroundInSignalSample) return true;
    return false;
}

// added plain tree creating stuff
void TopAnalysisPlain::Begin(TTree* tree)
{
    TopAnalysis::Begin(tree);

    // plain tree stuff

    // must be no open output file and no existing plain trees at this point
    if(fileOut_)
    {
      printf("****** Error in TopAnalysisPlain::Begin() apparently open output file %s\n", fileOut_->GetName());
      throw std::logic_error("TopAnalysisPlain::Begin() apparently open output file");
    }
    if(vecTreePlain_.size())
    {
      printf("****** Error in TopAnalysisPlain::Begin() apparently existing plain trees\n");
      throw std::logic_error("TopAnalysisPlain::Begin() apparently existing plain trees");
    }

    // open new output file
    TString fileNameFull = common::assignFolder(dirOut_, this->channel(), this->systematic());
    fileNameFull += this->outputFilename();
    printf("****** Opening file for PlainTree %s\n", fileNameFull.Data());
    fileOut_ = TFile::Open(fileNameFull, "recreate");
    if(!fileOut_ || fileOut_->IsZombie())
    {
      printf("****** Error: cannot create file %s\n", fileNameFull.Data());
      throw std::runtime_error("TopAnalysisPlain: cannot create file");
    }
    fileOut_->cd();

    // create plain tree setups if not provided
    //if(VecTreePlainSetup.size() == 0)
    //{
      std::vector<TreePlainSetup> vecSetups;
      // gen level
      // if(this->outputFilename().Contains("ttbarsignalplustau"))
      //if(this->isTopSignal())
      if(fillGenTree_){
        TreePlainSetup setup;
        setup.isGen = (this->isMC());
        setup.name = "plainTree_gen_step0";
        setup.title = "plainTree_gen_step0";
        setup.fillRecoTree = false;
        setup.fillSimpleRecoTree = false;
        setup.fillGenTree = true;
        setup.isTTbarSample = isTTbarSample_;
        setup.isSystSample = isSystSample_;
        setup.isSystVariation = isSystVariation_;
        vecSetups.push_back(setup);
      }

      // reco level
      if(fillRecoTree_){
        TreePlainSetup setup;
        setup.isGen = (this->isMC());
        setup.name = "plainTree_rec_step8";
        setup.title = "plainTree_rec_step8";
        setup.fillRecoTree = true;
        setup.fillSimpleRecoTree = false;
        setup.fillGenTree = false;
        setup.isTTbarSample = isTTbarSample_;
        setup.isSystSample = isSystSample_;
        setup.isSystVariation = isSystVariation_;
        vecSetups.push_back(setup);
      }
      if(fillSimpleRecoTree_){
        TreePlainSetup setup;
        setup.isGen = (this->isMC());
        // setup.name = "plainTree_recSimple_step8";
        // setup.title = "plainTree_recSimple_step8";
        setup.name = "plainTree_recSimple_step7L";
        setup.title = "plainTree_recSimple_step7L";
        setup.fillRecoTree = false;
        setup.fillSimpleRecoTree = true;
        setup.fillGenTree = false;
        setup.isTTbarSample = isTTbarSample_;
        setup.isSystSample = isSystSample_;
        setup.isSystVariation = isSystVariation_;
        vecSetups.push_back(setup);
      }
    //}

    // create plain trees
    //const long maxTreeSize = 6e9; // 6 GB
    // for(const auto& setup : vecSetups)
    for(std::vector<TreePlainSetup>::iterator setup=vecSetups.begin(); setup!=vecSetups.end(); ++setup)    //const auto TreePlainSetup& setup : vecSetups
    {
      printf("****** Creating TreePlain %s\n", setup->name.Data());
      const TString title = setup->name;
      // TreePlain* tree = new TreePlain(setup.Name, title, setup.IsGen);
      TreePlain* tree = new TreePlain(*setup);
      tree->InitWriteTree();
      //tree->TTreePtr->SetAutoSave(-maxTreeSize);
      //tree->TTreePtr->SetMaxTreeSize(maxTreeSize);
      vecTreePlain_.push_back(tree);
    }
}

// added plain tree filling stuff
void TopAnalysisPlain::fillAll(const std::string& selectionStep,
                             const EventMetadata& eventMetadata,
                             const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
			     const TopGenObjects& topGenObjects, const TopPseudoObjects& topPseudoObjects,
                             const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                             const LooseKinRecoSolution& looseKinRecoSolution,
                             const ttbar::GenObjectIndices& genObjectIndices, const ttbar::RecoObjectIndices& recoObjectIndices,
                             const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                             const double& defaultWeight,const double& defaultWeightLooseKinReco)const
{
    TopAnalysis::fillAll(selectionStep,
                            eventMetadata,
                            recoObjects, commonGenObjects,
                            topGenObjects, topPseudoObjects,
                            kinematicReconstructionSolutions,
                            looseKinRecoSolution,
                            genObjectIndices, recoObjectIndices,
                            genLevelWeights, recoLevelWeights,
                            defaultWeight,defaultWeightLooseKinReco);

    // OZ plain tree
    if(vecTreePlain_.size())
    {
      //this->GetPdfWeightsEntry(eventCounter_);


      const int& leadLepIdx = recoObjectIndices.leadingLeptonIndex_;
      const int& nLeadLepIdx = recoObjectIndices.nLeadingLeptonIndex_;
      const VLV& leptons = *recoObjects.allLeptons_;
      const std::vector<int>& pdgids = *recoObjects.lepPdgId_;
      const std::vector<float>& etas = *recoObjects.lepSCEta_;

      AnalysisBaseEntry entry;
      entry.analysisConfig = &analysisConfig_;
      entry.eventCounter = &eventCounter_;
      // entry.analysisBase = (TopAnalysisPlain*) this;
      entry.analysisBase = this;
      entry.eventMetadata = &eventMetadata;
      entry.commonGenObjects = &commonGenObjects;
      entry.recoObjects = &recoObjects;
      entry.topGenObjects = &topGenObjects;
      entry.topPseudoObjects = &topPseudoObjects;
      entry.kinematicReconstructionSolutions = &kinematicReconstructionSolutions;
      entry.looseKinRecoSolution = &looseKinRecoSolution;
      entry.recoObjectIndices = &recoObjectIndices;
      entry.genObjectIndices = &genObjectIndices;
      entry.genLevelWeights = &genLevelWeights;
      entry.weight = &defaultWeight;
      entry.weightLooseKinReco = &defaultWeightLooseKinReco;

      entry.jetIndices = &recoObjectIndices.jetIndices_;
      entry.jets = recoObjects.jets_;
      entry.jetFlavour = commonGenObjects.jetFlavour_;
      entry.btagDiscriminators = recoObjects.jetBtags_;

      entry.entryNumber = *entry.eventCounter;
      if(entry.analysisBase->b_eventNumber->GetEntryNumber() != *entry.eventCounter){
        entry.entryNumber = *entry.eventCounter;
      }else{
        entry.entryNumber = this->b_eventNumber->GetEntryNumber();
      }
      // float deltaRJetLep = this->analysisConfig_.selections().genDeltaRLeptonJetCut_;

      //correcting shift by -1 in entry value to make the resulting event_id consistent with the outputs (i.e. counting) from (in) TopAnalysis.cc
      entry.entryNumber--;

      for(auto tree : vecTreePlain_)
      {
        TString selectionStepTString(selectionStep);
        TString tempTreeNameTString(tree->TTreePtr->GetName());
        if(tempTreeNameTString.EndsWith(selectionStepTString)){


          if(tree->fillRecoTree_ || tree->fillSimpleRecoTree_){
            entry.lepSF[0] = this->weightLeptonSF(leadLepIdx, nLeadLepIdx, leptons, pdgids, etas, "LEPT_UP");
            entry.lepSF[1] = this->weightLeptonSF(leadLepIdx, nLeadLepIdx, leptons, pdgids, etas);
            entry.lepSF[2] = this->weightLeptonSF(leadLepIdx, nLeadLepIdx, leptons, pdgids, etas, "LEPT_DOWN");
              
            entry.lepSF[3] = this->weightLeptonSF(leadLepIdx, nLeadLepIdx, leptons, pdgids, etas, "ELE_ID_UP");
            entry.lepSF[4] = this->weightLeptonSF(leadLepIdx, nLeadLepIdx, leptons, pdgids, etas, "ELE_ID_DOWN");
            entry.lepSF[5] = this->weightLeptonSF(leadLepIdx, nLeadLepIdx, leptons, pdgids, etas, "ELE_RECO_UP");
            entry.lepSF[6] = this->weightLeptonSF(leadLepIdx, nLeadLepIdx, leptons, pdgids, etas, "ELE_RECO_DOWN");

            entry.lepSF[7] = this->weightLeptonSF(leadLepIdx, nLeadLepIdx, leptons, pdgids, etas, "MUON_ID_UP");
            entry.lepSF[8] = this->weightLeptonSF(leadLepIdx, nLeadLepIdx, leptons, pdgids, etas, "MUON_ID_DOWN");
            entry.lepSF[9] = this->weightLeptonSF(leadLepIdx, nLeadLepIdx, leptons, pdgids, etas, "MUON_ISO_UP");
            entry.lepSF[10] = this->weightLeptonSF(leadLepIdx, nLeadLepIdx, leptons, pdgids, etas, "MUON_ISO_DOWN");
            
            entry.kinRecoSF[0] = this->weightKinReco_Up();
            entry.kinRecoSF[1] = this->weightKinReco();
            entry.kinRecoSF[2] = this->weightKinReco_Dn();

            if(analysisConfig_.corrections().triggerSFReadoutMode_== 0){
                entry.trigSF[0] = 1.;
                entry.trigSF[1] = 1.;
                entry.trigSF[2] = 1.;
                // entry.trigSF[3] = 1.;
                // entry.trigSF[4] = 1.;
            }
            else if(analysisConfig_.corrections().triggerSFReadoutMode_== 2){
                entry.trigSF[0] = this->weightTriggerSF_Up(leadLepIdx, nLeadLepIdx, leptons, pdgids);
                entry.trigSF[1] = this->weightTriggerSF(leadLepIdx, nLeadLepIdx, leptons, pdgids);
                entry.trigSF[2] = this->weightTriggerSF_Dn(leadLepIdx, nLeadLepIdx, leptons, pdgids);
                // entry.trigSF[3] = this->weightTriggerSF_EtaUp(leadLepIdx, nLeadLepIdx, leptons, pdgids);
                // entry.trigSF[4] = this->weightTriggerSF_EtaDn(leadLepIdx, nLeadLepIdx, leptons, pdgids);
            }
            else {
                std::cerr << "TopAnalysisPlain::fillAll -->> " << analysisConfig_.corrections().triggerSFReadoutMode_ << " not implemented for plainTree production" << std::endl;
            }

            tree->SetBtagScaleFactors_Var_AllSyst(mapAllSystBTagScaleFactors_);
            entry.nominalBtagSF = entry.analysisBase->weightBtagSF(*entry.jetIndices, *entry.jets, *entry.jetFlavour, *entry.btagDiscriminators);

            entry.l1prefiringW[0] = entry.analysisBase->weightL1Prefiring_Up(recoObjects);
            entry.l1prefiringW[1] = entry.analysisBase->weightL1Prefiring(recoObjects);
            entry.l1prefiringW[2] = entry.analysisBase->weightL1Prefiring_Dn(recoObjects);
          }
          entry.pileupW[0] = weightPileup_Up(entry.entryNumber);
          entry.pileupW[1] = weightPileup(entry.entryNumber);
          entry.pileupW[2] = weightPileup_Dn(entry.entryNumber);
          // tree->Fill(&entry, deltaRJetLep);
          tree->Fill(&entry);
        }
      }
    }
}

// added plain tree closing stuff
void TopAnalysisPlain::Terminate()
{
    // Produce b-tag efficiencies if required for given correction mode
    this->produceBtagEfficiencies();

    // write plain trees
    TString outName;
    if(fileOut_)
    {
      fileOut_->cd();
      for(auto tree : vecTreePlain_)
      {
        printf("****** Storing TreePlain %s with %lld entries in %s\n", tree->TTreePtr->GetName(), tree->TTreePtr->GetEntries(), fileOut_->GetName());
        tree->TTreePtr->Write("",kWriteDelete);
        outName =  tree->TTreePtr->GetName();
      }
      h_weightedEvents->Write();
      TObjString(Channel::convert(channel_)).Write("channelName");
      TObjString(systematic_.name()).Write("systematicsName");
      TObjString(samplename_).Write("sampleName");
      TObjString(isTopSignal_ ? "1" : "0").Write("isSignal");
      TObjString(isHiggsSignal_ ? "1" : "0").Write("isHiggsSignal");
      TObjString(isMC_ ? "1" : "0").Write("isMC");

      // const double globalNormalisationFactor = trueLevelWeightSum_ != 0. ? trueLevelNoRenormalisationWeightSum_/trueLevelWeightSum_ : -999.;
      if (this->storeTrueLevelWeightSumHists())
      {
        std::unique_ptr<TH1D> trueLevelNoRenormalisationWeightSum_hist = std::make_unique<TH1D>(  TH1D("trueLevelNoRenormalisationWeightSum_hist","trueLevelNoRenormalisationWeightSum_hist",1,0.,1.) )  ;
        std::unique_ptr<TH1D> trueLevelWeightSum_hist = std::make_unique<TH1D>(  TH1D("trueLevelWeightSum_hist","trueLevelWeightSum_hist",1,0.,1.) );

        trueLevelNoRenormalisationWeightSum_hist->Fill(0.5,trueLevelNoRenormalisationWeightSum_);
        trueLevelWeightSum_hist->Fill(0.5,trueLevelWeightSum_);
        trueLevelNoRenormalisationWeightSum_hist->Write();
        trueLevelWeightSum_hist->Write();
      }

      if(this->readInAllSyst())
      {
        fileOut_->mkdir("normalization_histograms");
        fileOut_->cd("normalization_histograms");
        for(auto const& x : mapAllSystTrueLevelWeightSum_){
          TString name = "trueLevelWeightSum_hist_"+Systematic::convertType(std::get<0>(x.first))+Systematic::convertVariation(std::get<1>(x.first));
          if(std::get<0>(x.first)==Systematic::pdf || std::get<0>(x.first)==Systematic::psScaleWeight){
            TString addition;
            addition.Form("_%i", std::get<2>(x.first));
            name=name+addition;
          }
          std::unique_ptr<TH1D>  trueLevelWeightSum_hist_temp= std::make_unique<TH1D>(  TH1D(name, name, 1,0.,1.) ) ;
          trueLevelWeightSum_hist_temp->Fill(0.5, x.second);
          trueLevelWeightSum_hist_temp->Write();

        }
        for(auto const& x : mapAllSystTrueLevelNoRenormalisationWeightSum_)
        {
          TString name = "trueLevelNoRenormalisationWeightSum_hist"+Systematic::convertType(std::get<0>(x.first))+Systematic::convertVariation(std::get<1>(x.first));
          std::unique_ptr<TH1D>  trueLevelNoRenormalisationWeightSum_hist_temp = std::make_unique<TH1D>(  TH1D(name, name, 1,0.,1.) ) ;
          trueLevelNoRenormalisationWeightSum_hist_temp->Fill(0.5, x.second);
          trueLevelNoRenormalisationWeightSum_hist_temp->Write();
        }
      }

      trueLevelWeightSum_ = 0.;
      trueLevelNoRenormalisationWeightSum_ = 0.;

      fileOut_->Close();
      fileOut_ = NULL;
      for(auto tree : vecTreePlain_)
        delete tree;
      vecTreePlain_.clear();
    }

    // Do everything needed for TTree
    for(TreeHandlerBase* treeHandler : v_treeHandler_){
        if(treeHandler){
            // Produce and write tree
            treeHandler->writeTrees(this->outputFilename(), this->channel(), this->systematic());
            //treeHandler->writeTrees(fOutput);

            // Cleanup
            treeHandler->clear();
        }
    }
    // Defaults from AnalysisBase
    AnalysisBase::Terminate();
}
