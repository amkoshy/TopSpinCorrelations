#include "MiniTreeTopAnalysis.h"
#include "AnalyzerBaseClass.h"
#include "TreeHandlerBase.h"
#include "ScaleFactors.h"
#include "AnalysisConfig.h"
#include <TTree.h>
#include <TH1.h>
#include "../../common/include/KinematicReconstructionSolution.h"
#include "../../common/include/KinematicReconstruction.h"
#include "../../common/include/LooseKinRecoSolution.h"
#include "../../common/include/LooseKinReco.h"



void MiniTreeTopAnalysis::SetTriggerScaleFactors_Up(const TriggerScaleFactors* const scaleFactors)
{
    triggerScaleFactors_Up_ = scaleFactors;
}
void MiniTreeTopAnalysis::SetTriggerScaleFactors_Dn(const TriggerScaleFactors* const scaleFactors)
{
    triggerScaleFactors_Dn_ = scaleFactors;
}

void MiniTreeTopAnalysis::SetBtagScaleFactors_Var_AllSyst(const std::map<std::tuple<Systematic::Type, Systematic::Variation, int>, BtagScaleFactors*>* mapAllSystBTagScaleFactors)
{
    mapAllSystBTagScaleFactors_ = mapAllSystBTagScaleFactors;
}

  // kinematicReconstruction as dummy to keep same structure
void MiniTreeTopAnalysis::SetKinematicReconstruction_Up(const KinematicReconstruction* const kinematicReconstruction,
                                              const KinematicReconstructionScaleFactors* const kinematicReconstructionScaleFactors)
{
    (void) kinematicReconstruction;
    kinematicReconstructionScaleFactors_Up_ = kinematicReconstructionScaleFactors;
}
void MiniTreeTopAnalysis::SetKinematicReconstruction_Dn(const KinematicReconstruction* const kinematicReconstruction,
                                              const KinematicReconstructionScaleFactors* const kinematicReconstructionScaleFactors)
{
    (void) kinematicReconstruction;
    kinematicReconstructionScaleFactors_Dn_ = kinematicReconstructionScaleFactors;
}

  // LooseKinReco as dummy to keep same structure
void MiniTreeTopAnalysis::SetLooseKinReco_Up(const LooseKinReco* const looseKinReco,
					  const LooseKinRecoScaleFactors* const looseKinRecoScaleFactors)
{
    (void) looseKinReco;
    looseKinRecoScaleFactors_Up_ = looseKinRecoScaleFactors;
}
void MiniTreeTopAnalysis::SetLooseKinReco_Dn(const LooseKinReco* const looseKinReco,
					  const LooseKinRecoScaleFactors* const looseKinRecoScaleFactors)
{
    (void) looseKinReco;
    looseKinRecoScaleFactors_Dn_ = looseKinRecoScaleFactors;
}


void MiniTreeTopAnalysis::SetRunViaTau(const bool runViaTau)
{
    runViaTau_ = runViaTau;
    // if(runViaTau) this->SetTopSignal(false);
}

void MiniTreeTopAnalysis::SetTreeOutputDirectory(const TString& treeOutputDirectory)
{
   dirOut_ = treeOutputDirectory;
}


void MiniTreeTopAnalysis::SetTTBarSample(const bool& isTTBarMC)
{
  isTTbarSample_=isTTBarMC;
}

void MiniTreeTopAnalysis::SetSystSample(const bool& isSystMC)
{
  isSystSample_=isSystMC;
}
void MiniTreeTopAnalysis::SetSystVarSample(const bool& isSystVar)
{
  isSystVariation_=isSystVar;
}


void MiniTreeTopAnalysis::Begin(TTree* tree)
{
    TopAnalysis::Begin(tree);


    // must be no open output file and no existing trees at this point
    if(fileOut_)
    {
      printf("****** Error in MiniTreeTopAnalysis::Begin() apparently open output file %s\n", fileOut_->GetName());
      throw std::logic_error("MiniTreeTopAnalysis::Begin() apparently open output file");
    }
    if(vecMiniTree_.size())
    {
      printf("****** Error in MiniTreeTopAnalysis::Begin() apparently existing trees\n");
      throw std::logic_error("MiniTreeTopAnalysis::Begin() apparently existing trees");
    }

    // open new output file
    TString fileNameFull = common::assignFolder(dirOut_, this->channel(), this->systematic());
    fileNameFull += this->outputFilename();
    printf("****** Opening file for Tree %s\n", fileNameFull.Data());
    fileOut_ = TFile::Open(fileNameFull, "recreate");
    if(!fileOut_ || fileOut_->IsZombie())
    {
      printf("****** Error: cannot create file %s\n", fileNameFull.Data());
      throw std::runtime_error("MiniTreeTopAnalysis: cannot create file");
    }
    fileOut_->cd();

    // create tree setups if not provided
      std::vector<TreeSetup> vecSetups;
        TreeSetup setup;
        setup.isGen = (this->isMC());
        setup.name = "miniTree";
        setup.title = "miniTree";
        setup.isTTbarSample = isTTbarSample_;
        setup.isSystSample = isSystSample_;
        setup.isSystVariation = isSystVariation_;
        vecSetups.push_back(setup);


    // create trees
    for(std::vector<TreeSetup>::iterator setup=vecSetups.begin(); setup!=vecSetups.end(); ++setup)    //const auto TreeSetup& setup : vecSetups
    {
      printf("****** Creating %s\n", setup->name.Data());
      const TString title = setup->name;
      MiniTree* tree = new MiniTree(*setup);
      tree->InitWriteTree();
      vecMiniTree_.push_back(tree);
    }
}

// added tree filling stuff
void MiniTreeTopAnalysis::fillAll(const std::string& selectionStep,
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
    if(vecMiniTree_.size())
    {
        AnalysisEntry entry;
        entry.eventCounter = &eventCounter_;
        entry.eventCounterCorrected = *entry.eventCounter -1;
        entry.eventCounter = &entry.eventCounterCorrected;
        entry.analysisBase = this;
        entry.eventMetadata = &eventMetadata;
        entry.commonGenObjects = &commonGenObjects;
        entry.analysisConfig = &analysisConfig_;
        entry.recoObjects = &recoObjects;
        entry.topGenObjects = &topGenObjects;
        entry.topPseudoObjects = &topPseudoObjects;
        entry.genObjectIndices = &genObjectIndices;
        entry.genLevelWeights = &genLevelWeights;
        entry.weight = defaultWeight;


        entry.entryNumber = *entry.eventCounter;
        if(entry.analysisBase->b_eventNumber->GetEntryNumber() != *entry.eventCounter){
            entry.entryNumber = *entry.eventCounter;
        }else{
            entry.entryNumber = this->b_eventNumber->GetEntryNumber();
        }

        for(auto tree : vecMiniTree_)
        {
            TString selectionStepTString(selectionStep);
            if(selectionStepTString.Contains("lhAllGenLevel")){

                tree->SetTriggerScaleFactors_Up(triggerScaleFactors_Up_);
                tree->SetTriggerScaleFactors_Dn(triggerScaleFactors_Dn_);

                tree->SetBtagScaleFactors_Var_AllSyst(mapAllSystBTagScaleFactors_);

                tree->SetKinematicReconstruction_Up(kinematicReconstruction_Up_, kinematicReconstructionScaleFactors_Up_);
                tree->SetKinematicReconstruction_Dn(kinematicReconstruction_Dn_, kinematicReconstructionScaleFactors_Dn_);

                tree->SetLooseKinReco_Up(looseKinReco_Up_, looseKinRecoScaleFactors_Up_);
                tree->SetLooseKinReco_Dn(looseKinReco_Dn_, looseKinRecoScaleFactors_Dn_);

                entry.pileupW_Up = weightPileup_Up(entry.entryNumber);
                entry.pileupW = weightPileup(entry.entryNumber);
                entry.pileupW_Down = weightPileup_Dn(entry.entryNumber);

                tree->Fill(&entry);
            }
        }
    }
}

// added tree closing stuff
void MiniTreeTopAnalysis::Terminate()
{
    // Produce b-tag efficiencies if required for given correction mode
    this->produceBtagEfficiencies();

    // write trees
    TString outName;
    if(fileOut_)
    {
      fileOut_->cd();
      for(auto tree : vecMiniTree_)
      {
        printf("****** Storing %s with %lld entries in %s\n", tree->TTreePtr->GetName(), tree->TTreePtr->GetEntries(), fileOut_->GetName());
        tree->TTreePtr->Write("",kWriteDelete);
        outName =  tree->TTreePtr->GetName();
      }

      // Write additional information into file
      h_weightedEvents->Write();
      TObjString(Channel::convert(channel_)).Write("channelName");
      TObjString(systematic_.name()).Write("systematicsName");
      TObjString(samplename_).Write("sampleName");
      TObjString(isTopSignal_ ? "1" : "0").Write("isSignal");
      TObjString(isHiggsSignal_ ? "1" : "0").Write("isHiggsSignal");
      TObjString(isMC_ ? "1" : "0").Write("isMC");

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
      for(auto tree : vecMiniTree_)
        delete tree;
      vecMiniTree_.clear();
    }

    // Do everything needed for TTree
    for(TreeHandlerBase* treeHandler : v_treeHandler_){
        if(treeHandler){
            // Produce and write tree
            treeHandler->writeTrees(this->outputFilename(), this->channel(), this->systematic());
            // Cleanup
            treeHandler->clear();
        }
         std::cout<<"Created: \033[1;1m"<<outName<<"\033[0m\n\n";
    }

    // Defaults from AnalysisBase
    AnalysisBase::Terminate();
}
