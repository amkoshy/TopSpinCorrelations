#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TProof.h>
#include <TSelector.h>
#include <TObjString.h>
#include <TChain.h>
#include <TH1.h>
#include <Rtypes.h>

#include "TopAnalysis.h"
#include "TopAnalysisPlain.h"
#include "MiniTreeTopAnalysis.h"
#include "AnalysisConfig.h"
#include "analysisHelpers.h"
#include "AnalyzerBaseClass.h"
#include "AnalyzerControlPlots.h"
#include "AnalyzerBoostedTop.h"
#include "AnalyzerSpinCorr.h"
#include "AnalyzerGenObjs.h"
#include "AnalyzerDoubleDiffXS.h"
#include "AnalyzerKinReco.h"
#include "AnalyzerLooseKinReco.h"
#include "AnalyzerKinRecoQualityStudies.h"
#include "TreeHandlerBase.h"
#include "TreeHandlerTTBar.h"
#include "TreeHandlerSpinCorr.h"
#include "TreePlain.h"
#include "MiniTree.h"
#include "TreeHandlerBoostedTop.h"
#include "utils.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/utils.h"
#include "../../common/include/CommandLineParameters.h"
#include "../../common/include/KinematicReconstruction.h"
#include "../../common/include/LooseKinReco.h"
#include "../../common/include/ScaleFactors.h"





void load_Analysis(const TString& validFilenamePattern,
                   const Channel::Channel& channel,
                   const Systematic::Systematic& systematic,
                   const std::vector<AnalysisMode::AnalysisMode>& v_analysisMode,
		           const int specific_PDF,
		           const int specific_PS,
                   const int dy,
                   const TString& closure,
                   const float& slope,
                   const Long64_t& maxEvents,
                   const Long64_t& skipEvents,
                   const bool runsignalviatau,
                   const bool useTopMassSetup,
                   const bool dryRun,
                   const std::string& commandline
                  )
{
    std::cout<<std::endl;

    // Read analysis config from text file
    const AnalysisConfig analysisConfig;
    const AnalysisConfig::Corrections& corrections = analysisConfig.corrections();
    //analysisConfig.print(true);

    // Set analysis year suffix from era
    std::string yearSuffix = "";
    if (analysisConfig.general().era_ == Era::run2_13tev_2016) yearSuffix = "_2016";
    else if (analysisConfig.general().era_ == Era::run2_13tev_2017) yearSuffix = "_2017";
    else if (analysisConfig.general().era_ == Era::run2_13tev_2018) yearSuffix = "_2018";

    else if (analysisConfig.general().era_ == Era::run2_UL_13tev_2016preVFP) yearSuffix = "_2016preVFP";
    else if (analysisConfig.general().era_ == Era::run2_UL_13tev_2016postVFP) yearSuffix = "_2016postVFP";
    else if (analysisConfig.general().era_ == Era::run2_UL_13tev_2017) yearSuffix = "_2017";
    else if (analysisConfig.general().era_ == Era::run2_UL_13tev_2018) yearSuffix = "_2018";

    else std::cout << "Era-to-year conversion is not implemented. Proceeding with the setup without assigning dedicated output folders for each year.\n" << std::endl;

    // Set up the channels to run over
    std::vector<Channel::Channel> channels;
    if(channel != Channel::undefined) channels.push_back(channel);
    else channels = Channel::realChannels;

    // Set up kinematic reconstruction
    const KinematicReconstruction* kinematicReconstruction(0);
    kinematicReconstruction = new KinematicReconstruction(analysisConfig.general().era_, 1, true);

    const LooseKinReco* looseKinReco(0);
    looseKinReco = new LooseKinReco(analysisConfig.general().era_);

    // Set up kinematic reconstruction scale factors (null-pointer means no application)

    KinematicReconstructionScaleFactors* kinematicReconstructionScaleFactors(0);
    kinematicReconstructionScaleFactors = new KinematicReconstructionScaleFactors(analysisConfig.general().era_, channels, systematic);
    if(systematic.type()==Systematic::kin && !kinematicReconstructionScaleFactors){
      std::cout<<"Systematic for kinematic reconstruction requested, but scale factors not applied"
	       <<"\nStop running of analysis --- this is NOT an error, just avoiding double running\n"<<std::endl;
      exit(1);
    }
    LooseKinRecoScaleFactors* looseKinRecoScaleFactors(0);
    looseKinRecoScaleFactors = new LooseKinRecoScaleFactors(analysisConfig.general().era_, channels, systematic);
    if(systematic.type()==Systematic::kin && !looseKinRecoScaleFactors){
      std::cout<<"Systematic for loose kinematic reconstruction requested, but scale factors not applied"
	       <<"\nStop running of analysis --- this is NOT an error, just avoiding double running\n"<<std::endl;
      exit(1);
    }

    // Set up lepton efficiency scale factors
    LeptonScaleFactors* leptonScaleFactors(0);
    if(corrections.combinedSF_== true){
      leptonScaleFactors = new LeptonScaleFactors(corrections.electronSFInputFile_, corrections.muonSFInputFile_, corrections.leptonSFReadoutMode_, systematic);
    }else{
      if(corrections.combinedSF_==false){
        lepSF::SFHisto electronSFID_(corrections.SF_electron_ID_File_, corrections.SF_electron_ID_Histo_,
          corrections.SF_electron_ID_Format_, "","");
        lepSF::SFHisto electronSFReco_(corrections.SF_electron_Reco_File_,
          corrections.SF_electron_Reco_Histo_, corrections.SF_electron_Reco_Format_, "", "");
        lepSF::SFHisto electronSFReco_LowPt_(corrections.SF_electron_Reco_File_LowPt_,
          corrections.SF_electron_Reco_Histo_LowPt_, corrections.SF_electron_Reco_Format_LowPt_, "", "");
        lepSF::SFHisto muonSFID_(corrections.SF_muon_ID_File_, corrections.SF_muon_ID_Histo_,
          corrections.SF_muon_ID_Format_, corrections.SF_muon_ID_Histo_Stat_, corrections.SF_muon_ID_Histo_Sys_);
        lepSF::SFHisto muonSFISO_(corrections.SF_muon_ISO_File_, corrections.SF_muon_ISO_Histo_,
          corrections.SF_muon_ISO_Format_, corrections.SF_muon_ISO_Histo_Stat_, corrections.SF_muon_ISO_Histo_Sys_);

        leptonScaleFactors = new LeptonScaleFactors(electronSFID_, electronSFReco_, electronSFReco_LowPt_,
          muonSFID_, muonSFISO_, systematic);
      }
    }
    leptonScaleFactors->setDYExtrapolationUncFactors(corrections.extra_muon_ISO_DY_extrapolation_unc_, corrections.extra_electron_ID_DY_extrapolation_unc_);

    // Set up trigger efficiency scale factors
    TriggerScaleFactors* triggerScaleFactors(0);
    if(corrections.triggerSFReadoutMode_==1){
      triggerScaleFactors = new TriggerScaleFactors(corrections.triggerSFInputSuffix_, channels, corrections.triggerSFReadoutMode_, systematic);
    }
    else if(corrections.triggerSFReadoutMode_==2){
	       triggerSF::SFHisto eeSF_(corrections.trigger_DoubleElectron_SF_File_, corrections.trigger_DoubleElectron_SF_Histogram_,
	          corrections.trigger_DoubleElectron_SF_Histogram_Format_, "","");
        triggerSF::SFHisto emuSF_(corrections.trigger_ElectronMuon_SF_File_,
          corrections.trigger_ElectronMuon_SF_Histogram_, corrections.trigger_ElectronMuon_SF_Histogram_Format_, "", "");
        triggerSF::SFHisto mumuSF_(corrections.trigger_DoubleMuon_SF_File_,
          corrections.trigger_DoubleMuon_SF_Histogram_, corrections.trigger_DoubleMuon_SF_Histogram_Format_, "", "");

	      triggerScaleFactors = new TriggerScaleFactors(eeSF_, emuSF_, mumuSF_, corrections.triggerSFReadoutMode_, systematic, false);
      }
    else if (corrections.triggerSFReadoutMode_== 0){
    	std::cout << "** Warning: triggerSFReadoutMode_ = 0, probably due to no definition of this parameter in the config file. " << std::endl
    			  << "   TriggerSF weights will be taken as 1.0" << std::endl << std::endl;
        }
    else{
    	std::cerr << "triggerSFReadoutMode_ = " << corrections.triggerSFReadoutMode_ << " not implemented!" << std::endl
    	          << "Please use one of these options: 1, 2, or 0." << std::endl;
    	std::exit(1);

    }

    // Set up JER correction and systematic scale factors
    const JetEnergyResolutionScaleFactors*    jetEnergyResolutionScaleFactors   (nullptr);
    if(corrections.jet_JERHybridMethod_PtResolutionFile_ == ""){
        jetEnergyResolutionScaleFactors = new JetEnergyResolutionScaleFactors(corrections.jetCorrectionMode_, corrections.jerUncertaintySourceName_, systematic, corrections.jerSFBinning_);
    }else{
        const std::string jet_PtRes_File = common::DATA_PATH_COMMON()+"/"+corrections.jet_JERHybridMethod_PtResolutionFile_;

        jetEnergyResolutionScaleFactors    = new JetEnergyResolutionScaleFactors(corrections.jetCorrectionMode_, corrections.jerUncertaintySourceName_, systematic,
                                                                                 corrections.jet_JERHybridMethod_jetConeRadius_, jet_PtRes_File, corrections.jerSFBinning_);
    }
    // Set up JES correction and systematic scale factors
    std::string input_jesUncertaintySourceFile = "";
    if (corrections.useReducedJESSet_){
      if (corrections.jesUncertaintyReducedSourceFile_ == "") {
        std::cout<<"No jesUncertaintyReducedSourceFile have been provided in the config file !!!"<<std::endl;
        exit(1);
      }
      input_jesUncertaintySourceFile = corrections.jesUncertaintyReducedSourceFile_;
    }
    else{
      input_jesUncertaintySourceFile = corrections.jesUncertaintySourceFile_;
    }
    const JetEnergyScaleScaleFactors* jetEnergyScaleScaleFactors = new JetEnergyScaleScaleFactors(corrections.jetCorrectionMode_, input_jesUncertaintySourceFile, systematic, corrections.useReducedJESSet_);
    const JetEnergyScaleScaleFactors* jetEnergyScaleScaleFactors_L1only = new JetEnergyScaleScaleFactors(corrections.jetCorrectionMode_, input_jesUncertaintySourceFile, systematic, corrections.useReducedJESSet_);

    // Set up top-pt reweighting scale factors (null-pointer means no application)
    const TopPtScaleFactors* topPtScaleFactors(0);

    if (analysisConfig.corrections().topPtReweighting_) topPtScaleFactors = new TopPtScaleFactors(systematic);
    else if(systematic.type()==Systematic::topPt) topPtScaleFactors = new TopPtScaleFactors(Systematic::Systematic("Nominal"));

    if(systematic.type()==Systematic::topPt && !topPtScaleFactors){
        std::cout<<"Systematic for top-pt reweighting requested, but scale factors not applied"
                 <<"\nStop running of analysis --- this is NOT an error, just avoiding double running\n"<<std::endl;
        exit(1);
    }

    if(systematic.type()==Systematic::topPt){
        if (analysisConfig.corrections().topPtReweighting_){
            if (systematic.variation()==Systematic::undefinedVariation){
                std::cerr << "Inconsistent systematic parameter: " << systematic.name() << ".\n You are using top pT reweighting. For consistent treatment please use TOP_PT_UP and TOP_PT_DOWN.\n";
                std::exit(1);
            }
        }
        else if(systematic.variation()==Systematic::up || systematic.variation()==Systematic::down){
            std::cerr << "Inconsistent systematic parameter: " << systematic.name() << ".\n You are NOT using top pT reweighting. For consistent treatment please use TOP_PT.\n";
            std::exit(1);
        }
    }



    // Vector for setting up all analysers
    std::vector<AnalyzerBaseClass*> v_analyzer;

    // Histograms for Gen-Objects
    AnalyzerGenObjs* analyzerGenObjs(nullptr);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::gen) != v_analysisMode.end())
    {
      analyzerGenObjs = new AnalyzerGenObjs({"0", "1", "2", "3", "4", "5", "6", "7", "7L"});
      v_analyzer.push_back(analyzerGenObjs);
    }

    // Set up event yield histograms
    AnalyzerEventYields* analyzerEventYields(0);
    analyzerEventYields = new AnalyzerEventYields({"1", "2", "3", "4", "5", "6", "7", "8", "7L"});
    v_analyzer.push_back(analyzerEventYields);

    // Set up Drell-Yan scaling histograms
    AnalyzerDyScaling* analyzerDyScaling(0);
    analyzerDyScaling = new AnalyzerDyScaling({"4", "5", "6", "7", "8"}, "5");
    v_analyzer.push_back(analyzerDyScaling);

    // Set up basic histograms
    AnalyzerControlPlots* analyzerControlPlots(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::cp) != v_analysisMode.end()){
      //      analyzerControlPlots = new AnalyzerControlPlots({"3", "4", "4zWindow", "5", "5zWindow", "6", "6zWindow", "7", "7zWindow", "7L", "8", "8zWindow"});
      analyzerControlPlots = new AnalyzerControlPlots({"3", "4", "5", "6", "7", "8"});
        v_analyzer.push_back(analyzerControlPlots);
    }
    // Set up dda histograms
    AnalyzerDoubleDiffXS* analyzerDoubleDiffXS(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::dda) != v_analysisMode.end()){
        analyzerDoubleDiffXS = new AnalyzerDoubleDiffXS({"0","8"});
        v_analyzer.push_back(analyzerDoubleDiffXS);
    }

    // Set up boosted top histograms
    AnalyzerBoostedTop* analyzerBoostedTop(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::btop) != v_analysisMode.end()){
        analyzerBoostedTop = new AnalyzerBoostedTop({"0","8"});
        v_analyzer.push_back(analyzerBoostedTop);
    }

    // ajeeta                                                                                                                                                                                                            
    // Set up spin correlation histograms                                                                                                                                                                                
    AnalyzerSpinCorr* analyzerSpinCorr(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::spinCorr) != v_analysisMode.end()){
      analyzerSpinCorr = new AnalyzerSpinCorr(analysisConfig,analysisConfig.selections().leptonPtCut_, analysisConfig.selections().leptonEtaCut_, analysisConfig.selections().jetPtCut_, analysisConfig.selections().jetEtaCut_, analysisConfig.selections().genDeltaRLeptonJetCut_, {"8", "8_0to450_ttbarmass", "8_450to600_ttbarmass", "8_600to800_ttbarmass", "8_800toinf_ttbarmass", "8_m1tomhalf_topscatteringangle", "8_mhalfto0_topscatteringangle", "8_0tophalf_topscatteringangle", "8_phalftop1_topscatteringangle", "8_0_extrajets", "8_1_extrajets", "8_2_extrajets",  "8_3ormore_extrajets"});

      v_analyzer.push_back(analyzerSpinCorr);
    }

    // Set up KinReco histograms
    AnalyzerKinReco* analyzerKinReco(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::kinReco) != v_analysisMode.end()){
      analyzerKinReco = new AnalyzerKinReco(analysisConfig.general().era_, {"0","7","8"});
      v_analyzer.push_back(analyzerKinReco);
    }

    // Set up loose KinReco histograms
    AnalyzerLooseKinReco* analyzerLooseKinReco(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::looseKinReco) != v_analysisMode.end()){
      analyzerLooseKinReco = new AnalyzerLooseKinReco(analysisConfig.general().era_, {"7","8"});
      v_analyzer.push_back(analyzerLooseKinReco);
    }

    // Set up histograms for kin reco quality studies
    AnalyzerKinRecoQualityStudies* analyzerKinRecoQualityStudies(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::kinRecoQualityStudies) != v_analysisMode.end()){
      analyzerKinRecoQualityStudies = new AnalyzerKinRecoQualityStudies(analysisConfig.general().era_, {"7","8"});
      v_analyzer.push_back(analyzerKinRecoQualityStudies);
    }

    // Vector setting up all tree handlers for MVA input variables
    std::vector<TreeHandlerBase*> v_treeHandler;

    // Set up production of TTree for ttbar analysis
    TreeHandlerTTBar* treeHandlerTTBar(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::ddaTree) != v_analysisMode.end()){
        treeHandlerTTBar = new TreeHandlerTTBar("ddaInput", {"0","8"});
        v_treeHandler.push_back(treeHandlerTTBar);
    }

    // ajeeta                                                                                                                                                                                                            
    // Set up production of TTree for spin correlation analysis                                                                                                                                                          
    //TH1* wtedEvts = dynamic_cast<Th1*>(file->Get("EventsBeforeSelection/weightedEvents"));                                                                                                                              
    TreeHandlerSpinCorr* treeHandlerSpinCorr(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::spinCorr) != v_analysisMode.end()){

      if (yearSuffix == "_2016") treeHandlerSpinCorr = new TreeHandlerSpinCorr("spinCorrInput_2016", {"0","8"});
      if (yearSuffix == "_2016preVFP") treeHandlerSpinCorr = new TreeHandlerSpinCorr("spinCorrInput_2016preVFP", {"0","8"});
      if (yearSuffix == "_2016postVFP") treeHandlerSpinCorr = new TreeHandlerSpinCorr("spinCorrInput_2016postVFP", {"0","8"});
      if (yearSuffix == "_2017") treeHandlerSpinCorr = new TreeHandlerSpinCorr("spinCorrInput_2017", {"0","8"});
      if (yearSuffix == "_2018") treeHandlerSpinCorr = new TreeHandlerSpinCorr("spinCorrInput_2018", {"0","8"});

      v_treeHandler.push_back(treeHandlerSpinCorr);
    }

    // Set up production of TTree for Boosted Top analysis
    TreeHandlerBoostedTop* treeHandlerBoostedTop(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::btopTree) != v_analysisMode.end()){
        treeHandlerBoostedTop = new TreeHandlerBoostedTop("btopInput", {"0","8"});
        v_treeHandler.push_back(treeHandlerBoostedTop);
    }

    // Set up production of plain TTree (contains particle 4-momenta)
    bool isPlainTreeAnalysis2D = false;
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::plainTree) != v_analysisMode.end()){
      isPlainTreeAnalysis2D = true;
    }

    // Set up production of flat miniTrees
    bool isMiniTreeAnalysis = false;
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::miniTree) != v_analysisMode.end()){
      isMiniTreeAnalysis = true;
    }

    // Set up the analysis
    TopAnalysis* selector = NULL;
    // Use derived analyser class for 2D analysis with plain tree
    const TString outputPlainTreeFolder = "plainTree" + yearSuffix;
    const TString outputMiniTreeFolder = "miniTree" + yearSuffix;
    if (isPlainTreeAnalysis2D) {
      selector = new TopAnalysisPlain(analysisConfig);
      selector->SetPlainTreeOutputDirectory(outputPlainTreeFolder);
      selector->SetTrueLevelWeightSumHists(true);
      // selector->SetReadInAllSyst(true, -1, -1, 1);
    } else {
      if(isMiniTreeAnalysis){
        selector = new MiniTreeTopAnalysis(analysisConfig);
        selector->SetTreeOutputDirectory(outputMiniTreeFolder);
        selector->SetTrueLevelWeightSumHists(true);
        // selector->SetReadInAllSyst(true, -1, -1, 1);
      }else{
        selector = new TopAnalysis(analysisConfig);
      }
    }

    const std::string outputBaseFolder = "selectionRoot" + yearSuffix;
    const char* outputBase = outputBaseFolder.c_str();
    selector->SetDryRun(dryRun, commandline);
    selector->SetAnalysisOutputBase(outputBase);

    selector->SetKinematicReconstruction(kinematicReconstruction, kinematicReconstructionScaleFactors);

    selector->SetLooseKinReco(looseKinReco, looseKinRecoScaleFactors);

    selector->SetLeptonScaleFactors(leptonScaleFactors);
    selector->SetTriggerScaleFactors(triggerScaleFactors);
    selector->SetJetEnergyResolutionScaleFactors(jetEnergyResolutionScaleFactors);
    selector->SetJetEnergyScaleScaleFactors(jetEnergyScaleScaleFactors);
    selector->SetJetEnergyScaleScaleFactors_L1only(jetEnergyScaleScaleFactors_L1only);
    selector->SetTopPtScaleFactors(topPtScaleFactors);
    selector->SetAllAnalyzers(v_analyzer);
    selector->SetAllTreeHandlers(v_treeHandler);

    // Access selectionList containing all input sample nTuples
    const std::string selectionListName = "data/selectionList" + yearSuffix + ".txt";
    std::ifstream infile(selectionListName);
    if(!infile.good()){
        std::cerr << "Error! Please check the " << selectionListName << " file!\n" << std::endl;
        exit(773);
    }

    // Loop over all input files
    int filecounter = 0;
    while(!infile.eof()){

        // Access nTuple input filename from selectionList
        TString filename;
        infile>>filename;
        if(filename=="" || filename[0]=='#') continue; //empty or commented line? --> skip
        if(!filename.Contains(validFilenamePattern)) continue;
        std::cout << "\nPROCESSING File "<< ++filecounter << " (" << filename << ") from " << selectionListName << "\n" << std::endl;

        // Access the basic filename by stripping off the folders
        TString filenameBase(filename);
        if(filenameBase.Contains('/')){
            Ssiz_t last = filenameBase.Last('/');
            filenameBase = filenameBase.Data() + last + 1;
        }

        // Open nTuple file
        TFile *file = TFile::Open(filename);
	if(file->IsZombie()){
	  std::cerr<<"ERROR! Cannot open nTuple file with name: "<<filename<<std::endl;
	  exit(853);
	}
    

        // Check whether nTuple can be found
        TTree* tree = dynamic_cast<TTree*>(file->Get("writeNTuple/NTuple"));
        if(!tree){
            std::cerr<<"ERROR! TTree (=nTuple) not found in file!\n"<<std::endl;
            exit(854);
        }

        // Access information about samples, stored in the nTuples
        TObjString* channel_from_file = dynamic_cast<TObjString*>(file->Get("writeNTuple/channelName"));
        TObjString* systematics_from_file = dynamic_cast<TObjString*>(file->Get("writeNTuple/systematicsName"));
        TObjString* samplename = dynamic_cast<TObjString*>(file->Get("writeNTuple/sampleName"));
        TObjString* o_isSignal = dynamic_cast<TObjString*>(file->Get("writeNTuple/isSignal"));
        TObjString* o_isHiggsSignal = dynamic_cast<TObjString*>(file->Get("writeNTuple/isHiggsSignal"));
        TObjString* o_isMC = dynamic_cast<TObjString*>(file->Get("writeNTuple/isMC"));
        TH1* weightedEvents = dynamic_cast<TH1*>(file->Get("EventsBeforeSelection/weightedEvents"));
        if(!channel_from_file || !systematics_from_file || !o_isSignal || !o_isMC || !samplename){
            std::cerr<<"Error: info about sample missing!"<<std::endl;
            exit(855);
        }
        const Channel::Channel channelFromFile = Channel::convert(channel_from_file->GetString());
        const Systematic::Systematic systematicFromFile = Systematic::Systematic(systematics_from_file->GetString());

        // Configure information about samples
        const bool isTopSignal = o_isSignal->GetString() == "1";
        const bool isMC = o_isMC->GetString() == "1";
        const bool isHiggsSignal(o_isHiggsSignal && o_isHiggsSignal->GetString()=="1");
        const bool isTtbarV(samplename->GetString()=="ttbarw" || samplename->GetString()=="ttbarz");
        const bool isTtbarbg(samplename->GetString()=="ttbarbg");
        const bool isSingletop(samplename->GetString()=="singletop");
        const bool isDY(samplename->GetString()=="dy50inf" || samplename->GetString()=="dy1050");

        // Checks avoiding running on ill-defined configurations
        if(!isMC && systematic.type()!=Systematic::undefinedType){
            std::cout<<"Sample is DATA, so not running again for systematic variation\n";
            continue;
        }

        bool withPDFWeights = isTopSignal || isTtbarbg;
        bool withPSWeights = isTopSignal || isTtbarbg || isSingletop;
        bool withMEWeights = isTopSignal || isTtbarbg || isSingletop;
        bool withBFragWeights = isTopSignal || isTtbarbg;
        int nPDFWeights = -1;
        int nPSWeights = -1;
        if(isMiniTreeAnalysis||isPlainTreeAnalysis2D){
            if(withPSWeights || withPDFWeights)
            {
                TH1* pdfWeights = dynamic_cast<TH1*>(file->Get("EventsBeforeSelection/pdfEventWeights"));
                if(!pdfWeights){
                    std::cerr<<"ERROR in load_Analysis()! Cannot find histogram pdfEventWeights\n...break\n"<<std::endl;
                    exit(831);
                }
                TH1* psWeights = dynamic_cast<TH1*>(file->Get("EventsBeforeSelection/psEventWeights"));
                if(!psWeights){
                    std::cerr<<"ERROR in load_Analysis()! Cannot find histogram psEventWeights\n...break\n"<<std::endl;
                    exit(831);
                }
                nPDFWeights = utils::helper::GetNumberOfNonZeroBins(*pdfWeights);
                nPSWeights = utils::helper::GetNumberOfNonZeroBins(*psWeights);
                selector->SetReadInAllSyst(true, withPDFWeights, withMEWeights, withBFragWeights, withPSWeights,
                                           utils::helper::GetNumberOfNonZeroBins(*pdfWeights),
                                           utils::helper::GetNumberOfNonZeroBins(*psWeights), true);
            }else{
                selector->SetReadInAllSyst(true, withPDFWeights, withMEWeights, withBFragWeights, withPSWeights, -1, -1, true);
            }
        }

	//Check that histograms for pdf and PS weights exists and that systematic ID is ok
	if(systematic.type() == Systematic::pdf){
            if(!(isTopSignal || (isTtbarbg && useTopMassSetup)) || isHiggsSignal || isTtbarV || systematicFromFile.type() != Systematic::nominal ||
               !(isTopSignal || (isTtbarbg && isPlainTreeAnalysis2D)) || !(isTopSignal || (isTtbarbg && isMiniTreeAnalysis))){
                std::cout<<"Sample is not ttbar dilepton or not nominal, so not running PDF variation\n";
                continue;
            }
            TH1* pdfWeights = dynamic_cast<TH1*>(file->Get("EventsBeforeSelection/pdfEventWeights"));
            if(!pdfWeights){
                std::cerr<<"ERROR in load_Analysis()! Cannot find histogram pdfEventWeights\n...break\n"<<std::endl;
                exit(831);
            }
            if(specific_PDF >= (int) utils::helper::GetNumberOfNonZeroBins(*pdfWeights) || pdfWeights->GetBinContent(specific_PDF+1)<=0.){
                std::cout<<"ID specified for PDF variation is above number of variations stored in sample, so not running\n";
                continue;
            }
        }
        if((systematic.type() == Systematic::psScale) || systematic.type() == Systematic::psISRScale || systematic.type() == Systematic::psFSRScale || systematic.type() == Systematic::psScaleWeight){
            if(systematicFromFile.type() != Systematic::nominal){
                std::cout<<"Sample is not ttbar dilepton or not nominal, so not running PS variation\n";
                continue;
            }
            TH1* psWeights = dynamic_cast<TH1*>(file->Get("EventsBeforeSelection/psEventWeights"));
            if(!psWeights){
                std::cerr<<"ERROR in load_Analysis()! Cannot find histogram psEventWeights\n...break\n"<<std::endl;
                exit(831);
            }
            if(specific_PS>=(int) utils::helper::GetNumberOfNonZeroBins(*psWeights) || psWeights->GetBinContent(specific_PS+1)<=0.){
                std::cout<<"ID specified for PS variation is above number of variations stored in sample, so not running\n";
                continue;
            }
        }

        if(isMC){
        // Set up pileup reweighter
        const std::string pileup_MC_File = (corrections.pileup_MC_File_ == "@{INPUT_FILE}") ? std::string(filename) : common::DATA_PATH_COMMON()+"/"+corrections.pileup_MC_File_;


        const PileupScaleFactors* const pileupScaleFactors = new PileupScaleFactors(corrections.pileup_DATA_File_, corrections.pileup_DATA_sysUP_File_,
                                                                                    corrections.pileup_DATA_sysDN_File_, corrections.pileup_DATA_Histogram_,
                                                                                    pileup_MC_File,corrections.pileup_MC_Histogram_, systematic);


        selector->SetPileupScaleFactors(pileupScaleFactors);
        }

        // JetVetoMaps
	utils::JetVetoMap* jetVetoMap(nullptr);
        if((corrections.jetVetoMap_File_ + corrections.jetVetoMap_Histogram_ + corrections.jetVetoMap_Format_) != ""){
          jetVetoMap = new utils::JetVetoMap("JetVetoMap");
          jetVetoMap->set_verbose(false);
          const std::string jetVetoMap_File = common::DATA_PATH_COMMON()+"/"+corrections.jetVetoMap_File_;
          jetVetoMap->initMap(jetVetoMap_File, corrections.jetVetoMap_Histogram_, corrections.jetVetoMap_Format_);
          selector->SetJetVetoMap(jetVetoMap);
        }

        // Is the channel given in the file? This is true only for data which is preselected due to trigger,
        // and guarantees that only the proper channel is processed (undefined - usual channel for MC).
        // However, for single lepton datasets (se & smu) the analysis should be processed in ee, emu, mumu modes.
        if(channelFromFile != Channel::undefined && channelFromFile != Channel::se && channelFromFile != Channel::smu){
            channels.clear();
            channels.push_back(channelFromFile);
        }

        // Extract JES correction file (via string parsing)
        TObjArray* tokenL1  = TString(corrections.jesL1CorrectionDataFile_).Tokenize(";");
        TObjArray* tokenL2  = TString(corrections.jesL2CorrectionDataFile_).Tokenize(";");
        TObjArray* tokenL3  = TString(corrections.jesL3CorrectionDataFile_).Tokenize(";");
        TObjArray* tokenL2L3  = TString(corrections.jesL2L3CorrectionDataFile_).Tokenize(";");

        // declare JES MC correction file names
        std::string jesL1CorrectionDataFile;
        std::string jesL2CorrectionDataFile;
        std::string jesL3CorrectionDataFile;
        std::string jesL2L3CorrectionDataFile;

        // Check analysis era and that there is more than one file in the list
        if(tokenL1->GetSize() > 1 && tokenL2->GetSize() > 1 && tokenL3->GetSize() > 1
           && tokenL2L3->GetSize() > 1 && (analysisConfig.general().era_ == Era::run2_13tev_2016_25ns || analysisConfig.general().era_ == Era::run2_13tev_2016 || analysisConfig.general().era_ == Era::run2_UL_13tev_2016preVFP || analysisConfig.general().era_ == Era::run2_UL_13tev_2016postVFP)) {

          // Set jes corretion file based on run period
          if(filenameBase.Contains("run2016B.root") || filenameBase.Contains("run2016C.root") || filenameBase.Contains("run2016D.root")) {
            jesL1CorrectionDataFile = (((TObjString*)tokenL1->At(0))->GetString()).Data();
            jesL2CorrectionDataFile = (((TObjString*)tokenL2->At(0))->GetString()).Data();
            jesL3CorrectionDataFile = (((TObjString*)tokenL3->At(0))->GetString()).Data();
            jesL2L3CorrectionDataFile = (((TObjString*)tokenL2L3->At(0))->GetString()).Data();
          }
          else if(filenameBase.Contains("run2016E.root") || filenameBase.Contains("run2016F.root") || filenameBase.Contains("run2016F1.root")) {
            jesL1CorrectionDataFile = (((TObjString*)tokenL1->At(1))->GetString()).Data();
            jesL2CorrectionDataFile = (((TObjString*)tokenL2->At(1))->GetString()).Data();
            jesL3CorrectionDataFile = (((TObjString*)tokenL3->At(1))->GetString()).Data();
            jesL2L3CorrectionDataFile = (((TObjString*)tokenL2L3->At(1))->GetString()).Data();
          }
          else if(filenameBase.Contains("run2016G.root") || filenameBase.Contains("run2016F2.root") || filenameBase.Contains("run2016H.root")) {
            jesL1CorrectionDataFile = (((TObjString*)tokenL1->At(2))->GetString()).Data();
            jesL2CorrectionDataFile = (((TObjString*)tokenL2->At(2))->GetString()).Data();
            jesL3CorrectionDataFile = (((TObjString*)tokenL3->At(2))->GetString()).Data();
            jesL2L3CorrectionDataFile = (((TObjString*)tokenL2L3->At(2))->GetString()).Data();
          }
        }
        else if(tokenL1->GetSize() > 1 && tokenL2->GetSize() > 1 && tokenL3->GetSize() > 1
		&& tokenL2L3->GetSize() > 1 && (analysisConfig.general().era_ == Era::run2_13tev_2017 || analysisConfig.general().era_ == Era::run2_UL_13tev_2017)) {

		  // Set jes corretion file based on run period
		  if(filenameBase.Contains("run2017B.root")) {
			jesL1CorrectionDataFile = (((TObjString*)tokenL1->At(0))->GetString()).Data();
			jesL2CorrectionDataFile = (((TObjString*)tokenL2->At(0))->GetString()).Data();
			jesL3CorrectionDataFile = (((TObjString*)tokenL3->At(0))->GetString()).Data();
			jesL2L3CorrectionDataFile = (((TObjString*)tokenL2L3->At(0))->GetString()).Data();
		  }
		  else if(filenameBase.Contains("run2017C.root")) {
			jesL1CorrectionDataFile = (((TObjString*)tokenL1->At(1))->GetString()).Data();
			jesL2CorrectionDataFile = (((TObjString*)tokenL2->At(1))->GetString()).Data();
			jesL3CorrectionDataFile = (((TObjString*)tokenL3->At(1))->GetString()).Data();
			jesL2L3CorrectionDataFile = (((TObjString*)tokenL2L3->At(1))->GetString()).Data();
		  }
		  else if(filenameBase.Contains("run2017D.root")) {
			jesL1CorrectionDataFile = (((TObjString*)tokenL1->At(2))->GetString()).Data();
			jesL2CorrectionDataFile = (((TObjString*)tokenL2->At(2))->GetString()).Data();
			jesL3CorrectionDataFile = (((TObjString*)tokenL3->At(2))->GetString()).Data();
			jesL2L3CorrectionDataFile = (((TObjString*)tokenL2L3->At(2))->GetString()).Data();
		  }
		  else if(filenameBase.Contains("run2017E.root")) {
			jesL1CorrectionDataFile = (((TObjString*)tokenL1->At(3))->GetString()).Data();
			jesL2CorrectionDataFile = (((TObjString*)tokenL2->At(3))->GetString()).Data();
			jesL3CorrectionDataFile = (((TObjString*)tokenL3->At(3))->GetString()).Data();
			jesL2L3CorrectionDataFile = (((TObjString*)tokenL2L3->At(3))->GetString()).Data();
		  }
		  else if(filenameBase.Contains("run2017F.root")) {
			jesL1CorrectionDataFile = (((TObjString*)tokenL1->At(4))->GetString()).Data();
			jesL2CorrectionDataFile = (((TObjString*)tokenL2->At(4))->GetString()).Data();
			jesL3CorrectionDataFile = (((TObjString*)tokenL3->At(4))->GetString()).Data();
			jesL2L3CorrectionDataFile = (((TObjString*)tokenL2L3->At(4))->GetString()).Data();
		  }
	}
        else if(tokenL1->GetSize() > 1 && tokenL2->GetSize() > 1 && tokenL3->GetSize() > 1
		&& tokenL2L3->GetSize() > 1 && (analysisConfig.general().era_ == Era::run2_13tev_2018 || analysisConfig.general().era_ == Era::run2_UL_13tev_2018)) {

        		  // Set jes corretion file based on run period
        		  if(filenameBase.Contains("run2018A.root")) {
        			jesL1CorrectionDataFile = (((TObjString*)tokenL1->At(0))->GetString()).Data();
        			jesL2CorrectionDataFile = (((TObjString*)tokenL2->At(0))->GetString()).Data();
        			jesL3CorrectionDataFile = (((TObjString*)tokenL3->At(0))->GetString()).Data();
        			jesL2L3CorrectionDataFile = (((TObjString*)tokenL2L3->At(0))->GetString()).Data();
        		  }
        		  else if(filenameBase.Contains("run2018B.root")) {
        			jesL1CorrectionDataFile = (((TObjString*)tokenL1->At(1))->GetString()).Data();
        			jesL2CorrectionDataFile = (((TObjString*)tokenL2->At(1))->GetString()).Data();
        			jesL3CorrectionDataFile = (((TObjString*)tokenL3->At(1))->GetString()).Data();
        			jesL2L3CorrectionDataFile = (((TObjString*)tokenL2L3->At(1))->GetString()).Data();
        		  }
        		  else if(filenameBase.Contains("run2018C.root")) {
        			jesL1CorrectionDataFile = (((TObjString*)tokenL1->At(2))->GetString()).Data();
        			jesL2CorrectionDataFile = (((TObjString*)tokenL2->At(2))->GetString()).Data();
        			jesL3CorrectionDataFile = (((TObjString*)tokenL3->At(2))->GetString()).Data();
        			jesL2L3CorrectionDataFile = (((TObjString*)tokenL2L3->At(2))->GetString()).Data();
        		  }
        		  else if(filenameBase.Contains("run2018D.root")) {
        			jesL1CorrectionDataFile = (((TObjString*)tokenL1->At(3))->GetString()).Data();
        			jesL2CorrectionDataFile = (((TObjString*)tokenL2->At(3))->GetString()).Data();
        			jesL3CorrectionDataFile = (((TObjString*)tokenL3->At(3))->GetString()).Data();
        			jesL2L3CorrectionDataFile = (((TObjString*)tokenL2L3->At(3))->GetString()).Data();
        		  }
        		}
        else {
          // Use single file jes corrections
          jesL1CorrectionDataFile = corrections.jesL1CorrectionDataFile_;
          jesL2CorrectionDataFile = corrections.jesL2CorrectionDataFile_;
          jesL3CorrectionDataFile = corrections.jesL3CorrectionDataFile_;
          jesL2L3CorrectionDataFile = corrections.jesL2L3CorrectionDataFile_;
        }

        // Configure JES correction
        if(analysisConfig.general().era_==Era::run1_8tev) jetEnergyScaleScaleFactors->configure("", "", "", "", isMC);
        else {
            if(!isMC) {
            	jetEnergyScaleScaleFactors->configure(jesL1CorrectionDataFile, jesL2CorrectionDataFile, jesL3CorrectionDataFile, jesL2L3CorrectionDataFile , isMC);
            	jetEnergyScaleScaleFactors_L1only->configure(jesL1CorrectionDataFile);
            }
            else {
            	jetEnergyScaleScaleFactors->configure(corrections.jesL1CorrectionMcFile_, corrections.jesL2CorrectionMcFile_, corrections.jesL3CorrectionMcFile_, "", isMC);
            	jetEnergyScaleScaleFactors_L1only->configure(corrections.jesL1CorrectionMcFile_);
            }
        }


        // If no systematic is specified, read it from the file and use this (used for systematic variations of signal samples, and for nominal)
        const Systematic::Systematic selectedSystematic = systematic.type()==Systematic::undefinedType ? systematicFromFile : systematic;

        // If for specific systematic variations the nominal btagging efficiencies should be used
        const Systematic::Systematic systematicForBtagEfficiencies = (selectedSystematic.type()==Systematic::pdf || selectedSystematic.type()==Systematic::psScaleWeight || selectedSystematic.type()==Systematic::closure) ? Systematic::nominalSystematic() :  selectedSystematic;

        // Set up btag efficiency scale factors
        // This has to be done only after potentially setting systematic from file, since it is varied with signal systematics
        const std::string inputBTagEffFolder = "BTagEff" + yearSuffix;
        const std::string outputBTagEffFolder = outputBaseFolder + "/BTagEff" + yearSuffix;
        BtagScaleFactors* btagScaleFactors = new BtagScaleFactors(inputBTagEffFolder, outputBTagEffFolder,
                                          corrections.btagEfficiencyFilename_, corrections.btagEfficiencyFilename_, std::string(filenameBase),
                                          analysisConfig.general().era_, corrections.btagHeavyFlavourFile_, corrections.btagLightFlavourFile_,
                                          corrections.btagSFInputFile_, channels, systematicForBtagEfficiencies, corrections.btagCorrectionMode_,
                                          analysisConfig.selections().btagAlgorithm_, analysisConfig.selections().btagWorkingPoint_, corrections.btagSFCorrectionViaEff_);

        // Configure selector
        selector->SetTopSignal(isTopSignal);
        selector->SetHiggsSignal(isHiggsSignal);
        selector->SetMC(isMC);
        selector->SetWeightedEvents(weightedEvents);
        selector->SetSamplename(samplename->GetString());
        selector->SetGeneratorBools(samplename->GetString(), selectedSystematic);
        selector->SetSystematic(selectedSystematic);
        selector->SetBtagScaleFactors(btagScaleFactors);

        TriggerScaleFactors* triggerScaleFactors_Up(0);
        TriggerScaleFactors* triggerScaleFactors_Dn(0);

      	KinematicReconstructionScaleFactors* kinematicReconstructionScaleFactors_Up(0);
      	KinematicReconstructionScaleFactors* kinematicReconstructionScaleFactors_Dn(0);

      	LooseKinRecoScaleFactors* looseKinRecoScaleFactors_Up(0);
      	LooseKinRecoScaleFactors* looseKinRecoScaleFactors_Dn(0);

        // BtagScaleFactors* btagScaleFactors_Up(0);
        // BtagScaleFactors* btagScaleFactors_Dn(0);
        // BtagScaleFactors* btagScaleFactors_LjetUp(0);
        // BtagScaleFactors* btagScaleFactors_LjetDn(0);
        std::map<std::tuple<Systematic::Type, Systematic::Variation, int>, BtagScaleFactors*>* mapAllSystBTagScaleFactors(0);

        if(isPlainTreeAnalysis2D||isMiniTreeAnalysis)
        {
          std::cout<<"=============================================================================="<<std::endl;
          std::cout<<"Begin preparation of all instances relevant for plainTree/miniTree production:"<<std::endl;
          Systematic::Systematic tempSys(Systematic::undefinedSystematic());
          if(triggerScaleFactors){
              triggerSF::SFHisto eeSF_(corrections.trigger_DoubleElectron_SF_File_, corrections.trigger_DoubleElectron_SF_Histogram_,
                  corrections.trigger_DoubleElectron_SF_Histogram_Format_, "","");
              triggerSF::SFHisto emuSF_(corrections.trigger_ElectronMuon_SF_File_,
               corrections.trigger_ElectronMuon_SF_Histogram_, corrections.trigger_ElectronMuon_SF_Histogram_Format_, "", "");
              triggerSF::SFHisto mumuSF_(corrections.trigger_DoubleMuon_SF_File_,
               corrections.trigger_DoubleMuon_SF_Histogram_, corrections.trigger_DoubleMuon_SF_Histogram_Format_, "", "");

              tempSys = Systematic::Systematic("TRIG_UP");
              triggerScaleFactors_Up = new TriggerScaleFactors(eeSF_, emuSF_, mumuSF_, corrections.triggerSFReadoutMode_, tempSys, false);
              tempSys = Systematic::Systematic("TRIG_DOWN");
              triggerScaleFactors_Dn = new TriggerScaleFactors(eeSF_, emuSF_, mumuSF_, corrections.triggerSFReadoutMode_, tempSys, false);
          }

      	  tempSys = Systematic::Systematic("KIN_DOWN");
      	  kinematicReconstructionScaleFactors_Up = new KinematicReconstructionScaleFactors(analysisConfig.general().era_, channels, tempSys);
      	  tempSys = Systematic::Systematic("KIN_DOWN");
      	  kinematicReconstructionScaleFactors_Dn = new KinematicReconstructionScaleFactors(analysisConfig.general().era_, channels, tempSys);

      	  tempSys = Systematic::Systematic("KIN_DOWN");
      	  looseKinRecoScaleFactors_Up = new LooseKinRecoScaleFactors(analysisConfig.general().era_, channels, tempSys);
      	  tempSys = Systematic::Systematic("KIN_DOWN");
      	  looseKinRecoScaleFactors_Dn = new LooseKinRecoScaleFactors(analysisConfig.general().era_, channels, tempSys);


          mapAllSystBTagScaleFactors = new std::map<std::tuple<Systematic::Type, Systematic::Variation, int>, BtagScaleFactors*>;
          std::vector<std::string> listOfSysts;
          std::string fixedSysts[30] = {"BTAG_UP", "BTAG_DOWN", "BTAG_LJET_UP", "BTAG_LJET_DOWN", 
					"BTAG_CORR_UP", "BTAG_CORR_DOWN", "BTAG_LJET_CORR_UP", "BTAG_LJET_CORR_DOWN",
					"BTAG_UNCORR_UP", "BTAG_UNCORR_DOWN", "BTAG_LJET_UNCORR_UP", "BTAG_LJET_UNCORR_DOWN",
					"PU_UP", "PU_DOWN",
                                        "TRIG_UP", "TRIG_DOWN", "ELE_ID_UP", "ELE_ID_DOWN", "ELE_RECO_UP", "ELE_RECO_DOWN",
                                        "MUON_ID_UP", "MUON_ID_DOWN", "MUON_ISO_UP", "MUON_ISO_DOWN",
                                        "LEPT_UP","LEPT_DOWN",
                                        "KIN_UP","KIN_DOWN",
                                        "L1PREFIRING_UP", "L1PREFIRING_DOWN"};
          listOfSysts.insert(listOfSysts.end(), std::begin(fixedSysts), std::end(fixedSysts));
          if(withPDFWeights){
              listOfSysts.push_back("PDF_ALPHAS_UP");
              listOfSysts.push_back("PDF_ALPHAS_DOWN");
          }
          if(withMEWeights){
              listOfSysts.push_back("MERENSCALE_UP");
              listOfSysts.push_back("MERENSCALE_DOWN");
              listOfSysts.push_back("MEFACSCALE_UP");
              listOfSysts.push_back("MEFACSCALE_DOWN");
              listOfSysts.push_back("MESCALE_UP");
              listOfSysts.push_back("MESCALE_DOWN");
          }
          if(withBFragWeights){
              listOfSysts.push_back("BFRAG_UP");
              listOfSysts.push_back("BFRAG_DOWN");
              listOfSysts.push_back("BFRAG_PETERSON");
              listOfSysts.push_back("BFRAG_CENTRAL");
              listOfSysts.push_back("BSEMILEP_UP");
              listOfSysts.push_back("BSEMILEP_DOWN");
          }
          for(auto sysName: listOfSysts){
            tempSys = Systematic::Systematic(sysName);
            std::tuple<Systematic::Type, Systematic::Variation, int> tempTuple{tempSys.type(), tempSys.variation(), tempSys.variationNumber()};
            (*mapAllSystBTagScaleFactors)[tempTuple] = new BtagScaleFactors(inputBTagEffFolder, outputBTagEffFolder,
                                                                          corrections.btagEfficiencyFilename_,corrections.btagEfficiencyFilename_, std::string(filenameBase),
                                                                          analysisConfig.general().era_, corrections.btagHeavyFlavourFile_, corrections.btagLightFlavourFile_, corrections.btagSFInputFile_,
                                                                          channels, tempSys, corrections.btagCorrectionMode_,
                                                                          analysisConfig.selections().btagAlgorithm_, analysisConfig.selections().btagWorkingPoint_, corrections.btagSFCorrectionViaEff_, true);
          }
          if(nPDFWeights>0 && withPDFWeights){ // disabled for the moment
            for(int specific_PDF=0; specific_PDF<nPDFWeights; specific_PDF++){
              Systematic::Variation variation = !specific_PDF ? Systematic::central : specific_PDF%2 ? Systematic::up : Systematic::down;
              int variationNumber = (specific_PDF+1)/2;
              tempSys = Systematic::Systematic(Systematic::pdf, variation, variationNumber);
              std::tuple<Systematic::Type, Systematic::Variation, int> tempTuplepdf{tempSys.type(), tempSys.variation(), variationNumber};
              // (*mapAllSystBTagScaleFactors)[tempTuplepdf] = new BtagScaleFactors(inputBTagEffFolder, outputBTagEffFolder,
              //                                                                   corrections.btagEfficiencyFilename_,corrections.btagEfficiencyFilename_, std::string(filenameBase),
              //                                                                   analysisConfig.general().era_, corrections.btagHeavyFlavourFile_, corrections.btagLightFlavourFile_, corrections.btagSFInputFile_,
              //                                                                   channels, tempSys, corrections.btagCorrectionMode_,
              //                                                                   analysisConfig.selections().btagAlgorithm_, analysisConfig.selections().btagWorkingPoint_, corrections.btagSFCorrectionViaEff_, true);
            }
          }
          if(nPSWeights>0 && withPSWeights){
            // PS Weights: Up and Down Variations
            // REF: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopModGen#Event_Generation
            const std::vector<int> v_psScaleUp = {2, 3, 6, 7, 10, 11};
            for(int specific_PS=0; specific_PS<nPSWeights; specific_PS++){
                int variationNumber(-1);
                if(specific_PS < 14){
                    Systematic::Variation variation = (std::find(v_psScaleUp.begin(), v_psScaleUp.end(), specific_PS) != v_psScaleUp.end()) ? Systematic::up : Systematic::down;
                    if(specific_PS == 2 || specific_PS==4){ variationNumber = 2; }
                    else if(specific_PS == 3 || specific_PS==5){ variationNumber = 3; }
                    else if(specific_PS == 6 || specific_PS==8){ variationNumber = 4; }
                    else if(specific_PS == 7 || specific_PS==9){ variationNumber = 5; }
                    else if(specific_PS == 10 || specific_PS==12){ variationNumber = 6; }
                    else if(specific_PS == 11 || specific_PS==13){ variationNumber = 7; }
                    else if(specific_PS == 0){ variationNumber = 0; variation = Systematic::central;}//nominal PS weight
                    else if(specific_PS == 1){ variationNumber = 1; variation = Systematic::central;}//nominal PS weight
                    else {
                        std::cerr<<"ERROR in load_Analysis(), initializing of BTagScaleFactors for all systematics! PS Systematic ID is not valid: "<<specific_PS<<"\n...break\n"<<std::endl;
                        exit(1);
                    }
                    tempSys = Systematic::Systematic(Systematic::psScaleWeight, variation, variationNumber);
                }else{
                    const Systematic::Variation variation = !specific_PS ? Systematic::central : specific_PS%2 ? Systematic::down : Systematic::up;
                    variationNumber = (specific_PS+2)/2;
                    tempSys = Systematic::Systematic(Systematic::psScaleWeight, variation, variationNumber);
                }
              std::tuple<Systematic::Type, Systematic::Variation, int> tempTuplePs{tempSys.type(), tempSys.variation(), variationNumber};
              //consider only 4,5 up and down at the moment (default weights)
              if(variationNumber==4 || variationNumber==5){
                  (*mapAllSystBTagScaleFactors)[tempTuplePs] = new BtagScaleFactors(inputBTagEffFolder, outputBTagEffFolder,
                                                                                corrections.btagEfficiencyFilename_,corrections.btagEfficiencyFilename_, std::string(filenameBase),
                                                                                analysisConfig.general().era_, corrections.btagHeavyFlavourFile_, corrections.btagLightFlavourFile_, corrections.btagSFInputFile_,
                                                                                channels, tempSys, corrections.btagCorrectionMode_,
                                                                                analysisConfig.selections().btagAlgorithm_, analysisConfig.selections().btagWorkingPoint_, corrections.btagSFCorrectionViaEff_, true);
              }
            }
          }

          selector->SetBtagScaleFactors_Var_AllSyst(mapAllSystBTagScaleFactors);

          // tempSys = Systematic::Systematic("BTAG_UP");
          // btagScaleFactors_Up = new BtagScaleFactors(inputBTagEffFolder, outputBTagEffFolder,
          //                                   corrections.btagEfficiencyFilename_,corrections.btagEfficiencyFilename_, std::string(filenameBase),
          //                                   analysisConfig.general().era_, corrections.btagHeavyFlavourFile_, corrections.btagLightFlavourFile_, corrections.btagSFInputFile_,
          //                                   channels, tempSys, corrections.btagCorrectionMode_,
          //                                   analysisConfig.selections().btagAlgorithm_, analysisConfig.selections().btagWorkingPoint_, corrections.btagSFCorrectionViaEff_);
          // tempSys = Systematic::Systematic("BTAG_DOWN");
          // btagScaleFactors_Dn = new BtagScaleFactors(inputBTagEffFolder, outputBTagEffFolder,
          //                                   corrections.btagEfficiencyFilename_,corrections.btagEfficiencyFilename_, std::string(filenameBase),
          //                                   analysisConfig.general().era_, corrections.btagHeavyFlavourFile_, corrections.btagLightFlavourFile_, corrections.btagSFInputFile_,
          //                                   channels, tempSys, corrections.btagCorrectionMode_,
          //                                   analysisConfig.selections().btagAlgorithm_, analysisConfig.selections().btagWorkingPoint_, corrections.btagSFCorrectionViaEff_);
          // tempSys = Systematic::Systematic("BTAG_LJET_UP");
          // btagScaleFactors_LjetUp = new BtagScaleFactors(inputBTagEffFolder, outputBTagEffFolder,
          //                                   corrections.btagEfficiencyFilename_,corrections.btagEfficiencyFilename_, std::string(filenameBase),
          //                                   analysisConfig.general().era_, corrections.btagHeavyFlavourFile_, corrections.btagLightFlavourFile_, corrections.btagSFInputFile_,
          //                                   channels, tempSys, corrections.btagCorrectionMode_,
          //                                   analysisConfig.selections().btagAlgorithm_, analysisConfig.selections().btagWorkingPoint_, corrections.btagSFCorrectionViaEff_);
          // tempSys = Systematic::Systematic("BTAG_LJET_DOWN");
          // btagScaleFactors_LjetDn = new BtagScaleFactors(inputBTagEffFolder, outputBTagEffFolder,
          //                                   corrections.btagEfficiencyFilename_,corrections.btagEfficiencyFilename_, std::string(filenameBase),
          //                                   analysisConfig.general().era_, corrections.btagHeavyFlavourFile_, corrections.btagLightFlavourFile_, corrections.btagSFInputFile_,
          //                                   channels, tempSys, corrections.btagCorrectionMode_,
          //                                   analysisConfig.selections().btagAlgorithm_, analysisConfig.selections().btagWorkingPoint_, corrections.btagSFCorrectionViaEff_);
          if(isMC){
            const std::string pileup_MC_File = (corrections.pileup_MC_File_ == "@{INPUT_FILE}") ? std::string(filename) : common::DATA_PATH_COMMON()+"/"+corrections.pileup_MC_File_;
            tempSys = Systematic::Systematic("PU_UP");
            const PileupScaleFactors* const pileupScaleFactors_Up = new PileupScaleFactors(corrections.pileup_DATA_File_, corrections.pileup_DATA_sysUP_File_,
                                                                                        corrections.pileup_DATA_sysDN_File_, corrections.pileup_DATA_Histogram_,
                                                                                        pileup_MC_File,corrections.pileup_MC_Histogram_, tempSys);
            tempSys = Systematic::Systematic("PU_DOWN");
            const PileupScaleFactors* const pileupScaleFactors_Dn = new PileupScaleFactors(corrections.pileup_DATA_File_, corrections.pileup_DATA_sysUP_File_,
                                                                                        corrections.pileup_DATA_sysDN_File_, corrections.pileup_DATA_Histogram_,
                                                                                        pileup_MC_File,corrections.pileup_MC_Histogram_, tempSys);
            selector->SetPileupScaleFactors_Up(pileupScaleFactors_Up);
            selector->SetPileupScaleFactors_Dn(pileupScaleFactors_Dn);
          }
          if(triggerScaleFactors){
              selector->SetTriggerScaleFactors_Up(triggerScaleFactors_Up);
              selector->SetTriggerScaleFactors_Dn(triggerScaleFactors_Dn);
          }

      	  selector->SetKinematicReconstruction_Up(kinematicReconstruction, kinematicReconstructionScaleFactors_Up);
      	  selector->SetKinematicReconstruction_Dn(kinematicReconstruction, kinematicReconstructionScaleFactors_Dn);

      	  selector->SetLooseKinReco_Up(looseKinReco, looseKinRecoScaleFactors_Up);
      	  selector->SetLooseKinReco_Dn(looseKinReco, looseKinRecoScaleFactors_Dn);

          // selector->SetBtagScaleFactors_Up(btagScaleFactors_Up);
          // selector->SetBtagScaleFactors_Dn(btagScaleFactors_Dn);
          // selector->SetBtagScaleFactors_LjetUp(btagScaleFactors_LjetUp);
          // selector->SetBtagScaleFactors_LjetDn(btagScaleFactors_LjetDn);

          selector->SetRecoTree(true);
          selector->SetSimpleRecoTree(true);
          selector->SetGenTree(isMC);

          selector->SetTTBarSample( (filename.Contains("ttbarsignal")) || (filename.Contains("ttbarbg")) || (filename.Contains("ttbarH125")) );
          bool isNotNominal = !(selectedSystematic.type() == Systematic::nominal || selectedSystematic.type() == Systematic::cp5) || (filename.Contains("ttbar") && filename.Contains("madgraph"))
                              || (filename.Contains("ttbar") && filename.Contains("amcat"));
          selector->SetSystSample(isNotNominal);
          selector->SetSystVarSample(isNotNominal);
        }


        // Configure usage of a special trigger menu, if needed for a particular dataset
        selector->SetSpecialTriggerMenu(filenameBase);

        // Loop over channels and run selector
        for(const auto& selectedChannel : channels){

            // Set the channel
            const TString channelName = Channel::convert(selectedChannel);
            TString outputfilename = filenameBase.BeginsWith(channelName+"_") ? filenameBase : channelName+"_"+filenameBase;

            if(triggerScaleFactors) triggerScaleFactors->prepareChannel(selectedChannel);

            btagScaleFactors->prepareChannel(selectedChannel);


      	    if(kinematicReconstructionScaleFactors) kinematicReconstructionScaleFactors->prepareChannel(selectedChannel);

      	    if(looseKinRecoScaleFactors) looseKinRecoScaleFactors->prepareChannel(selectedChannel);


            selector->SetChannel(selectedChannel);
            if(channelFromFile == Channel::se || channelFromFile == Channel::smu) selector->SetChannelFromFile(channelFromFile); // set channel from file, once SetChannel() is done

            if(isPlainTreeAnalysis2D||isMiniTreeAnalysis){
              if(triggerScaleFactors){
                  triggerScaleFactors_Up->prepareChannel(selectedChannel);
                  triggerScaleFactors_Dn->prepareChannel(selectedChannel);
              }

      	      if(kinematicReconstructionScaleFactors_Up) kinematicReconstructionScaleFactors_Up->prepareChannel(selectedChannel);
      	      if(kinematicReconstructionScaleFactors_Dn) kinematicReconstructionScaleFactors_Dn->prepareChannel(selectedChannel);

      	      if(looseKinRecoScaleFactors_Up) looseKinRecoScaleFactors_Up->prepareChannel(selectedChannel);
      	      if(looseKinRecoScaleFactors_Dn) looseKinRecoScaleFactors_Dn->prepareChannel(selectedChannel);

              // btagScaleFactors_Up->prepareChannel(selectedChannel);
              // btagScaleFactors_Dn->prepareChannel(selectedChannel);
              // btagScaleFactors_LjetUp->prepareChannel(selectedChannel);
              // btagScaleFactors_LjetDn->prepareChannel(selectedChannel);
              for(std::map<std::tuple<Systematic::Type, Systematic::Variation, int>, BtagScaleFactors*>::iterator it = mapAllSystBTagScaleFactors->begin(); it!=mapAllSystBTagScaleFactors->end(); it++){
                  it->second->prepareChannel(selectedChannel);
              }
            }



            // Set up nTuple chain
            TChain chain("writeNTuple/NTuple");
            chain.Add(filename);
            // chain.SetProof(); //will work from 5.34 onwards
            std::cout<<"filename: "<<filename<<" channel: "<<selectedChannel<<std::endl;
            // Split Drell-Yan sample in decay modes ee, mumu, tautau
            selector->SetTrueLevelDYChannel(dy);
            if(dy){
                if(outputfilename.First("_dy") == kNPOS){
                    std::cerr << "DY variations must be run on DY samples!\n";
                    std::cerr << outputfilename << " must start with 'channel_dy'\n";
                    exit(1);
                }

                selector->SetDrellYan(true);
                const Channel::Channel dyChannel = dy == 11 ? Channel::ee : dy == 13 ? Channel::mumu : Channel::tautau;
                outputfilename.ReplaceAll("_dy", TString("_dy").Append(Channel::convert(dyChannel)));
            }
            else{
                selector->SetDrellYan(false);
            }

            // Run the selector
            selector->SetRunViaTau(0);
            selector->SetOutputfilename(outputfilename);
            selector->SetClosureTest(closure, slope);
            if(isTopSignal && !isTtbarV) selector->SetSampleForBtagEfficiencies(true);
            if (runsignalviatau){

                if (isTopSignal && !isTtbarV && closure == ""){
                    selector->SetRunViaTau(1);
                    outputfilename.ReplaceAll("signalplustau", "signalviatau");
                    selector->SetOutputfilename(outputfilename);
                }
                else{
                    std::cerr << "ERROR: inconsistent option --signalviatau\n";
                    std::exit(1);
                }
            }
            // For running on PDF systematics
            if(selectedSystematic.type() == Systematic::pdf){
                selector->SetPdfVariation(specific_PDF);
            }
    	    // For running on systematics from parton shower
    	    else if((selectedSystematic.type() == Systematic::psScale) || systematic.type() == Systematic::psISRScale || systematic.type() == Systematic::psFSRScale || systematic.type() == Systematic::psScaleWeight){
    	      selector->SetPSVariation(specific_PS);
    	    }
            // Reveighting to nominal CT14nlo
            else if(useTopMassSetup && (isTopSignal || isTtbarbg) && !isTtbarV && selectedSystematic.type()!=Systematic::amcatnlofxfx){
                selector->SetPdfVariation(0);
            }
            chain.Process(selector, "", maxEvents, skipEvents);
            selector->SetSampleForBtagEfficiencies(false);
            selector->SetPdfVariation(-1);
            selector->SetPSVariation(-1);

        }
        file->Close();
    }
}



/// All systematics allowed for analysis step
/// Only systematics which run on the nominal ntuples, e.g. pileup reweighting
namespace Systematic{
    const std::vector<Type> allowedSystematics = {
        nominal,
        pu, lept, trig, trigEta,
        ele, muon,
        eleID, eleIDStat, eleIDSyst,
        eleReco, eleRecoStat, eleRecoSyst,
        muonID, muonIDStat, muonIDSyst,
        muonIso, muonIsoStat, muonIsoSyst,
        eleScaleSyst, eleScaleGain, eleScaleStat, eleScaleEt,
        eleSmearingPhi, eleSmearingRho,
        eleScaleSmearing,
        muonScaleStat, muonScaleZpt, muonScaleEwk, muonScaleDeltaM, muonScaleEwk2, muonScale,
        jer, jes,
        jerEta0,jerEta1,
        unclustered,
        btag, btagPt, btagEta,
        btagLjet, btagLjetPt, btagLjetEta,
	btag_corr, btagLjet_corr, btag_uncorr, btagLjet_uncorr,

        kin,
        // lumi,
        cp5,
        topPt,
        match,
        meScale, meFacScale, meRenScale,
	    psISRScale, psFSRScale, psScale,
        alphasPdf, pdf, l1prefiring,
        bFrag, bFrag_central, bFrag_Peterson, bSemilep,
        psScaleWeight,

        mTop166, mTop169, mTop171, mTop173, mTop175, mTop178,

        // b-tag discriminator reweighting variations
        btagDiscrBpurity,
        btagDiscrLpurity,
        btagDiscrBstat1,
        btagDiscrBstat2,
        btagDiscrLstat1,
        btagDiscrLstat2,
        btagDiscrCerr1,
        btagDiscrCerr2,

        // independent variation of JES sources
        jesAbsoluteStat,
        jesAbsoluteScale,
        jesAbsoluteFlavMap,
        jesAbsoluteMPFBias,
        jesFragmentation,
        jesSinglePionECAL,
        jesSinglePionHCAL,
        jesFlavorQCD,
        jesFlavorZJet,
        jesFlavorRealistic,
        jesFlavorPhotonJet,
        jesFlavorPureGluon,
        jesFlavorPureQuark,
        jesFlavorPureCharm,
        jesFlavorPureBottom,
        jesTimePtEta,
        jesRelativeBal,
        jesRelativeJEREC1,
        jesRelativeJEREC2,
        jesRelativeJERHF,
        jesRelativePtBB,
        jesRelativePtEC1,
        jesRelativePtEC2,
        jesRelativePtHF,
        jesRelativeFSR,
        jesRelativeStatFSR,
        jesRelativeStatEC,
        jesRelativeStatHF,
        jesPileUpDataMC,
        jesPileUpPtRef,
        jesPileUpPtBB,
        jesPileUpPtEC1,
        jesPileUpPtEC2,
        jesPileUpPtHF,
        jesRelativeSample,
        jesUserDefinedHEM1516,

        //jes combined
        jesAbsolute,
        jesAbsoluteYear,
        jesBBEC1,
        jesBBEC1Year,
        jesEC2,
        jesEC2Year,
        jesHF,
        jesHFYear,

        closure

    };
}



int main(int argc, char** argv) {
    CLParameter<std::string> opt_f("f", "Restrict to filename pattern, e.g. ttbar", false, 1, 1);
    CLParameter<std::string> opt_c("c", "Specify a certain channel (ee, emu, mumu). No channel specified = run on all channels", false, 1, 1,
            common::makeStringCheck(Channel::convert(Channel::allowedChannelsAnalysis)));
    CLParameter<std::string> opt_s("s", "Run with a systematic that runs on the nominal ntuples, e.g. 'PDF', 'PU_UP' or 'TRIG_DOWN'", false, 1, 1,
            common::makeStringCheckBegin(Systematic::convertType(Systematic::allowedSystematics)));
    CLParameter<std::string> opt_mode("m", "Mode of analysis: control plots (cp), "
				      "double differential analysis (dda), "
				      "gen level plots (gen), "
				      "kin. reco. efficiency plots (kinReco), loose kin reco plots (looseKinReco), "
				      "kin reco plots for quality studies (kinRecoQualityStudies), "
				      "tree for 2d unfolding (ddaTree), plain tree (plainTree), flat mini trees (miniTree),"
                                      "Default is cp, dda, kinReco", false, 1, 100,
            common::makeStringCheck(AnalysisMode::convert(AnalysisMode::allowedAnalysisModes)));
    CLParameter<int> opt_pdfno("pdf", "Run a certain PDF systematic only, sets -s PDF. Use e.g. --pdf n, where n=0 is central, 1=variation 1 up, 2=1down, 3=2up, 4=2down, ...",
			       false, 1, 1,
	    [](int id){return id>=0;});
    CLParameter<int> opt_psno("ps", "ID of systematic variation for systematics requiring one, e.g. 'PS'", false, 1, 1,
	    [](int id){return id >= 0;});
    CLParameter<int> opt_dy("d", "Drell-Yan mode (11 for ee, 13 for mumu, 15 for tautau)", false, 1, 1,
            [](int dy){return dy == 11 || dy == 13 || dy == 15;});
    CLParameter<std::string> opt_closure("closure", "Enable the closure test. Valid: pttop|ytop|nominal", false, 1, 1,
            [](const std::string &c){return c == "pttop" || c == "ytop" || c == "nominal";});
    CLParameter<float> opt_closureSlope("slope", "Slope for closure test, use -0.01 to 0.01 for pt and -0.4 to 0.4 for ytop", false, 1, 1,
            [](float s){return std::abs(s) < 1;});
    CLParameter<Long64_t> opt_maxEvents("maxEvents", "Maximum number of events to process", false, 1, 1,
            [](const Long64_t mE){return mE > 0;});
    CLParameter<Long64_t> opt_skipEvents("skipEvents", "Number of events to be skipped", false, 1, 1,
            [](const Long64_t sE){return sE > 0;});
    CLParameter<std::string> opt_runsignalviatau("signalviatau", "Run over signal via tau", false, 0, 0);
    CLParameter<bool> opt_mass("topmass", "Use setup for top mass extraction analysis ", false, 0, 0);
    CLParameter<bool> opt_dryRun("dry", "Dry run, only print filenames and steering commands to text file", false, 0, 0);
    CLAnalyser::interpretGlobal(argc, argv);

    std::cout<<"\n"<<"--- Beginning setting up command line steering parameters\n";

    // Create command line string without dry mode
    std::string commandline;
    for(int iArg = 0; iArg < argc; ++iArg){
        std::stringstream ss_arg;
        ss_arg<<argv[iArg];
        const std::string arg = ss_arg.str();
        if(arg == "--dry") continue;
        commandline.append(arg).append(" ");
    }

    std::cout<<commandline<<std::endl;

    // // Set up dry run if chosen
    const bool dryRun = opt_dryRun.isSet();

    TString validFilenamePattern = opt_f.isSet() ? opt_f[0] : "";

    // Set up channel
    Channel::Channel channel(Channel::undefined);
    if(opt_c.isSet()) channel = Channel::convert(opt_c[0]);

    // Set up systematic
    Systematic::Systematic systematic(Systematic::undefinedSystematic());
    if(opt_s.isSet()) systematic = Systematic::Systematic(opt_s[0]);

    // Set up analysis mode
    std::vector<AnalysisMode::AnalysisMode> v_analysisMode({AnalysisMode::cp});

    if(opt_mode.isSet()) v_analysisMode = AnalysisMode::convert(opt_mode.getArguments());
    std::cout<<"\nRunning the following analysis modes:\n";
    for(const auto& analysisMode : v_analysisMode) std::cout<<AnalysisMode::convert(analysisMode)<<" , ";
    std::cout<<"\n\n";

    int dy = opt_dy.isSet() ? opt_dy[0] : 0;
    TString closure = opt_closure.isSet() ? opt_closure[0] : "";
    float slope = 0;
    if (closure != "" && closure != "nominal") {
        if (!opt_closureSlope.isSet()) {
            std::cerr << "closure test: need slope!\n"; exit(1);
        } else {
            slope = opt_closureSlope[0];
        }
    }


    // Set up systematic ID for systematics requiring one (e.g. PDF)
    int pdf_no(-1);
    int ps_no(-1);
    int variationNumber(-1);

    if(opt_pdfno.isSet()){
        pdf_no = opt_pdfno[0];
        if(systematic.type() == Systematic::pdf){
            const Systematic::Variation variation = !pdf_no ? Systematic::central : pdf_no%2 ? Systematic::up : Systematic::down;
            variationNumber = (pdf_no+1)/2;
            systematic = Systematic::Systematic(Systematic::pdf, variation, variationNumber);
        }
    }else{
        if(systematic.type() == Systematic::pdf){
            std::cerr<<"ERROR in load_Analysis executable! Systematic requires systematic ID, but none specified\n...break\n"<<std::endl;
            exit(1);
       }
   }
    // PS Weights: Up and Down Variations
    // REF: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopModGen#Event_Generation
    const std::vector<int> v_psScaleUp = {2, 3, 6, 7, 10, 11};

    if(opt_psno.isSet()){
        ps_no = opt_psno[0];
        if(systematic.type() == Systematic::psScale){
            if((ps_no < 2) || (ps_no > 14)){
                std::cerr << "ERROR in load_Analysis executable! Systematic ID not in allowed PS systematics\n...break\n" << std::endl;
                exit(1);
            }
            const Systematic::Variation variation = (std::find(v_psScaleUp.begin(), v_psScaleUp.end(), ps_no) != v_psScaleUp.end()) ? Systematic::up : Systematic::down;
            if(ps_no <  6){ variationNumber = 1; } // Reduced
            else if(ps_no < 10){ variationNumber = 2; } // Default
            else { variationNumber = 3; } // Conservative

            if (ps_no%2==0){ systematic = Systematic::Systematic(Systematic::psISRScale, variation, variationNumber); }
            else { systematic = Systematic::Systematic(Systematic::psFSRScale, variation, variationNumber); }
        }else if(systematic.type() == Systematic::psScaleWeight){
            if(ps_no < 14){
                Systematic::Variation variation = (std::find(v_psScaleUp.begin(), v_psScaleUp.end(), ps_no) != v_psScaleUp.end()) ? Systematic::up : Systematic::down;
                if(ps_no == 2 || ps_no==4){ variationNumber = 2; }
                else if(ps_no == 3 || ps_no==5){ variationNumber = 3; }
                else if(ps_no == 6 || ps_no==8){ variationNumber = 4; }
                else if(ps_no == 7 || ps_no==9){ variationNumber = 5; }
                else if(ps_no == 10 || ps_no==12){ variationNumber = 6; }
                else if(ps_no == 11 || ps_no==13){ variationNumber = 7; }
                else if(ps_no == 0){ variationNumber = 0; variation=Systematic::central;}//nominal PS weight
                else if(ps_no == 1){ variationNumber = 1; variation=Systematic::central;}//nominal PS weight
                else {
                    std::cerr<<"ERROR in load_Analysis executable! PS Systematic ID is not valid: "<<ps_no<<"\n...break\n"<<std::endl;
                    exit(1);
                }
                systematic = Systematic::Systematic(Systematic::psScaleWeight, variation, variationNumber);
            }else{
                const Systematic::Variation variation = !ps_no ? Systematic::central : ps_no%2 ? Systematic::down : Systematic::up;
                variationNumber = (ps_no+2)/2;
                systematic = Systematic::Systematic(Systematic::psScaleWeight, variation, variationNumber);
            }
        }else{
            std::cerr<<"ERROR in load_Analysis executable! Systematic is not of type requiring systematic ID, but one specified: "<<ps_no<<"\n...break\n"<<std::endl;
            exit(1);
        }
        std::cout<<"Systematic variation constructed from ID (ID, name): "<<ps_no<<" , "<<systematic.name()<<"\n";
    }else{
        if(systematic.type() == Systematic::psScale){
            std::cerr<<"ERROR in load_Analysis executable! Systematic requires systematic ID, but none specified\n...break\n"<<std::endl;
            exit(1);
        }
    }

    // Set up maximum number of events to process
    const Long64_t bigNumber(TChain::kBigNumber);
    const Long64_t maxEvents = dryRun ? 0 :
      opt_maxEvents.isSet() ? opt_maxEvents[0]
      : bigNumber;

    // Set up number of events to be skipped
    const Long64_t skipEvents = opt_skipEvents.isSet() ? opt_skipEvents[0] : 0;

    //     TProof* p = TProof::Open(""); // not before ROOT 5.34
    //    load_Analysis(validFilenamePattern, channel, systematic, v_analysisMode, pdf_no, ps_no, dy, closure, slope, maxEvents, skipEvents, opt_runsignalviatau.isSet(), opt_mass.isSet(),
    //                  dryRun, commandline);

    load_Analysis(validFilenamePattern, channel, systematic, v_analysisMode, pdf_no, ps_no, dy, closure, slope, maxEvents, skipEvents, opt_runsignalviatau.isSet(), opt_mass.isSet(),dryRun, commandline);

    //     delete p;
}
