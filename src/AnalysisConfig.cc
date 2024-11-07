#include <iostream>
#include <sstream>

#include "AnalysisConfig.h"
#include "../../common/include/TextfileReader.h"
#include "../../common/include/sampleHelpers.h"





AnalysisConfig::General::General():
era_(Era::undefined),
luminosity_(-999.),
luminosityUncertainty_(-999.)
{}



std::string AnalysisConfig::General::print(const bool screen)const
{
    std::ostringstream ss;
    ss<<"Beginning printing AnalysisConfig general\n"
      <<"\tera: "<<Era::convert(era_)<<"\n"
      <<"\tluminosity: "<<luminosity_<<"\n"
      <<"\tluminosity uncertainty: "<<luminosityUncertainty_<<"\n"
      <<"Finishing printing AnalysisConfig general\n\n";
    const std::string s = ss.str();
    if(screen) std::cout<<s;
    return s;
}



AnalysisConfig::Corrections::Corrections():
pileup_MC_File_(""),
pileup_MC_Histogram_(""),
pileup_DATA_File_(""),
pileup_DATA_Histogram_(""),
pileup_DATA_sysUP_File_(""),
pileup_DATA_sysDN_File_(""),
triggerSFInputSuffix_(""),
triggerSFReadoutMode_(0),
trigger_DoubleElectron_SF_File_(""),
trigger_DoubleElectron_SF_Histogram_(""),
trigger_DoubleElectron_SF_Histogram_Format_(""),
trigger_ElectronMuon_SF_File_(""),
trigger_ElectronMuon_SF_Histogram_(""),
trigger_ElectronMuon_SF_Histogram_Format_(""),
trigger_DoubleMuon_SF_File_(""),
trigger_DoubleMuon_SF_Histogram_(""),
trigger_DoubleMuon_SF_Histogram_Format_(""),
electronSFInputFile_(""),
muonSFInputFile_(""),
leptonSFReadoutMode_(0),
combinedSF_(true),
SF_electron_ID_File_(""),
SF_electron_ID_Histo_(""),
SF_electron_ID_Format_(""),
SF_electron_Reco_File_(""),
SF_electron_Reco_Histo_(""),
SF_electron_Reco_Format_(""),
SF_electron_Reco_File_LowPt_(""),
SF_electron_Reco_Histo_LowPt_(""),
SF_electron_Reco_Format_LowPt_(""),
SF_muon_ID_File_(""),
SF_muon_ID_Histo_(""),
SF_muon_ID_Histo_Sys_(""),
SF_muon_ID_Histo_Stat_(""),
SF_muon_ID_Format_(""),
SF_muon_ISO_File_(""),
SF_muon_ISO_Histo_(""),
SF_muon_ISO_Histo_Sys_(""),
SF_muon_ISO_Histo_Stat_(""),
SF_muon_ISO_Format_(""),
extra_muon_ISO_DY_extrapolation_unc_(0.),
extra_electron_ID_DY_extrapolation_unc_(0.),
jetVetoMap_File_(""),
jetVetoMap_Histogram_(""),
jetVetoMap_Format_(""),
applyMuonEnergyCorrections_    (false),
applyElectronEnergyCorrections_(false),
applyBJetEnergyRegression_(false),
propagateCorrectionsToMET_(false),
jetCorrectionMode_(0),
jerUncertaintySourceName_(""),
jerSFBinning_(""),
jet_JERHybridMethod_PtResolutionFile_(""),
jet_JERHybridMethod_jetConeRadius_(0.),
jesUncertaintySourceFile_(""),
jesUncertaintyReducedSourceFile_(""),
jesL1CorrectionDataFile_(""),
jesL2CorrectionDataFile_(""),
jesL3CorrectionDataFile_(""),
jesL2L3CorrectionDataFile_(""),
jesL1CorrectionMcFile_(""),
jesL2CorrectionMcFile_(""),
jesL3CorrectionMcFile_(""),
useReducedJESSet_(false),
btagEfficiencyFilename_(""),
btagCorrectionMode_(Btag::undefinedCorrectionMode),
btagHeavyFlavourFile_(""),
btagLightFlavourFile_(""),
btagSFInputFile_(""),
btagSFCorrectionViaEff_(false),
topPtReweighting_(false),
applyTopSlope2DGenLevelReweighting_(false)
{}



std::string AnalysisConfig::Corrections::print(const bool screen)const
{
    std::ostringstream ss;
    ss<<"Beginning printing AnalysisConfig corrections\n"
      <<"\tpileup data file: "<<pileup_DATA_File_<<"\n"
      <<"\tpileup data histogram name: "<<pileup_DATA_Histogram_<<"\n"
      <<"\tpileup data file sysUp: "<<pileup_DATA_sysUP_File_<<"\n"
      <<"\tpileup data file sysDn: "<<pileup_DATA_sysDN_File_<<"\n"
      <<"\tpileup MC file: "<<pileup_MC_File_<<"\n"
      <<"\tpileup MC histogram name: "<<pileup_MC_Histogram_<<"\n"
      <<"\ttrigger SF suffix: "<<triggerSFInputSuffix_<<"\n"
      <<"\ttrigger SF readout mode: "<<triggerSFReadoutMode_<<"\n"
      <<"\ttrigger DoubleElectron SF File: "<<trigger_DoubleElectron_SF_File_<<"\n"
      <<"\ttrigger DoubleElectron SF Histogram: "<<trigger_DoubleElectron_SF_Histogram_<<"\n"
      <<"\ttrigger DoubleElectron SF Histogram_Format: "<<trigger_DoubleElectron_SF_Histogram_Format_<<"\n"
      <<"\ttrigger ElectronMuon SF File: "<<trigger_ElectronMuon_SF_File_<<"\n"
      <<"\ttrigger ElectronMuon SF Histogram: "<<trigger_ElectronMuon_SF_Histogram_<<"\n"
      <<"\ttrigger ElectronMuon SF Histogram_Format: "<<trigger_ElectronMuon_SF_Histogram_Format_<<"\n"
      <<"\ttrigger DoubleMuon SF File: "<<trigger_DoubleMuon_SF_File_<<"\n"
      <<"\ttrigger DoubleMuon SF Histogram: "<<trigger_DoubleMuon_SF_Histogram_<<"\n"
      <<"\ttrigger DoubleMuon SF Histogram_Format: "<<trigger_DoubleMuon_SF_Histogram_Format_<<"\n"
      <<"\telectron SF file: "<<electronSFInputFile_<<"\n"
      <<"\tmuon SF file: "<<muonSFInputFile_<<"\n"
      <<"\tlepton SF readout mode: "<<leptonSFReadoutMode_<<"\n"
      <<"\tcombined lepton SF: "<<combinedSF_<<"\n"
      <<"\telectron ID SF file: "<<SF_electron_ID_File_<<"\n"
      <<"\telectron ID SF histogram name: "<<SF_electron_ID_Histo_<<"\n"
      <<"\telectron ID SF readout mode: "<<SF_electron_ID_Format_<<"\n"
      <<"\telectron Reconstruction Eff. SF file: "<<SF_electron_Reco_File_<<"\n"
      <<"\telectron Reconstruction Eff. SF histogram name: "<<SF_electron_Reco_Histo_<<"\n"
      <<"\telectron Reconstruction Eff. SF readout mode: "<<SF_electron_Reco_Format_<<"\n"
      <<"\telectron Reconstruction Eff. SF LowPt file: "<<SF_electron_Reco_File_LowPt_<<"\n"
      <<"\telectron Reconstruction Eff. SF LowPt histogram name: "<<SF_electron_Reco_Histo_LowPt_<<"\n"
      <<"\telectron Reconstruction Eff. SF LowPt readout mode: "<<SF_electron_Reco_Format_LowPt_<<"\n"
      <<"\telectron ID extra unc (DY extrapolation): "<<extra_electron_ID_DY_extrapolation_unc_<<"\n"
      <<"\tmuon ID SF file: "<<SF_muon_ID_File_<<"\n"
      <<"\tmuon ID SF histogram name: "<<SF_muon_ID_Histo_<<"\n"
      <<"\tmuon ID SF histogram name statErr: "<<SF_muon_ID_Histo_Stat_<<"\n"
      <<"\tmuon ID SF histogram name sysErr: "<<SF_muon_ID_Histo_Sys_<<"\n"
      <<"\tmuon ID SF readout mode: "<<SF_muon_ID_Format_<<"\n"
      <<"\tmuon ISO SF file: "<<SF_muon_ISO_File_<<"\n"
      <<"\tmuon ISO SF histogram name: "<<SF_muon_ISO_Histo_<<"\n"
      <<"\tmuon ISO SF histogram name statErr: "<<SF_muon_ISO_Histo_Stat_<<"\n"
      <<"\tmuon ISO SF histogram name sysErr: "<<SF_muon_ISO_Histo_Sys_<<"\n"
      <<"\tmuon ISO SF readout mode: "<<SF_muon_ISO_Format_<<"\n"
      <<"\tmuon ISO extra unc (DY extrapolation): "<<extra_muon_ISO_DY_extrapolation_unc_ <<"\n"
      << "\tJet VetoMap file: " << jetVetoMap_File_ << "\n"
      << "\tJet VetoMap histogram: " << jetVetoMap_Histogram_ << "\n"
      << "\tJet VetoMap format: " << jetVetoMap_Format_ << "\n"
      << "\tApply Muon Energy Corrections: "     << applyMuonEnergyCorrections_<<"\n"
      << "\n\tApply Electron Energy Corrections: " << applyElectronEnergyCorrections_<<"\n"
      << "\n\tApply BJet Energy Regression: " << applyBJetEnergyRegression_<<"\n"
      << "\n\tPropagate lepton and jet corrections to MET: " << propagateCorrectionsToMET_<<"\n"
      <<"\tjet correction mode: "<<jetCorrectionMode_<<"\n"
      <<"\tJER name: "<<jerUncertaintySourceName_<<"\n"
      <<"\tJER SF binning: "<<jerSFBinning_<<"\n"
      <<"\tJES uncertainty file: "<<jesUncertaintySourceFile_<<"\n"
      <<"\tReduced JES uncertainty file: "<<jesUncertaintyReducedSourceFile_<<"\n"
      <<"\tJet PtResolution file (JER Hybrid Method): "<<jet_JERHybridMethod_PtResolutionFile_<<"\n"
      <<"\tJet Cone Radius       (JER Hybrid Method): "<<jet_JERHybridMethod_jetConeRadius_<<"\n"
      <<"\tJES L1 Data SF file: "<<jesL1CorrectionDataFile_<<"\n"
      <<"\tJES L2 Data SF file: "<<jesL2CorrectionDataFile_<<"\n"
      <<"\tJES L3 Data SF file: "<<jesL3CorrectionDataFile_<<"\n"
      <<"\tJES L2L3 Data SF file: "<<jesL2L3CorrectionDataFile_<<"\n"
      <<"\tJES L1 MC SF file: "<<jesL1CorrectionMcFile_<<"\n"
      <<"\tJES L2 MC SF file: "<<jesL2CorrectionMcFile_<<"\n"
      <<"\tJES L3 MC SF file: "<<jesL3CorrectionMcFile_<<"\n"
      <<"\tb-tag efficiency filename: "<<btagEfficiencyFilename_<<"\n"
      <<"\tb-tag correction mode: "<<Btag::convertCorrectionMode(btagCorrectionMode_)<<"\n"
      <<"\tb-tag discriminator reweight HF file: "<<btagHeavyFlavourFile_<<"\n"
      <<"\tb-tag discriminator reweight LF file: "<<btagLightFlavourFile_<<"\n"
      <<"\tb-tag SF file: "<<btagSFInputFile_<<"\n"
      <<"\tb-tag scale factors from file correction requested : "<<btagSFCorrectionViaEff_<<"\n"
      <<"\tTop pt reweighting applied: "<<topPtReweighting_<<"\n"
      <<"\tTop pt 2d reweighting applied: "<<applyTopSlope2DGenLevelReweighting_<<"\n"
      <<"Finishing printing AnalysisConfig corrections\n\n";
    const std::string s = ss.str();
    if(screen) std::cout<<s;
    return s;
}



AnalysisConfig::Selections::Selections():
leptonEtaCut_(-999.),
leptonPtCut_(-999.),
leadingLeptonPtCut_(-999.),
muonIsoCut_(-999.),
jetEtaCut_(-999.),
jetPtCut_(-999.),
jetPFIDWorkingPoint_(-1.),
jetPileupIdWorkingPoint_(JetPileupId::undefinedWP),
jetPileupIdMinPt_(-1),
jetPileupIdMaxPt_(-1),
jetVetoMapCutValue_(999.),
deltaRLeptonJetCut_(-999.),
lead2JetPtCut_(-999.),
btagAlgorithm_(Btag::undefinedAlgorithm),
btagWorkingPoint_(Btag::undefinedWP),
metAlgorithm_(Met::undefinedAlgorithm),
metCut_(-999.),
genJetEtaCut_(-999.),
genJetPtCut_(-999.),
genDeltaRLeptonJetCut_(-999.),
vetoEventsWithMoreThan2Leptons_(0)
{}



std::string AnalysisConfig::Selections::print(const bool screen)const
{
    std::ostringstream ss;
    ss<<"Beginning printing AnalysisConfig selections\n"
      <<"\tlepton eta: "<<leptonEtaCut_<<"\n"
      <<"\tlepton pt: "<<leptonPtCut_<<"\n"
      <<"\tleading lepton pt: "<<leadingLeptonPtCut_<<"\n"
      <<"\tmuon iso: "<<muonIsoCut_<<"\n"
      <<"\tjet eta: "<<jetEtaCut_<<"\n"
      <<"\tjet pt: "<<jetPtCut_<<"\n"
      <<"\tjet PF-ID working point: " << jetPFIDWorkingPoint_ << "\n"
      <<"\tjet pileup ID working point: "<<JetPileupId::convertWorkingPoint(jetPileupIdWorkingPoint_)<<"\n"
      <<"\tjet pileup ID pT-range of applicability (min=" << jetPileupIdMinPt_ << ", max=" << jetPileupIdMaxPt_ << ")\n"
      <<"\tjet veto map cut value =" << jetVetoMapCutValue_ <<"\n"
      <<"\tdeltaR(lepton, jet): "<<deltaRLeptonJetCut_<<"\n"
      <<"\tleading 2 jet pt: "<<lead2JetPtCut_<<"\n"
      <<"\tb-tag algorithm: "<<Btag::convertAlgorithm(btagAlgorithm_)<<"\n"
      <<"\tb-tag working point: "<<Btag::convertWorkingPoint(btagWorkingPoint_)<<"\n"
      <<"\tMET algorithm: "<<Met::convertAlgorithm(metAlgorithm_)<<"\n"
      <<"\tMET et: "<<metCut_<<"\n"
      <<"\tgenJet eta: "<<genJetEtaCut_<<"\n"
      <<"\tgenJet pt: "<<genJetPtCut_<<"\n"
      <<"\tdeltaR(genJet, genLepton): "<<genDeltaRLeptonJetCut_<<"\n"
      <<"\tveto events with > 2 selected leptons: "<<vetoEventsWithMoreThan2Leptons_<<"\n"
      <<"Finishing printing AnalysisConfig selections\n\n";
    const std::string s = ss.str();
    if(screen) std::cout<<s;
    return s;
}



AnalysisConfig::SampleComposition::SampleComposition():
addSingleLeptonTriggerInfo_(false),
pseudodata_(0),
mergeLevel_(0),
pseudoTopMode_(0)
{}



std::string AnalysisConfig::SampleComposition::print(const bool screen)const
{
    std::ostringstream ss;
    ss<<"Beginning printing AnalysisConfig sample composition\n"
      <<"\tAddition of single lepton trigger info: "<<addSingleLeptonTriggerInfo_<<"\n"
      <<"\tPseudodata: "<<pseudodata_<<"\n"
      <<"\tMerge level: "<<mergeLevel_<<"\n"
      <<"\tPseudo-top mode: "<<pseudoTopMode_<<"\n"
      <<"Finishing printing AnalysisConfig sample composition\n\n";
    const std::string s = ss.str();
    if(screen) std::cout<<s;
    return s;
}



AnalysisConfig::PlotStyle::PlotStyle():
cmsLabel_(0)
{}



std::string AnalysisConfig::PlotStyle::print(const bool screen)const
{
    std::ostringstream ss;
    ss<<"Beginning printing AnalysisConfig plot style\n"
      <<"\tCMS label: "<<cmsLabel_<<"\n"
      <<"Finishing printing AnalysisConfig plot style\n\n";
    const std::string s = ss.str();
    if(screen) std::cout<<s;
    return s;
}



AnalysisConfig::AnalysisConfig(const std::string& configfilename)
{
    std::cout<<"--- Beginning setting up analysis config\n";

    TextfileReader textfileReader;

    // Read general info
    textfileReader.setStartMarker("[ general ]");
    textfileReader.setEndMarker("[ end - general ]");
    textfileReader.readFile(configfilename);
    general_.era_ = Era::convert(textfileReader.getValue<std::string>("era"));
    general_.luminosity_ = textfileReader.getValue<float>("luminosity");
    general_.luminosityUncertainty_ = textfileReader.getValue<float>("luminosityUncertainty");
    textfileReader.clear();

    // Read corrections
    textfileReader.setStartMarker("[ corrections ]");
    textfileReader.setEndMarker("[ end - corrections ]");
    textfileReader.readFile(configfilename);
    textfileReader.setRequireValues(false);
    corrections_.pileup_MC_File_ = textfileReader.getValue<std::string>("pileup_MC_File", "");
    corrections_.pileup_MC_Histogram_ = textfileReader.getValue<std::string>("pileup_MC_Histogram", "");
    corrections_.pileup_DATA_File_ = textfileReader.getValue<std::string>("pileup_DATA_File", "");
    corrections_.pileup_DATA_Histogram_ = textfileReader.getValue<std::string>("pileup_DATA_Histogram", "");
    corrections_.pileup_DATA_sysUP_File_ = textfileReader.getValue<std::string>("pileup_DATA_sysUP_File", "");
    corrections_.pileup_DATA_sysDN_File_ = textfileReader.getValue<std::string>("pileup_DATA_sysDN_File", "");
    corrections_.triggerSFInputSuffix_ = textfileReader.getValue<std::string>("triggerSFInputSuffix", "");
    corrections_.triggerSFReadoutMode_ = textfileReader.getValue<int>("triggerSFReadoutMode", 0);
    corrections_.trigger_DoubleElectron_SF_File_ = textfileReader.getValue<std::string>("trigger_DoubleElectron_SF_File", "");
    corrections_.trigger_DoubleElectron_SF_Histogram_ = textfileReader.getValue<std::string>("trigger_DoubleElectron_SF_Histogram", "");
    corrections_.trigger_DoubleElectron_SF_Histogram_Format_ = textfileReader.getValue<std::string>("trigger_DoubleElectron_SF_Histogram_Format", "");
    corrections_.trigger_ElectronMuon_SF_File_ = textfileReader.getValue<std::string>("trigger_ElectronMuon_SF_File", "");
    corrections_.trigger_ElectronMuon_SF_Histogram_ = textfileReader.getValue<std::string>("trigger_ElectronMuon_SF_Histogram", "");
    corrections_.trigger_ElectronMuon_SF_Histogram_Format_ = textfileReader.getValue<std::string>("trigger_ElectronMuon_SF_Histogram_Format", "");
    corrections_.trigger_DoubleMuon_SF_File_ = textfileReader.getValue<std::string>("trigger_DoubleMuon_SF_File", "");
    corrections_.trigger_DoubleMuon_SF_Histogram_ = textfileReader.getValue<std::string>("trigger_DoubleMuon_SF_Histogram", "");
    corrections_.trigger_DoubleMuon_SF_Histogram_Format_ = textfileReader.getValue<std::string>("trigger_DoubleMuon_SF_Histogram_Format", "");
    corrections_.electronSFInputFile_ = textfileReader.getValue<std::string>("electronSFInputFile", "");
    corrections_.muonSFInputFile_ = textfileReader.getValue<std::string>("muonSFInputFile", "");
    corrections_.leptonSFReadoutMode_ = textfileReader.getValue<int>("leptonSFReadoutMode", 0);
    corrections_.combinedSF_ = textfileReader.getValue<bool>("combinedSF", true);
    corrections_.SF_electron_ID_File_ = textfileReader.getValue<std::string>("SF_electron_ID_File", "");
    corrections_.SF_electron_ID_Histo_ = textfileReader.getValue<std::string>("SF_electron_ID_Histo", "");
    corrections_.SF_electron_ID_Format_ = textfileReader.getValue<std::string>("SF_electron_ID_Format", "");
    corrections_.SF_electron_Reco_File_ = textfileReader.getValue<std::string>("SF_electron_Reco_File", "");
    corrections_.SF_electron_Reco_Histo_ = textfileReader.getValue<std::string>("SF_electron_Reco_Histo", "");
    corrections_.SF_electron_Reco_Format_ = textfileReader.getValue<std::string>("SF_electron_Reco_Format", "");
    corrections_.SF_electron_Reco_File_LowPt_ = textfileReader.getValue<std::string>("SF_electron_Reco_File_LowPt", "");
    corrections_.SF_electron_Reco_Histo_LowPt_ = textfileReader.getValue<std::string>("SF_electron_Reco_Histo_LowPt", "");
    corrections_.SF_electron_Reco_Format_LowPt_ = textfileReader.getValue<std::string>("SF_electron_Reco_Format_LowPt", "");
    corrections_.extra_electron_ID_DY_extrapolation_unc_ = textfileReader.getValue<float>("extra_electron_ID_DY_extrapolation_unc", 0.);
    corrections_.SF_muon_ID_File_ = textfileReader.getValue<std::string>("SF_muon_ID_File", "");
    corrections_.SF_muon_ID_Histo_ = textfileReader.getValue<std::string>("SF_muon_ID_Histo", "");
    corrections_.SF_muon_ID_Histo_Sys_ = textfileReader.getValue<std::string>("SF_muon_ID_Histo_Sys", "");
    corrections_.SF_muon_ID_Histo_Stat_ = textfileReader.getValue<std::string>("SF_muon_ID_Histo_Stat", "");
    corrections_.SF_muon_ID_Format_ = textfileReader.getValue<std::string>("SF_muon_ID_Format", "");
    corrections_.SF_muon_ISO_File_ = textfileReader.getValue<std::string>("SF_muon_ISO_File", "");
    corrections_.SF_muon_ISO_Histo_ = textfileReader.getValue<std::string>("SF_muon_ISO_Histo", "");
    corrections_.SF_muon_ISO_Format_ = textfileReader.getValue<std::string>("SF_muon_ISO_Format", "");
    corrections_.SF_muon_ISO_Histo_Sys_ = textfileReader.getValue<std::string>("SF_muon_ISO_Histo_Sys", "");
    corrections_.SF_muon_ISO_Histo_Stat_ = textfileReader.getValue<std::string>("SF_muon_ISO_Histo_Stat", "");
    corrections_.extra_muon_ISO_DY_extrapolation_unc_ = textfileReader.getValue<float>("extra_muon_ISO_DY_extrapolation_unc", 0.);
    corrections_.jetVetoMap_File_ = textfileReader.getValue<std::string>("jetVetoMap_File", "");
    corrections_.jetVetoMap_Histogram_ = textfileReader.getValue<std::string>("jetVetoMap_Histogram", "");
    corrections_.jetVetoMap_Format_ = textfileReader.getValue<std::string>("jetVetoMap_Format", "");

    // Lepton JEC ---
    corrections_.applyMuonEnergyCorrections_     = textfileReader.getValue<bool>("applyMuonEnergyCorrections"    , false);
    corrections_.applyElectronEnergyCorrections_ = textfileReader.getValue<bool>("applyElectronEnergyCorrections", false);
    corrections_.applyBJetEnergyRegression_ = textfileReader.getValue<bool>("applyBJetEnergyRegression", false);
    corrections_.propagateCorrectionsToMET_ = textfileReader.getValue<bool>("propagateCorrectionsToMET", false);
    corrections_.jetCorrectionMode_ = textfileReader.getValue<int>("jetCorrectionMode", 0);
    corrections_.jerUncertaintySourceName_ = textfileReader.getValue<std::string>("jerUncertaintySourceName", "");
    corrections_.jerSFBinning_ = textfileReader.getValue<std::string>("jerSFBinning", "Eta");
    corrections_.  jet_JERHybridMethod_PtResolutionFile_ = textfileReader.getValue<std::string>("jet_JERHybridMethod_PtResolutionFile", "");
    corrections_.  jet_JERHybridMethod_jetConeRadius_    = textfileReader.getValue<float>     ("jet_JERHybridMethod_jetConeRadius"   , 0.4);
    corrections_.jesUncertaintySourceFile_ = textfileReader.getValue<std::string>("jesUncertaintySourceFile", "");
    corrections_.jesUncertaintyReducedSourceFile_ = textfileReader.getValue<std::string>("jesUncertaintyReducedSourceFile", "");
    corrections_.jesL1CorrectionDataFile_ = textfileReader.getValue<std::string>("jesL1CorrectionDataFile", "");
    corrections_.jesL2CorrectionDataFile_ = textfileReader.getValue<std::string>("jesL2CorrectionDataFile", "");
    corrections_.jesL3CorrectionDataFile_ = textfileReader.getValue<std::string>("jesL3CorrectionDataFile", "");
    corrections_.jesL2L3CorrectionDataFile_ = textfileReader.getValue<std::string>("jesL2L3CorrectionDataFile", "");
    corrections_.jesL1CorrectionMcFile_ = textfileReader.getValue<std::string>("jesL1CorrectionMcFile", "");
    corrections_.jesL2CorrectionMcFile_ = textfileReader.getValue<std::string>("jesL2CorrectionMcFile", "");
    corrections_.jesL3CorrectionMcFile_ = textfileReader.getValue<std::string>("jesL3CorrectionMcFile", "");
    corrections_.useReducedJESSet_ = textfileReader.getValue<bool>("useReducedJESSet", false);
    corrections_.btagEfficiencyFilename_ = textfileReader.getValue<std::string>("btagEfficiencyFilename", "ttbarsignalplustau_fromDilepton.root");
    corrections_.btagCorrectionMode_ = Btag::convertCorrectionMode(textfileReader.getValue<std::string>("btagCorrectionMode", "noCorrection"));
    corrections_.btagHeavyFlavourFile_ = textfileReader.getValue<std::string>("btagHeavyFlavourFile", "");
    corrections_.btagLightFlavourFile_ = textfileReader.getValue<std::string>("btagLightFlavourFile", "");
    corrections_.btagSFInputFile_ = textfileReader.getValue<std::string>("btagSFInputFile", "");
    corrections_.btagSFCorrectionViaEff_ = textfileReader.getValue<bool>("btagSFCorrectionViaEff", false);
    corrections_.topPtReweighting_ = textfileReader.getValue<bool>("topPtReweighting", false);
    corrections_.applyTopSlope2DGenLevelReweighting_ = textfileReader.getValue<bool>("applyTopSlope2DGenLevelReweighting", false);
    if (corrections_.topPtReweighting_ && corrections_.applyTopSlope2DGenLevelReweighting_) {
      std::cerr << "=====>>>>> IO ERROR: Only one pt reweighting method should be switched on. Please check flags in config file" << std::endl;
      exit(1);
    }
    textfileReader.setRequireValues(true);
    textfileReader.clear();

    // Read object selection
    textfileReader.setStartMarker("[ objectSelection ]");
    textfileReader.setEndMarker("[ end - objectSelection ]");
    textfileReader.readFile(configfilename);
    selections_.leptonEtaCut_ = textfileReader.getValue<float>("leptonEtaCut");
    selections_.leptonPtCut_ = textfileReader.getValue<float>("leptonPtCut");
    selections_.leadingLeptonPtCut_ = textfileReader.getValue<float>("leadingLeptonPtCut");
    selections_.muonIsoCut_ = textfileReader.getValue<float>("muonIsoCut");
    selections_.jetEtaCut_ = textfileReader.getValue<float>("jetEtaCut");
    selections_.jetPtCut_ = textfileReader.getValue<float>("jetPtCut");
    selections_.jetPFIDWorkingPoint_ = textfileReader.getValue<int>("jetPFIDWorkingPoint");
    selections_.jetPileupIdWorkingPoint_ = JetPileupId::convertWorkingPoint(textfileReader.getValue<std::string>("jetPileupIdWorkingPoint"));
    selections_.jetPileupIdMinPt_ = textfileReader.getValue<float>("jetPileupIdMinPt", -1);
    selections_.jetPileupIdMaxPt_ = textfileReader.getValue<float>("jetPileupIdMaxPt", -1);
    selections_.jetVetoMapCutValue_ = textfileReader.getValue<float>("jetVetoMapCutValue", 999.);
    selections_.deltaRLeptonJetCut_ = textfileReader.getValue<float>("deltaRLeptonJetCut");
    selections_.lead2JetPtCut_ = textfileReader.getValue<float>("lead2JetPtCut");
    selections_.btagAlgorithm_ = Btag::convertAlgorithm(textfileReader.getValue<std::string>("btagAlgorithm"));
    selections_.btagWorkingPoint_ = Btag::convertWorkingPoint(textfileReader.getValue<std::string>("btagWorkingPoint"));
    selections_.metAlgorithm_ = Met::convertAlgorithm(textfileReader.getValue<std::string>("metAlgorithm"));
    selections_.metCut_ = textfileReader.getValue<float>("metCut");
    selections_.genJetEtaCut_ = textfileReader.getValue<float>("genJetEtaCut");
    selections_.genJetPtCut_ = textfileReader.getValue<float>("genJetPtCut");
    selections_.genDeltaRLeptonJetCut_ = textfileReader.getValue<float>("genDeltaRLeptonJetCut");
    textfileReader.setRequireValues(false);
    selections_.vetoEventsWithMoreThan2Leptons_ = textfileReader.getValue<int>("vetoEventsWithMoreThan2Leptons", 0);
    textfileReader.setRequireValues(true);
    textfileReader.clear();

    // Read sample composition
    textfileReader.setStartMarker("[ sampleComposition ]");
    textfileReader.setEndMarker("[ end - sampleComposition ]");
    textfileReader.readFile(configfilename);
    sampleComposition_.addSingleLeptonTriggerInfo_ = textfileReader.getValue<bool>("addSingleLeptonTriggerInfo");
    sampleComposition_.pseudodata_ = textfileReader.getValue<int>("pseudodata");
    sampleComposition_.mergeLevel_ = textfileReader.getValue<int>("mergeLevel");
    textfileReader.setRequireValues(false);
    sampleComposition_.pseudoTopMode_ = textfileReader.getValue<int>("pseudoTopMode", 0);
    textfileReader.setRequireValues(true);
    textfileReader.clear();

    // Read plot style
    textfileReader.setStartMarker("[ plotStyle ]");
    textfileReader.setEndMarker("[ end - plotStyle ]");
    textfileReader.readFile(configfilename);
    plotStyle_.cmsLabel_ = textfileReader.getValue<int>("cmsLabel");
    textfileReader.clear();

    std::cout<<"=== Finishing setting up analysis config\n\n";
}



std::string AnalysisConfig::print(const bool screen)const
{
    std::ostringstream ss;
    ss<<general_.print();
    ss<<corrections_.print();
    ss<<selections_.print();
    ss<<sampleComposition_.print();
    ss<<plotStyle_.print();
    const std::string s = ss.str();
    if(screen) std::cout<<s;
    return s;
}
