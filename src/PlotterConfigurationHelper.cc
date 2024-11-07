#include <map>
#include <utility>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <vector>

#include <TMath.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TString.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include "TFractionFitter.h"
#include "TreePlain.h"
#include "utils.h"
#include <TChain.h>
#include "THStack.h"

#include "../../common/include/RootFileReader.h"

#include "PlotterConfigurationHelper.h"
#include "PlotterDiffXSec.h"

/// Put the variation of each systematics one after each other, starting from the UP variation.
/// NOT valid example: MATCH_UP, MASS_DOWN, MASS_UP
/// It is important to keep "BTAG_PT" and "BTAG_ETA" in consecutive order for BCJets and LJets, in example: "BTAG_PT_UP", "BTAG_PT_DOWN", "BTAG_ETA_UP", "BTAG_ETA_DOWN"


const std::map<TString, TString> VectorOfValidSystematicsCommon
{
{"JER", "ExpSys"},
/*
//{"JEREta0", "ExpSys"},
//{"JEREta1", "ExpSys"},
*/
{"UNCLUSTERED", "ExpSys"},
{"PU", "ExpSys"},
{"TRIG", "ExpSys"},
/*
//  {"TRIG_ETA", "ExpSys"},
//  {"LEPT", "ExpSys"},
//  {"ELE", "ExpSys"},
//  {"MUON", "ExpSys"},
*/
{"ELE_ID", "ExpSys"},
{"ELE_RECO", "ExpSys"},
//  {"ELE_SCALESMEARING", "ExpSys"},
{"ELE_SCALE_SYST", "ExpSys"},
{"ELE_SCALE_GAIN", "ExpSys"},
{"ELE_SCALE_STAT", "ExpSys"},
{"MUON_ID", "ExpSys"},
{"MUON_ISO", "ExpSys"},
{"MUON_SCALE", "ExpSys"},
{"DY", "ExpSys"},
{"BG", "ExpSys"},
//  {"TW", "ExpSys"},
//  {"LUMI", "ExpSys"},
//  {"KIN", "ExpSys"},
{"L1PREFIRING", "ExpSys"},
{"BTAG_CORR", "ExpSys"},
{"BTAG_LJET_CORR", "ExpSys"},
{"BTAG_UNCORR", "ExpSys"},
{"BTAG_LJET_UNCORR", "ExpSys"},
/*
// {"BTAG_PT", "ExpSys"},
// {"BTAG_ETA", "ExpSys"},
// {"BTAG_LJET_PT", "ExpSys"},
// {"BTAG_LJET_ETA", "ExpSys"},
// {"BTAG_BEFF", "ExpSys"},
// {"BTAG_CEFF", "ExpSys"},
// {"BTAG_LEFF", "ExpSys"},
*/
{"JES", "JESSys"},
{"JESAbsoluteStat", "JESSysSplit"}, // this JES source should be always the first (for unfolding with SVD, but keep like this for now)
{"JESAbsoluteScale", "JESSysSplit"},
//{"JESAbsoluteFlavMap", "JESSysSplit"},
{"JESAbsoluteMPFBias","JESSysSplit"},
//{"JESFlavorQCD", "JESSysSplit"},
{"JESFlavorPureBottom", "JESSysSplit"},
{"JESFlavorPureCharm", "JESSysSplit"},
{"JESFlavorPureGluon", "JESSysSplit"},
{"JESFlavorPureQuark", "JESSysSplit"},
{"JESFlavorPhotonJet", "JESSysSplit"},
{"JESFlavorRealistic", "JESSysSplit"},
{"JESFlavorZJet", "JESSysSplit"},
{"JESFragmentation", "JESSysSplit"},
{"JESPileUpDataMC", "JESSysSplit"},
{"JESPileUpPtRef", "JESSysSplit"},
{"JESPileUpPtEC1", "JESSysSplit"},
{"JESPileUpPtEC2", "JESSysSplit"},
{"JESPileUpPtHF", "JESSysSplit"},
{"JESPileUpPtBB", "JESSysSplit"},
{"JESRelativeBal", "JESSysSplit"},
{"JESRelativeJEREC1", "JESSysSplit"},
{"JESRelativeJEREC2", "JESSysSplit"},
{"JESRelativeJERHF", "JESSysSplit"},
{"JESRelativePtBB", "JESSysSplit"},
{"JESRelativePtEC1", "JESSysSplit"},
{"JESRelativePtEC2", "JESSysSplit"},
{"JESRelativePtHF", "JESSysSplit"},
{"JESRelativeFSR", "JESSysSplit"},
{"JESRelativeStatFSR", "JESSysSplit"},
{"JESRelativeStatEC", "JESSysSplit"},
{"JESRelativeStatHF", "JESSysSplit"},
{"JESSinglePionECAL", "JESSysSplit"},
{"JESSinglePionHCAL", "JESSysSplit"},
{"JESTimePtEta", "JESSysSplit"},

{"UETUNE", "ModISys"},
{"MATCH", "ModISys"},
{"MASS", "ModISys"},

//{"TOT_COLORREC", "ModISys"},

{"MESCALE", "ModWSys"},  {"MEFACSCALE", "ModWSys"},  {"MERENSCALE", "ModWSys"},
//{"TOT_SCALE", "ModWSys"},
{"BSEMILEP", "ModWSys"},
{"PDF_ALPHAS", "ModWSys"},
//{"BFRAG", "ModWSys"},
//{"TOT_BFRAG", "ModWSys"},

{"PDF_1", "ModWSys"},  {"PDF_2", "ModWSys"},
{"PDF_3", "ModWSys"},  {"PDF_4", "ModWSys"},
{"PDF_5", "ModWSys"},  {"PDF_6", "ModWSys"},
{"PDF_7", "ModWSys"},  {"PDF_8", "ModWSys"},
{"PDF_9", "ModWSys"},  {"PDF_10", "ModWSys"},
{"PDF_11", "ModWSys"},  {"PDF_12", "ModWSys"},
{"PDF_13", "ModWSys"},  {"PDF_14", "ModWSys"},
{"PDF_15", "ModWSys"},  {"PDF_16", "ModWSys"},
{"PDF_17", "ModWSys"},  {"PDF_18", "ModWSys"},
{"PDF_19", "ModWSys"},  {"PDF_20", "ModWSys"},
{"PDF_21", "ModWSys"},  {"PDF_22", "ModWSys"},
{"PDF_23", "ModWSys"},  {"PDF_24", "ModWSys"},
{"PDF_25", "ModWSys"},  {"PDF_26", "ModWSys"},
{"PDF_27", "ModWSys"},  {"PDF_28", "ModWSys"},
{"PDF_29", "ModWSys"},  {"PDF_30", "ModWSys"},
{"PDF_31", "ModWSys"},  {"PDF_32", "ModWSys"},
{"PDF_33", "ModWSys"},  {"PDF_34", "ModWSys"},
{"PDF_35", "ModWSys"},  {"PDF_36", "ModWSys"},
{"PDF_37", "ModWSys"},  {"PDF_38", "ModWSys"},
{"PDF_39", "ModWSys"},  {"PDF_40", "ModWSys"},
{"PDF_41", "ModWSys"},  {"PDF_42", "ModWSys"},
{"PDF_43", "ModWSys"},  {"PDF_44", "ModWSys"},
{"PDF_45", "ModWSys"},  {"PDF_46", "ModWSys"},
{"PDF_47", "ModWSys"},  {"PDF_48", "ModWSys"},
{"PDF_49", "ModWSys"},  {"PDF_50", "ModWSys"},

//{"PDF", "ModWSys"},
//{"TOT_PDF", "ModWSys"},  

//  {"PSSCALE_WEIGHT_2", "ModWSys"}, //ISR Reduced
//  {"PSSCALE_WEIGHT_3", "ModWSys"}, //FSR Reduced
//  {"PSSCALE_WEIGHT_4", "ModWSys"}, //ISR
//  {"PSSCALE_WEIGHT_5", "ModWSys"}, //FSR
{"PSSCALE_WEIGHT_6", "ModWSys"}, //ISR Conservative
{"PSSCALE_WEIGHT_7", "ModWSys"}, //FSR Conservative
//  {"PSSCALE_WEIGHT_8", "ModWSys"}, //fsr_G2GG_muR
//  {"PSSCALE_WEIGHT_9", "ModWSys"}, //fsr_G2QQ_muR
//  {"PSSCALE_WEIGHT_10", "ModWSys"}, //fsr_Q2QG_muR
//  {"PSSCALE_WEIGHT_11", "ModWSys"}, //fsr_X2XG_muR
//  {"PSSCALE_WEIGHT_12", "ModWSys"}, //fsr_G2GG_cNS
//  {"PSSCALE_WEIGHT_13", "ModWSys"}, //fsr_G2QQ_cNS
//  {"PSSCALE_WEIGHT_14", "ModWSys"}, //fsr_Q2QG_cNS
//  {"PSSCALE_WEIGHT_15", "ModWSys"}, //fsr_X2XG_cNS
//  {"PSSCALE_WEIGHT_16", "ModWSys"}, //isr_G2GG_muR
//  {"PSSCALE_WEIGHT_17", "ModWSys"}, //isr_G2QQ_muR
//  {"PSSCALE_WEIGHT_18", "ModWSys"}, //isr_Q2QG_muR
//  {"PSSCALE_WEIGHT_19", "ModWSys"}, //isr_X2XG_muR
//  {"PSSCALE_WEIGHT_20", "ModWSys"}, //isr_G2GG_cNS
//  {"PSSCALE_WEIGHT_21", "ModWSys"}, //isr_G2QQ_cNS
//  {"PSSCALE_WEIGHT_22", "ModWSys"}, //isr_Q2QG_cNS
//  {"PSSCALE_WEIGHT_23", "ModWSys"}, //isr_X2XG_cNS


};

const std::map<TString, TString> VectorOfValidSystematicsCommon_noUpDown
{

//{"BFRAG_PETERSON", "ModWSys"},
{"PDF_0_CENTRAL", "ModWSys"},

//  {"BFRAG_CENTRAL", "ModWSys"},

//  {"MCATNLO", "AltTheorySys"},
//  {"POWHEG", "AltTheorySys"},
{"MADGRAPHMLM", "AltTheorySys"}, {"AMCATNLOFXFX", "AltTheorySys"},  {"POWHEGV2HERWIG", "AltTheorySys"},

{"ERDON", "ModISys"},  {"ERDONRETUNE", "ModISys"},  {"GLUONMOVETUNE", "ModISys"},

{"TOP_PT", "ModWSys"}
};

// Year exclusive systematics.
const std::map<TString, TString> VectorOfValidSystematicsOnly2016
{
//  {"PSISRSCALE", "ModISys"},  // M2T4 setup (for now)
//  {"PSFSRSCALE", "ModISys"},  // M2T4 setup (for now)
//  {"PSISRSCALE_2", "ModWSys"},// CP5 setup
//  {"PSFSRSCALE_2", "ModWSys"} // CP5 setup
};

const std::map<TString, TString> VectorOfValidSystematicsOnly2016_noUpDown
{
    
};

const std::map<TString, TString> VectorOfValidSystematicsOnly2016preVFP
{

};

const std::map<TString, TString> VectorOfValidSystematicsOnly2016preVFP_noUpDown
{
    
};

const std::map<TString, TString> VectorOfValidSystematicsOnly2016postVFP
{

};

const std::map<TString, TString> VectorOfValidSystematicsOnly2016postVFP_noUpDown
{
    
};

const std::map<TString, TString> VectorOfValidSystematicsOnly2017
{

};

const std::map<TString, TString> VectorOfValidSystematicsOnly2017_noUpDown
{

};

const std::map<TString, TString> VectorOfValidSystematicsOnly2018
{

};

const std::map<TString, TString> VectorOfValidSystematicsOnly2018_noUpDown
{

};

const std::map<TString, TString> VectorOfValidSystematicsOnlyFullRun2UL
{

};

// Exclude systematics from VectorOfValidSystematicsCommon for a specific year.
const std::map<TString, TString> VectorOfExcludedSystematicsOnly2016
{
    // {"MASS", "ModISys"},
    // {"ERDONRETUNE", "ModISys"},
    // {"GLUONMOVETUNE", "ModISys"}
};
const std::map<TString, TString> VectorOfExcludedSystematicsOnly2016preVFP
{

};

const std::map<TString, TString> VectorOfExcludedSystematicsOnly2016postVFP
{

};

const std::map<TString, TString> VectorOfExcludedSystematicsOnly2017
{

};

const std::map<TString, TString> VectorOfExcludedSystematicsOnly2018
{
    
};

const std::map<TString, TString> VectorOfExcludedSystematicsOnlyFullRun2UL
{

};

// List of systematics to be included via option from command line
const std::map<TString, TString> VectorOfValidSystematicsAddons
{
  {"ELE", "ExpSys"},
  {"MUON", "ExpSys"}
};

const std::map<TString, TString> VectorOfValidSystematicsAddons_noUpDown
{
};

const std::map<TString, TString> VectorOfExcludedSysts_ttmd_only{
    // {"UETUNE", "ModISys"},
    // {"PU", "ExpSys"}
};

const std::map<TString, std::vector<TString>> EnvelopesMap{ //Currently only used in ttmd code

    {"TOT_COLORREC",
        {"ERDON", "ERDONRETUNE", "GLUONMOVETUNE"}
    },
      //    {"TOT_BFRAG",
      //        {"BFRAG_PETERSON", "BFRAG_UP", "BFRAG_DOWN"}
      //    },
    {"TOT_SCALE",
        {"MESCALE_UP", "MESCALE_DOWN", "MEFACSCALE_UP", "MEFACSCALE_DOWN", "MERENSCALE_UP", "MERENSCALE_DOWN"}
    },
    {"TOT_PDF",
        { "PDF_0_CENTRAL",
	  "PDF_1_UP", "PDF_1_DOWN",
	  "PDF_2_UP", "PDF_2_DOWN",
	  "PDF_3_UP", "PDF_3_DOWN",
	  "PDF_4_UP", "PDF_4_DOWN",
	  "PDF_5_UP", "PDF_5_DOWN",
	  "PDF_6_UP", "PDF_6_DOWN",
	  "PDF_7_UP", "PDF_7_DOWN",
	  "PDF_8_UP", "PDF_8_DOWN",
	  "PDF_9_UP", "PDF_9_DOWN",
	  "PDF_10_UP", "PDF_10_DOWN",
	  "PDF_11_UP", "PDF_11_DOWN",
	  "PDF_12_UP", "PDF_12_DOWN",
	  "PDF_13_UP", "PDF_13_DOWN",
	  "PDF_14_UP", "PDF_14_DOWN",
	  "PDF_15_UP", "PDF_15_DOWN",
	  "PDF_16_UP", "PDF_16_DOWN",
	  "PDF_17_UP", "PDF_17_DOWN",
	  "PDF_18_UP", "PDF_18_DOWN",
	  "PDF_19_UP", "PDF_19_DOWN",
	  "PDF_20_UP", "PDF_20_DOWN",
	  "PDF_21_UP", "PDF_21_DOWN",
	  "PDF_22_UP", "PDF_22_DOWN",
	  "PDF_23_UP", "PDF_23_DOWN",
	  "PDF_24_UP", "PDF_24_DOWN",
	  "PDF_25_UP", "PDF_25_DOWN",
	  "PDF_26_UP", "PDF_26_DOWN",
	  "PDF_27_UP", "PDF_27_DOWN",
	  "PDF_28_UP", "PDF_28_DOWN",
	  "PDF_29_UP", "PDF_29_DOWN",
	  "PDF_30_UP", "PDF_30_DOWN",
	  "PDF_31_UP", "PDF_31_DOWN",
	  "PDF_32_UP", "PDF_32_DOWN",
	  "PDF_33_UP", "PDF_33_DOWN",
	  "PDF_34_UP", "PDF_34_DOWN",
	  "PDF_35_UP", "PDF_35_DOWN",
	  "PDF_36_UP", "PDF_36_DOWN",
	  "PDF_37_UP", "PDF_37_DOWN",
	  "PDF_38_UP", "PDF_38_DOWN",
	  "PDF_39_UP", "PDF_39_DOWN",
	  "PDF_40_UP", "PDF_40_DOWN",
	  "PDF_41_UP", "PDF_41_DOWN",
	  "PDF_42_UP", "PDF_42_DOWN",
	  "PDF_43_UP", "PDF_43_DOWN",
	  "PDF_44_UP", "PDF_44_DOWN",
	  "PDF_45_UP", "PDF_45_DOWN",
	  "PDF_46_UP", "PDF_46_DOWN",
	  "PDF_47_UP", "PDF_47_DOWN",
	  "PDF_48_UP", "PDF_48_DOWN",
	  "PDF_49_UP", "PDF_49_DOWN",
	  "PDF_50_UP", "PDF_50_DOWN"}
    }

};

//references
//prescriptions exit for JES, JER, Lumi. No recommendaion from BTV POG yet so following arguments regarding correlations from:
//https://indico.cern.ch/event/853327/contributions/3623147/attachments/1935532/3207317/ttx_October2019_Ranken.pdf (slide 7,8,11)
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETRun2Corrections#Uncertainty_correlations_among_y (UNCLUSTERED, JER, JES)
//https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#Run2_JER_uncertainty_correlation (JER)
//https://docs.google.com/spreadsheets/d/1JZfk78_9SD225bcUuTWVo4i02vwI5FfeVKH-dwzUdhM/edit#gid=1345121349 (JES sources - split)
//https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#LumiComb (LUMI)
const std::map<TString, std::vector<double>> YearToYearCorrMap{
    //correlation factors follow the notation: {1617, 1618, 1718}
    {"JER", {0.,0.,0.}},
      // {"JEREta0", {0.,0.,0.}},
      // {"JEREta1", {0.,0.,0.}},
    {"UNCLUSTERED", {0.,0.,0.}},
    {"LUMI", {0.3,0.3,0.3}},
    {"BTAG", {0.5,0.5,0.5}},
    {"BTAG_LJET", {0.5,0.5,0.5}},
      //    {"BTAG_PT", {0.5,0.5,0.5}},
      //    {"BTAG_ETA", {0.5,0.5,0.5}},
      //    {"BTAG_LJET_PT", {0.5,0.5,0.5}},
      //    {"BTAG_LJET_ETA", {0.5,0.5,0.5}},
      //    {"BTAG_BEFF", {0.5,0.5,0.5}},
      //    {"BTAG_CEFF", {0.5,0.5,0.5}},
      //    {"BTAG_LEFF", {0.5,0.5,0.5}},
    {"JESAbsoluteStat", {0.,0.,0.}},
    {"JESTimePtEta", {0.,0.,0.}},
    {"JESRelativeBal", {0.5,0.5,0.5}},
    {"JESRelativeJEREC1", {0.,0.,0.}},
    {"JESRelativeJEREC2", {0.,0.,0.}},
    {"JESRelativeJERHF", {0.5,0.5,0.5}},
    {"JESRelativePtBB", {0.5,0.5,0.5}},
    {"JESRelativePtEC1", {0.,0.,0.}},
    {"JESRelativePtEC2", {0.,0.,0.}},
    {"JESRelativePtHF", {0.5,0.5,0.5}},
    {"JESRelativeFSR", {0.5,0.5,0.5}},
    {"JESRelativeStatFSR", {0.,0.,0.}},
    {"JESRelativeStatEC", {0.,0.,0.}},
    {"JESRelativeStatHF", {0.,0.,0.}},
    {"JESPileUpDataMC", {0.5,0.5,0.5}},
    {"JESPileUpPtRef", {0.5,0.5,0.5}},
    {"JESPileUpPtEC1", {0.5,0.5,0.5}},
    {"JESPileUpPtEC2", {0.5,0.5,0.5}},
    {"JESPileUpPtHF", {0.5,0.5,0.5}},
    {"JESPileUpPtBB", {0.5,0.5,0.5}}
};

std::pair<double,double> PlotterConfigurationHelper::get_lumi_and_err(TString in_year){

    std::pair<double,double> lumi_info;
    double y_lumi, y_lumiError;

    if (in_year == "2016preVFP") {
        y_lumi = 19500.;
         y_lumiError = 0.012;
    } else if (in_year == "2016postVFP") {
        y_lumi = 16810.;
         y_lumiError = 0.012;
    } else if (in_year == "2016") {
        y_lumi = 36310.;
         y_lumiError = 0.012;
    } else if (in_year == "2017") {
        y_lumi = 41480.;
         y_lumiError = 0.023;
    } else if (in_year == "2018") {
        y_lumi = 59830.;
        y_lumiError = 0.025;
    } else if (in_year == "fullRun2UL") {
        y_lumi = 137620.;
        y_lumiError = 0.016;
    } else {
        std::cerr << "ERROR in PlotterConfigurationHelper::PlotterConfigurationHelper() : year is not supported " << std::endl;
        exit(1);
    }

    lumi_info.first = y_lumi;
    lumi_info.second = y_lumiError;

    return lumi_info;

}

void PlotterConfigurationHelper::get_year_from_file(const TString& filename){

    if(filename == "")
      return;

    lumi_year.Clear();

    if (filename.Contains("2016preVFP")) lumi_year = "2016preVFP";
    else if (filename.Contains("2016postVFP")) lumi_year = "2016postVFP";
    else if (filename.Contains("2016")) lumi_year = "2016";
    else if (filename.Contains("2017")) lumi_year = "2017";
    else if (filename.Contains("2018")) lumi_year = "2018";
    else if (filename.Contains("fullRun2UL")) lumi_year = "fullRun2UL";
    else {
      std::cerr << "ERROR in PlotterConfigurationHelper::get_year_from_file() : no year string in" << filename  << std::endl;
      exit(1);
    }

}

double PlotterConfigurationHelper::get_year_corr_rescale_from_file(const TString& complete_Var_name, const TString& filepath){

    double yearToYearCorrRescale = 1.;

    if (filepath.Contains("run"))
      return yearToYearCorrRescale = 1.;

    TString sys = VarToSys(complete_Var_name);
    bool isBkg = false;
    bool isTTBkg = false;

    if (filepath.Contains("dyee") || filepath.Contains("dymumu") || filepath.Contains("dytautau") || filepath.Contains("wtolnu")
	                          || filepath.Contains("wwtoall") || filepath.Contains("wztoall") || filepath.Contains("zztoall")
	                          || filepath.Contains("ttbarW") || filepath.Contains("ttbarZ") || filepath.Contains("ttgjets")
	                          || filepath.Contains("single"))
      isBkg = true;

    if (filepath.Contains("ttbarbgviatau") || filepath.Contains("ttbarbg_") || filepath.Contains("ttbarbg."))
      isTTBkg = true;

    std::vector<TString> list_years_Corr;

    if (complete_Var_name.Contains("16")) {
      list_years_Corr.push_back("2016");
    }
    if (complete_Var_name.Contains("17")) {
      list_years_Corr.push_back("2017");
    }
    if (complete_Var_name.Contains("18")) {
      list_years_Corr.push_back("2018");
    }

    std::vector<TString> other_years;
    if (year_ == "fullRun2UL"  && yearCorr_ != "" && IsYearToYearCorr(year_, yearCorr_, complete_Var_name) && list_years_Corr.size() != 0) {
      other_years = this->years_list;
      for(auto iYearCorr = list_years_Corr.begin(); iYearCorr != list_years_Corr.end(); ++iYearCorr) {
	int idx;
	auto it = std::find(other_years.begin(), other_years.end(), *iYearCorr);
	if (it != other_years.end())
	  {
	    idx = std::distance(other_years.begin(), it);
	    other_years.erase(other_years.begin()+idx);
	  }
      }
    }

    get_year_from_file(filepath);

    if ((filepath.Contains("Nominal") && complete_Var_name == "Nominal") || yearCorr_ == "") {
      //case 1
      //no correlation treatment should be performed for these two cases: running Nominal or not running with --doCorrYear. In the latter case all systematics are treated as 100%        correlated between years by default
      yearToYearCorrRescale = 1.;
    }
    else if (yearCorr_ != "" && (std::find(other_years.begin(), other_years.end(), lumi_year) != other_years.end())) {
      //case 2
      //sample dir should be nominal for years not contained in name of systematic. I.e. for BTAG_1617 one needs to take nominal samples for 2018
      if (filepath.Contains("Nominal")) {
	yearToYearCorrRescale = 1.;
      }
      else {
        std::cerr << "ERROR in PlotterConfigurationHelper::get_year_corr_rescale_from_file() : year variation for year = " << lumi_year << " shouldn't be performed but the path doesn't refer to Nominal" << std::endl;
	std::cerr << "file path = " << filepath << std::endl;
	exit(1);
      }
    }
    else if (yearCorr_ != "" && !filepath.Contains("Nominal") && ((IsAltTheorSys(sys) && (isTTBkg || isBkg)) || ((IsModISys(sys) || IsModWSys(sys)) && isBkg))) {
      //case 3
      //if the sample dir is not nominal but the systematic variation is alternative theory and the sample is not signal the don't recale
      //if the sample dir is not nominal but the systematic variation is not experimental and the sample is bg then  don't rescale
      yearToYearCorrRescale = 1.;
    }
    else {
      yearToYearCorrRescale = GetYearToYearCorrRescale(year_, yearCorr_, complete_Var_name);
    }

    return yearToYearCorrRescale;
}

TString PlotterConfigurationHelper::getCorrespondingNominalDataset(const TString& datasetPath, const TString& systematic, const TString& channel){

    TString histoListNameNominal = baseFileList + "HistoFileList_Nominal_"+channel+".txt";
    std::ifstream FileList(histoListNameNominal);
    if (FileList.fail()) {
        std::cerr << "Error reading " << histoListNameNominal << std::endl;
        exit(1);
    }

    TString filePatternTTbar = "";
    if (datasetPath.Contains("ttbarsignalplustau_fromDilepton"))
      filePatternTTbar = "ttbarsignalplustau_fromDilepton";
    else if (datasetPath.Contains("ttbarsignalviatau_fromDilepton"))
      filePatternTTbar = "ttbarsignalviatau_fromDilepton";
    else if (datasetPath.Contains("ttbarsignalplustau"))
      filePatternTTbar = "ttbarsignalplustau";
    else if (datasetPath.Contains("ttbarsignalviatau"))
      filePatternTTbar = "ttbarsignalviatau";
    else if (datasetPath.Contains("ttbarbg_fromDilepton"))
      filePatternTTbar = "ttbarbg_fromDilepton";
    else if (datasetPath.Contains("ttbarbg_fromHadronic"))
      filePatternTTbar = "ttbarbg_fromHadronic";
    else if (datasetPath.Contains("ttbarbg_fromLjets"))
      filePatternTTbar = "ttbarbg_fromLjets";
    else if (datasetPath.Contains("ttbarbgviatau_fromDilepton"))
      filePatternTTbar = "ttbarbgviatau_fromDilepton";
    else if (datasetPath.Contains("ttbarbgviatau"))
      filePatternTTbar = "ttbarbgviatau";
    else if (datasetPath.Contains("ttbarbg"))
      filePatternTTbar = "ttbarbg";

    get_year_from_file(datasetPath);
    TString originalDatasetName = ExtractFileNameFromPath(datasetPath,false);
    TString correspondingNominalPath;
    TString datasetPathNominal;
    while (!FileList.eof()) {
        FileList >> datasetPathNominal;
	if (datasetPathNominal.Contains(originalDatasetName) && datasetPathNominal.Contains(lumi_year)) {
	  correspondingNominalPath = datasetPathNominal;
	  break;
	}
	else if (IsModISys(VarToSys(systematic)) || IsAltTheorSys(VarToSys(systematic))) {
	  if ((filePatternTTbar != "" && datasetPathNominal.Contains(filePatternTTbar)) && datasetPathNominal.Contains(lumi_year)) {
	    correspondingNominalPath = datasetPathNominal;
	    break;
	  }
	}
    }
    FileList.close();

    if (correspondingNominalPath == "") {
      std::cerr << "Error: no corresponding dataset in nominal filelist for the following: " << std::endl;
      std::cerr << "nominal file list = " << histoListNameNominal << std::endl;
      std::cerr << "systematic = " << systematic << std::endl;
      std::cerr << "dataset = " << originalDatasetName << std::endl;
      exit(1);
    }

    return correspondingNominalPath;

}

void PlotterConfigurationHelper::doYearToYearCorrRescaleSanityCheck(const TString& Var, const TString dirInput){

    TString temp_var = Var;
    RemoveYearCorr(year_, yearCorr_, temp_var);

    std::vector<TString> vars_to_fill;

    std::vector<TString> yearCorrList = GetYearCorrList(year_, yearCorr_, temp_var, IsYearToYearCorr(year_, yearCorr_, Var));

    //if (yearCorrList.size() == 0)
    //  return;

    TString var_16;
    TString var_17;
    TString var_18;
    TString var_161718;

    if(yearCorrList.size() != 0) {
      var_16 = temp_var;
      AppendYearCorr("fullRun2UL", "16", var_16, true);
      var_17 = temp_var;
      AppendYearCorr("fullRun2UL", "17", var_17, true);
      var_18 = temp_var;
      AppendYearCorr("fullRun2UL", "18", var_18, true);
      var_161718 = temp_var;
      AppendYearCorr("fullRun2UL", "161718", var_161718, true);

      if (std::find(begin(yearCorrList), end(yearCorrList), var_16) == end(yearCorrList)) return;
      if (std::find(begin(yearCorrList), end(yearCorrList), var_17) == end(yearCorrList)) return;
      if (std::find(begin(yearCorrList), end(yearCorrList), var_18) == end(yearCorrList)) return;
      if (std::find(begin(yearCorrList), end(yearCorrList), var_161718) == end(yearCorrList)) return;
      
      vars_to_fill.insert(vars_to_fill.end(), yearCorrList.begin(), yearCorrList.end());
      vars_to_fill.insert(vars_to_fill.begin(), "Nominal");}
    
    else {
      if (Var == "Nominal") return;
      vars_to_fill.push_back("Nominal");
      vars_to_fill.push_back(Var);
    }

    for (unsigned int i = 1; i < 5; i++) {

      std::map<TString, TH1D*> test_var;
      TString sampleType;

      if (i == 1)
	sampleType = "Rec";
      else if (i == 2)
	sampleType = "TTBkg";
      else if (i == 3)
	sampleType = "Bkg";
      else if (i == 4)
	sampleType = "Gen";

      for(const TString& var : vars_to_fill){
	const TString& dirTree = GetVarDirPlot(var);
	std::vector<TH1D*> readCP = PlotterDiffXSec::ReadCP(dirInput + "/" + dirTree + "/cp.txt");
	//printf("integral = %f \n", h->Integral());
	test_var[var] = (TH1D*) readCP[i]->Clone();

      }

      if(yearCorrList.size() != 0) {
	TH1D* h_Nom = NULL;
	TH1D* h_var_16 = NULL;
	TH1D* h_var_17 = NULL;
	TH1D* h_var_18 = NULL;
	TH1D* h_var_161718 = NULL;
	for(const TString& var : vars_to_fill){

	  TString label = GetLabelYearCorr(var);

	  if (var == "Nominal") {
	    h_Nom = (TH1D*) test_var.at(var)->Clone();
	  }
	  else if (label == "16") {
	    h_var_16 = (TH1D*) test_var.at(var)->Clone();
	  }
	  else if (label == "17") {
	    h_var_17 = (TH1D*) test_var.at(var)->Clone();
	  }
	  else if (label == "18") {
	    h_var_18 = (TH1D*) test_var.at(var)->Clone();
	  }
	  else if (label == "161718") {
	    h_var_161718 = (TH1D*) test_var.at(var)->Clone();
	  }
	}

	if (!h_Nom || !h_var_16 || !h_var_17 || !h_var_18 || !h_var_161718) {
	  std::cout << "PlotterConfigurationHelper::doYearToYearCorrRescaleSanityCheck: one or more variation histograms aren't set properly \n";
	  std::exit(1);
	}

	TH1D* hDiff_var_16 = (TH1D*) h_var_16->Clone();
	hDiff_var_16->Add(h_Nom, -1.0);

	TH1D* hDiff_var_17 = (TH1D*) h_var_17->Clone();
	hDiff_var_17->Add(h_Nom, -1.0);

	TH1D* hDiff_var_18 = (TH1D*) h_var_18->Clone();
	hDiff_var_18->Add(h_Nom, -1.0);

	TH1D* hDiff_var_161718 = (TH1D*) h_var_161718->Clone();
	hDiff_var_161718->Add(h_Nom, -1.0);

	if (hDiff_var_16->GetEntries() == 0 || hDiff_var_17->GetEntries() == 0 || hDiff_var_18->GetEntries() == 0 || hDiff_var_161718->GetEntries() == 0) {
	  printf("%s variation is not performed for %s \n", Var.Data(), sampleType.Data());
	  continue;
	}

	printf("h_Nom = %f \n", h_Nom->Integral());
	printf("h_var_16->Integral() = %f \n", h_var_16->Integral());
	printf("h_var_17->Integral() = %f \n", h_var_17->Integral());
	printf("h_var_18->Integral() = %f \n", h_var_18->Integral());
	printf("h_var_161718->Integral() = %f \n", h_var_161718->Integral());

	double rescale_var_16 = GetYearToYearCorrRescale("fullRun2UL", "all", var_16);
	double rescale_var_17 = GetYearToYearCorrRescale("fullRun2UL", "all", var_17 );
	double rescale_var_18 = GetYearToYearCorrRescale("fullRun2UL", "all", var_18);
	double rescale_var_161718 = GetYearToYearCorrRescale("fullRun2UL", "all", var_161718);

	printf("before \n");
	printf("hDiff_var_16->Integral() = %f \n", hDiff_var_16->Integral());
	printf("hDiff_var_17->Integral() = %f \n", hDiff_var_17->Integral());
	printf("hDiff_var_18->Integral() = %f \n", hDiff_var_18->Integral());
	printf("hDiff_var_161718->Integral() = %f \n", hDiff_var_161718->Integral());

	hDiff_var_16->Scale(1./rescale_var_16);
	hDiff_var_17->Scale(1./rescale_var_17);
	hDiff_var_18->Scale(1./rescale_var_18);
	TH1D* hDiff_var_161718_clone = (TH1D*) hDiff_var_161718->Clone();
	hDiff_var_161718_clone->Scale(1./rescale_var_161718);

	printf("after \n");
	printf("hDiff_var_16->Integral() = %f \n", hDiff_var_16->Integral());
	printf("hDiff_var_17->Integral() = %f \n", hDiff_var_17->Integral());
	printf("hDiff_var_18->Integral() = %f \n", hDiff_var_18->Integral());
	printf("hDiff_var_161718->Integral() = %f \n", hDiff_var_161718_clone->Integral());

	TH1D* h_add_16_17_18 = (TH1D*) hDiff_var_16->Clone();
	h_add_16_17_18->Add(hDiff_var_16, 1.0);
	h_add_16_17_18->Add(hDiff_var_17, 1.0);
	h_add_16_17_18->Add(hDiff_var_18, 1.0);

	printf("//********************* %s, %s ******************************** \n", Var.Data(), sampleType.Data());
	double no_rescale_161718 = h_add_16_17_18->Integral();
	double with_rescale_161718 = hDiff_var_161718->Integral();
	printf("variation 161718 with no rescale = %f for %s \n", no_rescale_161718, sampleType.Data());
	printf("variation 161718 with rescale = %f for %s \n", with_rescale_161718, sampleType.Data());
	printf("161718 rescale factor extracted from histogram = %f for %s \n", with_rescale_161718/no_rescale_161718, sampleType.Data());
	printf("161718 rescale factor applied at TreePlain level = %f for %s \n", rescale_var_161718, sampleType.Data());
	printf("//************************************************************ \n");
      }
      else {
	TH1D* h_Nom = NULL;
	TH1D* h_var = NULL;
	for(const TString& var : vars_to_fill){
	  
	  if(var == "Nominal")
	    h_Nom = (TH1D*) test_var.at(var)->Clone();
	  else
	    h_var = (TH1D*) test_var.at(var)->Clone();
	}
	
	TH1D* hDiff_var = (TH1D*) h_var->Clone();
	hDiff_var->Add(h_Nom, -1.0);
	
	if (hDiff_var->GetEntries() == 0) {
	  printf("%s variation is not performed for %s \n", Var.Data(), sampleType.Data());
	  continue;
	  
	}
      }
    }
}

PlotterConfigurationHelper::PlotterConfigurationHelper(RootFileReader* rootFileReader, const TString& year, const TString& yearCorr, bool isDYScale, bool setDiffXSec, bool flagProcSyst):
energy(13.0), // centre-of-mass energy in TeV
lumi(1.), // data luminosity in pb-1
lumiError(100.0), // in percentages
topxsec(830.91), // in pb, taken from taken from private calculation using Top++ 2.0 with NNPDF31_nnlo_hessian_pdfas (LHAPDF  306000) / 16.04.2021
bRatio{0.01147,0.01130,0.02277,0.04554},//ttbar->dilepton branching ratios for [ee, mumu, emu, combined] not including decays via tau
bRatioError(0.015), // in percentages
massUncReductionFactor(1.0), // relative factor reducing absolute uncertainty due to top quark mass to match +-1GeV variation (assuming usage of 169.5 (down) and 175.5 (up) samples)
fsrScaleUncReductionFactor(TMath::Sqrt(2.0)), // relative factor reducing absolute uncertainty due to PS FSR scale
revokeAntiQuantityCombination(false), //Prevents construction of one (averaged, if applicable) histogram out of two objects : quantity and anti-quantity, if allowed for given quantities (e.g. top quark & top anti-quark)
drawTotalMCErrorForCP(true), //Draw error band on control plot as obtained from total variations in all MCs (instead of standalone variations in ttbar)
drawPlotRatio(true),
cmsPrelimLabelMode(1), // 0 - CMS, 1 - CMS Preliminary, 2 - CMS Private Work
usePoissonSmearedPseudoData(false),//instead of real data, use pseudo-data defined as poisson-smeared sum of all MCs (requires relevant lists of input files; see README.txt for more details)
estimateSysUncFromPseudoData(!usePoissonSmearedPseudoData ? false : false),// Setting to "true" triggers estimation of syst. unc. from 'varied pseudo-data', aka estimation of expected sys. uncertainty (i.e. in case for each variation: using always nominal MC simulation for signal & bkgd estimations, while 'varied pseudo-data' is constructed out of varied MC simulations). If needed, use together with 'usePoissonSmearedPseudoData == 1' and create special filelists beforehand via: scripts/mk_inputForVariedPseudoData_13TeV.sh and scripts/mk_HistoFileListForVariedPseudoData_13TeV.sh
baseFileList("FileLists_"+year+"/"),
outputPlotsFolder("Plots_"+year+"/"),
doDYScaleFromFit(isDYScale && true),
fileReader(rootFileReader),
year_(year),
yearCorr_(yearCorr),
isSingleLeptonDatasetIncluded(true),// supported only for 2016 data
doDYScale(isDYScale),
flagUseZpeak(doDYScaleFromFit && true)
{

    setYearOptions(setDiffXSec, flagProcSyst);

    makeControlPlotEventCountTable = true;
    EventCountTableWhenHist = "HypTTBarpT";

}

PlotterConfigurationHelper::PlotterConfigurationHelper(const TString& year, const TString& yearCorr, bool setDiffXSec, bool flagProcSyst):
year_(year),
yearCorr_(yearCorr)
{
    setYearOptions(setDiffXSec, flagProcSyst);
}

void PlotterConfigurationHelper::setYearOptions(bool setDiffXSec, bool flagProcSyst){
    // Filling vector with systematics excluding the ones with no UP/DOWN variations, and "Nominal" (no UP and Down is added).
    if (setDiffXSec){
        fillVectorOfValidSystematics(_vectorOfSystematics, year_, yearCorr_, true, false, false, false, false, false, false, "all", true);
        // Filling vector with systematics including only the ones with no UP/DOWN variations.
        fillVectorOfValidSystematics(_vectorOfSystematics_noUpDown, year_, yearCorr_, true, false, false, false, false, true, true, "all");
        // Filling vector with systematics excluding "Nominal" (no UP and Down is added).
        fillVectorOfValidSystematics(_vectorOfSystematics_all, year_, yearCorr_, true, false, false, false, false, true, false, "all");
    }
    else{
        fillVectorOfValidSystematics(_vectorOfSystematics, year_, yearCorr_, false, false, false, false, false, false, false, "all", true);
        // Filling vector with systematics including only the ones with no UP/DOWN variations.
        fillVectorOfValidSystematics(_vectorOfSystematics_noUpDown, year_, yearCorr_, false, false, false, false, false, true, true, "all");
        // Filling vector with systematics excluding "Nominal" (no UP and Down is added).
        fillVectorOfValidSystematics(_vectorOfSystematics_all, year_, yearCorr_, false, false, false, false, false, true, false, "all");
    }

    std::pair<double,double> lumi_and_error = get_lumi_and_err(year_);
    lumi = lumi_and_error.first;
    lumiError = lumi_and_error.second;

    if (year_ == "2016preVFP") {

        energy_label = "13 TeV";
        years_list = {year_};

    } else if (year_ == "2016postVFP") {

        energy_label = "13 TeV";
        years_list = {year_};

    } else if (year_ == "2016") {

        energy_label = "13 TeV";
        years_list = {"2016preVFP", "2016postVFP"};

    } else if (year_ == "2017") {

        energy_label = "13 TeV";
        years_list = {year_};

    } else if (year_ == "2018") {

        energy_label = "13 TeV";
        years_list = {year_};

    } else if (year_ == "fullRun2UL") {

        energy_label = "13 TeV";
        years_list = {"2016preVFP", "2016postVFP", "2017", "2018"};

    } else {

        std::cerr << "ERROR in PlotterConfigurationHelper::PlotterConfigurationHelper() : year is not supported " << std::endl;
        exit(1);
    }

    lumi_fb << std::fixed;
    if (year_ == "fullRun2UL") lumi_fb << std::setprecision(0);
    else lumi_fb << std::setprecision(1);
    lumi_fb << lumi/1000;

    if(setDiffXSec){
      this->runningDiffXSec = true;
      setPlotterDiffXSecVars();
    }
    if(flagProcSyst)
      this->flagProcSyst = true;
    else
      this->flagProcSyst = false;

    this->_envelopesMap = EnvelopesMap;

    gSysJES = 1; // 1 for JES , 2 for JES split
}

void PlotterConfigurationHelper::setPlotterDiffXSecVars(){
    // Plotting and output options
    ttmd_createEventSummaryTable = true;
    ttmd_plot_uncertainties_summary = true;
    ttmd_plot_only_Absolute_and_Regularized = false;
    ttmd_plot_only_Normalized_and_Regularized = false;
    ttmd_plot_only_Regularized = true;
    ttmd_plot_comparison_between_channels = false;
    ttmd_plot_comparison_between_years = false;
    ttmd_plot_systs_mig_matrix_comparison_between_years = false;
    ttmd_include_alternative_theory = true;
    ttmd_skip_extra_plots_for_systematics = true;

    // corrections options
    ttmd_apply_trueLevelWeights = true;

    gBaseDir = "../";
    gPlainTreeVersion = year_; // updated simple kin reco: configuration for the initial setup

    longTitleCrossSection = false;

    gTTBARDebug = 0;

    fileList_prefix_plainTrees = "PlainTreeFileList_";
    fileList_format_plainTrees = ".txt";

    // Output Paths
    analysisSuffix = "140119"; // Analysis suffix (Usefull for storing different version of the plots)
    diLeptonic_dir = gBaseDir + "diLeptonic/";
    analysis_folder = "ttmd_analysis_" + year_;
    // gTreeDir  = diLeptonic_dir + "mergedPlainTree_" + year_; // directory with plainTree: configuration for the initial setup
    gAnalDir = diLeptonic_dir + analysis_folder; // directory with final plots
    gPseudoDir = gAnalDir + "cover/"; // directory with toys
    outDir = ""; // runtime active directory
    gCorrData = gAnalDir + "corrdata/";
    gForAN = gAnalDir + "forAN/";
    gMGTabsDir = gAnalDir + "/mgtabs_" + year_;
    gNPCorrSubdir = gAnalDir + "/npcorr_" + year_;

    gBTagEffData = diLeptonic_dir + "BTagEff_"+year_; // directory with BTagEffs

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //int gFlagMG = 1; // 1 Nov17, 3 Feb18
    gFlagMG = 4; // Mar18 -> Sep18
    //int gFlagMG = 5; // May18: mu = H/4, reduced number of histograms
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //
    gPlainTreeName = "plainTree_rec_step8";
    gPlainTreeName0 = "ttBar_treeVariables_step0";
    gPlainTreeName0Plain = "plainTree_gen_step0";
    //
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    gXsecTtbarTwiki = topxsec;// / (0.02277);
    gBrrTtbarToEMu =bRatio[2];//0.02277;
    gBrrTtbarToEMuUncRel =bRatioError;
    //double gXsecTtbarTwiki = 2031.76;// / (0.02277);
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    gPlotPixelSmall = 500;
    gPlotPixelLarge = 750;
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    gSysPDF = 0; //enable inclusion of different PDF systematics (off=0 and on=1, this was set to 1 in implementation from TOP-18-004)
    gOnlyDominantSyst = 0; // 0 skip nothing, 1 skip non-dominant, 2 skip all
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    gVarBkg = 0.3;
    gVarDY = 0.3;
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //std::vector<TString> gvKRName = { "full", "vis", "visE", "ll", "ll+MET", "MW" };
    //std::vector<TString> gvKRSuffix = { "", "nocor", "corcol", "corll", "cormw" };
    //std::vector<bool> gvIsSimple = { 0, 1, 1, 1, 1, 0 };
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    gUnfUncRspTotal = 1;

    _gVFilesData_pattern = {"run"};
    // TODO: Update for handling alternative theory files, check if works due to _alternativeThoery file name!
    _ttbar_signal_pattern = {"ttbarsignalplustau", "ttbarsignalviatau"};
    _ttbar_bgviatau_pattern = {"ttbarbgviatau"};
    _ttbar_bkg_pattern = {"ttbarbg_", "ttbarbg."};
    _DY_pattern = {"dyee", "dymumu", "dytautau"};
    _gVFilesMCbkg_pattern = {"single", "wtolnu","wwtoall","wztoall","zztoall","ttbarW","ttbarZ","ttgjets"};
    _gVFilesMCbkg_pattern.insert(_gVFilesMCbkg_pattern.end(), _DY_pattern.begin(), _DY_pattern.end());
    _gVFilesMCbkg_pattern_ID = {{"single", 4}, {"1050", 9}, {"50inf", 5}, {"wtolnu", 6}, {"wwtoall", 7}, {"wztoall", 7}, {"zztoall", 7}, {"ttbarW", 8}, {"ttbarZ", 8}, {"ttgjets", 8}};

    //Alternative Theory included in plots
    if (ttmd_include_alternative_theory) {
        //sample production folder name
        ttmd_alternative_MC_samples = { "AMCATNLOFXFX", "MADGRAPHMLM" };
        //Name for legends in plots
        ttmd_alternative_MC_samples_legend_names = { "FXFX+PYT", "MLM+PYT" };
    }
    else{
        ttmd_alternative_MC_samples = {};
        ttmd_alternative_MC_samples_legend_names = {};
    }

    std::cout << std::endl;
}

void PlotterConfigurationHelper::fillVectorOfValidSystematics(std::vector<TString>& vect_all, const TString& year, const TString& yearCorr, const bool runningDiffXSec,
                                                              const bool useAdditionalSystematics, const bool fill_as_Up_Down,
                                                              const bool include_all_as_syst, const bool include_Nominal_as_syst,
                                                              const bool include_vars_with_noUpDown, const bool only_include_vars_with_noUpDown,
                                                              TString sys_type, const bool printExclusionWarning){

    if (sys_type != "all" && sys_type != "ExpSys" && !sys_type.BeginsWith("JESSys") && sys_type != "ModWSys" && sys_type != "ModISys" && sys_type != "AltTheorySys"){
      throw std::logic_error(TString::Format("Unsupported systematic type in PlotterConfigurationHelper::fillVectorOfValidSystematics (sys_type = %s)", sys_type.Data()));
    }

    vect_all.clear();

    if (include_vars_with_noUpDown == false and only_include_vars_with_noUpDown == true) {
        std::cerr << "PlotterConfigurationHelper::fillVectorOfValidSystematics -->> Assertion error in the flags!" << std::endl;
    }

    TString up_label = "";
    TString down_label = "";

    if (fill_as_Up_Down){
        up_label = "_UP";
        down_label = "_DOWN";
    }
    if (include_Nominal_as_syst && only_include_vars_with_noUpDown == false){
        vect_all.push_back("Nominal");
    }
    if (include_all_as_syst && only_include_vars_with_noUpDown == false){
        vect_all.push_back("all");
    }

    std::map<TString, TString> VectorOfExcludedSystematics;
    std::map<TString, TString> VectorOfExcludedSystematicsOnlyYear;
    std::map<TString, TString> VectorOfValidSystematicsOnlyYear;
    std::map<TString, TString> VectorOfValidSystematicsOnlyYear_noUpDown;

    //Setting Year dependent vectors
    if (year=="2016preVFP"){
        VectorOfValidSystematicsOnlyYear = VectorOfValidSystematicsOnly2016preVFP;
        VectorOfExcludedSystematicsOnlyYear = VectorOfExcludedSystematicsOnly2016preVFP;
        VectorOfValidSystematicsOnlyYear_noUpDown = VectorOfValidSystematicsOnly2016preVFP_noUpDown;
    }
    else if (year=="2016postVFP"){
        VectorOfValidSystematicsOnlyYear = VectorOfValidSystematicsOnly2016postVFP;
        VectorOfExcludedSystematicsOnlyYear = VectorOfExcludedSystematicsOnly2016postVFP;
        VectorOfValidSystematicsOnlyYear_noUpDown = VectorOfValidSystematicsOnly2016postVFP_noUpDown;
    }
    else if (year=="2016"){
        VectorOfValidSystematicsOnlyYear = VectorOfValidSystematicsOnly2016;
        VectorOfExcludedSystematicsOnlyYear = VectorOfExcludedSystematicsOnly2016;
        VectorOfValidSystematicsOnlyYear_noUpDown = VectorOfValidSystematicsOnly2016_noUpDown;
    }
    else if (year=="2017"){
        VectorOfValidSystematicsOnlyYear = VectorOfValidSystematicsOnly2017;
        VectorOfExcludedSystematicsOnlyYear = VectorOfExcludedSystematicsOnly2017;
        VectorOfValidSystematicsOnlyYear_noUpDown = VectorOfValidSystematicsOnly2017_noUpDown;
    }
    else if (year=="2018"){
        VectorOfValidSystematicsOnlyYear = VectorOfValidSystematicsOnly2018;
        VectorOfExcludedSystematicsOnlyYear = VectorOfExcludedSystematicsOnly2018;
        VectorOfValidSystematicsOnlyYear_noUpDown = VectorOfValidSystematicsOnly2018_noUpDown;
    }
    else if (year=="fullRun2UL"){
        VectorOfValidSystematicsOnlyYear = VectorOfValidSystematicsOnlyFullRun2UL;
        VectorOfExcludedSystematicsOnlyYear = VectorOfExcludedSystematicsOnlyFullRun2UL;
    }
    else{
        std::cerr << "ERROR in PlotterConfigurationHelper::fillVectorOfValidSystematics() : year is not supported " << std::endl;
    }

    //Print warning if systematics to be excluded were found
    if (VectorOfExcludedSystematicsOnlyYear.size()>0){
        if (printExclusionWarning){
            std::cout << " \n	WARNING! These systematics have been excluded for " + year + ":" << std::endl;
            for (auto s: VectorOfExcludedSystematicsOnlyYear) std::cout << "		" << s.first << std::endl;
            std::cout << std::endl;
        }
        VectorOfExcludedSystematics.insert(VectorOfExcludedSystematicsOnlyYear.begin(), VectorOfExcludedSystematicsOnlyYear.end());
    }

    if (VectorOfExcludedSysts_ttmd_only.size()>0 && runningDiffXSec){
        if (printExclusionWarning){
            std::cout << " \n	WARNING! These systematics have been excluded for DiffXSec analysis:" << std::endl;
            for (auto s: VectorOfExcludedSysts_ttmd_only) std::cout << "		" << s.first << std::endl;
            std::cout << std::endl;
        }
        VectorOfExcludedSystematics.insert(VectorOfExcludedSysts_ttmd_only.begin(), VectorOfExcludedSysts_ttmd_only.end());
    }

    //Adding the common systematics excluding the excluded Systematics in VectorOfExcludedSystematics (systs without _UP and _DOWN).
    if (include_vars_with_noUpDown){
        for (auto s: VectorOfValidSystematicsCommon_noUpDown) {
                //Setting the type of the sys in the loop to all if sys_type == "all"
                TString this_var_type = "all";
                if (sys_type != "all") this_var_type = s.second;
                if (std::find(VectorOfExcludedSystematics.begin(), VectorOfExcludedSystematics.end(), s) == VectorOfExcludedSystematics.end() && this_var_type == sys_type) {
		  std::vector<TString> yearCorrList = GetYearCorrList(year, yearCorr, s.first, IsYearToYearCorr(year, yearCorr, s.first));
		    if (yearCorrList.size() != 0)
		      vect_all.insert(vect_all.end(), yearCorrList.begin(), yearCorrList.end());
		    else
		      vect_all.push_back(s.first);
		}
        }
    }

    if (only_include_vars_with_noUpDown == false) {
        //Adding the common systematics excluding the excluded Systematics in VectorOfExcludedSystematics (_UP and _DOWN).
        for (auto s: VectorOfValidSystematicsCommon) {
                //Setting the type of the sys in the loop to all if sys_type == "all"
                TString this_var_type = "all";
                if (sys_type != "all") this_var_type = s.second;

                if (std::find(VectorOfExcludedSystematics.begin(), VectorOfExcludedSystematics.end(), s) == VectorOfExcludedSystematics.end() && this_var_type == sys_type) {
		    std::vector<TString> yearCorrList = GetYearCorrList(year, yearCorr, s.first, IsYearToYearCorr(year, yearCorr, s.first));
                    if (fill_as_Up_Down){
		        if (yearCorrList.size() != 0) {
			  for(auto it = yearCorrList.begin(); it != yearCorrList.end(); ++it) {
			    vect_all.push_back(*it + up_label);
			    vect_all.push_back(*it + down_label);
			  }
			}
			else {
			  vect_all.push_back(s.first + up_label);
			  vect_all.push_back(s.first + down_label);
			}
                    }
                    else {
		        if (yearCorrList.size() != 0) {
			  vect_all.insert(vect_all.end(), yearCorrList.begin(), yearCorrList.end());
			}
			else
			  vect_all.push_back(s.first);
                    }
		}
        }
    }

    //Adding specific year systematics.
    if (include_vars_with_noUpDown){
      for (auto s: VectorOfValidSystematicsOnlyYear_noUpDown) {
        //Setting the type of the sys in the loop to all if sys_type == "all"
        TString this_var_type = "all";
        if (sys_type != "all") this_var_type = s.second;
        if (std::find(VectorOfExcludedSystematics.begin(), VectorOfExcludedSystematics.end(), s) == VectorOfExcludedSystematics.end() && this_var_type == sys_type){
	  std::vector<TString> yearCorrList = GetYearCorrList(year, yearCorr, s.first, IsYearToYearCorr(year, yearCorr, s.first));
	  if (yearCorrList.size() != 0)
	    vect_all.insert(vect_all.end(), yearCorrList.begin(), yearCorrList.end());
	  else
	    vect_all.push_back(s.first);
	}
      }
    }

    if (only_include_vars_with_noUpDown == false) {
        for (auto s: VectorOfValidSystematicsOnlyYear) {
            //Setting the type of the sys in the loop to all if sys_type == "all"
            TString this_var_type = "all";
            if (sys_type != "all") this_var_type = s.second;
	    std::vector<TString> yearCorrList = GetYearCorrList(year, yearCorr, s.first, IsYearToYearCorr(year, yearCorr, s.first));
	    if (fill_as_Up_Down && this_var_type == sys_type){
	      if (yearCorrList.size() != 0) {
		for(auto it = yearCorrList.begin(); it != yearCorrList.end(); ++it) {
		  vect_all.push_back(*it + up_label);
		  vect_all.push_back(*it + down_label);
		}
	      }
	      else {
		vect_all.push_back(s.first + up_label);
		vect_all.push_back(s.first + down_label);
	      }
	    }
	    else if (this_var_type == sys_type) {
	      if (yearCorrList.size() != 0) {
		vect_all.insert(vect_all.end(), yearCorrList.begin(), yearCorrList.end());
	      }
	      else
		vect_all.push_back(s.first);
	    }
        }
    }

    //Adding additional systematics.
    if (useAdditionalSystematics) {
        if (only_include_vars_with_noUpDown == false) {
            for (auto s: VectorOfValidSystematicsAddons) {
                //Setting the type of the sys in the loop to all if sys_type == "all"
                TString this_var_type = "all";
                if (sys_type != "all") this_var_type = s.second;
                if (std::find(begin(VectorOfValidSystematicsCommon), end(VectorOfValidSystematicsCommon), s) == end(VectorOfValidSystematicsCommon)) {
		  std::vector<TString> yearCorrList = GetYearCorrList(year, yearCorr, s.first, IsYearToYearCorr(year, yearCorr, s.first));
		  if (fill_as_Up_Down && this_var_type == sys_type){
		    if (yearCorrList.size() != 0) {
		      for(auto it = yearCorrList.begin(); it != yearCorrList.end(); ++it) {
			vect_all.push_back(*it + up_label);
			vect_all.push_back(*it + down_label);
		      }
		    }
		    else {
		      vect_all.push_back(s.first + up_label);
		      vect_all.push_back(s.first + down_label);
		    }
		  }
		  else if (this_var_type == sys_type) {
		    if (yearCorrList.size() != 0) {
		      vect_all.insert(vect_all.end(), yearCorrList.begin(), yearCorrList.end());
		    }
		    else
		      vect_all.push_back(s.first);
		  }
                }
            }
        }
        if (include_vars_with_noUpDown){
            for (auto s: VectorOfValidSystematicsAddons_noUpDown) {
                //Setting the type of the sys in the loop to all if sys_type == "all"
                TString this_var_type = "all";
                if (sys_type != "all") this_var_type = s.second;
                if (std::find(begin(VectorOfValidSystematicsCommon), end(VectorOfValidSystematicsCommon), s) == end(VectorOfValidSystematicsCommon) && this_var_type == sys_type) {
		  std::vector<TString> yearCorrList = GetYearCorrList(year, yearCorr, s.first, IsYearToYearCorr(year, yearCorr, s.first));
		  if (yearCorrList.size() != 0)
		    vect_all.insert(vect_all.end(), yearCorrList.begin(), yearCorrList.end());
		  else
		    vect_all.push_back(s.first);
		}
            }
        }

    }

}

double PlotterConfigurationHelper::SampleXSection(const TString& filename)
{
    //13 TeV MC cross sections taken from:
    //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
    if (filename.Contains("run"))                   {if (usePoissonSmearedPseudoData) return topxsec;
                                                     else return 1.;}
    else if (filename.Contains("ttbarsignal") && filename.Contains("AMCATNLOFXFX"))   {return topxsec;}// updated PDG 2020 with Dilepton BR Corrections
    else if (filename.Contains("ttbarsignal") && filename.Contains("POWHEGV2HERWIG"))   {return topxsec;}// updated PDG 2020 with Dilepton BR Corrections
    else if (filename.Contains("ttbarsignal") && filename.Contains("MADGRAPHMLM"))   {return topxsec * 0.10706;}// updated PDG 2020 with Dilepton BR Corrections
    else if (filename.Contains("ee_ttbarsignalplustau_fromDilepton"))   {return topxsec * 0.10706 * 0.964261576;}// updated PDG 2020 with Dilepton BR Corrections
    else if (filename.Contains("emu_ttbarsignalplustau_fromDilepton"))  {return topxsec * 0.10706 * 0.957058875;}// updated PDG 2020 with Dilepton BR Corrections
    else if (filename.Contains("mumu_ttbarsignalplustau_fromDilepton")) {return topxsec * 0.10706 * 0.949909976;}// updated PDG 2020 with Dilepton BR Corrections
    else if (filename.Contains("ee_ttbarsignalviatau_fromDilepton"))    {return topxsec * 0.10706 * 1.029827957;}// updated PDG 2020 with Dilepton BR Corrections
    else if (filename.Contains("emu_ttbarsignalviatau_fromDilepton"))   {return topxsec * 0.10706 * 1.026209047;}// updated PDG 2020 with Dilepton BR Corrections
    else if (filename.Contains("mumu_ttbarsignalviatau_fromDilepton"))  {return topxsec * 0.10706 * 1.022670477;}// updated PDG 2020 with Dilepton BR Corrections
    else if (filename.Contains("ttbarbg_fromDilepton"))       {return topxsec * 0.10706;}// updated PDG 2020
    else if (filename.Contains("ttbarbg_fromLjets"))          {return topxsec * 0.44113;}// updated PDG 2020
    else if (filename.Contains("ttbarbg_fromHadronic"))       {return topxsec * 0.45441;}// updated PDG 2020

    else if (filename.Contains("ttbar") &&
            !filename.Contains("ttbarW") &&
            !filename.Contains("ttbarZ"))             {return topxsec;}// updated
    else if (filename.Contains("single") &&
	     //NoFullyHadronic
	     filename.Contains("tw"))                 {return 35.85*0.54559;}// updated
             //Inclusive
    //	     filename.Contains("tw"))                 {return 35.85;}// updated
    else if (filename.Contains("single") &&
	     filename.Contains("_s"))                 {return 10.32;}// updated
    else if (filename.Contains("singletop") &&
	     filename.Contains("_t"))                 {return 136.02;}// updated
    else if (filename.Contains("singleantitop") &&
	     filename.Contains("_t"))                 {return 80.95;}// updated
    else if (filename.Contains("ww"))                 {return 118.7;}// updated
    else if (filename.Contains("wz"))                 {return 47.13;}// updated 04.09.17: NLO from MCFM
    else if (filename.Contains("zz"))                 {return 16.523;}// updated 04.09.17: NLO from MCFM
    else if (filename.Contains("1050"))               {return 22635.1;}// from NNLO FEWZ 3.1

    // k factors from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#notes
    else if (filename.Contains("50inf_ht0040to0070"))    {return 310.7*1.23;}// updated 13.05.20 / GenXSecAnalyzer
    else if (filename.Contains("50inf_ht0070to0100"))    {return 169.9*1.23;}// updated 22.05.19
    else if (filename.Contains("50inf_ht0100to0200"))    {return 147.40*1.23;}// updated 22.05.19
    else if (filename.Contains("50inf_ht0200to0400"))    {return 40.99*1.23;}// updated 22.05.19
    else if (filename.Contains("50inf_ht0400to0600"))    {return 5.678*1.23;}// updated 22.05.19
    else if (filename.Contains("50inf_ht0600to0800"))    {return 1.367*1.23;}//////TODO
    else if (filename.Contains("50inf_ht0800to1200"))    {return 0.6304*1.23;}// updated 22.05.19
    else if (filename.Contains("50inf_ht1200to2500"))    {return 0.1514*1.23;}// updated 22.05.19
    else if (filename.Contains("50inf_ht2500toINFT"))    {return 0.003565*1.23;}// updated 22.05.19

    // Splitting the 50inf half between the amcatnlofxfx and madgraphmlm
    else if (filename.Contains("50inf_amcatnlofxfx"))               {return 0.5 * 3.*2075.14;}// updated 23.05.19: NNLO FEWZ 3.1.b2 NNPDF31_nnlo_as_0118_luxqed
    else if (filename.Contains("50inf_madgraphmlm"))               {return 0.5 * 3.*2075.14;}// updated 23.05.19: NNLO FEWZ 3.1.b2 NNPDF31_nnlo_as_0118_luxqed

    else if (filename.Contains("wtolnu"))             {return 61526.7;}// updated 22.05.19
    else if (filename.Contains("ttgjets"))            {return 3.697;}// updated 23.05.19 NLO

    else if (filename.Contains("ttbarWjetstolnu"))    {return 0.2043;}// updated: nlo
    else if (filename.Contains("ttbarWjetstoqq"))     {return 0.4062;}// updated: nlo
    else if (filename.Contains("ttbarZtollnunu"))     {return 0.2529;}// updated: nlo
    else if (filename.Contains("ttbarZtoqq"))         {return 0.5297;}// updated: nlo

    return -1;
}


void PlotterConfigurationHelper::fillLegendColorDataset(const TString& fileListName, std::vector<TString>& legends, std::vector<int>& colors, std::vector<TString>& dataset)
{
    std::cout << "reading " << fileListName << std::endl;
    std::ifstream FileList(fileListName);
    if (FileList.fail()) {
        std::cerr << "Error reading " << fileListName << std::endl;
        exit(1);
    }

    dataset.clear();
    legends.clear();
    colors.clear();

    TString filename;
    while (!FileList.eof()) {
        FileList >> filename;
        if (filename == "") continue;//skip empty lines
        dataset.push_back(filename);
        if (filename.Contains("run"))               {legends.push_back("Data"); colors.push_back(kBlack);}
        else if (filename.Contains("ttbarsignal"))  {legends.push_back("t#bar{t} signal"); colors.push_back(kRed+1);}
        else if (filename.Contains("ttbarbg"))      {legends.push_back("t#bar{t} other"); colors.push_back(kRed-7);}
        else if (filename.Contains("single"))       {legends.push_back("Single t"); colors.push_back(kMagenta);}
        else if (filename.Contains("ww") ||
                 filename.Contains("wz") ||
                 filename.Contains("zz"))           {legends.push_back("Diboson"); colors.push_back(10);}
        //else if (filename.Contains("dytautau"))     {legends.push_back("Z / #gamma* #rightarrow #tau#tau"); colors.push_back(kAzure+8);}
        else if (filename.Contains("dymumu") ||
                 filename.Contains("dyee") ||
                 filename.Contains("dytautau"))     {legends.push_back("Z+jets"); colors.push_back(kAzure-2);}
        else if (filename.Contains("wtolnu"))       {legends.push_back("W+jets"); colors.push_back(kGreen-3);}
        //else if (filename.Contains("qcd"))          {legends.push_back("QCD multijet"); colors.push_back(kYellow);}
        else if (filename.Contains("ttbarZ") ||
                 filename.Contains("ttbarW") ||
				 //temporal position for ttgjets
				 filename.Contains("ttgjets"))       {legends.push_back("t#bar{t}+Z/W"); colors.push_back(kOrange-2);}
        else{
        	std::cerr << "File name " << filename << "couldn't be classified into legends!!" << std::endl;
        }
    }
    FileList.close();
}


double PlotterConfigurationHelper::CalcLumiWeight(const TString& WhichSample, const TString in_year)
{
    if (WhichSample.Contains("run") && !usePoissonSmearedPseudoData) return 1;
    double lumiWeight = 0.;
    double lumi_to_use = 0.;

    //    std::cout << WhichSample << std::endl;

    //Load lumi for input year or use default from PlotterConfigurationHelper() if no year is entered
    if (in_year == "use initialized"){
        lumi_to_use = lumi;
    }
    else{
        std::pair<double,double> lumi_and_error = get_lumi_and_err(in_year);
        lumi_to_use = lumi_and_error.first;
    }

    if (WhichSample != "") {
        double XSection = SampleXSection(WhichSample);
        if (XSection <= 0.) {
        	std::cout << WhichSample << " sample XSection is <0. Can't calculate luminosity weight!! returning" << std::endl;
            return 0;
        }
        //From 'filename' get the number of weighted (MC weights) event processed.
        const TH1 *h_NrOfEvts = fileReader->Get<TH1>(WhichSample, "weightedEvents");
        double NrOfEvts = h_NrOfEvts->GetBinContent(1);
        if (!WhichSample.Contains("run")) lumiWeight = lumi_to_use * XSection / NrOfEvts;
        if (WhichSample.Contains("run") && usePoissonSmearedPseudoData) lumiWeight = 2 * lumi_to_use * XSection / NrOfEvts; //if ttbarsignlaplustau is merged with ttbarbgviatau --> 2
        //std::cout << WhichSample << " " << lumiWeight << " xsec:" << XSection << std::endl;
        //std::cout << WhichSample <<", lumiw=" << lumiWeight <<", lumi=" << lumi <<", XSection=" << XSection <<", NrOfEvts=" << NrOfEvts << "\n";
    }

    if (lumiWeight == 0 && WhichSample != "") {
        std::cout << WhichSample << " has lumi weight 0\n";
    }

    return lumiWeight;
}


void PlotterConfigurationHelper::ApplyFlatWeights(TH1* varhists, const double weight)
{
    if (weight == 0) std::cout << "Warning: the weight your applying is 0. This will remove your distribution." << std::endl;
    varhists->Scale(weight);
}


void PlotterConfigurationHelper::ApplyFlatWeights(TH2* varhists, const double weight)
{
    if (weight == 0) std::cout<<"Warning: the weight your applying is 0. This will remove your distribution."<<std::endl;
    varhists->Scale(weight);
}

void PlotterConfigurationHelper::RescaleHistYearToYearCorr(TH1* varHist, TH1* nominalHist, const double year_corr_rescale_from_file)
{
  if (year_corr_rescale_from_file != 1.) {
    varHist->Add(nominalHist, -1.0);
    varHist->Scale(year_corr_rescale_from_file);
    varHist->Add(nominalHist, 1.0);
  }
}

void PlotterConfigurationHelper::RescaleHistYearToYearCorr(TH2* varHist, TH2* nominalHist, const double year_corr_rescale_from_file)
{
  if (year_corr_rescale_from_file != 1.) {
    varHist->Add(nominalHist, -1.0);
    varHist->Scale(year_corr_rescale_from_file);
    varHist->Add(nominalHist, 1.0);
  }
}

void PlotterConfigurationHelper::DYScaleFactor(TString SpecialComment,std::vector<double>& DYScale, TString name)
{
    DYScale = {1.,1.,1.,1.};

    if (!doDYScale || usePoissonSmearedPseudoData) return;

    TString nameAppendix = "";
    if (!SpecialComment.BeginsWith("_post") &&  SpecialComment != "Standard" ) {
        std::cout << "\n\n*******************************************************************" << std::endl;
        std::cout << "ERROR: When calculating the DY Scale factor you must specify in which step you want to calculate the DY SF:" << std::endl;
        std::cout << " '_postZcut', '_post2jets', '_postMET', '_post1btag', '_postKinReco' or 'Standard' = _postKinReco" << std::endl;
        std::cout << "*******************************************************************\n\n" << std::endl;
        exit(444);
    }
    if (SpecialComment.BeginsWith("_post")) {
        if (SpecialComment.EqualTo("_postZcut")) nameAppendix = "_step4";
        if (SpecialComment.EqualTo("_post2jets")) nameAppendix = "_step5";
        if (SpecialComment.EqualTo("_postMET")) nameAppendix = "_step6";
        if (SpecialComment.EqualTo("_post1btag")) nameAppendix = "_step7";
        if (SpecialComment.EqualTo("_postKinReco")) nameAppendix = "_step8";
    } else if (SpecialComment == "Standard") {
        nameAppendix = "_step5";
    }

    TString latexLabel = "";
    if (nameAppendix.EqualTo("_step4")) latexLabel = "Z-peak";
    if (nameAppendix.EqualTo("_step5")) latexLabel = "2 jets";
    if (nameAppendix.EqualTo("_step6")) latexLabel = "$E_\\text{T}^\\text{miss}$";
    if (nameAppendix.EqualTo("_step7")) latexLabel = "b-tag";
    if (nameAppendix.EqualTo("_step8")) latexLabel = "kin. fit";

    std::cout << "\n\nBegin DYSCALE FACTOR calculation at selection step " << nameAppendix << std::endl;

    //std::vector<TString> Vec_Files = InputNominalFileList("combined");
	std::vector<TString> Vec_Files;
	std::vector<TString> Vec_Files_root_names;
	TString histoListName = baseFileList + "HistoFileList_Nominal_combined.txt";
	std::ifstream FileList(histoListName);
	TString filename;
	while (!FileList.eof()) {
		FileList >> filename;
		//std::cout << filename << std::endl;
		TString nameRootFile = ExtractFileNameFromPath(filename,false);
		Vec_Files.push_back(filename);
		Vec_Files_root_names.push_back(nameRootFile);
	}

    if (Vec_Files.size() < 1) {
        std::cout << "WARNING(in DYScaleFactor)!!! No datasets available to calculate DY SF. EXITING!!"<<std::endl;
        return;
    }

    double NoutEEDYMC = 0, NinEEDYMC = 0, NoutMuMuDYMC = 0, NinMuMuDYMC = 0;//number of events in/out of z-veto region for the DY MC
    double NoutEMuDYMC = 0;//same as above, but for EMU channel, needed for correct combination of SFs for the "combined" channel
    double NinEE = 0, NinMuMu = 0, NinEMu = 0;//number of events in z-veto region for data
    double NinEEloose = 0, NinMuMuloose = 0;//number of data events in Z-Veto region with MET cut
    double NinEEMC = 0, NinMuMuMC = 0;//all other MC events

    for (size_t i = 0; i < Vec_Files.size(); i++) {
        double LumiWeight;
	if (year_ == "fullRun2UL" || year_ == "2016") {
	  get_year_from_file(Vec_Files.at(i));
	  LumiWeight = CalcLumiWeight(Vec_Files.at(i), lumi_year);
	}
	else
	  LumiWeight = CalcLumiWeight(Vec_Files.at(i));
        double allWeights = LumiWeight;//calculate here all the flat-weights we apply: Lumi*others*...
        if (Vec_Files_root_names.at(i).Contains("ee_") || Vec_Files_root_names.at(i).Contains("mumu_")) {
	    if (Vec_Files_root_names.at(i).Contains("run")) {
                TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString(TString("dyScaling_Zh1").Append(nameAppendix)).Append("zWindow"));
                TH1D *htemp1 = fileReader->GetClone<TH1D>(Vec_Files.at(i), "dyScaling_Looseh1");
                ApplyFlatWeights(htemp, allWeights);
                ApplyFlatWeights(htemp1, allWeights);
                if (Vec_Files_root_names.at(i).Contains("ee_")) {
		    printf("Vec_Files.at(i) = %s \n", Vec_Files.at(i).Data());
                    NinEE += htemp->Integral();
                    NinEEloose += htemp1->Integral();
                }
                if (Vec_Files_root_names.at(i).Contains("mumu_")) {
                    NinMuMu += htemp->Integral();
                    NinMuMuloose += htemp1->Integral();
                }
                delete htemp;
                delete htemp1;
            } else if (Vec_Files_root_names.at(i).Contains("dy")) {
                if (Vec_Files_root_names.at(i).Contains("50inf") || Vec_Files_root_names.at(i).Contains("50inf_ht")) {
                    TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString(TString("dyScaling_Zh1").Append(nameAppendix)).Append("zWindow"));
                    TH1D *htemp1 = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString("dyScaling_TTh1").Append(nameAppendix));
                    ApplyFlatWeights(htemp, LumiWeight);
                    ApplyFlatWeights(htemp1, LumiWeight);
                    if (Vec_Files_root_names.at(i).Contains("ee_")) {
                        NinEEDYMC += htemp->Integral();
                        NoutEEDYMC += htemp1->Integral();
                    }
                    if (Vec_Files_root_names.at(i).Contains("mumu_")) {
                        NinMuMuDYMC += htemp->Integral();
                        NoutMuMuDYMC += htemp1->Integral();
                    }
                    delete htemp;
                    delete htemp1;
                } else {
                    TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString("dyScaling_TTh1").Append(nameAppendix));
                    ApplyFlatWeights(htemp, LumiWeight);
                    if (Vec_Files_root_names.at(i).Contains("ee_")) NoutEEDYMC += htemp->Integral();
                    if (Vec_Files_root_names.at(i).Contains("mumu_")) NoutMuMuDYMC += htemp->Integral();
                    delete htemp;
                }
            } else {
                TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString(TString("dyScaling_Zh1").Append(nameAppendix)).Append("zWindow"));
                ApplyFlatWeights(htemp, LumiWeight);
                if (Vec_Files_root_names.at(i).Contains("ee_")) NinEEMC += htemp->Integral();
                if (Vec_Files_root_names.at(i).Contains("mumu_")) NinMuMuMC += htemp->Integral();
                delete htemp;
            }
        }

        if (Vec_Files_root_names.at(i).Contains("emu_")) {
            if (Vec_Files_root_names.at(i).Contains("run")) {
                TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString(TString("dyScaling_Zh1").Append(nameAppendix)).Append("zWindow"));
                ApplyFlatWeights(htemp, LumiWeight);
                NinEMu += htemp->Integral();
                delete htemp;
            }
            else if (Vec_Files_root_names.at(i).Contains("dy")) {
                if (Vec_Files_root_names.at(i).Contains("50inf") || Vec_Files_root_names.at(i).Contains("50inf_ht")) {
                    TH1D *htemp1 = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString("dyScaling_TTh1").Append(nameAppendix));
                    ApplyFlatWeights(htemp1, LumiWeight);
                    NoutEMuDYMC += htemp1->Integral();
                    delete htemp1;
                } else {
                    TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString("dyScaling_TTh1").Append(nameAppendix));
                    ApplyFlatWeights(htemp, LumiWeight);
                    NoutEMuDYMC += htemp->Integral();
                    delete htemp;
                }
            }
        }
    }

    double kee = sqrt(NinEEloose / NinMuMuloose);
    double kmumu = sqrt(NinMuMuloose / NinEEloose);
    double RoutinEE = NoutEEDYMC / NinEEDYMC;
    double RoutinMuMu = NoutMuMuDYMC / NinMuMuDYMC;
    double NoutMCEE = RoutinEE * (NinEE - 0.5 * NinEMu * kee);
    double NoutMCMuMu = RoutinMuMu * (NinMuMu - 0.5 * NinEMu * kmumu);

    double DYSFEE = NoutMCEE / NoutEEDYMC;
    double DYSFMuMu = NoutMCMuMu / NoutMuMuDYMC;
    double DYSFEMu = std::sqrt(DYSFEE * DYSFMuMu);
    double DYSFDilept = (DYSFEE * NoutEEDYMC + DYSFMuMu * NoutMuMuDYMC + DYSFEMu * NoutEMuDYMC) / (NoutEEDYMC + NoutMuMuDYMC + NoutEMuDYMC);

    std::cout << std::endl;
    std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;
    std::cout << "Calculation of DY Scale Factors for '" << name << "'  at selection step " << nameAppendix << std::endl;
    std::cout << "DYSFEE:                 " << DYSFEE << std::endl;
    std::cout << "DYSFMuMu:               " << DYSFMuMu << std::endl;
    std::cout << "DYSFEMu:                " << DYSFEMu << std::endl;
    std::cout << "DYSFDilept:             " << DYSFDilept << std::endl;
    std::cout << "NinEEloose:             " << NinEEloose << std::endl;
    std::cout << "NinMMloose:             " << NinMuMuloose << std::endl;
    std::cout << "kee:                    " << kee << " +- "<< kee * 0.5 * TMath::Sqrt(1. / NinMuMuloose + 1. / NinEEloose) << std::endl;
    std::cout << "kmumu:                  " << kmumu << " +- "<< kmumu * 0.5 * TMath::Sqrt(1. / NinMuMuloose + 1. / NinEEloose) << std::endl;
    std::cout << "Rout/Rin ee:            " << RoutinEE << std::endl;
    std::cout << "Rout/Rin Mumu:          " << RoutinMuMu << std::endl;
    std::cout << "Est. From Data(ee):     " << NoutMCEE << std::endl;
    std::cout << "Est. From Data(mumu):   " << NoutMCMuMu << std::endl;
    std::cout << "Est. From MC(ee):       " << NoutEEDYMC << std::endl;
    std::cout << "Est. From MC(mumu):     " << NoutMuMuDYMC << std::endl;
    std::cout << "Est. From MC(emu):      " << NoutEMuDYMC << std::endl;
    std::cout << "      " << latexLabel << " &  " << round(DYSFEE * 1000)/1000 << " & " << round(DYSFMuMu * 1000)/1000 << " & " << round(DYSFEMu * 1000)/1000 << " & " << round(DYSFDilept * 1000)/1000 << " \\\\" << std::endl;
    std::cout << "\\hline "+year_+" "+latexLabel+" & \\ee & \\mumu & \\emu \\\\" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$\\text{\\Zjets}_{MC}$ & " << NoutEEDYMC << " & " << NoutMuMuDYMC << " & " << NoutEMuDYMC << " \\\\" << std::endl;
    std::cout << "$\\text{\\Zjets}_{DATA}$ & " << NoutMCEE << " & " << NoutMCMuMu << " & -\\\\" << std::endl;
    std::cout << "$N^{l^{+}l^{-}}_{in,~Data}$ & " << NinEE << " & " << NinMuMu << " & " << NinEMu << " \\\\" << std::endl;
    std::cout << "$R_{out/in}$ & " << RoutinEE << " & " << RoutinMuMu << " & -\\\\" << std::endl;
    std::cout << "${\\cal C}_{\\rm{Z+jets}}$ & " << round(DYSFEE * 1000)/1000 << " & " << round(DYSFMuMu * 1000)/1000 << " & " << round(DYSFEMu * 1000)/1000 << " \\\\" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;
    std::cout << std::endl;

    if (DYSFEE < 0.2 || DYSFMuMu < 0.2) {
        std::cout << "The DY SF is too low (below 0.2). Something is probably wrong.\n";
        std::exit(1);
    }


    DYScale.at(0) = DYSFEE;
    DYScale.at(1) = DYSFMuMu;
    DYScale.at(2) = DYSFEMu;
    DYScale.at(3) = DYSFDilept;


    std::cout << "End DYSCALE FACTOR calculation\n" << std::endl;
}


void PlotterConfigurationHelper::DYScaleFactorFromFit(TString SpecialComment,std::vector<double>& DYScale, TString name)
{
    TH1::AddDirectory(kFALSE);
    DYScale = {1.,1.,1.,1.};

    if (!doDYScale || usePoissonSmearedPseudoData) return;

    TString nameAppendix = "";
    if (!SpecialComment.BeginsWith("_post") &&  SpecialComment != "Standard" ) {
        std::cout << "\n\n*******************************************************************" << std::endl;
        std::cout << "ERROR: When calculating the DY Scale factor you must specify in which step you want to calculate the DY SF:" << std::endl;
        std::cout << " '_postZcut', '_post2jets', '_postMET', '_post1btag', '_postKinReco' or 'Standard' = _post2jets" << std::endl;
        std::cout << "*******************************************************************\n\n" << std::endl;
        exit(444);
    }
    if (SpecialComment.BeginsWith("_post")) {
        if (SpecialComment.EqualTo("_postZcut")) nameAppendix = "_step4";
        if (SpecialComment.EqualTo("_post2jets")) nameAppendix = "_step5";
        if (SpecialComment.EqualTo("_postMET")) nameAppendix = "_step6";
        if (SpecialComment.EqualTo("_post1btag")) nameAppendix = "_step7";
        if (SpecialComment.EqualTo("_postKinReco")) nameAppendix = "_step8";
    } else if (SpecialComment == "Standard") {
        nameAppendix = "_step5";
    }

    TString latexLabel = "";
    if (nameAppendix.EqualTo("_step4")) latexLabel = "Z-peak";
    if (nameAppendix.EqualTo("_step5")) latexLabel = "2 jets";
    if (nameAppendix.EqualTo("_step6")) latexLabel = "$E_\\text{T}^\\text{miss}$";
    if (nameAppendix.EqualTo("_step7")) latexLabel = "b-tag";
    if (nameAppendix.EqualTo("_step8")) latexLabel = "kin. fit";
    
    std::cout << "\n\nBegin DYSCALE FACTOR calculation at selection step " << nameAppendix << std::endl;

    if (!flagUseZpeak)
      std::cout << "Calculating DY scale factors from fit of whole dilepton mass distribution" << std::endl;
    else
      std::cout << "Calculating DY scale factors from fit of exclusion region in dilepton mass distribution" << std::endl;

    std::vector<TString> Vec_Files;
    std::vector<TString> Vec_Files_root_names;
    TString histoListName = baseFileList + "HistoFileList_Nominal_combined.txt";
    std::ifstream FileList(histoListName);
    TString filename;
    while (!FileList.eof()) {
      FileList >> filename;
      //      std::cout << filename << std::endl;
      TString nameRootFile = ExtractFileNameFromPath(filename,false);
      Vec_Files.push_back(filename);
      Vec_Files_root_names.push_back(nameRootFile);
    }

    if (Vec_Files.size() < 1) {
        std::cout << "WARNING(in DYScaleFactor)!!! No datasets available to calculate DY SF. EXITING!!"<<std::endl;
        return;
    }

    TH1D *inputData_ee = nullptr;
    TH1D *templateDY_ee = nullptr;
    TH1D *templateOtherMC_ee = nullptr;

    TH1D *inputData_mumu = nullptr;
    TH1D *templateDY_mumu = nullptr;
    TH1D *templateOtherMC_mumu = nullptr;

    TH1D *inputData_emu = nullptr;
    //uncomment to perform fit for emu
    //TH1D *templateDY_emu = nullptr;
    //TH1D *templateOtherMC_emu = nullptr;

    for (size_t i = 0; i < Vec_Files.size(); i++) {
      double LumiWeight;

      //      std::cout << Vec_Files.at(i) << std::endl;

      if (year_ == "fullRun2UL" || year_ == "2016") {
	get_year_from_file(Vec_Files.at(i));
	LumiWeight = CalcLumiWeight(Vec_Files.at(i), lumi_year);
      }
      else
	LumiWeight = CalcLumiWeight(Vec_Files.at(i));
      double allWeights = LumiWeight;//calculate here all the flat-weights we apply: Lumi*others*...
      if (Vec_Files_root_names.at(i).Contains("ee_") || Vec_Files_root_names.at(i).Contains("mumu_") || Vec_Files_root_names.at(i).Contains("emu_")) {
	if (Vec_Files_root_names.at(i).Contains("run")) {
	  TH1D *inputData;
	  if (!flagUseZpeak)
	    inputData = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString(TString("dyScaling_Allh1")).Append(nameAppendix));
	  else
	    inputData = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString(TString("dyScaling_Zh1").Append(nameAppendix)).Append("zWindow"));
    	  ApplyFlatWeights(inputData, allWeights);
	  if (!inputData_ee && Vec_Files_root_names.at(i).Contains("ee_")) {
    	    inputData_ee = (TH1D*) inputData->Clone();
	    inputData_ee->SetDirectory(0);
    	  }
	  else if (Vec_Files_root_names.at(i).Contains("ee_")) {
    	    inputData_ee->Add(inputData);
    	  }
	  if (!inputData_mumu && Vec_Files_root_names.at(i).Contains("mumu_")) {
    	    inputData_mumu = (TH1D*) inputData->Clone();
	    inputData_mumu->SetDirectory(0);
    	  }
	  else if (Vec_Files_root_names.at(i).Contains("mumu_")) {
    	    inputData_mumu->Add(inputData);
    	  }
	  if (!inputData_emu && Vec_Files_root_names.at(i).Contains("emu_")) {
    	    inputData_emu = (TH1D*) inputData->Clone();
	    inputData_emu->SetDirectory(0);
    	  }
	  else if (Vec_Files_root_names.at(i).Contains("emu_")) {
    	    inputData_emu->Add(inputData);
    	  }
    	  delete inputData;
	} else if (Vec_Files_root_names.at(i).Contains("dy")) {
	  TH1D *templateDY;
	  if (!flagUseZpeak)
	    templateDY = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString(TString("dyScaling_Allh1")).Append(nameAppendix));
	  else
	    templateDY = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString(TString("dyScaling_Zh1").Append(nameAppendix)).Append("zWindow"));
    	  ApplyFlatWeights(templateDY, LumiWeight);
	  if (!templateDY_ee && Vec_Files_root_names.at(i).Contains("ee_")) {
    	    templateDY_ee = (TH1D*) templateDY->Clone();
	    templateDY_ee->SetDirectory(0);
    	  }
	  else if (Vec_Files_root_names.at(i).Contains("ee_")) {
    	    templateDY_ee->Add(templateDY);
    	  }
	  if (!templateDY_mumu && Vec_Files_root_names.at(i).Contains("mumu_")) {
    	    templateDY_mumu = (TH1D*) templateDY->Clone();
	    templateDY_mumu->SetDirectory(0);
    	  }
	  else if (Vec_Files_root_names.at(i).Contains("mumu_")) {
    	    templateDY_mumu->Add(templateDY);
    	  }
	  //uncomment to perform fit for emu
	  //if (!templateDY_emu && Vec_Files_root_names.at(i).Contains("emu_")) {
    	  //  templateDY_emu = (TH1D*) templateDY->Clone();
	  //  templateDY_emu->SetDirectory(0);
    	  //}
	  //else if (Vec_Files_root_names.at(i).Contains("emu_")) {
    	  //  templateDY_emu->Add(templateDY);
	  //}
	  delete templateDY;
	} else {
	  TH1D *templateOtherMC;
	  if (!flagUseZpeak)
	    templateOtherMC = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString(TString("dyScaling_Allh1")).Append(nameAppendix));
	  else
	    templateOtherMC = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString(TString("dyScaling_Zh1").Append(nameAppendix)).Append("zWindow"));
    	  ApplyFlatWeights(templateOtherMC, LumiWeight);
	  if (!templateOtherMC_ee && Vec_Files_root_names.at(i).Contains("ee_")) {
    	    templateOtherMC_ee = (TH1D*) templateOtherMC->Clone();
	    templateOtherMC_ee->SetDirectory(0);
    	  }
	  else if (Vec_Files_root_names.at(i).Contains("ee_")) {
    	    templateOtherMC_ee->Add(templateOtherMC);
    	  }
	  if (!templateOtherMC_mumu && Vec_Files_root_names.at(i).Contains("mumu_")) {
    	    templateOtherMC_mumu = (TH1D*) templateOtherMC->Clone();
	    templateOtherMC_mumu->SetDirectory(0);
    	  }
	  else if (Vec_Files_root_names.at(i).Contains("mumu_")) {
    	    templateOtherMC_mumu->Add(templateOtherMC);
    	  }
	  //uncomment to perform fit for emu
	  //if (!templateOtherMC_emu && Vec_Files_root_names.at(i).Contains("emu_")) {
    	  //  templateOtherMC_emu = (TH1D*) templateOtherMC->Clone();
	  //  templateOtherMC_emu->SetDirectory(0);
    	  //}
	  //else if (Vec_Files_root_names.at(i).Contains("emu_")) {
    	  //  templateOtherMC_emu->Add(templateOtherMC);
	  //}
    	  delete templateOtherMC;
    	}
      }
    } //end loop over Vec_Files


    TObjArray *mc_ee = new TObjArray(2);
    mc_ee->Add(templateDY_ee);
    mc_ee->Add(templateOtherMC_ee);

    TObjArray *mc_mumu = new TObjArray(2);
    mc_mumu->Add(templateDY_mumu);
    mc_mumu->Add(templateOtherMC_mumu);

    //uncomment to perform fit for emu
    //TObjArray *mc_emu = new TObjArray(2);
    //mc_emu->Add(templateDY_emu);
    //mc_emu->Add(templateOtherMC_emu);

    // initialise fits
    TFractionFitter* fit_ee   = new TFractionFitter(inputData_ee, mc_ee, "V");
    TFractionFitter* fit_mumu = new TFractionFitter(inputData_mumu, mc_mumu, "V");
    //TFractionFitter* fit_emu  = new TFractionFitter(inputData_emu, mc_emu); //uncomment to perform fit for emu

    // constrain fraction 0 to be between 0 and 1
    fit_ee->Constrain(0,0.9,1.0);
    fit_mumu->Constrain(0,0.9,1.0);
    //fit_emu->Constrain(0,0.0,1.0); //uncomment to perform fit for emu

    // constrain fraction 1 to be between 0 and 1
    fit_ee->Constrain(1,0.01,0.07);
    fit_mumu->Constrain(1,0.01,0.07);
    //fit_emu->Constrain(1,0.0,1.0); //uncomment to perform fit for emu

    // perform the fit
    Int_t status_ee = fit_ee->Fit();
    Int_t status_mumu = fit_mumu->Fit();
    //Int_t status_emu = fit_emu->Fit(); //uncomment to perform fit for emu

    std::cout << "fit status of ee channel: " << status_ee << std::endl;
    std::cout << "fit status of mumu channel: " << status_mumu << std::endl;
    //std::cout << "fit status of emu channel: " << status_emu << std::endl; //uncomment to perform fit for emu

    // check on fit status - ee channel
    if (status_ee == 0) {

      TCanvas* c = new TCanvas("h_fit_status_ee", "h_fit_status_ee", 1000, 800);
      TH1F* plot_result_ee = (TH1F*) fit_ee->GetPlot();

      THStack *hs = new THStack("hs", "ee channel");
      gStyle->SetErrorX(0);

      inputData_ee->SetMarkerStyle(20);
      inputData_ee->SetMarkerColor(kBlack);
      inputData_ee->SetStats(false);
      hs->Add(inputData_ee, "P");

      plot_result_ee->SetLineColor(kRed);
      hs->Add(plot_result_ee);

      templateDY_ee->SetMarkerStyle(2);
      templateDY_ee->SetMarkerColor(kBlue);
      templateDY_ee->SetMarkerSize(3);
      hs->Add(templateDY_ee, "P");

      templateOtherMC_ee->SetMarkerStyle(2);
      templateOtherMC_ee->SetMarkerColor(kGreen);
      templateOtherMC_ee->SetMarkerSize(3);
      hs->Add(templateOtherMC_ee, "P");

      hs->Draw("nostack");
      hs->GetXaxis()->SetLimits(56,126);
      hs->GetXaxis()->SetTitle("m_{ll}");

      TLegend* legend = new TLegend(.1,.7,.3,.9);
      legend->SetFillColor(0);

      legend->AddEntry(plot_result_ee, "data predicted (after fit)", "L");
      legend->AddEntry(inputData_ee, "data", "P");
      legend->AddEntry(templateDY_ee, "Drell-Yan template", "P");
      legend->AddEntry(templateOtherMC_ee, "MC other template", "P");
      legend->Draw();

      c->Print(outputPlotsFolder + "fit_status_ee" + ".pdf");
      legend->Clear();
      c->Clear();
      delete plot_result_ee;
      delete hs;
    }
    // check on fit status - mumu channel
    if (status_mumu == 0) {

      TCanvas* c = new TCanvas("h_fit_status_mumu", "h_fit_status_mumu", 1000, 800);
      TH1F* plot_result_mumu = (TH1F*) fit_mumu->GetPlot();

      THStack *hs = new THStack("hs", "#mu#mu channel");
      gStyle->SetErrorX(0);

      inputData_mumu->SetMarkerStyle(20);
      inputData_mumu->SetMarkerColor(kBlack);
      inputData_mumu->SetStats(false);
      hs->Add(inputData_mumu, "P");

      plot_result_mumu->SetLineColor(kRed);
      hs->Add(plot_result_mumu);

      templateDY_mumu->SetMarkerStyle(2);
      templateDY_mumu->SetMarkerColor(kBlue);
      templateDY_mumu->SetMarkerSize(3);
      hs->Add(templateDY_mumu, "P");

      templateOtherMC_mumu->SetMarkerStyle(2);
      templateOtherMC_mumu->SetMarkerColor(kGreen);
      templateOtherMC_mumu->SetMarkerSize(3);
      hs->Add(templateOtherMC_mumu, "P");

      hs->Draw("nostack");
      hs->GetXaxis()->SetLimits(56,126);
      hs->GetXaxis()->SetTitle("m_{ll}");

      TLegend* legend = new TLegend(.1,.7,.3,.9);
      legend->SetFillColor(0);

      legend->AddEntry(plot_result_mumu, "data predicted (after fit)", "L");
      legend->AddEntry(inputData_mumu, "data", "P");
      legend->AddEntry(templateDY_mumu, "Drell-Yan template", "P");
      legend->AddEntry(templateOtherMC_mumu, "MC other template", "P");
      legend->Draw();

      c->Print(outputPlotsFolder + "fit_status_mumu" + ".pdf");
      legend->Clear();
      c->Clear();
      delete plot_result_mumu;
      delete hs;
    }

    double DY_frac_ee; double DY_frac_ee_error;
    double DY_frac_mumu; double DY_frac_mumu_error;
    //double DY_frac_emu; double DY_frac_emu_error; //uncomment to perform fit for emu

    fit_ee->GetResult(0, DY_frac_ee, DY_frac_ee_error);
    fit_mumu->GetResult(0, DY_frac_mumu, DY_frac_mumu_error);
    //fit_emu->GetResult(0, DY_frac_emu, DY_frac_emu_error); //uncomment to perform fit for emu

    //ee
    int dataBins_ee = inputData_ee->GetNbinsX();
    double data_ee = inputData_ee->Integral(0,dataBins_ee+1);
    double data_ee_err = std::sqrt(data_ee);

    double DY_data_fit_ee = DY_frac_ee*data_ee;
    double DY_data_fit_ee_err = DY_data_fit_ee*std::sqrt(TMath::Power((DY_frac_ee_error/DY_frac_ee),2)+TMath::Power((data_ee_err/data_ee),2));

    int DYBins_ee = templateDY_ee->GetNbinsX();
    double DY_MC_ee = templateDY_ee->Integral(0,DYBins_ee+1);
    double DY_MC_ee_err = std::sqrt(DY_MC_ee);


    //mumu
    int dataBins_mumu = inputData_mumu->GetNbinsX();
    double data_mumu = inputData_mumu->Integral(0,dataBins_mumu+1);
    double data_mumu_err = std::sqrt(data_mumu);

    double DY_data_fit_mumu = DY_frac_mumu*data_mumu;
    double DY_data_fit_mumu_err = DY_data_fit_mumu*std::sqrt(TMath::Power((DY_frac_mumu_error/DY_frac_mumu),2)+TMath::Power((data_mumu_err/data_mumu),2));

    int DYBins_mumu = templateDY_mumu->GetNbinsX();
    double DY_MC_mumu = templateDY_mumu->Integral(0,DYBins_mumu+1);
    double DY_MC_mumu_err = std::sqrt(DY_MC_mumu);


    //emu
    int dataBins_emu = inputData_emu->GetNbinsX();
    double data_emu = inputData_emu->Integral(0,dataBins_emu+1);
    double data_emu_err = std::sqrt(data_emu);

    //uncomment to perform fit for emu
    //double DY_data_fit_emu = DY_frac_emu*data_emu;
    //double DY_data_fit_emu_err = DY_data_fit_emu*std::sqrt(TMath::Power((DY_frac_emu_error/DY_frac_emu),2)+TMath::Power((data_emu_err/data_emu),2));

    //int DYBins_emu = templateDY_emu->GetNbinsX();
    //double DY_MC_emu = templateDY_emu->Integral(0,DYBins_emu+1);
    //double DY_MC_emu_err = std::sqrt(DY_MC_emu);

    //scale factor calculations
    double DYSFEE = DY_data_fit_ee/DY_MC_ee;
    double DYSFEE_error = DYSFEE*std::sqrt(TMath::Power((DY_data_fit_ee_err/DY_data_fit_ee),2)+TMath::Power((DY_MC_ee_err/DY_MC_ee),2));

    double DYSFMuMu = DY_data_fit_mumu/DY_MC_mumu;
    double DYSFMuMu_error = DYSFMuMu*std::sqrt(TMath::Power((DY_data_fit_mumu_err/DY_data_fit_mumu),2)+TMath::Power((DY_MC_mumu_err/DY_MC_mumu),2));

    //uncomment to perform fit for emu
    //double DYSFEMu = DY_data_fit_emu/DY_MC_emu;
    //double DYSFEMu_error = DYSFEMu*std::sqrt(TMath::Power((DY_data_fit_emu_err/DY_data_fit_emu),2)+TMath::Power((DY_MC_emu_err/DY_MC_emu),2));

    //comment to perform fit for emu
    double DYSFEMu = std::sqrt(DYSFEE*DYSFMuMu);
    double DYSFEMu_error = 2*DYSFEMu*std::sqrt(TMath::Power((DYSFEE_error/DYSFEE),2)+TMath::Power((DYSFMuMu_error/DYSFMuMu),2));

    double DYSFDilept = (DYSFEE * data_ee + DYSFMuMu * data_mumu + DYSFEMu * data_emu) / (data_ee + data_mumu + data_emu);
    double DYSFDilept_error = DYSFDilept*std::sqrt(TMath::Power((DYSFEE_error/DYSFEE),2)+2*TMath::Power((data_ee_err/data_ee),2)+
						   TMath::Power((DYSFMuMu_error/DYSFMuMu),2)+2*TMath::Power((data_mumu_err/data_mumu),2)+
						   TMath::Power((DYSFEMu_error/DYSFEMu),2)+2*TMath::Power((data_emu_err/data_emu),2));

    std::cout << std::endl;
    std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;
    std::cout << "Calculation of DY Scale Factors for '" << name << "'  at selection step " << nameAppendix << std::endl;
    std::cout << "DYSFEE:                 " << DYSFEE << std::endl;
    std::cout << "DYSFEE error:           " << DYSFEE_error << std::endl;
    std::cout << "DYSFMuMu:               " << DYSFMuMu << std::endl;
    std::cout << "DYSFMuMu error:         " << DYSFMuMu_error << std::endl;
    std::cout << "DYSFEMu:                " << DYSFEMu << std::endl;
    std::cout << "DYSFEMu error:          " << DYSFEMu_error << std::endl;
    std::cout << "DYSFDilept:             " << DYSFDilept << std::endl;
    std::cout << "DYSFDilept error:       " << DYSFDilept_error << std::endl;
    std::cout << "      " << latexLabel << " &  " << round(DYSFEE * 1000)/1000 << " $\\pm$ " << round(DYSFEE_error * 1000)/1000 << " & " << round(DYSFMuMu * 1000)/1000 << " $\\pm$ " << round(DYSFMuMu_error * 1000)/1000 << " & " << round(DYSFEMu * 1000)/1000 << " $\\pm$ " << round(DYSFEMu_error * 1000)/1000 << " & " << round(DYSFDilept * 1000)/1000 << " $\\pm$ " << round(DYSFDilept_error * 1000)/1000 << " \\\\" << std::endl;
    std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;
    std::cout << std::endl;

    if (DYSFEE < 0.2 || DYSFMuMu < 0.2) {
        std::cout << "The DY SF is too low (below 0.2). Something is probably wrong.\n";
        std::exit(1);
    }


    DYScale.at(0) = DYSFEE;
    DYScale.at(1) = DYSFMuMu;
    DYScale.at(2) = DYSFEMu;
    DYScale.at(3) = DYSFDilept;


    delete inputData_ee;
    delete templateDY_ee;
    delete templateOtherMC_ee;

    delete inputData_mumu;
    delete templateDY_mumu;
    delete templateOtherMC_mumu;

    //uncomment to perform fit for emu
    //delete inputData_emu;
    //delete templateDY_emu;
    //delete templateOtherMC_emu;

    std::cout << "End DYSCALE FACTOR calculation\n" << std::endl;
}


/////////////////////////////////////////////////////////////
//                                                         //
// General settings and methods for TUnfold implementation //
//                                                         //
/////////////////////////////////////////////////////////////

 std::vector<double> PlotterConfigurationHelper::load_ttmd_DY_SFs(TString year){
    // "DYSFEE", "DYSFMuMu", "DYSFEMu", "diLepton"
    // scale factors are computed from fit to Z-exclusion region
    std::vector<double> SFs;
    if (year == "2016preVFP"){
      SFs = {1.0, 1.0, 1.0, 1.0};
    }
    else if (year == "2016postVFP"){
      SFs = {1.0, 1.0, 1.0, 1.0};
    }
    else if (year == "2016"){
      SFs = {1.0, 1.0, 1.0, 1.0};
    }
    else if (year == "2017"){
      SFs = {1.0, 1.0, 1.0, 1.0};
    }
    else if (year == "2018"){
      SFs = {1.0, 1.0, 1.0, 1.0};
    } else {
        std::cerr << "ERROR in static std::vector<double> load_ttmd_DY_SFs(TString year) --->>> year not implemented!" << std::endl;
        exit(1);
    }

    return SFs;
}

double PlotterConfigurationHelper::get_DY_SF(TString year, TString channel){
    std::vector<double> ttmd_DY_SFs =  load_ttmd_DY_SFs(year);
    double SF = -99.99;
    if (channel == "ee"){
        SF = ttmd_DY_SFs.at(0);
    }
    else if (channel == "mumu"){
        SF = ttmd_DY_SFs.at(1);
    }
    else if (channel == "emu"){
        SF = ttmd_DY_SFs.at(2);
    }
    else if (channel == "ll"){
        SF = ttmd_DY_SFs.at(3);
    }
    else{
        std::cerr << "ERROR in PlotterConfigurationHelper::get_DY_SF --->>> channel not implemented!" << std::endl;
        exit(1);
    }

    return SF;

}


void PlotterConfigurationHelper::DefaultSetup()
{
    // flush
    std::cout.setf(std::ios::unitbuf);
    setbuf(stdout, NULL);

    gTTBARDebug = 1;
    utils::ttmd::SetOZStyle();
    gStyle->SetPaintTextFormat("0.2f");
    TH1::SetDefaultSumw2();
}

double PlotterConfigurationHelper::GetTrueLevelRenormalisationWeight(const TString& WhichSample, const TString& varName){
    double renormWeight = -999.99;

    TString syst_weight_hist_branch_name = "normalization_histograms";
    TString weight_sum_hist_base_name = "trueLevelWeightSum_hist";
    TString no_weight_sum_hist_base_name = "trueLevelNoRenormalisationWeightSum_hist";
    TString weight_sum_hist_name = "";
    TString no_weight_sum_hist_name = "";

    //set histogram names for the specific systematic
    if (varName.BeginsWith("PU_") || varName.BeginsWith("MESCALE_") ||
        varName.BeginsWith("MEFACSCALE_") || varName.BeginsWith("MERENSCALE_")||
        varName.BeginsWith("BFRAG_") || varName.BeginsWith("BSEMILEP_") || varName.BeginsWith("PDF_") || varName.BeginsWith("PSSCALE_WEIGHT_")){
        
        TString varNameHist = varName;
        if (varName.BeginsWith("PSSCALE_WEIGHT_")) {
            varNameHist = varNameHist.ReplaceAll("_4_UP", "_UP_4");
            varNameHist = varNameHist.ReplaceAll("_5_UP", "_UP_5");
            varNameHist = varNameHist.ReplaceAll("_4_DOWN", "_DOWN_4");
            varNameHist = varNameHist.ReplaceAll("_5_DOWN", "_DOWN_5");
        }
        weight_sum_hist_name = syst_weight_hist_branch_name + "/" + weight_sum_hist_base_name + "_" + varNameHist;
    }
    else{
        weight_sum_hist_name = weight_sum_hist_base_name;
    }

    if (varName.BeginsWith("PU_")) no_weight_sum_hist_name = syst_weight_hist_branch_name + "/" + no_weight_sum_hist_base_name + varName;
    else no_weight_sum_hist_name = no_weight_sum_hist_base_name;

    // std::cout << varName << ": " << weight_sum_hist_name << "  " << no_weight_sum_hist_name << std::endl;

    const TH1D *h_trueLevelWeightSum_hist = fileReader->Get<TH1D>(WhichSample, weight_sum_hist_name);
    const TH1D *h_trueLevelNoRenormalisationWeightSum_hist = fileReader->Get<TH1D>(WhichSample, no_weight_sum_hist_name);
    const double trueLevelWeightSum_ = h_trueLevelWeightSum_hist->GetBinContent(1);
    const double trueLevelNoRenormalisationWeightSum_ = h_trueLevelNoRenormalisationWeightSum_hist->GetBinContent(1);
    renormWeight = trueLevelWeightSum_ != 0. ? trueLevelNoRenormalisationWeightSum_/trueLevelWeightSum_ : -999.;
    renormWeight = trueLevelNoRenormalisationWeightSum_/trueLevelWeightSum_;

    if(renormWeight <= 0.) {
        renormWeight=1.;
        std::cout<<"    *** WARNING: True lvl. renormalisation weight < 0!! Value has been set to 1.0!" << std::endl;
    }
    else std::cout<<"True lvl. renormalisation weight: "<<renormWeight<<std::endl;

    return renormWeight;

}

// ======Porting of ttmd_settings========
void PlotterConfigurationHelper::fillPlainTreesFileListSampleVectors(TString in_year, std::vector<TString>& pattern_vector, std::vector<TString>& output_vector, TString syst, TString channel){

    output_vector.clear();
    TString fileListName = "";
    TString baseFileList_plainTrees = "FileLists_"+in_year+"/";
    if (syst!="" && channel!=""){
        fileListName = this->fileList_prefix_plainTrees + syst + "_" + channel + this->fileList_format_plainTrees;
    }
    else {
        //Look for a Nominal FileList
        bool file_found = false;
        std::vector<TString> ch_to_check{"mumu","ee","emu"};
        for (TString i_ch: ch_to_check){
            fileListName = this->fileList_prefix_plainTrees + "Nominal_" + i_ch + this->fileList_format_plainTrees;
            std::ifstream infile(baseFileList_plainTrees + fileListName);
            if (infile.good()) {
                file_found = true;
                break;
            }
        }
        if (!file_found){
            std::cerr << "ERROR in PlotterConfigurationHelper::fillPlainTreesFileListSampleVectors -->> No Nominal fileList found for any channel!!" << std::endl;
            exit(1);
        }
    }
    TString fileListPath = baseFileList_plainTrees + fileListName;

    //Printing some info
    std::cout << std::endl;
    std::cout << "Reading patterns: ";
    for (TString pattern: pattern_vector) std::cout << pattern << " ";
    std::cout << std::endl << " from file: " + fileListPath << std::endl;

    std::ifstream FileList(fileListPath);
    if (!FileList.is_open()) {
        std::cerr << "ERROR in PlotterConfigurationHelper::FillPlainTreesFileListSampleVectors -->> " << fileListPath << " failed to open!" << std::endl;
        exit(1);
    }
    TString filename;
    while (!FileList.eof()) {
        FileList >> filename;
        filename = ExtractFileNameFromPath(filename,true);
        for (TString pattern: pattern_vector){
            if (filename.Contains(pattern)){
                std::cout << "      " + pattern + " <-- file: " + filename << std::endl;
                output_vector.push_back(filename);
            }

        }
    }
    std::cout << std::endl;

}

const TString PlotterConfigurationHelper::ExtractFileNameFromPath(TString path, bool exclude_channel_label){
    TString filename = "";
    std::vector<TString> channels = {"mumu_", "ee_", "emu_"};
    TObjArray* path_array  = path.Tokenize("/");
    if (!path_array->IsEmpty()){
        filename = ((TObjString*)path_array->Last())->GetString();
        if (exclude_channel_label) {
            for (TString ch: channels){
                if (filename.Contains(ch)){
                    filename.ReplaceAll(ch, "");
                }
            }
         }
    }

    delete path_array;
    return filename;
}

void PlotterConfigurationHelper::fillPlainTreesSampleIdVector(std::vector<TString>& fileListVector, std::vector<std::pair<TString, int>>& pattern_vector, std::vector<std::pair<TString, int>>& output_vector){
    output_vector.clear();
    for (TString file: fileListVector){
        for (std::pair<TString, int> pattern_id: pattern_vector){
            if (file.Contains(pattern_id.first)) output_vector.push_back({file, pattern_id.second});
        }
    }
}

const std::vector<TString>& PlotterConfigurationHelper::GetJECSrc(const bool fill_Up_Down)
{
    //vect_all, year, useAdditionalSystematics, fill_Up_Down, include_all, include_Nominal, include_noUpDown, only_noUpDown, sys_type
    // Fill list of JES systematics excluding "Nominal" for TUnfold implementation
    static std::vector<TString> fill_vect;
    fill_vect.clear();
    if (gSysJES == 1) fillVectorOfValidSystematics(fill_vect, year_, yearCorr_, true, false, fill_Up_Down, false, false, false, false, "JESSys");
    if (gSysJES == 2) fillVectorOfValidSystematics(fill_vect, year_, yearCorr_, true, false, fill_Up_Down, false, false, false, false, "JESSysSplit");

    static std::vector<TString> src;
    src.clear();
    src.insert(src.end(), fill_vect.begin(), fill_vect.end());

    return src;
}


const std::vector<TString>& PlotterConfigurationHelper::GetExpSys(const bool fill_Up_Down)
{

    //vect_all, year, useAdditionalSystematics, fill_Up_Down, include_all, include_Nominal, include_noUpDown, only_noUpDown, sys_type
    // Fill list of experimental systematics excluding "Nominal" for TUnfold implementation
    static std::vector<TString> fill_vect;
    fill_vect.clear();
    fillVectorOfValidSystematics(fill_vect, year_, yearCorr_, true, false, fill_Up_Down, false, false, false, false, "ExpSys");

    static std::vector<TString> syst;
    syst.clear();
    if(syst.size() == 0)
    {
      syst.insert(syst.end(), fill_vect.begin(), fill_vect.end());//

      for(auto& src : GetJECSrc(fill_Up_Down))
        if(!src.Contains("Total"))
          syst.push_back(src);
    }

    return syst;

}

const std::vector<TString>& PlotterConfigurationHelper::GetModPDFEig()
{
    static std::vector<TString> syst;
    if(syst.size() == 0)
    {
      for(int i = 1; i <= 28; i++)
        syst.push_back(TString::Format("PDF_CT14_eig%d", i));
    }
    return syst;
}

bool PlotterConfigurationHelper::IsModPDFEig(const TString& syst)
{
    auto vSyst = GetModPDFEig();
    return (std::find(vSyst.begin(), vSyst.end(), syst) != vSyst.end());
}

const std::vector<TString>& PlotterConfigurationHelper::GetModWSys(const bool fill_Up_Down)
{
    //vect_all, year, useAdditionalSystematics, fill_Up_Down, include_all, include_Nominal, include_noUpDown, only_noUpDown, sys_type
    // Fill list of ModWSys systematics excluding "Nominal" for TUnfold implementation
    static std::vector<TString> fill_vect;
    fill_vect.clear();
    fillVectorOfValidSystematics(fill_vect, year_, yearCorr_, true, false, fill_Up_Down, false, false, true, false, "ModWSys");

    static std::vector<TString> syst;
    syst.clear();
    if(syst.size() == 0)
    {
      syst.insert(syst.end(), fill_vect.begin(), fill_vect.end());//
      if(gSysPDF)
      {
        auto v2 = GetModPDFEig();
        syst.insert(syst.end(), v2.begin(), v2.end());
        // and CT14 as variation
        syst.push_back("PDF_CT14_AS");
      }
    }

    return syst;
}

const std::vector<TString>& PlotterConfigurationHelper::GetModISys(const bool fill_Up_Down)
{
    //vect_all, year, useAdditionalSystematics, fill_Up_Down, include_all, include_Nominal, include_noUpDown, only_noUpDown, sys_type
    // Fill list of ModISys systematics excluding "Nominal" for TUnfold implementation
    static std::vector<TString> fill_vect;
    fill_vect.clear();
    fillVectorOfValidSystematics(fill_vect, year_, yearCorr_, true, false, fill_Up_Down, false, false, true, false, "ModISys");

    static std::vector<TString> syst;
    syst.clear();
    if(syst.size() == 0)
    {
      syst.insert(syst.end(), fill_vect.begin(), fill_vect.end());//
    }

    return syst;
}

const std::vector<TString>& PlotterConfigurationHelper::GetAltTheorSys(const bool fill_Up_Down)
{

    //vect_all, year, useAdditionalSystematics, fill_Up_Down, include_all, include_Nominal, include_noUpDown, only_noUpDown, sys_type
    // Fill list of alternative theory systematics excluding "Nominal" for TUnfold implementation
    static std::vector<TString> fill_vect;
    fill_vect.clear();
    fillVectorOfValidSystematics(fill_vect, year_, yearCorr_, true, false, fill_Up_Down, false, false, true, false, "AltTheorySys");

    static std::vector<TString> syst;
    syst.clear();
    if(syst.size() == 0)
    {
      syst.insert(syst.end(), fill_vect.begin(), fill_vect.end());//
    }

    return syst;
}

const std::vector<TString>& PlotterConfigurationHelper::GetAllSys(const bool fill_Up_Down)
{
    static std::vector<TString> result;
    result.clear();

    if(result.size() == 0)
    {
      auto v1 = GetExpSys(fill_Up_Down);
      result.insert(result.end(), v1.begin(), v1.end());
      auto v2 = GetModWSys(fill_Up_Down);
      result.insert(result.end(), v2.begin(), v2.end());
      auto v3 = GetModISys(fill_Up_Down);
      result.insert(result.end(), v3.begin(), v3.end());
    }

    return result;
}

const std::vector<TString>& PlotterConfigurationHelper::GetAllVar(const bool fill_Up_Down)
{
    static std::vector<TString> result;
    result.clear();

    if(result.size() == 0)
    {
      auto syst = GetAllSys(fill_Up_Down);
      result.insert(result.end(), syst.begin(), syst.end());
      result.insert(result.begin(), "Nominal");
    }

    return result;
}

bool PlotterConfigurationHelper::IsExpSys(const TString& syst)
{
    auto vSyst = GetExpSys(false);
    return (std::find(vSyst.begin(), vSyst.end(), syst) != vSyst.end());
}


bool PlotterConfigurationHelper::IsModWSys(const TString& syst)
{
    if (syst == "TOT_BFRAG" || syst == "TOT_SCALE") return true;
    else{
        auto vSyst = GetModWSys(false);
        return (std::find(vSyst.begin(), vSyst.end(), syst) != vSyst.end());
    }
}

bool PlotterConfigurationHelper::IsModISys(const TString& syst)
{
    if (syst == "TOT_COLORREC" || syst == "CR") return true;
    else{
        auto vSyst = GetModISys(false);
        return (std::find(vSyst.begin(), vSyst.end(), syst) != vSyst.end());
    }
}

bool PlotterConfigurationHelper::IsAltTheorSys(const TString& syst)
{
    auto vSyst = GetAltTheorSys(false);
    return (std::find(vSyst.begin(), vSyst.end(), syst) != vSyst.end());
}

bool PlotterConfigurationHelper::IsNoUpDown(const TString& syst)
{
    if (syst == "Nominal") return true;
    else return (std::find(_vectorOfSystematics_noUpDown.begin(), _vectorOfSystematics_noUpDown.end(), syst) != _vectorOfSystematics_noUpDown.end());
}

bool PlotterConfigurationHelper::IsEnvelope(const TString& sys){
    bool result = false;
    auto it = this->_envelopesMap.find(sys);
    if (it != this->_envelopesMap.end()) result = true;

    return result;
}

bool PlotterConfigurationHelper::IsYearToYearCorr(const TString& year, const TString& yearCorr, const TString& sys){
    bool result = false;

    TString temp_sys = sys;
    RemoveYearCorr(year, yearCorr, temp_sys);

    const char* s = temp_sys.Data();
    const int n = temp_sys.Length();
    if (temp_sys.Contains("_UP"))
      temp_sys = TString(s, n - 3);
    else if (temp_sys.Contains("_DOWN"))
      temp_sys = TString(s, n - 5);

    auto it = YearToYearCorrMap.find(temp_sys);
    if (it != YearToYearCorrMap.end()) result = true;

    return result;
}

const std::vector<TString> PlotterConfigurationHelper::GetYearCorrList(const TString& year, const TString& yearCorr, const TString& sys, const bool isYearToYearCorr){

  std::vector<TString> yearCorrList;
if (year == "fullRun2UL" && isYearToYearCorr) {

    TString labelYearCorr;
    TString sysTmp16;
    TString sysTmp17;
    TString sysTmp18;
    TString sysTmp1617;
    TString sysTmp1618;
    TString sysTmp1718;
    TString sysTmp161718;

    if (yearCorr == "2016") {
      sysTmp16 = sys;
      labelYearCorr = "16";
      AppendYearCorr(year, labelYearCorr, sysTmp16, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp16) != 0.)
	yearCorrList.push_back(sysTmp16);

      sysTmp1617 = sys;
      labelYearCorr = "1617";
      AppendYearCorr(year, labelYearCorr, sysTmp1617, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp1617) != 0.)
	yearCorrList.push_back(sysTmp1617);

      sysTmp1618 = sys;
      labelYearCorr = "1618";
      AppendYearCorr(year, labelYearCorr, sysTmp1618, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp1618) != 0.)
	yearCorrList.push_back(sysTmp1618);

      //in case two years are fully correlated but partially correlated with the third year
      sysTmp161718 = sys;
      labelYearCorr = "161718";
      AppendYearCorr(year, labelYearCorr, sysTmp161718, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp161718) != 0.)
	yearCorrList.push_back(sysTmp161718);
    }
    else if (yearCorr == "2017") {
      sysTmp17 = sys;
      labelYearCorr = "17";
      AppendYearCorr(year, labelYearCorr, sysTmp17, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp17) != 0.)
	yearCorrList.push_back(sysTmp17);

      sysTmp1617 = sys;
      labelYearCorr = "1617";
      AppendYearCorr(year, labelYearCorr, sysTmp1617, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp1617) != 0.)
	yearCorrList.push_back(sysTmp1617);

      sysTmp1718 = sys;
      labelYearCorr = "1718";
      AppendYearCorr(year, labelYearCorr, sysTmp1718, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp1718) != 0.)
	yearCorrList.push_back(sysTmp1718);

      //in case two years are fully correlated but partially correlated with the third year
      sysTmp161718 = sys;
      labelYearCorr = "161718";
      AppendYearCorr(year, labelYearCorr, sysTmp161718, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp161718) != 0.)
	yearCorrList.push_back(sysTmp161718);
    }
    else if (yearCorr == "2018") {
      sysTmp18 = sys;
      labelYearCorr = "18";
      AppendYearCorr(year, labelYearCorr, sysTmp18, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp18) != 0.)
	yearCorrList.push_back(sysTmp18);

      sysTmp1618 = sys;
      labelYearCorr = "1618";
      AppendYearCorr(year, labelYearCorr, sysTmp1618, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp1618) != 0.)
	yearCorrList.push_back(sysTmp1618);

      sysTmp1718 = sys;
      labelYearCorr = "1718";
      AppendYearCorr(year, labelYearCorr, sysTmp1718, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp1718) != 0.)
	yearCorrList.push_back(sysTmp1718);

      //in case two years are fully correlated but partially correlated with the third year
      sysTmp161718 = sys;
      labelYearCorr = "161718";
      AppendYearCorr(year, labelYearCorr, sysTmp161718, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp161718) != 0.)
	yearCorrList.push_back(sysTmp161718);
    }
    else if (yearCorr == "all") {

      //2016
      sysTmp16 = sys;
      labelYearCorr = "16";
      AppendYearCorr(year, labelYearCorr, sysTmp16, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp16) != 0.)
	yearCorrList.push_back(sysTmp16);

      //2017
      sysTmp17 = sys;
      labelYearCorr = "17";
      AppendYearCorr(year, labelYearCorr, sysTmp17, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp17) != 0.)
	yearCorrList.push_back(sysTmp17);

      //2018
      sysTmp18 = sys;
      labelYearCorr = "18";
      AppendYearCorr(year, labelYearCorr, sysTmp18, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp18) != 0.)
	yearCorrList.push_back(sysTmp18);

      //1617
      sysTmp1617 = sys;
      labelYearCorr = "1617";
      AppendYearCorr(year, labelYearCorr, sysTmp1617, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp1617) != 0.)
	yearCorrList.push_back(sysTmp1617);

      //1618
      sysTmp1618 = sys;
      labelYearCorr = "1618";
      AppendYearCorr(year, labelYearCorr, sysTmp1618, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp1618) != 0.)
	yearCorrList.push_back(sysTmp1618);

      //1718
      sysTmp1718 = sys;
      labelYearCorr = "1718";
      AppendYearCorr(year, labelYearCorr, sysTmp1718, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp1718) != 0.)
	yearCorrList.push_back(sysTmp1718);

      //in case two years are fully correlated but partially correlated with the third year
      sysTmp161718 = sys;
      labelYearCorr = "161718";
      AppendYearCorr(year, labelYearCorr, sysTmp161718, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp161718) != 0.)
	yearCorrList.push_back(sysTmp161718);
    }
    else if (yearCorr == "only16") {
      sysTmp16 = sys;
      labelYearCorr = "16";
      AppendYearCorr(year, labelYearCorr, sysTmp16, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp16) != 0.) yearCorrList.push_back(sysTmp16);
    }
    else if (yearCorr == "only17") {
      sysTmp17 = sys;
      labelYearCorr = "17";
      AppendYearCorr(year, labelYearCorr, sysTmp17, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp17) != 0.) yearCorrList.push_back(sysTmp17);
    }
    else if (yearCorr == "only18") {
      sysTmp18 = sys;
      labelYearCorr = "18";
      AppendYearCorr(year, labelYearCorr, sysTmp18, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp18) != 0.) yearCorrList.push_back(sysTmp18);
    }
    else if (yearCorr == "only161718") {
      sysTmp161718 = sys;
      labelYearCorr = "161718";
      AppendYearCorr(year, labelYearCorr, sysTmp161718, isYearToYearCorr);
      if (GetYearToYearCorrRescale(year, yearCorr, sysTmp161718) != 0.) yearCorrList.push_back(sysTmp161718);
    }
  }

  return yearCorrList;
}

void PlotterConfigurationHelper::AppendYearCorr(const TString& year, const TString& labelYearCorr, TString& sys, const bool isYearToYearCorr, const bool noUpDown){
  if (year == "fullRun2UL" && labelYearCorr != "" && isYearToYearCorr) {
    TString sysNoUpDown = sys;

    const char* s = sys.Data();
    const int n = sys.Length();
    if (sys.Contains("_UP"))
      sysNoUpDown = TString(s, n - 3);
    else if (sys.Contains("_DOWN"))
      sysNoUpDown = TString(s, n - 5);

    if (sys.Contains("_UP") && !noUpDown)
      sys = sysNoUpDown + "_" + labelYearCorr + "_UP";
    else if (sys.Contains("_DOWN") && !noUpDown)
      sys = sysNoUpDown + "_" + labelYearCorr + "_DOWN";
    else
      sys = sysNoUpDown + "_" + labelYearCorr;
  }
}

void PlotterConfigurationHelper::RemoveYearCorr(const TString& year, const TString& yearCorr, TString& sys, const bool noUpDown){
  if (year == "fullRun2UL" && yearCorr != "") {
    TString sysNoUpDown = sys;

    const char* s = sys.Data();
    const int n = sys.Length();
    if (sys.Contains("_UP"))
      sysNoUpDown = TString(s, n - 3);
    else if (sys.Contains("_DOWN"))
      sysNoUpDown = TString(s, n - 5);

    TString sysNoYear = sysNoUpDown;
    const char* s_NoUpDown = sysNoUpDown.Data();
    const int n_NoUpDown = sysNoUpDown.Length();

    if(TString(s_NoUpDown + n_NoUpDown - 6, 6) == "161718") {
      sysNoYear = TString(s_NoUpDown, n_NoUpDown - 7);
    }
    else if(TString(s_NoUpDown + n_NoUpDown - 4, 4) == "1617" || TString(s_NoUpDown + n_NoUpDown - 4, 4) == "1618" || TString(s_NoUpDown + n_NoUpDown - 4, 4) == "1718") {
      sysNoYear = TString(s_NoUpDown, n_NoUpDown - 5);
    }
    else if(TString(s_NoUpDown + n_NoUpDown - 2, 2) == "16" || TString(s_NoUpDown + n_NoUpDown - 2, 2) == "17" || TString(s_NoUpDown + n_NoUpDown - 2, 2) == "18") {
      sysNoYear = TString(s_NoUpDown, n_NoUpDown - 3);
    }


    if (sys.Contains("_UP") && !noUpDown)
      sys = sysNoYear + "_UP";
    else if (sys.Contains("_DOWN") && !noUpDown)
      sys = sysNoYear + "_DOWN";
    else
      sys = sysNoYear;
  }
}

const TString PlotterConfigurationHelper::GetLabelYearCorr(const TString& var){

  if(var.Contains("161718"))
    return "161718";
  else if(var.Contains("1617"))
    return "1617";
  else if(var.Contains("1618"))
    return "1618";
  else if(var.Contains("1718"))
    return "1718";
  else if(var.Contains("16"))
    return "16";
  else if(var.Contains("17"))
    return "17";
  else if(var.Contains("18"))
    return "18";
  else
    return "";
}

const std::vector<TString> PlotterConfigurationHelper::GetEnvelopeVars(const TString& sys, bool get_UD_vars){
    std::vector<TString> result = {};
    auto it = this->_envelopesMap.find(sys);
    if (it != this->_envelopesMap.end()) result = it->second;

    std::vector<TString> out_vector={};
    if (get_UD_vars) out_vector = result;
    else {
      for (auto var_syst: result){
        if (var_syst.Contains("_UP")) var_syst.ReplaceAll("_UP", "");
        else if (var_syst.Contains("_DOWN")) var_syst.ReplaceAll("_DOWN", "");
        
        if (std::find(out_vector.begin(), out_vector.end(), var_syst) == out_vector.end()) out_vector.push_back(var_syst);
      }
    }

    return out_vector;
}

const std::pair<TString,TString> PlotterConfigurationHelper::GetTotEnvelopeVars(const TString& sys){
    std::pair<TString,TString> up_down;
    auto it = this->_envelopesMap.find(sys);
    if (it != this->_envelopesMap.end()) {
        up_down.first = sys + "_UP";
        up_down.second = sys + "_DOWN";
    }
    else{
        std::cerr << "ERROR in PlotterConfigurationHelper::GetTotEnvelopeVars(): This in not an envelope name <-- " << sys << std::endl;
    }

    return up_down;
}

const std::vector<TString> PlotterConfigurationHelper::GetEnvelopesList(const bool fill_tot_up_down_variations){
    std::vector<TString> result = {};

    for (auto envelope_info: this->_envelopesMap){
        if (!fill_tot_up_down_variations)  result.push_back(envelope_info.first);
        else{
            result.push_back(envelope_info.first + "_UP");
            result.push_back(envelope_info.first + "_DOWN");
        }
    }
    return result;
}

const TString PlotterConfigurationHelper::GetVarEnvelope(const TString& var){
    TString result = "";
    for (auto i_envelope: this->_envelopesMap){
        for (auto i_var: i_envelope.second){
            if (i_var == var) {
                result = i_envelope.first;
                break;
            }
        }
    }
    return result;
}

 bool PlotterConfigurationHelper::IsVarInEnvelope(const TString& var, const TString& envelope){
    if (!IsEnvelope(envelope)) std::cerr << "ERROR in PlotterConfigurationHelper::GetTotEnvelopeVars(): This in not an envelope name <-- " << envelope << std::endl;

    std::vector<TString> vars_list = GetEnvelopeVars(envelope);

    return !(std::find(begin(vars_list), end(vars_list), var) == end(vars_list));

}

const std::vector<TString> PlotterConfigurationHelper::FoldSystematicVarsListIntoEnvelops(const std::vector<TString>& list_of_vars, const bool list_with_variations_included){
    TString temp_envelope;
    std::vector<TString> result_systs;
    result_systs.clear();

    for (auto i_var: list_of_vars){
        if (list_with_variations_included) {
            temp_envelope = GetVarEnvelope(i_var);
        }else{
            std::vector<TString> list_of_vars_to_check = {i_var};
            list_of_vars_to_check.push_back(i_var + "_UP");
            list_of_vars_to_check.push_back(i_var + "_DOWN");
            for (auto var_to_check: list_of_vars_to_check) {
                temp_envelope = GetVarEnvelope(var_to_check);
                if (temp_envelope != "") break;
            }
        }

        if (temp_envelope == "" && std::find(begin(result_systs), end(result_systs), i_var) == end(result_systs)) result_systs.push_back(i_var);
        else if (std::find(begin(result_systs), end(result_systs), temp_envelope) == end(result_systs)){
            result_systs.push_back(temp_envelope);
        }
    }

    return result_systs;

}

const std::vector<TString> PlotterConfigurationHelper::GetSysVars(const TString& sys, bool getEnvelopesTotVar)
{
  std::vector<TString> result;

  //check for possible envelope and get all vars if it is the case
  if (this->IsEnvelope(sys) && !getEnvelopesTotVar) return GetEnvelopeVars(sys);

  // special cases
  if(sys == "MASS_CONSTMC")
  {
    result.push_back("MASS_CONSTMC");
  }
  else if(sys == "REWTOPPT")
  {
    result.push_back("REWTOPPT");
  }
  // else if(sys.BeginsWith("PDF_CT14_eig"))
  // {
  //   int n = TString(sys.Data() + 12).Atoi();
  //   result.push_back(TString::Format("PDF_%d", 101 + 1 + (n - 1) * 2));
  //   result.push_back(TString::Format("PDF_%d", 101 + 1 + (n - 1) * 2 + 1));
  // }
  // else if(sys == "PDF_CT14_AS")
  // {
  //   result.push_back("PDF_158");
  //   result.push_back("PDF_159");
  // }
  // else if(sys.BeginsWith("PDF_NNPDF31_eig"))
  // {
  //   int n = TString(sys.Data() + 15).Atoi();
  //   result.push_back(TString::Format("PDF_%d", 2*n - 1));
  //   result.push_back(TString::Format("PDF_%d", 2*n));
  // }

  else if(sys == "POWHEGV2HERWIG")
    result.push_back("POWHEGV2HERWIG");
  else if(sys == "AMCATNLOFXFX")
    result.push_back("AMCATNLOFXFX");

  else if(this->IsNoUpDown(sys)) result.push_back(sys);

  // otherwise up-down
  else
  {
    result.push_back(sys + "_UP");
    result.push_back(sys + "_DOWN");
  }

  return result;
}

const TString PlotterConfigurationHelper::VarToSys(const TString& var, int& sign)
{
    sign = 0;
    if(var == "Nominal")
      return var;

    if (flagProcSyst) {
      //Check if var is contained in some envelope
      TString envelope_name = GetVarEnvelope(var);
      if (envelope_name!="") return envelope_name;
    }

    const char* s = var.Data();
    const int n = var.Length();
    if(TString(s + n - 3, 3) == "_UP")
    {
      sign = 1;
      return TString(s, n - 3);
    }
    else if(TString(s + n - 5, 5) == "_DOWN")
    {
      sign = -1;
      return TString(s, n - 5);
    }
    else if(var.Contains("PDF_") && var.Contains("CENTRAL"))
      return var;
    else if(var == "MASS_CONSTMC")
      return "MASS_CONSTMC";
    else if(var == "REWTOPPT")
      return "REWTOPPT";
    // else if(var == "PDF_101")
    //   return "PDF_NNPDF30";
    // else if(TString(var.Data(), 4) == "PDF_")
    // {
    //   int n = TString(var.Data() + 4).Atoi();
    //   if(n < 1 || n > 100)
    //     throw std::logic_error(TString::Format("Unsupported n = %d (var = %s)", n, var.Data()));
    //   int eig = (n-1) / 2 + 1;
    //   if(((n-1) % 2) == 0)
    //     sign = -1;
    //   else
    //     sign = +1;
    //   return TString::Format("PDF_NNPDF31_eig%d", eig);
    // }

    else if(var == "POWHEGV2HERWIG")
      return "POWHEGV2HERWIG";
    else if(var == "AMCATNLOFXFX")
      return "AMCATNLOFXFX";
    else if(var == "MADGRAPHMLM")
      return "MADGRAPHMLM";
    else if (var == "ERDON")
      return "ERDON";
    else if (var == "ERDONRETUNE")
      return "ERDONRETUNE";
    else if (var == "GLUONMOVETUNE")
      return "GLUONMOVETUNE";
    else if (var == "BFRAG_PETERSON")
      return "BFRAG_PETERSON";
    else if (var == "TOP_PT")
      return "TOP_PT";

    throw std::logic_error(TString::Format("Error in VarToSys(): unknown variation %s.", var.Data()).Data());
}

const TString PlotterConfigurationHelper::VarToSys(const TString& var)
{
    int sign = 0;
    return VarToSys(var, sign);
}

const TString PlotterConfigurationHelper::GetVarDirTree(const TString& var)
{
    TString oldVar = var;
    RemoveYearCorr(year_, yearCorr_, oldVar);

    if(oldVar.Contains("BG") || oldVar.Contains("LUMI") || oldVar.Contains("DY")) return "Nominal";
    return oldVar;
}

const TString PlotterConfigurationHelper::GetVarWeight(const TString& var){

    TString oldVar = var;
    RemoveYearCorr(year_, yearCorr_, oldVar);

    TString sys = VarToSys(var);

    if(IsExpSys(sys) || IsModWSys(sys) || IsModISys(sys)){
      if(!sys.BeginsWith("JES") && !sys.BeginsWith("JER") && !sys.BeginsWith("UNCLUSTERED") && !sys.BeginsWith("BG") && !sys.BeginsWith("LUMI")
	                        && !sys.BeginsWith("DY") && !sys.BeginsWith("CR") && !sys.BeginsWith("ERDON") && !sys.BeginsWith("ERDONRETUNE") && !sys.BeginsWith("GLUONMOVETUNE")
	                        && !sys.BeginsWith("MASS") && !sys.BeginsWith("MATCH") && !sys.BeginsWith("UETUNE") &&  !sys.BeginsWith("MUON_SCALE")
                                && !sys.BeginsWith("ELE_SCALESMEARING")
	 )

	return oldVar;
      else return "Nominal";
    }
    return "Nominal";
}

const TString PlotterConfigurationHelper::GetVarDirTreePlain(const TString& var, const bool flagWeights)
{

    TString oldVar = var;
    RemoveYearCorr(year_, yearCorr_, oldVar);

    const TString& sys = VarToSys(var);

    if(sys.Contains("DY") || sys.Contains("BG") || sys.Contains("LUMI"))
      return "Nominal";
    if(sys.Contains("REWTOPPT"))
      return "Nominal";
    if(flagWeights)
    {
      if(GetVarWeight(var) == "Nominal")
        return oldVar;
      else
        return "Nominal";
      /*if(IsModISys(var))
        return oldVar;
      else
        return "Nominal";*/
    }
    else
      return oldVar;
}

const TString PlotterConfigurationHelper::GetVarDirPlot(const TString& var)
{
    const TString& sys = VarToSys(var);
    //if(sys == "BR")
    if(sys.Contains("BR"))
      return "Nominal";
    return var;
}

TString PlotterConfigurationHelper::GetShortVarName(const TString& var)
{
    TString shortStr;


    TString labelYearCorr = GetLabelYearCorr(var);

    //if(IsUpDownVar(var))
    if(var.Contains("_UP") || var.Contains("_DOWN"))
      {
      TString str = VarToSys(var);
      RemoveYearCorr(year_, yearCorr_, str);

      if(str == "BTAG") shortStr = "b";
      else if(str == "BTAG_ETA") shortStr = "bet";
      else if(str == "BTAG_LJET") shortStr = "blj";
      else if(str == "BTAG_LJET_ETA") shortStr = "bljet";
      else if(str == "BTAG_LJET_PT") shortStr = "bljpt";
      else if(str == "BTAG_PT") shortStr = "bpt";
      //else if(str == "JER") shortStr = "jer";
      //else if(str == "JES") shortStr = "jes";
      else if(str == "KIN") shortStr = "kin";
      else if(str == "LEPT") shortStr = "lept";
      else if(str == "PU") shortStr = "pu";
      else if(str == "TRIG") shortStr = "tr";
      else if(str == "TRIG_ETA") shortStr = "treta";
      else if(str == "MEFACSCALE") shortStr = "mf";
      else if(str == "MERENSCALE") shortStr = "mr";
      else if(str == "MESCALE") shortStr = "mu";
      else if(str == "TOT_SCALE") shortStr = "tot_scale";
      else if(str == "PDF_ALPHAS") shortStr = "pdfalphas";
      else if(str == "MASS") shortStr = "mt";
      else if(str == "MATCH") shortStr = "match";
      else if(str == "PSFSRSCALE" || "PSSCALE_WEIGHT_5") shortStr = "psfsr";
      else if(str == "PSISRSCALE" || "PSSCALE_WEIGHT_4") shortStr = "psisr";
      else if(str == "PSFSRSCALE_2") shortStr = "psfsr2";
      else if(str == "PSISRSCALE_2") shortStr = "psisr2";
      else if(str == "UETUNE") shortStr = "ue";
      else if(str == "BSEMILEP") shortStr = "bslmp";
      else if(str == "BG") shortStr = "bkg";
      else if(str == "DY") shortStr = "dy";
      else if(str == "LUMI") shortStr = "lumi";
      else if(str == "BR") shortStr = "br";
      else if(str == "TOT_BFRAG") shortStr = "bfrag";
      else if(var == "BFRAG_PETERSON") shortStr = "bfragP";
      else if(var == "BFRAG_CENTRAL") shortStr = "bfragC";
      else if(var == "BFRAG_UP") shortStr = "bfragU";
      else if(var == "BFRAG_DOWN") shortStr = "bfragD";
      else if(str == "L1PREFIRING") shortStr = "l1pref";

      // CT14 PDF, as
      else if(str == "PDF_CT14_AS") shortStr = "pdfas";
      else if(TString(str.Data(), 12) == "PDF_CT14_eig") shortStr = "pdfe";
      // NNPDF30
      else if(str == "PDF_NNPDF30") shortStr = "pdfnn30";

      else if(str == "POWHEGV2HERWIG") shortStr = "pwhrw";
      else if(str == "AMCATNLOFXFX") shortStr = "amcat";

      // Nov17 splitting of JES sources
      else if(str == "JER") shortStr = "jer";
      else if(str == "UNCLUSTERED") shortStr = "unmet";
      else if(str == "METJES") shortStr = "metjs";
      else
      {
        for(unsigned int s  = 0; s < GetJECSrc(false).size(); s++)
        {
          const TString& src = GetJECSrc(false)[s];
          if(!src.Contains("Total"))
            if(src == str)
              shortStr = TString::Format("jes%2d", s);
        }
      }

      if(shortStr == "")
        throw std::logic_error(TString::Format("Error in GetShortVarName(): unknown variation name %s\n", var.Data()).Data());

      const char* s = var.Data();
      const int n = var.Length();
      if(TString(s + n - 3, 3) == "_UP")
	if (labelYearCorr != "")
	  shortStr += labelYearCorr + "U";
	else
	  shortStr += "U";
      else if(TString(s + n - 5, 5) == "_DOWN")
	if (labelYearCorr != "")
	  shortStr += labelYearCorr + "D";
	else
	  shortStr += "D";
      // PDF, as are named not with _UP, _DOWN
      else if(shortStr.BeginsWith("pdf"))
      {
        int n = TString(str.Data() + 4).Atoi();
        if(((n - 102) % 2) == 0)
	  if (labelYearCorr != "")
	    shortStr += labelYearCorr + "D";
	  else
	    shortStr += "D";
        else
	  if (labelYearCorr != "")
	    shortStr += labelYearCorr + "U";
	  else
	    shortStr += "U";
      }
      else
        throw std::logic_error(TString::Format("Error in GetShortVarName(): not _UP/_DOWN variation %s\n", var.Data()).Data());
    }

    else if(var == "BFRAG_UP") {
      if (labelYearCorr != "")
	shortStr = "bfrag" + labelYearCorr + "U";
      else
	shortStr = "bfragU";
    }
    else if(var == "BFRAG_DOWN") {
      if (labelYearCorr != "")
	shortStr = "bfrag" + labelYearCorr + "D";
      else
	shortStr = "bfragD";
    }
    else if(var == "BFRAG_PETERSON") {
      if (labelYearCorr != "")
	shortStr = "bfrag" + labelYearCorr + "P";
      else
	shortStr = "bfragP";
    }
    else if(var == "BFRAG_CENTRAL") {
      if (labelYearCorr != "")
	shortStr = "bfrag" + labelYearCorr + "C";
      else
	shortStr = "bfragC";
    }
    else if(var == "ERDON") {
      if (labelYearCorr != "")
	shortStr = "cr" + labelYearCorr + "1";
      else
	shortStr = "cr1";
    }
    else if(var == "ERDONRETUNE") {
      if (labelYearCorr != "")
	shortStr = "cr" + labelYearCorr + "2";
      else
	shortStr = "cr2";
    }
    else if(var == "GLUONMOVETUNE") {
      if (labelYearCorr != "")
	shortStr = "cr" + labelYearCorr + "3";
      else
	shortStr = "cr3";
    }
    else
      throw std::logic_error(TString::Format("Error in GetShortVarName(): unknown variation name %s\n", var.Data()).Data());

    return shortStr;
}

TString PlotterConfigurationHelper::GetShortVarNameForTex(const TString& var)
{
    TString shortStr;

    //if(IsUpDownVar(var))
    if(var.Contains("_UP") || var.Contains("_DOWN"))
      {
      TString str = VarToSys(var);

      if(str == "BTAG") shortStr = "b-tagging";
      else if(str == "BTAG_LJET") shortStr = "b-tagging (light jets)";
      else if(str == "BTAG_LJET_ETA") shortStr = "b-tagging (light jets in eta)";
      else if(str == "BTAG_LJET_PT") shortStr = "b-tagging (light jets in pt)";
      else if(str == "BTAG_ETA") shortStr = "b-tagging (eta)";
      else if(str == "BTAG_PT") shortStr = "b-tagging (pt)";
      //else if(str == "JER") shortStr = "JER";
      //else if(str == "JES") shortStr = "JES";
      else if(str == "LEPT") shortStr = "lepton ID/ISO";
      else if(str == "PU") shortStr = "pileup";
      else if(str == "TRIG") shortStr = "trigger";
      else if(str == "TRIG_ETA") shortStr = "trigger ($\\eta$)";
      else if(str == "MESCALE") shortStr = "$\\mu_{r,f}$";
      else if(str == "MASS") shortStr = "$m_{\\rm t}^{\\rm MC}$";
      else if(str == "MATCH") shortStr = "$h_{\\rm damp}$";
      else if(str == "PSFSRSCALE" || str == "PSFSRSCALE_2" || "PSSCALE_WEIGHT_5") shortStr = "PS FSR";
      else if(str == "PSISRSCALE" || str == "PSISRSCALE_2" || "PSSCALE_WEIGHT_4") shortStr = "PS ISR";
      else if(str == "UETUNE") shortStr = "UE tune";
      else if(str == "BSEMILEP") shortStr = "branching ratio $B\\to \\mu$";
      else if(str == "BG") shortStr = "non-\\ttbar background";
      else if(str == "DY") shortStr = "DY background";
      else if(str == "KIN") shortStr = "kinematic reconstruction";
      else if(str == "LUMI") shortStr = "luminosity";
      else if(str == "BR") shortStr = "branching ratio $\\ttbar \\to \\ell\\ell$";
      else if (str == "TOT_COLORREC") shortStr = "color reconnection";
      else if (str == "TOT_SCALE") shortStr = "$\\mu_{r,f}$";
      else if (str == "TOT_BFRAG") shortStr = "fragmentation";
      else if (str == "PDF_ALPHAS") shortStr = "$\\alpha_{\\rm s}$";
      else if(var == "BFRAG_UP") shortStr = "fragmentation up";
      else if(var == "BFRAG_DOWN") shortStr = "fragmentation down";
      else if(var == "BFRAG_PETERSON") shortStr = "fragmentation PETERSON";
      else if(var == "BFRAG_CENTRAL") shortStr = "fragmentation CENTRAL";
      else if(str == "L1PREFIRING") shortStr = "l1 prefiring";

      // CT14 PDF, as
      else if(str == "PDF_CT14_AS") shortStr = "$\\alpha_{\\rm s}$";
      else if(TString(str.Data(), 9) == "PDF_CT14_") shortStr = TString::Format("PDF eigenvector %d", (TString(str.Data() + 12).Atoi()));

      // Nov17 splitting of JES sources
      else if(str == "JER") shortStr = "JER";
      //else if(str == "UNCLUSTERED") shortStr = "$\\ETmiss$";
      else if(str == "UNCLUSTERED") shortStr = "$\\ptvecmiss$";
      else
      {
        for(unsigned int s  = 0; s < GetJECSrc(false).size(); s++)
        {
          const TString& src = GetJECSrc(false)[s];
          if(!src.Contains("Total"))
            if(src == str)
              //shortStr = TString::Format("JES%d", s);
              shortStr = src;
        }
      }

      if(shortStr == "")
        throw std::logic_error(TString::Format("Error in GetShortVarNameForTex(): unknown variation name %s\n", var.Data()).Data());

      const char* s = var.Data();
      const int n = var.Length();
      if(TString(s + n - 3, 3) == "_UP")
        shortStr += " v1";
      else if(TString(s + n - 5, 5) == "_DOWN")
        shortStr += " v2";
      // PDF, as are named not with _UP, _DOWN
      else if(shortStr.BeginsWith("PDF") || shortStr == "$\\alpha_{\\rm s}$")
      {
        int n = TString(str.Data() + 4).Atoi();
        if(((n - 102) % 2) == 0)
          shortStr += " v2";
        else
          shortStr += " v1";
      }
      else
        throw std::logic_error(TString::Format("Error in GetShortVarNameForTex(): not _UP/_DOWN variation %s\n", var.Data()).Data());
    }

    else if(var == "BFRAG_UP") shortStr = "fragmentation $b\\to B$ v1";
    else if(var == "BFRAG_DOWN") shortStr = "fragmentation $b\\to B$ v2";
    else if(var == "BFRAG_PETERSON") shortStr = "fragmentation $b\\to B$ v3";
    else if(var == "BFRAG_CENTRAL") shortStr = "fragmentation $b\\to B$ v4";
    else if(var == "ERDON") shortStr = "colour reconnection v1";
    else if(var == "ERDONRETUNE") shortStr = "colour reconnection v2";
    else if(var == "GLUONMOVETUNE") shortStr = "colour reconnection v3";
    else
      throw std::logic_error(TString::Format("Error in GetShortVarNameForTex(): unknown variation name %s\n", var.Data()).Data());

    return shortStr;
}

double PlotterConfigurationHelper::GetYearToYearCorrRescale(const TString& year, const TString& yearCorr, const TString& complete_Var_name)
{

    TString sysTmp = complete_Var_name;
    RemoveYearCorr(year, yearCorr, sysTmp, true);
    if(!IsYearToYearCorr(year, yearCorr, sysTmp))
      return 1.0;

    //Warning: Don't use this function blindly. Year-to-year correlation handling should be worked out based on the systematic source on a case by case basis
    //->the following is subject to specific interpretations of POG recommendations
    //->please modify accordingly

    std::vector<double> corrFactors = YearToYearCorrMap.at(sysTmp);

    //corrFactors[0] = 1617, corrFactors[1] = 1618, corrFactors[2] = 1718
    if(complete_Var_name.Contains("16") && !complete_Var_name.Contains("17") && !complete_Var_name.Contains("18")) {
      if(corrFactors[2] == 1. && corrFactors[0] == corrFactors[1])
	return TMath::Sqrt(1-corrFactors[0]);
      else if(corrFactors[0] == corrFactors[1] && corrFactors[0] == corrFactors[2])
	return TMath::Sqrt(1-corrFactors[0]);
      else if((corrFactors[0]+corrFactors[1])<=1.)
	return TMath::Sqrt(1-corrFactors[0]-corrFactors[1]);
    }
    else if(complete_Var_name.Contains("17") && !complete_Var_name.Contains("16") && !complete_Var_name.Contains("18")) {
      if(corrFactors[1] == 1. && corrFactors[0] == corrFactors[2])
	return TMath::Sqrt(1-corrFactors[0]);
      else if(corrFactors[0] == corrFactors[1] && corrFactors[0] == corrFactors[2])
	return TMath::Sqrt(1-corrFactors[0]);
      else if((corrFactors[0]+corrFactors[2])<=1.)
	return TMath::Sqrt(1-corrFactors[0]-corrFactors[2]);
    }
    else if(complete_Var_name.Contains("18") && !complete_Var_name.Contains("16") && !complete_Var_name.Contains("17")) {
      if(corrFactors[0] == 1. && corrFactors[1] == corrFactors[2])
	return TMath::Sqrt(1-corrFactors[1]);
      else if(corrFactors[0] == corrFactors[1] && corrFactors[0] == corrFactors[2])
	return TMath::Sqrt(1-corrFactors[0]);
      else if((corrFactors[1]+corrFactors[2])<=1.)
	  return TMath::Sqrt(1-corrFactors[1]-corrFactors[2]);
    }
    else if(complete_Var_name.Contains("1617") && !complete_Var_name.Contains("161718")) {
      if(corrFactors[0] == 1. && corrFactors[1] == corrFactors[2])
	return TMath::Sqrt(1-corrFactors[1]);
      else if(corrFactors[0] == corrFactors[1] && corrFactors[0] == corrFactors[2])
	return 0.;
      else if((corrFactors[0]+corrFactors[1])<=1. && (corrFactors[0]+corrFactors[2])<=1.)
	return TMath::Sqrt(corrFactors[0]);
    }
    else if(complete_Var_name.Contains("1618") && !complete_Var_name.Contains("161718")) {
      if(corrFactors[1] == 1. && corrFactors[0] == corrFactors[2])
	return TMath::Sqrt(1-corrFactors[0]);
      else if(corrFactors[0] == corrFactors[1] && corrFactors[0] == corrFactors[2])
	return 0.;
      else if((corrFactors[0]+corrFactors[1])<=1. && (corrFactors[1]+corrFactors[2])<=1.)
	return TMath::Sqrt(corrFactors[1]);
    }
    else if(complete_Var_name.Contains("1718") && !complete_Var_name.Contains("161718")) {
      if(corrFactors[2] == 1. && corrFactors[0] == corrFactors[1])
	return TMath::Sqrt(1-corrFactors[0]);
      else if(corrFactors[0] == corrFactors[1] && corrFactors[0] == corrFactors[2])
	return 0.;
      else if((corrFactors[0]+corrFactors[2])<=1. && (corrFactors[1]+corrFactors[2])<=1.)
	return TMath::Sqrt(corrFactors[2]);
    }
    else if(complete_Var_name.Contains("161718")) {
      if(corrFactors[0] == 1. && corrFactors[1] == corrFactors[2])
	return TMath::Sqrt(corrFactors[1]);
      else if(corrFactors[1] == 1. && corrFactors[0] == corrFactors[2])
	return TMath::Sqrt(corrFactors[0]);
      else if(corrFactors[2] == 1. && corrFactors[0] == corrFactors[1])
	return TMath::Sqrt(corrFactors[0]);
      else if(corrFactors[0] == corrFactors[1] && corrFactors[0] == corrFactors[2])
	return TMath::Sqrt(corrFactors[0]);
    }


    if(!((corrFactors[0]+corrFactors[1])<=1.) && !((corrFactors[0]+corrFactors[2])<=1.) && !((corrFactors[1]+corrFactors[2])<=1.) && !(corrFactors[0] == 1. && corrFactors[1] == corrFactors[2]) && !(corrFactors[1] == 1. && corrFactors[0] == corrFactors[2]) && !(corrFactors[2] == 1. && corrFactors[0] == corrFactors[1])) {
      std::cerr << "ERROR in PlotterConfigurationHelper::GetYearToYearCorrRescale() : correlation factors for " <<
	"1617 = " << corrFactors[0] << ", 1618 = " << corrFactors[1] << " and 1718 = " << corrFactors[2] << " are not allowed for " << complete_Var_name << std::endl;
      exit(1);
    }
    else
      return 0.;
}

double PlotterConfigurationHelper::GetSysRescale(const TString& sys)
{

    TString sysTmp = sys;
    RemoveYearCorr(year_, yearCorr_, sysTmp);

    if(sysTmp == "MASS")
      return 1.0 / massUncReductionFactor;
    /*if(sysTmp == "PSFSRSCALE" || sysTmp == "PSFSRSCALE_2"){
      static double val = 1.0 / TMath::Sqrt(2.0);
      return val;
    }*/
    if(IsModPDFEig(sysTmp)) // rescale CT14 unc. CL 90% -> 68%
      return 1.0 / 1.64;
    else
      return 1.0;
}

double PlotterConfigurationHelper::GetVarRescaleSample(const TString& var, TString mode, TString i_year)
{

    double lumi_error_for_year;
    if (i_year == "") {
      lumi_error_for_year = lumiError;
    }
    else {
      std::pair<double,double> lumi_error_for_year_map = get_lumi_and_err(i_year);
      lumi_error_for_year = lumi_error_for_year_map.second;
    }

    int sign = 0;
    double output_factor = 1.0;
    const TString& sys = VarToSys(var, sign);
    bool is_no_up_down = IsNoUpDown(var);
    bool is_envelope = IsEnvelope(sys);
    if(sign != 1 && sign != -1 && !is_no_up_down && !is_envelope) throw std::logic_error(TString::Format("Error in GetVarRescaleSample(): sign = %d (should be +-1), var = %s\n", sign, var.Data()));
    // for default mode
    if(sys.Contains("LUMI")){
      output_factor = 1.0 + sign * lumi_error_for_year;
    }
    //for bkg mode
    else if(mode == "bkg"){
      if (sys.Contains("BG")) output_factor = 1.0 + sign * gVarBkg;
      else if (sys.Contains("DY")) output_factor = 1.0 + sign * gVarDY;
    }
    else output_factor = 1.0;

    return output_factor;
}

// plain tree file suffices
TString PlotterConfigurationHelper::GetVarPlainTreeFileSuffix(const TString& var, const TString& complete_Var_name, const TString& use_year)
{
    if (var == "Nominal")
      return "";

    (void) use_year; //to supress warning
    /*
    //Handling different file name in 2016 when using fullRun2 (PS..._2_..)
    if (use_year == "2016" && var.Contains("PSISRSCALE_2") && var.Contains("_UP")) return "_psisrscaleup";
    else if (use_year == "2016" && var.Contains("PSISRSCALE_2") && var.Contains("_DOWN")) return "_psisrscaledown";
    else if (use_year == "2016" && var.Contains("PSFSRSCALE_2") && var.Contains("_UP")) return "_psfsrscaleup";
    else if (use_year == "2016" && var.Contains("PSFSRSCALE_2") && var.Contains("_DOWN")) return "_psfsrscaledown";*/

    if(!IsModISys(VarToSys(complete_Var_name)) && !IsAltTheorSys(complete_Var_name))
      return "";

    if(var.Contains("MATCH") && var.Contains("_UP"))
      return "_matchup";
    if(var.Contains("MATCH") && var.Contains("_DOWN"))
      return "_matchdown";
    if(var.Contains("MASS") && var.Contains("_UP"))
      return "_175_massup";
    if(var.Contains("MASS") && var.Contains("_DOWN"))
      return "_169_massdown";
    /*
    if(var.Contains("PSISRSCALE") && var.Contains("_UP"))
      return "_psisrscaleup";
    if(var.Contains("PSISRSCALE") && var.Contains("_DOWN"))
      return "_psisrscaledown";
    if(var.Contains("PSFSRSCALE") && var.Contains("_UP"))
      return "_psfsrscaleup";
    if(var.Contains("PSFSRSCALE") && var.Contains("_DOWN"))
      return "_psfsrscaledown";
    */
    if(var.Contains("UETUNE") && var.Contains("_UP"))
      return "_uetuneup";
    if(var.Contains("UETUNE") && var.Contains("_DOWN"))
      return "_uetunedown";
    if(var.Contains("ERDONRETUNE"))
      return "_erdonretune";
    if(var.Contains("ERDON"))
      return "_erdon";
    if(var.Contains("GLUONMOVETUNE"))
      return "_gluonmovetune";

    if(var.Contains("POWHEGV2HERWIG"))
      return "_powhegv2Herwig_part1";
    if(var.Contains("AMCATNLOFXFX"))
      return "_amcatnlofxfx";
    if(var.Contains("MADGRAPHMLM"))
      return "_madgraphmlm";

    if(var.Contains("MASS_CONSTMC") || (var.Contains("MASS_KINRECO") && var.Contains("_UP")) || (var.Contains("MASS_KINRECO") && var.Contains("_DOWN")) || (var.Contains("MASS_KINRECO3GEV") && var.Contains("_UP")) || (var.Contains("MASS_KINRECO3GEV") && var.Contains("_DOWN")) || var.Contains("REWTOPPT"))
      return "";

    throw std::logic_error(TString::Format("Error in GetPlainTreeFileSuffix(): unknown variation %s", var.Data()).Data());
}

TString PlotterConfigurationHelper::GetSysBkgDir(const TString& sys)
{
    if(IsExpSys(sys))
      return sys;
    else
      return "Nominal";
}

std::vector<TString> PlotterConfigurationHelper::GetChannelsPlainTree(const TString& ch)
{
    static std::map<TString, std::vector<TString> > mCh;
    if(mCh.size() == 0)
    {
      mCh.insert(std::pair<TString, std::vector<TString> >("ee", { "ee" }));
      mCh.insert(std::pair<TString, std::vector<TString> >("mumu", { "mumu" }));
      mCh.insert(std::pair<TString, std::vector<TString> >("emu", { "emu" }));
      mCh.insert(std::pair<TString, std::vector<TString> >("ll", { "emu", "ee", "mumu" }));
    }
    assert(mCh.find(ch) != mCh.end());
    return mCh[ch];
}

double PlotterConfigurationHelper::GetChannelBR(const TString& ch)
{
    static std::map<TString, double> mCh;
    if(mCh.size() == 0)
    {
      mCh.insert(std::pair<TString, double>("ee", gBrrTtbarToEMu / 2.0));
      mCh.insert(std::pair<TString, double>("mumu", gBrrTtbarToEMu / 2.0));
      mCh.insert(std::pair<TString, double>("emu", gBrrTtbarToEMu));
      mCh.insert(std::pair<TString, double>("ll", 2 * gBrrTtbarToEMu));
    }
    assert(mCh.find(ch) != mCh.end());
    return mCh[ch];
}

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>> files >>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
TreePlain* PlotterConfigurationHelper::GetTreeDataPlain(const TString& ch, const TString& dirSyst, const TString& wSyst)
{
    // channels
    std::vector<TString> vChDirs = GetChannelsPlainTree(ch);
    // prepare chain
    TreePlain* tree = new TreePlain(gPlainTreeNamePlain, false, wSyst);

    //Taking the same file names for all systematics
    TString def_syst= "";
    TString def_ch = "";

    for (TString i_year: years_list){
      std::vector<TString> gVFilesData;
      assert(gVFilesData.size() == 0);

      // Filling the sample vector
      fillPlainTreesFileListSampleVectors(i_year,_gVFilesData_pattern, gVFilesData,def_syst,def_ch);

      TString i_TreeDir = diLeptonic_dir + "mergedPlainTree_" + i_year;

      for(unsigned int i = 0; i < gVFilesData.size(); i++){
	for(const auto& chDir : vChDirs){
	  if( (chDir == "ee" && gVFilesData[i].Contains("smu_")) || (chDir == "mumu" && gVFilesData[i].Contains("se_"))) continue;
	  tree->Add(i_TreeDir + "/" + dirSyst + "/" + chDir + "/" + chDir + "_" + gVFilesData[i]);
	}
      }

    }

    tree->InitReading();

    return tree;
}

TreePlain* PlotterConfigurationHelper::GetTreeMCBgrPlain(const TString& complete_Var_name, const TString& ch, const TString& dirSyst, const TString& wSyst)
{
    // channels
    std::vector<TString> vChDirs = GetChannelsPlainTree(ch);

    TString var_without_yearCorrLabel = complete_Var_name;
    if (yearCorr_ != "") RemoveYearCorr(year_, yearCorr_, var_without_yearCorrLabel);

    TString sys = VarToSys(complete_Var_name);

    double bkg_rescale = 1.0;
    if (complete_Var_name.Contains("DY_") || complete_Var_name.Contains("BG_")) {
      bkg_rescale = GetVarRescaleSample(complete_Var_name, "bkg");
    }

    // prepare chain
    TreePlain* tree = new TreePlain(gPlainTreeNamePlain, false, wSyst);

    //Taking the same file names for all systematics
    TString def_syst= "";
    TString def_ch = "";

    std::vector<TString> list_years_Corr;

    //check which years should be varied
    if (complete_Var_name.Contains("16")) {
      list_years_Corr.push_back("2016");
    }
    if (complete_Var_name.Contains("17")) {
      list_years_Corr.push_back("2017");
    }
    if (complete_Var_name.Contains("18")) {
      list_years_Corr.push_back("2018");
    }

    std::vector<TString> other_years;
    if ((year_ == "fullRun2UL") && yearCorr_ != "" && IsYearToYearCorr(year_, yearCorr_, complete_Var_name) && list_years_Corr.size() != 0) {
      other_years = this->years_list;
      for(auto iYearCorr = list_years_Corr.begin(); iYearCorr != list_years_Corr.end(); ++iYearCorr) {
	int idx;
	auto it = std::find(other_years.begin(), other_years.end(), *iYearCorr);
	if (it != other_years.end())
	  {
	    idx = std::distance(other_years.begin(), it);
	    other_years.erase(other_years.begin()+idx);
	  }
      }
    }

    for (TString i_year: this->years_list) {

      TString dirSystTmp = dirSyst;
      double rescaleTmp = GetVarRescaleSample(complete_Var_name, "", i_year);
      // bool apply_tlrw_to_i_year = true;

      if (yearCorr_ != "" && IsYearToYearCorr(year_, yearCorr_, complete_Var_name) && (std::find(other_years.begin(), other_years.end(), i_year) != other_years.end())) {
          dirSystTmp = "Nominal";
          // apply_tlrw_to_i_year = false;
      }


      double year_to_year_Corr_rescale = 1.;
      bool do_year_to_year_Corr_rescale = false;
      if ((dirSystTmp == "Nominal" && complete_Var_name == "Nominal") || yearCorr_ == "") {
	//case 1
	//no correlation treatment should be performed for these two cases: running Nominal or not running with --doCorrYear. In the latter case all systematics are treated as 100%        correlated between years by default
	year_to_year_Corr_rescale = 1.;
	do_year_to_year_Corr_rescale = false;
      }
      else if (yearCorr_ != "" && dirSystTmp == "Nominal" && (std::find(other_years.begin(), other_years.end(), i_year) != other_years.end())) {
	//case 2
	//sample dir should be nominal for years not contained in name of systematic. I.e. for BTAG_1617 one needs to take nominal samples for 2018
	year_to_year_Corr_rescale = 1.;
	do_year_to_year_Corr_rescale = false;
	rescaleTmp = 1.;
      }
      else if (yearCorr_ != "" && !IsExpSys(sys)) {
	//case 3
	//if sys isn't experimental then don't rescale bg samples
	year_to_year_Corr_rescale = 1.;
	do_year_to_year_Corr_rescale = false;
      }
      else {
	//case 4
        year_to_year_Corr_rescale = GetYearToYearCorrRescale(year_, yearCorr_, complete_Var_name);
	do_year_to_year_Corr_rescale = true;
      }

      //for systematics stored as weights in the nominal samples one needs to set doWeightVar To false for case 2
      bool doWeightVar = false;
      if ((yearCorr_ != "" && wSyst != "" && wSyst != "Nominal" && do_year_to_year_Corr_rescale) || (yearCorr_ == "" && wSyst != "" && wSyst != "Nominal"))
	doWeightVar = true;

      TString i_TreeDir = diLeptonic_dir + "mergedPlainTree_" + i_year;

      std::vector<TString> gVFilesMCbkg_noID;
      assert(gVFilesMCbkg_noID.size() == 0);
      std::vector<std::pair<TString, int> > gVFilesMCbkg;
      assert(gVFilesMCbkg.size() == 0);

      fillPlainTreesFileListSampleVectors(i_year,_gVFilesMCbkg_pattern, gVFilesMCbkg_noID,def_syst,def_ch);
      fillPlainTreesSampleIdVector(gVFilesMCbkg_noID,_gVFilesMCbkg_pattern_ID,gVFilesMCbkg);

      for(unsigned int i = 0; i < gVFilesMCbkg.size(); i++){
	for(const auto& chDir : vChDirs){

	  TString dirSystTmpTEMP = dirSystTmp;
	  bool doWeightVarTEMP = doWeightVar;
	  bool do_year_to_year_Corr_rescaleTEMP = do_year_to_year_Corr_rescale;
	  double year_to_year_Corr_rescaleTEMP = year_to_year_Corr_rescale;

	  //temporary fix for wrong handling of certain systematics that should have zero variation. TODO: fix problem at file production level
	  if ((chDir == "ee" && complete_Var_name.Contains("MUON"))
	      || (chDir == "mumu" && complete_Var_name.Contains("ELE"))
	      || (year_ == "2018" && complete_Var_name.Contains("L1PREFIRING")))
	    {
	      dirSystTmpTEMP = "Nominal";
	      doWeightVarTEMP = false;
	      do_year_to_year_Corr_rescaleTEMP = false;
	      year_to_year_Corr_rescaleTEMP = 1.;
	    }

	  TString nameFile = i_TreeDir + "/" + dirSystTmpTEMP + "/" + chDir + "/" + chDir + "_" + gVFilesMCbkg[i].first;
	  TString nameFile_nominal = i_TreeDir + "/Nominal/" + chDir + "/" + chDir + "_" + gVFilesMCbkg[i].first;

	  double w = CalcLumiWeight(nameFile,i_year) * rescaleTmp;
	  double w_nominal = CalcLumiWeight(nameFile_nominal,i_year);
        
        bool running_dy = false;
        if (gVFilesMCbkg[i].second == 9 || gVFilesMCbkg[i].second == 5) running_dy = true;

	  // if (ttmd_apply_trueLevelWeights && apply_tlrw_to_i_year && !running_dy) w *= GetTrueLevelRenormalisationWeight(nameFile,var_without_yearCorrLabel);

	  //Handling DY scale factor and systemtic
	  if (running_dy){
	    if (doDYScale) {
	      w *= get_DY_SF(i_year, ch);
	      w_nominal *= get_DY_SF(i_year, ch);
	      std::cout << "(** DY scale Factor = " << get_DY_SF(i_year, ch) << " **)" ;
	    }
	    if (complete_Var_name.Contains("DY_") && (do_year_to_year_Corr_rescaleTEMP || yearCorr_ == "")){
	      w *= bkg_rescale;
	      std::cout << "(rescaling to " << complete_Var_name << " with weight = " << bkg_rescale << ")";
	    }
	    std::cout << " -->> ";
	  }
	  //Handling no DY background systemtic
	  else if (complete_Var_name.Contains("BG_") && (do_year_to_year_Corr_rescaleTEMP || yearCorr_ == "")) {
	    w *= bkg_rescale;
	    std::cout << "(** rescaling to " << complete_Var_name << " with weight = " << bkg_rescale << " **) -->> ";
	  }

	  printf("Applying year-to-year correlation rescale factor: %f \n",year_to_year_Corr_rescaleTEMP);
	  w *= year_to_year_Corr_rescaleTEMP;

	  int sampleType = gVFilesMCbkg[i].second;
	  tree->Add(nameFile, w, doWeightVarTEMP, sampleType);
	  //make sure that the correlation rescale factor is only applied to the variation itself and not to nominal
	  if(year_to_year_Corr_rescaleTEMP != 1.) {
	    double w_nominal_rescaled = w_nominal-year_to_year_Corr_rescaleTEMP*w_nominal;
	    tree->Add(nameFile_nominal, w_nominal_rescaled, false, sampleType);
	  }
	}
      }
    }

    tree->InitReading();

    return tree;
}

TreePlain* PlotterConfigurationHelper::GetTreeMCBttPlain(const TString& complete_Var_name, const TString& ch, const TString& dirSyst, const TString& wSyst, const bool flagMode)
{
    // channels
    std::vector<TString> vChDirs = GetChannelsPlainTree(ch);

    TString var_without_yearCorrLabel = complete_Var_name;
    if (yearCorr_ != "") RemoveYearCorr(year_, yearCorr_, var_without_yearCorrLabel);

    TString sys = VarToSys(complete_Var_name);
    TString suffix = GetVarPlainTreeFileSuffix(dirSyst, complete_Var_name);
    TString dir = (suffix != "" && !IsModISys(sys) && flagMode) ? "Nominal" : dirSyst;

    TreePlain* tree = new TreePlain(gPlainTreeNamePlain, false, wSyst);

    //Taking the same file names for all systematics
    TString def_syst= "";
    TString def_ch = "";

    std::vector<TString> list_years_Corr;

    //check which years should be varied
    if (complete_Var_name.Contains("16")) {
      list_years_Corr.push_back("2016");
    }
    if (complete_Var_name.Contains("17")) {
      list_years_Corr.push_back("2017");
    }
    if (complete_Var_name.Contains("18")) {
      list_years_Corr.push_back("2018");
    }

    std::vector<TString> other_years;
    if (year_ == "fullRun2UL" && yearCorr_ != "" && IsYearToYearCorr(year_, yearCorr_, complete_Var_name) && list_years_Corr.size() != 0) {
      other_years = this->years_list;
      for(auto iYearCorr = list_years_Corr.begin(); iYearCorr != list_years_Corr.end(); ++iYearCorr) {
	int idx;
	auto it = std::find(other_years.begin(), other_years.end(), *iYearCorr);
	if (it != other_years.end())
	  {
	    idx = std::distance(other_years.begin(), it);
	    other_years.erase(other_years.begin()+idx);
	  }
      }
    }

    for (TString i_year: this->years_list){

      TString dirSystTmp = dirSyst;
      TString dirTmp = dir;
      double rescaleTmp = GetVarRescaleSample(complete_Var_name, "", i_year);
      bool apply_tlrw_to_i_year = true;

      if (yearCorr_ != "" && IsYearToYearCorr(year_, yearCorr_, complete_Var_name) && (std::find(other_years.begin(), other_years.end(), i_year) != other_years.end())) {
          dirSystTmp = "Nominal";
          dirTmp = "Nominal";
          apply_tlrw_to_i_year = false;
      }

      double year_to_year_Corr_rescale = 1.;
      bool do_year_to_year_Corr_rescale = false;
      if ((dirSystTmp == "Nominal" && complete_Var_name == "Nominal") || yearCorr_ == "") {
	//case 1
	//no correlation treatment should be performed for these two cases: running Nominal or not running with --doCorrYear. In the latter case all systematics are treated as 100%        correlated between years by default
	year_to_year_Corr_rescale = 1.;
	do_year_to_year_Corr_rescale = false;
      }
      else if (yearCorr_ != "" && dirSystTmp == "Nominal" && (std::find(other_years.begin(), other_years.end(), i_year) != other_years.end())) {
	//case 2
	//sample dir should be nominal for years not contained in name of systematic. I.e. for BTAG_1617 one needs to take nominal samples for 2018
	year_to_year_Corr_rescale = 1.;
	do_year_to_year_Corr_rescale = false;
	rescaleTmp = 1.;
      }
      else {
	//case 3
        year_to_year_Corr_rescale = GetYearToYearCorrRescale(year_, yearCorr_, complete_Var_name);
	do_year_to_year_Corr_rescale = true;
      }

      //for systematics stored as weights in the nominal samples one needs to set doWeightVar to false for case 2
      bool doWeightVar = false;
      if ((yearCorr_ != "" && wSyst != "" && wSyst != "Nominal" && do_year_to_year_Corr_rescale) || (yearCorr_ == "" && wSyst != "" && wSyst != "Nominal"))
	doWeightVar = true;

      suffix = GetVarPlainTreeFileSuffix(dirSystTmp, complete_Var_name, i_year);
      TString i_TreeDir = diLeptonic_dir + "mergedPlainTree_" + i_year;

      std::vector<TString> ttbar_bgviatau;
      assert(ttbar_bgviatau.size() == 0);
      std::vector<TString> ttbar_bkg;
      assert(ttbar_bkg.size() == 0);
      std::vector<TString> gVFilesMCtt;
      assert(gVFilesMCtt.size() == 0);
      std::vector<TString> gVFilesMCtt_nominal ;
      assert(gVFilesMCtt_nominal.size() == 0);

      fillPlainTreesFileListSampleVectors(i_year, _ttbar_bgviatau_pattern, ttbar_bgviatau, def_syst,def_ch);
      fillPlainTreesFileListSampleVectors(i_year, _ttbar_bkg_pattern, ttbar_bkg, def_syst,def_ch);

      gVFilesMCtt = ttbar_bgviatau;
      //including dedicated bkg samples
      for (auto bkg: ttbar_bkg) {
	gVFilesMCtt.push_back(bkg);
      }

      gVFilesMCtt_nominal = gVFilesMCtt;

      for(unsigned int i = 0; i < gVFilesMCtt.size(); i++) {
	if (suffix != "" && IsModISys(sys) && (gVFilesMCtt[i].Contains("_PSweights.root") || gVFilesMCtt[i].Contains("_TuneCP5.root"))){
	  if (gVFilesMCtt[i].Contains("_PSweights.root")) gVFilesMCtt[i].ReplaceAll("_PSweights.root", suffix + ".root");
	  else if (gVFilesMCtt[i].Contains("_TuneCP5.root")) gVFilesMCtt[i].ReplaceAll("_TuneCP5.root", suffix + ".root");
	}
	else{
	  gVFilesMCtt[i].ReplaceAll(".root", suffix + ".root");
	}
      }

      for(unsigned int i = 0; i < gVFilesMCtt.size(); i++){
	for(const auto& chDir : vChDirs){

	  TString dirTmpTEMP = dirTmp;
	  bool doWeightVarTEMP = doWeightVar;
	  bool do_year_to_year_Corr_rescaleTEMP = do_year_to_year_Corr_rescale;
	  double year_to_year_Corr_rescaleTEMP = year_to_year_Corr_rescale;

	  //temporary fix for wrong handling of certain systematics that should have zero variation. TODO: fix problem at file production level
	  if ((chDir == "ee" && complete_Var_name.Contains("MUON"))
	      || (chDir == "mumu" && complete_Var_name.Contains("ELE"))
	      || (year_ == "2018" && complete_Var_name.Contains("L1PREFIRING")))
	    {
	      dirTmpTEMP = "Nominal";
	      doWeightVarTEMP = false;
	      do_year_to_year_Corr_rescaleTEMP = false;
	      year_to_year_Corr_rescaleTEMP = 1.;
	    }

	  TString fileName = i_TreeDir + "/" + dirTmpTEMP + "/" + chDir + "/" + chDir + "_" + gVFilesMCtt[i];
	  TString fileName_nominal = i_TreeDir + "/Nominal/" + chDir + "/" + chDir + "_" + gVFilesMCtt_nominal[i];
	  double w = CalcLumiWeight(fileName,i_year) * rescaleTmp;
	  double w_nominal = CalcLumiWeight(fileName_nominal,i_year);
	  printf("Applying year-to-year correlation rescale factor: %f \n",year_to_year_Corr_rescaleTEMP);
	  w *= year_to_year_Corr_rescaleTEMP;

	  if (ttmd_apply_trueLevelWeights && apply_tlrw_to_i_year) w *= GetTrueLevelRenormalisationWeight(fileName,var_without_yearCorrLabel);

	  tree->Add(fileName, w, doWeightVarTEMP);
	  //printf("adding %s %d\n", TString(gTreeDir + "/" + dirTmpTEMP + "/emu/emu_" + gVFilesMCtt[i]).Data(), tree->TChainPtr->GetEntries());
	  //make sure that the correlation rescale factor is only applied to the variation itself and not to nominal
	  if(year_to_year_Corr_rescaleTEMP != 1.) {
	    double w_nominal_rescaled = w_nominal-year_to_year_Corr_rescaleTEMP*w_nominal;
	    tree->Add(fileName_nominal, w_nominal_rescaled, false);
	  }
	}
      }
    }

    tree->InitReading();

    return tree;
}

TreePlain* PlotterConfigurationHelper::GetTreeMCRecPlain(const TString& complete_Var_name, const TString& ch, const TString& dirSyst, const TString& wSyst, const bool flagMode, const bool DoParticle)
{
    // channels
    std::vector<TString> vChDirs = GetChannelsPlainTree(ch);

    TString var_without_yearCorrLabel = complete_Var_name;
    if (yearCorr_ != "") RemoveYearCorr(year_, yearCorr_, var_without_yearCorrLabel);

    TString sys = VarToSys(complete_Var_name);
    TString suffix = GetVarPlainTreeFileSuffix(dirSyst, complete_Var_name);
    TString dir = (suffix != "" && !IsModISys(sys) && flagMode) ? "Nominal" : dirSyst;

    TreePlain* tree = new TreePlain(gPlainTreeNamePlain, false, wSyst, DoParticle);

    //Taking the same file names for all systematics
    TString def_syst= "";
    TString def_ch = "";

    std::vector<TString> list_years_Corr;

    //check which years should be varied
    if (complete_Var_name.Contains("16")) {
      list_years_Corr.push_back("2016");
    }
    if (complete_Var_name.Contains("17")) {
      list_years_Corr.push_back("2017");
    }
    if (complete_Var_name.Contains("18")) {
      list_years_Corr.push_back("2018");
    }

    std::vector<TString> other_years;
    if ((year_ == "fullRun2UL") && yearCorr_ != "" && IsYearToYearCorr(year_, yearCorr_, complete_Var_name) && list_years_Corr.size() != 0) {
      other_years = this->years_list;
      for(auto iYearCorr = list_years_Corr.begin(); iYearCorr != list_years_Corr.end(); ++iYearCorr) {
	int idx;
	auto it = std::find(other_years.begin(), other_years.end(), *iYearCorr);
	if (it != other_years.end())
	  {
	    idx = std::distance(other_years.begin(), it);
	    other_years.erase(other_years.begin()+idx);
	  }
      }
    }

    for (TString i_year: this->years_list) {

      TString dirSystTmp = dirSyst;
      TString dirTmp = dir;
      double rescaleTmp = GetVarRescaleSample(complete_Var_name, "", i_year);
      bool apply_tlrw_to_i_year = true;
      if (yearCorr_ != "" && IsYearToYearCorr(year_, yearCorr_, complete_Var_name) && (std::find(other_years.begin(), other_years.end(), i_year) != other_years.end())) {
          dirSystTmp = "Nominal";
          dirTmp = "Nominal";
          apply_tlrw_to_i_year = false;
      }

      double year_to_year_Corr_rescale = 1.;
      bool do_year_to_year_Corr_rescale = false;
      if ((dirSystTmp == "Nominal" && complete_Var_name == "Nominal") || yearCorr_ == "") {
	//case 1
	//no correlation treatment should be performed for these two cases: running Nominal or not running with --doCorrYear. In the latter case all systematics are treated as 100%        correlated between years by default
	year_to_year_Corr_rescale = 1.;
	do_year_to_year_Corr_rescale = false;
      }
      else if (yearCorr_ != "" && dirSystTmp == "Nominal" && (std::find(other_years.begin(), other_years.end(), i_year) != other_years.end())) {
	//case 2
	//sample dir should be nominal for years not contained in name of systematic. I.e. for BTAG_1617 one needs to take nominal samples for 2018
	year_to_year_Corr_rescale = 1.;
	do_year_to_year_Corr_rescale = false;
	rescaleTmp = 1.;
      }
      else {
	//case 3
	year_to_year_Corr_rescale = GetYearToYearCorrRescale(year_, yearCorr_, complete_Var_name);
	do_year_to_year_Corr_rescale = true;
      }


      //for systematics stored as weights in the nominal samples one needs to set doWeightVar to false for case 2
      bool doWeightVar = false;
      if ((yearCorr_ != "" && wSyst != "" && wSyst != "Nominal" && do_year_to_year_Corr_rescale) || (yearCorr_ == "" && wSyst != "" && wSyst != "Nominal"))
	doWeightVar = true;

      suffix = GetVarPlainTreeFileSuffix(dirSystTmp, complete_Var_name, i_year);
      std::vector<TString> gVFilesMCRec;
      std::vector<TString> gVFilesMCRec_nominal;
      assert(gVFilesMCRec.size()==0);
      assert(gVFilesMCRec_nominal.size()==0);

      fillPlainTreesFileListSampleVectors(i_year,_ttbar_signal_pattern,gVFilesMCRec,def_syst,def_ch);
      gVFilesMCRec_nominal = gVFilesMCRec;
      for(unsigned int i = 0; i < gVFilesMCRec.size(); i++) {
	if (suffix != "" && IsModISys(sys) && (gVFilesMCRec[i].Contains("_PSweights.root") || gVFilesMCRec[i].Contains("_TuneCP5.root"))){
	  if (gVFilesMCRec[i].Contains("_PSweights.root")) gVFilesMCRec[i].ReplaceAll("_PSweights.root", suffix + ".root");
	  else if (gVFilesMCRec[i].Contains("_TuneCP5.root")) gVFilesMCRec[i].ReplaceAll("_TuneCP5.root", suffix + ".root");
	}
	else{
	  gVFilesMCRec[i].ReplaceAll(".root", suffix + ".root");
	}
      }

      TString i_TreeDir = diLeptonic_dir + "mergedPlainTree_" + i_year;
      //year_to_year_Corr_rescale = 1.0000001;

      for (unsigned int i = 0; i < gVFilesMCRec.size(); i++) {
	for(const auto& chDir : vChDirs){

	  TString dirTmpTEMP = dirTmp;
	  bool doWeightVarTEMP = doWeightVar;
	  bool do_year_to_year_Corr_rescaleTEMP = do_year_to_year_Corr_rescale;
	  double year_to_year_Corr_rescaleTEMP = year_to_year_Corr_rescale;

	  //temporary fix for wrong handling of certain systematics that should have zero variation. TODO: fix problem at file production level
	  if ((chDir == "ee" && complete_Var_name.Contains("MUON"))
	      || (chDir == "mumu" && complete_Var_name.Contains("ELE"))
	      || (year_ == "2018" && complete_Var_name.Contains("L1PREFIRING")))
	    {
	      dirTmpTEMP  = "Nominal";
	      doWeightVarTEMP = false;
	      do_year_to_year_Corr_rescaleTEMP = false;
	      year_to_year_Corr_rescaleTEMP = 1.;
	    }

	  TString fileName = i_TreeDir + "/" + dirTmpTEMP  + "/" + chDir + "/" + chDir + "_" + gVFilesMCRec[i];
	  TString fileName_nominal = i_TreeDir + "/Nominal/" + chDir + "/" + chDir + "_" + gVFilesMCRec_nominal[i];
	  double w = CalcLumiWeight(fileName,i_year) * rescaleTmp;
	  double w_nominal = CalcLumiWeight(fileName_nominal,i_year);
	  printf("Applying year-to-year correlation rescale factor: %f \n",year_to_year_Corr_rescaleTEMP);

	  w *= year_to_year_Corr_rescaleTEMP;

	  if (ttmd_apply_trueLevelWeights && apply_tlrw_to_i_year) w *= GetTrueLevelRenormalisationWeight(fileName,var_without_yearCorrLabel);

	  tree->Add(fileName, w, doWeightVarTEMP);
	  //printf("GetTreeMCRec(): w = %f = %f %f / %ld\n", w, gXsecTtbarTwiki, lumi, nevt);
	  //make sure that the correlation rescale factor is only applied to the variation itself and not to nominal
	  if(year_to_year_Corr_rescaleTEMP != 1.) {
	    double w_nominal_rescaled = w_nominal-year_to_year_Corr_rescaleTEMP*w_nominal;
	    tree->Add(fileName_nominal, w_nominal_rescaled, false);
	  }
	}
      }
    }

    tree->InitReading();

    return tree;
}

TreePlain* PlotterConfigurationHelper::GetTreeMCGenPlain(const TString& complete_Var_name, const TString& ch, const TString& dirSyst, const TString& wSyst, const bool flagMode, const bool DoParticle)
{
    // channels
    std::vector<TString> vChDirs = GetChannelsPlainTree(ch);

    TString var_without_yearCorrLabel = complete_Var_name;
    if (yearCorr_ != "") RemoveYearCorr(year_, yearCorr_, var_without_yearCorrLabel);

    TString sys = VarToSys(complete_Var_name);
    TString varWeightName = IsModWSys(sys) || sys.Contains("PU") ? wSyst : "Nominal";
    TString suffix = GetVarPlainTreeFileSuffix(dirSyst, complete_Var_name);
    bool isAltTheory = IsAltTheorSys(sys);

    TString dir = (suffix != "" && !IsModISys(sys) && !isAltTheory && flagMode) ? "Nominal" : dirSyst;
    TreePlain* tree = new TreePlain(gPlainTreeName0Plain, true, varWeightName, DoParticle);

    //Taking the same file names for all systematics
    TString def_syst= "";
    TString def_ch = "";

    std::vector<TString> list_years_Corr;

    if (complete_Var_name.Contains("16")) {
      list_years_Corr.push_back("2016");
    }
    if (complete_Var_name.Contains("17")) {
      list_years_Corr.push_back("2017");
    }
    if (complete_Var_name.Contains("18")) {
      list_years_Corr.push_back("2018");
    }

    std::vector<TString> other_years;
    if ((year_ == "fullRun2UL") && yearCorr_ != "" && IsYearToYearCorr(year_, yearCorr_, complete_Var_name) && list_years_Corr.size() != 0) {
      other_years = this->years_list;
      for(auto iYearCorr = list_years_Corr.begin(); iYearCorr != list_years_Corr.end(); ++iYearCorr) {
	int idx;
	auto it = std::find(other_years.begin(), other_years.end(), *iYearCorr);
	if (it != other_years.end())
	  {
	    idx = std::distance(other_years.begin(), it);
	    other_years.erase(other_years.begin()+idx);
	  }
      }
    }

    for (TString i_year: this->years_list){

      TString dirSystTmp = dirSyst;
      TString dirTmp = dir;
      double rescaleTmp = 1.;
      if(!complete_Var_name.Contains("LUMI"))
	rescaleTmp = GetVarRescaleSample(complete_Var_name, "", i_year);
      bool apply_tlrw_to_i_year = true; // apply true lvl ren. weights only to the years with the variation.

      if (yearCorr_ != "" && IsYearToYearCorr(year_, yearCorr_, complete_Var_name) && (std::find(other_years.begin(), other_years.end(), i_year) != other_years.end())) {
          dirSystTmp = "Nominal";
          dirTmp = "Nominal";
          apply_tlrw_to_i_year = false;
      }

      double year_to_year_Corr_rescale = 1.;
      bool do_year_to_year_Corr_rescale = false;
      if ((dirSystTmp == "Nominal" && complete_Var_name == "Nominal") || yearCorr_ == "") {
	//case 1
	//no correlation treatment should be performed for these two cases: running Nominal or not running with --doCorrYear. In the latter case all systematics are treated as 100%        correlated between years by default
	year_to_year_Corr_rescale = 1.;
	do_year_to_year_Corr_rescale = false;
      }
      else if (yearCorr_ != "" && dirSystTmp == "Nominal" && (std::find(other_years.begin(), other_years.end(), i_year) != other_years.end())) {
	//case 2
	//sample dir should be nominal for years not contained in name of systematic. I.e. for BTAG_1617 one needs to take nominal samples for 2018
	year_to_year_Corr_rescale = 1.;
	do_year_to_year_Corr_rescale = false;
	rescaleTmp = 1.;
      }
      else {
	//case 3
        year_to_year_Corr_rescale = GetYearToYearCorrRescale(year_, yearCorr_, complete_Var_name);
	do_year_to_year_Corr_rescale = true;
      }

      //for systematics stored as weights in the nominal samples one needs to set doWeightVar to false for case 2
      bool doWeightVar = false;
      if ((yearCorr_ != "" && varWeightName != "" && varWeightName != "Nominal" && do_year_to_year_Corr_rescale) || (yearCorr_ == "" && wSyst != "" && wSyst != "Nominal"))
	doWeightVar = true;

      suffix = GetVarPlainTreeFileSuffix(dirSystTmp, complete_Var_name, i_year);
      std::vector<TString> gVFilesMCGen;
      std::vector<TString> gVFilesMCGen_nominal;
      assert(gVFilesMCGen.size()==0);
      assert(gVFilesMCGen_nominal.size()==0);

      if (isAltTheory) {
	if (i_year == "2016preVFP") {
	  gVFilesMCGen.push_back("ttbarsignalplustau_fromDilepton" + suffix + ".root");
	  gVFilesMCGen_nominal.push_back("ttbarsignalplustau_fromDilepton_TuneCP5.root");
	}
	if (i_year == "2016postVFP") {
	  gVFilesMCGen.push_back("ttbarsignalplustau_fromDilepton" + suffix + ".root");
	  gVFilesMCGen_nominal.push_back("ttbarsignalplustau_fromDilepton_TuneCP5.root");
	}
	else if (i_year == "2017") {
	  gVFilesMCGen.push_back("ttbarsignalplustau" + suffix + ".root");
	  gVFilesMCGen_nominal.push_back("ttbarsignalplustau_fromDilepton_PSweights.root");
	}
	else {
	  gVFilesMCGen.push_back("ttbarsignalplustau" + suffix + ".root");
	  gVFilesMCGen_nominal.push_back("ttbarsignalplustau_fromDilepton.root");
	}
      }
      else {
	fillPlainTreesFileListSampleVectors(i_year,_ttbar_signal_pattern,gVFilesMCGen,def_syst,def_ch);
	gVFilesMCGen_nominal = gVFilesMCGen;
	for(unsigned int i = 0; i < gVFilesMCGen.size(); i++) {
	  if (suffix != "" && IsModISys(sys) && (gVFilesMCGen[i].Contains("_PSweights.root") || gVFilesMCGen[i].Contains("_TuneCP5.root"))){
	    if (gVFilesMCGen[i].Contains("_PSweights.root")) gVFilesMCGen[i].ReplaceAll("_PSweights.root", suffix + ".root");
	    else if (gVFilesMCGen[i].Contains("_TuneCP5.root")) gVFilesMCGen[i].ReplaceAll("_TuneCP5.root", suffix + ".root");
	  }
	  else{
	    gVFilesMCGen[i].ReplaceAll(".root", suffix + ".root");
	  }
	}
      }

      TString i_TreeDir = diLeptonic_dir + "mergedPlainTree_" + i_year;

      for(const auto& chDir : vChDirs){

	TString dirTmpTEMP = dirTmp;
	bool doWeightVarTEMP = doWeightVar;
	bool do_year_to_year_Corr_rescaleTEMP = do_year_to_year_Corr_rescale;
        double year_to_year_Corr_rescaleTEMP = year_to_year_Corr_rescale;

	//temporary fix for wrong handling of certain systematics that should have zero variation. TODO: fix problem at file production level
	if ((chDir == "ee" && complete_Var_name.Contains("MUON"))
	    || (chDir == "mumu" && complete_Var_name.Contains("ELE"))
	    || (year_ == "2018" && complete_Var_name.Contains("L1PREFIRING")))
	  {
	    dirTmpTEMP = "Nominal";
	    doWeightVarTEMP = false;
	    do_year_to_year_Corr_rescaleTEMP = false;
	    year_to_year_Corr_rescaleTEMP = 1.;
	  }

	for (unsigned int i = 0; i < gVFilesMCGen.size(); i++){
	  TString fileName = i_TreeDir + "/" + dirTmpTEMP + "/" + chDir + "/" + chDir + "_" + gVFilesMCGen[i];
	  TString fileName_nominal = i_TreeDir + "/Nominal/" + chDir + "/" + chDir + "_" + gVFilesMCGen_nominal[i];
	  double w = CalcLumiWeight(fileName,i_year) * rescaleTmp;
	  double w_nominal = CalcLumiWeight(fileName_nominal,i_year);
	  printf("Applying year-to-year correlation rescale factor: %f \n",year_to_year_Corr_rescaleTEMP);
	  w *= year_to_year_Corr_rescaleTEMP;

	  if (ttmd_apply_trueLevelWeights && apply_tlrw_to_i_year) w *= GetTrueLevelRenormalisationWeight(fileName, var_without_yearCorrLabel);

	  tree->Add(fileName, w, doWeightVarTEMP);
	  //make sure that the correlation rescale factor is only applied to the variation itself and not to nominal
	  if(year_to_year_Corr_rescaleTEMP != 1.) {
	    double w_nominal_rescaled = w_nominal-year_to_year_Corr_rescaleTEMP*w_nominal;
	    tree->Add(fileName_nominal, w_nominal_rescaled, false);
	  }
	}
      }
    }

    tree->InitReading();

    return tree;
}
