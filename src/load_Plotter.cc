#include <vector>
#include <iostream>
#include <set>
#include <cassert>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>

#include "../../common/include/CommandLineParameters.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/utils.h"

#include "../../common/include/RootFileReader.h"

#include "PlotterCP.h"
#include "MiniTreeAnalyzer.h"
#include "MiniTreeAppender.h"
#include "HistoListReader.h"
#include "PlotterConfigurationHelper.h"
// #include "sampleHelpers.h"
// #include "utils.h"

// namespace Systematic{


#include "ttmd_declareVars.h"
#include "ttmd_declareVars1D.h"
#include "ttmd_declareVars1D_particle.h"
#include "ttmd_declareVars2D.h"
#include "ttmd_declareVars3D.h"
#include "ttmd_declareVarsCP.h"
#include "PlotterDiffXSec.h"


using namespace std;

void load_Plotter(bool doControlPlots,
           std::vector<std::string> plots,
           std::string year,
           std::string yearCorr,
           std::vector<std::string> systematics,
           std::vector<std::string> channels,
           const bool drawUncBand,
           const bool mergeEnvelopeVariations,
           const bool useEnvelopesInSystematicsList
           )
{


  for (auto channel : channels) {
    std::string histoCPList = std::string("data/binnings/HistoList_CP_") + year + std::string("_") + channel;
    HistoListReader histoList(doControlPlots ? histoCPList.c_str() : "HistoList_13TeV"); //FIXME: doControlPlots always true, syntax kept to facilitate addition of new files for readout
    if (histoList.IsZombie()) exit(12);
    for (auto it = histoList.begin(); it != histoList.end(); ++it) {
        const PlotProperties& p = it->second;
        std::cout << std::endl << "--- checking " << p.name << std::endl;
        bool found = false;
        for (auto plot : plots) {
            if (plot.size() && plot[0] == '+') {
                if (p.name.CompareTo(&plot[1], TString::kIgnoreCase) == 0) {
                    found = true;
                    break;
                }
            } else if (p.name.Contains(plot, TString::kIgnoreCase)) {
                found = true;
                break;
            }
        }
        if (!found) continue;

        // Create PlotterCP and configure
        PlotterCP h_generalPlot(year, yearCorr);

        TString outpath = "";
        h_generalPlot.setOutpath(outpath);

        // Do fit in the ratio plot, if ratio plot done
        h_generalPlot.DoFitInRatio(0);

        // Communicate use of CP plotting only
        if (doControlPlots) h_generalPlot.OnlyCPSetup();

        h_generalPlot.setUseEnvelopesInSystematicsList(useEnvelopesInSystematicsList);

        h_generalPlot.setOptions(p.name, p.specialComment, p.ytitle, p.xtitle,
                                 p.rebin, p.do_dyscale, p.logX, p.logY, p.ymin, p.ymax, p.xmin, p.xmax,
                                 p.bins, p.xbinbounds, p.bincenters);

        h_generalPlot.DYScaleFactor(p.specialComment);


            for (auto systematic : systematics) {
                h_generalPlot.setDrawUncBand(drawUncBand);
                h_generalPlot.setMergeEnvelopeVariations(mergeEnvelopeVariations);
                h_generalPlot.plotWriterWrapper(channel, systematic);
            }

    }

  }
}

void load_ttmd(std::string year, std::string yearCorr,
           int modeKR,
           int coverIter,
           int coverMode,
           bool flagProcSyst,
           std::vector<TString> unfVars,
           std::vector<TString> chs,
	   TString ttmdDim,
           TString ttmdType
           )
{//FIXME make it prettier and more adaptable once we have ported the rest
    // 0: fr, 1: oz
    int mode = 1;
    bool applyDYScale = true;
    RootFileReader* fileReader;
    fileReader = RootFileReader::getInstance();
    PlotterConfigurationHelper* plotterConfigHelper = new PlotterConfigurationHelper(fileReader, year, yearCorr, applyDYScale, 1, flagProcSyst);

    plotterConfigHelper->DefaultSetup();
    // optional analysis suffix for plots directory
    TString analysisSuffix = plotterConfigHelper->analysisSuffix; // configuration for the initial setup
    // std::cout<<plotterConfigHelper->gFlagMG<<std::endl;
    PlotterDiffXSec* an = new PlotterDiffXSec(year,plotterConfigHelper, mode, analysisSuffix, modeKR);
    //an->MaxEvents = 50000;

    //bool tmpAll = 0;
    // bool flagCP       = 0;
    // bool flagCP25     = 0 || tmpAll;
    // bool flagPSE      = 0 || tmpAll;
    // bool flagXSec82D  = 0 || tmpAll;
    // bool flagXSec81D  = 0 || tmpAll;
    // bool flagXSec1D   = 0 || tmpAll;
    // bool flagXSec3D   = 0 || tmpAll;
    bool flag171018ExtraPlots = 1; // 25.10.18

    //an->ModeKinRecoTree = 0;
    //if(modeKR >= 0)
    //  an->ModeKinRecoTree = modeKR;

    an->vCh = chs;

    std::vector<TString> vVar;
    //std::vector<TString> vVar = { "JESFragmentation_UP", "JESFragmentation_DOWN" };
    for(auto& var : vVar)
    {
      an->VVarStoreCP.push_back(var);
      an->VVarStorePSE.push_back(var);
      an->VVarStoreXSec.push_back(var);
    }


    if(flag171018ExtraPlots)
    {

        // cp
        if (ttmdDim == "all" || ttmdDim == "cp"){
            if(an->ModeKinRecoTree == 0) an->AddXSec(DeclareXSecAllCP(plotterConfigHelper));
            else if(an->ModeKinRecoTree == 9){
                //an->AddXSec(DeclareXSecMttPSE25(an->ModeKinRecoTree));
                //an->AddXSec(DeclareXSecYttPSE25(an->ModeKinRecoTree));
                an->AddXSec(DeclareXSecNjCP(plotterConfigHelper,"def", 0.4, 30.0, 4));
                //an->AddXSec(DeclareXSecNjCP(plotterConfigHelper,"def", 0.4, 50.0, 4));
                //an->AddXSec(DeclareXSecNjCP(plotterConfigHelper,"def", 0.4, 100.0, 4));
            }
        }

        // 1d
        if (ttmdDim == "all" || ttmdDim == "1d"){
            if(an->ModeKinRecoTree == 0){

	      if (ttmdType == "all" || ttmdType == "parton"){
                an->AddXSec(DeclareXSecPtt1D(plotterConfigHelper));
                an->AddXSec(DeclareXSecYt1D(plotterConfigHelper));
                an->AddXSec(DeclareXSecPtat1D(plotterConfigHelper));
                an->AddXSec(DeclareXSecYat1D(plotterConfigHelper));
                an->AddXSec(DeclareXSecYtt1D(plotterConfigHelper));
                an->AddXSec(DeclareXSecMtt1D(plotterConfigHelper));
                an->AddXSec(DeclareXSecPttt1D_new(plotterConfigHelper));
                an->AddXSec(DeclareXSecDphitt1D_new(plotterConfigHelper));
                an->AddXSec(DeclareXSecDytt1D(plotterConfigHelper));
                an->AddXSec(DeclareXSecRatioPttMtt(plotterConfigHelper));
                an->AddXSec(DeclareXSecRatioPtttMtt(plotterConfigHelper));
		an->AddXSec(DeclareXSecRecoPartonMomFraction(plotterConfigHelper));
                an->AddXSec(DeclareXSecRecoAntipartonMomFraction(plotterConfigHelper));

		//an->AddXSec(DeclareXSecNj(plotterConfigHelper,"def", 0.4, 30.0, 2));
                //an->AddXSec(DeclareXSecNj(plotterConfigHelper,"def", 0.4, 30.0, 3));
                //an->AddXSec(DeclareXSecYtLead1D(plotterConfigHelper));
                //an->AddXSec(DeclareXSecYtNLead1D(plotterConfigHelper));
                //an->AddXSec(DeclareXSecPttLead1D(plotterConfigHelper));
                //an->AddXSec(DeclareXSecPttNLead1D(plotterConfigHelper));
                //an->AddXSec(DeclareXSecPttTTRestFrame1D(plotterConfigHelper));
                //an->AddXSec(DeclareXSecPtatTTRestFramE1d(plotterConfigHelper));
	      }
	      else if (ttmdType == "all" || ttmdType == "particle"){
		an->AddXSec(DeclareXSecPtt1D_particle(plotterConfigHelper));
		an->AddXSec(DeclareXSecYt1D_particle(plotterConfigHelper));
		an->AddXSec(DeclareXSecPtat1D_particle(plotterConfigHelper));
		an->AddXSec(DeclareXSecYat1D_particle(plotterConfigHelper));
		an->AddXSec(DeclareXSecYtt1DParticle(plotterConfigHelper));
		an->AddXSec(DeclareXSecMtt1DParticle(plotterConfigHelper));
		an->AddXSec(DeclareXSecPttt1D_new_particle(plotterConfigHelper));
		an->AddXSec(DeclareXSecDphitt1D_new_particle(plotterConfigHelper));
		an->AddXSec(DeclareXSecDytt1D_particle(plotterConfigHelper));
                an->AddXSec(DeclareXSecRatioPttMtt_particle(plotterConfigHelper));
                an->AddXSec(DeclareXSecRatioPtttMtt_particle(plotterConfigHelper));
                an->AddXSec(DeclareXSecRecoPartonMomFraction_particle(plotterConfigHelper));
                an->AddXSec(DeclareXSecRecoAntipartonMomFraction_particle(plotterConfigHelper));

		an->AddXSec(DeclareXSecPtbLead1D_particle(plotterConfigHelper));
		an->AddXSec(DeclareXSecPtbNLead1D_particle(plotterConfigHelper));
		an->AddXSec(DeclareXSecPtl1D_particle(plotterConfigHelper));
		an->AddXSec(DeclareXSecRatioPtbLeadt_particle(plotterConfigHelper));
		an->AddXSec(DeclareXSecRatioPtbNLeadt_particle(plotterConfigHelper));
		an->AddXSec(DeclareXSecRatioPtlt_particle(plotterConfigHelper));
		an->AddXSec(DeclareXSecMllbbmet1D_particle(plotterConfigHelper));
		an->AddXSec(DeclareXSecMllbb1D_particle(plotterConfigHelper));
		an->AddXSec(DeclareXSecMbb1D_particle(plotterConfigHelper));
		an->AddXSec(DeclareXSecMll1D_particle(plotterConfigHelper));

		//an->AddXSec(DeclareXSecYtLead1D_particle(plotterConfigHelper));
		//an->AddXSec(DeclareXSecYtNLead1D_particle(plotterConfigHelper));
		//an->AddXSec(DeclareXSecPttLead1D_particle(plotterConfigHelper));
		//an->AddXSec(DeclareXSecPttNLead1D_particle(plotterConfigHelper));
		//an->AddXSec(DeclareXSecPttTTRestFrame1D_particle(plotterConfigHelper));
		//an->AddXSec(DeclareXSecPtatTTRestFrame1D_particle(plotterConfigHelper));
		//an->AddXSec(DeclareXSecPtal1D_particle(plotterConfigHelper));
		//an->AddXSec(DeclareXSecPtlLead1D_particle(plotterConfigHelper));
		//an->AddXSec(DeclareXSecPtlNLead1D_particle(plotterConfigHelper));
		//an->AddXSec(DeclareXSecEtal1D_particle(plotterConfigHelper));
		//an->AddXSec(DeclareXSecEtaal1D_particle(plotterConfigHelper));
		//an->AddXSec(DeclareXSecEtalLead1D_particle(plotterConfigHelper));
		//an->AddXSec(DeclareXSecEtalNLead1D_particle(plotterConfigHelper));
		//an->AddXSec(DeclareXSecDphill1D_new_particle(plotterConfigHelper));
		//an->AddXSec(DeclareXSecDEtall1D_particle(plotterConfigHelper));
		//an->AddXSec(DeclareXSecPtll1D_new_particle(plotterConfigHelper));
		//an->AddXSec(DeclareXSecEtabLead1D_particle(plotterConfigHelper));
		//an->AddXSec(DeclareXSecEtabNLead1D_particle(plotterConfigHelper));
		//an->AddXSec(DeclareXSecPtbb1D_new_particle(plotterConfigHelper));
	      }

            }
            else if(an->ModeKinRecoTree == 9){
	      if (ttmdType == "all" || ttmdType == "parton"){
                an->AddXSec(DeclareXSecYtt1D(plotterConfigHelper, an->ModeKinRecoTree));
                an->AddXSec(DeclareXSecMtt1D(plotterConfigHelper, an->ModeKinRecoTree));
                an->AddXSec(DeclareXSecPttt1D_new(plotterConfigHelper, an->ModeKinRecoTree));
                an->AddXSec(DeclareXSecMtt1D_threshold(plotterConfigHelper, an->ModeKinRecoTree));
	      }
	      else if (ttmdType == "all" || ttmdType == "particle"){
		an->AddXSec(DeclareXSecYtt1DParticle(plotterConfigHelper, an->ModeKinRecoTree));
		an->AddXSec(DeclareXSecMtt1DParticle(plotterConfigHelper, an->ModeKinRecoTree));
		an->AddXSec(DeclareXSecPttt1D_new_particle(plotterConfigHelper, an->ModeKinRecoTree));
	      }
            }

        }
        // 1d extra jets
        if (ttmdDim == "all" || ttmdDim == "nj" || ttmdDim == "multidim"){
            if(an->ModeKinRecoTree == 0){
              an->AddXSec(DeclareXSecNj_simple(plotterConfigHelper,"def", 0.4, 30.0, 4));
              // an->AddXSec(DeclareXSecNj_simple(plotterConfigHelper,"def", 0.4, 40.0, 4));
              an->AddXSec(DeclareXSecNj_simple(plotterConfigHelper,"def", 0.4, 50.0, 4));
              // an->AddXSec(DeclareXSecNj_simple(plotterConfigHelper,"def", 0.4, 75.0, 4));
              an->AddXSec(DeclareXSecNj_simple(plotterConfigHelper,"def", 0.4, 100.0, 3)); // nj reduced to 3 due to bkg subtraction issue in 2017 (ee channel)
              // an->AddXSec(DeclareXSecNj_simple(plotterConfigHelper,"def", 0.4, 150.0, 3));
            }
        }


        // 2d
        const int binning_option = 0; // 1->>TOP-18-004, 2->>Detc. Levl. bin number reduced
        if (ttmdDim == "all" || ttmdDim == "2d" || ttmdDim == "multidim"){
            // xsecs for full and loose KR
            an->AddXSec(DeclareXSecNjMtt(plotterConfigHelper,"def", 0.4, 30.0, 3, an->ModeKinRecoTree));
            an->AddXSec(DeclareXSecNjYtt(plotterConfigHelper,"def", 0.4, 30.0, 3, an->ModeKinRecoTree));
            an->AddXSec(DeclareXSecNjPttt(plotterConfigHelper,"def", 0.4, 30.0, 3, an->ModeKinRecoTree));
            an->AddXSec(DeclareXSecYttPtttXSec8(plotterConfigHelper, binning_option, an->ModeKinRecoTree));
            an->AddXSec(DeclareXSecMttYttXSec8(plotterConfigHelper, binning_option, an->ModeKinRecoTree));
            an->AddXSec(DeclareXSecMttPtttXSec8(plotterConfigHelper, binning_option, an->ModeKinRecoTree));

            // xsecs only possible with full KR
            if(an->ModeKinRecoTree == 0){
                an->AddXSec(DeclareXSecMttPttXSec8(plotterConfigHelper, binning_option));
                an->AddXSec(DeclareXSecYtPttXSec8(plotterConfigHelper, binning_option));
                an->AddXSec(DeclareXSecMttDetattXSec8(plotterConfigHelper, binning_option));
                an->AddXSec(DeclareXSecMttYtXSec8(plotterConfigHelper, binning_option));
                an->AddXSec(DeclareXSecMttDphittXSec8(plotterConfigHelper, binning_option));
                an->AddXSec(DeclareXSecPttPtttXSec8(plotterConfigHelper, binning_option));

                an->AddXSec(DeclareXSecNjDetatt(plotterConfigHelper,"def", 0.4, 30.0, 3));
                an->AddXSec(DeclareXSecNjPtt(plotterConfigHelper,"def", 0.4, 30.0, 3));
                an->AddXSec(DeclareXSecNjYt(plotterConfigHelper,"def", 0.4, 30.0, 3));

            }
        }

        // 3d
        if (ttmdDim == "all" || ttmdDim == "3d" || ttmdDim == "multidim"){
            // xsecs for full and loose KR
            an->AddXSec(DeclareXSecNjMttYtt(plotterConfigHelper,"def", 0.4, 30.0, 2, 3, 4, an->ModeKinRecoTree)); //TOP-18-004
            an->AddXSec(DeclareXSecNjMttYtt(plotterConfigHelper,"def", 0.4, 30.0, 3, 3, 4, an->ModeKinRecoTree)); //TOP-18-004
            // an->AddXSec(DeclareXSecNjMttPttt(plotterConfigHelper,"def", 0.4, 30.0, 3, 3, 3, an->ModeKinRecoTree));
            an->AddXSec(DeclareXSecPtttMttYtt(plotterConfigHelper, an->ModeKinRecoTree));

            if(an->ModeKinRecoTree == 0){
                // an->AddXSec(DeclareXSecNjMttDeta(plotterConfigHelper,"def", 0.4, 30.0, 3, 3, 3, an->ModeKinRecoTree));
                an->AddXSec(DeclareXSecPttMttPttt(plotterConfigHelper));
            }
            else if(an->ModeKinRecoTree == 9){
                an->AddXSec(DeclareXSecPttt2BMttYtt(plotterConfigHelper,50.0, an->ModeKinRecoTree));
            }

            if (year == "fullRun2")
                an->AddXSec(DeclareXSecNjMttYtt(plotterConfigHelper,"def", 0.4, 30.0, 4, 3, 4, an->ModeKinRecoTree)); //Not available for single years. Only working with fullRun2!
            }
    }

    if(coverMode >= 0){
      for(auto& ch : an->vCh) an->CoverTest(ch, coverMode, coverIter);
    }

    // unfolding
    std::vector<TString> vars;
    if(unfVars.size() == 1 && unfVars[0] == "all") {
        vars = plotterConfigHelper->GetAllVar(true);
        std::cout << " **** All systematics will be processed: ";
        for (TString sys: vars) std::cout << sys << " ";
        std::cout << "\n \n";
    }
    else if (unfVars.size() == 0) vars.clear();
    else
      vars.insert(vars.end(), unfVars.begin(), unfVars.end());

    for(auto var : vars) {
      std::cout << " **** Processing: " << var << " ****"<< std::endl << std::endl;
      an->SetVar(var);
      for(auto& ch : an->vCh) an->Analyse(ch);
    }

    if(flagProcSyst){
      an->ProcessSystematics();
    }

}

void load_MiniTreeAnalysis(
           TString year,
           TString filename,
           TString systematic_string,
           TString channel_string,
           bool doDNNAppend
           )
{
    if(!doDNNAppend){
        MiniTreeAnalyzer l_Analyzer(year, filename, systematic_string, channel_string);
        l_Analyzer.Analyze();
        l_Analyzer.WriteOutput();
    }else{
        MiniTreeAppender l_Appender(year, filename, systematic_string, channel_string);
        l_Appender.Append();
    }
}

/**
 * Helper function to create a function which checks if a string found is in the
 * passed vector of string.
 *
 * @param allowed a vector of allowed strings (char*s)
 * @return a function taking a std::string and returning a bool
 */
std::function<bool(const std::string &s)> makeStringChecker(const std::vector<const char *> allowed)
{
    return [allowed](const std::string &test) {
        return std::find(begin(allowed), end(allowed), test) != end(allowed);
    };
}

bool areSystematicsConsistent(const std::vector<std::string> systematics,const std::vector<const char *> allowed)
{
    for (auto syst : systematics)
        if (std::find(begin(allowed), end(allowed), syst) == end(allowed)) return false;
    return true;
}
bool areSystematicsConsistent(const std::vector<std::string> systematics,const std::vector<TString> allowed)
{
    for (auto syst : systematics)
        if (std::find(begin(allowed), end(allowed), syst) == end(allowed)) return false;
    return true;
}

int main(int argc, char** argv)
{
    std::vector<TString> VectorOfValidSystematics;

    CLParameter<std::string> opt_type("t", "cp=contol plots, ttmd=ttmd unfolding, mt=MiniTree analysis", true, 1, 1, makeStringChecker({"cp","ttmd","mt"}));
    CLParameter<std::string> opt_year("y", "Specify year, valid: 2016, 2017, 2018, fullRun2UL", true, 1, 1, makeStringChecker({"2016", "2016preVFP", "2016postVFP", "2017", "2018", "fullRun2UL"}));
    CLParameter<std::string> opt_plots("p", "Name (pattern) of plot; multiple patterns possible; use '+Name' to match name exactly ", false, 1, 100);
    CLParameter<std::string> opt_channel("c", "Specify channel(s), valid: emu, ee, mumu, combined. Default: all channels ", false, 1, 4, makeStringChecker({"ee", "emu", "mumu", "combined"}));
    CLParameter<std::string> opt_sys("s", "Systematic variation - default is Nominal, use 'all' for all ", false, 1, 250);
    CLParameter<bool> opt_band("b", "If existing, draw uncertainty band if '-t cp' ", false, 0, 0);
    CLParameter<bool> opt_merge("m", "Merge scale variations if '-t cp' ", false, 0, 0);
    CLParameter<bool> opt_moresyst("moresyst", "Use setup with extended list of systematic variations ", false, 0, 0);
    //FIXME descriptions and types
    CLParameter<std::string> opt_ttmd_kinReco("kr", "Type of kinematic reconstruction for ttmd mode, default is -1; options are: 9 ", false, 1, 1, makeStringChecker({"9", "0"}));
    CLParameter<std::string> opt_ttmd_coverIter("coverIter", "Cover Iter?, default is -1; options are: ???", false, 1, 1);
    CLParameter<std::string> opt_ttmd_coverMode("coverMode", "Cover Test Mode?, default is -1; options are: ???", false, 1, 1);
    CLParameter<bool> opt_ttmd_doSyst("ttmdSyst", "If existing, process systs if '-t ttmd' ", false, 0, 0);
    CLParameter<std::string> opt_ttmd_dim("ttmdDim", "Run only over 1d, 2d, 3d, multidim, cp, all - default is all", false, 1, 1, makeStringChecker({"1d", "2d", "3d", "nj", "multidim", "cp", "all"}));
    CLParameter<std::string> opt_ttmd_type("ttmdType", "Run over parton, particle or all - default is all", false, 1, 1, makeStringChecker({"parton", "particle", "all"}));
    CLParameter<std::string> opt_mt_sys("ls", "Systematic variation - default is Nominal, use 'all' for all  if '-t mt' ", false, 1, 100, common::makeStringCheckBegin(Systematic::convertType(Systematic::allowedSystematics)));
    CLParameter<std::string> opt_mt_ch("lc", "Specify channel, valid: emu, ee, mumu. Default: emu, if '-t mt' ", false, 1, 4, common::makeStringCheck(Channel::convert(Channel::allowedChannelsAnalysis)));
    CLParameter<std::string> opt_mt_file("lf", "Specify filename, e.g. ttbarsignalplustau.root. Necessary if '-t mt' ", false, 1, 250);
    CLParameter<bool> opt_mt_append("lAppend", "Do miniTreeAppending of DNN score if '-t mt' ", false, 0, 0);
    CLParameter<std::string> opt_doCorrYear("doCorrYear", "treat year to year correlations for systematics when running fullRun2. Specify doCorrYear, valid: 2016, 2017, 2018, only16, only17, only18, only161718, all", false, 1, 1, makeStringChecker({"2016", "2017", "2018", "only16", "only17", "only18", "only161718", "all"}));

    CLAnalyser::interpretGlobal(argc, argv);

    std::string year; //FIXME: add functionality to handle all years, i.e. RunII, in a single plot
    if (opt_year.isSet()) year = opt_year.getArguments()[0];
    std::cout << "Processing setup for year: " << year << std::endl;

    bool doUnfolding=false;
    bool doControlPlots=false;
    bool doMiniTreeAnalysis=false;
    if (opt_type.isSet()){
      doMiniTreeAnalysis = opt_type[0] == "mt";
      doControlPlots = opt_type[0] == "cp";
      doUnfolding = opt_type[0] == "ttmd";
    }
    bool drawUncBand = (opt_band.isSet() && doControlPlots);

    std::string yearCorr;
    if ((year == "fullRun2UL" && opt_sys.isSet() && opt_doCorrYear.isSet()) || (drawUncBand && year == "fullRun2UL")) {
      if (opt_doCorrYear.isSet()) yearCorr = opt_doCorrYear.getArguments()[0];
      for (auto sys: opt_sys.getArguments()) {
	if (PlotterConfigurationHelper::IsYearToYearCorr(year, yearCorr, sys) || opt_sys.getArguments()[0] == "all")
	  std::cout << "Treating systematic " << sys << " " << "as uncorrelated or partially correlated for " << yearCorr << std::endl;
      }
    }
    else if (year == "fullRun2UL" && opt_ttmd_doSyst.isSet() && opt_doCorrYear.isSet()) {
	yearCorr = opt_doCorrYear.getArguments()[0];
    }
    else if (year == "fullRun2UL" && opt_sys.isSet() && !opt_doCorrYear.isSet()) {
      for (auto sys: opt_sys.getArguments()) {
	if (PlotterConfigurationHelper::IsYearToYearCorr(year, yearCorr, sys) || sys == "all") {
	  //should be changed to exit under these circumstances as soon as year to year corr handling for cp is implemented
	  if (opt_sys.getArguments()[0] == "all")
	    std::cerr << "WARNING: " << "some systematics are (un)correlated with the other years. Please use option --doCorrYear" << std::endl;
	  else
	    std::cerr << "WARNING: " << sys << " is (un)correlated with the other years. Please use option --doCorrYear" << std::endl;
	}
      }
    }
    else if (year != "fullRun2UL" && opt_doCorrYear.isSet()) {
      std::cerr << "ERROR: please only use option --doCorrYear with fullRun2UL" << std::endl;
      exit(1);
    }
    else if (!opt_sys.isSet() && !opt_ttmd_doSyst.isSet() && opt_doCorrYear.isSet() && !drawUncBand) {
      std::cerr << "ERROR: please specify a systematic using the -s option which you want to treat as uncorrelated or partially correlated with the other years " << std::endl;
      exit(1);
    }

    std::vector<std::string> channels {"emu", "ee", "mumu", "combined"};
    if (opt_channel.isSet()) channels = opt_channel.getArguments();
    if(doControlPlots||doUnfolding){
        std::cout << "Processing channels: ";
        for (auto ch: channels) std::cout << ch << " ";
        std::cout << "\n";
    }

    PlotterConfigurationHelper::fillVectorOfValidSystematics(VectorOfValidSystematics, year, yearCorr, doUnfolding, opt_moresyst.isSet());

    std::vector<std::string> systematics;
    if (opt_sys.isSet()) {
      for (auto sys: opt_sys.getArguments()) {
	std::vector<TString> yearCorrList = PlotterConfigurationHelper::GetYearCorrList(year, yearCorr, sys, PlotterConfigurationHelper::IsYearToYearCorr(year, yearCorr, sys));
	if (yearCorrList.size() != 0)
	  systematics.insert(systematics.end(), yearCorrList.begin(), yearCorrList.end());
	else
	  systematics.push_back(sys);
      }
      if (!areSystematicsConsistent(systematics, VectorOfValidSystematics)) {
	std::cerr << "ERROR: selected systematic(s) not allowed: " << std::endl;
	for (auto sys: systematics)
	  std::cerr << sys << std::endl;
	exit(1);
      }
    }
    
    if(doControlPlots){
        PlotterConfigurationHelper *dummy_config = new PlotterConfigurationHelper(year, yearCorr);
        std::vector<std::string> systematics {"Nominal"};
        if (opt_sys.isSet()) {
	  if (opt_sys.getArguments()[0] == "all") {
                systematics.clear();
                std::vector<TString> allSysExp = dummy_config->GetExpSys(true);
                std::vector<TString> allSysModW = dummy_config->GetModWSys(true);
                std::vector<TString> allSysModI = dummy_config->GetModISys(true);
                systematics.insert(systematics.end(), allSysExp.begin(), allSysExp.end());
                systematics.insert(systematics.end(), allSysModW.begin(), allSysModW.end());
                systematics.insert(systematics.end(), allSysModI.begin(), allSysModI.end());
                systematics.insert(systematics.begin(), "Nominal");
                /*for (TString syst: VectorOfValidSystematics) {
                    if (syst != "all" && syst != "POWHEGHERWIG") systematics.push_back(std::string(syst));
                }*/
	  }
	  else {
	    systematics.clear();
	    for (auto sys: opt_sys.getArguments()) {
	      std::vector<TString> yearCorrList = PlotterConfigurationHelper::GetYearCorrList(year, yearCorr, sys, PlotterConfigurationHelper::IsYearToYearCorr(year, yearCorr, sys));
	      if (yearCorrList.size() != 0)
		systematics.insert(systematics.end(), yearCorrList.begin(), yearCorrList.end());
	      else
		systematics.push_back(sys);
	    }
	  }
        }

        std::vector<std::string> plots {""};
        if (opt_plots.isSet()) plots = opt_plots.getArguments();
        bool mergeEnvelopeVariations = (opt_merge.isSet() && ((doControlPlots && !drawUncBand))); // FIXME: syntax with double brackets kept to facilitate addition of new options
        bool useEnvelopesInSystematicsList = (opt_merge.isSet() && drawUncBand);

        if (mergeEnvelopeVariations && doControlPlots) {
            systematics.clear();
            std::vector<TString> envelopes_variations = dummy_config->GetEnvelopesList(true);
            systematics.insert(systematics.begin(), envelopes_variations.begin(), envelopes_variations.end());
        }

        std::cout << "Processing systematics (use >>-s all<< to process all known systematics): ";
        for (std::string sys: systematics) std::cout << sys << " ";
        std::cout << "\n";

        load_Plotter(doControlPlots, plots, year, yearCorr, systematics, channels, drawUncBand, mergeEnvelopeVariations, useEnvelopesInSystematicsList);
    }

    // ====================================================
    // ======== Porting of ttmd_run_unfolding.cxx =========
    // ====================================================
    int argModeKR = -1;
    int argCoverIter = -1;
    int argCoverMode = -1;
    bool argFlagProcSyst = 0;
    TString ttmdDim = "all";
    TString ttmdType = "all";
    std::vector<TString> argUnfVar;

    if (opt_ttmd_kinReco.isSet()) argModeKR = stoi(opt_ttmd_kinReco.getArguments()[0]);
    if (opt_ttmd_coverIter.isSet()) argCoverIter = stoi(opt_ttmd_coverIter.getArguments()[0]);
    if (opt_ttmd_coverMode.isSet()) argCoverMode = stoi(opt_ttmd_coverMode.getArguments()[0]);
    if (opt_ttmd_dim.isSet()) ttmdDim = opt_ttmd_dim.getArguments()[0];
    if (opt_ttmd_type.isSet()) ttmdType = opt_ttmd_type.getArguments()[0];
    if (opt_ttmd_doSyst.isSet()) argFlagProcSyst = opt_ttmd_doSyst.isSet();
    //if (opt_ttmd_unfVar.isSet()) argUnfVar = (opt_ttmd_unfVar.getArguments()[0]);
    if (opt_sys.isSet()) argUnfVar.insert(argUnfVar.end(), systematics.begin(), systematics.end());
    else{
      //in case no opt_sys argument has been specified but opt_ttmd_doSyst is set don't default to Nominal
      if (argFlagProcSyst == 1) argUnfVar.clear();
      else argUnfVar.push_back("Nominal");
    }

    //FIXME decide what to use combined or ll?
    std::vector<TString> ttmdChs;
    for (auto ch: channels){
      if(ch=="combined"){
          ttmdChs.push_back("ll");
      }else{
          ttmdChs.push_back(ch);
      }
    }

    if (doUnfolding){
      if(argModeKR==-1){
        std::cerr << "ERROR: selected kinReco not allowed. " << std::endl;
        exit(1);
      }
      cout<<"Proccess ttmd Unfolding with following settings:"<<endl;
      cout<<"Kinematic reconstruction mode: "<<argModeKR<<endl;
      cout<<"Process xsecs with dimension(s): "<<ttmdDim<<endl;
      cout<<"Process xsecs with type(s): "<<ttmdType<<endl;
      cout<<"argCoverIter: "<<argCoverIter<<endl;
      cout<<"Cover Test mode: "<<argCoverMode<<endl;
      cout<<"Process systematics? : "<<argFlagProcSyst<<endl;
      for (auto unfVar: argUnfVar)
	cout<<"Unfolding variation(s): "<<unfVar<<endl;
      cout<<"channels: ";
      for (auto ch: ttmdChs){
        cout<<ch<<" ";
      }
      cout<<endl;

      load_ttmd(year, yearCorr, argModeKR, argCoverIter, argCoverMode, argFlagProcSyst, argUnfVar, ttmdChs, ttmdDim, ttmdType);
    }

    if(doMiniTreeAnalysis){
        TString filename("");
        if (opt_mt_file.isSet()){
            filename = opt_mt_file.getArguments()[0];
        }
        if (filename ==""){
            std::cerr << "ERROR: selected file is empty. " << std::endl;
            exit(1);
        }
        std::cout << "Processing file: ";
        std::cout << filename << " ";
        std::cout << "\n";
        TString channel("emu");
        if (opt_mt_ch.isSet()){
            channel = opt_mt_ch.getArguments()[0];
        }
        std::cout << "Processing channel: ";
        std::cout << channel << " ";
        std::cout << "\n";
        TString systematic("Nominal");
        if (opt_mt_sys.isSet()) {
            systematic = opt_mt_sys.getArguments()[0];
        }
        std::cout << "Processing systematic: ";
        std::cout << systematic << " ";
        std::cout << "\n";

        bool doAppend = opt_mt_append.isSet();
        load_MiniTreeAnalysis(year, filename, systematic, channel, doAppend);
    }
}
