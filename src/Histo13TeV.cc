#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TMath.h>

#include <vector>
#include <iostream>
#include <set>

#include "plotterclass13TeV.h"
#include "HistoListReader.h"
#include "UsefulTools13TeV.h"
#include "../../common/include/CommandLineParameters.h"
#include "../../common/include/sampleHelpers.h"

using namespace std;

void Histo13TeV(bool doControlPlots, bool doUnfold, bool doDiffXSPlotOnly,
           std::vector<std::string> plots, 
           std::vector<std::string> systematics, 
           std::vector<std::string> channels,
           const bool drawUncBand,
           const bool mergeEnvelopeVariations,
           const std::string& unfoldingType = "svd", // OZ 4.01.2017 for TUnfold (default SVD, if unfoldingType = "tunfold" then TUfold will be used)
           const bool useTopMassSetup = false
           )
{

    if (useTopMassSetup) std::cout<<"Using setup for top mass extraction analysis"<<std::endl;
    else std::cout<<"Using setup for differential cross-section analysis"<<std::endl;

    HistoListReader histoList(doControlPlots ? "HistoList_control_13TeV" : "HistoList_13TeV");
    if (histoList.IsZombie()) exit(12);
    for (auto it = histoList.begin(); it != histoList.end(); ++it) {
        const PlotProperties& p = it->second;
        std::cout << "checking " << p.name << std::endl;
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

        // Create Plotter13TeV 
        Plotter13TeV h_generalPlot;
        
        /////////////////////////////////////////////////////
        /////////   UNFOLDING OPTIONS     ///////////////////
        /////////////////////////////////////////////////////

        TString outpath = "";
        h_generalPlot.UnfoldingOptions(doUnfold);
        h_generalPlot.SetOutpath("");

        /// Do fit in the ratio plot, if ratio plot done
        h_generalPlot.DoFitInRatio(0);

        // Communicate use of CP plotting only
        if(doControlPlots) h_generalPlot.OnlyCPSetup();

        // Setup for top mass extraction analysis
        if (useTopMassSetup) h_generalPlot.SetTopMassSetup();

        /////////////////////////////////////////////////////
        /////////////////////////////////////////////////////
        ///////////////////////////////////////////////////// 
        h_generalPlot.setOptions(p.name,p.specialComment,p.ytitle,p.xtitle, 
                                    p.rebin, p.do_dyscale, p.logX, p.logY, 
                                    p.ymin, p.ymax, p.xmin, p.xmax, p.bins, p.xbinbounds, p.bincenters);
        h_generalPlot.DYScaleFactor(p.specialComment);
        //need preunfolding for ALL channels before unfolding!!
        for (auto channel : channels) {
            for (auto systematic : systematics) {
                h_generalPlot.setDrawUncBand(drawUncBand);
                h_generalPlot.setMergeEnvelopeVariations(mergeEnvelopeVariations);
                if (!doDiffXSPlotOnly) h_generalPlot.preunfolding(channel, systematic, unfoldingType);
            }
        }
        if (doUnfold) {
            for (auto systematic:systematics){
                for (auto channel:channels) {
                    h_generalPlot.unfolding(channel, systematic, unfoldingType);
                }
            }
        std::cout << "Done with the unfolding\n";
        }
        if (doDiffXSPlotOnly) {
            for (auto channel:channels){
                std::vector<TString> syst (systematics.begin(), systematics.end());
                std::vector<TString> valid_sys;
                for (size_t sys=0; sys<syst.size(); ++sys){
                    ///Ugly method to use the systematic convention used up to now
                    if(syst.at(sys) == TString("Nominal") || syst.at(sys).Contains("_DOWN") || syst.at(sys).Contains("HERWIG") || syst.at(sys).Contains("PERUGIA11") || syst.at(sys).Contains("SPINCORR") ){continue;}

                    if (!syst.at(sys).Contains("POWHEG") && !syst.at(sys).Contains("MCATNLO") && !syst.at(sys).Contains("BFRAG_PETERSON") &&
                        !syst.at(sys).Contains("ERDON") && !syst.at(sys).Contains("ERDONRETUNE") && !syst.at(sys).Contains("GLUONMOVETUNE") ) {valid_sys.push_back(syst.at(sys).Remove(syst.at(sys).Length()-2, 2));}
                    else {valid_sys.push_back(syst.at(sys));};
                };
                h_generalPlot.PlotDiffXSec(channel, valid_sys);
            }
        }
    }
}

/**
 * Helper function to create a function which checks if a string found is in the
 * passed vector of string.
 * 
 * @param allowed a vector of allowed strings (char*s)
 * @return a function taking a std::string and returning a bool
 */
std::function<bool(const std::string &s)> makeStringChecker(const std::vector<const char *> allowed) {
    return [allowed](const std::string &test) {
        return std::find(begin(allowed), end(allowed), test) != end(allowed);
    };
}

bool areSystematicsConsistent(const std::vector<std::string> systematics,const std::vector<const char *> allowed){
    for (auto syst : systematics)
        if (std::find(begin(allowed), end(allowed), syst) == end(allowed)) return false;
    return true;
}

int main(int argc, char** argv) {

    std::vector<const char*> VectorOfValidSystematics;
    UsefulTools13TeV::fillVectorOfValidSystematics(VectorOfValidSystematics,true);

    CLParameter<std::string> opt_type("t", "cp|unfold|plot - required, cp=contol plots, unfold, or only plot diffXS", true, 1, 1,
        makeStringChecker({"cp", "unfold", "plot"}));
    CLParameter<std::string> opt_plots("p", "Name (pattern) of plot; multiple patterns possible; use '+Name' to match name exactly", false, 1, 100);
    CLParameter<std::string> opt_channel("c", "Specify channel(s), valid: emu, ee, mumu, combined. Default: all channels", false, 1, 4,
        makeStringChecker({"ee", "emu", "mumu", "combined"}));
    CLParameter<std::string> opt_sys("s", "Systematic variation - default is Nominal, use 'all' for all", false, 1, 250,
        makeStringChecker(VectorOfValidSystematics));
    CLParameter<bool> opt_band("b", "If existing, draw uncertainty band if '-t cp' ", false, 0, 0);
    CLParameter<bool> opt_merge("m", "Merge scale variations if '-t cp' ", false, 0, 0);
    CLParameter<bool> opt_mass("topmass", "Use setup for top mass extraction analysis ", false, 0, 0);

    // OZ 4.01.2017 for TUnfold
    CLParameter<std::string> opt_unfolding("u", "Specify unfolding, valid: svd, tunfold. Default: svd", false, 1, 1,
        makeStringChecker({"svd", "tunfold"}));

    CLAnalyser::interpretGlobal(argc, argv);

    UsefulTools13TeV::fillVectorOfValidSystematics(VectorOfValidSystematics,opt_mass.isSet());
    if( !areSystematicsConsistent(opt_sys.getArguments(),VectorOfValidSystematics) ){
        std::cerr << "ERROR: selected systematic(s) not allowed. " << std::endl;
        exit(1);
    }

    // OZ 4.01.2017 for TUnfold
    std::string unfoldingType = "svd";
    if (opt_unfolding.isSet()) unfoldingType = opt_unfolding.getArguments()[0];
    //printf("OZ unfoldingType: %s\n", unfoldingType.c_str());
    
    std::vector<std::string> channels { "emu", "ee", "mumu", "combined" };
    if (opt_channel.isSet()) channels = opt_channel.getArguments();
    std::cout << "Processing channels: "; 
    for (auto ch: channels) std::cout << ch << " "; std::cout << "\n";
        
    std::vector<std::string> systematics { "Nominal" };
    if (opt_sys.isSet()) {
        systematics = opt_sys.getArguments();
        if (systematics[0] == "all") {
            systematics.clear();
            for (string syst: VectorOfValidSystematics) {
                if (syst != "all" && syst != "POWHEGHERWIG") systematics.push_back(syst);
            }
        }
    }
    std::cout << "Processing systematics (use >>-s all<< to process all knwon systematics): "; 
    for (auto sys: systematics) std::cout << sys << " "; std::cout << "\n";

    std::vector<std::string> plots { "" };
    if (opt_plots.isSet()) plots = opt_plots.getArguments();

    bool doControlPlots = opt_type[0] == "cp";
    bool doDiffXSPlotOnly = opt_type[0] == "plot";
    bool doUnfold = opt_type[0] == "unfold";
    bool drawUncBand = (opt_band.isSet() && doControlPlots);
    bool mergeEnvelopeVariations = (opt_merge.isSet() && ((doControlPlots && !drawUncBand) || doDiffXSPlotOnly));

    if (mergeEnvelopeVariations && doControlPlots){
        systematics.clear();
        systematics.push_back("TOT_SCALE_UP");
        systematics.push_back("TOT_SCALE_DOWN");
        systematics.push_back("TOT_BFRAG_UP");
        systematics.push_back("TOT_BFRAG_DOWN");
        systematics.push_back("TOT_COLORREC_UP");
        systematics.push_back("TOT_COLORREC_DOWN");
    }

    Histo13TeV(doControlPlots, doUnfold, doDiffXSPlotOnly, plots, systematics, channels, drawUncBand, mergeEnvelopeVariations ,unfoldingType, opt_mass.isSet());
}
