#include <fstream>
#include <iostream>
#include <cstdio>
#include <sstream>
#include <cmath>
#include <memory>

#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TExec.h>
#include <TStyle.h>
#include <TMath.h>
#include <TROOT.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TList.h>
#include <TGraphAsymmErrors.h>
#include <TError.h>
#include <TRandom3.h>
#include<TMatrixD.h>

#include "DilepSVDFunctions.h"
#include "plotterclass13TeV.h"
#include "UsefulTools13TeV.h"
#include "ttbarUtils.h"
#include "../../common/include/RootFileReader.h"
#include "../../common/include/plotterUtils.h"

// OZ 4.01.2017 for TUnfold
#include "PlotterForTUnfold.h"
#include "Samples.h"

using namespace std;

const bool Plotter13TeV::doClosureTest = 0; //Signal is MC - so add BG to it and dont do DY scaling (please don't use with "usePoissonSmearedPseudoData == 1" unless you really know what you are doing)
const bool Plotter13TeV::usePoissonSmearedPseudoData = 0; //Use in case you need pseudo-data defined as: poisson smeared sum of all MCs (please look in README file)
const bool doClosureTestWithLowStat = 0;
const bool performUnfoldingInCombinedChannel = false; //Perform unfolding based on the input from combined channel (result will be stored in "combinedUnfolded", instead of calculating results from weighted mean of ee/emu/mumu)
const bool drawTotalMCErrorForCP = true; //Draw error band on control plot as obtained from total variations in all MCs (instead of standalone variations in ttbar)
const bool revokeAntiQuantityCombination = true; //Prevents construction of one (averaged, if applicable) histogram out of two objects : quantity and anti-quantity, if allowed for given quantities (e.g. top quark & top anti-quark)
const bool evaluateTotCovMatricesForSVDUnfolding = false; //Allow evaluation of tot. cov matrices (don't use blindly, it is systematic-type dependent). For now, configured ONLY for use with SVD in ee, emu, mumu, combined-UNFOLDED.

double BranchingFraction[4] = { 0.01147,0.01130,0.02277,0.04554 };//[ee, mumu, emu, combined] not including tau //FIXME: put to proper place
//double BranchingFraction[4]={0.01582, 0.01573, 0.03155, 0.06310};//[ee, mumu, emu, combined] including tau into electron/muon - update needed

void Plotter13TeV::setDrawUncBand(bool drawUncBand)
{
    drawUncBand_ = drawUncBand;
}

void Plotter13TeV::UnfoldingOptions(bool doSVD)
{
    doUnfolding = doSVD;
    drawNLOCurves = true; // boolean to draw/not-draw extra theory curves in the Diff.XSection plots
    absoluteXsec = false; // boolean to calculate/not calculate absolute xsec
    particleXsec = false; // boolean to unfold/not unfold results back to particle (aka pseudo top) level

    if (particleXsec) {
        visTheoryPrefix = "VisPseudo";
        respMatrixPrefix = "PseudoReco";
    }
    else {
        visTheoryPrefix = "VisGen";
        respMatrixPrefix = "GenReco";
    }

    drawPlotRatio = true;
    drawSmoothMadgraph = false;
    drawMCATNLO = true;
    drawKidonakis = false;
    drawAhrens = false;
    drawPOWHEG = true;
    drawPOWHEGHERWIG = false;
    drawPERUGIA11 = false;
    drawMadScaleMatching = false;
    drawMadMass = false; //FIXME: update needed to support functionality in absxsec: update GetNLOCurveMass in a similar way to GetNLOCurve, also propagate further treatment for rebinned histograms

    //FIXME: move out options below to the separate new method : PlotterOptions()

    // Setting to "true" triggers estimation of syst. unc. from 'varied pseudo-data', aka estimation of expected sys. uncertainty
    // (i.e. in case for each variation: using always nominal MC simulation for signal & bkgd estimations, while 'varied pseudo-data' is constructed out of varied MC simulations).
    // If needed, use together with 'usePoissonSmearedPseudoData == 1' and create special filelists beforehand via: scripts/mk_inputForVariedPseudoData_13TeV.sh and scripts/mk_HistoFileListForVariedPseudoData_13TeV.sh
    estimateSysUncFromPseudoData = 0;

    lumiError = 0.025;
    BRError = 0.015;
    massUncReductionFactor = 3.0; // divide resulting mass uncertainty by 'factor' needed to match +-1GeV variation in top-mass (e.g. here we assume that 169.5 (down), 175.5 (up) samples are used)
    scaleUncReductionFactor = TMath::Sqrt(2.0); // applied only to FSR uncertainties, as recommended for now

    cmsPrelimLabelMode = 1;

}

void Plotter13TeV::setMergeEnvelopeVariations(bool mergeEnvelopeVariations)
{
    mergeEnvelopeVariationsForBandsInCPs_ = mergeEnvelopeVariations;
    mergeScaleVariations_ = mergeEnvelopeVariations;
    mergeBTagVariations_ = mergeEnvelopeVariations;//FIXME: steer proper variable or create smth similar
    mergeTriggerVariations_ = mergeEnvelopeVariations;//FIXME: steer proper variable or create smth similar
    mergeBFragVariations_ = mergeEnvelopeVariations;//FIXME: steer proper variable or create smth similar
    mergeColorRecVariations_ = mergeEnvelopeVariations;//FIXME: steer proper variable or create smth similar 
}

void Plotter13TeV::DoFitInRatio(bool doFit)
{
    doFit_ = doFit;
}


void Plotter13TeV::SetOutpath(TString path)
{
    outpath = path;
}


void Plotter13TeV::unfolding(TString channel, TString systematic, const TString& unfoldingType/* = "svd"*/) {

    std::cout << "\n\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;
    std::cout << "Starting Calculation of Diff. X Section for '" << name << "' in Channel '" << channel << "' and Systematic '" << systematic << "':" << std::endl;
    std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;

    if (channel == "combined" && !performUnfoldingInCombinedChannel)
    {
        ifstream ResultsEE("UnfoldingResults/" + systematic + "/ee/" + name + "Results.txt");
        ifstream ResultsEMu("UnfoldingResults/" + systematic + "/emu/" + name + "Results.txt");
        ifstream ResultsMuMu("UnfoldingResults/" + systematic + "/mumu/" + name + "Results.txt");

        if (!ResultsEE.is_open()) {
            Plotter13TeV::preunfolding("ee", systematic, unfoldingType);
            CalcDiffXSec("ee", systematic, unfoldingType);
        }
        else {
            ResultsEE.close();
        }
        if (!ResultsEMu.is_open()) {
            Plotter13TeV::preunfolding("emu", systematic, unfoldingType);
            CalcDiffXSec("emu", systematic, unfoldingType);
        }
        else {
            ResultsEMu.close();
        }
        if (!ResultsMuMu.is_open()) {
            Plotter13TeV::preunfolding("mumu", systematic, unfoldingType);
            CalcDiffXSec("mumu", systematic, unfoldingType);
        }
        else {
            ResultsMuMu.close();
        }
    }

    CalcDiffXSec(channel, systematic, unfoldingType);

    std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;
    std::cout << "Finished Calculation of Diff. X Section for '" << name << "' in Channel '" << channel << "' and Systematic '" << systematic << "':" << std::endl;
    std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;

}

void Plotter13TeV::preunfolding(TString Channel, TString Systematic, const TString& unfoldingType/* = "svd"*/)
{
    write(Channel, Systematic, unfoldingType);
}


void Plotter13TeV::DYScaleFactor(TString SpecialComment) {

    DYScale = { 1,1,1,1 };
    usefulTools13TeV->DYScaleFactor(SpecialComment, DYScale, name);
}

Plotter13TeV::Plotter13TeV()
{
    gErrorIgnoreLevel = 1001;
    name = "defaultName";
    specialComment = "Standard";
    rangemin = 0;
    rangemax = 3;
    YAxis = "N_{events}";
    initialized = false;
    datafiles = 0;

    channelLabel.insert(channelLabel.begin(), 4, "");

    outpath = "";
    outpathPlots = "Plots";
    outpathPlots = "UnfoldingResults";
    subfolderChannel = "";
    subfolderSpecial = "";

    fileReader = RootFileReader::getInstance();

    useTopMassSetup_ = false;
    isCP_ = false;

}

void Plotter13TeV::setOptions(TString name_, TString specialComment_, TString YAxis_, TString XAxis_, int rebin_, bool doDYScale_, bool logX_, bool logY_, double ymin_, double ymax_, double rangemin_, double rangemax_, int bins_, std::vector<double> XAxisbins_, std::vector<double> XAxisbinCenters_)
{
    name = name_; //Variable name you want to plot
    specialComment = specialComment_;
    YAxis = YAxis_; //Y-axis title
    XAxis = XAxis_; //X-axis title
    rebin = rebin_; //Nr. of bins to be merged together
    doDYScale = doDYScale_; //Apply DY scale factor?
    logX = logX_; //Draw X-axis in Log scale
    logY = logY_; //Draw Y-axis in Log scale
    ymin = ymin_; //Min. value in Y-axis
    ymax = ymax_; //Max. value in Y-axis
    rangemin = rangemin_; //Min. value in X-axis
    rangemax = rangemax_; //Max. value in X-axis
    bins = bins_; //Number of bins to plot
    XAxisbins.clear();
    XAxisbins = XAxisbins_; // Bins edges=bins+1
    XAxisbinCenters.clear();
    XAxisbinCenters = XAxisbinCenters_; //Central point for BinCenterCorrection=bins

    //Modify the X/Y-axis labels
    if (XAxis.Contains("band#bar{b}")) {//Histogram naming convention has to be smarter
        XAxis.ReplaceAll("band#bar{b}", 11, "b and #bar{b}", 13);
    }
    if (XAxis.Contains("tand#bar{t}")) {//Histogram naming convention has to be smarter
        XAxis.ReplaceAll("tand#bar{t}", 11, "t and #bar{t}", 13);
    }
    if (XAxis.Contains("l^{+}andl^{-}")) {//Histogram naming convention has to be smarter
        XAxis.ReplaceAll("l^{+}andl^{-}", 13, "l^{+} and l^{-}", 15);
    }
    if (YAxis.Contains("Toppairs")) {
        YAxis.ReplaceAll("Toppairs", 8, "Top-quark pairs", 15);
    }
    if (YAxis.Contains("Topquarks")) {
        YAxis.ReplaceAll("Topquarks", 9, "Top quarks", 10);
    }
    if (YAxis.Contains("Numberof")) {
        YAxis.ReplaceAll("Numberof", 8, "Number of ", 10);
    }

    DYScale.insert(DYScale.begin(), 4, 1.);//Initialize the DY scale-factor to (1., 1., 1., 1.)

    usefulTools13TeV = new UsefulTools13TeV(fileReader, doClosureTest, doDYScale);
    this->energy = usefulTools13TeV->energy;
    this->lumi = usefulTools13TeV->lumi;
    this->topxsec = usefulTools13TeV->topxsec;

}


void Plotter13TeV::setDataSet(std::vector<TString> dataset_, std::vector<double> scales_, std::vector<TString> legends_, std::vector<int> colors_, TString DYEntry_)
{
    dataset.clear();
    scales.clear();
    legends.clear();
    legendsSyst.clear();
    colors.clear();
    dataset = dataset_;
    scales = scales_;
    legends = legends_;
    colors = colors_;
    DYEntry = DYEntry_;
}

void Plotter13TeV::setDataSet(TString mode, TString Systematic)
{
    initialized = false;
    legendsSyst.clear();

    if (channelLabel.size() < 4) { channelLabel.insert(channelLabel.begin(), 4, ""); }

    if (mode == "ee") { channelType = 0;channelLabel.at(0) = "ee"; }
    if (mode == "mumu") { channelType = 1;channelLabel.at(1) = "#mu#mu"; }
    if (mode == "emu") { channelType = 2;channelLabel.at(2) = "e#mu"; }
    if (mode == "combined") { channelType = 3;channelLabel.at(3) = "Dilepton"; }

    if (!isCP_ && !particleXsec) {
        if (mode == "ee") { channelType = 0;channelLabel.at(0) = "ee: parton"; }
        if (mode == "mumu") { channelType = 1;channelLabel.at(1) = "#mu#mu: parton"; }
        if (mode == "emu") { channelType = 2;channelLabel.at(2) = "e#mu: parton"; }
        if (mode == "combined") { channelType = 3;channelLabel.at(3) = "Dilepton: parton"; }
    }

    if (!isCP_ && particleXsec) {
        if (mode == "ee") { channelType = 0;channelLabel.at(0) = "ee: particle"; }
        if (mode == "mumu") { channelType = 1;channelLabel.at(1) = "#mu#mu: particle"; }
        if (mode == "emu") { channelType = 2;channelLabel.at(2) = "e#mu: particle"; }
        if (mode == "combined") { channelType = 3;channelLabel.at(3) = "Dilepton: particle"; }
    }

    // Set dataset specific subfolders
    outpathPlots = "./Plots";
    outpathResults = "./UnfoldingResults";
    subfolderChannel = mode;
    subfolderChannel.Prepend("/");
    subfolderSpecial = "";
    if (specialComment.CompareTo("Standard") != 0) {
        //subfolderSpecial = specialComment.Prepend("/");
    }

    DYEntry = "Z+jets";

    if (Systematic.Contains("DY_") || Systematic.Contains("BG_")) { Systematic = "Nominal"; }//We just need to vary the nominal DY and BG systematics

    TString histoListName = "FileLists/HistoFileList_" + Systematic + "_" + mode + ".txt";
    if (estimateSysUncFromPseudoData) histoListName = "PseudoFileLists/HistoFileList_" + Systematic + "_" + mode + ".txt";
    std::cout << "reading " << histoListName << std::endl;

    //FIXME: 
    /// This block can be deleted if "datafiles" not using
    // Counting number of data files
    ifstream FileList(histoListName);
    if (FileList.fail()) {
        std::cerr << "Error reading " << histoListName << std::endl;
        exit(1);
    }
    TString filename;
    datafiles = 0;

    while (!FileList.eof()) {
        FileList >> filename;
        if (filename == "") { continue; }//Skip empty lines
        if (filename.Contains("run")) { datafiles++; }

    }
    FileList.close();
    /// This block can be deleted if "datafiles" not using

    //Fill 
    UsefulTools13TeV::fillLegendColorDataset(histoListName, legends, colors, dataset);

    if (usePoissonSmearedPseudoData) {

        if (!estimateSysUncFromPseudoData) std::cout << "Constructing pseudo-data as sum of 'nominal' ttbar and bkgd simulations" << std::endl;
        else std::cout << "Constructing pseudo-data as sum of 'varied' ttbar and bkgd simulations" << std::endl;
        TString nomFileListName = "FileLists/HistoFileList_Nominal_" + mode + ".txt";
        if (estimateSysUncFromPseudoData) nomFileListName = "FileLists/HistoFileList_" + Systematic + "_" + mode + ".txt";
        std::cout << "--- reading " << nomFileListName << " for pseudo-data" << std::endl;

        std::ifstream nomFileList(nomFileListName);
        if (nomFileList.fail()) {
            std::cerr << " Error reading " << nomFileListName << std::endl;
            exit(1);
        }

        TString nomFileName;
        pseudoDataset.clear();
        while (!nomFileList.eof()) {
            nomFileList >> nomFileName;
            if (nomFileName == "") { continue; }//Skip empty lines
            pseudoDataset.push_back(nomFileName);
        }
        nomFileList.close();

    }

}

bool Plotter13TeV::fillHisto()
{
    TH1::AddDirectory(kFALSE);
    if (initialized) { return true; }
    hists.clear();
    histsTtBgrOutOfSpace.clear();
    histsForPseudoData.clear();
    for (unsigned int i = 0; i < dataset.size(); i++) {
        //         std::cout << i << ": " << dataset.at(i) << std::endl;


        TH1D* hist = fileReader->GetClone<TH1D>(dataset.at(i), name, true);
        TH1D* histTtBgrOutOfSpace = NULL;
        TH1D* histForPseudoData = NULL;
        if (!hist) return false;
        if (usePoissonSmearedPseudoData) {
            if (pseudoDataset.size() != dataset.size()) {
                std::cout << "FileLists for nominal dataset and pseudo-data contain different number of input .root files!" << std::endl;
                return false;
            }
            histForPseudoData = fileReader->GetClone<TH1D>(pseudoDataset.at(i), name, true);
            if (!histForPseudoData) return false;
        }
        if (particleXsec && dataset.at(i).Contains("ttbarsignalplustau")) {
            histTtBgrOutOfSpace = fileReader->GetClone<TH1D>(dataset.at(i), name + "_OutOfSpace", true);
            if (!histTtBgrOutOfSpace) return false;
        }

        if (!revokeAntiQuantityCombination && !name.Contains("_step") && !name.Contains("bcp_") && !name.Contains("Lead") && !name.EndsWith("bkr") && !name.EndsWith("akr") && !name.EndsWith("afterbtag")
            && (name.Contains("Lepton") || name.Contains("BJet") || name.Contains("Top")))
        {
            TString stemp = name;
            if (name.Contains("Lepton")) { stemp.ReplaceAll("Lepton", 6, "AntiLepton", 10); }
            else if (name.Contains("BJet")) { stemp.ReplaceAll("BJet", 4, "AntiBJet", 8); }
            else if (name.Contains("Top")) { stemp.ReplaceAll("Top", 3, "AntiTop", 7); }
            const TH1D* other = fileReader->Get<TH1D>(dataset.at(i), stemp);
            if (other) hist->Add(other);
            else std::cerr << "Cannot find corresponding anti-quantity histogram: " << stemp << std::endl;

            if (particleXsec && dataset.at(i).Contains("ttbarsignalplustau")) {
                const TH1D* otherTtBkgd = fileReader->Get<TH1D>(dataset.at(i), stemp + "_OutOfSpace");
                if (otherTtBkgd) histTtBgrOutOfSpace->Add(otherTtBkgd);
                else std::cerr << "Cannot find corresponding anti-quantity histogram for tt bkgd: " << stemp + "_OutOfSpace" << std::endl;
            }

            if (usePoissonSmearedPseudoData) {
                const TH1D* otherForPseudoData = fileReader->Get<TH1D>(pseudoDataset.at(i), stemp);
                if (otherForPseudoData) histForPseudoData->Add(otherForPseudoData);
                else std::cerr << "Cannot find corresponding anti-quantity histogram for pseudo-data: " << stemp << std::endl;
            }
        }

        //Rescaling to the data luminosity
        double LumiWeight = usefulTools13TeV->CalcLumiWeight(dataset.at(i));
        //std::cout << "File " << dataset.at(i) << " has weight " << LumiWeight << "\n";

        usefulTools13TeV->ApplyFlatWeights(hist, LumiWeight);

        common::setHHStyle(*gStyle);
        hists.push_back(*hist);
        delete hist;

        if (particleXsec && dataset.at(i).Contains("ttbarsignalplustau")) {
            usefulTools13TeV->ApplyFlatWeights(histTtBgrOutOfSpace, LumiWeight);
            common::setHHStyle(*gStyle);
            if (dataset.at(i).Contains("ee_ttbarsignalplustau")) histTtBgrOutOfSpace->SetName("ee_" + (TString)histTtBgrOutOfSpace->GetName());
            if (dataset.at(i).Contains("emu_ttbarsignalplustau")) histTtBgrOutOfSpace->SetName("emu_" + (TString)histTtBgrOutOfSpace->GetName());
            if (dataset.at(i).Contains("mumu_ttbarsignalplustau")) histTtBgrOutOfSpace->SetName("mumu_" + (TString)histTtBgrOutOfSpace->GetName());
            histsTtBgrOutOfSpace.push_back(*histTtBgrOutOfSpace);
            delete histTtBgrOutOfSpace;
        }

        if (usePoissonSmearedPseudoData) {
            //Rescaling to the data luminosity
            double LumiWeightForPseudoData = usefulTools13TeV->CalcLumiWeight(pseudoDataset.at(i));
            usefulTools13TeV->ApplyFlatWeights(histForPseudoData, LumiWeightForPseudoData);
            common::setHHStyle(*gStyle);
            histsForPseudoData.push_back(*histForPseudoData);
        }
        delete histForPseudoData;

    }
    if (doClosureTest && !usePoissonSmearedPseudoData) {
        for (size_t i = 1; i < dataset.size(); ++i) {
            if (!dataset.at(i).Contains("ttbarsignalplustau")) {
                hists[0].Add(&hists[i]);
            }
        }
    }

    if (usePoissonSmearedPseudoData) {
        if (channelType != 3) {// for ee, emu or mumu channel
            for (size_t i = 1; i < pseudoDataset.size(); ++i) {
                if (!pseudoDataset.at(i).Contains("ttbarsignalplustau") && !pseudoDataset.at(i).Contains("ttbarbgviatau") && !pseudoDataset.at(i).Contains("run")) {
                    hists[0].Add(&histsForPseudoData[i]);
                }
            }
            Plotter13TeV::performPoissonSmearing(&hists[0]);
        }
        else if (channelType == 3) {
            //Indeces must match to the ordering of data files in filelist {syst}_combined.txt
            TString chToInclude[3] = { "ee_", "emu_", "mumu_" };
            for (int ch_it = 0; ch_it < 3; ++ch_it) {
                for (size_t i = ch_it + 1; i < pseudoDataset.size(); ++i) { //skipping pseudo-data files (using +1) to which MC is added
                    if (!pseudoDataset.at(i).Contains("ttbarsignalplustau") && !pseudoDataset.at(i).Contains("ttbarbgviatau") && !pseudoDataset.at(i).Contains("run") && pseudoDataset.at(i).Contains(chToInclude[ch_it])) {
                        hists[ch_it].Add(&histsForPseudoData[i]);
                    }
                }
                Plotter13TeV::performPoissonSmearing(&hists[ch_it]);
            }

        }
    }
    initialized = true;
    return true;
}


void Plotter13TeV::addAndDelete_or_Assign(TH1*& addToThis, TH1* addThis, bool removeNegExpectation, bool isUnfolding)
{
    // removeNegExpectation:
    // Here it is used for binning of diff. xsecs (important)
    if (removeNegExpectation && isUnfolding) removeNegativeExpectationsInBins(addThis);

    if (!addToThis) addToThis = addThis;
    else {
        addToThis->Add(addThis);
        delete addThis;
    }
}


void Plotter13TeV::removeNegativeExpectationsInBins(TH1*& thisHisto)
{
    // removeNegExpectation:
    // Needed for exclusion of bins with negative entries due to samples with negative weights (eg, aMC@NLO)
    // Needed in order to avoid negative expectation for one process and exclude its impact on other processes in the histogram stack
    // Safe to use with backgrounds, but don't use for ttbar signal, otherwise resp. will be screwed -> ensure to have enough statistics in the bin for diff. meas. 

    for (int ibin = 0; ibin <= thisHisto->GetNbinsX(); ++ibin) {
        //Check visible bins + underflow
        if (thisHisto->GetBinContent(ibin) < 0.0) {
            thisHisto->SetBinContent(ibin, 0.0);// set prediction to 0.0
            thisHisto->SetBinError(ibin, 1.0);// set error to 1 in order to match expectation of +/-1 event 
        }
    }
    if (thisHisto->GetBinContent(thisHisto->GetNbinsX() + 1) < 0.0) {
        thisHisto->SetBinContent(thisHisto->GetNbinsX() + 1, 0.0);// set prediction to 0.0
        thisHisto->SetBinError(thisHisto->GetNbinsX() + 1, 1.0);// set error to 1 in order to match expectation of +/-1 event 
    }
}

void Plotter13TeV::removeNegativeExpectationsInBinsDouble(TH1D*& thisHisto)
{
    // removeNegExpectation:
    // Needed for exclusion of bins with negative entries due to samples with negative weights (eg, aMC@NLO)
    // Needed in order to avoid negative expectation for one process and exclude its impact on other processes in the histogram stack
    // Safe to use with backgrounds, but don't use for ttbar signal, otherwise resp. will be screwed -> ensure to have enough statistics in the bin for diff. meas. 

    for (int ibin = 0; ibin <= thisHisto->GetNbinsX(); ++ibin) {
        //Check visible bins + underflow
        if (thisHisto->GetBinContent(ibin) < 0.0) {
            thisHisto->SetBinContent(ibin, 0.0);// set prediction to 0.0
            thisHisto->SetBinError(ibin, 1.0);// set error to 1 in order to match expectation of +/-1 event 
        }
    }
    if (thisHisto->GetBinContent(thisHisto->GetNbinsX() + 1) < 0.0) {
        thisHisto->SetBinContent(thisHisto->GetNbinsX() + 1, 0.0);// set prediction to 0.0
        thisHisto->SetBinError(thisHisto->GetNbinsX() + 1, 1.0);// set error to 1 in order to match expectation of +/-1 event 
    }
}


void Plotter13TeV::write(TString Channel, TString Systematic, const TString& unfoldingType/* = "svd"*/) // do scaling, stacking, legending, and write in file 
{
    // OZ 4.01.2017 for TUnfold fine binning is needed
    std::vector<double> XAxisbins = this->XAxisbins;
    int bins = this->bins;
    if (unfoldingType == "tunfold")
    {
        printf("Plotter13TeV::write(): reading fine binning for TUnfold\n");
        // PlotterForTUnfold is needed only to read the fine binning (should be provided in ControlCardsTUnfold/ directory)
        Samples samples;
        PlotterForTUnfold plotter(samples, 0.0, 0.0);
        std::vector<TString> v_plotName;
        v_plotName.push_back(name);
        v_plotName.push_back("1d");
        plotter.setOptions(v_plotName);
        XAxisbins = plotter.GetFineBinning()[0];
        bins = XAxisbins.size() - 1;
        printf("Plotter13TeV::write(): bins = %d from %f tp %f for %s\n", bins, XAxisbins[0], XAxisbins[bins], name.Data());
    }

    setDataSet(Channel, Systematic);
    if (!fillHisto()) return;
    if (hists.size() == 0) {
        std::cerr << "***ERROR! No histograms available! " << Channel << "/" << Systematic << std::endl;
        exit(11);
    }

    auto c = make_shared<TCanvas>("", "");
    auto stack = make_shared<THStack>("def", "def");
    TLegend* leg = new TLegend();

    std::vector<TH1*> drawhists;
    drawhists.resize(hists.size());

    std::stringstream ss;
    ss << DYScale[channelType];
    TString scale;
    scale = (TString)ss.str();
    int legchange = 0;
    leg->Clear();
    c->Clear();
    c->SetName("");

    c->SetTitle("");
    TH1* aDataHist = NULL;
    TH1* aBgrHist = NULL;
    TH1* aTtBgrHist = NULL;
    TH1* aTtBgrOutOfSpaceHist = NULL;
    TH1* aRecHist = NULL;
    TH1* aGenHist = NULL;
    TH1* aPseudoHist = NULL;
    TH1* aRespHist = NULL;
    double ttbarSignalLumiWeight = 0.0;

    double Xbins[XAxisbins.size()];
    TString newname = name;

    if (name.Contains("Hyp")) {//Histogram naming convention has to be smarter
        newname.ReplaceAll("Hyp", 3, "", 0);
    }

    bool init = false;

    for (unsigned int i = 0; i < XAxisbins.size();i++) { Xbins[i] = XAxisbins[i]; }

    std::unique_ptr<TH1> sumttbar;
    std::unique_ptr<TH1> allttbar;

    // OZ 18.01.2017 determine fine binning
    // The code below is needed if one wants to estimate possible "fine" (detector level) binning needed for TUnfold.
    // Fine bin edges need to match existing ones in the stored histograms; 
    // statistics should be gaussian in each fine bin (i.e. at least ~ 30 events);
    // this will be done automatically using available statistics in data.
    // To determine fine binning, uncomment code below and run.
    // NB hists[1] has larger statistics than hists[0]
    /*if(unfoldingType == "svd")
    {
      TH1D* hCoarse = new TH1D("hCoarse", "hCoarse", bins, &XAxisbins[0]);
      printf("VAR: %s\n", name.Data());
      TString message;
      std::vector<double> fineBinning = PlotterForTUnfold::DetermineFineBinning(hCoarse, &hists[1], &message);
      TString fileName = TString::Format("ControlCardsTUnfold/%s.txt", name.Data());
      FILE* fout = fopen(fileName.Data(), "w");
      fprintf(fout, "%s", message.Data());
      fclose(fout);
      printf("TUnfold control card for var %s written in %s\n", name.Data(), fileName.Data());
      delete hCoarse;
      //return;
      //exit(1);
    }*/

    for (unsigned int i = 0; i < hists.size(); i++) { // prepare histos and leg
        drawhists[i] = (TH1D*)hists[i].Clone();//rebin and scale the histograms
        if (rebin > 1) drawhists[i]->Rebin(rebin);

        //In case you would like to produce histogram with inclusive last bin for jet multiplicity
        /*if(name.Contains("HypjetMulti")){
            double incJetContent = 0.0;
            incJetContent = drawhists[i]->GetBinContent(7) + drawhists[i]->GetBinContent(8) + drawhists[i]->GetBinContent(9) + drawhists[i]->GetBinContent(10);
            drawhists[i]->SetBinContent(7, incJetContent);
            drawhists[i]->SetBinContent(8, 0.0);
            drawhists[i]->SetBinContent(9, 0.0);
            drawhists[i]->SetBinContent(10, 0.0);
        }*/

        if (XAxisbins.size() > 1) drawhists[i] = drawhists[i]->Rebin(bins, "tmp", Xbins);
        setStyle(drawhists[i], i, true);

        // Needed for exclusion of bins with negative entries due to samples with negative weights (eg, aMC@NLO)
        // Done for binning in CPs (important)
        // Needed in order to avoid negative expectation for one process and exclude its impact on other processes in the histogram stack
        if (legends.at(i) != "Data" && !doUnfolding) removeNegativeExpectationsInBins(drawhists[i]);

        if (legends.at(i) == "t#bar{t} signal") {
            if (sumttbar.get()) sumttbar->Add(drawhists[i]);
            else sumttbar = std::unique_ptr<TH1>{ static_cast<TH1*>(drawhists[i]->Clone()) };
        }
        if (legends.at(i) == "t#bar{t} signal" || legends.at(i) == "t#bar{t} other") {
            if (allttbar.get()) allttbar->Add(drawhists[i]);
            else allttbar = std::unique_ptr<TH1>{ static_cast<TH1*>(drawhists[i]->Clone()) };
        }
    }
    for (unsigned int i = 0; i < hists.size(); ++i) { // prepare histos and leg
        //         std::cout << "Legend ["<<i<<"] = " << legends.at(i) << std::endl;
        if (legends.at(i) != "Data") {
            //      drawhists[i]->Scale(12.1/5.1);

            if (XAxisbins.size() > 1) {//only distributions we want to unfold will have a binning vector
                //if(XAxisbins.size()>1||1){//only distributions we want to unfold will have a binning vector //to compare
                if (legends.at(i) == "t#bar{t} signal" && doUnfolding) {
                    TString ftemp = dataset.at(i);
                    ttbarSignalLumiWeight = usefulTools13TeV->CalcLumiWeight(dataset.at(i));
                    if (!init) {
                        aRespHist = fileReader->GetClone<TH2>(ftemp, respMatrixPrefix + newname);
                        aGenHist = fileReader->GetClone<TH1D>(ftemp, "VisGen" + newname);
                        if (particleXsec) {
                            aPseudoHist = fileReader->GetClone<TH1D>(ftemp, "VisPseudo" + newname);
                            aTtBgrOutOfSpaceHist = fileReader->GetClone<TH1D>(ftemp, "Hyp" + newname + "_OutOfSpace");
                        }
                        if (!revokeAntiQuantityCombination && !newname.Contains("Lead") && (newname.Contains("Lepton") || newname.Contains("Top") || newname.Contains("BJet"))) {
                            aRespHist->Add(fileReader->Get<TH2>(ftemp, respMatrixPrefix + "Anti" + newname));
                            aGenHist->Add(fileReader->Get<TH1D>(ftemp, "VisGenAnti" + newname));
                            if (particleXsec) {
                                aPseudoHist->Add(fileReader->Get<TH1D>(ftemp, "VisPseudoAnti" + newname));
                                aTtBgrOutOfSpaceHist->Add(fileReader->Get<TH1D>(ftemp, "HypAnti" + newname + "_OutOfSpace"));
                            }
                        }
                        init = true;
                    }
                    else {//account for more than one signal histogram
                        aRespHist->Add(fileReader->Get<TH2>(ftemp, respMatrixPrefix + newname));
                        aGenHist->Add(fileReader->Get<TH1D>(ftemp, "VisGen" + newname));
                        if (particleXsec) {
                            aPseudoHist->Add(fileReader->Get<TH1D>(ftemp, "VisPseudo" + newname));
                            aTtBgrOutOfSpaceHist->Add(fileReader->Get<TH1D>(ftemp, "Hyp" + newname + "_OutOfSpace"));
                        }
                        if (!revokeAntiQuantityCombination && !newname.Contains("Lead") && (newname.Contains("Lepton") || newname.Contains("Top") || newname.Contains("BJet"))) {
                            aRespHist->Add(fileReader->Get<TH2>(ftemp, respMatrixPrefix + "Anti" + newname));
                            aGenHist->Add(fileReader->Get<TH1D>(ftemp, "VisGenAnti" + newname));
                            if (particleXsec) {
                                aPseudoHist->Add(fileReader->Get<TH1D>(ftemp, "VisPseudoAnti" + newname));
                                aTtBgrOutOfSpaceHist->Add(fileReader->Get<TH1D>(ftemp, "HypAnti" + newname + "_OutOfSpace"));
                            }
                        }
                    }
                }
                //std::cout<<"Legend: "<<legends.at(i)<<std::endl;
                if (legends.at(i) == "t#bar{t} signal") {
                    addAndDelete_or_Assign(aRecHist, drawhists[i]->Rebin(bins, "aRecHist", Xbins), 1, doUnfolding);
                }
                else if (legends.at(i) == "t#bar{t} other") {//IMPORTANT: TTbar Other are added to the ttbarbackground histogram AND the Background Hist gram
                    addAndDelete_or_Assign(aTtBgrHist, drawhists[i]->Rebin(bins, "aTtBgrHist", Xbins), 1, doUnfolding);
                    addAndDelete_or_Assign(aBgrHist, drawhists[i]->Rebin(bins, "aTtBgrHist", Xbins), 1, doUnfolding);
                }
                else if ((legends.at(i) == DYEntry)) {
                    drawhists[i]->Scale(DYScale.at(channelType));

                    //Here we take into account the systematic shifts needed for DY systematic because it only modifies the nominal dataset
                    if (Systematic == "DY_UP") {
                        drawhists[i]->Scale(1.3);
                    }
                    if (Systematic == "DY_DOWN") {
                        drawhists[i]->Scale(0.7);
                    }
                    addAndDelete_or_Assign(aBgrHist, drawhists[i]->Rebin(bins, "aBgrHist", Xbins), 1, doUnfolding);
                }
                else {
                    //Here we take into account the systematic shifts needed for BG systematic because it only modifies the nominal dataset
                    if (Systematic == "BG_UP") {
                        drawhists[i]->Scale(1.3);
                    }
                    if (Systematic == "BG_DOWN") {
                        drawhists[i]->Scale(0.7);
                    }
                    addAndDelete_or_Assign(aBgrHist, drawhists[i]->Rebin(bins, "aBgrHist", Xbins), 1, doUnfolding);
                }
            }
            else if ((legends.at(i) == DYEntry) && DYScale.at(channelType) != 1) {
                drawhists[i]->Scale(DYScale.at(channelType));
                //std::cout<<"DY SF propagated to CP: "<<DYScale.at(channelType)<<"in channel: "<<channelType<<std::endl;
            }

            if (i > 0) {
                if (legends.at(i) != legends.at(i - 1)) {
                    legchange = i;
                    if ((legends.at(i) == DYEntry) && DYScale.at(channelType) != 1) leg->AddEntry(drawhists[i], legends.at(i), "f");
                    else leg->AddEntry(drawhists[i], legends.at(i), "f");
                }
                else {
                    drawhists[legchange]->Add(drawhists[i]);
                }
            }
            if (i != (hists.size() - 1)) {
                if (legends.at(i) != legends.at(i + 1)) {
                    drawhists[i]->SetLineColor(1);
                }
            }            
else {
                drawhists[i]->SetLineColor(1);
            }
            if (legends.at(i) != legends.at(i - 1)) {
                drawhists[i]->SetLineColor(1);
                if (!legends.at(i).Contains("QCD") || addQCDToControlPlot()) { stack->Add(drawhists[i]); };

            }
        }
        else {
            if (i == 0) leg->AddEntry(drawhists[i], legends.at(i), "pe");
            if (i > 0) {
                if (legends.at(i) != legends.at(i - 1)) {
                    leg->AddEntry(drawhists[i], legends.at(i), "pe");
                }
                if (legends.at(i) == legends.at(0)) {
                    drawhists[0]->Add(drawhists[i]);
                }
            }
        }
    }

    if (XAxisbins.size() > 1 && doUnfolding) {//only distributions we want to unfold will have a binning vector
        aDataHist = drawhists[0]->Rebin(bins, "aDataHist", Xbins);
        aRespHist->Scale(ttbarSignalLumiWeight);
        aGenHist->Scale(ttbarSignalLumiWeight);

        if (doClosureTestWithLowStat && doClosureTest) Plotter13TeV::performPoissonSmearing(aDataHist);
        if (particleXsec) {
            aPseudoHist->Scale(ttbarSignalLumiWeight);
            //treat out-of-space bkgd as aRecoHist histogram
            usefulTools13TeV->ApplyFlatWeights(aTtBgrOutOfSpaceHist, ttbarSignalLumiWeight);

            aTtBgrOutOfSpaceHist = aTtBgrOutOfSpaceHist->Rebin(bins, "aTtBgrOutOfSpaceHist", Xbins);
            removeNegativeExpectationsInBins(aTtBgrOutOfSpaceHist);
            aRecHist->Add(aTtBgrOutOfSpaceHist, -1.0);
            aTtBgrHist->Add(aTtBgrOutOfSpaceHist);
            aBgrHist->Add(aTtBgrOutOfSpaceHist);
        }

        TString outdir = ttbar::assignFolder("preunfolded", Channel, Systematic);
        TFile* f15 = new TFile(outdir.Copy() + name + "_UnfoldingHistos.root", "RECREATE");
        aDataHist->Write("aDataHist"); delete aDataHist;
        aTtBgrHist->Write("aTtBgrHist"); delete aTtBgrHist;
        aBgrHist->Write("aBgrHist"); delete aBgrHist;
        aGenHist->Write("aGenHist"); delete aGenHist;
        if (particleXsec) {
            aPseudoHist->Write("aPseudoHist"); delete aPseudoHist;
            aTtBgrOutOfSpaceHist->Write("aTtBgrOutOfSpaceHist"); delete aTtBgrOutOfSpaceHist;
        }
        aRespHist->Write("aRespHist"); delete aRespHist;
        aRecHist->Write("aRecHist"); delete aRecHist;

        f15->Close();
        delete f15;
    }
    if (doUnfolding) return;
    TLegend* leg1 = (TLegend*)leg->Clone("leg1");
    TLegend* leg2 = (TLegend*)leg->Clone("leg2");
    setControlPlotLegendStyle(drawhists, legends, leg, leg1, leg2);

    if (name.Contains("HypjetMultiXSec")) {

        double InclusiveXsectionWrite[4], InclusiveXsectionStatErrorWrite[4];
        CalcXSec(dataset, pseudoDataset, InclusiveXsectionWrite, InclusiveXsectionStatErrorWrite, Systematic, "");

        Plotter13TeV::MakeTable(Channel, Systematic);
        if (channelType == 3) Plotter13TeV::PlotXSec(Channel);
    }

    drawhists[0]->SetMinimum(ymin);

    if (rangemin != 0 || rangemax != 0) { drawhists[0]->SetAxisRange(rangemin, rangemax, "X"); }

    if (logY)c->SetLogy();
    if (ymax == 0) {
        if (logY) { drawhists[0]->SetMaximum(18 * drawhists[0]->GetBinContent(drawhists[0]->GetMaximumBin())); }
        else { drawhists[0]->SetMaximum(1.5 * drawhists[0]->GetBinContent(drawhists[0]->GetMaximumBin())); }
    }
    else { drawhists[0]->SetMaximum(ymax); }

    if (name.Contains("HypTopRapidity") || name.Contains("HypTTBarRapidity") || (name.Contains("HypAntiTopRapidity") && revokeAntiQuantityCombination)) drawhists[0]->GetXaxis()->SetNdivisions(511);

    drawhists[0]->GetXaxis()->SetNoExponent(kTRUE);
    TGaxis::SetMaxDigits(4);


    //Removal of extra ticks in JetMult plots
    if (name.Contains("HypJetMultpt")) {
        drawhists[0]->GetXaxis()->SetNdivisions(drawhists[0]->GetNbinsX(), 0, 0, 1);

        TString TitBin = "";
        for (int bin = 1; bin <= drawhists[0]->GetNbinsX(); bin++) {
            if (bin == drawhists[0]->GetNbinsX()) { TitBin += "#geq"; TitBin += drawhists[0]->GetBinCenter(bin); drawhists[0]->GetXaxis()->SetBinLabel(bin, TitBin); }
            else {
                TitBin += drawhists[0]->GetBinCenter(bin);
                drawhists[0]->GetXaxis()->SetBinLabel(bin, TitBin);
            }
            TitBin = "";
        }
    }

    /*if(name.Contains("HypjetMulti")) {
        drawhists[0]->GetXaxis()->SetNdivisions(drawhists[0]->GetNbinsX(),0,0, 1);

        TString TitBin = "";
        for(int bin = 1; bin <= drawhists[0]->GetNbinsX(); bin++) {
            if( bin == 7) {TitBin += "#geq"; TitBin += drawhists[0]->GetBinCenter(bin); drawhists[0]->GetXaxis()->SetBinLabel(bin,TitBin);}
            else{TitBin += drawhists[0]->GetBinCenter(bin);
            drawhists[0]->GetXaxis()->SetBinLabel(bin,TitBin);
            }
            TitBin  = "";
        }
    }*/



    //Add the binwidth to the yaxis in yield plots

    TString ytitle = drawhists[0]->GetYaxis()->GetTitle();
    double binwidth = drawhists[0]->GetXaxis()->GetBinWidth(1);
    std::ostringstream width;
    width << binwidth;

    if (name.Contains("DeltaR") || name.Contains("Rapidity") || name.Contains("Eta") || name.Contains("Phi") || name.Contains("Fraction") || (name.Contains("Mass") && name.Contains("1st"))) { ytitle.Append(" / ").Append(width.str()); }
    else if ((name.Contains("pT", TString::kIgnoreCase) && !name.Contains("JetMult")) || (name.Contains("Mass", TString::kIgnoreCase) && !name.Contains("1st")) || name.Contains("MET") || name.Contains("HT")) { ytitle.Append(" / ").Append(width.str()).Append(" GeV"); };
    drawhists[0]->GetYaxis()->SetTitle(ytitle);
    drawhists[0]->Draw("e1");
    gStyle->SetEndErrorSize(0);

    stack->Draw("same HIST");

    TH1* stacksum = common::summedStackHisto(stack.get());
    TH1* uncBand = nullptr, * uncBandPlot = nullptr;

    if (mergeEnvelopeVariationsForBandsInCPs_) {

        sumttbar->SetName(name + "_signalmc");
        allttbar->SetName(name + "_allttbar");
        stacksum->SetName(name + "_allmc");

        GetEnvelopeVariationsForBandsInCPs(sumttbar.get(), Systematic, Channel);
        GetEnvelopeVariationsForBandsInCPs(allttbar.get(), Systematic, Channel);
        GetEnvelopeVariationsForBandsInCPs(stacksum, Systematic, Channel);

    }

    if (drawUncBand_) {

        if (drawTotalMCErrorForCP) uncBand = dynamic_cast<TH1*>(stacksum->Clone("uncBand"));//common::summedStackHisto(stack.get());
        else uncBand = dynamic_cast<TH1*>(allttbar->Clone("uncBand"));//common::summedStackHisto(stack.get()); 

        getSignalUncertaintyBand(uncBand, Channel);
        uncBand->SetFillStyle(3354);
        uncBand->SetFillColor(kBlack);
        uncBand->SetLineColor(kBlack);
        gStyle->SetHatchesLineWidth(1);
        gStyle->SetHatchesSpacing(0.8);
        uncBand->SetMarkerStyle(0);

        uncBandPlot = dynamic_cast<TH1*> (uncBand->Clone("uncBandPlot"));
        for (int i = 0; i <= stacksum->GetNbinsX(); i++) {
            uncBandPlot->SetBinContent(i, stacksum->GetBinContent(i));
        }
        uncBandPlot->Draw("same,e2");
        leg->AddEntry(uncBand, "Uncertainty", "f");
        if (leg2) {
            leg2->AddEntry(uncBand, "Uncertainty", "f");
            // stupid resizing of the legend to have same size if leg1 and leg2 have different number of entries
            const float y1 = leg2->GetY1NDC();
            const float y2 = leg2->GetY2NDC();
            const float deltaY = std::fabs(y2 - y1);
            const int nentriesLeg1 = leg1->GetNRows();
            const int nentriesLeg2 = leg2->GetNRows();
            leg2->SetY1NDC(y2 - 1. * nentriesLeg2 / nentriesLeg1 * deltaY);
        }
    }


    gPad->RedrawAxis();
    TExec* setex1{ new TExec("setex1","gStyle->SetErrorX(0.5)") };//this is frustrating and stupid but apparently necessary...
    setex1->Draw();
    if (drawUncBand_)uncBandPlot->Draw("same,e2");
    TExec* setex2{ new TExec("setex2","gStyle->SetErrorX(0.)") };
    setex2->Draw();
    drawhists[0]->Draw("same,e1");

    DrawCMSLabels(cmsPrelimLabelMode);
    DrawDecayChLabel(channelLabel[channelType]);

    if (name.Contains("JetMult")) {
        TString legtit = "";
        if (name.Contains("pt60")) legtit += "p_{T}^{jet}> 60 GeV";
        else if (name.Contains("pt100")) legtit += "p_{T}^{jet}> 100 GeV";
        else legtit += "p_{T}^{jet}> 30 GeV, |#eta^{jet}| < 2.4";     //###########
        leg->SetHeader(legtit);
    }

    //    leg->Draw("SAME");
    if (leg1) leg1->Draw("SAME");
    if (leg2) leg2->Draw("SAME");

    if (drawPlotRatio) {
        double yminCP_ = 0.49, ymaxCP_ = 1.51;
        yRangeControlPlotRatio(yminCP_, ymaxCP_);
        common::drawRatio(drawhists[0], stacksum, uncBand, yminCP_, ymaxCP_, doFit_);
    }

    // Create Directory for Output Plots 
    TString outdir = ttbar::assignFolder(outpathPlots, Channel, Systematic);
    c->Print(outdir.Copy() + name + ".eps");
    c->Print(outdir.Copy() + name + ".pdf");

    // Get the ratio plot from the canvas
    TPad* tmpPad = dynamic_cast<TPad*>(c->GetPrimitive("rPad"));
    TH1* ratio = nullptr;
    if (tmpPad) ratio = dynamic_cast<TH1*>(tmpPad->GetPrimitive("ratio"));

    //save Canvas AND sources in a root file
    TFile out_root(outdir.Copy() + name + "_source.root", "RECREATE");
    drawhists[0]->Write(name + "_data");
    sumttbar->Write(name + "_signalmc");
    allttbar->Write(name + "_allttbar");
    stacksum->SetName(name);
    stacksum->Write(name + "_allmc");

    if (ratio && ratio->GetEntries())ratio->Write(name + "_ratio");
    c->Write(name + "_canvas");
    out_root.Close();
    c->Clear();

    if (useTopMassSetup_) {

        TH1D* mt_ttsignal = 0;
        TH1D* mt_ttother = 0;
        TH1D* mt_zjets = 0;
        TH1D* mt_wjets = 0;
        TH1D* mt_singlet = 0;
        TH1D* mt_diboson = 0;
        TH1D* mt_ttV = 0;

        TFile mt_out_root(outdir.Copy() + name + "_topmass.root", "RECREATE");

        for (unsigned int iHW = 0; iHW < drawhists.size(); ++iHW) {
            if (iHW > 0) {
                if (legends.at(iHW) == legends.at(iHW - 1)) continue;
            }

            if (legends.at(iHW) == "t#bar{t} signal") {
                if (!mt_ttsignal) mt_ttsignal = (TH1D*)drawhists.at(iHW)->Clone();
                else mt_ttsignal->Add(drawhists.at(iHW));
            }
            else if (legends.at(iHW) == "t#bar{t} other") {
                if (!mt_ttother) mt_ttother = (TH1D*)drawhists.at(iHW)->Clone();
                else mt_ttother->Add(drawhists.at(iHW));
            }
            else if (legends.at(iHW) == "Z+jets") {
                if (!mt_zjets) mt_zjets = (TH1D*)drawhists.at(iHW)->Clone();
                else mt_zjets->Add(drawhists.at(iHW));
            }
            else if (legends.at(iHW) == "W+jets") {
                if (!mt_wjets) mt_wjets = (TH1D*)drawhists.at(iHW)->Clone();
                else mt_wjets->Add(drawhists.at(iHW));
            }
            else if (legends.at(iHW) == "Single t") {
                if (!mt_singlet) mt_singlet = (TH1D*)drawhists.at(iHW)->Clone();
                else mt_singlet->Add(drawhists.at(iHW));
            }
            else if (legends.at(iHW) == "Diboson") {
                if (!mt_diboson) mt_diboson = (TH1D*)drawhists.at(iHW)->Clone();
                else mt_diboson->Add(drawhists.at(iHW));
            }
            else if (legends.at(iHW) == "t#bar{t}+Z/W") {
                if (!mt_ttV) mt_ttV = (TH1D*)drawhists.at(iHW)->Clone();
                else mt_ttV->Add(drawhists.at(iHW));
            }
            else if (legends.at(iHW) != "Data") {
                std::cerr << "ERROR: Legend entry " << legends.at(iHW) << " not linked to any output histogram. Maybe a new process was added?" << std::endl;
                exit(1);
            }
        }

        if (!mt_ttsignal || !mt_ttother || !mt_zjets || !mt_wjets || !mt_singlet || !mt_diboson || !mt_ttV) {
            std::cerr << "ERROR: one or more histograms for top mass extraction are not initialized" << std::endl;
            exit(1);
        }

        drawhists.at(0)->Write("mt_" + name + "_data");
        mt_ttsignal->Write("mt_" + name + "_ttsignal");
        mt_ttother->Write("mt_" + name + "_ttother");
        mt_singlet->Write("mt_" + name + "_singlet");
        mt_ttV->Write("mt_" + name + "_ttV");
        mt_zjets->Write("mt_" + name + "_zjets");
        mt_diboson->Write("mt_" + name + "_diboson");
        mt_wjets->Write("mt_" + name + "_wjets");

        delete mt_ttsignal;
        delete mt_ttother;
        delete mt_zjets;
        delete mt_wjets;
        delete mt_singlet;
        delete mt_diboson;
        delete mt_ttV;
    }

    for (TH1* h : drawhists) delete h;
}

void Plotter13TeV::setStyle(TH1* hist, unsigned int i, bool isControlPlot)
{
    TString phasespace = ""; //in case we need to specify vis/full
    hist->SetFillColor(colors[i]);
    hist->SetLineColor(colors[i]);
    hist->SetLineWidth(1);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.08);
    hist->GetYaxis()->SetTitleOffset(1.7);
    hist->GetXaxis()->SetLabelOffset(0.007);
    //hist->GetXaxis()->SetLabelSize(0.04);// 0.05 for HypjetMult or HypJetMultpt30 or 0.04 for all other Control Plots 
    hist->GetYaxis()->SetLabelOffset(0.007);

    if (legends.at(i) == "Data") {
        hist->SetFillColor(0);
        hist->SetMarkerStyle(20);
        hist->SetMarkerSize(1.2);
        hist->SetLineWidth(2);
        if ((name.Contains("pT", TString::kIgnoreCase) || name.Contains("Mass", TString::kIgnoreCase)) &&
            (!name.Contains("1st") && !name.Contains("Rapidity", TString::kIgnoreCase) && !name.Contains("Eta", TString::kIgnoreCase) && !name.Contains("Phi", TString::kIgnoreCase) && !name.Contains("JetMult") && !name.Contains("Fraction")
                && !name.Contains("Multiplicity", TString::kIgnoreCase))) {
            hist->GetXaxis()->SetTitle(XAxis + " #left[GeV#right]");
            hist->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d" + XAxis + "}" + " #left[GeV^{-1}#right]");
            if (absoluteXsec) hist->GetYaxis()->SetTitle("#frac{d#sigma^{" + phasespace + "}}{d" + XAxis + "}" + " #left[pb/GeV#right]");
        }
        else {
            hist->GetXaxis()->SetTitle(XAxis);
            hist->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d" + XAxis + "}");
            if (absoluteXsec) hist->GetYaxis()->SetTitle("#frac{d#sigma^{" + phasespace + "}}{d" + XAxis + "}" + " #left[pb#right]");
        }
        if (isControlPlot) hist->GetYaxis()->SetTitle(YAxis);
    }
}


void Plotter13TeV::PlotXSec(TString Channel) {

    TH1::AddDirectory(kFALSE);

    std::vector<TString> vec_systematic{ "MASS_", "MESCALE_", "MEFACSCALE_", "MERENSCALE_", "PSISRSCALE_", "PSFSRSCALE_", "BSEMILEP_", "BFRAG_", "BFRAG_PETERSON_", "ERDON_", "ERDONRETUNE_", "GLUONMOVETUNE_", "UETUNE_", "MATCH_", "PDF_ALPHAS_", "HAD_", "BTAG_", "BTAG_LJET_", "KIN_", "LEPT_", "TRIG_", "BG_", "DY_", "PU_", "JER_", "JES_" /*, "PDF_", "MOD_", */ };//For the time being until all systematics are finalished
    std::vector<TString> vec_channel{ "ee","mumu","emu","combined" };

    double BR_Error = BRError;
    double Lumi_Error = lumiError;

    double InclusiveXsectionPlot[4] = { 0. }, InclusiveXsectionStatErrorPlot[4] = { 0. }, InclusiveXsectionSysErrorPlot[4] = { 0. }, InclusiveXsectionTotalErrorPlot[4] = { 0. };
    for (int j = 0; j < (int)vec_channel.size(); j++) {
        TString outdir = ttbar::assignFolder(outpathPlots, vec_channel.at(j), TString("FinalResults"));
        ifstream SysResultsList("Plots/Nominal/" + vec_channel.at(j) + "/InclusiveXSec.txt");
        TString DUMMY;
        SysResultsList >> DUMMY >> DUMMY >> DUMMY >> DUMMY >> DUMMY >> InclusiveXsectionPlot[j] >> DUMMY >> InclusiveXsectionStatErrorPlot[j];
        SysResultsList.close();

        std::ofstream OutputFile(outdir.Copy() + "InclusiveXSecResultLateX.txt", std::ofstream::trunc);
        OutputFile << "Inclusive XSection Numerical Results for channel " << vec_channel.at(j) << std::endl;

        double syst_square_for_channel = 0.0;
        double tot_scale_variation = 0.0; // the 'envelope' (in this case out of integrated bin) type treatment is applied
        double tot_bfrag_variation = 0.0; // the 'envelope' (in this case out of integrated bin) type treatment is applied
        double tot_colorrec_variation = 0.0; // the 'envelope' (in this case out of integrated bin) type treatment is applied

        for (int i = 0; i < (int)vec_systematic.size(); i++) {

            ifstream SysUP, SysDOWN;

            if (vec_systematic.at(i) != "HAD_" && vec_systematic.at(i) != "MOD_" && vec_systematic.at(i) != "BFRAG_PETERSON_" &&
                vec_systematic.at(i) != "ERDON_" && vec_systematic.at(i) != "ERDONRETUNE_" && vec_systematic.at(i) != "GLUONMOVETUNE_") {
                SysUP.open("Plots/" + vec_systematic.at(i) + "UP/" + vec_channel.at(j) + "/InclusiveXSec.txt");
                SysDOWN.open("Plots/" + vec_systematic.at(i) + "DOWN/" + vec_channel.at(j) + "/InclusiveXSec.txt");
                if (!SysUP.is_open() || !SysDOWN.is_open()) continue;
            }
            else if (vec_systematic.at(i) == "HAD_") {
                SysUP.open("Plots/MCATNLO/" + vec_channel.at(j) + "/InclusiveXSec.txt");
                SysDOWN.open("Plots/MCATNLO/" + vec_channel.at(j) + "/InclusiveXSec.txt");
            }
            else if (vec_systematic.at(i) == "MOD_") {
                SysUP.open("Plots/POWHEG/" + vec_channel.at(j) + "/InclusiveXSec.txt");
                SysDOWN.open("Plots/POWHEG/" + vec_channel.at(j) + "/InclusiveXSec.txt");
            }
            else if (vec_systematic.at(i) == "BFRAG_PETERSON_") {
                SysUP.open("Plots/BFRAG_PETERSON/" + vec_channel.at(j) + "/InclusiveXSec.txt");
                SysDOWN.open("Plots/BFRAG_PETERSON/" + vec_channel.at(j) + "/InclusiveXSec.txt");
            }
            else if (vec_systematic.at(i) == "ERDON_") {
                SysUP.open("Plots/ERDON/" + vec_channel.at(j) + "/InclusiveXSec.txt");
                SysDOWN.open("Plots/ERDON/" + vec_channel.at(j) + "/InclusiveXSec.txt");
            }
            else if (vec_systematic.at(i) == "ERDONRETUNE_") {
                SysUP.open("Plots/ERDONRETUNE/" + vec_channel.at(j) + "/InclusiveXSec.txt");
                SysDOWN.open("Plots/ERDONRETUNE/" + vec_channel.at(j) + "/InclusiveXSec.txt");
            }
            else if (vec_systematic.at(i) == "GLUONMOVETUNE_") {
                SysUP.open("Plots/GLUONMOVETUNE/" + vec_channel.at(j) + "/InclusiveXSec.txt");
                SysDOWN.open("Plots/GLUONMOVETUNE/" + vec_channel.at(j) + "/InclusiveXSec.txt");
            }

            double VarUp = 0, VarDown = 0, StatErrUp = 0, StatErrDown = 0;

            SysUP >> DUMMY >> DUMMY >> DUMMY >> DUMMY >> DUMMY >> VarUp >> DUMMY >> StatErrUp;
            SysDOWN >> DUMMY >> DUMMY >> DUMMY >> DUMMY >> DUMMY >> VarDown >> DUMMY >> StatErrDown;
            SysUP.close();
            SysDOWN.close();

            //systematic error in %
            double sys_err = (TMath::Abs(InclusiveXsectionPlot[j] - VarUp) + TMath::Abs(InclusiveXsectionPlot[j] - VarDown)) * 0.5 / InclusiveXsectionPlot[j];
            if (vec_systematic.at(i) == "MASS_") sys_err = sys_err / massUncReductionFactor;
            if (vec_systematic.at(i) == "PSFSRSCALE_") sys_err = sys_err / scaleUncReductionFactor;
            if (vec_systematic.at(i) == "MESCALE_" || vec_systematic.at(i) == "MEFACSCALE_" || vec_systematic.at(i) == "MERENSCALE_" ||
                vec_systematic.at(i) == "PSISRSCALE_" || vec_systematic.at(i) == "PSFSRSCALE_") {
                if (sys_err > tot_scale_variation) tot_scale_variation = sys_err;
                OutputFile << vec_systematic.at(i) << " not propagated to total error directly (but included to related envelope and thus to total error) (%): " << setprecision(3) << sys_err * 100 << std::endl;
                continue;
            }
            if (vec_systematic.at(i) == "BFRAG_" || vec_systematic.at(i) == "BFRAG_PETERSON_") {
                if (sys_err > tot_bfrag_variation) tot_bfrag_variation = sys_err;
                OutputFile << vec_systematic.at(i) << " not propagated to total error directly (but included to related envelope and thus to total error) (%): " << setprecision(3) << sys_err * 100 << std::endl;
                continue;
            }
            if (vec_systematic.at(i) == "ERDON_" || vec_systematic.at(i) == "ERDONRETUNE_" || vec_systematic.at(i) == "GLUONMOVETUNE_") {
                if (sys_err > tot_colorrec_variation) tot_colorrec_variation = sys_err;
                OutputFile << vec_systematic.at(i) << " not propagated to total error directly (but included to related envelope and thus to total error) (%): " << setprecision(3) << sys_err * 100 << std::endl;
                continue;
            }
            if (vec_systematic.at(i) == "HAD_" || vec_systematic.at(i) == "MOD_") {
                OutputFile << vec_systematic.at(i) << " not propagated to total error (%): " << setprecision(3) << sys_err * 100 << std::endl;
                continue;
            }

            syst_square_for_channel += sys_err * sys_err;
            OutputFile << vec_systematic.at(i) << " (%): " << setprecision(3) << sys_err * 100 << std::endl;
        }
        // Consider tot_scale uncertainty, since it was not added
        syst_square_for_channel += tot_scale_variation * tot_scale_variation;
        syst_square_for_channel += tot_bfrag_variation * tot_bfrag_variation;
        syst_square_for_channel += tot_colorrec_variation * tot_colorrec_variation;

        OutputFile << "TOT_SCALE_ (%): " << setprecision(3) << tot_scale_variation * 100 << std::endl;
        OutputFile << "TOT_BFRAG_ (%): " << setprecision(3) << tot_bfrag_variation * 100 << std::endl;
        OutputFile << "TOT_COLORREC_ (%): " << setprecision(3) << tot_colorrec_variation * 100 << std::endl;
        OutputFile << "BranchingRatio (%): " << setprecision(3) << BR_Error * 100 << std::endl;
        OutputFile << "Luminosity (%): " << setprecision(3) << Lumi_Error * 100 << std::endl;

        InclusiveXsectionSysErrorPlot[j] = TMath::Sqrt(syst_square_for_channel + BR_Error * BR_Error + Lumi_Error * Lumi_Error);

        InclusiveXsectionTotalErrorPlot[j] = sqrt(InclusiveXsectionStatErrorPlot[j] * InclusiveXsectionStatErrorPlot[j] +
            InclusiveXsectionPlot[j] * InclusiveXsectionSysErrorPlot[j] * InclusiveXsectionPlot[j] * InclusiveXsectionSysErrorPlot[j]
        );

        OutputFile << "\n\n*******************************************************************************\n\n";
        OutputFile << " InclXsec[pb]     Stat.[pb]    Syst.[pb]   Total[pb]" << std::endl;
        OutputFile << setprecision(6) << InclusiveXsectionPlot[j] << " +- " << setprecision(3) << InclusiveXsectionStatErrorPlot[j] << " +- " << setprecision(4) << InclusiveXsectionSysErrorPlot[j] * InclusiveXsectionPlot[j] << " +- " << setprecision(4) << InclusiveXsectionTotalErrorPlot[j] << std::endl;
        OutputFile.close();
    }

    // measured results with statistical error
    Double_t mx[] = { 0.50,       1.50,       2.50,       3.50 };
    Double_t mexl[] = { 0.00,       0.00,       0.00,       0.00 };
    Double_t mexh[] = { 0.00,       0.00,       0.00,       0.00 };

    TGraphAsymmErrors* mplot = new TGraphAsymmErrors(4, mx, InclusiveXsectionPlot, mexl, mexh, InclusiveXsectionStatErrorPlot, InclusiveXsectionStatErrorPlot);
    mplot->SetMarkerStyle(20);
    mplot->GetYaxis()->SetNoExponent(kTRUE);
    mplot->SetMarkerColor(kBlack);
    mplot->SetMarkerSize(1.5);
    mplot->SetLineColor(kBlack);

    TGraphAsymmErrors* mplotwithsys = new TGraphAsymmErrors(4, mx, InclusiveXsectionPlot, mexl, mexh, InclusiveXsectionTotalErrorPlot, InclusiveXsectionTotalErrorPlot);
    mplotwithsys->SetMarkerStyle(20);
    mplotwithsys->SetMarkerColor(kBlack);
    mplotwithsys->SetMarkerSize(1.5);
    mplotwithsys->SetLineColor(kBlack);

    // kidonakis
    Double_t kidonmean = 234;
    Double_t kidonx[] = { -0.5,     0.5,   1.5,     2.5,     3.5,     4.5 };
    Double_t kidony[] = { kidonmean,kidonmean,kidonmean,kidonmean,kidonmean,kidonmean };
    Double_t kidonexl[] = { .4,    .4,      .5,      .5,      .5,      .5 };
    Double_t kidonexh[] = { .5,    .5,      .5,      .5,      .4,      .4 };
    Double_t kidoneyl[] = { 15.6,    15.6,    15.6,  15.6,    15.6,    15.6 };
    Double_t kidoneyh[] = { 13.9,    13.9,    13.9,  13.9,    13.9,    13.9 };

    TGraphAsymmErrors* kidonplot = new TGraphAsymmErrors(6, kidonx, kidony, kidonexl, kidonexh, kidoneyl, kidoneyh);
    kidonplot->SetLineColor(kGreen + 1);
    kidonplot->SetLineWidth(4);
    kidonplot->SetFillColor(kGreen + 1);
    kidonplot->SetFillStyle(3004);

    // Full NNLO
    Double_t nnlomean = 245.794;
    Double_t errDown = 10.656;
    Double_t errUp = 8.652;
    Double_t nnlox[] = { -0.5,     0.5,   1.5,     2.5,     3.5,     4.5 };
    Double_t nnloy[] = { nnlomean,nnlomean,nnlomean,nnlomean,nnlomean,nnlomean };
    Double_t nnloexl[] = { .4,    .4,      .5,      .5,      .5,      .5 };
    Double_t nnloexh[] = { .5,    .5,      .5,      .5,      .4,      .4 };
    Double_t nnloeyl[] = { errDown, errDown, errDown, errDown, errDown, errDown };
    Double_t nnloeyh[] = { errUp,   errUp,   errUp,   errUp,   errUp,   errUp };

    TGraphAsymmErrors* nnloplot = new TGraphAsymmErrors(6, nnlox, nnloy, nnloexl, nnloexh, nnloeyl, nnloeyh);
    nnloplot->SetLineColor(kGreen + 1);
    nnloplot->SetLineWidth(4);
    nnloplot->SetFillColor(kGreen + 1);
    nnloplot->SetFillStyle(3004);

    // TopLHC working group, prescription for m_top = 172.5 GeV
    //   https://indico.cern.ch/getFile.py/access?contribId=4&sessionId=1&resId=0&materialId=slides&confId=280522

    Double_t toplhcwgmean = 252.89;
    Double_t toplhcwgDown = 15.313;
    Double_t toplhcwgUp = 16.266;
    Double_t toplhcwgx[] = { -0.5,     0.5,   1.5,     2.5,     3.5,     4.5 };
    Double_t toplhcwgy[] = { toplhcwgmean,toplhcwgmean,toplhcwgmean,toplhcwgmean,toplhcwgmean,toplhcwgmean };
    Double_t toplhcwgexl[] = { .4,    .4,      .5,      .5,      .5,      .5 };
    Double_t toplhcwgexh[] = { .5,    .5,      .5,      .5,      .4,      .4 };
    Double_t toplhcwgeyl[] = { toplhcwgDown, toplhcwgDown, toplhcwgDown, toplhcwgDown, toplhcwgDown, toplhcwgDown };
    Double_t toplhcwgeyh[] = { toplhcwgUp,   toplhcwgUp,   toplhcwgUp,   toplhcwgUp,   toplhcwgUp,   toplhcwgUp };

    TGraphAsymmErrors* toplhcwgplot = new TGraphAsymmErrors(6, toplhcwgx, toplhcwgy, toplhcwgexl, toplhcwgexh, toplhcwgeyl, toplhcwgeyh);
    toplhcwgplot->SetLineColor(kGreen + 1);
    toplhcwgplot->SetLineWidth(4);
    toplhcwgplot->SetFillColor(kGreen + 1);
    toplhcwgplot->SetFillStyle(3004);

    // mcfm
    Double_t mcfmmean = 225.197;
    Double_t mcfmx[] = { -0.5,     0.5,     1.5,     2.5,     3.5,     4.5 };
    Double_t mcfmy[] = { mcfmmean,mcfmmean,mcfmmean,mcfmmean,mcfmmean,mcfmmean };
    Double_t mcfmexl[] = { .4,      .4,      .5,      .5,      .5,      .5 };
    Double_t mcfmexh[] = { .5,      .5,      .5,      .5,      .4,      .4 };
    Double_t mcfmeyl[] = { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0 };
    Double_t mcfmeyh[] = { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0 };

    TGraphAsymmErrors* mcfmplot = new TGraphAsymmErrors(6, mcfmx, mcfmy, mcfmexl, mcfmexh, mcfmeyl, mcfmeyh);
    mcfmplot->SetLineColor(kBlue + 1);
    mcfmplot->SetLineWidth(4);
    mcfmplot->SetFillColor(kBlue + 1);
    mcfmplot->SetFillStyle(3005);

    TH1F* framehist = new TH1F("framehist", "", 4, 0., 4.);
    framehist->SetMinimum(100);
    framehist->SetMaximum(380);
    framehist->GetXaxis()->SetTickLength(0);
    framehist->GetXaxis()->SetBinLabel(1, "");
    framehist->GetXaxis()->SetBinLabel(2, "");
    framehist->GetXaxis()->SetBinLabel(3, "");
    framehist->GetYaxis()->SetTitle("#sigma [pb]");
    framehist->GetYaxis()->CenterTitle(kTRUE);
    framehist->GetYaxis()->SetNoExponent(kTRUE);

    TPaveText* box1 = new TPaveText(0.25, 0.20, 0.33, 0.30, "NDC");
    box1->SetFillColor(10);
    box1->SetTextSize(0.04);
    box1->AddText("ee");

    TPaveText* box2 = new TPaveText(0.44, 0.20, 0.52, 0.30, "NDC");
    box2->SetFillColor(10);
    box2->SetTextSize(0.04);
    box2->AddText("#mu#mu");

    TPaveText* box3 = new TPaveText(0.62, 0.20, 0.72, 0.30, "NDC");
    box3->SetFillColor(10);
    box3->SetTextSize(0.04);
    box3->AddText("e#mu");

    TPaveText* box4 = new TPaveText(0.82, 0.20, 0.90, 0.30, "NDC");
    box4->SetFillColor(10);
    box4->SetTextSize(0.04);
    box4->AddText("combined");

    TLegend* leg = new TLegend(0.56, 0.65, 0.89, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextFont(62);
    leg->SetTextSize(0.03);
    leg->AddEntry(mplot, "Measurements", "p");
    leg->AddEntry(mcfmplot, "MCFM #otimes CTQE66M", "lf");
    //    leg->AddEntry( kidonplot,    "Kidonakis #otimes MSTW2008 NNLO",     "lf" );
    //    leg->AddEntry( nnloplot,    "NNLO #otimes MSTW2008 NNLO",     "lf" );
    leg->AddEntry(toplhcwgplot, "TOP LHC WG", "lf");

    TCanvas* c = new TCanvas("plot", "plot", 1200, 800);
    framehist->Draw();
    mcfmplot->Draw("C,2,SAME");
    //kidonplot->Draw("C,2,SAME");
    //nnloplot->Draw("C,2,SAME");
    toplhcwgplot->Draw("C,2,SAME");
    gStyle->SetEndErrorSize(8);
    mplot->Draw("p,SAME");
    mplotwithsys->Draw("p,SAME,Z");
    leg->Draw("SAME");

    box1->Draw("SAME");
    box2->Draw("SAME");
    box3->Draw("SAME");
    box4->Draw("SAME");

    TString outdir = ttbar::assignFolder(outpathPlots, Channel, TString("FinalResults"));
    c->Print(outdir.Copy() + "InclusiveXSec.eps");
    c->Print(outdir.Copy() + "InclusiveXSec.C");
    c->Clear();
    delete c;

}


void Plotter13TeV::MakeTable(TString Channel, TString Systematic) {

    TH1D* numhists5[hists.size()];
    TH1D* numhists6[hists.size()];
    TH1D* numhists7[hists.size()];
    TH1D* numhists8[hists.size()];
    TH1D* numhists9[hists.size()];

    TH1D* numhistsForPseudoData5[histsForPseudoData.size()];
    TH1D* numhistsForPseudoData6[histsForPseudoData.size()];
    TH1D* numhistsForPseudoData7[histsForPseudoData.size()];
    TH1D* numhistsForPseudoData8[histsForPseudoData.size()];
    TH1D* numhistsForPseudoData9[histsForPseudoData.size()];

    TString histoForYieldTable;
    if (usePoissonSmearedPseudoData) histoForYieldTable = "basic_jet_multiplicity";
    else histoForYieldTable = "events_weighted";

    for (unsigned int i = 0; i < dataset.size(); i++) {

        TH1D* temp_hist5 = fileReader->GetClone<TH1D>(dataset[i], histoForYieldTable + "_step4");
        TH1D* temp_hist6 = fileReader->GetClone<TH1D>(dataset[i], histoForYieldTable + "_step5");
        TH1D* temp_hist7 = fileReader->GetClone<TH1D>(dataset[i], histoForYieldTable + "_step6");
        TH1D* temp_hist8 = fileReader->GetClone<TH1D>(dataset[i], histoForYieldTable + "_step7");
        TH1D* temp_hist9 = fileReader->GetClone<TH1D>(dataset[i], histoForYieldTable + "_step8");

        double LumiWeight = usefulTools13TeV->CalcLumiWeight(dataset.at(i));
        usefulTools13TeV->ApplyFlatWeights(temp_hist5, LumiWeight);
        usefulTools13TeV->ApplyFlatWeights(temp_hist6, LumiWeight);
        usefulTools13TeV->ApplyFlatWeights(temp_hist7, LumiWeight);
        usefulTools13TeV->ApplyFlatWeights(temp_hist8, LumiWeight);
        usefulTools13TeV->ApplyFlatWeights(temp_hist9, LumiWeight);

        removeNegativeExpectationsInBinsDouble(temp_hist5);
        removeNegativeExpectationsInBinsDouble(temp_hist6);
        removeNegativeExpectationsInBinsDouble(temp_hist7);
        removeNegativeExpectationsInBinsDouble(temp_hist8);
        removeNegativeExpectationsInBinsDouble(temp_hist9);

        numhists5[i] = temp_hist5;
        numhists6[i] = temp_hist6;
        numhists7[i] = temp_hist7;
        numhists8[i] = temp_hist8;
        numhists9[i] = temp_hist9;

    }

    if (usePoissonSmearedPseudoData) {

        for (unsigned int i = 0; i < pseudoDataset.size(); i++) {

            TH1D* temp_hist5_pd = fileReader->GetClone<TH1D>(pseudoDataset[i], histoForYieldTable + "_step4");
            TH1D* temp_hist6_pd = fileReader->GetClone<TH1D>(pseudoDataset[i], histoForYieldTable + "_step5");
            TH1D* temp_hist7_pd = fileReader->GetClone<TH1D>(pseudoDataset[i], histoForYieldTable + "_step6");
            TH1D* temp_hist8_pd = fileReader->GetClone<TH1D>(pseudoDataset[i], histoForYieldTable + "_step7");
            TH1D* temp_hist9_pd = fileReader->GetClone<TH1D>(pseudoDataset[i], histoForYieldTable + "_step8");

            double LumiWeight_pd = usefulTools13TeV->CalcLumiWeight(pseudoDataset.at(i));
            usefulTools13TeV->ApplyFlatWeights(temp_hist5_pd, LumiWeight_pd);
            usefulTools13TeV->ApplyFlatWeights(temp_hist6_pd, LumiWeight_pd);
            usefulTools13TeV->ApplyFlatWeights(temp_hist7_pd, LumiWeight_pd);
            usefulTools13TeV->ApplyFlatWeights(temp_hist8_pd, LumiWeight_pd);
            usefulTools13TeV->ApplyFlatWeights(temp_hist9_pd, LumiWeight_pd);

            removeNegativeExpectationsInBinsDouble(temp_hist5_pd);
            removeNegativeExpectationsInBinsDouble(temp_hist6_pd);
            removeNegativeExpectationsInBinsDouble(temp_hist7_pd);
            removeNegativeExpectationsInBinsDouble(temp_hist8_pd);
            removeNegativeExpectationsInBinsDouble(temp_hist9_pd);

            numhistsForPseudoData5[i] = temp_hist5_pd;
            numhistsForPseudoData6[i] = temp_hist6_pd;
            numhistsForPseudoData7[i] = temp_hist7_pd;
            numhistsForPseudoData8[i] = temp_hist8_pd;
            numhistsForPseudoData9[i] = temp_hist9_pd;

        }

        if (Channel != "combined") {
            for (size_t i = 1; i < pseudoDataset.size(); ++i) {
                if (!pseudoDataset.at(i).Contains("ttbarsignalplustau") && !pseudoDataset.at(i).Contains("ttbarbgviatau") && !pseudoDataset.at(i).Contains("run")) {
                    numhists5[0]->Add(numhistsForPseudoData5[i]);
                    numhists6[0]->Add(numhistsForPseudoData6[i]);
                    numhists7[0]->Add(numhistsForPseudoData7[i]);
                    numhists8[0]->Add(numhistsForPseudoData8[i]);
                    numhists9[0]->Add(numhistsForPseudoData9[i]);
                }
            }
            Plotter13TeV::performPoissonSmearing(numhists5[0]);
            Plotter13TeV::performPoissonSmearing(numhists6[0]);
            Plotter13TeV::performPoissonSmearing(numhists7[0]);
            Plotter13TeV::performPoissonSmearing(numhists8[0]);
            Plotter13TeV::performPoissonSmearing(numhists9[0]);
        }
        else if (Channel == "combined") {
            //Indeces must match to the ordering of data files in filelist {syst}_combined.txt
            TString chToInclude[3] = { "ee_", "emu_", "mumu_" };
            for (int ch_it = 0; ch_it < 3; ++ch_it) {
                for (size_t i = ch_it + 1; i < pseudoDataset.size(); ++i) {//skipping pseudo-data files (using +1) to which MC is added
                    if (!pseudoDataset.at(i).Contains("ttbarsignalplustau") && !pseudoDataset.at(i).Contains("ttbarbgviatau") && !pseudoDataset.at(i).Contains("run") && pseudoDataset.at(i).Contains(chToInclude[ch_it])) {
                        numhists5[ch_it]->Add(numhistsForPseudoData5[i]);
                        numhists6[ch_it]->Add(numhistsForPseudoData6[i]);
                        numhists7[ch_it]->Add(numhistsForPseudoData7[i]);
                        numhists8[ch_it]->Add(numhistsForPseudoData8[i]);
                        numhists9[ch_it]->Add(numhistsForPseudoData9[i]);
                    }
                }
                Plotter13TeV::performPoissonSmearing(numhists5[ch_it]);
                Plotter13TeV::performPoissonSmearing(numhists6[ch_it]);
                Plotter13TeV::performPoissonSmearing(numhists7[ch_it]);
                Plotter13TeV::performPoissonSmearing(numhists8[ch_it]);
                Plotter13TeV::performPoissonSmearing(numhists9[ch_it]);
            }
        }
    }
    else {
        for (unsigned int i = 0; i < hists.size(); i++) { // prepare histos and leg
            if (legends.at(i) == DYEntry) {
                //numhists5[i]->Scale(DYScale[channelType]);//DYscale not applied in step5 and 6?
                numhists6[i]->Scale(DYScale[channelType]);
                numhists7[i]->Scale(DYScale.at(channelType));
                numhists8[i]->Scale(DYScale.at(channelType));
                numhists9[i]->Scale(DYScale.at(channelType));
            }
        }
    }

    ////////////////////////////Make output for tables
    double tmp_num5 = 0;
    double tmp_num6 = 0;
    double tmp_num7 = 0;
    double tmp_num8 = 0;
    double tmp_num9 = 0;

    TString outdir = ttbar::assignFolder(outpathPlots, Channel, Systematic);
    ofstream EventFile5; EventFile5.open(outdir.Copy() + "Events5.txt");
    ofstream EventFile6; EventFile6.open(outdir.Copy() + "Events6.txt");
    ofstream EventFile7; EventFile7.open(outdir.Copy() + "Events7.txt");
    ofstream EventFile8; EventFile8.open(outdir.Copy() + "Events8.txt");
    ofstream EventFile9; EventFile9.open(outdir.Copy() + "Events9.txt");

    double bg_num5 = 0;
    double bg_num6 = 0;
    double bg_num7 = 0;
    double bg_num8 = 0;
    double bg_num9 = 0;

    for (unsigned int i = 0; i < hists.size(); i++) {
        if (usePoissonSmearedPseudoData) {
            tmp_num5 += numhists5[i]->Integral(0, numhists5[i]->GetNbinsX() + 1);// considering over/underflow bins while using "basic_jet_multiplicity" histogram
            tmp_num6 += numhists6[i]->Integral(0, numhists6[i]->GetNbinsX() + 1);
            tmp_num7 += numhists7[i]->Integral(0, numhists7[i]->GetNbinsX() + 1);
            tmp_num8 += numhists8[i]->Integral(0, numhists8[i]->GetNbinsX() + 1);
            tmp_num9 += numhists9[i]->Integral(0, numhists9[i]->GetNbinsX() + 1);
        }
        else {
            tmp_num5 += numhists5[i]->Integral();
            tmp_num6 += numhists6[i]->Integral();
            tmp_num7 += numhists7[i]->Integral();
            tmp_num8 += numhists8[i]->Integral();
            tmp_num9 += numhists9[i]->Integral();
        }

        if (i == (hists.size() - 1)) {
            EventFile5 << legends.at(i) << ": " << tmp_num5 << std::endl;
            EventFile6 << legends.at(i) << ": " << tmp_num6 << std::endl;
            EventFile7 << legends.at(i) << ": " << tmp_num7 << std::endl;
            EventFile8 << legends.at(i) << ": " << tmp_num8 << std::endl;
            EventFile9 << legends.at(i) << ": " << tmp_num9 << std::endl;
            bg_num5 += tmp_num5;
            bg_num6 += tmp_num6;
            bg_num7 += tmp_num7;
            bg_num8 += tmp_num8;
            bg_num9 += tmp_num9;
            tmp_num5 = 0;
            tmp_num6 = 0;
            tmp_num7 = 0;
            tmp_num8 = 0;
            tmp_num9 = 0;
        }
        else if (legends.at(i) != legends.at(i + 1)) {
            EventFile5 << legends.at(i) << ": " << tmp_num5 << std::endl;
            EventFile6 << legends.at(i) << ": " << tmp_num6 << std::endl;
            EventFile7 << legends.at(i) << ": " << tmp_num7 << std::endl;
            EventFile8 << legends.at(i) << ": " << tmp_num8 << std::endl;
            EventFile9 << legends.at(i) << ": " << tmp_num9 << std::endl;
            if (legends.at(i) != "Data") {
                bg_num5 += tmp_num5;
                bg_num6 += tmp_num6;
                bg_num7 += tmp_num7;
                bg_num8 += tmp_num8;
                bg_num9 += tmp_num9;
            }
            tmp_num5 = 0;
            tmp_num6 = 0;
            tmp_num7 = 0;
            tmp_num8 = 0;
            tmp_num9 = 0;
        }
    }
    EventFile5 << "Total background: " << bg_num5 << std::endl;
    EventFile5.close();
    EventFile6 << "Total background: " << bg_num6 << std::endl;
    EventFile6.close();
    EventFile7 << "Total background: " << bg_num7 << std::endl;
    EventFile7.close();
    EventFile8 << "Total background: " << bg_num8 << std::endl;
    EventFile8.close();
    EventFile9 << "Total background: " << bg_num9 << std::endl;
    EventFile9.close();
    std::cout << "\nEvent yields saved in " << outdir << std::endl;
}

double Plotter13TeV::CalcXSec(std::vector<TString> datasetVec, std::vector<TString> pseudoDatasetVec, double InclusiveXsectionVec[4], double InclusiveXsectionStatErrorVec[4], TString Systematic, TString Shift) {

    double NrOfEvts_VisGen_afterSelection_noweight = 0, NrOfEvts_VisGen_afterSelection = 0;
    double NrOfEvts_afterSelection_noweight = 0, NrOfEvts_afterSelection = 0;
    double NrOfEvts_Gen_afterRecoSelection_noweight = 0, NrOfEvts_Gen_afterRecoSelection = 0;
    double NrOfEvts = 0;

    TH1D* numhists[hists.size()];
    double numbers[5] = { 0., 0., 0., 0., 0. };//[0]=data, [1]=Signal, [2]Signal(only lumi & PU weights), [3]ttbar background, [4]background(non-ttbar)
    double error_numbers[5] = { 0., 0., 0., 0., 0. };//Square of error: [0]=data, [1]=Signal, [2]Signal(only lumi & PU weights), [3]ttbar background, [4]background(non-ttbar)
    //     double TTbarBGnum =0;

    if (Systematic.Contains("UP")) { Shift = "Up"; }
    if (Systematic.Contains("DOWN")) { Shift = "Down"; }

    TString histoForInclXSec;
    if (usePoissonSmearedPseudoData) histoForInclXSec = "basic_jet_multiplicity_step8";// you might want to use: "events_weighted_step8" or "HypjetMultiXSec" or even "..._step7"
    else histoForInclXSec = "events_weighted_step8"; //check FIXME below, before changing it 


    for (unsigned int i = 0; i < datasetVec.size(); i++) {
        TH1D* hist = fileReader->GetClone<TH1D>(datasetVec[i], histoForInclXSec);
        usefulTools13TeV->ApplyFlatWeights(hist, usefulTools13TeV->CalcLumiWeight(datasetVec.at(i)));
        removeNegativeExpectationsInBinsDouble(hist); //FIXME: fine to do this, as long as one-bin histograms is used (like "events_weighted"), otherwise an integral should be checked instead       
        numhists[i] = hist;

    }

    if (usePoissonSmearedPseudoData) {

        TH1D* numhistsForPseudoData[histsForPseudoData.size()];

        for (unsigned int i = 0; i < pseudoDatasetVec.size(); i++) {
            TH1D* histForPseudoData = fileReader->GetClone<TH1D>(pseudoDatasetVec[i], histoForInclXSec);
            usefulTools13TeV->ApplyFlatWeights(histForPseudoData, usefulTools13TeV->CalcLumiWeight(pseudoDatasetVec.at(i)));
            removeNegativeExpectationsInBinsDouble(histForPseudoData);
            numhistsForPseudoData[i] = histForPseudoData;

        }

        if (channelType != 3) {// for ee, emu or mumu channel
            for (size_t i = 1; i < pseudoDatasetVec.size(); ++i) {
                if (!pseudoDatasetVec.at(i).Contains("ttbarsignalplustau") && !pseudoDatasetVec.at(i).Contains("ttbarbgviatau") && !pseudoDatasetVec.at(i).Contains("run")) {
                    numhists[0]->Add(numhistsForPseudoData[i]);
                }
            }
            Plotter13TeV::performPoissonSmearing(numhists[0]);
        }
        else if (channelType == 3) {// for combined channel
            //Indeces must match to the ordering of data files in filelist {syst}_combined.txt
            TString chToInclude[3] = { "ee_", "emu_", "mumu_" };
            for (int ch_it = 0; ch_it < 3; ++ch_it) {
                for (size_t i = ch_it + 1; i < pseudoDatasetVec.size(); ++i) { //skipping pseudo-data files (using +1) to which MC is added
                    if (!pseudoDatasetVec.at(i).Contains("ttbarsignalplustau") && !pseudoDatasetVec.at(i).Contains("ttbarbgviatau") && !pseudoDatasetVec.at(i).Contains("run") && pseudoDatasetVec.at(i).Contains(chToInclude[ch_it])) {
                        numhists[ch_it]->Add(numhistsForPseudoData[i]);
                    }
                }
                Plotter13TeV::performPoissonSmearing(numhists[ch_it]);
            }
        }
    }

    for (unsigned int i = 0; i < hists.size(); i++) { // prepare histos and leg 
        if (legends.at(i) == "Data") {
            if (usePoissonSmearedPseudoData) {
                numbers[0] += numhists[i]->Integral(0, numhists[i]->GetNbinsX() + 1);// includes under/overflow bins
                for (int j = 0; j <= numhists[i]->GetNbinsX() + 1; ++j) {
                    error_numbers[0] += TMath::Sqrt(numhists[i]->GetBinContent(j)) * TMath::Sqrt(numhists[i]->GetBinContent(j));//Not GetBinError() since MC stat will be used, but we treat "Data" 
                }
            }
            else {
                numbers[0] += numhists[i]->Integral();
                error_numbers[0] += numhists[i]->GetBinError(2) * numhists[i]->GetBinError(2); //This bin selection is hardcoded please change it if changes when filling in Analysis.C
            }
        }
        else if (legends.at(i) == "t#bar{t} signal") {
            if (usePoissonSmearedPseudoData) {
                numbers[1] += numhists[i]->Integral(0, numhists[i]->GetNbinsX() + 1);// includes under/overflow bins
                for (int j = 0; j <= numhists[i]->GetNbinsX() + 1; ++j) {
                    error_numbers[1] += numhists[i]->GetBinError(j) * numhists[i]->GetBinError(j);
                }
            }
            else {
                numbers[1] += numhists[i]->Integral();
                error_numbers[1] += numhists[i]->GetBinError(2) * numhists[i]->GetBinError(2); //This bin selection is hardcoded please change it if changes when filling in Analysis.C
            }

            TH1D* GenPlot = fileReader->GetClone<TH1D>(datasetVec.at(i), "GenAll");
            TH1D* GenPlot_noweight = fileReader->GetClone<TH1D>(datasetVec.at(i), "GenAll_noweight");
            TH1D* VisGenPlot = fileReader->GetClone<TH1D>(datasetVec.at(i), "VisGenAll");
            TH1D* VisGenPlot_noweight = fileReader->GetClone<TH1D>(datasetVec.at(i), "VisGenAll_noweight");
            TH1D* RecoGenPlot = fileReader->GetClone<TH1D>(datasetVec.at(i), "GenAll_RecoCuts");
            TH1D* RecoGenPlot_noweight = fileReader->GetClone<TH1D>(datasetVec.at(i), "GenAll_RecoCuts_noweight");
            TH1* h_NrOfEvts = fileReader->GetClone<TH1>(datasetVec.at(i), "weightedEvents");

            double LumiWeight = usefulTools13TeV->CalcLumiWeight(datasetVec.at(i));
            usefulTools13TeV->ApplyFlatWeights(GenPlot, LumiWeight);
            usefulTools13TeV->ApplyFlatWeights(GenPlot_noweight, LumiWeight);
            usefulTools13TeV->ApplyFlatWeights(VisGenPlot, LumiWeight);
            usefulTools13TeV->ApplyFlatWeights(VisGenPlot_noweight, LumiWeight);
            usefulTools13TeV->ApplyFlatWeights(RecoGenPlot, LumiWeight);
            usefulTools13TeV->ApplyFlatWeights(RecoGenPlot_noweight, LumiWeight);
            usefulTools13TeV->ApplyFlatWeights(h_NrOfEvts, LumiWeight);

            NrOfEvts = h_NrOfEvts->GetBinContent(1);
            NrOfEvts_afterSelection += GenPlot->Integral();
            NrOfEvts_afterSelection_noweight += GenPlot_noweight->Integral();
            NrOfEvts_VisGen_afterSelection += VisGenPlot->Integral();
            NrOfEvts_Gen_afterRecoSelection += RecoGenPlot->Integral();
            NrOfEvts_Gen_afterRecoSelection_noweight += RecoGenPlot_noweight->Integral();
            NrOfEvts_VisGen_afterSelection_noweight += VisGenPlot_noweight->Integral();

            numbers[2] += GenPlot->Integral();
            error_numbers[2] += GenPlot->GetBinError(18) * GenPlot->GetBinError(18); //This bin selection is hardcoded please change it if changes when filling in Analysis.C

        }
        else if (legends.at(i) == "t#bar{t} other") {
            if (usePoissonSmearedPseudoData) {
                numbers[3] += numhists[i]->Integral(0, numhists[i]->GetNbinsX() + 1);// includes under/overflow bins
                for (int j = 0; j <= numhists[i]->GetNbinsX() + 1; ++j) {
                    error_numbers[3] += numhists[i]->GetBinError(j) * numhists[i]->GetBinError(j);
                }
            }
            else {
                numbers[3] += numhists[i]->Integral();
                error_numbers[3] += numhists[i]->GetBinError(2) * numhists[i]->GetBinError(2); //This bin selection is hardcoded please change it if changes when filling in Analysis.C
            }
        }
        else {
            if ((legends.at(i) == DYEntry)) {
                numhists[i]->Scale(DYScale.at(channelType));
            }
            if ((legends.at(i) == DYEntry) && Systematic.Contains("DY_") && Shift == "Up") {
                numhists[i]->Scale(1.3);
            }
            if ((legends.at(i) == DYEntry) && Systematic.Contains("DY_") && Shift == "Down") {
                numhists[i]->Scale(0.7);
            }
            if (Systematic.Contains("BG_") && Shift == "Up" && legends.at(i) != "t#bar{t} other" && legends.at(i) != DYEntry) {
                numhists[i]->Scale(1.3);
            }
            if (Systematic.Contains("BG_") && Shift == "Down" && legends.at(i) != "t#bar{t} other" && legends.at(i) != DYEntry) {
                numhists[i]->Scale(0.7);
            }

            if (usePoissonSmearedPseudoData) {
                numbers[4] += numhists[i]->Integral(0, numhists[i]->GetNbinsX() + 1);// includes under/overflow bins
                for (int j = 0; j <= numhists[i]->GetNbinsX() + 1; ++j) {
                    error_numbers[4] += numhists[i]->GetBinError(j) * numhists[i]->GetBinError(j);
                }
            }
            else {
                numbers[4] += numhists[i]->Integral();
                error_numbers[4] += numhists[i]->GetBinError(2) * numhists[i]->GetBinError(2); //This bin selection is hardcoded please change it if changes when filling in Analysis.C
            }
        }
    }
    ////////////////////////////Make output for tables

    double tmp_num = 0;

    ofstream EventFile, XSecFile;
    TString outdir = ttbar::assignFolder(outpathPlots, subfolderChannel.Copy().Remove(0, 1), Systematic);
    EventFile.open(outdir.Copy() + "Events.txt");
    XSecFile.open(outdir.Copy() + "InclusiveXSec.txt");

    double bg_num = 0;
    for (unsigned int i = 0; i < hists.size(); i++) {
        if (usePoissonSmearedPseudoData) {
            tmp_num += numhists[i]->Integral(0, numhists[i]->GetNbinsX() + 1);
        }
        else tmp_num += numhists[i]->Integral();

        if (i == (hists.size() - 1)) {
            EventFile << legends.at(i) << ": " << tmp_num << std::endl;
            bg_num += tmp_num;
            tmp_num = 0;
        }
        else if (legends.at(i) != legends.at(i + 1)) {
            EventFile << legends.at(i) << ": " << tmp_num << std::endl;
            if (legends.at(i) != "Data")bg_num += tmp_num;
            tmp_num = 0;
        }
    }
    EventFile << "Total MCs: " << bg_num << std::endl;
    EventFile << "\nDataEvents= " << numbers[0] << std::endl;
    EventFile << "SignalReco= " << numbers[1] << std::endl;
    EventFile << "SignalGen = " << numbers[2] << std::endl;
    EventFile << "ttbar bags= " << numbers[3] << std::endl;
    EventFile << "All Backgd= " << numbers[4] << std::endl;
    EventFile << "Efficiency= " << (numbers[1] / numbers[2]) << std::endl;
    EventFile << "BrancRatio= " << BranchingFraction[channelType] << std::endl;
    EventFile << "Total Gen Events (no weights)= " << NrOfEvts << std::endl;
    EventFile << "Gen Events after Selection (no weights)= " << NrOfEvts_afterSelection_noweight << std::endl;
    EventFile << "Visible Gen Events after Selection (no weights)= " << NrOfEvts_VisGen_afterSelection_noweight << std::endl;
    EventFile << "Acceptance= " << NrOfEvts_afterSelection_noweight / NrOfEvts << std::endl;
    EventFile << "Visible Acceptance= " << NrOfEvts_VisGen_afterSelection_noweight / NrOfEvts << std::endl;
    EventFile << "------------------------------------------------------------------------------" << std::endl;
    EventFile << "Efficiency and acceptance definitions as proposed by the TopXSection conveners\n" << std::endl;
    EventFile << "N_rec = " << numbers[1] << std::endl;
    EventFile << "N_gen (with cuts at parton level) = " << NrOfEvts_VisGen_afterSelection << std::endl;
    EventFile << "N_gen (with cuts at parton level, no weights) = " << NrOfEvts_VisGen_afterSelection_noweight << std::endl;
    EventFile << "N_gen (with cuts at reco level) = " << NrOfEvts_Gen_afterRecoSelection << std::endl;
    EventFile << "N_gen (with cuts at reco level, no weights) = " << NrOfEvts_Gen_afterRecoSelection_noweight << std::endl;
    EventFile << "N_gen = " << numbers[2] << std::endl;
    EventFile << "\nEfficiency = N_rec / N_gen (with cuts at parton level) = " << numbers[1] / NrOfEvts_VisGen_afterSelection << std::endl;
    EventFile << "Efficiency = N_rec / N_gen (with cuts at parton level && noweights) = " << numbers[1] / NrOfEvts_VisGen_afterSelection_noweight << std::endl;
    EventFile << "\nEfficiency' = N_rec / N_gen (with cuts at reco level) = " << numbers[1] / NrOfEvts_Gen_afterRecoSelection << std::endl;
    EventFile << "Efficiency' = N_rec / N_gen (with cuts at reco level && noweights) = " << numbers[1] / NrOfEvts_Gen_afterRecoSelection_noweight << std::endl;
    EventFile << "\nAcceptance = N_gen (with cuts at parton level) / N_gen = " << NrOfEvts_VisGen_afterSelection / numbers[2] << std::endl;
    EventFile << "Acceptance = N_gen (with cuts at parton level && noweights) / N_gen = " << NrOfEvts_VisGen_afterSelection_noweight / numbers[2] << std::endl;
    EventFile << "Eff * Acc = " << numbers[1] / numbers[2] << std::endl;
    EventFile << "------------------------------------------------------------------------------" << std::endl;

    // Acceptance driven correction needed for avoiding dependence from N^gen_dilepton/N^gen_all(inclusive) ratio specific for each MC sample dataset while measure with pseudo-data based on ttbar signal
    double correctionForInclXSecWhilePseudoData = 1.;
    if (usePoissonSmearedPseudoData) {
        correctionForInclXSecWhilePseudoData = BranchingFraction[channelType] / (NrOfEvts_afterSelection_noweight / NrOfEvts);
    }

    double xsec = ((numbers[0] - numbers[4]) * (numbers[1] / (numbers[1] + numbers[3]))) / ((numbers[1] / numbers[2]) * BranchingFraction[channelType] * lumi);
    double xsecstaterror = TMath::Sqrt(error_numbers[0]) * (numbers[1] / (numbers[1] + numbers[3])) / ((numbers[1] / numbers[2]) * BranchingFraction[channelType] * lumi);

    if (usePoissonSmearedPseudoData) {
        xsec = correctionForInclXSecWhilePseudoData * xsec;
        xsecstaterror = correctionForInclXSecWhilePseudoData * xsecstaterror;
    }

    if (channelType != 3) {
        InclusiveXsectionVec[channelType] = xsec;
        InclusiveXsectionStatErrorVec[channelType] = xsecstaterror;
    }
    else {
        TString eefilename = "Plots/" + Systematic + "/ee/InclusiveXSec.txt";
        TString mumufilename = "Plots/" + Systematic + "/mumu/InclusiveXSec.txt";
        TString emufilename = "Plots/" + Systematic + "/emu/InclusiveXSec.txt";

        //check the existence of the file
        if (gSystem->AccessPathName(eefilename) || gSystem->AccessPathName(emufilename) || gSystem->AccessPathName(mumufilename)) {
            std::cout << "WARNING (in CalcXSec)!!" << std::endl;
            std::cout << "One of the input files you use for the combined XSection measurement doesn't exist!!\nExiting!!" << std::endl;
            exit(888);
        }

        ifstream ResultsEE(eefilename);
        ifstream ResultsEMu(emufilename);
        ifstream ResultsMuMu(mumufilename);

        TString Dummy = "";

        ResultsEE >> Dummy >> Dummy >> Dummy >> Dummy >> Dummy >> InclusiveXsectionVec[0] >> Dummy >> InclusiveXsectionStatErrorVec[0];
        ResultsMuMu >> Dummy >> Dummy >> Dummy >> Dummy >> Dummy >> InclusiveXsectionVec[1] >> Dummy >> InclusiveXsectionStatErrorVec[1];
        ResultsEMu >> Dummy >> Dummy >> Dummy >> Dummy >> Dummy >> InclusiveXsectionVec[2] >> Dummy >> InclusiveXsectionStatErrorVec[2];

        ResultsEE.close(); ResultsEMu.close(); ResultsMuMu.close();

        InclusiveXsectionVec[channelType] = (InclusiveXsectionVec[0] / (InclusiveXsectionStatErrorVec[0] * InclusiveXsectionStatErrorVec[0])
            + InclusiveXsectionVec[1] / (InclusiveXsectionStatErrorVec[1] * InclusiveXsectionStatErrorVec[1])
            + InclusiveXsectionVec[2] / (InclusiveXsectionStatErrorVec[2] * InclusiveXsectionStatErrorVec[2])
            ) /
            (1 / (InclusiveXsectionStatErrorVec[0] * InclusiveXsectionStatErrorVec[0])
                + 1 / (InclusiveXsectionStatErrorVec[1] * InclusiveXsectionStatErrorVec[1])
                + 1 / (InclusiveXsectionStatErrorVec[2] * InclusiveXsectionStatErrorVec[2])
                );

        InclusiveXsectionStatErrorVec[channelType] = 1 / (TMath::Sqrt(
            (1 / (InclusiveXsectionStatErrorVec[0] * InclusiveXsectionStatErrorVec[0]))
            + (1 / (InclusiveXsectionStatErrorVec[1] * InclusiveXsectionStatErrorVec[1]))
            + (1 / (InclusiveXsectionStatErrorVec[2] * InclusiveXsectionStatErrorVec[2]))
        ));
    }

    EventFile << "XSection  = " << InclusiveXsectionVec[channelType] << std::endl;
    EventFile << "XSecStaErr= " << InclusiveXsectionStatErrorVec[channelType] << std::endl;
    EventFile.close();
    XSecFile << "Systematic: " << Systematic << " Channel: " << subfolderChannel << " InclXSection: " << InclusiveXsectionVec[channelType] << " AbsStatError: " << InclusiveXsectionStatErrorVec[channelType] << std::endl;
    XSecFile.close();
    std::cout << "\nInclusive XSection information saved in: " << outdir << std::endl;
    return xsec;
}

int Plotter13TeV::CalcDiffXSec(TString Channel, TString Systematic, const TString& unfoldingType/* = "svd"*/) {

    double Xbins[XAxisbins.size()];
    double binWidth[XAxisbinCenters.size()];
    for (unsigned int i = 0; i < XAxisbins.size();i++) { Xbins[i] = XAxisbins[i]; }
    double GenSignalSum[XAxisbinCenters.size()];

    TString ftemp = "preunfolded/" + Systematic + "/" + Channel + "/" + name + "_UnfoldingHistos.root";
    if (gSystem->AccessPathName(ftemp)) {
        std::cout << "WARNING (in CalcDiffXSec)!!" << std::endl;
        std::cout << "File: " << ftemp << " doesn't exist" << std::endl;
        std::cout << "One of the input files you use for the combined XSection measurement doesn't exist!!" << std::endl;
        return -1;
    }
    TH1D* theDataHist = fileReader->GetClone<TH1D>(ftemp, "aDataHist");
    TH1D* theBgrHist = fileReader->GetClone<TH1D>(ftemp, "aBgrHist");
    TH1D* theTtBgrHist = fileReader->GetClone<TH1D>(ftemp, "aTtBgrHist");
    TH1D* theRecHist = fileReader->GetClone<TH1D>(ftemp, "aRecHist");
    TH1D* theGenHist = NULL;
    if (particleXsec) theGenHist = fileReader->GetClone<TH1D>(ftemp, "aPseudoHist");
    else theGenHist = fileReader->GetClone<TH1D>(ftemp, "aGenHist");
    TH2D* theRespHist = fileReader->GetClone<TH2D>(ftemp, "aRespHist");

    std::unique_ptr<TH1> theGenHistRebinned{ theGenHist->Rebin(bins,"aDataHist",Xbins) };
    for (Int_t bin = 0; bin < bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
        GenSignalSum[bin] = theGenHistRebinned->GetBinContent(bin + 1);
    }

    double GenDiffXSecVec[4][bins];
    double DiffXSecVec[4][bins];
    double DiffXSecStatErrorVec[4][bins];

    if (Channel == "ee") { channelType = 0; }
    if (Channel == "mumu") { channelType = 1; }
    if (Channel == "emu") { channelType = 2; }
    if (Channel == "combined") { channelType = 3; }

    if (doUnfolding == true && (channelType != 3 || performUnfoldingInCombinedChannel))//do the unfolding only in the individual channels: ee, emu, mumu; for "combined" channel dedicated option should be propagated
    {
        // OZ 4.01.2017 modified to support both SVD and TUnfold
        double UnfoldingResult[XAxisbinCenters.size()];
        double UnfoldingError[XAxisbinCenters.size()];
        if (unfoldingType == "svd")
        {

            // SVD Helper Class
            DilepSVDFunctions mySVDFunctions;
            mySVDFunctions.SetOutputPath(outpath);

            // Binning
            double* theBins = Xbins;
            int numberBins = bins;

            // Names and Labels
            TString channelLabelStr(channelLabel[channelType]);
            TString theChannelName = Channel;
            TString theParticleName = "";
            if (name.Contains("Lepton")) {
                theParticleName = "Leptons";
                if (revokeAntiQuantityCombination && !name.Contains("Anti")) theParticleName = "Lepton";
                if (revokeAntiQuantityCombination && name.Contains("Anti")) theParticleName = "AntiLepton";
            }
            if (name.Contains("LLBar")) theParticleName = "LepPair";
            if (name.Contains("Top")) {
                theParticleName = "TopQuarks";
                if (revokeAntiQuantityCombination && !name.Contains("Anti")) theParticleName = "TopQuark";
                if (revokeAntiQuantityCombination && name.Contains("Anti")) theParticleName = "TopAntiQuark";
            }
            if (name.Contains("TTBar")) theParticleName = "TtBar";
            if (name.Contains("BBBar")) theParticleName = "BBbar";
            if (name.Contains("BJet")) {
                theParticleName = "BJets";
                if (revokeAntiQuantityCombination && !name.Contains("Anti")) theParticleName = "BJet";
                if (revokeAntiQuantityCombination && name.Contains("Anti")) theParticleName = "AntiBJet";
            }
            if (name.Contains("JetMult")) theParticleName = "Jets";
            if (name.Contains("ExtraJet")) theParticleName = "ExtraJets";
            if (name.Contains("DeltaRJet12")) theParticleName = "DeltaR";
            TString theQuantityName = "";
            if (name.Contains("pT")) theQuantityName = "Pt";
            if (name.Contains("Eta")) theQuantityName = "Eta";
            if (name.Contains("Rapidity")) theQuantityName = "Rapidity";
            if (name.Contains("Mass")) theQuantityName = "Mass";
            if (name.Contains("JetMult")) theQuantityName = "Mult";
            if (name.Contains("ExtraJetEta")) theQuantityName = "Eta";
            if (name.Contains("ExtraJetpT")) theQuantityName = "Pt";
            if (name.Contains("DeltaRJet12")) theQuantityName = "DeltaR";
            TString theSpecialPostfix = "";
            //if (name.Contains("Lead")) theSpecialPostfix = name;
            theSpecialPostfix = name;
            if (specialComment.CompareTo("Standard") != 0) {
                //theSpecialPostfix = specialComment;
            }

            double totalDataEventsNom[1] = { TopSVDFunctions::SVD_Integral1D(theDataHist, 0, false) };
            double totalBgrEventsNom[1] = { TopSVDFunctions::SVD_Integral1D(theBgrHist, 0, false) };
            double totalTtBgrEventsNom[1] = { TopSVDFunctions::SVD_Integral1D(theTtBgrHist, 0, false) };
            double totalRecEventsNom[1] = { TopSVDFunctions::SVD_Integral1D(theRecHist, 0, false) };
            double totalGenEventsNom[1] = { TopSVDFunctions::SVD_Integral1D(theGenHist, 0, true) };

            // UNFOLDING 
            // Retrieve a histogram with the unfolded quantities.
            // Note: The unfolded histogram has additional side bins!
            // Keep this in mind when accessing bin content via indices
            TH1D* unfoldedDistribution = NULL;
            TH1D* unfoldedDistributionNormalized = NULL;
            int numSystematics = 0;
            mySVDFunctions.SVD_DoUnfold(
                theDataHist,
                theBgrHist,
                theTtBgrHist,
                theGenHist,
                theRecHist,
                theRespHist,
                totalDataEventsNom,
                totalBgrEventsNom,
                totalTtBgrEventsNom,
                totalRecEventsNom,
                totalGenEventsNom,
                theBins, numberBins,
                unfoldedDistribution,
                unfoldedDistributionNormalized,
                numSystematics,
                theChannelName, theParticleName, theQuantityName, theSpecialPostfix, "");

            // Make a vector from the result
            for (size_t i = 0; i < XAxisbinCenters.size(); i++) {
                UnfoldingResult[i] = unfoldedDistributionNormalized->GetBinContent(i + 2);//account for extra row in SVD unfolding
                UnfoldingError[i] = unfoldedDistributionNormalized->GetBinError(i + 2); //absolute bin error
                if (absoluteXsec) {
                    UnfoldingResult[i] = unfoldedDistribution->GetBinContent(i + 2) / lumi;//account for extra row in SVD unfolding
                    UnfoldingError[i] = unfoldedDistribution->GetBinError(i + 2) / lumi;
                }
                //UnfoldingResult[i] = unfoldedDistribution->GetBinContent(i+2);//account for extra row in SVD unfolding
                //UnfoldingError[i] = unfoldedDistribution->GetBinError(i+2);
            }
        } // end SVD
        else if (unfoldingType == "tunfold") // OZ 4.01.2017
        {
            int debug = 0; // more printouts 
            int doBBB = 0; // do bin-by-bin unfolding for comparison
            int doSimpleTUnfold = 0; // do simplified unregularizefd unfolding for comparison

            // create plotter
            Samples samples;
            PlotterForTUnfold plotter(samples, lumi, 0.0);
            std::vector<TString> v_plotName;
            v_plotName.push_back("1d");
            v_plotName.push_back(name);
            plotter.setOptions(v_plotName);
            plotter.SetPlotsFolder(Channel::convert(Channel), Systematic::Systematic(Systematic));

            // set generator distributions (needed for comparison)
            TH2D* hPlotGen = PlotterForTUnfold::H1dto2d((TH1D*)(theGenHistRebinned.get()));
            plotter.SetHistGen(hPlotGen);

            // set regularisation method (set PlotterForTunfold.cc for more details)
            TString regMethod = "rho";
            //TString regMethod = "0";
            plotter.setRegMethod(regMethod);

            // see PlotterForTunfold.cc if you need options below
            //plotter.setTauTest(runTauTest);
            //plotter.setCT(runCT);
            //if(var!=-999)plotter.setVar(var);

            // rebin response matrix, adding events which are generated but not reconstructed
            std::vector<double> fineBins = plotter.GetFineBinning()[1];
            TH2D* theRespHist_TUnfold = new TH2D("theRespHist_TUnfold", "theRespHist_TUnfold", bins, Xbins, fineBins.size() - 1, &fineBins[0]);
            for (int i = 0; i <= theRespHist->GetNbinsX() + 1; i++)
            {
                // determine bin index in re-binned histogram
                double y = theRespHist->GetXaxis()->GetBinCenter(i);
                int reBinY = theRespHist_TUnfold->GetYaxis()->FindBin(y);
                // check that this bin is fully inside one of the bins in re-binned histogram
                double yL = theRespHist->GetXaxis()->GetBinLowEdge(i);
                double yR = theRespHist->GetXaxis()->GetBinUpEdge(i);
                if (reBinY != 0 && reBinY != theRespHist_TUnfold->GetYaxis()->GetNbins() && ((theRespHist_TUnfold->GetYaxis()->FindBin(yL) != reBinY && !PlotterForTUnfold::IsEqual(yL, theRespHist_TUnfold->GetYaxis()->GetBinLowEdge(reBinY))) || (theRespHist_TUnfold->GetYaxis()->FindBin(yR) != reBinY && !PlotterForTUnfold::IsEqual(yR, theRespHist_TUnfold->GetYaxis()->GetBinUpEdge(reBinY)))))
                {
                    printf("Error in Plotter13TeV::CalcDiffXSec(): can not rebin 2D responce histogram, incompatible binning scheme for detector level\n");
                    printf("input bin [%f %f] -> [%f %f]\n", yL, yR, theRespHist_TUnfold->GetYaxis()->GetBinLowEdge(reBinY), theRespHist_TUnfold->GetYaxis()->GetBinUpEdge(reBinY));
                    exit(1);
                }
                for (int j = 0; j <= theRespHist->GetNbinsY() + 1; j++)
                {
                    // determine bin index in re-binned histogram
                    double x = theRespHist->GetYaxis()->GetBinCenter(j);
                    int reBinX = theRespHist_TUnfold->GetXaxis()->FindBin(x);
                    // check that this bin is fully inside one of the bins in re-binned histogram
                    if (i == 1)
                    {
                        double xL = theRespHist->GetYaxis()->GetBinLowEdge(j);
                        double xR = theRespHist->GetYaxis()->GetBinUpEdge(j);
                        if (reBinX != 0 && reBinX != theRespHist_TUnfold->GetXaxis()->GetNbins() && ((theRespHist_TUnfold->GetXaxis()->FindBin(xL) != reBinX && !PlotterForTUnfold::IsEqual(xL, theRespHist_TUnfold->GetXaxis()->GetBinLowEdge(reBinX))) || (theRespHist_TUnfold->GetXaxis()->FindBin(xR) != reBinX && !PlotterForTUnfold::IsEqual(xR, theRespHist_TUnfold->GetXaxis()->GetBinUpEdge(reBinX)))))
                        {
                            printf("Error in Plotter13TeV::CalcDiffXSec(): can not rebin 2D responce histogram, incompatible binning scheme for generator level\n");
                            printf("input bin [%f %f] -> [%f %f]\n", xL, xR, theRespHist_TUnfold->GetXaxis()->GetBinLowEdge(reBinX), theRespHist_TUnfold->GetXaxis()->GetBinUpEdge(reBinX));
                            exit(1);
                        }
                    }
                    // bin content and unc from input histogram
                    double thisBinContent = theRespHist->GetBinContent(i, j);
                    double thisBinUnc = theRespHist->GetBinError(i, j);
                    // existing bin content and unc in re-binned histogram
                    double prevBinConent = theRespHist_TUnfold->GetBinContent(reBinX, reBinY);
                    double prevBinUnc = theRespHist_TUnfold->GetBinError(reBinX, reBinY);
                    // calculate updated bin content and unc in re-binned histogram
                    double newBinContent = thisBinContent + prevBinConent;
                    double newBinUnc = sqrt(thisBinUnc * thisBinUnc + prevBinUnc * prevBinUnc);
                    theRespHist_TUnfold->SetBinContent(reBinX, reBinY, newBinContent);
                    theRespHist_TUnfold->SetBinError(reBinX, reBinY, newBinUnc);
                    //theRespHist_TUnfold->SetBinError(reBinX, reBinY, 0.0);
                }
            }
            // add generated but not reconstructed events to underflow bins of response matrix
            std::unique_ptr<TH1D> hRecCoarse(new TH1D("hRecCoarse", "hRecCoarse", bins, Xbins));
            TH1D* biniHist = theRespHist_TUnfold->ProjectionY(); // this is reco level fine binning
            // discard original content of under/overflow bin on both rec and gen level
            theRespHist_TUnfold->ClearUnderflowAndOverflow();
            biniHist->ClearUnderflowAndOverflow();
            theGenHistRebinned->ClearUnderflowAndOverflow();
            for (int i = 0; i <= biniHist->GetNbinsX() + 1; i++)
                hRecCoarse->Fill(biniHist->GetXaxis()->GetBinCenter(i), biniHist->GetBinContent(i));
            for (int bx = 0; bx <= theRespHist_TUnfold->GetNbinsX() + 1; bx++)
            {
                double reco = 0.0;
                for (int by = 0; by <= theRespHist_TUnfold->GetNbinsY() + 1; by++)
                {
                    reco += theRespHist_TUnfold->GetBinContent(bx, by);
                }
                double gen = theGenHistRebinned->GetBinContent(bx);
                theRespHist_TUnfold->SetBinContent(bx, 0, gen - reco);
                theRespHist_TUnfold->SetBinError(bx, 0, 0.0); // uncertainty is neglected
            }

            // background handling
            TH1D* theDataHistNoBkgr = new TH1D(*theDataHist); // this will be background free spectrum to unfold
            bool beForgiving = false;
            int numHist = 1;
            //theDataHistNoBkgr->Print("all");
            //theBgrHist->Print("all");
            //theTtBgrHist->Print("all");
            //biniHist->Print("all");
            //theDataHist->Print("all");
            // discard original content of under/overflow bin
            theDataHistNoBkgr->ClearUnderflowAndOverflow();
            theDataHist->ClearUnderflowAndOverflow();
            theBgrHist->ClearUnderflowAndOverflow();
            theTtBgrHist->ClearUnderflowAndOverflow();
            biniHist->ClearUnderflowAndOverflow();
            TopSVDFunctions::SVD_BackgrHandling(theDataHistNoBkgr, theBgrHist, theTtBgrHist, biniHist, theDataHist, beForgiving, numHist);
            if (debug)
            {
                printf("Total number of data events: %f with uo %f\n", theDataHist->Integral(), theDataHist->Integral(0, theDataHist->GetNbinsX() + 1));
                printf("Total number of data events no background : %f with uo %f\n", theDataHistNoBkgr->Integral(), theDataHistNoBkgr->Integral(0, theDataHistNoBkgr->GetNbinsX() + 1));
                printf("Total number of background events: %f with uo %f\n", theBgrHist->Integral(), theBgrHist->Integral(0, theBgrHist->GetNbinsX() + 1));
                printf("Total number of tt background events: %f with uo %f\n", theTtBgrHist->Integral(), theTtBgrHist->Integral(0, theTtBgrHist->GetNbinsX() + 1));
                printf("Total number of reconstructed events: %f with uo %f\n", theRespHist_TUnfold->Integral(0, theRespHist_TUnfold->GetNbinsX() + 1, 1, theRespHist_TUnfold->GetNbinsY()), theRespHist_TUnfold->Integral(0, theRespHist_TUnfold->GetNbinsX() + 1, 0, theRespHist_TUnfold->GetNbinsY() + 1));
                printf("Total number of generated events: %f with uo %f\n", theRespHist_TUnfold->Integral(1, theRespHist_TUnfold->GetNbinsX(), 0, theRespHist_TUnfold->GetNbinsY() + 1), theRespHist_TUnfold->Integral(0, theRespHist_TUnfold->GetNbinsX() + 1, 0, theRespHist_TUnfold->GetNbinsY() + 1));
            }

            // do TUnfold
            plotter.RunUnfolding(theRespHist_TUnfold, theDataHistNoBkgr, NULL, NULL);
            TH2D* histUf2d = (TH2D*)(absoluteXsec ? plotter.GetUnfoldedData() : plotter.GetUnfoldedDataNorm()); // this histogram is hold by the plotter
            TH1D* histUf = PlotterForTUnfold::H2dto1d(histUf2d);
            if (debug)
            {
                printf("TUnfold done, printing output histogram:\n");
                histUf->Print("all");
            }

            // the code below does TUnfold, but simpler than it is done above 
            // (e.g. no regularisation; see PlotterForTUnfold::DoSimpleUnfolding() for more details), 
            // therefore could be useful for cross checks (but gives no plots etc.)
            if (doSimpleTUnfold)
            {
                TH1* histUfOZ = PlotterForTUnfold::DoSimpleUnfolding(theDataHistNoBkgr, theRespHist_TUnfold);
                if (debug)
                {
                    printf("PlotterForTUnfold::DoSimpleUnfolding() done, printing output histogram:\n");
                    histUfOZ->Print("all");
                }
            }

            // need to adjust channel numbering
            int channelForTUnfold = channelType;
            if (channelForTUnfold == 2)
                channelForTUnfold = 1;
            else if (channelForTUnfold == 1)
                channelForTUnfold = 2;
            plotter.writePlotXSec(channelForTUnfold);
            //plotter.writePlotEPSAllBins(); // anyway purity-stability-efficiency plots will be empty, get them from SVD

            // make a vector from the result
            for (size_t i = 0; i < XAxisbinCenters.size(); i++)
            {
                UnfoldingResult[i] = histUf->GetBinContent(i + 1);
                UnfoldingError[i] = histUf->GetBinError(i + 1);
                if (absoluteXsec)
                {
                    UnfoldingResult[i] = histUf->GetBinContent(i + 1) / lumi;
                    UnfoldingError[i] = histUf->GetBinError(i + 1) / lumi;
                }
            }

            delete hPlotGen;
            delete theRespHist_TUnfold;
            delete theDataHistNoBkgr;
            delete biniHist;
            delete histUf;

            // as a cross check, do bin-by-bin unfolding (denoted BBB)
            if (doBBB)
            {
                // rebin response matrix (correct uncertainties are not needed)
                std::vector<double> fineBins = plotter.GetFineBinning()[0];
                std::unique_ptr<TH2D> theRespHist_BBB(new TH2D("theRespHist_BBB", "theRespHist_BBB", bins, Xbins, bins, Xbins));
                for (int i = 0; i <= theRespHist->GetNbinsX() + 1; i++)
                    for (int j = 0; j <= theRespHist->GetNbinsY() + 1; j++)
                        theRespHist_BBB->Fill(theRespHist->GetYaxis()->GetBinCenter(j), theRespHist->GetXaxis()->GetBinCenter(i), theRespHist->GetBinContent(i, j));

                // rebin data and background matrices
                std::unique_ptr<TH1D> theDataHistRebinned((TH1D*)theDataHist->Rebin(bins, "aDataHist", Xbins));
                std::unique_ptr<TH1D> theBgrHistRebinned((TH1D*)theBgrHist->Rebin(bins, "theBgrHist", Xbins));
                std::unique_ptr<TH1D> theTtBgrHistRebinned((TH1D*)theTtBgrHist->Rebin(bins, "theTtBgrHist", Xbins));

                // background handling
                std::unique_ptr<TH1D> theDataHistNoBkgr_BBB(new TH1D(*(static_cast<TH1D*>(theDataHistRebinned.get()))));
                TH1D* theDataHistNoBkgr_BBB_Ptr = theDataHistNoBkgr_BBB.get();
                std::unique_ptr<TH1D> biniHist_BBB{ (TH1D*)theRespHist_BBB->ProjectionY() };
                bool beForgiving = false;
                int numHist = 1;
                TopSVDFunctions::SVD_BackgrHandling(theDataHistNoBkgr_BBB_Ptr, theBgrHistRebinned.get(), theTtBgrHistRebinned.get(), biniHist_BBB.get(), theDataHistRebinned.get(), beForgiving, numHist);
                if (debug)
                {
                    printf("BBB Total number of data events: %f with uo %f\n", theDataHistNoBkgr_BBB->Integral(), theDataHistNoBkgr_BBB->Integral(0, theDataHistNoBkgr_BBB->GetNbinsX() + 1));
                    printf("BBB Total number of background events: %f with uo %f\n", theBgrHistRebinned->Integral(), theBgrHistRebinned->Integral(0, theBgrHistRebinned->GetNbinsX() + 1));
                    printf("BBB Total number of tt background events: %f with uo %f\n", theTtBgrHistRebinned->Integral(), theTtBgrHistRebinned->Integral(0, theTtBgrHistRebinned->GetNbinsX() + 1));
                    printf("BBB Total number of reconstructed events: %f with uo %f\n", theRespHist_BBB->Integral(0, theRespHist_BBB->GetNbinsX() + 1, 1, theRespHist_BBB->GetNbinsY()), theRespHist_BBB->Integral(0, theRespHist_BBB->GetNbinsX() + 1, 0, theRespHist_BBB->GetNbinsY() + 1));
                    //printf("BBB Total number of generated events: %f with uo %f\n", theRespHist_BBB->Integral(), theRespHist_BBB->Integral(0, theRespHist_BBB->GetNbinsX() + 1, 1, theRespHist_BBB->GetNbinsY()));
                }

                // calculate BBB (takes into account only data stat. uncertainties)
                std::unique_ptr<TH1> bbb(new TH1D("bbb", "bbb", bins, Xbins));
                for (int b = 1; b <= theDataHistNoBkgr_BBB->GetNbinsX(); b++)
                {
                    double data = theDataHistNoBkgr_BBB->GetBinContent(b);
                    double reco = 0.0;
                    for (int bb = 0; bb <= theRespHist_BBB->GetNbinsX() + 1; bb++)
                        reco += theRespHist_BBB->GetBinContent(bb, b);
                    double gen = theGenHistRebinned->GetBinContent(b);
                    double eff = reco / gen;
                    //printf("BBB[%d] data = %f reco = %f gen = %f eff = %f\n", b, data, reco, gen, eff);
                    //bbb->SetBinContent(b, data / eff / (bbb->GetBinLowEdge(b+1) - bbb->GetBinLowEdge(b))); // if needed to divide by bin widths
                    bbb->SetBinContent(b, data / eff);
                    bbb->SetBinError(b, bbb->GetBinContent(b) / TMath::Sqrt(theDataHistNoBkgr_BBB->GetBinContent(b)));
                }
                bbb->Scale(1.0 / bbb->Integral());
                //bbb->Scale(1.0 / bbb->Integral("width"));
                if (debug)
                {
                    printf("BBB unfolding done, printing output histogram:\n");
                    bbb->Print("all");
                    // print BBB results to file
                    gSystem->mkdir("UnfoldingResults/" + Systematic + "/" + Channel, true);
                    TString ResultsFilestring = TString::Format("%s%s/%s/%s", outpathResults.Data(), subfolderSpecial.Data(), Systematic.Data(), Channel.Data());
                    if (channelType == 3 && performUnfoldingInCombinedChannel) ResultsFilestring.Append("Unfolded");
                    ResultsFilestring.Append("/" + name + "ResultsBBBTUnfold.txt");
                    ofstream ResultsFile;
                    ResultsFile.open(ResultsFilestring.Data());
                    for (int b = 1; b <= bbb->GetNbinsX(); b++)
                    {
                        double bw = bbb->GetBinLowEdge(b + 1) - bbb->GetBinLowEdge(b);
                        ResultsFile << "XAxisbinCenters[bin]: " << bbb->GetBinCenter(b) << " bin: " << bbb->GetBinLowEdge(b) << " to " << bbb->GetBinLowEdge(b + 1) << " DiffXsec: " << bbb->GetBinContent(b) / bw << " StatError: " << bbb->GetBinError(b) / bw << std::endl;
                    }
                    ResultsFile.close();
                }
            } // end BBB unfolding cross check
        } // end TUnfold
        else
        {
            printf("Error in Plotter13TeV::CalcDiffXSec(): unsupported unfolding = %s\n", unfoldingType.Data());
        }

        SignalEventswithWeight = 0;
        // CROSS SECTION CALCULATION
        for (Int_t i = 0; i < bins; ++i) {
            SignalEventswithWeight += GenSignalSum[i];
        }

        for (Int_t i = 0; i < bins; ++i) {
            //      if(Channel!="combined"){
            binWidth[i] = Xbins[i + 1] - Xbins[i];
            DiffXSecVec[channelType][i] = UnfoldingResult[i] / (binWidth[i]);
            DiffXSecStatErrorVec[channelType][i] = UnfoldingError[i] / (binWidth[i]); // absolute statistical error
            GenDiffXSecVec[channelType][i] = (GenSignalSum[i] * topxsec) / (SignalEventswithWeight * binWidth[i]);//DIRTY (signal*topxsec)/(total events*binwidth)

            if (absoluteXsec && !particleXsec) { //extrapolate non-fiducial absolute measurements accordingly to BR
                DiffXSecVec[channelType][i] = DiffXSecVec[channelType][i] / BranchingFraction[channelType];
                DiffXSecStatErrorVec[channelType][i] = DiffXSecStatErrorVec[channelType][i] / BranchingFraction[channelType];
            }

            if (name.Contains("Lepton") || name.Contains("Top") || name.Contains("BJet")) {
                GenDiffXSecVec[channelType][i] = GenDiffXSecVec[channelType][i] / 2.;
            }
        }
    }

    if (doUnfolding && channelType == 3 && !performUnfoldingInCombinedChannel) {//for 'combined' channel: do an statistical combination of the the 3 independent channels

        TString eefilename = "UnfoldingResults/" + Systematic + "/ee/" + name + "Results.txt";
        TString emufilename = "UnfoldingResults/" + Systematic + "/emu/" + name + "Results.txt";
        TString mumufilename = "UnfoldingResults/" + Systematic + "/mumu/" + name + "Results.txt";

        //check the existence of the file
        if (gSystem->AccessPathName(eefilename) || gSystem->AccessPathName(emufilename) || gSystem->AccessPathName(mumufilename)) {
            std::cout << "WARNING (in CalcDiffXSec)!!" << std::endl;
            std::cout << "One of the input files you use for the combined XSection measurement doesn't exist!!" << std::endl;
            return -1;
        }

        ifstream ResultsEE(eefilename);
        ifstream ResultsEMu(emufilename);
        ifstream ResultsMuMu(mumufilename);
        double perChannelDiffXSecPlot[3][bins];      //perChannelDiffXSecPlot[channel][bin]
        double perChannelDiffXSecStatError[3][bins]; //perChannelDiffXSecStatError[channel][bin]
        double perChannelGenDiffXSec[3][bins];       //perChannelGenDiffXSec[channel][bin]
        TString Dummy = "";
        for (Int_t bin = 0; bin < bins; bin++) {//Retrieve arrays for plotting
            ResultsEE >> Dummy >> XAxisbinCenters[bin] >> Dummy >> Xbins[bin] >> Dummy >> Xbins[bin + 1] >> Dummy >> perChannelDiffXSecPlot[0][bin] >> Dummy >> perChannelDiffXSecStatError[0][bin] >> Dummy >> perChannelGenDiffXSec[0][bin];
            ResultsMuMu >> Dummy >> XAxisbinCenters[bin] >> Dummy >> Xbins[bin] >> Dummy >> Xbins[bin + 1] >> Dummy >> perChannelDiffXSecPlot[1][bin] >> Dummy >> perChannelDiffXSecStatError[1][bin] >> Dummy >> perChannelGenDiffXSec[1][bin];
            ResultsEMu >> Dummy >> XAxisbinCenters[bin] >> Dummy >> Xbins[bin] >> Dummy >> Xbins[bin + 1] >> Dummy >> perChannelDiffXSecPlot[2][bin] >> Dummy >> perChannelDiffXSecStatError[2][bin] >> Dummy >> perChannelGenDiffXSec[2][bin];
        }
        ResultsEE.close(); ResultsEMu.close(); ResultsMuMu.close();

        //for gen level distribution
        for (Int_t i = 0; i < bins; ++i) {
            SignalEventswithWeight += GenSignalSum[i];
        }

        //do the actual combined Diff.XSection calculation
        for (int i = 0; i < bins; i++) {
            binWidth[i] = Xbins[i + 1] - Xbins[i];
            for (int j = 0; j < 3; j++) {//check if any stat error is 0, in this case set their contribution to 0!!
                if (perChannelDiffXSecStatError[j][i] == 0) {
                    perChannelDiffXSecStatError[j][i] = 1e100;
                    perChannelDiffXSecPlot[j][i] = 0;
                }
            }
            DiffXSecVec[channelType][i] = (perChannelDiffXSecPlot[0][i] / (perChannelDiffXSecStatError[0][i] * perChannelDiffXSecStatError[0][i])
                + perChannelDiffXSecPlot[1][i] / (perChannelDiffXSecStatError[1][i] * perChannelDiffXSecStatError[1][i])
                + perChannelDiffXSecPlot[2][i] / (perChannelDiffXSecStatError[2][i] * perChannelDiffXSecStatError[2][i])) /
                ((1 / (perChannelDiffXSecStatError[0][i] * perChannelDiffXSecStatError[0][i])) +
                    (1 / (perChannelDiffXSecStatError[1][i] * perChannelDiffXSecStatError[1][i])) +
                    (1 / (perChannelDiffXSecStatError[2][i] * perChannelDiffXSecStatError[2][i])));
            DiffXSecStatErrorVec[channelType][i] = 1 / (TMath::Sqrt((1 / (perChannelDiffXSecStatError[0][i] * perChannelDiffXSecStatError[0][i])) +
                (1 / (perChannelDiffXSecStatError[1][i] * perChannelDiffXSecStatError[1][i])) +
                (1 / (perChannelDiffXSecStatError[2][i] * perChannelDiffXSecStatError[2][i]))));

            if (absoluteXsec && particleXsec) {// only for fiducial cross sections, where results are not supposed to be extrapolated
                DiffXSecVec[channelType][i] = perChannelDiffXSecPlot[0][i] + perChannelDiffXSecPlot[1][i] + perChannelDiffXSecPlot[2][i]; // Just linear sum of entries w/o extrapolation with BR
                DiffXSecStatErrorVec[channelType][i] = TMath::Sqrt((perChannelDiffXSecStatError[0][i] * perChannelDiffXSecStatError[0][i]) +
                    (perChannelDiffXSecStatError[1][i] * perChannelDiffXSecStatError[1][i]) +
                    (perChannelDiffXSecStatError[2][i] * perChannelDiffXSecStatError[2][i])); // Just sum of squares w/o extrapolation with BR
            }

            GenDiffXSecVec[channelType][i] = (GenSignalSum[i] * topxsec) / (SignalEventswithWeight * binWidth[i]);//DIRTY (signal*topxsec)/(total events*binwidth)
        }
    }

    ofstream ResultsFile, ResultsLateX;

    if (channelType == 3 && performUnfoldingInCombinedChannel) gSystem->mkdir("UnfoldingResults/" + Systematic + "/" + Channel + "Unfolded", true);
    else gSystem->mkdir("UnfoldingResults/" + Systematic + "/" + Channel, true);

    string ResultsFilestring = outpathResults.Data();
    ResultsFilestring.append(subfolderSpecial.Data());
    ResultsFilestring.append("/");
    ResultsFilestring.append(Systematic);
    ResultsFilestring.append("/");
    ResultsFilestring.append(Channel);
    if (channelType == 3 && performUnfoldingInCombinedChannel) ResultsFilestring.append("Unfolded");
    ResultsFilestring.append("/");
    ResultsFilestring.append(name);
    ResultsFilestring.append("Results.txt");
    ResultsFile.open(ResultsFilestring.c_str());


    string ResultsFilestringLatex = outpathPlots.Data();
    ResultsFilestringLatex.append(subfolderChannel.Data());
    if (channelType == 3 && performUnfoldingInCombinedChannel) ResultsFilestringLatex.append("Unfolded");
    ResultsFilestringLatex.append(subfolderSpecial.Data());
    ResultsFilestringLatex.append("/");
    ResultsFilestringLatex.append(name);
    ResultsFilestringLatex.append("ResultsLaTeX.txt");
    ResultsLateX.open(ResultsFilestringLatex.c_str());
    ResultsLateX << "Bin Center & Bin & 1/#sigma d#sigma/dX & stat(%) & syst(%) & total(%)" << std::endl;
    for (Int_t bin = 0; bin < bins; bin++) {//condense matrices to arrays for plotting
        ResultsFile << "XAxisbinCenters[bin]: " << XAxisbinCenters[bin] << " bin: " << Xbins[bin] << " to " << Xbins[bin + 1] << " DiffXsec: " << DiffXSecVec[channelType][bin] << " StatError: " << DiffXSecStatErrorVec[channelType][bin] << " GenDiffXsec: " << GenDiffXSecVec[channelType][bin] << std::endl;
    }
    ResultsFile.close();
    ResultsLateX.close();

    //clean up
    delete theDataHist;
    delete theBgrHist;
    delete theTtBgrHist;
    delete theRecHist;
    delete theGenHist;
    delete theRespHist;

    Plotter13TeV::PlotSingleDiffXSec(Channel, Systematic);

    return 1;
}

void Plotter13TeV::PlotDiffXSec(TString Channel, std::vector<TString>vec_systematic) {

    setDataSet(Channel, "Nominal");
    if (!fillHisto()) return;

    if (Channel == "ee") { channelType = 0; }
    if (Channel == "mumu") { channelType = 1; }
    if (Channel == "emu") { channelType = 2; }
    if (Channel == "combined") { channelType = 3; }

    TH1::AddDirectory(kFALSE);
    TGaxis::SetMaxDigits(2);

    double Xbins[XAxisbins.size()];
    for (unsigned int i = 0; i < XAxisbins.size();i++) { Xbins[i] = XAxisbins[i]; }
    double binCenters[XAxisbinCenters.size()];
    if (drawSmoothMadgraph) {
        for (unsigned int i = 0; i < XAxisbinCenters.size();i++) {
            binCenters[i] = XAxisbinCenters[i];
        }
    }
    else {
        for (unsigned int i = 0; i < XAxisbinCenters.size();i++) {
            binCenters[i] = XAxisbins[i] + (XAxisbins[i + 1] - XAxisbins[i]) / 2;
        }
    }

    double DataSum[XAxisbinCenters.size()];
    double GenSignalSum[XAxisbinCenters.size()];
    double BGSum[XAxisbinCenters.size()];
    bool init = false;
    TH1* varhists[hists.size()];

    TString newname = name;
    if (name.Contains("Hyp")) {//Histogram naming convention has to be smarter
        newname.ReplaceAll("Hyp", 3, "", 0);
    }

    TString ftemp = "preunfolded/Nominal/" + Channel + "/" + name + "_UnfoldingHistos.root";

    //     theDataHist =  fileReader->GetClone<TH1D>(ftemp, "aDataHist");
    //     theBgrHist =  fileReader->GetClone<TH1D>(ftemp, "aBgrHist");
    //     theTtBgrHist =  fileReader->GetClone<TH1D>(ftemp, "aTtBgrHist");
    //     RecoPlot =  fileReader->GetClone<TH1D>(ftemp, "aRecHist");
    TH1* GenPlotTheory = NULL;
    if (particleXsec) GenPlotTheory = fileReader->GetClone<TH1D>(ftemp, "aPseudoHist");
    else GenPlotTheory = fileReader->GetClone<TH1D>(ftemp, "aGenHist");
    TH1* GenPlotTheoryForEffAndAccCalculation = fileReader->GetClone<TH1D>(ftemp, "aGenHist");
    TH2* genReco2d = fileReader->GetClone<TH2D>(ftemp, "aRespHist");

    for (unsigned int i = 0; i < hists.size(); i++) {
        varhists[i] = hists[i].Rebin(bins, "varhists", Xbins);
        setStyle(varhists[i], i);
    }

    if (particleXsec) {
        TH1* histsOOS[histsTtBgrOutOfSpace.size()];
        for (unsigned int j = 0; j < histsTtBgrOutOfSpace.size(); j++) {
            TString tmp = histsTtBgrOutOfSpace[j].GetName();
            histsOOS[j] = histsTtBgrOutOfSpace[j].Rebin(bins, tmp, Xbins);
            removeNegativeExpectationsInBins(histsOOS[j]);
            for (unsigned int i = 0; i < hists.size(); i++) {
                if (legends.at(i) == "t#bar{t} signal" && dataset.at(i).Contains("emu_ttbarsignalplustau") && tmp.BeginsWith("emu_")) varhists[i]->Add(histsOOS[j], -1.0);
                if (legends.at(i) == "t#bar{t} other" && dataset.at(i).Contains("emu_ttbarbgviatau") && tmp.BeginsWith("emu_")) varhists[i]->Add(histsOOS[j]);
                if (legends.at(i) == "t#bar{t} signal" && dataset.at(i).Contains("ee_ttbarsignalplustau") && tmp.BeginsWith("ee_")) varhists[i]->Add(histsOOS[j], -1.0);
                if (legends.at(i) == "t#bar{t} other" && dataset.at(i).Contains("ee_ttbarbgviatau") && tmp.BeginsWith("ee_")) varhists[i]->Add(histsOOS[j]);
                if (legends.at(i) == "t#bar{t} signal" && dataset.at(i).Contains("mumu_ttbarsignalplustau") && tmp.BeginsWith("mumu_")) varhists[i]->Add(histsOOS[j], -1.0);
                if (legends.at(i) == "t#bar{t} other" && dataset.at(i).Contains("mumu_ttbarbgviatau") && tmp.BeginsWith("mumu_")) varhists[i]->Add(histsOOS[j]);
            }
        }
    }

    std::unique_ptr<TH1> GenPlot{ GenPlotTheory->Rebin(bins,"genplot",Xbins) };
    std::unique_ptr<TH1> GenPlotForEffAndAccCalculation{ GenPlotTheoryForEffAndAccCalculation->Rebin(bins,"genplotForEffAndAcc",Xbins) };

    THStack* stack = new THStack("def", "def");
    TLegend* leg = new TLegend();
    int legchange = 0;
    std::vector<TH1*> varhistsPlotting;
    varhistsPlotting.resize(hists.size());

    for (unsigned int i = 0; i < hists.size(); i++) { // prepare histos and leg
        setStyle(varhists[i], i);
        varhistsPlotting[i] = (TH1*)varhists[i]->Clone();
        if (legends.at(i) != "Data") {
            if ((legends.at(i) == DYEntry)) {
                varhists[i]->Scale(DYScale.at(channelType));
                varhistsPlotting[i]->Scale(DYScale.at(channelType));
            }

            if (i != (hists.size() - 1)) {
                if (legends.at(i) != legends.at(i + 1)) {
                    //std::cout<<legends.at(i)<<std::endl;
                    varhistsPlotting[i]->SetLineColor(1);
                }
            }
            else {
                varhistsPlotting[i]->SetLineColor(1);
            }

            if (legends.at(i) != legends.at(i - 1)) {
                varhistsPlotting[i]->SetLineColor(1);
                stack->Add(varhistsPlotting[i]);
            }
            if (i > 0) {
                if (legends.at(i) != legends.at(i - 1)) {
                    legchange = i;
                    if ((legends.at(i) == DYEntry) && DYScale.at(channelType) != 1) {
                        leg->AddEntry(varhistsPlotting[i], legends.at(i), "f");
                    }
                    else leg->AddEntry(varhistsPlotting[i], legends.at(i), "f");
                }
                else { varhistsPlotting[legchange]->Add(varhistsPlotting[i]); }
            }
        }
        else { if (i == 0) leg->AddEntry(varhistsPlotting[i], legends.at(i), "pe"); }
    }

    ///////////////////////////////////
    //purity and stability plots as taken from CombinedCrossSection...

    TH1* genHist = (TH1*)GenPlotForEffAndAccCalculation->Clone();
    TH1* fidHist = NULL;
    if (particleXsec) fidHist = (TH1*)GenPlot->Clone();
    TH1* genRecHist = new TH1D("", "", bins, Xbins);
    int intbinsX[XAxisbins.size()];
    int intbinsY[XAxisbins.size()];

    // fill the elements of the main diagonal of the 2d hist into binned 1D histogram
    for (unsigned int i = 0; i < XAxisbins.size(); ++i) {
        intbinsX[i] = genReco2d->GetXaxis()->FindBin(Xbins[i] + 0.001);
        intbinsY[i] = genReco2d->GetYaxis()->FindBin(Xbins[i] + 0.001);

        if (i > 0) {
            genRecHist->SetBinContent(i, ((TH2D*)genReco2d)->Integral(intbinsX[i - 1], intbinsX[i] - 1, intbinsY[i - 1], intbinsY[i] - 1));
        }
    }

    TH1* genPseHist = ((TH2D*)genReco2d)->ProjectionY();
    TH1* recPseHist = ((TH2D*)genReco2d)->ProjectionX();

    TH1* genBinHist = genPseHist->Rebin(bins, "genBinHist", Xbins);
    TH1* recBinHist = recPseHist->Rebin(bins, "recBinHist", Xbins);

    genRecHist->SetBinContent(0, 0);
    genRecHist->SetBinContent(bins + 1, 0);
    genBinHist->SetBinContent(0, 0);
    genBinHist->SetBinContent(bins + 1, 0);
    recBinHist->SetBinContent(0, 0);
    recBinHist->SetBinContent(bins + 1, 0);
    genHist->SetBinContent(0, 0);
    genHist->SetBinContent(bins + 1, 0);
    if (particleXsec) fidHist->SetBinContent(0, 0);
    if (particleXsec) fidHist->SetBinContent(bins + 1, 0);

    // this is realy ugly but necessary:
    // As it seems, somewhere a double is tranformed into a float so that
    // efficiencies can be larger than 1.
    for (Int_t i = 1; i <= genRecHist->GetNbinsX(); ++i) {
        if (genRecHist->GetBinContent(i) > recBinHist->GetBinContent(i)) {
            genRecHist->SetBinContent(i, recBinHist->GetBinContent(i));
            std::cout << "WARNING in PlotDifferentialCrossSections: number of events generated and reconstructed in bin" << i
                << " = " << genRecHist->GetBinContent(i) << " is larger than number of reconstructed events in that bin"
                << " = " << recBinHist->GetBinContent(i) << std::endl;
        }
        if (genRecHist->GetBinContent(i) > genBinHist->GetBinContent(i)) {
            genRecHist->SetBinContent(i, genBinHist->GetBinContent(i));
            std::cout << "WARNING in PlotDifferentialCrossSections: number of events generated and reconstructed in bin " << i
                << " is larger than number of generated events in that bin" << std::endl;
        }
    }
    // efficiency, purity, stability
    TGraphAsymmErrors* grE; // for efficiency
    TGraphAsymmErrors* grP; // for purity
    TGraphAsymmErrors* grS; // for stability
    TGraphAsymmErrors* grA = NULL; // for acceptance   OZ 15.01.2017 initialized with NULL to eliminate warning 'may be used uninitialized'

    Bool_t effAboveOne = false;
    if (particleXsec) {
        for (Int_t i = 0; i <= fidHist->GetNbinsX() + 1; ++i) {
            if (fidHist->GetBinContent(i) < recBinHist->GetBinContent(i)) effAboveOne = true;
        }
    }

    // efficiency
    if (particleXsec && effAboveOne) grE = new TGraphAsymmErrors(recBinHist, fidHist, "pois");
    else if (particleXsec && !effAboveOne) grE = new TGraphAsymmErrors(recBinHist, fidHist);
    else grE = new TGraphAsymmErrors(recBinHist, genHist);
    grE->SetMinimum(0);
    grE->SetMaximum(1);
    grE->SetLineColor(8);
    grE->SetLineWidth(2);
    grE->SetMarkerSize(2);
    grE->SetMarkerStyle(21);
    grE->SetMarkerColor(8);

    // purity
    grP = new TGraphAsymmErrors(genRecHist, recBinHist);
    grP->SetLineColor(4);
    grP->SetLineWidth(2);
    grP->SetMarkerSize(2);
    grP->SetMarkerStyle(23);
    grP->SetMarkerColor(4);

    // stability
    grS = new TGraphAsymmErrors(genRecHist, genBinHist);
    grS->SetLineColor(2);
    grS->SetLineWidth(2);
    grS->SetMarkerSize(2);
    grS->SetMarkerStyle(22);
    grS->SetMarkerColor(2);

    // acceptance
    if (particleXsec) {
        grA = new TGraphAsymmErrors(fidHist, genHist);
        grA->SetLineColor(kViolet + 2);
        grA->SetLineWidth(2);
        grA->SetMarkerSize(2);
        grA->SetMarkerStyle(20);
        grA->SetMarkerColor(kViolet + 2);
    }


    grE->GetXaxis()->SetTitle(XAxis);
    TCanvas* cESP = new TCanvas("ESP", "ESP");

    // this is a dummy to get the x axis range corrct

    recBinHist->Reset();
    recBinHist->Draw();
    recBinHist->SetMaximum(1.);
    if (effAboveOne) recBinHist->SetMaximum(1.5);
    recBinHist->GetXaxis()->SetNoExponent(kTRUE);
    if (name.Contains("pT") || name.Contains("Mass")) {
        recBinHist->GetXaxis()->SetTitle(XAxis.Copy().Prepend("Reconstructed ").Append(" #left[GeV#right]"));
        if (name.Contains("Rapidity")) recBinHist->GetXaxis()->SetTitle(XAxis.Copy().Prepend("Reconstructed "));
    }
    else recBinHist->GetXaxis()->SetTitle(XAxis.Copy().Prepend("Reconstructed "));
    DrawCMSLabels(cmsPrelimLabelMode);
    DrawDecayChLabel(channelLabel[channelType]);
    grE->Draw("P,SAME");
    grP->Draw("P,SAME");
    grS->Draw("P,SAME");
    //if(particleXsec) grA->Draw("P,SAME"); 

    TLegend* leg3 = new TLegend();
    leg3->AddEntry(grE, "Efficiency", "p");
    leg3->AddEntry(grP, "Purity", "p");
    leg3->AddEntry(grS, "Stability", "p");
    //if(particleXsec) leg3->AddEntry(grA, "Particle acceptance", "p" );     
    setResultLegendStyle(leg3, 0);
    leg3->Draw("SAME");

    TString outdir;
    if (Channel == "combined" && performUnfoldingInCombinedChannel) outdir = ttbar::assignFolder(outpathPlots, "combinedUnfolded", TString("FinalResults"));
    else outdir = ttbar::assignFolder(outpathPlots, Channel, TString("FinalResults"));

    cESP->Print(outdir.Copy() + "ESP_" + name + ".eps");
    cESP->Print(outdir.Copy() + "ESP_" + name + ".pdf");
    cESP->Clear();
    delete cESP;

    init = false;
    for (unsigned int hist = 0; hist < hists.size(); hist++) {
        if (legends.at(hist) == "Data") {
            for (Int_t bin = 0; bin < bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
                DataSum[bin] += varhists[hist]->GetBinContent(bin + 1);
            }
        }
        else if ((legends.at(hist) == "t#bar{t} signal") && init == false) {
            signalHist = hist;
            init = true;
            for (Int_t bin = 0; bin < bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
                GenSignalSum[bin] += GenPlot->GetBinContent(bin + 1);
            }
        }
        else {
            for (Int_t bin = 0; bin < bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
                BGSum[bin] += varhists[hist]->GetBinContent(bin + 1);
            }
        }
    }
    double totalDataSum = 0;
    double GenDiffXSecPlot[XAxisbinCenters.size()];
    for (Int_t bin = 0; bin < bins; ++bin) {
        totalDataSum += DataSum[bin];
    }
    TH1* h_DiffXSec = (TH1D*)varhists[0]->Clone(); h_DiffXSec->Reset();
    TH1* h_GenDiffXSec = (TH1D*)varhists[0]->Clone(); h_GenDiffXSec->Reset();

    double xsecnorm = 1.0; // needed for abs. cross xsecs to normalize theory to data 
    //The systematic array is filled in the order in which the Stack is filled
    double DiffXSecPlot[XAxisbinCenters.size()];
    double DiffXSecStatErrorPlot[XAxisbinCenters.size()];
    double DiffXSecSysErrorbySysPlot[XAxisbinCenters.size()][(int)vec_systematic.size()];
    double DiffXSecSysErrorPlot[XAxisbinCenters.size()];
    double DiffXSecTotalErrorPlot[XAxisbinCenters.size()];

    double ModelSysPlot[XAxisbinCenters.size()];
    double ExpSysPlot[XAxisbinCenters.size()];

    if (evaluateTotCovMatricesForSVDUnfolding) CalcTotCovMtrx(Channel, name, vec_systematic);

    //Read absolute systematic uncertainty of each bin from file
    for (int Syst = 0; Syst < (int)vec_systematic.size(); Syst++) {
        if (vec_systematic.at(Syst).Contains("HERWIG")) continue;
        if (mergeScaleVariations_ && (vec_systematic.at(Syst).Contains("MEFACSCALE_") || vec_systematic.at(Syst).Contains("MERENSCALE_") || vec_systematic.at(Syst).Contains("PSISRSCALE_") ||
            vec_systematic.at(Syst).Contains("PSFSRSCALE_"))) continue;
        if (mergeBTagVariations_ && (vec_systematic.at(Syst) == "BTAG_PT_" || vec_systematic.at(Syst) == "BTAG_ETA_" || vec_systematic.at(Syst) == "BTAG_LJET_PT_" || vec_systematic.at(Syst) == "BTAG_LJET_ETA_")) continue;
        if (mergeTriggerVariations_ && (vec_systematic.at(Syst) == "TRIG_ETA_")) continue;
        if (mergeBFragVariations_ && (vec_systematic.at(Syst) == "BFRAG_PETERSON")) continue;
        if (mergeColorRecVariations_ && (vec_systematic.at(Syst) == "ERDONRETUNE" || vec_systematic.at(Syst) == "GLUONMOVETUNE")) continue;

        TString sysup = vec_systematic.at(Syst) + "UP";
        TString sysdown = vec_systematic.at(Syst) + "DOWN";
        if (vec_systematic.at(Syst) == "POWHEG")
        {
            sysup = "POWHEG";
            sysdown = "Nominal";
            vec_systematic.at(Syst) = "MOD_";
        };
        if (vec_systematic.at(Syst) == "MCATNLO")
        {
            sysup = "MCATNLO";
            sysdown = "Nominal";
            vec_systematic.at(Syst) = "HAD_";
        };
        if (vec_systematic.at(Syst) == "BFRAG_PETERSON")
        {
            sysup = "BFRAG_PETERSON";
            sysdown = "Nominal";
            vec_systematic.at(Syst) = "BFRAG_PETERSON_";
        };
        if (vec_systematic.at(Syst) == "ERDON")
        {
            sysup = "ERDON";
            sysdown = "Nominal";
            vec_systematic.at(Syst) = "ERDON_";
        };
        if (vec_systematic.at(Syst) == "ERDONRETUNE")
        {
            sysup = "ERDONRETUNE";
            sysdown = "Nominal";
            vec_systematic.at(Syst) = "ERDONRETUNE_";
        };
        if (vec_systematic.at(Syst) == "GLUONMOVETUNE")
        {
            sysup = "GLUONMOVETUNE";
            sysdown = "Nominal";
            vec_systematic.at(Syst) = "GLUONMOVETUNE_";
        };
        if (mergeScaleVariations_ && vec_systematic.at(Syst) == "MESCALE_")
        {
            sysup = "Nominal";
            sysdown = "Nominal";
            vec_systematic.at(Syst) = "TOT_SCALE_";
        };
        if (mergeBTagVariations_ && vec_systematic.at(Syst) == "BTAG_")
        {
            sysup = "Nominal";
            sysdown = "Nominal";
            vec_systematic.at(Syst) = "TOT_BTAG_";
        };
        if (mergeBTagVariations_ && vec_systematic.at(Syst) == "BTAG_LJET_")
        {
            sysup = "Nominal";
            sysdown = "Nominal";
            vec_systematic.at(Syst) = "TOT_BTAG_LJET_";
        };
        if (mergeTriggerVariations_ && vec_systematic.at(Syst) == "TRIG_")
        {
            sysup = "Nominal";
            sysdown = "Nominal";
            vec_systematic.at(Syst) = "TOT_TRIG_";
        };
        if (mergeBFragVariations_ && vec_systematic.at(Syst) == "BFRAG_")
        {
            sysup = "Nominal";
            sysdown = "Nominal";
            vec_systematic.at(Syst) = "TOT_BFRAG_";
        };
        if (mergeColorRecVariations_ && vec_systematic.at(Syst) == "ERDON_")
        {
            sysup = "Nominal";
            sysdown = "Nominal";
            vec_systematic.at(Syst) = "TOT_COLORREC_";
        };

        ifstream SysResultsList;
        if (Channel == "combined" && performUnfoldingInCombinedChannel) SysResultsList.open("UnfoldingResults/" + vec_systematic.at(Syst) + "/" + Channel + "Unfolded/" + name + "Results.txt");
        else SysResultsList.open("UnfoldingResults/" + vec_systematic.at(Syst) + "/" + Channel + "/" + name + "Results.txt");
        if (!SysResultsList.is_open()) {
            /// Apply symmetrization of errors only if the file doesn't exist!!
            /// If file exists means that the error is set by hand!!
            if ((mergeScaleVariations_ && vec_systematic.at(Syst) == "TOT_SCALE_") ||
                (mergeBTagVariations_ && vec_systematic.at(Syst) == "TOT_BTAG_") ||
                (mergeBTagVariations_ && vec_systematic.at(Syst) == "TOT_BTAG_LJET_") ||
                (mergeTriggerVariations_ && vec_systematic.at(Syst) == "TOT_TRIG_") ||
                (mergeBFragVariations_ && vec_systematic.at(Syst) == "TOT_BFRAG_") ||
                (mergeColorRecVariations_ && vec_systematic.at(Syst) == "TOT_COLORREC_")) Plotter13TeV::GetEnvelopeForUnfoldedResults(Channel, vec_systematic.at(Syst), name);
            else Plotter13TeV::CalcUpDownDifference(Channel, sysup, sysdown, name);
        }
        SysResultsList.close();
        if (Channel == "combined" && performUnfoldingInCombinedChannel) SysResultsList.open("UnfoldingResults/" + vec_systematic.at(Syst) + "/" + Channel + "Unfolded/" + name + "Results.txt");
        else SysResultsList.open("UnfoldingResults/" + vec_systematic.at(Syst) + "/" + Channel + "/" + name + "Results.txt");
        for (Int_t bin = 0; bin < bins; ++bin) {
            TString DUMMY;
            SysResultsList >> DUMMY >> XAxisbinCenters[bin] >> DUMMY >> Xbins[bin] >> DUMMY >> Xbins[bin + 1] >> DUMMY >> DiffXSecSysErrorbySysPlot[bin][Syst];
        }
        SysResultsList.close();
    }
    TString Dummy;
    //Read central and abolute statistical uncertainty values from Nominal
    ifstream ResultsList;
    if (Channel == "combined" && performUnfoldingInCombinedChannel) ResultsList.open("UnfoldingResults/Nominal/" + Channel + "Unfolded/" + name + "Results.txt");
    else ResultsList.open("UnfoldingResults/Nominal/" + Channel + "/" + name + "Results.txt");
    for (Int_t bin = 0; bin < bins; bin++) {//Retrieve arrays for plotting
        ResultsList >> Dummy >> XAxisbinCenters[bin] >> Dummy >> Xbins[bin] >> Dummy >> Xbins[bin + 1] >> Dummy >> DiffXSecPlot[bin] >> Dummy >> DiffXSecStatErrorPlot[bin] >> Dummy >> GenDiffXSecPlot[bin];
        h_DiffXSec->SetBinContent(bin + 1, DiffXSecPlot[bin]);
        h_DiffXSec->SetBinError(bin + 1, DiffXSecStatErrorPlot[bin]);
        h_GenDiffXSec->SetBinContent(bin + 1, GenDiffXSecPlot[bin]);
    }

    double TotalVisXSection = 1.; // this can currently be set to 1. Re-normalization, if needed, is maintained through SVD functions. Can be changed in case of need.

    if (!absoluteXsec) TotalVisXSection = h_DiffXSec->Integral("width");
    else {
        xsecnorm = h_DiffXSec->Integral("width");
        if ((name.Contains("Lepton") || name.Contains("Top") || name.Contains("BJet")) && !name.Contains("Lead") && !revokeAntiQuantityCombination) {
            TotalVisXSection = 2.0;
            xsecnorm = 0.5 * xsecnorm;
        }

    }
    h_DiffXSec->Scale(1 / TotalVisXSection);

    for (Int_t bin = 0; bin < bins; bin++) {//condense matrices to arrays for plotting
        double syst_square = 0;
        if (absoluteXsec) syst_square = lumiError * lumiError;
        if (absoluteXsec && !particleXsec) syst_square += BRError * BRError;
        DiffXSecSysErrorPlot[bin] = 0.;
        ExpSysPlot[bin] = 0.;
        ModelSysPlot[bin] = 0.;
        for (int syst = 0; syst < (int)vec_systematic.size(); syst++) {
            if (mergeScaleVariations_ && (vec_systematic.at(syst).Contains("MEFACSCALE_") || vec_systematic.at(syst).Contains("MERENSCALE_") || vec_systematic.at(syst).Contains("PSISRSCALE_") ||
                vec_systematic.at(syst).Contains("PSFSRSCALE_"))) continue;
            if (mergeBTagVariations_ && (vec_systematic.at(syst) == "BTAG_PT_" || vec_systematic.at(syst) == "BTAG_ETA_" || vec_systematic.at(syst) == "BTAG_LJET_PT_" || vec_systematic.at(syst) == "BTAG_LJET_ETA_")) continue;
            if (mergeTriggerVariations_ && (vec_systematic.at(syst) == "TRIG_ETA_")) continue;
            if (mergeBFragVariations_ && (vec_systematic.at(syst) == "BFRAG_PETERSON")) continue;
            if (mergeColorRecVariations_ && (vec_systematic.at(syst) == "ERDONRETUNE" || vec_systematic.at(syst) == "GLUONMOVETUNE")) continue;


            syst_square += DiffXSecSysErrorbySysPlot[bin][syst] * DiffXSecSysErrorbySysPlot[bin][syst];
            if (vec_systematic.at(syst).Contains("JER_") || vec_systematic.at(syst).Contains("JES_") || vec_systematic.at(syst).Contains("PU_") ||
                vec_systematic.at(syst).Contains("DY_") || vec_systematic.at(syst).Contains("BG_") || vec_systematic.at(syst).Contains("TRIG_") ||
                vec_systematic.at(syst).Contains("LEPT_") || vec_systematic.at(syst).Contains("BTAG_") || vec_systematic.at(syst).Contains("KIN_")) {
                ExpSysPlot[bin] += DiffXSecSysErrorbySysPlot[bin][syst] * DiffXSecSysErrorbySysPlot[bin][syst];
            }
            else {
                ModelSysPlot[bin] += DiffXSecSysErrorbySysPlot[bin][syst] * DiffXSecSysErrorbySysPlot[bin][syst];
            }
        }
        if (absoluteXsec) ExpSysPlot[bin] += lumiError * lumiError;
        if (absoluteXsec && !particleXsec) ExpSysPlot[bin] += BRError * BRError;
        DiffXSecSysErrorPlot[bin] += syst_square;
        ExpSysPlot[bin] = sqrt(ExpSysPlot[bin]);
        ModelSysPlot[bin] = sqrt(ModelSysPlot[bin]);
        DiffXSecStatErrorPlot[bin] = DiffXSecStatErrorPlot[bin] / TotalVisXSection;
        DiffXSecPlot[bin] = DiffXSecPlot[bin] / TotalVisXSection;
        DiffXSecSysErrorPlot[bin] = sqrt(DiffXSecSysErrorPlot[bin]) * DiffXSecPlot[bin]; //absolute systematic error in bin 'bin'
        DiffXSecTotalErrorPlot[bin] = sqrt(DiffXSecSysErrorPlot[bin] * DiffXSecSysErrorPlot[bin] + DiffXSecStatErrorPlot[bin] * DiffXSecStatErrorPlot[bin]);//absolute total error
    }

    //The Markus plots
    TCanvas* c10 = new TCanvas("Markus", "Markus");
    THStack* SystHists = new THStack("MSTACK", "MSTACK");
    TLegend* leg10 = new TLegend(0.20, 0.65, 0.45, 0.90);

    FILE* systfile;
    systfile = fopen(outdir.Copy() + newname + "_SystematicsLaTeX.txt", "w");


    for (int systs = 0; systs < (int)vec_systematic.size(); systs++) {
        if (vec_systematic.at(systs) == "BTAG_ETA_" || vec_systematic.at(systs) == "BTAG_LJET_ETA_") { continue; }//Skip the BTAG_ETA systematic because it's added in quadrature to BTAG_PT
        if (mergeScaleVariations_ && (vec_systematic.at(systs).Contains("MEFACSCALE_") || vec_systematic.at(systs).Contains("MERENSCALE_") || vec_systematic.at(systs).Contains("PSISRSCALE_") ||
            vec_systematic.at(systs).Contains("PSFSRSCALE_"))) continue;
        if (mergeBTagVariations_ && (vec_systematic.at(systs) == "BTAG_PT_" || vec_systematic.at(systs) == "BTAG_ETA_" || vec_systematic.at(systs) == "BTAG_LJET_PT_" || vec_systematic.at(systs) == "BTAG_LJET_ETA_")) continue;
        if (mergeTriggerVariations_ && (vec_systematic.at(systs) == "TRIG_ETA_")) continue;
        if (mergeBFragVariations_ && (vec_systematic.at(systs) == "BFRAG_PETERSON")) continue;
        if (mergeColorRecVariations_ && (vec_systematic.at(systs) == "ERDONRETUNE" || vec_systematic.at(systs) == "GLUONMOVETUNE")) continue;

        TH1D* systtemp = (TH1D*)varhists[0]->Clone();
        systtemp->Reset();
        double TotalSyst = 0.0, TotalSqSyst = 0.0;
        //double AvgSyst, SqAvgSys; // OZ 15.01.2017 commented out to eliminate 'set but not used' warning

        TString systNameTemp = vec_systematic.at(systs);
        systNameTemp.Remove(systNameTemp.Length() - 1, 1);

        // Fill header of txt file, fill only for the first iteration systematic
        if (systs < 1) {
            fprintf(systfile, "Source ");
            for (Int_t bin = 0; bin < bins; bin++) {//condense matrices to arrays for plotting
                if (name.Contains("Rapidity") || name.Contains("Eta") || name.Contains("Phi") || name.Contains("DeltaR")) {
                    fprintf(systfile, "&  [%1.2f, %1.2f]", Xbins[bin], Xbins[bin + 1]);
                }
                else if (name.Contains("JetMultp")) {
                    fprintf(systfile, "&  %1.0f", Xbins[bin] + 0.5);
                }
                else { fprintf(systfile, "&  [%.0f, %.0f]", Xbins[bin], Xbins[bin + 1]); }
            }
            fprintf(systfile, "\\\\ \\hline \n");
        }

        for (Int_t bin = 0; bin < bins; bin++) {//condense matrices to arrays for plotting
            if (vec_systematic.at(systs) == "BTAG_PT_" || vec_systematic.at(systs) == "BTAG_LJET_PT_") {
                DiffXSecSysErrorbySysPlot[bin][systs] = TMath::Sqrt(
                    (DiffXSecSysErrorbySysPlot[bin][systs] * DiffXSecSysErrorbySysPlot[bin][systs])
                    + (DiffXSecSysErrorbySysPlot[bin][systs + 1] * DiffXSecSysErrorbySysPlot[bin][systs + 1])
                );
            }

            systtemp->SetBinContent(bin + 1, (DiffXSecSysErrorbySysPlot[bin][systs] * DiffXSecSysErrorbySysPlot[bin][systs]));
            if (bin == 0) {
                fprintf(systfile, "%s   ", systNameTemp.Data());
            }
            fprintf(systfile, " &  %2.2f  ", TMath::Sqrt(systtemp->GetBinContent(bin + 1)) * 100);
            if (bin > 0 && bin < bins - 1) {//Exclude the 2 side bins
                TotalSyst = TotalSyst + TMath::Sqrt(systtemp->GetBinContent(bin + 1));
                TotalSqSyst = TotalSqSyst + systtemp->GetBinContent(bin + 1);
            }
        }
        //AvgSyst=TotalSyst/(bins-2); // OZ 15.01.2017 commented out to eliminate 'set but not used' warning
        //SqAvgSys=TMath::Sqrt(TotalSqSyst/(bins-2)); // OZ 15.01.2017 commented out to eliminate 'set but not used' warning
        fprintf(systfile, "\\\\ \n");
        //         fprintf(systfile, "\\\\ Lin.Avg.(%%)= %.5f  Quad.Avg.(%%)=%.5f\n", 100*AvgSyst, 100*SqAvgSys);
        int legColorId = (int)vec_systematic.size() - systs + 1;
        if (legColorId > 10) legColorId += 27;
        systtemp->SetFillColor(legColorId);
        SystHists->Add(systtemp);

        if (vec_systematic.at(systs) == "BTAG_PT_") {
            leg10->AddEntry(systtemp, "BTAG_PT/ETA", "f");
        }
        else if (vec_systematic.at(systs) == "BTAG_LJET_PT_") {
            leg10->AddEntry(systtemp, "BTAG_LJET_PT/ETA", "f");
        }
        else leg10->AddEntry(systtemp, systNameTemp, "f");
    }

    SystHists->Draw();
    fclose(systfile);

    if (vec_systematic.size() > 0) {
        if (name.Contains("pT") || name.Contains("Mass")) {
            SystHists->GetHistogram()->GetXaxis()->SetTitle(XAxis.Copy().Append(" #left[GeV#right]"));
            if (name.Contains("Rapidity")) SystHists->GetHistogram()->GetXaxis()->SetTitle(XAxis);
        }
        else  SystHists->GetHistogram()->GetXaxis()->SetTitle(XAxis);
        SystHists->GetHistogram()->GetYaxis()->SetTitle("#sum #left( #frac{#Delta #sigma}{#sigma} #right)^{2}");
        SystHists->GetXaxis()->SetNoExponent(kTRUE);

        leg10->SetFillColor(0);
        leg10->Draw("SAME");
        c10->Print(outdir.Copy() + "MSP_" + name + ".eps");
        c10->Print(outdir.Copy() + "MSP_" + name + ".pdf");
        delete leg10;
    }
    delete c10;

    //The Experimental/Model/Statistical plot
    TCanvas* c11 = new TCanvas("EMS", "EMS");
    TH1D* ExpHist = (TH1D*)varhists[0]->Clone();    ExpHist->Reset();
    TH1D* ModelHist = (TH1D*)varhists[0]->Clone();  ModelHist->Reset();
    TH1D* StatHist = (TH1D*)varhists[0]->Clone();   StatHist->Reset();
    TH1D* TotalHist = (TH1D*)varhists[0]->Clone();  TotalHist->Reset();
    for (Int_t bin = 0; bin < bins; bin++) {//condense matrices to arrays for plotting
        ExpHist->SetBinContent(bin + 1, 100 * ExpSysPlot[bin]);
        ModelHist->SetBinContent(bin + 1, 100 * ModelSysPlot[bin]);
        StatHist->SetBinContent(bin + 1, 100 * DiffXSecStatErrorPlot[bin] / DiffXSecPlot[bin]);
        TotalHist->SetBinContent(bin + 1, 100 * DiffXSecTotalErrorPlot[bin] / DiffXSecPlot[bin]);
    }
    TotalHist->SetMinimum(0.);
    TotalHist->GetYaxis()->SetTitle("#frac{#Delta#sigma}{#sigma} [%]");
    TotalHist->SetLineColor(1);
    ExpHist->SetLineColor(kRed);
    StatHist->SetLineColor(kGreen);
    ModelHist->SetLineColor(kBlue);
    TLegend* leg11 = new TLegend(0.65, 0.60, 0.90, 0.85);
    leg11->SetFillColor(0);
    leg11->AddEntry(ExpHist->Clone(), "Experimental Uncertainty", "l");
    leg11->AddEntry(StatHist->Clone(), "Statistical Uncertainty", "l");
    leg11->AddEntry(ModelHist->Clone(), "Model Uncertainty", "l");
    leg11->AddEntry(TotalHist->Clone(), "Total Uncertainty", "l");
    TotalHist->Draw();ModelHist->Draw("SAME");ExpHist->Draw("SAME");StatHist->Draw("SAME");
    leg11->Draw("SAME");
    TotalHist->GetXaxis()->SetNoExponent(kTRUE);
    c11->Print(outdir.Copy() + "SEM_" + name + ".eps");
    c11->Print(outdir.Copy() + "SEM_" + name + ".pdf");
    c11->Clear();
    delete c11;

    Double_t mexl[XAxisbinCenters.size()];
    Double_t mexh[XAxisbinCenters.size()];
    for (unsigned int j = 0; j < XAxisbinCenters.size();j++) { mexl[j] = 0;mexh[j] = 0; }
    TGraphAsymmErrors* tga_DiffXSecPlot = new TGraphAsymmErrors(bins, binCenters, DiffXSecPlot, mexl, mexh, DiffXSecStatErrorPlot, DiffXSecStatErrorPlot);
    tga_DiffXSecPlot->SetMarkerStyle(1);
    tga_DiffXSecPlot->SetMarkerColor(kBlack);
    tga_DiffXSecPlot->SetMarkerSize(0.8);
    tga_DiffXSecPlot->SetLineWidth(2);
    tga_DiffXSecPlot->SetLineColor(kBlack);

    TGraphAsymmErrors* tga_DiffXSecPlotwithSys = new TGraphAsymmErrors(bins, binCenters, DiffXSecPlot, mexl, mexh, DiffXSecTotalErrorPlot, DiffXSecTotalErrorPlot);
    tga_DiffXSecPlotwithSys->SetMarkerStyle(20);
    tga_DiffXSecPlotwithSys->SetMarkerColor(kBlack);
    tga_DiffXSecPlotwithSys->SetMarkerSize(0.8);
    tga_DiffXSecPlotwithSys->SetLineWidth(2);
    tga_DiffXSecPlotwithSys->SetLineColor(kBlack);

    //GenPlotTheory->Scale(topxsec/(SignalEventswithWeight*GenPlotTheory->GetBinWidth(1)));
    GenPlot->Scale(topxsec / (SignalEventswithWeight * GenPlot->GetBinWidth(1)));
    if ((name.Contains("Lepton") || name.Contains("Top") || name.Contains("BJet")) && !name.Contains("Lead") && !revokeAntiQuantityCombination) {
        GenPlotTheory->Scale(1. / 2.);
    }
    //     GenPlotTheory->Rebin(2);
    GenPlotTheory->Scale(GenPlotTheory->Integral("width"));
    h_GenDiffXSec->Scale(h_GenDiffXSec->Integral("width"));

    TH1* madgraphhist = 0, * mcnlohist = 0, * mcnlohistup = 0, * mcnlohistdown = 0, * powheghist = 0, * powhegHerwighist = 0, * perugia11hist = 0;
    TH1* mcnlohistnorm = 0;
    TGraph* mcatnloBand = 0;

    TH1* mcnlohistnormBinned = 0, * mcnlohistupBinned = 0;
    TH1* mcnlohistdownBinned = 0, * mcnlohistBinned = 0;
    TH1* madgraphhistBinned = 0, * powheghistBinned = 0, * powhegHerwighistBinned = 0;
    TH1* perugia11histBinned = 0;

    TH1* Kidoth1_Binned = 0, * Ahrensth1_Binned = 0;

    TH1* madup = 0, * maddown = 0, * matchup = 0, * matchdown = 0, * match2up = 0, * match2down = 0;
    TH1* madupBinned = 0, * maddownBinned = 0, * matchupBinned = 0, * matchdownBinned = 0, * match2upBinned = 0, * match2downBinned = 0;

    madgraphhist = GetNloCurve(newname, "Nominal");
    madgraphhistBinned = RescaleNLOCurveAndPrepareBinnedCurve(madgraphhist, bins, Xbins, xsecnorm, BranchingFraction[channelType]);

    bool canDrawMCATNLO = true;
    if (drawNLOCurves && drawMCATNLO) {
        mcnlohist = GetNloCurve(newname, "MCATNLO");
        mcnlohistBinned = RescaleNLOCurveAndPrepareBinnedCurve(mcnlohist, bins, Xbins, xsecnorm, BranchingFraction[channelType]);

        if (name.Contains("LeptonpT")) { mcnlohistnorm = GetNloCurve("Leptons", "Pt", "MCatNLO"); }//temprorary until I change the naming convention in the root file
        else if (name.Contains("LeptonEta")) { mcnlohistnorm = GetNloCurve("Leptons", "Eta", "MCatNLO"); }
        else if (name.Contains("LLBarpT")) { mcnlohistnorm = GetNloCurve("LepPair", "Pt", "MCatNLO"); }
        else if (name.Contains("LLBarMass")) { mcnlohistnorm = GetNloCurve("LepPair", "Mass", "MCatNLO"); }
        else if (name.Contains("ToppT")) { mcnlohistnorm = GetNloCurve("TopQuarks", "Pt", "MCatNLO"); }
        else if (name.Contains("TopRapidity")) { mcnlohistnorm = GetNloCurve("TopQuarks", "Rapidity", "MCatNLO"); }
        else if (name.Contains("TTBarpT")) { mcnlohistnorm = GetNloCurve("TtBar", "Pt", "MCatNLO"); }
        else if (name.Contains("TTBarRapidity")) { mcnlohistnorm = GetNloCurve("TtBar", "Rapidity", "MCatNLO"); }
        else if (name.Contains("TTBarMass")) { mcnlohistnorm = GetNloCurve("TtBar", "Mass", "MCatNLO"); }
        else if (name.Contains("BJetpT")) { mcnlohistnorm = GetNloCurve("Jets", "Pt", "MCatNLO"); }
        else if (name.Contains("BJetEta")) { mcnlohistnorm = GetNloCurve("Jets", "Eta", "MCatNLO"); }
        //else if(name.Contains("LeptonBJetMass")){mcnlohistnorm = GetNloCurve("Jets","Eta","MCatNLO");}

        if (mcnlohistnorm) {
            mcnlohistnormBinned = mcnlohistnorm->Rebin(bins, "genBinHistNorm", Xbins);

            if (name.Contains("LeptonpT")) { mcnlohistup = GetNloCurve("Leptons", "Pt", "MCNLOup"); }//temprorary until I change the naming convention in the root file
            else if (name.Contains("LeptonEta")) { mcnlohistup = GetNloCurve("Leptons", "Eta", "MCNLOup"); }
            else if (name.Contains("LLBarpT")) { mcnlohistup = GetNloCurve("LepPair", "Pt", "MCNLOup"); }
            else if (name.Contains("LLBarMass")) { mcnlohistup = GetNloCurve("LepPair", "Mass", "MCNLOup"); }
            else if (name.Contains("ToppT")) { mcnlohistup = GetNloCurve("TopQuarks", "Pt", "MCNLOup"); }
            else if (name.Contains("TopRapidity")) { mcnlohistup = GetNloCurve("TopQuarks", "Rapidity", "MCNLOup"); }
            else if (name.Contains("TTBarpT")) { mcnlohistup = GetNloCurve("TtBar", "Pt", "MCNLOup"); }
            else if (name.Contains("TTBarRapidity")) { mcnlohistup = GetNloCurve("TtBar", "Rapidity", "MCNLOup"); }
            else if (name.Contains("TTBarMass")) { mcnlohistup = GetNloCurve("TtBar", "Mass", "MCNLOup"); }
            else if (name.Contains("BJetpT")) { mcnlohistup = GetNloCurve("Jets", "Pt", "MCNLOup"); }
            else if (name.Contains("BJetEta")) { mcnlohistup = GetNloCurve("Jets", "Eta", "MCNLOup"); }
            mcnlohistupBinned = mcnlohistup->Rebin(bins, "genBinHist", Xbins);


            if (name.Contains("LeptonpT")) { mcnlohistdown = GetNloCurve("Leptons", "Pt", "MCNLOdown"); }//temprorary until I change the naming convention in the root file
            else if (name.Contains("LeptonEta")) { mcnlohistdown = GetNloCurve("Leptons", "Eta", "MCNLOdown"); }
            else if (name.Contains("LLBarpT")) { mcnlohistdown = GetNloCurve("LepPair", "Pt", "MCNLOdown"); }
            else if (name.Contains("LLBarMass")) { mcnlohistdown = GetNloCurve("LepPair", "Mass", "MCNLOdown"); }
            else if (name.Contains("ToppT")) { mcnlohistdown = GetNloCurve("TopQuarks", "Pt", "MCNLOdown"); }
            else if (name.Contains("TopRapidity")) { mcnlohistdown = GetNloCurve("TopQuarks", "Rapidity", "MCNLOdown"); }
            else if (name.Contains("TTBarpT")) { mcnlohistdown = GetNloCurve("TtBar", "Pt", "MCNLOdown"); }
            else if (name.Contains("TTBarRapidity")) { mcnlohistdown = GetNloCurve("TtBar", "Rapidity", "MCNLOdown"); }
            else if (name.Contains("TTBarMass")) { mcnlohistdown = GetNloCurve("TtBar", "Mass", "MCNLOdown"); }
            else if (name.Contains("BJetpT")) { mcnlohistdown = GetNloCurve("Jets", "Pt", "MCNLOdown"); }
            else if (name.Contains("BJetEta")) { mcnlohistdown = GetNloCurve("Jets", "Eta", "MCNLOdown"); }
            mcnlohistdownBinned = mcnlohistdown->Rebin(bins, "genBinHist", Xbins);

            for (Int_t bin = 0; bin < bins; bin++) {
                mcnlohistupBinned->SetBinContent(bin + 1, mcnlohistupBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / mcnlohistup->GetBinWidth(1)));
                mcnlohistdownBinned->SetBinContent(bin + 1, mcnlohistdownBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / mcnlohistdown->GetBinWidth(1)));
                mcnlohistnormBinned->SetBinContent(bin + 1, mcnlohistnormBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / mcnlohistnorm->GetBinWidth(1)));
            }
            mcnlohistupBinned->Scale(xsecnorm / mcnlohistnormBinned->Integral("width"));
            mcnlohistdownBinned->Scale(xsecnorm / mcnlohistnormBinned->Integral("width"));
            mcnlohistnormBinned->Scale(xsecnorm / mcnlohistnormBinned->Integral("width"));

            for (Int_t bin = 0; bin < bins; bin++) {
                mcnlohistupBinned->SetBinContent(bin + 1, (mcnlohistupBinned->GetBinContent(bin + 1) / mcnlohistnormBinned->GetBinContent(bin + 1)) * mcnlohistBinned->GetBinContent(bin + 1));
                mcnlohistdownBinned->SetBinContent(bin + 1, (mcnlohistdownBinned->GetBinContent(bin + 1) / mcnlohistnormBinned->GetBinContent(bin + 1)) * mcnlohistBinned->GetBinContent(bin + 1));
            }

            //Uncertainty band for MC@NLO
            Double_t x[bins];
            Double_t xband[2 * bins];
            Double_t errup[bins];
            Double_t errdn[bins];
            Double_t errorband[2 * bins];

            for (Int_t j = 0; j < bins; j++) {
                x[j] = mcnlohistBinned->GetBinCenter(j + 1);
                errup[j] = (mcnlohistupBinned->GetBinContent(j + 1) / mcnlohistnormBinned->GetBinContent(j + 1)) * mcnlohistBinned->GetBinContent(j + 1);
                errdn[j] = (mcnlohistdownBinned->GetBinContent(j + 1) / mcnlohistnormBinned->GetBinContent(j + 1)) * mcnlohistBinned->GetBinContent(j + 1);

                xband[j] = x[j];
                errorband[j] = errdn[j]; //lower band
                xband[2 * bins - j - 1] = x[j];
                errorband[2 * bins - j - 1] = errup[j]; //upper band
            }

            mcatnloBand = new TGraph(2 * bins, xband, errorband);
            mcatnloBand->SetFillColor(kGray);
            mcatnloBand->SetFillStyle(1001);
            mcatnloBand->SetLineColor(kBlue);
            mcatnloBand->SetLineWidth(2);
            mcatnloBand->SetLineStyle(5);
            canDrawMCATNLO = false;
        }
        else {
            std::cout << "\n*************************\nMC@NLO Curve not available!\n**********************\n";
            canDrawMCATNLO = false;
        }
    }
    if (drawNLOCurves && drawPOWHEG) {
        powheghist = GetNloCurve(newname, "POWHEG");
        powheghistBinned = RescaleNLOCurveAndPrepareBinnedCurve(powheghist, bins, Xbins, xsecnorm, BranchingFraction[channelType]);
    }
    if (drawNLOCurves && drawPOWHEGHERWIG) {
        powhegHerwighist = GetNloCurve(newname, "POWHEGHERWIG");
        powhegHerwighistBinned = RescaleNLOCurveAndPrepareBinnedCurve(powhegHerwighist, bins, Xbins, xsecnorm, BranchingFraction[channelType]);
    }

    if (drawNLOCurves && drawPERUGIA11) {
        perugia11hist = GetNloCurve(newname, "PERUGIA11");
        perugia11histBinned = RescaleNLOCurveAndPrepareBinnedCurve(perugia11hist, bins, Xbins, xsecnorm, BranchingFraction[channelType]);
    }
    if (drawNLOCurves && drawKidonakis &&
        (name == "HypToppT" || name == "HypTopRapidity") && !name.Contains("Lead") && !name.Contains("RestFrame") && !revokeAntiQuantityCombination) {
        TString kidoFile = ttbar::DATA_PATH_DILEPTONIC() + "/kidonakisNNLO_8TeV.root";
        if (name.Contains("ToppT"))             Kidoth1_Binned = fileReader->GetClone<TH1>(kidoFile, "topPt");
        else if (name.Contains("TopRapidity"))  Kidoth1_Binned = fileReader->GetClone<TH1>(kidoFile, "topY");
        for (int iter = 0; iter <= Kidoth1_Binned->GetNbinsX(); iter++) Kidoth1_Binned->SetBinError(iter, 1e-10);
        Kidoth1_Binned->Scale(1. / Kidoth1_Binned->Integral("width"));
    }

    if (drawNLOCurves && drawAhrens && (name == "HypTTBarMass" || name == "HypTTBarpT")) {
        TString ahrensFile = ttbar::DATA_PATH_DILEPTONIC() + "/ahrensNNLL_8TeV.root";
        if (name == "HypTTBarMass")    Ahrensth1_Binned = fileReader->GetClone<TH1>(ahrensFile, "ttbarM");
        else if (name == "HypTTBarpT") Ahrensth1_Binned = fileReader->GetClone<TH1>(ahrensFile, "ttbarPt");
        for (int iter = 0; iter <= Ahrensth1_Binned->GetNbinsX(); iter++) Ahrensth1_Binned->SetBinError(iter, 1e-10);
        Ahrensth1_Binned->Scale(1. / Ahrensth1_Binned->Integral("width"));
    }
    else { drawAhrens = 0; }

    if (drawMadScaleMatching) {
        //    if(drawNLOCurves && drawMadScaleMatching){
        madup = GetNloCurve(newname, "SCALE_UP");
        madupBinned = RescaleNLOCurveAndPrepareBinnedCurve(madup, bins, Xbins, xsecnorm, BranchingFraction[channelType]);

        maddown = GetNloCurve(newname, "SCALE_DOWN");
        maddownBinned = RescaleNLOCurveAndPrepareBinnedCurve(maddown, bins, Xbins, xsecnorm, BranchingFraction[channelType]);

        matchup = GetNloCurve(newname, "MATCH_UP");
        matchupBinned = RescaleNLOCurveAndPrepareBinnedCurve(matchup, bins, Xbins, xsecnorm, BranchingFraction[channelType]);

        matchdown = GetNloCurve(newname, "MATCH_DOWN");
        matchdownBinned = RescaleNLOCurveAndPrepareBinnedCurve(matchdown, bins, Xbins, xsecnorm, BranchingFraction[channelType]);

    }

    if (drawNLOCurves && drawMadMass) {
        madup = GetNloCurveMass(newname, "MASS_UP", "173.5");
        madup->Scale(xsecnorm / madup->Integral("width"));
        madupBinned = madup->Rebin(bins, "madupplot", Xbins);
        for (Int_t bin = 0; bin < bins; bin++) {//condense matrices to arrays for plotting
            madupBinned->SetBinContent(bin + 1, madupBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / madup->GetBinWidth(1)));
        }
        madupBinned->Scale(xsecnorm / madupBinned->Integral("width"));

        maddown = GetNloCurveMass(newname, "MASS_DOWN", "171.5");
        maddown->Scale(xsecnorm / maddown->Integral("width"));
        maddownBinned = maddown->Rebin(bins, "maddownplot", Xbins);
        for (Int_t bin = 0; bin < bins; bin++) {//condense matrices to arrays for plotting
            maddownBinned->SetBinContent(bin + 1, maddownBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / maddown->GetBinWidth(1)));
        }
        maddownBinned->Scale(xsecnorm / maddownBinned->Integral("width"));

        matchup = GetNloCurveMass(newname, "MASS_UP", "175.5");
        matchup->Scale(xsecnorm / matchup->Integral("width"));
        matchupBinned = matchup->Rebin(bins, "matchupplot", Xbins);
        for (Int_t bin = 0; bin < bins; bin++) {//condense matrices to arrays for plotting
            matchupBinned->SetBinContent(bin + 1, matchupBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / matchup->GetBinWidth(1)));
        }
        matchupBinned->Scale(xsecnorm / matchupBinned->Integral("width"));

        matchdown = GetNloCurveMass(newname, "MASS_DOWN", "169.5");
        matchdown->Scale(xsecnorm / matchdown->Integral("width"));
        matchdownBinned = matchdown->Rebin(bins, "matchdownplot", Xbins);
        for (Int_t bin = 0; bin < bins; bin++) {//condense matrices to arrays for plotting
            matchdownBinned->SetBinContent(bin + 1, matchdownBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / matchdown->GetBinWidth(1)));
        }
        matchdownBinned->Scale(xsecnorm / matchdownBinned->Integral("width"));

        match2up = GetNloCurveMass(newname, "MASS_UP", "178.5");
        double match2scale = 1. / match2up->Integral("width");
        match2up->Scale(match2scale);
        match2upBinned = match2up->Rebin(bins, "match2upplot", Xbins);
        for (Int_t bin = 0; bin < bins; bin++) {//condense matrices to arrays for plotting
            match2upBinned->SetBinContent(bin + 1, match2upBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / match2up->GetBinWidth(1)));
        }
        match2upBinned->Scale(xsecnorm / match2upBinned->Integral("width"));

        match2down = GetNloCurveMass(newname, "MASS_DOWN", "166.5");
        match2down->Scale(xsecnorm / match2down->Integral("width"));
        match2downBinned = match2down->Rebin(bins, "match2downplot", Xbins);
        for (Int_t bin = 0; bin < bins; bin++) {//condense matrices to arrays for plotting
            match2downBinned->SetBinContent(bin + 1, match2downBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / match2down->GetBinWidth(1)));
        }
        match2downBinned->Scale(xsecnorm / match2downBinned->Integral("width"));

    }

    TCanvas* c = new TCanvas("DiffXS", "DiffXS");
    if (logY) {
        c->SetLogy();
    }
    h_DiffXSec->SetMarkerStyle(tga_DiffXSecPlotwithSys->GetMarkerStyle());
    h_DiffXSec->SetMarkerSize(tga_DiffXSecPlotwithSys->GetMarkerSize());
    h_DiffXSec->SetMarkerColor(tga_DiffXSecPlotwithSys->GetMarkerColor());
    //MCFMHist->SetMarkerStyle(2);
    if (ymax != 0) {
        if (logY) {
            madgraphhistBinned->SetMaximum(18 * madgraphhistBinned->GetBinContent(madgraphhistBinned->GetMaximumBin()));
        }
        else { madgraphhistBinned->SetMaximum(1.5 * madgraphhistBinned->GetBinContent(madgraphhistBinned->GetMaximumBin())); }
    }
    madgraphhistBinned->GetXaxis()->SetNoExponent(kTRUE);
    if (name.Contains("Rapidity") || name.Contains("Eta")) { madgraphhistBinned->GetYaxis()->SetNoExponent(kTRUE); }
    if (name.Contains("HypJetMultpt")) {
        //madgraphhistBinned->GetXaxis()->SetNdivisions(madgraphhistBinned->GetNbinsX(),0,0, 1);
        TString TitBin = "";
        for (int bin = 1; bin <= madgraphhistBinned->GetNbinsX(); bin++) {
            if (bin == madgraphhistBinned->GetNbinsX()) { TitBin += "#geq"; TitBin += madgraphhistBinned->GetBinCenter(bin); madgraphhistBinned->GetXaxis()->SetBinLabel(bin, TitBin); }
            else {
                TitBin += madgraphhistBinned->GetBinCenter(bin);
                madgraphhistBinned->GetXaxis()->SetBinLabel(bin, TitBin);
            }
            TitBin = "";
        }
    }


    if (ymax != 0) madgraphhistBinned->SetMaximum(ymax);
    if (ymin != 0) madgraphhistBinned->SetMinimum(ymin);
    gStyle->SetEndErrorSize(8);
    if (drawNLOCurves && drawMCATNLO && canDrawMCATNLO) {
        //    mcatnloBand->Draw("same, F");
        mcnlohistupBinned->SetFillColor(kGray);
        mcnlohistupBinned->SetLineColor(kGray);
        mcnlohistupBinned->Draw("same");
        mcnlohistdownBinned->SetLineColor(10);
        mcnlohistdownBinned->SetFillColor(10);
        mcnlohistdownBinned->Draw("same");
    }
    GenPlotTheory->SetLineColor(kRed + 1);
    GenPlotTheory->SetLineWidth(2);
    GenPlotTheory->SetLineStyle(1);
    h_GenDiffXSec->SetLineColor(GenPlotTheory->GetLineColor());
    h_GenDiffXSec->SetLineStyle(GenPlotTheory->GetLineStyle());

    // Plot statistical and totl error of full result
    TGraphAsymmErrors* ratio_stat = 0, * ratio_tota = 0;
    if (tga_DiffXSecPlot) ratio_stat = (TGraphAsymmErrors*)tga_DiffXSecPlot->Clone("ratio_stat");
    if (tga_DiffXSecPlotwithSys) ratio_tota = (TGraphAsymmErrors*)tga_DiffXSecPlotwithSys->Clone("ratio_tota");

    if (ratio_stat) {
        ratio_stat->SetFillStyle(1001);
        ratio_stat->SetFillColor(kGray + 1);
        ratio_stat->SetLineColor(0);
        for (Int_t iter = 0; iter < tga_DiffXSecPlot->GetN(); iter++)
        {
            double binWidth = (XAxisbins[iter + 1] - XAxisbins[iter]) / 2;
            double x = tga_DiffXSecPlot->GetX()[iter];
            double y_ratio_stat = tga_DiffXSecPlot->GetY()[iter] / tga_DiffXSecPlot->GetY()[iter];
            double abserr_ratio_stat = y_ratio_stat - std::abs(tga_DiffXSecPlot->GetErrorY(iter) - tga_DiffXSecPlot->GetY()[iter]) / tga_DiffXSecPlot->GetY()[iter];
            ratio_stat->SetPoint(iter, x, y_ratio_stat);
            ratio_stat->SetPointError(iter, binWidth, binWidth, abserr_ratio_stat, abserr_ratio_stat);
        }
    }
    if (ratio_tota) {
        ratio_tota->SetFillStyle(1001);
        ratio_tota->SetFillColor(kOrange - 4);
        ratio_tota->SetLineColor(0);
        for (Int_t iter = 0; iter < tga_DiffXSecPlotwithSys->GetN(); iter++)
        {
            double binWidth = (XAxisbins[iter + 1] - XAxisbins[iter]) / 2;
            double x = tga_DiffXSecPlotwithSys->GetX()[iter];
            double y_ratio_tota = tga_DiffXSecPlotwithSys->GetY()[iter] / tga_DiffXSecPlotwithSys->GetY()[iter];
            double abserr_ratio_tota = y_ratio_tota - std::abs(tga_DiffXSecPlotwithSys->GetErrorY(iter) - tga_DiffXSecPlotwithSys->GetY()[iter]) / tga_DiffXSecPlotwithSys->GetY()[iter];
            ratio_tota->SetPoint(iter, x, y_ratio_tota);
            ratio_tota->SetPointError(iter, binWidth, binWidth, abserr_ratio_tota, abserr_ratio_tota);
        }
    }

    TLegend* leg2 = new TLegend();
    if (doClosureTest || usePoissonSmearedPseudoData) {
        leg2->AddEntry(h_DiffXSec, "Pseudo-Data", "p");
    }
    else {
        leg2->AddEntry(h_DiffXSec, "Data", "p");
    }

    setTheoryStyleAndFillLegend(h_DiffXSec, "data");
    setTheoryStyleAndFillLegend(madgraphhist, "madgraph");
    setTheoryStyleAndFillLegend(madgraphhistBinned, "madgraph", leg2);
    madgraphhistBinned->GetXaxis()->SetTitle(varhists[0]->GetXaxis()->GetTitle());
    madgraphhistBinned->GetYaxis()->SetTitle(varhists[0]->GetYaxis()->GetTitle());
    madgraphhistBinned->Draw();

    gStyle->SetErrorX(0.5);
    if (drawNLOCurves && drawMCATNLO) {
        setTheoryStyleAndFillLegend(mcnlohist, "mcatnloherwig");
        setTheoryStyleAndFillLegend(mcnlohistBinned, "mcatnloherwig", leg2);
        mcnlohistBinned->Draw("SAME");
    }
    if (drawNLOCurves && drawPOWHEG) {
        setTheoryStyleAndFillLegend(powheghist, "powhegpythia");
        setTheoryStyleAndFillLegend(powheghistBinned, "powhegpythia", leg2);
        powheghistBinned->Draw("SAME");
    }
    if (drawNLOCurves && drawPOWHEGHERWIG) {
        setTheoryStyleAndFillLegend(powhegHerwighist, "powhegherwig");
        setTheoryStyleAndFillLegend(powhegHerwighistBinned, "powhegherwig", leg2);
        powhegHerwighistBinned->Draw("SAME");
    }
    if (drawNLOCurves && drawPERUGIA11) {
        setTheoryStyleAndFillLegend(perugia11hist, "perugia11");
        setTheoryStyleAndFillLegend(perugia11histBinned, "perugia11", leg2);
        perugia11histBinned->Draw("SAME");
    }
    if (drawNLOCurves && drawKidonakis &&
        (name == "HypToppT" || name == "HypTopRapidity") &&
        !name.Contains("Lead") && !name.Contains("RestFrame") && !revokeAntiQuantityCombination) {
        setTheoryStyleAndFillLegend(Kidoth1_Binned, "kidonakis", leg2);
        Kidoth1_Binned->Draw("SAME");
    }
    if (drawNLOCurves && drawAhrens && (name == "HypTTBarMass" || name == "HypTTBarpT"))
    {
        setTheoryStyleAndFillLegend(Ahrensth1_Binned, "ahrens", leg2);
        Ahrensth1_Binned->Draw("SAME");
    }
    if (drawNLOCurves && (drawMadScaleMatching || drawMadMass)) {
        if (drawMadScaleMatching) {
            setTheoryStyleAndFillLegend(madupBinned, "scaleup", leg2);
            setTheoryStyleAndFillLegend(maddownBinned, "scaledown", leg2);
            setTheoryStyleAndFillLegend(matchupBinned, "matchup", leg2);
            setTheoryStyleAndFillLegend(matchdownBinned, "matchdown", leg2);
        }
        if (drawMadMass) {
            setTheoryStyleAndFillLegend(madupBinned, "mass173.5", leg2);
            setTheoryStyleAndFillLegend(maddownBinned, "mass171.5", leg2);
            setTheoryStyleAndFillLegend(matchupBinned, "mass175.5", leg2);
            setTheoryStyleAndFillLegend(matchdownBinned, "mass169.5", leg2);
            setTheoryStyleAndFillLegend(match2upBinned, "mass178.5", leg2);
            setTheoryStyleAndFillLegend(match2downBinned, "mass166.5", leg2);
        }
        madupBinned->Draw("SAME");
        maddownBinned->Draw("SAME");
        matchupBinned->Draw("SAME");
        matchdownBinned->Draw("SAME");
        match2upBinned->Draw("SAME");
        match2downBinned->Draw("SAME");

        ofstream OutputFileXSec;
        if (Channel == "combined" && performUnfoldingInCombinedChannel) OutputFileXSec.open(string("Plots/" + Channel + "Unfolded/" + name + "DiffXsecMass.txt"));
        else OutputFileXSec.open(string("Plots/" + Channel + "/" + name + "DiffXsecMass.txt"));

        for (int i = 1; i < madupBinned->GetNbinsX(); i++) {
            //OutputFileXSec<<"Nominal "<<"Mass 181 GeV" << " Mass 175 GeV"<< "Mass 169 GeV" << "Mass 163 GeV"<<endl;                             OutputFileXSec<<h_DiffXSec->GetBinContent(i)<< " "<<tga_DiffXSecPlot->GetErrorY(i-1)<<" "<<tga_DiffXSecPlotwithSys->GetErrorY(i-1)<< " "<<h|
        }
        OutputFileXSec.close();
    }

    if (drawNLOCurves) {
        //if (drawPOWHEG && powheghist->GetEntries())                                                     leg2->AddEntry(powheghistBinned, "tt+1jet m=178.5 GeV","l");
        //if (drawPOWHEGHERWIG && powhegHerwighist->GetEntries())                                         leg2->AddEntry(powhegHerwighistBinned, "tt+1jet m=172.5 GeV","l");
    }

    if (drawSmoothMadgraph && !revokeAntiQuantityCombination) {
        if (name.Contains("HypTTBarMass")) {
            GenPlotTheory->Rebin(15);
            GenPlotTheory->Scale(1. / GenPlotTheory->Integral("width"));
            TH1D* SmoothMadgraph = (TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(20, "R");
            SmoothMadgraph->Draw("SAME, L");
        }
        else if (name.Contains("HypBJetpTNLead")) {
            GenPlotTheory->Rebin(5);
            GenPlotTheory->Scale(1. / GenPlotTheory->Integral("width"));
            GenPlotTheory->Draw("same,c");
        }
        else if (name.Contains("HypBJetEtaLead") || name.Contains("HypBJetEtaNLead") || name.Contains("HypBJetpTLead")) {
            GenPlotTheory->Rebin(2);
            GenPlotTheory->Scale(1. / GenPlotTheory->Integral("width"));
            TH1D* SmoothMadgraph = (TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(10);
            SmoothMadgraph->Draw("SAME, L");
        }
        else if (name.Contains("HypLeptonEtaNLead")) {
            GenPlotTheory->Rebin(2);
            GenPlotTheory->Scale(1. / GenPlotTheory->Integral("width"));
            TH1D* SmoothMadgraph = (TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(5);
            SmoothMadgraph->Draw("SAME, L");
        }
        else if (name.Contains("HypLeptonEta") || name.Contains("HypTopRapidity") ||
            name.Contains("HypBJetEta") || name.Contains("HypTTBarRapidity")) {
            TH1D* SmoothMadgraph = (TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(10);
            SmoothMadgraph->Draw("SAME, L");
        }
        else if (name.Contains("HypToppT")) {
            GenPlotTheory->Rebin(6);
            GenPlotTheory->Scale(1. / GenPlotTheory->Integral("width"));
            TH1D* SmoothMadgraph = (TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(10);
            SmoothMadgraph->Draw("SAME, L");
        }
        else if (name.Contains("HypLeptonpTLead")) {
            GenPlotTheory->Rebin(2);
            GenPlotTheory->Scale(1. / GenPlotTheory->Integral("width"));
            TH1D* SmoothMadgraph = (TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(4);
            SmoothMadgraph->Draw("SAME, L");
        }
        else if (name.Contains("HypLeptonpTNLead")) {
            GenPlotTheory->Rebin(3);
            GenPlotTheory->Scale(1. / GenPlotTheory->Integral("width"));
            GenPlotTheory->Draw("same,c");
        }
        else if (name.Contains("HypTTBarpT")) {
            GenPlotTheory->Rebin(5);
            GenPlotTheory->Scale(1. / GenPlotTheory->Integral("width"));
            GenPlotTheory->Draw("same,c");
        }
        else if (name.Contains("HypLLBarpT")) {
            GenPlotTheory->Rebin(2);
            GenPlotTheory->Scale(1. / GenPlotTheory->Integral("width"));
            GenPlotTheory->Draw("same,c");
        }
        else if (name.Contains("HypLLBarMass")) {
            GenPlotTheory->Rebin(4);
            GenPlotTheory->Scale(1. / GenPlotTheory->Integral("width"));
            TH1D* SmoothMadgraph = (TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(2);
            SmoothMadgraph->Draw("same,l");
        }
        else { GenPlotTheory->Draw("SAME,C"); } //### 150512 ### 
    }

    //     madgraphhistBinned->Draw("SAME");
    DrawCMSLabels(cmsPrelimLabelMode);
    DrawDecayChLabel(channelLabel[channelType]);

    if (name.Contains("JetMult")) {
        TString legtit = "";
        if (name.Contains("pt60")) legtit += "p_{T}^{jet}> 60 GeV";
        else if (name.Contains("pt100")) legtit += "p_{T}^{jet}> 100 GeV";
        else legtit += "p_{T}^{jet}> 30 GeV, |#eta^{jet}| < 2.4";     //###########
        leg2->SetHeader(legtit);

    }

    setResultLegendStyle(leg2);
    leg2->Draw("same");


    if (drawNLOCurves && drawKidonakis && (name.Contains("ToppT") || name.Contains("TopRapidity")) && !name.Contains("Lead") && !name.Contains("RestFrame") && !revokeAntiQuantityCombination) {
        DrawLabel("(arXiv:1210.7813)", leg2->GetX1NDC() + 0.06, leg2->GetY1NDC() - 0.025, leg2->GetX2NDC(), leg2->GetY1NDC(), 12, 0.025);
    }
    if (drawNLOCurves && drawAhrens && (name == "HypTTBarMass" || name == "HypTTBarpT")) {
        DrawLabel("(arXiv:1003.5827)", leg2->GetX1NDC() + 0.06, leg2->GetY1NDC() - 0.025, leg2->GetX2NDC(), leg2->GetY1NDC(), 12, 0.025);
    }

    //     madgraphhistBinned->Draw("SAME");
    gStyle->SetEndErrorSize(10);
    tga_DiffXSecPlot->Draw("p, SAME");
    tga_DiffXSecPlotwithSys->Draw("p, SAME, Z");
    gPad->RedrawAxis();

    if (drawPlotRatio) {
        TH1D* tmpKido = 0, * tmpAhrens = 0;
        if (Kidoth1_Binned) {
            /// Idiot definition of temporary histogram for Kidonakis due to the larger number of bins in raw histogram
            tmpKido = (TH1D*)h_DiffXSec->Clone();
            tmpKido->SetLineColor(Kidoth1_Binned->GetLineColor());
            tmpKido->SetLineStyle(Kidoth1_Binned->GetLineStyle());
            tmpKido->SetLineWidth(Kidoth1_Binned->GetLineWidth());
            for (int i = 0; i < (int)tmpKido->GetNbinsX() + 2; i++) { tmpKido->SetBinContent(i, Kidoth1_Binned->GetBinContent(Kidoth1_Binned->FindBin(tmpKido->GetBinCenter(i)))); };
        }
        if (Ahrensth1_Binned) {
            /// Idiot definition of temporary histogram for Ahrens due to the larger number of bins in raw histogram
            tmpAhrens = (TH1D*)h_DiffXSec->Clone();
            tmpAhrens->SetLineColor(Ahrensth1_Binned->GetLineColor());
            tmpAhrens->SetLineStyle(Ahrensth1_Binned->GetLineStyle());
            tmpAhrens->SetLineWidth(Ahrensth1_Binned->GetLineWidth());
            for (int i = 0; i < (int)tmpAhrens->GetNbinsX() + 2; i++) { tmpAhrens->SetBinContent(i, Ahrensth1_Binned->GetBinContent(Ahrensth1_Binned->FindBin(tmpAhrens->GetBinCenter(i)))); };
        }

        double yminRatio, ymaxRatio;
        setResultRatioRanges(yminRatio, ymaxRatio);
        common::drawRatioXSEC(h_DiffXSec, madgraphhistBinned, ratio_stat, ratio_tota, powheghistBinned, mcnlohistBinned, tmpKido, tmpAhrens, powhegHerwighistBinned, perugia11histBinned, yminRatio, ymaxRatio);
        c->Update();
        c->Modified();
    };

    if (drawNLOCurves && drawMadScaleMatching) common::drawRatioXSEC(h_DiffXSec, madgraphhistBinned, 0, 0, madupBinned, maddownBinned, matchupBinned, matchdownBinned, 0, 0, 0.4, 1.6);
    if (drawNLOCurves && drawMadMass) common::drawRatioXSEC(madgraphhistBinned, h_DiffXSec, ratio_stat, 0, madupBinned, maddownBinned, matchupBinned, matchdownBinned, match2upBinned, match2downBinned, 0.4, 1.6);


    c->Print(outdir.Copy() + "DiffXS_" + name + ".eps");
    c->Print(outdir.Copy() + "DiffXS_" + name + ".pdf");
    //c->Print(outdir.Copy()+"DiffXS_"+name+".C");
    TFile out_source(outdir.Copy() + "DiffXS_" + name + "_source.root", "RECREATE");
    c->Write("canvas");
    tga_DiffXSecPlot->Write("data_staterror_only");
    tga_DiffXSecPlotwithSys->Write("data");
    h_GenDiffXSec->Write("mc");
    out_source.Close();
    delete c;
    gStyle->SetEndErrorSize(0);


    FILE* file;
    file = fopen(outdir.Copy() + name + "_Chi2Values.txt", "w");
    fprintf(file, "Variable: %s  Channel: %s \n", name.Data(), subfolderChannel.Copy().Remove(0, 1).Data());
    fprintf(file, "Theory & $\\chi^{2}/ndof$ \\\\ \n");
    fprintf(file, "\\hline \n");
    fprintf(file, "MadGraph+Pythia & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlotwithSys, madgraphhistBinned));
    if (drawNLOCurves && drawPOWHEG && powheghistBinned && powheghistBinned->GetEntries()) {
        fprintf(file, "PowHeg+Pythia & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlotwithSys, powheghistBinned));
    }
    if (drawNLOCurves && drawPOWHEGHERWIG && powhegHerwighistBinned && powhegHerwighistBinned->GetEntries()) {
        fprintf(file, "PowHeg+Herwig & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlotwithSys, powhegHerwighistBinned));
    }
    if (drawNLOCurves && drawMCATNLO && mcnlohistBinned && mcnlohistBinned->GetEntries()) {
        fprintf(file, "MC\\@NLO+Herwig & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlotwithSys, mcnlohistBinned));
    }
    if (drawNLOCurves && drawKidonakis && Kidoth1_Binned && (name.Contains("ToppT") || name.Contains("TopRapidity")) && !name.Contains("Lead") && !name.Contains("RestFrame") && !revokeAntiQuantityCombination) {
        fprintf(file, "Approx. NNLO & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlotwithSys, Kidoth1_Binned));
    }
    if (drawAhrens && Ahrensth1_Binned && (name.Contains("TTBarpT") || name.Contains("TTBarMass"))) {
        fprintf(file, "NLO+NNLL & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlotwithSys, Ahrensth1_Binned));
    }
    if (drawNLOCurves && drawMadScaleMatching) {
        fprintf(file, "Q^{2} Up & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlotwithSys, madupBinned));
        fprintf(file, "Q^{2} Down & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlotwithSys, maddownBinned));
        fprintf(file, "ME/PS Up & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlotwithSys, matchupBinned));
        fprintf(file, "ME/PS Down & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlotwithSys, matchdownBinned));
    }

    PrintResultTotxtFile(Channel, binCenters, tga_DiffXSecPlot, tga_DiffXSecPlotwithSys);
    //PrintResultPlusGenTotxtFile(Channel, binCenters, tga_DiffXSecPlot, tga_DiffXSecPlotwithSys, madgraphhistBinned, powheghistBinned, powhegHerwighistBinned, mcnlohistBinned);


    TCanvas* c1 = new TCanvas("DiffXS", "DiffXS");
    TH1* stacksum = common::summedStackHisto(stack);

    for (unsigned int i = 1; i < hists.size(); i++) { // sum all data plots to first histogram
        if (legends.at(i) == legends.at(0)) {
            varhists[0]->Add(varhists[i]);
        }
    }
    TH1D* syshist = 0;
    syshist = (TH1D*)stacksum->Clone();
    double lumierr = lumiError;
    double brerr = BRError;
    //stat uncertainty::make a function 
    for (Int_t i = 0; i <= syshist->GetNbinsX(); ++i) {
        Double_t binc = 0;
        binc += stacksum->GetBinContent(i);
        syshist->SetBinContent(i, binc);
        // calculate uncertainty: lumi uncertainty
        Double_t binerr2 = binc * binc * lumierr * lumierr;
        Double_t topunc = 0; // uncertainty on top xsec
        Double_t bin_brerr2 = binc * binc * brerr * brerr;

        double topxsecErr2 = 2.2 * 2.2 + 11.6 * 11.6;

        double topRelUnc = TMath::Sqrt(topxsecErr2) / topxsec;
        //Functionality for multiple signal histograms
        topunc += varhists[signalHist]->GetBinContent(i) * topRelUnc;
        binerr2 += (topunc * topunc);
        if (absoluteXsec && !particleXsec) binerr2 += bin_brerr2;
        syshist->SetLineColor(1);
        syshist->SetBinError(i, TMath::Sqrt(binerr2));
    }

    setControlPlotLegendStyle(varhistsPlotting, legends, leg);
    syshist->SetFillStyle(3004);
    syshist->SetFillColor(kBlack);
    //leg->AddEntry( syshist, "Uncertainty", "f" );


    varhists[0]->SetMaximum(1.5 * varhists[0]->GetBinContent(varhists[0]->GetMaximumBin()));

    varhists[0]->SetMinimum(0);
    varhists[0]->GetYaxis()->SetTitle("events");
    varhists[0]->GetXaxis()->SetNoExponent(kTRUE);
    varhists[0]->Draw("e");

    //Add the binwidth to the yaxis in yield plots
    TString ytitle = TString(varhists[0]->GetYaxis()->GetTitle()).Copy();
    double binwidth = varhists[0]->GetXaxis()->GetBinWidth(1);
    std::ostringstream width;
    width << binwidth;
    if (name.Contains("Rapidity") || name.Contains("Eta") || name.Contains("Fraction")) { ytitle.Append(" / ").Append(width.str()); }
    else if (name.Contains("pT", TString::kIgnoreCase) || name.Contains("Mass", TString::kIgnoreCase) || name.Contains("MET") || name.Contains("HT")) { ytitle.Append(" / ").Append(width.str()).Append(" GeV"); };
    varhists[0]->GetYaxis()->SetTitle(ytitle);

    stack->Draw("same HIST");

    //Only necessary if we want error bands

    /*    TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0.5)");//this is frustrating and stupid but apparently necessary...
    setex1->Draw();
    syshist->Draw("same,E2");
    TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.)");
    setex2->Draw();*/
    varhists[0]->Draw("same, e1"); //############
    //varhists[0]->Draw("same, e"); 
    DrawCMSLabels(cmsPrelimLabelMode);
    DrawDecayChLabel(channelLabel[channelType]);
    leg->Draw("SAME");
    gPad->RedrawAxis();

    c1->Print(outdir.Copy() + "preunfolded_" + name + ".eps");
    c1->Print(outdir.Copy() + "preunfolded_" + name + ".pdf");
    TFile out_root(outdir.Copy() + "preunfolded_" + name + "_source.root", "RECREATE");

    varhists[0]->Write(name + "_data");
    stacksum->Write(name + "_allmc");
    c1->Write(name + "_canvas");
    out_root.Close();
    c1->Clear();
    delete c1;
    delete stacksum;
    for (unsigned int i = 0; i < hists.size(); i++) {
        delete varhists[i];
        //         delete varhistsPlotting.at(i);
    }
}


void Plotter13TeV::PlotSingleDiffXSec(TString Channel, TString Systematic) {

    setDataSet(Channel, Systematic);
    if (!fillHisto()) return;

    if (Channel == "ee") { channelType = 0; }
    if (Channel == "mumu") { channelType = 1; }
    if (Channel == "emu") { channelType = 2; }
    if (Channel == "combined") { channelType = 3; }

    TH1::AddDirectory(kFALSE);
    TGaxis::SetMaxDigits(2);

    double Xbins[XAxisbins.size()];
    double binCenters[XAxisbinCenters.size()];

    for (unsigned int i = 0; i < XAxisbins.size();i++) { Xbins[i] = XAxisbins[i]; }
    for (unsigned int i = 0; i < XAxisbinCenters.size();i++) { binCenters[i] = XAxisbins[i] + (XAxisbins[i + 1] - XAxisbins[i]) / 2; }

    double DataSum[XAxisbinCenters.size()];
    double GenSignalSum[XAxisbinCenters.size()];
    double BGSum[XAxisbinCenters.size()];

    TH1* varhists[hists.size()];
    TString newname = name;
    if (name.Contains("Hyp")) {//Histogram naming convention has to be smarter
        newname.ReplaceAll("Hyp", 3, "", 0);
    }

    TString ftemp = "preunfolded/" + Systematic + "/" + Channel + "/" + name + "_UnfoldingHistos.root";
    TH1* GenPlotTheory = NULL;
    if (particleXsec) GenPlotTheory = fileReader->GetClone<TH1D>(ftemp, "aPseudoHist");
    else GenPlotTheory = fileReader->GetClone<TH1D>(ftemp, "aGenHist");
    TH1* GenPlotTheoryForEffAndAccCalculation = fileReader->GetClone<TH1D>(ftemp, "aGenHist");
    TH2* genReco2d = fileReader->GetClone<TH2D>(ftemp, "aRespHist");

    for (unsigned int i = 0; i < hists.size(); i++) {
        varhists[i] = hists[i].Rebin(bins, "varhists", Xbins);
        setStyle(varhists[i], i);
    }

    if (particleXsec) {
        TH1* histsOOS[histsTtBgrOutOfSpace.size()];
        for (unsigned int j = 0; j < histsTtBgrOutOfSpace.size(); j++) {
            TString tmp = histsTtBgrOutOfSpace[j].GetName();
            histsOOS[j] = histsTtBgrOutOfSpace[j].Rebin(bins, tmp, Xbins);
            removeNegativeExpectationsInBins(histsOOS[j]);
            for (unsigned int i = 0; i < hists.size(); i++) {
                if (legends.at(i) == "t#bar{t} signal" && dataset.at(i).Contains("emu_ttbarsignalplustau") && tmp.BeginsWith("emu_")) varhists[i]->Add(histsOOS[j], -1.0);
                if (legends.at(i) == "t#bar{t} other" && dataset.at(i).Contains("emu_ttbarbgviatau") && tmp.BeginsWith("emu_")) varhists[i]->Add(histsOOS[j]);
                if (legends.at(i) == "t#bar{t} signal" && dataset.at(i).Contains("ee_ttbarsignalplustau") && tmp.BeginsWith("ee_")) varhists[i]->Add(histsOOS[j], -1.0);
                if (legends.at(i) == "t#bar{t} other" && dataset.at(i).Contains("ee_ttbarbgviatau") && tmp.BeginsWith("ee_")) varhists[i]->Add(histsOOS[j]);
                if (legends.at(i) == "t#bar{t} signal" && dataset.at(i).Contains("mumu_ttbarsignalplustau") && tmp.BeginsWith("mumu_")) varhists[i]->Add(histsOOS[j], -1.0);
                if (legends.at(i) == "t#bar{t} other" && dataset.at(i).Contains("mumu_ttbarbgviatau") && tmp.BeginsWith("mumu_")) varhists[i]->Add(histsOOS[j]);
            }
        }
    }

    std::unique_ptr<TH1> GenPlot{ GenPlotTheory->Rebin(bins,"genplot",Xbins) };
    std::unique_ptr<TH1> GenPlotForEffAndAccCalculation{ GenPlotTheoryForEffAndAccCalculation->Rebin(bins,"genplotForEffAndAcc",Xbins) };

    THStack* stack = new THStack("def", "def");
    TLegend* leg = new TLegend();
    int legchange = 0;
    std::vector<TH1*> varhistsPlotting;
    varhistsPlotting.resize(hists.size());


    for (unsigned int i = 0; i < hists.size(); i++) { // prepare histos and leg
        setStyle(varhists[i], i);
        varhistsPlotting[i] = (TH1*)varhists[i]->Clone();
        if (legends.at(i) != "Data") {
            if ((legends.at(i) == DYEntry)) {
                varhists[i]->Scale(DYScale.at(channelType));
                varhistsPlotting[i]->Scale(DYScale.at(channelType));
            }

            if (i != (hists.size() - 1)) {
                if (legends.at(i) != legends.at(i + 1)) { varhistsPlotting[i]->SetLineColor(1); }
            }
            else {
                varhistsPlotting[i]->SetLineColor(1);
            }

            if (legends.at(i) != legends.at(i - 1)) {
                varhistsPlotting[i]->SetLineColor(1);
                stack->Add(varhistsPlotting[i]);
            }
            if (i > 0) {
                if (legends.at(i) != legends.at(i - 1)) {
                    legchange = i;
                    if ((legends.at(i) == DYEntry) && DYScale.at(channelType) != 1) {
                        leg->AddEntry(varhistsPlotting[i], legends.at(i), "f");
                    }
                    else {
                        leg->AddEntry(varhistsPlotting[i], legends.at(i), "f");
                    }
                }
                else {
                    varhistsPlotting[legchange]->Add(varhistsPlotting[i]);
                }
            }
        }
        else {
            if (i == 0) leg->AddEntry(varhistsPlotting[i], legends.at(i), "pe");
        }
    }

    ///////////////////////////////////
    //purity and stability plots as taken from CombinedCrossSection...

    TH1* genHist = (TH1*)GenPlotForEffAndAccCalculation->Clone();
    TH1* fidHist = NULL;
    if (particleXsec) fidHist = (TH1*)GenPlot->Clone();
    TH1* genRecHist = new TH1D("", "", bins, Xbins);
    int intbinsX[XAxisbins.size()];
    int intbinsY[XAxisbins.size()];

    // fill the elements of the main diagonal of the 2d hist into binned 1D histogram
    for (unsigned int i = 0; i < XAxisbins.size(); ++i) {
        intbinsX[i] = genReco2d->GetXaxis()->FindBin(Xbins[i] + 0.001);
        intbinsY[i] = genReco2d->GetYaxis()->FindBin(Xbins[i] + 0.001);

        if (i > 0) genRecHist->SetBinContent(i,
            ((TH2D*)genReco2d)->Integral(intbinsX[i - 1], intbinsX[i] - 1,
                intbinsY[i - 1], intbinsY[i] - 1)
        );

    }

    TH1* genPseHist = ((TH2D*)genReco2d)->ProjectionY();
    TH1* recPseHist = ((TH2D*)genReco2d)->ProjectionX();

    TH1* genBinHist = genPseHist->Rebin(bins, "genBinHist", Xbins);
    TH1* recBinHist = recPseHist->Rebin(bins, "recBinHist", Xbins);

    genRecHist->SetBinContent(0, 0);
    genRecHist->SetBinContent(bins + 1, 0);
    genBinHist->SetBinContent(0, 0);
    genBinHist->SetBinContent(bins + 1, 0);
    recBinHist->SetBinContent(0, 0);
    recBinHist->SetBinContent(bins + 1, 0);
    genHist->SetBinContent(0, 0);
    genHist->SetBinContent(bins + 1, 0);
    if (particleXsec) fidHist->SetBinContent(0, 0);
    if (particleXsec) fidHist->SetBinContent(bins + 1, 0);

    // this is realy ugly but necessary:
    // As it seems, somewhere a double is tranformed into a float so that
    // efficiencies can be larger than 1.
    for (Int_t i = 1; i <= genRecHist->GetNbinsX(); ++i) {
        if (genRecHist->GetBinContent(i) > recBinHist->GetBinContent(i)) {
            genRecHist->SetBinContent(i, recBinHist->GetBinContent(i));
            std::cout << "WARNING in PlotDifferentialCrossSections: number of events generated and reconstructed in bin" << i
                << " = " << genRecHist->GetBinContent(i) << " is larger than number of reconstructed events in that bin"
                << " = " << recBinHist->GetBinContent(i) << std::endl;
        }
        if (genRecHist->GetBinContent(i) > genBinHist->GetBinContent(i)) {
            genRecHist->SetBinContent(i, genBinHist->GetBinContent(i));
            std::cout << "WARNING in PlotDifferentialCrossSections: number of events generated and reconstructed in bin " << i
                << " is larger than number of generated events in that bin" << std::endl;
        }
    }

    // efficiency, purity, stability
    TGraphAsymmErrors* grE; // for efficiency
    TGraphAsymmErrors* grP; // for purity
    TGraphAsymmErrors* grS; // for stability
    TGraphAsymmErrors* grA = NULL; // for acceptance   // OZ 15.01.2017 initialized with NULL to eliminate warning 'may be used uninitialized'

    Bool_t effAboveOne = false;
    if (particleXsec) {
        for (Int_t i = 0; i <= fidHist->GetNbinsX() + 1; ++i) {
            if (fidHist->GetBinContent(i) < recBinHist->GetBinContent(i)) effAboveOne = true;
        }
    }

    // efficiency
    if (particleXsec && effAboveOne) grE = new TGraphAsymmErrors(recBinHist, fidHist, "pois");
    else if (particleXsec && !effAboveOne) grE = new TGraphAsymmErrors(recBinHist, fidHist);
    else grE = new TGraphAsymmErrors(recBinHist, genHist);
    grE->SetMinimum(0);
    grE->SetMaximum(1);
    grE->SetLineColor(8);
    grE->SetLineWidth(2);
    grE->SetMarkerSize(2);
    grE->SetMarkerStyle(21);
    grE->SetMarkerColor(8);

    // purity
    grP = new TGraphAsymmErrors(genRecHist, recBinHist);
    grP->SetLineColor(4);
    grP->SetLineWidth(2);
    grP->SetMarkerSize(2);
    grP->SetMarkerStyle(23);
    grP->SetMarkerColor(4);

    // stability
    grS = new TGraphAsymmErrors(genRecHist, genBinHist);
    grS->SetLineColor(2);
    grS->SetLineWidth(2);
    grS->SetMarkerSize(2);
    grS->SetMarkerStyle(22);
    grS->SetMarkerColor(2);

    // acceptance
    if (particleXsec) {
        grA = new TGraphAsymmErrors(fidHist, genHist);
        grA->SetLineColor(kViolet + 2);
        grA->SetLineWidth(2);
        grA->SetMarkerSize(2);
        grA->SetMarkerStyle(20);
        grA->SetMarkerColor(kViolet + 2);
    }


    grE->GetXaxis()->SetTitle(XAxis);
    TCanvas* cESP = new TCanvas("ESP", "ESP");

    // this is a dummy to get the x axis range corrct

    recBinHist->Reset();
    recBinHist->Draw();
    recBinHist->SetMaximum(1.);
    if (effAboveOne) recBinHist->SetMaximum(1.5);
    recBinHist->GetXaxis()->SetNoExponent(kTRUE);
    if (name.Contains("pT") || name.Contains("Mass")) {
        recBinHist->GetXaxis()->SetTitle(XAxis.Copy().Prepend("Reconstructed ").Append(" #left[GeV#right]"));
        if (name.Contains("Rapidity")) recBinHist->GetXaxis()->SetTitle(XAxis.Copy().Prepend("Reconstructed "));
    }
    else recBinHist->GetXaxis()->SetTitle(XAxis.Copy().Prepend("Reconstructed "));
    DrawCMSLabels(cmsPrelimLabelMode);
    DrawDecayChLabel(channelLabel[channelType]);
    grE->Draw("P,SAME");
    grP->Draw("P,SAME");
    grS->Draw("P,SAME");
    //if(particleXsec) grA->Draw("P,SAME"); 

    TLegend* leg3 = new TLegend();
    leg3->AddEntry(grE, "Efficiency", "p");
    leg3->AddEntry(grP, "Purity", "p");
    leg3->AddEntry(grS, "Stability", "p");
    //if(particleXsec) leg3->AddEntry(grA, "Particle acceptance", "p" );   
    setResultLegendStyle(leg3, 0);
    leg3->Draw("SAME");

    TString outdir;
    if (Channel == "combined" && performUnfoldingInCombinedChannel) outdir = ttbar::assignFolder(outpathPlots, "combinedUnfolded", Systematic);
    else outdir = ttbar::assignFolder(outpathPlots, Channel, Systematic);

    cESP->Print(outdir.Copy() + "ESP_" + name + ".eps");
    cESP->Print(outdir.Copy() + "ESP_" + name + ".pdf");
    cESP->Clear();
    delete cESP;

    bool init = false;
    for (unsigned int hist = 0; hist < hists.size(); hist++) {
        if (legends.at(hist) == "Data") {
            for (Int_t bin = 0; bin < bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
                DataSum[bin] += varhists[hist]->GetBinContent(bin + 1);
            }
        }
        else if ((legends.at(hist) == "t#bar{t} signal") && init == false) {
            signalHist = hist;
            init = true;
            for (Int_t bin = 0; bin < bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
                GenSignalSum[bin] += GenPlot->GetBinContent(bin + 1);
            }
        }
        else {
            for (Int_t bin = 0; bin < bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
                BGSum[bin] += varhists[hist]->GetBinContent(bin + 1);
            }
        }
    }
    double totalDataSum = 0;
    double GenDiffXSecPlot[XAxisbinCenters.size()];
    for (Int_t bin = 0; bin < bins; ++bin) {
        totalDataSum += DataSum[bin];
    }

    TH1* h_DiffXSec = (TH1D*)varhists[0]->Clone(); h_DiffXSec->Reset();
    TH1* h_GenDiffXSec = (TH1D*)varhists[0]->Clone(); h_GenDiffXSec->Reset();

    double DiffXSecPlot[XAxisbinCenters.size()];
    double DiffXSecStatErrorPlot[XAxisbinCenters.size()];

    TString Dummy;
    //Read central and absolute statistical uncertainty values from Nominal
    ifstream ResultsList;
    if (Channel == "combined" && performUnfoldingInCombinedChannel) ResultsList.open("UnfoldingResults/" + Systematic + "/" + Channel + "Unfolded/" + name + "Results.txt");
    else ResultsList.open("UnfoldingResults/" + Systematic + "/" + Channel + "/" + name + "Results.txt");
    if (!ResultsList.is_open())
    {
        std::cout << "WARNING (in PlotSingleDiffXSec): File is not open.\nFix this. \nEXITING!!" << std::endl;
        exit(123);
    }
    for (Int_t bin = 0; bin < bins; bin++) {//Retrieve arrays for plotting
        ResultsList >> Dummy >> XAxisbinCenters[bin] >> Dummy >> Xbins[bin] >> Dummy >> Xbins[bin + 1] >> Dummy >> DiffXSecPlot[bin] >> Dummy >> DiffXSecStatErrorPlot[bin] >> Dummy >> GenDiffXSecPlot[bin];
        h_DiffXSec->SetBinContent(bin + 1, DiffXSecPlot[bin]);
        h_DiffXSec->SetBinError(bin + 1, DiffXSecStatErrorPlot[bin]);
        h_GenDiffXSec->SetBinContent(bin + 1, GenDiffXSecPlot[bin]);
    }

    //double TotalVisXSection = 1.; // this can currently be set to 1. Re-normalization, if needed, is maintained through SVD functions. Can be changed in case of need.

    //h_DiffXSec->Scale(1/TotalVisXSection);

    double xsecnorm = 1.0;
    double TotalVisXSection = 1.; // this can currently be set to 1. Re-normalization, if needed, is maintained through SVD functions. Can be changed in case of need.

    if (!absoluteXsec) TotalVisXSection = h_DiffXSec->Integral("width");
    else {
        xsecnorm = h_DiffXSec->Integral("width");
        if ((name.Contains("Lepton") || name.Contains("Top") || name.Contains("BJet")) && !name.Contains("Lead") && !revokeAntiQuantityCombination) {
            TotalVisXSection = 2.0;
            xsecnorm = 0.5 * xsecnorm;
        }
    }
    h_DiffXSec->Scale(1 / TotalVisXSection);


    for (Int_t bin = 0; bin < bins; bin++) {//Retrieve arrays for plotting
        DiffXSecPlot[bin] = h_DiffXSec->GetBinContent(bin + 1);
        DiffXSecStatErrorPlot[bin] = h_DiffXSec->GetBinError(bin + 1);
    }

    Double_t mexl[XAxisbinCenters.size()];
    Double_t mexh[XAxisbinCenters.size()];
    for (unsigned int j = 0; j < XAxisbinCenters.size();j++) { mexl[j] = 0;mexh[j] = 0; }
    TGraphAsymmErrors* tga_DiffXSecPlot = new TGraphAsymmErrors(bins, binCenters, DiffXSecPlot, mexl, mexh, DiffXSecStatErrorPlot, DiffXSecStatErrorPlot);
    tga_DiffXSecPlot->SetMarkerStyle(1);
    //tga_DiffXSecPlot->SetMarkerStyle(20);
    tga_DiffXSecPlot->SetMarkerColor(kBlack);
    tga_DiffXSecPlot->SetMarkerSize(1);
    tga_DiffXSecPlot->SetLineWidth(2);
    tga_DiffXSecPlot->SetLineColor(kBlack);

    GenPlot->Scale(topxsec / (SignalEventswithWeight * GenPlot->GetBinWidth(1)));
    if ((name.Contains("Lepton") || name.Contains("Top") || name.Contains("BJet")) &&
        !name.Contains("Lead") && !revokeAntiQuantityCombination) {
        GenPlotTheory->Scale(1. / 2.);
    }

    GenPlotTheory->Scale(1. / GenPlotTheory->Integral("width"));
    h_GenDiffXSec->Scale(1. / h_GenDiffXSec->Integral("width"));

    TH1* madgraphhist = 0, * mcnlohist = 0, * mcnlohistup = 0, * mcnlohistdown = 0, * powheghist = 0, * powhegHerwighist = 0, * perugia11hist = 0;
    TH1* mcnlohistnorm = 0;
    TGraph* mcatnloBand = 0;

    TH1* madgraphhistBinned = 0, * mcnlohistnormBinned = 0, * mcnlohistupBinned = 0;
    TH1* mcnlohistdownBinned = 0, * mcnlohistBinned = 0;
    TH1* powheghistBinned = 0, * powhegHerwighistBinned = 0;
    TH1* perugia11histBinned = 0;

    TH1* Kidoth1_Binned = 0, * Ahrensth1_Binned = 0;

    TH1* madup = 0, * maddown = 0, * matchup = 0, * matchdown = 0, * match2up = 0, * match2down = 0;
    TH1* madupBinned = 0, * maddownBinned = 0, * matchupBinned = 0, * matchdownBinned = 0, * match2upBinned = 0, * match2downBinned = 0;

    bool canDrawMCATNLO = true;
    if (drawNLOCurves && drawMCATNLO) {
        mcnlohist = GetNloCurve(newname, "MCATNLO");
        mcnlohistBinned = RescaleNLOCurveAndPrepareBinnedCurve(mcnlohist, bins, Xbins, xsecnorm, BranchingFraction[channelType]);

        if (name.Contains("LeptonpT")) { mcnlohistnorm = GetNloCurve("Leptons", "Pt", "MCatNLO"); }//temprorary until I change the naming convention in the root file
        else if (name.Contains("LeptonEta")) { mcnlohistnorm = GetNloCurve("Leptons", "Eta", "MCatNLO"); }
        else if (name.Contains("LLBarpT")) { mcnlohistnorm = GetNloCurve("LepPair", "Pt", "MCatNLO"); }
        else if (name.Contains("LLBarMass")) { mcnlohistnorm = GetNloCurve("LepPair", "Mass", "MCatNLO"); }
        else if (name.Contains("ToppT")) { mcnlohistnorm = GetNloCurve("TopQuarks", "Pt", "MCatNLO"); }
        else if (name.Contains("TopRapidity")) { mcnlohistnorm = GetNloCurve("TopQuarks", "Rapidity", "MCatNLO"); }
        else if (name.Contains("TTBarpT")) { mcnlohistnorm = GetNloCurve("TtBar", "Pt", "MCatNLO"); }
        else if (name.Contains("TTBarRapidity")) { mcnlohistnorm = GetNloCurve("TtBar", "Rapidity", "MCatNLO"); }
        else if (name.Contains("TTBarMass")) { mcnlohistnorm = GetNloCurve("TtBar", "Mass", "MCatNLO"); }
        else if (name.Contains("BJetpT")) { mcnlohistnorm = GetNloCurve("Jets", "Pt", "MCatNLO"); }
        else if (name.Contains("BJetEta")) { mcnlohistnorm = GetNloCurve("Jets", "Eta", "MCatNLO"); }

        if (mcnlohistnorm) {
            mcnlohistnormBinned = mcnlohistnorm->Rebin(bins, "genBinHistNorm", Xbins);

            if (name.Contains("LeptonpT")) { mcnlohistup = GetNloCurve("Leptons", "Pt", "MCNLOup"); }//temprorary until I change the naming convention in the root file
            else if (name.Contains("LeptonEta")) { mcnlohistup = GetNloCurve("Leptons", "Eta", "MCNLOup"); }
            else if (name.Contains("LLBarpT")) { mcnlohistup = GetNloCurve("LepPair", "Pt", "MCNLOup"); }
            else if (name.Contains("LLBarMass")) { mcnlohistup = GetNloCurve("LepPair", "Mass", "MCNLOup"); }
            else if (name.Contains("ToppT")) { mcnlohistup = GetNloCurve("TopQuarks", "Pt", "MCNLOup"); }
            else if (name.Contains("TopRapidity")) { mcnlohistup = GetNloCurve("TopQuarks", "Rapidity", "MCNLOup"); }
            else if (name.Contains("TTBarpT")) { mcnlohistup = GetNloCurve("TtBar", "Pt", "MCNLOup"); }
            else if (name.Contains("TTBarRapidity")) { mcnlohistup = GetNloCurve("TtBar", "Rapidity", "MCNLOup"); }
            else if (name.Contains("TTBarMass")) { mcnlohistup = GetNloCurve("TtBar", "Mass", "MCNLOup"); }
            else if (name.Contains("BJetpT")) { mcnlohistup = GetNloCurve("Jets", "Pt", "MCNLOup"); }
            else if (name.Contains("BJetEta")) { mcnlohistup = GetNloCurve("Jets", "Eta", "MCNLOup"); }
            mcnlohistupBinned = mcnlohistup->Rebin(bins, "genBinHist", Xbins);


            if (name.Contains("LeptonpT")) { mcnlohistdown = GetNloCurve("Leptons", "Pt", "MCNLOdown"); }//temprorary until I change the naming convention in the root file
            else if (name.Contains("LeptonEta")) { mcnlohistdown = GetNloCurve("Leptons", "Eta", "MCNLOdown"); }
            else if (name.Contains("LLBarpT")) { mcnlohistdown = GetNloCurve("LepPair", "Pt", "MCNLOdown"); }
            else if (name.Contains("LLBarMass")) { mcnlohistdown = GetNloCurve("LepPair", "Mass", "MCNLOdown"); }
            else if (name.Contains("ToppT")) { mcnlohistdown = GetNloCurve("TopQuarks", "Pt", "MCNLOdown"); }
            else if (name.Contains("TopRapidity")) { mcnlohistdown = GetNloCurve("TopQuarks", "Rapidity", "MCNLOdown"); }
            else if (name.Contains("TTBarpT")) { mcnlohistdown = GetNloCurve("TtBar", "Pt", "MCNLOdown"); }
            else if (name.Contains("TTBarRapidity")) { mcnlohistdown = GetNloCurve("TtBar", "Rapidity", "MCNLOdown"); }
            else if (name.Contains("TTBarMass")) { mcnlohistdown = GetNloCurve("TtBar", "Mass", "MCNLOdown"); }
            else if (name.Contains("BJetpT")) { mcnlohistdown = GetNloCurve("Jets", "Pt", "MCNLOdown"); }
            else if (name.Contains("BJetEta")) { mcnlohistdown = GetNloCurve("Jets", "Eta", "MCNLOdown"); }
            mcnlohistdownBinned = mcnlohistdown->Rebin(bins, "genBinHist", Xbins);

            for (Int_t bin = 0; bin < bins; bin++) {
                mcnlohistupBinned->SetBinContent(bin + 1, mcnlohistupBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / mcnlohistup->GetBinWidth(1)));
                mcnlohistdownBinned->SetBinContent(bin + 1, mcnlohistdownBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / mcnlohistdown->GetBinWidth(1)));
                mcnlohistnormBinned->SetBinContent(bin + 1, mcnlohistnormBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / mcnlohistnorm->GetBinWidth(1)));
            }
            mcnlohistupBinned->Scale(xsecnorm / mcnlohistnormBinned->Integral("width"));
            mcnlohistdownBinned->Scale(xsecnorm / mcnlohistnormBinned->Integral("width"));
            mcnlohistnormBinned->Scale(xsecnorm / mcnlohistnormBinned->Integral("width"));

            for (Int_t bin = 0; bin < bins; bin++) {
                mcnlohistupBinned->SetBinContent(bin + 1, (mcnlohistupBinned->GetBinContent(bin + 1) / mcnlohistnormBinned->GetBinContent(bin + 1)) * mcnlohistBinned->GetBinContent(bin + 1));
                mcnlohistdownBinned->SetBinContent(bin + 1, (mcnlohistdownBinned->GetBinContent(bin + 1) / mcnlohistnormBinned->GetBinContent(bin + 1)) * mcnlohistBinned->GetBinContent(bin + 1));
            }

            //Uncertainty band for MC@NLO
            Double_t x[bins];
            Double_t xband[2 * bins];
            Double_t errup[bins];
            Double_t errdn[bins];
            Double_t errorband[2 * bins];

            for (Int_t j = 0; j < bins; j++) {
                x[j] = mcnlohistBinned->GetBinCenter(j + 1);
                errup[j] = (mcnlohistupBinned->GetBinContent(j + 1) / mcnlohistnormBinned->GetBinContent(j + 1)) * mcnlohistBinned->GetBinContent(j + 1);
                errdn[j] = (mcnlohistdownBinned->GetBinContent(j + 1) / mcnlohistnormBinned->GetBinContent(j + 1)) * mcnlohistBinned->GetBinContent(j + 1);

                xband[j] = x[j];
                errorband[j] = errdn[j]; //lower band
                xband[2 * bins - j - 1] = x[j];
                errorband[2 * bins - j - 1] = errup[j]; //upper band
            }

            mcatnloBand = new TGraph(2 * bins, xband, errorband);
            mcatnloBand->SetFillColor(kGray);
            mcatnloBand->SetFillStyle(1001);
            mcatnloBand->SetLineColor(kBlue);
            mcatnloBand->SetLineWidth(2);
            mcatnloBand->SetLineStyle(5);
            canDrawMCATNLO = false;
        }
        else {
            std::cout << "\n*************************\nMC@NLO Curve not available!\n**********************\n";
            canDrawMCATNLO = false;
        }
    }

    madgraphhist = GetNloCurve(newname, "Nominal");
    madgraphhistBinned = RescaleNLOCurveAndPrepareBinnedCurve(madgraphhist, bins, Xbins, xsecnorm, BranchingFraction[channelType]);

    if (drawNLOCurves && drawPOWHEG) {
        powheghist = GetNloCurve(newname, "POWHEG");
        powheghistBinned = RescaleNLOCurveAndPrepareBinnedCurve(powheghist, bins, Xbins, xsecnorm, BranchingFraction[channelType]);
    }
    if (drawNLOCurves && drawPOWHEGHERWIG) {
        powhegHerwighist = GetNloCurve(newname, "POWHEGHERWIG");
        powhegHerwighistBinned = RescaleNLOCurveAndPrepareBinnedCurve(powhegHerwighist, bins, Xbins, xsecnorm, BranchingFraction[channelType]);
    }
    if (drawNLOCurves && drawPERUGIA11) {
        perugia11hist = GetNloCurve(newname, "PERUGIA11");
        perugia11histBinned = RescaleNLOCurveAndPrepareBinnedCurve(perugia11hist, bins, Xbins, xsecnorm, BranchingFraction[channelType]);
    }
    if (drawNLOCurves && drawKidonakis &&
        (name == "HypToppT" || name == "HypTopRapidity") && !name.Contains("Lead") && !name.Contains("RestFrame") && !revokeAntiQuantityCombination) {
        TString kidoFile = ttbar::DATA_PATH_DILEPTONIC() + "/kidonakisNNLO_8TeV.root";
        if (name.Contains("ToppT"))             Kidoth1_Binned = fileReader->GetClone<TH1>(kidoFile, "topPt");
        else if (name.Contains("TopRapidity"))  Kidoth1_Binned = fileReader->GetClone<TH1>(kidoFile, "topY");
        for (int iter = 0; iter <= Kidoth1_Binned->GetNbinsX(); iter++) Kidoth1_Binned->SetBinError(iter, 1e-10);
        Kidoth1_Binned->Scale(1. / Kidoth1_Binned->Integral("width"));
    }
    if (drawNLOCurves && drawAhrens && (name == "HypTTBarMass" || name == "HypTTBarpT")) {
        TString ahrensFile = ttbar::DATA_PATH_DILEPTONIC() + "/ahrensNNLL_8TeV.root";
        if (name == "HypTTBarMass")    Ahrensth1_Binned = fileReader->GetClone<TH1>(ahrensFile, "ttbarM");
        else if (name == "HypTTBarpT") Ahrensth1_Binned = fileReader->GetClone<TH1>(ahrensFile, "ttbarPt");
        for (int iter = 0; iter <= Ahrensth1_Binned->GetNbinsX(); iter++) Ahrensth1_Binned->SetBinError(iter, 1e-10);
        Ahrensth1_Binned->Scale(1. / Ahrensth1_Binned->Integral("width"));
    }
    else { drawAhrens = 0; }

    if (drawMadScaleMatching) {
        madup = GetNloCurve(newname, "SCALE_UP");
        madupBinned = RescaleNLOCurveAndPrepareBinnedCurve(madup, bins, Xbins, xsecnorm, BranchingFraction[channelType]);

        maddown = GetNloCurve(newname, "SCALE_DOWN");
        maddownBinned = RescaleNLOCurveAndPrepareBinnedCurve(maddown, bins, Xbins, xsecnorm, BranchingFraction[channelType]);

        matchup = GetNloCurve(newname, "MATCH_UP");
        matchupBinned = RescaleNLOCurveAndPrepareBinnedCurve(matchup, bins, Xbins, xsecnorm, BranchingFraction[channelType]);

        matchdown = GetNloCurve(newname, "MATCH_DOWN");
        matchdownBinned = RescaleNLOCurveAndPrepareBinnedCurve(matchdown, bins, Xbins, xsecnorm, BranchingFraction[channelType]);

    }

    if (drawNLOCurves && drawMadMass) {
        madup = GetNloCurveMass(newname, "MASS_UP", "173.5");
        madup->Scale(1. / madup->Integral("width"));
        madupBinned = madup->Rebin(bins, "madupplot", Xbins);
        for (Int_t bin = 0; bin < bins; bin++) {//condense matrices to arrays for plotting
            madupBinned->SetBinContent(bin + 1, madupBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / madup->GetBinWidth(1)));
        }
        madupBinned->Scale(xsecnorm / madupBinned->Integral("width"));

        maddown = GetNloCurveMass(newname, "MASS_DOWN", "171.5");
        maddown->Scale(1. / maddown->Integral("width"));
        maddownBinned = maddown->Rebin(bins, "maddownplot", Xbins);
        for (Int_t bin = 0; bin < bins; bin++) {//condense matrices to arrays for plotting
            maddownBinned->SetBinContent(bin + 1, maddownBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / maddown->GetBinWidth(1)));
        }
        maddownBinned->Scale(xsecnorm / maddownBinned->Integral("width"));

        matchup = GetNloCurveMass(newname, "MASS_UP", "175.5");
        matchup->Scale(1. / matchup->Integral("width"));
        matchupBinned = matchup->Rebin(bins, "matchupplot", Xbins);
        for (Int_t bin = 0; bin < bins; bin++) {//condense matrices to arrays for plotting
            matchupBinned->SetBinContent(bin + 1, matchupBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / matchup->GetBinWidth(1)));
        }
        matchupBinned->Scale(xsecnorm / matchupBinned->Integral("width"));

        matchdown = GetNloCurveMass(newname, "MASS_DOWN", "169.5");
        matchdown->Scale(1. / matchdown->Integral("width"));
        matchdownBinned = matchdown->Rebin(bins, "matchdownplot", Xbins);
        for (Int_t bin = 0; bin < bins; bin++) {//condense matrices to arrays for plotting
            matchdownBinned->SetBinContent(bin + 1, matchdownBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / matchdown->GetBinWidth(1)));
        }
        matchdownBinned->Scale(xsecnorm / matchdownBinned->Integral("width"));

        match2up = GetNloCurveMass(newname, "MASS_UP", "178.5");
        double match2scale = 1. / match2up->Integral("width");
        match2up->Scale(match2scale);
        match2upBinned = match2up->Rebin(bins, "match2upplot", Xbins);
        for (Int_t bin = 0; bin < bins; bin++) {//condense matrices to arrays for plotting
            match2upBinned->SetBinContent(bin + 1, match2upBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / match2up->GetBinWidth(1)));
        }
        match2upBinned->Scale(xsecnorm / match2upBinned->Integral("width"));

        match2down = GetNloCurveMass(newname, "MASS_DOWN", "166.5");
        double match2downscale = 1. / match2down->Integral("width");
        match2down->Scale(match2downscale);
        match2downBinned = match2down->Rebin(bins, "match2downplot", Xbins);
        for (Int_t bin = 0; bin < bins; bin++) {//condense matrices to arrays for plotting
            match2downBinned->SetBinContent(bin + 1, match2downBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / match2down->GetBinWidth(1)));
        }
        match2downBinned->Scale(xsecnorm / match2downBinned->Integral("width"));

    }

    TCanvas* c = new TCanvas("DiffXS", "DiffXS");
    if (logY) {
        c->SetLogy();
    }
    h_DiffXSec->SetMarkerStyle(20);
    if (ymax != 0) {
        if (logY) {
            madgraphhistBinned->SetMaximum(18 * madgraphhistBinned->GetBinContent(madgraphhistBinned->GetMaximumBin()));
        }
        else {
            madgraphhistBinned->SetMaximum(1.5 * madgraphhistBinned->GetBinContent(madgraphhistBinned->GetMaximumBin()));
        }
    }
    madgraphhistBinned->GetXaxis()->SetNoExponent(kTRUE);
    if (name.Contains("Rapidity") || name.Contains("Eta")) { madgraphhistBinned->GetYaxis()->SetNoExponent(kTRUE); }
    if (ymax != 0) madgraphhistBinned->SetMaximum(ymax);
    if (ymin != 0) madgraphhistBinned->SetMinimum(ymin);

    gStyle->SetEndErrorSize(8);
    if (drawNLOCurves && drawMCATNLO && canDrawMCATNLO) {
        mcnlohistupBinned->SetFillColor(kGray);
        mcnlohistupBinned->SetLineColor(kGray);
        mcnlohistupBinned->Draw("same");
        mcnlohistdownBinned->SetLineColor(10);
        mcnlohistdownBinned->SetFillColor(10);
        mcnlohistdownBinned->Draw("same");
    }
    GenPlotTheory->SetLineColor(kRed + 1);
    GenPlotTheory->SetLineWidth(2);
    GenPlotTheory->SetLineStyle(1);
    h_GenDiffXSec->SetLineColor(GenPlotTheory->GetLineColor());
    h_GenDiffXSec->SetLineStyle(GenPlotTheory->GetLineStyle());

    // Plot statistical and totl error of full result
    TGraphAsymmErrors* ratio_stat = 0;
    if (tga_DiffXSecPlot) ratio_stat = (TGraphAsymmErrors*)tga_DiffXSecPlot->Clone("ratio_stat");

    if (ratio_stat) {
        ratio_stat->SetFillStyle(1001);
        ratio_stat->SetFillColor(kGray + 1);
        ratio_stat->SetLineColor(0);
        for (Int_t iter = 0; iter < tga_DiffXSecPlot->GetN(); iter++)
        {
            double binWidth = (XAxisbins[iter + 1] - XAxisbins[iter]) / 2;
            double x = tga_DiffXSecPlot->GetX()[iter];
            double y_ratio_stat = tga_DiffXSecPlot->GetY()[iter] / tga_DiffXSecPlot->GetY()[iter];
            double abserr_ratio_stat = y_ratio_stat - std::abs(tga_DiffXSecPlot->GetErrorY(iter) - tga_DiffXSecPlot->GetY()[iter]) / tga_DiffXSecPlot->GetY()[iter];
            ratio_stat->SetPoint(iter, x, y_ratio_stat);
            ratio_stat->SetPointError(iter, binWidth, binWidth, abserr_ratio_stat, abserr_ratio_stat);
        }
    }

    TLegend* leg2 = new TLegend();
    leg2->SetHeader(Systematic);
    if (doClosureTest || usePoissonSmearedPseudoData) {
        leg2->AddEntry(h_DiffXSec, "Pseudo-Data", "p");
    }
    else {
        leg2->AddEntry(h_DiffXSec, "Data", "p");
    }

    setTheoryStyleAndFillLegend(h_DiffXSec, "data");
    setTheoryStyleAndFillLegend(madgraphhist, "madgraph");
    setTheoryStyleAndFillLegend(madgraphhistBinned, "madgraph", leg2);
    madgraphhistBinned->GetXaxis()->SetTitle(varhists[0]->GetXaxis()->GetTitle());
    madgraphhistBinned->GetYaxis()->SetTitle(varhists[0]->GetYaxis()->GetTitle());
    madgraphhistBinned->Draw();

    TH1* realTruthBinned = nullptr;
    if (!usePoissonSmearedPseudoData && doClosureTest && dataset.at(0).Contains("fakerun") && Channel != "combined")
    {
        newname.ReplaceAll("Hyp", visTheoryPrefix);
        newname.Prepend(visTheoryPrefix);
        TH1* realTruth = fileReader->GetClone<TH1>(dataset.at(0), newname);
        if (!revokeAntiQuantityCombination && !newname.Contains("Lead") && (newname.Contains("Lepton") || newname.Contains("Top") || newname.Contains("BJet")))
        {
            //loop over anti-particle histograms and add them
            TString antiName = newname.ReplaceAll(visTheoryPrefix, visTheoryPrefix + "Anti");
            realTruth->Add(fileReader->Get<TH1>(dataset.at(0), antiName));
        }
        // Rebin histogram to final binning
        realTruth->Scale(1. / realTruth->Integral("width"));
        realTruthBinned = realTruth->Rebin(bins, "realTruthBinned", Xbins);
        for (Int_t bin = 0; bin < bins; bin++) {
            realTruthBinned->SetBinContent(bin + 1, realTruthBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / realTruth->GetBinWidth(1)));
        }
        realTruthBinned->Scale(1. / realTruthBinned->Integral("width"));

        // Histogram style
        realTruthBinned->SetLineColor(kRed - 7);
        realTruthBinned->SetLineStyle(2);
        realTruthBinned->SetLineWidth(2);
        realTruthBinned->Draw("SAME][");
        leg2->AddEntry(realTruthBinned, "Simu. Reweighted", "l");
    }

    gStyle->SetErrorX(0.5);
    if (drawNLOCurves && drawMCATNLO) {
        setTheoryStyleAndFillLegend(mcnlohist, "mcatnloherwig");
        setTheoryStyleAndFillLegend(mcnlohistBinned, "mcatnloherwig", leg2);
        mcnlohistBinned->Draw("SAME");
    }
    if (drawNLOCurves && drawPOWHEG) {
        setTheoryStyleAndFillLegend(powheghist, "powhegpythia");
        setTheoryStyleAndFillLegend(powheghistBinned, "powhegpythia", leg2);
        powheghistBinned->Draw("SAME");
    }
    if (drawNLOCurves && drawPOWHEGHERWIG) {
        setTheoryStyleAndFillLegend(powhegHerwighist, "powhegherwig");
        setTheoryStyleAndFillLegend(powhegHerwighistBinned, "powhegherwig", leg2);
        powhegHerwighistBinned->Draw("SAME");
    }
    if (drawNLOCurves && drawPERUGIA11) {
        setTheoryStyleAndFillLegend(perugia11hist, "perugia11");
        setTheoryStyleAndFillLegend(perugia11histBinned, "perugia11", leg2);
        perugia11histBinned->Draw("SAME");
    }
    if (drawNLOCurves && drawKidonakis &&
        (name == "HypToppT" || name == "HypTopRapidity") &&
        !name.Contains("Lead") && !name.Contains("RestFrame") && !revokeAntiQuantityCombination) {
        setTheoryStyleAndFillLegend(Kidoth1_Binned, "kidonakis", leg2);
        Kidoth1_Binned->Draw("SAME");
    }
    if (drawNLOCurves && drawAhrens && (name == "HypTTBarMass" || name == "HypTTBarpT"))
    {
        setTheoryStyleAndFillLegend(Ahrensth1_Binned, "ahrens", leg2);
        Ahrensth1_Binned->Draw("SAME");
    }
    if (drawNLOCurves && (drawMadScaleMatching || drawMadMass)) {
        if (drawMadScaleMatching) {
            setTheoryStyleAndFillLegend(madupBinned, "scaleup", leg2);
            setTheoryStyleAndFillLegend(maddownBinned, "scaledown", leg2);
            setTheoryStyleAndFillLegend(matchupBinned, "matchup", leg2);
            setTheoryStyleAndFillLegend(matchdownBinned, "matchdown", leg2);
        }
        if (drawMadMass) {
            setTheoryStyleAndFillLegend(madupBinned, "mass173.5", leg2);
            setTheoryStyleAndFillLegend(maddownBinned, "mass171.5", leg2);
            setTheoryStyleAndFillLegend(matchupBinned, "mass175.5", leg2);
            setTheoryStyleAndFillLegend(matchdownBinned, "mass169.5", leg2);
            setTheoryStyleAndFillLegend(match2upBinned, "mass178.5", leg2);
            setTheoryStyleAndFillLegend(match2downBinned, "mass166.5", leg2);
        }
        madupBinned->Draw("SAME");
        maddownBinned->Draw("SAME");
        matchupBinned->Draw("SAME");
        matchdownBinned->Draw("SAME");
        match2upBinned->Draw("SAME");
        match2downBinned->Draw("SAME");

        ofstream OutputFileXSec;
        if (Channel == "combined" && performUnfoldingInCombinedChannel) OutputFileXSec.open(string("Plots/" + Channel + "Unfolded/" + name + "DiffXsecMass.txt"));
        else OutputFileXSec.open(string("Plots/" + Channel + "/" + name + "DiffXsecMass.txt"));

        for (int i = 1; i < madupBinned->GetNbinsX(); i++) {
            //OutputFileXSec<<"Nominal "<<"Mass 181 GeV" << " Mass 175 GeV"<< "Mass 169 GeV" << "Mass 163 GeV"<<endl;                             OutputFileXSec<<h_DiffXSec->GetBinContent(i)<< " "<<tga_DiffXSecPlot->GetErrorY(i-1)<<" "<<tga_DiffXSecPlotwithSys->GetErrorY(i-1)<< " "<<h|
        }
        OutputFileXSec.close();
    }

    //     madgraphhistBinned->Draw("SAME");
    DrawCMSLabels(cmsPrelimLabelMode);
    DrawDecayChLabel(channelLabel[channelType]);

    if (drawNLOCurves) {
        //if (drawPOWHEG && powheghist->GetEntries())                                                     leg2->AddEntry(powheghistBinned, "tt+1jet m=178.5 GeV","l");
        //if (drawPOWHEGHERWIG && powhegHerwighist->GetEntries())                                         leg2->AddEntry(powhegHerwighistBinned, "tt+1jet m=172.5 GeV","l");
    }


    if (name.Contains("JetMult")) {
        TString legtit = "";
        if (name.Contains("pt60")) legtit += "p_{T}^{jet}> 60 GeV";
        else if (name.Contains("pt100")) legtit += "p_{T}^{jet}> 100 GeV";
        else legtit += "p_{T}^{jet}> 30 GeV, |#eta^{jet}| < 2.4";     //###########
        leg2->SetHeader(legtit);
    }

    setResultLegendStyle(leg2);
    leg2->Draw("same");

    if (drawNLOCurves && drawKidonakis && (name.Contains("ToppT") || name.Contains("TopRapidity")) && !name.Contains("Lead") && !name.Contains("RestFrame") && !revokeAntiQuantityCombination) {
        DrawLabel("(arXiv:1210.7813)", leg2->GetX1NDC() + 0.06, leg2->GetY1NDC() - 0.025, leg2->GetX2NDC(), leg2->GetY1NDC(), 12, 0.025);
    }
    if (drawNLOCurves && drawAhrens && (name == "HypTTBarMass" || name == "HypTTBarpT")) {
        DrawLabel("(arXiv:1003.5827)", leg2->GetX1NDC() + 0.06, leg2->GetY1NDC() - 0.025, leg2->GetX2NDC(), leg2->GetY1NDC(), 12, 0.025);
    }

    //     madgraphhistBinned->Draw("same");
    gStyle->SetEndErrorSize(10);
    tga_DiffXSecPlot->Draw("p, SAME");
    ///Stupid clone to be able to see the stat error bars
    TGraphAsymmErrors* tga_DiffXSecPlotClone = (TGraphAsymmErrors*)tga_DiffXSecPlot->Clone();
    tga_DiffXSecPlotClone->SetMarkerStyle(20);
    tga_DiffXSecPlotClone->Draw("p, SAME");
    gPad->RedrawAxis();

    if (drawPlotRatio) {
        TH1D* tmpKido = 0, * tmpAhrens = 0;
        if (Kidoth1_Binned) {
            /// Idiot definition of temporary histogram for Kidonakis due to the larger number of bins in raw histogram
            tmpKido = (TH1D*)h_DiffXSec->Clone();
            tmpKido->SetLineColor(Kidoth1_Binned->GetLineColor());
            tmpKido->SetLineStyle(Kidoth1_Binned->GetLineStyle());
            tmpKido->SetLineWidth(Kidoth1_Binned->GetLineWidth());
            for (int i = 0; i < (int)tmpKido->GetNbinsX() + 2; i++) { tmpKido->SetBinContent(i, Kidoth1_Binned->GetBinContent(Kidoth1_Binned->FindBin(tmpKido->GetBinCenter(i)))); };
        }
        if (Ahrensth1_Binned) {
            /// Idiot definition of temporary histogram for Ahrens due to the larger number of bins in raw histogram
            tmpAhrens = (TH1D*)h_DiffXSec->Clone();
            tmpAhrens->SetLineColor(Ahrensth1_Binned->GetLineColor());
            tmpAhrens->SetLineStyle(Ahrensth1_Binned->GetLineStyle());
            tmpAhrens->SetLineWidth(Ahrensth1_Binned->GetLineWidth());
            for (int i = 0; i < (int)tmpAhrens->GetNbinsX() + 2; i++) { tmpAhrens->SetBinContent(i, Ahrensth1_Binned->GetBinContent(Ahrensth1_Binned->FindBin(tmpAhrens->GetBinCenter(i)))); };
        }

        double yminRatio, ymaxRatio;
        setResultRatioRanges(yminRatio, ymaxRatio);
        //if(!usePoissonSmearedPseudoData && doClosureTest) { common::drawRatioXSEC(h_DiffXSec, realTruthBinned, ratio_stat , 0,madgraphhistBinned, 0,
        //                                         0, 0, 0,0, 0.4, 1.6);
        if (!usePoissonSmearedPseudoData && doClosureTest) {
            common::drawRatioXSEC(h_DiffXSec, madgraphhistBinned, ratio_stat, 0, realTruthBinned, 0,
                mcnlohistBinned, tmpKido, tmpAhrens, powhegHerwighistBinned, 0.4, 1.6);
        }
        else { //common::drawRatioXSEC(h_DiffXSec, madgraphhistBinned, ratio_stat,0,powheghistBinned, mcnlohistBinned, tmpKido, tmpAhrens, powhegHerwighistBinned, perugia11histBinned,0.4, 1.6);
            if (drawNLOCurves && drawMadScaleMatching) common::drawRatioXSEC(h_DiffXSec, madgraphhistBinned, 0, 0, madupBinned, maddownBinned, matchupBinned, matchdownBinned, 0, 0, 0.4, 1.6);
            else if (drawNLOCurves && drawMadMass) common::drawRatioXSEC(madgraphhistBinned, h_DiffXSec, 0, 0, madupBinned, maddownBinned, matchupBinned, matchdownBinned, match2upBinned, match2downBinned, 0.4, 1.6);
            else common::drawRatioXSEC(h_DiffXSec, madgraphhistBinned, ratio_stat, 0, powheghistBinned, mcnlohistBinned, tmpKido, tmpAhrens, powhegHerwighistBinned, perugia11histBinned, yminRatio, ymaxRatio);
        };
    };
    c->Update();
    c->Modified();

    c->Print(outdir.Copy() + "DiffXS_" + name + ".eps");
    c->Print(outdir.Copy() + "DiffXS_" + name + ".pdf");
    TFile out_source(outdir.Copy() + "DiffXS_" + name + "_source.root", "RECREATE");
    c->Write("canvas");
    tga_DiffXSecPlot->Write("data_staterror_only");
    h_GenDiffXSec->Write("mc");
    out_source.Close();
    delete c;
    gStyle->SetEndErrorSize(0);

    //    double result_Integral = Plotter13TeV::CalculateIntegral(tga_DiffXSecPlot, Xbins);
    FILE* file;
    file = fopen(outdir.Copy() + name + "_Chi2Values.txt", "w");
    fprintf(file, "Variable: %s  Channel: %s \n", name.Data(), subfolderChannel.Copy().Remove(0, 1).Data());
    fprintf(file, "Theory & $\\chi^{2}/ndof$ \\\\ \n");
    fprintf(file, "\\hline \n");
    fprintf(file, "MadGraph+Pythia & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlot, madgraphhistBinned));
    if (drawNLOCurves && drawPOWHEG && powheghistBinned && powheghistBinned->GetEntries()) {
        fprintf(file, "PowHeg+Pythia & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlot, powheghistBinned));
    }
    if (drawNLOCurves && drawPOWHEGHERWIG && powhegHerwighistBinned && powhegHerwighistBinned->GetEntries()) {
        fprintf(file, "PowHeg+Herwig & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlot, powhegHerwighistBinned));
    }
    if (drawNLOCurves && drawMCATNLO && mcnlohistBinned && mcnlohistBinned->GetEntries()) {
        fprintf(file, "MC\\@NLO+Herwig & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlot, mcnlohistBinned));
    }
    if (drawNLOCurves && drawKidonakis && Kidoth1_Binned && (name.Contains("ToppT") || name.Contains("TopRapidity")) && !name.Contains("Lead") && !name.Contains("RestFrame") && !revokeAntiQuantityCombination) {
        fprintf(file, "Approx. NNLO & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlot, Kidoth1_Binned));
    }
    if (drawAhrens && Ahrensth1_Binned && (name.Contains("TTBarpT") || name.Contains("TTBarMass"))) {
        fprintf(file, "NLO+NNLL & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlot, Ahrensth1_Binned));
    }
    if (drawNLOCurves && drawMadScaleMatching) {
        fprintf(file, "Q^{2} Up & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlot, madupBinned));
        fprintf(file, "Q^{2} Down & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlot, maddownBinned));
        fprintf(file, "ME/PS Up & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlot, matchupBinned));
        fprintf(file, "ME/PS Down & %3.2f \\\\ \n", Plotter13TeV::GetChi2(tga_DiffXSecPlot, matchdownBinned));
    }

    TCanvas* c1 = new TCanvas("DiffXS", "DiffXS");
    TH1* stacksum = common::summedStackHisto(stack);

    for (unsigned int i = 1; i < hists.size(); i++) { // sum all data plots to first histogram
        if (legends.at(i) == legends.at(0)) {
            varhists[0]->Add(varhists[i]);
        }
    }
    TH1D* syshist = 0;
    syshist = (TH1D*)stacksum->Clone();
    double lumierr = lumiError;
    double brerr = BRError;
    //stat uncertainty::make a function 
    for (Int_t i = 0; i <= syshist->GetNbinsX(); ++i) {
        Double_t binc = 0;
        binc += stacksum->GetBinContent(i);
        syshist->SetBinContent(i, binc);
        // calculate uncertainty: lumi uncertainty
        Double_t binerr2 = binc * binc * lumierr * lumierr;
        Double_t topunc = 0; // uncertainty on top xsec
        Double_t bin_brerr2 = binc * binc * brerr * brerr;

        double topxsecErr2 = 2.2 * 2.2 + 11.6 * 11.6;

        double topRelUnc = TMath::Sqrt(topxsecErr2) / topxsec;
        //Functionality for multiple signal histograms
        topunc += varhists[signalHist]->GetBinContent(i) * topRelUnc;
        binerr2 += (topunc * topunc);
        if (absoluteXsec && !particleXsec) binerr2 += bin_brerr2;
        syshist->SetLineColor(1);
        syshist->SetBinError(i, TMath::Sqrt(binerr2));
    }

    syshist->SetFillStyle(3004);
    syshist->SetFillColor(kBlack);
    //leg->AddEntry( syshist, "Uncertainty", "f" );


    varhists[0]->SetMaximum(1.5 * varhists[0]->GetBinContent(varhists[0]->GetMaximumBin()));

    varhists[0]->SetMinimum(0);
    varhists[0]->GetYaxis()->SetTitle("events");
    varhists[0]->GetXaxis()->SetNoExponent(kTRUE);
    varhists[0]->Draw("e");

    //Add the binwidth to the yaxis in yield plots
    TString ytitle = TString(varhists[0]->GetYaxis()->GetTitle()).Copy();
    double binwidth = varhists[0]->GetXaxis()->GetBinWidth(1);
    std::ostringstream width;
    width << binwidth;
    if (name.Contains("Rapidity") || name.Contains("Eta")) { ytitle.Append(" / ").Append(width.str()); }
    else if (name.Contains("pT", TString::kIgnoreCase) || name.Contains("Mass", TString::kIgnoreCase) || name.Contains("MET") || name.Contains("HT")) { ytitle.Append(" / ").Append(width.str()).Append(" GeV"); };
    varhists[0]->GetYaxis()->SetTitle(ytitle);

    stack->Draw("same HIST");

    //Only necessary if we want error bands

    varhists[0]->Draw("same, e1");
    DrawCMSLabels(cmsPrelimLabelMode);
    DrawDecayChLabel(channelLabel[channelType]);
    setControlPlotLegendStyle(varhistsPlotting, legends, leg);
    leg->Draw("SAME");
    gPad->RedrawAxis();

    c1->Print(outdir.Copy() + "preunfolded_" + name + ".eps");
    c1->Print(outdir.Copy() + "preunfolded_" + name + ".pdf");
    TFile out_root(outdir.Copy() + "preunfolded_" + name + "_source.root", "RECREATE");

    varhists[0]->Write(name + "_data");
    stacksum->Write(name + "_allmc");
    c1->Write(name + "_canvas");
    out_root.Close();
    c1->Clear();
    delete c1;
    delete stacksum;
    for (unsigned int i = 0; i < hists.size(); i++) {
        delete varhists[i];
        //         delete varhistsPlotting.at(i);
    }
}



// get generator cross section curve for NLO prediction
TH1* Plotter13TeV::GetNloCurve(const char* particle, const char* quantity, const char* generator) {

    TH1::AddDirectory(kFALSE);
    TString histname;
    if (strcmp(particle, "TopQuarks") == 0 || strcmp(particle, "TtBar") == 0) {
        histname = "total_";
    }
    else {
        histname = "visible_";
    }
    histname.Append(particle);
    histname.Append(quantity);
    histname.Append("_");
    histname.Append(generator);


    TString filename;
    if (strcmp(generator, "Powheg") == 0) { filename = "selectionRoot/Nominal/emu/ttbarsignalplustau_powheg.root"; }
    else if (strcmp(generator, "MCatNLO") == 0) { filename = ttbar::DATA_PATH_DILEPTONIC() + "/MCatNLO_status3_v20120729.root"; }
    else if (strcmp(generator, "MCNLOup") == 0) { filename = ttbar::DATA_PATH_DILEPTONIC() + "/MCatNLO_Uncert_Up_status3_v20120729.root"; }
    else if (strcmp(generator, "MCNLOdown") == 0) { filename = ttbar::DATA_PATH_DILEPTONIC() + "/MCatNLO_Uncert_Down_status3_v20120729.root"; }

    TH1* hist = fileReader->GetClone<TH1>(filename, histname, true);
    if (hist) {
        TH1* weight = fileReader->GetClone<TH1>(filename, TString("total_LeptonsPt_").Append(generator).Data(), true);
        if (!weight) {
            std::cerr << "WARNING in GetNloCurve: histogram to extract original number of events could not be opened! No weighting applied!" << std::endl;
        }
        return hist;
    }
    std::cerr << "WARNING in GetNloCurve: input file could not been opened! Returning dummy!" << std::endl;
    hist = new TH1D();
    return hist; //I'd prefer to return nullptr
}

TH1* Plotter13TeV::GetNloCurve(TString NewName, TString Generator) {

    TString filename;
    if (Generator == "Nominal") {
        filename = "_ttbarsignalplustau.root";
    }
    else if (Generator == "MCATNLO") {
        filename = "_ttbarsignalplustau_powhegv2Herwig.root";
    }
    else if (Generator == "POWHEG") {
        filename = "_ttbarsignalplustau_fromDilepton_amcatnlofxfx.root";
    }
    else if (Generator == "POWHEGHERWIG") {
        filename = "_ttbarsignalplustau_madgraphmlm.root";
    }
    else if (Generator == "PERUGIA11") {
        filename = "_ttbarsignalplustau_Perugia11.root";
    }
    else if (Generator == "MATCH_UP") {
        filename = "_ttbarsignalplustau_matchup.root";
    }
    else if (Generator == "MATCH_DOWN") {
        filename = "_ttbarsignalplustau_matchdown.root";
    }
    else if (Generator == "SCALE_UP") {
        filename = "_ttbarsignalplustau_scaleup.root";
    }
    else if (Generator == "SCALE_DOWN") {
        filename = "_ttbarsignalplustau_scaledown.root";
    }
    else {
        std::cerr << "Unknown Generator!\n";
        std::exit(2);
    }

    const static std::vector<TString> channelName{ "ee", "mumu", "emu" };
    std::vector<TString> files;
    assert(channelType >= 0);
    assert(channelType <= 3);
    for (int i = 0; i <= 2; ++i) {
        if (channelType == i || channelType == 3) {
            files.push_back("selectionRoot/" + Generator + "/" + channelName.at(i) + "/" + channelName.at(i) + filename);
            std::cout << "Getting NLO curve from: " << files.at(files.size() - 1) << std::endl;
        }
    }

    TString histname(visTheoryPrefix + NewName);
    TH1* hist = fileReader->GetClone<TH1>(files.at(0), histname);
    for (size_t i = 1; i < files.size(); ++i) {
        hist->Add(fileReader->Get<TH1>(files.at(i), histname));
    }

    if (!revokeAntiQuantityCombination && !NewName.Contains("Lead")
        && (NewName.Contains("Lepton") || NewName.Contains("Top") || NewName.Contains("BJet")))
    {
        //loop over anti-particle histograms and add them
        TString antiName(visTheoryPrefix + "Anti" + NewName);
        for (const auto& file : files) {
            hist->Add(fileReader->Get<TH1>(file, antiName));
        }
    }

    if (absoluteXsec) hist->Scale(usefulTools13TeV->CalcLumiWeight(files.at(0)) / lumi);

    hist->SetName(Generator + "plot");

    return hist;
}

TH1* Plotter13TeV::RescaleNLOCurveAndPrepareBinnedCurve(TH1* inputhist, int bins, double Xbins[], double xsecnorm, double BrFrac) {

    double inputscale = 1. / inputhist->Integral("width");
    TString inputname = inputhist->GetName();
    if (!absoluteXsec) inputhist->Scale(inputscale);
    TH1* inputhistBinned = inputhist->Rebin(bins, inputname, Xbins);
    for (Int_t bin = 0; bin < bins; bin++) {
        if (!absoluteXsec) inputhistBinned->SetBinContent(bin + 1, inputhistBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin]) / inputhist->GetBinWidth(1)));
        else inputhistBinned->SetBinContent(bin + 1, inputhistBinned->GetBinContent(bin + 1) / ((Xbins[bin + 1] - Xbins[bin])));
        inputhistBinned->SetBinError(bin + 1, 1e-10);
    }
    if (!absoluteXsec) inputhistBinned->Scale(xsecnorm / inputhistBinned->Integral("width"));
    if (absoluteXsec && !particleXsec) inputhistBinned->Scale(1. / BrFrac);

    return inputhistBinned;
}

TH1* Plotter13TeV::GetNloCurveMass(TString NewName, TString Generator, TString Mass) {

    TString filename;
    if (Generator == "MASS_UP" && Mass == "173.5") {
        filename = "_ttbarsignalplustau_massup.root";
    }
    else if (Generator == "MASS_UP" && Mass == "175.5") {
        filename = "_ttbarsignalplustau_175_massup.root";
    }
    else if (Generator == "MASS_UP" && Mass == "178.5") {
        filename = "_ttbarsignalplustau_178_massup.root";
    }
    else if (Generator == "MASS_DOWN" && Mass == "171.5") {
        filename = "_ttbarsignalplustau_massdown.root";
    }
    else if (Generator == "MASS_DOWN" && Mass == "169.5") {
        filename = "_ttbarsignalplustau_169_massdown.root";
    }
    else if (Generator == "MASS_DOWN" && Mass == "166.5") {
        filename = "_ttbarsignalplustau_166_massdown.root";
    }
    else {
        std::cerr << "Unknown Generator!\n";
        std::exit(2);
    }

    const static std::vector<TString> channelName{ "ee", "mumu", "emu" };
    std::vector<TString> files;
    assert(channelType >= 0);
    assert(channelType <= 3);
    for (int i = 0; i <= 2; ++i) {
        if (channelType == i || channelType == 3) {
            files.push_back("selectionRoot/" + Generator + "/" + channelName.at(i) + "/" + channelName.at(i) + filename);
            cout << "Getting NLO curve from: " << files.at(files.size() - 1) << endl;
        }
    }

    TString histname(visTheoryPrefix + NewName);
    TH1* hist = fileReader->GetClone<TH1>(files.at(0), histname);
    for (size_t i = 1; i < files.size(); ++i) {
        hist->Add(fileReader->Get<TH1>(files.at(i), histname));
    }

    if (!revokeAntiQuantityCombination && !NewName.Contains("Lead")
        && (NewName.Contains("Lepton") || NewName.Contains("Top") || NewName.Contains("BJet")))
    {
        //loop over anti-particle histograms and add them
        TString antiName(visTheoryPrefix + "Anti" + NewName);
        for (const auto& file : files) {
            hist->Add(fileReader->Get<TH1>(file, antiName));
        }
    }
    return hist;
}


TH1F* Plotter13TeV::ConvertGraphToHisto(TGraphErrors* pGraph) {
    // takes data from a graph, determines binning and fills data into histogram
    Int_t NPoints = pGraph->GetN();
    Double_t BinLimits[NPoints + 1];
    // sort graph
    pGraph->Sort();
    // determine lower limit of histogram: half the distance to next point
    Double_t x0, x1, y;
    pGraph->GetPoint(0, x0, y);
    pGraph->GetPoint(1, x1, y);
    Double_t Distance = TMath::Abs(x0 - x1);
    BinLimits[0] = x0 - Distance / 2.;
    // now set upper limits for all the other points
    for (Int_t k = 0; k < NPoints - 1;k++) {
        pGraph->GetPoint(k, x0, y);
        pGraph->GetPoint(k + 1, x1, y);
        Distance = TMath::Abs(x0 - x1);
        BinLimits[k + 1] = x0 + Distance / 2.;
    }
    // for the last point set upper limit similar to first point:
    pGraph->GetPoint(NPoints - 2, x0, y);
    pGraph->GetPoint(NPoints - 1, x1, y);
    Distance = TMath::Abs(x0 - x1);
    BinLimits[NPoints] = x1 + Distance / 2.;
    // now we know the binning and can create the histogram:
    TString Name = "ConvertedHisto";
    // make name unique 
    Name += rand();
    TH1F* ThisHist = new TH1F(Name, "Converted Histogram", NPoints, BinLimits);
    // now fill the histogram
    for (Int_t i = 0; i < pGraph->GetN();i++) {
        Double_t x2, y2;
        pGraph->GetPoint(i, x2, y2);
        ThisHist->SetBinContent(i + 1, y2);
    }
    return ThisHist;
}

double Plotter13TeV::GetChi2(TGraphAsymmErrors* data, TH1* mc) {

    double chi2 = 0.0;

    for (int i = 0; i < (int)data->GetN(); ++i) {
        if (data->GetErrorYhigh(i) == 0 || data->GetErrorYlow(i) == 0) {
            std::cout << "When calculating the Chi2 the DATA TGraph has error 0 in bin " << i << std::endl;
            //exit(42);
            return 0;
        }
        double dataMinusMC = data->GetY()[i] - mc->GetBinContent(mc->FindBin(data->GetX()[i]));
        double dataError = (data->GetErrorYhigh(i) + data->GetErrorYlow(i)) / 2;
        chi2 += dataMinusMC * dataMinusMC / (dataError * dataError);
    }
    return chi2 / (data->GetN() - 1);
}


//TH1F* Plotter13TeV::reBinTH1FIrregularNewBinning(TH1F *histoOldBinning, const std::vector<double> &vecBinning, TString plotname, bool rescale=1){
TH1F* Plotter13TeV::reBinTH1FIrregularNewBinning(TH1F* histoOldBinning, TString plotname, bool rescale) {
    //  This function rebins a histogram using a variable binning
    // 
    //  (1) It is not required to have an equidistant binning.
    //  (2) Any type of ROOT-histgramme can be used, so the template 
    //      arguments should be 
    //      (a) histoT = TH1D,   TH1F,  ....
    //      (b) varT   = double, float, ....
    //  
    //  modified quantities: none
    //  used functions:      none
    //  used enumerators:    none
    //  
    //  "histoOldBinning":   plot to be re-binned
    //  "vecBinning":        vector containing all bin edges 
    //                       from xaxis.min to xaxis.max
    //  "rescale":           rescale the rebinned histogramme
    //                       (applicable for cross-section e.g.) 
    //  std::cout << std::endl;
    //  std::cout << std::endl;
    //  std::cout << "asdfasdfasdfasdfasdf hallo user! " << plotname << " " << rescale << std::endl;
    //  std::cout << "histoOldBinning = ";
    //  for ( int i = 0 ; i < histoOldBinning->GetXaxis()->GetNbins() + 1; i++ ) std::cout << " " << histoOldBinning->GetXaxis()->GetBinLowEdge(i+1);
    //  std::cout << std::endl;
    //  std::cout << std::endl;
    //  std::cout << std::endl;

    unsigned int vecIndex = 0;

    // fill vector into array to use appropriate constructor of TH1-classes
    const double* binArray = XAxisbins.data();

    // create histo with new binning
    TH1F* histoNewBinning = new TH1F("histoNewBinning" + plotname, "histoNewBinning" + plotname, XAxisbins.size() - 1, binArray);

    // fill contents of histoOldBinning into histoNewBinning and rescale if applicable
    for (vecIndex = 0; vecIndex < XAxisbins.size() - 1; vecIndex++) {

        double lowEdge = XAxisbins[vecIndex];
        if (plotname == "topPt" && vecIndex == 0 && lowEdge == 0.0) lowEdge += 10;  // adhoc fix to compensate for minimum top-Pt cut in NNLO curve
        double highEdge = XAxisbins[vecIndex + 1];
        double newBinWidth = highEdge - lowEdge;
        double newBinCenter = 0.5 * (highEdge + lowEdge);
        double binSum = 0.0;

        for (int j = 1; j < histoOldBinning->GetNbinsX(); j++) {

            double oldBin = histoOldBinning->GetBinCenter(j);

            if ((oldBin >= lowEdge) && (oldBin < highEdge)) {
                if (rescale) binSum += histoOldBinning->GetBinContent(j) * histoOldBinning->GetBinWidth(j);
                else         binSum += histoOldBinning->GetBinContent(j);
            }
        }

        if (rescale) histoNewBinning->Fill(newBinCenter, binSum / newBinWidth);
        else histoNewBinning->Fill(newBinCenter, binSum);
    }

    return (TH1F*)histoNewBinning->Clone();
}


void Plotter13TeV::setResultLegendStyle(TLegend* leg, const bool result)
{
    double x1 = 0.560, y1 = 0.655; //### orig
    //double x1 = 0.540, y1 = 0.685; //###
    double height = 0.175, width = 0.275;
    if (result) {
        if (name.Contains("Eta")) {
            x1 = 0.4;
            y1 = 0.39;
        }
        if (name.Contains("TopRapidity") || name.Contains("HypTTBarDeltaRapidity")) {
            x1 = 0.4;
            y1 = 0.41;
            height += 0.045;
        }
        if (name.Contains("DeltaPhi")) {
            x1 = 0.4;
            y1 = 0.56;
        }
        if (name.Contains("TTBarRapidity") || name.Contains("TTBarpT") || name.Contains("TTBarMass") || name.Contains("ToppT")) {
            height += 0.045;
            y1 -= 0.045;
        }
        if (name.Contains("TopRapidityAbs")) {
            x1 = 0.560, y1 = 0.655;
            height = 0.175, width = 0.275;
            height += 0.045;
            y1 -= 0.045;
        }
        if (name.Contains("JetMultpt")) { //########
            x1 = 0.540, y1 = 0.685;
            height = 0.175, width = 0.275;
            height += 0.045;
            y1 -= 0.045;
        }

    }
    leg->SetX1NDC(x1);
    leg->SetY1NDC(y1);
    leg->SetX2NDC(x1 + width);
    leg->SetY2NDC(y1 + height);

    leg->SetTextFont(42);
    leg->SetTextAlign(12);
    leg->SetTextSize(0.035);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
}


void Plotter13TeV::setControlPlotLegendStyle(std::vector< TH1* > drawhists, std::vector< TString > legends, TLegend* leg, TLegend* leg1, TLegend* leg2)
{
    //this appendix  to the legend aligns vertically all entries in the legend
    std::string  appendix = "                                                                                    ";

    //hardcoded ControlPlot legend
    std::vector<TString> OrderedLegends;
    TString pseudoData13TeVLeg = "Ps.-Data";
    OrderedLegends.push_back("Data");
    OrderedLegends.push_back("t#bar{t} signal"); //###
    OrderedLegends.push_back("t#bar{t} other");  //###  
    OrderedLegends.push_back("Single t");        //###
    OrderedLegends.push_back("W+jets");          //###
    //OrderedLegends.push_back("Z / #gamma* #rightarrow ee/#mu#mu/#tau#tau");
    OrderedLegends.push_back("Z+jets");          //###
    //OrderedLegends.push_back("Z / #gamma* #rightarrow #tau#tau");
    OrderedLegends.push_back("t#bar{t}+Z/W");
    OrderedLegends.push_back("Diboson");
    if (this->addQCDToControlPlot()) OrderedLegends.push_back("QCD multijet");

    leg->Clear();
    if (leg1) leg1->Clear();
    if (leg2) leg2->Clear();
    for (size_t i = 0; i < OrderedLegends.size(); ++i) {
        for (size_t j = 0; j < drawhists.size(); ++j) {
            if (OrderedLegends.at(i) == legends.at(j)) {
                if (OrderedLegends.at(i) == "Data") {
                    if (usePoissonSmearedPseudoData) {
                        leg->AddEntry(drawhists.at(j), pseudoData13TeVLeg + appendix, "pe");
                        if (leg1)leg1->AddEntry(drawhists.at(j), pseudoData13TeVLeg + appendix, "pe");
                    }
                    else {
                        leg->AddEntry(drawhists.at(j), OrderedLegends.at(i) + appendix, "pe");
                        if (leg1)leg1->AddEntry(drawhists.at(j), OrderedLegends.at(i) + appendix, "pe");
                    }
                    break;
                }
                else {
                    leg->AddEntry(drawhists.at(j), OrderedLegends.at(i) + appendix, "f");
                    if (leg1 && i < 4) leg1->AddEntry(drawhists.at(j), OrderedLegends.at(i) + appendix, "f");
                    if (leg2 && i >= 4) leg2->AddEntry(drawhists.at(j), OrderedLegends.at(i) + appendix, "f");
                    break;
                }
            }
        }
    }
    //coordinates for legend without splitting
    double x1 = 1.00 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.25 + 0.005;
    double y1 = 1.00 - gStyle->GetPadTopMargin() - gStyle->GetTickLength() - 0.05 - 0.03 * leg->GetNRows();
    double x2 = 1.00 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.025;
    double y2 = 1.00 - gStyle->GetPadTopMargin() - 0.8 * gStyle->GetTickLength();

    leg->SetX1NDC(x1);
    leg->SetY1NDC(y1);
    leg->SetX2NDC(x2);
    leg->SetY2NDC(y2);

    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextAlign(12);

    if (!leg1) return;
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.03);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextAlign(12);

    // Define shifts of legends
    double xShift = 0.8 * (x2 - x1);
    //double yShift=0.06; //### orig
    //double yShift=0.0; //###
    double yShift = 0.03; //###
    double y1mod = y1 + 0.45 * (y2 - y1);

    // coordinates for splitted legends
    leg1->SetX1NDC(x1 - xShift);
    leg1->SetY1NDC(y1mod - yShift);
    leg1->SetX2NDC(x2 - xShift);
    leg1->SetY2NDC(y2 - yShift);

    if (!leg2) return;
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.03);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetTextAlign(12);

    // coordinates for splitted legends
    leg2->SetX1NDC(x1 + 0.01);
    leg2->SetY1NDC(y1mod - yShift);
    leg2->SetX2NDC(x2 + 0.01);
    leg2->SetY2NDC(y2 - yShift);

}

void Plotter13TeV::DrawLabel(TString text, const double x1, const double y1, const double x2, const double y2, int centering, double textSize) {
    //Function to add Kidonakis references to DiffXSection plots' legends
    TPaveText* label = new TPaveText(x1, y1, x2, y2, "br NDC");
    label->AddText(text);
    label->SetFillStyle(0);
    label->SetBorderSize(0);

    if (textSize != 0) label->SetTextSize(textSize);
    label->SetTextAlign(centering);
    label->Draw("same");
}


// Draw label for Decay Channel in upper left corner of plot
void Plotter13TeV::DrawDecayChLabel(TString decaychannel, double textSize) {

    TPaveText* decch = new TPaveText();

    decch->AddText(decaychannel);

    //For old and semi-official CMS label
    decch->SetX1NDC(gStyle->GetPadLeftMargin() + gStyle->GetTickLength());
    decch->SetY1NDC(1.0 - gStyle->GetPadTopMargin() - gStyle->GetTickLength() - 0.05);
    decch->SetX2NDC(gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15);
    decch->SetY2NDC(1.0 - gStyle->GetPadTopMargin() - gStyle->GetTickLength());


    //For CMS official label: channel name at upper right corner 
    //decch->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.55   );
    //decch->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );
    //decch->SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.60 );
    //decch->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength()        );    

    decch->SetFillStyle(0);
    decch->SetBorderSize(0);
    if (textSize != 0) decch->SetTextSize(textSize);
    decch->SetTextAlign(12);
    decch->Draw("same");
}

// // Draw official labels (CMS Preliminary, luminosity and CM energy) above plot - OLD STYLE
// void Plotter13TeV::DrawCMSLabels(int cmsprelim, double energy, double textSize) {

//     const char *text;
//     if(cmsprelim ==2 ) {//Private work for PhDs students
//         text = "Private Work, %2.f pb^{-1} at #sqrt{s} = %2.f TeV";
//     } else if (cmsprelim==1) {//CMS preliminary label
//         text = "CMS Preliminary, %2.f pb^{-1} at #sqrt{s} = %2.f TeV";
//     } else {//CMS label
//         text = "CMS, %2.f pb^{-1} at #sqrt{s} = %2.f TeV";
//     }

//     TPaveText *label = new TPaveText();
//     label->SetX1NDC(gStyle->GetPadLeftMargin());
//     label->SetY1NDC(1.0-gStyle->GetPadTopMargin());
//     label->SetX2NDC(1.0-gStyle->GetPadRightMargin());
//     label->SetY2NDC(1.0);
//     label->SetTextFont(42);
//     label->AddText(Form(text, lumi, energy));
//     label->SetFillStyle(0);
//     label->SetBorderSize(0);
//     if (textSize!=0) label->SetTextSize(textSize);
//     label->SetTextAlign(32);
//     label->Draw("same");
// }

// Draw official labels (CMS Preliminary, luminosity and CM energy) above plot
void Plotter13TeV::DrawCMSLabels(int cmsprelim, double textSize) {

    const char* text = "%2.1f fb^{-1} (%2.f TeV)";
    //   const char *text = "%2.f pb^{-1} (%2.f TeV)";

    TPaveText* label = new TPaveText();
    label->SetX1NDC(gStyle->GetPadLeftMargin());
    //label->SetY1NDC(1.0-gStyle->GetPadTopMargin()); //## orig
    label->SetY1NDC(0.98 - gStyle->GetPadTopMargin());
    label->SetX2NDC(1.0 - gStyle->GetPadRightMargin());
    //label->SetY2NDC(1.0); //## orig
    label->SetY2NDC(0.98); //##
    label->SetTextFont(42);
    label->AddText(Form(text, lumi / 1000, energy)); //## orig
    //label->AddText(Form(text, lumi, energy));
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    if (textSize != 0) label->SetTextSize(textSize);
    label->SetTextAlign(32);
    label->Draw("same");

    TPaveText* cms = new TPaveText();
    cms->AddText("CMS");

    //Official
    //cms->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
    //cms->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );
    //cms->SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
    //cms->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength()        );

    //Semi-official
    cms->SetX1NDC(0.471 - gStyle->GetPadLeftMargin());
    cms->SetY1NDC(0.98 - gStyle->GetPadTopMargin());
    cms->SetX2NDC(1.0 - gStyle->GetPadRightMargin());
    cms->SetY2NDC(0.98);

    //std::cout << "############# " << gStyle->GetPadLeftMargin() << " " << gStyle->GetPadTopMargin() << " " <<  gStyle->GetPadRightMargin() << std::endl;

    cms->SetFillStyle(0);
    cms->SetBorderSize(0);
    if (textSize != 0) cms->SetTextSize(textSize * 1.1);
    cms->SetTextAlign(12);
    cms->SetTextFont(61);
    cms->Draw("same");

    if (cmsprelim > 0) {
        TPaveText* extra = new TPaveText();
        if (cmsprelim == 2) { extra->AddText("Private Work"); }
        else { extra->AddText("Preliminary"); }

        //Official
        //extra->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
        //extra->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.10 );
        //extra->SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
        //extra->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );

        //Semi-official
        extra->SetX1NDC(0.47 - gStyle->GetPadLeftMargin() + gStyle->GetTickLength());
        extra->SetY1NDC(0.975 - gStyle->GetPadTopMargin());
        extra->SetX2NDC(1.0 - gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 1.50);
        extra->SetY2NDC(0.975);

        extra->SetFillStyle(0);
        extra->SetBorderSize(0);
        if (textSize != 0) extra->SetTextSize(textSize);
        extra->SetTextAlign(12);
        extra->SetTextFont(52);
        extra->Draw("same");
    }
}


void Plotter13TeV::PrintResultTotxtFile(TString channel, double binCenters[], TGraphAsymmErrors* tga_DiffXSecPlot, TGraphAsymmErrors* tga_DiffXSecPlotwithSys)
{
    // Keeping this functionality for the case that further scripts rely on this output format (from 8 TeV)

    FILE* file;
    if (channel == "combined" && performUnfoldingInCombinedChannel) file = fopen("Plots/FinalResults/" + channel + "Unfolded/" + name + "LaTeX.txt", "w");
    else file = fopen("Plots/FinalResults/" + channel + "/" + name + "LaTeX.txt", "w");
    fprintf(file, "Variable: %s Channel: %s \n", name.Data(), channelLabel.at(channelType).Data());
    fprintf(file, "BinCenter & LowXbinEdge  &  HighXbinEdge  &   DiffXSec  &  StatError(\\%%)  & SystError(\\%%)  & TotalError(\\%%) \\\\ \n");
    fprintf(file, "\\hline \n");
    for (int i = 0; i < (int)tga_DiffXSecPlot->GetN(); i++) {
        double DiffXSec = tga_DiffXSecPlot->GetY()[i];
        double RelStatErr = 0, RelSysErr = 0, RelTotErr = 0;
        if (DiffXSec != 0.0) {
            RelStatErr = 100 * (tga_DiffXSecPlot->GetErrorY(i)) / DiffXSec;
            RelTotErr = 100 * (tga_DiffXSecPlotwithSys->GetErrorY(i)) / DiffXSec;
            if (RelTotErr >= RelStatErr) RelSysErr = TMath::Sqrt(RelTotErr * RelTotErr - RelStatErr * RelStatErr);
        }
        fprintf(file, "$%5.2f$ & $%5.2f$ to $%5.2f$ & %5.2e & %2.1f & %2.1f & %2.1f \\\\ \n", binCenters[i], XAxisbins.at(i), XAxisbins.at(i + 1), DiffXSec, RelStatErr, RelSysErr, RelTotErr);
    }
    fclose(file);
}

void Plotter13TeV::PrintResultPlusGenTotxtFile(TString channel, double binCenters[], TGraphAsymmErrors* tga_DiffXSecPlot, TGraphAsymmErrors* tga_DiffXSecPlotwithSys, TH1* mcPH, TH1* mcAMC, TH1* mcMLM, TH1* mcHpp)
{
    FILE* file;
    if (channel == "combined" && performUnfoldingInCombinedChannel) file = fopen("Plots/FinalResults/" + channel + "Unfolded/" + name + "WithGeneratorsLaTeX.txt", "w");
    else file = fopen("Plots/FinalResults/" + channel + "/" + name + "WithGeneratorsLaTeX.txt", "w");
    fprintf(file, "Variable: %s Channel: %s \n", name.Data(), channelLabel.at(channelType).Data());
    fprintf(file, "BinCenter & LowXbinEdge  &  HighXbinEdge  &   DiffXSec  &  StatError(\\%%)  & SystError(\\%%)  & TotalError(\\%%)& Pv2+P8 (\\%%) & aMCNLO+P8 (\\%%) & MGMLM+P8 (\\%%) & Pv2+H++(\\%%) \\\\ \n");
    fprintf(file, "\\hline \n");
    for (int i = 0; i < (int)tga_DiffXSecPlot->GetN(); i++) {
        double DiffXSec = tga_DiffXSecPlot->GetY()[i];
        double RelStatErr = 0, RelSysErr = 0, RelTotErr = 0;
        double vmcPH = mcPH->GetBinContent(mcPH->FindBin(tga_DiffXSecPlot->GetX()[i]));
        double vmcAMC = mcAMC->GetBinContent(mcAMC->FindBin(tga_DiffXSecPlot->GetX()[i]));
        double vmcMLM = mcMLM->GetBinContent(mcMLM->FindBin(tga_DiffXSecPlot->GetX()[i]));
        double vmcHpp = mcHpp->GetBinContent(mcHpp->FindBin(tga_DiffXSecPlot->GetX()[i]));
        if (DiffXSec != 0.0) {
            RelStatErr = 100 * (tga_DiffXSecPlot->GetErrorY(i)) / DiffXSec;
            RelTotErr = 100 * (tga_DiffXSecPlotwithSys->GetErrorY(i)) / DiffXSec;
            if (RelTotErr >= RelStatErr) RelSysErr = TMath::Sqrt(RelTotErr * RelTotErr - RelStatErr * RelStatErr);
        }
        fprintf(file, "$ %5.2f $ & $ %5.2f $ to $ %5.2f $ & %5.2e & %2.1f & %2.1f & %2.1f & %5.2e & %5.2e & %5.2e & %5.2e \\\\ \n", binCenters[i], XAxisbins.at(i), XAxisbins.at(i + 1), DiffXSec, RelStatErr, RelSysErr, RelTotErr, vmcPH, vmcAMC, vmcMLM, vmcHpp);
    }
    fclose(file);
}

void Plotter13TeV::CalcTotCovMtrx(TString Channel, TString Variable, std::vector<TString>vec_systematic) {

    //This function only works for standalone unfoldings (ee, emu, mumu, combined-unfolded)
    if (Channel == "combined" && !performUnfoldingInCombinedChannel) return;
    if (Channel == "combined" && performUnfoldingInCombinedChannel) Channel += "Unfolded";

    int debug = 0; // toggle on only if silent == 1
    int silent = 1;
    int mtrxplots = 0;

    ifstream NominalFile;
    string nominalfilename = string("UnfoldingResults/Nominal/" + Channel + "/" + Variable + "Results.txt");
    NominalFile.open(nominalfilename);
    if (!NominalFile.is_open()) {
        std::cout << "Nominal: " << nominalfilename << std::endl;
        std::cout << "CalcTotCovMtrx(): the input file cannot be opened. Exit!" << std::endl;
        exit(478);
    }

    std::vector<double> BinCenters, BinLowEdge, BinHigEdge, BinWidth;
    std::vector<double> CentralValueNom, StatErrorNom;
    TString Dummy = "";
    double dummy = 0;

    TString appendix = "";
    if (absoluteXsec) appendix += "absolute";
    else appendix += "normalized";
    if (particleXsec) appendix += ", particle";
    else appendix += ", parton";

    if (!silent) {
        std::cout << "Results Table for measurement: " << Variable << " in " << Channel << " channel, " << appendix << std::endl;
        std::cout << "Bin";
        for (Int_t bin = 0; bin < bins; bin++) std::cout << " " << bin + 1;
        std::cout << std::endl;
        std::cout << "Nominal central value";
    }

    for (Int_t bin = 0; bin < bins; bin++) {
        double XAxisbinCenters_Nom = 0, XLowEdge_Nom = 0, XHigEdge_Nom = 0;
        double CentralValue_Nom = 0, StatError_Nom = 0;
        NominalFile >> Dummy >> XAxisbinCenters_Nom >> Dummy >> XLowEdge_Nom >> Dummy >> XHigEdge_Nom >> Dummy >> CentralValue_Nom >> Dummy >> StatError_Nom >> Dummy >> dummy;

        BinCenters.push_back(XAxisbinCenters_Nom);
        BinLowEdge.push_back(XLowEdge_Nom);
        BinHigEdge.push_back(XHigEdge_Nom);
        BinWidth.push_back(XHigEdge_Nom - XLowEdge_Nom);
        CentralValueNom.push_back(CentralValue_Nom);
        StatErrorNom.push_back(StatError_Nom);
        if (!silent) std::cout << " " << CentralValue_Nom;
        if (debug) std::cout << Variable << " Bin: " << bin << " " << XAxisbinCenters_Nom << " " << XLowEdge_Nom << " " << XHigEdge_Nom << " " << CentralValue_Nom << " " << StatError_Nom << std::endl;
    }
    if (!silent) std::cout << std::endl;
    NominalFile.close();

    if (!silent) {
        std::cout << "Stat. unc. value";
        for (Int_t bin = 0; bin < bins; bin++) std::cout << " " << StatErrorNom.at(bin) / CentralValueNom.at(bin);
        std::cout << std::endl;
    }

    std::vector<TString> suffix;
    suffix.push_back("UP"); suffix.push_back("DOWN");

    std::vector<TString> fullVectorSystematic;
    for (int syst = 0; syst < (int)vec_systematic.size(); syst++) {
        for (int sfx = 0; sfx < (int)suffix.size(); sfx++) {
            TString systName = vec_systematic.at(syst);
            if (vec_systematic.at(syst) == "BFRAG_PETERSON" || vec_systematic.at(syst) == "ERDON" || vec_systematic.at(syst) == "ERDONRETUNE" || vec_systematic.at(syst) == "GLUONMOVETUNE" ||
                vec_systematic.at(syst) == "POWHEG" || vec_systematic.at(syst) == "MCATNLO") {
                if (suffix.at(sfx) == "DOWN") continue;
            }
            else systName += suffix.at(sfx);
            fullVectorSystematic.push_back(systName);
        }
    }

    if (absoluteXsec) {
        fullVectorSystematic.push_back("LUMI_UP");
        fullVectorSystematic.push_back("LUMI_DOWN");
        if (!particleXsec) {
            fullVectorSystematic.push_back("DILEPBR_UP");
            fullVectorSystematic.push_back("DILEPBR_DOWN");
        }
    }

    std::vector<std::vector<double>> systMatrix((int)fullVectorSystematic.size(), vector<double>(bins));

    for (int syst = 0; syst < (int)fullVectorSystematic.size(); syst++) {
        TString systName = fullVectorSystematic.at(syst);
        ifstream SystFile;
        if (!systName.Contains("LUMI") && !systName.Contains("DILEPBR")) {
            string systfilename = string("UnfoldingResults/" + systName + "/" + Channel + "/" + Variable + "Results.txt");
            SystFile.open(systfilename);
            if (!SystFile.is_open()) {
                std::cout << "Syst: " << systfilename << std::endl;
                std::cout << "CalcTotCovMtrx(): the input file cannot be opened. Exit!" << std::endl;
                exit(479);
            }
        }

        if (systName == "POWHEG") systName = "MOD";
        if (systName == "MCATNLO") systName = "HAD";

        if (!silent) std::cout << systName;
        for (Int_t bin = 0; bin < bins; bin++) {
            double CentralValue_Sys = 0;
            double SystVariationPerBin = 0;
            if (!systName.Contains("LUMI") && !systName.Contains("DILEPBR")) {
                SystFile >> Dummy >> dummy >> Dummy >> dummy >> Dummy >> dummy >> Dummy >> CentralValue_Sys >> Dummy >> dummy >> Dummy >> dummy;
                SystVariationPerBin = (CentralValue_Sys - CentralValueNom.at(bin)) / CentralValueNom.at(bin);
            }
            if (systName == "LUMI_UP") SystVariationPerBin = lumiError;
            if (systName == "LUMI_DOWN") SystVariationPerBin = (-1.) * lumiError;
            if (systName == "DILEPBR_UP") SystVariationPerBin = BRError;
            if (systName == "DILEPBR_DOWN") SystVariationPerBin = (-1.) * BRError;
            if (systName.Contains("MASS_")) SystVariationPerBin /= massUncReductionFactor;
            if (systName.Contains("PSFSRSCALE_")) SystVariationPerBin /= scaleUncReductionFactor;
            systMatrix.at(syst).at(bin) = SystVariationPerBin;
            if (!silent) std::cout << " " << SystVariationPerBin;
        }
        if (!silent) std::cout << std::endl;
        if (!systName.Contains("LUMI") && !systName.Contains("DILEPBR")) SystFile.close();
    }

    std::vector<std::vector<double>> systCovMatrix(bins, vector<double>(bins));
    if (!silent) std::cout << "Syst. covariance matrix for measurement: " << Variable << " in " << Channel << " channel, " << appendix << std::endl;

    for (int syst = 0; syst < (int)fullVectorSystematic.size(); syst++) {
        if ((int)fullVectorSystematic.size() - syst == 1 && !silent) {
            std::cout << "Bin";
            for (Int_t bin = 0; bin < bins; bin++) std::cout << " " << bin + 1;
            std::cout << std::endl;
        }

        TString systName = fullVectorSystematic.at(syst);
        Int_t variationFactor = 1; // factor depending on a number of variations considered for the envelope construction (see implementation of: Plotter13TeV::GetEnvelopeForUnfoldedResults())
        // or reflecting the factor with which the variation is considered ( 1 - for standalone (e.g MCATNLO), 2 - for up/down (e.g. JER))
// most usual case for UP/DOWN variation
        if (systName.Contains("_UP") || systName.Contains("_DOWN")) variationFactor = 2;
        // now consider envelope cases, if used
        if (mergeScaleVariations_ && (systName.Contains("MESCALE_") || systName.Contains("MEFACSCALE_") || systName.Contains("MERENSCALE_") ||
            systName.Contains("PSISRSCALE_") || systName.Contains("PSFSRSCALE_"))) variationFactor = 10;
        else if (mergeBTagVariations_ && (systName.Contains("BTAG_"))) variationFactor = 6;
        else if (mergeTriggerVariations_ && (systName.Contains("TRIG_"))) variationFactor = 4;
        else if (mergeBFragVariations_ && (systName.Contains("BFRAG_"))) variationFactor = 3;
        else if (mergeColorRecVariations_ && (systName.Contains("ERDON") || systName.Contains("ERDONRETUNE") || systName.Contains("GLUONMOVETUNE"))) variationFactor = 3;
        if (debug) std::cout << systName << " " << variationFactor << std::endl;

        for (Int_t ibin = 0; ibin < bins; ibin++) {
            if ((int)fullVectorSystematic.size() - syst == 1 && !silent) std::cout << ibin + 1;
            for (Int_t jbin = 0; jbin < bins; jbin++) {
                systCovMatrix.at(ibin).at(jbin) += systMatrix.at(syst).at(ibin) * systMatrix.at(syst).at(jbin) / variationFactor;
                if ((int)fullVectorSystematic.size() - syst == 1 && !silent) std::cout << " " << systCovMatrix.at(ibin).at(jbin);
            }
            if ((int)fullVectorSystematic.size() - syst == 1 && !silent) std::cout << std::endl;
        }
    }

    std::vector<std::vector<double>> systCorrMatrix(bins, vector<double>(bins));
    if (!silent) {
        std::cout << "Syst. correlation matrix for measurement: " << Variable << " in " << Channel << " channel, " << appendix << std::endl;
        std::cout << "Bin";
        for (Int_t bin = 0; bin < bins; bin++) std::cout << " " << bin + 1;
        std::cout << std::endl;
    }
    for (Int_t ibin = 0; ibin < bins; ibin++) {
        if (!silent) std::cout << ibin + 1;
        for (Int_t jbin = 0; jbin < bins; jbin++) {
            systCorrMatrix.at(ibin).at(jbin) = systCovMatrix.at(ibin).at(jbin) / std::sqrt(systCovMatrix.at(ibin).at(ibin) * systCovMatrix.at(jbin).at(jbin));
            if (!silent) std::cout << " " << systCorrMatrix.at(ibin).at(jbin);
        }
        if (!silent) std::cout << std::endl;
    }

    // Read-out unfolding (stat.) covariance matrix, and construct the one that is in compatible units with "syst. covariance"
    TString theParticleName = "";
    if (name.Contains("Lepton")) {
        theParticleName = "Leptons";
        if (revokeAntiQuantityCombination && !name.Contains("Anti")) theParticleName = "Lepton";
        if (revokeAntiQuantityCombination && name.Contains("Anti")) theParticleName = "AntiLepton";
    }
    if (name.Contains("LLBar")) theParticleName = "LepPair";
    if (name.Contains("Top")) {
        theParticleName = "TopQuarks";
        if (revokeAntiQuantityCombination && !name.Contains("Anti")) theParticleName = "TopQuark";
        if (revokeAntiQuantityCombination && name.Contains("Anti")) theParticleName = "TopAntiQuark";
    }
    if (name.Contains("TTBar")) theParticleName = "TtBar";
    if (name.Contains("BBBar")) theParticleName = "BBbar";
    if (name.Contains("BJet")) {
        theParticleName = "BJets";
        if (revokeAntiQuantityCombination && !name.Contains("Anti")) theParticleName = "BJet";
        if (revokeAntiQuantityCombination && name.Contains("Anti")) theParticleName = "AntiBJet";
    }
    if (name.Contains("JetMult")) theParticleName = "Jets";
    if (name.Contains("ExtraJet")) theParticleName = "ExtraJets";
    if (name.Contains("DeltaRJet12")) theParticleName = "DeltaR";

    TString theQuantityName = "";
    if (name.Contains("pT")) theQuantityName = "Pt";
    if (name.Contains("Eta")) theQuantityName = "Eta";
    if (name.Contains("Rapidity")) theQuantityName = "Rapidity";
    if (name.Contains("Mass")) theQuantityName = "Mass";
    if (name.Contains("JetMult")) theQuantityName = "Mult";
    if (name.Contains("ExtraJetEta")) theQuantityName = "Eta";
    if (name.Contains("ExtraJetpT")) theQuantityName = "Pt";
    if (name.Contains("DeltaRJet12")) theQuantityName = "DeltaR";

    TString theChannelName = Channel;
    if (theChannelName.Contains("combined")) theChannelName = "combined";

    TString svdFileName = "SVD/Unfolding_" + theChannelName + "_" + theParticleName + "_" + theQuantityName + "_" + Variable + ".root";
    std::ifstream svdFile(svdFileName);
    if (!svdFile.is_open()) {
        std::cout << "SVD File '" << svdFileName << "'\n cannot be opened" << std::endl;
        std::cout << "Exiting!!" << std::endl;
        exit(12);
    }
    svdFile.close();

    TFile svdTFile(svdFileName);
    TString svdHistoName = "";
    if (absoluteXsec) svdHistoName = "SVD_" + theChannelName + "_" + theParticleName + "_" + theQuantityName + "_" + Variable + "_STATCOV";
    else svdHistoName = "SVD_" + theChannelName + "_" + theParticleName + "_" + theQuantityName + "_" + Variable + "_STATCOVNORM";
    TH2D* statCovHisto = static_cast<TH2D*>(svdTFile.Get(svdHistoName));
    const int svdBins = statCovHisto->GetNbinsX() - 2;
    if (debug) std::cout << "SVD bins: " << svdBins << std::endl;
    if (bins != svdBins) {
        std::cout << "Binning in the distribution of input variable " << Variable << " and the propagated SVD file don't match!" << std::endl;
        std::cout << "Exiting!!" << std::endl;
        exit(441);
    }

    std::vector<std::vector<double>> statCovMatrix(bins, vector<double>(bins));
    for (Int_t ibin = 0; ibin < bins; ibin++) {
        for (Int_t jbin = 0; jbin < bins; jbin++) {
            statCovMatrix.at(ibin).at(jbin) = statCovHisto->GetBinContent(ibin + 2, jbin + 2) / (BinWidth.at(ibin) * BinWidth.at(jbin) * CentralValueNom.at(ibin) * CentralValueNom.at(jbin));
            if (absoluteXsec) {
                statCovMatrix.at(ibin).at(jbin) /= (lumi * lumi);
                if (!particleXsec) statCovMatrix.at(ibin).at(jbin) /= (BranchingFraction[channelType] * BranchingFraction[channelType]);
            }
        }
        if (debug) std::cout << "Check normalization values for stat. cov. mtrx.: " << ibin << " " << BinWidth.at(ibin) << " " << CentralValueNom.at(ibin) << std::endl;
    }
    if (svdTFile.IsOpen()) svdTFile.Close();

    // Evaluate total covariance and correlation matrices + stat. correlation matrix

    std::vector<std::vector<double>> totCovMatrix(bins, vector<double>(bins));
    for (Int_t ibin = 0; ibin < bins; ibin++) {
        for (Int_t jbin = 0; jbin < bins; jbin++) {
            totCovMatrix.at(ibin).at(jbin) = statCovMatrix.at(ibin).at(jbin) + systCovMatrix.at(ibin).at(jbin);
        }
    }
    std::vector<std::vector<double>> totCorrMatrix(bins, vector<double>(bins));
    std::vector<std::vector<double>> statCorrMatrix(bins, vector<double>(bins));
    for (Int_t ibin = 0; ibin < bins; ibin++) {
        for (Int_t jbin = 0; jbin < bins; jbin++) {
            totCorrMatrix.at(ibin).at(jbin) = totCovMatrix.at(ibin).at(jbin) / std::sqrt(totCovMatrix.at(ibin).at(ibin) * totCovMatrix.at(jbin).at(jbin));
            statCorrMatrix.at(ibin).at(jbin) = statCovMatrix.at(ibin).at(jbin) / std::sqrt(statCovMatrix.at(ibin).at(ibin) * statCovMatrix.at(jbin).at(jbin));
        }
    }

    // Plots the matrices via TMtrxD functionality

    if (mtrxplots) {
        TMatrixD systCovTMtrxD(bins, bins);
        TMatrixD systCorrTMtrxD(bins, bins);
        TMatrixD statCovTMtrxD(bins, bins);
        TMatrixD totCovTMtrxD(bins, bins);
        TMatrixD totCorrTMtrxD(bins, bins);
        for (Int_t ibin = 0; ibin < bins; ibin++) {
            for (Int_t jbin = 0; jbin < bins; jbin++) {
                systCovTMtrxD[ibin][jbin] = systCovMatrix.at(ibin).at(jbin);
                systCorrTMtrxD[ibin][jbin] = systCorrMatrix.at(ibin).at(jbin);
                statCovTMtrxD[ibin][jbin] = statCovMatrix.at(ibin).at(jbin);
                totCovTMtrxD[ibin][jbin] = totCovMatrix.at(ibin).at(jbin);
                totCorrTMtrxD[ibin][jbin] = totCorrMatrix.at(ibin).at(jbin);
            }
        }
        if (debug) {
            systCovTMtrxD.Print();
            std::cout << "Syst. covariance matrix determinant : " << systCovTMtrxD.Determinant() << std::endl;
            statCovTMtrxD.Print();
            std::cout << "Stat. covariance matrix determinant : " << statCovTMtrxD.Determinant() << std::endl;
            totCovTMtrxD.Print();
            std::cout << "Total covariance matrix determinant : " << totCovTMtrxD.Determinant() << std::endl;
        }

        TString outdir = ttbar::assignFolder(outpathPlots, TString(Channel + "/UnfMtrxAndTabs"), TString("FinalResults"));
        TCanvas* c = new TCanvas("Matrix", "Matrix", 1200, 1200);
        systCovTMtrxD.Draw("colz,text");
        c->Print(outdir.Copy() + Variable + "_systCovMtrx.pdf");
        c->Clear();
        systCorrTMtrxD.Draw("colz,text");
        c->Print(outdir.Copy() + Variable + "_systCorrMtrx.pdf");
        c->Clear();
        statCovTMtrxD.Draw("colz,text");
        c->Print(outdir.Copy() + Variable + "_statCovMtrx.pdf");
        c->Clear();
        totCovTMtrxD.Draw("colz,text");
        c->Print(outdir.Copy() + Variable + "_totCovMtrx.pdf");
        c->Clear();
        totCorrTMtrxD.Draw("colz,text");
        c->Print(outdir.Copy() + Variable + "_totCorrMtrx.pdf");
        c->Clear();
        delete c;
    }

    //Printing the matrices to .txt-files for analysis and as the latex-table

    //ofstream resTableFile;
    ofstream statCovMtrxFile;
    ofstream systCovMtrxFile;
    ofstream totCovMtrxFile;
    ofstream totCorrMtrxFile;

    //ofstream resTableLatex;
    ofstream statCorrMtrxLatex;
    ofstream totCorrMtrxLatex;

    statCovMtrxFile.open(ttbar::assignFolder(outpathPlots, TString(Channel + "/UnfMtrxAndTabs"), TString("FinalResults")) + Variable + "_statCovMtrxFile.txt");
    systCovMtrxFile.open(ttbar::assignFolder(outpathPlots, TString(Channel + "/UnfMtrxAndTabs"), TString("FinalResults")) + Variable + "_systCovMtrxFile.txt");
    totCovMtrxFile.open(ttbar::assignFolder(outpathPlots, TString(Channel + "/UnfMtrxAndTabs"), TString("FinalResults")) + Variable + "_totCovMtrxFile.txt");
    totCorrMtrxFile.open(ttbar::assignFolder(outpathPlots, TString(Channel + "/UnfMtrxAndTabs"), TString("FinalResults")) + Variable + "_totCorrMtrxFile.txt");

    statCorrMtrxLatex.open(ttbar::assignFolder(outpathPlots, TString(Channel + "/UnfMtrxAndTabs"), TString("FinalResults")) + Variable + "_statCorrMtrxLatexForMykola.txt");
    totCorrMtrxLatex.open(ttbar::assignFolder(outpathPlots, TString(Channel + "/UnfMtrxAndTabs"), TString("FinalResults")) + Variable + "_totCorrMtrxLatexForMykola.txt");

    if (!statCovMtrxFile.is_open() || !systCovMtrxFile.is_open() || !totCovMtrxFile.is_open() || !totCorrMtrxFile.is_open() ||
        !statCorrMtrxLatex.is_open() || !totCorrMtrxLatex.is_open()) {
        std::cout << "The output file cannot be opened. Exiting!!" << std::endl;
        exit(435);
    }

    fillUnfMtrxAndTabsToFile(statCovMtrxFile, statCovMatrix, bins, "Stat. covariance matrix", Variable, Channel, 0);
    fillUnfMtrxAndTabsToFile(systCovMtrxFile, systCovMatrix, bins, "Syst. covariance matrix", Variable, Channel, 0);
    fillUnfMtrxAndTabsToFile(totCovMtrxFile, totCovMatrix, bins, "Total covariance matrix", Variable, Channel, 0);
    fillUnfMtrxAndTabsToFile(totCorrMtrxFile, totCorrMatrix, bins, "Total correlation matrix", Variable, Channel, 0);

    TString latexApp = "";
    if (absoluteXsec) latexApp += "Abs";
    else latexApp += "Norm";
    if (particleXsec) latexApp += "Fid";
    else latexApp += "Full";

    fillUnfMtrxAndTabsToFile(statCorrMtrxLatex, statCorrMatrix, bins, "Stat. correlation matrix", Variable, TString((TString)("statCorrMtrx") + "-" + latexApp + "-" + Channel), 1);
    fillUnfMtrxAndTabsToFile(totCorrMtrxLatex, totCorrMatrix, bins, "Total correlation matrix", Variable, TString((TString)("totCorrMtrx") + "-" + latexApp + "-" + Channel), 1);

    statCovMtrxFile.close();
    systCovMtrxFile.close();
    totCovMtrxFile.close();
    totCorrMtrxFile.close();

    statCorrMtrxLatex.close();
    totCorrMtrxLatex.close();

}

void Plotter13TeV::fillUnfMtrxAndTabsToFile(std::ofstream& ifile, const std::vector<std::vector<double>>& matrix, int bins, TString Type, TString Variable, TString Channel, bool latex) {

    // purely for use in CalcTotCovMtrx()
    // Latex formatting is intended for use by Mykola for his PhD work (don't submit to git the reconfigured latex commands and other new commands for own use)

    TString caption, label, tabular = "";

    if (latex) {
        caption += "\\caption{" + Type + " " + Channel + " " + Variable + ".}";
        label += "\\label{table:" + Channel + "-" + Variable + "}";
        tabular += "\\begin{tabular}{c |";
        for (Int_t bin = 0; bin < bins; bin++) tabular += " c";
        tabular += "}";

        ifile << "\\begin{table}" << std::endl;
        ifile << "\\scriptsize" << std::endl;
        ifile << "\\aboverulesep=0ex" << std::endl;
        ifile << "\\belowrulesep=0ex" << std::endl;
        ifile << "\\renewcommand{\\arraystretch}{1.1}" << std::endl;
        ifile << caption << std::endl;
        ifile << "\\centering" << std::endl;
        ifile << label << std::endl;
        ifile << tabular << std::endl;
        ifile << "\\toprule" << std::endl;
    }
    else ifile << Type << " for " << Variable << " in " << Channel << " channel" << std::endl;
    ifile << "Bin";
    for (Int_t bin = 0; bin < bins; bin++) {
        if (latex) ifile << " & " << bin + 1;
        else ifile << " " << bin + 1;
    }
    if (latex) ifile << " \\\\";
    ifile << std::endl;
    if (latex) ifile << "\\midrule" << std::endl;

    for (Int_t ibin = 0; ibin < bins; ibin++) {
        ifile << ibin + 1;
        for (Int_t jbin = 0; jbin < bins; jbin++) {
            if (latex) {
                if (jbin < ibin) ifile << " & " << " ";
                else if (matrix.at(ibin).at(jbin) > 0.) ifile << " & +" << std::fixed << std::setprecision(1) << matrix.at(ibin).at(jbin) * 100.;
                else ifile << " & " << std::fixed << std::setprecision(1) << matrix.at(ibin).at(jbin) * 100.;
            }
            else ifile << " " << matrix.at(ibin).at(jbin);
        }
        if (latex) ifile << " \\\\";
        ifile << std::endl;
    }

    if (latex) {
        ifile << "\\bottomrule" << std::endl;
        ifile << "\\end{tabular}" << std::endl;
        ifile << "\\end{table}" << std::endl;
    }
}


void Plotter13TeV::CalcUpDownDifference(TString Channel, TString Syst_Up, TString Syst_Down, TString Variable) {

    ///Function to get the error of a certain systematic: Sqrt(0.5*(Err_Up*Err_Up + Err_Down*Err_Down))/Nominal

    ifstream NominalFile, UpFile, DownFile;
    string nominalfilename = string("UnfoldingResults/Nominal/" + Channel + "/" + Variable + "Results.txt");
    string upfilename = string("UnfoldingResults/" + Syst_Up + "/" + Channel + "/" + Variable + "Results.txt");
    string downfilename = string("UnfoldingResults/" + Syst_Down + "/" + Channel + "/" + Variable + "Results.txt");

    if (Channel == "combined" && performUnfoldingInCombinedChannel) {
        nominalfilename = string("UnfoldingResults/Nominal/" + Channel + "Unfolded/" + Variable + "Results.txt");
        upfilename = string("UnfoldingResults/" + Syst_Up + "/" + Channel + "Unfolded/" + Variable + "Results.txt");
        downfilename = string("UnfoldingResults/" + Syst_Down + "/" + Channel + "Unfolded/" + Variable + "Results.txt");
    }

    NominalFile.open(nominalfilename);
    UpFile.open(upfilename);
    DownFile.open(downfilename);

    std::vector<double> BinCenters, BinLowEdge, BinHigEdge;
    std::vector<double> RelativeError;
    TString Dummy = "";
    double dummy = 0;

    if (!NominalFile.is_open() || !UpFile.is_open() || !DownFile.is_open()) {
        std::cout << "Nominal: " << nominalfilename << std::endl;
        std::cout << "Sys Up : " << upfilename << std::endl;
        std::cout << "Sys Dow: " << downfilename << std::endl;
        std::cout << "The input file cannot be opened. Exiting!!" << std::endl;
        exit(433);
    }
    for (Int_t bin = 0; bin < bins; bin++) {

        double XAxisbinCenters_Nom = 0, XLowEdge_Nom = 0, XHigEdge_Nom = 0;
        double XAxisbinCenters_Up = 0, XLowEdge_Up = 0, XHigEdge_Up = 0;
        double XAxisbinCenters_Down = 0, XLowEdge_Down = 0, XHigEdge_Down = 0;
        double CentralValue_Nom = 0, CentralValue_Up = 0, CentralValue_Down = 0;

        NominalFile >> Dummy >> XAxisbinCenters_Nom >> Dummy >> XLowEdge_Nom >> Dummy >> XHigEdge_Nom >> Dummy >> CentralValue_Nom >> Dummy >> dummy >> Dummy >> dummy;
        UpFile >> Dummy >> XAxisbinCenters_Up >> Dummy >> XLowEdge_Up >> Dummy >> XHigEdge_Up >> Dummy >> CentralValue_Up >> Dummy >> dummy >> Dummy >> dummy;
        DownFile >> Dummy >> XAxisbinCenters_Down >> Dummy >> XLowEdge_Down >> Dummy >> XHigEdge_Down >> Dummy >> CentralValue_Down >> Dummy >> dummy >> Dummy >> dummy;

        BinCenters.push_back(XAxisbinCenters_Up);
        BinLowEdge.push_back(XLowEdge_Up);
        BinHigEdge.push_back(XHigEdge_Up);

        if (CentralValue_Nom != 0) {
            double up = std::fabs(CentralValue_Up - CentralValue_Nom);
            double down = std::fabs(CentralValue_Down - CentralValue_Nom);

            if (Syst_Up.Contains("POWHEG") || Syst_Down.Contains("POWHEG"))
            {
                up = std::fabs(CentralValue_Up - CentralValue_Down);
                down = up;
            }
            if (Syst_Up.Contains("MCATNLO") || Syst_Down.Contains("MCATNLO"))
            {
                up = std::fabs(CentralValue_Up - CentralValue_Down);
                down = up;
            }
            double rel_err = 0.5 * (up + down) / CentralValue_Nom;
            if (Syst_Up.Contains("MASS") && Syst_Down.Contains("MASS")) rel_err = rel_err / massUncReductionFactor;
            if (Syst_Up.Contains("PSFSRSCALE") && Syst_Down.Contains("PSFSRSCALE")) rel_err = rel_err / scaleUncReductionFactor;
            RelativeError.push_back(rel_err);
        }
    }
    NominalFile.close(); UpFile.close(); DownFile.close();

    if (Syst_Up == "POWHEG" && Syst_Down == "Nominal")
    {
        Syst_Up = "MOD_";
    }
    else if (Syst_Up == "MCATNLO" && Syst_Down == "Nominal")
    {
        Syst_Up = "HAD_";
    }
    else if (Syst_Up == "BFRAG_PETERSON" && Syst_Down == "Nominal")
    {
        Syst_Up = "BFRAG_PETERSON_";
    }
    else if (Syst_Up == "ERDON" && Syst_Down == "Nominal")
    {
        Syst_Up = "ERDON_";
    }
    else if (Syst_Up == "ERDONRETUNE" && Syst_Down == "Nominal")
    {
        Syst_Up = "ERDONRETUNE_";
    }
    else if (Syst_Up == "GLUONMOVETUNE" && Syst_Down == "Nominal")
    {
        Syst_Up = "GLUONMOVETUNE_";
    }
    else {
        Syst_Up.Remove(Syst_Up.Length() - 2, 2);
    }

    ofstream SystematicRelError;
    if (Channel == "combined" && performUnfoldingInCombinedChannel) SystematicRelError.open(ttbar::assignFolder("UnfoldingResults", "combinedUnfolded", Syst_Up) + Variable + "Results.txt");
    else SystematicRelError.open(ttbar::assignFolder("UnfoldingResults", Channel, Syst_Up) + Variable + "Results.txt");

    if (!SystematicRelError.is_open()) {
        std::cout << "The output file cannot be opened. Exiting!!" << std::endl;
        exit(434);
    }
    for (int bin = 0; bin < (int)RelativeError.size(); bin++) {
        SystematicRelError << "XAxisbinCenters[bin]: " << BinCenters.at(bin) << " bin: " << BinLowEdge.at(bin) << " to " << BinHigEdge.at(bin) << " SystematicRelError: " << RelativeError.at(bin) << std::endl;
    }
    SystematicRelError.close();

}


void Plotter13TeV::GetEnvelopeForUnfoldedResults(TString Channel, TString SystName, TString Variable) {

    if (SystName != "TOT_SCALE_" && SystName != "TOT_BTAG_" && SystName != "TOT_BTAG_LJET_" && SystName != "TOT_TRIG_" && SystName != "TOT_BFRAG_" && SystName != "TOT_COLORREC_") {
        std::cout << "Systematics name: " << SystName << std::endl;
        std::cout << "The calculation of envelope is not supported. Exiting!!" << std::endl;
        exit(433);
    }

    std::vector<TString> syst;

    if (SystName == "TOT_SCALE_") syst = { "MESCALE_UP", "MESCALE_DOWN", "MEFACSCALE_UP", "MEFACSCALE_DOWN", "MERENSCALE_UP", "MERENSCALE_DOWN",
    "PSISRSCALE_UP", "PSISRSCALE_DOWN", "PSFSRSCALE_UP", "PSFSRSCALE_DOWN" };
    if (SystName == "TOT_BTAG_") syst = { "BTAG_UP", "BTAG_DOWN", "BTAG_PT_UP", "BTAG_PT_DOWN", "BTAG_ETA_UP", "BTAG_ETA_DOWN" };
    if (SystName == "TOT_BTAG_LJET_") syst = { "BTAG_LJET_UP", "BTAG_LJET_DOWN", "BTAG_LJET_PT_UP", "BTAG_LJET_PT_DOWN", "BTAG_LJET_ETA_UP", "BTAG_LJET_ETA_DOWN" };
    if (SystName == "TOT_TRIG_") syst = { "TRIG_UP", "TRIG_DOWN", "TRIG_ETA_UP", "TRIG_ETA_DOWN" };
    if (SystName == "TOT_BFRAG_") syst = { "BFRAG_UP", "BFRAG_DOWN", "BFRAG_PETERSON" };
    if (SystName == "TOT_COLORREC_") syst = { "ERDON", "ERDONRETUNE", "GLUONMOVETUNE" };


    ///Function to get the error of a certain systematic: Sqrt(0.5*(Err_Up*Err_Up + Err_Down*Err_Down))/Nominal

    std::vector<double> BinCenters, BinLowEdge, BinHigEdge, CentralValues_Nom, CentralValues_Up, CentralValues_Down;
    std::vector<double> RelativeError;
    TString Dummy = "";
    double dummy = 0;

    ifstream NominalFile;
    string nominalfilename = string("UnfoldingResults/Nominal/" + Channel + "/" + Variable + "Results.txt");
    if (Channel == "combined" && performUnfoldingInCombinedChannel) {
        nominalfilename = string("UnfoldingResults/Nominal/" + Channel + "Unfolded/" + Variable + "Results.txt");
    }

    NominalFile.open(nominalfilename);
    if (!NominalFile.is_open()) {
        std::cout << "Nominal: " << nominalfilename << std::endl;
        std::cout << "The input file cannot be opened. Exiting!!" << std::endl;
        exit(433);
    }
    for (Int_t bin = 0; bin < bins; bin++) {
        double XAxisbinCenters_Nom = 0, XLowEdge_Nom = 0, XHigEdge_Nom = 0, CentralValue_Nom = 0;
        NominalFile >> Dummy >> XAxisbinCenters_Nom >> Dummy >> XLowEdge_Nom >> Dummy >> XHigEdge_Nom >> Dummy >> CentralValue_Nom >> Dummy >> dummy >> Dummy >> dummy;

        BinCenters.push_back(XAxisbinCenters_Nom);
        BinLowEdge.push_back(XLowEdge_Nom);
        BinHigEdge.push_back(XHigEdge_Nom);
        CentralValues_Nom.push_back(CentralValue_Nom);
    }
    NominalFile.close();

    CentralValues_Up = CentralValues_Nom;
    CentralValues_Down = CentralValues_Nom;

    for (size_t iter = 0; iter < syst.size(); iter++) {

        ifstream VarFile;
        string varfilename = string("UnfoldingResults/" + syst.at(iter) + "/" + Channel + "/" + Variable + "Results.txt");
        if (Channel == "combined" && performUnfoldingInCombinedChannel) {
            varfilename = string("UnfoldingResults/" + syst.at(iter) + "/" + Channel + "Unfolded/" + Variable + "Results.txt");
        }

        VarFile.open(varfilename);
        if (!VarFile.is_open()) {
            std::cout << "Variation: " << varfilename << std::endl;
            std::cout << "The input file cannot be opened. Exiting!!" << std::endl;
            exit(433);
        }
        for (Int_t bin = 0; bin < bins; bin++) {
            double CentralValue_Var = 0;
            VarFile >> Dummy >> dummy >> Dummy >> dummy >> Dummy >> dummy >> Dummy >> CentralValue_Var >> Dummy >> dummy >> Dummy >> dummy;

            if (syst.at(iter) == "PSFSRSCALE_UP" || syst.at(iter) == "PSFSRSCALE_DOWN") CentralValue_Var = CentralValues_Nom.at(bin) + (CentralValue_Var - CentralValues_Nom.at(bin)) / scaleUncReductionFactor;
            if (CentralValue_Var < CentralValues_Down.at(bin)) CentralValues_Down.at(bin) = CentralValue_Var;
            if (CentralValue_Var > CentralValues_Up.at(bin)) CentralValues_Up.at(bin) = CentralValue_Var;
        }
        VarFile.close();

    } //end syst loop

    for (Int_t bin = 0; bin < bins; bin++) {

        if (CentralValues_Nom.at(bin) != 0) {
            double up = std::fabs(CentralValues_Up.at(bin) - CentralValues_Nom.at(bin));
            double down = std::fabs(CentralValues_Down.at(bin) - CentralValues_Nom.at(bin));
            double rel_err = 0.5 * (up + down) / CentralValues_Nom.at(bin);
            RelativeError.push_back(rel_err);
        }
    }

    ofstream SystematicRelError;
    if (Channel == "combined" && performUnfoldingInCombinedChannel) SystematicRelError.open(ttbar::assignFolder("UnfoldingResults", "combinedUnfolded", SystName) + Variable + "Results.txt");
    else SystematicRelError.open(ttbar::assignFolder("UnfoldingResults", Channel, SystName) + Variable + "Results.txt");

    if (!SystematicRelError.is_open()) {
        std::cout << "The output file cannot be opened. Exiting!!" << std::endl;
        exit(434);
    }
    for (int bin = 0; bin < (int)RelativeError.size(); bin++) {
        SystematicRelError << "XAxisbinCenters[bin]: " << BinCenters.at(bin) << " bin: " << BinLowEdge.at(bin) << " to " << BinHigEdge.at(bin) << " SystematicRelError: " << RelativeError.at(bin) << std::endl;
    }
    SystematicRelError.close();

    return;
}



double Plotter13TeV::CalculateIntegral(TGraphAsymmErrors* tga_DiffXSecPlot, double Xbins[])
{
    double tmp_integral = 0;
    for (int i = 0; i < tga_DiffXSecPlot->GetN(); i++)
    {
        tmp_integral += tga_DiffXSecPlot->GetY()[i] * (Xbins[i + 1] - Xbins[i]);
    }
    std::cout << "Integral = " << tmp_integral << std::endl;
    return tmp_integral;
}



void Plotter13TeV::setTheoryStyleAndFillLegend(TH1* histo, TString theoryName, TLegend* leg) {

    histo->GetXaxis()->SetTitleOffset(1.08);
    histo->GetXaxis()->SetTitleSize(0.05);
    histo->GetXaxis()->SetLabelFont(42);
    histo->GetXaxis()->SetLabelOffset(0.007);
    histo->GetXaxis()->SetLabelSize(0.04);

    histo->GetYaxis()->SetTitleOffset(1.7);
    histo->GetYaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetLabelFont(42);
    histo->GetYaxis()->SetLabelOffset(0.007);
    histo->GetYaxis()->SetLabelSize(0.04);

    histo->GetXaxis()->SetLabelSize(0.04);//0.05 for njets; rest: 0.04


    histo->SetLineWidth(2);
    if (theoryName != "data") {
        histo->SetMarkerSize(0);
        histo->SetMarkerStyle(1);
    }

    if (theoryName == "madgraph") {
        histo->SetLineColor(kRed + 1);
        histo->SetLineStyle(1);
        //if(leg) leg->AddEntry(histo, "PowhegV2+P8",  "l");
        if (leg) leg->AddEntry(histo, "Powheg+Pythia8", "l");
    }
    if (theoryName == "powhegpythia") {
        histo->SetLineColor(kGreen + 1);
        histo->SetLineStyle(7);
        //if(leg) leg->AddEntry(histo, "aMC@NLOFxFx+P8",  "l");
        if (leg) leg->AddEntry(histo, "aMC@NLO+Pythia8", "l");
    }
    if (theoryName == "powhegherwig") {
        histo->SetLineColor(kBlue);
        histo->SetLineStyle(4);
        //if(leg) leg->AddEntry(histo, "MadGraphMLM+P8",  "l");
        if (leg) leg->AddEntry(histo, "MadGraph+Pythia8", "l");
    }
    if (theoryName == "mcatnloherwig") {
        histo->SetLineColor(kViolet - 5);
        histo->SetLineStyle(5);
        //if(leg) leg->AddEntry(histo, "PowhegV2+H++",  "l");
        if (leg) leg->AddEntry(histo, "Powheg+Herwig++", "l");
    }
    if (theoryName == "ahrens") {
        histo->SetLineColor(kViolet - 6);
        histo->SetLineStyle(6);
        if (leg) leg->AddEntry(histo, "NLO+NNLL", "l");
    }
    if (theoryName == "kidonakis") {
        histo->SetLineColor(kViolet - 6);
        histo->SetLineStyle(2);
        if (leg) leg->AddEntry(histo, "Approx. NNLO", "l");
    }
    if (theoryName == "matchup" || theoryName == "mass175.5") {
        histo->SetLineStyle(7);
        histo->SetLineColor(kPink - 7);
        if (leg && theoryName == "matchup") leg->AddEntry(histo, "Matching up", "l");
        if (leg && theoryName == "mass175.5") leg->AddEntry(histo, "Mass = 175.5 GeV", "l");
    }
    if (theoryName == "matchdown" || theoryName == "mass169.5") {
        histo->SetLineStyle(7);
        histo->SetLineColor(kRed - 7);
        if (leg && theoryName == "matchdown") leg->AddEntry(histo, "Matching down", "l");
        if (leg && theoryName == "mass169.5")   leg->AddEntry(histo, "Mass = 169.5 GeV", "l");
    }
    if (theoryName == "scaleup" || theoryName == "mass173.5") {
        histo->SetLineStyle(2);
        histo->SetLineColor(kAzure + 2);
        if (leg && theoryName == "scaleup") leg->AddEntry(histo, "4*Q^{2}", "l");
        if (leg && theoryName == "mass173.5") leg->AddEntry(histo, "Mass = 173.5 GeV", "l");
    }
    if (theoryName == "scaledown" || theoryName == "mass171.5") {
        histo->SetLineStyle(2);
        histo->SetLineColor(8);
        if (leg && theoryName == "scaledown") leg->AddEntry(histo, "Q^{2}/4", "l");
        if (leg && theoryName == "mass171.5")   leg->AddEntry(histo, "Mass = 171.5 GeV", "l");
    }
    if (theoryName == "mass178.5") {
        histo->SetLineStyle(3);
        histo->SetLineColor(42);
        if (leg && theoryName == "matchup") leg->AddEntry(histo, "Matching up", "l");
        if (leg && theoryName == "mass178.5") leg->AddEntry(histo, "Mass = 178.5 GeV", "l");
    }
    if (theoryName == "mass166.5") {
        histo->SetLineStyle(3);
        histo->SetLineColor(48);
        if (leg && theoryName == "matchdown") leg->AddEntry(histo, "Matching down", "l");
        if (leg && theoryName == "mass166.5")   leg->AddEntry(histo, "Mass = 166.5 GeV", "l");
    }

}

bool Plotter13TeV::addQCDToControlPlot()const
{
    if ((name.Contains("_step") && !name.Contains("7") && !name.Contains("8")) ||
        (name.Contains("events_") && !name.Contains("7") && !name.Contains("8")) ||
        (!name.Contains("Hyp") && (name.Contains("jetHT") || name.Contains("jetpT"))) ||
        name == "MET" || name.Contains("_noBTag") || name.Contains("_diLep") || name.Contains("Electron") || name.Contains("Muon")
        )
    {
        return 0;//1 in case you want QCD plotted
    }
    return 0;
}

void Plotter13TeV::GetEnvelopeVariationsForBandsInCPs(TH1* varHist, const TString Systematic, const TString channel)
{

    if (!varHist) return;

    std::vector<TString> syst;

    if (!Systematic.Contains("TOT_SCALE_") && !Systematic.Contains("TOT_BFRAG_") && !Systematic.Contains("TOT_COLORREC_")) {
        std::cout << "Systematics name: " << Systematic << std::endl;
        std::cout << "The calculation of envelope is not supported. Exiting!!" << std::endl;
        exit(433);
    }

    if (Systematic.Contains("TOT_SCALE_")) syst = { "MESCALE_UP", "MESCALE_DOWN", "MEFACSCALE_UP", "MEFACSCALE_DOWN", "MERENSCALE_UP", "MERENSCALE_DOWN",
    "PSISRSCALE_UP", "PSISRSCALE_DOWN", "PSFSRSCALE_UP", "PSFSRSCALE_DOWN" };
    if (Systematic.Contains("TOT_BFRAG_")) syst = { "BFRAG_UP", "BFRAG_DOWN", "BFRAG_PETERSON" };
    if (Systematic.Contains("TOT_COLORREC_")) syst = { "ERDON", "ERDONRETUNE", "GLUONMOVETUNE" };

    unsigned int nbins = varHist->GetNbinsX();

    TString file_nom = TString("Plots/Nominal/" + channel + "/" + name + "_source.root");
    ifstream inputFileStream(file_nom);
    if (inputFileStream.fail()) { return; }
    inputFileStream.close();

    TH1D* tmp_nom = nullptr;
    tmp_nom = fileReader->GetClone<TH1D>(file_nom, varHist->GetName(), 1);
    tmp_nom->SetDirectory(0);

    for (size_t iter = 0; iter < syst.size(); iter++) {

        TString file = TString("Plots/" + syst.at(iter) + "/" + channel + "/" + name + "_source.root");
        ifstream inputFileStream(file);

        if (inputFileStream.fail()) { continue; }
        inputFileStream.close();

        TH1D* tmp = nullptr;

        TString thisName = varHist->GetName();
        tmp = fileReader->GetClone<TH1D>(file, thisName, 1);

        tmp->SetDirectory(0);

        for (size_t ibin = 1; ibin <= nbins; ++ibin) {

            double tmp_bin = tmp->GetBinContent(ibin);
            if (syst.at(iter).Contains("PSFSRSCALE_")) tmp_bin = tmp_nom->GetBinContent(ibin) + (tmp_bin - tmp_nom->GetBinContent(ibin)) / scaleUncReductionFactor;

            if (Systematic.Contains("_DOWN")) {
                if (tmp_bin < varHist->GetBinContent(ibin))
                    varHist->SetBinContent(ibin, tmp_bin);
            }

            else {
                if (tmp_bin > varHist->GetBinContent(ibin))
                    varHist->SetBinContent(ibin, tmp_bin);
            }

        } //end bin loop

        delete tmp;

    } //end syst loop

    delete tmp_nom;

    return;
}


void Plotter13TeV::getSignalUncertaintyBand(TH1* uncBand, TString channel_)
{
    if (!uncBand)  return;
    std::vector<TString> syst{ "MASS_", "BSEMILEP_", "UETUNE_", "MATCH_", "PDF_ALPHAS_", "JES_", "JER_", "PU_", "LEPT_", "TRIG_",
                               "BTAG_", "BTAG_LJET_", "TOT_SCALE_", "TOT_BFRAG_", "TOT_COLORREC_" /*, "HAD_", "MOD_", "BTAG_PT_", "BTAG_ETA_", "BTAG_LJET_PT_", "BTAG_LJET_ETA_",*/ };

    double norm_Events = uncBand->Integral(-1e6, 1e6);
    int nbins = uncBand->GetNbinsX();
    std::vector<double> vec_varup(nbins, 0), vec_vardown(nbins, 0);
    for (size_t iter = 0; iter < syst.size(); iter++)
    {
        TString file_up = "", file_do = "";
        if (syst.at(iter) == "HAD_") {
            file_up = TString("Plots/MCATNLO/" + channel_ + "/" + name + "_source.root");
            file_do = TString("Plots/MCATNLO/" + channel_ + "/" + name + "_source.root");
        }
        else if (syst.at(iter) == "MOD_") {
            file_up = TString("Plots/POWHEG/" + channel_ + "/" + name + "_source.root");
            file_do = TString("Plots/POWHEG/" + channel_ + "/" + name + "_source.root");
        }
        else {
            file_up = TString("Plots/" + syst.at(iter) + "UP/" + channel_ + "/" + name + "_source.root");
            file_do = TString("Plots/" + syst.at(iter) + "DOWN/" + channel_ + "/" + name + "_source.root");
        }
        ifstream inputFileStream(file_up);
        if (inputFileStream.fail()) { continue; }
        inputFileStream.close();
        ifstream inputFileStream2(file_do);
        if (inputFileStream2.fail()) { continue; }
        inputFileStream2.close();

        // This lines crashes the code, some probles arises form the HistoListReader class
        TH1D* tmpUp = nullptr;
        TH1D* tmpDo = nullptr;

        if (drawTotalMCErrorForCP) {
            tmpUp = fileReader->GetClone<TH1D>(file_up, name + "_allmc", 1);
            tmpDo = fileReader->GetClone<TH1D>(file_do, name + "_allmc", 1);
        }

        else {
            tmpUp = fileReader->GetClone<TH1D>(file_up, name + "_allttbar", 1);
            tmpDo = fileReader->GetClone<TH1D>(file_do, name + "_allttbar", 1);
        }

        if (!tmpUp && tmpDo) { delete tmpDo; continue; }
        if (tmpUp && !tmpDo) { delete tmpUp; continue; }
        if (!tmpUp && !tmpDo) continue;

        if (nbins != tmpUp->GetNbinsX() || nbins != tmpDo->GetNbinsX()) {
            continue;
        }
        tmpDo->SetDirectory(0);
        tmpUp->SetDirectory(0);

        double upIntegral = tmpUp->Integral(-1e6, 1e6);
        double doIntegral = tmpDo->Integral(-1e6, 1e6);

        double factorVar = 1.0;
        if (syst.at(iter) == "MASS_") factorVar = massUncReductionFactor;

        tmpUp->Scale(norm_Events / upIntegral);
        tmpDo->Scale(norm_Events / doIntegral);
        for (Int_t nbin = 0; nbin < nbins; nbin++)
        {
            double binContent = uncBand->GetBinContent(nbin + 1);
            double rel_diffup = std::abs(tmpUp->GetBinContent(nbin + 1) - binContent) / binContent;
            double rel_diffdo = std::abs(tmpDo->GetBinContent(nbin + 1) - binContent) / binContent;
            if (binContent < 1e-6) {
                rel_diffdo = 0;
                rel_diffup = 0;
            }
            vec_varup.at(nbin) += (rel_diffup / factorVar) * (rel_diffup / factorVar);
            vec_vardown.at(nbin) += (rel_diffdo / factorVar) * (rel_diffdo / factorVar);
        }
        factorVar = 1.0;
        delete tmpDo; delete tmpUp;
    }
    for (size_t iter = 0; iter < vec_varup.size(); iter++)
    {
        double centralValue = uncBand->GetBinContent(iter + 1);
        uncBand->SetBinError(iter + 1, centralValue * 0.5 * (std::sqrt(vec_vardown.at(iter)) + std::sqrt(vec_varup.at(iter))));
    }
}

void Plotter13TeV::performPoissonSmearing(TH1* hist)
{

    int seed = 4357;
    if (doClosureTestWithLowStat) seed = 0;

    TRandom3* randClosureTest = new TRandom3(seed);
    for (int binIndex = 0; binIndex <= hist->GetNbinsX() + 1; ++binIndex) { // "0" and "+1" needed for under/overflow bins
        double binContent = hist->GetBinContent(binIndex);
        double newContent = randClosureTest->PoissonD(binContent);
        hist->SetBinContent(binIndex, newContent);
    }
    delete randClosureTest;
}


void Plotter13TeV::yRangeControlPlotRatio(double& yminCP_, double& ymaxCP_)const
{
    yminCP_ = 0.51;
    ymaxCP_ = 1.49;

    if (name.Contains("ToppT")) { yminCP_ = 0.51; ymaxCP_ = 1.49; }
    if (name.Contains("TopRapidity")) { yminCP_ = 0.51; ymaxCP_ = 1.49; }
    if (name.Contains("TTBarDeltaRapidity")) { yminCP_ = 0.51; ymaxCP_ = 1.49; }

    if (name.Contains("TTBarpT")) { yminCP_ = 0.51; ymaxCP_ = 1.49; }
    if (name.Contains("TTBarMass")) { yminCP_ = 0.51; ymaxCP_ = 1.49; }
    if (name.Contains("TTBarRapidity")) { yminCP_ = 0.51; ymaxCP_ = 1.49; }

    if (name.Contains("LeptonpT")) { yminCP_ = 0.51; ymaxCP_ = 1.49; }
    if (name.Contains("LeptonEta")) { yminCP_ = 0.51; ymaxCP_ = 1.49; }

    if (name.Contains("LLBarMass")) { yminCP_ = 0.51; ymaxCP_ = 1.49; }
    if (name.Contains("LLBarpT")) { yminCP_ = 0.51; ymaxCP_ = 1.49; }
    if (name.Contains("LLBarDPhi")) { yminCP_ = 0.51; ymaxCP_ = 1.49; }

    if (name.Contains("BJetpT")) { yminCP_ = 0.51; ymaxCP_ = 1.49; }
    if (name.Contains("BJetEta")) { yminCP_ = 0.51; ymaxCP_ = 1.49; }

    if (name.Contains("BBBarMass")) { yminCP_ = 0.51; ymaxCP_ = 1.49; }
    if (name.Contains("BBBarpT")) { yminCP_ = 0.51; ymaxCP_ = 1.49; }

    if (name.Contains("LeptonBjetMass")) { yminCP_ = 0.51; ymaxCP_ = 1.49; }
}



void Plotter13TeV::setResultRatioRanges(double& yminRes_, double& ymaxRes_)const
{
    yminRes_ = 0.51;
    ymaxRes_ = 1.49;

    if (name.Contains("TopPartonFraction")) { yminRes_ = 0.75; ymaxRes_ = 1.35; return; }
    if (name.Contains("ToppTTTRestFrame")) { yminRes_ = 0.70; ymaxRes_ = 1.65; return; }
    if (name.Contains("ToppTNLead")) { yminRes_ = 0.70; ymaxRes_ = 1.65; return; }
    if (name.Contains("ToppTLead")) { yminRes_ = 0.70; ymaxRes_ = 1.65; return; }
    if (name.Contains("ToppT")) { yminRes_ = 0.51; ymaxRes_ = 1.49; return; }

    if (name.Contains("TopRapidity")) { yminRes_ = 0.71; ymaxRes_ = 1.29; return; }

    if (name.Contains("TTBarDeltaRapidity")) { yminRes_ = 0.51; ymaxRes_ = 1.49; return; }
    if (name.Contains("TTBarDeltaPhi")) { yminRes_ = 0.85; ymaxRes_ = 1.25; return; }

    if (name.Contains("TTBarpT")) { yminRes_ = 0.51; ymaxRes_ = 1.49; return; }
    if (name.Contains("TTBarMass")) { yminRes_ = 0.51; ymaxRes_ = 1.49; return; }
    if (name.Contains("TTBarRapidity")) { yminRes_ = 0.71; ymaxRes_ = 1.29; return; }

    if (name.Contains("LeptonpT")) { yminRes_ = 0.71; ymaxRes_ = 1.29; return; }
    if (name.Contains("LeptonEta")) { yminRes_ = 0.21; ymaxRes_ = 1.79; return; }

    if (name.Contains("LLBarMass")) { yminRes_ = 0.85; ymaxRes_ = 1.25; return; }
    if (name.Contains("LLBarpT")) { yminRes_ = 0.85; ymaxRes_ = 1.25; return; }
    if (name.Contains("LLBarDPhi")) { yminRes_ = 0.81; ymaxRes_ = 1.19; return; }
    if (name.Contains("BBBarDPhi")) { yminRes_ = 0.71; ymaxRes_ = 1.29; return; }

    if (name.Contains("BJetpT")) { yminRes_ = 0.21; ymaxRes_ = 1.79; return; }
    if (name.Contains("BJetEta")) { yminRes_ = 0.21; ymaxRes_ = 1.79; return; }

    if (name.Contains("BBBarMass")) { yminRes_ = 0.90; ymaxRes_ = 1.14; return; }
    if (name.Contains("BBBarpT")) { yminRes_ = 0.60; ymaxRes_ = 1.49; return; }

    if (name.Contains("LeptonBjetMass")) { yminRes_ = 0.70; ymaxRes_ = 1.45; return; }
}
