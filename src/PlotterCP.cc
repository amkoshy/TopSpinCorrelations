#include <fstream>
#include <iostream>
#include <cstdio>
#include <sstream>
#include <cmath>
#include <memory>

#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TExec.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TList.h>
#include <THStack.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TError.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TRandom3.h>
#include <TMatrixD.h>

#include "../../common/include/RootFileReader.h"
#include "../../common/include/plotterUtils.h"

#include "DilepSVDFunctions.h"
#include "PlotterCP.h"
#include "PlotterConfigurationHelper.h"
#include "utils.h"
#include "Samples.h"

using namespace std;


PlotterCP::PlotterCP(TString analysis_year, TString treat_correlated_year)
{
    year = analysis_year;
    yearCorr = treat_correlated_year;

    outpath = "";
    outpathPlots = "Plots_"+year+"/";
    subfolderChannel = "";
    subfolderSpecial = "";
    gErrorIgnoreLevel = 1001;
    name = "defaultName";
    specialComment = "Standard";
    rangemin = 0;
    rangemax = 3;
    YAxis = "N_{events}";
    initialized = false;
    datafiles = 0;
    channelLabel.insert(channelLabel.begin(), 4, "");
    fileReader = RootFileReader::getInstance();


    isCP_ = false;
}


void PlotterCP::setOptions(TString name_, TString specialComment_, TString YAxis_, TString XAxis_,
                           int rebin_, bool doDYScale_, bool logX_, bool logY_, double ymin_, double ymax_, double rangemin_, double rangemax_,
                           int bins_, std::vector<double> XAxisbins_, std::vector<double> XAxisbinCenters_)
{
    name = name_;//Variable name you want to plot
    specialComment = specialComment_;
    YAxis = YAxis_;//Y-axis title
    XAxis = XAxis_;//X-axis title
    rebin = rebin_;//Nr. of bins to be merged together
    doDYScale = doDYScale_;//Boolean to control application of DY SFs
    logX = logX_;//Draw X-axis in Log scale
    logY = logY_;//Draw Y-axis in Log scale
    ymin = ymin_;//Min. value in Y-axis
    ymax = ymax_;//Max. value in Y-axis
    rangemin = rangemin_;//Min. value in X-axis
    rangemax = rangemax_;//Max. value in X-axis
    bins = bins_;//Number of bins to plot
    XAxisbins.clear();
    XAxisbins = XAxisbins_;// Bins edges = bins+1
    XAxisbinCenters.clear();
    XAxisbinCenters = XAxisbinCenters_;//Central point for BinCenterCorrection = bins

    //Modify the X/Y-axis labels
    if (XAxis.Contains("band#bar{b}")) {
        XAxis.ReplaceAll("band#bar{b}", 11, "b and #bar{b}", 13);
    }
    if (XAxis.Contains("tand#bar{t}")) {
        XAxis.ReplaceAll("tand#bar{t}", 11, "t and #bar{t}", 13);
    }
    if (XAxis.Contains("l^{+}andl^{-}")) {
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

    DYScale.insert(DYScale.begin(), 4, 1.);//Set DY scale factor array to (1., 1., 1., 1.)

    //Get configurations from the helper class
    configHelper = new PlotterConfigurationHelper(fileReader, year, yearCorr, doDYScale);

    this->energy = configHelper->energy;
    this->lumi = configHelper->lumi;
    this->lumi_fb = configHelper->lumi_fb.str();
    this->lumiError = configHelper->lumiError;
    this->topxsec = configHelper->topxsec;
    std::copy(std::begin(configHelper->bRatio), std::end(configHelper->bRatio), std::begin(this->bRatio));
    this->bRatioError = configHelper->bRatioError;

    this->revokeAntiQuantityCombination = configHelper->revokeAntiQuantityCombination;
    this->usePoissonSmearedPseudoData = configHelper->usePoissonSmearedPseudoData;
    this->estimateSysUncFromPseudoData = configHelper->estimateSysUncFromPseudoData;

    this->drawTotalMCErrorForCP = configHelper->drawTotalMCErrorForCP;
    this->massUncReductionFactor = configHelper->massUncReductionFactor;
    this->fsrScaleUncReductionFactor = configHelper->fsrScaleUncReductionFactor;

    this->cmsPrelimLabelMode = configHelper->cmsPrelimLabelMode;
    this->drawPlotRatio = configHelper->drawPlotRatio;

    std::vector<TString> allSysExp = configHelper->GetExpSys(false);
    std::vector<TString> allSysModW = configHelper->GetModWSys(false);
    std::vector<TString> allSysModI = configHelper->GetModISys(false);
    this->uncertaintyBandSystematics.clear();
    this->uncertaintyBandSystematics.insert(this->uncertaintyBandSystematics.end(), allSysExp.begin(), allSysExp.end());
    this->uncertaintyBandSystematics.insert(this->uncertaintyBandSystematics.end(), allSysModW.begin(), allSysModW.end());
    this->uncertaintyBandSystematics.insert(this->uncertaintyBandSystematics.end(), allSysModI.begin(), allSysModI.end());

    if (useEnvelopesInSystematicsList_) {
        std::vector<TString> temp_uncertaintyBandSystematics = this->uncertaintyBandSystematics;
        this->uncertaintyBandSystematics.clear();
        this->uncertaintyBandSystematics = this->configHelper->FoldSystematicVarsListIntoEnvelops(temp_uncertaintyBandSystematics);
    }

}


void PlotterCP::setDataSet(std::vector<TString> dataset_, std::vector<double> scales_, std::vector<TString> legends_, std::vector<int> colors_, TString DYEntry_)
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


void PlotterCP::setDataSet(TString mode, TString Systematic)
{
    initialized=false;
    legendsSyst.clear();

    if (channelLabel.size() < 4) channelLabel.insert(channelLabel.begin(), 4, "");
    if (mode == "ee") {channelType=0; channelLabel.at(0) = "ee";}
    if (mode == "mumu") {channelType=1; channelLabel.at(1) = "#mu#mu";}
    if (mode == "emu") {channelType=2; channelLabel.at(2) = "e#mu";}
    if (mode == "combined") {channelType=3; channelLabel.at(3) = "Dilepton";}

    // Set dataset specific subfolders
    outpathPlots = configHelper->outputPlotsFolder;
    subfolderChannel = mode;
    subfolderChannel.Prepend("/");
    subfolderSpecial = "";

    DYEntry = "Z+jets";

    if (Systematic.Contains("DY_") || Systematic.Contains("BG_") || Systematic.Contains("TW_")) Systematic = "Nominal";//need to vary the nominal DY and BG systematics

    TString tempSystematic = Systematic;
    if (Systematic.Contains("161718")) {
      tempSystematic = Systematic;
      configHelper->RemoveYearCorr(year, yearCorr, tempSystematic);
    }

    if (mode == "ee" && tempSystematic.Contains("MUON"))
      tempSystematic = "Nominal";

    if (mode == "mumu" && tempSystematic.Contains("ELE"))
      tempSystematic = "Nominal";

    TString histoListName = configHelper->baseFileList + "HistoFileList_" + tempSystematic + "_" + mode + ".txt";
    if (estimateSysUncFromPseudoData) histoListName = "PseudoFileLists/HistoFileList_" + tempSystematic + "_" + mode + ".txt";
    std::cout << "reading " << histoListName << std::endl;

    // Readout check
    ifstream FileList(histoListName);
    if (FileList.fail()) {
        std::cerr << "Error reading " << histoListName << std::endl;
        exit(1);
    }
    FileList.close();

    PlotterConfigurationHelper::fillLegendColorDataset(histoListName, legends, colors, dataset);

    if (usePoissonSmearedPseudoData) {
        if (!estimateSysUncFromPseudoData) std::cout << "Constructing pseudo-data as sum of 'nominal' ttbar and bkgd simulations" << std::endl;
        else std::cout << "Constructing pseudo-data as sum of 'varied' ttbar and bkgd simulations" << std::endl;

        TString nomFileListName = configHelper->baseFileList + "HistoFileList_Nominal_" + mode + ".txt";
        if (estimateSysUncFromPseudoData) nomFileListName = configHelper->baseFileList + "HistoFileList_" + tempSystematic + "_" + mode + ".txt";
        std::cout << "--- reading " << nomFileListName << " for pseudo-data " << std::endl;

        std::ifstream nomFileList(nomFileListName);
        if (nomFileList.fail()) {
            std::cerr << " Error reading " << nomFileListName << std::endl;
            exit(1);
        }

        TString nomFileName;
        pseudoDataset.clear();
        while(!nomFileList.eof()) {
            nomFileList >> nomFileName;
            if (nomFileName == "") continue;//skip empty lines
            pseudoDataset.push_back(nomFileName);
        }
        nomFileList.close();
    }
}

bool PlotterCP::fillHisto(TString Channel, TString Systematic)
{
    TH1::AddDirectory(kFALSE);
    if (initialized) return true;
    hists.clear();
    histsForPseudoData.clear();

    for (unsigned int i = 0; i < dataset.size(); i++) {
        TH1D *hist = fileReader->GetClone<TH1D>(dataset.at(i), name, true);
        if (!hist) return false;
        TH1D *histForPseudoData = NULL;
        if (usePoissonSmearedPseudoData) {
            if (pseudoDataset.size() != dataset.size()) {
                std::cout << "FileLists for nominal dataset and pseudo-data contain different number of input .root files!" << std::endl;
                return false;
            }
            histForPseudoData = fileReader->GetClone<TH1D>(pseudoDataset.at(i), name, true);
            if (!histForPseudoData) return false;
        }

	TString stemp = "";
        if (!revokeAntiQuantityCombination &&
            !name.Contains("_step") && !name.Contains("bcp_") && !name.Contains("Lead") && !name.EndsWith("bkr") && !name.EndsWith("akr") && !name.EndsWith("afterbtag") &&
            (name.Contains("Lepton") || name.Contains("BJet") || name.Contains("Top")))
        {
            stemp = name;
	    //            if (name.Contains("Lepton")) stemp.ReplaceAll("Lepton", 6, "AntiLepton", 10);
	    //            else if (name.Contains("BJet")) stemp.ReplaceAll("BJet", 4, "AntiBJet", 8);
	    //            else if (name.Contains("Top")) stemp.ReplaceAll("Top", 3, "AntiTop", 7);

            const TH1D *other = fileReader->Get<TH1D>(dataset.at(i), stemp);
            if (other) hist->Add(other);
            else std::cerr << "Cannot find corresponding anti-quantity histogram: " << stemp<<std::endl;

            if (usePoissonSmearedPseudoData) {
                const TH1D *otherForPseudoData = fileReader->Get<TH1D>(pseudoDataset.at(i), stemp);
                if (otherForPseudoData) histForPseudoData->Add(otherForPseudoData);
                else std::cerr << "Cannot find corresponding anti-quantity histogram for pseudo-data: " << stemp << std::endl;
            }
        }

        //Rescaling to the data luminosity
	double LumiWeight;
	double year_corr_rescale_from_file = 1.;
	if (year == "fullRun2UL" || year == "2016") {
	  year_corr_rescale_from_file = configHelper->get_year_corr_rescale_from_file(Systematic, dataset.at(i));
	  configHelper->get_year_from_file(dataset.at(i));
	  LumiWeight = configHelper->CalcLumiWeight(dataset.at(i), configHelper->lumi_year);
	}
	else
	  LumiWeight = configHelper->CalcLumiWeight(dataset.at(i));
	configHelper->ApplyFlatWeights(hist, LumiWeight);
	if (year_corr_rescale_from_file != 1.) {
	  TString datasetNominal = configHelper->getCorrespondingNominalDataset(dataset.at(i), Systematic, Channel);
	  TH1D* histNominal = fileReader->GetClone<TH1D>(datasetNominal, name, true);
	  if (!histNominal) return false;
	  if (stemp != "") {
	    const TH1D *otherNominal = fileReader->Get<TH1D>(datasetNominal, stemp);
	    if (otherNominal) histNominal->Add(otherNominal);
	  }
	  double LumiWeightNominal = configHelper->CalcLumiWeight(datasetNominal, configHelper->lumi_year);
	  configHelper->ApplyFlatWeights(histNominal, LumiWeightNominal);
	  configHelper->RescaleHistYearToYearCorr(hist, histNominal, year_corr_rescale_from_file);
	}

        common::setHHStyle(*gStyle);
        hists.push_back(*hist);
        delete hist;
        if (usePoissonSmearedPseudoData) {
	    double LumiWeightForPseudoData;
	    double year_corr_rescale_from_file_pseudo = 1.;
	    if (year == "fullRun2UL" || year == "2016") {
	      year_corr_rescale_from_file_pseudo = configHelper->get_year_corr_rescale_from_file(Systematic, pseudoDataset.at(i));
	      configHelper->get_year_from_file(pseudoDataset.at(i));
	      LumiWeightForPseudoData = configHelper->CalcLumiWeight(pseudoDataset.at(i), configHelper->lumi_year);
	    }
	    else
	      LumiWeightForPseudoData = configHelper->CalcLumiWeight(pseudoDataset.at(i));
            configHelper->ApplyFlatWeights(histForPseudoData, LumiWeightForPseudoData);
	    if (year_corr_rescale_from_file != 1.) {
	      TString pseudoDatasetNominal = configHelper->getCorrespondingNominalDataset(pseudoDataset.at(i), Systematic, Channel);
	      TH1D* histForPseudoDataNominal = fileReader->GetClone<TH1D>(pseudoDatasetNominal, name, true);
	      if (!histForPseudoDataNominal) return false;
	      double LumiWeightForPseudoDataNominal = configHelper->CalcLumiWeight(pseudoDatasetNominal, configHelper->lumi_year);
	      configHelper->ApplyFlatWeights(histForPseudoDataNominal, LumiWeightForPseudoDataNominal);
	      configHelper->RescaleHistYearToYearCorr(histForPseudoData, histForPseudoDataNominal, year_corr_rescale_from_file_pseudo);
	    }
            histsForPseudoData.push_back(*histForPseudoData);
        }
        delete histForPseudoData;
    }

    if (usePoissonSmearedPseudoData) {
        if (channelType != 3) {// for ee, emu or mumu channel
            for (size_t i = 1; i < pseudoDataset.size(); ++i) {
	      if (!pseudoDataset.at(i).Contains("ttbarsignalplustau") && !pseudoDataset.at(i).Contains("ttbarsignalviatau") && !pseudoDataset.at(i).Contains("ttbarbgviatau") && !pseudoDataset.at(i).Contains("run")) {
                    hists[0].Add(&histsForPseudoData[i]);
                }
            }
            PlotterCP::performPoissonSmearing(&hists[0]);
        } else if (channelType == 3) {
            //Indeces must match to the order of data files in filelist {syst}_combined.txt
            TString chToInclude[3] = {"ee_", "emu_", "mumu_"};
            for (int ch_it = 0; ch_it < 3; ++ch_it) {
                for (size_t i = ch_it + 1; i < pseudoDataset.size(); ++i) {//skipping pseudo-data files (moving iterator by +1 when initialized) to which MC is added
		  if (!pseudoDataset.at(i).Contains("ttbarsignalplustau") && !pseudoDataset.at(i).Contains("ttbarsignalviatau") && !pseudoDataset.at(i).Contains("ttbarbgviatau") &&
                        !pseudoDataset.at(i).Contains("run") && pseudoDataset.at(i).Contains(chToInclude[ch_it]))
                    {
                        hists[ch_it].Add(&histsForPseudoData[i]);
                    }
                }
                PlotterCP::performPoissonSmearing(&hists[ch_it]);
            }
        }
    }

    initialized=true;
    return true;
}


void PlotterCP::DYScaleFactor(TString SpecialComment)
{
    DYScale = {1., 1., 1., 1.};
    if (!configHelper->doDYScaleFromFit)
      configHelper->DYScaleFactor(SpecialComment, DYScale, name);
    else
      configHelper->DYScaleFactorFromFit(SpecialComment, DYScale, name);
}


void PlotterCP::plotWriterWrapper(TString Channel, TString Systematic)
{
    writeCP(Channel, Systematic);
}


void PlotterCP::writeCP(TString Channel, TString Systematic) // do scaling, stacking, legending, and write in file
{
    setDataSet(Channel, Systematic);
    if (!fillHisto(Channel, Systematic)) return;
    if (hists.size() == 0) {
        std::cerr << "***ERROR! No histograms available! " << Channel << "/" << Systematic << std::endl;
        exit(11);
    }
    std::vector<TH1*> drawhists;
    drawhists.resize(hists.size());

    auto c = make_shared<TCanvas>("", "");
    c->Clear();
    c->SetName("c");
    c->SetTitle("");

    auto stack = make_shared<THStack>("def", "def");

    TLegend *leg = new TLegend();
    // TLegend *leg = new TLegend(0.1,0.1,0.2,0.2,"NDC");
    leg->Clear();
    int legchange=0;

    std::vector<double> XAxisbins = this->XAxisbins;//FIXME: former functionality for readout of x-binning vector for unfolding; check whether works for CPs, e.g. to introduce non-equidistant binning
    int bins = this->bins;
    double Xbins[XAxisbins.size()];
    for (unsigned int i = 0; i < XAxisbins.size(); i++) Xbins[i]=XAxisbins[i];

    TString newname = name;
    if(name.Contains("Hyp")) newname.ReplaceAll("Hyp", 3, "", 0);

    std::unique_ptr<TH1> sumttbar;
    std::unique_ptr<TH1> allttbar;

    if (hists.size() != legends.size()){
    	std::cerr << "Size of legends don't match with the histograms. This is usually related with a sample not included in the legend classification." << std::endl;
    	exit(11);
    }

    // Rebin and manipulate histograms per source
    for (unsigned int i = 0; i < hists.size(); i++) {
        drawhists[i] = (TH1D*) hists[i].Clone();
        if (rebin > 1) drawhists[i]->Rebin(rebin);

        //In case you would like to produce histogram with the inclusive last bin
        if (name.Contains("HypjetMulti_diLep") || name.Contains("HypjetMulti_noBTag") || name.Contains("HypjetMulti_step7")) {
            double incJetContent = 0.0;
            incJetContent = drawhists[i]->GetBinContent(7) + drawhists[i]->GetBinContent(8) + drawhists[i]->GetBinContent(9) + drawhists[i]->GetBinContent(10);
            drawhists[i]->SetBinContent(7, incJetContent);
            drawhists[i]->SetBinContent(8, 0.0);
            drawhists[i]->SetBinContent(9, 0.0);
            drawhists[i]->SetBinContent(10, 0.0);
        }
        else if (name.Contains("HypBjetMulti_noBTag") || name.Contains("HypBjetMulti_step7") || name.Contains("HypexjetMulti")) {
            double incJetContent = 0.0;
            for (int kk = 6; kk < drawhists[i]->GetNbinsX(); kk++) {
                incJetContent += drawhists[i]->GetBinContent(kk+1);
                drawhists[i]->SetBinContent(kk+1, 0.0);
            }
            incJetContent += drawhists[i]->GetBinContent(6);
            drawhists[i]->SetBinContent(6, incJetContent);
        }
        else if (name.Contains("basic_jet_multiplicity_step4")) {
            double incJetContent = 0.0;
            for (int kk = 10; kk < drawhists[i]->GetNbinsX(); kk++) {
                incJetContent += drawhists[i]->GetBinContent(kk+1);
                drawhists[i]->SetBinContent(kk+1, 0.0);
            }
            incJetContent += drawhists[i]->GetBinContent(10);
            drawhists[i]->SetBinContent(10, incJetContent);
        }

        if (XAxisbins.size() > 1) drawhists[i] = drawhists[i]->Rebin(bins, "tmp", Xbins);
        if (legends.at(i) != "Data") removeNegativeExpectations(drawhists[i]);
        setStyle(drawhists[i], i, true);

        if (legends.at(i) == "t#bar{t} signal") {
            if (sumttbar.get()) sumttbar->Add(drawhists[i]);
            else sumttbar = std::unique_ptr<TH1>{static_cast<TH1*>(drawhists[i]->Clone())};
        }
        if (legends.at(i) == "t#bar{t} signal" || legends.at(i) == "t#bar{t} other") {
            if (allttbar.get()) allttbar->Add(drawhists[i]);
            else allttbar = std::unique_ptr<TH1>{static_cast<TH1*>(drawhists[i]->Clone())};
        }
    }

    // Operate histograms stack, assign legend colours and rescale DY
    for (unsigned int i = 0; i < hists.size(); ++i) {
        if (legends.at(i) != "Data") {

	    // Normalization uncertainties rescaling
	    //      const bool is_tw_bkg = legends.at(i) == "Single t";
	    const bool is_dy_bkg = legends.at(i) == DYEntry;
	    const bool is_ttbar_sig = (legends.at(i) == "t#bar{t} signal");
	    if (is_dy_bkg && Systematic == "DY_UP") drawhists[i]->Scale(1.2);
	    else if (is_dy_bkg && Systematic == "DY_DOWN") drawhists[i]->Scale(0.8);
	    //     	    else if (is_tw_bkg && Systematic == "TW_UP") drawhists[i]->Scale(1.3);
	    //	    else if (is_tw_bkg && Systematic == "TW_DOWN") drawhists[i]->Scale(0.7);
	    //	    else if (Systematic.Contains("BG_UP") && !is_ttbar_sig && !is_dy_bkg && !is_tw_bkg) drawhists[i]->Scale(1.3);
	    //	    else if (Systematic.Contains("BG_DOWN") && !is_ttbar_sig && !is_dy_bkg && !is_tw_bkg) drawhists[i]->Scale(0.7);
	    else if (Systematic.Contains("BG_UP") && !is_ttbar_sig && !is_dy_bkg) drawhists[i]->Scale(1.3);
	    else if (Systematic.Contains("BG_DOWN") && !is_ttbar_sig && !is_dy_bkg) drawhists[i]->Scale(0.7);

            if ((legends.at(i) == DYEntry) && DYScale.at(channelType) != 1) drawhists[i]->Scale(DYScale.at(channelType));

            if (i > 0) {
                if (legends.at(i) != legends.at(i-1)) {
                    legchange = i;
                    if ((legends.at(i) == DYEntry) && DYScale.at(channelType) != 1) leg->AddEntry(drawhists[i], legends.at(i), "f");
                    else leg->AddEntry(drawhists[i], legends.at(i), "f");
                } else {
                    drawhists[legchange]->Add(drawhists[i]);
                }
            }

            if (i != (hists.size()-1)) {
                if (legends.at(i) != legends.at(i+1)) {
                    drawhists[i]->SetLineColor(1);
                }
            } else {
                drawhists[i]->SetLineColor(1);
            }

            if (legends.at(i) != legends.at(i-1)) {
                drawhists[i]->SetLineColor(1);
                if (!legends.at(i).Contains("QCD") || addQCDToControlPlot()) stack->Add(drawhists[i]);
            }
        } else {
            if (i==0) leg->AddEntry(drawhists[i], legends.at(i), "pe");
            if (i>0) {
                if (legends.at(i) != legends.at(i-1)) {
                    leg->AddEntry(drawhists[i], legends.at(i), "pe");
                }
                if (legends.at(i) == legends.at(0)) {
                    drawhists[0]->Add(drawhists[i]);
                }
            }
        }
    }

    // Calculation of inclusive cross section and preparation of event table is steered via "HypjetMultiXSec" histogram
    if (name.Contains("HypjetMultiXSec")) {
        double InclusiveXsectionWrite[4], InclusiveXsectionStatErrorWrite[4];
        CalcXSec(dataset, pseudoDataset, InclusiveXsectionWrite, InclusiveXsectionStatErrorWrite, Systematic, "");
        if (channelType==3) PlotterCP::PlotXSec(Channel);
        if (configHelper->makeControlPlotEventCountTable) {
            PlotterCP::MakeTable(Channel, Systematic);
        }
    }

    if (name.Contains(configHelper->EventCountTableWhenHist) && configHelper->makeControlPlotEventCountTable) {
        PlotterCP::MakeTable(Channel, Systematic);
    }

    // Configure legends
    TLegend *leg1 = (TLegend*)leg->Clone("leg1");
    TLegend *leg2 = (TLegend*)leg->Clone("leg2");
    TLegend *leg3 = (TLegend*)leg->Clone("leg3");
    setControlPlotLegendStyle(drawhists, legends, leg, leg1, leg2, leg3);

    // Configure ranges and axes in plots
    drawhists[0]->SetMinimum(ymin);
    if (rangemin != 0 || rangemax != 0) drawhists[0]->SetAxisRange(rangemin, rangemax, "X");
    if (logY) c->SetLogy();
    if(name.Contains("_vs_TTBarMass") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets")) c->SetGridx();
    if (ymax==0) {
        if (logY) drawhists[0]->SetMaximum(18  * drawhists[0]->GetBinContent(drawhists[0]->GetMaximumBin()));
        else drawhists[0]->SetMaximum(1.5 * drawhists[0]->GetBinContent(drawhists[0]->GetMaximumBin()));
    } else {
        double current_plot_ymax = ymax;
        drawhists[0]->SetMaximum(current_plot_ymax);
    }
    if (name.Contains("HypTopRapidity") || name.Contains("HypTTBarRapidity") || (name.Contains("HypAntiTopRapidity") && revokeAntiQuantityCombination)) drawhists[0]->GetXaxis()->SetNdivisions(511);
    // Jason
    if(name.Contains("_LLBar") || 
       name.Contains("_Lepton") || 
       name.Contains("_AntiLepton") || 
       name.Contains("_TTBarMass") ||
       name.Contains("_ScatteringAngle_TTBarFrame") ) {
      drawhists[0]->GetXaxis()->SetNdivisions(-20604);
    }
    drawhists[0]->GetXaxis()->SetNoExponent(kTRUE);
    TGaxis::SetMaxDigits(4);

    //Removal of extra ticks in plots
    if (name.Contains("HypJetMultpt")) {
        drawhists[0]->GetXaxis()->SetNdivisions(drawhists[0]->GetNbinsX(), 0, 0, 1);
        TString TitBin = "";
        for (int bin = 1; bin <= drawhists[0]->GetNbinsX(); bin++) {
            if (bin == drawhists[0]->GetNbinsX()) {
                TitBin += "#geq"; TitBin += drawhists[0]->GetBinCenter(bin);
                drawhists[0]->GetXaxis()->SetBinLabel(bin, TitBin);
            } else {
                TitBin += drawhists[0]->GetBinCenter(bin);
                drawhists[0]->GetXaxis()->SetBinLabel(bin, TitBin);
            }
            TitBin = "";
        }
    }
    else if (name.Contains("HypjetMulti_diLep") || name.Contains("HypjetMulti_noBTag") || name.Contains("HypjetMulti_step7")) {
        drawhists[0]->GetXaxis()->SetNdivisions(drawhists[0]->GetNbinsX(), 0, 0, 1);
        TString TitBin = "";
        for (int bin = 1; bin <= drawhists[0]->GetNbinsX(); bin++) {
            if( bin == 7) {
                TitBin += "#geq"; TitBin += drawhists[0]->GetBinCenter(bin);
                drawhists[0]->GetXaxis()->SetBinLabel(bin, TitBin);
            } else {
                TitBin += drawhists[0]->GetBinCenter(bin);
                drawhists[0]->GetXaxis()->SetBinLabel(bin, TitBin);
            }
            TitBin = "";
        }
    }
    else if (name.Contains("HypBjetMulti_noBTag") || name.Contains("HypBjetMulti_step7")) {
        drawhists[0]->GetXaxis()->SetNdivisions(drawhists[0]->GetNbinsX(), 0, 0, 1);
        TString TitBin = "";
        for (int bin = 1; bin <= drawhists[0]->GetNbinsX(); bin++) {
            if (bin == 6) {
                TitBin += "#geq"; TitBin += drawhists[0]->GetBinCenter(bin);
                drawhists[0]->GetXaxis()->SetBinLabel(bin, TitBin);
            } else {
                TitBin += drawhists[0]->GetBinCenter(bin);
                drawhists[0]->GetXaxis()->SetBinLabel(bin, TitBin);
            }
            TitBin = "";
        }
    }
    else if (name.Contains("basic_jet_multiplicity_step4")) {
        drawhists[0]->GetXaxis()->SetNdivisions(drawhists[0]->GetNbinsX(), 0, 0, 1);
        TString TitBin = "";
        for (int bin = 1; bin <= drawhists[0]->GetNbinsX(); bin++) {
            if (bin == 10) {
                TitBin += "#geq"; TitBin += drawhists[0]->GetBinCenter(bin);
                drawhists[0]->GetXaxis()->SetBinLabel(bin, TitBin);
            } else {
                TitBin += drawhists[0]->GetBinCenter(bin);
                drawhists[0]->GetXaxis()->SetBinLabel(bin, TitBin);
            }
            TitBin = "";
        }
    }

    //Add the binwidth to the yaxis in yield plots
    TString ytitle = drawhists[0]->GetYaxis()->GetTitle();
    double binwidth = drawhists[0]->GetXaxis()->GetBinWidth(1);
    std::ostringstream width;
    width<<binwidth;

    if (name.Contains("DeltaR") || name.Contains("Rapidity") || name.Contains("Eta") ||
        name.Contains("Phi") || name.Contains("Fraction") || (name.Contains("Mass") && name.Contains("1st")))
    {
      //        ytitle.Append(" / ").Append(width.str());
        ytitle = ytitle;
    }

    else if(
	    (name.Contains("_vs_TTBarMass", TString::kIgnoreCase ) ||
             name.Contains("_vs_ToppT", TString::kIgnoreCase ) ||
             name.Contains("_vs_ScatteringAngle", TString::kIgnoreCase ) ||
             name.Contains("_AntiLepton", TString::kIgnoreCase) || 
             name.Contains("_Lepton", TString::kIgnoreCase) || 
             name.Contains("_LLBar", TString::kIgnoreCase) || 

             name.Contains("Bj", TString::kIgnoreCase) || 
             name.Contains("Bk", TString::kIgnoreCase) || 
             name.Contains("Br", TString::kIgnoreCase) || 
             name.Contains("Bq", TString::kIgnoreCase) || 
             name.Contains("Bn", TString::kIgnoreCase) || 

             name.Contains("Ckk", TString::kIgnoreCase) || 
             name.Contains("Crr", TString::kIgnoreCase) || 
             name.Contains("Cnn", TString::kIgnoreCase) || 
             name.Contains("Ckj", TString::kIgnoreCase) || 
             name.Contains("Crq", TString::kIgnoreCase) || 

             name.Contains("Chan", TString::kIgnoreCase) || 
             name.Contains("Csca", TString::kIgnoreCase) || 
             name.Contains("Ctra", TString::kIgnoreCase) || 
             name.Contains("CkjL", TString::kIgnoreCase) || 
             name.Contains("CrqL", TString::kIgnoreCase) || 

             name.Contains("CPrk", TString::kIgnoreCase) || 
             name.Contains("CMrk", TString::kIgnoreCase) || 
             name.Contains("CPnr", TString::kIgnoreCase) || 
             name.Contains("CMnr", TString::kIgnoreCase) || 
             name.Contains("CPnk", TString::kIgnoreCase) || 
             name.Contains("CMnk", TString::kIgnoreCase) || 
             name.Contains("CPrj", TString::kIgnoreCase) || 
             name.Contains("CMrj", TString::kIgnoreCase) || 

             name.Contains("CrkP", TString::kIgnoreCase) || 
             name.Contains("CrkM", TString::kIgnoreCase) || 
             name.Contains("CnrP", TString::kIgnoreCase) || 
             name.Contains("CnrM", TString::kIgnoreCase) || 
             name.Contains("CnkP", TString::kIgnoreCase) || 
             name.Contains("CnkM", TString::kIgnoreCase) || 

             name.Contains("DPhi", TString::kIgnoreCase) || 
             name.Contains("cHel", TString::kIgnoreCase) || 
             name.Contains("cLab", TString::kIgnoreCase))
	    )
    {
        ytitle = ytitle;
    }

    else if(
	    (name.Contains("pT", TString::kIgnoreCase) && !name.Contains("JetMult")) || 
	    (name.Contains("Mass", TString::kIgnoreCase) && !name.Contains("1st")) || 
	    name.Contains("MET") || 
	    name.Contains("HT") && 
	    (!name.Contains("1st") && 
             !name.Contains("Rapidity", TString::kIgnoreCase) && 
             !name.Contains("Eta", TString::kIgnoreCase) && 
             !name.Contains("Phi", TString::kIgnoreCase) && 
             !name.Contains("JetMult") && 
             !name.Contains("Fraction") && 
             !name.Contains("Multiplicity", TString::kIgnoreCase))
	    )
    {
      //        ytitle.Append(" / ").Append(width.str()).Append(" GeV");
        ytitle = ytitle;

    }

    else {
      //      ytitle.Append(" / ").Append(width.str());
        ytitle = ytitle;
    }

    drawhists[0]->GetYaxis()->SetTitle(ytitle);
    drawhists[0]->Draw("e1");
    gStyle->SetEndErrorSize(0);

    // Draw histogram stack
    stack->Draw("same HIST");

    // Configure and draw uncertainty band
    TH1 *stacksum = common::summedStackHisto(stack.get());
    TH1 *uncBand = nullptr, *uncBandPlot = nullptr;
    std::map<TString, std::vector<double>> varUnc;

    if (mergeEnvelopeVariationsForBandsInCPs_) {
        sumttbar->SetName(name + "_signalmc");
        allttbar->SetName(name + "_allttbar");
        stacksum->SetName(name + "_allmc");
        GetEnvelopeVariationsForBandsInCPs(sumttbar.get(), Systematic, Channel);
        GetEnvelopeVariationsForBandsInCPs(allttbar.get(), Systematic, Channel);
        GetEnvelopeVariationsForBandsInCPs(stacksum, Systematic, Channel);
    }

    if (drawUncBand_) {
        if (drawTotalMCErrorForCP) uncBand = dynamic_cast<TH1*>(stacksum->Clone("uncBand"));
        else uncBand = dynamic_cast<TH1*>(allttbar->Clone("uncBand"));

        getSignalUncertaintyBand(uncBand, varUnc, Channel);
        uncBand->SetFillStyle(3354);
        uncBand->SetFillColor(kBlack);
        uncBand->SetLineColor(kBlack);
        uncBand->SetMarkerStyle(0);
        gStyle->SetHatchesLineWidth(1);
        gStyle->SetHatchesSpacing(0.8);

        uncBandPlot = dynamic_cast<TH1*> (uncBand->Clone("uncBandPlot"));
        for (int i = 0; i <= stacksum->GetNbinsX(); i++) {
            uncBandPlot->SetBinContent(i, stacksum->GetBinContent(i));
        }
        uncBandPlot->Draw("same,e2");
        leg->AddEntry(uncBand, "Uncertainty", "f");
        if (leg3) {
            leg3->AddEntry(uncBand, "Uncertainty", "f");
            const float y1 = leg3->GetY1NDC();
            const float y2 = leg3->GetY2NDC();
            const float deltaY = std::fabs(y2 - y1);
            const int nentriesLeg2 = leg2->GetNRows();
            const int nentriesLeg3 = leg3->GetNRows();
            //leg3->SetY1NDC(y2 - 1. * nentriesLeg3 / nentriesLeg2 * deltaY);
            leg3->SetY1NDC(0.7665217);
        }
    }

    // Redraw data, draw labels and legend etc.
    gPad->RedrawAxis();
    TExec *setex1 = new TExec("setex1", "gStyle->SetErrorX(0.5)");
    setex1->Draw();
    if (drawUncBand_) uncBandPlot->Draw("same,e2");
    TExec *setex2 = new TExec("setex2", "gStyle->SetErrorX(0.)");
    setex2->Draw();
    drawhists[0]->Draw("same,e1");

    drawCMSLabels(cmsPrelimLabelMode);
    drawDecayChLabel(channelLabel[channelType]);

    if (name.Contains("_step8_") && ( name.Contains("_topscatteringangle")  ||  name.Contains("_ttbarmass") ) )     draw2DLabels();
    else if (name.Contains("_vs_TTBarMass") || name.Contains("_vs_ScatteringAngle_TTBarFrame")  ||  name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets") )      draw2DLabels();

    if (name.Contains("JetMult")) {
        TString legtit = "";
        if (name.Contains("pt60")) legtit += "p_{T}^{jet} > 60 GeV";
        else if (name.Contains("pt100")) legtit += "p_{T}^{jet} > 100 GeV";
        else legtit += "p_{T}^{jet} > 30 GeV, |#eta^{jet}| < 2.4";
        leg->SetHeader(legtit);
    }

    if (leg1) leg1->Draw("same");
    if (leg2) leg2->Draw("same");
    if (leg3) leg3->Draw("same");


    if (drawPlotRatio) {
        double yminCP_ = 0.49, ymaxCP_ = 1.51;
        yRangeControlPlotRatio(yminCP_, ymaxCP_);
        common::drawRatio(drawhists[0], stacksum, uncBand, yminCP_, ymaxCP_, doFit_, name);
    }

    // Create directory for output plots
    TString outdir = utils::assignFolder(this->outpathPlots, Channel, Systematic);
    //    c->Print(outdir.Copy() + name + ".eps");
    //    c->Print(outdir.Copy() + name + ".C");
    c->Print(outdir.Copy() + name + ".pdf");

    // Get the ratio plot from the canvas
    TPad *tmpPad = dynamic_cast<TPad*>(c->GetPrimitive("rPad"));
    TH1 *ratio = nullptr;
    if (tmpPad) ratio = dynamic_cast<TH1*>(tmpPad->GetPrimitive("ratio"));

    //save canvas and sources in a root file
    TFile out_root(outdir.Copy() + name + "_source.root", "RECREATE");
    drawhists[0]->Write(name + "_data");
    sumttbar->Write(name + "_signalmc");
    allttbar->Write(name + "_allttbar");
    stacksum->SetName(name);
    stacksum->Write(name + "_allmc");

    if (ratio && ratio->GetEntries()) ratio->Write("ratio");
    c->Write(name + "_canvas");
    c->Clear();

    if (drawUncBand_) {
        plotUncertaintyBandContributions(varUnc, Channel, uncBand);
    }

    //

    TH1D * mt_ttsignal = 0 ;
    TH1D * mt_ttother  = 0 ;
    TH1D * mt_zjets    = 0 ;
    TH1D * mt_wjets    = 0 ;
    TH1D * mt_singlet  = 0 ;
    TH1D * mt_diboson  = 0 ;
    TH1D * mt_ttV      = 0 ;

    for ( unsigned int iHW = 0; iHW < drawhists.size(); ++iHW ){
      if ( iHW > 0 ){
	if ( legends.at(iHW) == legends.at(iHW-1) ) continue;
      }

      if (legends.at(iHW)=="t#bar{t} signal"){
	if (!mt_ttsignal) mt_ttsignal = (TH1D*) drawhists.at(iHW)->Clone();
	else mt_ttsignal->Add(drawhists.at(iHW));
      }
      else if (legends.at(iHW)=="t#bar{t} other"){
	if (!mt_ttother) mt_ttother = (TH1D*) drawhists.at(iHW)->Clone();
	else mt_ttother->Add(drawhists.at(iHW));
      }
      else if (legends.at(iHW)=="Z+jets"){
	if (!mt_zjets) mt_zjets = (TH1D*) drawhists.at(iHW)->Clone();
	else mt_zjets->Add(drawhists.at(iHW));
      }
      else if (legends.at(iHW)=="W+jets"){
	if (!mt_wjets) mt_wjets = (TH1D*) drawhists.at(iHW)->Clone();
	else mt_wjets->Add(drawhists.at(iHW));
      }
      else if (legends.at(iHW)=="Single t"){
	if (!mt_singlet) mt_singlet = (TH1D*) drawhists.at(iHW)->Clone();
	else mt_singlet->Add(drawhists.at(iHW));
      }
      else if (legends.at(iHW)=="Diboson"){
	if (!mt_diboson) mt_diboson = (TH1D*) drawhists.at(iHW)->Clone();
	else mt_diboson->Add(drawhists.at(iHW));
      }
      else if (legends.at(iHW)=="t#bar{t}+Z/W"){
	if (!mt_ttV) mt_ttV = (TH1D*) drawhists.at(iHW)->Clone();
	else mt_ttV->Add(drawhists.at(iHW));
      }
      else if (legends.at(iHW)!="Data") {
	std::cerr << "ERROR: Legend entry " << legends.at(iHW) << " not linked to any output histogram. Maybe a new process was added?" << std::endl;
	exit(1);
      }
    }

    if ( !mt_ttsignal || !mt_ttother || !mt_zjets || !mt_wjets || !mt_singlet || !mt_diboson || !mt_ttV ){
      std::cerr << "ERROR: one or more histograms are not initialized" << std::endl;
      exit(1);
    }
    
    //    drawhists.at(0)->Write(name+"_data");
    mt_ttsignal->Write(name+"_ttbarsignal");
    mt_ttother->Write(name+"_ttbarother");
    mt_singlet->Write(name+"_singlet");
    mt_ttV->Write(name+"_ttV");
    mt_zjets->Write(name+"_zjets");
    mt_diboson->Write(name+"_diboson");
    mt_wjets->Write(name+"_wjets");

    delete mt_ttsignal;
    delete mt_ttother;
    delete mt_zjets;
    delete mt_wjets;
    delete mt_singlet;
    delete mt_diboson;
    delete mt_ttV;

    //

    out_root.Close();

    for (TH1* h : drawhists) delete h;
}


double PlotterCP::CalcXSec(std::vector<TString> datasetVec, std::vector<TString> pseudoDatasetVec, double InclusiveXsectionVec[4], double InclusiveXsectionStatErrorVec[4], TString Systematic, TString Shift)
{
    double NrOfEvts_VisGen_afterSelection_noweight = 0, NrOfEvts_VisGen_afterSelection = 0;
    double NrOfEvts_afterSelection_noweight = 0, NrOfEvts_afterSelection = 0;
    double NrOfEvts_Gen_afterRecoSelection_noweight = 0, NrOfEvts_Gen_afterRecoSelection = 0;
    double NrOfEvts = 0;

    TH1D *numhists[hists.size()];
    double numbers[5] = {0., 0., 0., 0., 0.};//[0]=data, [1]=Signal, [2]Signal(only lumi & PU weights), [3]ttbar background, [4]non-ttbar background
    double error_numbers[5] = {0., 0., 0., 0., 0.};//square of error: [0]=data, [1]=Signal, [2]Signal(only lumi & PU weights), [3]ttbar background, [4]non-ttbar background

    if (Systematic.Contains("UP"))  Shift = "Up";
    if (Systematic.Contains("DOWN")) Shift = "Down";

    TString histoForInclXSec;
    if (usePoissonSmearedPseudoData) histoForInclXSec = "basic_jet_multiplicity_step8";// you may use other histograms, e.g.: "events_weighted_step8" or "HypjetMultiXSec" or even "..._step7"
    else histoForInclXSec = "events_weighted_step8";//check FIXME below (in 5 lines), before changing it

    TString Channel;
    if (channelType == 0) Channel = "ee";
    if (channelType == 1) Channel = "mumu";
    if (channelType == 2) Channel = "emu";
    if (channelType == 3) Channel = "combined";

    for (unsigned int i = 0; i < datasetVec.size(); i++) {
        TH1D *hist = fileReader->GetClone<TH1D>(datasetVec[i], histoForInclXSec);
	double LumiWeight;
	double year_corr_rescale_from_file = 1.;
	if (year == "fullRun2UL" || year == "2016") {
	  year_corr_rescale_from_file = configHelper->get_year_corr_rescale_from_file(Systematic, datasetVec.at(i));
	  configHelper->get_year_from_file(datasetVec.at(i));
	  LumiWeight = configHelper->CalcLumiWeight(datasetVec.at(i), configHelper->lumi_year);
	  configHelper->ApplyFlatWeights(hist, LumiWeight);
	  if (year_corr_rescale_from_file != 1.) {
	    TString datasetNominal = configHelper->getCorrespondingNominalDataset(datasetVec[i], Systematic, Channel);
	    TH1D* histNominal = fileReader->GetClone<TH1D>(datasetNominal, histoForInclXSec);
	    if (!histNominal) return false;
	    double LumiWeightNominal = configHelper->CalcLumiWeight(datasetNominal, configHelper->lumi_year);
	    configHelper->ApplyFlatWeights(histNominal, LumiWeightNominal);
	    configHelper->RescaleHistYearToYearCorr(hist, histNominal, year_corr_rescale_from_file);
	  }
	}
	else
	  configHelper->ApplyFlatWeights(hist, configHelper->CalcLumiWeight(datasetVec.at(i)));
        removeNegativeExpectations(hist);//FIXME: fine to do this as long as a one-bin histograms is used (like "events_weighted"), otherwise the integral must be operated
        numhists[i] = hist;
    }

    if (usePoissonSmearedPseudoData) {
        TH1D *numhistsForPseudoData[histsForPseudoData.size()];

        for (unsigned int i = 0; i < pseudoDatasetVec.size(); i++) {
            TH1D *histForPseudoData = fileReader->GetClone<TH1D>(pseudoDatasetVec[i], histoForInclXSec);
	    double LumiWeightPseudoData;
	    double year_corr_rescale_from_file_pseudo = 1.;
	    if (year == "fullRun2UL" || year == "2016") {
	      year_corr_rescale_from_file_pseudo = configHelper->get_year_corr_rescale_from_file(Systematic, pseudoDatasetVec.at(i));
	      configHelper->get_year_from_file(pseudoDatasetVec.at(i));
	      LumiWeightPseudoData = configHelper->CalcLumiWeight(pseudoDatasetVec.at(i), configHelper->lumi_year);
	      configHelper->ApplyFlatWeights(histForPseudoData, LumiWeightPseudoData);
	      if (year_corr_rescale_from_file_pseudo != 1.) {
		TString pseudoDatasetNominal = configHelper->getCorrespondingNominalDataset(pseudoDatasetVec.at(i), Systematic, Channel);
		TH1D* histForPseudoDataNominal = fileReader->GetClone<TH1D>(pseudoDatasetNominal, histoForInclXSec);
		if (!histForPseudoDataNominal) return false;
		double LumiWeightPseudoDataNominal = configHelper->CalcLumiWeight(pseudoDatasetNominal, configHelper->lumi_year);
		configHelper->ApplyFlatWeights(histForPseudoDataNominal, LumiWeightPseudoDataNominal);
		configHelper->RescaleHistYearToYearCorr(histForPseudoData, histForPseudoDataNominal, year_corr_rescale_from_file_pseudo);
	      }
	    }
	    else
	      configHelper->ApplyFlatWeights(histForPseudoData, configHelper->CalcLumiWeight(pseudoDatasetVec.at(i)));
            removeNegativeExpectations(histForPseudoData);
            numhistsForPseudoData[i] = histForPseudoData;
        }

        if (channelType != 3) {// for ee, emu or mumu channel
            for (size_t i = 1; i < pseudoDatasetVec.size(); ++i) {
	      if (!pseudoDatasetVec.at(i).Contains("ttbarsignalplustau") && !pseudoDatasetVec.at(i).Contains("ttbarsignalviatau") && !pseudoDatasetVec.at(i).Contains("ttbarbgviatau") && !pseudoDatasetVec.at(i).Contains("run")) {
                    numhists[0]->Add(numhistsForPseudoData[i]);
                }
            }
            PlotterCP::performPoissonSmearing(numhists[0]);
        } else if (channelType == 3) {// for combined channel
            //Indeces must match to the order of data files in filelist {syst}_combined.txt
            TString chToInclude[3] = {"ee_", "emu_", "mumu_"};
            for (int ch_it = 0; ch_it < 3; ++ch_it) {
                for (size_t i = ch_it + 1; i < pseudoDatasetVec.size(); ++i) {//skipping pseudo-data files (moving iterator by +1 when initialized) to which MC is added
		  if (!pseudoDatasetVec.at(i).Contains("ttbarsignalplustau") && !pseudoDatasetVec.at(i).Contains("ttbarsignalviatau") && !pseudoDatasetVec.at(i).Contains("ttbarbgviatau") &&
                        !pseudoDatasetVec.at(i).Contains("run") && pseudoDatasetVec.at(i).Contains(chToInclude[ch_it]))
                    {
                        numhists[ch_it]->Add(numhistsForPseudoData[i]);
                    }
                }
                PlotterCP::performPoissonSmearing(numhists[ch_it]);
            }
        }
    }

    for (unsigned int i = 0; i < hists.size(); i++) {
        if (legends.at(i) == "Data") {
            if (usePoissonSmearedPseudoData) {
                numbers[0] += numhists[i]->Integral(0, numhists[i]->GetNbinsX()+1);//including under/overflow bins
                for (int j = 0; j <= numhists[i]->GetNbinsX()+1; ++j) {
                    error_numbers[0] += TMath::Sqrt(numhists[i]->GetBinContent(j)) * TMath::Sqrt(numhists[i]->GetBinContent(j));//"Data"-like approach
                }
            } else {
                numbers[0] += numhists[i]->Integral();
                error_numbers[0] += numhists[i]->GetBinError(2) * numhists[i]->GetBinError(2);//this bin selection is hardcoded; please change it if changes when filling in Analysis.C
            }
        } else if (legends.at(i) == "t#bar{t} signal") {
            if (usePoissonSmearedPseudoData) {
                numbers[1] += numhists[i]->Integral(0, numhists[i]->GetNbinsX()+1);//including under/overflow bins
                for (int j = 0; j <= numhists[i]->GetNbinsX()+1; ++j) {
                    error_numbers[1] += numhists[i]->GetBinError(j) * numhists[i]->GetBinError(j);
                }
            } else {
                numbers[1] += numhists[i]->Integral();
                error_numbers[1] += numhists[i]->GetBinError(2) * numhists[i]->GetBinError(2);//this bin selection is hardcoded; please change it if changes when filling in Analysis.C
            }

            TH1D *GenPlot = fileReader->GetClone<TH1D>(datasetVec.at(i), "GenAll");
            TH1D *GenPlot_noweight = fileReader->GetClone<TH1D>(datasetVec.at(i), "GenAll_noweight");
            TH1D *VisGenPlot = fileReader->GetClone<TH1D>(datasetVec.at(i), "VisGenAll");
            TH1D *VisGenPlot_noweight = fileReader->GetClone<TH1D>(datasetVec.at(i), "VisGenAll_noweight");
            TH1D *RecoGenPlot = fileReader->GetClone<TH1D>(datasetVec.at(i), "GenAll_RecoCuts");
            TH1D *RecoGenPlot_noweight = fileReader->GetClone<TH1D>(datasetVec.at(i), "GenAll_RecoCuts_noweight");
            TH1 *h_NrOfEvts = fileReader->GetClone<TH1>(datasetVec.at(i), "weightedEvents");

	    double LumiWeight;
	    double year_corr_rescale_from_file = 1.;
	    if (year == "fullRun2UL" || year == "2016") {
	      year_corr_rescale_from_file = configHelper->get_year_corr_rescale_from_file(Systematic, datasetVec.at(i));
	      configHelper->get_year_from_file(datasetVec.at(i));
	      LumiWeight = configHelper->CalcLumiWeight(datasetVec.at(i), configHelper->lumi_year);
	    }
	    else
	      LumiWeight = configHelper->CalcLumiWeight(datasetVec.at(i));

	    configHelper->ApplyFlatWeights(GenPlot, LumiWeight);
	    configHelper->ApplyFlatWeights(GenPlot_noweight, LumiWeight);
	    configHelper->ApplyFlatWeights(VisGenPlot, LumiWeight);
	    configHelper->ApplyFlatWeights(VisGenPlot_noweight, LumiWeight);
	    configHelper->ApplyFlatWeights(RecoGenPlot, LumiWeight);
	    configHelper->ApplyFlatWeights(RecoGenPlot_noweight, LumiWeight);
	    configHelper->ApplyFlatWeights(h_NrOfEvts, LumiWeight);

	    if (year_corr_rescale_from_file != 1.) {
	      TString datasetNominal = configHelper->getCorrespondingNominalDataset(datasetVec.at(i), Systematic, Channel);

	      TH1D *GenPlot_nominal = fileReader->GetClone<TH1D>(datasetNominal, "GenAll");
	      TH1D *GenPlot_noweight_nominal = fileReader->GetClone<TH1D>(datasetNominal, "GenAll_noweight");
	      TH1D *VisGenPlot_nominal = fileReader->GetClone<TH1D>(datasetNominal, "VisGenAll");
	      TH1D *VisGenPlot_noweight_nominal = fileReader->GetClone<TH1D>(datasetNominal, "VisGenAll_noweight");
	      TH1D *RecoGenPlot_nominal = fileReader->GetClone<TH1D>(datasetNominal, "GenAll_RecoCuts");
	      TH1D *RecoGenPlot_noweight_nominal = fileReader->GetClone<TH1D>(datasetNominal, "GenAll_RecoCuts_noweight");
	      TH1 *h_NrOfEvts_nominal = fileReader->GetClone<TH1>(datasetNominal, "weightedEvents");

	      if (!GenPlot_nominal) return false;
	      if (!GenPlot_noweight_nominal) return false;
	      if (!VisGenPlot_nominal) return false;
	      if (!VisGenPlot_noweight_nominal) return false;
	      if (!RecoGenPlot_nominal) return false;
	      if (!RecoGenPlot_noweight_nominal) return false;
	      if (!h_NrOfEvts_nominal) return false;

	      double LumiWeightNominal = configHelper->CalcLumiWeight(datasetNominal, configHelper->lumi_year);

	      configHelper->ApplyFlatWeights(GenPlot_nominal, LumiWeightNominal);
	      configHelper->ApplyFlatWeights(GenPlot_noweight_nominal, LumiWeightNominal);
	      configHelper->ApplyFlatWeights(VisGenPlot_nominal, LumiWeightNominal);
	      configHelper->ApplyFlatWeights(VisGenPlot_noweight_nominal, LumiWeightNominal);
	      configHelper->ApplyFlatWeights(RecoGenPlot_nominal, LumiWeightNominal);
	      configHelper->ApplyFlatWeights(RecoGenPlot_noweight_nominal, LumiWeightNominal);
	      configHelper->ApplyFlatWeights(h_NrOfEvts_nominal, LumiWeightNominal);

	      configHelper->RescaleHistYearToYearCorr(GenPlot, GenPlot_nominal, year_corr_rescale_from_file);
	      configHelper->RescaleHistYearToYearCorr(GenPlot_noweight, GenPlot_noweight_nominal, year_corr_rescale_from_file);
	      configHelper->RescaleHistYearToYearCorr(VisGenPlot, VisGenPlot_nominal, year_corr_rescale_from_file);
	      configHelper->RescaleHistYearToYearCorr(VisGenPlot_noweight, VisGenPlot_noweight_nominal, year_corr_rescale_from_file);
	      configHelper->RescaleHistYearToYearCorr(RecoGenPlot, RecoGenPlot_nominal, year_corr_rescale_from_file);
	      configHelper->RescaleHistYearToYearCorr(RecoGenPlot_noweight, RecoGenPlot_noweight_nominal, year_corr_rescale_from_file);
	      configHelper->RescaleHistYearToYearCorr(h_NrOfEvts, h_NrOfEvts_nominal, year_corr_rescale_from_file);
	    }

            NrOfEvts = h_NrOfEvts->GetBinContent(1);
            NrOfEvts_afterSelection += GenPlot->Integral();
            NrOfEvts_afterSelection_noweight += GenPlot_noweight->Integral();
            NrOfEvts_VisGen_afterSelection += VisGenPlot->Integral();
            NrOfEvts_Gen_afterRecoSelection += RecoGenPlot->Integral();
            NrOfEvts_Gen_afterRecoSelection_noweight += RecoGenPlot_noweight->Integral();
            NrOfEvts_VisGen_afterSelection_noweight += VisGenPlot_noweight->Integral();

            numbers[2] += GenPlot->Integral();
            error_numbers[2] += GenPlot->GetBinError(18) * GenPlot->GetBinError(18);//this bin selection is hardcoded; please change it if changes when filling in Analysis.C
        } else if (legends.at(i) == "t#bar{t} other") {
            if (usePoissonSmearedPseudoData) {
                numbers[3] += numhists[i]->Integral(0, numhists[i]->GetNbinsX()+1);//including under/overflow bins
                for (int j = 0; j <= numhists[i]->GetNbinsX()+1; ++j) {
                    error_numbers[3] += numhists[i]->GetBinError(j) * numhists[i]->GetBinError(j);
                }
            } else {
                numbers[3] += numhists[i]->Integral();
                error_numbers[3] += numhists[i]->GetBinError(2) * numhists[i]->GetBinError(2);//this bin selection is hardcoded; please change it if changes when filling in Analysis.C
            }
        } else {
            if ((legends.at(i) == DYEntry)) numhists[i]->Scale(DYScale.at(channelType));
            if ((legends.at(i) == DYEntry) && Systematic.Contains("DY_") && Shift == "Up") numhists[i]->Scale(1.3);
            if ((legends.at(i) == DYEntry) && Systematic.Contains("DY_") && Shift == "Down") numhists[i]->Scale(0.7);
            if (Systematic.Contains("BG_") && Shift=="Up" && legends.at(i)!= "t#bar{t} other" && legends.at(i) != DYEntry) numhists[i]->Scale(1.3);
            if (Systematic.Contains("BG_") && Shift=="Down" && legends.at(i)!= "t#bar{t} other" && legends.at(i) != DYEntry) numhists[i]->Scale(0.7);

            if (usePoissonSmearedPseudoData) {
                numbers[4] += numhists[i]->Integral(0, numhists[i]->GetNbinsX()+1);//including under/overflow bins
                for (int j = 0; j <= numhists[i]->GetNbinsX()+1; ++j) {
                    error_numbers[4] += numhists[i]->GetBinError(j) * numhists[i]->GetBinError(j);
                }
            } else {
                numbers[4] += numhists[i]->Integral();
                error_numbers[4] += numhists[i]->GetBinError(2) * numhists[i]->GetBinError(2);//this bin selection is hardcoded; please change it if changes when filling in Analysis.C
            }
        }
    }
    ////////////////////////////Make output for tables
    double tmp_num = 0;

    ofstream EventFile, XSecFile;
    TString outdir = utils::assignFolder(this->outpathPlots, subfolderChannel.Copy().Remove(0,1), Systematic);
    EventFile.open(outdir.Copy() + "Events.txt");
    XSecFile.open(outdir.Copy() + "InclusiveXSec.txt");

    double bg_num = 0;
    for (unsigned int i = 0; i < hists.size(); i++) {
        if (usePoissonSmearedPseudoData) {
            tmp_num += numhists[i]->Integral(0,numhists[i]->GetNbinsX()+1);
        } else {
            tmp_num += numhists[i]->Integral();
        }

        if (i == (hists.size()-1)) {
            EventFile << legends.at(i) << ": " << tmp_num << std::endl;
            bg_num += tmp_num;
            tmp_num = 0;
        } else if (legends.at(i) != legends.at(i+1)) {
            EventFile << legends.at(i) << ": " << tmp_num << std::endl;
            if (legends.at(i) != "Data") bg_num += tmp_num;
            tmp_num = 0;
        }
    }
    EventFile << "Total MCs: "      << bg_num << std::endl;
    EventFile << "\nDataEvents= "   << numbers[0] << std::endl;
    EventFile << "SignalReco= "     << numbers[1] << std::endl;
    EventFile << "SignalGen = "     << numbers[2] << std::endl;
    EventFile << "ttbar bags= "     << numbers[3] << std::endl;
    EventFile << "All Backgd= "     << numbers[4] << std::endl;
    EventFile << "Efficiency= "     << (numbers[1] / numbers[2]) << std::endl;
    EventFile << "BrancRatio= "     << bRatio[channelType] << std::endl;
    EventFile << "Total Gen Events (no weights)= "                      << NrOfEvts << std::endl;
    EventFile << "Gen Events after Selection (no weights)= "            << NrOfEvts_afterSelection_noweight << std::endl;
    EventFile << "Visible Gen Events after Selection (no weights)= "    << NrOfEvts_VisGen_afterSelection_noweight << std::endl;
    EventFile << "Acceptance= "                                         << NrOfEvts_afterSelection_noweight / NrOfEvts << std::endl;
    EventFile << "Visible Acceptance= "                                 << NrOfEvts_VisGen_afterSelection_noweight / NrOfEvts << std::endl;
    EventFile << "------------------------------------------------------------------------------" << std::endl;
    EventFile << "Efficiency and acceptance definitions as proposed by the TopXSection conveners\n" << std::endl;
    EventFile << "N_rec = "                                         << numbers[1] << std::endl;
    EventFile << "N_gen (with cuts at parton level) = "             << NrOfEvts_VisGen_afterSelection << std::endl;
    EventFile << "N_gen (with cuts at parton level, no weights) = " << NrOfEvts_VisGen_afterSelection_noweight << std::endl;
    EventFile << "N_gen (with cuts at reco level) = "               << NrOfEvts_Gen_afterRecoSelection << std::endl;
    EventFile << "N_gen (with cuts at reco level, no weights) = "   << NrOfEvts_Gen_afterRecoSelection_noweight << std::endl;
    EventFile << "N_gen = "                                         << numbers[2] << std::endl;
    EventFile << "\nEfficiency = N_rec / N_gen (with cuts at parton level) = "              << numbers[1] / NrOfEvts_VisGen_afterSelection << std::endl;
    EventFile << "Efficiency = N_rec / N_gen (with cuts at parton level && noweights) = "   << numbers[1] / NrOfEvts_VisGen_afterSelection_noweight << std::endl;
    EventFile << "\nEfficiency' = N_rec / N_gen (with cuts at reco level) = "               << numbers[1] / NrOfEvts_Gen_afterRecoSelection << std::endl;
    EventFile << "Efficiency' = N_rec / N_gen (with cuts at reco level && noweights) = "    << numbers[1] / NrOfEvts_Gen_afterRecoSelection_noweight << std::endl;
    EventFile << "\nAcceptance = N_gen (with cuts at parton level) / N_gen = "              << NrOfEvts_VisGen_afterSelection/numbers[2] << std::endl;
    EventFile << "Acceptance = N_gen (with cuts at parton level && noweights) / N_gen = "   << NrOfEvts_VisGen_afterSelection_noweight/numbers[2] << std::endl;
    EventFile << "Eff * Acc = " << numbers[1] / numbers[2] << std::endl;
    EventFile << "------------------------------------------------------------------------------" << std::endl;

    // Acceptance driven correction needed to avoiding dependence from N^gen_dilepton/N^gen_all(inclusive) ratio specific for each MC sample dataset when using pseudo-data
    double correctionForInclXSecWhilePseudoData = 1.;
    if (usePoissonSmearedPseudoData) correctionForInclXSecWhilePseudoData = bRatio[channelType] / (NrOfEvts_afterSelection_noweight / NrOfEvts);

    double xsec = ((numbers[0] - numbers[4]) * (numbers[1] / (numbers[1] + numbers[3]))) / ((numbers[1] / numbers[2]) * bRatio[channelType] * lumi);
    double xsecstaterror = TMath::Sqrt(error_numbers[0]) * (numbers[1] / (numbers[1] + numbers[3])) / ((numbers[1] / numbers[2]) * bRatio[channelType] * lumi);

    if (usePoissonSmearedPseudoData) {
        xsec = correctionForInclXSecWhilePseudoData * xsec;
        xsecstaterror = correctionForInclXSecWhilePseudoData * xsecstaterror;
    }

    if (channelType != 3) {
        InclusiveXsectionVec[channelType] = xsec;
        InclusiveXsectionStatErrorVec[channelType] = xsecstaterror;
    } else {
        TString eefilename = this->outpathPlots + "/" + Systematic + "/ee/InclusiveXSec.txt";
        TString mumufilename = this->outpathPlots + "/" + Systematic + "/mumu/InclusiveXSec.txt";
        TString emufilename = this->outpathPlots + "/" + Systematic + "/emu/InclusiveXSec.txt";

        // Check the existence of the file
        if (gSystem->AccessPathName(eefilename) || gSystem->AccessPathName(emufilename) || gSystem->AccessPathName(mumufilename)) {
            std::cout << "WARNING (in CalcXSec)!!" << std::endl;
            std::cout << "One of the input files you use for the combined XSection measurement doesn't exist!!\nExiting!!" << std::endl;
            exit(888);
        }

        TString Dummy="";
        ifstream ResultsEE(eefilename);
        ifstream ResultsEMu(emufilename);
        ifstream ResultsMuMu(mumufilename);
        ResultsEE   >> Dummy >> Dummy >> Dummy >> Dummy >> Dummy >> InclusiveXsectionVec[0] >> Dummy >> InclusiveXsectionStatErrorVec[0];
        ResultsMuMu >> Dummy >> Dummy >> Dummy >> Dummy >> Dummy >> InclusiveXsectionVec[1] >> Dummy >> InclusiveXsectionStatErrorVec[1];
        ResultsEMu  >> Dummy >> Dummy >> Dummy >> Dummy >> Dummy >> InclusiveXsectionVec[2] >> Dummy >> InclusiveXsectionStatErrorVec[2];

        ResultsEE.close(); ResultsEMu.close(); ResultsMuMu.close();

        InclusiveXsectionVec[channelType] = (InclusiveXsectionVec[0] / (InclusiveXsectionStatErrorVec[0] * InclusiveXsectionStatErrorVec[0]) +
                                             InclusiveXsectionVec[1] / (InclusiveXsectionStatErrorVec[1] * InclusiveXsectionStatErrorVec[1]) +
                                             InclusiveXsectionVec[2] / (InclusiveXsectionStatErrorVec[2] * InclusiveXsectionStatErrorVec[2])
                                            ) /
                                            (1 / (InclusiveXsectionStatErrorVec[0] * InclusiveXsectionStatErrorVec[0]) +
                                             1 / (InclusiveXsectionStatErrorVec[1] * InclusiveXsectionStatErrorVec[1]) +
                                             1 / (InclusiveXsectionStatErrorVec[2] * InclusiveXsectionStatErrorVec[2])
                                            );

        InclusiveXsectionStatErrorVec[channelType] = 1 / (TMath::Sqrt(
                                            (1 / (InclusiveXsectionStatErrorVec[0] * InclusiveXsectionStatErrorVec[0])) +
                                            (1 / (InclusiveXsectionStatErrorVec[1] * InclusiveXsectionStatErrorVec[1])) +
                                            (1 / (InclusiveXsectionStatErrorVec[2] * InclusiveXsectionStatErrorVec[2]))
                                                                     ));
    }

    EventFile << "XSection  = " << InclusiveXsectionVec[channelType] << std::endl;
    EventFile << "XSecStaErr= " << InclusiveXsectionStatErrorVec[channelType] << std::endl;
    EventFile.close();
    XSecFile << "Systematic: " << Systematic << " Channel: " << subfolderChannel << " InclXSection: "
        << InclusiveXsectionVec[channelType] << " AbsStatError: " << InclusiveXsectionStatErrorVec[channelType] << std::endl;
    XSecFile.close();
    std::cout << "\nInclusive XSection information saved in: " << outdir << std::endl;
    return xsec;
}


void PlotterCP::PlotXSec(TString Channel)
{
    TH1::AddDirectory(kFALSE);
    //Using systematics defined in PlotterConfigurationHelper
    std::vector<TString> syst = this->uncertaintyBandSystematics;

    std::vector<TString> vec_systematic;
    for (int i = 0; i < (int) syst.size(); i++) {
        if (configHelper->IsNoUpDown(syst.at(i))) vec_systematic.push_back(syst.at(i));
        else vec_systematic.push_back(syst.at(i) + "_");
    }

    //Uncomment this for customized systematics
    /*std::vector<TString> vec_systematic {
          "MASS_", "MESCALE_", "MEFACSCALE_", "MERENSCALE_", "PSISRSCALE_2_", "PSFSRSCALE_2_", "BSEMILEP_", "BFRAG_", "BFRAG_PETERSON_",
       // "ERDON_", "ERDONRETUNE_", "GLUONMOVETUNE_", "UETUNE_", "MATCH_", "PDF_ALPHAS_", "PDF_",
          "UETUNE_", "MATCH_", "PDF_ALPHAS_", "PDF_",
          "PDF_ALPHAS_", "PDF_",
          "UNCLUSTERED_", "BTAG_", "BTAG_LJET_", "KIN_", "LEPT_", "L1PREFIRING_", "TRIG_", "BG_", "DY_", "PU_", "JER_",
          "JESAbsoluteStat_", "JESAbsoluteScale_", "JESAbsoluteMPFBias_", "JESFragmentation_", "JESSinglePionECAL_", "JESSinglePionHCAL_",
          "JESFlavorQCD_", "JESTimePtEta_", "JESRelativeBal_", "JESRelativeJEREC1_", "JESRelativePtBB_", "JESRelativePtEC1_", "JESRelativeFSR_",
          "JESRelativeStatFSR_", "JESRelativeStatEC_", "JESPileUpDataMC_", "JESPileUpPtRef_", "JESPileUpPtEC1_", "JESPileUpPtBB_"
        //Not needed: , "MOD_", "HAD_", "JES_", "JESAbsoluteFlavMap_", "JESRelativeJEREC2_", "JESRelativeJERHF_", "JESRelativePtEC2_",
        //"JESRelativePtHF_", "JESRelativeStatHF_", "JESPileUpPtEC2_", "JESPileUpPtHF_",

    };*/

    std::vector<TString> vec_channel {"ee", "mumu", "emu", "combined"};

    double BR_Error = bRatioError;
    double Lumi_Error = lumiError;

    double InclusiveXsectionPlot[4] = {0.}, InclusiveXsectionStatErrorPlot[4] = {0.}, InclusiveXsectionSysErrorPlot[4] = {0.}, InclusiveXsectionTotalErrorPlot[4] = {0.};
    std::vector<TString> syst_files_not_found;

    for (int j = 0; j < (int)vec_channel.size(); j++) {
        TString outdir = utils::assignFolder(this->outpathPlots, vec_channel.at(j), TString("FinalResults"));
        ifstream SysResultsList(this->outpathPlots + "/Nominal/" + vec_channel.at(j) + "/InclusiveXSec.txt");
        TString DUMMY;
        SysResultsList >> DUMMY >> DUMMY >> DUMMY >> DUMMY >> DUMMY >> InclusiveXsectionPlot[j] >> DUMMY >> InclusiveXsectionStatErrorPlot[j];
        SysResultsList.close();

        std::ofstream OutputFile(outdir.Copy() + "InclusiveXSecResultLateX.txt", std::ofstream::trunc);
        OutputFile << "Inclusive XSection Numerical Results for channel " << vec_channel.at(j) << std::endl;

        //FIXME: propagate handle for envelopes
        double syst_square_for_channel = 0.0;
        double tot_scale_variation = 0.0; // the 'envelope' (in this case out of integrated bin) type treatment is applied
        double tot_bfrag_variation = 0.0; // the 'envelope' (in this case out of integrated bin) type treatment is applied
        double tot_colorrec_variation = 0.0; // the 'envelope' (in this case out of integrated bin) type treatment is applied
        double tot_jes_variation_squared = 0.0; // the quadratic sum of all 27 sources (in this case out of integrated bin) type treatment is applied
	std::map<TString, double> tot_year_to_year_corr_variations_squared; // the quadratic sum of all year-to-year correlated variations
	double tot_jer_variation_squared = 0.0; // the quadratic sum of all JER sources
	double tot_lept_variation_squared = 0.0; // the quadratic sum of all LEPT sources

        for (int i = 0; i < (int) vec_systematic.size(); i++) {
            ifstream SysUP, SysDOWN;

            if (vec_systematic.at(i) == "HAD_") {
                SysUP.open(this->outpathPlots + "/MCATNLO/" + vec_channel.at(j) + "/InclusiveXSec.txt");
                SysDOWN.open(this->outpathPlots + "/MCATNLO/" + vec_channel.at(j) + "/InclusiveXSec.txt");
            } else if (vec_systematic.at(i) == "MOD_") {
                SysUP.open(this->outpathPlots + "/POWHEG/" + vec_channel.at(j) + "/InclusiveXSec.txt");
                SysDOWN.open(this->outpathPlots + "/POWHEG/" + vec_channel.at(j) + "/InclusiveXSec.txt");
            } else if (configHelper->IsNoUpDown(vec_systematic.at(i))) {
                SysUP.open(this->outpathPlots + "/" + vec_systematic.at(i) +"/" + vec_channel.at(j) + "/InclusiveXSec.txt");
                SysDOWN.open(this->outpathPlots + "/" + vec_systematic.at(i) + "/" + vec_channel.at(j) + "/InclusiveXSec.txt");
            } else if (vec_systematic.at(i) == "POWHEG") {
                continue;
            }
            else {
                SysUP.open(this->outpathPlots + "/" + vec_systematic.at(i) + "UP/" + vec_channel.at(j) + "/InclusiveXSec.txt");
                SysDOWN.open(this->outpathPlots + "/" + vec_systematic.at(i) + "DOWN/" + vec_channel.at(j) + "/InclusiveXSec.txt");
                //if (!SysUP.is_open() || !SysDOWN.is_open()) continue;
            }

            if (!SysUP.is_open()) {
                syst_files_not_found.push_back(this->outpathPlots + "/" + vec_systematic.at(i) + "UP/" + vec_channel.at(j) + "/InclusiveXSec.txt");
                continue;
            }
            if (!SysDOWN.is_open()) {
                syst_files_not_found.push_back(this->outpathPlots + "/" + vec_systematic.at(i) + "DOWN/" + vec_channel.at(j) + "/InclusiveXSec.txt");
                continue;
            }

            double VarUp = 0, VarDown = 0, StatErrUp = 0, StatErrDown = 0;

            SysUP >> DUMMY >> DUMMY >> DUMMY >> DUMMY >> DUMMY >> VarUp >> DUMMY >> StatErrUp;
            SysDOWN >> DUMMY >> DUMMY >> DUMMY >> DUMMY >> DUMMY >> VarDown >> DUMMY >> StatErrDown;
            SysUP.close();
            SysDOWN.close();

            // Systematic error in %
            double sys_err = (TMath::Abs(InclusiveXsectionPlot[j] - VarUp) + TMath::Abs(InclusiveXsectionPlot[j] - VarDown)) * 0.5 / InclusiveXsectionPlot[j];
            if (vec_systematic.at(i).Contains("MASS")) sys_err = sys_err / massUncReductionFactor;
            // if (vec_systematic.at(i).Contains("PSFSRSCALE")) sys_err = sys_err / fsrScaleUncReductionFactor;
            if (vec_systematic.at(i) == "MESCALE_" || vec_systematic.at(i) == "MEFACSCALE_" || vec_systematic.at(i) == "MERENSCALE_") {
               if (sys_err > tot_scale_variation) tot_scale_variation = sys_err;
               OutputFile << "==> " << vec_systematic.at(i) << " not propagated to total error directly (but included to related envelope and thus to total error) (%): "
                    << setprecision(3) << sys_err * 100 << std::endl;
               continue;
            }
            if (vec_systematic.at(i) == "BFRAG_" || vec_systematic.at(i) == "BFRAG_PETERSON" || vec_systematic.at(i) == "BFRAG_CENTRAL") {
               if (sys_err > tot_bfrag_variation) tot_bfrag_variation = sys_err;
               OutputFile << "==> " << vec_systematic.at(i) << " not propagated to total error directly (but included to related envelope and thus to total error) (%): "
                    << setprecision(3) << sys_err * 100 << std::endl;
               continue;
            }
            if (vec_systematic.at(i) == "ERDON" || vec_systematic.at(i) == "ERDONRETUNE" || vec_systematic.at(i) == "GLUONMOVETUNE") {
               if (sys_err > tot_colorrec_variation) tot_colorrec_variation = sys_err;
               OutputFile << "==> " << vec_systematic.at(i) << " not propagated to total error directly (but included to related envelope and thus to total error) (%): "
                    << setprecision(3) << sys_err * 100 << std::endl;
               continue;
            }
            if (vec_systematic.at(i) == "HAD_" || vec_systematic.at(i) == "MOD_") {
               OutputFile << "==> " << vec_systematic.at(i) << " not propagated to total error (%): " << setprecision(3) << sys_err * 100 << std::endl;
               continue;
            }

            if (vec_systematic.at(i) == "JES_") {
               if (this->configHelper->gSysJES == 2){
                   OutputFile << "==> " << vec_systematic.at(i) << " the cumulative jes variation is not propagated to total error, quoted here for the reference (%): "
                              << setprecision(3) << sys_err * 100 << std::endl;
                    continue;
               }
               else if (this->configHelper->gSysJES == 1){
                   OutputFile << "==> " << vec_systematic.at(i) <<
                    " total jes source, quoted here for the reference (%): "
                        << setprecision(3) << sys_err * 100 << std::endl;
                   tot_jes_variation_squared += sys_err * sys_err;
                    continue;
               }
            }
            if (vec_systematic.at(i).Contains("JES") && !(vec_systematic.at(i) == "JES_") && this->configHelper->gSysJES == 2) {
               OutputFile << "==> " << vec_systematic.at(i) <<
                " individual jes source propagated to total jes error via quadratic sum (and thus to total error), quoted here for the reference (%): "
                    << setprecision(3) << sys_err * 100 << std::endl;
               tot_jes_variation_squared += sys_err * sys_err;
               continue;
            }


	    if(vec_systematic.at(i).Contains("JER"))
	      tot_jer_variation_squared += sys_err * sys_err;

	    if(vec_systematic.at(i).Contains("ELE") || vec_systematic.at(i).Contains("MUON"))
	      tot_lept_variation_squared += sys_err * sys_err;

	    TString temp_sys = vec_systematic.at(i);
	    const char* s = temp_sys.Data();
	    const int n = temp_sys.Length();
	    if(TString(s + n - 1, 1) == "_")
	      temp_sys = TString(s, n - 1);
	    if (PlotterConfigurationHelper::IsYearToYearCorr(year, yearCorr, temp_sys) && yearCorr != "") {
	       PlotterConfigurationHelper::RemoveYearCorr(year, yearCorr, temp_sys);
	       if (!configHelper->IsNoUpDown(temp_sys))
		 temp_sys = temp_sys + "_";
	       if (tot_year_to_year_corr_variations_squared.find(temp_sys) != tot_year_to_year_corr_variations_squared.end())
		 tot_year_to_year_corr_variations_squared[temp_sys] += sys_err * sys_err;
	       else
		 tot_year_to_year_corr_variations_squared[temp_sys] = sys_err * sys_err;
            }


            syst_square_for_channel += sys_err * sys_err;
            OutputFile << vec_systematic.at(i) << " (%): " << setprecision(3) << sys_err * 100 << std::endl;
        }

	if(yearCorr != "") {
	  OutputFile << "==> " << " total variations for systematics with correlations among years" << std::endl;
	  OutputFile << "***" << std::endl;
	  for (auto tot_var: tot_year_to_year_corr_variations_squared) {
	    OutputFile << tot_var.first << " (%): " << setprecision(3) << TMath::Sqrt(tot_var.second) * 100 << std::endl;
	  }
	  OutputFile << "***" << std::endl;
	}

	// Consider the relevant sources, which were not added yet
        syst_square_for_channel += tot_scale_variation * tot_scale_variation;
        syst_square_for_channel += tot_bfrag_variation * tot_bfrag_variation;
        syst_square_for_channel += tot_colorrec_variation * tot_colorrec_variation;
        syst_square_for_channel += tot_jes_variation_squared;

	OutputFile << "TOT_LEPT (%): "      << setprecision(3) << TMath::Sqrt(tot_lept_variation_squared) * 100 << std::endl;
	OutputFile << "TOT_JER (%): "       << setprecision(3) << TMath::Sqrt(tot_jer_variation_squared) * 100 << std::endl;
        OutputFile << "TOT_JES_ (%): "      << setprecision(3) << TMath::Sqrt(tot_jes_variation_squared) * 100 << std::endl;
        OutputFile << "TOT_SCALE_ (%): "    << setprecision(3) << tot_scale_variation * 100 << std::endl;
        OutputFile << "TOT_BFRAG_ (%): "    << setprecision(3) << tot_bfrag_variation * 100 << std::endl;
        OutputFile << "TOT_COLORREC_ (%): " << setprecision(3) << tot_colorrec_variation * 100 << std::endl;
        OutputFile << "BranchingRatio (%): "<< setprecision(3) << BR_Error * 100 << std::endl;
        OutputFile << "Luminosity (%): "    << setprecision(3) << Lumi_Error * 100 << std::endl;

        InclusiveXsectionSysErrorPlot[j] = TMath::Sqrt(syst_square_for_channel + BR_Error * BR_Error + Lumi_Error * Lumi_Error);

        InclusiveXsectionTotalErrorPlot[j] = sqrt(InclusiveXsectionStatErrorPlot[j] * InclusiveXsectionStatErrorPlot[j] +
                                                  InclusiveXsectionPlot[j] * InclusiveXsectionSysErrorPlot[j] * InclusiveXsectionPlot[j] * InclusiveXsectionSysErrorPlot[j]
                                                 );

        OutputFile << "\n\n*******************************************************************************\n\n";
        OutputFile << " InclXsec[pb]     Stat.[pb]    Syst.[pb]   Total[pb]"<<std::endl;
        OutputFile << setprecision(6) << InclusiveXsectionPlot[j] << " +- " << setprecision(3) << InclusiveXsectionStatErrorPlot[j]
            << " +- " << setprecision(4) << InclusiveXsectionSysErrorPlot[j] * InclusiveXsectionPlot[j]
                << " +- " << setprecision(4) << InclusiveXsectionTotalErrorPlot[j] << std::endl;
        OutputFile.close();
    }

    if (syst_files_not_found.size()>0){
        std::cout << "      WARNING: " << syst_files_not_found.size() << " file(s) not found for Inclusive Xsection plotting." << std::endl
                  << "          Files not found: ";
        for (auto file: syst_files_not_found) std::cout << file << " ";
        std::cout << endl << endl;
    }

    // measured results with statistical error
    Double_t mx[]   = {      0.50,       1.50,       2.50,       3.50};
    Double_t mexl[] = {      0.00,       0.00,       0.00,       0.00};
    Double_t mexh[] = {      0.00,       0.00,       0.00,       0.00};

    TGraphAsymmErrors *mplot = new TGraphAsymmErrors(4, mx, InclusiveXsectionPlot, mexl, mexh, InclusiveXsectionStatErrorPlot, InclusiveXsectionStatErrorPlot);
    mplot->SetMarkerStyle(20);
    mplot->GetYaxis()->SetNoExponent(kTRUE);
    mplot->SetMarkerColor(kBlack);
    mplot->SetMarkerSize(0.2);
    mplot->SetLineColor(kBlack);

    TGraphAsymmErrors *mplotwithsys = new TGraphAsymmErrors(4, mx, InclusiveXsectionPlot, mexl, mexh, InclusiveXsectionTotalErrorPlot, InclusiveXsectionTotalErrorPlot);
    mplotwithsys->SetMarkerStyle(20);
    mplotwithsys->SetMarkerColor(kBlack);
    mplotwithsys->SetMarkerSize(1.5);
    mplotwithsys->SetLineColor(kBlack);

    // top-16-005
    Double_t cmsemu15x[]   = {    2.65};
    Double_t cmsemu15y[]   = {    815.};
    Double_t cmsemu15exl[] = {      0.0};
    Double_t cmsemu15exh[] = {      0.0};
    Double_t cmsemu15eyl_stat[] = {    9.};
    Double_t cmsemu15eyh_stat[] = {    9.};
    Double_t cmsemu15eyl_stat_syst[] = {    43.43};
    Double_t cmsemu15eyh_stat_syst[] = {    43.43};

    TGraphAsymmErrors *cmsemu15plot = new TGraphAsymmErrors(6, cmsemu15x, cmsemu15y, cmsemu15exl, cmsemu15exh, cmsemu15eyl_stat, cmsemu15eyh_stat);
    cmsemu15plot->SetMarkerStyle(21);
    cmsemu15plot->GetYaxis()->SetNoExponent(kTRUE);
    cmsemu15plot->SetMarkerColor(kRed+1);
    cmsemu15plot->SetMarkerSize(1.5);
    cmsemu15plot->SetLineColor(kRed+1);

    TGraphAsymmErrors *cmsemu15plot_stat_syst = new TGraphAsymmErrors(6, cmsemu15x, cmsemu15y, cmsemu15exl, cmsemu15exh, cmsemu15eyl_stat_syst, cmsemu15eyh_stat_syst);
    cmsemu15plot_stat_syst->SetMarkerStyle(21);
    cmsemu15plot_stat_syst->GetYaxis()->SetNoExponent(kTRUE);
    cmsemu15plot_stat_syst->SetMarkerColor(kRed+1);
    cmsemu15plot_stat_syst->SetMarkerSize(1.5);
    cmsemu15plot_stat_syst->SetLineColor(kRed+1);

    // NNLO + NNLL : arXiv:1112.5675; arXiv:1303.6254.
    Double_t nnlomean = 831.76;
    Double_t errDown = 45.63; // scale + pdf + alpha_s
    Double_t errUp   = 40.25; // scale + pdf + alpha_s
    Double_t nnlox[]   = {    -0.5,     0.5,   1.5,     2.5,     3.5,     4.5};
    Double_t nnloy[]   = {nnlomean,nnlomean,nnlomean,nnlomean,nnlomean,nnlomean};
    Double_t nnloexl[] = {      .4,    .4,      .5,      .5,      .5,      .5};
    Double_t nnloexh[] = {      .5,    .5,      .5,      .5,      .4,      .4};
    Double_t nnloeyl[] = { errDown, errDown, errDown, errDown, errDown, errDown};
    Double_t nnloeyh[] = {   errUp,   errUp,   errUp,   errUp,   errUp,   errUp};

    TGraphAsymmErrors *nnloplot = new TGraphAsymmErrors(6, nnlox, nnloy, nnloexl, nnloexh, nnloeyl, nnloeyh);
    nnloplot->SetLineColor(kGreen-7);
    nnloplot->SetLineWidth(4);
    nnloplot->SetFillColor(kGreen-7);
    nnloplot->SetFillStyle(3013);

    TH1F* framehist = new TH1F("framehist", "", 4, 0., 4.);
    framehist->SetMinimum(700);
    framehist->SetMaximum(1000);
    framehist->GetXaxis()->SetTickLength(0);
    framehist->GetXaxis()->SetBinLabel(1, "");
    framehist->GetXaxis()->SetBinLabel(2, "");
    framehist->GetXaxis()->SetBinLabel(3, "");
    framehist->GetYaxis()->SetTitle("#sigma_{t#bar{t}}^{incl} [pb]");
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

    TLegend* leg =  new TLegend(0.37, 0.7, 0.62, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    TString thiswork_label = "This work," + this->lumi_fb + "fb^{-1}";
    leg->AddEntry(mplotwithsys, thiswork_label, "p");
    leg->AddEntry(cmsemu15plot, "CMS, 2.2 fb^{-1} [EPJC 77:172 (2017)]", "p");
    leg->AddEntry(nnloplot, "NNLO+NNLL, m_{t} = 172.5 GeV [PRL 110:252004 (2013)]", "lf");

    TCanvas* c = new TCanvas("plot", "plot", 1200, 800);
    framehist->Draw();
    box1->Draw("SAME");
    box2->Draw("SAME");
    box3->Draw("SAME");
    box4->Draw("SAME");
    nnloplot->Draw("C,2,SAME");

    gStyle->SetEndErrorSize(8);
    cmsemu15plot->Draw("p,SAME");
    cmsemu15plot_stat_syst->Draw("p,SAME,Z");
    mplot->Draw("p,SAME");
    mplotwithsys->Draw("p,SAME,Z");
    leg ->Draw("SAME");

    TPaveText *label = new TPaveText();
    label->SetX1NDC(1.0 - gStyle->GetPadRightMargin() - 0.2);
    label->SetY1NDC(0.98 - gStyle->GetPadTopMargin());
    label->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - 0.03);
    label->SetY2NDC(0.98);
    label->SetTextFont(42);
    label->AddText("13 TeV");
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->SetTextSize(0.04);
    label->SetTextAlign(32);
    label->Draw("SAME");

    TString outdir = utils::assignFolder(this->outpathPlots, Channel, TString("FinalResults"));
    c->Print(outdir.Copy() + "InclusiveXSec.eps");
    c->Print(outdir.Copy() + "InclusiveXSec.C");
    c->Print(outdir.Copy() + "InclusiveXSec.pdf");
    c->Print(outdir.Copy() + "InclusiveXSec.root");
    c->Clear();
    delete c;
}


void PlotterCP::MakeTable(TString Channel, TString Systematic)
{
    setDataSet(Channel, Systematic);
    TH1D *numhists5[hists.size()];
    TH1D *numhists6[hists.size()];
    TH1D *numhists7[hists.size()];
    TH1D *numhists8[hists.size()];
    TH1D *numhists9[hists.size()];
    TH1D *numhists7L[hists.size()];

    TH1D *numhistsForPseudoData5[histsForPseudoData.size()];
    TH1D *numhistsForPseudoData6[histsForPseudoData.size()];
    TH1D *numhistsForPseudoData7[histsForPseudoData.size()];
    TH1D *numhistsForPseudoData8[histsForPseudoData.size()];
    TH1D *numhistsForPseudoData9[histsForPseudoData.size()];
    TH1D *numhistsForPseudoData7L[histsForPseudoData.size()];

    TString histoForYieldTable;
    if (usePoissonSmearedPseudoData) histoForYieldTable = "basic_jet_multiplicity";
    else histoForYieldTable = "events_weighted";

    for (unsigned int i = 0; i < dataset.size(); i++) {
        TH1D *temp_hist5 = fileReader->GetClone<TH1D>(dataset[i], histoForYieldTable + "_step4");
        TH1D *temp_hist6 = fileReader->GetClone<TH1D>(dataset[i], histoForYieldTable + "_step5");
        TH1D *temp_hist7 = fileReader->GetClone<TH1D>(dataset[i], histoForYieldTable + "_step6");
        TH1D *temp_hist8 = fileReader->GetClone<TH1D>(dataset[i], histoForYieldTable + "_step7");
        TH1D *temp_hist9 = fileReader->GetClone<TH1D>(dataset[i], histoForYieldTable + "_step8");
        TH1D *temp_hist7L = fileReader->GetClone<TH1D>(dataset[i], histoForYieldTable + "_step7L");

	double LumiWeight;
	double year_corr_rescale_from_file = 1.;
	if (year == "fullRun2UL" || year == "2016") {
	  year_corr_rescale_from_file = configHelper->get_year_corr_rescale_from_file(Systematic, dataset.at(i));
	  configHelper->get_year_from_file(dataset.at(i));
	  LumiWeight = configHelper->CalcLumiWeight(dataset.at(i), configHelper->lumi_year);
	}
	else
	  LumiWeight = configHelper->CalcLumiWeight(dataset.at(i));

        configHelper->ApplyFlatWeights(temp_hist5, LumiWeight);
        configHelper->ApplyFlatWeights(temp_hist6, LumiWeight);
        configHelper->ApplyFlatWeights(temp_hist7, LumiWeight);
        configHelper->ApplyFlatWeights(temp_hist8, LumiWeight);
        configHelper->ApplyFlatWeights(temp_hist9, LumiWeight);
        configHelper->ApplyFlatWeights(temp_hist7L, LumiWeight);

	//printf("dataset.at(%d) = %s \n", i, dataset.at(i).Data());
	//printf("year_corr_rescale_from_file = %f \n", year_corr_rescale_from_file);

	if (year_corr_rescale_from_file != 1.) {
	  TString datasetNominal = configHelper->getCorrespondingNominalDataset(dataset.at(i), Systematic, Channel);

	  //printf("datasetNominal = %s \n", datasetNominal.Data());

	  TH1D *temp_hist5_nominal = fileReader->GetClone<TH1D>(datasetNominal, histoForYieldTable + "_step4");
	  TH1D *temp_hist6_nominal = fileReader->GetClone<TH1D>(datasetNominal, histoForYieldTable + "_step5");
	  TH1D *temp_hist7_nominal = fileReader->GetClone<TH1D>(datasetNominal, histoForYieldTable + "_step6");
	  TH1D *temp_hist8_nominal = fileReader->GetClone<TH1D>(datasetNominal, histoForYieldTable + "_step7");
	  TH1D *temp_hist9_nominal = fileReader->GetClone<TH1D>(datasetNominal, histoForYieldTable + "_step8");
	  TH1D *temp_hist7L_nominal = fileReader->GetClone<TH1D>(datasetNominal, histoForYieldTable + "_step7L");

	  double LumiWeightNominal = configHelper->CalcLumiWeight(datasetNominal, configHelper->lumi_year);

	  configHelper->ApplyFlatWeights(temp_hist5_nominal, LumiWeightNominal);
	  configHelper->ApplyFlatWeights(temp_hist6_nominal, LumiWeightNominal);
	  configHelper->ApplyFlatWeights(temp_hist7_nominal, LumiWeightNominal);
	  configHelper->ApplyFlatWeights(temp_hist8_nominal, LumiWeightNominal);
	  configHelper->ApplyFlatWeights(temp_hist9_nominal, LumiWeightNominal);
	  configHelper->ApplyFlatWeights(temp_hist7L_nominal, LumiWeightNominal);

	  configHelper->RescaleHistYearToYearCorr(temp_hist5, temp_hist5_nominal, year_corr_rescale_from_file);
	  configHelper->RescaleHistYearToYearCorr(temp_hist6, temp_hist6_nominal, year_corr_rescale_from_file);
	  configHelper->RescaleHistYearToYearCorr(temp_hist7, temp_hist7_nominal, year_corr_rescale_from_file);
	  configHelper->RescaleHistYearToYearCorr(temp_hist8, temp_hist8_nominal, year_corr_rescale_from_file);
	  configHelper->RescaleHistYearToYearCorr(temp_hist9, temp_hist9_nominal, year_corr_rescale_from_file);
	  configHelper->RescaleHistYearToYearCorr(temp_hist7L, temp_hist7L_nominal, year_corr_rescale_from_file);
	}

        removeNegativeExpectations(temp_hist5);
        removeNegativeExpectations(temp_hist6);
        removeNegativeExpectations(temp_hist7);
        removeNegativeExpectations(temp_hist8);
        removeNegativeExpectations(temp_hist9);
        removeNegativeExpectations(temp_hist7L);

        numhists5[i] = temp_hist5;
        numhists6[i] = temp_hist6;
        numhists7[i] = temp_hist7;
        numhists8[i] = temp_hist8;
        numhists9[i] = temp_hist9;
        numhists7L[i] = temp_hist7L;
    }

    if (usePoissonSmearedPseudoData) {
        for (unsigned int i = 0; i < pseudoDataset.size(); i++) {
            TH1D *temp_hist5_pd = fileReader->GetClone<TH1D>(pseudoDataset[i], histoForYieldTable + "_step4");
            TH1D *temp_hist6_pd = fileReader->GetClone<TH1D>(pseudoDataset[i], histoForYieldTable + "_step5");
            TH1D *temp_hist7_pd = fileReader->GetClone<TH1D>(pseudoDataset[i], histoForYieldTable + "_step6");
            TH1D *temp_hist8_pd = fileReader->GetClone<TH1D>(pseudoDataset[i], histoForYieldTable + "_step7");
            TH1D *temp_hist9_pd = fileReader->GetClone<TH1D>(pseudoDataset[i], histoForYieldTable + "_step8");
            TH1D *temp_hist7L_pd = fileReader->GetClone<TH1D>(pseudoDataset[i], histoForYieldTable + "_step7L");

            double LumiWeight_pd;
	    double year_corr_rescale_from_file_pd = 1.;
	    if (year == "fullRun2UL" || year == "2016") {
	      year_corr_rescale_from_file_pd = configHelper->get_year_corr_rescale_from_file(Systematic, pseudoDataset.at(i));
	      configHelper->get_year_from_file(pseudoDataset.at(i));
	      LumiWeight_pd = configHelper->CalcLumiWeight(pseudoDataset.at(i), configHelper->lumi_year);
	    }
	    else
	      LumiWeight_pd = configHelper->CalcLumiWeight(pseudoDataset.at(i));

            configHelper->ApplyFlatWeights(temp_hist5_pd, LumiWeight_pd);
            configHelper->ApplyFlatWeights(temp_hist6_pd, LumiWeight_pd);
            configHelper->ApplyFlatWeights(temp_hist7_pd, LumiWeight_pd);
            configHelper->ApplyFlatWeights(temp_hist8_pd, LumiWeight_pd);
            configHelper->ApplyFlatWeights(temp_hist9_pd, LumiWeight_pd);
            configHelper->ApplyFlatWeights(temp_hist7L_pd, LumiWeight_pd);

	    if (year_corr_rescale_from_file_pd != 1.) {
	      TString pseudoDatasetNominal = configHelper->getCorrespondingNominalDataset(pseudoDataset[i], Systematic, Channel);

	      TH1D *temp_hist5_nominal_pd = fileReader->GetClone<TH1D>(pseudoDatasetNominal, histoForYieldTable + "_step4");
	      TH1D *temp_hist6_nominal_pd = fileReader->GetClone<TH1D>(pseudoDatasetNominal, histoForYieldTable + "_step5");
	      TH1D *temp_hist7_nominal_pd = fileReader->GetClone<TH1D>(pseudoDatasetNominal, histoForYieldTable + "_step6");
	      TH1D *temp_hist8_nominal_pd = fileReader->GetClone<TH1D>(pseudoDatasetNominal, histoForYieldTable + "_step7");
	      TH1D *temp_hist9_nominal_pd = fileReader->GetClone<TH1D>(pseudoDatasetNominal, histoForYieldTable + "_step8");
	      TH1D *temp_hist7L_nominal_pd = fileReader->GetClone<TH1D>(pseudoDatasetNominal, histoForYieldTable + "_step7L");

	      double LumiWeightNominal_pd = configHelper->CalcLumiWeight(pseudoDatasetNominal, configHelper->lumi_year);

	      configHelper->ApplyFlatWeights(temp_hist5_nominal_pd, LumiWeightNominal_pd);
	      configHelper->ApplyFlatWeights(temp_hist6_nominal_pd, LumiWeightNominal_pd);
	      configHelper->ApplyFlatWeights(temp_hist7_nominal_pd, LumiWeightNominal_pd);
	      configHelper->ApplyFlatWeights(temp_hist8_nominal_pd, LumiWeightNominal_pd);
	      configHelper->ApplyFlatWeights(temp_hist9_nominal_pd, LumiWeightNominal_pd);
	      configHelper->ApplyFlatWeights(temp_hist7L_nominal_pd, LumiWeightNominal_pd);

	      configHelper->RescaleHistYearToYearCorr(temp_hist5_pd, temp_hist5_nominal_pd, year_corr_rescale_from_file_pd );
	      configHelper->RescaleHistYearToYearCorr(temp_hist6_pd, temp_hist6_nominal_pd, year_corr_rescale_from_file_pd );
	      configHelper->RescaleHistYearToYearCorr(temp_hist7_pd, temp_hist7_nominal_pd, year_corr_rescale_from_file_pd );
	      configHelper->RescaleHistYearToYearCorr(temp_hist8_pd, temp_hist8_nominal_pd, year_corr_rescale_from_file_pd );
	      configHelper->RescaleHistYearToYearCorr(temp_hist9_pd, temp_hist9_nominal_pd, year_corr_rescale_from_file_pd );
	      configHelper->RescaleHistYearToYearCorr(temp_hist7L_pd, temp_hist7L_nominal_pd, year_corr_rescale_from_file_pd );
	    }

            removeNegativeExpectations(temp_hist5_pd);
            removeNegativeExpectations(temp_hist6_pd);
            removeNegativeExpectations(temp_hist7_pd);
            removeNegativeExpectations(temp_hist8_pd);
            removeNegativeExpectations(temp_hist9_pd);
            removeNegativeExpectations(temp_hist7L_pd);

            numhistsForPseudoData5[i] = temp_hist5_pd;
            numhistsForPseudoData6[i] = temp_hist6_pd;
            numhistsForPseudoData7[i] = temp_hist7_pd;
            numhistsForPseudoData8[i] = temp_hist8_pd;
            numhistsForPseudoData9[i] = temp_hist9_pd;
            numhistsForPseudoData7L[i] = temp_hist7L_pd;
        }

        if (Channel != "combined") {
            for (size_t i = 1; i < pseudoDataset.size(); ++i) {
	      if (!pseudoDataset.at(i).Contains("ttbarsignalplustau") && !pseudoDataset.at(i).Contains("ttbarsignalviatau") && !pseudoDataset.at(i).Contains("ttbarbgviatau") && !pseudoDataset.at(i).Contains("run")) {
                    numhists5[0]->Add(numhistsForPseudoData5[i]);
                    numhists6[0]->Add(numhistsForPseudoData6[i]);
                    numhists7[0]->Add(numhistsForPseudoData7[i]);
                    numhists8[0]->Add(numhistsForPseudoData8[i]);
                    numhists9[0]->Add(numhistsForPseudoData9[i]);
                    numhists7L[0]->Add(numhistsForPseudoData7L[i]);
                }
            }
            PlotterCP::performPoissonSmearing(numhists5[0]);
            PlotterCP::performPoissonSmearing(numhists6[0]);
            PlotterCP::performPoissonSmearing(numhists7[0]);
            PlotterCP::performPoissonSmearing(numhists8[0]);
            PlotterCP::performPoissonSmearing(numhists9[0]);
            PlotterCP::performPoissonSmearing(numhists7L[0]);
        } else if (Channel == "combined") {
            //Indeces must match to the order of data files in filelist {syst}_combined.txt
            TString chToInclude[3] = {"ee_", "emu_", "mumu_"};
            for (int ch_it = 0; ch_it < 3; ++ch_it) {
                for (size_t i = ch_it + 1; i < pseudoDataset.size(); ++i) {//skipping pseudo-data files (moving iterator by +1 when initialized) to which MC is added
		  if (!pseudoDataset.at(i).Contains("ttbarsignalplustau") && !pseudoDataset.at(i).Contains("ttbarsignalviatau") && !pseudoDataset.at(i).Contains("ttbarbgviatau") &&
                        !pseudoDataset.at(i).Contains("run") && pseudoDataset.at(i).Contains(chToInclude[ch_it]))
                    {
                        numhists5[ch_it]->Add(numhistsForPseudoData5[i]);
                        numhists6[ch_it]->Add(numhistsForPseudoData6[i]);
                        numhists7[ch_it]->Add(numhistsForPseudoData7[i]);
                        numhists8[ch_it]->Add(numhistsForPseudoData8[i]);
                        numhists9[ch_it]->Add(numhistsForPseudoData9[i]);
                        numhists7L[ch_it]->Add(numhistsForPseudoData7L[i]);
                    }
                }
                PlotterCP::performPoissonSmearing(numhists5[ch_it]);
                PlotterCP::performPoissonSmearing(numhists6[ch_it]);
                PlotterCP::performPoissonSmearing(numhists7[ch_it]);
                PlotterCP::performPoissonSmearing(numhists8[ch_it]);
                PlotterCP::performPoissonSmearing(numhists9[ch_it]);
                PlotterCP::performPoissonSmearing(numhists7L[ch_it]);
            }
        }
    } else {
        for (unsigned int i = 0; i < hists.size(); i++) {
            if (legends.at(i) == DYEntry) {
                //numhists5[i]->Scale(DYScale[channelType]);//Has to be consistent with the step where the application of DY SF starts
                numhists6[i]->Scale(DYScale[channelType]);
                numhists7[i]->Scale(DYScale.at(channelType));
                numhists8[i]->Scale(DYScale.at(channelType));
                numhists9[i]->Scale(DYScale.at(channelType));
                numhists7L[i]->Scale(DYScale.at(channelType));
            }
        }
    }

    ////////////////////////////Make output for tables
    double tmp_num5 = 0;
    double tmp_num6 = 0;
    double tmp_num7 = 0;
    double tmp_num8 = 0;
    double tmp_num9 = 0;
    double tmp_num7L = 0;

    TString outdir = utils::assignFolder(this->outpathPlots, Channel, Systematic);
    ofstream EventFile5; EventFile5.open(outdir.Copy() + "Events5.txt");
    ofstream EventFile6; EventFile6.open(outdir.Copy() + "Events6.txt");
    ofstream EventFile7; EventFile7.open(outdir.Copy() + "Events7.txt");
    ofstream EventFile8; EventFile8.open(outdir.Copy() + "Events8.txt");
    ofstream EventFile9; EventFile9.open(outdir.Copy() + "Events9.txt");
    ofstream EventFile7L; EventFile7L.open(outdir.Copy() + "Events7L.txt");

    double bg_num5 = 0;
    double bg_num6 = 0;
    double bg_num7 = 0;
    double bg_num8 = 0;
    double bg_num9 = 0;
    double bg_num7L = 0;

    for (unsigned int i = 0; i < hists.size(); i++) {
        if (usePoissonSmearedPseudoData) {
            tmp_num5 += numhists5[i]->Integral(0, numhists5[i]->GetNbinsX()+1);//considering use of over/underflow bins
            tmp_num6 += numhists6[i]->Integral(0, numhists6[i]->GetNbinsX()+1);
            tmp_num7 += numhists7[i]->Integral(0, numhists7[i]->GetNbinsX()+1);
            tmp_num8 += numhists8[i]->Integral(0, numhists8[i]->GetNbinsX()+1);
            tmp_num9 += numhists9[i]->Integral(0, numhists9[i]->GetNbinsX()+1);
            tmp_num7L += numhists7L[i]->Integral(0, numhists7L[i]->GetNbinsX()+1);
        } else {
            tmp_num5 += numhists5[i]->Integral();
            tmp_num6 += numhists6[i]->Integral();
            tmp_num7 += numhists7[i]->Integral();
            tmp_num8 += numhists8[i]->Integral();
            tmp_num9 += numhists9[i]->Integral();
            tmp_num7L += numhists7L[i]->Integral();
        }

        if (i == (hists.size()-1)) {
            EventFile5 << legends.at(i) << ": " << tmp_num5 << std::endl;
            EventFile6 << legends.at(i) << ": " << tmp_num6 << std::endl;
            EventFile7 << legends.at(i) << ": " << tmp_num7 << std::endl;
            EventFile8 << legends.at(i) << ": " << tmp_num8 << std::endl;
            EventFile9 << legends.at(i) << ": " << tmp_num9 << std::endl;
            EventFile7L << legends.at(i) << ": " << tmp_num7L << std::endl;
            bg_num5 += tmp_num5;
            bg_num6 += tmp_num6;
            bg_num7 += tmp_num7;
            bg_num8 += tmp_num8;
            bg_num9 += tmp_num9;
            bg_num7L += tmp_num7L;
            tmp_num5 = 0;
            tmp_num6 = 0;
            tmp_num7 = 0;
            tmp_num8 = 0;
            tmp_num9 = 0;
            tmp_num7L = 0;
        } else if (legends.at(i) != legends.at(i+1)) {
            EventFile5 << legends.at(i) << ": " << tmp_num5 << std::endl;
            EventFile6 << legends.at(i) << ": " << tmp_num6 << std::endl;
            EventFile7 << legends.at(i) << ": " << tmp_num7 << std::endl;
            EventFile8 << legends.at(i) << ": " << tmp_num8 << std::endl;
            EventFile9 << legends.at(i) << ": " << tmp_num9 << std::endl;
            EventFile7L << legends.at(i) << ": " << tmp_num7L << std::endl;
            if (legends.at(i) != "Data") {
                bg_num5 += tmp_num5;
                bg_num6 += tmp_num6;
                bg_num7 += tmp_num7;
                bg_num8 += tmp_num8;
                bg_num9 += tmp_num9;
                bg_num7L += tmp_num7L;
            }
            tmp_num5 = 0;
            tmp_num6 = 0;
            tmp_num7 = 0;
            tmp_num8 = 0;
            tmp_num9 = 0;
            tmp_num7L = 0;
        }
    }
    EventFile5 << "Total background: " << bg_num5 << std::endl;
    EventFile5.close();
    EventFile6 << "Total background: " << bg_num6 << std::endl;
    EventFile6.close();
    EventFile7 << "Total background: " << bg_num7 << std::endl;
    EventFile7.close();
    EventFile7L << "Total background: " << bg_num7L << std::endl;
    EventFile7L.close();
    EventFile8 << "Total background: " << bg_num8 << std::endl;
    EventFile8.close();
    EventFile9 << "Total background: " << bg_num9 << std::endl;
    EventFile9.close();
    std::cout << "\nEvent yields saved in " << outdir << std::endl;
}


bool PlotterCP::addQCDToControlPlot()const
{
    if ((name.Contains("_step") && !name.Contains("7") && !name.Contains("8")) ||
        (name.Contains("events_") && !name.Contains("7") && !name.Contains("8")) ||
        (!name.Contains("Hyp") && (name.Contains("jetHT") || name.Contains("jetpT"))) ||
        name == "MET" || name.Contains("_noBTag") || name.Contains("_diLep") || name.Contains("Electron") || name.Contains("Muon"))
       {
            return 0;//for now, permanently switched off, use: 1 in case you want QCD plotted
       }
    return 0;
}


void PlotterCP::GetEnvelopeVariationsForBandsInCPs(TH1 * varHist, const TString Systematic, const TString channel)
{
    if (!varHist) return;

    TString systematic_name = utils::GetNameFromUpDownVariation(Systematic);

    if (!this->configHelper->IsEnvelope(systematic_name)) {
        std::cout << "Systematics name: " << Systematic << std::endl;
        std::cout << "The calculation of envelope is not supported. Exiting!!" << std::endl;
        exit(433);
    }

    std::vector<TString> syst = this->configHelper->GetEnvelopeVars(systematic_name);

    /*
    if (Systematic.Contains("TOT_SCALE_")) syst = {"MESCALE_UP", "MESCALE_DOWN", "MEFACSCALE_UP", "MEFACSCALE_DOWN", "MERENSCALE_UP", "MERENSCALE_DOWN"};
    else if (Systematic.Contains("TOT_BFRAG_")) syst = {"BFRAG_UP", "BFRAG_DOWN", "BFRAG_PETERSON"};
    // if (Systematic.Contains("TOT_COLORREC_")) syst = {"ERDON", "ERDONRETUNE", "GLUONMOVETUNE"};
    else if (Systematic.Contains("TOT_COLORREC_")) syst = {"ERDON"};
    //FIXME: merging envelope for b-tagging and trigger variations needed?*/

    unsigned int nbins = varHist->GetNbinsX();

    TString file_nom = TString(this->outpathPlots + "Nominal/" + channel + "/" + name + "_source.root");
    ifstream inputFileStream(file_nom);
    if (inputFileStream.fail()) return;
    inputFileStream.close();

    TH1D *tmp_nom = nullptr;
    tmp_nom = fileReader->GetClone<TH1D>(file_nom, varHist->GetName(), 1);
    tmp_nom->SetDirectory(0);

    for (size_t iter = 0; iter < syst.size(); iter++) {
        TString file = TString(this->outpathPlots + syst.at(iter) + "/" + channel + "/" + name + "_source.root");
        std::cout << "Adding to " << Systematic << ": " << file << std::endl;
        ifstream inputFileStream(file);

        if (inputFileStream.fail()) {
            std::cout << "** WARNING ** File " << file << " not oppened!" << std::endl;
            continue;
        }
        inputFileStream.close();

        TH1D *tmp = nullptr;

        TString thisName = varHist->GetName();
        tmp = fileReader->GetClone<TH1D>(file, thisName, 1);

        tmp->SetDirectory(0);

        for (size_t ibin = 1; ibin <= nbins; ++ibin) {
            double tmp_bin = tmp->GetBinContent(ibin);
            // if (syst.at(iter).Contains("PSFSRSCALE_")) tmp_bin = tmp_nom->GetBinContent(ibin) + (tmp_bin - tmp_nom->GetBinContent(ibin)) / fsrScaleUncReductionFactor;
            if (Systematic.Contains("_DOWN")) {
                if (tmp_bin < varHist->GetBinContent(ibin)) varHist->SetBinContent(ibin, tmp_bin);
            } else {
                if (tmp_bin > varHist->GetBinContent(ibin)) varHist->SetBinContent(ibin, tmp_bin);
            }
        }

        delete tmp;
    }

    delete tmp_nom;

    return;
}


void PlotterCP::getSignalUncertaintyBand(TH1* uncBand, std::map<TString, std::vector<double> > &varUnc, TString channel_)
{
  if(!uncBand)  return;

  //Using systematics defined in PlotterConfigurationHelper
  //std::vector<TString> syst = this->configHelper->_vectorOfSystematics_all;
  std::vector<TString> syst;
  syst = this->uncertaintyBandSystematics;

  std::cout << "Using these systematics for Signal Uncertainty Band: ";
  for (auto syst_i: syst) std::cout << syst_i << "	";
  std::cout << std::endl;

  //Uncomment this for customized systematics
  /*std::vector<TString> syst {
      // "MASS_", "BSEMILEP_", "UETUNE_", "MATCH_", "PSISRSCALE_", "PSFSRSCALE_", "PDF_ALPHAS_",
         "MASS_", "BSEMILEP_", "UETUNE_", "MATCH_", "PSISRSCALE_2_", "PSFSRSCALE_2_", "PDF_ALPHAS_",
         "JESAbsoluteStat_", "JESAbsoluteScale_", "JESAbsoluteMPFBias_", "JESFragmentation_",
         "JESSinglePionECAL_", "JESSinglePionHCAL_", "JESFlavorQCD_", "JESTimePtEta_",
         "JESRelativeBal_", "JESRelativeJEREC1_", "JESRelativePtBB_", "JESRelativePtEC1_",
         "JESRelativeFSR_", "JESRelativeStatFSR_", "JESRelativeStatEC_", "JESPileUpDataMC_",
         "JESPileUpPtRef_", "JESPileUpPtEC1_", "JESPileUpPtBB_",
         "JER_", "PU_", "LEPT_", "L1PREFIRING_", "TRIG_", "BTAG_", "BTAG_LJET_", "UNCLUSTERED_",
      // "MESCALE_", "MEFACSCALE_", "MERENSCALE_","BFRAG_",
         "TOT_SCALE_", "TOT_BFRAG_", "TOT_COLORREC_"
  };*/

  std::vector<TString> syst_files_not_found;

  // List of systematics that could be included if needed:
  // "HAD_", "MOD_", "BTAG_PT_", "BTAG_ETA_", "BTAG_LJET_PT_", "BTAG_LJET_ETA_",

    int nbins = uncBand->GetNbinsX();
    std::vector<double> vec_varup(nbins, 0), vec_vardown(nbins, 0);
    for (size_t iter = 0; iter < syst.size(); iter++) {
        if (syst.at(iter).Contains("LUMI")) continue; 
        TString file_up = "", file_do = "";
        // Handling samples without UP/DOWN
        //TODO: Create a general function for this
        if (syst.at(iter) == "HAD") {
            file_up = TString(this->outpathPlots + "/AMCATNLOFXFX/" + channel_ + "/" + name + "_source.root");
            file_do = TString(this->outpathPlots + "/AMCATNLOFXFX/" + channel_ + "/" + name + "_source.root");
        } else if (syst.at(iter) == "TOP_PT") {
            file_up = TString(this->outpathPlots + "/TOP_PT/" + channel_ + "/" + name + "_source.root");
            file_do = TString(this->outpathPlots + "/TOP_PT/" + channel_ + "/" + name + "_source.root");
        } else if (syst.at(iter) == "MOD") {
            file_up = TString(this->outpathPlots + "/POWHEG/" + channel_ + "/" + name + "_source.root");
            file_do = TString(this->outpathPlots + "/POWHEG/" + channel_ + "/" + name + "_source.root");
        } else if (configHelper->IsNoUpDown(syst.at(iter))) {
            file_up = TString(this->outpathPlots + "/" + syst.at(iter)+ "/" + channel_ + "/" + name + "_source.root");
            file_do = TString(this->outpathPlots + "/" + syst.at(iter)+ "/" + channel_ + "/" + name + "_source.root");
        } else if (syst.at(iter) == "POWHEG") {
            continue;
        } else {
            file_up = TString(this->outpathPlots + syst.at(iter) + "_UP/" + channel_ + "/" + name + "_source.root");
            file_do = TString(this->outpathPlots + syst.at(iter) + "_DOWN/" + channel_ + "/" + name + "_source.root");
        }

        ifstream inputFileStream(file_up);
        if (inputFileStream.fail()) {
            syst_files_not_found.push_back(file_up);
            continue;
        }
        inputFileStream.close();
        ifstream inputFileStream2(file_do);
        if (inputFileStream2.fail()) {
            syst_files_not_found.push_back(file_do);
            continue;
        }
        inputFileStream2.close();

        // This lines crashes the code, some probles arises form the HistoListReader class
        TH1D *tmpUp = nullptr;
        TH1D *tmpDo = nullptr;

        if (drawTotalMCErrorForCP) {
            tmpUp = fileReader->GetClone<TH1D>(file_up, name + "_allmc", 1);
            tmpDo = fileReader->GetClone<TH1D>(file_do, name + "_allmc", 1);
        } else {
            tmpUp = fileReader->GetClone<TH1D>(file_up, name + "_allttbar", 1);
            tmpDo = fileReader->GetClone<TH1D>(file_do, name + "_allttbar", 1);
        }

        if (!tmpUp && tmpDo) {delete tmpDo; continue;}
        if (tmpUp && !tmpDo) {delete tmpUp; continue;}
        if (!tmpUp && !tmpDo) continue;

        if (nbins != tmpUp->GetNbinsX() || nbins != tmpDo->GetNbinsX()) continue;
        tmpDo->SetDirectory(0);
        tmpUp->SetDirectory(0);

       //        if(syst.at(iter).Contains("PSSCALE_WEIGHT")){
       //        double norm_Events = uncBand->Integral(-1e6, 1e6);
       //        double upIntegral = tmpUp->Integral(-1e6, 1e6);
       //        double doIntegral = tmpDo->Integral(-1e6, 1e6);
       //        tmpUp->Scale(norm_Events/upIntegral);
       //        tmpDo->Scale(norm_Events/doIntegral);

        double factorVar = 1.0;
        if(syst.at(iter).Contains("MASS")) factorVar = massUncReductionFactor;
        // if(syst.at(iter).Contains("PSFSRSCALE_2")) factorVar = fsrScaleUncReductionFactor;

        double temp_uncert_u = 0.0;
        double temp_uncert_d = 0.0;
        std::vector<double> all_temp_uncert_u, all_temp_uncert_d;
        for (Int_t nbin = 0; nbin < nbins; nbin++) {
            double binContent = uncBand->GetBinContent(nbin+1);

	    double rel_diffup = (tmpUp->GetBinContent(nbin+1) - binContent) / binContent;
	    double rel_diffdo = (tmpDo->GetBinContent(nbin+1) - binContent) / binContent;
            if (binContent < 1e-6) {
                rel_diffdo = 0;
                rel_diffup = 0;
            }

            double mod_diff_u = std::abs(rel_diffup) / factorVar;
            double mod_diff_d = std::abs(rel_diffdo) / factorVar;
            vec_varup.at(nbin) += (mod_diff_u) * (mod_diff_u);
            vec_vardown.at(nbin) += (mod_diff_d) * (mod_diff_d);

            temp_uncert_u += (mod_diff_u) * (mod_diff_u);
            temp_uncert_d += (mod_diff_d) * (mod_diff_d);

            all_temp_uncert_u.push_back(rel_diffup);
            all_temp_uncert_d.push_back(rel_diffdo);

	}

        std::cout << syst.at(iter) << "(%): UP=" << (std::sqrt(temp_uncert_u)/std::sqrt(nbins))*100 << "; DOWN=" << (std::sqrt(temp_uncert_d)/std::sqrt(nbins))*100 << std::endl;
        if (configHelper->IsNoUpDown(syst.at(iter))) varUnc.insert(std::pair<TString, std::vector<double> >(syst.at(iter), all_temp_uncert_u));
        else {
            varUnc.insert(std::pair<TString, std::vector<double>>(syst.at(iter)+"_UP", all_temp_uncert_u));
            varUnc.insert(std::pair<TString, std::vector<double>>(syst.at(iter)+"_DOWN", all_temp_uncert_d));
        }

        factorVar = 1.0;
        delete tmpDo; delete tmpUp;
    }
    double lumiError2 = lumiError*lumiError;
    double bRatioError2 = bRatioError*bRatioError;
    for (size_t iter = 0; iter < vec_varup.size(); iter++) {
        double centralValue = uncBand->GetBinContent(iter+1);
        // Adding lumi uncertainty here
        vec_varup.at(iter) += lumiError2;
        vec_varup.at(iter) += bRatioError2;
        vec_vardown.at(iter) += lumiError2;
        vec_vardown.at(iter) += bRatioError2;
        
        uncBand->SetBinError(iter+1, centralValue * 0.5 * (std::sqrt(vec_vardown.at(iter)) + std::sqrt(vec_varup.at(iter))));
    }

    //Print warning if files weren't opened
    if (syst_files_not_found.size()>0){
        std::cout << "      WARNING: " << syst_files_not_found.size() << " file(s) not found for uncertainty band plotting." << std::endl
                  << "          Files not found: ";
        for (auto file: syst_files_not_found) std::cout << file << " ";
        std::cout << endl << endl;
    }

}


void PlotterCP::plotUncertaintyBandContributions(std::map<TString, std::vector<double> > &varUnc, const TString channel_, TH1* refHisto, const std::vector<TString>& syst_list, const TString& fileNameSufix) {
    TString out_base_path = this->outpathPlots + "/Nominal/" + channel_ + "/" + name + "_uncBreakdown" + fileNameSufix;

    double uncRange = .01;
    // uncRange = 7;
    // double offbin = 0.5;
    // double eps = 1e-4;
    double minY = -1 * uncRange;
    double maxY = uncRange;
    int nbins = varUnc.begin()->second.size();
    int widthX = 800;
    int widthY = 600;
    double lineWidth = 1.0;

    TCanvas* cc = new TCanvas(name + "_uncSummaryCanvas" + fileNameSufix, name + "_uncSummaryCanvas" + fileNameSufix, widthX, widthY);
    cc->cd();

    TPad* p1 = new TPad(name + "_uncSummaryPad", name + "_uncSummaryPad", 0.0, 0.0, 0.75, 1.0);
    p1->SetTopMargin(0.10);
    p1->SetBottomMargin(0.15);
    p1->SetLeftMargin(0.18);
    p1->SetRightMargin(0.05);
    TPad* p2 = new TPad(name + "_LegendPad", name + "_LegendPad", 0.75, 0.0, 1.0, 1.0);
    p2->SetTopMargin(0.025);
    p2->SetBottomMargin(0.025);
    p2->SetLeftMargin(0.0);
    p2->SetRightMargin(0.025);

    //    double legX0 = 0.825;
    //    double legX1 = 0.99;
    double legX0 = 0.0;
    double legX1 = 0.99;


    TLegend* leg = new TLegend(legX0, 0.12, legX1, 0.94);
    leg->SetNColumns(1);
    leg->SetFillStyle(0);
    leg->SetEntrySeparation(0.0);
    leg->SetColumnSeparation(0.0);
    leg->SetTextAlign(12);
    leg->SetBorderSize(0);
    leg->SetLineWidth(2);
    leg->SetTextSize(0.07);

    p1->cd();

    TH1 *hrange = nullptr;
    hrange = dynamic_cast<TH1*> (refHisto->Clone(name + "_uncSummaryHRange" + fileNameSufix));
    for(int b = 0; b < refHisto->GetNbinsX(); b++) {
        hrange->SetBinContent(b+1,0);
        hrange->SetBinError(b+1, 0);
    }
    hrange->GetXaxis()->SetTitle(XAxis);
    hrange->GetXaxis()->SetTickLength(0.04);
    hrange->GetYaxis()->SetTitle("(Nominal-Syst)/Nominal");
    hrange->GetYaxis()->SetTitleSize(0.04);
    hrange->GetYaxis()->SetTitleOffset(1.9);
    hrange->GetXaxis()->SetTitleOffset(1.1);
    hrange->GetYaxis()->SetLabelSize(0.04);
    hrange->SetLineColor(1);
    hrange->SetFillColor(0);
    hrange->GetYaxis()->SetNdivisions(510);
    if(name.Contains("_LLBar") || 
       name.Contains("_Lepton") || 
       name.Contains("_AntiLepton") || 
       name.Contains("_TTBarMass") ||
       name.Contains("_ScatteringAngle_TTBarFrame") ) {
      hrange->GetXaxis()->SetNdivisions(-20604);
    }
    hrange->GetYaxis()->SetRangeUser(minY, maxY);
    hrange->SetMinimum(minY);
    hrange->SetMaximum(maxY);
    if (rangemin != 0 || rangemax != 0) hrange->SetAxisRange(rangemin, rangemax, "X");

    hrange->Draw();
    if(name.Contains("_vs_TTBarMass") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets")) {
      p1->SetGridx();
      draw2DLabels(0.825);
    }    

    // fill plots here
    int id_style = 0;
    int i_color = 1;
    int lineStyles[6] = {kSolid, kDashed, kDotted, kDashDotted, 2, 4};
    int step_for_style = 8;

    const std::vector<TString> *syst = nullptr;
    //    if (syst_list.size() == 0) syst = &this->uncertaintyBandSystematicsShapeOnly;
    if (syst_list.size() == 0) syst = &this->uncertaintyBandSystematics;
    else syst = &syst_list;

    for (auto const& i_syst : *syst){
        TString leg_mod_syst = i_syst;
        if (i_syst == "ELE_SCALESMEARING_REDUCED_noPhi") leg_mod_syst = "ELE_SCALE_(noPhi)";
        else if (i_syst == "ELE_SCALESMEARING_REDUCED_Phi") leg_mod_syst = "ELE_SCALE_(0.5Phi)";
        else if (i_syst == "ELE_SCALESMEARING_REDUCED_noPhi_noRho") leg_mod_syst = "ELE_SCALE_(noPhi,noRho)";
        else if (i_syst == "ELE_SCALESMEARING") leg_mod_syst = "ELE_SCALE";
        else if (i_syst == "PSSCALE_WEIGHT_2") leg_mod_syst = "ISR Reduced";
        else if (i_syst == "PSSCALE_WEIGHT_3") leg_mod_syst = "FSR Reduced";
        else if (i_syst == "PSSCALE_WEIGHT_4") leg_mod_syst = "ISR Default";
        else if (i_syst == "PSSCALE_WEIGHT_5") leg_mod_syst = "FSR Default";
        else if (i_syst == "PSSCALE_WEIGHT_6") leg_mod_syst = "ISR Conservative";
        else if (i_syst == "PSSCALE_WEIGHT_7") leg_mod_syst = "FSR Conservative";
        else if (i_syst == "PSSCALE_WEIGHT_8") leg_mod_syst = "fsr_G2GG_muR";
        else if (i_syst == "PSSCALE_WEIGHT_9") leg_mod_syst = "fsr_G2QQ_muR";
        else if (i_syst == "PSSCALE_WEIGHT_10") leg_mod_syst = "fsr_Q2QG_muR";
        else if (i_syst == "PSSCALE_WEIGHT_11") leg_mod_syst = "fsr_X2XG_muR";
        else if (i_syst == "PSSCALE_WEIGHT_12") leg_mod_syst = "fsr_G2GG_cNS";
        else if (i_syst == "PSSCALE_WEIGHT_13") leg_mod_syst = "fsr_G2QQ_cNS";
        else if (i_syst == "PSSCALE_WEIGHT_14") leg_mod_syst = "fsr_Q2QG_cNS";
        else if (i_syst == "PSSCALE_WEIGHT_15") leg_mod_syst = "fsr_X2XG_cNS";
        else if (i_syst == "PSSCALE_WEIGHT_16") leg_mod_syst = "isr_G2GG_muR";
        else if (i_syst == "PSSCALE_WEIGHT_17") leg_mod_syst = "isr_G2QQ_muR";
        else if (i_syst == "PSSCALE_WEIGHT_18") leg_mod_syst = "isr_Q2QG_muR";
        else if (i_syst == "PSSCALE_WEIGHT_19") leg_mod_syst = "isr_X2XG_muR";
        else if (i_syst == "PSSCALE_WEIGHT_20") leg_mod_syst = "isr_G2GG_cNS";
        else if (i_syst == "PSSCALE_WEIGHT_21") leg_mod_syst = "isr_G2QQ_cNS";
        else if (i_syst == "PSSCALE_WEIGHT_22") leg_mod_syst = "isr_Q2QG_cNS";
        else if (i_syst == "PSSCALE_WEIGHT_23") leg_mod_syst = "isr_X2XG_cNS";

        else if (i_syst == "BFRAG_PETERSON") leg_mod_syst = "BFRAG_P";
        else if (i_syst == "ELE_ID_new") leg_mod_syst = "ELE_ID (with DY extr.)";
        else if (i_syst == "MUON_ISO_new") leg_mod_syst = "MUON_ISO (with DY extr.)";

        int nvar = 1;
        for (auto & varName : configHelper->GetSysVars(i_syst, true)){
            std::vector<double> bin_vals = varUnc[varName];
            TH1 *hDashedSyst = nullptr;
            hDashedSyst = dynamic_cast<TH1*> (refHisto->Clone(name + "_unc_"+varName + fileNameSufix));

	    //	    std::cout << "varName " << varName << std::endl;
            for(int b = 0; b < nbins; b++) {
                hDashedSyst->SetBinContent(b+1, bin_vals.at(b));
            }
            hDashedSyst->SetLineColor(utils::ttmd::GetDistinctColor(i_color));
            hDashedSyst->SetLineStyle(lineStyles[id_style]);
            hDashedSyst->SetLineWidth(lineWidth);
            hDashedSyst->SetFillColor(0);
	    if(name.Contains("_vs_TTBarMass") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets")) {
	      if (hDashedSyst->GetMaximum() > hrange->GetMaximum()) hrange->SetMaximum(ceil(hDashedSyst->GetMaximum()*100)/100 + 0.05);
	    }
	    else {
	      if (hDashedSyst->GetMaximum() > hrange->GetMaximum()) hrange->SetMaximum(ceil(1.2*hDashedSyst->GetMaximum()*100)/100);
	    }
	    if (hDashedSyst->GetMinimum() < hrange->GetMinimum()) hrange->SetMinimum(floor(1.2*hDashedSyst->GetMinimum()*100)/100);
            hDashedSyst->Draw("hist0 same");
            if (nvar==1) leg->AddEntry(hDashedSyst, leg_mod_syst, "l");

            nvar++;
        }

        if (i_color > step_for_style) {
            i_color = 1;
            id_style++;
        } else i_color++;
    }

    drawCMSLabels(cmsPrelimLabelMode);
    //    drawDecayChLabel(channelLabel[channelType]);

    cc->cd();
    p1->Draw();

    p2->cd();
    leg->Draw();
    cc->cd();
    p2->Draw();
    
    cc->Draw();

    utils::ttmd::SaveCanvas(cc, out_base_path);
    cc->Print(this->outpathPlots + "/Nominal/" + channel_ + "/" + name + "_uncBreakdown" + ".C");
}


void PlotterCP::removeNegativeExpectations(TH1* thisHisto)
{
    // Remove negative expectations in bins for samples with negative weights (e.g. aMC@NLO)
    // Needed to avoid negative expectations for a process and exclude its impact on other processes in the histogram stack
    // Has to be used only after the binning is finalized
    // Recommended to be used only for backgrounds (instead, when using for ttbar signal, ensure to have enough statistics per bin of differential measurement; for now, ttbar from aMC@NLO used only for comparisons)

    for (int ibin = 0; ibin <= thisHisto->GetNbinsX(); ++ibin) {
        //Check visible bins + underflow
        if(thisHisto->GetBinContent(ibin) < 0.0) {
           thisHisto->SetBinContent(ibin, 0.0);// set prediction to 0.0
           thisHisto->SetBinError(ibin, 1.0);// set error to 1 in order to match expectation of +/-1 event
        }
    }

    if(thisHisto->GetBinContent(thisHisto->GetNbinsX()+1) < 0.0) {
        thisHisto->SetBinContent(thisHisto->GetNbinsX()+1, 0.0);// set prediction to 0.0
        thisHisto->SetBinError(thisHisto->GetNbinsX()+1, 1.0);// set error to 1 in order to match expectation of +/-1 event
    }
}


void PlotterCP::performPoissonSmearing(TH1 *hist)
{
    int seed = 4357;

    TRandom3 * randClosureTest = new TRandom3(seed);
    for (int binIndex = 0; binIndex <= hist->GetNbinsX()+1; ++binIndex) {// "0" and "+1" needed for under- and overflow bins
        double binContent = hist->GetBinContent(binIndex);
        double newContent = randClosureTest->PoissonD(binContent);
        hist->SetBinContent(binIndex, newContent);
    }
    delete randClosureTest;
}


void PlotterCP::addAndDelete_or_Assign(TH1*& addToThis, TH1* addThis)
{
    if (!addToThis) addToThis = addThis;
    else {
        addToThis->Add(addThis);
        delete addThis;
    }
}


void PlotterCP::setStyle(TH1 *hist, unsigned int i, bool isControlPlot)
{
    hist->SetFillColor(colors[i]);
    hist->SetLineColor(colors[i]);
    hist->SetLineWidth(1);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetTitleOffset(1.08);
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->GetXaxis()->SetLabelOffset(0.007);
    if (name.Contains("HypBjetMulti_noBTag") || name.Contains("basic_jet_multiplicity_step4") || name.Contains("HypJetMultpt30")) hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelOffset(0.007);

    // Jason

    if(name.Contains("_vs_TTBarMass")) hist->GetYaxis()->SetNdivisions(504); 
    else if(name.Contains("_vs_ScatteringAngle_TTBarFrame")) hist->GetYaxis()->SetNdivisions(504); 
    else if(name.Contains("_vs_ToppT")) hist->GetYaxis()->SetNdivisions(504); 
    
    else if(name.Contains("_TTBarMass")) hist->GetYaxis()->SetNdivisions(504); 
    else if(name.Contains("_TTBarpT")) hist->GetYaxis()->SetNdivisions(504); 
    else if(name.Contains("_TTBarRapidity")) hist->GetYaxis()->SetNdivisions(504); 

    else if(name.Contains("_AntiLepton")) hist->GetYaxis()->SetNdivisions(504); 
    else if(name.Contains("_Lepton")) hist->GetYaxis()->SetNdivisions(504); 
    else if(name.Contains("_LLBar")) hist->GetYaxis()->SetNdivisions(504); 


    if (legends.at(i) == "Data") {
        hist->SetFillColor(0);
        hist->SetMarkerStyle(20);
        hist->SetMarkerSize(1.2);
        hist->SetLineWidth(2);
	if ((name.Contains("pT", TString::kIgnoreCase) || name.Contains("Mass", TString::kIgnoreCase)) &&
	    (!name.Contains("1st") && 
	     !name.Contains("Rapidity", TString::kIgnoreCase) && 
	     !name.Contains("Eta", TString::kIgnoreCase) && 
	     !name.Contains("Phi", TString::kIgnoreCase) && 
	     !name.Contains("JetMult") && 
             !name.Contains("Fraction") && 
             !name.Contains("Multiplicity", TString::kIgnoreCase) && 

             !name.Contains("_AntiLepton", TString::kIgnoreCase) && 
             !name.Contains("_Lepton", TString::kIgnoreCase) && 
             !name.Contains("_LLBar", TString::kIgnoreCase) && 

             !name.Contains("_vs_ScatteringAngle", TString::kIgnoreCase ) && 
             !name.Contains("_vs_TTBarMass", TString::kIgnoreCase ) && 
             !name.Contains("_vs_ToppT", TString::kIgnoreCase ))
	    )
        {
	    hist->GetXaxis()->SetTitle(XAxis+" #left[GeV#right]");
        } else {
            hist->GetXaxis()->SetTitle(XAxis);
        }
        if (isControlPlot) hist->GetYaxis()->SetTitle(YAxis);
    }
}


void PlotterCP::setControlPlotLegendStyle(std::vector< TH1* > drawhists, std::vector< TString > legends, TLegend* leg, TLegend *leg1, TLegend *leg2, TLegend *leg3)
{
    //this appendix  to the legend aligns vertically all entries in the legend
    std::string  appendix = "                                                                                    ";

    //hardcoded controlplot legend
    std::vector<TString> OrderedLegends;
    TString pseudoData13TeVLeg = "Ps.-Data";
    OrderedLegends.push_back("Data");
    OrderedLegends.push_back("t#bar{t} signal");
    OrderedLegends.push_back("t#bar{t} other");
    OrderedLegends.push_back("Single t");
    OrderedLegends.push_back("W+jets");
    OrderedLegends.push_back("Z+jets");
    OrderedLegends.push_back("t#bar{t}+Z/W");
    OrderedLegends.push_back("Diboson");
    if (this->addQCDToControlPlot()) OrderedLegends.push_back("QCD multijet");

    leg->Clear();
    if (leg1) leg1->Clear();
    if (leg2) leg2->Clear();
    if (leg3) leg3->Clear();
    for (size_t i = 0; i < OrderedLegends.size(); ++i) {
        for (size_t j = 0; j < drawhists.size(); ++j) {
            if (OrderedLegends.at(i) == legends.at(j)) {
                if (OrderedLegends.at(i) == "Data") {
                    if (usePoissonSmearedPseudoData) {
                        leg->AddEntry(drawhists.at(j), pseudoData13TeVLeg+appendix, "pe");
                        if (leg1) leg1->AddEntry(drawhists.at(j), pseudoData13TeVLeg + appendix, "pe");
                    } else {
                        leg->AddEntry(drawhists.at(j), OrderedLegends.at(i)+appendix, "pe");
                        if (leg1) leg1->AddEntry(drawhists.at(j), OrderedLegends.at(i) + appendix, "pe");
                    }
                    break;
                } else {
                    leg->AddEntry(drawhists.at(j), OrderedLegends.at(i) + appendix, "f");
                    if (leg1 && i < 3) leg1->AddEntry(drawhists.at(j), OrderedLegends.at(i) + appendix, "f");
                    else if (leg2 && i >= 3 && i < 6) leg2->AddEntry(drawhists.at(j), OrderedLegends.at(i) + appendix, "f");
                    else if (leg3 && i >= 6) leg3->AddEntry(drawhists.at(j), OrderedLegends.at(i) + appendix, "f");
                    break;
                }
            }
        }
    }
    //coordinates for legend without splitting
    float x1 = 1.00 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.25 + 0.005;
    float y1 = 1.00 - gStyle->GetPadTopMargin() - gStyle->GetTickLength() - 0.05 - 0.03 * leg->GetNRows();
    float x2 = 1.00 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.025;
    float y2 = 1.00 - gStyle->GetPadTopMargin() - 0.8 * gStyle->GetTickLength();

    leg->SetX1NDC(x1);
    leg->SetY1NDC(y1);
    leg->SetX2NDC(x2);
    leg->SetY2NDC(y2);

    leg->SetTextFont(42);
    leg->SetTextSize(0.0244);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextAlign(12);

    if (!leg1) return;
    leg1->Draw("same");
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.0244);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextAlign(12);

    // Define shifts of legends
    float xShift = 0.8 * (x2 - x1);
    float yShift = 0.03;
    float y1mod = y1 + 0.45 * (y2 - y1);

    // coordinates for splitted legend 1
    gPad->Update();
    //    leg1->SetX1NDC(x1 - xShift);
    //    leg1->SetY1NDC(y1mod - yShift);
    //    leg1->SetX2NDC(x2 - xShift);
    //    leg1->SetY2NDC(y2 - yShift);
    leg1->SetX1NDC(0.397);
    leg1->SetY1NDC(0.7665217);
    leg1->SetX2NDC(0.567);
    leg1->SetY2NDC(0.8865217);

    if (!leg2) return;
    leg2->Draw("same");
    gPad->Update();
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.0244);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetTextAlign(12);

    // coordinates for splitted legend 2
    //    leg2->SetX1NDC(x1 + 0.01);
    //    leg2->SetY1NDC(y1mod - yShift);
    //    leg2->SetX2NDC(x2 + 0.01);
    //    leg2->SetY2NDC(y2 - yShift);

    leg2->SetX1NDC(0.593);
    leg2->SetY1NDC(0.7665217);
    leg2->SetX2NDC(0.763);
    leg2->SetY2NDC(0.8865217);


    if (!leg3) return;
    leg3->Draw("same");
    gPad->Update();
    leg3->SetTextFont(42);
    leg3->SetTextSize(0.0244);
    leg3->SetFillStyle(0);
    leg3->SetBorderSize(0);
    leg3->SetTextAlign(12);

    // coordinates for splitted legend 3
    //    leg3->SetX1NDC(x1 + 0.01);
    //    leg3->SetY1NDC(y1mod - yShift);
    //    leg3->SetX2NDC(x2 + 0.01);
    //    leg3->SetY2NDC(y2 - yShift);
    leg3->SetX1NDC(0.779);
    leg3->SetY1NDC(0.8065217);
    leg3->SetX2NDC(0.959);
    leg3->SetY2NDC(0.8865217);


}


void PlotterCP::drawDecayChLabel(TString decaychannel, double textSize)
{
    TPaveText *decch = new TPaveText();
    decch->AddText(decaychannel);
    decch->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
    decch->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );
    decch->SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
    decch->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength()        );
    decch->SetFillStyle(0);
    decch->SetBorderSize(0);
    if (textSize != 0) decch->SetTextSize(textSize);
    decch->SetTextAlign(12);
    if (decaychannel != "Dilepton") decch->Draw("same");
}


void PlotterCP::drawCMSLabels(int cmsprelim, double textSize)
{

    const char *text = (year == "fullRun2UL") ? "%.0f fb^{-1} (%2.f TeV)" : "%2.1f fb^{-1} (%2.f TeV)";

    TPaveText *label = new TPaveText();
    label->SetX1NDC(gStyle->GetPadLeftMargin() + 0.03);
    label->SetY1NDC(0.98 - gStyle->GetPadTopMargin());
    label->SetX2NDC(1.0 - gStyle->GetPadRightMargin() + 0.03);
    label->SetY2NDC(0.98);
    label->SetTextFont(42);
    label->AddText(Form(text, lumi/1000, energy));
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    if (textSize!=0) label->SetTextSize(0.045);
    label->SetTextAlign(32);
    label->Draw("same");

    TPaveText *cms = new TPaveText();
    cms->AddText("CMS");
    cms->SetX1NDC(0.437 - gStyle->GetPadLeftMargin());
    cms->SetY1NDC(0.98 - gStyle->GetPadTopMargin());
    cms->SetX2NDC(1.0 - gStyle->GetPadRightMargin());
    cms->SetY2NDC(0.98);
    cms->SetFillStyle(0);
    cms->SetBorderSize(0);
    if (textSize != 0) cms->SetTextSize(0.050);
    cms->SetTextAlign(12);
    cms->SetTextFont(61);
    cms->Draw("same");

    if (cmsprelim > 0) {
        TPaveText *extra = new TPaveText();
        if (cmsprelim == 2) extra->AddText("Private Work");
        else extra->AddText("Preliminary");
        extra->SetX1NDC(0.46 - gStyle->GetPadLeftMargin() + gStyle->GetTickLength());
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


void PlotterCP::draw2DLabels(double y1_NDC)
{

  TString hist_name = name;
  TString dimlabeltext = "";
  TString dimlabel_1text = "";
  TString dimlabel_2text = "";
  TString dimlabel_3text = "";
  TString dimlabel_4text = "";

  if (hist_name.Contains("_step8_0to450_ttbarmass")) dimlabeltext = "0 < m_{t#bar{t}} #leq 450";
  else if (hist_name.Contains("_step8_450to550_ttbarmass")) dimlabeltext = "450 < m_{t#bar{t}} #leq 550";
  else if (hist_name.Contains("_step8_550to800_ttbarmass")) dimlabeltext = "550 < m_{t#bar{t}} #leq 800"; 
  else if (hist_name.Contains("_step8_800toinf_ttbarmass")) dimlabeltext = "m_{t#bar{t}} > 800"; 
  else if (hist_name.Contains("_step8_m1tomhalf_topscatteringangle")) dimlabeltext = "-1.0 #leq cos#Theta #leq -0.5";
  else if (hist_name.Contains("_step8_mhalfto0_topscatteringangle")) dimlabeltext = "-0.5 < cos#Theta #leq 0.0";
  else if (hist_name.Contains("_step8_0tophalf_topscatteringangle")) dimlabeltext = "0.0 < cos#Theta #leq 0.5"; 
  else if (hist_name.Contains("_step8_phalftop1_topscatteringangle")) dimlabeltext = "0.5 < cos#Theta #leq 1.0";

  else if (hist_name.Contains("_vs_TTBarMass")) {
    dimlabel_1text = "0 < m_{t#bar{t}} #leq 450";
    dimlabel_2text = "450 < m_{t#bar{t}} #leq 600";
    dimlabel_3text = "600 < m_{t#bar{t}} #leq 800"; 
    dimlabel_4text = "m_{t#bar{t}} > 800"; 
  }
  else if (hist_name.Contains("_vs_ScatteringAngle_TTBarFrame")) {
    dimlabel_1text = "-1.0 #leq cos#Theta #leq -0.5";
    dimlabel_2text = "-0.5 < cos#Theta #leq 0.0";
    dimlabel_3text = "0.0 < cos#Theta #leq 0.5"; 
    dimlabel_4text = "0.5 < cos#Theta #leq 1.0";
  }
  else if (hist_name.Contains("_vs_ToppT")) {
    dimlabel_1text = "0 < p_{T}^{t} #leq 80";
    dimlabel_2text = "80 < p_{T}^{t} #leq 150";
    dimlabel_3text = "150 < p_{T}^{t} #leq 250"; 
    dimlabel_4text = "p_{T}^{t} > 250"; 
  }
  else if (hist_name.Contains("_vs_ExtraJets")) {
    dimlabel_1text = "N_{jets} = 2";
    dimlabel_2text = "N_{jets} = 3";
    dimlabel_3text = "N_{jets} = 4"; 
    dimlabel_4text = "N_{jets} #geq 5"; 
  }



  TPaveText *dimlabel = new TPaveText();
  dimlabel->AddText(dimlabeltext);
  dimlabel->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() );
  dimlabel->SetY1NDC(1.0 - 1.2*gStyle->GetPadTopMargin()  - 2*gStyle->GetTickLength() );
  dimlabel->SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
  dimlabel->SetY2NDC(1.0 - 1.2*gStyle->GetPadTopMargin()  - 2*gStyle->GetTickLength() - 0.05 );
  dimlabel->SetFillStyle(0);
  dimlabel->SetBorderSize(0);
  dimlabel->SetTextSize(0.040);
  dimlabel->SetTextAlign(22);
  dimlabel->Draw("same");

  TPaveText *dimlabel_1 = new TPaveText();
  dimlabel_1->AddText(dimlabel_1text);
  dimlabel_1->SetX1NDC(0.18);
  dimlabel_1->SetY1NDC(y1_NDC);
  dimlabel_1->SetX2NDC(0.3725);
  dimlabel_1->SetY2NDC(y1_NDC + 0.05);
  dimlabel_1->SetFillStyle(0);
  dimlabel_1->SetBorderSize(0);
  dimlabel_1->SetTextSize(0.025);
  dimlabel_1->SetTextFont(62);
  dimlabel_1->SetTextAlign(22);
  dimlabel_1->Draw("same");

  TPaveText *dimlabel_2 = new TPaveText();
  dimlabel_2->AddText(dimlabel_2text);
  dimlabel_2->SetX1NDC(0.3725);
  dimlabel_2->SetY1NDC(y1_NDC);
  dimlabel_2->SetX2NDC(0.565);
  dimlabel_2->SetY2NDC(y1_NDC + 0.05);
  dimlabel_2->SetFillStyle(0);
  dimlabel_2->SetBorderSize(0);
  dimlabel_2->SetTextSize(0.025);
  dimlabel_2->SetTextFont(62);  
  dimlabel_2->SetTextAlign(22);
  dimlabel_2->Draw("same");

  TPaveText *dimlabel_3 = new TPaveText();
  dimlabel_3->AddText(dimlabel_3text);
  dimlabel_3->SetX1NDC(0.565);
  dimlabel_3->SetY1NDC(y1_NDC);
  dimlabel_3->SetX2NDC(0.7575);
  dimlabel_3->SetY2NDC(y1_NDC + 0.05);
  dimlabel_3->SetFillStyle(0);
  dimlabel_3->SetBorderSize(0);
  dimlabel_3->SetTextSize(0.025);
  dimlabel_3->SetTextFont(62);
  dimlabel_3->SetTextAlign(22);
  dimlabel_3->Draw("same");

  TPaveText *dimlabel_4 = new TPaveText();
  dimlabel_4->AddText(dimlabel_4text);
  dimlabel_4->SetX1NDC(0.7575);
  dimlabel_4->SetY1NDC(y1_NDC);
  dimlabel_4->SetX2NDC(0.95);
  dimlabel_4->SetY2NDC(y1_NDC + 0.05);
  dimlabel_4->SetFillStyle(0);
  dimlabel_4->SetBorderSize(0);
  dimlabel_4->SetTextSize(0.025);
  dimlabel_4->SetTextFont(62);
  dimlabel_4->SetTextAlign(22);
  dimlabel_4->Draw("same");

}

void PlotterCP::yRangeControlPlotRatio(double &yminCP_, double &ymaxCP_) const
{
  //    yminCP_ = 0.841;  // 0.9 to 1.1
  //    ymaxCP_ = 1.159;  

  //    yminCP_ = 0.682;  // 0.8 to 1.2
  //    ymaxCP_ = 1.318;

    yminCP_ = 0.682;
    ymaxCP_ = 1.318;

    if (name.Contains("vertMulti"))    {yminCP_ = 0.5;    ymaxCP_ = 1.5;}


    if (name.Contains("LeptonB") || name.Contains("LLBarC"))    {yminCP_ = 0.841;    ymaxCP_ = 1.159;}
}

