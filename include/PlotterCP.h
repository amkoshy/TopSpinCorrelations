#ifndef PlotterCP_h
#define PlotterCP_h

#include <vector>
#include <set>

class TLegend;
class TString;
class TH1;
class TH1F;
class TH1D;
class TGraphErrors;
class TGraphAsymmErrors;

class RootFileReader;
class PlotterConfigurationHelper;



class PlotterCP
{
  public:

    PlotterCP(TString analysis_year, TString treat_correlated_year);
    void setOptions(TString name_, TString specialComment_, TString YAxis_, TString XAxis_, int rebin_, bool doDYScale_, bool logX_, bool logY_, double ymin_, double ymax_, double rangemin_,
                    double rangemax_, int bins_, std::vector<double> XAxisbins_, std::vector<double> XAxisbinCenters_);
    void setDataSet(std::vector<TString>, std::vector<double>, std::vector<TString>, std::vector<int>, TString);
    void setDataSet(TString, TString);
    bool fillHisto(TString Channel, TString Systematic);
    void DYScaleFactor(TString);
    void plotWriterWrapper(TString Channel="", TString Systematic="");
    void writeCP(TString, TString);
    double CalcXSec(std::vector<TString> datasetVec, std::vector<TString> pseudoDatasetVec , double InclusiveXsectionVec[4],double InclusiveXsectionStatErrorVec[4], TString Systematic, TString Shift);
    void PlotXSec(TString Channel);
    void MakeTable(TString Channel, TString Systematic);
    void setStyle(TH1*, unsigned int, bool = false);
    void setControlPlotLegendStyle(std::vector<TH1*> drawhists, std::vector<TString> legends, TLegend *leg, TLegend *leg1 = nullptr, TLegend *leg2 = nullptr, TLegend *leg3 = nullptr);

    void DoFitInRatio(bool doFit) {doFit_ = doFit;}
    void setMergeEnvelopeVariations(bool mergeEnvelopeVariations) {mergeEnvelopeVariationsForBandsInCPs_ = mergeEnvelopeVariations;}
    void setUseEnvelopesInSystematicsList(bool useEnvelopesInSystematicsList) {useEnvelopesInSystematicsList_ = useEnvelopesInSystematicsList;}
    void setOutpath(TString path) {outpath = path;}
    void setDrawUncBand(bool drawUncBand) {drawUncBand_ = drawUncBand;}
    void OnlyCPSetup(){isCP_=true;};
    void SetOutpath(TString path);

  private:

    /// boolean to decide to add or not the QCD background to the control plot
    bool addQCDToControlPlot()const;

    void GetEnvelopeVariationsForBandsInCPs( TH1 *varHist , const TString Systematic , const TString channel );

    /// Derive tt signal model unceratinty band
    void getSignalUncertaintyBand(TH1 *uncBand, std::map<TString, std::vector<double> > &varUnc, TString channel_);
    void plotUncertaintyBandContributions(std::map<TString, std::vector<double> > &varUnc, const TString channel_, TH1* refHisto, const std::vector<TString>& syst_list = {}, const TString& fileNameSufix="");


    void removeNegativeExpectations(TH1* thisHisto);

    /// Perform Poisson smearing of the histogram
    void performPoissonSmearing(TH1 *hist);

    ///add addThis to addToThis (and delete it afterwards) - or assign it it to addToThis if addToThis is nullptr.
    void   addAndDelete_or_Assign(TH1*& addToThis, TH1* addThis);

    void drawDecayChLabel(TString decaychannel="", double textSize=0.04);

    void drawCMSLabels(int cmsprelim=true, double textSize=0.04);

    void draw2DLabels(double y1_NDC = 0.71);

    /// Control plot ratio yaxis rangemax
    void yRangeControlPlotRatio(double &yminCP_, double &ymaxCP_)const;

    TString year;
    TString yearCorr;
    TString name;
    TString specialComment;
    int bins, datafiles, rebin;
    double rangemin, rangemax, ymin, ymax;

    std::vector<TString> dataset;
    std::vector<TString> pseudoDataset;
    std::vector<double> scales;
    std::vector<TString> legends, legendsSyst;
    std::vector<int> colors;
    std::vector<double> XAxisbins, XAxisbinCenters;

    std::vector<double> DYScale;
    TString DYEntry;
    TString YAxis;
    TString XAxis;
    TString channel;
    int channelType; //0=ee 1==mumu 2==emu 3==combined  
    bool initialized, logX, logY, doDYScale;
    std::vector<TString> channelLabel;

    std::vector<TH1D> hists;
    std::vector<TH1D> histsForPseudoData;

    RootFileReader *fileReader;
    PlotterConfigurationHelper *configHelper;
    double energy, lumi, lumiError, topxsec, bRatio[4], bRatioError;
    TString lumi_fb;
    std::vector<TString> uncertaintyBandSystematics, uncertaintyBandSystematicsShapeOnly;
    double massUncReductionFactor;
    double fsrScaleUncReductionFactor;
    bool drawPlotRatio;
    int cmsPrelimLabelMode;
    bool drawTotalMCErrorForCP;
    bool revokeAntiQuantityCombination;
    bool doFit_;
    bool drawUncBand_;
    bool mergeEnvelopeVariationsForBandsInCPs_, useEnvelopesInSystematicsList_;
    bool isCP_;

    bool usePoissonSmearedPseudoData;
    // For estimation of syst. uncertainty from pseudo-data
    bool estimateSysUncFromPseudoData;

    TString outpath;
    TString outpathPlots;
    TString subfolderChannel;
    TString subfolderSpecial;
};

#endif
