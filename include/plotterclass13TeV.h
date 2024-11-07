#ifndef plotterclass13TeV_h
#define plotterclass13TeV_h

#include <vector>
#include <set>

class TGraphErrors;
class TLegend;
class TGraphAsymmErrors;
class TH1;
class TString;
class TH1F;
class TH1D;

class RootFileReader;
class UsefulTools13TeV;





class Plotter13TeV {

public:
    Plotter13TeV();
    void   setOptions(TString name_, TString specialComment_, TString YAxis_, TString XAxis_, int rebin_, bool doDYScale_, bool logX_, bool logY_, double ymin_, double ymax_, double rangemin_, double rangemax_, int bins_, std::vector<double> XAxisbins_, std::vector<double> XAxisbinCenters_);
    void   setDataSet(std::vector<TString>, std::vector<double>, std::vector<TString>, std::vector<int>, TString);
    void   setDataSet(TString, TString);
    bool   fillHisto();
    void   setStyle(TH1*, unsigned int, bool = false);
    void   unfolding(TString channel, TString systematic, const TString& unfoldingType = "svd"); // OZ 4.01.2017 modified for TUnfold
    void   preunfolding(TString Channel="", TString Systematic="", const TString& unfoldingType = "svd"); // OZ 4.01.2017 modified for TUnfold
    void   DoFitInRatio(bool doFit = 0);
    void   setDrawUncBand(bool drawUncBand);
    void   setMergeEnvelopeVariations(bool mergeEnvelopeVariations);
    void   SetTopMassSetup(){useTopMassSetup_=true;};
    void   OnlyCPSetup(){isCP_=true;};
    
    ///add addThis to addToThis (and delete it afterwards) - or assign it it to addToThis if addToThis is nullptr.
    void   addAndDelete_or_Assign(TH1*& addToThis, TH1* addThis, bool removeNegExpectation = 0, bool isUnfolding = 0);
    void   removeNegativeExpectationsInBins(TH1*& thisHisto);
    void   removeNegativeExpectationsInBinsDouble(TH1D*& thisHisto);
    void   write(TString, TString, const TString& unfoldingType = "svd"); // OZ 4.01.2017 modified for TUnfold
    
    double CalcXSec(std::vector<TString> datasetVec, std::vector<TString> pseudoDatasetVec , double InclusiveXsectionVec[4],double InclusiveXsectionStatErrorVec[4], TString Systematic, TString Shift);
    void MakeTable(TString Channel, TString Systematic);


    void PlotXSec(TString Channel);
    //  void CalcDiffXSec(TH1* varhists[], TH1* RecoPlot, TH1* GenPlot, TH2* genReco2d, double DiffXSecVec[4][10], double DiffXSecStatErrorVec[4][10]); 
    int CalcDiffXSec(TString, TString, const TString& unfoldingType = "svd"); // OZ 4.01.2017 modified for TUnfold
    void PlotDiffXSec(TString, std::vector<TString>);
    void PlotSingleDiffXSec(TString, TString);

    void DYScaleFactor(TString);

    void PrintResultTotxtFile(TString, double[], TGraphAsymmErrors *, TGraphAsymmErrors *);
    void PrintResultPlusGenTotxtFile (TString , double[], TGraphAsymmErrors *, TGraphAsymmErrors *, TH1 *, TH1 *,TH1 *,TH1 * );
    void CalcTotCovMtrx ( TString Channel, TString Variable, std::vector<TString>vec_systematic, double[], double[], double[], double[], double[]);
    void fillUnfMtrxAndTabsToFile (std::ofstream& ifile, const std::vector<std::vector<double>>& matrix, int bins, TString Type, TString Variable, TString Channel, bool latex, bool noverbosity = 0);
    void CalcUpDownDifference ( TString Channel, TString Syst_Up, TString Syst_Down, TString Variable);
    void GetEnvelopeForUnfoldedResults ( TString Channel, TString SystName, TString Variable);

    TH1* GetNloCurve(const char *particle, const char *quantity, const char *generator);
    TH1* GetNloCurve(TString NewName, TString Generator);
    TH1* GetNloCurveMass(TString NewName, TString Generator, TString Mass);
    TH1* RescaleNLOCurveAndPrepareBinnedCurve(TH1 *inputhist , int bins, double Xbins[], double xsecnorm, double BrFrac);
    TH1F* ConvertGraphToHisto(TGraphErrors *pGraph);
    double GetChi2 (TGraphAsymmErrors *data, TH1 *mc);
    TH1F* reBinTH1FIrregularNewBinning(TH1F *histoOldBinning, TString plotname, bool rescale);

    void UnfoldingOptions(bool doSVD);
    void SetOutpath(TString path); 
    void DrawLabel(TString text, const double x1, const double y1, const double x2, const double y2, int centering, double textSize);

    double CalculateIntegral(TGraphAsymmErrors *tga_DiffXSecPlot, double Xbins[]);
    

private:

    TString name;

    TString specialComment;
    int bins, datafiles, rebin;
    double rangemin, rangemax, ymin, ymax;

    std::vector<TString> dataset, datasetUp, datasetDown;
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

    std::vector<TH1D> hists;
    std::vector<TH1D> histsTtBgrOutOfSpace;
    std::vector<TH1D> histsForPseudoData;
    std::vector<TH1D> systhistsUp;
    std::vector<TH1D> systhistsDown;

    bool doFit_;
    bool initialized, logX, logY, doDYScale;
    double energy, lumi, topxsec;
    int signalHist;
    
    std::vector<TString> channelLabel;

    double SignalEventswithWeight;
    
    bool doUnfolding; 
    bool doSystematics;
    bool drawSmoothMadgraph, drawPlotRatio;
    bool drawNLOCurves, drawMadSpinCorr, drawMCATNLO, drawKidonakis, drawAhrens;
    bool drawPOWHEG, drawPOWHEGHERWIG, drawPERUGIA11;
    bool drawMadMass, drawMadScaleMatching;

    // Booleans and strings configuring diff. cross section measurements
    bool absoluteXsec;  
    bool fullPhaseXsec;
    bool particleXsec;
    TString visTheoryPrefix;
    TString respMatrixPrefix;
  
    // For estimation of syst. uncertainty from pseudo-data
    bool estimateSysUncFromPseudoData;
    double lumiError;
    double BRError;
    double massUncReductionFactor;
    double scaleUncReductionFactor;
    TString outpath;
    TString outpathPlots;
    TString outpathResults;
    TString subfolderChannel;
    TString subfolderSpecial;

    static const bool doClosureTest;
    static const bool usePoissonSmearedPseudoData;
    RootFileReader *fileReader;
    UsefulTools13TeV *usefulTools13TeV;
    void DrawDecayChLabel(TString decaychannel="", double textSize=0.04);
    void DrawCMSLabels(int cmsprelim=true, double textSize=0.04);
    // Mode defining plotting of preliminary labels, to be used in DrawCMSLabels() 
    int cmsPrelimLabelMode;
    
    /// Define members and enums for theory curves style properties
    void setTheoryStyleAndFillLegend(TH1 *histo, TString theoryName, TLegend *leg = 0);

    /// Set style of result and control plot legend
    void setResultLegendStyle(TLegend *leg, const bool result = 1);
    void setControlPlotLegendStyle(std::vector<TH1*> drawhists, std::vector<TString> legends, TLegend *leg, TLegend *leg1 = nullptr, TLegend *leg2 = nullptr);

    /// boolean to decide to add or not the QCD background to the control plot
    bool addQCDToControlPlot()const;

    /// Derive tt signal model unceratinty band
    void getSignalUncertaintyBand(TH1 *uncBand, TString channel_);
    void GetEnvelopeVariationsForBandsInCPs( TH1 *varHist , const TString Systematic , const TString channel );
    bool drawUncBand_;
    bool mergeEnvelopeVariationsForBandsInCPs_;
    bool mergeScaleVariations_;
    bool mergeBTagVariations_;
    bool mergeTriggerVariations_;
    bool mergeBFragVariations_;
    bool mergeColorRecVariations_;
    
    /// Perform Poisson smearing of the histogram
    void performPoissonSmearing(TH1 *hist);

    /// Control plot ratio yaxis rangemax
    void yRangeControlPlotRatio(double &yminCP_, double &ymaxCP_)const;

    /// Set yaxis range of ratio plot of the result
    void setResultRatioRanges(double &ymin, double &ymax)const;

    // optional setup for top mass extraction analysis
    bool useTopMassSetup_;

    // communicate use of CP plotting only
    bool isCP_;
};

#endif
