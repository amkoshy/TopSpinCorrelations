// 
// This file taken from Ievgen Korol's code (/nfs/dust/cms/user/korol/CMSSW_5_3_27) on 14.01.2017 as used for 2D analysis TOP-14-013
// with purpose to run TUnfold. Modifed to work for 1D.
//
// Modified to work in 1D: "fake" first dimension with 1 bin has to be provided, see example in Plotter13TeV::CalcDiffXSec()
//
// Oleksandr Zenaiev (oleksandr.zenaiev@desy.de), 14.01.2017
//
#ifndef PlotterForTUnfold_h
#define PlotterForTUnfold_h

#include <Rtypes.h>
#include <TRandom3.h>

class TCanvas;
class TLegend;
class TH1;
class TH2;
class Output;
class TDirectory;

#include "TUnfoldBinning.h"
#include "TUnfoldDensity.h"
#include "SamplesFwd.h"
#include "Samples.h"
#include "Sample.h"

class PlotterForTUnfold{
    
public:
    
    /// Constructor
    PlotterForTUnfold(const Samples& samples, double lumi, double lumiUnc);
    
    /// Destructor
    ~PlotterForTUnfold(){};
    
    /// Set options specific for variables according to list in NameList
    void setOptions(const std::vector<TString>);
    
    void setRegMethod(TString method);
    void setTauTest(bool isTau);
    void setCT(bool isCT);
    void setVar(double);
    
    /// Produce the plots for histogram under consideration from all samples
    void producePlots();
    //void printf(const char* arg1, const char* Data);
    
    // OZ 4.01.2017 get TUnfold results
    TH1* GetUnfoldedData() { return unfoldedData_; };
    TH1* GetUnfoldedDataNorm() { return unfoldedDataNorm_; };
    // OZ 15.01.2017 get fine binning
    std::vector<std::vector<Double_t> >  GetFineBinning() { return v_fineBins_; }

private:
    
    void clearMemory();
    void clearMemoryPerSystematicChannel();
    
    ///Prepare Histograms
     void prepareHistograms(const std::vector<Sample>& v_sample);
     
     /// Write canvas
     void writeCanvas(TCanvas* cnavas, const TString name);
     
     const Samples samples_;
     double lumi_ = -999;
     double lumi0_;
     const double lumiUnc_;
     
     double topxsec_;
     double nDataInc_;
     double nBgrInc_;    //backgrounds without tt-bar
     
     double nRecoTtAllInc_;
     double nRecoSigInc_;
     double nRecoAllMC_;
     double nRecoTtOtherInc_;
     double nGenSigInc_;
     double nRecoSingleTop_;
     double nRecoWPlusJets_;
     double nRecoZee_;
     double nRecoZtautau_;
     double nRecoTtPlusZWG_;
     double nRecoDiboson_;
     double nRecoQCD_;
     
     std::vector<TString> v_plotName_;
     std::vector<TString> v_plotTitle_;
     std::vector<TString> v_plotUnits_;
     std::vector<TString> v_plotYunitscp_;
     std::vector<int> v_plotYlog_;
     
     TString plotsFolder_;
     TString tauFolder_;
     TString ctFolder_;
     TString systematicName_;
     
     
     std::vector<TH1D* > v_histOnReco_;
     std::vector<TH1D* > v_histOnGen_;
     
    ///Control plots
    std::vector<int> v_cpNBins_;
    std::vector<Double_t> v_R1_;
    std::vector<Double_t> v_R2_;
    std::vector<Double_t> v_ymax_;
    std::vector<Int_t> v_leg2_;
    std::vector<Int_t> v_eban_;
    std::vector<std::vector<TH1D* > > vv_SampleHist_; // v-nD_ , v-sample , Hist 
    std::vector<TH1* >  v_histRecoAllBins_;
    std::vector<std::vector<std::vector<TH1D* > > > vvv_SampleHist_; // v-nD_ , v-bin , v-sample , Hist    
    std::vector<std::vector<TH1D* > > vv_UnderflowSampleHist_; //v-nD_ , v-sample , Hist 
    std::vector<std::vector<TH1D* > > vv_OverflowSampleHist_;
    void writePlotCP(const std::vector<Sample>& v_sample,const std::vector<TH1D* >& v_SampleHist,const int& ind, const int binNum = -999);
    void writePlotCPAllBins(const std::vector<Sample>& v_sample ,const std::vector<TH1* >& v_SampleHist);
    
    ///Unfolding Options
    TString regMethod_;
    bool isTauTest_;
    double tau_;
    double logTau_;
    bool isCT_;
  
    double var_;
    bool isVar_;
    
    
    ///Unfolding
    std::vector<std::vector<Double_t> >  v_coarseBins_;
    std::vector<std::vector<Double_t> >  v_fineBins_;
    std::vector<int> v_uReco_;
    std::vector<int> v_oReco_;
    std::vector<int> v_uTrue_;
    std::vector<int> v_oTrue_;
    TUnfoldBinning* detectorBinning_;
    TUnfoldBinning* generatorBinning_;
    TUnfoldBinning* detectorDistribution_;
    TUnfoldBinning* generatorDistribution_;
    TH2* histMigration_;
    void cleanMigrationMatrix(TH2*& histMigration);
    std::vector<double> v_covX_;
    std::vector<double> vX_;
    std::vector<double> v_covXnorm_;
    std::vector<double> vXnorm_;
    std::vector<double> vXnormMC_;
    std::vector<double> vXmc_;
    std::vector<double> v_covY_;
    std::vector<double> v_dY_;
    std::vector<double> v_covYnorm_;
    std::vector<double> v_dYnorm_;
    std::vector<double> vYmc_;
    std::vector<double> vYnormMC_;
    
    
    TH1* histBgrUo_;
    TH1* histBgrUoScaled_;
    TH1* histBgr_;
    TH1* histBgrScaled_;
    TH1* histData_;
    TH1* histDataScaled_;
    
    TH1* histFsigReco_;
    TH1* histTtAllReco_;
    TH1* histTtSigReco_;
    TH1* histAllReco_;
    TH1* histAllRecoScaled_;
    
    TH1* histBgrUoGenBin_;
    TH1* histBgrUoScaledGenBin_;
    TH1* histBgrGenBin_;
    TH1* histBgrScaledGenBin_;
    TH1* histDataGenBin_;
    TH1* histDataScaledGenBin_;
    
    TH1* histFsigRecoGenBin_;
    TH1* histTtAllRecoGenBin_;
    TH1* histTtSigRecoGenBin_;
    
    
    TH1* unfoldedData_;
    TH1* unfoldedDataNorm_;
    TH1* histUnfGenBin_;
    TH1* histTauLog_;
    TH1* histTauLogRhoScan_;
    TH1* histTauLogLCurve_;
    
    TH2* rhoIJtotal_;
    
    TRandom3* r3_;
    
    int genBin_(const std::vector<float>& val);
    int recoBin_(const std::vector<float>& val);
    void runUnfolding(const TH2* histMigration,const TH1* histInput,
                      const TH1* histBgr,const TH1* histBgrUo);
    
    ///Closure test
    TH1* histSR_;
    TH1* histSG_;
    void performPoissonSmearing(TH1* histIn,TH1*& histOut);
    std::vector<TH1* > v_histGenRewUp_;
    std::vector<TH1* > v_histGenRewDown_;
    std::vector<TH1* > v_histRecoRewUp_;
    std::vector<TH1* > v_histRecoRewDown_;
    
    TH1* histGenPtRewUp_;
    TH1* histGenPtRewDown_;
    TH1* histGenRapidityRewUp_;
    TH1* histGenRapidityRewDown_;
    TH1* histRecoPtRewUp_;
    TH1* histRecoPtRewDown_;
    TH1* histRecoRapidityRewUp_;
    TH1* histRecoRapidityRewDown_;
    double rewTopPtUp(double pt);
    double rewTopPtDown(double pt);
    double rewTopRapidityUp(double y);
    double rewTopRapidityDown(double y);
    
    TH1* histGen_;
    TH1* histGenNorm_;
    
    TH1* histPurityAllBins_;
    TH1* histStabilityAllBins_;
    std::vector<double> vectorPurityAllBins_;
    std::vector<double> vectorStabilityAllBins_;
    TH1* histRecoGenAllBins_;
    TH1* histEffAllBins_;
    
    void writePlotEPSAllBins();
    TH1D* plotNiceTitel(double yAxisRangeMin,double yAxisRangeMax,TString yTitle="");
    TH1D* plotNiceTitelReco(double yAxisRangeMin,double yAxisRangeMax);
    void writePlotXSec(const TH1* hDataNorm,const TH1* hMCNorm, const TH1* hData,const TH1* hMC, const int channel);
    std::vector<double> v_BR_;
    
    /// Tree branches
    Int_t entry_,entry0_;
    Float_t eventWeight_, trueLevelWeight_, trueLevelWeight0_;
    Int_t isKinReco_;
    bool isSignal_;
    std::vector<float> branchVals_;
    std::vector<float> branchValsGen_;
    std::vector<float> branchValsGen0_;
    
    
    bool isPDF_ ;
    
    void moveLegend(TLegend*& legend,double moveX,double moveY,double plusWx,double plusWy);

  public:
    // OZ 17.01.2017 1D to 2D (with 1 bin in 1st dim) histogram converter
    static TH2D* H1dto2d(const TH1D* hIn);

    // OZ 17.01.2017 2D (with 1 bin in 1st dim) to 1D histogram converter (overlow underflow in 1st dimension is lost)
    static TH1D* H2dto1d(const TH2D* hIn);
    
    // OZ 17.01.2017 set generator level histogram
    void SetHistGen(TH2D* histGen);
  
    // OZ 15.01.2017 method to determine fine binning automatically
    static std::vector<double> DetermineFineBinning(const TH1* hCoarse, const TH1* hInput, TString* message = NULL);

    // OZ 15.01.2017 check equality of float numbers with tolerance 
    static bool IsEqual(const double val1, const double val2, const double eps = 1e-6);
    
    // OZ 18.01.2017 public method to run TUnfold
    void RunUnfolding(const TH2* histMigration,const TH1* histInput,const TH1* histBgr,const TH1* histBgrUo)
    {
      runUnfolding(histMigration, histInput, histBgr, histBgrUo);
    }
    
    // OZ 18.01.2017 simplified arguments
    void writePlotXSec(const int channel) 
    {
      writePlotXSec(unfoldedDataNorm_, histGenNorm_, unfoldedData_, histGen_, channel);
    }

    // OZ 18.01.2017 set output folder
    void SetPlotsFolder(const Channel::Channel& channel, const Systematic::Systematic& systematic)
    {
      plotsFolder_ = common::assignFolder("TUnfold",channel,systematic,v_plotName_.at(0)+"-"+v_plotName_.at(1));
    }

    // OZ 15.01.2017 simplified unregularized unfolding, useful as cross check
    static TH1* DoSimpleUnfolding(const TH1* hData, const TH2* hResp);
};




#endif







