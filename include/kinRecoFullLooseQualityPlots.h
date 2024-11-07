#ifndef kinRecoFullLooseQualityPlots_h
#define kinRecoFullLooseQualityPlots_h

#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TH2.h>
#include <TH1.h>
#include <Rtypes.h>
#include <TFile.h>
#include <iostream>
#include <vector>
#include <TGraphAsymmErrors.h>

class PrepareQualityPlots
{

 public:
  PrepareQualityPlots(const char *titleTrueVsDeltaTrueReco, const char *titleProjYTrueVsDeltaTrueReco, const char *titleGenVsRMS, const char *titleGenVsMean);
  ~PrepareQualityPlots(){}

  void prepareQualityPlots(const TString& fileDir, const TString& fileName, const TString& histoName, const int& Xnb, const float& Xr1, const float& Xr2);

  TH2D* getHistoTrueVsDeltaTrueReco(){return h_TrueVsDeltaTrueReco_;};
  TH1D* getHistoProjYTrueVsDeltaTrueReco(){return h_projYTrueVsDeltaTrueReco_;};
  TH1D* getHistoGenVsRMS(){return h_GenVsRMS_;};
  TH1D* getHistoGenVsMean(){return h_GenVsMean_;};

 private:

  const char *titleTrueVsDeltaTrueReco_;
  const char *titleProjYTrueVsDeltaTrueReco_;
  const char *titleGenVsRMS_;
  const char *titleGenVsMean_;

  TH2D *h_TrueVsDeltaTrueReco_;
  TH1D *h_projYTrueVsDeltaTrueReco_;
  TH1D *h_GenVsRMS_;
  TH1D *h_GenVsMean_;

};

class DrawQualityPlots
{

 public:
  DrawQualityPlots(const TString& plotsDir, const TString& plotsBaseName, const TString& ch, const TString& year);
  ~DrawQualityPlots(){}

  void drawQualityPlots(TH2D *h_TrueVsDeltaTrueReco, TH1D *h_projYTrueVsDeltaTrueReco, TH1D *h_GenVsRMS, TH1D *h_GenVsMean);
  TGraphAsymmErrors* drawAsGraph(const TH1D* h, const bool flagUnc = true, const double offset = 0.5, const TString& option = "ZP0") const;
  void rebin2D(TH1 *h, Int_t ngx, Int_t ngy) const;
  void calculatePSE(const float& rebin) const;
  void drawMigrationMatrix(const TString& varXTitle, const float& rebin) const;

  TH2D* getHistoRsp(){ return hResp_; };
  TH1D* getHistoGen(){ return hGen_; };
  TH2D* setHistoRsp(const TString& fileDir, const TString& fileName, const TString& histoName)
  { 
    TFile rootFile(fileDir + fileName);
    hResp_ = (TH2D*)rootFile.Get(histoName);
    return hResp_; 
  };
  TH1D* setHistoGen(const TString& fileDir, const TString& fileName, const TString& histoName)
  { 
    TFile rootFile(fileDir + fileName);
    hGen_ = (TH1D*)rootFile.Get(histoName);
    return hGen_; 
  };

 private:
  const TString plotsDir_;
  const TString plotsBaseName_;
  const TString channel_;
  const TString year_;
  TH2D *hResp_;
  TH1D *hGen_;

};

class CompareQualityPlots
{

 public:
  CompareQualityPlots(const TString& compareDir, const TString& compareBaseName, const TString& legend1, const TString& legend2, const TString& ch, const TString& year);
  ~CompareQualityPlots(){}

  void comparePlots(const TString& filePath1, const TString& histName1, const TString& filePath2, const TString& histName2);

 private:
  const TString compareDir_;
  const TString compareBaseName_;
  const TString channel_;
  const TString year_;
  const TString legend1_;
  const TString legend2_;

};

#endif
