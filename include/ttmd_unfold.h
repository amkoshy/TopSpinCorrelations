#ifndef ZUnfold_h
#define ZUnfold_h

#include <TUnfoldDensity.h>

#include "PlotterConfigurationHelper.h"
class TString;
class TH1D;
class TH2D;
class UnfoldingXSecHandler;
class UnfoldingIOHandler;

class ZUnfold
{
  public:

    TString Remark;
    TString Suffix;
    int MarkerStyle;
    int Color;
    bool FlagTreatStatisticsSeriously = false;
    TUnfoldDensity* tunf = NULL;
    const UnfoldingIOHandler* _unfoldingIOHandler = NULL;
    TH1D* hInput = NULL;
    TH2D* hRsp = NULL;
    TH1D* hUnf = NULL;
    TH2D* hEmatTotal = NULL;
    TH2D* hCormatTotal = NULL;
    TH1D* hUnfN = NULL;
    TH2D* hEmatTotalN = NULL;
    TH2D* hCormatTotalN = NULL;
    TH1D* hEffBBB = NULL;
    TH1D* hGen = NULL;
    TH1D* hGenN = NULL;

    double Chi2Data;
    double Chi2Reg;
    int Dof;
    double Chi2;
    double Chi2N;
    double TauBest;

    int flagUFRec;
    int flagOFRec;
    int flagUFGen;
    int flagOFGen;

    TUnfold::EConstraint constraintMode;
    TUnfold::ERegMode regMode;
    TUnfoldDensity::EDensityMode densityFlags;

    // detailed steering for regularisation
    const char *REGDISTR;
    const char *REGAXISSTEER;

    enum ZUnfoldReg
    {
      RegNo,
      RegRhoAvg,
      RegRhoMax,
      RegLCurve
    };
    ZUnfoldReg regChoice;
    double tau;

    ZUnfold(PlotterConfigurationHelper *configHelper);

    int GetDataChi2Dim() const;
    void RunUnfolding(UnfoldingIOHandler* unfoldingIOHandler, const TString& ch, const TString& fileName = ""); // empty fileName no plots
    void PlotInput(const TString& fileName);
    void SetHistoStyle(TH1* h);

  private:

    PlotterConfigurationHelper *configHelper_;

    void RunTUnfold(const UnfoldingIOHandler* unfoldingIOHandler, const TString& fileNameBase = ""); // empty fileName no plots
    void RunBBB(UnfoldingIOHandler* unfoldingIOHandler, const int flagEff = 1);
    void MakeNormalisation();
    void CheckStatisticsIsGaussian(const TH1* h);
    void CalculateChi2();
    TH1D* GetEfficiency(UnfoldingIOHandler* unfoldingIOHandler, const  int flag = 1); // flag = 0: gen_i & rec / gen_i, flag = 1: rec_i / gen_i
};

#endif
