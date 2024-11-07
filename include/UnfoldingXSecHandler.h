#ifndef UnfoldingXSecHandler_h
#define UnfoldingXSecHandler_h

#include <vector>
#include <TUnfoldBinning.h>
#include <TF2.h>
#include <TF3.h>
#include "ttmd_vars.h"
#include "PlotterConfigurationHelper.h"
#include "UnfoldingIOHandler.h"

#include <exception>
#include <cassert>

class TF1;
class TString;
class TH1;
class TH1D;
class TH2;
class TH2D;
class TreePlain;
class ZMadGraph;
class ZPredSetup;


class UnfoldingXSecHandler
{
  private:

    PlotterConfigurationHelper *configHelper_;

    //bool flagApplyCut = false;

    virtual std::vector<double> GetPredictionImplPlain(const UnfoldingIOHandler* unfoldingIOHandler, ZMadGraph* mg, const ZPredSetup& predSetup) const;
    virtual std::vector<double> GetPredictionImplIntegrate(const UnfoldingIOHandler* unfoldingIOHandler, ZMadGraph* mg, const ZPredSetup& predSetup) const;
    virtual std::vector<double> GetPredictionImplPlainFast(const UnfoldingIOHandler* unfoldingIOHandler, ZMadGraph* mg, const ZPredSetup& predSetup) const;
    virtual std::vector<double> GetPredictionImplIntegrateFast(const UnfoldingIOHandler* unfoldingIOHandler, ZMadGraph* mg, const ZPredSetup& predSetup) const;

    static TH1D* ReadPredictionHistoSum(ZMadGraph* mg, const ZPredSetup& predSetup, const std::vector<int>& vHistoID, const int order);
    static float* ReadPredictionHistoSumFast(ZMadGraph* mg, const ZPredSetup& predSetup, const std::vector<int>& vHistoID, const int& oneBin, const int order, int* nbinsPtr = NULL);

  public:

    void CheckResponseNormalisation(const UnfoldingIOHandler* unfoldingIOHandler);

    // Rew. functions
    static double d1_st1(double* xx, double* par);
    static double d2_st1_mult(double* xx, double* par);
    static double d3_st1_mult(double* xx, double* par);
    static TF1* GetRewFunction(const UnfoldingIOHandler* unfoldingIOHandler, const std::vector<std::vector<double> >& vvPar);
    static TF1* GetRewFunction(const UnfoldingIOHandler* unfoldingIOHandler, const int grade);

    std::vector<double> GetPrediction(const UnfoldingIOHandler* unfoldingIOHandler, ZMadGraph* mg, const ZPredSetup& predSetup, const bool flagNorm = false, const bool flagDivideBW = false) const;
    TH1D* GetPredictionHisto(const UnfoldingIOHandler* unfoldingIOHandler, ZMadGraph* mg, const ZPredSetup& predSetup, const bool flagNorm = false, const bool flagDivideBW = false) const;
};

#endif
