#ifndef PlotterDiffXSec_h
#define PlotterDiffXSec_h

#include "PlotterConfigurationHelper.h"

#include <cassert>
#include <vector>

#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TTree.h>
#include <TLegend.h>
#include <TPolyLine.h>

#include "UnfoldingXSecHandler.h"
#include "UnfoldingIOHandler.h"
#include "ttmd_unfold.h"
#include "UnfoldingClosureTestsHandler.h"
#include "ttmd_measur.h"
#include "ttmd_ZMadGraph.h"

class UnfoldingXSecHandler;
class UnfoldingIOHandler;
class TString;
class TH1D;
class ZUnfold;

class PlotterDiffXSec
{
  private:

    bool zFlagDataFilled = false;
    bool zFlagMCFilled = false;
    bool zFlagMCBkgTTbarFilled = false;
    bool zFlagMCBkgFilled = false;
    TString zVarLast;
    TString year;
    PlotterConfigurationHelper *configHelper;

    static std::vector<std::vector<TH1D*> > CreateCPHistos(const UnfoldingIOHandler* unfoldingIOHandler, const int binning, const std::vector<TH1D*> vh);
    static std::vector<std::pair<TH1D*, TH1D*> > CreateLastDimUDHistograms(const UnfoldingIOHandler* unfoldingIOHandler, const std::pair<TH1D*, TH1D*>& pairUD, const int binning, const int flagNorm);
    static std::vector<std::vector<std::pair<TH1D*, TH1D*> > > CreateLastDimUDHistograms(const UnfoldingIOHandler* unfoldingIOHandler, const std::vector<std::pair<TH1D*, TH1D*> >& pairUD, const int binning, const int flagNorm);

    static void DrawBinLabel(const UnfoldingIOHandler* unfoldingIOHandler, const int i, const double axisX, const double histoX);
    static void DrawBinLabelFull(const UnfoldingIOHandler* unfoldingIOHandler, const int i, const double axisX, const double histoX, TString plot_type);
    static std::vector<TPolyLine*> CreatePolylineBins(const UnfoldingIOHandler* unfoldingIOHandler, const int nbins, const double minY, const double maxY);
    static void DrawCommonXTitle(const double x0, const double y0, const double x1, const double y1, const TH1* h);
    static void DrawBinLabelFullVertical(const UnfoldingIOHandler* unfoldingIOHandler, const int i/*, const double axisX, const double histoX*/);

    static void PrintHistBins(const TH1* h, FILE* f);
    static std::vector<float> ReadHistoBins(FILE* f);

    static void PrintXsecAllAlg(std::vector<ZUnfold*> vAlg, FILE* fout, const int typeXSec);

    static void PrintTexTabSystHeader(FILE* f, const int flagLong, const int nBins, const TString& strXSecLabel, const TString& strXSecLabel0, const TString& strXSec, const TString& strExtraBins);
    static void PrintTexTabSystFooter(FILE* f, const int flagLong);

  public:

    PlotterDiffXSec(TString analysis_year,PlotterConfigurationHelper *configHelper_, int mode = 1, const TString& analysisSuffix = "", int ModeKinRecoTree = 0);
    PlotterDiffXSec();
    ~PlotterDiffXSec(){};

    std::vector<UnfoldingIOHandler*> vXSec;
    std::vector<TString> vCh;
    TString current_channel;
    long MaxEvents;
    TString Var;
    //TString Dir;
    const int zMode;
    const TString AnalysisSuffix;
    const int ModeKinRecoTree; // 0 default full kin reco, 9 loose kin reco
    // bool FlagXFitterWrite = true;
    bool FlagSkipExtraPlots;

    std::vector<TString> VVarStoreCP = { "Nominal", "REWTOPPT", "MASS_UP", "MASS_DOWN", "MASS_CONSTMC",
                                         "MASS_KINRECO3GEV_DOWN", "MASS_KINRECO3GEV_UP", "MASS_KINRECO_DOWN", "MASS_KINRECO_UP" };
    std::vector<TString> VVarStorePSE = { "Nominal" };
    std::vector<TString> VVarStoreXSec = { "Nominal" };

    static constexpr int PlotPixelX = 1000;
    static constexpr int PlotPixelY = 750;
    static constexpr int NDiv = 305;
    static int FlagPaper; // 0 nothing, 1 paper, 2 pas
    static bool FlagAnalysisNote;

    struct Params
    {
        int Binning = 0;
        bool LogX = false;
        bool LogY = false;
        double RatMin = 0.65;
        double RatMax = 1.35;
        double YMin = -1.0;
        double YMax = -1.0;
        int DivideBW = 1;
    };

    // *****************************************************************************************
    // *************************** CONTROL PLOT ************************************************
    // *****************************************************************************************
    struct CPParams : public Params
    {
        //int Binning = 0;
    };

    void ControlPlot(const UnfoldingIOHandler* unfoldingIOHandler, const TString& fileName, CPParams pars);
    void ControlPlot(const UnfoldingIOHandler* unfoldingIOHandler, const std::vector<TH1D*> vh, const TString& fileName, CPParams pars)
    {
      const std::pair<TH1D*, TH1D*> totUDZero = std::pair<TH1D*, TH1D*>(NULL, NULL);
      ControlPlot(unfoldingIOHandler, vh, totUDZero, fileName, pars);
    };
    void ControlPlot(const UnfoldingIOHandler* unfoldingIOHandler, const std::vector<TH1D*> vh, const std::pair<TH1D*, TH1D*> vMCTotUD, const TString& fileName, CPParams pars);

    // *****************************************************************************************
    // *************************** CROSS SECTION ***********************************************
    // *****************************************************************************************
    struct XSecParams : public Params
    {
        XSecParams(const UnfoldingIOHandler* unfoldingIOHandler = NULL)
        {
          if(unfoldingIOHandler && unfoldingIOHandler->Dim() == 3)
          {
            RatMin = 0.5;
            RatMax = 1.5;
          }
          //Binning = 1;
          DivideBW = unfoldingIOHandler ? unfoldingIOHandler->DivideBW : true;
          if(unfoldingIOHandler && unfoldingIOHandler->DoRj)
          {
            //YMin = 0.3;
            //YMax = 0.8;
            YMin = 0.0;
            YMax = 1.0;
            RatMin = 0.8;
            RatMax = 1.2;
          }
        }
    };
    void XSec(const UnfoldingIOHandler* unfoldingIOHandler, const std::vector<TH1D*> hDat,
                     std::vector<TH1D*> hTh, const TString& fileName, XSecParams pars);
    void XSec(const UnfoldingIOHandler* unfoldingIOHandler, const std::vector<TH1D*> hDat, const std::vector<std::pair<TH1D*, TH1D*> > vDataTotUD,
                        std::vector<TH1D*> hTh, const TString& fileName, XSecParams pars, const std::pair<TH1D*, TH1D*> vMCTotUD = std::pair<TH1D*, TH1D*>(NULL, NULL), bool draw_channel_label=true);

    // *****************************************************************************************
    // *************************** VARIATIONS **************************************************
    // *****************************************************************************************
    static void PlotVars(const UnfoldingIOHandler* unfoldingIOHandler,
                         const std::vector<std::vector<TH1D*> >& vhVarDat,
                         const std::vector<TH1D*>& vhVarGen,
                         const std::vector<TH1D*>& vhVarRec,
                         const TString& fileName, const int flagAbs, const int flagEig = 0);
    static void PlotVars(const UnfoldingIOHandler* unfoldingIOHandler,
                         const std::vector<TH1D*>& vhVarGen,
                         const std::vector<TH1D*>& vhVarRec,
                         const TString& fileName, const int flagAbs, const int flagEig = 0);

    // *****************************************************************************************
    // *************************** UNCERTAINTY SUMMARY *****************************************
    // *****************************************************************************************
    struct UncSummaryParams : public Params
    {
        int NLegendRows = 1;
        bool FlagAbs = 0;
        double Rat = -1.0;
    };

    std::vector<TString> GetUncSummaryList(const std::vector<TString>& selectVar, TString& year);

    static void PlotUncSummary(const UnfoldingIOHandler* unfoldingIOHandler, const std::vector<std::vector<TH1D*> >& vhMarkers,
                               const std::vector<TH1D*>& vhTotal, const std::vector<TH1D*>& vhDashed,
                               const TString& fileName, const UncSummaryParams& params);

    void DrawCMSTitle(int flag, int dim, TString plot_type);

    void DrawChannelLabel(TString channel, int flag, int dim, TString plot_type, double x0Leg, double x1Leg, double y1Leg);

    static void SetYMinMax(std::vector<std::vector<TH1D*> >& hh, Params* pars);

    void AddXSec(UnfoldingIOHandler* unfoldingIOHandler);

    void AddXSec(std::vector<UnfoldingIOHandler*> unfoldingIOHandler);

    void SetVar(const TString& var);

    static std::vector<TString> CreateSuffVector();

    std::vector<ZUnfold*> CreateZUnfoldVector(const int n = 3);

    void LoopDataTreePlain(const TString& ch);

    void StandaloneLoopMCTreeGenPlain(const TString& ch, const int binning, TH1D* h, UnfoldingIOHandler* unfoldingIOHandler, const TString fileIn = "");

    void StandaloneLoopMCTreeGenPlain(const TString& ch, const int binning, std::vector<TH1D*> vH, std::vector<UnfoldingIOHandler*> vXSec, const TString fileIn = "", const TString weight = "");

    void LoopMCTreePlain(const TString& ch);

    void LoopMCBkgTreePlain(const TString& ch);

    void LoopMCBkgTTbarTreePlain(const TString& ch);

    void LoopAllTreesPlain(const TString& ch);

    std::vector<TH1D*> GetNPCorrHistos(std::vector<UnfoldingIOHandler*> vXSec, const TString& inputSuffix, const TString weight = "");

    void Analyse(const TString& ch);

    void StoreAndPlotNPCorr(const std::vector<std::vector<std::vector<TH1D*> > > vvvH, const TString& dirNPCorr);

    void StoreAndPlotNPCorrAll4(const std::vector<std::vector<std::vector<TH1D*> > > vvvH, const TString& dirNPCorr);

    void StoreAndPlotNPCorrHad2(const std::vector<std::vector<std::vector<TH1D*> > > vvvH, const TString& dirNPCorr);

    void StoreAndPlotNPCorrMPI3(const std::vector<std::vector<std::vector<TH1D*> > > vvvH, const TString& dirNPCorr);

    void StoreAndPlotNPCorrDef2(const std::vector<std::vector<std::vector<TH1D*> > > vvvH, const TString& dirNPCorr);

    void StoreAndPlotNPCorrCheck1(const std::vector<std::vector<std::vector<TH1D*> > > vvvH, const TString& dirNPCorr);

    void StoreAndPlotNPCorrVars(const std::vector<std::vector<std::vector<TH1D*> > > vvvH, const TString& dirNPCorr);

    void StoreAndPlotNPCorrVarsPS(const std::vector<std::vector<std::vector<TH1D*> > > vvvH, const TString& dirNPCorr);

    void StoreAndPlotNPCorrVarsME(const std::vector<std::vector<std::vector<TH1D*> > > vvvH, const TString& dirNPCorr);

    void CalculateNPCorr();

    void DoXSec(UnfoldingIOHandler* unfoldingIOHandler, const TString& ch, const std::vector<ZUnfold*> vAlg, const TString& dirOut, XSecParams pars,
                const bool flagStoreExtraPlots = true);

    void StoreRecUnf(const UnfoldingIOHandler* unfoldingIOHandler, const ZUnfold* unf, const bool flagNorm, const TString& dirOut);

    std::vector<std::pair<TString, std::vector<TString> > > CreateSysToPlotVector(const bool flagSingleVars, const std::vector<TString>& vSys);

    void PlotAllVars(const UnfoldingIOHandler* unfoldingIOHandler,
                     const std::vector<ZUnfold*> vUnf,
                     const std::map<TString, TH1D*>& mhVarGen,
                     const std::map<TString, TH1D*>& mhVarRec,
                     std::vector<std::map<TString, TH1D*> > vmhVarDat,
                     std::vector<std::pair<TString, std::vector<TString> > > sysToPlot,
                     const bool flagPlotGenRecOnly,
                     const TString& sufXSec,
                     const TString& dirOutBase
                     );


    double SetScales(ZPredSetup& pred, const std::vector<TString>& vStrMu, const int mu);

    void PlotScaleVars(ZMadGraph* mg, ZPredSetup pred, const UnfoldingIOHandler* unfoldingIOHandler, const bool flagNorm, const std::vector<TString>& vStrMu, const TString& fileName);

    std::vector<TH1D*> GetVHMtDep(const ZMadGraph* mg, const ZMeasur* measXSec);

    void ProcessSystematics();

    void PlotMGTheory(ZMadGraph* mg, const TString fileName, const std::vector<ZPredSetup>& vPredSetup, const std::vector<TString> vTitle,
                      const int dof, const TString genTitle, const UnfoldingIOHandler* unfoldingIOHandler, ZMeasur* measXSec, std::pair<TH1D*, TH1D*> hXSecToPlotUD, TH1D* hXSecGenNom,
                      const bool flagPlotPaper = false);

    // (with several sets of data for mt dependence)
    void PlotMGTheory(ZMadGraph* mg, const TString fileName, const std::vector<ZPredSetup>& vPredSetup, const std::vector<TString> vTitle,
                      const int dof, const TString genTitle, const UnfoldingIOHandler* unfoldingIOHandler, std::vector<TH1D*> vhData, const TH2D* hCov, std::pair<TH1D*, TH1D*> hXSecToPlotUD, TH1D* hXSecGenNom,
                      const bool flagPlotPaper = false);

    // simplified version with theory histograms provided directly
    void PlotMGTheory(const TString fileName, const std::vector<TH1D*>& vPredH,
                      const TString genTitle, const UnfoldingIOHandler* unfoldingIOHandler, ZMeasur* measXSec, std::pair<TH1D*, TH1D*> hXSecToPlotUD, TH1D* hXSecGenNom);


    ZMeasur* FillMeas(const std::map<TString, TH1D*>& mapVar, const std::vector<TString>& selectVar, const TH2D* hStatCov = NULL, const bool flagShapeOnly = false);

    std::vector<TH1D*> GetVectorUncUD(const ZMeasur& meas, const TString& type);

    std::vector<TH1D*> PlotJESSubtotalUnc(const UnfoldingIOHandler* unfoldingIOHandler, const std::map<TString, TH1D*>& vmhVarDat, const TString& nameSubtotal,
                                          const std::vector<TString>& vJESsrc, std::vector<TH1D*> vhUncUnf,
                                          const TString& fileNameBase, const UncSummaryParams& pars, const bool flagPlotSrc = true);

    void CoverTest(const TString& ch, const int type, const int niter = -1);

    void CoverTestWrite(std::vector<ZUnfold*>& vUnf, const int niter, const TString suffix, const int flagRew, const int seed, const TString& ch);

    void CoverTestRead(std::vector<ZUnfold*>& vUnf, const int niter, const TString suffix, const TString& ch);

    static int NDigFloatXSec;
    static int NDigFloatStat;
    static int NDigFloatSyst;
    static int NDigFloatCorr;
    static int NDigFloatVar;

    static void PrintXsecAllAlg(std::vector<ZUnfold*> vAlg, const TString& fileName);

    // nominal x-section and covariance matrix
    static void PrintDataXsec(const TH1D* vH, const TH1D* vT, const TString& fileName);
    static void PrintDataCorMatrix(const TH2D* hCor, const TString& fileName);
    static std::pair<std::pair<TH1D*, TH1D*>, TH2D*> ReadDataXsecAndCov(const TString& fileNameBase);
    //                          hDat  hGen               Stat            systU   systD
    static std::pair<std::pair<TH1D*, TH1D*>, std::pair<TH1D*, std::pair<TH1D*, TH1D*>>> ReadDataXsecWithSystAndStats(const TString& fileNameBase, bool get_total_uncert_UD=false);

    // detector level
    static std::pair<std::pair<TH1D*, TH1D*>, TH2D*> ReadRecUnf(const UnfoldingIOHandler* unfoldingIOHandler, const TString& dirName);

     // with total u/d systematics
    static void PrintDataXsec(const TH1D* vH, const std::pair<TH1D*, TH1D*> hud, const TH1D* vT, const TString& fileName);

    // LaTex table with binning
    static void PrintTexBinning(const UnfoldingIOHandler* unfoldingIOHandler, const TString& fileName);
    static void PrintTexTabXSec(const UnfoldingIOHandler* unfoldingIOHandler, const ZMeasur* meas, const TString& fileName);
    static void PrintTexTabCorr(const UnfoldingIOHandler* unfoldingIOHandler, const ZMeasur* meas, const TString& fileName);
    static void PrintTexTabSyst(PlotterConfigurationHelper *configHelper__, const UnfoldingIOHandler* unfoldingIOHandler, const ZMeasur* meas,
                                std::vector<TString> &systematics_list, TString fileName);

    // ...
    static void PrintDataXsecFinal(const TString& fileNameIn, const TString& fileNameOut, const TString& fileNameOutSys, const std::vector<std::vector<std::vector<TH1D*> > >& vvvSysAll, const std::vector<TString>& vShortName, const std::vector<std::vector<TH1D*> >& vvSysTot, const int a);

    // control plot
    static void PrintCP(const std::vector<TH1D*> vh, const TString& fileName);
    static void PrintCPwithMCunc(const std::vector<TH1D*> vh, const std::pair<TH1D*, TH1D*> pairUD, const TString& fileName);
    static std::vector<TH1D*> ReadCP(const TString& fileName);

    // all variations
    static void PrintEigenvectors(const TH1D* hNom, const std::pair<TH1D*, TH1D*> hUD, const std::map<std::pair<TString, TString>, TH1D*>& mapVars, const TString& fileName);

    void CalculatePSE(const UnfoldingIOHandler* unfoldingIOHandler, const TString& fileName, const bool flagOrigBinning1D = true);

    TH2D* DrawMigrationMatrix(const UnfoldingIOHandler* unfoldingIOHandler, const TString& fileName = "", const bool flagOrigBinning1D = true) const;
    void PrintMigrationMatrix(const UnfoldingIOHandler* unfoldingIOHandler, const TString& fileName = "");
    void PlotMigrationMatrixFromPlainVector(std::vector<float> plain_matrix, const UnfoldingIOHandler* unfoldingIOHandler, const TString& fileName = "", std::vector<int> bin_info = {}, const bool flagOrigBinning1D=true);
    void PlotMultipleMigrationMatrixFromPlainVectors(std::map<TString, std::vector<float>> plain_matrices, const UnfoldingIOHandler* unfoldingIOHandler, const TString& fileName, std::vector<float> extra_bin_info = {}, TString extra_bin_info_label="", TString plot_title="", const bool flagOrigBinning1D=true);
    std::vector<float> ReadMigrationMatrix(const TString& fileName, const bool convert_to_percent=false);
};

#endif
