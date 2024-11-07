#ifndef ZMADGRAPH_H
#define ZMADGRAPH_H

#include <TString.h>
#include <TFile.h>
#include <TH1D.h>
#include <vector>
#include <map>
#include <cassert>

#include "PlotterConfigurationHelper.h"

class TH2D;
class UnfoldingXSecHandler;
class UnfoldingIOHandler;
class ZMeasur;
class TFile;
class TextFile;

class ZPredSetup
{
  public:
    enum ZPredSetupPDF
    {
      CT14nlo = 10,//13100,
      CT14nlo_as_0111 = 67,//13158,
      CT14nlo_as_0113 = 68,//13160,
      CT14nlo_as_0116 = 69,//13163,
      CT14nlo_as_0118 = 70,//13165,
      CT14nlo_as_0120 = 71,//13167,
      CT14nlo_as_0123 = 72,//13170,
      MMHT2014nlo68cl = 73,//25100,
      MMHT2014nlo_asmzlargerange = 124,//25270,
      HERAPDF20_NLO_EIG = 145,//61100,
      HERAPDF20_NLO_VAR = 174,//61130,
      HERAPDF20_NLO_ALPHAS_110 = 188,//61710,
      HERAPDF20_NLO_ALPHAS_113 = 189,//61713,
      HERAPDF20_NLO_ALPHAS_116 = 190,//61716,
      HERAPDF20_NLO_ALPHAS_118 = 191,//61718,
      HERAPDF20_NLO_ALPHAS_122 = 192,//61722,
      HERAPDF20_NLO_ALPHAS_126 = 193,//61726,
      NNPDF31_nlo_as_0116 = 194,//318900,
      NNPDF31_nlo_as_0118 = 195,//303400,
      NNPDF31_nlo_as_0120 = 296,//319100,
      ABMP16als114_5_nlo = 297,//43114,
      ABMP16als118_5_nlo = 327,//43118,
      ABMP16als122_5_nlo = 357,//43122,
      JR14NLO08VF = 387,//81200,
      CJ15nlo = 426,//44001,
      //
      Undefined = -1
    };

    static TString GetLHAPDFName(const ZPredSetupPDF pdf)
    {
      static std::vector<TString> v;
      if(v.size() == 0)
      {
        v = std::vector<TString>(500);
        v[ZPredSetupPDF::CT14nlo] = "CT14nlo";
        v[ZPredSetupPDF::CT14nlo_as_0111] = "CT14nlo_as_0111";
        v[ZPredSetupPDF::CT14nlo_as_0113] = "CT14nlo_as_0113";
        v[ZPredSetupPDF::CT14nlo_as_0116] = "CT14nlo_as_0116";
        v[ZPredSetupPDF::CT14nlo_as_0118] = "CT14nlo_as_0118";
        v[ZPredSetupPDF::CT14nlo_as_0120] = "CT14nlo_as_0120";
        v[ZPredSetupPDF::CT14nlo_as_0123] = "CT14nlo_as_0123";
        v[ZPredSetupPDF::MMHT2014nlo68cl] = "MMHT2014nlo68cl";
        v[ZPredSetupPDF::MMHT2014nlo_asmzlargerange] = "MMHT2014nlo_asmzlargerange";
        v[ZPredSetupPDF::HERAPDF20_NLO_EIG] = "HERAPDF20_NLO_EIG";
        v[ZPredSetupPDF::HERAPDF20_NLO_VAR] = "HERAPDF20_NLO_VAR";
        v[ZPredSetupPDF::HERAPDF20_NLO_ALPHAS_110] = "HERAPDF20_NLO_ALPHAS_110";
        v[ZPredSetupPDF::HERAPDF20_NLO_ALPHAS_113] = "HERAPDF20_NLO_ALPHAS_113";
        v[ZPredSetupPDF::HERAPDF20_NLO_ALPHAS_116] = "HERAPDF20_NLO_ALPHAS_116";
        v[ZPredSetupPDF::HERAPDF20_NLO_ALPHAS_118] = "HERAPDF20_NLO_ALPHAS_118";
        v[ZPredSetupPDF::HERAPDF20_NLO_ALPHAS_122] = "HERAPDF20_NLO_ALPHAS_122";
        v[ZPredSetupPDF::HERAPDF20_NLO_ALPHAS_126] = "HERAPDF20_NLO_ALPHAS_126";
        v[ZPredSetupPDF::NNPDF31_nlo_as_0116] = "NNPDF31_nlo_as_0116";
        v[ZPredSetupPDF::NNPDF31_nlo_as_0118] = "NNPDF31_nlo_as_0118";
        v[ZPredSetupPDF::NNPDF31_nlo_as_0120] = "NNPDF31_nlo_as_0120";
        v[ZPredSetupPDF::ABMP16als114_5_nlo] = "ABMP16als114_5_nlo";
        v[ZPredSetupPDF::ABMP16als118_5_nlo] = "ABMP16als118_5_nlo";
        v[ZPredSetupPDF::ABMP16als122_5_nlo] = "ABMP16als122_5_nlo";
        v[ZPredSetupPDF::JR14NLO08VF] = "JR14NLO08VF";
        v[ZPredSetupPDF::CJ15nlo] = "CJ15nlo";
      }
      assert(pdf >= 0 && pdf < v.size());
      return v[pdf];
    }

    static int GetNScales(const double mr, const double mf);

    ZPredSetupPDF PDF = Undefined;
    int PDFmem = -1;
    double Mur = -1.0;
    double Muf = -1.0;
    double Mt = -1.0;

    int NMu = -1;

    // optional
    TString Title;
    int MarkerStyle = -1;
    int MarkerColor = -1;

    ZPredSetup() {}

    ZPredSetup(const ZPredSetupPDF& pdf, const int mem, const double mur, const double muf, const double mt,
               const TString title = "", const int style = 1, const int color = 1)
    {
      PDF = pdf;
      PDFmem = mem;
      Mur = mur;
      Muf = muf;
      Mt = mt;
      Title = title;
      MarkerStyle = style;
      MarkerColor = color;
      SetNMu();
    }

    static double GetNMu(const double mur, const double muf)
    {
      if(mur == 1.0 && muf == 1.0)
        return 0;
      else if(mur == 0.5 && muf == 1.0)
        return 1;
      else if(mur == 2.0 && muf == 1.0)
        return 2;
      else if(mur == 1.0 && muf == 0.5)
        return 3;
      else if(mur == 1.0 && muf == 2.0)
        return 4;
      else if(mur == 0.5 && muf == 0.5)
        return 5;
      else if(mur == 2.0 && muf == 2.0)
        return 6;
      // 16.07.2018
      else if(mur == 0.0 && muf == 12.0)
        return 7;
      else if(mur == 0.0 && muf == 13.0)
        return 8;
      else
        throw std::logic_error("Error: should not be here");
    }

    void SetNMu()
    {
      NMu = GetNMu(Mur, Muf);
    }
};

class ZPredHisto
{
  public:
    ZPredHisto(const ZPredSetup& predSetup, const int histoID, const int order)
    {
      PredSetup = predSetup;
      HistoID = histoID;
      Order = order;
    }

    int Order = -1;
    int HistoID = -1;
    ZPredSetup PredSetup;

    bool operator <(const ZPredHisto& rhs) const
    {
      bool result = false;

      result = this->HistoID < rhs.HistoID;
      if(result)
        return result;

      result = this->Order < rhs.Order;
      if(result)
        return result;

      result = this->PredSetup.Mt < rhs.PredSetup.Mt;
      if(result)
        return result;

      result = this->PredSetup.Muf < rhs.PredSetup.Muf;
      if(result)
        return result;

      result = this->PredSetup.Mur < rhs.PredSetup.Mur;
      if(result)
        return result;

      result = this->PredSetup.PDFmem < rhs.PredSetup.PDFmem;
      if(result)
        return result;

      result = this->PredSetup.PDF < rhs.PredSetup.PDF;
      if(result)
        return result;
    }
};

#define ZMadGraph_MAXHISTONBINS 10
#define ZMadGraph_MAXNMT 5
#define ZMadGraph_MAXNORDER 3
#define ZMadGraph_MAXNSCALES 9

class ZMadGraph
{
  private:
    double timeReadFiles = 0.0;
    double timeReadPred = 0.0;
    double timeProcess = 0.0;
    double timeAddVar = 0.0;

    PlotterConfigurationHelper *configHelper_;

    std::vector<std::vector<TFile*> > vvFile;
    std::vector<std::vector<std::vector<TH1D*> > > vvvHistos;
    std::vector<std::vector<TString> > vvDirAG;
    std::map<TString, TextFile*> mFileAG;
    //std::vector<std::vector<std::vector<std::vector<std::vector<TH1D*> > > > > vvvvvHistosAG; // [order] [mt] [scales] [PDF] [histoID]
    TH1D* vvvvvHistosAG[ZMadGraph_MAXNORDER][ZMadGraph_MAXNMT][ZMadGraph_MAXNSCALES][500][500]; // [order] [mt] [scales] [PDF] [histoID]
    float v5PredAG[ZMadGraph_MAXNORDER][ZMadGraph_MAXNMT][ZMadGraph_MAXNSCALES][500][500][ZMadGraph_MAXHISTONBINS]; // [order] [mt] [scales] [PDF] [histoID]
    char v5PredAGnbins[ZMadGraph_MAXNORDER][ZMadGraph_MAXNMT][ZMadGraph_MAXNSCALES][500][500]; // [order] [mt] [scales] [PDF] [histoID]

    //std::vector<std::vector<std::vector<std::map<const UnfoldingIOHandler*, TH1D*   > > > > vvvmPredNom;
    //std::vector<std::vector<std::vector<std::map<const UnfoldingIOHandler*, TH2D*   > > > > vvvmPredCov;
    std::vector<std::vector<std::vector<std::map<const UnfoldingIOHandler*, ZMeasur*> > > > vvvmPredMeas;

    void ResizeStoredPredictions();

    ZMeasur*& PredMeas(const UnfoldingIOHandler* unfoldingIOHandler, const int npdf, const int nmu, const int nmt)
    {
      return vvvmPredMeas[npdf][nmu][nmt][unfoldingIOHandler];
    }

  public:

    ZMadGraph(PlotterConfigurationHelper *configHelper, const bool flagAG = true, const bool flagUseNomPDFCovMat = true);
    ~ZMadGraph();

    bool FlagAG = 1;
    bool FlagUseNomPDFCovMat = 0;

    const int NWeightOffset = 10000;
    const int NMaxHistoID = 500;
    const int NOrder = 3;

    const double MtNom = 172.5;

    //std::vector<double> VMt = { 172.5, 167.5, 177.5 };
    //std::vector<TString> VMtName = { "mt1725", "mt1675", "mt1775" };
    const std::vector<double> VMt = { 172.5, 167.5, 177.5, 170.0, 175.0 };
    const std::vector<TString> VMtName = { "mt1725", "mt1675", "mt1775", "mt1700", "mt1750" };
    const int NMtNom = 0;

    //std::vector<int> VMu = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    const int NPDF = 500;

    int GetNMt(const double mt) const;

    TString DirNPCorr;

    static int GetOffsetHistoPtjThreshold(PlotterConfigurationHelper *configHelper__, const double pt);
    static int GetOffsetHistoPtttThreshold(PlotterConfigurationHelper *configHelper__, const double pt);

    const TH1D* ReadPred(const ZPredHisto& predHisto);
    float* ReadPredFast(const ZPredHisto& predHisto, int* nbinsPtr = NULL);
    std::vector<TH1D*> GetPredictionWithUnc(const UnfoldingIOHandler* unfoldingIOHandler, const ZPredSetup& predNom, const std::vector<ZPredSetup>& vPredUnc, const TString& uncType, ZMeasur* th, const bool flagNorm, const bool flagDivideBW, const double rescale = 1.0);

    /*TH1D*& PredNom(const UnfoldingIOHandler* unfoldingIOHandler, const int npdf, const int nmu, const int nmt)
    {
      return vvvmPredNom[npdf][nmu][nmt][unfoldingIOHandler];
    }

    TH2D*& PredCov(const UnfoldingIOHandler* unfoldingIOHandler, const int npdf, const int nmu, const int nmt)
    {
      return vvvmPredCov[npdf][nmu][nmt][unfoldingIOHandler];
    }*/

    ZMeasur*& PredMeas(const UnfoldingIOHandler* unfoldingIOHandler, const std::pair<int, int> pdf, const std::pair<double, double> mu, const double mt)
    {
      const int npdf = pdf.first + pdf.second;
      const int nmu = ZPredSetup::GetNMu(mu.first, mu.second);
      const int nmt = GetNMt(mt);
      return vvvmPredMeas[npdf][nmu][nmt][unfoldingIOHandler];
    }

    std::vector<ZMeasur*> PredMeasVec(const UnfoldingIOHandler* unfoldingIOHandler, const std::pair<int, int> pdf, const std::pair<double, double> mu, const std::vector<double>& vMt);
    std::vector<TH2D*> PredCovVec(const UnfoldingIOHandler* unfoldingIOHandler, const std::pair<int, int> pdf, const std::pair<double, double> mu, const std::vector<double>& vMt);

    int GetFullHistoID(const ZPredHisto& predHisto) const; // as used in MG

    std::vector<ZPredSetup> GetPredGroup(const TString& name);
    ZPredSetup GetPred(const TString& name);

    ZPredSetup::ZPredSetupPDF GetPDFForCovMatrix(ZPredSetup::ZPredSetupPDF pdfIn);
};

#endif // ZMADGRAPH_H
