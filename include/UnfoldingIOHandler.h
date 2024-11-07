#ifndef UnfoldingIOHandler_h
#define UnfoldingIOHandler_h

#include <vector>
#include <TUnfoldBinning.h>
#include "ttmd_vars.h"
#include <PlotterConfigurationHelper.h>

class TF1;
class TString;
class TH1;
class TH1D;
class TH2;
class TH2D;
class TreePlain;
class ZMadGraph;
class ZPredSetup;
//class ZMeasur;

/*#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "vars.h"
#include <TUnfoldBinning.h>
//#include "mystaff.cxx"
#include <TPaveText.h>
#include <TGraphErrors.h>
#include "framework.h"
#include "settings.h"
#include "TRandom3.h"
#include "TF1.h"*/

class UnfoldingIOHandler
{
  private:
    std::vector<ZVar*> zVVar;
    std::vector<TH1D*> zHDat;
    std::vector<TH1D*> zHDnb;
    std::vector<TH1D*> zHRec;
    std::vector<TH1D*> zHGen;
    std::vector<TH1D*> zHBgr;
    std::vector<TH1D*> zHBtt;
    std::vector<TH1D*> zHOutOfSpace;
    std::vector<TH2D*> zHRsp;
    std::vector<double> zRecWeight;
    TF1* zFRew;
    std::vector<std::pair<ZVar*, TF1*> >* zvVFRew;
    //TString zDir;

    PlotterConfigurationHelper *configHelper_;

    bool flagApplyCut = false;

    void FillTH1(TH1* h, const TreePlain* tree, const int binning, double w);

    // use weight from treeX, reweight treeY
    void FillTH2(TH2* h, const TreePlain* treeX, const int binningX, const TreePlain* treeY, const int binningY, const double w);

    void FillTH2UnderflowX(TH2* h, const TreePlain* treeY, const int binningY, const double w);

  public:

    TString Expression;
    TString Suffix;
    UnfoldingIOHandler(const std::vector<ZVar*>& argVVar, PlotterConfigurationHelper *configHelper, const TString& suffix = "");
    UnfoldingIOHandler(const UnfoldingIOHandler& old, PlotterConfigurationHelper *configHelper);

    int DoXSec = 0; // 0 nothing, 1 cross-section level only, 2 detector level only, 3 both cross-section and detector level
    //bool DoXSecFineDetOnly = false; // only detector level distributions: do not unfold
    bool DoYields = false;
    bool DoPSE = false;
    bool DoParticle = false;
    //bool DoEfficiency = false;
    //bool DoPurStab = false;
    bool DoRj = false;
    bool DivideBW = true;
    std::vector<int> VSkip;

    int NPertOrders = 1; // need multiple as order predictions
    bool FlagNeedNP = false;

    float RatMin = +1001.0;
    float RatMax = -1001.0;

    // ZVar* Var(const int v) const { return zVVar[v]; }
    ZVar* Var(const int v) const { return zVVar[v]; }
    std::vector<ZVar*> Vars() const { return zVVar; }

    void SetReweightFunction(TF1* argF) { zFRew = argF; }
    const TF1* GetReweightFunction() const { return zFRew; }
    void AddVarReweight(ZVar* var, TF1* f)
    {
      if(!zvVFRew)
        zvVFRew = new std::vector<std::pair<ZVar*, TF1*> >;
      zvVFRew->push_back(std::pair<ZVar*, TF1*>(new ZVar(*var), f));
    }

    //TString Dir() const { return (zDir != "") ? zDir : "plots/"; }
    //void SetDir(const TString& str) { zDir = str; }

    int Dim() const { return zVVar.size(); }
    int NBinsC() const { return zVVar[0]->NBins(0); }
    int NBinsF() const { return zVVar[0]->NBins(1); }
    const std::vector<double>& BinsC() const { return zVVar[0]->Bins[0]; }
    const std::vector<double>& BinsF() const { return zVVar[0]->Bins[1]; }
    int NBinnings() const { return zVVar[0]->NBinnings(); }

    std::vector<TUnfoldBinning*> zTUnfBinning;
    std::vector<TUnfoldBinning*> zTUnfBinningNode;

    const TH1D* HDat(const int binning) const { return zHDat[binning]; }
    const TH1D* HDatF() const { return zHDat[1]; }
    const TH1D* HDatC() const { return zHDat[0]; }
    const TH1D* HDnb(const int binning) const { return zHDnb[binning]; }
    const TH1D* HDnbF() const { return zHDnb[1]; }
    const TH1D* HDnbC() const { return zHDnb[0]; }
    const std::vector<TH1D*>& HRec() const { return zHRec; }
    const TH1D* HRec(const int binning) const { return zHRec[binning]; }
    const TH1D* HRecF() const { return zHRec[1]; }
    const TH1D* HRecC() const { return zHRec[0]; }
    const TH1D* HGen(const int binning) const { return zHGen[binning]; }
    const TH1D* HGenF() const { return zHGen[1]; }
    const TH1D* HGenC() const { return zHGen[0]; }
    const TH1D* HBgr(const int binning) const { return zHBgr[binning]; }
    const TH1D* HBgrF() const { return zHBgr[1]; }
    const TH1D* HBgrC() const { return zHBgr[0]; }
    const TH1D* HBtt(const int binning) const { return zHBtt[binning]; }
    const TH1D* HBttF() const { return zHBtt[1]; }
    const TH1D* HBttC() const { return zHBtt[0]; }
    const TH1D* HOutOfSpace(const int binning) const { return zHOutOfSpace[binning]; }
    const TH1D* HOutOfSpaceF() const { return zHOutOfSpace[1]; }
    const TH1D* HOutOfSpaceC() const { return zHOutOfSpace[0]; }
    const TH2D* HRspC() const { return zHRsp[0]; }
    const TH2D* HRspTUnfold() const { return zHRsp[1]; }

    double GetRecWeight(const int sample) const { return zRecWeight[sample]; }
    void PrintRecWeights(const TString &ch, TString Var, int kin_reco_mode) const;

    void SetPseudoData(const std::vector<TH1D*>& v);
    void GeneratePseudoData(const int binning, const TH1D* hGen, const TH2D* hRsp, const TH1D* hRefUnc = NULL, const bool flagFluct = true);

    TUnfoldBinning* GetTUnfoldBinningF() const { return zTUnfBinning[1]; }
    TUnfoldBinning* GetTUnfoldBinningC() const { return zTUnfBinning[0]; }

    TH2D* HRspTUnfold11() const;

    void Init();

    void Reset(const TString& opt = "drgob");

    int NBinsXSec(const int binning) const;

    double Reweight(const TreePlain* tree);

    double GetBin(const TreePlain* tree, const int binning) const;

    int GetBin(const double* val, const int binning) const;

    void FillData(const TreePlain* tree, const double rew = 1.0);

    void FillBgr(const TreePlain* tree, const double rew = 1.0);

    void FillBtt(const TreePlain* tree, const double rew = 1.0);

    void FillOutOfSpace(const TreePlain* tree, const double rew = 1.0);

    void FillMC(const TreePlain* treeRec, const TreePlain* treeGen, const double rew = 1.0);

    void StandaloneFillMCGen(const TreePlain* treeGen, const int binning, TH1D* h, const double rew = 1.0);

    void AddOutOfSpaceToTTother();

    void SubtrackBackground(const bool flagDo = true);

    void FillLastDimHisto(TH1D* hFill, const TH1D* hGlobal, const bool flagNorm = 0) const;

    std::vector<TH1D*> CreateLastDimHistos(const int binning, const int dim = 0, std::vector<double> varArray = {}) const;

    TH1D* GetPlainHisto(const int binning = 0) const { return (TH1D*) zTUnfBinning[binning]->CreateHistogram(""); }

    void InitTree(TreePlain* tree, bool flagNoRew = false);
    //void InitTreeForReweighting(TreePlain* tree);
    void ClearTree(TreePlain* tree, bool flagNoRew = false);

    std::vector<std::pair<std::vector<int>, int> > VMGHisto;
    std::vector<std::pair<std::vector<int>, int> > VMGHistoMinus;
    bool FlagMGHistoIntegrate = false;
    std::vector<int> vMGHistoOneBinOnly;

    // TODO: move to MadGraph class?
    void PutPredictionPttt(double ptThreshold);
    void PutPredictionPtttIntegrate(double ptThreshold);

    static std::pair<std::vector<double>,std::vector<double>> GetBinningInfoFromFile(const TString binningFileName);
};

#endif
