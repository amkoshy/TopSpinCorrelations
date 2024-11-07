#ifndef UnfoldingClosureTestsHandler_h
#define UnfoldingClosureTestsHandler_h

#include <vector>
#include <TString.h>
#include <PlotterConfigurationHelper.h>
class UnfoldingXSecHandler;
class UnfoldingIOHandler;
class TH1D;
class TH2D;
class ZUnfold;
class TVirtualPad;

class UnfoldingClosureTestsHandler
{
private:
  const UnfoldingIOHandler* _unfoldingIOHandler;
  const int _nalg;
  const int _niter;
  PlotterConfigurationHelper *configHelper_;

  std::vector<TH1D*> SetupHistogramsChi2(const int nalg/*, const int niter*/, const int nbins);
  std::vector<std::vector<TH1D*> > SetupHistogramsResiduals(const int nalg/*, const int niter*/, const int nbins);
  std::vector<std::vector<TH1D*> > SetupHistogramsPulls(const int nalg/*, const int niter*/, const int nbins);

  std::pair<TH1D*, TH2D*> GetAverageAndCov(const std::vector<TH1D*>& vH);

  void PlotVHInBins(const TString& prop, const std::vector<std::vector<TH1D*> >& vvH, const TString& fileName);
  TVirtualPad* PlotVH(const TString& prop, const std::vector<TH1D*>& vH, const TString& fileName = "", TVirtualPad* cParent = NULL, TH2D* hrIn = NULL);

  void PlotPulResSummary(const TString& prop, const std::vector<std::vector<TH1D*> >& vvH, std::vector<ZUnfold*> vUnf, const TString& fileName);

  void Plot2Cov(TH2D* covEstim, TH2D* covSample, const TString& fileName);

public:
  UnfoldingClosureTestsHandler(const UnfoldingIOHandler* unfoldingIOHandler, const int nalg, const int niter, PlotterConfigurationHelper *configHelper):
   _unfoldingIOHandler(unfoldingIOHandler),
   _nalg(nalg),
   _niter(niter),
   configHelper_(configHelper)
   {
   }

  void ReadSample(const TString& rootFileName, std::vector<ZUnfold*> vUnf, const TString& dirOut);
  void WriteSample(const TString& rootFileName, std::vector<ZUnfold*> vUnf, UnfoldingIOHandler* unfoldingIOHandler, const TH1D* hGen, const TString& ch, const int seed = -1);

  TString GetRootFileNameToStore(UnfoldingIOHandler* unfoldingIOHandler, const TString& suffix, const int niter, const TString& ch, const int modeKinRecoTree);
};

#endif
