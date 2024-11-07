#ifndef TTBAR_treePseudo_h
#define TTBAR_treePseudo_h

#include <TChain.h>
#include <TFile.h>
class ZUnfold;
class TH1D;
class TH2D;

class ZTreePseudo: public TChain
{
  private:
    TFile* _f;
    int _nalg;
    int _nbins;
    int _nbinsUndiag;
    int _niter;
    std::vector<int> _chi2Dim;
    std::vector<int> _isCov;
    std::vector<std::vector<float> > _brChi2;
    std::vector<float> _brTau;
    std::vector<std::vector<float> > _brHUnf;
    std::vector<std::vector<float> > _brHCovDiag;
    std::vector<std::vector<float> > _brHCovUndiag;
    std::vector<TBranch*> _branch_tau;
    
  public:
    int MaxEvents;
  
    ZTreePseudo(const std::vector<ZUnfold*>& vUnf, const int nbins, const TString fileName);
    ~ZTreePseudo();

    float Tau(const int a) { return _brTau[a]; }
    float Chi2A(const int a) { return _brChi2[a][0]; }
    float Chi2L(const int a) { return _brChi2[a][1]; }    
    float Unf(const int a, const int b) { return _brHUnf[a][b]; }
    float CovDiag(const int a, const int b) { return _brHCovDiag[a][b]; }
    float CovUndiag(const int a, const int b) { return _brHCovUndiag[a][b]; }
    int Chi2Dim(const int a) { return _chi2Dim[a]; }

    void HUnfCov(const int a, TH1D* hunf, TH2D* hcov);
    
    std::pair<TH1D*, TH2D*> HUnfCov(const int a);
    
    TH1D* HGen() { return (TH1D*) _f->Get("hgen"); }
    TH1D* HGenRsp() { return (TH1D*) _f->Get("hgenRsp"); }
    
    void GetPulResMinMax(std::vector<double>& vMinPul, std::vector<double>& vMaxPul, std::vector<double>& vMinRes, std::vector<double>& vMaxRes);

    std::vector<double> GetChi2ALMinMax();

    std::vector<double> GetTauMinMax();
  
    std::pair<std::vector<std::vector<TH1D*> >, std::vector<std::vector<TH2D*> > > ReadAllHUnfCov();
};

#endif
