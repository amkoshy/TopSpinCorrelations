#include "ttmd_treePseudo.h"
#include "ttmd_unfold.h"
#include "utils.h"
#include <TChain.h>

ZTreePseudo::ZTreePseudo(const std::vector<ZUnfold*>& vUnf, const int nbins, const TString fileName): TChain("treeUnf")
{
  // input file
  this->Add(fileName);
  _niter = this->GetEntries();
  _f = TFile::Open(fileName);
  
  // number of algorithms and xsec bins
  _nalg = vUnf.size();
  _nbins = nbins;
  _nbinsUndiag = nbins * (nbins - 1) / 2;
  
  // resize arrays for branches, initialise everything with nan
  _brChi2.resize(_nalg);
  _brTau.resize(_nalg, std::numeric_limits<double>::quiet_NaN());
  _brHUnf.resize(_nalg);
  _brHCovDiag.resize(_nalg);
  _brHCovUndiag.resize(_nalg);
  for(int a = 0; a < _nalg; a++)
  {
    _brChi2[a].resize(2, std::numeric_limits<double>::quiet_NaN());
    _brHUnf[a].resize(nbins, std::numeric_limits<double>::quiet_NaN());
    _brHCovDiag[a].resize(nbins, std::numeric_limits<double>::quiet_NaN());
    _brHCovUndiag[a].resize(_nbinsUndiag, std::numeric_limits<double>::quiet_NaN());
  }
  _chi2Dim.resize(_nalg, 0);
  _isCov.resize(_nalg, 0);
  _branch_tau.resize(_nalg, NULL);
  
  // setup branches
  for(int a = 0; a < _nalg; a++)
  {
    _chi2Dim[a] = vUnf[a]->GetDataChi2Dim();
    if(_chi2Dim[a] > 0)
      this->SetBranchAddress(TString::Format("chi2%s", vUnf[a]->Remark.Data()), &_brChi2[a][0]);
    if(_chi2Dim[a] == 2)
      this->SetBranchAddress(TString::Format("tau%s", vUnf[a]->Remark.Data()), &_brTau[a], &_branch_tau[a]);
    this->SetBranchAddress(TString::Format("hunf%s", vUnf[a]->Remark.Data()), &_brHUnf[a][0]);
    this->SetBranchAddress(TString::Format("hcovdiag%s", vUnf[a]->Remark.Data()), &_brHCovDiag[a][0]);
    if(vUnf[a]->Remark != "BBB")
    {
      _isCov[a] = 1;
      this->SetBranchAddress(TString::Format("hcovundiag%s", vUnf[a]->Remark.Data()), &_brHCovUndiag[a][0]);
    }
  }
}
    
ZTreePseudo::~ZTreePseudo()
{
  _f->Close();
}

void ZTreePseudo::HUnfCov(const int a, TH1D* hunf, TH2D* hcov)
{
  if(_isCov[a])
    utils::ttmd::FArraysToHCov(_nbins, &_brHUnf[a][0], &_brHCovDiag[a][0], &_brHCovUndiag[a][0], hunf, hcov);
  else
    utils::ttmd::FArraysToHCov(_nbins, &_brHUnf[a][0], &_brHCovDiag[a][0], NULL, hunf, hcov);
}
    
std::pair<TH1D*, TH2D*> ZTreePseudo::HUnfCov(const int a)
{
  TH1D* hunf = new TH1D("", "", _nbins, 0.5, _nbins + 0.5);
  TH2D* hcov = new TH2D("", "", _nbins, 0.5, _nbins + 0.5, _nbins, 0.5, _nbins + 0.5);
  HUnfCov(a, hunf, hcov);
  return std::pair<TH1D*, TH2D*>(hunf, hcov);
}
    
void ZTreePseudo::GetPulResMinMax(std::vector<double>& vMinPul, std::vector<double>& vMaxPul, std::vector<double>& vMinRes, std::vector<double>& vMaxRes)
{
  //printf("%d %d\n", _nalg, _niter);
  for(int t = 1; t < _niter; t++)
  {
    this->GetEntry(t);
    for(int u = 0; u < _nalg; u++)
    {
      /*// chi2
      if(tree.Chi2Dim(u) >= 1)
      {
    chi2Amin = std::min(chi2Amin, tree.Chi2A(u));
    chi2Amax = std::max(chi2Amax, tree.Chi2A(u));
      }
      if(tree.Chi2Dim(u) == 2)
      {
    chi2Lmin = std::min(chi2Amin, tree.Chi2L(u));
    chi2Lmax = std::max(chi2Amax, tree.Chi2L(u));
      }*/

      // pul and res
      for(int b = 0; b < _nbins; b++)
      {
        double res = this->Unf(u, b) - this->HGen()->GetBinContent(b + 1);
        double resRel = res / this->HGen()->GetBinContent(b + 1);
        vMinRes[b] = std::min(vMinRes[b], resRel);
        vMaxRes[b] = std::max(vMaxRes[b], resRel);
        double pul = res / TMath::Sqrt(this->CovDiag(u, b));
        vMinPul[b] = std::min(vMinPul[b], pul);
        vMaxPul[b] = std::max(vMaxPul[b], pul);
        //printf("%f %f\n", res, pul);
      }
    }
  }
  
  // res, pul total min max
  vMaxRes[_nbins] = *std::max_element(vMaxRes.begin(), vMaxRes.end());
  vMinRes[_nbins] = *std::min_element(vMinRes.begin(), vMinRes.end());
  vMaxPul[_nbins] = *std::max_element(vMaxPul.begin(), vMaxPul.end());
  vMinPul[_nbins] = *std::min_element(vMinPul.begin(), vMinPul.end());
  
  for(int b = 0; b <= _nbins; b++)
  {
    printf("bin %2d  res %e %e pul %e %e\n", b, vMinRes[b], vMaxRes[b], vMinPul[b], vMaxPul[b]);
  }
  //throw;
}

std::vector<double> ZTreePseudo::GetChi2ALMinMax()
{
  float chi2Amin = 1.0e10;
  float chi2Amax = 0.0;
  float chi2Lmin = 1.0e10;
  float chi2Lmax = 0.0;
  for(int t = 1; t < _niter; t++)
  {
    this->GetEntry(t);
    for(int u = 0; u < _nalg; u++)
    {
      if(this->Chi2Dim(u) >= 1)
      {
        chi2Amin = std::min(chi2Amin, this->Chi2A(u));
        chi2Amax = std::max(chi2Amax, this->Chi2A(u));
      }
      if(this->Chi2Dim(u) == 2)
      {
        chi2Lmin = std::min(chi2Amin, this->Chi2L(u));
        chi2Lmax = std::max(chi2Amax, this->Chi2L(u));
      }
    }
  }
  printf("ZTreePseudo::GetChi2ALMinMax() chi2A [%e %e] chi2L [%e %e]\n", chi2Amin, chi2Amax, chi2Lmin, chi2Lmax);
  std::vector<double> vec = {chi2Amin, chi2Amax, chi2Lmin, chi2Lmax};
  return vec;
}

std::vector<double> ZTreePseudo::GetTauMinMax()
{
  float taumin = 1.0e10;
  float taumax = 0.0;
  for(int u = 0; u < _nalg; u++)
    if(_branch_tau[u])
    {
      printf("branch name %s\n", _branch_tau[u]->GetName());
      taumin = std::min(taumin, (float) GetMinimum(_branch_tau[u]->GetName()));
      taumax = std::max(taumax, (float) GetMaximum(_branch_tau[u]->GetName()));
    }

  printf("ZTreePseudo::GetTauMinMax() [%e %e]\n", taumin, taumax);
  std::vector<double> vec = { taumin, taumax };
  return vec;
}
  
  
std::pair<std::vector<std::vector<TH1D*> >, std::vector<std::vector<TH2D*> > > ZTreePseudo::ReadAllHUnfCov()
{
  std::vector<std::vector<TH1D*> > vvHUnf = utils::ttmd::CreateVectorPtr2D<TH1D>(_nalg, _niter);
  std::vector<std::vector<TH2D*> > vvHCov = utils::ttmd::CreateVectorPtr2D<TH2D>(_nalg, _niter);
  for(int t = 0; t < _niter; t++)
  {
    this->GetEntry(t);
    for(int u = 0; u < _nalg; u++)
    {
      //printf(" t = %d u = %d\n", t, u);
      std::pair<TH1D*, TH2D*> hunfcov = this->HUnfCov(u);
      vvHUnf[u][t] = hunfcov.first;
      vvHCov[u][t] = hunfcov.second;
    }
  }
  return std::pair<std::vector<std::vector<TH1D*> >, std::vector<std::vector<TH2D*> > >(vvHUnf, vvHCov);
}
