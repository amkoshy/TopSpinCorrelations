#include "ttmd_measur.h"
#include "utils.h"
#include "PlotterConfigurationHelper.h"
//#include "ttmd_settings.h"
#include <TMath.h>
#include <cassert>

int ZMeasur::nCallAddVar = 0;

ZMeasur::ZMeasur(PlotterConfigurationHelper *configHelper):
configHelper_(configHelper)
{
  zHNom = NULL;
  zHTotalUD.first = zHTotalUD.second = NULL;
  zHSystUD.first = zHSystUD.second = NULL;
  zHCovSys = NULL;
  zHCovStat = NULL;
  zHCovFull = NULL;
  zHCovSysEnvelope = NULL;
  zHCovFullEnvelope = NULL;
}

ZMeasur::~ZMeasur()
{
  if(zHNom)
    delete zHNom;
  if(zHCovSys)
    delete zHCovSys;
  if(zHCovSysEnvelope)
    delete zHCovSysEnvelope;
  if(zHCovFull)
    delete zHCovFull;
  if(zHCovFullEnvelope)
    delete zHCovFullEnvelope;
  for(auto var : zHEigens)
    delete var.second;
  for(auto sys : zHSys)
  {
    delete sys.second.first;
    delete sys.second.second;
  }
}

void ZMeasur::SetNom(const TH1D* h) { zHNom = new TH1D(*h);}

void ZMeasur::SetCovStatMatrix(const TH2D* h) { zHCovStat = new TH2D(*h); }

void ZMeasur::AddVarDiff(const TString& var, TH1D* h, int nVarInSys, TString sys)
{
  if(sys == "")
  {
    sys = configHelper_->VarToSys(var);
    assert(nVarInSys == 0 || configHelper_->GetSysVars(sys).size() == nVarInSys);
  }
  else
    assert(nVarInSys);

  // add to systematics (envelope)
  // TODO check this, adding now without thinking
  if(DoEnvelope)
    AddVarToSys(sys, h);

  // if stat. cov. matrix is set already, enable "full" quadrature calculation
  if(DoEigenvectors)
  {
    if(nVarInSys == 0)
      nVarInSys = configHelper_->GetSysVars(sys).size();
    //else
      //printf("Warning in AddVarDiff(): explicitly provided nVarInSys = %d but DoQuadrature = false\n", nVarInSys);
    if(zHEigens.find(var) != zHEigens.end())
      throw std::logic_error(TString::Format("Error in AddVar(): var = %s already set\n", var.Data()).Data());
    // 20.04.18 do not clone
    //auto added = zHEigens.insert(std::pair<TString, TH1D*>(var, (TH1D*) h->Clone()));
    //added.first->second->Scale(1.0 / TMath::Sqrt(nVarInSys));
    h->Scale(1.0 / TMath::Sqrt(nVarInSys));
    zHEigens.insert(std::pair<TString, TH1D*>(var, h));
  }
}

void ZMeasur::AddVar(const TString& var, const TH1D* h, int nVarInSys, TString sys)
{
  //nCallAddVar++;
  //if(nCallAddVar % 100 == 0)
  //  printf("nCallAddVar: %d\n", nCallAddVar);
  assert(zHNom);
  TH1D* hDiff = (TH1D*) h->Clone();
  if(FlagShapeOnly)
    hDiff->Scale(zHNom->Integral() / hDiff->Integral());
  hDiff->Add(zHNom, -1.0);
  this->AddVarDiff(var, hDiff, nVarInSys, sys);
}

void ZMeasur::AddRep(const int nrep, const TH1D* h)
{
  //if(zHRep.size() <= nrep)
  //  zHRep.resize(nrep);
  assert(!zHRep[nrep]);
  zHRep[nrep] = (TH1D*) h->Clone();
}

void ZMeasur::Process()
{
  // MC replicas
  if(zHRep.size())
  {
    assert(zHRep.size() > 1);
    int nbins = zHNom->GetNbinsX();
    size_t nrep = zHRep.size();

    // mean and RMS
    // TODO check RMS
    std::vector<double> mean(nbins);
    for(int b = 0; b < nbins; b++)
    {
      for(size_t r = 0; r < nrep; r++)
        mean[b] += zHRep[r]->GetBinContent(b + 1);
      mean[b] /= nrep;
    }
    std::vector<double> rms(nbins);
    for(int b = 0; b < nbins; b++)
    {
      for(size_t r = 0; r < nrep; r++)
        rms[b] += TMath::Power(zHRep[r]->GetBinContent(b + 1) - mean[b], 2.0);
      rms[b] = TMath::Sqrt(rms[b] / nrep);
      //rms[b] = TMath::Sqrt(rms[b] * nrep / (nrep - 1));
    }

    // uncertainties
    zHTotalUD.first = new TH1D(*zHNom);
    zHTotalUD.first->Reset();
    zHTotalUD.second = new TH1D(*zHNom);
    zHTotalUD.second->Reset();
    for(int b = 0; b < nbins; b++)
    {
      zHTotalUD.first->SetBinContent(b + 1, rms[b]);
      zHTotalUD.second->SetBinContent(b + 1, -1.0 * rms[b]);
    }
    zHSystUD.first = new TH1D(*zHTotalUD.first);
    zHSystUD.second = new TH1D(*zHTotalUD.second);

    // covariance matrix
    // TODO check
    assert(zHCovStat);
    zHCovFull = new TH2D(*zHCovStat);
    for(int b1 = 0; b1 < nbins; b1++)
    {
      for(int b2 = 0; b2 < nbins; b2++)
      {
        double val = 0.0;
        for(size_t r = 0; r < zHRep.size(); r++)
          val += (zHRep[r]->GetBinContent(b1 + 1) - mean[b1]) * (zHRep[r]->GetBinContent(b2 + 1) - mean[b2]);
        val /= nrep;
        //val *=  (1.0 * nrep / (nrep - 1));
        zHCovFull->SetBinContent(b1 + 1, b2 + 1, val);
      }
    }

    zHCovSys = zHCovFull;
  }
  else
  {
    const int nbins = zHNom->GetNbinsX();
    //double* bins = zHNom->GetXaxis()->GetXbins()->GetArray();

    // clear target histograms if Process us called many times
    if(zHTotalUD.first)
      delete zHTotalUD.first;
    if(zHTotalUD.second)
      delete zHTotalUD.second;
    if(zHSystUD.first)
      delete zHSystUD.first;
    if(zHSystUD.second)
      delete zHSystUD.second;
    if(zHCovSys)
      delete zHCovSys;
    if(zHCovFull)
      delete zHCovFull;
    if(zHCovSysEnvelope)
      delete zHCovSysEnvelope;
    if(zHCovFullEnvelope)
      delete zHCovFullEnvelope;

    // sum envelopes in quadrature to get total up and down uncertainties
    zHTotalUD.first = new TH1D(*zHNom);
    zHTotalUD.first->Reset();
    zHSystUD.first = new TH1D(*zHTotalUD.first);
    zHTotalUD.second = new TH1D(*zHNom);
    zHTotalUD.second->Reset();
    zHSystUD.second = new TH1D(*zHTotalUD.second);
    for(auto v : zHSys)
      utils::ttmd::AddVarInQuadrature(zHSystUD, v.second);

    if(DoEigenvectors)
    {
      if(!zHCovStat)
        throw std::logic_error("Error in ZMeasur::Process(): zHCovStat = NULL, DoQuadrature = true\n");
      // build covariance matrix from vars
      zHCovSys = new TH2D(*zHCovStat);
      zHCovSys->Reset();
      for(int i = 0; i < nbins; i++)
        for(int j = 0; j <= i; j++)
        {
          double val = zHCovSys->GetBinContent(i + 1, j + 1);
          for(auto v : zHEigens)
            val += v.second->GetBinContent(i + 1) * v.second->GetBinContent(j + 1);
          zHCovSys->SetBinContent(i + 1, j + 1, val);
          zHCovSys->SetBinContent(j + 1, i + 1, val);
        }
      zHCovFull = new TH2D(*zHCovSys);
      zHCovFull->Add(zHCovStat);

      //build "alternative" covariance matrix from envelopes
      zHCovSysEnvelope = new TH2D(*zHCovStat);
      zHCovSysEnvelope->Reset();
      for(int i = 0; i < nbins; i++)
        for(int j = 0; j <= i; j++)
        {
          double val = zHCovSysEnvelope->GetBinContent(i + 1, j + 1);
          for(auto v : zHSys)
          {
            TH1D* hu = v.second.first;
            TH1D* hd = v.second.second;
            val += (hu->GetBinContent(i + 1) * hu->GetBinContent(j + 1) + hd->GetBinContent(i + 1) * hd->GetBinContent(j + 1)) / 2.0;
          }
          zHCovSysEnvelope->SetBinContent(i + 1, j + 1, val);
          zHCovSysEnvelope->SetBinContent(j + 1, i + 1, val);
        }
      zHCovFullEnvelope = new TH2D(*zHCovSysEnvelope);
      zHCovFullEnvelope->Add(zHCovStat);
    }

    // total uncertainties (stat + systUD)
    for(int i = 0; i < nbins; i++)
    {
      double stat2 = 0.0;
      if(zHCovStat)
        stat2 = zHCovStat->GetBinContent(i + 1, i + 1);
      zHTotalUD.first->SetBinContent(i + 1, TMath::Sqrt(stat2 + TMath::Power(zHSystUD.first->GetBinContent(i + 1), 2.0)));
      zHTotalUD.second->SetBinContent(i + 1, -1 * TMath::Sqrt(stat2 + TMath::Power(zHSystUD.second->GetBinContent(i + 1), 2.0)));
    }
  }
}

TH1D* ZMeasur::GetNom() const
{
  assert(zHNom);
  return zHNom;
}

/*const TH1D*ZMeasur::GetVar(const TString& varName)
{
  assert(this->zHSys)
}*/

std::pair<TH1D*, TH1D*> ZMeasur::GetSystUD() const
{
  assert(zHSystUD.first);
  assert(zHSystUD.second);
  return zHSystUD;
}

std::pair<TH1D*, TH1D*> ZMeasur::GetTotalUD() const {
  assert(zHTotalUD.first);
  assert(zHTotalUD.second);
  return zHTotalUD;
}

TH2D* ZMeasur::GetCovFullMatrix() const
{
  assert(zHCovFull);
  return zHCovFull;
}

void ZMeasur::SetCovFullMatrix(TH2D* h)
{
  zHCovFull = h;
}

TH2D*ZMeasur::GetCovFullMatrixEnvelope() const
{
  assert(zHCovFullEnvelope);
  return zHCovFullEnvelope;
}

double ZMeasur::Chi2(const TH1D* h, int skip/* = 0*/) const
{
  double chi2 = utils::ttmd::Chi2(zHNom, zHCovFull, h, skip);
  return chi2;
}

double ZMeasur::Chi2(const TH1D* h, const std::vector<int>& vSkip) const
{
  double chi2 = utils::ttmd::Chi2(zHNom, zHCovFull, h, vSkip);
  return chi2;
}

void ZMeasur::AddVarToSys(const TString& sys, const TH1D* h)
{
  auto it = zHSys.find(sys);
  if(it == zHSys.end())
  {
    std::pair<TH1D*, TH1D*> hud;
    hud.first = new TH1D(*h);
    hud.first->Reset();
    hud.second = new TH1D(*h);
    hud.second->Reset();
    auto added = zHSys.insert(std::pair<TString, std::pair<TH1D*, TH1D*> >(sys, hud));
    it = added.first;
  }
  utils::ttmd::AddVarToEnvelope(it->second, h);
}
