#include "UnfoldingXSecHandler.h"
#include "ttmd_measur.h"
#include "ttmd_ZMadGraph.h"
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "ttmd_vars.h"
#include "utils.h"
#include <TPaveText.h>
#include <TGraphErrors.h>
#include "PlotterConfigurationHelper.h"
#include "UnfoldingIOHandler.h"
#include "TRandom3.h"
#include "TF1.h"
#include <cmath>

using namespace utils;


void UnfoldingXSecHandler::CheckResponseNormalisation(const UnfoldingIOHandler* unfoldingIOHandler)
{
  TH2D* hProb = utils::ttmd::MakeProbMatrix(unfoldingIOHandler->HRspTUnfold());
  double max = 0.0;
  for(int i = 1; i <= hProb->GetNbinsX(); i++)
  {
    double expect = 0.0;
    for(int j = 0; j <= hProb->GetNbinsY() + 1; j++)
      expect += unfoldingIOHandler->HGenC()->GetBinContent(j) * hProb->GetBinContent(i, j);
    double content = unfoldingIOHandler->HRecF()->GetBinContent(i);
    double diff = (content - expect) / expect;
    if(TMath::Abs(diff) > max)
      max = TMath::Abs(diff);
    //printf("bin %3d content = %6.1f expect %6.1f diff = %6.1e\n", i, content, expect, diff);
  }
  printf("UnfoldingXSecHandler::CheckResponseNormalisation() max relative deviation = %.1e\n", max);
  delete hProb;
}

std::vector<double> UnfoldingXSecHandler::GetPredictionImplPlain(const UnfoldingIOHandler* unfoldingIOHandler, ZMadGraph* mg, const ZPredSetup& predSetup) const
{
  const int nbins = unfoldingIOHandler->HGenC()->GetNbinsX();
  std::vector<double> result(nbins);

  assert(unfoldingIOHandler->VMGHisto.size());
  int bin = 0;
  for(size_t i = 0; i < unfoldingIOHandler->VMGHisto.size(); i++)
  {
    TH1D* h = ReadPredictionHistoSum(mg, predSetup, unfoldingIOHandler->VMGHisto[i].first, unfoldingIOHandler->VMGHisto[i].second);
    for(int j = 0; j < h->GetNbinsX(); j++)
    {
      double val = h->GetBinContent(j + 1);
      result[bin] = val;
      bin++;
    }
  }
  assert(bin == nbins);

  if(unfoldingIOHandler->VMGHistoMinus.size())
  {
    if(unfoldingIOHandler->VMGHistoMinus.size() == unfoldingIOHandler->VMGHisto.size())
    {
      int bin = 0;
      for(size_t i = 0; i < unfoldingIOHandler->VMGHistoMinus.size(); i++)
      {
        if(unfoldingIOHandler->VMGHistoMinus[i].first.size() == 0)
        {
          // skip bins, incrementing bin counter
          TH1D* h = ReadPredictionHistoSum(mg, predSetup, unfoldingIOHandler->VMGHisto[i].first, unfoldingIOHandler->VMGHisto[i].second);
          for(int j = 0; j < h->GetNbinsX(); j++)
            bin++;
        }
        else
        {
          TH1D* hMinus = ReadPredictionHistoSum(mg, predSetup, unfoldingIOHandler->VMGHistoMinus[i].first, unfoldingIOHandler->VMGHistoMinus[i].second);
          for(int j = 0; j < hMinus->GetNbinsX(); j++)
          {
            double val = hMinus->GetBinContent(j + 1);
            result[bin] -= val;
            bin++;
          }
        }
      }
      assert(bin == nbins);
    }
    else
      throw std::logic_error(TString::Format("Error: VHistoOrderMadGraph2 = %zu != %zu = VHistoOrderMadGraphMinus2", unfoldingIOHandler->VMGHisto.size(), unfoldingIOHandler->VMGHistoMinus.size()));
  }

  return result;
}

std::vector<double> UnfoldingXSecHandler::GetPredictionImplIntegrate(const UnfoldingIOHandler* unfoldingIOHandler, ZMadGraph* mg, const ZPredSetup& predSetup) const
{
  const int nbins = unfoldingIOHandler->HGenC()->GetNbinsX();
  std::vector<double> result(nbins);

  assert(unfoldingIOHandler->VMGHisto.size());
  int bin = 0;
  for(size_t i = 0; i < unfoldingIOHandler->VMGHisto.size(); i++)
  {
    TH1D* h = ReadPredictionHistoSum(mg, predSetup, unfoldingIOHandler->VMGHisto[i].first, unfoldingIOHandler->VMGHisto[i].second);
    double val = h->Integral();
    result[bin] = val;
    bin++;
  }
  assert(bin == nbins);

  if(unfoldingIOHandler->VMGHistoMinus.size())
  {
    if(unfoldingIOHandler->VMGHistoMinus.size() == unfoldingIOHandler->VMGHisto.size())
    {
      int bin = 0;
      for(size_t i = 0; i < unfoldingIOHandler->VMGHistoMinus.size(); i++)
      {
        if(unfoldingIOHandler->VMGHistoMinus[i].first.size() == 0)
        {
          // skip bins, incrementing bin counter
          //TH1D* h = ReadPredictionHistoSum(mg, predSetup, unfoldingIOHandler->VMGHisto[i].first, unfoldingIOHandler->VMGHisto[i].second);
          bin++;
        }
        else
        {
          TH1D* hMinus = ReadPredictionHistoSum(mg, predSetup, unfoldingIOHandler->VMGHistoMinus[i].first, unfoldingIOHandler->VMGHistoMinus[i].second);
          double val = hMinus->Integral();
          result[bin] -= val;
          bin++;
        }
      }
      assert(bin == nbins);
    }
    else
      throw std::logic_error(TString::Format("Error: unfoldingIOHandler->VMGHisto = %zu != %zu = unfoldingIOHandler->VMGHistoMinus", unfoldingIOHandler->VMGHisto.size(), unfoldingIOHandler->VMGHistoMinus.size()));
  }

  return result;
}

std::vector<double> UnfoldingXSecHandler::GetPredictionImplPlainFast(const UnfoldingIOHandler* unfoldingIOHandler, ZMadGraph* mg, const ZPredSetup& predSetup) const
{
  const int nbins = unfoldingIOHandler->HGenC()->GetNbinsX();
  std::vector<double> result(nbins);

  assert(unfoldingIOHandler->VMGHisto.size());
  int bin = 0;
  for(size_t i = 0; i < unfoldingIOHandler->VMGHisto.size(); i++)
  {
    // 25.06.18 sanity check
    //assert(unfoldingIOHandler->VMGHisto[i].first.size() == 1);
    int nbinsCurrent = 0;
    float* pred = ReadPredictionHistoSumFast(mg, predSetup, unfoldingIOHandler->VMGHisto[i].first, unfoldingIOHandler->vMGHistoOneBinOnly[i], unfoldingIOHandler->VMGHisto[i].second, &nbinsCurrent);
    assert(nbinsCurrent);
    for(int j = 0; j < nbinsCurrent; j++)
    {
      result[bin] = pred[j];
      bin++;
    }
  }
  assert(bin == nbins);

  if(unfoldingIOHandler->VMGHistoMinus.size())
  {
    if(unfoldingIOHandler->VMGHistoMinus.size() == unfoldingIOHandler->VMGHisto.size())
    {
      int bin = 0;
      for(size_t i = 0; i < unfoldingIOHandler->VMGHistoMinus.size(); i++)
      {
        if(unfoldingIOHandler->VMGHistoMinus[i].first.size() == 0)
        {
          // skip bins, incrementing bin counter
          int nbinsCurrent = 0;
          //float* pred = ReadPredictionHistoSumFast(mg, predSetup, unfoldingIOHandler->VMGHisto[i].first, unfoldingIOHandler->vMGHistoOneBinOnly[i], unfoldingIOHandler->VMGHisto[i].second, &nbinsCurrent);
          assert(nbinsCurrent);
          bin += nbinsCurrent;
        }
        else
        {
          // 25.06.18 sanity check
          //assert(unfoldingIOHandler->VMGHistoMinus[i].first.size() == 1);
          int nbinsCurrent = 0;
          float* pred = ReadPredictionHistoSumFast(mg, predSetup, unfoldingIOHandler->VMGHistoMinus[i].first, unfoldingIOHandler->vMGHistoOneBinOnly[i], unfoldingIOHandler->VMGHistoMinus[i].second, &nbinsCurrent);
          assert(nbinsCurrent);
          for(int j = 0; j < nbinsCurrent; j++)
          {
            result[bin] -= pred[j];
            bin++;
          }
        }
      }
      assert(bin == nbins);
    }
    else
      throw std::logic_error(TString::Format("Error: VHistoOrderMadGraph2 = %zu != %zu = VHistoOrderMadGraphMinus2", unfoldingIOHandler->VMGHisto.size(), unfoldingIOHandler->VMGHistoMinus.size()));
  }

  return result;
}

std::vector<double> UnfoldingXSecHandler::GetPredictionImplIntegrateFast(const UnfoldingIOHandler* unfoldingIOHandler, ZMadGraph* mg, const ZPredSetup& predSetup) const
{
  const int nbins = unfoldingIOHandler->HGenC()->GetNbinsX();
  std::vector<double> result(nbins);

  assert(unfoldingIOHandler->VMGHisto.size());
  int bin = 0;
  for(size_t i = 0; i < unfoldingIOHandler->VMGHisto.size(); i++)
  {
    int nbinsCurrent = 0;
    float* res = ReadPredictionHistoSumFast(mg, predSetup, unfoldingIOHandler->VMGHisto[i].first, unfoldingIOHandler->vMGHistoOneBinOnly[i], unfoldingIOHandler->VMGHisto[i].second, &nbinsCurrent);
    assert(nbinsCurrent);
    double val = 0.0;
    for(int b = 0; b < nbinsCurrent; b++)
      val += res[b];
    result[bin] = val;
    bin++;
  }
  assert(bin == nbins);

  if(unfoldingIOHandler->VMGHistoMinus.size())
  {
    if(unfoldingIOHandler->VMGHistoMinus.size() == unfoldingIOHandler->VMGHisto.size())
    {
      int bin = 0;
      for(size_t i = 0; i < unfoldingIOHandler->VMGHistoMinus.size(); i++)
      {
        if(unfoldingIOHandler->VMGHistoMinus[i].first.size() == 0)
        {
          // skip bins, incrementing bin counter
          //TH1D* h = ReadPredictionHistoSum(mg, predSetup, unfoldingIOHandler->VMGHisto[i].first, unfoldingIOHandler->VMGHisto[i].second);
          bin++;
        }
        else
        {
          int nbinsCurrent = 0;
          float* pred = ReadPredictionHistoSumFast(mg, predSetup, unfoldingIOHandler->VMGHistoMinus[i].first, unfoldingIOHandler->vMGHistoOneBinOnly[i], unfoldingIOHandler->VMGHistoMinus[i].second, &nbinsCurrent);
          assert(nbinsCurrent);
          double val = 0.0;
          for(int b = 0; b < nbinsCurrent; b++)
            val += pred[b];
          result[bin] -= val;
          bin++;
        }
      }
      assert(bin == nbins);
    }
    else
      throw std::logic_error(TString::Format("Error: unfoldingIOHandler->VMGHisto = %zu != %zu = unfoldingIOHandler->VMGHistoMinus", unfoldingIOHandler->VMGHisto.size(), unfoldingIOHandler->VMGHistoMinus.size()));
  }

  return result;
}

TH1D* UnfoldingXSecHandler::ReadPredictionHistoSum(ZMadGraph* mg, const ZPredSetup& predSetup, const std::vector<int>& vHistoID, const int order)
{
  TH1D* h = NULL;
  for(size_t ii = 0; ii < vHistoID.size(); ii++)
  {
    ZPredHisto predHisto(predSetup, vHistoID[ii], order);
    const TH1D* hCurrent = mg->ReadPred(predHisto);
    if(!h)
      h = (TH1D*) hCurrent->Clone();
    else
      h->Add(hCurrent);
  }
  //h->Print("all");
  return h;
}

float* UnfoldingXSecHandler::ReadPredictionHistoSumFast(ZMadGraph* mg, const ZPredSetup& predSetup, const std::vector<int>& vHistoID, const int& oneBin, const int order, int* nbinsPtr)
{
  static float res[ZMadGraph_MAXHISTONBINS]; // static in order to exist after leaving this scope
  int nbins = 0;
  //int nbinsLast = 0;
  for(size_t ii = 0; ii < vHistoID.size(); ii++)
  {
    ZPredHisto predHisto(predSetup, vHistoID[ii], order);
    const float* resCurrent = mg->ReadPredFast(predHisto, &nbins);
    assert(nbins);
    assert(ii == 0 || nbins == nbinsLast);
    if(ii == 0)
    {
      for(int b = 0; b < nbins; b++)
      {
        if(oneBin == -1 || oneBin == b)
          res[b] = resCurrent[b];
        else
          res[b] = 0.0;
      }
    }
    else
    {
      if(oneBin >= 0)
        res[oneBin] += resCurrent[oneBin];
      else
        for(int b = 0; b < nbins; b++)
          if(oneBin == -1 || oneBin == b)
            res[b] += resCurrent[b];
    }
    assert(nbinsLast = nbins);
  }
  if(nbinsPtr)
    *nbinsPtr = nbins;
  return res;
}

//void UnfoldingXSecHandler::PutPredictionPttt(double ptThreshold)
//{
//  int nptttAll = ZMadGraph::GetOffsetHistoPtttThreshold(configHelper_, -1.0);
//  int npttt = ZMadGraph::GetOffsetHistoPtttThreshold(configHelper_, ptThreshold);
//  std::vector<int> hid[2];
//  for(int p = 0; p < nptttAll; p++)
//  {
//    if(p < npttt)
//      hid[0].push_back(100 + 4 * (p + 1));
//    else
//      hid[1].push_back(100 + 4 * (p + 1));
//  }
//  // 1st bin pttt
//  for(int i = 1; i <= 4; i++)
//  {
//    std::vector<int> vh;
//    for(size_t ii = 0; ii < hid[0].size(); ii++)
//      vh.push_back(hid[0][ii] + i);
//    for(size_t ii = 0; ii < hid[1].size(); ii++)
//      vh.push_back(hid[1][ii] + i);
//    VMGHisto.push_back(std::pair<std::vector<int>, int>(vh, 0));
//    vh.clear();
//    for(size_t ii = 0; ii < hid[1].size(); ii++)
//      vh.push_back(hid[1][ii] + i);
//    VMGHistoMinus.push_back(std::pair<std::vector<int>, int>(vh, 1));
//  }
//  // 2nd bin pttt
//  for(int i = 1; i <= 4; i++)
//  {
//    std::vector<int> vh;
//    for(size_t ii = 0; ii < hid[1].size(); ii++)
//      vh.push_back(hid[1][ii] + i);
//    VMGHisto.push_back(std::pair<std::vector<int>, int>(vh, 1));
//    VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));
//  }
//}
//
//void UnfoldingXSecHandler::PutPredictionPtttIntegrate(double ptThreshold)
//{
//  int nptttAll = ZMadGraph::GetOffsetHistoPtttThreshold(configHelper_, -1.0);
//  int npttt = ZMadGraph::GetOffsetHistoPtttThreshold(configHelper_, ptThreshold);
//  std::vector<int> hid[2];
//  for(int p = 0; p < nptttAll; p++)
//  {
//    if(p < npttt)
//      hid[0].push_back(100 + 4 * (p + 1));
//    else
//      hid[1].push_back(100 + 4 * (p + 1));
//  }
//  // 1st bin pttt
//  std::vector<int> vh;
//  for(int i = 1; i <= 4; i++)
//  {
//    for(size_t ii = 0; ii < hid[0].size(); ii++)
//      vh.push_back(hid[0][ii] + i);
//    for(size_t ii = 0; ii < hid[1].size(); ii++)
//      vh.push_back(hid[1][ii] + i);
//  }
//  VMGHisto.push_back(std::pair<std::vector<int>, int>(vh, 0));
//  vh.clear();
//  for(int i = 1; i <= 4; i++)
//  {
//    for(size_t ii = 0; ii < hid[1].size(); ii++)
//      vh.push_back(hid[1][ii] + i);
//  }
//  VMGHistoMinus.push_back(std::pair<std::vector<int>, int>(vh, 1));
//  // 2nd bin pttt
//  vh.clear();
//  for(int i = 1; i <= 4; i++)
//  {
//    for(size_t ii = 0; ii < hid[1].size(); ii++)
//      vh.push_back(hid[1][ii] + i);
//  }
//  VMGHisto.push_back(std::pair<std::vector<int>, int>(vh, 1));
//  vh.clear();
//  VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));
//}

std::vector<double> UnfoldingXSecHandler::GetPrediction(const UnfoldingIOHandler* unfoldingIOHandler, ZMadGraph* mg, const ZPredSetup& predSetup, const bool flagNorm, const bool flagDivideBW) const
{
  std::vector<double> result;
  if(!(unfoldingIOHandler->FlagMGHistoIntegrate))
    //result = unfoldingIOHandler->GetPredictionImplPlain(unfoldingIOHandler, mg, predSetup);
    result = GetPredictionImplPlainFast(unfoldingIOHandler, mg, predSetup);
  else
    //result = unfoldingIOHandler->GetPredictionImplIntegrate(unfoldingIOHandler, mg, predSetup);
    result = GetPredictionImplIntegrateFast(unfoldingIOHandler, mg, predSetup);

  // consistency check
  for(const auto& val : result){
      if (val <= 0.0) {
          std::cerr << "ERROR in UnfoldingXSecHandler::GetPrediction -->> value <= 0 was found!!" << std::endl;
          exit(1);
      }
  }

  if(unfoldingIOHandler->FlagNeedNP)
  {
    TString fileName = mg->DirNPCorr + "/" + unfoldingIOHandler->Suffix + "/" + "npcorr.dat";
    TH1D* h = unfoldingIOHandler->GetPlainHisto();
    utils::ttmd::ReadHistoKFactor(h, fileName, unfoldingIOHandler->Dim() * 2 + 1 + 1, 1);
    const int nbins = unfoldingIOHandler->HGenC()->GetNbinsX();
    assert(h->GetNbinsX() == nbins);
    for(int b = 0; b < nbins; b++)
      result[b] *= h->GetBinContent(b + 1);
  }

  if(flagNorm)
  {
    double sum = 0.0;
    for(auto& val : result)
      sum += val;
    for(auto& val : result)
      val /= sum;
  }

  if(flagDivideBW)
  {
    ZVar* varLastDim = unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - 1);
    const std::vector<double>& binsLastDim = varLastDim->BinsC();
    const int nbinsLastDim = binsLastDim.size() - 1;
    for(size_t b = 0; b < result.size(); b++)
    {
      int binLastDim = b % nbinsLastDim;
      double binWidth = binsLastDim[binLastDim + 1] - binsLastDim[binLastDim];
      result[b] /= binWidth;
    }
  }

  return result;
}

TH1D* UnfoldingXSecHandler::GetPredictionHisto(const UnfoldingIOHandler* unfoldingIOHandler, ZMadGraph* mg, const ZPredSetup& predSetup, const bool flagNorm, const bool flagDivideBW) const
{
  TH1D* h = (TH1D*) unfoldingIOHandler->HGenC()->Clone();
  //return h;
  std::vector<double> vPred = GetPrediction(unfoldingIOHandler, mg, predSetup, flagNorm, flagDivideBW);
  for(int b = 0; b < h->GetNbinsX(); b++)
    h->SetBinContent(b + 1, vPred[b]);
  if(unfoldingIOHandler->DoRj)
    utils::ttmd::DoRjHisto(h);
  return h;
}

double UnfoldingXSecHandler::d1_st1(double* xx, double* par){
  // parameters
  double min = par[0];
  double max = par[1];
  double a = par[2];
  double b = par[3];
  double c = par[4];
  double d = par[5];
  double e = par[6];
  double f = par[7];

  // variable
  double x = (xx[0] - min) / (max - min);

  const double eps = 1e-6;
  x = std::max(x, eps);
  x = std::min(x, 1.0 - eps);

  // function
  double eval = a * TMath::Power(x, b) * TMath::Power(1 - x, c) * (1 + d * x + e * x * x + f * x * x * x);
  //printf("x = %f  d1_st1(%f) = %f\n", x, xx[0], eval);
  return eval;
}

double UnfoldingXSecHandler::d2_st1_mult(double* xx, double* par)
{
  double v1 = d1_st1(xx, par);
  double v2 = d1_st1(xx + 1, par + 8);
  return v1 * v2;
}

double UnfoldingXSecHandler::d3_st1_mult(double* xx, double* par)
{
  double v1 = d1_st1(xx, par);
  double v2 = d1_st1(xx + 1, par + 8);
  double v3 = d1_st1(xx + 2, par + 16);
  return v1 * v2 * v3;
}

TF1* UnfoldingXSecHandler::GetRewFunction(const UnfoldingIOHandler* unfoldingIOHandler, const std::vector<std::vector<double> >& vvPar)
{
  assert(unfoldingIOHandler->Dim() <= vvPar.size());
  const int dim = unfoldingIOHandler->Dim();
  TF1* fRew = NULL;
  if(dim == 1)
    fRew = new TF1("", UnfoldingXSecHandler::d1_st1, 0.0, 1.0, 8);
  else if(dim == 2)
    fRew = new TF2("", UnfoldingXSecHandler::d2_st1_mult, 0.0, 1.0, 0.0, 1.0, 16);
  else if(dim == 3)
    fRew = new TF3("", UnfoldingXSecHandler::d3_st1_mult, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 24);
  else
    throw std::logic_error(TString::Format("Error: unfoldingIOHandler->Dim() = %d\n", unfoldingIOHandler->Dim()).Data());

  double width = 1.0;
  for(int d = 0; d < unfoldingIOHandler->Dim(); d++)
  {
    const ZVar* var = unfoldingIOHandler->Var(d);
    double xmin = var->BinsC()[0];
    double xmax = var->BinsC()[var->BinsC().size() - 1];
    printf("dim: %d  xmin, xmax: %f %f\n", d, xmin, xmax);
    width *= (xmax - xmin);
    //xmin -= ratioExpand * (xmax - xmin);
    //xmax += ratioExpand * (xmax - xmin);
    const int offset = d * 8;
    fRew->SetParameter(0 + offset, xmin);
    fRew->SetParameter(1 + offset, xmax);
    fRew->SetParameter(2 + offset, 1.0);
    fRew->SetParameter(3 + offset, vvPar[d][0]);
    fRew->SetParameter(4 + offset, vvPar[d][1]);
    fRew->SetParameter(5 + offset, vvPar[d][2]);
    fRew->SetParameter(6 + offset, vvPar[d][3]);
    fRew->SetParameter(7 + offset, vvPar[d][4]);
  }

  double norm = unfoldingIOHandler->HGenC()->Integral();
  double integral = -1.0;
  if(dim == 1)
    integral = fRew->Integral(fRew->GetParameter(0), fRew->GetParameter(1));
  else if(dim == 2)
    integral = ((TF2*)fRew)->Integral(fRew->GetParameter(0), fRew->GetParameter(1), fRew->GetParameter(8), fRew->GetParameter(9));
  else if(dim == 3)
    integral = ((TF3*)fRew)->Integral(fRew->GetParameter(0), fRew->GetParameter(1), fRew->GetParameter(8), fRew->GetParameter(9), fRew->GetParameter(16), fRew->GetParameter(17));
  fRew->SetParameter(2, norm / width / integral);
  return fRew;
}

TF1* UnfoldingXSecHandler::GetRewFunction(const UnfoldingIOHandler* unfoldingIOHandler, const int grade)
{
  static std::vector<std::vector<std::vector<double> > > vvvPar; // [grade][dim][bcdef]
  if(vvvPar.size() == 0)
  {
    //vvvPar.push_back({ { 0.02, 0.003, 0.03, -0.01, -0.003}, { 0.005,  0.004, 0.001, -0.04, -0.02}, { 0.002,  0.003, 0.003, -0.01, 0.02} });
    for(int i = 0; i < 4; i++)
    {
      double k = -1.0;
      if(i == 0)
        k = 1.0;
      else if(i == 1)
        k = 5.0;
      else if(i == 2)
        k = 15.0;
      else if(i == 3)
        k = 50.0;
      assert(k > 0.0);
      vvvPar.push_back({ { 0.020 * k, 0.003 * k, 0.030 * k, -0.01 * k, -0.003 * k},
                         { 0.025 * k, 0.004 * k, 0.001 * k, -0.03 * k,  0.04 * k},
                         { 0.002 * k, 0.003 * k, 0.015 * k, -0.01 * k,  0.02 * k}
                       });
    }
    //vvvPar.push_back({ { 0.55,  0.08,  0.9, -0.2, -0.08}, { 0.10,  0.02, 0.3, -1.0, -0.30}, { 0.10, 0.06,  0.1, 0.08, -0.03} });
  }
  assert(grade <= vvvPar.size());
  TF1* f = GetRewFunction(unfoldingIOHandler, vvvPar[grade - 1]);
  return f;
}
