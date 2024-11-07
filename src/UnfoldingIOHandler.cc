#include <fstream>
#include "UnfoldingIOHandler.h"
//#include "measur.h"
#include "ttmd_ZMadGraph.h"
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TObjString.h>
#include "ttmd_vars.h"
#include "utils.h"
#include <TPaveText.h>
#include <TGraphErrors.h>
#include "ttmd_framework.h"
#include "PlotterConfigurationHelper.h"
//#include "ttmd_settings.h"
#include "TRandom3.h"
#include "TF1.h"
#include <cmath>
#include "UnfoldingXSecHandler.h"

using namespace utils;

UnfoldingIOHandler::UnfoldingIOHandler(const std::vector<ZVar*>& argVVar, PlotterConfigurationHelper *configHelper, const TString& suffix/* = ""*/)
{
  zVVar = argVVar;
  zFRew = NULL;
  zvVFRew = NULL;
  if(suffix != "")
    Suffix = suffix;
  for(unsigned int i = 0; i < argVVar.size(); i++)
    Expression += argVVar[i]->Expression;
  configHelper_= configHelper;
}

UnfoldingIOHandler::UnfoldingIOHandler(const UnfoldingIOHandler& old, PlotterConfigurationHelper *configHelper)
{
  Suffix = old.Suffix;
  //zVVar = old.zVVar;
  zVVar.resize(old.zVVar.size());
  for(size_t v = 0; v < zVVar.size(); v++)
    zVVar[v] = old.zVVar[v]->clone();
  Expression = old.Expression;
  zFRew = NULL;
  zvVFRew = NULL;
  DoXSec = old.DoXSec;
  DoYields = old.DoYields;
  DoPSE = old.DoPSE;
  FlagNeedNP = old.FlagNeedNP;
  NPertOrders = old.NPertOrders;
  configHelper_= configHelper;
}

void UnfoldingIOHandler::PrintRecWeights(const TString& ch, TString in_Var, int kin_reco_mode) const
{
  double dat = this->GetRecWeight(0);
  double mcsig = this->GetRecWeight(1);
  double mcoth = this->GetRecWeight(2);
  double mcbgr = this->GetRecWeight(3);
  double mcbgr_singletop = this->GetRecWeight(4);
  //double mcbgr_zjets = this->GetRecWeight(5);
  double mcbgr_zjets1050 = this->GetRecWeight(9);
  double mcbgr_zjets50inf = this->GetRecWeight(5);
  double mccbgr_zjets = mcbgr_zjets1050 + mcbgr_zjets50inf;
  double mcbgr_wjets = this->GetRecWeight(6);
  double mcbgr_diboson = this->GetRecWeight(7);
  double mcbgr_ttbarWZg = this->GetRecWeight(8);
  //double check = mcbgr_singletop + mcbgr_zjets + mcbgr_wjets + mcbgr_diboson + mcbgr_ttbarWZg;
  double mctot = mcsig + mcoth + mcbgr;
  printf("Data %.0f [MC sum %.0f ttbar sig %.0f(%.1f%%) ttbar oth %.0f(%.1f%%) bgr %.0f(%.1f%%)]\n",
         dat, mctot, mcsig, mcsig / mctot * 100, mcoth, mcoth / mctot * 100, mcbgr, mcbgr / mctot * 100);
  //printf("bgr [singletop %.0f(%.1f%%) zjets %.0f(%.1f%%) wjets %.0f(%.1f%%) diboson %.0f(%.1f%%) ttbarWZg %.0f(%.1f%%)]\n",
  //       mcbgr_singletop, mcbgr_singletop / dat * 100, mcbgr_zjets, mcbgr_zjets / dat * 100,
  //       mcbgr_wjets, mcbgr_wjets / dat * 100, mcbgr_diboson, mcbgr_diboson / dat * 100, mcbgr_ttbarWZg, mcbgr_ttbarWZg / dat * 100);
  printf("bgr [singletop %.0f(%.1f%%) zjets1050 %.0f(%.1f%%) zjets50inf %.0f(%.1f%%) wjets %.0f(%.1f%%) diboson %.0f(%.1f%%) ttbarWZg %.0f(%.1f%%)]\n",
         mcbgr_singletop, mcbgr_singletop / dat * 100, mcbgr_zjets1050, mcbgr_zjets1050 / dat * 100, mcbgr_zjets50inf, mcbgr_zjets50inf / dat * 100,
         mcbgr_wjets, mcbgr_wjets / dat * 100, mcbgr_diboson, mcbgr_diboson / dat * 100, mcbgr_ttbarWZg, mcbgr_ttbarWZg / dat * 100);


  //summary txt file for Nominal
  if (in_Var == "Nominal"){
      std::map<TString, double> samples_info_map {
          {"Data", dat},
          {"ttbarWZg", mcbgr_ttbarWZg},
          {"Z+", mccbgr_zjets},
          {"Diboson", mcbgr_diboson},
          {"W+", mcbgr_wjets},
          {"Single", mcbgr_singletop},
          {"other", mcoth},
          {"signal", mcsig},
          {"Total", mctot}
      };

      TString outdir = this->configHelper_->gAnalDir + "/";
      TString txt_suffix = "";
      if (kin_reco_mode == 0) txt_suffix = "9";
      else if (kin_reco_mode == 9) txt_suffix = "7L";
      else {
          std::cerr << "ERROR in UnfoldingIOHandler::PrintRecWeights ---->>> Kin. Reco. mode " << kin_reco_mode << " not implemented in this function!" << std::endl;
          exit(1);
      }
      TString txt_file_name = ch + "_eventCount"+ txt_suffix +".txt";
      std::ofstream txt_EventCount_file;
      utils::ttmd::CheckFile(outdir);
      txt_EventCount_file.open(outdir + txt_file_name);

      for (auto sample_info: samples_info_map){
          txt_EventCount_file << sample_info.first << ": " << sample_info.second << std::endl;
      }
      txt_EventCount_file.close();
      std::cout << "Event yields saved in " << outdir + txt_file_name << std::endl;
  }


}

void UnfoldingIOHandler::SetPseudoData(const std::vector<TH1D*>& v)
{
  size_t NBinnings_unsigned = NBinnings();
  if(v.size() != NBinnings_unsigned)
  {
    printf("Error in UnfoldingIOHandler::SetPseudoData() inconsistent number of binnings %ld for input vector, %d for var\n", v.size(), NBinnings());
    throw;
  }
  if(zHDat.size() != 0)
  {
    zHDat.clear();
    printf("Warning in UnfoldingIOHandler::SetPseudoData() zHDat.size() = %ld will be overwritten\n", zHDat.size());
  }
  for(unsigned int i = 0; i < v.size(); i++)
    zHDat.push_back((TH1D*) v[i]->Clone());
}

void UnfoldingIOHandler::GeneratePseudoData(const int binning, const TH1D* hGen, const TH2D* hRsp,
                               const TH1D* hRefUnc/* = NULL*/, const bool flagFluct/* = true*/)
{
  /*if(zHDat.size() != 0)
  {
    zHDat.clear();
    printf("Warning in UnfoldingIOHandler::GeneratePseudoData() zHDat.size() = %ld will be overwritten\n", zHDat.size());
  }*/

  TH2D* hProb = utils::ttmd::MakeProbMatrix(hRsp);
  TH1D* hRec = (TH1D*)zTUnfBinning[binning]->CreateHistogram(TString::Format("h_rec_binning%d", binning));
  TRandom3& rnd = ZRandom::GetInstance();
  //if(gTTBARDebug)
  //  printf("Info in UnfoldingIOHandler::GeneratePseudoDataF() start seed %d\n", rnd.GetSeed());
  for(int i = 1; i <= hProb->GetNbinsX(); i++)
  {
    double expect = 0.0;
    for(int j = 0; j <= hProb->GetNbinsY() + 1; j++)
      expect += hGen->GetBinContent(j) * hProb->GetBinContent(i, j);
    double fluct = 0.0;
    if(flagFluct)
    {
      double mean = hRefUnc ? hRefUnc->GetBinContent(i) : expect;
      fluct = rnd.PoissonD(mean) - mean;
      //if(i == 2) printf("expect = %f  mean = %f  fluct = %f\n", expect, mean, fluct);
    }
    hRec->SetBinContent(i, expect + fluct);
    hRec->SetBinError(i, TMath::Sqrt(expect + fluct));
  }
  //zHDat.push_back(hRec);
  if(zHDat[binning])
    delete zHDat[binning];
  zHDat[binning] = hRec;
  //hRec->Print("all");
  //if(gTTBARDebug)
  //  printf("Info in UnfoldingIOHandler::GeneratePseudoDataF() finish next seed %d\n", rnd.GetSeed());
  delete hProb;
  //delete hRec;
}

TH2D* UnfoldingIOHandler::HRspTUnfold11() const
{
  TH2D* hRsp = (TH2D*) zHRsp[2]->Clone();
  for(int b = 0; b <= hRsp->GetNbinsY() + 1; b++)
  {
    // check that events in reco are the same as column integrals in response matrix
    /*double rec1 = zHRec[1]->GetBinContent(b);
    double rec2 = zHRsp[0]->Integral(1, zHRsp[0]->GetNbinsX(), b, b);
    if(IsEqual(rec1, rec2))
    {
    printf("mismatch rec1 %f != %f rec2, rel. diff. %f\n", rec1, rec2, (rec1 - rec2) / rec1);
    }*/
    double genANDrec = hRsp->Integral(0, hRsp->GetNbinsX() + 1, b, b);
    // neglect uncertainty
    double genANDrecUnc = 0.0;//hRsp->IntegralError(0, hRsp->GetNbinsX() + 1, b, b);
    double genNOTrec = zHGen[0]->GetBinContent(b) - genANDrec;
    // assume that gen = reco + genNOTrec -> unc(genNOTrec) = sqrt(unc(gen)^2-unc(rec)^2)
    double genNOTrecUnc = 0.0;
    if(zHGen[0]->GetBinError(b) > zHRec[0]->GetBinError(b))
      genNOTrecUnc = TMath::Sqrt(TMath::Power(zHGen[0]->GetBinError(b), 2.0) - TMath::Power(genANDrecUnc, 2.0));
    double oldContent = hRsp->GetBinContent(0, b);
    double oldUnc = hRsp->GetBinError(0, b);
    double newContent = oldContent + genNOTrec;
    // assume that rec and genNOTrec are uncorrelated
    double newUnc = TMath::Sqrt(TMath::Power(oldUnc, 2.0) + TMath::Power(genNOTrecUnc, 2.0));
    //if(b == 1) printf("bx = 1  %f = %f - %f + %f\n", newContent, zHGen[1]->GetBinContent(b), zHRec[1]->GetBinContent(b), hRsp->GetBinContent(0, b));
    hRsp->SetBinContent(0, b, newContent);
    hRsp->SetBinError(0, b, newUnc);
  }
  return hRsp;
}

void UnfoldingIOHandler::Init()
{
  // 25.06.18 imlemented option to take one bin only from all hisotgrams
  // put -1 if it is empty (init() call must be after providing MG histograms)
  if(vMGHistoOneBinOnly.size() == 0)
    vMGHistoOneBinOnly = std::vector<int>(VMGHisto.size(), -1);

  // determine number of binnings, check that it is consistent for all vars
  int nBinnings0 = 0;
  for(int v = 0; v < Dim(); v++)
  {
    if(v == 0)
      nBinnings0 = zVVar[v]->NBinnings();
    else if(nBinnings0 != zVVar[v]->NBinnings())
    {
      printf("Error: inconsistent number of binnings %d for 0 var, %d for %d var\n", nBinnings0, zVVar[v]->NBinnings(), v);
      throw;
    }
  }

  // create TUnfoldBinning instance for each binning
  // need to clear if this object was created by copying
  zTUnfBinning.clear();
  zTUnfBinningNode.clear();
  for(int i = 0; i < nBinnings0; i++)
  {
    zTUnfBinning.push_back(new TUnfoldBinning(TString::Format("binning%d", i)));
    TUnfoldBinning* last = zTUnfBinning.back();
    zTUnfBinningNode.push_back(last->AddBinning("dummy"));
    TUnfoldBinning* ptr = zTUnfBinningNode.back();
    // backward order of variables to keep last dimension variable smooth
    //for(int v = 0; v < Dim(); v++)
    for(int v = Dim() - 1; v >= 0; v--)
    {
      const std::vector<double>& bins = zVVar[v]->Bins[i];
      ptr->AddAxis(zVVar[v]->Expression, bins.size() - 1, &bins[0], false, false);
    }
  }

  // create 1D histograms
  //TString hNameBase = "h";
  //for(int v = 0; v < Dim(); v++)
  //  hNameBase += "_" + zVVar[v]->Expression;
  zHDat.clear();
  zHRec.clear();
  zHGen.clear();
  zHBgr.clear();
  zHBtt.clear();
  zHOutOfSpace.clear();
  for(int i = 0; i < NBinnings(); i++)
  {
    // data
    zHDat.push_back((TH1D*)zTUnfBinning[i]->CreateHistogram(""));
    // MC signal reco
    zHRec.push_back((TH1D*)zTUnfBinning[i]->CreateHistogram(""));
    // MC signal gen
    zHGen.push_back((TH1D*)zTUnfBinning[i]->CreateHistogram(""));
    // MC non-ttbar background
    zHBgr.push_back((TH1D*)zTUnfBinning[i]->CreateHistogram(""));
    // MC ttbar background
    zHBtt.push_back((TH1D*)zTUnfBinning[i]->CreateHistogram(""));
    // MC ttbar out of space background
    zHOutOfSpace.push_back((TH1D*)zTUnfBinning[i]->CreateHistogram(""));
    /*// data
    zHDat.push_back((TH1D*)zTUnfBinning[i]->CreateHistogram(TString::Format(hNameBase + "_dat_binning%d", i)));
    // MC signal reco
    zHRec.push_back((TH1D*)zTUnfBinning[i]->CreateHistogram(TString::Format(hNameBase + "_rec_binning%d", i)));
    // MC signal gen
    zHGen.push_back((TH1D*)zTUnfBinning[i]->CreateHistogram(TString::Format(hNameBase + "_gen_binning%d", i)));
    // MC non-ttbar background
    zHBgr.push_back((TH1D*)zTUnfBinning[i]->CreateHistogram(TString::Format(hNameBase + "_bgr_binning%d", i)));
    // MC ttbar background
    zHBtt.push_back((TH1D*)zTUnfBinning[i]->CreateHistogram(TString::Format(hNameBase + "_btt_binning%d", i)));*/
  }

  // create TH2 response matrix
  //zHRsp.push_back(TUnfoldBinning::CreateHistogramOfMigrations(zTUnfBinning[0], zTUnfBinning[0], hNameBase + "_rsp_v0"));
  zHRsp.clear();
  zHRsp.push_back(TUnfoldBinning::CreateHistogramOfMigrations(zTUnfBinning[0], zTUnfBinning[0], ""));
  if(NBinnings() > 1)
  {
    //zHRsp.push_back(TUnfoldBinning::CreateHistogramOfMigrations(zTUnfBinning[1], zTUnfBinning[0], hNameBase + "_rsp_v1"));
    //zHRsp.push_back(TUnfoldBinning::CreateHistogramOfMigrations(zTUnfBinning[1], zTUnfBinning[0], hNameBase + "_rsp_v11"));
    zHRsp.push_back(TUnfoldBinning::CreateHistogramOfMigrations(zTUnfBinning[1], zTUnfBinning[0], ""));
    zHRsp.push_back(TUnfoldBinning::CreateHistogramOfMigrations(zTUnfBinning[1], zTUnfBinning[0], ""));
  }

  // create vector with weights
  zRecWeight.resize(11);

  // set suffix if not set already
  if(Suffix == "")
    Suffix = Expression;
}

void UnfoldingIOHandler::Reset(const TString& opt)
{
  for(int i = 0; i < NBinnings(); i++)
  {
    // data
    if(opt.First('d') != kNPOS)
    {
      zHDat[i]->Reset(); // "d"
      zRecWeight[0] = 0;
    }
    // MC signal reco
    if(opt.First('r') != kNPOS)
    {
      zHRec[i]->Reset(); // "r"
      zRecWeight[1] = 0;
    }
    // MC signal gen
    if(opt.First('g') != kNPOS)
    {
      zHGen[i]->Reset(); // "g"
    }
    // MC non-ttbar background
    if(opt.First('b') != kNPOS)
    {
      zHBgr[i]->Reset(); // "b"
      for(unsigned int a = 3; a < zRecWeight.size(); a++)
        zRecWeight[a] = 0;
    }
    // MC ttbar background
    if(opt.First('o') != kNPOS)
    {
      zHBtt[i]->Reset(); // "o"
      zRecWeight[2] = 0;
    }
    // MC ttbar out of space background
    if(opt.First('o') != kNPOS && opt.First('p') != kNPOS && opt.First('s') != kNPOS)
    {
      zHOutOfSpace[i]->Reset(); // "ops"
      zRecWeight[10] = 0;
    }
  }

  // TH2 response matrix
  // "rg"
  if(opt.First('r') != kNPOS && opt.First('g') != kNPOS)
  {
    for(unsigned int i = 0; i < zHRsp.size(); i++)
      if(zHRsp[i])
        zHRsp[i]->Reset();
  }
}

int UnfoldingIOHandler::NBinsXSec(const int binning) const
{
  TH1D* h = (TH1D*)zTUnfBinning[binning]->CreateHistogram("");
  int n = h->GetNbinsX();
  delete h;
  return n;
}

double UnfoldingIOHandler::Reweight(const TreePlain* tree)
{
  // 01.05.18 now this is fake
  //return 1.0;

  double wTotal = 1.0;
  if(zvVFRew)
  {
    double w1 = 1.0;
    for(auto& it : *zvVFRew)
    {
      double v1 = it.first->Get(tree);
      double v[1] = {v1};
      double w = it.second->EvalPar(v);
      w1 *= w;
    }
    wTotal *= w1;
  }
  if(zFRew)
  {
    double w2 = 1.0;
    if(Dim() == 1)
    {
      double v1 = zVVar[0]->Get(tree);
      double v[1] = {v1};
      w2 = zFRew->EvalPar(v);
      if(w2 != w2 || std::isinf(w2))
      {
        printf("reweight(%f) = %f, reset to 1.0\n", v1, zFRew->EvalPar(v));
        w2 = 1.0;
      }
    }
    else if(Dim() == 2)
    {
      double v1 = zVVar[0]->Get(tree);
      double v2 = zVVar[1]->Get(tree);
      double v[2] = {v1, v2};
      //zFRew->Print("all");
      w2 = zFRew->EvalPar(v);
      //printf("reweight(%f %f) = %f, reset to 1.0\n", v1, v2, zFRew->EvalPar(v));
      if(w2 != w2 || std::isinf(w2))
      {
        //printf("reweight(%f %f) = %f, reset to 1.0\n", v1, v2, zFRew->EvalPar(v));
        w2 = 1.0;
      }
    }
    else if(Dim() == 3)
    {
      double v1 = zVVar[0]->Get(tree);
      double v2 = zVVar[1]->Get(tree);
      double v3 = zVVar[2]->Get(tree);
      double v[3] = {v1, v2, v3};
      w2 = zFRew->EvalPar(v);
      //w2 = 1.0;
      if(w2 != w2 || std::isinf(w2))
      {
        printf("reweight(%f %f %f) = %f, reset to 1.0\n", v1, v2, v3, zFRew->EvalPar(v));
        w2 = 1.0;
      }
    }
    else
    {
      printf("Error in Reweight(): Dim > 3 not implemented\n");
      throw;
    }
    wTotal *= w2;
  }
  return wTotal;
}

// TODO check should it return int
double UnfoldingIOHandler::GetBin(const TreePlain* tree, const int binning) const
{
  int bin = -1;
  if(Dim() == 1)
  {
    double v1 = zVVar[0]->Get(tree);
    bin = zTUnfBinningNode[binning]->GetGlobalBinNumber(v1);
    //printf("%e  %d\n", v1, bin);
  }
  else if(Dim() == 2)
  {
    //double v1 = zVVar[0]->Get(tree);
    //double v2 = zVVar[1]->Get(tree);
    double v2 = zVVar[0]->Get(tree);
    double v1 = zVVar[1]->Get(tree);
    bin = zTUnfBinningNode[binning]->GetGlobalBinNumber(v1, v2);
  }
  else if(Dim() == 3)
  {
    //double v1 = zVVar[0]->Get(tree);
    //double v2 = zVVar[1]->Get(tree);
    //double v3 = zVVar[2]->Get(tree);
    double v3 = zVVar[0]->Get(tree);
    double v2 = zVVar[1]->Get(tree);
    double v1 = zVVar[2]->Get(tree);
    //printf("%f  %f  %f\n", v1, v2, v3);
    bin = zTUnfBinningNode[binning]->GetGlobalBinNumber(v1, v2, v3);
  }
  else
  {
    printf("Error in GetBin(): Dim > 3 not implemented\n");
    throw;
  }
  //if(bin == 0) printf("bin = 0  var = %f  tree = %s\n", zVVar[0]->Get(tree), tree->GetName());
  return bin;
}

int UnfoldingIOHandler::GetBin(const double* val, const int binning) const
{
  int bin = zTUnfBinningNode[binning]->GetGlobalBinNumber(val);
  return bin;
}

void UnfoldingIOHandler::FillTH1(TH1* h, const TreePlain* tree, const int binning, double w)
{
  int bin = GetBin(tree, binning);
  // ignore bins not in histo range
  if(bin < 1 || bin > h->GetNbinsX())
    return;
  h->Fill(bin, w);
}

// use weight from treeX, reweight treeY
void UnfoldingIOHandler::FillTH2(TH2* h, const TreePlain* treeX, const int binningX, const TreePlain* treeY, const int binningY, const double w)
{
  int binX = GetBin(treeX, binningX);
  if(binX < 1 || binX > h->GetNbinsX())
    return;
  int binY = GetBin(treeY, binningY);
  // ignore bins not in histo range
  if(binY < 1 || binY > h->GetNbinsY())
    return;
  //printf("binX = %d binY = %d w = %f\n", binX, binY, w);
  h->Fill(binX, binY, w);
}

void UnfoldingIOHandler::FillTH2UnderflowX(TH2* h, const TreePlain* treeY, const int binningY, const double w)
{
  int binY = GetBin(treeY, binningY);
  // ignore bins not in histo range
  if(binY < 1 || binY > h->GetNbinsY())
    return;
  h->Fill(0.0, (double)binY, w);
}

void UnfoldingIOHandler::FillData(const TreePlain* tree, const double rew)
{
  double w = tree->Weight(flagApplyCut) * rew;
  zRecWeight[0] += w;
  for(int i = 0; i < NBinnings(); i++)
    FillTH1(zHDat[i], tree, i, w);
}

void UnfoldingIOHandler::FillBgr(const TreePlain* tree, const double rew)
{
  double w = tree->Weight(flagApplyCut) * rew;
  zRecWeight[3] += w;
  zRecWeight[tree->SampleType()] += w;
  for(int i = 0; i < NBinnings(); i++)
    FillTH1(zHBgr[i], tree, i, w);
}

void UnfoldingIOHandler::FillBtt(const TreePlain* tree, const double rew)
{
  double w = tree->Weight(flagApplyCut) * rew;
  zRecWeight[2] += w;
  for(int i = 0; i < NBinnings(); i++)
    FillTH1(zHBtt[i], tree, i, w);
}

void UnfoldingIOHandler::FillOutOfSpace(const TreePlain* tree, const double rew)
{
  double w = tree->Weight(flagApplyCut) * rew;
  zRecWeight[10] += w;
  for(int i = 0; i < NBinnings(); i++)
    FillTH1(zHOutOfSpace[i], tree, i, w);
}

void UnfoldingIOHandler::FillMC(const TreePlain* treeRec, const TreePlain* treeGen, const double rew)
{
  double wRec = (treeRec ? (treeRec->Weight(flagApplyCut)) : 1.0) * rew;
  double wGen = (treeGen ? (treeGen->Weight(flagApplyCut)) : 1.0) * rew;
  if(treeRec)
    zRecWeight[1] += wRec;
  if(treeGen && treeRec)
  {
    //printf("GEN and REC\n");
    FillTH2(zHRsp[0], treeRec, 0, treeGen, 0, wRec);
    if(NBinnings() > 1)
    {
      FillTH2(zHRsp[1], treeRec, 1, treeGen, 0, wRec);
      FillTH2(zHRsp[2], treeRec, 1, treeGen, 0, wRec);
    }
    // suppress reco-gen weight in underflow bin
    FillTH2UnderflowX(zHRsp[0], treeGen, 0,  1 * wGen);
    FillTH2UnderflowX(zHRsp[0], treeGen, 0, -1 * wRec);
    if(NBinnings() > 1)
    {
      FillTH2UnderflowX(zHRsp[1], treeGen, 0,  1 * wGen);
      FillTH2UnderflowX(zHRsp[1], treeGen, 0, -1 * wRec);
    }
    for(int i = 0; i < NBinnings(); i++)
    {
      FillTH1(zHGen[i], treeGen, i, wGen);
      FillTH1(zHRec[i], treeRec, i, wRec);
    }
  }
  else if(treeGen)
  {
    //printf("GEN\n");
    FillTH2UnderflowX(zHRsp[0], treeGen, 0, wGen);
    if(NBinnings() > 1)
      FillTH2UnderflowX(zHRsp[1], treeGen, 0, wGen);
    for(int i = 0; i < NBinnings(); i++)
      FillTH1(zHGen[i], treeGen, i, wGen);
  }
  else if(treeRec)
  {
    //printf("Event rec but not gen, should not happen\n");
    //throw;
    //printf("REC\n");
    for(int i = 0; i < NBinnings(); i++)
      FillTH1(zHRec[i], treeRec, i, wRec);
  }
}

void UnfoldingIOHandler::StandaloneFillMCGen(const TreePlain* treeGen, const int binning, TH1D* h, const double rew)
{
  double w = treeGen->Weight(flagApplyCut) * rew;
  FillTH1(h, treeGen, binning, w);
}

void UnfoldingIOHandler::AddOutOfSpaceToTTother()
{

  for(int i = 0; i < NBinnings(); i++)
    zHBtt[i]->Add(zHOutOfSpace[i], 1.0);
}

void UnfoldingIOHandler::SubtrackBackground(const bool flagDo/* = true*/)
{
  if(zHDnb.size() != 0)
  {
    //printf("Warning in ZXsec::SubtrackBackground() non empty hDnb.size() = %ld, clearing\n", zHDnb.size());
    /*for(auto& h : zHDnb)
      if(h)
        delete h;*/
    zHDnb.clear();
  }

  for(int i = 0; i < NBinnings(); i++)
  {
    //zHDat[i]->Print("all");
    if(flagDo)
      zHDnb.push_back(SubtractBkg(zHDat[i], zHBgr[i], zHBtt[i], zHRec[i]));
    else
      zHDnb.push_back((TH1D*) zHDat[i]->Clone());
  }

  int b = 0;
  if(configHelper_->gTTBARDebug && 0)
  {
    printf("Total number of data events: %f with uo %f\n", zHDat[b]->Integral(), zHDat[b]->Integral(0, zHDat[b]->GetNbinsX() + 1));
    printf("Total number of data events no background : %f with uo %f\n", zHDnb[b]->Integral(), zHDnb[b]->Integral(0, zHDnb[b]->GetNbinsX() + 1));
    printf("Total number of non-background events: %f with uo %f\n", zHBgr[b]->Integral(), zHBgr[b]->Integral(0, zHBgr[b]->GetNbinsX() + 1));
    printf("Total number of tt background events: %f with uo %f\n", zHBtt[b]->Integral(), zHBtt[b]->Integral(0, zHBtt[b]->GetNbinsX() + 1));
    if (zHOutOfSpace[b])
      printf("Total number of out of space background events: %f with uo %f\n", zHOutOfSpace[b]->Integral(), zHOutOfSpace[b]->Integral(0, zHOutOfSpace[b]->GetNbinsX() + 1));
    printf("Total number of reconstructed events: %f with uo %f\n", zHRec[b]->Integral(), zHRec[b]->Integral(0, zHRec[b]->GetNbinsX() + 1));
    printf("Total number of generated events: %f with uo %f\n", zHGen[b]->Integral(), zHGen[b]->Integral(0, zHGen[b]->GetNbinsX() + 1));
  }
}

void UnfoldingIOHandler::FillLastDimHisto(TH1D* hFill, const TH1D* hGlobal, const bool flagNorm/* = 0*/) const
{
  for(int b = 0; b < hFill->GetNbinsX(); b++)
  {
    int gb = (int) (hFill->GetBinContent(b + 1) + 0.5);
    double bw = flagNorm ? (hFill->GetBinLowEdge(b + 2) - hFill->GetBinLowEdge(b + 1)) : 1.0;
    hFill->SetBinContent(b + 1, hGlobal->GetBinContent(gb) / bw);
    hFill->SetBinError(b + 1, hGlobal->GetBinError(gb) / bw);
  }
}

std::vector<TH1D*> UnfoldingIOHandler::CreateLastDimHistos(const int binning, const int dim/* = 0*/, std::vector<double> varArray/* = {}*/) const
{
  //printf("called CreateLastDimHistos binning = %d dim = %d\n", binning, dim);
  if(dim < 0 || dim >= Dim())
  {
    printf("Error in CreateLastDimHistos() invalid dim = %d for Dim() = %d\n", dim, Dim());
    throw;
  }
  std::vector<TH1D*> out;
  const ZVar* var = zVVar[dim];
  const std::vector<double>& bins = var->Bins[binning];
  if(dim == Dim() - 1) // last dimension
  {
    TH1D* h = var->CreateTH1D(binning, "");
    for(int b = 1; b <= h->GetNbinsX(); b++)
    {
      double value = bins[b - 1] + 0.5 * (bins[b] - bins[b - 1]);
      varArray.push_back(value);
      //printf("varArray size = %d: %f %f\n", varArray.size(), varArray[0], varArray[1]);
      std::vector<double>  reversedValArray = varArray;
      std::reverse(reversedValArray.begin(), reversedValArray.end());
      //int gbin = zTUnfBinningNode[binning]->GetGlobalBinNumber(&varArray[0]);
      int gbin = zTUnfBinningNode[binning]->GetGlobalBinNumber(&reversedValArray[0]);
      varArray.pop_back();
      h->SetBinContent(b, gbin);
    }
    out.push_back(h);
  }
  else // iterative call
  {
    for(unsigned int dd = 0; dd < bins.size() - 1; dd++)
    {
      double value = bins[dd] + 0.5 * (bins[dd + 1] - bins[dd]);
      varArray.push_back(value);
      std::vector<TH1D*> ret = CreateLastDimHistos(binning, dim + 1, varArray);
      varArray.pop_back();
      for(unsigned int i = 0; i < ret.size(); i++)
      {
        TString nameExt = TString::Format("d%db%d", dim, dd);
        TString name = TString::Format("%s_%s", ret[i]->GetName(), nameExt.Data());
        ret[i]->SetName(name);
        TString titleExt = var->GetIneq(binning, dd);
        TString title = ret[i]->GetTitle();
        if(title != TString(""))
          title += ", ";
        title += titleExt;
        ret[i]->SetTitle(title);
      }
      out.insert(out.end(), ret.begin(), ret.end());
    }
  }
  return out;
}

void UnfoldingIOHandler::InitTree(TreePlain* tree, bool flagNoRew)
{
  // 13.06.18 flagNoRew is used to skip initialising ptt variable in reco level tree for pT(t) reweighting
  // (no such variable in simplified kin reco)
  flagApplyCut = false;
  for(auto var : zVVar)
  {
    var->RegisterTree(tree);
    if(var->ExpressionCut != "")
      flagApplyCut = true;
  }
  if(zvVFRew && !flagNoRew)
  {
    for(auto& it : *zvVFRew)
    {
      ZVar* var = it.first;
      var->RegisterTree(tree);
    }
  }
}

//void UnfoldingIOHandler::InitTreeForReweighting(TreePlain* tree)
//{
//}

void UnfoldingIOHandler::ClearTree(TreePlain* tree, bool flagNoRew)
{
  // 13.06.18 (same as above) flagNoRew is used to skip initialising ptt variable in reco level tree for pT(t) reweighting
  // (no such variable in simplified kin reco)
  flagApplyCut = false;
  for(auto var : zVVar)
    var->DeregisterTree(tree);
  if(zvVFRew && !flagNoRew)
  {
    for(auto& it : *zvVFRew)
    {
      ZVar* var = it.first;
      var->DeregisterTree(tree);
    }
  }
}

void UnfoldingIOHandler::PutPredictionPttt(double ptThreshold)
{
  int nptttAll = ZMadGraph::GetOffsetHistoPtttThreshold(configHelper_, -1.0);
  int npttt = ZMadGraph::GetOffsetHistoPtttThreshold(configHelper_, ptThreshold);
  std::vector<int> hid[2];
  for(int p = 0; p < nptttAll; p++)
  {
    if(p < npttt)
      hid[0].push_back(100 + 4 * (p + 1));
    else
      hid[1].push_back(100 + 4 * (p + 1));
  }
  // 1st bin pttt
  for(int i = 1; i <= 4; i++)
  {
    std::vector<int> vh;
    for(size_t ii = 0; ii < hid[0].size(); ii++)
      vh.push_back(hid[0][ii] + i);
    for(size_t ii = 0; ii < hid[1].size(); ii++)
      vh.push_back(hid[1][ii] + i);
    VMGHisto.push_back(std::pair<std::vector<int>, int>(vh, 0));
    vh.clear();
    for(size_t ii = 0; ii < hid[1].size(); ii++)
      vh.push_back(hid[1][ii] + i);
    VMGHistoMinus.push_back(std::pair<std::vector<int>, int>(vh, 1));
  }
  // 2nd bin pttt
  for(int i = 1; i <= 4; i++)
  {
    std::vector<int> vh;
    for(size_t ii = 0; ii < hid[1].size(); ii++)
      vh.push_back(hid[1][ii] + i);
    VMGHisto.push_back(std::pair<std::vector<int>, int>(vh, 1));
    VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));
  }
}

void UnfoldingIOHandler::PutPredictionPtttIntegrate(double ptThreshold)
{
  int nptttAll = ZMadGraph::GetOffsetHistoPtttThreshold(configHelper_, -1.0);
  int npttt = ZMadGraph::GetOffsetHistoPtttThreshold(configHelper_, ptThreshold);
  std::vector<int> hid[2];
  for(int p = 0; p < nptttAll; p++)
  {
    if(p < npttt)
      hid[0].push_back(100 + 4 * (p + 1));
    else
      hid[1].push_back(100 + 4 * (p + 1));
  }
  // 1st bin pttt
  std::vector<int> vh;
  for(int i = 1; i <= 4; i++)
  {
    for(size_t ii = 0; ii < hid[0].size(); ii++)
      vh.push_back(hid[0][ii] + i);
    for(size_t ii = 0; ii < hid[1].size(); ii++)
      vh.push_back(hid[1][ii] + i);
  }
  VMGHisto.push_back(std::pair<std::vector<int>, int>(vh, 0));
  vh.clear();
  for(int i = 1; i <= 4; i++)
  {
    for(size_t ii = 0; ii < hid[1].size(); ii++)
      vh.push_back(hid[1][ii] + i);
  }
  VMGHistoMinus.push_back(std::pair<std::vector<int>, int>(vh, 1));
  // 2nd bin pttt
  vh.clear();
  for(int i = 1; i <= 4; i++)
  {
    for(size_t ii = 0; ii < hid[1].size(); ii++)
      vh.push_back(hid[1][ii] + i);
  }
  VMGHisto.push_back(std::pair<std::vector<int>, int>(vh, 1));
  vh.clear();
  VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));
}

std::pair<std::vector<double>,std::vector<double>> UnfoldingIOHandler::GetBinningInfoFromFile(const TString binningFileName){
    std::vector<double> rec_lvl_bins = {};
    std::vector<double> det_lvl_bins = {};
    TString label_for_rec_bins = "rec_level";
    TString label_for_det_bins = "det_level";

    std::ifstream input_file;
    input_file.open(binningFileName);
    std::string temp_line = "";

    if (input_file.is_open()){
        while ( getline (input_file,temp_line) ){
            TString tstring_line = temp_line;
            //skipping comments
            if (tstring_line.BeginsWith("#")) continue;
            // std::cout << tstring_line << std::endl;

            //Spliting label and data
            TObjArray* label_and_info  = tstring_line.Tokenize(":");
            TString data_label = (((TObjString*)label_and_info->At(0))->GetString()).Data();
            TString data_in_tstring_format = (((TObjString*)label_and_info->At(1))->GetString()).Data();

            //Checking for rec. lvl. binning
            if (data_label == label_for_rec_bins) {
                rec_lvl_bins = utils::GetVectorFromTStringData(data_in_tstring_format);
            }
            //Checking for det. lvl. binning
            else if ( data_label == label_for_det_bins){
                det_lvl_bins = utils::GetVectorFromTStringData(data_in_tstring_format);
            }
        }
        input_file.close();
    }
    else{
        std::cerr << "ERROR in UnfoldingIOHandler::GetBinningInfoFromFile -->> File not found: " << binningFileName << std::endl;
        exit(1);
    }

    return std::pair<std::vector<double>,std::vector<double>>(rec_lvl_bins,det_lvl_bins);
}
