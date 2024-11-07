#include <fstream>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <sstream>
#include <cmath>
#include <memory>
#include <cassert>

#include <TTree.h>
#include <TLegend.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TPaveText.h>

#include "PlotterDiffXSec.h"

#include "PlotterConfigurationHelper.h"
#include "utils.h"
#include "ttmd_declareVars.h"
#include "UnfoldingXSecHandler.h"
#include "ttmd_unfold.h"
#include "UnfoldingClosureTestsHandler.h"
#include "UnfoldingIOHandler.h"
#include "ttmd_measur.h"
#include "ttmd_ZMadGraph.h"

using namespace std;

int PlotterDiffXSec::NDigFloatXSec = 5;
int PlotterDiffXSec::NDigFloatStat = 4;
int PlotterDiffXSec::NDigFloatSyst = 3;
int PlotterDiffXSec::NDigFloatVar  = 2;
int PlotterDiffXSec::NDigFloatCorr = 3;

PlotterDiffXSec::PlotterDiffXSec(TString analysis_year,PlotterConfigurationHelper *configHelper_ , int mode, const TString& analysisSuffix, int ModeKinRecoTree) :
configHelper(configHelper_),
zMode(mode),
AnalysisSuffix(analysisSuffix),
ModeKinRecoTree(ModeKinRecoTree)
{
  year = analysis_year;
  current_channel = "";

  MaxEvents = 0;
  zVarLast = "Nominal";
  //DirMCSyst = "Nominal";
  //DirMCBkgSyst = "Nominal";

  FlagSkipExtraPlots = configHelper->ttmd_skip_extra_plots_for_systematics;

  //if(mode == 0)
  //  {
  //    configHelper->gTreeDir  = configHelper->gBaseDir + "../plainTree-fr/" + configHelper->gPlainTreeVersion;
  //    configHelper->gAnalDir = configHelper->gAnalDir + "-fr";
  //  }
  if(mode == 1)
    {
      std::cout << "Running on analysis mode 1" << std::endl;
    }
  else
    printf("Warning: unknown analysis mode %d\n", zMode);

  if(AnalysisSuffix != "")
    configHelper->gAnalDir = configHelper->gAnalDir + "/" + AnalysisSuffix;
  else
    configHelper->gAnalDir = configHelper->gAnalDir + "/def";

  if(ModeKinRecoTree == 0)
    configHelper->gPlainTreeNamePlain = "plainTree_rec_step8";
  else if(ModeKinRecoTree == 9)
    configHelper->gPlainTreeNamePlain = "plainTree_recSimple_step7L";
    //configHelper->gPlainTreeNamePlain = "plainTree_recSimple_step8"; //old name for loose kin reco plain trees to run on older plain tree versions
  else
    throw std::logic_error(TString::Format("Error: unsupported ModeKinRecoTree = %d", ModeKinRecoTree).Data());
}

PlotterDiffXSec::PlotterDiffXSec():
zMode(0),
AnalysisSuffix(""),
ModeKinRecoTree(0)
{}

int PlotterDiffXSec::FlagPaper = 2;
bool PlotterDiffXSec::FlagAnalysisNote = false;

std::vector<std::vector<TH1D*> > PlotterDiffXSec::CreateCPHistos(const UnfoldingIOHandler* unfoldingIOHandler, const int binning, const std::vector<TH1D*> vh)
{
  std::vector<std::vector<TH1D*> > ret;
  std::vector<TH1D*> hh = unfoldingIOHandler->CreateLastDimHistos(binning);
  for(unsigned int i = 0; i < hh.size(); i++)
  {
    ret.push_back(std::vector<TH1D*> {});

    // non-ttbar
    ret[i].push_back((TH1D*)hh[i]->Clone());
    ret[i].back()->SetTitle("Non-t#bar{t}");
    ret[i].back()->SetFillColor(4);
    ret[i].back()->SetLineColor(1);
    unfoldingIOHandler->FillLastDimHisto(ret[i].back(), vh[3]);
    
    // ttbar other
    ret[i].push_back((TH1D*)hh[i]->Clone());
    ret[i].back()->SetTitle("t#bar{t} other");
    //ret[i].back()->SetFillColor(6);
    ret[i].back()->SetFillColor(kMagenta - 7);
    ret[i].back()->SetLineColor(1);
    unfoldingIOHandler->FillLastDimHisto(ret[i].back(), vh[2]);
    ret[i].back()->Add(ret[i][ret[i].size() - 2]);
    
    // ttbar signal
    ret[i].push_back((TH1D*)hh[i]->Clone());
    ret[i].back()->SetTitle("t#bar{t} signal");
    ret[i].back()->SetFillColor(2);
    ret[i].back()->SetLineColor(1);
    unfoldingIOHandler->FillLastDimHisto(ret[i].back(), vh[1]);
    ret[i].back()->Add(ret[i][ret[i].size() - 2]);

    // data
    ret[i].push_back((TH1D*)hh[i]->Clone());
    //ret[i].back()->SetTitle("Data");
    ret[i].back()->SetMarkerStyle(20);
    ret[i].back()->SetMarkerColor(1);
    ret[i].back()->SetMarkerSize(0.75);
    ret[i].back()->SetLineColor(1);
    unfoldingIOHandler->FillLastDimHisto(ret[i].back(), vh[0]);
  }
  
  for(unsigned int i = 0; i < hh.size(); i++)
    delete hh[i];
  
  return ret;
}

void PlotterDiffXSec::SetYMinMax(std::vector<std::vector<TH1D*> >& hh, Params* pars)
{
  if(pars->YMax < 0.0)
  {
    for(unsigned int i = 0; i < hh.size(); i++)
      for(unsigned int j = 0; j < hh[i].size(); j++)
        if(pars->YMax < hh[i][j]->GetMaximum())
          pars->YMax = hh[i][j]->GetMaximum();
    pars->YMax *= (pars->LogY) ? 2.5 : 1.3;
  }
  if(pars->YMin < 0.0)
  {
    if(pars->LogY)
    {
      pars->YMin = pars->YMax;
      for(unsigned int i = 0; i < hh.size(); i++)
        for(unsigned int j = 0; j < hh[i].size(); j++)
          if(pars->YMin > hh[i][j]->GetMinimum())
            pars->YMin = hh[i][j]->GetMinimum();
      pars->YMin /= 2.0;
    }
    else
      pars->YMin = 0.0;
  }
}

std::vector<std::pair<TH1D*, TH1D*> > PlotterDiffXSec::CreateLastDimUDHistograms(const UnfoldingIOHandler* unfoldingIOHandler, const std::pair<TH1D*, TH1D*>& pairUD, const int binning, const int flagNorm)
{
  std::vector<std::pair<TH1D*, TH1D*> > result;
  auto vh = unfoldingIOHandler->CreateLastDimHistos(binning);
  for(unsigned int i = 0; i < vh.size(); i++)
  {
    auto h1 = (TH1D*) vh[i]->Clone();
    auto h2 = (TH1D*) vh[i]->Clone();
    result.push_back(std::make_pair(h1, h2));
    unfoldingIOHandler->FillLastDimHisto(result[i].first, pairUD.first, flagNorm);
    unfoldingIOHandler->FillLastDimHisto(result[i].second, pairUD.second, flagNorm);
  }
  return result;
}

std::vector<std::vector<std::pair<TH1D*, TH1D*> > > PlotterDiffXSec::CreateLastDimUDHistograms(const UnfoldingIOHandler* unfoldingIOHandler, const std::vector<std::pair<TH1D*, TH1D*> >& pairUD, const int binning, const int flagNorm)
{
  std::vector<std::vector<std::pair<TH1D*, TH1D*> > > result;
  auto vh = unfoldingIOHandler->CreateLastDimHistos(binning);
  for(unsigned int i = 0; i < vh.size(); i++)
  {
    result.push_back(std::vector<std::pair<TH1D*, TH1D*> >());
    for(unsigned int j = 0; j < pairUD.size(); j++)
    {
      auto h1 = (TH1D*) vh[i]->Clone();
      auto h2 = (TH1D*) vh[i]->Clone();
      result[i].push_back(std::make_pair(h1, h2));
      unfoldingIOHandler->FillLastDimHisto(result[i][j].first, pairUD[j].first, flagNorm);
      unfoldingIOHandler->FillLastDimHisto(result[i][j].second, pairUD[j].second, flagNorm);
    }
  }
  return result;
}

void PlotterDiffXSec::DrawBinLabel(const UnfoldingIOHandler* unfoldingIOHandler, const int i, const double axisX, const double histoX)
{
  return;
  {
    //TPaveText* pt = new TPaveText(0.05 + ((i == 0) ? 0.32 * ((axisX + histoX) / histoX - 1.0) : 0.0), 0.80, 0.95, 0.90, "ndc");
    TPaveText* pt = new TPaveText(0.05 + ((i == 0) ? 0.37 * ((axisX + histoX) / histoX - 1.0) : 0.0), 0.80, 0.95, 0.90, "ndc");
    pt->SetTextSize(25);
    if(unfoldingIOHandler->Dim() == 3)
      pt->SetTextSize(15);
    pt->SetTextFont(63);
    pt->SetBorderSize(0);
    int divide = 1;
    for(int d = 1; d < unfoldingIOHandler->Dim(); d++)
    {
      int bin = (i / divide) % (unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - d - 1)->BinsC().size() - 1);
      divide = unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - d - 1)->BinsC().size() - 1;
      bin += 1;
      pt->AddText(TString::Format("%s #%d", unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - d - 1)->Title.Data(), bin));
    }
    pt->SetFillStyle(0);
    pt->SetFillColor(0);
    pt->Draw();
  }
}

void PlotterDiffXSec::DrawBinLabelFullVertical(const UnfoldingIOHandler* unfoldingIOHandler, const int i/*, const double axisX, const double histoX*/)
{
  bool flagLargest = (unfoldingIOHandler->Dim() == 3) && (unfoldingIOHandler->Var(0)->Bins[0].size() >= 4);
  //TPaveText* pt = new TPaveText(0.05 + ((i == 0) ? 0.32 * ((axisX + histoX) / histoX - 1.0) : 0.0), 0.80, 0.95, 0.90, "ndc");
  //TPaveText* pt = new TPaveText(0.05 + ((i == 0) ? (flagLargest ? 0.38 : 0.47) * ((axisX + histoX) / histoX - 1.0) : 0.0), 0.80, 0.98, 0.90, "ndc");
  TPaveText* pt = new TPaveText(0.52, 0.40, 0.98, 0.90, "ndc");
  //pt->SetTextSize(17);
  pt->SetTextSize(24);
  //if(flagLargest)
  //  pt->SetTextFont(43);
  //else
    pt->SetTextFont(63);
  pt->SetBorderSize(0);
  if(unfoldingIOHandler->Dim() == 3)
  {
    //pt->SetTextSize(12);
    pt->SetTextSize(20);
    if(flagLargest)
    {
      //pt->SetTextFont(83);
      //pt->SetTextSize(8);
      pt->SetTextSize(20);
    }
  }
  //pt->SetTextFont(73);
  //pt->SetTextSize(0.04);
  int divide = 1;
  //pt->SetTextAngle(90);
  TString strFull;
  for(int d = 1; d < unfoldingIOHandler->Dim(); d++)
  {
    int bin = (i / divide) % (unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - d - 1)->BinsC().size() - 1);
    divide = unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - d - 1)->BinsC().size() - 1;
    bin += 1;
    //pt->AddText(TString::Format("%s #%d", unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - d - 1)->Title.Data(), bin));
    TString str = unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - d - 1)->GetIneq(0, bin - 1);
    // some customising
    if(str == "-0.5 < N_{jet} < 0.5") str = "N_{jet} = 0";
    else if(str == "0.5 < N_{jet} < 1.5") str = "N_{jet} = 1";
    else if(str == "0.5 < N_{jet} < 8.5") str = "N_{jet} > 0";
    else if(str == "1.5 < N_{jet} < 8.5") str = "N_{jet} > 1";
    else if(str == "1.5 < N_{jet} < 2.5") str = "N_{jet} = 2";
    else if(str == "2.5 < N_{jet} < 8.5") str = "N_{jet} > 2";

    str.ReplaceAll(" ", "");
    strFull += str;
    if(d == 1)
      strFull += ",  ";
    //printf("str = %s\n", str.Data());
    //if(d == 1)
    //  tt->SetTextAngle(90);
  }
  auto tt = pt->AddText(strFull);
  tt->SetTextAngle(90);
  pt->SetFillStyle(0);
  pt->SetFillColor(0);
  //pt->SetTextAngle(90);
  pt->Draw();
}

void PlotterDiffXSec::DrawBinLabelFull(const UnfoldingIOHandler* unfoldingIOHandler, const int i, const double axisX, const double histoX, TString plot_type)
{
  if (plot_type != "xsec" && plot_type != "cp") {
      std::cerr << "ERROR in PlotterDiffXSec::DrawBinLabelFull -->> plot_type not implemented: " << plot_type << std::endl;
      exit(1);
  }
  // 14.11.18 split long line
  bool flagSplitLine = 1;
  // 7.11.18
  //DrawBinLabelFullVertical(unfoldingIOHandler, i, axisX, histoX);
  //return;
  //
  bool flagLargest = (unfoldingIOHandler->Dim() == 3) && (unfoldingIOHandler->Var(0)->Bins[0].size() >= 4);
  bool flag4njbins = (unfoldingIOHandler->Dim() == 3) && (unfoldingIOHandler->Var(0)->Bins[0].size() == 5);
  int totalNumberOfLabeledBins = 1;
  for (int i=0; i<unfoldingIOHandler->Dim()-1; i++){
    totalNumberOfLabeledBins = totalNumberOfLabeledBins * (unfoldingIOHandler->Var(i)->Bins[0].size()-1);
  }
  //TPaveText* pt = new TPaveText(0.05 + ((i == 0) ? 0.32 * ((axisX + histoX) / histoX - 1.0) : 0.0), 0.80, 0.95, 0.90, "ndc");
  //TPaveText* pt = new TPaveText(0.05 + ((i == 0) ? (flagLargest ? 0.38 : 0.47) * ((axisX + histoX) / histoX - 1.0) : 0.0), 0.80, 0.98, 0.90, "ndc");
  TPaveText* pt = new TPaveText(0.05 + ((i == 0) ? (flagLargest ? 0.38 : 0.47) * ((axisX + histoX) / histoX - 1.0) + 0.10 : -0.02), 0.80, 0.98, 0.90, "ndc");
  //pt->SetBorderSize(2);
  if(flagSplitLine)
  {
    pt->SetY1(0.69);
    pt->SetY2(0.89);
    if(unfoldingIOHandler->Dim() == 3)
    {
      pt->SetY1(0.67);
      pt->SetY2(0.87);
      if(unfoldingIOHandler->Var(0)->Bins[0].size() != 4) // only for nj2mttytt
      {
        pt->SetX1(pt->GetX1() - 0.05);
        pt->SetX2(pt->GetX2() - 0.05);
      }
    }
  }
  if(unfoldingIOHandler->Dim() == 2)
  {
    pt->SetY2(pt->GetY2() + 0.05);
    pt->SetY1(pt->GetY1() + 0.05);
  }

  if (plot_type == "xsec"){
    pt->SetTextSize(28);
  }
  else if (plot_type == "cp"){
    pt->SetTextSize(18);
  }
  else {
     pt->SetTextSize(28);
  }

  //if(flagLargest && !flagSplitLine)
  //  pt->SetTextFont(43);
  //else
  //  pt->SetTextFont(63);
  pt->SetTextFont(43);
  pt->SetBorderSize(0);
  if(unfoldingIOHandler->Dim() == 3)
  {
    //pt->SetTextSize(12);
    pt->SetTextSize(20);
    if (totalNumberOfLabeledBins>=12) pt->SetTextSize(9);
    else if(flagLargest)
    {
      //pt->SetTextFont(83);
      //pt->SetTextSize(8);
      if (flag4njbins) {
          if (plot_type == "xsec") pt->SetTextSize(16);
          else pt->SetTextSize(9);
      }
      else {
          if (plot_type == "xsec") pt->SetTextSize(14);
          else pt->SetTextSize(10);
      }
    }
  }
  //pt->SetTextFont(73);
  //pt->SetTextSize(0.04);
  int divide = 1;
  for(int d = 1; d < unfoldingIOHandler->Dim(); d++)
  {
    int bin = (i / divide) % (unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - d - 1)->BinsC().size() - 1);
    divide = unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - d - 1)->BinsC().size() - 1;
    bin += 1;
    //pt->AddText(TString::Format("%s #%d", unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - d - 1)->Title.Data(), bin));
    TString str = unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - d - 1)->GetIneq(0, bin - 1);
    // some customising
    // some customising
    if(str == "-0.5 < N_{jet} < 0.5") str = "N_{jet} = 0";
    else if(str == "0.5 < N_{jet} < 1.5") str = "N_{jet} = 1";
    else if(str == "0.5 < N_{jet} < 8.5") str = "N_{jet} > 0";
    else if(str == "1.5 < N_{jet} < 8.5") str = "N_{jet} > 1";
    else if(str == "1.5 < N_{jet} < 2.5") str = "N_{jet} = 2";
    else if(str == "2.5 < N_{jet} < 8.5") str = "N_{jet} > 2";

    str.ReplaceAll(" ", "");
    // 15.01.19 fix y(t) -> |y(t)|
    //if(str.Contains("y(t)"))
    //  str.ReplaceAll("y(t)", "|y(t)|");
    // 14.11.18
    if(unfoldingIOHandler->Dim() == 3 && flagSplitLine && ((d == 1) || totalNumberOfLabeledBins>=12))
    {
      int n_comp_operator = TString(str.Data()).CountChar('<');
      if (n_comp_operator == 0 || str.Contains("N_{jet}")) {
        pt->AddText(str);
      }
      else{
        if (n_comp_operator >= 1) {
        TString str1 = TString(str.Data(), str.Last('<'));
        pt->AddText(str1);
        }
        if (n_comp_operator >= 2) {
          TString str2 = TString(str.Data() + str.Last('<'));
          pt->AddText(str2);
        }
      }
      
      //pt->SetTextAlign(31);
      if (totalNumberOfLabeledBins>=12) {
         if (plot_type == "xsec") pt->SetTextSize(16);
         else pt->SetTextSize(9);
      }
      else if (flag4njbins) {
          if (plot_type == "xsec") pt->SetTextSize(16);
          else pt->SetTextSize(9);
      }
      else pt->SetTextSize(16);
    }
    else
      pt->AddText(str);
  }
  if(unfoldingIOHandler->Dim() == 3 && plot_type == "xsec") pt->SetTextAlign(31);
  else if (unfoldingIOHandler->Dim() == 3 && plot_type == "cp") pt->SetTextAlign(11);
  pt->SetFillStyle(0);
  pt->SetFillColor(0);
  pt->Draw();
}

void PlotterDiffXSec::DrawCMSTitle(int flag, int dim, TString plot_type) // 1 paper, 2 pas or work in progress
{
  TString cms_suffix_label = "";
  if (flag == 1 ) cms_suffix_label = "";
  // else if (flag == 2) cms_suffix_label = "Preliminary";
  else if (flag == 2) cms_suffix_label = "Work in progress";
  
  if (plot_type != "xsec" && plot_type != "cp" && plot_type != "pse") {
    std::cerr << "ERROR in PlotterDiffXSec::DrawCMSTitle -->> plot_type not implemented: " << plot_type << std::endl;
    exit(1);
  }

  if(flag == 3) return;

  TPaveText* pt = new TPaveText(dim == 2 ? 0.11 : 0.09, 0.96, dim == 2 ? 0.18 : 0.145, 0.99, "ndc");
  pt->SetTextSize(40);
  
  if(dim == 1 && plot_type != "cp") {
    if (FlagAnalysisNote) {
      pt = new TPaveText(0.20, 0.96, 0.27, 0.98, "ndc");
      pt->SetTextSize(35);
    }
    else {
      pt = new TPaveText(0.11, 0.96, 0.18, 0.99, "ndc");
      pt->SetTextSize(37);
    }
  }
  else if (dim == 1 && plot_type == "cp") {
    if (FlagAnalysisNote) {
      pt = new TPaveText(0.21, 0.96, 0.28, 0.98, "ndc");
      pt->SetTextSize(35);
    }
    else {
      pt = new TPaveText(0.18, 0.95, 0.25, 0.99, "ndc");
      pt->SetTextSize(37);
    }
  }
  else if (dim > 1 && plot_type == "cp") {
    pt = new TPaveText(0.18, 0.94, 0.23, 0.98, "ndc");
    pt->SetTextSize(33);
  }
  else if (plot_type == "pse") {
    pt = new TPaveText(0.21, 0.96, 0.28, 0.98, "ndc");
    pt->SetTextSize(35);
  }

  pt->SetTextFont(63);
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->AddText("CMS");
  pt->Draw();

  // add preliminary if needed
  pt = new TPaveText(dim == 2 || dim == 3 ? 0.24 : 0.185, 0.95, dim == 2 || dim == 3 ? 0.285 : 0.255, 0.99, "ndc");
  pt->SetTextSize(30);
  if(dim == 1 && plot_type != "cp") {
    if (FlagAnalysisNote) {
      pt = new TPaveText(0.34, 0.95, 0.48, 0.98, "ndc");
      pt->SetTextSize(27);
    }
    else {
      pt = new TPaveText(0.24, 0.95, 0.38, 0.99, "ndc");
      pt->SetTextSize(30);
    }
  }
  else if (dim == 1 && plot_type == "cp") {
    if (FlagAnalysisNote) {
      pt = new TPaveText(0.35, 0.95, 0.49, 0.98, "ndc");
      pt->SetTextSize(27);
    }
    else {
      pt = new TPaveText(0.18, 0.94, 0.57, 0.99, "ndc");
      pt->SetTextSize(30);
    }
  }
  else if(dim > 1 && plot_type == "cp") {
    pt = new TPaveText(0.2, 0.93, 0.53, 0.98, "ndc");
    pt->SetTextSize(33);
  }
  else if (plot_type == "pse") {
    pt = new TPaveText(0.35, 0.95, 0.49, 0.98, "ndc");
    pt->SetTextSize(27);
  }

  pt->SetTextFont(53);
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->AddText(cms_suffix_label);
  pt->Draw();

  // lumi and energy
  pt = new TPaveText(0.69, 0.95, 0.81, 0.99, "ndc");
  pt->SetTextSize(28);
  if(dim == 1 && plot_type != "cp") {
    if (FlagAnalysisNote) {
      pt = new TPaveText(0.80, 0.95, 0.88, 0.99, "ndc");
      pt->SetTextSize(33);
    }
    else {
      pt = new TPaveText(0.61, 0.95, 0.69, 0.99, "ndc");
      pt->SetTextSize(35);
    }
  }
  else if (dim == 1 && plot_type == "cp") {
    if (FlagAnalysisNote) {
      pt = new TPaveText(0.80, 0.95, 0.88, 0.99, "ndc");
      pt->SetTextSize(33);
    }
    else {
      pt = new TPaveText(0.61, 0.94, 0.71, 0.99, "ndc");
      pt->SetTextSize(35);
    }
  }
  else if(dim > 1 && plot_type == "cp") {
    pt = new TPaveText(0.61, 0.93, 0.71, 0.98, "ndc");
    pt->SetTextSize(32);
  }
  else if (plot_type == "pse") {
    pt = new TPaveText(0.80, 0.95, 0.88, 0.99, "ndc");
    pt->SetTextSize(33);
  }
  pt->SetTextFont(43);
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);

  TString label_for_lumi_and_energy = this->configHelper->lumi_fb.str() +  " fb^{-1} (" + this->configHelper->energy_label + ")";
  pt->AddText(label_for_lumi_and_energy);
  pt->Draw();
}

void PlotterDiffXSec::DrawChannelLabel(TString channel, int flag, int dim, TString plot_type, double x0Leg, double x1Leg, double y1Leg){ // 1 paper, 2 pas or work in progress

  if (plot_type != "xsec" && plot_type != "cp") {
    std::cerr << "ERROR in PlotterDiffXSec::DrawChannelLabel -->> plot_type not implemented: " << plot_type << std::endl;
    exit(1);
  }
  if (flag > 2) {
    std::cerr << "ERROR in PlotterDiffXSec::DrawChannelLabel -->> plot type flag not implemented: " << flag << std::endl;
    exit(1);
  }

  TString channel_label = "";
  if (channel == "ll") channel_label = "dilepton";
  else channel_label = channel;

  TPaveText* pt = new TPaveText();
  if (dim > 1){
    if (plot_type == "xsec") pt = new TPaveText(x0Leg * 1.01 - 0.01 - 0.07, y1Leg * 0.99 - 0.04, x1Leg * 0.99 + 0.01, y1Leg * 0.99, "ndc");
    else if (plot_type == "cp") pt = new TPaveText(x0Leg * 1.01, y1Leg * 0.99 - 0.02, x1Leg * 0.99 + 0.01, y1Leg * 1.05, "ndc");
  }
  else if (dim == 1){
    if (FlagAnalysisNote) {
      if (plot_type == "xsec") pt = new TPaveText(0.15, y1Leg * 0.98 - 0.05, 0.15+0.25, y1Leg * 1.05 - 0.05, "ndc");
      else if (plot_type == "cp") pt = new TPaveText(0.15, y1Leg * 0.98 - 0.05, 0.15+0.25, y1Leg * 1.05 - 0.05, "ndc");
    }
    else {
      if (plot_type == "xsec") pt = new TPaveText(x0Leg * 1.01, y1Leg * 0.99 - 0.02, x1Leg * 0.99 + 0.01, y1Leg * 1.05, "ndc");
      else if (plot_type == "cp") pt = new TPaveText(x0Leg * 1.01, y1Leg * 0.99 - 0.02, x1Leg * 0.99 + 0.01, y1Leg * 1.05, "ndc");
    }
  }
  else pt = new TPaveText(x0Leg * 1.01 - 0.01 - 0.07, y1Leg * 0.99 - 0.04, x1Leg * 0.99 + 0.01, y1Leg * 0.99, "ndc");

  pt->SetTextSize(33);
  pt->SetTextFont(63);
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->AddText(channel_label);
  pt->Draw();
}

std::vector<TPolyLine*> PlotterDiffXSec::CreatePolylineBins(const UnfoldingIOHandler* unfoldingIOHandler, const int nbins, const double minY, const double maxY)
{
  std::vector<TPolyLine*> plineBins;
  for(int b = 0; b < nbins; b++)
  {
    plineBins.push_back(new TPolyLine(2, &(std::vector<double>{0.5 + b, 0.5 + b}[0]), &(std::vector<double>{minY, maxY}[0])));
    plineBins[b]->SetLineStyle(3);
    plineBins[b]->SetLineWidth(1);
    int nbinsLastDim = unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - 1)->BinsC().size() - 1;
    if(b != 0 && b % nbinsLastDim == 0)
    {
      plineBins[b]->SetLineWidth(2);
      plineBins[b]->SetLineStyle(2);
      // for dim = 3
      int nbinsLastDim2 = unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - 2)->BinsC().size() - 1;
      if(unfoldingIOHandler->Dim() == 3 && ((b / nbinsLastDim) % nbinsLastDim2) == 0)
      {
        plineBins[b]->SetLineWidth(3);
        plineBins[b]->SetLineStyle(9);
      }
    }
  }
  return plineBins;
}

void PlotterDiffXSec::DrawCommonXTitle(const double x0, const double y0, const double x1, const double y1, const TH1* h)
{
  //printf("x0, y0, x1, y1: %f %f %f %f\n", x0, y0, x1, y1);
  TPaveText* ptXTitle = new TPaveText(x0, y0, x1, y1, "ndc");
  ptXTitle->SetBorderSize(0);
  //ptXTitle->SetBorderSize(1);
  //ptXTitle->SetLineColor(2);
  ptXTitle->SetTextAlign(31);
  ptXTitle->SetMargin(0.0);
  TString title = h->GetXaxis()->GetTitle();
  //if(title == "y(t#bar{t})")
  //  title = "|y(t#bar{t})|";
  //if(title == "y(t)")
  //  title = "|y(t)|";
  ptXTitle->AddText(title);
  ptXTitle->SetTextFont(h->GetXaxis()->GetTitleFont());
  ptXTitle->SetTextSize(h->GetXaxis()->GetTitleSize());
  //ptXTitle->AddText("dupa");
  ptXTitle->SetFillStyle(0);
  ptXTitle->SetFillColor(0);
  ptXTitle->Draw();
}

void PlotterDiffXSec::ControlPlot(const UnfoldingIOHandler* unfoldingIOHandler, const TString& fileName, CPParams pars)
{
  //const auto& vh = { unfoldingIOHandler->HDat(pars.Binning), unfoldingIOHandler->HRec(pars.Binning), unfoldingIOHandler->HBtt(pars.Binning), unfoldingIOHandler->HBgr(pars.Binning) };
  std::vector<TH1D*> vh;
  vh.push_back((TH1D*) unfoldingIOHandler->HDat(pars.Binning)->Clone());
  vh.push_back((TH1D*) unfoldingIOHandler->HRec(pars.Binning)->Clone());
  vh.push_back((TH1D*) unfoldingIOHandler->HBtt(pars.Binning)->Clone());
  vh.push_back((TH1D*) unfoldingIOHandler->HBgr(pars.Binning)->Clone());
  ControlPlot(unfoldingIOHandler, vh, fileName, pars);
}

void PlotterDiffXSec::ControlPlot(const UnfoldingIOHandler* unfoldingIOHandler, const std::vector<TH1D*> vh, const std::pair<TH1D*, TH1D*> vMCTotUD, const TString& fileName, CPParams pars)
{

  double marL = 0.02;
  double marR = 0.02;
  double marB = 0.02;
  double marT = 0.02;
  double axisX = 0.10;
  double axisY = 0.10;
  double ratNorm = 0.3;
  double ratLegend = 0.20;
  double markerSize = 1.5;
  if(FlagPaper && !FlagAnalysisNote)
  {
    axisX = 0.17;
    axisY = 0.14;
    marL = 0.002;
    marR = 0.002;
    marB = 0.002;
    marT = 0.02;
    ratLegend = 0.21;
  }
  else if(FlagAnalysisNote && unfoldingIOHandler->Dim() == 1)
  {
    marT = 0.007;
    ratLegend = 0.15;
  }

  double scaleAxisFontSize = 1.25;
  if(FlagPaper)
    scaleAxisFontSize = 1.6;
  double scaleOffsetsY = 1.0;

  if(unfoldingIOHandler->Dim() == 2)
  {
    scaleOffsetsY = 3.0;
  }
  if(unfoldingIOHandler->Dim() == 3)
  {
    //scaleAxisFontSize = 1.0;
    scaleOffsetsY = 5.0;
    if(FlagPaper)
      scaleOffsetsY = 4.0;
  }
  if(unfoldingIOHandler->Dim() == 1 && FlagPaper && !FlagAnalysisNote)
    scaleOffsetsY = 0.94;
  else if(FlagAnalysisNote && unfoldingIOHandler->Dim() == 1)
    scaleOffsetsY = 1.5;

  utils::ttmd::SetOZStyle();
  TCanvas* c;
  if (unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote)
    {
      c = new TCanvas("cp", "", 800, PlotPixelY);
    }
  else
    {
     c = new TCanvas("cp", "", PlotPixelX, PlotPixelY);
    }
  std::vector<std::vector<TH1D*> > hh = CreateCPHistos(unfoldingIOHandler, pars.Binning, vh);
  int nhh = hh.size();

  // create u/d histograms
  std::vector<std::pair<TH1D*, TH1D*> > hhMCTotUD;
  bool flagDrawBand = false;
  if(vMCTotUD.first && vMCTotUD.second)
  {
    flagDrawBand = true;
    hhMCTotUD = CreateLastDimUDHistograms(unfoldingIOHandler, vMCTotUD, pars.Binning, 0);
  }

  const ZVar* var = unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - 1);
  double histoX;
  if(unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote)
    histoX = (1.0 - marL - marR - axisX) / nhh;
  else
    histoX = (1.0 - marL - marR - axisX - ratLegend) / nhh;
  double histoY = (1.0 - marB - marT - axisY) / (1.0 + ratNorm);
  double histoYNorm = histoY * ratNorm;
  double x0 = marL;
  double x1 = x0 + histoX + axisX;
  double y0Norm = marB;
  double y1Norm = marB + axisY + histoYNorm;
  double y0 = y1Norm;
  double y1 = y0 + histoY;
  double epsAxis = 1e-4;

  double x0Leg = 1.0 - marB - ratLegend;
  double x1Leg = 1.0 - marB;
  if(unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote) {
    x0Leg = 1.0 - marL - ratLegend;
    x1Leg = 1.0 - marL;
  }
  double y0Leg;
  if(unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote)
    y0Leg = y0 + histoY * 0.4;
  else
    y0Leg = y0 + histoY * 0.5;
  double y1Leg = y0 + histoY * 1.0 - 0.05;

  SetYMinMax(hh, &pars);
  double xMin = hh[0][0]->GetBinLowEdge(1);
  double xMax = hh[0][0]->GetBinLowEdge(hh[0][0]->GetNbinsX() + 1);

  TH2D* hrangeOrig = new TH2D("", "", 1, xMin + epsAxis, xMax - epsAxis, 1, pars.YMin + epsAxis, pars.YMax - epsAxis);
  utils::ttmd::SetHistoAxisFonts(hrangeOrig);
  hrangeOrig->GetXaxis()->SetTitle(var->GetXTitle());
  hrangeOrig->GetYaxis()->SetTitle("Events");
  if(FlagPaper && TString(hrangeOrig->GetYaxis()->GetTitle()) == "Events" && !unfoldingIOHandler->Suffix.BeginsWith("nj"))
  {
    double width = hh[0][0]->GetBinLowEdge(2) - hh[0][0]->GetBinLowEdge(1);
    int ndig = unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - 1)->NDigits;
    if(unfoldingIOHandler->Suffix == "ptt-pse25")
      ndig = 0;
    TString str = TString::Format("Events / %.*f", ndig, width);
    if(unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - 1)->Units != "")
      str += " " + unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - 1)->Units;
    hrangeOrig->GetYaxis()->SetTitle(str);
  }
  hrangeOrig->GetYaxis()->SetTitleOffset(1.4 * scaleOffsetsY);
  utils::ttmd::ScaleAxisFonts(hrangeOrig->GetXaxis(), scaleAxisFontSize);
  utils::ttmd::ScaleAxisFonts(hrangeOrig->GetYaxis(), scaleAxisFontSize);
  int ndiv = unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - 1)->NDiv;
  if(ndiv < 0)
    ndiv = NDiv;
  hrangeOrig->GetXaxis()->SetNdivisions(ndiv);
  hrangeOrig->GetYaxis()->SetNdivisions(NDiv);
  if(FlagPaper)
  {
    hrangeOrig->GetYaxis()->SetNoExponent(1);
    hrangeOrig->GetYaxis()->SetTitleOffset(2.2 * scaleOffsetsY);
    // special case is 3x3 bins and 2x3x4
    if((unfoldingIOHandler->Dim() == 2 && unfoldingIOHandler->Var(0)->Bins[0].size() == 4 && unfoldingIOHandler->Var(1)->Bins[0].size() == 4) ||
       (unfoldingIOHandler->Dim() == 3 && unfoldingIOHandler->Var(0)->Bins[0].size() == 3 && unfoldingIOHandler->Var(1)->Bins[0].size() == 4 && unfoldingIOHandler->Var(2)->Bins[0].size() == 5))
      hrangeOrig->GetYaxis()->SetTitleOffset(1.8 * scaleOffsetsY);
    //hrange->SetExponentOffset(-10.0, -10.0, "xy");
  }

  if(unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote) {
    hrangeOrig->GetYaxis()->SetNoExponent(0);
  }

  TH2D* hrangeNormOrig = (TH2D*) hrangeOrig->Clone();
  hrangeNormOrig->SetBins(1, xMin + epsAxis, xMax - epsAxis, 1, pars.RatMin + epsAxis, pars.RatMax - epsAxis);
  hrangeNormOrig->GetYaxis()->SetTitle("Ratio ");
  hrangeNormOrig->GetXaxis()->SetLabelOffset(0.02);
  hrangeNormOrig->GetXaxis()->SetTitleOffset(3.2);
  if(FlagPaper)
    hrangeNormOrig->GetXaxis()->SetTitleOffset(3.8);
  hrangeNormOrig->GetYaxis()->SetLabelOffset(0.01 * scaleOffsetsY);
  hrangeNormOrig->GetYaxis()->SetTitleOffset(1.4 * scaleOffsetsY);
  if(FlagPaper)
  {
    hrangeNormOrig->GetYaxis()->SetTitleOffset(2.2 * scaleOffsetsY);
    // special case is 3x3 bins and 2x3x4
    if((unfoldingIOHandler->Dim() == 2 && unfoldingIOHandler->Var(0)->Bins[0].size() == 4 && unfoldingIOHandler->Var(1)->Bins[0].size() == 4) ||
       (unfoldingIOHandler->Dim() == 3 && unfoldingIOHandler->Var(0)->Bins[0].size() == 3 && unfoldingIOHandler->Var(1)->Bins[0].size() == 4 && unfoldingIOHandler->Var(2)->Bins[0].size() == 5))
      hrangeNormOrig->GetYaxis()->SetTitleOffset(1.8 * scaleOffsetsY);
  }

  TLegend* legend = new TLegend();
  // 1D - setup
  if(unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote) {
    legend = new TLegend(x0Leg * 0.86, y0Leg * 1.12, x1Leg * 0.97, y1Leg * 0.97);
    legend->SetMargin(0.17);
    legend->SetTextFont(42);
  }
  else {
    legend = new TLegend(x0Leg * 1.01, y0Leg * 1.01, x1Leg * 0.99, y1Leg * 0.99);
    legend->SetMargin(0.25);
    legend->SetTextSize(30);
    legend->SetTextFont(63);

    if(FlagPaper)
      {
	legend->SetTextSize(38);
	legend->SetTextFont(43);
	legend->SetY1(0.55);
      }
  }

  if(unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote && unfoldingIOHandler->Suffix == "dphitt-xsec8")
    legend = new TLegend(x0Leg * 0.80, y0Leg * 1.12, x1Leg * 0.92, y1Leg * 0.97);
  else if(unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote && unfoldingIOHandler->Suffix == "ytLead-xsec8")
    legend = new TLegend(x0Leg * 0.86, y0Leg * 1.15, x1Leg * 0.98, y1Leg * 1.00);
  else if(unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote && unfoldingIOHandler->Suffix == "ytNLead-xsec8")
    legend = new TLegend(x0Leg * 0.86, y0Leg * 1.15, x1Leg * 0.98, y1Leg * 1.00);

  std::vector<TPad*> vPadToDraw;

  // for further drawing of x title spanned over all pads
  double x0NormTotal = 0.0;
  double x1NormTotal = 0.0;

  for(unsigned int i = 0; i < hh.size(); i++)
  {
    if(i == 0)
    {
      if(FlagPaper)
        legend->AddEntry(hh[i][hh[i].size() - 1], "Data", "p");
      else
        legend->AddEntry(hh[i][hh[i].size() - 1], "Data", "pe");
    }

    if(i > 0)
      x1 += histoX;
    TPad* pad = new TPad("", "", x0, y0, x1, y1);
    //printf("pad: %f %f %f %f\n", x0, y0, x1, y1);
    pad->SetMargin((i == 0) ? axisX / (x1 - x0) : 0.0, 0.0, 0.0, 0.07);
    //printf("TPad %f %f %f %f\n", x0, y0, x1, y1);
    pad->cd();
    pad->SetLogx(pars.LogX);
    pad->SetLogy(pars.LogY);
    TH2D* hrange = (TH2D*) hrangeOrig->Clone();
    if(i != 0)
    {
      hrange->GetYaxis()->SetTitle("");
    }
    hrange->GetXaxis()->SetTitle("");
    hrange->Draw();
    for(int j = hh[i].size() - 2; j >= 0; j--)
    {
      hh[i][j]->Draw("hist0 same");
      if(i == 0)
      {
        legend->AddEntry(hh[i][j], hh[i][j]->GetTitle(), "f");
      }
    }
    // draw uncertainty band
    if(flagDrawBand){
          TGraphAsymmErrors* uncert_band_ptr = utils::ttmd::DrawBand(hhMCTotUD[i], hh[i][2]);
          if(i == 0) legend->AddEntry(uncert_band_ptr, "Syst. unc.", "f");
          uncert_band_ptr->SetFillStyle(3354);
          uncert_band_ptr->SetFillColor(12);
          uncert_band_ptr->SetLineColor(12);
          uncert_band_ptr->SetMarkerStyle(0);
    }
    hh[i][hh[i].size() - 1]->SetMarkerSize(markerSize);
    utils::ttmd::DrawAsGraph(hh[i][hh[i].size() - 1]);
    hrange->Draw("axis same");

    if(unfoldingIOHandler->Dim() > 1)
    {
      if(FlagPaper == 0)
        DrawBinLabel(unfoldingIOHandler, i, axisX, histoX);
      else
        DrawBinLabelFull(unfoldingIOHandler, i, axisX, histoX, "cp");
    }

    TPad* padNorm = new TPad("", "", x0, y0Norm, x1, y1Norm);
    if(i == 0)
      x0NormTotal = x0;
    if(i == (hh.size() - 1))
      x1NormTotal = x1;
    padNorm->SetMargin((i == 0) ? axisX / (x1 - x0) : 0.0, 0.0, axisY / (y1Norm - y0Norm), 0.0);
    padNorm->cd();
    TH1D* hRat = (TH1D*) hh[i][hh[i].size() - 1]->Clone();
    TH1D* hUnity = (TH1D*) hRat->Clone();
    for(int b = 0; b < hUnity->GetNbinsX(); b++)
    {
      hUnity->SetBinContent(b + 1, 1.0);
      hUnity->SetBinError(b + 1, 0.0);
    }
    for(int b = 0; b < hRat->GetNbinsX(); b++)
    {
      double val = hRat->GetBinContent(b + 1) / hh[i][hh[i].size() - 2]->GetBinContent(b + 1);
      double err = hRat->GetBinError(b + 1) / hh[i][hh[i].size() - 2]->GetBinContent(b + 1);
      hRat->SetBinContent(b + 1, val);
      hRat->SetBinError(b + 1, err);
      hRat->SetMarkerStyle(20);
      hRat->SetMarkerColor(1);
      //hRat->SetMarkerSize(0.75);
      hRat->SetMarkerSize(markerSize);
      hRat->SetLineColor(1);
    }
    TH2D* hrangeNorm = (TH2D*) hrangeNormOrig->Clone();
    if(i != 0)
    {
      hrangeNorm->GetYaxis()->SetTitle("");
    }
    if(i != (hh.size() -1 ))
    {
      hrangeNorm->GetXaxis()->SetTitle("");
    }
    // disable x axis title always: will be drawn later spanned over multiple bins
    hrangeNorm->GetXaxis()->SetTitle("");
    hrangeNorm->Draw();
    hUnity->Draw("h same");
    if(vMCTotUD.first && vMCTotUD.second)
    {
      std::pair<TH1D*, TH1D*> hhMCTotUDRat;
      hhMCTotUDRat.first = (TH1D*) hhMCTotUD[i].first;
      hhMCTotUDRat.second = (TH1D*) hhMCTotUD[i].second;
      hhMCTotUDRat.first->Divide(hh[i][hh[i].size() - 2]);
      hhMCTotUDRat.second->Divide(hh[i][hh[i].size() - 2]);
      TGraphAsymmErrors* ratio_unc_band = utils::ttmd::DrawBand(hhMCTotUDRat, hUnity);
      ratio_unc_band->SetFillStyle(3554);
      ratio_unc_band->SetFillColor(12);
      ratio_unc_band->SetLineColor(12);
      ratio_unc_band->SetMarkerStyle(0);

    }
    utils::ttmd::DrawAsGraph(hRat);
    hrangeNorm->Draw("axis same");

    c->cd();
    vPadToDraw.push_back(padNorm);
    vPadToDraw.push_back(pad);
    x0 = x1;
  }

  c->cd();
  for(int i = vPadToDraw.size() - 1; i >= 0; i--)
    vPadToDraw[i]->Draw();
  legend->Draw();
  // draw bottom right label with xaxis title
  if(FlagPaper)
    DrawCommonXTitle(x0NormTotal, y0Norm - 0.02, x1NormTotal, y0Norm + axisY * 0.6, hrangeNormOrig);
  else
    DrawCommonXTitle(x0NormTotal, y0Norm, x1NormTotal, y0Norm + axisY * 0.6, hrangeNormOrig);

  // draw CMS title if needed
  if(FlagPaper > 0) DrawCMSTitle(FlagPaper, unfoldingIOHandler->Dim(),"cp");

  if (current_channel!="") DrawChannelLabel(current_channel,FlagPaper,unfoldingIOHandler->Dim(),"cp", x0Leg, x1Leg, y1Leg);

  utils::ttmd::SaveCanvas(c, fileName);
  delete c;
}

void PlotterDiffXSec::XSec(const UnfoldingIOHandler* unfoldingIOHandler, const std::vector<TH1D*> hDat,
                    std::vector<TH1D*> hTh, const TString& fileName, XSecParams pars)
{
  const std::vector<std::pair<TH1D*, TH1D*> > vDataTotUD;
  XSec(unfoldingIOHandler, hDat, vDataTotUD, hTh, fileName, pars);
}

void PlotterDiffXSec::XSec(const UnfoldingIOHandler* unfoldingIOHandler, const std::vector<TH1D*> hDat, const std::vector<std::pair<TH1D*, TH1D*> > vDataTotUD,
                    std::vector<TH1D*> hTh, const TString& fileName, XSecParams pars, const std::pair<TH1D*, TH1D*> vMCTotUD, bool draw_channel_label)
{
  double markerSize = 0.5;
  if(hDat.size() == 1)
    markerSize = 1.0;

  bool all_channel_comparison_plot = fileName.Contains("all_channel_comparison");
  bool all_years_comparison_plot = fileName.Contains("all_year_comparison");
  bool keep_input_style = all_channel_comparison_plot || all_years_comparison_plot;

  std::vector<double> offset = {};
  if (keep_input_style) offset = {0.35, 0.45, 0.55, 0.65};
  else offset = {0.55, 0.35, 0.45, 0.65};

  std::vector<TH1D*> hh = unfoldingIOHandler->CreateLastDimHistos(0);
  const ZVar* var = unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - 1);
  std::vector<std::vector<TH1D*> > vhD;
  std::vector<std::vector<TH1D*> > vhT;
  int iNorm = -1; // numebr of reference theory histogram
  for(unsigned int i = 0; i < hh.size(); i++)
  {
    vhD.push_back(std::vector<TH1D*>{});
    for(unsigned int t = 0; t < hDat.size(); t++)
    {
      TH1D* h = (TH1D*) hh[i]->Clone();
      unfoldingIOHandler->FillLastDimHisto(h, hDat[t], pars.DivideBW);
      vhD[i].push_back(h);
      vhD[i].back()->SetLineColor(hDat[t]->GetLineColor());
      vhD[i].back()->SetMarkerColor(hDat[t]->GetMarkerColor());
      vhD[i].back()->SetMarkerStyle(hDat[t]->GetMarkerStyle());
      vhD[i].back()->SetMarkerSize(markerSize);
      vhD[i].back()->SetTitle(hDat[t]->GetTitle());
    }
    vhT.push_back(std::vector<TH1D*>{});
    for(unsigned int t = 0; t < hTh.size(); t++)
    {
      TH1D* h = (TH1D*) hh[i]->Clone();
      unfoldingIOHandler->FillLastDimHisto(h, hTh[t], pars.DivideBW);
      h->SetLineColor(hTh[t]->GetLineColor());
      h->SetLineStyle(hTh[t]->GetLineStyle());
      vhT[i].push_back(h);
      //vhT[i].back()->SetLineColor(t + 1);
      vhT[i].back()->SetTitle(hTh[t]->GetTitle());
      if(iNorm < 0 && h->GetLineColor() != 0)
        iNorm = t;
    }
  }
  // std::cout << "************************************ REFERENCE " << iNorm << std::endl;
  assert(iNorm >= 0);

  // create u/d histograms
  bool flagDrawSyst = false;
  if(vDataTotUD.size())
    flagDrawSyst = true;
  std::vector<std::vector<std::pair<TH1D*, TH1D*> > > hhDataTotUD;
  if(flagDrawSyst)
    hhDataTotUD = CreateLastDimUDHistograms(unfoldingIOHandler, vDataTotUD, 0, pars.DivideBW);

  double marL = 0.02;
  //if(FlagPaper)
  //  marL = 0.01;
  double marR = 0.02;
  double marB = 0.02;
  double marT = 0.02;
  double axisX = 0.10;
  double axisY = 0.10;
  double ratNorm = 0.3;
  double ratLegend = 0.20;
  if(unfoldingIOHandler->Suffix.Contains("ratio") && unfoldingIOHandler->Dim() == 1)
    {
      pars.RatMin = 0.8;
      pars.RatMax = 1.2;
    }

  if(FlagPaper && unfoldingIOHandler->Dim() != 1)
  {
    //axisX = 0.12;
    axisX = 0.08;
    if(unfoldingIOHandler->Dim() == 2)
      axisX = 0.107;
    if(unfoldingIOHandler->Dim() == 2)
      axisY = 0.12;
    else
      axisY = 0.10;
    marL = 0.005;
    marR = 0.005;
    marB = 0.002;
    marT = 0.007;
    ratLegend = 0.17;
  }
  else if(unfoldingIOHandler->Dim() == 1 && !FlagAnalysisNote)
    {
      marL = 0.005;
      marT = 0.007;
    }
  else if(unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote)
    {
      marL = 0.005;
      marT = 0.007;
      ratLegend = 0.17;
    }

  bool flagTransparent = 0;

  double scaleAxisFontSize = 1.25;
  if(FlagPaper && unfoldingIOHandler->Dim() != 1)
    scaleAxisFontSize = 1.45 * 1.20;
  double scaleOffsetsY = 1.0;
  if(unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote)
  {
    scaleOffsetsY = 1.05;
  }
  if(unfoldingIOHandler->Dim() == 2)
  {
    //scaleOffsetsY = 3.0;
    scaleOffsetsY = 2.05;
    if(FlagPaper){
      if (fileName.Contains("nj")) scaleOffsetsY = 1.5;
      else scaleOffsetsY = 1.9;
    }
  }
  if(unfoldingIOHandler->Dim() == 3)
  {
    //scaleAxisFontSize = 1.0;
    scaleOffsetsY = 5.0;
    if(FlagPaper)
    {
      //scaleOffsetsY = 4.0;
      if (fileName.Contains("pttmttpttt")) scaleOffsetsY = 1.55;
      else scaleOffsetsY = 2.45;
    }
  }

  utils::ttmd::SetOZStyle();
  TCanvas* c;
  if (unfoldingIOHandler->Dim() == 1 && !FlagAnalysisNote)
    {
      c = new TCanvas("cp", "", PlotPixelX, PlotPixelY);
    }
  else if (unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote)
    {
      c = new TCanvas("cp", "", 800, PlotPixelY);
    }
  else
    {
      c = new TCanvas("cp", "", FlagPaper == 0 ? PlotPixelX : 1400, PlotPixelY);
    }
  //c->SetMargin(0.001, 0.001, 0.001, 0.001);
  //c->SetFillColor(0);
  //c->SetFillStyle(0);
  if(flagTransparent)
  {
    c->SetFillStyle(4000);
    c->SetFillColor(0);
    c->SetFillStyle(0);
    c->SetFrameFillColor(0);
    c->SetFrameFillStyle(0);
  }
  int nhh = hh.size();

  // create u/d histograms
  std::vector<std::pair<TH1D*, TH1D*> > hhMCTotUD;
  bool flagDrawBand = false;
  if(vMCTotUD.first && vMCTotUD.second)
  {
    flagDrawBand = true;
    hhMCTotUD = CreateLastDimUDHistograms(unfoldingIOHandler, vMCTotUD, pars.Binning, pars.DivideBW);
  }

  double histoX;
  if(unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote)
    histoX = (1.0 - marL - marR - axisX) / nhh;
  else
    histoX = (1.0 - marL - marR - axisX - ratLegend) / nhh;
  double histoY = (1.0 - marB - marT - axisY) / (1.0 + ratNorm);
  double histoYNorm = histoY * ratNorm;
  double x0 = marL;
  double x1 = x0 + histoX + axisX;
  double y0Norm = marB;
  //if(FlagPaper)
  //  y0Norm -= 0.01;
  double y1Norm = marB + axisY + histoYNorm;
  double y0 = y1Norm;
  double y1 = y0 + histoY;
  //double epsAxis = 1e-5;

  double x0Leg = 1.0 - marB - ratLegend;
  double x1Leg = 1.0 - marB;
  if(!FlagPaper || unfoldingIOHandler->Dim() == 1)
    x0Leg += 0.005;
  else if(unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote) {
    x0Leg = 1.0 - marL - ratLegend;
    x1Leg = 1.0 - marL;
  }
  double y0Leg;
  if(unfoldingIOHandler->Dim() == 1)
    y0Leg = y0 + histoY * 0.5;
  else
    y0Leg = y0 + histoY * (FlagPaper ? 0.4 : 0.5);
  double y1Leg = y0 + histoY * 1.0 - 0.05;

  SetYMinMax(vhD, &pars);
  double xMin = hh[0]->GetBinLowEdge(1);
  double xMax = hh[0]->GetBinLowEdge(hh[0]->GetNbinsX() + 1);

  double epsXAxis = (xMax - xMin) / 1e5;
  double epsYAxis = (pars.YMax - pars.YMin) / 1e4;
  if(unfoldingIOHandler->Suffix.Contains("ratio") && fileName.Contains("-abs") && unfoldingIOHandler->Dim() == 1)
    {
      epsYAxis = (pars.YMax - pars.YMin) / 1e8;
    }
  TH2D* hrange = new TH2D("", "", 1, xMin + epsXAxis, xMax - epsXAxis, 1, pars.YMin + epsYAxis, pars.YMax - epsYAxis);
  utils::ttmd::SetHistoAxisFonts(hrange);
  hrange->GetXaxis()->SetTitle(var->GetXTitle());
  if(fileName.Contains("-abs"))
    hrange->GetYaxis()->SetTitle(var->GetYTitleXSec());
  else
    hrange->GetYaxis()->SetTitle(var->GetYTitleXSecNorm());
  if(unfoldingIOHandler->DoRj)
    hrange->GetYaxis()->SetTitle("R_{j}");
  hrange->GetYaxis()->SetTitleOffset(1.4 * scaleOffsetsY);
  hrange->GetYaxis()->SetLabelOffset(0.01 * scaleOffsetsY);
  utils::ttmd::ScaleAxisFonts(hrange->GetXaxis(), scaleAxisFontSize);
  utils::ttmd::ScaleAxisFonts(hrange->GetYaxis(), scaleAxisFontSize);
  int ndiv = unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - 1)->NDiv;
  if(ndiv < 0)
    ndiv = NDiv;
  hrange->GetXaxis()->SetNdivisions(ndiv);
  hrange->GetYaxis()->SetNdivisions(NDiv);
  if(FlagPaper && unfoldingIOHandler->Dim() != 1)
  {
    hrange->GetYaxis()->SetNoExponent(1);
    hrange->GetYaxis()->SetTitleOffset(2.2 * scaleOffsetsY);
    // special case is 3x3 bins and 2x3x4
    if((unfoldingIOHandler->Dim() == 2 && unfoldingIOHandler->Var(0)->Bins[0].size() == 4 && unfoldingIOHandler->Var(1)->Bins[0].size() == 4) ||
       (unfoldingIOHandler->Dim() == 3 && unfoldingIOHandler->Var(0)->Bins[0].size() == 3 && unfoldingIOHandler->Var(1)->Bins[0].size() == 4 && unfoldingIOHandler->Var(2)->Bins[0].size() == 5))
      hrange->GetYaxis()->SetTitleOffset(1.8 * scaleOffsetsY);
    //hrange->SetExponentOffset(-10.0, -10.0, "xy");
  }

  TH2D* hrangeNorm = (TH2D*) hrange->Clone();
  hrangeNorm->SetBins(1, xMin + epsXAxis, xMax - epsXAxis, 1, pars.RatMin + epsYAxis, pars.RatMax - epsYAxis);
  hrangeNorm->GetYaxis()->SetTitle("Ratio ");
  hrangeNorm->GetXaxis()->SetLabelOffset(0.02);
  hrangeNorm->GetXaxis()->SetTitleOffset(1.2);
  if(FlagPaper)
  {
    hrangeNorm->GetXaxis()->SetTitleOffset(1.15);
    //printf("set offset %f for %s\n", hrangeNorm->GetXaxis()->GetTitleOffset(), hrangeNorm->GetXaxis()->GetTitle());
  }
  hrangeNorm->GetYaxis()->SetLabelOffset(0.01 * scaleOffsetsY);
  hrangeNorm->GetYaxis()->SetTitleOffset(1.4 * scaleOffsetsY);
  if(FlagPaper && unfoldingIOHandler->Dim() != 1)
  {
    hrangeNorm->GetYaxis()->SetTitleOffset(2.2 * scaleOffsetsY);
    // special case is 3x3 bins and 2x3x4
    if((unfoldingIOHandler->Dim() == 2 && unfoldingIOHandler->Var(0)->Bins[0].size() == 4 && unfoldingIOHandler->Var(1)->Bins[0].size() == 4) ||
       (unfoldingIOHandler->Dim() == 3 && unfoldingIOHandler->Var(0)->Bins[0].size() == 3 && unfoldingIOHandler->Var(1)->Bins[0].size() == 4 && unfoldingIOHandler->Var(2)->Bins[0].size() == 5))
      hrangeNorm->GetYaxis()->SetTitleOffset(1.8 * scaleOffsetsY);
  }

  TLegend* legend = new TLegend();
  // 1D - setup
  if(unfoldingIOHandler->Dim() == 1)
    if(FlagAnalysisNote) legend = new TLegend(x0Leg * 0.85, y0Leg * 1.12, x1Leg * 0.97, y1Leg * 0.97);
    else legend = new TLegend(x0Leg * 1.01 - 0.01, y0Leg * 1.01, x1Leg * 0.99 + 0.01, y1Leg * 0.99);
  // 2D & 3D - setup
  else if(unfoldingIOHandler->Dim() > 1) legend = new TLegend(x0Leg * 1.01 - 0.01, y0Leg * 1.01 - 0.04, x1Leg * 0.99 + 0.01, y1Leg * 0.99 - 0.04);
  else {
      std::cerr << "ERROR in PlotterDiffXSec::XSec -->> Dimension of the xsec is not allowed!" << std::endl;
      exit(1);
  }

  //legend->SetFillStyle(1);
  //legend->SetFillColor(4);
  legend->SetMargin(0.15);
  //legend->SetTextSize(30);
  //legend->SetTextFont(63);
  legend->SetTextFont(62);
  legend->SetTextSize(legend->GetTextSize() * 4.0);

  if(unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote && unfoldingIOHandler->Suffix == "dphitt-xsec8")
    legend = new TLegend(x0Leg * 0.80, y0Leg * 1.12, x1Leg * 0.92, y1Leg * 0.97);
  else if(unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote && unfoldingIOHandler->Suffix == "ytLead-xsec8")
    legend = new TLegend(x0Leg * 0.86, y0Leg * 1.15, x1Leg * 0.98, y1Leg * 1.00);
  else if(unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote && unfoldingIOHandler->Suffix == "ytNLead-xsec8")
    legend = new TLegend(x0Leg * 0.86, y0Leg * 1.15, x1Leg * 0.98, y1Leg * 1.00);
  
  if(FlagPaper){
    legend->SetTextFont(42);
    legend->SetMargin(0.17);
    legend->SetTextSize(legend->GetTextSize() * 7.0);
    if (current_channel!="" && draw_channel_label) DrawChannelLabel(current_channel,FlagPaper,unfoldingIOHandler->Dim(),"xsec", x0Leg, x1Leg, y1Leg);
  }

  std::vector<TPad*> vPadToDraw;

  // for further drawing of x title spanned over all pads
  double x0NormTotal = 0.0;
  double x1NormTotal = 0.0;

  for(unsigned int i = 0; i < hh.size(); i++)
  {
    if(i == 0)
    {
      for(unsigned int j = 0; j < vhD[i].size(); j++)
      {
        if(FlagPaper)
          legend->AddEntry(vhD[i][j], vhD[i][j]->GetTitle(), "p");
        else
          legend->AddEntry(vhD[i][j], vhD[i][j]->GetTitle(), "pe");
      }
    }

    if(i > 0)
      x1 += histoX;

    TPad* pad = new TPad("", "", x0, y0, x1, y1);
    if(flagTransparent)
    {
      pad->SetFillStyle(4000);
      pad->SetFillColor(0);
      pad->SetFillStyle(0);
      pad->SetFrameFillColor(0);
      pad->SetFrameFillStyle(0);
    }
    pad->SetMargin((i == 0) ? axisX / (x1 - x0) : 0.0, 0.0, 0.0, 0.07);
    pad->cd();
    pad->SetLogx(pars.LogX);
    pad->SetLogy(pars.LogY);
    //if(FlagPaper)
    //  pad->SetLogy(1);
    //hrange->Draw();
    TH2D* hrangeClone = (TH2D*) hrange->Clone();
    hrangeClone->GetXaxis()->SetTitle("");
    if(i != 0)
    {
      hrangeClone->GetYaxis()->SetTitle("");
    }
    hrangeClone->GetXaxis()->SetTitle("");
    hrangeClone->Draw();
    for(unsigned int t = 0; t < vhT[i].size(); t++)
    {
      vhT[i][t]->SetLineWidth(hTh[t]->GetLineWidth());
      //vhT[i][t]->SetLineColor(kRed);
      vhT[i][t]->Draw("hist same");
      if(i == 0)
      {
        if(FlagPaper && FlagPaper != 3)
        {
          TString str = vhT[i][t]->GetTitle();
          if(str != "SKIP")
          {
            if(str.Contains("["))
            {
              TString str0 = TString(str.Data(), str.First('[') - 1);
              //printf("str = %s\n", str.Data());
              TString str1 = TString(str.Data() + str.First('[') + 1);
              //printf("str1 = %s\n", str1.Data());
              if(str1.Contains("(")) // PDFs
              {
                str0 += ", #chi^{2}=" + TString(str1.Data(), str1.First(')') + 1);
              }
              else
                str0 += ", #chi^{2}=" + TString(str1.Data(), str1.First('/'));
              str = str0;
              legend->AddEntry(vhT[i][t], str, "l");
            }
            else if(str.Contains("NLO CT14"))
            {
              if(fileName.Contains("mtdep"))
              {
                TString str1 = "NLO CT14";
                TString str2 = str.Data() + 8;
                legend->AddEntry(vhT[i][t], str1 + str2, "l");
              }
              else
              {
                if(!TString(str.Data() + 9).BeginsWith("#alpha"))
                {
                  TString str1 = "NLO CT14";
                  legend->AddEntry(vhT[i][t], str1, "l");
                  TString str2 = str.Data() + 8;
                  legend->AddEntry(vhT[i][t], str2, "l");
                }
                else
                {
                  TString str1 = "NLO CT14";
                  legend->AddEntry(vhT[i][t], str1, "l");
                  TString str2 = "#alpha_{s}=0.118";
                  legend->AddEntry(vhT[i][t], str2, "l");
                  TString str3 = "m_{t}^{pole}=";
                  legend->AddEntry(vhT[i][t], str3, "l");
                }
              }
            }
            else if(str == "NLO")
            {
              TString str1 = "NLO";
              legend->AddEntry(vhT[i][t], str1, "l");
              TString str2 = "m_{t}^{pole} = 172.5 GeV";
              legend->AddEntry(vhT[i][t], str2, "l");
            }
            else
              legend->AddEntry(vhT[i][t], str, "l");
          }
        }
        else {
            //legend->AddEntry(vhT[i][t], vhT[i][t]->GetTitle(), "l");
            TString str = vhT[i][t]->GetTitle();
            if(str.Contains("[")){
                TString str0 = TString(str.Data(), str.First('[') - 1);
                //printf("str = %s\n", str.Data());
                TString str1 = TString(str.Data() + str.First('[') + 1);
                //printf("str1 = %s\n", str1.Data());
                if(str1.Contains("(")){ // PDFs
                    str0 += ", #chi^{2}=" + TString(str1.Data(), str1.First(')') + 1);
                }
                else str0 += ", #chi^{2}=" + TString(str1.Data(), str1.First('/'));
                str = str0;
                legend->AddEntry(vhT[i][t], str, "l");
            }
            else legend->AddEntry(vhT[i][t], vhT[i][t]->GetTitle(), "l");
            }
      }
    }
    // draw uncertainty band
    if(flagDrawBand)
    {
        TGraphAsymmErrors* xsec_uncert_band_ptr = utils::ttmd::DrawBand(hhMCTotUD[i], vhT[i][0]);
        if(i == 0) legend->AddEntry(xsec_uncert_band_ptr, TString(vhT[i][0]->GetTitle(), TString(vhT[i][0]->GetTitle()).First(' ')) + " unc.", "f");
        xsec_uncert_band_ptr->SetFillStyle(1001);
        xsec_uncert_band_ptr->SetFillColor(kRed-9);
        xsec_uncert_band_ptr->SetLineColor(kRed-9);
        xsec_uncert_band_ptr->SetLineWidth(0);
        xsec_uncert_band_ptr->SetMarkerStyle(0);
    }

    for(unsigned int t = 0; t < vhT[i].size(); t++){
    vhT[i][t]->Draw("hist same");
    }

    for(unsigned int j = 0; j < vhD[i].size(); j++)
    {
      if(flagDrawSyst && !keep_input_style)
      {
        vhD[i][j]->SetMarkerStyle(20);
        vhD[i][j]->SetMarkerColor(1);
        utils::ttmd::DrawAsGraph(vhD[i][j], 1, offset[j], "P0");
        utils::ttmd::DrawAsGraph(hhDataTotUD[i][j], vhD[i][j], 1, offset[j], "ZP0");
      }
      else
        utils::ttmd::DrawAsGraph(vhD[i][j], 1, offset[j], "ZP0");
    }
    hrangeClone->Draw("axis same");
    
    if(unfoldingIOHandler->Dim() > 1)
    {
      if(FlagPaper == 0)
        DrawBinLabel(unfoldingIOHandler, i, axisX, histoX);
      else
        DrawBinLabelFull(unfoldingIOHandler, i, axisX, histoX, "xsec");
    }

    TPad* padNorm = new TPad("", "", x0, y0Norm, x1, y1Norm);
    if(i == 0)
      x0NormTotal = x0;
    if(i == (hh.size() - 1))
      x1NormTotal = x1;
    if(flagTransparent)
    {
      padNorm->SetFillStyle(4000);
      padNorm->SetFillColor(0);
      padNorm->SetFillStyle(0);
      padNorm->SetFrameFillColor(0);
      padNorm->SetFrameFillStyle(0);
    }
    padNorm->SetMargin((i == 0) ? axisX / (x1 - x0) : 0.0, 0.0, axisY / (y1Norm - y0Norm), 0.0);
    padNorm->cd();
    // clone and disable some axis
    //hrangeNorm->Draw();
    TH2D* hrangeNormClone = (TH2D*) hrangeNorm->Clone();
    if(i != 0)
    {
      hrangeNormClone->GetYaxis()->SetTitle("");
    }
    if(i != (hh.size() -1 ))
    {
      hrangeNormClone->GetXaxis()->SetTitle("");
    }
    // disable x axis title always: will be drawn later spanned over multiple bins
    hrangeNormClone->GetXaxis()->SetTitle("");
    hrangeNormClone->Draw();

    std::vector<TH1D*> vhT_hRat_vector;
    for(unsigned int j = 0; j < vhT[i].size(); j++)
    {
      TH1D* hRat = (TH1D*) vhT[i][j]->Clone();
      for(int b = 0; b < hRat->GetNbinsX(); b++)
      {
        hRat->SetBinContent(b + 1, hRat->GetBinContent(b + 1) / vhT[i][iNorm]->GetBinContent(b + 1));
        hRat->SetBinError(b + 1, hRat->GetBinError(b + 1) / vhT[i][iNorm]->GetBinContent(b + 1));
      }
      hRat->SetMarkerStyle(vhT[i][j]->GetMarkerStyle());
      hRat->SetLineColor(vhT[i][j]->GetLineColor());
      hRat->SetMarkerColor(vhT[i][j]->GetMarkerColor());
      hRat->SetMarkerSize(markerSize);
      hRat->SetLineWidth(hTh[j]->GetLineWidth());
      vhT_hRat_vector.push_back(hRat);
      // hRat->Draw("hist0 same");
    }
    for(unsigned int j = 0; j < vhD[i].size(); j++)
    {
      TH1D* hRat = (TH1D*) vhD[i][j]->Clone();
      TH1D* hUnity = (TH1D*) hRat->Clone();
      hUnity->SetLineColor(2);
      for(int b = 0; b < hUnity->GetNbinsX(); b++)
      {
        hUnity->SetBinContent(b + 1, 1.0);
        hUnity->SetBinError(b + 1, 0.0);
        //hUnity->SetLineColor(kRed);
        hRat->SetBinContent(b + 1, hRat->GetBinContent(b + 1) / vhT[i][iNorm]->GetBinContent(b + 1));
        hRat->SetBinError(b + 1, hRat->GetBinError(b + 1) / vhT[i][iNorm]->GetBinContent(b + 1));
      }
      hRat->SetMarkerStyle(vhD[i][j]->GetMarkerStyle());
      hRat->SetLineColor(vhD[i][j]->GetLineColor());
      hRat->SetMarkerColor(vhD[i][j]->GetMarkerColor());
      hRat->SetMarkerSize(markerSize);

      /*
      if(!FlagPaper || (!fileName.Contains("xsec-fit-") && !fileName.Contains("xsec-topptproblem-")))
        hUnity->Draw("h same");*/

      if(vMCTotUD.first && vMCTotUD.second)
      {
        std::pair<TH1D*, TH1D*> hhMCTotUDRat;
        hhMCTotUDRat.first = (TH1D*) hhMCTotUD[i].first;
        hhMCTotUDRat.second = (TH1D*) hhMCTotUD[i].second;
        hhMCTotUDRat.first->Divide(vhT[i][iNorm]);
        hhMCTotUDRat.second->Divide(vhT[i][iNorm]);
        TGraphAsymmErrors* xsec_ratio_uncert_band_ptr = utils::ttmd::DrawBand(hhMCTotUDRat, hUnity);
        xsec_ratio_uncert_band_ptr->SetFillStyle(1001);
        xsec_ratio_uncert_band_ptr->SetFillColor(kRed-9);
        xsec_ratio_uncert_band_ptr->SetLineColor(kRed-9);
        xsec_ratio_uncert_band_ptr->SetMarkerStyle(0);
      }
      // 13.12.2018 suppress unity black line for plots with "NLO fit" (which have green ratio histogram)
      if(!FlagPaper || (!fileName.Contains("xsec-fit-") && !fileName.Contains("xsec-topptproblem-")))
          hUnity->Draw("h same");

      for (auto i_vhT_hRat: vhT_hRat_vector) i_vhT_hRat->Draw("hist0 same");

      if(flagDrawSyst)
      {
        utils::ttmd::DrawAsGraph(hRat, 1, offset[j], "P0");
        std::pair<TH1D*, TH1D*> hhRatTotUD = std::make_pair<TH1D*, TH1D*>(NULL, NULL);
        hhRatTotUD.first = (TH1D*) hhDataTotUD[i][j].first->Clone();
        hhRatTotUD.first->Divide(vhT[i][iNorm]);
        hhRatTotUD.second = (TH1D*) hhDataTotUD[i][j].second->Clone();
        hhRatTotUD.second->Divide(vhT[i][iNorm]);
        utils::ttmd::DrawAsGraph(hhRatTotUD, hRat, 1, offset[j], "ZP0");
      }
      else
        utils::ttmd::DrawAsGraph(hRat, 1, offset[j], "ZP0");
    }
    hrangeNormClone->Draw("axis same");

    c->cd();
    vPadToDraw.push_back(padNorm);
    vPadToDraw.push_back(pad);
    x0 = x1;
  }

  c->cd();
  for(int i = vPadToDraw.size() - 1; i >= 0; i--)
    vPadToDraw[i]->Draw();
  if(FlagPaper && vhT[0].size() > 5)
    legend->SetY1(0.35);
  if(FlagPaper && fileName.Contains("xsecmtdep"))
    legend->SetY1(legend->GetY1() - 0.05);
  if(FlagPaper && fileName.Contains("topptproblem"))
    legend->SetY1(legend->GetY1() - 0.10);
  legend->Draw();
  // draw bottom right label with xaxis title
  if(FlagPaper)
  {
    if(unfoldingIOHandler->Dim() == 2)
      DrawCommonXTitle(x0NormTotal, y0Norm + 0.01, x1NormTotal, y0Norm + axisY * 0.3, hrangeNorm);
    else
      DrawCommonXTitle(x0NormTotal, y0Norm - 0.01, x1NormTotal, y0Norm + axisY * 0.3, hrangeNorm);
  }
  else
    DrawCommonXTitle(x0NormTotal, y0Norm, x1NormTotal, y0Norm + axisY * 0.6, hrangeNorm);

  // draw CMS title if needed
  if(FlagPaper > 0)
    DrawCMSTitle(FlagPaper, unfoldingIOHandler->Dim(),"xsec");

  utils::ttmd::SaveCanvas(c, fileName);
  delete c;
}

void PlotterDiffXSec::PlotVars(const UnfoldingIOHandler* unfoldingIOHandler, const std::vector<TH1D*>& vhVarGen, const std::vector<TH1D*>& vhVarRec, const TString& fileName, const int flagAbs, const int flagEig/* = 0*/)
{
  std::vector<std::vector<TH1D*> > vhVarDat;
  PlotVars(unfoldingIOHandler, vhVarDat, vhVarGen, vhVarRec, fileName, flagAbs, flagEig);
}

void PlotterDiffXSec::PlotVars(const UnfoldingIOHandler* unfoldingIOHandler, const std::vector<std::vector<TH1D*> >& vhVarDat, const std::vector<TH1D*>& vhVarGen, const std::vector<TH1D*>& vhVarRec, const TString& fileName, const int flagAbs, const int flagEig/* = 0*/)
{
  // label for two line types
  std::pair<TString, TString> pairVarType = flagEig ? std::make_pair("up", "dn") : std::make_pair("gen", "rec");

  double lineWidth = 1.0;
  double sizeMarker = 1.0;
  double eps = 1e-4;
  double deltaY = flagAbs ? 0.2 : 0.1;
  double nomY = 1.0;
  double minY = nomY - deltaY;
  double maxY = nomY + deltaY;
  double legY0 = 0.8;
  
  int nalgs = vhVarDat.size();
  bool isData = (nalgs > 0);
  int nMarkers = nalgs * (vhVarRec.size() - 1);
  double offsetMarkerWidth = 0.6 / nMarkers;
  double offsetMarkerStart = 0.2;

  TCanvas* cc = new TCanvas("", "", 800, 600);
  TPad* c = new TPad("", "", 0.0, 0.0, 1.0, legY0 - 0.01);
  c->SetTopMargin(0.0);
  c->SetRightMargin(0.03);
  int nbins = unfoldingIOHandler->HGenC()->GetNbinsX();

  double offbin = 0.5;
  TH2D* hrange = new TH2D("", "", 1, offbin, nbins + offbin, 1, minY + eps, maxY - eps);
  //TString xTitle = ZVar::TitleWithExtensionCombined(unfoldingIOHandler->Vars());
  TString xTitle = ZVar::TitleFull(unfoldingIOHandler->Vars());
  hrange->GetXaxis()->SetTitle(TString::Format("Bin %s", xTitle.Data()));
  hrange->GetXaxis()->SetNdivisions(216);
  hrange->GetXaxis()->SetTickLength(0.0);
  TString yTitle = "Ratio to nominal";
  if(!FlagPaper)
  {
    if(isData){
      if(flagAbs)
        yTitle += " (abs. data)";
      else
        yTitle += " (nor. data)";
    }
  }
  hrange->GetYaxis()->SetTitle(yTitle);
  hrange->GetYaxis()->SetTitleOffset(1.42);
  std::vector<double> vecyy(2, nomY);
  TPolyLine plineUnity(2, &(std::vector<double>{offbin, nbins + offbin}[0]), &vecyy[0]);
  std::vector<TPolyLine*> plineBins = CreatePolylineBins(unfoldingIOHandler, nbins, minY, maxY);

  TLegend* leg = new TLegend(FlagPaper ? 0.08 : 0.00, legY0, 0.96 * (0.7 + vhVarDat.size() * 0.10), 0.98);
  leg->SetFillColor(0);
  leg->SetNColumns(3 + nalgs);
  leg->SetEntrySeparation(0.03);
  leg->SetColumnSeparation(0.03);
  leg->SetMargin(0.33);
  //leg->SetTextSize(leg->GetTextSize() * 2.0);
  leg->SetTextAlign(12);

  c->cd();
  hrange->Draw();
  for(unsigned int v = 1; v < vhVarRec.size(); v++)
  {
    int color = utils::ttmd::GetDistinctColor(v);
    // gen
    TH1D* hGenRat = new TH1D(*vhVarGen[v]);
    hGenRat->Divide(vhVarGen[0]);
    hGenRat->SetLineColor(color);
    hGenRat->SetLineWidth(lineWidth);
    hGenRat->Draw("hist0 same");
    // rec
    TH1D* hRecRat = new TH1D(*vhVarRec[v]);
    hRecRat->Divide(vhVarRec[0]);
    hRecRat->SetLineColor(color);
    hRecRat->SetLineWidth(lineWidth);
    hRecRat->SetLineStyle(2);
    hRecRat->Draw("hist0 same");
    // legend entries
    //size_t posBracket = TString(vhVarRec[v]->GetTitle()).First('[');
    Ssiz_t posBracket = TString(vhVarRec[v]->GetTitle()).First('[');

    // if '[' is present in histo title, then it is considered as pair of variations,
    // otherwise only one variation will be added to the legend
    const bool flagTwoVars = (posBracket != TString::kNPOS);
    // 16.07.18 if "[[" is present in histo title, then it is considered as single variation; one of '[' is omitted in the legend
    const bool flagSingleVar = flagTwoVars && *(vhVarRec[v]->GetTitle() + posBracket + 1) == '[';
    TString header = flagTwoVars ? TString(hRecRat->GetTitle(), posBracket) : TString(hRecRat->GetTitle());
    //printf("title: %s %d %d %d\n", vhVarRec[v]->GetTitle(), posBracket, flagTwoVars, flagSingleVar);
    if(flagTwoVars)
    {
      leg->AddEntry((TObject*)NULL, header, "");
      if(!flagSingleVar)
      {
        leg->AddEntry(hGenRat, TString::Format("%s%s", pairVarType.first.Data(), hGenRat->GetTitle() + posBracket), "l");
        leg->AddEntry(hRecRat, TString::Format("%s%s", pairVarType.second.Data(), hRecRat->GetTitle() + posBracket), "l");
      }
      else
      {
        leg->AddEntry(hGenRat, TString::Format("%s", hGenRat->GetTitle() + posBracket + 1), "l");
        leg->AddEntry((TObject*)(NULL), "", "");
      }
    }
    else
      leg->AddEntry(hGenRat, hRecRat->GetTitle(), "l");
    // data if present
    for(int a = 0; a < nalgs; a++)
    {
      TH1D* hRat = new TH1D(*vhVarDat[a][v]);
      hRat->Divide(vhVarDat[a][0]);
      int colorGraph = color + 3 * (a - (nalgs - 1) / 2);
      hRat->SetLineColor(colorGraph);
      hRat->SetMarkerColor(colorGraph);
      hRat->SetMarkerSize(sizeMarker);
      hRat->SetMarkerStyle(24 + a);
      utils::ttmd::DrawAsGraph(hRat, 0, offsetMarkerStart + offsetMarkerWidth * ((v - 1) * nalgs + a));
      leg->AddEntry(hRat, vhVarDat[a][v]->GetTitle(), "p");
    }
  }

  c->cd();
  plineUnity.Draw();
  //if(isData)
  for(auto pl : plineBins)
    pl->Draw();
  hrange->Draw("axis same");

  cc->cd();
  c->Draw();
  leg->Draw();
  utils::ttmd::SaveCanvas(cc, fileName);
}

std::vector<TString> PlotterDiffXSec::GetUncSummaryList(const std::vector<TString>& selectVar, TString& year)
{

    if (year=="fullRun2") {
        std::vector<TString> vType;
        for(std::vector<TString>::const_iterator it = selectVar.begin(); it != selectVar.end(); it++) {
          bool isYearToYearCorr = PlotterConfigurationHelper::IsYearToYearCorr("fullRun2", "all", *it);
          std::vector<TString> yearCorrList = PlotterConfigurationHelper::GetYearCorrList("fullRun2", this->configHelper->GetCorrYearValue(), *it, isYearToYearCorr);
          if (yearCorrList.size() != 0)
            vType.insert(vType.end(), yearCorrList.begin(), yearCorrList.end());
          else
            vType.push_back(*it);
        }
        return vType;
    }
  else return selectVar;
}

void PlotterDiffXSec::PlotUncSummary(const UnfoldingIOHandler* unfoldingIOHandler, const std::vector<std::vector<TH1D*> >& vhMarkers,
                              const std::vector<TH1D*>& vhTotal, const std::vector<TH1D*>& vhDashed,
			      const TString& fileName, const UncSummaryParams& params)
{
  assert(vhTotal.size() == 2);
  assert(vhDashed.size() == 4 || vhDashed.size() == 2 || vhDashed.size() == 0);
  for(unsigned int i = 0; i < vhMarkers.size(); i++)
    assert(vhMarkers[i].size() == 2);

  double lineWidth = 1.0;
  double sizeMarker = 1.0;
  int lineColor = 1;
  double eps = 1e-4;
  double minY = -1 * params.Rat;
  double maxY = params.Rat;
  if(maxY < 0.0) // use default RatMin and RatMax
  {
    if(params.FlagAbs)
    {
      minY = -20.0;
      maxY = +20.0;
    }
    else
    {
      minY = -20.0;
      maxY = +20.0;
    }
    if(fileName.Contains("jes-"))
    {
      minY = -15.0;
      maxY = +15.0;

    }
    if(unfoldingIOHandler->Dim() == 3)
    {
      if(unfoldingIOHandler->Var(0)->Bins[0].size() == 3 && unfoldingIOHandler->Var(1)->Bins[0].size() == 4)
      {
        // njmttytt-b2-mtt3
        minY *= 1;
        maxY *= 1;
      }
      else
      {
        minY *= 2;
        maxY *= 2;
      }
    }
  }
  double legY0 = 0.9;

  bool flagPlotMarkers = (vhMarkers.size() > 0);
  int nMarkers = flagPlotMarkers ? vhMarkers.size() : 0;
  double offsetMarkerWidth = flagPlotMarkers ? (0.6 / (nMarkers - 1)) : 0;
  double offsetMarkerStart = (flagPlotMarkers && vhMarkers[0].size() > 1) ? 0.2 : 0.5;

  bool flagPlotDashed = (vhDashed.size() > 0);

  int widthX = 800;
  int widthY = 600;
  if(unfoldingIOHandler->Dim() == 3)
    widthX = 1200;
  TCanvas* cc = new TCanvas("", "", widthX, widthY);
  TPad* c = new TPad("", "", 0.0, 0.0, 1.0, legY0 - 0.01);
  c->SetTopMargin(0.0);
  c->SetBottomMargin(FlagPaper ? 0.18 : 0.13);
  c->SetLeftMargin(0.10);
  c->SetRightMargin(0.03);
  int nbins = unfoldingIOHandler->HGenC()->GetNbinsX();

  double offbin = 0.5;
  TH2D* hrange = new TH2D("", "", 1, offbin, nbins + offbin, 1, minY + eps, maxY - eps);
  //TString xTitle = ZVar::TitleWithExtensionCombined(unfoldingIOHandler->Vars());
  TString xTitle = ZVar::TitleFull(unfoldingIOHandler->Vars());
  /*if(FlagPaper && unfoldingIOHandler->Dim() == 3)
  {
    int nb = unfoldingIOHandler->Var(0)->Bins[0].size() - 1;
    if(nb == 2)
      xTitle.ReplaceAll("N_{jet}", "N_{jet}^{0,1+}");
    else if(nb == 3)
      xTitle.ReplaceAll("N_{jet}", "N_{jet}^{0,1,2+}");
  }*/
  hrange->GetXaxis()->SetTitle(TString::Format("Bin %s", xTitle.Data()));
  hrange->GetXaxis()->SetNdivisions(216);
  hrange->GetXaxis()->SetTickLength(0.0);
  TString yTitle = "Uncertainty [%]";
  if(!FlagPaper)
  {
    if(params.FlagAbs)
      yTitle += " (abs. data)";
    else
      yTitle += " (nor. data)";
  }
  hrange->GetYaxis()->SetTitle(yTitle);
  hrange->GetYaxis()->SetTitleOffset(1.0);
  hrange->GetXaxis()->SetTitleOffset(1.1);
  utils::ttmd::ScaleHistoFonts(hrange, 1.30);
  if(FlagPaper)
  {
    utils::ttmd::ScaleHistoFonts(hrange, 1.45);
    hrange->GetYaxis()->SetTitleOffset(0.6);
  }
  std::vector<double> vecyy(2, 0.0);
  TPolyLine plineUnity(2, &(std::vector<double>{offbin, nbins + offbin}[0]), &vecyy[0]);
  std::vector<TPolyLine*> plineBins = CreatePolylineBins(unfoldingIOHandler, nbins, minY, maxY);

  double legWidthMin = 0.6;
  double legWidthMax = 0.95;
  double legOffsettoX0 = 0.0;
  double legWidth = (vhMarkers.size() == 0) ? legWidthMin : std::min((legWidthMin + 0.2 * vhMarkers.size()), legWidthMax);
  double legX0 = 0.50 - legOffsettoX0 - legWidth / 2.0;
  double legX1 = 0.50 - legOffsettoX0 + legWidth / 2.0;
  TLegend* leg = new TLegend(FlagPaper ? legX0 + 0.08 : legX0, legY0, legX1, 0.98);
  leg->SetFillColor(0);
  leg->SetNColumns((3 + vhMarkers.size()) / params.NLegendRows);
  leg->SetEntrySeparation(0.0);
  leg->SetColumnSeparation(0.0);
  leg->SetMargin(0.20);
  leg->SetTextAlign(12);

  c->cd();
  hrange->Draw();
  for(unsigned int v = 0; v < vhTotal.size(); v++)
  {
    // gen
    TH1D* hTotal = new TH1D(*vhTotal[v]);
    hTotal->SetLineColor(lineColor);
    hTotal->SetLineWidth(lineWidth * 2);
    hTotal->Draw("hist0 same");
    if(v == 0)
      leg->AddEntry(hTotal, vhTotal[0]->GetTitle(), "l");
    // rec
    if(flagPlotDashed)
    {
      for(unsigned int vd = 0; vd < vhDashed.size() / vhTotal.size(); vd++)
      {
        TH1D* hDashed = new TH1D(*vhDashed[v + vd * 2]);
        hDashed->SetLineColor(lineColor);
        hDashed->SetLineWidth(lineWidth);
        if(vd == 0)
          hDashed->SetLineStyle(2);
        else if(vd == 1)
        {
          hDashed->SetLineStyle(3);
          //hDashed->SetLineStyle(1);
          //hDashed->SetLineWidth(1);
        }
        hDashed->Draw("hist0 same");
        if(v == 0)
          leg->AddEntry(hDashed, vhDashed[vd * 2]->GetTitle(), "l");
      }
    }
    // markers if present
    for(unsigned int a = 0; a < vhMarkers.size(); a++)
    {
      int colorGraph = utils::ttmd::GetDistinctColor(a);
      TH1D* hGraph = new TH1D(*vhMarkers[a][v]);
      hGraph->SetLineColor(colorGraph);
      hGraph->SetMarkerColor(colorGraph);
      hGraph->SetMarkerSize(sizeMarker);
      hGraph->SetMarkerStyle(vhMarkers[a][0]->GetMarkerStyle());
      //double offset = 0.5;
      //if(a > 0)
      //  offset = 0.3 + 0.2 * (a - 1);
      utils::ttmd::DrawAsGraph(hGraph, 0, offsetMarkerStart + offsetMarkerWidth * a);
      //utils::ttmd::DrawAsGraph(hGraph, 0, offsetMarkerStart + offsetMarkerWidth * ((v / 2) * vhMarkers.size() + a));
      //utils::ttmd::DrawAsGraph(hGraph, 0, 0.5, "ZP0C");
      //utils::ttmd::DrawAsGraph(hGraph, 0, 0.5, "");
      if(v == 0)
        leg->AddEntry(hGraph, vhMarkers[a][v]->GetTitle(), "p");
    }
  }

  c->cd();
  plineUnity.Draw();
  //if(isData)
  for(auto pl : plineBins)
    pl->Draw();
  hrange->Draw("axis same");

  cc->cd();
  c->Draw();
  leg->Draw();
  utils::ttmd::SaveCanvas(cc, fileName);
}


void PlotterDiffXSec::AddXSec(UnfoldingIOHandler* unfoldingIOHandler)
{
  //gTreeDir  = gBaseDir + "../plainTree-oz/" + gPlainTreeVersion;
  //gPlotsDir = gBaseDir + "./plots-oz";
  vXSec.push_back(unfoldingIOHandler);
}

void PlotterDiffXSec::AddXSec(std::vector<UnfoldingIOHandler*> vXSec_input)
{
  for (auto xsec: vXSec_input){
    AddXSec(xsec);
  }
}

void PlotterDiffXSec::SetVar(const TString& var)
{
  Var = var;
}

std::vector<TString> PlotterDiffXSec::CreateSuffVector()
{
  return std::vector<TString> { "nor", "abs" };
}

std::vector<ZUnfold*> PlotterDiffXSec::CreateZUnfoldVector(const int n)
{
  std::vector<ZUnfold*> vUnf;
  // BBB
  vUnf.push_back(new ZUnfold(configHelper));
  vUnf.back()->Remark = "BBB";
  vUnf.back()->Suffix = "bbb";
  vUnf.back()->MarkerStyle = 25;
  vUnf.back()->Color = 1;
  if(n == 1)
    return vUnf;
  // no reg
  vUnf.push_back(new ZUnfold(configHelper));
  vUnf.back()->regChoice = ZUnfold::RegNo;
  vUnf.back()->Remark = "NoReg";
  vUnf.back()->Suffix = "unr";
  vUnf.back()->MarkerStyle = 26;
  vUnf.back()->Color = 4;
  if(n == 2)
    return vUnf;
  // reg min rho average
  vUnf.push_back(new ZUnfold(configHelper));
  vUnf.back()->Remark = "MinRhoAvg";
  vUnf.back()->Suffix = "reg";
  vUnf.back()->MarkerStyle = 24;
  vUnf.back()->Color = 2;
  if(n == 3)
    return vUnf;
  // BBB0
  vUnf.push_back(new ZUnfold(configHelper));
  vUnf.back()->Remark = "BBB0";
  vUnf.back()->Suffix = "bbb0";
  vUnf.back()->MarkerStyle = 25;
  vUnf.back()->Color = 6;
  if(n == 4)
    return vUnf;

  return vUnf;
}


// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>> data >>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void PlotterDiffXSec::LoopDataTreePlain(const TString& ch)
{
  // reset
  for(unsigned int x = 0; x < vXSec.size(); x++)
    vXSec[x]->Reset("d");
  // tree
  TString dir = "Nominal";
  // TODO find better way to set directory (only these variations apply in data)
  if(Var == "MASS_KINRECO_UP" || Var == "MASS_KINRECO_DOWN" || Var == "MASS_KINRECO3GEV_UP" || Var == "MASS_KINRECO3GEV_DOWN")
    dir = Var;
  //if(Var == "MASS_UP")
  //  dir = "MASS_KINRECO3GEV_UP";
  //if(Var == "MASS_DOWN")
  //  dir = "MASS_KINRECO3GEV_DOWN";
  TreePlain* tree = configHelper->GetTreeDataPlain(ch, dir);
  if(configHelper->gTTBARDebug)
    printf("Data events: %lld\n", tree->TTreePtr->GetEntries());
  for(unsigned int x = 0; x < vXSec.size(); x++)
    vXSec[x]->InitTree(tree, true);
  for(int e = 0; e < (MaxEvents ? MaxEvents : tree->TTreePtr->GetEntries()); e++)
    {
      tree->TTreePtr->GetEntry(e);
      for(unsigned int x = 0; x < vXSec.size(); x++)
	vXSec[x]->FillData(tree);
    }
  for(unsigned int x = 0; x < vXSec.size(); x++)
    vXSec[x]->ClearTree(tree, true);
  delete tree;
}

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>> MC signal gen only >>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void PlotterDiffXSec::StandaloneLoopMCTreeGenPlain(const TString& ch, const int binning, TH1D* h, UnfoldingIOHandler* unfoldingIOHandler, const TString fileIn)
{
  h->Reset();
  bool flagAddfiducialInfo = false;
  if (unfoldingIOHandler->DoParticle)
      flagAddfiducialInfo = true;

  TreePlain* treeGen = NULL;
  if(fileIn == "")
    treeGen = configHelper->GetTreeMCGenPlain(Var, ch, Var, "Nominal", true, flagAddfiducialInfo);
  else
    {
      // manual setup (used for npcorr)
      treeGen = new TreePlain(configHelper->gPlainTreeName0Plain, true);
      treeGen->Add(fileIn);
      treeGen->InitReading();
    }
  if(configHelper->gTTBARDebug)
    {
      printf("MC events GEN: %lld\n", treeGen->TTreePtr->GetEntries());
    }
  unfoldingIOHandler->InitTree(treeGen);
  for(int e0 = 0; e0 < (MaxEvents ? MaxEvents : treeGen->TTreePtr->GetEntries()); e0++)
    {
      //if(e0 % 1000 == 0) printf("e0 = %d\n", e0);
      treeGen->TTreePtr->GetEntry(e0);
      double rew = unfoldingIOHandler->Reweight(treeGen);
      if(unfoldingIOHandler->DoParticle) {
	if(treeGen->Vars.inFiducialPhaseSpace)
	  unfoldingIOHandler->StandaloneFillMCGen(treeGen, binning, h, rew);
      }
      else
	unfoldingIOHandler->StandaloneFillMCGen(treeGen, binning, h, rew);
    }
  unfoldingIOHandler->ClearTree(treeGen);
  delete treeGen;
}

void PlotterDiffXSec::StandaloneLoopMCTreeGenPlain(const TString& ch, const int binning, std::vector<TH1D*> vH, std::vector<UnfoldingIOHandler*> vXSec, const TString fileIn, const TString weight)
{
  assert(vH.size());
  assert(vH.size() == vXSec.size());
  for(auto& h : vH)
    h->Reset();
  bool flagAddfiducialInfo = false;
  for(unsigned int x = 0; x < vXSec.size(); x++) {
    if (vXSec[x]->DoParticle)
      flagAddfiducialInfo = true;
  }

  TreePlain* treeGen = NULL;
  if(fileIn == "")
    treeGen = configHelper->GetTreeMCGenPlain(Var, ch, "Nominal", "Nominal", true, flagAddfiducialInfo);
  else
    {
      // manual setup (used for npcorr)
      treeGen = new TreePlain(configHelper->gPlainTreeName0Plain, true, weight);
      treeGen->Add(fileIn);
      treeGen->InitReading();
    }
  if(configHelper->gTTBARDebug)
    {
      printf("MC events GEN: %lld\n", treeGen->TTreePtr->GetEntries());
    }
  for(auto& unfoldingIOHandler : vXSec)
    unfoldingIOHandler->InitTree(treeGen);
  std::vector<double> vRew(vXSec.size());
  for(int e0 = 0; e0 < (MaxEvents ? MaxEvents : treeGen->TTreePtr->GetEntries()); e0++)
    {
      for(unsigned int x = 0; x < vXSec.size(); x++)
	vRew[x] = vXSec[x]->Reweight(treeGen);
      //if(e0 % 100000 == 0)
      //  printf("e0 = %d\n", e0);
      treeGen->TTreePtr->GetEntry(e0);
      for(size_t x = 0; x < vH.size(); x++) {
	if (vXSec[x]->DoParticle) {
	  if(treeGen->Vars.inFiducialPhaseSpace)
	    vXSec[x]->StandaloneFillMCGen(treeGen, binning, vH[x], vRew[x]);
	}
	else
	  vXSec[x]->StandaloneFillMCGen(treeGen, binning, vH[x], vRew[x]);
      }
    }
  for(auto& unfoldingIOHandler : vXSec)
    unfoldingIOHandler->ClearTree(treeGen);
  delete treeGen;
}

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>> MC signal >>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void PlotterDiffXSec::LoopMCTreePlain(const TString& ch)
{
  // reset
  bool flagAddfiducialInfo = false;
  for(unsigned int x = 0; x < vXSec.size(); x++) {
    vXSec[x]->Reset("rg");
    vXSec[x]->Reset("ops");
    if (vXSec[x]->DoParticle)
      flagAddfiducialInfo = true;
  }

  // tree
  TString dir = configHelper->GetVarDirTreePlain(Var, zMode);
  TString wSyst = zMode ? configHelper->GetVarWeight(Var) : "Nominal";
  TreePlain* treeRec = configHelper->GetTreeMCRecPlain(Var, ch, dir, wSyst, zMode, flagAddfiducialInfo);
  TreePlain* treeGen = configHelper->GetTreeMCGenPlain(Var, ch, dir, wSyst, zMode, flagAddfiducialInfo);
  if(configHelper->gTTBARDebug)
    {
      printf("MC events GEN: %lld\n", treeGen->TTreePtr->GetEntries());
      printf("MC events REC: %lld\n", treeRec->TTreePtr->GetEntries());
    }
  for(unsigned int x = 0; x < vXSec.size(); x++)
    {
      if(vXSec[x]->DoPSE || vXSec[x]->DoXSec)
        {
          vXSec[x]->InitTree(treeGen);
          //vXSec[x]->InitTreeForReweighting(treeGen);
        }
      vXSec[x]->InitTree(treeRec, true);
    }
  int e = 0;
  bool flagRecUsed = true;
  std::vector<double> vRew(vXSec.size());
  for(int e0 = 0; e0 < (MaxEvents ? MaxEvents : treeGen->TTreePtr->GetEntries()); e0++)
    //for(int e0 = 0; e0 < 30; e0++)
    {
      //treeGen->Vars.b_eventCounter->GetEntry(e0);
      treeGen->TTreePtr->GetEntry(e0);
      for(unsigned int x = 0; x < vXSec.size(); x++)
	vRew[x] = vXSec[x]->Reweight(treeGen);
      if(flagRecUsed)
        {
          //treeRec->Vars.b_eventCounter->GetEntry(e);
          treeRec->TTreePtr->GetEntry(e);
          e++;
          flagRecUsed = false;
        }
      //printf("%d %d\n", tree->entry, tree0->entry);

      if(treeRec->Vars.eventCounter == treeGen->Vars.eventCounter) // rec level
        {
          //treeGen->TTreePtr->GetEntry(e0);
          //treeRec->Vars.b_eventCounter->GetEntry(e);
          flagRecUsed = true;

	  for(unsigned int x = 0; x < vXSec.size(); x++) {
	    if (vXSec[x]->DoParticle) {
	      if(treeRec->Vars.inFiducialPhaseSpace) {
		if(vXSec[x]->DoPSE || vXSec[x]->DoXSec)
		  vXSec[x]->FillMC(treeRec, treeGen, vRew[x]);
		else
		  vXSec[x]->FillMC(treeRec, NULL, vRew[x]);
	      }
	      else if(!treeRec->Vars.inFiducialPhaseSpace) {
		//fill out of space background for particle-level cross sections
		vXSec[x]->FillOutOfSpace(treeRec);
	      }
	    }
	    else {
	      if(vXSec[x]->DoPSE || vXSec[x]->DoXSec)
		vXSec[x]->FillMC(treeRec, treeGen, vRew[x]);
	      else
		vXSec[x]->FillMC(treeRec, NULL, vRew[x]);
	    }
	  }
	}
      else // fill gen level only
        {
	  for(unsigned int x = 0; x < vXSec.size(); x++) {
	    //handle underflow for events that are in the fiducial phase space at generator-level but not reconstructed at detector-level
	    if (vXSec[x]->DoParticle) {
	      if((vXSec[x]->DoPSE || vXSec[x]->DoXSec) && treeGen->Vars.inFiducialPhaseSpace)
		vXSec[x]->FillMC(NULL, treeGen, vRew[x]);
	    }
	    else {
	      if(vXSec[x]->DoPSE || vXSec[x]->DoXSec)
		vXSec[x]->FillMC(NULL, treeGen, vRew[x]);
	    }
	  }
	}
    }
  for(unsigned int x = 0; x < vXSec.size(); x++)
    {
      vXSec[x]->ClearTree(treeRec, true);
      if(vXSec[x]->DoPSE || vXSec[x]->DoXSec)
	vXSec[x]->ClearTree(treeGen);
    }
  delete treeRec;
  delete treeGen;
}

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>> MC background >>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void PlotterDiffXSec::LoopMCBkgTreePlain(const TString& ch)
{
  // reset
  for(unsigned int x = 0; x < vXSec.size(); x++)
    vXSec[x]->Reset("b");
  // tree
  TString dir = "Nominal";
  TString wSyst = "Nominal";
  if(configHelper->IsExpSys(configHelper->VarToSys(Var)))
    {
      dir = configHelper->GetVarDirTreePlain(Var, zMode);
      wSyst = zMode ? configHelper->GetVarWeight(Var) : "Nominal";
    }
  TreePlain* treeBgr = configHelper->GetTreeMCBgrPlain(Var, ch, dir, wSyst);
  if(configHelper->gTTBARDebug)
    printf("MC events BGR: %lld\n", treeBgr->TTreePtr->GetEntries());
  for(unsigned int x = 0; x < vXSec.size(); x++)
    vXSec[x]->InitTree(treeBgr, true);
  for(int e = 0; e < (MaxEvents ? MaxEvents : treeBgr->TTreePtr->GetEntries()); e++)
    {
      treeBgr->TTreePtr->GetEntry(e);
      for(unsigned int x = 0; x < vXSec.size(); x++)
	vXSec[x]->FillBgr(treeBgr);
    }
  for(unsigned int x = 0; x < vXSec.size(); x++)
    vXSec[x]->ClearTree(treeBgr, true);
  delete treeBgr;
}

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>> MC background ttbar >>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void PlotterDiffXSec::LoopMCBkgTTbarTreePlain(const TString& ch)
{
  // reset
  for(unsigned int x = 0; x < vXSec.size(); x++)
    vXSec[x]->Reset("o");
  // tree
  TString dir = configHelper->GetVarDirTreePlain(Var, zMode);
  TString wSyst = zMode ? configHelper->GetVarWeight(Var) : "Nominal";
  TreePlain* treeBtt = configHelper->GetTreeMCBttPlain(Var, ch, dir, wSyst, zMode);
  if(configHelper->gTTBARDebug)
    printf("MC events BGR ttbar: %lld\n", treeBtt->TTreePtr->GetEntries());
  for(unsigned int x = 0; x < vXSec.size(); x++)
    vXSec[x]->InitTree(treeBtt, true);
  //std::vector<double> vRew(vXSec.size());
  for(int e = 0; e < (MaxEvents ? MaxEvents : treeBtt->TTreePtr->GetEntries()); e++)
    {
      // 13.06.18 disabled reweighting (used for pT(t)): do not know how to reweight reco level based on gen level variable here
      //for(int x = 0; x < vXSec.size(); x++)
      //  vRew[x] = vXSec[x]->Reweight(treeBtt);
      treeBtt->TTreePtr->GetEntry(e);
      for(unsigned int x = 0; x < vXSec.size(); x++)
	//vXSec[x]->FillBtt(treeBtt, vRew[x]);
	vXSec[x]->FillBtt(treeBtt);
    }
  for(unsigned int x = 0; x < vXSec.size(); x++)
    vXSec[x]->ClearTree(treeBtt, true);
  delete treeBtt;
}

void PlotterDiffXSec::LoopAllTreesPlain(const TString& ch)
{
  if(!zFlagDataFilled || true)
    {
      LoopDataTreePlain(ch);
      zFlagDataFilled = 1;
    }
  if(!zFlagMCFilled || true)
    {
      LoopMCTreePlain(ch);
      zFlagMCFilled = 1;
    }
  if(!zFlagMCBkgTTbarFilled || true)
    {
      LoopMCBkgTTbarTreePlain(ch);
      zFlagMCBkgTTbarFilled = 1;
    }
  if(!zFlagMCBkgFilled || configHelper->IsExpSys(configHelper->VarToSys(zVarLast)) || configHelper->IsExpSys(configHelper->VarToSys(Var)) || true)
    {
      LoopMCBkgTreePlain(ch);
      zFlagMCBkgFilled = 1;
    }
  // store current variation
  zVarLast = Var;
  // print number of selected events in different samples
  if(vXSec.size())
    vXSec[0]->PrintRecWeights(ch, Var, ModeKinRecoTree);
}

std::vector<TH1D*> PlotterDiffXSec::GetNPCorrHistos(std::vector<UnfoldingIOHandler*> vXSec, const TString& inputSuffix, const TString weight)
{
  std::vector<TH1D*> vH(vXSec.size());
  for(size_t h = 0; h < vH.size(); h++)
    vH[h] = (TH1D*) vXSec[h]->HGenC()->Clone();
  TString fileIn = configHelper->gBaseDir + "../plainTree-oz/npcorr/emu_ttbarsignalplustau_" + inputSuffix + ".root";
  const TString ch = ""; // dummy
  StandaloneLoopMCTreeGenPlain(ch, 0, vH, vXSec, fileIn, weight);
  return vH;
}


void PlotterDiffXSec::Analyse(const TString& ch)
{
  current_channel = ch;
  if(Var == "BR_UP" || Var == "BR_DOWN")
    return;
  printf("Analyse()  var = %s  ch = %s\n", Var.Data(), ch.Data());

  //// pT(t) reweighting
  //if(Var == "REWTOPPT")
  //  {
  //    // https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting (accessed 01.05.18)
  //    TF1* f = new TF1("", "TMath::Sqrt(TMath::Exp(0.0615 - 0.0005 * x))", 0.0, 1000.0);
  //    ZVar* pt = new ZVarPtt();
  //    ZVar* pta = new ZVarPtat();
  //    for(UnfoldingIOHandler* unfoldingIOHandler : vXSec)
  //      {
  //        unfoldingIOHandler->AddVarReweight(pt, f);
  //        unfoldingIOHandler->AddVarReweight(pta, (TF1*)f->Clone());
  //      }
  //  }

  // gen only run for alternative theory MC samples
  if(configHelper->IsAltTheorSys(Var))
    {
      for(unsigned int x = 0; x < vXSec.size(); x++)
        {
          UnfoldingIOHandler* unfoldingIOHandler = vXSec[x];
          if(ModeKinRecoTree == 0)
            configHelper->outDir = configHelper->gAnalDir + "/" + ch + "/" + unfoldingIOHandler->Suffix + "/" + Var;
          else
            configHelper->outDir = configHelper->gAnalDir + "/" + ch + TString::Format("-kr%d/", ModeKinRecoTree) + unfoldingIOHandler->Suffix + "/" + Var;
          TH1D* hGen = unfoldingIOHandler->GetPlainHisto();
          StandaloneLoopMCTreeGenPlain(ch, 0, hGen, unfoldingIOHandler);
          //hGen->Print("all");
          std::vector<TH1D*> vh;
          vh.push_back(hGen);
          TString fileName = configHelper->outDir + "/gen.txt";
          PrintCP(vh, fileName);
        }
      return;
    }

  LoopAllTreesPlain(ch);
  //throw;

  //std::vector<ZUnfold*> vUnf = CreateZUnfoldVector(1);
  std::vector<ZUnfold*> vUnf = CreateZUnfoldVector(3);
  //std::vector<ZUnfold*> vUnf = CreateZUnfoldVector(4);

  for(unsigned int x = 0; x < vXSec.size(); x++)
    {
      UnfoldingIOHandler* unfoldingIOHandler = vXSec[x];

      if(unfoldingIOHandler->DoParticle)
	unfoldingIOHandler->AddOutOfSpaceToTTother();

      if(ModeKinRecoTree == 0)
	configHelper->outDir = configHelper->gAnalDir + "/" + ch + "/" + unfoldingIOHandler->Suffix + "/" + Var;
      else
	configHelper->outDir = configHelper->gAnalDir + "/" + ch + TString::Format("-kr%d/", ModeKinRecoTree) + unfoldingIOHandler->Suffix + "/" + Var;
      //ClearDir(gOutDir);

      // control plot
      if(unfoldingIOHandler->DoYields)
        {
          std::vector<TH1D*> vh;
          vh.push_back((TH1D*) unfoldingIOHandler->HDatC()->Clone());
          vh.push_back((TH1D*) unfoldingIOHandler->HRecC()->Clone());
          vh.back()->SetTitle("MCsig");
          vh.push_back((TH1D*) unfoldingIOHandler->HBttC()->Clone());
          vh.back()->SetTitle("MCbtt");
          vh.push_back((TH1D*) unfoldingIOHandler->HBgrC()->Clone());
          vh.back()->SetTitle("MCbkg");
          vh.push_back((TH1D*) unfoldingIOHandler->HGenC()->Clone());
          vh.back()->SetTitle("MCgen");
          PrintCP(vh, configHelper->outDir + "/cp.txt");
          //const std::vector<TH1D*> vhRead = ReadCP(outDir + "/cp.txt");

          if(!FlagSkipExtraPlots || std::find(VVarStoreCP.begin(), VVarStoreCP.end(), Var) != VVarStoreCP.end())
	    {
	      CPParams pars;
	      if(unfoldingIOHandler->RatMin < +1000.0)
		pars.RatMin = unfoldingIOHandler->RatMin;
	      if(unfoldingIOHandler->RatMax > -1000.0)
		pars.RatMax = unfoldingIOHandler->RatMax;
	      //ControlPlot(unfoldingIOHandler, vhRead, outDir + "/cp", pars);
	      ControlPlot(unfoldingIOHandler, configHelper->outDir + "/cp", pars);
	    }
        }

      // calculate PSE
      if(!FlagSkipExtraPlots || std::find(VVarStorePSE.begin(), VVarStorePSE.end(), Var) != VVarStorePSE.end())
        {
          if(unfoldingIOHandler->DoPSE)
	    {
	      CalculatePSE(unfoldingIOHandler, configHelper->outDir + "/pse");
	      DrawMigrationMatrix(unfoldingIOHandler, configHelper->outDir + "/mig");
	    }
        }
      if(unfoldingIOHandler->DoPSE) PrintMigrationMatrix(unfoldingIOHandler, configHelper->outDir + "/mig.txt");

      // calculate x-section
      if(unfoldingIOHandler->DoXSec)
        {
          XSecParams pars(unfoldingIOHandler);
          //pars.RatMin = 0.5;
          //pars.RatMax = 3.5;
          //pars.YMax = 0.0075;
          //vUnf[1]->Color = 0;
          const bool flagStoreExtraPlots = !FlagSkipExtraPlots || std::find(VVarStoreXSec.begin(), VVarStoreXSec.end(), Var) != VVarStoreXSec.end();
          unfoldingIOHandler->SubtrackBackground();
          DoXSec(unfoldingIOHandler, ch, vUnf, configHelper->outDir, pars, flagStoreExtraPlots);
        }
      else
        {
          // print rec level only for xsec with same reco and gen binning
          bool flagPrintRec = 1;
          for(int d = 0; d < unfoldingIOHandler->Dim(); d++)
	    {
	      ZVar* var = unfoldingIOHandler->Var(d);
	      if(var->Bins.size() != 2 || var->Bins[0].size() != var->Bins[1].size())
		{
		  flagPrintRec = 0;
		  break;
		}
	    }
          if(flagPrintRec)
	    {
	      unfoldingIOHandler->SubtrackBackground();
	      StoreRecUnf(unfoldingIOHandler, NULL, 0, configHelper->outDir);
	      StoreRecUnf(unfoldingIOHandler, NULL, 1, configHelper->outDir);
	    }
        }
      // prtint tex table with binning
      if(Var == "Nominal")
        {
          TString fileName = configHelper->outDir + "/bins.tex";
          PrintTexBinning(unfoldingIOHandler, fileName);
        }
    }

  for(auto alg : vUnf)
    delete alg;
  configHelper->outDir = "";
}

void PlotterDiffXSec::StoreAndPlotNPCorr(const std::vector<std::vector<std::vector<TH1D*> > > vvvH, const TString& dirNPCorr)
{
  const TString baseFileName = "npcorr";

  // dummy for storing
  //std::pair<TH1D*, TH1D*> hUDEmpty = std::pair<TH1D*, TH1D*>(NULL, NULL);
  const std::map<std::pair<TString, TString>, TH1D*> mapVarsEmpty;

  for(unsigned int x = 0; x < vXSec.size(); x++)
    {
      // histograms to plot
      std::vector<TH1D*> vHPlot(2);
      // 1
      TH1D* hUnity = utils::ttmd::CloneTH1D(vvvH[0][0][x], 1.0);
      vHPlot[0] = hUnity;
      vHPlot[0]->SetTitle("");
      // full
      vHPlot[1] = (TH1D*) vvvH[0][0][x]->Clone();
      // TODO vvvH[1][4][x] includes finite top width effect: but should it be included?
      //vHPlot[1]->Divide(  vvvH[1][4][x]); // with finite width
      vHPlot[1]->Divide(  vvvH[1][3][x]);
      vHPlot[1]->SetTitle("NP full");

      // make plot
      TString fileName = dirNPCorr + "/" + vXSec[x]->Suffix + "/" + baseFileName;
      PlotVars(vXSec[x], vHPlot, vHPlot, fileName, true, false);

      // store txt file
      fileName += ".dat";
      std::vector<TH1D*> vHStore(vHPlot.begin() + 1, vHPlot.end());
      // ZIOxFitter::PrintXFitterDat(vXSec[x], vHStore, hUDEmpty, mapVarsEmpty, fileName, true);
      printf("Stored %s\n", fileName.Data());
    }
}

void PlotterDiffXSec::StoreAndPlotNPCorrAll4(const std::vector<std::vector<std::vector<TH1D*> > > vvvH, const TString& dirNPCorr)
    {
      const TString baseFileName = "npcorr4";

      // dummy for storing
      //std::pair<TH1D*, TH1D*> hUDEmpty = std::pair<TH1D*, TH1D*>(NULL, NULL);
      const std::map<std::pair<TString, TString>, TH1D*> mapVarsEmpty;

      for(unsigned int x = 0; x < vXSec.size(); x++)
      {
        // histograms to plot
        std::vector<TH1D*> vHPlot(5);
        // 1
        TH1D* hUnity = utils::ttmd::CloneTH1D(vvvH[0][0][x], 1.0);
        vHPlot[0] = hUnity;
        vHPlot[0]->SetTitle("");
        // full
        vHPlot[1] = (TH1D*) vvvH[0][0][x]->Clone();
        vHPlot[1]->Divide(  vvvH[1][3][x]);
        vHPlot[1]->SetTitle("NP full");
        // had
        vHPlot[2] = (TH1D*) vvvH[0][0][x]->Clone();
        vHPlot[2]->Divide(  vvvH[0][2][x]);
        vHPlot[2]->SetTitle("NP had.");
        // MPI
        vHPlot[3] = (TH1D*) vvvH[0][0][x]->Clone();
        vHPlot[3]->Divide(  vvvH[0][1][x]);
        vHPlot[3]->SetTitle("NP MPI");
        // definition
        vHPlot[4] = (TH1D*) vvvH[0][3][x]->Clone();
        vHPlot[4]->Divide(  vvvH[1][3][x]);
        vHPlot[4]->SetTitle("NP def.");

        // make plot
        TString fileName = dirNPCorr + "/" + vXSec[x]->Suffix + "/" + baseFileName;
        PlotVars(vXSec[x], vHPlot, vHPlot, fileName, true, false);

        // store txt file
        fileName += ".dat";
        std::vector<TH1D*> vHStore(vHPlot.begin() + 1, vHPlot.end());
        // ZIOxFitter::PrintXFitterDat(vXSec[x], vHStore, hUDEmpty, mapVarsEmpty, fileName, true);
        printf("Stored %s\n", fileName.Data());
      }
    }

void PlotterDiffXSec::StoreAndPlotNPCorrHad2(const std::vector<std::vector<std::vector<TH1D*> > > vvvH, const TString& dirNPCorr)
    {
      const TString baseFileName = "npcorr-had2";

      // dummy for storing
      //std::pair<TH1D*, TH1D*> hUDEmpty = std::pair<TH1D*, TH1D*>(NULL, NULL);
      const std::map<std::pair<TString, TString>, TH1D*> mapVarsEmpty;

      for(unsigned int x = 0; x < vXSec.size(); x++)
      {
        // histograms to plot
        std::vector<TH1D*> vHPlot(3);
        // 1
        TH1D* hUnity = utils::ttmd::CloneTH1D(vvvH[0][0][x], 1.0);
        vHPlot[0] = hUnity;
        vHPlot[0]->SetTitle("");
        // had (1)
        vHPlot[1] = (TH1D*) vvvH[0][0][x]->Clone();
        vHPlot[1]->Divide(  vvvH[0][2][x]);
        vHPlot[1]->SetTitle("had. def-nom/def-nohad");

        // had (2)
        vHPlot[2] = (TH1D*) vvvH[0][1][x]->Clone();
        vHPlot[2]->Divide(  vvvH[0][3][x]);
        vHPlot[2]->SetTitle("had. def-nompi/def-nonp");

        // make plot
        TString fileName = dirNPCorr + "/" + vXSec[x]->Suffix + "/" + baseFileName;
        PlotVars(vXSec[x], vHPlot, vHPlot, fileName, true, false);

        // store txt file
        fileName += ".dat";
        std::vector<TH1D*> vHStore(vHPlot.begin() + 1, vHPlot.end());
        // ZIOxFitter::PrintXFitterDat(vXSec[x], vHStore, hUDEmpty, mapVarsEmpty, fileName, true);
        printf("Stored %s\n", fileName.Data());
      }
    }

void PlotterDiffXSec::StoreAndPlotNPCorrMPI3(const std::vector<std::vector<std::vector<TH1D*> > > vvvH, const TString& dirNPCorr)
    {
      const TString baseFileName = "npcorr-mpi3";

      // dummy for storing
      //std::pair<TH1D*, TH1D*> hUDEmpty = std::pair<TH1D*, TH1D*>(NULL, NULL);
      const std::map<std::pair<TString, TString>, TH1D*> mapVarsEmpty;

      for(unsigned int x = 0; x < vXSec.size(); x++)
      {
        // histograms to plot
        std::vector<TH1D*> vHPlot(4);
        // 1
        TH1D* hUnity = utils::ttmd::CloneTH1D(vvvH[0][0][x], 1.0);
        vHPlot[0] = hUnity;
        vHPlot[0]->SetTitle("");
        // MPI(1)
        vHPlot[1] = (TH1D*) vvvH[0][0][x]->Clone();
        vHPlot[1]->Divide(  vvvH[0][1][x]);
        vHPlot[1]->SetTitle("MPI def-nom/def-nompi");
        // MPI(2)
        vHPlot[2] = (TH1D*) vvvH[0][2][x]->Clone();
        vHPlot[2]->Divide(  vvvH[0][3][x]);
        vHPlot[2]->SetTitle("MPI def-nohad/def-nonp");
        // MPI(3)
        vHPlot[3] = (TH1D*) vvvH[1][2][x]->Clone();
        vHPlot[3]->Divide(  vvvH[1][3][x]);
        vHPlot[3]->SetTitle("MPI parnt-nohad/parnt-nonp");

        // make plot
        TString fileName = dirNPCorr + "/" + vXSec[x]->Suffix + "/" + baseFileName;
        PlotVars(vXSec[x], vHPlot, vHPlot, fileName, true, false);

        // store txt file
        fileName += ".dat";
        std::vector<TH1D*> vHStore(vHPlot.begin() + 1, vHPlot.end());
        // ZIOxFitter::PrintXFitterDat(vXSec[x], vHStore, hUDEmpty, mapVarsEmpty, fileName, true);
        printf("Stored %s\n", fileName.Data());
      }
    }

void PlotterDiffXSec::StoreAndPlotNPCorrDef2(const std::vector<std::vector<std::vector<TH1D*> > > vvvH, const TString& dirNPCorr)
    {
      const TString baseFileName = "npcorr-def2";

      // dummy for storing
      //std::pair<TH1D*, TH1D*> hUDEmpty = std::pair<TH1D*, TH1D*>(NULL, NULL);
      const std::map<std::pair<TString, TString>, TH1D*> mapVarsEmpty;

      for(unsigned int x = 0; x < vXSec.size(); x++)
      {
        // histograms to plot
        std::vector<TH1D*> vHPlot(3);
        // 1
        TH1D* hUnity = utils::ttmd::CloneTH1D(vvvH[0][0][x], 1.0);
        vHPlot[0] = hUnity;
        vHPlot[0]->SetTitle("");
        // definition(1)
        vHPlot[1] = (TH1D*) vvvH[0][3][x]->Clone();
        vHPlot[1]->Divide(  vvvH[1][3][x]);
        vHPlot[1]->SetTitle("def. def-nonp/parnt-nonp");
        // definition(2)
        vHPlot[2] = (TH1D*) vvvH[0][2][x]->Clone();
        vHPlot[2]->Divide(  vvvH[1][2][x]);
        vHPlot[2]->SetTitle("def. def-nohad/parnt-nohad");

        // make plot
        TString fileName = dirNPCorr + "/" + vXSec[x]->Suffix + "/" + baseFileName;
        PlotVars(vXSec[x], vHPlot, vHPlot, fileName, true, false);

        // store txt file
        fileName += ".dat";
        std::vector<TH1D*> vHStore(vHPlot.begin() + 1, vHPlot.end());
        // ZIOxFitter::PrintXFitterDat(vXSec[x], vHStore, hUDEmpty, mapVarsEmpty, fileName, true);
        printf("Stored %s\n", fileName.Data());
      }
    }

void PlotterDiffXSec::StoreAndPlotNPCorrCheck1(const std::vector<std::vector<std::vector<TH1D*> > > vvvH, const TString& dirNPCorr)
    {
      const TString baseFileName = "npcorr-check1";

      // dummy for storing
      //std::pair<TH1D*, TH1D*> hUDEmpty = std::pair<TH1D*, TH1D*>(NULL, NULL);
      const std::map<std::pair<TString, TString>, TH1D*> mapVarsEmpty;

      for(unsigned int x = 0; x < vXSec.size(); x++)
      {
        // histograms to plot
        std::vector<TH1D*> vHPlot(2);
        // 1
        TH1D* hUnity = utils::ttmd::CloneTH1D(vvvH[0][0][x], 1.0);
        vHPlot[0] = hUnity;
        vHPlot[0]->SetTitle("");
        // check
        vHPlot[1] = (TH1D*) vvvH[1][3][x]->Clone();
        vHPlot[1]->Divide(  vvvH[1][4][x]);
        vHPlot[1]->SetTitle("NP check parnt-nonp/parnt-nonpd");

        // make plot
        TString fileName = dirNPCorr + "/" + vXSec[x]->Suffix + "/" + baseFileName;
        PlotVars(vXSec[x], vHPlot, vHPlot, fileName, true, false);

        // store txt file
        fileName += ".dat";
        std::vector<TH1D*> vHStore(vHPlot.begin() + 1, vHPlot.end());
        // ZIOxFitter::PrintXFitterDat(vXSec[x], vHStore, hUDEmpty, mapVarsEmpty, fileName, true);
        printf("Stored %s\n", fileName.Data());
      }
    }

void PlotterDiffXSec::StoreAndPlotNPCorrVars(const std::vector<std::vector<std::vector<TH1D*> > > vvvH, const TString& dirNPCorr)
    {
      const TString baseFileName = "npcorr-vars";

      // dummy for storing
      //std::pair<TH1D*, TH1D*> hUDEmpty = std::pair<TH1D*, TH1D*>(NULL, NULL);
      const std::map<std::pair<TString, TString>, TH1D*> mapVarsEmpty;

      for(unsigned int x = 0; x < vXSec.size(); x++)
      {
        // histograms to plot
        std::vector<TH1D*> vHPlot(7);
        // 1
        TH1D* hUnity = utils::ttmd::CloneTH1D(vvvH[0][0][x], 1.0);
        vHPlot[0] = hUnity;
        vHPlot[0]->SetTitle("");
        // nom
        vHPlot[1] = (TH1D*) vvvH[0][0][x]->Clone();
        //vHPlot[1]->Divide(  vvvH[1][4][x]);
        vHPlot[1]->Divide(  vvvH[1][3][x]);
        vHPlot[1]->SetTitle("NP nom");
        // UETUNE_DOWN
        vHPlot[2] = (TH1D*) vvvH[0][5][x]->Clone();
        vHPlot[2]->Scale(vvvH[0][0][x]->Integral() / vvvH[0][5][x]->Integral());
        //vHPlot[2]->Divide(  vvvH[1][4][x]);
        vHPlot[2]->Divide(  vvvH[1][3][x]);
        vHPlot[2]->SetTitle("NP UETUNE_DOWN");
        // UETUNE_UP
        vHPlot[3] = (TH1D*) vvvH[0][6][x]->Clone();
        vHPlot[3]->Scale(vvvH[0][0][x]->Integral() / vvvH[0][6][x]->Integral());
        //vHPlot[3]->Divide(  vvvH[1][4][x]);
        vHPlot[3]->Divide(  vvvH[1][3][x]);
        vHPlot[3]->SetTitle("NP UETUNE_UP");
        // ERDON
        vHPlot[4] = (TH1D*) vvvH[0][7][x]->Clone();
        vHPlot[4]->Scale(vvvH[0][0][x]->Integral() / vvvH[0][7][x]->Integral());
        //vHPlot[4]->Divide(  vvvH[1][4][x]);
        vHPlot[4]->Divide(  vvvH[1][3][x]);
        vHPlot[4]->SetTitle("NP ERDON");
        // ERDONRETUNE
        vHPlot[5] = (TH1D*) vvvH[0][8][x]->Clone();
        vHPlot[5]->Scale(vvvH[0][0][x]->Integral() / vvvH[0][8][x]->Integral());
        //vHPlot[5]->Divide(  vvvH[1][4][x]);
        vHPlot[5]->Divide(  vvvH[1][3][x]);
        vHPlot[5]->SetTitle("NP ERDONRETUNE");
        // GLUONMOVETUNE
        vHPlot[6] = (TH1D*) vvvH[0][9][x]->Clone();
        vHPlot[6]->Scale(vvvH[0][0][x]->Integral() / vvvH[0][9][x]->Integral());
        //vHPlot[6]->Divide(  vvvH[1][4][x]);
        vHPlot[6]->Divide(  vvvH[1][3][x]);
        vHPlot[6]->SetTitle("NP GLUONMOVETUNE");

        // make plot
        TString fileName = dirNPCorr + "/" + vXSec[x]->Suffix + "/" + baseFileName;
        PlotVars(vXSec[x], vHPlot, vHPlot, fileName, true, false);

        // store txt file
        fileName += ".dat";
        std::vector<TH1D*> vHStore(vHPlot.begin() + 1, vHPlot.end());
        // ZIOxFitter::PrintXFitterDat(vXSec[x], vHStore, hUDEmpty, mapVarsEmpty, fileName, true);
        printf("Stored %s\n", fileName.Data());
      }
    }

void PlotterDiffXSec::StoreAndPlotNPCorrVarsPS(const std::vector<std::vector<std::vector<TH1D*> > > vvvH, const TString& dirNPCorr)
    {
      const TString baseFileName = "npcorr-vars-ps";

      // dummy for storing
      //std::pair<TH1D*, TH1D*> hUDEmpty = std::pair<TH1D*, TH1D*>(NULL, NULL);
      const std::map<std::pair<TString, TString>, TH1D*> mapVarsEmpty;

      for(unsigned int x = 0; x < vXSec.size(); x++)
      {
        // histograms to plot
        std::vector<TH1D*> vHPlot(8);
        // 1
        TH1D* hUnity = utils::ttmd::CloneTH1D(vvvH[0][0][x], 1.0);
        vHPlot[0] = hUnity;
        vHPlot[0]->SetTitle("");
        // nom
        vHPlot[1] = (TH1D*) vvvH[0][0][x]->Clone();
        //vHPlot[1]->Divide(  vvvH[1][4][x]);
        vHPlot[1]->Divide(  vvvH[1][3][x]);
        vHPlot[1]->SetTitle("NP nom");
        // PSFSRSCALE_UP
        int nVar = 10;
        int offset = 6;
        vHPlot[2] = (TH1D*) vvvH[0][nVar][x]->Clone();
        //vHPlot[2]->Scale(vvvH[0][0][x]->Integral() / vvvH[0][nVar][x]->Integral());
        vHPlot[2]->Divide(  vvvH[1][nVar + offset][x]);
        //vHPlot[2]->Scale(vvvH[1][nVar + offset][x]->Integral() / vvvH[1][3][x]->Integral());
        utils::ttmd::ScaleRelativeHistoVar(vHPlot[2], 1.0 / TMath::Sqrt(2.0), vHPlot[1]);
        vHPlot[2]->SetTitle("NP PSFSRSCALE_UP");
        // PSFSRSCALE_DOWN
        nVar = 11;
        vHPlot[3] = (TH1D*) vvvH[0][nVar][x]->Clone();
        //vHPlot[3]->Scale(vvvH[0][0][x]->Integral() / vvvH[0][nVar][x]->Integral());
        vHPlot[3]->Divide(  vvvH[1][nVar + offset][x]);
        utils::ttmd::ScaleRelativeHistoVar(vHPlot[3], 1.0 / TMath::Sqrt(2.0), vHPlot[1]);
        //vHPlot[3]->Scale(vvvH[1][nVar + offset][x]->Integral() / vvvH[1][3][x]->Integral());
        vHPlot[3]->SetTitle("NP PSFSRSCALE_DOWN");
        // PSISRSCALE_UP
        nVar = 12;
        vHPlot[4] = (TH1D*) vvvH[0][nVar][x]->Clone();
        //vHPlot[4]->Scale(vvvH[0][0][x]->Integral() / vvvH[0][nVar][x]->Integral());
        vHPlot[4]->Divide(  vvvH[1][nVar + offset][x]);
        //vHPlot[4]->Scale(vvvH[1][nVar + offset][x]->Integral() / vvvH[1][3][x]->Integral());
        vHPlot[4]->SetTitle("NP PSISRSCALE_UP");
        // PSISRSCALE_DOWN
        nVar = 13;
        vHPlot[5] = (TH1D*) vvvH[0][nVar][x]->Clone();
        //vHPlot[5]->Scale(vvvH[0][0][x]->Integral() / vvvH[0][nVar][x]->Integral());
        vHPlot[5]->Divide(  vvvH[1][nVar + offset][x]);
        //vHPlot[5]->Scale(vvvH[1][nVar + offset][x]->Integral() / vvvH[1][3][x]->Integral());
        vHPlot[5]->SetTitle("NP PSISRSCALE_DOWN");
        // MATCH_UP
        nVar = 14;
        vHPlot[6] = (TH1D*) vvvH[0][nVar][x]->Clone();
        //vHPlot[6]->Scale(vvvH[0][0][x]->Integral() / vvvH[0][nVar][x]->Integral());
        vHPlot[6]->Divide(  vvvH[1][nVar + offset][x]);
        //vHPlot[6]->Scale(vvvH[1][nVar + offset][x]->Integral() / vvvH[1][3][x]->Integral());
        vHPlot[6]->SetTitle("NP MATCH_UP");
        // MATCH_DOWN
        nVar = 15;
        vHPlot[7] = (TH1D*) vvvH[0][nVar][x]->Clone();
        //vHPlot[7]->Scale(vvvH[0][0][x]->Integral() / vvvH[0][nVar][x]->Integral());
        vHPlot[7]->Divide(  vvvH[1][nVar + offset][x]);
        //vHPlot[7]->Scale(vvvH[1][nVar + offset][x]->Integral() / vvvH[1][3][x]->Integral());
        vHPlot[7]->SetTitle("NP MATCH_DOWN");

        // make plot
        TString fileName = dirNPCorr + "/" + vXSec[x]->Suffix + "/" + baseFileName;
        PlotVars(vXSec[x], vHPlot, vHPlot, fileName, true, false);

        // store txt file
        fileName += ".dat";
        std::vector<TH1D*> vHStore(vHPlot.begin() + 1, vHPlot.end());
        // ZIOxFitter::PrintXFitterDat(vXSec[x], vHStore, hUDEmpty, mapVarsEmpty, fileName, true);
        printf("Stored %s\n", fileName.Data());
      }
    }

void PlotterDiffXSec::StoreAndPlotNPCorrVarsME(const std::vector<std::vector<std::vector<TH1D*> > > vvvH, const TString& dirNPCorr)
    {
      const TString baseFileName = "npcorr-vars-me";

      // dummy for storing
      //std::pair<TH1D*, TH1D*> hUDEmpty = std::pair<TH1D*, TH1D*>(NULL, NULL);
      const std::map<std::pair<TString, TString>, TH1D*> mapVarsEmpty;

      for(unsigned int x = 0; x < vXSec.size(); x++)
      {
        // histograms to plot
        std::vector<TH1D*> vHPlot(8);
        // 1
        TH1D* hUnity = utils::ttmd::CloneTH1D(vvvH[0][0][x], 1.0);
        vHPlot[0] = hUnity;
        vHPlot[0]->SetTitle("");
        // nom
        vHPlot[1] = (TH1D*) vvvH[0][0][x]->Clone();
        //vHPlot[1]->Divide(  vvvH[1][4][x]);
        vHPlot[1]->Divide(  vvvH[1][3][x]);
        vHPlot[1]->SetTitle("NP nom");
        // MERENSCALE_UP
        int nVar = 22;
        int offset = 6;
        vHPlot[2] = (TH1D*) vvvH[0][nVar][x]->Clone();
        vHPlot[2]->Divide(  vvvH[1][nVar + offset][x]);
        utils::ttmd::ScaleRelativeHistoVar(vHPlot[2], 1.0 / TMath::Sqrt(2.0), vHPlot[1]);
        vHPlot[2]->SetTitle("NP MERENSCALE_UP");
        // MERENSCALE_DOWN
        nVar = 23;
        vHPlot[3] = (TH1D*) vvvH[0][nVar][x]->Clone();
        vHPlot[3]->Divide(  vvvH[1][nVar + offset][x]);
        utils::ttmd::ScaleRelativeHistoVar(vHPlot[3], 1.0 / TMath::Sqrt(2.0), vHPlot[1]);
        vHPlot[3]->SetTitle("NP MERENSCALE_DOWN");
        // MEFACSCALE_UP
        nVar = 24;
        vHPlot[4] = (TH1D*) vvvH[0][nVar][x]->Clone();
        vHPlot[4]->Divide(  vvvH[1][nVar + offset][x]);
        vHPlot[4]->SetTitle("NP MEFACSCALE_UP");
        // MEFACSCALE_DOWN
        nVar = 25;
        vHPlot[5] = (TH1D*) vvvH[0][nVar][x]->Clone();
        vHPlot[5]->Divide(  vvvH[1][nVar + offset][x]);
        vHPlot[5]->SetTitle("NP MEFACSCALE_DOWN");
        // MESCALE_UP
        nVar = 26;
        vHPlot[6] = (TH1D*) vvvH[0][nVar][x]->Clone();
        vHPlot[6]->Divide(  vvvH[1][nVar + offset][x]);
        vHPlot[6]->SetTitle("NP MESCALE_UP");
        // MESCALE_DOWN
        nVar = 27;
        vHPlot[7] = (TH1D*) vvvH[0][nVar][x]->Clone();
        vHPlot[7]->Divide(  vvvH[1][nVar + offset][x]);
        vHPlot[7]->SetTitle("NP MESCALE_DOWN");

        // make plot
        TString fileName = dirNPCorr + "/" + vXSec[x]->Suffix + "/" + baseFileName;
        PlotVars(vXSec[x], vHPlot, vHPlot, fileName, true, false);

        // store txt file
        fileName += ".dat";
        std::vector<TH1D*> vHStore(vHPlot.begin() + 1, vHPlot.end());
        // ZIOxFitter::PrintXFitterDat(vXSec[x], vHStore, hUDEmpty, mapVarsEmpty, fileName, true);
        printf("Stored %s\n", fileName.Data());
      }
    }

void PlotterDiffXSec::CalculateNPCorr()
    {
      printf("CalculateNPCorr()  var = %s\n", Var.Data());
      //this->MaxEvents = 50000;

      //TString ch = "";
      const TString dirNPCorr = configHelper->gAnalDir + "/" + configHelper->gNPCorrSubdir;
      //std::pair<TH1D*, TH1D*> hUDEmpty = std::pair<TH1D*, TH1D*>(NULL, NULL);
      //const std::map<std::pair<TString, TString>, TH1D*> mapVarsEmpty;

      // default jet definition (particle level isolated from b,l)
      std::vector<UnfoldingIOHandler*>& vXSecDef = vXSec;

      // parton no top definition
      std::vector<UnfoldingIOHandler*> vXSecParnt = vXSec;
      for(unsigned int x = 0; x < vXSec.size(); x++)
      {
        vXSecParnt[x] = new UnfoldingIOHandler(*vXSecDef[x]);
        for(ZVar* var : vXSecParnt[x]->Vars())
        {
          //var->Expression.ReplaceAll("njdef", "njparnt");
          //var->ExpressionCut.ReplaceAll("njdef", "njparnt");
          var->Expression.ReplaceAll("def", "parnt");
          var->ExpressionCut.ReplaceAll("def", "parnt");
        }
        vXSecParnt[x]->Init();
      }

      // some notations
      // jet definition
      //   [0] def: default
      //   [1] parnt: parton no top decay products
      // MC sample
      //   [0] nom: nominal
      //   [1] nompi: no MPI
      //   [2] nohad: no hadronisation
      //   [3] nonp: no hadronisation & no MPI
      //   [4] nonpd: stable top & no hadronisation & no MPI
      //   [5] nom_uetuneup: variation
      std::vector<std::vector<std::vector<TH1D*> > > vvvH = utils::ttmd::CreateVectorPtr3D<TH1D>(2, 5 + 5, vXSecDef.size());
      std::vector<TString> vMCName = { "nom", "nompi", "nohad", "nonp", "nonpd",
                                       "nom_uetunedown", "nom_uetuneup", "nom_erdon", "nom_erdonretune", "nom_gluonmovetune" };

      // 2.09.18 PS variations
      bool flagPSvars = 1;
      if(flagPSvars)
      {
        vvvH = utils::ttmd::CreateVectorPtr3D<TH1D>(2, 5 + 5 + 12 + 12, vXSecDef.size());
        vMCName = { "nom", "nompi", "nohad", "nonp", "nonpd",
                                               "nom_uetunedown", "nom_uetuneup", "nom_erdon", "nom_erdonretune", "nom_gluonmovetune",
                  "nom_psfsrscaleup", "nom_psfsrscaledown", "nom_psisrscaleup", "nom_psisrscaledown", "nom_matchup", "nom_matchdown",
                  "nonp_psfsrscaleup", "nonp_psfsrscaledown", "nonp_psisrscaleup", "nonp_psisrscaledown", "nonp_matchup", "nonp_matchdown",
                  "nom_MERENSCALE_DOWN", "nom_MERENSCALE_UP", "nom_MEFACSCALE_DOWN", "nom_MEFACSCALE_UP", "nom_MESCALE_DOWN", "nom_MESCALE_UP",
        "nonp_MERENSCALE_DOWN", "nonp_MERENSCALE_UP", "nonp_MEFACSCALE_DOWN", "nonp_MEFACSCALE_UP", "nonp_MESCALE_DOWN", "nonp_MESCALE_UP"};
      }

      for(size_t n = 0; n < vMCName.size(); n++)
      {
        //if(n != 0 && n != 4 && n < 5) continue;
        TString nameBase = vMCName[n];
        TString weight = "";
        if(vMCName[n].Contains("SCALE"))
        {
          // "SCALE" in name is present for ME scale vars
          // they could start only with nom_ or _nonp
          weight = vMCName[n];
          if(weight.BeginsWith("nom_"))
          {
            nameBase = "nom_meweights";
            weight.ReplaceAll("nom_", "");
          }
          else if(weight.BeginsWith("nonp_"))
          {
            nameBase = "nonp_meweights";
            weight.ReplaceAll("nonp_", "");
          }
          else
            assert(0);
        }
        if(vMCName[n] != "nonpd" && !vMCName[n].BeginsWith("nonp_"))
          vvvH[0][n] = GetNPCorrHistos(vXSecDef, nameBase, weight);
        if(vMCName[n] != "nompi" && vMCName[n] != "nom" && !vMCName[n].BeginsWith("nom_"))
          vvvH[1][n] = GetNPCorrHistos(vXSecParnt, nameBase, weight);
        /*if(vMCName[n] == "nonp" || vMCName[n] == "nonpd")
        {
          printf("name: %s\n", vMCName[n].Data());
          if(vvvH[0][n][0])
            vvvH[0][n][0]->Print("all");
          if(vvvH[1][n][0])
            vvvH[1][n][0]->Print("all");
        }*/
      }

      // normalised input histograms
      // TODO this way correctins are valid only for normalised cross sections
      std::vector<std::vector<std::vector<TH1D*> > > vvvHNorm = utils::ttmd::CreateVectorPtr3D<TH1D>(vvvH.size(), vvvH[0].size(), vvvH[0][0].size());
      for(size_t i = 0; i < vvvH.size(); i++)
      {
        for(size_t ii = 0; ii < vvvH[0].size(); ii++)
        {
          for(size_t iii = 0; iii < vvvH[0][0].size(); iii++)
          {
            if(vvvH[i][ii][iii])
            {
              vvvHNorm[i][ii][iii] = new TH1D(*vvvH[i][ii][iii]);
              vvvHNorm[i][ii][iii]->Scale(1.0 / vvvHNorm[i][ii][iii]->Integral());
            }
          }
        }
      }

      // plots and txt output
      StoreAndPlotNPCorr(vvvH, dirNPCorr);
      StoreAndPlotNPCorrAll4(vvvH, dirNPCorr);
      StoreAndPlotNPCorrHad2(vvvH, dirNPCorr);
      StoreAndPlotNPCorrMPI3(vvvH, dirNPCorr);
      StoreAndPlotNPCorrDef2(vvvH, dirNPCorr);
      StoreAndPlotNPCorrCheck1(vvvH, dirNPCorr);
      StoreAndPlotNPCorrVars(vvvH, dirNPCorr);
      if(flagPSvars)
      {
        StoreAndPlotNPCorrVarsPS(vvvHNorm, dirNPCorr);
        // ME scale vars
        StoreAndPlotNPCorrVarsME(vvvHNorm, dirNPCorr);
      }
    }

void PlotterDiffXSec::DoXSec(UnfoldingIOHandler* unfoldingIOHandler, const TString& ch, const std::vector<ZUnfold*> vAlg, const TString& dirOut, XSecParams pars,
                const bool flagStoreExtraPlots)
    {
      printf("DoXSec()  unfoldingIOHandler = %s  var = %s  ch = %s\n", unfoldingIOHandler->Suffix.Data(), Var.Data(), ch.Data());

      int nbins = unfoldingIOHandler->HGenC()->GetNbinsX();
      std::vector<TH1D*> vDtAbs;
      std::vector<TH1D*> vDtNor;
      std::vector<TH1D*> vThAbs;
      std::vector<TH1D*> vThNor;

      // for table with unfolding chi2, dof and tau
      TString strTab1 = "";
      TString strTab2 = "";
      for(unsigned int a = 0; a < vAlg.size(); a++)
      {
        // 1.10.18 detector level distributions only
        if(unfoldingIOHandler->DoXSec > 1)
        {
          StoreRecUnf(unfoldingIOHandler, NULL, 0, dirOut);
          StoreRecUnf(unfoldingIOHandler, NULL, 1, dirOut);
          break;
        }

        // unfold data
        // supress low stat. warnings for ee, mm channels
        vAlg[a]->FlagTreatStatisticsSeriously = false;
        if(ch == "ee" || ch == "mumu")
          vAlg[a]->FlagTreatStatisticsSeriously = false;
        TString fileNameBase;
        if(flagStoreExtraPlots)
          fileNameBase = dirOut + "/" + vAlg[a]->Suffix + "-";
        // do not unfold total cross section, only BBB
        if(unfoldingIOHandler->HGen(0)->GetNbinsX() == 1 && !vAlg[a]->Remark.Contains("BBB"))
          continue;
        vAlg[a]->RunUnfolding(unfoldingIOHandler, ch, fileNameBase);
        if(flagStoreExtraPlots)
        {
          if(!vAlg[a]->Remark.Contains("BBB"))
          {
            vAlg[a]->PlotInput(dirOut + "/tunf-" + vAlg[a]->Suffix);
          }
        }

        // 20.07.18 table for AN
        if(vAlg[a]->regChoice == ZUnfold::RegNo)
          strTab1 = TString::Format(" %d & %.1f &", vAlg[a]->Dof, vAlg[a]->Chi2Data);
        else if(vAlg[a]->regChoice == ZUnfold::RegRhoAvg)
          strTab2 = TString::Format(" %.1f + %.1f & %.2e \\\\", vAlg[a]->Chi2Data - vAlg[a]->Chi2Reg, vAlg[a]->Chi2Reg, vAlg[a]->TauBest);

        // 15.07.2018 save unfolding input for reco level studies
        // do this only for reg. unfolding
        if(vAlg[a]->regChoice == ZUnfold::RegNo)
        //if(vAlg[a]->regChoice != ZUnfold::RegNo)
        {
          StoreRecUnf(unfoldingIOHandler, vAlg[a], 0, dirOut);
          StoreRecUnf(unfoldingIOHandler, vAlg[a], 1, dirOut);
        }

        // MC (once per x-section)
        if(a == 0)
        {
          vThAbs.push_back((TH1D*) vAlg[a]->hGen->Clone());
          vThAbs.back()->SetTitle("POWHEG+PYTHIA");
          vThAbs.back()->SetLineColor(kBlack);
          vThNor.push_back((TH1D*) vThAbs.back()->Clone());
          vThNor.back()->Scale(1.0 / vThAbs.back()->Integral());
        }

        // push and print absolute x-section
        vDtAbs.push_back(vAlg[a]->hUnf);
        vDtAbs.back()->SetTitle(TString::Format("%s [%.0f/%d]", vAlg[a]->Suffix.Data(), vAlg[a]->Chi2, nbins));
        PrintDataXsec(vDtAbs.back(), vThAbs[0], dirOut + "/xsec-abs-" + vAlg[a]->Suffix + "-vec.txt");
        PrintDataCorMatrix(vAlg[a]->hCormatTotal, dirOut + "/xsec-abs-" + vAlg[a]->Suffix + "-cor.txt");

        // push and print normalised x-section
        vDtNor.push_back(vAlg[a]->hUnfN);
        vDtNor.back()->SetTitle(TString::Format("%s [%.0f/%d]", vAlg[a]->Suffix.Data(), vAlg[a]->Chi2N, nbins - 1));
        PrintDataXsec(vDtNor.back(), vThNor[0], dirOut + "/xsec-nor-" + vAlg[a]->Suffix + "-vec.txt");
        PrintDataCorMatrix(vAlg[a]->hCormatTotalN, dirOut + "/xsec-nor-" + vAlg[a]->Suffix + "-cor.txt");

        //ReadDataXsecAndCov(dirOut + "/xsec-nor-" + vAlg[a]->Suffix);
        //ReadDataXsecAndCov(dirOut + "/xsec-abs-" + vAlg[a]->Suffix);

      }

      // store for AN table
      TString fileName = dirOut + "/datatest.txt";
      utils::ttmd::CheckFile(fileName);
      FILE* f = fopen(fileName, "w");
      fprintf(f, "%s%s\n", strTab1.Data(), strTab2.Data());
      fclose(f);
      printf("Stored %s\n", fileName.Data());

      // plot x-sections
      if(unfoldingIOHandler->DoXSec == 1 || unfoldingIOHandler->DoXSec == 3)
      {
        if(flagStoreExtraPlots)
        {
          XSec(unfoldingIOHandler, vDtAbs, vThAbs, dirOut + "/xsec-abs", pars);
          XSec(unfoldingIOHandler, vDtNor, vThNor, dirOut + "/xsec-nor", pars);
        }

        // print summary x-sections tables
        PrintXsecAllAlg(vAlg, dirOut + "/xsec-all.txt");
      }
    }

void PlotterDiffXSec::StoreRecUnf(const UnfoldingIOHandler* unfoldingIOHandler, const ZUnfold* unf, const bool flagNorm, const TString& dirOut)
    {
      TString flagNormStr = flagNorm ? "-nor" : "-abs";

      const TH1D* hGen = new TH1D(*unfoldingIOHandler->HGen(0));
      TH2D* hRspNorm = DrawMigrationMatrix(unfoldingIOHandler);
      //TH1D* hInput = new TH1D(*unf->hInput);
      TH1D* hInput = (TH1D*) unfoldingIOHandler->HDnbC()->Clone();
      TH2D* hInputNCov = utils::ttmd::MakeCovTH2(hInput);
      if(flagNorm)
      {
        utils::ttmd::Normalise(hInput, hInputNCov);
        //hGen->Scale(1.0 / hGen->Integral());
      }

      //printf("StoreRecUnf **** Gen level: **** ");
      //hGen->Print("all");
      //printf("StoreRecUnf **** Reco level distribution: **** ");
      //hInput->Print("all");
      //printf("StoreRecUnf **** Response matrix: **** ");
      //hRspNorm->Print("all");

      TH1D* hPred = new TH1D(*hInput);
      for(int br = 1; br <= hPred->GetNbinsX(); br++)
      {
        double sum = 0.0;
        for(int bg = 1; bg <= hGen->GetNbinsX(); bg++)
          sum += hGen->GetBinContent(bg) * hRspNorm->GetBinContent(br, bg);
        hPred->SetBinContent(br, sum);
      }
      if(flagNorm)
        hPred->Scale(1.0 / hPred->Integral());
      //printf("StoreRecUnf **** Reco level prediction: **** ");
      //hPred->Print("all");
      double chi2Rec = utils::ttmd::Chi2(hInput, hInputNCov, hPred, flagNorm);
      printf("StoreRecUnf chi2Rec (norm = %d, dof = %d) = %.2f\n", flagNorm, hInput->GetNbinsX() - flagNorm, chi2Rec);

      // store txt
      // data
      TString fileName = dirOut + "/rec" + flagNormStr + "/data.txt";
      utils::ttmd::CheckFile(fileName);
      FILE* f = fopen(fileName.Data(), "w");
      for(int b = 1; b <= hInput->GetNbinsX(); b++)
        fprintf(f, "%.5e %.6e\n", hInput->GetBinContent(b), hInput->GetBinError(b));
      fclose(f);
      printf("StoreRecUnf Stored %s\n", fileName.Data());
      // response
      fileName = dirOut + "/rec" + flagNormStr + "/resp.txt";
      utils::ttmd::CheckFile(fileName);
      f = fopen(fileName.Data(), "w");
      for(int b = 1; b <= hRspNorm->GetNbinsY(); b++)
      {
        for(int bb = 1; bb <= hRspNorm->GetNbinsX(); bb++)
        {
          if(bb > 1)
            fprintf(f, " ");
          fprintf(f, "%.5e", hRspNorm->GetBinContent(bb, b));
        }
        fprintf(f, "\n");
      }
      fclose(f);
      printf("StoreRecUnf Stored %s\n", fileName.Data());
      // gen
      fileName = dirOut + "/rec" + flagNormStr + "/gen.txt";
      utils::ttmd::CheckFile(fileName);
      f = fopen(fileName.Data(), "w");
      for(int b = 1; b <= hGen->GetNbinsX(); b++)
        fprintf(f, "%.5e\n", hGen->GetBinContent(b));
      fclose(f);
      printf("StoreRecUnf Stored %s\n", fileName.Data());
      // pred
      fileName = dirOut + "/rec" + flagNormStr + "/pred.txt";
      utils::ttmd::CheckFile(fileName);
      f = fopen(fileName.Data(), "w");
      for(int b = 1; b <= hPred->GetNbinsX(); b++)
        fprintf(f, "%.5e\n", hPred->GetBinContent(b));
      fclose(f);
      printf("StoreRecUnf Stored %s\n", fileName.Data());
      // chi2
      fileName = dirOut + "/rec" + flagNormStr + "/chi2.txt";
      utils::ttmd::CheckFile(fileName);
      f = fopen(fileName.Data(), "w");
      fprintf(f, "%.5e %d %d\n", chi2Rec, hInput->GetNbinsX() - flagNorm, flagNorm);
      fclose(f);
      printf("StoreRecUnf Stored %s\n", fileName.Data());

      // unfolded
      if(unf)
      {
        TH1D* hUnf = new TH1D(*unf->hUnf);
        TH2D* hUnfCov = new TH2D(*unf->hEmatTotal);
        TH1D* hUnfGen = new TH1D(*unf->hGen);
        if(flagNorm)
        {
          utils::ttmd::Normalise(hUnf, hUnfCov);
          hUnfGen->Scale(1.0 / hUnfGen->Integral());
        }
        double chi2Unf = utils::ttmd::Chi2(hUnf, hUnfCov, hUnfGen, flagNorm);
        printf("StoreRecUnf chi2Unf (norm = %d, dof = %d) = %.2f\n", flagNorm, hInput->GetNbinsX() - flagNorm, chi2Unf);
        // unf
        fileName = dirOut + "/unf" + flagNormStr + "/data.txt";
        utils::ttmd::CheckFile(fileName);
        f = fopen(fileName.Data(), "w");
        for(int b = 1; b <= hUnf->GetNbinsX(); b++)
          fprintf(f, "%.5e\n", hUnf->GetBinContent(b));
        fclose(f);
        printf("StoreRecUnf Stored %s\n", fileName.Data());
        // unf cov
        fileName = dirOut + "/unf" + flagNormStr + "/datacov.txt";
        utils::ttmd::CheckFile(fileName);
        f = fopen(fileName.Data(), "w");
        for(int b = 1; b <= hUnfCov->GetNbinsY(); b++)
        {
          for(int bb = 1; bb <= hUnfCov->GetNbinsX(); bb++)
          {
            if(bb > 1)
              fprintf(f, " ");
            fprintf(f, "%+.5e", hUnfCov->GetBinContent(bb, b));
          }
          fprintf(f, "\n");
        }
        fclose(f);
        printf("StoreRecUnf Stored %s\n", fileName.Data());
        // pred
        fileName = dirOut + "/unf" + flagNormStr + "/pred.txt";
        utils::ttmd::CheckFile(fileName);
        f = fopen(fileName.Data(), "w");
        for(int b = 1; b <= hUnfGen->GetNbinsX(); b++)
          fprintf(f, "%.5e\n", hUnfGen->GetBinContent(b));
        fclose(f);
        printf("StoreRecUnf Stored %s\n", fileName.Data());
        // chi2
        fileName = dirOut + "/unf" + flagNormStr + "/chi2.txt";
        utils::ttmd::CheckFile(fileName);
        f = fopen(fileName.Data(), "w");
        fprintf(f, "%.5e %d %d\n", chi2Unf, hInput->GetNbinsX() - flagNorm, flagNorm);
        fclose(f);
        printf("StoreRecUnf Stored %s\n", fileName.Data());

        // chi2 from TUnfold
        fileName = dirOut + "/unf" + flagNormStr + "/all-chi2-fromUnf.txt";
        utils::ttmd::CheckFile(fileName);
        f = fopen(fileName.Data(), "w");
        double Chi2Data = unf->Chi2Data;
        double Chi2Reg = unf->Chi2Reg;
        int Dof = unf->Dof;
        fprintf(f, "%.5e %.5e %d\n", Chi2Data, Chi2Reg, Dof);
        fclose(f);
        printf("Unf chi2 Stored %s\n", fileName.Data());
      }
    }

std::vector<std::pair<TString, std::vector<TString> > > PlotterDiffXSec::CreateSysToPlotVector(const bool flagSingleVars, const std::vector<TString>& vSys)
    {
      const TString dirSingleVars = "single";
      std::vector<std::pair<TString, std::vector<TString> > > sysToPlot;
      for(auto sys : vSys)
      {
        //if(sys != "CR") continue;
        // sys variations
        std::vector<TString> vars = { "Nominal" };
        std::vector<TString> sysVars = configHelper->GetSysVars(sys);
        vars.insert(vars.end(), sysVars.begin(), sysVars.end());
        sysToPlot.push_back(std::make_pair(sys, vars));
        // single variations
        if(flagSingleVars)
        {
          for(unsigned int v = 1; v < vars.size(); v++)
          {
            const TString& var = vars[v];
            //if(!var.Contains("GLUONMOVETUNE")) continue;
            std::vector<TString> varsSingle = { "Nominal", var };
            sysToPlot.push_back(std::make_pair(dirSingleVars + "/" + var, varsSingle));
          }
        }
      }

      if(configHelper->gOnlyDominantSyst != 2) // if not skipping all
      {
        // all scale variations
        std::vector<TString> vars = { "Nominal", "MEFACSCALE_UP", "MEFACSCALE_DOWN", "MERENSCALE_UP", "MERENSCALE_DOWN", "MESCALE_UP", "MESCALE_DOWN" };
        sysToPlot.push_back(std::make_pair("ALL_MESCALES", vars));
      }

      return sysToPlot;
    }

void PlotterDiffXSec::PlotAllVars(const UnfoldingIOHandler* unfoldingIOHandler,
                     const std::vector<ZUnfold*> vUnf,
                     const std::map<TString, TH1D*>& mhVarGen,
                     const std::map<TString, TH1D*>& mhVarRec,
                     std::vector<std::map<TString, TH1D*> > vmhVarDat,
                     std::vector<std::pair<TString, std::vector<TString> > > sysToPlot,
                     const bool flagPlotGenRecOnly,
                     const TString& sufXSec,
                     const TString& dirOutBase
                     )
    {
      std::vector<TH1D*> vhVarGen;
      std::vector<TH1D*> vhVarGenRescaled;
      std::vector<TH1D*> vhVarRec;
      std::vector<TH1D*> vhVarRecRescaled;
      std::vector<std::vector<TH1D*> > vvhVarDat;
      vvhVarDat.resize(vUnf.size());
      for(const std::pair<TString, std::vector<TString> >& sys : sysToPlot)
      {
        TString baseName = sys.first;
        const std::vector<TString>& vars = sys.second;
        vhVarGen.resize(vars.size());
        vhVarGenRescaled.resize(vars.size());
        vhVarRec.resize(vars.size());
        vhVarRecRescaled.resize(vars.size());
        for(auto& vec : vvhVarDat)
          vec.resize(vars.size());
        for(unsigned int v = 0; v < vars.size(); v++)
        {
          const TString& var = vars[v];
          assert( v > 0 || var == "Nominal");
          // gen
          assert(mhVarGen.at(var) != NULL);
          vhVarGen[v] = (TH1D*) mhVarGen.at(var)->Clone();
          //if(v > 0)
          //  utils::ttmd::ScaleRelativeHistoVar(vhVarGen[v], configHelper->GetSysRescale(var), vhVarGen[0]);
          vhVarGen[v]->SetTitle(TString::Format("%s[I=%.3f]", var.Data(), vhVarGen[v]->Integral() / vhVarGen[0]->Integral()));
          vhVarGenRescaled[v] = (TH1D*) vhVarGen[v]->Clone();
          double scaleToNorm = vhVarGen[0]->Integral() / vhVarGen[v]->Integral();
          vhVarGenRescaled[v]->Scale(scaleToNorm);
          vhVarGenRescaled[v]->SetTitle(var);
          // rec
          assert(mhVarRec.at(var));
          vhVarRec[v] = (TH1D*) mhVarRec.at(var)->Clone();
          //if(v > 0)
          //  utils::ttmd::ScaleRelativeHistoVar(vhVarRec[v], configHelper->GetSysRescale(var), vhVarRec[0]);
          vhVarRec[v]->SetTitle(TString::Format("%s[I=%.3f]", var.Data(), vhVarRec[v]->Integral() / vhVarRec[0]->Integral()));
          vhVarRecRescaled[v] = (TH1D*) vhVarRec[v]->Clone();
          vhVarRecRescaled[v]->Scale(scaleToNorm);
          vhVarRecRescaled[v]->SetTitle(TString::Format("%s[I=%.3f]", var.Data(), vhVarRecRescaled[v]->Integral() / vhVarRecRescaled[0]->Integral()));
          // data
          if(unfoldingIOHandler->DoXSec)
          {
            for(unsigned int a = 0; a < vUnf.size(); a++)
            {
              assert(vmhVarDat[a][var]);
              vvhVarDat[a][v] = (TH1D*) vmhVarDat[a][var]->Clone();
              //if(v > 0)
              //  utils::ttmd::ScaleRelativeHistoVar(vvhVarDat[a][v], configHelper->GetSysRescale(var), vvhVarDat[a][0]);
              if(sufXSec == "abs")
                vvhVarDat[a][v]->SetTitle(TString::Format("%s[I=%.3f,~%.1f%%]", vUnf[a]->Suffix.Data(), vvhVarDat[a][v]->Integral() / vvhVarDat[a][0]->Integral(), 100 * utils::ttmd::GetMeanDeviationFromRef(vvhVarDat[a][v], vvhVarDat[a][0])));
              else
                vvhVarDat[a][v]->SetTitle(TString::Format("%s[~%.1f%%]", vUnf[a]->Suffix.Data(), 100 * utils::ttmd::GetMeanDeviationFromRef(vvhVarDat[a][v], vvhVarDat[a][0])));
            }
          }
        }
        if(flagPlotGenRecOnly)
        {
          TString fileName = dirOutBase + "/" + baseName + "-genorig";
          PlotVars(unfoldingIOHandler, vhVarGen, vhVarRec, fileName, 1);
          fileName = dirOutBase + "/" + baseName;
          PlotVars(unfoldingIOHandler, vhVarGenRescaled, vhVarRecRescaled, fileName, 0);
        }
        if(unfoldingIOHandler->DoXSec)
        {
          int flagAbs = (sufXSec == "abs");
          TString fileName = dirOutBase + "/" + baseName + "-xsec-" + sufXSec + "-all";
          PlotVars(unfoldingIOHandler, vvhVarDat, vhVarGenRescaled, vhVarRecRescaled, fileName, flagAbs);
          for(unsigned int a = 0; a < vUnf.size(); a++)
          {
            TString fileName = dirOutBase + "/" + baseName + "-xsec-" + sufXSec + "-" + vUnf[a]->Suffix;
            PlotVars(unfoldingIOHandler, { vvhVarDat[a] }, vhVarGenRescaled, vhVarRecRescaled, fileName, flagAbs);
          }
        }
      }
    }


double PlotterDiffXSec::SetScales(ZPredSetup& pred, const std::vector<TString>& vStrMu, const int mu)
    {
      double coef = ((mu % 2) == 0) ? 2.0 : 0.5;
      if(vStrMu[mu] == "#mu_{r}")
      {
        pred.Mur = coef;
        pred.Muf = 1.0;
      }
      else if(vStrMu[mu] == "#mu_{f}")
      {
        pred.Mur = 1.0;
        pred.Muf = coef;
      }
      else if(vStrMu[mu] == "#mu_{r,f}")
      {
        pred.Mur = pred.Muf = coef;
      }
      // 16.07.2018 support alternative scales
      else if(vStrMu[mu] == "#mu_{r,f} = H/2")
      {
        pred.Mur = 0.0;
        pred.Muf = 12.0;
      }
      else if(vStrMu[mu] == "#mu_{r,f} = H''/2")
      {
        pred.Mur = 0.0;
        pred.Muf = 13.0;
      }
      else
        throw std::logic_error("Error: should not be here");
      pred.SetNMu();
      return coef;
    }

void PlotterDiffXSec::PlotScaleVars(ZMadGraph* mg, ZPredSetup pred, const UnfoldingIOHandler* unfoldingIOHandler, const bool flagNorm, const std::vector<TString>& vStrMu, const TString& fileName)
    {
      UnfoldingXSecHandler* unfoldingXSecHandler = new UnfoldingXSecHandler();
      TH1D* hNom = unfoldingXSecHandler->GetPredictionHisto(unfoldingIOHandler, mg, pred, flagNorm);
      std::vector<TH1D*> vHMuU(5);
      std::vector<TH1D*> vHMuD(5);
      if(vStrMu.size() == 6)
      {
        vHMuU.resize(4);
        vHMuD.resize(4);
      }
      vHMuU[0] = vHMuD[0] = hNom;
      for(size_t mu = 0; mu < vStrMu.size(); mu++)
      {
        SetScales(pred, vStrMu, mu);
        TH1D* hVar = unfoldingXSecHandler->GetPredictionHisto(unfoldingIOHandler, mg, pred, flagNorm);
        hVar->SetTitle(TString::Format("%s [~%.1f%%]", vStrMu[mu].Data(), 100 * utils::ttmd::GetMeanDeviationFromRef(hVar, hNom)));
        if(mu < 6)
        {
          if((mu % 2) == 0)
            vHMuU[mu / 2 + 1] = hVar;
          else
            vHMuD[mu / 2 + 1] = hVar;
        }
        else if(mu == 6)
        {
          hVar->SetTitle(TString(hVar->GetTitle()).ReplaceAll("[~", "[[~"));
          vHMuU[4] = vHMuD[4] = hVar;
        }
        else
          assert(0);
      }
      PlotVars(unfoldingIOHandler, vHMuU, vHMuD, fileName, true, true);
    }

std::vector<TH1D*> PlotterDiffXSec::GetVHMtDep(const ZMadGraph* mg, const ZMeasur* measXSec)
    {
      std::vector<TH1D*> vhMtDep(mg->VMt.size());
      for(size_t m = 0; m < mg->VMt.size(); m++)
      //for(size_t m = 0; m < 3; m++)
      {
        double mt = mg->VMt[m];
        double diffMt = mt - mg->VMt[mg->NMtNom];
        TH1D* hNom = measXSec->GetNom();
        TH1D* hMtVar = (TH1D*) hNom->Clone();
        std::map<TString, TH1D*> map = measXSec->GetEigenvectors();
        if(diffMt < 0.0)
          hMtVar->Add(map["MASS_DOWN"]);
        else if(diffMt > 0.0)
          hMtVar->Add(map["MASS_UP"]);
        else
          hMtVar = (TH1D*) hNom->Clone();
        //printf("mt = %f  before\n", mt);
        //hMtVar->Print("all");
        if(diffMt != 0.0)
          utils::ttmd::ScaleRelativeHistoVar(hMtVar, TMath::Abs(diffMt), hNom);
        // reset relative stat. unc. taken from nominal
        for(int b = 1; b <= hNom->GetNbinsX(); b++)
          hMtVar->SetBinError(b, hMtVar->GetBinContent(b) * (hNom->GetBinError(b) / hNom->GetBinContent(b)));
        //hMtVar->SetBinError(b, 0.0);
        /*hMtVar->Add(hNom, -1.0);
        hMtVar->Scale(diffMt);
        hMtVar->Add(hNom, +1.0);
        hMtVar->Scale(1.0 / hMtVar->Integral());*/
        vhMtDep[m] = hMtVar;
        //printf("mt = %f  after\n", mt);
        //vhMtDep[m]->Print("all");
      }
      return vhMtDep;
    }

void PlotterDiffXSec::ProcessSystematics()
    {
      //bool flagCPShapeOnly = 1;
      const TString nameMeasDir = "all";
      const TString nameVarDir = "var";
      const TString nameUncDir = "unc";
      // const TString nameXFitterDir = "xfitter";
      const bool flagUseApplGrid = 1;
      const bool flagUseNomPDFCovMat = 0;
      const TString nameThMGDir = flagUseApplGrid ? "mg" : "mg-flagAG0";

      const bool flagOnlyAbsReg     = this->configHelper->ttmd_plot_only_Absolute_and_Regularized;
      const bool flagOnlyNorReg     = this->configHelper->ttmd_plot_only_Normalized_and_Regularized;
      const bool flagOnlyReg     = this->configHelper->ttmd_plot_only_Regularized;
      const bool flagAllChCompar    = this->configHelper->ttmd_plot_comparison_between_channels;
      const bool flagAllYearCompar    = this->configHelper->ttmd_plot_comparison_between_years;
      const bool flagAllYearMigMatrixSystsCompar    = this->configHelper->ttmd_plot_systs_mig_matrix_comparison_between_years;

      const bool flagPlotSysVars    = 0;
      const bool flagPlotSingleVars = 0;
      std::vector<TString> vSys = {}; // syst. variations to plot (if empty, plot all)
      //      /std::vector<TString> vSys = { "MASS", "MASS_CONSTMC", "MASS_KINRECO", "MASS_MCKINRECO"};
      //std::vector<TString> vSys = { "MASS", "MASS_CONSTMC", "MASS_KINRECO", "REWTOPPT"};
      //std::vector<TString> vSys = { "MASS", "MASS_CONSTMC", "REWTOPPT" };
      //vSys.push_back("MASS_KINRECO");
      //vSys.push_back("JESSinglePionECAL");
      //vSys.push_back("JESSinglePionHCAL");
      //vSys.push_back("JESFragmentation");
      const bool flagPlotUncSummary = this->configHelper->ttmd_plot_uncertainties_summary;
      //const bool flagPlotPDFEig     = 0;
      const bool flagPlotJESSrc     = 0;
      const bool flagPlotFOQCD      = 0;//TOP-18-004: 1
      const bool flagPlotFOQCDNorRegOnly = 0;//TOP-18-004: 1
      const bool flagPlotFOQCDAlsoAbsReg = 0;//TOP-18-004: 1
      const bool flagPlotPaper      = 1;
      //const bool flagPlotPas      = 1;
      const bool flagDoRj = 0;
      const std::vector<int> vSkipRj = {  };
      //std::vector<int> vSkipRj = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
      //const std::vector<int> vSkipRj = { 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 };
      //const std::vector<int> vSkipRj = { -1 };
      //std::vector<TString> vAltMC = { "POWHEGV2HERWIG", "AMCATNLOFXFX", "MADGRAPHMLM" };
      //std::vector<TString> vAltMCTitle = { "POW+HER", "FXFX+PYT", "MLM+PYT" };
      //std::vector<TString> vAltMC = { "MASS_UP", "MASS_DOWN" };
      //std::vector<TString> vAltMCTitle = { "MASS_UP", "MASS_DOWN" };
      //std::vector<TString> vAltMC = { "POWHEGV2HERWIG", "AMCATNLOFXFX" };
      //std::vector<TString> vAltMCTitle = { "POW+HER", "FXFX+PYT" };
      //std::vector<TString> vAltMCTitle = { "POW+HER", "MG5+PYT" };
      std::vector<TString> vAltMC = this->configHelper->ttmd_alternative_MC_samples;
      std::vector<TString> vAltMCTitle = this->configHelper->ttmd_alternative_MC_samples_legend_names;
      const int nmu = 8; // 1 nominal + 6 up-down + 1 alternative

      configHelper->gOnlyDominantSyst = 0;
      // set proper JES unc. flag --- FIXME: removes jes from plotting in some cases - perhaps not needed functionality
      //if(configHelper->gOnlyDominantSyst > 0)
      //  configHelper->gSysJES = 1;
      //else
      //  configHelper->gSysJES = 2;

      assert(flagPlotSysVars == 0 || flagOnlyNorReg == 0);

      std::pair<std::pair<TH1D*, TH1D*>, TH2D*> readXSec;
      std::vector<TH1D*> readCP;
      std::pair<std::pair<TH1D*, TH1D*>, TH2D*> readRecUnf;

      std::vector<TString> vSuf = CreateSuffVector();
      std::vector<ZUnfold*> vUnf = CreateZUnfoldVector();

      
      if(vSys.size() == 0)
        vSys = configHelper->GetAllSys(false);

      //std::vector<TString> allSysNormalization = {"LUMI", "BR"};
      //vSys.insert(vSys.end(), allSysNormalization.begin(), allSysNormalization.end());

      std::cout << "Systematics to be plotted: " << std::endl;
      for (auto syst_name: vSys) std::cout << syst_name << " ";
      std::cout << std::endl << std::endl;

      std::vector<std::pair<TString, std::vector<TString> > > sysToPlot = CreateSysToPlotVector(flagPlotSingleVars, vSys);

      std::vector<TString> allSys = configHelper->GetAllSys(false);
      allSys.insert(allSys.begin(), "Nominal");
      allSys = this->configHelper->FoldSystematicVarsListIntoEnvelops(allSys,false);
      std::cout << "Systematics included in allSys (used for plotting and systematics summary): " << std::endl;
      for (auto syst_name: allSys) std::cout << syst_name << " ";
      std::cout << std::endl << std::endl;

      std::vector<TString> allSysExp = configHelper->GetExpSys(false);
      allSysExp.insert(allSysExp.begin(), "Nominal");
      allSysExp = this->configHelper->FoldSystematicVarsListIntoEnvelops(allSysExp,false);
      std::cout << "Systematics included in allSysExp (used for systematics summary): " << std::endl;
      for (auto syst_name: allSysExp) std::cout << syst_name << " ";
      std::cout << std::endl << std::endl;

      std::vector<TString> allSysModW = configHelper->GetModWSys(false);
      allSysModW.insert(allSysModW.begin(), "Nominal");
      allSysModW = this->configHelper->FoldSystematicVarsListIntoEnvelops(allSysModW,false);
      std::cout << "Systematics included in allSysModW (used for systematics summary): " << std::endl;
      for (auto syst_name: allSysModW) std::cout << syst_name << " ";
      std::cout << std::endl << std::endl;

      std::vector<TString> allSysModI = configHelper->GetModISys(false);
      allSysModI.insert(allSysModI.begin(), "Nominal");
      allSysModI = this->configHelper->FoldSystematicVarsListIntoEnvelops(allSysModI,false);
      std::cout << "Systematics included in allSysModI (used for systematics summary): " << std::endl;
      for (auto syst_name: allSysModI) std::cout << syst_name << " ";
      std::cout << std::endl << std::endl;

      const auto& vars_to_fill = configHelper->GetAllVar(true);
      std::cout << "Variations to be processed into measurament: " << std::endl;
      for (auto syst_name: vars_to_fill) std::cout << syst_name << " ";
      std::cout << std::endl << std::endl;

      ZMadGraph* mg = NULL;
      if(flagPlotFOQCD)
      {
        mg = new ZMadGraph(configHelper,flagUseApplGrid, flagUseNomPDFCovMat);
        mg->DirNPCorr = configHelper->gAnalDir + "/" + configHelper->gNPCorrSubdir;
      }

      for(unsigned int x = 0; x < vXSec.size(); x++)
      {
        UnfoldingIOHandler* unfoldingIOHandler = vXSec[x];
        if(flagDoRj)
        {
          unfoldingIOHandler->DoRj = 1;
          unfoldingIOHandler->VSkip = vSkipRj;
          unfoldingIOHandler->DivideBW = 0;
          //unfoldingIOHandler->Var(2)->Get
        }

        // store x-sections for different channels, sufficex and algrorithms; map: suffix [alg+xsec], channel
        std::map<TString, std::vector<TH1D*> > mhAllSufXSec;
        std::map<TString, std::vector<std::pair<TH1D*, TH1D*> > > mhAllSufXSecUD;
        std::map<TString, std::vector<TH1D*> > mhAllSufXSecGen;
        std::map<TString, std::vector<TString> > mhAllSufXSecTitleChi2StatUncOnly;

        for(size_t c = 0; c < vCh.size(); c++)
        {
          const TString& ch = vCh[c];
          current_channel = ch;
          //TString dir = gAnalDir + "/" + ch + "/" + unfoldingIOHandler->Suffix;
          TString dir;
          if(ModeKinRecoTree == 0)
            dir = configHelper->gAnalDir + "/" + ch + "/" + unfoldingIOHandler->Suffix;
          else
            dir = configHelper->gAnalDir + "/" + ch + TString::Format("-kr%d/", ModeKinRecoTree) + unfoldingIOHandler->Suffix;
          const TString dirInput = dir;
          if(unfoldingIOHandler->DoRj)
            dir += "-rj";

          std::map<TString, TH1D*> mhVarGen;
          std::map<TString, TH1D*> mhVarRec;
          std::map<TString, TH1D*> mhVarRecWithBkg;
          std::vector<TH1D*> vCPNom; // for CP plot with systematics uncertainty band
          TH1D* hRecNom = NULL;
          TH1D* hGenNom = NULL;
          bool flagYieldsDone = 0;
          bool flagCPDrawn = 0;
          // detector level
          std::map<TString, TH2D*> mhVarDetResp;
          std::map<TString, TH1D*> mhVarDetDat;
          //std::map<TString, TH1D*> mhVarDetGen;
          TH2D* hDetRespNom = NULL;
          //TH1D* hDetDatNom = NULL;
          //TH1D* hDetGenNom = NULL;

          for(size_t s = 0; s < vSuf.size(); s++)
          {
            const TString suf = vSuf[s];
            if(flagOnlyNorReg && suf != "nor")
              continue;
            if(flagOnlyAbsReg && suf != "abs")
              continue;
            std::vector<std::map<TString, TH1D*> > vmhVarDat;
            vmhVarDat.resize(vUnf.size());
            int dof = unfoldingIOHandler->HGenC()->GetNbinsX();
            bool flagNorm = false;
            if(suf == "nor")
            {
              flagNorm = true;
              dof -= 1;
            }
            if(unfoldingIOHandler->DoRj)
            {
              assert(dof == 24);
              dof = 12;
            }

            std::vector<TH1D*> vXSecToPlot;
            std::vector<std::pair<TH1D*, TH1D*> > vXSecToPlotUD;
            TH2D* hXSecCovNom = NULL;
            TH1D* hXSecGenNom = NULL;
            for(size_t a = 0; a < vUnf.size(); a++){
                if(flagOnlyReg && vUnf[a]->Suffix != "reg")
                    continue;
                if(flagOnlyNorReg && vUnf[a]->Suffix != "reg")
                    continue;
                if(flagOnlyAbsReg && vUnf[a]->Suffix != "reg")
                    continue;
                vmhVarDat.push_back(std::map<TString, TH1D*>());
                TH1D* hXSecNom = NULL;

                /*for(const TString& var : vars_to_fill){
                    const TString& sys1 = configHelper->VarToSys(var);
                    printf("var = %s\n", var.Data());
                    printf("sys = %s\n", sys1.Data());
                }*/

                for(const TString& var : vars_to_fill){
                    // std::cout << "Processing: " << unfoldingIOHandler->Suffix << " <<--" << var << std::endl;
                    const TString& sys = configHelper->VarToSys(var);
                    //int nVarInSys = configHelper->GetSysVars(configHelper->VarToSys(var)).size();
                    const TString& dirTree = configHelper->GetVarDirPlot(var);
                    // read yields (once per each x-section)
                    //if(unfoldingIOHandler->DoYields && a == 0 && s == 0)
                    if(unfoldingIOHandler->DoYields && !flagYieldsDone)
                    {
                      readCP = ReadCP(dirInput + "/" + dirTree + "/cp.txt");
                      TH1D* hRec = (TH1D*) readCP[1]->Clone();
                      hRec->Add(readCP[2]);
                      hRec->Add(readCP[3]);
                      TH1D*& hGen = readCP[4];
                      if(var == "Nominal")
                      {
                        vCPNom = readCP;
                        hRecNom = hRec;
                        hGenNom = hGen;
                      }
                      else
                      {
                        utils::ttmd::ScaleRelativeHistoVar(hGen, configHelper->GetSysRescale(sys), hGenNom);
                        utils::ttmd::ScaleRelativeHistoVar(hRec, configHelper->GetSysRescale(sys), hRecNom);
                      }

		      //reverse engineer readCP histograms for a systematic with 16, 17, 18, 1617, 1618 and 1718 variations to check year-by-year rescaling is correct
		      //configHelper->doYearToYearCorrRescaleSanityCheck(var, dirInput);

                      // store rec in map
                      assert(mhVarRec.find(var) == mhVarRec.end());
                      //assert(mhVarRecWithBkg.find(var) == mhVarRecWithBkg.end());
                      mhVarRec[var] = new TH1D(*hRec);
                      //mhVarRecWithBkg[var] = new TH1D(*hRec);
                      // store gen in map
                      assert(mhVarGen.find(var) == mhVarGen.end());
                      mhVarGen[var] = new TH1D(*hGen);
                    }

                    // x-section
                    if(unfoldingIOHandler->DoXSec && unfoldingIOHandler->DoXSec != 2)
                    {
                      readXSec = ReadDataXsecAndCov(dirInput + "/" + dirTree + "/xsec-" + suf + "-" + vUnf[a]->Suffix);
                      TH1D*& hXSec = readXSec.first.first;
                      // 5.11.18 for KR Rj
                      if(flagDoRj)
                      {
                        if(var == "Nominal")
                        {
                          //hXSec->Print("all");
                          utils::ttmd::DoRjHisto(hXSec, readXSec.second);
                          //utils::ttmd::DoRjHisto(hXSec);
                          //for(int bx = 1; bx <= hXSec->GetNbinsX(); bx++)
                          //  for(int by = 1; by <= hXSec->GetNbinsX(); by++)
                          //    hXSecCovNom->SetBinContent(bx, by, (bx == by) ? (hXSec->GetBinError(bx) * hXSec->GetBinError(bx)) : 0);
                          //hXSec->Print("all");
                          //throw;
                          utils::ttmd::DoRjHisto(readXSec.first.second);
                        }
                        else
                          utils::ttmd::DoRjHisto(hXSec);
                      }
                      if(var == "Nominal")
                      {
                        hXSecNom = hXSec;
                        hXSecCovNom = readXSec.second;
                        hXSecGenNom = readXSec.first.second;
                      }
                      else
                      {
                        utils::ttmd::ScaleRelativeHistoVar(hXSec, configHelper->GetSysRescale(sys), hXSecNom);
                        if(suf == "abs" && !flagDoRj)
                        {
                          if(var == "LUMI_UP")
                            hXSec->Scale(1.0 / (1.0 + configHelper->lumiError));
                          else if(var == "LUMI_DOWN")
                            hXSec->Scale(1.0 / (1.0 - configHelper->lumiError));
                          if(var == "BR_UP")
                            hXSec->Scale(1.0 / (1.0 + configHelper->gBrrTtbarToEMuUncRel));
                          else if(var == "BR_DOWN")
                            hXSec->Scale(1.0 / (1.0 - configHelper->gBrrTtbarToEMuUncRel));
                        }
                      }
                      // store data in map
                      assert(vmhVarDat[a].find(var) == vmhVarDat[a].end());
                      //printf("var = %s %f\n", var.Data(), hXSec->GetBinContent(14));
                      //hXSec->Print("all");
                      vmhVarDat[a][var] = new TH1D(*hXSec);
                    }

                    // detector level and response matrix
                    if(unfoldingIOHandler->DoXSec > 1 && !flagYieldsDone)
                    {
                      const TString& dirTree = configHelper->GetVarDirPlot(var);
                      readRecUnf = ReadRecUnf(unfoldingIOHandler, dirInput + "/" + dirTree + "/rec-abs");
                      //TH1D*& hDet = readRecUnf.first.first;
                      if(var == "Nominal")
                      {
                        //hDetDatNom = hDet;
                        hDetRespNom = readRecUnf.second;
                        //hDetGenNom = readRecUnf.first.second;
                      }
                      else
                      {
                        /*utils::ttmd::ScaleRelativeHistoVar(hDet, configHelper->GetSysRescale(sys), hDetDatNom);
                        if(suf == "abs")
                        {
                          if(var == "LUMI_UP")
                            hDet->Scale(1.0 / (1.0 + lumiError));
                          else if(var == "LUMI_DOWN")
                            hDet->Scale(1.0 / (1.0 - lumiError));
                          if(var == "BR_UP")
                            hDet->Scale(1.0 / (1.0 + gBrrTtbarToEMuUncRel));
                          else if(var == "BR_DOWN")
                            hDet->Scale(1.0 / (1.0 - gBrrTtbarToEMuUncRel));
                        }*/
                      }
                      // store data in map
                      //if(var == "Nominal")
                      {
                        assert(mhVarDetResp.find(var) == mhVarDetResp.end());
                        mhVarDetResp[var] = new TH2D(*readRecUnf.second);
                        assert(mhVarDetDat.find(var) == mhVarDetDat.end());
                        mhVarDetDat[var] = new TH1D(*readRecUnf.first.first);
                      }
                    }
              }
              // TODO: this is dangerous
              if(unfoldingIOHandler->DoYields && !flagYieldsDone)
                flagYieldsDone = 1;

              // draw control plot with systematics
              //if(unfoldingIOHandler->DoYields && a == 0 && s == 0)
              //if(unfoldingIOHandler->DoYields && a == (vUnf.size() - 1) && s == (vSuf.size() - 1))

              if(unfoldingIOHandler->DoYields && !flagCPDrawn)
              {
                flagCPDrawn = 1;
                CPParams pars;
                // shape only variations
                ZMeasur* measRec = FillMeas(mhVarRec, allSys, NULL, true);
                ControlPlot(unfoldingIOHandler, vCPNom, measRec->GetTotalUD(), dir + "/" + nameMeasDir + "/cp", pars);
                PrintCPwithMCunc(vCPNom, measRec->GetTotalUD(), dir + "/" + nameMeasDir + "/cp.txt");
                if(flagPlotPaper && flagNorm)
                {
                  FlagPaper = 1;
                  ControlPlot(unfoldingIOHandler, vCPNom, measRec->GetTotalUD(), dir + "/" + nameMeasDir + "/cp-pap", pars);
                  FlagPaper = 2;
                  ControlPlot(unfoldingIOHandler, vCPNom, measRec->GetTotalUD(), dir + "/" + nameMeasDir + "/cp-pas", pars);
                  FlagPaper = 0;
                }
                // absolute variations
                ZMeasur* measRecAbs = FillMeas(mhVarRec, allSys, NULL, false);
                FlagPaper = 2;
                ControlPlot(unfoldingIOHandler, vCPNom, measRecAbs->GetTotalUD(), dir + "/" + nameMeasDir + "/cp-abs", pars);
                PrintCPwithMCunc(vCPNom, measRecAbs->GetTotalUD(), dir + "/" + nameMeasDir + "/cp-abs.txt");
              }

              // detector level
              if(unfoldingIOHandler->DoXSec == 2 || unfoldingIOHandler->DoXSec == 3)
              {
                std::map<TString, TH1D*> mDetConv;
                TH1D* hGenNom = mhVarGen["Nominal"];
                TH1D* hConvNom = utils::ttmd::Convolute(hGenNom, mhVarDetResp["Nominal"]);
                const TH1D* hDetNom = mhVarDetDat["Nominal"];
                for(std::map<TString, TH2D*>::iterator it = mhVarDetResp.begin(); it != mhVarDetResp.end(); it++)
                {
                  const TString& var = it->first;
                  const TH2D* hResp = it->second;
                  TH1D* hConv = utils::ttmd::Convolute(hGenNom, hResp);
                  //hConv->Print("all");
                  //hConvNom->Print("all");
                  //hDetNom->Print("all");
                  hConv->Add(hConvNom, -1.0);
                  hConv->Add(hDetNom);
                  if(flagNorm)
                    hConv->Scale(1.0 / hConv->Integral());
                  assert(mDetConv.find(var) == mDetConv.end());
                  mDetConv[var] = hConv;
                }
                //TH2D* hDetCov = new TH2D(*hXSecCovNom);
                TH2D* hDetCov = new TH2D("", "", hGenNom->GetNbinsX(), hGenNom->GetXaxis()->GetXbins()->GetArray(), hGenNom->GetNbinsX(), hGenNom->GetXaxis()->GetXbins()->GetArray());
                for(int bx = 0; bx < hDetCov->GetNbinsX(); bx++)
                  hDetCov->SetBinContent(bx + 1, bx + 1, TMath::Power(hDetNom->GetBinError(bx + 1), 2.0));
                utils::ttmd::Normalise(new TH1D(*hDetNom), hDetCov);
                ZMeasur* measXSecDet = FillMeas(mDetConv, allSys, hDetCov);
                TH1D* hDetPred = utils::ttmd::Convolute(hGenNom, hDetRespNom);
                //hDetRespNom->Print("all");
                //hGenNom->Print("all");
                //hDetPred->Print("all");
                if(flagNorm)
                  hDetPred->Scale(1.0 / hDetPred->Integral());
                measXSecDet->GetCovFullMatrix()->Print("all");
                TH1D* h = new TH1D(*measXSecDet->GetNom());
                h->Add(hDetPred, -1.0);
                h->Print("all");
                {
                  TCanvas* c = new TCanvas("", "", 600, 600);
                  measXSecDet->GetNom()->Draw();
                  hDetPred->SetLineColor(2);
                  hDetPred->Draw("same");
                  utils::ttmd::SaveCanvas(c, "tmpplot");
                }
                //double chi2 = utils::ttmd::Chi2(hDetPred, measXSecDet->GetCovFullMatrix(), measXSecDet->GetNom(), 0);
                double chi2 = measXSecDet->Chi2(hDetPred, flagNorm);
                printf("detector chi2/dof = %.2f/%d\n", chi2, dof);
              }

              // create x-section plot with systematics
              if(unfoldingIOHandler->DoXSec == 1 || unfoldingIOHandler->DoXSec == 3)
              {
                // make txt files
                ZMeasur* measXSec = FillMeas(vmhVarDat[a], allSys, hXSecCovNom);
                int skipBin = (suf == "nor");
                double chi2 = 0.0;
                if(flagDoRj)
                  chi2 = measXSec->Chi2(hXSecGenNom, vSkipRj);
                else
                  chi2 = measXSec->Chi2(hXSecGenNom, skipBin);
                printf("unfoldingIOHandler = %s  ALG = %s  SUF = %s  TOTAL CHI2 = %.1f\n", unfoldingIOHandler->Suffix.Data(), vUnf[a]->Suffix.Data(), suf.Data(), chi2);
                PrintDataXsec(measXSec->GetNom(), measXSec->GetSystUD(), hXSecGenNom, dir + "/" + nameMeasDir + "/xsec-" + suf + "-" + vUnf[a]->Suffix + "-vec.txt");
                TH2D* hCorStat = utils::ttmd::MakeCorrMatrix(utils::ttmd::MakeCorrMatrix(measXSec->GetCovStatMatrix()));
                PrintDataCorMatrix(hCorStat, dir + "/" + nameMeasDir + "/xsec-" + suf + "-" + vUnf[a]->Suffix + "-cor.txt");
                TH2D* hCorFull = utils::ttmd::MakeCorrMatrix(utils::ttmd::MakeCorrMatrix(measXSec->GetCovFullMatrix()));
                PrintDataCorMatrix(hCorFull, dir + "/" + nameMeasDir + "/xsec-" + suf + "-" + vUnf[a]->Suffix + "-corfull.txt");
                TH2D* hCorFullEnv = utils::ttmd::MakeCorrMatrix(utils::ttmd::MakeCorrMatrix(measXSec->GetCovFullMatrixEnvelope()));
                PrintDataCorMatrix(hCorFullEnv, dir + "/" + nameMeasDir + "/xsec-" + suf + "-" + vUnf[a]->Suffix + "-corfullenv.txt");

                // make txt file with all variations
                const std::map<TString, TH1D*>& mapAllSysVars = measXSec->GetEigenvectors();
                std::map<std::pair<TString, TString>, TH1D*> mapAllSysVarsForPrint;
                for(const auto& item : mapAllSysVars)
                {
                  const auto& var = item.first;
                  std::pair<TString, TString> pairNames = std::pair<TString, TString>(var, configHelper->GetShortVarName(var));
                  mapAllSysVarsForPrint[pairNames] = item.second;
                }
                PrintEigenvectors(measXSec->GetNom(), measXSec->GetTotalUD(), mapAllSysVarsForPrint,  dir + "/" + nameMeasDir + "/xsec-" + suf + "-" + vUnf[a]->Suffix + "-eig.txt");

                // make xFitter files
                //if(unfoldingIOHandler->Dim() == 2)
                // {
                //   TString nameFileXFitter = "xfitter-" + unfoldingIOHandler->Suffix + "-" + suf + "-" + vUnf[a]->Suffix;
                //   ZIOxFitter::PrintXFitterDat(unfoldingIOHandler, (TH1D*) measXSec->GetNom()->Clone(), measXSec->GetTotalUD(), mapAllSysVarsForPrint, dir + "/" + nameXFitterDir + "/" + nameFileXFitter + ".dat");
                //   ZIOxFitter::PrintXFitterCor(unfoldingIOHandler, utils::ttmd::MakeCorrMatrix(measXSec->GetCovStatMatrix()), dir + "/" + nameXFitterDir + "/" + nameFileXFitter + "__" + nameFileXFitter + ".dat");
                // }

                // make plot
                vXSecToPlot.push_back((TH1D*) measXSec->GetNom()->Clone());
                if(flagDoRj)
                  chi2 = utils::ttmd::Chi2(vXSecToPlot.back(), measXSec->GetCovFullMatrix(), hXSecGenNom, vSkipRj);
                else
                  chi2 = utils::ttmd::Chi2(vXSecToPlot.back(), measXSec->GetCovFullMatrix(), hXSecGenNom, (suf == "nor"));
                vUnf[a]->SetHistoStyle(vXSecToPlot.back());
                vXSecToPlot.back()->SetTitle(TString::Format("%s [%.0f/%d]", vUnf[a]->Suffix.Data(), chi2, dof));
                vXSecToPlotUD.push_back(measXSec->GetTotalUD());
                //hXSecGenNom->SetTitle("POWHEG+PYTHIA");
                //hXSecGenNom->SetTitle("POW+PYT");
                hXSecGenNom->SetTitle(TString::Format("%s [%.0f/%d]", "POW+PYT", chi2, dof));
                hXSecGenNom->SetLineColor(kBlack);
                hXSecGenNom->SetLineWidth(2);
                //XSec(unfoldingIOHandler, { vXSecToPlot.back() }, { vXSecToPlotUD.back() }, { hXSecGenNom }, dir + "/" + nameMeasDir + "/xsec-" + suf + "-" + vUnf[a]->Suffix, XSecParams(xsec));
                std::vector<TH1D*> vhXSecGenNom;
                hXSecGenNom->SetLineColor(utils::ttmd::GetDistinctColor(1));
                vhXSecGenNom.push_back(hXSecGenNom);
                //if(a == 1) printf("integral %f\n", hXSecGenNom->Integral());
                // new data histogram
                std::vector<TH1D*> vXSecToPlotCopy;
                vXSecToPlotCopy.push_back((TH1D*)vXSecToPlot.back()->Clone());
                vXSecToPlotCopy.back()->SetTitle(TString::Format("data, %s", vUnf[a]->Suffix.Data()));
                vXSecToPlotCopy.back()->SetLineColor(1);
                vXSecToPlotCopy.back()->SetMarkerColor(1);
                // alternative MC
                TString fileName = dir + "/all/chi2-allmc-"+ suf + ".txt";
                FILE* fout = fopen(fileName.Data(), "w");
                fprintf(fout, "%s %.2f %d\n", "POWHEGV2PYTHIA", chi2, dof);
                for(size_t mc = 0; mc < vAltMC.size(); mc++)
                {
                  TH1D* h = unfoldingIOHandler->GetPlainHisto();
                  TString fileName;
                  // 27.07.18 tmp solution to read MASS_UP, MASS_DOWN predictions from cp.txt
                  int column = 2;
                  TString fileNameBase = "gen.txt";
                  if(vAltMC[mc] == "MASS_UP" || vAltMC[mc] == "MASS_DOWN")
                  {
                    column = 7;
                    fileNameBase = "cp.txt";
                  }
                  if(ModeKinRecoTree == 0)
                    fileName = configHelper->gAnalDir + "/" + ch + "/" + unfoldingIOHandler->Suffix + "/" + vAltMC[mc] + "/" + fileNameBase;
                  else
                    fileName = configHelper->gAnalDir + "/" + ch + TString::Format("-kr%d/", ModeKinRecoTree) + unfoldingIOHandler->Suffix + "/" + vAltMC[mc] + "/" + fileNameBase;
                  utils::ttmd::ReadHistoKFactor(h, fileName, column, 3, flagNorm);
                  if(flagDoRj)
                    utils::ttmd::DoRjHisto(h);
                  // TODO achieve proper MC normalisation (POW+PYT now does not match NNLO x-sec)
                  if(!flagNorm)
                    h->Scale(hXSecGenNom->Integral() / h->Integral());
                  //if(a == 1) printf("integral mc = %d %f\n", mc, h->Integral());
                  h->SetLineColor(utils::ttmd::GetDistinctColor(mc + 2));
                  double chi2 = 0.0;
                  if(flagDoRj)
                    chi2 = utils::ttmd::Chi2(vXSecToPlotCopy.back(), measXSec->GetCovFullMatrix(), h, vSkipRj);
                  else
                    chi2 = utils::ttmd::Chi2(vXSecToPlotCopy.back(), measXSec->GetCovFullMatrix(), h, flagNorm);
                  h->SetTitle(TString::Format("%s [%.0f/%d]", vAltMCTitle[mc].Data(), chi2, dof));
                  h->SetLineWidth(2);
                  vhXSecGenNom.push_back(h);
                  fprintf(fout, "%s %.2f %d\n", vAltMC[mc].Data(), chi2, dof);
                }
                fclose(fout);
                //printf("nom: %f\n", measXSec->GetNom()->GetBinContent(14));
                //printf("ud: %f %f\n", vXSecToPlotUD.back().first->GetBinContent(14), vXSecToPlotUD.back().second->GetBinContent(13));
                XSec(unfoldingIOHandler, { vXSecToPlotCopy.back() }, { vXSecToPlotUD.back() }, vhXSecGenNom, dir + "/" + nameMeasDir + "/xsec-" + suf + "-" + vUnf[a]->Suffix, XSecParams(unfoldingIOHandler));

                // 17.10.18 plot with POW-PYT unc. band
                TH2D* hGenCov = new TH2D(*hXSecCovNom);
                std::map<TString, TH1D*> mhVarGenCopy;
                for(std::map<TString, TH1D*>::iterator it = mhVarGen.begin(); it != mhVarGen.end(); it++)
                {
                  TH1D* h = new TH1D(*it->second);
                  // only normalisation uncertainties
                  h->Scale(vhXSecGenNom[0]->Integral() / h->Integral());
                  mhVarGenCopy[it->first] = h;
                }
                ZMeasur* measGen = FillMeas(mhVarGenCopy, allSys, hGenCov);
                std::pair<TH1D*, TH1D*> measGenUD = measGen->GetTotalUD();
                // XSec(unfoldingIOHandler, { vXSecToPlotCopy.back() }, { vXSecToPlotUD.back() }, vhXSecGenNom, dir + "/" + nameMeasDir + "/xsec-bands-" + suf + "-" + vUnf[a]->Suffix, XSecParams(unfoldingIOHandler), measGenUD);

                // 26.10.18 paper plot with POW-PYT unc. band
                //if(flagPlotPaper && flagNorm)
                if(flagPlotPaper)
                {
                  std::vector<TH1D*> vhXSecGenNomCopy;
                  for(size_t i = 0; i < vhXSecGenNom.size(); i++)
                  {
                    vhXSecGenNomCopy.push_back(vhXSecGenNom[i]);
                    vhXSecGenNomCopy.back()->SetLineStyle(utils::ttmd::GetDistinctLineStyle(i));
                  }
                  TH1D* hd = new TH1D(*vXSecToPlotCopy.back());
                  //hd->SetTitle("Data");
                  hd->SetTitle(TString::Format("Data, dof=%d", dof));

                  hd->SetTitle(TString::Format("Data, dof=%d", dof));
                  FlagPaper = 1;
                  TString extraSuffix = "pap";
                  XSecParams pars = XSecParams(unfoldingIOHandler);
                  //pars.LogY = 1;
                  XSec(unfoldingIOHandler, { hd }, { vXSecToPlotUD.back() }, vhXSecGenNom, dir + "/" + nameMeasDir + "/xsec-bands-" + suf + "-" + vUnf[a]->Suffix + "-" + extraSuffix, pars, measGenUD);
                  FlagPaper = 2;
                  extraSuffix = "pas";
                  XSec(unfoldingIOHandler, { hd }, { vXSecToPlotUD.back() }, vhXSecGenNom, dir + "/" + nameMeasDir + "/xsec-bands-" + suf + "-" + vUnf[a]->Suffix + "-" + extraSuffix, pars, measGenUD);
                  FlagPaper = 0;
                  //hd->SetTitle("Data");
                  hd->SetTitle(TString::Format("Data, dof=%d", dof));
                }

                if(flagPlotPaper && flagNorm == 1 && vUnf[a]->Suffix == "reg" && (unfoldingIOHandler->Dim() == 2 || unfoldingIOHandler->Dim() == 3))
                {
                  TString fileName = dir + "/" + nameMeasDir + "/tab-xsec-" + suf + "-" + vUnf[a]->Suffix + ".tex";
                  PrintTexTabXSec(unfoldingIOHandler, measXSec, fileName);
                  fileName = dir + "/" + nameMeasDir + "/tab-corr-" + suf + "-" + vUnf[a]->Suffix + ".tex";
                  PrintTexTabCorr(unfoldingIOHandler, measXSec, fileName);
                  fileName = dir + "/" + nameMeasDir + "/tab-syst-" + suf + "-" + vUnf[a]->Suffix + ".tex";
                  // PrintTexTabSyst(configHelper ,unfoldingIOHandler, measXSec, allSys, fileName);
                }

                // store xsec for further plotting vs different channels
                TString nameSufAlg = suf + "-" + vUnf[a]->Suffix;
                mhAllSufXSec[nameSufAlg].push_back((TH1D*) vXSecToPlot.back()->Clone());
                mhAllSufXSecUD[nameSufAlg].push_back(vXSecToPlotUD.back());
                mhAllSufXSecGen[nameSufAlg].push_back((TH1D*) hXSecGenNom->Clone());
                // this is title with chi2 stat. unc. only
                double chi2StatOnly = 0.0;
                if(flagDoRj)
                  chi2StatOnly = utils::ttmd::Chi2(vXSecToPlot.back(), measXSec->GetCovStatMatrix(), hXSecGenNom, vSkipRj);
                else
                  chi2StatOnly = utils::ttmd::Chi2(vXSecToPlot.back(), measXSec->GetCovStatMatrix(), hXSecGenNom, (suf == "nor"));
                TString titleChi2StatOnly = TString::Format("%s [%.0f/%d]", vUnf[a]->Suffix.Data(), chi2StatOnly, vXSecToPlot.back()->GetNbinsX());
                mhAllSufXSecTitleChi2StatUncOnly[nameSufAlg].push_back(titleChi2StatOnly);

                // TODO move code below to appropriate place
                // fixed order QCD (MadGraph)
                if(unfoldingIOHandler->VMGHisto.size() && flagPlotFOQCD && (!flagPlotFOQCDNorRegOnly || (suf == "nor" && vUnf[a]->Suffix == "reg")
                //if(unfoldingIOHandler->VMGHisto.size() && (suf == "abs" && vUnf[a]->Suffix == "reg") && ( 1
                                                              || (flagPlotFOQCDAlsoAbsReg && suf == "abs" && vUnf[a]->Suffix == "reg")) )
                {
                  //double mt = 172.5;
                  ZPredSetup predNom(ZPredSetup::CT14nlo, 0, 1.0, 1.0, mg->MtNom);
                  ZPredSetup predEmpty;
                  ZPredSetup pred;
                  //std::vector<TString> vStrMu = {"#mu_{r}", "#mu_{r}", "#mu_{f}", "#mu_{f}", "#mu_{r,f}", "#mu_{r,f}"};
                  std::vector<TString> vStrMu = {"#mu_{r}", "#mu_{r}", "#mu_{f}", "#mu_{f}", "#mu_{r,f}", "#mu_{r,f}", "#mu_{r,f} = H/2"};
                  //std::vector<TString> vStrMu = {"#mu_{r}", "#mu_{r}", "#mu_{f}", "#mu_{f}", "#mu_{r,f}", "#mu_{r,f}", "#mu_{r,f} = H/2", "#mu_{r,f} = H''/2"};
                  // 16.07.18 mu12 not yet available for tt2j
                  bool flagSkipMu12 = 0;
                  // 26.10.18 available now
                  //if(unfoldingIOHandler->NPertOrders == 3)
                  //  flagSkipMu12 = 1;
                  if(flagSkipMu12)
                    vStrMu.resize(6);

                  // read predictions and PDF unc.
                  for(size_t mt = 0; mt < mg->VMt.size(); mt++)
                  {
                    for(size_t mu = 0; mu < nmu; mu++)
                    {
                      if(flagSkipMu12 && mu == 7)
                        continue;
                      // if using PDF cov. matrix from nominal setup, skip all variations
                      //bool flagNoCov = (flagUseNomPDFCovMat && (mu != 0 || mt != 0));
                      //if(flagUseNomPDFCovMat && (mu != 0 || mt != 0))
                      //  continue;
                      // if not ApplGrid mode, skip all scales variations
                      if(!flagUseApplGrid && mu != 0)
                        continue;
                      // prediction with varied scales
                      ZPredSetup pred;
                      if(mu == 0)
                        pred.Mur = pred.Muf = 1.0;
                      else
                        SetScales(pred, vStrMu, mu - 1);
                    }
                  }
                }

                // summary plots of uncertainties
                if(flagPlotUncSummary)
                {
                  // total
                  std::vector<TH1D*> vhUncTotal = GetVectorUncUD(*measXSec, "total");
                  std::vector<TH1D*> vhUncUnf = GetVectorUncUD(*measXSec, "stat");
                  std::vector<std::vector<TH1D*> > vhUncSyst;
                  vhUncSyst.push_back(GetVectorUncUD(*measXSec, "syst"));
                  vhUncSyst.back()[0]->SetMarkerStyle(31);
                  vhUncSyst.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], allSysExp), "exp"));
                  vhUncSyst.back()[0]->SetMarkerStyle(24);
                  vhUncSyst.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], allSysModW), "mod (w)"));
                  vhUncSyst.back()[0]->SetMarkerStyle(26);
                  vhUncSyst.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], allSysModI), "mod (i)"));
                  vhUncSyst.back()[0]->SetMarkerStyle(32);
                  TString fileName = dir + "/" + nameUncDir + "/total-" + suf + "-" + vUnf[a]->Suffix;
                  UncSummaryParams pars;
                  pars.FlagAbs = (suf == "abs");
                  PlotUncSummary(unfoldingIOHandler, vhUncSyst, vhUncTotal, vhUncUnf, fileName, pars);

                  // summary plots of experimental uncertainties
                  std::vector<std::vector<TH1D*> > vhUncExp;
                  // std::vector<TString> vBTAG = {"Nominal", "BTAG", "BTAG_ETA", "BTAG_LJET", "BTAG_LJET_ETA", "BTAG_LJET_PT", "BTAG_PT"};
                  std::vector<TString> vBTAG = {"Nominal", "BTAG", "BTAG_LJET"};
                  std::vector<TString> vLEPT = {"Nominal", "ELE_ID", "ELE_RECO", "ELE_SCALESMEARING", "MUON_ID", "MUON_ISO", "MUON_SCALE"};
                  std::vector<TH1D*> vhUncLEPT_UD = GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList(vLEPT, year)), "LEPT");
                  //std::vector<TString> vBTAG = {"Nominal", "BTAG", "BTAG_LJET"};
                  vhUncExp.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList(vBTAG, year)), "BTAG"));
                  vhUncExp.back()[0]->SetMarkerStyle(24);
                  // 1.10.18 METJES was just a cross check: it is small
                  //std::vector<TString> srcJETMET = { "Nominal", "JER", "UNCLUSTERED", "METJES" };
                  //std::vector<TString> srcJETMET = { "Nominal", "JER", "UNCLUSTERED" };
                  std::vector<TString> srcJETMET_noJES = { "Nominal", "JEREta0", "JEREta1", "UNCLUSTERED" };
                  std::vector<TString> srcJETMET = srcJETMET_noJES;
                  //std::vector<TString> srcJETMET = { "Nominal", "JEREta0", "JEREta1" };

                  std::vector<TString> vJES = {"Nominal"};
                  for(auto& src : configHelper->GetJECSrc(false)){
                    if(!src.Contains("Total")){
                      srcJETMET.push_back(src);
                      vJES.push_back(src);
                    }
                  }
                  
                  std::vector<TH1D*> vhUncJES_UD = GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList(vJES, year)), "JES");
                  std::vector<TH1D*> vhUncJETMET_UD = GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList(srcJETMET, year)), "JETMET");
                  vhUncExp.push_back(vhUncJETMET_UD);
                  vhUncExp.back()[0]->SetMarkerStyle(25);
                  if(ModeKinRecoTree == 0)
                  {
                    vhUncExp.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList({"Nominal", "KIN"}, year)), "KIN"));
                    vhUncExp.back()[0]->SetMarkerStyle(26);
                  }
                  vhUncExp.push_back(vhUncLEPT_UD);
                  vhUncExp.back()[0]->SetMarkerStyle(32);
                  if (configHelper->IsExpSys("PU")){
                    vhUncExp.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList({"Nominal", "PU"}, year)), "PU"));
                    vhUncExp.back()[0]->SetMarkerStyle(27);
                  }
                  vhUncExp.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList({"Nominal", "LUMI"}, year)), "LUMI"));
                  vhUncExp.back()[0]->SetMarkerStyle(29);
                  //vhUncExp.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], {"Nominal", "TRIG", "TRIG_ETA", "L1PREFIRING"}), "TRIG"));
                  vhUncExp.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList({"Nominal", "TRIG"}, year)), "TRIG"));
                  vhUncExp.back()[0]->SetMarkerStyle(28);
                  vhUncExp.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList({"Nominal", "L1PREFIRING"}, year)), "L1PREF"));
                  vhUncExp.back()[0]->SetMarkerStyle(32);
                  vhUncExp.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList({"Nominal", "BG", "DY"}, year)), "BG"));
                  vhUncExp.back()[0]->SetMarkerStyle(30);
                  /*if(suf == "abs" && !flagDoRj){
                    vhUncExp.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], {"Nominal", "LUMI", "BR"}), "NORM"));
                    vhUncExp.back()[0]->SetMarkerStyle(31);
                  }*/
                  fileName = dir + "/" + nameUncDir + "/exp-" + suf + "-" + vUnf[a]->Suffix;
                  pars.NLegendRows = 2;
                  PlotUncSummary(unfoldingIOHandler, vhUncExp, vhUncSyst[1], {}, fileName, pars);

                  // summary plots of mod (w) uncertainties
                  /*
                  TString isr_folder, fsr_folder;
                  if (year=="2016"){
                      isr_folder = "PSISRSCALE";
                      fsr_folder = "PSFSRSCALE";
                  }
                  else{
                      isr_folder = "PSISRSCALE_2";
                      fsr_folder = "PSFSRSCALE_2";
                  } */
                  TString isr_folder = "PSSCALE_WEIGHT_4";
                  TString fsr_folder = "PSSCALE_WEIGHT_5";
                  std::vector<std::vector<TH1D*> > vhUncModW;
                  
                  std::vector<TString> vScale = this->configHelper->GetEnvelopeVars("TOT_SCALE", false);
                  vScale.insert(vScale.begin(), "Nominal");
                  // std::vector<TString> vScale = {"Nominal", "MESCALE", "MEFACSCALE", "MERENSCALE"};
                  std::vector<TH1D*> Scale_UD = GetVectorUncUD(*FillMeas(vmhVarDat[a], vScale), "TOT_SCALE");
                  vhUncModW.push_back(Scale_UD);
                  vhUncModW.back()[0]->SetMarkerStyle(24);

                  std::vector<TString> vBFrag = this->configHelper->GetEnvelopeVars("TOT_BFRAG", false);
                  vBFrag.insert(vBFrag.begin(), "Nominal");
                  std::vector<TH1D*> BFrag_UD = GetVectorUncUD(*FillMeas(vmhVarDat[a], vBFrag), "TOT_BFRAG");
                  vhUncModW.push_back(BFrag_UD);
                  vhUncModW.back()[0]->SetMarkerStyle(25);

                  vhUncModW.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList({"Nominal", "BSEMILEP"}, year)), "BSEMILEP"));
                  vhUncModW.back()[0]->SetMarkerStyle(26);
                  vhUncModW.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList({"Nominal", "PDF_ALPHAS"}, year)), "PDF_ALPHAS"));
                  vhUncModW.back()[0]->SetMarkerStyle(27);
                  vhUncModW.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList({"Nominal", fsr_folder}, year)), "PSFSRSCALE"));
                  vhUncModW.back()[0]->SetMarkerStyle(28);
                  vhUncModW.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList({"Nominal", isr_folder}, year)), "PSISRSCALE"));
                  vhUncModW.back()[0]->SetMarkerStyle(32);
                  /*std::vector<TString> srcPDF = configHelper->GetModPDFEig();
                  srcPDF.insert(srcPDF.begin(), "Nominal");
                  std::vector<TH1D*> vhPDFTot = GetVectorUncUD(*FillMeas(vmhVarDat[a], srcPDF), "PDF_CT14");
                  vhUncModW.push_back(vhPDFTot);
                  vhUncModW.back()[0]->SetMarkerStyle(32);
                  std::vector<TString> srcPDFAS = { "Nominal", "PDF_CT14_AS" };
                  vhUncModW.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], srcPDFAS), "PDF_CT14_AS"));
                  vhUncModW.back()[0]->SetMarkerStyle(27);*/
                  fileName = dir + "/" + nameUncDir + "/modw-" + suf + "-" + vUnf[a]->Suffix;
                  pars.NLegendRows = 2;
                  pars.Rat = 8.0;
                  PlotUncSummary(unfoldingIOHandler, vhUncModW, vhUncSyst[2], {}, fileName, pars);

                  // summary plots of mod (i) uncertainties
                  std::vector<std::vector<TH1D*> > vhUncModI;
                  // Checking if "MASS" was excluded
                  if (configHelper->IsModISys("MASS")){
                    vhUncModI.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList({"Nominal", "MASS"}, year)), "MASS"));
                    vhUncModI.back()[0]->SetMarkerStyle(24);
                  }
                  vhUncModI.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList({"Nominal", "MATCH"}, year)), "MATCH"));
                  vhUncModI.back()[0]->SetMarkerStyle(25);
                  
                  /*
                  vhUncModI.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList({"Nominal", fsr_folder})), "PSFSRSCALE"));
                  vhUncModI.back()[0]->SetMarkerStyle(26);
                  vhUncModI.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList({"Nominal", isr_folder})), "PSISRSCALE"));
                  vhUncModI.back()[0]->SetMarkerStyle(32);*/
                  std::vector<TString> vCR = this->configHelper->GetEnvelopeVars("TOT_COLORREC", false);
                  vCR.insert(vCR.begin(), "Nominal");
                  std::vector<TH1D*> CR_UD = GetVectorUncUD(*FillMeas(vmhVarDat[a], {"Nominal", "TOT_COLORREC"}), "CR");
                  vhUncModI.push_back(CR_UD);
                  vhUncModI.back()[0]->SetMarkerStyle(27);

                  // Checking if "UETUNE" was excluded
                bool includeUETune = false;
                for (auto modISys: allSysModI) {
                  if (modISys.Contains("UETUNE")) includeUETune = true;
                }
                  if (includeUETune){
                      vhUncModI.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], GetUncSummaryList({"Nominal", "UETUNE"}, year)), "UETUNE"));
                      vhUncModI.back()[0]->SetMarkerStyle(28);
                  }
                  fileName = dir + "/" + nameUncDir + "/modi-" + suf + "-" + vUnf[a]->Suffix;
                  pars.NLegendRows = 2;
                  pars.Rat = 12.0;
                  PlotUncSummary(unfoldingIOHandler, vhUncModI, vhUncSyst[3], {}, fileName, pars);

                  // summary plots of ALL_MESCALE uncertainties
                  std::vector<std::vector<TH1D*> > vhUncAllScales;
                  int marker_style_id = 23;
                  for (auto src: vScale) {
                    if (src == "Nominal") continue;
                    std::vector<TString> summary_list = GetUncSummaryList({src}, year);
                    summary_list.insert(summary_list.begin(), "Nominal");
                    vhUncAllScales.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], summary_list), src));
                    vhUncAllScales.back()[0]->SetMarkerStyle(marker_style_id);
                  }
                  fileName = dir + "/" + nameUncDir + "/allscale-" + suf + "-" + vUnf[a]->Suffix;
                  pars.NLegendRows = 1;
                  pars.Rat = 8.0;
                  PlotUncSummary(unfoldingIOHandler, vhUncAllScales, Scale_UD, {}, fileName, pars);

                  // summary plots of COLOR_RECONECTION (CR) uncertainties
                  std::vector<std::vector<TH1D*> > vhUncAllColorRec;
                  marker_style_id = 23;
                  for (auto src: vCR) {
                    if (src == "Nominal") continue;
                    std::vector<TString> summary_list = GetUncSummaryList({src}, year);
                    summary_list.insert(summary_list.begin(), "Nominal");
                    vhUncAllColorRec.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], summary_list), src));
                    vhUncAllColorRec.back()[0]->SetMarkerStyle(marker_style_id);
                  }
                  fileName = dir + "/" + nameUncDir + "/all_CR-" + suf + "-" + vUnf[a]->Suffix;
                  pars.NLegendRows = 1;
                  pars.Rat = 8.0;
                  PlotUncSummary(unfoldingIOHandler, vhUncAllColorRec, CR_UD, {}, fileName, pars);

                  // summary plots of LEPT uncertainties
                  std::vector<std::vector<TH1D*> > vhUncLepton;
                  marker_style_id = 23;
                  for (auto src: vLEPT) {
                    if (src == "Nominal") continue;
                    std::vector<TString> summary_list = GetUncSummaryList({src}, year);
                    summary_list.insert(summary_list.begin(), "Nominal");
                    vhUncLepton.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], summary_list), src));
                    vhUncLepton.back()[0]->SetMarkerStyle(marker_style_id);
                    marker_style_id++;
                  }
                  fileName = dir + "/" + nameUncDir + "/lepton-" + suf + "-" + vUnf[a]->Suffix;
                  pars.NLegendRows = 2;
                  pars.Rat = 15.0;
                  PlotUncSummary(unfoldingIOHandler, vhUncLepton, vhUncLEPT_UD, {}, fileName, pars);

                  // summary plots of BFRAG uncertainties
                  std::vector<std::vector<TH1D*> > vhUncBFrag;
                  marker_style_id = 23;
                  for (auto src: vBFrag) {
                    if (src == "Nominal") continue;
                    std::vector<TString> summary_list = GetUncSummaryList({src}, year);
                    summary_list.insert(summary_list.begin(), "Nominal");
                    vhUncBFrag.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], summary_list), src));
                    vhUncBFrag.back()[0]->SetMarkerStyle(marker_style_id);
                    marker_style_id++;
                  }
                  fileName = dir + "/" + nameUncDir + "/all_bfrag-" + suf + "-" + vUnf[a]->Suffix;
                  pars.NLegendRows = 1;
                  pars.Rat = 4.0;
                  PlotUncSummary(unfoldingIOHandler, vhUncBFrag, BFrag_UD, {}, fileName, pars);

                  // summary plots of JETMET uncertainties
                  std::vector<std::vector<TH1D*> > vhUncJETMETsources;
                  marker_style_id = 23;
                  for (auto src: srcJETMET_noJES) {
                    if (src == "Nominal") continue;
                    std::vector<TString> summary_list = GetUncSummaryList({src}, year);
                    summary_list.insert(summary_list.begin(), "Nominal");
                    vhUncJETMETsources.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], summary_list), src));
                    vhUncJETMETsources.back()[0]->SetMarkerStyle(marker_style_id);
                    marker_style_id++;
                  }
                  vhUncJETMETsources.push_back(vhUncJES_UD);
                  vhUncJETMETsources.back()[0]->SetMarkerStyle(marker_style_id);
                  
                  fileName = dir + "/" + nameUncDir + "/JETMET-" + suf + "-" + vUnf[a]->Suffix;
                  pars.NLegendRows = 1;
                  pars.Rat = 20.0;
                  PlotUncSummary(unfoldingIOHandler, vhUncJETMETsources, vhUncJETMET_UD, {}, fileName, pars);



                  // ************* summary plots of JEC uncertainties *************
                  if(flagPlotJESSrc)
                  {
                    fileName = dir + "/" + nameUncDir + "/JESXXX-" + suf + "-" + vUnf[a]->Suffix;
                    pars.NLegendRows = 2;
                    //pars.Rat = 15.0;
                    pars.Rat = -1;

                    // SubTotalPileUp
                    std::vector<TString> vJESsrc = { "PileUpPtRef", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF", "PileUpDataMC" };
                    std::vector<TH1D*> hJESSumPileUp = PlotJESSubtotalUnc(unfoldingIOHandler, vmhVarDat[a], "SubTotalPileUp", vJESsrc, vhUncUnf, fileName, pars);

                    // SubTotalRelative
                    vJESsrc = { "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB", "RelativePtEC1", "RelativePtEC2",
                                "RelativePtHF", "RelativeBal", "RelativeFSR", "RelativeStatFSR", "RelativeStatEC", "RelativeStatHF"};
                    pars.NLegendRows = 3;
                    std::vector<TH1D*> hJESSumRelative = PlotJESSubtotalUnc(unfoldingIOHandler, vmhVarDat[a], "SubTotalRelative", vJESsrc, vhUncUnf, fileName, pars);

                    // SubTotalPt
                    vJESsrc = { "Fragmentation", "SinglePionECAL", "SinglePionHCAL" };
                    pars.NLegendRows = 2;
                    std::vector<TH1D*> hJESSumPt = PlotJESSubtotalUnc(unfoldingIOHandler, vmhVarDat[a], "SubTotalPt", vJESsrc, vhUncUnf, fileName, pars);

                    // SubTotalScale
                    vJESsrc = { "AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias" };
                    pars.NLegendRows = 2;
                    std::vector<TH1D*> hJESSumScale = PlotJESSubtotalUnc(unfoldingIOHandler, vmhVarDat[a], "SubTotalScale", vJESsrc, vhUncUnf, fileName, pars);

                    // SubTotalAbsolute
                    pars.NLegendRows = 2;
                    vJESsrc = { "SubTotalPt", "SubTotalScale" };
                    PlotJESSubtotalUnc(unfoldingIOHandler, vmhVarDat[a], "SubTotalAbsolute", vJESsrc, vhUncUnf, fileName, pars);

                    // SubTotalMC
                    vJESsrc.clear();
                    pars.NLegendRows = 1;
                    //pars.Rat = 20.0;
                    for(auto& src : configHelper->GetJECSrc(false))
                      if(!src.BeginsWith("JESPileUpPt") && !src.Contains("Total"))
                        vJESsrc.push_back(TString(src.Data() + 3));
                    PlotJESSubtotalUnc(unfoldingIOHandler, vmhVarDat[a], "SubTotalMC", vJESsrc, vhUncUnf, fileName, pars, false);

                    // TotalNoFlavor
                    vJESsrc.clear();
                    for(auto& src : configHelper->GetJECSrc(false))
                      if(!src.BeginsWith("JESFlavor") && !src.Contains("Total"))
                        vJESsrc.push_back(TString(src.Data() + 3));
                    PlotJESSubtotalUnc(unfoldingIOHandler, vmhVarDat[a], "TotalNoFlavor", vJESsrc, vhUncUnf, fileName, pars, false);

                    // TotalNoTime
                    vJESsrc.clear();
                    for(auto& src : configHelper->GetJECSrc(false))
                      if(!src.BeginsWith("JESTime") && !src.Contains("Total"))
                        vJESsrc.push_back(TString(src.Data() + 3));
                    PlotJESSubtotalUnc(unfoldingIOHandler, vmhVarDat[a], "TotalNoTime", vJESsrc, vhUncUnf, fileName, pars, false);

                    // TotalNoFlavorNoTime
                    vJESsrc.clear();
                    for(auto& src : configHelper->GetJECSrc(false))
                      if(!src.BeginsWith("JESFlavor") && !src.BeginsWith("JESTime") && !src.Contains("Total"))
                        vJESsrc.push_back(TString(src.Data() + 3));
                    PlotJESSubtotalUnc(unfoldingIOHandler, vmhVarDat[a], "TotalNoFlavorNoTime", vJESsrc, vhUncUnf, fileName, pars, false);

                    // Total
                    vJESsrc.clear();
                    for(auto& src : configHelper->GetJECSrc(false))
                      if(!src.Contains("Total"))
                        vJESsrc.push_back(TString(src.Data() + 3));
                    std::vector<TH1D*> hJESSumTotal = PlotJESSubtotalUnc(unfoldingIOHandler, vmhVarDat[a], "Total", vJESsrc, vhUncUnf, fileName, pars, false);

                    // Total vs historical "JES"
                    std::vector<TH1D*> uncTotalOneVar = GetVectorUncUD(*FillMeas(vmhVarDat[a], {"Nominal", "JES"}), "JES one variation");
                    fileName = dir + "/" + nameUncDir + "/jes-Total-" + suf + "-" + vUnf[a]->Suffix;
                    uncTotalOneVar.insert(uncTotalOneVar.end(), vhUncUnf.begin(), vhUncUnf.end());
                    PlotUncSummary(unfoldingIOHandler, {}, hJESSumTotal, uncTotalOneVar, fileName, pars);

                    // Plot all contributions
                    std::vector<std::vector<TH1D*> > vhJESContib;
                    vhJESContib.push_back(hJESSumPileUp);
                    vhJESContib.back()[0]->SetMarkerStyle(24);
                    vhJESContib.push_back(hJESSumRelative);
                    vhJESContib.back()[0]->SetMarkerStyle(25);
                    vhJESContib.push_back(hJESSumPt);
                    vhJESContib.back()[0]->SetMarkerStyle(26);
                    vhJESContib.push_back(hJESSumScale);
                    vhJESContib.back()[0]->SetMarkerStyle(32);
                    vhJESContib.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], {"Nominal", "JESFlavorQCD"}), "FlavorQCD"));
                    vhJESContib.back()[0]->SetMarkerStyle(27);
                    vhJESContib.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], {"Nominal", "JESTimePtEta"}), "TimePtEta"));
                    vhJESContib.back()[0]->SetMarkerStyle(28);
                    vhJESContib.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], {"Nominal", "JER"}), "JER"));
                    vhJESContib.back()[0]->SetMarkerStyle(30);
                    vhJESContib.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], {"Nominal", "UNCLUSTERED"}), "Unclustered MET"));
                    vhJESContib.back()[0]->SetMarkerStyle(3);
                    //vhJESContib.push_back(GetVectorUncUD(*FillMeas(vmhVarDat[a], {"Nominal", "METJES"}), "MET JES"));
                    //vhJESContib.back()[0]->SetMarkerStyle(5);

                    //vJESsrc = { "Nominal", "JER", "UNCLUSTERED", "METJES" };
                    vJESsrc = { "Nominal", "JER", "UNCLUSTERED" };
                    for(auto& src : configHelper->GetJECSrc(false))
                      if(!src.Contains("Total"))
                        vJESsrc.push_back(src);
                    std::vector<TH1D*> uncTotal = GetVectorUncUD(*FillMeas(vmhVarDat[a], vJESsrc), "JETMETTotal");
                    std::vector<TH1D*> uncTotalOrig = GetVectorUncUD(*FillMeas(vmhVarDat[a], {"Nominal", "JESTotal"}), "JESTotal-src.");

                    fileName = dir + "/" + nameUncDir + "/jes-all-" + suf + "-" + vUnf[a]->Suffix;
                    pars.NLegendRows = 3;
                    //pars.Rat = 20.0;
                    //uncTotal.insert(uncTotal.end(), vhUncUnf.begin(), vhUncUnf.end());
                    PlotUncSummary(unfoldingIOHandler, vhJESContib, uncTotal, {}, fileName, pars);
                  }
                }
              }
            }

            // draw all algorithms x-section
            if(unfoldingIOHandler->DoXSec && !flagOnlyNorReg)
              XSec(unfoldingIOHandler, vXSecToPlot, vXSecToPlotUD, { hXSecGenNom }, dir + "/" + nameMeasDir + "/xsec-" + suf, XSecParams(unfoldingIOHandler));

            // draw systematics
            if(flagPlotSysVars)
              if(unfoldingIOHandler->DoYields)
                PlotAllVars(unfoldingIOHandler, vUnf, mhVarGen, mhVarRec, vmhVarDat, sysToPlot, (s == 0), suf, dir + "/" + nameVarDir);
          }
        } // end of channel loop

        // draw all channel x-section
        if(flagAllChCompar)
        {
          if(unfoldingIOHandler->DoXSec){
            std::map<TString, std::vector<TH1D*> > mhAllSufXSec_by_year;
            std::map<TString, std::vector<std::pair<TH1D*, TH1D*> > > mhAllSufXSecUD_by_year;
            std::map<TString, std::vector<TH1D*> > mhAllSufXSecGen_by_year;

            std::vector<int> vMarkers = {24, 25, 26, 32};
            std::vector<int> vColors = {4, 8, kOrange - 4, 2};

            std::vector<TString> channels_list = { "ll", "ee", "emu", "mumu"};
            std::vector<TString> undfolding_opts = {"-nor-reg", "-abs-reg"};
            TString plots_version = configHelper->analysisSuffix;

            for (TString undfolding_opt: undfolding_opts){
                  TString out_dir = "";
                  if(ModeKinRecoTree == 0) out_dir = configHelper->gAnalDir + "/all_channel_comparison" + "/" + unfoldingIOHandler->Suffix;
                  else out_dir = configHelper->gAnalDir + "/all_channel_comparison" + "/" + TString::Format("-kr%d", ModeKinRecoTree) + "-" + unfoldingIOHandler->Suffix;

                  std::vector<TH1D*> vector_hDat;
                  std::vector<TH1D*> vector_hGen;
                  std::vector<std::pair<TH1D*, TH1D*> > vector_TotUD;

                  unsigned int ch_counter = 0;

                  for (TString i_ch: channels_list){
                    TString in_dir = "";
                    if(ModeKinRecoTree == 0) in_dir = configHelper->diLeptonic_dir + "ttmd_analysis_" + year + "/" + plots_version + "/" + i_ch + "/" + unfoldingIOHandler->Suffix;
                    else in_dir = configHelper->diLeptonic_dir + "ttmd_analysis_" + year + "/" + plots_version + "/" + i_ch + TString::Format("-kr%d/", ModeKinRecoTree) + "/" + unfoldingIOHandler->Suffix;
                    in_dir = in_dir + "/all/xsec" + undfolding_opt;
                    // std::cout << "LOADING: " << in_dir << std::endl;
                    std::pair<std::pair<TH1D*, TH1D*>, std::pair<TH1D*, std::pair<TH1D*, TH1D*>>> readXSec_with_systs = ReadDataXsecWithSystAndStats(in_dir, true);
                    TH1D* hDat = readXSec_with_systs.first.first;
                    TH1D* hGen = readXSec_with_systs.first.second;
                    hDat->SetTitle(i_ch);
                    hGen->SetTitle("POW+PYT (" + i_ch + ")");

                    hDat->SetMarkerStyle(vMarkers[ch_counter]);
                    hDat->SetMarkerColor(vColors[ch_counter]);
                    hDat->SetLineColor(vColors[ch_counter]);

                    hGen->SetLineColor(vColors[ch_counter]);

                    vector_hDat.push_back(hDat);
                    vector_hGen.push_back(hGen);
                    vector_TotUD.push_back(readXSec_with_systs.second.second);
                    ch_counter++;
                  }
                  FlagPaper = 1;
                  XSec(unfoldingIOHandler, vector_hDat, vector_TotUD, vector_hGen, out_dir + undfolding_opt + "-pap", XSecParams(unfoldingIOHandler), std::pair<TH1D*, TH1D*>(NULL, NULL), false);
                  FlagPaper = 2;
                  XSec(unfoldingIOHandler, vector_hDat, vector_TotUD, vector_hGen, out_dir + undfolding_opt + "-pas", XSecParams(unfoldingIOHandler), std::pair<TH1D*, TH1D*>(NULL, NULL), false);
            }

            // TString dir = "";
            // if(ModeKinRecoTree == 0) dir = configHelper->gAnalDir + "/" + "all_channel_comparison" + "/" + unfoldingIOHandler->Suffix;
            // else dir = configHelper->gAnalDir + "/" + "all_channel_comparison" + "/" + TString::Format("-kr%d", ModeKinRecoTree) + "-" + unfoldingIOHandler->Suffix;
            // for(auto it : mhAllSufXSec)
            // {
            //   const TString suf = it.first;
            //   TString fileName = dir + "/" + nameMeasDir + "/xsec-" + suf;
            //   std::vector<int> vMarkers = {24, 25, 26, 32};
            //   std::vector<int> vColors = {4, kOrange - 6, 2, 1};
            //   for(size_t c = 0; c < vCh.size(); c++)
            //   {
            //     mhAllSufXSec[suf][c]->SetMarkerStyle(vMarkers[c]);
            //     mhAllSufXSec[suf][c]->SetMarkerColor(vColors[c]);
            //     mhAllSufXSec[suf][c]->SetLineColor(vColors[c]);
            //     TString oldTitle = mhAllSufXSec[suf][c]->GetTitle();
            //     // replace 1st word, keep [chi2/dof]
            //     TString oldTitleChi2 = oldTitle.Data() + oldTitle.First(' ') + 1;
            //     mhAllSufXSec[suf][c]->SetTitle(vCh[c] + " " + oldTitleChi2);
            //   }
            //   //XSec(unfoldingIOHandler, mhAllSufXSec[suf], mhAllSufXSecUD[suf], mhAllSufXSecGen[suf], fileName, XSecParams(unfoldingIOHandler));
            //   XSec(unfoldingIOHandler, mhAllSufXSec[suf], mhAllSufXSecUD[suf], { mhAllSufXSecGen[suf][0] }, fileName, XSecParams(unfoldingIOHandler));
            //   // with stat. unc. only
            //   fileName = dir + "/" + nameMeasDir + "/xsec-statonly-" + suf;
            //   std::vector<std::pair<TH1D*, TH1D*> > uncEmpty;
            //   for(size_t c = 0; c < vCh.size(); c++)
            //   {
            //     TString oldTitle = mhAllSufXSecTitleChi2StatUncOnly[suf][c];
            //     // replace 1st word, keep [chi2/dof]
            //     TString oldTitleChi2 = oldTitle.Data() + oldTitle.First(' ') + 1;
            //     mhAllSufXSec[suf][c]->SetTitle(vCh[c] + " " + oldTitleChi2);
            //   }
            //   //XSec(unfoldingIOHandler, mhAllSufXSec[suf], uncEmpty, mhAllSufXSecGen[suf], fileName, XSecParams(unfoldingIOHandler));
            //   XSec(unfoldingIOHandler, mhAllSufXSec[suf], uncEmpty, { mhAllSufXSecGen[suf][0] }, fileName, XSecParams(unfoldingIOHandler));
            //   // check that gen. xsec is consistent in all channels (differences smaller than 1% or data unc.)
            //   for(int b = 1; b <= mhAllSufXSecGen[suf][0]->GetNbinsX(); b++)
            //   {
            //     assert(vCh.size() == 0 || mhAllSufXSecGen[suf][0]->GetBinContent(b) > 1e-10);
            //     for(size_t c = 1; c < vCh.size(); c++)
            //     {
            //       //double diff = (mhAllSufXSecGen[suf][c]->GetBinContent(b) - mhAllSufXSecGen[suf][0]->GetBinContent(b)) / mhAllSufXSecGen[suf][0]->GetBinContent(b);
            //       assert(diff < 1e-3 || diff < mhAllSufXSec[suf][0]->GetBinError(b));
            //     }
            //   }
            // }


	  }
	}
        
        // draw all year comparison with fullRun2
        if (flagAllYearCompar && year == "fullRun2" && unfoldingIOHandler->Suffix != "njmttyttdefIso04Pt30-b4-mtt3"){
            if(unfoldingIOHandler->DoXSec){
                std::map<TString, std::vector<TH1D*> > mhAllSufXSec_by_year;
                std::map<TString, std::vector<std::pair<TH1D*, TH1D*> > > mhAllSufXSecUD_by_year;
                std::map<TString, std::vector<TH1D*> > mhAllSufXSecGen_by_year;
                //std::map<TString, std::vector<TString> > mhAllSufXSecTitleChi2StatUncOnly_by_year;

                std::vector<int> vMarkers = {24, 25, 26, 32};
                std::vector<int> vColors = {4, 8, kOrange - 4, 2};

                std::vector<TString> years_list = configHelper->years_list;
                years_list.insert(years_list.begin(),"fullRun2");
                std::vector<TString> channels_to_compare = vCh;
                std::vector<TString> undfolding_opts = {"-nor-reg", "-abs-reg"};
                TString plots_version = configHelper->analysisSuffix;

                for (TString undfolding_opt: undfolding_opts){
                  for (TString i_ch: channels_to_compare){
                    TString out_dir = "";
                    if(ModeKinRecoTree == 0) out_dir = configHelper->gAnalDir + "/all_year_comparison" + "/" + i_ch + "-" + unfoldingIOHandler->Suffix;
                    else out_dir = configHelper->gAnalDir + "/all_year_comparison" + "/" + i_ch + TString::Format("-kr%d", ModeKinRecoTree) + "-" + unfoldingIOHandler->Suffix;
                    //std::cout << "SAVING in : " << out_dir << std::endl;

                    std::vector<TH1D*> vector_hDat;
                    std::vector<TH1D*> vector_hGen;
                    std::vector<std::pair<TH1D*, TH1D*> > vector_TotUD;

                    unsigned int year_counter = 0;
                    for (TString i_year: years_list){
                        TString in_dir = "";
                        if(ModeKinRecoTree == 0) in_dir = configHelper->diLeptonic_dir + "ttmd_analysis_" + i_year + "/" + plots_version + "/" + i_ch + "/" + unfoldingIOHandler->Suffix;
                        else in_dir = configHelper->diLeptonic_dir + "ttmd_analysis_" + i_year + "/" + plots_version + "/" + i_ch + TString::Format("-kr%d/", ModeKinRecoTree) + "/" + unfoldingIOHandler->Suffix;
                        in_dir = in_dir + "/all/xsec" + undfolding_opt;
                        // std::cout << "LOADING: " << in_dir << std::endl;
                        std::pair<std::pair<TH1D*, TH1D*>, std::pair<TH1D*, std::pair<TH1D*, TH1D*>>> readXSec_with_systs = ReadDataXsecWithSystAndStats(in_dir, true);
                        TH1D* hDat = readXSec_with_systs.first.first;
                        TH1D* hGen = readXSec_with_systs.first.second;
                        hDat->SetTitle(i_year);
                        hGen->SetTitle("POW+PYT (" + i_year + ")");

                        hDat->SetMarkerStyle(vMarkers[year_counter]);
                        hDat->SetMarkerColor(vColors[year_counter]);
                        hDat->SetLineColor(vColors[year_counter]);

                        hGen->SetLineColor(vColors[year_counter]);

                        vector_hDat.push_back(hDat);
                        vector_hGen.push_back(hGen);
                        vector_TotUD.push_back(readXSec_with_systs.second.second);
                        year_counter++;
                    }
                    FlagPaper = 1;
                    XSec(unfoldingIOHandler, vector_hDat, vector_TotUD, vector_hGen, out_dir + undfolding_opt + "-pap", XSecParams(unfoldingIOHandler));
                    FlagPaper = 2;
                    XSec(unfoldingIOHandler, vector_hDat, vector_TotUD, vector_hGen, out_dir + undfolding_opt + "-pas", XSecParams(unfoldingIOHandler));
                  }

                }

            }
        }

        // draw migration matrix years comparison with fullRun2
        if (flagAllYearMigMatrixSystsCompar && year == "fullRun2" && unfoldingIOHandler->DoXSec){
            std::vector<TString> years_list = configHelper->years_list;
            unsigned int N_years = years_list.size();
            years_list.insert(years_list.begin(),"fullRun2");
            std::vector<TString> channels_to_compare = vCh;

            bool doIndependentVarsPlots = true;
            bool doIndependentSystsPlots = false;
            bool doSummarySystsPlots = false;

            TString plots_version = configHelper->analysisSuffix;
            TString file_prefix = "/mig.txt";

            for (TString i_ch: channels_to_compare){
                TString out_dir = "";
                if(ModeKinRecoTree == 0) out_dir = configHelper->gAnalDir + "/all_year_mig_matrix_comparison" + "/" + i_ch + "-" + unfoldingIOHandler->Suffix;
                else out_dir = configHelper->gAnalDir + "/all_year_mig_matrix_comparison" + "/" + i_ch + TString::Format("-kr%d", ModeKinRecoTree) + "-" + unfoldingIOHandler->Suffix;
                //std::cout << "SAVING in : " << out_dir << std::endl;
                //        syst               year      matrix-plain
                std::map<TString, std::map<TString, std::vector<float>>> ByYearBySystMigMatricesPlain;
                //        syst    diff-matrix-plain
                std::map<TString, std::vector<float>> BySystDiffMeanPlain;
                //        syst    max-diff-matrix-plain
                std::map<TString, std::vector<float>> BySystDiffMaxPlain;
                //        syst          stat-name  stat-value
                std::map<TString, std::map<TString,float>> BySystStats;

                for (TString i_from_list_var: vars_to_fill){
                    std::map<TString, std::vector<float>> by_year_info, by_year_diff;
                    for (TString i_year: years_list){
                        TString in_base_dir = "";
                        TString i_var = i_from_list_var;
                        if (i_year == "2016" && (i_var.Contains("PSFSRSCALE") || i_var.Contains("PSISRSCALE"))) {
                            if (i_var.Contains("_UP")) i_var = i_var.ReplaceAll("_2_UP", "_UP");
                            else i_var = i_var.ReplaceAll("_2_DOWN", "_DOWN");
                        }
                        if(ModeKinRecoTree == 0) in_base_dir = configHelper->diLeptonic_dir + "ttmd_analysis_" + i_year + "/" + plots_version + "/" + i_ch + "/" + unfoldingIOHandler->Suffix;
                        else in_base_dir = configHelper->diLeptonic_dir + "ttmd_analysis_" + i_year + "/" + plots_version + "/" + i_ch + TString::Format("-kr%d/", ModeKinRecoTree) + "/" + unfoldingIOHandler->Suffix;
                        in_base_dir = in_base_dir + "/" + i_var + file_prefix;
                        //std::cout << "LOADING: " << in_base_dir << std::endl;
                        std::vector<float> i_plain_matrix = ReadMigrationMatrix(in_base_dir, true);
                        by_year_info.insert(std::pair<TString, std::vector<float>>(i_year,i_plain_matrix));
                    }
                    ByYearBySystMigMatricesPlain.insert(std::pair<TString, std::map<TString, std::vector<float>>> (i_from_list_var, by_year_info));

                    if (doIndependentVarsPlots) {
                        std::vector<float> DiffMeanPlainVector = {};
                        std::vector<float> DiffPlainVector = {};
                        std::vector<float> MaxDiffPlainVector = {};
                        std::vector<float> refVector = {};
                        auto it_ref = by_year_info.find("fullRun2");
                        if (it_ref != by_year_info.end()) refVector = it_ref->second;

                        for (auto i_year_info: by_year_info){
                            if (i_year_info.first != "fullRun2"){
                                for (unsigned int i_bin=0; i_bin<refVector.size(); i_bin++){
                                    float i_diff = TMath::Abs(i_year_info.second.at(i_bin) - refVector.at(i_bin));
                                    if (DiffMeanPlainVector.size() != refVector.size()) DiffMeanPlainVector.push_back(0.0);
                                    DiffMeanPlainVector.at(i_bin) += i_diff;
                                    if (MaxDiffPlainVector.size() != refVector.size()) MaxDiffPlainVector.push_back(-999.99);
                                    if (i_diff > MaxDiffPlainVector.at(i_bin)) MaxDiffPlainVector.at(i_bin) = i_diff;

                                    if (DiffPlainVector.size() != refVector.size()) DiffPlainVector.push_back(0.0);
                                    DiffPlainVector.at(i_bin) = i_diff;
                                }
                                by_year_diff.insert(std::pair<TString, std::vector<float>>(i_year_info.first,DiffPlainVector));
                                DiffPlainVector.clear();
                            }
                        }

                        //Computing mean
                        for (unsigned int i_bin=0; i_bin<DiffMeanPlainVector.size(); i_bin++) {
                            DiffMeanPlainVector.at(i_bin) = DiffMeanPlainVector.at(i_bin)/N_years;
                        }

                        // PlotMigrationMatrixFromPlainVector(DiffMeanPlainVector,unfoldingIOHandler,out_dir+"-"+i_from_list_var+"-mean_diff");
                        PlotMultipleMigrationMatrixFromPlainVectors(by_year_info,unfoldingIOHandler,out_dir+"-"+i_from_list_var+"-all_vals",refVector, "fullRun2 value (in %)");
                        PlotMultipleMigrationMatrixFromPlainVectors(by_year_diff,unfoldingIOHandler,out_dir+"-"+i_from_list_var+"-all_diff",MaxDiffPlainVector, " max. diff. with fullRun2 (in %)");

                    }

                    by_year_info.clear();
                }


                if (doSummarySystsPlots){
                    // Data for all systs
                    std::vector<float> AllSystDiffMeanPlainVector = {};
                    std::vector<float> AllSystMaxDiffPlainVector = {};
                    std::vector<int> AllSystMaxDiffYearVector = {};
                    unsigned int all_N_vectors = 0;

                    for (auto i_syst: allSys){
                        std::vector<TString> i_vars = configHelper->GetSysVars(i_syst);

                        std::vector<float> DiffMeanPlainVector = {};
                        std::vector<float> MaxDiffPlainVector = {};
                        std::vector<float> refVector = {};
                        unsigned int N_vectors = 0;

                        for (auto i_var: i_vars){
                            auto it_ref = ByYearBySystMigMatricesPlain.find(i_var);
                            if (it_ref == ByYearBySystMigMatricesPlain.end()) std::cerr << "ERROR in PlotterDifdXsec(migration matrices stduies) <--- syst variation not found in ByYearBySystMigMatricesPlain: " << i_var << std::endl;

                            auto it_ref_year = it_ref->second.find("fullRun2");
                            if (it_ref_year != it_ref->second.end()) refVector = it_ref_year->second;

                            for (auto i_year_info: it_ref->second){
                                if (i_year_info.first != "fullRun2"){
                                    for (unsigned int i_bin=0; i_bin<refVector.size(); i_bin++){
                                        float i_diff = TMath::Abs(i_year_info.second.at(i_bin) - refVector.at(i_bin));

                                        if (DiffMeanPlainVector.size() != refVector.size()) DiffMeanPlainVector.push_back(0.0);
                                        DiffMeanPlainVector.at(i_bin) += i_diff;

                                        if (MaxDiffPlainVector.size() != refVector.size()) MaxDiffPlainVector.push_back(-999.99);
                                        if (i_diff > MaxDiffPlainVector.at(i_bin)) MaxDiffPlainVector.at(i_bin) = i_diff;

                                        //Filling All systs vectors
                                        if (AllSystDiffMeanPlainVector.size() != refVector.size()) AllSystDiffMeanPlainVector.push_back(0.0);
                                        AllSystDiffMeanPlainVector.at(i_bin) += i_diff;

                                        if (AllSystMaxDiffPlainVector.size() != refVector.size()) {
                                            AllSystMaxDiffPlainVector.push_back(-999.99);
                                            AllSystMaxDiffYearVector.push_back(9999);
                                        }
                                        if (i_diff > AllSystMaxDiffPlainVector.at(i_bin)) {
                                            AllSystMaxDiffPlainVector.at(i_bin) = i_diff;
                                            AllSystMaxDiffYearVector.at(i_bin) = i_year_info.first.Atoi();
                                        }

                                    }
                                    N_vectors++;
                                    all_N_vectors++;
                                }
                            }
                        }


                        if (doIndependentSystsPlots){
                            //Computing mean
                            for (unsigned int i_bin=0; i_bin<DiffMeanPlainVector.size(); i_bin++) {
                                DiffMeanPlainVector.at(i_bin) = DiffMeanPlainVector.at(i_bin)/N_vectors;
                            }
                            PlotMigrationMatrixFromPlainVector(DiffMeanPlainVector,unfoldingIOHandler,out_dir+"-"+i_syst+"-mean_diff");
                            PlotMigrationMatrixFromPlainVector(MaxDiffPlainVector,unfoldingIOHandler,out_dir+"-"+i_syst+"-max_diff");

                        }

                    }

                    //Computing mean
                    for (unsigned int i_bin=0; i_bin<AllSystDiffMeanPlainVector.size(); i_bin++) {
                        AllSystDiffMeanPlainVector.at(i_bin) = AllSystDiffMeanPlainVector.at(i_bin)/all_N_vectors;
                    }
                    PlotMigrationMatrixFromPlainVector(AllSystDiffMeanPlainVector,unfoldingIOHandler,out_dir+"-allSysts-mean_diff");
                    //PlotMigrationMatrixFromPlainVector(AllSystMaxDiffPlainVector,unfoldingIOHandler,out_dir+"-allSysts-max_diff", AllSystMaxDiffYearVector);
                    PlotMigrationMatrixFromPlainVector(AllSystMaxDiffPlainVector,unfoldingIOHandler,out_dir+"-allSysts-max_diff");
                }
            }
        }

      }
      
      if(mg)
        delete mg;
    }


void PlotterDiffXSec::PlotMGTheory(ZMadGraph* mg, const TString fileName, const std::vector<ZPredSetup>& vPredSetup, const std::vector<TString> vTitle,
                      const int dof, const TString genTitle, const UnfoldingIOHandler* unfoldingIOHandler, ZMeasur* measXSec, std::pair<TH1D*, TH1D*> hXSecToPlotUD, TH1D* hXSecGenNom,
			       const bool flagPlotPaper)
    {
      UnfoldingXSecHandler* unfoldingXSecHandler = new UnfoldingXSecHandler();
      assert(vPredSetup.size() == vTitle.size());

      TH1D* hXSecToPlot = (TH1D*) measXSec->GetNom()->Clone();
      hXSecToPlot->SetMarkerStyle(24);

      //hXSecToPlot->SetMarkerStyle(24);
      hXSecToPlot->SetMarkerColor(1);
      hXSecToPlot->SetLineColor(1);
      hXSecToPlot->SetTitle(TString::Format("Data, dof=%d", dof));

      std::vector<TH1D*> vhXSecGen;
      //vhXSecGen.push_back(hXSecGenNom);

      bool flagNorm = ( (hXSecGenNom->GetNbinsX() - 1) == dof );
      for(unsigned int i = 0; i < vPredSetup.size(); i++)
      {
        //TH1D* h = unfoldingXSecHandler->GetPredictionHisto(unfoldingIOHandler, mg, vPredSetup[i], flagNorm, flagNorm);
        TH1D* h = unfoldingXSecHandler->GetPredictionHisto(unfoldingIOHandler, mg, vPredSetup[i], flagNorm, false);
        if(dof > 0)
        {
          double chi2 = 0.0;
          if(unfoldingIOHandler->DoRj)
            chi2 = utils::ttmd::Chi2(hXSecToPlot, measXSec->GetCovFullMatrix(), h, unfoldingIOHandler->VSkip);
          else
            chi2 = utils::ttmd::Chi2(hXSecToPlot, measXSec->GetCovFullMatrix(), h, flagNorm);
          h->SetTitle(TString::Format("%s [%.0f/%d]", vTitle[i].Data(), chi2, dof));
        }
        h->SetLineColor(utils::ttmd::GetDistinctColor(i));
        //h->SetLineStyle(i + 1);
        if(flagPlotPaper)
          h->SetLineStyle(utils::ttmd::GetDistinctLineStyle(i));
        h->SetLineWidth(2);
        vhXSecGen.push_back(h);
      }

      TH1D* hXSecGen0 = (TH1D*) hXSecGenNom->Clone();
      hXSecGen0->SetTitle(genTitle);
      hXSecGen0->SetLineColor(0);
      hXSecGen0->Scale(1e8);
      //vhXSecGen.push_back(hXSecGen0);
      vhXSecGen.insert(vhXSecGen.begin(), hXSecGen0);

      XSec(unfoldingIOHandler, { hXSecToPlot }, { hXSecToPlotUD }, vhXSecGen, fileName, XSecParams(unfoldingIOHandler));
      //printf("integral = %f\n", hXSecToPlot->Integral());
      //printf("integral = %f\n", vhXSecGen[0]->Integral());
    }

    // (with several sets of data for mt dependence)
void PlotterDiffXSec::PlotMGTheory(ZMadGraph* mg, const TString fileName, const std::vector<ZPredSetup>& vPredSetup, const std::vector<TString> vTitle,
                      const int dof, const TString genTitle, const UnfoldingIOHandler* unfoldingIOHandler, std::vector<TH1D*> vhData, const TH2D* hCov, std::pair<TH1D*, TH1D*> hXSecToPlotUD, TH1D* hXSecGenNom, const bool flagPlotPaper)
    {

      UnfoldingXSecHandler* unfoldingXSecHandler = new UnfoldingXSecHandler();

      assert(vPredSetup.size() == vTitle.size());
      assert(vhData.size() >= vPredSetup.size());

      std::vector<TH1D*> vhXSecGen;
      //vhXSecGen.push_back(hXSecGenNom);

      bool flagNorm = ( (hXSecGenNom->GetNbinsX() - 1) == dof );
      for(unsigned int i = 0; i < vPredSetup.size(); i++)
      {
        //TH1D* h = unfoldingXSecHandler->GetPredictionHisto(unfoldingIOHandler, mg, vPredSetup[i], flagNorm, flagNorm);
        TH1D* h = unfoldingXSecHandler->GetPredictionHisto(unfoldingIOHandler, mg, vPredSetup[i], flagNorm, false);
        if(dof > 0)
        {
          double chi2 = 0.0;
          if(unfoldingIOHandler->DoRj)
            chi2 = utils::ttmd::Chi2(vhData[i], hCov, h, unfoldingIOHandler->VSkip);
          else
            chi2 = utils::ttmd::Chi2(vhData[i], hCov, h, true);
          if(flagPlotPaper > 0 && fileName.Contains("mtdep"))
            h->SetTitle(TString::Format("%s", vTitle[i].Data()));
          else
            h->SetTitle(TString::Format("%s [%.0f/%d]", vTitle[i].Data(), chi2, dof));
        }
        h->SetLineColor(utils::ttmd::GetDistinctColor(i));
        //h->SetLineStyle(i + 1);
        if(flagPlotPaper)
          h->SetLineStyle(utils::ttmd::GetDistinctLineStyle(i));
        h->SetLineWidth(2);
        vhXSecGen.push_back(h);
        vhData[i]->SetMarkerStyle(utils::ttmd::GetDistinctMarkerStyle(i));
        vhData[i]->SetMarkerColor(h->GetLineColor());
        vhData[i]->SetLineColor(h->GetLineColor());
      }

      TH1D* hXSecGen0 = (TH1D*) hXSecGenNom->Clone();
      hXSecGen0->SetTitle(genTitle);
      hXSecGen0->SetLineColor(0);
      hXSecGen0->Scale(1e8);
      //vhXSecGen.push_back(hXSecGen0);
      vhXSecGen.insert(vhXSecGen.begin(), hXSecGen0);

      std::vector<std::pair<TH1D*, TH1D*> > vhXSecToPlotUD(vhData.size(), hXSecToPlotUD);

      XSec(unfoldingIOHandler, vhData, vhXSecToPlotUD, vhXSecGen, fileName, XSecParams(unfoldingIOHandler));
      //printf("integral = %f\n", hXSecToPlot->Integral());
      //printf("integral = %f\n", vhXSecGen[0]->Integral());
    }

    // simplified version with theory histograms provided directly
void PlotterDiffXSec::PlotMGTheory(const TString fileName, const std::vector<TH1D*>& vPredH,
                      const TString genTitle, const UnfoldingIOHandler* unfoldingIOHandler, ZMeasur* measXSec, std::pair<TH1D*, TH1D*> hXSecToPlotUD, TH1D* hXSecGenNom)
    {
      TH1D* hXSecToPlot = (TH1D*) measXSec->GetNom()->Clone();
      hXSecToPlot->SetMarkerStyle(24);

      //hXSecToPlot->SetMarkerStyle(24);
      hXSecToPlot->SetMarkerColor(1);
      hXSecToPlot->SetLineColor(1);
      TString str = vPredH.back()->GetTitle();
      TString str1 = str.Data() + str.First('/') + 1;
      int dof = TString(str1, str1.First(']')).Atoi();
      hXSecToPlot->SetTitle(TString::Format("Data, dof=%d", dof));

      std::vector<TH1D*> vhXSecGen;
      //vhXSecGen.push_back(hXSecGenNom);

      for(unsigned int i = 0; i < vPredH.size(); i++)
      {
        TH1D* h = (TH1D*) vPredH[i]->Clone();
        h->SetLineColor(utils::ttmd::GetDistinctColor(i));
        //h->SetLineStyle(i + 1);
        h->SetLineWidth(2);
        vhXSecGen.push_back(h);
      }

      TH1D* hXSecGen0 = (TH1D*) hXSecGenNom->Clone();
      hXSecGen0->SetTitle(genTitle);
      hXSecGen0->SetLineColor(0);
      hXSecGen0->Scale(1e8);
      //vhXSecGen.push_back(hXSecGen0);
      vhXSecGen.insert(vhXSecGen.begin(), hXSecGen0);

      XSec(unfoldingIOHandler, { hXSecToPlot }, { hXSecToPlotUD }, vhXSecGen, fileName, XSecParams(unfoldingIOHandler));
      //printf("integral = %f\n", hXSecToPlot->Integral());
      //printf("integral = %f\n", vhXSecGen[0]->Integral());
    }

ZMeasur* PlotterDiffXSec::FillMeas(const std::map<TString, TH1D*>& mapVar, const std::vector<TString>& selectVar, const TH2D* hStatCov, const bool flagShapeOnly)
    {
      ZMeasur* meas = new ZMeasur(configHelper);
      if(hStatCov)
        meas->DoEigenvectors = true;
      meas->FlagShapeOnly = flagShapeOnly;

      for(std::vector<TString>::const_iterator it = selectVar.begin(); it != selectVar.end(); it++)
      {
        const auto& sys = *it;
        //if(configHelper->IsSkipped(sys, 1))
        //  continue;
        for(const auto& var : configHelper->GetSysVars(sys))
        {
          auto it = mapVar.find(var);
          if (it == mapVar.end()) {
              std::cerr << "ERROR in PlotterDiffXSec::FillMeas(): This var is not included in the provided map <-- " << var << " (comming from: " << sys << ")" << std::endl;
              exit(1);
          }
          if(var == "Nominal")
            meas->SetNom(mapVar.at(var));
          else
            meas->AddVar(var, mapVar.at(var));
        }
      }

      if(hStatCov)
        meas->SetCovStatMatrix(utils::ttmd::MakeCovMatrixFromCorr(hStatCov, meas->GetNom()));
      meas->Process();

      return meas;
    }

std::vector<TH1D*> PlotterDiffXSec::GetVectorUncUD(const ZMeasur& meas, const TString& type)
    {
      std::vector<TH1D*> vhUnc;
      if(type == "total")
      {
        vhUnc.push_back((TH1D*) meas.GetTotalUD().first->Clone());
        vhUnc.push_back((TH1D*) meas.GetTotalUD().second->Clone());
        vhUnc[0]->Divide(meas.GetNom());
        vhUnc[1]->Divide(meas.GetNom());
        vhUnc[0]->Scale(100.0);
        vhUnc[1]->Scale(100.0);
      }
      else if(type == "stat")
      {
        vhUnc.push_back((TH1D*) meas.GetSystUD().first->Clone());
        vhUnc.push_back((TH1D*) meas.GetSystUD().second->Clone());
        for(int b = 0; b < vhUnc[0]->GetNbinsX(); b++)
        {
          vhUnc[0]->SetBinContent(b + 1, 100 * TMath::Sqrt(meas.GetCovStatMatrix()->GetBinContent(b + 1, b + 1)));
          vhUnc[1]->SetBinContent(b + 1, -100 *  TMath::Sqrt(meas.GetCovStatMatrix()->GetBinContent(b + 1, b + 1)));
        }
        vhUnc[0]->Divide(meas.GetNom());
        vhUnc[1]->Divide(meas.GetNom());
      }
      else
        //else if(type == "syst")
      {
        vhUnc.push_back((TH1D*) meas.GetSystUD().first->Clone());
        vhUnc.push_back((TH1D*) meas.GetSystUD().second->Clone());
        vhUnc[0]->Divide(meas.GetNom());
        vhUnc[1]->Divide(meas.GetNom());
        vhUnc[0]->Scale(100.0);
        vhUnc[1]->Scale(100.0);
      }
      //else
      //  throw std::logic_error(TString::Format("Error in GetVectorUncUD(): unknown type = %s", type.Data()));
      double dev = utils::ttmd::GetMeanDeviationFromRef(vhUnc[0], vhUnc[1], NULL);
      vhUnc[0]->SetTitle(TString::Format("%s [~%.1f%%]", type.Data(), dev));
      return vhUnc;
    }

std::vector<TH1D*> PlotterDiffXSec::PlotJESSubtotalUnc(const UnfoldingIOHandler* unfoldingIOHandler, const std::map<TString, TH1D*>& vmhVarDat, const TString& nameSubtotal,
                                          const std::vector<TString>& vJESsrc, std::vector<TH1D*> vhUncUnf,
					  const TString& fileNameBase, const UncSummaryParams& pars, const bool flagPlotSrc)
    {
      std::vector<int> vMarkerStyle = { 0, 24, 25, 26, 32, 27, 28, 30, 31, 2, 3, 33, 34 }; // 0 is fake for 0th entry (Nominal)

      std::vector<TString> vJES = vJESsrc;
      // add prefix "JES"
      for(unsigned int v = 0; v < vJES.size(); v++)
        vJES[v] = "JES" + vJES[v];
      // insert reference
      vJES.insert(vJES.begin(), "Nominal");

      std::vector<TH1D*> hJESSumOrig = GetVectorUncUD(*FillMeas(vmhVarDat, { "Nominal", "JES" + nameSubtotal }), nameSubtotal + "-src.");
      std::vector<TH1D*> hJESSumCheck = GetVectorUncUD(*FillMeas(vmhVarDat, vJES), nameSubtotal);
      std::vector<std::vector<TH1D*> > vhUncJESSubTotalPileUp;
      for(unsigned int v = 1; v < vJES.size(); v++)
      {
        // skip prefix "JES"
        TString name = TString(vJES[v].Data() + 3);
        vhUncJESSubTotalPileUp.push_back(GetVectorUncUD(*FillMeas(vmhVarDat, {"Nominal", vJES[v]}), name));
        vhUncJESSubTotalPileUp.back()[0]->SetMarkerStyle(vMarkerStyle[v]);
      }
      hJESSumOrig.insert(hJESSumOrig.end(), vhUncUnf.begin(), vhUncUnf.end());
      TString fileName = fileNameBase;
      fileName.ReplaceAll("JESXXX", "JES" + nameSubtotal);
      fileName.ReplaceAll("JES", "jes-");
      if(flagPlotSrc)
      {
        PlotUncSummary(unfoldingIOHandler, vhUncJESSubTotalPileUp, hJESSumCheck, hJESSumOrig, fileName, pars);
      }
      else
      {
        PlotUncSummary(unfoldingIOHandler, {}, hJESSumCheck, hJESSumOrig, fileName, pars);
      }

      return hJESSumCheck;
    }

void PlotterDiffXSec::CoverTest(const TString& ch, const int type, const int niter)
    {
      const bool flagWrite = 1;
      const bool flagRead = 1;

      std::vector<ZUnfold*> vUnf = CreateZUnfoldVector();
      this->Var = "Nominal";

      const int seed = 19;
      if(type == 0)
      {
        const TString suffix = "covertest-pseudo";
        //const int niter = 5;
        //const int niter = 1500;
        //const int niter = 150;
        //const int niter = 1000;
        assert(niter >= 0);

        const int flagRew = 0;
        if(flagWrite)
          CoverTestWrite(vUnf, niter, suffix, flagRew, seed, ch);
        if(flagRead)
          CoverTestRead(vUnf, niter, suffix, ch);
      }
      else
      {
        const int niter = 0;
        const int flagRew = type;
        const TString suffix = TString::Format("covertest-rew%d", flagRew);
        if(flagWrite)
          CoverTestWrite(vUnf, niter, suffix, flagRew, seed, ch);
        if(flagRead)
          CoverTestRead(vUnf, niter, suffix, ch);
      }

      // delete unfolding
      for(unsigned int u = 0; u < vUnf.size(); u++)
        delete vUnf[u];
    }

void PlotterDiffXSec::CoverTestWrite(std::vector<ZUnfold*>& vUnf, const int niter, const TString suffix, const int flagRew, const int seed, const TString& ch)
    {
      // TODO remove duplication (see Analyse())
      //if(ModeKinRecoTree == 0)
      //  configHelper->gPlainTreeNamePlain = "plainTree_rec_step8";
      //else if(ModeKinRecoTree >= 1 || ModeKinRecoTree <= 9)
      //	configHelper->gPlainTreeNamePlain = "plainTree_recSimple_step8";
      //  //configHelper->gPlainTreeNamePlain = "plainTree_recSimple_step7L";
      //else
      //  throw std::logic_error(TString::Format("Error: unsupported ModeKinRecoTree = %d", ModeKinRecoTree).Data());
      //const TString ch = "em";
      //MaxEvents = 30000;
      LoopDataTreePlain(ch);
      LoopMCTreePlain(ch);

      for(unsigned int x = 0; x < vXSec.size(); x++)
      {
        UnfoldingIOHandler* unfoldingIOHandler = vXSec[x];
        TH1D* hGenRew = (TH1D*) unfoldingIOHandler->HGenC()->Clone();
        //hGenRef->Print("all");

        // reweighting function
        if(flagRew > 0)
        {
          // reweighting function
          TF1* fRew = UnfoldingXSecHandler::GetRewFunction(unfoldingIOHandler, flagRew);
          unfoldingIOHandler->SetReweightFunction(fRew);

          // re-normalise reweighting function
          TH1D* hGenRefCopy = (TH1D*) hGenRew->Clone();
          const double norm0 = hGenRew->Integral();
          //printf(">>>>> reweighting 1\n");
          StandaloneLoopMCTreeGenPlain(ch, 0, hGenRew, unfoldingIOHandler);
          //hGenRef->Print("all");
          const double norm1 = hGenRew->Integral();
          // set new normalisation parameter to get the same total x-section
          printf("old norm par: %f\n", fRew->GetParameter(2));
          fRew->SetParameter(2, fRew->GetParameter(2) * norm0 / norm1);
          printf("new norm par: %f\n", fRew->GetParameter(2));
          //printf(">>>>> reweighting 2\n");
          StandaloneLoopMCTreeGenPlain(ch, 0, hGenRew, unfoldingIOHandler);
          //hGenRef->Print("all");
          const double norm2 = hGenRew->Integral();
          printf("norm0 = %f norm1 = %f norm2 = %f\n", norm0, norm1, norm2);

          printf("hGenRef before reweighting\n");
          hGenRefCopy->Print("all");
          hGenRefCopy->Divide(hGenRew);
          printf("hGenRef after reweighting\n");
          hGenRew->Print("all");
          printf("Reweighting ratio\n");
          hGenRefCopy->Print("all");
        }

        // generate pseudo data
        UnfoldingClosureTestsHandler unfoldingClosureTestsHandler(unfoldingIOHandler, vUnf.size(), niter, configHelper);
        // TString fileName = UnfoldingClosureTestsHandler::GetRootFileNameToStore(unfoldingIOHandler, suffix, niter, ch, ModeKinRecoTree);
        TString fileName = unfoldingClosureTestsHandler.GetRootFileNameToStore(unfoldingIOHandler, suffix, niter, ch, ModeKinRecoTree);
        unfoldingClosureTestsHandler.WriteSample(fileName, vUnf, unfoldingIOHandler, hGenRew, ch, seed);
        delete hGenRew;
      }
    }

void PlotterDiffXSec::CoverTestRead(std::vector<ZUnfold*>& vUnf, const int niter, const TString suffix, const TString& ch)
    {
      for(unsigned int x = 0; x < vXSec.size(); x++)
      {
        UnfoldingIOHandler* unfoldingIOHandler = vXSec[x];
        UnfoldingClosureTestsHandler unfoldingClosureTestsHandler(unfoldingIOHandler, vUnf.size(), niter, configHelper);
        // TString fileName = UnfoldingClosureTestsHandler::GetRootFileNameToStore(unfoldingIOHandler, suffix, niter, ch, ModeKinRecoTree);
        TString fileName = unfoldingClosureTestsHandler.GetRootFileNameToStore(unfoldingIOHandler, suffix, niter, ch, ModeKinRecoTree);
        TString chDir = ch;
        if(ModeKinRecoTree)
          chDir += TString::Format("-kr%d/", ModeKinRecoTree);
        configHelper->outDir = configHelper->gAnalDir + "/" + chDir + "/" + suffix + "/" + unfoldingIOHandler->Suffix;
        unfoldingClosureTestsHandler.ReadSample(fileName, vUnf, configHelper->outDir);
        configHelper->outDir = "";
      }
    }

void PlotterDiffXSec::PrintXsecAllAlg(std::vector<ZUnfold*> vAlg, const TString& fileName)
{
  utils::ttmd::CheckFile(fileName);
  FILE* fout = fopen(fileName.Data(), "w");
  PrintXsecAllAlg(vAlg, fout, 0);
  fprintf(fout, "************************************************************************\n");
  PrintXsecAllAlg(vAlg, fout, 1);
  fclose(fout);
}

void PlotterDiffXSec::PrintXsecAllAlg(std::vector<ZUnfold*> vAlg, FILE* fout, const int typeXSec)
{
  if(typeXSec != 0 && typeXSec != 1)
    throw std::logic_error(TString::Format("Error in PrintXsecAllAlg(): unsupported typeXSec = %d", typeXSec).Data());

  if(vAlg.size() == 0)
    throw std::logic_error("Error in PrintXsecAllAlg(): empty vAlg");

  if(typeXSec == 0)
    fprintf(fout, " *** ABSOLUTE X-SECTION *** \n");
  else if(typeXSec == 1)
    fprintf(fout, " *** NORMALISED X-SECTION *** \n");

  int nbins = vAlg[0]->hUnf->GetNbinsX();
  if(typeXSec == 0)
  {
    fprintf(fout, "INTEGRATED XSEC ");
    for(unsigned int a = 0; a < vAlg.size(); a++)
    {
      TH1D* h = (typeXSec == 0) ? vAlg[a]->hUnf : vAlg[a]->hUnfN;
      if(!h)
        continue; // this happens for total xsec not unfolded
      double xsec = -1.0;
      double unc = -1.0;
      xsec = h->IntegralAndError(1, nbins, unc);
      fprintf(fout, "%s = %9.5e [%.1f%%] ", vAlg[a]->Suffix.Data(), xsec, unc / xsec * 100);
    }
    fprintf(fout, "\n");
  }

  for(int b = 0; b < nbins; b++)
  {
    double eff = vAlg[0]->hEffBBB->GetBinContent(b + 1) * 100.0;
    fprintf(fout, "bin [%2d]  EFF = %5.1f%%  XSEC ", b + 1, eff);
    for(unsigned int a = 0; a < vAlg.size(); a++)
    {
      TH1D* h = (typeXSec == 0) ? vAlg[a]->hUnf : vAlg[a]->hUnfN;
      if(!h)
        continue; // this happens for total xsec not unfolded
      double cs = h->GetBinContent(b + 1);
      double unc = h->GetBinError(b + 1) / cs * 100.0;
      fprintf(fout, "%s = %9.5e [%.1f%%] ", vAlg[a]->Suffix.Data(), cs, unc);
    }
    fprintf(fout, "\n");
  }
}

void PlotterDiffXSec::PrintDataXsec(const TH1D* hD, const TH1D* vT, const TString& fileName)
{
  //printf("PrintDataXsec fileName = %s\n", fileName.Data());
  utils::ttmd::CheckFile(fileName);
  FILE* f = fopen(fileName.Data(), "w");
  PrintHistBins(hD, f);
  fprintf(f, "%8s%12s%8s%12s\n", "hypbin", "xsec", "stat", "genMC");
  for(int b = 0; b < hD->GetNbinsX(); b++)
  {
    fprintf(f, "%8d", b);
    fprintf(f, "%12.*e", NDigFloatXSec, hD->GetBinContent(b + 1));
    fprintf(f, "%8.*f", NDigFloatStat, hD->GetBinError(b + 1) / hD->GetBinContent(b + 1) * 100.0);
    fprintf(f, "%12.*e", NDigFloatXSec, vT->GetBinContent(b + 1));
    fprintf(f, "\n");
  }
  fclose(f);
}

void PlotterDiffXSec::PrintDataXsec(const TH1D* hD, const std::pair<TH1D*, TH1D*> hud, const TH1D* vT, const TString& fileName)
{
  //printf("PrintDataXsec fileName = %s\n", fileName.Data());
  utils::ttmd::CheckFile(fileName);
  FILE* f = fopen(fileName.Data(), "w");
  PrintHistBins(hD, f);
  fprintf(f, "%8s%12s%8s%8s%8s%12s\n", "hypbin", "xsec", "stat", "systU", "systD", "genMC");
  for(int b = 0; b < hD->GetNbinsX(); b++)
  {
    fprintf(f, "%8d", b);
    fprintf(f, "%12.*e", NDigFloatXSec, hD->GetBinContent(b + 1));
    fprintf(f, "%8.*f", NDigFloatStat, hD->GetBinError(b + 1) / hD->GetBinContent(b + 1) * 100.0);
    fprintf(f, "%8.*f", NDigFloatSyst, hud.first->GetBinContent(b + 1) / hD->GetBinContent(b + 1) * 100.0);
    fprintf(f, "%8.*f", NDigFloatSyst, -1 * hud.second->GetBinContent(b + 1) / hD->GetBinContent(b + 1) * 100.0);
    fprintf(f, "%12.*e", NDigFloatXSec, vT->GetBinContent(b + 1));
    fprintf(f, "\n");
  }
  fclose(f);
}

void PlotterDiffXSec::PrintTexBinning(const UnfoldingIOHandler* unfoldingIOHandler, const TString& fileName)
{
  utils::ttmd::CheckFile(fileName);
  FILE* f = fopen(fileName.Data(), "w");
  const int dim = unfoldingIOHandler->Dim();
  const TH1D* hNom = unfoldingIOHandler->GetPlainHisto();

  // header
  fprintf(f, "\\begin{tabular}{|");
  for(int d = 0; d < dim; d++)
    fprintf(f, "c|");
  fprintf(f, "c|}\n");
  fprintf(f, "\\hline\n");
  for(int d = 0; d < dim; d++)
  {
    TString title = utils::ttmd::TexAdoptVarTitle(unfoldingIOHandler->Var(d)->Title);
    fprintf(f, "$%s$ ", title.Data());
    if(unfoldingIOHandler->Var(d)->Units != "")
      fprintf(f, "[%s${}^{-1}$] ", unfoldingIOHandler->Var(d)->Units.Data());
    fprintf(f, " & ");
  }
  fprintf(f, "bin\\\\\n");
  fprintf(f, "\\hline\n");

  // main part
  for(int b = 0; b < hNom->GetNbinsX(); b++)
  {
    int bin = b;
    int nbins = hNom->GetNbinsX();
    for(int d = 0; d < unfoldingIOHandler->Dim(); d++)
    {
      int nbinsDim = unfoldingIOHandler->Var(d)->BinsC().size() - 1;
      int norm = nbins / nbinsDim;
      int localBin = bin / norm;
      nbins /= nbinsDim;
      bin -= localBin * norm;
      int ndig = unfoldingIOHandler->Var(d)->NDigits;
      const double varmin = unfoldingIOHandler->Var(d)->BinsC()[localBin];
      const double varmax = unfoldingIOHandler->Var(d)->BinsC()[localBin + 1];
      fprintf(f, "$%10.*f$--$%10.*f$ & ", ndig, varmin, ndig, varmax);
    }
    fprintf(f, "%8d \\\\\n", b + 1);
    /*double bw = 1.0;
    for(const auto& h : vHNom)
      fprintf(f, "%12.*e", NDigFloatXSec, h->GetBinContent(b + 1) / bw);
    if(!flagPlain)
    {
      fprintf(f, "%8.*f", NDigFloatStat, hNom->GetBinError(b + 1) / hNom->GetBinContent(b + 1) * 100.0);
      if(hUD.first)
        fprintf(f, "%8.*f", NDigFloatSyst, hUD.first->GetBinContent(b + 1) / hNom->GetBinContent(b + 1) * 100.0);
      if(hUD.second)
        fprintf(f, "%8.*f", NDigFloatSyst, -1 * hUD.second->GetBinContent(b + 1) / hNom->GetBinContent(b + 1) * 100.0);
      for(std::map<std::pair<TString, TString>, TH1D*>::const_iterator it = mapVars.begin(); it != mapVars.end(); it++)
        fprintf(f, "%8.*f", NDigFloatSyst, it->second->GetBinContent(b + 1) / hNom->GetBinContent(b + 1) * 100);
    }
    fprintf(f, "\n");*/
  }

  // footer
  fprintf(f, "\\hline\n");
  fprintf(f, "\\end{tabular}\n");
  fclose(f);
}

void PlotterDiffXSec::PrintTexTabXSec(const UnfoldingIOHandler* unfoldingIOHandler, const ZMeasur* meas, const TString& fileName)
{
  assert(unfoldingIOHandler->Dim() == 2 || unfoldingIOHandler->Dim() == 3);
  const bool flagLessLines = 1;
  utils::ttmd::CheckFile(fileName);
  FILE* f = fopen(fileName.Data(), "w");
  const int dim = unfoldingIOHandler->Dim();
  const TH1D* hNom = unfoldingIOHandler->GetPlainHisto();
  TString strExtraBins;
  TString strXSec = utils::ttmd::TexXsecTitle(unfoldingIOHandler->Suffix, strExtraBins);
  strXSec.ReplaceAll("p_{T}", "p_{\rm T}");
  TString strXSecLabel = utils::ttmd::TexXsecLabel(strXSec, unfoldingIOHandler->Var(0)->Bins[0].size() - 1);

  // NP
  TH1D* hNP = NULL;
  if(unfoldingIOHandler->FlagNeedNP)
  {
    TString fileNameNP = TString(fileName.Data(), fileName.Last('/')) + "/../../../npcorr/" + unfoldingIOHandler->Suffix + "/" + "npcorr.dat";
    hNP = unfoldingIOHandler->GetPlainHisto();
    utils::ttmd::ReadHistoKFactor(hNP, fileNameNP, unfoldingIOHandler->Dim() * 2 + 1 + 1, 1);
    //const int nbins = unfoldingIOHandler->HGenC()->GetNbinsX();
    assert(hNP->GetNbinsX() == nbins);
    //hNP->Print("all");
  }

  // header
  fprintf(f, "\\begin{table}[b]\n");
  fprintf(f, "\\begin{center}\n");
  TString strExtraNP;
  if(unfoldingIOHandler->FlagNeedNP)
    strExtraNP = ", and NP corrections (see Section~\\ref{sec:asmt})";
  fprintf(f, "\\caption{The measured %s cross sections%s, along with their relative statistical and systematic uncertainties%s.}\n", strXSec.Data(), strExtraBins.Data(), strExtraNP.Data());
  fprintf(f, "\\label{tab:%s_xsec}\n", strXSecLabel.Data());
  if(hNom->GetNbinsX() == 36)
    fprintf(f, "\\footnotesize\n");
  else
    fprintf(f, "\\footnotesize\n");
  fprintf(f, "\\renewcommand*{\\arraystretch}{1.55}\n");
  fprintf(f, "\\tabcolsep2.0mm\n");
  fprintf(f, "\\begin{tabular}{");
  if(!flagLessLines)
    fprintf(f, "|");
  for(int d = 0; d < dim; d++)
  {
    if(flagLessLines)
      fprintf(f, "c");
    else
      fprintf(f, "c|");
  }
  if(flagLessLines)
    fprintf(f, "cccc");
  else
    fprintf(f, "c|c|c|c|");
  if(unfoldingIOHandler->FlagNeedNP)
  {
    if(flagLessLines)
      fprintf(f, "c");
    else
      fprintf(f, "c|");
  }
  fprintf(f, "}\n");
  if(!flagLessLines)
    fprintf(f, "\\hline\n");
  for(int d = 0; d < dim; d++)
  {
    TString title = utils::ttmd::TexAdoptVarTitle(unfoldingIOHandler->Var(d)->Title);
    fprintf(f, "$%s$ ", title.Data());
    if(unfoldingIOHandler->Var(d)->Units != "")
      fprintf(f, "[{\\rm %s}] ", unfoldingIOHandler->Var(d)->Units.Data());
    fprintf(f, " & ");
  }
  TString title = utils::ttmd::TexAdoptVarYAxisTitle(unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - 1)->GetYTitleXSec());
  fprintf(f, "$%s$ & ", title.Data());
  fprintf(f, "stat. [\\%%] & ");
  fprintf(f, "syst. [\\%%] & ");
  if(unfoldingIOHandler->FlagNeedNP)
    fprintf(f, "{$\\cal{C_{\\text{NP}}}$} & ");
  fprintf(f, "bin\\\\\n");
  fprintf(f, "\\hline\n");

  // main part
  for(int b = 0; b < hNom->GetNbinsX(); b++)
  {
    int bin = b;
    int nbins = hNom->GetNbinsX();
    for(int d = 0; d < unfoldingIOHandler->Dim(); d++)
    {
      int nbinsDim = unfoldingIOHandler->Var(d)->BinsC().size() - 1;
      int norm = nbins / nbinsDim;
      int localBin = bin / norm;
      nbins /= nbinsDim;
      bin -= localBin * norm;
      int ndig = unfoldingIOHandler->Var(d)->NDigits;
      const double varmin = unfoldingIOHandler->Var(d)->BinsC()[localBin];
      const double varmax = unfoldingIOHandler->Var(d)->BinsC()[localBin + 1];
      char charBinLabel[256];
      sprintf(charBinLabel, "$%.*f$--$%.*f$ & ", ndig, varmin, ndig, varmax);
      TString strBinLabel = utils::ttmd::TexAdoptBinLabel(TString(charBinLabel));
      fprintf(f, "%20s", strBinLabel.Data());
    }
    double xsecVal = meas->GetNom()->GetBinContent(b + 1);
    int nbinsLastDim = unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - 1)->Bins[0].size() - 1;
    double binWidth = unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - 1)->Bins[0][b % nbinsLastDim + 1] - unfoldingIOHandler->Var(unfoldingIOHandler->Dim() - 1)->Bins[0][b % nbinsLastDim];
    double xsecValBW = xsecVal / binWidth;
    fprintf(f, "$%s$ & ", utils::ttmd::GetDigitTimesFormat(xsecValBW, 3).Data());
    double stat = meas->GetNom()->GetBinError(b + 1) / xsecVal * 100.0;
    fprintf(f, "$%5.1f$ & ", stat);
    double systU = meas->GetSystUD().first->GetBinContent(b + 1) / xsecVal * 100.0;
    double systD = meas->GetSystUD().second->GetBinContent(b + 1) / xsecVal * 100.0;
    fprintf(f, "${}_{%+5.1f}^{%+5.1f}$ & ", systD, systU);
    if(unfoldingIOHandler->FlagNeedNP)
      fprintf(f, "$%.3f$ & ", hNP->GetBinContent(b + 1));
    fprintf(f, "$%8d$ \\\\\n", b + 1);
  }

  // footer
  if(!flagLessLines)
    fprintf(f, "\\hline\n");
  fprintf(f, "\\end{tabular}\n");
  fprintf(f, "\\end{center}\n");
  fprintf(f, "\\end{table}\n");

  fclose(f);
  printf("Stored %s\n", fileName.Data());
}

void PlotterDiffXSec::PrintTexTabCorr(const UnfoldingIOHandler* unfoldingIOHandler, const ZMeasur* meas, const TString& fileName)
{
  assert(unfoldingIOHandler->Dim() == 2 || unfoldingIOHandler->Dim() == 3);
  //const bool flagLessLines = 1;
  utils::ttmd::CheckFile(fileName);
  FILE* f = fopen(fileName.Data(), "w");
  //const int dim = unfoldingIOHandler->Dim();
  const TH1D* hNom = unfoldingIOHandler->GetPlainHisto();
  TString strExtraBins;
  TString strXSec = utils::ttmd::TexXsecTitle(unfoldingIOHandler->Suffix, strExtraBins);
  TString strXSecLabel = utils::ttmd::TexXsecLabel(strXSec, unfoldingIOHandler->Var(0)->Bins[0].size() - 1);
  int flagLong = 0;
  if(hNom->GetNbinsX() == 24)
    flagLong = 1;
  else if(hNom->GetNbinsX() == 36)
    flagLong = 2;

  // header
  fprintf(f, "\\begin{table}[b]\n");
  fprintf(f, "\\begin{center}\n");
  fprintf(f, "\\caption{The correlation matrix of statistical uncertainties for the measured %s cross sections%s. The values are expressed as percentages. For bin indices see Table~\\ref{tab:%s_xsec}.}\n", strXSec.Data(), strExtraBins.Data(), strXSecLabel.Data());
  fprintf(f, "\\label{tab:%s_corr}\n", strXSecLabel.Data());
  //fprintf(f, "\\scriptsize\n");
  if(flagLong == 2)
  {
    fprintf(f, "\\tiny\n");
    fprintf(f, "\\tabcolsep0.3mm\n");
  }
  else if(flagLong == 1)
  {
    fprintf(f, "\\scriptsize\n");
    fprintf(f, "\\tabcolsep0.75mm\n");
  }
  else
  {
    fprintf(f, "\\footnotesize\n");
    fprintf(f, "\\tabcolsep2.0mm\n");
  }
  fprintf(f, "\\renewcommand*{\\arraystretch}{1.45}\n");
  //if(flagLong)
  //  fprintf(f, "\\tabcolsep0.17mm\n");
  //else
  //  fprintf(f, "\\tabcolsep0.5mm\n");
  fprintf(f, "\\begin{tabular}{l|");
  //for(int b = 0; b < hNom->GetNbinsX(); b++)
  for(int b = 1; b < hNom->GetNbinsX(); b++)
    fprintf(f, "r");
  fprintf(f, "}\n");
  fprintf(f, "bin ");
  //for(int b = 0; b < hNom->GetNbinsX(); b++)
  for(int b = 1; b < hNom->GetNbinsX(); b++)
    fprintf(f, "& $%d$ ", b + 1);
  fprintf(f, "\\\\ \n");
  fprintf(f, "\\hline \n");

  // main part
  //for(int b = 0; b < hNom->GetNbinsX(); b++)
  for(int b = 0; b < hNom->GetNbinsX() - 1; b++)
  {
    fprintf(f, "$%d$ ", b + 1);
    //for(int bb = 0; bb < hNom->GetNbinsX(); bb++)
    for(int bb = 1; bb < hNom->GetNbinsX(); bb++)
    {
      //if(bb < b)
      if(bb <= b)
        fprintf(f, "& ");
      else
      {
        double corr = meas->GetCovStatMatrix()->GetBinContent(b + 1, bb + 1) / meas->GetNom()->GetBinError(b + 1) / meas->GetNom()->GetBinError(bb + 1) * 100.0;
        //corr /= 100.0;
        //fprintf(f, "& $%+.2f$ ", corr);
        char buf[1024];
        sprintf(buf, "%+.0f", corr);
        TString bufStr = buf;
        if(bufStr == "+0" || bufStr == "-0")
          bufStr = "0";
        fprintf(f, "& $%s$ ", bufStr.Data());
        //fprintf(f, "& $%.0f$ ", corr);
      }
    }
    fprintf(f, "\\\\ \n");
  }

  // footer
  fprintf(f, "\\end{tabular}\n");
  fprintf(f, "\\end{center}\n");
  fprintf(f, "\\end{table}\n");

  fclose(f);
  printf("Stored %s\n", fileName.Data());
}

void PlotterDiffXSec::PrintTexTabSystHeader(FILE* f, const int flagLong, const int nBins, const TString& strXSecLabel, const TString& strXSecLabel0, const TString& strXSec, const TString& strExtraBins)
{
  if(flagLong == 2)
    fprintf(f, "\\begin{sidewaystable}[b]\n");
  else
    fprintf(f, "\\begin{table}[b]\n");
  fprintf(f, "\\begin{center}\n");
  if(strXSecLabel == strXSecLabel0)
    fprintf(f, "\\caption{Sources and values of the relative systematic uncertainties in percent of the measured %s cross sections%s. For bin indices see Table~\\ref{tab:%s_xsec}.}\n", strXSec.Data(), strExtraBins.Data(), strXSecLabel.Data());
  else
    fprintf(f, "\\caption{Table~\\ref{tab:%s_syst} continued.}\n", strXSecLabel0.Data());
  fprintf(f, "\\label{tab:%s_syst}\n", strXSecLabel.Data());
  fprintf(f, "\\renewcommand*{\\arraystretch}{1.75}\n");
  fprintf(f, "\\tiny\n");
  if(flagLong == 2)
    fprintf(f, "\\tabcolsep0.75mm\n");
  else if(flagLong == 1)
    fprintf(f, "\\tabcolsep0.70mm\n");
  else
    fprintf(f, "\\tabcolsep2.0mm\n");
  fprintf(f, "\\begin{tabular}{l|\n");
  for(int b = 0; b < nBins; b++)
    fprintf(f, "r");
  fprintf(f, "}\n");
  fprintf(f, "source / bin ");
  for(int b = 0; b < nBins; b++)
    fprintf(f, "& $%d$ ", b + 1);
  fprintf(f, "\\\\ \n");
}

void PlotterDiffXSec::PrintTexTabSystFooter(FILE* f, const int flagLong)
{
  fprintf(f, "\\end{tabular}\n");
  fprintf(f, "\\end{center}\n");
  if(flagLong == 2)
    fprintf(f, "\\end{sidewaystable}\n");
  else
    fprintf(f, "\\end{table}\n");
  fclose(f);
}

void PlotterDiffXSec::PrintTexTabSyst(PlotterConfigurationHelper *configHelper__, const UnfoldingIOHandler* unfoldingIOHandler, const ZMeasur* meas,
                                      std::vector<TString>& systematics_list, TString fileName)
{
  assert(unfoldingIOHandler->Dim() == 2 || unfoldingIOHandler->Dim() == 3);
  //const bool flagLessLines = 1;
  utils::ttmd::CheckFile(fileName);
  FILE* f = fopen(fileName.Data(), "w");
  const TString fileName0 = fileName;
  //const int dim = unfoldingIOHandler->Dim();
  const TH1D* hNom = unfoldingIOHandler->GetPlainHisto();
  int flagLong = 0;
  if(hNom->GetNbinsX() == 24)
    flagLong = 1;
  else if(hNom->GetNbinsX() == 36)
    flagLong = 2;

  TString strExtraBins;
  TString strXSec = utils::ttmd::TexXsecTitle(unfoldingIOHandler->Suffix, strExtraBins);
  TString strXSecLabel = utils::ttmd::TexXsecLabel(strXSec, unfoldingIOHandler->Var(0)->Bins[0].size() - 1);
  TString strXSecLabel0 = strXSecLabel;

  // header
  PrintTexTabSystHeader(f, flagLong, hNom->GetNbinsX(), strXSecLabel, strXSecLabel0, strXSec, strExtraBins);

  // main part
  std::map<TString, TH1D*> mSystOrig = meas->GetEigenvectors();
  std::vector<std::pair<TString, TH1D*> > vSystOrigFinal;

  for(size_t i = 0; i < systematics_list.size(); i++)
  {
      for (auto var: configHelper__->GetSysVars(systematics_list[i])){
          if (var != "Nominal"){
              vSystOrigFinal.push_back(std::pair<TString, TH1D*>(var, mSystOrig[var]));
              mSystOrig.erase(var);
          }
      }
  }

  for(auto it = mSystOrig.begin(); it != mSystOrig.end(); it++)
    printf("%s\n", it->first.Data());
  assert(mSystOrig.size() == 0);
  //for(auto it = vSystOrigFinal.begin(); it != vSystOrigFinal.end(); it++)
  //  printf("%s\n", it->first.Data());
  TString varNameFinalLast;
  int nSources = 0;
  int nSourcesMax = 26;
  if(flagLong == 2)
    nSourcesMax = 16;

  for(size_t i = 0; i < vSystOrigFinal.size(); i++)
  {
    TString varName = configHelper__->GetShortVarNameForTex(vSystOrigFinal[i].first);
    //printf("varName = %s\n", varName.Data());
    varName.ReplaceAll(" v", "|");
    TString varNameFinal = TString(varName.Data(), varName.First('|'));
    // 20.11.2018 final adjustment of systematic names
    varNameFinal.ReplaceAll("JESFlavor", "JESFlavour");
    varNameFinal.ReplaceAll("br(", "BR(");
    if(varNameFinal != varNameFinalLast)
    {
      nSources++;
      if(nSources % nSourcesMax == 0)
      {
        // start new table
        PrintTexTabSystFooter(f, flagLong);
        printf("Stored %s\n", fileName.Data());
        fileName = fileName0;
        fileName.ReplaceAll(".tex", TString::Format("-%d.tex", nSources / nSourcesMax));
        utils::ttmd::CheckFile(fileName);
        f = fopen(fileName.Data(), "w");
        TString strXSecLabelNew = strXSecLabel + TString::Format("_%d", nSources / nSourcesMax);
        PrintTexTabSystHeader(f, flagLong, hNom->GetNbinsX(), strXSecLabelNew, strXSecLabel, strXSec, strExtraBins);
      }
      fprintf(f, "\\hline \n");
      int nVar = 2;
      if(varName.BeginsWith("BFRAG"))
        nVar = 4;
      if(varName.BeginsWith("CR"))
        nVar = 3;
      fprintf(f, "\\multirow{%d}{*}{%s}", nVar, varNameFinal.Data());
      varNameFinalLast = varNameFinal;
    }
    for(int bb = 0; bb < hNom->GetNbinsX(); bb++)
    {
      double var = vSystOrigFinal[i].second->GetBinContent(bb + 1) / meas->GetNom()->GetBinContent(bb + 1) * 100.0;
      char buf[1024];
      sprintf(buf, "%+.1f", var);
      TString bufStr = buf;
      if(bufStr == "+0.0" || bufStr == "-0.0")
        bufStr = "0.0";
      fprintf(f, "& $%s$ ", bufStr.Data());
    }
    fprintf(f, "\\\\ \n");
  }


  printf("nSources: %d\n", nSources);

  // footer
  PrintTexTabSystFooter(f, flagLong);
  printf("Stored %s\n", fileName.Data());
}

void PlotterDiffXSec::PrintDataCorMatrix(const TH2D* hCor, const TString& fileName)
{
  // correlation matrix
  utils::ttmd::CheckFile(fileName);
  FILE* f = fopen(fileName.Data(), "w");
  fprintf(f, "%8s", "hyperbin");
  for(int h = 0; h < hCor->GetNbinsX(); h++)
    fprintf(f, "%8d", h);
  fprintf(f, "\n");
  for(int hh = 0; hh < hCor->GetNbinsX(); hh++)
  {
    fprintf(f, "%8d", hh);
    for(int h = 0; h < hCor->GetNbinsX(); h++)
    {
      if(h >= hh)
        fprintf(f, "%8.*f", NDigFloatCorr, hCor->GetBinContent(h + 1, hh + 1) * 100.0);
      else
        fprintf(f, "%8s", "");
    }
    fprintf(f, "\n");
  }
  fclose(f);
}

void PlotterDiffXSec::PrintCP(const std::vector<TH1D*> vh, const TString& fileName)
{
  assert(vh.size());
  utils::ttmd::CheckFile(fileName);
  FILE* f = fopen(fileName.Data(), "w");
  PrintHistBins(vh[0], f);
  fprintf(f, "%8s%12s%11s", "hyperbin", "data", "dataUnc");
  for(unsigned int i = 1; i < vh.size(); i++)
    fprintf(f, "%12s", vh[i]->GetTitle());
  fprintf(f, "\n");
  for(int b = 0; b < vh[0]->GetNbinsX(); b++)
  {
    fprintf(f, "%8d", b);
    for(unsigned int i = 0; i < vh.size(); i++)
    {
      fprintf(f, "%12.1f", vh[i]->GetBinContent(b + 1));
      if(i == 0)
        fprintf(f, "%11.1f", vh[i]->GetBinError(b + 1));
    }
    fprintf(f, "\n");
  }
  fclose(f);
}

void PlotterDiffXSec::PrintCPwithMCunc(const std::vector<TH1D*> vh, const std::pair<TH1D*, TH1D*> pairUD, const TString& fileName)
{
  std::vector<TH1D*> vhExt = vh;
  vhExt.push_back(pairUD.first);
  vhExt.back()->SetTitle("MCuncU");
  vhExt.push_back(pairUD.second);
  vhExt.back()->SetTitle("MCuncD");
  PrintCP(vhExt, fileName);
}

std::vector<TH1D*> PlotterDiffXSec::ReadCP(const TString& fileName)
{
  utils::ttmd::FileExists(fileName, 1);
  FILE* f = fopen(fileName.Data(), "r");
  const int maxBuf = 10000;
  std::vector<char> buf;
  buf.resize(maxBuf + 1);

  // read bins, create histograms
  std::vector<float> bins = ReadHistoBins(f);
  int nbins = bins.size() - 2;
  std::vector<TH1D*> vh;
  const int ncol = 5;
  for(int c = 0; c < ncol; c++)
    vh.push_back(new TH1D("", "", nbins, &bins[1]));

  // read historgam content
  if(!fgets(&buf[0], maxBuf, f)) utils::ttmd::ErFgets(); // skip line with column labels
  if(!fgets(&buf[0], maxBuf, f)) utils::ttmd::ErFgets();
  int counter = 0;
  while(!feof(f))
  {
    int hbin = -1;
    int offset = 0;
    int readN = 0;
    float ev = -1.0;
    // read bin
    sscanf(&buf[0], "%d%n", &hbin, &readN);
    for(int c = 0; c < ncol; c++)
    {
      offset += readN;
      sscanf(&buf[0] + offset, "%f%n", &ev, &readN);
      vh[c]->SetBinContent(hbin + 1, ev);
      if(c == 0)
      {
        // read data uncertainty
        offset += readN;
        sscanf(&buf[0] + offset, "%f%n", &ev, &readN);
        vh[0]->SetBinError(hbin + 1, ev);
      }
    }
    counter++;
    if(!fgets(&buf[0], maxBuf, f)) break;
  }
  assert(counter == nbins);
  fclose(f);
  return vh;
}

void PlotterDiffXSec::PrintEigenvectors(const TH1D* hNom, const std::pair<TH1D*, TH1D*> hUD, const std::map<std::pair<TString, TString>, TH1D*>& mapVars, const TString& fileName)
{
  utils::ttmd::FileExists(fileName);
  FILE* f = fopen(fileName.Data(), "w");

  // header line
  fprintf(f, "%8s%12s%8s%8s%8s", "hypbin", "xsec", "stat", "totU", "totD");
  for(std::map<std::pair<TString, TString>, TH1D*>::const_iterator it = mapVars.begin(); it != mapVars.end(); it++)
    fprintf(f, "%8s", it->first.second.Data());
  fprintf(f, "\n");

  // main part
  for(int b = 0;b < hNom->GetNbinsX(); b++)
  {
    fprintf(f, "%8d", b);
    fprintf(f, "%12.*e", NDigFloatXSec, hNom->GetBinContent(b + 1));
    fprintf(f, "%8.*f", NDigFloatStat, hNom->GetBinError(b + 1) / hNom->GetBinContent(b + 1) * 100.0);
    fprintf(f, "%8.*f", NDigFloatSyst, hUD.first->GetBinContent(b + 1) / hNom->GetBinContent(b + 1) * 100.0);
    fprintf(f, "%8.*f", NDigFloatSyst, -1 * hUD.second->GetBinContent(b + 1) / hNom->GetBinContent(b + 1) * 100.0);
    for(std::map<std::pair<TString, TString>, TH1D*>::const_iterator it = mapVars.begin(); it != mapVars.end(); it++)
      fprintf(f, "%8.*f", NDigFloatSyst, it->second->GetBinContent(b + 1) / hNom->GetBinContent(b + 1) * 100);
    fprintf(f, "\n");
  }

  fclose(f);
}

void PlotterDiffXSec::PrintDataXsecFinal(const TString& fileNameIn, const TString& fileNameOut, const TString& fileNameOutSys, const std::vector<std::vector<std::vector<TH1D*> > >& vvvSysAll, const std::vector<TString>& vShortName, const std::vector<std::vector<TH1D*> >& vvSysTot, const int a)
{
  utils::ttmd::CheckFile(fileNameIn, 0);
  FILE* fin = fopen(fileNameIn.Data(), "r");
  utils::ttmd::CheckFile(fileNameOut);
  FILE* fout = fopen(fileNameOut.Data(), "w");
  fprintf(fout, "%8s%8s%8s%8s%8s%12s%8s%8s%8s%12s\n", "bin1", "bin2", "bin3", "bin4", "hbin", "xsec", "stat", "sysU", "sysD", "genMC");
  utils::ttmd::CheckFile(fileNameOutSys);
  FILE* foutSys = fopen(fileNameOutSys.Data(), "w");
  fprintf(foutSys, "%8s%8s%8s%8s%8s%12s%8s%8s%8s", "bin1", "bin2", "bin3", "bin4", "hbin", "xsec", "stat", "sysU", "sysD");
  for(unsigned int s = 0; s < vShortName.size(); s++)
    fprintf(foutSys, "%7s", vShortName[s].Data());
  fprintf(foutSys, "\n");
  const int nbuf = 1023;
  std::vector<char> buf;
  buf.reserve(nbuf + 1);
  if(!fgets(&buf[0], nbuf, fin)) utils::ttmd::ErFgets();
  if(!fgets(&buf[0], nbuf, fin)) utils::ttmd::ErFgets();
  //int hbin = 0;
  while(!feof(fin))
  {
    float bin1, bin2, bin3, bin4, cs, stat, gen;
    int hbin;
    sscanf(&buf[0], "%f%f%f%f%d%f%f%*f%*f%f", &bin1, &bin2, &bin3, &bin4, &hbin, &cs, &stat, &gen);
    //hbin++;
    //sscanf(&buf[0], "%f%f%f%f%*d%f%f%*f%*f%f", &bin1, &bin2, &bin3, &bin4, &cs, &stat, &gen);
    double statRel = stat / cs * 100.0;
    double sysU = vvSysTot[0][a]->GetBinContent(hbin) / cs * 100;
    double sysD = vvSysTot[1][a]->GetBinContent(hbin) / cs * 100;
    fprintf(fout, "%8.2f%8.2f%8.2f%8.2f%8d%12.4e%8.2f%8.2f%8.2f%12.4e\n", bin1, bin2, bin3, bin4, hbin, cs, statRel, sysU, sysD, gen);
    // full systematics
    fprintf(foutSys, "%8.2f%8.2f%8.2f%8.2f%8d%12.4e%8.2f%8.2f%8.2f", bin1, bin2, bin3, bin4, hbin, cs, statRel, sysU, sysD);
    for(unsigned int s = 0; s < vvvSysAll.size(); s++)
    {
      fprintf(foutSys, "%7.2f", vvvSysAll[s][0][a]->GetBinContent(hbin) / cs * 100);
      fprintf(foutSys, "%7.2f", vvvSysAll[s][1][a]->GetBinContent(hbin) / cs * 100);
    }
    fprintf(foutSys, "\n");
    if(!fgets(&buf[0], nbuf, fin)) break;
  }
  fclose(fin);
  fclose(fout);
  fclose(foutSys);
}

std::pair<std::pair<TH1D*, TH1D*>, TH2D*> PlotterDiffXSec::ReadDataXsecAndCov(const TString& fileNameBase)
{
  // vector
  TString fileName = fileNameBase + "-vec.txt";
  utils::ttmd::FileExists(fileName, 1);
  FILE* f = fopen(fileName.Data(), "r");
  const int maxBuf = 10000;
  std::vector<char> buf;
  buf.resize(maxBuf + 1);

  // read bins, create histograms
  std::vector<float> bins = ReadHistoBins(f);
  int nbins = bins.size() - 2;
  TH1D* hDat = new TH1D("", "", nbins, &bins[1]);
  TH1D* hGen = new TH1D("", "", nbins, &bins[1]);
  TH2D* hCor = new TH2D("", "", nbins, &bins[1], nbins, &bins[1]);

  // read x-sections
  float cs = -1.0;
  float unc = 0.0;
  float gen = -1.0;
  if(!fgets(&buf[0], maxBuf, f)) utils::ttmd::ErFgets(); // skip line with column labels
  if(!fgets(&buf[0], maxBuf, f)) utils::ttmd::ErFgets();
  int hbin = -1;
  int counter = 0;
  while(!feof(f))
  {
    sscanf(&buf[0], "%d%f%f%f", &hbin, &cs, &unc, &gen);
    hDat->SetBinContent(hbin + 1, cs);
    hDat->SetBinError(hbin + 1, unc * cs / 100);
    hGen->SetBinContent(hbin + 1, gen);
    counter++;
    if(!fgets(&buf[0], maxBuf, f)) break;
  }
  if(nbins != counter)
    throw std::runtime_error(TString::Format("Error in ReadDataXsec(): nbins = %d counter = %d reading %s\n", nbins, counter, fileName.Data()).Data());
  fclose(f);

  // correlation matrix
  fileName = fileNameBase + "-cor.txt";
  if(utils::ttmd::CheckFile(fileName))
    hCor = NULL;
  else
  {
    f = fopen(fileName.Data(), "r");
    int binRead = -1;
    float cor = -1.0;
    std::vector<int> vHyperbins(nbins, -1);
    if(!fgets(&buf[0], maxBuf, f)) utils::ttmd::ErFgets();
    int readN = 0;
    sscanf(&buf[0], "%*s%n", &readN);
    int offset = readN;
    for(int h = 0; h < nbins; h++)
    {
      sscanf(&buf[0] + offset, "%d %n", &vHyperbins[h], &readN);
      offset += readN;
    }
    if(!fgets(&buf[0], maxBuf, f)) utils::ttmd::ErFgets();
    counter = 0;
    while(!feof(f))
    {
      sscanf(&buf[0], "%d %n", &binRead, &readN);
      offset = readN;
      for(int h = 0; h < nbins; h++)
      {
        if(h >= binRead)
        {
          sscanf(&buf[0] + offset, "%f %n", &cor, &readN);
          hCor->SetBinContent(binRead + 1, vHyperbins[h] + 1, cor / 100.0);
          hCor->SetBinContent(vHyperbins[h] + 1, binRead + 1, cor / 100.0);
          offset += readN;
        }
      }
      counter++;
      if(!fgets(&buf[0], maxBuf, f)) break;
    }
    if(nbins != counter)
      throw std::runtime_error(TString::Format("Error in ReadDataXsec(): nbins = %d counter = %d reading %s\n", nbins, counter, fileName.Data()).Data());
    fclose(f);

    //hDat->Print("all");
    //hGen->Print("all");
    //hCor->Print("all");
    //throw;

    // chi2
    //int skip = fileNameBase.Contains("-abs") ? 0 : 1;
    //double chi2 = Chi2(hDat, utils::ttmd::MakeCovMatrixFromCorr(hCor, hDat), hGen, skip);
    //printf("ReadDataXsec %s chi2 = %.1f\n", fileNameBase.Data(), chi2);
  }

  // return results
  std::pair<TH1D*, TH1D*> first = std::pair<TH1D*, TH1D*>(hDat, hGen);
  std::pair<std::pair<TH1D*, TH1D*>, TH2D*> result = std::pair<std::pair<TH1D*, TH1D*>, TH2D*>(first, hCor);
  return result;
}

std::pair<std::pair<TH1D*, TH1D*>, std::pair<TH1D*, std::pair<TH1D*, TH1D*>>> PlotterDiffXSec::ReadDataXsecWithSystAndStats(const TString& fileNameBase, bool get_total_uncert_UD){
    // vector
    TString fileName = fileNameBase + "-vec.txt";
    utils::ttmd::FileExists(fileName, 1);
    FILE* f = fopen(fileName.Data(), "r");
    const int maxBuf = 10000;
    std::vector<char> buf;
    buf.resize(maxBuf + 1);

    // read bins, create histograms
    std::vector<float> bins = ReadHistoBins(f);
    int nbins = bins.size() - 2;
    TH1D* hDat = new TH1D("", "", nbins, &bins[1]);
    TH1D* hStat = new TH1D("", "", nbins, &bins[1]);
    TH1D* hSyst_U = new TH1D("", "", nbins, &bins[1]);
    TH1D* hSyst_D = new TH1D("", "", nbins, &bins[1]);
    TH1D* hGen = new TH1D("", "", nbins, &bins[1]);

    // read x-sections
    float cs = -1.0;
    float stat_unc = 0.0;
    float syst_unc_u, syst_unc_d  = 0.0;
    float gen = -1.0;
    if(!fgets(&buf[0], maxBuf, f)) utils::ttmd::ErFgets(); // skip line with column labels
    if(!fgets(&buf[0], maxBuf, f)) utils::ttmd::ErFgets();
    int hbin = -1;
    int counter = 0;
    while(!feof(f)){
      sscanf(&buf[0], "%d%f%f%f%f%f", &hbin, &cs, &stat_unc, &syst_unc_u, &syst_unc_d, &gen);
      //std::cout << hbin << " " << cs << " " << stat_unc << " " << syst_unc_u << " " << syst_unc_d << " " << gen << std::endl;
      if (get_total_uncert_UD){
          float tot_uncertainty_u = TMath::Sqrt(stat_unc*stat_unc + syst_unc_u*syst_unc_u);
          float tot_uncertainty_d = TMath::Sqrt(stat_unc*stat_unc + syst_unc_d*syst_unc_d);
          syst_unc_u = tot_uncertainty_u;
          syst_unc_d = tot_uncertainty_d;
      }
      hDat->SetBinContent(hbin + 1, cs);
      hDat->SetBinError(hbin + 1, stat_unc * cs / 100);
      hGen->SetBinContent(hbin + 1, gen);
      hStat->SetBinContent(hbin + 1, stat_unc * cs / 100);
      hSyst_U->SetBinContent(hbin + 1, syst_unc_u * cs / 100);
      hSyst_D->SetBinContent(hbin + 1, -1 * syst_unc_d * cs / 100);

      counter++;
      if(!fgets(&buf[0], maxBuf, f)) break;
    }
    if(nbins != counter)
      throw std::runtime_error(TString::Format("Error in ReadDataXsec(): nbins = %d counter = %d reading %s\n", nbins, counter, fileName.Data()).Data());
    fclose(f);

    std::pair<TH1D*, TH1D*> data_and_mc = std::pair<TH1D*, TH1D*>(hDat, hGen);
    std::pair<TH1D*, TH1D*> syst_uncertainties = std::pair<TH1D*, TH1D*>(hSyst_U, hSyst_D);
    std::pair<TH1D*, std::pair<TH1D*, TH1D*>> uncertainties = std::pair<TH1D*, std::pair<TH1D*, TH1D*>>(hStat, syst_uncertainties);
    std::pair<std::pair<TH1D*, TH1D*>, std::pair<TH1D*, std::pair<TH1D*, TH1D*>>> result;
    result.first = data_and_mc;
    result.second = uncertainties;

    return result;

}

std::pair<std::pair<TH1D*, TH1D*>, TH2D*> PlotterDiffXSec::ReadRecUnf(const UnfoldingIOHandler* unfoldingIOHandler, const TString& dirName)
{
  // detector data and generator
  TString fileName = dirName + "/data.txt";
  utils::ttmd::FileExists(fileName, 1);
  FILE* f = fopen(fileName.Data(), "r");
  TString fileNameGen = dirName + "/gen.txt";
  utils::ttmd::FileExists(fileNameGen, 1);
  FILE* fGen = fopen(fileNameGen.Data(), "r");
  const int maxBuf = 10000;
  std::vector<char> buf;
  buf.resize(maxBuf + 1);
  printf("reading %s\n", fileNameGen.Data());

  // read bins, create histograms
  //std::vector<double> bins = unfoldingIOHandler->BinsF();

  TH1D* hDet = unfoldingIOHandler->GetPlainHisto();
  TH1D* hGen = unfoldingIOHandler->GetPlainHisto();
  int nbins = hDet->GetNbinsX();
  //assert(bins.size() == nbins);
  TH2D* hResp = new TH2D("", "", nbins, 0.0, nbins, nbins, 0.0, nbins);

  // read x-sections
  float det = -1.0;
  float unc = -1.0;
  float gen = -1.0;
  //int hbin = -1;
  int counter = 0;
  while(true)
  {
    if(!fgets(&buf[0], maxBuf, f)) break;
    sscanf(&buf[0], "%f%f", &det, &unc);
    hDet->SetBinContent(counter + 1, det);
    hDet->SetBinError(counter + 1, unc);
    // generator
    fgets(&buf[0], maxBuf, fGen);
    if(feof(f))
      break;
    gen = -999.9;
    int ret = sscanf(&buf[0], "%f", &gen);
    if(ret != 1 || gen == -999.9)
      break;
    hGen->SetBinContent(counter + 1, gen);
    counter++;
  }
  if(nbins != counter)
    throw std::runtime_error(TString::Format("Error in ReadRecUnf(): nbins = %d counter = %d reading %s\n", nbins, counter, fileName.Data()).Data());
  fclose(f);
  fclose(fGen);

  // response matrix
  fileName = dirName + "/resp.txt";
  utils::ttmd::FileExists(fileName, 1);
  f = fopen(fileName.Data(), "r");
  //int binRead = -1;
  float resp = -1.0;
  //std::vector<int> vHyperbins(nbins, -1);
  int readN = 0;
  int offset = 0;
  counter = 0;
  while(!feof(f))
  {
    if(!fgets(&buf[0], maxBuf, f)) break;
    counter++;
    offset = 0;
    //binRead = 0;
    for(int h = 0; h < nbins; h++)
    {
      int bx = counter;
      int by = h + 1;
      //if(h >= binRead)
      {
        sscanf(&buf[0] + offset, "%f %n", &resp, &readN);
        //assert(resp == 0.0);
        hResp->SetBinContent(bx, by, resp);
        //hResp->SetBinContent(by, bx, resp);
        offset += readN;
        //by++;
        //binRead++;
      }
    }
  }
  //hResp->Print("all");
  //throw;
  if(nbins != counter)
    throw std::runtime_error(TString::Format("Error in ReadRecUnf(): nbins = %d counter = %d reading %s\n", nbins, counter, fileName.Data()).Data());
  fclose(f);

  std::pair<std::pair<TH1D*, TH1D*>, TH2D*> result = std::pair<std::pair<TH1D*, TH1D*>, TH2D*>(std::pair<TH1D*, TH1D*>(hDet, hGen), hResp);
  return result;
}

void PlotterDiffXSec::PrintHistBins(const TH1* h, FILE* f)
{
  fprintf(f, "nbins %d bins", h->GetNbinsX());
  for(int b = 0; b < h->GetNbinsX() + 2; b++)
    fprintf(f, " %e", h->GetBinLowEdge(b));
  fprintf(f, "\n");
}

std::vector<float> PlotterDiffXSec::ReadHistoBins(FILE* f)
{
  const int maxBuf = 10000;
  std::vector<char> buf;
  buf.resize(maxBuf + 1);
  std::vector<float> bins;
  int nbins = -1;
  if(!fgets(&buf[0], maxBuf, f)) utils::ttmd::ErFgets();
  int readN = 0;
  sscanf(&buf[0], "%*s%d%*s%n", &nbins, &readN);
  if(nbins == -1)
    return bins;
  int counter = 0;
  int offset = readN;
  bins.resize(nbins + 2);
  while(readN > 1)
  {
    if(sscanf(&buf[0] + offset, "%f%n", &bins[counter], &readN) != 1)
      break;
    offset += readN;
    counter++;
  }
  if(nbins + 2 != counter)
    bins.clear();
  return bins;
}

void PlotterDiffXSec::CalculatePSE(const UnfoldingIOHandler* unfoldingIOHandler, const TString& fileName, const bool flagOrigBinning1D/* = true*/)
{
  double markerSize = 1.5;
  int pixel = configHelper->gPlotPixelSmall;

  const TH2D* hResp = unfoldingIOHandler->HRspC();
  const TH1D* hGen = unfoldingIOHandler->HGenC();
  //hResp->Print("all");
  int n = hResp->GetNbinsX();
  //TH1D* hp = new TH1D("", "", n, 0.0, n);
  //TH1D* hs = new TH1D("", "", n, 0.0, n);
  //TH1D* he = new TH1D("", "", n, 0.0, n);
  const bool binOrig = flagOrigBinning1D && unfoldingIOHandler->Dim() == 1;
  TH1D* hp = (TH1D*) unfoldingIOHandler->GetTUnfoldBinningC()->CreateHistogram("", binOrig);
  TH1D* hs = (TH1D*) unfoldingIOHandler->GetTUnfoldBinningC()->CreateHistogram("", binOrig);
  TH1D* he = (TH1D*) unfoldingIOHandler->GetTUnfoldBinningC()->CreateHistogram("", binOrig);
  TH1D* he1 = (TH1D*) unfoldingIOHandler->GetTUnfoldBinningC()->CreateHistogram("", binOrig);
  //TH1D* he2 = (TH1D*) unfoldingIOHandler->GetTUnfoldBinningC()->CreateHistogram("", binOrig);
  for(int b = 1; b <= n; b++)
  {
    double GenInBinAndRecInBin = hResp->GetBinContent(b, b);
    double RecInBin = hResp->Integral(b, b, 1, n);
    double GenInBinAndRec = hResp->Integral(1, n, b, b);
    double GenInBinAll = hGen->GetBinContent(b);
    hp->SetBinContent(b, GenInBinAndRecInBin / RecInBin);
    hs->SetBinContent(b, GenInBinAndRecInBin / GenInBinAndRec);
    // 23.07.17 significant differences between efficiency in 1D abd 2D
    he->SetBinContent(b, GenInBinAndRec / GenInBinAll);
    he1->SetBinContent(b, RecInBin / GenInBinAll);
    //he2->SetBinContent(b, GenInBinAndRecInBin / GenInBinAll);
  }
  //hp->Print("all");
  //hs->Print("all");
  //he->Print("all");

  TCanvas* c;
  if (unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote)
    {
      c = new TCanvas("cp", "", 800, PlotPixelY);
      c->SetMargin(0.12, 0.02, 0.12, 0.05);
    }
  else
    {
      c = new TCanvas("", "", pixel, pixel);
    }
  c->cd();
  TLegend* leg = new TLegend(0.15, 0.67, 0.40, 0.85);
  leg->SetTextFont(62);
  TH2D* hr;
  if (unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote) {
    leg = new TLegend(0.15, 0.75, 0.40, 0.90);
    hr = new TH2D("", "", 1, he->GetBinLowEdge(1), he->GetBinLowEdge(n + 1), 1, 0.0, 1.5);
  }
  else
    hr = new TH2D("", "", 1, he->GetBinLowEdge(1), he->GetBinLowEdge(n + 1), 1, 0.0, 1.0);
  if(binOrig)
    hr->GetXaxis()->SetTitle(unfoldingIOHandler->Var(0)->GetXTitle());
  else
    hr->GetXaxis()->SetTitle("Bin");
  hr->Draw();
  hp->SetMarkerColor(4);
  hp->SetMarkerStyle(23);
  hp->SetMarkerSize(markerSize);
  leg->AddEntry(hp, "Purity", "p");
  utils::ttmd::DrawAsGraph(hp);
  hs->SetMarkerColor(2);
  hs->SetMarkerStyle(22);
  hs->SetMarkerSize(markerSize);
  leg->AddEntry(hs, "Stability", "p");
  utils::ttmd::DrawAsGraph(hs);
  he->SetMarkerColor(8);
  he->SetMarkerStyle(20);
  he->SetMarkerSize(markerSize);
  leg->AddEntry(he, "Efficiency", "p");
  utils::ttmd::DrawAsGraph(he);
  // 23.07.17
  he1->SetMarkerColor(kGreen + 4);
  he1->SetMarkerStyle(24);
  he1->SetMarkerSize(markerSize);
  leg->AddEntry(he1, "Efficiency (R/G)", "p");
  utils::ttmd::DrawAsGraph(he1);
  //
  leg->Draw();
  if (unfoldingIOHandler->Dim() == 1 && FlagAnalysisNote) {
    TLine *cutOffLine = new TLine(he->GetBinLowEdge(1),1.,he->GetBinLowEdge(n + 1),1.);
    cutOffLine->SetLineColor(kRed);
    cutOffLine->SetLineWidth(1);
    cutOffLine->Draw("SAME");
    DrawCMSTitle(2, unfoldingIOHandler->Dim(), "pse");
  }
  utils::ttmd::SaveCanvas(c, fileName);
}

TH2D* PlotterDiffXSec::DrawMigrationMatrix(const UnfoldingIOHandler* unfoldingIOHandler, const TString& fileName, const bool flagOrigBinning1D) const
{
  int pixel = configHelper->gPlotPixelSmall;

  const TH2D* hResp = unfoldingIOHandler->HRspC();
  //const TH1D* hGen = unfoldingIOHandler->HGenC();
  int n = hResp->GetNbinsX();

  TH2D* hMig = (TH2D*) hResp->Clone();
  const bool binOrig = flagOrigBinning1D && unfoldingIOHandler->Dim() == 1;
  TString units = binOrig ? unfoldingIOHandler->Var(0)->GetXTitle() : "Bin";
  hMig->GetXaxis()->SetTitle(units + ", rec. level");
  hMig->GetYaxis()->SetTitle(units + ", gen. level");

  for(int by = 1; by <= n; by++)
  {
    //double gen = hGen->GetBinContent(by);
    //double GenAndRec = hResp->Integral(1, n, by, by);
    double GenAndRec = hResp->Integral(0, n, by, by);
    for(int bx = 1; bx <= n; bx++)
    for(int bx = 1; bx <= n; bx++)
    {
      double val = hResp->GetBinContent(bx, by);
      hMig->SetBinContent(bx, by, val / GenAndRec);
    }
  }

  if(fileName != "")
  {
    TCanvas c("", "", pixel, pixel);
    c.SetMargin(0.10, 0.15, 0.10, 0.05);
    // decrease font if there are many bins
    double scaleFont = 1.0;
    if(hMig->GetNbinsX() > 40)
      scaleFont = 0.37;
    else if(hMig->GetNbinsX() > 20)
      scaleFont = 0.60;
    hMig->SetMarkerSize(hMig->GetMarkerSize() * scaleFont);
    if(fileName.Contains("njmttyttfine"))
      hMig->Draw("colz");
    else
      hMig->Draw("text colz");
    //hMig->Draw("colz");
    utils::ttmd::SaveCanvas(&c, fileName);
  }

  //hResp->Print("all");
  return hMig;
}

void PlotterDiffXSec::PrintMigrationMatrix(const UnfoldingIOHandler* unfoldingIOHandler, const TString& fileName){
  utils::ttmd::CheckFile(fileName);
  ofstream mig_output_file;
  mig_output_file.open(fileName);
  if (mig_output_file.is_open()){
      const TH2D* hResp = unfoldingIOHandler->HRspC();
      int n = hResp->GetNbinsX();
      double mig_val = -999.9;

      for(int by = 1; by <= n; by++){
        double GenAndRec = hResp->Integral(0, n, by, by);
        for(int bx = 1; bx <= n; bx++){
          double val = hResp->GetBinContent(bx, by);
          mig_val = val / GenAndRec;
          mig_output_file << mig_val << endl;
        }
      }
      mig_output_file.close();
  }
  else{
      std::cerr << "ERROR in PlotterDiffXSec::PrintMigrationMatrix -->> Can't create output file: " << fileName << std::endl;
      exit(1);
  }
}


std::vector<float> PlotterDiffXSec::ReadMigrationMatrix(const TString& fileName, const bool convert_to_percent){
  ifstream mig_input_file;
  mig_input_file.open(fileName);

  double k_factor = 1.0;
  if (convert_to_percent) k_factor = 100;

  std::vector<float> plain_mig_matrix = {};

  if (mig_input_file.is_open()){
      float bin_content = -999.99;
      while (!mig_input_file.eof()) {
         mig_input_file >> bin_content;
         plain_mig_matrix.push_back(bin_content*k_factor);
      }
      mig_input_file.close();

      //for (auto i_val: plain_mig_matrix) std::cout << i_val << " ";
      //std::cout<< std::endl;
  }
  else{
      std::cerr << "ERROR in PlotterDiffXSec::ReadMigrationMatrix -->> Can't read input file: " << fileName << std::endl;
      exit(1);
  }

  return plain_mig_matrix;

}

void PlotterDiffXSec::PlotMigrationMatrixFromPlainVector(std::vector<float> plain_matrix, const UnfoldingIOHandler* unfoldingIOHandler, const TString& fileName, std::vector<int> bin_info, const bool flagOrigBinning1D){
    int pixel = configHelper->gPlotPixelSmall;

    bool plotting_info = (bin_info.size()>0);

    const TH2D* hResp = unfoldingIOHandler->HRspC();
    //const TH1D* hGen = unfoldingIOHandler->HGenC();
    int n = hResp->GetNbinsX();

    TH2D* hMig_info = (TH2D*) hResp->Clone();
    if (plotting_info) {
        assert(bin_info.size() == plain_matrix.size());
    }

    TH2D* hMig = (TH2D*) hResp->Clone();
    const bool binOrig = flagOrigBinning1D && unfoldingIOHandler->Dim() == 1;
    TString units = binOrig ? unfoldingIOHandler->Var(0)->GetXTitle() : "Bin";
    hMig->GetXaxis()->SetTitle(units + ", rec. level");
    hMig->GetYaxis()->SetTitle(units + ", gen. level");

    unsigned int i_bin_plain = 0;

    for(int by = 1; by <= n; by++)
    {
      for(int bx = 1; bx <= n; bx++)
      {
        if (plotting_info) hMig_info->SetBinContent(bx, by, bin_info.at(i_bin_plain));
        hMig->SetBinContent(bx, by, plain_matrix.at(i_bin_plain));
        i_bin_plain++;
      }
    }

    if(fileName != "")
    {
      TCanvas c("", "", pixel, pixel);
      c.SetMargin(0.10, 0.15, 0.10, 0.05);
      // decrease font if there are many bins
      double scaleFont = 1.0;
      if(hMig->GetNbinsX() > 40)
        scaleFont = 0.37;
      else if(hMig->GetNbinsX() > 20)
        scaleFont = 0.60;
      hMig->SetMarkerSize(hMig->GetMarkerSize() * scaleFont);
      if(fileName.Contains("njmttyttfine"))
        hMig->Draw("colz");
      else
        hMig->Draw("text colz");
      //hMig->Draw("colz");
      if (plotting_info) {
          hMig_info->SetBarOffset(0.2);
          hMig_info->Draw("TEXT SAME");
      }
      utils::ttmd::SaveCanvas(&c, fileName);
    }

}

void PlotterDiffXSec::PlotMultipleMigrationMatrixFromPlainVectors(std::map<TString, std::vector<float>> plain_matrices, const UnfoldingIOHandler* unfoldingIOHandler, const TString& fileName, std::vector<float> extra_bin_info, TString extra_bin_info_label, TString plot_title,const bool flagOrigBinning1D){
    int pixel = configHelper->gPlotPixelSmall;

    TString colbar_title = extra_bin_info_label;

    assert(plain_matrices.size()==3);
    bool include_extra_bin_info = (extra_bin_info.size()>0);

    const TH2D* hResp = unfoldingIOHandler->HRspC();
    //const TH1D* hGen = unfoldingIOHandler->HGenC();
    int n = hResp->GetNbinsX();

    TH2D* hMig_1 = (TH2D*) hResp->Clone();
    TH2D* hMig_2 = (TH2D*) hResp->Clone();
    TH2D* hMig_3 = (TH2D*) hResp->Clone();
    TH2D* hMig_info = (TH2D*) hResp->Clone();

    std::vector<TString> plain_matrix_names;
    std::vector<std::vector<float>> plain_matrices_data;

    for (auto it: plain_matrices){
        plain_matrices_data.push_back(it.second);
        plain_matrix_names.push_back(it.first);
    }

    const bool binOrig = flagOrigBinning1D && unfoldingIOHandler->Dim() == 1;
    TString units = binOrig ? unfoldingIOHandler->Var(0)->GetXTitle() : "Bin";

    unsigned int i_bin_plain = 0;

    for(int by = 1; by <= n; by++)
    {
      for(int bx = 1; bx <= n; bx++)
      {
        hMig_1->SetBinContent(bx, by, plain_matrices_data.at(0).at(i_bin_plain));
        hMig_2->SetBinContent(bx, by, plain_matrices_data.at(1).at(i_bin_plain));
        hMig_3->SetBinContent(bx, by, plain_matrices_data.at(2).at(i_bin_plain));
        if (include_extra_bin_info) hMig_info->SetBinContent(bx, by, extra_bin_info.at(i_bin_plain));
        i_bin_plain++;
      }
    }

    if(fileName != "")
    {
      TCanvas c("", "", pixel, pixel);
      c.SetMargin(0.12, 0.15, 0.12, 0.05);

      // decrease font if there are many bins
      double scaleFont = 1.0;
      if(hMig_1->GetNbinsX() > 40)
        scaleFont = 0.30;
      else if(hMig_1->GetNbinsX() > 20)
        scaleFont = 0.55;

      hMig_1->SetMarkerSize(hMig_1->GetMarkerSize() * scaleFont);
      hMig_2->SetMarkerSize(hMig_1->GetMarkerSize() * scaleFont);
      hMig_3->SetMarkerSize(hMig_1->GetMarkerSize() * scaleFont);

      if (include_extra_bin_info){
          gStyle->SetPalette(kBird);
          hMig_info->Draw("COLZ");
          hMig_info->GetZaxis()->SetTitle(extra_bin_info_label);
          hMig_info->GetZaxis()->SetTitleOffset(1.2);
          hMig_info->GetZaxis()->CenterTitle();
          hMig_1->Draw("TEXT SAME");
          hMig_info->GetXaxis()->SetTitle(units + ", rec. level");
          hMig_info->GetYaxis()->SetTitle(units + ", gen. level");
          hMig_info->SetTitle(plot_title);
      }
      else{
          hMig_1->Draw("TEXT");
          hMig_1->GetXaxis()->SetTitle(units + ", rec. level");
          hMig_1->GetYaxis()->SetTitle(units + ", gen. level");
          hMig_1->SetTitle(plot_title);
      }


      hMig_1->SetBarOffset(0.2);
      //hMig_1->SetTitle(plain_matrix_names.at(0));
      //hMig_1->SetMarkerColor(kRed);
      hMig_2->Draw("TEXT SAME");
      //hMig_2->SetMarkerColor(kBlue);
      //hMig_2->SetTitle(plain_matrix_names.at(1));
      hMig_3->SetBarOffset(-0.2);
      //hMig_3->SetMarkerColor(kSpring);
      //hMig_3->SetTitle(plain_matrix_names.at(2));
      hMig_3->Draw("TEXT SAME");

      auto legend = new TLegend(0.008, 0.12, 0.11, 0.008);
      legend->AddEntry(hMig_1, plain_matrix_names.at(0),"");
      legend->AddEntry(hMig_2, plain_matrix_names.at(1),"");
      legend->AddEntry(hMig_3, plain_matrix_names.at(2),"");
      legend->SetFillColor(kBlue-10);
      legend->SetTextSize(0.03);
      legend->Draw();
      utils::ttmd::SaveCanvas(&c, fileName);
    }

}
