#ifndef declareVars3D_h
#define declareVars3D_h

#include "ttmd_defineVars.h"
#include "UnfoldingIOHandler.h"
#include "utils.h"
#include "ttmd_ZMadGraph.h"
#include <cassert>

// ******************************************************************************************
// ************************* M(ttbar),y(ttbar) **********************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecNjMttYtt(PlotterConfigurationHelper *configHelper, const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const int njNBins = 3, const int mttNBins = 4, const int yttNBins = 4, const int mode = 0)
{
  ZVarNj* varNj = new ZVarNj(def, iso, pt);
  varNj->NDigits = 1;
  if(njNBins == 2)
    varNj->Bins.push_back({-0.5, 0.5, 8.5});
  else if(njNBins == 3)
    varNj->Bins.push_back({-0.5, 0.5, 1.5, 8.5});
  else if(njNBins == 4)
    varNj->Bins.push_back({-0.5, 0.5, 1.5, 2.5, 8.5});
  else
    throw std::logic_error(TString::Format("Error in DeclareXSecNjMttYtt(): njNBins = %d not supported", njNBins));

  if(njNBins == 4) varNj->Bins.push_back({-0.5, 0.5, 1.5, 2.5, 8.5});
  else varNj->Bins.push_back({-0.5, 0.5, 1.5, 8.5});

  ZVar* varMtt = new ZVarMtt(mode);
  if(mttNBins == 4)
    varMtt->Bins.push_back({300, 400, 500, 650, 1500});
  else if(mttNBins == 3)
    varMtt->Bins.push_back({300, 400, 500, 1500});
  else if(mttNBins == 8)
    varMtt->Bins.push_back({300, 380, 400, 425, 500, 550, 650, 800, 1500});
  else
    throw std::logic_error(TString::Format("Error in DeclareXSecNjMttYtt(): mttNBins = %d not supported", mttNBins));
  if(mode == 0)
    // varMtt->Bins.push_back({300,373,386,400,420,442,467,500,540,592,650,740,1500});
    varMtt->Bins.push_back({300,373,410,430,460,500,550,620,750,1500});
  else if(mode == 6 || mode == 9)
    varMtt->Bins.push_back({0,340,375,400,430,460,495,540,580,630,700,800,2500});
  else if(mode == 2)
    varMtt->Bins.push_back({0,310,340,370,430,460,500,550,620,750,2500});
  else
    assert(0);

  ZVarYtt* varYtt = new ZVarYtt(mode);
  varYtt->NDiv = 303;
  if(yttNBins == 4)
    varYtt->Bins.push_back({0.0, 0.35, 0.75, 1.15, 2.5});
  else if(yttNBins == 8)
    varYtt->Bins.push_back({0.0, 0.15, 0.35, 0.55, 0.75, 0.95, 1.15, 1.35, 2.5});
  else
    throw std::logic_error(TString::Format("Error in DeclareXSecNjMttYtt(): yttNBins = %d not supported", yttNBins));
  if(mode == 0)
    varYtt->Bins.push_back({0.0,0.08,0.16,0.25,0.35,0.47,0.61,0.75,0.92,1.15,1.41,2.5});
  else if(mode == 6 || mode == 9)
    varYtt->Bins.push_back({0.0,0.09,0.18,0.26,0.36,0.47,0.59,0.71,0.84,1.03,1.25,2.5});
  else if(mode == 2)
    varYtt->Bins.push_back({0.0,0.06,0.14,0.22,0.30,0.38,0.47,0.56,0.64,0.74,0.85,2.5});
  else
    assert(0);

  TString name = utils::ttmd::GetNameExtraJetVar("njmttytt", def, iso, pt);
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varNj, varMtt, varYtt}, configHelper, name + TString::Format("-b%d", njNBins));
  if(mttNBins == 3)
    unfoldingIOHandler->Suffix +=  TString::Format("-mtt%d", mttNBins);
  if(mttNBins == 8 && yttNBins == 8)
    unfoldingIOHandler->Suffix +=  TString::Format("-mtt%dytt%d", mttNBins, yttNBins);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->FlagNeedNP = false; // FIXME : for readout of npcorr; if no npcorr available - switch off.
  unfoldingIOHandler->NPertOrders = njNBins;

  if((njNBins == 2 || njNBins == 3) && (mttNBins == 3 || mttNBins == 4) && yttNBins == 4)
  {
    int offset = ZMadGraph::GetOffsetHistoPtjThreshold(configHelper, pt);
    // these histograms are provided only fot pt(jet) > 30 GeV, i.e. offset = 1
    if(offset != 1)
      throw std::runtime_error(TString::Format("Error: pt(jet) = %f not supported", pt).Data());

    /*// tt0j
    for(int i = 1; i <= 4; i++)
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({30 + i}, 0));
    for(int i = 1; i <= 4; i++)
      unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({30 + i}, 1));

    // tt1j
    for(int i = 1; i <= 4; i++)
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({30 + i}, 1));
    if(njNBins == 3)
      for(int i = 1; i <= 4; i++)
        unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({30 + i}, 2));
    else
      for(int i = 1; i <= 4; i++)
        unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));

    // tt2j (if 3 nj bins)
    if(njNBins == 3)
    {
      for(int i = 1; i <= 4; i++)
        unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({30 + i}, 2));
      for(int i = 1; i <= 4; i++)
        unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));
    }*/

    // tt0j
    unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({30 + 1}, 0));
    unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({30 + 2}, 0));
    if(mttNBins == 4)
    {
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({30 + 3}, 0));
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({30 + 4}, 0));
    }
    else if(mttNBins == 3)
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({30 + 3, 30 + 4}, 0));
    //
    unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({30 + 1}, 1));
    unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({30 + 2}, 1));
    if(mttNBins == 4)
    {
      unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({30 + 3}, 1));
      unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({30 + 4}, 1));
    }
    else if(mttNBins == 3)
      unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({30 + 3, 30 + 4}, 1));

    // tt1j
    unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({30 + 1}, 1));
    unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({30 + 2}, 1));
    if(mttNBins == 4)
    {
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({30 + 3}, 1));
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({30 + 4}, 1));
    }
    else if(mttNBins == 3)
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({30 + 3, 30 + 4}, 1));
    //
    if(njNBins == 3)
    {
      unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({30 + 1}, 2));
      unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({30 + 2}, 2));
      if(mttNBins == 4)
      {
        unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({30 + 3}, 2));
        unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({30 + 4}, 2));
      }
      else if(mttNBins == 3)
        unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({30 + 3, 30 + 4}, 2));
    }
    else
    {
      unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));
      unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));
      if(mttNBins == 4)
      {
        unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));
        unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));
      }
      else if(mttNBins == 3)
        unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));
    }

    // tt2j (if 3 nj bins)
    if(njNBins == 3)
    {
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({30 + 1}, 2));
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({30 + 2}, 2));
      if(mttNBins == 4)
      {
        unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({30 + 3}, 2));
        unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({30 + 4}, 2));
      }
      else if(mttNBins == 3)
        unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({30 + 3, 30 + 4}, 2));
      //
      unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));
      unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));
      if(mttNBins == 4)
      {
        unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));
        unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));
      }
      else if(mttNBins == 3)
        unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));
    }
  }
  unfoldingIOHandler->Init();

  return unfoldingIOHandler;
}

// fine binning detector level (1st and 2nd mtt bins combined)
UnfoldingIOHandler* DeclareXSecNjMttYttFineDet(PlotterConfigurationHelper *configHelper, const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const int njNBins = 3, const int nBinsTotal = 242)
{
  ZVarNj* varNj = new ZVarNj(def, iso, pt);
  varNj->NDigits = 1;
  if(njNBins == 2)
    varNj->Bins.push_back({-0.5, 0.5, 8.5});
  else if(njNBins == 3)
    varNj->Bins.push_back({-0.5, 0.5, 1.5, 8.5});
  else
    throw std::logic_error(TString::Format("Error in DeclareXSecNjMttYtt(): njNBins = %d not supported", njNBins));
  varNj->Bins.push_back(varNj->Bins.back());

  // mode 9 is default for this cross section
  ZVar* varMtt = new ZVarMtt(9);
  if(nBinsTotal == 242)
    varMtt->Bins.push_back({0,375,400,430,460,495,540,580,630,700,800,2500});
  else if(nBinsTotal == 264)
    varMtt->Bins.push_back({0,340,375,400,430,460,495,540,580,630,700,800,2500});
  else
    throw std::logic_error(TString::Format("Error in DeclareXSecNjMttYttFineDet(): nBinsTotal = %d not supported", nBinsTotal));
  if(nBinsTotal == 242)
    varMtt->Bins.push_back({0,375,400,430,460,495,540,580,630,700,800,2500});
  else if(nBinsTotal == 264)
    varMtt->Bins.push_back({0,340,375,400,430,460,495,540,580,630,700,800,2500});

  ZVarYtt* varYtt = new ZVarYtt(9);
  varYtt->NDiv = 303;
  varYtt->Bins.push_back({0.0,0.09,0.18,0.26,0.36,0.47,0.59,0.71,0.84,1.03,1.25,2.5});
  varYtt->Bins.push_back({0.0,0.09,0.18,0.26,0.36,0.47,0.59,0.71,0.84,1.03,1.25,2.5});

  TString name = utils::ttmd::GetNameExtraJetVar(TString::Format("njmttyttfine%d", nBinsTotal), def, iso, pt);
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varNj, varMtt, varYtt}, configHelper, name + TString::Format("-b%d", njNBins));
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  //unfoldingIOHandler->DoXSec = true;
  //unfoldingIOHandler->DoXSecFineDetOnly = true;
  unfoldingIOHandler->DoXSec = 2;
  /*unfoldingIOHandler->FlagNeedNP = true;
  unfoldingIOHandler->NPertOrders = njNBins;

  if(njNBins == 2 || njNBins == 3)
  {
    int offset = ZMadGraph::GetOffsetHistoPtjThreshold(pt);
    // these histograms are provided only fot pt(jet) > 30 GeV, i.e. offset = 1
    if(offset != 1)
      throw std::runtime_error(TString::Format("Error: pt(jet) = %f not supported", pt).Data());

    // tt0j
    if(nBinsTotal == 242)
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({411 + 1, 411 + 2}, 0));
    else if(nBinsTotal == 264)
    {
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({411 + 1}, 0));
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({411 + 2}, 0));
    }
    for(int yttbin = 3; yttbin <= 12; yttbin++)
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({411 + yttbin}, 0));

    // tt1j
    if(nBinsTotal == 242)
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({411 + 1, 411 + 2}, 1));
    else if(nBinsTotal == 264)
    {
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({411 + 1}, 1));
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({411 + 2}, 1));
    }
    for(int yttbin = 3; yttbin <= 12; yttbin++)
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({411 + yttbin}, 1));

    // tt2j (if 3 nj bins)
    if(njNBins == 3)
    {
      if(nBinsTotal == 242)
        unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({411 + 1, 411 + 2}, 2));
      else if(nBinsTotal == 264)
      {
        unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({411 + 1}, 2));
        unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({411 + 2}, 2));
      }
      for(int yttbin = 3; yttbin <= 12; yttbin++)
        unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({411 + yttbin}, 2));
      //
      if(nBinsTotal == 242)
        unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({}, 0));
      else if(nBinsTotal == 264)
      {
        unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({}, 0));
        unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({}, 0));
      }
      for(int yttbin = 3; yttbin <= 12; yttbin++)
        unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({}, 0));
    }
  }*/
  unfoldingIOHandler->Init();

  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecPtttMttYtt(PlotterConfigurationHelper *configHelper, const int mode = 0)
{
  ZVarPttt* varPttt = new ZVarPttt(mode);
  varPttt->Title = "p_{T}(t#bar{t})";
  varPttt->TitleTxt = "pT(ttbar)";
  varPttt->Units = "GeV";
  varPttt->NDigits = 0;
  varPttt->NDiv = 303;
  varPttt->Bins.push_back({0, 40, 120, 500});
  varPttt->Bins.push_back({0,20,40,50,65,90,120,160,500});

  ZVar* varMtt = new ZVarMtt(mode);
  //varMtt->Expression = "mtt";
  varMtt->Title = "M(t#bar{t})";
  varMtt->Units = "GeV";
  varMtt->NDigits = 0;
  varMtt->Bins.push_back({340, 400, 500, 650, 1500});
  //varMtt->Bins.push_back({340,366,373,380,386,393,400,410,420,431,442,454,467,482,500,519,540,563,592,620,650,685,740,850,1500});
  varMtt->Bins.push_back({340,380,400,428,460,500,560,650,740,1500});

  ZVarYtt* varYtt = new ZVarYtt(mode);
  //varYtt->Expression = "ytt";
  varYtt->Title = "y(t#bar{t})";
  varYtt->Units = "";
  varYtt->NDigits = 2;
  varYtt->NDiv = 303;
  varYtt->Bins.push_back({0.0, 0.35, 0.75, 1.15, 2.5});
  //varYtt->Bins.push_back({0.0,0.04,0.08,0.12,0.16,0.2,0.25,0.3,0.35,0.41,0.47,0.54,0.61,0.68,0.75,0.83,0.92,1.02,1.15,1.28,1.41,2.5});
  varYtt->Bins.push_back({0.0,0.15,0.25,0.35,0.47,0.61,0.75,0.92,1.15,1.30,2.5});

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varPttt, varMtt, varYtt}, configHelper, "ptttmttytt");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->NPertOrders = 2;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecPttt2BMttYtt(PlotterConfigurationHelper *configHelper, const double ptThreshold = 30.0, const int mode = 0)
{
  ZVarPttt* varPttt = new ZVarPttt(mode);
  varPttt->NDiv = 303;
  varPttt->Bins.push_back({0, ptThreshold, 1000});
  //varPttt->Bins.push_back({0,20,40,50,65,80,100,150,1000}); //8
  varPttt->Bins.push_back({0,33,62,100,2000}); //8

  ZVar* varMtt = new ZVarMtt(mode);
  varMtt->Bins.push_back({300, 400, 500, 650, 1500});
  if(mode == 0)
    varMtt->Bins.push_back({300,373,386,400,420,442,467,500,540,592,650,740,1500});
    //varMtt->Bins.push_back({340,380,400,440,500,560,650,740,1500}); //9
  else if(mode == 6 || mode == 9)
    varMtt->Bins.push_back({0,340,375,400,430,460,495,540,580,630,700,800,2500});
  else if(mode == 2)
    varMtt->Bins.push_back({0,310,340,370,430,460,500,550,620,750,2500});
  else
    assert(0);

  ZVarYtt* varYtt = new ZVarYtt(mode);
  varYtt->NDiv = 303;
  varYtt->Bins.push_back({0.0, 0.35, 0.75, 1.15, 2.5});
  if(mode == 0)
    //varYtt->Bins.push_back({0.0,0.15,0.25,0.35,0.47,0.61,0.75,0.92,1.15,1.35,2.5}); //10
    varYtt->Bins.push_back({0.0,0.08,0.16,0.25,0.35,0.47,0.61,0.75,0.92,1.15,1.41,2.5});
  else if(mode == 6 || mode == 9)
    varYtt->Bins.push_back({0.0,0.09,0.18,0.26,0.36,0.47,0.59,0.71,0.84,1.03,1.25,2.5});
  else if(mode == 2)
    varYtt->Bins.push_back({0.0,0.06,0.14,0.22,0.30,0.38,0.47,0.56,0.64,0.74,0.85,2.5});
  else
    assert(0);

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varPttt, varMtt, varYtt}, configHelper, TString::Format("pttt%.0fmttytt", ptThreshold));
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  // prediction
  unfoldingIOHandler->PutPredictionPttt(ptThreshold);

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecPttMttPttt(PlotterConfigurationHelper *configHelper)
{
  ZVarPtt* varPtt = new ZVarPtt;
  varPtt->Title = "p_{T}(t)";
  varPtt->TitleTxt = "pT(t)";
  varPtt->Expression = "ptt";
  varPtt->Units = "GeV";
  varPtt->NDigits = 0;
  varPtt->NDiv = 303;
  varPtt->Bins.push_back({0, 100, 600});
  varPtt->Bins.push_back({0,64,104,141,210,600});
  
  ZVar* varMtt = new ZVarMtt(0);
  //varMtt->Expression = "mtt";
  varMtt->Title = "M(t#bar{t})";
  varMtt->Units = "GeV";
  varMtt->NDigits = 0;
  varMtt->Bins.push_back({340, 450, 1500});
  //varMtt->Bins.push_back({340,366,373,380,386,393,400,410,420,431,442,454,467,482,500,519,540,563,592,620,650,685,740,850,1500});
  varMtt->Bins.push_back({340,400,460,560,740,1500});

  ZVarPttt* varPttt = new ZVarPttt;
  varPttt->Title = "p_{T}(t#bar{t})";
  varPttt->TitleTxt = "pT(ttbar)";
  varPttt->Expression = "pttt";
  varPttt->Units = "GeV";
  varPttt->NDigits = 0;
  varPttt->NDiv = 303;
  varPttt->Bins.push_back({0, 50, 500});
  varPttt->Bins.push_back({0,33,62,100,500});

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varPtt, varMtt, varPttt}, configHelper, "pttmttpttt");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->NPertOrders = 2;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecNjRhojYttj(PlotterConfigurationHelper *configHelper, const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const TString& leadStr = "2")
{
  ZVar* varNj = new ZVarNj(def, iso, pt);
  varNj->Bins.push_back({});
  int nbins = leadStr.Atoi() + 1;
  for(int b = 0; b < nbins; b++)
    varNj->Bins.back().push_back(-0.5 + b);
  varNj->Bins.back().push_back(8.5);
  varNj->Bins.push_back({-0.5, 0.5, 1.5, 8.5});

  ZVar* varRhoj = new ZVarRhoj(def, iso, pt, leadStr);
  varRhoj->Bins.push_back({0.0, 0.4, 0.5, 0.6, 1.0});
  varRhoj->Bins.push_back({0.0,0.35,0.42,0.46,0.5,0.54,0.58,0.62,0.67,0.73,1.0});

  ZVar* varYttj = new ZVarYttj(def, iso, pt, leadStr);
  varYttj->Bins.push_back({0.0, 0.3, 0.7, 1.2, 2.5});
  varYttj->Bins.push_back({0.0,0.15,0.3,0.5,0.7,0.9,1.15,1.45,2.5});
  varYttj->NDiv = 303;

  TString name = utils::ttmd::GetNameExtraJetVar("njrhojyttj", def, iso, pt, leadStr);
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varNj, varRhoj, varYttj}, configHelper, name + TString::Format("-b%d", nbins));
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->FlagNeedNP = true;
  if(nbins == 2)
    unfoldingIOHandler->NPertOrders = 2;
  else if(nbins == 3)
    unfoldingIOHandler->NPertOrders = 3;
  else
    unfoldingIOHandler->NPertOrders = 1;

  if(nbins == 2 || nbins == 3)
  {
    int offset = ZMadGraph::GetOffsetHistoPtjThreshold(configHelper, pt);
    int offset0 = 0;
    int hid0 = 340;
    int hid0Minus = 440;

    // tt0j
    for(int i = 1; i <= 4; i++)
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({hid0 + i * 10 + offset0}, 0));
    for(int i = 1; i <= 4; i++)
      unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({hid0Minus + i * 10 + offset}, 1));

    // tt1j
    for(int i = 1; i <= 4; i++)
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({hid0 + i * 10 + offset}, 1));
    if(nbins == 3)
      for(int i = 1; i <= 4; i++)
        unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({hid0Minus + i * 10 + offset}, 2));
    else
      for(int i = 1; i <= 4; i++)
        unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));

    // tt2j (if 3 nj bins)
    if(nbins == 3)
    {
      for(int i = 1; i <= 4; i++)
        unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({hid0 + i * 10 + offset}, 2));
      for(int i = 1; i <= 4; i++)
        unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));
    }
  }
  unfoldingIOHandler->Init();

  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecNjMttDeta(PlotterConfigurationHelper *configHelper, const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const int njNBins = 3, const int mttNBins = 4, const int detaNBins = 3, const int mode = 0){
    ZVarNj* varNj = new ZVarNj(def, iso, pt);
    varNj->NDigits = 1;
    if(njNBins == 2)
      varNj->Bins.push_back({-0.5, 0.5, 8.5});
    else if(njNBins == 3)
      varNj->Bins.push_back({-0.5, 0.5, 1.5, 8.5});
    else
      throw std::logic_error(TString::Format("Error in DeclareXSecNjMttDeta(): njNBins = %d not supported", njNBins));
    varNj->Bins.push_back({-0.5, 0.5, 1.5, 8.5});

    ZVar* varMtt = new ZVarMtt(mode);
    if(mttNBins == 4)
      varMtt->Bins.push_back({300, 400, 500, 650, 1500});
    else if(mttNBins == 3)
      varMtt->Bins.push_back({300, 400, 500, 1500});
    else if(mttNBins == 8)
      varMtt->Bins.push_back({300, 380, 400, 425, 500, 550, 650, 800, 1500});
    else
      throw std::logic_error(TString::Format("Error in DeclareXSecNjMttDeta(): mttNBins = %d not supported", mttNBins));
    if(mode == 0)
      // varMtt->Bins.push_back({300,373,386,400,420,442,467,500,540,592,650,740,1500});
      varMtt->Bins.push_back({300,373,410,430,460,500,550,620,750,1500});
    else if(mode == 6 || mode == 9)
      varMtt->Bins.push_back({0,340,375,400,430,460,495,540,580,630,700,800,2500});
    else if(mode == 2)
      varMtt->Bins.push_back({0,310,340,370,430,460,500,550,620,750,2500});
    else
      assert(0);

    ZVarDetatt* varDetatt = new ZVarDetatt;
    varDetatt->Title = "#Delta#eta(t,#bar{t})";
    varDetatt->TitleTxt = "delta_eta(t,tbar)";
    varDetatt->Expression = "detatt";
    varDetatt->Units = "";
    varDetatt->NDigits = 1;
    if (detaNBins==3) varDetatt->Bins.push_back({0.0, 0.4, 1.2, 6.0});
    else throw std::logic_error(TString::Format("Error in DeclareXSecNjMttDeta(): detaNBins = %d not supported", detaNBins));
    varDetatt->Bins.push_back({0.0,0.25,0.55,0.85,1.2,1.65,2.5,6.0});
    //varDetatt->Bins.push_back({0.0,0.12,0.25,0.4,0.55,0.70,0.85,1.0,1.2,1.4,1.65,2.0,2.5,6.0});

    TString name = utils::ttmd::GetNameExtraJetVar("njmttdeta", def, iso, pt);
    UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varNj, varMtt, varDetatt}, configHelper, name + TString::Format("-b%d", njNBins));

    unfoldingIOHandler->DoYields = true;
    unfoldingIOHandler->DoPSE = true;
    unfoldingIOHandler->DoXSec = true;
    unfoldingIOHandler->FlagNeedNP = false; // FIXME : for readout of npcorr; if no npcorr available - switch off.
    unfoldingIOHandler->NPertOrders = njNBins;

    unfoldingIOHandler->Init();

    return unfoldingIOHandler;

}

UnfoldingIOHandler* DeclareXSecNjMttPttt(PlotterConfigurationHelper *configHelper, const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const int njNBins = 3, const int mttNBins = 4, const int ptttNBins = 3, const int mode = 0){
    ZVarNj* varNj = new ZVarNj(def, iso, pt);
    varNj->NDigits = 1;
    if(njNBins == 2)
      varNj->Bins.push_back({-0.5, 0.5, 8.5});
    else if(njNBins == 3)
      varNj->Bins.push_back({-0.5, 0.5, 1.5, 8.5});
    else
      throw std::logic_error(TString::Format("Error in DeclareXSecNjMttPttt(): njNBins = %d not supported", njNBins));
    varNj->Bins.push_back({-0.5, 0.5, 1.5, 8.5});

    ZVar* varMtt = new ZVarMtt(mode);
    if(mttNBins == 4)
      varMtt->Bins.push_back({300, 400, 500, 650, 1500});
    else if(mttNBins == 3)
      varMtt->Bins.push_back({300, 400, 500, 1500});
    else if(mttNBins == 8)
      varMtt->Bins.push_back({300, 380, 400, 425, 500, 550, 650, 800, 1500});
    else
      throw std::logic_error(TString::Format("Error in DeclareXSecNjMttPttt(): mttNBins = %d not supported", mttNBins));
    if(mode == 0)
      // varMtt->Bins.push_back({300,373,386,400,420,442,467,500,540,592,650,740,1500});
      varMtt->Bins.push_back({300,373,410,430,460,500,550,620,750,1500});
    else if(mode == 6 || mode == 9)
      varMtt->Bins.push_back({0,340,375,400,430,460,495,540,580,630,700,800,2500});
    else if(mode == 2)
      varMtt->Bins.push_back({0,310,340,370,430,460,500,550,620,750,2500});
    else
      assert(0);

    ZVarPttt* varPttt = new ZVarPttt(mode);
    varPttt->Title = "p_{T}(t#bar{t})";
    varPttt->TitleTxt = "pT(ttbar)";
    varPttt->Units = "GeV";
    varPttt->NDigits = 0;
    varPttt->NDiv = 303;
    if (ptttNBins == 4) varPttt->Bins.push_back({0, 30, 75, 150, 500});
    else if (ptttNBins == 3) varPttt->Bins.push_back({0, 40, 120, 500});
    else throw std::logic_error(TString::Format("Error in DeclareXSecNjMttPttt(): ptttNBins = %d not supported", ptttNBins));
    varPttt->Bins.push_back({0,20,40,50,65,90,120,160,500});
    //varPttt->Bins.push_back({0,7,10,13,16,19,22,25,28,32,36,40,45,50,55,62,70,79,89,100,112,129,150,195,500});
    //varPttt->Bins.push_back({0,10,16,22,28,36,45,55,70,89,112,150,500});


    TString name = utils::ttmd::GetNameExtraJetVar("njmttpttt", def, iso, pt);
    UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varNj, varMtt, varPttt}, configHelper, name + TString::Format("-b%d", njNBins));

    unfoldingIOHandler->DoYields = true;
    unfoldingIOHandler->DoPSE = true;
    unfoldingIOHandler->DoXSec = true;
    unfoldingIOHandler->FlagNeedNP = false; // FIXME : for readout of npcorr; if no npcorr available - switch off.
    unfoldingIOHandler->NPertOrders = njNBins;

    unfoldingIOHandler->Init();

    return unfoldingIOHandler;

}

#endif
