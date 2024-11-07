#ifndef declareVars_h
#define declareVars_h

#include "ttmd_defineVars.h"
#include "ttmd_declareVars3D.h"
#include "UnfoldingIOHandler.h"
#include "utils.h"
#include "ttmd_ZMadGraph.h"
#include <cassert>

// ******************************************************************************************
// ************************* y(ttbar) *******************************************************
// ******************************************************************************************
// UnfoldingIOHandler* DeclareXSecYttSigned()
// {
//   ZVar* var = new ZVarYttSigned;
//   var->Bins.push_back({-2.5, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.5});
//   var->Bins.push_back({-2.5, -1.9, -1.7, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9, 2.5});
//   var->Title = "y(t#bar{t})";
//   var->Units = "";
//
//   UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, "ytts-xsec");
//   unfoldingIOHandler->DoYields = true;
//   unfoldingIOHandler->DoPSE = true;
//   unfoldingIOHandler->DoXSec = true;
//   unfoldingIOHandler->Init();
//   return unfoldingIOHandler;
// }

UnfoldingIOHandler* DeclareXSecYttXSec8(PlotterConfigurationHelper *configHelper,const int mode = 0, const int yttNBins = 10, const int nsubbins = 1)
{
  assert(yttNBins == 4 || yttNBins == 10);
  ZVar* var = new ZVarYtt(mode);
  if(yttNBins == 4)
    var->Bins.push_back({0.0, 0.35, 0.75, 1.15, 2.5});
  else
    var->Bins.push_back({0.0, 0.1, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25, 1.5, 1.75, 2.50});
  if(nsubbins != 1)
    var->Bins.back() = utils::ttmd::MakeSubbins(var->Bins.back(), nsubbins);
  if(mode == 0 || mode == 6 || mode == 9)
    var->Bins.push_back({0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.70,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.80,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.00,1.02,1.04,1.06,1.08,1.10,1.12,1.14,1.16,1.18,1.20,1.22,1.24,1.26,1.28,1.30,1.32,1.34,1.36,1.38,1.40,1.45,1.50,1.55,1.60,1.65,1.70,1.80,1.90,2.0,2.5});
  else if(mode == 2)
    var->Bins.push_back({0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.70,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.80,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.00,1.02,1.04,1.06,1.08,1.10,1.12,1.14,1.16,1.18,1.20,1.22,1.24,1.26,1.28,1.30,1.32,1.34,1.36,1.38,1.40,1.45,1.50,1.60,2.5});
  else assert(0);

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ytt-xsec8");
  if(yttNBins == 4)
    unfoldingIOHandler->Suffix +=  TString::Format("-ytt%d", yttNBins);
  if(nsubbins != 1)
    unfoldingIOHandler->Suffix += TString::Format("-nsubb%d", nsubbins);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  if(nsubbins == 1)
  {
    if(yttNBins == 4)
    {
      unfoldingIOHandler->FlagMGHistoIntegrate = 1;
      unfoldingIOHandler->vMGHistoOneBinOnly = {0, 1, 2, 3};
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({31,32,33,34}, 0) );
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({31,32,33,34}, 0) );
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({31,32,33,34}, 0) );
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({31,32,33,34}, 0) );
    }
    else if(yttNBins == 8)
      unfoldingIOHandler->VMGHisto = { std::pair<std::vector<int>, int>({4}, 0) };
  }
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

/*UnfoldingIOHandler* DeclareXSecYtt2j(double ptMin = 30.0)
{
  ZVar* var = new ZVarYtt2j(ptMin);
  var->Bins.push_back({0.0, 0.1, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25, 1.5, 1.75, 2.50});
  var->Bins.push_back({0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.70,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.80,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.00,1.02,1.04,1.06,1.08,1.10,1.12,1.14,1.16,1.18,1.20,1.22,1.24,1.26,1.28,1.30,1.32,1.34,1.36,1.38,1.40,1.45,1.50,1.55,1.60,1.65,1.70,1.80,1.90,2.0,2.5});

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, TString::Format("yttj%.0f", ptMin));
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}*/

// UnfoldingIOHandler* DeclareXSecYttCP50(const int mode = 0)
// {
//   ZVarYttSigned* var = new ZVarYttSigned(mode);
//   var->Bins.push_back({});
//   var->Bins.back().push_back(-2.5);
//   var->Bins.back().push_back(-2.15);
//   double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), -2.0, 40, 0.1);
//   assert(IsEqual(last, 2.0));
//   var->Bins.back().push_back(2.0);
//   var->Bins.back().push_back(2.15);
//   var->Bins.back().push_back(2.5);
//   var->Title = "y(t#bar{t})";
//   var->Units = "";
//
//   UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, "ytt-cp50");
//   unfoldingIOHandler->DoYields = true;
//   unfoldingIOHandler->DoPSE = true;
//   unfoldingIOHandler->Init();
//   return unfoldingIOHandler;
// }

UnfoldingIOHandler* DeclareXSecYttPSE25(PlotterConfigurationHelper *configHelper,const int mode = 0)
{
  ZVarYttSigned* var = new ZVarYttSigned(mode);
  var->Bins.push_back({});
  double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), -2.5, 25, 0.2);
  var->Bins.back().push_back(last);
  assert(IsEqual(var->Bins.back().back(), 2.5));
  var->Title = "y(t#bar{t})";
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ytt-pse25");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* M(ttbar) *******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecMttXSec8(PlotterConfigurationHelper *configHelper,const int mode = 0, const int mttNBins = 8, const int nsubbins = 1/*, const bool flagSameRecoBins = 0*/)
{
  assert(mttNBins == 3 || mttNBins == 4 || mttNBins == 8);
  ZVar* var = new ZVarMtt(mode);
  if(mttNBins == 4)
    var->Bins.push_back({300, 400, 500, 650, 1500});
  else if(mttNBins == 3)
    var->Bins.push_back({300, 400, 500, 1500});
  else
    var->Bins.push_back({300.0, 370.0, 395.0, 425.0, 460.0, 500.0, 570.0, 650.0, 750.0, 950.0, 1500.0});
  if(nsubbins != 1)
    var->Bins.back() = utils::ttmd::MakeSubbins(var->Bins.back(), nsubbins);
  if(mode == 0)
    var->Bins.push_back({300,360,365,370,375,380,385,390,395,400,405,410,415,420,425,430,435,440,445,450,455,460,465,470,475,480,485,490,495,500,505,510,515,520,525,530,535,540,545,550,555,560,565,570,575,580,590,600,610,620,630,640,650,660,670,680,700,720,740,760,780,800,850,900,950,1000,1100,1200,1300,1500});
  else if(mode == 2 || mode == 6 || mode == 9)
  {
    var->Bins.push_back({});
    //var->Bins.back().push_back(0.0);
    //var->Bins.back().push_back(200.0);
    //var->Bins.back().push_back(225.0);
    //var->Bins.back().push_back(240.0);
    //double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 250.0, 5, 10.0);
    //double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 300.0, 100, 5.0);
    //last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 800.0, 20, 10.0);
    //last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 1000.0, 10, 20.0);
    //last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 1200.0, 4, 50.0);
    //last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 1400.0, 2, 100.0);
    var->Bins.back().push_back(1625.0);
    var->Bins.back().push_back(1850.0);
    var->Bins.back().push_back(2500.0);
  }
  else
    assert(0);

  //if(flagSameRecoBins)
  //  var->Bins[1] = var->Bins[0];

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "mtt-xsec8");
  if(mttNBins == 3 || mttNBins == 4)
    unfoldingIOHandler->Suffix +=  TString::Format("-mtt%d", mttNBins);
  if(nsubbins != 1)
    unfoldingIOHandler->Suffix += TString::Format("-nsubb%d", nsubbins);
  //if(flagSameRecoBins)
  //  unfoldingIOHandler->Suffix += TString::Format("-eqgenrec");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  if(nsubbins == 1)
  {
    if(mttNBins == 3)
    {
      unfoldingIOHandler->FlagMGHistoIntegrate = 1;
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({31}, 0) );
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({32}, 0) );
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({33, 34}, 0) );
    }
    else if(mttNBins == 4)
    {
      unfoldingIOHandler->FlagMGHistoIntegrate = 1;
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({31}, 0) );
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({32}, 0) );
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({33}, 0) );
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({34}, 0) );
    }
    else if(mttNBins == 8)
      unfoldingIOHandler->VMGHisto = { std::pair<std::vector<int>, int>({5}, 0) };
  }
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

/*UnfoldingIOHandler* DeclareXSecMtt2j(double ptMin = 30.0)
{
  ZVar* var = new ZVarMtt2j(ptMin);
  var->Bins.push_back({340.0, 370.0, 395.0, 425.0, 460.0,500.0, 570.0, 650.0, 750.0, 950.0, 1500.0});
  var->Bins.push_back({340,360,365,370,375,380,385,390,395,400,405,410,415,420,425,430,435,440,445,450,455,460,465,470,475,480,485,490,495,500,505,510,515,520,525,530,535,540,545,550,555,560,565,570,575,580,590,600,610,620,630,640,650,660,670,680,700,720,740,760,780,800,850,900,950,1000,1100,1200,1300,1500});
  var->Title = TString::Format("M(t#bar{t}2j)_{p_{T}(j) > %.0f GeV}", ptMin);
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, TString::Format("mtt2j%.0f", ptMin));
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecMttaj()
{
  ZVar* var = new ZVarMttaj;
  var->Bins.push_back({340.0, 370.0, 395.0, 425.0, 460.0, 500.0, 570.0, 650.0, 750.0, 950.0, 1500.0});
  var->Bins.push_back({340,360,365,370,375,380,385,390,395,400,405,410,415,420,425,430,435,440,445,450,455,460,465,470,475,480,485,490,495,500,505,510,515,520,525,530,535,540,545,550,555,560,565,570,575,580,590,600,610,620,630,640,650,660,670,680,700,720,740,760,780,800,850,900,950,1000,1100,1200,1300,1500});

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, "mttaj");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}*/

// UnfoldingIOHandler* DeclareXSecMttCP50(const int mode = 0)
// {
//   ZVar* var = new ZVarMtt(mode);
//   var->Bins.push_back({});
//   var->Bins.back().push_back(340.0);
//   for(int b = 1; b < 15; b++)
//     var->Bins.back().push_back(var->Bins.back().back() + 8.0 * (1 + b / 12.0));
//   for(int b = 1; b < 15; b++)
//     var->Bins.back().push_back(var->Bins.back().back() + 18.0 * (1 + b / 22.0));
//   for(int b = 0; b < 15; b++)
//     var->Bins.back().push_back(var->Bins.back().back() + 30.0 * (1 + b / 5.5));
//   var->Bins.back().push_back(2000.0);
//   var->Title = "M(t#bar{t})";
//   var->Units = "GeV";
//   var->NDiv = 505;
//
//   UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, "mtt-cp50");
//   unfoldingIOHandler->DoYields = true;
//   unfoldingIOHandler->DoPSE = true;
//   unfoldingIOHandler->Init();
//   return unfoldingIOHandler;
// }

UnfoldingIOHandler* DeclareXSecMttPSE25(PlotterConfigurationHelper *configHelper,const int mode = 0, const int flagPaper = 0)
{
  ZVarMtt* var = new ZVarMtt(mode);
  var->Bins.push_back({});
  if(mode == 0 && flagPaper == 0)
  {
    double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 300.0, 25, 60.0);
    var->Bins.back().push_back(last);
    assert(IsEqual(var->Bins.back().back(), 1800.0));
  }
  else if(mode == 2 || mode == 6 || mode == 9 || flagPaper == 1)
  {
    //double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 0.0, 25, 100.0);
    double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 200.0, 18, 100.0);
    var->Bins.back().push_back(last);
    //assert(IsEqual(var->Bins.back().back(), 2500.0));
    assert(IsEqual(var->Bins.back().back(), 2000.0));
  }
  else
    assert(0);
  var->Title = "M(t#bar{t})";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "mtt-pse25");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

/*UnfoldingIOHandler* DeclareXSecMtt2jCP50(double ptMin = 30.0)
{
  ZVar* var = new ZVarMtt2j(ptMin);
  var->Bins.push_back({});
  var->Bins.back().push_back(340.0);
  for(int b = 1; b < 15; b++)
    var->Bins.back().push_back(var->Bins.back().back() + 8.0 * (1 + b / 12.0));
  for(int b = 1; b < 15; b++)
    var->Bins.back().push_back(var->Bins.back().back() + 18.0 * (1 + b / 22.0));
  for(int b = 0; b < 15; b++)
    var->Bins.back().push_back(var->Bins.back().back() + 30.0 * (1 + b / 5.5));
  var->Bins.back().push_back(2000.0);
  var->NDiv = 505;

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, TString::Format("mtt2j%.0f-cp50", ptMin));
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecMtt2jPSE25(double ptMin = 30.0)
{
  ZVarMtt2j* var = new ZVarMtt2j(ptMin);
  var->Bins.push_back({});
  double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 340.0, 25, 60.0);
  var->Bins.back().push_back(last);
  assert(IsEqual(var->Bins.back().back(), 1840.0));

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, TString::Format("mtt2j%.0f-pse25", ptMin));
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}*/

// ******************************************************************************************
// ************************* N(jets) ********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecNj(PlotterConfigurationHelper *configHelper,const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const int nbins = 2)
{
  assert(nbins >= 2 && nbins < 10);
  ZVar* var = new ZVarNj(def, iso, pt);
  var->Bins.push_back({});
  for(int b = 0; b < nbins; b++)
    var->Bins.back().push_back(-0.5 + b);
  var->Bins.back().push_back(8.5);
  if(pt == 100.0)
    var->Bins.push_back({-0.5, 0.5, 1.5, 2.5, 3.5, 8.5});
  else if(pt == 150.0)
    var->Bins.push_back({-0.5, 0.5, 1.5, 2.5, 8.5});
  else
    var->Bins.push_back({-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 8.5});

  TString name = var->Expression;
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, name + TString::Format("-b%d", nbins));
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
    unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({1}, 0));
    unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({300 + offset}, 1));
    unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({300 + offset}, 1));
    if(nbins == 3)
      unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({300 + offset}, 2));
    else
      unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));
    if(nbins == 3)
    {
      unfoldingIOHandler->VMGHisto.push_back(std::pair<std::vector<int>, int>({300 + offset}, 2));
      unfoldingIOHandler->VMGHistoMinus.push_back(std::pair<std::vector<int>, int>({}, 0));
    }
  }
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecNjCP(PlotterConfigurationHelper *configHelper,const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const int nLastBinBoundary = 4)
{
  ZVar* var = new ZVarNj(def, iso, pt);
  var->Bins.push_back({});
  for(int b = 0; b <= nLastBinBoundary + 1; b++)
    var->Bins.back().push_back(-0.5 + b);
  assert(var->Bins.back().back() < 8.0);
  var->Bins.back().push_back(8.5);

  TString name = var->Expression;
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, name + TString::Format("-cp%d", nLastBinBoundary));
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecNj_insideTTbarCP(PlotterConfigurationHelper *configHelper,const double iso = 0.4, const double pt = 30.0, const int nLastBinBoundary = 4)
{
  ZVar* var = new ZVarNj_inside(iso, pt);
  var->Bins.push_back({});
  for(int b = 0; b <= nLastBinBoundary + 1; b++)
    var->Bins.back().push_back(-0.5 + b);
  assert(var->Bins.back().back() < 8.0);
  var->Bins.back().push_back(8.5);

  TString name = var->Expression;
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, name + TString::Format("-cp%d", nLastBinBoundary));
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecNjPSE(PlotterConfigurationHelper *configHelper,const TString& def = "def", const double iso = 0.4, const double pt = 30.0)
{
  ZVar* var = new ZVarNj(def, iso, pt);
  var->Bins.push_back({});
  for(int b = 0; b <= 9; b++)
    var->Bins.back().push_back(-0.5 + b);
  assert(IsEqual(var->Bins.back().back(), 8.5));

  TString name = var->Expression;
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, name + "-pse");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* rho(ttj) *******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecRhoj(PlotterConfigurationHelper *configHelper,const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const TString& leadStr = "2",
                       const int binningScheme = 1, const bool flagAndNj1 = false, const int mode = 0)
{
  assert(flagAndNj1); // otherwise it cannot be compared to NLO
  ZVar* var = new ZVarRhoj(def, iso, pt, leadStr);

  if(flagAndNj1)
  {
    assert(leadStr == "1");
    var->ExpressionCut = utils::ttmd::GetNameExtraJetVar("nj", def, iso, pt);
    var->TitleExtension += " & N_{j} > 0";
  }

  // TOP-13-006: [0.0,0.2,0.3,0.45,0.6,1.0]
  if(binningScheme == 8)
    var->Bins.push_back({0.0, 0.2, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 1.0});
  else if(binningScheme == 4)
    var->Bins.push_back({0.0, 0.4, 0.5, 0.6, 1.0});
  else
    throw std::logic_error(TString::Format("Error in DeclareXSecRhoj(): binningScheme = %d not supported", binningScheme));

  var->Bins.push_back({0.0, 0.13, 0.16, 0.18, 0.20});
  for(int b = 0; b < 65; b++)
    var->Bins.back().push_back(0.20 + 0.01 * (b + 1));
  //printf("%f\n", var->Bins.back().back());
  assert(TMath::Abs(var->Bins.back().back() - 0.85) < 1e-6);
  var->Bins.back().push_back(0.88);
  var->Bins.back().push_back(0.91);
  //var->Bins.back().push_back(0.95);
  var->Bins.back().push_back(1.0);
  if(mode == 2)
  {
    var->Bins.back().push_back(1.1);
    var->Bins.back().push_back(1.2);
    var->Bins.back().push_back(1.3);
    //var->Bins.back().push_back(1.5);
    //var->Bins.back().push_back(1.75);
    var->Bins.back().push_back(2.5);
  }

  TString name = var->Expression;
  if(flagAndNj1)
    name += "-andNj1";
  if(binningScheme == 4)
    name += TString::Format("-b%d", binningScheme);
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, name);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->FlagNeedNP = true;
  unfoldingIOHandler->NPertOrders = 2;

  if(flagAndNj1)
  {
    //throw std::logic_error("Not implemented");
    int nptj = ZMadGraph::GetOffsetHistoPtjThreshold(configHelper, pt);
    if(binningScheme == 4)
    {
      unfoldingIOHandler->FlagMGHistoIntegrate = true;
      for(int i = 1; i <= 4; i++)
      {
        std::vector<int> vh = {340 + i * 10 + nptj};
        unfoldingIOHandler->VMGHisto.push_back(std::make_pair(vh, 1));
      }
    }
    else if(binningScheme == 8)
    {
      std::vector<int> vh = {330 + nptj};
      unfoldingIOHandler->VMGHisto.push_back(std::make_pair(vh, 1));
    }
  }
  else
  {
    throw std::logic_error("Not implemented");
  }
  unfoldingIOHandler->Init();

  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecRhojCP50(PlotterConfigurationHelper *configHelper,const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const TString& leadStr = "2", const bool flagAndNj1 = false)
{
  ZVar* var = new ZVarRhoj(def, iso, pt, leadStr);
  if(flagAndNj1)
  {
    var->ExpressionCut = utils::ttmd::GetNameExtraJetVar("nj", def, iso, pt);
    var->TitleExtension += " & N_{j} > 0";
  }
  var->Bins.push_back({});
  for(int b = 0; b < 51; b++)
    var->Bins.back().push_back(0.02 * b);
  assert(var->Bins.back().back() == 1.0);

  TString name = var->Expression;
  if(flagAndNj1)
    name += "-andNj1";
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, name + "-cp50");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecRhojPSE25(PlotterConfigurationHelper *configHelper,const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const TString& leadStr = "2", const bool flagAndNj1 = false)
{
  ZVar* var = new ZVarRhoj(def, iso, pt, leadStr);
  if(flagAndNj1)
  {
    var->ExpressionCut = utils::ttmd::GetNameExtraJetVar("nj", def, iso, pt);
    var->TitleExtension += " & N_{j} > 0";
  }
  var->Bins.push_back({});
  for(int b = 0; b < 26; b++)
    var->Bins.back().push_back(0.04 * b);
  assert(var->Bins.back().back() == 1.0);

  TString name = var->Expression;
  if(flagAndNj1)
    name += "-andNj1";
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, name + "-pse25");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* y(ttj) *********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecYttj(PlotterConfigurationHelper *configHelper,const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const TString& leadStr = "2")
{
  ZVar* var = new ZVarYttj(def, iso, pt, leadStr);
  var->Bins.push_back({0.0, 0.1, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25, 1.5, 1.75, 2.50});
  var->Bins.push_back({0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.70,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.80,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.00,1.02,1.04,1.06,1.08,1.10,1.12,1.14,1.16,1.18,1.20,1.22,1.24,1.26,1.28,1.30,1.32,1.34,1.36,1.38,1.40,1.45,1.50,1.55,1.60,1.65,1.70,1.80,1.90,2.0,2.5});

  TString name = var->Expression;
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, name);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecYttjCP50(PlotterConfigurationHelper *configHelper,const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const TString& leadStr = "2")
{
  ZVar* var = new ZVarYttjSigned(def, iso, pt, leadStr);
  var->Bins.push_back({});
  var->Bins.back().push_back(-2.5);
  var->Bins.back().push_back(-2.15);
  //double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), -2.0, 40, 0.1);
  assert(IsEqual(utils::ttmd::DVectorAppendNWidth(var->Bins.back(), -2.0, 40, 0.1), 2.0));
  var->Bins.back().push_back(2.0);
  var->Bins.back().push_back(2.15);
  var->Bins.back().push_back(2.5);

  TString name = var->Expression;
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, name + "-cp50");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecYttjPSE25(PlotterConfigurationHelper *configHelper,const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const TString& leadStr = "2")
{
  ZVar* var = new ZVarYttj(def, iso, pt, leadStr);
  var->Bins.push_back({});
  for(int b = 0; b < 26; b++)
    var->Bins.back().push_back(0.1 * b);
  assert(IsEqual(var->Bins.back().back(), 2.5));

  TString name = var->Expression;
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, name + "-pse25");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* y(t) ***********************************************************
// ******************************************************************************************

/*UnfoldingIOHandler* DeclareXSecYtLCP50()
{
  ZVarYtLSigned* var = new ZVarYtLSigned;
  var->Bins.push_back({});
  var->Bins.back().push_back(-2.5);
  var->Bins.back().push_back(-2.15);
  double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), -2.0, 40, 0.1);
  assert(IsEqual(last, 2.0));
  var->Bins.back().push_back(2.0);
  var->Bins.back().push_back(2.15);
  var->Bins.back().push_back(2.5);
  var->Title = "y(t_{L}})";
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, "ytl-cp50");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecYtLPSE25()
{
  ZVarYtLSigned* var = new ZVarYtLSigned;
  var->Bins.push_back({});
  double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), -2.5, 25, 0.2);
  var->Bins.back().push_back(last);
  assert(IsEqual(var->Bins.back().back(), 2.5));
  var->Title = "y(t_{L}})";
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, "ytl-pse25");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecYtSCP50()
{
  ZVarYtSSigned* var = new ZVarYtSSigned;
  var->Bins.push_back({});
  var->Bins.back().push_back(-2.5);
  var->Bins.back().push_back(-2.15);
  double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), -2.0, 40, 0.1);
  assert(IsEqual(last, 2.0));
  var->Bins.back().push_back(2.0);
  var->Bins.back().push_back(2.15);
  var->Bins.back().push_back(2.5);
  var->Title = "y(t_{S}})";
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, "yts-cp50");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecYtSPSE25()
{
  ZVarYtSSigned* var = new ZVarYtSSigned;
  var->Bins.push_back({});
  double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), -2.5, 25, 0.2);
  var->Bins.back().push_back(last);
  assert(IsEqual(var->Bins.back().back(), 2.5));
  var->Title = "y(t_{S}})";
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, "yts-pse25");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}*/

// UnfoldingIOHandler* DeclareXSecYtCP50()
// {
//   ZVarYtSigned* var = new ZVarYtSigned;
//   var->Bins.push_back({});
//   var->Bins.back().push_back(-2.5);
//   var->Bins.back().push_back(-2.15);
//   double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), -2.0, 40, 0.1);
//   assert(IsEqual(last, 2.0));
//   var->Bins.back().push_back(2.0);
//   var->Bins.back().push_back(2.15);
//   var->Bins.back().push_back(2.5);
//   var->Title = "y(t)";
//   var->Units = "";
//
//   UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, "yt-cp50");
//   unfoldingIOHandler->DoYields = true;
//   unfoldingIOHandler->DoPSE = true;
//   unfoldingIOHandler->Init();
//   return unfoldingIOHandler;
// }

UnfoldingIOHandler* DeclareXSecYtPSE25(PlotterConfigurationHelper *configHelper)
{
  ZVarYtSigned* var = new ZVarYtSigned;
  var->Bins.push_back({});
  double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), -2.5, 25, 0.2);
  var->Bins.back().push_back(last);
  assert(IsEqual(var->Bins.back().back(), 2.5));
  var->Title = "y(t)";
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "yt-pse25");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecYtXSec8(PlotterConfigurationHelper *configHelper)
{
  ZVarYt* var = new ZVarYt;
  var->Bins.push_back({0.0, 0.12, 0.3, 0.5, 0.70, 0.9, 1.15, 1.40, 1.65, 1.9, 2.50});
  var->Bins.push_back({0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.70,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.80,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.00,1.02,1.04,1.06,1.08,1.10,1.12,1.14,1.16,1.18,1.20,1.22,1.24,1.26,1.28,1.30,1.32,1.34,1.36,1.38,1.40,1.45,1.50,1.55,1.60,1.65,1.70,1.80,1.90,2.0,2.5});
  var->Title = "|y(t)|";
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "yt-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* pT(t) **********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecPttXSec8(PlotterConfigurationHelper *configHelper)
{
  ZVarPtt* var = new ZVarPtt;
  var->Bins.push_back({0, 35, 65, 95, 130, 165, 195, 240, 310, 400, 600});
  var->Bins.push_back({0,15,22,30,36,40,45,51,55,59,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,260,270,280,290,300,310,320,330,350,370,390,420,450,490,530,600});
  var->Title = "p_{T}(t)";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ptt-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->VMGHisto = { std::pair<std::vector<int>, int>({2}, 0) };
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// UnfoldingIOHandler* DeclareXSecPttCP50()
// {
//   ZVarPtt* var = new ZVarPtt;
//   var->Bins.push_back({});
//   var->Bins.back().push_back(0.0);
//   for(int b = 1; b < 15; b++)
//     var->Bins.back().push_back(var->Bins.back().back() + 8.0 * (1 + b / 30.0));
//   for(int b = 1; b < 15; b++)
//     var->Bins.back().push_back(var->Bins.back().back() + 12.0 * (1 + b / 30.0));
//   for(int b = 0; b < 15; b++)
//     var->Bins.back().push_back(var->Bins.back().back() + 20.0 * (1 + b / 19.0));
//   var->Bins.back().push_back(800.0);
//   var->Title = "p_{T}(t)";
//   var->Units = "GeV";
//
//   UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, "ptt-cp50");
//   unfoldingIOHandler->DoYields = true;
//   unfoldingIOHandler->DoPSE = true;
//   unfoldingIOHandler->Init();
//   return unfoldingIOHandler;
// }

UnfoldingIOHandler* DeclareXSecPttPSE25(PlotterConfigurationHelper *configHelper)
{
  ZVarPtt* var = new ZVarPtt;
  var->Bins.push_back({});
  double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 0.0, 25, 30.0);
  var->Bins.back().push_back(last);
  assert(IsEqual(var->Bins.back().back(), 750.0));
  var->Title = "p_{T}(t)";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ptt-pse25");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* pT(ttbar) ******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecPttt1D(PlotterConfigurationHelper *configHelper)
{
  ZVarPttt* var = new ZVarPttt;
  var->Title = "p_{T}(t#bar{t})";
  var->TitleTxt = "pT(ttbar)";
  var->Expression = "pttt";
  var->Units = "GeV";
  var->NDigits = 0;
  var->NDiv = 505;
  var->Bins.push_back({0, 15, 30, 50, 70, 100, 140, 200, 300, 450, 700});
  var->Bins.push_back({0,5,7,10,13,16,19,22,25,28,30,32,36,40,45,50,55,62,70,75,79,89,100,112,129,140,160,180,200,230,260,300,350,430,550,700});

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "pttt-1d");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecPtttPSE25(PlotterConfigurationHelper *configHelper)
{
  ZVarPttt* var = new ZVarPttt;
  var->Title = "p_{T}(t#bar{t})";
  var->TitleTxt = "pT(ttbar)";
  var->Expression = "pttt";
  var->Units = "GeV";
  var->Bins.push_back({});
  double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 0.0, 25, 30.0);
  var->Bins.back().push_back(last);
  assert(IsEqual(var->Bins.back().back(), 750));
  var->NDigits = 0;
  var->NDiv = 505;

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "pttt-pse25");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// UnfoldingIOHandler* DeclareXSecPtttCP50()
// {
//   ZVarPttt* var = new ZVarPttt;
//   var->Title = "p_{T}(t#bar{t})";
//   var->TitleTxt = "pT(ttbar)";
//   var->Expression = "pttt";
//   var->Units = "GeV";
//   var->Bins.push_back({});
//   var->Bins.back().push_back(0.0);
//   for(int b = 1; b < 15; b++)
//     var->Bins.back().push_back(var->Bins.back().back() + 4.0 * (1 + b / 15.0));
//   for(int b = 1; b < 15; b++)
//     var->Bins.back().push_back(var->Bins.back().back() + 8.0 * (1 + b / 15.0));
//   for(int b = 0; b < 16; b++)
//     var->Bins.back().push_back(var->Bins.back().back() + 17.0 * (1 + b / 12.0));
//   var->Bins.back().push_back(750.0);
//   var->NDigits = 0;
//   var->NDiv = 505;
//
//   UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, "pttt-cp50");
//   unfoldingIOHandler->DoYields = true;
//   unfoldingIOHandler->DoPSE = true;
//   unfoldingIOHandler->Init();
//   return unfoldingIOHandler;
// }


UnfoldingIOHandler* DeclareXSecPttt2B(PlotterConfigurationHelper *configHelper,const double ptThreshold = 30.0)
{
  // pttt with 2 bins
  std::vector<double> binsPttt2BRec = {0,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,
                                       56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,92,94,96,98,100,102,104,106,108,110,112,115,118,121,124,127,130,134,138,
                                       142,146,150,155,161,165,170,175,180,185,190,195,200,210,220,230,240,250,260,280,300,330,360,400,450,500};
  ZVarPttt* varPttt = new ZVarPttt;
  varPttt->Title = "p_{T}(t#bar{t})";
  varPttt->TitleTxt = "pT(ttbar)";
  varPttt->Expression = "pttt";
  varPttt->Units = "GeV";
  varPttt->NDigits = 0;
  varPttt->NDiv = 303;
  varPttt->Bins.push_back({0, ptThreshold, 500});
  varPttt->Bins.push_back(binsPttt2BRec);

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varPttt}, configHelper, TString::Format("pttt%.0f", ptThreshold));
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->NPertOrders = 2;

  // prediction
  unfoldingIOHandler->PutPredictionPtttIntegrate(ptThreshold);
  unfoldingIOHandler->FlagMGHistoIntegrate = true;
  unfoldingIOHandler->Init();

  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* delta_phi(t,tbar) **********************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecDphitt1D(PlotterConfigurationHelper *configHelper)
{
  ZVarDphitt* var = new ZVarDphitt;
  var->Title = "#Delta#phi(t,#bar{t})";
  var->TitleTxt = "delta_phi(t,tbar)";
  var->Expression = "dphitt";
  var->Units = "rad";
  var->NDigits = 2;
  var->NDiv = 505;
  var->Bins.push_back({0, 1.0, 1.7, 2.2, 2.6, 2.9, 3.05, 3.1415});
  var->Bins.push_back({0,0.25,0.5,0.75,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.65,2.70,2.75,2.80,2.85,2.9,2.95,3.0,3.05,3.08,3.10,3.1415});

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "dphitt-1d");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecDphittPSE25(PlotterConfigurationHelper *configHelper)
{
  ZVarDphitt* var = new ZVarDphitt;
  var->Title = "#Delta#phi(t,#bar{t})";
  var->TitleTxt = "delta_phi(t,tbar)";
  var->Expression = "dphitt";
  var->Units = "rad";
  var->NDigits = 2;
  var->NDiv = 505;
  var->Bins.push_back({});
  double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 0.0, 25, (TMath::Pi() + 1e-5) / 25.0);
  var->Bins.back().push_back(last);
  assert(IsEqual(var->Bins.back().back(), (TMath::Pi() + 1e-5)));

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "dphitt-pse25");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// UnfoldingIOHandler* DeclareXSecDphittCP50()
// {
//   ZVarDphitt* var = new ZVarDphitt;
//   var->Title = "#Delta#phi(t,#bar{t})";
//   var->TitleTxt = "delta_phi(t,tbar)";
//   var->Expression = "dphitt";
//   var->Units = "rad";
//   var->NDigits = 2;
//   var->NDiv = 505;
//   var->Bins.push_back({});
//   var->Bins.back().push_back(0.0);
//   for(int b = 1; b < 15; b++)
//     var->Bins.back().push_back(var->Bins.back().back() + 0.17 * (1 - b / 30.0));
//   for(int b = 1; b < 15; b++)
//     var->Bins.back().push_back(var->Bins.back().back() + 0.085 * (1 - b / 30.0));
//   for(int b = 0; b < 16; b++)
//     var->Bins.back().push_back(var->Bins.back().back() + 0.04 * (1 - b / 26.0));
//   var->Bins.back().push_back(TMath::Pi() + 1e-5);
//
//   UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, "dphitt-cp50");
//   unfoldingIOHandler->DoYields = true;
//   unfoldingIOHandler->DoPSE = true;
//   unfoldingIOHandler->Init();
//   return unfoldingIOHandler;
// }

// ******************************************************************************************
// ************************* delta_eta(t,tbar) **********************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecDetatt1D(PlotterConfigurationHelper *configHelper)
{
  ZVarDetatt* var = new ZVarDetatt;
  var->Title = "#Delta#eta(t,#bar{t})";
  var->TitleTxt = "delta_eta(t,tbar)";
  var->Expression = "detatt";
  var->Units = "";
  var->NDigits = 2;
  var->NDiv = 505;
  var->Bins.push_back({0.0, 0.2, 0.5, 0.8, 1.2, 1.7, 2.5, 3.5, 5.0, 7.0});
  var->Bins.push_back({0.0,0.03,0.06,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.5,3.0,3.5,4.0,4.5,5.0,6.0,7.0});

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "detatt-1d");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecDetattPSE25(PlotterConfigurationHelper *configHelper)
{
  ZVarDetatt* var = new ZVarDetatt;
  var->Title = "#Delta#eta(t,#bar{t})";
  var->TitleTxt = "delta_eta(t,tbar)";
  var->Expression = "detatt";
  var->Units = "";
  var->NDigits = 1;
  var->NDiv = 505;
  var->Bins.push_back({});
  double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 0.0, 25, 7.5 / 25);
  var->Bins.back().push_back(last);
  assert(IsEqual(var->Bins.back().back(), 7.5));

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "detatt-pse25");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}


// UnfoldingIOHandler* DeclareXSecDetattCP50()
// {
//   ZVarDetatt* var = new ZVarDetatt;
//   var->Title = "#Delta#eta(t,#bar{t})";
//   var->TitleTxt = "delta_eta(t,tbar)";
//   var->Expression = "detatt";
//   var->Units = "";
//   var->NDigits = 1;
//   var->NDiv = 505;
//   var->Bins.push_back({});
//   var->Bins.back().push_back(0.0);
//   for(int b = 1; b < 15; b++)
//     var->Bins.back().push_back(var->Bins.back().back() + 0.08 * (1 + b / 30.0));
//   for(int b = 1; b < 15; b++)
//     var->Bins.back().push_back(var->Bins.back().back() + 0.12 * (1 + b / 30.0));
//   for(int b = 0; b < 15; b++)
//     var->Bins.back().push_back(var->Bins.back().back() + 0.20 * (1 + b / 19.0));
//   var->Bins.back().push_back(8.0);
//
//   UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, "detatt-cp50");
//   unfoldingIOHandler->DoYields = true;
//   unfoldingIOHandler->DoPSE = true;
//   unfoldingIOHandler->Init();
//   return unfoldingIOHandler;
// }

// ******************************************************************************************
// ************************** total *********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecTotal(PlotterConfigurationHelper *configHelper, const int mode = 0)
{
  ZVar* var = new ZVarYtt(mode);
  var->Bins.push_back({0.0, 10.0});
  var->Bins.push_back({0.0, 10.0});
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "total");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* N jets *********************************************************
// ******************************************************************************************
/*
UnfoldingIOHandler* DeclareXSecNj3(const double ptMin = 30.0)
{
  ZVarNj* var = new ZVarNj(ptMin);
  var->Bins.push_back({-0.5, 0.5, 1.5, 8.5});
  //var->Bins.push_back({1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5});
  if(ptMin != 150.0)
    var->Bins.push_back({-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 8.5});
  else
    var->Bins.push_back({-0.5, 0.5, 1.5, 2.5, 3.5, 8.5});

  TString dir = "nj-3b";
  if(ptMin != 30.0)
    dir = TString::Format("nj%.0f-3b", ptMin);
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, dir);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecNj2(const double ptMin = 30.0)
{
  ZVarNj* var = new ZVarNj(ptMin);
  var->Bins.push_back({-0.5, 0.5, 8.5});
  if(ptMin != 150.0)
    var->Bins.push_back({-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 8.5});
  else
    var->Bins.push_back({-0.5, 0.5, 1.5, 2.5, 3.5, 8.5});

  TString dir = "nj-2b";
  if(ptMin != 30.0)
    dir = TString::Format("nj%.0f-2b", ptMin);
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, dir);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecNj6(const double ptMin = 30.0)
{
  ZVarNj* var = new ZVarNj(ptMin);
  var->Bins.push_back({-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 8.5});
  var->Bins.push_back({-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 8.5});

  TString dir = "nj-6b";
  if(ptMin != 30.0)
    dir = TString::Format("nj%.0f-6b", ptMin);
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, dir);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecNjCP(const double ptMin = 30.0)
{
  ZVarNj* var = new ZVarNj(ptMin);
  var->Bins.push_back({-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 8.5});

  TString dir = "nj-cp";
  if(ptMin != 30.0)
    dir = TString::Format("nj%.0f-cp", ptMin);
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, dir);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecNjCP2(const double ptMin = 30.0)
{
  ZVarNj* var = new ZVarNj(ptMin);
  var->Bins.push_back({-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5});

  TString dir = "nj-cp";
  if(ptMin != 30.0)
    dir = TString::Format("nj%.0f-cp", ptMin);
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, dir);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecNjPSE(const double ptMin = 30.0)
{
  ZVarNj* var = new ZVarNj(ptMin);
  var->Bins.push_back({-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5});

  TString dir = "nj-pse";
  if(ptMin != 30.0)
    dir = TString::Format("nj%.0f-pse", ptMin);
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, dir);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}
*/
// ******************************************************************************************
// ************************* "converters" ***************************************************
// ******************************************************************************************

// UnfoldingIOHandler* DeclareXSecCp25FromCP50(UnfoldingIOHandler* unfoldingIOHandler0)
// {
//   ZVar* var0 = unfoldingIOHandler0->Var(0);
//
//   // skip jet multiplicity
//   if(unfoldingIOHandler0->Expression == "njet")
//     return unfoldingIOHandler0;
//
//   std::vector<double> bins;
//   for(int b = 0; b < var0->BinsC().size(); b = b + 2)
//     bins.push_back(var0->BinsC()[b]);
//   auto var = var0;
//   var->Bins.clear();
//   var->Bins.push_back(bins);
//
//   var->Title = var0->Title;
//   var->Units = var0->Units;
//
//   UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, unfoldingIOHandler0->Suffix.ReplaceAll("cp50", "cp25"));
//   unfoldingIOHandler->DoYields = unfoldingIOHandler0->DoYields;
//   unfoldingIOHandler->DoPSE = unfoldingIOHandler0->DoPSE;
//   unfoldingIOHandler->Init();
//
//   delete unfoldingIOHandler0;
//   return unfoldingIOHandler;
// }

// ******************************************************************************************
// ************************* summary methods ************************************************
// ******************************************************************************************

// std::vector<UnfoldingIOHandler*> DeclareXSecAllCP50()
// {
//   std::vector<UnfoldingIOHandler*> vec;
//   // vec.push_back(DeclareXSecMttCP50());
//   // vec.push_back(DeclareXSecYttCP50());
//   // vec.push_back(DeclareXSecPttCP50());
//   // vec.push_back(DeclareXSecYtCP50());
//   //vec.push_back(DeclareXSecYtLCP50());
//   //vec.push_back(DeclareXSecYtSCP50());
//   // vec.push_back(DeclareXSecPtttCP50());
//   // vec.push_back(DeclareXSecDetattCP50());
//   // vec.push_back(DeclareXSecDphittCP50());
//   vec.push_back(DeclareXSecNjCP("def", 0.4, 30, 4));
//   vec.push_back(DeclareXSecNjCP("def", 0.4, 30, 7));
//   //vec.push_back(DeclareXSecNjCP2());
//   return vec;
// }

// std::vector<UnfoldingIOHandler*> DeclareXSecAllCP25()
// {
//   std::vector<UnfoldingIOHandler*> vec;
//   auto vunfoldingIOHandler0 = DeclareXSecAllCP50();
//   //for_each (vunfoldingIOHandler0.begin(), vunfoldingIOHandler0.end(), DeclareXSecCp25FromCP50);
//   for(auto unfoldingIOHandler0 : vunfoldingIOHandler0)
//     vec.push_back(DeclareXSecCp25FromCP50(unfoldingIOHandler0));
//   return vec;
// }

std::vector<UnfoldingIOHandler*> DeclareXSecAllPSE25(PlotterConfigurationHelper *configHelper)
{
  std::vector<UnfoldingIOHandler*> vec;
  vec.push_back(DeclareXSecMttPSE25(configHelper));
  vec.push_back(DeclareXSecYttPSE25(configHelper));
  vec.push_back(DeclareXSecPttPSE25(configHelper));
  vec.push_back(DeclareXSecYtPSE25(configHelper));
  //vec.push_back(DeclareXSecYtLPSE25());
  //vec.push_back(DeclareXSecYtSPSE25());
  vec.push_back(DeclareXSecPtttPSE25(configHelper));
  vec.push_back(DeclareXSecDetattPSE25(configHelper));
  vec.push_back(DeclareXSecDphittPSE25(configHelper));
  return vec;
}

// std::vector<UnfoldingIOHandler*> DeclareXSecAll81D()
// {
//   std::vector<UnfoldingIOHandler*> vec;
//   vec.push_back(DeclareXSecMttXSec8());
//   vec.push_back(DeclareXSecYttXSec8());
//   vec.push_back(DeclareXSecYtXSec8());
//   vec.push_back(DeclareXSecPttXSec8());
//   return vec;
// }

// std::vector<UnfoldingIOHandler*> DeclareXSecAll1D()
// {
//   std::vector<UnfoldingIOHandler*> vec;
//   vec.push_back(DeclareXSecPttt1D());
//   vec.push_back(DeclareXSecDetatt1D());
//   vec.push_back(DeclareXSecDphitt1D());
//   vec.push_back(DeclareXSecNj("def", 0.4, 30, 1));
//   vec.push_back(DeclareXSecNj("def", 0.4, 30, 2));
//   vec.push_back(DeclareXSecNj("def", 0.4, 30, 6));
//   return vec;
// }

/*std::vector<UnfoldingIOHandler*> DeclareXSecAll3D()
{
  std::vector<UnfoldingIOHandler*> vec;
  vec.push_back(DeclareXSecNjMttYtt());
  vec.push_back(DeclareXSecNj2MttYtt());
  return vec;
}*/

// std::vector<UnfoldingIOHandler*> DeclareXSecAllPttt2B()
// {
//   std::vector<UnfoldingIOHandler*> vec;
//   vec.push_back(DeclareXSecPttt2B(30.0));
//   vec.push_back(DeclareXSecPttt2B(40.0));
//   vec.push_back(DeclareXSecPttt2B(50.0));
//   vec.push_back(DeclareXSecPttt2B(75.0));
//   vec.push_back(DeclareXSecPttt2B(100.0));
//   vec.push_back(DeclareXSecPttt2B(150.0));
//   return vec;
// }

// std::vector<UnfoldingIOHandler*> DeclareXSecAllPttt2BMttYtt()
// {
//   std::vector<UnfoldingIOHandler*> vec;
//   vec.push_back(DeclareXSecPttt2BMttYtt(30.0));
//   vec.push_back(DeclareXSecPttt2BMttYtt(40.0));
//   vec.push_back(DeclareXSecPttt2BMttYtt(50.0));
//   vec.push_back(DeclareXSecPttt2BMttYtt(75.0));
//   vec.push_back(DeclareXSecPttt2BMttYtt(100.0));
//   vec.push_back(DeclareXSecPttt2BMttYtt(150.0));
//   return vec;
// }

#endif
