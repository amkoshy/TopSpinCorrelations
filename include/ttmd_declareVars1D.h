#ifndef ttmd_declareVars1D_h
#define ttmd_declareVars1D_h

#include "ttmd_defineVars.h"
#include "UnfoldingIOHandler.h"
#include "utils.h"
#include "ttmd_ZMadGraph.h"
#include <cassert>

// ****************************************************************************************************
// ************************* parton level  ************************************************************
// ****************************************************************************************************

// ******************************************************************************************
// ************************* y(ttbar) *******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecYtt1D(PlotterConfigurationHelper *configHelper, const int mode = 0)
{
  ZVarYttSigned* var = new ZVarYttSigned(mode);
  var->Bins.push_back({-2.6, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.6}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.6, -1.8, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.6}); // coarse binning at reco level. Works for all three years and all systematics (including 2016 M2T4). Now also works for loose kin reco
  //var->Bins.push_back({-2.6, -1.9, -1.7, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9, 2.6}); //fine binning at reco level
  var->Title = "y(t#bar{t})";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ytts-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* M(ttbar) *******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecMtt1D(PlotterConfigurationHelper *configHelper, const int mode = 0)
{
  ZVarMtt* var = new ZVarMtt(mode);
  var->Bins.push_back({300, 380, 470, 620, 820, 1100, 1500, 2500}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({300, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 590, 630, 670, 750, 850, 1000, 1500, 2500}); // coarse binning at reco level. Works for all three years and all systematics (including 2016 M2T4)
  //var->Bins.push_back({300,360,365,370,375,380,385,390,395,400,405,410,415,420,425,430,435,440,445,450,455,460,465,470,475,480,485,490,495,500,505,510,515,520,525,530,535,540,545,550,555,560,565,570,575,580,590,600,610,620,630,640,650,660,670,680,700,720,740,760,780,800,850,900,950,1000,1100,1200,1300,1500,2500}); //fine binning at reco level 

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "mtt-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ****************************************************************************************************
// ************************* M(ttbar) threshold *******************************************************
// ****************************************************************************************************

UnfoldingIOHandler* DeclareXSecMtt1D_threshold(PlotterConfigurationHelper *configHelper, const int mode = 9)
{
  ZVarMtt* var = new ZVarMtt(mode);
  var->Bins.push_back({200, 260, 330, 410, 500}); //binning according to resolution gen level
  var->Bins.push_back({200, 250, 280, 310, 340, 370, 400, 420, 440, 460, 480, 500});

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "mtt-threshold-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* pT(ttbar) ******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecPttt1D_new(PlotterConfigurationHelper *configHelper, const int mode = 0)
{
  ZVarPttt* var = new ZVarPttt(mode);
  var->Bins.push_back({0, 40, 100, 200, 310, 420, 570}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0, 10, 20, 30, 40, 50, 80, 100, 120, 140, 160, 180, 200, 230, 260, 310, 360, 420, 570}); // coarse binning at reco level. Works for all three years and all systematics (including 2016 M2T4)
  //var->Bins.push_back({0,5,7,10,13,16,19,22,25,28,30,32,36,40,45,50,55,62,70,75,79,89,100,112,129,140,160,180,200,230,260,300,350,430,570}); //fine binning at reco level  
  var->Title = "p_{T}(t#bar{t})";
  var->TitleTxt = "pT(ttbar)";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "pttt-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* y(t) ***********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecYt1D(PlotterConfigurationHelper *configHelper)
{
  ZVarYtSigned* var = new ZVarYtSigned;
  var->Bins.push_back({-2.6, -1.8, -1.35, -0.90, -0.45, 0.0, 0.45, 0.9, 1.35, 1.8, 2.6}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.6, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.6, -0.45, -0.3, -0.2, 0.0, 0.2, 0.3, 0.45, 0.6, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.6}); // coarse binning at reco level. Works for all three years and all systematics (including 2016 M2T4)
  //var->Bins.push_back({-2.6, -1.9, -1.7, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9, 2.6}); //fine binning at reco level (doesn't work for UETUNE_DOWN when running 2016 M2T4)
  var->Title = "y(t)";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "yts-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* y(tbar) ********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecYat1D(PlotterConfigurationHelper *configHelper)
{
  ZVarYatSigned* var = new ZVarYatSigned;
  var->Bins.push_back({-2.6, -1.8, -1.35, -0.90, -0.45, 0.0, 0.45, 0.9, 1.35, 1.8, 2.6}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.6, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.6, -0.45, -0.3, -0.2, 0.0, 0.2, 0.3, 0.45, 0.6, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.6}); // coarse binning at reco level. Works for all three years and all systematics (including 2016 M2T4)
  //var->Bins.push_back({-2.6, -1.9, -1.7, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9, 2.6}); //fine binning at reco level (doesn't work for UETUNE_DOWN when running 2016 M2T4)
  var->Title = "y(#bar{t})";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "yats-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* y(t) leading ***************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecYtLead1D(PlotterConfigurationHelper *configHelper)
{
  ZVarYtLead* var = new ZVarYtLead;
  var->Bins.push_back({-2.6, -1.65, -1.1, -0.55, 0.0, 0.55, 1.1, 1.65, 2.6}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.6, -1.75, -1.5, -1.25, -1.0, -0.75, -0.6, -0.45, -0.3, -0.2, 0.0, 0.2, 0.3, 0.45, 0.6, 0.75, 1.0, 1.25, 1.5, 1.75, 2.6});
  var->Title = "y(t) leading";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ytLead-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* y(t) trailing **************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecYtNLead1D(PlotterConfigurationHelper *configHelper)
{
  ZVarYtNLead* var = new ZVarYtNLead;
  var->Bins.push_back({-2.6, -1.65, -1.1, -0.55, 0.0, 0.55, 1.1, 1.65, 2.6}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.6, -1.75, -1.5, -1.25, -1.0, -0.75, -0.6, -0.45, -0.3, -0.2, 0.0, 0.2, 0.3, 0.45, 0.6, 0.75, 1.0, 1.25, 1.5, 1.75, 2.6});
  var->Title = "y(t) trailing";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ytNLead-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* pT(t) **********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecPtt1D(PlotterConfigurationHelper *configHelper)
{
  ZVarPtt* var = new ZVarPtt;
  var->Bins.push_back({0, 65, 125, 200, 290, 400, 550}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 230, 250, 270, 290, 310, 330, 370, 430, 490, 550}); // coarse binning at reco level. Works for all three years and all systematics (including 2016 M2T4)
  //var->Bins.push_back({0,15,22,30,36,40,45,51,55,59,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,260,270,280,290,300,310,320,330,350,370,390,420,450,490,550}); //fine binning at reco level  
  var->Title = "p_{T}(t)";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ptt-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* pT(tbar) *******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecPtat1D(PlotterConfigurationHelper *configHelper)
{
  ZVarPtat* var = new ZVarPtat;
  var->Bins.push_back({0, 65, 125, 200, 290, 400, 550}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 230, 250, 270, 290, 310, 330, 370, 430, 490, 550}); // coarse binning at reco level. Works for all three years and all systematics (including 2016 M2T4)
  //var->Bins.push_back({0,15,22,30,36,40,45,51,55,59,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,260,270,280,290,300,310,320,330,350,370,390,420,450,490,550}); //fine binning at reco level  
  var->Title = "p_{T}(#bar{t})";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ptat-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// **************************************************************************************************
// ************************* pT(t) leading **********************************************************
// **************************************************************************************************

UnfoldingIOHandler* DeclareXSecPttLead1D(PlotterConfigurationHelper *configHelper)
{
  ZVarPttLead* var = new ZVarPttLead;
  var->Bins.push_back({0, 65, 125, 200, 290, 400, 550}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 230, 250, 270, 290, 310, 330, 370, 430, 490, 550});
  var->Title = "p_{T}(t) leading";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "pttLead-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ***************************************************************************************************
// ************************* pT(t) trailing **********************************************************
// ***************************************************************************************************

UnfoldingIOHandler* DeclareXSecPttNLead1D(PlotterConfigurationHelper *configHelper)
{
  ZVarPttNLead* var = new ZVarPttNLead;
  var->Bins.push_back({0, 65, 125, 200, 290, 400, 550}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 230, 250, 270, 290, 310, 330, 370, 430, 490, 550});
  var->Title = "p_{T}(t) trailing";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "pttNLead-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// *****************************************************************************************************
// ************************* pT(t) rest frame **********************************************************
// *****************************************************************************************************

UnfoldingIOHandler* DeclareXSecPttTTRestFrame1D(PlotterConfigurationHelper *configHelper)
{
  ZVarPttTTRestFrame* var = new ZVarPttTTRestFrame;
  var->Bins.push_back({0, 65, 125, 200, 290, 400, 550}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 230, 250, 270, 290, 310, 330, 370, 430, 490, 550});
  var->Title = "p_{T}(t) (t#bar{t} RF)";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "pttTTRestFrame-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ********************************************************************************************************
// ************************* pT(tbar) rest frame **********************************************************
// ********************************************************************************************************

UnfoldingIOHandler* DeclareXSecPtatTTRestFrame1D(PlotterConfigurationHelper *configHelper)
{
  ZVarPtatTTRestFrame* var = new ZVarPtatTTRestFrame;
  var->Bins.push_back({0, 65, 125, 200, 290, 400, 550}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 230, 250, 270, 290, 310, 330, 370, 430, 490, 550});
  var->Title = "p_{T}(#bar{t}) (t#bar{t} RF)";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ptatTTRestFrame-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}


// ******************************************************************************************
// ************************* delta_phi(t,tbar) **********************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecDphitt1D_new(PlotterConfigurationHelper *configHelper)
{
  ZVarDphitt* var = new ZVarDphitt;
  var->Bins.push_back({0, 1.57, 2.67, 3.02, 3.142}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0,0.25,0.75,1.0,1.25,1.5,1.8,2.5,2.6,2.65,2.70,2.75,2.80,2.85,2.9,2.95,3.0,3.05,3.08,3.10,3.1415});
  var->Title = "#Delta#phi(t,#bar{t})";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "dphitt-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}


// ******************************************************************************************
// ************************* delta_y(t,tbar) ************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecDytt1D(PlotterConfigurationHelper *configHelper)
{
  ZVarDytt* var = new ZVarDytt;
  var->Bins.push_back({-2.6, -1.4, -0.9, -0.4, 0.0, 0.4, 0.9, 1.4, 2.6}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.6, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.6, -0.45, -0.3, -0.2, 0.0, 0.2, 0.3, 0.45, 0.6, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.6});
  var->Title = "#Delta y(t,#bar{t})";
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "dytt-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* Ptt/Mtt ********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecRatioPttMtt(PlotterConfigurationHelper *configHelper)
{
  ZVarRatioPttMtt* var = new ZVarRatioPttMtt;
  var->Bins.push_back({0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.0});
  var->Bins.push_back({0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.8, 1.0});
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ratioPttMtt-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* Pttt/Mtt *******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecRatioPtttMtt(PlotterConfigurationHelper *configHelper)
{
  ZVarRatioPtttMtt* var = new ZVarRatioPtttMtt;
  var->Bins.push_back({0., 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.5});
  var->Bins.push_back({0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 1.2, 1.5});
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ratioPtttMtt-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* RecoPartonMomFraction  *******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecRecoPartonMomFraction(PlotterConfigurationHelper *configHelper)
{
  ZVarRecoPartonMomFraction* var = new ZVarRecoPartonMomFraction;
  var->Bins.push_back({-3.0, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, 0.});
  var->Bins.push_back({-3.0, -2.6,  -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.3, 0.});
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "RecoPartonMomFraction-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* RecoAntipartonMomFraction *******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecRecoAntipartonMomFraction(PlotterConfigurationHelper *configHelper)
{
  ZVarRecoAntipartonMomFraction* var = new ZVarRecoAntipartonMomFraction;
  var->Bins.push_back({-3.0, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, 0.});
  var->Bins.push_back({-3.0, -2.6, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.3, 0.});
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "RecoAntipartonMomFraction-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* Extra Jets *****************************************************
// ******************************************************************************************
UnfoldingIOHandler* DeclareXSecNj_simple(PlotterConfigurationHelper *configHelper,const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const int nLastBinBoundary = 4)
{
  ZVar* var = new ZVarNj(def, iso, pt);
  var->Bins.push_back({});
  for(int b = 0; b <= nLastBinBoundary + 1; b++)
    var->Bins.back().push_back(-0.5 + b);
  assert(var->Bins.back().back() < 8.0);
  var->Bins.back().push_back(8.5);

  TString name = utils::ttmd::GetNameExtraJetVar("nj", def, iso, pt);
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, name + TString::Format("-xsec%d", nLastBinBoundary));
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

#endif
