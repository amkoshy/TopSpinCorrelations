#ifndef ttmd_declareVars1D_particle_h
#define ttmd_declareVars1D_particle_h

#include "ttmd_defineVars.h"
#include "UnfoldingIOHandler.h"
#include "utils.h"
#include "ttmd_ZMadGraph.h"
#include <cassert>

// ****************************************************************************************************
// ************************* particle level  ************************************************************
// ****************************************************************************************************

// ******************************************************************************************
// ************************* y(ttbar) *******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecYtt1DParticle(PlotterConfigurationHelper *configHelper, const int mode = 0)
{
  ZVarYttSigned* var = new ZVarYttSigned(mode, true);
  //var->Bins.push_back({-2.6, -2.1, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.1, 2.6}); // updated 30.07.20
  //var->Bins.push_back({-2.6, -2.2, -1.8, -1.6, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.2, 2.6}); // updated 30.07.20
  var->Bins.push_back({-2.6, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.6}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.6, -1.8, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.6});
  var->Title = "y(t#bar{t})";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ytts-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* M(ttbar) *******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecMtt1DParticle(PlotterConfigurationHelper *configHelper, const int mode = 0)
{
  ZVarMtt* var = new ZVarMtt(mode, true);
  //var->Bins.push_back({300, 340, 380, 425, 470, 545, 620, 720, 820, 1100, 1500, 2500}); // updated 30.07.20
  //var->Bins.push_back({300, 330, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 590, 610, 630, 650, 670, 710, 750, 800, 850, 1000, 1500, 2500}); // updated 30.07.20
  var->Bins.push_back({300, 380, 470, 620, 820, 1100, 1500, 2500}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({300, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 590, 630, 670, 750, 850, 1000, 1500, 2500});

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "mtt-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* pT(ttbar) ******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecPttt1D_new_particle(PlotterConfigurationHelper *configHelper, const int mode = 0)
{
  ZVarPttt* var = new ZVarPttt(mode, true);
  //var->Bins.push_back({0, 20, 40, 70, 100, 150, 200, 310, 420, 570}); // updated 30.07.20
  //var->Bins.push_back({0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 65, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 230, 260, 310, 360, 420, 570}); // updated 30.07.20
  var->Bins.push_back({0, 40, 100, 200, 310, 420, 570}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0, 10, 20, 30, 40, 50, 80, 100, 120, 140, 160, 180, 200, 230, 260, 310, 360, 420, 570});
  var->Title = "p_{T}(t#bar{t})";
  var->TitleTxt = "pT(ttbar)";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "pttt-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* y(t) ***********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecYt1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarYtSigned* var = new ZVarYtSigned(true);
  //var->Bins.push_back({-2.6, -2.2, -1.8, -1.575, -1.35, -1.225, -0.9, -0.675, -0.45, -0.225, 0.0, 0.225, 0.45, 0.675, 0.9, 1.225, 1.35, 1.575, 1.8, 2.2, 2.6}); // updated 30.07.20
  //var->Bins.push_back({-2.6, -2.3, -2.0, -1.875, -1.75, -1.625, -1.5, -1.375, -1.25, -1.125, -1.0, -0.875, -0.75, -0.675, -0.6, -0.525, -0.45, -0.375, -0.3, -0.25, -0.2, -0.1, 0.0, 0.1, 0.2, 0.25, 0.3, 0.375, 0.45, 0.525, 0.6, 0.675, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2.0, 2.3, 2.6}); // updated 30.07.20
  var->Bins.push_back({-2.6, -1.8, -1.35, -0.90, -0.45, 0.0, 0.45, 0.9, 1.35, 1.8, 2.6}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.6, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.6, -0.45, -0.3, -0.2, 0.0, 0.2, 0.3, 0.45, 0.6, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.6});
  var->Title = "y(t)";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "yts-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* y(tbar) ********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecYat1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarYatSigned* var = new ZVarYatSigned(true);
  //var->Bins.push_back({-2.6, -2.2, -1.8, -1.575, -1.35, -1.225, -0.9, -0.675, -0.45, -0.225, 0.0, 0.225, 0.45, 0.675, 0.9, 1.225, 1.35, 1.575, 1.8, 2.2, 2.6}); // updated 30.07.20
  //var->Bins.push_back({-2.6, -2.3, -2.0, -1.875, -1.75, -1.625, -1.5, -1.375, -1.25, -1.125, -1.0, -0.875, -0.75, -0.675, -0.6, -0.525, -0.45, -0.375, -0.3, -0.25, -0.2, -0.1, 0.0, 0.1, 0.2, 0.25, 0.3, 0.375, 0.45, 0.525, 0.6, 0.675, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2.0, 2.3, 2.6}); //updated 30.07.20
  var->Bins.push_back({-2.6, -1.8, -1.35, -0.90, -0.45, 0.0, 0.45, 0.9, 1.35, 1.8, 2.6}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.6, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.6, -0.45, -0.3, -0.2, 0.0, 0.2, 0.3, 0.45, 0.6, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.6});
  var->Title = "y(#bar{t})";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "yats-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* y(t) leading ***************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecYtLead1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarYtLead* var = new ZVarYtLead(true);
  //var->Bins.push_back({-2.6, -2.125, -1.65, -1.375, -1.1, -0.825, -0.55, -0.275, 0.0, 0.275, 0.55, 0.825, 1.1, 1.375, 1.65, 2.125, 2.6}); // updated 30.07.20
  //var->Bins.push_back({-2.6, -2.175, -1.75, -1.625, -1.5, -1.375, -1.25, -1.125, -1.0, -0.875, -0.75, -0.675, -0.6, -0.525, -0.45, -0.375, -0.3, -0.25, -0.2, -0.1, 0.0, 0.1, 0.2, 0.25, 0.3, 0.375, 0.45, 0.525, 0.6, 0.675, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 2.175, 2.6}); // updated 30.07.20
  var->Bins.push_back({-2.6, -1.65, -1.1, -0.55, 0.0, 0.55, 1.1, 1.65, 2.6}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.6, -1.75, -1.5, -1.25, -1.0, -0.75, -0.6, -0.45, -0.3, -0.2, 0.0, 0.2, 0.3, 0.45, 0.6, 0.75, 1.0, 1.25, 1.5, 1.75, 2.6});
  var->Title = "y(t) leading";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ytLead-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* y(t) trailing **************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecYtNLead1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarYtNLead* var = new ZVarYtNLead(true);
  //var->Bins.push_back({-2.6, -2.125, -1.65, -1.375, -1.1, -0.825, -0.55, -0.275, 0.0, 0.275, 0.55, 0.825, 1.1, 1.375, 1.65, 2.125, 2.6}); // updated 30.07.20
  //var->Bins.push_back({-2.6, -2.175, -1.75, -1.625, -1.5, -1.375, -1.25, -1.125, -1.0, -0.875, -0.75, -0.675, -0.6, -0.525, -0.45, -0.375, -0.3, -0.25, -0.2, -0.1, 0.0, 0.1, 0.2, 0.25, 0.3, 0.375, 0.45, 0.525, 0.6, 0.675, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 2.175, 2.6}); // updated 30.07.20
  var->Bins.push_back({-2.6, -1.65, -1.1, -0.55, 0.0, 0.55, 1.1, 1.65, 2.6}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.6, -1.75, -1.5, -1.25, -1.0, -0.75, -0.6, -0.45, -0.3, -0.2, 0.0, 0.2, 0.3, 0.45, 0.6, 0.75, 1.0, 1.25, 1.5, 1.75, 2.6});
  var->Title = "y(t) trailing";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ytNLead-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* pT(t) **********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecPtt1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarPtt* var = new ZVarPtt(true);
  //var->Bins.push_back({0, 32.5, 65, 95, 125, 162.5, 200, 290, 400, 550}); //updated 30.07.20
  //var->Bins.push_back({0, 7.5, 15, 22.5, 30, 37.5, 45, 52.5, 60, 67.5, 75, 82.5, 90, 97.5, 105, 112.5, 120, 127.5, 135, 142.5, 150, 157.5, 165, 172.5, 180, 187.5, 195, 202.5, 210, 230, 250, 270, 290, 310, 330, 370, 430, 490, 550}); //updated 30.07.20
  var->Bins.push_back({0, 65, 125, 200, 290, 400, 550}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 230, 250, 270, 290, 310, 330, 370, 430, 490, 550});
  var->Title = "p_{T}(t)";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ptt-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* pT(tbar) *******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecPtat1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarPtat* var = new ZVarPtat(true);
  //var->Bins.push_back({0, 32.5, 65, 95, 125, 162.5, 200, 290, 400, 550});// updated 30.07.20
  //var->Bins.push_back({0, 7.5, 15, 22.5, 30, 37.5, 45, 52.5, 60, 67.5, 75, 82.5, 90, 97.5, 105, 112.5, 120, 127.5, 135, 142.5, 150, 157.5, 165, 172.5, 180, 187.5, 195, 202.5, 210, 230, 250, 270, 290, 310, 330, 370, 430, 490, 550}); // updated 30.07.20
  var->Bins.push_back({0, 65, 125, 200, 290, 400, 550}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 230, 250, 270, 290, 310, 330, 370, 430, 490, 550});
  var->Title = "p_{T}(#bar{t})";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ptat-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// **************************************************************************************************
// ************************* pT(t) leading **********************************************************
// **************************************************************************************************

UnfoldingIOHandler* DeclareXSecPttLead1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarPttLead* var = new ZVarPttLead(true);
  //var->Bins.push_back({0, 32.5, 65, 95, 125, 162.5, 200, 290, 400, 550}); // updated 30.07.20
  //var->Bins.push_back({0, 7.5, 15, 22.5, 30, 37.5, 45, 52.5, 60, 67.5, 75, 82.5, 90, 97.5, 105, 112.5, 120, 127.5, 135, 142.5, 150, 157.5, 165, 172.5, 180, 187.5, 195, 202.5, 210, 230, 250, 270, 290, 310, 330, 370, 430, 490, 550}); // updated 30.07.20
  var->Bins.push_back({0, 65, 125, 200, 290, 400, 550}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 230, 250, 270, 290, 310, 330, 370, 430, 490, 550});
  var->Title = "p_{T}(t) leading";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "pttLead-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ***************************************************************************************************
// ************************* pT(t) trailing **********************************************************
// ***************************************************************************************************

UnfoldingIOHandler* DeclareXSecPttNLead1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarPttNLead* var = new ZVarPttNLead(true);
  //var->Bins.push_back({0, 32.5, 65, 95, 125, 162.5, 200, 290, 400, 550}); //updated 30.07.20
  //var->Bins.push_back({0, 7.5, 15, 22.5, 30, 37.5, 45, 52.5, 60, 67.5, 75, 82.5, 90, 97.5, 105, 112.5, 120, 127.5, 135, 142.5, 150, 157.5, 165, 172.5, 180, 187.5, 195, 202.5, 210, 230, 250, 270, 290, 310, 330, 370, 430, 490, 550}); // updated 30.07.20
  var->Bins.push_back({0, 65, 125, 200, 290, 400, 550}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 230, 250, 270, 290, 310, 330, 370, 430, 490, 550});
  var->Title = "p_{T}(t) trailing";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "pttNLead-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// *****************************************************************************************************
// ************************* pT(t) rest frame **********************************************************
// *****************************************************************************************************

UnfoldingIOHandler* DeclareXSecPttTTRestFrame1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarPttTTRestFrame* var = new ZVarPttTTRestFrame(true);
  //var->Bins.push_back({0, 32.5, 65, 95, 125, 162.5, 200, 290, 400, 550}); // updated 30.07.20
  //var->Bins.push_back({0, 7.5, 15, 22.5, 30, 37.5, 45, 52.5, 60, 67.5, 75, 82.5, 90, 97.5, 105, 112.5, 120, 127.5, 135, 142.5, 150, 157.5, 165, 172.5, 180, 187.5, 195, 202.5, 210, 230, 250, 270, 290, 310, 330, 370, 430, 490, 550}); //updated 30.07.20
  var->Bins.push_back({0, 65, 125, 200, 290, 400, 550}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 230, 250, 270, 290, 310, 330, 370, 430, 490, 550});
  var->Title = "p_{T}(t) (t#bar{t} RF)";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "pttTTRestFrame-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ********************************************************************************************************
// ************************* pT(tbar) rest frame **********************************************************
// ********************************************************************************************************

UnfoldingIOHandler* DeclareXSecPtatTTRestFrame1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarPtatTTRestFrame* var = new ZVarPtatTTRestFrame(true);
  //var->Bins.push_back({0, 32.5, 65, 95, 125, 162.5, 200, 290, 400, 550}); // updated 30.07.20
  //var->Bins.push_back({0, 7.5, 15, 22.5, 30, 37.5, 45, 52.5, 60, 67.5, 75, 82.5, 90, 97.5, 105, 112.5, 120, 127.5, 135, 142.5, 150, 157.5, 165, 172.5, 180, 187.5, 195, 202.5, 210, 230, 250, 270, 290, 310, 330, 370, 430, 490, 550}); // updated 30.07.20
  var->Bins.push_back({0, 65, 125, 200, 290, 400, 550}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 230, 250, 270, 290, 310, 330, 370, 430, 490, 550});
  var->Title = "p_{T}(#bar{t}) (t#bar{t} RF)";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ptatTTRestFrame-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}


// ******************************************************************************************
// ************************* delta_phi(t,tbar) **********************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecDphitt1D_new_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarDphitt* var = new ZVarDphitt(true);
  var->Bins.push_back({0, 1.57, 2.67, 3.02, 3.142}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0,0.25,0.75,1.0,1.25,1.5,1.8,2.5,2.6,2.65,2.70,2.75,2.80,2.85,2.9,2.95,3.0,3.05,3.08,3.10,3.1415});
  var->Title = "#Delta#phi(t,#bar{t})";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "dphitt-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}


// ******************************************************************************************
// ************************* delta_y(t,tbar) ************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecDytt1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarDytt* var = new ZVarDytt(true);
  var->Bins.push_back({-2.6, -1.4, -0.9, -0.4, 0.0, 0.4, 0.9, 1.4, 2.6}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.6, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.6, -0.45, -0.3, -0.2, 0.0, 0.2, 0.3, 0.45, 0.6, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.6});
  var->Title = "#Delta y(t,#bar{t})";
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "dytt-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* pT(l) **********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecPtl1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarPtl* var = new ZVarPtl();
  var->Bins.push_back({20, 40, 70, 120, 180, 400}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({20, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 230, 250, 270, 290, 310, 330, 370, 400});
  var->Title = "p_{T}(l)";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ptl-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* pT(lbar) *******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecPtal1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarPtal* var = new ZVarPtal();
  var->Bins.push_back({20, 40, 70, 120, 180, 400}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({20, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 230, 250, 270, 290, 310, 330, 370, 400});
  var->Title = "p_{T}(#bar{l})";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ptal-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// **************************************************************************************************
// ************************* pT(l) leading **********************************************************
// **************************************************************************************************

UnfoldingIOHandler* DeclareXSecPtlLead1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarPtlLead* var = new ZVarPtlLead();
  var->Bins.push_back({20, 40, 70, 120, 180, 400}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({20, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 230, 250, 270, 290, 310, 330, 370, 400});
  var->Title = "p_{T}(l) leading";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ptlLead-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ***************************************************************************************************
// ************************* pT(l) trailing **********************************************************
// ***************************************************************************************************

UnfoldingIOHandler* DeclareXSecPtlNLead1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarPtlNLead* var = new ZVarPtlNLead();

  var->Title = "p_{T}(l) trailing";
  var->Units = "GeV";
  var->Bins.push_back({20, 35, 50, 90, 140, 400}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({20, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 230, 250, 270, 290, 310, 330, 370, 400});
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ptlNLead-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* eta(l) ***********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecEtal1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarEtal* var = new ZVarEtal();
  var->Bins.push_back({-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.4, -2.25, -2.1, -1.95, -1.80, -1.65, -1.5, -1.35, -1.2, -1.1, -0.9, -0.75, -0.6, -0.45, -0.3, -0.15, 0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.1, 1.2, 1.35, 1.5, 1.65, 1.8, 1.95, 2.1, 2.25, 2.4});
  var->Title = "#eta(l)";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "etal-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* eta(lbar) ********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecEtaal1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarEtaal* var = new ZVarEtaal();
  var->Bins.push_back({-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.4, -2.25, -2.1, -1.95, -1.80, -1.65, -1.5, -1.35, -1.2, -1.1, -0.9, -0.75, -0.6, -0.45, -0.3, -0.15, 0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.1, 1.2, 1.35, 1.5, 1.65, 1.8, 1.95, 2.1, 2.25, 2.4});
  var->Title = "#eta(#bar{l})";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "etaal-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* eta(l) leading ***************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecEtalLead1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarEtalLead* var = new ZVarEtalLead();
  var->Bins.push_back({-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.4, -2.25, -2.1, -1.95, -1.80, -1.65, -1.5, -1.35, -1.2, -1.1, -0.9, -0.75, -0.6, -0.45, -0.3, -0.15, 0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.1, 1.2, 1.35, 1.5, 1.65, 1.8, 1.95, 2.1, 2.25, 2.4});
  var->Title = "#eta(l) leading";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "etalLead-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* eta(l) trailing **************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecEtalNLead1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarEtalNLead* var = new ZVarEtalNLead();
  var->Bins.push_back({-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.4, -2.25, -2.1, -1.95, -1.80, -1.65, -1.5, -1.35, -1.2, -1.1, -0.9, -0.75, -0.6, -0.45, -0.3, -0.15, 0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.1, 1.2, 1.35, 1.5, 1.65, 1.8, 1.95, 2.1, 2.25, 2.4});
  var->Title = "#eta(l) trailing";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "etalNLead-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* delta_phi(l,lbar) **********************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecDphill1D_new_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarDphill* var = new ZVarDphill();
  var->Bins.push_back({0, 0.4, 0.78, 1.14, 1.48, 1.8, 2.1, 2.38, 2.64, 2.89, 3.142}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0, 0.2, 0.4, 0.59, 0.78, 0.96, 1.14, 1.31, 1.48, 1.64, 1.8, 1.95, 2.1, 2.24, 2.38, 2.51, 2.64, 2.77, 2.89, 3.016, 3.142});
  var->Title = "#Delta#phi(l,#bar{l})";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "dphill-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}


// ******************************************************************************************
// ************************* delta_eta(l,lbar) ************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecDEtall1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarDetall* var = new ZVarDetall();
  var->Bins.push_back({-2.4, -1.7, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.7, 2.4}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.4, -2.1, -1.7, -1.45, -1.2, -1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.45, 1.7, 2.1, 2.4});
  var->Title = "#Delta #eta(l,#bar{l})";
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "detall-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* pT(llbar) ******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecPtll1D_new_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarPtll* var = new ZVarPtll();
  var->Bins.push_back({0, 10, 20, 40, 60, 100, 150, 400}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 125, 150, 275, 400});
  var->Title = "p_{T}(l#bar{l})";
  var->TitleTxt = "pT(llbar)";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ptll-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* M(ll) **********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecMll1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarMll* var = new ZVarMll();
  var->Bins.push_back({20, 30, 50, 76, 106, 130, 170, 260, 650});
  var->Bins.push_back({20, 25, 30, 40, 50, 63, 76, 91, 106, 118, 130, 150, 170, 215, 260, 455, 650});

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "mll-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* M(llbb) **********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecMllbb1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarMllbb* var = new ZVarMllbb();
  var->Bins.push_back({90, 200, 300, 400, 510, 630, 760, 1000, 2500}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({90, 145, 200, 250, 300, 350, 400, 455, 510, 570, 630, 695, 760, 930, 1000, 1500, 2500});

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "mllbb-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* M(llbbmet) **********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecMllbbmet1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarMllbbmet* var = new ZVarMllbbmet();
  var->Bins.push_back({90, 200, 300, 400, 510, 630, 760, 1000, 2500}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({90, 145, 200, 250, 300, 350, 400, 455, 510, 570, 630, 695, 760, 930, 1000, 1500, 2500});

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "mllbbmet-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// **************************************************************************************************
// ************************* pT(b) leading **********************************************************
// **************************************************************************************************

UnfoldingIOHandler* DeclareXSecPtbLead1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarPtbLead* var = new ZVarPtbLead();
  var->Bins.push_back({30, 60, 95, 150, 230, 500}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({30, 45, 60, 77.5, 95, 122.5, 150, 190, 230, 365, 500});
  var->Title = "p_{T}(b) leading";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ptbLead-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ***************************************************************************************************
// ************************* pT(b) trailing **********************************************************
// ***************************************************************************************************

UnfoldingIOHandler* DeclareXSecPtbNLead1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarPtbNLead* var = new ZVarPtbNLead();
  var->Bins.push_back({30, 45, 70, 110, 170, 500}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({30, 37.5, 45, 57.5, 70, 90, 110, 140, 170, 335, 500});
  var->Title = "p_{T}(b) trailing";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ptbNLead-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* eta(b) leading ***************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecEtabLead1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarEtabLead* var = new ZVarEtabLead();
  var->Bins.push_back({-2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4});
  var->Title = "#eta(b) leading";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "etabLead-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* eta(b) trailing **************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecEtabNLead1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarEtabNLead* var = new ZVarEtabNLead();
  var->Bins.push_back({-2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4});
  var->Title = "#eta(b) trailing";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "etabNLead-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* pT(bb) ******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecPtbb1D_new_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarPtbb* var = new ZVarPtbb();
  var->Bins.push_back({0, 30, 60, 100, 180, 400}); //binning at gen level corresponds to TOP-17-014
  var->Bins.push_back({0, 15, 30, 45, 60, 80, 100, 140, 180, 290, 400});
  var->Title = "p_{T}(b#bar{b})";
  var->TitleTxt = "pT(bbbar)";
  var->Units = "GeV";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ptbbar-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* M(bb) **********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecMbb1D_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarMbb* var = new ZVarMbb();
  var->Bins.push_back({0, 60, 120, 240, 650});
  var->Bins.push_back({0, 30, 60, 90, 120, 180, 240, 445, 650});

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "mbb-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* Ptt/Mtt ********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecRatioPttMtt_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarRatioPttMtt* var = new ZVarRatioPttMtt(true);
  var->Bins.push_back({0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.0});
  var->Bins.push_back({0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.8, 1.0});
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ratioPttMtt-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* Pttt/Mtt *******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecRatioPtttMtt_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarRatioPtttMtt* var = new ZVarRatioPtttMtt(true);
  var->Bins.push_back({0., 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.5});
  var->Bins.push_back({0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 1.2, 1.5});
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ratioPtttMtt-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* RecoPartonMomFraction  *******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecRecoPartonMomFraction_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarRecoPartonMomFraction* var = new ZVarRecoPartonMomFraction(true);
  var->Bins.push_back({-3.0, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, 0.});
  var->Bins.push_back({-3.0, -2.6,  -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.3, 0.});
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "RecoPartonMomFraction-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* RecoAntipartonMomFraction *******************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecRecoAntipartonMomFraction_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarRecoAntipartonMomFraction* var = new ZVarRecoAntipartonMomFraction(true);
  var->Bins.push_back({-3.0, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, 0.});
  var->Bins.push_back({-3.0, -2.6, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.3, 0.});
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "RecoAntipartonMomFraction-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* Pt(b leading)/Pt(t) ********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecRatioPtbLeadt_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarRatioPtbLeadt* var = new ZVarRatioPtbLeadt;
  var->Bins.push_back({0., 0.2, 0.3, 0.4, 0.6, 0.8, 1.2, 2.0});
  var->Bins.push_back({0., 0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.6, 2.0});
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ratioPtbLeadt-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* Pt(b trailing)/Pt(t) ********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecRatioPtbNLeadt_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarRatioPtbNLeadt* var = new ZVarRatioPtbNLeadt;
  var->Bins.push_back({0., 0.2, 0.3, 0.4, 0.6, 0.8, 1.2, 2.0});
  var->Bins.push_back({0., 0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.6, 2.0});
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ratioPtbNLeadt-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// ******************************************************************************************
// ************************* Pt(l)/Pt(t) ********************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecRatioPtlt_particle(PlotterConfigurationHelper *configHelper)
{
  ZVarRatioPtlt* var = new ZVarRatioPtlt;
  var->Bins.push_back({0., 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.5});
  var->Bins.push_back({0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 1.2, 1.5});
  var->Units = "";

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var}, configHelper, "ratioPtlt-particle-xsec8");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->DoParticle = true;

  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

#endif
