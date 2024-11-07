#ifndef DECLAREVARSCP_H
#define DECLAREVARSCP_H

#include "ttmd_defineVars.h"
#include "UnfoldingIOHandler.h"
#include "utils.h"
#include <cassert>

UnfoldingIOHandler* DeclareXSecNvtx(PlotterConfigurationHelper *configHelper)
{
  ZVar* var = new ZVarNvtx;
  var->Bins.push_back({});
  double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 0, 50, 1);
  var->Bins.back().push_back(last);
  assert(IsEqual(var->Bins.back().back(), 50));

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var},configHelper, "nvtx-cp");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->RatMin = 0.0;
  unfoldingIOHandler->RatMax = 2.0;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// b-tagged jet multiplicity
UnfoldingIOHandler* DeclareXSecNbj(PlotterConfigurationHelper *configHelper)
{
  ZVar* var = new ZVarNbj;
  var->Bins.push_back({});
  double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 0, 6, 1);
  var->Bins.back().push_back(last);
  assert(IsEqual(var->Bins.back().back(), 6));

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var},configHelper, "nbj-cp");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

// all jet multiplicity (but with jet-lepton cleaning cut)
UnfoldingIOHandler* DeclareXSecNaj(PlotterConfigurationHelper *configHelper)
{
  ZVar* var = new ZVarNaj;
  var->Bins.push_back({});
  double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 0, 9, 1);
  var->Bins.back().push_back(last);
  assert(IsEqual(var->Bins.back().back(), 9));

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var},configHelper, "naj-cp");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecMet(PlotterConfigurationHelper *configHelper)
{
  ZVar* var = new ZVarMet;
  var->Bins.push_back({});
  double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 0, 30, 10.0);
  var->Bins.back().push_back(last);
  assert(IsEqual(var->Bins.back().back(), 300.0));

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var},configHelper, "met-cp");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecMll(PlotterConfigurationHelper *configHelper)
{
  ZVar* var = new ZVarMll;
  var->Bins.push_back({});
  double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 0, 30, 10.0);
  var->Bins.back().push_back(last);
  assert(IsEqual(var->Bins.back().back(), 300.0));

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var},configHelper, "mll-cp");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecPtDaug(PlotterConfigurationHelper *configHelper, const TString daug)
{
  assert(daug == "j1" || daug == "j2");

  ZVar* var = new ZVarPtDaug(daug);
  var->Bins.push_back({});
  double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), 0, 50, 5.0);
  var->Bins.back().push_back(last);
  assert(IsEqual(var->Bins.back().back(), 250.0));

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var},configHelper, "pt" + daug + "-cp");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecEtaDaug(PlotterConfigurationHelper *configHelper, const TString daug)
{
  assert(daug == "j1" || daug == "j2");

  ZVar* var = new ZVarEtaDaug(daug);
  var->Bins.push_back({});
  double last = utils::ttmd::DVectorAppendNWidth(var->Bins.back(), -2.5, 50, 0.1);
  var->Bins.back().push_back(last);
  assert(IsEqual(var->Bins.back().back(), 2.5));

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({var},configHelper, "eta" + daug + "-cp");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

std::vector<UnfoldingIOHandler*> DeclareXSecAllCP(PlotterConfigurationHelper *configHelper)
{
  std::vector<UnfoldingIOHandler*> vec;
  vec.push_back(DeclareXSecNvtx(configHelper));
  vec.push_back(DeclareXSecNbj(configHelper));
  vec.push_back(DeclareXSecNaj(configHelper));
  vec.push_back(DeclareXSecNjCP(configHelper));
  vec.push_back(DeclareXSecNj_insideTTbarCP(configHelper));
  vec.push_back(DeclareXSecMet(configHelper));
  vec.push_back(DeclareXSecMll(configHelper));
  vec.push_back(DeclareXSecPtDaug(configHelper,"j1"));
  vec.push_back(DeclareXSecEtaDaug(configHelper,"j1"));
  vec.push_back(DeclareXSecPtDaug(configHelper,"j2"));
  vec.push_back(DeclareXSecEtaDaug(configHelper,"j2"));
  vec.push_back(DeclareXSecPtDaug(configHelper,"ej1"));
  vec.push_back(DeclareXSecEtaDaug(configHelper,"ej1"));
  vec.push_back(DeclareXSecPtDaug(configHelper,"ej2"));
  vec.push_back(DeclareXSecEtaDaug(configHelper,"ej2"));
  vec.push_back(DeclareXSecPttPSE25(configHelper));
  vec.push_back(DeclareXSecYtPSE25(configHelper));
  vec.push_back(DeclareXSecPtttPSE25(configHelper));
  vec.push_back(DeclareXSecMttPSE25(configHelper));
  vec.push_back(DeclareXSecYttPSE25(configHelper));
  return vec;
}

#endif // DECLAREVARSCP_H
