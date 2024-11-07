#ifndef declareVars2D_h
#define declareVars2D_h

#include "ttmd_defineVars.h"
#include "UnfoldingIOHandler.h"
#include "utils.h"
#include <cassert>

// ******************************************************************************************
// ************************* M(ttbar),y(ttbar) **********************************************
// ******************************************************************************************

/*UnfoldingIOHandler* DeclareXSecMttYtt()
{
  ZVar* varMtt = new ZVarMtt;
  //varMtt->Expression = "mtt";
  varMtt->Title = "M(t#bar{t})";
  varMtt->Units = "GeV";
  varMtt->NDigits = 0;
  varMtt->Bins.push_back({340, 400, 500, 650, 1500});
  varMtt->Bins.push_back({340,366,373,380,386,393,400,410,420,431,442,454,467,482,500,519,540,563,592,620,650,685,740,850,1500});

  ZVarYtt* varYtt = new ZVarYtt;
  //varYtt->Expression = "ytt";
  varYtt->Title = "y(t#bar{t})";
  varYtt->Units = "";
  varYtt->NDigits = 2;
  varYtt->Bins.push_back({0.0, 0.35, 0.75, 1.15, 2.5});
  varYtt->Bins.push_back({0.0,0.04,0.08,0.12,0.16,0.2,0.25,0.3,0.35,0.41,0.47,0.54,0.61,0.68,0.75,0.83,0.92,1.02,1.15,1.28,1.41,2.5});

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varMtt, varYtt}, "mttytt");
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}*/

// ******************************************************************************************
// ************************* TOP-14-013 *****************************************************
// ******************************************************************************************

UnfoldingIOHandler* DeclareXSecYtPttXSec8(PlotterConfigurationHelper *configHelper, int binning_option = 0)
{
  ZVarYtt* varYtt = new ZVarYtt;
  varYtt->Title = "|y(t)|";
  varYtt->TitleTxt = "y(t)";
  varYtt->Expression = "yt";
  varYtt->Units = "";
  varYtt->NDigits = 2;
  varYtt->Bins.push_back({0.0, 0.35, 0.85, 1.45, 2.5});
  // Fine binning
  if (binning_option == 0) varYtt->Bins.push_back({0.0,0.05,0.10,0.16,0.22,0.28,0.35,0.42,0.49,0.57,0.66,0.75,0.85,0.94,1.03,1.12,1.22,1.33,1.45,1.7,2.5});
  // Thick binning
  else if (binning_option == 1) varYtt->Bins.push_back({0.0,0.10,0.22,0.35,0.49,0.66,0.85,1.03,1.22,1.45,1.7,2.5});

  ZVar* varPtt = new ZVarPtt;
  varPtt->Title = "p_{T}(t)";
  varPtt->TitleTxt = "pT(t)";
  varPtt->Expression = "ptt";
  varPtt->Units = "GeV";
  varPtt->NDigits = 0;
  varPtt->Bins.push_back({0, 80, 150, 250, 600});
  if (binning_option == 0) varPtt->Bins.push_back({0,22,36,51,64,72,80,93,104,115,125,133,141,150,170,190,220,250,290,340,600});
  else if (binning_option == 1) varPtt->Bins.push_back({0,22,36,51,64,72,80,93,104,115,125,133,141,150,170,190,220,250,290,340,600});

  TString xsec_name = "ytptt-xsec8";
  if(binning_option!=0) xsec_name = xsec_name + "-binOpt" + std::to_string(binning_option);

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varYtt, varPtt}, configHelper, xsec_name);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->VMGHisto = { std::pair<std::vector<int>, int>({11}, 0),
                     std::pair<std::vector<int>, int>({12}, 0),
                     std::pair<std::vector<int>, int>({13}, 0),
                     std::pair<std::vector<int>, int>({14}, 0) };
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

ZVar* DecalteVarMttXSec8(const int mode = 0, int binning_option = 0)
{
  ZVar* varMtt = new ZVarMtt(mode);
  varMtt->Bins.push_back({300, 400, 500, 650, 1500});
  if(mode == 0 && binning_option == 0) varMtt->Bins.push_back({300,366,373,380,386,393,400,410,420,431,442,454,467,482,500,519,540,563,592,620,650,685,740,850,1500});
  else if (mode == 0 && binning_option == 1) varMtt->Bins.push_back({300,373,386,400,420,442,467,500,540,592,650,740,1500});
  else if((mode == 6 || mode == 9) && binning_option == 0)
    varMtt->Bins.push_back({0,315,330,345,360,372,382,400,420,440,460,480,500,520,540,570,600,650,700,780,900,2500});
  else if(mode == 2 && binning_option == 0)
    varMtt->Bins.push_back({0,290,310,330,345,360,380,400,420,440,460,480,500,520,540,570,600,650,720,850,2500});
  else
    assert(0);
  return varMtt;
}

ZVarNj* DeclareVarNj(const int njNBins, const TString& def = "def", const double iso = 0.4, const double pt = 30.0){
    ZVarNj* varNj = new ZVarNj(def, iso, pt);
    varNj->NDigits = 1;

    TString binning_filename = "";

    if(njNBins == 2) binning_filename = "data/binnings/multiDiffXSecs/nj2.txt";
    else if(njNBins == 3) binning_filename = "data/binnings/multiDiffXSecs/nj3.txt";
    else if(njNBins == 4) binning_filename = "data/binnings/multiDiffXSecs/nj4.txt";
    else throw std::logic_error(TString::Format("Error in DeclareXSecNjMttYtt(): njNBins = %d not supported", njNBins));

    std::pair<std::vector<double>,std::vector<double>> binning_info = UnfoldingIOHandler::GetBinningInfoFromFile(binning_filename);
    varNj->Bins.push_back(binning_info.first);
    varNj->Bins.push_back(binning_info.second);

    /*if(njNBins == 2)
      varNj->Bins.push_back({-0.5, 0.5, 8.5});
    else if(njNBins == 3)
      varNj->Bins.push_back({-0.5, 0.5, 1.5, 8.5});
    else if(njNBins == 4)
      varNj->Bins.push_back({-0.5, 0.5, 1.5, 2.5, 8.5});
    else
      throw std::logic_error(TString::Format("Error in DeclareXSecNjMttYtt(): njNBins = %d not supported", njNBins));

    if(njNBins == 4) varNj->Bins.push_back({-0.5, 0.5, 1.5, 2.5, 8.5});
    else varNj->Bins.push_back({-0.5, 0.5, 1.5, 8.5});*/

    return varNj;
}

UnfoldingIOHandler* DeclareXSecMttYtXSec8(PlotterConfigurationHelper *configHelper, int binning_option = 0)
{
  ZVar* varMtt = DecalteVarMttXSec8(0, binning_option);

  ZVarYt* varYt = new ZVarYt;
  varYt->Title = "|y(t)|";
  varYt->TitleTxt = "y(t)";
  varYt->Expression = "yt";
  varYt->Units = "";
  varYt->NDigits = 2;
  varYt->NDiv = 303;
  varYt->Bins.push_back({0.0, 0.35, 0.85, 1.45, 2.5});
  if (binning_option == 0) varYt->Bins.push_back({0.0,0.05,0.10,0.16,0.22,0.28,0.35,0.42,0.49,0.57,0.66,0.75,0.85,0.94,1.03,1.12,1.22,1.33,1.45,1.7,2.5}); //TOP-18-004
  // small number of bins in rec level for testing
  else if (binning_option == 1) varYt->Bins.push_back({0.0,0.10,0.16,0.28,0.42,0.57,0.75,0.94,1.12,1.33,1.7,2.5});
  else {
      std::cerr << "ERROR in DeclareXSecMttYtXSec8 -->> binning_option not implemented: " << binning_option << std::endl;
      exit(1);
  }

  TString xsec_name = "mttyt-xsec8";
  if(binning_option!=0) xsec_name = xsec_name + "-binOpt" + std::to_string(binning_option);

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varMtt, varYt}, configHelper, xsec_name);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->VMGHisto = { std::pair<std::vector<int>, int>({21}, 0),
                     std::pair<std::vector<int>, int>({22}, 0),
                     std::pair<std::vector<int>, int>({23}, 0),
                     std::pair<std::vector<int>, int>({24}, 0) };
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecMttYttXSec8(PlotterConfigurationHelper *configHelper, const int binning_opt = 0, const int mode = 0, const int nsubbins = 1)
{
  ZVar* varMtt = DecalteVarMttXSec8(mode, binning_opt);
  if(nsubbins != 1)
    varMtt->Bins[0] = utils::ttmd::MakeSubbins(varMtt->Bins[0], nsubbins);

  ZVarYtt* varYtt = new ZVarYtt(mode);
  varYtt->NDiv = 303;
  varYtt->Bins.push_back({0.0, 0.35, 0.75, 1.15, 2.5});
  if(nsubbins != 1)
    varYtt->Bins.back() = utils::ttmd::MakeSubbins(varYtt->Bins.back(), nsubbins);
  if(mode == 0 || mode == 6 || mode == 9){
    if (binning_opt == 0) varYtt->Bins.push_back({0.0,0.04,0.08,0.12,0.16,0.2,0.25,0.3,0.35,0.41,0.47,0.54,0.61,0.68,0.75,0.83,0.92,1.02,1.15,1.28,1.41,2.5});
    else if (binning_opt == 1) varYtt->Bins.push_back({0.0,0.07,0.15,0.24,0.34,0.46,0.60,0.74,0.91,1.10,2.5});
  }
  else assert(0);
  TString xsec_name = "mttytt";
  if(mode!=0) xsec_name = xsec_name + "-mode" + std::to_string(mode);
  if(binning_opt!=0) xsec_name = xsec_name + "-binOpt" + std::to_string(binning_opt);
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varMtt, varYtt}, configHelper, xsec_name);
  if(nsubbins != 1)
    unfoldingIOHandler->Suffix += TString::Format("-nsubb%d", nsubbins);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  if(nsubbins == 1)
    unfoldingIOHandler->VMGHisto = { std::pair<std::vector<int>, int>({31}, 0),
                       std::pair<std::vector<int>, int>({32}, 0),
                       std::pair<std::vector<int>, int>({33}, 0),
                       std::pair<std::vector<int>, int>({34}, 0) };
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecMttDetattXSec8(PlotterConfigurationHelper *configHelper, int binning_option = 0)
{
  ZVar* varMtt = DecalteVarMttXSec8(0, binning_option);

  ZVarDetatt* varDetatt = new ZVarDetatt;
  varDetatt->Title = "#Delta#eta(t,#bar{t})";
  varDetatt->TitleTxt = "delta_eta(t,tbar)";
  varDetatt->Expression = "detatt";
  varDetatt->Units = "";
  varDetatt->NDigits = 1;
  varDetatt->Bins.push_back({0.0, 0.4, 1.2, 6.0});
  if (binning_option == 0) varDetatt->Bins.push_back({0.0,0.12,0.25,0.4,0.55,0.70,0.85,1.0,1.2,1.4,1.65,2.0,2.5,6.0});
  else if (binning_option == 1) varDetatt->Bins.push_back({0.0,0.25,0.55,0.85,1.2,1.65,2.5,6.0});

  TString xsec_name = "mttdetatt";
  if(binning_option!=0) xsec_name = xsec_name + "-binOpt" + std::to_string(binning_option);

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varMtt, varDetatt}, configHelper, xsec_name);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->VMGHisto = { std::pair<std::vector<int>, int>({41}, 0),
                     std::pair<std::vector<int>, int>({42}, 0),
                     std::pair<std::vector<int>, int>({43}, 0),
                     std::pair<std::vector<int>, int>({44}, 0) };
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecMttPtttXSec8(PlotterConfigurationHelper *configHelper, int binning_option = 0, const int mode = 0)
{
  ZVar* varMtt = DecalteVarMttXSec8(mode, binning_option);

  ZVarPttt* varPttt = new ZVarPttt(mode);
  varPttt->Title = "p_{T}(t#bar{t})";
  varPttt->TitleTxt = "pT(ttbar)";
  varPttt->Units = "GeV";
  varPttt->NDigits = 0;
  varPttt->NDiv = 303;
  varPttt->Bins.push_back({0, 30, 75, 150, 500});
  if(binning_option == 0) varPttt->Bins.push_back({0,7,10,13,16,19,22,25,28,32,36,40,45,50,55,62,70,79,89,100,112,129,150,195,500});
  else if(binning_option == 1) varPttt->Bins.push_back({0,10,16,22,28,36,45,55,70,89,112,150,500});

  TString xsec_name = "mttpttt-xsec8";
  if(binning_option!=0) xsec_name = xsec_name + "-binOpt" + std::to_string(binning_option);

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varMtt, varPttt}, configHelper, xsec_name);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->VMGHisto = { std::pair<std::vector<int>, int>({61}, 0),
                     std::pair<std::vector<int>, int>({62}, 0),
                     std::pair<std::vector<int>, int>({63}, 0),
                     std::pair<std::vector<int>, int>({64}, 0) };
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecYttPtttXSec8(PlotterConfigurationHelper *configHelper, int binning_option = 0, const int mode = 0, const int nsubbins = 1)
{
  ZVarYtt* varYtt = new ZVarYtt(mode);
  varYtt->NDiv = 303;
  varYtt->Bins.push_back({0.0, 0.35, 0.75, 1.15, 2.5});
  if(nsubbins != 1)
    varYtt->Bins.back() = utils::ttmd::MakeSubbins(varYtt->Bins.back(), nsubbins);
  if(mode == 0 || mode == 6 || mode == 9){
    if (binning_option == 0) varYtt->Bins.push_back({0.0,0.04,0.08,0.12,0.16,0.2,0.25,0.3,0.35,0.41,0.47,0.54,0.61,0.68,0.75,0.83,0.92,1.02,1.15,1.28,1.41,2.5});
    else if (binning_option == 1) varYtt->Bins.push_back({0.0,0.07,0.15,0.24,0.34,0.46,0.60,0.74,0.91,1.10,2.5});
  }
  else assert(0);

  ZVarPttt* varPttt = new ZVarPttt(mode);
  varPttt->Title = "p_{T}(t#bar{t})";
  varPttt->TitleTxt = "pT(ttbar)";
  varPttt->Units = "GeV";
  varPttt->NDigits = 0;
  varPttt->NDiv = 303;
  varPttt->Bins.push_back({0, 30, 75, 150, 500});
  if(binning_option == 0) varPttt->Bins.push_back({0,7,10,13,16,19,22,25,28,32,36,40,45,50,55,62,70,79,89,100,112,129,150,195,500});
  else if(binning_option == 1) varPttt->Bins.push_back({0,10,16,22,28,36,45,55,70,89,112,150,500});

  TString xsec_name = "yttpttt";
  if(binning_option!=0) xsec_name = xsec_name + "-binOpt" + std::to_string(binning_option);

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varYtt, varPttt}, configHelper, xsec_name);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  // unfoldingIOHandler->VMGHisto = { std::pair<std::vector<int>, int>({61}, 0),
  //                    std::pair<std::vector<int>, int>({62}, 0),
  //                    std::pair<std::vector<int>, int>({63}, 0),
  //                    std::pair<std::vector<int>, int>({64}, 0) };
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecPttPtttXSec8(PlotterConfigurationHelper *configHelper, int binning_option = 0)
{
  ZVar* varPtt = new ZVarPtt;
  varPtt->Title = "p_{T}(t)";
  varPtt->TitleTxt = "pT(t)";
  varPtt->Expression = "ptt";
  varPtt->Units = "GeV";
  varPtt->NDigits = 0;
  varPtt->Bins.push_back({0, 80, 150, 250, 600});
  if (binning_option == 0) varPtt->Bins.push_back({0,22,36,51,64,72,80,93,104,115,125,133,141,150,170,190,220,250,290,340,600});
  else if (binning_option == 1) varPtt->Bins.push_back({0,22,36,51,64,72,80,93,104,115,125,133,141,150,170,190,220,250,290,340,600});

  ZVarPttt* varPttt = new ZVarPttt(0);
  varPttt->Title = "p_{T}(t#bar{t})";
  varPttt->TitleTxt = "pT(ttbar)";
  varPttt->Expression = "pttt";
  varPttt->Units = "GeV";
  varPttt->NDigits = 0;
  varPttt->NDiv = 303;
  varPttt->Bins.push_back({0, 30, 75, 150, 500});
  if(binning_option == 0) varPttt->Bins.push_back({0,7,10,13,16,19,22,25,28,32,36,40,45,50,55,62,70,79,89,100,112,129,150,195,500});
  else if(binning_option == 1) varPttt->Bins.push_back({0,10,16,22,28,36,45,55,70,89,112,150,500});

  TString xsec_name = "pttpttt";
  if(binning_option!=0) xsec_name = xsec_name + "-binOpt" + std::to_string(binning_option);

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varPtt, varPttt}, configHelper, xsec_name);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecMttPttXSec8(PlotterConfigurationHelper *configHelper, int binning_option = 0)
{
  ZVar* varMtt = new ZVarMtt(0);
  varMtt->Title = "M(t#bar{t})";
  varMtt->TitleTxt = "M(ttbar)";
  varMtt->Expression = "mtt";
  varMtt->Units = "GeV";
  varMtt->NDigits = 0;
  //varMtt->Bins.push_back({340, 440, 550, 750, 1500});
  //varMtt->Bins.push_back({340, 400, 500, 650, 1500});
  //varMtt->Bins.push_back({340, 450, 600, 1500});
  //varMtt->Bins.push_back({340,373,386,400,420,442,467,500,540,592,650,740,1500});
  varMtt->Bins.push_back({300, 450, 600, 1500});
  //varMtt->Bins.push_back({300,373,386,400,420,442,467,500,540,592,650,740,1500});
  if (binning_option == 0) varMtt->Bins.push_back({300,382,400,420,442,467,500,540,592,650,740,1500});
  else if(binning_option == 1) varMtt->Bins.push_back({300,400,442,500,592,740,1500});

  ZVarPtt* varPtt = new ZVarPtt;
  varPtt->Title = "p_{T}(t)";
  varPtt->TitleTxt = "pT(t)";
  varPtt->Expression = "ptt";
  varPtt->Units = "GeV";
  varPtt->NDigits = 0;
  varPtt->NDiv = 303;
  //varPtt->Bins.push_back({0, 80, 150, 220, 600});
  varPtt->Bins.push_back({0, 100, 180, 600});
  if(binning_option == 0) varPtt->Bins.push_back({0,36,64,80,104,125,141,165,210,600});
  else if(binning_option == 1) varPtt->Bins.push_back({0,80,125,141,165,210,600});

  TString xsec_name = "mttptt-xsec8";
  if(binning_option!=0) xsec_name = xsec_name + "-binOpt" + std::to_string(binning_option);

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varMtt, varPtt}, configHelper, xsec_name);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecPtej1PttXSec8(PlotterConfigurationHelper *configHelper, int binning_option = 0)
{
  ZVar* Ptej1 = new ZVarPtDaug("ej1");
  Ptej1->Title = "p_{T}^{lead}(extra jet)";
  Ptej1->TitleTxt = "pT(leading extra jet)";
  Ptej1->Units = "GeV";
  Ptej1->NDigits = 0;
  Ptej1->Bins.push_back({0, 100, 180, 600});
  Ptej1->Bins.push_back({0,80,125,141,165,210,600});

  ZVarPtt* varPtt = new ZVarPtt;
  varPtt->Title = "p_{T}(t)";
  varPtt->TitleTxt = "pT(t)";
  varPtt->Expression = "ptt";
  varPtt->Units = "GeV";
  varPtt->NDigits = 0;
  varPtt->NDiv = 303;
  varPtt->Bins.push_back({0, 100, 180, 600});
  if(binning_option == 0) varPtt->Bins.push_back({0,36,64,80,104,125,141,165,210,600});
  else if(binning_option == 1) varPtt->Bins.push_back({0,80,125,141,165,210,600});

  TString xsec_name = "ptej1ptt";
  if(binning_option!=0) xsec_name = xsec_name + "-binOpt" + std::to_string(binning_option);

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({Ptej1, varPtt}, configHelper, xsec_name);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecMttDphittXSec8(PlotterConfigurationHelper *configHelper, int binning_option = 0)
{
  ZVar* varMtt = DecalteVarMttXSec8(0, binning_option);

  ZVarDphitt* varDphitt = new ZVarDphitt;
  varDphitt->Title = "#Delta#phi(t,#bar{t})";
  varDphitt->TitleTxt = "delta_phi(t,tbar)";
  varDphitt->Expression = "dphitt";
  varDphitt->Units = "rad";
  varDphitt->NDigits = 2;
  varDphitt->Bins.push_back({0, 2.2, 2.95, 3.1415});
  varDphitt->Bins.push_back({0,2.2,2.5,2.67,2.8,2.85,2.9,2.95,3.0,3.05,3.09,3.1415});

  TString xsec_name = "mttdphitt-xsec8";
  if(binning_option!=0) xsec_name = xsec_name + "-binOpt" + std::to_string(binning_option);

  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varMtt, varDphitt}, configHelper, xsec_name);
  unfoldingIOHandler->DoYields = true;
  unfoldingIOHandler->DoPSE = true;
  unfoldingIOHandler->DoXSec = true;
  unfoldingIOHandler->VMGHisto = { std::pair<std::vector<int>, int>({51}, 0),
                     std::pair<std::vector<int>, int>({52}, 0),
                     std::pair<std::vector<int>, int>({53}, 0),
                     std::pair<std::vector<int>, int>({54}, 0) };
  unfoldingIOHandler->Init();
  return unfoldingIOHandler;
}

UnfoldingIOHandler* DeclareXSecNjRhoj(PlotterConfigurationHelper *configHelper, const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const TString& leadStr = "2")
{
  ZVar* varNj = new ZVarNj(def, iso, pt);
  varNj->Bins.push_back({});
  int nbins = leadStr.Atoi() + 1;
  assert(nbins == 2 || nbins == 3);
  for(int b = 0; b < nbins; b++)
    varNj->Bins.back().push_back(-0.5 + b);
  varNj->Bins.back().push_back(8.5);
  varNj->Bins.push_back({-0.5, 0.5, 1.5, 8.5});

  ZVar* varRhoj = new ZVarRhoj(def, iso, pt, leadStr);
  varRhoj->Title += ";N_{j} > 0";
  varRhoj->Bins.push_back({0.0, 0.4, 0.5, 0.6, 1.0});
  varRhoj->Bins.push_back({0.0,0.35,0.42,0.46,0.5,0.54,0.58,0.62,0.67,0.73,1.0});

  TString name = utils::ttmd::GetNameExtraJetVar("njrhoj", def, iso, pt, leadStr);
  UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varNj, varRhoj}, configHelper, name);
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
    //throw std::logic_error("Error: not implemented");
    unfoldingIOHandler->FlagMGHistoIntegrate = true;
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


// 2D vs. extra jets

UnfoldingIOHandler* DeclareXSecNjDetatt(PlotterConfigurationHelper *configHelper,const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const int njNBins = 3){

    ZVarNj* varNj = DeclareVarNj(njNBins, def, iso, pt);

    ZVarDetatt* varDetatt = new ZVarDetatt;
    varDetatt->Title = "#Delta#eta(t,#bar{t})";
    varDetatt->TitleTxt = "delta_eta(t,tbar)";
    varDetatt->Expression = "detatt";
    varDetatt->Units = "";
    varDetatt->NDigits = 1;
    varDetatt->Bins.push_back({0.0, 0.4, 1.2, 6.0});
    // varDetatt->Bins.push_back({0.0,0.12,0.25,0.4,0.55,0.70,0.85,1.0,1.2,1.4,1.65,2.0,2.5,6.0});
    varDetatt->Bins.push_back({0.0,0.25,0.55,0.85,1.2,1.65,2.5,6.0});

    TString name = utils::ttmd::GetNameExtraJetVar("njdeta", def, iso, pt);
    UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varNj, varDetatt}, configHelper, name + TString::Format("-b%d", njNBins));

    unfoldingIOHandler->DoYields = true;
    unfoldingIOHandler->DoPSE = true;
    unfoldingIOHandler->DoXSec = true;
    unfoldingIOHandler->NPertOrders = njNBins;

    unfoldingIOHandler->Init();
    return unfoldingIOHandler;

}

UnfoldingIOHandler* DeclareXSecNjPtt(PlotterConfigurationHelper *configHelper,const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const int njNBins = 3){

    ZVarNj* varNj = DeclareVarNj(njNBins, def, iso, pt);

    ZVarPtt* varPtt = new ZVarPtt;
    varPtt->Title = "p_{T}(t)";
    varPtt->TitleTxt = "pT(t)";
    varPtt->Expression = "ptt";
    varPtt->Units = "GeV";
    varPtt->NDigits = 0;
    varPtt->NDiv = 303;
    //varPtt->Bins.push_back({0, 80, 150, 220, 600});
    varPtt->Bins.push_back({0, 100, 180, 600});
    //varPtt->Bins.push_back({0,36,64,80,104,125,141,165,210,600});
    varPtt->Bins.push_back({0,80,125,141,165,210,600});

    TString name = utils::ttmd::GetNameExtraJetVar("njptt", def, iso, pt);
    UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varNj, varPtt}, configHelper, name + TString::Format("-b%d", njNBins));

    unfoldingIOHandler->DoYields = true;
    unfoldingIOHandler->DoPSE = true;
    unfoldingIOHandler->DoXSec = true;
    unfoldingIOHandler->NPertOrders = njNBins;

    unfoldingIOHandler->Init();
    return unfoldingIOHandler;

}

UnfoldingIOHandler* DeclareXSecNjYt(PlotterConfigurationHelper *configHelper,const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const int njNBins = 3){

  ZVarNj* varNj = DeclareVarNj(njNBins, def, iso, pt);

  ZVarYt* varYt = new ZVarYt;
  varYt->Title = "|y(t)|";
  varYt->TitleTxt = "y(t)";
  varYt->Expression = "yt";
  varYt->Units = "";
  varYt->NDigits = 2;
  varYt->NDiv = 303;
  varYt->Bins.push_back({0.0, 0.35, 0.85, 1.45, 2.5});
  varYt->Bins.push_back({0.0,0.05,0.10,0.16,0.22,0.28,0.35,0.42,0.49,0.57,0.66,0.75,0.85,0.94,1.03,1.12,1.22,1.33,1.45,1.7,2.5}); //TOP-18-004
  // varYt->Bins.push_back({0.0,0.10,0.16,0.28,0.42,0.57,0.75,0.94,1.12,1.33,1.7,2.5});

    TString name = utils::ttmd::GetNameExtraJetVar("njyt", def, iso, pt);
    UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varNj, varYt}, configHelper, name + TString::Format("-b%d", njNBins));

    unfoldingIOHandler->DoYields = true;
    unfoldingIOHandler->DoPSE = true;
    unfoldingIOHandler->DoXSec = true;
    unfoldingIOHandler->NPertOrders = njNBins;

    unfoldingIOHandler->Init();
    return unfoldingIOHandler;

}

UnfoldingIOHandler* DeclareXSecNjYtt(PlotterConfigurationHelper *configHelper,const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const int njNBins = 3, const int mode = 0){

    ZVarNj* varNj = DeclareVarNj(njNBins, def, iso, pt);

    ZVarYtt* varYtt = new ZVarYtt(mode);
    varYtt->NDiv = 303;
    varYtt->Bins.push_back({0.0, 0.35, 0.75, 1.15, 2.5});
    varYtt->Bins.push_back({0.0,0.04,0.08,0.12,0.16,0.2,0.25,0.3,0.35,0.41,0.47,0.54,0.61,0.68,0.75,0.83,0.92,1.02,1.15,1.28,1.41,2.5});

    TString name = utils::ttmd::GetNameExtraJetVar("njytt", def, iso, pt);
    UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varNj, varYtt}, configHelper, name + TString::Format("-b%d", njNBins));

    unfoldingIOHandler->DoYields = true;
    unfoldingIOHandler->DoPSE = true;
    unfoldingIOHandler->DoXSec = true;
    unfoldingIOHandler->NPertOrders = njNBins;

    unfoldingIOHandler->Init();
    return unfoldingIOHandler;

}

UnfoldingIOHandler* DeclareXSecNjPttt(PlotterConfigurationHelper *configHelper,const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const int njNBins = 3, const int mode = 0){

    ZVarNj* varNj = DeclareVarNj(njNBins, def, iso, pt);

    ZVarPttt* varPttt = new ZVarPttt(mode);
    varPttt->Title = "p_{T}(t#bar{t})";
    varPttt->TitleTxt = "pT(ttbar)";
    varPttt->Units = "GeV";
    varPttt->NDigits = 0;
    varPttt->NDiv = 303;
    varPttt->Bins.push_back({0, 30, 75, 130, 500});
    //varPttt->Bins.push_back({0,7,10,13,16,19,22,25,28,32,36,40,45,50,55,62,70,79,89,100,112,129,150,195,500});
    varPttt->Bins.push_back({0,10,16,22,28,36,45,55,70,89,112,150,500});

    TString name = utils::ttmd::GetNameExtraJetVar("njpttt", def, iso, pt);
    UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varNj, varPttt}, configHelper, name + TString::Format("-b%d", njNBins));

    unfoldingIOHandler->DoYields = true;
    unfoldingIOHandler->DoPSE = true;
    unfoldingIOHandler->DoXSec = true;
    unfoldingIOHandler->NPertOrders = njNBins;

    unfoldingIOHandler->Init();
    return unfoldingIOHandler;

}

UnfoldingIOHandler* DeclareXSecNjMtt(PlotterConfigurationHelper *configHelper,const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const int njNBins = 3, const int mode = 0){

    ZVarNj* varNj = DeclareVarNj(njNBins, def, iso, pt);

    ZVar* varMtt = new ZVarMtt(mode);
    varMtt->Title = "M(t#bar{t})";
    varMtt->TitleTxt = "M(ttbar)";
    varMtt->Units = "GeV";
    varMtt->NDigits = 0;

    //varMtt->Bins.push_back({300, 400, 500, 650, 1500});
    //varMtt->Bins.push_back({300,373,386,400,420,442,467,500,540,592,650,740,1500});

    std::pair<std::vector<double>,std::vector<double>> mtt_binning_info = UnfoldingIOHandler::GetBinningInfoFromFile("data/binnings/multiDiffXSecs/mtt.txt");
    varMtt->Bins.push_back(mtt_binning_info.first);
    varMtt->Bins.push_back(mtt_binning_info.second);

    TString name = utils::ttmd::GetNameExtraJetVar("njmtt", def, iso, pt);
    UnfoldingIOHandler* unfoldingIOHandler = new UnfoldingIOHandler({varNj, varMtt}, configHelper, name + TString::Format("-b%d", njNBins));

    unfoldingIOHandler->DoYields = true;
    unfoldingIOHandler->DoPSE = true;
    unfoldingIOHandler->DoXSec = true;
    unfoldingIOHandler->NPertOrders = njNBins;

    unfoldingIOHandler->Init();
    return unfoldingIOHandler;

}


// ******************************************************************************************
// ************************* summary methods ************************************************
// ******************************************************************************************

std::vector<UnfoldingIOHandler*> DeclareXSecAll82D(PlotterConfigurationHelper *configHelper)
{
  std::vector<UnfoldingIOHandler*> vec;
  vec.push_back(DeclareXSecYtPttXSec8(configHelper));
  vec.push_back(DeclareXSecMttYtXSec8(configHelper));
  vec.push_back(DeclareXSecMttYttXSec8(configHelper));
  vec.push_back(DeclareXSecMttDetattXSec8(configHelper));
  vec.push_back(DeclareXSecMttPtttXSec8(configHelper));
  vec.push_back(DeclareXSecMttDphittXSec8(configHelper));
  return vec;
}
#endif
