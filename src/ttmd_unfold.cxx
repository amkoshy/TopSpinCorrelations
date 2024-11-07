#include "ttmd_unfold.h"

#include "UnfoldingXSecHandler.h"
#include "utils.h"
#include "PlotterConfigurationHelper.h"
#include "UnfoldingIOHandler.h"
#include "PlotterDiffXSec.h"
#include "TUnfoldDensity.h"

#include <TLegend.h>

ZUnfold::ZUnfold(PlotterConfigurationHelper *configHelper):
configHelper_(configHelper)
{
  tunf = NULL;
  regChoice = RegRhoAvg;
  tau = -1.0;

  flagUFRec = 0;
  flagOFRec = 0;
  flagUFGen = 0;
  flagOFGen = 0;

  constraintMode= TUnfold::kEConstraintArea;
  //constraintMode= TUnfold::kEConstraintNone;

  regMode = TUnfold::kRegModeCurvature;
  //regMode = TUnfold::kRegModeDerivative;
  //regMode = TUnfold::kRegModeSize;
  //regMode = TUnfold::kRegModeNone;

  //densityFlags = TUnfoldDensity::kDensityModeBinWidth;
  densityFlags = TUnfoldDensity::kDensityModeNone;

  // detailed steering for regularisation
  REGDISTR = NULL;
  REGAXISSTEER = NULL;

  // temporary
  hEffBBB = NULL;
}

int ZUnfold::GetDataChi2Dim() const
{
  if(Remark == "BBB")
    return 0;
  else if(Remark == "NoReg")
    return 1;
  else if(Remark == "MinRhoAvg")
    return 2;
  else
  {
    printf("Error in ZUnfold::GetDataChi2Dim(): Remark = %s not implemented\n", Remark.Data());
    throw;
  }
}

void ZUnfold::RunTUnfold(const UnfoldingIOHandler* unfoldingIOHandler, const TString& fileNameBase/* = ""*/)
{
  _unfoldingIOHandler = unfoldingIOHandler;
  bool flagSameBins = false;
  if(_unfoldingIOHandler->NBinnings() == 1)
  {
    printf("Warning in RunTUnfold(): no fine detector binning, using coarse binning instead\n");
    flagSameBins = true;
  }

  TUnfoldBinning* genBinning = unfoldingIOHandler->GetTUnfoldBinningC();
  TUnfoldBinning* detBinning = flagSameBins ? genBinning : unfoldingIOHandler->GetTUnfoldBinningF();
  //detBinning->PrintStream(std::cout);
  //genBinning->PrintStream(std::cout);
  if(regChoice == RegNo && tau < 0)
    tau = 0.0;

  //TH2D* hRsp = (TH2D*) (flagSameBins ? unfoldingIOHandler->HRspC() : unfoldingIOHandler->HRspTUnfold())->Clone();
  hRsp = new TH2D(*(flagSameBins ? unfoldingIOHandler->HRspC() : unfoldingIOHandler->HRspTUnfold()));
  //hRsp->Print("all");
  //if(hInput)
  //  delete hInput;
  hInput = (TH1D*) (flagSameBins ? unfoldingIOHandler->HDnbC() : unfoldingIOHandler->HDnbF())->Clone();
  CheckStatisticsIsGaussian(hInput);
  //if(!tunf)
  //printf("hInput size: %d\n", hInput->GetNbinsX());
  //printf("hRsp size: %d %d\n", hRsp->GetNbinsX(), hRsp->GetNbinsY());
  tunf = new TUnfoldDensity(hRsp, TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags, genBinning, detBinning, REGDISTR, REGAXISSTEER);
  //if(tunf->SetInput(hInput) >= 10000)
  if(tunf->SetInput(hInput, 1.0) >= 10000)
  {
    printf("Unfolding result may be wrong\n");
    throw "SetInput >= 10000";
  }

  if(regChoice == RegNo)
    tunf->DoUnfold(tau);
  else if(regChoice == RegRhoAvg || regChoice == RegRhoMax)
  {
    int nScan = 100;
    TSpline *rhoLogTau = NULL;
    double minTau = 1e-6;
    double maxTau = 1e-0;
    TGraph *lCurve = NULL;
    TSpline *lCurveX = NULL;
    TSpline *lCurveY = NULL;
    const char *SCAN_DISTRIBUTION = NULL;
    const char *SCAN_AXISSTEERING = NULL;
    TUnfoldDensity::EScanTauMode scanTauMode;
    if(regChoice == RegRhoAvg)
      scanTauMode = TUnfoldDensity::kEScanTauRhoAvg;
    else if(regChoice == RegRhoMax)
      scanTauMode = TUnfoldDensity::kEScanTauRhoMax;
    int iBest = tunf->ScanTau(nScan, minTau, maxTau, &rhoLogTau, scanTauMode, SCAN_DISTRIBUTION, SCAN_AXISSTEERING, &lCurve, &lCurveX, &lCurveY);

    // create graphs with one point to visualize best choice of tau
    double t[1], rho[1], x[1], y[1];
    rhoLogTau->GetKnot(iBest, t[0], rho[0]);
    lCurve->GetPoint(iBest, x[0], y[0]);
    TGraph* bestRhoLogTau = new TGraph(1, t, rho);
    //TGraph* bestLCurve = new TGraph(1, x, y);
    double* tAll = new double[nScan];
    double* rhoAll = new double[nScan];
    for(int i = 0; i < nScan; i++)
      rhoLogTau->GetKnot(i, tAll[i], rhoAll[i]);
    //TGraph* knots = new TGraph(nScan, tAll, rhoAll);
    tau = TMath::Power(10.0, t[0]);
    printf("best tau = %.5e\n", tau);
    //tau = 1000000.0;
    tunf->DoUnfold(tau);

    // plot graphs
    if(fileNameBase != "")
    {
      TCanvas canvas("", "", configHelper_->gPlotPixelSmall,configHelper_-> gPlotPixelSmall);
      //canvas.Divide(1, 1);
      //canvas.cd();
      // (4) scan of correlation vs tau
      canvas.cd();
      TH2D* hr = new TH2D("", "", 1, TMath::Log10(minTau), TMath::Log10(maxTau), 1, 0, 1);
      hr->GetXaxis()->SetTitle("lg(#tau)");
      if(regChoice == RegRhoAvg)
        hr->GetYaxis()->SetTitle("Average global cor. coef.");
      else if(regChoice == RegRhoMax)
        hr->GetYaxis()->SetTitle("Maximum global cor. coef.");
      hr->GetXaxis()->SetTitleOffset(1.10);
      hr->GetYaxis()->SetTitleOffset(1.25);
      hr->Draw();
      //knots->Draw("*");
      bestRhoLogTau->SetMarkerColor(kRed);
      bestRhoLogTau->Draw("same *");
      rhoLogTau->Draw("same");
      // (5) global correlation coefficients for the distributions
      //   used during the scan
      /*canvas.cd(2);
      //histCorrCoeff->Draw("BOX");
      TH1 *histGlobalCorrScan=unfold.GetRhoIJtotal
       ("histGlobalCorrScan");
      histGlobalCorrScan->Draw("HIST");
      // (6) L-curve
      canvas.cd(3);
      lCurve=0;
      iBest=unfold.ScanLcurve(nScan,1e-6,1e-1,&lCurve,&lCurveX,&lCurveY);
      lCurve->GetPoint(iBest,x[0],y[0]);
      lCurveX->GetKnot(iBest,x[0],y[0]);
      printf("best tau (l-curve) = %f\n", TMath::Power(10.0, x[0]));
      bestLCurve=new TGraph(1,x,y);
      lCurve->Draw("AL");
      bestLCurve->SetMarkerColor(kRed);
      bestLCurve->Draw("*");*/
      utils::ttmd::SaveCanvas(&canvas, fileNameBase + "scanrho");

      // right tau to file
      utils::ttmd::CheckFile(fileNameBase + "tau.txt");
      FILE* fout = fopen(fileNameBase + "tau.txt", "w");
      fprintf(fout, "%.3e", tau);
      TauBest = tau;
      fclose(fout);
    }
  }

  // get results
  printf("chi2 = %.1f + %.1f / %d\n", tunf->GetChi2A(), tunf->GetChi2L(), tunf->GetNdf());
  Chi2Data = tunf->GetChi2A() + tunf->GetChi2L();
  Chi2Reg = tunf->GetChi2L();
  Dof = tunf->GetNdf();
  hUnf = (TH1D*) tunf->GetOutput("", NULL, NULL, NULL, false);
  //printf("hUnf size: %d\n", hUnf->GetNbinsX());
  if(configHelper_->gUnfUncRspTotal == 1)
    hEmatTotal = (TH2D*) tunf->GetEmatrixTotal("");
  else if(configHelper_->gUnfUncRspTotal == 0)
    hEmatTotal = (TH2D*) tunf->GetEmatrixInput("");
  else
    assert(0);
  for(int b = 0; b < hUnf->GetNbinsX(); b++)
    hUnf->SetBinError(b + 1, TMath::Sqrt(hEmatTotal->GetBinContent(b + 1, b + 1)));
}

void ZUnfold::MakeNormalisation()
{
  // calculate normalised cross section
  //if(hUnfN)
  //  delete hUnfN;
  hUnfN = (TH1D*) hUnf->Clone();
  //if(hEmatTotalN)
  //  delete hEmatTotalN;
  hEmatTotalN = (TH2D*) hEmatTotal->Clone();
  utils::ttmd::Normalise(hUnfN, hEmatTotalN);
  //if(hCormatTotalN)
  //  delete hCormatTotalN;
  hCormatTotalN = utils::ttmd::MakeCorrMatrix(hEmatTotalN);
  //if(hGenN)
  //  delete hGenN;
  hGenN = (TH1D*) hGen->Clone();
  hGenN->Scale(1.0 / hGen->Integral());
}

void ZUnfold::CalculateChi2()
{
  this->Chi2 = utils::ttmd::Chi2(hUnf, hEmatTotal, hGen);
  this->Chi2N = utils::ttmd::Chi2(hUnfN, hEmatTotalN, hGenN, 1);
}

TH1D* ZUnfold::GetEfficiency(UnfoldingIOHandler* unfoldingIOHandler, const int flag/* = 1*/)
{
  TH1D* hEff = (TH1D*) unfoldingIOHandler->HRecC()->Clone();

  if(flag == 1)
    hEff->Divide(unfoldingIOHandler->HGenC());
  else if(flag == 0)
  {
    hEff->Reset();
    const TH2D* hRsp = unfoldingIOHandler->HRspC();
    int nbins = hEff->GetNbinsX();
    for(int i = 0; i < nbins; i++)
    {
      double val = 0.0;
      for(int j = 0; j < nbins; j++)
        val += hRsp->GetBinContent(j + 1, i + 1);
      hEff->SetBinContent(i + 1, val);
    }
    hEff->Divide(unfoldingIOHandler->HGenC());
  }
  else
    throw std::logic_error(TString::Format("Error in GetEfficiency(): unknown flag = %d", flag));

  for(int b = 0; b < hEff->GetNbinsX(); b++)
    hEff->SetBinError(b + 1, 0.0);

  return hEff;
}

void ZUnfold::PlotInput(const TString& fileName)
{
  if(!hInput || !tunf)
  {
    printf("Error in ZUnfold::PlotInput() no input hInput = %p, tunf = %p, call RunTUnfold()\n", (void *)hInput, (void *)tunf);
    throw;
  }

  double ratMin = 0.6;
  double ratMax = 1.4;
  double scaleAxisFontSize = 5.0;
  if(_unfoldingIOHandler->Dim() == 1)
    scaleAxisFontSize = 1.35;
  double eps = 1e-3;
  int div = 304;

  // folded back
  TH1D *hFold = (TH1D*) tunf->GetFoldedOutput("", NULL, NULL, NULL, false);
  utils::ttmd::SetErrorsToZero(hFold);

  // prepare histograms
  PlotterDiffXSec plotter;
  const ZVar* var = _unfoldingIOHandler->Var(_unfoldingIOHandler->Dim() - 1);
  const int binning = (_unfoldingIOHandler->NBinnings() == 1) ? 0 : 1;
  std::vector<TH1D*> hhData = _unfoldingIOHandler->CreateLastDimHistos(binning);
  std::vector<TH1D*> hhFold = _unfoldingIOHandler->CreateLastDimHistos(binning);
  std::vector<TH1D*> hhFoldRatio;
  std::vector<TH1D*> hhDataRatio;
  std::vector<TH2D*> hRange;
  std::vector<TH2D*> hRangeRatio;
  int npads = hhData.size();
  double max = 0.0;
  for(int b = 1; b <= hInput->GetNbinsX(); b++)
  {
    max = std::max(max, hInput->GetBinContent(b));
    max = std::max(max, hFold->GetBinContent(b));
  }
  max *= 1.1;

  // plot
  TCanvas c("", "", configHelper_->gPlotPixelLarge, configHelper_->gPlotPixelLarge);
  TCanvas cR("", "", configHelper_->gPlotPixelLarge, configHelper_->gPlotPixelLarge);
  if(_unfoldingIOHandler->Dim() == 1)
  {
    c.SetMargin(0.15, 0.05, 0.15, 0.05);
    cR.SetMargin(0.15, 0.05, 0.15, 0.05);
  }
  else
  {
    c.SetMargin(0.30, 0.10, 0.30, 0.10);
    cR.SetMargin(0.30, 0.10, 0.30, 0.10);
  }
  int w = -1;
  int h = -1;
  if(_unfoldingIOHandler->Dim() == 1)
  {
    w = 1;
    h = 1;
    c.Divide(1);
    cR.Divide(1);
  }
  else
  {
    utils::ttmd::DivideSquareWH(&c, npads + 1, w, h);
    c.DivideSquare(npads + 1, 0.005, 0.00);
    cR.DivideSquare(npads + 1, 0.005, 0.00);
    //c.DivideSquare(npads + 1);
  }
  TObject* objDataAsGraph = NULL;

  TH2D hr("", "", 1, hhFold[0]->GetBinLowEdge(1) + eps, hhFold[0]->GetBinLowEdge(hhFold[0]->GetNbinsX() + 1) - eps, 1, eps, max - eps);
  hr.GetXaxis()->SetTitle(var->GetXTitle());
  hr.GetYaxis()->SetTitle("Events");
  hr.GetYaxis()->SetTitleOffset(1.1);
  utils::ttmd::ScaleHistoFonts(&hr, scaleAxisFontSize);
  if(_unfoldingIOHandler->Dim() != 1)
  {
    hr.GetXaxis()->SetNdivisions(div);
    hr.GetYaxis()->SetNdivisions(div);
  }

  TH2D hrR("", "", 1, hhFold[0]->GetBinLowEdge(1) + eps, hhFold[0]->GetBinLowEdge(hhFold[0]->GetNbinsX() + 1) - eps, 1, ratMin + eps, ratMax - eps);
  hrR.GetXaxis()->SetTitle(var->GetXTitle());
  hrR.GetYaxis()->SetTitle("Ratio");
  hrR.GetXaxis()->SetNdivisions(div);
  hrR.GetYaxis()->SetNdivisions(div);
  hrR.GetYaxis()->SetTitleOffset(1.1);
  utils::ttmd::ScaleHistoFonts(&hrR, scaleAxisFontSize);

  for(int p = npads - 1; p >= 0; p--)
  {
    _unfoldingIOHandler->FillLastDimHisto(hhData[p], hInput);
    _unfoldingIOHandler->FillLastDimHisto(hhFold[p], hFold);
    utils::ttmd::SetErrorsToZero(hhFold[p]);
    c.cd(p + 1);
    hRange.push_back((TH2D*) hr.Clone());
    utils::ttmd::AdjustDividedCanvasHistograms(hRange.back(), p, npads, w/*, h*/);
    hRange.back()->Draw();
    hhFold[p]->SetLineColor(4);
    hhFold[p]->Draw("hist same");
    hhData[p]->SetMarkerSize(0.5);
    hhData[p]->SetMarkerStyle(20.0);
    hhData[p]->SetMarkerColor(1);
    hhData[p]->SetLineColor(1);
    objDataAsGraph = utils::ttmd::DrawAsGraph(hhData[p]);
    hhFold[p]->Draw("hist same");
    hRange.back()->Draw("axis same");

    // ratio
    hRangeRatio.push_back((TH2D*) hrR.Clone());
    hhDataRatio.push_back((TH1D*)hhData[p]->Clone());
    hhFoldRatio.push_back((TH1D*)hhFold[p]->Clone());
    hhDataRatio.back()->Divide(hhFold[p]);
    hhFoldRatio.back()->Divide(hhFold[p]);
    cR.cd(p + 1);
    utils::ttmd::AdjustDividedCanvasHistograms(hRangeRatio.back(), p, npads, w/*, h*/);
    hRangeRatio.back()->Draw();
    hhFoldRatio.back()->Draw("hist same");
    utils::ttmd::DrawAsGraph(hhDataRatio.back());
    hRangeRatio.back()->Draw("axis same");
  }

  // legend
  TLegend* leg = NULL;
  if(_unfoldingIOHandler->Dim() == 1)
    leg = new TLegend(0.65, 0.65, 0.89, 0.89, "", "NDC");
  else
    leg = new TLegend(0.1, 0.35, 0.9, 0.95, "", "NDC");
  leg->SetTextFont(62);
  leg->SetTextSize(1.5 * leg->GetTextSize());
  leg->SetFillStyle(0);
  leg->AddEntry(objDataAsGraph, "Data", "pe");
  leg->AddEntry(hhFold[0], "Fold", "l");
  if(_unfoldingIOHandler->Dim() == 1)
    c.cd(1);
  else
    c.cd(npads + 1);
  leg->Draw();
  if(_unfoldingIOHandler->Dim() == 1)
    cR.cd(1);
  else
    cR.cd(npads + 1);
  leg->Draw();

  utils::ttmd::SaveCanvas(&c, fileName);
  utils::ttmd::SaveCanvas(&cR, fileName + "-rat");

  for(int p = 0; p < npads; p++)
  {
    delete hhData[p];
    delete hhFold[p];
    delete hRange[p];
    delete hhDataRatio[p];
    delete hhFoldRatio[p];
    delete hRangeRatio[p];
  }
}

void ZUnfold::RunBBB(UnfoldingIOHandler* unfoldingIOHandler, const int flagEff/* = 1*/)
{
  TH1D* hEff = this->GetEfficiency(unfoldingIOHandler, flagEff);
  //if(hUnf)
  //  delete hUnf;
  hUnf = (TH1D*) unfoldingIOHandler->HDnbC()->Clone();
  hUnf->Divide(hEff);
  hEmatTotal = utils::ttmd::MakeCovTH2(hUnf);
  hEffBBB = hEff;
}

void ZUnfold::RunUnfolding(UnfoldingIOHandler* unfoldingIOHandler, const TString& ch, const TString& fileName/* = ""*/)
{
  //if(hGen)
  //  delete hGen;
  hGen = (TH1D*) unfoldingIOHandler->HGenC()->Clone();

  if(Remark == "BBB")
    RunBBB(unfoldingIOHandler, 1);
  else if(Remark == "BBB0")
    RunBBB(unfoldingIOHandler, 0);
  else
    RunTUnfold(unfoldingIOHandler, fileName);

  // absolute scaling
  if (!unfoldingIOHandler->DoParticle) {
    double br = configHelper_->GetChannelBR(ch);
    hGen->Scale(1.0 / configHelper_->lumi / br);
    this->hUnf->Scale(1.0 / configHelper_->lumi / br);
    hEmatTotal->Scale(TMath::Power(1.0 / configHelper_->lumi / br, 2.0));
  }
  else if (unfoldingIOHandler->DoParticle) {
    hGen->Scale(1.0 / configHelper_->lumi);
    this->hUnf->Scale(1.0 / configHelper_->lumi);
    hEmatTotal->Scale(TMath::Power(1.0 / configHelper_->lumi, 2.0));
  }

  //if(hCormatTotal)
  //  delete hCormatTotal;
  hCormatTotal = utils::ttmd::MakeCorrMatrix(hEmatTotal);
  SetHistoStyle(this->hUnf);
  MakeNormalisation();
  CalculateChi2();
}

void ZUnfold::CheckStatisticsIsGaussian(const TH1* h)
{
  const double gausMin = 30;
  double minContent = 1e10;
  double minEntries = 1e10;
  for(int b = 1; b < h->GetNbinsX(); b++)
  {
    // check bin content
    double content = h->GetBinContent(b);
    // check number of entries assuming Poisson events with equal weights
    double unc = h->GetBinError(b);
    double entr = TMath::Power(content / unc, 2.0);
    if(content < minContent)
    {
      minContent = content;
      minEntries = entr;
    }
    // issue warning
    if(content < gausMin || content < gausMin)
    {
      // report warning only if corresponding flag is true, unless there are 0 events
      if(FlagTreatStatisticsSeriously || content <= 0.0)
        printf("Warning ZUnfold::CheckStatisticsIsGaussian() bin %d contains %.0f (estimate %.0f entries) [below %.0f] histogram name = %s title = %s\n", b, content, entr, gausMin, h->GetName(), h->GetTitle());
      if(FlagTreatStatisticsSeriously && content <= 0.0)
      {
        h->Print("all");
        throw std::runtime_error("Error: too low statistics. Stop.");
      }
    }
  }
  //if(gTTBARDebug)
  printf("Info ZUnfold::CheckStatisticsIsGaussian() min bin content = %.0f (estimate %.0f entries) histogram name = %s title = %s\n", minContent, minEntries, h->GetName(), h->GetTitle());
}

void ZUnfold::SetHistoStyle(TH1* h)
{
  h->SetTitle(Suffix);
  h->SetLineColor(Color);
  h->SetMarkerColor(Color);
  h->SetMarkerStyle(MarkerStyle);
}
