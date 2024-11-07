#include "UnfoldingClosureTestsHandler.h"

#include "utils.h"
#include "UnfoldingXSecHandler.h"
#include "PlotterDiffXSec.h"
#include "ttmd_treePseudo.h"
#include "PlotterConfigurationHelper.h"
#include "UnfoldingIOHandler.h"
#include "ttmd_unfold.h"

#include <TLegend.h>

void UnfoldingClosureTestsHandler::WriteSample(const TString& rootFileName, std::vector<ZUnfold*> vUnf,
                               UnfoldingIOHandler* unfoldingIOHandler, const TH1D* hGen, const TString& ch, const int seed/* = -1*/)
{
  printf("ZClosureTest::WriteSample() fileName = %s\n", rootFileName.Data());

  //const TString ch = "emu";

  // set seed
  if(seed >= 0)
  {
    printf("ZClosureTest::WriteSample() seed = %d\n", seed);
    ZRandom::GetInstance().SetSeed(seed);
  }

  // open file, setup tree
  utils::ttmd::CheckFile(rootFileName);
  TFile* fout = TFile::Open(rootFileName, "recreate");
  fout->cd();
  TTree* tree = new TTree("treeUnf", "");

  // gen histo
  //TH1D* hGenStore = (TH1D*) vUnf[0]->hGen->Clone();
  //TH1D* hGenStore = (TH1D*) hGen->Clone();
  //hGenStore->Scale(_unfoldingIOHandler->HGenC()->Integral() / hGenStore->Integral());

  // if there is reweight function in provided unfoldingIOHandler, write original gen distribution used to build response matrix
  if(unfoldingIOHandler->GetReweightFunction())
  {
    TH1D* hgenRsp = (TH1D*) _unfoldingIOHandler->HGenC()->Clone();
    //TH1D* hgenRsp = (TH1D*) hGen->Clone();
    //hgenRsp->Scale(hGen->Integral() / hgenRsp->Integral());
    double br = configHelper_->GetChannelBR(ch);
    hgenRsp->Scale(1.0 / configHelper_->lumi / br);
    hgenRsp->SetName("hgenRsp");
    hgenRsp->Write();
  }

  // branch with data chi2
  int nalg = vUnf.size();
  assert(nalg);
  float brChi2[nalg][2];
  int chi2Dim[nalg];
  float brTau[nalg];
  for(int a = 0; a < nalg; a++)
  {
    chi2Dim[a] = vUnf[a]->GetDataChi2Dim();
    if(chi2Dim[a] > 0)
    tree->Branch(TString::Format("chi2%s", vUnf[a]->Remark.Data()), brChi2[a], TString::Format("chi2%s[%d]/F", vUnf[a]->Remark.Data(), chi2Dim[a]));
    if(chi2Dim[a] == 2)
    tree->Branch(TString::Format("tau%s", vUnf[a]->Remark.Data()), &brTau[a], TString::Format("tau%s/F", vUnf[a]->Remark.Data()));
  }

  // branch with histograms
  //int nbins = unfoldingIOHandler->NBinsC();
  const int nbins = unfoldingIOHandler->HGenC()->GetNbinsX();
  int nbinsUndiag = nbins * (nbins - 1) / 2;
  float brHUnf[nalg][nbins];
  float brHCovDiag[nalg][nbins];
  float brHCovUndiag[nalg][nbinsUndiag];
  for(int a = 0; a < nalg; a++)
  {
    tree->Branch(TString::Format("hunf%s", vUnf[a]->Remark.Data()), brHUnf[a], TString::Format("hunf%s[%d]/F", vUnf[a]->Remark.Data(), nbins));
    tree->Branch(TString::Format("hcovdiag%s", vUnf[a]->Remark.Data()), brHCovDiag[a], TString::Format("hcovdiag%s[%d]/F", vUnf[a]->Remark.Data(), nbins));
    if(vUnf[a]->Remark != "BBB")
      tree->Branch(TString::Format("hcovundiag%s", vUnf[a]->Remark.Data()), brHCovUndiag[a], TString::Format("hcovundiag%s[%d]/F", vUnf[a]->Remark.Data(), nbinsUndiag));
  }

  // iterate
  for(int t = 0; t < _niter + 1; t++)
  {
    printf("*** ZClosureTest::WriteSample()  ITERATION = %d  unfoldingIOHandler = %s ***\n", t, _unfoldingIOHandler->Suffix.Data());

    // generate pseudodata
    unfoldingIOHandler->GeneratePseudoData(0, hGen, unfoldingIOHandler->HRspC(), NULL, t);
    unfoldingIOHandler->GeneratePseudoData(1, hGen, unfoldingIOHandler->HRspTUnfold(), NULL, t);
    unfoldingIOHandler->SubtrackBackground(0);

    // unfold pseudodata
    for(int u = 0; u < nalg; u++)
    {
      ZUnfold* unf = vUnf[u];
      unf->RunUnfolding(unfoldingIOHandler, ch);

      // fill data chi2 branch
      if(chi2Dim[u] >= 1)
        brChi2[u][0] = vUnf[u]->tunf->GetChi2A();
      if(chi2Dim[u] == 2)
      {
        brChi2[u][1] = vUnf[u]->tunf->GetChi2L();
        brTau[u] = vUnf[u]->tau;
      }

      // fill histo branch
      // TODO optimise
      utils::ttmd::HCovToFArrays(vUnf[u]->hUnf, vUnf[u]->hEmatTotal, brHUnf[u], brHCovDiag[u], brHCovUndiag[u]);

      if(t == 0 && u == 0)
      {
      }
    }

    // fill tree
    tree->Fill();
  }

  // write gen histo
  TH1D* hGenStore = (TH1D*) hGen->Clone();
  double br = configHelper_->GetChannelBR(ch);
  hGenStore->Scale(1.0 / configHelper_->lumi / br);
  //TH1D* hGenStore = (TH1D*) hGen->Clone();
  //hGenStore->Scale(_unfoldingIOHandler->HGenC()->Integral() / hGenStore->Integral());
  hGenStore->SetName("hgen");
  hGenStore->Write();

  // save ROOT file
  fout->cd();
  tree->Write();

  fout->Close();
}

void UnfoldingClosureTestsHandler::ReadSample(const TString& rootFileName, std::vector<ZUnfold*> vUnf, const TString& dirOut)
{
  //int niter = _niter + 1;
  int nalg = vUnf.size();
  //int nbins = _unfoldingIOHandler->NBinsC();
  // TODO find better solution to get number of bins
  TH1D* hTmp = _unfoldingIOHandler->GetPlainHisto();
  int nbins = hTmp->GetNbinsX();
  delete hTmp;

  // read tree and gen histo
  printf("ReadSample(): rootFileName = %s\n", rootFileName.Data());
  ZTreePseudo tree(vUnf, nbins, rootFileName);
  if((_niter + 1) != tree.GetEntries())
    printf("Warning in ZClosureTest::ClosureTest2() (_niter + 1) = %d != %lld = tree->GetEvents()\n", (_niter + 1), tree.GetEntries());
  TH1D* hGen = tree.HGen();
  TH1D* hGenRsp = tree.HGenRsp();

  // structures with stored results
  std::pair<std::vector<std::vector<TH1D*> >, std::vector<std::vector<TH2D*> > > pairHUnfCov = tree.ReadAllHUnfCov();
  std::vector<std::vector<TH1D*> > vvHUnf = pairHUnfCov.first;
  std::vector<std::vector<TH2D*> > vvHCov = pairHUnfCov.second;

  // averaged covariance
  std::vector<TH2D*> avCov(nalg);
  for(int u = 0; u < nalg; u++)
  {
    avCov[u] = (TH2D*) vvHCov[u][0]->Clone();
    for(int t = 1; t < _niter; t++)
      avCov[u]->Add(vvHCov[u][t]);
    avCov[u]->Scale(1.0 / _niter);
  }

  // unfolded and covariance from sample (^)
  std::pair<TH1D*, TH2D*> avHAndCov[nalg];
  for(int u = 0; u < nalg; u++)
    avHAndCov[u] = UnfoldingClosureTestsHandler::GetAverageAndCov(vvHUnf[u]);

  // compare covariance averaged vs from sample
  for(int u = 0; u < nalg; u++)
    this->Plot2Cov(avCov[u], avHAndCov[u].second, dirOut + "/cov-" + vUnf[u]->Suffix);

  // determine x axis boundaries
  /*std::vector<double> vMinMaxChi2AL = tree.GetChi2ALMinMax();
  std::vector<double> vMinMaxTau = tree.GetTauMinMax();
  std::vector<double> vMinPul(nbins + 1);
  std::vector<double> vMaxPul(nbins + 1);
  std::vector<double> vMinRes(nbins + 1);
  std::vector<double> vMaxRes(nbins + 1);
  tree.GetPulResMinMax(vMinPul, vMaxPul, vMinRes, vMaxRes);*/

  // chi2 ^ vs gen
  std::vector<TH1D*> vHChi2 = UnfoldingClosureTestsHandler::SetupHistogramsChi2(nalg/*, _niter*/, nbins);
  for(int u = 0; u < nalg; u++)
  {
    vUnf[u]->SetHistoStyle(vHChi2[u]);
    for(int t = 1; t < _niter; t++)
    {
      //vvHUnf[u][t]->Print("all");
      //hGen->Print("all");
      double chi2 = utils::ttmd::Chi2(vvHUnf[u][t], vvHCov[u][t], hGen);
      //printf("chi2 = %f\n", chi2);
      vHChi2[u]->Fill(chi2);
    }
    //vHChi2[u]->Print("all");
  }
  UnfoldingClosureTestsHandler::PlotVH("chi", vHChi2, dirOut + "/chi2");
  utils::ttmd::ClearVectorWithContent(vHChi2);

  // residuals and pulls in bins
  std::vector<std::vector<TH1D*> > vvHRes = UnfoldingClosureTestsHandler::SetupHistogramsResiduals(nalg/*, _niter*/, nbins);
  std::vector<std::vector<TH1D*> > vvHPul = UnfoldingClosureTestsHandler::SetupHistogramsPulls(nalg/*, _niter*/, nbins);
  for(int u = 0; u < nalg; u++)
  {
    for(int t = 1; t < _niter; t++)
    {
      for(int b = 0; b < nbins; b++)
      {
        double res = vvHUnf[u][t]->GetBinContent(b + 1) - hGen->GetBinContent(b + 1);
        double resRel = res / hGen->GetBinContent(b + 1);
        vvHRes[u][b]->Fill(resRel);
        double pul = res / TMath::Sqrt(vvHCov[u][t]->GetBinContent(b + 1, b + 1));
        vvHPul[u][b]->Fill(pul);
        if(t == 1)
        {
          vUnf[u]->SetHistoStyle(vvHRes[u][b]);
          vUnf[u]->SetHistoStyle(vvHPul[u][b]);
        }
      }
    }
  }
  UnfoldingClosureTestsHandler::PlotVHInBins("res", vvHRes, dirOut + "/res");
  UnfoldingClosureTestsHandler::PlotVHInBins("pul", vvHPul, dirOut + "/pul");
  PlotPulResSummary("res", vvHRes, vUnf, dirOut + "/res-summary");
  PlotPulResSummary("pul", vvHPul, vUnf, dirOut + "/pul-summary");

  // plot cross sections from 0th iterations (no fluctuations)
  int dof = nbins;
  std::vector<TH1D*> vDt;
  std::vector<TH1D*> vTh;
  //hGen->SetTitle("POWHEG+PYTHIA");
  hGen->SetLineColor(kBlack);
  hGen->SetTitle("true");
  printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> hGen\n");
  hGen->Print("all");
  vTh.push_back(hGen);
  if(hGenRsp)
  {
    // renormalise
    //hGenRsp->Scale(hGen->Integral() / hGenRsp->Integral());
    printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> hGenRsp\n");
    hGenRsp->Print("all");
    //hGenRsp->SetTitle("MC resp. matr.");
    hGenRsp->SetTitle("MC");
    hGenRsp->SetLineColor(kOrange + 10);
    hGenRsp->SetLineStyle(3);
    vTh.push_back(hGenRsp);
  }
  for(int u = 0; u < nalg; u++)
  {
    double chi2 = utils::ttmd::Chi2(vvHUnf[u][0], vvHCov[u][0], hGen);
    vUnf[u]->SetHistoStyle(vvHUnf[u][0]);
    vvHUnf[u][0]->SetTitle(TString::Format("%s [%.0f/%d]", vUnf[u]->Remark.Data(), chi2, dof));
    vDt.push_back(vvHUnf[u][0]);
  }
  TString csFileName = dirOut + "/xsec0";
  PlotterDiffXSec::XSecParams pars;
  pars.RatMin = 0.6;
  pars.RatMax = 1.4;
  PlotterDiffXSec* plotterDiffXSec = new PlotterDiffXSec();
  plotterDiffXSec->XSec(_unfoldingIOHandler, vDt, vTh, csFileName, pars);
}

TString UnfoldingClosureTestsHandler::GetRootFileNameToStore(UnfoldingIOHandler* unfoldingIOHandler, const TString& suffix, const int niter, const TString& ch, const int modeKinRecoTree)
{
  TString chDir = ch;
  if(modeKinRecoTree)
    chDir += TString::Format("-kr%d/", modeKinRecoTree);
  TString fileName = TString::Format("%s/%s/%s/%s", configHelper_->gPseudoDir.Data(), chDir.Data(), unfoldingIOHandler->Suffix.Data(), suffix.Data());
  if(niter > 0)
    fileName += TString::Format("-i%d", niter);
  fileName += TString::Format(".root");
  return fileName;
}

void UnfoldingClosureTestsHandler::Plot2Cov(TH2D* covEstim, TH2D* covSample, const TString& fileName)
{
  TH2D* corrEstim = utils::ttmd::MakeCorrMatrix(covEstim);
  TH2D* corrSample = utils::ttmd::MakeCorrMatrix(covSample);

  // decrease font if there are many bins
  double scaleFont = 1.0;
  if(corrEstim->GetNbinsX() > 40)
    scaleFont = 0.3;
  else if(corrEstim->GetNbinsX() > 20)
    scaleFont = 0.6;
  corrEstim->SetMarkerSize(corrEstim->GetMarkerSize() * scaleFont);
  corrSample->SetMarkerSize(corrSample->GetMarkerSize() * scaleFont);

  TCanvas c("", "", configHelper_->gPlotPixelSmall, configHelper_->gPlotPixelSmall);
  c.cd();
  gStyle->SetPaintTextFormat("4.2f");
  corrEstim->Draw("text");
  corrSample->SetMarkerColor(2);
  corrSample->Draw("text same");
  utils::ttmd::SaveCanvas(&c, fileName);
}

std::pair<TH1D*, TH2D*> UnfoldingClosureTestsHandler::GetAverageAndCov(const std::vector<TH1D*>& vH)
{
  // check input
  if(vH.size() < 1)
  {
    printf("Error in ZClosureTest::GetAverageAndCov() input vH.size() = %ld\n", vH.size());
    throw;
  }
  const int niter = vH.size();

  // create new histos
  int nbinsXSec = vH[0]->GetNbinsX();
  TH1D* avH = (TH1D*) vH[0]->Clone();
  TH2D* avHCov = new TH2D("", "", nbinsXSec, 0.0, nbinsXSec, nbinsXSec, 0.0, nbinsXSec);

  // calculate average histogram
  // skip 1st entry: zero fluctuation
  for(int i = 1; i < niter; i++)
    avH->Add(vH[i]);
  avH->Scale(1.0 / niter);

  // calculate covariance
  for(int bx = 0; bx <= avHCov->GetNbinsX() + 1; bx++)
    for(int by = 0; by <= avHCov->GetNbinsX() + 1; by++)
      for(int i = 1; i < niter; i++)
        avHCov->SetBinContent(bx, by, avHCov->GetBinContent(bx, by) + (vH[i]->GetBinContent(bx) - avH->GetBinContent(bx)) * (vH[i]->GetBinContent(by) - avH->GetBinContent(by)));
  avHCov->Scale(1.0 / niter);

  // return result
  return std::pair<TH1D*, TH2D*>(avH, avHCov);
}

std::vector<TH1D*> UnfoldingClosureTestsHandler::SetupHistogramsChi2(const int nalg/*, const int niter*/, const int nbinsXSec)
{
  // determine binning
  //int nbinHChi2 = niter / 20;
  double minHChi2 = 0.0;
  double maxHChi2 = nbinsXSec * 4;
  //if(nbinHChi2 < 2)
  //  nbinHChi2 = 2;
  //nbinHChi2 = 50;
  //if(nbinHChi2 > (maxHChi2 - minHChi2))
  //  nbinHChi2 = maxHChi2 - minHChi2;
  int nbinHChi2 = maxHChi2 - minHChi2;

  std::vector<TH1D*> vHChi2;
  for(int u = 0; u < nalg; u++)
    vHChi2.push_back(new TH1D("", "", nbinHChi2, minHChi2, maxHChi2));
  return vHChi2;
}

std::vector<std::vector<TH1D*> > UnfoldingClosureTestsHandler::SetupHistogramsResiduals(const int nalg/*, const int niter*/, const int nbinsXSec)
{
  // steering
  const double eps = 1e-3;
  const double min = -0.15 + eps;
  const double max = 0.15 - eps;
  //int nbins = niter / 20;
  //if(nbins < 2)
  //  nbins = 2;
  //if(nbins > 50)
  //  nbins = 50;
  //nbins = 20;
  int nbins = 50;

  std::vector<std::vector<TH1D*> > vvH;
  for(int i = 0; i < nalg; i++)
  {
    vvH.push_back(std::vector<TH1D*>());
    for(int j = 0; j < nbinsXSec; j++)
      vvH.back().push_back(new TH1D("", "", nbins, min, max));
  }
  return vvH;
}

std::vector<std::vector<TH1D*> > UnfoldingClosureTestsHandler::SetupHistogramsPulls(const int nalg/*, const int niter*/, const int nbinsXSec)
{
  // steering
  const double eps = 1e-3;
  const double min = -3.0 + eps;
  const double max = 3.0 - eps;
  //int nbins = niter / 20;
  //if(nbins < 2)
  //  nbins = 2;
  //if(nbins > 50)
  //  nbins = 50;
  //nbins = 20;
  int nbins = 50;

  std::vector<std::vector<TH1D*> > vvH;
  for(int i = 0; i < nalg; i++)
  {
    vvH.push_back(std::vector<TH1D*>());
    for(int j = 0; j < nbinsXSec; j++)
      vvH.back().push_back(new TH1D("", "", nbins, min, max));
  }
  return vvH;
}

void UnfoldingClosureTestsHandler::PlotVHInBins(const TString& prop, const std::vector<std::vector<TH1D*> >& vvH, const TString& fileName)
{
  //printf("PlotVHInBins prop = %s vvH.size() = %ld vvH[0].size() = %ld\n", prop.Data(), vvH.size(), vvH[0].size());

  // x, y boundaries
  std::pair<double, double> minmax = utils::ttmd::GetMinMaxFromTH1DVector2D(vvH);
  double xmin = vvH[0][0]->GetBinLowEdge(1);
  double xmax = vvH[0][0]->GetBinLowEdge(vvH[0][0]->GetNbinsX() + 1);
  double ymin = 0.0;
  double ymax = minmax.second;
  const double scaleMax = 1.6;
  TH2D hr("", "", 1, xmin, xmax, 1, ymin + 1e-3, ymax * scaleMax);
  double scaleAxisFontSize = 3.0;
  const int div = 3;
  utils::ttmd::ScaleHistoFonts(&hr, scaleAxisFontSize);
  hr.GetXaxis()->SetNdivisions(div);
  hr.GetYaxis()->SetNdivisions(div);

  int nalg = vvH.size();
  int nbins = vvH[0].size();
  std::vector<TH1D*> vhsum;

  TCanvas c("", "", configHelper_->gPlotPixelLarge, configHelper_->gPlotPixelLarge);
  c.SetMargin(0.30, 0.10, 0.30, 0.10);
  int w = -1;
  int h = -1;
  utils::ttmd::DivideSquareWH(&c, nbins, w, h);
  c.DivideSquare(nbins, 0.005, 0.00);

  for(int b = 0; b < nbins; b++)
  {
    std::vector<TH1D*> vh;
    for(int a = 0; a < nalg; a++)
      vh.push_back(vvH[a][b]);
    c.cd(b + 1);
    TH2D* hrClone = (TH2D*)hr.Clone();
    utils::ttmd::AdjustDividedCanvasHistograms(hrClone, b, nbins, w, h);
    PlotVH(prop, vh, "", gPad, hrClone);
  }
  // plot in bins
  utils::ttmd::SaveCanvas(&c, fileName + "-bins");
  // plot total
  std::vector<TH1D*> vh;
  for(int a = 0; a < nalg; a++)
  {
    vh.push_back((TH1D*)vvH[a][0]->Clone());
    for(int b = 1; b < nbins; b++)
      vh.back()->Add(vvH[a][b]);
  }
  PlotVH(prop, vh, fileName + "-all");
}

TVirtualPad* UnfoldingClosureTestsHandler::PlotVH(const TString& prop, const std::vector<TH1D*>& vH, const TString& fileName/* = ""*/,
                                  TVirtualPad* cParent/* = NULL*/, TH2D* hrIn/* = NULL*/)
{
  //printf("PlotVH prop = %s vH.size() = %ld\n", prop.Data(), vH.size());
  const double markerSize = 0.5;
  int nsets = vH[0]->Integral(0, vH[0]->GetNbinsX() + 1);

  TH2D* hr = NULL;
  if(hrIn)
    hr = hrIn;
  else
  {
    std::pair<double, double> minmax = utils::ttmd::GetMinMaxFromTH1DVector(vH);
    double xmin = vH[0]->GetBinLowEdge(1);
    double xmax = vH[0]->GetBinLowEdge(vH[0]->GetNbinsX() + 1);
    double ymin = 0.0;
    double ymax = minmax.second;
    const double scaleMax = 1.6;
    hr = new TH2D("", "", 1, xmin, xmax, 1, ymin, ymax * scaleMax);
  }
  double xmin = hr->GetXaxis()->GetBinLowEdge(1);
  double xmax = hr->GetXaxis()->GetBinLowEdge(hr->GetNbinsX() + 1);
  hr->GetYaxis()->SetTitle("Entries");
  if(prop == "chi")
    hr->GetXaxis()->SetTitle("#chi^{2}");
  else if(prop == "pul")
    hr->GetXaxis()->SetTitle("Pull");
  else if(prop == "res")
    hr->GetXaxis()->SetTitle("Residual");
  else
  {
    printf("Error in ZClosureTest::PlotVH() unknown prop = %s\n", prop.Data());
    throw;
  }

  TVirtualPad* c = NULL;
  if(cParent)
    c = cParent;
  else
    c = new TCanvas("", "", configHelper_->gPlotPixelSmall, configHelper_->gPlotPixelSmall);
  c->cd();
  hr->Draw();

  // draw input histograms
  for(int u = 0; u < _nalg; u++)
  {
    vH[u]->SetMarkerSize(markerSize);
    vH[u]->Draw("hist sames");
    //if(prop == "chi")
    //  vH[u]->Print("all");
  }

  // draw expectation
  TH1D* hEx = NULL;
  int dof = _unfoldingIOHandler->HGenC()->GetNbinsX();
  if(prop == "chi" || prop == "pul")
  {
    int nbins = vH[0]->GetNbinsX();
    hEx = new TH1D("", "", nbins, xmin, xmax);
    TF1 *f = NULL;
    if(prop == "chi")
    {
      f = new TF1("f", "ROOT::Math::chisquared_pdf(x, [0])", xmin, xmax);
      f->SetParameter(0, dof);
      //printf("INTEGRAL = %f\n", f->Mean( 0.0, 100.0));
    }
    else if(prop == "pul")
      f = new TF1("f", "ROOT::Math::gaussian_pdf(x)", xmin, xmax);
    for(int b = 1; b <= nbins; b++)
      hEx->SetBinContent(b, f->Integral(hEx->GetBinLowEdge(b), hEx->GetBinLowEdge(b + 1)) * nsets);
    //printf("hEx->Integral() = %f\n", hEx->Integral());
    hEx->SetLineColor(kGray + 2);
    hEx->SetMarkerColor(kGray + 2);
    hEx->SetMarkerStyle(20);
    hEx->SetMarkerSize(0.5);
    hEx->SetLineStyle(3);
    hEx->Draw("same");
    utils::ttmd::DrawAsGraph(hEx, 0, 0.5, "pc");
    //hEx->Print("all");
  }

  // draw legend
  TLegend* leg = new TLegend(0.40, 0.65, (!cParent) ? 0.88 : 0.96, (!cParent) ? 0.88 : 0.96);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetMargin(0.1);
  //leg->AddEntry((TObject*) NULL, TString::Format("%d pseudo sets", nsets), "");
  for(unsigned int a = 0; a < vH.size(); a++)
  {
    int ndig = (prop == "res") ? 3 : 1;
    TString str;
    str = TString::Format("%s mean = %.*f RMS = %.*f", vH[a]->GetTitle(), ndig, vH[a]->GetMean(), ndig, vH[a]->GetRMS());
    leg->AddEntry(vH[a], str, "l");
  }
  if(prop == "chi")
    leg->AddEntry(hEx, TString::Format("%s mean = %.1f RMS = %.1f", "expect", (double)dof, TMath::Sqrt((double)2 * dof)), "lp");
  if(prop == "pul")
    leg->AddEntry(hEx, TString::Format("%s mean = %.1f RMS = %.1f", "expect", 0.0, 1.0), "lp");
  leg->Draw();

  hr->Draw("axis same");
  if(fileName != "")
    utils::ttmd::SaveCanvas(c, fileName);
  return c;
}

void UnfoldingClosureTestsHandler::PlotPulResSummary(const TString& prop, const std::vector<std::vector<TH1D*> >& vvH, std::vector<ZUnfold*> vUnf, const TString& fileName)
{
  double ymin = 0.0;
  double ymax = 0.0;
  if(prop == "pul")
  {
    ymin = -3.0;
    ymax = 4.5;
  }
  if(prop == "res")
  {
    ymin = -0.05;
    ymax = 0.15;
  }
  const double markerSize = 0.5;
  //int nalg = vUnf.size();
  int nbins = vvH[0].size();
  double xmin = 0.5;
  double xmax = nbins + 0.5;
  TH2D hr("", "", 1, xmin, xmax, 1, ymin, ymax);
  hr.GetXaxis()->SetTitle("bin");
  if(prop == "pul")
    hr.GetYaxis()->SetTitle("pull mean/RMS ");
  if(prop == "res")
    hr.GetYaxis()->SetTitle("residual mean/RMS ");
  TCanvas c("", "", configHelper_->gPlotPixelSmall, configHelper_->gPlotPixelSmall);
  c.cd();
  hr.Draw();
  TH1D hpulEmpty("", "", nbins, xmin, xmax);
  TH1D hresEmpty("", "", nbins, xmin, xmax);
  TLegend leg(0.15, 0.70, 0.85, 0.88);
  leg.SetNColumns(2);
  leg.AddEntry((TObject*)(NULL), "Mean", "");
  leg.AddEntry((TObject*)(NULL), "RMS", "");
  for(unsigned int u = 0; u < vUnf.size(); u++)
  {
    TH1D* hpul = (TH1D*) hpulEmpty.Clone();
    TH1D* hres = (TH1D*) hresEmpty.Clone();
    for(int b = 0; b < nbins; b++)
    {
      hres->SetBinContent(b + 1, vvH[u][b]->GetMean());
      hpul->SetBinContent(b + 1, vvH[u][b]->GetRMS());
    }
    vUnf[u]->SetHistoStyle(hpul);
    vUnf[u]->SetHistoStyle(hres);
    hpul->SetMarkerSize(markerSize);
    hres->SetMarkerSize(markerSize);
    hres->SetMarkerStyle(hres->GetMarkerStyle() - 4);
    leg.AddEntry(hres, vUnf[u]->Remark, "p");
    leg.AddEntry(hpul, vUnf[u]->Remark, "p");
    utils::ttmd::DrawAsGraph(hpul);
    utils::ttmd::DrawAsGraph(hres);
  }
  TH1D hpulExpect("", "", nbins, xmin, xmax);
  TH1D hresExpect("", "", nbins, xmin, xmax);
  for(int b = 0; b < nbins; b++)
  {
    hresEmpty.SetBinContent(b + 1, 0.0);
    hpulEmpty.SetBinContent(b + 1, 1.0);
  }
  hpulEmpty.SetLineColor(1);
  hresEmpty.SetLineColor(1);
  hpulEmpty.Draw("hist0 same");
  hresEmpty.Draw("hist0 same");
  leg.Draw();
  utils::ttmd::SaveCanvas(&c, fileName);
}
