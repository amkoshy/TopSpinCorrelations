#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TChain.h>
#include <TCanvas.h>
#include <iostream>
#include "TLorentzVector.h"
#include "TMath.h"
#include <vector>
#include <TLegend.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include "TStyle.h"
#include <string>
#include "TLatex.h"
#include "TPaveText.h"
#include "THStack.h"
#include <sstream>
#include <iomanip>
#include "Math/QuantFuncMathCore.h"
#include "../../common/include/CommandLineParameters.h"
#include "TSystem.h"


 
using namespace std;

// Code to make histograms for trigger efficiencies. To run in root:
// ex: root -l -b makehistograms.cc++(\"2017_UL\")
// ex: root -l -b makehistograms.cc++'("2017_UL")'


string era_;


/// Functions

double CalcLumiWeight(TFile* fMC_, double XSection_, string era_) {

  double lumi = 0;
  
  if (era_ == "2016")  lumi = 35922.0;
  if (era_ == "2016preVFP_UL")  lumi = 19500.0;
  if (era_ == "2016postVFP_UL")  lumi = 16800.0;
  else if (era_ == "2016preVFP_UL")  lumi = 19667.81;
  else if (era_ == "2016postVFP_UL")  lumi = 16977.7;
  else if (era_ == "2017" || era_ == "2017_UL")  lumi = 41479.68;
  else if (era_ == "2018" || era_ == "2018_UL")  lumi = 59832.48;

  TH1D *h_NrOfEvts = (TH1D*)fMC_->Get("weightedEvents");  
  
  double NrOfEvts = h_NrOfEvts->GetBinContent(1);
  //  double NrOfEvts = h_NrOfEvts->GetEntries();

  double lumiWeight = lumi * XSection_ / NrOfEvts;


  std::cout << ", lumiw=" << lumiWeight <<", lumi=" << lumi <<", XSection=" << XSection_ <<", NrOfEvts=" << NrOfEvts << "\n";
 

  return lumiWeight;

}


void ApplyFlatWeights(TH1* hist, const double weight) {
  
  if (weight == 0) std::cout << "Warning: the weight your applying is 0. This will remove your distribution." << std::endl;
  hist->Scale(weight);

}


void setPad(TPad *p, bool hasGridX_ = true) {

  p->SetLeftMargin(0.12);
  p->SetRightMargin(0.05);
  p->SetBottomMargin(0.2);
  p->SetTopMargin(0.1);

  if (hasGridX_) p->SetGridx();
  p->SetGridy();
  p->SetTickx(1);
  p->SetTicky(1);

  p->Range(80,0.75,220,1.25);
  
  p->Draw();
  p->cd();

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("3.3f");

}


void set2DPad(TPad *p) {

  p->SetLeftMargin(0.1);
  p->SetRightMargin(0.125);
  p->SetBottomMargin(0.15);
  p->SetTopMargin(0.1);
  p->SetTickx(1);
  p->SetTicky(1);
  p->Draw();
  p->cd();

  gStyle->SetOptStat(0);

}


void draw_CMS_lumi_ch(string era, string channel_) {

  TPaveText *pt_CMS = new TPaveText(0.13,0.9,0.23,1.0,"brNDC");
  pt_CMS->SetName("");
  pt_CMS->SetBorderSize(0);
  pt_CMS->SetFillStyle(0);
  pt_CMS->SetTextAlign(12);
  pt_CMS->SetTextFont(61);
  pt_CMS->SetTextSize(0.044);
  pt_CMS->AddText("CMS");
  pt_CMS->Draw("SAME"); 

  TPaveText *pt_Prelim = new TPaveText(0.23,0.895,0.33,0.995,"brNDC");
  pt_Prelim->SetName("");
  pt_Prelim->SetBorderSize(0);
  pt_Prelim->SetFillStyle(0);
  pt_Prelim->SetTextAlign(12);
  pt_Prelim->SetTextFont(52);
  pt_Prelim->SetTextSize(0.04);
  pt_Prelim->AddText("Preliminary");
  pt_Prelim->Draw("SAME"); 

  TPaveText *pt_Lumi = new TPaveText(0.76,0.9,0.96,1.0,"brNDC");
  pt_Lumi->SetName("");
  pt_Lumi->SetBorderSize(0);
  pt_Lumi->SetFillStyle(0);
  pt_Lumi->SetTextAlign(32);
  pt_Lumi->SetTextFont(42);
  pt_Lumi->SetTextSize(0.04);
  if (era == "2017") pt_Lumi->AddText("41.53 fb^{-1} (13 TeV)");
  else if (era == "2018" ) pt_Lumi->AddText("59.97 fb^{-1} (13 TeV)");
  else if (era == "2016" ) pt_Lumi->AddText("35.9 fb^{-1} (13 TeV)");

  else if (era == "2016preVFP_UL" ) pt_Lumi->AddText("19.7 fb^{-1} (13 TeV)");
  else if (era == "2016postVFP_UL" ) pt_Lumi->AddText("17.0 fb^{-1} (13 TeV)");
  else if (era == "2017_UL" ) pt_Lumi->AddText("41.5 fb^{-1} (13 TeV)");
  else if (era == "2018_UL" ) pt_Lumi->AddText("59.8 fb^{-1} (13 TeV)");
  pt_Lumi->Draw("SAME"); 

  TLatex tex_ch;
  tex_ch.SetTextAngle(0);
  tex_ch.SetTextColor(kBlack);
  tex_ch.SetTextAlign(31);
  tex_ch.SetTextSize(0.06);
  if (channel_ == "ee") {
    tex_ch.DrawLatexNDC(0.25,0.8,"ee");
  } else if (channel_ == "emu") {
    tex_ch.DrawLatexNDC(0.25,0.8,"e#mu");
  } else if (channel_ == "mumu") {
    tex_ch.DrawLatexNDC(0.25,0.8,"#mu#mu");
  }

}

void PropagateKatzStatError(TH1F *h_x, TH1F *h_m, TH1F *h_y, TH1F *h_n, TH1F *h_t, TGraphAsymmErrors *g_t  ) {


  int nbins = h_t->GetNbinsX();;

  for (Int_t nbin = 0; nbin < nbins; nbin++) {

    double x = h_x->GetBinContent(nbin+1);
    double m = h_m->GetBinContent(nbin+1);
    double y = h_y->GetBinContent(nbin+1);
    double n = h_n->GetBinContent(nbin+1);
    double t = h_t->GetBinContent(nbin+1);

     if ( x == m ) x = x - 0.5;
     if ( y == n ) y = y - 0.5;
    
    double var = (1/x) - (1/m) + (1/y) - (1/n);

    // 68-95-99.7 rule: 1 sigma  
    double alpha = 0.3173;

    double chi =  ROOT::Math::normal_quantile(1 - alpha/2,1);
    
    double t_yl = t * std::exp( -1 * chi * std::sqrt(var));
    double t_yh = t * std::exp( +1 * chi * std::sqrt(var));

    double t_xl = g_t->GetErrorXlow(nbin);
    double t_xh = g_t->GetErrorXhigh(nbin);

    double e = (std::abs(t - t_yl) + std::abs(t - t_yh))/2;

    h_t->SetBinError(nbin+1, e );


    //    g_t->SetPointY(nbin,h_t->GetBinError(nbin+1));
    g_t->SetPointError(nbin,t_xl,t_xh,std::abs(t - t_yl),std::abs(t - t_yh));

  }

}


void PropagateKatzStatError2D(TH2F *h2D_x, TH2F *h2D_m, TH2F *h2D_y, TH2F *h2D_n, TH2F *h2D_t) {


  int nbins_X = h2D_t->GetNbinsX();
  int nbins_Y = h2D_t->GetNbinsY();

    for (Int_t nbin_X = 0; nbin_X < nbins_X; nbin_X++) {

      for (Int_t nbin_Y = 0; nbin_Y < nbins_Y; nbin_Y++) {
	
	double x = h2D_x->GetBinContent(nbin_X+1,nbin_Y+1);
	double m = h2D_m->GetBinContent(nbin_X+1,nbin_Y+1);
	double y = h2D_y->GetBinContent(nbin_X+1,nbin_Y+1);
	double n = h2D_n->GetBinContent(nbin_X+1,nbin_Y+1);
	double t = h2D_t->GetBinContent(nbin_X+1,nbin_Y+1);
	
	if ( x == m ) x = x - 0.5;
	if ( y == n ) y = y - 0.5;
    
	double var = (1/x) - (1/m) + (1/y) - (1/n);

	// 68-95-99.7 rule: 1 sigma
	double alpha = 0.3173;

	double chi =  ROOT::Math::normal_quantile(1 - alpha/2,1);
    
	double t_yl = t * std::exp( -1 * chi * std::sqrt(var));
	double t_yh = t * std::exp( +1 * chi * std::sqrt(var));


	double e = (std::abs(t - t_yl) + std::abs(t - t_yh))/2;

	h2D_t->SetBinError(nbin_X+1,nbin_Y+1, e );

      }

    }

}

void drawGraph(TH1F *h1, TH1F *h2, TGraphAsymmErrors *g1, TGraphAsymmErrors *g2, string yAxisTitle_, string channel_, double xMin_, double xMax_, string era) {


  g1->SetLineWidth(2);
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(1.2);
  g1->SetMarkerColor(kAzure);
  g1->SetLineColor(kAzure);
  g2->SetLineWidth(2);
  g2->SetMarkerStyle(20);
  g2->SetMarkerSize(1.2);
  g2->SetMarkerColor(kRed);
  g2->SetLineColor(kRed);

  h1->SetLineWidth(2);
  h1->SetLineColor(kAzure);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1.2);
  h1->SetMarkerColor(kAzure);

  h1->GetXaxis()->SetLabelSize(0.);
  h1->GetYaxis()->SetLabelSize(0.03);
  h1->GetYaxis()->SetLabelOffset(0.01);

  h1->GetYaxis()->SetTitleSize(0.04);
  h1->GetYaxis()->SetTitleOffset(1.12);

  h1->GetYaxis()->SetTitle((yAxisTitle_).c_str());
  h1->GetXaxis()->SetLimits(xMin_, xMax_); 
  if (yAxisTitle_ == "Efficiency") {
    h1->SetMinimum(0.8);
    h1->SetMaximum(1.1);
  }
  if (yAxisTitle_ == "Entries") {
    h1->GetYaxis()->SetTitleOffset(1.5);
  }

  h2->SetLineWidth(2);
  h2->SetLineColor(kRed);
  h2->SetMarkerStyle(20);
  h2->SetMarkerSize(1.2);
  h2->SetMarkerColor(kRed);

  h1->Draw("X");
  h2->Draw("XSAME");
   g1->Draw("EPSAME");
   g2->Draw("EPSAME");

   TH1F *uncBand_ = new TH1F();
   uncBand_->SetFillColor(1);
   uncBand_->SetFillStyle(3354);


  TLegend *leg = new TLegend(0.68,0.78,0.93,0.88);  
  if (yAxisTitle_ == "Efficiency") {
    //    leg->AddEntry(h1,Form("Data (%i)", era), "lep");
    leg->AddEntry(h1,"Data", "lep");
    leg->AddEntry(h2,"MC", "lep");
    leg->AddEntry(uncBand_,"Uncertainty", "f");
  } else {
    leg->AddEntry(h1,"Dilepton Triggers", "lep");
    leg->AddEntry(h2,"Cross Triggers", "lep");
  }
  leg->SetBorderSize(0);
	leg->SetFillColor(0);
  leg->Draw("SAME");

  draw_CMS_lumi_ch(era, channel_);

}


void drawGraph(TH1F *h1, TH1F *h2, string yAxisTitle_, string channel_, double xMin_, double xMax_, string era) {

  h1->SetLineWidth(2);
  h1->SetLineColor(kAzure);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1.2);
  h1->SetMarkerColor(kAzure);

  h1->GetXaxis()->SetLabelSize(0.);
  h1->GetYaxis()->SetLabelSize(0.03);
  h1->GetYaxis()->SetLabelOffset(0.01);

  h1->GetYaxis()->SetTitleSize(0.04);
  h1->GetYaxis()->SetTitleOffset(1.12);

  h1->GetYaxis()->SetTitle((yAxisTitle_).c_str());
  h1->GetXaxis()->SetLimits(xMin_, xMax_); 
  if (yAxisTitle_ == "Efficiency") {
    h1->SetMinimum(0.8);
    h1->SetMaximum(1.1);
  }
  if (yAxisTitle_ == "Entries") {
    h1->GetYaxis()->SetTitleOffset(1.5);
  }

  h1->Draw("");

  h2->SetLineWidth(2);
  h2->SetLineColor(kRed);
  h2->SetMarkerStyle(20);
  h2->SetMarkerSize(1.2);
  h2->SetMarkerColor(kRed);
  h2->Draw("SAME");

  TLegend *leg = new TLegend(0.68,0.78,0.93,0.88);  
  if (yAxisTitle_ == "Efficiency") {
    //    leg->AddEntry(h1,Form("Data (%i)", era), "lep");
    leg->AddEntry(h1,"Data", "lep");
    leg->AddEntry(h2,"MC", "lep");
  } else {
    leg->AddEntry(h1,"Dilepton Triggers", "lep");
    leg->AddEntry(h2,"Cross Triggers", "lep");
  }
  leg->SetBorderSize(0);
	leg->SetFillColor(0);
  leg->Draw("SAME");

  draw_CMS_lumi_ch(era, channel_);

}


void drawGraph_alpha(TH1F *g1, TH1F *g2, TH1F *g3, string channel_, double xMin_, double xMax_, string era) {

  g1->SetLineWidth(2);
  g1->SetLineColor(kAzure);
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(1.2);
  g1->SetMarkerColor(kAzure);

  g1->GetXaxis()->SetLabelSize(0.);
  g1->GetYaxis()->SetLabelSize(0.03);
  g1->GetYaxis()->SetLabelOffset(0.01);

  g1->GetYaxis()->SetTitleSize(0.04);
  g1->GetYaxis()->SetTitleOffset(0.99);
  g1->GetYaxis()->SetTitle("Events");

  g1->GetXaxis()->SetLimits(xMin_, xMax_); 
  g1->SetMinimum(0.0);
  g1->SetMaximum(1.5);
  g1->Draw();

  g2->SetLineWidth(2);
  g2->SetLineColor(kRed);
  g2->SetMarkerStyle(20);
  g2->SetMarkerSize(1.2);
  g2->SetMarkerColor(kRed);
  g2->Draw("SAME");

  g3->SetLineWidth(2);
  g3->SetLineColor(kGreen);
  g3->SetMarkerStyle(20);
  g3->SetMarkerSize(1.2);
  g3->SetMarkerColor(kGreen);
  g3->Draw("SAME");

  TLegend *leg = new TLegend(0.68,0.78,0.93,0.88);  
  leg->AddEntry(g1, "MET", "lep");
  leg->AddEntry(g2, "Dilepton", "lep");
  leg->AddEntry(g3, "MET+Dilepton", "lep");
  leg->SetBorderSize(0);
	leg->SetFillColor(0);
  leg->Draw("SAME");

  draw_CMS_lumi_ch(era, channel_);

}


void DrawUncertainties(TH1F *fullsyst, TH1F *stat, TH1F *njets, TH1F *nvtx, TH1F *era_dependency, string channel_, string xAxisTitle_, double xMin_, double xMax_, string era) {

  fullsyst->SetLineWidth(2);
  fullsyst->SetLineColor(kAzure);
  fullsyst->SetMarkerStyle(20);
  fullsyst->SetMarkerSize(1.2);
  fullsyst->SetMarkerColor(kAzure);

  fullsyst->GetXaxis()->SetLabelSize(0.03);
  fullsyst->GetXaxis()->SetTitleSize(0.04);
  fullsyst->GetXaxis()->SetTitleOffset(1.2);


  if (channel_ == "emu" && xAxisTitle_ == "leading lepton pT (GeV)") {
    fullsyst->GetXaxis()->SetTitle(("Electron pT (GeV)"));
  }
  else if (channel_ == "emu" && xAxisTitle_ == "subleading lepton pT (GeV)") {
    fullsyst->GetXaxis()->SetTitle(("Muon pT (GeV)"));
  }
  else if (channel_ == "emu" && xAxisTitle_ == "leading lepton #eta") {
    fullsyst->GetXaxis()->SetTitle(("Electron #eta"));
  }
  else if (channel_ == "emu" && xAxisTitle_ == "subleading lepton #eta") {
    fullsyst->GetXaxis()->SetTitle(("Muon #eta"));
  }
  else {
    fullsyst->GetXaxis()->SetTitle((xAxisTitle_).c_str());
  }

  fullsyst->GetYaxis()->SetLabelSize(0.03);
  fullsyst->GetYaxis()->SetLabelOffset(0.01);

  fullsyst->GetYaxis()->SetTitleSize(0.04);
  fullsyst->GetYaxis()->SetTitleOffset(1.25);
  fullsyst->GetYaxis()->SetTitle("Uncertainty (%)");

  fullsyst->GetXaxis()->SetLimits(xMin_, xMax_); 


  float AxisMin = 0.0;
  float AxisMax = 0.01;

  int nbins_X = fullsyst->GetNbinsX();

  TH1F *totalstatsyst = (TH1F*)fullsyst->Clone("h_totalstatsyst");

  for (Int_t nbin_X = 0; nbin_X < nbins_X; nbin_X++) {
    double statbinContent = stat->GetBinContent(nbin_X+1);
    double fullsystbinContent = fullsyst->GetBinContent(nbin_X+1);
    double totalstatsystbinContent = std::sqrt((statbinContent)*(statbinContent) + (fullsystbinContent)*(fullsystbinContent));
    totalstatsyst->SetBinContent(nbin_X+1, totalstatsystbinContent );

    if (totalstatsystbinContent > AxisMax) AxisMax = totalstatsystbinContent;   

  }

    
  fullsyst->SetMinimum(AxisMin);
  fullsyst->SetMaximum(1.25*AxisMax);

  fullsyst->SetFillStyle(0);
  fullsyst->Draw("Hist");

  stat->SetLineWidth(2);
  stat->SetLineColor(kRed);
  stat->SetMarkerStyle(20);
  stat->SetMarkerSize(1.2);
  stat->SetMarkerColor(kRed);
  stat->SetFillStyle(0);
  stat->Draw("SAME,Hist");

  totalstatsyst->SetLineWidth(2);
  totalstatsyst->SetLineColor(kBlack);
  totalstatsyst->SetMarkerStyle(20);
  totalstatsyst->SetMarkerSize(1.2);
  totalstatsyst->SetMarkerColor(kBlack);
  totalstatsyst->SetFillStyle(0);
  totalstatsyst->Draw("SAME,Hist");

  njets->SetLineWidth(2);
  njets->SetLineColor(kGreen);
  njets->SetMarkerStyle(20);
  njets->SetMarkerSize(1.2);
  njets->SetMarkerColor(kGreen);
  njets->SetFillStyle(0);
  njets->Draw("SAME,Hist");

  nvtx->SetLineWidth(2);
  nvtx->SetLineColor(kCyan);
  nvtx->SetMarkerStyle(20);
  nvtx->SetMarkerSize(1.2);
  nvtx->SetMarkerColor(kCyan);
  nvtx->Draw("SAME,Hist");

  era_dependency->SetLineWidth(2);
  era_dependency->SetLineColor(kMagenta);
  era_dependency->SetMarkerStyle(20);
  era_dependency->SetMarkerSize(1.2);
  era_dependency->SetMarkerColor(kMagenta);
  era_dependency->SetFillStyle(0);
  era_dependency->Draw("SAME,Hist");

  TLegend *leg = new TLegend(0.3,0.78,0.5,0.88);  
  leg->AddEntry(fullsyst, "Total Systematic", "lep");
  leg->AddEntry(totalstatsyst, "Total Systematic + Statistical", "lep");
  leg->AddEntry(stat, "Total Statistical", "lep");
  leg->AddEntry(njets, "Number of Jets", "lep");
  leg->AddEntry(nvtx, "Number of Vertices", "lep");
  leg->AddEntry(era_dependency, "Era Dependency", "lep");
  leg->SetBorderSize(0);
	leg->SetFillColor(0);
  leg->Draw("SAME");

  draw_CMS_lumi_ch(era, channel_);

}



void drawRatio(TH1F *h, TGraphAsymmErrors *g, string channel_, string xAxisTitle_, double xMin_, double xMax_, string yAxisTitle_, double yMin_, double yMax_) {
  
  h->SetMinimum(yMin_); 
  h->SetMaximum(yMax_);
  h->GetXaxis()->SetLimits(xMin_, xMax_); 

  h->GetYaxis()->SetNdivisions(505);
  h->GetXaxis()->SetNdivisions(510);

  h->GetXaxis()->SetLabelSize(0.1);
  h->GetXaxis()->SetTitleSize(0.1);
  h->GetXaxis()->SetTitleOffset(0.95);

  if (channel_ == "emu" && xAxisTitle_ == "leading lepton pT (GeV)") {
    h->GetXaxis()->SetTitle(("Electron pT (GeV)"));
  }
  else if (channel_ == "emu" && xAxisTitle_ == "subleading lepton pT (GeV)") {
    h->GetXaxis()->SetTitle(("Muon pT (GeV)"));
  }
  else if (channel_ == "emu" && xAxisTitle_ == "leading lepton #eta") {
    h->GetXaxis()->SetTitle(("Electron #eta"));
  }
  else if (channel_ == "emu" && xAxisTitle_ == "subleading lepton #eta") {
    h->GetXaxis()->SetTitle(("Muon #eta"));
  }
  else {
    h->GetXaxis()->SetTitle((xAxisTitle_).c_str());
  }

  h->GetYaxis()->SetLabelSize(0.1);
  h->GetYaxis()->SetTitleSize(0.12);
  h->GetYaxis()->SetTitleOffset(0.37);
  h->GetYaxis()->SetTitle((yAxisTitle_).c_str());

  h->SetLineWidth(2);
  h->SetLineColor(kAzure);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.2);
  h->SetMarkerColor(kAzure);
  h->Draw("X");

  g->SetLineWidth(2);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.2);
  g->SetMarkerColor(kAzure);
  g->SetLineColor(kAzure);

  g->Draw("EPSAME");
    
  g->Fit("pol0", "W");
  TF1 *func = g->GetFunction("pol0");
  Double_t p0 = func->GetParameter(0);
  Double_t e0 = func->GetParError(0);
  // cout << p0 << " +/- " << e0 << endl;

  TLatex tex;
  tex.SetTextAngle(0);
  tex.SetTextColor(kBlack);
  tex.SetTextAlign(31);
  tex.SetTextSize(0.12);
  if (yAxisTitle_ == "Alpha") {
    tex.DrawLatexNDC(0.5,0.75,Form("%.4f #pm %.4f", p0, e0));
  } else {
    tex.DrawLatexNDC(0.5,0.75,Form("%.3f #pm %.3f", p0, e0));
  }

  TLine *l = new TLine(xMin_, p0, xMax_, p0); 
  l->SetLineWidth(2);
  l->SetLineColor(kRed);
  l->Draw("SAME");
    

}


void drawRatio(TH1F *h, string channel_, string xAxisTitle_, double xMin_, double xMax_, string yAxisTitle_, double yMin_, double yMax_) {
  
  h->SetMinimum(yMin_); 
  h->SetMaximum(yMax_);
  h->GetXaxis()->SetLimits(xMin_, xMax_); 

  h->GetYaxis()->SetNdivisions(505);
  h->GetXaxis()->SetNdivisions(510);

  h->GetXaxis()->SetLabelSize(0.1);
  h->GetXaxis()->SetTitleSize(0.1);
  h->GetXaxis()->SetTitleOffset(0.95);

  if (channel_ == "emu" && xAxisTitle_ == "leading lepton pT (GeV)") {
    h->GetXaxis()->SetTitle(("Electron pT (GeV)"));
  }
  else if (channel_ == "emu" && xAxisTitle_ == "subleading lepton pT (GeV)") {
    h->GetXaxis()->SetTitle(("Muon pT (GeV)"));
  }
  else if (channel_ == "emu" && xAxisTitle_ == "leading lepton #eta") {
    h->GetXaxis()->SetTitle(("Electron #eta"));
  }
  else if (channel_ == "emu" && xAxisTitle_ == "subleading lepton #eta") {
    h->GetXaxis()->SetTitle(("Muon #eta"));
  }
  else {
    h->GetXaxis()->SetTitle((xAxisTitle_).c_str());
  }

  h->GetYaxis()->SetLabelSize(0.1);
  h->GetYaxis()->SetTitleSize(0.12);
  h->GetYaxis()->SetTitleOffset(0.37);
  h->GetYaxis()->SetTitle((yAxisTitle_).c_str());

  h->SetLineWidth(2);
  h->SetLineColor(kAzure);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.2);
  h->SetMarkerColor(kAzure);
  h->Draw("");

  
  h->Fit("pol0");
  TF1 *func = h->GetFunction("pol0");
  Double_t p0 = func->GetParameter(0);
  Double_t e0 = func->GetParError(0);
  // cout << p0 << " +/- " << e0 << endl;

  TLatex tex;
  tex.SetTextAngle(0);
  tex.SetTextColor(kBlack);
  tex.SetTextAlign(31);
  tex.SetTextSize(0.12);
  if (yAxisTitle_ == "Alpha") {
    tex.DrawLatexNDC(0.5,0.75,Form("%.4f #pm %.4f", p0, e0));
  } else {
    tex.DrawLatexNDC(0.5,0.75,Form("%.3f #pm %.3f", p0, e0));
  }

  TLine *l = new TLine(xMin_, p0, xMax_, p0); 
  l->SetLineWidth(2);
  l->SetLineColor(kRed);
  l->Draw("SAME");
    

}


void drawMultivariables(TH2F *h, string channel_, string multivariablesxAxisTitles_, string multivariablesyAxisTitles_, string era) {

  gStyle->SetPalette(1);

  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.035);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleSize(0.035);
  h->GetXaxis()->SetTitleOffset(1.3);

  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.035);
  h->GetYaxis()->SetLabelOffset(0.01);
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleSize(0.035);
  h->GetYaxis()->SetTitleOffset(1.3);

  h->GetZaxis()->SetLabelFont(42);
  h->GetZaxis()->SetLabelSize(0.035);
  h->GetZaxis()->SetLabelOffset(0.005);
  h->GetZaxis()->SetTitleFont(42);
  h->GetZaxis()->SetTitleSize(0.035);

  if (channel_ == "emu" && multivariablesxAxisTitles_ == "leading lepton pT (GeV)" && multivariablesyAxisTitles_ == "subleading lepton pT (GeV)") {
    h->GetXaxis()->SetTitle(("Electron pT (GeV)"));
    h->GetYaxis()->SetTitle(("Muon pT (GeV)"));
    h->SetMarkerSize(0.7);
  }
  else if (channel_ == "emu" && multivariablesxAxisTitles_ == "leading lepton #eta" && multivariablesyAxisTitles_ == "subleading lepton #eta") {
    h->GetXaxis()->SetTitle(("Electron #eta"));
    h->GetYaxis()->SetTitle(("Muon #eta"));
  }
  else {
    h->GetXaxis()->SetTitle((multivariablesxAxisTitles_).c_str());
    h->GetYaxis()->SetTitle((multivariablesyAxisTitles_).c_str());
  }

  if (multivariablesxAxisTitles_ == "leading lepton pT (GeV)") {

    float AxisMin = 0.99;
    float AxisMax = 0.01;

    int nbins_X = h->GetNbinsX();
    int nbins_Y = h->GetNbinsY();

    for (Int_t nbin_X = 0; nbin_X < nbins_X; nbin_X++) {

      for (Int_t nbin_Y = 0; nbin_Y < nbins_Y; nbin_Y++) {

	double binContent = h->GetBinContent(nbin_X+1,nbin_Y+1);

	if (binContent < AxisMin && binContent > 0.0) AxisMin = binContent;
	if (binContent > AxisMax) AxisMax = binContent;

      }
    }
    
    h->SetMinimum(AxisMin);
    h->SetMaximum(AxisMax);
  }


  //  h->Draw("colz TEXT error");
  h->Draw("colz TEXT");

  draw_CMS_lumi_ch(era, channel_);

}

void AddSystError_njet_nvtx(TH2F *h2D_SF, TH2F *h2D_SF_njets_high, TH2F *h2D_SF_njets_low, TH2F *h2D_SF_nvtx_high, TH2F *h2D_SF_nvtx_low) {

  int nbins_X = h2D_SF->GetNbinsX();
  int nbins_Y = h2D_SF->GetNbinsY();
  
  std::vector<std::vector<double>> vec_var(nbins_X, std::vector<double> (nbins_Y, 0)); 

  for (Int_t nbin_X = 0; nbin_X < nbins_X; nbin_X++) {

    for (Int_t nbin_Y = 0; nbin_Y < nbins_Y; nbin_Y++) {

      double binContent = h2D_SF->GetBinContent(nbin_X+1,nbin_Y+1);
      double binError = h2D_SF->GetBinError(nbin_X+1,nbin_Y+1);
      if (binContent > 1e-6) vec_var[nbin_X][nbin_Y] += (binError/binContent) * (binError/binContent);
      else vec_var[nbin_X][nbin_Y] += 0*0;

      //      cout << "Check Stat Uncertainty: " << binContent*std::sqrt(vec_var[nbin_X][nbin_Y]) << endl;
      
      double rel_diff_njets_high;
      double rel_diff_njets_low;
      double rel_diff_nvtx_high;
      double rel_diff_nvtx_low;
      
      if (h2D_SF_njets_high->GetBinContent(nbin_X+1,nbin_Y+1) < 1e-6) rel_diff_njets_high = 0;
      else rel_diff_njets_high = std::fabs(h2D_SF_njets_high->GetBinContent(nbin_X+1,nbin_Y+1) - binContent) / binContent;
      if (h2D_SF_njets_low->GetBinContent(nbin_X+1,nbin_Y+1) < 1e-6) rel_diff_njets_low = 0;
      else rel_diff_njets_low = std::fabs(h2D_SF_njets_low->GetBinContent(nbin_X+1,nbin_Y+1) - binContent) / binContent;
      if (h2D_SF_nvtx_high->GetBinContent(nbin_X+1,nbin_Y+1) < 1e-6) rel_diff_nvtx_high = 0;
      else rel_diff_nvtx_high = std::fabs(h2D_SF_nvtx_high->GetBinContent(nbin_X+1,nbin_Y+1) - binContent) / binContent;
      if (h2D_SF_nvtx_low->GetBinContent(nbin_X+1,nbin_Y+1) < 1e-6) rel_diff_nvtx_low = 0;
      else rel_diff_nvtx_low = std::fabs(h2D_SF_nvtx_low->GetBinContent(nbin_X+1,nbin_Y+1) - binContent) / binContent;

      if (binContent < 1e-6) {
	rel_diff_njets_high = 0;
	rel_diff_njets_low = 0;
	rel_diff_nvtx_high = 0;
	rel_diff_nvtx_low = 0;
      }

      double rel_diff_njets = std::max(rel_diff_njets_high,rel_diff_njets_low);
      double rel_diff_nvtx = std::max(rel_diff_nvtx_high,rel_diff_nvtx_low);
      //      cout << "rel_diff_njets: " << rel_diff_njets << "rel_diff_nvtx " << rel_diff_nvtx << endl;


      vec_var[nbin_X][nbin_Y] += (rel_diff_njets) * (rel_diff_njets);
      vec_var[nbin_X][nbin_Y] += (rel_diff_nvtx) * (rel_diff_nvtx);

      h2D_SF->SetBinError(nbin_X+1, nbin_Y+1, binContent*std::sqrt(vec_var[nbin_X][nbin_Y]) );
      //      cout << "Checking Final Uncertainty: " << h2D_SF->GetBinError(nbin_X+1,nbin_Y+1) << endl;

    }

  }

}


void AddSystError_EraRun(TH2F *h2D_SF, vector<TH2F*> v_2D_EraRuns,  vector<float> v_2D_EraRuns_Scale, TH2F *h2D_SF_Errors) {

  TH2F *h2D_SF_lumiweighted = (TH2F*)v_2D_EraRuns[0]->Clone("h2D_SF_lumiweighted");
  h2D_SF_lumiweighted->Scale(v_2D_EraRuns_Scale[0]);

  for (UInt_t i = 1; i < v_2D_EraRuns.size(); i++) {	
    v_2D_EraRuns[i]->Scale(v_2D_EraRuns_Scale[i]);
    h2D_SF_lumiweighted->Add(v_2D_EraRuns[i]);
  }

  int nbins_X = h2D_SF->GetNbinsX();
  int nbins_Y = h2D_SF->GetNbinsY();
  
  std::vector<std::vector<double>> vec_var(nbins_X, std::vector<double> (nbins_Y, 0)); 

  for (Int_t nbin_X = 0; nbin_X < nbins_X; nbin_X++) {

    for (Int_t nbin_Y = 0; nbin_Y < nbins_Y; nbin_Y++) {

      double binContent = h2D_SF->GetBinContent(nbin_X+1,nbin_Y+1);

      double binContent_lumiweighted = h2D_SF_lumiweighted->GetBinContent(nbin_X+1,nbin_Y+1);

      double rel_diff;

      if (binContent_lumiweighted < 1e-6) rel_diff = 0;
      else rel_diff = std::fabs(binContent - binContent_lumiweighted) / binContent;

      double binError = h2D_SF->GetBinError(nbin_X+1,nbin_Y+1);
      h2D_SF_Errors->SetBinError(nbin_X+1, nbin_Y+1, 0 );

      if (binContent > 1e-6) vec_var[nbin_X][nbin_Y] += (binError/binContent) * (binError/binContent);
      else vec_var[nbin_X][nbin_Y] += 0*0;


      /*
      for (UInt_t i = 0; i < v_2D_EraRuns.size(); i++) {

	if (rel_diff < std::fabs(v_2D_EraRuns.at(i)->GetBinContent(nbin_X+1,nbin_Y+1) - binContent) / binContent) rel_diff = std::fabs(v_2D_EraRuns.at(i)->GetBinContent(nbin_X+1,nbin_Y+1) - binContent) / binContent;
      
      }
      */

      if (binContent < 1e-6) rel_diff = 0;

      vec_var[nbin_X][nbin_Y] += (rel_diff) * (rel_diff) ;

      
      h2D_SF->SetBinError(nbin_X+1, nbin_Y+1, binContent*std::sqrt(vec_var[nbin_X][nbin_Y]) );

      h2D_SF_Errors->SetBinContent(nbin_X+1, nbin_Y+1, binContent*std::sqrt(vec_var[nbin_X][nbin_Y]) );



    }

  }

}


void CalcUncBand_njet_nvtx(TH1F *h_SF, TH1F *h_SF_njets_high, TH1F *h_SF_njets_low, TH1F *h_SF_nvtx_high, TH1F *h_SF_nvtx_low, TH1F *uncBand, TH1F *uncBand_rel_diff_stat, TH1F *uncBand_rel_diff_njets, TH1F *uncBand_rel_diff_nvtx ) {

  int nbins_X = h_SF->GetNbinsX();
  
  std::vector<double> vec_var(nbins_X, 0); 

  for (Int_t nbin_X = 0; nbin_X < nbins_X; nbin_X++) {

      double binContent = h_SF->GetBinContent(nbin_X+1);
      double binError = h_SF->GetBinError(nbin_X+1);

      double rel_diff_stat = binError / binContent;      

      double rel_diff_njets_high;
      double rel_diff_njets_low;
      double rel_diff_nvtx_high;
      double rel_diff_nvtx_low;
      
      if (h_SF_njets_high->GetBinContent(nbin_X+1) < 1e-6) rel_diff_njets_high = 0;
      else rel_diff_njets_high = std::fabs(h_SF_njets_high->GetBinContent(nbin_X+1) - binContent) / binContent;
      if (h_SF_njets_low->GetBinContent(nbin_X+1) < 1e-6) rel_diff_njets_low = 0;
      else rel_diff_njets_low = std::fabs(h_SF_njets_low->GetBinContent(nbin_X+1) - binContent) / binContent;
      if (h_SF_nvtx_high->GetBinContent(nbin_X+1) < 1e-6) rel_diff_nvtx_high = 0;
      else rel_diff_nvtx_high = std::fabs(h_SF_nvtx_high->GetBinContent(nbin_X+1) - binContent) / binContent;
      if (h_SF_nvtx_low->GetBinContent(nbin_X+1) < 1e-6) rel_diff_nvtx_low = 0;
      else rel_diff_nvtx_low = std::fabs(h_SF_nvtx_low->GetBinContent(nbin_X+1) - binContent) / binContent;

      if (binContent < 1e-6) {
	rel_diff_stat = 0;

	rel_diff_njets_high = 0;
	rel_diff_njets_low = 0;
	rel_diff_nvtx_high = 0;
	rel_diff_nvtx_low = 0;
      }


      double rel_diff_njets = std::max(rel_diff_njets_high,rel_diff_njets_low);
      double rel_diff_nvtx = std::max(rel_diff_nvtx_high,rel_diff_nvtx_low);
      //      cout << "rel_diff_njets: " << rel_diff_njets << "rel_diff_nvtx " << rel_diff_nvtx << endl;

      vec_var[nbin_X] += (rel_diff_njets) * (rel_diff_njets);
      vec_var[nbin_X] += (rel_diff_nvtx) * (rel_diff_nvtx);

      uncBand_rel_diff_stat->SetBinContent(nbin_X+1, rel_diff_stat );
      uncBand_rel_diff_njets->SetBinContent(nbin_X+1, rel_diff_njets );
      uncBand_rel_diff_nvtx->SetBinContent(nbin_X+1, rel_diff_nvtx );

      uncBand_rel_diff_stat->SetBinError(nbin_X+1, 0 );
      uncBand_rel_diff_njets->SetBinError(nbin_X+1, 0 );
      uncBand_rel_diff_nvtx->SetBinError(nbin_X+1, 0 );

      uncBand->SetBinError(nbin_X+1,  binContent*std::sqrt(vec_var[nbin_X]) );
      //      cout << "Checking Final Uncertainty: " << uncBand->GetBinError(nbin_X+1) << endl;
  }

}



void CalcUncBand_EraRun(TH1F *h_SF, vector<TH1F*> v_SF_EraRuns, vector<float> v_SF_EraRuns_Scale, TH1F *uncBand, TH1F *uncBand_rel_diff_era, TH1F *uncBand_rel_diff_Full) {

  TH1F *h_SF_lumiweighted = (TH1F*)v_SF_EraRuns[0]->Clone("h_SF_lumiweighted");
  h_SF_lumiweighted->Scale(v_SF_EraRuns_Scale[0]);

  for (UInt_t i = 1; i < v_SF_EraRuns.size(); i++) {	
    v_SF_EraRuns[i]->Scale(v_SF_EraRuns_Scale[i]);
    h_SF_lumiweighted->Add(v_SF_EraRuns[i]);
  }

  int nbins_X = h_SF->GetNbinsX();
  
  std::vector<double> vec_var(nbins_X, 0); 

  for (Int_t nbin_X = 0; nbin_X < nbins_X; nbin_X++) {

    double binContent = h_SF->GetBinContent(nbin_X+1);
    double uncBandError = uncBand->GetBinError(nbin_X+1);

    double binContent_lumiweighted = h_SF_lumiweighted->GetBinContent(nbin_X+1);

    double rel_diff;
    
    if (binContent_lumiweighted < 1e-6) rel_diff = 0;
    else rel_diff = std::fabs(binContent - binContent_lumiweighted) / binContent;

    double uncBandRelDiff = uncBandError/binContent;


    /*
    for (UInt_t i = 0; i < v_SF_EraRuns.size(); i++) {

      if (rel_diff < std::fabs(v_SF_EraRuns.at(i)->GetBinContent(nbin_X+1) - binContent) / binContent) rel_diff = std::fabs(v_SF_EraRuns.at(i)->GetBinContent(nbin_X+1) - binContent) / binContent;
      
    }
    */

    if (binContent < 1e-6) {
      rel_diff = 0;
      uncBandRelDiff = 0;
    }

    vec_var[nbin_X] += (rel_diff) * (rel_diff) + (uncBandRelDiff) * (uncBandRelDiff);


    uncBand_rel_diff_era->SetBinContent(nbin_X+1,  rel_diff );
    uncBand_rel_diff_Full->SetBinContent(nbin_X+1,  std::sqrt(vec_var[nbin_X]) );
    uncBand_rel_diff_era->SetBinError(nbin_X+1,  0 );
    uncBand_rel_diff_Full->SetBinError(nbin_X+1,  0 );

    uncBand->SetBinError(nbin_X+1,  binContent*std::sqrt(vec_var[nbin_X]) );

  }

}




/// Main Function                                                                                                                                                               

void triggereffmakehistograms(string era_){

  TH1::SetDefaultSumw2();

  if (era_ == "2017") {
    gSystem->Exec("rm triggereff/output_hists_Data_All_2017.root");
    gSystem->Exec("hadd triggereff/output_hists_Data_All_2017.root triggereff/output_hists_Data_run2017B_2017.root triggereff/output_hists_Data_run2017C_2017.root triggereff/output_hists_Data_run2017D_2017.root triggereff/output_hists_Data_run2017E_2017.root triggereff/output_hists_Data_run2017F_2017.root");
  }

  if (era_ == "2016preVFP_UL") {
    gSystem->Exec("rm triggereff/output_hists_Data_All_2016preVFP_UL.root");
    gSystem->Exec("hadd triggereff/output_hists_Data_All_2016preVFP_UL.root triggereff/output_hists_Data_run2016_ULB_2016preVFP_UL.root triggereff/output_hists_Data_run2016_ULC_2016preVFP_UL.root triggereff/output_hists_Data_run2016_ULD_2016preVFP_UL.root triggereff/output_hists_Data_run2016_ULE_2016preVFP_UL.root triggereff/output_hists_Data_run2016_ULF1_2016preVFP_UL.root");
  }

  if (era_ == "2016postVFP_UL") {
    gSystem->Exec("rm triggereff/output_hists_Data_All_2016postVFP_UL.root");
    gSystem->Exec("hadd triggereff/output_hists_Data_All_2016postVFP_UL.root triggereff/output_hists_Data_run2016_ULF2_2016postVFP_UL.root triggereff/output_hists_Data_run2016_ULG_2016postVFP_UL.root triggereff/output_hists_Data_run2016_ULH_2016postVFP_UL.root");
  }

  if (era_ == "2017_UL") {
    gSystem->Exec("rm triggereff/output_hists_Data_All_2017_UL.root");
    gSystem->Exec("hadd triggereff/output_hists_Data_All_2017_UL.root triggereff/output_hists_Data_run2017_ULB_2017_UL.root triggereff/output_hists_Data_run2017_ULC_2017_UL.root triggereff/output_hists_Data_run2017_ULD_2017_UL.root triggereff/output_hists_Data_run2017_ULE_2017_UL.root triggereff/output_hists_Data_run2017_ULF_2017_UL.root");
  }

  if (era_ == "2018_UL") {
    gSystem->Exec("rm triggereff/output_hists_Data_All_2018_UL.root");
    gSystem->Exec("hadd triggereff/output_hists_Data_All_2018_UL.root triggereff/output_hists_Data_run2018_ULA_2018_UL.root triggereff/output_hists_Data_run2018_ULB_2018_UL.root triggereff/output_hists_Data_run2018_ULC_2018_UL.root triggereff/output_hists_Data_run2018_ULD_2018_UL.root");
  }

  TFile* fdata = new TFile(("triggereff/output_hists_Data_All_"+era_+".root").c_str(),"READ");

  vector<TFile*> fMC;
  vector<double> fMC_LumiWeight;

  double topxsec = 831.76;

  TFile* fMC_ttbarsignalplustau_fromDilepton = new TFile(("triggereff/output_hists_MC_ttbarsignalplustau_fromDilepton_"+era_+".root").c_str(),"READ");
  fMC.push_back(fMC_ttbarsignalplustau_fromDilepton);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_ttbarsignalplustau_fromDilepton, topxsec * 0.10706, era_));

  /*
  TFile* fMC_ttbarbg_fromLjets = new TFile(("triggereff/output_hists_MC_ttbarbg_fromLjets_"+era_+".root").c_str(),"READ");
  fMC.push_back(fMC_ttbarbg_fromLjets);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_ttbarbg_fromLjets, topxsec * 0.4380, era_));

  TFile* fMC_ttbarbg_fromHadronic = new TFile(("triggereff/output_hists_MC_ttbarbg_fromHadronic_"+era_+".root").c_str(),"READ");
  fMC.push_back(fMC_ttbarbg_fromHadronic);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_ttbarbg_fromHadronic, topxsec * 0.4570, era_));

  TFile* fMC_singleantitop_tw = new TFile(("triggereff/output_hists_MC_singleantitop_tw_"+era_+".root").c_str(),"READ");
  fMC.push_back(fMC_singleantitop_tw);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_singleantitop_tw, 35.85, era_));

  TFile* fMC_singletop_tw = new TFile(("triggereff/output_hists_MC_singletop_tw_"+era_+".root").c_str(),"READ");
  fMC.push_back(fMC_singletop_tw);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_singletop_tw, 35.85, era_));

  TFile* fMC_dy1050 = new TFile(("triggereff/output_hists_MC_dy1050_"+era_+".root").c_str(),"READ");
  fMC.push_back(fMC_dy1050);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_dy1050, 16.523, era_));

  TFile* fMC_dy50inf = new TFile(("triggereff/output_hists_MC_dy50inf_"+era_+".root").c_str(),"READ");
  fMC.push_back(fMC_dy50inf);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_dy50inf, 3.* 2075.14, era_));
  */

  vector<TFile*> fdata_EraRuns;
  vector<string> EraRuns;
  vector<float> EraRuns_Scale;

  if (era_ == "2016") {
    EraRuns = {"B", "C", "D", "E", "F", "G", "H"};
    for (UInt_t i = 0; i < EraRuns.size(); i++) {
      TFile* fdata_tmp = new TFile(("triggereff/output_hists_Data_run"+era_+EraRuns[i]+"_"+era_+".root").c_str(),"READ");
      fdata_EraRuns.push_back(fdata_tmp);
    }
    EraRuns_Scale.push_back(0.160);
    EraRuns_Scale.push_back(0.072);
    EraRuns_Scale.push_back(0.118);
    EraRuns_Scale.push_back(0.112);
    EraRuns_Scale.push_back(0.087);
    EraRuns_Scale.push_back(0.211);
    EraRuns_Scale.push_back(0.241);
  }  else if (era_ == "2017") {
    EraRuns = {"B", "C", "D", "E", "F"};
    for (UInt_t i = 0; i < EraRuns.size(); i++) {
      TFile* fdata_tmp = new TFile(("triggereff/output_hists_Data_run"+era_+EraRuns[i]+"_"+era_+".root").c_str(),"READ");
      fdata_EraRuns.push_back(fdata_tmp);
    }
    EraRuns_Scale.push_back(0.115);
    EraRuns_Scale.push_back(0.232);
    EraRuns_Scale.push_back(0.102);
    EraRuns_Scale.push_back(0.224);
    EraRuns_Scale.push_back(0.326);
  }  else if (era_ == "2018") {
    EraRuns = {"A", "B", "C", "D",};
    for (UInt_t i = 0; i < EraRuns.size(); i++) {
      TFile* fdata_tmp = new TFile(("triggereff/output_hists_Data_run"+era_+EraRuns[i]+"_"+era_+".root").c_str(),"READ");
      fdata_EraRuns.push_back(fdata_tmp);
    }
    EraRuns_Scale.push_back(0.233);
    EraRuns_Scale.push_back(0.118);
    EraRuns_Scale.push_back(0.116);
    EraRuns_Scale.push_back(0.532);
  }  else if (era_ == "2016preVFP_UL") {
    EraRuns = {"B", "C", "D", "E", "F1"};
    for (UInt_t i = 0; i < EraRuns.size(); i++) {
      TFile* fdata_tmp = new TFile(("triggereff/output_hists_Data_run2016_UL"+EraRuns[i]+"_"+era_+".root").c_str(),"READ");
      fdata_EraRuns.push_back(fdata_tmp);
    }
    EraRuns_Scale.push_back(0.296394305);
    EraRuns_Scale.push_back(0.133278918);
    EraRuns_Scale.push_back(0.21792112);
    EraRuns_Scale.push_back(0.20673243);
    EraRuns_Scale.push_back(0.145673226);
  }  else if (era_ == "2016postVFP_UL") {
    EraRuns = {"F2", "G", "H"};
    for (UInt_t i = 0; i < EraRuns.size(); i++) {
      TFile* fdata_tmp = new TFile(("triggereff/output_hists_Data_run2016_UL"+EraRuns[i]+"_"+era_+".root").c_str(),"READ");
      fdata_EraRuns.push_back(fdata_tmp);
    }
    EraRuns_Scale.push_back(0.034416982);
    EraRuns_Scale.push_back(0.450783111);
    EraRuns_Scale.push_back(0.514799907);
  }  else if (era_ == "2017_UL") {
    EraRuns = {"B", "C", "D", "E", "F"};
    for (UInt_t i = 0; i < EraRuns.size(); i++) {
      TFile* fdata_tmp = new TFile(("triggereff/output_hists_Data_run"+era_+EraRuns[i]+"_"+era_+".root").c_str(),"READ");
      fdata_EraRuns.push_back(fdata_tmp);
    }
    EraRuns_Scale.push_back(0.115800593);
    EraRuns_Scale.push_back(0.230812526);
    EraRuns_Scale.push_back(0.102406592);
    EraRuns_Scale.push_back(0.224557685);
    EraRuns_Scale.push_back(0.326422605);
  }  else if (era_ == "2018_UL") {
    EraRuns = {"A", "B", "C", "D",};
    for (UInt_t i = 0; i < EraRuns.size(); i++) {
      TFile* fdata_tmp = new TFile(("triggereff/output_hists_Data_run"+era_+EraRuns[i]+"_"+era_+".root").c_str(),"READ");
      fdata_EraRuns.push_back(fdata_tmp);
    }
    EraRuns_Scale.push_back(0.234448169);
    EraRuns_Scale.push_back(0.11810563);
    EraRuns_Scale.push_back(0.115302214);
    EraRuns_Scale.push_back(0.532143987);
  }


  
  string output_dir = era_;
  system(("mkdir -p triggereff/"+output_dir).c_str());
  system(("mkdir -p triggereff/"+output_dir+"/effData/").c_str());
  system(("mkdir -p triggereff/"+output_dir+"/effMC/").c_str());
  system(("mkdir -p triggereff/"+output_dir+"/alpha/").c_str());
  TFile *output_hists = new TFile(("triggereff/"+output_dir+"/results"+output_dir+".root").c_str(),"RECREATE");

  vector<string> channels = {"ee", "emu", "mumu"};
  vector<string> selections = {"afterSelections", "dilepTriggers", "crossTriggers", "crossAndDilepTriggers"};

  /* sidereal
  vector<string> variables = {"met", "njets", "nvtx", "nvtx_good", "lepApt", "lepBpt", "lepAeta", "lepBeta", "lepAphi", "lepBphi", "lepApt_eeIn", "lepBpt_eeIn", "lepApt_eeSplit", "lepBpt_eeSplit", "lepApt_eeOut", "lepBpt_eeOut", "lepApt_emuIn", "lepBpt_emuIn", "lepApt_emuOut", "lepBpt_emuOut", "sidereal"}; 
  vector<string> xAxisTitles = {"Missing Transverse Energy (GeV)", "number of jets", "number of vertex", "number of good vertex", "leading lepton pT (GeV)", "subleading lepton pT (GeV)", "leading lepton #eta", "subleading lepton #eta", "leading lepton #phi", "subleading lepton #phi", "leading lepton pT (GeV)", "subleading lepton pT (GeV)", "leading lepton pT (GeV)", "subleading lepton pT (GeV)", "leading lepton pT (GeV)", "subleading lepton pT (GeV)", "leading lepton pT (GeV)", "subleading lepton pT (GeV)", "leading lepton pT (GeV)", "subleading lepton pT (GeV)", "sidereal hour"};
  vector<string> xAxisTitles_emu = {"Missing Transverse Energy (GeV)", "number of jets", "number of vertex", "number of good vertex", "Electron pT (GeV)", "Muon pT (GeV)", "Electron #eta", "Muon #eta", "Electron #phi", "Muon #phi", "Electron pT (GeV)", "Muon pT (GeV)", "Electron pT (GeV)", "Muon pT (GeV)", "Electron pT (GeV)", "Muon pT (GeV)", "Electron pT (GeV)", "Muon pT (GeV)", "Electron pT (GeV)", "Muon pT (GeV)", "sidereal hour"};
  vector<double> xMins = {100., 0., 0., 0., 20., 15., -2.5, -2.5, -4., -4., 20., 15., 20., 15., 20., 15., 20., 15., 20., 15., 0.};
  vector<double> xMaxs = {200., 7., 60., 60., 500., 500., 2.5, 2.5, 4., 4., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500., 24.};
  */

  vector<string> variables = {"met", "njets", "nvtx", "nvtx_good", "lepApt", "lepBpt", "lepAeta", "lepBeta", "lepAphi", "lepBphi", "lepApt_eeIn", "lepBpt_eeIn", "lepApt_eeSplit", "lepBpt_eeSplit", "lepApt_eeOut", "lepBpt_eeOut", "lepApt_emuIn", "lepBpt_emuIn", "lepApt_emuOut", "lepBpt_emuOut"}; 
  vector<string> xAxisTitles = {"Missing Transverse Energy (GeV)", "number of jets", "number of vertex", "number of good vertex", "leading lepton pT (GeV)", "subleading lepton pT (GeV)", "leading lepton #eta", "subleading lepton #eta", "leading lepton #phi", "subleading lepton #phi", "leading lepton pT (GeV)", "subleading lepton pT (GeV)", "leading lepton pT (GeV)", "subleading lepton pT (GeV)", "leading lepton pT (GeV)", "subleading lepton pT (GeV)", "leading lepton pT (GeV)", "subleading lepton pT (GeV)", "leading lepton pT (GeV)", "subleading lepton pT (GeV)"};
  vector<string> xAxisTitles_emu = {"Missing Transverse Energy (GeV)", "number of jets", "number of vertex", "number of good vertex", "Electron pT (GeV)", "Muon pT (GeV)", "Electron #eta", "Muon #eta", "Electron #phi", "Muon #phi", "Electron pT (GeV)", "Muon pT (GeV)", "Electron pT (GeV)", "Muon pT (GeV)", "Electron pT (GeV)", "Muon pT (GeV)", "Electron pT (GeV)", "Muon pT (GeV)", "Electron pT (GeV)", "Muon pT (GeV)"};
  vector<double> xMins = {100., 0., 0., 0., 20., 15., -2.5, -2.5, -4., -4., 20., 15., 20., 15., 20., 15., 20., 15., 20., 15.};
  vector<double> xMaxs = {200., 7., 60., 60., 500., 500., 2.5, 2.5, 4., 4., 500., 500., 500., 500., 500., 500., 500., 500., 500., 500.};


  vector<string> multivariables = {"lepABpt", "lepABeta", "lepABeta_onebin", "lepABpt_eeIn", "lepABpt_eeOut", "lepABpt_eeSplit", "lepABpt_emuIn", "lepABpt_emuOut"} ;
  vector<string> multivariablesxAxisTitles = {"leading lepton pT (GeV)", "leading lepton #eta", "leading lepton #eta", "leading lepton pT (GeV)", "leading lepton pT (GeV)", "leading lepton pT (GeV)", "leading lepton pT (GeV)", "leading lepton pT (GeV)"};
  vector<string> multivariablesyAxisTitles = {"subleading lepton pT (GeV)", "subleading lepton #eta", "subleading lepton #eta", "subleading lepton pT (GeV)", "subleading lepton pT (GeV)", "subleading lepton pT (GeV)", "subleading lepton pT (GeV)", "subleading lepton pT (GeV)"};

  vector<string> systematics = {"nominal", "low2_njets", "high2_njets",  "low3_njets", "high3_njets",  "low4_njets", "high4_njets", "low20_nvtx", "high20_nvtx", "low25_nvtx", "high25_nvtx", "low30_nvtx", "high30_nvtx", "low35_nvtx", "high35_nvtx", "low40_nvtx", "high40_nvtx", "low45_nvtx", "high45_nvtx", "low50_nvtx", "high50_nvtx"};


  vector<string> njetsystematics = {"low2_njets", "high2_njets",  "low3_njets", "high3_njets",  "low4_njets", "high4_njets"};
  vector<string> nvtxsystematics = {"low20_nvtx", "high20_nvtx", "low25_nvtx", "high25_nvtx", "low30_nvtx", "high30_nvtx", "low35_nvtx", "high35_nvtx", "low40_nvtx", "high40_nvtx", "low45_nvtx", "high45_nvtx", "low50_nvtx", "high50_nvtx"};


  for (UInt_t i = 0; i < channels.size(); i++) {

    cout << channels[i] << endl;

    for (UInt_t j = 0; j < variables.size(); j++) {

      cout << variables[j] << endl;

      if(channels[i] != "ee") {
	if(variables[j] == "lepApt_eeIn" || variables[j] == "lepBpt_eeIn" || variables[j] == "lepApt_eeOut" || variables[j] == "lepBpt_eeOut" || variables[j] == "lepApt_eeSplit" || variables[j] == "lepBpt_eeSplit" ) continue;
      }

      if(channels[i] != "emu") {
	if(variables[j] == "lepApt_emuIn" || variables[j] == "lepBpt_emuIn" || variables[j] == "lepApt_emuOut" || variables[j] == "lepBpt_emuOut") continue;
      }



      for (UInt_t k = 0; k < systematics.size(); k++) {

	cout << systematics[k] << endl;

	if(systematics[k]== "low2_njets" || systematics[k]== "high2_njets" || systematics[k]== "low4_njets" || systematics[k]== "high4_njets" || systematics[k]== "low5_njets" || systematics[k]== "high5_njets") continue;
	if(systematics[k]== "low50_nvtx" || systematics[k]== "high50_nvtx") continue;

	if(era_ == "2016" || era_ == "2016preVFP_UL" || era_ == "2016postVFP_UL") {
	  if(systematics[k]=="low35_nvtx") continue;
	  if(systematics[k]=="high35_nvtx") continue;
	  if(systematics[k]=="low40_nvtx") continue;
	  if(systematics[k]=="high40_nvtx") continue;
	  if(systematics[k]=="low45_nvtx") continue;
	  if(systematics[k]=="high45_nvtx") continue;
	  if(systematics[k]=="low50_nvtx") continue;
	  if(systematics[k]=="high50_nvtx") continue;	  
	}

	TH1F *h_crossTriggers_data = (TH1F*)fdata->Get((channels[i]+"_"+variables[j]+"_"+systematics[k]+"_"+selections[2]).c_str());  
	TH1F *h_crossAndDilepTriggers_data = (TH1F*)fdata->Get((channels[i]+"_"+variables[j]+"_"+systematics[k]+"_"+selections[3]).c_str()); 
	h_crossTriggers_data->Sumw2();
	h_crossAndDilepTriggers_data->Sumw2(); 

	//	THStack *hs_crossTriggers_MC = new THStack((channels[i]+"_"+variables[j]+"_"+systematics[k]+"_"+selections[2]+"_stack").c_str(),"");
	//	THStack *hs_crossAndDilepTriggers_MC = new THStack((channels[i]+"_"+variables[j]+"_"+systematics[k]+"_"+selections[3]+"_stack").c_str(),"");

	TH1F *h_crossTriggers_MC = (TH1F*)fMC[0]->Get((channels[i]+"_"+variables[j]+"_"+systematics[k]+"_"+selections[2]).c_str());
	h_crossTriggers_MC->Sumw2();
	//	h_crossTriggers_MC->Scale(fMC_LumiWeight[0]);
	//	hs_crossTriggers_MC->Add(h_crossTriggers_MC);


	TH1F *h_crossAndDilepTriggers_MC = (TH1F*)fMC[0]->Get((channels[i]+"_"+variables[j]+"_"+systematics[k]+"_"+selections[3]).c_str());;
	h_crossAndDilepTriggers_MC->Sumw2();
	//	h_crossAndDilepTriggers_MC->Scale(fMC_LumiWeight[0]);;	
	//	hs_crossAndDilepTriggers_MC->Add(h_crossAndDilepTriggers_MC);


	/*
	for (UInt_t m = 1; m < 7; m++) {

	  cout << m << endl;

	  TH1F *h_crossTriggers_MC_tmp = (TH1F*)fMC[m]->Get((channels[i]+"_"+variables[j]+"_"+systematics[k]+"_"+selections[2]).c_str());
	  h_crossTriggers_MC_tmp->Sumw2();
	  h_crossTriggers_MC_tmp->Scale(fMC_LumiWeight[m]);
	  h_crossTriggers_MC->Add(h_crossTriggers_MC_tmp);
	  hs_crossTriggers_MC->Add(h_crossTriggers_MC_tmp);

	  TH1F *h_crossAndDilepTriggers_MC_tmp = (TH1F*)fMC[m]->Get((channels[i]+"_"+variables[j]+"_"+systematics[k]+"_"+selections[3]).c_str());
	  h_crossAndDilepTriggers_MC_tmp->Sumw2();
	  h_crossAndDilepTriggers_MC_tmp->Scale(fMC_LumiWeight[m]);;	
	  h_crossAndDilepTriggers_MC->Add(h_crossAndDilepTriggers_MC_tmp);
	  hs_crossAndDilepTriggers_MC->Add(h_crossAndDilepTriggers_MC_tmp);

	}
	*/


	///TO GET EFF and SF
	TH1F *h_eff_data = (TH1F*)h_crossAndDilepTriggers_data->Clone(("h_eff_data_"+channels[i]+"_"+variables[j]+"_"+systematics[k]).c_str());
	h_eff_data->Divide(h_crossAndDilepTriggers_data,h_crossTriggers_data,1,1,"B");
	TH1F *h_eff_MC = (TH1F*)h_crossAndDilepTriggers_MC->Clone(("h_eff_MC_"+channels[i]+"_"+variables[j]+"_"+systematics[k]).c_str());
	h_eff_MC->Divide(h_crossAndDilepTriggers_MC,h_crossTriggers_MC,1,1,"B");


      	TGraphAsymmErrors *g_eff_data = new TGraphAsymmErrors(h_crossAndDilepTriggers_data, h_crossTriggers_data, "cp");
      	TGraphAsymmErrors *g_eff_MC = new TGraphAsymmErrors(h_crossAndDilepTriggers_MC, h_crossTriggers_MC, "cp");

	TH1F *h_SF = (TH1F*)h_eff_data->Clone(("h_SF_"+channels[i]+"_"+variables[j]+"_"+systematics[k]).c_str()); 
	h_SF->Divide(h_eff_MC);

	TGraphAsymmErrors *g_SF = new TGraphAsymmErrors(h_SF);

	PropagateKatzStatError(h_crossAndDilepTriggers_data,h_crossTriggers_data,h_crossAndDilepTriggers_MC,h_crossTriggers_MC,h_SF,g_SF);


	//	hs_crossAndDilepTriggers_MC->Write();
	//	hs_crossTriggers_MC->Write();
	h_eff_data->Write();
	h_eff_MC->Write();
	h_SF->Write();

	//	if( k == 0 ) {
	  TCanvas *c_eff = new TCanvas(("c_eff_"+channels[i]+"_"+variables[j]+"_"+systematics[k]).c_str(), "", 1000., 1000.);
	  TPad *p_eff = new TPad(("p_eff_"+channels[i]+"_"+variables[j]+"_"+systematics[k]).c_str(), "", 0, 0.15, 1, 1.0); 
	  setPad(p_eff);
	  drawGraph(h_eff_data, h_eff_MC, g_eff_data, g_eff_MC, "Efficiency", channels[i], xMins[j], xMaxs[j], era_);
	  c_eff->cd();
	  TPad *p_ratio_eff = new TPad(("p_ratio_eff_"+channels[i]+"_"+variables[j]+"_"+systematics[k]).c_str(), "", 0, 0.05, 1, 0.3);  
	  setPad(p_ratio_eff, false);
	  drawRatio(h_SF, g_SF, channels[i], xAxisTitles[j], xMins[j], xMaxs[j], "Scale Factor", 0.95, 1.05);
	  system(("mkdir -p triggereff/"+output_dir+"/"+systematics[k]+"/").c_str());

	  c_eff->SaveAs(("triggereff/"+output_dir+"/"+systematics[k]+"/g_"+variables[j]+"_"+channels[i]+".pdf").c_str());  
	  c_eff->SaveAs(("triggereff/"+output_dir+"/"+systematics[k]+"/g_"+variables[j]+"_"+channels[i]+".C").c_str()); 
	  
	  c_eff->Write();

	  TCanvas *c_eff_data = new TCanvas(("c_eff_data_"+channels[i]+"_"+variables[j]+"_"+systematics[k]).c_str(), "", 1000., 1000.);
	  TPad *p_eff_data = new TPad(("p_eff_data_"+channels[i]+"_"+variables[j]+"_"+systematics[k]).c_str(), "", 0, 0.15, 1, 1.0); 
	  setPad(p_eff_data);
	  drawGraph(h_crossAndDilepTriggers_data,h_crossTriggers_data,"Entries",channels[i], xMins[j], xMaxs[j], era_);
	  c_eff_data->cd();
	  TPad *p_ratio_eff_data = new TPad(("p_ratio_eff_data_"+channels[i]+"_"+variables[j]+"_"+systematics[k]).c_str(), "", 0, 0.05, 1, 0.3);  
	  setPad(p_ratio_eff_data, false);
	  drawRatio(h_eff_data, g_eff_data, channels[i], xAxisTitles[j], xMins[j], xMaxs[j], "Data Efficiency", 0.95, 1.05);
	  system(("mkdir -p triggereff/"+output_dir+"/effData/"+systematics[k]+"/").c_str());

	  c_eff_data->SaveAs(("triggereff/"+output_dir+"/effData/"+systematics[k]+"/"+"g_"+variables[j]+"_"+channels[i]+"_data.pdf").c_str());      
	  c_eff_data->SaveAs(("triggereff/"+output_dir+"/effData/"+systematics[k]+"/"+"g_"+variables[j]+"_"+channels[i]+"_data.C").c_str());      

	  c_eff_data->Write();
	  
	  TCanvas *c_eff_MC = new TCanvas(("c_eff_MC_"+channels[i]+"_"+variables[j]+"_"+systematics[k]).c_str(), "", 1000., 1000.);
	  TPad *p_eff_MC = new TPad(("p_eff_MC_"+channels[i]+"_"+variables[j]+"_"+systematics[k]).c_str(), "", 0, 0.15, 1, 1.0); 
	  setPad(p_eff_MC);
	  drawGraph(h_crossAndDilepTriggers_MC,h_crossTriggers_MC,"Entries",channels[i], xMins[j], xMaxs[j], era_);
	  c_eff_MC->cd();
	  TPad *p_ratio_eff_MC = new TPad(("p_ratio_eff_MC_"+channels[i]+"_"+variables[j]+"_"+systematics[k]).c_str(), "", 0, 0.05, 1, 0.3);  
	  setPad(p_ratio_eff_MC, false);
	  drawRatio(h_eff_MC, g_eff_MC, channels[i], xAxisTitles[j], xMins[j], xMaxs[j], "MC Efficiency", 0.95, 1.05);
	  system(("mkdir -p triggereff/"+output_dir+"/effMC/"+systematics[k]+"/").c_str());

	  c_eff_MC->SaveAs(("triggereff/"+output_dir+"/effMC/"+systematics[k]+"/"+"g_"+variables[j]+"_"+channels[i]+"_MC.pdf").c_str());      
	  c_eff_MC->SaveAs(("triggereff/"+output_dir+"/effMC/"+systematics[k]+"/"+"g_"+variables[j]+"_"+channels[i]+"_MC.C").c_str());      

	  c_eff_MC->Write();
	  //	}


	/// 


	///TO GET ALPHA 
	if( k == 0 ) {

	  TH1F *h_afterSelections_MC = (TH1F*)fMC[0]->Get((channels[i]+"_"+variables[j]+"_"+systematics[k]+"_"+selections[0]).c_str());      
	  TH1F *h_dilepTriggers_MC = (TH1F*)fMC[0]->Get((channels[i]+"_"+variables[j]+"_"+systematics[k]+"_"+selections[1]).c_str());
      
	  TH1F *h_eff_dilepTriggers_MC = (TH1F*)h_dilepTriggers_MC->Clone(("h_eff_dilepTriggers_MC_"+channels[i]+"_"+variables[j]+"_"+systematics[k]).c_str());
	  h_eff_dilepTriggers_MC->Divide(h_dilepTriggers_MC,h_afterSelections_MC,1,1,"B");
	  TH1F *h_eff_crossTriggers_MC = (TH1F*)h_crossTriggers_MC->Clone(("h_eff_crossTriggers_MC_"+channels[i]+"_"+variables[j]+"_"+systematics[k]).c_str());
	  h_eff_crossTriggers_MC->Divide(h_crossTriggers_MC,h_afterSelections_MC,1,1,"B");
	  TH1F *h_eff_crossAndDilepTriggers_MC = (TH1F*)h_crossAndDilepTriggers_MC->Clone(("h_eff_crossAndDilepTriggers_MC_"+channels[i]+"_"+variables[j]+"_"+systematics[k]).c_str());
	  h_eff_crossAndDilepTriggers_MC->Divide(h_crossAndDilepTriggers_MC,h_afterSelections_MC,1,1,"B");
	  TH1F *h_crossAndDilepTriggers_overCrossTriggers_MC = (TH1F*)h_crossAndDilepTriggers_MC->Clone(("h_eff_crossAndDilepTriggers_MC_"+channels[i]+"_"+variables[j]).c_str());
	  h_crossAndDilepTriggers_overCrossTriggers_MC->Divide(h_crossAndDilepTriggers_MC,h_crossTriggers_MC,1,1,"B");

	  TH1F *h_alpha = (TH1F*)h_eff_dilepTriggers_MC->Clone(("h_alpha_"+channels[i]+"_"+variables[j]).c_str());
	  h_alpha->Divide(h_crossAndDilepTriggers_overCrossTriggers_MC);
	  
	  h_eff_dilepTriggers_MC->Write();
	  h_eff_crossTriggers_MC->Write();
	  h_eff_crossAndDilepTriggers_MC->Write();
	  h_crossAndDilepTriggers_overCrossTriggers_MC->Write();
	  h_alpha->Write();

	  TCanvas *c_alpha = new TCanvas(("c_alpha_"+channels[i]+"_"+variables[j]).c_str(), "", 1000., 1000.);
	  TPad *p_alpha = new TPad(("p_alpha_"+channels[i]+"_"+variables[j]).c_str(), "", 0, 0.15, 1, 1.0); 
	  setPad(p_alpha);
	  drawGraph_alpha(h_eff_crossTriggers_MC, h_eff_dilepTriggers_MC, h_eff_crossAndDilepTriggers_MC, channels[i], xMins[j], xMaxs[j], era_);
	  c_alpha->cd();
	  TPad *p_ratio_alpha = new TPad(("p_ratio_alpha_"+channels[i]+"_"+variables[j]).c_str(), "", 0, 0.05, 1, 0.3);  
	  setPad(p_ratio_alpha, false);
	  drawRatio(h_alpha, channels[i], xAxisTitles[j], xMins[j], xMaxs[j], "Alpha", 0.98, 1.02);
	  system(("mkdir -p triggereff/"+output_dir+"/alpha/"+systematics[k]+"/").c_str());
	  c_alpha->SaveAs(("triggereff/"+output_dir+"/alpha/"+systematics[k]+"/"+"g_"+variables[j]+"_"+channels[i]+"_alpha.pdf").c_str());
	  //	  c_alpha->SaveAs(("triggereff/"+output_dir+"/alpha/"+systematics[k]+"/"+"g_"+variables[j]+"_"+channels[i]+"_alpha.C").c_str());
	  c_alpha->Write();
	}
	///

      }

    }

    ///TO GET 2D SF
    for (UInt_t j = 0; j < multivariables.size(); j++) {

      for (UInt_t k = 0; k < systematics.size(); k++) {
	
	TH2F *h2D_crossTriggers_data = (TH2F*)fdata->Get((channels[i]+"_"+multivariables[j]+"_"+systematics[k]+"_"+selections[2]).c_str());  
	TH2F *h2D_crossAndDilepTriggers_data = (TH2F*)fdata->Get((channels[i]+"_"+multivariables[j]+"_"+systematics[k]+"_"+selections[3]).c_str()); 
	h2D_crossTriggers_data->Sumw2();
	h2D_crossAndDilepTriggers_data->Sumw2(); 

	TH2F *h2D_crossTriggers_MC = (TH2F*)fMC[0]->Get((channels[i]+"_"+multivariables[j]+"_"+systematics[k]+"_"+selections[2]).c_str());
	TH2F *h2D_crossAndDilepTriggers_MC = (TH2F*)fMC[0]->Get((channels[i]+"_"+multivariables[j]+"_"+systematics[k]+"_"+selections[3]).c_str());


	/*
	if(multivariables[j] == "lepABpt_eeIn"  || multivariables[j] == "lepABpt_eeOut" ||  multivariables[j] == "lepABpt_eeSplit" ||  multivariables[j] == "lepABpt_emuIn" ||  multivariables[j] == "lepABpt_emuOut" ) {
	  h2D_crossTriggers_data->Rebin2D(2,2);
	  h2D_crossAndDilepTriggers_data->Rebin2D(2,2);
	  h2D_crossTriggers_MC->Rebin2D(2,2);
	  h2D_crossAndDilepTriggers_MC->Rebin2D(2,2);
	}
	else if(multivariables[j] == "lepABpt_eeOut" ||  multivariables[j] == "lepABpt_eeSplit") {
	  h2D_crossTriggers_data->Rebin2D(2,2);
	  h2D_crossAndDilepTriggers_data->Rebin2D(2,2);
	  h2D_crossTriggers_MC->Rebin2D(2,2);
	  h2D_crossAndDilepTriggers_MC->Rebin2D(2,2);
	}
	*/


	TH2F *h2D_eff_data = (TH2F*)h2D_crossAndDilepTriggers_data->Clone(("h2D_eff_data_"+channels[i]+"_"+multivariables[j]+"_"+systematics[k]).c_str());
	h2D_eff_data->Divide(h2D_crossAndDilepTriggers_data,h2D_crossTriggers_data,1,1,"B");
	//	h2D_eff_data->SaveAs(("triggereff/"+output_dir+"/effData/"+systematics[k]+"/h2D_"+multivariables[j]+"_"+systematics[k]+"_"+channels[i]+".C").c_str());

	TH2F *h2D_eff_MC = (TH2F*)h2D_crossAndDilepTriggers_MC->Clone(("h2D_eff_MC_"+channels[i]+"_"+multivariables[j]+"_"+systematics[k]).c_str());
	h2D_eff_MC->Divide(h2D_crossAndDilepTriggers_MC,h2D_crossTriggers_MC,1,1,"B");
	//	h2D_eff_MC->SaveAs(("triggereff/"+output_dir+"/effMC/"+systematics[k]+"/h2D_"+multivariables[j]+"_"+systematics[k]+"_"+channels[i]+".C").c_str());
	
	TH2F *h2D_SF = (TH2F*)h2D_eff_data->Clone(("h2D_SF_"+channels[i]+"_"+multivariables[j]+"_"+systematics[k]).c_str());
	h2D_SF->Divide(h2D_eff_MC);


	PropagateKatzStatError2D(h2D_crossAndDilepTriggers_data,h2D_crossTriggers_data,h2D_crossAndDilepTriggers_MC,h2D_crossTriggers_MC,h2D_SF);

	TCanvas *c_2D = new TCanvas(("c_2D_"+channels[i]+"_"+multivariables[j]+"_"+systematics[k]).c_str(), "", 1200., 800.);
	TPad *p_2D = new TPad(("p_2D_"+channels[i]+"_"+multivariables[j]+"_"+systematics[k]).c_str(), "", 0, 0, 1, 1); 
	set2DPad(p_2D);
	drawMultivariables(h2D_SF, channels[i], multivariablesxAxisTitles[j], multivariablesyAxisTitles[j], era_);
	h2D_SF->Write();
	if(systematics[k] == "nominal") {
	  c_2D->SaveAs(("triggereff/"+output_dir+"/"+systematics[k]+"/h2D_"+multivariables[j]+"_"+channels[i]+".pdf").c_str());
	  //	  c_2D->SaveAs(("triggereff/"+output_dir+"/"+systematics[k]+"/h2D_"+multivariables[j]+"_"+channels[i]+".C").c_str());
	}
	
	if(era_ == "2016" || era_ == "2016preVFP_UL" || era_ == "2016postVFP_UL") {
	  if(systematics[k]=="low20_nvtx" || systematics[k]=="low3_njets" || systematics[k]=="high20_nvtx" || systematics[k]=="high3_njets") c_2D->SaveAs(("triggereff/"+output_dir+"/"+systematics[k]+"/h2D_"+multivariables[j]+"_"+channels[i]+"_"+systematics[k]+".pdf").c_str());
	}

	if(era_ == "2017" || era_ == "2017_UL") {
	  if(systematics[k]=="low30_nvtx" || systematics[k]=="low3_njets" || systematics[k]=="high30_nvtx" || systematics[k]=="high3_njets") c_2D->SaveAs(("triggereff/"+output_dir+"/"+systematics[k]+"/h2D_"+multivariables[j]+"_"+channels[i]+"_"+systematics[k]+".pdf").c_str());
	}

	if(era_ == "2018" || era_ == "2018_UL") {
	  if(systematics[k]=="low30_nvtx" || systematics[k]=="low3_njets" || systematics[k]=="high30_nvtx" || systematics[k]=="high3_njets") c_2D->SaveAs(("triggereff/"+output_dir+"/"+systematics[k]+"/h2D_"+multivariables[j]+"_"+channels[i]+"_"+systematics[k]+".pdf").c_str());
	}

	c_2D->Write();


	TCanvas *c_2D_data = new TCanvas(("c_2D_data_"+channels[i]+"_"+multivariables[j]+"_"+systematics[k]).c_str(), "", 1200., 800.);
	TPad *p_2D_data = new TPad(("p_2D_data_"+channels[i]+"_"+multivariables[j]+"_"+systematics[k]).c_str(), "", 0, 0, 1, 1); 
	set2DPad(p_2D_data);
	drawMultivariables(h2D_eff_data, channels[i], multivariablesxAxisTitles[j], multivariablesyAxisTitles[j], era_);
	h2D_eff_data->Write();
	if(systematics[k] == "nominal") {
	  c_2D_data->SaveAs(("triggereff/"+output_dir+"/effData/"+systematics[k]+"/h2D_"+multivariables[j]+"_"+channels[i]+"_data.pdf").c_str());
	}
	
	if(era_ == "2016" || era_ == "2016preVFP_UL" || era_ == "2016postVFP_UL") {
	  if(systematics[k]=="low20_nvtx" || systematics[k]=="low3_njets" || systematics[k]=="high20_nvtx" || systematics[k]=="high3_njets") c_2D_data->SaveAs(("triggereff/"+output_dir+"/effData/"+systematics[k]+"/h2D_"+multivariables[j]+"_"+channels[i]+"_"+systematics[k]+"_data.pdf").c_str());
	}

	if(era_ == "2017" || era_ == "2017_UL") {
	  if(systematics[k]=="low30_nvtx" || systematics[k]=="low3_njets" || systematics[k]=="high30_nvtx" || systematics[k]=="high3_njets") c_2D_data->SaveAs(("triggereff/"+output_dir+"/effData/"+systematics[k]+"/h2D_"+multivariables[j]+"_"+channels[i]+"_"+systematics[k]+"_data.pdf").c_str());
	}

	if(era_ == "2018" || era_ == "2018_UL") {
	  if(systematics[k]=="low30_nvtx" || systematics[k]=="low3_njets" || systematics[k]=="high30_nvtx" || systematics[k]=="high3_njets") c_2D_data->SaveAs(("triggereff/"+output_dir+"/effData/"+systematics[k]+"/h2D_"+multivariables[j]+"_"+channels[i]+"_"+systematics[k]+"_data.pdf").c_str());
	}

	c_2D_data->Write();



	TCanvas *c_2D_MC = new TCanvas(("c_2D_MC_"+channels[i]+"_"+multivariables[j]+"_"+systematics[k]).c_str(), "", 1200., 800.);
	TPad *p_2D_MC = new TPad(("p_2D_MC_"+channels[i]+"_"+multivariables[j]+"_"+systematics[k]).c_str(), "", 0, 0, 1, 1); 
	set2DPad(p_2D_MC);
	drawMultivariables(h2D_eff_MC, channels[i], multivariablesxAxisTitles[j], multivariablesyAxisTitles[j], era_);
	h2D_eff_MC->Write();
	//	if(systematics[k] == "nominal") {
	  c_2D_MC->SaveAs(("triggereff/"+output_dir+"/effMC/"+systematics[k]+"/h2D_"+multivariables[j]+"_"+channels[i]+"_MC.pdf").c_str());
	  //	}
	
	if(era_ == "2016" || era_ == "2016preVFP_UL" || era_ == "2016postVFP_UL") {
	  if(systematics[k]=="low20_nvtx" || systematics[k]=="low3_njets" || systematics[k]=="high20_nvtx" || systematics[k]=="high3_njets") c_2D_MC->SaveAs(("triggereff/"+output_dir+"/effMC/"+systematics[k]+"/h2D_"+multivariables[j]+"_"+channels[i]+"_"+systematics[k]+"_MC.pdf").c_str());
	}

	if(era_ == "2017" || era_ == "2017_UL") {
	  if(systematics[k]=="low30_nvtx" || systematics[k]=="low3_njets" || systematics[k]=="high30_nvtx" || systematics[k]=="high3_njets") c_2D_MC->SaveAs(("triggereff/"+output_dir+"/effMC/"+systematics[k]+"/h2D_"+multivariables[j]+"_"+channels[i]+"_"+systematics[k]+"_MC.pdf").c_str());
	}

	if(era_ == "2018" || era_ == "2018_UL") {
	  if(systematics[k]=="low30_nvtx" || systematics[k]=="low3_njets" || systematics[k]=="high30_nvtx" || systematics[k]=="high3_njets") c_2D_MC->SaveAs(("triggereff/"+output_dir+"/effMC/"+systematics[k]+"/h2D_"+multivariables[j]+"_"+channels[i]+"_"+systematics[k]+"_MC.pdf").c_str());
	}

	c_2D_MC->Write();
      }

    }
    ///

    ///To Get Era Run Dependence

    cout << "fdata_EraRuns.size " << fdata_EraRuns.size() << endl;

    for (UInt_t j = 0; j < fdata_EraRuns.size(); j++) {

      cout << "EraRun "   << era_  <<  EraRuns[j] << endl;

      for (UInt_t k = 0; k < variables.size(); k++) {

	TH1F *h_crossTriggers_data = (TH1F*)fdata_EraRuns.at(j)->Get((channels[i]+"_"+variables[k]+"_"+systematics[0]+"_"+selections[2]).c_str());  
	TH1F *h_crossAndDilepTriggers_data = (TH1F*)fdata_EraRuns.at(j)->Get((channels[i]+"_"+variables[k]+"_"+systematics[0]+"_"+selections[3]).c_str()); 
	h_crossTriggers_data->Sumw2();
	h_crossAndDilepTriggers_data->Sumw2(); 
	
	TH1F *h_crossTriggers_MC = (TH1F*)fMC[0]->Get((channels[i]+"_"+variables[k]+"_"+systematics[0]+"_"+selections[2]).c_str());
	TH1F *h_crossAndDilepTriggers_MC = (TH1F*)fMC[0]->Get((channels[i]+"_"+variables[k]+"_"+systematics[0]+"_"+selections[3]).c_str());


	TH1F *h_eff_data = (TH1F*)h_crossAndDilepTriggers_data->Clone(("h_eff_data_"+channels[i]+"_"+variables[k]+"_"+systematics[0]+"_"+EraRuns[j]).c_str());
	h_eff_data->Divide(h_crossAndDilepTriggers_data,h_crossTriggers_data,1,1,"B");
	TH1F *h_eff_MC = (TH1F*)h_crossAndDilepTriggers_MC->Clone(("h_eff_MC_"+channels[i]+"_"+variables[k]+"_"+systematics[0]).c_str());
	h_eff_MC->Divide(h_crossAndDilepTriggers_MC,h_crossTriggers_MC,1,1,"B");

	TH1F *h_SF = (TH1F*)h_eff_data->Clone(("h_SF_"+channels[i]+"_"+variables[k]+"_"+systematics[0]+"_"+EraRuns[j]).c_str()); 
	h_SF->Divide(h_eff_MC);

	h_SF->Write();


      }

      for (UInt_t k = 0; k < multivariables.size(); k++) {

	TH2F *h2D_crossTriggers_data = (TH2F*)fdata_EraRuns.at(j)->Get((channels[i]+"_"+multivariables[k]+"_"+systematics[0]+"_"+selections[2]).c_str());  
	TH2F *h2D_crossAndDilepTriggers_data = (TH2F*)fdata_EraRuns.at(j)->Get((channels[i]+"_"+multivariables[k]+"_"+systematics[0]+"_"+selections[3]).c_str()); 
	h2D_crossTriggers_data->Sumw2();
	h2D_crossAndDilepTriggers_data->Sumw2(); 

	TH2F *h2D_crossTriggers_MC = (TH2F*)fMC[0]->Get((channels[i]+"_"+multivariables[k]+"_"+systematics[0]+"_"+selections[2]).c_str());
	TH2F *h2D_crossAndDilepTriggers_MC = (TH2F*)fMC[0]->Get((channels[i]+"_"+multivariables[k]+"_"+systematics[0]+"_"+selections[3]).c_str());


	TH2F *h2D_eff_data = (TH2F*)h2D_crossAndDilepTriggers_data->Clone(("h2D_eff_data_"+channels[i]+"_"+multivariables[k]+"_"+systematics[0]+"_"+EraRuns[j]).c_str());
	h2D_eff_data->Divide(h2D_crossAndDilepTriggers_data,h2D_crossTriggers_data,1,1,"B");

	TH2F *h2D_eff_MC = (TH2F*)h2D_crossAndDilepTriggers_MC->Clone(("h2D_eff_MC_"+channels[i]+"_"+multivariables[k]+"_"+systematics[0]).c_str());
	h2D_eff_MC->Divide(h2D_crossAndDilepTriggers_MC,h2D_crossTriggers_MC,1,1,"B");
	
	/*
	if(multivariables[k] == "lepABpt_eeIn"  || multivariables[k] == "lepABpt_eeOut" ||  multivariables[k] == "lepABpt_eeSplit" ||  multivariables[k] == "lepABpt_emuIn" ||  multivariables[k] == "lepABpt_emuOut" ) {
	  h2D_eff_data->Rebin2D(2,2);
	  h2D_eff_MC->Rebin2D(2,2);
	}
	*/

	TH2F *h2D_SF = (TH2F*)h2D_eff_data->Clone(("h2D_SF_"+channels[i]+"_"+multivariables[k]+"_"+systematics[0]+"_"+EraRuns[j]).c_str());
	h2D_SF->Divide(h2D_eff_MC);

	PropagateKatzStatError2D(h2D_crossAndDilepTriggers_data,h2D_crossTriggers_data,h2D_crossAndDilepTriggers_MC,h2D_crossTriggers_MC,h2D_SF);

	h2D_SF->Write();

      }

    }

    ///



  }
  
  output_hists->Close();




  ///

  TFile fSF(("triggereff/"+era_+"/results"+era_+".root").c_str(),"READ");

  TFile *intermediate_results = new TFile(("triggereff/"+output_dir+"/intermediate_results"+output_dir+".root").c_str(),"RECREATE");

  for (UInt_t i = 0; i < channels.size(); i++) {

    for (UInt_t j = 0; j < multivariables.size(); j++) {

      TH2F *h2D_SF_nominal = (TH2F*)fSF.Get(("h2D_SF_"+channels[i]+"_"+multivariables[j]+"_"+systematics[0]).c_str());

      /// To add nvtx/njet systematic uncertainties to statistical uncertainties for 2D histograms

      for (UInt_t k = 0; k < (nvtxsystematics.size())/2; k++) {


	if(nvtxsystematics[2*k]== "low50_nvtx" || nvtxsystematics[2*k+1]== "high50_nvtx") continue;

	if(era_ == "2016" || era_ == "2016preVFP_UL" || era_ == "2016postVFP_UL") {	
	  if(nvtxsystematics[2*k]=="low35_nvtx") continue;
	  if(nvtxsystematics[2*k+1]=="high35_nvtx") continue;
	  if(nvtxsystematics[2*k]=="low40_nvtx") continue;
	  if(nvtxsystematics[2*k+1]=="high40_nvtx") continue;
	  if(nvtxsystematics[2*k]=="low45_nvtx") continue;
	  if(nvtxsystematics[2*k+1]=="high45_nvtx") continue;
	  if(nvtxsystematics[2*k]=="low50_nvtx") continue;
	  if(nvtxsystematics[2*k+1]=="high50_nvtx") continue;
	}

	for (UInt_t l = 0; l < (njetsystematics.size())/2; l++) {  

	if(njetsystematics[2*l]== "low2_njets" || njetsystematics[2*l+1]== "high2_njets" || njetsystematics[2*l]== "low4_njets" || njetsystematics[2*l+1]== "high4_njets" || njetsystematics[2*l]== "low5_njets" || njetsystematics[2*l+1]== "high5_njets") continue;


	  TH2F *h2D_SF = (TH2F*)h2D_SF_nominal->Clone(("h2D_SF_"+channels[i]+"_"+multivariables[j]+"_njet_nvtx").c_str());
	  TH2F *h2D_SF_njets_high = (TH2F*)fSF.Get(("h2D_SF_"+channels[i]+"_"+multivariables[j]+"_"+njetsystematics[2*l+1]).c_str());
	  TH2F *h2D_SF_njets_low = (TH2F*)fSF.Get(("h2D_SF_"+channels[i]+"_"+multivariables[j]+"_"+njetsystematics[2*l]).c_str());
	  TH2F *h2D_SF_nvtx_high = (TH2F*)fSF.Get(("h2D_SF_"+channels[i]+"_"+multivariables[j]+"_"+nvtxsystematics[2*k+1]).c_str());
	  TH2F *h2D_SF_nvtx_low = (TH2F*)fSF.Get(("h2D_SF_"+channels[i]+"_"+multivariables[j]+"_"+nvtxsystematics[2*k]).c_str());
       
	  AddSystError_njet_nvtx(h2D_SF, h2D_SF_njets_high, h2D_SF_njets_low, h2D_SF_nvtx_high, h2D_SF_nvtx_low);
	  //	  h2D_SF->SaveAs((era_+"/h2D_SF_"+channels[i]+"_"+multivariables[j]+"_"+njetsystematics[2*l+1]+"_"+nvtxsystematics[2*k+1]+"_FullSystError.C").c_str());

	  if(era_ == "2016" || era_ == "2016preVFP_UL" || era_ == "2016postVFP_UL") {
	    if(nvtxsystematics[2*k]=="low20_nvtx" && njetsystematics[2*l]=="low3_njets") h2D_SF->Write();
	  }

	  if(era_ == "2017" || era_ == "2017_UL") {
	    if(nvtxsystematics[2*k]=="low30_nvtx" && njetsystematics[2*l]=="low3_njets") h2D_SF->Write();
	  }

	  if(era_ == "2018" || era_ == "2018_UL") {
	    if(nvtxsystematics[2*k]=="low30_nvtx" && njetsystematics[2*l]=="low3_njets") h2D_SF->Write();
	  }
	   
	}

      }

      ///

      /// To additionally add EraRun systematic uncertainty to the bin errors and make bin error plots for 2D histograms

      TH2F *h2D_SF_njet_nvtx = (TH2F*)intermediate_results->Get(("h2D_SF_"+channels[i]+"_"+multivariables[j]+"_njet_nvtx").c_str());
      TH2F *h2D_SF = (TH2F*)h2D_SF_njet_nvtx->Clone(("h2D_SF_"+channels[i]+"_"+multivariables[j]+"_FullError").c_str());
      TH2F *h2D_SF_Errors = (TH2F*)h2D_SF_njet_nvtx->Clone(("h2D_SF_"+channels[i]+"_"+multivariables[j]+"_BinErrors").c_str());

      vector<TH2F*> v_2D_EraRuns;

      for (UInt_t k = 0; k < EraRuns.size(); k++) {

	cout << "EraRun "   << era_  <<  EraRuns[k] << endl;
       
	TH2F *h_2D_EraRun = (TH2F*)fSF.Get(("h2D_SF_"+channels[i]+"_"+multivariables[j]+"_"+systematics[0]+"_"+EraRuns[k]).c_str());
	v_2D_EraRuns.push_back(h_2D_EraRun);

      }

      AddSystError_EraRun(h2D_SF,v_2D_EraRuns,EraRuns_Scale,h2D_SF_Errors);
      
      TCanvas *c_2D = new TCanvas(("c_2D_"+channels[i]+"_"+multivariables[j]+"_BinErrors").c_str(), "", 1200., 800.);
      TPad *p_2D = new TPad(("p_2D_"+channels[i]+"_"+multivariables[j]+"_BinErrors").c_str(), "", 0, 0, 1, 1); 
      set2DPad(p_2D);
      drawMultivariables(h2D_SF_Errors, channels[i], multivariablesxAxisTitles[j], multivariablesyAxisTitles[j], era_);
      h2D_SF->Write();
      c_2D->SaveAs(("triggereff/"+output_dir+"/h2D_"+multivariables[j]+"_"+channels[i]+"_BinErrors.pdf").c_str());
      c_2D->SaveAs(("triggereff/"+output_dir+"/h2D_"+multivariables[j]+"_"+channels[i]+"_BinErrors.C").c_str());
      
      h2D_SF->Write();
      h2D_SF_Errors->Write();

    ///

    }


    for (UInt_t j = 0; j < variables.size(); j++) {

      if(channels[i] != "ee") {
	if(variables[j] == "lepApt_eeIn" || variables[j] == "lepBpt_eeIn" || variables[j] == "lepApt_eeOut" || variables[j] == "lepBpt_eeOut" || variables[j] == "lepApt_eeSplit" || variables[j] == "lepBpt_eeSplit" ) continue;
      }

      if(channels[i] != "emu") {
	if(variables[j] == "lepApt_emuIn" || variables[j] == "lepBpt_emuIn" || variables[j] == "lepApt_emuOut" || variables[j] == "lepBpt_emuOut") continue;
      }



      TH1F *h_SF_nominal = (TH1F*)fSF.Get(("h_SF_"+channels[i]+"_"+variables[j]+"_"+systematics[0]).c_str());	

      /// To add nvtx/njet systematic uncertainty band to 1D histograms

      for (UInt_t k = 0; k < (nvtxsystematics.size())/2; k++) {

	if(nvtxsystematics[2*k]== "low50_nvtx" || nvtxsystematics[2*k+1]== "high50_nvtx") continue;

	if(era_ == "2016" || era_ == "2016preVFP_UL" || era_ == "2016postVFP_UL") {	
	  if(nvtxsystematics[2*k]=="low35_nvtx") continue;
	  if(nvtxsystematics[2*k+1]=="high35_nvtx") continue;
	  if(nvtxsystematics[2*k]=="low40_nvtx") continue;
	  if(nvtxsystematics[2*k+1]=="high40_nvtx") continue;
	  if(nvtxsystematics[2*k]=="low45_nvtx") continue;
	  if(nvtxsystematics[2*k+1]=="high45_nvtx") continue;
	  if(nvtxsystematics[2*k]=="low50_nvtx") continue;
	  if(nvtxsystematics[2*k+1]=="high50_nvtx") continue;
	}

	for (UInt_t l = 0; l < (njetsystematics.size())/2; l++) {  

	if(njetsystematics[2*l]== "low2_njets" || njetsystematics[2*l+1]== "high2_njets" || njetsystematics[2*l]== "low4_njets" || njetsystematics[2*l+1]== "high4_njets" || njetsystematics[2*l]== "low5_njets" || njetsystematics[2*l+1]== "high5_njets") continue;

	  TH1F *h_crossTriggers_data = (TH1F*)fdata->Get((channels[i]+"_"+variables[j]+"_"+systematics[0]+"_"+selections[2]).c_str());  
	  TH1F *h_crossAndDilepTriggers_data = (TH1F*)fdata->Get((channels[i]+"_"+variables[j]+"_"+systematics[0]+"_"+selections[3]).c_str()); 
	  //	  h_crossTriggers_data->Sumw2();
	  //	  h_crossAndDilepTriggers_data->Sumw2(); 
	  
	  TH1F *h_crossTriggers_MC = (TH1F*)fMC[0]->Get((channels[i]+"_"+variables[j]+"_"+systematics[0]+"_"+selections[2]).c_str());
	  TH1F *h_crossAndDilepTriggers_MC = (TH1F*)fMC[0]->Get((channels[i]+"_"+variables[j]+"_"+systematics[0]+"_"+selections[3]).c_str());

	  TH1F *h_eff_data = (TH1F*)h_crossAndDilepTriggers_data->Clone(("h_eff_data_"+channels[i]+"_"+variables[j]+"_"+systematics[0]).c_str());
	  h_eff_data->Divide(h_crossAndDilepTriggers_data,h_crossTriggers_data,1,1,"B");
	  TH1F *h_eff_MC = (TH1F*)h_crossAndDilepTriggers_MC->Clone(("h_eff_MC_"+channels[i]+"_"+variables[j]+"_"+systematics[0]).c_str());
	  h_eff_MC->Divide(h_crossAndDilepTriggers_MC,h_crossTriggers_MC,1,1,"B");


	  TGraphAsymmErrors *g_eff_data = new TGraphAsymmErrors(h_crossAndDilepTriggers_data, h_crossTriggers_data, "cp");
	  TGraphAsymmErrors *g_eff_MC = new TGraphAsymmErrors(h_crossAndDilepTriggers_MC, h_crossTriggers_MC, "cp");
	  TH1F *h_SF = (TH1F*)h_SF_nominal->Clone(("h_"+channels[i]+"_"+variables[j]+"_"+njetsystematics[2*l+1]+"_"+nvtxsystematics[2*k+1]+"_SystUncBand").c_str());
	  TGraphAsymmErrors *g_SF = new TGraphAsymmErrors(h_SF);
	  PropagateKatzStatError(h_crossAndDilepTriggers_data,h_crossTriggers_data,h_crossAndDilepTriggers_MC,h_crossTriggers_MC,h_SF,g_SF);

	  TH1F *h_SF_njets_high = (TH1F*)fSF.Get(("h_SF_"+channels[i]+"_"+variables[j]+"_"+njetsystematics[2*l+1]).c_str());
	  TH1F *h_SF_njets_low = (TH1F*)fSF.Get(("h_SF_"+channels[i]+"_"+variables[j]+"_"+njetsystematics[2*l]).c_str());
	  TH1F *h_SF_nvtx_high = (TH1F*)fSF.Get(("h_SF_"+channels[i]+"_"+variables[j]+"_"+nvtxsystematics[2*k+1]).c_str());
	  TH1F *h_SF_nvtx_low = (TH1F*)fSF.Get(("h_SF_"+channels[i]+"_"+variables[j]+"_"+nvtxsystematics[2*k]).c_str());

	  TH1F *uncBand_rel_diff_stat = (TH1F*)h_SF->Clone(("h_SF_"+channels[i]+"_"+variables[j]+"_stat_rel_diff_UncBand").c_str()); 
	  TH1F *uncBand_rel_diff_njets = (TH1F*)h_SF->Clone(("h_SF_"+channels[i]+"_"+variables[j]+"_njets_rel_diff_UncBand").c_str()); 
	  TH1F *uncBand_rel_diff_nvtx = (TH1F*)h_SF->Clone(("h_SF_"+channels[i]+"_"+variables[j]+"_nvtx_rel_diff_UncBand").c_str()); 
	  TH1F *uncBand = (TH1F*)h_SF->Clone(("h_SF_"+channels[i]+"_"+variables[j]+"_njetnvtxUncBand").c_str()); 

	  CalcUncBand_njet_nvtx(h_SF, h_SF_njets_high, h_SF_njets_low, h_SF_nvtx_high, h_SF_nvtx_low, uncBand, uncBand_rel_diff_stat, uncBand_rel_diff_njets, uncBand_rel_diff_nvtx);

	  TCanvas *c_eff = new TCanvas(("c_eff_"+channels[i]+"_"+variables[j]+"_"+njetsystematics[2*l+1]+"_"+nvtxsystematics[2*k+1]+"_SystUncBand").c_str(), "", 1000., 1000.);
	  TPad *p_eff = new TPad(("p_eff_"+channels[i]+"_"+variables[j]+"_"+njetsystematics[2*l+1]+"_"+nvtxsystematics[2*k+1]+"_SystUncBand").c_str(), "", 0, 0.15, 1, 1.0); 
	  setPad(p_eff);
	  drawGraph(h_eff_data, h_eff_MC,g_eff_data, g_eff_MC, "Efficiency", channels[i], xMins[j], xMaxs[j], era_);
	  c_eff->cd();
	  TPad *p_ratio_eff = new TPad(("p_ratio_eff_"+channels[i]+"_"+variables[j]+"_"+njetsystematics[2*l+1]+"_"+nvtxsystematics[2*k+1]+"_SystUncBand").c_str(), "", 0, 0.05, 1, 0.3);  
	  setPad(p_ratio_eff, false);
	  drawRatio(h_SF, g_SF, channels[i], xAxisTitles[j], xMins[j], xMaxs[j], "Scale Factor", 0.95, 1.05);

	  uncBand->SetFillColor(1);
	  uncBand->SetFillStyle(3354);
	  uncBand->SetMarkerStyle(0);
	  uncBand->Draw("same,e2");


	  if(era_ == "2016" || era_ == "2016preVFP_UL" || era_ == "2016postVFP_UL") {
	    if(njetsystematics[2*l+1] == "high3_njets" && nvtxsystematics[2*k+1] == "high20_nvtx") {
	      //	      c_eff->SaveAs(("triggereff/"+output_dir+"/g_"+channels[i]+"_"+variables[j]+"_SystUncBand"+".pdf").c_str());  
	      uncBand->Write();  
	      uncBand_rel_diff_stat->Write();  
	      uncBand_rel_diff_njets->Write();  
	      uncBand_rel_diff_nvtx->Write();  
	    }
	  }
	  else if(era_ == "2017" || era_ == "2017_UL") {	
	    if(njetsystematics[2*l+1] == "high3_njets" && nvtxsystematics[2*k+1] == "high30_nvtx") {
	      //	      c_eff->SaveAs(("triggereff/"+output_dir+"/g_"+channels[i]+"_"+variables[j]+"_SystUncBand"+".pdf").c_str());  
	      uncBand->Write();  
	      uncBand_rel_diff_stat->Write();  
	      uncBand_rel_diff_njets->Write();  
	      uncBand_rel_diff_nvtx->Write();  
	    }
	  }
	  else if(era_ == "2018" || era_ == "2018_UL") {	
	    if(njetsystematics[2*l+1] == "high3_njets" && nvtxsystematics[2*k+1] == "high30_nvtx") {
	      //	      c_eff->SaveAs(("triggereff/"+output_dir+"/g_"+channels[i]+"_"+variables[j]+"_SystUncBand"+".pdf").c_str());  
	      uncBand->Write();  
	      uncBand_rel_diff_stat->Write();  
	      uncBand_rel_diff_njets->Write();  
	      uncBand_rel_diff_nvtx->Write();  
	    }
	  }

	  //	  c_eff->SaveAs(("triggereff/"+output_dir+"/g_"+channels[i]+"_"+variables[j]+"_"+njetsystematics[2*l+1]+"_"+nvtxsystematics[2*k+1]+"_SystUncBand"+".pdf").c_str());  
	  //	  c_eff->SaveAs(("triggereff/"+output_dir+"/g_"+channels[i]+"_"+variables[j]+"_"+njetsystematics[2*l+1]+"_"+nvtxsystematics[2*k+1]+"_SystUncBand"+".C").c_str()); 
	  //	  c_eff->Write();
	}

      }
      

      /// To add EraRun systematic uncertainty band to 1D histograms

      TH1F *h_crossTriggers_data = (TH1F*)fdata->Get((channels[i]+"_"+variables[j]+"_"+systematics[0]+"_"+selections[2]).c_str());  
      TH1F *h_crossAndDilepTriggers_data = (TH1F*)fdata->Get((channels[i]+"_"+variables[j]+"_"+systematics[0]+"_"+selections[3]).c_str()); 
      //	  h_crossTriggers_data->Sumw2();
      //	  h_crossAndDilepTriggers_data->Sumw2(); 
      
      TH1F *h_crossTriggers_MC = (TH1F*)fMC[0]->Get((channels[i]+"_"+variables[j]+"_"+systematics[0]+"_"+selections[2]).c_str());
      TH1F *h_crossAndDilepTriggers_MC = (TH1F*)fMC[0]->Get((channels[i]+"_"+variables[j]+"_"+systematics[0]+"_"+selections[3]).c_str());

      TH1F *h_eff_data = (TH1F*)h_crossAndDilepTriggers_data->Clone(("h_eff_data_"+channels[i]+"_"+variables[j]+"_"+systematics[0]).c_str());
      h_eff_data->Divide(h_crossAndDilepTriggers_data,h_crossTriggers_data,1,1,"B");
      TH1F *h_eff_MC = (TH1F*)h_crossAndDilepTriggers_MC->Clone(("h_eff_MC_"+channels[i]+"_"+variables[j]+"_"+systematics[0]).c_str());
      h_eff_MC->Divide(h_crossAndDilepTriggers_MC,h_crossTriggers_MC,1,1,"B");
      
      vector<TH1F*> v_SF_EraRuns;

      TH1F *h_SF = (TH1F*)h_SF_nominal->Clone(("h_"+channels[i]+"_"+variables[j]+"_EraRun_SystUncBand").c_str());
      TGraphAsymmErrors *g_eff_data = new TGraphAsymmErrors(h_crossAndDilepTriggers_data, h_crossTriggers_data, "cp");
      TGraphAsymmErrors *g_eff_MC = new TGraphAsymmErrors(h_crossAndDilepTriggers_MC, h_crossTriggers_MC, "cp");
      TGraphAsymmErrors *g_SF = new TGraphAsymmErrors(h_SF);
      PropagateKatzStatError(h_crossAndDilepTriggers_data,h_crossTriggers_data,h_crossAndDilepTriggers_MC,h_crossTriggers_MC,h_SF,g_SF);

      TH1F *uncBand_njetnvtx = (TH1F*)intermediate_results->Get(("h_SF_"+channels[i]+"_"+variables[j]+"_njetnvtxUncBand").c_str());
      TH1F *uncBand_rel_diff_stat = (TH1F*)intermediate_results->Get(("h_SF_"+channels[i]+"_"+variables[j]+"_stat_rel_diff_UncBand").c_str());
      TH1F *uncBand_rel_diff_njets = (TH1F*)intermediate_results->Get(("h_SF_"+channels[i]+"_"+variables[j]+"_njets_rel_diff_UncBand").c_str());
      TH1F *uncBand_rel_diff_nvtx = (TH1F*)intermediate_results->Get(("h_SF_"+channels[i]+"_"+variables[j]+"_nvtx_rel_diff_UncBand").c_str());


      TH1F *uncBand_rel_diff_era = (TH1F*)uncBand_njetnvtx->Clone(("h_SF_"+channels[i]+"_"+variables[j]+"_era_rel_diff_UncBand").c_str());
      TH1F *uncBand = (TH1F*)uncBand_njetnvtx->Clone(("h_SF_"+channels[i]+"_"+variables[j]+"_Full_UncBand").c_str());
      TH1F *uncBand_rel_diff_Full = (TH1F*)uncBand_njetnvtx->Clone(("h_SF_"+channels[i]+"_"+variables[j]+"_Full_rel_diff_UncBand").c_str());

      for (UInt_t k = 0; k < EraRuns.size(); k++) {
	
	TH1F *h_SF_EraRun = (TH1F*)fSF.Get(("h_SF_"+channels[i]+"_"+variables[j]+"_nominal_"+EraRuns[k]).c_str());
	v_SF_EraRuns.push_back(h_SF_EraRun);

      }
      
      CalcUncBand_EraRun(h_SF, v_SF_EraRuns, EraRuns_Scale, uncBand, uncBand_rel_diff_era, uncBand_rel_diff_Full);
      
      TCanvas *c_eff = new TCanvas(("c_eff_"+channels[i]+"_"+variables[j]+"_FullSystUncBand").c_str(), "", 1000., 1000.);
      TPad *p_eff = new TPad(("p_eff_"+channels[i]+"_"+variables[j]+"_FullSystUncBand").c_str(), "", 0, 0.15, 1, 1.0); 
      setPad(p_eff);
      drawGraph(h_eff_data, h_eff_MC,g_eff_data, g_eff_MC, "Efficiency", channels[i], xMins[j], xMaxs[j], era_);
      c_eff->cd();
      TPad *p_ratio_eff = new TPad(("p_ratio_eff_"+channels[i]+"_"+variables[j]+"_FullSystUncBand").c_str(), "", 0, 0.05, 1, 0.3);  
      setPad(p_ratio_eff, false);
      drawRatio(h_SF, g_SF, channels[i], xAxisTitles[j], xMins[j], xMaxs[j], "Scale Factor", 0.95, 1.05);

      uncBand->SetFillColor(1);
      uncBand->SetFillStyle(3354);
      uncBand->SetMarkerStyle(0);
      uncBand->Draw("same,e2");
      uncBand->Write();
      uncBand_rel_diff_era->Write();
      uncBand_rel_diff_Full->Write();

      c_eff->SaveAs(("triggereff/"+output_dir+"/g_"+channels[i]+"_"+variables[j]+"_FullSystUncBand"+".pdf").c_str());  

      TCanvas *c_unc = new TCanvas(("c_unc_"+channels[i]+"_"+variables[j]).c_str(), "", 900., 600.);
      TPad *p_unc = new TPad(("p_unc_"+channels[i]+"_"+variables[j]).c_str(), "", 0, 0, 1, 1.0); 
      setPad(p_unc);
      DrawUncertainties(uncBand_rel_diff_Full, uncBand_rel_diff_stat, uncBand_rel_diff_njets, uncBand_rel_diff_nvtx, uncBand_rel_diff_era, channels[i], xAxisTitles[j], xMins[j], xMaxs[j], era_);

      c_unc->SaveAs(("triggereff/"+output_dir+"/g_"+channels[i]+"_"+variables[j]+"_ErrorsBreakdown"+".pdf").c_str());  
      c_unc->SaveAs(("triggereff/"+output_dir+"/g_"+channels[i]+"_"+variables[j]+"_ErrorsBreakdown"+".C").c_str());  


    ///
   
    }




  }

  

  intermediate_results->Close();

  TFile intermediate_results_(("triggereff/"+output_dir+"/intermediate_results"+output_dir+".root").c_str(),"READ");
  TFile results_(("triggereff/"+output_dir+"/results"+output_dir+".root").c_str(),"READ");
  TFile *final_results = new TFile(("triggereff/"+output_dir+"/TriggerSF_"+era_+".root").c_str(),"RECREATE");

  vector<string> multivariables_final = {"lepABpt", "lepABpt_eeIn", "lepABpt_eeOut", "lepABpt_eeSplit", "lepABpt_emuIn", "lepABpt_emuOut"};

  //  vector<string> variables_final = {"sidereal"};

  for (UInt_t i = 0; i < channels.size(); i++) {

    for (UInt_t j = 0; j < multivariables_final.size(); j++) {


      if(channels[i] == "ee") {
	if(multivariables_final[j] == "lepABpt_emuIn" || multivariables_final[j] == "lepABpt_emuOut") continue;
      }
      if(channels[i] == "emu") {
	if(multivariables_final[j] == "lepABpt_eeIn" || multivariables_final[j] == "lepABpt_eeOut" || multivariables_final[j] == "lepABpt_eeSplit") continue;
      }
      if(channels[i] == "mumu") {
	if(multivariables_final[j] == "lepABpt_emuIn" || multivariables_final[j] == "lepABpt_emuOut") continue;
	if(multivariables_final[j] == "lepABpt_eeIn" || multivariables_final[j] == "lepABpt_eeOut" || multivariables_final[j] == "lepABpt_eeSplit") continue;
      }

      TH2F *h2D_SF = (TH2F*)intermediate_results_.Get(("h2D_SF_"+channels[i]+"_"+multivariables_final[j]+"_FullError").c_str());
      h2D_SF->Write();
    }
    /* sidereal
    for (UInt_t j = 0; j < variables_final.size(); j++) {
      TH1F *h_SF = (TH1F*)results_.Get(("h_SF_"+channels[i]+"_"+variables_final[j]+"_nominal").c_str());
      h_SF->Write();
      TH1F *h_SF_Full_UncBand = (TH1F*)intermediate_results_.Get(("h_SF_"+channels[i]+"_"+variables_final[j]+"_Full_UncBand").c_str());
      h_SF_Full_UncBand->Write();
    }
    */
  }


  final_results->Close();

}


int main(int argc, char** argv) {

    CLParameter<std::string> opt_y("y", "Era", false, 1, 1);
    CLAnalyser::interpretGlobal(argc, argv);
    
    triggereffmakehistograms(opt_y[0]);

}
