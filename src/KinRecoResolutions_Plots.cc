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

// To run in root, example:
// root -l -q -b src/KinRecoResolutions_Plots.cc++'("prompt")'

/// Functions

double CalcLumiWeight(TFile* fMC_, double XSection_, string era_) {

  double lumi = 0;
  
  if (era_ == "2016preVFP")  lumi = 19500;
  else if (era_ == "2016postVFP")  lumi = 16810;
  else if (era_ == "2017")  lumi = 41480;
  else if (era_ == "2018")  lumi = 59830;

  TH1D *h_NrOfEvts = (TH1D*)fMC_->Get("weightedEvents");  
  
  double NrOfEvts = h_NrOfEvts->GetBinContent(1);
  //  double NrOfEvts = h_NrOfEvts->GetEntries();

  double lumiWeight = lumi * XSection_ / NrOfEvts;


  std::cout << "lumiw=" << lumiWeight <<", lumi=" << lumi <<", XSection=" << XSection_ <<", NrOfEvts=" << NrOfEvts << "\n";
 

  return lumiWeight;

}


void setPad(TPad *p, bool hasGridX_ = true) {

  p->SetLeftMargin(0.15);
  p->SetRightMargin(0.05);
  p->SetBottomMargin(0.1);
  p->SetTopMargin(0.1);

  if (hasGridX_) p->SetGridx();
  p->SetGridy();
  p->SetTickx(1);
  p->SetTicky(1);

  //p->Range(80,0.75,220,1.25);
  
  p->Draw();
  p->cd();

  gStyle->SetPaintTextFormat("3.3f");

}


void set2DPad(TPad *p) {

  p->SetLeftMargin(0.125);
  p->SetRightMargin(0.125);
  p->SetBottomMargin(0.125);
  p->SetTopMargin(0.1);
  p->SetTickx(1);
  p->SetTicky(1);

  p->SetLogz();

  p->Draw();
  p->cd();

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

}


/// Main Function                                                                                                                                                               

void KinRecoResolutions_Plots(std::string prompt_or_viatau_){

  system("mkdir -p KinRecoResolutions/Plots");

  // Combine luminosity scaled histograms from the four eras

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  vector<TFile*> fMC;
  vector<double> fMC_LumiWeight;

  double topxsec = 830.91;

  if (prompt_or_viatau_ == "prompt") {

    TFile* fMC_prompt_2016preVFP = new TFile("KinRecoResolutions/output_hists_prompt_2016preVFP.root","READ");
    fMC.push_back(fMC_prompt_2016preVFP);  
    fMC_LumiWeight.push_back(CalcLumiWeight(fMC_prompt_2016preVFP, topxsec * 0.10706 * 0.957072325, "2016preVFP"));

    TFile* fMC_prompt_2016postVFP = new TFile("KinRecoResolutions/output_hists_prompt_2016postVFP.root","READ");
    fMC.push_back(fMC_prompt_2016postVFP);  
    fMC_LumiWeight.push_back(CalcLumiWeight(fMC_prompt_2016postVFP, topxsec * 0.10706 * 0.957072325, "2016postVFP"));

    TFile* fMC_prompt_2017 = new TFile("KinRecoResolutions/output_hists_prompt_2017.root","READ");
    fMC.push_back(fMC_prompt_2017);  
    fMC_LumiWeight.push_back(CalcLumiWeight(fMC_prompt_2017, topxsec * 0.10706 * 0.957072325, "2017"));

    TFile* fMC_prompt_2018 = new TFile("KinRecoResolutions/output_hists_prompt_2018.root","READ");
    fMC.push_back(fMC_prompt_2018);  
    fMC_LumiWeight.push_back(CalcLumiWeight(fMC_prompt_2018, topxsec * 0.10706 * 0.957072325, "2018"));

  }

  else if (prompt_or_viatau_ == "viatau") {

    TFile* fMC_viatau_2016preVFP = new TFile("KinRecoResolutions/output_hists_viatau_2016preVFP.root","READ");
    fMC.push_back(fMC_viatau_2016preVFP);  
    fMC_LumiWeight.push_back(CalcLumiWeight(fMC_viatau_2016preVFP, topxsec * 0.10706 * 1.026252753, "2016preVFP"));

    TFile* fMC_viatau_2016postVFP = new TFile("KinRecoResolutions/output_hists_viatau_2016postVFP.root","READ");
    fMC.push_back(fMC_viatau_2016postVFP);  
    fMC_LumiWeight.push_back(CalcLumiWeight(fMC_viatau_2016postVFP, topxsec * 0.10706 * 1.026252753, "2016postVFP"));

    TFile* fMC_viatau_2017 = new TFile("KinRecoResolutions/output_hists_viatau_2017.root","READ");
    fMC.push_back(fMC_viatau_2017);  
    fMC_LumiWeight.push_back(CalcLumiWeight(fMC_viatau_2017, topxsec * 0.10706 * 1.026252753, "2017"));

    TFile* fMC_viatau_2018 = new TFile("KinRecoResolutions/output_hists_viatau_2018.root","READ");
    fMC.push_back(fMC_viatau_2018);  
    fMC_LumiWeight.push_back(CalcLumiWeight(fMC_viatau_2018, topxsec * 0.10706 * 1.026252753, "2018"));

  }


  TFile *output_hists = new TFile(("KinRecoResolutions/Plots/hists_"+prompt_or_viatau_+".root").c_str(),"RECREATE");


  vector<string> variables = {"ttbar_mass_residual", "ttbar_pT_residual", "ttbar_rapidity_residual", "top_pT_residual", "top_rapidity_residual", "top_phi_residual", "top_scatteringangle_ttbarframe_residual"};
  vector<string> variablesXAxisTitles = {"M(t#bar{t})^{gen} - M(t#bar{t})^{rec} [GeV]", "p_{T}(t#bar{t})^{gen} - p_{T}(t#bar{t})^{rec} [GeV]", "y(t#bar{t})^{gen} - y(t#bar{t})^{rec}", "p_{T}(t)^{gen} - p_{T}(t)^{rec} [GeV]", "y(t)^{gen} - y(t)^{rec}", "#phi(t)^{gen} - #phi(t)^{rec}", "cos#Theta(t)^{gen} - cos#Theta(t)^{rec}"};
  vector<string> variablesYAxisTitles = {"Top Quark/Anti-Quark Pairs", "Top Quark/Anti-Quark Pairs", "Top Quark/Anti-Quark Pairs", "Top Quarks", "Top Quarks", "Top Quarks", "Top Quarks"};

  vector<string> multivariables = {"ttbar_mass_genreco", "ttbar_pT_genreco", "ttbar_rapidity_genreco", "top_pT_genreco", "top_rapidity_genreco", "top_phi_genreco", "top_scatteringangle_ttbarframe_genreco"};
  vector<string> multivariablesXAxisTitles = {"M(t#bar{t})^{gen} [GeV]", "p_{T}(t#bar{t})^{gen} [GeV]", "y(t#bar{t})^{gen}", "p_{T}(t)^{gen} [GeV]", "y(t)^{gen}", "#phi(t)^{gen}", "cos#Theta(t)^{gen}"};
  vector<string> multivariablesYAxisTitles = {"M(t#bar{t})^{rec} [GeV]", "p_{T}(t#bar{t})^{rec} [GeV]", "y(t#bar{t})^{rec}", "p_{T}(t)^{rec} [GeV]", "y(t)^{rec}", "#phi(t)^{rec}", "cos#Theta(t)^{rec}"};

  vector<string> multiresidualvariables = {"ttbar_mass_multiresidual", "ttbar_pT_multiresidual", "ttbar_rapidity_multiresidual", "top_pT_multiresidual", "top_rapidity_multiresidual", "top_phi_multiresidual", "top_scatteringangle_ttbarframe_multiresidual"};
  vector<string> multiresidualvariablesXAxisTitles = {"M(t#bar{t})^{gen} [GeV]", "p_{T}(t#bar{t})^{gen} [GeV]", "y(t#bar{t})^{gen}", "p_{T}(t)^{gen} [GeV]", "y(t)^{gen}", "#phi(t)^{gen}", "cos#Theta(t)^{gen}"};
  vector<string> multiresidualvariablesYAxisTitles = {"Mean & Sigma [M(t#bar{t})^{gen} - M(t#bar{t})^{rec}] [GeV]", "Mean & Sigma [p_{T}(t#bar{t})^{gen} - p_{T}(t#bar{t})^{rec}] [GeV]", "Mean & Sigma [y(t#bar{t})^{gen} - y(t#bar{t})^{rec}]", "Mean & Sigma [p_{T}(t)^{gen} - p_{T}(t)^{rec}] [GeV]", "Mean & Sigma [y(t)^{gen} - y(t)^{rec}]", "Mean & Sigma [#phi(t)^{gen} - #phi(t)^{rec}]", "Mean & Sigma [cos#Theta(t)^{gen} - cos#Theta(t)^{rec}]"};
  vector<double> multiresidualvariablesMins = {-200, -50, -1.0, -100, -1.2, -0.5, -0.5};
  vector<double> multiresidualvariablesMaxs = {600, 100, 1.0, 200, 1.2, 0.5, 0.5};
  

  for (UInt_t i = 0; i < variables.size(); i++) {

    std::cout << variables[i] << std::endl;

    TH1D *h_var = (TH1D*)fMC[0]->Get(variables[i].c_str());
    h_var->Sumw2();
    h_var->Scale(fMC_LumiWeight[0]);

    TH1D *hs_var = (TH1D*)h_var->Clone(variables[i].c_str());

    for (UInt_t j = 1; j < fMC.size(); j++) {

      h_var = (TH1D*)fMC[j]->Get((variables[i]).c_str());
      h_var->Sumw2();
      h_var->Scale(fMC_LumiWeight[j]);

      hs_var->Add(h_var);

    }

    hs_var->Write();

  }


  for (UInt_t i = 0; i < multivariables.size(); i++) {

    std::cout << multivariables[i] << std::endl;

    TH2D *h_multivar = (TH2D*)fMC[0]->Get((multivariables[i]).c_str());
    h_multivar->Sumw2();
    h_multivar->Scale(fMC_LumiWeight[0]);

    TH2D *hs_multivar = (TH2D*)h_multivar->Clone(multivariables[i].c_str());

    for (UInt_t j = 1; j < fMC.size(); j++) {

      h_multivar = (TH2D*)fMC[j]->Get((multivariables[i]).c_str());
      h_multivar->Sumw2();
      h_multivar->Scale(fMC_LumiWeight[j]);

      hs_multivar->Add(h_multivar);

    }

    hs_multivar->Write();

  }


  for (UInt_t i = 0; i < multiresidualvariables.size(); i++) {

    std::cout << multiresidualvariables[i] << std::endl;

    TH2D *h_multiresidualvar = (TH2D*)fMC[0]->Get((multiresidualvariables[i]).c_str());
    h_multiresidualvar->Sumw2();
    h_multiresidualvar->Scale(fMC_LumiWeight[0]);

    TH2D *hs_multiresidualvar = (TH2D*)h_multiresidualvar->Clone(multiresidualvariables[i].c_str());

    for (UInt_t j = 1; j < fMC.size(); j++) {

      h_multiresidualvar = (TH2D*)fMC[j]->Get((multiresidualvariables[i]).c_str());
      h_multiresidualvar->Sumw2();
      h_multiresidualvar->Scale(fMC_LumiWeight[j]);

      hs_multiresidualvar->Add(h_multiresidualvar);

    }

    hs_multiresidualvar->Write();

    hs_multiresidualvar->FitSlicesY();

    TH1D *h_multiresidualvar_0 = (TH1D*)gDirectory->Get((multiresidualvariables[i]+"_0").c_str());
    h_multiresidualvar_0->Write();
    TH1D *h_multiresidualvar_1 = (TH1D*)gDirectory->Get((multiresidualvariables[i]+"_1").c_str());
    h_multiresidualvar_1->Write();
    TH1D *h_multiresidualvar_2 = (TH1D*)gDirectory->Get((multiresidualvariables[i]+"_2").c_str());
    h_multiresidualvar_2->Write();
  }

  output_hists->Close();


  // Make plots from luminosity scaled histograms

  TFile *hists = new TFile(("KinRecoResolutions/Plots/hists_"+prompt_or_viatau_+".root").c_str(),"READ");

  TFile *output_plots = new TFile(("KinRecoResolutions/Plots/plots_"+prompt_or_viatau_+".root").c_str(),"RECREATE");

  for (UInt_t i = 0; i < variables.size(); i++) {

    std::cout << variables[i] << std::endl;

    TH1D *h_var = (TH1D*)hists->Get((variables[i]).c_str());
    h_var->GetXaxis()->SetTitle(variablesXAxisTitles[i].c_str());
    h_var->GetYaxis()->SetTitle(variablesYAxisTitles[i].c_str());

    TCanvas *c_var = new TCanvas(("c_var_"+variables[i]).c_str(), "", 800., 800.);
    c_var->cd();
    TPad *p_var = new TPad(("p_var_"+variables[i]).c_str(), "", 0, 0, 1, 1); 
    setPad(p_var);

    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat("mruo");

    double mean = h_var->GetBinCenter( h_var->GetMaximumBin() );
    double sigma = h_var->GetRMS();

    TF1 *f1 = new TF1("f1","gaus", mean - sigma, mean + sigma);

    h_var->Fit("f1", "RMQ");

    h_var->GetXaxis()->SetNdivisions(510);

    h_var->GetXaxis()->SetLabelSize(0.03);
    h_var->GetXaxis()->SetTitleSize(0.04);
    h_var->GetXaxis()->SetTitleOffset(1.0);

    h_var->GetYaxis()->SetNdivisions(510);
    h_var->GetYaxis()->SetMaxDigits(3);

    h_var->GetYaxis()->SetLabelSize(0.03);
    h_var->GetYaxis()->SetLabelOffset(0.01);
    h_var->GetYaxis()->SetTitleSize(0.04);
    h_var->GetYaxis()->SetTitleOffset(1.4);

    h_var->SetLineWidth(2);
    h_var->SetLineColor(kAzure);
    h_var->SetMarkerStyle(20);
    h_var->SetMarkerSize(1.2);
    h_var->SetMarkerColor(kAzure);

    h_var->Draw();

    TLatex cmslabel;
    cmslabel.SetNDC();
    cmslabel.DrawLatex(0.225, 0.94, "#font[61]{CMS}   #scale[0.76]{#font[52]{Work in Progress}}   #scale[0.76]{#font[52]{Simulation}}");

    TLatex label;
    label.SetNDC();
    label.DrawLatex(0.175, 0.85, ("#scale[0.76]{#font[61]{"+prompt_or_viatau_+"}}").c_str());

    c_var->SaveAs(("KinRecoResolutions/Plots/"+variables[i]+"_"+prompt_or_viatau_+".pdf").c_str());
    c_var->SaveAs(("KinRecoResolutions/Plots/"+variables[i]+"_"+prompt_or_viatau_+".C").c_str());
    c_var->Write();
  }

  for (UInt_t i = 0; i < multivariables.size(); i++) {

    std::cout << multivariables[i] << std::endl;

    TH2D *h_multivar = (TH2D*)hists->Get((multivariables[i]).c_str());
    h_multivar->GetXaxis()->SetTitle(multivariablesXAxisTitles[i].c_str());
    h_multivar->GetYaxis()->SetTitle(multivariablesYAxisTitles[i].c_str());

    TCanvas *c_multivar = new TCanvas(("c_multivar_"+multivariables[i]).c_str(), "", 800., 800.);
    c_multivar->cd();
    TPad *p_multivar = new TPad(("p_multivar_"+multivariables[i]).c_str(), "", 0, 0, 1, 1); 
    set2DPad(p_multivar);

    gStyle->SetOptTitle(0);

    h_multivar->GetXaxis()->SetNdivisions(510);
    //h_multivar->GetXaxis()->SetMaxDigits(3);
    h_multivar->GetYaxis()->SetNdivisions(510);
    h_multivar->GetYaxis()->SetMaxDigits(3);

    h_multivar->GetXaxis()->SetLabelFont(42);
    h_multivar->GetXaxis()->SetLabelSize(0.035);
    h_multivar->GetXaxis()->SetLabelOffset(0.01);
    h_multivar->GetXaxis()->SetTitleFont(42);
    h_multivar->GetXaxis()->SetTitleSize(0.040);
    h_multivar->GetXaxis()->SetTitleOffset(1.3);

    h_multivar->GetYaxis()->SetLabelFont(42);
    h_multivar->GetYaxis()->SetLabelSize(0.035);
    h_multivar->GetYaxis()->SetLabelOffset(0.01);
    h_multivar->GetYaxis()->SetTitleFont(42);
    h_multivar->GetYaxis()->SetTitleSize(0.045);
    h_multivar->GetYaxis()->SetTitleOffset(1.3);

    h_multivar->GetZaxis()->SetLabelFont(42);
    h_multivar->GetZaxis()->SetLabelSize(0.035);
    h_multivar->GetZaxis()->SetLabelOffset(0.005);
    h_multivar->GetZaxis()->SetTitleFont(42);
    h_multivar->GetZaxis()->SetTitleSize(0.045);

    h_multivar->SetMarkerSize(0.7);

    h_multivar->Draw("colz");

    TLatex cmslabel;
    cmslabel.SetNDC();
    cmslabel.DrawLatex(0.2, 0.925, "#font[61]{CMS}   #scale[0.76]{#font[52]{Work in Progress}}   #scale[0.76]{#font[52]{Simulation}}");

    TLatex label;
    label.SetNDC();
    label.DrawLatex(0.175, 0.85, ("#scale[0.76]{#font[61]{"+prompt_or_viatau_+"}}").c_str());

    c_multivar->SaveAs(("KinRecoResolutions/Plots/"+multivariables[i]+"_"+prompt_or_viatau_+".pdf").c_str());
    c_multivar->SaveAs(("KinRecoResolutions/Plots/"+multivariables[i]+"_"+prompt_or_viatau_+".C").c_str());
    c_multivar->Write();
  }

  for (UInt_t i = 0; i < multiresidualvariables.size(); i++) {

    std::cout << multiresidualvariables[i] << std::endl;

    TH1D *h_multiresidualvar_0 = (TH1D*)hists->Get((multiresidualvariables[i]+"_0").c_str());
    TH1D *h_multiresidualvar_1 = (TH1D*)hists->Get((multiresidualvariables[i]+"_1").c_str());
    TH1D *h_multiresidualvar_2 = (TH1D*)hists->Get((multiresidualvariables[i]+"_2").c_str());
    h_multiresidualvar_1->GetXaxis()->SetTitle(multiresidualvariablesXAxisTitles[i].c_str());
    h_multiresidualvar_1->GetYaxis()->SetTitle(multiresidualvariablesYAxisTitles[i].c_str());

    TCanvas *c_multiresidualvar = new TCanvas(("c_multiresidualvar_"+multiresidualvariables[i]).c_str(), "", 800., 800.);
    c_multiresidualvar->cd();
    TPad *p_multiresidualvar = new TPad(("p_multiresidualvar_"+multiresidualvariables[i]).c_str(), "", 0, 0, 1, 1); 
    setPad(p_multiresidualvar);

    gStyle->SetOptTitle(0);

    h_multiresidualvar_1->GetXaxis()->SetNdivisions(510);

    h_multiresidualvar_1->GetXaxis()->SetLabelSize(0.03);
    h_multiresidualvar_1->GetXaxis()->SetTitleSize(0.04);
    h_multiresidualvar_1->GetXaxis()->SetTitleOffset(1.0);

    h_multiresidualvar_1->GetYaxis()->SetNdivisions(510);
    h_multiresidualvar_1->GetYaxis()->SetMaxDigits(3);

    h_multiresidualvar_1->GetYaxis()->SetLabelSize(0.03);
    h_multiresidualvar_1->GetYaxis()->SetLabelOffset(0.01);
    h_multiresidualvar_1->GetYaxis()->SetTitleSize(0.04);
    h_multiresidualvar_1->GetYaxis()->SetTitleOffset(1.4);

    h_multiresidualvar_1->SetLineWidth(2);
    h_multiresidualvar_1->SetLineColor(kAzure);
    h_multiresidualvar_1->SetMarkerStyle(20);
    h_multiresidualvar_1->SetMarkerSize(1.2);
    h_multiresidualvar_1->SetMarkerColor(kAzure);

    h_multiresidualvar_2->SetLineWidth(2);
    h_multiresidualvar_2->SetLineColor(kRed);
    h_multiresidualvar_2->SetMarkerStyle(20);
    h_multiresidualvar_2->SetMarkerSize(1.2);
    h_multiresidualvar_2->SetMarkerColor(kRed);

    h_multiresidualvar_1->SetMinimum(multiresidualvariablesMins[i]);
    h_multiresidualvar_1->SetMaximum(multiresidualvariablesMaxs[i]);

    h_multiresidualvar_1->Draw();
    h_multiresidualvar_2->Draw("SAME");

    TLegend *leg = new TLegend(0.65,0.775,0.75,0.875);  
    leg->AddEntry(h_multiresidualvar_1, "Mean", "lep");
    leg->AddEntry(h_multiresidualvar_2, "Sigma", "lep");
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->Draw("SAME");

    TLatex cmslabel;
    cmslabel.SetNDC();
    cmslabel.DrawLatex(0.2, 0.925, "#font[61]{CMS}   #scale[0.76]{#font[52]{Work in Progress}}   #scale[0.76]{#font[52]{Simulation}}");

    TLatex label;
    label.SetNDC();
    label.DrawLatex(0.175, 0.85, ("#scale[0.76]{#font[61]{"+prompt_or_viatau_+"}}").c_str());

    c_multiresidualvar->SaveAs(("KinRecoResolutions/Plots/"+multiresidualvariables[i]+"_"+prompt_or_viatau_+".pdf").c_str());
    c_multiresidualvar->SaveAs(("KinRecoResolutions/Plots/"+multiresidualvariables[i]+"_"+prompt_or_viatau_+".C").c_str());
    c_multiresidualvar->Write();
  }


  output_plots->Close();

}

int main(int argc, char** argv) {

  CLParameter<std::string> opt_t("t", "Prompt or via Tau", false, 1, 1);
  CLAnalyser::interpretGlobal(argc, argv);


  KinRecoResolutions_Plots(opt_t[0]);

}
