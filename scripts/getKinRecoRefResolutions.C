#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TH2.h>
#include <TPaveText.h>
#include <TText.h>
#include <iostream>

void setStyle(TStyle* myStyle)
{
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetCanvasColor(kWhite);
  myStyle->SetCanvasDefH(600); //Height of canvas
  myStyle->SetCanvasDefW(600); //Width of canvas
  myStyle->SetCanvasDefX(0);   //Position on screen
  myStyle->SetCanvasDefY(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetPadColor(kWhite);
  myStyle->SetPadGridX(false);
  myStyle->SetPadGridY(false);
  myStyle->SetGridColor(0);
  myStyle->SetGridStyle(3);
  myStyle->SetGridWidth(1);
  myStyle->SetFrameBorderMode(0);
  myStyle->SetFrameBorderSize(1);
  myStyle->SetFrameFillColor(0);
  myStyle->SetFrameFillStyle(0);
  myStyle->SetFrameLineColor(1);
  myStyle->SetFrameLineStyle(1);
  myStyle->SetFrameLineWidth(1);
  myStyle->SetPadTopMargin(0.1);
  myStyle->SetPadBottomMargin(0.12);
  myStyle->SetPadLeftMargin(0.11);
  myStyle->SetPadRightMargin(0.15);
  myStyle->SetPaperSize(20.,20.);
  myStyle->SetOptStat(0);
  //myStyle->SetTitleX(0.25); // to set the x start position in NDC
  myStyle->SetTitleY(0.96); // to set the y position in NDC [0,1]
  //myStyle->SetTitleW(0.5); //to set the title box width
  //myStyle->SetTitleH(0.05);
}

void getReso(TString histoName, TString year)
{

  TString yearsuffix = year;
  year = year.ReplaceAll("UL","");

  TString outDir = "kinRecoRefResolutions/" + year + "/";

  if(histoName.Contains("mbl_true") || histoName.Contains("W_mass")) histoName += "_step0";
  else histoName += "_step7";

  TFile *emu_file = new TFile("./selectionRoot_" + year + "/Nominal/emu/emu_ttbarsignalplustau_fromDilepton_"+yearsuffix+".root");
  if(!emu_file) {
    std::cout<<"File in channel emu doesn't exist"<<std::endl;
    return;
  }

  setStyle(gStyle);

  TH1D *hist = (TH1D*) emu_file->Get(histoName);
  hist->SetDirectory(0);
  emu_file->Close();

  TFile *ee_file = new TFile("./selectionRoot_" + year + "/Nominal/ee/ee_ttbarsignalplustau_fromDilepton_"+yearsuffix+".root");
  if(!ee_file) {
    std::cout<<"File in channel ee doesn't exist"<<std::endl;
    return;
  }

  hist->Add((TH1D*) ee_file->Get(histoName));
  ee_file->Close();

  TFile *mumu_file = new TFile("./selectionRoot_" + year + "/Nominal/mumu/mumu_ttbarsignalplustau_fromDilepton_"+yearsuffix+".root");
  if(!mumu_file) {
    std::cout<<"File in channel mumu doesn't exist"<<std::endl;
    return;
  }

  hist->Add((TH1D*) mumu_file->Get(histoName));
  mumu_file->Close();        

  if(histoName.Contains("d_angle_jet")) hist->GetXaxis()->SetRangeUser(0.0, 0.2);
  if(histoName.Contains("d_angle_lep")) hist->GetXaxis()->SetRangeUser(0.0, 0.002);
  if(histoName.Contains("fE_jet")) hist->GetXaxis()->SetRangeUser(0.0, 2.8);
  if(histoName.Contains("fE_lep")) hist->GetXaxis()->SetRangeUser(0.86, 1.24);
  if(histoName.Contains("W_mass")) hist->GetXaxis()->SetRangeUser(66.0, 94.0);
  if(histoName.Contains("mbl_true")) hist->GetXaxis()->SetRangeUser(0.0, 180.0);
    
  if(histoName.Contains("mbl_true")) hist->Rebin(8);     

  /*TPaveText *text = new TPaveText(0.55,0.80,0.8, 0.9, "NDC");
    text->SetFillColor(0);
    text->SetFillStyle(0);
    text->SetBorderSize(0);*/


  /*TString channelLatex = "";
    if(channel == "mumu"){channelLatex = "#mu#mu channel";}
    else if (channel == "emu"){channelLatex = "e#mu channel";}
    else if (channel == "ee"){channelLatex = "ee channel";}
    text->AddText(channelLatex);*/

  gSystem->mkdir(outDir);
  TCanvas *c = new TCanvas("", "", 0,0, 1000, 750);

  //gStyle->SetPaintTextFormat("4.2f");

  if(histoName.Contains("d_angle_lep") || histoName.Contains("d_angle_jet")) TGaxis::SetMaxDigits(3);

  TString xtitle = "";
  if(histoName.Contains("d_angle_jet")) xtitle = "#alpha [rad]";
  if(histoName.Contains("d_angle_lep")) xtitle = "#alpha [rad]";
  if(histoName.Contains("fE_jet")) xtitle = "E^{true}_{jet}/E^{reco}_{jet}";
  if(histoName.Contains("fE_lep")) xtitle = "E^{true}_{l}/E^{reco}_{l}";
  if(histoName.Contains("W_mass")) xtitle = "m_{W} [GeV]";
  if(histoName.Contains("mbl_true")) xtitle = "m_{lb-system} [GeV]";
  TString ytitle = "a.u.";
  TString title = "";
  if(histoName.Contains("d_angle_jet")) title = "Jet smearing angle";
  if(histoName.Contains("d_angle_lep")) title = "Lepton smearing angle";
  if(histoName.Contains("fE_jet")) title = "Jet energy correction factor";
  if(histoName.Contains("fE_lep")) title = "Lepton energy correction factor";
  if(histoName.Contains("W_mass")) title = "W boson mass";
  if(histoName.Contains("mbl_true")) title = "Mass of lepton-b-jet system";
  hist->GetXaxis()->SetTitle(xtitle);
  hist->GetYaxis()->SetTitle(ytitle);
  hist->GetYaxis()->SetTitleOffset(1.05*hist->GetYaxis()->GetTitleOffset());
  hist->GetXaxis()->SetTitleOffset(0.95*hist->GetXaxis()->GetTitleOffset());
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->SetTitle(title);
  hist->SetLineColor(kRed+1);
  hist->SetLineWidth(3);
  hist->SetLineStyle(1);
  hist->SetMarkerColor(kRed+1);
  /*beff->GetZaxis()->SetRangeUser(0.6,0.9);
    beff->Draw("colz,text,E");
    text->Draw("same");
    c->SetLogx();*/
  TH1F *h2 = (TH1F*) hist->Clone();
  //hist->Draw();
  hist->Draw("c hist same");
  //h2->SetLineColor(kBlue+1);
  //h2->Draw("C same");
  c->Update();
  c->Print(outDir.Copy()+histoName+".eps");
  c->Print(outDir.Copy()+histoName+".pdf");
  c->Print(outDir.Copy()+histoName+".C");
  //c->Print(outDir.Copy()+histoName+".root");

  delete c;

}




void getKinRecoRefResolutions(TString year)
{

  TString histos[6] = {"KinReco_d_angle_jet", "KinReco_d_angle_lep", "KinReco_fE_jet", "KinReco_fE_lep", "KinReco_W_mass", "KinReco_mbl_true"};
  //TString histos[2] = {"KinReco_W_mass", "KinReco_mbl_true"};
  for (unsigned int iter =0; iter < 6; iter ++){ // set to 6
    getReso(histos[iter], year);
  };
}
