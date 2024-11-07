
void Linearity_summary()
{
    gStyle->SetPadLeftMargin(0.08);
    gStyle->SetPadRightMargin(0.02);
    //gStyle->SetTitleSize(0.08, "XYZ");
    //gStyle->SetLabelSize(0.07, "XYZ");
    gStyle->SetTickLength(0.01, "XYZ");
    gStyle->SetOptTitle(0);
    gStyle->SetPadGridY(1);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetOptStat("");


  // TString Filename1 = "Linearity_summary_correctBFF.root";
  //TString Filename1 = "Linearity_summary_reweightbiastest_BFFonly.root";
  //TString Filename2 = "Linearity_summary_defaultTUR.root";
  //TString Filename2 = "Linearity_summary_reweightbiastest.root";
  //TString Filename2 = "Linearity_summary_reweightbiastest_BFFonly.root";

  //TString Filename2 = "Linearity_summary_reweightbiastest_mx3.root";
  
  //TString Filename2 = "Linearity_summary_1stbinonlyWlin.root";
  //TString Filename2 = "Linearity_summary_3rdbinonlylin.root";
  //TString Filename2 = "Linearity_summary_RWwide3.root";
  // TString Filename2 = "Linearity_summary_righthalfonlylin.root";
  

  TString Filename1 = "TUnfoldResults_lin_16_simpleplotsforAN/Nominal/mumu/Linearity_summary.root";
  TString Filename2 = "~/DESY_test/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/TUnfoldResults_lin_120/Nominal/mumu/Linearity_summary.root";

	TFile *File1 = new TFile(Filename1.Data(), "READ");
	TFile *File2 = new TFile(Filename2.Data(), "READ");

	TH1D* hAfbLinearityMaxBiasFrac_rebinnedA_1 = (TH1D*)File1->Get("hAfbLinearityMaxBiasFrac_rebinnedA");
	TH1D* hAfbLinearityMaxBiasFracStat_rebinnedA_1 = (TH1D*)File1->Get("hAfbLinearityMaxBiasFracStat_rebinnedA");
	TH1D* hAfbLinearityMaxBiasSlope_rebinnedA_1 = (TH1D*)File1->Get("hAfbLinearityMaxBiasSlope_rebinnedA");
	TH1D* hAfbLinearityMaxBiasFracMax_rebinnedA_1 = (TH1D*)File1->Get("hAfbLinearityMaxBiasFracMax_rebinnedA");
	TH1D* hAfbLinearityMaxBiasFracStatMax_rebinnedA_1 = (TH1D*)File1->Get("hAfbLinearityMaxBiasFracStatMax_rebinnedA");
	TH1D* hAfbLinearityMaxBiasSlopeMax_rebinnedA_1 = (TH1D*)File1->Get("hAfbLinearityMaxBiasSlopeMax_rebinnedA");
	TH1D* hAfbLinearityAvgBiasFrac_rebinnedA_1 = (TH1D*)File1->Get("hAfbLinearityAvgBiasFrac_rebinnedA");
	TH1D* hAfbLinearityAvgBiasFracStat_rebinnedA_1 = (TH1D*)File1->Get("hAfbLinearityAvgBiasFracStat_rebinnedA");
	TH1D* hAfbLinearityAvgBiasSlope_rebinnedA_1 = (TH1D*)File1->Get("hAfbLinearityAvgBiasSlope_rebinnedA");
	TH1D* hAfbLinearityMaxBiasFrac_rebinnedB_1 = (TH1D*)File1->Get("hAfbLinearityMaxBiasFrac_rebinnedB");
	TH1D* hAfbLinearityMaxBiasFracStat_rebinnedB_1 = (TH1D*)File1->Get("hAfbLinearityMaxBiasFracStat_rebinnedB");
	TH1D* hAfbLinearityMaxBiasSlope_rebinnedB_1 = (TH1D*)File1->Get("hAfbLinearityMaxBiasSlope_rebinnedB");
	TH1D* hAfbLinearityMaxBiasFracMax_rebinnedB_1 = (TH1D*)File1->Get("hAfbLinearityMaxBiasFracMax_rebinnedB");
	TH1D* hAfbLinearityMaxBiasFracStatMax_rebinnedB_1 = (TH1D*)File1->Get("hAfbLinearityMaxBiasFracStatMax_rebinnedB");
	TH1D* hAfbLinearityMaxBiasSlopeMax_rebinnedB_1 = (TH1D*)File1->Get("hAfbLinearityMaxBiasSlopeMax_rebinnedB");
	TH1D* hAfbLinearityAvgBiasFrac_rebinnedB_1 = (TH1D*)File1->Get("hAfbLinearityAvgBiasFrac_rebinnedB");
	TH1D* hAfbLinearityAvgBiasFracStat_rebinnedB_1 = (TH1D*)File1->Get("hAfbLinearityAvgBiasFracStat_rebinnedB");
	TH1D* hAfbLinearityAvgBiasSlope_rebinnedB_1 = (TH1D*)File1->Get("hAfbLinearityAvgBiasSlope_rebinnedB");

	TH1D* hAfbLinearityMaxBiasFrac_rebinnedA_2 = (TH1D*)File2->Get("hAfbLinearityMaxBiasFrac_rebinnedA");
	TH1D* hAfbLinearityMaxBiasFracStat_rebinnedA_2 = (TH1D*)File2->Get("hAfbLinearityMaxBiasFracStat_rebinnedA");
	TH1D* hAfbLinearityMaxBiasSlope_rebinnedA_2 = (TH1D*)File2->Get("hAfbLinearityMaxBiasSlope_rebinnedA");
	TH1D* hAfbLinearityMaxBiasFracMax_rebinnedA_2 = (TH1D*)File2->Get("hAfbLinearityMaxBiasFracMax_rebinnedA");
	TH1D* hAfbLinearityMaxBiasFracStatMax_rebinnedA_2 = (TH1D*)File2->Get("hAfbLinearityMaxBiasFracStatMax_rebinnedA");
	TH1D* hAfbLinearityMaxBiasSlopeMax_rebinnedA_2 = (TH1D*)File2->Get("hAfbLinearityMaxBiasSlopeMax_rebinnedA");
	TH1D* hAfbLinearityAvgBiasFrac_rebinnedA_2 = (TH1D*)File2->Get("hAfbLinearityAvgBiasFrac_rebinnedA");
	TH1D* hAfbLinearityAvgBiasFracStat_rebinnedA_2 = (TH1D*)File2->Get("hAfbLinearityAvgBiasFracStat_rebinnedA");
	TH1D* hAfbLinearityAvgBiasSlope_rebinnedA_2 = (TH1D*)File2->Get("hAfbLinearityAvgBiasSlope_rebinnedA");
	TH1D* hAfbLinearityMaxBiasFrac_rebinnedB_2 = (TH1D*)File2->Get("hAfbLinearityMaxBiasFrac_rebinnedB");
	TH1D* hAfbLinearityMaxBiasFracStat_rebinnedB_2 = (TH1D*)File2->Get("hAfbLinearityMaxBiasFracStat_rebinnedB");
	TH1D* hAfbLinearityMaxBiasSlope_rebinnedB_2 = (TH1D*)File2->Get("hAfbLinearityMaxBiasSlope_rebinnedB");
	TH1D* hAfbLinearityMaxBiasFracMax_rebinnedB_2 = (TH1D*)File2->Get("hAfbLinearityMaxBiasFracMax_rebinnedB");
	TH1D* hAfbLinearityMaxBiasFracStatMax_rebinnedB_2 = (TH1D*)File2->Get("hAfbLinearityMaxBiasFracStatMax_rebinnedB");
	TH1D* hAfbLinearityMaxBiasSlopeMax_rebinnedB_2 = (TH1D*)File2->Get("hAfbLinearityMaxBiasSlopeMax_rebinnedB");
	TH1D* hAfbLinearityAvgBiasFrac_rebinnedB_2 = (TH1D*)File2->Get("hAfbLinearityAvgBiasFrac_rebinnedB");
	TH1D* hAfbLinearityAvgBiasFracStat_rebinnedB_2 = (TH1D*)File2->Get("hAfbLinearityAvgBiasFracStat_rebinnedB");
	TH1D* hAfbLinearityAvgBiasSlope_rebinnedB_2 = (TH1D*)File2->Get("hAfbLinearityAvgBiasSlope_rebinnedB");

	TH1D* hAfbLinearityMaxBiasFrac_rebinnedA_2c = (TH1D*) hAfbLinearityMaxBiasFrac_rebinnedA_2->Clone();
	TH1D* hAfbLinearityMaxBiasFracStat_rebinnedA_2c = (TH1D*) hAfbLinearityMaxBiasFracStat_rebinnedA_2->Clone();
	TH1D* hAfbLinearityMaxBiasSlope_rebinnedA_2c = (TH1D*) hAfbLinearityMaxBiasSlope_rebinnedA_2->Clone();
	TH1D* hAfbLinearityMaxBiasFracMax_rebinnedA_2c = (TH1D*) hAfbLinearityMaxBiasFracMax_rebinnedA_2->Clone();
	TH1D* hAfbLinearityMaxBiasFracStatMax_rebinnedA_2c = (TH1D*) hAfbLinearityMaxBiasFracStatMax_rebinnedA_2->Clone();
	TH1D* hAfbLinearityMaxBiasSlopeMax_rebinnedA_2c = (TH1D*) hAfbLinearityMaxBiasSlopeMax_rebinnedA_2->Clone();
	TH1D* hAfbLinearityAvgBiasFrac_rebinnedA_2c = (TH1D*) hAfbLinearityAvgBiasFrac_rebinnedA_2->Clone();
	TH1D* hAfbLinearityAvgBiasFracStat_rebinnedA_2c = (TH1D*) hAfbLinearityAvgBiasFracStat_rebinnedA_2->Clone();
	TH1D* hAfbLinearityAvgBiasSlope_rebinnedA_2c = (TH1D*) hAfbLinearityAvgBiasSlope_rebinnedA_2->Clone();
	TH1D* hAfbLinearityMaxBiasFrac_rebinnedB_2c = (TH1D*) hAfbLinearityMaxBiasFrac_rebinnedB_2->Clone();
	TH1D* hAfbLinearityMaxBiasFracStat_rebinnedB_2c = (TH1D*) hAfbLinearityMaxBiasFracStat_rebinnedB_2->Clone();
	TH1D* hAfbLinearityMaxBiasSlope_rebinnedB_2c = (TH1D*) hAfbLinearityMaxBiasSlope_rebinnedB_2->Clone();
	TH1D* hAfbLinearityMaxBiasFracMax_rebinnedB_2c = (TH1D*) hAfbLinearityMaxBiasFracMax_rebinnedB_2->Clone();
	TH1D* hAfbLinearityMaxBiasFracStatMax_rebinnedB_2c = (TH1D*) hAfbLinearityMaxBiasFracStatMax_rebinnedB_2->Clone();
	TH1D* hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c = (TH1D*) hAfbLinearityMaxBiasSlopeMax_rebinnedB_2->Clone();
	TH1D* hAfbLinearityAvgBiasFrac_rebinnedB_2c = (TH1D*) hAfbLinearityAvgBiasFrac_rebinnedB_2->Clone();
	TH1D* hAfbLinearityAvgBiasFracStat_rebinnedB_2c = (TH1D*) hAfbLinearityAvgBiasFracStat_rebinnedB_2->Clone();
	TH1D* hAfbLinearityAvgBiasSlope_rebinnedB_2c = (TH1D*) hAfbLinearityAvgBiasSlope_rebinnedB_2->Clone();

	hAfbLinearityMaxBiasFrac_rebinnedA_2c->Add(hAfbLinearityMaxBiasFrac_rebinnedA_1,-1.);
	hAfbLinearityMaxBiasFracStat_rebinnedA_2c->Add(hAfbLinearityMaxBiasFracStat_rebinnedA_1,-1.);
	hAfbLinearityMaxBiasSlope_rebinnedA_2c->Add(hAfbLinearityMaxBiasSlope_rebinnedA_1,-1.);
	hAfbLinearityMaxBiasFracMax_rebinnedA_2c->Add(hAfbLinearityMaxBiasFracMax_rebinnedA_1,-1.);
	hAfbLinearityMaxBiasFracStatMax_rebinnedA_2c->Add(hAfbLinearityMaxBiasFracStatMax_rebinnedA_1,-1.);
	hAfbLinearityMaxBiasSlopeMax_rebinnedA_2c->Add(hAfbLinearityMaxBiasSlopeMax_rebinnedA_1,-1.);
	hAfbLinearityAvgBiasFrac_rebinnedA_2c->Add(hAfbLinearityAvgBiasFrac_rebinnedA_1,-1.);
	hAfbLinearityAvgBiasFracStat_rebinnedA_2c->Add(hAfbLinearityAvgBiasFracStat_rebinnedA_1,-1.);
	hAfbLinearityAvgBiasSlope_rebinnedA_2c->Add(hAfbLinearityAvgBiasSlope_rebinnedA_1,-1.);
	hAfbLinearityMaxBiasFrac_rebinnedB_2c->Add(hAfbLinearityMaxBiasFrac_rebinnedB_1,-1.);
	hAfbLinearityMaxBiasFracStat_rebinnedB_2c->Add(hAfbLinearityMaxBiasFracStat_rebinnedB_1,-1.);
	hAfbLinearityMaxBiasSlope_rebinnedB_2c->Add(hAfbLinearityMaxBiasSlope_rebinnedB_1,-1.);
	hAfbLinearityMaxBiasFracMax_rebinnedB_2c->Add(hAfbLinearityMaxBiasFracMax_rebinnedB_1,-1.);
	hAfbLinearityMaxBiasFracStatMax_rebinnedB_2c->Add(hAfbLinearityMaxBiasFracStatMax_rebinnedB_1,-1.);
	hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->Add(hAfbLinearityMaxBiasSlopeMax_rebinnedB_1,-1.);
	hAfbLinearityAvgBiasFrac_rebinnedB_2c->Add(hAfbLinearityAvgBiasFrac_rebinnedB_1,-1.);
	hAfbLinearityAvgBiasFracStat_rebinnedB_2c->Add(hAfbLinearityAvgBiasFracStat_rebinnedB_1,-1.);
	hAfbLinearityAvgBiasSlope_rebinnedB_2c->Add(hAfbLinearityAvgBiasSlope_rebinnedB_1,-1.);


  TH1D* hAfbLinearityMaxBiasFrac_rebinnedA_2c_residual = (TH1D*) hAfbLinearityMaxBiasFrac_rebinnedA_2c->Clone();
  TH1D* hAfbLinearityMaxBiasFracStat_rebinnedA_2c_residual = (TH1D*) hAfbLinearityMaxBiasFracStat_rebinnedA_2c->Clone();
  TH1D* hAfbLinearityMaxBiasSlope_rebinnedA_2c_residual = (TH1D*) hAfbLinearityMaxBiasSlope_rebinnedA_2c->Clone();
  TH1D* hAfbLinearityMaxBiasFracMax_rebinnedA_2c_residual = (TH1D*) hAfbLinearityMaxBiasFracMax_rebinnedA_2c->Clone();
  TH1D* hAfbLinearityMaxBiasFracStatMax_rebinnedA_2c_residual = (TH1D*) hAfbLinearityMaxBiasFracStatMax_rebinnedA_2c->Clone();
  TH1D* hAfbLinearityMaxBiasSlopeMax_rebinnedA_2c_residual = (TH1D*) hAfbLinearityMaxBiasSlopeMax_rebinnedA_2c->Clone();
  TH1D* hAfbLinearityAvgBiasFrac_rebinnedA_2c_residual = (TH1D*) hAfbLinearityAvgBiasFrac_rebinnedA_2c->Clone();
  TH1D* hAfbLinearityAvgBiasFracStat_rebinnedA_2c_residual = (TH1D*) hAfbLinearityAvgBiasFracStat_rebinnedA_2c->Clone();
  TH1D* hAfbLinearityAvgBiasSlope_rebinnedA_2c_residual = (TH1D*) hAfbLinearityAvgBiasSlope_rebinnedA_2c->Clone();
  TH1D* hAfbLinearityMaxBiasFrac_rebinnedB_2c_residual = (TH1D*) hAfbLinearityMaxBiasFrac_rebinnedB_2c->Clone();
  TH1D* hAfbLinearityMaxBiasFracStat_rebinnedB_2c_residual = (TH1D*) hAfbLinearityMaxBiasFracStat_rebinnedB_2c->Clone();
  TH1D* hAfbLinearityMaxBiasSlope_rebinnedB_2c_residual = (TH1D*) hAfbLinearityMaxBiasSlope_rebinnedB_2c->Clone();
  TH1D* hAfbLinearityMaxBiasFracMax_rebinnedB_2c_residual = (TH1D*) hAfbLinearityMaxBiasFracMax_rebinnedB_2c->Clone();
  TH1D* hAfbLinearityMaxBiasFracStatMax_rebinnedB_2c_residual = (TH1D*) hAfbLinearityMaxBiasFracStatMax_rebinnedB_2c->Clone();
  TH1D* hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c_residual = (TH1D*) hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->Clone();
  TH1D* hAfbLinearityAvgBiasFrac_rebinnedB_2c_residual = (TH1D*) hAfbLinearityAvgBiasFrac_rebinnedB_2c->Clone();
  TH1D* hAfbLinearityAvgBiasFracStat_rebinnedB_2c_residual = (TH1D*) hAfbLinearityAvgBiasFracStat_rebinnedB_2c->Clone();
  TH1D* hAfbLinearityAvgBiasSlope_rebinnedB_2c_residual = (TH1D*) hAfbLinearityAvgBiasSlope_rebinnedB_2c->Clone();

  for (int i_bin = 1; i_bin <= hAfbLinearityMaxBiasFrac_rebinnedA_2c->GetNbinsX(); ++i_bin)
  {
    hAfbLinearityMaxBiasFrac_rebinnedA_2c_residual->SetBinContent(i_bin,hAfbLinearityMaxBiasFrac_rebinnedA_2c->GetBinContent(i_bin)+1);
    hAfbLinearityMaxBiasFracStat_rebinnedA_2c_residual->SetBinContent(i_bin,hAfbLinearityMaxBiasFracStat_rebinnedA_2c->GetBinContent(i_bin)+1);
    hAfbLinearityMaxBiasSlope_rebinnedA_2c_residual->SetBinContent(i_bin,hAfbLinearityMaxBiasSlope_rebinnedA_2c->GetBinContent(i_bin)+1);
    hAfbLinearityMaxBiasFracMax_rebinnedA_2c_residual->SetBinContent(i_bin,hAfbLinearityMaxBiasFracMax_rebinnedA_2c->GetBinContent(i_bin)+1);
    hAfbLinearityMaxBiasFracStatMax_rebinnedA_2c_residual->SetBinContent(i_bin,hAfbLinearityMaxBiasFracStatMax_rebinnedA_2c->GetBinContent(i_bin)+1);
    hAfbLinearityMaxBiasSlopeMax_rebinnedA_2c_residual->SetBinContent(i_bin,hAfbLinearityMaxBiasSlopeMax_rebinnedA_2c->GetBinContent(i_bin)+1);
    hAfbLinearityAvgBiasFrac_rebinnedA_2c_residual->SetBinContent(i_bin,hAfbLinearityAvgBiasFrac_rebinnedA_2c->GetBinContent(i_bin)+1);
    hAfbLinearityAvgBiasFracStat_rebinnedA_2c_residual->SetBinContent(i_bin,hAfbLinearityAvgBiasFracStat_rebinnedA_2c->GetBinContent(i_bin)+1);
    hAfbLinearityAvgBiasSlope_rebinnedA_2c_residual->SetBinContent(i_bin,hAfbLinearityAvgBiasSlope_rebinnedA_2c->GetBinContent(i_bin)+1);
    hAfbLinearityMaxBiasFrac_rebinnedB_2c_residual->SetBinContent(i_bin,hAfbLinearityMaxBiasFrac_rebinnedB_2c->GetBinContent(i_bin)+1);
    hAfbLinearityMaxBiasFracStat_rebinnedB_2c_residual->SetBinContent(i_bin,hAfbLinearityMaxBiasFracStat_rebinnedB_2c->GetBinContent(i_bin)+1);
    hAfbLinearityMaxBiasSlope_rebinnedB_2c_residual->SetBinContent(i_bin,hAfbLinearityMaxBiasSlope_rebinnedB_2c->GetBinContent(i_bin)+1);
    hAfbLinearityMaxBiasFracMax_rebinnedB_2c_residual->SetBinContent(i_bin,hAfbLinearityMaxBiasFracMax_rebinnedB_2c->GetBinContent(i_bin)+1);
    hAfbLinearityMaxBiasFracStatMax_rebinnedB_2c_residual->SetBinContent(i_bin,hAfbLinearityMaxBiasFracStatMax_rebinnedB_2c->GetBinContent(i_bin)+1);
    hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c_residual->SetBinContent(i_bin,hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->GetBinContent(i_bin)+1);
    hAfbLinearityAvgBiasFrac_rebinnedB_2c_residual->SetBinContent(i_bin,hAfbLinearityAvgBiasFrac_rebinnedB_2c->GetBinContent(i_bin)+1);
    hAfbLinearityAvgBiasFracStat_rebinnedB_2c_residual->SetBinContent(i_bin,hAfbLinearityAvgBiasFracStat_rebinnedB_2c->GetBinContent(i_bin)+1);
    hAfbLinearityAvgBiasSlope_rebinnedB_2c_residual->SetBinContent(i_bin,hAfbLinearityAvgBiasSlope_rebinnedB_2c->GetBinContent(i_bin)+1);
  }

      TCanvas *c2Afbr = new TCanvas("c2Afbr", "RegBias",1024,1280);
      gStyle->SetPadLeftMargin(0.08);

      c2Afbr->Divide(1,3);

      c2Afbr->cd(1);
      gPad->SetLogy(0);
      hAfbLinearityMaxBiasSlopeMax_rebinnedB_1->SetLineColor(2);
      hAfbLinearityMaxBiasSlopeMax_rebinnedB_1->GetXaxis()->SetTitle("Variable");
      hAfbLinearityMaxBiasSlopeMax_rebinnedB_1->GetYaxis()->SetTitle("Measured coefficient linearity slope");
      if(hAfbLinearityMaxBiasSlopeMax_rebinnedB_1->GetMaximum()<hAfbLinearityMaxBiasSlopeMax_rebinnedA_1->GetMaximum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB_1->SetMaximum(hAfbLinearityMaxBiasSlopeMax_rebinnedA_1->GetMaximum());
      if(hAfbLinearityMaxBiasSlopeMax_rebinnedB_1->GetMinimum()>hAfbLinearityMaxBiasSlopeMax_rebinnedA_1->GetMinimum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB_1->SetMinimum(hAfbLinearityMaxBiasSlopeMax_rebinnedA_1->GetMinimum());
      if(hAfbLinearityMaxBiasSlopeMax_rebinnedB_1->GetMaximum()<hAfbLinearityAvgBiasSlope_rebinnedA_1->GetMaximum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB_1->SetMaximum(hAfbLinearityAvgBiasSlope_rebinnedA_1->GetMaximum());
      if(hAfbLinearityMaxBiasSlopeMax_rebinnedB_1->GetMinimum()>hAfbLinearityAvgBiasSlope_rebinnedA_1->GetMinimum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB_1->SetMinimum(hAfbLinearityAvgBiasSlope_rebinnedA_1->GetMinimum());
      if(hAfbLinearityMaxBiasSlopeMax_rebinnedB_1->GetMaximum()<hAfbLinearityAvgBiasSlope_rebinnedB_1->GetMaximum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB_1->SetMaximum(hAfbLinearityAvgBiasSlope_rebinnedB_1->GetMaximum());
      if(hAfbLinearityMaxBiasSlopeMax_rebinnedB_1->GetMinimum()>hAfbLinearityAvgBiasSlope_rebinnedB_1->GetMinimum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB_1->SetMinimum(hAfbLinearityAvgBiasSlope_rebinnedB_1->GetMinimum());
      hAfbLinearityMaxBiasSlopeMax_rebinnedB_1->SetMinimum(0.99);
      hAfbLinearityMaxBiasSlopeMax_rebinnedB_1->Draw();
      hAfbLinearityMaxBiasSlope_rebinnedB_1->Draw("sames");
      hAfbLinearityAvgBiasSlope_rebinnedB_1->SetLineColor(3);
      hAfbLinearityAvgBiasSlope_rebinnedB_1->Draw("sames");
      hAfbLinearityMaxBiasSlopeMax_rebinnedA_1->SetLineColor(2);
      hAfbLinearityMaxBiasSlopeMax_rebinnedA_1->SetLineStyle(2);
      hAfbLinearityMaxBiasSlopeMax_rebinnedA_1->Draw("sames");
      hAfbLinearityMaxBiasSlope_rebinnedA_1->SetLineStyle(2);
      hAfbLinearityMaxBiasSlope_rebinnedA_1->Draw("sames");
      hAfbLinearityAvgBiasSlope_rebinnedA_1->SetLineColor(3);
      hAfbLinearityAvgBiasSlope_rebinnedA_1->SetLineStyle(2);
      hAfbLinearityAvgBiasSlope_rebinnedA_1->Draw("sames");


      c2Afbr->cd(2);
      gPad->SetLogy(0);
      hAfbLinearityMaxBiasSlopeMax_rebinnedB_2->SetLineColor(2);
      hAfbLinearityMaxBiasSlopeMax_rebinnedB_2->GetXaxis()->SetTitle("Variable");
      hAfbLinearityMaxBiasSlopeMax_rebinnedB_2->GetYaxis()->SetTitle("Measured coefficient linearity slope");
      if(hAfbLinearityMaxBiasSlopeMax_rebinnedB_2->GetMaximum()<hAfbLinearityMaxBiasSlopeMax_rebinnedA_2->GetMaximum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB_2->SetMaximum(hAfbLinearityMaxBiasSlopeMax_rebinnedA_2->GetMaximum());
      if(hAfbLinearityMaxBiasSlopeMax_rebinnedB_2->GetMinimum()>hAfbLinearityMaxBiasSlopeMax_rebinnedA_2->GetMinimum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB_2->SetMinimum(hAfbLinearityMaxBiasSlopeMax_rebinnedA_2->GetMinimum());
      if(hAfbLinearityMaxBiasSlopeMax_rebinnedB_2->GetMaximum()<hAfbLinearityAvgBiasSlope_rebinnedA_2->GetMaximum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB_2->SetMaximum(hAfbLinearityAvgBiasSlope_rebinnedA_2->GetMaximum());
      if(hAfbLinearityMaxBiasSlopeMax_rebinnedB_2->GetMinimum()>hAfbLinearityAvgBiasSlope_rebinnedA_2->GetMinimum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB_2->SetMinimum(hAfbLinearityAvgBiasSlope_rebinnedA_2->GetMinimum());
      if(hAfbLinearityMaxBiasSlopeMax_rebinnedB_2->GetMaximum()<hAfbLinearityAvgBiasSlope_rebinnedB_2->GetMaximum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB_2->SetMaximum(hAfbLinearityAvgBiasSlope_rebinnedB_2->GetMaximum());
      if(hAfbLinearityMaxBiasSlopeMax_rebinnedB_2->GetMinimum()>hAfbLinearityAvgBiasSlope_rebinnedB_2->GetMinimum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB_2->SetMinimum(hAfbLinearityAvgBiasSlope_rebinnedB_2->GetMinimum());
      hAfbLinearityMaxBiasSlopeMax_rebinnedB_2->Draw();
      hAfbLinearityMaxBiasSlope_rebinnedB_2->Draw("sames");
      hAfbLinearityAvgBiasSlope_rebinnedB_2->SetLineColor(3);
      hAfbLinearityAvgBiasSlope_rebinnedB_2->Draw("sames");
      hAfbLinearityMaxBiasSlopeMax_rebinnedA_2->SetLineColor(2);
      hAfbLinearityMaxBiasSlopeMax_rebinnedA_2->SetLineStyle(2);
      hAfbLinearityMaxBiasSlopeMax_rebinnedA_2->Draw("sames");
      hAfbLinearityMaxBiasSlope_rebinnedA_2->SetLineStyle(2);
      hAfbLinearityMaxBiasSlope_rebinnedA_2->Draw("sames");
      hAfbLinearityAvgBiasSlope_rebinnedA_2->SetLineColor(3);
      hAfbLinearityAvgBiasSlope_rebinnedA_2->SetLineStyle(2);
      hAfbLinearityAvgBiasSlope_rebinnedA_2->Draw("sames");



      c2Afbr->cd(3);
      gPad->SetLogy(0);
      hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->SetLineColor(2);
      hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->GetXaxis()->SetTitle("Variable");
      hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->GetYaxis()->SetTitle("Measured coefficient linearity slope");
      if(hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->GetMaximum()<hAfbLinearityMaxBiasSlopeMax_rebinnedA_2c->GetMaximum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->SetMaximum(hAfbLinearityMaxBiasSlopeMax_rebinnedA_2c->GetMaximum());
      if(hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->GetMinimum()>hAfbLinearityMaxBiasSlopeMax_rebinnedA_2c->GetMinimum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->SetMinimum(hAfbLinearityMaxBiasSlopeMax_rebinnedA_2c->GetMinimum());
      if(hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->GetMaximum()<hAfbLinearityAvgBiasSlope_rebinnedA_2c->GetMaximum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->SetMaximum(hAfbLinearityAvgBiasSlope_rebinnedA_2c->GetMaximum());
      if(hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->GetMinimum()>hAfbLinearityAvgBiasSlope_rebinnedA_2c->GetMinimum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->SetMinimum(hAfbLinearityAvgBiasSlope_rebinnedA_2c->GetMinimum());
      if(hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->GetMaximum()<hAfbLinearityAvgBiasSlope_rebinnedB_2c->GetMaximum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->SetMaximum(hAfbLinearityAvgBiasSlope_rebinnedB_2c->GetMaximum());
      if(hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->GetMinimum()>hAfbLinearityAvgBiasSlope_rebinnedB_2c->GetMinimum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->SetMinimum(hAfbLinearityAvgBiasSlope_rebinnedB_2c->GetMinimum());
      hAfbLinearityMaxBiasSlopeMax_rebinnedB_2c->Draw();
      hAfbLinearityMaxBiasSlope_rebinnedB_2c->Draw("sames");
      hAfbLinearityAvgBiasSlope_rebinnedB_2c->SetLineColor(3);
      hAfbLinearityAvgBiasSlope_rebinnedB_2c->Draw("sames");
      hAfbLinearityMaxBiasSlopeMax_rebinnedA_2c->SetLineColor(2);
      hAfbLinearityMaxBiasSlopeMax_rebinnedA_2c->SetLineStyle(2);
      hAfbLinearityMaxBiasSlopeMax_rebinnedA_2c->Draw("sames");
      hAfbLinearityMaxBiasSlope_rebinnedA_2c->SetLineStyle(2);
      hAfbLinearityMaxBiasSlope_rebinnedA_2c->Draw("sames");
      hAfbLinearityAvgBiasSlope_rebinnedA_2c->SetLineColor(3);
      hAfbLinearityAvgBiasSlope_rebinnedA_2c->SetLineStyle(2);
      hAfbLinearityAvgBiasSlope_rebinnedA_2c->Draw("sames");


      c2Afbr->Print("Linearity_slope_summary.pdf");




      TCanvas *c2Afbr_residual = new TCanvas("c2Afbr_residual", "RegBias_residual",1024,500);
      hAfbLinearityMaxBiasSlope_rebinnedB_2->SetMinimum(0.99);
      hAfbLinearityMaxBiasSlope_rebinnedB_2->GetXaxis()->SetTitle("Variable");
      hAfbLinearityMaxBiasSlope_rebinnedB_2->GetYaxis()->SetTitle("Measured coefficient linearity slope");
      hAfbLinearityMaxBiasSlope_rebinnedB_2->Draw();
      hAfbLinearityMaxBiasSlope_rebinnedB_2c_residual->SetLineStyle(2);
      hAfbLinearityMaxBiasSlope_rebinnedB_2c_residual->Draw("same");
      c2Afbr_residual->Print("Linearity_slope_summary_residual.pdf");





  TH1D* hAfbConstantMaxBiasSlope_rebinnedA_1 = (TH1D*)File1->Get("hAfbConstantMaxBiasSlope_rebinnedA");
  TH1D* hAfbConstantAvgBiasSlope_rebinnedA_1 = (TH1D*)File1->Get("hAfbConstantAvgBiasSlope_rebinnedA");
  TH1D* hAfbConstantMaxBiasSlope_rebinnedB_1 = (TH1D*)File1->Get("hAfbConstantMaxBiasSlope_rebinnedB");
  TH1D* hAfbConstantAvgBiasSlope_rebinnedB_1 = (TH1D*)File1->Get("hAfbConstantAvgBiasSlope_rebinnedB");

  TH1D* hAfbConstantMaxBiasSlope_rebinnedA_2 = (TH1D*)File2->Get("hAfbConstantMaxBiasSlope_rebinnedA");
  TH1D* hAfbConstantAvgBiasSlope_rebinnedA_2 = (TH1D*)File2->Get("hAfbConstantAvgBiasSlope_rebinnedA");
  TH1D* hAfbConstantMaxBiasSlope_rebinnedB_2 = (TH1D*)File2->Get("hAfbConstantMaxBiasSlope_rebinnedB");
  TH1D* hAfbConstantAvgBiasSlope_rebinnedB_2 = (TH1D*)File2->Get("hAfbConstantAvgBiasSlope_rebinnedB");

  TH1D* hAfbConstantMaxBiasSlope_rebinnedA_2c = (TH1D*) hAfbConstantMaxBiasSlope_rebinnedA_2->Clone();
  TH1D* hAfbConstantAvgBiasSlope_rebinnedA_2c = (TH1D*) hAfbConstantAvgBiasSlope_rebinnedA_2->Clone();
  TH1D* hAfbConstantMaxBiasSlope_rebinnedB_2c = (TH1D*) hAfbConstantMaxBiasSlope_rebinnedB_2->Clone();
  TH1D* hAfbConstantAvgBiasSlope_rebinnedB_2c = (TH1D*) hAfbConstantAvgBiasSlope_rebinnedB_2->Clone();

  hAfbConstantMaxBiasSlope_rebinnedA_2c->Add(hAfbConstantMaxBiasSlope_rebinnedA_1,-1.);
  hAfbConstantAvgBiasSlope_rebinnedA_2c->Add(hAfbConstantAvgBiasSlope_rebinnedA_1,-1.);
  hAfbConstantMaxBiasSlope_rebinnedB_2c->Add(hAfbConstantMaxBiasSlope_rebinnedB_1,-1.);
  hAfbConstantAvgBiasSlope_rebinnedB_2c->Add(hAfbConstantAvgBiasSlope_rebinnedB_1,-1.);




      TCanvas *c2Afbrc = new TCanvas("c2Afbrc", "RegBias",1024,1280);
      gStyle->SetPadLeftMargin(0.08);

      c2Afbrc->Divide(1,3);

      c2Afbrc->cd(1);
      gPad->SetLogy(0);
      hAfbConstantMaxBiasSlope_rebinnedB_1->SetLineColor(2);
      hAfbConstantMaxBiasSlope_rebinnedB_1->GetXaxis()->SetTitle("Variable");
      hAfbConstantMaxBiasSlope_rebinnedB_1->GetYaxis()->SetTitle("Measured coefficient bias (offset)");
      if(hAfbConstantMaxBiasSlope_rebinnedB_1->GetMaximum()<hAfbConstantMaxBiasSlope_rebinnedA_1->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedB_1->SetMaximum(hAfbConstantMaxBiasSlope_rebinnedA_1->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_1->GetMinimum()>hAfbConstantMaxBiasSlope_rebinnedA_1->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedB_1->SetMinimum(hAfbConstantMaxBiasSlope_rebinnedA_1->GetMinimum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_1->GetMaximum()<hAfbConstantAvgBiasSlope_rebinnedA_1->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedB_1->SetMaximum(hAfbConstantAvgBiasSlope_rebinnedA_1->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_1->GetMinimum()>hAfbConstantAvgBiasSlope_rebinnedA_1->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedB_1->SetMinimum(hAfbConstantAvgBiasSlope_rebinnedA_1->GetMinimum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_1->GetMaximum()<hAfbConstantAvgBiasSlope_rebinnedB_1->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedB_1->SetMaximum(hAfbConstantAvgBiasSlope_rebinnedB_1->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_1->GetMinimum()>hAfbConstantAvgBiasSlope_rebinnedB_1->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedB_1->SetMinimum(hAfbConstantAvgBiasSlope_rebinnedB_1->GetMinimum());
      hAfbConstantMaxBiasSlope_rebinnedB_1->Draw();
      hAfbConstantAvgBiasSlope_rebinnedB_1->SetLineColor(3);
      hAfbConstantAvgBiasSlope_rebinnedB_1->Draw("sames");
      hAfbConstantMaxBiasSlope_rebinnedA_1->SetLineColor(2);
      hAfbConstantMaxBiasSlope_rebinnedA_1->SetLineStyle(2);
      hAfbConstantMaxBiasSlope_rebinnedA_1->Draw("sames");
      hAfbConstantAvgBiasSlope_rebinnedA_1->SetLineColor(3);
      hAfbConstantAvgBiasSlope_rebinnedA_1->SetLineStyle(2);
      hAfbConstantAvgBiasSlope_rebinnedA_1->Draw("sames");


      c2Afbrc->cd(2);
      gPad->SetLogy(0);
      hAfbConstantMaxBiasSlope_rebinnedB_2->SetLineColor(2);
      hAfbConstantMaxBiasSlope_rebinnedB_2->GetXaxis()->SetTitle("Variable");
      hAfbConstantMaxBiasSlope_rebinnedB_2->GetYaxis()->SetTitle("Measured coefficient bias (offset)");
      if(hAfbConstantMaxBiasSlope_rebinnedB_2->GetMaximum()<hAfbConstantMaxBiasSlope_rebinnedA_2->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedB_2->SetMaximum(hAfbConstantMaxBiasSlope_rebinnedA_2->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_2->GetMinimum()>hAfbConstantMaxBiasSlope_rebinnedA_2->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedB_2->SetMinimum(hAfbConstantMaxBiasSlope_rebinnedA_2->GetMinimum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_2->GetMaximum()<hAfbConstantAvgBiasSlope_rebinnedA_2->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedB_2->SetMaximum(hAfbConstantAvgBiasSlope_rebinnedA_2->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_2->GetMinimum()>hAfbConstantAvgBiasSlope_rebinnedA_2->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedB_2->SetMinimum(hAfbConstantAvgBiasSlope_rebinnedA_2->GetMinimum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_2->GetMaximum()<hAfbConstantAvgBiasSlope_rebinnedB_2->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedB_2->SetMaximum(hAfbConstantAvgBiasSlope_rebinnedB_2->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_2->GetMinimum()>hAfbConstantAvgBiasSlope_rebinnedB_2->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedB_2->SetMinimum(hAfbConstantAvgBiasSlope_rebinnedB_2->GetMinimum());
      hAfbConstantMaxBiasSlope_rebinnedB_2->Draw();
      hAfbConstantAvgBiasSlope_rebinnedB_2->SetLineColor(3);
      hAfbConstantAvgBiasSlope_rebinnedB_2->Draw("sames");
      hAfbConstantMaxBiasSlope_rebinnedA_2->SetLineColor(2);
      hAfbConstantMaxBiasSlope_rebinnedA_2->SetLineStyle(2);
      hAfbConstantMaxBiasSlope_rebinnedA_2->Draw("sames");
      hAfbConstantAvgBiasSlope_rebinnedA_2->SetLineColor(3);
      hAfbConstantAvgBiasSlope_rebinnedA_2->SetLineStyle(2);
      hAfbConstantAvgBiasSlope_rebinnedA_2->Draw("sames");



      c2Afbrc->cd(3);
      gPad->SetLogy(0);
      hAfbConstantMaxBiasSlope_rebinnedB_2c->SetLineColor(2);
      hAfbConstantMaxBiasSlope_rebinnedB_2c->GetXaxis()->SetTitle("Variable");
      hAfbConstantMaxBiasSlope_rebinnedB_2c->GetYaxis()->SetTitle("Measured coefficient bias (offset)");
      if(hAfbConstantMaxBiasSlope_rebinnedB_2c->GetMaximum()<hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedB_2c->SetMaximum(hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_2c->GetMinimum()>hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedB_2c->SetMinimum(hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMinimum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_2c->GetMaximum()<hAfbConstantAvgBiasSlope_rebinnedA_2c->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedB_2c->SetMaximum(hAfbConstantAvgBiasSlope_rebinnedA_2c->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_2c->GetMinimum()>hAfbConstantAvgBiasSlope_rebinnedA_2c->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedB_2c->SetMinimum(hAfbConstantAvgBiasSlope_rebinnedA_2c->GetMinimum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_2c->GetMaximum()<hAfbConstantAvgBiasSlope_rebinnedB_2c->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedB_2c->SetMaximum(hAfbConstantAvgBiasSlope_rebinnedB_2c->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_2c->GetMinimum()>hAfbConstantAvgBiasSlope_rebinnedB_2c->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedB_2c->SetMinimum(hAfbConstantAvgBiasSlope_rebinnedB_2c->GetMinimum());
      hAfbConstantMaxBiasSlope_rebinnedB_2c->Draw();
      hAfbConstantAvgBiasSlope_rebinnedB_2c->SetLineColor(3);
      hAfbConstantAvgBiasSlope_rebinnedB_2c->Draw("sames");
      hAfbConstantMaxBiasSlope_rebinnedA_2c->SetLineColor(2);
      hAfbConstantMaxBiasSlope_rebinnedA_2c->SetLineStyle(2);
      hAfbConstantMaxBiasSlope_rebinnedA_2c->Draw("sames");
      hAfbConstantAvgBiasSlope_rebinnedA_2c->SetLineColor(3);
      hAfbConstantAvgBiasSlope_rebinnedA_2c->SetLineStyle(2);
      hAfbConstantAvgBiasSlope_rebinnedA_2c->Draw("sames");


      c2Afbrc->Print("Linearity_constant_summary.pdf");




      double dC = 0.1;
      hAfbLinearityMaxBiasSlope_rebinnedB_2c->Scale(dC);
      hAfbLinearityAvgBiasSlope_rebinnedB_2c->Scale(dC);
      hAfbLinearityMaxBiasSlope_rebinnedA_2c->Scale(dC);
      hAfbLinearityAvgBiasSlope_rebinnedA_2c->Scale(dC);



      TCanvas *c2Afbrc2 = new TCanvas("c2Afbrc2", "RegBias",1024,1280);
      gStyle->SetPadLeftMargin(0.10);

      c2Afbrc2->Divide(1,3);

      c2Afbrc2->cd(1);
      gPad->SetLogy(0);
      hAfbConstantMaxBiasSlope_rebinnedB_2c->SetLineColor(2);
      hAfbConstantMaxBiasSlope_rebinnedB_2c->GetXaxis()->SetTitle("Variable");
      hAfbConstantMaxBiasSlope_rebinnedB_2c->GetYaxis()->SetTitle("Measured coefficient bias (#Deltaoffset)");
      if(hAfbConstantMaxBiasSlope_rebinnedB_2c->GetMaximum()<hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedB_2c->SetMaximum(hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_2c->GetMinimum()>hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedB_2c->SetMinimum(hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMinimum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_2c->GetMaximum()<hAfbConstantAvgBiasSlope_rebinnedA_2c->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedB_2c->SetMaximum(hAfbConstantAvgBiasSlope_rebinnedA_2c->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_2c->GetMinimum()>hAfbConstantAvgBiasSlope_rebinnedA_2c->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedB_2c->SetMinimum(hAfbConstantAvgBiasSlope_rebinnedA_2c->GetMinimum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_2c->GetMaximum()<hAfbConstantAvgBiasSlope_rebinnedB_2c->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedB_2c->SetMaximum(hAfbConstantAvgBiasSlope_rebinnedB_2c->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedB_2c->GetMinimum()>hAfbConstantAvgBiasSlope_rebinnedB_2c->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedB_2c->SetMinimum(hAfbConstantAvgBiasSlope_rebinnedB_2c->GetMinimum());
      hAfbConstantMaxBiasSlope_rebinnedB_2c->Draw();
      hAfbConstantAvgBiasSlope_rebinnedB_2c->SetLineColor(3);
      hAfbConstantAvgBiasSlope_rebinnedB_2c->Draw("sames");
      hAfbConstantMaxBiasSlope_rebinnedA_2c->SetLineColor(2);
      hAfbConstantMaxBiasSlope_rebinnedA_2c->SetLineStyle(2);
      hAfbConstantMaxBiasSlope_rebinnedA_2c->Draw("sames");
      hAfbConstantAvgBiasSlope_rebinnedA_2c->SetLineColor(3);
      hAfbConstantAvgBiasSlope_rebinnedA_2c->SetLineStyle(2);
      hAfbConstantAvgBiasSlope_rebinnedA_2c->Draw("sames");



      c2Afbrc2->cd(2);
      gPad->SetLogy(0);
      hAfbLinearityMaxBiasSlope_rebinnedB_2c->SetLineColor(2);
      hAfbLinearityMaxBiasSlope_rebinnedB_2c->GetXaxis()->SetTitle("Variable");
      hAfbLinearityMaxBiasSlope_rebinnedB_2c->GetYaxis()->SetTitle("Measured coefficient bias (0.1 #Deltaslope)");
      if(hAfbLinearityMaxBiasSlope_rebinnedB_2c->GetMaximum()<hAfbLinearityMaxBiasSlope_rebinnedA_2c->GetMaximum()) hAfbLinearityMaxBiasSlope_rebinnedB_2c->SetMaximum(hAfbLinearityMaxBiasSlope_rebinnedA_2c->GetMaximum());
      if(hAfbLinearityMaxBiasSlope_rebinnedB_2c->GetMinimum()>hAfbLinearityMaxBiasSlope_rebinnedA_2c->GetMinimum()) hAfbLinearityMaxBiasSlope_rebinnedB_2c->SetMinimum(hAfbLinearityMaxBiasSlope_rebinnedA_2c->GetMinimum());
      if(hAfbLinearityMaxBiasSlope_rebinnedB_2c->GetMaximum()<hAfbLinearityAvgBiasSlope_rebinnedA_2c->GetMaximum()) hAfbLinearityMaxBiasSlope_rebinnedB_2c->SetMaximum(hAfbLinearityAvgBiasSlope_rebinnedA_2c->GetMaximum());
      if(hAfbLinearityMaxBiasSlope_rebinnedB_2c->GetMinimum()>hAfbLinearityAvgBiasSlope_rebinnedA_2c->GetMinimum()) hAfbLinearityMaxBiasSlope_rebinnedB_2c->SetMinimum(hAfbLinearityAvgBiasSlope_rebinnedA_2c->GetMinimum());
      if(hAfbLinearityMaxBiasSlope_rebinnedB_2c->GetMaximum()<hAfbLinearityAvgBiasSlope_rebinnedB_2c->GetMaximum()) hAfbLinearityMaxBiasSlope_rebinnedB_2c->SetMaximum(hAfbLinearityAvgBiasSlope_rebinnedB_2c->GetMaximum());
      if(hAfbLinearityMaxBiasSlope_rebinnedB_2c->GetMinimum()>hAfbLinearityAvgBiasSlope_rebinnedB_2c->GetMinimum()) hAfbLinearityMaxBiasSlope_rebinnedB_2c->SetMinimum(hAfbLinearityAvgBiasSlope_rebinnedB_2c->GetMinimum());
      hAfbLinearityMaxBiasSlope_rebinnedB_2c->Draw("hist");
      hAfbLinearityAvgBiasSlope_rebinnedB_2c->SetLineColor(3);
      hAfbLinearityAvgBiasSlope_rebinnedB_2c->Draw("hist sames");
      hAfbLinearityMaxBiasSlope_rebinnedA_2c->SetLineColor(2);
      hAfbLinearityMaxBiasSlope_rebinnedA_2c->SetLineStyle(2);
      hAfbLinearityMaxBiasSlope_rebinnedA_2c->Draw("hist sames");
      hAfbLinearityAvgBiasSlope_rebinnedA_2c->SetLineColor(3);
      hAfbLinearityAvgBiasSlope_rebinnedA_2c->SetLineStyle(2);
      hAfbLinearityAvgBiasSlope_rebinnedA_2c->Draw("hist sames");




      c2Afbrc2->cd(3);
      gPad->SetLogy(0);

      hAfbConstantMaxBiasSlope_rebinnedA_2c->SetLineColor(2);
      hAfbConstantMaxBiasSlope_rebinnedA_2c->SetLineStyle(2);
      hAfbConstantMaxBiasSlope_rebinnedA_2c->GetXaxis()->SetTitle("Variable");
      hAfbConstantMaxBiasSlope_rebinnedA_2c->GetYaxis()->SetTitle("Measured coefficient bias (#Deltaoffset or 0.1 #Deltaslope)");

      if(hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMaximum()<hAfbConstantMaxBiasSlope_rebinnedB_2c->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedA_2c->SetMaximum(hAfbConstantMaxBiasSlope_rebinnedB_2c->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMinimum()>hAfbConstantMaxBiasSlope_rebinnedB_2c->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedA_2c->SetMinimum(hAfbConstantMaxBiasSlope_rebinnedB_2c->GetMinimum());
      if(hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMaximum()<hAfbConstantAvgBiasSlope_rebinnedA_2c->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedA_2c->SetMaximum(hAfbConstantAvgBiasSlope_rebinnedA_2c->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMinimum()>hAfbConstantAvgBiasSlope_rebinnedA_2c->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedA_2c->SetMinimum(hAfbConstantAvgBiasSlope_rebinnedA_2c->GetMinimum());
      if(hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMaximum()<hAfbConstantAvgBiasSlope_rebinnedB_2c->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedA_2c->SetMaximum(hAfbConstantAvgBiasSlope_rebinnedB_2c->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMinimum()>hAfbConstantAvgBiasSlope_rebinnedB_2c->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedA_2c->SetMinimum(hAfbConstantAvgBiasSlope_rebinnedB_2c->GetMinimum());

      if(hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMaximum()<hAfbLinearityMaxBiasSlope_rebinnedA_2c->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedA_2c->SetMaximum(hAfbLinearityMaxBiasSlope_rebinnedA_2c->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMinimum()>hAfbLinearityMaxBiasSlope_rebinnedA_2c->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedA_2c->SetMinimum(hAfbLinearityMaxBiasSlope_rebinnedA_2c->GetMinimum());
      if(hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMaximum()<hAfbLinearityMaxBiasSlope_rebinnedB_2c->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedA_2c->SetMaximum(hAfbLinearityMaxBiasSlope_rebinnedB_2c->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMinimum()>hAfbLinearityMaxBiasSlope_rebinnedB_2c->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedA_2c->SetMinimum(hAfbLinearityMaxBiasSlope_rebinnedB_2c->GetMinimum());
      if(hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMaximum()<hAfbLinearityAvgBiasSlope_rebinnedA_2c->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedA_2c->SetMaximum(hAfbLinearityAvgBiasSlope_rebinnedA_2c->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMinimum()>hAfbLinearityAvgBiasSlope_rebinnedA_2c->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedA_2c->SetMinimum(hAfbLinearityAvgBiasSlope_rebinnedA_2c->GetMinimum());
      if(hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMaximum()<hAfbLinearityAvgBiasSlope_rebinnedB_2c->GetMaximum()) hAfbConstantMaxBiasSlope_rebinnedA_2c->SetMaximum(hAfbLinearityAvgBiasSlope_rebinnedB_2c->GetMaximum());
      if(hAfbConstantMaxBiasSlope_rebinnedA_2c->GetMinimum()>hAfbLinearityAvgBiasSlope_rebinnedB_2c->GetMinimum()) hAfbConstantMaxBiasSlope_rebinnedA_2c->SetMinimum(hAfbLinearityAvgBiasSlope_rebinnedB_2c->GetMinimum());

      hAfbConstantMaxBiasSlope_rebinnedA_2c->Draw();
      hAfbConstantMaxBiasSlope_rebinnedB_2c->SetLineColor(2);
      hAfbConstantMaxBiasSlope_rebinnedB_2c->Draw("sames");
      hAfbConstantAvgBiasSlope_rebinnedB_2c->SetLineColor(3);
      hAfbConstantAvgBiasSlope_rebinnedB_2c->Draw("sames");
      hAfbConstantAvgBiasSlope_rebinnedA_2c->SetLineColor(3);
      hAfbConstantAvgBiasSlope_rebinnedA_2c->SetLineStyle(2);
      hAfbConstantAvgBiasSlope_rebinnedA_2c->Draw("sames");

      hAfbLinearityMaxBiasSlope_rebinnedB_2c->SetLineColor(kRed-1); 
      hAfbLinearityMaxBiasSlope_rebinnedB_2c->Draw("hist sames"); 
      hAfbLinearityAvgBiasSlope_rebinnedB_2c->SetLineColor(kGreen-1); 
      hAfbLinearityAvgBiasSlope_rebinnedB_2c->Draw("hist sames"); 
      hAfbLinearityMaxBiasSlope_rebinnedA_2c->SetLineColor(kRed-1); 
      hAfbLinearityMaxBiasSlope_rebinnedA_2c->SetLineStyle(2); 
      hAfbLinearityMaxBiasSlope_rebinnedA_2c->Draw("hist sames"); 
      hAfbLinearityAvgBiasSlope_rebinnedA_2c->SetLineColor(kGreen-1); 
      hAfbLinearityAvgBiasSlope_rebinnedA_2c->SetLineStyle(2); 
      hAfbLinearityAvgBiasSlope_rebinnedA_2c->Draw("hist sames"); 


      c2Afbrc2->Print("Linearity_difference_summary.pdf");



}