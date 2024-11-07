#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TH2.h>
#include <TH1.h>
#include <Rtypes.h>
#include <TFile.h>
#include <iostream> 
#include <vector>
#include <TGraphAsymmErrors.h>
#include "utils.h"

#include "kinRecoFullLooseQualityPlots.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/utils.h"
#include "../../common/include/CommandLineParameters.h"

PrepareQualityPlots::PrepareQualityPlots(const char *titleTrueVsDeltaTrueReco, const char *titleProjYTrueVsDeltaTrueReco, const char *titleGenVsRMS, const char *titleGenVsMean):
titleTrueVsDeltaTrueReco_(titleTrueVsDeltaTrueReco),
titleProjYTrueVsDeltaTrueReco_(titleProjYTrueVsDeltaTrueReco),
titleGenVsRMS_(titleGenVsRMS),
titleGenVsMean_(titleGenVsMean),
h_TrueVsDeltaTrueReco_(0),
h_projYTrueVsDeltaTrueReco_(0),
h_GenVsRMS_(0),
h_GenVsMean_(0)
{
}

void PrepareQualityPlots::prepareQualityPlots(const TString& fileDir, const TString& fileName, const TString& histoName, const int& Xnb, const float& Xr1, const float& Xr2)
{
  TH1::AddDirectory(kFALSE);
   
  float dXbin;

  TFile rootFile(fileDir + fileName);
  h_TrueVsDeltaTrueReco_ = (TH2D*)rootFile.Get(histoName);
  h_TrueVsDeltaTrueReco_->SetTitle(titleTrueVsDeltaTrueReco_);
  h_TrueVsDeltaTrueReco_->GetXaxis()->SetTitleOffset(1.20);
  h_TrueVsDeltaTrueReco_->GetYaxis()->SetTitleOffset(1.80);
  h_TrueVsDeltaTrueReco_->SetStats(0);
  h_TrueVsDeltaTrueReco_->SetDirectory(0);

  h_projYTrueVsDeltaTrueReco_= h_TrueVsDeltaTrueReco_->ProjectionY("h_py",1,h_TrueVsDeltaTrueReco_->GetXaxis()->GetNbins(),"");
  h_projYTrueVsDeltaTrueReco_->SetTitle(titleProjYTrueVsDeltaTrueReco_);
  h_projYTrueVsDeltaTrueReco_->SetTitleOffset(2.0);
  h_projYTrueVsDeltaTrueReco_->GetXaxis()->SetTitleOffset(1.20);
  h_projYTrueVsDeltaTrueReco_->GetYaxis()->SetTitleOffset(1.80);
  h_projYTrueVsDeltaTrueReco_->SetDirectory(0);

  dXbin=(Xr2-Xr1)/((float)(Xnb));

  h_GenVsRMS_ = new TH1D();
  h_GenVsRMS_->SetBins(Xnb,Xr1,Xr2);
  h_GenVsRMS_->SetTitleOffset(2.0);
  h_GenVsRMS_->GetXaxis()->SetTitleOffset(1.20);
  h_GenVsRMS_->GetYaxis()->SetTitleOffset(1.30);
  h_GenVsRMS_->SetTitle(titleGenVsRMS_); 
  h_GenVsRMS_->SetStats(0);

  h_GenVsMean_ = new TH1D();
  h_GenVsMean_->SetBins(Xnb,Xr1,Xr2);
  h_GenVsMean_->SetTitleOffset(2.0);
  h_GenVsMean_->GetXaxis()->SetTitleOffset(1.20);
  h_GenVsMean_->GetYaxis()->SetTitleOffset(1.30);
  h_GenVsMean_->SetTitle(titleGenVsMean_);
  h_GenVsMean_->SetStats(0);
    
  for(int i=0;i<Xnb;i++)
    {
      h_GenVsRMS_->SetBinContent(i+1,(h_TrueVsDeltaTrueReco_->ProjectionY("_py",h_TrueVsDeltaTrueReco_->GetXaxis()->FindFixBin(Xr1+i*dXbin) ,h_TrueVsDeltaTrueReco_->GetXaxis()->FindFixBin(Xr1+(i+1)*dXbin),""))->GetRMS());
      h_GenVsRMS_->SetBinError(i+1,(h_TrueVsDeltaTrueReco_->ProjectionY("_py",h_TrueVsDeltaTrueReco_->GetXaxis()->FindFixBin(Xr1+i*dXbin) ,h_TrueVsDeltaTrueReco_->GetXaxis()->FindFixBin(Xr1+(i+1)*dXbin),""))->GetRMSError());

      h_GenVsMean_->SetBinContent(i+1,(h_TrueVsDeltaTrueReco_->ProjectionY("_py",h_TrueVsDeltaTrueReco_->GetXaxis()->FindFixBin(Xr1+i*dXbin) ,h_TrueVsDeltaTrueReco_->GetXaxis()->FindFixBin(Xr1+(i+1)*dXbin),""))->GetMean());
      h_GenVsMean_->SetBinError(i+1,(h_TrueVsDeltaTrueReco_->ProjectionY("_py",h_TrueVsDeltaTrueReco_->GetXaxis()->FindFixBin(Xr1+i*dXbin) ,h_TrueVsDeltaTrueReco_->GetXaxis()->FindFixBin(Xr1+(i+1)*dXbin),""))->GetMeanError()); 
      
    }

  h_GenVsRMS_->SetDirectory(0);
  h_GenVsMean_->SetDirectory(0);

  return;
}

DrawQualityPlots::DrawQualityPlots(const TString& plotsDir, const TString& plotsBaseName, const TString& ch, const TString& year):
plotsDir_(plotsDir),
plotsBaseName_(plotsBaseName),
channel_(ch),
year_(year),
hResp_(0),
hGen_(0)
{
}

void DrawQualityPlots::drawQualityPlots(TH2D *h_TrueVsDeltaTrueReco, TH1D *h_projYTrueVsDeltaTrueReco, TH1D *h_GenVsRMS, TH1D *h_GenVsMean)
{

  TString outDir = utils::assignFolder(plotsDir_, channel_, "Nominal");

  TCanvas *canvas = new TCanvas("canvas","",585, 246, 863, 837);
  canvas->Divide(2,2);
  canvas->GetPad(3)->SetGridy(1);
  canvas->GetPad(4)->SetGridy(1);

  canvas->cd(1);
  h_projYTrueVsDeltaTrueReco->Draw();

  canvas->cd(2);
  h_TrueVsDeltaTrueReco->Draw("colz");

  canvas->cd(3);
  canvas->SetGridy();
  h_GenVsRMS->Draw("ep"); 

  canvas->cd(4);
  canvas->SetGridy();
  h_GenVsMean->Draw("ep");
  
  canvas->Print(outDir + plotsBaseName_ + "_" + year_ + ".pdf");
  canvas->Clear();

  TCanvas *c = new TCanvas("c","",585, 246, 863, 837);
  c->SetGridy(1);

  TFile out_root(outDir + plotsBaseName_ + "_" + year_ + "_source.root", "RECREATE");

  h_projYTrueVsDeltaTrueReco->Draw();
  h_projYTrueVsDeltaTrueReco->Write(plotsBaseName_ + "_" + year_ + "_projy");
  c->Print(outDir + plotsBaseName_ + "_" + year_ + "_projy.pdf");

  h_TrueVsDeltaTrueReco->Draw("colz");
  h_TrueVsDeltaTrueReco->Write(plotsBaseName_ + "_" + year_ + "_TrueVsDeltaTrueReco");
  c->Print(outDir + plotsBaseName_ + "_" + year_ + "_TrueVsDeltaTrueReco.pdf");

  h_GenVsRMS->Draw("ep"); 
  h_GenVsRMS->Write(plotsBaseName_ + "_" + year_ + "_RMS");
  c->Print(outDir + plotsBaseName_ + "_" + year_ + "_RMS.pdf");

  h_GenVsMean->Draw("ep");
  h_GenVsMean->Write(plotsBaseName_ + "_" + year_ + "_mean");
  c->Print(outDir + plotsBaseName_ + "_" + year_ + "_mean.pdf");

  out_root.Close();
  c->Clear();

  return;
}



TGraphAsymmErrors* DrawQualityPlots::drawAsGraph(const TH1D* h, const bool flagUnc/* = true*/, const double offset/* = 0.5*/, const TString& option/* = "ZP0"*/) const
{
  TGraphAsymmErrors* g = new TGraphAsymmErrors;
  for(int b = 0; b < h->GetNbinsX(); b++)
    {
      double x = h->GetBinLowEdge(b + 1) + offset * h->GetBinWidth(b + 1);
      double y = h->GetBinContent(b + 1);
      double uncLowY = h->GetBinError(b + 1);
      double uncHighY = h->GetBinError(b + 1);
      g->SetPoint(b, x, y);
      g->SetPointError(b, 0.0, 0.0, uncLowY, uncHighY);
    }
  g->SetLineColor(h->GetLineColor());
  g->SetMarkerColor(h->GetMarkerColor());
  g->SetMarkerStyle(h->GetMarkerStyle());
  g->SetLineStyle(h->GetLineStyle());
  g->SetMarkerSize(h->GetMarkerSize());
  TString drawOption = option;
  if(!flagUnc)
    drawOption += "X";
  g->Draw(drawOption);
  return g;
}



void DrawQualityPlots::rebin2D(TH1 *h, Int_t ngx, Int_t ngy) const
{
  //Rebin 2-d histogram h, grouping ngx bins together along X
  //and ngy bins together along Y
  //NB: this macro ignores histogram errors if defined
   
  //make a clone of h
  TH1 *hold = (TH1*)h->Clone();
  hold->SetDirectory(0);

  int nbinsx = hold->GetXaxis()->GetNbins();
  int nbinsy = hold->GetYaxis()->GetNbins();
  float xmin  = hold->GetXaxis()->GetXmin();
  float xmax  = hold->GetXaxis()->GetXmax();
  float ymin  = hold->GetYaxis()->GetXmin();
  float ymax  = hold->GetYaxis()->GetXmax();
  int nx = nbinsx/ngx;
  int ny = nbinsy/ngy;
  h->SetBins (nx,xmin,xmax,ny,ymin,ymax);

  //loop on all bins to reset contents and errors
  double cu;
  float bx,by;
  int ix,iy,ibin,bin,binx,biny;
  for (biny=1;biny<=nbinsy;biny++) {
    for (binx=1;binx<=nbinsx;binx++) {
      ibin = h->GetBin(binx,biny);
      h->SetBinContent(ibin,0);
    }
  }
  //loop on all bins and refill
  for (biny=1;biny<=nbinsy;biny++) {
    by  = hold->GetYaxis()->GetBinCenter(biny);
    iy  = h->GetYaxis()->FindBin(by);
    for (binx=1;binx<=nbinsx;binx++) {
      bx = hold->GetXaxis()->GetBinCenter(binx);
      ix  = h->GetXaxis()->FindBin(bx);
      bin = hold->GetBin(binx,biny);
      ibin= h->GetBin(ix,iy);
      cu  = hold->GetBinContent(bin);
      h->AddBinContent(ibin,cu);
    }
  }
  delete hold;          
}



void DrawQualityPlots::calculatePSE(const float& rebin) const
{

  TString outDir = utils::assignFolder(plotsDir_, channel_, "Nominal");

  TH2D *hResp = hResp_;
  TH1D *hGen = hGen_;
  this->rebin2D(hResp_, rebin, rebin);
  hGen->Rebin(rebin);


  double markerSize = 1.5;
  int pixel = 500;
                                                                         
  int n = hResp->GetNbinsX();

  TH1D* hp = new TH1D("", "", n, 0.0, n); 
  TH1D* hs = new TH1D("", "", n, 0.0, n); 
  TH1D* he = new TH1D("", "", n, 0.0, n); 

  for(int b = 1; b <= n; b++)
    {
      double GenInBinAndRecInBin = hResp->GetBinContent(b, b);
      double RecInBin = hResp->Integral(b, b, 1, n);
      double GenInBinAndRec = hResp->Integral(1, n, b, b);
      double GenInBinAll = hGen->GetBinContent(b);
      //std::cout<< "b = " << b <<std::endl;
      //std::cout<< "GenInBinAndRec/GenInBinAll = " << GenInBinAndRec / GenInBinAll <<std::endl;
      hp->SetBinContent(b, GenInBinAndRecInBin / RecInBin);
      hs->SetBinContent(b, GenInBinAndRecInBin / GenInBinAndRec);
      he->SetBinContent(b, GenInBinAndRec / GenInBinAll);
    }

  TCanvas c("", "", pixel, pixel);
  c.cd();
  TLegend* leg = new TLegend(0.15, 0.67, 0.40, 0.85);
  leg->SetTextFont(62);
  TH2D* hr = new TH2D("", "", 1, he->GetBinLowEdge(1), he->GetBinLowEdge(n + 1), 1, 0.0, 1.0);
    hr->GetXaxis()->SetTitle("Bin");
  hr->SetStats(0);
  hr->Draw();
  hp->SetMarkerColor(4);
  hp->SetMarkerStyle(23);
  hp->SetMarkerSize(markerSize);
  leg->AddEntry(hp, "Purity", "p");
  this->drawAsGraph(hp);
  hs->SetMarkerColor(2);
  hs->SetMarkerStyle(22);
  hs->SetMarkerSize(markerSize);
  leg->AddEntry(hs, "Stability", "p");
  this->drawAsGraph(hs);
  he->SetMarkerColor(8);
  he->SetMarkerStyle(20);
  he->SetMarkerSize(markerSize);
  leg->AddEntry(he, "Efficiency", "p");
  this->drawAsGraph(he);
                                                                                                                                                                                   
  leg->Draw();
  c.Print(outDir + "pse.pdf");

  return;
}



void DrawQualityPlots::drawMigrationMatrix(const TString& varXTitle, const float& rebin) const
{

  TString outDir = utils::assignFolder(plotsDir_, channel_, "Nominal");

  TH2D *hResp = hResp_;
  this->rebin2D(hResp_, rebin, rebin);

  int pixel = 500;
 
  int n = hResp->GetNbinsX();

  TH2D* hMig = (TH2D*) hResp->Clone();

  hMig->SetTitle(""); 
  hMig->GetXaxis()->SetTitle(varXTitle + ", rec. level");
  hMig->GetYaxis()->SetTitle(varXTitle + ", gen. level");

  for(int by = 1; by <= n; by++)
    {
      double GenAndRec = hResp->Integral(1, n, by, by);
      for(int bx = 1; bx <= n; bx++)
	for(int bx = 1; bx <= n; bx++)
	  {
	    double val = hResp->GetBinContent(bx, by);
	    hMig->SetBinContent(bx, by, val / GenAndRec);
	  }
    }

  TCanvas c("", "", pixel, pixel);
  c.SetMargin(0.10, 0.15, 0.10, 0.05);
  // decrease font if there are many bins
  double scaleFont = 1.0;
  if(hMig->GetNbinsX() > 40)
    scaleFont = 0.37;
  else if(hMig->GetNbinsX() > 20)
    scaleFont = 0.60;
  hMig->SetMarkerSize(hMig->GetMarkerSize() * scaleFont);
  hMig->SetStats(0);
  hMig->Draw("text colz");
  c.Print(outDir + "mig.pdf");

  return;
}



CompareQualityPlots::CompareQualityPlots(const TString& compareDir, const TString& compareBaseName, const TString& legend1, const TString& legend2, const TString& ch, const TString& year):
compareDir_(compareDir),
compareBaseName_(compareBaseName),
channel_(ch),
year_(year),
legend1_(legend1),
legend2_(legend2)
{
}

void CompareQualityPlots::comparePlots(const TString& filePath1, const TString& histName1, const TString& filePath2, const TString& histName2)
{

  TFile rootFile1(filePath1);
  TH1D *h_GenVsRMS_1 = (TH1D*)rootFile1.Get(histName1+"_RMS");
  TH1D *h_GenVsMean_1 = (TH1D*)rootFile1.Get(histName1+"_mean");

  TFile rootFile2(filePath2);
  TH1D *h_GenVsRMS_2 = (TH1D*)rootFile2.Get(histName2+"_RMS");
  TH1D *h_GenVsMean_2 = (TH1D*)rootFile2.Get(histName2+"_mean");

  TString outDir = utils::assignFolder(compareDir_, channel_, "Nominal");

  TCanvas *c = new TCanvas("c","",585, 246, 863, 837);
  c->SetGridy(1);

  //create the legend
  TLegend* leg = new TLegend(0.15, 0.67, 0.40, 0.85, NULL,"brNDC");
  leg->SetTextFont(62);
  leg->SetBorderSize(1);
  leg->SetLineColor(0);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  
  //create a text box
  TLatex * tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.038);
  tex->SetTextFont(42);

  h_GenVsRMS_1->SetLineColor(kBlue);
  h_GenVsRMS_1->Draw();
  h_GenVsRMS_2->SetLineColor(kRed);
  h_GenVsRMS_2->Draw("same");
  leg->AddEntry(h_GenVsRMS_1, legend1_ ,"l");
  leg->AddEntry(h_GenVsRMS_2, legend2_ ,"l");
  leg->Draw();
  c->Print(outDir + compareBaseName_ + "_" + year_ + "_GenVsRMS.pdf");

  h_GenVsMean_1->SetLineColor(kBlue);
  h_GenVsMean_1->Draw();
  h_GenVsMean_2->SetLineColor(kRed);
  h_GenVsMean_2->Draw("same");
  leg->Draw();
  c->Print(outDir + compareBaseName_ + "_" + year_ + "_GenVsMean.pdf");

  c->Clear();

  delete h_GenVsRMS_1;
  delete h_GenVsRMS_2;
  delete h_GenVsMean_1;
  delete h_GenVsMean_2;


  return;
}


int main(int argc, char** argv) {

  CLParameter<std::string> opt_c("c", "Specify a certain channel (ee, emu, mumu). No channel specified = run on all channels", false, 1, 1, common::makeStringCheck(Channel::convert(Channel::allowedChannelsAnalysis)));
  CLParameter<std::string> opt_year("y", "Specify year, valid: 2016preVFP, 2016postVFP, 2017, 2018 ", true, 1, 1, common::makeStringCheck({"2016preVFP", "2016postVFP", "2017", "2018"}));
  CLParameter<std::string> opt_kr_flag("kr_flag", "Specify kin reco type, valid options: full (0), loose (9), both (-1). By default both are done", true, 1, 1, common::makeStringCheck({"0", "9", "-1"}));

  CLAnalyser::interpretGlobal(argc, argv);

  std::string year; 
  if (opt_year.isSet()) year = opt_year.getArguments()[0];
  std::cout << "Processing setup for year: " << year << std::endl;

  std::vector<std::string> channels {"emu", "ee", "mumu"};
  if (opt_c.isSet()) channels = opt_c.getArguments();
  std::cout << "Processing channels: ";
  for (auto ch: channels) std::cout << ch << " ";
  std::cout << "\n";

  int kr_flag = -1; 
  if (opt_kr_flag.isSet()) kr_flag = stoi(opt_kr_flag.getArguments()[0]);

  TH1::AddDirectory(kFALSE);

  for (auto ch: channels){

    std::string input_dir = "./selectionRoot_" + year + "/Nominal/" + ch + "/";

    std::cout << "input_dir: " << input_dir << std::endl;

    std::string input_file;

    if (year == "2016preVFP")
      input_file =  ch + "_ttbarsignalplustau_fromDilepton_2016ULpreVFP.root";
    else if (year == "2016postVFP")
      input_file =  ch + "_ttbarsignalplustau_fromDilepton_2016ULpostVFP.root";
    else if (year == "2017")
      input_file =  ch + "_ttbarsignalplustau_fromDilepton_2017UL.root";
    else if (year == "2018")
      input_file =  ch + "_ttbarsignalplustau_fromDilepton_2018UL.root";

    std::cout << "input_file: " << input_file << std::endl;

    std::string output_dir = "./QualityPlots_" +year +"/";

    if ((kr_flag == 0 || kr_flag == -1)) {

      //--------------------------------------- Pt Top from full kin reco ---------------------------------------------       

      const char *titleTrueVsDeltaTrueReco_ptTop_full =";p_{t, true}^{top}, GeV;p_{t, true}^{top} - p_{t, reco}^{top}, GeV";
      const char *titleProjYTrueVsDeltaTrueReco_ptTop_full ="p_{t, true}^{top} - p_{t, reco}^{top};#Deltap_{t}^{top}, Gev;Nentries";
      const char *titleGenVsRMS_ptTop_full ="RMS (p_{t, true}^{top} - p_{t, reco}^{top}) vs p_{t, true}^{top}; p_{t, true}^{top}, GeV;RMS";
      const char *titleGenVsMean_ptTop_full ="Mean (p_{t, true}^{top} - p_{t, reco}^{top}) vs p_{t, true}^{top}; p_{t, true}^{top}, GeV;Mean";
      
      PrepareQualityPlots fullKinRecoPtTop(titleTrueVsDeltaTrueReco_ptTop_full, titleProjYTrueVsDeltaTrueReco_ptTop_full, titleGenVsRMS_ptTop_full, titleGenVsMean_ptTop_full);
      
      fullKinRecoPtTop.prepareQualityPlots(input_dir, input_file, "kinRecoQualityStudies_GenVsDeltaTrueReco_ToppT_step8", 10, 30., 500.);
      
      TH2D *h_TrueVsDeltaTrueReco_ptTop_full = fullKinRecoPtTop.getHistoTrueVsDeltaTrueReco();
      h_TrueVsDeltaTrueReco_ptTop_full->SetDirectory(0);
      TH1D *h_projYTrueVsDeltaTrueReco_ptTop_full = fullKinRecoPtTop.getHistoProjYTrueVsDeltaTrueReco();
      h_projYTrueVsDeltaTrueReco_ptTop_full->Rebin(4);
      h_projYTrueVsDeltaTrueReco_ptTop_full->GetXaxis()->SetRangeUser(-500,500);
      h_projYTrueVsDeltaTrueReco_ptTop_full->SetDirectory(0);
      TH1D *h_GenVsRMS_ptTop_full = fullKinRecoPtTop.getHistoGenVsRMS();
      h_GenVsRMS_ptTop_full->GetYaxis()->SetRangeUser(40,130);
      h_GenVsRMS_ptTop_full->SetDirectory(0);
      TH1D *h_GenVsMean_ptTop_full = fullKinRecoPtTop.getHistoGenVsMean();
      h_GenVsMean_ptTop_full->GetYaxis()->SetRangeUser(-30,100);
      h_GenVsMean_ptTop_full->SetDirectory(0);
      
      DrawQualityPlots drawFullKinRecoPtTop(output_dir + "fullKinReco/ptTop", "dPtTopFull", ch, year);
      drawFullKinRecoPtTop.drawQualityPlots(h_TrueVsDeltaTrueReco_ptTop_full, h_projYTrueVsDeltaTrueReco_ptTop_full, h_GenVsRMS_ptTop_full, h_GenVsMean_ptTop_full);
      
      //set response and gen histograms for PSE studies and draw migration matrix
      drawFullKinRecoPtTop.setHistoRsp(input_dir, input_file, "kinRecoQualityStudies_RecoVsGenToppT_step8");
      drawFullKinRecoPtTop.setHistoGen(input_dir, input_file, "kinRecoQualityStudies_GenToppT_step8");
      drawFullKinRecoPtTop.calculatePSE(25);
      drawFullKinRecoPtTop.drawMigrationMatrix("top pt", 25);
      
      delete h_TrueVsDeltaTrueReco_ptTop_full;
      delete h_projYTrueVsDeltaTrueReco_ptTop_full;
      delete h_GenVsRMS_ptTop_full;
      delete h_GenVsMean_ptTop_full;
      
      //--------------------------------------- Eta Top from full kin reco ---------------------------------------------
      
      const char *titleTrueVsDeltaTrueReco_etaTop_full =";y_{true}^{top} ;y_{true}^{top} - y_{reco}^{top}";
      const char *titleProjYTrueVsDeltaTrueReco_etaTop_full ="y_{true}^{top} - y_{reco}^{top};#Deltay^{top} ;Nentries";
      const char *titleGenVsRMS_etaTop_full ="RMS (y_{true}^{top} - y_{reco}^{top}) vs y_{true}^{top}; y_{true}^{top} ;RMS";
      const char *titleGenVsMean_etaTop_full ="Mean (y_{true}^{top} - y_{reco}^{top}) vs y_{true}^{top}; y_{true}^{top} ;Mean";
      
      PrepareQualityPlots fullKinRecoEtaTop(titleTrueVsDeltaTrueReco_etaTop_full, titleProjYTrueVsDeltaTrueReco_etaTop_full, titleGenVsRMS_etaTop_full, titleGenVsMean_etaTop_full);
      
      fullKinRecoEtaTop.prepareQualityPlots(input_dir, input_file, "kinRecoQualityStudies_GenVsDeltaTrueReco_TopRapidity_step8", 11, -2.5, 2.5);
      
      TH2D *h_TrueVsDeltaTrueReco_etaTop_full = fullKinRecoEtaTop.getHistoTrueVsDeltaTrueReco();
      h_TrueVsDeltaTrueReco_etaTop_full->GetXaxis()->SetRangeUser(-3.5,3.5);
      h_TrueVsDeltaTrueReco_etaTop_full->GetYaxis()->SetRangeUser(-3.5,3.5);
      h_TrueVsDeltaTrueReco_etaTop_full->SetDirectory(0);
      TH1D *h_projYTrueVsDeltaTrueReco_etaTop_full = fullKinRecoEtaTop.getHistoProjYTrueVsDeltaTrueReco();
      h_projYTrueVsDeltaTrueReco_etaTop_full->Rebin(2);
      h_projYTrueVsDeltaTrueReco_etaTop_full->SetDirectory(0);
      TH1D *h_GenVsRMS_etaTop_full = fullKinRecoEtaTop.getHistoGenVsRMS();
      h_GenVsRMS_etaTop_full->GetYaxis()->SetRangeUser(0,0.7);
      h_GenVsRMS_etaTop_full->SetDirectory(0);
      TH1D *h_GenVsMean_etaTop_full = fullKinRecoEtaTop.getHistoGenVsMean();
      h_GenVsMean_etaTop_full->GetYaxis()->SetRangeUser(-1,1);
      h_GenVsMean_etaTop_full->SetDirectory(0);
      
      DrawQualityPlots drawFullKinRecoEtaTop(output_dir + "fullKinReco/etaTop", "dEtaTopFull", ch, year);
      drawFullKinRecoEtaTop.drawQualityPlots(h_TrueVsDeltaTrueReco_etaTop_full, h_projYTrueVsDeltaTrueReco_etaTop_full, h_GenVsRMS_etaTop_full, h_GenVsMean_etaTop_full);
      
      //set response and gen histograms for PSE studies and draw migration matrix
      drawFullKinRecoEtaTop.setHistoRsp(input_dir, input_file, "kinRecoQualityStudies_RecoVsGenTopRapidity_step8");
      drawFullKinRecoEtaTop.setHistoGen(input_dir, input_file, "kinRecoQualityStudies_GenTopRapidity_step8");
      drawFullKinRecoEtaTop.calculatePSE(5);
      drawFullKinRecoEtaTop.drawMigrationMatrix("top rapidity", 5);
      
      delete h_TrueVsDeltaTrueReco_etaTop_full;
      delete h_projYTrueVsDeltaTrueReco_etaTop_full;
      delete h_GenVsRMS_etaTop_full;
      delete h_GenVsMean_etaTop_full;
      
      //--------------------------------------- Pt TTbar from full kin reco ---------------------------------------------       
      
      const char *titleTrueVsDeltaTrueReco_ptTTbar_full =";p_{t, true}^{t#bart}, GeV;p_{t, true}^{t#bart} - p_{t, reco}^{t#bart}, GeV";
      const char *titleProjYTrueVsDeltaTrueReco_ptTTbar_full ="p_{t, true}^{t#bart} - p_{t, reco}^{t#bart};#Deltap_{t}^{t#bart}, Gev;Nentries";
      const char *titleGenVsRMS_ptTTbar_full ="RMS (p_{t, true}^{t#bart} - p_{t, reco}^{t#bart}) vs p_{t, true}^{t#bart}; p_{t, true}^{t#bart}, GeV;RMS";
      const char *titleGenVsMean_ptTTbar_full ="Mean (p_{t, true}^{t#bart} - p_{t, reco}^{t#bart}) vs p_{t, true}^{t#bart}; p_{t, true}^{t#bart}, GeV;Mean";
      
      PrepareQualityPlots fullKinRecoPtTTbar(titleTrueVsDeltaTrueReco_ptTTbar_full, titleProjYTrueVsDeltaTrueReco_ptTTbar_full, titleGenVsRMS_ptTTbar_full, titleGenVsMean_ptTTbar_full);
      
      fullKinRecoPtTTbar.prepareQualityPlots(input_dir, input_file, "kinRecoQualityStudies_GenVsDeltaTrueReco_TTBarpT_step8", 20, 60., 1000.);
      
      TH2D *h_TrueVsDeltaTrueReco_ptTTbar_full = fullKinRecoPtTTbar.getHistoTrueVsDeltaTrueReco();
      h_TrueVsDeltaTrueReco_ptTTbar_full->SetDirectory(0);
      TH1D *h_projYTrueVsDeltaTrueReco_ptTTbar_full = fullKinRecoPtTTbar.getHistoProjYTrueVsDeltaTrueReco();
      h_projYTrueVsDeltaTrueReco_ptTTbar_full->Rebin(4);
      h_projYTrueVsDeltaTrueReco_ptTTbar_full->GetXaxis()->SetRangeUser(-500,500);
      h_projYTrueVsDeltaTrueReco_ptTTbar_full->SetDirectory(0);
      TH1D *h_GenVsRMS_ptTTbar_full = fullKinRecoPtTTbar.getHistoGenVsRMS();
      h_GenVsRMS_ptTTbar_full->GetYaxis()->SetRangeUser(40,130);
      h_GenVsRMS_ptTTbar_full->SetDirectory(0);
      TH1D *h_GenVsMean_ptTTbar_full = fullKinRecoPtTTbar.getHistoGenVsMean();
      h_GenVsMean_ptTTbar_full->GetYaxis()->SetRangeUser(-30,100);
      h_GenVsMean_ptTTbar_full->SetDirectory(0);
      
      DrawQualityPlots drawFullKinRecoPtTTbar(output_dir + "fullKinReco/ptTTbar", "dPtTTbarFull", ch, year);
      drawFullKinRecoPtTTbar.drawQualityPlots(h_TrueVsDeltaTrueReco_ptTTbar_full, h_projYTrueVsDeltaTrueReco_ptTTbar_full, h_GenVsRMS_ptTTbar_full, h_GenVsMean_ptTTbar_full);
      
      //set response and gen histograms for PSE studies and draw migration matrix
      drawFullKinRecoPtTTbar.setHistoRsp(input_dir, input_file, "kinRecoQualityStudies_RecoVsGenTTBarpT_step8");
      drawFullKinRecoPtTTbar.setHistoGen(input_dir, input_file, "kinRecoQualityStudies_GenTTBarpT_step8");
      drawFullKinRecoPtTTbar.calculatePSE(10);
      drawFullKinRecoPtTTbar.drawMigrationMatrix("ttbar pt", 10);
      
      delete h_TrueVsDeltaTrueReco_ptTTbar_full;
      delete h_projYTrueVsDeltaTrueReco_ptTTbar_full;
      delete h_GenVsRMS_ptTTbar_full;
      delete h_GenVsMean_ptTTbar_full;
      
      //--------------------------------------- Eta TTbar from full kin reco ---------------------------------------------
      
      const char *titleTrueVsDeltaTrueReco_etaTTbar_full =";y_{true}^{t#bart} ;y_{true}^{t#bart} - y_{reco}^{t#bart}";
      const char *titleProjYTrueVsDeltaTrueReco_etaTTbar_full ="y_{true}^{t#bart} - y_{reco}^{t#bart};#Deltay^{t#bart} ;Nentries";
      const char *titleGenVsRMS_etaTTbar_full ="RMS (y_{true}^{t#bart} - y_{reco}^{t#bart}) vs y_{true}^{t#bart}; y_{true}^{t#bart} ;RMS";
      const char *titleGenVsMean_etaTTbar_full ="Mean (y_{true}^{t#bart} - y_{reco}^{t#bart}) vs y_{true}^{t#bart}; y_{true}^{t#bart} ;Mean";
      
      PrepareQualityPlots fullKinRecoEtaTTbar(titleTrueVsDeltaTrueReco_etaTTbar_full, titleProjYTrueVsDeltaTrueReco_etaTTbar_full, titleGenVsRMS_etaTTbar_full, titleGenVsMean_etaTTbar_full);
      
      fullKinRecoEtaTTbar.prepareQualityPlots(input_dir, input_file, "kinRecoQualityStudies_GenVsDeltaTrueReco_TTBarRapidity_step8", 22, -5., 5.);
      
      TH2D *h_TrueVsDeltaTrueReco_etaTTbar_full = fullKinRecoEtaTTbar.getHistoTrueVsDeltaTrueReco();
      h_TrueVsDeltaTrueReco_etaTTbar_full->GetXaxis()->SetRangeUser(-3.5,3.5);
      h_TrueVsDeltaTrueReco_etaTTbar_full->GetYaxis()->SetRangeUser(-3.5,3.5);
      h_TrueVsDeltaTrueReco_etaTTbar_full->SetDirectory(0);
      TH1D *h_projYTrueVsDeltaTrueReco_etaTTbar_full = fullKinRecoEtaTTbar.getHistoProjYTrueVsDeltaTrueReco();
      h_projYTrueVsDeltaTrueReco_etaTTbar_full->Rebin(2);
      h_projYTrueVsDeltaTrueReco_etaTTbar_full->SetDirectory(0);
      TH1D *h_GenVsRMS_etaTTbar_full = fullKinRecoEtaTTbar.getHistoGenVsRMS();
      h_GenVsRMS_etaTTbar_full->GetYaxis()->SetRangeUser(0,0.7);
      h_GenVsRMS_etaTTbar_full->SetDirectory(0);
      TH1D *h_GenVsMean_etaTTbar_full = fullKinRecoEtaTTbar.getHistoGenVsMean();
      h_GenVsMean_etaTTbar_full->GetYaxis()->SetRangeUser(-1,1);
      h_GenVsMean_etaTTbar_full->SetDirectory(0);
      
      DrawQualityPlots drawFullKinRecoEtaTTbar(output_dir + "fullKinReco/etaTTbar", "dEtaTTbarFull", ch, year);
      drawFullKinRecoEtaTTbar.drawQualityPlots(h_TrueVsDeltaTrueReco_etaTTbar_full, h_projYTrueVsDeltaTrueReco_etaTTbar_full, h_GenVsRMS_etaTTbar_full, h_GenVsMean_etaTTbar_full);
      
      //set response and gen histograms for PSE studies and draw migration matrix
      drawFullKinRecoEtaTTbar.setHistoRsp(input_dir, input_file, "kinRecoQualityStudies_RecoVsGenTTBarRapidity_step8");
      drawFullKinRecoEtaTTbar.setHistoGen(input_dir, input_file, "kinRecoQualityStudies_GenTTBarRapidity_step8");
      drawFullKinRecoEtaTTbar.calculatePSE(5);
      drawFullKinRecoEtaTTbar.drawMigrationMatrix("ttbar rapidity", 5);
      
      delete h_TrueVsDeltaTrueReco_etaTTbar_full;
      delete h_projYTrueVsDeltaTrueReco_etaTTbar_full;
      delete h_GenVsRMS_etaTTbar_full;
      delete h_GenVsMean_etaTTbar_full;
      
      //--------------------------------------- Mass TTbar from full kin reco ---------------------------------------------
      
      const char *titleTrueVsDeltaTrueReco_mTTbar_full =";M_{true}^{t#bart}, GeV;M_{true}^{t#bart} - M_{reco}^{t#bart}, GeV";
      const char *titleProjYTrueVsDeltaTrueReco_mTTbar_full ="M_{true}^{t#bart} - M_{reco}^{t#bart};#DeltaM^{t#bart}, Gev;Nentries";
      const char *titleGenVsRMS_mTTbar_full ="RMS (M_{true}^{t#bart} - M_{reco}^{t#bart}) vs M_{true}^{t#bart}; M_{true}^{t#bart}, GeV; RMS";
      const char *titleGenVsMean_mTTbar_full ="Mean (M_{true}^{t#bart} - M_{reco}^{t#bart}) vs M_{true}^{t#bart}; M_{true}^{t#bart}, GeV;Mean";
      
      PrepareQualityPlots fullKinRecoMTTbar(titleTrueVsDeltaTrueReco_mTTbar_full, titleProjYTrueVsDeltaTrueReco_mTTbar_full, titleGenVsRMS_mTTbar_full, titleGenVsMean_mTTbar_full);
      
      fullKinRecoMTTbar.prepareQualityPlots(input_dir, input_file, "kinRecoQualityStudies_GenVsDeltaTrueReco_TTBarMass_step8", 14, 350., 1500.);
      
      TH2D *h_TrueVsDeltaTrueReco_mTTbar_full = fullKinRecoMTTbar.getHistoTrueVsDeltaTrueReco();
      h_TrueVsDeltaTrueReco_mTTbar_full->GetXaxis()->SetRangeUser(350,1500);
      h_TrueVsDeltaTrueReco_mTTbar_full->GetYaxis()->SetRangeUser(-1200,1200);
      h_TrueVsDeltaTrueReco_mTTbar_full->SetDirectory(0);
      TH1D *h_projYTrueVsDeltaTrueReco_mTTbar_full = fullKinRecoMTTbar.getHistoProjYTrueVsDeltaTrueReco();
      h_projYTrueVsDeltaTrueReco_mTTbar_full->Rebin(8);
      h_projYTrueVsDeltaTrueReco_mTTbar_full->SetDirectory(0);
      TH1D *h_GenVsRMS_mTTbar_full = fullKinRecoMTTbar.getHistoGenVsRMS();
      h_GenVsRMS_mTTbar_full->GetYaxis()->SetRangeUser(0,300);
      h_GenVsRMS_mTTbar_full->SetDirectory(0);
      TH1D *h_GenVsMean_mTTbar_full = fullKinRecoMTTbar.getHistoGenVsMean();
      h_GenVsMean_mTTbar_full->GetYaxis()->SetRangeUser(-40,230);
      h_GenVsMean_mTTbar_full->SetDirectory(0);
      
      DrawQualityPlots drawFullKinRecoMTTbar(output_dir + "fullKinReco/massTTbar", "dMTTbarFull", ch, year);
      drawFullKinRecoMTTbar.drawQualityPlots(h_TrueVsDeltaTrueReco_mTTbar_full, h_projYTrueVsDeltaTrueReco_mTTbar_full, h_GenVsRMS_mTTbar_full, h_GenVsMean_mTTbar_full);
      
      //set response and gen histograms for PSE studies and draw migration matrix
      drawFullKinRecoMTTbar.setHistoRsp(input_dir, input_file, "kinRecoQualityStudies_RecoVsGenTTBarMass_step8");
      drawFullKinRecoMTTbar.setHistoGen(input_dir, input_file, "kinRecoQualityStudies_GenTTBarMass_step8");
      drawFullKinRecoMTTbar.calculatePSE(15);
      drawFullKinRecoMTTbar.drawMigrationMatrix("ttbar mass", 15);
      
      delete h_TrueVsDeltaTrueReco_mTTbar_full;
      delete h_projYTrueVsDeltaTrueReco_mTTbar_full;
      delete h_GenVsRMS_mTTbar_full;
      delete h_GenVsMean_mTTbar_full;
    
    }

    if ((kr_flag == 9 || kr_flag == -1)) {

      //--------------------------------------- Pt TTbar from loose kin reco ---------------------------------------------       
      
      const char *titleTrueVsDeltaTrueReco_ptTTbar_loose =";p_{t, true}^{t#bart}, GeV;p_{t, true}^{t#bart} - p_{t, reco}^{t#bart}, GeV";
      const char *titleProjYTrueVsDeltaTrueReco_ptTTbar_loose ="p_{t, true}^{t#bart} - p_{t, reco}^{t#bart};#Deltap_{t}^{t#bart}, Gev;Nentries";
      const char *titleGenVsRMS_ptTTbar_loose ="RMS (p_{t, true}^{t#bart} - p_{t, reco}^{t#bart}) vs p_{t, true}^{t#bart}; p_{t, true}^{t#bart}, GeV;RMS";
      const char *titleGenVsMean_ptTTbar_loose ="Mean (p_{t, true}^{t#bart} - p_{t, reco}^{t#bart}) vs p_{t, true}^{t#bart}; p_{t, true}^{t#bart}, GeV;Mean";
      
      PrepareQualityPlots looseKinRecoPtTTbar(titleTrueVsDeltaTrueReco_ptTTbar_loose, titleProjYTrueVsDeltaTrueReco_ptTTbar_loose, titleGenVsRMS_ptTTbar_loose, titleGenVsMean_ptTTbar_loose);
      
      looseKinRecoPtTTbar.prepareQualityPlots(input_dir, input_file, "kinRecoQualityStudies_LooseGenVsDeltaTrueReco_TTbarpT_step7L", 20, 60., 1000.);
      
      TH2D *h_TrueVsDeltaTrueReco_ptTTbar_loose = looseKinRecoPtTTbar.getHistoTrueVsDeltaTrueReco();
      h_TrueVsDeltaTrueReco_ptTTbar_loose->SetDirectory(0);
      TH1D *h_projYTrueVsDeltaTrueReco_ptTTbar_loose = looseKinRecoPtTTbar.getHistoProjYTrueVsDeltaTrueReco();
      h_projYTrueVsDeltaTrueReco_ptTTbar_loose->Rebin(4);
      h_projYTrueVsDeltaTrueReco_ptTTbar_loose->GetXaxis()->SetRangeUser(-500,500);
      h_projYTrueVsDeltaTrueReco_ptTTbar_loose->SetDirectory(0);
      TH1D *h_GenVsRMS_ptTTbar_loose = looseKinRecoPtTTbar.getHistoGenVsRMS();
      h_GenVsRMS_ptTTbar_loose->GetYaxis()->SetRangeUser(40,130);
      h_GenVsRMS_ptTTbar_loose->SetDirectory(0);
      TH1D *h_GenVsMean_ptTTbar_loose = looseKinRecoPtTTbar.getHistoGenVsMean();
      h_GenVsMean_ptTTbar_loose->GetYaxis()->SetRangeUser(-30,100);
      h_GenVsMean_ptTTbar_loose->SetDirectory(0);
      
      DrawQualityPlots drawLooseKinRecoPtTTbar(output_dir + "looseKinReco/ptTTbar", "dPtTTbarLoose", ch, year);
      drawLooseKinRecoPtTTbar.drawQualityPlots(h_TrueVsDeltaTrueReco_ptTTbar_loose, h_projYTrueVsDeltaTrueReco_ptTTbar_loose, h_GenVsRMS_ptTTbar_loose, h_GenVsMean_ptTTbar_loose);
      
      //set response and gen histograms for PSE studies and draw migration matrix
      drawLooseKinRecoPtTTbar.setHistoRsp(input_dir, input_file, "kinRecoQualityStudies_RecoVsGenTTBarpT_step7L");
      drawLooseKinRecoPtTTbar.setHistoGen(input_dir, input_file, "kinRecoQualityStudies_GenTTBarpT_step7L");
      drawLooseKinRecoPtTTbar.calculatePSE(10);
      drawLooseKinRecoPtTTbar.drawMigrationMatrix("loose ttbar pt", 10);
      
      delete h_TrueVsDeltaTrueReco_ptTTbar_loose;
      delete h_projYTrueVsDeltaTrueReco_ptTTbar_loose;
      delete h_GenVsRMS_ptTTbar_loose;
      delete h_GenVsMean_ptTTbar_loose;
      
      //--------------------------------------- Eta TTbar from loose kin reco ---------------------------------------------
      
      const char *titleTrueVsDeltaTrueReco_etaTTbar_loose =";y_{true}^{t#bart} ;y_{true}^{t#bart} - y_{reco}^{t#bart}";
      const char *titleProjYTrueVsDeltaTrueReco_etaTTbar_loose ="y_{true}^{t#bart} - y_{reco}^{t#bart};#Deltay^{t#bart} ;Nentries";
      const char *titleGenVsRMS_etaTTbar_loose ="RMS (y_{true}^{t#bart} - y_{reco}^{t#bart}) vs y_{true}^{t#bart}; y_{true}^{t#bart} ;RMS";
      const char *titleGenVsMean_etaTTbar_loose ="Mean (y_{true}^{t#bart} - y_{reco}^{t#bart}) vs y_{true}^{t#bart}; y_{true}^{t#bart} ;Mean";
      
      PrepareQualityPlots looseKinRecoEtaTTbar(titleTrueVsDeltaTrueReco_etaTTbar_loose, titleProjYTrueVsDeltaTrueReco_etaTTbar_loose, titleGenVsRMS_etaTTbar_loose, titleGenVsMean_etaTTbar_loose);
      
      looseKinRecoEtaTTbar.prepareQualityPlots(input_dir, input_file, "kinRecoQualityStudies_LooseGenVsDeltaTrueReco_TTbarRapidity_step7L", 22, -5., 5.);
      
      TH2D *h_TrueVsDeltaTrueReco_etaTTbar_loose = looseKinRecoEtaTTbar.getHistoTrueVsDeltaTrueReco();
      h_TrueVsDeltaTrueReco_etaTTbar_loose->GetXaxis()->SetRangeUser(-3.5,3.5);
      h_TrueVsDeltaTrueReco_etaTTbar_loose->GetYaxis()->SetRangeUser(-3.5,3.5);
      h_TrueVsDeltaTrueReco_etaTTbar_loose->SetDirectory(0);
      TH1D *h_projYTrueVsDeltaTrueReco_etaTTbar_loose = looseKinRecoEtaTTbar.getHistoProjYTrueVsDeltaTrueReco();
      h_projYTrueVsDeltaTrueReco_etaTTbar_loose->Rebin(2);
      h_projYTrueVsDeltaTrueReco_etaTTbar_loose->SetDirectory(0);
      TH1D *h_GenVsRMS_etaTTbar_loose = looseKinRecoEtaTTbar.getHistoGenVsRMS();
      h_GenVsRMS_etaTTbar_loose->GetYaxis()->SetRangeUser(0,0.7);
      h_GenVsRMS_etaTTbar_loose->SetDirectory(0);
      TH1D *h_GenVsMean_etaTTbar_loose = looseKinRecoEtaTTbar.getHistoGenVsMean();
      h_GenVsMean_etaTTbar_loose->GetYaxis()->SetRangeUser(-1,1);
      h_GenVsMean_etaTTbar_loose->SetDirectory(0);
      
      DrawQualityPlots drawLooseKinRecoEtaTTbar(output_dir + "looseKinReco/etaTTbar", "dEtaTTbarLoose", ch, year);
      drawLooseKinRecoEtaTTbar.drawQualityPlots(h_TrueVsDeltaTrueReco_etaTTbar_loose, h_projYTrueVsDeltaTrueReco_etaTTbar_loose, h_GenVsRMS_etaTTbar_loose, h_GenVsMean_etaTTbar_loose);
      
      //set response and gen histograms for PSE studies and draw migration matrix
      drawLooseKinRecoEtaTTbar.setHistoRsp(input_dir, input_file, "kinRecoQualityStudies_RecoVsGenTTBarRapidity_step7L");
      drawLooseKinRecoEtaTTbar.setHistoGen(input_dir, input_file, "kinRecoQualityStudies_GenTTBarRapidity_step7L");
      drawLooseKinRecoEtaTTbar.calculatePSE(5);
      drawLooseKinRecoEtaTTbar.drawMigrationMatrix("loose ttbar rapidity", 5);
      
      delete h_TrueVsDeltaTrueReco_etaTTbar_loose;
      delete h_projYTrueVsDeltaTrueReco_etaTTbar_loose;
      delete h_GenVsRMS_etaTTbar_loose;
      delete h_GenVsMean_etaTTbar_loose;
      
      //--------------------------------------- Mtt TTbar from loose kin reco ---------------------------------------------
      
      const char *titleTrueVsDeltaTrueReco_mTTbar_loose =";M_{true}^{t#bart}, GeV;M_{true}^{t#bart} - M_{reco}^{t#bart}, GeV";
      const char *titleProjYTrueVsDeltaTrueReco_mTTbar_loose ="M_{true}^{t#bart} - M_{reco}^{t#bart};#DeltaM^{t#bart}, Gev;Nentries";
      const char *titleGenVsRMS_mTTbar_loose ="RMS (M_{true}^{t#bart} - M_{reco}^{t#bart}) vs M_{true}^{t#bart}; M_{true}^{t#bart}, GeV; RMS";
      const char *titleGenVsMean_mTTbar_loose ="Mean (M_{true}^{t#bart} - M_{reco}^{t#bart}) vs M_{true}^{t#bart}; M_{true}^{t#bart}, GeV;Mean";
      
      PrepareQualityPlots looseKinRecoMTTbar(titleTrueVsDeltaTrueReco_mTTbar_loose, titleProjYTrueVsDeltaTrueReco_mTTbar_loose, titleGenVsRMS_mTTbar_loose, titleGenVsMean_mTTbar_loose);
      
      looseKinRecoMTTbar.prepareQualityPlots(input_dir, input_file, "kinRecoQualityStudies_LooseGenVsDeltaTrueReco_TTBarMass_step7L", 14, 350., 1500.);
      
      TH2D *h_TrueVsDeltaTrueReco_mTTbar_loose = looseKinRecoMTTbar.getHistoTrueVsDeltaTrueReco();
      h_TrueVsDeltaTrueReco_mTTbar_loose->GetXaxis()->SetRangeUser(350,1500);
      h_TrueVsDeltaTrueReco_mTTbar_loose->GetYaxis()->SetRangeUser(-1200,1200);
      h_TrueVsDeltaTrueReco_mTTbar_loose->SetDirectory(0);
      TH1D *h_projYTrueVsDeltaTrueReco_mTTbar_loose = looseKinRecoMTTbar.getHistoProjYTrueVsDeltaTrueReco();
      h_projYTrueVsDeltaTrueReco_mTTbar_loose->Rebin(8);
      h_projYTrueVsDeltaTrueReco_mTTbar_loose->SetDirectory(0);
      TH1D *h_GenVsRMS_mTTbar_loose = looseKinRecoMTTbar.getHistoGenVsRMS();
      h_GenVsRMS_mTTbar_loose->GetYaxis()->SetRangeUser(0,500);
      h_GenVsRMS_mTTbar_loose->SetDirectory(0);
      TH1D *h_GenVsMean_mTTbar_loose = looseKinRecoMTTbar.getHistoGenVsMean();
      h_GenVsMean_mTTbar_loose->GetYaxis()->SetRangeUser(-40,230);
      h_GenVsMean_mTTbar_loose->SetDirectory(0);
      
      DrawQualityPlots drawLooseKinRecoMTTbar(output_dir + "looseKinReco/massTTbar", "dMTTbarLoose", ch, year);
      drawLooseKinRecoMTTbar.drawQualityPlots(h_TrueVsDeltaTrueReco_mTTbar_loose, h_projYTrueVsDeltaTrueReco_mTTbar_loose, h_GenVsRMS_mTTbar_loose, h_GenVsMean_mTTbar_loose);
      
      //set response and gen histograms for PSE studies and draw migration matrix
      drawLooseKinRecoMTTbar.setHistoRsp(input_dir, input_file, "kinRecoQualityStudies_RecoVsGenTTBarMass_step7L");
      drawLooseKinRecoMTTbar.setHistoGen(input_dir, input_file, "kinRecoQualityStudies_GenTTBarMass_step7L");
      drawLooseKinRecoMTTbar.calculatePSE(15);
      drawLooseKinRecoMTTbar.drawMigrationMatrix("loose ttbar mass", 15);

      //t bins = h_GenVsRMS_mTTbar_loose->GetNbinsX();
      //uble binContent;
      //uble binCenter;
      //r(int i=0;i<bins;i++)
      //{
      //  binContent = h_GenVsRMS_mTTbar_loose->GetBinContent(i+1);
      //  binCenter = h_GenVsRMS_mTTbar_loose->GetBinCenter(i+1);
      //  std::cout << "bin " << i+1 << ", binCenter = " << binCenter << ", binContent = " << binContent << std::endl;
      //}
      
      delete h_TrueVsDeltaTrueReco_mTTbar_loose;
      delete h_projYTrueVsDeltaTrueReco_mTTbar_loose;
      delete h_GenVsRMS_mTTbar_loose;
      delete h_GenVsMean_mTTbar_loose;

    }

    if (kr_flag == -1) {
    
      //-------------------------- make comparisons between full and loose quality plots (only TTbar) ---------------------
      
      CompareQualityPlots comparePtTTbar(output_dir + "compareFullToLoose/ptTTbar", "dPtTTbar", "full kin reco", "loose kin reco", ch, year);
      comparePtTTbar.comparePlots(output_dir + "fullKinReco/ptTTbar/Nominal/" + ch + "/dPtTTbarFull_" + year + "_source.root", "dPtTTbarFull_" + year, output_dir + "looseKinReco/ptTTbar/Nominal/" + ch + "/dPtTTbarLoose_" + year + "_source.root", "dPtTTbarLoose_" + year);
      
      CompareQualityPlots compareEtaTTbar(output_dir + "compareFullToLoose/etaTTbar", "dEtaTTbar", "full kin reco", "loose kin reco", ch, year);
      compareEtaTTbar.comparePlots(output_dir + "fullKinReco/etaTTbar/Nominal/" + ch + "/dEtaTTbarFull_" + year + "_source.root", "dEtaTTbarFull_" + year, output_dir + "looseKinReco/etaTTbar/Nominal/" + ch + "/dEtaTTbarLoose_" + year + "_source.root", "dEtaTTbarLoose_" + year);
      
      CompareQualityPlots compareMttTTbar(output_dir + "compareFullToLoose/massTTbar", "dMTTbar", "full kin reco", "loose kin reco", ch, year);
      compareMttTTbar.comparePlots(output_dir + "fullKinReco/massTTbar/Nominal/" + ch + "/dMTTbarFull_" + year + "_source.root", "dMTTbarFull_" + year, output_dir + "looseKinReco/massTTbar/Nominal/" + ch + "/dMTTbarLoose_" + year + "_source.root", "dMTTbarLoose_" + year);
      
    }
    
    //--------------------------------------- loose kin reco solution plots ---------------------------------------------  
    
    //TString output_folder = utils::assignFolder((output_dir + "looseKinReco/solutions").data(), ch, "Nominal");


    // ************** probe ttbar mass onset ****************
    //TFile rootFile((input_dir + input_file).data());
    //TH1D *LooseKinRecoTTbarM_step7L = (TH1D*)rootFile.Get("kinRecoQualityStudies_LooseKinRecoTTbarM_step7L");
    //
    //int bins_mttbar = LooseKinRecoTTbarM_step7L->GetNbinsX();
    //double binContent_mttbar;
    //double binCenter_mttbar;
    //for(int i=0;i<bins_mttbar;i++)
    //  {
    //	binContent_mttbar = LooseKinRecoTTbarM_step7L->GetBinContent(i+1);
    //	binCenter_mttbar = LooseKinRecoTTbarM_step7L->GetBinCenter(i+1);
    //	std::cout << "bin " << i+1 << ", binCenter_mttbar = " << binCenter_mttbar << ", binContent_mttbar = " << binContent_mttbar << std::endl;
    //  }

    //TH1D *LooseKinRecoTTbarRapidity_step7L = (TH1D*)rootFile.Get("kinRecoQualityStudies_LooseKinRecoTTbarRapidity_step7L");
    //TH1D *LooseKinRecoTTbarPt_step7L = (TH1D*)rootFile.Get("kinRecoQualityStudies_LooseKinRecoTTbarPt_step7L");
    //TH1D *LooseWWmass_step7L = (TH1D*)rootFile.Get("kinRecoQualityStudies_LooseWWmass_step7L");
    //TH1D *LooseMinMlb_step7L = (TH1D*)rootFile.Get("kinRecoQualityStudies_LooseMinMlb_step7L");
    //
    //TCanvas *canvas = new TCanvas("canvas","",585, 246, 863, 837);
    //canvas->Divide(3,3);
    //
    //canvas->cd(1);
    //LooseKinRecoTTbarM_step7L->Draw();
    //
    //canvas->cd(2);
    //LooseKinRecoTTbarRapidity_step7L->Draw();
    //
    //canvas->cd(3);
    //LooseKinRecoTTbarPt_step7L->Draw(); 
    //
    //canvas->cd(4);
    //LooseWWmass_step7L->Draw();
    //
    //canvas->cd(5);
    //LooseMinMlb_step7L->Draw();
    //
    //canvas->Print(output_folder + "looseSolutionOverview.pdf");
    //
    //delete LooseKinRecoTTbarM_step7L;
    //delete LooseKinRecoTTbarRapidity_step7L;
    //delete LooseKinRecoTTbarPt_step7L; 
    //delete LooseWWmass_step7L;
    //delete LooseMinMlb_step7L;

  } //end channel loop

}
