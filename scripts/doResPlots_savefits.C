#include <sstream>

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 3)
{
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}


double CalcLumiWeight(TFile* fMC_, double XSection_, string era_) {

  double lumi = 0;
  
  if (era_ == "2016preVFP")  lumi = 19500;
  else if (era_ == "2016postVFP")  lumi = 16810;
  else if (era_ == "2017")  lumi = 41480;
  else if (era_ == "2018")  lumi = 59830;

  TH1D *h_NrOfEvts = (TH1D*)fMC_->Get("hNrOfEvts");  
  
  double NrOfEvts = h_NrOfEvts->GetBinContent(1);
  //  double NrOfEvts = h_NrOfEvts->GetEntries();

  double lumiWeight = lumi * XSection_ / NrOfEvts;


  std::cout << "lumiw=" << lumiWeight <<", lumi=" << lumi <<", XSection=" << XSection_ <<", NrOfEvts=" << NrOfEvts << "\n";
 

  return lumiWeight;

}


void fitPeak(TH1 *h1, double& mean, double& meanerr, double& sigma, double& sigmaerr)
{

   mean = 0.;
   meanerr = 0.;
   sigma = 0.;
   sigmaerr = 0.;

   if(h1->GetMean() != 0) {

     mean = h1->GetBinCenter( h1->GetMaximumBin() );
     sigma = h1->GetRMS();

     TF1 *f1 = new TF1("f1","gaus", mean - sigma, mean + sigma);

     h1->Fit("f1", "RMQ");
     //mean = f1->GetParameter(1);
     meanerr = f1->GetParError(1);
     sigma = f1->GetParameter(2);
     sigmaerr = f1->GetParError(2);

     //std::cout<<"fitted1: "<<mean<<" "<<meanerr<<" "<<sigma<<" "<<sigmaerr<<std::endl;
     //std::cout<<"ratio1: "<<mean/h1->GetMean()<<" "<<meanerr/h1->GetMeanError()<<" "<<sigma/h1->GetRMS()<<" "<<sigmaerr/h1->GetRMSError()<<std::endl;

     if( sigmaerr == 0. || sigmaerr > 10* h1->GetRMSError() )
     {
       mean = h1->GetMean();
       meanerr = h1->GetMeanError();
       sigma = h1->GetRMS();
       sigmaerr = h1->GetRMSError();
     }
     else {
       TF1 *f2 = new TF1("f2","gaus", mean - sigma, mean + sigma);

       h1->Fit("f2", "RMQ");
       mean = f2->GetParameter(1);
       meanerr = f2->GetParError(1);
       sigma = f2->GetParameter(2);
       sigmaerr = f2->GetParError(2);
     }

     if( sigmaerr == 0. || sigmaerr > 16* h1->GetRMSError() )
     {
       mean = h1->GetMean();
       meanerr = h1->GetMeanError();
       sigma = h1->GetRMS();
       sigmaerr = h1->GetRMSError();
     }
     else if( sigma > 1.1*h1->GetRMS() )
     {
       mean = h1->GetMean();
       meanerr = h1->GetMeanError();
       sigma = h1->GetRMS();
       sigmaerr = h1->GetRMSError();
     }

   }

/*

*/
   //std::cout<<"fitted2: "<<mean<<" "<<meanerr<<" "<<sigma<<" "<<sigmaerr<<std::endl;
   //std::cout<<"ratio2: "<<mean/h1->GetMean()<<" "<<meanerr/h1->GetMeanError()<<" "<<sigma/h1->GetRMS()<<" "<<sigmaerr/h1->GetRMSError()<<std::endl;
}


void drawResPlots(int& i_canvas,TCanvas *c1,TCanvas *c2,TCanvas *c3,TH1 *hresolution_var,TH1 *hresmean_var,TH1 *hresRMS_var,TH1 *hresMeanFit_var,TH1 *hresSigmaFit_var,TH1 *hresMode_var,TH1 *hbinwidth_var)
{
  double mean,meanerr,sigma,sigmaerr;

  c1->cd(i_canvas);
  TPad *p1 = new TPad("p1", "", 0, 0, 1, 1); 
  p1->SetGridx();
  p1->SetGridy();
  p1->SetTickx(1);
  p1->SetTicky(1);
  p1->SetTopMargin(0.1);
  p1->Draw();
  p1->cd();
  hresolution_var->SetLineWidth(1);
  hresolution_var->SetLineColor(kAzure);
  hresolution_var->SetMarkerStyle(20);
  hresolution_var->SetMarkerSize(0.6);
  hresolution_var->SetMarkerColor(kAzure);
  hresolution_var->SetMaximum(hresolution_var->GetMaximum()*1.3);
  hresolution_var->Draw();
  fitPeak(hresolution_var,mean,meanerr,sigma,sigmaerr);
  double SF_resolution_var = sigma/hresolution_var->GetRMS();
  std::cout<<"Resolution SF for "<<hresolution_var->GetName()<<" = "<<SF_resolution_var<<std::endl;
  c2->cd(i_canvas);
  TPad *p2 = new TPad("p2", "", 0, 0, 1, 1); 
  p2->SetGridx();
  p2->SetGridy();
  p2->SetTickx(1);
  p2->SetTicky(1);
  p2->Draw();
  p2->cd();
  hbinwidth_var->SetMaximum(hbinwidth_var->GetMaximum()*1.5);
  if(hbinwidth_var->GetMaximum() < hresRMS_var->GetMaximum() ) hbinwidth_var->SetMaximum( 1.3*hresRMS_var->GetMaximum() );
  if(hbinwidth_var->GetMaximum() < hresSigmaFit_var->GetMaximum() ) hbinwidth_var->SetMaximum( 1.3*hresSigmaFit_var->GetMaximum() );
  hbinwidth_var->SetMinimum(0);
  hbinwidth_var->SetLineStyle(2);
  hbinwidth_var->SetLineWidth(2);
  hbinwidth_var->Draw();
  hresRMS_var->SetLineColor(3);
  hresRMS_var->SetLineWidth(2);
  hresRMS_var->SetMarkerColor(3);
  hresRMS_var->Draw("same");
  TH1D* hresRMS_var_scaled = (TH1D*) hresRMS_var->Clone();
  hresRMS_var_scaled->Scale(SF_resolution_var);
  hresRMS_var_scaled->SetLineColor(99);
  hresRMS_var_scaled->SetMarkerColor(99);
  //hresRMS_var_scaled->Draw("same");
  hresSigmaFit_var->SetLineColor(2);
  hresSigmaFit_var->SetLineWidth(2);
  hresSigmaFit_var->SetMarkerColor(2);
  hresSigmaFit_var->Draw("same");
  c3->cd(i_canvas);
  TPad *p3 = new TPad("p3", "", 0, 0, 1, 1); 
  p3->SetGridx();
  p3->SetGridy();
  p3->SetTickx(1);
  p3->SetTicky(1);
  p3->Draw();
  p3->cd();
  hresmean_var->SetLineColor(3);
  hresmean_var->Draw();
  hresMode_var->SetLineWidth(1);
  //hresMode_var->SetLineStyle(3);
  hresMode_var->SetLineColor(38);
  hresMode_var->Draw("same");
  hresMeanFit_var->SetLineColor(1);
  hresMeanFit_var->Draw("same");
  ++i_canvas;

}


void fillmeanandRMShistos(TH2 *h2, TH1 *hmean, TH1 *hRMS, TH1 *hmeanfit, TH1 *hsigmafit, TH1 *hresMode, TH1 *hbinwidth, TString obsname)
{

  int nbins = h2->GetNbinsY();
  double binwidth = hmean->GetBinWidth(1);
  double binboundary = hmean->GetBinCenter(1) - 0.5*binwidth;

/*
  for (int i = 3; i <= nbins; ++i)
  {
    TH1D* htemp = (TH1D*) h2->ProjectionX("_px",i,i);



    std::cout<<"raw: "<<htemp->GetMean()<<" "<<htemp->GetMeanError()<<" "<<htemp->GetRMS()<<" "<<htemp->GetRMSError()<<std::endl;
  }
*/

  TCanvas *c = new TCanvas("c", "Resolution_bins",1600,1600);
  c->Divide(3,4);

  for (int i = 1; i <= nbins; ++i)
  {
    TH1D* htemp = (TH1D*) h2->ProjectionX(Form("%s_px_bin%i",h2->GetName(),i),i,i);
    c->cd(i);
    TPad *pad = new TPad(("p"+std::to_string(i)).c_str(), "", 0, 0, 1, 1); 
    pad->SetGridx();
    pad->SetGridy();
    pad->SetTickx(1);
    pad->SetTicky(1);
    pad->SetTopMargin(0.175);
    pad->Draw();
    pad->cd();

    gStyle->SetPaintTextFormat("3.3f");
    
    htemp->SetLineWidth(1);
    htemp->SetLineColor(kAzure);
    htemp->SetMarkerStyle(20);
    htemp->SetMarkerSize(0.6);
    htemp->SetMarkerColor(kAzure);

    htemp->GetXaxis()->SetNdivisions(510);
    htemp->GetXaxis()->SetTitle(to_string_with_precision(binboundary) + " < " + obsname + " < " + to_string_with_precision(binboundary + binwidth) + " (reco-gen)");

    binwidth = hmean->GetBinWidth(i);

    binboundary = hmean->GetBinCenter(i) - 0.5*binwidth;
    binboundary = binboundary + binwidth;

    htemp->SetMaximum(htemp->GetMaximum()*1.3);

    htemp->Draw();

    hmean->SetBinContent(i, htemp->GetMean());
    hmean->SetBinError(i, htemp->GetMeanError());
    hRMS->SetBinContent(i, htemp->GetRMS());
    hRMS->SetBinError(i, htemp->GetRMSError());
    hresMode->SetBinContent(i, htemp->GetBinCenter( htemp->GetMaximumBin() ) );
    hresMode->SetBinError(i, 0.);
    hbinwidth->SetBinContent(i, 2.*h2->GetYaxis()->GetBinWidth(i));
    hbinwidth->SetBinError(i, 0.);

    double mean;
    double meanerr;
    double sigma;
    double sigmaerr;

    fitPeak(htemp,mean,meanerr,sigma,sigmaerr);

    hmeanfit->SetBinContent(i, mean);
    hmeanfit->SetBinError(i, meanerr);
    hsigmafit->SetBinContent(i, sigma);
    hsigmafit->SetBinError(i, sigmaerr);

  }

  c->Print( Form("UnfoldingHistosResolutionPlots/%s_Resolution_bins.pdf",h2->GetName()) );
  c->Print( Form("UnfoldingHistosResolutionPlots/%s_Resolution_bins.root",h2->GetName()) );
  c->Print( Form("UnfoldingHistosResolutionPlots/%s_Resolution_bins.C",h2->GetName()) );
  c->Close();

}


TString ObsName(TString VariableName1, bool islatex){

  TString quantity;

  if(islatex) {
    if(VariableName1.Contains("pt")) quantity = "p_{T}";
    else if(VariableName1.Contains("mass")) quantity = "M";
    else if(VariableName1.Contains("rapidity")) quantity =  "\\eta";
    else if(VariableName1.Contains("delta_phi")) quantity = "|\\Delta\\phi_{\\ell\\ell}|";
    else if(VariableName1.Contains("costheta")) quantity = "\\cos\\theta";
    else if(VariableName1.Contains("1n")) quantity = "\\cos\\theta_{1}^{n}";
    else if(VariableName1.Contains("2n")) quantity = "\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("1r")) quantity = "\\cos\\theta_{1}^{r}";
    else if(VariableName1.Contains("2r")) quantity = "\\cos\\theta_{2}^{r}";
    else if(VariableName1.Contains("1k")) quantity = "\\cos\\theta_{1}^{k}";
    else if(VariableName1.Contains("2k")) quantity = "\\cos\\theta_{2}^{k}";
    else if(VariableName1.Contains("1j")) quantity = "\\cos\\theta_{1}^{k*}";
    else if(VariableName1.Contains("2j")) quantity = "\\cos\\theta_{2}^{k*}";
    else if(VariableName1.Contains("1q")) quantity = "\\cos\\theta_{1}^{r*}";
    else if(VariableName1.Contains("2q")) quantity = "\\cos\\theta_{2}^{r*}";
    else if(VariableName1.Contains("Pnn")) quantity = "\\cos\\theta_{1}^{n}+\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("Mnn")) quantity = "\\cos\\theta_{1}^{n}-\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("Prr")) quantity = "\\cos\\theta_{1}^{r}+\\cos\\theta_{2}^{r}";
    else if(VariableName1.Contains("Mrr")) quantity = "\\cos\\theta_{1}^{r}-\\cos\\theta_{2}^{r}";
    else if(VariableName1.Contains("Pkk")) quantity = "\\cos\\theta_{1}^{k}+\\cos\\theta_{2}^{k}";
    else if(VariableName1.Contains("Mkk")) quantity = "\\cos\\theta_{1}^{k}-\\cos\\theta_{2}^{k}";
    else if(VariableName1.Contains("Pjj")) quantity = "\\cos\\theta_{1}^{k*}+\\cos\\theta_{2}^{k*}";
    else if(VariableName1.Contains("Mjj")) quantity = "\\cos\\theta_{1}^{k*}-\\cos\\theta_{2}^{k*}";
    else if(VariableName1.Contains("Pqq")) quantity = "\\cos\\theta_{1}^{r*}+\\cos\\theta_{2}^{r*}";
    else if(VariableName1.Contains("Mqq")) quantity = "\\cos\\theta_{1}^{r*}-\\cos\\theta_{2}^{r*}";
    else if(VariableName1.Contains("c_Prk")) quantity = "\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{k}+\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{r}";
    else if(VariableName1.Contains("c_Mrk")) quantity = "\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{k}-\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{r}";
    else if(VariableName1.Contains("c_Pnr")) quantity = "\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{r}+\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("c_Mnr")) quantity = "\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{r}-\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("c_Pnk")) quantity = "\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{k}+\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("c_Mnk")) quantity = "\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{k}-\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("c_han")) quantity = "+\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{k}-\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{r}-\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("c_sca")) quantity = "-\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{k}+\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{r}-\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("c_tra")) quantity = "-\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{k}-\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{r}+\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("c_kjL")) quantity = "-\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{k*}-\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{r}-\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("c_rqL")) quantity = "-\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{k}-\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{r*}-\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("c_rkP")) quantity = "-\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{k}-\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{r}-\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("c_rkM")) quantity = "-\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{k}+\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{r}-\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("c_nrP")) quantity = "-\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{r}-\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{n}-\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{k}";
    else if(VariableName1.Contains("c_nrM")) quantity = "-\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{r}+\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{n}-\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{k}";
    else if(VariableName1.Contains("c_nkP")) quantity = "-\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{k}-\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{n}-\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{r}";
    else if(VariableName1.Contains("c_nkM")) quantity = "-\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{k}+\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{n}-\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{r}";
    else if(VariableName1.Contains("c_kk")) quantity = "\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{k}";
    else if(VariableName1.Contains("c_kn")) quantity = "\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("c_kr")) quantity = "\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{r}";
    else if(VariableName1.Contains("c_nk")) quantity = "\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{k}";
    else if(VariableName1.Contains("c_nn")) quantity = "\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("c_nr")) quantity = "\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{r}";
    else if(VariableName1.Contains("c_rk")) quantity = "\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{k}";
    else if(VariableName1.Contains("c_rn")) quantity = "\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("c_rr")) quantity = "\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{r}";
    else if(VariableName1.Contains("c_kj")) quantity = "\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{k*}";
    else if(VariableName1.Contains("c_jk")) quantity = "\\cos\\theta_{1}^{k*}\\cos\\theta_{2}^{k}";
    else if(VariableName1.Contains("cHel")) quantity = "\\cos\\varphi";
    else if(VariableName1.Contains("cLab")) quantity = "\\cos\\varphi_{\\mathrm{lab}}";
    else if(VariableName1.Contains("kNorm")) quantity = "\\sin\\theta";
    else if(VariableName1.Contains("rNorm")) quantity = "\\sin\\theta";

    else quantity="";
  }
  else {
    if(VariableName1.Contains("pt")) quantity = "p_{T}";
    else if(VariableName1.Contains("mass")) quantity = "M";
    else if(VariableName1.Contains("rapidity")) quantity =  "#eta";
    else if(VariableName1.Contains("delta_phi")) quantity = "|#Delta#phi_{ll}|";
    else if(VariableName1.Contains("costheta")) quantity = "cos#theta";
    else if(VariableName1.Contains("1n")) quantity = "cos#theta_{1}^{n}";
    else if(VariableName1.Contains("2n")) quantity = "cos#theta_{2}^{n}";
    else if(VariableName1.Contains("1r")) quantity = "cos#theta_{1}^{r}";
    else if(VariableName1.Contains("2r")) quantity = "cos#theta_{2}^{r}";
    else if(VariableName1.Contains("1k")) quantity = "cos#theta_{1}^{k}";
    else if(VariableName1.Contains("2k")) quantity = "cos#theta_{2}^{k}";
    else if(VariableName1.Contains("1j")) quantity = "cos#theta_{1}^{k*}";
    else if(VariableName1.Contains("2j")) quantity = "cos#theta_{2}^{k*}";
    else if(VariableName1.Contains("1q")) quantity = "cos#theta_{1}^{r*}";
    else if(VariableName1.Contains("2q")) quantity = "cos#theta_{2}^{r*}";
    else if(VariableName1.Contains("Pnn")) quantity = "cos#theta_{1}^{n}+cos#theta_{2}^{n}";
    else if(VariableName1.Contains("Mnn")) quantity = "cos#theta_{1}^{n}-cos#theta_{2}^{n}";
    else if(VariableName1.Contains("Prr")) quantity = "cos#theta_{1}^{r}+cos#theta_{2}^{r}";
    else if(VariableName1.Contains("Mrr")) quantity = "cos#theta_{1}^{r}-cos#theta_{2}^{r}";
    else if(VariableName1.Contains("Pkk")) quantity = "cos#theta_{1}^{k}+cos#theta_{2}^{k}";
    else if(VariableName1.Contains("Mkk")) quantity = "cos#theta_{1}^{k}-cos#theta_{2}^{k}";
    else if(VariableName1.Contains("Pjj")) quantity = "cos#theta_{1}^{k*}+cos#theta_{2}^{k*}";
    else if(VariableName1.Contains("Mjj")) quantity = "cos#theta_{1}^{k*}-cos#theta_{2}^{k*}";
    else if(VariableName1.Contains("Pqq")) quantity = "cos#theta_{1}^{r*}+cos#theta_{2}^{r*}";
    else if(VariableName1.Contains("Mqq")) quantity = "cos#theta_{1}^{r*}-cos#theta_{2}^{r*}";
    else if(VariableName1.Contains("c_Prk")) quantity = "cos#theta_{1}^{r}cos#theta_{2}^{k}+cos#theta_{1}^{k}cos#theta_{2}^{r}";
    else if(VariableName1.Contains("c_Mrk")) quantity = "cos#theta_{1}^{r}cos#theta_{2}^{k}-cos#theta_{1}^{k}cos#theta_{2}^{r}";
    else if(VariableName1.Contains("c_Pnr")) quantity = "cos#theta_{1}^{n}cos#theta_{2}^{r}+cos#theta_{1}^{r}cos#theta_{2}^{n}";
    else if(VariableName1.Contains("c_Mnr")) quantity = "cos#theta_{1}^{n}cos#theta_{2}^{r}-cos#theta_{1}^{r}cos#theta_{2}^{n}";
    else if(VariableName1.Contains("c_Pnk")) quantity = "cos#theta_{1}^{n}cos#theta_{2}^{k}+cos#theta_{1}^{k}cos#theta_{2}^{n}";
    else if(VariableName1.Contains("c_Mnk")) quantity = "cos#theta_{1}^{n}cos#theta_{2}^{k}-cos#theta_{1}^{k}cos#theta_{2}^{n}";
    else if(VariableName1.Contains("c_han")) quantity = "+cos#theta_{1}^{k}cos#theta_{2}^{k}-cos#theta_{1}^{r}cos#theta_{2}^{r}-cos#theta_{1}^{n}cos#theta_{2}^{n}";
    else if(VariableName1.Contains("c_sca")) quantity = "-cos#theta_{1}^{k}cos#theta_{2}^{k}+cos#theta_{1}^{r}cos#theta_{2}^{r}-cos#theta_{1}^{n}cos#theta_{2}^{n}";
    else if(VariableName1.Contains("c_tra")) quantity = "-cos#theta_{1}^{k}cos#theta_{2}^{k}-cos#theta_{1}^{r}cos#theta_{2}^{r}+cos#theta_{1}^{n}cos#theta_{2}^{n}";
    else if(VariableName1.Contains("c_kjL")) quantity = "-cos#theta_{1}^{k}cos#theta_{2}^{k*}-cos#theta_{1}^{r}cos#theta_{2}^{r}-cos#theta_{1}^{n}cos#theta_{2}^{n}";
    else if(VariableName1.Contains("c_rqL")) quantity = "-cos#theta_{1}^{k}cos#theta_{2}^{k}-cos#theta_{1}^{r}cos#theta_{2}^{r*}-cos#theta_{1}^{n}cos#theta_{2}^{n}";
    else if(VariableName1.Contains("c_rkP")) quantity = "-cos#theta_{1}^{r}cos#theta_{2}^{k}-cos#theta_{1}^{k}cos#theta_{2}^{r}-cos#theta_{1}^{n}cos#theta_{2}^{n}";
    else if(VariableName1.Contains("c_rkM")) quantity = "-cos#theta_{1}^{r}cos#theta_{2}^{k}+cos#theta_{1}^{k}cos#theta_{2}^{r}-cos#theta_{1}^{n}cos#theta_{2}^{n}";
    else if(VariableName1.Contains("c_nrP")) quantity = "-cos#theta_{1}^{n}cos#theta_{2}^{r}-cos#theta_{1}^{r}cos#theta_{2}^{n}-cos#theta_{1}^{k}cos#theta_{2}^{k}";
    else if(VariableName1.Contains("c_nrM")) quantity = "-cos#theta_{1}^{n}cos#theta_{2}^{r}+cos#theta_{1}^{r}cos#theta_{2}^{n}-cos#theta_{1}^{k}cos#theta_{2}^{k}";
    else if(VariableName1.Contains("c_nkP")) quantity = "-cos#theta_{1}^{n}cos#theta_{2}^{k}-cos#theta_{1}^{k}cos#theta_{2}^{n}-cos#theta_{1}^{r}cos#theta_{2}^{r}";
    else if(VariableName1.Contains("c_nkM")) quantity = "-cos#theta_{1}^{n}cos#theta_{2}^{k}+cos#theta_{1}^{k}cos#theta_{2}^{n}-cos#theta_{1}^{r}cos#theta_{2}^{r}";
    else if(VariableName1.Contains("c_kk")) quantity = "cos#theta_{1}^{k}cos#theta_{2}^{k}";
    else if(VariableName1.Contains("c_kn")) quantity = "cos#theta_{1}^{k}cos#theta_{2}^{n}";
    else if(VariableName1.Contains("c_kr")) quantity = "cos#theta_{1}^{k}cos#theta_{2}^{r}";
    else if(VariableName1.Contains("c_nk")) quantity = "cos#theta_{1}^{n}cos#theta_{2}^{k}";
    else if(VariableName1.Contains("c_nn")) quantity = "cos#theta_{1}^{n}cos#theta_{2}^{n}";
    else if(VariableName1.Contains("c_nr")) quantity = "cos#theta_{1}^{n}cos#theta_{2}^{r}";
    else if(VariableName1.Contains("c_rk")) quantity = "cos#theta_{1}^{r}cos#theta_{2}^{k}";
    else if(VariableName1.Contains("c_rn")) quantity = "cos#theta_{1}^{r}cos#theta_{2}^{n}";
    else if(VariableName1.Contains("c_rr")) quantity = "cos#theta_{1}^{r}cos#theta_{2}^{r}";
    else if(VariableName1.Contains("c_kj")) quantity = "cos#theta_{1}^{k}cos#theta_{2}^{k*}";
    else if(VariableName1.Contains("c_jk")) quantity = "cos#theta_{1}^{k*}cos#theta_{2}^{k}";
    else if(VariableName1.Contains("cHel")) quantity = "cos#varphi";
    else if(VariableName1.Contains("cLab")) quantity = "cos#varphi_{lab}";
    else if(VariableName1.Contains("kNorm")) quantity = "sin#theta";
    else if(VariableName1.Contains("rNorm")) quantity = "sin#theta";

    else quantity="";
  }

  return quantity;

}


void doResPlots_savefits()
{

  // To run in root:
  //root -b -q macros/doResPlots_savefits.C

  vector<TFile*> fMC;
  vector<double> fMC_LumiWeight;

  double topxsec = 830.91;

  // prompt Unfolding Histos
  TFile* fMC_ee_prompt_2016preVFP = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2016preVFP/Nominal/ee/histosTUnfold_ee_ttbarsignalplustau_fromDilepton_2016ULpreVFP.root","READ");
  fMC.push_back(fMC_ee_prompt_2016preVFP);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_ee_prompt_2016preVFP, topxsec * 0.10706 * 0.964261576, "2016preVFP"));

  TFile* fMC_ee_prompt_2016postVFP = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2016postVFP/Nominal/ee/histosTUnfold_ee_ttbarsignalplustau_fromDilepton_2016ULpostVFP.root","READ");
  fMC.push_back(fMC_ee_prompt_2016postVFP);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_ee_prompt_2016postVFP, topxsec * 0.10706 * 0.964261576, "2016postVFP"));

  TFile* fMC_ee_prompt_2017 = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2017/Nominal/ee/histosTUnfold_ee_ttbarsignalplustau_fromDilepton_2017UL.root","READ");
  fMC.push_back(fMC_ee_prompt_2017);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_ee_prompt_2017, topxsec * 0.10706 * 0.964261576, "2017"));

  TFile* fMC_ee_prompt_2018 = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2018/Nominal/ee/histosTUnfold_ee_ttbarsignalplustau_fromDilepton_2018UL.root","READ");
  fMC.push_back(fMC_ee_prompt_2018);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_ee_prompt_2018, topxsec * 0.10706 * 0.964261576, "2018"));

  TFile* fMC_emu_prompt_2016preVFP = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2016preVFP/Nominal/emu/histosTUnfold_emu_ttbarsignalplustau_fromDilepton_2016ULpreVFP.root","READ");
  fMC.push_back(fMC_emu_prompt_2016preVFP);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_emu_prompt_2016preVFP, topxsec * 0.10706 * 0.957058875, "2016preVFP"));

  TFile* fMC_emu_prompt_2016postVFP = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2016postVFP/Nominal/emu/histosTUnfold_emu_ttbarsignalplustau_fromDilepton_2016ULpostVFP.root","READ");
  fMC.push_back(fMC_emu_prompt_2016postVFP);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_emu_prompt_2016postVFP, topxsec * 0.10706 * 0.957058875, "2016postVFP"));

  TFile* fMC_emu_prompt_2017 = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2017/Nominal/emu/histosTUnfold_emu_ttbarsignalplustau_fromDilepton_2017UL.root","READ");
  fMC.push_back(fMC_emu_prompt_2017);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_emu_prompt_2017, topxsec * 0.10706 * 0.957058875, "2017"));

  TFile* fMC_emu_prompt_2018 = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2018/Nominal/emu/histosTUnfold_emu_ttbarsignalplustau_fromDilepton_2018UL.root","READ");
  fMC.push_back(fMC_emu_prompt_2018);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_emu_prompt_2018, topxsec * 0.10706 * 0.957058875, "2018"));

  TFile* fMC_mumu_prompt_2016preVFP = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2016preVFP/Nominal/mumu/histosTUnfold_mumu_ttbarsignalplustau_fromDilepton_2016ULpreVFP.root","READ");
  fMC.push_back(fMC_mumu_prompt_2016preVFP);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_mumu_prompt_2016preVFP, topxsec * 0.10706 * 0.949909976, "2016preVFP"));

  TFile* fMC_mumu_prompt_2016postVFP = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2016postVFP/Nominal/mumu/histosTUnfold_mumu_ttbarsignalplustau_fromDilepton_2016ULpostVFP.root","READ");
  fMC.push_back(fMC_mumu_prompt_2016postVFP);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_mumu_prompt_2016postVFP, topxsec * 0.10706 * 0.949909976, "2016postVFP"));

  TFile* fMC_mumu_prompt_2017 = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2017/Nominal/mumu/histosTUnfold_mumu_ttbarsignalplustau_fromDilepton_2017UL.root","READ");
  fMC.push_back(fMC_mumu_prompt_2017);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_mumu_prompt_2017, topxsec * 0.10706 * 0.949909976, "2017"));

  TFile* fMC_mumu_prompt_2018 = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2018/Nominal/mumu/histosTUnfold_mumu_ttbarsignalplustau_fromDilepton_2018UL.root","READ");
  fMC.push_back(fMC_mumu_prompt_2018);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_mumu_prompt_2018, topxsec * 0.10706 * 0.949909976, "2018"));
  //

  // viatau Unfolding Histos
  TFile* fMC_ee_viatau_2016preVFP = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2016preVFP/Nominal/ee/histosTUnfold_ee_ttbarsignalviatau_fromDilepton_2016ULpreVFP.root","READ");
  fMC.push_back(fMC_ee_viatau_2016preVFP);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_ee_viatau_2016preVFP, topxsec * 0.10706 * 1.029827957, "2016preVFP"));

  TFile* fMC_ee_viatau_2016postVFP = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2016postVFP/Nominal/ee/histosTUnfold_ee_ttbarsignalviatau_fromDilepton_2016ULpostVFP.root","READ");
  fMC.push_back(fMC_ee_viatau_2016postVFP);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_ee_viatau_2016postVFP, topxsec * 0.10706 * 1.029827957, "2016postVFP"));

  TFile* fMC_ee_viatau_2017 = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2017/Nominal/ee/histosTUnfold_ee_ttbarsignalviatau_fromDilepton_2017UL.root","READ");
  fMC.push_back(fMC_ee_viatau_2017);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_ee_viatau_2017, topxsec * 0.10706 * 1.029827957, "2017"));

  TFile* fMC_ee_viatau_2018 = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2018/Nominal/ee/histosTUnfold_ee_ttbarsignalviatau_fromDilepton_2018UL.root","READ");
  fMC.push_back(fMC_ee_viatau_2018);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_ee_viatau_2018, topxsec * 0.10706 * 1.029827957, "2018"));

  TFile* fMC_emu_viatau_2016preVFP = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2016preVFP/Nominal/emu/histosTUnfold_emu_ttbarsignalviatau_fromDilepton_2016ULpreVFP.root","READ");
  fMC.push_back(fMC_emu_viatau_2016preVFP);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_emu_viatau_2016preVFP, topxsec * 0.10706 * 1.026209047, "2016preVFP"));

  TFile* fMC_emu_viatau_2016postVFP = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2016postVFP/Nominal/emu/histosTUnfold_emu_ttbarsignalviatau_fromDilepton_2016ULpostVFP.root","READ");
  fMC.push_back(fMC_emu_viatau_2016postVFP);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_emu_viatau_2016postVFP, topxsec * 0.10706 * 1.026209047, "2016postVFP"));

  TFile* fMC_emu_viatau_2017 = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2017/Nominal/emu/histosTUnfold_emu_ttbarsignalviatau_fromDilepton_2017UL.root","READ");
  fMC.push_back(fMC_emu_viatau_2017);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_emu_viatau_2017, topxsec * 0.10706 * 1.026209047, "2017"));

  TFile* fMC_emu_viatau_2018 = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2018/Nominal/emu/histosTUnfold_emu_ttbarsignalviatau_fromDilepton_2018UL.root","READ");
  fMC.push_back(fMC_emu_viatau_2018);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_emu_viatau_2018, topxsec * 0.10706 * 1.026209047, "2018"));

  TFile* fMC_mumu_viatau_2016preVFP = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2016preVFP/Nominal/mumu/histosTUnfold_mumu_ttbarsignalviatau_fromDilepton_2016ULpreVFP.root","READ");
  fMC.push_back(fMC_mumu_viatau_2016preVFP);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_mumu_viatau_2016preVFP, topxsec * 0.10706 * 1.022670477, "2016preVFP"));

  TFile* fMC_mumu_viatau_2016postVFP = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2016postVFP/Nominal/mumu/histosTUnfold_mumu_ttbarsignalviatau_fromDilepton_2016ULpostVFP.root","READ");
  fMC.push_back(fMC_mumu_viatau_2016postVFP);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_mumu_viatau_2016postVFP, topxsec * 0.10706 * 1.022670477, "2016postVFP"));

  TFile* fMC_mumu_viatau_2017 = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2017/Nominal/mumu/histosTUnfold_mumu_ttbarsignalviatau_fromDilepton_2017UL.root","READ");
  fMC.push_back(fMC_mumu_viatau_2017);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_mumu_viatau_2017, topxsec * 0.10706 * 1.022670477, "2017"));

  TFile* fMC_mumu_viatau_2018 = new TFile("/depot/cms/top/bakshi3/TopSpinCorr_Run2_alt_version/CMSSW_10_6_29/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos_2018/Nominal/mumu/histosTUnfold_mumu_ttbarsignalviatau_fromDilepton_2018UL.root","READ");
  fMC.push_back(fMC_mumu_viatau_2018);  
  fMC_LumiWeight.push_back(CalcLumiWeight(fMC_mumu_viatau_2018, topxsec * 0.10706 * 1.022670477, "2018"));
  //


  // Binnings
  const int gen_ntop_ptbin = 7-1;
  const int reco_ntop_ptbin = 14-2;
  double gen_top_ptbins[gen_ntop_ptbin+1] = {0.0, 65.0, 125.0, 200.0, 290.0, 400.0, 1200.0};
  double reco_top_ptbins[reco_ntop_ptbin+1] = {0.0, 32.0, 65.0, 95.0, 125.0, 160.0, 200.0, 245.0, 290.0, 340.0, 400.0, 600.0, 1200.0};

  const int gen_nl_ptbin = 7-1;
  const int reco_nl_ptbin = 14-2;
  double gen_l_ptbins[gen_nl_ptbin+1] = {20.0, 40.0, 70.0, 120.0, 180.0, 400.0, 1200.0};
  double reco_l_ptbins[reco_nl_ptbin+1] = {20.0, 30.0, 40.0, 55.0, 70.0, 95.0, 120.0, 150.0, 180.0, 280.0, 400.0, 600.0, 1200.0};

  const int gen_nlbar_ptbin = 7-1;
  const int reco_nlbar_ptbin = 14-2;
  double gen_lbar_ptbins[gen_nlbar_ptbin+1] = {20.0, 40.0, 70.0, 120.0, 180.0, 400.0, 1200.0};
  double reco_lbar_ptbins[reco_nlbar_ptbin+1] = {20.0, 30.0, 40.0, 55.0, 70.0, 95.0, 120.0, 150.0, 180.0, 280.0, 400.0, 600.0, 1200.0};


  const int gen_nttbar_ptbin = 6-1;
  const int reco_nttbar_ptbin = 12-2;
  double gen_ttbar_ptbins[gen_nttbar_ptbin+1] = {0.0, 30.0, 80.0, 170.0, 300.0, 1200.0};
  double reco_ttbar_ptbins[reco_nttbar_ptbin+1] = {0.0, 15.0, 30.0, 55.0, 80.0, 120.0, 170.0, 225.0, 300.0, 500.0, 1200.0};

  const int gen_nttbar_massbin = 5-1;
  const int reco_nttbar_massbin = 10-2;
  double gen_ttbar_massbins[gen_nttbar_massbin+1] = {300.0, 450.0, 600.0, 800.0, 2000.0};
  double reco_ttbar_massbins[reco_nttbar_massbin+1] = {300.0, 375.0, 450.0, 525.0, 600.0, 700.0, 800.0, 1400.0, 2000.0};

  const int gen_ntop_scatteringangle_ttbarframebin = 5-1;
  const int reco_ntop_scatteringangle_ttbarframebin = 10-2;
  double gen_top_scatteringangle_ttbarframebins[gen_nttbar_massbin+1] = {-1.0, -0.5, 0.0, 0.5, 1.0};
  double reco_top_scatteringangle_ttbarframebins[reco_nttbar_massbin+1] = {-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0};

  const int gen_nllbar_ptbin = 9-1;
  const int reco_nllbar_ptbin = 18-2;
  double gen_llbar_ptbins[gen_nllbar_ptbin+1] = {0.0, 10.0, 20.0, 40.0, 60.0, 100.0, 150.0, 400.0, 457.0};
  double reco_llbar_ptbins[reco_nllbar_ptbin+1] = {0.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 125.0, 150.0, 275.0, 400.0, 455.0, 457.0};

  const int gen_nllbar_massbin = 10-1;
  const int reco_nllbar_massbin = 20-2;
  double gen_llbar_massbins[gen_nllbar_massbin+1] = {20.0, 30.0, 50.0, 76.0, 106.0, 130.0, 170.0, 260.0, 400.0, 447.5};
  double reco_llbar_massbins[reco_nllbar_massbin+1] = {20.0, 25.0, 30.0, 40.0, 50.0, 63.0, 76.0, 90.0, 106.0, 118.0, 130.0, 150.0, 170.0, 215.0, 260.0, 320.0, 400.0, 445.0, 447.5};

  const int gen_ndphibin = 9-1;
  const int reco_ndphibin = 18-2;
  double gen_dphibins[gen_ndphibin+1] = {0.0, 1.0, 1.6, 2.0, 2.4, 2.6, 2.8, 3.0, 3.2};
  double reco_dphibins[reco_ndphibin+1] = {0.0, 0.5, 1.0, 1.3, 1.6, 1.8, 2.0, 2.2, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2};

  const int gen_ndetabin = 9-1;
  const int reco_ndetabin = 18-2;
  double gen_detabins[gen_ndetabin+1] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0, 5.0, 10.0};
  double reco_detabins[reco_ndetabin+1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0};

  const int gen_netabin = 11-1;
  const int reco_netabin = 22-2;
  double gen_etabins[gen_netabin+1] = {-3.2, -2.0, -1.4, -0.8, -0.4, 0.0, 0.4, 0.8, 1.4, 2.0, 3.2};
  double reco_etabins[reco_netabin+1] = {-3.2, -2.6, -2.0, -1.7, -1.4, -1.1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.1, 1.4, 1.7, 2.0, 2.6, 3.2};



  const int gen_ncosbin = 6;
  const int reco_ncosbin = 12;
  double gen_cosbins[gen_ncosbin+1] = {-1.0, -2./3., -1./3., 0., 1./3., 2./3., 1.0};
  double reco_cosbins[reco_ncosbin+1] = {-1.0, -5./6., -2./3., -1./2., -1./3., -1./6., 0., 1./6., 1./3., 1./2., 2./3., 5./6., 1.0};

  //const int gen_nbbin = 11-1;
  //const int reco_nbbin = 22-2;
  //double gen_bbins[gen_nbbin+1] = {-2.0, -1.4, -1.0, -0.6, -0.3, 0.0, 0.3, 0.6, 1.0, 1.4, 2.0};
  //double reco_bbins[reco_nbbin+1] = {-2.0, -1.7, -1.4, -1.2, -1.0, -0.8, -0.6, -0.45, -0.3, -0.15, 0.0, 0.15, 0.3, 0.45, 0.6, 0.8,  1.0, 1.2,  1.4, 1.7,  2.0};

  //const int gen_nbbin = 8;
  //const int reco_nbbin = 16;
  //double gen_bbins[gen_nbbin+1] = {-2.0, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 2.0};
  //double reco_bbins[reco_nbbin+1] = {-2.0, -1.6, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0};

  const int gen_nbbin = 6;
  const int reco_nbbin = 12;
  double gen_bbins[gen_nbbin+1] = {-1.0, -2./3., -1./3., 0., 1./3., 2./3., 1.0};
  double reco_bbins[reco_nbbin+1] = {-1.0, -5./6., -2./3., -1./2., -1./3., -1./6., 0., 1./6., 1./3., 1./2., 2./3., 5./6., 1.0};

  //const int gen_nbbin = 8;
  //const int reco_nbbin = 16;
  //double gen_bbins[gen_nbbin+1] = {-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0};
  //double reco_bbins[reco_nbbin+1] = {-2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0};

  //const int gen_ncbin = 9-1;
  //const int reco_ncbin = 18-2;
  //double gen_cbins[gen_ncbin+1] = {-1.0, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 1.0};
  //double reco_cbins[reco_ncbin+1] = {-1.0, -0.8, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0};

  //const int gen_ncbin = 6;
  //const int reco_ncbin = 12;
  //double gen_cbins[gen_ncbin+1] = {-1.0, -0.6, -0.3, 0.0, 0.3, 0.6, 1.0};
  //double reco_cbins[reco_ncbin+1] = {-1.0, -0.8, -0.6, -0.45, -0.3, -0.15, 0.0, 0.15, 0.3, 0.45, 0.6, 0.8, 1.0};

  const int gen_ncbin = 6;
  const int reco_ncbin = 12;
  double gen_cbins[gen_ncbin+1] = {-1.0, -2./3., -1./3., 0., 1./3., 2./3., 1.0};
  double reco_cbins[reco_ncbin+1] = {-1.0, -5./6., -2./3., -1./2., -1./3., -1./6., 0., 1./6., 1./3., 1./2., 2./3., 5./6., 1.0};

  //const int gen_ncijbin = 6;
  //const int reco_ncijbin = 12;
  //double gen_cijbins[gen_ncijbin+1] = {-1.0, -0.6, -0.3, 0.0, 0.3, 0.6, 1.0};
  //double reco_cijbins[reco_ncijbin+1] = {-1.0, -0.8, -0.6, -0.45, -0.3, -0.15, 0.0, 0.15, 0.3, 0.45, 0.6, 0.8, 1.0};



  const int gen_nllbar_delta_phibin = 6;
  const int reco_nllbar_delta_phibin = 12;
  double gen_llbar_delta_phibins[gen_nllbar_delta_phibin+1] = {0., 1.*TMath::Pi()/6., 2.*TMath::Pi()/6., 3.*TMath::Pi()/6., 4.*TMath::Pi()/6., 5.*TMath::Pi()/6., TMath::Pi()};
  double reco_llbar_delta_phibins[reco_nllbar_delta_phibin+1] = {0., 1.*TMath::Pi()/12., 2.*TMath::Pi()/12., 3.*TMath::Pi()/12., 4.*TMath::Pi()/12., 5.*TMath::Pi()/12., 6.*TMath::Pi()/12., 7.*TMath::Pi()/12., 8.*TMath::Pi()/12., 9.*TMath::Pi()/12., 10.*TMath::Pi()/12., 11.*TMath::Pi()/12., TMath::Pi()};

  //


  //
  string variables[33] = {"b1n", "b2n", "b1r", "b2r", "b1k", "b2k", "b1j", "b2j", "b1q", "b2q", "c_kk", "c_rr", "c_nn", "c_Prk", "c_Mrk", "c_Pnr", "c_Mnr", "c_Pnk", "c_Mnk", "c_han", "c_sca", "c_tra", "c_kjL", "c_rqL", "c_rkP", "c_rkM", "c_nrP", "c_nrM", "c_nkP", "c_nkM", "ll_cHel", "ttbar_mass", "top_scatteringangle_ttbarframe"};

  int reco_nbin[33] = {reco_nbbin, reco_nbbin, reco_nbbin, reco_nbbin, reco_nbbin, reco_nbbin, reco_nbbin, reco_nbbin, reco_nbbin, reco_nbbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncosbin, reco_nttbar_massbin, reco_ntop_scatteringangle_ttbarframebin};

  double* reco_bins[33] = {reco_bbins, reco_bbins, reco_bbins, reco_bbins, reco_bbins, reco_bbins, reco_bbins, reco_bbins, reco_bbins, reco_bbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins, reco_cosbins, reco_ttbar_massbins, reco_top_scatteringangle_ttbarframebins};
  //
  

  // Adding together eras and channels for each variable
  system("mkdir -p UnfoldingHistosResolutionPlots");

  TFile *output_hists = new TFile("UnfoldingHistosResolutionPlots/resolutionhists_fullRun2UL.root","RECREATE");

  for (UInt_t i = 0; i < sizeof(variables)/sizeof(variables[0]); i++) {

    std::cout <<  ("hresolutionbins_"+variables[i]).c_str() << std::endl;

    TH2D *h_var = (TH2D*)fMC[0]->Get(( "hresolutionbins_"+variables[i]).c_str());
    h_var->Scale(fMC_LumiWeight[0]);

    TH2D *hs_var = (TH2D*)h_var->Clone(( "hresolutionbins_"+variables[i]).c_str());

    for (UInt_t j = 1; j < fMC.size(); j++) {

      h_var = (TH2D*)fMC[j]->Get(( "hresolutionbins_"+variables[i]).c_str());
      h_var->Scale(fMC_LumiWeight[j]);

      hs_var->Add(h_var);

    }

    hs_var->Write();

  }

  output_hists->Close();

  TString Filename = "UnfoldingHistosResolutionPlots/resolutionhists_fullRun2UL.root";
  TFile *File = new TFile(Filename.Data(), "READ");


  // polarizations
  string pol_variables[10] = {"b1n", "b2n", "b1r", "b2r", "b1k", "b2k", "b1j", "b2j", "b1q", "b2q"};

  int pol_reco_nbin[10] = {reco_nbbin, reco_nbbin, reco_nbbin, reco_nbbin, reco_nbbin, reco_nbbin, reco_nbbin, reco_nbbin, reco_nbbin, reco_nbbin};

  double* pol_reco_bins[10] = {reco_bbins, reco_bbins, reco_bbins, reco_bbins, reco_bbins, reco_bbins, reco_bbins, reco_bbins, reco_bbins, reco_bbins};

  gStyle->SetOptStat("MR");
  gStyle->SetOptFit(111);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.13);
  gStyle->SetStatY(1.0);
  
  TCanvas *c1 = new TCanvas("c1", "pol_Resolution",1600,1600);
  c1->Divide(2,5);

  TCanvas *c2 = new TCanvas("c2", "pol_Resolution_bins",1600,1600);
  c2->Divide(2,5);
  
  TCanvas *c3 = new TCanvas("c3", "pol_Bias_bins",1600,1600);
  c3->Divide(2,5);
  
  int i_canvas = 1;

  for (UInt_t i = 0; i < sizeof(pol_variables)/sizeof(pol_variables[0]); i++) {

    std::cout << pol_variables[i] << std::endl;

    TH2D* hresolutionbins = (TH2D*)File->Get(("hresolutionbins_"+pol_variables[i]).c_str());

    hresolutionbins->Rebin2D(1,4);

    //    hresolutionbins->GetXaxis()->SetTitle(ObsName(pol_variables[i],0)+" (reco-gen)");
    hresolutionbins->GetYaxis()->SetTitle("Entries");
    
    TH1D* hresolution = (TH1D*) hresolutionbins->ProjectionX();

    TH1D *hresmean;
    TH1D *hresRMS;
    TH1D *hresMeanFit;
    TH1D *hresSigmaFit;
    TH1D *hresMode;
    TH1D *hbinwidth;

    TH1::SetDefaultSumw2();

    hresmean = new TH1D(("hresmean_"+pol_variables[i]).c_str(), "; (gen);mean(reco-gen)",pol_reco_nbin[i],pol_reco_bins[i]);  
    hresRMS = new TH1D(("hresRMS_"+pol_variables[i]).c_str(), "; (gen);RMS(reco-gen)",pol_reco_nbin[i],pol_reco_bins[i]); 
    hresMeanFit = new TH1D(("hresMeanFit_"+pol_variables[i]).c_str(), "; (gen);MeanFit(reco-gen)",pol_reco_nbin[i],pol_reco_bins[i]);  
    hresSigmaFit = new TH1D(("hresSigmaFit_"+pol_variables[i]).c_str(), "; (gen);SigmaFit(reco-gen)",pol_reco_nbin[i],pol_reco_bins[i]);  
    hresMode = new TH1D(("hresMode_"+pol_variables[i]).c_str(), "; (gen);Mode(reco-gen)",pol_reco_nbin[i],pol_reco_bins[i]);  
    hbinwidth = new TH1D(("hbinwidth_"+pol_variables[i]).c_str(), Form(";%s (gen);Estimated resolution for %s",(ObsName(pol_variables[i],0)).Data(),(ObsName(pol_variables[i],0)).Data() ),pol_reco_nbin[i],pol_reco_bins[i]);  

    fillmeanandRMShistos(hresolutionbins, hresmean, hresRMS, hresMeanFit, hresSigmaFit, hresMode, hbinwidth, ObsName(pol_variables[i],0));

    drawResPlots(i_canvas, c1, c2, c3, hresolution, hresmean, hresRMS, hresMeanFit, hresSigmaFit, hresMode, hbinwidth);

  }

  c1->Print("UnfoldingHistosResolutionPlots/pol_Resolution.pdf");
  c1->Print("UnfoldingHistosResolutionPlots/pol_Resolution.root");
  c1->Print("UnfoldingHistosResolutionPlots/pol_Resolution.C");
  c1->Close();
    
  gStyle->SetOptStat("");
  
  c2->Print("UnfoldingHistosResolutionPlots/pol_Resolution_bins.pdf");
  c2->Print("UnfoldingHistosResolutionPlots/pol_Resolution_bins.root");
  c2->Print("UnfoldingHistosResolutionPlots/pol_Resolution_bins.C");
  c2->Close();
  
  c3->Print("UnfoldingHistosResolutionPlots/pol_Bias_bins.pdf");
  c3->Print("UnfoldingHistosResolutionPlots/pol_Bias_bins.root");
  c3->Print("UnfoldingHistosResolutionPlots/pol_Bias_bins.C");
  c3->Close();
  //

  // spin correlations`
  string corr_variables[10] = {"ll_cHel", "c_kk", "c_rr", "c_nn", "c_Prk", "c_Mrk", "c_Pnr", "c_Mnr", "c_Pnk", "c_Mnk"};

  int corr_reco_nbin[10] = {reco_ncosbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin};

  double* corr_reco_bins[10] = {reco_cosbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins};

  gStyle->SetOptStat("MR");
  gStyle->SetOptFit(111);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.13);
  gStyle->SetStatY(1.0);
  
  TCanvas *c4 = new TCanvas("c4", "corr_Resolution",1600,1600);
  c4->Divide(2,5);

  TCanvas *c5 = new TCanvas("c5", "corr_Resolution_bins",1600,1600);
  c5->Divide(2,5);
  
  TCanvas *c6 = new TCanvas("c6", "corr_Bias_bins",1600,1600);
  c6->Divide(2,5);
  
  i_canvas = 1;

  for (UInt_t i = 0; i < sizeof(corr_variables)/sizeof(corr_variables[0]); i++) {

    std::cout << corr_variables[i] << std::endl;

    TH2D* hresolutionbins = (TH2D*)File->Get(("hresolutionbins_"+corr_variables[i]).c_str());

    hresolutionbins->Rebin2D(1,4);

    hresolutionbins->GetXaxis()->SetTitle(ObsName(corr_variables[i],0)+" (reco-gen)");
    hresolutionbins->GetYaxis()->SetTitle("Entries");
    
    TH1D* hresolution = (TH1D*) hresolutionbins->ProjectionX();

    TH1D *hresmean;
    TH1D *hresRMS;
    TH1D *hresMeanFit;
    TH1D *hresSigmaFit;
    TH1D *hresMode;
    TH1D *hbinwidth;

    TH1::SetDefaultSumw2();

    hresmean = new TH1D(("hresmean_"+corr_variables[i]).c_str(), "; (gen);mean(reco-gen)",corr_reco_nbin[i],corr_reco_bins[i]);  
    hresRMS = new TH1D(("hresRMS_"+corr_variables[i]).c_str(), "; (gen);RMS(reco-gen)",corr_reco_nbin[i],corr_reco_bins[i]); 
    hresMeanFit = new TH1D(("hresMeanFit_"+corr_variables[i]).c_str(), "; (gen);MeanFit(reco-gen)",corr_reco_nbin[i],corr_reco_bins[i]);  
    hresSigmaFit = new TH1D(("hresSigmaFit_"+corr_variables[i]).c_str(), "; (gen);SigmaFit(reco-gen)",corr_reco_nbin[i],corr_reco_bins[i]);  
    hresMode = new TH1D(("hresMode_"+corr_variables[i]).c_str(), "; (gen);Mode(reco-gen)",corr_reco_nbin[i],corr_reco_bins[i]);  
    hbinwidth = new TH1D(("hbinwidth_"+corr_variables[i]).c_str(), Form(";%s (gen);Estimated resolution for %s",(ObsName(corr_variables[i],0)).Data(),(ObsName(corr_variables[i],0)).Data() ),corr_reco_nbin[i],corr_reco_bins[i]);  

    fillmeanandRMShistos(hresolutionbins, hresmean, hresRMS, hresMeanFit, hresSigmaFit, hresMode, hbinwidth, ObsName(corr_variables[i],0));

    drawResPlots(i_canvas, c4, c5, c6, hresolution, hresmean, hresRMS, hresMeanFit, hresSigmaFit, hresMode, hbinwidth);

  }

  c4->Print("UnfoldingHistosResolutionPlots/corr_Resolution.pdf");
  c4->Print("UnfoldingHistosResolutionPlots/corr_Resolution.root");
  c4->Print("UnfoldingHistosResolutionPlots/corr_Resolution.C");
  c4->Close();
    
  gStyle->SetOptStat("");
  
  c5->Print("UnfoldingHistosResolutionPlots/corr_Resolution_bins.pdf");
  c5->Print("UnfoldingHistosResolutionPlots/corr_Resolution_bins.root");
  c5->Print("UnfoldingHistosResolutionPlots/corr_Resolution_bins.C");
  c5->Close();
  
  c6->Print("UnfoldingHistosResolutionPlots/corr_Bias_bins.pdf");
  c6->Print("UnfoldingHistosResolutionPlots/corr_Bias_bins.root");
  c6->Print("UnfoldingHistosResolutionPlots/corr_Bias_bins.C");
  c6->Close();
  //

  // spin correlation linear combinations
  string lin_variables[11] = {"c_han", "c_sca", "c_tra", "c_kjL", "c_rqL", "c_rkP", "c_rkM", "c_nrP", "c_nrM", "c_nkP", "c_nkM"};

  int lin_reco_nbin[11] = {reco_ncosbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin, reco_ncbin};

  double* lin_reco_bins[11] = {reco_cosbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins, reco_cbins};

  gStyle->SetOptStat("MR");
  gStyle->SetOptFit(111);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.13);
  gStyle->SetStatY(1.0);
  
  TCanvas *c7 = new TCanvas("c7", "lin_Resolution",1600,1600);
  c7->Divide(2,6);

  TCanvas *c8 = new TCanvas("c8", "lin_Resolution_bins",1600,1600);
  c8->Divide(2,6);
  
  TCanvas *c9 = new TCanvas("c9", "lin_Bias_bins",1600,1600);
  c9->Divide(2,6);
  
  i_canvas = 1;

  for (UInt_t i = 0; i < sizeof(lin_variables)/sizeof(lin_variables[0]); i++) {

    std::cout << lin_variables[i] << std::endl;

    TH2D* hresolutionbins = (TH2D*)File->Get(("hresolutionbins_"+lin_variables[i]).c_str());

    hresolutionbins->Rebin2D(1,4);

    hresolutionbins->GetXaxis()->SetTitle(ObsName(lin_variables[i],0)+" (reco-gen)");
    hresolutionbins->GetYaxis()->SetTitle("Entries");
    
    TH1D* hresolution = (TH1D*) hresolutionbins->ProjectionX();

    TH1D *hresmean;
    TH1D *hresRMS;
    TH1D *hresMeanFit;
    TH1D *hresSigmaFit;
    TH1D *hresMode;
    TH1D *hbinwidth;

    TH1::SetDefaultSumw2();

    hresmean = new TH1D(("hresmean_"+lin_variables[i]).c_str(), "; (gen);mean(reco-gen)",lin_reco_nbin[i],lin_reco_bins[i]);  
    hresRMS = new TH1D(("hresRMS_"+lin_variables[i]).c_str(), "; (gen);RMS(reco-gen)",lin_reco_nbin[i],lin_reco_bins[i]); 
    hresMeanFit = new TH1D(("hresMeanFit_"+lin_variables[i]).c_str(), "; (gen);MeanFit(reco-gen)",lin_reco_nbin[i],lin_reco_bins[i]);  
    hresSigmaFit = new TH1D(("hresSigmaFit_"+lin_variables[i]).c_str(), "; (gen);SigmaFit(reco-gen)",lin_reco_nbin[i],lin_reco_bins[i]);  
    hresMode = new TH1D(("hresMode_"+lin_variables[i]).c_str(), "; (gen);Mode(reco-gen)",lin_reco_nbin[i],lin_reco_bins[i]);  
    hbinwidth = new TH1D(("hbinwidth_"+lin_variables[i]).c_str(), Form(";%s (gen);Estimated resolution for %s",(ObsName(lin_variables[i],0)).Data(),(ObsName(lin_variables[i],0)).Data() ),lin_reco_nbin[i],lin_reco_bins[i]);  

    fillmeanandRMShistos(hresolutionbins, hresmean, hresRMS, hresMeanFit, hresSigmaFit, hresMode, hbinwidth, ObsName(lin_variables[i],0));

    drawResPlots(i_canvas, c7, c8, c9, hresolution, hresmean, hresRMS, hresMeanFit, hresSigmaFit, hresMode, hbinwidth);

  }

  c7->Print("UnfoldingHistosResolutionPlots/lin_Resolution.pdf");
  c7->Print("UnfoldingHistosResolutionPlots/lin_Resolution.root");
  c7->Print("UnfoldingHistosResolutionPlots/lin_Resolution.C");
  c7->Close();
    
  gStyle->SetOptStat("");
  
  c8->Print("UnfoldingHistosResolutionPlots/lin_Resolution_bins.pdf");
  c8->Print("UnfoldingHistosResolutionPlots/lin_Resolution_bins.root");
  c8->Print("UnfoldingHistosResolutionPlots/lin_Resolution_bins.C");
  c8->Close();
  
  c9->Print("UnfoldingHistosResolutionPlots/lin_Bias_bins.pdf");
  c9->Print("UnfoldingHistosResolutionPlots/lin_Bias_bins.root");
  c9->Print("UnfoldingHistosResolutionPlots/lin_Bias_bins.C");
  c9->Close();
  //


  // 2D variables
  string TwoD_variables[2] = {"ttbar_mass", "top_scatteringangle_ttbarframe"};

  int TwoD_reco_nbin[2] = {reco_nttbar_massbin, reco_ntop_scatteringangle_ttbarframebin};

  double* TwoD_reco_bins[2] = {reco_ttbar_massbins, reco_top_scatteringangle_ttbarframebins};

  gStyle->SetOptStat("MR");
  gStyle->SetOptFit(111);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.13);
  gStyle->SetStatY(1.0);
  
  TCanvas *c10 = new TCanvas("c10", "TwoD_Resolution",1600,400);
  c10->Divide(2,1);

  TCanvas *c11 = new TCanvas("c11", "TwoD_Resolution_bins",1600,400);
  c11->Divide(2,1);
  
  TCanvas *c12 = new TCanvas("c12", "TwoD_Bias_bins",1600,400);
  c12->Divide(2,1);
  
  i_canvas = 1;

  for (UInt_t i = 0; i < sizeof(TwoD_variables)/sizeof(TwoD_variables[0]); i++) {

    std::cout << TwoD_variables[i] << std::endl;

    TH2D* hresolutionbins = (TH2D*)File->Get(("hresolutionbins_"+TwoD_variables[i]).c_str());

    hresolutionbins->Rebin2D(1,2);

    hresolutionbins->GetXaxis()->SetTitle(ObsName(TwoD_variables[i],0)+" (reco-gen)");
    hresolutionbins->GetYaxis()->SetTitle("Entries");
    
    TH1D* hresolution = (TH1D*) hresolutionbins->ProjectionX();

    TH1D *hresmean;
    TH1D *hresRMS;
    TH1D *hresMeanFit;
    TH1D *hresSigmaFit;
    TH1D *hresMode;
    TH1D *hbinwidth;

    TH1::SetDefaultSumw2();

    hresmean = new TH1D(("hresmean_"+TwoD_variables[i]).c_str(), "; (gen);mean(reco-gen)",TwoD_reco_nbin[i],TwoD_reco_bins[i]);  
    hresRMS = new TH1D(("hresRMS_"+TwoD_variables[i]).c_str(), "; (gen);RMS(reco-gen)",TwoD_reco_nbin[i],TwoD_reco_bins[i]); 
    hresMeanFit = new TH1D(("hresMeanFit_"+TwoD_variables[i]).c_str(), "; (gen);MeanFit(reco-gen)",TwoD_reco_nbin[i],TwoD_reco_bins[i]);  
    hresSigmaFit = new TH1D(("hresSigmaFit_"+TwoD_variables[i]).c_str(), "; (gen);SigmaFit(reco-gen)",TwoD_reco_nbin[i],TwoD_reco_bins[i]);  
    hresMode = new TH1D(("hresMode_"+TwoD_variables[i]).c_str(), "; (gen);Mode(reco-gen)",TwoD_reco_nbin[i],TwoD_reco_bins[i]);  
    hbinwidth = new TH1D(("hbinwidth_"+TwoD_variables[i]).c_str(), Form(";%s (gen);Estimated resolution for %s",(ObsName(TwoD_variables[i],0)).Data(),(ObsName(TwoD_variables[i],0)).Data() ),TwoD_reco_nbin[i],TwoD_reco_bins[i]);  

    fillmeanandRMShistos(hresolutionbins, hresmean, hresRMS, hresMeanFit, hresSigmaFit, hresMode, hbinwidth, ObsName(TwoD_variables[i],0));

    drawResPlots(i_canvas, c10, c11, c12, hresolution, hresmean, hresRMS, hresMeanFit, hresSigmaFit, hresMode, hbinwidth);

  }

  c10->Print("UnfoldingHistosResolutionPlots/TwoD_Resolution.pdf");
  c10->Print("UnfoldingHistosResolutionPlots/TwoD_Resolution.root");
  c10->Print("UnfoldingHistosResolutionPlots/TwoD_Resolution.C");
  c10->Close();
    
  gStyle->SetOptStat("");
  
  c11->Print("UnfoldingHistosResolutionPlots/TwoD_Resolution_bins.pdf");
  c11->Print("UnfoldingHistosResolutionPlots/TwoD_Resolution_bins.root");
  c11->Print("UnfoldingHistosResolutionPlots/TwoD_Resolution_bins.C");
  c11->Close();
  
  c12->Print("UnfoldingHistosResolutionPlots/TwoD_Bias_bins.pdf");
  c12->Print("UnfoldingHistosResolutionPlots/TwoD_Bias_bins.root");
  c12->Print("UnfoldingHistosResolutionPlots/TwoD_Bias_bins.C");
  c12->Close();
  //


}
