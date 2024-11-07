#include <TH1D.h>
#include <iostream>
//#include "../unfolding/TopSVDFunctions.h"

#ifndef TTBAR_FRAMEWORK_H
#define TTBAR_FRAMEWORK_H

// BACKGROUND HANDLING
// The raw data are given plus the simulated
// backgrounds from MC. The result is returned in
// the placeholder parameter 'dataHist'.
// Notice:
// (1) Errors from background will be neglected! 
// (2) If 'ttbgrInputHist' is given, it will be removed by
//     a multiplicative approach that is xSec independent!   
// (3) If 'bgrInputHist' is given, it will be removed from the
//     data by substraction.
// (4) If both 'bgrInputHist' and 'ttbgrInputHist' are given,
//     the difference between 'bgrInputHist' and 'ttbgrInputHist'
//     will be substracted from the data and the background given
//     int 'ttbgrInputHist' will be treated as described in (2).
// (5) If neither 'bgrInputHist' nor 'ttbgrInputHist' are given,
//     no background will be considered.
class TopSVDFunctions_ttmd
{
  public:
    // Background
    static void SVD_BackgrHandling(TH1D*& dataHist, TH1D* bgrHist, TH1D* ttbgrHist, TH1D* biniHist, TH1D* rawHist, bool beForgiving, int numHist = 1);
};

void TopSVDFunctions_ttmd::SVD_BackgrHandling(TH1D*& dataHist, TH1D* bgrHist, TH1D* ttbgrHist, TH1D* biniHist, TH1D* rawHist, bool beForgiving, int numHist)
{ 
  //beForgiving = 1;
  using namespace std;

  int flag_verbose = 1;
  double tolerance_factor = 20.0;


  // Steer background handling
  bool doBgr = false;
  if ( bgrHist != NULL ) doBgr = true;
  bool doTtBgr = false;
  if ( ttbgrHist != NULL ) doTtBgr = true;
  bool flagAddBkgStatUnc = 1;

  // Nbins
  int nbins = rawHist->GetNbinsX();

  // Loop over all histograms
  for ( int h = 0 ; h < numHist ; h++ ) {

    // loop bins
    for ( int i = 1 ; i <= nbins ; i++) {

      // get bin and errors for data
      // Notice: Only ONE histogram will be used here!
      double value_data = rawHist->GetBinContent(i);
      double err_data = rawHist->GetBinError(i);
      // Get background value to be substracted!
      double value_bgr = 0.;
      if ( doBgr == true   ) value_bgr = (bgrHist+h)->GetBinContent(i);
      // if ttbar signal fraction method is used: consider only non ttbar BG in check
      double relevant_bgr = value_bgr;
      if ( doTtBgr == true ) relevant_bgr-=(ttbgrHist+h)->GetBinContent(i);
      // check for data to be smaller as BG that will be subtracted
      if(value_data<relevant_bgr){
        if ( beForgiving == false ) {
          // 18.02.18 adding different channels: low statistics, do reasonable check that the difference is within 4 data stat. unc.
          // TODO tune factor below?
	  // or data = 0 and background is small enough
      double tolerance = tolerance_factor * TMath::Sqrt(value_data);
          if((relevant_bgr - value_data) < tolerance || (value_data == 0.0 && relevant_bgr < 1.0))
          {
            printf("Warning in TopSVDFunctions::SVD_BackgrHandling(): bin %d (relevant_bgr[%.1f] - value_data[%.1f]) < tolerance[%.1f]\n", i, relevant_bgr, value_data, tolerance);
            // below via label relevant_bgr is not relevant
            value_bgr = value_data;
            if(doTtBgr)
              value_bgr -= (ttbgrHist+h)->GetBinContent(i);
            goto labelMakeDiffNull;
          }
          std::cout << "ERROR in TopSVDFunctions::SVD_BackgrHandling: " << std::endl;
          std::cout << "N_MC BG > N_data in bin " << i << std::endl;
          for ( int k = 1 ; k <= nbins ; k++ ) {
            std::cout << "   Bin " << k << ":  Data=" << rawHist->GetBinContent(k) << ",    Bgr=";
            double relevant_bgr_bin =(bgrHist+h)->GetBinContent(k);
            if ( doTtBgr == true ) relevant_bgr_bin-=(ttbgrHist+h)->GetBinContent(k);
            std::cout << relevant_bgr_bin;
            if ( relevant_bgr > value_data ) {
              cout << " !@#$%^&*! " << endl;
            }
            else {
              cout << endl;
            }
          }
          exit(1);
        } else {
          labelMakeDiffNull : ;
          relevant_bgr = value_data;
        }
      }

      // Debug output
      if(flag_verbose>=3){
        std::cout << std::endl << "bin " << i << std::endl;
        std::cout << "value_data: " << value_data << std::endl;
        std::cout << "value_bgr: "  << value_bgr << std::endl;
        std::cout << "err_data: "   << err_data << std::endl;
      }


      // Get background scale factor!
      double sigFrac = 1.;
      if ( doTtBgr == true ) {
        double value_ttSig=(biniHist+h)->GetBinContent(i);
        double value_ttBgr=(ttbgrHist+h)->GetBinContent(i);
        if(value_ttBgr>value_bgr){
          if ( beForgiving == false ) {
            std::cout << "ERROR in TopSVDFunctions::SVD_BackgrHandling: " << std::endl;
            std::cout << "The TtBar-Background is larger than the complete Background!" << std::endl;
            std::cout << "This happens in bin " << i << " (range " << rawHist->GetBinLowEdge(i) << ",";
            std::cout << rawHist->GetBinLowEdge(i+1) << " )" << std::endl;
            std::cout << "    TTBar-Bgr: " << value_ttBgr << std::endl;
            std::cout << "    Bgr:       " << value_bgr << std::endl;
            rawHist->Print("all");
            exit(0);
          } else {
            value_ttBgr=value_bgr;
          }
        }

        // don't subtract ttbar BG !!!!!!!
        value_bgr-=value_ttBgr;

        // calculate signal fraction
        if ( value_ttSig+value_ttBgr > 0. ) {
          sigFrac=value_ttSig/(value_ttSig+value_ttBgr);
        }

        // debug output
        if(flag_verbose>=3){
          std::cout << "value_ttBgr: " << value_ttBgr << std::endl;
          std::cout << "value_bgr(no ttbar): " << value_bgr << std::endl;
          std::cout << "value_ttSig: " << value_ttSig << std::endl;
          std::cout << "sigFrac: " << sigFrac << std::endl;
        }

      }

      // New Values
      double value_new = (value_data - value_bgr)*sigFrac;
      double err_new = err_data*sigFrac;
      if(flagAddBkgStatUnc && doBgr)
      {
        double err_bkg = bgrHist->GetBinError(i);
        double err_new_withbkg = TMath::Sqrt(err_data * err_data + err_bkg * err_bkg) * sigFrac;
        //printf("error increase with bkg: %.2f\n", err_new_withbkg / err_new);
        err_new = err_new_withbkg;
      }


      // debug output
      if(flag_verbose>=3){
        std::cout << "data_value_new: " << value_new << std::endl;
        std::cout << "data_err_new: "   << err_new   << std::endl;
      }

      // Save new values
      (dataHist+h)->SetBinContent(i, value_new);
      (dataHist+h)->SetBinError(i, err_new);

    }
  }
}

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

TH1D* SubtractBkg(const TH1D* dataHist, const TH1D* nonttbgrHist, const TH1D* ttbgrHist, const TH1D* ttHist)
{
  TH1D* hBgrTotal = (TH1D*) nonttbgrHist->Clone();
  // 18.02.18 implementing all dilepton channels: make sure all bins are non-negative
  for(int b = 0; b <= hBgrTotal->GetNbinsX(); b++)
    if(hBgrTotal->GetBinContent(b) < 0.0)
    {
      printf("Warning in SubtractBkg(): negative bin content (this is not-tt background) hBgrTotal->GetBinContent(%d) = %.1f -> resetting to 0\n", b, hBgrTotal->GetBinContent(b));
      hBgrTotal->SetBinContent(b, 0.0);
    }
  hBgrTotal->Add(ttbgrHist);
  TH1D* hDatNoBgr = (TH1D*) dataHist->Clone();
  bool beForgiving = false;
  int numHist = 1;

  // copy and later delete some histograms to keep the arguments const
  TH1D* ttbgrHistCopy = (TH1D*) ttbgrHist->Clone();
  TH1D* ttHistCopy = (TH1D*) ttHist->Clone();
  TH1D* dataHistCopy = (TH1D*) dataHist->Clone();

  TopSVDFunctions_ttmd::SVD_BackgrHandling(hDatNoBgr, hBgrTotal, ttbgrHistCopy, ttHistCopy, dataHistCopy, beForgiving, numHist);

  delete ttbgrHistCopy;
  delete ttHistCopy;
  delete hBgrTotal;
  delete dataHistCopy;

  return hDatNoBgr;
}

#endif
