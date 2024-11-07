/******************************************************************************************************
 * File        : MatrixUnfControl.C 
 * Authors     : Ajeeta Khatiwada, Andreas Jung, Jacob Linacre
 * Description : Reads in all different matrices from bg, signal and the data vector; applies unfolding
 ******************************************************************************************************/

#include "stdInc.h"
#include "TRandom3.h"
#include "THStack.h"
#include <vector>
#include <TSystem.h>
#include <iostream>
#include <fstream>
#include "TPaveText.h"
#include "TGaxis.h"
#include "TDecompSVD.h"

// For the unfolding!
#include "TUnfoldDensity.h"
#include "MatrixUnf.h"
#include "MatrixUnfControl.h"
#include "../../common/include/plotterUtils.h"

// Amandeep : For reading in XML binning
#include <TDOMParser.h>
#include <TXMLDocument.h>
#include "TUnfoldBinningXML.h"
// End

TRandom3 gMyRandom;

MatrixUnfControl::MatrixUnfControl()
{
  //fMCFile=NULL;
  fMCFiles=NULL;
  fBgFiles=NULL;
  fDataFiles=NULL;
  fSysFiles=NULL;  
  //fXSectionSig=-1.0;
  fDataArray=NULL;
  fRecArray=NULL;
  fGenArray=NULL;
  fVisGenArray=NULL;
  fMigMatArray=NULL;
  fBgArray=NULL;
  fBgSamplesArray=NULL;
  fMigMatSysArray=NULL;
  fTrueDataArray=NULL;
  fEra="";
  fLumi=-1.0;
  fTauMin=-1.0;
  fTauMax=-1.0;
  fTauSys=-1.0;
  fMinRhoAVG=kFALSE;
  fClosureTest = kFALSE;  
  fSelfConsistency = kFALSE; 
  fEnsembleTest = kFALSE;
  fLinearityTest = kFALSE;
  fAllowRebinOption=kFALSE;
  fCombinedSumptt=kFALSE;
  fConstrainTau=kFALSE;
  fSubtrBgMCFromData=kFALSE;
  fDoSubtrBg=kFALSE;
  fUseAreaConstrain=kFALSE;
  fUseBiasInUnfold=kFALSE;
  fRegMode=-1.0;
  fsetSyst=kFALSE;
}

// ***********************************************************
// DoAllUnfold :
// Sets up data, reco, gen and response matrices for unfolding
// Then calls the unfolding method from MatrixUnf.cc 
// ***********************************************************

void MatrixUnfControl::DoAllUnfold(TString PREUNFFILE, TString  TUNFFILE, TString CHAN){

  //---------------------- For the output --------------------------------//

  TString SummaryFilename           = TUNFFILE + "Linearity_summary.root";
  if(fEnsembleTest) SummaryFilename = TUNFFILE + "Bootstrap_deltapull_histos.root";
  TFile *summaryunfoutfile          = new TFile(SummaryFilename.Data(),"RECREATE");
  std::cout << " summaryunfoutfile" << std::endl;
  summaryunfoutfile->ls();
  summaryunfoutfile->cd();
  
  std::cout << "fVariableNames.size(): " << fVariableNames.size() << std::endl;

  // Number of points to scan beween coefficient = C_SM-0.5 and C_SM+0.5 
  // when running linearity tests (should be an odd number)
  const int nLinearity = 21;
  
  // Amandeep : Changing from 3000 to 10
  // Number of pseudoexperiments to create/read, when running pull tests
  const int nPE = 1;
  
  // Whether to create simple PEs from 1D poisson variations (false) or read in bootstrapped histograms (true)
  const bool doBootstrap = kTRUE;

  const int nVar         = fVariableNames.size();
  const int nbinsrhoi = 24; // Amandeep : Hardcoded here, what is the spirit of variable ?
  //const int nbinsrhoi    = 6;
  const int nbinsrn      = nbinsrhoi*nVar;

  // Vectors to store ensemble results for all variables
  std::vector<std::vector<std::vector<float>>> CoefPE_rebinnedA(nbinsrhoi, std::vector<std::vector<float>>(nVar, std::vector<float>(nPE)));
  std::vector<std::vector<std::vector<float>>> BinPE_rebinnedA(nbinsrhoi , std::vector<std::vector<float>>(nVar, std::vector<float>(nPE)));
  std::vector<std::vector<std::vector<float>>> CoefPE_rebinnedB(nbinsrhoi, std::vector<std::vector<float>>(nVar, std::vector<float>(nPE)));
  std::vector<std::vector<std::vector<float>>> BinPE_rebinnedB(nbinsrhoi , std::vector<std::vector<float>>(nVar, std::vector<float>(nPE)));

  TH1D* hAllVarUnf_rebinnedA;
  TH1D* hAllVarUnf_rebinnedB;

  TH1D* hist_acombfactor_rebinnedA[nbinsrhoi/2];
  TH1D* hist_acombfactor_rebinnedB[nbinsrhoi/2];

  if(doBootstrap && fEnsembleTest) {
    hAllVarUnf_rebinnedA = new TH1D("hAllVarUnf_rebinnedA", "hAllVarUnf_rebinnedA", nbinsrn,0,nbinsrn);
    hAllVarUnf_rebinnedB = new TH1D("hAllVarUnf_rebinnedB", "hAllVarUnf_rebinnedB", nbinsrn,0,nbinsrn);  
    
    for (int i_a = 0; i_a < nbinsrhoi/2; ++i_a)
    {
      hist_acombfactor_rebinnedA[i_a] = new TH1D(Form("hist_acombfactor_rebinnedA%d",i_a), Form("hist_acombfactor_rebinnedA%d",i_a), nVar, 0, nVar);
      hist_acombfactor_rebinnedB[i_a] = new TH1D(Form("hist_acombfactor_rebinnedB%d",i_a), Form("hist_acombfactor_rebinnedB%d",i_a), nVar, 0, nVar);
    }
  }

  if(fLinearityTest && !fEnsembleTest){
    // Amandeep : 
    // TODO : This is hardcoded
    // MatrixUnfControl::CreateLinearitySummaryHistos((int) nbinsrhoi / 6);
    MatrixUnfControl::CreateLinearitySummaryHistos((int) nbinsrhoi);
    // End
  }
  
  int nbgfiles  = 0;
  int ndatafiles= 0;
  int nMCfiles  = 0;
  int nSysFiles = 0;

  if(fBgFiles)   nbgfiles   = fBgFiles->GetEntries();
  if(fDataFiles) ndatafiles = fDataFiles->GetEntries();
  if(fMCFiles)   nMCfiles   = fMCFiles->GetEntries();
  if(fSysFiles) {
    nSysFiles  = fSysFiles->GetEntries();
    nSysFiles /= (nMCfiles > 0 ? nMCfiles : 1);
  }

  // Looping over the variables
  for(unsigned int i=0; i < fVariableNames.size(); i++){

    if(fVariableNames[i]=="") continue;  // histograms for TUnfold

    int n_sub_bins = 2;

    if(fAllowRebinOption) {
      std::ifstream in(Form("binning/%s_rebin.txt", fVariableNames[i].Data()));
      in >> n_sub_bins;
      in.close();
      std::cout << "Rebin: " << n_sub_bins << std::endl;
    }
    else n_sub_bins = 2;

    const int rebinfine = n_sub_bins;

    TDOMParser parser;
    Int_t error = parser.ParseFile(Form("binning/" + fVariableNames[i] + "_binning.xml"));
    
    if (error) std::cout << "error=" << error << " from TDOMParser\n";

    TXMLDocument const *XMLdocument     = parser.GetXMLDocument();
    TUnfoldBinningXML *generatorBinning = TUnfoldBinningXML::ImportXML(XMLdocument, "generator");
    TUnfoldBinningXML *detectorBinning  = TUnfoldBinningXML::ImportXML(XMLdocument, "detector");

    TString PreUnfoldedFilename = PREUNFFILE+fVariableNames[i]+".root";
    TFile *preunfoutfile        = new TFile(PreUnfoldedFilename.Data(),"RECREATE");
    // TFile *preunfoutfile     = new TFile(PreUnfoldedFilename.Data(),"READ");
    printf("Open file root-File %s\n\n", PreUnfoldedFilename.Data());                                                                             
    
    if(preunfoutfile==NULL) exit(0);                                                                                                                                           
       
    // Write the preunfolded object to file
    printf("\n All preunfolded - Writing output file ...\n"); 
    preunfoutfile->cd();
     
    ((TH1D*) ((*fDataArray)[i]))->Write();
    ((TH1D*) ((*fGenArray)[i]))->Write();
    ((TH1D*) ((*fVisGenArray)[i]))->Write();
    ((TH1D*) ((*fRecArray)[i]))->Write();
    ((TH1D*) ((*fBgArray)[i]))->Write();
    ((TH2D*) ((*fMigMatArray)[i]))->Write();
    for(int j=0; j<nSysFiles; j++){ ((TH2D*) ((*fMigMatSysArray)[i*nSysFiles+j]))->Write();}
    if(fClosureTest) ((TH1D*) ((*fTrueDataArray)[i]))->Write();
    
    preunfoutfile->Write();
    //preunfoutfile->Close();
    delete preunfoutfile;

    TH1D* hData;
    TH1D* hRecMC;
    TH1D* hGenMC;
    TH1D* hVisGenMC;
    TH2D* hMigMat;
    TH1D* hBgSubtrData;
    TH1D* hBgMC;
    TH1D* hTrueData;

    std::vector<TH1D*> vecBgMC;

    std::vector<TH2D*>vecMigMatSys;
    
    TString TUnfoldResultsFilename = Form(TUNFFILE + fVariableNames[i] + ".root");
    printf("Open file root-File %s\n\n", TUnfoldResultsFilename.Data());
    TFile *outfile = new TFile(TUnfoldResultsFilename.Data(),"RECREATE");                                                                                 
    if(outfile==NULL) exit(0);                                                                                                                                           
    outfile->cd();
    
    // Amandeep : Modifying since there is no preunfolded
    // Original :
    hData     = (TH1D*)((TH1D*)fDataArray->At(i))->Clone(Form("%s_norm",((TH1D*)fDataArray->At(i))->GetName()));
    hRecMC    = (TH1D*)((TH1D*)fRecArray->At(i))->Clone(Form("%s_norm",((TH1D*)fRecArray->At(i))->GetName()));
    hGenMC    = (TH1D*)((TH1D*)fGenArray->At(i))->Clone(Form("%s_norm",((TH1D*)fGenArray->At(i))->GetName()));
    hVisGenMC = (TH1D*)((TH1D*)fVisGenArray->At(i))->Clone(Form("%s_norm",((TH1D*)fVisGenArray->At(i))->GetName()));
    hMigMat   = (TH2D*)((TH2D*)fMigMatArray->At(i))->Clone(Form("%s_norm",((TH2D*)fMigMatArray->At(i))->GetName()));
    hBgMC     = (TH1D*)((TH1D*)fBgArray->At(i))->Clone(Form("%s_norm",((TH1D*)fBgArray->At(i))->GetName()));

    // Modified :
    // hData     = (TH1D*)((TH1D*)preunfoutfile->Get(Form("%sData", fVariableNames[i])))->Clone(Form("%sData_norm", fVariableNames[i]));
    // hRecMC    = (TH1D*)((TH1D*)preunfoutfile->Get(Form("%sReco", fVariableNames[i])))->Clone(Form("%sReco_norm", fVariableNames[i]));
    // hGenMC    = (TH1D*)((TH1D*)preunfoutfile->Get(Form("%sGen" , fVariableNames[i])))->Clone(Form("%sGen_norm", fVariableNames[i]);
    // hVisGenMC = (TH1D*)((TH1D*)preunfoutfile->Get(Form("%sVisGen" , fVariableNames[i])))->Clone(Form("%sVisGen_norm", fVariableNames[i]));
    // hMigMat   = (TH2D*)((TH2D*)preunfoutfile->Get(Form("%sRespMat", fVariableNames[i])))->Clone(Form("%sRespMat_norm", fVariableNames[i]));
    // hBgMC     = (TH1D*)((TH1D*)preunfoutfile->Get(Form("%sRecoBg" , fVariableNames[i])))->Clone(Form("%sRecoBg_norm", fVariableNames[i]));

    TString quantity = MatrixUnf::ObsName(fVariableNames[i],kFALSE);

    hRecMC->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    hGenMC->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    hVisGenMC->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    hMigMat->GetXaxis()->SetTitle(Form("%s (Gen)",quantity.Data()));
    hMigMat->GetYaxis()->SetTitle(Form("%s (Reco)",quantity.Data()));
    hBgMC->GetXaxis()->SetTitle(Form("%s",quantity.Data()));

    for(int j=0; j < nbgfiles; j++){
      // TH1D* hBgSample = (TH1D*)((TH1D*)fBgSamplesArray->At(i*nbgfiles+j))->Clone(Form("%s_norm",((TH1D*)fBgSamplesArray->At(i*nbgfiles+j))->GetName()));
      TH1D* hBgSample    = (TH1D*)fBgSamplesArray->At(i*nbgfiles+j);
      vecBgMC.push_back(hBgSample);
      std::cout << "Bg sample name :: " << j << " " << hBgSample->GetName() << " " << vecBgMC.at(j)->Integral() << std::endl;
    }

    // Check the sum of vecBgMC is consistent with the total, hBgMC
    double bkg_integral_check = 0;
    for(int j=0; j < nbgfiles; j++){
      bkg_integral_check += vecBgMC.at(j)->Integral();
    }
    if( fabs( 1. - bkg_integral_check/hBgMC->Integral() ) > 1e-10 ) {
      std::cout << "Error in background subtraction." << std::endl;
      exit(0);
    }

    for(int j=0; j < nSysFiles; j++){
      TH2D* hMigMatSys = (TH2D*)((TH2D*)fMigMatSysArray->At(i*nSysFiles+j))->Clone(Form("%s_norm",((TH2D*)fMigMatSysArray->At(i*nSysFiles+j))->GetName()));
      vecMigMatSys.push_back(hMigMatSys);
    }    

    // Amandeep : What does this do then ?
    if(fClosureTest){
      hTrueData = (TH1D*)((TH1D*)fTrueDataArray->At(i))->Clone(Form("%s_norm",((TH1D*)fTrueDataArray->At(i))->GetName()));
    }
    
    // Define bg subtracted data (or ensemble data) histos
    double tau_min,tau_max;
    
    if(!fConstrainTau){
      tau_min=fTauMin;
      tau_max=fTauMax;
    }
    else{
      tau_min=fTauSys;
      tau_max=fTauSys+(5e-10);
    }

    if(tau_min<0.0 || tau_max<0.0 || tau_min>tau_max){
      std::cout<<"Set different range for tau!!"<<std::endl;
      exit(0);
    }

    Bool_t is2DHist      = kFALSE;
    MatrixUnf *pttMatUnf = new MatrixUnf(fVariableNames[i],fVariableNames[i]+"(gen)[GeV]",fVariableNames[i],fVariableNames[i]+"(rec)[GeV]","Entries",Form("binning/%s_binning.xml",fVariableNames[i].Data()),is2DHist,fUseBiasInUnfold, CHAN, fsetSyst, hGenMC, rebinfine);

    if(fSelfConsistency) {
      // *******
      // CLOSURE
      // *******

      TH1D* hPseudoData = (TH1D*)hRecMC->Clone(Form("hPseudoData_%s",fVariableNames[i].Data()));

      if(fSubtrBgMCFromData || fDoSubtrBg) {
        hPseudoData->Add(hBgMC);
      }

      for(int x_bin=1; x_bin<=hPseudoData->GetNbinsX(); x_bin++){
        hPseudoData->SetBinError(x_bin,TMath::Sqrt(hPseudoData->GetBinContent(x_bin)));
      }

      if(fSubtrBgMCFromData && !fDoSubtrBg) {
        hPseudoData->Add(hBgMC,-1.); // TODO: should really use regularBgSubtMethod here
      } 

      unfoldData(pttMatUnf, fVariableNames[i], hPseudoData, hMigMat, hGenMC, hVisGenMC, hRecMC, vecBgMC, tau_min,tau_max, kFALSE, CHAN, vecMigMatSys);
    } 
    
    else {
      // **************
      // DATA UNFOLDING
      // **************

      if(fSubtrBgMCFromData && !fDoSubtrBg){
        std::cout<<"Subtracting background"<<std::endl;

        regularBgSubtMethod(hBgSubtrData, hBgMC, hData);

        unfoldData(pttMatUnf, fVariableNames[i], hBgSubtrData, hMigMat, hGenMC, hVisGenMC, hRecMC, vecBgMC, tau_min,tau_max, kFALSE, CHAN, vecMigMatSys);

      } 
      else { // if subtrBgMCFromData finished
        unfoldData(pttMatUnf, fVariableNames[i], hData, hMigMat, hGenMC, hVisGenMC, hRecMC, vecBgMC, tau_min,tau_max, kFALSE, CHAN, vecMigMatSys);
      }
    }
    
    // No self consistency check
    error = parser.ParseFile(Form("binning/" + fVariableNames[i] + "_binning_rebinnedB.xml"));
    if (error) std::cout << "error=" << error << " from TDOMParser\n";

    XMLdocument = parser.GetXMLDocument();
    TUnfoldBinningXML *generatorRebinnedBinning = TUnfoldBinningXML::ImportXML(XMLdocument, "generator_rebinnedB");
    TUnfoldBinningXML *detectorRebinnedBinning  = TUnfoldBinningXML::ImportXML(XMLdocument, "detector_rebinnedB");

    if (!generatorRebinnedBinning)
      std::cout << "could not read 'generator_rebinnedB' binning\n";

    if (!detectorRebinnedBinning)
      std::cout << "could not read 'detector_rebinnedB' binning\n";

    const TUnfoldBinning *binning = generatorRebinnedBinning->FindNode("ttbargen_rebinnedB");
    const int numDimensions       = binning->GetDistributionDimension();
    
    // ***********************************************************************    
    // ***********************************************************************
    // The remainder of this loop is for the optional linearity and pull tests
    // ***********************************************************************
    // ***********************************************************************
    
    const int unf_nbins               = hGenMC->GetNbinsX();
    const int unf_nbins_rebinnedS     = unf_nbins / pow(rebinfine, pttMatUnf->GetDim());
    const int numBinsPerCoefficient   = 6;
    const int numRebinnedCoefficients = unf_nbins_rebinnedS / numBinsPerCoefficient;
    const int numCoefficients         = unf_nbins / (numBinsPerCoefficient * rebinfine);

    TH1D *hGenMC_rebinnedS = (TH1D*)hGenMC->Clone(Form("hGenMC_rebinnedS_%s",fVariableNames[i].Data()));
    
    if (fAllowRebinOption) {
      MatrixUnf::rebinMultidimensionalInput(hGenMC_rebinnedS,
                                            generatorRebinnedBinning->FindNode("ttbargen_rebinnedB"),
                                            generatorBinning->FindNode("ttbargen"),
                                            rebinfine);
      //hGenMC_rebinnedS->Rebin(rebinfine);
    }


    // ***************
    // Linearity Tests
    // ***************

    if(fLinearityTest && !fEnsembleTest){
      std::cout<<"Performing linearity test"<<std::endl;

      TH1D* hLinearity[nLinearity];
      TH1D* hLinearityGen[nLinearity];
      TH1D* hLinearityGen_rebinnedS[nLinearity];
      TH1D* hAfbLinearityGen_rebinnedS[nLinearity];
      double delta_coefficient[nLinearity];

      // create histograms to store linearity results for each variable

      TH1D *hLinearityDBias = (TH1D*)hGenMC->Clone(Form("hLinearityDBias_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityCBias = (TH1D*)hGenMC->Clone(Form("hLinearityCBias_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDInjectedShape = (TH1D*)hGenMC->Clone(Form("hLinearityDInjectedShape_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDMeasuredShape = (TH1D*)hGenMC->Clone(Form("hLinearityDMeasuredShape_%s",fVariableNames[i].Data())); 
      TH1D *hLinearitySlope          = (TH1D*)hGenMC->Clone(Form("hLinearitySlope_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDBiasOverStatErr  = (TH1D*)hGenMC->Clone(Form("hLinearityDBiasOverStatErr_%s",fVariableNames[i].Data()));
      TH1D *hLinearityCBiasOverStatErr  = (TH1D*)hGenMC->Clone(Form("hLinearityCBiasOverStatErr_%s",fVariableNames[i].Data()));
      TH1D *hLinearityDBiasFractional   = (TH1D*)hGenMC->Clone(Form("hLinearityDBiasFractional_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityCBiasFractional   = (TH1D*)hGenMC->Clone(Form("hLinearityCBiasFractional_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDBiasMax          = (TH1D*)hGenMC->Clone(Form("hLinearityDBiasMax_%s",fVariableNames[i].Data())); 
      //unused TH1D *hLinearityCBiasMax = (TH1D*)hGenMC->Clone(Form("hLinearityCBiasMax_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDInjectedShapeMax = (TH1D*)hGenMC->Clone(Form("hLinearityDInjectedShapeMax_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDMeasuredShapeMax = (TH1D*)hGenMC->Clone(Form("hLinearityDMeasuredShapeMax_%s",fVariableNames[i].Data())); 
      TH1D *hLinearitySlopeMax = (TH1D*)hGenMC->Clone(Form("hLinearitySlopeMax_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDBiasOverStatErrMax = (TH1D*)hGenMC->Clone(Form("hLinearityDBiasOverStatErrMax_%s",fVariableNames[i].Data()));
      //unused TH1D *hLinearityCBiasOverStatErrMax = (TH1D*)hGenMC->Clone(Form("hLinearityCBiasOverStatErrMax_%s",fVariableNames[i].Data()));
      TH1D *hLinearityDBiasFractionalMax  = (TH1D*)hGenMC->Clone(Form("hLinearityDBiasFractionalMax_%s",fVariableNames[i].Data())); 
      //unused TH1D *hLinearityCBiasFractionalMax = (TH1D*)hGenMC->Clone(Form("hLinearityCBiasFractionalMax_%s",fVariableNames[i].Data())); 


      TH1D *hLinearityDBias_rebinnedA = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDBias_rebinnedA_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityCBias_rebinnedA = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityCBias_rebinnedA_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDInjectedShape_rebinnedA = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDInjectedShape_rebinnedA_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDMeasuredShape_rebinnedA = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDMeasuredShape_rebinnedA_%s",fVariableNames[i].Data())); 
      TH1D *hLinearitySlope_rebinnedA          = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearitySlope_rebinnedA_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDBiasOverStatErr_rebinnedA  = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDBiasOverStatErr_rebinnedA_%s",fVariableNames[i].Data()));
      TH1D *hLinearityCBiasOverStatErr_rebinnedA  = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityCBiasOverStatErr_rebinnedA_%s",fVariableNames[i].Data()));
      TH1D *hLinearityDBiasFractional_rebinnedA   = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDBiasFractional_rebinnedA_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityCBiasFractional_rebinnedA   = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityCBiasFractional_rebinnedA_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDBiasMax_rebinnedA          = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDBiasMax_rebinnedA_%s",fVariableNames[i].Data())); 
      //unused TH1D *hLinearityCBiasMax_rebinnedA = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityCBiasMax_rebinnedA_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDInjectedShapeMax_rebinnedA = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDInjectedShapeMax_rebinnedA_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDMeasuredShapeMax_rebinnedA = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDMeasuredShapeMax_rebinnedA_%s",fVariableNames[i].Data())); 
      TH1D *hLinearitySlopeMax_rebinnedA = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearitySlopeMax_rebinnedA_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDBiasOverStatErrMax_rebinnedA = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDBiasOverStatErrMax_rebinnedA_%s",fVariableNames[i].Data()));
      //unused TH1D *hLinearityCBiasOverStatErrMax_rebinnedA = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityCBiasOverStatErrMax_rebinnedA_%s",fVariableNames[i].Data()));
      TH1D *hLinearityDBiasFractionalMax_rebinnedA  = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDBiasFractionalMax_rebinnedA_%s",fVariableNames[i].Data()));
      //unused TH1D *hLinearityCBiasFractionalMax_rebinnedA = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityCBiasFractionalMax_rebinnedA_%s",fVariableNames[i].Data()));


      TH1D *hLinearityDBias_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDBias_rebinnedB_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityCBias_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityCBias_rebinnedB_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDInjectedShape_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDInjectedShape_rebinnedB_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDMeasuredShape_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDMeasuredShape_rebinnedB_%s",fVariableNames[i].Data())); 
      TH1D *hLinearitySlope_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearitySlope_rebinnedB_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDBiasOverStatErr_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDBiasOverStatErr_rebinnedB_%s",fVariableNames[i].Data()));
      TH1D *hLinearityCBiasOverStatErr_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityCBiasOverStatErr_rebinnedB_%s",fVariableNames[i].Data()));
      TH1D *hLinearityDBiasFractional_rebinnedB  = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDBiasFractional_rebinnedB_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityCBiasFractional_rebinnedB  = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityCBiasFractional_rebinnedB_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDBiasMax_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDBiasMax_rebinnedB_%s",fVariableNames[i].Data())); 
      //unused TH1D *hLinearityCBiasMax_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityCBiasMax_rebinnedB_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDInjectedShapeMax_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDInjectedShapeMax_rebinnedB_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDMeasuredShapeMax_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDMeasuredShapeMax_rebinnedB_%s",fVariableNames[i].Data())); 
      TH1D *hLinearitySlopeMax_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearitySlopeMax_rebinnedB_%s",fVariableNames[i].Data())); 
      TH1D *hLinearityDBiasOverStatErrMax_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDBiasOverStatErrMax_rebinnedB_%s",fVariableNames[i].Data()));
      //unused TH1D *hLinearityCBiasOverStatErrMax_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityCBiasOverStatErrMax_rebinnedB_%s",fVariableNames[i].Data()));
      TH1D *hLinearityDBiasFractionalMax_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityDBiasFractionalMax_rebinnedB_%s",fVariableNames[i].Data()));
      //unused TH1D *hLinearityCBiasFractionalMax_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hLinearityCBiasFractionalMax_rebinnedB_%s",fVariableNames[i].Data()));

      TH1D *hAfbLinearityDBias_rebinnedA = new TH1D(Form("hAfbLinearityDBias_rebinnedA_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityCBias_rebinnedA = new TH1D(Form("hAfbLinearityCBias_rebinnedA_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityDInjectedShape_rebinnedA = new TH1D(Form("hAfbLinearityDInjectedShape_rebinnedA_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityDMeasuredShape_rebinnedA = new TH1D(Form("hAfbLinearityDMeasuredShape_rebinnedA_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearitySlope_rebinnedA = new TH1D(Form("hAfbLinearitySlope_rebinnedA_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityDBiasOverStatErr_rebinnedA = new TH1D(Form("hAfbLinearityDBiasOverStatErr_rebinnedA_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityCBiasOverStatErr_rebinnedA = new TH1D(Form("hAfbLinearityCBiasOverStatErr_rebinnedA_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityDBiasFractional_rebinnedA  = new TH1D(Form("hAfbLinearityDBiasFractional_rebinnedA_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityCBiasFractional_rebinnedA  = new TH1D(Form("hAfbLinearityCBiasFractional_rebinnedA_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityDBiasMax_rebinnedA = new TH1D(Form("hAfbLinearityDBiasMax_rebinnedA_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      // unused TH1D *hAfbLinearityCBiasMax_rebinnedA = (TH1D*)hGenMC_rebinnedS->Clone(Form("hAfbLinearityCBiasMax_rebinnedA_%s",fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityDInjectedShapeMax_rebinnedA = new TH1D(Form("hAfbLinearityDInjectedShapeMax_rebinnedA_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityDMeasuredShapeMax_rebinnedA = new TH1D(Form("hAfbLinearityDMeasuredShapeMax_rebinnedA_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearitySlopeMax_rebinnedA = new TH1D(Form("hAfbLinearitySlopeMax_rebinnedA_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityDBiasOverStatErrMax_rebinnedA = new TH1D(Form("hAfbLinearityDBiasOverStatErrMax_rebinnedA_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      // unused TH1D *hAfbLinearityCBiasOverStatErrMax_rebinnedA = (TH1D*)hGenMC_rebinnedS->Clone(Form("hAfbLinearityCBiasOverStatErrMax_rebinnedA_%s",fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityDBiasFractionalMax_rebinnedA = new TH1D(Form("hAfbLinearityDBiasFractionalMax_rebinnedA_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      // unused TH1D *hAfbLinearityCBiasFractionalMax_rebinnedA = (TH1D*)hGenMC_rebinnedS->Clone(Form("hAfbLinearityCBiasFractionalMax_rebinnedA_%s",fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);

      TH1D *hAfbLinearityDBias_rebinnedB = new TH1D(Form("hAfbLinearityDBias_rebinnedB_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityCBias_rebinnedB = new TH1D(Form("hAfbLinearityCBias_rebinnedB_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityDInjectedShape_rebinnedB = new TH1D(Form("hAfbLinearityDInjectedShape_rebinnedB_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityDMeasuredShape_rebinnedB = new TH1D(Form("hAfbLinearityDMeasuredShape_rebinnedB_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearitySlope_rebinnedB = new TH1D(Form("hAfbLinearitySlope_rebinnedB_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityDBiasOverStatErr_rebinnedB = new TH1D(Form("hAfbLinearityDBiasOverStatErr_rebinnedB_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityCBiasOverStatErr_rebinnedB = new TH1D(Form("hAfbLinearityCBiasOverStatErr_rebinnedB_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityDBiasFractional_rebinnedB  = new TH1D(Form("hAfbLinearityDBiasFractional_rebinnedB_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityCBiasFractional_rebinnedB  = new TH1D(Form("hAfbLinearityCBiasFractional_rebinnedB_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityDBiasMax_rebinnedB = new TH1D(Form("hAfbLinearityDBiasMax_rebinnedB_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      // unused TH1D *hAfbLinearityCBiasMax_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hAfbLinearityCBiasMax_rebinnedB_%s",fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityDInjectedShapeMax_rebinnedB = new TH1D(Form("hAfbLinearityDInjectedShapeMax_rebinnedB_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityDMeasuredShapeMax_rebinnedB = new TH1D(Form("hAfbLinearityDMeasuredShapeMax_rebinnedB_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearitySlopeMax_rebinnedB = new TH1D(Form("hAfbLinearitySlopeMax_rebinnedB_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityDBiasOverStatErrMax_rebinnedB = new TH1D(Form("hAfbLinearityDBiasOverStatErrMax_rebinnedB_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      // unused TH1D *hAfbLinearityCBiasOverStatErrMax_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hAfbLinearityCBiasOverStatErrMax_rebinnedB_%s",fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      TH1D *hAfbLinearityDBiasFractionalMax_rebinnedB = new TH1D(Form("hAfbLinearityDBiasFractionalMax_rebinnedB_%s", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      // unused TH1D *hAfbLinearityCBiasFractionalMax_rebinnedB = (TH1D*)hGenMC_rebinnedS->Clone(Form("hAfbLinearityCBiasFractionalMax_rebinnedB_%s",fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);

      //clear results vectors to prevent crash when looping over variables
      fPseudoUnfResults.clear();
      fPseudoUnfResultsRebinnedA.clear();
      fPseudoUnfResultsRebinnedB.clear();
      fPseudoUnfAfbResultsRebinnedA.clear();
      fPseudoUnfAfbResultsRebinnedB.clear();
      fPseudoUnfOtherResultsRebinnedA.clear();
      fPseudoUnfOtherResultsRebinnedB.clear();

      // Create and unfold linearity test histos
      TH2D *hLinearityInputMig = (TH2D*)hMigMat->Clone(Form("hLinearityInputMig_%s",fVariableNames[i].Data()));

      for (int i_linearity = 0; i_linearity < nLinearity; i_linearity++) {
        double wLinearity = nLinearity-1;
        delta_coefficient[i_linearity] = (i_linearity - wLinearity/2. )/wLinearity;
        MakeLinearityHistos(hLinearityInputMig, hLinearity[i_linearity], hLinearityGen[i_linearity], fVariableNames[i], delta_coefficient[i_linearity], generatorRebinnedBinning, generatorBinning, numCoefficients, rebinfine, kFALSE);
        unfoldData(pttMatUnf, fVariableNames[i]+Form("_linearity%d",i_linearity), hLinearity[i_linearity], hMigMat, hGenMC, hVisGenMC, hRecMC, vecBgMC, tau_min,tau_max, kTRUE, CHAN, vecMigMatSys);

        hLinearityGen_rebinnedS[i_linearity] = (TH1D*)hLinearityGen[i_linearity]->Clone(Form("%s_linearity%d_ReweightedGen_rebinnedS",fVariableNames[i].Data(),i_linearity));
        if (fAllowRebinOption) {
          MatrixUnf::rebinMultidimensionalInput(hLinearityGen_rebinnedS[i_linearity],
                                    generatorRebinnedBinning->FindNode("ttbargen_rebinnedB"),
                                    generatorBinning->FindNode("ttbargen"),
                                    rebinfine);
          //hLinearityGen_rebinnedS[i_linearity]->Rebin(rebinfine);
        }

        // Create gen-level hiso for asymmetries and coefficient
        hAfbLinearityGen_rebinnedS[i_linearity] = new TH1D(Form("%s_Afblinearity%d_ReweightedGen_rebinnedS", fVariableNames[i].Data(), i_linearity),
                                                          hLinearityGen_rebinnedS[i_linearity]->GetTitle(),
                                                          (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2,
                                                          0,
                                                          1);
        //hAfbLinearityGen_rebinnedS[i_linearity] = (TH1D *)hLinearityGen_rebinnedS[i_linearity]->Clone(Form("%s_Afblinearity%d_ReweightedGen_rebinnedS", fVariableNames[i].Data(), i_linearity));
        MatrixUnf::fillAsymsAndCoefficient(hLinearityGen_rebinnedS[i_linearity], nullptr, hGenMC_rebinnedS, fVariableNames[i], hAfbLinearityGen_rebinnedS[i_linearity], nullptr, true, true, numRebinnedCoefficients);


        // Convert asymmetry to coefficient
        for (int i_coefficient = 0; i_coefficient < numRebinnedCoefficients; i_coefficient++) {
          double offset = i_coefficient * (numBinsPerCoefficient + 2);
          std::vector factor_in = {1., 1., 1.};
          double Afactor = MatrixUnf::AtoCfactor(fVariableNames[i], hLinearityGen_rebinnedS[i_linearity], factor_in, numRebinnedCoefficients, i_coefficient);
          hAfbLinearityGen_rebinnedS[i_linearity]->SetBinContent(offset + unf_nbins_rebinnedS, Afactor * hAfbLinearityGen_rebinnedS[i_linearity]->GetBinContent(offset + unf_nbins_rebinnedS) );
          hAfbLinearityGen_rebinnedS[i_linearity]->SetBinError(offset + unf_nbins_rebinnedS, Afactor * hAfbLinearityGen_rebinnedS[i_linearity]->GetBinError(offset + unf_nbins_rebinnedS) );
          fPseudoUnfAfbResultsRebinnedA[i_linearity]->SetBinContent(offset + unf_nbins_rebinnedS, Afactor * fPseudoUnfAfbResultsRebinnedA[i_linearity]->GetBinContent(offset + unf_nbins_rebinnedS) );
          fPseudoUnfAfbResultsRebinnedA[i_linearity]->SetBinError(offset + unf_nbins_rebinnedS, Afactor * fPseudoUnfAfbResultsRebinnedA[i_linearity]->GetBinError(offset + unf_nbins_rebinnedS) );
          fPseudoUnfAfbResultsRebinnedB[i_linearity]->SetBinContent(offset + unf_nbins_rebinnedS, Afactor * fPseudoUnfAfbResultsRebinnedB[i_linearity]->GetBinContent(offset + unf_nbins_rebinnedS) );
          fPseudoUnfAfbResultsRebinnedB[i_linearity]->SetBinError(offset + unf_nbins_rebinnedS, Afactor * fPseudoUnfAfbResultsRebinnedB[i_linearity]->GetBinError(offset + unf_nbins_rebinnedS) );
        }
      }
      std::cout<<"Finished creating linearity test input distributions"<<std::endl;

      // Fill linearity plots
      fillLinearityGraphs(fPseudoUnfResults, 1, false, "", nLinearity, i, unf_nbins, hLinearityGen, delta_coefficient, hLinearityDBiasFractional, hLinearityDBiasFractionalMax, hLinearityDBias, hLinearityDBiasMax, hLinearityDInjectedShape, hLinearityDInjectedShapeMax, hLinearityDMeasuredShape, hLinearityDMeasuredShapeMax, hLinearitySlope, hLinearitySlopeMax, hLinearityDBiasOverStatErr, hLinearityDBiasOverStatErrMax, hLinearityCBias, hLinearityCBiasFractional, hLinearityCBiasOverStatErr, hLinearityMaxBiasFrac, hLinearityMaxBiasFracStat, hLinearityMaxBiasSlope, hLinearityMaxBiasFracMax, hLinearityMaxBiasFracStatMax, hLinearityMaxBiasSlopeMax, hLinearityAvgBiasFrac, hLinearityAvgBiasFracStat, hLinearityAvgBiasSlope, nullptr, nullptr, numCoefficients);
      fillLinearityGraphs(fPseudoUnfResultsRebinnedA, rebinfine, false, "_rebinnedA", nLinearity, i, unf_nbins_rebinnedS, hLinearityGen_rebinnedS, delta_coefficient, hLinearityDBiasFractional_rebinnedA, hLinearityDBiasFractionalMax_rebinnedA, hLinearityDBias_rebinnedA, hLinearityDBiasMax_rebinnedA, hLinearityDInjectedShape_rebinnedA, hLinearityDInjectedShapeMax_rebinnedA, hLinearityDMeasuredShape_rebinnedA, hLinearityDMeasuredShapeMax_rebinnedA, hLinearitySlope_rebinnedA, hLinearitySlopeMax_rebinnedA, hLinearityDBiasOverStatErr_rebinnedA, hLinearityDBiasOverStatErrMax_rebinnedA, hLinearityCBias_rebinnedA, hLinearityCBiasFractional_rebinnedA, hLinearityCBiasOverStatErr_rebinnedA, hLinearityMaxBiasFrac_rebinnedA, hLinearityMaxBiasFracStat_rebinnedA, hLinearityMaxBiasSlope_rebinnedA, hLinearityMaxBiasFracMax_rebinnedA, hLinearityMaxBiasFracStatMax_rebinnedA, hLinearityMaxBiasSlopeMax_rebinnedA, hLinearityAvgBiasFrac_rebinnedA, hLinearityAvgBiasFracStat_rebinnedA, hLinearityAvgBiasSlope_rebinnedA, nullptr, nullptr, numRebinnedCoefficients);
      fillLinearityGraphs(fPseudoUnfResultsRebinnedB, rebinfine, false, "_rebinnedB", nLinearity, i, unf_nbins_rebinnedS, hLinearityGen_rebinnedS, delta_coefficient, hLinearityDBiasFractional_rebinnedB, hLinearityDBiasFractionalMax_rebinnedB, hLinearityDBias_rebinnedB, hLinearityDBiasMax_rebinnedB, hLinearityDInjectedShape_rebinnedB, hLinearityDInjectedShapeMax_rebinnedB, hLinearityDMeasuredShape_rebinnedB, hLinearityDMeasuredShapeMax_rebinnedB, hLinearitySlope_rebinnedB, hLinearitySlopeMax_rebinnedB, hLinearityDBiasOverStatErr_rebinnedB, hLinearityDBiasOverStatErrMax_rebinnedB, hLinearityCBias_rebinnedB, hLinearityCBiasFractional_rebinnedB, hLinearityCBiasOverStatErr_rebinnedB, hLinearityMaxBiasFrac_rebinnedB, hLinearityMaxBiasFracStat_rebinnedB, hLinearityMaxBiasSlope_rebinnedB, hLinearityMaxBiasFracMax_rebinnedB, hLinearityMaxBiasFracStatMax_rebinnedB, hLinearityMaxBiasSlopeMax_rebinnedB, hLinearityAvgBiasFrac_rebinnedB, hLinearityAvgBiasFracStat_rebinnedB, hLinearityAvgBiasSlope_rebinnedB, nullptr, nullptr, numRebinnedCoefficients);

      bool calculateCoefficients = true;
      if(calculateCoefficients){
        fillLinearityGraphs(fPseudoUnfAfbResultsRebinnedA, rebinfine, true, "_rebinnedA", nLinearity, i, unf_nbins_rebinnedS, hAfbLinearityGen_rebinnedS, delta_coefficient, hAfbLinearityDBiasFractional_rebinnedA, hAfbLinearityDBiasFractionalMax_rebinnedA, hAfbLinearityDBias_rebinnedA, hAfbLinearityDBiasMax_rebinnedA, hAfbLinearityDInjectedShape_rebinnedA, hAfbLinearityDInjectedShapeMax_rebinnedA, hAfbLinearityDMeasuredShape_rebinnedA, hAfbLinearityDMeasuredShapeMax_rebinnedA, hAfbLinearitySlope_rebinnedA, hAfbLinearitySlopeMax_rebinnedA, hAfbLinearityDBiasOverStatErr_rebinnedA, hAfbLinearityDBiasOverStatErrMax_rebinnedA, hAfbLinearityCBias_rebinnedA, hAfbLinearityCBiasFractional_rebinnedA, hAfbLinearityCBiasOverStatErr_rebinnedA, hAfbLinearityMaxBiasFrac_rebinnedA, hAfbLinearityMaxBiasFracStat_rebinnedA, hAfbLinearityMaxBiasSlope_rebinnedA, hAfbLinearityMaxBiasFracMax_rebinnedA, hAfbLinearityMaxBiasFracStatMax_rebinnedA, hAfbLinearityMaxBiasSlopeMax_rebinnedA, hAfbLinearityAvgBiasFrac_rebinnedA, hAfbLinearityAvgBiasFracStat_rebinnedA, hAfbLinearityAvgBiasSlope_rebinnedA, hAfbConstantMaxBiasSlope_rebinnedA, hAfbConstantAvgBiasSlope_rebinnedA,numRebinnedCoefficients);
        fillLinearityGraphs(fPseudoUnfAfbResultsRebinnedB, rebinfine, true, "_rebinnedB", nLinearity, i, unf_nbins_rebinnedS, hAfbLinearityGen_rebinnedS, delta_coefficient, hAfbLinearityDBiasFractional_rebinnedB, hAfbLinearityDBiasFractionalMax_rebinnedB, hAfbLinearityDBias_rebinnedB, hAfbLinearityDBiasMax_rebinnedB, hAfbLinearityDInjectedShape_rebinnedB, hAfbLinearityDInjectedShapeMax_rebinnedB, hAfbLinearityDMeasuredShape_rebinnedB, hAfbLinearityDMeasuredShapeMax_rebinnedB, hAfbLinearitySlope_rebinnedB, hAfbLinearitySlopeMax_rebinnedB, hAfbLinearityDBiasOverStatErr_rebinnedB, hAfbLinearityDBiasOverStatErrMax_rebinnedB, hAfbLinearityCBias_rebinnedB, hAfbLinearityCBiasFractional_rebinnedB, hAfbLinearityCBiasOverStatErr_rebinnedB, hAfbLinearityMaxBiasFrac_rebinnedB, hAfbLinearityMaxBiasFracStat_rebinnedB, hAfbLinearityMaxBiasSlope_rebinnedB, hAfbLinearityMaxBiasFracMax_rebinnedB, hAfbLinearityMaxBiasFracStatMax_rebinnedB, hAfbLinearityMaxBiasSlopeMax_rebinnedB, hAfbLinearityAvgBiasFrac_rebinnedB, hAfbLinearityAvgBiasFracStat_rebinnedB, hAfbLinearityAvgBiasSlope_rebinnedB, hAfbConstantMaxBiasSlope_rebinnedB, hAfbConstantAvgBiasSlope_rebinnedB, numRebinnedCoefficients);
      }


      // plot linearity plots
      bool dosimpleplots = kTRUE;

      TCanvas *c1 = new TCanvas("c1", "RegBias", 1280, 1024);
      gStyle->SetPadLeftMargin(0.15);

      c1->Divide(2, 2);
      c1->cd(1);
      hLinearityDInjectedShapeMax->SetLineColor(2);
      if (dosimpleplots) hLinearityDInjectedShapeMax->SetLineWidth(0);
      hLinearityDInjectedShapeMax->Draw();
      hLinearityDMeasuredShapeMax_rebinnedA->SetLineColor(3);
      // if(dosimpleplots) hLinearityDMeasuredShapeMax_rebinnedA->SetLineWidth(0);
      hLinearityDMeasuredShapeMax_rebinnedA->SetLineStyle(2);
      hLinearityDMeasuredShapeMax_rebinnedA->Draw("sames");
      hLinearityDMeasuredShapeMax->SetLineColor(3);
      // if(dosimpleplots) hLinearityDMeasuredShapeMax->SetLineWidth(0);
      hLinearityDMeasuredShapeMax->Draw("sames");
      hLinearityDInjectedShape->Draw("sames");
      hLinearityDInjectedShapeMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hLinearityDInjectedShapeMax_rebinnedA->SetLineWidth(0);
      hLinearityDInjectedShapeMax_rebinnedA->SetLineStyle(2);
      hLinearityDInjectedShapeMax_rebinnedA->Draw("sames");
      hLinearityDInjectedShape_rebinnedA->SetLineStyle(2);
      hLinearityDInjectedShape_rebinnedA->Draw("sames");
      c1->cd(2);
      hLinearitySlopeMax->SetLineColor(2);
      if (dosimpleplots) hLinearitySlopeMax->SetLineWidth(0);
      hLinearitySlopeMax->Draw();
      hLinearitySlope->Draw("sames");
      hLinearitySlopeMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hLinearitySlopeMax_rebinnedA->SetLineWidth(0);
      hLinearitySlopeMax_rebinnedA->SetLineStyle(2);
      hLinearitySlopeMax_rebinnedA->Draw("sames");
      hLinearitySlope_rebinnedA->SetLineStyle(2);
      hLinearitySlope_rebinnedA->Draw("sames");
      // hLinearityDBiasMax->SetLineColor(2);
      // if(dosimpleplots) hLinearityDBiasMax->SetLineWidth(0);
      // hLinearityDBiasMax->Draw();
      // hLinearityDBias->Draw("sames");
      // hLinearityDBiasMax_rebinnedA->SetLineColor(2);
      // if(dosimpleplots) hLinearityDBiasMax_rebinnedA->SetLineWidth(0);
      // hLinearityDBiasMax_rebinnedA->SetLineStyle(2);
      // hLinearityDBiasMax_rebinnedA->Draw("sames");
      // hLinearityDBias_rebinnedA->SetLineStyle(2);
      // hLinearityDBias_rebinnedA->Draw("sames");
      c1->cd(3);
      hLinearityDBiasOverStatErrMax->SetLineColor(2);
      if (dosimpleplots) hLinearityDBiasOverStatErrMax->SetLineWidth(0);
      hLinearityDBiasOverStatErrMax->Draw();
      hLinearityDBiasOverStatErr->Draw("sames");
      hLinearityDBiasOverStatErrMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hLinearityDBiasOverStatErrMax_rebinnedA->SetLineWidth(0);
      hLinearityDBiasOverStatErrMax_rebinnedA->SetLineStyle(2);
      hLinearityDBiasOverStatErrMax_rebinnedA->Draw("sames");
      hLinearityDBiasOverStatErr_rebinnedA->SetLineStyle(2);
      hLinearityDBiasOverStatErr_rebinnedA->Draw("sames");
      c1->cd(4);
      hLinearityDBiasFractionalMax->SetLineColor(2);
      if (dosimpleplots) hLinearityDBiasFractionalMax->SetLineWidth(0);
      hLinearityDBiasFractionalMax->Draw();
      hLinearityDBiasFractional->Draw("sames");
      hLinearityDBiasFractionalMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hLinearityDBiasFractionalMax_rebinnedA->SetLineWidth(0);
      hLinearityDBiasFractionalMax_rebinnedA->SetLineStyle(2);
      hLinearityDBiasFractionalMax_rebinnedA->Draw("sames");
      hLinearityDBiasFractional_rebinnedA->SetLineStyle(2);
      hLinearityDBiasFractional_rebinnedA->Draw("sames");

      c1->Print(Form("%sPlots/Regularisation_bias_%s.pdf", TUNFFILE.Data(), fVariableNames[i].Data()));
      c1->Print(Form("%sPlots/Regularisation_bias_%s.root", TUNFFILE.Data(), fVariableNames[i].Data()));
      c1->Print(Form("%sPlots/Regularisation_bias_%s.C", TUNFFILE.Data(), fVariableNames[i].Data()));
      c1->Close();

      TCanvas *c1r = new TCanvas("c1r", "RegBias", 1280, 1024);
      gStyle->SetPadLeftMargin(0.15);

      c1r->Divide(2, 2);
      c1r->cd(1);
      hLinearityDInjectedShapeMax_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hLinearityDInjectedShapeMax_rebinnedB->SetLineWidth(0);
      if (hLinearityDInjectedShapeMax_rebinnedB->GetMaximum() < hLinearityDInjectedShapeMax_rebinnedA->GetMaximum()) hLinearityDInjectedShapeMax_rebinnedB->SetMaximum(hLinearityDInjectedShapeMax_rebinnedA->GetMaximum());
      if (hLinearityDInjectedShapeMax_rebinnedB->GetMinimum() > hLinearityDInjectedShapeMax_rebinnedA->GetMinimum()) hLinearityDInjectedShapeMax_rebinnedB->SetMinimum(hLinearityDInjectedShapeMax_rebinnedA->GetMinimum());
      hLinearityDInjectedShapeMax_rebinnedB->Draw();
      hLinearityDMeasuredShapeMax_rebinnedA->SetLineColor(3);
      // if(dosimpleplots) hLinearityDMeasuredShapeMax_rebinnedA->SetLineWidth(0);
      hLinearityDMeasuredShapeMax_rebinnedA->SetLineStyle(2);
      hLinearityDMeasuredShapeMax_rebinnedA->Draw("sames");
      hLinearityDMeasuredShapeMax_rebinnedB->SetLineColor(3);
      // if(dosimpleplots) hLinearityDMeasuredShapeMax_rebinnedB->SetLineWidth(0);
      hLinearityDMeasuredShapeMax_rebinnedB->Draw("sames");
      hLinearityDInjectedShape_rebinnedB->Draw("sames");
      hLinearityDInjectedShapeMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hLinearityDInjectedShapeMax_rebinnedA->SetLineWidth(0);
      hLinearityDInjectedShapeMax_rebinnedA->SetLineStyle(2);
      hLinearityDInjectedShapeMax_rebinnedA->Draw("sames");
      hLinearityDInjectedShape_rebinnedA->SetLineStyle(2);
      hLinearityDInjectedShape_rebinnedA->Draw("sames");
      c1r->cd(2);
      hLinearitySlopeMax_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hLinearitySlopeMax_rebinnedB->SetLineWidth(0);
      if (hLinearitySlopeMax_rebinnedB->GetMaximum() < hLinearitySlopeMax_rebinnedA->GetMaximum()) hLinearitySlopeMax_rebinnedB->SetMaximum(hLinearitySlopeMax_rebinnedA->GetMaximum());
      if (hLinearitySlopeMax_rebinnedB->GetMinimum() > hLinearitySlopeMax_rebinnedA->GetMinimum()) hLinearitySlopeMax_rebinnedB->SetMinimum(hLinearitySlopeMax_rebinnedA->GetMinimum());
      if (hLinearitySlopeMax_rebinnedB->GetMaximum() < (1 + 0.15 * (1 - hLinearitySlopeMax_rebinnedB->GetMinimum()))) hLinearitySlopeMax_rebinnedB->SetMaximum(1 + 0.16 * (1 - hLinearitySlopeMax_rebinnedB->GetMinimum()));
      if (hLinearitySlopeMax_rebinnedB->GetMinimum() > (1 + 0.15 * (1 - hLinearitySlopeMax_rebinnedB->GetMaximum()))) hLinearitySlopeMax_rebinnedB->SetMinimum(1 + 0.16 * (1 - hLinearitySlopeMax_rebinnedB->GetMaximum()));
      hLinearitySlopeMax_rebinnedB->Draw();
      hLinearitySlope_rebinnedB->Draw("sames");
      hLinearitySlopeMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hLinearitySlopeMax_rebinnedA->SetLineWidth(0);
      hLinearitySlopeMax_rebinnedA->SetLineStyle(2);
      hLinearitySlopeMax_rebinnedA->Draw("sames");
      hLinearitySlope_rebinnedA->SetLineStyle(2);
      hLinearitySlope_rebinnedA->Draw("sames");

      // hLinearityDBiasMax_rebinnedB->SetLineColor(2);
      // if(dosimpleplots) hLinearityDBiasMax_rebinnedB->SetLineWidth(0);
      // if(hLinearityDBiasMax_rebinnedB->GetMaximum()<hLinearityDBiasMax_rebinnedA->GetMaximum()) hLinearityDBiasMax_rebinnedB->SetMaximum(hLinearityDBiasMax_rebinnedA->GetMaximum());
      // if(hLinearityDBiasMax_rebinnedB->GetMinimum()>hLinearityDBiasMax_rebinnedA->GetMinimum()) hLinearityDBiasMax_rebinnedB->SetMinimum(hLinearityDBiasMax_rebinnedA->GetMinimum());
      // hLinearityDBiasMax_rebinnedB->Draw();
      // hLinearityDBias_rebinnedB->Draw("sames");
      // hLinearityDBiasMax_rebinnedA->SetLineColor(2);
      // if(dosimpleplots) hLinearityDBiasMax_rebinnedA->SetLineWidth(0);
      // hLinearityDBiasMax_rebinnedA->SetLineStyle(2);
      // hLinearityDBiasMax_rebinnedA->Draw("sames");
      // hLinearityDBias_rebinnedA->SetLineStyle(2);
      // hLinearityDBias_rebinnedA->Draw("sames");

      c1r->cd(3);
      hLinearityDBiasOverStatErrMax_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hLinearityDBiasOverStatErrMax_rebinnedB->SetLineWidth(0);
      if (hLinearityDBiasOverStatErrMax_rebinnedB->GetMaximum() < hLinearityDBiasOverStatErrMax_rebinnedA->GetMaximum()) hLinearityDBiasOverStatErrMax_rebinnedB->SetMaximum(hLinearityDBiasOverStatErrMax_rebinnedA->GetMaximum());
      if (hLinearityDBiasOverStatErrMax_rebinnedB->GetMinimum() > hLinearityDBiasOverStatErrMax_rebinnedA->GetMinimum()) hLinearityDBiasOverStatErrMax_rebinnedB->SetMinimum(hLinearityDBiasOverStatErrMax_rebinnedA->GetMinimum());
      hLinearityDBiasOverStatErrMax_rebinnedB->Draw();
      hLinearityDBiasOverStatErr_rebinnedB->Draw("sames");
      hLinearityDBiasOverStatErrMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hLinearityDBiasOverStatErrMax_rebinnedA->SetLineWidth(0);
      hLinearityDBiasOverStatErrMax_rebinnedA->SetLineStyle(2);
      hLinearityDBiasOverStatErrMax_rebinnedA->Draw("sames");
      hLinearityDBiasOverStatErr_rebinnedA->SetLineStyle(2);
      hLinearityDBiasOverStatErr_rebinnedA->Draw("sames");
      c1r->cd(4);
      hLinearityDBiasFractionalMax_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hLinearityDBiasFractionalMax_rebinnedB->SetLineWidth(0);
      if (hLinearityDBiasFractionalMax_rebinnedB->GetMaximum() < hLinearityDBiasFractionalMax_rebinnedA->GetMaximum()) hLinearityDBiasFractionalMax_rebinnedB->SetMaximum(hLinearityDBiasFractionalMax_rebinnedA->GetMaximum());
      if (hLinearityDBiasFractionalMax_rebinnedB->GetMinimum() > hLinearityDBiasFractionalMax_rebinnedA->GetMinimum()) hLinearityDBiasFractionalMax_rebinnedB->SetMinimum(hLinearityDBiasFractionalMax_rebinnedA->GetMinimum());
      hLinearityDBiasFractionalMax_rebinnedB->Draw();
      hLinearityDBiasFractional_rebinnedB->Draw("sames");
      hLinearityDBiasFractionalMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hLinearityDBiasFractionalMax_rebinnedA->SetLineWidth(0);
      hLinearityDBiasFractionalMax_rebinnedA->SetLineStyle(2);
      hLinearityDBiasFractionalMax_rebinnedA->Draw("sames");
      hLinearityDBiasFractional_rebinnedA->SetLineStyle(2);
      hLinearityDBiasFractional_rebinnedA->Draw("sames");

      c1r->Print(Form("%sPlots/Regularisation_bias_rebinned_%s.pdf", TUNFFILE.Data(), fVariableNames[i].Data()));
      c1r->Print(Form("%sPlots/Regularisation_bias_rebinned_%s.root", TUNFFILE.Data(), fVariableNames[i].Data()));
      c1r->Print(Form("%sPlots/Regularisation_bias_rebinned_%s.C", TUNFFILE.Data(), fVariableNames[i].Data()));
      c1r->Close();

      TCanvas *c1c = new TCanvas("c1c", "RegBias", 1280, 1024);
      gStyle->SetPadLeftMargin(0.15);

      c1c->Divide(2, 2);
      c1c->cd(1);
      hLinearityCBias->Draw();
      hLinearityCBias_rebinnedA->SetLineStyle(2);
      hLinearityCBias_rebinnedA->Draw("sames");
      hLinearityCBias_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hLinearityCBias_rebinnedB->SetLineWidth(0);
      hLinearityCBias_rebinnedB->Draw("sames");
      c1c->cd(2);
      if (hLinearityDBias->GetMaximum() < hLinearityCBias->GetMaximum()) hLinearityDBias->SetMaximum(hLinearityCBias->GetMaximum());
      if (hLinearityDBias->GetMinimum() > hLinearityCBias->GetMinimum()) hLinearityDBias->SetMinimum(hLinearityCBias->GetMinimum());
      hLinearityDBias->SetLineColor(kBlue - 1);
      hLinearityDBias->Draw();
      hLinearityDBias_rebinnedA->SetLineStyle(2);
      hLinearityDBias_rebinnedA->SetLineColor(kBlue - 1);
      hLinearityDBias_rebinnedA->Draw("sames");
      hLinearityDBias_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hLinearityDBias_rebinnedB->SetLineWidth(0);
      hLinearityDBias_rebinnedB->SetLineColor(kRed - 1);
      hLinearityDBias_rebinnedB->Draw("sames");
      hLinearityCBias->Draw("sames");
      hLinearityCBias_rebinnedA->SetLineStyle(2);
      hLinearityCBias_rebinnedA->Draw("sames");
      hLinearityCBias_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hLinearityCBias_rebinnedB->SetLineWidth(0);
      hLinearityCBias_rebinnedB->Draw("sames");
      c1c->cd(3);
      hLinearityCBiasOverStatErr->Draw();
      hLinearityCBiasOverStatErr_rebinnedA->SetLineStyle(2);
      hLinearityCBiasOverStatErr_rebinnedA->Draw("sames");
      hLinearityCBiasOverStatErr_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hLinearityCBiasOverStatErr_rebinnedB->SetLineWidth(0);
      hLinearityCBiasOverStatErr_rebinnedB->Draw("sames");
      c1c->cd(4);
      hLinearityCBiasFractional->Draw();
      hLinearityCBiasFractional_rebinnedA->SetLineStyle(2);
      hLinearityCBiasFractional_rebinnedA->Draw("sames");
      hLinearityCBiasFractional_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hLinearityCBiasFractional_rebinnedB->SetLineWidth(0);
      hLinearityCBiasFractional_rebinnedB->Draw("sames");

      c1c->Print(Form("%sPlots/Regularisation_bias_constant_%s.pdf", TUNFFILE.Data(), fVariableNames[i].Data()));
      c1c->Print(Form("%sPlots/Regularisation_bias_constant_%s.root", TUNFFILE.Data(), fVariableNames[i].Data()));
      c1c->Print(Form("%sPlots/Regularisation_bias_constant_%s.C", TUNFFILE.Data(), fVariableNames[i].Data()));
      c1c->Close();

      TCanvas *c1Afbr = new TCanvas("c1Afbr", "RegBias", 1280, 1024);
      gStyle->SetPadLeftMargin(0.15);

      c1Afbr->Divide(2, 2);
      c1Afbr->cd(1);
      hAfbLinearityDInjectedShapeMax_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hAfbLinearityDInjectedShapeMax_rebinnedB->SetLineWidth(0);
      if (hAfbLinearityDInjectedShapeMax_rebinnedB->GetMaximum() < hAfbLinearityDInjectedShapeMax_rebinnedA->GetMaximum()) hAfbLinearityDInjectedShapeMax_rebinnedB->SetMaximum(hAfbLinearityDInjectedShapeMax_rebinnedA->GetMaximum());
      if (hAfbLinearityDInjectedShapeMax_rebinnedB->GetMinimum() > hAfbLinearityDInjectedShapeMax_rebinnedA->GetMinimum()) hAfbLinearityDInjectedShapeMax_rebinnedB->SetMinimum(hAfbLinearityDInjectedShapeMax_rebinnedA->GetMinimum());
      hAfbLinearityDInjectedShapeMax_rebinnedB->Draw();
      hAfbLinearityDMeasuredShapeMax_rebinnedA->SetLineColor(3);
      
      // if(dosimpleplots) hAfbLinearityDMeasuredShapeMax_rebinnedA->SetLineWidth(0);
      hAfbLinearityDMeasuredShapeMax_rebinnedA->SetLineStyle(2);
      hAfbLinearityDMeasuredShapeMax_rebinnedA->Draw("sames");
      hAfbLinearityDMeasuredShapeMax_rebinnedB->SetLineColor(3);
      
      // if(dosimpleplots) hAfbLinearityDMeasuredShapeMax_rebinnedB->SetLineWidth(0);
      hAfbLinearityDMeasuredShapeMax_rebinnedB->Draw("sames");
      hAfbLinearityDInjectedShape_rebinnedB->Draw("sames");
      hAfbLinearityDInjectedShapeMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hAfbLinearityDInjectedShapeMax_rebinnedA->SetLineWidth(0);
      hAfbLinearityDInjectedShapeMax_rebinnedA->SetLineStyle(2);
      hAfbLinearityDInjectedShapeMax_rebinnedA->Draw("sames");
      hAfbLinearityDInjectedShape_rebinnedA->SetLineStyle(2);
      hAfbLinearityDInjectedShape_rebinnedA->Draw("sames");
      c1Afbr->cd(2);
      hAfbLinearitySlopeMax_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hAfbLinearitySlopeMax_rebinnedB->SetLineWidth(0);
      if (hAfbLinearitySlopeMax_rebinnedB->GetMaximum() < hAfbLinearitySlopeMax_rebinnedA->GetMaximum()) hAfbLinearitySlopeMax_rebinnedB->SetMaximum(hAfbLinearitySlopeMax_rebinnedA->GetMaximum());
      if (hAfbLinearitySlopeMax_rebinnedB->GetMinimum() > hAfbLinearitySlopeMax_rebinnedA->GetMinimum()) hAfbLinearitySlopeMax_rebinnedB->SetMinimum(hAfbLinearitySlopeMax_rebinnedA->GetMinimum());
      if (hAfbLinearitySlopeMax_rebinnedB->GetMaximum() < (1 + 0.15 * (1 - hAfbLinearitySlopeMax_rebinnedB->GetMinimum()))) hAfbLinearitySlopeMax_rebinnedB->SetMaximum(1 + 0.16 * (1 - hAfbLinearitySlopeMax_rebinnedB->GetMinimum()));
      if (hAfbLinearitySlopeMax_rebinnedB->GetMinimum() > (1 + 0.15 * (1 - hAfbLinearitySlopeMax_rebinnedB->GetMaximum()))) hAfbLinearitySlopeMax_rebinnedB->SetMinimum(1 + 0.16 * (1 - hAfbLinearitySlopeMax_rebinnedB->GetMaximum()));
      hAfbLinearitySlopeMax_rebinnedB->Draw();
      hAfbLinearitySlope_rebinnedB->Draw("sames");
      hAfbLinearitySlopeMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hAfbLinearitySlopeMax_rebinnedA->SetLineWidth(0);
      hAfbLinearitySlopeMax_rebinnedA->SetLineStyle(2);
      hAfbLinearitySlopeMax_rebinnedA->Draw("sames");
      hAfbLinearitySlope_rebinnedA->SetLineStyle(2);
      hAfbLinearitySlope_rebinnedA->Draw("sames");

      // hAfbLinearityDBiasMax_rebinnedB->SetLineColor(2);
      // if(dosimpleplots) hAfbLinearityDBiasMax_rebinnedB->SetLineWidth(0);
      // if(hAfbLinearityDBiasMax_rebinnedB->GetMaximum()<hAfbLinearityDBiasMax_rebinnedA->GetMaximum()) hAfbLinearityDBiasMax_rebinnedB->SetMaximum(hAfbLinearityDBiasMax_rebinnedA->GetMaximum());
      // if(hAfbLinearityDBiasMax_rebinnedB->GetMinimum()>hAfbLinearityDBiasMax_rebinnedA->GetMinimum()) hAfbLinearityDBiasMax_rebinnedB->SetMinimum(hAfbLinearityDBiasMax_rebinnedA->GetMinimum());
      // hAfbLinearityDBiasMax_rebinnedB->Draw();
      // hAfbLinearityDBias_rebinnedB->Draw("sames");
      // hAfbLinearityDBiasMax_rebinnedA->SetLineColor(2);
      // if(dosimpleplots) hAfbLinearityDBiasMax_rebinnedA->SetLineWidth(0);
      // hAfbLinearityDBiasMax_rebinnedA->SetLineStyle(2);
      // hAfbLinearityDBiasMax_rebinnedA->Draw("sames");
      // hAfbLinearityDBias_rebinnedA->SetLineStyle(2);
      // hAfbLinearityDBias_rebinnedA->Draw("sames");

      c1Afbr->cd(3);
      hAfbLinearityDBiasOverStatErrMax_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hAfbLinearityDBiasOverStatErrMax_rebinnedB->SetLineWidth(0);
      if (hAfbLinearityDBiasOverStatErrMax_rebinnedB->GetMaximum() < hAfbLinearityDBiasOverStatErrMax_rebinnedA->GetMaximum()) hAfbLinearityDBiasOverStatErrMax_rebinnedB->SetMaximum(hAfbLinearityDBiasOverStatErrMax_rebinnedA->GetMaximum());
      if (hAfbLinearityDBiasOverStatErrMax_rebinnedB->GetMinimum() > hAfbLinearityDBiasOverStatErrMax_rebinnedA->GetMinimum()) hAfbLinearityDBiasOverStatErrMax_rebinnedB->SetMinimum(hAfbLinearityDBiasOverStatErrMax_rebinnedA->GetMinimum());
      hAfbLinearityDBiasOverStatErrMax_rebinnedB->Draw();
      hAfbLinearityDBiasOverStatErr_rebinnedB->Draw("sames");
      hAfbLinearityDBiasOverStatErrMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hAfbLinearityDBiasOverStatErrMax_rebinnedA->SetLineWidth(0);
      hAfbLinearityDBiasOverStatErrMax_rebinnedA->SetLineStyle(2);
      hAfbLinearityDBiasOverStatErrMax_rebinnedA->Draw("sames");
      hAfbLinearityDBiasOverStatErr_rebinnedA->SetLineStyle(2);
      hAfbLinearityDBiasOverStatErr_rebinnedA->Draw("sames");
      c1Afbr->cd(4);
      hAfbLinearityDBiasFractionalMax_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hAfbLinearityDBiasFractionalMax_rebinnedB->SetLineWidth(0);
      if (hAfbLinearityDBiasFractionalMax_rebinnedB->GetMaximum() < hAfbLinearityDBiasFractionalMax_rebinnedA->GetMaximum()) hAfbLinearityDBiasFractionalMax_rebinnedB->SetMaximum(hAfbLinearityDBiasFractionalMax_rebinnedA->GetMaximum());
      if (hAfbLinearityDBiasFractionalMax_rebinnedB->GetMinimum() > hAfbLinearityDBiasFractionalMax_rebinnedA->GetMinimum()) hAfbLinearityDBiasFractionalMax_rebinnedB->SetMinimum(hAfbLinearityDBiasFractionalMax_rebinnedA->GetMinimum());
      hAfbLinearityDBiasFractionalMax_rebinnedB->Draw();
      hAfbLinearityDBiasFractional_rebinnedB->Draw("sames");
      hAfbLinearityDBiasFractionalMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hAfbLinearityDBiasFractionalMax_rebinnedA->SetLineWidth(0);
      hAfbLinearityDBiasFractionalMax_rebinnedA->SetLineStyle(2);
      hAfbLinearityDBiasFractionalMax_rebinnedA->Draw("sames");
      hAfbLinearityDBiasFractional_rebinnedA->SetLineStyle(2);
      hAfbLinearityDBiasFractional_rebinnedA->Draw("sames");

      c1Afbr->Print(Form("%sPlots/Regularisation_bias_Afb_rebinned_%s.pdf", TUNFFILE.Data(), fVariableNames[i].Data()));
      c1Afbr->Print(Form("%sPlots/Regularisation_bias_Afb_rebinned_%s.root", TUNFFILE.Data(), fVariableNames[i].Data()));
      c1Afbr->Print(Form("%sPlots/Regularisation_bias_Afb_rebinned_%s.C", TUNFFILE.Data(), fVariableNames[i].Data()));
      c1Afbr->Close();

      TCanvas *c1Afbc = new TCanvas("c1Afbc", "RegBias", 1280, 1024);
      gStyle->SetPadLeftMargin(0.15);

      c1Afbc->Divide(2, 2);
      c1Afbc->cd(1);
      hAfbLinearityCBias_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hAfbLinearityCBias_rebinnedB->SetLineWidth(0);
      hAfbLinearityCBias_rebinnedB->Draw();
      hAfbLinearityCBias_rebinnedA->SetLineStyle(2);
      hAfbLinearityCBias_rebinnedA->Draw("sames");
      c1Afbc->cd(2);
      if (hAfbLinearityDBias_rebinnedB->GetMaximum() < hAfbLinearityCBias_rebinnedB->GetMaximum()) hAfbLinearityDBias_rebinnedB->SetMaximum(hAfbLinearityCBias_rebinnedB->GetMaximum());
      if (hAfbLinearityDBias_rebinnedB->GetMinimum() > hAfbLinearityCBias_rebinnedB->GetMinimum()) hAfbLinearityDBias_rebinnedB->SetMinimum(hAfbLinearityCBias_rebinnedB->GetMinimum());
      hAfbLinearityDBias_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hAfbLinearityDBias_rebinnedB->SetLineWidth(0);
      hAfbLinearityDBias_rebinnedB->SetLineColor(kRed - 1);
      hAfbLinearityDBias_rebinnedB->Draw();
      hAfbLinearityDBias_rebinnedA->SetLineStyle(2);
      hAfbLinearityDBias_rebinnedA->SetLineColor(kBlue - 1);
      hAfbLinearityDBias_rebinnedA->Draw("sames");
      hAfbLinearityCBias_rebinnedA->SetLineStyle(2);
      hAfbLinearityCBias_rebinnedA->Draw("sames");
      hAfbLinearityCBias_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hAfbLinearityCBias_rebinnedB->SetLineWidth(0);
      hAfbLinearityCBias_rebinnedB->Draw("sames");
      c1Afbc->cd(3);
      hAfbLinearityCBiasOverStatErr_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hAfbLinearityCBiasOverStatErr_rebinnedB->SetLineWidth(0);
      hAfbLinearityCBiasOverStatErr_rebinnedB->Draw();
      hAfbLinearityCBiasOverStatErr_rebinnedA->SetLineStyle(2);
      hAfbLinearityCBiasOverStatErr_rebinnedA->Draw("sames");
      c1Afbc->cd(4);
      hAfbLinearityCBiasFractional_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hAfbLinearityCBiasFractional_rebinnedB->SetLineWidth(0);
      hAfbLinearityCBiasFractional_rebinnedB->Draw();
      hAfbLinearityCBiasFractional_rebinnedA->SetLineStyle(2);
      hAfbLinearityCBiasFractional_rebinnedA->Draw("sames");

      c1Afbc->Print(Form("%sPlots/Regularisation_bias_Afb_constant_%s.pdf", TUNFFILE.Data(), fVariableNames[i].Data()));
      c1Afbc->Print(Form("%sPlots/Regularisation_bias_Afb_constant_%s.root", TUNFFILE.Data(), fVariableNames[i].Data()));
      c1Afbc->Print(Form("%sPlots/Regularisation_bias_Afb_constant_%s.C", TUNFFILE.Data(), fVariableNames[i].Data()));
      c1Afbc->Close();

      if (i == fVariableNames.size() - 1) gStyle->SetPadBottomMargin(0.12);
      std::vector<TString> CoefNames2;
      for (unsigned int iVar = 0; iVar < fVariableNames.size(); ++iVar) {
          for (unsigned int iCoef = 0; iCoef < numRebinnedCoefficients; iCoef++) {
              CoefNames2.push_back(MatrixUnf::CoefName(fVariableNames[iVar], kFALSE));
          }
      }

      TCanvas *c2 = new TCanvas("c2", "RegBias", 1024, 1280);
      gStyle->SetPadLeftMargin(0.08);

      c2->Divide(1, 3);
      c2->cd(2);
      gPad->SetLogy();
      hLinearityMaxBiasFracMax->SetLineColor(2);
      if (dosimpleplots) hLinearityMaxBiasFracMax->SetLineWidth(0);
      hLinearityMaxBiasFracMax->GetXaxis()->SetTitle("Variable");
      // hLinearityMaxBiasFracMax->GetYaxis()->SetTitle("Max/Avg #DeltaFractionalBias / #DeltaCoef.");
      hLinearityMaxBiasFracMax->GetYaxis()->SetTitle("Max and avg |fractional bias| per unit injected coef.");
      hLinearityMaxBiasFracMax->Draw();

      hLinearityMaxBiasFracMax->SetTitle("");
      c2->Update();
      hLinearityMaxBiasFracMax->GetXaxis()->SetTitle("");
      hLinearityMaxBiasFracMax->GetXaxis()->SetTitleOffset(2.0);
      hLinearityMaxBiasFracMax->GetXaxis()->SetTitleSize(0.04);
      hLinearityMaxBiasFracMax->GetXaxis()->SetLabelSize(0.03);
      // hLinearityMaxBiasFracMax->GetYaxis()->SetTitleOffset(1.8);
      hLinearityMaxBiasFracMax->GetYaxis()->SetTitleSize(0.04);
      hLinearityMaxBiasFracMax->GetYaxis()->SetLabelSize(0.04);
      for (unsigned int i_bin = 1; i_bin <= CoefNames2.size(); i_bin++) {
          hLinearityMaxBiasFracMax->GetXaxis()->SetBinLabel(i_bin, CoefNames2.at(i_bin - 1).Data());
      }
      hLinearityMaxBiasFracMax->GetXaxis()->LabelsOption("v");
      // hLinearityMaxBiasFracMax->GetYaxis()->LabelsOption("v");
      // DrawCMSLabels(0.04);
      if (CHAN != "combined") DrawDecayChLabel(CHAN, 0.04);

      hLinearityMaxBiasFrac->Draw("sames");
      hLinearityAvgBiasFrac->SetLineColor(3);
      // if(dosimpleplots) hLinearityAvgBiasFrac->SetLineWidth(0);
      hLinearityAvgBiasFrac->Draw("sames");
      hLinearityMaxBiasFracMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hLinearityMaxBiasFracMax_rebinnedA->SetLineWidth(0);
      hLinearityMaxBiasFracMax_rebinnedA->SetLineStyle(2);
      hLinearityMaxBiasFracMax_rebinnedA->Draw("sames");
      hLinearityMaxBiasFrac_rebinnedA->SetLineStyle(2);
      hLinearityMaxBiasFrac_rebinnedA->Draw("sames");
      hLinearityAvgBiasFrac_rebinnedA->SetLineColor(3);
      // if(dosimpleplots) hLinearityAvgBiasFrac_rebinnedA->SetLineWidth(0);
      hLinearityAvgBiasFrac_rebinnedA->SetLineStyle(2);
      hLinearityAvgBiasFrac_rebinnedA->Draw("sames");
      c2->cd(3);
      gPad->SetLogy();
      hLinearityMaxBiasFracStatMax->SetLineColor(2);
      if (dosimpleplots) hLinearityMaxBiasFracStatMax->SetLineWidth(0);
      hLinearityMaxBiasFracStatMax->GetXaxis()->SetTitle("Variable");
      // hLinearityMaxBiasFracStatMax->GetYaxis()->SetTitle("Max/Avg #Delta(Bias/StatUnc.) / #DeltaCoef.");
      hLinearityMaxBiasFracStatMax->GetYaxis()->SetTitle("Max and average |bias/#sigma_{stat}| per unit injected coef.");
      hLinearityMaxBiasFracStatMax->Draw();

      hLinearityMaxBiasFracStatMax->SetTitle("");
      c2->Update();
      hLinearityMaxBiasFracStatMax->GetXaxis()->SetTitle("");
      hLinearityMaxBiasFracStatMax->GetXaxis()->SetTitleOffset(2.0);
      hLinearityMaxBiasFracStatMax->GetXaxis()->SetTitleSize(0.04);
      hLinearityMaxBiasFracStatMax->GetXaxis()->SetLabelSize(0.03);
      // hLinearityMaxBiasFracStatMax->GetYaxis()->SetTitleOffset(1.8);
      hLinearityMaxBiasFracStatMax->GetYaxis()->SetTitleSize(0.04);
      hLinearityMaxBiasFracStatMax->GetYaxis()->SetLabelSize(0.04);
      for (unsigned int i_bin = 1; i_bin <= CoefNames2.size(); i_bin++) {
          hLinearityMaxBiasFracStatMax->GetXaxis()->SetBinLabel(i_bin, CoefNames2.at(i_bin - 1).Data());
      }
      hLinearityMaxBiasFracStatMax->GetXaxis()->LabelsOption("v");
      // hLinearityMaxBiasFracStatMax->GetYaxis()->LabelsOption("v");
      // DrawCMSLabels(0.04);
      if (CHAN != "combined") DrawDecayChLabel(CHAN, 0.04);

      hLinearityMaxBiasFracStat->Draw("sames");
      hLinearityAvgBiasFracStat->SetLineColor(3);
      // if(dosimpleplots) hLinearityAvgBiasFracStat->SetLineWidth(0);
      hLinearityAvgBiasFracStat->Draw("sames");
      hLinearityMaxBiasFracStatMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hLinearityMaxBiasFracStatMax_rebinnedA->SetLineWidth(0);
      hLinearityMaxBiasFracStatMax_rebinnedA->SetLineStyle(2);
      hLinearityMaxBiasFracStatMax_rebinnedA->Draw("sames");
      hLinearityMaxBiasFracStat_rebinnedA->SetLineStyle(2);
      hLinearityMaxBiasFracStat_rebinnedA->Draw("sames");
      hLinearityAvgBiasFracStat_rebinnedA->SetLineColor(3);
      // if(dosimpleplots) hLinearityAvgBiasFracStat_rebinnedA->SetLineWidth(0);
      hLinearityAvgBiasFracStat_rebinnedA->SetLineStyle(2);
      hLinearityAvgBiasFracStat_rebinnedA->Draw("sames");
      c2->cd(1);
      gPad->SetLogy();
      hLinearityMaxBiasSlopeMax->SetLineColor(2);
      if (dosimpleplots) hLinearityMaxBiasSlopeMax->SetLineWidth(0);
      hLinearityMaxBiasSlopeMax->GetXaxis()->SetTitle("Variable");
      hLinearityMaxBiasSlopeMax->GetYaxis()->SetTitle("Maximum and average |slope-1| of linearity plot");
      hLinearityMaxBiasSlopeMax->Draw();

      hLinearityMaxBiasSlopeMax->SetTitle("");
      c2->Update();
      hLinearityMaxBiasSlopeMax->GetXaxis()->SetTitle("");
      hLinearityMaxBiasSlopeMax->GetXaxis()->SetTitleOffset(2.0);
      hLinearityMaxBiasSlopeMax->GetXaxis()->SetTitleSize(0.04);
      hLinearityMaxBiasSlopeMax->GetXaxis()->SetLabelSize(0.03);
      // hLinearityMaxBiasSlopeMax->GetYaxis()->SetTitleOffset(1.8);
      hLinearityMaxBiasSlopeMax->GetYaxis()->SetTitleSize(0.04);
      hLinearityMaxBiasSlopeMax->GetYaxis()->SetLabelSize(0.04);
      for (unsigned int i_bin = 1; i_bin <= CoefNames2.size(); i_bin++) {
          hLinearityMaxBiasSlopeMax->GetXaxis()->SetBinLabel(i_bin, CoefNames2.at(i_bin - 1).Data());
      }
      hLinearityMaxBiasSlopeMax->GetXaxis()->LabelsOption("v");
      // hLinearityMaxBiasSlopeMax->GetYaxis()->LabelsOption("v");
      // DrawCMSLabels(0.04);
      if (CHAN != "combined") DrawDecayChLabel(CHAN, 0.04);

      hLinearityMaxBiasSlope->Draw("sames");
      hLinearityAvgBiasSlope->SetLineColor(3);
      // if(dosimpleplots) hLinearityAvgBiasSlope->SetLineWidth(0);
      hLinearityAvgBiasSlope->Draw("sames");
      hLinearityMaxBiasSlopeMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hLinearityMaxBiasSlopeMax_rebinnedA->SetLineWidth(0);
      hLinearityMaxBiasSlopeMax_rebinnedA->SetLineStyle(2);
      hLinearityMaxBiasSlopeMax_rebinnedA->Draw("sames");
      hLinearityMaxBiasSlope_rebinnedA->SetLineStyle(2);
      hLinearityMaxBiasSlope_rebinnedA->Draw("sames");
      hLinearityAvgBiasSlope_rebinnedA->SetLineColor(3);
      // if(dosimpleplots) hLinearityAvgBiasSlope_rebinnedA->SetLineWidth(0);
      hLinearityAvgBiasSlope_rebinnedA->SetLineStyle(2);
      hLinearityAvgBiasSlope_rebinnedA->Draw("sames");

      c2->Print(Form("%sPlots/Regularisation_bias_summary.pdf", TUNFFILE.Data()));
      c2->Print(Form("%sPlots/Regularisation_bias_summary.root", TUNFFILE.Data()));
      c2->Print(Form("%sPlots/Regularisation_bias_summary.C", TUNFFILE.Data()));
      c2->Close();

      TCanvas *c2r = new TCanvas("c2r", "RegBias", 1024, 1280);
      gStyle->SetPadLeftMargin(0.08);

      c2r->Divide(1, 3);
      c2r->cd(2);
      gPad->SetLogy();
      hLinearityMaxBiasFracMax_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hLinearityMaxBiasFracMax_rebinnedB->SetLineWidth(0);
      hLinearityMaxBiasFracMax_rebinnedB->GetXaxis()->SetTitle("Variable");
      // hLinearityMaxBiasFracMax_rebinnedB->GetYaxis()->SetTitle("Max/Avg #DeltaFractionalBias / #DeltaCoef.");
      hLinearityMaxBiasFracMax_rebinnedB->GetYaxis()->SetTitle("Max and avg |fractional bias| per unit injected coef.");
      if (hLinearityMaxBiasFracMax_rebinnedB->GetMaximum() < hLinearityMaxBiasFracMax_rebinnedA->GetMaximum()) hLinearityMaxBiasFracMax_rebinnedB->SetMaximum(hLinearityMaxBiasFracMax_rebinnedA->GetMaximum());

      //if (hLinearityMaxBiasFracMax_rebinnedB->GetMinimum() > hLinearityMaxBiasFracMax_rebinnedA->GetMinimum()) hLinearityMaxBiasFracMax_rebinnedB->SetMinimum(hLinearityMaxBiasFracMax_rebinnedA->GetMinimum());
      if (hLinearityMaxBiasFracMax_rebinnedB->GetMinimum() > hLinearityAvgBiasFrac_rebinnedA->GetMinimum()) hLinearityMaxBiasFracMax_rebinnedB->SetMinimum(hLinearityAvgBiasFrac_rebinnedA->GetMinimum());
      //if (hLinearityMaxBiasFracMax_rebinnedB->GetMinimum() > hLinearityMaxBiasFrac_rebinnedA->GetMinimum()) hLinearityMaxBiasFracMax_rebinnedB->SetMinimum(hLinearityMaxBiasFrac_rebinnedA->GetMinimum());

      hLinearityMaxBiasFracMax_rebinnedB->Draw();

      hLinearityMaxBiasFracMax_rebinnedB->SetTitle("");
      c2r->Update();
      hLinearityMaxBiasFracMax_rebinnedB->GetXaxis()->SetTitle("");
      hLinearityMaxBiasFracMax_rebinnedB->GetXaxis()->SetTitleOffset(2.0);
      hLinearityMaxBiasFracMax_rebinnedB->GetXaxis()->SetTitleSize(0.04);
      hLinearityMaxBiasFracMax_rebinnedB->GetXaxis()->SetLabelSize(0.03);
      // hLinearityMaxBiasFracMax_rebinnedB->GetYaxis()->SetTitleOffset(1.8);
      hLinearityMaxBiasFracMax_rebinnedB->GetYaxis()->SetTitleSize(0.04);
      hLinearityMaxBiasFracMax_rebinnedB->GetYaxis()->SetLabelSize(0.04);
      for (unsigned int i_bin = 1; i_bin <= CoefNames2.size(); i_bin++) {
          hLinearityMaxBiasFracMax_rebinnedB->GetXaxis()->SetBinLabel(i_bin, CoefNames2.at(i_bin - 1).Data());
      }
      hLinearityMaxBiasFracMax_rebinnedB->GetXaxis()->LabelsOption("v");
      // hLinearityMaxBiasFracMax_rebinnedB->GetYaxis()->LabelsOption("v");
      // DrawCMSLabels(0.04);
      if (CHAN != "combined") DrawDecayChLabel(CHAN, 0.04);

      hLinearityMaxBiasFrac_rebinnedB->Draw("sames");
      hLinearityAvgBiasFrac_rebinnedB->SetLineColor(3);
      // if(dosimpleplots) hLinearityAvgBiasFrac_rebinnedB->SetLineWidth(0);
      hLinearityAvgBiasFrac_rebinnedB->Draw("sames");
      hLinearityMaxBiasFracMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hLinearityMaxBiasFracMax_rebinnedA->SetLineWidth(0);
      hLinearityMaxBiasFracMax_rebinnedA->SetLineStyle(2);
      hLinearityMaxBiasFracMax_rebinnedA->Draw("sames");
      hLinearityMaxBiasFrac_rebinnedA->SetLineStyle(2);
      hLinearityMaxBiasFrac_rebinnedA->Draw("sames");
      hLinearityAvgBiasFrac_rebinnedA->SetLineColor(3);
      // if(dosimpleplots) hLinearityAvgBiasFrac_rebinnedA->SetLineWidth(0);
      hLinearityAvgBiasFrac_rebinnedA->SetLineStyle(2);
      hLinearityAvgBiasFrac_rebinnedA->Draw("sames");
      c2r->cd(3);
      gPad->SetLogy();
      hLinearityMaxBiasFracStatMax_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hLinearityMaxBiasFracStatMax_rebinnedB->SetLineWidth(0);
      hLinearityMaxBiasFracStatMax_rebinnedB->GetXaxis()->SetTitle("Variable");
      // hLinearityMaxBiasFracStatMax_rebinnedB->GetYaxis()->SetTitle("Max/Avg #Delta(Bias/StatUnc.) / #DeltaCoef.");
      hLinearityMaxBiasFracStatMax_rebinnedB->GetYaxis()->SetTitle("Max and average |bias/#sigma_{stat}| per unit injected coef.");
      if (hLinearityMaxBiasFracStatMax_rebinnedB->GetMaximum() < hLinearityMaxBiasFracStatMax_rebinnedA->GetMaximum()) hLinearityMaxBiasFracStatMax_rebinnedB->SetMaximum(hLinearityMaxBiasFracStatMax_rebinnedA->GetMaximum());
      if (hLinearityMaxBiasFracStatMax_rebinnedB->GetMinimum() > hLinearityAvgBiasFracStat_rebinnedA->GetMinimum()) hLinearityMaxBiasFracStatMax_rebinnedB->SetMinimum(hLinearityAvgBiasFracStat_rebinnedA->GetMinimum());
      hLinearityMaxBiasFracStatMax_rebinnedB->Draw();

      hLinearityMaxBiasFracStatMax_rebinnedB->SetTitle("");
      c2r->Update();
      hLinearityMaxBiasFracStatMax_rebinnedB->GetXaxis()->SetTitle("");
      hLinearityMaxBiasFracStatMax_rebinnedB->GetXaxis()->SetTitleOffset(2.0);
      hLinearityMaxBiasFracStatMax_rebinnedB->GetXaxis()->SetTitleSize(0.04);
      hLinearityMaxBiasFracStatMax_rebinnedB->GetXaxis()->SetLabelSize(0.03);
      // hLinearityMaxBiasFracStatMax_rebinnedB->GetYaxis()->SetTitleOffset(1.8);
      hLinearityMaxBiasFracStatMax_rebinnedB->GetYaxis()->SetTitleSize(0.04);
      hLinearityMaxBiasFracStatMax_rebinnedB->GetYaxis()->SetLabelSize(0.04);
      for (unsigned int i_bin = 1; i_bin <= CoefNames2.size(); i_bin++) {
          hLinearityMaxBiasFracStatMax_rebinnedB->GetXaxis()->SetBinLabel(i_bin, CoefNames2.at(i_bin - 1).Data());
      }
      hLinearityMaxBiasFracStatMax_rebinnedB->GetXaxis()->LabelsOption("v");
      // hLinearityMaxBiasFracStatMax_rebinnedB->GetYaxis()->LabelsOption("v");
      // DrawCMSLabels(0.04);
      if (CHAN != "combined") DrawDecayChLabel(CHAN, 0.04);

      hLinearityMaxBiasFracStat_rebinnedB->Draw("sames");
      hLinearityAvgBiasFracStat_rebinnedB->SetLineColor(3);
      // if(dosimpleplots) hLinearityAvgBiasFracStat_rebinnedB->SetLineWidth(0);
      hLinearityAvgBiasFracStat_rebinnedB->Draw("sames");
      hLinearityMaxBiasFracStatMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hLinearityMaxBiasFracStatMax_rebinnedA->SetLineWidth(0);
      hLinearityMaxBiasFracStatMax_rebinnedA->SetLineStyle(2);
      hLinearityMaxBiasFracStatMax_rebinnedA->Draw("sames");
      hLinearityMaxBiasFracStat_rebinnedA->SetLineStyle(2);
      hLinearityMaxBiasFracStat_rebinnedA->Draw("sames");
      hLinearityAvgBiasFracStat_rebinnedA->SetLineColor(3);
      // if(dosimpleplots) hLinearityAvgBiasFracStat_rebinnedA->SetLineWidth(0);
      hLinearityAvgBiasFracStat_rebinnedA->SetLineStyle(2);
      hLinearityAvgBiasFracStat_rebinnedA->Draw("sames");
      c2r->cd(1);
      gPad->SetLogy();
      hLinearityMaxBiasSlopeMax_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hLinearityMaxBiasSlopeMax_rebinnedB->SetLineWidth(0);
      hLinearityMaxBiasSlopeMax_rebinnedB->GetXaxis()->SetTitle("Variable");
      hLinearityMaxBiasSlopeMax_rebinnedB->GetYaxis()->SetTitle("Maximum and average |slope-1| of linearity plot");
      if (hLinearityMaxBiasSlopeMax_rebinnedB->GetMaximum() < hLinearityMaxBiasSlopeMax_rebinnedA->GetMaximum()) hLinearityMaxBiasSlopeMax_rebinnedB->SetMaximum(hLinearityMaxBiasSlopeMax_rebinnedA->GetMaximum());
      if (hLinearityMaxBiasSlopeMax_rebinnedB->GetMinimum() > hLinearityAvgBiasSlope_rebinnedA->GetMinimum()) hLinearityMaxBiasSlopeMax_rebinnedB->SetMinimum(hLinearityAvgBiasSlope_rebinnedA->GetMinimum());
      hLinearityMaxBiasSlopeMax_rebinnedB->Draw();

      hLinearityMaxBiasSlopeMax_rebinnedB->SetTitle("");
      c2r->Update();
      hLinearityMaxBiasSlopeMax_rebinnedB->GetXaxis()->SetTitle("");
      hLinearityMaxBiasSlopeMax_rebinnedB->GetXaxis()->SetTitleOffset(2.0);
      hLinearityMaxBiasSlopeMax_rebinnedB->GetXaxis()->SetTitleSize(0.04);
      hLinearityMaxBiasSlopeMax_rebinnedB->GetXaxis()->SetLabelSize(0.03);
      // hLinearityMaxBiasSlopeMax_rebinnedB->GetYaxis()->SetTitleOffset(1.8);
      hLinearityMaxBiasSlopeMax_rebinnedB->GetYaxis()->SetTitleSize(0.04);
      hLinearityMaxBiasSlopeMax_rebinnedB->GetYaxis()->SetLabelSize(0.04);
      for (unsigned int i_bin = 1; i_bin <= CoefNames2.size(); i_bin++) {
          hLinearityMaxBiasSlopeMax_rebinnedB->GetXaxis()->SetBinLabel(i_bin, CoefNames2.at(i_bin - 1).Data());
      }
      hLinearityMaxBiasSlopeMax_rebinnedB->GetXaxis()->LabelsOption("v");
      // hLinearityMaxBiasSlopeMax_rebinnedB->GetYaxis()->LabelsOption("v");
      // DrawCMSLabels(0.04);
      if (CHAN != "combined") DrawDecayChLabel(CHAN, 0.04);

      hLinearityMaxBiasSlope_rebinnedB->Draw("sames");
      hLinearityAvgBiasSlope_rebinnedB->SetLineColor(3);
      // if(dosimpleplots) hLinearityAvgBiasSlope_rebinnedB->SetLineWidth(0);
      hLinearityAvgBiasSlope_rebinnedB->Draw("sames");
      hLinearityMaxBiasSlopeMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hLinearityMaxBiasSlopeMax_rebinnedA->SetLineWidth(0);
      hLinearityMaxBiasSlopeMax_rebinnedA->SetLineStyle(2);
      hLinearityMaxBiasSlopeMax_rebinnedA->Draw("sames");
      hLinearityMaxBiasSlope_rebinnedA->SetLineStyle(2);
      hLinearityMaxBiasSlope_rebinnedA->Draw("sames");
      hLinearityAvgBiasSlope_rebinnedA->SetLineColor(3);
      // if(dosimpleplots) hLinearityAvgBiasSlope_rebinnedA->SetLineWidth(0);
      hLinearityAvgBiasSlope_rebinnedA->SetLineStyle(2);
      hLinearityAvgBiasSlope_rebinnedA->Draw("sames");

      c2r->Print(Form("%sPlots/Regularisation_bias_rebinned_summary.pdf", TUNFFILE.Data()));
      c2r->Print(Form("%sPlots/Regularisation_bias_rebinned_summary.root", TUNFFILE.Data()));
      c2r->Print(Form("%sPlots/Regularisation_bias_rebinned_summary.C", TUNFFILE.Data()));
      c2r->Close();

      TCanvas *c2Afbr = new TCanvas("c2Afbr", "RegBias", 1024, 900);
      // TCanvas *c2Afbr = new TCanvas("c2Afbr", "RegBias",1024,1280);
      gStyle->SetPadLeftMargin(0.08);

      c2Afbr->Divide(1, 2);
      // c2Afbr->Divide(1,3);
      c2Afbr->cd(2);
      // gPad->SetLogy(0);
      // hAfbLinearityMaxBiasFracMax_rebinnedB->SetLineColor(2);
      // if(dosimpleplots) hAfbLinearityMaxBiasFracMax_rebinnedB->SetLineWidth(0);
      // hAfbLinearityMaxBiasFracMax_rebinnedB->GetXaxis()->SetTitle("Variable");
      // hAfbLinearityMaxBiasFracMax_rebinnedB->GetYaxis()->SetTitle("#DeltaBias / #DeltaCoef.");
      // hAfbLinearityMaxBiasFracMax_rebinnedB->GetYaxis()->SetTitle("Measured coef. bias per unit injected coef.");
      // if(hAfbLinearityMaxBiasFracMax_rebinnedB->GetMaximum()<hAfbLinearityMaxBiasFracMax_rebinnedA->GetMaximum()) hAfbLinearityMaxBiasFracMax_rebinnedB->SetMaximum(hAfbLinearityMaxBiasFracMax_rebinnedA->GetMaximum());
      // if(hAfbLinearityMaxBiasFracMax_rebinnedB->GetMinimum()>hAfbLinearityMaxBiasFracMax_rebinnedA->GetMinimum()) hAfbLinearityMaxBiasFracMax_rebinnedB->SetMinimum(hAfbLinearityMaxBiasFracMax_rebinnedA->GetMinimum());
      // hAfbLinearityMaxBiasFracMax_rebinnedB->Draw();
      // hAfbLinearityMaxBiasFrac_rebinnedB->Draw("sames");
      // hAfbLinearityAvgBiasFrac_rebinnedB->SetLineColor(3);
      // if(dosimpleplots) hAfbLinearityAvgBiasFrac_rebinnedB->SetLineWidth(0);
      // hAfbLinearityAvgBiasFrac_rebinnedB->Draw("sames");
      // hAfbLinearityMaxBiasFracMax_rebinnedA->SetLineColor(2);
      // if(dosimpleplots) hAfbLinearityMaxBiasFracMax_rebinnedA->SetLineWidth(0);
      // hAfbLinearityMaxBiasFracMax_rebinnedA->SetLineStyle(2);
      // hAfbLinearityMaxBiasFracMax_rebinnedA->Draw("sames");
      // hAfbLinearityMaxBiasFrac_rebinnedA->SetLineStyle(2);
      // hAfbLinearityMaxBiasFrac_rebinnedA->Draw("sames");
      // hAfbLinearityAvgBiasFrac_rebinnedA->SetLineColor(3);
      // if(dosimpleplots) hAfbLinearityAvgBiasFrac_rebinnedA->SetLineWidth(0);
      // hAfbLinearityAvgBiasFrac_rebinnedA->SetLineStyle(2);
      // hAfbLinearityAvgBiasFrac_rebinnedA->Draw("sames");
      // c2Afbr->cd(3);
      gPad->SetLogy(0);
      hAfbLinearityMaxBiasFracStatMax_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hAfbLinearityMaxBiasFracStatMax_rebinnedB->SetLineWidth(0);
      hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetXaxis()->SetTitle("Variable");
      // hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetYaxis()->SetTitle("#Delta(Bias/StatUnc.) / #DeltaCoef.");
      hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetYaxis()->SetTitle("Measured coef. bias/#sigma_{stat} per unit injected coef.");
      if (hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetMaximum() < hAfbLinearityMaxBiasFracStatMax_rebinnedA->GetMaximum()) hAfbLinearityMaxBiasFracStatMax_rebinnedB->SetMaximum(hAfbLinearityMaxBiasFracStatMax_rebinnedA->GetMaximum());
      if (hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetMinimum() > hAfbLinearityMaxBiasFracStat_rebinnedA->GetMinimum()) hAfbLinearityMaxBiasFracStatMax_rebinnedB->SetMinimum(hAfbLinearityMaxBiasFracStat_rebinnedA->GetMinimum());
      // if(hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetMaximum()>0 && hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetMinimum()>-hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetMaximum()/8. ) hAfbLinearityMaxBiasFracStatMax_rebinnedB->SetMinimum(-hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetMaximum()/8.);
      hAfbLinearityMaxBiasFracStatMax_rebinnedB->Draw();

      hAfbLinearityMaxBiasFracStatMax_rebinnedB->SetTitle("");
      c2Afbr->Update();
      hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetXaxis()->SetTitle("");
      hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetXaxis()->SetTitleOffset(2.0);
      hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetXaxis()->SetTitleSize(0.04);
      hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetXaxis()->SetLabelSize(0.03);
      // hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetYaxis()->SetTitleOffset(1.8);
      hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetYaxis()->SetTitleSize(0.04);
      hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetYaxis()->SetLabelSize(0.04);
      for (unsigned int i_bin = 1; i_bin <= CoefNames2.size(); i_bin++) {
          hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetXaxis()->SetBinLabel(i_bin, CoefNames2.at(i_bin - 1).Data());
      }
      hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetXaxis()->LabelsOption("v");
      // hAfbLinearityMaxBiasFracStatMax_rebinnedB->GetYaxis()->LabelsOption("v");
      // DrawCMSLabels(0.04);
      if (CHAN != "combined") DrawDecayChLabel(CHAN, 0.04);

      hAfbLinearityMaxBiasFracStat_rebinnedB->Draw("sames");
      hAfbLinearityAvgBiasFracStat_rebinnedB->SetLineColor(3);
      if (dosimpleplots) hAfbLinearityAvgBiasFracStat_rebinnedB->SetLineWidth(0);
      hAfbLinearityAvgBiasFracStat_rebinnedB->Draw("sames");
      hAfbLinearityMaxBiasFracStatMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hAfbLinearityMaxBiasFracStatMax_rebinnedA->SetLineWidth(0);
      hAfbLinearityMaxBiasFracStatMax_rebinnedA->SetLineStyle(2);
      hAfbLinearityMaxBiasFracStatMax_rebinnedA->Draw("sames");
      hAfbLinearityMaxBiasFracStat_rebinnedA->SetLineStyle(2);
      hAfbLinearityMaxBiasFracStat_rebinnedA->Draw("sames");
      hAfbLinearityAvgBiasFracStat_rebinnedA->SetLineColor(3);
      if (dosimpleplots) hAfbLinearityAvgBiasFracStat_rebinnedA->SetLineWidth(0);
      hAfbLinearityAvgBiasFracStat_rebinnedA->SetLineStyle(2);
      hAfbLinearityAvgBiasFracStat_rebinnedA->Draw("sames");
      c2Afbr->cd(1);
      gPad->SetLogy(0);
      hAfbLinearityMaxBiasSlopeMax_rebinnedB->SetLineColor(2);
      if (dosimpleplots) hAfbLinearityMaxBiasSlopeMax_rebinnedB->SetLineWidth(0);
      hAfbLinearityMaxBiasSlopeMax_rebinnedB->GetXaxis()->SetTitle("Variable");
      hAfbLinearityMaxBiasSlopeMax_rebinnedB->GetYaxis()->SetTitle("Measured coefficient linearity slope");
      if (hAfbLinearityMaxBiasSlopeMax_rebinnedB->GetMaximum() < hAfbLinearityMaxBiasSlopeMax_rebinnedA->GetMaximum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB->SetMaximum(hAfbLinearityMaxBiasSlopeMax_rebinnedA->GetMaximum());
      if (hAfbLinearityMaxBiasSlopeMax_rebinnedB->GetMinimum() > hAfbLinearityMaxBiasSlope_rebinnedA->GetMinimum()) hAfbLinearityMaxBiasSlopeMax_rebinnedB->SetMinimum(hAfbLinearityMaxBiasSlope_rebinnedA->GetMinimum());
      if (hAfbLinearityMaxBiasSlopeMax_rebinnedB->GetMaximum() > 1 && hAfbLinearityMaxBiasSlopeMax_rebinnedB->GetMinimum() > 1. - (hAfbLinearityMaxBiasSlopeMax_rebinnedB->GetMaximum() - 1.) / 8.) hAfbLinearityMaxBiasSlopeMax_rebinnedB->SetMinimum(1. - (hAfbLinearityMaxBiasSlopeMax_rebinnedB->GetMaximum() - 1.) / 8.);
      hAfbLinearityMaxBiasSlopeMax_rebinnedB->Draw();

      hAfbLinearityMaxBiasSlopeMax_rebinnedB->SetTitle("");
      c2Afbr->Update();
      hAfbLinearityMaxBiasSlopeMax_rebinnedB->GetXaxis()->SetTitle("");
      hAfbLinearityMaxBiasSlopeMax_rebinnedB->GetXaxis()->SetTitleOffset(2.0);
      hAfbLinearityMaxBiasSlopeMax_rebinnedB->GetXaxis()->SetTitleSize(0.04);
      hAfbLinearityMaxBiasSlopeMax_rebinnedB->GetXaxis()->SetLabelSize(0.03);
      // hAfbLinearityMaxBiasSlopeMax_rebinnedB->GetYaxis()->SetTitleOffset(1.8);
      hAfbLinearityMaxBiasSlopeMax_rebinnedB->GetYaxis()->SetTitleSize(0.04);
      hAfbLinearityMaxBiasSlopeMax_rebinnedB->GetYaxis()->SetLabelSize(0.04);
      for (unsigned int i_bin = 1; i_bin <= CoefNames2.size(); i_bin++) {
          hAfbLinearityMaxBiasSlopeMax_rebinnedB->GetXaxis()->SetBinLabel(i_bin, CoefNames2.at(i_bin - 1).Data());
      }
      hAfbLinearityMaxBiasSlopeMax_rebinnedB->GetXaxis()->LabelsOption("v");
      // hAfbLinearityMaxBiasSlopeMax_rebinnedB->GetYaxis()->LabelsOption("v");
      // DrawCMSLabels(0.04);
      if (CHAN != "combined") DrawDecayChLabel(CHAN, 0.04);

      hAfbLinearityMaxBiasSlope_rebinnedB->Draw("sames");
      hAfbLinearityAvgBiasSlope_rebinnedB->SetLineColor(3);
      if (dosimpleplots) hAfbLinearityAvgBiasSlope_rebinnedB->SetLineWidth(0);
      hAfbLinearityAvgBiasSlope_rebinnedB->Draw("sames");
      hAfbLinearityMaxBiasSlopeMax_rebinnedA->SetLineColor(2);
      if (dosimpleplots) hAfbLinearityMaxBiasSlopeMax_rebinnedA->SetLineWidth(0);
      hAfbLinearityMaxBiasSlopeMax_rebinnedA->SetLineStyle(2);
      hAfbLinearityMaxBiasSlopeMax_rebinnedA->Draw("sames");
      hAfbLinearityMaxBiasSlope_rebinnedA->SetLineStyle(2);
      hAfbLinearityMaxBiasSlope_rebinnedA->Draw("sames");
      hAfbLinearityAvgBiasSlope_rebinnedA->SetLineColor(3);
      if (dosimpleplots) hAfbLinearityAvgBiasSlope_rebinnedA->SetLineWidth(0);
      hAfbLinearityAvgBiasSlope_rebinnedA->SetLineStyle(2);
      hAfbLinearityAvgBiasSlope_rebinnedA->Draw("sames");

      c2Afbr->Print(Form("%sPlots/Regularisation_bias_Afb_rebinned_summary.pdf", TUNFFILE.Data()));
      c2Afbr->Print(Form("%sPlots/Regularisation_bias_Afb_rebinned_summary.root", TUNFFILE.Data()));
      c2Afbr->Print(Form("%sPlots/Regularisation_bias_Afb_rebinned_summary.C", TUNFFILE.Data()));
      c2Afbr->Close();
    }


    // **************
    // Ensemble Tests
    // **************

    if(fEnsembleTest){
      std::cout << "Performing ensemble test" << std::endl;
      
      TH1D *hPseudoUnfDistPerBin[unf_nbins];
      TH1D *hPseudoUnfPullPerBin[unf_nbins];
      TH1D *hPseudoUnfrebinnedAPullPerBin[unf_nbins_rebinnedS];
      TH1D *hPseudoUnfrebinnedAErrorPerBin[unf_nbins_rebinnedS];
      TH1D *hPseudoUnfrebinnedAAfbPullPerBin[unf_nbins_rebinnedS];
      TH1D *hPseudoUnfrebinnedAAfbErrorPerBin[unf_nbins_rebinnedS];
      TH1D *hPseudoUnfrebinnedAOtherPerBin[unf_nbins_rebinnedS];
      TH1D *hPseudoUnfrebinnedBPullPerBin[unf_nbins_rebinnedS];
      TH1D *hPseudoUnfrebinnedBErrorPerBin[unf_nbins_rebinnedS];
      TH1D *hPseudoUnfrebinnedBAfbPullPerBin[unf_nbins_rebinnedS];
      TH1D *hPseudoUnfrebinnedBAfbErrorPerBin[unf_nbins_rebinnedS];
      TH1D *hPseudoUnfrebinnedBOtherPerBin[unf_nbins_rebinnedS];

      TH1D *hPseudoUnfResults;
      TH1D *hPseudoUnfPullResults;
      TH1D *hPseudoUnfMeans;
      TH1D *hPseudoUnfWidths;

      TH1D *hPseudoUnfPulls;
      TH1D *hPseudoUnfPullWidths;
      TH1D *hPseudoUnfrebinnedAPullResults;
      TH1D *hPseudoUnfrebinnedAPulls;
      TH1D *hPseudoUnfrebinnedAPullWidths;
      TH1D *hPseudoUnfrebinnedAAfbPullResults;
      TH1D *hPseudoUnfrebinnedAAfbPulls;
      TH1D *hPseudoUnfrebinnedAAfbPullWidths;
      TH1D *hPseudoUnfrebinnedBPullResults;
      TH1D *hPseudoUnfrebinnedBPulls;
      TH1D *hPseudoUnfrebinnedBPullWidths;
      TH1D *hPseudoUnfrebinnedBAfbPullResults;
      TH1D *hPseudoUnfrebinnedBAfbPulls;
      TH1D *hPseudoUnfrebinnedBAfbPullWidths;

      TH1D *hPseudoInput;

      hPseudoUnfMeans  = (TH1D*)hGenMC->Clone(Form("hPseudoUnfMeans_%s",fVariableNames[i].Data()));
      hPseudoUnfWidths = (TH1D*)hGenMC->Clone(Form("hPseudoUnfWidths_%s",fVariableNames[i].Data()));
      hPseudoUnfPulls  = (TH1D*)hGenMC->Clone(Form("hPseudoUnfPulls_%s",fVariableNames[i].Data()));
      hPseudoUnfResults     = (TH1D*)hGenMC->Clone(Form("hPseudoUnfResults_%s",fVariableNames[i].Data()));
      hPseudoUnfPullResults = (TH1D*)hGenMC->Clone(Form("hPseudoUnfPullResults_%s",fVariableNames[i].Data()));
      hPseudoUnfPullWidths  = (TH1D*)hGenMC->Clone(Form("hPseudoUnfPullWidths_%s",fVariableNames[i].Data()));

      hPseudoUnfrebinnedAPullResults = (TH1D*)hGenMC_rebinnedS->Clone(Form("hPseudoUnfrebinnedAPullResults_%s",fVariableNames[i].Data()));
      hPseudoUnfrebinnedAPulls       = (TH1D*)hGenMC_rebinnedS->Clone(Form("hPseudoUnfrebinnedAPulls_%s",fVariableNames[i].Data()));
      hPseudoUnfrebinnedAPullWidths  = (TH1D*)hGenMC_rebinnedS->Clone(Form("hPseudoUnfrebinnedAPullWidths_%s",fVariableNames[i].Data()));
      hPseudoUnfrebinnedAAfbPullResults = (TH1D*)hGenMC_rebinnedS->Clone(Form("hPseudoUnfrebinnedAAfbPullResults_%s",fVariableNames[i].Data()));
      hPseudoUnfrebinnedAAfbPulls       = (TH1D*)hGenMC_rebinnedS->Clone(Form("hPseudoUnfrebinnedAAfbPulls_%s",fVariableNames[i].Data()));
      hPseudoUnfrebinnedAAfbPullWidths  = (TH1D*)hGenMC_rebinnedS->Clone(Form("hPseudoUnfrebinnedAAfbPullWidths_%s",fVariableNames[i].Data()));

      hPseudoUnfrebinnedBPullResults = (TH1D*)hGenMC_rebinnedS->Clone(Form("hPseudoUnfrebinnedBPullResults_%s",fVariableNames[i].Data()));
      hPseudoUnfrebinnedBPulls       = (TH1D*)hGenMC_rebinnedS->Clone(Form("hPseudoUnfrebinnedBPulls_%s",fVariableNames[i].Data()));
      hPseudoUnfrebinnedBPullWidths  = (TH1D*)hGenMC_rebinnedS->Clone(Form("hPseudoUnfrebinnedBPullWidths_%s",fVariableNames[i].Data()));
      hPseudoUnfrebinnedBAfbPullResults = (TH1D*)hGenMC_rebinnedS->Clone(Form("hPseudoUnfrebinnedBAfbPullResults_%s",fVariableNames[i].Data()));
      hPseudoUnfrebinnedBAfbPulls       = (TH1D*)hGenMC_rebinnedS->Clone(Form("hPseudoUnfrebinnedBAfbPulls_%s",fVariableNames[i].Data()));
      hPseudoUnfrebinnedBAfbPullWidths  = (TH1D*)hGenMC_rebinnedS->Clone(Form("hPseudoUnfrebinnedBAfbPullWidths_%s",fVariableNames[i].Data()));

      // Clear results vectors to prevent crash when looping over variables                                                                                                                                                                                                                                               
      fPseudoUnfResults.clear();
      fPseudoUnfResultsRebinnedA.clear();
      fPseudoUnfResultsRebinnedB.clear();
      fPseudoUnfAfbResultsRebinnedA.clear();
      fPseudoUnfAfbResultsRebinnedB.clear();
      fPseudoUnfOtherResultsRebinnedA.clear();
      fPseudoUnfOtherResultsRebinnedB.clear();

      // We don't currently vary the background in the PEs 
      // so set background stat uncertainty to zero so it doesn't bias the pull widths
      for(int x_bin=1; x_bin <= hBgMC->GetNbinsX(); x_bin++){
        hBgMC->SetBinError(x_bin,0);
        for(int j=0; j < nbgfiles; j++) vecBgMC.at(j)->SetBinError(x_bin,0);
      }

      // Separately create the first "pseudoexperiment" in the vector to find the central value before Poisson fluctuations
      if(fSelfConsistency) {
        hPseudoInput = (TH1D*)hRecMC->Clone(Form("hPseudoInput_%s",fVariableNames[i].Data()));

        if(fSubtrBgMCFromData || fDoSubtrBg) {
          hPseudoInput->Add(hBgMC);
        }

        if(!doBootstrap) {
          for(int x_bin=1; x_bin<=hPseudoInput->GetNbinsX(); x_bin++){
            hPseudoInput->SetBinError(x_bin,TMath::Sqrt(hPseudoInput->GetBinContent(x_bin)));
          }
        }

      }
      else hPseudoInput = (TH1D*)hData->Clone(Form("hPseudoInput_%s",fVariableNames[i].Data()));

      TH1D* hPseudoInput_centralvalue = (TH1D*) hPseudoInput->Clone( Form("%s_centralvalue",hPseudoInput->GetName()) );
      
      if(fSubtrBgMCFromData && !fDoSubtrBg) {
        hPseudoInput_centralvalue->Add(hBgMC, -1.); 
        // TODO: should really use regularBgSubtMethod here
      }

      unfoldData(pttMatUnf, fVariableNames[i]+"_pseudo_centralvalue", hPseudoInput_centralvalue, hMigMat, hGenMC, hVisGenMC, hRecMC, vecBgMC, tau_min,tau_max, kTRUE, CHAN, vecMigMatSys);  

      // *************
      // Bootstrapping
      // *************
       
      if(doBootstrap) {
        // Fill single histogram with results of all variables
        for (int i_bin = 0; i_bin < unf_nbins_rebinnedS; ++i_bin)
        {
          hAllVarUnf_rebinnedA->SetBinContent(1 + i_bin + unf_nbins_rebinnedS*i, fPseudoUnfResultsRebinnedA[0]->GetBinContent(i_bin+1));
          hAllVarUnf_rebinnedA->SetBinError(1   + i_bin + unf_nbins_rebinnedS*i, fPseudoUnfResultsRebinnedA[0]->GetBinError(i_bin+1));
          hAllVarUnf_rebinnedB->SetBinContent(1 + i_bin + unf_nbins_rebinnedS*i, fPseudoUnfResultsRebinnedB[0]->GetBinContent(i_bin+1));
          hAllVarUnf_rebinnedB->SetBinError(1   + i_bin + unf_nbins_rebinnedS*i, fPseudoUnfResultsRebinnedB[0]->GetBinError(i_bin+1));
        }

        // Fill factors used in combining binned asymmetries into coefficients
        // TODO: Fix this for multidimensional shit
        std::vector<double> factor {1.,1.,1.};

        MatrixUnf::AtoCfactor(fVariableNames[i], hGenMC_rebinnedS, factor, 1, 0);

        // Amandeep : What do these do ?
        hist_acombfactor_rebinnedA[0]->SetBinContent(i+1, (fPseudoUnfAfbResultsRebinnedA[0]->GetBinContent(0))/factor[0] );
        hist_acombfactor_rebinnedA[1]->SetBinContent(i+1, (fPseudoUnfAfbResultsRebinnedA[0]->GetBinContent(unf_nbins_rebinnedS+1))/factor[1] );
        hist_acombfactor_rebinnedA[2]->SetBinContent(i+1, (1. - fPseudoUnfAfbResultsRebinnedA[0]->GetBinContent(0) - fPseudoUnfAfbResultsRebinnedA[0]->GetBinContent(unf_nbins_rebinnedS+1))/factor[2] );

        hist_acombfactor_rebinnedB[0]->SetBinContent(i+1, (fPseudoUnfAfbResultsRebinnedB[0]->GetBinContent(0))/factor[0] );
        hist_acombfactor_rebinnedB[1]->SetBinContent(i+1, (fPseudoUnfAfbResultsRebinnedB[0]->GetBinContent(unf_nbins_rebinnedS+1))/factor[1] );
        hist_acombfactor_rebinnedB[2]->SetBinContent(i+1, (1. - fPseudoUnfAfbResultsRebinnedB[0]->GetBinContent(0) - fPseudoUnfAfbResultsRebinnedB[0]->GetBinContent(unf_nbins_rebinnedS+1))/factor[2] );
      }

      // gMyRandom.SetSeed(101);
      gMyRandom.SetSeed(101+i); // to make sure the PEs are uncorrelated for each variable

      for(int i_pseudo=0; i_pseudo < nPE; i_pseudo++){

        if(!doBootstrap) {
          TH1D* hPseudo;
          MakePseudoHistos(hPseudoInput, hPseudo, fVariableNames[i]);
          if(fSubtrBgMCFromData && !fDoSubtrBg) {
            hPseudo->Add(hBgMC,-1.); // TODO: should really use regularBgSubtMethod here
          }
          unfoldData(pttMatUnf, fVariableNames[i]+Form("_pseudo%d",i_pseudo), hPseudo, hMigMat, hGenMC, hVisGenMC, hRecMC, vecBgMC, tau_min,tau_max, kTRUE, CHAN, vecMigMatSys);
          delete hPseudo;
        }

        else {
          TH1D* hDataBootstrap;

          if(fSelfConsistency) {
            hDataBootstrap = (TH1D*) ((TH1D*)((TFile*)fMCFiles->At(0))->Get(Form("hrecoBootstrap_%s%i",fVariableNames[i].Data(),i_pseudo)))->Clone();
            delete (TH1D*)((TFile*)fMCFiles->At(0))->Get(Form("hrecoBootstrap_%s%i",fVariableNames[i].Data(),i_pseudo));
            
            // Amandeep : Commenting to experiment with external lumi scaling 
            //hDataBootstrap->Scale(getLumiWeight( (TFile*) fMCFiles->At(0), fXSectionMC.at(0), fLumiMC.at(0)));
            // End

            // hDataBootstrap = (TH1D*) ((TH1D*)((TFile*)fMCFiles->At(0))->Get(Form("hrecoBootstrap_%s%i",fVariableNames[i].Data(),i_pseudo)))->Clone();
            // std::cout << i_pseudo << " " << ((TH1D*)((TFile*)fMCFiles->At(0))->Get(Form("hrecoBootstrap_%s%i",fVariableNames[i].Data(),i_pseudo)))->GetName() << std::endl;
            for(int j=1; j < nMCfiles; j++){
              TH1D* hDataBootstraptemp = (TH1D*)((TFile*)fMCFiles->At(j))->Get(Form("hrecoBootstrap_%s%i",fVariableNames[i].Data(),i_pseudo));
              // Amandeep : Commenting to experiment with external lumi scaling 
              //hDataBootstraptemp->Scale(getLumiWeight( (TFile*) fMCFiles->At(j), fXSectionMC.at(j), fLumiMC.at(j)));
              // End
              hDataBootstrap->Add(hDataBootstraptemp);
              delete hDataBootstraptemp;
            }
          }

          else {      
            hDataBootstrap = (TH1D*) ((TH1D*)((TFile*)fDataFiles->At(0))->Get(Form("hrecoBootstrap_%s%i",fVariableNames[i].Data(),i_pseudo)))->Clone();
            delete (TH1D*)((TFile*)fDataFiles->At(0))->Get(Form("hrecoBootstrap_%s%i",fVariableNames[i].Data(),i_pseudo));
            
            // std::cout << i_pseudo << " " << ((TH1D*)((TFile*)fDataFiles->At(0))->Get(Form("hrecoBootstrap_%s%i",fVariableNames[i].Data(),i_pseudo)))->GetName() << std::endl;
            
            for(int j=1; j < ndatafiles; j++){
              hDataBootstrap->Add((TH1D*)((TFile*)fDataFiles->At(j))->Get(Form("hrecoBootstrap_%s%i",fVariableNames[i].Data(),i_pseudo)));
              delete (TH1D*)((TFile*)fDataFiles->At(j))->Get(Form("hrecoBootstrap_%s%i",fVariableNames[i].Data(),i_pseudo));
            }
          }

          if(fSelfConsistency && (fSubtrBgMCFromData || fDoSubtrBg) ) {
            hDataBootstrap->Add(hBgMC);
          }
          if(fSubtrBgMCFromData && !fDoSubtrBg) {
            hDataBootstrap->Add(hBgMC,-1.); // TODO: should really use regularBgSubtMethod here
          }

          unfoldData(pttMatUnf, fVariableNames[i]+Form("_pseudo%d",i_pseudo), hDataBootstrap, hMigMat, hGenMC, hVisGenMC, hRecMC, vecBgMC, tau_min,tau_max, kTRUE, CHAN, vecMigMatSys);
          delete hDataBootstrap;
        }
      }

      std::cout<<"Finished creating ensembles"<<std::endl;
      
      /*
      for(int j=0;j<ndatafiles;j++){
        ((TFile*)fDataFiles->At(j))->Close("R");
        ((TFile*)fDataFiles->At(j))->ReOpen("READ");
      }
      */

      // Create gen-level histo for asymmetries and coefficient
      // TH1D* hAfbGenMC_rebinnedS = (TH1D*)hGenMC_rebinnedS->Clone(Form("%s_AfbGen_rebinnedS",fVariableNames[i].Data()));
      
      TH1D *hAfbGenMC_rebinnedS    = new TH1D(Form("%s_AfbGen_rebinnedS", fVariableNames[i].Data()), "", (numRebinnedCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      MatrixUnf::fillAsymsAndCoefficient(hGenMC_rebinnedS, nullptr, hGenMC_rebinnedS, fVariableNames[i], hAfbGenMC_rebinnedS, nullptr, true, true, numRebinnedCoefficients);

      double nsigmas = 5.;
      
      for(int i_bin = 0; i_bin<=unf_nbins+1; i_bin++){
        std::cout << i_bin << " " << unf_nbins << std::endl;
        std::cout << fPseudoUnfResults.size()  << std::endl;
        std::cout << fPseudoUnfResults.at(0)->GetName() << " " << fPseudoUnfResults.at(1)->GetName() << std::endl;
        
        double bincontent = fPseudoUnfResults[0]->GetBinContent(i_bin);
        double binerror   = fPseudoUnfResults[0]->GetBinError(i_bin);
        
        std::cout << bincontent << " " << binerror << " " << hGenMC->GetBinContent(i_bin) << " " << hGenMC->GetBinError(i_bin) << std::endl;

        hPseudoUnfDistPerBin[i_bin] = new TH1D(Form("hPseudoUnfDistPerBin_%d",i_bin), Form("hPseudoUnfDistPerBin_%d",i_bin), 1000,bincontent-nsigmas*binerror, bincontent+nsigmas*binerror);
        hPseudoUnfPullPerBin[i_bin] = new TH1D(Form("hPseudoUnfPullPerBin_%d",i_bin), Form("hPseudoUnfPullPerBin_%d",i_bin), 1000, -nsigmas, nsigmas);
        
        if(i_bin<=unf_nbins_rebinnedS+1) {
          double OtherMean_rebinnedA  = 0;
          double OtherError_rebinnedA = 0;
          double ErrorMean_rebinnedA  = 0;
          double ErrorError_rebinnedA = 0;
          double AfbErrorMean_rebinnedA  = 0;
          double AfbErrorError_rebinnedA = 0;

          double OtherMean_rebinnedB  = 0;
          double OtherError_rebinnedB = 0;
          double ErrorMean_rebinnedB  = 0;
          double ErrorError_rebinnedB = 0;
          double AfbErrorMean_rebinnedB  = 0;
          double AfbErrorError_rebinnedB = 0;

          MatrixUnfControl::CalculateHistRange(fPseudoUnfOtherResultsRebinnedA, i_bin, OtherMean_rebinnedA, OtherError_rebinnedA, kFALSE);
          MatrixUnfControl::CalculateHistRange(fPseudoUnfResultsRebinnedA, i_bin, ErrorMean_rebinnedA, ErrorError_rebinnedA, kTRUE);
          MatrixUnfControl::CalculateHistRange(fPseudoUnfAfbResultsRebinnedA, i_bin, AfbErrorMean_rebinnedA, AfbErrorError_rebinnedA, kTRUE);

          MatrixUnfControl::CalculateHistRange(fPseudoUnfOtherResultsRebinnedB, i_bin, OtherMean_rebinnedB, OtherError_rebinnedB, kFALSE);
          MatrixUnfControl::CalculateHistRange(fPseudoUnfResultsRebinnedB, i_bin, ErrorMean_rebinnedB, ErrorError_rebinnedB, kTRUE);
          MatrixUnfControl::CalculateHistRange(fPseudoUnfAfbResultsRebinnedB, i_bin, AfbErrorMean_rebinnedB, AfbErrorError_rebinnedB, kTRUE);

          hPseudoUnfrebinnedAPullPerBin[i_bin]     = new TH1D(Form("hPseudoUnfrebinnedAPullPerBin_%d",i_bin), Form("hPseudoUnfrebinnedAPullPerBin_%d",i_bin), 1000,-nsigmas, nsigmas);
          hPseudoUnfrebinnedAErrorPerBin[i_bin]    = new TH1D(Form("hPseudoUnfrebinnedAErrorPerBin_%d",i_bin), Form("hPseudoUnfrebinnedAErrorPerBin_%d",i_bin), 1000,ErrorMean_rebinnedA-ErrorError_rebinnedA*nsigmas, ErrorMean_rebinnedA+ErrorError_rebinnedA*nsigmas);
          hPseudoUnfrebinnedAAfbPullPerBin[i_bin]  = new TH1D(Form("hPseudoUnfrebinnedAAfbPullPerBin_%d",i_bin), Form("hPseudoUnfrebinnedAAfbPullPerBin_%d",i_bin), 1000,-nsigmas, nsigmas);
          hPseudoUnfrebinnedAAfbErrorPerBin[i_bin] = new TH1D(Form("hPseudoUnfrebinnedAAfbErrorPerBin_%d",i_bin), Form("hPseudoUnfrebinnedAAfbErrorPerBin_%d",i_bin), 1000,AfbErrorMean_rebinnedA-AfbErrorError_rebinnedA*nsigmas, AfbErrorMean_rebinnedA+AfbErrorError_rebinnedA*nsigmas);
          hPseudoUnfrebinnedAOtherPerBin[i_bin]    = new TH1D(Form("hPseudoUnfrebinnedAOtherPerBin_%d",i_bin), Form("hPseudoUnfrebinnedAOtherPerBin_%d",i_bin), 1000,OtherMean_rebinnedA-nsigmas*OtherError_rebinnedA, OtherMean_rebinnedA+nsigmas*OtherError_rebinnedA);
          hPseudoUnfrebinnedBPullPerBin[i_bin]     = new TH1D(Form("hPseudoUnfrebinnedBPullPerBin_%d",i_bin), Form("hPseudoUnfrebinnedBPullPerBin_%d",i_bin), 1000,-nsigmas, nsigmas);
          hPseudoUnfrebinnedBErrorPerBin[i_bin]    = new TH1D(Form("hPseudoUnfrebinnedBErrorPerBin_%d",i_bin), Form("hPseudoUnfrebinnedBErrorPerBin_%d",i_bin), 1000,ErrorMean_rebinnedB-ErrorError_rebinnedB*nsigmas, ErrorMean_rebinnedB+ErrorError_rebinnedB*nsigmas);
          hPseudoUnfrebinnedBAfbPullPerBin[i_bin]  = new TH1D(Form("hPseudoUnfrebinnedBAfbPullPerBin_%d",i_bin), Form("hPseudoUnfrebinnedBAfbPullPerBin_%d",i_bin), 1000,-nsigmas, nsigmas);
          hPseudoUnfrebinnedBAfbErrorPerBin[i_bin] = new TH1D(Form("hPseudoUnfrebinnedBAfbErrorPerBin_%d",i_bin), Form("hPseudoUnfrebinnedBAfbErrorPerBin_%d",i_bin), 1000,AfbErrorMean_rebinnedB-AfbErrorError_rebinnedB*nsigmas, AfbErrorMean_rebinnedB+AfbErrorError_rebinnedB*nsigmas);
          hPseudoUnfrebinnedBOtherPerBin[i_bin]    = new TH1D(Form("hPseudoUnfrebinnedBOtherPerBin_%d",i_bin), Form("hPseudoUnfrebinnedBOtherPerBin_%d",i_bin), 1000,OtherMean_rebinnedB-nsigmas*OtherError_rebinnedB, OtherMean_rebinnedB+nsigmas*OtherError_rebinnedB);
        }

        // Exclude 0 from the loop because this is the central value result without Poisson fluctuations   
        for(unsigned int i_pseudo=1; i_pseudo<fPseudoUnfResults.size(); i_pseudo++){
          
          hPseudoUnfDistPerBin[i_bin]->Fill(fPseudoUnfResults[i_pseudo]->GetBinContent(i_bin));

          // Calculate the pull distribution wrt the known central value, so the pull mean is 0 by construction
          // unless doing fSelfConsistency, in which case use the known gen-level value 
          if (fSelfConsistency) hPseudoUnfPullPerBin[i_bin]->Fill( fPseudoUnfResults[i_pseudo]->GetBinError(i_bin) > 0 ? (fPseudoUnfResults[i_pseudo]->GetBinContent(i_bin) - hGenMC->GetBinContent(i_bin) )/ fPseudoUnfResults[i_pseudo]->GetBinError(i_bin) : 0 );
          else hPseudoUnfPullPerBin[i_bin]->Fill( fPseudoUnfResults[i_pseudo]->GetBinError(i_bin) > 0 ? (fPseudoUnfResults[i_pseudo]->GetBinContent(i_bin) - fPseudoUnfResults[0]->GetBinContent(i_bin) )/ fPseudoUnfResults[i_pseudo]->GetBinError(i_bin) : 0 );

          if(i_bin<=unf_nbins_rebinnedS+1) {
            if(fSelfConsistency) {
              hPseudoUnfrebinnedAPullPerBin[i_bin]->Fill( fPseudoUnfResultsRebinnedA[i_pseudo]->GetBinError(i_bin) > 0 ? (fPseudoUnfResultsRebinnedA[i_pseudo]->GetBinContent(i_bin) - hGenMC_rebinnedS->GetBinContent(i_bin) )/ fPseudoUnfResultsRebinnedA[i_pseudo]->GetBinError(i_bin) : 0 );
              hPseudoUnfrebinnedBPullPerBin[i_bin]->Fill( fPseudoUnfResultsRebinnedB[i_pseudo]->GetBinError(i_bin) > 0 ? (fPseudoUnfResultsRebinnedB[i_pseudo]->GetBinContent(i_bin) - hGenMC_rebinnedS->GetBinContent(i_bin) )/ fPseudoUnfResultsRebinnedB[i_pseudo]->GetBinError(i_bin) : 0 );
              hPseudoUnfrebinnedAAfbPullPerBin[i_bin]->Fill( fPseudoUnfAfbResultsRebinnedA[i_pseudo]->GetBinError(i_bin) > 0 ? (fPseudoUnfAfbResultsRebinnedA[i_pseudo]->GetBinContent(i_bin) - hAfbGenMC_rebinnedS->GetBinContent(i_bin) )/ fPseudoUnfAfbResultsRebinnedA[i_pseudo]->GetBinError(i_bin) : 0 );
              hPseudoUnfrebinnedBAfbPullPerBin[i_bin]->Fill( fPseudoUnfAfbResultsRebinnedB[i_pseudo]->GetBinError(i_bin) > 0 ? (fPseudoUnfAfbResultsRebinnedB[i_pseudo]->GetBinContent(i_bin) - hAfbGenMC_rebinnedS->GetBinContent(i_bin) )/ fPseudoUnfAfbResultsRebinnedB[i_pseudo]->GetBinError(i_bin) : 0 );
            }
            else {
              hPseudoUnfrebinnedAPullPerBin[i_bin]->Fill( fPseudoUnfResultsRebinnedA[i_pseudo]->GetBinError(i_bin) > 0 ? (fPseudoUnfResultsRebinnedA[i_pseudo]->GetBinContent(i_bin) - fPseudoUnfResultsRebinnedA[0]->GetBinContent(i_bin) )/ fPseudoUnfResultsRebinnedA[i_pseudo]->GetBinError(i_bin) : 0 );
              hPseudoUnfrebinnedBPullPerBin[i_bin]->Fill( fPseudoUnfResultsRebinnedB[i_pseudo]->GetBinError(i_bin) > 0 ? (fPseudoUnfResultsRebinnedB[i_pseudo]->GetBinContent(i_bin) - fPseudoUnfResultsRebinnedB[0]->GetBinContent(i_bin) )/ fPseudoUnfResultsRebinnedB[i_pseudo]->GetBinError(i_bin) : 0 );
              hPseudoUnfrebinnedAAfbPullPerBin[i_bin]->Fill( fPseudoUnfAfbResultsRebinnedA[i_pseudo]->GetBinError(i_bin) > 0 ? (fPseudoUnfAfbResultsRebinnedA[i_pseudo]->GetBinContent(i_bin) - fPseudoUnfAfbResultsRebinnedA[0]->GetBinContent(i_bin) )/ fPseudoUnfAfbResultsRebinnedA[i_pseudo]->GetBinError(i_bin) : 0 );
              hPseudoUnfrebinnedBAfbPullPerBin[i_bin]->Fill( fPseudoUnfAfbResultsRebinnedB[i_pseudo]->GetBinError(i_bin) > 0 ? (fPseudoUnfAfbResultsRebinnedB[i_pseudo]->GetBinContent(i_bin) - fPseudoUnfAfbResultsRebinnedB[0]->GetBinContent(i_bin) )/ fPseudoUnfAfbResultsRebinnedB[i_pseudo]->GetBinError(i_bin) : 0 );
            }
            hPseudoUnfrebinnedAErrorPerBin[i_bin]->Fill( fPseudoUnfResultsRebinnedA[i_pseudo]->GetBinError(i_bin) );
            hPseudoUnfrebinnedBErrorPerBin[i_bin]->Fill( fPseudoUnfResultsRebinnedB[i_pseudo]->GetBinError(i_bin) );
            hPseudoUnfrebinnedAAfbErrorPerBin[i_bin]->Fill( fPseudoUnfAfbResultsRebinnedA[i_pseudo]->GetBinError(i_bin) );
            hPseudoUnfrebinnedBAfbErrorPerBin[i_bin]->Fill( fPseudoUnfAfbResultsRebinnedB[i_pseudo]->GetBinError(i_bin) );
            hPseudoUnfrebinnedAOtherPerBin[i_bin]->Fill( fPseudoUnfOtherResultsRebinnedA[i_pseudo]->GetBinContent(i_bin) );
            hPseudoUnfrebinnedBOtherPerBin[i_bin]->Fill( fPseudoUnfOtherResultsRebinnedB[i_pseudo]->GetBinContent(i_bin) );

            if(doBootstrap && i_bin>0 && i_bin<=unf_nbins_rebinnedS) {
              BinPE_rebinnedA[i_bin-1][i][i_pseudo-1]  = fPseudoUnfResultsRebinnedA[i_pseudo]->GetBinError(i_bin) > 0    ? (fPseudoUnfResultsRebinnedA[i_pseudo]->GetBinContent(i_bin)    - fPseudoUnfResultsRebinnedA[0]->GetBinContent(i_bin) )/ fPseudoUnfResultsRebinnedA[i_pseudo]->GetBinError(i_bin) : 0;
              CoefPE_rebinnedA[i_bin-1][i][i_pseudo-1] = fPseudoUnfAfbResultsRebinnedA[i_pseudo]->GetBinError(i_bin) > 0 ? (fPseudoUnfAfbResultsRebinnedA[i_pseudo]->GetBinContent(i_bin) - fPseudoUnfAfbResultsRebinnedA[0]->GetBinContent(i_bin) )/ fPseudoUnfAfbResultsRebinnedA[i_pseudo]->GetBinError(i_bin) : 0;

              BinPE_rebinnedB[i_bin-1][i][i_pseudo-1]  = fPseudoUnfResultsRebinnedB[i_pseudo]->GetBinError(i_bin) > 0    ? (fPseudoUnfResultsRebinnedB[i_pseudo]->GetBinContent(i_bin)    - fPseudoUnfResultsRebinnedB[0]->GetBinContent(i_bin) )/ fPseudoUnfResultsRebinnedB[i_pseudo]->GetBinError(i_bin) : 0;
              CoefPE_rebinnedB[i_bin-1][i][i_pseudo-1] = fPseudoUnfAfbResultsRebinnedB[i_pseudo]->GetBinError(i_bin) > 0 ? (fPseudoUnfAfbResultsRebinnedB[i_pseudo]->GetBinContent(i_bin) - fPseudoUnfAfbResultsRebinnedB[0]->GetBinContent(i_bin) )/ fPseudoUnfAfbResultsRebinnedB[i_pseudo]->GetBinError(i_bin) : 0;
            }
          }
        }

        hPseudoUnfDistPerBin[i_bin]->BufferEmpty(1);
        hPseudoUnfDistPerBin[i_bin]->Sumw2();

        hPseudoUnfPullPerBin[i_bin]->BufferEmpty(1);
        hPseudoUnfPullPerBin[i_bin]->Sumw2();

        std::cout << "Unfolded results for central value: "         << fVariableNames[i].Data() << " " << i_bin << " bin ----> Value: "<<fPseudoUnfResults[0]->GetBinContent(i_bin)<<" +- "<<fPseudoUnfResults[0]->GetBinError(i_bin)<<std::endl;
        std::cout << "Unfolded results from Pseudo experiment for " << fVariableNames[i].Data() << " " << i_bin << " bin ----> Mean: "<<hPseudoUnfDistPerBin[i_bin]->GetMean()<<" +- "<<hPseudoUnfDistPerBin[i_bin]->GetMeanError()<<" and RMS: "<<hPseudoUnfDistPerBin[i_bin]->GetRMS()<<" +- "<<hPseudoUnfDistPerBin[i_bin]->GetRMSError()<<std::endl;
        std::cout << "Unfolded pulls from Pseudo experiment for   " << fVariableNames[i].Data() << " " << i_bin << " bin ----> Pull Mean: "<<hPseudoUnfPullPerBin[i_bin]->GetMean()<<" +- "<<hPseudoUnfPullPerBin[i_bin]->GetMeanError()<<" and Pull RMS: "<<hPseudoUnfPullPerBin[i_bin]->GetRMS()<<" +- "<<hPseudoUnfPullPerBin[i_bin]->GetRMSError()<<std::endl;

        hPseudoUnfResults->SetBinContent(i_bin, hPseudoUnfDistPerBin[i_bin]->GetMean());
        hPseudoUnfResults->SetBinError(i_bin, hPseudoUnfDistPerBin[i_bin]->GetRMS());
        hPseudoUnfPullResults->SetBinContent(i_bin, hPseudoUnfPullPerBin[i_bin]->GetMean());
        hPseudoUnfPullResults->SetBinError(i_bin, hPseudoUnfPullPerBin[i_bin]->GetRMS());

        hPseudoUnfMeans->SetBinContent(i_bin, hPseudoUnfDistPerBin[i_bin]->GetMean());
        hPseudoUnfMeans->SetBinError(i_bin, hPseudoUnfDistPerBin[i_bin]->GetMeanError());
        hPseudoUnfWidths->SetBinContent(i_bin, hPseudoUnfDistPerBin[i_bin]->GetRMS());
        hPseudoUnfWidths->SetBinError(i_bin, hPseudoUnfDistPerBin[i_bin]->GetRMSError());

        hPseudoUnfPulls->SetBinContent(i_bin, hPseudoUnfPullPerBin[i_bin]->GetMean());
        hPseudoUnfPulls->SetBinError(i_bin, hPseudoUnfPullPerBin[i_bin]->GetMeanError());
        hPseudoUnfPullWidths->SetBinContent(i_bin, hPseudoUnfPullPerBin[i_bin]->GetRMS());
        hPseudoUnfPullWidths->SetBinError(i_bin, hPseudoUnfPullPerBin[i_bin]->GetRMSError());

        hPseudoUnfDistPerBin[i_bin]->Rebin(40);
        hPseudoUnfPullPerBin[i_bin]->Rebin(40);

        hPseudoUnfDistPerBin[i_bin]->Write();
        hPseudoUnfPullPerBin[i_bin]->Write();

        delete hPseudoUnfDistPerBin[i_bin];
        delete hPseudoUnfPullPerBin[i_bin];

        if(i_bin<=unf_nbins_rebinnedS+1) {
          hPseudoUnfrebinnedAPullPerBin[i_bin]->BufferEmpty(1);
          hPseudoUnfrebinnedAPullPerBin[i_bin]->Sumw2();

          hPseudoUnfrebinnedAAfbPullPerBin[i_bin]->BufferEmpty(1);
          hPseudoUnfrebinnedAAfbPullPerBin[i_bin]->Sumw2();

          hPseudoUnfrebinnedAErrorPerBin[i_bin]->BufferEmpty(1);
          hPseudoUnfrebinnedAErrorPerBin[i_bin]->Sumw2();

          hPseudoUnfrebinnedAAfbErrorPerBin[i_bin]->BufferEmpty(1);
          hPseudoUnfrebinnedAAfbErrorPerBin[i_bin]->Sumw2();

          hPseudoUnfrebinnedAOtherPerBin[i_bin]->BufferEmpty(1);
          hPseudoUnfrebinnedAOtherPerBin[i_bin]->Sumw2();

          hPseudoUnfrebinnedBPullPerBin[i_bin]->BufferEmpty(1);
          hPseudoUnfrebinnedBPullPerBin[i_bin]->Sumw2();

          hPseudoUnfrebinnedBAfbPullPerBin[i_bin]->BufferEmpty(1);
          hPseudoUnfrebinnedBAfbPullPerBin[i_bin]->Sumw2();

          hPseudoUnfrebinnedBErrorPerBin[i_bin]->BufferEmpty(1);
          hPseudoUnfrebinnedBErrorPerBin[i_bin]->Sumw2();

          hPseudoUnfrebinnedBAfbErrorPerBin[i_bin]->BufferEmpty(1);
          hPseudoUnfrebinnedBAfbErrorPerBin[i_bin]->Sumw2();

          hPseudoUnfrebinnedBOtherPerBin[i_bin]->BufferEmpty(1);
          hPseudoUnfrebinnedBOtherPerBin[i_bin]->Sumw2();

          std::cout << "Unfolded (rebinnedA) pulls from Pseudo experiment for "    << fVariableNames[i].Data() << " " << i_bin << " bin ----> Mean: "<< hPseudoUnfrebinnedAPullPerBin[i_bin]->GetMean()<<" +- "<<hPseudoUnfrebinnedAPullPerBin[i_bin]->GetMeanError()<<" and RMS: "<<hPseudoUnfrebinnedAPullPerBin[i_bin]->GetRMS()<<" +- "<<hPseudoUnfrebinnedAPullPerBin[i_bin]->GetRMSError()<<std::endl;
          std::cout << "Unfolded (rebinnedA) Afb pulls from Pseudo experiment for "<< fVariableNames[i].Data() << " " << i_bin << " bin ----> Mean: "<< hPseudoUnfrebinnedAAfbPullPerBin[i_bin]->GetMean()<<" +- "<<hPseudoUnfrebinnedAAfbPullPerBin[i_bin]->GetMeanError()<<" and RMS: "<<hPseudoUnfrebinnedAAfbPullPerBin[i_bin]->GetRMS()<<" +- "<<hPseudoUnfrebinnedAAfbPullPerBin[i_bin]->GetRMSError()<<std::endl;
          std::cout << "Unfolded (rebinnedB) pulls from Pseudo experiment for "    << fVariableNames[i].Data() << " " << i_bin << " bin ----> Mean: "<< hPseudoUnfrebinnedBPullPerBin[i_bin]->GetMean()<<" +- "<<hPseudoUnfrebinnedBPullPerBin[i_bin]->GetMeanError()<<" and RMS: "<<hPseudoUnfrebinnedBPullPerBin[i_bin]->GetRMS()<<" +- "<<hPseudoUnfrebinnedBPullPerBin[i_bin]->GetRMSError()<<std::endl;
          std::cout << "Unfolded (rebinnedB) Afb pulls from Pseudo experiment for "<< fVariableNames[i].Data() << " " << i_bin << " bin ----> Mean: "<< hPseudoUnfrebinnedBAfbPullPerBin[i_bin]->GetMean()<<" +- "<<hPseudoUnfrebinnedBAfbPullPerBin[i_bin]->GetMeanError()<<" and RMS: "<<hPseudoUnfrebinnedBAfbPullPerBin[i_bin]->GetRMS()<<" +- "<<hPseudoUnfrebinnedBAfbPullPerBin[i_bin]->GetRMSError()<<std::endl;

          hPseudoUnfrebinnedAPullResults->SetBinContent(i_bin, hPseudoUnfrebinnedAPullPerBin[i_bin]->GetMean());
          hPseudoUnfrebinnedAPullResults->SetBinError(i_bin, hPseudoUnfrebinnedAPullPerBin[i_bin]->GetRMS());
          hPseudoUnfrebinnedAPulls->SetBinContent(i_bin, hPseudoUnfrebinnedAPullPerBin[i_bin]->GetMean());
          hPseudoUnfrebinnedAPulls->SetBinError(i_bin, hPseudoUnfrebinnedAPullPerBin[i_bin]->GetMeanError());
          hPseudoUnfrebinnedAPullWidths->SetBinContent(i_bin, hPseudoUnfrebinnedAPullPerBin[i_bin]->GetRMS());
          hPseudoUnfrebinnedAPullWidths->SetBinError(i_bin, hPseudoUnfrebinnedAPullPerBin[i_bin]->GetRMSError());
          hPseudoUnfrebinnedAAfbPullResults->SetBinContent(i_bin, hPseudoUnfrebinnedAAfbPullPerBin[i_bin]->GetMean());
          hPseudoUnfrebinnedAAfbPullResults->SetBinError(i_bin, hPseudoUnfrebinnedAAfbPullPerBin[i_bin]->GetRMS());
          hPseudoUnfrebinnedAAfbPulls->SetBinContent(i_bin, hPseudoUnfrebinnedAAfbPullPerBin[i_bin]->GetMean());
          hPseudoUnfrebinnedAAfbPulls->SetBinError(i_bin, hPseudoUnfrebinnedAAfbPullPerBin[i_bin]->GetMeanError());
          hPseudoUnfrebinnedAAfbPullWidths->SetBinContent(i_bin, hPseudoUnfrebinnedAAfbPullPerBin[i_bin]->GetRMS());
          hPseudoUnfrebinnedAAfbPullWidths->SetBinError(i_bin, hPseudoUnfrebinnedAAfbPullPerBin[i_bin]->GetRMSError());

          hPseudoUnfrebinnedBPullResults->SetBinContent(i_bin, hPseudoUnfrebinnedBPullPerBin[i_bin]->GetMean());
          hPseudoUnfrebinnedBPullResults->SetBinError(i_bin, hPseudoUnfrebinnedBPullPerBin[i_bin]->GetRMS());
          hPseudoUnfrebinnedBPulls->SetBinContent(i_bin, hPseudoUnfrebinnedBPullPerBin[i_bin]->GetMean());
          hPseudoUnfrebinnedBPulls->SetBinError(i_bin, hPseudoUnfrebinnedBPullPerBin[i_bin]->GetMeanError());
          hPseudoUnfrebinnedBPullWidths->SetBinContent(i_bin, hPseudoUnfrebinnedBPullPerBin[i_bin]->GetRMS());
          hPseudoUnfrebinnedBPullWidths->SetBinError(i_bin, hPseudoUnfrebinnedBPullPerBin[i_bin]->GetRMSError());
          hPseudoUnfrebinnedBAfbPullResults->SetBinContent(i_bin, hPseudoUnfrebinnedBAfbPullPerBin[i_bin]->GetMean());
          hPseudoUnfrebinnedBAfbPullResults->SetBinError(i_bin, hPseudoUnfrebinnedBAfbPullPerBin[i_bin]->GetRMS());
          hPseudoUnfrebinnedBAfbPulls->SetBinContent(i_bin, hPseudoUnfrebinnedBAfbPullPerBin[i_bin]->GetMean());
          hPseudoUnfrebinnedBAfbPulls->SetBinError(i_bin, hPseudoUnfrebinnedBAfbPullPerBin[i_bin]->GetMeanError());
          hPseudoUnfrebinnedBAfbPullWidths->SetBinContent(i_bin, hPseudoUnfrebinnedBAfbPullPerBin[i_bin]->GetRMS());
          hPseudoUnfrebinnedBAfbPullWidths->SetBinError(i_bin, hPseudoUnfrebinnedBAfbPullPerBin[i_bin]->GetRMSError());

          hPseudoUnfrebinnedAPullPerBin[i_bin]->Rebin(40);
          hPseudoUnfrebinnedAErrorPerBin[i_bin]->Rebin(40);
          hPseudoUnfrebinnedAAfbPullPerBin[i_bin]->Rebin(40);
          hPseudoUnfrebinnedAAfbErrorPerBin[i_bin]->Rebin(40);
          hPseudoUnfrebinnedAOtherPerBin[i_bin]->Rebin(10);
          hPseudoUnfrebinnedBPullPerBin[i_bin]->Rebin(40);
          hPseudoUnfrebinnedBErrorPerBin[i_bin]->Rebin(40);
          hPseudoUnfrebinnedBAfbPullPerBin[i_bin]->Rebin(40);
          hPseudoUnfrebinnedBAfbErrorPerBin[i_bin]->Rebin(40);
          hPseudoUnfrebinnedBOtherPerBin[i_bin]->Rebin(10);

          hPseudoUnfrebinnedAPullPerBin[i_bin]->Write();
          hPseudoUnfrebinnedAErrorPerBin[i_bin]->Write();
          hPseudoUnfrebinnedAAfbPullPerBin[i_bin]->Write();
          hPseudoUnfrebinnedAAfbErrorPerBin[i_bin]->Write();
          hPseudoUnfrebinnedAOtherPerBin[i_bin]->Write();
          hPseudoUnfrebinnedBPullPerBin[i_bin]->Write();
          hPseudoUnfrebinnedBErrorPerBin[i_bin]->Write();
          hPseudoUnfrebinnedBAfbPullPerBin[i_bin]->Write();
          hPseudoUnfrebinnedBAfbErrorPerBin[i_bin]->Write();
          hPseudoUnfrebinnedBOtherPerBin[i_bin]->Write();

          delete hPseudoUnfrebinnedAPullPerBin[i_bin];
          delete hPseudoUnfrebinnedAErrorPerBin[i_bin];
          delete hPseudoUnfrebinnedAAfbPullPerBin[i_bin];
          delete hPseudoUnfrebinnedAAfbErrorPerBin[i_bin];
          delete hPseudoUnfrebinnedAOtherPerBin[i_bin];
          delete hPseudoUnfrebinnedBPullPerBin[i_bin];
          delete hPseudoUnfrebinnedBErrorPerBin[i_bin];
          delete hPseudoUnfrebinnedBAfbPullPerBin[i_bin];
          delete hPseudoUnfrebinnedBAfbErrorPerBin[i_bin];
          delete hPseudoUnfrebinnedBOtherPerBin[i_bin];
        }

      }

      hPseudoUnfResults->Write();
      hPseudoUnfPullResults->Write();
      hPseudoUnfMeans->Write();
      hPseudoUnfWidths->Write();

      hPseudoUnfPulls->Write();
      hPseudoUnfPullWidths->Write();
      hPseudoUnfrebinnedAPullResults->Write();
      hPseudoUnfrebinnedAPulls->Write();
      hPseudoUnfrebinnedAPullWidths->Write();
      hPseudoUnfrebinnedAAfbPullResults->Write();
      hPseudoUnfrebinnedAAfbPulls->Write();
      hPseudoUnfrebinnedAAfbPullWidths->Write();
      hPseudoUnfrebinnedBPullResults->Write();
      hPseudoUnfrebinnedBPulls->Write();
      hPseudoUnfrebinnedBPullWidths->Write();
      hPseudoUnfrebinnedBAfbPullResults->Write();
      hPseudoUnfrebinnedBAfbPulls->Write();
      hPseudoUnfrebinnedBAfbPullWidths->Write();


      delete hPseudoUnfResults;
      delete hPseudoUnfPullResults;
      delete hPseudoUnfMeans;
      delete hPseudoUnfWidths;
      delete hPseudoUnfPulls;
      delete hPseudoUnfPullWidths;
      delete hPseudoUnfrebinnedAPullResults;
      delete hPseudoUnfrebinnedAPulls;
      delete hPseudoUnfrebinnedAPullWidths;
      delete hPseudoUnfrebinnedAAfbPullResults;
      delete hPseudoUnfrebinnedAAfbPulls;
      delete hPseudoUnfrebinnedAAfbPullWidths;
      delete hPseudoUnfrebinnedBPullResults;
      delete hPseudoUnfrebinnedBPulls;
      delete hPseudoUnfrebinnedBPullWidths;
      delete hPseudoUnfrebinnedBAfbPullResults;
      delete hPseudoUnfrebinnedBAfbPulls;
      delete hPseudoUnfrebinnedBAfbPullWidths;

    }


    // Write the MatrixUnf object to file                                                                                                                                         
    printf("\n All unfolded - Writing output file ...\n");
    outfile->cd();
    outfile->Write();
    // outfile->Close();
    delete outfile;

    delete pttMatUnf;
    // return 1;
  }
  // Variable list ends here
  // Total number of events in data before bg subtraction


  // CalculateCorrelationMatrices for all variables in the variable list
  if(doBootstrap && fEnsembleTest) {

    printf("\n Calculating correlation matrices for the %d variables ...\n",nVar);

    // Calculate and fill correlation matrices
    TString BootstrapCorrelationMatricesFilename = TUNFFILE + "BootstrapCorrelationMatrices.root";
    TFile  *BootstrapCorrelationMatricesFile     = new TFile(BootstrapCorrelationMatricesFilename.Data(), "RECREATE");
    BootstrapCorrelationMatricesFile->cd();

    TH1D *hAllVarUnfClone_rebinnedA     = (TH1D*)hAllVarUnf_rebinnedA->Clone("hAllVarUnfClone_rebinnedA");
    TH1D *hAllVarUnfClone_rebinnedB     = (TH1D*)hAllVarUnf_rebinnedB->Clone("hAllVarUnfClone_rebinnedB");

    TH1D *hAllVarUnfCorrected_rebinnedA = (TH1D*)hAllVarUnf_rebinnedA->Clone("hAllVarUnfCorrected_rebinnedA");
    TH1D *hAllVarUnfCorrected_rebinnedB = (TH1D*)hAllVarUnf_rebinnedB->Clone("hAllVarUnfCorrected_rebinnedB");

    std::vector<double> pullwidthBins_rebinnedA, pullwidthVBins_rebinnedA, pullwidthBins_rebinnedB, pullwidthVBins_rebinnedB;

    // Amandeep : Commenting to experiment
    summaryunfoutfile->cd();
    MatrixUnfControl::calculateCorrelationMatrices(nVar, nbinsrhoi, nPE, kFALSE, "rebinnedA", BootstrapCorrelationMatricesFile, hAllVarUnfClone_rebinnedA, hist_acombfactor_rebinnedA, BinPE_rebinnedA, CoefPE_rebinnedA, pullwidthBins_rebinnedA, pullwidthVBins_rebinnedA);
    summaryunfoutfile->cd();
    MatrixUnfControl::calculateCorrelationMatrices(nVar, nbinsrhoi, nPE, kFALSE, "rebinnedB", BootstrapCorrelationMatricesFile, hAllVarUnfClone_rebinnedB, hist_acombfactor_rebinnedB, BinPE_rebinnedB, CoefPE_rebinnedB, pullwidthBins_rebinnedB, pullwidthVBins_rebinnedB);
    
    // Repeat, this time correcting pulls using the measured pull widths 
    // (from half the width of the bin_i + bin_i distribution)
    summaryunfoutfile->cd();
    MatrixUnfControl::calculateCorrelationMatrices(nVar, nbinsrhoi, nPE, kTRUE, "rebinnedA", BootstrapCorrelationMatricesFile, hAllVarUnfCorrected_rebinnedA, hist_acombfactor_rebinnedA, BinPE_rebinnedA, CoefPE_rebinnedA, pullwidthBins_rebinnedA, pullwidthVBins_rebinnedA);
    summaryunfoutfile->cd();
    MatrixUnfControl::calculateCorrelationMatrices(nVar, nbinsrhoi, nPE, kTRUE, "rebinnedB", BootstrapCorrelationMatricesFile, hAllVarUnfCorrected_rebinnedB, hist_acombfactor_rebinnedB, BinPE_rebinnedB, CoefPE_rebinnedB, pullwidthBins_rebinnedB, pullwidthVBins_rebinnedB);
    
    BootstrapCorrelationMatricesFile->Write();
    delete BootstrapCorrelationMatricesFile;

    delete hAllVarUnf_rebinnedA;
    delete hAllVarUnf_rebinnedB;

    for (int i_a = 0; i_a < nbinsrhoi/2; ++i_a)
     {
       delete hist_acombfactor_rebinnedA[i_a];
       delete hist_acombfactor_rebinnedB[i_a];
     }
  }

  printf("\n Writing summary output file ...\n");
  summaryunfoutfile->Write();
  //summaryunfoutfile->Close();
  delete summaryunfoutfile;
}


// *********************************
// FillCorrelationMatrixAllVariables
// *********************************

void MatrixUnfControl::FillCorrelationMatrixAllVariables(int nbins, TString namestring, std::vector<std::vector<std::vector<TH1D*>>>hPseudoUnfDeltapull, TH2D* &h_RhoVV, TMatrixD &m_RhoVV){

  double RhoVV[3][nbins][nbins] = {0};
  double RhoVVerr[3][nbins][nbins] = {0};

  double oneoversumoneoversigmasq=0.;
  double oneoversigmasq[2]={0.};

  h_RhoVV    = new TH2D(Form("h_Rho%s",namestring.Data()), Form("h_Rho%s",namestring.Data()), nbins, 0, nbins, nbins, 0, nbins);
  // m_RhoVV = new TMatrixD(nbins,nbins);

  // Amandeep : Filling the final matrix

  for(int i=0; i < nbins; i++){    // For the case of the full bin-by-bin matrix this is just nbinsrhoi * nVar 
    for(int j=0; j < nbins; j++){  // For the case of the full bin-by-bin matrix this is just nbinsrhoi * nVar 
      for(int l=0; l < 3; l++)   { // +/- 1 sigma and combined
        
        if(l<2){
          // sigma^2 = 2 +/- 2*rho_ij => rho_ij = (-/+)(1 - (sigma^2)/2 )
          RhoVV[l][j][i]    = (l==0?1.:-1.) * (1.-0.5*pow(hPseudoUnfDeltapull[l][j][i]->GetRMS(),2));  
          RhoVVerr[l][j][i] = hPseudoUnfDeltapull[l][j][i]->GetRMSError()*hPseudoUnfDeltapull[l][j][i]->GetRMS();
        }
        else{
          // combine + and - results by weighted average
          oneoversigmasq[0] = RhoVVerr[0][j][i] > 0 ? pow(RhoVVerr[0][j][i],-2) : 1e100;
          oneoversigmasq[1] = RhoVVerr[1][j][i] > 0 ? pow(RhoVVerr[1][j][i],-2) : 1e100;
          oneoversumoneoversigmasq = 1./(oneoversigmasq[0]+oneoversigmasq[1]);

          RhoVV[2][j][i]    = (RhoVV[0][j][i]*oneoversigmasq[0] + RhoVV[1][j][i]*oneoversigmasq[1]) * oneoversumoneoversigmasq;
          RhoVVerr[2][j][i] = sqrt(oneoversumoneoversigmasq);

          h_RhoVV->SetBinContent( j+1, i+1, RhoVV[2][j][i] );
          h_RhoVV->SetBinError( j+1, i+1, RhoVVerr[2][j][i]);
          m_RhoVV(j,i) = RhoVV[2][j][i];
        }   
      }
    }
  }
}


// ****************************
// calculateCorrelationMatrices
// ****************************

void MatrixUnfControl::calculateCorrelationMatrices(int nVar, int nbinsrhoi, int nPE, bool docorrectpullwidth, TString rebinnedstring, TFile* BootstrapCorrelationMatricesFile, TH1D* hAllVarUnf, TH1D* hist_acombfactor[], std::vector<std::vector<std::vector<float>>> &BinPE, std::vector<std::vector<std::vector<float>>> &CoefPE, std::vector<double> &pullwidthBins, std::vector<double> &pullwidthVBins) {

  const int nbinsrn = nbinsrhoi * nVar;

  std::vector<TString> mnames;
  mnames.push_back("Afb0");
  mnames.push_back("Afb1");
  mnames.push_back("Afb2");
  mnames.push_back("BinSum");
  mnames.push_back("Coef");
  mnames.push_back("Afb");

  // Amandeep : Commenting to experiment
  // TString nametag = docorrectpullwidth ? "Corrected_"+rebinnedstring : "_"+rebinnedstring;
  TString nametag    =  "_" + rebinnedstring;

  if(!docorrectpullwidth) {
    pullwidthBins.clear();
    pullwidthVBins.clear();
  }

  if(docorrectpullwidth) {
    for(int j=0; j < nbinsrn; j++) {
      for (int iPE = 0; iPE < nPE; ++iPE) {

        BinPE[j%nbinsrhoi][j/nbinsrhoi][iPE]  /= pullwidthBins[j];
        CoefPE[j%nbinsrhoi][j/nbinsrhoi][iPE] /= pullwidthVBins[j];
      }

      hAllVarUnf->SetBinError(j+1, hAllVarUnf->GetBinError(j+1)*pullwidthBins[j]);
    }
  }

  // Create and fill histograms for the bootstrap ensemble results
  std::vector<std::vector<std::vector<std::vector<TH1D*>>>> hPseudoUnfDeltapull(nbinsrhoi, std::vector<std::vector<std::vector<TH1D*>>>(2, std::vector<std::vector<TH1D*>>(nbinsrn, std::vector<TH1D*>(nbinsrn))));
  std::vector<std::vector<std::vector<TH1D*>>> hPseudoUnfDeltapullBins(2 , std::vector<std::vector<TH1D*>>(nbinsrn, std::vector<TH1D*>(nbinsrn)));
  std::vector<std::vector<std::vector<TH1D*>>> hPseudoUnfDeltapullVBins(2, std::vector<std::vector<TH1D*>>(nbinsrn, std::vector<TH1D*>(nbinsrn)));

  // Also directly calculate Pearson correlation coefficients for comparison with the 'FillCorrelationMatrixAllVariables' method. 
  // The results are basically the same, but prefer the 'FillCorrelationMatrixAllVariables' method because it also returns an estimate of the uncertainty (and also seems to be slightly more precise).
  std::vector<std::vector<std::vector<std::vector<double>>>> rPearson(nbinsrhoi, std::vector<std::vector<std::vector<double>>>(2, std::vector<std::vector<double>>(nbinsrn, std::vector<double>(nbinsrn))));
  std::vector<std::vector<std::vector<double>>> rPearsonBins(2, std::vector<std::vector<double>>(nbinsrn, std::vector<double>(nbinsrn)));
  std::vector<std::vector<std::vector<double>>> rPearsonVBins(2, std::vector<std::vector<double>>(nbinsrn, std::vector<double>(nbinsrn)));


  double nsigmas = 10.;

  for(int k=0; k < nbinsrhoi; k++){ // Over the bins
    for(int i=0; i < nVar; i++)   { // Over the variables
      for(int j=0; j < nVar; j++) { // Over the variables
        for(int l=0; l < 2; l++)  { // +/- 1 sigma

          hPseudoUnfDeltapull[k][l][j][i] = new TH1D(Form("hPseudoUnfDeltapull%s%s_%s_%s_%s",nametag.Data(), l>0?"Plus":"Minus", mnames[k].Data(), fVariableNames[j].Data(), fVariableNames[i].Data()), Form("hPseudoUnfDeltapull%s%s_%s_%s_%s",nametag.Data(), l>0?"Plus":"Minus", mnames[k].Data(), fVariableNames[j].Data(), fVariableNames[i].Data()), 200, -nsigmas, nsigmas);

          for (int iPE = 0; iPE < nPE; ++iPE)
          {
            hPseudoUnfDeltapull[k][l][j][i]->Fill(CoefPE[k][j][iPE] + (l==0?-1.:1.)*CoefPE[k][i][iPE]);
            rPearson[k][l][j][i] += (CoefPE[k][j][iPE]*(l==0?-1.:1.)*CoefPE[k][i][iPE]);
          }
          rPearson[k][l][j][i] /= nPE;
        }
      }
    }
  }

  for(int i=0; i < nbinsrn; i++){
    for(int j=0; j < nbinsrn; j++){
      for(int l=0; l < 2; l++){

        hPseudoUnfDeltapullBins[l][j][i]  = new TH1D(Form("hPseudoUnfDeltapullBins%s%s_%s%d_%s%d" ,nametag.Data(), l>0?"Plus":"Minus", fVariableNames[j/nbinsrhoi].Data(), j%nbinsrhoi, fVariableNames[i/nbinsrhoi].Data(), i%nbinsrhoi), Form("hPseudoUnfDeltapullBins%s%s_%s%d_%s%d" , nametag.Data(), l>0?"Plus":"Minus", fVariableNames[j/nbinsrhoi].Data(), j%nbinsrhoi, fVariableNames[i/nbinsrhoi].Data(),i%nbinsrhoi), 200, -nsigmas, nsigmas);
        hPseudoUnfDeltapullVBins[l][j][i] = new TH1D(Form("hPseudoUnfDeltapullVBins%s%s_%s%d_%s%d",nametag.Data(), l>0?"Plus":"Minus", fVariableNames[j/nbinsrhoi].Data(), j%nbinsrhoi, fVariableNames[i/nbinsrhoi].Data(), i%nbinsrhoi), Form("hPseudoUnfDeltapullVBins%s%s_%s%d_%s%d", nametag.Data(), l>0?"Plus":"Minus", fVariableNames[j/nbinsrhoi].Data(), j%nbinsrhoi, fVariableNames[i/nbinsrhoi].Data(),i%nbinsrhoi), 200, -nsigmas, nsigmas);

        for (int iPE = 0; iPE < nPE; ++iPE)
        {
          hPseudoUnfDeltapullBins[l][j][i]->Fill(BinPE[j%nbinsrhoi][j/nbinsrhoi][iPE]   + (l==0?-1.:1.)*BinPE[i%nbinsrhoi][i/nbinsrhoi][iPE]);
          hPseudoUnfDeltapullVBins[l][j][i]->Fill(CoefPE[j%nbinsrhoi][j/nbinsrhoi][iPE] + (l==0?-1.:1.)*CoefPE[i%nbinsrhoi][i/nbinsrhoi][iPE]);
          rPearsonBins[l][j][i] += (BinPE[j%nbinsrhoi][j/nbinsrhoi][iPE] * (l==0?-1.:1.) * BinPE[i%nbinsrhoi][i/nbinsrhoi][iPE]);
          rPearsonVBins[l][j][i]+= (CoefPE[j%nbinsrhoi][j/nbinsrhoi][iPE]* (l==0?-1.:1.) * CoefPE[i%nbinsrhoi][i/nbinsrhoi][iPE]);
        }
        rPearsonBins[l][j][i] /= nPE;
        rPearsonVBins[l][j][i]/= nPE;

        // Amandeep : adding for experimentation
        std::cout << "rPearsonBins :: "  << rPearsonBins[l][j][i]  << std::endl;
        std::cout << "rPearsonVBins :: " << rPearsonVBins[l][j][i] << std::endl;
        // End
      }
    }

    if(!docorrectpullwidth) {
      pullwidthBins.push_back(hPseudoUnfDeltapullBins[1][i][i]->GetRMS()/2.);
      pullwidthVBins.push_back(hPseudoUnfDeltapullVBins[1][i][i]->GetRMS()/2.);
    }
  }

  BootstrapCorrelationMatricesFile->cd();
  
  TH2D* h_RhoVV[nbinsrhoi];

  for (int k = 0; k < nbinsrhoi; ++k)
  {
    TMatrixD m_RhoVV(nVar,nVar);
    MatrixUnfControl::FillCorrelationMatrixAllVariables(nVar, Form("VV%s_%s",nametag.Data(), mnames[k].Data()), hPseudoUnfDeltapull.at(k), h_RhoVV[k], m_RhoVV);
  }
  std::cout << "Pearson (above)" << std::endl;
  std::cout << "DeltaPullWidth (below):" << std::endl;
  for(int i=0; i < nVar ;i++){
    for(int j=0; j < nVar; j++){
      std::cout << i+1 << " " << j+1 << " " << h_RhoVV[4]->GetBinContent(j+1,i+1)<< std::endl;
    }
  }

  TH2D* h_RhoBB;
  TMatrixD m_RhoBB(nbinsrn,nbinsrn);
  MatrixUnfControl::FillCorrelationMatrixAllVariables(nbinsrn, Form("BB%s",nametag.Data()), hPseudoUnfDeltapullBins, h_RhoBB, m_RhoBB);

  TH2D* h_RhoVBVB;
  TMatrixD m_RhoVBVB(nbinsrn,nbinsrn);
  MatrixUnfControl::FillCorrelationMatrixAllVariables(nbinsrn, Form("VBVB%s",nametag.Data()), hPseudoUnfDeltapullVBins, h_RhoVBVB, m_RhoVBVB);

  // Cross-check that we get the same (to first order) correlation matrix for each of the measured variables 
  // when we calculate it starting from the correlation matrix for all the measured bins, m_RhoBB

  TMatrixD m_AFB(nVar,nVar);
  std::vector<double> afbvar;
  std::vector<double> afberrvar;
  MatrixUnf::GetAfbCorrMAllVars(hAllVarUnf, m_RhoBB, afbvar, afberrvar, m_AFB);

  TH2D* h_RhoAA = new TH2D(Form("h_RhoAA%s",nametag.Data()), Form("h_RhoAA%s",nametag.Data()), nVar, 0, nVar, nVar, 0, nVar);
  for(int i=0; i < nVar; i++){
    for(int j=0; j < nVar; j++){
      h_RhoAA->SetBinContent( j+1, i+1, m_AFB(j,i) );
    }
  }

  TMatrixD m_BinSum(nVar,nVar);
  std::vector<double> binsum;
  std::vector<double> binsumerr;
  MatrixUnf::GetBinSumCorrMAllVars(hAllVarUnf, m_RhoBB, binsum, binsumerr, m_BinSum);

  TH2D* h_RhoBSBS = new TH2D(Form("h_RhoBSBS%s",nametag.Data()), Form("h_RhoBSBS%s",nametag.Data()), nVar, 0, nVar, nVar, 0, nVar);
  for(int i=0; i < nVar; i++){
    for(int j=0; j < nVar; j++){
      h_RhoBSBS->SetBinContent( j+1, i+1, m_BinSum(j,i) );
    }
  }

  TMatrixD m_AFBB(nVar*nbinsrhoi/2,nVar*nbinsrhoi/2);
  TMatrixD m_C(nVar,nVar);
  std::vector<std::vector<double>> afbvecvec;
  std::vector<std::vector<double>> afberrvecvec;
  MatrixUnf::GetAfbBinsCorrMAllVars(hAllVarUnf, m_RhoBB, afbvecvec, afberrvecvec, m_AFBB, hist_acombfactor, m_C);

  TH2D* h_RhoCC = new TH2D(Form("h_RhoCC%s",nametag.Data()), Form("h_RhoCC%s",nametag.Data()), nVar, 0, nVar, nVar, 0, nVar);
  for(int i=0; i < nVar; i++){
    for(int j=0; j < nVar; j++){
      h_RhoCC->SetBinContent( j+1, i+1, m_C(j,i) );
    }
  }

  TH2D* h_RhoABAB = new TH2D(Form("h_RhoABAB%s",nametag.Data()), Form("h_RhoABAB%s",nametag.Data()), nVar*nbinsrhoi/2, 0, nVar*nbinsrhoi/2, nVar*nbinsrhoi/2, 0, nVar*nbinsrhoi/2);
  for(int i=0; i < nVar*nbinsrhoi/2; i++){
    for(int j=0; j < nVar*nbinsrhoi/2; j++){
      h_RhoABAB->SetBinContent( j+1, i+1, m_AFBB(j,i) );
    }
  }

  TH2D* h_RhoAiAi[nbinsrhoi/2];
  for(int k=0;k<nbinsrhoi/2;k++){
    h_RhoAiAi[k] = new TH2D(Form("h_RhoA%dA%d%s",k,k,nametag.Data()), Form("h_RhoA%dA%d%s",k,k,nametag.Data()), nVar, 0, nVar, nVar, 0, nVar);
    for(int i=0; i < nVar; i++){
      for(int j=0; j < nVar; j++){
        h_RhoAiAi[k]->SetBinContent( j+1, i+1, m_AFBB(j*nbinsrhoi/2+k,i*nbinsrhoi/2+k) );
      }
    }
  }
}


// **********
// unfoldData
// **********

Int_t MatrixUnfControl::unfoldData(MatrixUnf* pttMatUnf, TString NameTUnfObjX,TH1D* input, TH2D* MigMat, TH1D* gen, TH1D* visgen, TH1D* reco, std::vector<TH1D*> vecbg, Double_t tauMin, Double_t tauMax, Bool_t ensTest, TString Chan, std::vector<TH2D*>vecMigMatSys) 
{
  //Modify this later based on whether the MigMat is filled correctly
  /*
  input = (TH1D*) MigMat->ProjectionY();
  reco  = (TH1D*) MigMat->ProjectionY();
  gen   = (TH1D*) MigMat->ProjectionX();

  TH1D* bg = (TH1D*) vecbg.at(0)->Clone();
  for (int i = 1; i < vecbg.size(); ++i)
  {
    bg->Add(vecbg.at(i));
  }
  */

  double Afb, AfbErr;
  MatrixUnf::GetAfb(input, Afb, AfbErr);
  std::cout << "Input AFB: " << Afb << " +/- " << AfbErr << std::endl;

  pttMatUnf->SetName(NameTUnfObjX);
  // Fill the input histos
  if(!ensTest) pttMatUnf->SetInputDists(gen, visgen, reco, vecbg, reco, MigMat, vecMigMatSys);
  if(!ensTest) pttMatUnf->createInputResults(ensTest); 
  pttMatUnf->prepRunUnfolding(input, fRegMode, 0, tauMin, tauMax, fLumi, fDoSubtrBg, fUseAreaConstrain, ensTest, fMinRhoAVG, fPseudoUnfResults, fPseudoUnfResultsRebinnedA, fPseudoUnfResultsRebinnedB, fPseudoUnfAfbResultsRebinnedA, fPseudoUnfAfbResultsRebinnedB, fPseudoUnfOtherResultsRebinnedA, fPseudoUnfOtherResultsRebinnedB);
  if(!ensTest) pttMatUnf->AddObjectsToFile();
  
  return 1;
}


// ***************
// FillHistoArrays 
// ***************

// Get the histograms stored in a file and store the pointer in a TObjArray for further use
// also pass a string here so that that string can be appended in front of the name

bool MatrixUnfControl::FillHistoArrays(TString CHAN){

  TString bgfracplotsfile = "BGfracplots/"+CHAN+"/";
  gSystem->mkdir(bgfracplotsfile.Data(), true);

  int nbgfiles=0;
  if(fBgFiles) nbgfiles = fBgFiles->GetEntries();
  int ndatafiles=0;
  if(fDataFiles) ndatafiles = fDataFiles->GetEntries();
  int nMCfiles=0;
  // Amandeep : Here the spirit of the variable is count the number of signal files 
  if(fMCFiles) nMCfiles = fMCFiles->GetEntries();
  int nsysfiles=0;
  if(fSysFiles) {
    nsysfiles = fSysFiles->GetEntries();
    // Amandeep : Maybe this needs to be modified for 2 files ?
    if( nsysfiles % (nMCfiles>0?nMCfiles:1) != 0 ) std::cout<<"MatrixUnfControl::FillHistoArrays: Error, didnt find the same number of systematics for each channel."<<std::endl;
    nsysfiles /= (nMCfiles>0?nMCfiles:1);
  }
 
  if(fDataArray==NULL){
    fDataArray = new TObjArray();
  }
  else{
    fDataArray->Clear();
  }
  
  if(fRecArray==NULL){
    fRecArray = new TObjArray();
  }
  else{
    fRecArray->Clear();
  }
  
  if(fGenArray==NULL){
    fGenArray = new TObjArray();
  }
  else{
    fGenArray->Clear();
  }

  if(fVisGenArray==NULL){
    fVisGenArray = new TObjArray();
  }
  else{
    fVisGenArray->Clear();
  }
  
  if(fMigMatArray==NULL){
    fMigMatArray = new TObjArray();
  }
  else{
    fMigMatArray->Clear();
  }
  
  if(fBgArray==NULL){
    fBgArray = new TObjArray();
  }
  else{
    fBgArray->Clear();
  }
  if(fBgSamplesArray==NULL){
    fBgSamplesArray = new TObjArray();
  }
  else{
    fBgSamplesArray->Clear();
  }

  if(fMigMatSysArray==NULL){
    fMigMatSysArray = new TObjArray();
  }
  else{
    fMigMatSysArray->Clear();
  }

  if(fTrueDataArray==NULL){
    fTrueDataArray = new TObjArray();
  }
  else{
    fTrueDataArray->Clear();
  }
  
  std::cout << "fVariableNames.size() :: " << fVariableNames.size() << std::endl;
  
  // Amandeep : Over all variables
  for(unsigned int i=0; i < fVariableNames.size(); i++)
  {
 
    std::cout << "Currently working with variable :: " << fVariableNames[i].Data() << std::endl;

    // ****
    // DATA
    // ****

    if(fDataFiles == NULL && !fSelfConsistency) {
      // Require 2 root files for 2015 data and more for 2016
      std::cout<<"require 2 root files for 2015 and more for 2016"<<std::endl;
      return kFALSE;
    }
    if(fMCFiles == NULL){
      std::cout<<"No MC file(s) provided"<<std::endl;
      return kFALSE;
    }

    // No requirement in number of background files because should have option to not use any background at all
    
    std::vector<TString>channels;
    channels.push_back("ee");
    channels.push_back("emu");
    channels.push_back("mumu");
    channels.push_back("combined");

    // Amandeep : kMaxNumDataFiles is set to 100 in the .h file
    TH1D* htempdata[kMaxNumDataFiles+1];
    TH1D* hdatabgsub_channels[4];

    double data_yields[4]       = {0,0,0,0};
    double bkgsubdata_yields[4] = {0,0,0,0};
    
    if(fDataFiles) 
    {
      for (int j = 0; j < 4; ++j)
      {
        hdatabgsub_channels[j]= (TH1D*)((TH1D*)((TFile*)fDataFiles->At(0))->Get(Form("hreco_%s",fVariableNames[i].Data())))->Clone(Form("%sDataBgSub_%s",fVariableNames[i].Data(),channels.at(j).Data()));
        hdatabgsub_channels[j]->Reset();
      }
      htempdata[0]= (TH1D*)((TH1D*)((TFile*)fDataFiles->At(0))->Get(Form("hreco_%s",fVariableNames[i].Data())))->Clone(Form("%sData",fVariableNames[i].Data()));
      htempdata[0]->Reset();

      for(int j=0; j < ndatafiles; j++){
      
        htempdata[j+1]= (TH1D*) ((TFile*)fDataFiles->At(j))->Get(Form("hreco_%s",fVariableNames[i].Data()));
        std::cout << fDataFiles->At(j)->GetName() << " Integral: " << htempdata[j+1]->Integral() << std::endl;
        if( TString(fDataFiles->At(j)->GetName()).Contains("histosTUnfold_ee_") )   {hdatabgsub_channels[0]->Add(htempdata[j+1]); data_yields[0]+=htempdata[j+1]->Integral();}
        if( TString(fDataFiles->At(j)->GetName()).Contains("histosTUnfold_emu_") )  {hdatabgsub_channels[1]->Add(htempdata[j+1]); data_yields[1]+=htempdata[j+1]->Integral();}
        if( TString(fDataFiles->At(j)->GetName()).Contains("histosTUnfold_mumu_") ) {hdatabgsub_channels[2]->Add(htempdata[j+1]); data_yields[2]+=htempdata[j+1]->Integral();}
        hdatabgsub_channels[3]->Add(htempdata[j+1]);
        data_yields[3]+=htempdata[j+1]->Integral();

        htempdata[0]->Add(htempdata[j+1]);
      }
      std::cout << "data yields ee emu mumu combined: " << data_yields[0] << " "<< data_yields[1] << " " << data_yields[2]<<" "<<data_yields[3]<<std::endl;
      fDataArray->AddLast(htempdata[0]);
    }


    // ***************
    // MC signal files
    // ***************

    // Reco
    TH1D* htempreco[nMCfiles+1];
    // Gen
    TH1D* htempgen[nMCfiles+1];
    // VisGen
    TH1D* htempvisgen[nMCfiles+1];
    // MigMat  
    TH2D* htemprecoVsgen[nMCfiles+1];
    
    TH1D* hrecoplustau_channels[4];
    double sig_yields[4] = {0,0,0,0};
    
    if(fMCFiles) 
    {

      htempreco[0]  = (TH1D*)((TH1D*)((TFile*)fMCFiles->At(0))->Get(Form("hreco_%s",fVariableNames[i].Data())))->Clone(Form("%sReco",fVariableNames[i].Data()));
      htempreco[0]->Reset();
      // std::cout<<fMCFiles->At(0)->GetName()<<" Integral: "<<htempreco[0]->Integral()<<" xsec: "<<fXSectionMC.at(0)<<std::endl;

      for (int j = 0; j < 4; ++j)
      {
        hrecoplustau_channels[j]= (TH1D*)((TH1D*)((TFile*)fMCFiles->At(0))->Get(Form("hreco_%s",fVariableNames[i].Data())))->Clone(Form("%sReco_%s",fVariableNames[i].Data(),channels.at(j).Data()));
        hrecoplustau_channels[j]->Reset();
      }

      htempgen[0]    = (TH1D*)((TH1D*)((TFile*)fMCFiles->At(0))->Get(Form("hgen_%s",fVariableNames[i].Data())))->Clone(Form("%sGen",fVariableNames[i].Data()));
      htempgen[0]->Reset();

      htempvisgen[0] = (TH1D*)((TH1D*)((TFile*)fMCFiles->At(0))->Get(Form("hvisgen_%s",fVariableNames[i].Data())))->Clone(Form("%sVisGen",fVariableNames[i].Data()));
      htempvisgen[0]->Reset();

      for (int j=0; j < nMCfiles; j++)
      {
        TString MCFileName = fMCFiles->At(j)->GetName();

        htempreco[j+1] = (TH1D*)((TFile*)fMCFiles->At(j))->Get(Form("hreco_%s",fVariableNames[i].Data()));
        // Amandeep : Commenting to experiment with external lumi scaling 
        //htempreco[j+1]->Scale(getLumiWeight( (TFile*) fMCFiles->At(j), fXSectionMC.at(j), fLumiMC.at(j)));
        // End
        std::cout << MCFileName << " Integral: " << htempreco[j+1]->Integral() << " xsec: " << fXSectionMC.at(j) << std::endl;
        
        if( MCFileName.Contains("histosTUnfold_ee_") )   {hrecoplustau_channels[0]->Add(htempreco[j+1]); sig_yields[0]+=htempreco[j+1]->Integral();}
        if( MCFileName.Contains("histosTUnfold_emu_") )  {hrecoplustau_channels[1]->Add(htempreco[j+1]); sig_yields[1]+=htempreco[j+1]->Integral();}
        if( MCFileName.Contains("histosTUnfold_mumu_") ) {hrecoplustau_channels[2]->Add(htempreco[j+1]); sig_yields[2]+=htempreco[j+1]->Integral();}
        
        hrecoplustau_channels[3]->Add(htempreco[j+1]);
        sig_yields[3]+=htempreco[j+1]->Integral();

        htempreco[0]->Add(htempreco[j+1]);

        htempgen[j+1]    = (TH1D*)((TFile*)fMCFiles->At(j))->Get(Form("hgen_%s",fVariableNames[i].Data()));
        // Amandeep : Commenting to experiment with external lumi scaling 
        //htempgen[j+1]->Scale(getLumiWeight( (TFile*) fMCFiles->At(j), fXSectionMC.at(j), fLumiMC.at(j)));
        // End
        htempgen[0]->Add(htempgen[j+1]);

        htempvisgen[j+1] = (TH1D*)((TFile*)fMCFiles->At(j))->Get(Form("hvisgen_%s",fVariableNames[i].Data()));
        // Amandeep : Commenting to experiment with external lumi scaling  
        //htempvisgen[j+1]->Scale(getLumiWeight( (TFile*) fMCFiles->At(j), fXSectionMC.at(j), fLumiMC.at(j)));
        // End
        htempvisgen[0]->Add(htempvisgen[j+1]);
      }

      std::cout << "sig yields ee emu mumu combined: " << sig_yields[0] << " " << sig_yields[1] << " " << sig_yields[2] << " " << sig_yields[3] << std::endl;

      fRecArray->AddLast(htempreco[0]);
      fGenArray->AddLast(htempgen[0]);
      fVisGenArray->AddLast(htempvisgen[0]);
    }

    // ***********
    // BACKGROUNDS
    // ***********

    // Amandeep : Maybe this is even right ?
    // kMaxNumBgFiles is defined in the .h to be 100 

    TH1D* htempbg[kMaxNumBgFiles + 1];
    TH1D* htaubg_channels[4];
    TH1D* hMCtaubg_channels[4];
    TH1D* htWfrac_channels[4];
    TH1D* htaufrac_channels[4];

    double bkg_yields[4] = {0,0,0,0};
    double tau_yields[4] = {0,0,0,0};
    
    // Amandeep : commenting since we don't have taus as bgs
    // double tau_SF[4]={0,0,0,0};
    // End

    double channel_SF[4] = {1,1,1,1};
    double taufrac_SF[4] = {1,1,1,1};

    if(fBgFiles) {
      htempbg[0]= (TH1D*)((TH1D*)((TFile*)fBgFiles->At(0))->Get(Form("hreco_%s",fVariableNames[i].Data())))->Clone(Form("%sRecoBg",fVariableNames[i].Data()));
      htempbg[0]->Reset();
      //std::cout<<fBgFiles->At(0)->GetName()<<" Integral: "<<htempbg[0]->Integral()<<" xsec: "<<fXSectionBg.at(0)<<std::endl;

      for (int j = 0; j < 4; ++j)
      {
        htWfrac_channels[j]= (TH1D*)((TH1D*)((TFile*)fBgFiles->At(0))->Get(Form("hreco_%s",fVariableNames[i].Data())))->Clone(Form("%stWfrac_%s",fVariableNames[i].Data(),channels.at(j).Data()));
        htWfrac_channels[j]->Reset();
      }
      
      std::vector<TString> taubgname;
      std::vector<int> itau;
      //std::vector<TH1D*> tauhist;

      for (int j=0; j < nbgfiles; j++)
      {
        TString BgFileName = fBgFiles->At(j)->GetName();

        htempbg[j+1]       = (TH1D*)((TH1D*)((TFile*) fBgFiles->At(j))->Get(Form("hreco_%s",fVariableNames[i].Data())))->Clone(Form("%sRecoBg_%s",fVariableNames[i].Data(),BgFileName.Data()));
        // Amandeep : Commenting to experiment with external lumi scaling 
        //double lumiweight  = getLumiWeight(((TFile*)fBgFiles->At(j)), fXSectionBg[j],  fLumiBg[j]);
        //htempbg[j+1]->Scale(lumiweight);
        // End

        double bkgintegral, bkgerror, bkgnentries, bkgneffentries;
        bkgintegral    = htempbg[j+1]->IntegralAndError(1, htempbg[j+1]->GetNbinsX(), bkgerror);
        bkgnentries    = htempbg[j+1]->GetEntries();
        bkgneffentries = htempbg[j+1]->GetEffectiveEntries();

        // std::cout<< BgFileName << " Integral: " << bkgintegral << " +/- " << bkgerror << " xsec: " << fXSectionBg.at(j) << " lumiweight: " << lumiweight << " bkgnentries: " << bkgnentries << " bkgneffentries: " << bkgneffentries << std::endl;
        
        if( BgFileName.Contains("histosTUnfold_ee_") )   {hdatabgsub_channels[0]->Add(htempbg[j+1],-1.); bkg_yields[0] += bkgintegral;}
        if( BgFileName.Contains("histosTUnfold_emu_") )  {hdatabgsub_channels[1]->Add(htempbg[j+1],-1.); bkg_yields[1] += bkgintegral;}
        if( BgFileName.Contains("histosTUnfold_mumu_") ) {hdatabgsub_channels[2]->Add(htempbg[j+1],-1.); bkg_yields[2] += bkgintegral;}
        
        hdatabgsub_channels[3]->Add(htempbg[j+1],-1.);
        bkg_yields[3] += bkgintegral;
        htempbg[0]->Add(htempbg[j+1]);

        if( BgFileName.Contains("single") ) {
          if( BgFileName.Contains("histosTUnfold_ee_") )   {htWfrac_channels[0]->Add(htempbg[j+1]);}
          if( BgFileName.Contains("histosTUnfold_emu_") )  {htWfrac_channels[1]->Add(htempbg[j+1]);}
          if( BgFileName.Contains("histosTUnfold_mumu_") ) {htWfrac_channels[2]->Add(htempbg[j+1]);}
          htWfrac_channels[3]->Add(htempbg[j+1]);
        }

      }

      for (int j = 0; j < 4; ++j)
      {
        for (int i_bin = 1; i_bin <= hdatabgsub_channels[j]->GetNbinsX(); ++i_bin)
        {
          if(hdatabgsub_channels[j]->GetBinContent( i_bin )<0) {
            std::cout<<"belowzerobin: "<<i_bin<<" "<<hdatabgsub_channels[j]->GetBinContent( i_bin )<<std::endl;
            hdatabgsub_channels[j]->SetBinContent( i_bin, 0 );
          }
        }
      }

      // Amandeep : changing to kFALSE to deal with tau stuff 
      // Original : 
      // bool doDataFractionalTauSub = kTRUE;
      // bool useCombinedTauFrac     = kTRUE;

      // Modified :
      bool doDataFractionalTauSub = kFALSE;
      bool useCombinedTauFrac     = kFALSE;
      // End

      // Commenting the block below and this is effectively deadcode now
      if(doDataFractionalTauSub) {

        // std::cout<<"tau yields ee emu mumu combined: "<<tau_yields[0]<<" "<<tau_yields[1]<<" "<<tau_yields[2]<<" "<<tau_yields[3]<<std::endl;
        // std::cout<<"bkg yields ee emu mumu combined: "<<bkg_yields[0]<<" "<<bkg_yields[1]<<" "<<bkg_yields[2]<<" "<<bkg_yields[3]<<std::endl;

        // for (int j = 3; j >= 0; --j)
        // {
        //   hMCtaubg_channels[j]=(TH1D*) htaubg_channels[j]->Clone(Form("%s_MCtaubg_%s",fVariableNames[i].Data(),channels.at(j).Data()));
        //   hrecoplustau_channels[j]->Add(hMCtaubg_channels[j]);
        //   std::cout<<"SimpleIntegral: "<<hdatabgsub_channels[j]->Integral()*hMCtaubg_channels[j]->Integral()/hrecoplustau_channels[j]->Integral()<<std::endl;

        //   if(useCombinedTauFrac) taufrac_SF[j] =  hMCtaubg_channels[j]->Integral() == 0 ? 0 : hMCtaubg_channels[j]->Integral()/hrecoplustau_channels[j]->Integral();

        //   htaubg_channels[j]->Divide(hrecoplustau_channels[j]);
        //   //htaubg_channels[j]->Print("all");

        //   htaufrac_channels[j]=(TH1D*) htaubg_channels[j]->Clone(Form("%s_TauFraction_%s",fVariableNames[i].Data(),channels.at(j).Data()));

        //   if(useCombinedTauFrac && j!=3) htaubg_channels[j] = (TH1D*) htaufrac_channels[3]->Clone(htaubg_channels[j]->GetName());
        //   htaubg_channels[j]->Multiply(hdatabgsub_channels[j]);
        //   //htaubg_channels[j]->Print("all");

        //   htWfrac_channels[j]->Divide(hrecoplustau_channels[j]);
        // }


        // if(useCombinedTauFrac) {

        //   for (int j = 0; j < 4; ++j) {
        //     taufrac_SF[j]/=taufrac_SF[3];
        //   }

        //   //check the tau fractions are compatible between channels
        //   TCanvas *ctaufrac = new TCanvas("ctaufrac", "ctaufrac");
        //   for (int j = 3; j >= 0; --j)
        //   {

        //     htaufrac_channels[j]->SetLineColor(4-j);
        //     if(j==3) htaufrac_channels[j]->Draw();
        //     else {
        //       htaufrac_channels[j]->Draw("same");
        //       htaufrac_channels[j]->Chi2Test(htaufrac_channels[3],"WW CHI2 P");
        //       htaufrac_channels[j]->KolmogorovTest(htaufrac_channels[3],"D");
        //     }


        //   }
        //   ctaufrac->SaveAs(Form("%staufrac_%s.pdf",bgfracplotsfile.Data(),fVariableNames[i].Data()));
        //   for (int j = 3; j >= 0; --j)
        //   {

        //     htWfrac_channels[j]->SetLineColor(4-j);
        //     if(j==3) htWfrac_channels[j]->Draw();
        //     else {
        //       htWfrac_channels[j]->Draw("same");
        //       htWfrac_channels[j]->Chi2Test(htWfrac_channels[3],"WW CHI2 P");
        //       htWfrac_channels[j]->KolmogorovTest(htWfrac_channels[3],"D");
        //     }


        //   }
        //   gStyle->SetOptStat(0);
        //   ctaufrac->SaveAs(Form("%stWfrac_%s.pdf",bgfracplotsfile.Data(),fVariableNames[i].Data()));

        //   htaufrac_channels[3]->GetYaxis()->SetRangeUser(0,0.16);
        //   htaufrac_channels[3]->GetYaxis()->SetTitle("Fraction of events");
        //   htaufrac_channels[3]->Draw();
        //   htWfrac_channels[3]->SetLineColor(2);
        //   htWfrac_channels[3]->Draw("same");
        //   gStyle->SetOptStat(0);
        //   ctaufrac->SaveAs(Form("%stWtaufrac_%s.pdf",bgfracplotsfile.Data(),fVariableNames[i].Data()));
        //   gStyle->SetOptStat(1111);

        //   TCanvas *ctaufracfrac = new TCanvas("ctaufracfrac", "ctaufracfrac");
        //   gPad->SetLogy();
        //   for (int j = 2; j >= 0; --j)
        //   {

        //     htaufrac_channels[j]->Divide(htaufrac_channels[3]);
        //     htaufrac_channels[j]->SetMaximum(2.);
        //     htaufrac_channels[j]->SetMinimum(0.5);
        //     if(j==2) htaufrac_channels[j]->Draw();
        //     else {
        //       htaufrac_channels[j]->Draw("same");
        //     }
        //     htaufrac_channels[j]->Fit("pol0");
        //     //if(htaufrac_channels[j]->GetFunction("pol0")) taufrac_SF[j] = htaufrac_channels[j]->GetFunction("pol0")->GetParameter(0);

        //   }
        //   ctaufracfrac->SaveAs(Form("%staufracfrac_%s.pdf",bgfracplotsfile.Data(),fVariableNames[i].Data()));
        //   std::cout<<"taufrac_SF: "<<taufrac_SF[0]<<" "<<taufrac_SF[1]<<" "<<taufrac_SF[2]<<" "<<taufrac_SF[3]<<std::endl;
        // }


        // for (int j = 3; j >= 0; --j)
        // {
        //   if(useCombinedTauFrac && j!=3) htaubg_channels[j]->Scale( taufrac_SF[j] );

        //   //set uncertainties based on the MC stat uncertainties, because the tau MC stat uncertainties give approximately the stat uncertainty on the tau fraction (temporary solution : it would be better to let the tau subtraction vary as a systematic)
        //   for (int i_bin = 1; i_bin <= htaubg_channels[j]->GetNbinsX(); ++i_bin)
        //   {
        //     double binSF = htaubg_channels[j]->GetBinContent(i_bin) / hMCtaubg_channels[j]->GetBinContent(i_bin);
        //     if( !(binSF>0.7 && binSF<1/0.7) ) {
        //       std::cout<<"Warning, inconsistent data/MC ratio for ttbarinctau. binSF: "<<i_bin<<" "<<binSF<<" "<<hMCtaubg_channels[j]->GetBinError(i_bin)<<", new binSF: ";
        //       binSF = binSF<1?0.7:1/0.7;
        //       std::cout<<binSF<<std::endl;
        //     }
        //     htaubg_channels[j]->SetBinError( i_bin, hMCtaubg_channels[j]->GetBinError(i_bin) * binSF );
        //   }


        //   tau_yields[j]=htaubg_channels[j]->Integral();
        //   bkg_yields[j]+=tau_yields[j];
        //   bkgsubdata_yields[j] = data_yields[j] - bkg_yields[j];
        //   hdatabgsub_channels[j]->Add(htaubg_channels[j],-1.);
        // }


        // std::cout<<"final tau yields ee emu mumu combined: "<<tau_yields[0]<<" "<<tau_yields[1]<<" "<<tau_yields[2]<<" "<<tau_yields[3]<<std::endl;
        // std::cout<<"final bkg yields ee emu mumu combined: "<<bkg_yields[0]<<" "<<bkg_yields[1]<<" "<<bkg_yields[2]<<" "<<bkg_yields[3]<<std::endl;
        // std::cout<<"bkgsubdata yields ee emu mumu combined: "<<bkgsubdata_yields[0]<<" "<<bkgsubdata_yields[1]<<" "<<bkgsubdata_yields[2]<<" "<<bkgsubdata_yields[3]<<std::endl;
        // std::cout<<"bkgsubdata yields ee emu mumu combined: "<<hdatabgsub_channels[0]->Integral()<<" "<<hdatabgsub_channels[1]->Integral()<<" "<<hdatabgsub_channels[2]->Integral()<<" "<<hdatabgsub_channels[3]->Integral()<<std::endl;


        // for (unsigned int i_file = 0; i_file < itau.size(); ++i_file)
        // {
        //   int i_chan = 3;
        
        //   if( taubgname.at(i_file).Contains("histosTUnfold_ee_") ) i_chan = 0;
        //   if( taubgname.at(i_file).Contains("histosTUnfold_emu_") ) i_chan = 1;
        //   if( taubgname.at(i_file).Contains("histosTUnfold_mumu_") ) i_chan = 2;

        //   htempbg[ itau.at(i_file) ] = htaubg_channels[i_chan];

        //   htempbg[0]->Add( htaubg_channels[i_chan] );

        // }
      }
      // End

      else{

        // Amandeep : Commenting since we don't hav tau bgs
        // for (int i_chan = 0; i_chan < 4; ++i_chan)
        // {
        //   if( sig_yields[i_chan] > 0 ) tau_SF[i_chan] = ( data_yields[i_chan] - bkg_yields[i_chan] ) / ( sig_yields[i_chan] + tau_yields[i_chan] );
        // }
        // End

        std::cout<<"bkg yields ee emu mumu combined: "<<bkg_yields[0]<<" "<<bkg_yields[1]<<" "<<bkg_yields[2]<<" "<<bkg_yields[3]<<std::endl;
        std::cout<<"tau yields ee emu mumu combined: "<<tau_yields[0]<<" "<<tau_yields[1]<<" "<<tau_yields[2]<<" "<<tau_yields[3]<<std::endl;
        // std::cout<<"tau SF ee emu mumu combined: "<<tau_SF[0]<<" "<<tau_SF[1]<<" "<<tau_SF[2]<<" "<<tau_SF[3]<<std::endl;


        for (unsigned int i_file = 0; i_file < itau.size(); ++i_file)
        {
          // Amandeep : Revisit this logic
          int i_chan = 3;
        
          if( taubgname.at(i_file).Contains("histosTUnfold_ee_") )   i_chan = 0;
          if( taubgname.at(i_file).Contains("histosTUnfold_emu_") )  i_chan = 1;
          if( taubgname.at(i_file).Contains("histosTUnfold_mumu_") ) i_chan = 2;

          // std::cout<<"Scaling "<< taubgname.at(i_file) <<" by "<<tau_SF[i_chan]<<" :"<<std::endl;
          // htempbg[ itau.at(i_file) ]->Scale( tau_SF[i_chan] );

          std::cout<< taubgname.at(i_file) << " Integral: " << htempbg[ itau.at(i_file) ]->Integral() << " itau: " << itau.at(i_file) << std::endl;
          bkg_yields[i_chan] += htempbg[ itau.at(i_file) ]->Integral();

          if(i_chan!=3) bkg_yields[3]+=htempbg[ itau.at(i_file) ]->Integral();
          htempbg[0]->Add(htempbg[ itau.at(i_file) ]);

          // htempbg[ itau.at(i_file) ]->Print("all");
          std::cout << "Integral: "<< htempbg[ itau.at(i_file)]->Integral() << std::endl;
        }


        for (int i_chan = 0; i_chan < 4; ++i_chan)
        {
          // tau_yields[i_chan]    *= tau_SF[i_chan];
          bkgsubdata_yields[i_chan] = data_yields[i_chan] - bkg_yields[i_chan];
        }

        std::cout<<"final tau yields ee emu mumu combined: "<<tau_yields[0]<<" "<<tau_yields[1]<<" "<<tau_yields[2]<<" "<<tau_yields[3]<<std::endl;
        std::cout<<"final bkg yields ee emu mumu combined: "<<bkg_yields[0]<<" "<<bkg_yields[1]<<" "<<bkg_yields[2]<<" "<<bkg_yields[3]<<std::endl;
        std::cout<<"bkgsubdata yields ee emu mumu combined: "<<bkgsubdata_yields[0]<<" "<<bkgsubdata_yields[1]<<" "<<bkgsubdata_yields[2]<<" "<<bkgsubdata_yields[3]<<std::endl;
      }

      // Amandeep : nMCfiles was hardcoded here | Replace with channel.compare('combined') ?
      if( CHAN.CompareTo("combined") == 0 ) 
      { 
        CalculateChannelFractions(bkgsubdata_yields, sig_yields, channel_SF);
      }

      fBgArray->AddLast(htempbg[0] );

      for(int j=0;j<nbgfiles;j++)
      {
        fBgSamplesArray->AddLast(htempbg[j+1] );
      }
    }

    bool doCorrectChannelFractions = kTRUE;

    if(fMCFiles) 
    {
      htemprecoVsgen[0] = (TH2D*)((TH2D*)((TFile*)fMCFiles->At(0))->Get(Form("hrecoVsgen_%s",fVariableNames[i].Data())))->Clone(Form("%sRespMat",fVariableNames[i].Data()));
      htemprecoVsgen[0]->Reset();

      for(int j=0; j < nMCfiles; j++)
      {
        TString MCFileName  = fMCFiles->At(j)->GetName();

        htemprecoVsgen[j+1] = (TH2D*)((TFile*)fMCFiles->At(j))->Get(Form("hrecoVsgen_%s",fVariableNames[i].Data()));    
        // Amandeep : Commenting to experiment with external lumi scaling 
        //htemprecoVsgen[j+1]->Scale(getLumiWeight( (TFile*) fMCFiles->At(j), fXSectionMC.at(j), fLumiMC.at(j)));
        // End

        // Amandeep : nMCfiles was hardcoded here | Replace with channel.compare('combined')
        if( doCorrectChannelFractions && (!fSelfConsistency) && CHAN.CompareTo("combined") == 0 ) {
          CorrectChannelFractions(htemprecoVsgen[j+1], MCFileName, channel_SF);
        }
        htemprecoVsgen[0]->Add(htemprecoVsgen[j+1]);
      }
      fMigMatArray->AddLast(htemprecoVsgen[0]);
    }

    // ***********
    // SYSTEMATICS
    // ***********
    
    // Same for systematics matrix

    if(fsetSyst)
    {
      // Migration Matrix 
      // Amandeep : implementing this a vector breaks 
      // Original :
      TH2D *htemp4[kMaxNumSysFiles];
      // Modified :
      // TH2D *htemp4[nsysfiles];

      // Over the systematics 
      for(int j=0; j < nsysfiles; j++){
        std::cout << "Filling histograms for systematics" << std::endl;

        // Amandeep : nMCfiles was hardcoded here | Replace with channel.compare('combined')
        if( CHAN.CompareTo("combined") == 0 ) 
        {
          // re-evaluate the SFs for each systematic variation
          double syst_yields[4] = {0,0,0,0};

          for(int k=0; k < nMCfiles; k++)
          {
            TString SysFileName = fSysFiles->At(j*nMCfiles+k)->GetName();
            TH1D* hrecosysttemp = (TH1D*) ((TFile*)fSysFiles->At(j*nMCfiles+k))->Get(Form("hreco_%s",fVariableNames[i].Data()));
            // Amandeep : Commenting to experiment with external lumi scaling 
            //hrecosysttemp->Scale(getLumiWeight( (TFile*) fSysFiles->At(j*nMCfiles+k), fXSectionSys.at(j*nMCfiles+k), fLumiSys.at(j*nMCfiles+k)));
            // End

            if( SysFileName.Contains("histosTUnfold_ee_") )   syst_yields[0]+=hrecosysttemp->Integral();
            if( SysFileName.Contains("histosTUnfold_emu_") )  syst_yields[1]+=hrecosysttemp->Integral();
            if( SysFileName.Contains("histosTUnfold_mumu_") ) syst_yields[2]+=hrecosysttemp->Integral();
            syst_yields[3]+=hrecosysttemp->Integral();
          }

          if(fSelfConsistency) CalculateChannelFractions(sig_yields, syst_yields, channel_SF);
          else CalculateChannelFractions(bkgsubdata_yields, syst_yields, channel_SF);

        }

        htemp4[j] = (TH2D*)((TH2D*) ((TFile*)fSysFiles->At(j*nMCfiles))->Get(Form("hrecoVsgen_%s",fVariableNames[i].Data())))->Clone(Form("%sRespMatSys_%s",fVariableNames[i].Data(), fVectorOfValidSystematics[j].c_str()));
        htemp4[j]->Reset();

        // Over the channels/ files ??
        for(int k=0; k < nMCfiles; k++) {
          TString SysFileName = fSysFiles->At(j*nMCfiles+k)->GetName();          
          
          // Amandeep : Adding this here to debug
          std::cout << i  << " , " << j << " , " << k << std::endl;
          std::cout << "Sysfilename :: " << SysFileName << ", " << nMCfiles  << " , " << nsysfiles << " , " << j*nMCfiles+k << std::endl;
          // End         
          
          TH2D* htemp4temp    = (TH2D*) ((TFile*)fSysFiles->At(j*nMCfiles+k))->Get(Form("hrecoVsgen_%s",fVariableNames[i].Data()));
          
          // Amandeep : Commenting to experiment with external lumi scaling 
          //htemp4temp->Scale(getLumiWeight( (TFile*) fSysFiles->At(j*nMCfiles+k), fXSectionSys.at(j*nMCfiles+k), fLumiSys.at(j*nMCfiles+k)));
          // End

          // Amandeep : nMCfiles was hardcoded here | Replace with channel.compare('combined')
          if( doCorrectChannelFractions && CHAN.CompareTo("combined") == 0 ) {
            // Not sure about the math here, need to follow the logic fully
            CorrectChannelFractions(htemp4temp, SysFileName, channel_SF);
          }
          // This adds all channel files of a particular syst to the htemp
          htemp4[j]->Add(htemp4temp);
        }
        fMigMatSysArray->AddLast(htemp4[j]); //consists of more entries than fMigMatArray (nsysfiles*nvariables)
      }  
    }

    if(fClosureTest){
     TH1D* htemptruedata[kMaxNumDataFiles];
     htemptruedata[0]= (TH1D*)((TH1D*)((TFile*)fDataFiles->At(0))->Get(Form("hgen_%s",fVariableNames[i].Data())))->Clone(Form("%sTrueData",fVariableNames[i].Data()));
     for(int j = 1; j < ndatafiles; j++)
     {
       htemptruedata[j]= (TH1D*) ((TFile*)fDataFiles->At(j))->Get(Form("hgen_%s",fVariableNames[i].Data()));
       htemptruedata[0]->Add(htemptruedata[j]);
     }
     fTrueDataArray->AddLast(htemptruedata[0]);
   }
    
  }
  return kTRUE; 
}


// *************************
// CalculateChannelFractions
// *************************

// Amandeep : Need to ensure this makes sense in our new scenario with taus as signal
void MatrixUnfControl::CalculateChannelFractions(double bkgsubdata_yields[], double sig_yields[], double channel_SF[]){

  for (int i_chan = 0; i_chan < 4; ++i_chan)
  {
    if( sig_yields[i_chan] > 0 ) channel_SF[i_chan] = ( bkgsubdata_yields[i_chan]/(bkgsubdata_yields[0]+bkgsubdata_yields[1]+bkgsubdata_yields[2]) ) / ( sig_yields[i_chan]/(sig_yields[0]+sig_yields[1]+sig_yields[2]) );
    else channel_SF[i_chan] = 0;
  }

  std::cout<<"channel SF ee emu mumu combined: "<<channel_SF[0]<<" "<<channel_SF[1]<<" "<<channel_SF[2]<<" "<<channel_SF[3]<<std::endl;

}


// ***********************
// CorrectChannelFractions
// ***********************

// Amandeep : Need to ensure this makes sense in our new scenario with taus as signal
void MatrixUnfControl::CorrectChannelFractions(TH2D*& htemprecoVsgen, TString MCFileName, double channel_SF[]){

  std::cout << "Correcting channel yield fractions using the following SFs (ee, emu, mumu, combined): " << channel_SF[0] << " " << channel_SF[1] << " " << channel_SF[2] << " " << channel_SF[3] << std::endl;

  // Check integral and error
  double integralN, integralNerror, integralD, integralDerror;

  integralD = htemprecoVsgen->IntegralAndError(1, htemprecoVsgen->GetNbinsX(), 0, htemprecoVsgen->GetNbinsY()+1, integralDerror);
  integralN = htemprecoVsgen->IntegralAndError(1, htemprecoVsgen->GetNbinsX(), 1, htemprecoVsgen->GetNbinsY()  , integralNerror);
  
  std::cout << MCFileName << " IntegralDB: " << integralD << " +/- " << integralDerror << " IntegralNB: " << integralN << " +/- " << integralNerror << std::endl;

  double tempSF=1;
  if( MCFileName.Contains("histosTUnfold_ee_") )   tempSF = channel_SF[0];
  if( MCFileName.Contains("histosTUnfold_emu_") )  tempSF = channel_SF[1];
  if( MCFileName.Contains("histosTUnfold_mumu_") ) tempSF = channel_SF[2];

  double delta_uf[htemprecoVsgen->GetNbinsX()] = {0};
  double deltaerr2_uf[htemprecoVsgen->GetNbinsX()] = {0};


  for(Int_t ri =1; ri <= htemprecoVsgen->GetNbinsX(); ri++)  {
    for(Int_t rj =1; rj <= htemprecoVsgen->GetNbinsY(); rj++) { 
      delta_uf[ri-1]     += (1.-tempSF)*htemprecoVsgen->GetBinContent(ri,rj);
      deltaerr2_uf[ri-1] += (1.-tempSF*tempSF)*htemprecoVsgen->GetBinError(ri,rj)*htemprecoVsgen->GetBinError(ri,rj);
      htemprecoVsgen->SetBinContent(ri,rj,tempSF*htemprecoVsgen->GetBinContent(ri,rj));
      htemprecoVsgen->SetBinError(ri,rj,tempSF*htemprecoVsgen->GetBinError(ri,rj));
    }
  }

  for(Int_t ri =1; ri <= htemprecoVsgen->GetNbinsX(); ri++)  { 
      htemprecoVsgen->SetBinContent(ri,0,delta_uf[ri-1]+htemprecoVsgen->GetBinContent(ri,0));
      htemprecoVsgen->SetBinError(ri,0,sqrt(deltaerr2_uf[ri-1]+htemprecoVsgen->GetBinError(ri,0)*htemprecoVsgen->GetBinError(ri,0)) );
  }

  integralD = htemprecoVsgen->IntegralAndError(1,htemprecoVsgen->GetNbinsX(),0,htemprecoVsgen->GetNbinsY()+1,integralDerror);
  integralN = htemprecoVsgen->IntegralAndError(1,htemprecoVsgen->GetNbinsX(),1,htemprecoVsgen->GetNbinsY(),integralNerror);
  std::cout << MCFileName << " IntegralDA: " << integralD << " +/- " << integralDerror << " IntegralNA: " << integralN << " +/- " << integralNerror << std::endl;

}


void MatrixUnfControl::sigfracBgSubtMethod(TH1D*& data_subt, TH1D* sig, TH1D* bg, TH1D* data){
  int sigfrac;
  sigfrac=sig->Integral()/(sig->Integral()+bg->Integral());
  data_subt = (TH1D*) ((TH1D*) data)->Clone(Form("%s_subtr",((TH1D*) data)->GetName()) );
  data_subt->Scale(sigfrac);
}


void MatrixUnfControl::regularBgSubtMethod(TH1D*& data_subt, TH1D* bg, TH1D* data){

  data_subt = (TH1D*) ((TH1D*) data)->Clone(Form("%s_subtr",((TH1D*) data)->GetName()) );
  
  for(int i=1; i<=data->GetNbinsX(); i++){
    if(bg->GetBinContent(i)>data->GetBinContent(i)) {
      data_subt->SetBinContent(i,0);
      std::cout<<"MatrixUnfControl::regularBgSubtMethod: setting bin content to zero "<<i<<std::endl;
    }
    else data_subt->SetBinContent(i,data->GetBinContent(i) - bg->GetBinContent(i));
    data_subt->SetBinError(i, sqrt(data->GetBinError(i)*data->GetBinError(i) + bg->GetBinError(i)*bg->GetBinError(i))); 
  }
}


Double_t MatrixUnfControl::getLumiWeight(TFile* f, Double_t XSection, Double_t Lumi){

  TH1D *hNrOfEvts   = (TH1D*) f->Get("hNrOfEvts");
  Double_t nrOfEvts = hNrOfEvts->GetBinContent(1);

  Double_t lumiweight = Lumi*XSection/(nrOfEvts); 
  return lumiweight; 
}


// *******
// AddFile
// *******

// Takes in file, xsection and a fileptr 
// where we append the files from the filelist

bool MatrixUnfControl::AddFile(TString filename, double xsec_norm, TObjArray *&arr, std::vector<double>&norm_list, std::vector<double>&lumi_list) {
  if(arr==nullptr){
    arr = new TObjArray();
    norm_list.clear();
  }

  double lumi_to_use;
  if      (filename.Contains("2016preVFP"))  lumi_to_use = 19668.0;
  else if (filename.Contains("2016postVFP")) lumi_to_use = 16778.0;
  else if (filename.Contains("2016")) lumi_to_use = 36446.0;
  else if (filename.Contains("2017"))     lumi_to_use = 41480.0;
  else if (filename.Contains("2018"))     lumi_to_use = 59830.0;
  else if (filename.Contains("fullRun2")) lumi_to_use = 137650.0;
  else {
    std::cerr << "ERROR in MatrixUnfControl::AddFile() : no year string in" << filename  << std::endl;
    exit(1);
  }

  std::cout<<"Opening file: "<<filename.Data()<<std::endl;
  TFile *file = new TFile(filename.Data(), "READ");
  if(file==0){
    std::cout<<"couldnt open "<<filename<<std::endl;
    return 0;
  }
  else{
    //pointers to the files are added in an TObjArray for background and data files that have more than one sample files
    arr->AddLast(file);
    norm_list.push_back(xsec_norm);
    lumi_list.push_back(lumi_to_use);
  }
  return kTRUE;
}

bool MatrixUnfControl::SetTauOptions(bool constrainTau,float tausys, double taumin, double taumax, bool minrhoavg){
  fConstrainTau = constrainTau;
  fTauSys       = tausys;
  // Change the values from 1e-20 (small regularization) 
  // to 1(completely regularized), 0 searches for itself
  fTauMin    = taumin;
  fTauMax    = taumax;
  fMinRhoAVG = minrhoavg;
  return kTRUE;
}

bool MatrixUnfControl::SetOptions(std::string era, double lumi, bool closuretest, bool selfconsistency, bool ensembletest, bool linearitytest,  bool allowrebin, bool combinedsumptt, bool subtrbgmcfromdata, bool dosubtrbg, bool usebiasintunfold, int useareaconstrain, int regmode, bool setSyst){
  // ************************************************************************
  // these are the flags that can be tweaked during the running of the script
  // ************************************************************************
  
  //used for normalization

  fEra  = era;
  fLumi = lumi;
  fClosureTest     = closuretest;
  fSelfConsistency = selfconsistency;
  fEnsembleTest    = ensembletest;
  fLinearityTest   = linearitytest;
  fAllowRebinOption  = allowrebin;
  fCombinedSumptt    = combinedsumptt;
  fUseBiasInUnfold   = usebiasintunfold;
  fUseAreaConstrain  = useareaconstrain;
  fSubtrBgMCFromData = subtrbgmcfromdata; // Manual background subtraction
  fDoSubtrBg = dosubtrbg;                 // TUnfold background subtraction
  // change the fRegMode - see MatrixUnf.C
  fRegMode = regmode; // 1=size,2=derivative,3=curvature,4=mixed,5=densityCurvature
  fsetSyst = setSyst;

  return kTRUE;
}

// For single variable
bool MatrixUnfControl::AddUnfoldVariable(std::string varname){
  // In both cases the varname is of type string
  fVariableNames.push_back(varname);
  return kTRUE;
}

// For entire list
bool MatrixUnfControl::AddUnfoldVariables(TString varlist){
  // Read fVariableNames from this file
  std::ifstream readvar;
  readvar.open(varlist.Data());
  if(!readvar.is_open()){
    std::cout<<"The file with the list of variables is not open!"<<std::endl;
    exit(0);
  }

  while(!readvar.eof()){
    std::string varname;
    getline(readvar, varname);
    if(varname=="") continue;
    // In both cases the varname is of type string
    fVariableNames.push_back(varname);
  }
  return kTRUE;
}

bool MatrixUnfControl::AddSystematics(std::vector<std::string> vecinputsyst, std::vector<std::string> vecsyst){
  fVectorOfInputSystematics = vecinputsyst;
  fVectorOfValidSystematics = vecsyst;
}

bool MatrixUnfControl::MakePseudoHistos(TH1D *inputhisto, TH1D*& pseudohisto, TString varname){

  pseudohisto= (TH1D*)inputhisto->Clone(Form("hpseudo_%s_norm",varname.Data()));
  
  for(int x_bin=1; x_bin<=inputhisto->GetNbinsX(); x_bin++){
    
    float mu = inputhisto->GetBinContent(x_bin);
    float u = gMyRandom.PoissonD(mu);
    
    if(u>0) {
      pseudohisto->SetBinContent(x_bin, u);
      pseudohisto->SetBinError(x_bin,TMath::Sqrt(u));
    }
    if(u<0) {
         pseudohisto->SetBinContent(x_bin, 0);
         pseudohisto->SetBinError(x_bin, 0);
    }
  }
  return kTRUE;
}

/**
 * @brief Create the linearity histogram inputs to be used to check for bias in the unfolding method
 * 
 * @param inputmighisto The input nominal migration matrix that will be reweighted with a linear slope
 * @param linearityrecohisto The linearly reweighted reco-level distribution that will be fed as input in the unfolding procedure to test
 *                           for bias in the unfolding procedure. This is simply the projection of the reweighted migration matrix onto the y-axis.
 * @param linearitygenhisto  The linearly reweighted gen-level distribution that will be used to compare the unfolded linearly reweighted 
 *                           distribution against. This is simply the projection of the reweighted migration matrix onto the x-axis.
 * @param varname The name of the variable that is being linearly reweighted. Used to determine how much "slope" to inject for a change in the variable.
 * @param deltaCoefficient The change in the coefficient that is being added to the distribution. 
 * @param rebinnedBinning The rebinned binning of the unfolded/gen distribution
 * @param originalBinning The fine binning of the unfolded/gen distribution
 * @param numCoefficients The number of coefficients within the distribution
 * @param rebin The amount to rebin the distribution by when performing the linear reweighting. This is only used if we use the coarse binning to reweight by
 * @param reweightusingwidebincentres Whether to perform the linear reweighting with the coarse binning
 * @param binToReweight What bin to reweight. Default is -1 which means reweight the entire distribution
 * @return true Always returns true
 * @return false Never returns false
 */
bool MatrixUnfControl::MakeLinearityHistos(TH2D *inputmighisto, TH1D *&linearityrecohisto, TH1D *&linearitygenhisto, TString varname, float deltaCoefficient, const TUnfoldBinning* rebinnedBinning, const TUnfoldBinning* originalBinning, int numCoefficients, int rebin, bool reweightusingwidebincentres, int binToReweight) {
  std::cout << "Making the linearity histograms" << std::endl;
  std::cout << "Initial deltaCoefficient = " << deltaCoefficient << std::endl;
  // because the reweighting in this function changes the value of the (Cij+-Cji) coefficient by 2*deltaCoefficient
  if (MatrixUnf::isCPMtype(varname)) deltaCoefficient /= 2.;                                                                 
  
  if (MatrixUnf::isCtype(varname) || MatrixUnf::isCPMtype(varname) || MatrixUnf::isDtype(varname) || MatrixUnf::isLinearCombtype(varname)) deltaCoefficient *= -1.;  // sign convention for spin correlation coefficients

  TH2D *weightedmighisto = (TH2D *)inputmighisto->Clone(Form("hlinearity_ReweightedMig_%s_%f", varname.Data(), deltaCoefficient));
  // weightedmighisto->Sumw2();
  TH1D *inputgenhisto = (TH1D *)inputmighisto->ProjectionX("_px");
  // TH1D* inputrecohisto = (TH1D*) inputmighisto->ProjectionY("_py",1,NbinsY);
  TH1D *rebinnedgenhisto = (TH1D *)inputgenhisto->Clone(Form("hlinearity_Gen_rebinned_%s_%f", varname.Data(), deltaCoefficient));
  if (rebin > 1) {
    MatrixUnf::rebinMultidimensionalInput(rebinnedgenhisto, rebinnedBinning->FindNode("ttbargen_rebinnedB"), originalBinning->FindNode("ttbargen"), rebin);    
  }

  // Set up number of bin information
  int dim    = originalBinning->GetDistributionDimension();
  int NbinsX = inputmighisto->GetNbinsX();
  int NbinsY = inputmighisto->GetNbinsY();
  int numGenBinsPerCoefficient  = NbinsX / numCoefficients;
  int numRecoBinsPerCoefficient = NbinsY / numCoefficients;

  int numLinearityCoefficients = 0;
  if (reweightusingwidebincentres) {
    numLinearityCoefficients = numCoefficients / pow(rebin, dim - 1);
  } 
  else {
    numLinearityCoefficients = numCoefficients;
  }

  // grab the bin edges
  const TVectorD *vec_rebinnedBinEdges = rebinnedBinning->FindNode("ttbargen_rebinnedB")->GetDistributionBinning(0);
  const double *rebinnedBinEdges       = vec_rebinnedBinEdges->GetMatrixArray();
  const TVectorD *vec_originalBinEdges = originalBinning->FindNode("ttbargen")->GetDistributionBinning(0);
  const double *binEdges = vec_originalBinEdges->GetMatrixArray();
  
  double centre = (binEdges[numGenBinsPerCoefficient] + binEdges[0]) / 2;
  double width  =  binEdges[numGenBinsPerCoefficient] - binEdges[0];
  
  std::cout << "Center of distribution = " << centre << ", and width = " << width << std::endl;
  
  for (int i = 0; i < numLinearityCoefficients; i++) {
    std::cout << "Linear reweighting of " << i << "th-variable portion of the distribution" << std::endl;
    int binOffset = i * numGenBinsPerCoefficient;
    int rebinnedBinOffset   = i * numGenBinsPerCoefficient / rebin;
    double inputgenintegral = inputgenhisto->Integral(binOffset + 1, binOffset + numGenBinsPerCoefficient);

    std::cout << "Total number of events in this portion of the distribution: " << inputgenintegral << std::endl;
    // double inputrecointegral = inputrecohisto->Integral(i * numRecoBinsPerCoefficient + 1, (i + 1) * numRecoBinsPerCoefficient);

    for (int x_bin = 0; x_bin < numGenBinsPerCoefficient; x_bin++) {
      int binIndex      = binOffset + x_bin + 1;
      double binCenter  = (binEdges[x_bin] + binEdges[x_bin + 1]) / 2;
      double binContent = inputgenhisto->GetBinContent(binIndex);
      double symmetricBinContent = inputgenhisto->GetBinContent(binOffset + numGenBinsPerCoefficient - x_bin);
      double binWidth   = binEdges[x_bin + 1] - binEdges[x_bin];

      // use the same slope for all fine bins in the same large bin
      if (reweightusingwidebincentres) {
        binIndex   = rebinnedBinOffset + (x_bin / rebin) + 1;
        binCenter  = rebinnedgenhisto->GetXaxis()->GetBinCenter(binIndex);
        binContent = rebinnedgenhisto->GetBinContent(binIndex);
        symmetricBinContent = rebinnedgenhisto->GetBinContent(rebinnedBinOffset + (numGenBinsPerCoefficient - x_bin + 1) / rebin + 1);
        binWidth   = rebinnedBinEdges[binIndex] - rebinnedBinEdges[binIndex - 1];
      }

      // linear reweighting based on bin centre (to avoid bias from dependence of the A and S matrices on the slope, so that we are only testing for bias from the regularisation)
      double binslope = deltaCoefficient * (binCenter - centre) / (width / 2.);

      // double addedslope = binslope*inputgenintegral*binWidth/width;   //add a constant slope function (deltaCoefficient*x) to the distribution (so that regularisation of curvature should be unbiased in default TUnfold)
      // if(reweightusingwidebincentres) addedslope = binslope*inputgenintegral*rebinnedgenhisto->GetXaxis()->GetBinWidth((x_bin-1)/rebin + 1)/width;   //add a constant slope function (deltaCoefficient*x) to the distribution (so that regularisation of curvature should be unbiased in default TUnfold)

      // Inject shape using the expected change in functional form from a change in coefficient

      double addedslope = 0.;
      if (MatrixUnf::isBPMtype(varname)) { 
        //  For the b_P/M* variables a change in sum/difference of coefficients always leads to the addition of 
        // (constant*x)*(2 - abs(x)) to the functional form:
        addedslope = deltaCoefficient * inputgenintegral / 8. * (2. - fabs(binCenter)) * binCenter * binWidth;
      } 
      else if (MatrixUnf::isCtype(varname) || MatrixUnf::isCPMtype(varname) || MatrixUnf::isBtype(varname) || MatrixUnf::isDtype(varname) || MatrixUnf::isLinearCombtype(varname)) {
        //  Add (deltaCoefficient*x)*f_sym(x) 
        // (Realistic NP case for the change in shape for the other spin density matrix variables)
        addedslope = binslope * (binContent + symmetricBinContent) / 2.;  
      } else if (MatrixUnf::isOtherSymType(varname)) {
        // For other variables where we assume the NP shows up roughly as the addition of a linear shape 
        // (ith_variable.e. dist_injected = dist_orig + constant*x)
        addedslope = binslope * inputgenintegral * binWidth / width;
      } else {
        // Don't symmetrise the SM shape for the remaining (non-spin) variables, ith_variable.e. add (deltaCoefficient*x)*f(x)  (For other variables in general the expected functional form of the change in shape is unknown. 
        // This injection of shape has no physical motivation, and is chosen to result in "unbiased" regularisation (as measured in the linearity tests) when the unsymmetrised BinFactorFunction is used)
        addedslope = binslope * (binContent);  
      }

      std::cout << x_bin << ": added slope = " << addedslope << ", binslope = " << binslope << ", binContent = " << binContent << ", binCenter = " << binCenter << std::endl;

      double weight = 1. + addedslope / binContent;
      //std::cout<<"NbinsX,GetNbinsY,centre,width,inputgenintegral,binslope,addedslope,weight: "<<NbinsX<<" "<<NbinsY<<" "<<centre<<" "<<width<<" "<<inputgenintegral<<" "<<binslope<<" "<<addedslope<<" "<<weight<<std::endl;

      for (int y_bin = 0; y_bin <= NbinsY; y_bin++) {
        double u = inputmighisto->GetBinContent(binOffset + x_bin + 1, y_bin);
        double delU = inputmighisto->GetBinError(binOffset + x_bin + 1, y_bin);
        if (binToReweight < 1  // default binToReweight = -1
          || (x_bin - 1) / rebin == (binToReweight - 1)
          || (binToReweight == 998 && x_bin <= NbinsX / 2)
          || (binToReweight == 999 && x_bin > NbinsX / 2)) { 
          u *= weight;
          delU *= weight;
        }

        if (u > 0) {
          weightedmighisto->SetBinContent(binOffset + x_bin + 1, y_bin, u);
          weightedmighisto->SetBinError(binOffset + x_bin + 1, y_bin, delU);
        }
        if (u < 0) {
          weightedmighisto->SetBinContent(x_bin, y_bin, 0);
          weightedmighisto->SetBinError(x_bin, y_bin, 0);
        }
      }
    }
  }
  
  linearityrecohisto = (TH1D *)weightedmighisto->ProjectionY("_py", 1, NbinsY, "e");
  linearitygenhisto = (TH1D *)weightedmighisto->ProjectionX("_px", 0, -1, "e");

  // double outputgenintegral = linearitygenhisto->Integral();
  // double outputrecointegral = linearityrecohisto->Integral();
  // std::cout<<"LinIntegrals: "<<outputgenintegral/inputgenintegral<<" "<<outputrecointegral/inputrecointegral<<std::endl;

  // linearityrecohisto->Scale(inputgenintegral/outputgenintegral);
  // linearitygenhisto->Scale(inputgenintegral/outputgenintegral);

  // linearityrecohisto->Scale(inputrecointegral/outputrecointegral);
  // linearitygenhisto->Scale(inputrecointegral/outputrecointegral);

  std::cout << "Made one linearity histo" << std::endl;
  std::cout << "linearityrecohisto->GetName(): " << linearityrecohisto->GetName() << std::endl;

  return kTRUE;
}


void MatrixUnfControl::MakeAllPlots(TString TUNFFILE, TString  PLOTFILE, TString CHAN, bool doOptimiseTotalUnc){

  // Note that this function also calculates the systematics and fits the anomalous couplings

  common::setHHStyle(*gStyle);

  bool divBinWidth = kTRUE; //apply bin width correction to unfolded plots, covariance matrices, and yaml tables (note I only made certain to apply this to the plots for the paper/hepdata)

  TString SystFilename = TUNFFILE+"Systematics.root";
  TFile *systoutfile;
  if(fsetSyst) systoutfile = new TFile(SystFilename.Data(),"RECREATE");

  std::vector<TString>backgroundnames;
  backgroundnames.push_back("DYeemm");
  backgroundnames.push_back("DYtautau");
  //backgroundnames.push_back("ttbarbgviatau");
  backgroundnames.push_back("singletop");
  backgroundnames.push_back("ww");
  backgroundnames.push_back("wz");
  backgroundnames.push_back("zz");
  backgroundnames.push_back("ttbarW");
  backgroundnames.push_back("ttbarZ");
  backgroundnames.push_back("ttbarbg");
  backgroundnames.push_back("wtolnu");
  backgroundnames.push_back("other");

  // Amandeep : Here it seems to add backgrounds to systematics if we want background subtraction
  std::vector<std::string> tempVectorOfInputSystematics = fVectorOfInputSystematics;
  if(fDoSubtrBg) {
    for (unsigned int ibkg = 0; ibkg < backgroundnames.size(); ++ibkg)
    {
      tempVectorOfInputSystematics.push_back(backgroundnames.at(ibkg).Data());
    }
  }

  std::vector<TString>BKGnames;
  BKGnames.push_back("BKGDYeemm");
  BKGnames.push_back("BKGDYtautau");
  //BKGnames.push_back("BKGttbar");
  BKGnames.push_back("BKGsingletop");
  BKGnames.push_back("BKGww");
  BKGnames.push_back("BKGwz");
  BKGnames.push_back("BKGzz");
  BKGnames.push_back("BKGttbarW");
  BKGnames.push_back("BKGttbarZ");
  BKGnames.push_back("BKGttbarother");
  BKGnames.push_back("BKGwtolnu");
  BKGnames.push_back("BKGother");

  std::vector<std::string> tempVectorOfInputSystematics_BKGnames = fVectorOfInputSystematics;
  if(fDoSubtrBg) {
    for (unsigned int ibkg = 0; ibkg < BKGnames.size(); ++ibkg)
    {
      tempVectorOfInputSystematics_BKGnames.push_back(BKGnames.at(ibkg).Data());
    }
  }

  std::vector<TString>flatsystnames;
  flatsystnames.push_back("Lumi");
  flatsystnames.push_back("BR");

  bool doFlatSysts = kTRUE;
  if(doFlatSysts) {
    for (unsigned int iflat = 0; iflat < flatsystnames.size(); ++iflat)
    {
      tempVectorOfInputSystematics.push_back(flatsystnames.at(iflat).Data());
      tempVectorOfInputSystematics_BKGnames.push_back(flatsystnames.at(iflat).Data());
    }
  }

  const int nsyst = tempVectorOfInputSystematics.size();
  // TODO: Hard-coded number of bins
  // can't really change this without reworking all of the multidimensional arrays into vectors
  const int nbinsrhoi_temp = 24;
  const int numCoefficients = nbinsrhoi_temp / 6;
  const int nbinstotal = fVariableNames.size() * nbinsrhoi_temp;
  const int nbinstotal_norm = fVariableNames.size() * (nbinsrhoi_temp - 1);
  
  // Amandeep : Why are chi2 numbers are hardcoded ?
  double systTable[4][fVariableNames.size()][numCoefficients][nsyst + 6] = {0.};
  double coefTable[4][fVariableNames.size()][numCoefficients]            = {0.};
  std::vector<double> fSMTable;
  std::vector<double> fSMTableErr;
  std::vector<double> fSMTableStatErr;
  std::vector<double> fSMTableSystErr;
  std::vector<TString> fSMTableNames;
  double MCcoefTable[fVariableNames.size()][numCoefficients][3] = {0.};
  double AMCATNLOcoefTable[fVariableNames.size()][numCoefficients][3] = {0.};
  double binsumsystTable[4][fVariableNames.size()][numCoefficients][nsyst + 6] = {0.};
  double binsumTable[4][fVariableNames.size()][numCoefficients] = {0.};
  double Chi2values[2][2][fVariableNames.size()][8]    = {0.};
  double CoefChi2values[3][fVariableNames.size()][8]   = {0.};
  double AllVarsChi2values[3][fVariableNames.size()][8] = {0.};

  std::vector<TString>combsystnames_temp; //categories for combining systematic subcomponents
  combsystnames_temp.push_back("JES");
  combsystnames_temp.push_back("PDF");
  combsystnames_temp.push_back("BKG");

  std::vector<TString>combsystnames;

  bool addedcombtolist[combsystnames_temp.size()] = {0};
  for (unsigned int syst = 0; syst < tempVectorOfInputSystematics_BKGnames.size(); ++syst) {

    if(tempVectorOfInputSystematics_BKGnames[syst] == "AMCATNLOFXFX" || tempVectorOfInputSystematics_BKGnames[syst] == "POWHEGV2HERWIG") continue; //systematics not to include in deltaSyst plots or table of combined systs

    bool addtolist = kTRUE;
    int i_comb = -1;
    for (unsigned int i_cs = 0; i_cs < combsystnames_temp.size(); ++i_cs) {
      if( (TString(tempVectorOfInputSystematics_BKGnames[syst]).BeginsWith( combsystnames_temp.at(i_cs).Data() ) )  ) {
        addtolist = kFALSE;
        i_comb = i_cs;
        break;
      }
    }
    if(addtolist) {
      combsystnames.push_back(tempVectorOfInputSystematics_BKGnames[syst].c_str());
      std::cout<<tempVectorOfInputSystematics_BKGnames[syst].c_str()<<std::endl;
    }
    else if ( addedcombtolist[i_comb] == kFALSE ) {
      combsystnames.push_back( combsystnames_temp.at(i_comb) );
      std::cout<<combsystnames_temp.at(i_comb)<<std::endl;
      addedcombtolist[i_comb] = kTRUE;
    }
  }

  double combsystTable[4][fVariableNames.size()][numCoefficients][combsystnames.size()] = {0.};
  double binsumcombsystTable[4][fVariableNames.size()][numCoefficients][combsystnames.size()] = {0.};


  TH1D* hdeltasquaredStat_AllVar[4];

  // Amandeep : coefficient values for all spincorr variables ?
  // Shouldn't be hardcoded 0-21 since order of spin corr variables can change
  double NLOWvalues[38]       = {4.0e-3, 4.0e-3, 1.6e-3, 1.6e-3, 5.7e-3, 5.7e-3, 1e-3  , 1e-3  , 1e-3  , 1e-3  , 0.331, 0.071, 0.326, -0.206, 0, 1.06e-3, 0, 2.15e-3, 0, -0.242666667, 0.180510 , 0.107535, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  double NLOWvaluesUnCorr[38] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0.128998 , 0.200713, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  // array to use in anomalous couplings fits. Values from https://arxiv.org/abs/1508.05271
  const int nFit = 7;
  double mutTableMulti[nFit][38] = {
    {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.917,2.475,2.025,0.740,0.,0.,0.,0.,0.,-1.805666667,0.,0., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // mu_t
    {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-1.218,-0.697,-0.0799,-0.306,0.,0.,0.,0.,0.,(1.218+0.697+0.0799)/3.,0.,0., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // cVV
    {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-4.143,0.,-0.800,0.,0.,0., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // d_t
    {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-1.226,0.,-2.157,0.,0.,0., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // c--
    {1.607/2.,1.607/2.,0.210/2.,0.210/2.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // cVA
    {0.,0.,0.,0.,0.,0.,0.883/2.,0.883/2.,0.792/2.,0.792/2.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // cAV
    {0.,0.,0.,0.,4.876/2.,-4.876/2.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} // c-+
  };

  TMatrixD Measvec_allVars(nbinstotal,1);
  TMatrixD Theoryvec_allVars(nbinstotal,1);
  // TMatrixD AMCATNLOvec_allVars(nbinstotal,1);
  // TMatrixD POWHEGV2HERWIGvec_allVars(nbinstotal,1);
  // TMatrixD TOP_PTvec_allVars(nbinstotal,1);
  TMatrixD NLOWvec_allVars(nbinstotal,1);
  TMatrixD NLOWuncorrvec_allVars(nbinstotal,1);
  TMatrixD NLOW3vec_allVars(nbinstotal,1);
  TMatrixD NLOW4vec_allVars(nbinstotal,1);

  std::vector<TMatrixD> NLOWvecMulti_allVars;
  for (int iM = 0; iM < nFit; ++iM)
  {
    TMatrixD NLOWvecMulti_allVars_temp(nbinstotal,1);
    NLOWvecMulti_allVars.push_back(NLOWvecMulti_allVars_temp);
  }


  TMatrixD Measvec_norm_allVars(nbinstotal_norm,1);
  TMatrixD Theoryvec_norm_allVars(nbinstotal_norm,1);
  // TMatrixD AMCATNLOvec_norm_allVars(nbinstotal_norm,1);
  // TMatrixD POWHEGV2HERWIGvec_norm_allVars(nbinstotal_norm,1);
  // TMatrixD TOP_PTvec_norm_allVars(nbinstotal_norm,1);
  TMatrixD NLOWvec_norm_allVars(nbinstotal_norm,1);
  TMatrixD NLOWuncorrvec_norm_allVars(nbinstotal_norm,1);
  TMatrixD NLOW3vec_norm_allVars(nbinstotal_norm,1);
  TMatrixD NLOW4vec_norm_allVars(nbinstotal_norm,1);

  std::vector<TMatrixD> NLOWvecMulti_norm_allVars;
  for (int iM = 0; iM < nFit; ++iM)
  {
    TMatrixD NLOWvecMulti_norm_allVars_temp(nbinstotal_norm,1);
    NLOWvecMulti_norm_allVars.push_back(NLOWvecMulti_norm_allVars_temp);
  }


  TMatrixD Meas_Coefs(fVariableNames.size(), numCoefficients);
  TMatrixD Theory_Coefs(fVariableNames.size(), numCoefficients);
  TMatrixD AMCATNLO_Coefs(fVariableNames.size(), numCoefficients);
  TMatrixD POWHEGV2HERWIG_Coefs(fVariableNames.size(), numCoefficients);
  TMatrixD TOP_PT_Coefs(fVariableNames.size(), numCoefficients);
  TMatrixD NLOW_Coefs(fVariableNames.size(), numCoefficients);
  TMatrixD NLOWuncorr_Coefs(fVariableNames.size(), numCoefficients);
  TMatrixD NLOW3_Coefs(fVariableNames.size(), numCoefficients);
  TMatrixD NLOW4_Coefs(fVariableNames.size(), numCoefficients);

  TMatrixD Meas_Coefs_abs(fVariableNames.size(), numCoefficients);

  TMatrixD Meas_BinSums(fVariableNames.size(), numCoefficients);



  TH2D* hTUnfoldDataStatCorrM_AllVar[4];
  TH2D* hTUnfoldTotalStatCorrM_AllVar[4];

  for (int i_n = 0; i_n < 4; ++i_n)
  {
    hTUnfoldDataStatCorrM_AllVar[i_n] = new TH2D(Form("hTUnfoldDataStatCorrM_AllVar_%d",i_n), Form("hTUnfoldDataStatCorrM_AllVar_%d",i_n), nbinstotal,0,nbinstotal, nbinstotal,0,nbinstotal);
    hTUnfoldDataStatCorrM_AllVar[i_n]->SetDirectory(0);

    hTUnfoldTotalStatCorrM_AllVar[i_n] = new TH2D(Form("hTUnfoldTotalStatCorrM_AllVar_%d",i_n), Form("hTUnfoldTotalStatCorrM_AllVar_%d",i_n), nbinstotal,0,nbinstotal, nbinstotal,0,nbinstotal);
    hTUnfoldTotalStatCorrM_AllVar[i_n]->SetDirectory(0);
  }

  TH1D* hist_acombfactor_allVar[4][nbinsrhoi_temp/2];

  for (int i_a = 0; i_a < nbinsrhoi_temp/2; ++i_a)
  {
    for (int i_n = 0; i_n < 4; ++i_n)
    {
      hist_acombfactor_allVar[i_n][i_a] = new TH1D(Form("hist_acombfactor_allVar%d_%d",i_a,i_n), Form("hist_acombfactor_allVar%d_%d",i_a,i_n), fVariableNames.size(),0,fVariableNames.size());
      hist_acombfactor_allVar[i_n][i_a]->SetDirectory(0);
    }
  }


  TH1D* hpowheg_AllVar[2];
  TH1D* hmcatnlo_AllVar[2];

  for (int inorm = 0; inorm <= 1; ++inorm)
  {
    hpowheg_AllVar[inorm]  = new TH1D(Form("hpowhegAllVar%s",(inorm?"Norm":"")), Form("hpowhegAllVar%s",(inorm?"Norm":"")), fVariableNames.size()*nbinsrhoi_temp,0,fVariableNames.size()*nbinsrhoi_temp);
    hpowheg_AllVar[inorm]->SetDirectory(0);
    hmcatnlo_AllVar[inorm] = new TH1D(Form("hmcatnloAllVar%s",(inorm?"Norm":"")), Form("hmcatnloAllVar%s",(inorm?"Norm":"")), fVariableNames.size()*nbinsrhoi_temp,0,fVariableNames.size()*nbinsrhoi_temp);
    hmcatnlo_AllVar[inorm]->SetDirectory(0);
  }

  TH1D* hbinwidth_AllVar = new TH1D(Form("hbinwidth_AllVar"), Form("hbinwidth_AllVar"), fVariableNames.size()*nbinsrhoi_temp,0,fVariableNames.size()*nbinsrhoi_temp);
  hbinwidth_AllVar->SetDirectory(0);

  TH1D* hbinlowedge_AllVar = new TH1D(Form("hbinlowedge_AllVar"), Form("hbinlowedge_AllVar"), fVariableNames.size()*nbinsrhoi_temp,0,fVariableNames.size()*nbinsrhoi_temp);
  hbinlowedge_AllVar->SetDirectory(0);

  std::vector<const TUnfoldBinning*> generatorBinningsRebinned;
  std::vector<const TUnfoldBinning*> generatorBinnings;
  std::vector<const TUnfoldBinning*> detectorBinningsRebinned;
  std::vector<const TUnfoldBinning*> detectorBinnings;
  for (unsigned int ithVariable = 0; ithVariable < fVariableNames.size(); ithVariable++) {
    TDOMParser parser;
    Int_t error = parser.ParseFile(Form("binning/" + fVariableNames[ithVariable] + "_binning_rebinnedB.xml"));
    if(error) std::cout << "error=" << error << " from TDOMParser\n";

    TXMLDocument const* XMLdocument     = parser.GetXMLDocument();
    TUnfoldBinningXML* generatorBinning = TUnfoldBinningXML::ImportXML(XMLdocument,"generator_rebinnedB");
    TUnfoldBinningXML* detectorBinning = TUnfoldBinningXML::ImportXML(XMLdocument,"detector_rebinnedB");

    generatorBinningsRebinned.push_back(generatorBinning->FindNode("ttbargen_rebinnedB"));
    detectorBinningsRebinned.push_back(detectorBinning->FindNode("ttbarreco_rebinnedB"));

    error = parser.ParseFile(Form("binning/" + fVariableNames[ithVariable] + "_binning.xml"));
    if(error) std::cout << "error=" << error << " from TDOMParser\n";

    XMLdocument     = parser.GetXMLDocument();
    TUnfoldBinningXML* generatorBinning_fine = TUnfoldBinningXML::ImportXML(XMLdocument,"generator");
    TUnfoldBinningXML* detectorBinning_fine = TUnfoldBinningXML::ImportXML(XMLdocument,"detector");

    generatorBinnings.push_back(generatorBinning_fine->FindNode("ttbargen"));
    detectorBinnings.push_back(detectorBinning_fine->FindNode("ttbarreco"));
  }


  for(unsigned int i = 0 ; i<fVariableNames.size(); i++){
    if(fVariableNames[i] == "") continue;  // histograms for TUnfold

    // Amandeep : For extracting the number of dimensions
    // From AJ
    int n_sub_bins = 2;
    if (fAllowRebinOption) {
        std::ifstream in(Form("binning/%s_rebin.txt", fVariableNames[i]));
        in >> n_sub_bins;
        in.close();
        std::cout << "Rebin: " << n_sub_bins << std::endl;
    } else
        n_sub_bins = 2;

    const int rebinfine = n_sub_bins;

    TDOMParser parser;
    Int_t error = parser.ParseFile("binning/" + fVariableNames[i] + "_binning_rebinnedB.xml");
    if(error) std::cout << "error=" << error << " from TDOMParser\n";

    TXMLDocument const* XMLdocument     = parser.GetXMLDocument();
    TUnfoldBinningXML* generatorBinning = TUnfoldBinningXML::ImportXML(XMLdocument,"generator_rebinnedB");

    if(!generatorBinning)  std::cout<<"could not read 'generator_rebinnedB' binning\n";

    const TUnfoldBinning *binning = generatorBinning->FindNode("ttbargen_rebinnedB");
    const int numDimensions       = binning->GetDistributionDimension();
    std::cout << numDimensions << " " << rebinfine << std::endl;

    error = parser.ParseFile(Form("binning/" + fVariableNames[i] + "_binning.xml"));
    if(error) std::cout << "error=" << error << " from TDOMParser\n";

    XMLdocument     = parser.GetXMLDocument();
    TUnfoldBinningXML* generatorBinning_fine = TUnfoldBinningXML::ImportXML(XMLdocument,"generator");
    TUnfoldBinningXML* detectorBinning_fine = TUnfoldBinningXML::ImportXML(XMLdocument,"detector");

    if(!generatorBinning_fine) std::cout<<"could not read 'generator' binning\n";
    if(!detectorBinning_fine)  std::cout<<"could not read 'detector' binning\n";

    const TVectorD *binning1DTVectorD = binning->GetDistributionBinning(0);
    binning1DTVectorD->Print();

    if (numDimensions > 1) {
      const TVectorD *binning2DTVectorD = binning->GetDistributionBinning(1);
      binning2DTVectorD->Print();
    }

    // End 

    gStyle->SetPadLeftMargin(0.11);

    //TString TUnfoldResultsFilename = Form(TUNFFILE + fVariableNames[i] + "_constrainedtau%d_regMode%d_closure%d_selfconsistency%d_minrhoavg%d.root", fConstrainTau, fRegMode, fClosureTest, fSelfConsistency, fMinRhoAVG);
    TString TUnfoldResultsFilename   = Form(TUNFFILE + fVariableNames[i] + ".root");
    
    printf("Open file root-File %s\n\n", TUnfoldResultsFilename.Data());
    TFile *unffile = new TFile(TUnfoldResultsFilename.Data(),"READ");                                                                                 
    if(unffile==NULL) exit(0);                                                                                                                                          

    
    TString quantity     = MatrixUnf::ObsName(fVariableNames[i],kFALSE);
    TString quantitycoef = MatrixUnf::CoefName(fVariableNames[i],kFALSE);

    TCanvas *c1 = new TCanvas("c1", "c1");
    TH1D* tunfresult = (TH1D*)unffile->Get(fVariableNames[i]+"TUnfResult");
    tunfresult->SetLineColor(kRed);      
    tunfresult->SetLineWidth(2);      
    tunfresult->SetMarkerColor(kRed);      
    tunfresult->SetMarkerStyle(20);      
    tunfresult->SetMarkerSize(0.5);
    tunfresult->SetTitle("");      
    tunfresult->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    tunfresult->GetXaxis()->SetNdivisions(-20604);
    if(fVariableNames[i].Contains("_mttbar")) {
      c1->SetGridx();
      tunfresult->GetXaxis()->SetNdivisions(-20608);
    }
    SetAxis(fVariableNames[i], tunfresult->GetXaxis());
    //tunfresult->SetMinimum(0);
    tunfresult->SetMaximum(tunfresult->GetMaximum()*1.4);
    tunfresult->GetYaxis()->SetTitle("Unfolded (uncorrected) number of events");
    tunfresult->GetXaxis()->SetTitleOffset(1.36);
    tunfresult->GetYaxis()->SetTitleOffset(1.6);
    tunfresult->Draw();

    Draw2DLabels(fVariableNames[i]);
    DrawCMSLabels(0.04, fLumi);
    if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

    c1->SaveAs(Form("%s/UnfoldedResultsUncorr_%s.pdf",PLOTFILE.Data(),fVariableNames[i].Data()));
    c1->SaveAs(Form("%s/UnfoldedResultsUncorr_%s.root",PLOTFILE.Data(),fVariableNames[i].Data()));
    c1->SaveAs(Form("%s/UnfoldedResultsUncorr_%s.C",PLOTFILE.Data(),fVariableNames[i].Data()));

    TCanvas *c1A = new TCanvas("c1A", "c1A");
    TH1D* tunfresult_rebinnedA = (TH1D*)unffile->Get(fVariableNames[i]+"TUnfResult_rebinnedA");
    tunfresult_rebinnedA->SetLineColor(kRed);      
    tunfresult_rebinnedA->SetLineWidth(2);      
    tunfresult_rebinnedA->SetMarkerColor(kRed);      
    tunfresult_rebinnedA->SetMarkerStyle(20);      
    tunfresult_rebinnedA->SetMarkerSize(0.5);
    tunfresult_rebinnedA->SetTitle("");      
    tunfresult_rebinnedA->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    tunfresult_rebinnedA->GetXaxis()->SetNdivisions(-20604);

    if(fVariableNames[i].Contains("_mttbar")) {
      tunfresult_rebinnedA->GetXaxis()->SetNdivisions(-20604);
      c1A->SetGridx();
    }
    
    SetAxis(fVariableNames[i], tunfresult_rebinnedA->GetXaxis());
    //tunfresult_rebinnedA->SetMinimum(0);
    tunfresult_rebinnedA->SetMaximum(tunfresult_rebinnedA->GetMaximum()*1.4);
    tunfresult_rebinnedA->GetYaxis()->SetTitle("Unfolded (uncorrected) number of events");
    tunfresult_rebinnedA->GetXaxis()->SetTitleOffset(1.36);
    tunfresult_rebinnedA->GetYaxis()->SetTitleOffset(1.6);
    tunfresult_rebinnedA->Draw();

    Draw2DLabels(fVariableNames[i]);
    DrawCMSLabels(0.04, fLumi);
    if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

    c1A->SaveAs(Form("%s/UnfoldedResults_rebinnedA_%s.pdf",PLOTFILE.Data(),fVariableNames[i].Data()));
    c1A->SaveAs(Form("%s/UnfoldedResults_rebinnedA_%s.root",PLOTFILE.Data(),fVariableNames[i].Data()));
    c1A->SaveAs(Form("%s/UnfoldedResults_rebinnedA_%s.C",PLOTFILE.Data(),fVariableNames[i].Data()));

    TCanvas *c1B = new TCanvas("c1B", "c1B");
    TH1D* tunfresult_rebinnedB = (TH1D*)unffile->Get(fVariableNames[i]+"TUnfResult_rebinnedB");
    tunfresult_rebinnedB->SetLineColor(kRed);      
    tunfresult_rebinnedB->SetLineWidth(2);      
    tunfresult_rebinnedB->SetMarkerColor(kRed);      
    tunfresult_rebinnedB->SetMarkerStyle(20);      
    tunfresult_rebinnedB->SetMarkerSize(0.5);
    tunfresult_rebinnedB->SetTitle("");      
    tunfresult_rebinnedB->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    tunfresult_rebinnedB->GetXaxis()->SetNdivisions(-20604);
    if(fVariableNames[i].Contains("_mttbar")) {
      tunfresult_rebinnedB->GetXaxis()->SetNdivisions(-20604);
      c1B->SetGridx();
    }
    SetAxis(fVariableNames[i], tunfresult_rebinnedB->GetXaxis());
    tunfresult_rebinnedB->SetMaximum(tunfresult_rebinnedB->GetMaximum()*1.05);
    tunfresult_rebinnedB->SetMinimum(tunfresult_rebinnedB->GetMinimum()*0.95);
    tunfresult_rebinnedB->GetYaxis()->SetTitle("Unfolded (uncorrected) number of events");
    tunfresult_rebinnedB->GetXaxis()->SetTitleOffset(1.36);
    tunfresult_rebinnedB->GetYaxis()->SetTitleOffset(1.6);
    tunfresult_rebinnedB->Draw();
    Draw2DLabels(fVariableNames[i]);
    DrawCMSLabels(0.04, fLumi);
    if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);
    c1B->SaveAs(Form("%s/UnfoldedResults_rebinnedB_%s.pdf",PLOTFILE.Data(),fVariableNames[i].Data()));

    TCanvas *c1A_Afb = new TCanvas("c1A_Afb", "c1A_Afb");
    TH1D* Afbresult_rebinnedA = (TH1D*)unffile->Get(fVariableNames[i]+"AfbResult_rebinnedA");
    Afbresult_rebinnedA->SetLineColor(kRed);      
    Afbresult_rebinnedA->SetLineWidth(2);      
    Afbresult_rebinnedA->SetMarkerColor(kRed);      
    Afbresult_rebinnedA->SetMarkerStyle(20);      
    Afbresult_rebinnedA->SetMarkerSize(0.5);
    Afbresult_rebinnedA->SetTitle("");      
    Afbresult_rebinnedA->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, A_{total} (right)",quantity.Data(),quantitycoef.Data()));
    Afbresult_rebinnedA->GetYaxis()->SetTitle("Unfolded asymmetries and coefficients");
    Afbresult_rebinnedA->GetXaxis()->SetTitleOffset(1.36);
    Afbresult_rebinnedA->GetYaxis()->SetTitleOffset(1.6);
    Afbresult_rebinnedA->Draw();

    Draw2DLabels(fVariableNames[i]);
    DrawCMSLabels(0.04, fLumi);
    if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

    c1A_Afb->SaveAs(Form("%s/UnfoldedAfbResults_rebinnedA_%s.pdf",PLOTFILE.Data(),fVariableNames[i].Data()));

    TCanvas *c1B_Afb = new TCanvas("c1B_Afb", "c1B_Afb");
    TH1D* Afbresult_rebinnedB = (TH1D*)unffile->Get(fVariableNames[i]+"AfbResult_rebinnedB");
    Afbresult_rebinnedB->SetLineColor(kRed);      
    Afbresult_rebinnedB->SetLineWidth(2);      
    Afbresult_rebinnedB->SetMarkerColor(kRed);      
    Afbresult_rebinnedB->SetMarkerStyle(20);      
    Afbresult_rebinnedB->SetMarkerSize(0.5);
    Afbresult_rebinnedB->SetTitle("");      
    Afbresult_rebinnedB->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, A_{total} (right)",quantity.Data(),quantitycoef.Data()));
    Afbresult_rebinnedB->GetYaxis()->SetTitle("Unfolded asymmetries and coefficients");
    Afbresult_rebinnedB->GetXaxis()->SetTitleOffset(1.36);
    Afbresult_rebinnedB->GetYaxis()->SetTitleOffset(1.6);
    Afbresult_rebinnedB->Draw();

    Draw2DLabels(fVariableNames[i]);
    DrawCMSLabels(0.04, fLumi);
    if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

    c1B_Afb->SaveAs(Form("%s/UnfoldedAfbResults_rebinnedB_%s.pdf",PLOTFILE.Data(),fVariableNames[i].Data()));


    TCanvas *c3 = new TCanvas("c3", "c3");
    TH1D* inputfoldedback = (TH1D*)unffile->Get(fVariableNames[i]+"InputFoldedBack_rebinnedA");
    inputfoldedback->SetLineColor(kRed);      
    inputfoldedback->SetLineWidth(2);      
    inputfoldedback->SetMarkerColor(kRed);      
    inputfoldedback->SetMarkerStyle(20);      
    inputfoldedback->SetMarkerSize(0.5);
    inputfoldedback->SetTitle("");      
    inputfoldedback->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    inputfoldedback->GetXaxis()->SetNdivisions(-20604);
    if(fVariableNames[i].Contains("_mttbar")) {
      inputfoldedback->GetXaxis()->SetNdivisions(-20608);
      c3->SetGridx();
    }
    SetAxis(fVariableNames[i], inputfoldedback->GetXaxis());
    inputfoldedback->SetMinimum(0);
    inputfoldedback->SetMaximum(inputfoldedback->GetMaximum()*1.4);
    inputfoldedback->GetYaxis()->SetTitle("Number of events");
    inputfoldedback->GetXaxis()->SetTitleOffset(1.36);
    inputfoldedback->GetYaxis()->SetTitleOffset(1.6);
    inputfoldedback->Draw();

    // TODO: Fix me
    // TH1D* aequimeasvec = (TH1D*)unffile->Get(fVariableNames[i]+"OrigInput_rebinnedB");
    // aequimeasvec->SetLineColor(kBlack);      
    // aequimeasvec->SetLineWidth(1);      
    // aequimeasvec->SetMarkerColor(kBlack);      
    // aequimeasvec->SetMarkerStyle(5);      
    // aequimeasvec->SetTitle("");      
    // aequimeasvec->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    // aequimeasvec->GetYaxis()->SetTitle("Number of events");
    // aequimeasvec->GetXaxis()->SetTitleOffset(1.36);
    // aequimeasvec->GetYaxis()->SetTitleOffset(1.6);
    // aequimeasvec->Draw("same");

    TLegend *l3 = new TLegend(0.8-gStyle->GetPadRightMargin(),0.875-gStyle->GetPadTopMargin(),1.0-gStyle->GetPadRightMargin(),1.0-gStyle->GetPadTopMargin());
    l3->AddEntry(inputfoldedback,"Folded Back");
    // l3->AddEntry(aequimeasvec, "Measured");
    l3->SetFillColor(kWhite);
    l3->SetLineColor(kBlack);
    l3->Draw("same");

    Draw2DLabels(fVariableNames[i]);
    DrawCMSLabels(0.04, fLumi);
    if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

    c3->SaveAs(Form("%s/MeasuredVsFoldedBack_%s.pdf" ,PLOTFILE.Data(),fVariableNames[i].Data()));
    c3->SaveAs(Form("%s/MeasuredVsFoldedBack_%s.root",PLOTFILE.Data(),fVariableNames[i].Data()));
    c3->SaveAs(Form("%s/MeasuredVsFoldedBack_%s.C"   ,PLOTFILE.Data(),fVariableNames[i].Data()));

    TCanvas *c4 = new TCanvas("c4", "c4");
    TGraph* rhoavgcurve = (TGraph*)unffile->Get(fVariableNames[i]+" ScanOfRhoAvg2");
    rhoavgcurve->SetLineColor(kBlack);      
    rhoavgcurve->SetLineWidth(2);      
    rhoavgcurve->SetMarkerColor(kBlack);      
    rhoavgcurve->SetMarkerStyle(21);      

    TH1F* favg = new TH1F("favg", "favg", 100, TMath::MinElement(rhoavgcurve->GetN(),rhoavgcurve->GetX()), TMath::MaxElement(rhoavgcurve->GetN(),rhoavgcurve->GetX()));
    favg->SetTitle("");
    favg->GetXaxis()->SetTitle("Log10(#tau)");
    favg->GetYaxis()->SetTitle("#rho_{avg}");
    favg->GetYaxis()->SetTitleOffset(1.0);
    favg->GetYaxis()->SetRangeUser(TMath::MinElement(rhoavgcurve->GetN(),rhoavgcurve->GetY())-0.05,TMath::MaxElement(rhoavgcurve->GetN(),rhoavgcurve->GetY())+0.05);

    rhoavgcurve->SetHistogram(favg);
    rhoavgcurve->Draw("AL");

    TGraph* rhoavg =(TGraph*)unffile->Get(fVariableNames[i]+"RhoAvg");
    rhoavg->SetLineColor(kRed);      
    rhoavg->SetLineWidth(2);      
    rhoavg->SetMarkerColor(kRed);      
    rhoavg->SetMarkerStyle(20);      
    rhoavg->Draw("P same");

    c4->SaveAs(Form("%s/AvgRho_%s.pdf",PLOTFILE.Data(),fVariableNames[i].Data()));
    c4->SaveAs(Form("%s/AvgRho_%s.root",PLOTFILE.Data(),fVariableNames[i].Data()));
    c4->SaveAs(Form("%s/AvgRho_%s.C",PLOTFILE.Data(),fVariableNames[i].Data()));

    std::ofstream intau(Form("%s/OptimalTauAndAvgRho_%s.txt",PLOTFILE.Data(),fVariableNames[i].Data()));
    intau <<fVariableNames[i].Data() << std::endl;
    intau <<"Optimal Tau: "<<std::pow(10,rhoavg->GetX()[0]) <<std::endl;
    intau <<"Average Rho: "<<rhoavg->GetY()[0] <<std::endl;
    intau.close();

    // Amandeep : : For consistency in style with other plots.
    // gStyle->SetPalette(1);
    // End
    gStyle->SetPadRightMargin(0.13);
    gStyle->SetPalette(1);

    TCanvas *c5 = new TCanvas("c5", "c5");
    TH2D* resmat = (TH2D*)unffile->Get(fVariableNames[i]+"ResMat");
    resmat->SetTitle("");
    resmat->GetXaxis()->SetTitle(Form("%s (Gen)",quantity.Data()));
    resmat->GetYaxis()->SetTitle(Form("%s (Reco)",quantity.Data()));
    resmat->GetXaxis()->SetTitleOffset(1.45);
    resmat->GetYaxis()->SetTitleOffset(1.55);
    resmat->GetXaxis()->SetNdivisions(-20604);
    resmat->GetYaxis()->SetNdivisions(-20604);
    if(fVariableNames[i].Contains("_mttbar")) {
      resmat->GetXaxis()->SetNdivisions(-20608);
      resmat->GetYaxis()->SetNdivisions(-20616);
      resmat->GetXaxis()->SetLabelSize(0);
      resmat->GetYaxis()->SetLabelSize(0);
      c5->SetGridx();
      c5->SetGridy();
    }

    SetAxis(fVariableNames[i], resmat->GetXaxis());
    SetAxis(fVariableNames[i], resmat->GetYaxis());
    gPad->SetLogz();
    resmat->Draw("colz");

    DrawResMat2DAxisLabels(fVariableNames[i]);
    Draw2DLabels(fVariableNames[i], 0.85, 0.019);
    Draw2DLabelsVertical(fVariableNames[i], 0.82, 0.019);
    DrawCMSLabels(0.04, fLumi);
    if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

    c5->SaveAs(Form("%s/ResponseMatrix_%s.pdf",PLOTFILE.Data(),fVariableNames[i].Data()));
    c5->SaveAs(Form("%s/ResponseMatrix_%s.root",PLOTFILE.Data(),fVariableNames[i].Data()));
    c5->SaveAs(Form("%s/ResponseMatrix_%s.C",PLOTFILE.Data(),fVariableNames[i].Data()));

    /*
    TCanvas *c6 = new TCanvas("c6", "c6");
    TH2D* aequiresmat = (TH2D*)unffile->Get(fVariableNames[i]+"AequiResMat");
    aequiresmat->SetTitle("");
    aequiresmat->GetXaxis()->SetTitle(Form("Equidistance bins %s (Gen)",quantity.Data()));
    aequiresmat->GetYaxis()->SetTitle(Form("Equidistance bins %s (Reco)",quantity.Data()));
    aequiresmat->Draw("colz");
    DrawCMSLabels(0.04);
    if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

    c6->SaveAs(Form("%s/EquidistantResMat_%s.pdf",PLOTFILE.Data(),fVariableNames[i].Data()));
    */

    gStyle->SetPadRightMargin(0.05);

    // Amandeep : Commenting to try and fix combined channel plotting
    TCanvas *c7 = new TCanvas("c7", "c7");
    TH1D* gen_norm = (TH1D*)unffile->Get(fVariableNames[i]+"Gen_norm");
    gen_norm->SetLineColor(kRed);      
    gen_norm->SetLineWidth(2);      
    gen_norm->SetMarkerColor(kRed);      
    gen_norm->SetMarkerStyle(20);      
    gen_norm->SetMarkerSize(0.5);
    gen_norm->SetTitle("");      
    gen_norm->GetXaxis()->SetTitle(Form("%s (Gen)",quantity.Data()));
    gen_norm->GetXaxis()->SetNdivisions(-20604);
    if(fVariableNames[i].Contains("_mttbar")) {
      gen_norm->GetXaxis()->SetNdivisions(-20608);
      c7->SetGridx();
    }
    SetAxis(fVariableNames[i], gen_norm->GetXaxis());
    gen_norm->SetMinimum(0);
    gen_norm->SetMaximum(gen_norm->GetMaximum()*1.4);
    gen_norm->GetYaxis()->SetTitle("Number of events");
    gen_norm->GetXaxis()->SetTitleOffset(1.36);
    gen_norm->GetYaxis()->SetTitleOffset(1.6);
    gen_norm->Draw();

    Draw2DLabels(fVariableNames[i]);
    DrawCMSLabels(0.04, fLumi);
    if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

    c7->SaveAs(Form("%s/Gen_%s.pdf",PLOTFILE.Data(),fVariableNames[i].Data()));
    c7->SaveAs(Form("%s/Gen_%s.root",PLOTFILE.Data(),fVariableNames[i].Data()));
    c7->SaveAs(Form("%s/Gen_%s.C",PLOTFILE.Data(),fVariableNames[i].Data()));

    TCanvas *c8 = new TCanvas("c8", "c8");
    TH1D* visgen_norm = (TH1D*)unffile->Get(fVariableNames[i]+"VisGen_norm");
    visgen_norm->SetLineColor(kRed);      
    visgen_norm->SetLineWidth(2);      
    visgen_norm->SetMarkerColor(kRed);      
    visgen_norm->SetMarkerStyle(20);      
    visgen_norm->SetMarkerSize(0.5);
    visgen_norm->SetTitle("");      
    visgen_norm->GetXaxis()->SetTitle(Form("%s (VisGen)",quantity.Data()));
    visgen_norm->GetXaxis()->SetNdivisions(-20604);
    if(fVariableNames[i].Contains("_mttbar")) {
      visgen_norm->GetXaxis()->SetNdivisions(-20608);
      c8->SetGridx();
    }
    SetAxis(fVariableNames[i], visgen_norm->GetXaxis());
    visgen_norm->SetMinimum(0);
    visgen_norm->SetMaximum(visgen_norm->GetMaximum()*1.4);
    visgen_norm->GetYaxis()->SetTitle("Number of events");
    visgen_norm->GetXaxis()->SetTitleOffset(1.36);
    visgen_norm->GetYaxis()->SetTitleOffset(1.6);
    visgen_norm->Draw();

    Draw2DLabels(fVariableNames[i]);
    DrawCMSLabels(0.04);
    if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

    c8->SaveAs(Form("%s/VisGen_%s.pdf",PLOTFILE.Data(),fVariableNames[i].Data()));

    TCanvas *c9 = new TCanvas("c9", "c9");
    TH1D* reco_norm = (TH1D*)unffile->Get(fVariableNames[i]+"Reco_norm");
    reco_norm->SetLineColor(kBlack);      
    reco_norm->SetLineWidth(2);      
    reco_norm->SetMarkerColor(kBlack);      
    reco_norm->SetMarkerStyle(21);      
    reco_norm->SetMarkerSize(0.5);
    reco_norm->SetTitle("");      
    reco_norm->GetXaxis()->SetTitle(Form("%s (Reco)",quantity.Data()));
    reco_norm->GetXaxis()->SetNdivisions(-20604);
    if(fVariableNames[i].Contains("_mttbar")) {
      reco_norm->GetXaxis()->SetNdivisions(-20616);
      reco_norm->GetXaxis()->SetLabelSize(0.02);
      c9->SetGridx();
    } 
    SetAxis(fVariableNames[i], reco_norm->GetXaxis());
    reco_norm->SetMinimum(0);
    reco_norm->SetMaximum(reco_norm->GetMaximum()*1.4);
    reco_norm->GetYaxis()->SetTitle("Number of events");
    reco_norm->GetXaxis()->SetTitleOffset(1.36);
    reco_norm->GetYaxis()->SetTitleOffset(1.6);
    reco_norm->Draw("");

    Draw2DLabels(fVariableNames[i]);
    DrawCMSLabels(0.04, fLumi);
    if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

    c9->SaveAs(Form("%s/Reco_%s.pdf",PLOTFILE.Data(),fVariableNames[i].Data()));
    c9->SaveAs(Form("%s/Reco_%s.root",PLOTFILE.Data(),fVariableNames[i].Data()));
    c9->SaveAs(Form("%s/Reco_%s.C",PLOTFILE.Data(),fVariableNames[i].Data()));

    TCanvas *c10 = new TCanvas("c10", "c10");
    TH1D* purity = (TH1D*)unffile->Get(fVariableNames[i]+"MigPurFromMatHist");
    purity->SetLineColor(kGreen+2);      
    purity->SetLineWidth(2);      
    purity->SetMarkerColor(kGreen+2);      
    purity->SetMarkerStyle(22);      
    purity->SetMarkerSize(0.5);
    purity->SetTitle("");      
    purity->GetXaxis()->SetTitle(Form("%s (Reco)",quantity.Data()));
    purity->GetYaxis()->SetTitle("Purity");
    purity->GetXaxis()->SetTitleOffset(1.36);
    purity->GetYaxis()->SetTitleOffset(1.6);

    TCanvas *c11 = new TCanvas("c11", "c11");
    TH1D* stability = (TH1D*)unffile->Get(fVariableNames[i]+"StabFromMatHist");
    stability->SetLineColor(kBlue-4);      
    stability->SetLineWidth(2);      
    stability->SetMarkerColor(kBlue-4);      
    stability->SetMarkerStyle(23);      
    stability->SetMarkerSize(0.5);
    stability->SetTitle("");      
    stability->GetXaxis()->SetTitle(Form("%s (Gen)",quantity.Data()));
    stability->GetYaxis()->SetTitle("Stability");
    stability->GetXaxis()->SetTitleOffset(1.36);
    stability->GetYaxis()->SetTitleOffset(1.6);

    TCanvas *c12 = new TCanvas("c12", "c12");
    TH1D* toteff = (TH1D*)unffile->Get(fVariableNames[i]+"TotEffiHist");
    toteff->SetLineColor(kBlack);      
    toteff->SetLineWidth(2);      
    toteff->SetMarkerColor(kBlack);      
    toteff->SetMarkerStyle(21);      
    toteff->SetMarkerSize(0.5);
    toteff->SetTitle("");      
    toteff->GetXaxis()->SetTitle(Form("%s (Gen)",quantity.Data()));
    toteff->GetXaxis()->SetNdivisions(-20604);
    if(fVariableNames[i].Contains("_mttbar")) {
      toteff->GetXaxis()->SetNdivisions(-20608);
      toteff->GetXaxis()->SetLabelSize(0);
      c12->SetGridx();
    }
    SetAxis(fVariableNames[i], toteff->GetXaxis());
    toteff->SetMinimum(0);
    toteff->SetMaximum(toteff->GetMaximum()*1.4);
    toteff->GetYaxis()->SetTitle("Efficiency X Acceptance");
    toteff->GetXaxis()->SetTitleOffset(1.36);
    toteff->GetYaxis()->SetTitleOffset(1.6);
    toteff->Draw("");

    DrawTotEff2DAxisLabels(fVariableNames[i]);
    Draw2DLabels(fVariableNames[i]);
    DrawCMSLabels(0.04, fLumi);
    if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

    c12->SaveAs(Form("%s/TotEff_%s.pdf" ,PLOTFILE.Data(),fVariableNames[i].Data()));
    c12->SaveAs(Form("%s/TotEff_%s.root",PLOTFILE.Data(),fVariableNames[i].Data()));
    c12->SaveAs(Form("%s/TotEff_%s.C"   ,PLOTFILE.Data(),fVariableNames[i].Data()));

    TCanvas *c13 = new TCanvas("c13", "c13");
    TH1D* eff = (TH1D*)unffile->Get(fVariableNames[i]+"EffiHist");
    eff->SetLineColor(kBlack);      
    eff->SetLineWidth(2);      
    eff->SetMarkerColor(kBlack);      
    eff->SetMarkerStyle(21);      
    eff->SetMarkerSize(0.5);
    eff->SetTitle("");      
    eff->GetXaxis()->SetTitle(Form("%s (Gen)",quantity.Data()));
    eff->GetXaxis()->SetNdivisions(-20604);
    if(fVariableNames[i].Contains("_mttbar")) {
      eff->GetXaxis()->SetNdivisions(-20608);
      c13->SetGridx();
    }
    SetAxis(fVariableNames[i], eff->GetXaxis());
    eff->SetMinimum(0);
    eff->SetMaximum(eff->GetMaximum()*1.4);
    eff->GetYaxis()->SetTitle("Efficiency");
    eff->GetXaxis()->SetTitleOffset(1.36);
    eff->GetYaxis()->SetTitleOffset(1.6);
    eff->Draw("");

    Draw2DLabels(fVariableNames[i]);
    DrawCMSLabels(0.04, fLumi);
    if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

    c13->SaveAs(Form("%s/Eff_%s.pdf",PLOTFILE.Data(),fVariableNames[i].Data()));

    TCanvas *c14 = new TCanvas("c14", "c14");
    TH1D* accep = (TH1D*)unffile->Get(fVariableNames[i]+"AccepHist");
    accep->SetLineColor(kBlack);      
    accep->SetLineWidth(2);      
    accep->SetMarkerColor(kBlack);      
    accep->SetMarkerStyle(21);      
    accep->SetMarkerSize(0.5);
    accep->SetTitle("");      
    accep->GetXaxis()->SetTitle(Form("%s (Gen)",quantity.Data()));
    accep->GetXaxis()->SetNdivisions(-20604);
    if(fVariableNames[i].Contains("_mttbar")) {
      accep->GetXaxis()->SetNdivisions(-20608);
      c14->SetGridx();
    }
    SetAxis(fVariableNames[i], accep->GetXaxis());
    accep->SetMinimum(0);
    accep->SetMaximum(accep->GetMaximum()*1.4);
    accep->GetYaxis()->SetTitle("Acceptance");
    accep->GetXaxis()->SetTitleOffset(1.36);
    accep->GetYaxis()->SetTitleOffset(1.6);
    accep->Draw("");

    Draw2DLabels(fVariableNames[i]);
    DrawCMSLabels(0.04, fLumi);
    if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

    c14->SaveAs(Form("%s/Accep_%s.pdf",PLOTFILE.Data(),fVariableNames[i].Data()));

    TCanvas *c14comb = new TCanvas("c14comb", "c14comb");
    purity->GetYaxis()->SetTitle("");    
    stability->GetYaxis()->SetTitle("");    
    purity->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    purity->GetYaxis()->SetRangeUser(0, 0.7);    
    purity->GetXaxis()->SetNdivisions(-20604);
    stability->GetXaxis()->SetNdivisions(-20604);
    if(fVariableNames[i].Contains("_mttbar")) {
      purity->GetXaxis()->SetNdivisions(-20604);
      stability->GetXaxis()->SetNdivisions(-20604);
      purity->GetXaxis()->SetLabelSize(0);
      stability->GetXaxis()->SetLabelSize(0);
      c14comb->SetGridx();
    }
    SetAxis(fVariableNames[i], purity->GetXaxis());
    SetAxis(fVariableNames[i], stability->GetXaxis());
    purity->SetMinimum(0);
    purity->SetMaximum(purity->GetMaximum()*1.4);
    stability->SetMinimum(0);
    stability->SetMaximum(stability->GetMaximum()*1.4);
    if(MatrixUnf::isOtherSymType(fVariableNames[i])) purity->GetYaxis()->SetRangeUser(0.9, 1.1);    
    purity->Draw("");
    stability->Draw("same");

    TLegend *l14comb=new TLegend(0.8-gStyle->GetPadRightMargin(),0.875-gStyle->GetPadTopMargin(),1.0-gStyle->GetPadRightMargin(),1.0-gStyle->GetPadTopMargin());
    l14comb->AddEntry(purity, "Purity");
    l14comb->AddEntry(stability, "Stability");
    l14comb->SetFillColor(kWhite);
    l14comb->SetLineColor(kBlack);
    l14comb->Draw("same");

    DrawPurStab2DAxisLabels(fVariableNames[i]);
    Draw2DLabels(fVariableNames[i]);
    DrawCMSLabels(0.04, fLumi);
    if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

    c14comb->SaveAs(Form("%s/PurStab_%s.pdf",PLOTFILE.Data(),fVariableNames[i].Data()));
    c14comb->SaveAs(Form("%s/PurStab_%s.root",PLOTFILE.Data(),fVariableNames[i].Data()));
    c14comb->SaveAs(Form("%s/PurStab_%s.C",PLOTFILE.Data(),fVariableNames[i].Data()));


    if(!fsetSyst) continue;

    // *****************************************************************
    // Calculate systematic and combined (stat+syst) covariance matrices
    // *****************************************************************

    std::vector<TMatrixD> TotalCovMats;
    std::vector<TMatrixD> StatCovMats;

    TString NameTag = "_rebinnedA";

    for (int i_tag = 0; i_tag < 2; ++i_tag)
    {
      if(i_tag==1) NameTag = "_rebinnedB";

      TH2D* htotalsystcovmat[2];
      for (int inorm = 0; inorm <= 1; ++inorm)
      {
        // Amandeep : For consistency in style with other plots.
        // gStyle->SetPalette(1);
        // End
        gStyle->SetPadLeftMargin(0.11);
        gStyle->SetPadRightMargin(0.13);

        TCanvas *c15 = new TCanvas("c15", "c15");
        TH2D* hstatcovmat = (TH2D*)unffile->Get(fVariableNames[i]+"Ematrix"+(inorm?"Norm":"Cor")+NameTag);
        hstatcovmat->SetTitle("");
        hstatcovmat->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
        hstatcovmat->GetYaxis()->SetTitle(Form("%s",quantity.Data()));
	hstatcovmat->GetXaxis()->SetTitleSize(0.035);
	hstatcovmat->GetXaxis()->SetTitleOffset(1.65);
	hstatcovmat->GetYaxis()->SetTitleSize(0.035);
	hstatcovmat->GetYaxis()->SetTitleOffset(1.65);

        hstatcovmat->GetXaxis()->SetNdivisions(-20604);
        hstatcovmat->GetYaxis()->SetNdivisions(-20604);
	hstatcovmat->GetZaxis()->SetLabelSize(0.03);

        if(fVariableNames[i].Contains("_mttbar")) {
          hstatcovmat->GetXaxis()->SetNdivisions(-20604);
          hstatcovmat->GetYaxis()->SetNdivisions(-20604);
          hstatcovmat->GetXaxis()->SetLabelSize(0);
          hstatcovmat->GetYaxis()->SetLabelSize(0);
          c15->SetGridx();
          c15->SetGridy();
        }

        SetAxis(fVariableNames[i], hstatcovmat->GetXaxis());
        SetAxis(fVariableNames[i], hstatcovmat->GetYaxis());
        hstatcovmat->Draw("colz");

	DrawCov2DAxisLabels(fVariableNames[i]);
        Draw2DLabels(fVariableNames[i], 0.85, 0.019);
        Draw2DLabelsVertical(fVariableNames[i], 0.82, 0.019);
        DrawCMSLabels(0.04, fLumi);
        if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

        c15->SaveAs(Form("%s/StatCovMatrix%s%s_%s.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data(),fVariableNames[i].Data()));
        c15->SaveAs(Form("%s/StatCovMatrix%s%s_%s.C",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data(),fVariableNames[i].Data()));

        gStyle->SetPadRightMargin(0.05);


        const int nbinsrhoi = hstatcovmat->GetNbinsX();

        if(i==0) {
          hdeltasquaredStat_AllVar[i_tag*2+inorm] = new TH1D(Form("hdeltasquaredStat_AllVar%s%s",(inorm?"Norm":""),NameTag.Data()), Form("hdeltasquaredStat_AllVar%s%s",(inorm?"Norm":""),NameTag.Data()), fVariableNames.size()*nbinsrhoi,0,fVariableNames.size()*nbinsrhoi);
          hdeltasquaredStat_AllVar[i_tag*2+inorm]->SetDirectory(0);
        }

        // Amandeep : Why is the MC stat uncertainty labelled RespMatSys_Uncorr_CovMatSys ? 
        TH2D* hMCstatcovmat = (TH2D*)unffile->Get(fVariableNames[i]+"RespMatSys_Uncorr_CovMatSys"+(inorm?"Norm":"")+NameTag);
        TH2D* hBkgstatcovmat[backgroundnames.size()];
        TH2D* hBkgstatcovmatTotal;

        TMatrixD TotalCovMat         = MatrixUnf::TH2toTMatrixD(hstatcovmat);
        TMatrixD DataStatOnlyCorrMat = MatrixUnf::TH2toTMatrixD(hstatcovmat);
        TMatrixD TotalStatCorrMat    = MatrixUnf::TH2toTMatrixD(hstatcovmat);
        
        // Amandeep : Disabling to gauge the impact of MC stat cov matrix
        TotalCovMat      += MatrixUnf::TH2toTMatrixD(hMCstatcovmat);
        TotalStatCorrMat += MatrixUnf::TH2toTMatrixD(hMCstatcovmat);
        // End

        for (int i_bin = 1; i_bin <= nbinsrhoi; ++i_bin){
          // Amandeep : Disabling to gauge the impact of MC stat cov matrix
          // Original :
          hdeltasquaredStat_AllVar[i_tag*2+inorm]->SetBinContent(i_bin+i*nbinsrhoi, ( hstatcovmat->GetBinContent(i_bin,i_bin) + hMCstatcovmat->GetBinContent(i_bin,i_bin) ) );
          // Modified :
          // hdeltasquaredStat_AllVar[i_tag*2+inorm]->SetBinContent(i_bin+i*nbinsrhoi, ( hstatcovmat->GetBinContent(i_bin,i_bin) ) );
          // End
        }

        // Loop over all the input systematics
        if(fDoSubtrBg) {
          for (unsigned int ibkg = 0; ibkg < backgroundnames.size(); ++ibkg)
          {
            hBkgstatcovmat[ibkg] = (TH2D*)unffile->Get(fVariableNames[i]+"BkgSysUncorr_"+backgroundnames.at(ibkg)+"_CovMatSys"+(inorm?"Norm":"")+NameTag);
            if(ibkg==0) hBkgstatcovmatTotal = (TH2D*) hBkgstatcovmat[ibkg]->Clone(hBkgstatcovmat[ibkg]->GetName()+TString("Total"));
            else hBkgstatcovmatTotal->Add(hBkgstatcovmat[ibkg]);
            DataStatOnlyCorrMat -= MatrixUnf::TH2toTMatrixD(hBkgstatcovmat[ibkg]);
          }

        }

        DataStatOnlyCorrMat.NormByDiag(TMatrixDDiag(DataStatOnlyCorrMat));
        TotalStatCorrMat.NormByDiag(TMatrixDDiag(TotalStatCorrMat));

        for (int i_bin = 1; i_bin <= nbinsrhoi; ++i_bin)
        {
          for (int j_bin = 1; j_bin <= nbinsrhoi; ++j_bin)
          {
            hTUnfoldDataStatCorrM_AllVar[i_tag*2+inorm]->SetBinContent(i_bin+i*nbinsrhoi, j_bin+i*nbinsrhoi, DataStatOnlyCorrMat(i_bin-1,j_bin-1) );
            hTUnfoldTotalStatCorrM_AllVar[i_tag*2+inorm]->SetBinContent(i_bin+i*nbinsrhoi, j_bin+i*nbinsrhoi, TotalStatCorrMat(i_bin-1,j_bin-1) );
          }
        }


        TCanvas *c16 = new TCanvas(Form("c16%s",(inorm?"Norm":"")), Form("c16%s",(inorm?"Norm":"")));
        //gPad->SetLogy();
        int nxchecksys = 2;
        const int nleg = 4;
        double ff[nleg];
        TLegend *l16[nleg];
        int nperleg     = (nsyst+1-nxchecksys)/nleg;
        int remainder   = (nsyst+1-nxchecksys)%nleg;
        double padwidth = 1.0 - gStyle->GetPadLeftMargin() - gStyle->GetPadRightMargin();

        for (int ileg = 0; ileg < nleg; ++ileg)
        {
          ff[ileg]  = double(nperleg + 0.2 + (remainder>ileg?1:0) ) / double(nperleg + 0.2 + (remainder>0?1:0) );
          l16[ileg] = new TLegend(gStyle->GetPadLeftMargin()+ileg*padwidth/nleg + 0.02, 1.0 - gStyle->GetPadTopMargin() - ff[ileg]*0.5/nleg - 0.05, gStyle->GetPadLeftMargin()+(ileg+1)*padwidth/nleg - 0.02, 1.0 - gStyle->GetPadTopMargin() - 0.03);
        }

        TH1D* hsysttemp[nsyst];

        TH1D *htotsyst_reldelta = (TH1D*)((TH1D*)unffile->Get(fVariableNames[i]+"TUnfResultCor"+NameTag)->Clone(fVariableNames[i]+"_TotDeltaSyst"+(inorm?"Norm":"")+"Unfolded"+NameTag));
        htotsyst_reldelta->Reset();
        TH1D* hcombsystdelta[ combsystnames.size() ];

        std::vector<TH2D*>hcomponentsystcovmat;

        TH2D* hcovtemp[nsyst];
        TH2D* hcombsystcovmat[ combsystnames.size() ];

        int legcounter = 0;

        for(int syst = 0; syst < nsyst; syst++){

          // Amandeep : This where the systematics are computed
          hsysttemp[syst] = (TH1D*)(MatrixUnfControl::GetEnvelopeForUnfoldedResults(CHAN, tempVectorOfInputSystematics[syst], fVariableNames[i], unffile, inorm, NameTag));
          if(hsysttemp[syst] == NULL) hsysttemp[syst] = (TH1D*)(MatrixUnfControl::CalcUpDownDifference(CHAN, tempVectorOfInputSystematics[syst], fVariableNames[i], unffile, inorm, NameTag));

          // Using GetCombinedCovarianceMatrix function that calculates a fully consistent syst covariance matrix. 
          // Do this by calculating the correlation matrix for each source (from the sum in quadrature) 
          // then multiply each matrix element by sigma_i*sigma_j 
          // calculated from the GetEnvelopeForUnfoldedResults or CalcUpDownDifference.
          
          hcovtemp[syst] = (TH2D*) MatrixUnfControl::GetCombinedCovarianceMatrix(CHAN, tempVectorOfInputSystematics[syst], fVariableNames[i], unffile, hsysttemp[syst], inorm, NameTag);
          
          if( tempVectorOfInputSystematics[syst] != "AMCATNLOFXFX" && tempVectorOfInputSystematics[syst] != "POWHEGV2HERWIG"  ) hcomponentsystcovmat.push_back( hcovtemp[syst] );
          if(syst==0) {
            htotalsystcovmat[inorm] = (TH2D*) hcovtemp[syst]->Clone(Form("%s_Total",hcovtemp[syst]->GetName()));

            for (unsigned int i_cs  = 0; i_cs < combsystnames.size(); ++i_cs) {
              hcombsystcovmat[i_cs] = (TH2D*) hcovtemp[syst]->Clone(Form("%s_%s",hcovtemp[syst]->GetName(),combsystnames.at(i_cs).Data()));
              hcombsystcovmat[i_cs]->Reset();

              hcombsystdelta[i_cs] = (TH1D*) hsysttemp[syst]->Clone(Form("%s_%s",hsysttemp[syst]->GetName(),combsystnames.at(i_cs).Data()));
              hcombsystdelta[i_cs]->Reset();
            }
          }

          else if( tempVectorOfInputSystematics[syst] != "AMCATNLOFXFX" && tempVectorOfInputSystematics[syst] != "POWHEGV2HERWIG"  ) htotalsystcovmat[inorm]->Add( hcovtemp[syst] );

          for (unsigned int i_cs = 0; i_cs < combsystnames.size(); ++i_cs) {
            if( (TString(tempVectorOfInputSystematics_BKGnames[syst]).BeginsWith( combsystnames.at(i_cs).Data() ) )  ) hcombsystcovmat[i_cs]->Add( hcovtemp[syst] );
          }

          hsysttemp[syst]->SetLineColor(syst+1-9*(syst/9));
          hsysttemp[syst]->SetLineWidth(2);
          hsysttemp[syst]->SetMarkerStyle(1);
          hsysttemp[syst]->SetLineStyle(syst/9+1);
          //hsysttemp[syst]->GetYaxis()->SetRangeUser(0,0.2);
          hsysttemp[syst]->SetTitle("");
          hsysttemp[syst]->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
          hsysttemp[syst]->GetYaxis()->SetTitle("(Nominal-Syst)/Nominal");
          hsysttemp[syst]->GetXaxis()->SetTitleOffset(1.36);
          hsysttemp[syst]->GetYaxis()->SetTitleOffset(1.6);

          for(int ibin = 1; ibin<=htotsyst_reldelta->GetNbinsX(); ibin++){

            TString SystName = tempVectorOfInputSystematics[syst];
            double scalingfactor = 1.;
            // Amandeep : Jason said to change this to 1
            if(SystName=="MASS") scalingfactor=1.;               // convert from +/-3 GeV to +/-1 GeV
            if(SystName.Contains("PDF")) scalingfactor=10.;      // convert from sum in quadrature to RMS (for 100 replicas)
            if(SystName.Contains("PDF_ALPHAS")) scalingfactor=1.;
            //std::cout<<"test: "<<SystName<<" scalingfactor: "<<scalingfactor<<std::endl;

            if( tempVectorOfInputSystematics[syst] != "AMCATNLOFXFX" && tempVectorOfInputSystematics[syst] != "POWHEGV2HERWIG"  )  htotsyst_reldelta->SetBinContent(ibin, TMath::Sqrt(hsysttemp[syst]->GetBinContent(ibin)*hsysttemp[syst]->GetBinContent(ibin)/scalingfactor/scalingfactor + htotsyst_reldelta->GetBinContent(ibin)*htotsyst_reldelta->GetBinContent(ibin)));

            for (unsigned int i_cs = 0; i_cs < combsystnames.size(); ++i_cs) {
              if( (TString(tempVectorOfInputSystematics_BKGnames[syst]).BeginsWith( combsystnames.at(i_cs).Data() ) )  ) hcombsystdelta[i_cs]->SetBinContent(ibin, TMath::Sqrt(hsysttemp[syst]->GetBinContent(ibin)*hsysttemp[syst]->GetBinContent(ibin)/scalingfactor/scalingfactor + hcombsystdelta[i_cs]->GetBinContent(ibin)*hcombsystdelta[i_cs]->GetBinContent(ibin)));
            }


          }

          c16->cd();
          hsysttemp[syst]->GetXaxis()->SetNdivisions(-20604);
          htotsyst_reldelta->GetXaxis()->SetNdivisions(-20604);
          SetAxis(fVariableNames[i],  hsysttemp[syst]->GetXaxis());
          SetAxis(fVariableNames[i], htotsyst_reldelta->GetXaxis());

          if(syst==0) {
            hsysttemp[syst]->SetMinimum(0);
            hsysttemp[syst]->SetMaximum(hsysttemp[syst]->GetMaximum()*1.5);
	    if(fVariableNames[i].Contains("_mttbar")) {
	      hsysttemp[syst]->GetXaxis()->SetLabelSize(0);
	    }
            hsysttemp[syst]->Draw("hist");
          }

          else if( tempVectorOfInputSystematics[syst] != "AMCATNLOFXFX" && tempVectorOfInputSystematics[syst] != "POWHEGV2HERWIG"  ) hsysttemp[syst]->Draw("histsame");
          
          if( tempVectorOfInputSystematics[syst] != "AMCATNLOFXFX" && tempVectorOfInputSystematics[syst] != "POWHEGV2HERWIG"  ) {
            l16[legcounter%nleg]->AddEntry(hsysttemp[syst],tempVectorOfInputSystematics[syst].c_str());
            ++legcounter;
          }
          if(syst==(nsyst-1)) {
            htotsyst_reldelta->SetLineColor(12);
            htotsyst_reldelta->SetLineWidth(2);
            htotsyst_reldelta->SetLineStyle(1);
            htotsyst_reldelta->SetMarkerStyle(1);
            hsysttemp[0]->GetYaxis()->SetRangeUser(0,htotsyst_reldelta->GetMaximum()*1.29);
            htotsyst_reldelta->Draw("histsame");
            l16[(syst+1-nxchecksys)%nleg]->AddEntry(htotsyst_reldelta,"Total Syst");
          }

        }

        StatCovMats.push_back(TotalCovMat);
        TotalCovMat+=MatrixUnf::TH2toTMatrixD(htotalsystcovmat[inorm]);
        TotalCovMats.push_back(TotalCovMat);
            
        c16->cd();
        for (int ileg = 0; ileg < nleg; ++ileg)
        {
          // l16[ileg]->SetFillColor(kWhite);
	  l16[ileg]->SetBorderSize(0);
          // l16[ileg]->SetFillStyle(0);
          // l16[ileg]->SetLineColor(kBlack);
          // l16[ileg]->SetLineWidth(1);
          l16[ileg]->Draw("same");
        }

	DrawPurStab2DAxisLabels(fVariableNames[i]);
	Draw2DLabels(fVariableNames[i], 0.65);
        DrawCMSLabels(0.04, fLumi);
        if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

        c16->SaveAs(Form("%s/deltaSyst%s%s_%s.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data(),fVariableNames[i].Data()));
        hsysttemp[0]->GetYaxis()->SetRangeUser(1e-4,1);
	hsysttemp[0]->SetMinimum(.00001);
	hsysttemp[0]->SetMaximum(100);
        gPad->SetLogy();
        c16->SaveAs(Form("%s/deltaSystlog%s%s_%s.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data(),fVariableNames[i].Data()));


        TCanvas *c16comb = new TCanvas(Form("c16comb%s",(inorm?"Norm":"")), Form("c16comb%s",(inorm?"Norm":"")));
        if(inorm) nxchecksys=2;
        else nxchecksys=0;
        TLegend *l16comb[nleg];
        nperleg = (combsystnames.size()+1-nxchecksys)/nleg;
        remainder = (combsystnames.size()+1-nxchecksys)%nleg;

        for (int ileg = 0; ileg < nleg; ++ileg)
        {
          ff[ileg] = double(nperleg + 0.2 + (remainder>ileg?1:0) ) / double(nperleg + 0.2 + (remainder>0?1:0) );
          l16comb[ileg]=new TLegend(gStyle->GetPadLeftMargin()+ileg*padwidth/nleg + 0.02, 1.0 - gStyle->GetPadTopMargin() - ff[ileg]*0.5/nleg - 0.05, gStyle->GetPadLeftMargin()+(ileg+1)*padwidth/nleg - 0.02, 1.0 - gStyle->GetPadTopMargin() - 0.03);
        }

        c16comb->cd();
        legcounter=0;

        for (unsigned int i_cs = 0; i_cs < combsystnames.size(); ++i_cs) {

          if( !inorm || (combsystnames[i_cs] != "Lumi" && combsystnames[i_cs] != "BR") ) { 

            hcombsystdelta[i_cs]->SetLineColor(i_cs+1-9*(i_cs/9));
            hcombsystdelta[i_cs]->SetLineWidth(2);
            hcombsystdelta[i_cs]->SetMarkerStyle(1);
            hcombsystdelta[i_cs]->SetLineStyle(i_cs/9+1);
            hcombsystdelta[i_cs]->SetTitle("");
	    hcombsystdelta[i_cs]->GetXaxis()->SetNdivisions(-20604);
	    SetAxis(fVariableNames[i], hcombsystdelta[i_cs]->GetXaxis());
            hcombsystdelta[i_cs]->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
            hcombsystdelta[i_cs]->GetYaxis()->SetTitle("(Nominal-Syst)/Nominal");
            hcombsystdelta[i_cs]->GetXaxis()->SetTitleOffset(1.36);
            hcombsystdelta[i_cs]->GetYaxis()->SetTitleOffset(1.6);

	    if(i_cs==0 && fVariableNames[i].Contains("_mttbar")) {
	      hcombsystdelta[i_cs]->GetXaxis()->SetLabelSize(0);
	    }

            if(i_cs==0) hcombsystdelta[i_cs]->Draw("hist");
            else hcombsystdelta[i_cs]->Draw("histsame");

            l16comb[legcounter%nleg]->AddEntry(hcombsystdelta[i_cs],combsystnames[i_cs].Data());
            ++legcounter;

          }
        }

        hcombsystdelta[0]->GetYaxis()->SetRangeUser(0,htotsyst_reldelta->GetMaximum()*1.4);
        htotsyst_reldelta->Draw("histsame");
        l16comb[(combsystnames.size()-nxchecksys)%nleg]->AddEntry(htotsyst_reldelta,"Total Syst");

        for (int ileg = 0; ileg < nleg; ++ileg)
        {
          // l16comb[ileg]->SetFillColor(kWhite);
	  l16comb[ileg]->SetBorderSize(0);
          // l16comb[ileg]->SetFillStyle(0);
          // l16comb[ileg]->SetLineColor(kBlack);
          // l16comb[ileg]->SetLineWidth(1);

          l16comb[ileg]->SetTextFont(42);
          l16comb[ileg]->SetTextAlign(12);
          l16comb[ileg]->SetTextSize(0.0125);

          l16comb[ileg]->Draw("same");
        }

	if(fVariableNames[i].Contains("_mttbar")) {
	  c16comb->SetGridx();
	}

	DrawPurStab2DAxisLabels(fVariableNames[i]);
        Draw2DLabels(fVariableNames[i], 0.65);
        DrawCMSLabels(0.04, fLumi);
        if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

        c16comb->SaveAs(Form("%s/deltaSystCombined%s%s_%s.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data(),fVariableNames[i].Data()));
        c16comb->SaveAs(Form("%s/deltaSystCombined%s%s_%s.root",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data(),fVariableNames[i].Data()));
        c16comb->SaveAs(Form("%s/deltaSystCombined%s%s_%s.C",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data(),fVariableNames[i].Data()));
        hcombsystdelta[0]->GetYaxis()->SetRangeUser(1e-4,1);
        gPad->SetLogy();
	hcombsystdelta[0]->SetMinimum(.00001);
	hcombsystdelta[0]->SetMaximum(100);
        c16comb->SaveAs(Form("%s/deltaSystCombinedlog%s%s_%s.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data(),fVariableNames[i].Data()));
        c16comb->SaveAs(Form("%s/deltaSystCombinedlog%s%s_%s.root",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data(),fVariableNames[i].Data()));
        c16comb->SaveAs(Form("%s/deltaSystCombinedlog%s%s_%s.C",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data(),fVariableNames[i].Data()));


        TCanvas *c16stack = new TCanvas(Form("c16stack%s",(inorm?"Norm":"")), Form("c16stack%s",(inorm?"Norm":"")),500,500);
        if(inorm) nxchecksys=2;
        else nxchecksys=0;
        TLegend *l16stack[nleg];
        nperleg   = (combsystnames.size()-nxchecksys)/nleg;
        remainder = (combsystnames.size()-nxchecksys)%nleg;

        for (int ileg = 0; ileg < nleg; ++ileg)
        {
          ff[ileg] = double(nperleg + 0.2 + (remainder>ileg?1:0) ) / double(nperleg + 0.2 + (remainder>0?1:0) );
          l16stack[ileg]=new TLegend(gStyle->GetPadLeftMargin()+ileg*padwidth/nleg,1.0 - gStyle->GetPadTopMargin() - ff[ileg]*0.4/nleg, gStyle->GetPadLeftMargin()+(ileg+1)*padwidth/nleg,1.0 - gStyle->GetPadTopMargin());
        }

        THStack *hstack = new THStack(Form("hstack%s",(inorm?"Norm":"")), Form("hstack%s",(inorm?"Norm":"")));
        TH1D* hvartemp[combsystnames.size()];


        //for(unsigned int i_cs = 0;i_cs<combsystnames.size();i_cs++){
        for(int i_cs = combsystnames.size()-1;i_cs>=0;--i_cs){

          if( !inorm || (combsystnames[i_cs] != "Lumi" && combsystnames[i_cs] != "BR") ) { 
            /*
            TH1D* hvartemp[i_cs] = new TH1D(Form("%s_Variance",hcombsystdelta[i_cs]->GetName()), Form("%s_Variance",hcombsystdelta[i_cs]->GetName()), hcombsystdelta[i_cs]->GetNbinsX(),0,hcombsystdelta[i_cs]->GetNbinsX());

            for(int ibin = 1; ibin<=hvartemp[i_cs]->GetNbinsX(); ibin++){
              hvartemp[i_cs]->SetBinContent(ibin, hcombsystdelta[i_cs]->GetBinContent(ibin)*hcombsystdelta[i_cs]->GetBinContent(ibin));
            }
            */

            c16stack->cd();
            hvartemp[i_cs] = (TH1D*) hcombsystdelta[i_cs]->Clone(Form("%s_Variance",hcombsystdelta[i_cs]->GetName()));
            hvartemp[i_cs]->Reset();
            hvartemp[i_cs]->Divide(hcombsystdelta[i_cs],htotsyst_reldelta);
            hvartemp[i_cs]->Multiply(hvartemp[i_cs]);

            for(int ibin = 1; ibin<=hvartemp[i_cs]->GetNbinsX(); ibin++){
              //hvartemp[i_cs]->SetBinContent(ibin, hcombsystdelta[i_cs]->GetBinContent(ibin)*hcombsystdelta[i_cs]->GetBinContent(ibin));
              hvartemp[i_cs]->SetBinError(ibin, 0);
            }

            //hvartemp[i_cs]->SetFillColor(i_cs-10*(i_cs/10));
            hvartemp[i_cs]->SetFillColor(i_cs+1-9*(i_cs/9));
            //hvartemp[i_cs]->SetFillColor(i_cs);
            hvartemp[i_cs]->SetLineWidth(1);
            hvartemp[i_cs]->SetLineColor(1);
            hvartemp[i_cs]->SetLineStyle(1);
            hvartemp[i_cs]->SetMarkerStyle(1);
            hvartemp[i_cs]->GetXaxis()->SetTitleOffset(1.36);
            hvartemp[i_cs]->GetYaxis()->SetTitleOffset(1.6);
            hvartemp[i_cs]->Draw("hist");
            hstack->Add(hvartemp[i_cs]);
            //l16stack[i_cs%nleg]->AddEntry(hvartemp[i_cs],combsystnames[i_cs].Data());
          }

        }

        legcounter=0;
        for(unsigned int i_cs = 0;i_cs<combsystnames.size();i_cs++){
          if( !inorm || (combsystnames[i_cs] != "Lumi" && combsystnames[i_cs] != "BR") ) { 
            l16stack[legcounter%nleg]->AddEntry(hvartemp[i_cs],combsystnames[i_cs].Data());
            ++legcounter;
          }
        }

        c16stack->cd();
        hstack->Draw("");

        if(fVariableNames[i].Contains("_mttbar")) {
          hstack->GetXaxis()->SetNdivisions(-20604);
          c16stack->SetGridx();
        }

        hstack->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
        hstack->GetYaxis()->SetTitle("#sum #Delta#sigma^{2} / #sigma^{2}");
        hstack->GetXaxis()->SetTitleSize(0.03);
        hstack->GetXaxis()->SetLabelSize(0.03);
        hstack->GetYaxis()->SetTitleSize(0.03);
        hstack->GetYaxis()->SetLabelSize(0.03);
        hstack->GetXaxis()->SetTitleOffset(1.36);
        hstack->GetYaxis()->SetTitleOffset(1.6);
        hstack->SetMaximum( hstack->GetMaximum()*1.4);
        for (int ileg = 0; ileg < nleg; ++ileg)
        {
          l16stack[ileg]->SetFillColor(kWhite);
          //l16stack[ileg]->SetFillStyle(0);
          l16stack[ileg]->SetLineColor(kBlack);
          //l16stack[ileg]->SetLineWidth(1);

          l16stack[ileg]->SetTextFont(42);
          l16stack[ileg]->SetTextAlign(12);
          if(inorm) l16stack[ileg]->SetTextSize(0.021);
          else l16stack[ileg]->SetTextSize(0.019);

          l16stack[ileg]->Draw("same");
        }

	Draw2DLabels(fVariableNames[i]);
        DrawCMSLabels(0.03, fLumi);
        if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.03);

        c16stack->SaveAs(Form("%s/deltaSystStack%s%s_%s.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data(),fVariableNames[i].Data()));
        c16stack->SaveAs(Form("%s/deltaSystStack%s%s_%s.root",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data(),fVariableNames[i].Data()));
        c16stack->SaveAs(Form("%s/deltaSystStack%s%s_%s.C",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data(),fVariableNames[i].Data()));


        gStyle->SetPadRightMargin(0.13);

        TCanvas *c17t = new TCanvas("c17t", "c17t");
        htotalsystcovmat[inorm]->SetTitle("");
        htotalsystcovmat[inorm]->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
        htotalsystcovmat[inorm]->GetYaxis()->SetTitle(Form("%s",quantity.Data()));
	htotalsystcovmat[inorm]->GetXaxis()->SetTitleSize(0.035);
	htotalsystcovmat[inorm]->GetXaxis()->SetTitleOffset(1.65);
	htotalsystcovmat[inorm]->GetYaxis()->SetTitleSize(0.035);
	htotalsystcovmat[inorm]->GetYaxis()->SetTitleOffset(1.65);
        htotalsystcovmat[inorm]->GetXaxis()->SetNdivisions(-20604);
        htotalsystcovmat[inorm]->GetYaxis()->SetNdivisions(-20604);
	htotalsystcovmat[inorm]->GetZaxis()->SetLabelSize(0.03);
        if(fVariableNames[i].Contains("_mttbar")) {
          htotalsystcovmat[inorm]->GetXaxis()->SetNdivisions(-20604);
          htotalsystcovmat[inorm]->GetYaxis()->SetNdivisions(-20604);
          htotalsystcovmat[inorm]->GetXaxis()->SetLabelSize(0);
          htotalsystcovmat[inorm]->GetYaxis()->SetLabelSize(0);
          c17t->SetGridx();
          c17t->SetGridy();
        }
        SetAxis(fVariableNames[i], htotalsystcovmat[inorm]->GetXaxis());
        SetAxis(fVariableNames[i], htotalsystcovmat[inorm]->GetYaxis());
        htotalsystcovmat[inorm]->Draw("colz");

	DrawCov2DAxisLabels(fVariableNames[i]);
	Draw2DLabels(fVariableNames[i], 0.85, 0.019);
	Draw2DLabelsVertical(fVariableNames[i], 0.82, 0.019);
        DrawCMSLabels(0.04, fLumi);
        if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

        c17t->SaveAs(Form("%s/TotalSystCovMatrix%s%s_%s.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data(),fVariableNames[i].Data()));
        c17t->SaveAs(Form("%s/TotalSystCovMatrix%s%s_%s.C",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data(),fVariableNames[i].Data()));

        systoutfile->cd();
        htotalsystcovmat[inorm]->Write();

        common::setHHStyle(*gStyle);

        // ************
        // Coefficients
        // ************

        // Amandeep : Hardcoded here
        // TODO: fix these... hardcoded for 1 variable currently
        if (nbinsrhoi % 6 == 0) {
          int numBinsPerCoefficient = 6;
          int numCoefficients = (int)nbinsrhoi / numBinsPerCoefficient;
          TH1D *tunfresult_temp =
              (TH1D *)((TH1D *)unffile->Get(fVariableNames[i] + "TUnfResultCor" + NameTag))
                  ->Clone(fVariableNames[i] + "TUnfResult_temp");
          if (inorm > 0)
            tunfresult_temp->Scale(1. / tunfresult_temp->Integral());
          TH1D *hGenMC_rebinnedS = (TH1D *)((TH1D *)unffile->Get(fVariableNames[i] + "Gen_norm_rebinnedB"))
                                       ->Clone(Form("hGenMC_rebinnedS_%s", fVariableNames[i].Data()));

          Double_t CoefStatOnly[numCoefficients], CoefErrOptStatOnly[numCoefficients],
              CoefErrTestStatOnly[numCoefficients];
          Double_t Coef[numCoefficients], CoefErrOpt[numCoefficients], CoefErrTest[numCoefficients];
          Double_t Afb[numCoefficients], AfbErr[numCoefficients];
          std::vector<double> Afbvec, AfbErrvec, Coefvec, a_optimal_recoStatOnly, a_optimal_reco, a_optimal_reco_temp;
          TMatrixD m_AFB(nbinsrhoi / 2, nbinsrhoi / 2);
          TMatrixD m_C(nbinsrhoi / 2, nbinsrhoi / 2);
          TMatrixD statcovmat = MatrixUnf::TH2toTMatrixD(hstatcovmat);
          bool isCoefficient = kFALSE;
          MatrixUnf::GetAfbWithCorrelationsBinByBin(
              tunfresult_temp, statcovmat, Afbvec, AfbErrvec, m_AFB, numCoefficients);
          isCoefficient = MatrixUnf::GetOptimalCoefficient(Afbvec,
                                                           m_AFB,
                                                           hGenMC_rebinnedS,
                                                           fVariableNames[i],
                                                           CoefStatOnly,
                                                           CoefErrOptStatOnly,
                                                           CoefErrTestStatOnly,
                                                           Coefvec,
                                                           m_C,
                                                           a_optimal_recoStatOnly,
                                                           numCoefficients);
          if (!isCoefficient) {
            MatrixUnf::GetAfbWithCorrelations(tunfresult_temp, statcovmat, Afb, AfbErr, numCoefficients);
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              CoefStatOnly[ithCoefficient] = Afb[ithCoefficient];
              CoefErrTestStatOnly[ithCoefficient] = AfbErr[ithCoefficient];
              CoefErrOptStatOnly[ithCoefficient] = AfbErr[ithCoefficient];
            }
          }
          MatrixUnf::GetAfbWithCorrelationsBinByBin(
              tunfresult_temp, TotalCovMat, Afbvec, AfbErrvec, m_AFB, numCoefficients);
          isCoefficient = MatrixUnf::GetOptimalCoefficient(Afbvec,
                                                           m_AFB,
                                                           hGenMC_rebinnedS,
                                                           fVariableNames[i],
                                                           Coef,
                                                           CoefErrOpt,
                                                           CoefErrTest,
                                                           Coefvec,
                                                           m_C,
                                                           a_optimal_reco,
                                                           numCoefficients);
          if (!isCoefficient) {
            MatrixUnf::GetAfbWithCorrelations(tunfresult_temp, TotalCovMat, Afb, AfbErr, numCoefficients);
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              Coef[ithCoefficient] = Afb[ithCoefficient];
              CoefErrTest[ithCoefficient] = AfbErr[ithCoefficient];
              CoefErrOpt[ithCoefficient] = AfbErr[ithCoefficient];
            }
          }

          if (!isCoefficient) {
            // to reconstruct the inclusive asymmetry from the binned asymmetries, the appropriate factors are the fractions of events in each bin
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              int binOffset = ithCoefficient * numBinsPerCoefficient;
              int asymmetryOffset = ithCoefficient * numBinsPerCoefficient / 2;
              double integral = tunfresult_temp->Integral(1 + binOffset, binOffset + numBinsPerCoefficient);
              for (int i_bin = 0; i_bin < numBinsPerCoefficient / 2; ++i_bin) {
                double a_temp = (tunfresult_temp->GetBinContent(binOffset + i_bin + 1) +
                                 tunfresult_temp->GetBinContent(binOffset + numBinsPerCoefficient - i_bin)) /
                                integral;
                a_optimal_recoStatOnly.at(asymmetryOffset + i_bin) = a_temp;
                a_optimal_reco.at(asymmetryOffset + i_bin) = a_temp;
              }
            }
          }

          if (doOptimiseTotalUnc) {
            a_optimal_reco_temp = a_optimal_reco;
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              coefTable[i_tag * 2 + inorm][i][ithCoefficient] = Coef[ithCoefficient];
            }
          } else {
            a_optimal_reco_temp = a_optimal_recoStatOnly;
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              coefTable[i_tag * 2 + inorm][i][ithCoefficient] = CoefStatOnly[ithCoefficient];
            }
          }

          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            int binOffset = ithCoefficient * numBinsPerCoefficient;
            binsumTable[i_tag * 2 + inorm][i][ithCoefficient] =
                tunfresult_temp->Integral(binOffset + 1, numBinsPerCoefficient + binOffset);
          }
          double binsum_temp[numCoefficients], binsumerrsysttotal[numCoefficients], binsumerrstat[numCoefficients],
              binsumerrMCstat[numCoefficients], binsumerrBkgstatTotal[numCoefficients], binsumerrtotal[numCoefficients];

          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            int asymmetryOffset = ithCoefficient * numBinsPerCoefficient / 2;
            std::vector<double> factor{1., 1., 1.};
            MatrixUnf::AtoCfactor(
                fVariableNames[i], hGenMC_rebinnedS, factor, numCoefficients, ithCoefficient);

            hist_acombfactor_allVar[i_tag * 2 + inorm][asymmetryOffset + 0]->SetBinContent(
                i + 1, (a_optimal_reco_temp.at(asymmetryOffset + 0)) / factor[0]);
            hist_acombfactor_allVar[i_tag * 2 + inorm][asymmetryOffset + 1]->SetBinContent(
                i + 1, (a_optimal_reco_temp.at(asymmetryOffset + 1)) / factor[1]);
            hist_acombfactor_allVar[i_tag * 2 + inorm][asymmetryOffset + 2]->SetBinContent(
                i + 1,
                (1. - a_optimal_reco_temp.at(asymmetryOffset + 0) - a_optimal_reco_temp.at(asymmetryOffset + 1)) /
                    factor[2]);
          }

          for (int syst = 0; syst < nsyst; syst++) {
            double Cerrsyst[numCoefficients];
            calculateCoefficientUncertainty(hcovtemp[syst],
                                            tunfresult_temp,
                                            hGenMC_rebinnedS,
                                            a_optimal_reco_temp,
                                            fVariableNames[i],
                                            binsum_temp,
                                            numCoefficients,
                                            Cerrsyst);
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              std::cout << "Cerr at " << ithCoefficient << " is " << Cerrsyst[ithCoefficient] << std::endl;
            }
            printf("%s Coef%s%s %s Uncertainty: ",
                   fVariableNames[i].Data(),
                   NameTag.Data(),
                   (inorm ? "Norm" : ""),
                   tempVectorOfInputSystematics[syst].c_str());
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              systTable[i_tag * 2 + inorm][i][ithCoefficient][syst] = Cerrsyst[ithCoefficient];
              binsumsystTable[i_tag * 2 + inorm][i][ithCoefficient][syst] = binsum_temp[ithCoefficient];
            }

            delete hcovtemp[syst];
          }

          for (unsigned int i_cs = 0; i_cs < combsystnames.size(); i_cs++) {
            double Cerrsyst[numCoefficients];
            calculateCoefficientUncertainty(hcombsystcovmat[i_cs],
                                            tunfresult_temp,
                                            hGenMC_rebinnedS,
                                            a_optimal_reco_temp,
                                            fVariableNames[i],
                                            binsum_temp,
                                            numCoefficients,
                                            Cerrsyst);
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              printf("%s Coef%s%s %s Combined Uncertainty: %f \n",
                     fVariableNames[i].Data(),
                     NameTag.Data(),
                     (inorm ? "Norm" : ""),
                     combsystnames[i_cs].Data(),
                     Cerrsyst[ithCoefficient]);
              combsystTable[i_tag * 2 + inorm][i][ithCoefficient][i_cs] = Cerrsyst[ithCoefficient];
              binsumcombsystTable[i_tag * 2 + inorm][i][ithCoefficient][i_cs] = binsum_temp[ithCoefficient];
            }

            delete hcombsystcovmat[i_cs];
          }
          double Cerrsysttotal[numCoefficients], Cerrstat[numCoefficients], CerrMCstat[numCoefficients];

          calculateCoefficientUncertainty(htotalsystcovmat[inorm],
                                          tunfresult_temp,
                                          hGenMC_rebinnedS,
                                          a_optimal_reco_temp,
                                          fVariableNames[i],
                                          binsumerrsysttotal,
                                          numCoefficients,
                                          Cerrsysttotal);
          calculateCoefficientUncertainty(hstatcovmat,
                                          tunfresult_temp,
                                          hGenMC_rebinnedS,
                                          a_optimal_reco_temp,
                                          fVariableNames[i],
                                          binsumerrstat,
                                          numCoefficients,
                                          Cerrstat);
          calculateCoefficientUncertainty(hMCstatcovmat,
                                          tunfresult_temp,
                                          hGenMC_rebinnedS,
                                          a_optimal_reco_temp,
                                          fVariableNames[i],
                                          binsumerrMCstat,
                                          numCoefficients,
                                          CerrMCstat);
          // printf("%s Total Coef%s%s Uncertainties: %f %f %f \n", fVariableNames[i].Data(),NameTag.Data(),(inorm?"Norm":""), Cerrstat, CerrMCst at, Cerrsysttotal );
          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            systTable[i_tag * 2 + inorm][i][ithCoefficient][nsyst] = Cerrsysttotal[ithCoefficient];
            systTable[i_tag * 2 + inorm][i][ithCoefficient][nsyst + 1] = Cerrstat[ithCoefficient];
            systTable[i_tag * 2 + inorm][i][ithCoefficient][nsyst + 2] = CerrMCstat[ithCoefficient];
            systTable[i_tag * 2 + inorm][i][ithCoefficient][nsyst + 4] =
                sqrt(Cerrstat[ithCoefficient] * Cerrstat[ithCoefficient] +
                     CerrMCstat[ithCoefficient] * CerrMCstat[ithCoefficient]);
            systTable[i_tag * 2 + inorm][i][ithCoefficient][nsyst + 5] =
                sqrt(Cerrstat[ithCoefficient] * Cerrstat[ithCoefficient] +
                     CerrMCstat[ithCoefficient] * CerrMCstat[ithCoefficient] +
                     Cerrsysttotal[ithCoefficient] * Cerrsysttotal[ithCoefficient]);

            binsumsystTable[i_tag * 2 + inorm][i][ithCoefficient][nsyst] =
                binsumerrsysttotal[ithCoefficient];
            binsumsystTable[i_tag * 2 + inorm][i][ithCoefficient][nsyst + 1] = binsumerrstat[ithCoefficient];
            binsumsystTable[i_tag * 2 + inorm][i][ithCoefficient][nsyst + 2] =
                binsumerrMCstat[ithCoefficient];
            binsumsystTable[i_tag * 2 + inorm][i][ithCoefficient][nsyst + 4] =
                sqrt(binsumerrstat[ithCoefficient] * binsumerrstat[ithCoefficient] +
                     binsumerrMCstat[ithCoefficient] * binsumerrMCstat[ithCoefficient]);
            binsumsystTable[i_tag * 2 + inorm][i][ithCoefficient][nsyst + 5] =
                sqrt(binsumerrstat[ithCoefficient] * binsumerrstat[ithCoefficient] +
                     binsumerrMCstat[ithCoefficient] * binsumerrMCstat[ithCoefficient] +
                     binsumerrsysttotal[ithCoefficient] * binsumerrsysttotal[ithCoefficient]);
          }

          // note, the BG stat uncertainty is already included in hstatcovmat ( from GetEmatrix() ) - just printed here for info.
          double CerrBkgstat[backgroundnames.size()][numCoefficients];
          double CerrBkgstatTotal[numCoefficients];
          if (fDoSubtrBg) {
            for (unsigned int ibkg = 0; ibkg < backgroundnames.size(); ++ibkg) {
              double dummyCoefUnc[numCoefficients];
              calculateCoefficientUncertainty(hBkgstatcovmat[ibkg],
                                              tunfresult_temp,
                                              hGenMC_rebinnedS,
                                              a_optimal_reco_temp,
                                              fVariableNames[i],
                                              binsum_temp,
                                              numCoefficients,
                                              dummyCoefUnc);

              for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
                CerrBkgstat[ibkg][ithCoefficient] = dummyCoefUnc[ithCoefficient];
              }
            }
            calculateCoefficientUncertainty(hBkgstatcovmatTotal,
                                            tunfresult_temp,
                                            hGenMC_rebinnedS,
                                            a_optimal_reco_temp,
                                            fVariableNames[i],
                                            binsumerrBkgstatTotal,
                                            numCoefficients,
                                            CerrBkgstatTotal);
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              systTable[i_tag * 2 + inorm][i][ithCoefficient][nsyst + 3] = CerrBkgstatTotal[ithCoefficient];
              systTable[i_tag * 2 + inorm][i][ithCoefficient][nsyst + 1] =
                  sqrt(Cerrstat[ithCoefficient] * Cerrstat[ithCoefficient] -
                       CerrBkgstatTotal[ithCoefficient] * CerrBkgstatTotal[ithCoefficient]);

              binsumsystTable[i_tag * 2 + inorm][i][ithCoefficient][nsyst + 3] =
                  binsumerrBkgstatTotal[ithCoefficient];
              binsumsystTable[i_tag * 2 + inorm][i][ithCoefficient][nsyst + 1] =
                  sqrt(binsumerrstat[ithCoefficient] * binsumerrstat[ithCoefficient] -
                       binsumerrBkgstatTotal[ithCoefficient] * binsumerrBkgstatTotal[ithCoefficient]);
            }
          }

          if (fDoSubtrBg) {
            printf("%s Total Coef%s%s BkgStatUncertainties: ",
                   fVariableNames[i].Data(),
                   NameTag.Data(),
                   (inorm ? "Norm" : ""));
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              printf("For the %dth coefficient:", ithCoefficient);
              for (unsigned int ibkg = 0; ibkg < backgroundnames.size(); ++ibkg) {
                printf(" %s: %f ", backgroundnames.at(ibkg).Data(), CerrBkgstat[ibkg][ithCoefficient]);
              }
              printf(" Total: %f \n", CerrBkgstatTotal[ithCoefficient]);
            }
          }
          double Cerrtotal[numCoefficients];
          calculateCoefficientUncertainty(
              MatrixUnf::TMatrixDtoTH2(TotalCovMat, Form("hTotalCovMat%s%s", (inorm ? "Norm" : ""), NameTag.Data())),
              tunfresult_temp,
              hGenMC_rebinnedS,
              a_optimal_reco_temp,
              fVariableNames[i],
              binsumerrtotal,
              numCoefficients,
              Cerrtotal);

          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            printf("%s Total %dth Coef%s%s Uncertainties: %f %f %f %f %g %f %f %f %f %f %f %f %f %f %f %f %f \n",
                   fVariableNames[i].Data(),
                   ithCoefficient,
                   NameTag.Data(),
                   (inorm ? "Norm" : ""),
                   Cerrstat[ithCoefficient],
                   CerrMCstat[ithCoefficient],
                   Cerrsysttotal[ithCoefficient],
                   Cerrtotal[ithCoefficient],
                   Cerrtotal[ithCoefficient] * Cerrtotal[ithCoefficient] -
                       Cerrsysttotal[ithCoefficient] * Cerrsysttotal[ithCoefficient] -
                       CerrMCstat[ithCoefficient] * CerrMCstat[ithCoefficient] -
                       Cerrstat[ithCoefficient] * Cerrstat[ithCoefficient],
                   a_optimal_reco[ithCoefficient * 3 + 0],
                   a_optimal_reco[ithCoefficient * 3 + 1],
                   a_optimal_reco[ithCoefficient * 3 + 2],
                   Coef[ithCoefficient],
                   CoefErrOpt[ithCoefficient],
                   CoefErrTest[ithCoefficient],
                   a_optimal_recoStatOnly[ithCoefficient * 3 + 0],
                   a_optimal_recoStatOnly[ithCoefficient * 3 + 1],
                   a_optimal_recoStatOnly[ithCoefficient * 3 + 2],
                   CoefStatOnly[ithCoefficient],
                   CoefErrOptStatOnly[ithCoefficient],
                   CoefErrTestStatOnly[ithCoefficient]);
          }

          delete tunfresult_temp;
          delete hGenMC_rebinnedS;
        }
      } //inorm
    } //i_tag

    // Now set up the theoretical predictions for fitting
    common::setHHStyle(*gStyle);

    TH1D* theoryxsectemp = (TH1D*)unffile->Get(fVariableNames[i] + "TheoryXsec");
    
    // Amandeep : Theory xsection hardcoded here ?
    std::cout << "Correcting theory xsec: " << 830.91/theoryxsectemp->Integral() << std::endl;
    theoryxsectemp->Scale(830.91/theoryxsectemp->Integral());
    // End

    /*
    theoryxsec->SetLineColor(kRed);
    theoryxsec->SetLineWidth(2);
    theoryxsec->SetMarkerColor(kRed);
    theoryxsec->SetMarkerStyle(24);
    theoryxsec->SetMarkerSize(0.6);
    theoryxsec->SetTitle("");
    */

    // Amandeep : For debugging
    if     (MatrixUnf::isCtype(fVariableNames[i]))   std::cout << fVariableNames[i] << " is a Ctype() variable " << std::endl;
    else if(MatrixUnf::isCPMtype(fVariableNames[i])) std::cout << fVariableNames[i] << " is a CPMtype() variable " << std::endl;
    else if(MatrixUnf::isBtype(fVariableNames[i]))   std::cout << fVariableNames[i] << " is a Btype() variable " << std::endl;
    else if(MatrixUnf::isDtype(fVariableNames[i]))   std::cout << fVariableNames[i] << " is a Dtype() variable " << std::endl;
    else if(MatrixUnf::isOtherSymType(fVariableNames[i]))   std::cout << fVariableNames[i] << " is a isOtherSymType() variable " << std::endl;   
    else if(MatrixUnf::isLinearCombtype(fVariableNames[i])) std::cout << fVariableNames[i] << " is a isLinearCombtype() variable " << std::endl;
    // End

    // Amandeep : This is nbins_final 24(2D) or 6(1D)
    const int nbinsrhoi = theoryxsectemp->GetNbinsX();

    // Amandeep : Should this be physical binwidth ?
    for (int i_bin = 1; i_bin <= nbinsrhoi; ++i_bin)
    {
      hbinwidth_AllVar->SetBinContent(i_bin+i*nbinsrhoi, theoryxsectemp->GetBinWidth(i_bin) );
      hbinwidth_AllVar->SetBinError(i_bin+i*nbinsrhoi, 0. );

      hbinlowedge_AllVar->SetBinContent(i_bin+i*nbinsrhoi, theoryxsectemp->GetXaxis()->GetBinLowEdge(i_bin) );
      hbinlowedge_AllVar->SetBinError(i_bin+i*nbinsrhoi, 0. );
    }

    // Analytic functions
    TF1 *analyticdiffxsec;
    TF1 *analyticdiffxsec2;
    TF1 *analyticdiffxsec3;
    TF1 *analyticdiffxsec4;

    if     ( MatrixUnf::isCtype(fVariableNames[i]) )  analyticdiffxsec  = new TF1("analyticdiffxsec","(1-[0]*x)*log(abs(x))/(-2)", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) );
    else if( MatrixUnf::isCPMtype(fVariableNames[i])) analyticdiffxsec  = new TF1("analyticdiffxsec","(1-[0]*x/2)*acos(abs(x))/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 
    else if( MatrixUnf::isBtype(fVariableNames[i]) )  analyticdiffxsec  = new TF1("analyticdiffxsec","(1+[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 
    else if( MatrixUnf::isDtype(fVariableNames[i]) )  analyticdiffxsec  = new TF1("analyticdiffxsec","(1-[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 

    // Amandeep : For current variables that break 
    // Jacob suggested inserting a dummy BType functional form    
    else if( MatrixUnf::isOtherSymType(fVariableNames[i]) )  analyticdiffxsec  = new TF1("analyticdiffxsec","(1+[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) );     
    
    // Amandeep : For the new variables
    // There is an assumption here that since they are flat they should be Dtype
    else if( MatrixUnf::isLinearCombtype(fVariableNames[i]) )  analyticdiffxsec = new TF1("analyticdiffxsec","(1-[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 
    // End

    if     ( MatrixUnf::isCtype(fVariableNames[i]) )   analyticdiffxsec2 = new TF1("analyticdiffxsec2","(1-[0]*x)*log(abs(x))/(-2)", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) );
    else if( MatrixUnf::isCPMtype(fVariableNames[i]))  analyticdiffxsec2 = new TF1("analyticdiffxsec2","(1-[0]*x/2)*acos(abs(x))/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 
    else if( MatrixUnf::isBtype(fVariableNames[i]) )   analyticdiffxsec2 = new TF1("analyticdiffxsec2","(1+[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 
    else if( MatrixUnf::isDtype(fVariableNames[i]) )   analyticdiffxsec2 = new TF1("analyticdiffxsec2","(1-[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 

    // Amandeep : For current variables that break  
    // Jacob suggested inserting a dummy BType functional form 
    else if( MatrixUnf::isOtherSymType(fVariableNames[i]) )  analyticdiffxsec2  = new TF1("analyticdiffxsec2","(1+[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) );     
    
    // Amandeep : For the new variables
    // There is an assumption here that since they are flat they should be Dtype
    else if( MatrixUnf::isLinearCombtype(fVariableNames[i]) )  analyticdiffxsec2 = new TF1("analyticdiffxsec2","(1-[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 
    // End

    if     ( MatrixUnf::isCtype(fVariableNames[i]) )    analyticdiffxsec3 = new TF1("analyticdiffxsec3","(1-[0]*x)*log(abs(x))/(-2)", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) );
    else if( MatrixUnf::isCPMtype(fVariableNames[i]))   analyticdiffxsec3 = new TF1("analyticdiffxsec3","(1-[0]*x/2)*acos(abs(x))/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 
    else if( MatrixUnf::isBtype(fVariableNames[i]) )    analyticdiffxsec3 = new TF1("analyticdiffxsec3","(1+[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 
    else if( MatrixUnf::isDtype(fVariableNames[i]) )    analyticdiffxsec3 = new TF1("analyticdiffxsec3","(1-[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 

    // Amandeep : For current variables that break    
    // Jacob suggested inserting a dummy BType functional form 
    else if( MatrixUnf::isOtherSymType(fVariableNames[i]) )  analyticdiffxsec3  = new TF1("analyticdiffxsec3","(1+[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) );     
   
    // Amandeep : For the new variables
    // There is an assumption here that since they are flat they should be Dtype
    else if( MatrixUnf::isLinearCombtype(fVariableNames[i]) )  analyticdiffxsec3 = new TF1("analyticdiffxsec3","(1-[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 
    // End

    if     ( MatrixUnf::isCtype(fVariableNames[i]) )   analyticdiffxsec4 = new TF1("analyticdiffxsec4","(1-[0]*x)*log(abs(x))/(-2)", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) );
    else if( MatrixUnf::isCPMtype(fVariableNames[i]))  analyticdiffxsec4 = new TF1("analyticdiffxsec4","(1-[0]*x/2)*acos(abs(x))/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 
    else if( MatrixUnf::isBtype(fVariableNames[i]) )   analyticdiffxsec4 = new TF1("analyticdiffxsec4","(1+[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 
    else if( MatrixUnf::isDtype(fVariableNames[i]) )   analyticdiffxsec4 = new TF1("analyticdiffxsec4","(1-[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 

    // Amandeep : For current variables that break    
    // Jacob suggested inserting a dummy BType functional form 
    else if( MatrixUnf::isOtherSymType(fVariableNames[i]) )  analyticdiffxsec4  = new TF1("analyticdiffxsec4","(1+[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) );     
   
    // Amandeep : For the new variables
    // There is an assumption here that since they are flat they should be Dtype
    else if( MatrixUnf::isLinearCombtype(fVariableNames[i]) )  analyticdiffxsec4 = new TF1("analyticdiffxsec4","(1-[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 
    // End


    TF1 *analyticdiffxsecMulti[nFit];

    for (int iM = 0; iM < nFit; ++iM)
    {
      if     ( MatrixUnf::isCtype(fVariableNames[i]) )   analyticdiffxsecMulti[iM] = new TF1(Form("analyticdiffxsecMulti%i",iM),"(1-[0]*x)*log(abs(x))/(-2)", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) );
      else if( MatrixUnf::isCPMtype(fVariableNames[i]) ) analyticdiffxsecMulti[iM] = new TF1(Form("analyticdiffxsecMulti%i",iM),"(1-[0]*x/2)*acos(abs(x))/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 
      else if( MatrixUnf::isBtype(fVariableNames[i]) )   analyticdiffxsecMulti[iM] = new TF1(Form("analyticdiffxsecMulti%i",iM),"(1+[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 
      else if( MatrixUnf::isDtype(fVariableNames[i]) )   analyticdiffxsecMulti[iM] = new TF1(Form("analyticdiffxsecMulti%i",iM),"(1-[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 

      // Amandeep : For current variables that break    
      // Jacob suggested inserting a dummy BType functional form 
      else if( MatrixUnf::isOtherSymType(fVariableNames[i]) )  analyticdiffxsecMulti[iM]  = new TF1(Form("analyticdiffxsecMulti%i",iM),"(1+[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) );
    
      // Amandeep : For the new variables
      // There is an assumption here that since they are flat they should be Dtype
      else if( MatrixUnf::isLinearCombtype(fVariableNames[i]) ) analyticdiffxsecMulti[iM] = new TF1(Form("analyticdiffxsecMulti%i",iM),"(1-[0]*x)/2", theoryxsectemp->GetXaxis()->GetBinLowEdge(1), theoryxsectemp->GetXaxis()->GetBinUpEdge(theoryxsectemp->GetNbinsX()) ); 
      // End
    }

    std::cout << "Using coef = " << NLOWvalues[i] << std::endl;

    if( MatrixUnf::isCtype(fVariableNames[i]) || fVariableNames[i] == "c_Prk" || MatrixUnf::isDtype(fVariableNames[i]) || MatrixUnf::isOtherSymType(fVariableNames[i]) ) {
      fSMTable.push_back((coefTable[1][i][0] - NLOWvaluesUnCorr[i]) / (NLOWvalues[i] - NLOWvaluesUnCorr[i]));
      fSMTableErr.push_back(systTable[1][i][0][nsyst + 5] / fabs(NLOWvalues[i] - NLOWvaluesUnCorr[i]));      // total excluding theory uncertainties
      fSMTableStatErr.push_back(systTable[1][i][0][nsyst + 4] / fabs(NLOWvalues[i] - NLOWvaluesUnCorr[i]));  // stat
      fSMTableSystErr.push_back(systTable[1][i][0][nsyst] / fabs(NLOWvalues[i] - NLOWvaluesUnCorr[i]));      // syst
      fSMTableNames.push_back( MatrixUnf::CoefName(fVariableNames[i],kTRUE)  );
    }

    // Amandeep : adding here for checks
    std::cout << "fSM Table" << std::endl;
    for (double i : fSMTable)        { std::cout << "Coefficients :: "<< i << std::endl; }
    for (double i : fSMTableErr)     { std::cout << "Errors :: "      << i << std::endl; }
    for (double i : fSMTableStatErr) { std::cout << "Stat Errors :: " << i << std::endl; }
    for (double i : fSMTableSystErr) { std::cout << "Syst Errors :: " << i << std::endl; }
    for (auto i   : fSMTableNames)   { std::cout << "Element Name ::" << i << std::endl; }
    for (int iM = 0; iM < nFit; ++iM) {analyticdiffxsecMulti[iM]->Print();}
    // End

    const int Npx = 10000;
    analyticdiffxsec->FixParameter(0,NLOWvalues[i]);
    analyticdiffxsec->SetNpx(Npx*nbinsrhoi);
    TH1D* tempNLOWxsec = (TH1D*) analyticdiffxsec->GetHistogram();
    tempNLOWxsec->Rebin(Npx);
    //tempNLOWxsec->Scale(1./Npx/3.);
    tempNLOWxsec->Scale(1./tempNLOWxsec->Integral());

    analyticdiffxsec2->FixParameter(0,0);
    analyticdiffxsec2->SetNpx(Npx*nbinsrhoi);
    TH1D* tempNLOWxsecuncorr = (TH1D*) analyticdiffxsec2->GetHistogram();
    tempNLOWxsecuncorr->Rebin(Npx);
    //tempNLOWxsecuncorr->Scale(1./Npx/3.);
    tempNLOWxsecuncorr->Scale(1./tempNLOWxsecuncorr->Integral());

    // TODO: hard-coded to 0th coefficient
    analyticdiffxsec3->FixParameter(0, coefTable[0][i][0]);
    analyticdiffxsec3->SetNpx(Npx*nbinsrhoi);
    TH1D* tempNLOWxsec3 = (TH1D*) analyticdiffxsec3->GetHistogram();
    tempNLOWxsec3->Rebin(Npx);
    //tempNLOWxsec3->Scale(1./Npx/3.);
    tempNLOWxsec3->Scale(1./tempNLOWxsec3->Integral());

    // TODO: hard-coded to 0th coefficient
    analyticdiffxsec4->FixParameter(0, coefTable[1][i][0]);
    analyticdiffxsec4->SetNpx(Npx*nbinsrhoi);
    TH1D* tempNLOWxsec4 = (TH1D*) analyticdiffxsec4->GetHistogram();
    tempNLOWxsec4->Rebin(Npx);
    //tempNLOWxsec4->Scale(1./Npx/3.);
    tempNLOWxsec4->Scale(1./tempNLOWxsec4->Integral());


    // Amandeep : This is what causes the ll_clab and llbar_delta_phi to fail
    TH1D* tempNLOWxsecMulti[nFit];
    for (int iM = 0; iM < nFit; ++iM)
    {
      analyticdiffxsecMulti[iM]->FixParameter(0,NLOWvalues[i]+0.1*mutTableMulti[iM][i]);
      analyticdiffxsecMulti[iM]->SetNpx(Npx*nbinsrhoi);
      tempNLOWxsecMulti[iM] = ( (TH1D*) ( (analyticdiffxsecMulti[iM])->GetHistogram() ) );
      tempNLOWxsecMulti[iM]->Rebin(Npx);
      //tempNLOWxsecMulti[iM]->Scale(1./Npx/3.);
      tempNLOWxsecMulti[iM]->Scale(1./tempNLOWxsecMulti[iM]->Integral());
    }


    TH1D* NLOW_histo[12];
    std::vector<TString> TheoryDataFileList;

    if (fVariableNames[i].Contains("llbar_delta_phi") || fVariableNames[i].Contains("llbar_delta_eta"))
    {
      TheoryDataFileList.push_back("theory/lhc13_mu1m_dphill_spin.dat");
      TheoryDataFileList.push_back("theory/lhc13_mu2m_dphill_spin.dat");
      TheoryDataFileList.push_back("theory/lhc13_muhm_dphill_spin.dat");
      TheoryDataFileList.push_back("theory/lhc13_mu1m_dphill_nospin.dat");
      TheoryDataFileList.push_back("theory/lhc13_mu2m_dphill_nospin.dat");
      TheoryDataFileList.push_back("theory/lhc13_muhm_dphill_nospin.dat");

      for (int i_th = 0; i_th < 12; ++i_th) {
        NLOW_histo[i_th] = (TH1D*) theoryxsectemp->Clone(fVariableNames[i]+(i_th>=6?"LOxsec":"NLOWxsec")+TheoryDataFileList.at(i_th%6).Data());
        NLOW_histo[i_th]->Reset();
      }

      for (int i_th = 0; i_th < 6; ++i_th) {
        Float_t v0, v1, v2, v3, v4;
        Int_t ncols, nlines;
        nlines = 0;
        FILE *fp = fopen(TheoryDataFileList.at(i_th).Data(), "r");
        while (1)
        {
          //    \phi*6/pi    \phi_min               \phi_max              1/sig dsig/d\phi: LO QCD     NLO(expanded)
          ncols = fscanf(fp, "%f %f %f %f %f", &v0, &v1, &v2, &v3, &v4);
          if (ncols < 0) break;
          if(ncols==5) {
            printf("%s=%8f, v=%8f\n",fVariableNames[i].Data(), (v1+v2)/2., v4);
            NLOW_histo[i_th]->SetBinContent(nlines+1,v4);
            NLOW_histo[i_th+6]->SetBinContent(nlines+1,v3);
            nlines++;
          }
        }
        fclose (fp);
      }
    }


    if (fVariableNames[i].Contains("ll_cLab"))
    {

      TheoryDataFileList.push_back("theory/lhc13_mu1m_open_spin.dat");
      TheoryDataFileList.push_back("theory/lhc13_mu2m_open_spin.dat");
      TheoryDataFileList.push_back("theory/lhc13_muhm_open_spin.dat");
      TheoryDataFileList.push_back("theory/lhc13_mu1m_open_nospin.dat");
      TheoryDataFileList.push_back("theory/lhc13_mu2m_open_nospin.dat");
      TheoryDataFileList.push_back("theory/lhc13_muhm_open_nospin.dat");


      for (int i_th = 0; i_th < 12; ++i_th)
      {

        NLOW_histo[i_th] = (TH1D*) theoryxsectemp->Clone(fVariableNames[i]+(i_th>=6?"LOxsec":"NLOWxsec")+TheoryDataFileList.at(i_th%6).Data());
        NLOW_histo[i_th]->Reset();

      }

      for (int i_th = 0; i_th < 6; ++i_th)
      {

        Float_t v1, v2, v3, v4;
        double v3_rebinned=0;
        double v4_rebinned=0;
        Int_t ncols, nlines;
        nlines = 0;
        FILE *fp = fopen(TheoryDataFileList.at(i_th).Data(), "r");
        while (1)
        {
          //   cos\phi_ll(min)          cos\phi_ll(max)               LO QCD                 NLO(expanded)
          ncols = fscanf(fp, "%f %f %f %f", &v1, &v2, &v3, &v4);
          if (ncols < 0) break;

          if(ncols==4) {
            printf("%s=%8f, v=%8f\n",fVariableNames[i].Data(), (v1+v2)/2., v4);

            if(nlines%2==0) {
              v4_rebinned=v4;
              v3_rebinned=v3;
            }
            else{
              v4_rebinned+=v4;
              v3_rebinned+=v3;
              NLOW_histo[i_th]->SetBinContent(nlines/2+1,v4_rebinned/2.);
              NLOW_histo[i_th+6]->SetBinContent(nlines/2+1,v3_rebinned/2.);
            }
            nlines++;

          }

        }
        fclose (fp);

      //        NLOW_histo[i_th]->Rebin(2);
      //        NLOW_histo[i_th+6]->Rebin(2);
      //
      //        NLOW_histo[i_th]->Scale(0.5);
      //        NLOW_histo[i_th+6]->Scale(0.5);

      }
    }


    TH1D* NNLO_histo[7];
    std::vector<TString> NNLOFileList;

    if (fVariableNames[i].Contains("llbar_delta_phi") || fVariableNames[i].Contains("llbar_delta_eta"))
    {

      NNLOFileList.push_back("theory/NNLO-HT4-NNPDF31-dPhiLl_2.dat");
      for (int i_th = 0; i_th < 7; ++i_th)
      {
        NNLO_histo[i_th] = (TH1D*) theoryxsectemp->Clone( Form("%sNNLOxsec%i",fVariableNames[i].Data(),i_th) );
        NNLO_histo[i_th]->Reset();
      }

        Float_t v[16] = {0};
        Int_t ncols, nlines;
        nlines = 0;
        FILE *fp = fopen(NNLOFileList.at(0).Data(), "r");
        while (1)
        {
          //    \phi*6/pi    \phi_min               \phi_max              1/sig dsig/d\phi: LO QCD     NLO(expanded)
          ncols = fscanf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", &v[0], &v[1], &v[2], &v[3], &v[4], &v[5], &v[6], &v[7], &v[8], &v[9], &v[10], &v[11], &v[12], &v[13], &v[14], &v[15]);
          if (ncols < 0) break;

          if(ncols==16) {
            printf("%s=%8f, v=%8f\n",fVariableNames[i].Data(), v[0], v[1]);

            for (int i_th = 0; i_th < 7; ++i_th) {
              NNLO_histo[i_th]->SetBinContent(nlines+1,v[i_th+1]);
              NNLO_histo[i_th]->SetBinError(nlines+1,v[i_th+8]);
            }

            nlines++;
          }
        }
        fclose (fp);

        for (int i_th = 0; i_th < 7; ++i_th)
        {
          NNLO_histo[i_th]->Scale(1./NNLO_histo[i_th]->Integral());
        }
    }

    // TODO: hard-coded to 6 bins per coefficient that is extracted
    const int numBinsPerCoefficient = 6;
    const int numCoefficients       = nbinsrhoi / numBinsPerCoefficient;

    //TH1D *hAfbtheoryxsec = (TH1D *)theoryxsectemp->Clone(fVariableNames[i] + "Afbtheoryxsec");
    TH1D *hAfbtheoryxsec = new TH1D(fVariableNames[i] + "Afbtheoryxsec", "", (numCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
    MatrixUnf::fillAsymsAndCoefficient(theoryxsectemp, nullptr, theoryxsectemp, fVariableNames[i], hAfbtheoryxsec, nullptr, true, true, numCoefficients);

    double nominal_coef[numCoefficients];
    double nominal_coef_staterr[numCoefficients];
    for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
      std::cout << "Nominal(PP8) coef = " << hAfbtheoryxsec->GetBinContent(ithCoefficient * 8 + numBinsPerCoefficient - 1) << " +/- " << hAfbtheoryxsec->GetBinError(ithCoefficient * 8 + numBinsPerCoefficient - 1) << std::endl;
      nominal_coef[ithCoefficient] = hAfbtheoryxsec->GetBinContent(ithCoefficient * 8 + numBinsPerCoefficient - 1);
      nominal_coef_staterr[ithCoefficient] = hAfbtheoryxsec->GetBinError(ithCoefficient * 8 + numBinsPerCoefficient - 1);
    }

    for (int i_bin = 1; i_bin <= nbinsrhoi; ++i_bin)
    {
      hpowheg_AllVar[0]->SetBinContent(i_bin+i*nbinsrhoi, theoryxsectemp->GetBinContent(i_bin) );
      hpowheg_AllVar[0]->SetBinError(i_bin+i*nbinsrhoi, theoryxsectemp->GetBinError(i_bin) );

      hpowheg_AllVar[1]->SetBinContent(i_bin+i*nbinsrhoi, theoryxsectemp->GetBinContent(i_bin)/theoryxsectemp->Integral() );
      hpowheg_AllVar[1]->SetBinError(i_bin+i*nbinsrhoi, theoryxsectemp->GetBinError(i_bin)/theoryxsectemp->Integral() );
    }


    // Amandeep : The plotting fails since it can't retrieve RespMatSys_ for these variables
    std::cout << "Onto AMCATNLO histograms" << std::endl;
    /*TH1D *AMCATNLOxsectemptemp = (TH1D *)((TH2D *)unffile->Get(fVariableNames[i] + "RespMatSys_AMCATNLOFXFX"))->ProjectionX();
    for (int i_bin = 0; i_bin < AMCATNLOxsectemptemp->GetNbinsX(); ++i_bin) {
      if (!(AMCATNLOxsectemptemp->GetBinError(i_bin + 1) >= 0)) AMCATNLOxsectemptemp->SetBinError(i_bin + 1, sqrt(AMCATNLOxsectemptemp->GetBinContent(i_bin + 1)));
    }
    MatrixUnf::rebinMultidimensionalInput(AMCATNLOxsectemptemp, binning, generatorBinning_fine->FindNode("ttbargen"), rebinfine);
    TH1D *AMCATNLOxsectemp = new TH1D(AMCATNLOxsectemptemp->GetName(), AMCATNLOxsectemptemp->GetTitle(), nbinsrhoi, 0.5, nbinsrhoi * rebinfine + 0.5);
    for (int i_bin = 0; i_bin < AMCATNLOxsectemptemp->GetNbinsX(); ++i_bin) {
      AMCATNLOxsectemp->SetBinContent(i_bin + 1, AMCATNLOxsectemptemp->GetBinContent(i_bin + 1));
      AMCATNLOxsectemp->SetBinError(i_bin + 1, AMCATNLOxsectemptemp->GetBinError(i_bin + 1));
    }
    TH1D *hAfbAMCATNLOxsec = new TH1D(fVariableNames[i] + "AfbAMCATNLOxsec", "", (numCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
    MatrixUnf::fillAsymsAndCoefficient(AMCATNLOxsectemp, nullptr, theoryxsectemp, fVariableNames[i], hAfbAMCATNLOxsec, nullptr, true, true, numCoefficients);*/

    for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
      /*std::cout << "AMCATNLO coef = " << hAfbAMCATNLOxsec->GetBinContent(ithCoefficient * 8 + numBinsPerCoefficient - 1) << " +/- " << hAfbAMCATNLOxsec->GetBinError(ithCoefficient * 8 + numBinsPerCoefficient - 1) << std::endl;
      AMCATNLOcoefTable[i][ithCoefficient][0] = hAfbAMCATNLOxsec->GetBinContent(ithCoefficient * 8 + numBinsPerCoefficient - 1);
      // we only have the MC stat uncertainty for this sample
      AMCATNLOcoefTable[i][ithCoefficient][1] = hAfbAMCATNLOxsec->GetBinError(ithCoefficient * 8 + numBinsPerCoefficient - 1);
      AMCATNLOcoefTable[i][ithCoefficient][2] = 0;*/
      AMCATNLOcoefTable[i][ithCoefficient][0] = 0.0;
      AMCATNLOcoefTable[i][ithCoefficient][1] = 0.0;
      AMCATNLOcoefTable[i][ithCoefficient][2] = 0.0;
    }  

    for (int i_bin = 1; i_bin <= nbinsrhoi; ++i_bin) {
      /*hmcatnlo_AllVar[0]->SetBinContent(i_bin+i*nbinsrhoi, AMCATNLOxsectemp->GetBinContent(i_bin) );
      hmcatnlo_AllVar[0]->SetBinError(i_bin+i*nbinsrhoi, AMCATNLOxsectemp->GetBinError(i_bin) );

      hmcatnlo_AllVar[1]->SetBinContent(i_bin+i*nbinsrhoi, AMCATNLOxsectemp->GetBinContent(i_bin)/AMCATNLOxsectemp->Integral() );
      hmcatnlo_AllVar[1]->SetBinError(i_bin+i*nbinsrhoi, AMCATNLOxsectemp->GetBinError(i_bin)/AMCATNLOxsectemp->Integral() );*/

      hmcatnlo_AllVar[0]->SetBinContent(i_bin+i*nbinsrhoi, 0.0);
      hmcatnlo_AllVar[1]->SetBinContent(i_bin+i*nbinsrhoi, 0.0);
    }

    /*TH1D* POWHEGV2HERWIGxsectemp = (TH1D*) ((TH2D*)unffile->Get(fVariableNames[i]+"RespMatSys_POWHEGV2HERWIG"))->ProjectionX();
    MatrixUnf::rebinMultidimensionalInput(POWHEGV2HERWIGxsectemp, binning, generatorBinning_fine->FindNode("ttbargen"), rebinfine);
    TH1D *hAfbPOWHEGV2HERWIGxsec = new TH1D(fVariableNames[i]+"AfbPOWHEGV2HERWIGxsec", "", (numCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
    MatrixUnf::fillAsymsAndCoefficient(POWHEGV2HERWIGxsectemp, nullptr, theoryxsectemp, fVariableNames[i], hAfbPOWHEGV2HERWIGxsec, nullptr, true, true, numCoefficients);

    for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
      std::cout << "POWHEGV2HERWIG coef = " << hAfbPOWHEGV2HERWIGxsec->GetBinContent(ithCoefficient * 8 + numBinsPerCoefficient - 1) << " +/- " << hAfbPOWHEGV2HERWIGxsec->GetBinError(ithCoefficient * 8 + numBinsPerCoefficient - 1) << std::endl;
    }*/


    TH1D* TOP_PTxsectemp = (TH1D*) ((TH2D*)unffile->Get(fVariableNames[i]+"RespMatSys_TOP_PT"))->ProjectionX();
    MatrixUnf::rebinMultidimensionalInput(TOP_PTxsectemp, binning, generatorBinning_fine->FindNode("ttbargen"), rebinfine);
    TH1D* hAfbTOP_PTxsec = new TH1D(fVariableNames[i]+"AfbTOP_PTxsec", "", (numCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
    MatrixUnf::fillAsymsAndCoefficient(TOP_PTxsectemp, nullptr, theoryxsectemp, fVariableNames[i], hAfbTOP_PTxsec, nullptr, true, true, numCoefficients);

    for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
      std::cout << "TOP_PT coef = " << hAfbTOP_PTxsec->GetBinContent(ithCoefficient * 8 + numBinsPerCoefficient - 1) << " +/- " << hAfbTOP_PTxsec->GetBinError(ithCoefficient * 8 + numBinsPerCoefficient - 1) << std::endl;
    }


    // Calculate MC scale uncertainty
    std::vector<TString> ScaleEnvelopeList;
    std::vector<std::vector<double>> ScaleEnvelopeCoefficient;
    for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
      ScaleEnvelopeCoefficient.push_back({});
    }

    ScaleEnvelopeList.push_back("MESCALE_UP");
    ScaleEnvelopeList.push_back("MESCALE_DOWN");
    ScaleEnvelopeList.push_back("MEFACSCALE_UP");
    ScaleEnvelopeList.push_back("MEFACSCALE_DOWN");
    ScaleEnvelopeList.push_back("MERENSCALE_UP");
    ScaleEnvelopeList.push_back("MERENSCALE_DOWN");
    ScaleEnvelopeList.push_back("PSSCALE_WEIGHT_6_UP");
    ScaleEnvelopeList.push_back("PSSCALE_WEIGHT_6_DOWN");
    ScaleEnvelopeList.push_back("PSSCALE_WEIGHT_7_UP");
    ScaleEnvelopeList.push_back("PSSCALE_WEIGHT_7_DOWN");

    TH1D* ScaleEnvelopeHists[ScaleEnvelopeList.size()];
    TH1D* ScaleEnvelopeAfbHists[ScaleEnvelopeList.size()];

    for (unsigned int i_env = 0; i_env < ScaleEnvelopeList.size(); ++i_env) {
      ScaleEnvelopeHists[i_env] = (TH1D *)((TH2D *)unffile->Get(fVariableNames[i] + "RespMatSys_" + ScaleEnvelopeList[i_env].Data()))->ProjectionX();
      MatrixUnf::rebinMultidimensionalInput(ScaleEnvelopeHists[i_env], binning, generatorBinning_fine->FindNode("ttbargen"), rebinfine);
      //ScaleEnvelopeHists[i_env]->Rebin(ScaleEnvelopeHists[i_env]->GetNbinsX() / nbinsrhoi);

      //ScaleEnvelopeAfbHists[i_env] = (TH1D *)theoryxsectemp->Clone(fVariableNames[i] + "ScaleEnvelopeAfbHist_" + ScaleEnvelopeList[i_env].Data());
      ScaleEnvelopeAfbHists[i_env] = new TH1D(fVariableNames[i] + "ScaleEnvelopeAfbHist_" + ScaleEnvelopeList[i_env].Data(), "", (numCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
      MatrixUnf::fillAsymsAndCoefficient(ScaleEnvelopeHists[i_env], nullptr, theoryxsectemp, fVariableNames[i], ScaleEnvelopeAfbHists[i_env], nullptr, true, true, numCoefficients);

      for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
        std::cout << fVariableNames[i] << " ScaleEnvelope " << ScaleEnvelopeList[i_env] << " coef = " << ScaleEnvelopeAfbHists[i_env]->GetBinContent(ithCoefficient * 8 + numBinsPerCoefficient - 1) << " +/- " << ScaleEnvelopeAfbHists[i_env]->GetBinError(ithCoefficient * 8 + numBinsPerCoefficient - 1) << std::endl;

        if (!ScaleEnvelopeList[i_env].Contains("PS")) ScaleEnvelopeCoefficient[ithCoefficient].push_back(ScaleEnvelopeAfbHists[i_env]->GetBinContent(ithCoefficient * 8 + numBinsPerCoefficient - 1) - nominal_coef[ithCoefficient]);
      }
    }

    for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
      ScaleEnvelopeCoefficient[ithCoefficient].push_back(0.);

      std::sort(ScaleEnvelopeCoefficient[ithCoefficient].begin(), ScaleEnvelopeCoefficient[ithCoefficient].end());
      std::cout << fVariableNames[i] << " ScaleEnvelope: " << ScaleEnvelopeCoefficient[ithCoefficient].at(ScaleEnvelopeCoefficient[ithCoefficient].size() - 1) << " , " << ScaleEnvelopeCoefficient[ithCoefficient].at(0) << std::endl;

      MCcoefTable[i][ithCoefficient][0] = nominal_coef[ithCoefficient];
      // sum the scale envelope in quadrature with the MC stat uncertainty to define the final MCstat+scale uncertainties
      MCcoefTable[i][ithCoefficient][1] = sqrt(nominal_coef_staterr[ithCoefficient] * nominal_coef_staterr[ithCoefficient] + ScaleEnvelopeCoefficient[ithCoefficient].at(ScaleEnvelopeCoefficient[ithCoefficient].size() - 1) * ScaleEnvelopeCoefficient[ithCoefficient].at(ScaleEnvelopeCoefficient[ithCoefficient].size() - 1));
      MCcoefTable[i][ithCoefficient][2] = sqrt(nominal_coef_staterr[ithCoefficient] * nominal_coef_staterr[ithCoefficient] + ScaleEnvelopeCoefficient[ithCoefficient].at(0) * ScaleEnvelopeCoefficient[ithCoefficient].at(0));
    }

    std::vector<std::vector<double>> theoryCoefficient;
    for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
        theoryCoefficient.push_back({});
        theoryCoefficient[ithCoefficient].clear();
    }

    if (fVariableNames[i].Contains("ll_cLab") || fVariableNames[i].Contains("llbar_delta_phi") || fVariableNames[i].Contains("llbar_delta_eta")) {

      TH1D* theoryAfbHists[12];
      for (int i_th = 0; i_th < 12; ++i_th)
      {

        NLOW_histo[i_th]->Scale(1./NLOW_histo[i_th]->Integral());

        //theoryAfbHists[i_th] = (TH1D*) NLOW_histo[i_th]->Clone(fVariableNames[i]+"theoryAfbHist_"+(i_th>=6?"LOxsec":"NLOWxsec")+TheoryDataFileList[i_th%6].Data());
        theoryAfbHists[i_th] = new TH1D(fVariableNames[i] + "theoryAfbHist_" + (i_th >= 6 ? "LOxsec" : "NLOWxsec") + TheoryDataFileList[i_th % 6].Data(), "", (numCoefficients * (numBinsPerCoefficient + 2)) - 2, 0, 1);
        MatrixUnf::fillAsymsAndCoefficient(NLOW_histo[i_th], nullptr, theoryxsectemp, fVariableNames[i], theoryAfbHists[i_th], nullptr, true, true, 1);
        std::cout<<fVariableNames[i]<<" theory " << (i_th>=6?"LOxsec":"NLOWxsec") << TheoryDataFileList[i_th%6] << " coef = " << theoryAfbHists[i_th]->GetBinContent(nbinsrhoi-1)<<" +/- "<<theoryAfbHists[i_th]->GetBinError(nbinsrhoi-1)<<std::endl;

        for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
          std::cout << fVariableNames[i] << " theory " << (i_th >= 6 ? "LOxsec" : "NLOWxsec") << TheoryDataFileList[i_th % 6] << " coef = " << theoryAfbHists[i_th]->GetBinContent(ithCoefficient * 8 + numBinsPerCoefficient - 1) << " +/- " << theoryAfbHists[i_th]->GetBinError(ithCoefficient * 8 + numBinsPerCoefficient - 1) << std::endl;

          // TODO: fix this since it is always grabbing 0th coefficient but that's because NLO calculation only has a single NLO value currently
          theoryCoefficient[ithCoefficient].push_back(theoryAfbHists[i_th]->GetBinContent(0 * 8 + numBinsPerCoefficient - 1));
        }
      }
    }

    for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
      Meas_BinSums(i, ithCoefficient) = binsumTable[0][i][ithCoefficient];  // take the measured BinSum from the absolute diff xsec

      Meas_Coefs_abs(i, ithCoefficient) = coefTable[0][i][ithCoefficient];  // take the measured coefficient from the absolute diff xsec

      Meas_Coefs(i, ithCoefficient) = coefTable[1][i][ithCoefficient];  // take the measured coefficient from the normalised diff xsec
      Theory_Coefs(i, ithCoefficient) = hAfbtheoryxsec->GetBinContent(nbinsrhoi - 1);
      //AMCATNLO_Coefs(i, ithCoefficient) = hAfbAMCATNLOxsec->GetBinContent(nbinsrhoi - 1);
      AMCATNLO_Coefs(i, ithCoefficient) = 0.0;
      //POWHEGV2HERWIG_Coefs(i, ithCoefficient) = hAfbPOWHEGV2HERWIGxsec->GetBinContent(nbinsrhoi - 1);
      POWHEGV2HERWIG_Coefs(i, ithCoefficient) = 0.0;
      TOP_PT_Coefs(i, ithCoefficient) = hAfbTOP_PTxsec->GetBinContent(nbinsrhoi - 1);

      if (i < 35) {
        NLOW_Coefs(i, ithCoefficient) = NLOWvalues[i];
        NLOWuncorr_Coefs(i, ithCoefficient) = 0;
        NLOW3_Coefs(i, ithCoefficient) = coefTable[0][i][ithCoefficient];
        NLOW4_Coefs(i, ithCoefficient) = coefTable[1][i][ithCoefficient];
      } else {
        NLOW_Coefs(i, ithCoefficient) = theoryCoefficient[ithCoefficient][0];
        NLOWuncorr_Coefs(i, ithCoefficient) = theoryCoefficient[ithCoefficient][3];
        NLOW3_Coefs(i, ithCoefficient) = coefTable[0][i][ithCoefficient];
        NLOW4_Coefs(i, ithCoefficient) = coefTable[1][i][ithCoefficient];
      }
    }

    // Modified :
    // NLOW_Coefs(i,0) = NLOWvalues[i];
    // NLOWuncorr_Coefs(i,0) = 0;
    // NLOW3_Coefs(i,0) = coefTable[0][i];
    // NLOW4_Coefs(i,0) = coefTable[1][i];
    // End

    if(i==fVariableNames.size()-1) {
      std::cout<<"measured and correlated and uncorrelated NLOW coefficient vectors: "<<std::endl;
      Meas_Coefs.Print("f=%1.5g ");
      NLOW_Coefs.Print("f=%1.5g ");
      NLOWuncorr_Coefs.Print("f=%1.5g ");
      std::cout<<"measured BinSum vector: "<<std::endl;
      Meas_BinSums.Print("f=%1.5g "); 
    }

    for (int inorm = 0; inorm <= 1; ++inorm)
    {

      TCanvas *c2 = new TCanvas("c2", "c2");
      // TLegend *l2  = new TLegend();
      // TLegend *l2N = new TLegend();
      // TLegend *l2N = new TLegend(0.6 , 0.75,   0.7, 0.85,NULL,"brNDC");

      // Amandeep : Changing for 2D variables 
      float xlow_legend  = 0.2;
      float xhigh_legend = 0.3;
      //if (numDimensions == 2 )  xlow_legend = 0.56; xhigh_legend = 0.66;
      TLegend *l2  = new TLegend(xlow_legend, 0.775, xhigh_legend, 0.875,NULL,"brNDC");
      TLegend *l2N = new TLegend(xlow_legend, 0.775, xhigh_legend, 0.875,NULL,"brNDC");

      TH1D* theoryxsec = (TH1D*) theoryxsectemp->Clone(fVariableNames[i]+"TheoryXsec");

      if (numDimensions == 1) {
	      theoryxsec->GetXaxis()->SetNdivisions(-20604);
        theoryxsec->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
        theoryxsec->GetYaxis()->SetTitle(Form("#frac{d#sigma}{d%s} [pb]",quantity.Data()));
        if(inorm) theoryxsec->GetYaxis()->SetTitle(Form("#frac{1}{#sigma} #frac{d#sigma}{d%s}",quantity.Data()));

        if( MatrixUnf::isCtype(fVariableNames[i]) || MatrixUnf::isCPMtype(fVariableNames[i]) ) {
          theoryxsec->GetYaxis()->SetTitle(Form("#frac{d#sigma}{d%s} [pb]",quantity.Data()));
          if(inorm) theoryxsec->GetYaxis()->SetTitle(Form("#frac{1}{#sigma} #frac{d#sigma}{d%s}",quantity.Data()));
        }
      } else if (numDimensions == 2) {
        c2->SetGridx();
        theoryxsec->GetXaxis()->SetNdivisions(-20604);
        theoryxsec->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
        theoryxsec->GetYaxis()->SetTitleSize(0.03);
        theoryxsec->GetYaxis()->SetTitleOffset(2.8);
        if(fVariableNames[i].Contains("_mttbar")) {
          theoryxsec->GetYaxis()->SetTitle(Form("#frac{d^{2}#sigma}{d(%s) dm_{t#bar{t}}} [#frac{pb}{GeV}]",quantity.Data()));
          if(inorm) theoryxsec->GetYaxis()->SetTitle(Form("#frac{1}{#sigma} #frac{d^{2}#sigma}{d(%s) dm_{t#bar{t}}}",quantity.Data()));
          if( MatrixUnf::isCtype(fVariableNames[i]) || MatrixUnf::isCPMtype(fVariableNames[i]) ) {
            theoryxsec->GetYaxis()->SetTitle(Form("#frac{d^{2}#sigma}{d(%s) dm_{t#bar{t}}} [#frac{pb}{GeV}]",quantity.Data()));
            if(inorm) theoryxsec->GetYaxis()->SetTitle(Form("#frac{1}{#sigma} #frac{d^{2}#sigma}{d(%s) dm_{t#bar{t}}}",quantity.Data()));
          }
        }

        if( MatrixUnf::isCtype(fVariableNames[i]) || MatrixUnf::isCPMtype(fVariableNames[i]) ) {
          theoryxsec->GetYaxis()->SetTitle(Form("#frac{d^{2}#sigma}{d(%s) dm_{t#bar{t}}} [#frac{pb}{GeV}]",quantity.Data()));
          if(inorm) theoryxsec->GetYaxis()->SetTitle(Form("#frac{1}{#sigma} #frac{d^{2}#sigma}{d(%s) dm_{t#bar{t}}}",quantity.Data()));
        }
      } 
      // End

      /*TH1D* AMCATNLOxsec = (TH1D*) AMCATNLOxsectemp->Clone(fVariableNames[i]+"AMCATNLOxsec");
      AMCATNLOxsec->Scale(1./AMCATNLOxsec->Integral());
      AMCATNLOxsec->Scale(theoryxsec->Integral());

      TH1D* POWHEGV2HERWIGxsec = (TH1D*) POWHEGV2HERWIGxsectemp->Clone(fVariableNames[i]+"POWHEGV2HERWIGxsec");
      POWHEGV2HERWIGxsec->Scale(1./POWHEGV2HERWIGxsec->Integral());
      POWHEGV2HERWIGxsec->Scale(theoryxsec->Integral());*/

      TH1D* TOP_PTxsec = (TH1D*) TOP_PTxsectemp->Clone(fVariableNames[i]+"TOP_PTxsec");
      TOP_PTxsec->Scale(1./TOP_PTxsec->Integral());
      TOP_PTxsec->Scale(theoryxsec->Integral());

      TH1D* tunfresultcor = (TH1D*) ((TH1D*)unffile->Get(fVariableNames[i]+"TUnfResultCor_rebinnedA"))->Clone(fVariableNames[i]+"TUnfResultCor_rebinnedA_clone");
      tunfresultcor->SetLineColor(kBlack);
      tunfresultcor->SetLineWidth(2);
      tunfresultcor->SetMarkerColor(kBlack);
      tunfresultcor->SetMarkerSize(0.5);
      tunfresultcor->SetMarkerStyle(21);
      tunfresultcor->SetTitle("");
      tunfresultcor->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
      SetAxis(fVariableNames[i], tunfresultcor->GetXaxis());

      // tunfresultcor->GetYaxis()->SetTitle(Form("d#sigma/d%s [pb]",quantity.Data()));
      // tunfresultcorNorm->GetYaxis()->SetTitle(Form("1/#sigma d#sigma/d%s",quantity.Data()));


      TMatrixD Measvec   = MatrixUnf::TH1toTMatrixD(tunfresultcor);
      TMatrixD Theoryvec = MatrixUnf::TH1toTMatrixD(theoryxsec);
      // TMatrixD AMCATNLOvec = MatrixUnf::TH1toTMatrixD(AMCATNLOxsec);
      // TMatrixD POWHEGV2HERWIGvec = MatrixUnf::TH1toTMatrixD(POWHEGV2HERWIGxsec);
      // TMatrixD TOP_PTvec  = MatrixUnf::TH1toTMatrixD(TOP_PTxsec);

      TMatrixD NLOWvec       = MatrixUnf::TH1toTMatrixD(tempNLOWxsec);
      TMatrixD NLOWuncorrvec = MatrixUnf::TH1toTMatrixD(tempNLOWxsecuncorr);
      TMatrixD NLOW3vec      = MatrixUnf::TH1toTMatrixD(tempNLOWxsec3);
      TMatrixD NLOW4vec      = MatrixUnf::TH1toTMatrixD(tempNLOWxsec4);

      std::vector<TMatrixD> NLOWvecMulti;
      for (int iM = 0; iM < nFit; ++iM)
      {
        TMatrixD NLOWvecMultitemp = MatrixUnf::TH1toTMatrixD(tempNLOWxsecMulti[iM]);
        NLOWvecMulti.push_back(NLOWvecMultitemp);
      }

      // Amandeep : this might have had to do with certain lab frame variables 
      // if(i<20){
      if(i<39){
      }
      else {
        NLOWvec       = MatrixUnf::TH1toTMatrixD(NLOW_histo[0]);
        NLOWuncorrvec = MatrixUnf::TH1toTMatrixD(NLOW_histo[3]);
        NLOW3vec      = MatrixUnf::TH1toTMatrixD(NLOW_histo[0]);
        if (fVariableNames[i].Contains("llbar_delta_phi") || fVariableNames[i].Contains("llbar_delta_eta")) NLOW4vec = MatrixUnf::TH1toTMatrixD(NNLO_histo[0]);
        else NLOW4vec = MatrixUnf::TH1toTMatrixD(NLOW_histo[0]);

        for (int iM = 0; iM < nFit; ++iM) {
          NLOWvecMulti[iM] = MatrixUnf::TH1toTMatrixD(NLOW_histo[0]);
        }
      }

      NLOWvec       *= theoryxsectemp->Integral();
      NLOWuncorrvec *= theoryxsectemp->Integral();
      NLOW3vec *= theoryxsectemp->Integral();
      NLOW4vec *= theoryxsectemp->Integral();

      for (int iM = 0; iM < nFit; ++iM)
      {
        NLOWvecMulti[iM]*=theoryxsectemp->Integral();
      }

      const int binToIgnore_const = 1;

      if(inorm==0) {
        Theoryvec     *= tunfresultcor->Integral()/theoryxsec->Integral();
        // AMCATNLOvec       *= tunfresultcor->Integral()/theoryxsec->Integral();
        // POWHEGV2HERWIGvec *= tunfresultcor->Integral()/theoryxsec->Integral();
        // TOP_PTvec         *= tunfresultcor->Integral()/theoryxsec->Integral();
        NLOWvec       *= tunfresultcor->Integral()/theoryxsec->Integral();
        NLOWuncorrvec *= tunfresultcor->Integral()/theoryxsec->Integral();
        NLOW3vec      *= tunfresultcor->Integral()/theoryxsec->Integral();
        NLOW4vec      *= tunfresultcor->Integral()/theoryxsec->Integral();

        for (int iM = 0; iM < nFit; ++iM) {
          NLOWvecMulti[iM] *= tunfresultcor->Integral()/theoryxsec->Integral();
        }

        for (int i_bin = 0; i_bin < nbinsrhoi; ++i_bin)
        {
          Measvec_allVars(nbinsrhoi*i + i_bin,0)   = Measvec(i_bin,0);
          Theoryvec_allVars(nbinsrhoi*i + i_bin,0) = Theoryvec(i_bin,0);
          // AMCATNLOvec_allVars(nbinsrhoi*i + i_bin,0)       = AMCATNLOvec(i_bin,0);
          // POWHEGV2HERWIGvec_allVars(nbinsrhoi*i + i_bin,0) = POWHEGV2HERWIGvec(i_bin,0);
          // TOP_PTvec_allVars(nbinsrhoi*i + i_bin,0)         = TOP_PTvec(i_bin,0);
          NLOWvec_allVars(nbinsrhoi*i + i_bin,0)       = NLOWvec(i_bin,0);
          NLOWuncorrvec_allVars(nbinsrhoi*i + i_bin,0) = NLOWuncorrvec(i_bin,0);
          NLOW3vec_allVars(nbinsrhoi*i + i_bin,0)      = NLOW3vec(i_bin,0);
          NLOW4vec_allVars(nbinsrhoi*i + i_bin,0)      = NLOW4vec(i_bin,0);

          for (int iM = 0; iM < nFit; ++iM) {
            (NLOWvecMulti_allVars[iM])(nbinsrhoi*i + i_bin,0) = (NLOWvecMulti[iM])(i_bin,0);
          }

        }

      }
      else {
        Measvec   *= 1./tunfresultcor->Integral();
        Theoryvec *= 1./theoryxsec->Integral();
        // AMCATNLOvec*=1./theoryxsec->Integral();
        // POWHEGV2HERWIGvec*=1./theoryxsec->Integral();
        // TOP_PTvec*=1./theoryxsec->Integral();
        NLOWvec       *= 1./theoryxsec->Integral();
        NLOWuncorrvec *= 1./theoryxsec->Integral();
        NLOW3vec *= 1./theoryxsec->Integral();
        NLOW4vec *= 1./theoryxsec->Integral();

        for (int iM = 0; iM < nFit; ++iM)
        {
          NLOWvecMulti[iM]*=1./theoryxsec->Integral();
        }

        TMatrixD Measvec_temp = Measvec;
        Measvec.ResizeTo(nbinsrhoi-1,1);
        TMatrixD Theoryvec_temp = Theoryvec;
        Theoryvec.ResizeTo(nbinsrhoi-1,1);
        
        // TMatrixD AMCATNLOvec_temp = AMCATNLOvec;
        // AMCATNLOvec.ResizeTo(nbinsrhoi-1,1);
        // TMatrixD POWHEGV2HERWIGvec_temp = POWHEGV2HERWIGvec;
        // POWHEGV2HERWIGvec.ResizeTo(nbinsrhoi-1,1);
        // TMatrixD TOP_PTvec_temp = TOP_PTvec;
        // TOP_PTvec.ResizeTo(nbinsrhoi-1,1);
        
        TMatrixD NLOWvec_temp = NLOWvec;
        NLOWvec.ResizeTo(nbinsrhoi-1,1);
        TMatrixD NLOWuncorrvec_temp = NLOWuncorrvec;
        NLOWuncorrvec.ResizeTo(nbinsrhoi-1,1);
        TMatrixD NLOW3vec_temp = NLOW3vec;
        NLOW3vec.ResizeTo(nbinsrhoi-1,1);
        TMatrixD NLOW4vec_temp = NLOW4vec;
        NLOW4vec.ResizeTo(nbinsrhoi-1,1);

        std::vector<TMatrixD> NLOWvecMulti_temp;
        for (int iM = 0; iM < nFit; ++iM)
        {
          TMatrixD NLOWvecMulti_tempcopy = NLOWvecMulti[iM];
          NLOWvecMulti_temp.push_back(NLOWvecMulti_tempcopy);
          NLOWvecMulti[iM].ResizeTo(nbinsrhoi-1,1);
        }

        const int nbin_total  = (nbinsrhoi);
        const int binToIgnore = binToIgnore_const;

        int iShp = 0;

        for (int iAbs = 0; iAbs < nbin_total; ++iAbs) {
          if (iAbs % nbinsrhoi == binToIgnore) continue;

          Measvec(iShp,0)   = Measvec_temp(iAbs,0);
          Theoryvec(iShp,0) = Theoryvec_temp(iAbs,0);
          // AMCATNLOvec(iShp,0)       = AMCATNLOvec_temp(iAbs,0);
          // POWHEGV2HERWIGvec(iShp,0) = POWHEGV2HERWIGvec_temp(iAbs,0);
          // TOP_PTvec(iShp,0)  = TOP_PTvec_temp(iAbs,0);
          NLOWvec(iShp,0)       = NLOWvec_temp(iAbs,0);
          NLOWuncorrvec(iShp,0) = NLOWuncorrvec_temp(iAbs,0);
          NLOW3vec(iShp,0) = NLOW3vec_temp(iAbs,0);
          NLOW4vec(iShp,0) = NLOW4vec_temp(iAbs,0);

          for (int iM = 0; iM < nFit; ++iM)
          {
            (NLOWvecMulti[iM])(iShp,0) = (NLOWvecMulti_temp[iM])(iAbs,0);
          }

          ++iShp;
        }


        for (int i_bin = 0; i_bin < (nbinsrhoi-1); ++i_bin)
        {
          Measvec_norm_allVars((nbinsrhoi-1)*i + i_bin,0)   = Measvec(i_bin,0);
          Theoryvec_norm_allVars((nbinsrhoi-1)*i + i_bin,0) = Theoryvec(i_bin,0);

          // AMCATNLOvec_norm_allVars((nbinsrhoi-1)*i + i_bin,0)       = AMCATNLOvec(i_bin,0);
          // POWHEGV2HERWIGvec_norm_allVars((nbinsrhoi-1)*i + i_bin,0) = POWHEGV2HERWIGvec(i_bin,0);
          // TOP_PTvec_norm_allVars((nbinsrhoi-1)*i + i_bin,0) = TOP_PTvec(i_bin,0);

          NLOWvec_norm_allVars((nbinsrhoi-1)*i + i_bin,0)       = NLOWvec(i_bin,0);
          NLOWuncorrvec_norm_allVars((nbinsrhoi-1)*i + i_bin,0) = NLOWuncorrvec(i_bin,0);
          NLOW3vec_norm_allVars((nbinsrhoi-1)*i + i_bin,0) = NLOW3vec(i_bin,0);
          NLOW4vec_norm_allVars((nbinsrhoi-1)*i + i_bin,0) = NLOW4vec(i_bin,0);

          for (int iM = 0; iM < nFit; ++iM){
            (NLOWvecMulti_norm_allVars[iM])((nbinsrhoi-1)*i + i_bin,0) = (NLOWvecMulti[iM])(i_bin,0);
          }
        }
      }


      for (int i_unctype = 0; i_unctype <= 1; ++i_unctype) {

        TMatrixD CovMatInv = inorm ? (i_unctype ? StatCovMats.at(1) : TotalCovMats.at(1)) : (i_unctype ? StatCovMats.at(0) : TotalCovMats.at(0));
        if(inorm) {

          const int nbin_shape = (nbinsrhoi-1);
          const int nbin_total = (nbinsrhoi);

          std::vector<TMatrixD> covmats_5bins;

          for (int i_ignore = 0; i_ignore < nbinsrhoi; ++i_ignore)
          {

            const int binToIgnore = i_ignore;
            TMatrixD covmat_5bins(nbin_shape, nbin_shape);

            int iShpR = 0, iShpC = 0;

            for (int iAbsR = 0; iAbsR < nbin_total; ++iAbsR) {
              if (iAbsR % nbinsrhoi == binToIgnore) continue;

              iShpC = 0;
              for (int iAbsC = 0; iAbsC < nbin_total; ++iAbsC) {
                if (iAbsC % nbinsrhoi == binToIgnore) continue;

                covmat_5bins(iShpR, iShpC) = CovMatInv(iAbsR, iAbsC);
                ++iShpC;
              }
              ++iShpR;
            }
            covmats_5bins.push_back(covmat_5bins);
          }
          CovMatInv.ResizeTo(nbinsrhoi-1,nbinsrhoi-1);
          CovMatInv = covmats_5bins.at(binToIgnore_const);
        }


        TMatrixD CovMat = CovMatInv;
        CovMatInv.Invert();

        // Amandeep : Main unfolding Chi2
        TMatrixD dMeas  = Measvec-Theoryvec;
        TMatrixD dMeasT = dMeas;
        dMeasT          = dMeasT.T();
        TMatrixD MChi2  = dMeasT*CovMatInv*dMeas;

        std::cout<<fVariableNames[i]<<(inorm?"Norm":"")<<(i_unctype ? " Stat" : " Total" )<< " Chi2 for powheg = "<<MChi2(0,0)<<std::endl;
        Chi2values[inorm][i_unctype][i][0] = MChi2(0,0);


        //calculate chi2s for data vs various predictions
        if(kTRUE) {

          // TMatrixD dMeasAMCATNLO  = Measvec-AMCATNLOvec;
          // TMatrixD dMeasAMCATNLOT = dMeasAMCATNLO;
          // dMeasAMCATNLOT          = dMeasAMCATNLOT.T();
          // TMatrixD MChi2AMCATNLO  = dMeasAMCATNLOT*CovMatInv*dMeasAMCATNLO;

          // std::cout<<fVariableNames[i]<<(inorm?"Norm":"")<<(i_unctype ? " Stat" : " Total" )<< " Chi2 for AMCATNLO = "<<MChi2AMCATNLO(0,0)<<std::endl;
          // Chi2values[inorm][i_unctype][i][1] = MChi2AMCATNLO(0,0);


          // TMatrixD dMeasPOWHEGV2HERWIG  = Measvec-POWHEGV2HERWIGvec;
          // TMatrixD dMeasPOWHEGV2HERWIGT = dMeasPOWHEGV2HERWIG;
          // dMeasPOWHEGV2HERWIGT          = dMeasPOWHEGV2HERWIGT.T();
          // TMatrixD MChi2POWHEGV2HERWIG  = dMeasPOWHEGV2HERWIGT*CovMatInv*dMeasPOWHEGV2HERWIG;

          // std::cout<<fVariableNames[i]<<(inorm?"Norm":"")<<(i_unctype ? " Stat" : " Total" )<< " Chi2 for POWHEGV2HERWIG = "<<MChi2POWHEGV2HERWIG(0,0)<<std::endl;
          // Chi2values[inorm][i_unctype][i][2] = MChi2POWHEGV2HERWIG(0,0);


          // TMatrixD dMeasTOP_PT = Measvec-TOP_PTvec;
          // TMatrixD dMeasTOP_PTT = dMeasTOP_PT;
          // dMeasTOP_PTT = dMeasTOP_PTT.T();
          // TMatrixD MChi2TOP_PT = dMeasTOP_PTT*CovMatInv*dMeasTOP_PT;

          // std::cout<<fVariableNames[i]<<(inorm?"Norm":"")<<(i_unctype ? " Stat" : " Total" )<< " Chi2 for TOP_PT = "<<MChi2TOP_PT(0,0)<<std::endl;
          // Chi2values[inorm][i_unctype][i][3] = MChi2TOP_PT(0,0);


          TMatrixD dMeasNLOW  = Measvec-NLOWvec;
          TMatrixD dMeasNLOWT = dMeasNLOW;
          dMeasNLOWT          = dMeasNLOWT.T();
          TMatrixD MChi2NLOW  = dMeasNLOWT*CovMatInv*dMeasNLOW;

          std::cout<<fVariableNames[i]<<(inorm?"Norm":"")<<(i_unctype ? " Stat" : " Total" )<< " Chi2 for NLOW = "<<MChi2NLOW(0,0)<<std::endl;
          Chi2values[inorm][i_unctype][i][4] = MChi2NLOW(0,0);


          TMatrixD dMeasNLOWuncorr  = Measvec-NLOWuncorrvec;
          TMatrixD dMeasNLOWuncorrT = dMeasNLOWuncorr;
          dMeasNLOWuncorrT          = dMeasNLOWuncorrT.T();
          TMatrixD MChi2NLOWuncorr  = dMeasNLOWuncorrT*CovMatInv*dMeasNLOWuncorr;

          std::cout<<fVariableNames[i]<<(inorm?"Norm":"")<<(i_unctype ? " Stat" : " Total" )<< " Chi2 for NLOWuncorr = "<<MChi2NLOWuncorr(0,0)<<std::endl;
          Chi2values[inorm][i_unctype][i][5] = MChi2NLOWuncorr(0,0);


          TMatrixD dMeasNLOW3 = Measvec-NLOW3vec;
          TMatrixD dMeasNLOW3T = dMeasNLOW3;
          dMeasNLOW3T = dMeasNLOW3T.T();
          TMatrixD MChi2NLOW3 = dMeasNLOW3T*CovMatInv*dMeasNLOW3;

          std::cout<<fVariableNames[i]<<(inorm?"Norm":"")<<(i_unctype ? " Stat" : " Total" )<< " Chi2 for NLOW3 = "<<MChi2NLOW3(0,0)<<std::endl;
          Chi2values[inorm][i_unctype][i][6] = MChi2NLOW3(0,0);

          TMatrixD dMeasNLOW4  = Measvec-NLOW4vec;
          TMatrixD dMeasNLOW4T = dMeasNLOW4;
          dMeasNLOW4T          = dMeasNLOW4T.T();
          TMatrixD MChi2NLOW4  = dMeasNLOW4T*CovMatInv*dMeasNLOW4;

          std::cout<<fVariableNames[i]<<(inorm?"Norm":"")<<(i_unctype ? " Stat" : " Total" )<< " Chi2 for NLOW4 = "<<MChi2NLOW4(0,0)<<std::endl;
          Chi2values[inorm][i_unctype][i][7] = MChi2NLOW4(0,0);
        }


        // Analytically calculate Chi2 fit result for C (or f_SM in the case of the lab frame variables) by converting each bin into a measure of C (or f_SM), then combine the bins

        TMatrixD CMeasvec = Measvec;
        // Amandeep : Modifying this
        // TMatrixD CNLOWvec = i<20 ? NLOW4vec : NLOWvec; //Use NLOW4vec instead of NLOWvec to avoid dividing by zero below for the c_M variables
        TMatrixD CNLOWvec    = i<39 ? NLOW4vec : NLOWvec;
        // End
        CMeasvec -= NLOWuncorrvec;
        CNLOWvec -= NLOWuncorrvec;
        //std::cout<<"Checking CMeasvec:"<<std::endl;
        //CMeasvec.Print("f=%1.5g ");
        //CNLOWvec.Print("f=%1.5g ");
        if(i<39) CNLOWvec *= 1./(coefTable[1][i][0]-NLOWvaluesUnCorr[i]);
        //CNLOWvec.Print("f=%1.5g ");
        CMeasvec /= TMatrixDColumn(CNLOWvec,0);
        //CMeasvec.Print("f=%1.5g ");

        //CovMat.Print("f=%1.5g ");
        CovMat /= TMatrixDColumn(CNLOWvec,0);
        CovMat /= TMatrixDRow(CNLOWvec.T(),0);
        //CovMat.Print("f=%1.5g ");
        TMatrixD CovMat_copy = CovMat;
        CovMat.Invert();

        TMatrixD unit(Measvec.GetNrows(),1);
        for(int i_row=0;i_row<Measvec.GetNrows();i_row++) {
          unit(i_row,0) = 1;
        }


        TMatrixD Cweights = CovMat*unit;
        double sum_weight = 0.;
        for(int i_row=0;i_row<Measvec.GetNrows();i_row++) {
          sum_weight += Cweights(i_row,0);
        }
        Cweights *= 1./sum_weight;
        std::cout<<"CMeasvec weights:"<<std::endl;
        Cweights.Print("f=%1.5g ");


        //calculate Chi2 fit result 
        TMatrixD measCmat = (CMeasvec.T())*Cweights;
        std::cout<<"measCval: "<<measCmat(0,0)<<std::endl;
        //measCmat.Print("f=%1.6g ");

        TMatrixD CweightsT = Cweights;
        CweightsT = CweightsT.T();

        TMatrixD measCvar = CweightsT*CovMat_copy*Cweights;
        std::cout<<"measCvar: "<<sqrt(measCvar(0,0))<<std::endl;
        //measCvar.Print("f=%1.6g ");

      }

      if(inorm) tunfresultcor->Scale(1./tunfresultcor->Integral());
      for (int i_bin = 0; i_bin < tunfresultcor->GetNbinsX(); ++i_bin)
      {
        tunfresultcor->SetBinError( i_bin+1, sqrt(TotalCovMats.at(inorm)(i_bin,i_bin))  );
      }
      // Original :
      // if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(tunfresultcor);
      // Modified :
      if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(tunfresultcor, fVariableNames[i]);
      // End

      if(inorm) theoryxsec->Scale(1./theoryxsec->Integral());
      // Original :
      // if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(theoryxsec);
      // Modified :
      if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(theoryxsec, fVariableNames[i]);
      // End

      /*if(inorm) AMCATNLOxsec->Scale(1./AMCATNLOxsec->Integral());
      if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(AMCATNLOxsec, fVariableNames[i]);

      if(inorm) POWHEGV2HERWIGxsec->Scale(1./POWHEGV2HERWIGxsec->Integral());
      if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(POWHEGV2HERWIGxsec, fVariableNames[i]);*/

      if(inorm) TOP_PTxsec->Scale(1./TOP_PTxsec->Integral());
      if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(TOP_PTxsec, fVariableNames[i]);

      gSystem->mkdir("yaml", true);
      FILE *yamlfile = fopen(Form("yaml/predicted%sdiffxsec_%s.yaml",(inorm?"norm":""),fVariableNames[i].Data()), "w");
      fprintf(yamlfile, "dependent_variables:\n");
      if(inorm) fprintf(yamlfile, "- header: {name: '$\\frac{1}{\\sigma} \\frac{{d}\\sigma}{{d}%s}$'}\n", (MatrixUnf::ObsName(fVariableNames[i],kTRUE)).Data() );
      else fprintf(yamlfile, "- header: {name: '$\\frac{{d}\\sigma}{{d}%s}$ [pb]'}\n", (MatrixUnf::ObsName(fVariableNames[i],kTRUE)).Data() );
      fprintf(yamlfile, "  qualifiers:\n");
      fprintf(yamlfile, "  - {name: prediction, value: POWHEGV2}\n");
      fprintf(yamlfile, "  values:\n");
      for (int i_bin = 0; i_bin < nbinsrhoi; ++i_bin)
      {
        fprintf(yamlfile, "  - errors:\n");
        fprintf(yamlfile, "    - {label: stat, symerror: %2.6f}\n", theoryxsec->GetBinError(i_bin+1) );
        fprintf(yamlfile, "    value: %2.6f\n", theoryxsec->GetBinContent(i_bin+1) );
      }
      if(inorm) fprintf(yamlfile, "- header: {name: '$\\frac{1}{\\sigma} \\frac{{d}\\sigma}{{d}%s}$'}\n", (MatrixUnf::ObsName(fVariableNames[i],kTRUE)).Data() );
      else fprintf(yamlfile, "- header: {name: '$\\frac{{d}\\sigma}{{d}%s}$ [pb]'}\n", (MatrixUnf::ObsName(fVariableNames[i],kTRUE)).Data() );
      fprintf(yamlfile, "  qualifiers:\n");
      fprintf(yamlfile, "  - {name: prediction, value: MG5_aMC@NLO}\n");
      fprintf(yamlfile, "  values:\n");
      for (int i_bin = 0; i_bin < nbinsrhoi; ++i_bin)

      // {
      //   fprintf(yamlfile, "  - errors:\n");
      //   fprintf(yamlfile, "    - {label: stat, symerror: %2.6f}\n", AMCATNLOxsec->GetBinError(i_bin+1) );
      //   fprintf(yamlfile, "    value: %2.6f\n", AMCATNLOxsec->GetBinContent(i_bin+1) );
      // }

      /*
      if(inorm) fprintf(yamlfile, "- header: {name: '$\\frac{1}{\\sigma} \\frac{{d}\\sigma}{{d}%s}$'}\n", (MatrixUnf::ObsName(fVariableNames[i],kTRUE)).Data() );
      else fprintf(yamlfile, "- header: {name: '$\\frac{{d}\\sigma}{{d}%s}$ [pb]'}\n", (MatrixUnf::ObsName(fVariableNames[i],kTRUE)).Data() );
      fprintf(yamlfile, "  qualifiers:\n");
      fprintf(yamlfile, "  - {name: prediction, value: data}\n");
      fprintf(yamlfile, "  values:\n");
      for (int i_bin = 0; i_bin < nbinsrhoi; ++i_bin)
      {
        fprintf(yamlfile, "  - errors:\n");
        fprintf(yamlfile, "    - {label: total, symerror: %2.6f}\n", tunfresultcor->GetBinError(i_bin+1) );
        fprintf(yamlfile, "    value: %2.6f\n", tunfresultcor->GetBinContent(i_bin+1) );
      }
      */

      fprintf(yamlfile, "independent_variables:\n");
      fprintf(yamlfile, "- header: {name: '$%s$'}\n", (MatrixUnf::ObsName(fVariableNames[i],kTRUE)).Data() );
      fprintf(yamlfile, "  values:\n");

      for (int i_bin = 0; i_bin < nbinsrhoi; ++i_bin)
      {
        fprintf(yamlfile, "  - {high: %2.6f, low: %2.6f}\n", hbinlowedge_AllVar->GetBinContent(i*nbinsrhoi+i_bin+1)+hbinwidth_AllVar->GetBinContent(i*nbinsrhoi+i_bin+1), hbinlowedge_AllVar->GetBinContent(i*nbinsrhoi+i_bin+1) );
      }
      fclose(yamlfile);

      // for (int i_bin = 0; i_bin < nbinsrhoi; ++i_bin)
      // {
      //   AMCATNLOxsec->SetBinError(i_bin+1, 1e-10);
      //   POWHEGV2HERWIGxsec->SetBinError(i_bin+1, 1e-10);
      //   TOP_PTxsec->SetBinError(i_bin+1, 1e-10);
      // }

      /*
      tunfresultcor->GetYaxis()->SetRangeUser(theoryxsec->GetMaximum()*0.7,theoryxsec->GetMaximum()*1.3 );
      if( theoryxsec->GetMaximum()*0.7 > theoryxsec->GetMinimum() ) tunfresultcor->GetYaxis()->SetRangeUser(0,theoryxsec->GetMaximum()*1.3 );
      if(i==19) tunfresultcor->GetYaxis()->SetRangeUser(theoryxsec->GetMinimum()*0.8,theoryxsec->GetMaximum()*1.3 );
      */

      //if(i<20)       theoryxsec->GetYaxis()->SetRangeUser(theoryxsec->GetMinimum()*0.7,theoryxsec->GetMaximum()*1.3 );
      //else if(i==19) theoryxsec->GetYaxis()->SetRangeUser(theoryxsec->GetMinimum()*0.7,theoryxsec->GetMaximum()*1.3 );
      //else if(i>=20) theoryxsec->GetYaxis()->SetRangeUser(theoryxsec->GetMinimum()*0.7,theoryxsec->GetMaximum()*1.3 );
      
      //else if( theoryxsec->GetMaximum()*0.6 > theoryxsec->GetMinimum() ) theoryxsec->GetYaxis()->SetRangeUser(0,theoryxsec->GetMaximum()*1.45 );
      //else theoryxsec->GetYaxis()->SetRangeUser(theoryxsec->GetMinimum()*0.81,theoryxsec->GetMaximum()+0.19*theoryxsec->GetMinimum() );
      
      // Amandeep : Modifying for 2D plots
      // Original :
      // else if(inorm) {
      //   if( MatrixUnf::isCtype(fVariableNames[i]) )        theoryxsec->GetYaxis()->SetRangeUser(0.,1.57*1.05);
      //   else if( MatrixUnf::isCPMtype(fVariableNames[i]) ) theoryxsec->GetYaxis()->SetRangeUser(0.,1.04*1.05);
      //   else theoryxsec->GetYaxis()->SetRangeUser(0.46001,0.576);
      // }

      // Modified : Adjust based on Min and Max
      theoryxsec->GetYaxis()->SetRangeUser(theoryxsec->GetMinimum()*0.7,theoryxsec->GetMaximum()*1.3 );

      // Jason 
      if(fVariableNames[i].Contains("_mttbar")) {
	if(inorm) {
	  gPad->SetLogy();
	  theoryxsec->GetYaxis()->SetRangeUser(.000001,1.);
	}
	else{
	  gPad->SetLogy();
	  theoryxsec->GetYaxis()->SetRangeUser(.001,1000.);
	}
      }
      else {
	if(inorm) {
	  gPad->SetLogy();
	  theoryxsec->GetYaxis()->SetRangeUser(.01,1000.);
	}
	else{
	  gPad->SetLogy();
	  theoryxsec->GetYaxis()->SetRangeUser(1.,1000000.);
	}
      }
      // End
      
      
      //if(inorm) {
      //  if( MatrixUnf::isCtype(fVariableNames[i]) )        theoryxsec->GetYaxis()->SetRangeUser(theoryxsec->GetMinimum()*0.6,theoryxsec->GetMaximum()*1.3 );
      //  else if( MatrixUnf::isCPMtype(fVariableNames[i]) ) theoryxsec->GetYaxis()->SetRangeUser(theoryxsec->GetMinimum()*0.6,theoryxsec->GetMaximum()*1.3 );
      //  else theoryxsec->GetYaxis()->SetRangeUser(theoryxsec->GetMinimum()*0.7,theoryxsec->GetMaximum()*1.2 );
      //}
      // End
      
      //else if(i<19) {
      // if( theoryxsec->GetMaximum()*0.6 > theoryxsec->GetMinimum() ) 
      //     theoryxsec->GetYaxis()->SetRangeUser(0,theoryxsec->GetMaximum()*1.45 );
      //  else 
      //      theoryxsec->GetYaxis()->SetRangeUser(theoryxsec->GetMinimum()*0.81,theoryxsec->GetMaximum()+0.19*theoryxsec->GetMinimum() );
      //}

      // Amandeep :  For sure this is wrong for 2D since it tries to get bin width
      // Main question is how to define the stat and totalcov matrix for 2D
      // Does it still remain the same ?

      // Original :
      // if (divBinWidth) {
      //     for (int j = 0; j < nbinsrhoi; j++)
      //         for (int k = 0; k < nbinsrhoi; k++) {
      //             StatCovMats.at(inorm)(j, k)  /= tunfresultcor->GetBinWidth(j + 1) * tunfresultcor->GetBinWidth(k + 1);
      //             TotalCovMats.at(inorm)(j, k) /= tunfresultcor->GetBinWidth(j + 1) * tunfresultcor->GetBinWidth(k + 1);
      //         }
      // }

      // Modified :
      if (divBinWidth) {
          for (int j = 0; j < nbinsrhoi; j++)
              for (int k = 0; k < nbinsrhoi; k++) {
                  StatCovMats.at(inorm)(j, k)  /= GetPhysicalBinWidth(fVariableNames[i], j + 1) * GetPhysicalBinWidth(fVariableNames[i], k + 1);
                  TotalCovMats.at(inorm)(j, k) /= GetPhysicalBinWidth(fVariableNames[i], j + 1) * GetPhysicalBinWidth(fVariableNames[i], k + 1);
              }
      }
      // End


      // Plot statistical and total error of full result
      Double_t mexl[nbinsrhoi];
      Double_t mexh[nbinsrhoi];
      Double_t binWidths[nbinsrhoi];
      Double_t binCenters[nbinsrhoi];
      Double_t binValues[nbinsrhoi];
      Double_t binStatErrors[nbinsrhoi];
      Double_t binTotalErrors[nbinsrhoi];

      for (int j=0; j<nbinsrhoi;j++){
        mexl[j] = 0;
        mexh[j] = 0;
        binWidths[j]      = tunfresultcor->GetBinWidth(j+1);
        binCenters[j]     = tunfresultcor->GetBinCenter(j+1);
        binValues[j]      = tunfresultcor->GetBinContent(j+1);
        binStatErrors[j]  = sqrt(StatCovMats.at(inorm)(j,j));
        binTotalErrors[j] = sqrt(TotalCovMats.at(inorm)(j,j));

        std::cout << std::endl;
        std::cout << "binWidths for bin "  << j << " " << binWidths[j]  << std::endl;
        std::cout << "binValues for bin "  << j << " " << binValues[j]  << std::endl;
        std::cout << "binStatErrors for bin "  << j << " " << binStatErrors[j]  << std::endl;
        std::cout << std::endl;
      }

      TGraphAsymmErrors *tga_DiffXSecPlot        = new TGraphAsymmErrors(nbinsrhoi, binCenters, binValues, mexl, mexh, binStatErrors, binStatErrors);
      tga_DiffXSecPlot->SetMarkerStyle(1);
      tga_DiffXSecPlot->SetMarkerColor(kBlack);
      tga_DiffXSecPlot->SetMarkerSize(0.8);
      tga_DiffXSecPlot->SetLineWidth(2);
      tga_DiffXSecPlot->SetLineColor(kBlack);
     
      TGraphAsymmErrors *tga_DiffXSecPlotwithSys = new TGraphAsymmErrors(nbinsrhoi, binCenters, binValues, mexl, mexh, binTotalErrors, binTotalErrors);
      tga_DiffXSecPlotwithSys->SetMarkerStyle(20);
      tga_DiffXSecPlotwithSys->SetMarkerColor(kBlack);
      tga_DiffXSecPlotwithSys->SetMarkerSize(0.8);
      tga_DiffXSecPlotwithSys->SetLineWidth(2);
      tga_DiffXSecPlotwithSys->SetLineColor(kBlack);


      tunfresultcor->SetMarkerStyle(tga_DiffXSecPlotwithSys->GetMarkerStyle());
      tunfresultcor->SetMarkerSize(tga_DiffXSecPlotwithSys->GetMarkerSize());
      tunfresultcor->SetMarkerColor(tga_DiffXSecPlotwithSys->GetMarkerColor());


      // Plot statistical and total error of full result
      TGraphAsymmErrors *ratio_stat = 0, *ratio_tota = 0;
      if(tga_DiffXSecPlot) ratio_stat        = (TGraphAsymmErrors*)tga_DiffXSecPlot->Clone("ratio_stat");
      if(tga_DiffXSecPlotwithSys) ratio_tota = (TGraphAsymmErrors*)tga_DiffXSecPlotwithSys->Clone("ratio_tota");

      if(ratio_stat){
          ratio_stat->SetFillStyle(1001);
          ratio_stat->SetFillColor(kGray+1);
          ratio_stat->SetLineColor(0);

          std::cout << "binStatErrors " << std::endl;
          for (Int_t iter=0; iter < tga_DiffXSecPlot->GetN(); iter++)
          {
              // Amandeep : This is also wrong then 
              // Orignal  :
              // double binWidth = tunfresultcor->GetBinWidth(iter + 1); // 1
              // Modified : 
              // For 1D : pow(4,1), and for 2D pow(2,2)
              double binWidth    = binWidths[iter];
              // double binWidth = pow(rebinfine, numDimensions);
              // FIXME :  pow(rebinfine, numDimensions) doesn't seem to be working 
              // End
              double x                 = tga_DiffXSecPlot->GetX()[iter];
              double y_ratio_stat      = tga_DiffXSecPlot->GetY()[iter] / tga_DiffXSecPlot->GetY()[iter];
              double abserr_ratio_stat = y_ratio_stat - std::abs(tga_DiffXSecPlot->GetErrorY(iter) - tga_DiffXSecPlot->GetY()[iter]) / tga_DiffXSecPlot->GetY()[iter];
              ratio_stat->SetPoint(iter, x, y_ratio_stat);
              ratio_stat->SetPointError(iter, binWidth/2, binWidth/2, abserr_ratio_stat, abserr_ratio_stat);
              std::cout << abserr_ratio_stat  << std::endl;
          }
      }

      if(ratio_tota){
          ratio_tota->SetFillStyle(1001);
          ratio_tota->SetFillColor(kOrange-4);
          ratio_tota->SetLineColor(0);
          
          std::cout << "binTotalErrors " << std::endl;
          for (Int_t iter = 0; iter<tga_DiffXSecPlotwithSys->GetN(); iter++)
          {
              // Amandeep  
              // Orignal  :
              // double binWidth = tunfresultcor->GetBinWidth(iter + 1);
              // Modified :
              // For 1D : pow(4,1), and for 2D pow(2,2)
              double binWidth    = binWidths[iter];
              // double binWidth = pow(rebinfine, numDimensions);
              // FIXME :  pow(rebinfine, numDimensions) doesn't seem to be working 
              // End
              double x                 = tga_DiffXSecPlotwithSys->GetX()[iter];
              double y_ratio_tota      = tga_DiffXSecPlotwithSys->GetY()[iter] / tga_DiffXSecPlotwithSys->GetY()[iter];
              double abserr_ratio_tota = y_ratio_tota - std::abs(tga_DiffXSecPlotwithSys->GetErrorY(iter) - tga_DiffXSecPlotwithSys->GetY()[iter]) / tga_DiffXSecPlotwithSys->GetY()[iter];
              ratio_tota->SetPoint(iter, x, y_ratio_tota);
              ratio_tota->SetPointError(iter, binWidth/2, binWidth/2, abserr_ratio_tota, abserr_ratio_tota);
              std::cout << abserr_ratio_tota << std::endl;
          }
      }


      bool doPlotAllMCs = kFALSE;
      gStyle->SetErrorX(0.5);
      setTheoryStyleAndFillLegend(tunfresultcor, "data");
      l2->AddEntry(tga_DiffXSecPlotwithSys,"Unfolded data", "pE");
      setTheoryStyleAndFillLegend(theoryxsec, "nominal", l2);

      // Increase Y TitleOffset for Bs
      //if( MatrixUnf::isBtype(fVariableNames[i]) ) {
      //  theoryxsec->GetYaxis()->SetTitleOffset(1.77);
      //  gStyle->SetPadLeftMargin(0.19);
      //  gStyle->SetPadRightMargin(0.04);
      //}

      theoryxsec->GetYaxis()->SetLabelSize(0.04);
      theoryxsec->GetYaxis()->SetTitleSize(0.035);
      theoryxsec->GetYaxis()->SetTitleOffset(2.4);

      // if(doPlotAllMCs ) setTheoryStyleAndFillLegend(POWHEGV2HERWIGxsec, "powhegv2herwig", l2);
      // if(doPlotAllMCs || i>=0) setTheoryStyleAndFillLegend(AMCATNLOxsec, "amcatnlo", l2);
      // if(doPlotAllMCs ) setTheoryStyleAndFillLegend(TOP_PTxsec, "toppt", l2);
      theoryxsec->Draw();
      // if(doPlotAllMCs ) POWHEGV2HERWIGxsec->Draw("same");
      // if(doPlotAllMCs || i>=0) AMCATNLOxsec->Draw("same");
      // if(doPlotAllMCs ) TOP_PTxsec->Draw("same");


      TH1D* NLOWxsec;
      TH1D* NLOWxsecuncorr;
      TH1D* NNLOxsec;

      // Amandeep : Modifying this 
      // if(i<20) {
      if(i<39) {
        NLOWxsec       = (TH1D*) theoryxsec->Clone("NLOWxsec");
        NLOWxsec->Reset();
        NLOWxsecuncorr = (TH1D*) theoryxsec->Clone("NLOWxsecuncorr");
        NLOWxsecuncorr->Reset();

        for (int i_bin = 0; i_bin < nbinsrhoi; ++i_bin)
        {
          NLOWxsec->SetBinContent(i_bin+1, tempNLOWxsec->GetBinContent(i_bin+1));
          NLOWxsecuncorr->SetBinContent(i_bin+1, tempNLOWxsecuncorr->GetBinContent(i_bin+1));
          std::cout << "tempNLOWxsec->GetBinContent(i_bin+1) " << i_bin+1 <<  tempNLOWxsec->GetBinContent(i_bin+1) << std::endl;
          std::cout << "tempNLOWxsecuncorr->GetBinContent(i_bin+1) " << i_bin+1 << tempNLOWxsecuncorr->GetBinContent(i_bin+1) << std::endl;
          std::cout << std::endl;
        }
      }

      else {
        NLOWxsec       = (TH1D*) NLOW_histo[0]->Clone("NLOWxsec");
        NLOWxsecuncorr = (TH1D*) NLOW_histo[3]->Clone("NLOWxsecuncorr");
      }

      NLOWxsec->Scale(theoryxsectemp->Integral());
      NLOWxsecuncorr->Scale(theoryxsectemp->Integral());

      for (int i_bin = 0; i_bin < nbinsrhoi; ++i_bin)
      {
        // so the histogram is not drawn with vertical connecting lines
        NLOWxsec->SetBinError(i_bin+1, 1e-10);
        NLOWxsecuncorr->SetBinError(i_bin+1, 1e-10);
      }

      if(inorm) NLOWxsec->Scale(1./NLOWxsec->Integral());
      if(inorm) NLOWxsecuncorr->Scale(1./NLOWxsecuncorr->Integral());
      if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(NLOWxsec);
      if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(NLOWxsecuncorr);


      if (fVariableNames[i] != "llbar_delta_phi") {
        // Amandeep : Disabling for now
        setTheoryStyleAndFillLegend(NLOWxsec, "NLOW", l2N);
        if(i<10) setTheoryStyleAndFillLegend(NLOWxsecuncorr, "unpol", l2N);
        else setTheoryStyleAndFillLegend(NLOWxsecuncorr, "uncorr", l2N);
        NLOWxsecuncorr->Draw("SAME");
        NLOWxsec->Draw("SAME");
        setResultLegendStyleN(l2N, fVariableNames[i]);
        // End

      }
      else{

        NNLOxsec = (TH1D*) NNLO_histo[0]->Clone("NNLOxsec");
        NNLOxsec->Scale(theoryxsectemp->Integral());
        if(inorm) NNLOxsec->Scale(1./NNLOxsec->Integral());
        if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(NNLOxsec);
        // setTheoryStyleAndFillLegend(NLOWxsec, "NLOW", l2);
        // setTheoryStyleAndFillLegend(NLOWxsecuncorr, "uncorr", l2);
        // setTheoryStyleAndFillLegend(NNLOxsec, "NNLO", l2);

        // Amandeep : Disabling for now
        // NLOWxsecuncorr->Draw("SAME");
        // NLOWxsec->Draw("SAME");
        // NNLOxsec->Draw("SAME");
        // End

      }

      //tunfresultcor->Draw("same");
      gStyle->SetEndErrorSize(10);
      tga_DiffXSecPlot->Draw("p, SAME");
      tga_DiffXSecPlotwithSys->Draw("p, SAME, Z");
      //gPad->RedrawAxis();

      /*
      TLegend *l2=new TLegend(0.7,(i<20?0.7:0.75),0.9,0.9);
      l2->AddEntry(tunfresultcor,"Unfolded Data", "p");
      l2->AddEntry(theoryxsec, "POWHEG+Pythia8");
      if(i<20) l2->AddEntry(NLOWxsec, "NLO");
      if(i<10) l2->AddEntry(NLOWxsecuncorr, "Unpolarized");
      else if(i<20) l2->AddEntry(NLOWxsecuncorr, "No spin correlation");
      //l2->SetFillColor(kWhite);
      //l2->SetLineColor(kBlack);
      l2->SetTextFont(42);
      l2->SetTextAlign(12);
      l2->SetTextSize(0.035);
      l2->SetFillStyle(0);
      l2->SetBorderSize(0);
      */

      setResultLegendStyle(l2, fVariableNames[i]);
      l2->Draw("same");
      // if (fVariableNames[i] != "llbar_delta_phi") l2N->Draw("same");

      Draw2DLabels(fVariableNames[i]);
      DrawCMSLabels(0.045, fLumi);
      if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);


      // Amandeep : Modifying for 2D
      // Original :
      double ratioL = 0.81;
      double ratioU = 1.19;

      if(inorm) {
        // Consider different kinds of variables and num dimensions
        // 1D variables
        if ( MatrixUnf::isCtype(fVariableNames[i]) && numDimensions == 1) {
          ratioL = 0.8501;
          ratioU = 1.1499;
        }
        else if( MatrixUnf::isBtype(fVariableNames[i]) && numDimensions == 1) {
          ratioL = 0.9251;
          ratioU = 1.0749;
        }
        else if( MatrixUnf::isCPMtype(fVariableNames[i]) && numDimensions == 1) {
          ratioL = 0.91;
          ratioU = 1.09;
        }
        else if( MatrixUnf::isDtype(fVariableNames[i]) && numDimensions == 1) {
          ratioL = 0.9251;
          ratioU = 1.0749;
        }
        else if( MatrixUnf::isOtherSymType(fVariableNames[i]) && numDimensions == 1 ) {
          ratioL = 0.9251;
          ratioU = 1.0749;
        }

        // 2D variables
        else if( MatrixUnf::isCtype(fVariableNames[i]) && numDimensions == 2) {
          ratioL = 0.66;
          ratioU = 1.34;
        }
        else if( MatrixUnf::isBtype(fVariableNames[i]) && numDimensions == 2) {
          ratioL = 0.71;
          ratioU = 1.29;
        }
        else if( MatrixUnf::isCPMtype(fVariableNames[i]) && numDimensions == 2) {
          ratioL = 0.71;
          ratioU = 1.29;
        }
        // FIXME : Add others
      }

      else {
        
        // Consider different kinds of variables and num dimensions
        // 1D variables
        if( MatrixUnf::isCtype(fVariableNames[i]) && numDimensions == 1) {
          ratioL = 0.81;
          ratioU = 1.19;
        }
        else if( MatrixUnf::isBtype(fVariableNames[i]) && numDimensions == 1) {
          ratioL = 0.8501;
          ratioU = 1.1499;
        }
        else if( MatrixUnf::isCPMtype(fVariableNames[i]) && numDimensions == 1) {
          ratioL = 0.8501;
          ratioU = 1.1499;
        }
        else if( MatrixUnf::isDtype(fVariableNames[i]) && numDimensions == 1) {
          ratioL = 0.8501;
          ratioU = 1.1499;
        }
        else if( MatrixUnf::isOtherSymType(fVariableNames[i]) && numDimensions == 1 ) {
          ratioL = 0.8501;
          ratioU = 1.1499;
        }

        // 2D variables
        else if( MatrixUnf::isCtype(fVariableNames[i]) && numDimensions == 2) {
          ratioL = 0.71;
          ratioU = 1.29;
        }
        else if( MatrixUnf::isBtype(fVariableNames[i]) && numDimensions == 2) {
          ratioL = 0.71;
          ratioU = 1.29;
        }
        else if( MatrixUnf::isCPMtype(fVariableNames[i]) && numDimensions == 2) {
          ratioL = 0.71;
          ratioU = 1.29;
        }
        // FIXME : Add others
      }

      // Modified :
      // double ratioL = 0.51;
      // double ratioU = 1.49;
      // End

      // if(doPlotAllMCs || i>=0) common::drawRatioXSEC(tunfresultcor, theoryxsec, ratio_stat, ratio_tota, POWHEGV2HERWIGxsec, AMCATNLOxsec, TOP_PTxsec, NLOWxsecuncorr, NLOWxsec, nullptr, ratioL, ratioU);
      
      // Original
      // if (fVariableNames[i] == "llbar_delta_phi") common::drawRatioXSEC(tunfresultcor, theoryxsec, ratio_stat, ratio_tota, AMCATNLOxsec, NLOWxsecuncorr, NLOWxsec, NNLOxsec,nullptr,nullptr, ratioL, ratioU);
      // else if (fVariableNames[i] == "c_Prk" || fVariableNames[i] == "ll_cHel" || fVariableNames[i] == "ll_cLab") common::drawRatioXSEC(tunfresultcor, theoryxsec, ratio_stat, ratio_tota, AMCATNLOxsec, NLOWxsecuncorr, NLOWxsec, nullptr,nullptr,nullptr, ratioL, ratioU, kTRUE);
      // else common::drawRatioXSEC(tunfresultcor, theoryxsec, ratio_stat, ratio_tota, AMCATNLOxsec, NLOWxsecuncorr, NLOWxsec, nullptr,nullptr,nullptr, ratioL, ratioU);
      
      // Modified : Replace AMCATNLOxsec with nullptr
      // if      (fVariableNames[i] == "llbar_delta_phi") common::drawRatioXSEC(tunfresultcor, theoryxsec, ratio_stat, ratio_tota, nullptr, NLOWxsecuncorr, NLOWxsec, NNLOxsec,nullptr,nullptr, ratioL, ratioU);
      // else if (fVariableNames[i] == "c_Prk" || fVariableNames[i] == "ll_cHel" || fVariableNames[i] == "ll_cLab") common::drawRatioXSEC(tunfresultcor, theoryxsec, ratio_stat, ratio_tota, nullptr, NLOWxsecuncorr, NLOWxsec, nullptr,nullptr,nullptr, ratioL, ratioU, kTRUE);
      // else     common::drawRatioXSEC(tunfresultcor, theoryxsec, ratio_stat, ratio_tota, nullptr, nullptr, nullptr, nullptr,nullptr,nullptr, ratioL, ratioU);
      
      // Amandeep : ll_cHel breaks probably due to the custom TStyle, use this for everything (for now)
      common::drawRatioXSEC(tunfresultcor, theoryxsec, ratio_stat, ratio_tota, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, ratioL, ratioU, kFALSE, fVariableNames[i]);
      // End
      
      //else if(doPlotAllMCs || i>=0) common::drawRatioXSEC(tunfresultcor, theoryxsec, ratio_stat, ratio_tota, AMCATNLOxsec, NLOWxsecuncorr, NLOWxsec, nullptr,nullptr,nullptr, ratioL, ratioU);
      //else common::drawRatioXSEC(tunfresultcor, theoryxsec, ratio_stat, ratio_tota, NLOWxsecuncorr, NLOWxsec, nullptr,nullptr,nullptr, nullptr, ratioL, ratioU);
      
      c2->Update();
      c2->Modified();
      c2->Update();

      c2->SaveAs(Form("%s/UnfoldedResults%s_%s.pdf" ,PLOTFILE.Data(),inorm?"Norm":"",fVariableNames[i].Data()));
      c2->SaveAs(Form("%s/UnfoldedResults%s_%s.root",PLOTFILE.Data(),inorm?"Norm":"",fVariableNames[i].Data()));
      c2->SaveAs(Form("%s/UnfoldedResults%s_%s.C"   ,PLOTFILE.Data(),inorm?"Norm":"",fVariableNames[i].Data()));

      //revert change made for Bs
      if( MatrixUnf::isBtype(fVariableNames[i]) ) {
        gStyle->SetPadLeftMargin(0.18);
        gStyle->SetPadRightMargin(0.05);
      }


      //dPhi theory plot for talk
      if (i>=0) {
        bool drawscalevars = kTRUE;
        if(i<39 || !inorm) drawscalevars = kFALSE;
        TCanvas *c3 = new TCanvas("c3", "c3");
        TLegend *l3 = new TLegend();
        TLegend *l4 = new TLegend();

        // TH1D* PowhegScaleVarHists[ScaleEnvelopeList.size()-4];

        if(inorm) {
          setTheoryStyleAndFillLegend(theoryxsec, "nominal", l4);
        //   for (unsigned int i_env = 0; i_env < ScaleEnvelopeList.size()-4; ++i_env)
        //   {
        //     PowhegScaleVarHists[i_env] = (TH1D*) ScaleEnvelopeHists[i_env]->Clone();
        //     if(inorm) PowhegScaleVarHists[i_env]->Scale(1./PowhegScaleVarHists[i_env]->Integral());
        //     setTheoryStyleAndFillLegend(PowhegScaleVarHists[i_env], "POWHEGvars", l4);
        //     l4->AddEntry(PowhegScaleVarHists[i_env], ScaleEnvelopeList[i_env].Data(), "l");
        //   }
        }

        setTheoryStyleAndFillLegend(NLOWxsec, "NLOW", l3);
        setTheoryStyleAndFillLegend(NLOWxsecuncorr, "uncorr", l3);

        if(drawscalevars) {
          if (fVariableNames[i] == "llbar_delta_phi") {
            setTheoryStyleAndFillLegend(NNLOxsec, "NNLO", l3);
            setTheoryStyleAndFillLegend(NNLO_histo[1],"NNLO 2,2", nullptr);
            l3->AddEntry(NNLO_histo[1],"NNLO 2,2", "l");
            setTheoryStyleAndFillLegend(NNLO_histo[2],"NNLO .5,.5", nullptr);
            l3->AddEntry(NNLO_histo[2],"NNLO .5,.5", "l");
            setTheoryStyleAndFillLegend(NNLO_histo[3],"NNLO 1,.5", nullptr);
            l3->AddEntry(NNLO_histo[3],"NNLO 1,.5", "l");
            setTheoryStyleAndFillLegend(NNLO_histo[4],"NNLO .5,1", nullptr);
            l3->AddEntry(NNLO_histo[4],"NNLO .5,1", "l");
            setTheoryStyleAndFillLegend(NNLO_histo[5],"NNLO 2,1", nullptr);
            l3->AddEntry(NNLO_histo[5],"NNLO 2,1", "l");
            setTheoryStyleAndFillLegend(NNLO_histo[6],"NNLO 1,2", nullptr);
            l3->AddEntry(NNLO_histo[6],"NNLO 1,2", "l");
          }
          setTheoryStyleAndFillLegend(NLOW_histo[1],"NLOWscaleup", l3);
          setTheoryStyleAndFillLegend(NLOW_histo[2],"NLOWscaledown", l3);
          setTheoryStyleAndFillLegend(NLOW_histo[4],"uncorrscaleup", l3);
          setTheoryStyleAndFillLegend(NLOW_histo[5],"uncorrscaledown", l3);
        }

        NLOWxsecuncorr->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
        //NLOWxsecuncorr->GetYaxis()->SetTitle(Form("d#sigma/d%s [pb]",quantity.Data()));
        NLOWxsecuncorr->GetYaxis()->SetTitle(Form("#frac{d#sigma}{d%s} [pb]",quantity.Data()));
        if(inorm) NLOWxsecuncorr->GetYaxis()->SetTitle(Form("#frac{1}{#sigma} #frac{d#sigma}{d%s}",quantity.Data()));

        NLOWxsec->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
        //NLOWxsec->GetYaxis()->SetTitle(Form("d#sigma/d%s [pb]",quantity.Data()));
        NLOWxsec->GetYaxis()->SetTitle(Form("#frac{d#sigma}{d%s} [pb]",quantity.Data()));
        if(inorm) NLOWxsec->GetYaxis()->SetTitle(Form("#frac{1}{#sigma} #frac{d#sigma}{d%s}",quantity.Data()));


        NLOWxsec->SetLineStyle(1);
        NLOWxsec->SetLineWidth(2);
        NLOWxsecuncorr->SetLineStyle(1);
        NLOWxsecuncorr->SetLineWidth(2);


        for (int i_bin = 0; i_bin < nbinsrhoi; ++i_bin)
        {
          //so the histogram is  drawn with vertical connecting lines
          NLOWxsec->SetBinError(i_bin+1, 0);
          NLOWxsecuncorr->SetBinError(i_bin+1, 0);

          if(inorm) {

            // for (unsigned int i_env = 0; i_env < ScaleEnvelopeList.size()-4; ++i_env)
            // {
            //   PowhegScaleVarHists[i_env]->SetBinError(i_bin+1, 0);
            // }

          }
        }

        if (fVariableNames[i] == "llbar_delta_phi") {
          NLOWxsecuncorr->Draw();
          NLOWxsec->Draw("SAME");
        }
        else {
          NLOWxsec->Draw();
          NLOWxsecuncorr->Draw("SAME");
        }

         
        if(inorm) {
          theoryxsec->Draw("SAME");
          // for (unsigned int i_env = 0; i_env < ScaleEnvelopeList.size()-4; ++i_env)
          // {
          //   if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(PowhegScaleVarHists[i_env]);
          //   PowhegScaleVarHists[i_env]->Draw("SAME");
          // }
        }

        if(drawscalevars) {
          if (fVariableNames[i] == "llbar_delta_phi") {
            NNLOxsec->Draw("SAME");

            if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(NNLO_histo[1]);
            if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(NNLO_histo[2]);
            if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(NNLO_histo[3]);
            if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(NNLO_histo[4]);
            if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(NNLO_histo[5]);
            if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(NNLO_histo[6]);

            NNLO_histo[1]->Draw("SAME");
            NNLO_histo[2]->Draw("SAME");
            NNLO_histo[3]->Draw("SAME");
            NNLO_histo[4]->Draw("SAME");
            NNLO_histo[5]->Draw("SAME");
            NNLO_histo[6]->Draw("SAME");
          }

          if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(NLOW_histo[1]);
          if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(NLOW_histo[2]);
          if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(NLOW_histo[4]);
          if(divBinWidth) MatrixUnf::GetBinWidthCorrectedCrossSec(NLOW_histo[5]);

          NLOW_histo[1]->Draw("SAME");
          NLOW_histo[2]->Draw("SAME");
          NLOW_histo[4]->Draw("SAME");
          NLOW_histo[5]->Draw("SAME");
        }

        if(drawscalevars) setResultLegendStyle(l3, "large");
        else setResultLegendStyle(l3, "");
        setResultLegendStyleN(l4, "large");

        if(drawscalevars) l3->SetTextSize(0.02);
        else if(inorm) l3->SetTextSize(0.03);
        else l3->SetTextSize(0.045);
        l4->SetTextSize(0.025);
        l3->Draw("same");
        if(inorm) l4->Draw("same");

        // Plot statistical and total error of full result
        //TGraphAsymmErrors *ratio_stat = 0, *ratio_tota = 0;

        ratioL = 0.81;
        ratioU = 1.19;

        // Amandeep : Commenting for now
        // if(drawscalevars && (fVariableNames[i] == "llbar_delta_phi")) common::drawRatioXSEC(NLOWxsecuncorr, NLOWxsec, nullptr, nullptr, NLOW_histo[1], NLOW_histo[2], NLOW_histo[4], NLOW_histo[5],NNLO_histo[1],NNLO_histo[2], ratioL, ratioU, kFALSE, fVariableNames[i]);
        // else if(drawscalevars) common::drawRatioXSEC(NLOWxsecuncorr, NLOWxsec, nullptr, nullptr, NLOW_histo[1], NLOW_histo[2], NLOW_histo[4], NLOW_histo[5],nullptr,nullptr, ratioL, ratioU, kFALSE, fVariableNames[i]);
        // else if(inorm) common::drawRatioXSEC(NLOWxsecuncorr, NLOWxsec, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,nullptr, nullptr, ratioL, ratioU, kFALSE, fVariableNames[i]);
        // else common::drawRatioXSEC(NLOWxsecuncorr, NLOWxsec, nullptr, nullptr, nullptr,nullptr,nullptr, nullptr,nullptr,nullptr, ratioL, ratioU, kFALSE, fVariableNames[i]);
        common::drawRatioXSEC(NLOWxsecuncorr, NLOWxsec, nullptr, nullptr, nullptr,nullptr,nullptr, nullptr,nullptr,nullptr, ratioL, ratioU, kFALSE, fVariableNames[i]);


        c3->Update();
        c3->Modified();
        c3->Update();

        c3->SaveAs(Form("%s/%sTheory%s.pdf",PLOTFILE.Data(),fVariableNames[i].Data(),inorm?"Norm":""));
        c3->SaveAs(Form("%s/%sTheory%s.root",PLOTFILE.Data(),fVariableNames[i].Data(),inorm?"Norm":""));
        c3->SaveAs(Form("%s/%sTheory%s.C",PLOTFILE.Data(),fVariableNames[i].Data(),inorm?"Norm":""));

        delete c3;
        delete l3;
        delete l4;

      }

      delete theoryxsec;
      // delete AMCATNLOxsec;
      // delete POWHEGV2HERWIGxsec;
      // delete TOP_PTxsec;
      delete tunfresultcor;
      delete NLOWxsec;
      delete NLOWxsecuncorr;

      delete c2;
      delete l2;
      // delete l2N;


    }//inorm


    unffile->Close();
    delete unffile;

    
  }

  if(fsetSyst) {

    printf("\n Writing systematics output file ...\n");
    systoutfile->Write();
    //systoutfile->Close();
    delete systoutfile;


    for (int inorm = 0; inorm <= 1; ++inorm) {

      //chi2 table
      for (int i_unctype = 0; i_unctype <= 1; ++i_unctype)
      {

        printf(" Table of %s $\\chi^2$ values (5 ndf in all cases) (%s): \n",(i_unctype ? " Stat" : " Total" ),(inorm ? "Normalised" : "Absolute" ));

        printf("\\begin{tabular}{l | c c c c c c c c } \n");
        printf("\\hline \n");
        printf("\\textbf{Distribution} & \\small{\\Powhegvtwo P8} & \\small{\\MGaMCatNLO} & \\small{\\Powhegvtwo H++} & \\small{Top \\pt rw} & \\small{NLOW} & \\small{No SC/Pol.} & \\small{measA} & \\small{measN} \\\\ \n");
        for (unsigned int iVar = 0; iVar < fVariableNames.size(); ++iVar)
        {
          printf("$%s$ & $%2.1f$ & $%2.1f$ & $%2.1f$ & $%2.1f$ & $%2.1f$ & $%2.1f$ & $%2.1f$ & $%2.1f$  \\\\ \n",(MatrixUnf::CoefName(fVariableNames[iVar],kTRUE)).Data(), Chi2values[inorm][i_unctype][iVar][0], Chi2values[inorm][i_unctype][iVar][1], Chi2values[inorm][i_unctype][iVar][2], Chi2values[inorm][i_unctype][iVar][3], Chi2values[inorm][i_unctype][iVar][4], Chi2values[inorm][i_unctype][iVar][5], Chi2values[inorm][i_unctype][iVar][6], Chi2values[inorm][i_unctype][iVar][7]);
        }
        printf("\\hline \n");
        printf("\\end{tabular} \n \n");

      }

    }


    printf(" Table of $f_{SM}$ results excluding theory uncertainty");

    printf("\\begin{tabular}{l | c } \n");
    printf("\\hline \n");
    printf("\\textbf{Coefficient} & \\textbf{$f_{SM}$}  \\\\ \n");
    for (unsigned int iVar = 0; iVar < fSMTableNames.size(); ++iVar)
    {
      printf("$%s$ & $%2.5f \\pm %2.5f \\pm %2.5f (\\pm %2.5f)$  \\\\ \n",fSMTableNames.at(iVar).Data(), fSMTable.at(iVar), fSMTableStatErr.at(iVar), fSMTableSystErr.at(iVar), fSMTableErr.at(iVar));
    }
    printf("\\hline \n");
    printf("\\end{tabular} \n");



    std::vector<std::string> SystNames = tempVectorOfInputSystematics;
    SystNames.push_back("Total Syst.");
    SystNames.push_back("Stat.");
    SystNames.push_back("MC Stat.");
    SystNames.push_back("Bkg Stat.");
    SystNames.push_back("Total Stat.");
    SystNames.push_back("Total");


    std::vector<TString>CoefNames;
    std::vector<TString>CoefNames2;
    for (unsigned int iVar = 0; iVar < fVariableNames.size(); ++iVar)
    {
        for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
	  const TVectorD *vec_rebinnedBinEdges = fVariableNames[iVar].Contains("_mttbar") ? generatorBinningsRebinned[iVar]->GetDistributionBinning(1) : generatorBinningsRebinned[iVar]->GetDistributionBinning(0);
          const double *rebinnedBinEdges = vec_rebinnedBinEdges->GetMatrixArray();
          double lowerBoundary = rebinnedBinEdges[ithCoefficient];
          double upperBoundary = rebinnedBinEdges[ithCoefficient + 1];
          CoefNames.push_back(MatrixUnf::CoefName(fVariableNames[iVar], kTRUE, lowerBoundary, upperBoundary));
          CoefNames2.push_back(MatrixUnf::CoefName(fVariableNames[iVar], kFALSE, lowerBoundary, upperBoundary));
        }
    }


    // Print tables of results and systematics for the coefficients extracted from the rebinnedA normalised differential xsecs
    int nBvar = 0;
    int nCvar = 0;
    for (int ithVariableName = 0; ithVariableName < fVariableNames.size(); ithVariableName++) {
      if (MatrixUnf::isCtype(fVariableNames[ithVariableName].Data()) 
      || MatrixUnf::isCPMtype(fVariableNames[ithVariableName].Data())
      || MatrixUnf::isDtype(fVariableNames[ithVariableName].Data())
      || MatrixUnf::isOtherSymType(fVariableNames[ithVariableName].Data())) {
        nCvar += 1;
      } else if(MatrixUnf::isBtype(fVariableNames[ithVariableName].Data()) 
      || MatrixUnf::isBPMtype(fVariableNames[ithVariableName].Data())) {
        nBvar += 1;
      }
    }

    TString NameTag = "_rebinnedA";

    for (int i_tag = 0; i_tag < 2; ++i_tag)
    {
      if(i_tag==1) NameTag = "_rebinnedB";

      for (int inorm = 0; inorm <= 1; ++inorm)
      {

        printf(" Table of results for Coef%s%s: \n",NameTag.Data(),(inorm?"Norm":""));

        printf("\\begin{tabular}{l | c c c} \n");
        printf("\\hline \n");
        printf("\\textbf{Coefficient} & \\textbf{Measured} & \\textbf{\\Powhegvtwo} & \\textbf{\\MGaMCatNLO} & \\textbf{NLO calculation}  \\\\ \n");
        for (unsigned int iVar = 0; iVar < fVariableNames.size(); ++iVar)
        {
          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            printf("$%s$ & $%2.3f \\pm %2.3f$ & $%2.3f^{+%2.3f}_{-%2.3f}$ & $%2.3f \\pm %2.3f$ &  \\\\ \n", CoefNames.at(iVar * numCoefficients + ithCoefficient).Data(), coefTable[i_tag * 2 + inorm][iVar][ithCoefficient], systTable[i_tag * 2 + inorm][iVar][ithCoefficient][SystNames.size() - 1], MCcoefTable[iVar][ithCoefficient][0], MCcoefTable[iVar][ithCoefficient][1], MCcoefTable[iVar][ithCoefficient][2], AMCATNLOcoefTable[iVar][ithCoefficient][0], AMCATNLOcoefTable[iVar][ithCoefficient][1]);
          }
        }
        printf("\\hline \n");
        printf("\\end{tabular} \n");



        //printf(" Table of systematics %s%s: \n", NameTag.Data(),(inorm?"Norm":"") );
        printf(" Table of systematics for B Coef%s%s: \n",NameTag.Data(),(inorm?"Norm":""));

        printf("\\begin{tabular}{l");
        for (int iVar = 0; iVar < nBvar; ++iVar)
        {
          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            printf("c");
          }
        }
        printf("} \n");
        printf("\\hline \n");
        printf("\\textbf{Source} & \\multicolumn{%d}{c}{\\textbf{Uncertainty}} \\\\ \n", numCoefficients * nBvar);
        for (int iVar = 0; iVar < nBvar; ++iVar) {
          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            printf(" & $%s$", CoefNames.at(iVar * numCoefficients + ithCoefficient).Data());
          }
        }
        printf(" \\\\ \n");
        for (unsigned int isyst = 0; isyst < SystNames.size(); ++isyst)
        {
          printf("%s",SystNames[isyst].c_str());

          for (int iVar = 0; iVar < nBvar; ++iVar)
          {
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              printf(" & %2.3f", systTable[i_tag * 2 + inorm][iVar][ithCoefficient][isyst]);
            }
          }

          printf(" \\\\ \n");

          if(isyst == 8  || isyst == nsyst-3-backgroundnames.size()-flatsystnames.size() || isyst == nsyst-1-backgroundnames.size()-flatsystnames.size() || isyst == nsyst-1-flatsystnames.size() || int(isyst) == nsyst-1 || int(isyst) == nsyst || int(isyst) >= nsyst+3 )  printf("\\hline \n");

        }
        for (unsigned int isyst = 0; isyst < combsystnames.size(); ++isyst)
        {
          printf("%s",combsystnames[isyst].Data());

          for (int iVar = 0; iVar < nBvar; ++iVar)
          {
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              printf(" & %2.3f", combsystTable[i_tag * 2 + inorm][iVar][ithCoefficient][isyst]);
            }
          }
          printf(" \\\\ \n");

        }
        printf("\\end{tabular} \n");


        printf(" Table of systematics for C Coef%s%s: \n",NameTag.Data(),(inorm?"Norm":""));

        printf("\\begin{tabular}{l");
        for (int iVar = nBvar; iVar < nBvar + nCvar; ++iVar)
        {
          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            printf("c");
          }
        }
        printf("} \n");
        printf("\\hline \n");
        printf("\\textbf{Source} & \\multicolumn{%d}{c}{\\textbf{Uncertainty}} \\\\ \n",nCvar);
        for (int iVar = 0; iVar < nCvar; ++iVar) {
          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            printf(" & $%s$", CoefNames.at(iVar * numCoefficients + ithCoefficient).Data());
          }
        }
        printf(" \\\\ \n");
        for (unsigned int isyst = 0; isyst < SystNames.size(); ++isyst)
        {
          printf("%s",SystNames[isyst].c_str());

          for (int iVar = 0; iVar < nCvar; ++iVar) {
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              printf(" & %2.3f", systTable[i_tag * 2 + inorm][iVar][ithCoefficient][isyst]);
            }
          }

          printf(" \\\\ \n");

          if(isyst == 8  || isyst == nsyst-3-backgroundnames.size()-flatsystnames.size() || isyst == nsyst-1-backgroundnames.size()-flatsystnames.size() || isyst == nsyst-1-flatsystnames.size() || int(isyst) == nsyst-1 || int(isyst) == nsyst || int(isyst) >= nsyst+3 )  printf("\\hline \n");

        }
        for (unsigned int isyst = 0; isyst < combsystnames.size(); ++isyst)
        {
          printf("%s",combsystnames[isyst].Data());

          for (int iVar = 0; iVar < nCvar; ++iVar) {
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
                printf(" & %2.3f", combsystTable[i_tag * 2 + inorm][iVar][ithCoefficient][isyst]);
            }
          }
          printf(" \\\\ \n");

        }
        printf("\\end{tabular} \n");


        printf(" Table of results for BinSum%s%s: \n",NameTag.Data(),(inorm?"Norm":""));

        printf("\\begin{tabular}{l | c c c} \n");
        printf("\\hline \n");
        printf("\\textbf{BinSum} & \\textbf{Measured} & \\textbf{\\Powhegvtwo} & \\textbf{NLO calculation}  \\\\ \n");
        for (unsigned int iVar = 0; iVar < fVariableNames.size(); ++iVar)
        {
          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            printf("$%s$ & $%2.1f \\pm %2.1f$ &  &   \\\\ \n", CoefNames.at(iVar * numCoefficients + ithCoefficient).Data(), binsumTable[i_tag * 2 + inorm][iVar][numCoefficients], binsumsystTable[i_tag * 2 + inorm][iVar][numCoefficients][SystNames.size() - 1]);
          }
        }
        printf("\\hline \n");
        printf("\\end{tabular} \n");

        //printf(" Table of systematics %s%s: \n", NameTag.Data(),(inorm?"Norm":"") );
        printf(" Table of systematics for B BinSum%s%s: \n",NameTag.Data(),(inorm?"Norm":""));

        printf("\\begin{tabular}{l");
        for (int iVar = 0; iVar < nBvar; ++iVar)
        {
          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            printf("c");
          }
        }
        printf("} \n");
        printf("\\hline \n");
        printf("\\textbf{Source} & \\multicolumn{%d}{c}{\\textbf{Uncertainty}} \\\\ \n", numCoefficients * nBvar);
        for (int iVar = 0; iVar < nBvar; ++iVar) {
          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            printf(" & $%s$", CoefNames.at(iVar * numCoefficients + ithCoefficient).Data());
          }
        }
        printf(" \\\\ \n");
        for (unsigned int isyst = 0; isyst < SystNames.size(); ++isyst)
        {
          printf("%s",SystNames[isyst].c_str());

          for (int iVar = 0; iVar < nBvar; ++iVar)
          {
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              printf(" & %2.3f", binsumsystTable[i_tag * 2 + inorm][iVar][ithCoefficient][isyst]);
            }
          }

          printf(" \\\\ \n");

          if(isyst == 8  || isyst == nsyst-3-backgroundnames.size()-flatsystnames.size() || isyst == nsyst-1-backgroundnames.size()-flatsystnames.size() || isyst == nsyst-1-flatsystnames.size() || int(isyst) == nsyst-1 || int(isyst) == nsyst || int(isyst) >= nsyst+3 )  printf("\\hline \n");

        }
        printf("\\end{tabular} \n");


        printf(" Table of systematics for C BinSum%s%s: \n",NameTag.Data(),(inorm?"Norm":""));

        printf("\\begin{tabular}{l");
        for (int iVar = 0; iVar < nCvar; ++iVar) {
          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            printf("c");
          }
        }
        printf("} \n");
        printf("\\hline \n");
        printf("\\textbf{Source} & \\multicolumn{%d}{c}{\\textbf{Uncertainty}} \\\\ \n", numCoefficients * nCvar);
        for (int iVar = 0; iVar < nCvar; ++iVar) {
          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            printf(" & $%s$", CoefNames.at(iVar * numCoefficients + ithCoefficient).Data());
          }
        }
        printf(" \\\\ \n");
        for (unsigned int isyst = 0; isyst < SystNames.size(); ++isyst)
        {
          printf("%s",SystNames[isyst].c_str());

          for (int iVar = 0; iVar < nCvar; ++iVar) {
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              printf(" & %2.3f", binsumsystTable[i_tag * 2 + inorm][iVar][ithCoefficient][isyst]);
            }
          }

          printf(" \\\\ \n");

          if(isyst == 8  || isyst == nsyst-3-backgroundnames.size()-flatsystnames.size() || isyst == nsyst-1-backgroundnames.size()-flatsystnames.size() || isyst == nsyst-1-flatsystnames.size() || int(isyst) == nsyst-1 || int(isyst) == nsyst || int(isyst) >= nsyst+3 )  printf("\\hline \n");

        }
        printf("\\end{tabular} \n");

      }
    }

    // Write yaml table for rebinnedA coefficients
    gSystem->mkdir("yaml", true);
    for (int inorm = 0; inorm <= 1; ++inorm)
    {
      FILE *yamlfile_coef = fopen(Form("yaml/coefficients%s.yaml",(inorm?"norm":"")), "w");
      fprintf(yamlfile_coef, "dependent_variables:\n");
      if(inorm) fprintf(yamlfile_coef, "- header: {name: 'Coefficient'}\n");
      else fprintf(yamlfile_coef, "- header: {name: 'Coefficient'}\n");
      fprintf(yamlfile_coef, "  qualifiers:\n");
      fprintf(yamlfile_coef, "  - {name: RE, value: P P --> TOP TOPBAR X}\n");
      fprintf(yamlfile_coef, "  - {name: SQRT(S), units: GEV, value: '13000'}\n");
      fprintf(yamlfile_coef, "  values:\n");
      for (unsigned int iVar = 0; iVar < fVariableNames.size(); ++iVar)
      {
        for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
          fprintf(yamlfile_coef, "  - errors:\n");
          fprintf(yamlfile_coef, "    - {label: stat, symerror: %2.6f}\n", systTable[0 * 2 + inorm][iVar][ithCoefficient][SystNames.size() - 2]);
          fprintf(yamlfile_coef, "    - {label: sys, symerror: %2.6f}\n", systTable[0 * 2 + inorm][iVar][ithCoefficient][SystNames.size() - 6]);
          fprintf(yamlfile_coef, "    value: %2.6f\n", coefTable[0 * 2 + inorm][iVar][ithCoefficient]);
        }
      }
      fprintf(yamlfile_coef, "independent_variables:\n");
      fprintf(yamlfile_coef, "- header: {name: 'Coefficient'}\n");
      fprintf(yamlfile_coef, "  values:\n");
      for (unsigned int iVar = 0; iVar < fVariableNames.size(); ++iVar)
      {
        for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
          fprintf(yamlfile_coef, "  - {value: '$%s$'}\n", CoefNames.at(iVar * numCoefficients + ithCoefficient).Data());
        }
      }
      fclose(yamlfile_coef);

      if(inorm==1) {
        FILE *yamlfile_fSM = fopen(Form("yaml/fSM%s.yaml",(inorm?"norm":"")), "w");
        fprintf(yamlfile_fSM, "dependent_variables:\n");
        if(inorm) fprintf(yamlfile_fSM, "- header: {name: '$f_{\\mathrm{SM}}$'}\n");
        else fprintf(yamlfile_fSM, "- header: {name: '$f_{\\mathrm{SM}}$'}\n");
        fprintf(yamlfile_fSM, "  qualifiers:\n");
        fprintf(yamlfile_fSM, "  - {name: RE, value: P P --> TOP TOPBAR X}\n");
        fprintf(yamlfile_fSM, "  - {name: SQRT(S), units: GEV, value: '13000'}\n");
        fprintf(yamlfile_fSM, "  values:\n");
        for (unsigned int iVar = 0; iVar < fSMTableNames.size(); ++iVar)
        {
          fprintf(yamlfile_fSM, "  - errors:\n");
          fprintf(yamlfile_fSM, "    - {label: stat, symerror: %2.6f}\n", fSMTableStatErr.at(iVar) );
          fprintf(yamlfile_fSM, "    - {label: sys, symerror: %2.6f}\n", fSMTableSystErr.at(iVar) );
          fprintf(yamlfile_fSM, "    value: %2.6f\n", fSMTable.at(iVar) );
        }
        fprintf(yamlfile_fSM, "independent_variables:\n");
        fprintf(yamlfile_fSM, "- header: {name: 'Coefficient'}\n");
        fprintf(yamlfile_fSM, "  values:\n");
        for (unsigned int iVar = 0; iVar < fSMTableNames.size(); ++iVar)
        {
          fprintf(yamlfile_fSM, "  - {value: '$%s$'}\n", fSMTableNames[iVar].Data());
        }
        fclose(yamlfile_fSM);
      }

    }


    gStyle->SetPadLeftMargin(0.11);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadRightMargin(0.11);
    gStyle->SetPadTopMargin(0.09);
    gStyle->SetLabelOffset(0.006, "XY");

    // *************
    // AllVars plots
    // *************

    TFile *bootstrapfile = new TFile("BootstrapCorrelationMatrices.root","READ");
    TH1D* hnom_AllVar_copy[4];

    TString SystFilename_AllVars = TUNFFILE+"Systematics_AllVars.root";
    TFile *systoutfile_AllVars   = new TFile(SystFilename_AllVars.Data(),"RECREATE");
    systoutfile_AllVars->cd();

    NameTag = "_rebinnedA";

    std::vector<TMatrixD> SystCovMats_allVars;
    std::vector<TMatrixD> StatCovMats_allVars;
    std::vector<TMatrixD> TotalCovMats_allVars;
    std::vector<TMatrixD> StatCovMatsCoef_allVars;
    std::vector<TMatrixD> SystCovMatsCoef_allVars;
    std::vector<TMatrixD> TotalCovMatsCoef_allVars;
    std::vector<TMatrixD> StatCovMatsBinSum_allVars;
    std::vector<TMatrixD> SystCovMatsBinSum_allVars;
    std::vector<TMatrixD> TotalCovMatsBinSum_allVars;
    std::vector<TMatrixD> InvSystCovMats_allVars;
    std::vector<TMatrixD> InvStatCovMats_allVars;
    std::vector<TMatrixD> InvTotalCovMats_allVars;

    std::vector<TMatrixD> m_totalstatsystcovmat_AllVar_5bins;
    std::vector<TMatrixD> m_totalstatcovmat_AllVar_5bins;
    std::vector<TMatrixD> m_totalsystcovmat_AllVar_5bins;


    for (int i_tag = 0; i_tag < 2; ++i_tag)
    {
      if(i_tag==1) NameTag = "_rebinnedB";

      TH2D* htotalsystcovmat_AllVar[2];
      for (int inorm = 0; inorm <= 1; ++inorm)
      {
        
        // Loop over all the input systematics
        std::vector<TH2D*>hcomponentsystcovmat_AllVar;

        for(int syst = 0; syst<nsyst; syst++){
          TH1D* hnom_AllVar;
          TH1D* htot_reldelta_AllVar;

          MatrixUnfControl::GetDeltasAllVars(hnom_AllVar, htot_reldelta_AllVar, CHAN, tempVectorOfInputSystematics[syst], TUNFFILE, inorm, NameTag);
          TH2D* hcovtemp = (TH2D*) MatrixUnfControl::GetCombinedCovarianceMatrixAllVars(hnom_AllVar, htot_reldelta_AllVar, CHAN, tempVectorOfInputSystematics[syst], TUNFFILE, inorm, NameTag);

          int nbinsrhoi  = hcovtemp->GetNbinsX()/CoefNames2.size();

          if(syst==0) {
            hnom_AllVar_copy[i_tag*2+inorm] = (TH1D*) hnom_AllVar->Clone( Form("Unfolded_AllVar%s%s",(inorm?"Norm":""),NameTag.Data()) );
            //hnom_AllVar_copy[i_tag*2+inorm]->SetDirectory(0);
            hnom_AllVar_copy[i_tag*2+inorm]->SetTitle("");
            hnom_AllVar_copy[i_tag*2+inorm]->GetXaxis()->SetTitle("Bin of normalized differential cross section");
            hnom_AllVar_copy[i_tag*2+inorm]->GetXaxis()->SetTitleOffset(2.07);
            hnom_AllVar_copy[i_tag*2+inorm]->GetXaxis()->SetTickLength(0.);
            hnom_AllVar_copy[i_tag*2+inorm]->GetXaxis()->SetTitleSize(0.03);
            hnom_AllVar_copy[i_tag*2+inorm]->GetXaxis()->SetLabelSize(0.03);
            hnom_AllVar_copy[i_tag*2+inorm]->GetYaxis()->SetTitle(inorm?"Fraction of events":"d#sigma/dX [pb]");
            //hnom_AllVar_copy[i_tag*2+inorm]->GetYaxis()->SetTitleOffset(2.02);
            //hnom_AllVar_copy[i_tag*2+inorm]->GetYaxis()->SetTickLength(0.);
            hnom_AllVar_copy[i_tag*2+inorm]->GetYaxis()->SetTitleSize(0.03);
            hnom_AllVar_copy[i_tag*2+inorm]->GetYaxis()->SetLabelSize(0.03);

            for (unsigned int i_bin=1;i_bin<=CoefNames2.size()*nbinsrhoi;i_bin++) {
              if (i_bin%nbinsrhoi == 3) {
                hnom_AllVar_copy[i_tag*2+inorm]->GetXaxis()->SetBinLabel(i_bin,CoefNames2.at((i_bin-1)/nbinsrhoi).Data());
              }
              else {
                hnom_AllVar_copy[i_tag*2+inorm]->GetXaxis()->SetBinLabel(i_bin,"");
              }
            }
          }

          hcovtemp->SetTitle("");
          hcovtemp->GetXaxis()->SetTitle("Bin of normalized differential cross section");
          hcovtemp->GetXaxis()->SetTitleOffset(2.07);
          hcovtemp->GetXaxis()->SetTickLength(0.);
          hcovtemp->GetXaxis()->SetTitleSize(0.03);
          hcovtemp->GetXaxis()->SetLabelSize(0.03);
          hcovtemp->GetYaxis()->SetTitle("Bin of normalized differential cross section");
          hcovtemp->GetYaxis()->SetTitleOffset(2.02);
          hcovtemp->GetYaxis()->SetTickLength(0.);
          hcovtemp->GetYaxis()->SetTitleSize(0.03);
          hcovtemp->GetYaxis()->SetLabelSize(0.03);
          hcovtemp->GetZaxis()->SetTitleSize(0.03);
          hcovtemp->GetZaxis()->SetLabelSize(0.03);

          for (unsigned int i_bin=1; i_bin<=CoefNames2.size()*nbinsrhoi; i_bin++) {
            if (i_bin%nbinsrhoi == 3) {
              hcovtemp->GetXaxis()->SetBinLabel(i_bin,CoefNames2.at((i_bin-1)/nbinsrhoi).Data());
              hcovtemp->GetYaxis()->SetBinLabel(i_bin,CoefNames2.at((i_bin-1)/nbinsrhoi).Data());
            }
            else {
              hcovtemp->GetXaxis()->SetBinLabel(i_bin,"");
              hcovtemp->GetYaxis()->SetBinLabel(i_bin,"");
            }
          }

          if( tempVectorOfInputSystematics[syst] != "AMCATNLOFXFX" && tempVectorOfInputSystematics[syst] != "POWHEGV2HERWIG"  ) hcomponentsystcovmat_AllVar.push_back( hcovtemp );
          if(syst==0) { 
            htotalsystcovmat_AllVar[inorm] = (TH2D*) hcovtemp->Clone( Form("%sSystCovMatrix_AllVar%s%s","Total",(inorm?"Norm":""),NameTag.Data()) ); 
            htotalsystcovmat_AllVar[inorm]->SetDirectory(0);
          }
          else if( tempVectorOfInputSystematics[syst] != "AMCATNLOFXFX" && tempVectorOfInputSystematics[syst] != "POWHEGV2HERWIG"  ) htotalsystcovmat_AllVar[inorm]->Add( hcovtemp );

          calculateCoefficientSystCorrmAllVars(TUNFFILE, hcovtemp, hist_acombfactor_allVar[i_tag*2+inorm], NameTag, inorm);

          systoutfile_AllVars->cd();
          hcovtemp->Write();

          delete hnom_AllVar;
          delete htot_reldelta_AllVar;
        }

        int nbinsrhoi = htotalsystcovmat_AllVar[inorm]->GetNbinsX()/CoefNames2.size();

        TCanvas *cAllVars = new TCanvas(Form("cAllVars%s%s",(inorm?"Norm":""),NameTag.Data()), Form("cAllVars%s%s",(inorm?"Norm":""),NameTag.Data()));
        htotalsystcovmat_AllVar[inorm]->Draw("colz");
        cAllVars->Update();

        //TGaxis *labels = new TGaxis(0,cAllVars->GetUymin(),CoefNames2.size()*nbinsrhoi,cAllVars->GetUymin(),0,CoefNames2.size()*nbinsrhoi,ndivisions,"U");
        TF1 *f1=new TF1("f1","x",0,CoefNames2.size());
        TGaxis *labelsx  = new TGaxis(0,cAllVars->GetUymin(),CoefNames2.size()*nbinsrhoi,cAllVars->GetUymin(),"f1",CoefNames2.size(),"SU+");
        TGaxis *labelsx2 = new TGaxis(0,cAllVars->GetUymax(),CoefNames2.size()*nbinsrhoi,cAllVars->GetUymax(),"f1",CoefNames2.size(),"SU-");
        labelsx->SetTickSize(0.01);
        labelsx2->SetTickSize(0.01);
        labelsx->Draw();
        labelsx2->Draw();
        TGaxis *labelsy  = new TGaxis(cAllVars->GetUxmin(),0,cAllVars->GetUxmin(),CoefNames2.size()*nbinsrhoi,"f1",CoefNames2.size(),"SU-");
        TGaxis *labelsy2 = new TGaxis(cAllVars->GetUxmax(),0,cAllVars->GetUxmax(),CoefNames2.size()*nbinsrhoi,"f1",CoefNames2.size(),"SU+");
        labelsy->SetTickSize(0.01);
        labelsy2->SetTickSize(0.01);
        labelsy->Draw();
        labelsy2->Draw();
        DrawCMSLabels(0.04, fLumi);
        if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

        cAllVars->SaveAs(Form("%s/TotalSystCovMatrix%s%s_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalSystCovMatrix%s%s_AllVars.root",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalSystCovMatrix%s%s_AllVars.C",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));


        TH2D* h_RhoCCSyst = calculateCoefficientSystCorrmAllVars(TUNFFILE, htotalsystcovmat_AllVar[inorm], hist_acombfactor_allVar[i_tag*2+inorm], NameTag, inorm);
        h_RhoCCSyst->SetTitle("");
        h_RhoCCSyst->Draw("colz");
        cAllVars->Update();
        h_RhoCCSyst->GetXaxis()->SetTitle("Coefficient");
        h_RhoCCSyst->GetXaxis()->SetTitleOffset(2.07);
        h_RhoCCSyst->GetXaxis()->SetTitleSize(0.03);
        h_RhoCCSyst->GetXaxis()->SetLabelSize(0.03);
        h_RhoCCSyst->GetYaxis()->SetTitle("Coefficient");
        h_RhoCCSyst->GetYaxis()->SetTitleOffset(2.02);
        h_RhoCCSyst->GetYaxis()->SetTitleSize(0.03);
        h_RhoCCSyst->GetYaxis()->SetLabelSize(0.03);
        h_RhoCCSyst->GetZaxis()->SetTitleSize(0.03);
        h_RhoCCSyst->GetZaxis()->SetLabelSize(0.03);
        h_RhoCCSyst->GetZaxis()->SetRangeUser(-1,1);
        for (unsigned int i_bin=1;i_bin<=CoefNames2.size();i_bin++) {
          h_RhoCCSyst->GetXaxis()->SetBinLabel(i_bin,CoefNames2.at(i_bin-1).Data());
          h_RhoCCSyst->GetYaxis()->SetBinLabel(i_bin,CoefNames2.at(i_bin-1).Data());
        }
        h_RhoCCSyst->GetXaxis()->LabelsOption("v");
        //h_RhoCCSyst->GetYaxis()->LabelsOption("v");
        DrawCMSLabels(0.035);
        if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);
        cAllVars->SaveAs(Form("%s/TotalSystCorrMatrix%s%s_Coefficient_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalSystCorrMatrix%s%s_Coefficient_AllVars.root",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalSystCorrMatrix%s%s_Coefficient_AllVars.C",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));

        TH2D* h_RhoBSBSSyst = calculateBinSumSystCorrmAllVars(TUNFFILE, htotalsystcovmat_AllVar[inorm], hist_acombfactor_allVar[i_tag*2+inorm], NameTag, inorm);
        h_RhoBSBSSyst->SetTitle("");
        h_RhoBSBSSyst->Draw("colz");
        cAllVars->Update();
        h_RhoBSBSSyst->GetXaxis()->SetTitle("Bin Sum");
        h_RhoBSBSSyst->GetXaxis()->SetTitleOffset(2.07);
        h_RhoBSBSSyst->GetXaxis()->SetTitleSize(0.03);
        h_RhoBSBSSyst->GetXaxis()->SetLabelSize(0.03);
        h_RhoBSBSSyst->GetYaxis()->SetTitle("Bin Sum");
        h_RhoBSBSSyst->GetYaxis()->SetTitleOffset(2.02);
        h_RhoBSBSSyst->GetYaxis()->SetTitleSize(0.03);
        h_RhoBSBSSyst->GetYaxis()->SetLabelSize(0.03);
        h_RhoBSBSSyst->GetZaxis()->SetTitleSize(0.03);
        h_RhoBSBSSyst->GetZaxis()->SetLabelSize(0.02);
        //h_RhoBSBSSyst->GetZaxis()->SetRangeUser(-1,1);
        for (unsigned int i_bin=1;i_bin<=CoefNames2.size();i_bin++) {
          h_RhoBSBSSyst->GetXaxis()->SetBinLabel(i_bin,CoefNames2.at(i_bin-1).Data());
          h_RhoBSBSSyst->GetYaxis()->SetBinLabel(i_bin,CoefNames2.at(i_bin-1).Data());
        }
        h_RhoBSBSSyst->GetXaxis()->LabelsOption("v");
        //h_RhoBSBSSyst->GetYaxis()->LabelsOption("v");
        DrawCMSLabels(0.04, fLumi);
        if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);
        cAllVars->SaveAs(Form("%s/TotalSystCorrMatrix%s%s_BinSum_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalSystCorrMatrix%s%s_BinSum_AllVars.root",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalSystCorrMatrix%s%s_BinSum_AllVars.C",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));

        TMatrixD totalsystcorrmat_AllVar = MatrixUnf::TH2toTMatrixD(htotalsystcovmat_AllVar[inorm]);
        totalsystcorrmat_AllVar.NormByDiag(TMatrixDDiag(totalsystcorrmat_AllVar));
        TH2D* htotalsystcorrmat_AllVar   = MatrixUnf::TMatrixDtoTH2(totalsystcorrmat_AllVar,"temp");
        htotalsystcorrmat_AllVar->SetTitle("");
        htotalsystcorrmat_AllVar->Draw("colz");
        cAllVars->Update();
        htotalsystcorrmat_AllVar->GetXaxis()->SetTitle("Bin of normalized differential cross section");
        htotalsystcorrmat_AllVar->GetXaxis()->SetTitleOffset(2.07);
        htotalsystcorrmat_AllVar->GetXaxis()->SetTickLength(0.);
        htotalsystcorrmat_AllVar->GetXaxis()->SetTitleSize(0.03);
        htotalsystcorrmat_AllVar->GetXaxis()->SetLabelSize(0.03);
        htotalsystcorrmat_AllVar->GetYaxis()->SetTitle("Bin of normalized differential cross section");
        htotalsystcorrmat_AllVar->GetYaxis()->SetTitleOffset(2.02);
        htotalsystcorrmat_AllVar->GetYaxis()->SetTickLength(0.);
        htotalsystcorrmat_AllVar->GetYaxis()->SetTitleSize(0.03);
        htotalsystcorrmat_AllVar->GetYaxis()->SetLabelSize(0.03);
        htotalsystcorrmat_AllVar->GetZaxis()->SetTitleSize(0.03);
        htotalsystcorrmat_AllVar->GetZaxis()->SetLabelSize(0.03);
        htotalsystcorrmat_AllVar->GetZaxis()->SetRangeUser(-1,1);

        for (unsigned int i_bin=1;i_bin<=CoefNames2.size()*nbinsrhoi;i_bin++) {
          if (i_bin%nbinsrhoi == 3) {
            htotalsystcorrmat_AllVar->GetXaxis()->SetBinLabel(i_bin,CoefNames2.at((i_bin-1)/nbinsrhoi).Data());
            htotalsystcorrmat_AllVar->GetYaxis()->SetBinLabel(i_bin,CoefNames2.at((i_bin-1)/nbinsrhoi).Data());
          }
          else {
            htotalsystcorrmat_AllVar->GetXaxis()->SetBinLabel(i_bin,"");
            htotalsystcorrmat_AllVar->GetYaxis()->SetBinLabel(i_bin,"");
          }
        }
        labelsx->Draw();
        labelsx2->Draw();
        labelsy->Draw();
        labelsy2->Draw();
        DrawCMSLabels(0.035, fLumi);
        if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);
        cAllVars->SaveAs(Form("%s/TotalSystCorrMatrix%s%s_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalSystCorrMatrix%s%s_AllVars.root",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalSystCorrMatrix%s%s_AllVars.C",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));


        bool doCheckInverse = kTRUE;

        if(doCheckInverse) {

          std::cout<<"Inverting total systematic cov and corr matrices "<<(inorm?"Norm":"")<<NameTag<<std::endl;
          TMatrixD invtotalsystcorrmat_AllVar = MatrixUnf::TH2toTMatrixD(htotalsystcovmat_AllVar[inorm]);
          TMatrixD invtotalsystcovmat_AllVar  = MatrixUnf::TH2toTMatrixD(htotalsystcovmat_AllVar[inorm]);
          invtotalsystcorrmat_AllVar.NormByDiag(TMatrixDDiag(invtotalsystcorrmat_AllVar));

          SystCovMats_allVars.push_back(invtotalsystcovmat_AllVar);

          invtotalsystcovmat_AllVar.Invert();
          invtotalsystcorrmat_AllVar.Invert();

          InvSystCovMats_allVars.push_back(invtotalsystcovmat_AllVar);

          invtotalsystcorrmat_AllVar.Draw("colz");
          cAllVars->SaveAs(Form("%s/InvSystCorrMatrix%s%s_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
          invtotalsystcorrmat_AllVar.NormByDiag(TMatrixDDiag(invtotalsystcorrmat_AllVar));
          invtotalsystcorrmat_AllVar.Draw("colz");
          cAllVars->SaveAs(Form("%s/InvSystCorrMatrixNBD%s%s_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));

          invtotalsystcovmat_AllVar.Draw("colz");
          cAllVars->SaveAs(Form("%s/InvSystCovMatrix%s%s_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
          invtotalsystcovmat_AllVar.NormByDiag(TMatrixDDiag(invtotalsystcovmat_AllVar));
          invtotalsystcovmat_AllVar.Draw("colz");
          cAllVars->SaveAs(Form("%s/InvSystCovMatrixNBD%s%s_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));

        }

        delete htotalsystcorrmat_AllVar;

        // Amandeep : GetNormalisedCovarianceMatrixAllVars this the plotting code fails 
        // FIXME    : The function has hardcoded nbinsrhoi = 6; which needs to be changed

        TH2D* htemptotalstatcorrmat_AllVar   = (TH2D*)bootstrapfile->Get("h_RhoBBCorrected"+NameTag);
        TMatrixD temptotalstatcorrmat_AllVar = MatrixUnf::TH2toTMatrixD(htemptotalstatcorrmat_AllVar);
        
        // The diagonal should be 1s already, but just in case of rounding errors etc:
        temptotalstatcorrmat_AllVar.NormByDiag(TMatrixDDiag(temptotalstatcorrmat_AllVar));
        TVectorD deltaAllVar     = MatrixUnf::TH1toTVectorD(hdeltasquaredStat_AllVar[i_tag*2]);
        TVectorD deltaAllVarNorm = MatrixUnf::TH1toTVectorD(hdeltasquaredStat_AllVar[i_tag*2+1]);

        // Convert back to covariance matrix
        // Amandeep : NormByDiag does a(i,j)*sqrt(abs*(v(i)*v(j)))  
        // where v is the TVectorD that is passed; hdeltasquaredStat_AllVar in this case

        temptotalstatcorrmat_AllVar.NormByDiag(deltaAllVar,0);
        
        TH2D* htotalstatcovmat_AllVar = MatrixUnf::TMatrixDtoTH2(temptotalstatcorrmat_AllVar,Form("TotalStatCovMatrix_AllVar%s%s",(inorm?"Norm":""),NameTag.Data()));
        if(inorm) {
          MatrixUnf::GetNormalisedCovarianceMatrixAllVars(hnom_AllVar_copy[i_tag*2], htotalstatcovmat_AllVar, htotalstatcovmat_AllVar);

          bool doUseNormdeltaVec = kFALSE;
          if(doUseNormdeltaVec) {
            TMatrixD temptotalstatcovmat_AllVar = MatrixUnf::TH2toTMatrixD(htotalstatcovmat_AllVar);
            temptotalstatcovmat_AllVar.NormByDiag(TMatrixDDiag(temptotalstatcovmat_AllVar));
            temptotalstatcovmat_AllVar.NormByDiag(deltaAllVarNorm,0);
            delete htotalstatcovmat_AllVar;
            htotalstatcovmat_AllVar = MatrixUnf::TMatrixDtoTH2(temptotalstatcovmat_AllVar,Form("TotalStatCovMatrix_AllVar%s%s",(inorm?"Norm":""),NameTag.Data()));
            MatrixUnf::GetNormalisedCovarianceMatrixAllVars(hnom_AllVar_copy[i_tag*2+inorm], htotalstatcovmat_AllVar, htotalstatcovmat_AllVar);
          }

          // TUnfold matrices for testing
          TMatrixD TUnfoldTotalStatCovM_AllVar  = MatrixUnf::TH2toTMatrixD(hTUnfoldTotalStatCorrM_AllVar[i_tag*2+inorm]);
          TUnfoldTotalStatCovM_AllVar.NormByDiag(deltaAllVar,0);
          TH2D* hTUnfoldTotalStatCovM_AllVar    = MatrixUnf::TMatrixDtoTH2(TUnfoldTotalStatCovM_AllVar,Form("hTUnfoldTotalStatCovM_AllVar_%d",i_tag*2+inorm));
          MatrixUnf::GetNormalisedCovarianceMatrixAllVars(hnom_AllVar_copy[i_tag*2], hTUnfoldTotalStatCovM_AllVar, hTUnfoldTotalStatCovM_AllVar);
          TMatrixD TUnfoldTotalStatCorrM_AllVar = MatrixUnf::TH2toTMatrixD(hTUnfoldTotalStatCovM_AllVar);
          delete hTUnfoldTotalStatCovM_AllVar;
          TUnfoldTotalStatCorrM_AllVar.NormByDiag(TMatrixDDiag(TUnfoldTotalStatCorrM_AllVar));
          delete hTUnfoldTotalStatCorrM_AllVar[i_tag*2+inorm];
          hTUnfoldTotalStatCorrM_AllVar[i_tag*2+inorm] = MatrixUnf::TMatrixDtoTH2(TUnfoldTotalStatCorrM_AllVar,Form("hTUnfoldTotalStatCorrM_AllVar_%d",i_tag*2+inorm));

          TMatrixD TUnfoldDataStatCovM_AllVar  = MatrixUnf::TH2toTMatrixD(hTUnfoldDataStatCorrM_AllVar[i_tag*2+inorm]);
          TUnfoldDataStatCovM_AllVar.NormByDiag(deltaAllVar,0);
          TH2D* hTUnfoldDataStatCovM_AllVar    = MatrixUnf::TMatrixDtoTH2(TUnfoldDataStatCovM_AllVar,Form("hTUnfoldDataStatCovM_AllVar_%d",i_tag*2+inorm));
          MatrixUnf::GetNormalisedCovarianceMatrixAllVars(hnom_AllVar_copy[i_tag*2], hTUnfoldDataStatCovM_AllVar, hTUnfoldDataStatCovM_AllVar);
          TMatrixD TUnfoldDataStatCorrM_AllVar = MatrixUnf::TH2toTMatrixD(hTUnfoldDataStatCovM_AllVar);
          delete hTUnfoldDataStatCovM_AllVar;
          TUnfoldDataStatCorrM_AllVar.NormByDiag(TMatrixDDiag(TUnfoldDataStatCorrM_AllVar));
          delete hTUnfoldDataStatCorrM_AllVar[i_tag*2+inorm];
          hTUnfoldDataStatCorrM_AllVar[i_tag*2+inorm] = MatrixUnf::TMatrixDtoTH2(TUnfoldDataStatCorrM_AllVar,Form("hTUnfoldDataStatCorrM_AllVar_%d",i_tag*2+inorm));

        }

        TH2D* h_RhoCCStat   = calculateCoefficientSystCorrmAllVars(TUNFFILE, htotalstatcovmat_AllVar, hist_acombfactor_allVar[i_tag*2+inorm], NameTag, inorm);
        TH2D* h_RhoBSBSStat = calculateBinSumSystCorrmAllVars(TUNFFILE, htotalstatcovmat_AllVar, hist_acombfactor_allVar[i_tag*2+inorm], NameTag, inorm);

        TMatrixD totalstatsystcovmat_AllVar = MatrixUnf::TH2toTMatrixD(htotalstatcovmat_AllVar);
        totalstatsystcovmat_AllVar         += MatrixUnf::TH2toTMatrixD(htotalsystcovmat_AllVar[inorm]);        
        TH2D* htotalstatsystcovmat_AllVar   = MatrixUnf::TMatrixDtoTH2(totalstatsystcovmat_AllVar,Form("TotalStatSystCovMatrix_AllVar%s%s",(inorm?"Norm":""),NameTag.Data()));


        TH2D* h_RhoCCStatSyst   = calculateCoefficientSystCorrmAllVars(TUNFFILE, htotalstatsystcovmat_AllVar, hist_acombfactor_allVar[i_tag*2+inorm], NameTag, inorm);
        TH2D* h_RhoBSBSStatSyst = calculateBinSumSystCorrmAllVars(TUNFFILE, htotalstatsystcovmat_AllVar, hist_acombfactor_allVar[i_tag*2+inorm], NameTag, inorm);


        TMatrixD invtotalstatsystBScormat = MatrixUnf::TH2toTMatrixD(h_RhoBSBSStatSyst);
        TMatrixD invtotalstatBScormat     = MatrixUnf::TH2toTMatrixD(h_RhoBSBSStat);
        TMatrixD invtotalsystBScormat     = MatrixUnf::TH2toTMatrixD(h_RhoBSBSSyst);

        invtotalstatsystBScormat.Invert();
        invtotalstatBScormat.Invert();
        invtotalsystBScormat.Invert();


        if(doCheckInverse) {
          std::cout<<"Inverting total stat+syst cov and corr matrices "<<(inorm?"Norm":"")<<NameTag<<std::endl;

          TMatrixD invtotalstatsystcovmat_AllVar = MatrixUnf::TH2toTMatrixD(htotalstatcovmat_AllVar);
          TMatrixD invtotalstatcovmat_AllVar     = MatrixUnf::TH2toTMatrixD(htotalstatcovmat_AllVar);
          TMatrixD invtotalsystcovmat_AllVar     = MatrixUnf::TH2toTMatrixD(htotalsystcovmat_AllVar[inorm]);

          invtotalstatsystcovmat_AllVar += MatrixUnf::TH2toTMatrixD(htotalsystcovmat_AllVar[inorm]);

          TMatrixD invtotalstatsystcorrmat_AllVar = invtotalstatsystcovmat_AllVar;
          invtotalstatsystcorrmat_AllVar.NormByDiag(TMatrixDDiag(invtotalstatsystcorrmat_AllVar));


          TotalCovMats_allVars.push_back(invtotalstatsystcovmat_AllVar);
          StatCovMats_allVars.push_back(invtotalstatcovmat_AllVar);

          /*
          TMatrixD m_RhoCCStatSyst = MatrixUnf::TH2toTMatrixD( calculateCoefficientSystCorrmAllVars(TUNFFILE, htotalstatsystcovmat_AllVar, hist_acombfactor_allVar[i_tag*2+inorm], NameTag, inorm) );
          TMatrixD m_CovCCStatSyst = m_RhoCCStatSyst;
          */
          TMatrixD m_RhoCCStat = MatrixUnf::TH2toTMatrixD( h_RhoCCStat );
          TMatrixD m_CovCCStat = m_RhoCCStat;

          TMatrixD m_RhoCCSyst = MatrixUnf::TH2toTMatrixD( h_RhoCCSyst );
          TMatrixD m_CovCCSyst = m_RhoCCSyst;


          //TVectorD v_Coef(fVariableNames.size());
          TVectorD v_CoefStat(fVariableNames.size() * numCoefficients);
          TVectorD v_CoefSyst(fVariableNames.size() * numCoefficients);
          for (unsigned int i_var = 0; i_var < fVariableNames.size(); ++i_var) {
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              // v_Coef(i_var) = pow(systTable[i_tag*2+inorm][i_var][SystNames.size()-1],2);
              v_CoefStat(i_var * numCoefficients + ithCoefficient) = pow(systTable[i_tag * 2 + inorm][i_var][ithCoefficient][SystNames.size() - 2], 2);
              v_CoefSyst(i_var * numCoefficients + ithCoefficient) = pow(systTable[i_tag * 2 + inorm][i_var][ithCoefficient][SystNames.size() - 6], 2);
            }
          }
          //m_CovCCStatSyst.NormByDiag(v_Coef,0);
          m_CovCCStat.NormByDiag(v_CoefStat,0);
          m_CovCCSyst.NormByDiag(v_CoefSyst,0);


          TMatrixD m_CovCCStatSyst = m_CovCCStat+m_CovCCSyst;
          TMatrixD m_RhoCCStatSyst = m_CovCCStatSyst;
          m_RhoCCStatSyst.NormByDiag(TMatrixDDiag(m_RhoCCStatSyst));


          TotalCovMatsCoef_allVars.push_back(m_CovCCStatSyst);
          StatCovMatsCoef_allVars.push_back(m_CovCCStat);
          SystCovMatsCoef_allVars.push_back(m_CovCCSyst);

          TMatrixD m_RhoBSBSStat = MatrixUnf::TH2toTMatrixD( h_RhoBSBSStat );
          TMatrixD m_CovBSBSStat = m_RhoBSBSStat;

          TMatrixD m_RhoBSBSSyst = MatrixUnf::TH2toTMatrixD( h_RhoBSBSSyst );
          TMatrixD m_CovBSBSSyst = m_RhoBSBSSyst;


          //TVectorD v_BinSum(fVariableNames.size());
          TVectorD v_BinSumStat(fVariableNames.size() * numCoefficients);
          TVectorD v_BinSumSyst(fVariableNames.size() * numCoefficients);
          for (unsigned int i_var = 0; i_var < fVariableNames.size(); ++i_var) {
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              // v_BinSum(i_var) = pow(binsumsystTable[i_tag*2+inorm][i_var][SystNames.size()-1],2);
              v_BinSumStat(i_var * numCoefficients + ithCoefficient) = pow(binsumsystTable[i_tag * 2 + inorm][i_var][ithCoefficient][SystNames.size() - 2], 2);
              v_BinSumSyst(i_var * numCoefficients + ithCoefficient) = pow(binsumsystTable[i_tag * 2 + inorm][i_var][ithCoefficient][SystNames.size() - 6], 2);
            }
          }
          //m_CovBSBSStatSyst.NormByDiag(v_BinSum,0);
          m_CovBSBSStat.NormByDiag(v_BinSumStat,0);
          m_CovBSBSSyst.NormByDiag(v_BinSumSyst,0);


          TMatrixD m_CovBSBSStatSyst = m_CovBSBSStat+m_CovBSBSSyst;
          TMatrixD m_RhoBSBSStatSyst = m_CovBSBSStatSyst;
          m_RhoBSBSStatSyst.NormByDiag(TMatrixDDiag(m_RhoBSBSStatSyst));


          TotalCovMatsBinSum_allVars.push_back(m_CovBSBSStatSyst);
          StatCovMatsBinSum_allVars.push_back(m_CovBSBSStat);
          SystCovMatsBinSum_allVars.push_back(m_CovBSBSSyst);

          const int nbin_shape = fVariableNames.size()*(nbinsrhoi-1);
          const int nbin_total = fVariableNames.size()*(nbinsrhoi);

          std::vector<TMatrixD> totalstatsystcovmats_AllVar_5bins;
          std::vector<TMatrixD> totalstatcovmats_AllVar_5bins;
          std::vector<TMatrixD> totalsystcovmats_AllVar_5bins;

          for (int i_ignore = 0; i_ignore < nbinsrhoi; ++i_ignore)
          {

            const int binToIgnore = i_ignore;
            TMatrixD totalstatsystcovmat_AllVar_5bins(nbin_shape, nbin_shape);
            TMatrixD totalstatcovmat_AllVar_5bins(nbin_shape, nbin_shape);
            TMatrixD totalsystcovmat_AllVar_5bins(nbin_shape, nbin_shape);

            int iShpR = 0, iShpC = 0;

            for (int iAbsR = 0; iAbsR < nbin_total; ++iAbsR) {
              if (iAbsR % nbinsrhoi == binToIgnore) continue;

              iShpC = 0;
              for (int iAbsC = 0; iAbsC < nbin_total; ++iAbsC) {
                if (iAbsC % nbinsrhoi == binToIgnore) continue;

                totalstatsystcovmat_AllVar_5bins(iShpR, iShpC) = invtotalstatsystcovmat_AllVar(iAbsR, iAbsC);
                totalstatcovmat_AllVar_5bins(iShpR, iShpC) = invtotalstatcovmat_AllVar(iAbsR, iAbsC);
                totalsystcovmat_AllVar_5bins(iShpR, iShpC) = invtotalsystcovmat_AllVar(iAbsR, iAbsC);
                ++iShpC;
              }

              ++iShpR;
            }

            totalstatsystcovmats_AllVar_5bins.push_back(totalstatsystcovmat_AllVar_5bins);
            totalstatcovmats_AllVar_5bins.push_back(totalstatcovmat_AllVar_5bins);
            totalsystcovmats_AllVar_5bins.push_back(totalsystcovmat_AllVar_5bins);

          }

          invtotalstatsystcovmat_AllVar.Invert();
          invtotalstatsystcorrmat_AllVar.Invert();

          InvTotalCovMats_allVars.push_back(invtotalstatsystcovmat_AllVar);

          invtotalstatsystcorrmat_AllVar.Draw("colz");
          cAllVars->SaveAs(Form("%s/InvStatPlusSystCorrMatrix1%s%s_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
          invtotalstatsystcorrmat_AllVar.NormByDiag(TMatrixDDiag(invtotalstatsystcorrmat_AllVar));
          invtotalstatsystcorrmat_AllVar.Draw("colz");
          cAllVars->SaveAs(Form("%s/InvStatPlusSystCorrMatrixNBD1%s%s_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));

          invtotalstatsystcovmat_AllVar.Draw("colz");
          cAllVars->SaveAs(Form("%s/InvStatPlusSystCovMatrix1%s%s_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
          invtotalstatsystcovmat_AllVar.NormByDiag(TMatrixDDiag(invtotalstatsystcovmat_AllVar));
          invtotalstatsystcovmat_AllVar.Draw("colz");
          cAllVars->SaveAs(Form("%s/InvStatPlusSystCovMatrixNBD1%s%s_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));


          for (int iremove = 0; iremove < int(fVariableNames.size()); ++iremove)
          {
            for (int i_ignore = 0; i_ignore < nbinsrhoi; ++i_ignore)
            {

              TMatrixD invtotalstatsystcovmat_AllVar_5bins_copy = totalstatsystcovmats_AllVar_5bins.at(i_ignore);

              invtotalstatsystcovmat_AllVar_5bins_copy.ResizeTo(nbin_shape-iremove*(nbinsrhoi-1),nbin_shape-iremove*(nbinsrhoi-1));

              int iShpR = 0, iShpC = 0;

              for (int iAbsR = 0; iAbsR < nbin_shape; ++iAbsR) {
                if (iAbsR >= nbin_shape-iremove*(nbinsrhoi-1)) continue;
                //if (iAbsR < iremove*(nbinsrhoi-1)) continue;

                iShpC = 0;
                for (int iAbsC = 0; iAbsC < nbin_shape; ++iAbsC) {
                  if (iAbsC >= nbin_shape-iremove*(nbinsrhoi-1)) continue;
                  //if (iAbsC < iremove*(nbinsrhoi-1)) continue;

                  invtotalstatsystcovmat_AllVar_5bins_copy(iShpR, iShpC) = totalstatsystcovmats_AllVar_5bins.at(i_ignore)(iAbsR, iAbsC);
                  ++iShpC;
                }

                ++iShpR;
              }


              TDecompSVD svd_invtotalstatsystcovmat_AllVar_5bins(invtotalstatsystcovmat_AllVar_5bins_copy);
              svd_invtotalstatsystcovmat_AllVar_5bins.Decompose();
              std::cout<<" invtotalstatsystcovmat_AllVar_5bins iremove= "<<iremove<<" i_ignore= "<<i_ignore<<" condition number: "<<svd_invtotalstatsystcovmat_AllVar_5bins.Condition()<<" "<<svd_invtotalstatsystcovmat_AllVar_5bins.GetCondition()<<std::endl;

              if(i_tag==0 && inorm==1) {
                //global correlation coefficient
                TH1D* hrho = new TH1D(Form("hrho%s%s_%i_%i",(inorm?"Norm":"Abs"),NameTag.Data(),iremove,i_ignore), Form("hrho%s%s_%i_%i",(inorm?"Norm":"Abs"),NameTag.Data(),iremove,i_ignore), nbin_shape-iremove*(nbinsrhoi-1),0,nbin_shape-iremove*(nbinsrhoi-1));
                hrho->SetDirectory(0);

                TMatrixD CovMatInv = invtotalstatsystcovmat_AllVar_5bins_copy;
                CovMatInv.Invert();

                for (int i_bin=0;i_bin<nbin_shape-iremove*(nbinsrhoi-1);i_bin++) {
                  hrho->SetBinContent(i_bin+1, sqrt(1.-1/(CovMatInv(i_bin,i_bin)*invtotalstatsystcovmat_AllVar_5bins_copy(i_bin,i_bin))) );
                }

                TCanvas *crho = new TCanvas(Form("crho%s%s_%i_%i",(inorm?"Norm":"Abs"),NameTag.Data(),iremove,i_ignore), Form("crho%s%s_%i_%i",(inorm?"Norm":"Abs"),NameTag.Data(),iremove,i_ignore));
                hrho->Draw();

                crho->SaveAs(Form("%s/Rho%s%s_%i_%i.pdf",PLOTFILE.Data(),(inorm?"Norm":"Abs"),NameTag.Data(),iremove,i_ignore));

              }
            }
          }

          cAllVars->cd();

          const int binToIgnore_const = 1;
          TMatrixD invtotalstatsystcovmat_AllVar_5bins = totalstatsystcovmats_AllVar_5bins.at(binToIgnore_const);
          m_totalstatsystcovmat_AllVar_5bins.push_back( totalstatsystcovmats_AllVar_5bins.at(binToIgnore_const) );
          m_totalstatcovmat_AllVar_5bins.push_back( totalstatcovmats_AllVar_5bins.at(binToIgnore_const) );
          m_totalsystcovmat_AllVar_5bins.push_back( totalsystcovmats_AllVar_5bins.at(binToIgnore_const) );

          invtotalstatsystcovmat_AllVar_5bins.Invert();

          invtotalstatsystcovmat_AllVar_5bins.Draw("colz");
          cAllVars->SaveAs(Form("%s/InvStatPlusSystCovMatrix1%s%s_AllVars_5bins.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
          invtotalstatsystcovmat_AllVar_5bins.NormByDiag(TMatrixDDiag(invtotalstatsystcovmat_AllVar_5bins));
          invtotalstatsystcovmat_AllVar_5bins.Draw("colz");
          cAllVars->SaveAs(Form("%s/InvStatPlusSystCovMatrixNBD1%s%s_AllVars_5bins.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));


        }

        htotalstatcovmat_AllVar->SetTitle("");
        htotalstatcovmat_AllVar->Draw("colz");
        cAllVars->Update();
        htotalstatcovmat_AllVar->GetXaxis()->SetTitle("Bin of normalized differential cross section");
        htotalstatcovmat_AllVar->GetXaxis()->SetTitleOffset(2.07);
        htotalstatcovmat_AllVar->GetXaxis()->SetTickLength(0.);
        htotalstatcovmat_AllVar->GetXaxis()->SetTitleSize(0.03);
        htotalstatcovmat_AllVar->GetXaxis()->SetLabelSize(0.03);
        htotalstatcovmat_AllVar->GetYaxis()->SetTitle("Bin of normalized differential cross section");
        htotalstatcovmat_AllVar->GetYaxis()->SetTitleOffset(2.02);
        htotalstatcovmat_AllVar->GetYaxis()->SetTickLength(0.);
        htotalstatcovmat_AllVar->GetYaxis()->SetTitleSize(0.03);
        htotalstatcovmat_AllVar->GetYaxis()->SetLabelSize(0.03);
        htotalstatcovmat_AllVar->GetZaxis()->SetTitleSize(0.03);
        htotalstatcovmat_AllVar->GetZaxis()->SetLabelSize(0.03);
        for (unsigned int i_bin=1;i_bin<=CoefNames2.size()*nbinsrhoi;i_bin++) {
          if (i_bin%nbinsrhoi == 3) {
            htotalstatcovmat_AllVar->GetXaxis()->SetBinLabel(i_bin,CoefNames2.at((i_bin-1)/nbinsrhoi).Data());
            htotalstatcovmat_AllVar->GetYaxis()->SetBinLabel(i_bin,CoefNames2.at((i_bin-1)/nbinsrhoi).Data());
          }
          else {
            htotalstatcovmat_AllVar->GetXaxis()->SetBinLabel(i_bin,"");
            htotalstatcovmat_AllVar->GetYaxis()->SetBinLabel(i_bin,"");
          }
        }
        labelsx->Draw();
        labelsx2->Draw();
        labelsy->Draw();
        labelsy2->Draw();
        DrawCMSLabels(0.04, fLumi);
        if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

        cAllVars->SaveAs(Form("%s/TotalStatCovMatrix%s%s_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalStatCovMatrix%s%s_AllVars.root",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalStatCovMatrix%s%s_AllVars.C",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));



        if(i_tag==0) {

          //write yaml tables for rebinnedA results

          gSystem->mkdir("yaml", true);

          FILE *yamlfile_systcov = fopen(Form("yaml/systcov_%sdiffxsec_AllVars.yaml",(inorm?"norm":"")), "w");
          fprintf(yamlfile_systcov, "dependent_variables:\n");
          if(inorm) fprintf(yamlfile_systcov, "- header: {name: 'Systematic covariance for all normalized differential cross section bins'}\n");
          else fprintf(yamlfile_systcov, "- header: {name: 'Systematic covariance for all differential cross section bins [pb$^{2}$]'}\n");
          fprintf(yamlfile_systcov, "  qualifiers:\n");
          fprintf(yamlfile_systcov, "  - {name: RE, value: P P --> TOP TOPBAR X}\n");
          fprintf(yamlfile_systcov, "  - {name: SQRT(S), units: GEV, value: '13000'}\n");
          fprintf(yamlfile_systcov, "  values:\n");
          for (int j_bin = 0; j_bin < htotalsystcovmat_AllVar[inorm]->GetNbinsX(); ++j_bin) for (int i_bin = 0; i_bin < htotalsystcovmat_AllVar[inorm]->GetNbinsX(); ++i_bin)
          {
            //fprintf(yamlfile_systcov, "  - errors: []\n");
            fprintf(yamlfile_systcov, "  - {value: %2.10g}\n", htotalsystcovmat_AllVar[inorm]->GetBinContent(i_bin+1,j_bin+1)/(divBinWidth?(hbinwidth_AllVar->GetBinContent(i_bin+1)*hbinwidth_AllVar->GetBinContent(j_bin+1)):1.) );
          }
          fprintf(yamlfile_systcov, "independent_variables:\n");
          fprintf(yamlfile_systcov, "- header: {name: 'Bin'}\n");
          fprintf(yamlfile_systcov, "  values:\n");
          for (int i_row = 0; i_row < htotalsystcovmat_AllVar[inorm]->GetNbinsX(); ++i_row) for (unsigned int i_var = 0; i_var < fVariableNames.size(); ++i_var) for (int i_bin = 0; i_bin < nbinsrhoi; ++i_bin) {
            fprintf(yamlfile_systcov, "  - {value: '$%s$ bin %i'}\n", (MatrixUnf::ObsName(fVariableNames[i_var],kTRUE)).Data(), i_bin+1);
          }
          fprintf(yamlfile_systcov, "- header: {name: 'Bin'}\n");
          fprintf(yamlfile_systcov, "  values:\n");
          for (unsigned int i_var = 0; i_var < fVariableNames.size(); ++i_var) for (int i_bin = 0; i_bin < nbinsrhoi; ++i_bin) for (int i_row = 0; i_row < htotalsystcovmat_AllVar[inorm]->GetNbinsX(); ++i_row) {
            fprintf(yamlfile_systcov, "  - {value: '$%s$ bin %i'}\n", (MatrixUnf::ObsName(fVariableNames[i_var],kTRUE)).Data(), i_bin+1);
          }
          fclose(yamlfile_systcov);


          FILE *yamlfile_statcov = fopen(Form("yaml/statcov_%sdiffxsec_AllVars.yaml",(inorm?"norm":"")), "w");
          fprintf(yamlfile_statcov, "dependent_variables:\n");
          if(inorm) fprintf(yamlfile_statcov, "- header: {name: 'Statistical covariance for all normalized differential cross section bins'}\n");
          else fprintf(yamlfile_statcov, "- header: {name: 'Statistical covariance for all differential cross section bins [pb$^{2}$]'}\n");
          fprintf(yamlfile_statcov, "  qualifiers:\n");
          fprintf(yamlfile_statcov, "  - {name: RE, value: P P --> TOP TOPBAR X}\n");
          fprintf(yamlfile_statcov, "  - {name: SQRT(S), units: GEV, value: '13000'}\n");
          fprintf(yamlfile_statcov, "  values:\n");
          for (int j_bin = 0; j_bin < htotalstatcovmat_AllVar->GetNbinsX(); ++j_bin) for (int i_bin = 0; i_bin < htotalstatcovmat_AllVar->GetNbinsX(); ++i_bin)
          {
            //fprintf(yamlfile_statcov, "  - errors: []\n");
            fprintf(yamlfile_statcov, "  - {value: %2.10g}\n", htotalstatcovmat_AllVar->GetBinContent(i_bin+1,j_bin+1)/(divBinWidth?(hbinwidth_AllVar->GetBinContent(i_bin+1)*hbinwidth_AllVar->GetBinContent(j_bin+1)):1.) );
          }
          fprintf(yamlfile_statcov, "independent_variables:\n");
          fprintf(yamlfile_statcov, "- header: {name: 'Bin'}\n");
          fprintf(yamlfile_statcov, "  values:\n");
          for (int i_row = 0; i_row < htotalstatcovmat_AllVar->GetNbinsX(); ++i_row) for (unsigned int i_var = 0; i_var < fVariableNames.size(); ++i_var) for (int i_bin = 0; i_bin < nbinsrhoi; ++i_bin) {
            fprintf(yamlfile_statcov, "  - {value: '$%s$ bin %i'}\n", (MatrixUnf::ObsName(fVariableNames[i_var],kTRUE)).Data(), i_bin+1);
          }
          fprintf(yamlfile_statcov, "- header: {name: 'Bin'}\n");
          fprintf(yamlfile_statcov, "  values:\n");
          for (unsigned int i_var = 0; i_var < fVariableNames.size(); ++i_var) for (int i_bin = 0; i_bin < nbinsrhoi; ++i_bin) for (int i_row = 0; i_row < htotalstatcovmat_AllVar->GetNbinsX(); ++i_row) {
            fprintf(yamlfile_statcov, "  - {value: '$%s$ bin %i'}\n", (MatrixUnf::ObsName(fVariableNames[i_var],kTRUE)).Data(), i_bin+1);
          }
          fclose(yamlfile_statcov);


          for(unsigned int i=0;i<fVariableNames.size();i++){

            FILE *yamlfile = fopen(Form("yaml/%sdiffxsec_%s.yaml",(inorm?"norm":""),fVariableNames[i].Data()), "w");
            fprintf(yamlfile, "dependent_variables:\n");
            if(inorm) fprintf(yamlfile, "- header: {name: '$\\frac{1}{\\sigma} \\frac{{d}\\sigma}{{d}%s}$'}\n", (MatrixUnf::ObsName(fVariableNames[i],kTRUE)).Data() );
            else fprintf(yamlfile, "- header: {name: '$\\frac{{d}\\sigma}{{d}%s}$ [pb]'}\n", (MatrixUnf::ObsName(fVariableNames[i],kTRUE)).Data() );
            fprintf(yamlfile, "  qualifiers:\n");
            fprintf(yamlfile, "  - {name: RE, value: P P --> TOP TOPBAR X}\n");
            fprintf(yamlfile, "  - {name: SQRT(S), units: GEV, value: '13000'}\n");
            fprintf(yamlfile, "  values:\n");
            for (int i_bin = 0; i_bin < nbinsrhoi; ++i_bin)
            {
              fprintf(yamlfile, "  - errors:\n");
              fprintf(yamlfile, "    - {label: stat, symerror: %2.6f}\n", sqrt(htotalstatcovmat_AllVar->GetBinContent(i*nbinsrhoi+i_bin+1,i*nbinsrhoi+i_bin+1)/(divBinWidth?(hbinwidth_AllVar->GetBinContent(i*nbinsrhoi+i_bin+1)*hbinwidth_AllVar->GetBinContent(i*nbinsrhoi+i_bin+1)):1.)) );
              fprintf(yamlfile, "    - {label: sys, symerror: %2.6f}\n", sqrt(htotalsystcovmat_AllVar[inorm]->GetBinContent(i*nbinsrhoi+i_bin+1,i*nbinsrhoi+i_bin+1)/(divBinWidth?(hbinwidth_AllVar->GetBinContent(i*nbinsrhoi+i_bin+1)*hbinwidth_AllVar->GetBinContent(i*nbinsrhoi+i_bin+1)):1.)) );
              fprintf(yamlfile, "    value: %2.6f\n", hnom_AllVar_copy[i_tag*2+inorm]->GetBinContent(i*nbinsrhoi+i_bin+1)/(divBinWidth?(hbinwidth_AllVar->GetBinContent(i*nbinsrhoi+i_bin+1)):1.) );
            }
            fprintf(yamlfile, "independent_variables:\n");
            fprintf(yamlfile, "- header: {name: '$%s$'}\n", (MatrixUnf::ObsName(fVariableNames[i],kTRUE)).Data() );
            fprintf(yamlfile, "  values:\n");
            for (int i_bin = 0; i_bin < nbinsrhoi; ++i_bin)
            {
              fprintf(yamlfile, "  - {high: %2.6f, low: %2.6f}\n", hbinlowedge_AllVar->GetBinContent(i*nbinsrhoi+i_bin+1)+hbinwidth_AllVar->GetBinContent(i*nbinsrhoi+i_bin+1), hbinlowedge_AllVar->GetBinContent(i*nbinsrhoi+i_bin+1) );
            }
            fclose(yamlfile);
          }
        }


        h_RhoCCStat->SetTitle("");
        h_RhoCCStat->Draw("colz");
        cAllVars->Update();
        h_RhoCCStat->GetXaxis()->SetTitle("Coefficient");
        h_RhoCCStat->GetXaxis()->SetTitleOffset(2.07);
        h_RhoCCStat->GetXaxis()->SetTitleSize(0.03);
        h_RhoCCStat->GetXaxis()->SetLabelSize(0.03);
        h_RhoCCStat->GetYaxis()->SetTitle("Coefficient");
        h_RhoCCStat->GetYaxis()->SetTitleOffset(2.02);
        h_RhoCCStat->GetYaxis()->SetTitleSize(0.03);
        h_RhoCCStat->GetYaxis()->SetLabelSize(0.03);
        h_RhoCCStat->GetZaxis()->SetTitleSize(0.03);
        h_RhoCCStat->GetZaxis()->SetLabelSize(0.03);
        h_RhoCCStat->GetZaxis()->SetRangeUser(-1,1);
        for (unsigned int i_bin=1;i_bin<=CoefNames2.size();i_bin++) {
          h_RhoCCStat->GetXaxis()->SetBinLabel(i_bin,CoefNames2.at(i_bin-1).Data());
          h_RhoCCStat->GetYaxis()->SetBinLabel(i_bin,CoefNames2.at(i_bin-1).Data());
        }
        h_RhoCCStat->GetXaxis()->LabelsOption("v");
        //h_RhoCCStat->GetYaxis()->LabelsOption("v");
        DrawCMSLabels(0.035, fLumi);
        if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);
        cAllVars->SaveAs(Form("%s/TotalStatCorrMatrix%s%s_Coefficient_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalStatCorrMatrix%s%s_Coefficient_AllVars.root",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalStatCorrMatrix%s%s_Coefficient_AllVars.C",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));

        h_RhoBSBSStat->SetTitle("");
        h_RhoBSBSStat->Draw("colz");
        cAllVars->Update();
        h_RhoBSBSStat->GetXaxis()->SetTitle("Bin Sum");
        h_RhoBSBSStat->GetXaxis()->SetTitleOffset(2.07);
        h_RhoBSBSStat->GetXaxis()->SetTitleSize(0.03);
        h_RhoBSBSStat->GetXaxis()->SetLabelSize(0.03);
        h_RhoBSBSStat->GetYaxis()->SetTitle("Bin Sum");
        h_RhoBSBSStat->GetYaxis()->SetTitleOffset(2.02);
        h_RhoBSBSStat->GetYaxis()->SetTitleSize(0.03);
        h_RhoBSBSStat->GetYaxis()->SetLabelSize(0.03);
        h_RhoBSBSStat->GetZaxis()->SetTitleSize(0.03);
        h_RhoBSBSStat->GetZaxis()->SetLabelSize(0.03);
        //h_RhoBSBSStat->GetZaxis()->SetRangeUser(-1,1);
        for (unsigned int i_bin=1;i_bin<=CoefNames2.size();i_bin++) {
          h_RhoBSBSStat->GetXaxis()->SetBinLabel(i_bin,CoefNames2.at(i_bin-1).Data());
          h_RhoBSBSStat->GetYaxis()->SetBinLabel(i_bin,CoefNames2.at(i_bin-1).Data());
        }
        h_RhoBSBSStat->GetXaxis()->LabelsOption("v");
        //h_RhoBSBSStat->GetYaxis()->LabelsOption("v");
        DrawCMSLabels(0.04, fLumi);
        if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);
        cAllVars->SaveAs(Form("%s/TotalStatCorrMatrix%s%s_BinSum_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalStatCorrMatrix%s%s_BinSum_AllVars.root",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalStatCorrMatrix%s%s_BinSum_AllVars.C",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));


        h_RhoCCStatSyst->SetTitle("");
        h_RhoCCStatSyst->Draw("colz");
        cAllVars->Update();
        h_RhoCCStatSyst->GetXaxis()->SetTitle("Coefficient");
        h_RhoCCStatSyst->GetXaxis()->SetTitleOffset(2.07);
        h_RhoCCStatSyst->GetXaxis()->SetTitleSize(0.03);
        h_RhoCCStatSyst->GetXaxis()->SetLabelSize(0.03);
        h_RhoCCStatSyst->GetYaxis()->SetTitle("Coefficient");
        h_RhoCCStatSyst->GetYaxis()->SetTitleOffset(2.02);
        h_RhoCCStatSyst->GetYaxis()->SetTitleSize(0.03);
        h_RhoCCStatSyst->GetYaxis()->SetLabelSize(0.03);
        h_RhoCCStatSyst->GetZaxis()->SetTitleSize(0.03);
        h_RhoCCStatSyst->GetZaxis()->SetLabelSize(0.03);
        h_RhoCCStatSyst->GetZaxis()->SetRangeUser(-1,1);
        for (unsigned int i_bin=1;i_bin<=CoefNames2.size();i_bin++) {
          h_RhoCCStatSyst->GetXaxis()->SetBinLabel(i_bin,CoefNames2.at(i_bin-1).Data());
          h_RhoCCStatSyst->GetYaxis()->SetBinLabel(i_bin,CoefNames2.at(i_bin-1).Data());
        }
        h_RhoCCStatSyst->GetXaxis()->LabelsOption("v");
        //h_RhoCCStatSyst->GetYaxis()->LabelsOption("v");
        DrawCMSLabels(0.04, fLumi);
        if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);
        cAllVars->SaveAs(Form("%s/TotalStatSystCorrMatrix%s%s_Coefficient_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalStatSystCorrMatrix%s%s_Coefficient_AllVars.root",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalStatSystCorrMatrix%s%s_Coefficient_AllVars.C",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));

        h_RhoBSBSStatSyst->SetTitle("");
        h_RhoBSBSStatSyst->Draw("colz");
        cAllVars->Update();
        h_RhoBSBSStatSyst->GetXaxis()->SetTitle("Bin Sum");
        h_RhoBSBSStatSyst->GetXaxis()->SetTitleOffset(2.07);
        h_RhoBSBSStatSyst->GetXaxis()->SetTitleSize(0.03);
        h_RhoBSBSStatSyst->GetXaxis()->SetLabelSize(0.03);
        h_RhoBSBSStatSyst->GetYaxis()->SetTitle("Bin Sum");
        h_RhoBSBSStatSyst->GetYaxis()->SetTitleOffset(2.02);
        h_RhoBSBSStatSyst->GetYaxis()->SetTitleSize(0.03);
        h_RhoBSBSStatSyst->GetYaxis()->SetLabelSize(0.03);
        h_RhoBSBSStatSyst->GetZaxis()->SetTitleSize(0.03);
        h_RhoBSBSStatSyst->GetZaxis()->SetLabelSize(0.02);
        //h_RhoBSBSStatSyst->GetZaxis()->SetRangeUser(-1,1);
        for (unsigned int i_bin=1;i_bin<=CoefNames2.size();i_bin++) {
          h_RhoBSBSStatSyst->GetXaxis()->SetBinLabel(i_bin,CoefNames2.at(i_bin-1).Data());
          h_RhoBSBSStatSyst->GetYaxis()->SetBinLabel(i_bin,CoefNames2.at(i_bin-1).Data());
        }
        h_RhoBSBSStatSyst->GetXaxis()->LabelsOption("v");
        //h_RhoBSBSStatSyst->GetYaxis()->LabelsOption("v");
        DrawCMSLabels(0.04, fLumi);
        if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);
        cAllVars->SaveAs(Form("%s/TotalStatSystCorrMatrix%s%s_BinSum_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalStatSystCorrMatrix%s%s_BinSum_AllVars.root",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalStatSystCorrMatrix%s%s_BinSum_AllVars.C",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));


        invtotalstatsystBScormat.Draw("COLZ");
        cAllVars->SaveAs(Form("%s/InvStatSystCorrMatrix%s%s_BinSum_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));

        invtotalstatBScormat.Draw("COLZ");
        cAllVars->SaveAs(Form("%s/InvStatCorrMatrix%s%s_BinSum_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));

        invtotalsystBScormat.Draw("COLZ");
        cAllVars->SaveAs(Form("%s/InvSystCorrMatrix%s%s_BinSum_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));


        TH2D* htotalsystcovmatNorm_AllVar = (TH2D*) htotalsystcovmat_AllVar[inorm]->Clone(htotalsystcovmat_AllVar[inorm]->GetName()+TString("NormCalculated"));
        MatrixUnf::GetNormalisedCovarianceMatrixAllVars(hnom_AllVar_copy[i_tag*2+inorm], htotalsystcovmatNorm_AllVar, htotalsystcovmatNorm_AllVar);

        TH2D* htotalstatsystcovmatCoefficient_AllVar = MatrixUnf::TMatrixDtoTH2(TotalCovMatsCoef_allVars.at(i_tag*2+inorm),Form("TotalCovMatrix_Coefficients_AllVar%s%s",(inorm?"Norm":""),NameTag.Data()));
        TH2D* htotalstatcovmatCoefficient_AllVar = MatrixUnf::TMatrixDtoTH2(StatCovMatsCoef_allVars.at(i_tag*2+inorm),Form("StatCovMatrix_Coefficients_AllVar%s%s",(inorm?"Norm":""),NameTag.Data()));
        TH2D* htotalsystcovmatCoefficient_AllVar = MatrixUnf::TMatrixDtoTH2(SystCovMatsCoef_allVars.at(i_tag*2+inorm),Form("SystCovMatrix_Coefficients_AllVar%s%s",(inorm?"Norm":""),NameTag.Data()));

        TH1D* hCoefficients_AllVar = new TH1D(Form("Coefficients_%s%s",(inorm?"Norm":""),NameTag.Data()), Form("Coefficients_%s%s",(inorm?"Norm":""),NameTag.Data()), fVariableNames.size(),0,fVariableNames.size());
        for (unsigned int i_bin=0;i_bin<fVariableNames.size();i_bin++) {
          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            hCoefficients_AllVar->SetBinContent(i_bin * numCoefficients + ithCoefficient + 1, coefTable[i_tag * 2 + inorm][i_bin][ithCoefficient]);
          }
        }

        TH2D* htotalstatsystcovmatBinSum_AllVar = MatrixUnf::TMatrixDtoTH2(TotalCovMatsBinSum_allVars.at(i_tag*2+inorm),Form("TotalCovMatrix_BinSums_AllVar%s%s",(inorm?"Norm":""),NameTag.Data()));
        TH2D* htotalstatcovmatBinSum_AllVar = MatrixUnf::TMatrixDtoTH2(StatCovMatsBinSum_allVars.at(i_tag*2+inorm),Form("StatCovMatrix_BinSums_AllVar%s%s",(inorm?"Norm":""),NameTag.Data()));
        TH2D* htotalsystcovmatBinSum_AllVar = MatrixUnf::TMatrixDtoTH2(SystCovMatsBinSum_allVars.at(i_tag*2+inorm),Form("SystCovMatrix_BinSums_AllVar%s%s",(inorm?"Norm":""),NameTag.Data()));

        TH1D* hBinSums_AllVar = new TH1D(Form("BinSums_%s%s",(inorm?"Norm":""),NameTag.Data()), Form("BinSums_%s%s",(inorm?"Norm":""),NameTag.Data()), fVariableNames.size(),0,fVariableNames.size());
        for (unsigned int i_bin=0;i_bin<fVariableNames.size();i_bin++) {
          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            hBinSums_AllVar->SetBinContent(i_bin + 1, binsumTable[i_tag * 2 + inorm][i_bin][ithCoefficient]);
          }
        }


        systoutfile_AllVars->cd();
        htotalsystcovmat_AllVar[inorm]->Write();
        // h_RhoCCSyst->Write();
        if(!inorm) htotalsystcovmatNorm_AllVar->Write();
        htotalstatcovmat_AllVar->Write();
        htotalstatsystcovmat_AllVar->Write();
        // h_RhoCCStat->Write();
        hnom_AllVar_copy[i_tag*2+inorm]->Write();

        if(i_tag==0) hpowheg_AllVar[inorm]->Write();
        if(i_tag==0) hmcatnlo_AllVar[inorm]->Write();

        htotalsystcovmatCoefficient_AllVar->Write();
        htotalstatcovmatCoefficient_AllVar->Write();
        htotalstatsystcovmatCoefficient_AllVar->Write();
        hCoefficients_AllVar->Write();


        htotalsystcovmatBinSum_AllVar->Write();
        htotalstatcovmatBinSum_AllVar->Write();
        htotalstatsystcovmatBinSum_AllVar->Write();
        hBinSums_AllVar->Write();


        if(i_tag==0) {

          //write yaml tables for rebinnedA coefficient covariance matrices

          FILE *yamlfile_systcov_coef = fopen(Form("yaml/systcov_%scoef_AllVars.yaml",(inorm?"norm":"")), "w");
          fprintf(yamlfile_systcov_coef, "dependent_variables:\n");
          if(inorm) fprintf(yamlfile_systcov_coef, "- header: {name: 'Systematic covariance for all coefficients'}\n");
          else fprintf(yamlfile_systcov_coef, "- header: {name: 'Systematic covariance for all coefficients'}\n");
          fprintf(yamlfile_systcov_coef, "  qualifiers:\n");
          fprintf(yamlfile_systcov_coef, "  - {name: RE, value: P P --> TOP TOPBAR X}\n");
          fprintf(yamlfile_systcov_coef, "  - {name: SQRT(S), units: GEV, value: '13000'}\n");
          fprintf(yamlfile_systcov_coef, "  values:\n");
          for (int j_bin = 0; j_bin < htotalsystcovmatCoefficient_AllVar->GetNbinsX(); ++j_bin) for (int i_bin = 0; i_bin < htotalsystcovmatCoefficient_AllVar->GetNbinsX(); ++i_bin)
          {
            //fprintf(yamlfile_systcov_coef, "  - errors: []\n");
            fprintf(yamlfile_systcov_coef, "  - {value: %2.10g}\n", htotalsystcovmatCoefficient_AllVar->GetBinContent(i_bin+1,j_bin+1) );
          }
          fprintf(yamlfile_systcov_coef, "independent_variables:\n");
          fprintf(yamlfile_systcov_coef, "- header: {name: 'Coefficient'}\n");
          fprintf(yamlfile_systcov_coef, "  values:\n");
          for (int i_row = 0; i_row < htotalsystcovmatCoefficient_AllVar->GetNbinsX(); ++i_row) for (unsigned int i_var = 0; i_var < fVariableNames.size(); ++i_var) {
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              fprintf(yamlfile_systcov_coef, "  - {value: '$%s$'}\n", (CoefNames.at(i_var * numCoefficients + ithCoefficient)).Data());
            }
          }
          fprintf(yamlfile_systcov_coef, "- header: {name: 'Coefficient'}\n");
          fprintf(yamlfile_systcov_coef, "  values:\n");
          for (unsigned int i_var = 0; i_var < fVariableNames.size(); ++i_var) for (int i_row = 0; i_row < htotalsystcovmatCoefficient_AllVar->GetNbinsX(); ++i_row) {
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              fprintf(yamlfile_systcov_coef, "  - {value: '$%s$'}\n", (CoefNames.at(i_var * numCoefficients + ithCoefficient)).Data());
            }
          }
          fclose(yamlfile_systcov_coef);


          FILE *yamlfile_statcov_coef = fopen(Form("yaml/statcov_%scoef_AllVars.yaml",(inorm?"norm":"")), "w");
          fprintf(yamlfile_statcov_coef, "dependent_variables:\n");
          if(inorm) fprintf(yamlfile_statcov_coef, "- header: {name: 'Statistical covariance for all coefficients'}\n");
          else fprintf(yamlfile_statcov_coef, "- header: {name: 'Statistical covariance for all coefficients'}\n");
          fprintf(yamlfile_statcov_coef, "  qualifiers:\n");
          fprintf(yamlfile_statcov_coef, "  - {name: RE, value: P P --> TOP TOPBAR X}\n");
          fprintf(yamlfile_statcov_coef, "  - {name: SQRT(S), units: GEV, value: '13000'}\n");
          fprintf(yamlfile_statcov_coef, "  values:\n");
          for (int j_bin = 0; j_bin < htotalstatcovmatCoefficient_AllVar->GetNbinsX(); ++j_bin) for (int i_bin = 0; i_bin < htotalstatcovmatCoefficient_AllVar->GetNbinsX(); ++i_bin)
          {
            //fprintf(yamlfile_statcov_coef, "  - errors: []\n");
            fprintf(yamlfile_statcov_coef, "  - {value: %2.10g}\n", htotalstatcovmatCoefficient_AllVar->GetBinContent(i_bin+1,j_bin+1) );
          }
          fprintf(yamlfile_statcov_coef, "independent_variables:\n");
          fprintf(yamlfile_statcov_coef, "- header: {name: 'Coefficient'}\n");
          fprintf(yamlfile_statcov_coef, "  values:\n");
          for (int i_row = 0; i_row < htotalstatcovmatCoefficient_AllVar->GetNbinsX(); ++i_row) for (unsigned int i_var = 0; i_var < fVariableNames.size(); ++i_var) {
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              fprintf(yamlfile_statcov_coef, "  - {value: '$%s$'}\n", (CoefNames.at(i_var * numCoefficients + ithCoefficient)).Data());
            }
          }
          fprintf(yamlfile_statcov_coef, "- header: {name: 'Coefficient'}\n");
          fprintf(yamlfile_statcov_coef, "  values:\n");
          for (unsigned int i_var = 0; i_var < fVariableNames.size(); ++i_var) for (int i_row = 0; i_row < htotalstatcovmatCoefficient_AllVar->GetNbinsX(); ++i_row) {
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              fprintf(yamlfile_statcov_coef, "  - {value: '$%s$'}\n", (CoefNames.at(i_var * numCoefficients + ithCoefficient)).Data());
            }
          }
          fclose(yamlfile_statcov_coef);
        }

        TMatrixD totalstatcorrmat_AllVar = MatrixUnf::TH2toTMatrixD(htotalstatcovmat_AllVar);
        totalstatcorrmat_AllVar.NormByDiag(TMatrixDDiag(totalstatcorrmat_AllVar));
        TH2D* htotalstatcorrmat_AllVar = MatrixUnf::TMatrixDtoTH2(totalstatcorrmat_AllVar,"temp");
        htotalstatcorrmat_AllVar->SetTitle("");
        htotalstatcorrmat_AllVar->Draw("colz");
        cAllVars->Update();
        htotalstatcorrmat_AllVar->GetXaxis()->SetTitle("Bin of normalized differential cross section");
        htotalstatcorrmat_AllVar->GetXaxis()->SetTitleOffset(2.07);
        htotalstatcorrmat_AllVar->GetXaxis()->SetTickLength(0.);
        htotalstatcorrmat_AllVar->GetXaxis()->SetTitleSize(0.03);
        htotalstatcorrmat_AllVar->GetXaxis()->SetLabelSize(0.03);
        htotalstatcorrmat_AllVar->GetYaxis()->SetTitle("Bin of normalized differential cross section");
        htotalstatcorrmat_AllVar->GetYaxis()->SetTitleOffset(2.02);
        htotalstatcorrmat_AllVar->GetYaxis()->SetTickLength(0.);
        htotalstatcorrmat_AllVar->GetYaxis()->SetTitleSize(0.03);
        htotalstatcorrmat_AllVar->GetYaxis()->SetLabelSize(0.03);
        htotalstatcorrmat_AllVar->GetZaxis()->SetTitleSize(0.03);
        htotalstatcorrmat_AllVar->GetZaxis()->SetLabelSize(0.03);
        htotalstatcorrmat_AllVar->GetZaxis()->SetRangeUser(-1,1);
        for (unsigned int i_bin=1;i_bin<=CoefNames2.size()*nbinsrhoi;i_bin++) {
          if (i_bin%nbinsrhoi == 3) {
            htotalstatcorrmat_AllVar->GetXaxis()->SetBinLabel(i_bin,CoefNames2.at((i_bin-1)/nbinsrhoi).Data());
            htotalstatcorrmat_AllVar->GetYaxis()->SetBinLabel(i_bin,CoefNames2.at((i_bin-1)/nbinsrhoi).Data());
          }
          else {
            htotalstatcorrmat_AllVar->GetXaxis()->SetBinLabel(i_bin,"");
            htotalstatcorrmat_AllVar->GetYaxis()->SetBinLabel(i_bin,"");
          }
        }
        labelsx->Draw();
        labelsx2->Draw();
        labelsy->Draw();
        labelsy2->Draw();
        DrawCMSLabels(0.035, fLumi);
        if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);
        cAllVars->SaveAs(Form("%s/TotalStatCorrMatrix%s%s_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalStatCorrMatrix%s%s_AllVars.root",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        cAllVars->SaveAs(Form("%s/TotalStatCorrMatrix%s%s_AllVars.C",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        //delete htotalstatcorrmat_AllVar;

        htotalstatcorrmat_AllVar->Add(hTUnfoldDataStatCorrM_AllVar[i_tag*2+inorm],-1.);

        for (int i_bin = 1; i_bin <= htotalstatcorrmat_AllVar->GetNbinsX(); ++i_bin)
        {
          for (int j_bin = 1; j_bin <= htotalstatcorrmat_AllVar->GetNbinsX(); ++j_bin)
          {
            if( hTUnfoldDataStatCorrM_AllVar[i_tag*2]->GetBinContent(i_bin, j_bin) == 0 ) htotalstatcorrmat_AllVar->SetBinContent(i_bin, j_bin, 0);
          }
        }

        htotalstatcorrmat_AllVar->Draw("colz");
        cAllVars->SaveAs(Form("%s/TotalStatCorrMatrix_BSminusTUnfData%s%s_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));

        delete htotalstatcorrmat_AllVar;

        //TH2D* htotalstatcorrmat_AllVar_clone = (TH2D*) htotalstatcorrmat_AllVar->Clone();
        TH2D* htotalstatcorrmat_AllVar_clone = MatrixUnf::TMatrixDtoTH2(totalstatcorrmat_AllVar,"tempc");

        htotalstatcorrmat_AllVar_clone->SetTitle("");
        htotalstatcorrmat_AllVar_clone->Draw("colz");
        cAllVars->Update();
        htotalstatcorrmat_AllVar_clone->GetXaxis()->SetTitle("Bin of normalized differential cross section");
        htotalstatcorrmat_AllVar_clone->GetXaxis()->SetTitleOffset(2.07);
        htotalstatcorrmat_AllVar_clone->GetXaxis()->SetTickLength(0.);
        htotalstatcorrmat_AllVar_clone->GetXaxis()->SetTitleSize(0.03);
        htotalstatcorrmat_AllVar_clone->GetXaxis()->SetLabelSize(0.03);
        htotalstatcorrmat_AllVar_clone->GetYaxis()->SetTitle("Bin of normalized differential cross section");
        htotalstatcorrmat_AllVar_clone->GetYaxis()->SetTitleOffset(2.02);
        htotalstatcorrmat_AllVar_clone->GetYaxis()->SetTickLength(0.);
        htotalstatcorrmat_AllVar_clone->GetYaxis()->SetTitleSize(0.03);
        htotalstatcorrmat_AllVar_clone->GetYaxis()->SetLabelSize(0.03);
        htotalstatcorrmat_AllVar_clone->GetZaxis()->SetTitleSize(0.03);
        htotalstatcorrmat_AllVar_clone->GetZaxis()->SetLabelSize(0.03);
        htotalstatcorrmat_AllVar_clone->GetZaxis()->SetRangeUser(-1,1);
        for (unsigned int i_bin=1;i_bin<=CoefNames2.size()*nbinsrhoi;i_bin++) {
          if (i_bin%nbinsrhoi == 3) {
            htotalstatcorrmat_AllVar_clone->GetXaxis()->SetBinLabel(i_bin,CoefNames2.at((i_bin-1)/nbinsrhoi).Data());
            htotalstatcorrmat_AllVar_clone->GetYaxis()->SetBinLabel(i_bin,CoefNames2.at((i_bin-1)/nbinsrhoi).Data());
          }
          else {
            htotalstatcorrmat_AllVar_clone->GetXaxis()->SetBinLabel(i_bin,"");
            htotalstatcorrmat_AllVar_clone->GetYaxis()->SetBinLabel(i_bin,"");
          }
        }
        labelsx->Draw();
        labelsx2->Draw();
        labelsy->Draw();
        labelsy2->Draw();
        DrawCMSLabels(0.04, fLumi);
        if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);

        htotalstatcorrmat_AllVar_clone->Add(hTUnfoldTotalStatCorrM_AllVar[i_tag*2+inorm],-1.);

        for (int i_bin = 1; i_bin <= htotalstatcorrmat_AllVar_clone->GetNbinsX(); ++i_bin)
        {
          for (int j_bin = 1; j_bin <= htotalstatcorrmat_AllVar_clone->GetNbinsX(); ++j_bin)
          {
            if( hTUnfoldTotalStatCorrM_AllVar[i_tag*2]->GetBinContent(i_bin, j_bin) == 0 ) htotalstatcorrmat_AllVar_clone->SetBinContent(i_bin, j_bin, 0);
          }
        }

        htotalstatcorrmat_AllVar_clone->Draw("colz");
        cAllVars->SaveAs(Form("%s/TotalStatCorrMatrix_BSminusTUnfTotal%s%s_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));

        delete htotalstatcorrmat_AllVar_clone;

        /*
        htemptotalstatcorrmat_AllVar->SetTitle("");
        htemptotalstatcorrmat_AllVar->Draw("colz");
        cAllVars->Update();
        htemptotalstatcorrmat_AllVar->GetXaxis()->SetTitle("Bin of normalized differential cross section");
        htemptotalstatcorrmat_AllVar->GetXaxis()->SetTitleOffset(2.07);
        htemptotalstatcorrmat_AllVar->GetXaxis()->SetTickLength(0.);
        htemptotalstatcorrmat_AllVar->GetXaxis()->SetTitleSize(0.03);
        htemptotalstatcorrmat_AllVar->GetXaxis()->SetLabelSize(0.03);
        htemptotalstatcorrmat_AllVar->GetYaxis()->SetTitle("Bin of normalized differential cross section");
        htemptotalstatcorrmat_AllVar->GetYaxis()->SetTitleOffset(2.02);
        htemptotalstatcorrmat_AllVar->GetYaxis()->SetTickLength(0.);
        htemptotalstatcorrmat_AllVar->GetYaxis()->SetTitleSize(0.03);
        htemptotalstatcorrmat_AllVar->GetYaxis()->SetLabelSize(0.03);
        htemptotalstatcorrmat_AllVar->GetZaxis()->SetTitleSize(0.03);
        htemptotalstatcorrmat_AllVar->GetZaxis()->SetLabelSize(0.03);
        htemptotalstatcorrmat_AllVar->GetZaxis()->SetRangeUser(-1,1);
        for (unsigned int i_bin=1;i_bin<=CoefNames2.size()*nbinsrhoi;i_bin++) {
          if (i_bin%nbinsrhoi == 3) {
            htemptotalstatcorrmat_AllVar->GetXaxis()->SetBinLabel(i_bin,CoefNames2.at((i_bin-1)/nbinsrhoi).Data());
            htemptotalstatcorrmat_AllVar->GetYaxis()->SetBinLabel(i_bin,CoefNames2.at((i_bin-1)/nbinsrhoi).Data());
          }
          else {
            htemptotalstatcorrmat_AllVar->GetXaxis()->SetBinLabel(i_bin,"");
            htemptotalstatcorrmat_AllVar->GetYaxis()->SetBinLabel(i_bin,"");
          }
        }
        labelsx->Draw();
        labelsx2->Draw();
        labelsy->Draw();
        labelsy2->Draw();
        DrawCMSLabels(0.04);
        if(CHAN!="combined") DrawDecayChLabel(CHAN, 0.04);
        cAllVars->SaveAs(Form("%s/tempTotalStatCorrMatrix%s%s_AllVars.pdf",PLOTFILE.Data(),(inorm?"Norm":""),NameTag.Data()));
        delete htemptotalstatcorrmat_AllVar;
        */

        // Combine measurements b1+/-b2 etc. use coef correlation matrices 
        // h_RhoCCStat h_RhoCCSyst combined with total stat and syst uncertainties from the table.
        
        // coefTable[i_tag*2+inorm][iVar], systTable[i_tag*2+inorm][iVar][SystNames.size()-1], systTable[i_tag*2+inorm][iVar][SystNames.size()-2], systTable[i_tag*2+inorm][iVar][SystNames.size()-6]

        double bPM[nBvar * numCoefficients] = {0};
        double bPMStat[nBvar * numCoefficients] = {0};
        double bPMSyst[nBvar * numCoefficients] = {0};
        double bPMTot[nBvar * numCoefficients] = {0};
        double bPMStat2[nBvar * numCoefficients] = {0};
        double bPMSyst2[nBvar * numCoefficients] = {0};

        for (int iVar = 0; iVar < nBvar; ++iVar) {
          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            int oldIndex = 2 * (iVar / 2);
            int index = 2 * (iVar / 2) + ithCoefficient;
            int firstIndex = 2 * (iVar / 2) * numCoefficients + ithCoefficient;
            int secondIndex = (1 + (2 * (iVar / 2))) * numCoefficients + ithCoefficient;
            bPM[index] += coefTable[i_tag * 2 + inorm][iVar][ithCoefficient];
            bPMStat[index] += pow(systTable[i_tag * 2 + inorm][iVar][ithCoefficient][SystNames.size() - 2], 2);
            bPMStat[index] += h_RhoCCStat->GetBinContent(firstIndex + 1, secondIndex + 1) * systTable[i_tag * 2 + inorm][oldIndex][ithCoefficient][SystNames.size() - 2] * systTable[i_tag * 2 + inorm][oldIndex + 1][ithCoefficient][SystNames.size() - 2];
            bPMSyst[index] += pow(systTable[i_tag * 2 + inorm][iVar][ithCoefficient][SystNames.size() - 6], 2);
            bPMSyst[index] += h_RhoCCSyst->GetBinContent(firstIndex + 1, secondIndex + 1) * systTable[i_tag * 2 + inorm][oldIndex][ithCoefficient][SystNames.size() - 6] * systTable[i_tag * 2 + inorm][oldIndex + 1][ithCoefficient][SystNames.size() - 6];

            bPM[index + 1] += (iVar % 2 > 0 ? -1 : 1) * coefTable[i_tag * 2 + inorm][iVar][ithCoefficient];
            bPMStat[index + 1] += pow(systTable[i_tag * 2 + inorm][iVar][ithCoefficient][SystNames.size() - 2], 2);
            bPMStat[index + 1] -= h_RhoCCStat->GetBinContent(firstIndex + 1, secondIndex + 1) * systTable[i_tag * 2 + inorm][oldIndex][ithCoefficient][SystNames.size() - 2] * systTable[i_tag * 2 + inorm][oldIndex + 1][ithCoefficient][SystNames.size() - 2];
            bPMSyst[index + 1] += pow(systTable[i_tag * 2 + inorm][iVar][ithCoefficient][SystNames.size() - 6], 2);
            bPMSyst[index + 1] -= h_RhoCCSyst->GetBinContent(firstIndex + 1, secondIndex + 1) * systTable[i_tag * 2 + inorm][oldIndex][ithCoefficient][SystNames.size() - 6] * systTable[i_tag * 2 + inorm][oldIndex + 1][ithCoefficient][SystNames.size() - 6];

            bPMTot[index] += TotalCovMatsCoef_allVars.at(i_tag * 2 + inorm)(iVar, iVar);
            bPMTot[index] += TotalCovMatsCoef_allVars.at(i_tag * 2 + inorm)(index, index + 1);
            bPMTot[index + 1] += TotalCovMatsCoef_allVars.at(i_tag * 2 + inorm)(iVar, iVar);
            bPMTot[index + 1] -= TotalCovMatsCoef_allVars.at(i_tag * 2 + inorm)(index, index + 1);

            bPMStat2[index] += StatCovMatsCoef_allVars.at(i_tag * 2 + inorm)(iVar, iVar);
            bPMStat2[index] += StatCovMatsCoef_allVars.at(i_tag * 2 + inorm)(index, index + 1);
            bPMStat2[index + 1] += StatCovMatsCoef_allVars.at(i_tag * 2 + inorm)(iVar, iVar);
            bPMStat2[index + 1] -= StatCovMatsCoef_allVars.at(i_tag * 2 + inorm)(index, index + 1);

            bPMSyst2[index] += SystCovMatsCoef_allVars.at(i_tag * 2 + inorm)(iVar, iVar);
            bPMSyst2[index] += SystCovMatsCoef_allVars.at(i_tag * 2 + inorm)(index, index + 1);
            bPMSyst2[index + 1] += SystCovMatsCoef_allVars.at(i_tag * 2 + inorm)(iVar, iVar);
            bPMSyst2[index + 1] -= SystCovMatsCoef_allVars.at(i_tag * 2 + inorm)(index, index + 1);
          }
        }

        for (int iVar = 0; iVar < nBvar; ++iVar) {
          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            int oldIndex = 2 * (iVar / 2);
            int index = iVar * numCoefficients + ithCoefficient;
            printf("%s%s%s%s%s: %f +/- %f +/- %f ( +/- %f ) \n", fVariableNames[oldIndex].Data(), (iVar % 2 > 0 ? "-" : "+"), fVariableNames[oldIndex + 1].Data(), NameTag.Data(), (inorm ? "Norm" : ""), bPM[index], sqrt(bPMStat[index]), sqrt(bPMSyst[index]), sqrt(bPMStat[index] + bPMSyst[index]));
          }
        }

        for (int iVar = 0; iVar < nBvar; ++iVar) {
          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            int oldIndex = 2 * (iVar / 2);
            int index = iVar * numCoefficients + ithCoefficient;
            printf(" %s%s%s%s%s: %f +/- %f +/- %f ( +/- %f ) \n", fVariableNames[oldIndex].Data(), (iVar % 2 > 0 ? "-" : "+"), fVariableNames[oldIndex + 1].Data(), NameTag.Data(), (inorm ? "Norm" : ""), bPM[index], sqrt(bPMStat2[index]), sqrt(bPMSyst2[index]), sqrt(bPMTot[index]));
          }
        }

        // Write yaml table for rebinnedA bPM coefficients
        if(i_tag==0){
          FILE *yamlfile_coef = fopen(Form("yaml/bPMcoefficients%s.yaml",(inorm?"norm":"")), "w");
          fprintf(yamlfile_coef, "dependent_variables:\n");
          if(inorm) fprintf(yamlfile_coef, "- header: {name: 'Coefficient'}\n");
          else fprintf(yamlfile_coef, "- header: {name: 'Coefficient'}\n");
          fprintf(yamlfile_coef, "  qualifiers:\n");
          fprintf(yamlfile_coef, "  - {name: RE, value: P P --> TOP TOPBAR X}\n");
          fprintf(yamlfile_coef, "  - {name: SQRT(S), units: GEV, value: '13000'}\n");
          fprintf(yamlfile_coef, "  values:\n");
          for (int iVar = 0; iVar < nBvar; ++iVar)
          {
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              int index = iVar * numCoefficients + ithCoefficient;
              fprintf(yamlfile_coef, "  - errors:\n");
              fprintf(yamlfile_coef, "    - {label: stat, symerror: %2.6f}\n", sqrt(bPMStat[index]));
              fprintf(yamlfile_coef, "    - {label: sys, symerror: %2.6f}\n", sqrt(bPMSyst[index]));
              fprintf(yamlfile_coef, "    value: %2.6f\n", bPM[index]);
            }
          }
          fprintf(yamlfile_coef, "independent_variables:\n");
          fprintf(yamlfile_coef, "- header: {name: 'Coefficient'}\n");
          fprintf(yamlfile_coef, "  values:\n");
          for (int iVar = 0; iVar < nBvar; ++iVar)
          {
            for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
              int index = iVar * numCoefficients + ithCoefficient;                          
              fprintf(yamlfile_coef, "  - {value: '$%s%s%s$'}\n", (CoefNames.at(index)).Data(), (iVar % 2 > 0 ? "-" : "+"), (CoefNames.at(index)).Data());
            }
          }
          fclose(yamlfile_coef);
        }

        // compute D from trace of C matrix
        double Dcalc[numCoefficients];
        double DcalcStat[numCoefficients];
        double DcalcSyst[numCoefficients];
        for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
          Dcalc[ithCoefficient] = 0.0;
          DcalcStat[ithCoefficient] = 0.0;
          DcalcSyst[ithCoefficient] = 0.0;
        }

        for (int iVar = nBvar; iVar < nBvar + 3; ++iVar) {
          for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
            Dcalc[ithCoefficient] += coefTable[i_tag * 2 + inorm][iVar][ithCoefficient];
            DcalcStat[ithCoefficient] += pow(systTable[i_tag * 2 + inorm][iVar][ithCoefficient][SystNames.size() - 2], 2);
            DcalcSyst[ithCoefficient] += pow(systTable[i_tag * 2 + inorm][iVar][ithCoefficient][SystNames.size() - 2], 2);

            if (iVar == nBvar) {
              int firstIndex = iVar * numCoefficients + ithCoefficient;
              int secondIndex = (iVar + 1) * numCoefficients + ithCoefficient;
              int thirdIndex = (iVar + 2) * numCoefficients + ithCoefficient;

              DcalcStat[ithCoefficient] += 2. * h_RhoCCStat->GetBinContent(firstIndex + 1, secondIndex + 1) * systTable[i_tag * 2 + inorm][iVar][ithCoefficient][SystNames.size() - 2] * systTable[i_tag * 2 + inorm][iVar + 1][ithCoefficient][SystNames.size() - 2];
              DcalcStat[ithCoefficient] += 2. * h_RhoCCStat->GetBinContent(firstIndex + 1, thirdIndex + 1) * systTable[i_tag * 2 + inorm][iVar][ithCoefficient][SystNames.size() - 2] * systTable[i_tag * 2 + inorm][iVar + 2][ithCoefficient][SystNames.size() - 2];
              DcalcStat[ithCoefficient] += 2. * h_RhoCCStat->GetBinContent(secondIndex + 1, thirdIndex + 3) * systTable[i_tag * 2 + inorm][iVar + 1][ithCoefficient][SystNames.size() - 2] * systTable[i_tag * 2 + inorm][iVar + 2][ithCoefficient][SystNames.size() - 2];

              DcalcSyst[ithCoefficient] += 2. * h_RhoCCSyst->GetBinContent(firstIndex + 1, secondIndex + 1) * systTable[i_tag * 2 + inorm][iVar][ithCoefficient][SystNames.size() - 2] * systTable[i_tag * 2 + inorm][iVar + 1][ithCoefficient][SystNames.size() - 2];
              DcalcSyst[ithCoefficient] += 2. * h_RhoCCSyst->GetBinContent(firstIndex + 1, thirdIndex + 1) * systTable[i_tag * 2 + inorm][iVar][ithCoefficient][SystNames.size() - 2] * systTable[i_tag * 2 + inorm][iVar + 2][ithCoefficient][SystNames.size() - 2];
              DcalcSyst[ithCoefficient] += 2. * h_RhoCCSyst->GetBinContent(secondIndex + 1, thirdIndex + 1) * systTable[i_tag * 2 + inorm][iVar + 1][ithCoefficient][SystNames.size() - 2] * systTable[i_tag * 2 + inorm][iVar + 2][ithCoefficient][SystNames.size() - 2];
            }
          }
        }

        for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
          printf("%s%s%s: %f +/- %f +/- %f ( +/- %f ) \n", "Dcalc", NameTag.Data(), (inorm ? "Norm" : ""), -Dcalc[ithCoefficient] / 3., sqrt(DcalcStat[ithCoefficient]) / 3., sqrt(DcalcSyst[ithCoefficient]) / 3., sqrt(DcalcStat[ithCoefficient] + DcalcSyst[ithCoefficient]) / 3.);
        }
      } //inorm
    } //i_tag


    if(kTRUE){

      // Calculate chi2s and fit results before and after removing the lab frame observables (and then remove the other variables one at a time just for testing, similar to Afiq's evolution plots)
      // TString chi2label[3] = {"","nodPhill","TopRestOnly"};
      TString chi2label[22]   = {"","nodPhill","TopRestOnly","19","18","17","16","15","14","13","12","11","10","9","8","7","6","5","4","3","2","1"};
      //

      double graph_binsum[2][fVariableNames.size()]     = {0.};
      double graph_binsumunc[2][fVariableNames.size()]  = {0.};
      double graph_binsummean[2][fVariableNames.size()] = {0.};
      double graph_nVars[2][fVariableNames.size()]      = {0.};

      for (int iremove = 0; iremove < int(fVariableNames.size()); ++iremove)
      {

        TString NameTag = "_rebinnedA";

        for (int i_unctype = 0; i_unctype <= 2; ++i_unctype) {

          for (int inorm = 0; inorm <= 1; ++inorm)
          {
            // Amandeep : Hardcoded here
            int nbf = 24 - inorm;
            //int nbf = 6 - inorm;
            int nbt = inorm ? nbinstotal_norm : nbinstotal;

            TMatrixD StatCovMat = inorm ? m_totalstatcovmat_AllVar_5bins.at(1) : StatCovMats_allVars.at(0);
            TMatrixD SystCovMat = inorm ? m_totalsystcovmat_AllVar_5bins.at(1) : SystCovMats_allVars.at(0);
            // TMatrixD CovMat  = inorm ? (i_unctype ? (i_unctype > 1 ? m_totalsystcovmat_AllVar_5bins.at(1) : m_totalstatcovmat_AllVar_5bins.at(1) ) : m_totalstatsystcovmat_AllVar_5bins.at(1)) : (i_unctype ? (i_unctype > 1 ? SystCovMats_allVars.at(0) : StatCovMats_allVars.at(0) ) : TotalCovMats_allVars.at(0));

            TMatrixD TotalCovMat = StatCovMat + SystCovMat;

            TMatrixD CovMat = i_unctype ? (i_unctype > 1 ? SystCovMat : StatCovMat ) : TotalCovMat;

            CovMat.ResizeTo(nbt-iremove*nbf,nbt-iremove*nbf);

            TMatrixD CovMatInv = CovMat;
            CovMatInv.Invert();

            TMatrixD Measvec_allVars_copy   = inorm ? Measvec_norm_allVars : Measvec_allVars;
            TMatrixD Theoryvec_allVars_copy = inorm ? Theoryvec_norm_allVars : Theoryvec_allVars;
            // TMatrixD AMCATNLOvec_allVars_copy       = inorm ? AMCATNLOvec_norm_allVars : AMCATNLOvec_allVars;
            // TMatrixD POWHEGV2HERWIGvec_allVars_copy = inorm ? POWHEGV2HERWIGvec_norm_allVars : POWHEGV2HERWIGvec_allVars;
            // TMatrixD TOP_PTvec_allVars_copy  = inorm ? TOP_PTvec_norm_allVars : TOP_PTvec_allVars;
            TMatrixD NLOWvec_allVars_copy       = inorm ? NLOWvec_norm_allVars : NLOWvec_allVars;
            TMatrixD NLOWuncorrvec_allVars_copy = inorm ? NLOWuncorrvec_norm_allVars : NLOWuncorrvec_allVars;
            TMatrixD NLOW3vec_allVars_copy      = inorm ? NLOW3vec_norm_allVars : NLOW3vec_allVars;
            TMatrixD NLOW4vec_allVars_copy      = inorm ? NLOW4vec_norm_allVars : NLOW4vec_allVars;


            std::vector<TMatrixD> NLOWvecMulti_allVars_copy;
            for (int iM = 0; iM < nFit; ++iM)
            {
              if(inorm) NLOWvecMulti_allVars_copy.push_back( NLOWvecMulti_norm_allVars[iM] );
              else NLOWvecMulti_allVars_copy.push_back( NLOWvecMulti_allVars[iM] );

              (NLOWvecMulti_allVars_copy[iM]).ResizeTo(nbt-iremove*nbf,1);
            }

            TMatrixD NLOWvecMulti_allVars_Multi(nbt,nFit);

            for (int i_bin = 0; i_bin < nbt; ++i_bin)
            {
              for (int iM = 0; iM < nFit; ++iM)
              {
                if(inorm) NLOWvecMulti_allVars_Multi(i_bin,iM) = (NLOWvecMulti_norm_allVars[iM])(i_bin,0) - NLOWvec_allVars_copy(i_bin,0) ;
                else      NLOWvecMulti_allVars_Multi(i_bin,iM) = (NLOWvecMulti_allVars[iM])(i_bin,0)      - NLOWvec_allVars_copy(i_bin,0) ;

              }
            }

            NLOWvecMulti_allVars_Multi.ResizeTo(nbt-iremove*nbf,nFit);

            Measvec_allVars_copy.ResizeTo(nbt-iremove*nbf,1);
            Theoryvec_allVars_copy.ResizeTo(nbt-iremove*nbf,1);
            // AMCATNLOvec_allVars_copy.ResizeTo(nbt-iremove*nbf,1);
            // POWHEGV2HERWIGvec_allVars_copy.ResizeTo(nbt-iremove*nbf,1);
            // TOP_PTvec_allVars_copy.ResizeTo(nbt-iremove*nbf,1);
            NLOWvec_allVars_copy.ResizeTo(nbt-iremove*nbf,1);
            NLOWuncorrvec_allVars_copy.ResizeTo(nbt-iremove*nbf,1);
            NLOW3vec_allVars_copy.ResizeTo(nbt-iremove*nbf,1);
            NLOW4vec_allVars_copy.ResizeTo(nbt-iremove*nbf,1);

            // Global correlation coefficient
            TH1D* hrho = new TH1D(Form("hrho%s%s%s%s",(i_unctype ? (i_unctype > 1 ? "Syst" : "Stat") : "Total" ),(inorm?"Norm":""),NameTag.Data(),chi2label[iremove].Data()), Form("hrho%s%s%s%s",(i_unctype ? (i_unctype > 1 ? "Syst" : "Stat") : "Total" ),(inorm?"Norm":""),NameTag.Data(),chi2label[iremove].Data()), nbt-iremove*nbf,0,nbt-iremove*nbf);
            hrho->SetDirectory(0);

            for (int i_bin=0;i_bin<nbt-iremove*nbf;i_bin++) {
              hrho->SetBinContent(i_bin+1, sqrt(1.-1/(CovMatInv(i_bin,i_bin)*CovMat(i_bin,i_bin))) );
            }

            TMatrixD CovMatCoef = i_unctype ? (i_unctype > 1 ? SystCovMatsCoef_allVars.at(inorm) :  StatCovMatsCoef_allVars.at(inorm) ) : TotalCovMatsCoef_allVars.at(inorm);

            CovMatCoef.ResizeTo(nbt/nbf-iremove,nbt/nbf-iremove);

            TMatrixD CovMatCoefInv = CovMatCoef;
            CovMatCoefInv.Invert();

            TMatrixD Meas_Coefs_copy = inorm ? Meas_Coefs : Meas_Coefs_abs;

            Meas_Coefs_copy.ResizeTo(nbt/nbf-iremove,1);
            Theory_Coefs.ResizeTo(nbt/nbf-iremove,1);
            // AMCATNLO_Coefs.ResizeTo(nbt/nbf-iremove,1);
            // POWHEGV2HERWIG_Coefs.ResizeTo(nbt/nbf-iremove,1);
            // TOP_PT_Coefs.ResizeTo(nbt/nbf-iremove,1);
            NLOW_Coefs.ResizeTo(nbt/nbf-iremove,1);
            NLOWuncorr_Coefs.ResizeTo(nbt/nbf-iremove,1);
            NLOW3_Coefs.ResizeTo(nbt/nbf-iremove,1);
            NLOW4_Coefs.ResizeTo(nbt/nbf-iremove,1);


            TMatrixD CovMatBinSum = i_unctype ? (i_unctype > 1 ? SystCovMatsBinSum_allVars.at(0) :  StatCovMatsBinSum_allVars.at(0) ) : TotalCovMatsBinSum_allVars.at(0);

            CovMatBinSum.ResizeTo(nbt/nbf-iremove,nbt/nbf-iremove);

            TMatrixD CovMatBinSumInv = CovMatBinSum;
            CovMatBinSumInv.Invert();

            Meas_BinSums.ResizeTo(nbt/nbf-iremove,1);

            // Global correlation coefficient
            TH1D* hrhoCoef = new TH1D(Form("hrhoCoef%s%s%s%s",(i_unctype ? (i_unctype > 1 ? "Syst" : "Stat") : "Total" ),(inorm?"Norm":""),NameTag.Data(),chi2label[iremove].Data()), Form("hrhoCoef%s%s%s%s",(i_unctype ? (i_unctype > 1 ? "Syst" : "Stat") : "Total" ),(inorm?"Norm":""),NameTag.Data(),chi2label[iremove].Data()), nbt/nbf-iremove,0,nbt/nbf-iremove);
            hrhoCoef->SetDirectory(0);

            for (int i_bin=0;i_bin<nbt/nbf-iremove;i_bin++) {
              hrhoCoef->SetBinContent(i_bin+1, sqrt(1.-1/(CovMatCoefInv(i_bin,i_bin)*CovMatCoef(i_bin,i_bin))) );
            }

            /*
            CovMatCoef.NormByDiag(TMatrixDDiag(CovMatCoef));
            for (int i_bin=0;i_bin<nbt/nbf-iremove;i_bin++) {
              hrhoCoef->SetBinContent(i_bin+1, sqrt(1.-1/(CovMatCoefInv(i_bin,i_bin))) );
            }
            */

            //check closure
            TMatrixD CovMat_closure = CovMat;
            CovMat_closure*=CovMatInv;


            TCanvas *cAllVarsClosure = new TCanvas(Form("cAllVarsClosure%s%s%s%s",(i_unctype ? (i_unctype > 1 ? "Syst" : "Stat") : "Total" ),(inorm?"Norm":""),NameTag.Data(),chi2label[iremove].Data()), Form("cAllVarsClosure%s%s%s%s",(i_unctype ? (i_unctype > 1 ? "Syst" : "Stat") : "Total" ),(inorm?"Norm":""),NameTag.Data(),chi2label[iremove].Data()));
            hrho->Draw();
            cAllVarsClosure->SaveAs(Form("%s/Rho%s%s%s%s_AllVars.pdf",PLOTFILE.Data(),(i_unctype ? (i_unctype > 1 ? "Syst" : "Stat") : "Total" ),(inorm?"Norm":""),NameTag.Data(),chi2label[iremove].Data()));
            hrhoCoef->Draw();
            cAllVarsClosure->SaveAs(Form("%s/RhoCoef%s%s%s%s_AllVars.pdf",PLOTFILE.Data(),(i_unctype ? (i_unctype > 1 ? "Syst" : "Stat") : "Total" ),(inorm?"Norm":""),NameTag.Data(),chi2label[iremove].Data()));
            CovMat_closure.Draw("colz");
            cAllVarsClosure->SaveAs(Form("%s/ClosureMatrix%s%s%s%s_AllVars.pdf",PLOTFILE.Data(),(i_unctype ? (i_unctype > 1 ? "Syst" : "Stat") : "Total" ),(inorm?"Norm":""),NameTag.Data(),chi2label[iremove].Data()));


            // Calculate chi2s for data vs various predictions
            if(kTRUE) {

              TMatrixD dMeas  = Measvec_allVars_copy-Theoryvec_allVars_copy;
              TMatrixD dMeasT = dMeas;
              dMeasT          = dMeasT.T();
              TMatrixD MChi2  = dMeasT*CovMatInv*dMeas;

              std::cout<<"AllVars"<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<< " AllVarsChi2"<<chi2label[iremove]<<" for powheg = "<<MChi2(0,0)<<std::endl;
              AllVarsChi2values[i_unctype][iremove][0] = MChi2(0,0);


              // TMatrixD dMeasAMCATNLO  = Measvec_allVars_copy-AMCATNLOvec_allVars_copy;
              // TMatrixD dMeasAMCATNLOT = dMeasAMCATNLO;
              // dMeasAMCATNLOT          = dMeasAMCATNLOT.T();
              // TMatrixD MChi2AMCATNLO  = dMeasAMCATNLOT*CovMatInv*dMeasAMCATNLO;

              // std::cout<<"AllVars"<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<< " AllVarsChi2"<<chi2label[iremove]<<" for AMCATNLO = "<<MChi2AMCATNLO(0,0)<<std::endl;
              // AllVarsChi2values[i_unctype][iremove][1] = MChi2AMCATNLO(0,0);

              // TMatrixD dMeasPOWHEGV2HERWIG  = Measvec_allVars_copy-POWHEGV2HERWIGvec_allVars_copy;
              // TMatrixD dMeasPOWHEGV2HERWIGT = dMeasPOWHEGV2HERWIG;
              // dMeasPOWHEGV2HERWIGT          = dMeasPOWHEGV2HERWIGT.T();
              // TMatrixD MChi2POWHEGV2HERWIG  = dMeasPOWHEGV2HERWIGT*CovMatInv*dMeasPOWHEGV2HERWIG;

              // std::cout<<"AllVars"<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<< " AllVarsChi2"<<chi2label[iremove]<<" for POWHEGV2HERWIG = "<<MChi2POWHEGV2HERWIG(0,0)<<std::endl;
              // AllVarsChi2values[i_unctype][iremove][2] = MChi2POWHEGV2HERWIG(0,0);

              // TMatrixD dMeasTOP_PT  = Measvec_allVars_copy-TOP_PTvec_allVars_copy;
              // TMatrixD dMeasTOP_PTT = dMeasTOP_PT;
              // dMeasTOP_PTT          = dMeasTOP_PTT.T();
              // TMatrixD MChi2TOP_PT  = dMeasTOP_PTT*CovMatInv*dMeasTOP_PT;

              // std::cout<<"AllVars"<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<< " AllVarsChi2"<<chi2label[iremove]<<" for TOP_PT = "<<MChi2TOP_PT(0,0)<<std::endl;
              // AllVarsChi2values[i_unctype][iremove][3] = MChi2TOP_PT(0,0);

              TMatrixD dMeasNLOW  = Measvec_allVars_copy-NLOWvec_allVars_copy;
              TMatrixD dMeasNLOWT = dMeasNLOW;
              dMeasNLOWT          = dMeasNLOWT.T();
              TMatrixD MChi2NLOW  = dMeasNLOWT*CovMatInv*dMeasNLOW;

              std::cout<<"AllVars"<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<< " AllVarsChi2"<<chi2label[iremove]<<" for NLOW = "<<MChi2NLOW(0,0)<<std::endl;
              AllVarsChi2values[i_unctype][iremove][4] = MChi2NLOW(0,0);


              TMatrixD dMeasNLOWuncorr  = Measvec_allVars_copy-NLOWuncorrvec_allVars_copy;
              TMatrixD dMeasNLOWuncorrT = dMeasNLOWuncorr;
              dMeasNLOWuncorrT          = dMeasNLOWuncorrT.T();
              TMatrixD MChi2NLOWuncorr  = dMeasNLOWuncorrT*CovMatInv*dMeasNLOWuncorr;

              std::cout<<"AllVars"<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<< " AllVarsChi2"<<chi2label[iremove]<<" for NLOWuncorr = "<<MChi2NLOWuncorr(0,0)<<std::endl;
              AllVarsChi2values[i_unctype][iremove][5] = MChi2NLOWuncorr(0,0);


              TMatrixD dMeasNLOW3  = Measvec_allVars_copy-NLOW3vec_allVars_copy;
              TMatrixD dMeasNLOW3T = dMeasNLOW3;
              dMeasNLOW3T          = dMeasNLOW3T.T();
              TMatrixD MChi2NLOW3  = dMeasNLOW3T*CovMatInv*dMeasNLOW3;

              std::cout<<"AllVars"<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<< " AllVarsChi2"<<chi2label[iremove]<<" for NLOW3 = "<<MChi2NLOW3(0,0)<<std::endl;
              AllVarsChi2values[i_unctype][iremove][6] = MChi2NLOW3(0,0);

              TMatrixD dMeasNLOW4  = Measvec_allVars_copy-NLOW4vec_allVars_copy;
              TMatrixD dMeasNLOW4T = dMeasNLOW4;
              dMeasNLOW4T          = dMeasNLOW4T.T();
              TMatrixD MChi2NLOW4  = dMeasNLOW4T*CovMatInv*dMeasNLOW4;

              std::cout<<"AllVars"<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<< " AllVarsChi2"<<chi2label[iremove]<<" for NLOW4 = "<<MChi2NLOW4(0,0)<<std::endl;
              AllVarsChi2values[i_unctype][iremove][7] = MChi2NLOW4(0,0);
            }

            // Analytic solution for chi2 fit

            // first make each bin proportional to mu_t
            TMatrixD CdMeasvec = Measvec_allVars_copy;
            CdMeasvec         -= NLOWvec_allVars_copy;

            // 1D fits
            for (int iM = 0; iM < nFit; ++iM)
            {
              TMatrixD dcvec  = NLOWvecMulti_allVars_copy[iM];
              dcvec -= NLOWvec_allVars_copy;
              TMatrixD dcvecT = dcvec;
              dcvecT = dcvecT.T();

              TMatrixD temp_Av   = dcvecT*CovMatInv;
              TMatrixD temp_Avy  = temp_Av*CdMeasvec;
              TMatrixD temp_AvyN =  temp_Av*dcvec;

              double temp_solution    = temp_Avy(0,0)/temp_AvyN(0,0);
              double temp_solutionerr = sqrt(1./temp_AvyN(0,0));

              std::cout<<" "<<(inorm?"Norm":"")<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<<" "<<chi2label[iremove]<<" AllVars solution for coupling"<<iM<<": "<<temp_solution<<" +/- "<<temp_solutionerr<<std::endl;

              /*
              //print weight for the measured deviation in each bin (coupling_measured = sum[weight*deviation] )
              std::cout<<"effective weight for the measured deviation in each bin:"<<std::endl;
              temp_Av *= 1./temp_AvyN(0,0);
              TMatrixD AvT = temp_Av;
              AvT.T();
              AvT.Print("f=%1.5g ");

              //print contribution to the variance from each bin
              std::cout<<" contribution to the variance from each bin:"<<std::endl;
              TMatrixD dVar = CovMat*AvT;
              dVar *= TMatrixDColumn(AvT,0);
              dVar.Print("f=%1.5g ");

              //print contribution to the inverse variance from each bin
              std::cout<<" contribution to the inverse variance from each bin:"<<std::endl;
              TMatrixD dIVar = CovMatInv*dcvec;
              dIVar *= TMatrixDColumn(dcvec,0);
              dIVar.Print("f=%1.5g ");
              */

            }
    
            // Amandeep : Where is this used ?
            // multi-dimensional fit

            TMatrixD dmutvecMulti  = NLOWvecMulti_allVars_Multi;
            TMatrixD dmutvecMultiT = dmutvecMulti;
            dmutvecMultiT          = dmutvecMultiT.T();

            TMatrixD AvMulti   = dmutvecMultiT*CovMatInv;
            TMatrixD AvMultiy  = AvMulti*CdMeasvec;
            TMatrixD AvMultiyN = AvMulti*dmutvecMulti;
            AvMultiyN.Invert();
            TMatrixD solutionMulti = AvMultiyN*AvMultiy;

            std::cout<<" "<<(inorm?"Norm":"")<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<<" "<<chi2label[iremove]<<" Multiple coupling fit results:"<<std::endl;
            solutionMulti.Print("f=%1.5g ");
            AvMultiyN.Print("f=%1.5g ");

            bool goodsol = kTRUE;
            // if one or more of the couplings is not constrained by the fitted distributions there will be zeros in the result, so calculation of the correlation matrix (below) will fail
            for (int iM = 0; iM < nFit; ++iM)
            {
              if( solutionMulti(iM,0) == 0 ) goodsol = kFALSE; 
            }

            if(goodsol) {
              std::cout << " Correlation matrix:" << std::endl;
              AvMultiyN.NormByDiag(TMatrixDDiag(AvMultiyN));
              AvMultiyN.Print("f=%1.5g ");
            }
            else std::cout<<" One or more couplings not constrained by input distributions."<<std::endl;

            // Calculate chi2s for data vs various predictions
            if(kTRUE) {

              TMatrixD dMeas  = Meas_Coefs_copy-Theory_Coefs;
              TMatrixD dMeasT = dMeas;
              dMeasT          = dMeasT.T();
              TMatrixD MChi2  = dMeasT*CovMatCoefInv*dMeas;

              std::cout<<"AllCoefs"<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<< " AllCoefsChi2"<<chi2label[iremove]<<" for powheg = "<<MChi2(0,0)<<std::endl;
              CoefChi2values[i_unctype][iremove][0] = MChi2(0,0);


              // TMatrixD dMeasAMCATNLO = Meas_Coefs_copy-AMCATNLO_Coefs;
              // TMatrixD dMeasAMCATNLOT = dMeasAMCATNLO;
              // dMeasAMCATNLOT = dMeasAMCATNLOT.T();
              // TMatrixD MChi2AMCATNLO = dMeasAMCATNLOT*CovMatCoefInv*dMeasAMCATNLO;

              // std::cout<<"AllCoefs"<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<< " AllCoefsChi2"<<chi2label[iremove]<<" for AMCATNLO = "<<MChi2AMCATNLO(0,0)<<std::endl;
              // CoefChi2values[i_unctype][iremove][1] = MChi2AMCATNLO(0,0);

              // TMatrixD dMeasPOWHEGV2HERWIG = Meas_Coefs_copy-POWHEGV2HERWIG_Coefs;
              // TMatrixD dMeasPOWHEGV2HERWIGT = dMeasPOWHEGV2HERWIG;
              // dMeasPOWHEGV2HERWIGT = dMeasPOWHEGV2HERWIGT.T();
              // TMatrixD MChi2POWHEGV2HERWIG = dMeasPOWHEGV2HERWIGT*CovMatCoefInv*dMeasPOWHEGV2HERWIG;

              // std::cout<<"AllCoefs"<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<< " AllCoefsChi2"<<chi2label[iremove]<<" for POWHEGV2HERWIG = "<<MChi2POWHEGV2HERWIG(0,0)<<std::endl;
              // CoefChi2values[i_unctype][iremove][2] = MChi2POWHEGV2HERWIG(0,0);

              // TMatrixD dMeasTOP_PT = Meas_Coefs_copy-TOP_PT_Coefs;
              // TMatrixD dMeasTOP_PTT = dMeasTOP_PT;
              // dMeasTOP_PTT = dMeasTOP_PTT.T();
              // TMatrixD MChi2TOP_PT = dMeasTOP_PTT*CovMatCoefInv*dMeasTOP_PT;

              // std::cout<<"AllCoefs"<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<< " AllCoefsChi2"<<chi2label[iremove]<<" for TOP_PT = "<<MChi2TOP_PT(0,0)<<std::endl;
              // CoefChi2values[i_unctype][iremove][3] = MChi2TOP_PT(0,0);


              TMatrixD dMeasNLOW  = Meas_Coefs_copy-NLOW_Coefs;
              TMatrixD dMeasNLOWT = dMeasNLOW;
              dMeasNLOWT          = dMeasNLOWT.T();
              TMatrixD MChi2NLOW  = dMeasNLOWT*CovMatCoefInv*dMeasNLOW;

              std::cout<<"AllCoefs"<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<< " AllCoefsChi2"<<chi2label[iremove]<<" for NLOW = "<<MChi2NLOW(0,0)<<std::endl;
              CoefChi2values[i_unctype][iremove][4] = MChi2NLOW(0,0);


              TMatrixD dMeasNLOWuncorr  = Meas_Coefs_copy-NLOWuncorr_Coefs;
              TMatrixD dMeasNLOWuncorrT = dMeasNLOWuncorr;
              dMeasNLOWuncorrT          = dMeasNLOWuncorrT.T();
              TMatrixD MChi2NLOWuncorr  = dMeasNLOWuncorrT*CovMatCoefInv*dMeasNLOWuncorr;

              std::cout<<"AllCoefs"<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<< " AllCoefsChi2"<<chi2label[iremove]<<" for NLOWuncorr = "<<MChi2NLOWuncorr(0,0)<<std::endl;
              CoefChi2values[i_unctype][iremove][5] = MChi2NLOWuncorr(0,0);


              TMatrixD dMeasNLOW3  = Meas_Coefs_copy-NLOW3_Coefs;
              TMatrixD dMeasNLOW3T = dMeasNLOW3;
              dMeasNLOW3T          = dMeasNLOW3T.T();
              TMatrixD MChi2NLOW3  = dMeasNLOW3T*CovMatCoefInv*dMeasNLOW3;

              std::cout<<"AllCoefs"<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<< " AllCoefsChi2"<<chi2label[iremove]<<" for NLOW3 = "<<MChi2NLOW3(0,0)<<std::endl;
              CoefChi2values[i_unctype][iremove][6] = MChi2NLOW3(0,0);

              TMatrixD dMeasNLOW4  = Meas_Coefs_copy-NLOW4_Coefs;
              TMatrixD dMeasNLOW4T = dMeasNLOW4;
              dMeasNLOW4T          = dMeasNLOW4T.T();
              TMatrixD MChi2NLOW4  = dMeasNLOW4T*CovMatCoefInv*dMeasNLOW4;

              std::cout<<"AllCoefs"<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<< " AllCoefsChi2"<<chi2label[iremove]<<" for NLOW4 = "<<MChi2NLOW4(0,0)<<std::endl;
              CoefChi2values[i_unctype][iremove][7] = MChi2NLOW4(0,0);

            }

            // For testing only. For absolute differential cross sections, fit "combined" inclusive cross section to show this becomes biased 
            // when combining several strongly correlated measurements (since the absolute xsec is almost the same for each distribution)
            if(inorm==0) {

              double total_binsum = 0;

              for (unsigned int i_vec = 0; i_vec < fVariableNames.size()-iremove; ++i_vec) {
                total_binsum += Meas_BinSums(i_vec,0);
              }

              total_binsum /= double(fVariableNames.size()-iremove);


              const int nChi2 = fVariableNames.size();
              double a_chi2[nChi2];
              double a_BinSum[nChi2];

              // numerical chi2 fit, but could easily be replaced with analytical result
              for (int iChi = 0; iChi < nChi2; ++iChi)
              {

                double deltaSum=double(iChi-(nChi2-1)/2.)*(i_unctype?1:20.);

                TMatrixD Test_BinSums = Meas_BinSums;

                for (unsigned int i_vec = 0; i_vec < fVariableNames.size()-iremove; ++i_vec) {
                  Test_BinSums(i_vec,0) = total_binsum + deltaSum;
                }

                TMatrixD dMeasTest  = Meas_BinSums-Test_BinSums;
                TMatrixD dMeasTestT = dMeasTest;
                dMeasTestT          = dMeasTestT.T();
                TMatrixD MChi2Test  = dMeasTestT*CovMatBinSumInv*dMeasTest;

                //std::cout<<"AllBinSums"<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<< " AllBinSumsChi2"<<chi2label[iremove]<<" for deltaSum= "<<deltaSum<<" : "<<MChi2Test(0,0)<<std::endl;

                a_chi2[iChi]   = MChi2Test(0,0);
                a_BinSum[iChi] = total_binsum + deltaSum;
              }

              TGraph *Chi2_curve = new TGraph(nChi2,a_BinSum,a_chi2);
              Chi2_curve->SetName("Chi2vsBinSum");
              Chi2_curve->GetXaxis()->SetTitle("BinSum");
              Chi2_curve->GetYaxis()->SetTitle("#chi^{2}");

              double param[3]={0.};

              Chi2_curve->Fit("pol2");
              if(Chi2_curve->GetFunction("pol2")) param[0] = Chi2_curve->GetFunction("pol2")->GetParameter(0);
              if(Chi2_curve->GetFunction("pol2")) param[1] = Chi2_curve->GetFunction("pol2")->GetParameter(1);
              if(Chi2_curve->GetFunction("pol2")) param[2] = Chi2_curve->GetFunction("pol2")->GetParameter(2);

              std::cout<<"AllBinSums"<<(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" )<< " Chi2 curve result: "<<param[0]<<" "<<param[1]<<" "<<param[2]<<", "<<-param[1]/(2.*param[2])<<" +/- "<<1./sqrt(param[2])<<std::endl;

              graph_binsum[i_unctype][iremove]=-param[1]/(2.*param[2]);
              graph_binsumunc[i_unctype][iremove]=1./sqrt(param[2]);
              graph_nVars[i_unctype][iremove]=fVariableNames.size()-iremove;

              double binsummean = 0.0;
              for (unsigned int i_bin = 0; i_bin < fVariableNames.size() - iremove; ++i_bin) {
                for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
                  binsummean += binsumTable[0][i_bin][ithCoefficient];
                }
              }
              binsummean /= double((fVariableNames.size() - iremove) * numCoefficients);
              graph_binsummean[i_unctype][iremove] = binsummean;
            }

          }//inorm

        }//i_unctype

      }//iremove


      gStyle->SetPadLeftMargin(0.13);
      gStyle->SetPadRightMargin(0.05);
      gStyle->SetTitleOffset(1.2, "Y");


      for (int i_unctype = 0; i_unctype <= 1; ++i_unctype)
      {

        TCanvas *cBinsumUnc = new TCanvas(Form("cBinsumUnc%s",(i_unctype ? "Stat" : "Total" )), Form("cBinsumUnc%s",(i_unctype ? "Stat" : "Total" )));

        TGraph *UncvsnVars = new TGraph(fVariableNames.size(),graph_nVars[i_unctype],graph_binsumunc[i_unctype]);
        UncvsnVars->SetName("Uncertainty vs nVars");
        UncvsnVars->GetYaxis()->SetTitle("Cross section uncertainty [pb]");
        UncvsnVars->GetXaxis()->SetTitle("Number of fitted variables");

        UncvsnVars->Draw();
        cBinsumUnc->SaveAs(Form("%s/UncvsnVars%s.pdf",PLOTFILE.Data(),(i_unctype ? "Stat" : "Total" ) ));

        TCanvas *cBinsum = new TCanvas(Form("cBinsum%s",(i_unctype ? "Stat" : "Total" )), Form("cBinsum%s",(i_unctype ? "Stat" : "Total" )));

        TGraph *ValuevsnVars = new TGraph(fVariableNames.size(),graph_nVars[i_unctype],graph_binsum[i_unctype]);
        ValuevsnVars->SetName("Measured xsec vs nVars");
        ValuevsnVars->GetYaxis()->SetTitle("Cross section [pb]");
        ValuevsnVars->GetXaxis()->SetTitle("Number of fitted variables");

        ValuevsnVars->Draw();
        cBinsum->SaveAs(Form("%s/ValuevsnVars%s.pdf",PLOTFILE.Data(),(i_unctype ? "Stat" : "Total" ) ));

        TGraph *MeanvsnVars = new TGraph(fVariableNames.size(),graph_nVars[i_unctype],graph_binsummean[i_unctype]);
        MeanvsnVars->SetName("Measured xsec vs nVars");
        MeanvsnVars->GetYaxis()->SetTitle("Cross section [pb]");
        MeanvsnVars->GetXaxis()->SetTitle("Number of fitted variables");

        MeanvsnVars->SetLineColor(2);
        MeanvsnVars->SetMarkerColor(2);
        MeanvsnVars->SetMarkerSize(0);

        MeanvsnVars->Draw();
        cBinsum->SaveAs(Form("%s/ValueandMeanvsnVars%s.pdf",PLOTFILE.Data(),(i_unctype ? "Stat" : "Total" ) ));


      }

      // AllVars chi2 table
      for (int i_unctype = 0; i_unctype <= 2; ++i_unctype)
      {

        printf(" Table of %s $\\chi^2$ values for all bins at once: \n",(i_unctype ? (i_unctype > 1 ? " Syst" : " Stat") : " Total" ));

        printf("\\begin{tabular}{l | c c c c c c c c } \n");
        printf("\\hline \n");
        printf("\\textbf{Distribution} & \\small{\\Powhegvtwo P8} & \\small{\\MGaMCatNLO} & \\small{\\Powhegvtwo H++} & \\small{Top \\pt rw} & \\small{NLOW} & \\small{No SC/Pol.} & \\small{measA} & \\small{measN} \\\\ \n");
        for (unsigned int iVar = 0; iVar < fVariableNames.size(); ++iVar)
        {
          printf("$%s$ & $%2.1f$ & $%2.1f$ & $%2.1f$ & $%2.1f$ & $%2.1f$ & $%2.1f$ & $%2.1f$ & $%2.1f$  \\\\ \n",chi2label[iVar].Data(), AllVarsChi2values[i_unctype][iVar][0], AllVarsChi2values[i_unctype][iVar][1], AllVarsChi2values[i_unctype][iVar][2], AllVarsChi2values[i_unctype][iVar][3], AllVarsChi2values[i_unctype][iVar][4], AllVarsChi2values[i_unctype][iVar][5], AllVarsChi2values[i_unctype][iVar][6], AllVarsChi2values[i_unctype][iVar][7]);
        }
        printf("\\hline \n");
        printf("\\end{tabular} \n \n");

      }

      // Coef chi2 table
      for (int i_unctype = 0; i_unctype <= 1; ++i_unctype)
      {

        printf(" Table of %s $\\chi^2$ values for multiple coefficients fit together: \n",(i_unctype ? " Stat" : " Total" ));

        printf("\\begin{tabular}{l | c c c c c c c c } \n");
        printf("\\hline \n");
        printf("\\textbf{Distribution} & \\small{\\Powhegvtwo P8} & \\small{\\MGaMCatNLO} & \\small{\\Powhegvtwo H++} & \\small{Top \\pt rw} & \\small{NLOW} & \\small{No SC/Pol.} & \\small{measA} & \\small{measN} \\\\ \n");
        for (unsigned int iVar = 0; iVar < fVariableNames.size(); ++iVar)
        {
          printf("$%s$ & $%2.1f$ & $%2.1f$ & $%2.1f$ & $%2.1f$ & $%2.1f$ & $%2.1f$ & $%2.1f$ & $%2.1f$  \\\\ \n",chi2label[iVar].Data(), CoefChi2values[i_unctype][iVar][0], CoefChi2values[i_unctype][iVar][1], CoefChi2values[i_unctype][iVar][2], CoefChi2values[i_unctype][iVar][3], CoefChi2values[i_unctype][iVar][4], CoefChi2values[i_unctype][iVar][5], CoefChi2values[i_unctype][iVar][6], CoefChi2values[i_unctype][iVar][7]);
        }
        printf("\\hline \n");
        printf("\\end{tabular} \n \n");

      }
    }

    printf("\n Writing AllVars systematics output file ...\n");
    systoutfile_AllVars->Write();
    //systoutfile_AllVars->Close();
    delete systoutfile_AllVars;

  }
}

void MatrixUnfControl::calculateCoefficientUncertainty(TH2D *hsystcovmat, TH1D *tunfresult_temp, TH1D* genhistogram, std::vector<double> &a_optimal_reco, TString VariableName, double *binsum_temp, int numCoefficients, double Cerr[]) {
  int numBinsPerCoefficient = 6; // hardcoded
    
  int nbinsrhoi = hsystcovmat->GetNbinsX();

  TMatrixD m_TUnfCovMatSys(nbinsrhoi,nbinsrhoi);
  for (Int_t cmi = 0; cmi < nbinsrhoi; cmi++)
  {
    for (Int_t cmj = 0; cmj < nbinsrhoi; cmj++)
    {
      m_TUnfCovMatSys(cmi, cmj) = hsystcovmat->GetBinContent(cmi + 1, cmj + 1);
    }
  }

  TH1D* hist_acombfactor_single[nbinsrhoi/2];

  for (int i_a = 0; i_a < nbinsrhoi/2; ++i_a)
  {
    hist_acombfactor_single[i_a] = new TH1D(Form("hist_acombfactor_single%d",i_a), Form("hist_acombfactor_single%d",i_a), 1,0,1);
  }

  std::vector<std::vector<double>> afbvecvec;
  std::vector<std::vector<double>> afberrvecvec;
  for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
    std::vector<double> factor{1., 1., 1.};
    MatrixUnf::AtoCfactor(VariableName, genhistogram, factor, numCoefficients, ithCoefficient);
    Cerr[ithCoefficient] = 0.0;

    int binOffset = ithCoefficient * numBinsPerCoefficient;
    int asymmetryOffset = ithCoefficient * numBinsPerCoefficient / 2;

    // fill factors used in combining binned asymmetries into coefficients
    hist_acombfactor_single[asymmetryOffset + 0]->SetBinContent(1, (a_optimal_reco.at(asymmetryOffset + 0)) / factor[0]);
    hist_acombfactor_single[asymmetryOffset + 1]->SetBinContent(1, (a_optimal_reco.at(asymmetryOffset + 1)) / factor[1]);
    hist_acombfactor_single[asymmetryOffset + 2]->SetBinContent(1, (1. - a_optimal_reco.at(asymmetryOffset + 0) - a_optimal_reco.at(asymmetryOffset + 1)) / factor[2]);
  }

  //hist_acombfactor_single[0]->Print("all");
  //hist_acombfactor_single[1]->Print("all");
  //hist_acombfactor_single[2]->Print("all");


  //calculate the sys on coef etc
  /*
  TMatrixD m_AFB(1,1);
  std::vector<double> afbvar;
  std::vector<double> afberrvar;
  MatrixUnf::GetAfbCorrMAllVars(tunfresult_temp, m_TUnfCovMatSys, afbvar, afberrvar, m_AFB, kTRUE);
  //m_AFB.Print("f=%1.5g ");

  TMatrixD m_BinSum(1,1);
  std::vector<double> binsum;
  std::vector<double> binsumerr;
  MatrixUnf::GetBinSumCorrMAllVars(tunfresult_temp, m_TUnfCovMatSys, binsum, binsumerr, m_BinSum, kTRUE);
  //m_BinSum.Print("f=%1.5g ");
  */

  if (MatrixUnf::isOtherSymType(VariableName)) {
    TMatrixD m_AFB(numCoefficients, numCoefficients);
    std::vector<double> afbvar;
    std::vector<double> afberrvar;
    MatrixUnf::GetAfbCorrMAllVars(tunfresult_temp, m_TUnfCovMatSys, afbvar, afberrvar, m_AFB, kTRUE);
    for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
      Cerr[ithCoefficient] = m_AFB(ithCoefficient, ithCoefficient) > 0 ? sqrt(m_AFB(ithCoefficient, ithCoefficient)) : 0.0;
      std::cout << "Cerr at " << ithCoefficient << " is " << Cerr[ithCoefficient] << std::endl;
    }
  } else {
    TMatrixD m_AFBB(nbinsrhoi / 2, nbinsrhoi / 2);
    TMatrixD m_C(numCoefficients, numCoefficients);
    MatrixUnf::GetAfbBinsCorrMAllVars(tunfresult_temp, m_TUnfCovMatSys, afbvecvec, afberrvecvec, m_AFBB, hist_acombfactor_single, m_C, kTRUE, numCoefficients);
    // m_AFBB.Print("f=%1.5g ");
    // m_C.Print("f=%1.5g ");
    for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
      Cerr[ithCoefficient] = m_C(ithCoefficient, ithCoefficient) > 0 ? sqrt(m_C(ithCoefficient, ithCoefficient)) : 0.0;
      std::cout << "Cerr at " << ithCoefficient << " is " << Cerr[ithCoefficient] << std::endl;
    }
  }


  TMatrixD m_BinSum(numCoefficients, numCoefficients);
  std::vector<double> binsum;
  std::vector<double> binsumerr;
  MatrixUnf::GetBinSumCorrMAllVars(tunfresult_temp, m_TUnfCovMatSys, binsum, binsumerr, m_BinSum, kTRUE);
  for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
    binsum_temp[ithCoefficient] = m_BinSum(ithCoefficient, ithCoefficient) > 0 ? sqrt(m_BinSum(ithCoefficient, ithCoefficient)) : 0.0;
  }

  
  //print the correlation matrices to check the elements are all +/-1 both before and after rebinning
  //std::cout<<"corrm: "<<std::endl;
  m_TUnfCovMatSys.NormByDiag(TMatrixDDiag(m_TUnfCovMatSys));
  //m_TUnfCovMatSys.Print("f=%1.5g ");
  

  for (int ithCoefficient = 0; ithCoefficient < numCoefficients; ithCoefficient++) {
    int asymmetryOffset = ithCoefficient * numBinsPerCoefficient / 2;

    delete hist_acombfactor_single[asymmetryOffset + 0];
    delete hist_acombfactor_single[asymmetryOffset + 1];
    delete hist_acombfactor_single[asymmetryOffset + 2];
  }
}


TH2D* MatrixUnfControl::calculateCoefficientSystCorrmAllVars(TString TUNFFILE, TH2D* hsystcovmat, TH1D* hist_acombfactor[], TString NameTag, int inorm) {

  // Amandeep : Hardcoded here
  const int nbinsrhoi = 24;
  //const int nbinsrhoi    = 6;

  const int nbins = hsystcovmat->GetNbinsX();
  const int nVar  = fVariableNames.size();

  TMatrixD m_TUnfCovMatSys(nbins,nbins);

  for (Int_t cmi = 0; cmi < nbins; cmi++) {
    for (Int_t cmj = 0; cmj < nbins; cmj++) {
      m_TUnfCovMatSys(cmi, cmj) = hsystcovmat->GetBinContent(cmi + 1, cmj + 1);
    }
  }

  /*
    TH1D* hist_acombfactor[nbinsrhoi/2];
    for (int i_a = 0; i_a < nbinsrhoi/2; ++i_a) {
      hist_acombfactor[i_a] = new TH1D(Form("hist_acombfactor%d",i_a), Form("hist_acombfactor%d",i_a), nVar,0,nVar);
    }
  */

  TH1D* hnom_AllVar = new TH1D("hnom_AllVar", "hnom_AllVar", fVariableNames.size()*nbinsrhoi,0,fVariableNames.size()*nbinsrhoi);
  int bin_counter = 1;

  for(int i=0; i < nVar; i++)
  {

    TString VariableName = fVariableNames.at(i);

    if(VariableName=="") continue;

    TString TUnfoldResultsFilename = Form(TUNFFILE+VariableName+".root");
    TFile *unffile                 = new TFile(TUnfoldResultsFilename.Data(),"READ");
    if(unffile==NULL) exit(0);

    TH1D* tunfresult_temp = (TH1D*)((TH1D*)unffile->Get(VariableName+"TUnfResultCor"+NameTag))->Clone(VariableName+"TUnfResult_temp");
    if(inorm>0) tunfresult_temp->Scale(1./tunfresult_temp->Integral());
    
    /*
        TH1D* afbresult_temp = (TH1D*)((TH1D*)unffile->Get(VariableName+"AfbResult"+NameTag))->Clone(VariableName+"AfbResult_temp");

        // Fill factors used in combining binned asymmetries into coefficients
        TH1D *hGenMC_rebinnedS = (TH1D*)((TH1D*)unffile->Get(VariableName+"Gen_norm_rebinnedB"))->Clone(Form("hGenMC_rebinnedS_%s",VariableName.Data()));
        std::vector<double> factor {1.,1.,1.};
        MatrixUnf::AtoCfactor(VariableName, hGenMC_rebinnedS, factor);

        hist_acombfactor[0]->SetBinContent(i+1, (afbresult_temp->GetBinContent(0))/factor[0] );
        hist_acombfactor[1]->SetBinContent(i+1, (afbresult_temp->GetBinContent(nbinsrhoi+1))/factor[1] );
        hist_acombfactor[2]->SetBinContent(i+1, (1. - afbresult_temp->GetBinContent(0) - afbresult_temp->GetBinContent(nbinsrhoi+1))/factor[2] );
    */
    
    for(int ibins = 1; ibins<=tunfresult_temp->GetNbinsX(); ibins++){
      hnom_AllVar->SetBinContent(bin_counter,tunfresult_temp->GetBinContent(ibins));
      ++bin_counter;
    }

    delete tunfresult_temp;
    // delete afbresult_temp;
    // delete hGenMC_rebinnedS;

    unffile->Close();
    delete unffile;
  }


  // Calculate the sys on coef etc
  /*
    TMatrixD m_AFB(1,1);
    std::vector<double> afbvar;
    std::vector<double> afberrvar;
    MatrixUnf::GetAfbCorrMAllVars(tunfresult_temp, m_TUnfCovMatSys, afbvar, afberrvar, m_AFB, kTRUE);
    //m_AFB.Print("f=%1.5g ");

    TMatrixD m_BinSum(1,1);
    std::vector<double> binsum;
    std::vector<double> binsumerr;
    MatrixUnf::GetBinSumCorrMAllVars(tunfresult_temp, m_TUnfCovMatSys, binsum, binsumerr, m_BinSum, kTRUE);
    //m_BinSum.Print("f=%1.5g ");
  */

  TMatrixD m_AFBB(nVar*nbinsrhoi/2,nVar*nbinsrhoi/2);
  TMatrixD m_C(nVar,nVar);
  std::vector<std::vector<double>> afbvecvec;
  std::vector<std::vector<double>> afberrvecvec;
  MatrixUnf::GetAfbBinsCorrMAllVars(hnom_AllVar, m_TUnfCovMatSys, afbvecvec, afberrvecvec, m_AFBB, hist_acombfactor, m_C, kTRUE);
  //m_AFBB.Print("f=%1.5g ");
  //m_C.Print("f=%1.5g ");

  TVectorD diag_nozeros = TMatrixDDiag(m_C);
  for(int xbin=1; xbin<=nVar; xbin++)
  {
    if( diag_nozeros(xbin-1) == 0 ) diag_nozeros(xbin-1) = 1e-20;
  }
  m_C.NormByDiag(diag_nozeros);
  //std::cout<<"Coef sys corrm: "<<std::endl;
  //m_C.Print("f=%1.5g ");

  TH2D* h_RhoCC = new TH2D(Form("RhoCC_SystTotal%s%s",(inorm?"Norm":""),NameTag.Data()), Form("RhoCC_SystTotal%s%s",(inorm?"Norm":""),NameTag.Data()), nVar, 0, nVar, nVar, 0, nVar);
  for(int i=0; i < nVar; i++){
    for(int j=0; j < nVar; j++){
      h_RhoCC->SetBinContent( i+1, j+1, m_C(i,j) );
    }
  }

  // delete hist_acombfactor[0];
  // delete hist_acombfactor[1];
  // delete hist_acombfactor[2];

  delete hnom_AllVar;

  return h_RhoCC;

}


TH2D* MatrixUnfControl::calculateBinSumSystCorrmAllVars(TString TUNFFILE, TH2D* hsystcovmat, TH1D* hist_acombfactor[], TString NameTag, int inorm) {
  // Amandeep : Hardcoded here
  const int nbinsrhoi = 24;
  //const int nbinsrhoi    = 6;
  
  const int nbins = hsystcovmat->GetNbinsX();
  const int nVar  = fVariableNames.size();

  TMatrixD m_TUnfCovMatSys(nbins,nbins);

  for (Int_t cmi = 0; cmi < nbins; cmi++) {
    for (Int_t cmj = 0; cmj < nbins; cmj++) {
      m_TUnfCovMatSys(cmi, cmj) = hsystcovmat->GetBinContent(cmi + 1, cmj + 1);
    }
  }

  /*
    TH1D* hist_acombfactor[nbinsrhoi/2];
    for (int i_a = 0; i_a < nbinsrhoi/2; ++i_a) {
      hist_acombfactor[i_a] = new TH1D(Form("hist_acombfactor%d",i_a), Form("hist_acombfactor%d",i_a), nVar,0,nVar);
    }
  */

  TH1D* hnom_AllVar = new TH1D("hnom_AllVar", "hnom_AllVar", fVariableNames.size()*nbinsrhoi,0,fVariableNames.size()*nbinsrhoi);
  int bin_counter = 1;

  for(int i=0;i<nVar;i++)
  {

    TString VariableName = fVariableNames.at(i);

    if(VariableName=="") continue;

    TString TUnfoldResultsFilename = Form(TUNFFILE + VariableName + ".root");
    TFile *unffile                 = new TFile(TUnfoldResultsFilename.Data(), "READ");
    
    if(unffile==NULL) exit(0);

    TH1D* tunfresult_temp = (TH1D*)((TH1D*)unffile->Get(VariableName+"TUnfResultCor"+NameTag))->Clone(VariableName+"TUnfResult_temp");
    if(inorm>0) tunfresult_temp->Scale(1./tunfresult_temp->Integral());

    /*
        TH1D* afbresult_temp = (TH1D*)((TH1D*)unffile->Get(VariableName+"AfbResult"+NameTag))->Clone(VariableName+"AfbResult_temp");
        
        // Fill factors used in combining binned asymmetries into coefficients
        TH1D *hGenMC_rebinnedS = (TH1D*)((TH1D*)unffile->Get(VariableName+"Gen_norm_rebinnedB"))->Clone(Form("hGenMC_rebinnedS_%s",VariableName.Data()));
        std::vector<double> factor {1.,1.,1.};
        MatrixUnf::AtoCfactor(VariableName, hGenMC_rebinnedS, factor);

        hist_acombfactor[0]->SetBinContent(i+1, (afbresult_temp->GetBinContent(0))/factor[0] );
        hist_acombfactor[1]->SetBinContent(i+1, (afbresult_temp->GetBinContent(nbinsrhoi+1))/factor[1] );
        hist_acombfactor[2]->SetBinContent(i+1, (1. - afbresult_temp->GetBinContent(0) - afbresult_temp->GetBinContent(nbinsrhoi+1))/factor[2] );
    */

    for(int ibins = 1; ibins<=tunfresult_temp->GetNbinsX(); ibins++){

      hnom_AllVar->SetBinContent(bin_counter,tunfresult_temp->GetBinContent(ibins));
      ++bin_counter;

    }

    delete tunfresult_temp;
    // delete afbresult_temp;
    // delete hGenMC_rebinnedS;

    unffile->Close();
    delete unffile;
  }


  // Calculate the sys on coef etc
  /*
    TMatrixD m_AFB(1,1);
    std::vector<double> afbvar;
    std::vector<double> afberrvar;
    MatrixUnf::GetAfbCorrMAllVars(tunfresult_temp, m_TUnfCovMatSys, afbvar, afberrvar, m_AFB, kTRUE);
    //m_AFB.Print("f=%1.5g ");

    TMatrixD m_BinSum(1,1);
    std::vector<double> binsum;
    std::vector<double> binsumerr;
    MatrixUnf::GetBinSumCorrMAllVars(tunfresult_temp, m_TUnfCovMatSys, binsum, binsumerr, m_BinSum, kTRUE);
    //m_BinSum.Print("f=%1.5g ");
  */


  TMatrixD m_BinSum(nVar,nVar);
  std::vector<double> binsum;
  std::vector<double> binsumerr;
  MatrixUnf::GetBinSumCorrMAllVars(hnom_AllVar, m_TUnfCovMatSys, binsum, binsumerr, m_BinSum, kTRUE);
  // m_BinSum.Print("f=%1.5g ");
  
  TVectorD diag_nozeros = TMatrixDDiag(m_BinSum);
  
  for(int xbin=1; xbin<=nVar; xbin++) {
    if( diag_nozeros(xbin-1) == 0 ) diag_nozeros(xbin-1) = 1e-20;
  }

  m_BinSum.NormByDiag(diag_nozeros);

  // std::cout<<"Coef sys corrm: "<<std::endl;
  // m_BinSum.Print("f=%1.5g ");

  TH2D* h_RhoCC = new TH2D(Form("RhoBSBS_SystTotal%s%s",(inorm?"Norm":""),NameTag.Data()), Form("RhoBSBS_SystTotal%s%s",(inorm?"Norm":""),NameTag.Data()), nVar, 0, nVar, nVar, 0, nVar);
  for(int i=0; i < nVar; i++){
    for(int j=0; j < nVar; j++){
      h_RhoCC->SetBinContent( i+1, j+1, m_BinSum(i,j) );
    }
  }

  // delete hist_acombfactor[0];
  // delete hist_acombfactor[1];
  // delete hist_acombfactor[2];

  delete hnom_AllVar;

  return h_RhoCC;
}


void MatrixUnfControl::CalculateHistRange(std::vector<TH1D*>PseudoUnfResults, int i_bin, double &mean, double &rms, bool iserror, bool wholerange) {

  int n = PseudoUnfResults.size();
  mean  = 0;
  rms   = 0;

  if(wholerange) {
    double min = 1e100;
    double max = -1e100;
    for(int i_pseudo=1; i_pseudo<n; i_pseudo++){
      if (!iserror) {
        min = PseudoUnfResults[i_pseudo]->GetBinContent(i_bin) < min ? PseudoUnfResults[i_pseudo]->GetBinContent(i_bin) : min;
        max = PseudoUnfResults[i_pseudo]->GetBinContent(i_bin) > max ? PseudoUnfResults[i_pseudo]->GetBinContent(i_bin) : max;
      }
      else {
        min = PseudoUnfResults[i_pseudo]->GetBinError(i_bin) < min ? PseudoUnfResults[i_pseudo]->GetBinError(i_bin) : min;
        max = PseudoUnfResults[i_pseudo]->GetBinError(i_bin) > max ? PseudoUnfResults[i_pseudo]->GetBinError(i_bin) : max;
      }
    }

    mean = (min+max)/2.;
    rms = (max-min)/8.;
    if (rms <= 0) rms = 1e-10;

  }
  else{ 
    for(int i_pseudo=1; i_pseudo<n; i_pseudo++){
      if (!iserror) {
        mean += PseudoUnfResults[i_pseudo]->GetBinContent(i_bin);
        rms += pow(PseudoUnfResults[i_pseudo]->GetBinContent(i_bin),2);
      }
      else {
        mean += PseudoUnfResults[i_pseudo]->GetBinError(i_bin);
        rms += pow(PseudoUnfResults[i_pseudo]->GetBinError(i_bin),2);
      }
    }

    mean /= (n-1);
    rms /= (n-1);
    rms -= pow(mean,2);
    if (rms <= 0) rms = 1e-20;
    rms = sqrt(rms);
    //std::cout<<i_bin<<" mean: "<<mean<<" rms: "<<rms<<std::endl;
  }

}


/**
 * @brief Creates all of the linearity summary histograms.
 * 
 * @param numCoefficients How many coefficients per distribution.
 */
void MatrixUnfControl::CreateLinearitySummaryHistos(int numCoefficients) {
  // TODO: This is assuming there is the same number of variables being extracted per distribution.
  gStyle->SetPadLeftMargin(0.08);
  gStyle->SetPadRightMargin(0.02);
  // gStyle->SetTitleSize(0.08, "XYZ");
  // gStyle->SetLabelSize(0.07, "XYZ");
  gStyle->SetTickLength(0.01, "XYZ");
  gStyle->SetOptTitle(0);
  gStyle->SetPadGridY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat("");

  int totalNumCoefficients = numCoefficients * fVariableNames.size();

  hLinearityMaxBiasFrac = new TH1D("hLinearityMaxBiasFrac", "hLinearityMaxBiasFrac", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityMaxBiasFracStat = new TH1D("hLinearityMaxBiasFracStat", "hLinearityMaxBiasFracStat", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityMaxBiasSlope = new TH1D("hLinearityMaxBiasSlope", "hLinearityMaxBiasSlope", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityMaxBiasFracMax = new TH1D("hLinearityMaxBiasFracMax", "hLinearityMaxBiasFracMax", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityMaxBiasFracStatMax = new TH1D("hLinearityMaxBiasFracStatMax", "hLinearityMaxBiasFracStatMax", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityMaxBiasSlopeMax = new TH1D("hLinearityMaxBiasSlopeMax", "hLinearityMaxBiasSlopeMax", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityAvgBiasFrac = new TH1D("hLinearityAvgBiasFrac", "hLinearityAvgBiasFrac", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityAvgBiasFracStat = new TH1D("hLinearityAvgBiasFracStat", "hLinearityAvgBiasFracStat", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityAvgBiasSlope = new TH1D("hLinearityAvgBiasSlope", "hLinearityAvgBiasSlope", totalNumCoefficients, 0, totalNumCoefficients);

  hLinearityMaxBiasFrac_rebinnedA = new TH1D("hLinearityMaxBiasFrac_rebinnedA", "hLinearityMaxBiasFrac_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityMaxBiasFracStat_rebinnedA = new TH1D("hLinearityMaxBiasFracStat_rebinnedA", "hLinearityMaxBiasFracStat_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityMaxBiasSlope_rebinnedA = new TH1D("hLinearityMaxBiasSlope_rebinnedA", "hLinearityMaxBiasSlope_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityMaxBiasFracMax_rebinnedA = new TH1D("hLinearityMaxBiasFracMax_rebinnedA", "hLinearityMaxBiasFracMax_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityMaxBiasFracStatMax_rebinnedA = new TH1D("hLinearityMaxBiasFracStatMax_rebinnedA", "hLinearityMaxBiasFracStatMax_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityMaxBiasSlopeMax_rebinnedA = new TH1D("hLinearityMaxBiasSlopeMax_rebinnedA", "hLinearityMaxBiasSlopeMax_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityAvgBiasFrac_rebinnedA = new TH1D("hLinearityAvgBiasFrac_rebinnedA", "hLinearityAvgBiasFrac_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityAvgBiasFracStat_rebinnedA = new TH1D("hLinearityAvgBiasFracStat_rebinnedA", "hLinearityAvgBiasFracStat_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityAvgBiasSlope_rebinnedA = new TH1D("hLinearityAvgBiasSlope_rebinnedA", "hLinearityAvgBiasSlope_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);

  hLinearityMaxBiasFrac_rebinnedB = new TH1D("hLinearityMaxBiasFrac_rebinnedB", "hLinearityMaxBiasFrac_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityMaxBiasFracStat_rebinnedB = new TH1D("hLinearityMaxBiasFracStat_rebinnedB", "hLinearityMaxBiasFracStat_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityMaxBiasSlope_rebinnedB = new TH1D("hLinearityMaxBiasSlope_rebinnedB", "hLinearityMaxBiasSlope_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityMaxBiasFracMax_rebinnedB = new TH1D("hLinearityMaxBiasFracMax_rebinnedB", "hLinearityMaxBiasFracMax_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityMaxBiasFracStatMax_rebinnedB = new TH1D("hLinearityMaxBiasFracStatMax_rebinnedB", "hLinearityMaxBiasFracStatMax_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityMaxBiasSlopeMax_rebinnedB = new TH1D("hLinearityMaxBiasSlopeMax_rebinnedB", "hLinearityMaxBiasSlopeMax_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityAvgBiasFrac_rebinnedB = new TH1D("hLinearityAvgBiasFrac_rebinnedB", "hLinearityAvgBiasFrac_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityAvgBiasFracStat_rebinnedB = new TH1D("hLinearityAvgBiasFracStat_rebinnedB", "hLinearityAvgBiasFracStat_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);
  hLinearityAvgBiasSlope_rebinnedB = new TH1D("hLinearityAvgBiasSlope_rebinnedB", "hLinearityAvgBiasSlope_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);

  hAfbLinearityMaxBiasFrac_rebinnedA = new TH1D("hAfbLinearityMaxBiasFrac_rebinnedA", "hAfbLinearityMaxBiasFrac_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);
  hAfbLinearityMaxBiasFracStat_rebinnedA = new TH1D("hAfbLinearityMaxBiasFracStat_rebinnedA", "hAfbLinearityMaxBiasFracStat_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);
  hAfbLinearityMaxBiasSlope_rebinnedA = new TH1D("hAfbLinearityMaxBiasSlope_rebinnedA", "hAfbLinearityMaxBiasSlope_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);
  hAfbLinearityMaxBiasFracMax_rebinnedA = new TH1D("hAfbLinearityMaxBiasFracMax_rebinnedA", "hAfbLinearityMaxBiasFracMax_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);
  hAfbLinearityMaxBiasFracStatMax_rebinnedA = new TH1D("hAfbLinearityMaxBiasFracStatMax_rebinnedA", "hAfbLinearityMaxBiasFracStatMax_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);
  hAfbLinearityMaxBiasSlopeMax_rebinnedA = new TH1D("hAfbLinearityMaxBiasSlopeMax_rebinnedA", "hAfbLinearityMaxBiasSlopeMax_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);
  hAfbLinearityAvgBiasFrac_rebinnedA = new TH1D("hAfbLinearityAvgBiasFrac_rebinnedA", "hAfbLinearityAvgBiasFrac_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);
  hAfbLinearityAvgBiasFracStat_rebinnedA = new TH1D("hAfbLinearityAvgBiasFracStat_rebinnedA", "hAfbLinearityAvgBiasFracStat_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);
  hAfbLinearityAvgBiasSlope_rebinnedA = new TH1D("hAfbLinearityAvgBiasSlope_rebinnedA", "hAfbLinearityAvgBiasSlope_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);

  hAfbLinearityMaxBiasFrac_rebinnedB = new TH1D("hAfbLinearityMaxBiasFrac_rebinnedB", "hAfbLinearityMaxBiasFrac_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);
  hAfbLinearityMaxBiasFracStat_rebinnedB = new TH1D("hAfbLinearityMaxBiasFracStat_rebinnedB", "hAfbLinearityMaxBiasFracStat_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);
  hAfbLinearityMaxBiasSlope_rebinnedB = new TH1D("hAfbLinearityMaxBiasSlope_rebinnedB", "hAfbLinearityMaxBiasSlope_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);
  hAfbLinearityMaxBiasFracMax_rebinnedB = new TH1D("hAfbLinearityMaxBiasFracMax_rebinnedB", "hAfbLinearityMaxBiasFracMax_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);
  hAfbLinearityMaxBiasFracStatMax_rebinnedB = new TH1D("hAfbLinearityMaxBiasFracStatMax_rebinnedB", "hAfbLinearityMaxBiasFracStatMax_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);
  hAfbLinearityMaxBiasSlopeMax_rebinnedB = new TH1D("hAfbLinearityMaxBiasSlopeMax_rebinnedB", "hAfbLinearityMaxBiasSlopeMax_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);
  hAfbLinearityAvgBiasFrac_rebinnedB = new TH1D("hAfbLinearityAvgBiasFrac_rebinnedB", "hAfbLinearityAvgBiasFrac_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);
  hAfbLinearityAvgBiasFracStat_rebinnedB = new TH1D("hAfbLinearityAvgBiasFracStat_rebinnedB", "hAfbLinearityAvgBiasFracStat_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);
  hAfbLinearityAvgBiasSlope_rebinnedB = new TH1D("hAfbLinearityAvgBiasSlope_rebinnedB", "hAfbLinearityAvgBiasSlope_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);

  // unused hAfbConstantMaxBiasFrac_rebinnedA = new TH1D("hAfbConstantMaxBiasFrac_rebinnedA", "hAfbConstantMaxBiasFrac_rebinnedA", totalNumCoefficients,0,totalNumCoefficients);
  // unused hAfbConstantMaxBiasFracStat_rebinnedA = new TH1D("hAfbConstantMaxBiasFracStat_rebinnedA", "hAfbConstantMaxBiasFracStat_rebinnedA", totalNumCoefficients,0,totalNumCoefficients);
  hAfbConstantMaxBiasSlope_rebinnedA = new TH1D("hAfbConstantMaxBiasSlope_rebinnedA", "hAfbConstantMaxBiasSlope_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);
  // unused hAfbConstantMaxBiasFracMax_rebinnedA = new TH1D("hAfbConstantMaxBiasFracMax_rebinnedA", "hAfbConstantMaxBiasFracMax_rebinnedA", totalNumCoefficients,0,totalNumCoefficients);
  // unused hAfbConstantMaxBiasFracStatMax_rebinnedA = new TH1D("hAfbConstantMaxBiasFracStatMax_rebinnedA", "hAfbConstantMaxBiasFracStatMax_rebinnedA", totalNumCoefficients,0,totalNumCoefficients);
  // unused hAfbConstantMaxBiasSlopeMax_rebinnedA = new TH1D("hAfbConstantMaxBiasSlopeMax_rebinnedA", "hAfbConstantMaxBiasSlopeMax_rebinnedA", totalNumCoefficients,0,totalNumCoefficients);
  // unused hAfbConstantAvgBiasFrac_rebinnedA = new TH1D("hAfbConstantAvgBiasFrac_rebinnedA", "hAfbConstantAvgBiasFrac_rebinnedA", totalNumCoefficients,0,totalNumCoefficients);
  // unused hAfbConstantAvgBiasFracStat_rebinnedA = new TH1D("hAfbConstantAvgBiasFracStat_rebinnedA", "hAfbConstantAvgBiasFracStat_rebinnedA", totalNumCoefficients,0,totalNumCoefficients);
  hAfbConstantAvgBiasSlope_rebinnedA = new TH1D("hAfbConstantAvgBiasSlope_rebinnedA", "hAfbConstantAvgBiasSlope_rebinnedA", totalNumCoefficients, 0, totalNumCoefficients);

  // unused hAfbConstantMaxBiasFrac_rebinnedB = new TH1D("hAfbConstantMaxBiasFrac_rebinnedB", "hAfbConstantMaxBiasFrac_rebinnedB", totalNumCoefficients,0,totalNumCoefficients);
  // unused hAfbConstantMaxBiasFracStat_rebinnedB = new TH1D("hAfbConstantMaxBiasFracStat_rebinnedB", "hAfbConstantMaxBiasFracStat_rebinnedB", totalNumCoefficients,0,totalNumCoefficients);
  hAfbConstantMaxBiasSlope_rebinnedB = new TH1D("hAfbConstantMaxBiasSlope_rebinnedB", "hAfbConstantMaxBiasSlope_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);
  // unused hAfbConstantMaxBiasFracMax_rebinnedB = new TH1D("hAfbConstantMaxBiasFracMax_rebinnedB", "hAfbConstantMaxBiasFracMax_rebinnedB", totalNumCoefficients,0,totalNumCoefficients);
  // unused hAfbConstantMaxBiasFracStatMax_rebinnedB = new TH1D("hAfbConstantMaxBiasFracStatMax_rebinnedB", "hAfbConstantMaxBiasFracStatMax_rebinnedB", totalNumCoefficients,0,totalNumCoefficients);
  // unused hAfbConstantMaxBiasSlopeMax_rebinnedB = new TH1D("hAfbConstantMaxBiasSlopeMax_rebinnedB", "hAfbConstantMaxBiasSlopeMax_rebinnedB", totalNumCoefficients,0,totalNumCoefficients);
  // unused hAfbConstantAvgBiasFrac_rebinnedB = new TH1D("hAfbConstantAvgBiasFrac_rebinnedB", "hAfbConstantAvgBiasFrac_rebinnedB", totalNumCoefficients,0,totalNumCoefficients);
  // unused hAfbConstantAvgBiasFracStat_rebinnedB = new TH1D("hAfbConstantAvgBiasFracStat_rebinnedB", "hAfbConstantAvgBiasFracStat_rebinnedB", totalNumCoefficients,0,totalNumCoefficients);
  hAfbConstantAvgBiasSlope_rebinnedB = new TH1D("hAfbConstantAvgBiasSlope_rebinnedB", "hAfbConstantAvgBiasSlope_rebinnedB", totalNumCoefficients, 0, totalNumCoefficients);
}


void MatrixUnfControl::simpleFit(TGraph* graph, TH1D* hist, TH1D* histMax, int i_bin, double DeltaCtotal, double &param, double &parammax, bool isslope) {
    
  double mean_temp,slope_temp;
  int ifail;
  int npoints = graph->GetN();

  graph->LeastSquareLinearFit(1,mean_temp,slope_temp,ifail);

  if(!isslope) {
    param = slope_temp;

    double maxg = TMath::MaxElement(npoints,graph->GetY())-graph->GetY()[(npoints-1)/2];
    double ming = TMath::MinElement(npoints,graph->GetY())-graph->GetY()[(npoints-1)/2];
    parammax = std::max(fabs(maxg),fabs(ming)) * fabs(slope_temp)/(fabs(slope_temp)>0?slope_temp:1e-25) * 2./(DeltaCtotal);
  }
  else {
    param = mean_temp;

    double maxg = TMath::MaxElement(npoints,graph->GetY());
    double ming = TMath::MinElement(npoints,graph->GetY());
    parammax = mean_temp>1.?maxg:ming;
  }

  hist->SetBinContent(i_bin, param);
  hist->SetBinError(i_bin, 0.);

  histMax->SetBinContent(i_bin, parammax);
  histMax->SetBinError(i_bin, 0.);

}

/**
 * @brief 
 * 
 * @param rebin 
 * @param isasym 
 * @param nametag 
 * @param nLinearity 
 * @param i 
 * @param unf_nbins 
 * @param hLinearityGen 
 * @param delta_coefficient 
 * @param hLinearityDBiasFractional 
 * @param hLinearityDBiasFractionalMax 
 * @param hLinearityDBias 
 * @param hLinearityDBiasMax 
 * @param hLinearityDInjectedShape 
 * @param hLinearityDInjectedShapeMax 
 * @param hLinearityDMeasuredShape 
 * @param hLinearityDMeasuredShapeMax 
 * @param hLinearitySlope 
 * @param hLinearitySlopeMax 
 * @param hLinearityDBiasOverStatErr 
 * @param hLinearityDBiasOverStatErrMax 
 * @param hLinearityCBias 
 * @param hLinearityCBiasFractional 
 * @param hLinearityCBiasOverStatErr 
 * @param hLinearityMaxBiasFrac 
 * @param hLinearityMaxBiasFracStat 
 * @param hLinearityMaxBiasSlope 
 * @param hLinearityMaxBiasFracMax 
 * @param hLinearityMaxBiasFracStatMax 
 * @param hLinearityMaxBiasSlopeMax 
 * @param hLinearityAvgBiasFrac 
 * @param hLinearityAvgBiasFracStat 
 * @param hLinearityAvgBiasSlope 
 * @param hConstantMaxBiasSlope 
 * @param hConstantAvgBiasSlope 
 * @param numCoefficients 
 */
void MatrixUnfControl::fillLinearityGraphs(std::vector<TH1D*>PseudoUnfResults, int rebin, bool isasym, TString nametag, int nLinearity, int i, int unf_nbins, TH1D* hLinearityGen[], double delta_coefficient[], TH1D* hLinearityDBiasFractional, TH1D* hLinearityDBiasFractionalMax, TH1D* hLinearityDBias, TH1D* hLinearityDBiasMax, TH1D* hLinearityDInjectedShape, TH1D* hLinearityDInjectedShapeMax, TH1D* hLinearityDMeasuredShape, TH1D* hLinearityDMeasuredShapeMax, TH1D* hLinearitySlope, TH1D* hLinearitySlopeMax, TH1D* hLinearityDBiasOverStatErr, TH1D* hLinearityDBiasOverStatErrMax, TH1D* hLinearityCBias, TH1D* hLinearityCBiasFractional, TH1D* hLinearityCBiasOverStatErr,    TH1D* hLinearityMaxBiasFrac, TH1D* hLinearityMaxBiasFracStat, TH1D* hLinearityMaxBiasSlope, TH1D* hLinearityMaxBiasFracMax, TH1D* hLinearityMaxBiasFracStatMax, TH1D* hLinearityMaxBiasSlopeMax, TH1D* hLinearityAvgBiasFrac, TH1D* hLinearityAvgBiasFracStat, TH1D* hLinearityAvgBiasSlope, TH1D* hConstantMaxBiasSlope, TH1D* hConstantAvgBiasSlope, int numCoefficients) {
  hLinearityDBias->Reset();
  hLinearityCBias->Reset();
  hLinearityDInjectedShape->Reset();
  hLinearityDMeasuredShape->Reset();
  hLinearitySlope->Reset();
  hLinearityDBiasOverStatErr->Reset();
  hLinearityCBiasOverStatErr->Reset();
  hLinearityDBiasFractional->Reset();
  hLinearityCBiasFractional->Reset();
  hLinearityDBiasMax->Reset();
  //unused hLinearityCBiasMax->Reset();
  hLinearityDInjectedShapeMax->Reset();
  hLinearityDMeasuredShapeMax->Reset();
  hLinearitySlopeMax->Reset();
  hLinearityDBiasOverStatErrMax->Reset();
  //unused hLinearityCBiasOverStatErrMax->Reset();
  hLinearityDBiasFractionalMax->Reset();
  //unused hLinearityCBiasFractionalMax->Reset();

  int numBinsPerCoefficient = unf_nbins / numCoefficients;
	int variableOffset = i * numCoefficients;
	const char *varName = fVariableNames[i].Data();
	const char *name = nametag.Data();
	std::string asym = "";
	if (isasym) {
		asym = "Afb";
	}

  TGraph *gLinearityDBias[(numCoefficients * (numBinsPerCoefficient + 2)) - 2];
  TGraph *gLinearityDInjectedShape[(numCoefficients * (numBinsPerCoefficient + 2)) - 2];
  TGraph *gLinearityDMeasuredShape[(numCoefficients * (numBinsPerCoefficient + 2)) - 2];
  TGraph *gLinearitySlope[(numCoefficients * (numBinsPerCoefficient + 2)) - 2];
  TGraph *gLinearityForAfiq[(numCoefficients * (numBinsPerCoefficient + 2)) - 2];
  TGraph *gLinearityDBiasOverStatErr[(numCoefficients * (numBinsPerCoefficient + 2)) - 2];
  TGraph *gLinearityDBiasFractional[(numCoefficients * (numBinsPerCoefficient + 2)) - 2];

  double MaxBiasFrac[numCoefficients], MaxBiasFracStat[numCoefficients], MaxBiasSlope[numCoefficients];
  double MaxBiasFracMax[numCoefficients], MaxBiasFracStatMax[numCoefficients], MaxBiasSlopeMax[numCoefficients];
  double AvgBiasFrac[numCoefficients], AvgBiasFracStat[numCoefficients], AvgBiasSlope[numCoefficients];

  for (int i_coefficient = 0; i_coefficient < numCoefficients; i_coefficient++) {
    int binOffset = 0;
    int arrayOffset = 0;
    if (isasym) {
      binOffset = (numBinsPerCoefficient + 2) * i_coefficient;
      arrayOffset = (numBinsPerCoefficient + 2) * i_coefficient;
    } else {
      binOffset = numBinsPerCoefficient * i_coefficient;
      arrayOffset = numBinsPerCoefficient * i_coefficient;
    }

		MaxBiasFrac[i_coefficient] = -1.;
		MaxBiasFracStat[i_coefficient] = -1.;
		MaxBiasSlope[i_coefficient] = -1.;
		MaxBiasFracMax[i_coefficient] = -1.;
		MaxBiasFracStatMax[i_coefficient] = -1.;
		MaxBiasSlopeMax[i_coefficient] = -1.;
		AvgBiasFrac[i_coefficient] = 0.;
		AvgBiasFracStat[i_coefficient] = 0.;
		AvgBiasSlope[i_coefficient] = 0.;

    std::vector<int> SlopePointsToRemove;
    for(int i_bin = 1; i_bin <= numBinsPerCoefficient; i_bin++) {
      gLinearityDBias[arrayOffset + i_bin - 1] = new TGraph(nLinearity);
			gLinearityDBias[arrayOffset + i_bin - 1]->SetName(Form("%s_%sLinearityDBias%s_bin%d", varName, asym.c_str(), name, arrayOffset + i_bin));
			gLinearityDInjectedShape[arrayOffset + i_bin - 1] = new TGraph(nLinearity);
			gLinearityDInjectedShape[arrayOffset + i_bin - 1]->SetName(Form("%s_%sLinearityDInjectedShape%s_bin%d", varName, asym.c_str(), name, arrayOffset + i_bin));
			gLinearityDMeasuredShape[arrayOffset + i_bin - 1] = new TGraph(nLinearity);
			gLinearityDMeasuredShape[arrayOffset + i_bin - 1]->SetName(Form("%s_%sLinearityDMeasuredShape%s_bin%d", varName, asym.c_str(), name, arrayOffset + i_bin));
			gLinearitySlope[arrayOffset + i_bin - 1] = new TGraph(nLinearity);
			gLinearitySlope[arrayOffset + i_bin - 1]->SetName(Form("%s_%sLinearitySlope%s_bin%d", varName, asym.c_str(), name, arrayOffset + i_bin));
			gLinearityForAfiq[arrayOffset + i_bin - 1] = new TGraph(nLinearity);
			gLinearityForAfiq[arrayOffset + i_bin - 1]->SetName(Form("%s_%sLinearityForAfiq%s_bin%d", varName, asym.c_str(), name, arrayOffset + i_bin));
			gLinearityDBiasOverStatErr[arrayOffset + i_bin - 1] = new TGraph(nLinearity);
			gLinearityDBiasOverStatErr[arrayOffset + i_bin - 1]->SetName(Form("%s_%sLinearityDBiasOverStatErr%s_bin%d", varName, asym.c_str(), name, arrayOffset + i_bin));
			gLinearityDBiasFractional[arrayOffset + i_bin - 1] = new TGraph(nLinearity);
			gLinearityDBiasFractional[arrayOffset + i_bin - 1]->SetName(Form("%s_%sLinearityDBiasFractional%s_bin%d", varName, asym.c_str(), name, arrayOffset + i_bin));

			SlopePointsToRemove.clear();
			for (int i_linearity = 0; i_linearity < nLinearity; i_linearity++) {
				double deltaBias = PseudoUnfResults[i_linearity]->GetBinContent(binOffset + i_bin) - hLinearityGen[i_linearity]->GetBinContent(arrayOffset + i_bin);
				double deltaAfb  = PseudoUnfResults[i_linearity]->GetBinContent(binOffset + i_bin) - PseudoUnfResults[(nLinearity - 1) / 2]->GetBinContent(binOffset + i_bin);
				double afbError  = PseudoUnfResults[i_linearity]->GetBinError(binOffset + i_bin);
				double injectedShape = hLinearityGen[i_linearity]->GetBinContent(arrayOffset + i_bin) - hLinearityGen[(nLinearity - 1) / 2]->GetBinContent(arrayOffset + i_bin);

				gLinearityDBias[arrayOffset + i_bin - 1]->SetPoint(i_linearity, delta_coefficient[i_linearity], deltaBias);
				gLinearityDMeasuredShape[arrayOffset + i_bin - 1]->SetPoint(i_linearity, delta_coefficient[i_linearity], deltaAfb);
				gLinearityDBiasOverStatErr[arrayOffset + i_bin - 1]->SetPoint(i_linearity, delta_coefficient[i_linearity], deltaBias / afbError);
				gLinearityDInjectedShape[arrayOffset + i_bin - 1]->SetPoint(i_linearity, delta_coefficient[i_linearity], injectedShape);

				// as long as contents in the bin i_bin in the i_linearity variation changed by 0.001% fill plots with results
				if (fabs(1. - hLinearityGen[i_linearity]->GetBinContent(arrayOffset + i_bin) / hLinearityGen[(nLinearity - 1) / 2]->GetBinContent(arrayOffset + i_bin)) > 1e-5) {
					gLinearitySlope[arrayOffset + i_bin - 1]->SetPoint(i_linearity, delta_coefficient[i_linearity], deltaAfb / injectedShape);
					gLinearityForAfiq[arrayOffset + i_bin - 1]->SetPoint(i_linearity, hLinearityGen[i_linearity]->GetBinContent(arrayOffset + i_bin), PseudoUnfResults[i_linearity]->GetBinContent(binOffset + i_bin));
				} else
					SlopePointsToRemove.push_back(i_linearity);

				if (isasym)
					gLinearityDBiasFractional[arrayOffset + i_bin - 1]->SetPoint(i_linearity, delta_coefficient[i_linearity], deltaBias);
				else
					gLinearityDBiasFractional[arrayOffset + i_bin - 1]->SetPoint(i_linearity, delta_coefficient[i_linearity], hLinearityGen[i_linearity]->GetBinContent(arrayOffset + i_bin) > 0 ? deltaBias / hLinearityGen[i_linearity]->GetBinContent(arrayOffset + i_bin) : 1.);
			}

			for (int iremove = SlopePointsToRemove.size() - 1; iremove >= 0; --iremove) {
				gLinearitySlope[arrayOffset + i_bin - 1]->RemovePoint(SlopePointsToRemove[iremove]);
				gLinearityForAfiq[arrayOffset + i_bin - 1]->RemovePoint(SlopePointsToRemove[iremove]);
			}

			double param, parammax;

			MatrixUnfControl::simpleFit(gLinearityDBiasFractional[arrayOffset + i_bin - 1], hLinearityDBiasFractional, hLinearityDBiasFractionalMax, arrayOffset + i_bin, delta_coefficient[nLinearity - 1] - delta_coefficient[0], param, parammax, false);
			if (fabs(param) > MaxBiasFrac[i_coefficient]) MaxBiasFrac[i_coefficient] = fabs(param);
			if (fabs(parammax) > MaxBiasFracMax[i_coefficient]) MaxBiasFracMax[i_coefficient] = fabs(parammax);
			AvgBiasFrac[i_coefficient] += fabs(param);

			MatrixUnfControl::simpleFit(gLinearityDBiasOverStatErr[arrayOffset + i_bin - 1], hLinearityDBiasOverStatErr, hLinearityDBiasOverStatErrMax, arrayOffset + i_bin, delta_coefficient[nLinearity - 1] - delta_coefficient[0], param, parammax, false);
			if (fabs(param) > MaxBiasFracStat[i_coefficient]) MaxBiasFracStat[i_coefficient] = fabs(param);
			if (fabs(parammax) > MaxBiasFracStatMax[i_coefficient]) MaxBiasFracStatMax[i_coefficient] = fabs(parammax);
			AvgBiasFracStat[i_coefficient] += fabs(param);

			MatrixUnfControl::simpleFit(gLinearityDBias[arrayOffset + i_bin - 1], hLinearityDBias, hLinearityDBiasMax, arrayOffset + i_bin, delta_coefficient[nLinearity - 1] - delta_coefficient[0], param, parammax, false);

			MatrixUnfControl::simpleFit(gLinearityDInjectedShape[arrayOffset + i_bin - 1], hLinearityDInjectedShape, hLinearityDInjectedShapeMax, arrayOffset + i_bin, delta_coefficient[nLinearity - 1] - delta_coefficient[0], param, parammax, false);

			MatrixUnfControl::simpleFit(gLinearityDMeasuredShape[arrayOffset + i_bin - 1], hLinearityDMeasuredShape, hLinearityDMeasuredShapeMax, arrayOffset + i_bin, delta_coefficient[nLinearity - 1] - delta_coefficient[0], param, parammax, false);

			if (SlopePointsToRemove.size() < nLinearity) {
				MatrixUnfControl::simpleFit(gLinearitySlope[arrayOffset + i_bin - 1], hLinearitySlope, hLinearitySlopeMax, arrayOffset + i_bin, delta_coefficient[nLinearity - 1] - delta_coefficient[0], param, parammax, true);
				if (fabs(param - 1.) > MaxBiasSlope[i_coefficient]) MaxBiasSlope[i_coefficient] = fabs(param - 1.);
				if (fabs(parammax - 1.) > MaxBiasSlopeMax[i_coefficient]) MaxBiasSlopeMax[i_coefficient] = fabs(parammax - 1.);
				AvgBiasSlope[i_coefficient] += fabs(param - 1.);
			} else {
				hLinearitySlope->SetBinContent(arrayOffset + i_bin, 1.);
				hLinearitySlopeMax->SetBinContent(arrayOffset + i_bin, 1.);
			}

			hLinearityCBias->SetBinContent(arrayOffset + i_bin, gLinearityDBias[arrayOffset + i_bin - 1]->GetY()[(nLinearity - 1) / 2]);
			hLinearityCBias->SetBinError(arrayOffset + i_bin, 0.);

			hLinearityCBiasFractional->SetBinContent(arrayOffset + i_bin, gLinearityDBiasFractional[arrayOffset + i_bin - 1]->GetY()[(nLinearity - 1) / 2]);
			hLinearityCBiasFractional->SetBinError(arrayOffset + i_bin, 0.);

			hLinearityCBiasOverStatErr->SetBinContent(arrayOffset + i_bin, gLinearityDBiasOverStatErr[arrayOffset + i_bin - 1]->GetY()[(nLinearity - 1) / 2]);
			hLinearityCBiasOverStatErr->SetBinError(arrayOffset + i_bin, 0.);

			gLinearityDBias[arrayOffset + i_bin - 1]->Write();
			gLinearityDInjectedShape[arrayOffset + i_bin - 1]->Write();
			gLinearityDMeasuredShape[arrayOffset + i_bin - 1]->Write();
			gLinearitySlope[arrayOffset + i_bin - 1]->Write();
			gLinearityForAfiq[arrayOffset + i_bin - 1]->Write();
			gLinearityDBiasOverStatErr[arrayOffset + i_bin - 1]->Write();
			gLinearityDBiasFractional[arrayOffset + i_bin - 1]->Write(); 
    }

    AvgBiasFrac[i_coefficient] /= numBinsPerCoefficient;
		AvgBiasFracStat[i_coefficient] /= numBinsPerCoefficient;
		AvgBiasSlope[i_coefficient] /= numBinsPerCoefficient;


    if (isasym) {
      // hLinearityMaxBiasFrac->SetBinContent(i+1,MaxBiasFrac);
			// hLinearityMaxBiasFracMax->SetBinContent(i+1,MaxBiasFracMax);
			// hLinearityAvgBiasFrac->SetBinContent(i+1,AvgBiasFrac);
			// hLinearityMaxBiasFracStat->SetBinContent(i+1,MaxBiasFracStat);
			// hLinearityMaxBiasFracStatMax->SetBinContent(i+1,MaxBiasFracStatMax);
			// hLinearityAvgBiasFracStat->SetBinContent(i+1,AvgBiasFracStat);
			// hLinearityMaxBiasSlope->SetBinContent(i+1,MaxBiasSlope);
			// hLinearityMaxBiasSlopeMax->SetBinContent(i+1,MaxBiasSlopeMax);
			// hLinearityAvgBiasSlope->SetBinContent(i+1,AvgBiasSlope);

			hLinearityMaxBiasFrac->SetBinContent(variableOffset + i_coefficient + 1, hLinearityDBiasFractional->GetBinContent(binOffset + numBinsPerCoefficient - 1));
			hLinearityMaxBiasFracMax->SetBinContent(variableOffset + i_coefficient + 1, hLinearityDBiasFractionalMax->GetBinContent(binOffset + numBinsPerCoefficient - 1));
			hLinearityAvgBiasFrac->SetBinContent(variableOffset + i_coefficient + 1, hLinearityDBiasFractional->GetBinContent(binOffset + numBinsPerCoefficient));
			hLinearityMaxBiasFracStat->SetBinContent(variableOffset + i_coefficient + 1, hLinearityDBiasOverStatErr->GetBinContent(binOffset + numBinsPerCoefficient - 1));
			hLinearityMaxBiasFracStatMax->SetBinContent(variableOffset + i_coefficient + 1, hLinearityDBiasOverStatErrMax->GetBinContent(binOffset + numBinsPerCoefficient - 1));
			hLinearityAvgBiasFracStat->SetBinContent(variableOffset + i_coefficient + 1, hLinearityDBiasOverStatErr->GetBinContent(binOffset + numBinsPerCoefficient));
			hLinearityMaxBiasSlope->SetBinContent(variableOffset + i_coefficient + 1, hLinearitySlope->GetBinContent(binOffset + numBinsPerCoefficient - 1));
			hLinearityMaxBiasSlopeMax->SetBinContent(variableOffset + i_coefficient + 1, hLinearitySlopeMax->GetBinContent(binOffset + numBinsPerCoefficient - 1));
			hLinearityAvgBiasSlope->SetBinContent(variableOffset + i_coefficient + 1, hLinearitySlope->GetBinContent(binOffset + numBinsPerCoefficient));

			hConstantMaxBiasSlope->SetBinContent(variableOffset + i_coefficient + 1, hLinearityCBias->GetBinContent(binOffset + numBinsPerCoefficient - 1));
			hConstantAvgBiasSlope->SetBinContent(variableOffset + i_coefficient + 1, hLinearityCBias->GetBinContent(binOffset + numBinsPerCoefficient));
    } else {
      hLinearityMaxBiasFrac->SetBinContent(variableOffset + i_coefficient + 1, MaxBiasFrac[i_coefficient]);
			hLinearityMaxBiasFracMax->SetBinContent(variableOffset + i_coefficient + 1, MaxBiasFracMax[i_coefficient]);
			hLinearityAvgBiasFrac->SetBinContent(variableOffset + i_coefficient + 1, AvgBiasFrac[i_coefficient]);
			hLinearityMaxBiasFracStat->SetBinContent(variableOffset + i_coefficient + 1, MaxBiasFracStat[i_coefficient]);
			hLinearityMaxBiasFracStatMax->SetBinContent(variableOffset + i_coefficient + 1, MaxBiasFracStatMax[i_coefficient]);
			hLinearityAvgBiasFracStat->SetBinContent(variableOffset + i_coefficient + 1, AvgBiasFracStat[i_coefficient]);
			hLinearityMaxBiasSlope->SetBinContent(variableOffset + i_coefficient + 1, MaxBiasSlope[i_coefficient]);
			hLinearityMaxBiasSlopeMax->SetBinContent(variableOffset + i_coefficient + 1, MaxBiasSlopeMax[i_coefficient]);
			hLinearityAvgBiasSlope->SetBinContent(variableOffset + i_coefficient + 1, AvgBiasSlope[i_coefficient]);

			if (rebin > 1) {
				// scale the rebinned absolute histos so they can be drawn on the same axis
				hLinearityDBias->Scale(1. / rebin);
				hLinearityDInjectedShape->Scale(1. / rebin);
				hLinearityDMeasuredShape->Scale(1. / rebin);
				hLinearityDBiasMax->Scale(1. / rebin);
				hLinearityDInjectedShapeMax->Scale(1. / rebin);
				hLinearityDMeasuredShapeMax->Scale(1. / rebin);

				hLinearityCBias->Scale(1. / rebin);
				// unused hLinearityCBiasMax->Scale(1./rebin);
			}
    }
  }
  TString quantity = MatrixUnf::ObsName(fVariableNames[i], kFALSE);
  TString quantitycoef = MatrixUnf::CoefName(fVariableNames[i], kFALSE);

  hLinearityDBias->GetYaxis()->SetTitle("#DeltaBias / #Delta coef. (avg)");
  hLinearityDBias->GetXaxis()->SetTitle(Form("%s", quantity.Data()));
  hLinearityDInjectedShape->GetYaxis()->SetTitle("Injected shape, #Delta(bin value) / #Delta coef. (avg)");
  hLinearityDInjectedShape->GetXaxis()->SetTitle(Form("%s", quantity.Data()));
  hLinearityDMeasuredShape->GetYaxis()->SetTitle("Measured shape, #Delta(bin value) / #Delta coef. (avg)");
  hLinearityDMeasuredShape->GetXaxis()->SetTitle(Form("%s", quantity.Data()));
  hLinearitySlope->GetYaxis()->SetTitle("Linearity slope, #Delta(measured shape) / #Delta(injected shape) (avg)");
  hLinearitySlope->GetXaxis()->SetTitle(Form("%s", quantity.Data()));
  hLinearityDBiasOverStatErr->GetYaxis()->SetTitle("#Delta(bias/#sigma_{stat}) / #Delta coef. (avg)");
  hLinearityDBiasOverStatErr->GetXaxis()->SetTitle(Form("%s", quantity.Data()));

  hLinearityDBiasMax->GetYaxis()->SetTitle("#DeltaBias / #Delta coef.");
  hLinearityDBiasMax->GetXaxis()->SetTitle(Form("%s", quantity.Data()));
  hLinearityDInjectedShapeMax->GetYaxis()->SetTitleOffset(1.8);
  hLinearityDInjectedShapeMax->GetYaxis()->SetTitle("Injected shape, #Delta(bin value) / #Delta coef.");
  hLinearityDInjectedShapeMax->GetXaxis()->SetTitle(Form("%s", quantity.Data()));
  hLinearityDMeasuredShapeMax->GetYaxis()->SetTitle("Measured shape, #Delta(bin value) / #Delta coef.");
  hLinearityDMeasuredShapeMax->GetXaxis()->SetTitle(Form("%s", quantity.Data()));
  hLinearitySlopeMax->GetYaxis()->SetTitleOffset(1.5);
  hLinearitySlopeMax->GetYaxis()->SetTitle("Linearity slope, #Delta(measured shape) / #Delta(injected shape)");
  hLinearitySlopeMax->GetXaxis()->SetTitle(Form("%s", quantity.Data()));
  hLinearityDBiasOverStatErrMax->GetYaxis()->SetTitle("#Delta(bias/#sigma_{stat}) / #Delta coef.");
  hLinearityDBiasOverStatErrMax->GetXaxis()->SetTitle(Form("%s", quantity.Data()));

  hLinearityCBias->GetYaxis()->SetTitle("Bias (constant)");
  hLinearityCBias->GetXaxis()->SetTitle(Form("%s", quantity.Data()));
  hLinearityCBiasOverStatErr->GetYaxis()->SetTitle("Bias/StatUnc. (constant)");
  hLinearityCBiasOverStatErr->GetXaxis()->SetTitle(Form("%s", quantity.Data()));
  // unused hLinearityCBiasMax->GetYaxis()->SetTitle("Bias");
  // unused hLinearityCBiasMax->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
  // unused hLinearityCBiasOverStatErrMax->GetYaxis()->SetTitle("Bias/StatUnc.");
  // unused hLinearityCBiasOverStatErrMax->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
  // unused hLinearityCBiasFractionalMax->GetYaxis()->SetTitle("FractionalBias");
  // unused hLinearityCBiasFractionalMax->GetXaxis()->SetTitle(Form("%s",quantity.Data()));

  if (isasym) {
    hLinearityDBiasFractional->GetYaxis()->SetTitle("#DeltaBias / #Delta coef. (avg)");
    hLinearityDBiasFractional->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, C(A_{total}) (right)", quantity.Data(), quantitycoef.Data()));

    hLinearityDBiasFractionalMax->GetYaxis()->SetTitle("#DeltaBias / #Delta coef.");
    hLinearityDBiasFractionalMax->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, C(A_{total}) (right)", quantity.Data(), quantitycoef.Data()));

    hLinearityCBiasFractional->GetYaxis()->SetTitle("Bias (constant)");
    hLinearityCBiasFractional->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, C(A_{total}) (right)", quantity.Data(), quantitycoef.Data()));

    hLinearityDBias->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, C(A_{total}) (right)", quantity.Data(), quantitycoef.Data()));
    hLinearityDInjectedShape->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, C(A_{total}) (right)", quantity.Data(), quantitycoef.Data()));
    hLinearityDMeasuredShape->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, C(A_{total}) (right)", quantity.Data(), quantitycoef.Data()));
    hLinearitySlope->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, C(A_{total}) (right)", quantity.Data(), quantitycoef.Data()));
    hLinearityDBiasOverStatErr->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, C(A_{total}) (right)", quantity.Data(), quantitycoef.Data()));
    hLinearityDBiasMax->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, C(A_{total}) (right)", quantity.Data(), quantitycoef.Data()));
    hLinearityDInjectedShapeMax->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, C(A_{total}) (right)", quantity.Data(), quantitycoef.Data()));
    hLinearityDMeasuredShapeMax->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, C(A_{total}) (right)", quantity.Data(), quantitycoef.Data()));
    hLinearitySlopeMax->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, C(A_{total}) (right)", quantity.Data(), quantitycoef.Data()));
    hLinearityDBiasOverStatErrMax->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, C(A_{total}) (right)", quantity.Data(), quantitycoef.Data()));
    hLinearityCBias->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, C(A_{total}) (right)", quantity.Data(), quantitycoef.Data()));
    hLinearityCBiasOverStatErr->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, C(A_{total}) (right)", quantity.Data(), quantitycoef.Data()));
  } else {
    hLinearityDBiasFractional->GetYaxis()->SetTitle("#Delta(fractional bias) / #Delta coef. (avg)");
    hLinearityDBiasFractional->GetXaxis()->SetTitle(Form("%s", quantity.Data()));

    hLinearityDBiasFractionalMax->GetYaxis()->SetTitle("#Delta(fractional bias) / #Delta coef.");
    hLinearityDBiasFractionalMax->GetXaxis()->SetTitle(Form("%s", quantity.Data()));

    hLinearityCBiasFractional->GetYaxis()->SetTitle("Fractional bias (constant)");
    hLinearityCBiasFractional->GetXaxis()->SetTitle(Form("%s", quantity.Data()));
  }

  hLinearityDBiasFractionalMax->GetYaxis()->SetTitleOffset(1.5);
}


void MatrixUnfControl::doPlotPulls(TString TUNFFILE, TString  PLOTFILE){


  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.02);
  //gStyle->SetTitleSize(0.08, "XYZ");
  //gStyle->SetLabelSize(0.07, "XYZ");
  gStyle->SetTickLength(0.01, "XYZ");
  gStyle->SetOptTitle(0);
  gStyle->SetPadGridY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat("eMR");
  gStyle->SetOptFit(1111);
  gStyle->SetStatW(0.35);
  gStyle->SetStatH(0.13);
  gStyle->SetStatY(1.0);

  int nVar = fVariableNames.size();

  TCanvas *c1 = new TCanvas("c1", "PullRMS",1600,1600);
  c1->Divide(4,6);
  TH1D *hPseudoUnfPullResults[nVar];

  TCanvas *c2 = new TCanvas("c2", "Pulls",1600,1600);
  c2->Divide(4,6);
  TH1D *hPseudoUnfPulls[nVar];

  TCanvas *c3 = new TCanvas("c3", "PullWidths",1600,1600);
  c3->Divide(4,6);
  TH1D *hPseudoUnfPullWidths[nVar];

  TCanvas *c4 = new TCanvas("c4", "Means",1600,1600);
  c4->Divide(4,6);
  TH1D *hPseudoUnfMeans[nVar];

  TCanvas *c5 = new TCanvas("c5", "Widths",1600,1600);
  c5->Divide(4,6);
  TH1D *hPseudoUnfWidths[nVar];

  TCanvas *c6 = new TCanvas("c6", "PullRMS_rebinnedA",1600,1600);
  c6->Divide(4,6);
  TH1D *hPseudoUnfrebinnedAPullResults[nVar];

  TCanvas *c7 = new TCanvas("c7", "Pulls_rebinnedA",1600,1600);
  c7->Divide(4,6);
  TH1D *hPseudoUnfrebinnedAPulls[nVar];

  TCanvas *c8 = new TCanvas("c8", "PullWidths_rebinnedA",1600,1600);
  c8->Divide(4,6);
  TH1D *hPseudoUnfrebinnedAPullWidths[nVar];


  TCanvas *c9 = new TCanvas("c9", "AfbPullRMS_rebinnedA",1600,1600);
  c9->Divide(4,6);
  TH1D *hPseudoUnfrebinnedAAfbPullResults[nVar];

  TCanvas *c10 = new TCanvas("c10", "AfbPulls_rebinnedA",1600,1600);
  c10->Divide(4,6);
  TH1D *hPseudoUnfrebinnedAAfbPulls[nVar];

  TCanvas *c11 = new TCanvas("c11", "AfbPullWidths_rebinnedA",1600,1600);
  c11->Divide(4,6);
  TH1D *hPseudoUnfrebinnedAAfbPullWidths[nVar];



  TCanvas *c12 = new TCanvas("c12", "PullRMS_rebinnedB",1600,1600);
  c12->Divide(4,6);
  TH1D *hPseudoUnfrebinnedBPullResults[nVar];

  TCanvas *c13 = new TCanvas("c13", "Pulls_rebinnedB",1600,1600);
  c13->Divide(4,6);
  TH1D *hPseudoUnfrebinnedBPulls[nVar];

  TCanvas *c14 = new TCanvas("c14", "PullWidths_rebinnedB",1600,1600);
  c14->Divide(4,6);
  TH1D *hPseudoUnfrebinnedBPullWidths[nVar];


  TCanvas *c15 = new TCanvas("c15", "AfbPullRMS_rebinnedB",1600,1600);
  c15->Divide(4,6);
  TH1D *hPseudoUnfrebinnedBAfbPullResults[nVar];

  TCanvas *c16 = new TCanvas("c16", "AfbPulls_rebinnedB",1600,1600);
  c16->Divide(4,6);
  TH1D *hPseudoUnfrebinnedBAfbPulls[nVar];

  TCanvas *c17 = new TCanvas("c17", "AfbPullWidths_rebinnedB",1600,1600);
  c17->Divide(4,6);
  TH1D *hPseudoUnfrebinnedBAfbPullWidths[nVar];


  //TF1 *Chi2Func = new TF1("Chi2Func", "[0]*exp(-0.5*x)*x**(0.5*[1]-1)",0,200);
  //gamma distribution
  TF1 *Chi2Func = new TF1("Chi2Func", "[0]*exp(-x/[2])*x**([1]/[2]-1)",0,200);

  TCanvas *cChi2Unf_rebinnedA = new TCanvas("cChi2Unf_rebinnedA", "Chi2Unf_rebinnedA",1600,1600);
  cChi2Unf_rebinnedA->Divide(4,6);

  TCanvas *cChi2Unf_rebinnedB = new TCanvas("cChi2Unf_rebinnedB", "Chi2Unf_rebinnedB",1600,1600);
  cChi2Unf_rebinnedB->Divide(4,6);

  TCanvas *cChi2RB_rebinnedA = new TCanvas("cChi2RB_rebinnedA", "Chi2RB_rebinnedA",1600,1600);
  cChi2RB_rebinnedA->Divide(4,6);

  TCanvas *cChi2RB_rebinnedB = new TCanvas("cChi2RB_rebinnedB", "Chi2RB_rebinnedB",1600,1600);
  cChi2RB_rebinnedB->Divide(4,6);

  TCanvas *cChi2RB2_rebinnedA = new TCanvas("cChi2RB2_rebinnedA", "Chi2RB2_rebinnedA",1600,1600);
  cChi2RB2_rebinnedA->Divide(4,6);

  TCanvas *cChi2RB2_rebinnedB = new TCanvas("cChi2RB2_rebinnedB", "Chi2RB2_rebinnedB",1600,1600);
  cChi2RB2_rebinnedB->Divide(4,6);



  TH1::AddDirectory(kFALSE);

  std::cout << "Loading histos" << std::endl;
  for(unsigned int i=0; i < fVariableNames.size(); i++){
    if(fVariableNames[i]=="") continue;
    // Amandeep : Revisit this one
    if(i>=24) continue;

    TString quantity     = MatrixUnf::ObsName(fVariableNames[i],kFALSE);
    TString quantitycoef = MatrixUnf::CoefName(fVariableNames[i],kFALSE);

    TString TUnfoldResultsFilename = Form(TUNFFILE+fVariableNames[i]+".root");
    //TString Filename = "TUnfoldResults/"+SYST+"/"+CHAN+"/"+fVariableNames[i]+"_constrainedtau0_regMode3_closure1_selfconsistency1_minrhoavg1.root";
    printf("Open file root-File %s\n\n", TUnfoldResultsFilename.Data());
    TFile *File = new TFile(TUnfoldResultsFilename.Data(),"READ");                                                                                 
    if(File==NULL) exit(0);    

    c1->cd(i+1);
    //hPseudoUnfPullResults[i] = (TH1D*) ((TH1D*)File->Get(Form("hPseudoUnfPullResults_%s",fVariableNames[i].Data())))->Clone();
    hPseudoUnfPullResults[i] = (TH1D*)File->Get(Form("hPseudoUnfPullResults_%s",fVariableNames[i].Data()));
    //hPseudoUnfPullResults[i]->SetDirectory(0); // "detach" the histogram from the file
    hPseudoUnfPullResults[i]->GetYaxis()->SetTitle("Pull mean and width");
    hPseudoUnfPullResults[i]->GetYaxis()->SetTitleOffset(1.5);
    hPseudoUnfPullResults[i]->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    //hPseudoUnfPullResults[i]->Fit("pol0");
    hPseudoUnfPullResults[i]->Draw();

    c2->cd(i+1);
    hPseudoUnfPulls[i] = (TH1D*)File->Get(Form("hPseudoUnfPulls_%s",fVariableNames[i].Data()));
    //hPseudoUnfPulls[i]->SetDirectory(0); // "detach" the histogram from the file
    hPseudoUnfPulls[i]->GetYaxis()->SetTitle("Pull mean");
    hPseudoUnfPulls[i]->GetYaxis()->SetTitleOffset(1.5);
    hPseudoUnfPulls[i]->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    hPseudoUnfPulls[i]->Draw();
    hPseudoUnfPulls[i]->Fit("pol0");

    c3->cd(i+1);
    hPseudoUnfPullWidths[i] = (TH1D*)File->Get(Form("hPseudoUnfPullWidths_%s",fVariableNames[i].Data()));
    //hPseudoUnfPullWidths[i]->SetDirectory(0); // "detach" the histogram from the file
    hPseudoUnfPullWidths[i]->GetYaxis()->SetTitle("Pull width");
    hPseudoUnfPullWidths[i]->GetYaxis()->SetTitleOffset(1.5);
    hPseudoUnfPullWidths[i]->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    hPseudoUnfPullWidths[i]->SetMinimum(0.99*hPseudoUnfPullWidths[i]->GetMinimum(0.1));
    //hPseudoUnfPullWidths[i]->SetMaximum(1.01*hPseudoUnfPullWidths[i]->GetMaximum());
    hPseudoUnfPullWidths[i]->Draw();
    hPseudoUnfPullWidths[i]->Fit("pol0");

    c4->cd(i+1);
    hPseudoUnfMeans[i] = (TH1D*)File->Get(Form("hPseudoUnfMeans_%s",fVariableNames[i].Data()));
    //hPseudoUnfMeans[i]->SetDirectory(0); // "detach" the histogram from the file
    hPseudoUnfMeans[i]->GetYaxis()->SetTitle("Mean value");
    hPseudoUnfMeans[i]->GetYaxis()->SetTitleOffset(1.5);
    hPseudoUnfMeans[i]->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    hPseudoUnfMeans[i]->Draw();

    c5->cd(i+1);
    hPseudoUnfWidths[i] = (TH1D*)File->Get(Form("hPseudoUnfWidths_%s",fVariableNames[i].Data()));
    //hPseudoUnfWidths[i]->SetDirectory(0); // "detach" the histogram from the file
    hPseudoUnfWidths[i]->GetYaxis()->SetTitle("RMS value");
    hPseudoUnfWidths[i]->GetYaxis()->SetTitleOffset(1.5);
    hPseudoUnfWidths[i]->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    hPseudoUnfWidths[i]->Draw();


    c6->cd(i+1);
    hPseudoUnfrebinnedAPullResults[i] = (TH1D*)File->Get(Form("hPseudoUnfrebinnedAPullResults_%s",fVariableNames[i].Data()));
    //hPseudoUnfrebinnedAPullResults[i]->SetDirectory(0); // "detach" the histogram from the file
    hPseudoUnfrebinnedAPullResults[i]->GetYaxis()->SetTitle("Pull mean and width (rebinnedA)");
    hPseudoUnfrebinnedAPullResults[i]->GetYaxis()->SetTitleOffset(1.5);
    hPseudoUnfrebinnedAPullResults[i]->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    //hPseudoUnfrebinnedAPullResults[i]->Fit("pol0");
    hPseudoUnfrebinnedAPullResults[i]->Draw();

    c7->cd(i+1);
    hPseudoUnfrebinnedAPulls[i] = (TH1D*)File->Get(Form("hPseudoUnfrebinnedAPulls_%s",fVariableNames[i].Data()));
    //hPseudoUnfrebinnedAPulls[i]->SetDirectory(0); // "detach" the histogram from the file
    hPseudoUnfrebinnedAPulls[i]->GetYaxis()->SetTitle("Pull mean (rebinnedA)");
    hPseudoUnfrebinnedAPulls[i]->GetYaxis()->SetTitleOffset(1.5);
    hPseudoUnfrebinnedAPulls[i]->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    hPseudoUnfrebinnedAPulls[i]->Draw();
    hPseudoUnfrebinnedAPulls[i]->Fit("pol0");

    c8->cd(i+1);
    hPseudoUnfrebinnedAPullWidths[i] = (TH1D*)File->Get(Form("hPseudoUnfrebinnedAPullWidths_%s",fVariableNames[i].Data()));
    //hPseudoUnfrebinnedAPullWidths[i]->SetDirectory(0); // "detach" the histogram from the file
    hPseudoUnfrebinnedAPullWidths[i]->GetYaxis()->SetTitle("Pull width (rebinnedA)");
    hPseudoUnfrebinnedAPullWidths[i]->GetYaxis()->SetTitleOffset(1.5);
    hPseudoUnfrebinnedAPullWidths[i]->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    hPseudoUnfrebinnedAPullWidths[i]->SetMinimum(0.99*hPseudoUnfrebinnedAPullWidths[i]->GetMinimum(0.1));
    //hPseudoUnfrebinnedAPullWidths[i]->SetMaximum(1.01*hPseudoUnfrebinnedAPullWidths[i]->GetMaximum());
    hPseudoUnfrebinnedAPullWidths[i]->Draw();
    hPseudoUnfrebinnedAPullWidths[i]->Fit("pol0");

    c9->cd(i+1);
    hPseudoUnfrebinnedAAfbPullResults[i] = (TH1D*)File->Get(Form("hPseudoUnfrebinnedAAfbPullResults_%s",fVariableNames[i].Data()));
    //hPseudoUnfrebinnedAAfbPullResults[i]->SetDirectory(0); // "detach" the histogram from the file
    hPseudoUnfrebinnedAAfbPullResults[i]->GetYaxis()->SetTitle("Asymmetry pull mean and width (rebinnedA)");
    hPseudoUnfrebinnedAAfbPullResults[i]->GetYaxis()->SetTitleOffset(1.5);
    hPseudoUnfrebinnedAAfbPullResults[i]->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, A_{total} (right)",quantity.Data(),quantitycoef.Data()));
    //hPseudoUnfrebinnedAAfbPullResults[i]->Fit("pol0");
    hPseudoUnfrebinnedAAfbPullResults[i]->Draw();

    c10->cd(i+1);
    hPseudoUnfrebinnedAAfbPulls[i] = (TH1D*)File->Get(Form("hPseudoUnfrebinnedAAfbPulls_%s",fVariableNames[i].Data()));
    //hPseudoUnfrebinnedAAfbPulls[i]->SetDirectory(0); // "detach" the histogram from the file
    hPseudoUnfrebinnedAAfbPulls[i]->GetYaxis()->SetTitle("Asymmetry pull mean (rebinnedA)");
    hPseudoUnfrebinnedAAfbPulls[i]->GetYaxis()->SetTitleOffset(1.5);
    hPseudoUnfrebinnedAAfbPulls[i]->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, A_{total} (right)",quantity.Data(),quantitycoef.Data()));
    hPseudoUnfrebinnedAAfbPulls[i]->Draw();
    hPseudoUnfrebinnedAAfbPulls[i]->Fit("pol0");

    c11->cd(i+1);
    hPseudoUnfrebinnedAAfbPullWidths[i] = (TH1D*)File->Get(Form("hPseudoUnfrebinnedAAfbPullWidths_%s",fVariableNames[i].Data()));
    //hPseudoUnfrebinnedAAfbPullWidths[i]->SetDirectory(0); // "detach" the histogram from the file
    hPseudoUnfrebinnedAAfbPullWidths[i]->GetYaxis()->SetTitle("Asymmetry pull width (rebinnedA)");
    hPseudoUnfrebinnedAAfbPullWidths[i]->GetYaxis()->SetTitleOffset(1.5);
    hPseudoUnfrebinnedAAfbPullWidths[i]->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, A_{total} (right)",quantity.Data(),quantitycoef.Data()));
    hPseudoUnfrebinnedAAfbPullWidths[i]->SetMinimum(0.99*hPseudoUnfrebinnedAAfbPullWidths[i]->GetMinimum(0.1));
    //hPseudoUnfrebinnedAAfbPullWidths[i]->SetMaximum(1.01*hPseudoUnfrebinnedAAfbPullWidths[i]->GetMaximum());
    hPseudoUnfrebinnedAAfbPullWidths[i]->Draw();
    hPseudoUnfrebinnedAAfbPullWidths[i]->Fit("pol0");


    c12->cd(i+1);
    hPseudoUnfrebinnedBPullResults[i] = (TH1D*)File->Get(Form("hPseudoUnfrebinnedBPullResults_%s",fVariableNames[i].Data()));
    //hPseudoUnfrebinnedBPullResults[i]->SetDirectory(0); // "detach" the histogram from the file
    hPseudoUnfrebinnedBPullResults[i]->GetYaxis()->SetTitle("Pull mean and width (rebinnedB)");
    hPseudoUnfrebinnedBPullResults[i]->GetYaxis()->SetTitleOffset(1.5);
    hPseudoUnfrebinnedBPullResults[i]->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    //hPseudoUnfrebinnedBPullResults[i]->Fit("pol0");
    hPseudoUnfrebinnedBPullResults[i]->Draw();

    c13->cd(i+1);
    hPseudoUnfrebinnedBPulls[i] = (TH1D*)File->Get(Form("hPseudoUnfrebinnedBPulls_%s",fVariableNames[i].Data()));
    //hPseudoUnfrebinnedBPulls[i]->SetDirectory(0); // "detach" the histogram from the file
    hPseudoUnfrebinnedBPulls[i]->GetYaxis()->SetTitle("Pull mean (rebinnedB)");
    hPseudoUnfrebinnedBPulls[i]->GetYaxis()->SetTitleOffset(1.5);
    hPseudoUnfrebinnedBPulls[i]->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    hPseudoUnfrebinnedBPulls[i]->Draw();
    hPseudoUnfrebinnedBPulls[i]->Fit("pol0");

    c14->cd(i+1);
    hPseudoUnfrebinnedBPullWidths[i] = (TH1D*)File->Get(Form("hPseudoUnfrebinnedBPullWidths_%s",fVariableNames[i].Data()));
    //hPseudoUnfrebinnedBPullWidths[i]->SetDirectory(0); // "detach" the histogram from the file
    hPseudoUnfrebinnedBPullWidths[i]->GetYaxis()->SetTitle("Pull width (rebinnedB)");
    hPseudoUnfrebinnedBPullWidths[i]->GetYaxis()->SetTitleOffset(1.5);
    hPseudoUnfrebinnedBPullWidths[i]->GetXaxis()->SetTitle(Form("%s",quantity.Data()));
    hPseudoUnfrebinnedBPullWidths[i]->SetMinimum(0.99*hPseudoUnfrebinnedBPullWidths[i]->GetMinimum(0.1));
    //hPseudoUnfrebinnedBPullWidths[i]->SetMaximum(1.01*hPseudoUnfrebinnedBPullWidths[i]->GetMaximum());
    hPseudoUnfrebinnedBPullWidths[i]->Draw();
    hPseudoUnfrebinnedBPullWidths[i]->Fit("pol0");

    c15->cd(i+1);
    hPseudoUnfrebinnedBAfbPullResults[i] = (TH1D*)File->Get(Form("hPseudoUnfrebinnedBAfbPullResults_%s",fVariableNames[i].Data()));
    //hPseudoUnfrebinnedBAfbPullResults[i]->SetDirectory(0); // "detach" the histogram from the file
    hPseudoUnfrebinnedBAfbPullResults[i]->GetYaxis()->SetTitle("Asymmetry pull mean and width (rebinnedB)");
    hPseudoUnfrebinnedBAfbPullResults[i]->GetYaxis()->SetTitleOffset(1.5);
    hPseudoUnfrebinnedBAfbPullResults[i]->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, A_{total} (right)",quantity.Data(),quantitycoef.Data()));
    //hPseudoUnfrebinnedBAfbPullResults[i]->Fit("pol0");
    hPseudoUnfrebinnedBAfbPullResults[i]->Draw();

    c16->cd(i+1);
    hPseudoUnfrebinnedBAfbPulls[i] = (TH1D*)File->Get(Form("hPseudoUnfrebinnedBAfbPulls_%s",fVariableNames[i].Data()));
    //hPseudoUnfrebinnedBAfbPulls[i]->SetDirectory(0); // "detach" the histogram from the file
    hPseudoUnfrebinnedBAfbPulls[i]->GetYaxis()->SetTitle("Asymmetry pull mean (rebinnedB)");
    hPseudoUnfrebinnedBAfbPulls[i]->GetYaxis()->SetTitleOffset(1.5);
    hPseudoUnfrebinnedBAfbPulls[i]->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, A_{total} (right)",quantity.Data(),quantitycoef.Data()));
    hPseudoUnfrebinnedBAfbPulls[i]->Draw();
    hPseudoUnfrebinnedBAfbPulls[i]->Fit("pol0");

    c17->cd(i+1);
    hPseudoUnfrebinnedBAfbPullWidths[i] = (TH1D*)File->Get(Form("hPseudoUnfrebinnedBAfbPullWidths_%s",fVariableNames[i].Data()));
    //hPseudoUnfrebinnedBAfbPullWidths[i]->SetDirectory(0); // "detach" the histogram from the file
    hPseudoUnfrebinnedBAfbPullWidths[i]->GetYaxis()->SetTitle("Asymmetry pull width (rebinnedB)");
    hPseudoUnfrebinnedBAfbPullWidths[i]->GetYaxis()->SetTitleOffset(1.5);
    hPseudoUnfrebinnedBAfbPullWidths[i]->GetXaxis()->SetTitle(Form("A(|%s|) (left); bin sum, %s, A_{total} (right)",quantity.Data(),quantitycoef.Data()));
    hPseudoUnfrebinnedBAfbPullWidths[i]->SetMinimum(0.99*hPseudoUnfrebinnedBAfbPullWidths[i]->GetMinimum(0.1));
    //hPseudoUnfrebinnedBAfbPullWidths[i]->SetMaximum(1.01*hPseudoUnfrebinnedBAfbPullWidths[i]->GetMaximum());
    hPseudoUnfrebinnedBAfbPullWidths[i]->Draw();
    hPseudoUnfrebinnedBAfbPullWidths[i]->Fit("pol0");



    double p0,p1,p2;

    cChi2Unf_rebinnedA->cd(i+1);
    TH1D* hChi2UnfRebinnedA = (TH1D*)File->Get(Form("hPseudoUnfrebinnedAOtherPerBin_5"));
    hChi2UnfRebinnedA->GetYaxis()->SetTitle("Chi2Unf (rebinnedA)");
    hChi2UnfRebinnedA->GetYaxis()->SetTitleOffset(1.5);
    hChi2UnfRebinnedA->GetXaxis()->SetTitle(Form("#chi^2 (unf) for %s",quantity.Data()));
    p1=hChi2UnfRebinnedA->GetMean();
    p2=pow(hChi2UnfRebinnedA->GetRMS(),2)/hChi2UnfRebinnedA->GetMean();
    p2=2;
    p0=hChi2UnfRebinnedA->GetEntries()/(pow(p2,p1/p2)*TMath::Gamma(p1/p2));
    Chi2Func->SetParameter(0,p0);
    Chi2Func->SetParLimits(0,p0/3.3,p0*3.3);
    Chi2Func->SetParameter(1,p1);
    Chi2Func->SetParLimits(1,p1/3.3,p1*3.3);
    Chi2Func->FixParameter(2,2);
    //Chi2Func->SetParLimits(2,2/1.5,2*1.5);
    hChi2UnfRebinnedA->Fit(Chi2Func,"PEMRL");
    hChi2UnfRebinnedA->Draw();
    std::cout<<"Chi2Unf (rebinnedA) for "<<fVariableNames[i]<<": "<<hChi2UnfRebinnedA->GetMean()<<" +/- "<<hChi2UnfRebinnedA->GetMeanError()<<" , estimated ndf in unfolded distribution: "<<48-1-hChi2UnfRebinnedA->GetMean()<<std::endl;
    //std::cout<<"Chi2Unf (rebinnedA) for "<<fVariableNames[i]<<": "<<hChi2UnfRebinnedA->GetMean()<<" +/- "<<hChi2UnfRebinnedA->GetMeanError()<<" Fitted: "<<Chi2Func->GetParameter(1)<<" +/- "<<Chi2Func->GetParError(1)<<" , estimated ndf in unfolded distribution: "<<48-1-Chi2Func->GetParameter(1)<<std::endl;


    cChi2Unf_rebinnedB->cd(i+1);
    TH1D* hChi2UnfRebinnedB = (TH1D*)File->Get(Form("hPseudoUnfrebinnedBOtherPerBin_5"));
    hChi2UnfRebinnedB->GetYaxis()->SetTitle("Chi2Unf (rebinnedB)");
    hChi2UnfRebinnedB->GetYaxis()->SetTitleOffset(1.5);
    hChi2UnfRebinnedB->GetXaxis()->SetTitle(Form("#chi^2 (unf) for %s",quantity.Data()));
    p1=hChi2UnfRebinnedB->GetMean();
    p2=pow(hChi2UnfRebinnedB->GetRMS(),2)/hChi2UnfRebinnedB->GetMean();
    p2=2;
    p0=hChi2UnfRebinnedB->GetEntries()/(pow(p2,p1/p2)*TMath::Gamma(p1/p2));
    Chi2Func->SetParameter(0,p0);
    Chi2Func->SetParLimits(0,p0/3.3,p0*3.3);
    Chi2Func->SetParameter(1,p1);
    Chi2Func->SetParLimits(1,p1/3.3,p1*3.3);
    Chi2Func->FixParameter(2,2);
    //Chi2Func->SetParLimits(2,2/1.5,2*1.5);
    hChi2UnfRebinnedB->Fit(Chi2Func,"PEMRL");
    hChi2UnfRebinnedB->Draw();
    std::cout<<"Chi2Unf (rebinnedB) for "<<fVariableNames[i]<<": "<<hChi2UnfRebinnedB->GetMean()<<" +/- "<<hChi2UnfRebinnedB->GetMeanError()<<" , estimated ndf in unfolded distribution: "<<12-1-hChi2UnfRebinnedB->GetMean()<<std::endl;
    //std::cout<<"Chi2Unf (rebinnedB) for "<<fVariableNames[i]<<": "<<hChi2UnfRebinnedB->GetMean()<<" +/- "<<hChi2UnfRebinnedB->GetMeanError()<<" Fitted: "<<Chi2Func->GetParameter(1)<<" +/- "<<Chi2Func->GetParError(1)<<" , estimated ndf in unfolded distribution: "<<12-1-Chi2Func->GetParameter(1)<<std::endl;

    std::cout<<"dndf(Chi2Unf) for "<<fVariableNames[i]<<": "<<(48-1-hChi2UnfRebinnedA->GetMean()) - (12-1-hChi2UnfRebinnedB->GetMean())<<std::endl;



    cChi2RB_rebinnedA->cd(i+1);
    TH1D* hChi2RBRebinnedA = (TH1D*)File->Get(Form("hPseudoUnfrebinnedAOtherPerBin_7"));
    hChi2RBRebinnedA->GetYaxis()->SetTitle("Chi2RB (rebinnedA)");
    hChi2RBRebinnedA->GetYaxis()->SetTitleOffset(1.5);
    hChi2RBRebinnedA->GetXaxis()->SetTitle(Form("#chi^2 (RB) for %s",quantity.Data()));
    p1=hChi2RBRebinnedA->GetMean();
    p2=pow(hChi2RBRebinnedA->GetRMS(),2)/hChi2RBRebinnedA->GetMean();
    p0=hChi2RBRebinnedA->GetEntries()/(pow(p2,p1/p2)*TMath::Gamma(p1/p2));
    Chi2Func->SetParameter(0,p0);
    Chi2Func->SetParLimits(0,p0/3.3,p0*3.3);
    Chi2Func->SetParameter(1,p1);
    Chi2Func->SetParLimits(1,p1/3.3,p1*3.3);
    Chi2Func->SetParameter(2,p2);
    Chi2Func->SetParLimits(2,p2/3.3,p2*3.3);
    hChi2RBRebinnedA->Fit(Chi2Func,"PEMRL");
    hChi2RBRebinnedA->Draw();
    std::cout<<"Chi2RB (rebinnedA) for "<<fVariableNames[i]<<": "<<hChi2RBRebinnedA->GetMean()<<" +/- "<<hChi2RBRebinnedA->GetMeanError()<<" , estimated ndf in unfolded distribution: "<<12-1-hChi2RBRebinnedA->GetMean()<<std::endl;
    //std::cout<<"Chi2RB (rebinnedA) for "<<fVariableNames[i]<<": "<<hChi2RBRebinnedA->GetMean()<<" +/- "<<hChi2RBRebinnedA->GetMeanError()<<" Fitted: "<<Chi2Func->GetParameter(1)<<" +/- "<<Chi2Func->GetParError(1)<<" , estimated ndf in unfolded distribution: "<<12-1-Chi2Func->GetParameter(1)<<std::endl;


    cChi2RB_rebinnedB->cd(i+1);
    TH1D* hChi2RBRebinnedB = (TH1D*)File->Get(Form("hPseudoUnfrebinnedBOtherPerBin_7"));
    hChi2RBRebinnedB->GetYaxis()->SetTitle("Chi2RB (rebinnedB)");
    hChi2RBRebinnedB->GetYaxis()->SetTitleOffset(1.5);
    hChi2RBRebinnedB->GetXaxis()->SetTitle(Form("#chi^2 (RB) for %s",quantity.Data()));
    p1=hChi2RBRebinnedB->GetMean();
    p2=pow(hChi2RBRebinnedB->GetRMS(),2)/hChi2RBRebinnedB->GetMean();
    p0=hChi2RBRebinnedB->GetEntries()/(pow(p2,p1/p2)*TMath::Gamma(p1/p2));
    Chi2Func->SetParameter(0,p0);
    Chi2Func->SetParLimits(0,p0/3.3,p0*3.3);
    Chi2Func->SetParameter(1,p1);
    Chi2Func->SetParLimits(1,p1/3.3,p1*3.3);
    Chi2Func->SetParameter(2,p2);
    Chi2Func->SetParLimits(2,p2/3.3,p2*3.3);
    hChi2RBRebinnedB->Fit(Chi2Func,"PEMRL");
    hChi2RBRebinnedB->Draw();
    std::cout<<"Chi2RB (rebinnedB) for "<<fVariableNames[i]<<": "<<hChi2RBRebinnedB->GetMean()<<" +/- "<<hChi2RBRebinnedB->GetMeanError()<<" , estimated ndf in unfolded distribution: "<<12-1-hChi2RBRebinnedB->GetMean()<<std::endl;
    //std::cout<<"Chi2RB (rebinnedB) for "<<fVariableNames[i]<<": "<<hChi2RBRebinnedB->GetMean()<<" +/- "<<hChi2RBRebinnedB->GetMeanError()<<" Fitted: "<<Chi2Func->GetParameter(1)<<" +/- "<<Chi2Func->GetParError(1)<<" , estimated ndf in unfolded distribution: "<<12-1-Chi2Func->GetParameter(1)<<std::endl;

    std::cout<<"dndf(Chi2RB) for "<<fVariableNames[i]<<": "<<(12-1-hChi2RBRebinnedA->GetMean()) - (12-1-hChi2RBRebinnedB->GetMean())<<std::endl;




    cChi2RB2_rebinnedA->cd(i+1);
    TH1D* hChi2RB2RebinnedA = (TH1D*)File->Get(Form("hPseudoUnfrebinnedAOtherPerBin_6"));
    hChi2RB2RebinnedA->GetYaxis()->SetTitle("Chi2RB2 (rebinnedA)");
    hChi2RB2RebinnedA->GetYaxis()->SetTitleOffset(1.5);
    hChi2RB2RebinnedA->GetXaxis()->SetTitle(Form("#chi^2 (RB2) for %s",quantity.Data()));
    p1=hChi2RB2RebinnedA->GetMean();
    p2=pow(hChi2RB2RebinnedA->GetRMS(),2)/hChi2RB2RebinnedA->GetMean();
    p0=hChi2RB2RebinnedA->GetEntries()/(pow(p2,p1/p2)*TMath::Gamma(p1/p2));
    Chi2Func->SetParameter(0,p0);
    Chi2Func->SetParLimits(0,p0/10,p0*10);
    Chi2Func->SetParameter(1,p1);
    Chi2Func->SetParLimits(1,p1/10,p1*10);
    Chi2Func->SetParameter(2,p2);
    Chi2Func->SetParLimits(2,p2/10,p2*10);
    //hChi2RB2RebinnedA->Fit(Chi2Func,"PEMRL");
    hChi2RB2RebinnedA->Draw();
    std::cout<<"Chi2RB2 (rebinnedA) for "<<fVariableNames[i]<<": "<<hChi2RB2RebinnedA->GetMean()<<" +/- "<<hChi2RB2RebinnedA->GetMeanError()<<" , estimated ndf in unfolded distribution: "<<6-1-hChi2RB2RebinnedA->GetMean()<<std::endl;
    //std::cout<<"Chi2RB2 (rebinnedA) for "<<fVariableNames[i]<<": "<<hChi2RB2RebinnedA->GetMean()<<" +/- "<<hChi2RB2RebinnedA->GetMeanError()<<" Fitted: "<<Chi2Func->GetParameter(1)<<" +/- "<<Chi2Func->GetParError(1)<<" , estimated ndf in unfolded distribution: "<<6-1-Chi2Func->GetParameter(1)<<std::endl;


    cChi2RB2_rebinnedB->cd(i+1);
    TH1D* hChi2RB2RebinnedB = (TH1D*)File->Get(Form("hPseudoUnfrebinnedBOtherPerBin_6"));
    hChi2RB2RebinnedB->GetYaxis()->SetTitle("Chi2RB2 (rebinnedB)");
    hChi2RB2RebinnedB->GetYaxis()->SetTitleOffset(1.5);
    hChi2RB2RebinnedB->GetXaxis()->SetTitle(Form("#chi^2 (RB2) for %s",quantity.Data()));
    p1=hChi2RB2RebinnedB->GetMean();
    p2=pow(hChi2RB2RebinnedB->GetRMS(),2)/hChi2RB2RebinnedB->GetMean();
    p0=hChi2RB2RebinnedB->GetEntries()/(pow(p2,p1/p2)*TMath::Gamma(p1/p2));
    Chi2Func->SetParameter(0,p0);
    Chi2Func->SetParLimits(0,p0/10,p0*10);
    Chi2Func->SetParameter(1,p1);
    Chi2Func->SetParLimits(1,p1/10,p1*10);
    Chi2Func->SetParameter(2,p2);
    Chi2Func->SetParLimits(2,p2/10,p2*10);
    // hChi2RB2RebinnedB->Fit(Chi2Func,"PEMRL");
    hChi2RB2RebinnedB->Draw();

    // Amandeep : Hardcoded here
    std::cout<<"Chi2RB2 (rebinnedB) for "<<fVariableNames[i]<<": "<<hChi2RB2RebinnedB->GetMean()<<" +/- "<<hChi2RB2RebinnedB->GetMeanError()<<" , estimated ndf in unfolded distribution: "<<6-1-hChi2RB2RebinnedB->GetMean()<<std::endl;
    // std::cout<<"Chi2RB2 (rebinnedB) for "<<fVariableNames[i]<<": "<<hChi2RB2RebinnedB->GetMean()<<" +/- "<<hChi2RB2RebinnedB->GetMeanError()<<" Fitted: "<<Chi2Func->GetParameter(1)<<" +/- "<<Chi2Func->GetParError(1)<<" , estimated ndf in unfolded distribution: "<<6-1-Chi2Func->GetParameter(1)<<std::endl;
    // Amandeep : Hardcoded here
    std::cout<<"dndf(Chi2RB2) for "<<fVariableNames[i]<<": "<<(6-1-hChi2RB2RebinnedA->GetMean()) - (6-1-hChi2RB2RebinnedB->GetMean())<<std::endl;


    double pullwidthmean=0.;
    double pullmean=0.;
    double bincount=0.;
    //hPseudoUnfPullResults[i]->Print("all");
    for(int i_bin=1;i_bin<=hPseudoUnfPullResults[i]->GetNbinsX();i_bin++){
      if(hPseudoUnfPullResults[i]->GetBinContent(i_bin) != 0) {
        pullmean+=hPseudoUnfPullResults[i]->GetBinContent(i_bin);
        pullwidthmean+=hPseudoUnfPullResults[i]->GetBinError(i_bin);
        bincount+=1.;
      }
    }
    pullwidthmean/=bincount;
    pullmean/=bincount;
    std::cout<<"Average pull and pull width for "<<fVariableNames[i]<<": "<<pullmean<<" , "<<pullwidthmean<<std::endl;


    pullwidthmean=0.;
    pullmean=0.;
    bincount=0.;
    for(int i_bin=1;i_bin<=hPseudoUnfrebinnedAPullResults[i]->GetNbinsX();i_bin++){
      if(hPseudoUnfrebinnedAPullResults[i]->GetBinContent(i_bin) != 0) {
        pullmean+=hPseudoUnfrebinnedAPullResults[i]->GetBinContent(i_bin);
        pullwidthmean+=hPseudoUnfrebinnedAPullResults[i]->GetBinError(i_bin);
        bincount+=1.;
      }
    }
    pullwidthmean/=bincount;
    pullmean/=bincount;
    std::cout<<"Average pull and pull width (rebinnedA) for "<<fVariableNames[i]<<": "<<pullmean<<" , "<<pullwidthmean<<std::endl;


    pullwidthmean=0.;
    pullmean=0.;
    bincount=0.;
    for(int i_bin=1;i_bin<=hPseudoUnfrebinnedAAfbPullResults[i]->GetNbinsX();i_bin++){
      if(hPseudoUnfrebinnedAAfbPullResults[i]->GetBinContent(i_bin) != 0) {
        pullmean+=hPseudoUnfrebinnedAAfbPullResults[i]->GetBinContent(i_bin);
        pullwidthmean+=hPseudoUnfrebinnedAAfbPullResults[i]->GetBinError(i_bin);
        bincount+=1.;
      }
    }
    pullwidthmean/=bincount;
    pullmean/=bincount;
    std::cout<<"Average Afb pull and pull width (rebinnedA) for "<<fVariableNames[i]<<": "<<pullmean<<" , "<<pullwidthmean<<std::endl;




    pullwidthmean=0.;
    pullmean=0.;
    bincount=0.;
    for(int i_bin=1;i_bin<=hPseudoUnfrebinnedBPullResults[i]->GetNbinsX();i_bin++){
      if(hPseudoUnfrebinnedBPullResults[i]->GetBinContent(i_bin) != 0) {
        pullmean+=hPseudoUnfrebinnedBPullResults[i]->GetBinContent(i_bin);
        pullwidthmean+=hPseudoUnfrebinnedBPullResults[i]->GetBinError(i_bin);
        bincount+=1.;
      }
    }
    pullwidthmean/=bincount;
    pullmean/=bincount;
    std::cout<<"Average pull and pull width (rebinnedB) for "<<fVariableNames[i]<<": "<<pullmean<<" , "<<pullwidthmean<<std::endl;



    pullwidthmean=0.;
    pullmean=0.;
    bincount=0.;
    for(int i_bin=1;i_bin<=hPseudoUnfrebinnedBAfbPullResults[i]->GetNbinsX();i_bin++){
      if(hPseudoUnfrebinnedBAfbPullResults[i]->GetBinContent(i_bin) != 0) {
        pullmean+=hPseudoUnfrebinnedBAfbPullResults[i]->GetBinContent(i_bin);
        pullwidthmean+=hPseudoUnfrebinnedBAfbPullResults[i]->GetBinError(i_bin);
        bincount+=1.;
      }
    }
    pullwidthmean/=bincount;
    pullmean/=bincount;
    std::cout<<"Average Afb pull and pull width (rebinnedB) for "<<fVariableNames[i]<<": "<<pullmean<<" , "<<pullwidthmean<<std::endl;



    if(true) {
      TCanvas *cpull = new TCanvas("cpull", "Pulldist",1600,1200);
      cpull->Divide(3,2);
      for(int i_bin=1;i_bin<=hPseudoUnfrebinnedAPullResults[i]->GetNbinsX();i_bin++)
      {
        cpull->cd(i_bin);
        TH1D* hist = (TH1D*)File->Get(Form("hPseudoUnfrebinnedAPullPerBin_%d",i_bin));
        hist->GetYaxis()->SetTitle("Number of PEs");
        hist->GetYaxis()->SetTitleOffset(1.6);
        hist->GetXaxis()->SetTitle(Form("%s bin %d pull",quantity.Data(),i_bin));
        hist->Draw();
      }
      cpull->Print(Form("%s/PullDist_rebinnedA_%s.pdf",PLOTFILE.Data(),fVariableNames[i].Data()));
      cpull->Print(Form("%s/PullDist_rebinnedA_%s.root",PLOTFILE.Data(),fVariableNames[i].Data()));
      cpull->Print(Form("%s/PullDist_rebinnedA_%s.C",PLOTFILE.Data(),fVariableNames[i].Data()));
      delete cpull;
    }





    File->Close();
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  gStyle->SetStatW(0.26);
  gStyle->SetStatH(0.13);
  gStyle->SetStatY(1.0);

  c1->Print(Form("%s/PullsandWidths.pdf",PLOTFILE.Data()));
  c1->Print(Form("%s/PullsandWidths.root",PLOTFILE.Data()));
  c1->Print(Form("%s/PullsandWidths.C",PLOTFILE.Data()));
  c1->Close();

  c2->Print(Form("%s/PullMean.pdf",PLOTFILE.Data()));
  c2->Print(Form("%s/PullMean.root",PLOTFILE.Data()));
  c2->Print(Form("%s/PullMean.C",PLOTFILE.Data()));
  c2->Close();

  c3->Print(Form("%s/PullWidth.pdf",PLOTFILE.Data()));
  c3->Print(Form("%s/PullWidth.root",PLOTFILE.Data()));
  c3->Print(Form("%s/PullWidth.C",PLOTFILE.Data()));
  c3->Close();

  c4->Print(Form("%s/ValueMean.pdf",PLOTFILE.Data()));
  c4->Print(Form("%s/ValueMean.root",PLOTFILE.Data()));
  c4->Print(Form("%s/ValueMean.C",PLOTFILE.Data()));
  c4->Close();

  c5->Print(Form("%s/ValueRMS.pdf",PLOTFILE.Data()));
  c5->Print(Form("%s/ValueRMS.root",PLOTFILE.Data()));
  c5->Print(Form("%s/ValueRMS.C",PLOTFILE.Data()));
  c5->Close();

  c6->Print(Form("%s/PullsandWidthsrebinnedA.pdf",PLOTFILE.Data()));
  c6->Print(Form("%s/PullsandWidthsrebinnedA.root",PLOTFILE.Data()));
  c6->Print(Form("%s/PullsandWidthsrebinnedA.C",PLOTFILE.Data()));
  c6->Close();

  c7->Print(Form("%s/PullMeanrebinnedA.pdf",PLOTFILE.Data()));
  c7->Print(Form("%s/PullMeanrebinnedA.root",PLOTFILE.Data()));
  c7->Print(Form("%s/PullMeanrebinnedA.C",PLOTFILE.Data()));
  c7->Close();

  c8->Print(Form("%s/PullWidthrebinnedA.pdf",PLOTFILE.Data()));
  c8->Print(Form("%s/PullWidthrebinnedA.root",PLOTFILE.Data()));
  c8->Print(Form("%s/PullWidthrebinnedA.C",PLOTFILE.Data()));
  c8->Close();

  c9->Print(Form("%s/AfbPullsandWidthsrebinnedA.pdf",PLOTFILE.Data()));
  c9->Print(Form("%s/AfbPullsandWidthsrebinnedA.root",PLOTFILE.Data()));
  c9->Print(Form("%s/AfbPullsandWidthsrebinnedA.C",PLOTFILE.Data()));
  c9->Close();

  c10->Print(Form("%s/AfbPullMeanrebinnedA.pdf",PLOTFILE.Data()));
  c10->Print(Form("%s/AfbPullMeanrebinnedA.root",PLOTFILE.Data()));
  c10->Print(Form("%s/AfbPullMeanrebinnedA.C",PLOTFILE.Data()));
  c10->Close();

  c11->Print(Form("%s/AfbPullWidthrebinnedA.pdf",PLOTFILE.Data()));
  c11->Print(Form("%s/AfbPullWidthrebinnedA.root",PLOTFILE.Data()));
  c11->Print(Form("%s/AfbPullWidthrebinnedA.C",PLOTFILE.Data()));
  c11->Close();

  c12->Print(Form("%s/PullsandWidthsrebinnedB.pdf",PLOTFILE.Data()));
  c12->Print(Form("%s/PullsandWidthsrebinnedB.root",PLOTFILE.Data()));
  c12->Print(Form("%s/PullsandWidthsrebinnedB.C",PLOTFILE.Data()));
  c12->Close();

  c13->Print(Form("%s/PullMeanrebinnedB.pdf",PLOTFILE.Data()));
  c13->Print(Form("%s/PullMeanrebinnedB.root",PLOTFILE.Data()));
  c13->Print(Form("%s/PullMeanrebinnedB.C",PLOTFILE.Data()));
  c13->Close();

  c14->Print(Form("%s/PullWidthrebinnedB.pdf",PLOTFILE.Data()));
  c14->Print(Form("%s/PullWidthrebinnedB.root",PLOTFILE.Data()));
  c14->Print(Form("%s/PullWidthrebinnedB.C",PLOTFILE.Data()));
  c14->Close();

  c15->Print(Form("%s/AfbPullsandWidthsrebinnedB.pdf",PLOTFILE.Data()));
  c15->Print(Form("%s/AfbPullsandWidthsrebinnedB.root",PLOTFILE.Data()));
  c15->Print(Form("%s/AfbPullsandWidthsrebinnedB.C",PLOTFILE.Data()));
  c15->Close();

  c16->Print(Form("%s/AfbPullMeanrebinnedB.pdf",PLOTFILE.Data()));
  c16->Print(Form("%s/AfbPullMeanrebinnedB.root",PLOTFILE.Data()));
  c16->Print(Form("%s/AfbPullMeanrebinnedB.C",PLOTFILE.Data()));
  c16->Close();

  c17->Print(Form("%s/AfbPullWidthrebinnedB.pdf",PLOTFILE.Data()));
  c17->Print(Form("%s/AfbPullWidthrebinnedB.root",PLOTFILE.Data()));
  c17->Print(Form("%s/AfbPullWidthrebinnedB.C",PLOTFILE.Data()));
  c17->Close();


  cChi2Unf_rebinnedA->Print(Form("%s/Chi2Unf_rebinnedA.pdf",PLOTFILE.Data()));
  cChi2Unf_rebinnedA->Print(Form("%s/Chi2Unf_rebinnedA.root",PLOTFILE.Data()));
  cChi2Unf_rebinnedA->Print(Form("%s/Chi2Unf_rebinnedA.C",PLOTFILE.Data()));
  cChi2Unf_rebinnedA->Close();

  cChi2Unf_rebinnedB->Print(Form("%s/Chi2Unf_rebinnedB.pdf",PLOTFILE.Data()));
  cChi2Unf_rebinnedB->Print(Form("%s/Chi2Unf_rebinnedB.root",PLOTFILE.Data()));
  cChi2Unf_rebinnedB->Print(Form("%s/Chi2Unf_rebinnedB.C",PLOTFILE.Data()));
  cChi2Unf_rebinnedB->Close();


  cChi2RB_rebinnedA->Print(Form("%s/Chi2RB_rebinnedA.pdf",PLOTFILE.Data()));
  cChi2RB_rebinnedA->Print(Form("%s/Chi2RB_rebinnedA.root",PLOTFILE.Data()));
  cChi2RB_rebinnedA->Print(Form("%s/Chi2RB_rebinnedA.C",PLOTFILE.Data()));
  cChi2RB_rebinnedA->Close();

  cChi2RB_rebinnedB->Print(Form("%s/Chi2RB_rebinnedB.pdf",PLOTFILE.Data()));
  cChi2RB_rebinnedB->Print(Form("%s/Chi2RB_rebinnedB.root",PLOTFILE.Data()));
  cChi2RB_rebinnedB->Print(Form("%s/Chi2RB_rebinnedB.C",PLOTFILE.Data()));
  cChi2RB_rebinnedB->Close();

  cChi2RB2_rebinnedA->Print(Form("%s/Chi2RB2_rebinnedA.pdf",PLOTFILE.Data()));
  cChi2RB2_rebinnedA->Print(Form("%s/Chi2RB2_rebinnedA.root",PLOTFILE.Data()));
  cChi2RB2_rebinnedA->Print(Form("%s/Chi2RB2_rebinnedA.C",PLOTFILE.Data()));
  cChi2RB2_rebinnedA->Close();

  cChi2RB2_rebinnedB->Print(Form("%s/Chi2RB2_rebinnedB.pdf",PLOTFILE.Data()));
  cChi2RB2_rebinnedB->Print(Form("%s/Chi2RB2_rebinnedB.root",PLOTFILE.Data()));
  cChi2RB2_rebinnedB->Print(Form("%s/Chi2RB2_rebinnedB.C",PLOTFILE.Data()));
  cChi2RB2_rebinnedB->Close();

}


int main(int argc, char* argv[]) {


  if (argc != 13) {
    std::cout << "Wrong number of arguments provided, please call the tool with the following options:" << std::endl;
    std::cout << "./MatrixUnfControl <Era> <Variable_Name> <Systematics> <Channel> <SelfConsistency?> <UseBiasIntUnfold?> <AllowRebin?> <EnsembleTest?> <LinearityTest?> <SetSyst?> <plotoption?0:unfonly;1:unf+plots;2:resultplotsonly;3:pullplotsonly>" << std::endl;
    return 1;
  }

  std::string era(argv[1]);
  std::string variable_name(argv[2]);
  std::string systematics(argv[3]);
  std::string channel(argv[4]);

  bool selfconsistency(std::stoi(argv[5]));
  bool usebiasintunfold(std::stoi(argv[6]));
  bool allowrebin(std::stoi(argv[7]));
  bool ensembletest(std::stoi(argv[8]));
  bool linearitytest(std::stoi(argv[9]));
  bool setsyst(std::stoi(argv[10]));
  
  int bkgsuboption(std::stoi(argv[11]));
  int plotoption(std::stoi(argv[12]));

  std::cout << era << " " <<  variable_name << " " << systematics << " " << channel << " " << selfconsistency << " " << usebiasintunfold << " " << allowrebin << " " << ensembletest << " " << linearitytest << " " << setsyst << " " << bkgsuboption << " " << plotoption << std::endl;

  // Grab lumi information for each specific era
  float lumi = 0.0;
  if (era.compare("2016UL") == 0) {
    lumi = 36310.0;
  }
  else if (era.compare("2017UL") == 0) {
    lumi = 41480.0;
  }
  else if (era.compare("2018UL") == 0) {
    lumi = 59830.0;
  }
  else if (era.compare("fullRun2UL") == 0) {
    lumi = 137620.0;
  }
  else if (era.compare("2016ULpreVFP") == 0) {
    lumi = 19500.0;
  }
  else if (era.compare("2016ULpostVFP") == 0) {
    lumi = 16810.0;
  }
  else {
    std::cout << "No lumi information available for the given era. Check for typo in " << era << " or add lumi information in MatrixUnfControl::main." << std::endl;
    std::cout << "Current support for eras: 2017UL, 2018UL, 2016ULpreVFP, 2016ULpostVFP" << std::endl;
    return 1;
  }

  MatrixUnfControl mcontrol;
  gMyRandom.SetSeed(0);

  std::cout << "Setting tau options" << std::endl;

  mcontrol.SetTauOptions(kFALSE, //constrainTau,
    0,     // tausys
    1e-6,  // taumin
    1e-3,  // taumax
    kTRUE);// minRhoAVG

  std::cout << "Setting other options" << std::endl;
  
  mcontrol.SetOptions(era,// era
    lumi,                 // lumi
    kTRUE,                // closure test -- if self consistency is on, this doesn't matter
    selfconsistency,      // self consistency
    ensembletest,         // ensemble test
    linearitytest,        // linearity test
    allowrebin,           // allowrebin
    kFALSE,               // combinedsumptt
    (bkgsuboption == 2 ? kTRUE : kFALSE), // subtbgmcfromdata -- this subtracts bg mc from data before using TUnfold
    (bkgsuboption == 1 ? kTRUE : kFALSE), // dosubtrbg -- turn this on to provide bg distribution to be subtracted by TUnfold
    usebiasintunfold,     // usebiasintunfold
    kTRUE,                // useareaconstrain
    3,                    // regmode  1=size,2=derivative,3=curvature,4=mixed,5=densityCurvature
    setsyst);             // setsyst -- if on it will use TUnfoldSys class's AddSysError()

  const TString ERA = era; // 2017UL, 2018UL, 2016ULpreVFP, and 2016ULpostVFP are possible eras for now
  TString ERA_noUL  = ERA;
  ERA_noUL.ReplaceAll("UL","");

  const TString SYST = systematics; // Nominal or something else from the fVectorOfValidSystematics
  const TString CHAN = channel;     // ee, mumu, emu, or combined?

  // The cross sections are hard-coded at the moment

  // For the case of using all variables from list
  if (variable_name.compare("All") == 0){
    std::cout << "Adding variables for unfolding from TUnfoldVariablesList.txt" << std::endl;
    mcontrol.AddUnfoldVariables("TUnfoldVariablesList.txt");
  }

  // For the case of using just a single variable to parallelize
  else {
    std::cout << "Adding variable for unfolding from command line argument" << std::endl;
    mcontrol.AddUnfoldVariable(variable_name);
  }


  TString preunfoldedfile = "preunfolded_"    + ERA_noUL + "/" + SYST + "/" + CHAN + "/";
  TString tunfoldfile     = "TUnfoldResults_" + ERA_noUL + "/" + SYST + "/" + CHAN + "/";
  TString plotsfile       = "TUnfoldResults_" + ERA_noUL + "/" + SYST + "/" + CHAN + "/Plots";

  std::vector<std::string> VectorOfInputSystematics;
  std::vector<std::string> VectorOfValidSystematics;

  if(setsyst || ( plotoption > 0 && plotoption != 3 ) ){
    // Read the syst file and make a vector use the above vector below
    // Also save a copy of it as fVectorOfValidSystematics
    
    // Read VectorOfValidSystematics from this file
    std::ifstream readsyst("TUnfoldSystList.txt");
    if(!readsyst.is_open()){
      std::cout<<"The file with the list of systematics is not open!"<<std::endl;
      exit(0);
    }
    while(!readsyst.eof()){
      std::string systname;
      getline(readsyst, systname);
      if(systname=="") continue;
      // Save a copy of actual input systematics such as "PU", "KIN" 
      VectorOfInputSystematics.push_back(systname);

      // Make a copy of actual systematics such as "PU_UP", "PU_DOWN", "KIN_UP", "KIN_DOWN"  
      if((systname=="JER")            || (systname.find("JES") != std::string::npos)  || (systname=="PU") || (systname=="LEPT") || (systname=="KIN")   || 
         (systname=="BSEMILEP")       || (systname=="UNCLUSTERED")    || (systname=="UETUNE")     || (systname=="MASS")         || (systname=="MATCH") || 
         (systname=="PDF_ALPHAS")     || (systname=="L1PREFIRING")    || (systname=="ELE_ID")     || (systname=="ELE_RECO")         || (systname=="ELE_SCALESMEARING") || 
         (systname=="MUON_ID")        || (systname=="MUON_ISO")       || (systname=="MUON_SCALE") || (systname=="PSSCALE_WEIGHT_6") || (systname=="PSSCALE_WEIGHT_7")  || 
         (systname=="ELE_SCALE_GAIN") || (systname=="ELE_SCALE_STAT") || (systname=="ELE_SCALE_SYST") || (systname=="TRIG") ||
	       (systname=="BTAG_CORR")      || (systname=="BTAG_UNCORR")    || (systname=="BTAG_LJET_CORR") || (systname=="BTAG_LJET_UNCORR"))
      {
        VectorOfValidSystematics.push_back(systname+"_UP");
        VectorOfValidSystematics.push_back(systname+"_DOWN");
      }
      else if(systname=="SCALE"){
        VectorOfValidSystematics.push_back("MESCALE_UP");
        VectorOfValidSystematics.push_back("MESCALE_DOWN");
        VectorOfValidSystematics.push_back("MEFACSCALE_UP");
        VectorOfValidSystematics.push_back("MEFACSCALE_DOWN");
        VectorOfValidSystematics.push_back("MERENSCALE_UP");
        VectorOfValidSystematics.push_back("MERENSCALE_DOWN");
      }
      else if(systname=="BFRAG"){
        VectorOfValidSystematics.push_back(systname+"_UP");
        VectorOfValidSystematics.push_back(systname+"_DOWN");
        VectorOfValidSystematics.push_back(systname+"_PETERSON");
      }
      else if(systname=="COLORREC"){
        VectorOfValidSystematics.push_back("ERDON");
        VectorOfValidSystematics.push_back("ERDONRETUNE");
        VectorOfValidSystematics.push_back("GLUONMOVETUNE");
      }
      else VectorOfValidSystematics.push_back(systname);       
    }
    // Make a copy of VectorOfValidSystematics
    std::cout << "before add systematics" << std::endl;
    mcontrol.fVectorOfInputSystematics = VectorOfInputSystematics;
    std::cout << "inbetween add systematics" << std::endl;
    mcontrol.fVectorOfValidSystematics = VectorOfValidSystematics;

    // mcontrol.AddSystematics(VectorOfInputSystematics, VectorOfValidSystematics);
    std::cout << "after add systematics" << std::endl;
    }

  if(plotoption<2) {
    std::cout<<"adding samples and data"<<std::endl;
    std::ifstream inputsamples;

    inputsamples.open("TUnfoldHistoFileLists_Lumiweighted_" + ERA_noUL + "/TUnfoldHistoFileList_" + SYST + "_" + CHAN + ".txt");

    if(!inputsamples.is_open()){

      std::cout << "Error opening file TUnfoldHistoFileLists_Lumiweighted_" + ERA_noUL + "/TUnfoldHistoFileList_" + SYST + "_" + CHAN + ".txt\n";
      exit(0);
    }

    // Amandeep : Replace nMCfiles with nchannels 
    // that is the spirit of the variable, atleast in this scope
    int nchannels  = 0;
    if (channel.compare("combined") == 0) nchannels = 3;
    else nchannels = 1;
    // End

    double topxsec = 830.91;

    // Amandeep : Old correction
    // scale factor for the LO DY sample, so we can use both samples to increase statistics but maintain the same total DY yield
    // scale only the 50inf part because that's where the main discrepancy between the two sample yields occurs
    // double tempLODY50infSF = 1.187373682;

    double dySF(1.0); 

    if (era.compare("2018UL") == 0) {
      if (channel=="ee")            dySF = 1.12717;
      else if (channel=="emu")      dySF = 1.10941;
      else if (channel=="mumu")     dySF = 1.09194;
      else if (channel=="combined") dySF = 1.1027;
    } else if (era.compare("2017UL") == 0) {
      if (channel=="ee")            dySF = 1.17378;
      else if (channel=="emu")      dySF = 1.16718;
      else if (channel=="mumu")     dySF = 1.16061;
      else if (channel=="combined") dySF = 1.16462;
    } else if (era.compare("2016UL") == 0) {
      if (channel=="ee")            dySF = 1.18945;
      else if (channel=="emu")      dySF = 1.16376;
      else if (channel=="mumu")     dySF = 1.13863;
      else if (channel=="combined") dySF = 1.15266;
    } else if (era.compare("2016ULpostVFP") == 0) {
      if (channel=="ee")            dySF = 1.19939;
      else if (channel=="emu")      dySF = 1.1901;
      else if (channel=="mumu")     dySF = 1.18087;
      else if (channel=="combined") dySF = 1.18599;
    } else if (era.compare("2016ULpreVFP") == 0) {
      if (channel=="ee")            dySF = 1.1809;
      else if (channel=="emu")      dySF = 1.14192;
      else if (channel=="mumu")     dySF = 1.10422;
      else if (channel=="combined") dySF = 1.12537;
    } else if (era.compare("fullRun2UL") == 0) {
      if (channel=="ee")            dySF = 1.15707;
      else if (channel=="emu")      dySF = 1.14002;
      else if (channel=="mumu")     dySF = 1.12321;
      else if (channel=="combined") dySF = 1.13327;
    }

    while (!inputsamples.eof()) {
      std::string samplename;
      getline(inputsamples,samplename);

      if (samplename == "") continue; // skip empty lines
      if(samplename.find("run") != std::string::npos) mcontrol.AddDataFile(samplename.c_str(), 1);

      // Amandeep : Update new cross-sections
      // Note that samplename is a std::string and not TString object
      // hence we cannot use the Contains(<desired_string>) method, 
      // but have to use find(<desired_string>) != std::string::npos instead
      // Cross-sections borrowed from : https://gitlab.cern.ch/jthieman/TopAnalysis/-/blob/SpinCorr-Run2-UL/Configuration/analysis/diLeptonic/src/PlotterConfigurationHelper.cc#L1209
      
      else if (samplename.find("ee_ttbarsignalplustau_fromDilepton")  != std::string::npos) {mcontrol.AddMCFile(samplename.c_str(), topxsec * 0.10706 * 0.964261576); } // ++nMCfiles; } // updated PDG 2020 with Dilepton BR Corrections
      else if (samplename.find("emu_ttbarsignalplustau_fromDilepton") != std::string::npos) {mcontrol.AddMCFile(samplename.c_str(), topxsec * 0.10706 * 0.957058875); } // ++nMCfiles; } // updated PDG 2020 with Dilepton BR Corrections
      else if (samplename.find("mumu_ttbarsignalplustau_fromDilepton")!= std::string::npos) {mcontrol.AddMCFile(samplename.c_str(), topxsec * 0.10706 * 0.949909976); } // ++nMCfiles; } // updated PDG 2020 with Dilepton BR Corrections
      
      // Amandeep : ALSO SIGNAL MC
      // I don't increment the nmcfiles counter here, but these are also part of signal MC 
      // This is because nmcfiles == 1 corresponds to single channel L9337
      // However we need to use the signalviatau samples to fully consider the dilepton channel

      else if (samplename.find("ee_ttbarsignalviatau_fromDilepton")   != std::string::npos) {mcontrol.AddMCFile(samplename.c_str(), topxsec * 0.10706 * 1.029827957); }// updated PDG 2020 with Dilepton BR Corrections
      else if (samplename.find("emu_ttbarsignalviatau_fromDilepton")  != std::string::npos) {mcontrol.AddMCFile(samplename.c_str(), topxsec * 0.10706 * 1.026209047); }// updated PDG 2020 with Dilepton BR Corrections
      else if (samplename.find("mumu_ttbarsignalviatau_fromDilepton") != std::string::npos) {mcontrol.AddMCFile(samplename.c_str(), topxsec * 0.10706 * 1.022670477); }// updated PDG 2020 with Dilepton BR Corrections
      
      else if (samplename.find("bg_fromDilepton")!= std::string::npos) {mcontrol.AddBgFile(samplename.c_str(), topxsec * 0.10706);}// updated PDG 2020
      else if (samplename.find("fromLjets")!= std::string::npos)       {mcontrol.AddBgFile(samplename.c_str(), topxsec * 0.44113);}// updated PDG 2020
      else if (samplename.find("fromHadronic")!= std::string::npos)    {mcontrol.AddBgFile(samplename.c_str(), topxsec * 0.45441);}// updated PDG 2020

      else if ( samplename.find("ttbar") != std::string::npos &&
              !(samplename.find("ttbarW")!= std::string::npos) &&
              !(samplename.find("ttbarZ")!= std::string::npos) )      {mcontrol.AddBgFile(samplename.c_str(),topxsec);}// updated
      
      else if (samplename.find("single") != std::string::npos &&
               // NoFullyHadronic
               samplename.find("tw") != std::string::npos)            {mcontrol.AddBgFile(samplename.c_str(),35.85*0.54559);}// updated
              
      else if (samplename.find("single") != std::string::npos &&
               samplename.find("_s") != std::string::npos)                 {mcontrol.AddBgFile(samplename.c_str(), 10.32);}// updated
      else if (samplename.find("singletop") != std::string::npos &&
               samplename.find("_t") != std::string::npos)                 {mcontrol.AddBgFile(samplename.c_str(), 136.02);}// updated
      else if (samplename.find("singleantitop") != std::string::npos &&
               samplename.find("_t") != std::string::npos)                 {mcontrol.AddBgFile(samplename.c_str(), 80.95);}// updated
      
      else if (samplename.find("ww") != std::string::npos)                 {mcontrol.AddBgFile(samplename.c_str(), 118.7);}// updated
      else if (samplename.find("wz") != std::string::npos)                 {mcontrol.AddBgFile(samplename.c_str(), 47.13);}// updated 04.09.17: NLO from MCFM
      else if (samplename.find("zz") != std::string::npos)                 {mcontrol.AddBgFile(samplename.c_str(), 16.523);}// updated 04.09.17: NLO from MCFM
      else if (samplename.find("1050") != std::string::npos)               {mcontrol.AddBgFile(samplename.c_str(), 22635.1 * dySF);}// from NNLO FEWZ 3.1

      // Splitting the 50inf half between the amcatnlofxfx and madgraphmlm
      else if (samplename.find("50inf_amcatnlofxfx") != std::string::npos) {mcontrol.AddBgFile(samplename.c_str(), 0.5 * 3.*2075.14 * dySF);}// updated 23.05.19: NNLO FEWZ 3.1.b2 NNPDF31_nnlo_as_0118_luxqed
      else if (samplename.find("50inf_madgraphmlm")  != std::string::npos) {mcontrol.AddBgFile(samplename.c_str(), 0.5 * 3.*2075.14 * dySF);}// updated 23.05.19: NNLO FEWZ 3.1.b2 NNPDF31_nnlo_as_0118_luxqed

      else if (samplename.find("wtolnu")  != std::string::npos)            {mcontrol.AddBgFile(samplename.c_str(), 61526.7);}// updated 22.05.19

      else if (samplename.find("ttbarWjetstolnu") != std::string::npos)    {mcontrol.AddBgFile(samplename.c_str(), 0.2043);}// updated: nlo
      else if (samplename.find("ttbarWjetstoqq")  != std::string::npos)    {mcontrol.AddBgFile(samplename.c_str(), 0.4062);}// updated: nlo
      else if (samplename.find("ttbarZtollnunu")  != std::string::npos)    {mcontrol.AddBgFile(samplename.c_str(), 0.2529);}// updated: nlo
      else if (samplename.find("ttbarZtoqq") != std::string::npos)         {mcontrol.AddBgFile(samplename.c_str(), 0.5297);}// updated: nlo
      else mcontrol.AddBgFile(samplename.c_str(), topxsec);
      // End
    }


    std::vector<TString>channels;
    channels.push_back("ee");
    channels.push_back("emu");
    channels.push_back("mumu");
    
    // ***********
    // Systematics
    // ***********
    
    if(setsyst){
      std::cout<<"adding systematic samples"<<std::endl;
      // Over systematics
      for(unsigned int isyst = 0; isyst < VectorOfValidSystematics.size(); isyst++){ 
        std::cout<<"Adding "<<VectorOfValidSystematics[isyst]<<" systematics"<<std::endl;
        std::string systname = VectorOfValidSystematics[isyst]; 

        // Over channels
        for (int i_chan = 0; i_chan < nchannels; ++i_chan)
        {

          TString CHANstring;
          // Amandeep : replace nMCfiles with nchannels
          if(nchannels == 1) CHANstring = CHAN;
          else CHANstring = channels.at(i_chan);


          // fullRun2 combined
          if (ERA_noUL == "fullRun2") {
            if (CHANstring == "ee"){
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_fullRun2/"  + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalplustau_fromDilepton_fullRun2UL.root" , topxsec * 0.10706 * 0.964261576);
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_fullRun2/"  + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalviatau_fromDilepton_fullRun2UL.root"  , topxsec * 0.10706 * 1.029827957);
            }

            else if (CHANstring == "emu"){
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_fullRun2/"  + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalplustau_fromDilepton_fullRun2UL.root" , topxsec * 0.10706 * 0.957058875);
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_fullRun2/"  + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalviatau_fromDilepton_fullRun2UL.root"  , topxsec * 0.10706 * 1.026209047);
            }
            
            else if (CHANstring == "mumu"){
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_fullRun2/"  + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalplustau_fromDilepton_fullRun2UL.root" , topxsec * 0.10706 * 0.949909976);
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_fullRun2/"  + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalviatau_fromDilepton_fullRun2UL.root"  , topxsec * 0.10706 * 1.022670477);
            }
          }

          // 2016 combined
          else if (ERA_noUL == "2016") {
            if (CHANstring == "ee"){
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_2016preVFP/"  + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalplustau_fromDilepton_2016ULpreVFP.root" , topxsec * 0.10706 * 0.964261576);
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_2016preVFP/"  + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalviatau_fromDilepton_2016ULpreVFP.root"  , topxsec * 0.10706 * 1.029827957);
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_2016postVFP/" + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalplustau_fromDilepton_2016ULpostVFP.root", topxsec * 0.10706 * 0.964261576);
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_2016postVFP/" + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalviatau_fromDilepton_2016ULpostVFP.root" , topxsec * 0.10706 * 1.029827957);          
            }
            else if (CHANstring == "emu"){
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_2016preVFP/"  + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalplustau_fromDilepton_2016ULpreVFP.root" , topxsec * 0.10706 * 0.957058875);
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_2016preVFP/"  + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalviatau_fromDilepton_2016ULpreVFP.root"  , topxsec * 0.10706 * 1.026209047);
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_2016postVFP/" + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalplustau_fromDilepton_2016ULpostVFP.root", topxsec * 0.10706 * 0.957058875);
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_2016postVFP/" + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalviatau_fromDilepton_2016ULpostVFP.root" , topxsec * 0.10706 * 1.026209047);
            }
            else if (CHANstring == "mumu"){
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_2016preVFP/"  + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalplustau_fromDilepton_2016ULpreVFP.root" , topxsec * 0.10706 * 0.949909976);
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_2016preVFP/"  + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalviatau_fromDilepton_2016ULpreVFP.root"  , topxsec * 0.10706 * 1.022670477);
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_2016postVFP/" + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalplustau_fromDilepton_2016ULpostVFP.root", topxsec * 0.10706 * 0.949909976);
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_2016postVFP/" + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalviatau_fromDilepton_2016ULpostVFP.root" , topxsec * 0.10706 * 1.022670477);            
            }
          }

          // Default case
          else  {
            if (CHANstring == "ee"){
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_" + ERA_noUL + "/" + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalplustau_fromDilepton_" + ERA + ".root", topxsec * 0.10706 * 0.964261576);
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_" + ERA_noUL + "/" + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalviatau_fromDilepton_" + ERA + ".root", topxsec * 0.10706 * 1.029827957);
            }
            else if (CHANstring == "emu"){
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_" + ERA_noUL + "/" + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalplustau_fromDilepton_" + ERA + ".root", topxsec * 0.10706 * 0.957058875);
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_" + ERA_noUL + "/" + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalviatau_fromDilepton_" + ERA + ".root", topxsec * 0.10706 * 1.026209047);
            }
            else if (CHANstring == "mumu"){
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_" + ERA_noUL + "/" + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalplustau_fromDilepton_" + ERA + ".root", topxsec * 0.10706 * 0.949909976);
              mcontrol.AddSysFile("UnfoldingHistos_Lumiweighted_" + ERA_noUL + "/" + VectorOfValidSystematics[isyst] + "/" + CHANstring + "/" + "histosTUnfold_" + CHANstring + "_ttbarsignalviatau_fromDilepton_" + ERA + ".root", topxsec * 0.10706 * 1.022670477);
            }
          }
          // End
        }
      }
    }

    // All histo arrays for sig, bg, data obtained here
    std::cout<<"filling histogram arrays"<<std::endl;
    mcontrol.FillHistoArrays(CHAN);

    gSystem->mkdir(preunfoldedfile.Data(), true);
    gSystem->mkdir(tunfoldfile.Data(), true);
    gSystem->mkdir(plotsfile.Data(), true);
    mcontrol.DoAllUnfold(preunfoldedfile, tunfoldfile, CHAN);
  }
  if(plotoption>0) {
    bool doOptimiseTotalUnc = setsyst > 1 ? kTRUE : kFALSE;
    gSystem->mkdir(plotsfile.Data(), true);
    if(plotoption!=3) mcontrol.MakeAllPlots(tunfoldfile, plotsfile, CHAN, doOptimiseTotalUnc);
    if( (ensembletest && plotoption!=2) || plotoption==3) mcontrol.doPlotPulls(tunfoldfile, plotsfile);
  }
  return 1;
}


// ********************
// CalcUpDownDifference
// ********************

TH1D* MatrixUnfControl::CalcUpDownDifference ( TString Channel, TString SystName, TString Variable, TFile* unffile, Bool_t normXsec, TString NameTag ){
  std::cout << "Starting MatrixUnfControl::CalcUpDownDifference" << std::endl;
  TH1D *hnom, *deltarelavg, *hnom_norm;  
  hnom = (TH1D*) ((TH1D*)unffile->Get(Variable+"TUnfResultCor"+NameTag))->Clone("hnom");

  if(normXsec){
    hnom_norm = (TH1D*) hnom->Clone("hnom_norm");
    hnom_norm->Scale(1./hnom_norm->Integral());
  }
  std::cout<<"SystName: "<<SystName<<std::endl;

  // Amandeep : Samples with UP and DOWN variations
  // Jason mentioned not enveloping BTAG and BTAG_LJET
  if( (SystName=="JER") || (SystName.Contains("JES")) || (SystName=="PU")         || (SystName=="LEPT") || 
      (SystName=="KIN") || (SystName=="BSEMILEP")     || (SystName=="UNCLUSTERED")|| (SystName=="TRIG") ||
      (SystName=="UETUNE") || (SystName=="MASS")      || (SystName=="MATCH")      || (SystName=="PDF_ALPHAS") || 
      (SystName=="L1PREFIRING")    || (SystName=="ELE_ID")         || (SystName=="ELE_RECO")       || 
      (SystName=="ELE_SCALE_GAIN") || (SystName=="ELE_SCALE_SYST") || (SystName=="ELE_SCALE_STAT") || 
      (SystName=="BTAG_CORR")      || (SystName=="BTAG_UNCORR")    || (SystName=="BTAG_LJET_CORR") || (SystName=="BTAG_LJET_UNCORR") || 
      (SystName=="MUON_ID")        || (SystName=="MUON_ISO")       || (SystName=="MUON_SCALE")     || (SystName=="PSSCALE_WEIGHT_6") || 
      (SystName=="PSSCALE_WEIGHT_7"))
    {
    TH1D* deltaup   = (TH1D*) ((TH1D*)unffile->Get(Variable + "RespMatSys_" + SystName + "_UP_deltaSys"   + NameTag))->Clone("deltaup");
    TH1D* deltadown = (TH1D*) ((TH1D*)unffile->Get(Variable + "RespMatSys_" + SystName + "_DOWN_deltaSys" + NameTag))->Clone("deltadown");
    deltarelavg     = (TH1D*) hnom->Clone(Variable+"RespMatSys_" + SystName + "_deltaSysRel" + NameTag + (normXsec?"_Norm":""));

    if(normXsec) {
      deltaup->Add(hnom);
      deltadown->Add(hnom);
      deltaup->Scale(1./deltaup->Integral());
      deltadown->Scale(1./deltadown->Integral());
      deltaup->Add(hnom_norm,-1);
      deltadown->Add(hnom_norm,-1);
    }
 
    for(int ibin = 1; ibin<=deltarelavg->GetNbinsX(); ibin++){
      double bincontent_temp = hnom->GetBinContent(ibin);
      if(normXsec) bincontent_temp = hnom_norm->GetBinContent(ibin);
      deltarelavg->SetBinContent(ibin, (fabs(deltaup->GetBinContent(ibin) - deltadown->GetBinContent(ibin))*0.5)/bincontent_temp);
    }

  delete deltaup;
  delete deltadown;
  }

  // *******************
  // Alternative samples 
  // *******************

  else if((SystName=="MCATNLO") || (SystName=="POWHEG") || (SystName=="POWHEGHERWIG") || (SystName=="POWHEGV2") || (SystName=="POWHEGV2HERWIG") || (SystName=="AMCATNLOFXFX") || (SystName=="TOP_PT") || (SystName=="MADGRAPHMLM") || (SystName=="PERUGIA11") || (SystName=="SPINCORR") || (SystName.Contains("PDF")) ){
    deltarelavg = (TH1D*)((TH1D*)unffile->Get(Variable+"RespMatSys_"+SystName+"_deltaSys"+NameTag))->Clone(Variable+"RespMatSys_"+SystName+"_deltaSysRel"+NameTag+(normXsec?"_Norm":""));

    if(normXsec) {
      deltarelavg->Add(hnom);
      deltarelavg->Scale(1./deltarelavg->Integral());
      deltarelavg->Add(hnom_norm,-1);
    }

    for(int ibin = 1; ibin<=deltarelavg->GetNbinsX(); ibin++){
      double bincontent_temp = hnom->GetBinContent(ibin);
      if(normXsec) bincontent_temp = hnom_norm->GetBinContent(ibin);
      deltarelavg->SetBinContent(ibin, fabs(deltarelavg->GetBinContent(ibin)/bincontent_temp));
    }
  }

  // ***********
  // Backgrounds
  // ***********

  else if( (SystName=="DYeemm") || (SystName=="DYtautau") || (SystName=="singletop")  || (SystName=="ww") || (SystName=="wz") || (SystName=="zz") || (SystName=="ttbarW") || (SystName=="ttbarZ") || (SystName=="ttbarbg") || (SystName=="wtolnu") || (SystName=="other") ) {
    deltarelavg = (TH1D*)((TH1D*)unffile->Get(Variable + "BkgSys_" + SystName + "_deltaSys" + NameTag))->Clone(Variable + "BkgSys_" + SystName + "_deltaSysRel" + NameTag + (normXsec?"_Norm":""));

    if(normXsec) {
      deltarelavg->Add(hnom);
      deltarelavg->Scale(1./deltarelavg->Integral());
      deltarelavg->Add(hnom_norm,-1);
    }

    for(int ibin = 1; ibin<=deltarelavg->GetNbinsX(); ibin++){
      double bincontent_temp = hnom->GetBinContent(ibin);
      if(normXsec) bincontent_temp = hnom_norm->GetBinContent(ibin);
      deltarelavg->SetBinContent(ibin, fabs(deltarelavg->GetBinContent(ibin)/bincontent_temp));
    }
  }

  // ******************
  // Flat uncertainties
  // ******************

  else if( (SystName=="Lumi") || (SystName=="BR") ) {
    double flat_uncertainty = 0.025;
    if(SystName=="BR") flat_uncertainty = 0.015;

    deltarelavg = (TH1D*) hnom->Clone(Variable + "FlatSys_" + SystName + "_deltaSysRel" + NameTag + (normXsec?"_Norm":""));

    if(normXsec) {
      deltarelavg->Reset();
    }
    else {

      for(int ibin = 1; ibin<=deltarelavg->GetNbinsX(); ibin++){
        deltarelavg->SetBinContent(ibin, flat_uncertainty);
      }
    }

  }

  else 
  {
    std::cout<<"SystName is not one of the options defined in the CalcUpDownDifference"<<std::endl;
    return 0;
  }

  delete hnom;
  if(normXsec) delete hnom_norm; 
  return deltarelavg;
}


// *****************************
// GetEnvelopeForUnfoldedResults
// *****************************


TH1D* MatrixUnfControl::GetEnvelopeForUnfoldedResults( TString Channel, TString SystName, TString Variable, TFile* unffile, Bool_t normXsec, TString NameTag ){
  
  if(SystName!="SCALE" && SystName!="BFRAG" && SystName!="COLORREC") {
    std::cout<<"The calculation of envelope is not supported for "<<SystName<<". Exiting the function!!"<<std::endl;
    return 0;
  }
  
  std::vector<TString> syst;
  
  // Amandeep : removing the PS weights from here for ISR and FSR since they aren't properly defined 
  if(SystName=="SCALE")     syst = { "MESCALE_UP", "MESCALE_DOWN", "MEFACSCALE_UP", "MEFACSCALE_DOWN", "MERENSCALE_UP", "MERENSCALE_DOWN" };  
  if(SystName=="BFRAG")     syst = { "BFRAG_UP", "BFRAG_DOWN", "BFRAG_PETERSON" };
  if(SystName=="COLORREC")  syst = { "ERDON"   , "ERDONRETUNE", "GLUONMOVETUNE" };

  TH1D *htot_reldelta = (TH1D*)((TH1D*)unffile->Get(Variable + "TUnfResultCor" + NameTag)->Clone(Variable + "RespMatSys_" + SystName + "_deltaSysRel" + NameTag + (normXsec?"_Norm":"")));
  TH1D *hnom          = (TH1D*)((TH1D*)unffile->Get(Variable + "TUnfResultCor" + NameTag))->Clone("hnom");

  TH1D *hnom_norm;
  if(normXsec){
    hnom_norm = (TH1D*) hnom->Clone("hnom_norm");
    hnom_norm->Scale(1./hnom_norm->Integral());
  }

  std::vector<double> max_up;
  std::vector<double> max_down;
  
  for (size_t iter = 0; iter<syst.size(); iter++){

    double componentscalingfactor = 1.;
    //    if(SystName=="SCALE" && syst[iter].BeginsWith("PSFSR") ) componentscalingfactor = sqrt(2.); // applied only to FSR uncertainties

    TH1D *delta = (TH1D*) ((TH1D*)unffile->Get(Variable + "RespMatSys_" + syst[iter] + "_deltaSys" + NameTag))->Clone("delta");

    if(normXsec) {
      delta->Add(hnom);
      delta->Scale(1./delta->Integral());
      delta->Add(hnom_norm,-1);
    }

    delta->Scale(1./componentscalingfactor);

    for( int ibins = 1; ibins<=delta->GetNbinsX(); ibins++){
      if(iter==0) {
        max_up.push_back(0);
        max_down.push_back(0);
      }
      if(max_up[ibins-1]   < delta->GetBinContent(ibins)) max_up[ibins-1]   = delta->GetBinContent(ibins); 
      if(max_down[ibins-1] > delta->GetBinContent(ibins)) max_down[ibins-1] = delta->GetBinContent(ibins);
    }
    delete delta;
  }
  
  for(int ibins = 1; ibins<=htot_reldelta->GetNbinsX(); ibins++){
    double bincontent_temp = hnom->GetBinContent(ibins);
    if(normXsec) bincontent_temp = hnom_norm->GetBinContent(ibins);
    htot_reldelta->SetBinContent(ibins, (fabs(max_up[ibins-1]-max_down[ibins-1])*0.5)/bincontent_temp);
  }
  
  delete hnom;
  if(normXsec) delete hnom_norm; 
  std::cout<<"Finished MatrixUnfControl::GetEnvelopeForUnfoldedResults."<<std::endl;
  return htot_reldelta;
}


// ***************************
// GetCombinedCovarianceMatrix
// ***************************

TH2D* MatrixUnfControl::GetCombinedCovarianceMatrix( TString Channel, TString SystName, TString Variable, TFile* unffile, TH1D* htot_reldelta, Bool_t normXsec, TString NameTag ){

  // Calculate the correlation matrix for the sum of sub-components, then convert this to a covariance matrix using the envelope results
  std::cout<<"Starting MatrixUnfControl::GetCombinedCovarianceMatrix"<<std::endl;

  TH1D* htot_absdelta = (TH1D*) htot_reldelta->Clone("htot_absdelta");
  TH1D *hnom          = (TH1D*) ((TH1D*)unffile->Get(Variable + "TUnfResultCor" + NameTag))->Clone("hnom");

  TH1D *hnom_norm;
  if(normXsec){
    hnom_norm = (TH1D*) hnom->Clone("hnom_norm");
    hnom_norm->Scale(1./hnom_norm->Integral());
  }

  for(int ibins = 1; ibins<=htot_reldelta->GetNbinsX(); ibins++){
    if(normXsec) htot_absdelta->SetBinContent(ibins, htot_reldelta->GetBinContent(ibins)*hnom_norm->GetBinContent(ibins) );
    else         htot_absdelta->SetBinContent(ibins, htot_reldelta->GetBinContent(ibins)*hnom->GetBinContent(ibins) );
  }

  std::vector<TString> syst;
  bool isBkgSyst = kFALSE;
  
  if(SystName=="SCALE")    syst = { "MESCALE_UP", "MESCALE_DOWN", "MEFACSCALE_UP", "MEFACSCALE_DOWN", "MERENSCALE_UP", "MERENSCALE_DOWN" };
  if(SystName=="BFRAG")    syst = { "BFRAG_UP", "BFRAG_DOWN", "BFRAG_PETERSON" };
  if(SystName=="COLORREC") syst = { "ERDON", "ERDONRETUNE", "GLUONMOVETUNE" };

  // Amandeep : Samples with UP and DOWN variations 
  // Jason mentioned not enveloping BTAG and BTAG_LJET

  if( (SystName=="JER")  || (SystName.Contains("JES")) || (SystName=="PU")         || (SystName=="LEPT") || 
      (SystName=="KIN")  || (SystName=="BSEMILEP")     || (SystName=="UNCLUSTERED")|| (SystName=="TRIG") ||
      (SystName=="UETUNE") || (SystName=="MASS")       || (SystName=="MATCH")      || (SystName=="PDF_ALPHAS") || 
      (SystName=="L1PREFIRING")    || (SystName=="ELE_ID")         || (SystName=="ELE_RECO") || 
      (SystName=="ELE_SCALE_GAIN") || (SystName=="ELE_SCALE_SYST") || (SystName=="ELE_SCALE_STAT") || 
      (SystName=="BTAG_CORR")      || (SystName=="BTAG_UNCORR")    || (SystName=="BTAG_LJET_CORR") || (SystName=="BTAG_LJET_UNCORR") || 
      (SystName=="MUON_ID")        || (SystName=="MUON_ISO")       || (SystName=="MUON_SCALE")     || (SystName=="PSSCALE_WEIGHT_6") || 
      (SystName=="PSSCALE_WEIGHT_7")) 
      syst = {SystName+"_UP", SystName+"_DOWN"};
  
  else if(SystName.Contains("PDF")) syst = {SystName};

  // *******************
  // Alternative samples 
  // *******************

  if((SystName=="MCATNLO") || (SystName=="POWHEG") || (SystName=="POWHEGHERWIG") || (SystName=="POWHEGV2") || (SystName=="POWHEGV2HERWIG") || (SystName=="AMCATNLOFXFX") || (SystName=="TOP_PT") || (SystName=="MADGRAPHMLM") || (SystName=="PERUGIA11") || (SystName=="SPINCORR")) syst = {SystName};
  
  // ***********
  // Backgrounds
  // ***********

  if( (SystName=="DYeemm") || (SystName=="DYtautau") || (SystName=="singletop")  || (SystName=="ww") || (SystName=="wz") || (SystName=="zz") || (SystName=="ttbarW") || (SystName=="ttbarZ") || (SystName=="ttbarbg") || (SystName=="wtolnu") || (SystName=="other") ) {
    syst = {SystName};
    isBkgSyst = kTRUE;
  }

  double scalingfactor = 1.;
  if(SystName=="MASS")                scalingfactor=1.;  // convert from +/-3 GeV to +/-1 GeV
  if(SystName.Contains("PDF"))        scalingfactor=10.; // convert from sum in quadrature to RMS (for 100 replicas)
  if(SystName.Contains("PDF_ALPHAS")) scalingfactor=1.;

  const int xbins=hnom->GetNbinsX();
  const int ybins=hnom->GetNbinsX();

  TMatrixD systCov(xbins,ybins);
  // std::vector<std::vector<double>>systCov(xbins, std::vector<double>(ybins));

  for(unsigned int isyst = 0;isyst<syst.size();isyst++){

    double componentscalingfactor = 1.;
    // if(SystName=="SCALE" && syst[isyst].BeginsWith("PSFSR") ) componentscalingfactor=sqrt(2.); // applied only to FSR uncertainties

    TH1D *delta = (TH1D*) ((TH1D*)unffile->Get(Variable + (isBkgSyst?"BkgSys_":"RespMatSys_") + syst[isyst] + "_deltaSys" + NameTag))->Clone("delta");

    if(normXsec) {
      delta->Add(hnom);
      delta->Scale(1./delta->Integral());
      delta->Add(hnom_norm,-1);
    }

    for(int xbin=1; xbin<=xbins; xbin++){
      for(int ybin=1; ybin<=ybins; ybin++){
       systCov(xbin-1,ybin-1) += (delta->GetBinContent(xbin)*delta->GetBinContent(ybin)/componentscalingfactor/componentscalingfactor);
      }
    }

    delete delta;
  }

  // ******************
  // Flat uncertainties
  // ******************

  // Convert to correlation matrix (or fill fully-correlated correlation matrix for the flat systematics)
  if( (SystName=="Lumi") || (SystName=="BR") ) {
    for(int xbin=1; xbin<=xbins; xbin++){
      for(int ybin=1; ybin<=ybins; ybin++){
        systCov(xbin-1,ybin-1) = 1;
      }
    }
  }

  else {
    TVectorD diag_nozeros = TMatrixDDiag(systCov);
    for(int xbin=1; xbin<=xbins; xbin++){
      if( diag_nozeros(xbin-1) == 0 ) diag_nozeros(xbin-1) = 1e-20;
    }
    systCov.NormByDiag(diag_nozeros);
  }

  // Convert back to covariance matrix, but using the bin errors from the envelope method
  for(int xbin=1; xbin<=xbins; xbin++){
    for(int ybin=1; ybin<=ybins; ybin++){
      systCov(xbin-1,ybin-1)*= fabs( htot_absdelta->GetBinContent(xbin)*htot_absdelta->GetBinContent(ybin) );
    }
  }

  TH2D *systCovMatrix = new TH2D(Variable + "_systCovMatrix" + NameTag + (normXsec?"_Norm":""), Variable + "_systCovMatrix" + NameTag + (normXsec?"_Norm":""), xbins,0,xbins, xbins,0,xbins);
  // if(!normXsec) systCovMatrix = (TH2D*)((TH2D*)unffile->Get(Variable+(isBkgSyst?"BkgSys_":"RespMatSys_") + syst[0] + "_CovMatSys" + NameTag))->Clone(Variable + "_systCovMatrix" + NameTag + (normXsec?"_Norm":"")); 
  // else          systCovMatrix = (TH2D*)((TH2D*)unffile->Get(Variable+(isBkgSyst?"BkgSys_":"RespMatSys_") + syst[0] + "_CovMatSysNorm" + NameTag))->Clone(Variable + "_systCovMatrix" + NameTag + (normXsec?"_Norm":"")); 
  
  for(int xbin=1; xbin<=systCovMatrix->GetNbinsX(); xbin++){
    for(int ybin=1; ybin<=systCovMatrix->GetNbinsY(); ybin++){
      systCovMatrix->SetBinContent(xbin,ybin,systCov(xbin-1,ybin-1)/scalingfactor/scalingfactor);
    }
  }

  // Envelope method does not preserve covariance characteristics of a normalised differential xsec 
  // (e.g, the error on the sum of all bins should equal zero). 
  // Fix this by recalculating each bin as bin/(sum of all bins):
  if(normXsec) MatrixUnf::GetNormalisedCovarianceMatrix(hnom_norm, systCovMatrix, systCovMatrix);

  delete hnom;
  if(normXsec) delete hnom_norm;
  delete htot_absdelta;
  return systCovMatrix;
}


// ****************
// GetDeltasAllVars
// ****************

void MatrixUnfControl::GetDeltasAllVars(TH1D*& hnom_AllVar, TH1D*& htot_reldelta_AllVar, TString Channel, TString SystName, TString TUNFFILE, Bool_t normXsec, TString NameTag ){
  
  // Amandeep : Hardcoded here
  const int nbinsrhoi    = 24;
  // const int nbinsrhoi = 6;
  
  // This constructs a Nominal and Delta vector which is as wide as n_vars * n_bins
  hnom_AllVar          = new TH1D(Form("hnom_AllVar%s%s",(normXsec?"Norm":""),NameTag.Data()), Form("hnom_AllVar%s%s",(normXsec?"Norm":""),NameTag.Data()), fVariableNames.size()*nbinsrhoi,0,fVariableNames.size()*nbinsrhoi);
  htot_reldelta_AllVar = new TH1D(Form("htot_reldelta_AllVar%s%s",(normXsec?"Norm":""),NameTag.Data()), Form("htot_reldelta_AllVar%s%s",(normXsec?"Norm":""),NameTag.Data()), fVariableNames.size()*nbinsrhoi,0,fVariableNames.size()*nbinsrhoi);

  int bin_counter = 1;

  // Over the variables
  for(unsigned int i=0; i < fVariableNames.size(); i++)
  {

    TString Variable = fVariableNames.at(i);

    if(Variable=="") continue;

    TString TUnfoldResultsFilename = Form(TUNFFILE + Variable + ".root");
    TFile*  unffile                = new TFile(TUnfoldResultsFilename.Data(),"READ");
    if(unffile==NULL) exit(0);

    TH1D *htot_reldelta = (TH1D*)((TH1D*)unffile->Get(Variable + "TUnfResultCor" + NameTag)->Clone(Variable + "RespMatSys_" + SystName + "_deltaSysRel" + NameTag + (normXsec?"_Norm":"")));
    TH1D *hnom          = (TH1D*)((TH1D*)unffile->Get(Variable + "TUnfResultCor" + NameTag))->Clone("hnom");

    TH1D *hnom_norm;
    if(normXsec){
      hnom_norm = (TH1D*) hnom->Clone("hnom_norm");
      hnom_norm->Scale(1./hnom_norm->Integral());
    }


    if(SystName=="SCALE" || SystName=="BFRAG" || SystName=="COLORREC") {

      std::vector<TString> syst;

      if(SystName=="SCALE")    syst = { "MESCALE_UP", "MESCALE_DOWN", "MEFACSCALE_UP", "MEFACSCALE_DOWN", "MERENSCALE_UP", "MERENSCALE_DOWN" };
      if(SystName=="BFRAG")    syst = { "BFRAG_UP", "BFRAG_DOWN", "BFRAG_PETERSON" };
      if(SystName=="COLORREC") syst = { "ERDON", "ERDONRETUNE", "GLUONMOVETUNE" };

      std::vector<double> max_up;
      std::vector<double> max_down;

      for (size_t iter = 0; iter < syst.size(); iter++){

        double componentscalingfactor = 1.;
	      // if(SystName=="SCALE" && syst[iter].BeginsWith("PSFSR") ) componentscalingfactor = sqrt(2.); // applied only to FSR uncertainties

        // This is the fluctuation with dimension n_bins
        TH1D *delta = (TH1D*) ((TH1D*)unffile->Get(Variable+"RespMatSys_"+syst[iter]+"_deltaSys"+NameTag))->Clone("delta");

        if(normXsec) {
          delta->Add(hnom);
          delta->Scale(1./delta->Integral());
          delta->Add(hnom_norm,-1);
        }

        delta->Scale(1./componentscalingfactor);

        for( int ibins = 1; ibins<=delta->GetNbinsX(); ibins++){
          if(iter==0) {
            max_up.push_back(0);
            max_down.push_back(0);
          }
          if(max_up[ibins-1]   < delta->GetBinContent(ibins)) max_up[ibins-1]   = delta->GetBinContent(ibins); 
          if(max_down[ibins-1] > delta->GetBinContent(ibins)) max_down[ibins-1] = delta->GetBinContent(ibins);
        }
        delete delta;
      }

      for(int ibins = 1; ibins<=htot_reldelta->GetNbinsX(); ibins++){
        double bincontent_temp       = hnom->GetBinContent(ibins);
        if(normXsec) bincontent_temp = hnom_norm->GetBinContent(ibins);
        htot_reldelta->SetBinContent(ibins, (fabs(max_up[ibins-1]-max_down[ibins-1])*0.5)/bincontent_temp);
      }

    }

    // Amandeep : Samples with UP and DOWN variations 
    // Jason mentioned not enveloping BTAG and BTAG_LJET
    else if( 
      (SystName=="JER")  || (SystName.Contains("JES")) || (SystName=="PU")         || (SystName=="LEPT") || 
      (SystName=="KIN")  || (SystName=="BSEMILEP")     || (SystName=="UNCLUSTERED")|| (SystName=="TRIG") ||
      (SystName=="UETUNE") || (SystName=="MASS")       || (SystName=="MATCH")      || (SystName=="PDF_ALPHAS") || 
      (SystName=="L1PREFIRING")    || (SystName=="ELE_ID")         || (SystName=="ELE_RECO") || 
      (SystName=="ELE_SCALE_GAIN") || (SystName=="ELE_SCALE_SYST") || (SystName=="ELE_SCALE_STAT") || 
      (SystName=="BTAG_CORR")      || (SystName=="BTAG_UNCORR")    || (SystName=="BTAG_LJET_CORR") || (SystName=="BTAG_LJET_UNCORR") || 
      (SystName=="MUON_ID")        || (SystName=="MUON_ISO")       || (SystName=="MUON_SCALE")     || (SystName=="PSSCALE_WEIGHT_6") || 
      (SystName=="PSSCALE_WEIGHT_7"))
      {
      
      // This is the upward or downward fluctuation with dimension n_bins
      TH1D* deltaup   = (TH1D*) ((TH1D*)unffile->Get(Variable+"RespMatSys_"+SystName+"_UP_deltaSys"+NameTag))->Clone("deltaup");
      TH1D* deltadown = (TH1D*) ((TH1D*)unffile->Get(Variable+"RespMatSys_"+SystName+"_DOWN_deltaSys"+NameTag))->Clone("deltadown");

      if(normXsec) {
        deltaup->Add(hnom);
        deltadown->Add(hnom);
        deltaup->Scale(1./deltaup->Integral());
        deltadown->Scale(1./deltadown->Integral());
        deltaup->Add(hnom_norm,-1);
        deltadown->Add(hnom_norm,-1);
      }

      for(int ibin = 1; ibin<=htot_reldelta->GetNbinsX(); ibin++){
        double bincontent_temp       = hnom->GetBinContent(ibin);
        if(normXsec) bincontent_temp = hnom_norm->GetBinContent(ibin);
        htot_reldelta->SetBinContent(ibin, (fabs(deltaup->GetBinContent(ibin) - deltadown->GetBinContent(ibin))*0.5)/bincontent_temp);
      }

    delete deltaup;
    delete deltadown;
    }

    // *******************
    // Alternative samples 
    // *******************

    else if((SystName=="MCATNLO") || (SystName=="POWHEG") || (SystName=="POWHEGHERWIG") || (SystName=="POWHEGV2") || (SystName=="POWHEGV2HERWIG") || (SystName=="AMCATNLOFXFX") || (SystName=="TOP_PT") || (SystName=="MADGRAPHMLM") || (SystName=="PERUGIA11") || (SystName=="SPINCORR") || (SystName.Contains("PDF")) ){
      
      // This is the fluctuation with dimension n_bins
      htot_reldelta = (TH1D*)((TH1D*)unffile->Get(Variable+"RespMatSys_"+SystName+"_deltaSys"+NameTag))->Clone(Variable+"RespMatSys_"+SystName+"_deltaSysRel"+NameTag+(normXsec?"_Norm":""));

      if(normXsec) {
        htot_reldelta->Add(hnom);
        htot_reldelta->Scale(1./htot_reldelta->Integral());
        htot_reldelta->Add(hnom_norm,-1);
      }

      for(int ibin = 1; ibin<=htot_reldelta->GetNbinsX(); ibin++){
        double bincontent_temp       = hnom->GetBinContent(ibin);
        if(normXsec) bincontent_temp = hnom_norm->GetBinContent(ibin);
        htot_reldelta->SetBinContent(ibin, fabs(htot_reldelta->GetBinContent(ibin)/bincontent_temp));
      }
    }

    // ***********
    // Backgrounds
    // ***********

    else if( (SystName=="DYeemm") || (SystName=="DYtautau") || (SystName=="singletop")  || (SystName=="ww") || (SystName=="wz") || (SystName=="zz") || (SystName=="ttbarW") || (SystName=="ttbarZ") || (SystName=="ttbarbg") || (SystName=="wtolnu") || (SystName=="other") ) {
      
      // This is the fluctuation with dimension n_bins
      htot_reldelta = (TH1D*)((TH1D*)unffile->Get(Variable+"BkgSys_"+SystName+"_deltaSys"+NameTag))->Clone(Variable+"BkgSys_"+SystName+"_deltaSysRel"+NameTag+(normXsec?"_Norm":""));

      if(normXsec) {
        htot_reldelta->Add(hnom);
        htot_reldelta->Scale(1./htot_reldelta->Integral());
        htot_reldelta->Add(hnom_norm,-1);
      }

      for(int ibin = 1; ibin<=htot_reldelta->GetNbinsX(); ibin++){
        double bincontent_temp       = hnom->GetBinContent(ibin);
        if(normXsec) bincontent_temp = hnom_norm->GetBinContent(ibin);
        htot_reldelta->SetBinContent(ibin, fabs(htot_reldelta->GetBinContent(ibin)/bincontent_temp));
      }
    }
  
    // ******************
    // Flat uncertainties
    // ******************

    else if( (SystName=="Lumi") || (SystName=="BR") ) {
      double flat_uncertainty = 0.025;
      if(SystName=="BR") flat_uncertainty = 0.015;

      // This is the fluctuation with dimension n_bins
      htot_reldelta = (TH1D*) hnom->Clone(Variable+"FlatSys_"+SystName+"_deltaSysRel"+NameTag+(normXsec?"_Norm":""));

      if(normXsec) {
        htot_reldelta->Reset();
      }
      else {
        for(int ibin = 1; ibin<=htot_reldelta->GetNbinsX(); ibin++){
          htot_reldelta->SetBinContent(ibin, flat_uncertainty);
        }
      }
    }

    else
    {
      std::cout<<"SystName is not one of the options defined in MatrixUnfControl::GetDeltasAllVars"<<std::endl;
      return;
    }

    for(int ibins = 1; ibins<=hnom->GetNbinsX(); ibins++){
      double bincontent_temp       = hnom->GetBinContent(ibins);
      if(normXsec) bincontent_temp = hnom_norm->GetBinContent(ibins);

      hnom_AllVar->SetBinContent(bin_counter,bincontent_temp);
      htot_reldelta_AllVar->SetBinContent(bin_counter,htot_reldelta->GetBinContent(ibins));

      ++bin_counter;
    }

    delete hnom;
    if(normXsec) delete hnom_norm;
    delete htot_reldelta;

    unffile->Close();
    delete unffile;
  }
  std::cout<<"Finished MatrixUnfControl::GetDeltasAllVars."<<std::endl;
  return;
}


// **********************************
// GetCombinedCovarianceMatrixAllVars
// **********************************

TH2D* MatrixUnfControl::GetCombinedCovarianceMatrixAllVars(TH1D *hnom_AllVar, TH1D *htot_reldelta_AllVar, TString Channel, TString SystName, TString TUNFFILE, Bool_t normXsec, TString NameTag ){

  // Calculate the correlation matrix for the sum of sub-components, 
  // then convert this to a covariance matrix using the envelope results
  std::cout << "Starting MatrixUnfControl::GetCombinedCovarianceMatrixAllVars" << std::endl;

  // Amandeep : Hardcoded here
  const int nbinsrhoi    = 24;  
  // const int nbinsrhoi = 6;

  TH1D* htot_absdelta_AllVar = (TH1D*) htot_reldelta_AllVar->Clone("htot_absdelta_AllVar");

  for(int ibins = 1; ibins<=htot_reldelta_AllVar->GetNbinsX(); ibins++){
    htot_absdelta_AllVar->SetBinContent(ibins, htot_reldelta_AllVar->GetBinContent(ibins)*hnom_AllVar->GetBinContent(ibins) );
  }


  std::vector<TString> syst;
  bool isBkgSyst = kFALSE;

  if(SystName=="SCALE")    syst = { "MESCALE_UP", "MESCALE_DOWN", "MEFACSCALE_UP", "MEFACSCALE_DOWN", "MERENSCALE_UP", "MERENSCALE_DOWN" };
  if(SystName=="BFRAG")    syst = { "BFRAG_UP", "BFRAG_DOWN", "BFRAG_PETERSON" };
  if(SystName=="COLORREC") syst = { "ERDON", "ERDONRETUNE", "GLUONMOVETUNE" };

  // Amandeep : Samples with UP and DOWN variations 
  // Jason mentioned not enveloping BTAG and BTAG_LJET

  if( (SystName=="JER")    || (SystName.Contains("JES")) || (SystName=="PU")         || (SystName=="LEPT") || 
      (SystName=="KIN")    || (SystName=="BSEMILEP")     || (SystName=="UNCLUSTERED")|| 
      (SystName=="UETUNE") || (SystName=="MASS")         || (SystName=="MATCH")      || (SystName=="PDF_ALPHAS") || 
      (SystName=="L1PREFIRING")    || (SystName=="ELE_ID")         || (SystName=="ELE_RECO") || 
      (SystName=="ELE_SCALE_GAIN") || (SystName=="ELE_SCALE_SYST") || (SystName=="ELE_SCALE_STAT") || 
      (SystName=="BTAG_CORR")      || (SystName=="BTAG_UNCORR")    || (SystName=="BTAG_LJET_CORR") || (SystName=="BTAG_LJET_UNCORR") || 
      (SystName=="MUON_ID")        || (SystName=="MUON_ISO")       || (SystName=="MUON_SCALE")     || (SystName=="PSSCALE_WEIGHT_6") || 
      (SystName=="PSSCALE_WEIGHT_7")) 
      syst = {SystName+"_UP",SystName+"_DOWN"};
  
  else if(SystName.Contains("PDF")) syst = {SystName};

  // *******************
  // Alternative samples 
  // *******************

  if((SystName=="MCATNLO") || (SystName=="POWHEG") || (SystName=="POWHEGHERWIG") || (SystName=="POWHEGV2") || (SystName=="POWHEGV2HERWIG") || (SystName=="AMCATNLOFXFX") || (SystName=="TOP_PT") || (SystName=="MADGRAPHMLM") || (SystName=="PERUGIA11") || (SystName=="SPINCORR")) syst = {SystName};
  
  // ***********
  // Backgrounds
  // ***********

  if( (SystName=="DYeemm") || (SystName=="DYtautau") || (SystName=="singletop")  || (SystName=="ww") || (SystName=="wz") || (SystName=="zz") || (SystName=="ttbarW") || (SystName=="ttbarZ") || (SystName=="ttbarbg") || (SystName=="other") ) {
    syst = {SystName};
    isBkgSyst = kTRUE;
  }


  double scalingfactor = 1.;
  if(SystName=="MASS") scalingfactor = 1.;          // convert from +/-3 GeV to +/-1 GeV
  if(SystName.Contains("PDF")) scalingfactor = 10.; // convert from sum in quadrature to RMS (for 100 replicas)
  if(SystName.Contains("PDF_ALPHAS")) scalingfactor = 1.;

  const int xbins = htot_absdelta_AllVar->GetNbinsX();
  const int ybins = htot_absdelta_AllVar->GetNbinsX();

  TMatrixD systCov(xbins,ybins);

  // Over the systematics
  for(unsigned int isyst = 0;isyst < syst.size();isyst++){

    TH1D* hdelta_AllVar = new TH1D(Form("hdelta_AllVar%s%s",(normXsec?"Norm":""),NameTag.Data()), Form("hdelta_AllVar%s%s",(normXsec?"Norm":""),NameTag.Data()), fVariableNames.size()*nbinsrhoi,0,fVariableNames.size()*nbinsrhoi);

    double componentscalingfactor = 1.;
    // if(SystName=="SCALE" && syst[isyst].BeginsWith("PSFSR") ) componentscalingfactor=sqrt(2.); // applied only to FSR uncertainties

    int bin_counter = 1;

    // Over the variables
    for(unsigned int i=0; i < fVariableNames.size(); i++)
    {
      TString Variable = fVariableNames.at(i);

      if(Variable=="") continue;

      TString TUnfoldResultsFilename = Form(TUNFFILE+Variable+".root");
      TFile *unffile                 = new TFile(TUnfoldResultsFilename.Data(),"READ");

      if(unffile==NULL) exit(0);

      // This is the 1D vector with fluctuations of dimension n_bins
      TH1D *delta = (TH1D*) ((TH1D*)unffile->Get(Variable+(isBkgSyst?"BkgSys_":"RespMatSys_")+syst[isyst]+"_deltaSys"+NameTag))->Clone("delta");

      if(normXsec){
        TH1D *hnom      = (TH1D*) ((TH1D*)unffile->Get(Variable+"TUnfResultCor"+NameTag))->Clone("hnom");
        TH1D *hnom_norm = (TH1D*) hnom->Clone("hnom_norm");
        hnom_norm->Scale(1./hnom_norm->Integral());

        delta->Add(hnom);
        delta->Scale(1./delta->Integral());
        delta->Add(hnom_norm,-1);

       delete hnom;
       delete hnom_norm;
      }

      for(int ibins = 1; ibins <= delta->GetNbinsX(); ibins++){
        hdelta_AllVar->SetBinContent(bin_counter,delta->GetBinContent(ibins));
        ++bin_counter;
      }

      delete delta;
      unffile->Close();
      delete unffile;

    }

    for(int xbin=1; xbin <= xbins; xbin++){
      for(int ybin=1; ybin <= ybins; ybin++){
       systCov(xbin-1,ybin-1) += (hdelta_AllVar->GetBinContent(xbin)*hdelta_AllVar->GetBinContent(ybin)/componentscalingfactor/componentscalingfactor);
      }
    }

    delete hdelta_AllVar;
  }

  // ******************
  // Flat uncertainties
  // ******************
  
  // Convert to correlation matrix (or fill fully-correlated correlation matrix for the flat systematics)
  
  if( (SystName=="Lumi") || (SystName=="BR") ) {

    for(int xbin=1; xbin<=xbins; xbin++){
      for(int ybin=1; ybin<=ybins; ybin++){
        systCov(xbin-1,ybin-1) = 1;
      }
    }

  }
  else {
    TVectorD diag_nozeros = TMatrixDDiag(systCov);
    for(int xbin=1; xbin<=xbins; xbin++)
    {
      if( diag_nozeros(xbin-1) == 0 ) diag_nozeros(xbin-1) = 1e-20;
    }
    systCov.NormByDiag(diag_nozeros);
  }

  // Convert back to covariance matrix, but using the bin errors from the envelope method
  for(int xbin=1; xbin<=xbins; xbin++){
    for(int ybin=1; ybin<=ybins; ybin++){
      systCov(xbin-1,ybin-1)*= fabs( htot_absdelta_AllVar->GetBinContent(xbin)*htot_absdelta_AllVar->GetBinContent(ybin) );
    }
  }

  TH2D *systCovMatrix = new TH2D(Form("%sSystCovMatrix_AllVar%s%s",SystName.Data(),(normXsec?"Norm":""),NameTag.Data()), Form("%sSystCovMatrix_AllVar%s%s",SystName.Data(),(normXsec?"Norm":""),NameTag.Data()), fVariableNames.size()*nbinsrhoi,0,fVariableNames.size()*nbinsrhoi, fVariableNames.size()*nbinsrhoi,0,fVariableNames.size()*nbinsrhoi);
  for(int xbin=1; xbin<=systCovMatrix->GetNbinsX(); xbin++){
    for(int ybin=1; ybin<=systCovMatrix->GetNbinsY(); ybin++){
      systCovMatrix->SetBinContent(xbin,ybin,systCov(xbin-1,ybin-1)/scalingfactor/scalingfactor);
    }
  }

  // The envelope method does not preserve covariance characteristics of a normalised differential xsec 
  // (e.g, the error on the sum of all bins should equal zero). 
  // Fix this by recalculating each bin as bin/(sum of 6 bins):
  if(normXsec) MatrixUnf::GetNormalisedCovarianceMatrixAllVars(hnom_AllVar, systCovMatrix, systCovMatrix);

  delete htot_absdelta_AllVar;
  return systCovMatrix;
}



// Draw official labels (CMS Preliminary, luminosity and CM energy) above plot
void MatrixUnfControl::DrawCMSLabels(double textSize, double lumi, int cmsprelim, double energy) {

  const char *text = (fEra == "fullRun2UL") ? "%.0f fb^{-1} (%2.f TeV)" : "%2.1f fb^{-1} (%2.f TeV)";
  //const char *text = "%2.1f fb^{-1} (%2.f TeV)";
  
  TPaveText *label = new TPaveText();
  label->SetX1NDC(gStyle->GetPadLeftMargin());
  label->SetY1NDC(0.98-gStyle->GetPadTopMargin());
  label->SetX2NDC(1.0-gStyle->GetPadRightMargin() + 0.04);
  label->SetY2NDC(0.98); //##
  label->SetTextFont(42);
  label->AddText(Form(text, lumi/1000, energy)); //## orig
  label->SetFillStyle(0);
  label->SetBorderSize(0);
  if (textSize!=0) label->SetTextSize(textSize);
    label->SetTextAlign(32);
    label->Draw("same");
    
    TPaveText *cms = new TPaveText();
    cms->AddText("CMS");
    //Official
    //cms->SetX1NDC(	  gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
    //cms->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );
    //cms->SetX2NDC(	  gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
    //cms->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength()        );

    //Semi-official
    //cms->SetX1NDC(0.481 - gStyle->GetPadLeftMargin()  );
    //cms->SetX1NDC( gStyle->GetPadLeftMargin() + gStyle->GetTickLength()  );
    cms->SetX1NDC( gStyle->GetPadLeftMargin() - 0.02 );
    cms->SetY1NDC(0.98  - gStyle->GetPadTopMargin() + 0.005  );
    cms->SetX2NDC( gStyle->GetPadLeftMargin()  + 0.18 );
    cms->SetY2NDC(0.98 + 0.005                               );
    
    //std::cout << "############# " << gStyle->GetPadLeftMargin() << " " << gStyle->GetPadTopMargin() << " " <<  gStyle->GetPadRightMargin() << std::endl;
    
    cms->SetFillStyle(0);
    cms->SetBorderSize(0);
    if (textSize!=0) cms->SetTextSize(textSize*1.1);
    cms->SetTextAlign(22);
    cms->SetTextFont(61);
    cms->Draw("same");

    if(cmsprelim > 0) {
      TPaveText *extra = new TPaveText();
      if(cmsprelim == 2) { extra->AddText("Private Work"); }
      else { extra->AddText("Work in Progress"); }
      
      //Official
      //extra->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
      //extra->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.10 );
      //extra->SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
      //extra->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );

      //Semi-official
      //extra->SetX1NDC(0.47 - gStyle->GetPadLeftMargin() + gStyle->GetTickLength()     );
      extra->SetX1NDC( gStyle->GetPadLeftMargin() + .22 );
      extra->SetY1NDC(0.9775 - gStyle->GetPadTopMargin() );
      extra->SetX2NDC( gStyle->GetPadLeftMargin() + .38 );
      extra->SetY2NDC(0.9775            );

      extra->SetFillStyle(0);
      extra->SetBorderSize(0);
      if (textSize!=0) extra->SetTextSize(textSize);
      extra->SetTextAlign(22);
      extra->SetTextFont(52);
      extra->Draw("same");
    }
}


// Draw label for Decay Channel in upper left corner of plot
void MatrixUnfControl::DrawDecayChLabel(TString decaychannel, double textSize) {

    TPaveText *decch = new TPaveText();

    decch->AddText(decaychannel);

    decch->SetX1NDC(	  gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
    decch->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.02 );
    //decch->SetX2NDC(	  gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
    decch->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() + 0.04 );
    decch->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() );
    //decch->SetY1NDC(0.975 - gStyle->GetPadTopMargin() );
    //decch->SetY2NDC(0.975            );

    decch->SetFillStyle(0);
    decch->SetBorderSize(0);
    if (textSize!=0) decch->SetTextSize(textSize);
    decch->SetTextAlign(12);
    decch->Draw("same");
}

void MatrixUnfControl::setResultLegendStyle(TLegend *leg, TString name)
{
    double height = 0.095, width = 0.2;
    double x1 = 1-width-gStyle->GetPadRightMargin()-gStyle->GetTickLength()*0.9, y1 = 1-height-gStyle->GetPadTopMargin()-gStyle->GetTickLength()*0.75;

    // Move all legends to the left
    x1  = gStyle->GetPadLeftMargin()+gStyle->GetTickLength()*0.7;
    y1 -= height/2.;
    height += height/2.;

    if(name.Contains("llbar_delta_phi") || name.Contains("large")) {
      y1 -= 0.7*3.*height/2.;
      height += 0.7*3.*height/2.;
    }

    leg->SetX1NDC(x1);
    leg->SetY1NDC(y1);
    leg->SetX2NDC(x1 + width);
    leg->SetY2NDC(y1 + height);

    leg->SetTextFont(42);
    leg->SetTextAlign(12);
    leg->SetTextSize(0.035);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
}

void MatrixUnfControl::setResultLegendStyleN(TLegend *leg, TString name)
{
    double height = 0.095, width = 0.2;
    double x1 = 1-0.3-gStyle->GetPadRightMargin()-gStyle->GetTickLength(), y1 = 1-height-gStyle->GetPadTopMargin()-gStyle->GetTickLength()*0.75;

    if(name.Contains("large")) {

      y1 -= height/2.;
      height += height/2.;

      y1 -= 0.7*3.*height/2.;
      height += 0.7*3.*height/2.;
    }


    leg->SetX1NDC(x1);
    leg->SetY1NDC(y1);
    leg->SetX2NDC(x1 + width);
    leg->SetY2NDC(y1 + height);

    leg->SetTextFont(42);
    leg->SetTextAlign(12);
    leg->SetTextSize(0.035);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
}

void MatrixUnfControl::setTheoryStyleAndFillLegend(TH1* histo, TString theoryName, TLegend *leg){

    histo->GetXaxis()->SetTitleOffset(1.08);
    histo->GetXaxis()->SetTitleSize(0.05);
    histo->GetXaxis()->SetLabelFont(42);
    histo->GetXaxis()->SetLabelOffset(0.0075);
    histo->GetXaxis()->SetLabelSize(0.045);
    histo->GetXaxis()->SetTickLength(0.02);

    histo->GetYaxis()->SetTitleOffset(1.64);
    histo->GetYaxis()->SetTitleSize(0.049);
    histo->GetYaxis()->SetLabelFont(42);
    histo->GetYaxis()->SetLabelOffset(0.007);
    histo->GetYaxis()->SetLabelSize(0.045);

    histo->SetLineWidth(2);
    if(theoryName != "data"){
        histo->SetMarkerSize(0);
        histo->SetMarkerStyle(1);
    }

    if(theoryName == "nominal"){
        histo->SetLineColor(kRed+1);
        histo->SetMarkerColor(kRed+1);
        histo->SetLineStyle(1);
        //if(leg) leg->AddEntry(histo, "PowhegV2+P8",  "l");
        //if(leg) leg->AddEntry(histo, "Powhegv2 + Pythia8",  "l");
        if(leg) leg->AddEntry(histo, "POWHEGV2 + PYTHIA8",  "l");
    }
    if(theoryName == "POWHEGvars"){
        histo->SetLineColor(kRed+1);
        histo->SetMarkerColor(kRed+1);
        histo->SetLineStyle(1);
    }
    if(theoryName == "uncorr"){
        histo->SetLineColor(kGreen+1);
        histo->SetMarkerColor(kGreen+1);
        histo->SetLineStyle(7);
        if(leg) leg->AddEntry(histo, "NLO, uncorrelated",  "l");
    }
    if(theoryName == "uncorrscaleup"){
        histo->SetLineColor(kGreen+1);
        histo->SetMarkerColor(kGreen+1);
        histo->SetLineStyle(8);
        if(leg) leg->AddEntry(histo, "uncorr, scale up",  "l");
    }
    if(theoryName == "uncorrscaledown"){
        histo->SetLineColor(kGreen+1);
        histo->SetMarkerColor(kGreen+1);
        histo->SetLineStyle(9);
        if(leg) leg->AddEntry(histo, "uncorr, scale down",  "l");
    }
    if(theoryName == "unpol"){
        histo->SetLineColor(kGreen+1);
        histo->SetMarkerColor(kGreen+1);
        histo->SetLineStyle(7);
        if(leg) leg->AddEntry(histo, "NLO, unpolarized",  "l");
    }
    if(theoryName == "NLOW"){
        histo->SetLineColor(kBlue);
        histo->SetMarkerColor(kBlue);
        histo->SetLineStyle(2);
        //if(leg) leg->AddEntry(histo, "aMC@NLOFxFx+P8",  "l");
        if(leg) leg->AddEntry(histo, "NLO, SM",  "l");
    }
    if(theoryName == "NLOWscaleup"){
        histo->SetLineColor(kBlue);
        histo->SetMarkerColor(kBlue);
        histo->SetLineStyle(3);
        //if(leg) leg->AddEntry(histo, "aMC@NLOFxFx+P8",  "l");
        if(leg) leg->AddEntry(histo, "NLO, SM, scale up",  "l");
    }
    if(theoryName == "NLOWscaledown"){
        histo->SetLineColor(kBlue);
        histo->SetMarkerColor(kBlue);
        histo->SetLineStyle(4);
        //if(leg) leg->AddEntry(histo, "aMC@NLOFxFx+P8",  "l");
        if(leg) leg->AddEntry(histo, "NLO, SM, scale down",  "l");
    }
    if(theoryName.Contains("NNLO")){
        histo->SetLineColor(kBlack);
        histo->SetMarkerColor(kBlack);
        histo->SetLineStyle(3);
        //if(leg) leg->AddEntry(histo, "aMC@NLOFxFx+P8",  "l");
        if(leg) leg->AddEntry(histo, "NNLO, SM",  "l");
    }
    if(theoryName == "amcatnlo"){
        histo->SetLineColor(kViolet-2);
        histo->SetMarkerColor(kViolet-2);
        histo->SetLineStyle(5);
        //if(leg) leg->AddEntry(histo, "aMC@NLOFxFx + P8",  "l");
        if(leg) leg->AddEntry(histo, "MG5_aMC@NLO + PYTHIA8 [FxFx]",  "l");
    }

    if(theoryName == "toppt"){
        histo->SetLineColor(kRed-5);
        histo->SetMarkerColor(kRed-5);
        histo->SetLineStyle(3);
        //if(leg) leg->AddEntry(histo, "MadGraphMLM+P8",  "l");
        if(leg) leg->AddEntry(histo, "Top p_{T} reweighting",  "l");
    }
    if(theoryName == "powhegv2herwig"){
        histo->SetLineColor(kAzure+6);
        histo->SetMarkerColor(kAzure+6);
        histo->SetLineStyle(4);
        //if(leg) leg->AddEntry(histo, "PowhegV2+H++",  "l");
        if(leg) leg->AddEntry(histo, "Powhegv2 + Herwig++",  "l");
    }
    if(theoryName == "ahrens"){
        histo->SetLineColor(kViolet-6);
        histo->SetMarkerColor(kViolet-6);
        histo->SetLineStyle(6);
        if(leg) leg->AddEntry(histo, "NLO+NNLL",  "l");
    }
    if(theoryName == "kidonakis"){
        histo->SetLineColor(kViolet-6);
        histo->SetMarkerColor(kViolet-6);
        histo->SetLineStyle(2);
        if(leg) leg->AddEntry(histo, "Approx. NNLO",  "l");
    }
    if(theoryName == "matchup" || theoryName == "mass175.5"){
        histo->SetLineStyle(7);
        histo->SetLineColor(kPink-7);
        histo->SetMarkerColor(kPink-7);
        if(leg && theoryName == "matchup") leg->AddEntry(histo,"Matching up",  "l");
        if(leg && theoryName == "mass175.5") leg->AddEntry(histo,"Mass = 175.5 GeV",  "l");
    }
    if(theoryName == "matchdown" || theoryName == "mass169.5"){
        histo->SetLineStyle(7);
        histo->SetLineColor(kRed-7);
        histo->SetMarkerColor(kRed-7);
        if(leg && theoryName == "matchdown") leg->AddEntry(histo,"Matching down",  "l");
        if(leg && theoryName == "mass169.5")   leg->AddEntry(histo,"Mass = 169.5 GeV",  "l");
    }
    if(theoryName == "scaleup" || theoryName == "mass173.5"){
        histo->SetLineStyle(2);
        histo->SetLineColor(kAzure+2);
        histo->SetMarkerColor(kAzure+2);
        if(leg && theoryName == "scaleup") leg->AddEntry(histo,"4*Q^{2}",  "l");
        if(leg && theoryName == "mass173.5") leg->AddEntry(histo,"Mass = 173.5 GeV",  "l");
    }
    if(theoryName == "scaledown" || theoryName == "mass171.5"){
        histo->SetLineStyle(2);
        histo->SetLineColor(8);
        histo->SetMarkerColor(8);
        if(leg && theoryName == "scaledown") leg->AddEntry(histo,"Q^{2}/4",  "l");
        if(leg && theoryName == "mass171.5")   leg->AddEntry(histo,"Mass = 171.5 GeV",  "l");
    }
    if(theoryName == "mass178.5"){
        histo->SetLineStyle(3);
        histo->SetLineColor(42);
        histo->SetMarkerColor(42);
        if(leg && theoryName == "matchup") leg->AddEntry(histo,"Matching up",  "l");
        if(leg && theoryName == "mass178.5") leg->AddEntry(histo,"Mass = 178.5 GeV",  "l");
    }
    if(theoryName == "mass166.5"){
        histo->SetLineStyle(3);
        histo->SetLineColor(48);
        histo->SetMarkerColor(48);
        if(leg && theoryName == "matchdown") leg->AddEntry(histo,"Matching down",  "l");
        if(leg && theoryName == "mass166.5")   leg->AddEntry(histo,"Mass = 166.5 GeV",  "l");
    }
}

// Amandeep : Adding for generalized definition of bin width
// In the N-D case, tunfresult->GetBinWidth() doesn't quite work
// This loads in the final binning for the variable and then gives the correct size

// PS : This should be memoized, shouldn't have to open the XML file each time
Double_t MatrixUnfControl::GetPhysicalBinWidth(TString VariableName, int BinNumber)
{
  // Read in binning schemes in XML format
  TDOMParser parser;
  Int_t error = parser.ParseFile("binning/" + VariableName + "_binning_rebinnedB.xml");
  if(error) std::cout << "error=" << error << " from TDOMParser\n";

  TXMLDocument const* XMLdocument     = parser.GetXMLDocument();
  TUnfoldBinningXML* generatorBinning = TUnfoldBinningXML::ImportXML(XMLdocument,"generator_rebinnedB");

  if(!generatorBinning)  std::cout<<"could not read 'generator_rebinnedB' binning\n";

  const TUnfoldBinning *binning = generatorBinning->FindNode("ttbargen_rebinnedB");
  
  // GetBinSize() automatically gets the N-D size of the bin
  return binning->GetBinSize(BinNumber);
}

void MatrixUnfControl::SetAxis(TString VariableName, TAxis *Axis) {
  if(VariableName.Contains("llbar_delta_phi")) {
    Axis->Set(Axis->GetNbins(), 0.0, TMath::Pi());
  }
  else if(VariableName.Contains("llbar_delta_eta")) {
    Axis->Set(Axis->GetNbins(), 0.0, 5.0);
  }
  else {
    Axis->Set(Axis->GetNbins(), -1.0, 1.0);
  }
}

void MatrixUnfControl::Draw2DLabels(TString VariableName, double y1_NDC, double text_size)
{

  TString dimlabeltext = "";
  TString dimlabel_1text = "";
  TString dimlabel_2text = "";
  TString dimlabel_3text = "";
  TString dimlabel_4text = "";

  if (VariableName.Contains("_mttbar")) {
    dimlabel_1text = "0 < m_{t#bar{t}} #leq 450";
    dimlabel_2text = "450 < m_{t#bar{t}} #leq 600";
    dimlabel_3text = "600 < m_{t#bar{t}} #leq 800"; 
    dimlabel_4text = "m_{t#bar{t}} > 800"; 
  }
  else if (VariableName.Contains("_vs_ScatteringAngle_TTBarFrame")) {
    dimlabel_1text = "-1.0 #leq cos#Theta #leq -0.5";
    dimlabel_2text = "-0.5 < cos#Theta #leq 0.0";
    dimlabel_3text = "0.0 < cos#Theta #leq 0.5"; 
    dimlabel_4text = "0.5 < cos#Theta #leq 1.0";
  }
  else if (VariableName.Contains("_vs_ToppT")) {
    dimlabel_1text = "0 < p_{T}^{t} #leq 80";
    dimlabel_2text = "80 < p_{T}^{t} #leq 150";
    dimlabel_3text = "150 < p_{T}^{t} #leq 250"; 
    dimlabel_4text = "p_{T}^{t} > 250"; 
  }
  else if (VariableName.Contains("_vs_ExtraJets")) {
    dimlabel_1text = "N_{jets} = 2";
    dimlabel_2text = "N_{jets} = 3";
    dimlabel_3text = "N_{jets} = 4"; 
    dimlabel_4text = "N_{jets} #geq 5"; 
  }

  double bin2D = (1 - gStyle->GetPadLeftMargin() - gStyle->GetPadRightMargin())/4;

  TPaveText *dimlabel_1 = new TPaveText();
  dimlabel_1->AddText(dimlabel_1text);
  dimlabel_1->SetX1NDC(gStyle->GetPadLeftMargin() + 0*bin2D);
  dimlabel_1->SetY1NDC(y1_NDC);
  dimlabel_1->SetX2NDC(gStyle->GetPadLeftMargin() + 1*bin2D);
  dimlabel_1->SetY2NDC(y1_NDC + 0.05);
  dimlabel_1->SetFillStyle(0);
  dimlabel_1->SetBorderSize(0);
  dimlabel_1->SetTextSize(text_size);
  dimlabel_1->SetTextFont(62);
  dimlabel_1->SetTextAlign(22);
  dimlabel_1->Draw("same");

  TPaveText *dimlabel_2 = new TPaveText();
  dimlabel_2->AddText(dimlabel_2text);
  dimlabel_2->SetX1NDC(gStyle->GetPadLeftMargin() + 1*bin2D);
  dimlabel_2->SetY1NDC(y1_NDC);
  dimlabel_2->SetX2NDC(gStyle->GetPadLeftMargin() + 2*bin2D);
  dimlabel_2->SetY2NDC(y1_NDC + 0.05);
  dimlabel_2->SetFillStyle(0);
  dimlabel_2->SetBorderSize(0);
  dimlabel_2->SetTextSize(text_size);
  dimlabel_2->SetTextFont(62);
  dimlabel_2->SetTextAlign(22);
  dimlabel_2->Draw("same");

  TPaveText *dimlabel_3 = new TPaveText();
  dimlabel_3->AddText(dimlabel_3text);
  dimlabel_3->SetX1NDC(gStyle->GetPadLeftMargin() + 2*bin2D);
  dimlabel_3->SetY1NDC(y1_NDC);
  dimlabel_3->SetX2NDC(gStyle->GetPadLeftMargin() + 3*bin2D);
  dimlabel_3->SetY2NDC(y1_NDC + 0.05);
  dimlabel_3->SetFillStyle(0);
  dimlabel_3->SetBorderSize(0);
  dimlabel_3->SetTextSize(text_size);
  dimlabel_3->SetTextFont(62);
  dimlabel_3->SetTextAlign(22);
  dimlabel_3->Draw("same");

  TPaveText *dimlabel_4 = new TPaveText();
  dimlabel_4->AddText(dimlabel_4text);
  dimlabel_4->SetX1NDC(gStyle->GetPadLeftMargin() + 3*bin2D);
  dimlabel_4->SetY1NDC(y1_NDC);
  dimlabel_4->SetX2NDC(gStyle->GetPadLeftMargin() + 4*bin2D);
  dimlabel_4->SetY2NDC(y1_NDC + 0.05);
  dimlabel_4->SetFillStyle(0);
  dimlabel_4->SetBorderSize(0);
  dimlabel_4->SetTextSize(text_size);
  dimlabel_4->SetTextFont(62);
  dimlabel_4->SetTextAlign(22);
  dimlabel_4->Draw("same");

}

void MatrixUnfControl::Draw2DLabelsVertical(TString VariableName, double x1_NDC, double text_size)
{

  TString dimlabeltext = "";
  TString dimlabel_1text = "";
  TString dimlabel_2text = "";
  TString dimlabel_3text = "";
  TString dimlabel_4text = "";

  if (VariableName.Contains("_mttbar")) {
    dimlabel_1text = "450 #geq m_{t#bar{t}} > 0";
    dimlabel_2text = "600 #geq m_{t#bar{t}} > 450";
    dimlabel_3text = "800 #geq m_{t#bar{t}} > 600"; 
    dimlabel_4text = "800 < m_{t#bar{t}}"; 
  }
  else if (VariableName.Contains("_vs_ScatteringAngle_TTBarFrame")) {
    dimlabel_1text = "-0.5 #geq cos#Theta #geq -1.0";
    dimlabel_2text = "0.0 #geq cos#Theta > -0.5";
    dimlabel_3text = "0.5 #geq cos#Theta > 0.0"; 
    dimlabel_4text = "1.0 #geq cos#Theta > 0.5";
  }
  else if (VariableName.Contains("_vs_ToppT")) {
    dimlabel_1text = "80 #geq p_{T}^{t} > 0";
    dimlabel_2text = "150 #geq p_{T}^{t} > 80";
    dimlabel_3text = "250 #geq p_{T}^{t} > 150"; 
    dimlabel_4text = " 250 < p_{T}^{t}"; 
  }
  else if (VariableName.Contains("_vs_ExtraJets")) {
    dimlabel_1text = "2 = N_{jets}";
    dimlabel_2text = "3 = N_{jets}";
    dimlabel_3text = "4 = N_{jets}"; 
    dimlabel_4text = "5 #leq N_{jets}"; 
  }

  double bin2D = (1 - gStyle->GetPadBottomMargin() - gStyle->GetPadTopMargin())/4;

  TPaveText *dimlabel_1 = new TPaveText();
  TText *dimlabel_1Text = dimlabel_1->AddText(dimlabel_1text);
  dimlabel_1Text->SetTextAngle(-90);
  dimlabel_1->SetX1NDC(x1_NDC);
  dimlabel_1->SetY1NDC(gStyle->GetPadBottomMargin() + 0*bin2D);
  dimlabel_1->SetX2NDC(x1_NDC + 0.05);
  dimlabel_1->SetY2NDC(gStyle->GetPadBottomMargin() + 1*bin2D);
  dimlabel_1->SetFillStyle(0);
  dimlabel_1->SetBorderSize(0);
  dimlabel_1->SetTextSize(text_size);
  dimlabel_1->SetTextFont(62);
  dimlabel_1->SetTextAlign(22);
  dimlabel_1->Draw("same");

  TPaveText *dimlabel_2 = new TPaveText();
  TText *dimlabel_2Text = dimlabel_2->AddText(dimlabel_2text);
  dimlabel_2Text->SetTextAngle(-90);
  dimlabel_2->SetX1NDC(x1_NDC);
  dimlabel_2->SetY1NDC(gStyle->GetPadBottomMargin() + 1*bin2D);
  dimlabel_2->SetX2NDC(x1_NDC + 0.05);
  dimlabel_2->SetY2NDC(gStyle->GetPadBottomMargin() + 2*bin2D);
  dimlabel_2->SetFillStyle(0);
  dimlabel_2->SetBorderSize(0);
  dimlabel_2->SetTextSize(text_size);
  dimlabel_2->SetTextFont(62);
  dimlabel_2->SetTextAlign(22);
  dimlabel_2->Draw("same");

  TPaveText *dimlabel_3 = new TPaveText();
  TText *dimlabel_3Text = dimlabel_3->AddText(dimlabel_3text);
  dimlabel_3Text->SetTextAngle(-90);
  dimlabel_3->SetX1NDC(x1_NDC);
  dimlabel_3->SetY1NDC(gStyle->GetPadBottomMargin() + 2*bin2D);
  dimlabel_3->SetX2NDC(x1_NDC + 0.05);
  dimlabel_3->SetY2NDC(gStyle->GetPadBottomMargin() + 3*bin2D);
  dimlabel_3->SetFillStyle(0);
  dimlabel_3->SetBorderSize(0);
  dimlabel_3->SetTextSize(text_size);
  dimlabel_3->SetTextFont(62);
  dimlabel_3->SetTextAlign(22);
  dimlabel_3->Draw("same");

  TPaveText *dimlabel_4 = new TPaveText();
  TText *dimlabel_4Text = dimlabel_4->AddText(dimlabel_4text);
  dimlabel_4Text->SetTextAngle(-90);
  dimlabel_4->SetX1NDC(x1_NDC);
  dimlabel_4->SetY1NDC(gStyle->GetPadBottomMargin() + 3*bin2D);
  dimlabel_4->SetX2NDC(x1_NDC + 0.05);
  dimlabel_4->SetY2NDC(gStyle->GetPadBottomMargin() + 4*bin2D);
  dimlabel_4->SetFillStyle(0);
  dimlabel_4->SetBorderSize(0);
  dimlabel_4->SetTextSize(text_size);
  dimlabel_4->SetTextFont(62);
  dimlabel_4->SetTextAlign(22);
  dimlabel_4->Draw("same");

}

void MatrixUnfControl::DrawCov2DAxisLabels(TString name) {

  //  double bin2D_horizontal = (1 - gStyle->GetPadLeftMargin() - gStyle->GetPadRightMargin())/4;
  //  double bin2D_vertical = (1 - gStyle->GetPadBottomMargin() - gStyle->GetPadTopMargin())/4;

  if(name.Contains("llbar_delta_phi") && !(name.Contains("_mttbar") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets"))) {

    TString binlabels[] = {"0","#pi/4","#pi/2","3#pi/4","#pi"};
    TLatex label_h;
    label_h.SetTextSize(0.03);
    label_h.SetTextFont(42);
    label_h.SetTextAlign(22);
    TLatex label_v;
    label_v.SetTextSize(0.03);
    label_v.SetTextFont(42);
    label_v.SetTextAlign(22);

    for (Int_t k=0;k<=4;k++) {
      label_h.DrawLatexNDC(0.11125+k*0.19,0.13,binlabels[k]);
      label_v.DrawLatexNDC(0.0875,0.15+k*0.1875,binlabels[k]);
    }
  }
  else if(name.Contains("llbar_delta_phi") && (name.Contains("_mttbar") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets"))) {

    TString binlabels[] = {"0","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi"};

    TLatex label_h;
    label_h.SetTextSize(0.03);
    label_h.SetTextFont(42);
    label_h.SetTextAlign(22);
    TLatex label_v;
    label_v.SetTextSize(0.03);
    label_v.SetTextFont(42);
    label_v.SetTextAlign(22);

    for (Int_t k=0;k<=8;k++) {
      label_h.DrawLatexNDC(0.11125+k*0.095,0.13,binlabels[k]);
      label_v.DrawLatexNDC(0.0875,0.15+k*0.09375,binlabels[k]);
    }
   
  }
  else if(name.Contains("llbar_delta_eta") && (name.Contains("_mttbar") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets"))) {

    TString binlabels[] = {"0","2.5","5","2.5","5","2.5","5","2.5","5"};

    TLatex label_h;
    label_h.SetTextSize(0.03);
    label_h.SetTextFont(42);
    label_h.SetTextAlign(22);
    TLatex label_v;
    label_v.SetTextSize(0.03);
    label_v.SetTextFont(42);
    label_v.SetTextAlign(22);

    for (Int_t k=0;k<=8;k++) {
      label_h.DrawLatexNDC(0.11125+k*0.095,0.13,binlabels[k]);
      label_v.DrawLatexNDC(0.0875,0.15+k*0.09375,binlabels[k]);
    }
   
  }

  else if((name.Contains("b1") || name.Contains("b2") || name.Contains("c_") || name.Contains("ll_c")) && (name.Contains("_mttbar") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets"))) {

    TString binlabels[] = {"-1","0","1","0","1","0","1","0","1"};

    TLatex label_h;
    label_h.SetTextSize(0.03);
    label_h.SetTextFont(42);
    label_h.SetTextAlign(22);
    TLatex label_v;
    label_v.SetTextSize(0.03);
    label_v.SetTextFont(42);
    label_v.SetTextAlign(22);

    for (Int_t k=0;k<=8;k++) {
      label_h.DrawLatexNDC(0.11125+k*0.095,0.13,binlabels[k]);
      label_v.DrawLatexNDC(0.0875,0.15+k*0.09375,binlabels[k]);
    }
   
  }


}

void MatrixUnfControl::DrawResMat2DAxisLabels(TString name) {

  //  double bin2D_horizontal = (1 - gStyle->GetPadLeftMargin() - gStyle->GetPadRightMargin())/4;
  //  double bin2D_vertical = (1 - gStyle->GetPadBottomMargin() - gStyle->GetPadTopMargin())/4;

  if(name.Contains("llbar_delta_phi") && !(name.Contains("_mttbar") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets"))) {

    TString binlabels_h[] = {"0","#pi/4","#pi/2","3#pi/4","#pi"};
    TString binlabels_v[] = {"0","#pi/8","#pi/4","3#pi/8","#pi/2","5#pi/8","3#pi/4","7#pi/8","#pi"};
    TLatex label_h;
    label_h.SetTextSize(0.03);
    label_h.SetTextFont(42);
    label_h.SetTextAlign(22);
    TLatex label_v;
    label_v.SetTextSize(0.03);
    label_v.SetTextFont(42);
    label_v.SetTextAlign(22);

    for (Int_t k=0;k<=4;k++) {
      label_h.DrawLatexNDC(0.11125+k*0.19,0.13,binlabels_h[k]);
    }

    for (Int_t k=0;k<=8;k++) {
      label_v.DrawLatexNDC(0.0875,0.15+k*0.09375,binlabels_v[k]);
    }

  }
  else if(name.Contains("llbar_delta_phi") && (name.Contains("_mttbar") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets"))) {

    TString binlabels_h[] = {"0","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi"};
    TString binlabels_v[] = {"0","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi"};

    TLatex label_h;
    label_h.SetTextSize(0.03);
    label_h.SetTextFont(42);
    label_h.SetTextAlign(22);
    TLatex label_v;
    label_v.SetTextSize(0.03);
    label_v.SetTextFont(42);
    label_v.SetTextAlign(22);

    for (Int_t k=0;k<=16;k++) {
      label_h.DrawLatexNDC(0.11125+k*0.0475,0.13,binlabels_h[k]);
    }
    for (Int_t k=0;k<=32;k++) {
      label_v.DrawLatexNDC(0.0875,0.15+k*0.0234375,binlabels_v[k]);
    }
   
  }
  else if(name.Contains("llbar_delta_eta") && (name.Contains("_mttbar") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets"))) {

    TString binlabels_h[] = {"0","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5"};
    TString binlabels_v[] = {"0","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5"};

    TLatex label_h;
    label_h.SetTextSize(0.03);
    label_h.SetTextFont(42);
    label_h.SetTextAlign(22);
    TLatex label_v;
    label_v.SetTextSize(0.03);
    label_v.SetTextFont(42);
    label_v.SetTextAlign(22);

    for (Int_t k=0;k<=16;k++) {
      label_h.DrawLatexNDC(0.11125+k*0.0475,0.13,binlabels_h[k]);
    }
    for (Int_t k=0;k<=32;k++) {
      label_v.DrawLatexNDC(0.0875,0.15+k*0.0234375,binlabels_v[k]);
    }
   
  }

  else if((name.Contains("b1") || name.Contains("b2") || name.Contains("c_") || name.Contains("ll_c")) && (name.Contains("_mttbar") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets"))) {

    TString binlabels_h[] = {"-1","0","1","0","1","0","1","0","1","0","1","0","1","0","1","0","1"};
    TString binlabels_v[] = {"-1","0","1","0","1","0","1","0","1","0","1","0","1","0","1","0","1","0","1","0","1","0","1","0","1","0","1","0","1","0","1","0","1"};

    TLatex label_h;
    label_h.SetTextSize(0.03);
    label_h.SetTextFont(42);
    label_h.SetTextAlign(22);
    TLatex label_v;
    label_v.SetTextSize(0.02);
    label_v.SetTextFont(42);
    label_v.SetTextAlign(22);

    for (Int_t k=0;k<=16;k++) {
      label_h.DrawLatexNDC(0.11125+k*0.0475,0.13,binlabels_h[k]);
    }
    for (Int_t k=0;k<=32;k++) {
      label_v.DrawLatexNDC(0.09,0.15+k*0.0234375,binlabels_v[k]);
    }
   
  }


}

void MatrixUnfControl::DrawTotEff2DAxisLabels(TString name) {

  if(name.Contains("llbar_delta_phi") && !(name.Contains("_mttbar") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets"))) {

    TString binlabels[] = {"0","#pi/4","#pi/2","3#pi/4","#pi"};
    TLatex label;
    label.SetTextSize(0.03);
    label.SetTextFont(42);
    label.SetTextAlign(22);

    for (Int_t k=0;k<=4;k++) {
      label.DrawLatexNDC(0.18+k*0.1925,0.13,binlabels[k]);
    }
  }
  else if(name.Contains("llbar_delta_phi") && (name.Contains("_mttbar") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets"))) {

    TString binlabels[] = {"0","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi"};
    TLatex label;
    label.SetTextSize(0.03);
    label.SetTextFont(42);
    label.SetTextAlign(22);

    for (Int_t k=0;k<=16;k++) {
      label.DrawLatexNDC(0.11125+k*0.0525,0.13,binlabels[k]);
    }
   
  }
  else if(name.Contains("llbar_delta_eta") && (name.Contains("_mttbar") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets"))) {

    TString binlabels[] = {"0","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5","2.5","5"};
    TLatex label;
    label.SetTextSize(0.03);
    label.SetTextFont(42);
    label.SetTextAlign(22);

    for (Int_t k=0;k<=16;k++) {
      label.DrawLatexNDC(0.11125+k*0.0525,0.13,binlabels[k]);
    }
   
  }

  else if((name.Contains("b1") || name.Contains("b2") || name.Contains("c_") || name.Contains("ll_c")) && (name.Contains("_mttbar") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets"))) {

    TString binlabels[] = {"-1","0","1","0","1","0","1","0","1","0","1","0","1","0","1","0","1"};
    TLatex label;
    label.SetTextSize(0.03);
    label.SetTextFont(42);
    label.SetTextAlign(22);

    for (Int_t k=0;k<=16;k++) {
      label.DrawLatexNDC(0.11125+k*0.0525,0.13,binlabels[k]);
    }
   
  }


}


void MatrixUnfControl::DrawPurStab2DAxisLabels(TString name) {

  if(name.Contains("llbar_delta_phi") && !(name.Contains("_mttbar") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets"))) {

    TString binlabels[] = {"0","#pi/4","#pi/2","3#pi/4","#pi"};
    TLatex label;
    label.SetTextSize(0.03);
    label.SetTextFont(42);
    label.SetTextAlign(22);

    for (Int_t k=0;k<=4;k++) {
      label.DrawLatexNDC(0.18+k*0.1925,0.13,binlabels[k]);
    }
  }
  else if(name.Contains("llbar_delta_phi") && (name.Contains("_mttbar") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets"))) {

    TString binlabels[] = {"0","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi","#pi/2","#pi"};
    TLatex label;
    label.SetTextSize(0.03);
    label.SetTextFont(42);
    label.SetTextAlign(22);

    for (Int_t k=0;k<=8;k++) {
      label.DrawLatexNDC(0.11125+k*0.105,0.13,binlabels[k]);
    }
   
  }
  else if(name.Contains("llbar_delta_eta") && (name.Contains("_mttbar") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets"))) {

    TString binlabels[] = {"0","2.5","5","2.5","5","2.5","5","2.5","5"};
    TLatex label;
    label.SetTextSize(0.03);
    label.SetTextFont(42);
    label.SetTextAlign(22);

    for (Int_t k=0;k<=8;k++) {
      label.DrawLatexNDC(0.11125+k*0.105,0.13,binlabels[k]);
    }
   
  }

  else if((name.Contains("b1") || name.Contains("b2") || name.Contains("c_") || name.Contains("ll_c")) && (name.Contains("_mttbar") || name.Contains("_vs_ScatteringAngle_TTBarFrame") || name.Contains("_vs_ToppT") || name.Contains("_vs_ExtraJets"))) {

    TString binlabels[] = {"-1","0","1","0","1","0","1","0","1"};
    TLatex label;
    label.SetTextSize(0.03);
    label.SetTextFont(42);
    label.SetTextAlign(22);

    for (Int_t k=0;k<=8;k++) {
      label.DrawLatexNDC(0.11125+k*0.105,0.13,binlabels[k]);
    }
   
  }


}
// End
