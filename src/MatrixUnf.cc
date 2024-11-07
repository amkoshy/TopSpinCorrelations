/************************************************************
File:         MatrixUnf.cc
Author:       Ajeeta Khatiwada, Andreas Jung, Jacob Linacre 
Description:  Class for unfolding. 
************************************************************/

#include "stdInc.h"
#include "MatrixUnf.h"
#include "TUnfoldDensity.h"

// root stuff
#include <iostream>
#include <stdlib.h>
#include "TMinuit.h"
#include <TDOMParser.h>
#include <TXMLDocument.h>
#include "TUnfoldBinningXML.h"

using namespace std;

//global variables
bool dbgUnf = kTRUE;

Bool_t fsetBiasDisTUnfold;
/*******************************************************************************
 * constructors -- just create the needed histograms. 
 *
 *
 */
MatrixUnf::MatrixUnf() {   //std-constructor needed for loading the class in root 

}

/*******************************************************************************
 * MatrixUnf -- constructor when filling this class and calculation the results: 
 *
 * Parameters: 1. Variable Name; 
 *             2. x-Axis Title of the resulting histos
 *             3. y-Axis Title 
 *             4. File Name (in .xml format) containing the binning for this variable will loaded here.
 **/ 
MatrixUnf::MatrixUnf(TString VariableName1, TString xAxisTitle, TString VariableName2, TString yAxisTitle, TString zAxisTitle, TString binFileName, Bool_t InputIs2D, Bool_t setBiasDisTUnfold, TString channel, Bool_t setSyst, TH1D* GenVec, Int_t rebinfine) {
    
  fVariableName1 = VariableName1;
  fVariableName2 = VariableName2;
  fsetBiasDisTUnfold = setBiasDisTUnfold;
  fXAxisTitle = xAxisTitle; 
  fYAxisTitle = yAxisTitle; 
  fZAxisTitle = zAxisTitle; 
  fsetSyst    = setSyst;
  frebinfine  = rebinfine;

  // Set the branching ratio for different channels
  
  // Amandeep : Changing the BR since we now include tau decays
  // Original :
  // if(channel=="ee") fBR=0.01147;
  // else if(channel=="mumu") fBR=0.01130;
  // else if(channel=="emu") fBR=0.02277;
  // else if(channel=="combined") fBR=0.01147+0.01130+0.02277;
  
  // Modified : Adding tau contributions
  // Numbers from PDG W decay column in Jasons excel sheet
  if     (channel=="ee")   fBR = 0.01147 + 0.00476;
  else if(channel=="emu")  fBR = 0.02277 + 0.00935;
  else if(channel=="mumu") fBR = 0.01130 + 0.00460;
  else if(channel=="combined") fBR = 0.01147 + 0.01130 + 0.02277 + 0.00476 + 0.00460 + 0.00935;
  // End

  this->SetNameTitle(VariableName1+"Unfold", VariableName1+" "+VariableName2+" MatrixUnf Objects "); 

  //create TObjArrays: 
  printf("Creating first TObjetArray\n");
  fHistUnfArray = new TObjArray(35);  
  
  // TObjArray *HistArray2;
  printf("Creating second TObjArray\n"); 
  fResultUnfHists = new TObjArray(50); 
  
  fBkgHistArray = new TObjArray(200);

  fResSysHist      = new TObjArray(250);
  fAequiResSysHist = new TObjArray(250);
  fTUnfDeltaSys    = new TObjArray(250);
  fTUnfDeltaSys_rebinnedA  = new TObjArray(250);
  fTUnfDeltaSys_rebinnedB  = new TObjArray(250);
  fTUnfCovMatSys           = new TObjArray(250);
  fTUnfCovMatSys_rebinnedA = new TObjArray(250);
  fTUnfCovMatSys_rebinnedB = new TObjArray(250);
  fTUnfCovMatSysNorm       = new TObjArray(250);
  fTUnfCovMatSysNorm_rebinnedA = new TObjArray(250);
  fTUnfCovMatSysNorm_rebinnedB = new TObjArray(250);


  // Read binning schemes in XML format
  TDOMParser parser;
  Int_t error = parser.ParseFile(binFileName.Data());
  if(error) cout<<"error="<<error<<" from TDOMParser\n";

  TXMLDocument const *XMLdocument = parser.GetXMLDocument();
  detectorBinning  = TUnfoldBinningXML::ImportXML(XMLdocument,"detector");
  generatorBinning = TUnfoldBinningXML::ImportXML(XMLdocument,"generator");

  if(detectorBinning) {
    //detectorBinning->PrintStream(cout,1,1);
  } else {
    cout<<"could not read 'detector' binning\n";
  }
  if(generatorBinning) {
    //generatorBinning->PrintStream(cout,1,1);
  } else {
    cout<<"could not read 'generator' binning\n";
  }

  // pointers to various nodes in the binning scheme
  const TUnfoldBinning *ttbarreco = detectorBinning->FindNode("ttbarreco");
  const TUnfoldBinning *ttbargen  = generatorBinning->FindNode("ttbargen");

  dim = ttbargen->GetDistributionDimension();
  printf("MatrixUnf::MatrixUnf: Creating %d-dim Objects for variables: %s and %s\n", dim, VariableName1.Data(), VariableName2.Data());

  std::vector<const double*> vec_genBinEdges;
  std::vector<int> vec_numGenBins;
  fNBinsX = 1;
  fNBinsY = 1;
  for (int i = 0; i < ttbargen->GetDistributionDimension(); ++i) {
    const TVectorD* vec_binEdges = ttbargen->GetDistributionBinning(i);
    std::cout << "Adding bin edges for " << i << "th distribution" << std::endl;
    for (int j = 0; j < vec_binEdges->GetNoElements(); ++j) {
      std::cout << "Adding bin edge: " << vec_binEdges->GetMatrixArray()[j] << " ";
    }
    std::cout << std::endl;
    vec_genBinEdges.push_back(vec_binEdges->GetMatrixArray());
    fNBinsX *= vec_binEdges->GetNoElements() - 1;
    vec_numGenBins.push_back(vec_binEdges->GetNoElements() - 1);
  }

  std::vector<const double*> vec_recoBinEdges;
  for (int i = 0; i < ttbarreco->GetDistributionDimension(); ++i) {
    const TVectorD* vec_binEdges = ttbarreco->GetDistributionBinning(i);
    vec_recoBinEdges.push_back(vec_binEdges->GetMatrixArray());
    fNBinsY *= vec_binEdges->GetNoElements() - 1;
  }

  std::cout << "Number of bins on the X (GEN) = " << fNBinsX << ", and the number of bins on the Y (RECO) = " << fNBinsY << std::endl;

  // TODO: Move bin-dependent factors (BinFactorFunction) calculations to unfolding input creation as this is more appropriate
  // generator level binning with BFF
  generatorBinning_withBFF = new TUnfoldBinning("generatorBinning_withBFF");
  TUnfoldBinning* genDistribution_withBFF = generatorBinning_withBFF->AddBinning("ttbargen");
  for (int i = 0; i < ttbargen->GetDistributionDimension(); ++i) {
    genDistribution_withBFF->AddAxis(ttbargen->GetDistributionAxisLabel(i), vec_numGenBins[i], vec_genBinEdges[i], false, false);
  }
  
  // reweight bias distribution for systematic studies
  GenVecHist_forBias = (TH1D*) GenVec->Clone((Form("%s_forBias",GenVec->GetName())));
  if(freweight_bias_test) MatrixUnf::ReweightBiasTest(GenVecHist_forBias);

  //this function multiplies the distribution being unfolded such that any deviation from the SM will appear as a linear effect, so that regularising the curvature gives no bias.
  //e.g. for Ckk, dσ/dξ=-2*(1−Ckk*ξ)*ln|ξ|,  where  ξ=costheta1*costheta2, so we have to multiply by -1/ln|ξ| so that changes in Ckk give a linear effect.
  //if the functional form is unknown, instead this could be the function that makes the distribution ~linear, meaning deviations are more likely to be linear to first order.
  TVectorD *BFFvec = new TVectorD(fNBinsX);
  calculateBFF(GenVecHist_forBias, BFFvec, genDistribution_withBFF);

  //if( MatrixUnf::isBPMtype(fVariableName1) || MatrixUnf::isBtype(fVariableName1) || MatrixUnf::isDtype(fVariableName1) ) genDistribution->SetBinFactorFunction(1.0,BFFfunc); 
  genDistribution_withBFF->SetBinFactorFunction(1.0, (TF1*) BFFvec); 

  if (dbgUnf) 
    generatorBinning_withBFF->PrintStream(cout, 1, 1);

  // Create rebinned TUnfoldBinning objects
  std::vector<double*> vec_genBinEdges_rebinnedBefore;
  std::vector<int> vec_numGenRebinnedBins;
  for (int i = 0; i < ttbargen->GetDistributionDimension(); ++i) {
    int numRebinnedBins = (ttbargen->GetDistributionBinning(i)->GetNoElements() - 1) / frebinfine;
    vec_numGenRebinnedBins.push_back(numRebinnedBins);

    double* rebinnedGenBinEdges = new double[numRebinnedBins + 1];
    std::cout << "Adding bin edges for " << i << "th rebinnedB gen distribution" << std::endl;
    for (int j = 0; j <= numRebinnedBins; ++j) {
      rebinnedGenBinEdges[j] = vec_genBinEdges[i][j * frebinfine];
      std::cout << "Adding bin edge: " << vec_genBinEdges[i][j * frebinfine] << " ";
    }
    vec_genBinEdges_rebinnedBefore.push_back(rebinnedGenBinEdges);
  }

  std::vector<double*> vec_recoBinEdges_rebinnedBefore;
  std::vector<int> vec_numRecoRebinnedBins;
  for (int i = 0; i < ttbarreco->GetDistributionDimension(); ++i) {
    int numRebinnedBins = (ttbarreco->GetDistributionBinning(i)->GetNoElements() - 1) / frebinfine;
    vec_numRecoRebinnedBins.push_back(numRebinnedBins);

    double* rebinnedRecoBinEdges = new double[numRebinnedBins + 1];
    for (int j = 0; j <= numRebinnedBins; ++j) {
      rebinnedRecoBinEdges[j] = vec_recoBinEdges[i][j * frebinfine];
    }
    vec_recoBinEdges_rebinnedBefore.push_back(rebinnedRecoBinEdges);
  }

  // detector level binning
  detectorBinning_rebinnedB = new TUnfoldBinning("detector_rebinnedB");
  TUnfoldBinning* detectorDistribution_rebinnedB = detectorBinning_rebinnedB->AddBinning("ttbarreco_rebinnedB");
  for (int i = 0; i < ttbarreco->GetDistributionDimension(); ++i) {
    detectorDistribution_rebinnedB->AddAxis(ttbarreco->GetDistributionAxisLabel(i), vec_numRecoRebinnedBins[i], vec_recoBinEdges_rebinnedBefore[i], false, false);
  }

  // generator level binning
  generatorBinning_rebinnedB = new TUnfoldBinning("generator_rebinnedB");
  TUnfoldBinning* genDistribution_rebinnedB = generatorBinning_rebinnedB->AddBinning("ttbargen_rebinnedB");
  for (int i = 0; i < ttbargen->GetDistributionDimension(); ++i) {
    genDistribution_rebinnedB->AddAxis(ttbargen->GetDistributionAxisLabel(i), vec_numGenRebinnedBins[i], vec_genBinEdges_rebinnedBefore[i], false, false);
  }

  // set BinFactorFunction
  GenVecHist_forBias_rebinnedB = (TH1D*)GenVecHist_forBias->Clone((Form("%s_rebinnedB", GenVecHist_forBias->GetName())));
  if (frebinfine > 1) {
    MatrixUnf::rebinMultidimensionalInput(GenVecHist_forBias_rebinnedB,
                                        generatorBinning_rebinnedB->FindNode("ttbargen_rebinnedB"),
                                        generatorBinning->FindNode("ttbargen"),
                                        frebinfine);
  }

  TVectorD *BFFvec_rebinnedB = new TVectorD(fNBinsX / pow(frebinfine, dim));
  calculateBFF(GenVecHist_forBias_rebinnedB, BFFvec_rebinnedB, genDistribution_rebinnedB);
  //if( MatrixUnf::isBPMtype(fVariableName1) || MatrixUnf::isBtype(fVariableName1) || MatrixUnf::isDtype(fVariableName1) ) genDistribution_rebinnedB->SetBinFactorFunction(1.0,BFFfunc_rebinnedB); 
  genDistribution_rebinnedB->SetBinFactorFunction(1.0, (TF1*) BFFvec_rebinnedB); 
  if (dbgUnf)
    generatorBinning_rebinnedB->PrintStream(cout, 1, 1);

  // Amandeep : Dump the gen level binning to xml since it is local to the constructor scope
  std::ofstream xmlOut("binning/" + fVariableName1 + "_binning_rebinnedB.xml");
  TUnfoldBinningXML::ExportXML(*detectorBinning_rebinnedB , xmlOut, kTRUE,  kFALSE);
  TUnfoldBinningXML::ExportXML(*generatorBinning_rebinnedB, xmlOut, kFALSE, kTRUE);
  TUnfoldBinningXML::WriteDTD("binning/tunfoldbinning.dtd");
  xmlOut.close();
  // End

  // result histos will be created and filled when calculating the results not before. 
  printf("Done!\n"); 
  return; 
}

//destructor:
MatrixUnf::~MatrixUnf() { 
  
  printf("MatrixUnf::~MatrixUnf: Deleting Object in TObjArrays and the TObjArrays itself only. Nothing else!\n"); 

  fHistUnfArray->Clear(); 
  fResultUnfHists->Clear();
  fBkgHistArray->Clear();
  fResSysHist->Clear();
  fAequiResSysHist->Clear();
  fTUnfDeltaSys->Clear();
  fTUnfDeltaSys_rebinnedA->Clear();
  fTUnfDeltaSys_rebinnedB->Clear();
  fTUnfCovMatSys->Clear();
  fTUnfCovMatSys_rebinnedA->Clear();
  fTUnfCovMatSys_rebinnedB->Clear();
  fTUnfCovMatSysNorm->Clear();
  fTUnfCovMatSysNorm_rebinnedA->Clear();
  fTUnfCovMatSysNorm_rebinnedB->Clear();
  return; 
}

/*********************************************************/

void MatrixUnf::AddObjectsToFile() { 
  
  printf("MatrixUnf::AddObjectsToFile Adding ObjArrays to the Output File\n"); 
  fHistUnfArray->Write(); 
  //fResultUnfHists->Write();
  fBkgHistArray->Write();
  fResSysHist->Write();
  fAequiResSysHist->Write();
  fTUnfDeltaSys->Write();
  fTUnfDeltaSys_rebinnedA->Write();
  fTUnfDeltaSys_rebinnedB->Write();
  fTUnfCovMatSys->Write();
  fTUnfCovMatSys_rebinnedA->Write();
  fTUnfCovMatSys_rebinnedB->Write();
  fTUnfCovMatSysNorm->Write();
  fTUnfCovMatSysNorm_rebinnedA->Write();
  fTUnfCovMatSysNorm_rebinnedB->Write();

}

/*********************************************************

 * Dump -- std debug output of a root class. 
 *
 */
void MatrixUnf::DumpObjArr() { 

  printf("MatrixUnf::DumpObjArr\n"); 
  fHistUnfArray->Dump();  
  fResultUnfHists->Dump();
  fBkgHistArray->Dump();
  fResSysHist->Dump();
  fAequiResSysHist->Dump();
  fTUnfDeltaSys->Dump();
  fTUnfDeltaSys_rebinnedA->Dump();
  fTUnfDeltaSys_rebinnedB->Dump();
  fTUnfCovMatSys->Dump();
  fTUnfCovMatSys_rebinnedA->Dump();
  fTUnfCovMatSys_rebinnedB->Dump();
  fTUnfCovMatSysNorm->Dump();
  fTUnfCovMatSysNorm_rebinnedA->Dump();
  fTUnfCovMatSysNorm_rebinnedB->Dump();

}

/***************************************************
 * create Results -- after filling all the above 
 * histogram it is now possible to create and 
 * fill the results just by dividing!
 */ 
Bool_t MatrixUnf::createInputResults(Bool_t ensTest){    //create and fill all  result histos 
  
  Bool_t b0 = kTRUE;
  Bool_t b1 = kTRUE; 
  Bool_t b2 = kTRUE; 
  Bool_t b3 = kTRUE; 
  Bool_t b4 = kTRUE; 
  Bool_t b5 = kTRUE; 
  Bool_t b6 = kTRUE;
  Bool_t b7 = kTRUE;
  Bool_t b8 = kTRUE;

  b1 = createResMat();
  b3 = createAequiResMat();

  //don't create the supplementary results when running pseudoexperiments
  if(!ensTest) {
    std::cout << "Calculating Eff BbB Vec" << std::endl;
    b0 = createEffBbBVec();
    std::cout << "Create ResMatSys" << std::endl;
    b2 = createResMatSys();
    std::cout << "Create Aequi ResMat Sys" << std::endl;
    b4 = createAequiResMatSys(); //add in .h
    std::cout << "Create Aequi GenVec" << std::endl;
    b5 = createAequiGenVec(); // not needed here, allGen (0) directly set in setInputDist
    b7 = createTruthMigrationPurityStability();
    //b8 = createResolution();
  }

  if(( b0==kFALSE) || (b1==kFALSE) || (b2==kFALSE) || (b3==kFALSE) || (b4==kFALSE) || (b5==kFALSE) || (b6==kFALSE) || (b7==kFALSE) || (b8==kFALSE)) return kFALSE; 
  else return kTRUE; 
}


Bool_t MatrixUnf::createResMat() {
  
  TH2D* ResMat = (TH2D*) GetRecResHist()->Clone();
  ResMat->SetName(fVariableName1+"ResMat");
  TString Title = Form("ResMat for %s", fVariableName1.Data()); 
  ResMat->SetTitle(Title.Data()); 
  ResMat->GetXaxis()->SetTitle(fXAxisTitle);
  
  //if(dbgUnf) printf("MatrixUnf::createResMat: NO subtraction done!\n");
  
  fResultUnfHists->AddAt(ResMat,2); 
  return kTRUE; 
  
}

Bool_t MatrixUnf::createResMatSys() {
  TObjArray *sysarray = (TObjArray*) GetResSysHists()->Clone();
  for(Int_t k=0; k<sysarray->GetEntries(); k++){
    TH2D* ResMatSys = (TH2D*)((TH2D*)sysarray->At(k))->Clone();
    std::string str, str2; 
    str        = ResMatSys->GetName();
    size_t pos = str.find("_norm");
    str2       = str.substr(0,pos);
    //ResMatSys->SetName(Form("%sResMatSys%d",fVariableName1.Data(),k));
    ResMatSys->SetName(str2.c_str());
    //TString Title = Form("ResMatSys%d for %s", k, fVariableName1.Data()); 
    TString Title   = Form("%s for %s", ResMatSys->GetName(), fVariableName1.Data()); 
    ResMatSys->SetTitle(Title.Data()); 
    ResMatSys->GetXaxis()->SetTitle(fXAxisTitle);
  
    //if(dbgUnf) printf("MatrixUnf::createResMatSys: NO subtraction done!\n");
  }
  delete sysarray;
  return kTRUE; 
  
}

// Aequidistant measurement vector
Bool_t MatrixUnf::createAequiGenVec() {
  
  TH1D* GenVec_aequi = new TH1D("", "", GetGenVecHist()->GetNbinsX(), 0, GetGenVecHist()->GetNbinsX());
  GenVec_aequi->SetName(fVariableName1+"AequiGenVec");
  TString Title = Form("AequiGenVec for %s", fVariableName1.Data()); 
  GenVec_aequi->SetTitle(Title.Data()); 
  GenVec_aequi->GetXaxis()->SetTitle(fXAxisTitle);
  
  for(Int_t j =1; j <= GetGenVecHist()->GetNbinsX(); j++) {
    GenVec_aequi->SetBinContent(j,GetGenVecHist()->GetBinContent(j));
    GenVec_aequi->SetBinError(j, GetGenVecHist()->GetBinError(j));
    if(dbgUnf) printf("MatrixUnf::createAequiGenVec (%d): Gen %.2f +- %.2f\n",j, GenVec_aequi->GetBinContent(j), GenVec_aequi->GetBinError(j));
  }
  if(dbgUnf) printf("MatrixUnf::createAequiGenVec: AequiGenVec with %f entries, %d bins for Unfold... \n", GenVec_aequi->Integral(), (int) GenVec_aequi->GetNbinsX());
  
  fResultUnfHists->AddAt(GenVec_aequi,0);
  return kTRUE; 
}


// aequidistant measurement vector
Bool_t MatrixUnf::createAequiMeasVec(TH1D* MeasVec) {
  
  TH1D* MeasVec_aequi = new TH1D("", "", MeasVec->GetNbinsX(), 0, MeasVec->GetNbinsX());
  MeasVec_aequi->SetName(fVariableName1+"AequiMeasVec");
  TString Title = Form("AequiMeasVec for %s", fVariableName1.Data()); 
  MeasVec_aequi->SetTitle(Title.Data()); 
  MeasVec_aequi->GetXaxis()->SetTitle(fXAxisTitle);
  
  for(Int_t j =1; j <= MeasVec->GetNbinsX(); j++) {
    MeasVec_aequi->SetBinContent(j,MeasVec->GetBinContent(j));
    MeasVec_aequi->SetBinError(j, MeasVec->GetBinError(j));
    //if(dbgUnf) printf("MatrixUnf::createAequiMeasVec (%d): Meas %.2f +- %.2f\n",j, MeasVec_aequi->GetBinContent(j), MeasVec_aequi->GetBinError(j));
  }
  if(dbgUnf) printf("MatrixUnf::createAequiMeasVec: AequiMeasVec with %f entries, %d bins for Unfold... \n", MeasVec_aequi->Integral(), (int) MeasVec_aequi->GetNbinsX());
  
  fResultUnfHists->AddAt(MeasVec_aequi,3);
  return kTRUE; 
}

// aequidistant measurement vector
Bool_t MatrixUnf::createAequiMeasVecBgSubtracted(TH1D* MeasVec) {
  
  TH1D* MeasVec_aequi_bgSubtr = new TH1D("", "", MeasVec->GetNbinsX(), 0, MeasVec->GetNbinsX());
  MeasVec_aequi_bgSubtr->SetName(fVariableName1+"AequiMeasVec");
  TString Title = Form("AequiMeasVecBgSubtracted for %s", fVariableName1.Data()); 
  MeasVec_aequi_bgSubtr->SetTitle(Title.Data()); 
  MeasVec_aequi_bgSubtr->GetXaxis()->SetTitle(fXAxisTitle);
  Double_t checkEntries=0;
  
/*
  if(dbgUnf) printf("MatrixUnf::createAequiMeasVecBgSubtracted: Subtract input histo (%d bins) by bg histo (%d bins) with %.2f integral! Adding sidebands: %d to %d bins\n", MeasVec->GetNbinsX(),this->GetFullRecBgVec()->GetNbinsX(),MeasVec->Integral(),MeasVec->GetNbinsX(),MeasVec_aequi_bgSubtr->GetNbinsX());
  
  for(Int_t j =1; j <= MeasVec->GetNbinsX(); j++) {
    if(dbgUnf) printf("MatrixUnf::createAequiMeasVecBgSubtracted (%d): Meas %.2f +- %.2f - Bg %.2f +- %.2f = %.2f +- %.2f\n",j, MeasVec->GetBinContent(j),MeasVec->GetBinError(j), this->GetFullRecBgVec()->GetBinContent(j),this->GetFullRecBgVec()->GetBinError(j),  MeasVec->GetBinContent(j)-this->GetFullRecBgVec()->GetBinContent(j),sqrt(MeasVec->GetBinError(j)*MeasVec->GetBinError(j)+this->GetFullRecBgVec()->GetBinError(j)*this->GetFullRecBgVec()->GetBinError(j)) );
  }
*/

  for (int i = 0; i < this->GetBgHists()->GetEntries(); ++i)
  {
    MeasVec->Add( (TH1D*) this->GetBgHists()->At(i), -1.);
  }
  
  for(Int_t j =1; j <= MeasVec_aequi_bgSubtr->GetNbinsX(); j++) {
    MeasVec_aequi_bgSubtr->SetBinContent(j,MeasVec->GetBinContent(j));
    MeasVec_aequi_bgSubtr->SetBinError(j,MeasVec->GetBinError(j));
    checkEntries += MeasVec_aequi_bgSubtr->GetBinContent(j);
    
    // avoid negative entries
    if((MeasVec_aequi_bgSubtr->GetBinContent(j))<0.) {
      MeasVec_aequi_bgSubtr->SetBinContent(j,0.);
      MeasVec_aequi_bgSubtr->SetBinError(j,1.);
      if(dbgUnf)printf("MatrixUnf::createAequiMeasVecBgSubtracted (%d): Meas %.2f\n",j, MeasVec_aequi_bgSubtr->GetBinContent(j));
    }
    if(dbgUnf) printf("MatrixUnf::createAequiMeasVecBgSubtracted (%d): BgSubtr %.2f +- %.2f\n",j, MeasVec_aequi_bgSubtr->GetBinContent(j), MeasVec_aequi_bgSubtr->GetBinError(j));
  }
  if(dbgUnf) printf("MatrixUnf::createAequiMeasVecBgSubtracted: AequiMeasVecBgSubtracted with %d entries, %d bins for Unfold... \n", (Int_t) checkEntries,MeasVec_aequi_bgSubtr->GetNbinsX());
  
  fResultUnfHists->AddAt(MeasVec_aequi_bgSubtr,30);
  return kTRUE; 
}

// matrix aequidistant
Bool_t MatrixUnf::createAequiResMat() {
  
  if(dbgUnf) printf("MatrixUnf::createAequiResMat \n");
  
  TH2D* ResMat_aequi = new TH2D("", "", GetResMat()->GetNbinsX(), 0, GetResMat()->GetNbinsX(), GetResMat()->GetNbinsY(), 0, GetResMat()->GetNbinsY());
  printf("here \n");
  ResMat_aequi->SetName(fVariableName1+"AequiResMat");
  TString Title = Form("AequiResMat for %s", fVariableName1.Data()); 
  ResMat_aequi->SetTitle(Title.Data()); 
  ResMat_aequi->GetXaxis()->SetTitle(fXAxisTitle);
  ResMat_aequi->GetYaxis()->SetTitle(fYAxisTitle);
  ResMat_aequi->Sumw2();
  
  if(dbgUnf) printf("MatrixUnf::createAequiResMat: Copying bins... \n");
  
  // convert to efficiency corrected number of entries --> eff correction done at a later point
  // remember: does not make a difference as only global factor, migmat is anyway normalized to 1 internally
  
  for(Int_t i =1; i <= ResMat_aequi->GetNbinsX(); i++)  {
    for(Int_t j =0; j <= ResMat_aequi->GetNbinsY() + 1; j++) { //jacob: include the underflow bin, j=0, and also the overflow bin so that both underflow and overflow can be used for non-reconstructed events (as in default TUnfold)
      ResMat_aequi->SetBinContent(i,j,GetResMat()->GetBinContent(i,j));
      ResMat_aequi->SetBinError(i,j,GetResMat()->GetBinError(i,j));
      
      //if(dbgUnf) printf("MatrixUnf::createAequiResMat: %d-%d Entries:%.2e Error:%.2e rel. Error is: %.2f\n",i,j,ResMat_aequi->GetBinContent(i,j), ResMat_aequi->GetBinError(i,j), ResMat_aequi->GetBinError(i,j)/ResMat_aequi->GetBinContent(i,j));

      if( !(ResMat_aequi->GetBinContent(i,j) >= 0.) ) {
        printf("MatrixUnf::createAequiResMat: Warning, bin with negative content: %.3f\n", ResMat_aequi->GetBinContent(i,j) );
        ResMat_aequi->SetBinContent(i,j,0.01*TMath::Abs(GetResMat()->GetBinContent(i,j))); //just to avoid having exactly zero, which TUnfold doesn't seem to like when combined with a non-zero error
        //make sure the errors are conservative
        ResMat_aequi->SetBinError(i,j, sqrt(1.01*1.01*GetResMat()->GetBinContent(i,j)*GetResMat()->GetBinContent(i,j) + GetResMat()->GetBinError(i,j)*GetResMat()->GetBinError(i,j) ) );
      }

      if( !(ResMat_aequi->GetBinError(i,j) >= 0.) ) {
        printf("MatrixUnf::createAequiResMat: Error, bin with negative uncertainty: %.3f\n", ResMat_aequi->GetBinError(i,j) );
        ResMat_aequi->SetBinError(i,j, sqrt(fabs(ResMat_aequi->GetBinContent(i,j))) ); //temporary fix for AMCATNLO sample with a bin with nan error.
        printf("MatrixUnf::createAequiResMat: Set bin %d %d uncertainty to sqrt(content=%.3f): %.3f\n", i, j, ResMat_aequi->GetBinContent(i,j), ResMat_aequi->GetBinError(i,j) );
      }

    }
  }
  
  //check there is nothing in the other underflow and overflow bins. the gen-level underflow and overflow should be empty, and for reco-level we choose to put the non-reconstructed events only in the underflow bin (although putting them in the overflow bin should work too)
  for(Int_t i =0; i <= ResMat_aequi->GetNbinsX()+1; i++) {
    for(Int_t j =0; j <= ResMat_aequi->GetNbinsY()+1; j++) {
      if( (j == ResMat_aequi->GetNbinsY()+1 || i == ResMat_aequi->GetNbinsX()+1 || i == 0) && (GetResMat()->GetBinContent(i,j) != 0 || GetResMat()->GetBinError(i,j) != 0) ) printf("MatrixUnf::createAequiResMat: Warning, unexpected non-zero underflow or overflow bin.\n");
    }
  }

  if(dbgUnf) printf("MatrixUnf::createAequiResMat: Finished! \n");
  fResultUnfHists->AddAt(ResMat_aequi,4);
  return kTRUE; 
}

// matrix aequidistant //ajeeta
Bool_t MatrixUnf::createAequiResMatSys() {
  
  if(dbgUnf) printf("MatrixUnf::createAequiResMatSys \n");

  TObjArray *sysarray = (TObjArray*) GetResSysHists()->Clone();
  for(Int_t k=0; k<sysarray->GetEntries(); k++){

    TH2D* ResMatSys = (TH2D*)((TH2D*) ((*sysarray)[k]))->Clone();
    TH2D* ResMatSys_aequi = new TH2D("", "", ResMatSys->GetNbinsX(), 0, ResMatSys->GetNbinsX(), ResMatSys->GetNbinsY(), 0, ResMatSys->GetNbinsY());
    printf("here \n");

    std::string str, str2; 
    str  = ResMatSys->GetName();
    size_t pos = str.find("_norm");
    str2 = str.substr(0,pos);
    
    std::cout<<str2<<"_AequiResMatSys"<<std::endl;

    ResMatSys_aequi->SetName(Form("%s_AequiResMatSys",str2.c_str()));
    //ResMatSys_aequi->SetName(Form("%sAequiResMatSys%d", fVariableName1.Data(),k));
    //TString Title = Form("AequiResMatSys%d for %s",k, fVariableName1.Data()); 
    TString Title = Form("%s for %s",ResMatSys_aequi->GetName(), fVariableName1.Data()); 
    ResMatSys_aequi->SetTitle(Title.Data()); 
    ResMatSys_aequi->GetXaxis()->SetTitle(fXAxisTitle);
    ResMatSys_aequi->GetYaxis()->SetTitle(fYAxisTitle);
    ResMatSys_aequi->Sumw2();
  
    if(dbgUnf) printf("MatrixUnf::createAequiResMatSys: Copying bins... \n");
  
    // convert to efficiency corrected number of entries --> eff correction done at a later point
    // remember: does not make a difference as only global factor, migmat is anyway normalized to
    // 1 internally
  
    for(Int_t i =1; i <= ResMatSys_aequi->GetNbinsX(); i++)  {
      for(Int_t j =0; j <= ResMatSys_aequi->GetNbinsY() + 1; j++) { //jacob: include the underflow bin, j=0, and also the overflow bin so that both underflow and overflow can be used for non-reconstructed events (as in default TUnfold)
        ResMatSys_aequi->SetBinContent(i,j,ResMatSys->GetBinContent(i,j));
        ResMatSys_aequi->SetBinError(i,j,ResMatSys->GetBinError(i,j));
        //if(dbgUnf) printf("MatrixUnf::createAequiResMatSys: %d-%d Entries:%.2e Error:%.2e rel. Error is: %.2f\n",i,j,ResMatSys_aequi->GetBinContent(i,j), ResMatSys_aequi->GetBinError(i,j), ResMatSys_aequi->GetBinError(i,j)/ResMatSys_aequi->GetBinContent(i,j));

        if( !(ResMatSys_aequi->GetBinContent(i,j) >= 0.) ) {
          printf("MatrixUnf::createAequiResMatSys: Warning, bin with negative content: %.3f\n", ResMatSys_aequi->GetBinContent(i,j) );
          ResMatSys_aequi->SetBinContent(i,j,0.01*TMath::Abs(ResMatSys->GetBinContent(i,j))); //just to avoid having exactly zero, which TUnfold doesn't seem to like when combined with a non-zero error
          //make sure the errors are conservative
          ResMatSys_aequi->SetBinError(i,j, sqrt(1.01*1.01*ResMatSys->GetBinContent(i,j)*ResMatSys->GetBinContent(i,j) + ResMatSys->GetBinError(i,j)*ResMatSys->GetBinError(i,j) ) );
        }
        if( !(ResMatSys_aequi->GetBinError(i,j) >= 0.) ) {
          printf("MatrixUnf::createAequiResMatSys: Error, bin with negative uncertainty: %.3f\n", ResMatSys_aequi->GetBinError(i,j) );
          ResMatSys_aequi->SetBinError(i,j, sqrt(fabs(ResMatSys_aequi->GetBinContent(i,j))) ); //temporary fix for AMCATNLO sample with a bin with nan error.
          printf("MatrixUnf::createAequiResMatSys: Set bin %d %d uncertainty to sqrt(content=%.3f): %.3f\n", i, j, ResMatSys_aequi->GetBinContent(i,j), ResMatSys_aequi->GetBinError(i,j) );
        }

      }
    }
    
    //check there is nothing in the other underflow and overflow bins. the gen-level underflow and overflow should be empty, and for reco-level we choose to put the non-reconstructed events only in the underflow bin (although putting them in the overflow bin should work too)
    for(Int_t i =0; i <= ResMatSys_aequi->GetNbinsX()+1; i++) {
      for(Int_t j =0; j <= ResMatSys_aequi->GetNbinsY()+1; j++) {
        if( (j == ResMatSys_aequi->GetNbinsY()+1 || i == ResMatSys_aequi->GetNbinsX()+1 || i == 0) && (ResMatSys->GetBinContent(i,j) != 0 || ResMatSys->GetBinError(i,j) != 0) ) printf("MatrixUnf::createAequiResMatSys: Warning, unexpected non-zero underflow or overflow bin.\n");
      }
    }
    //ResMatSys_aequi->Write();
    if(dbgUnf) printf("MatrixUnf::createAequiResMatSys: Finished! \n");
    fAequiResSysHist->AddAt(ResMatSys_aequi,k); //put in a diff objarray ajeeta
  }
  return kTRUE; 
}

// Efficiency Histos from gen distributions: 
Bool_t MatrixUnf::createEffBbBVec() { 
  
  if(GetGenVecHist()->GetEntries()==0){ 
    printf("MatrixUnf::createEffBbBVec: GenHist for Variable %s not filled!!\n", fVariableName1.Data()); 
    return kFALSE; 
  }  
  
  //if(dbgUnf) printf("MatrixUnf::createEffBbBVec: Calculating Errors\n");
  
  // Total efficiency
  TH1D *TotEffiHist = NULL;
  // Acceptance
  TH1D *AccepHist = NULL;
  // Efficiency
  TH1D *EffiHist = NULL;  

  //if(dbgUnf) printf("MatrixUnf::createEffBbBVec: Creating the final result #Gen %.1f\n", GetGenVecHist()->GetEntries());
  
  if(GetGenVecHist() !=NULL) {
    TotEffiHist  = (TH1D*) GetGenVecHist()->Clone(fVariableName1+"TotEffiHist");
    TotEffiHist->Reset();
    TString Title = Form("Total-Efficiency in %s Vector (bin-by-bin)", fVariableName1.Data()); 
    TotEffiHist->SetTitle(Title.Data()); 
    TotEffiHist->GetXaxis()->SetTitle(fXAxisTitle);
    TotEffiHist->GetYaxis()->SetTitle("Acceptance #times Efficiency");
    TH1D* dummy = (TH1D*) GetResMat()->ProjectionX("_px",1,fNBinsY);
    
    for(Int_t i=1; i <= dummy->GetNbinsX();i++){
      if(GetGenVecHist()->GetBinContent(i)>0) TotEffiHist->SetBinContent(i, dummy->GetBinContent(i) / GetGenVecHist()->GetBinContent(i));
      if(dbgUnf) printf("MatrixUnf::createEffBbBVec: N(part.) fitted %.2f / Gentot %.0f --> EffBbBVector %.2f\n", dummy->GetBinContent(i), GetGenVecHist()->GetBinContent(i), (dummy->GetBinContent(i)/GetGenVecHist()->GetBinContent(i)));
    }
    
    //TotEffiHist->GetYaxis()->SetRangeUser(0.0,1.2); 
    TotEffiHist->Draw(); 
    
    delete dummy; 
  }

  if(GetVisGenVecHist() !=NULL) {
    EffiHist  = (TH1D*) GetGenVecHist()->Clone(fVariableName1+"EffiHist");
    EffiHist->Reset();
    TString Title = Form("Detector-Efficiency in %s Vector (bin-by-bin)", fVariableName1.Data()); 
    EffiHist->SetTitle(Title.Data()); 
    EffiHist->GetXaxis()->SetTitle(fXAxisTitle);
    EffiHist->GetYaxis()->SetTitle("Efficiency");

    TH1D* dummy1 = (TH1D*) GetResMat()->ProjectionX("_px",1,fNBinsY);
    
    for(Int_t i=1; i <= dummy1->GetNbinsX();i++){
      if(GetVisGenVecHist()->GetBinContent(i)>0) EffiHist->SetBinContent(i, dummy1->GetBinContent(i) / GetVisGenVecHist()->GetBinContent(i));
      if(dbgUnf) printf("MatrixUnf::createEffBbBVec: N(part.) fitted %.2f / Gentot %.0f --> EffBbBVector %.2f\n", dummy1->GetBinContent(i), GetVisGenVecHist()->GetBinContent(i), (dummy1->GetBinContent(i)/GetVisGenVecHist()->GetBinContent(i)));
    }
    
    //EffiHist->GetYaxis()->SetRangeUser(0.0,1.2); 
    EffiHist->Draw(); 
    
    delete dummy1; 
  }

  if(GetGenVecHist() !=NULL) {
    AccepHist  = (TH1D*) GetGenVecHist()->Clone(fVariableName1+"AccepHist");
    AccepHist->Reset();
    TString Title = Form("Acceptance in %s Vector (bin-by-bin)", fVariableName1.Data()); 
    AccepHist->SetTitle(Title.Data()); 
    AccepHist->GetXaxis()->SetTitle(fXAxisTitle);
    AccepHist->GetYaxis()->SetTitle("Acceptance");

    TH1D* dummy2 = (TH1D*) GetVisGenVecHist()->Clone(fVariableName1+"dummy2");
    
    for(Int_t i=1; i <= dummy2->GetNbinsX();i++){
      if(GetGenVecHist()->GetBinContent(i)>0) AccepHist->SetBinContent(i, dummy2->GetBinContent(i) / GetGenVecHist()->GetBinContent(i));
      if(dbgUnf) printf("MatrixUnf::createEffBbBVec: N(part.) fitted %.2f / Gentot %.0f --> EffBbBVector %.2f\n", dummy2->GetBinContent(i), GetGenVecHist()->GetBinContent(i), (dummy2->GetBinContent(i)/GetGenVecHist()->GetBinContent(i)));
    }
    
    // AccepHist->GetYaxis()->SetRangeUser(0.0,1.2); 
    AccepHist->Draw(); 
    
    delete dummy2; 
    
    if(dbgUnf) printf("MatrixUnf::createEffBbBVec: Finished! \n");
  }
  
  // Add to result TObjArray: 
  fResultUnfHists->AddAt(TotEffiHist,19); 
  fResultUnfHists->AddAt(EffiHist,1); 
  fResultUnfHists->AddAt(AccepHist,5); 
  return kTRUE; 
}



Bool_t MatrixUnf::SetInputDists(TH1D* GenVec, TH1D* VisGenVec, TH1D* RecVec, vector<TH1D*> BgVec, TH1D* MissRecVec, TH2D* ResMat, vector<TH2D*>vecResMatSys) {
  
  // put the non-aequidistant bin matrix in 
  fHistUnfArray->AddAt(GenVec, 0);
  fHistUnfArray->AddAt(ResMat, 1); 
  fHistUnfArray->AddAt(VisGenVec, 2);   
  fHistUnfArray->AddAt(RecVec, 3); 
  //fHistUnfArray->AddAt(BgVec, 12);
  fHistUnfArray->AddAt(MissRecVec, 14);

  for(unsigned int k=0; k<BgVec.size(); k++){
   fBkgHistArray->AddAt(BgVec.at(k), k);
  }

  //ajeeta -- add systematics res to a different Object Array 
  for(unsigned int k=0; k<vecResMatSys.size(); k++){
   fResSysHist->AddAt(vecResMatSys[k], k);  
  }
  if(dbgUnf)   printf("MatrixUnf::SetInputDists: Input hists are set!\n");
  //if(dbgUnf) printf("MatrixUnf::SetInputDists: gen=%f, rec=%f, bgRec=%f, missRec=%f, migMat=%f are set!\n", GenVec->Integral(), RecVec->Integral(), BgVec->Integral(), MissRecVec->Integral(), ResMat->Integral());
  if(dbgUnf)   printf("MatrixUnf::SetInputDists: gen=%f, rec=%f, missRec=%f, migMat=%f are set!\n", GenVec->Integral(), RecVec->Integral(), MissRecVec->Integral(), ResMat->Integral());

  return kTRUE;
}

/*
// Efficiency Histo from Matrix
Bool_t MatrixUnf::createEffFromMat() { 
  
  if(dbgUnf) printf("MatrixUnf::createEffFromMat: Calculation Eff\n"); 
  
  TString name = GetResMat()->GetName();
  
  // calculate the efficiencies from the matrix
  TH1D* GenCompHist = (TH1D*) this->GetAequiGen()->Clone(fVariableName1+"GenCompHist_eff");
  TH1D* RecCompHist = new TH1D(fVariableName1+"RecCompHist_eff", "", GetResMat()->GetNbinsY(), 0, GetResMat()->GetNbinsY());
  TH2D* dummyresmat = (TH2D*) GetResMat()->Clone("dummyresmat");

  GenCompHist->Sumw2();
  RecCompHist->Sumw2();

  TString Title = Form("Sum of generated columns in %s vs %s", fVariableName1.Data(), fVariableName2.Data()); 
  GenCompHist->SetTitle(Title.Data()); 
  GenCompHist->GetXaxis()->SetTitle(fXAxisTitle);

  if(frebinfine){
    GenCompHist->Rebin(frebinfine);
    RecCompHist->Rebin(frebinfine);
    dummyresmat->Rebin2D(frebinfine,frebinfine);
  }
  
  //for (Int_t m = 1; m <= GetResMat()->GetNbinsY(); m++){
  for (Int_t m = 1; m <= dummyresmat->GetNbinsY(); m++){
    Double_t rowSumN = 0.;
    //for (Int_t n=1; n <= GetResMat()->GetNbinsX(); n++){
    for (Int_t n=1; n <= dummyresmat->GetNbinsX(); n++){
      //rowSumN += GetResMat()->GetBinContent(n,m);
      rowSumN += dummyresmat->GetBinContent(n,m);
      if(dbgUnf) printf("rowSumEff(rec) %f for %d-%d\n",rowSumN, n, m);
    }
    RecCompHist->SetBinContent(m, rowSumN); 
    if(dbgUnf) printf("rowSumEff(rec) %f for %d\n",rowSumN, m);
  }
  
  TH1D* RecCompHist_save = (TH1D*) RecCompHist->Clone("RecCompHist_save");
  RecCompHist_save->Rebin(2);

  //total efficiency
  TH1D* RecEffFromMatrix = new TH1D(fVariableName1+"EffFromMatHist", "", GenCompHist->GetNbinsX(), 0, GenCompHist->GetNbinsX());
  TString Title1 = Form("Detector-Efficiency in %s from Matrix", fVariableName1.Data()); 
  RecEffFromMatrix->SetTitle(Title1.Data()); 
  RecEffFromMatrix->GetXaxis()->SetTitle(fYAxisTitle);
  RecEffFromMatrix->GetYaxis()->SetTitle("Total Efficiency");

  std::cout<<"++++++++++++++++++++++++++HERE+++++++++++++++++++++++++"<<std::endl;

  for (Int_t i = 1; i <= GenCompHist->GetNbinsX(); i++){
    if (GenCompHist->GetBinContent(i)>0) RecEffFromMatrix->SetBinContent(i, RecCompHist_save->GetBinContent(i) / GenCompHist->GetBinContent(i) );
    if(dbgUnf) printf("Rebin %d: Rec %.1f -- Gen %.1f --> Eff %.2f\n",i , RecCompHist_save->GetBinContent(i), GenCompHist->GetBinContent(i), RecEffFromMatrix->GetBinContent(i));
  }
  
  if(dbgUnf) printf("MatrixUnf::createEffFromMat: Creating the final result\n");
  
  RecEffFromMatrix->GetYaxis()->SetRangeUser(0.0,1.2); 
  RecEffFromMatrix->Draw(); 
  
  if(dbgUnf) printf("MatrixUnf::createEffFromMat: Finished! \n");

  //Add to result TObjArray:                                                                                    
  fResultUnfHists->AddAt(RecEffFromMatrix,6); 

  delete RecCompHist_save;
  delete dummyresmat;
  return kTRUE; 
}
*/

Bool_t MatrixUnf::createResolution() { 
  
  if(dbgUnf) printf("MatrixUnf::createResolution: Calculation of resolution using response Matrix\n"); 
  
  //TString name = GetResMat()->GetName();
  
  //TString Title = Form("(Rec-Gen)/Gen in %s vs %s", fVariableName1.Data(), fVariableName2.Data()); 
  
  TH2D* dummy = (TH2D*) GetResMat()->Clone("dummy");
  
  dummy->FitSlicesY(0,0,-1,0,"Q");

  printf("Get fitted slices: %s\n", Form("%s_0", dummy->GetName()));
  TH1D* fittedSlice = (TH1D*) gDirectory->Get(Form("%s_0;", dummy->GetName()));
  fittedSlice->GetYaxis()->SetTitle(Form("constant %s",fVariableName2.Data()));
  fittedSlice->GetXaxis()->SetTitle(fVariableName1.Data());
  fittedSlice->SetName(Form("%s_ConstGausFitRecVsGen", fVariableName1.Data()));
  
  printf("Get fitted slices: %s\n", Form("%s_1", dummy->GetName()));
  TH1D* fittedSlice1 = (TH1D*) gDirectory->Get(Form("%s_1;", dummy->GetName()));
  Double_t uppBoundY = fittedSlice1->GetXaxis()->GetBinLowEdge(fittedSlice1->GetNbinsX())+fittedSlice1->GetXaxis()->GetBinWidth(fittedSlice1->GetNbinsX());
  Double_t lowBoundY = fittedSlice1->GetXaxis()->GetBinLowEdge(1);
  fittedSlice1->GetYaxis()->SetRangeUser(lowBoundY, uppBoundY);
  fittedSlice1->GetYaxis()->SetTitle(Form("mean %s",fVariableName2.Data()));
  fittedSlice1->GetXaxis()->SetTitle(fVariableName1.Data());
  fittedSlice1->SetName(Form("%s_MeanGausFitRecVsGen", fVariableName1.Data()));

  printf("Get fitted slices: %s\n", Form("%s_2", dummy->GetName()));
  TH1D* fittedSlice2 = (TH1D*) gDirectory->Get(Form("%s_2;", dummy->GetName()));
  fittedSlice2->GetYaxis()->SetTitle(Form("StdDev %s",fVariableName2.Data()));
  fittedSlice2->GetXaxis()->SetTitle(fVariableName1.Data());
  fittedSlice2->SetName(Form("%s_SigmaGausFitRecVsGen", fVariableName1.Data()));

  printf("Get fitted slices: %s\n", Form("%s_chi2", dummy->GetName()));
  TH1D* fittedSlice3 = (TH1D*) gDirectory->Get(Form("%s_chi2;", dummy->GetName()));
  fittedSlice3->GetYaxis()->SetTitle(Form("chi2 %s",fVariableName2.Data()));
  fittedSlice3->GetXaxis()->SetTitle(fVariableName1.Data());
  fittedSlice3->SetName(Form("%s_chi2GausFitRecVsGen", fVariableName1.Data()));

  //Add to result TObjArray: 
  //Parameter 0 (constant) of a gauss fit of the 2D resolution projection along Y
  fResultUnfHists->AddAt(fittedSlice1,31);
  //Parameter 1 (mean) of a gauss fit of the 2D resolution projection along Y
  fResultUnfHists->AddAt(fittedSlice2,32);
  //Parameter 2 (stdev) of a gauss fit of the 2D resolution projection along Y
  fResultUnfHists->AddAt(fittedSlice3,33);
  delete dummy;

  return kTRUE; 


}
Bool_t MatrixUnf::createTruthMigrationPurityStability() { 
  
  if(dbgUnf) printf("MatrixUnf::createTruthMigrationPurityStability: Calculation of Migrations inside Matrix\n"); 

  TString name = GetResMat()->GetName();
  
  TH2D* ResMat_Rebin = (TH2D*) GetResMat()->Clone("ResMat_Rebin");

  TH1D* GenVecHist_Rebin = (TH1D*) GetGenVecHist()->Clone("GenVecHist_Rebin");

  //we want to measure the purity and stability in the final (wide) bins
    // ResMat_Rebin->Rebin2D(frebinfine, frebinfine);
  // ResMat_Rebin->RebinY(2);
  MatrixUnf::rebinMultidimMigrationMatrix(ResMat_Rebin, generatorBinning_rebinnedB->FindNode("ttbargen_rebinnedB"),
                                          detectorBinning_rebinnedB->FindNode("ttbarreco_rebinnedB"), generatorBinning->FindNode("ttbargen"),
                                          detectorBinning->FindNode("ttbarreco"), frebinfine);
  MatrixUnf::rebinMultidimMigrationMatrix(ResMat_Rebin, generatorBinning_rebinnedB->FindNode("ttbargen_rebinnedB"),
                                          detectorBinning_rebinnedB->FindNode("ttbarreco_rebinnedB"), 2, 1);
  
  // GenVecHist_Rebin->Rebin(frebinfine);
  MatrixUnf::rebinMultidimensionalInput(GenVecHist_Rebin, generatorBinning_rebinnedB->FindNode("ttbargen_rebinnedB"),
                                        generatorBinning->FindNode("ttbargen"), frebinfine);

  // calculate the efficiencies fromt the matrix -- same no of bins for both Gen and Rec
  TH1D* GenPurCompHist = new TH1D(fVariableName1+"GenPurCompHist", "", ResMat_Rebin->GetNbinsX(), 0, ResMat_Rebin->GetNbinsX());
  TH1D* RecPurCompHist = new TH1D(fVariableName1+"RecPurCompHist", "", ResMat_Rebin->GetNbinsX(), 0, ResMat_Rebin->GetNbinsX());
  TH1D* RecAndGenCompHist = new TH1D(fVariableName1+"RecAndGenCompHist", "", ResMat_Rebin->GetNbinsX(), 0, ResMat_Rebin->GetNbinsX());
  GenPurCompHist->Sumw2();
  RecPurCompHist->Sumw2();
  RecAndGenCompHist->Sumw2();

  TString Title = Form("Sum of generated columns in %s vs %s", fVariableName2.Data(), fVariableName1.Data()); 
  GenPurCompHist->SetTitle(Title.Data()); 
  GenPurCompHist->GetXaxis()->SetTitle(fXAxisTitle);

  TString Title1 = Form("Sum of reconstructed rows in %s vs %s", fVariableName2.Data(), fVariableName1.Data()); 
  RecPurCompHist->SetTitle(Title1.Data()); 
  RecPurCompHist->GetXaxis()->SetTitle(fXAxisTitle);



  // for the migrations only the diagonal element is important, but have to rebin first
  for (Int_t m = 1; m <= ResMat_Rebin->GetNbinsY(); m++){
    RecAndGenCompHist->SetBinContent(m, ResMat_Rebin->GetBinContent(m,m));
    //if(dbgUnf) printf("rowSumEff(rec) %f for %d-%d\n",ResMat_Rebin->GetBinContent(m,m), m, m);
  }
  

  //for each gen bin sum all the events reconstructed (excluding non-reconstructed events which are in the underflow bin)
  for (Int_t m = 1; m <= ResMat_Rebin->GetNbinsX(); m++){
    Double_t colSumN = 0.;
    for (Int_t n=1; n <= ResMat_Rebin->GetNbinsY();n++){
      colSumN += ResMat_Rebin->GetBinContent(m,n);
      //if(dbgUnf) printf("colSumEff(gen) %f for %d-%d\n",colSumN ,m,n);
    }
    GenPurCompHist->SetBinContent(m, colSumN);
    //if(dbgUnf) printf("colSumEff(gen) %f for %d\n",colSumN ,m);
  }

  //for each reco bin sum all the events generated
  for (Int_t n = 1; n <= ResMat_Rebin->GetNbinsY(); n++){
    Double_t rowSumM = 0.;
    for (Int_t m=1; m <= ResMat_Rebin->GetNbinsX();m++){
      rowSumM += ResMat_Rebin->GetBinContent(m,n);
      //if(dbgUnf) printf("rowSumEff(rec) %f for %d-%d\n",rowSumM ,m,n);
    }
    RecPurCompHist->SetBinContent(n, rowSumM);
    //if(dbgUnf) printf("rowSumEff(rec) %f for %d\n",rowSumM ,n);
  }
  
  TH1D* MigPurFromMatrix = (TH1D*) GenVecHist_Rebin->Clone(fVariableName1+"MigPurFromMatHist");
  MigPurFromMatrix->Reset();
  TString Title2 = Form("Migrations (Purity) in %s from Matrix", fVariableName1.Data()); 
  MigPurFromMatrix->SetTitle(Title2.Data()); 
  MigPurFromMatrix->GetXaxis()->SetTitle(fYAxisTitle);
  MigPurFromMatrix->GetYaxis()->SetTitle("Migrations");
  
  MigPurFromMatrix->SetBinContent(1,0);
  for (Int_t i = 1; i <= RecAndGenCompHist->GetNbinsX(); i++){
    if (RecPurCompHist->GetBinContent(i)>0) MigPurFromMatrix->SetBinContent(i, RecAndGenCompHist->GetBinContent(i) / RecPurCompHist->GetBinContent(i) );
    if(dbgUnf) printf("MatrixUnf::createTruthMigrationPurityStability: Bin %d: 'Purity' %.2f\n",i, MigPurFromMatrix->GetBinContent(i));
  }
  if(dbgUnf) printf("MatrixUnf::createTruthMigrationPurityStability: 'MeanPurity:' %.3f , 'Weighted:' %.3f \n", MigPurFromMatrix->Integral()/MigPurFromMatrix->GetNbinsX() ,  RecAndGenCompHist->Integral() / RecPurCompHist->Integral()  );
  
  //MigPurFromMatrix->GetYaxis()->SetRangeUser(0.0,1.2); 
  MigPurFromMatrix->Draw(); 
  
  if(dbgUnf) printf("MatrixUnf::createTruthMigrationPurityStability: Finished! \n");


  /******************************************/
  //stability calculation  --- 
  TH1D* StabFromMatrix = (TH1D*) GenVecHist_Rebin->Clone(fVariableName1+"StabFromMatHist");
  StabFromMatrix->Reset();
  TString Title3 = Form("Stability in %s from Matrix", fVariableName1.Data()); 
  StabFromMatrix->SetTitle(Title3.Data()); 
  StabFromMatrix->GetXaxis()->SetTitle(fYAxisTitle);
  StabFromMatrix->GetYaxis()->SetTitle("Stability");

  for (Int_t i = 1; i <= RecAndGenCompHist->GetNbinsX(); i++){
    if (GenPurCompHist->GetBinContent(i)>0) StabFromMatrix->SetBinContent(i, RecAndGenCompHist->GetBinContent(i) / GenPurCompHist->GetBinContent(i) );
    if(dbgUnf) printf("MatrixUnf::createTruthMigrationPurityStability: Bin %d: 'Stability' %.2f\n",i, StabFromMatrix->GetBinContent(i));
  }
  if(dbgUnf) printf("MatrixUnf::createTruthMigrationPurityStability: 'MeanStability:' %.3f , 'Weighted:' %.3f \n", StabFromMatrix->Integral()/StabFromMatrix->GetNbinsX() ,  RecAndGenCompHist->Integral() / GenPurCompHist->Integral()  );
  
  if(dbgUnf) printf("MatrixUnf::createTruthMigrationPurityStability: Creating the final result  sumOfWeights %.2e entries %.2e\n",GetResMat()->GetSumOfWeights(), GetResMat()->GetEntries());
  
  //StabFromMatrix->GetYaxis()->SetRangeUser(0.0,1.2); 
  StabFromMatrix->Draw(); 
  
  if(dbgUnf) printf("MatrixUnf::createTruthMigrationPurityStability: Finished! \n");

  //Add to result TObjArray: 
  fResultUnfHists->AddAt(MigPurFromMatrix,28);
  fResultUnfHists->AddAt(StabFromMatrix,29);
  delete ResMat_Rebin;
  return kTRUE; 
}


// ************
// Helper stuff
// ************

TH1D* MatrixUnf::calcDiffHisto(TH1D* in1, TH1D* in2) { 
  
  Float_t Diff1 =0;
  TH1D* out = (TH1D*) in1->Clone(fVariableName1+"EffDiffHisto");
  for (Int_t i=1;i<in1->GetNbinsX();i++){
    Diff1 = TMath::Abs(in2->GetBinContent(i) - in1->GetBinContent(i));
    if(in2->GetBinContent(i) > 0)out->SetBinContent(i, Diff1/(in2->GetBinContent(i)));
  }
  
  return out; 
  
}


Bool_t MatrixUnf::prepRunUnfolding(TH1D* measInput, Int_t RegMode, Int_t skipOuterBins, Double_t tauMin, Double_t tauMax, Double_t Lumi, Bool_t subtrBg, Int_t  useAreaConstrain, Bool_t ensTest, Bool_t minrhoavg, std::vector<TH1D*>&pseudoUnfResults, std::vector<TH1D*>&pseudoUnfResultsRebinnedA, std::vector<TH1D*>&pseudoUnfResultsRebinnedB, std::vector<TH1D*>&pseudoUnfAfbResultsRebinnedA, std::vector<TH1D*>&pseudoUnfAfbResultsRebinnedB, std::vector<TH1D*>&pseudoUnfOtherResultsRebinnedA, std::vector<TH1D*>&pseudoUnfOtherResultsRebinnedB) {

  // Do some Chi2 tests to compare the different results
  bool doChi2tests = kTRUE;

  // Copying a variable bin into a equidistance bins

  if(dbgUnf) printf("\nMatrixUnf::prepRunUnfolding: Preparing Unfolding...subtractBgMCFromInput: %d\n\n", subtrBg);
  
  
  // Create the aequidistant measurement vector
  if(subtrBg)this->createAequiMeasVecBgSubtracted( (TH1D*) measInput->Clone(Form("%s_bkgsub",measInput->GetName())) );
  this->createAequiMeasVec(measInput);
  
  if(dbgUnf) printf("\nMatrixUnf::prepRunUnfolding: Input vector created...\n");
  
  if(!subtrBg)if(this->GetAequiMeas()->GetEntries()<10){
      if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: ERROR too few entries (%f) to unfold: %s\n", this->GetAequiMeas()->GetEntries(), measInput->GetName());
      return 0;
    }
  if(subtrBg)if(this->GetAequiMeasBgSubtracted()->GetEntries()<10){
      if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: ERROR too few entries (%f) to unfold: %s\n", this->GetAequiMeasBgSubtracted()->GetEntries(), measInput->GetName());
      return 0;
  }
  
  if(dbgUnf && !subtrBg) printf("MatrixUnf::prepRunUnfolding: Response Matrix with %d gen and %d rec columns -- measured %d bins\n", this->GetAequiResMat()->GetNbinsX(), this->GetAequiResMat()->GetNbinsY(), this->GetAequiMeas()->GetNbinsX());
  if(dbgUnf &&  subtrBg) printf("MatrixUnf::prepRunUnfolding: Response Matrix with %d gen and %d rec columns -- measured %d bins\n", this->GetAequiResMat()->GetNbinsX(), this->GetAequiResMat()->GetNbinsY(), this->GetAequiMeasBgSubtracted()->GetNbinsX());


  
  // Declare REG and initialize to no regularisation
  TUnfold::ERegMode REG = TUnfold::kRegModeNone;
  switch (RegMode) {
  case 0:  //
    REG = TUnfold::kRegModeNone;
    break; 
  case 1:  //
    REG = TUnfold::kRegModeSize;
    break; 
  case 2:  //
    REG = TUnfold::kRegModeDerivative;
    break; 
  case 3:  //
    REG = TUnfold::kRegModeCurvature;
    break; 
  case 4:  //
    REG = TUnfold::kRegModeMixed;
    break; 
  }
  
  // declare AREA and initialize to no area constrain
  TUnfold::EConstraint AREA = TUnfold::kEConstraintNone;
  switch (useAreaConstrain) {
  case 0:  //
    AREA = TUnfold::kEConstraintNone;
    break; 
  case 1:  //
    AREA = TUnfold::kEConstraintArea;
    break; 
  }

  // TUnfoldDensity::EDensityMode  DENSITY = TUnfoldDensity::kDensityModeNone;
  // TUnfoldDensity::EDensityMode  DENSITY = TUnfoldDensity::kDensityModeBinWidth;
  // TUnfoldDensity::EDensityMode  DENSITY = TUnfoldDensity::kDensityModeUser;
  TUnfoldDensity::EDensityMode  DENSITY    = TUnfoldDensity::kDensityModeBinWidthAndUser;

  // older versions of root <6.10 don't know about the TUnfoldDensity::kDensityModes as used above
  // so use the line below instead ( kDensityModeNone =0, kDensityModeBinWidth =1, kDensityModeUser =2, kDensityModeBinWidthAndUser =3 )
  // TUnfoldDensity::EDensityMode  DENSITY = (TUnfoldDensity::EDensityMode) 3;

  TUnfoldDensity::EScanTauMode TAUMODE    = TUnfoldDensity::kEScanTauRhoAvg;
  // TUnfoldDensity::EScanTauMode TAUMODE = TUnfoldDensity::kEScanTauRhoAvgSys;
 
  // Create a copy of the response matrix to use for the unfolding
  // this can either be set to aequiresmat or regular resmat
  TH2D* ResMat = (TH2D*) this->GetAequiResMat()->Clone();

  // Set up the regularisation and the TUnfoldDensity object
  // if manually regularising always initialise to TUnfold::kRegModeNone 
  // because the regularisation is set via RegularizeBins or RegularizeDistribution 
  // (which add to, rather than override, whatever is set here)
  // TUnfoldDensity unfold(ResMat,TUnfold::kHistMapOutputHoriz, 
  // TUnfold::kRegModeNone, AREA, DENSITY,generatorBinning_withBFF,detectorBinning_withBFF);
  // otherwise it can all be set in this line:

  // Amandeep : Changing here so that it only uses the inner axis for regularization
  // Original :
  // TUnfoldDensity unfold(ResMat,TUnfold::kHistMapOutputHoriz, REG, AREA, DENSITY, generatorBinning_withBFF, detectorBinning, 0, "*[UOB]" );
  
  // Modified : 
  // Default case for 1D, refers to regularizing along all axes
  TString AxisSteering = "*[UOB]";

  // For 2D only regularize along the internal axis
  if(generatorBinning->FindNode("ttbargen")->GetDistributionDimension() == 2) {
    // FIXME : Add support for other types later
    // if      (MatrixUnf::isBtype(fVariableName1))   AxisSteering = fVariableName1(0,3) + "[UOB];mttbar[N]"; // b1k, b1r etc.
    // else if (MatrixUnf::isCtype(fVariableName1))   AxisSteering = fVariableName1(0,4) + "[UOB];mttbar[N]"; // c_kk, c_kr etc.
    // else if (MatrixUnf::isCPMtype(fVariableName1)) AxisSteering = fVariableName1(0,5) + "[UOB];mttbar[N]"; // c_Pnk, c_Mrk etc. 
    // FIXED : Pull from AJ
    AxisSteering = generatorBinning->FindNode("ttbargen")->GetDistributionAxisLabel(0) + "[UOB];" + generatorBinning->FindNode("ttbargen")->GetDistributionAxisLabel(1) + "[N]";
  }

  // FIXME : Generalize this maybe ttbarreco->GetAxisLabel(0)[UOB]; ttbarreco->GetAxisLabel(1)[N] ??
  TUnfoldDensity unfold(ResMat,TUnfold::kHistMapOutputHoriz, REG, AREA, DENSITY, generatorBinning_withBFF, detectorBinning, 0, AxisSteering);
  // End

  // Set regulatisation on a 1-dimensional curve
  // unfold.RegularizeBins(1,1,ResMat->GetNbinsX(),TUnfold::kRegModeCurvature);
  // unfold.RegularizeDistribution(REG, DENSITY, 0, "*[UOB]" ); //this gives the same result as RegularizeBins (except for redefining tau = 2*tau )

  // set the bias distribution, minimized is x-x0 where x0 is the bias-dist. 
  // TUnfold automatically uses the generated xsec for the x0 dist.

  if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: SetBiasDist: %d\n", fsetBiasDisTUnfold);
  
  Double_t scaleBias_temp   = 0.0;
  Double_t oneOverZeroError = 0.0;
  
  if (fsetBiasDisTUnfold && subtrBg)       scaleBias_temp = (this->GetAequiMeasBgSubtracted()->Integral()) / this->GetFullRecVec()->Integral() ;
  else if (fsetBiasDisTUnfold && !subtrBg) scaleBias_temp = (this->GetAequiMeas()->Integral()) / this->GetFullRecVec()->Integral() ;
  else scaleBias_temp = 0.0;
  
  Double_t const scaleBias = scaleBias_temp;
  
  // Set the input distribution
  unfold.SetInput(this->GetAequiMeas(), scaleBias, oneOverZeroError);

  // background subtraction using TUnfold

  std::vector<TString>backgroundnames;
  backgroundnames.push_back("DYeemm");
  backgroundnames.push_back("DYtautau");
  backgroundnames.push_back("singletop");
  backgroundnames.push_back("ww");
  backgroundnames.push_back("wz");
  backgroundnames.push_back("zz");
  backgroundnames.push_back("ttbarW");
  backgroundnames.push_back("ttbarZ");
  backgroundnames.push_back("ttbarbg");
  backgroundnames.push_back("wtolnu");
  backgroundnames.push_back("other");
  // backgroundnames.push_back("ttbarbgviatau");

  std::vector<float>backgrounduncertainties;
  backgrounduncertainties.push_back(0.3);
  backgrounduncertainties.push_back(0.3);
  backgrounduncertainties.push_back(0.3);
  backgrounduncertainties.push_back(0.3);
  backgrounduncertainties.push_back(0.3);
  backgrounduncertainties.push_back(0.3);
  backgrounduncertainties.push_back(0.3);
  backgrounduncertainties.push_back(0.3);
  backgrounduncertainties.push_back(0.3);
  backgrounduncertainties.push_back(0.3);
  backgrounduncertainties.push_back(0.3);
  // backgrounduncertainties.push_back(0.025);  // PDG values: Wtotau: 11.38 +/- 0.21, 1.85%;  Wtoe,mu: 10.67 +/- 0.16, 1.50% [but correlated with Wtotau?];  tautoe/mu: 17.6 +/- 0.04, 0.23%   sum in quad => 2.39%, round up to 2.5%.
  
  const int nCombinedBkgHists = backgroundnames.size();
  TH1D* CombinedBkgHists[nCombinedBkgHists];
  TH1D* CombinedBkgHists_rebinnedB[nCombinedBkgHists];

  if(subtrBg) {

    for (int i_C = 0; i_C < nCombinedBkgHists; ++i_C)
    {
      CombinedBkgHists[i_C] = (TH1D*) this->GetFullRecVec()->Clone(fVariableName1+"_"+backgroundnames.at(i_C)+"bkg");
      CombinedBkgHists[i_C]->Reset();
    }

    for (int i_bkg = 0; i_bkg < this->GetBgHists()->GetEntries(); ++i_bkg)
    {
      TString BkgName = ((TH1D*)this->GetBgHists()->At(i_bkg))->GetName();

      double bkgintegral, bkgerror, bkgnentries, bkgneffentries;
      bkgintegral    = ((TH1D*)this->GetBgHists()->At(i_bkg))->IntegralAndError(1, ((TH1D*)this->GetBgHists()->At(i_bkg))->GetNbinsX(), bkgerror);
      bkgnentries    = ((TH1D*)this->GetBgHists()->At(i_bkg))->GetEntries();
      bkgneffentries = ((TH1D*)this->GetBgHists()->At(i_bkg))->GetEffectiveEntries();

      if(BkgName.Contains("dyee") || BkgName.Contains("dymumu")) {
        if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Subtracting background %s. Type: %s, Uncertainty: %f, Integral: %f +/- %f, nEntries: %f, nEffEntries: %f\n", BkgName.Data(), backgroundnames[0].Data(), backgrounduncertainties[0], bkgintegral, bkgerror, bkgnentries, bkgneffentries);
        CombinedBkgHists[0]->Add( (TH1D*)this->GetBgHists()->At(i_bkg) );
      }

      else if(BkgName.Contains("dytautau")) {
        if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Subtracting background %s. Type: %s, Uncertainty: %f, Integral: %f +/- %f, nEntries: %f, nEffEntries: %f\n", BkgName.Data(), backgroundnames[1].Data(), backgrounduncertainties[1], bkgintegral, bkgerror, bkgnentries, bkgneffentries);
        CombinedBkgHists[1]->Add( (TH1D*)this->GetBgHists()->At(i_bkg) );
      }

      //      else if(BkgName.Contains("ttbarbgviatau")) {
      //        if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Subtracting background %s. Type: %s, Uncertainty: %f, Integral: %f +/- %f, nEntries: %f, nEffEntries: %f\n", BkgName.Data(), backgroundnames[2].Data(), backgrounduncertainties[2], bkgintegral, bkgerror, bkgnentries, bkgneffentries);
      //        CombinedBkgHists[2]->Add( (TH1D*)this->GetBgHists()->At(i_bkg) );
      //      }

      else if(BkgName.Contains("single")) {
        if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Subtracting background %s. Type: %s, Uncertainty: %f, Integral: %f +/- %f, nEntries: %f, nEffEntries: %f\n", BkgName.Data(), backgroundnames[2].Data(), backgrounduncertainties[2], bkgintegral, bkgerror, bkgnentries, bkgneffentries);
        CombinedBkgHists[2]->Add( (TH1D*)this->GetBgHists()->At(i_bkg) );
      }
      else if(BkgName.Contains("ww")) {
        if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Subtracting background %s. Type: %s, Uncertainty: %f, Integral: %f +/- %f, nEntries: %f, nEffEntries: %f\n", BkgName.Data(), backgroundnames[3].Data(), backgrounduncertainties[3], bkgintegral, bkgerror, bkgnentries, bkgneffentries);
        CombinedBkgHists[3]->Add( (TH1D*)this->GetBgHists()->At(i_bkg) );
      }
      else if(BkgName.Contains("wz")) {
        if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Subtracting background %s. Type: %s, Uncertainty: %f, Integral: %f +/- %f, nEntries: %f, nEffEntries: %f\n", BkgName.Data(), backgroundnames[4].Data(), backgrounduncertainties[4], bkgintegral, bkgerror, bkgnentries, bkgneffentries);
        CombinedBkgHists[4]->Add( (TH1D*)this->GetBgHists()->At(i_bkg) );
      }
      else if(BkgName.Contains("zz")) {
        if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Subtracting background %s. Type: %s, Uncertainty: %f, Integral: %f +/- %f, nEntries: %f, nEffEntries: %f\n", BkgName.Data(), backgroundnames[5].Data(), backgrounduncertainties[5], bkgintegral, bkgerror, bkgnentries, bkgneffentries);
        CombinedBkgHists[5]->Add( (TH1D*)this->GetBgHists()->At(i_bkg) );
      }
      else if(BkgName.Contains("ttbarW")) {
        if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Subtracting background %s. Type: %s, Uncertainty: %f, Integral: %f +/- %f, nEntries: %f, nEffEntries: %f\n", BkgName.Data(), backgroundnames[6].Data(), backgrounduncertainties[6], bkgintegral, bkgerror, bkgnentries, bkgneffentries);
        CombinedBkgHists[6]->Add( (TH1D*)this->GetBgHists()->At(i_bkg) );
      }
      else if(BkgName.Contains("ttbarZ")) {
        if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Subtracting background %s. Type: %s, Uncertainty: %f, Integral: %f +/- %f, nEntries: %f, nEffEntries: %f\n", BkgName.Data(), backgroundnames[7].Data(), backgrounduncertainties[7], bkgintegral, bkgerror, bkgnentries, bkgneffentries);
        CombinedBkgHists[7]->Add( (TH1D*)this->GetBgHists()->At(i_bkg) );
      }
      else if(BkgName.Contains("ttbarbg")) {
        if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Subtracting background %s. Type: %s, Uncertainty: %f, Integral: %f +/- %f, nEntries: %f, nEffEntries: %f\n", BkgName.Data(), backgroundnames[8].Data(), backgrounduncertainties[8], bkgintegral, bkgerror, bkgnentries, bkgneffentries);
        CombinedBkgHists[8]->Add( (TH1D*)this->GetBgHists()->At(i_bkg) );
      }
      else if(BkgName.Contains("wtolnu")) {
        if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Subtracting background %s. Type: %s, Uncertainty: %f, Integral: %f +/- %f, nEntries: %f, nEffEntries: %f\n", BkgName.Data(), backgroundnames[9].Data(), backgrounduncertainties[9], bkgintegral, bkgerror, bkgnentries, bkgneffentries);
        CombinedBkgHists[9]->Add( (TH1D*)this->GetBgHists()->At(i_bkg) );
      }

      else {
        if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Subtracting background %s. Type: %s, Uncertainty: %f, Integral: %f +/- %f, nEntries: %f, nEffEntries: %f\n", BkgName.Data(), backgroundnames[10].Data(), backgrounduncertainties[10], bkgintegral, bkgerror, bkgnentries, bkgneffentries);
        CombinedBkgHists[10]->Add( (TH1D*)this->GetBgHists()->At(i_bkg) );
      }
    }

    // Note that we are using a slightly modified version of DoBackgroundSubtraction in TUnfoldSysV17.C
    // where the background systematic assigned below is NOT included in the input covariance matrix fVyy. 
    // This is done to keep fVyy diagonal, so that large background systematics can't bias the minrhoavg scan.
    // note that this means that background systematics are excluded from TUnfold's GetEmatrix().

    for (int i_C = 0; i_C < nCombinedBkgHists; ++i_C)
    {
      unfold.SubtractBackground(CombinedBkgHists[i_C],backgroundnames[i_C],1.0,backgrounduncertainties[i_C]);
    }

  }


  // Set alternative bias distribution
  if(freweight_bias_test) unfold.SetBias(GenVecHist_forBias);

  if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Input is set - RegMode for Unfolding is %d, AreaConstrain is %d. scaleBias is %f.\n", RegMode,useAreaConstrain,scaleBias);

  // Now set up the rebinning in two ways: After (A) and Before (B) unfolding.

  // Create bin map for combining the fine bins into the final bins when rebinning (A)fter unfolding
  int binMap[fNBinsX+2];
  CreateBinMap(generatorBinning->FindNode("ttbargen"), frebinfine, binMap);

  if(dbgUnf) {
    cout<<"binMap: ";
    for (int i = 0; i < fNBinsX+2; ++i)  cout << binMap[i] <<",";
    cout<<endl;
  }

  // const int binMap[fNBinsX+2] = {-1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,-1};

  // Set up the alternative unfolding with inputs rebinned (B)efore unfolding
  TH2D* ResMat_rebinnedB     = (TH2D*) this->GetAequiResMat()->Clone((Form("%s_rebinnedB",this->GetAequiResMat()->GetName())));
  TH1D* AequiM_rebinnedB     = (TH1D*) this->GetAequiMeas()->Clone((Form("%s_rebinnedB",this->GetAequiMeas()->GetName())));
  TH1D* GenVecHist_rebinnedB = (TH1D*) this->GetGenVecHist()->Clone((Form("%s_rebinnedB",this->GetGenVecHist()->GetName())));

  // rebin input histos
  if(frebinfine>1) {
    MatrixUnf::rebinMultidimMigrationMatrix(ResMat_rebinnedB,
                                            generatorBinning_rebinnedB->FindNode("ttbargen_rebinnedB"),
                                            detectorBinning_rebinnedB->FindNode("ttbarreco_rebinnedB"),
                                            generatorBinning->FindNode("ttbargen"),
                                            detectorBinning->FindNode("ttbarreco"),
                                            frebinfine);

    MatrixUnf::rebinMultidimensionalInput(AequiM_rebinnedB,
                                        detectorBinning_rebinnedB->FindNode("ttbarreco_rebinnedB"),
                                        detectorBinning->FindNode("ttbarreco"),
                                        frebinfine);

    MatrixUnf::rebinMultidimensionalInput(GenVecHist_rebinnedB,
                                        generatorBinning_rebinnedB->FindNode("ttbargen_rebinnedB"),
                                        generatorBinning->FindNode("ttbargen"),
                                        frebinfine);
  }

  // Set up the TUnfoldDensity object
  
  // Amandeep : Modifying the axis steering to only regularize on the inner axis
  // Original :
  // TUnfoldDensity unfold_rebinnedB(ResMat_rebinnedB,TUnfold::kHistMapOutputHoriz, REG, AREA, DENSITY, generatorBinning_rebinnedB, detectorBinning_rebinnedB, 0, "*[UOB]" );
  // Modified :
  TUnfoldDensity unfold_rebinnedB(ResMat_rebinnedB,TUnfold::kHistMapOutputHoriz, REG, AREA, DENSITY, generatorBinning_rebinnedB, detectorBinning_rebinnedB, 0, AxisSteering );
  // End

  unfold_rebinnedB.SetInput(AequiM_rebinnedB, scaleBias, oneOverZeroError);
  
  if(subtrBg) {
    for (int i_C = 0; i_C < nCombinedBkgHists; ++i_C)
    {
      CombinedBkgHists_rebinnedB[i_C] = (TH1D*) CombinedBkgHists[i_C]->Clone(fVariableName1+"_"+backgroundnames.at(i_C)+"bkg_rebinnedB");
      if (frebinfine > 1) {
        MatrixUnf::rebinMultidimensionalInput(CombinedBkgHists_rebinnedB[i_C],
                                            detectorBinning_rebinnedB->FindNode("ttbarreco_rebinnedB"),
                                            detectorBinning->FindNode("ttbarreco"),
                                            frebinfine);
      }
      unfold_rebinnedB.SubtractBackground(CombinedBkgHists_rebinnedB[i_C],backgroundnames[i_C],1.0,backgrounduncertainties[i_C]);
      // delete CombinedBkgHists[i_C];
      // delete CombinedBkgHists_rebinnedB[i_C];
    }
  }

  if(freweight_bias_test) unfold_rebinnedB.SetBias(GenVecHist_forBias_rebinnedB);

  // ***********
  // SYSTEMATICS -- ajeeta
  // ***********

  if(fsetSyst)
  {
    TObjArray *dummysysarray = (TObjArray*) GetAequiResSysHists()->Clone();

    for(Int_t k=0; k<dummysysarray->GetEntries(); k++){
      TH2D *dummySys = (TH2D*)((*dummysysarray)[k])->Clone();
      TH2D *dummySys_rebinnedB = (TH2D*)((*dummysysarray)[k])->Clone();

      std::string str, str2; 
      str        = dummySys->GetName();
      size_t pos = str.find("_AequiResMatSys");
      str2       = str.substr(0,pos);
      std::cout<<str2<<"_Syserr"<<std::endl;

      unfold.AddSysError(dummySys,Form("%s_Syserr",str2.c_str()),TUnfold::kHistMapOutputHoriz, TUnfoldSys::kSysErrModeMatrix);

      // rebin histo
      if (frebinfine > 1) {
        MatrixUnf::rebinMultidimMigrationMatrix(dummySys_rebinnedB,
                                                generatorBinning_rebinnedB->FindNode("ttbargen_rebinnedB"),
                                                detectorBinning_rebinnedB->FindNode("ttbarreco_rebinnedB"),
                                                generatorBinning->FindNode("ttbargen"),
                                                detectorBinning->FindNode("ttbarreco"),
                                                frebinfine);
      }
      unfold_rebinnedB.AddSysError(dummySys_rebinnedB,Form("%s_Syserr",str2.c_str()),TUnfold::kHistMapOutputHoriz, TUnfoldSys::kSysErrModeMatrix);
    }
   delete dummysysarray;
  }

  // Regularisation strength optimisation
  // Set up histograms for the Rho_i vectors and covariance matrices after rebinning
  const int nbinsrhoi             = fNBinsX / pow(frebinfine, dim);     // number of bins after rebinning
  const int numBinsPerCoefficient = 6;                                  // number of bins in the spin variable's distribution
  const int numCoefficients       = nbinsrhoi / numBinsPerCoefficient;  // number of different phase space regions we are measuring the variable in

  TH1D *Rhoi_rebinnedA    = new TH1D(fVariableName1+"Rhoi_rebinnedA","Global Corr Coeff (rebinnedA)",nbinsrhoi,0,nbinsrhoi);
  TH2D *Ematrix_rebinnedA = new TH2D(fVariableName1+"Ematrix_rebinnedA","Covariance Matrix (rebinnedA)",nbinsrhoi,0,nbinsrhoi, nbinsrhoi,0,nbinsrhoi);
  TH1D *Rhoi_rebinnedB    = new TH1D(fVariableName1+"Rhoi_rebinnedB","Global Corr Coeff (rebinnedB)",nbinsrhoi,0,nbinsrhoi);
  TH2D *Ematrix_rebinnedB = new TH2D(fVariableName1+"Ematrix_rebinnedB","Covariance Matrix (rebinnedB)",nbinsrhoi,0,nbinsrhoi, nbinsrhoi,0,nbinsrhoi);

  //Initialise variables used in the regularisation strength optimisation
  Double_t tauMinOpt=1;
  Double_t rhoAvgMinOpt=999;
  Double_t tauMinOpt_rebinnedB=1;
  Double_t rhoAvgMinOpt_rebinnedB=999;

  if(minrhoavg) {

    if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Minimising rhoAvg for Unfolding \n");

    Double_t Afb[numCoefficients], AfbErr[numCoefficients], Coef[numCoefficients], CoefErrOpt[numCoefficients], CoefErrTest[numCoefficients];
    std::vector<double> Afbvec, AfbErrvec, Coefvec, a_optimal_reco;
    TMatrixD m_E_rebinnedA(nbinsrhoi, nbinsrhoi);
    TMatrixD m_AFB(nbinsrhoi/2, nbinsrhoi/2);
    TMatrixD m_C(nbinsrhoi/2, nbinsrhoi/2);
   
    int nIt1 = 0;
    const int nScan1   = 31;
    double logtaustart = -15.;

    Double_t xRAvg1[nScan1], yRAvg1[nScan1];
    Double_t xAFBerr1[nScan1], yAFBerr1[numCoefficients][nScan1];
    Double_t xCerr1[nScan1], yCerr1[numCoefficients][nScan1];
    Double_t xCerrOpt1[nScan1], yCerrOpt1[numCoefficients][nScan1];
    Double_t ya0_1[numCoefficients][nScan1], ya1_1[numCoefficients][nScan1], ya2_1[numCoefficients][nScan1];
    Double_t yCerr0_1[numCoefficients][nScan1], yCerr1_1[numCoefficients][nScan1], yCerr2_1[numCoefficients][nScan1];
    Double_t xRAvg1_rebinnedB[nScan1], yRAvg1_rebinnedB[nScan1];

    for (Int_t i=0; i < nScan1; i++) {
      Double_t logtau = logtaustart*(1.-i/double(nScan1-1));
      Double_t tauReg = TMath::Power(10,logtau);

      unfold.DoUnfold(tauReg);

      if(frebinfine > 1)  {
        unfold.GetRhoI(Rhoi_rebinnedA,binMap);
      }
      else {
        unfold.GetRhoI(Rhoi_rebinnedA,0);
      }

      double rhoaverage = Rhoi_rebinnedA->Integral()/nbinsrhoi;
      if(rhoaverage>0 && rhoaverage<rhoAvgMinOpt){
        rhoAvgMinOpt = rhoaverage;
        tauMinOpt    = tauReg;
      }

      unfold_rebinnedB.DoUnfold(tauReg);
      unfold_rebinnedB.GetRhoI(Rhoi_rebinnedB,0);
      double rhoaverage_rebinnedB = Rhoi_rebinnedB->Integral()/nbinsrhoi;

      // Amandeep :
      std::cout << "frebinfine :: " << frebinfine << " dim :: " << dim << " fNBinsX :: " << fNBinsX << std::endl;
      std::cout << "Rhoi_rebinnedB->Integral() :: " << Rhoi_rebinnedB->Integral() << " nbinsrhoi " << nbinsrhoi << std::endl;
      std::cout << "rhoaverage_rebinnedB :: " << rhoaverage_rebinnedB << " " << unfold_rebinnedB.GetRhoAvg() << " " << unfold_rebinnedB.GetScanVariable(TAUMODE,0,"*[UO]") << std::endl;
      // End

      if (rhoaverage_rebinnedB > 0 && rhoaverage_rebinnedB < rhoAvgMinOpt_rebinnedB){
        rhoAvgMinOpt_rebinnedB = rhoaverage_rebinnedB;
        tauMinOpt_rebinnedB    = tauReg;
      }

      // Make some plots showing the dependence of the stat uncertainty in the measured coefficient on tau
      if(!ensTest) {
        TH1D* TUnfResult_rebinnedA = (TH1D*) unfold.GetOutput(fVariableName1+"TUnfResult_rebinnedA",0,"ttbargen","*[UO]",0);

        if(frebinfine>1)  {
          unfold.GetEmatrix(Ematrix_rebinnedA,binMap);
        }
        else {
          unfold.GetEmatrix(Ematrix_rebinnedA,0);
        }

        if (frebinfine > 1) {
          MatrixUnf::rebinMultidimensionalInput(TUnfResult_rebinnedA,
                                              generatorBinning_rebinnedB->FindNode("ttbargen_rebinnedB"),
                                              generatorBinning->FindNode("ttbargen"),
                                              frebinfine);
        }

        for (Int_t cmi = 0; cmi < nbinsrhoi; cmi++)
        {
          for (Int_t cmj = 0; cmj < nbinsrhoi; cmj++)
          {
            m_E_rebinnedA(cmi, cmj) = Ematrix_rebinnedA->GetBinContent(cmi + 1, cmj + 1);
          }
        }

        std::cout << "Computing the Afb w/ correlations" << std::endl;
        MatrixUnf::GetAfbWithCorrelations(TUnfResult_rebinnedA, m_E_rebinnedA, Afb, AfbErr, numCoefficients);
        std::cout << "Computing the Afb w/ correlations bin by bin" << std::endl;
        MatrixUnf::GetAfbWithCorrelationsBinByBin(TUnfResult_rebinnedA, m_E_rebinnedA, Afbvec, AfbErrvec, m_AFB, numCoefficients);
        std::cout << "Computing the optimal coefficient correlations" << std::endl;
        bool isCoefficient = MatrixUnf::GetOptimalCoefficient(Afbvec,
                                                              m_AFB,
                                                              GenVecHist_rebinnedB,
                                                              fVariableName1,
                                                              Coef,
                                                              CoefErrOpt,
                                                              CoefErrTest,
                                                              Coefvec,
                                                              m_C,
                                                              a_optimal_reco,
                                                              numCoefficients);
        if (!isCoefficient) {
          for (int i_var = 0; i_var < numCoefficients; i_var++) {
            Coef[i_var] = Afb[i_var];
            CoefErrTest[i_var] = AfbErr[i_var];
            CoefErrOpt[i_var] = AfbErr[i_var];
            int offset = i_var * numBinsPerCoefficient;
            int AfbOffset = offset / 2;
            double integral = TUnfResult_rebinnedA->Integral(1 + offset, offset + numBinsPerCoefficient);
            for (int i_bin = 0; i_bin < numBinsPerCoefficient / 2; ++i_bin) {
              a_optimal_reco[AfbOffset + i_bin] = (TUnfResult_rebinnedA->GetBinContent(AfbOffset + i_bin + 1) +
                                      TUnfResult_rebinnedA->GetBinContent(AfbOffset + numBinsPerCoefficient - i_bin)) /
                                      integral;
            }
          }
        }

        // factor to convert inclusive asymmetry to the coefficient
        double Afactor = MatrixUnf::AtoCfactor(fVariableName1);

        for (int j = 0; j < numCoefficients; j++) {
          AfbErr[j] = AfbErr[j] * Afactor;
        }

        if (rhoaverage > 0 && !ensTest) {
          xRAvg1[nIt1] = TMath::Log10(unfold.GetTau());
          yRAvg1[nIt1] = rhoaverage;
          xRAvg1_rebinnedB[nIt1] = TMath::Log10(unfold_rebinnedB.GetTau());
          yRAvg1_rebinnedB[nIt1] = rhoaverage_rebinnedB > 0 ? rhoaverage_rebinnedB : 1;
          xAFBerr1[nIt1] = TMath::Log10(unfold.GetTau());
          xCerr1[nIt1] = TMath::Log10(unfold.GetTau());
          xCerrOpt1[nIt1] = TMath::Log10(unfold.GetTau());
          for (int j = 0; j < numCoefficients; j++) {
            int AfbOffset = j * numBinsPerCoefficient / 2;
            yAFBerr1[j][nIt1] = AfbErr[j];
            yCerr1[j][nIt1] = CoefErrTest[j];
            yCerrOpt1[j][nIt1] = CoefErrOpt[j];
            ya0_1[j][nIt1] = a_optimal_reco.at(AfbOffset + 0);
            ya1_1[j][nIt1] = a_optimal_reco.at(AfbOffset + 1);
            ya2_1[j][nIt1] = a_optimal_reco.at(AfbOffset + 2);
            yCerr0_1[j][nIt1] = sqrt(m_C(AfbOffset + 0, AfbOffset + 0));
            yCerr1_1[j][nIt1] = sqrt(m_C(AfbOffset + 1, AfbOffset + 1));
            if (numBinsPerCoefficient > 2)
              yCerr2_1[j][nIt1] = sqrt(m_C(AfbOffset + 2, AfbOffset + 2));
          }
          nIt1++;
        }
        delete TUnfResult_rebinnedA;
      }
    }

    double logtauoffset = TMath::Log10(tauMinOpt_rebinnedB) - TMath::Log10(tauMinOpt);
    int nIt2         = 0;
    const int nScan2 = 401;
    double scanwidth = -2.*logtaustart/double(nScan1-1);
    logtaustart      = TMath::Log10(tauMinOpt) - scanwidth/2.; 

    Double_t xRAvg2[nScan2],yRAvg2[nScan2];
    Double_t xRAvg2_rebinnedB[nScan2],yRAvg2_rebinnedB[nScan2];

    for (Int_t i=0;i<nScan2;i++) {
      Double_t logtau = logtaustart + scanwidth*i/double(nScan2-1);
      Double_t tauReg = TMath::Power(10,logtau);
      unfold.DoUnfold(tauReg);

      if(frebinfine>1) {
        unfold.GetRhoI(Rhoi_rebinnedA,binMap);
      }
      else {
        unfold.GetRhoI(Rhoi_rebinnedA,0);
      } 

      double rhoaverage = Rhoi_rebinnedA->Integral()/nbinsrhoi;
      if(rhoaverage>0 && rhoaverage<rhoAvgMinOpt){
        rhoAvgMinOpt = rhoaverage;
        tauMinOpt = tauReg;
      }

      tauReg=TMath::Power(10,logtau+logtauoffset);
      unfold_rebinnedB.DoUnfold(tauReg);
      unfold_rebinnedB.GetRhoI(Rhoi_rebinnedB,0);
      double rhoaverage_rebinnedB = Rhoi_rebinnedB->Integral()/nbinsrhoi;
      //cout<<"rhoaverage_rebinnedB: "<<rhoaverage_rebinnedB<<" "<<unfold_rebinnedB.GetRhoAvg()<<" "<<unfold_rebinnedB.GetScanVariable(TAUMODE,0,"*[UO]")<<endl;
      if(rhoaverage_rebinnedB>0 && rhoaverage_rebinnedB<rhoAvgMinOpt_rebinnedB){
        rhoAvgMinOpt_rebinnedB = rhoaverage_rebinnedB;
        tauMinOpt_rebinnedB = tauReg;
      }

      if(rhoaverage>0 && !ensTest){
        xRAvg2[nIt2]=TMath::Log10(unfold.GetTau());
        yRAvg2[nIt2]=rhoaverage;
        xRAvg2_rebinnedB[nIt2]=TMath::Log10(unfold_rebinnedB.GetTau());
        yRAvg2_rebinnedB[nIt2]=rhoaverage_rebinnedB;
        nIt2++;
      }

    }

      /*
      TSpline *rhoScan=0;
      unfold.ScanTau(50,10e-12,1.,&rhoScan, TAUMODE, 0, "*[UO]");
      tauMinOpt = unfold.GetTau();

      //now fine-tune the optimisation using smaller step size (may change to Minuit minimisation later)
      unfold.ScanTau(50,tauMinOpt/2.,tauMinOpt*2.,&rhoScan, TAUMODE, 0, "*[UO]");
      if( std::fabs( unfold.GetTau()/(tauMinOpt/2.) - 1. ) < 0.1 ||  std::fabs( unfold.GetTau()/(tauMinOpt*2.) - 1. ) < 0.1  ) printf("MatrixUnf::prepRunUnfolding: Error, tau scan zoomed in too far: %.3f %.3f \n", unfold.GetTau()/(tauMinOpt/2.) - 1. , unfold.GetTau()/(tauMinOpt*2.) - 1. );
      tauMinOpt = unfold.GetTau();
      //rhoAvgMinOpt = unfold.GetRhoAvg();
      rhoAvgMinOpt = unfold.GetScanVariable(TAUMODE,0,"*[UO]");
     */

      //TH1 * rihist = unfold.GetRhoIstatbgr("rihist", 0, 0, "*[UO]");
      //rihist->Print("all");

    if (!ensTest) {
      // fResultUnfHists->AddAt(rhoijAVG_curve, 22);

      TGraph* rhoijAVG_curve1 = new TGraph(nIt1, xRAvg1, yRAvg1);
      rhoijAVG_curve1->SetName(fVariableName1 + "ScanOfRhoAvg1");
      rhoijAVG_curve1->GetXaxis()->SetTitle("Log10(#tau)");
      rhoijAVG_curve1->GetYaxis()->SetTitle("Avg global #rho");
      rhoijAVG_curve1->Write();

      TGraph *AFBerr_curve1, *Cerr_curve1, *CerrOpt_curve1;
      TGraph *a0_curve1, *a1_curve1, *a2_curve1;
      TGraph *Cerr0_curve1, *Cerr1_curve1, *Cerr2_curve1;
      for (int j = 0; j < numCoefficients; j++) {
        AFBerr_curve1 = new TGraph(nIt1, xAFBerr1, yAFBerr1[j]);
        AFBerr_curve1->SetName(Form("%s %i ScanOfAFBerr1", fVariableName1.Data(), j));
        AFBerr_curve1->GetXaxis()->SetTitle("Log10(#tau)");
        AFBerr_curve1->GetYaxis()->SetTitle("AFBerr");
        AFBerr_curve1->Write();

        Cerr_curve1 = new TGraph(nIt1, xCerr1, yCerr1[j]);
        Cerr_curve1->SetName(Form("%s %i ScanOfCerr1", fVariableName1.Data(), j));
        Cerr_curve1->GetXaxis()->SetTitle("Log10(#tau)");
        Cerr_curve1->GetYaxis()->SetTitle("Cerr");
        Cerr_curve1->Write();

        CerrOpt_curve1 = new TGraph(nIt1, xCerrOpt1, yCerrOpt1[j]);
        CerrOpt_curve1->SetName(Form("%s %i ScanOfCerrOpt1", fVariableName1.Data(), j));
        CerrOpt_curve1->GetXaxis()->SetTitle("Log10(#tau)");
        CerrOpt_curve1->GetYaxis()->SetTitle("CerrOpt");
        CerrOpt_curve1->Write();

        a0_curve1 = new TGraph(nIt1, xCerr1, ya0_1[j]);
        a0_curve1->SetName(Form("%s %i ScanOfa0_1", fVariableName1.Data(), j));
        a0_curve1->GetXaxis()->SetTitle("Log10(#tau)");
        a0_curve1->GetYaxis()->SetTitle("a0");
        a0_curve1->Write();

        a1_curve1 = new TGraph(nIt1, xCerr1, ya1_1[j]);
        a1_curve1->SetName(Form("%s %i ScanOfa1_1", fVariableName1.Data(), j));
        a1_curve1->GetXaxis()->SetTitle("Log10(#tau)");
        a1_curve1->GetYaxis()->SetTitle("a1");
        a1_curve1->Write();

        a2_curve1 = new TGraph(nIt1, xCerr1, ya2_1[j]);
        a2_curve1->SetName(Form("%s %i ScanOfa2_1", fVariableName1.Data(), j));
        a2_curve1->GetXaxis()->SetTitle("Log10(#tau)");
        a2_curve1->GetYaxis()->SetTitle("a2");
        a2_curve1->Write();

        Cerr0_curve1 = new TGraph(nIt1, xCerr1, yCerr0_1[j]);
        Cerr0_curve1->SetName(Form("%s %i ScanOfCerr0_1", fVariableName1.Data(), j));
        Cerr0_curve1->GetXaxis()->SetTitle("Log10(#tau)");
        Cerr0_curve1->GetYaxis()->SetTitle("Cerr0");
        Cerr0_curve1->Write();

        Cerr1_curve1 = new TGraph(nIt1, xCerr1, yCerr1_1[j]);
        Cerr1_curve1->SetName(Form("%s %i ScanOfCerr1_1", fVariableName1.Data(), j));
        Cerr1_curve1->GetXaxis()->SetTitle("Log10(#tau)");
        Cerr1_curve1->GetYaxis()->SetTitle("Cerr1");
        Cerr1_curve1->Write();

        Cerr2_curve1 = new TGraph(nIt1, xCerr1, yCerr2_1[j]);
        Cerr2_curve1->SetName(Form("%s %i ScanOfCerr2_1", fVariableName1.Data(), j));
        Cerr2_curve1->GetXaxis()->SetTitle("Log10(#tau)");
        Cerr2_curve1->GetYaxis()->SetTitle("Cerr2");
        Cerr2_curve1->Write();
      }

      TGraph* rhoijAVG_curve2 = new TGraph(nIt2, xRAvg2, yRAvg2);
      rhoijAVG_curve2->SetName(Form("%s ScanOfRhoAvg2", fVariableName1.Data()));
      rhoijAVG_curve2->GetXaxis()->SetTitle("Log10(#tau)");
      rhoijAVG_curve2->GetYaxis()->SetTitle("Avg global #rho");
      rhoijAVG_curve2->Write();

      fResultUnfHists->AddAt(rhoijAVG_curve2, 22);

      TGraph* rhoijAVG_curve1_rebinnedB = new TGraph(nIt1, xRAvg1_rebinnedB, yRAvg1_rebinnedB);
      rhoijAVG_curve1_rebinnedB->SetName(fVariableName1 + "ScanOfRhoAvg1_rebinnedB");
      rhoijAVG_curve1_rebinnedB->GetXaxis()->SetTitle("Log10(#tau)");
      rhoijAVG_curve1_rebinnedB->GetYaxis()->SetTitle("Avg global #rho");
      rhoijAVG_curve1_rebinnedB->Write();

      TGraph* rhoijAVG_curve2_rebinnedB = new TGraph(nIt2, xRAvg2_rebinnedB, yRAvg2_rebinnedB);
      rhoijAVG_curve2_rebinnedB->SetName(fVariableName1 + "ScanOfRhoAvg2_rebinnedB");
      rhoijAVG_curve2_rebinnedB->GetXaxis()->SetTitle("Log10(#tau)");
      rhoijAVG_curve2_rebinnedB->GetYaxis()->SetTitle("Avg global #rho");
      rhoijAVG_curve2_rebinnedB->Write();

      // rhoScan->Write();
      // fResultUnfHists->AddAt(rhoScan, 22);

      // save minimum in rhoAvg corresponding to the one minimized using rhoAVG
      Double_t xR[1], yR[1];
      yR[0] = rhoAvgMinOpt;
      xR[0] = TMath::Log10(tauMinOpt);
      TGraph* rhoijAVG = new TGraph(1, xR, yR);
      rhoijAVG->SetName(fVariableName1 + "RhoAvg");
      rhoijAVG->GetXaxis()->SetTitle("Log10(#tau)");
      rhoijAVG->GetYaxis()->SetTitle("Avg global #rho");
      rhoijAVG->Write();

      fResultUnfHists->AddAt(rhoijAVG, 23);

      Double_t xR_rebinnedB[1], yR_rebinnedB[1];
      yR_rebinnedB[0] = rhoAvgMinOpt_rebinnedB;
      xR_rebinnedB[0] = TMath::Log10(tauMinOpt_rebinnedB);
      TGraph* rhoijAVG_rebinnedB = new TGraph(1, xR_rebinnedB, yR_rebinnedB);
      rhoijAVG_rebinnedB->SetName(fVariableName1 + "RhoAvg_rebinnedB");
      rhoijAVG_rebinnedB->GetXaxis()->SetTitle("Log10(#tau)");
      rhoijAVG_rebinnedB->GetYaxis()->SetTitle("Avg global #rho");
      rhoijAVG_rebinnedB->Write();

      delete rhoijAVG_curve1;
      delete AFBerr_curve1;
      delete Cerr_curve1;
      delete CerrOpt_curve1;

      delete a0_curve1;
      delete a1_curve1;
      delete a2_curve1;
      delete Cerr0_curve1;
      delete Cerr1_curve1;
      delete Cerr2_curve1;

      delete rhoijAVG_curve2;
      delete rhoijAVG;
      delete rhoijAVG_curve1_rebinnedB;
      delete rhoijAVG_curve2_rebinnedB;
      delete rhoijAVG_rebinnedB;
    }
  }
  
  
  // Only the minrhoavg regularisation strength optimisation method is compatible with using fine bins for the unfolding, so doScanLcurve is not fully supported
  bool doScanLcurve = false;
  if(doScanLcurve) {

    if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Scanning L-curve for Unfolding \n");
    
    Int_t nScan = 30;
    Int_t iBest = 1;
    Int_t iBest_rebinnedB = 1;
    TSpline *logTauX, *logTauY;
    TGraph *lCurve;
    
    //Int_t nBinsOfScanHists = 1000;
    //if(ensTest)nBinsOfScanHists=1;
    
    if(!minrhoavg){
      
      // reset tau for the best value, since the scan for average global correlations                            
      // minimizes the interval in tau (andy) 
      // do not want to reset to 0 but want to keep what is passed on in the arguments
      
      // this method scans the parameter tau and finds the kink in the L curve                                   
      // finally, the unfolding is done for the "best" choice of tau                                             
      iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
      iBest_rebinnedB=unfold_rebinnedB.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);

    }
    
    else{
      
      // minimize based on rho
      // this overwrites tauMin and tauMax based on the minimum found from minimization of rho AVG

      tauMin = 0.01*tauMinOpt;
      tauMax = 100.*tauMinOpt;
      
      iBest=unfold.ScanLcurve(nScan, tauMin, tauMax, &lCurve, &logTauX, &logTauY);

      tauMin = 0.01*tauMinOpt_rebinnedB;
      tauMax = 100.*tauMinOpt_rebinnedB;
      
      iBest_rebinnedB=unfold_rebinnedB.ScanLcurve(nScan, tauMin, tauMax, &lCurve, &logTauX, &logTauY);
      
    }

    if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: iBest after scanning of L-curve between %.3e and %.3e is %d with tau=%.3e and rho=%.3f , c.f. tauOpt=%.3e and rhoOpt=%.3f (ratios: %.3f, %.3f )\n", tauMin, tauMax, iBest, unfold.GetTau(),unfold.GetScanVariable(TAUMODE,0,"*[UO]"),tauMinOpt,rhoAvgMinOpt, unfold.GetTau()/tauMinOpt, unfold.GetScanVariable(TAUMODE,0,"*[UO]")/rhoAvgMinOpt );
    if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: ScanLcurve Chi2 of matrixTerm %.2f -- Chi2 of regTerm %.2f, NDF %d --> tot Chi2/n.d.f = %.2f\n", unfold.GetChi2A(), unfold.GetChi2L(), unfold.GetNdf(), ((unfold.GetChi2A()+unfold.GetChi2L())/((double) unfold.GetNdf())));

    if(!ensTest) {
      
      // save point corresponding to the kink in the L curve as TGraph 
      
      Double_t t[1],x[1],y[1];
      logTauX->GetKnot(iBest,t[0],x[0]);
      logTauY->GetKnot(iBest,t[0],y[0]);
      
      lCurve->SetName(fVariableName1+"lCurve");
      lCurve->GetYaxis()->SetTitle("log_{10}(curvature), L");
      lCurve->GetXaxis()->SetTitle("log_{10}(#chi^{2}) of A");
      lCurve->Write();
      
      // save point with smallest correlation as a graph
      Double_t logTau=TMath::Log10(t[0]);
      TGraph *lCurveUser = new TGraph(1,x,y);
      lCurveUser->SetName(fVariableName1+"lCurveUser");
      lCurveUser->GetYaxis()->SetTitle("log_{10}(curvature), L");
      lCurveUser->GetXaxis()->SetTitle("log_{10}(#chi^{2}) of A");
      lCurveUser->Write();
      
      
      TGraph *logTauXuser = new TGraph(1,t,x);
      logTauXuser->SetName(fVariableName1+"logTauXuser");
      logTauXuser->Write();
      
      // save minimum in rhoAvg corresponding to the one in the l-curve                                          
      Double_t xL[1],yL[1];
      yL[0] = unfold.GetScanVariable(TAUMODE,0,"*[UO]");
      xL[0] = TMath::Log10(unfold.GetTau());
      TGraph *bestLRhoAvg=new TGraph(1,xL,yL);
      bestLRhoAvg->SetName(fVariableName1+"bestLRhoAvg");
      bestLRhoAvg->Write();

      if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Adding lCurve to ObjArray \n"); 
      logTauX->Draw();
      fResultUnfHists->AddAt(logTauX, 7);
      
      logTauXuser->Draw();
      fResultUnfHists->AddAt(logTauXuser, 8);
      lCurve->Draw();
      fResultUnfHists->AddAt(lCurve, 9);
      lCurveUser->Draw();
      fResultUnfHists->AddAt(lCurveUser, 10);

      delete lCurveUser;
      delete logTauXuser;
      delete bestLRhoAvg;
    }
    
  }
  
  if(minrhoavg){
    //now do the unfolding with the rhoAvgMin optimal tau value
    unfold.DoUnfold(tauMinOpt);
    //unfold.DoUnfold(tauMinOpt*30.);
    //unfold.DoUnfold(1e-3);
    //if(dbgUnf) cout<<"ScanTau Rho Variables: "<<unfold.GetScanVariable(0,0,"*[UO]")<<" "<<unfold.GetScanVariable(1,0,"*[UO]")<<" "<<unfold.GetScanVariable(2,0,"*[UO]")<<" "<<unfold.GetScanVariable(3,0,"*[UO]")<<" "<<unfold.GetScanVariable(4,0,"*[UO]")<<" "<<unfold.GetScanVariable(5,0,"*[UO]")<<" "<<tauMinOpt<<" "<<unfold.GetRhoAvg()<<" "<<unfold.GetRhoMax()<<" "<<rhoAvgMinOpt<<endl; 
    if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Did the unfolding with tau = %.3e, Average Global Correlation = %.3f\n",unfold.GetTau(),rhoAvgMinOpt);
    if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: MinRhoAvg Chi2 of matrixTerm %.5f -- Chi2 of regTerm %.5f, NDF %d --> tot Chi2/n.d.f = %.5f\n", unfold.GetChi2A(), unfold.GetChi2L(), unfold.GetNdf(), ((unfold.GetChi2A()+unfold.GetChi2L())/((double) unfold.GetNdf())));

    unfold_rebinnedB.DoUnfold(tauMinOpt_rebinnedB);
    //unfold_rebinnedB.DoUnfold(tauMinOpt_rebinnedB*30.);
    //unfold_rebinnedB.DoUnfold(1e-4);
    if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Did the unfolding_rebinnedB with tau = %.3e, Average Global Correlation = %.3f\n",unfold_rebinnedB.GetTau(),rhoAvgMinOpt_rebinnedB);
    if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: MinRhoAvg_rebinnedB Chi2 of matrixTerm %.5f -- Chi2 of regTerm %.5f, NDF %d --> tot Chi2/n.d.f = %.5f\n", unfold_rebinnedB.GetChi2A(), unfold_rebinnedB.GetChi2L(), unfold_rebinnedB.GetNdf(), ((unfold_rebinnedB.GetChi2A()+unfold_rebinnedB.GetChi2L())/((double) unfold_rebinnedB.GetNdf())));
  }

  // get unfolded results
  TH1D* TUnfResult           = (TH1D*) unfold.GetOutput(fVariableName1+"TUnfResult",0,"ttbargen","*[UO]",0);
  TH1D* TUnfResult_rebinnedB = (TH1D*) unfold_rebinnedB.GetOutput(fVariableName1+"TUnfResult_rebinnedB",0,"ttbargen_rebinnedB","*[UO]",0);
  TH1D* TUnfResult_rebinnedA;

  // get error matrices (data stat errors only)
  TH2D *ErrMatrix = new TH2D(fVariableName1+"ErrMatrixData","Covariance Matrix",fNBinsX,0.,double(fNBinsX),fNBinsX,0.,double(fNBinsX));
  unfold.GetEmatrix(ErrMatrix);

  unfold_rebinnedB.GetEmatrix(Ematrix_rebinnedB);

  if(frebinfine>1)  {
    unfold.GetEmatrix(Ematrix_rebinnedA,binMap);
  }
  else {
    unfold.GetEmatrix(Ematrix_rebinnedA,0);
  }

  // calculate TUnfResult_rebinnedA (rebinned after unfolding), by rebinning TUnfResult
  TUnfResult_rebinnedA = (TH1D*) TUnfResult->Clone(fVariableName1+"TUnfResult_rebinnedA");

  if (frebinfine > 1) {
    MatrixUnf::rebinMultidimensionalInput(TUnfResult_rebinnedA,
                                        generatorBinning_rebinnedB->FindNode("ttbargen_rebinnedB"),
                                        generatorBinning->FindNode("ttbargen"),
                                        frebinfine);
  }

  for (Int_t cmi = 0; cmi < nbinsrhoi; cmi++)
  {
    TUnfResult_rebinnedA->SetBinError(cmi + 1, sqrt( Ematrix_rebinnedA->GetBinContent(cmi + 1, cmi + 1) ) );
    if(dbgUnf && !ensTest) printf("MatrixUnf::prepRunUnfolding: %s bin %i: TUnfRebinnedA: %.3f+-%.3f \n", fVariableName1.Data(), cmi+1, TUnfResult_rebinnedA->GetBinContent(cmi+1), TUnfResult_rebinnedA->GetBinError(cmi+1));
  }

  for(Int_t i =1; i <= TUnfResult_rebinnedB->GetNbinsX(); i++) {
    if(dbgUnf && !ensTest) printf("MatrixUnf::prepRunUnfolding: %s bin %i: TUnfRebinnedB: %.3f+-%.3f \n", fVariableName1.Data(), i, TUnfResult_rebinnedB->GetBinContent(i), TUnfResult_rebinnedB->GetBinError(i));
  }


  // Save output histos

  TString Tit1 = Form("Unfolded Cross section in %s (not corrected by lumi,binWidth,fBR)", fVariableName1.Data()); 
  TUnfResult->SetTitle(Tit1.Data()); 
  TUnfResult->GetXaxis()->SetTitle(fYAxisTitle);
  TString unfTit = Form("(uncor.) N(part.) in %s", fVariableName1.Data()); 
  TUnfResult->GetYaxis()->SetTitle(unfTit.Data());
  
  fResultUnfHists->AddAt(TUnfResult, 14);
  fResultUnfHists->AddAt(TUnfResult_rebinnedA, 35);
  fResultUnfHists->AddAt(TUnfResult_rebinnedB, 38);
  
  if(ensTest) {
    pseudoUnfResults.push_back(TUnfResult);
    pseudoUnfResultsRebinnedA.push_back(TUnfResult_rebinnedA);
    pseudoUnfResultsRebinnedB.push_back(TUnfResult_rebinnedB);
  }

  if(!ensTest){

    // Check closure (delta should be zero when unfolding using fSelfConsistency=true)
    for(Int_t i =1; i <= this->GetGenVecHist()->GetNbinsX(); i++) {
      if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: %s bin %i: Gen: %.3f+-%.3f \n", fVariableName1.Data(), i, this->GetGenVecHist()->GetBinContent(i), this->GetGenVecHist()->GetBinError(i));
      if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: %s bin %i: delta(Unfolded-Gen): %.2g \n", fVariableName1.Data(), i, TUnfResult->GetBinContent(i) - this->GetGenVecHist()->GetBinContent(i));
    }

  }

  // **************************************
  // Calculate asymmetries and coefficients
  // **************************************

  bool calculateCoefficients = true;
  if(calculateCoefficients) {

    // Calculate asymmetry both inclusively (GetAfbWithCorrelations) 
    // and for opposite pairs of bins (GetAfbWithCorrelationsBinByBin) in order to test the calculated covariance matrix
    
    // Store results in the AfbResult histogram (which is used in the pseudodata tests) 
    // where the left half contains the binned asymmetry results, 
    // and the last 3 bins of the right half contain the normalised sum of all bins, the measured coefficient, and the inclusive asymmetry.
    
    if(!ensTest) {
      std::cout << "filling gen asyms and coefficient " << std::endl;
      //TH1D* AfbResult_gen   = (TH1D*)GenVecHist_rebinnedB->Clone(fVariableName1 + "AfbResult_gen");
      TH1D* AfbResult_gen     = new TH1D(fVariableName1 + "AfbResult_gen", "", numCoefficients * (numBinsPerCoefficient + 2) + 20, 0, 1);
      //TH1D* OtherResult_gen = (TH1D*)GenVecHist_rebinnedB->Clone(fVariableName1 + "OtherResult_gen");
      TH1D* OtherResult_gen   = new TH1D(fVariableName1 + "OtherResult_gen", "", numCoefficients * (numBinsPerCoefficient + 2) + 20, 0, 1);
      MatrixUnf::fillAsymsAndCoefficient(GenVecHist_rebinnedB, nullptr, GenVecHist_rebinnedB, fVariableName1, AfbResult_gen, OtherResult_gen, ensTest, true, numCoefficients);
    }

    //TH1D* AfbResult_rebinnedA   = (TH1D*)TUnfResult_rebinnedA->Clone(fVariableName1 + "AfbResult_rebinnedA");
    TH1D* AfbResult_rebinnedA     = new TH1D(fVariableName1 + "AfbResult_rebinnedA", "", numCoefficients * (numBinsPerCoefficient + 2) + 20, 0, 1);
    //TH1D* OtherResult_rebinnedA = (TH1D*)TUnfResult_rebinnedA->Clone(fVariableName1 + "OtherResult_rebinnedA");
    TH1D* OtherResult_rebinnedA   = new TH1D(fVariableName1 + "OtherResult_rebinnedA", "", numCoefficients * (numBinsPerCoefficient + 2) + 20, 0, 1);
    MatrixUnf::fillAsymsAndCoefficient(TUnfResult_rebinnedA, Ematrix_rebinnedA, GenVecHist_rebinnedB, fVariableName1, AfbResult_rebinnedA, OtherResult_rebinnedA, ensTest, false, numCoefficients);
    OtherResult_rebinnedA->SetBinContent(0, TMath::Log10(tauMinOpt));
    OtherResult_rebinnedA->SetBinContent(numCoefficients * 2 + 1, rhoAvgMinOpt);


    //TH1D* AfbResult_rebinnedB = (TH1D*)TUnfResult_rebinnedB->Clone(fVariableName1 + "AfbResult_rebinnedB");
    TH1D* AfbResult_rebinnedB = new TH1D(fVariableName1 + "AfbResult_rebinnedB", "", numCoefficients * (numBinsPerCoefficient + 2) + 20, 0, 1);
    //TH1D* OtherResult_rebinnedB = (TH1D*)TUnfResult_rebinnedB->Clone(fVariableName1 + "OtherResult_rebinnedB");
    TH1D* OtherResult_rebinnedB = new TH1D(fVariableName1 + "OtherResult_rebinnedB", "", numCoefficients * (numBinsPerCoefficient + 2) + 20, 0, 1);
    MatrixUnf::fillAsymsAndCoefficient(TUnfResult_rebinnedB, Ematrix_rebinnedB, GenVecHist_rebinnedB, fVariableName1, AfbResult_rebinnedB, OtherResult_rebinnedB, ensTest, false, numCoefficients);
    OtherResult_rebinnedB->SetBinContent(0, TMath::Log10(tauMinOpt_rebinnedB));
    OtherResult_rebinnedB->SetBinContent(numCoefficients * 2 + 1, rhoAvgMinOpt_rebinnedB);

    if(doChi2tests) {
      printf("MatrixUnf::prepRunUnfolding: Doing Chi2 tests\n");
      //double chi2unfRB, chi2unfRB_rebinnedA, chi2unfRB_rebinnedB;
      std::vector<double> chi2unf;
      MatrixUnf::Chi2testsForPEs(unfold, unfold_rebinnedB, TUnfResult_rebinnedA, chi2unf);

      // note: using bin error here as a lazy way to store more variables without creating a new histogram
      OtherResult_rebinnedA->SetBinContent(numCoefficients * 2 + 5, chi2unf.at(0));
      OtherResult_rebinnedA->SetBinContent(numCoefficients * 2 + 4, chi2unf.at(3));
      OtherResult_rebinnedA->SetBinContent(numCoefficients * 2 + 3, unfold.GetChi2A());
      OtherResult_rebinnedA->SetBinError(numCoefficients * 2 + 5, chi2unf.at(1));
      OtherResult_rebinnedA->SetBinError(numCoefficients * 2 + 4, chi2unf.at(4));
      OtherResult_rebinnedA->SetBinError(numCoefficients * 2 + 3, unfold.GetChi2L());

      OtherResult_rebinnedB->SetBinContent(numCoefficients * 2 + 5, chi2unf.at(2));
      OtherResult_rebinnedB->SetBinContent(numCoefficients * 2 + 4, chi2unf.at(5));
      OtherResult_rebinnedB->SetBinContent(numCoefficients * 2 + 3, unfold_rebinnedB.GetChi2A());
      OtherResult_rebinnedB->SetBinError(numCoefficients * 2 + 5, 100 + chi2unf.at(0) - chi2unf.at(2));
      OtherResult_rebinnedB->SetBinError(numCoefficients * 2 + 4, 100 + chi2unf.at(3) - chi2unf.at(5));
      OtherResult_rebinnedB->SetBinError(numCoefficients * 2 + 3, unfold_rebinnedB.GetChi2L());
    }


    if(ensTest) {
      pseudoUnfAfbResultsRebinnedA.push_back(AfbResult_rebinnedA);
      pseudoUnfAfbResultsRebinnedB.push_back(AfbResult_rebinnedB);
      pseudoUnfOtherResultsRebinnedA.push_back(OtherResult_rebinnedA);
      pseudoUnfOtherResultsRebinnedB.push_back(OtherResult_rebinnedB);
    }

    fResultUnfHists->AddAt(AfbResult_rebinnedA, 36);
    fResultUnfHists->AddAt(AfbResult_rebinnedB, 39);

  }// calculateCoefficients


  // Extra results and tests
  if(!ensTest){

    TH1D *FoldedBack = (TH1D*) unfold.GetFoldedOutput(fVariableName1+"InputFoldedBack",0,"ttbarreco","*[UO]",1);
    fResultUnfHists->AddAt(FoldedBack, 11);

    // note, GetRhoIJ, GetEmatrix, and GetInputInverseEmatrix are the old TUnfold functions 
    // (so they only include data stat uncertainties)
    TH2D *Corrmatrix = new TH2D(fVariableName1+"Corrmatrix","Correlation Matrix",fNBinsX,0.,double(fNBinsX),fNBinsX,0.,double(fNBinsX));
    unfold.GetRhoIJ(Corrmatrix);
    fResultUnfHists->AddAt(Corrmatrix, 12);
     
    TH2D *invErrMatrixInputData = new TH2D(fVariableName1+"invErrMatrixInputData","Inverse Covariance Matrix (Errors from meas. Input)",fNBinsY,0.,double(fNBinsY),fNBinsY,0.,double(fNBinsY));
    unfold.GetInputInverseEmatrix(invErrMatrixInputData);
    fResultUnfHists->AddAt(invErrMatrixInputData, 13);


    TH1D *Rhoi = new TH1D(fVariableName1+"Rhoi","Global Corr Coeff",fNBinsX,0,fNBinsX);
    TH2D *InvE = new TH2D(fVariableName1+"InvE","Inverse E matrix",fNBinsX,0,fNBinsX, fNBinsX,0,fNBinsX);
    unfold.GetRhoI(Rhoi,0,InvE);

    TH2D *InvE_rebinnedA = new TH2D(fVariableName1+"InvE_rebinnedA","Inverse E matrix (rebinnedA)",nbinsrhoi,0,nbinsrhoi, nbinsrhoi,0,nbinsrhoi);
    TH2D *Corrmatrix_rebinnedA = new TH2D(fVariableName1+"Corrmatrix_rebinnedA","Correlation matrix (rebinnedA)",nbinsrhoi,0,nbinsrhoi, nbinsrhoi,0,nbinsrhoi);

    TH2D *InvE_rebinnedB = new TH2D(fVariableName1+"InvE_rebinnedB","Inverse E matrix (rebinnedB)",nbinsrhoi,0,nbinsrhoi, nbinsrhoi,0,nbinsrhoi);
    TH2D *Corrmatrix_rebinnedB = new TH2D(fVariableName1+"Corrmatrix_rebinnedB","Correlation matrix (rebinnedB)",nbinsrhoi,0,nbinsrhoi, nbinsrhoi,0,nbinsrhoi);

    if(frebinfine>1)  {
      unfold.GetRhoI(Rhoi_rebinnedA,binMap,InvE_rebinnedA);
      unfold.GetRhoIJ(Corrmatrix_rebinnedA,binMap);
    }
    else {
      unfold.GetRhoI(Rhoi_rebinnedA,0,InvE_rebinnedA);
      unfold.GetRhoIJ(Corrmatrix_rebinnedA,0);
    }

    unfold_rebinnedB.GetRhoI(Rhoi_rebinnedB,0,InvE_rebinnedB);
    unfold_rebinnedB.GetRhoIJ(Corrmatrix_rebinnedB,0);

    if(doChi2tests) MatrixUnf::Chi2tests(unfold, unfold_rebinnedB, TUnfResult, TUnfResult_rebinnedA, TUnfResult_rebinnedB);
  }


  // Only do corrections and systematics if not running pseudodata tests
  if(!ensTest){

    // Correct the result for binWidth, BR, Lumi & store it
    TH1D* TUnfResultCorXsec = (TH1D*)TUnfResult->Clone(fVariableName1+"TUnfResultCor");
    GetAbsoluteCrossSec(TUnfResultCorXsec, fBR, Lumi);

    TH1D* TUnfResultCorXsec_rebinnedA = (TH1D*)TUnfResult_rebinnedA->Clone(fVariableName1+"TUnfResultCor_rebinnedA");
    GetAbsoluteCrossSec(TUnfResultCorXsec_rebinnedA, fBR, Lumi);

    TH1D* TUnfResultCorXsec_rebinnedB = (TH1D*)TUnfResult_rebinnedB->Clone(fVariableName1+"TUnfResultCor_rebinnedB");
    GetAbsoluteCrossSec(TUnfResultCorXsec_rebinnedB, fBR, Lumi);


    TString Tit2    = Form("Unfolded Cross section in %s (corrected)", fVariableName1.Data()); 
    TString xsecTit = Form("d#sigma/d%s [nb]", fVariableName1.Data()); 

    TUnfResultCorXsec->SetTitle(Tit2.Data()); 
    TUnfResultCorXsec->GetXaxis()->SetTitle(fYAxisTitle);
    TUnfResultCorXsec->GetYaxis()->SetTitle(xsecTit.Data());

    TUnfResultCorXsec_rebinnedA->SetTitle(Tit2.Data()); 
    TUnfResultCorXsec_rebinnedA->GetXaxis()->SetTitle(fYAxisTitle);
    TUnfResultCorXsec_rebinnedA->GetYaxis()->SetTitle(xsecTit.Data());

    TUnfResultCorXsec_rebinnedB->SetTitle(Tit2.Data()); 
    TUnfResultCorXsec_rebinnedB->GetXaxis()->SetTitle(fYAxisTitle);
    TUnfResultCorXsec_rebinnedB->GetYaxis()->SetTitle(xsecTit.Data());

    fResultUnfHists->AddAt(TUnfResultCorXsec_rebinnedA, 15);

    // Get theory Xsec
    TH1D* TheoryXsec = (TH1D*) this->GetGenVecHist()->Clone(fVariableName1+"TheoryXsec");
    if (frebinfine > 1) {
      MatrixUnf::rebinMultidimensionalInput(TheoryXsec,
                                        generatorBinning_rebinnedB->FindNode("ttbargen_rebinnedB"),
                                        generatorBinning->FindNode("ttbargen"),
                                        frebinfine);
    }
    GetAbsoluteCrossSec(TheoryXsec, fBR, Lumi);
    TString Tit3 = Form("Theory Cross Section %s (corrected)", fVariableName1.Data()); 
    TheoryXsec->SetTitle(Tit3.Data()); 
    TheoryXsec->GetXaxis()->SetTitle(fYAxisTitle);
    TString ThxsecTit = Form("d#sigma/d%s [nb]", fVariableName1.Data()); 
    TheoryXsec->GetYaxis()->SetTitle(ThxsecTit.Data());
    fResultUnfHists->AddAt(TheoryXsec, 24);

    TH2D* ErrMatrixCorXsec         = (TH2D*) ErrMatrix->Clone( fVariableName1+"ErrMatrixDataCor" );
    TH2D* EmatrixCorXsec_rebinnedA = (TH2D*) Ematrix_rebinnedA->Clone( fVariableName1+"EmatrixCor_rebinnedA" );
    TH2D* EmatrixCorXsec_rebinnedB = (TH2D*) Ematrix_rebinnedB->Clone( fVariableName1+"EmatrixCor_rebinnedB" );

    GetAbsoluteCovarianceMatrix(TUnfResultCorXsec, ErrMatrixCorXsec, fBR, Lumi);
    GetAbsoluteCovarianceMatrix(TUnfResultCorXsec_rebinnedA, EmatrixCorXsec_rebinnedA, fBR, Lumi);
    GetAbsoluteCovarianceMatrix(TUnfResultCorXsec_rebinnedB, EmatrixCorXsec_rebinnedB, fBR, Lumi);

    TH2D *ErrMatrixNorm         = new TH2D(fVariableName1+"ErrMatrixDataNorm","Covariance Matrix", fNBinsX,0.,double(fNBinsX),fNBinsX,0.,double(fNBinsX));
    TH2D *EmatrixNorm_rebinnedA = new TH2D(fVariableName1+"EmatrixNorm_rebinnedA","Covariance matrixNorm (rebinnedA)",nbinsrhoi,0,nbinsrhoi, nbinsrhoi,0,nbinsrhoi);
    TH2D *EmatrixNorm_rebinnedB = new TH2D(fVariableName1+"EmatrixNorm_rebinnedB","Covariance matrixNorm (rebinnedB)",nbinsrhoi,0,nbinsrhoi, nbinsrhoi,0,nbinsrhoi);

    GetNormalisedCovarianceMatrix(TUnfResultCorXsec, ErrMatrixCorXsec, ErrMatrixNorm);
    GetNormalisedCovarianceMatrix(TUnfResultCorXsec_rebinnedA, EmatrixCorXsec_rebinnedA, EmatrixNorm_rebinnedA);
    GetNormalisedCovarianceMatrix(TUnfResultCorXsec_rebinnedB, EmatrixCorXsec_rebinnedB, EmatrixNorm_rebinnedB);


    fResultUnfHists->AddAt(ErrMatrixCorXsec, 13);
    fResultUnfHists->AddAt(EmatrixCorXsec_rebinnedA, 34);
    fResultUnfHists->AddAt(EmatrixCorXsec_rebinnedB, 37);

    // ***************************
    // Sytematics result -- ajeeta
    // ***************************

    if(fsetSyst){
      TObjArray *dummysysarray = (TObjArray*) GetAequiResSysHists()->Clone();
      //TObjArray *dummysysarray = (TObjArray*) GetResSysHists()->Clone();

      for(Int_t k=0; k<dummysysarray->GetEntries(); k++){

        //TH2D* dummyhisto=(TH2D*) ((*dummysysarray)[k])->Clone();
        std::string str, str2; 
        //str  = dummyhisto->GetName();
        str = ((*dummysysarray)[k])->GetName();
        size_t pos = str.find("_AequiResMatSys");
        str2 = str.substr(0,pos);
        std::cout<<"str2: "<<str2<<std::endl;

        TH1D* TUnfDeltaSys = (TH1D*) unfold.GetDeltaSysSource(Form("%s_Syserr", str2.c_str()),Form("%s_deltaSys",str2.c_str()),0,"ttbargen","*[UO]", 0);
        TH1D* TUnfDeltaSys_rebinnedA = (TH1D*) unfold.GetDeltaSysSource(Form("%s_Syserr", str2.c_str()),Form("%s_deltaSys_rebinnedA",str2.c_str()),0,"ttbargen","*[UO]", 0);
        if (frebinfine > 1) {
          MatrixUnf::rebinMultidimensionalInput(TUnfDeltaSys_rebinnedA,
                                              generatorBinning_rebinnedB->FindNode("ttbargen_rebinnedB"),
                                              generatorBinning->FindNode("ttbargen"),
                                              frebinfine);
        }
        TH1D* TUnfDeltaSys_rebinnedB = (TH1D*) unfold_rebinnedB.GetDeltaSysSource(Form("%s_Syserr", str2.c_str()),Form("%s_deltaSys_rebinnedB",str2.c_str()),0,"ttbargen_rebinnedB","*[UO]", 0);
        if(!TUnfDeltaSys_rebinnedB) {
          //in rare cases where the systematic variation is exactly zero, TUnfoldSys fails to add it to the list: https://root.cern.ch/doc/master/TUnfoldSys_8cxx_source.html#l00316
          //this is less likely to happen for _rebinnedA, so avoid the crash by taking that histo and resetting it to zero
          TUnfDeltaSys_rebinnedB = (TH1D*) unfold.GetDeltaSysSource(Form("%s_Syserr", str2.c_str()),Form("%s_deltaSys_rebinnedB",str2.c_str()),0,"ttbargen","*[UO]", 0);
          if (frebinfine > 1) {
            MatrixUnf::rebinMultidimensionalInput(TUnfDeltaSys_rebinnedB,
                                                generatorBinning_rebinnedB->FindNode("ttbargen_rebinnedB"),
                                                generatorBinning->FindNode("ttbargen"),
                                                frebinfine);
          }
          TUnfDeltaSys_rebinnedB->Reset();
        }

        TH2D* TUnfCovMatSys           = (TH2D*) ErrMatrix->Clone(Form("%s_CovMatSys",str2.c_str()));
        TH2D* TUnfCovMatSys_rebinnedA = (TH2D*) Ematrix_rebinnedA->Clone(Form("%s_CovMatSys_rebinnedA",str2.c_str()));
        TH2D* TUnfCovMatSys_rebinnedB = (TH2D*) Ematrix_rebinnedB->Clone(Form("%s_CovMatSys_rebinnedB",str2.c_str()));

        unfold.GetEmatrixSysSource(TUnfCovMatSys,Form("%s_Syserr",str2.c_str()));
        if(frebinfine>1){
         unfold.GetEmatrixSysSource(TUnfCovMatSys_rebinnedA,Form("%s_Syserr",str2.c_str()),binMap);
        }
        else unfold.GetEmatrixSysSource(TUnfCovMatSys_rebinnedA,Form("%s_Syserr",str2.c_str()));
        unfold_rebinnedB.GetEmatrixSysSource(TUnfCovMatSys_rebinnedB,Form("%s_Syserr",str2.c_str()));


        //apply corrections for absolute and normalised cross-sections

        GetAbsoluteCrossSec(TUnfDeltaSys, fBR, Lumi);
        GetAbsoluteCrossSec(TUnfDeltaSys_rebinnedA, fBR, Lumi);
        GetAbsoluteCrossSec(TUnfDeltaSys_rebinnedB, fBR, Lumi);

        fTUnfDeltaSys->AddAt(TUnfDeltaSys, k);
        fTUnfDeltaSys_rebinnedA->AddAt(TUnfDeltaSys_rebinnedA, k);
        fTUnfDeltaSys_rebinnedB->AddAt(TUnfDeltaSys_rebinnedB, k);
        //if(dbgUnf) cout<<"Added DeltaSys to the Obj array"<<std::endl;

        GetAbsoluteCovarianceMatrix(TUnfResultCorXsec, TUnfCovMatSys, fBR, Lumi);
        GetAbsoluteCovarianceMatrix(TUnfResultCorXsec_rebinnedA, TUnfCovMatSys_rebinnedA, fBR, Lumi);
        GetAbsoluteCovarianceMatrix(TUnfResultCorXsec_rebinnedB, TUnfCovMatSys_rebinnedB, fBR, Lumi);

        fTUnfCovMatSys->AddAt(TUnfCovMatSys, k);
        fTUnfCovMatSys_rebinnedA->AddAt(TUnfCovMatSys_rebinnedA, k);
        fTUnfCovMatSys_rebinnedB->AddAt(TUnfCovMatSys_rebinnedB, k);
        //if(dbgUnf) cout<<"Added covariance matrices to the Obj array"<<std::endl;

        TH2D* TUnfCovMatSysNorm =(TH2D*) ErrMatrix->Clone(Form("%s_CovMatSysNorm",str2.c_str()));
        TH2D* TUnfCovMatSysNorm_rebinnedA =(TH2D*) Ematrix_rebinnedA->Clone(Form("%s_CovMatSysNorm_rebinnedA",str2.c_str()));
        TH2D* TUnfCovMatSysNorm_rebinnedB =(TH2D*) Ematrix_rebinnedB->Clone(Form("%s_CovMatSysNorm_rebinnedB",str2.c_str()));

        GetNormalisedCovarianceMatrix(TUnfResultCorXsec, TUnfCovMatSys, TUnfCovMatSysNorm);
        GetNormalisedCovarianceMatrix(TUnfResultCorXsec_rebinnedA, TUnfCovMatSys_rebinnedA, TUnfCovMatSysNorm_rebinnedA);
        GetNormalisedCovarianceMatrix(TUnfResultCorXsec_rebinnedB, TUnfCovMatSys_rebinnedB, TUnfCovMatSysNorm_rebinnedB);

        fTUnfCovMatSysNorm->AddAt(TUnfCovMatSysNorm, k);
        fTUnfCovMatSysNorm_rebinnedA->AddAt(TUnfCovMatSysNorm_rebinnedA, k);
        fTUnfCovMatSysNorm_rebinnedB->AddAt(TUnfCovMatSysNorm_rebinnedB, k);


      }
      const int nsysfilled = dummysysarray->GetEntries();
      delete dummysysarray;

      // ************************
      // Background uncertainties
      // ************************

      vector<TMatrixD> TUnfCovMatBkgUncorr_rebinnedA_M;
      if(subtrBg) {
        for (int i_C = 0; i_C < nCombinedBkgHists; ++i_C)
        {

          //std::cout<<"i_C: "<<i_C<<" "<<nsysfilled<<std::endl;

          TH1D* TUnfDeltaBkg = (TH1D*) unfold.GetDeltaSysBackgroundScale(backgroundnames[i_C].Data(),Form("%sBkgSys_%s_deltaSys",fVariableName1.Data(),backgroundnames[i_C].Data()),0,"ttbargen","*[UO]",0);
          TH1D* TUnfDeltaBkg_rebinnedA = (TH1D*) unfold.GetDeltaSysBackgroundScale(backgroundnames[i_C].Data(),Form("%sBkgSys_%s_deltaSys_rebinnedA",fVariableName1.Data(),backgroundnames[i_C].Data()),0,"ttbargen","*[UO]",0);
          if (frebinfine > 1) {
            MatrixUnf::rebinMultidimensionalInput(TUnfDeltaBkg_rebinnedA,
                                                generatorBinning_rebinnedB->FindNode("ttbargen_rebinnedB"),
                                                generatorBinning->FindNode("ttbargen"),
                                                frebinfine);
          }
          TH1D* TUnfDeltaBkg_rebinnedB = (TH1D*) unfold_rebinnedB.GetDeltaSysBackgroundScale(backgroundnames[i_C].Data(),Form("%sBkgSys_%s_deltaSys_rebinnedB",fVariableName1.Data(),backgroundnames[i_C].Data()),0,"ttbargen_rebinnedB","*[UO]",0);


          TH2D* TUnfCovMatBkgSys           = (TH2D*) ErrMatrix->Clone(Form("%sBkgSys_%s_CovMatSys",fVariableName1.Data(),backgroundnames[i_C].Data()));
          TH2D* TUnfCovMatBkgSys_rebinnedA = (TH2D*) Ematrix_rebinnedA->Clone(Form("%sBkgSys_%s_CovMatSys_rebinnedA",fVariableName1.Data(),backgroundnames[i_C].Data()));
          TH2D* TUnfCovMatBkgSys_rebinnedB = (TH2D*) Ematrix_rebinnedB->Clone(Form("%sBkgSys_%s_CovMatSys_rebinnedB",fVariableName1.Data(),backgroundnames[i_C].Data()));

          unfold.GetEmatrixSysBackgroundScale(TUnfCovMatBkgSys,backgroundnames[i_C].Data());
          
          if(frebinfine>1){
           unfold.GetEmatrixSysBackgroundScale(TUnfCovMatBkgSys_rebinnedA,backgroundnames[i_C].Data(),binMap);
          }

          else unfold.GetEmatrixSysBackgroundScale(TUnfCovMatBkgSys_rebinnedA,backgroundnames[i_C].Data());
          unfold_rebinnedB.GetEmatrixSysBackgroundScale(TUnfCovMatBkgSys_rebinnedB,backgroundnames[i_C].Data());

          TH2D* TUnfCovMatBkgUncorr           = (TH2D*) unfold.GetEmatrixSysBackgroundUncorr(backgroundnames[i_C].Data(),Form("%sBkgSysUncorr_%s_CovMatSys",fVariableName1.Data(),backgroundnames[i_C].Data()),0,"ttbargen","*[UO]",1);
          TH2D* TUnfCovMatBkgUncorr_rebinnedA = (TH2D*) unfold.GetEmatrixSysBackgroundUncorr(backgroundnames[i_C].Data(),Form("%sBkgSysUncorr_%s_CovMatSys_rebinnedA",fVariableName1.Data(),backgroundnames[i_C].Data()),0,"ttbargen","*[UO]",1);
          if (frebinfine > 1) {
            MatrixUnf::rebinMultidimMigrationMatrix(TUnfCovMatBkgUncorr_rebinnedA,
                                                    generatorBinning_rebinnedB->FindNode("ttbargen_rebinnedB"),
                                                    generatorBinning_rebinnedB->FindNode("ttbargen_rebinnedB"),
                                                    generatorBinning->FindNode("ttbargen"),
                                                    generatorBinning->FindNode("ttbargen"),
                                                    frebinfine);
            // TUnfCovMatBkgUncorr_rebinnedA->Rebin2D(frebinfine, frebinfine);
          }
          TH2D* TUnfCovMatBkgUncorr_rebinnedB = (TH2D*) unfold_rebinnedB.GetEmatrixSysBackgroundUncorr(backgroundnames[i_C].Data(),Form("%sBkgSysUncorr_%s_CovMatSys_rebinnedB",fVariableName1.Data(),backgroundnames[i_C].Data()),0,"ttbargen_rebinnedB","*[UO]",1);

          // Apply corrections for absolute and normalised cross-sections
          GetAbsoluteCrossSec(TUnfDeltaBkg, fBR, Lumi);
          GetAbsoluteCrossSec(TUnfDeltaBkg_rebinnedA, fBR, Lumi);
          GetAbsoluteCrossSec(TUnfDeltaBkg_rebinnedB, fBR, Lumi);

          fTUnfDeltaSys->AddAt(TUnfDeltaBkg, i_C+nsysfilled);
          fTUnfDeltaSys_rebinnedA->AddAt(TUnfDeltaBkg_rebinnedA, i_C+nsysfilled);
          fTUnfDeltaSys_rebinnedB->AddAt(TUnfDeltaBkg_rebinnedB, i_C+nsysfilled);

          GetAbsoluteCovarianceMatrix(TUnfResultCorXsec, TUnfCovMatBkgSys, fBR, Lumi);
          GetAbsoluteCovarianceMatrix(TUnfResultCorXsec_rebinnedA, TUnfCovMatBkgSys_rebinnedA, fBR, Lumi);
          GetAbsoluteCovarianceMatrix(TUnfResultCorXsec_rebinnedB, TUnfCovMatBkgSys_rebinnedB, fBR, Lumi);

          fTUnfCovMatSys->AddAt(TUnfCovMatBkgSys, i_C+nsysfilled);
          fTUnfCovMatSys_rebinnedA->AddAt(TUnfCovMatBkgSys_rebinnedA, i_C+nsysfilled);
          fTUnfCovMatSys_rebinnedB->AddAt(TUnfCovMatBkgSys_rebinnedB, i_C+nsysfilled);

          TH2D* TUnfCovMatBkgSysNorm =(TH2D*) ErrMatrix->Clone(Form("%sBkgSys_%s_CovMatSysNorm",fVariableName1.Data(),backgroundnames[i_C].Data()));
          TH2D* TUnfCovMatBkgSysNorm_rebinnedA =(TH2D*) Ematrix_rebinnedA->Clone(Form("%sBkgSys_%s_CovMatSysNorm_rebinnedA",fVariableName1.Data(),backgroundnames[i_C].Data()));
          TH2D* TUnfCovMatBkgSysNorm_rebinnedB =(TH2D*) Ematrix_rebinnedB->Clone(Form("%sBkgSys_%s_CovMatSysNorm_rebinnedB",fVariableName1.Data(),backgroundnames[i_C].Data()));

          GetNormalisedCovarianceMatrix(TUnfResultCorXsec, TUnfCovMatBkgSys, TUnfCovMatBkgSysNorm);
          GetNormalisedCovarianceMatrix(TUnfResultCorXsec_rebinnedA, TUnfCovMatBkgSys_rebinnedA, TUnfCovMatBkgSysNorm_rebinnedA);
          GetNormalisedCovarianceMatrix(TUnfResultCorXsec_rebinnedB, TUnfCovMatBkgSys_rebinnedB, TUnfCovMatBkgSysNorm_rebinnedB);

          fTUnfCovMatSysNorm->AddAt(TUnfCovMatBkgSysNorm, i_C+nsysfilled);
          fTUnfCovMatSysNorm_rebinnedA->AddAt(TUnfCovMatBkgSysNorm_rebinnedA, i_C+nsysfilled);
          fTUnfCovMatSysNorm_rebinnedB->AddAt(TUnfCovMatBkgSysNorm_rebinnedB, i_C+nsysfilled);


          GetAbsoluteCovarianceMatrix(TUnfResultCorXsec, TUnfCovMatBkgUncorr, fBR, Lumi);
          GetAbsoluteCovarianceMatrix(TUnfResultCorXsec_rebinnedA, TUnfCovMatBkgUncorr_rebinnedA, fBR, Lumi);
          GetAbsoluteCovarianceMatrix(TUnfResultCorXsec_rebinnedB, TUnfCovMatBkgUncorr_rebinnedB, fBR, Lumi);

          fTUnfCovMatSys->AddAt(TUnfCovMatBkgUncorr, i_C+nCombinedBkgHists+nsysfilled);
          fTUnfCovMatSys_rebinnedA->AddAt(TUnfCovMatBkgUncorr_rebinnedA, i_C+nCombinedBkgHists+nsysfilled);
          fTUnfCovMatSys_rebinnedB->AddAt(TUnfCovMatBkgUncorr_rebinnedB, i_C+nCombinedBkgHists+nsysfilled);

          TH2D* TUnfCovMatBkgUncorrNorm           = (TH2D*) ErrMatrix->Clone(Form("%sBkgSysUncorr_%s_CovMatSysNorm",fVariableName1.Data(),backgroundnames[i_C].Data()));
          TH2D* TUnfCovMatBkgUncorrNorm_rebinnedA = (TH2D*) Ematrix_rebinnedA->Clone(Form("%sBkgSysUncorr_%s_CovMatSysNorm_rebinnedA",fVariableName1.Data(),backgroundnames[i_C].Data()));
          TH2D* TUnfCovMatBkgUncorrNorm_rebinnedB = (TH2D*) Ematrix_rebinnedB->Clone(Form("%sBkgSysUncorr_%s_CovMatSysNorm_rebinnedB",fVariableName1.Data(),backgroundnames[i_C].Data()));

          GetNormalisedCovarianceMatrix(TUnfResultCorXsec, TUnfCovMatBkgUncorr, TUnfCovMatBkgUncorrNorm);
          GetNormalisedCovarianceMatrix(TUnfResultCorXsec_rebinnedA, TUnfCovMatBkgUncorr_rebinnedA, TUnfCovMatBkgUncorrNorm_rebinnedA);
          GetNormalisedCovarianceMatrix(TUnfResultCorXsec_rebinnedB, TUnfCovMatBkgUncorr_rebinnedB, TUnfCovMatBkgUncorrNorm_rebinnedB);

          fTUnfCovMatSysNorm->AddAt(TUnfCovMatBkgUncorrNorm, i_C+nCombinedBkgHists+nsysfilled);
          fTUnfCovMatSysNorm_rebinnedA->AddAt(TUnfCovMatBkgUncorrNorm_rebinnedA, i_C+nCombinedBkgHists+nsysfilled);
          fTUnfCovMatSysNorm_rebinnedB->AddAt(TUnfCovMatBkgUncorrNorm_rebinnedB, i_C+nCombinedBkgHists+nsysfilled);


          TUnfCovMatBkgUncorr_rebinnedA_M.push_back(TH2toTMatrixD(TUnfCovMatBkgUncorr_rebinnedA));


        }
      }

      // Response matrix MC stat uncertainties
      TH2D* TUnfCovMatSysUncorr           = (TH2D*) unfold.GetEmatrixSysUncorr(Form("%sRespMatSys_Uncorr_CovMatSys",fVariableName1.Data()),0,"ttbargen","*[UO]",1);
      TH2D* TUnfCovMatSysUncorr_rebinnedA = (TH2D*) unfold.GetEmatrixSysUncorr(Form("%sRespMatSys_Uncorr_CovMatSys_rebinnedA",fVariableName1.Data()),0,"ttbargen","*[UO]",1);
      
      if (frebinfine > 1) {
        MatrixUnf::rebinMultidimMigrationMatrix(TUnfCovMatSysUncorr_rebinnedA,
                                                generatorBinning_rebinnedB->FindNode("ttbargen_rebinnedB"),
                                                generatorBinning_rebinnedB->FindNode("ttbargen_rebinnedB"),
                                                generatorBinning->FindNode("ttbargen"),
                                                generatorBinning->FindNode("ttbargen"),
                                                frebinfine);
      }
      TH2D* TUnfCovMatSysUncorr_rebinnedB = (TH2D*) unfold_rebinnedB.GetEmatrixSysUncorr(Form("%sRespMatSys_Uncorr_CovMatSys_rebinnedB",fVariableName1.Data()),0,"ttbargen_rebinnedB","*[UO]",1);

      // Apply corrections for absolute and normalised cross-sections
      GetAbsoluteCovarianceMatrix(TUnfResultCorXsec, TUnfCovMatSysUncorr, fBR, Lumi);
      GetAbsoluteCovarianceMatrix(TUnfResultCorXsec_rebinnedA, TUnfCovMatSysUncorr_rebinnedA, fBR, Lumi);
      GetAbsoluteCovarianceMatrix(TUnfResultCorXsec_rebinnedB, TUnfCovMatSysUncorr_rebinnedB, fBR, Lumi);

      fTUnfCovMatSys->AddAt(TUnfCovMatSysUncorr, 2*nCombinedBkgHists+nsysfilled);
      fTUnfCovMatSys_rebinnedA->AddAt(TUnfCovMatSysUncorr_rebinnedA, 2*nCombinedBkgHists+nsysfilled);
      fTUnfCovMatSys_rebinnedB->AddAt(TUnfCovMatSysUncorr_rebinnedB, 2*nCombinedBkgHists+nsysfilled);

      TH2D* TUnfCovMatSysUncorrNorm =(TH2D*) ErrMatrix->Clone(Form("%sRespMatSys_Uncorr_CovMatSysNorm",fVariableName1.Data()));
      TH2D* TUnfCovMatSysUncorrNorm_rebinnedA =(TH2D*) Ematrix_rebinnedA->Clone(Form("%sRespMatSys_Uncorr_CovMatSysNorm_rebinnedA",fVariableName1.Data()));
      TH2D* TUnfCovMatSysUncorrNorm_rebinnedB =(TH2D*) Ematrix_rebinnedB->Clone(Form("%sRespMatSys_Uncorr_CovMatSysNorm_rebinnedB",fVariableName1.Data()));

      GetNormalisedCovarianceMatrix(TUnfResultCorXsec, TUnfCovMatSysUncorr, TUnfCovMatSysUncorrNorm);
      GetNormalisedCovarianceMatrix(TUnfResultCorXsec_rebinnedA, TUnfCovMatSysUncorr_rebinnedA, TUnfCovMatSysUncorrNorm_rebinnedA);
      GetNormalisedCovarianceMatrix(TUnfResultCorXsec_rebinnedB, TUnfCovMatSysUncorr_rebinnedB, TUnfCovMatSysUncorrNorm_rebinnedB);

      fTUnfCovMatSysNorm->AddAt(TUnfCovMatSysUncorrNorm, 2*nCombinedBkgHists+nsysfilled);
      fTUnfCovMatSysNorm_rebinnedA->AddAt(TUnfCovMatSysUncorrNorm_rebinnedA, 2*nCombinedBkgHists+nsysfilled);
      fTUnfCovMatSysNorm_rebinnedB->AddAt(TUnfCovMatSysUncorrNorm_rebinnedB, 2*nCombinedBkgHists+nsysfilled);

    }
  }

  
  if(dbgUnf) printf("MatrixUnf::prepRunUnfolding: Finished! Adding results to ObjArray\n"); 
  
  // clean up some stuff
  if(ensTest) {
    delete Rhoi_rebinnedA;
    delete Rhoi_rebinnedB;

    delete ErrMatrix;
    delete Ematrix_rebinnedA;
    delete Ematrix_rebinnedB;

    delete this->GetAequiMeas();
    delete GenVecHist_rebinnedB;
    delete AequiM_rebinnedB;

    delete ResMat;
    delete ResMat_rebinnedB;
  }

  return kTRUE;
}

/**
 * @brief Creates a mapping of bins between two N dimensional histograms where one has been binned more finely than the
 *        other and both have been unraveled to be 1 dimensional histograms.
 * 
 * @param binning The TUnfoldBinning object that contains the information of how each axis is binned
 * @param n_sub_bins The division of binning between the two histograms. Assumed every axis is binned equally fine.
 * @param binMap The mapping of the bins between the histograms.
 */
void MatrixUnf::CreateBinMap(const TUnfoldBinning* binning, int n_sub_bins, int binMap[]) {
  // Grab the number of fine bins from binning object and total number of bins
  std::vector<int> numFineBins;
  int totalFineBins = 1;
  for (int i = 0; i < binning->GetDistributionDimension(); ++i) {
    const TVectorD* vec_binEdges = binning->GetDistributionBinning(i);
    numFineBins.push_back(vec_binEdges->GetNoElements() - 1);
    totalFineBins *= numFineBins[i];
  }

  // underflow and overflow map to underflow in new histogram
  binMap[0] = 0;
  binMap[totalFineBins + 1] = 0;

  // create the bin mapping
  for (int i = 0; i < totalFineBins; ++i) {
    // the index gets contributions from where you are at along each of the axes
    // the axes are unraveled though so need to do some fancy flooring to get correct mapping
    int index  = 0;
    for (int j = 0; j < numFineBins.size(); ++j) {
      // Calculate how many non-fine bins you needed to skip over to get to this fine bin index
      // for example, if you at a fine bin index corresponding to 0th bin on axis 0 and 1st bin on axis 1
      //      for non-fine binning then you had to first skip over numFineBins[0] bins to get to the 1st bin on axis 1
      int binOffset = 1;
      for (int k = 0; k < j; ++k) {
        binOffset *= numFineBins[k];
      }

      // Find the bin index independent of the other axes
      int binRepetitionAmount = 1;
      for (int k = 0; k <= j; ++k) {
        binRepetitionAmount *= numFineBins[k];
      }
      int jthAxisBinNumber = (i % binRepetitionAmount) / binOffset;

      index += (binOffset / pow(n_sub_bins, j)) * floor(floor(jthAxisBinNumber) / n_sub_bins);
    }
    binMap[i + 1] = index + 1;
  }
}

/**
 * @brief Compute the inclusive forward-backward asymmetry in the distribution
 * 
 * @param h The histogram/distribution to compute the forward-backward asymmetry
 * @param afb The forward-backward asymmetry
 * @param afberr The uncertainty in forward-backward asymmetry taking into account correlations
 */
void MatrixUnf::GetAfb(TH1D* h, Double_t& afb, Double_t& afberr) {
  Int_t nbins = h->GetNbinsX();
  Double_t event_minus;
  Double_t event_plus;
  Double_t event_total;
  Double_t event_plus_err;
  Double_t event_minus_err;

  // event_minus  = h-> IntegralAndError(0, nbins/2, event_plus_err,"");
  event_minus = h->IntegralAndError(1, nbins / 2, event_minus_err, "");
  // event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_minus_err,"");
  event_plus = h->IntegralAndError(nbins / 2 + 1, nbins, event_plus_err, "");
  event_total = event_plus + event_minus;

  // std::cout<<event_minus<<" "<<event_minus_err<<" "<<event_plus<<" "<<event_plus_err<<" "<<event_total<<std::endl;

  afb = (event_plus - event_minus) / (event_plus + event_minus);
  afberr = sqrt(4 *
                (event_plus * event_plus * event_minus_err * event_minus_err +
                 event_minus * event_minus * event_plus_err * event_plus_err) /
                (event_total * event_total * event_total * event_total));
}

/**
 * @brief Compute the inclusive forward-backward asymmetry of a distribution 
 * while taking into account the correlations between sources of uncertainty
 * 
 * @param histogram The histogram/distribution to compute the forward-backward asymmetry
 * @param covarianceM Covariance matrix for the uncertainties on the bin contents of the distribution
 * @param afb The forward-backward asymmetries
 * @param afberr The uncertainties in forward-backward asymmetry with taking into account correlations
 * @param numAsymmetries The number of asymmetries in the distribution that need to be computed 
 */
void MatrixUnf::GetAfbWithCorrelations(TH1* histogram, TMatrixD& covarianceM, Double_t* afb, Double_t* afberr, int numAsymmetries) {
  // Get histogram info
  const int nbins = histogram->GetNbinsX();
  const int numBinsPerAsymmetry = nbins / numAsymmetries;
  double n[nbins];
  for (int i = 0; i < nbins; i++) {
    n[i] = histogram->GetBinContent(i + 1);
  }

  // Setup Alpha Vector
  double alpha[numAsymmetries][numBinsPerAsymmetry];
  for (int j = 0; j < numAsymmetries; j++) {
    for (int i = 0; i < numBinsPerAsymmetry; i++) {
      if (i < numBinsPerAsymmetry / 2) {
        alpha[j][i] = -1;
      } else {
        alpha[j][i] = 1;
      }
    }
  }

  // Components of the error calculation
  double sum_n[numAsymmetries];
  double sum_alpha_n[numAsymmetries];
  for (int j = 0; j < numAsymmetries; j++) {
    int offset = numBinsPerAsymmetry * j;
    for (int i = 0; i < numBinsPerAsymmetry; i++) {
      // initialize to 0
      if (i == 0) {
        sum_n[j] = 0.0;
        sum_alpha_n[j] = 0.0;
      }

      sum_n[j] += n[offset + i];
      sum_alpha_n[j] += alpha[j][i] * n[offset + i];
    }
  }

  // partial derivative of the asymmetry with respect to that bin's contents
  double dfdn[nbins];
  for (int j = 0; j < numAsymmetries; j++) {
    int offset = numBinsPerAsymmetry * j;
    for (int i = 0; i < numBinsPerAsymmetry; i++) {
      dfdn[offset + i] = (alpha[j][i] * sum_n[j] - sum_alpha_n[j]) / pow(sum_n[j], 2);
    }
  }

  // Error Calculation
  for (int k = 0; k < numAsymmetries; k++) {
    int offset = numBinsPerAsymmetry * k;
    afberr[k] = 0.0;
    for (int i = 0; i < numBinsPerAsymmetry; i++) {
      for (int j = 0; j < numBinsPerAsymmetry; j++) {
        afberr[k] += covarianceM(offset + i, offset + j) * dfdn[offset + i] * dfdn[offset + j];
      }
    }
    afberr[k] = sqrt(afberr[k]);
  }

  // Calculate Afb
  for (int i = 0; i < numAsymmetries; i++) {
    afb[i] = sum_alpha_n[i] / sum_n[i];
  }
}


/**
 * @brief Compute the per-bin forward-backward asymmetry of a distribution. For example the first asymmetry
 * is between the first and last bin.
 * 
 * @param histogram The histogram/distribution to compute the forward-backward asymmetry
 * @param covarianceM Covariance matrix for the uncertainties on the bin contents of the distribution
 * @param myafb The forward-backward asymmetries
 * @param myerr The uncertainties in forward-backward asymmetry with taking into account correlations
 * @param AFBcovarianceM The covariance matrix for the asymmetries of the distribution. Essentially the corrected
 *   covariance matrix by the Jacobian for linear error propagation.
 * @param numAsymmetries The number of asymmetries in the distribution that need to be computed 
 */
void MatrixUnf::GetAfbWithCorrelationsBinByBin(TH1* histogram, TMatrixD& covarianceM, std::vector<double>& myafb, std::vector<double>& myerr, TMatrixD& AFBcovarianceM, int numAsymmetries) {
  myafb.clear();
  myerr.clear();

  // Get histogram info
  const int nbins = histogram->GetNbinsX();
  const int numBinsPerAsymmetry = nbins / numAsymmetries;

  Double_t afbbin[numAsymmetries][numBinsPerAsymmetry];
  Double_t afberrbin[numAsymmetries][numBinsPerAsymmetry];

  double n[numAsymmetries][numBinsPerAsymmetry];
  for (int j = 0; j < numAsymmetries; j++) {
    int offset = numBinsPerAsymmetry * j;
    for (int i = 0; i < numBinsPerAsymmetry; i++) {
      n[j][i] = histogram->GetBinContent(offset + i + 1);
    }
  }

  // Setup Alpha Vector
  double alpha[numAsymmetries][numBinsPerAsymmetry];
  for (int j = 0; j < numAsymmetries; j++) {
    for (int i = 0; i < numBinsPerAsymmetry; i++) {
      if (i < numBinsPerAsymmetry / 2) {
        alpha[j][i] = -1;
      } else {
        alpha[j][i] = 1;
      }
    }
  }

  // Components of the error calculation
  double sum_n[numAsymmetries][numBinsPerAsymmetry / 2];
  double sum_alpha_n[numAsymmetries][numBinsPerAsymmetry / 2];
  for (int j = 0; j < numAsymmetries; j++) {
    for (int i = 0; i < numBinsPerAsymmetry / 2; i++) {
      sum_n[j][i] = n[j][i] + n[j][numBinsPerAsymmetry - 1 - i];
      sum_alpha_n[j][i] = alpha[j][i] * n[j][i] + alpha[j][numBinsPerAsymmetry - 1 - i] * n[j][numBinsPerAsymmetry - 1 - i];
    }
  }

  double dfdn[numAsymmetries][numBinsPerAsymmetry];
  for (int j = 0; j < numAsymmetries; j++) {
    for (int i = 0; i < numBinsPerAsymmetry; i++) {
      int k;
      if (i < numBinsPerAsymmetry / 2)
        k = i;
      else
        k = numBinsPerAsymmetry - 1 - i;
      dfdn[j][i] = (alpha[j][i] * sum_n[j][k] - sum_alpha_n[j][k]) / pow(sum_n[j][k], 2);
    }
  }

  // Error Calculation
  for (int i = 0; i < nbins / 2; i++) {
    for (int j = 0; j < nbins / 2; j++) {
      AFBcovarianceM(i, j) = 0;
    }
  }

  for (int m = 0; m < numAsymmetries; m++) {
    int offset = numBinsPerAsymmetry * m;
    int AfbOffset = numBinsPerAsymmetry * m / 2;
    for (int k = 0; k < numBinsPerAsymmetry / 2; k++) {
      afberrbin[m][k] = 0.;
      for (int l = 0; l < numBinsPerAsymmetry / 2; l++) {
        for (int i = 0; i < numBinsPerAsymmetry; i++) {
          for (int j = 0; j < numBinsPerAsymmetry; j++) {
            if ((i == k || i == numBinsPerAsymmetry - 1 - k) && (j == l || j == numBinsPerAsymmetry - 1 - l)) {
              if (l == k)
                afberrbin[m][k] += covarianceM(offset + i, offset + j) * dfdn[m][i] * dfdn[m][j];
              AFBcovarianceM(AfbOffset + k, AfbOffset + l) += covarianceM(offset + i, offset + j) * dfdn[m][i] * dfdn[m][j];
            }
          }
        }
      }
      afberrbin[m][k] = sqrt(afberrbin[m][k]);
      afbbin[m][k] = sum_alpha_n[m][k] / sum_n[m][k];

      myafb.push_back(afbbin[m][k]);
      myerr.push_back(afberrbin[m][k]);
    }
  }
  // double afb = sum_alpha_n_total / sum_n_total;
  // std::cout<<"AFB = "<<afb<<std::endl;

  /*
      //check that we get the same result for the error propagation using matrix multiplication
      cout<<"AFBcovarianceM: "<<endl;
      AFBcovarianceM.Print("f=%1.5g ");

      TMatrixD m_dfdn(nbins2, nbins);
      for(int k=0;k<nbins2;k++){
        for(int i=0;i<nbins;i++){
          if( (i==k || i==nbins-1-k ) ) m_dfdn(k,i) = ( alpha[i] * sum_n[k] - sum_alpha_n[k] ) / pow(sum_n[k],2);
          else m_dfdn(k,i) = 0;
        }
      }
      TMatrixD m_dfdnT = m_dfdn;
      m_dfdnT = m_dfdnT.T();

      TMatrixD m_covf = m_dfdn*covarianceM*m_dfdnT;

      cout<<"Matrix version: "<<endl;
      m_covf.Print("f=%1.5g ");
    */
}



void MatrixUnf::GetAfbCorrMAllVars(TH1D* histAllVars, TMatrixD &corrM, std::vector<double> &myafb, std::vector<double> &myerr, TMatrixD &AFBcorrM, bool isCovInput){

  myafb.clear();
  myerr.clear();

  // Amandeep : Hardcoded here
  // const int nbinvar = 24;
  const int nbinvar = 6;
  const int nbins   = histAllVars->GetNbinsX();
  const int nvar    = nbins/nbinvar;

  double n[nbins];
  double ne[nbins];
  
  for(int i=0; i < nbins; i++){
    n[i]  = histAllVars->GetBinContent(i+1);
    ne[i] = histAllVars->GetBinError(i+1);
  }

  // Setup Alpha Vector
  double alpha[nbinvar];
  for(int i=0;i<nbinvar;i++) if(i < nbinvar/2 ){ alpha[i] = -1;}else{ alpha[i] = 1;}

  // Components of the error calculation
  double sum_n[nvar]={0};
  double sum_alpha_n[nvar]={0};

  for(int k=0;k<nvar;k++){
    for(int i=0;i<nbinvar;i++){
      sum_n[k] += n[i+k*nbinvar];
      sum_alpha_n[k] += alpha[i] * n[i+k*nbinvar];
    }
  }

  TMatrixD covarianceM(nbins, nbins);
  for(int j=0;j<nbins;j++){
    for(int i=0;i<nbins;i++){
      if(isCovInput) covarianceM(j,i) = corrM(j,i);
      else covarianceM(j,i) = corrM(j,i) * ne[j] * ne[i];
    }
  }

  TMatrixD m_dfdn(nvar, nbins);
  for(int k=0;k<nvar;k++){
    for(int i=0;i<nbins;i++){
      if( i/nbinvar == k ) m_dfdn(k,i) = ( alpha[i%nbinvar] * sum_n[k] - sum_alpha_n[k] ) / pow(sum_n[k],2);
      else m_dfdn(k,i) = 0;
    }
  }
  TMatrixD m_dfdnT = m_dfdn;
  m_dfdnT = m_dfdnT.T();

  AFBcorrM = m_dfdn*covarianceM*m_dfdnT;

  for(int k=0;k<nvar;k++){
    myafb.push_back(sum_alpha_n[k] / sum_n[k]);
    myerr.push_back( sqrt(AFBcorrM(k,k) ) );
    //cout<<"AFBV"<<k<<": "<<myafb[k]<<" "<<myerr[k]<<endl;
  }

  if(!isCovInput) AFBcorrM.NormByDiag(TMatrixDDiag(AFBcorrM));
  //cout<<"Correlation matrix for AFB of all variables: "<<endl;
  //AFBcorrM.Print("f=%1.5g ");

}


void MatrixUnf::GetAfbBinsCorrMAllVars(TH1D* histAllVars, TMatrixD &corrM, std::vector<std::vector<double>> &afbvecvec, std::vector<std::vector<double>> &errvecvec, TMatrixD &AFBBinscorrM, TH1D* hist_acombfactor[], TMatrixD &CcorrM, bool isCovInput, int numCoefficients){
  afbvecvec.clear();
  errvecvec.clear();

  // Get histAllVars info
  // Amandeep : Hardcoded here
  // const int nbinvar = 24;
  const int nbinvar = 6;
  const int nbins   = histAllVars->GetNbinsX();
  const int nvar    = nbins/nbinvar;
  const int nbins2  = nbinvar/2;

  double n[nbins];
  double ne[nbins];
  for(int i=0;i<nbins;i++){
    n[i]  = histAllVars->GetBinContent(i+1);
    ne[i] = histAllVars->GetBinError(i+1);
  }

  //Setup Alpha Vector
  double alpha[nbinvar];
  for(int i=0;i<nbinvar;i++) if(i < nbinvar/2 ){ alpha[i] = -1;}else{ alpha[i] = 1;}

  //Components of the error calculation
  double sum_n[nvar][nbins2]={0};
  double sum_alpha_n[nvar][nbins2]={0};
  double sum_n_total[nvar]={0};
  //double sum_alpha_n_total[nvar]={0};

  for(int k=0;k<nvar;k++){
    for(int i=0;i<nbins2;i++){
      sum_n[k][i] = n[i + k*nbinvar] + n[nbinvar-1-i + k*nbinvar];
      sum_alpha_n[k][i] = alpha[i] * n[i + k*nbinvar] + alpha[nbinvar-1-i] * n[nbinvar-1-i + k*nbinvar];
      sum_n_total[k] += sum_n[k][i];
      //sum_alpha_n_total[k] += sum_alpha_n[k][i];
    }
  }

  //cout<<"Creating covariance matrix"<<endl;
  TMatrixD covarianceM(nbins, nbins);
  for(int j=0;j<nbins;j++){
    for(int i=0;i<nbins;i++){
      if(isCovInput) covarianceM(j,i) = corrM(j,i);
      else covarianceM(j,i) = corrM(j,i) * ne[j] * ne[i];
    }
  }

  //cout<<"Creating dfdn vector"<<endl;
  TMatrixD m_dfdn(nvar*nbins2, nbins);
  for(int k=0;k<nvar*nbins2;k++){
    for(int i=0;i<nbins;i++){
      if( i/nbinvar == k/nbins2 && (k%nbins2 == i%nbinvar || k%nbins2 == nbinvar-1-i%nbinvar) ) m_dfdn(k,i) = ( alpha[i%nbinvar] * sum_n[k/nbins2][k%nbins2] - sum_alpha_n[k/nbins2][k%nbins2] ) / pow(sum_n[k/nbins2][k%nbins2],2);
      else m_dfdn(k,i) = 0;
    }
  }
  //cout<<"Creating transpoes of dfdn vector"<<endl;
  TMatrixD m_dfdnT = m_dfdn;
  m_dfdnT = m_dfdnT.T();

  //cout<<"Calculate AFB bins correlation matrix"<<endl;
  AFBBinscorrM = m_dfdn*covarianceM*m_dfdnT;

  for(int k=0;k<nvar;k++){
    std::vector<double> myafb;
    std::vector<double> myerr;
    for(int i=0;i<nbins2;i++){
      myafb.push_back(sum_alpha_n[k][i] / sum_n[k][i]);
      myerr.push_back( sqrt(AFBBinscorrM(k*nbins2+i,k*nbins2+i) ) );
      //cout<<"AFB"<<i<<"V"<<k<<": "<<myafb[i]<<" "<<myerr[i]<<endl;
    }
    afbvecvec.push_back(myafb);
    errvecvec.push_back(myerr);
  }


  TMatrixD m_dCdab(nvar,nvar*nbins2);
  for(int k=0;k<nvar;k++){
    for(int i=0;i<nvar*nbins2;i++){
      if( i/nbins2 == k ) m_dCdab(k,i) = hist_acombfactor[i % (nbins2 * numCoefficients)]->GetBinContent(k / numCoefficients + 1);
      else m_dCdab(k,i) = 0;
    }
  }
  TMatrixD m_dCdabT = m_dCdab;
  m_dCdabT = m_dCdabT.T();

  CcorrM = m_dCdab*AFBBinscorrM*m_dCdabT;

  for(int k=0;k<nvar;k++){
    //cout<<"CerrV"<<k<<": "<<sqrt(CcorrM(k,k))<<endl;
  }

  if(!isCovInput) AFBBinscorrM.NormByDiag(TMatrixDDiag(AFBBinscorrM));
  //cout<<"Correlation matrix for AFB bins of all variables: "<<endl;
  //AFBBinscorrM.Print("f=%1.5g ");
  if(!isCovInput) CcorrM.NormByDiag(TMatrixDDiag(CcorrM));
  //cout<<"Correlation matrix for extracted coefficient from all variables, calculated from AFBBinscorrM: "<<endl;
  //CcorrM.Print("f=%1.5g ");

}



void MatrixUnf::GetBinSumCorrMAllVars(TH1D* histAllVars, TMatrixD &corrM, std::vector<double> &binsum, std::vector<double> &binsumerr, TMatrixD &BinSumcorrM, bool isCovInput){

  binsum.clear();
  binsumerr.clear();

  // Get histAllVars info
  // Amandeep : Hardcoded here
  // const int nbinvar = 24;
  const int nbinvar = 6;
  const int nbins   = histAllVars->GetNbinsX();
  const int nvar    = nbins/nbinvar;

  double n[nbins];
  double ne[nbins];
  for(int i=0;i<nbins;i++){
    n[i] = histAllVars->GetBinContent(i+1);
    ne[i] = histAllVars->GetBinError(i+1);
  }

  double sum_n[nvar]={0};

  for(int k=0;k<nvar;k++){
    for(int i=0;i<nbinvar;i++){
      sum_n[k] += n[i+k*nbinvar];
    }
  }

  TMatrixD covarianceM(nbins, nbins);
  for(int j=0;j<nbins;j++){
    for(int i=0;i<nbins;i++){
      if(isCovInput) covarianceM(j,i) = corrM(j,i);
      else covarianceM(j,i) = corrM(j,i) * ne[j] * ne[i];
    }
  }

  TMatrixD m_dfdn(nvar, nbins);
  for(int k=0;k<nvar;k++){
    for(int i=0;i<nbins;i++){
      if( i/nbinvar == k ) m_dfdn(k,i) = 1.;
      else m_dfdn(k,i) = 0;
    }
  }
  TMatrixD m_dfdnT = m_dfdn;
  m_dfdnT = m_dfdnT.T();

  BinSumcorrM = m_dfdn*covarianceM*m_dfdnT;

  for(int k=0;k<nvar;k++){
    binsum.push_back(sum_n[k]);
    binsumerr.push_back( sqrt( BinSumcorrM(k,k)>0?BinSumcorrM(k,k):0 ) );
    //cout<<"BinSumV"<<k<<": "<<binsum[k]<<" "<<binsumerr[k]<<endl;
  }

  if(!isCovInput) BinSumcorrM.NormByDiag(TMatrixDDiag(BinSumcorrM));
  //cout<<"Correlation matrix for AFB of all variables: "<<endl;
  //BinSumcorrM.Print("f=%1.5g ");

}

/**
 * @brief Compute the per-bin forward-backward asymmetry of a distribution. For example the first asymmetry
 * is between the first and last bin.
 * 
 * @param histogram The histogram/distribution to compute the forward-backward asymmetry
 * @param myafb The forward-backward asymmetries
 * @param myerr The uncertainties in forward-backward asymmetry with taking into account correlations
 * @param numCoefficients The number of coefficients in the distribution that need to be computed/extracted 
 *   (via forward-backward asymmetries) 
 */
void MatrixUnf::GetAfbBinByBin(TH1D* histogram, std::vector<double>& myafb, std::vector<double>& myerr, int numCoefficients) {
  myafb.clear();
  myerr.clear();

  // Get histogram info
  const int nbins = histogram->GetNbinsX();
  const int numBinsPerCoefficient = nbins / numCoefficients;
  const int numAsymmetries = numBinsPerCoefficient / 2;
  for (int j = 0; j < numCoefficients; j++) {
    int offset = j * numBinsPerCoefficient;
    Double_t afb;
    Double_t afberr;

    Double_t event_minus;
    Double_t event_plus;
    Double_t event_minus_err;
    Double_t event_plus_err;
    Double_t event_total;

    for (int i = 0; i < numAsymmetries; i++) {
      event_minus = histogram->GetBinContent(offset + i + 1);
      event_plus = histogram->GetBinContent(offset + numBinsPerCoefficient - i);
      event_minus_err = histogram->GetBinError(offset + i + 1);
      event_plus_err = histogram->GetBinError(offset + numBinsPerCoefficient - i);
      event_total = event_plus + event_minus;

      afb = (event_plus - event_minus) / (event_plus + event_minus);
      afberr = sqrt(4 *
                    (event_plus * event_plus * event_minus_err * event_minus_err +
                    event_minus * event_minus * event_plus_err * event_plus_err) /
                    (event_total * event_total * event_total * event_total));

      myafb.push_back(afb);
      myerr.push_back(afberr);
    }
  }
}

void MatrixUnf::GetAbsoluteCrossSec(TH1D*& uncorrHist, Double_t BR, Double_t Lumi) {
  // correct the result for BR, Lumi & store it
  for(Int_t i =1; i <=uncorrHist->GetNbinsX(); i++)	{
    Double_t cF   = (BR*Lumi);
    //Double_t cF = (BR*Lumi*uncorrHist->GetBinWidth(i));
    uncorrHist->SetBinContent(i, uncorrHist->GetBinContent(i)/cF);
    uncorrHist->SetBinError(i, uncorrHist->GetBinError(i)/cF);
  }
} 


// Amandeep : Original
void MatrixUnf::GetBinWidthCorrectedCrossSec(TH1D*& uncorrHist) {
  // correct the result for bin width & store it
  for(Int_t i = 1; i <= uncorrHist->GetNbinsX(); i++) {
    Double_t cF = uncorrHist->GetBinWidth(i);
    uncorrHist->SetBinContent(i, uncorrHist->GetBinContent(i)/cF);
    uncorrHist->SetBinError(i, uncorrHist->GetBinError(i)/cF);
  }
} 


// Amandeep : Modifying this for N-D cross-sections
void MatrixUnf::GetBinWidthCorrectedCrossSec(TH1D*& uncorrHist, TString VariableName) {

  // Read in binning schemes in XML format
  TDOMParser parser;
  Int_t error = parser.ParseFile("binning/" + VariableName + "_binning_rebinnedB.xml");
  if(error) std::cout << "error=" << error << " from TDOMParser\n";

  TXMLDocument const* XMLdocument     = parser.GetXMLDocument();
  TUnfoldBinningXML* generatorBinning = TUnfoldBinningXML::ImportXML(XMLdocument,"generator_rebinnedB");

  if(generatorBinning)  generatorBinning->PrintStream(cout,1,1); 
  else {std::cout<<"could not read 'generator_rebinnedB' binning\n"; }

  const TUnfoldBinning *binning = generatorBinning->FindNode("ttbargen_rebinnedB");

  // GetBinSize() automatically gets the N-D size of the bin
  // Ignore underflow and overflow and scale bin contents and errors
  for(Int_t i =1; i <=uncorrHist->GetNbinsX(); i++) {
    Double_t cF = binning->GetBinSize(i);
    uncorrHist->SetBinContent(i, uncorrHist->GetBinContent(i)/cF);
    uncorrHist->SetBinError(i, uncorrHist->GetBinError(i)/cF);
  }

} 
// End


void MatrixUnf::GetAbsoluteCovarianceMatrix(TH1D* unfHist, TH2D*& covMat, Double_t BR, Double_t Lumi) {
  // correct the result for BR, Lumi & store it
  for (Int_t i = 1; i <= covMat->GetNbinsX(); i++) {
   for (Int_t j = 1; j <= covMat->GetNbinsY(); j++) {
    Double_t cFij   = pow(BR * Lumi,2);
    //Double_t cFij = pow(BR*Lumi,2)*unfHist->GetXaxis()->GetBinWidth(i)*unfHist->GetXaxis()->GetBinWidth(j);
    covMat->SetBinContent(i, j, covMat->GetBinContent(i,j)/cFij);
   }
  }
} 

TMatrixD MatrixUnf::TH2toTMatrixD(TH2* covMhist){

  const int nbins = covMhist->GetNbinsX();
  TMatrixD covarianceM(nbins, nbins);

  for(int i=0; i < nbins; i++) { 
    for(int j=0; j < nbins; j++) {
    covarianceM(i,j) = covMhist->GetBinContent(i+1,j+1);
    }
  }
  return covarianceM;
}

TH2D* MatrixUnf::TMatrixDtoTH2(TMatrixD &covM, TString name){

  const int nbins = covM.GetNrows();
  TH2D* covMhist = new TH2D(name, name, nbins, 0, nbins, nbins, 0, nbins);

  for(int i=0;i<nbins;i++)   for(int j=0;j<nbins;j++)  {
    covMhist->SetBinContent(i+1,j+1,covM(i,j));
  }

  return covMhist;
}

TVectorD MatrixUnf::TH1toTVectorD(TH1* deltahist){

  const int nbins = deltahist->GetNbinsX();
  TVectorD diag(nbins);

  for(int i=0;i<nbins;i++) {
    diag(i) = deltahist->GetBinContent(i+1);
  }

  return diag;
}

TMatrixD MatrixUnf::TH1toTMatrixD(TH1* deltahist){

  const int nbins = deltahist->GetNbinsX();
  TMatrixD diag(nbins,1);

  for(int i=0;i<nbins;i++) {
    diag(i,0) = deltahist->GetBinContent(i+1);
  }

  return diag;
}

void MatrixUnf::GetNormalisedCovarianceMatrix(TH1D* histogram, TH2D* covMhist, TH2D*& NormcovMhist){

  TMatrixD NormcovarianceM(NormcovMhist->GetNbinsX(), NormcovMhist->GetNbinsY());

  const int nbins  = histogram->GetNbinsX();
  double hintegral = histogram->Integral();
  //std::cout << "MatrixUnf::GetNormalisedCovarianceMatrix nbins: " << nbins << std::endl;
  //std::cout << "MatrixUnf::GetNormalisedCovarianceMatrix integral: " << hintegral << std::endl;

  double n[nbins];
  for(int i=0; i < nbins; i++){
    n[i]  = histogram->GetBinContent(i+1);
    n[i] /= hintegral;
  }

  double sum_n = 0.;
  for(int i=0;i<nbins;i++) sum_n += n[i] ;

  if (fabs(sum_n-1.)>2e-6) cout<<"WARNING: input distribution not normalised to unit area: "<<(sum_n-1.)<<endl;
  sum_n = 1.;

  // Amandeep : This is the Jacobian
  double dfdn[nbins][nbins];

  // And its calculation given that the operation is normalization
  for(int k=0; k < nbins; k++){ 
    for(int j=0; j < nbins; j++){
    if(k==j) dfdn[k][j]  = sum_n-n[k] / pow(sum_n,2) ;
    else     dfdn[k][j]  = -n[k] / pow(sum_n,2) ;
    NormcovarianceM(k,j) = 0;
    }
  }

  for(int k=0;k<nbins;k++)   for(int i=0;i<nbins;i++)   for(int j=0;j<nbins;j++)   for(int l=0;l<nbins;l++){
    // Amandeep : When normalized the uncertainty transforms as : sigma_norm = J * sigma_abs * J^T
    // where J is the Jacobian as computed above and sigma_abs is the original covariance matrix
    NormcovarianceM(k,l)   += covMhist->GetBinContent(i+1,j+1) * dfdn[k][i] * dfdn[l][j] / hintegral / hintegral;
  }

  for(int k=0;k<nbins;k++)   for(int j=0;j<nbins;j++)  {
    NormcovMhist->SetBinContent(k+1, j+1, NormcovarianceM(k,j));
  }
}


void MatrixUnf::GetNormalisedCovarianceMatrixAllVars(TH1D* histogram, TH2D* covMhist, TH2D*& NormcovMhist){

  TMatrixD NormcovarianceM(NormcovMhist->GetNbinsX(), NormcovMhist->GetNbinsY());

  const int nbins = histogram->GetNbinsX();
  //double hintegral=histogram->Integral();
  //std::cout<<"MatrixUnf::GetNormalisedCovarianceMatrix nbins: "<<nbins<<std::endl;
  //std::cout<<"MatrixUnf::GetNormalisedCovarianceMatrix integral: "<<hintegral<<std::endl;

  // Amandeep : Hardcoded here
  // const int nbinsrhoi = 24;
  const int nbinsrhoi    = 6;

  double n[nbins];
  double hintegral[nbins/nbinsrhoi] = {0.};
  for(int i=0;i<nbins;i++){
    n[i] = histogram->GetBinContent(i+1);
    if(i%nbinsrhoi==0) hintegral[i/nbinsrhoi] = histogram->Integral(nbinsrhoi*(i/nbinsrhoi)+1,nbinsrhoi*(i/nbinsrhoi)+nbinsrhoi);
    n[i] /= hintegral[i/nbinsrhoi];
  }

  double sum_n[nbins/nbinsrhoi] = {0.};
  for(int i=0;i<nbins;i++) sum_n[i/nbinsrhoi] += n[i] ;

  for (int i_var = 0; i_var < nbins/nbinsrhoi; ++i_var)
  {
    if (fabs(sum_n[i_var]-1.)>2e-6) cout<<"WARNING in GetNormalisedCovarianceMatrixAllVars: input distribution not normalised to unit area: "<<(sum_n[i_var]-1.)<<" "<<i_var<<endl;
    sum_n[i_var] = 1.;
  }

  double dfdn[nbins][nbins];

  for(int k=0;k<nbins;k++)   for(int j=0;j<nbins;j++)  {
    if(k==j) dfdn[k][j] = sum_n[k/nbinsrhoi]-n[k] / pow(sum_n[k/nbinsrhoi],2) ;
    else if(k/nbinsrhoi==j/nbinsrhoi) dfdn[k][j] = -n[k] / pow(sum_n[k/nbinsrhoi],2) ;
    else dfdn[k][j] = 0 ;
    NormcovarianceM(k,j) = 0;
  }

  for(int k=0;k<nbins;k++)   for(int i=0;i<nbins;i++)   for(int j=0;j<nbins;j++)   for(int l=0;l<nbins;l++)   {
    // NormcovarianceM(k,l) += covarianceM(i,j) * dfdn[k][i] * dfdn[l][j] / hintegral / hintegral;
    NormcovarianceM(k,l)    += covMhist->GetBinContent(i+1,j+1) * dfdn[k][i] * dfdn[l][j] / hintegral[i/nbinsrhoi] / hintegral[j/nbinsrhoi];
  }

  for(int k=0;k<nbins;k++)   for(int j=0;j<nbins;j++)  {
    NormcovMhist->SetBinContent(k+1, j+1, NormcovarianceM(k,j));
    //if(nbins<7) std::cout<<"NormcovMhist: "<<NormcovarianceM(k,j)<<std::endl;
  }
}

/**
 * @brief Extracts the optimal coefficient from the asymmetry somehow
 * 
 * @param Afbvec The vector of asymmetries in the distribution. In order from outward asymmetry to inward asymmetry, 
 *               i.e., first asymmetry is between first and last bin
 * @param m_AFB Covariance matrix for the asymmetries
 * @param genhistogram The generator-level distribution to extract the bin-by-bin asymmetries and compute asymmetry-to-coefficient conversion factors from
 * @param VariableName1 The name of the variable/coefficient that we are extracting
 * @param Coef The optimal coefficients that were extracted
 * @param CoefErrOpt Uncertainties on coefficients that were extracted
 * @param CoefErrTest Uncertainties on the generator-level coefficients that were extracted
 * @param Coefvec The bin-by-bin coefficients
 * @param m_C Covariance matrix for the coefficients converted from each asymmetry.
 * @param a_optimal_reco Relative weighting of the bin-by-bin coefficients for the optimal coefficient
 * @param numCoefficients The number of coefficients we will be extracting
 * @return true if asymmetry-to-coefficient conversion factors do exist
 * @return false if no asymmetry-to-coefficient conversion factors exist
 */
bool MatrixUnf::GetOptimalCoefficient(std::vector<double>& Afbvec,
                                      TMatrixD& m_AFB,
                                      TH1D* genhistogram,
                                      TString VariableName1,
                                      Double_t* Coef,
                                      Double_t* CoefErrOpt,
                                      Double_t* CoefErrTest,
                                      std::vector<double>& Coefvec,
                                      TMatrixD& m_C,
                                      std::vector<double>& a_optimal_reco,
                                      int numCoefficients) {
  Coefvec.clear();
  a_optimal_reco.clear();
  int nbinsrhoi = genhistogram->GetNbinsX();
  int numBinsPerCoefficient = nbinsrhoi / numCoefficients;

  std::vector<double> zerovec((nbinsrhoi + 1) / 2);
  a_optimal_reco = zerovec;
  Coefvec = zerovec;

  std::vector<double> Afbvec_gen, AfbErrvec_gen, CoefErrSqvec_gen;
  MatrixUnf::GetAfbBinByBin(genhistogram, Afbvec_gen, AfbErrvec_gen, numCoefficients);

  for (int i = 0; i < numCoefficients; i++) {
    int AfbOffset = i * numBinsPerCoefficient / 2;
    int offset = i * numBinsPerCoefficient;

    Coef[i] = 0.;
    CoefErrOpt[i] = 0.;
    CoefErrTest[i] = 0.;

    std::vector<double> factor{1., 1., 1.};
    bool factorsExist = MatrixUnf::AtoCfactor(VariableName1, genhistogram, factor, numCoefficients, i);

    if (!factorsExist) {
      if (numBinsPerCoefficient != 6) {
        cout << "Warning in MatrixUnf::GetOptimalCoefficient, only " << numBinsPerCoefficient
            << " bins and at least 6 are assumed. May cause problems when using a_optimal_reco later." << endl;
      }

      for (int j = 0; j < numCoefficients; j++) {
        Coef[j] = 0.0;
        CoefErrOpt[j] = 0.0;
        CoefErrTest[j] = 0.0;
      }
      return kFALSE;
    } else {
      std::vector<double> Ciivec{Afbvec[AfbOffset + 0] / factor[0], Afbvec[AfbOffset + 1] / factor[1], Afbvec[AfbOffset + 2] / factor[2]};
      double Ciimean = Afbvec[AfbOffset + 0] / factor[0] + Afbvec[AfbOffset + 1] / factor[1] + Afbvec[AfbOffset + 2] / factor[2];
      double Ciimeanerr = 0.;

      for (int iii = 0; iii < numBinsPerCoefficient / 2; iii++) {
        for (int j = 0; j < numBinsPerCoefficient / 2; j++) {
          double m_Cij = m_AFB(AfbOffset + iii, AfbOffset + j) / factor[iii] / factor[j];
          m_C(AfbOffset + iii, AfbOffset + j) = m_Cij;
          Ciimeanerr += m_Cij;
        }
      }
      Ciimeanerr = sqrt(Ciimeanerr);

      Ciimeanerr = Ciimeanerr / 3.;
      Ciimean = Ciimean / 3.;

      CoefErrSqvec_gen.clear();
      for (int iii = 0; iii < numBinsPerCoefficient / 2; iii++) {
        CoefErrSqvec_gen.push_back(pow(AfbErrvec_gen[AfbOffset + iii] / factor[iii], 2));
      }

      std::vector<double> a_optimal_gen{0, 0, 0};
      a_optimal_gen[0] = CoefErrSqvec_gen[1] * CoefErrSqvec_gen[2];
      a_optimal_gen[1] = CoefErrSqvec_gen[0] * CoefErrSqvec_gen[2];
      a_optimal_gen[2] = CoefErrSqvec_gen[0] * CoefErrSqvec_gen[1];

      double sum_a_optimal_gen = a_optimal_gen[0] + a_optimal_gen[1] + a_optimal_gen[2];

      a_optimal_gen[0] /= sum_a_optimal_gen;
      a_optimal_gen[1] /= sum_a_optimal_gen;
      a_optimal_gen[2] /= sum_a_optimal_gen;

      double Ciigenopt = Ciivec[0] * a_optimal_gen[0] + Ciivec[1] * a_optimal_gen[1] + Ciivec[2] * a_optimal_gen[2];
      double Ciigenopterr = 0.;

      for (int iii = 0; iii < numBinsPerCoefficient / 2; iii++) {
        for (int j = 0; j < numBinsPerCoefficient / 2; j++) {
          Ciigenopterr += m_C(AfbOffset + iii, AfbOffset + j) * a_optimal_gen[iii] * a_optimal_gen[j];
        }
      }
      Ciigenopterr = sqrt(Ciigenopterr);

      std::vector<double> a_optimal{0, 0, 0};
      a_optimal[0] = m_C(AfbOffset + 1, AfbOffset + 1) * m_C(AfbOffset + 2, AfbOffset + 2) - m_C(AfbOffset + 1, AfbOffset + 2) * m_C(AfbOffset + 1, AfbOffset + 2) - m_C(AfbOffset + 0, AfbOffset + 1) * (m_C(AfbOffset + 2, AfbOffset + 2) - m_C(AfbOffset + 1, AfbOffset + 2)) -
                    m_C(AfbOffset + 0, AfbOffset + 2) * (m_C(AfbOffset + 1, AfbOffset + 1) - m_C(AfbOffset + 1, AfbOffset + 2));
      a_optimal[1] = m_C(AfbOffset + 0, AfbOffset + 0) * m_C(AfbOffset + 2, AfbOffset + 2) - m_C(AfbOffset + 0, AfbOffset + 2) * m_C(AfbOffset + 0, AfbOffset + 2) - m_C(AfbOffset + 0, AfbOffset + 1) * (m_C(AfbOffset + 2, AfbOffset + 2) - m_C(AfbOffset + 0, AfbOffset + 2)) -
                    m_C(AfbOffset + 1, AfbOffset + 2) * (m_C(AfbOffset + 0, AfbOffset + 0) - m_C(AfbOffset + 0, AfbOffset + 2));
      a_optimal[2] = m_C(AfbOffset + 0, AfbOffset + 0) * m_C(AfbOffset + 1, AfbOffset + 1) - m_C(AfbOffset + 0, AfbOffset + 1) * m_C(AfbOffset + 0, AfbOffset + 1) - m_C(AfbOffset + 0, AfbOffset + 2) * (m_C(AfbOffset + 1, AfbOffset + 1) - m_C(AfbOffset + 0, AfbOffset + 1)) -
                    m_C(AfbOffset + 1, AfbOffset + 2) * (m_C(AfbOffset + 0, AfbOffset + 0) - m_C(AfbOffset + 0, AfbOffset + 1));

      double sum_a_optimal = a_optimal[0] + a_optimal[1] + a_optimal[2];

      a_optimal[0] /= sum_a_optimal;
      a_optimal[1] /= sum_a_optimal;
      a_optimal[2] /= sum_a_optimal;

      // don't allow combination with negative weights
      if (a_optimal[0] < 0 || a_optimal[1] < 0 || a_optimal[2] < 0) {
        // test if any two variables can be combined without negative weights. If so, do it; if not, use the single measurement with lowest variance.
        // an arbitrary starting value that will always be larger than the single or combined coefficient error^2
        double min_var = m_C(AfbOffset + 0, AfbOffset + 0) + m_C(AfbOffset + 1, AfbOffset + 1) + m_C(AfbOffset + 2, AfbOffset + 2);  
        double min_var_single = min_var;
        int i_pair_min = -1;
        int i_single_min = -1;
        std::vector<double> v_pair{0, 0, 0};

        for (int i_pair = 0; i_pair < 3; ++i_pair) {
          int j_pair = (i_pair + 1) % 3;
          // if the pair can be combined without one measurement having a negative weight
          if (m_C(AfbOffset + i_pair, AfbOffset + j_pair) <= m_C(AfbOffset + i_pair, AfbOffset + i_pair) 
           && m_C(AfbOffset + i_pair, AfbOffset + j_pair) <= m_C(AfbOffset + j_pair, AfbOffset + j_pair)) {  
            v_pair[i_pair] = (m_C(AfbOffset + i_pair, AfbOffset + i_pair) * m_C(AfbOffset + j_pair, AfbOffset + j_pair) - m_C(AfbOffset + i_pair, AfbOffset + j_pair) * m_C(AfbOffset + j_pair, AfbOffset + i_pair)) /
                            (m_C(AfbOffset + i_pair, AfbOffset + i_pair) + m_C(AfbOffset + j_pair, AfbOffset + j_pair) - m_C(AfbOffset + i_pair, AfbOffset + j_pair) - m_C(AfbOffset + j_pair, AfbOffset + i_pair));
            if (v_pair[i_pair] < min_var) {
              min_var = v_pair[i_pair];
              i_pair_min = i_pair;
            }
          }
          if (m_C(AfbOffset + i_pair, AfbOffset + i_pair) < min_var_single) {
            min_var_single = m_C(AfbOffset + i_pair, AfbOffset + i_pair);
            i_single_min = i_pair;
          }
        }
        // cout<<"Info in MatrixUnf::GetOptimalCoefficient, chose combination "<<i_pair_min<<" with lowest uncertainty^2 from "<<v_pair[0]<<" "<<v_pair[1]<<" "<<v_pair[2]<<endl;
        // cout<<"Info in MatrixUnf::GetOptimalCoefficient, chose measurement "<<i_single_min<<" with lowest uncertainty^2 from "<<m_C(0,0)<<" "<<m_C(1,1)<<" "<<m_C(2,2)<<endl;

        if (i_pair_min != -1) {
          int j_pair_min = (i_pair_min + 1) % 3;
          int k_pair_min = (i_pair_min + 2) % 3;
          a_optimal[i_pair_min] =
              (m_C(AfbOffset + j_pair_min, AfbOffset + j_pair_min) - m_C(AfbOffset + j_pair_min, AfbOffset + i_pair_min)) /
              (m_C(AfbOffset + i_pair_min, AfbOffset + i_pair_min) - 2. * m_C(AfbOffset + j_pair_min, AfbOffset + i_pair_min) + m_C(AfbOffset + j_pair_min, AfbOffset + j_pair_min));
          a_optimal[j_pair_min] = 1. - a_optimal[i_pair_min];
          a_optimal[k_pair_min] = 0;
        } else if (i_single_min != -1) {
          int j_single_min = (i_single_min + 1) % 3;
          int k_single_min = (i_single_min + 2) % 3;
          a_optimal[i_single_min] = 1;
          a_optimal[j_single_min] = 0;
          a_optimal[k_single_min] = 0;
        } else
          cout << "Error in MatrixUnf::GetOptimalCoefficient, no solution found: " << i_pair_min << " " << i_single_min
              << endl;
      }

      // cout<<VariableName1<<" a_optimal: "<<a_optimal[0]<<" "<<a_optimal[1]<<" "<<a_optimal[2]<<endl;

      // check again for negative weights (should never happen)
      if (a_optimal[0] < 0 || a_optimal[1] < 0 || a_optimal[2] < 0) {
        cout << "Warning in MatrixUnf::GetOptimalCoefficient, negative weights in combination: a_optimal[] = "
            << a_optimal[0] << " " << a_optimal[1] << " " << a_optimal[2] << endl;
        m_C.Print("f=%1.10g ");
        a_optimal[0] = a_optimal_gen[0];
        a_optimal[1] = a_optimal_gen[1];
        a_optimal[2] = a_optimal_gen[2];
      }

      double Cii_optimal = Ciivec[0] * a_optimal[0] + Ciivec[1] * a_optimal[1] + Ciivec[2] * a_optimal[2];
      double Ciierr_optimal = 0.;

      for (int iii = 0; iii < numBinsPerCoefficient / 2; iii++) {
        for (int j = 0; j < numBinsPerCoefficient / 2; j++) {
          Ciierr_optimal += m_C(AfbOffset + iii, AfbOffset + j) * a_optimal[iii] * a_optimal[j];
        }
      }
      Ciierr_optimal = sqrt(Ciierr_optimal);

      for (int iii = 0; iii < numBinsPerCoefficient / 2; iii++) {
        Coefvec[AfbOffset + iii] = Ciivec[iii];
        a_optimal_reco[AfbOffset + iii] = a_optimal[iii];
      }
      // CoefErrTest = Ciimeanerr;
      CoefErrTest[i] = Ciigenopterr;
      Coef[i] = Cii_optimal;
      CoefErrOpt[i] = Ciierr_optimal;
    }
  }
  return kTRUE;
}


TString MatrixUnf::CoefName(TString VariableName1, bool islatex, double lowerBoundary, double upperBoundary) {
  TString quantity;
  TString otherQuantity;

  int PRECISION_VAL = -1; // how many places after decimal to report

  if(VariableName1.Contains("delta_phi")) quantity = islatex? "A_{|\\Delta\\phi_{\\ell \\bar{\\ell}}|}" : "A_{|#Delta#phi_{l#bar{l}}|}";
  else if(VariableName1.Contains("delta_eta")) quantity = islatex? "A_{|\\Delta\\eta_{\\ell \\bar{\\ell}}|}" : "A_{|#Delta#eta_{l#bar{l}}|}";

  else if(VariableName1.Contains("1n")) quantity = "B_{1}^{n}";
  else if(VariableName1.Contains("2n")) quantity = "B_{2}^{n}";
  else if(VariableName1.Contains("1r")) quantity = "B_{1}^{r}";
  else if(VariableName1.Contains("2r")) quantity = "B_{2}^{r}";
  else if(VariableName1.Contains("1k")) quantity = "B_{1}^{k}";
  else if(VariableName1.Contains("2k")) quantity = "B_{2}^{k}";
  else if(VariableName1.Contains("1j")) quantity = "B_{1}^{k*}";
  else if(VariableName1.Contains("2j")) quantity = "B_{2}^{k*}";
  else if(VariableName1.Contains("1q")) quantity = "B_{1}^{r*}";
  else if(VariableName1.Contains("2q")) quantity = "B_{2}^{r*}";

  else if(VariableName1.Contains("Pnn")) quantity = "B_{1}^{n}+B_{2}^{n}";
  else if(VariableName1.Contains("Mnn")) quantity = "B_{1}^{n}-B_{2}^{n}";
  else if(VariableName1.Contains("Prr")) quantity = "B_{1}^{r}+B_{2}^{r}";
  else if(VariableName1.Contains("Mrr")) quantity = "B_{1}^{r}-B_{2}^{r}";
  else if(VariableName1.Contains("Pkk")) quantity = "B_{1}^{k}+B_{2}^{k}";
  else if(VariableName1.Contains("Mkk")) quantity = "B_{1}^{k}-B_{2}^{k}";
  else if(VariableName1.Contains("Pjj")) quantity = "B_{1}^{k*}+B_{2}^{k*}";
  else if(VariableName1.Contains("Mjj")) quantity = "B_{1}^{k*}-B_{2}^{k*}";
  else if(VariableName1.Contains("Pqq")) quantity = "B_{1}^{r*}+B_{2}^{r*}";
  else if(VariableName1.Contains("Mqq")) quantity = "B_{1}^{r*}-B_{2}^{r*}";

  else if(VariableName1.Contains("c_han")) quantity = islatex? "+ C_{kk} - C_{rr} - C_{nn}" : "+C_{kk}#kern[-0.6]{ }-#kern[-0.3]{ }C_{rr}#kern[-0.6]{ }-#kern[-0.3]{ }C_{nn}";
  else if(VariableName1.Contains("c_sca")) quantity = islatex? "- C_{kk} + C_{rr} - C_{nn}" : "-C_{kk}#kern[-0.6]{ }+#kern[-0.3]{ }C_{rr}#kern[-0.6]{ }-#kern[-0.3]{ }C_{nn}";
  else if(VariableName1.Contains("c_tra")) quantity = islatex? "- C_{kk} - C_{rr} + C_{nn}" : "-C_{kk}#kern[-0.6]{ }-#kern[-0.3]{ }C_{rr}#kern[-0.6]{ }+#kern[-0.3]{ }C_{nn}";
  else if(VariableName1.Contains("c_kjL")) quantity = islatex? "- C_{kk*} - C_{rr} - C_{nn}" : "-C_{kk*}#kern[-0.6]{ }-#kern[-0.3]{ }C_{rr}#kern[-0.6]{ }-#kern[-0.3]{ }C_{nn}";
  else if(VariableName1.Contains("c_rqL")) quantity = islatex? "- C_{kk} - C_{rr*} - C_{nn}" : "-C_{kk}#kern[-0.6]{ }-#kern[-0.3]{ }C_{rr*}#kern[-0.6]{ }-#kern[-0.3]{ }C_{nn}";

  else if(VariableName1.Contains("c_rkP")) quantity = islatex? "- C_{rk} - C_{kr} - C_{nn}" : "-C_{rk}#kern[-0.6]{ }-#kern[-0.3]{ }C_{kr}#kern[-0.6]{ }+#kern[-0.3]{ }C_{nn}";
  else if(VariableName1.Contains("c_rkM")) quantity = islatex? "- C_{rk} + C_{kr} - C_{nn}" : "-C_{rk}#kern[-0.6]{ }+#kern[-0.3]{ }C_{kr}#kern[-0.6]{ }-#kern[-0.3]{ }C_{nn}";
  else if(VariableName1.Contains("c_nrP")) quantity = islatex? "- C_{nr} - C_{rn} - C_{kk}" : "-C_{nr}#kern[-0.6]{ }-#kern[-0.3]{ }C_{rn}#kern[-0.6]{ }-#kern[-0.3]{ }C_{kk}";
  else if(VariableName1.Contains("c_nrM")) quantity = islatex? "- C_{nr} + C_{rn} - C_{kk}" : "-C_{nr}#kern[-0.6]{ }+#kern[-0.3]{ }C_{rn}#kern[-0.6]{ }+#kern[-0.3]{ }C_{kk}";
  else if(VariableName1.Contains("c_nkP")) quantity = islatex? "+ C_{nk} - C_{kn} - C_{rr}" : "-C_{nk}#kern[-0.6]{ }-#kern[-0.3]{ }C_{kn}#kern[-0.6]{ }-#kern[-0.3]{ }C_{rr}";
  else if(VariableName1.Contains("c_nkM")) quantity = islatex? "- C_{nk} + C_{kn} - C_{rr}" : "-C_{nk}#kern[-0.6]{ }+#kern[-0.3]{ }C_{kn}#kern[-0.6]{ }-#kern[-0.3]{ }C_{rr}";

  else if(VariableName1.Contains("c_kk")) quantity = "C_{kk}";
  else if(VariableName1.Contains("c_kn")) quantity = "C_{kn}";
  else if(VariableName1.Contains("c_kr")) quantity = "C_{kr}";
  else if(VariableName1.Contains("c_nk")) quantity = "C_{nk}";
  else if(VariableName1.Contains("c_nn")) quantity = "C_{nn}";
  else if(VariableName1.Contains("c_nr")) quantity = "C_{nr}";
  else if(VariableName1.Contains("c_rk")) quantity = "C_{rk}";
  else if(VariableName1.Contains("c_rn")) quantity = "C_{rn}";
  else if(VariableName1.Contains("c_rr")) quantity = "C_{rr}";
  else if(VariableName1.Contains("c_kj")) quantity = "C_{kk*}";
  else if(VariableName1.Contains("c_rq")) quantity = "C_{rr*}";

  else if(VariableName1.Contains("c_Prk")) quantity = islatex? "C_{rk}+C_{kr}" : "C_{rk}+#kern[-0.8]{ }C_{kr}";
  else if(VariableName1.Contains("c_Mrk")) quantity = islatex? "C_{rk}-C_{kr}" : "C_{rk}#kern[-0.6]{ }-#kern[-0.3]{ }C_{kr}";
  else if(VariableName1.Contains("c_Pnr")) quantity = islatex? "C_{nr}+C_{rn}" : "C_{nr}+#kern[-0.8]{ }C_{rn}";
  else if(VariableName1.Contains("c_Mnr")) quantity = islatex? "C_{nr}-C_{rn}" : "C_{nr}#kern[-0.6]{ }-#kern[-0.3]{ }C_{rn}";
  else if(VariableName1.Contains("c_Pnk")) quantity = islatex? "C_{nk}+C_{kn}" : "C_{nk}+#kern[-0.8]{ }C_{kn}";
  else if(VariableName1.Contains("c_Mnk")) quantity = islatex? "C_{nk}-C_{kn}" : "C_{nk}#kern[-0.6]{ }-#kern[-0.3]{ }C_{kn}";
  else if(VariableName1.Contains("c_Prj")) quantity = islatex? "C_{rk*}+C_{k*r}" : "C_{rk*}+#kern[-0.8]{ }C_{k*r}";
  else if(VariableName1.Contains("c_Mrj")) quantity = islatex? "C_{rk*}-C_{k*r}" : "C_{rk*}#kern[-0.6]{ }-#kern[-0.3]{ }C_{k*r}";

  else if(VariableName1.Contains("cHel")) quantity = "D";
  else if(VariableName1.Contains("cLab")) quantity = islatex? "A_{\\cos\\varphi}^{\\mathrm{lab}}" : "A_{cos#varphi}^{lab}";
  else if(VariableName1.Contains("kNorm")) quantity = islatex? "\\sin\\theta" : "sin#theta";
  else if(VariableName1.Contains("rNorm")) quantity = islatex? "\\sin\\theta" : "sin#theta";
  else quantity="";

  if (VariableName1.Contains("mttbar")) {
    otherQuantity = islatex ? "m_{t\\bar{t}}" : "m_{t#bar{t}}";
  } else {
    otherQuantity = "";
  }

  if (!otherQuantity.EqualTo("")) {
    if (quantity.Contains("_")) {
      quantity.Append("_{, ");
    } else {
      quantity.Append("_{");
    }
    if (lowerBoundary != -999) {
      std::string strLowerBoundary = std::to_string(lowerBoundary).substr(0, std::to_string(lowerBoundary).find(".") + PRECISION_VAL + 1);
      quantity.Append(strLowerBoundary);
      quantity.Append(" < ");
    }
    quantity.Append(otherQuantity);
    if (upperBoundary != -999) {
      quantity.Append(" < ");
      std::string strUpperBoundary = std::to_string(upperBoundary).substr(0, std::to_string(upperBoundary).find(".") + PRECISION_VAL + 1);
      quantity.Append(strUpperBoundary);
    }
    quantity.Append("}");
  }


  return quantity;

}

TString MatrixUnf::ObsName(TString VariableName1, bool islatex){

  TString quantity;

  if(islatex) {
    if(VariableName1.Contains("delta_phi")) quantity = "|\\Delta\\phi_{\\ell\\bar{\\ell}}|";
    else if(VariableName1.Contains("delta_eta")) quantity = "|\\Delta\\eta_{\\ell\\bar{\\ell}}|";

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
    else if(VariableName1.Contains("c_rq")) quantity = "\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{r*}";

    else if(VariableName1.Contains("c_Prk")) quantity = "\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{k}+\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{r}";
    else if(VariableName1.Contains("c_Mrk")) quantity = "\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{k}-\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{r}";
    else if(VariableName1.Contains("c_Pnr")) quantity = "\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{r}+\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("c_Mnr")) quantity = "\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{r}-\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("c_Pnk")) quantity = "\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{k}+\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("c_Mnk")) quantity = "\\cos\\theta_{1}^{n}\\cos\\theta_{2}^{k}-\\cos\\theta_{1}^{k}\\cos\\theta_{2}^{n}";
    else if(VariableName1.Contains("c_Prj")) quantity = "\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{k*}+\\cos\\theta_{1}^{k*}\\cos\\theta_{2}^{r}";
    else if(VariableName1.Contains("c_Mrj")) quantity = "\\cos\\theta_{1}^{r}\\cos\\theta_{2}^{k*}-\\cos\\theta_{1}^{k*}\\cos\\theta_{2}^{r}";

    else if(VariableName1.Contains("cHel")) quantity = "\\cos\\varphi";
    else if(VariableName1.Contains("cLab")) quantity = "\\cos\\varphi_{\\mathrm{lab}}";
    else if(VariableName1.Contains("kNorm")) quantity = "\\sin\\theta";
    else if(VariableName1.Contains("rNorm")) quantity = "\\sin\\theta";
    else quantity="";
  }
  else {
    if(VariableName1.Contains("delta_phi")) quantity = "|#Delta#phi_{l#bar{l}}|";
    else if(VariableName1.Contains("delta_eta")) quantity = "|#Delta#eta_{l#bar{l}}|";

    else if(VariableName1.Contains("1n")) quantity = "cos#kern[-0.8]{ }#theta_{1}^{n}";
    else if(VariableName1.Contains("2n")) quantity = "cos#kern[-0.8]{ }#theta_{2}^{n}";
    else if(VariableName1.Contains("1r")) quantity = "cos#kern[-0.8]{ }#theta_{1}^{r}";
    else if(VariableName1.Contains("2r")) quantity = "cos#kern[-0.8]{ }#theta_{2}^{r}";
    else if(VariableName1.Contains("1k")) quantity = "cos#kern[-0.8]{ }#theta_{1}^{k}";
    else if(VariableName1.Contains("2k")) quantity = "cos#kern[-0.8]{ }#theta_{2}^{k}";
    else if(VariableName1.Contains("1j")) quantity = "cos#kern[-0.8]{ }#theta_{1}^{k*}";
    else if(VariableName1.Contains("2j")) quantity = "cos#kern[-0.8]{ }#theta_{2}^{k*}";
    else if(VariableName1.Contains("1q")) quantity = "cos#kern[-0.8]{ }#theta_{1}^{r*}";
    else if(VariableName1.Contains("2q")) quantity = "cos#kern[-0.8]{ }#theta_{2}^{r*}";

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

    else if(VariableName1.Contains("c_han")) quantity = "+cos#kern[-0.9]{ }#theta_{1}^{k}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{k}#kern[-0.6]{ }-#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{r}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{r}#kern[-0.6]{ }-#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{n}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{n}";
    else if(VariableName1.Contains("c_sca")) quantity = "-cos#kern[-0.9]{ }#theta_{1}^{k}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{k}#kern[-0.6]{ }+#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{r}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{r}#kern[-0.6]{ }-#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{n}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{n}";
    else if(VariableName1.Contains("c_tra")) quantity = "-cos#kern[-0.9]{ }#theta_{1}^{k}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{k}#kern[-0.6]{ }-#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{r}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{r}#kern[-0.6]{ }+#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{n}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{n}";
    else if(VariableName1.Contains("c_kjL")) quantity = "-cos#kern[-0.9]{ }#theta_{1}^{k}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{k*}#kern[-0.6]{ }-#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{r}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{r}#kern[-0.6]{ }-#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{n}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{n}";
    else if(VariableName1.Contains("c_rqL")) quantity = "-cos#kern[-0.9]{ }#theta_{1}^{k}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{k}#kern[-0.6]{ }-#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{r}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{r*}#kern[-0.6]{ }-#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{n}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{n}";

    else if(VariableName1.Contains("c_rkP")) quantity = "-cos#kern[-0.9]{ }#theta_{1}^{r}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{k}#kern[-0.6]{ }-#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{k}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{r}#kern[-0.6]{ }-#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{n}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{n}";
    else if(VariableName1.Contains("c_rkM")) quantity = "-cos#kern[-0.9]{ }#theta_{1}^{r}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{k}#kern[-0.6]{ }+#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{k}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{r}#kern[-0.6]{ }-#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{n}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{n}";
    else if(VariableName1.Contains("c_nrP")) quantity = "-cos#kern[-0.9]{ }#theta_{1}^{n}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{r}#kern[-0.6]{ }-#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{r}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{n}#kern[-0.6]{ }-#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{k}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{k}";
    else if(VariableName1.Contains("c_nrM")) quantity = "-cos#kern[-0.9]{ }#theta_{1}^{n}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{r}#kern[-0.6]{ }+#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{r}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{n}#kern[-0.6]{ }-#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{k}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{k}";
    else if(VariableName1.Contains("c_nkP")) quantity = "-cos#kern[-0.9]{ }#theta_{1}^{n}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{k}#kern[-0.6]{ }-#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{k}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{n}#kern[-0.6]{ }-#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{r}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{r}";
    else if(VariableName1.Contains("c_nkM")) quantity = "-cos#kern[-0.9]{ }#theta_{1}^{n}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{k}#kern[-0.6]{ }+#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{k}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{n}#kern[-0.6]{ }-#kern[-0.4]{ }cos#kern[-0.9]{ }#theta_{1}^{r}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{r}";

    else if(VariableName1.Contains("c_kk")) quantity = "cos#kern[-0.9]{ }#theta_{1}^{k}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{k}";
    else if(VariableName1.Contains("c_kn")) quantity = "cos#kern[-0.9]{ }#theta_{1}^{k}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{n}";
    else if(VariableName1.Contains("c_kr")) quantity = "cos#kern[-0.9]{ }#theta_{1}^{k}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{r}";
    else if(VariableName1.Contains("c_nk")) quantity = "cos#kern[-0.9]{ }#theta_{1}^{n}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{k}";
    else if(VariableName1.Contains("c_nn")) quantity = "cos#kern[-0.9]{ }#theta_{1}^{n}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{n}";
    else if(VariableName1.Contains("c_nr")) quantity = "cos#kern[-0.9]{ }#theta_{1}^{n}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{r}";
    else if(VariableName1.Contains("c_rk")) quantity = "cos#kern[-0.9]{ }#theta_{1}^{r}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{k}";
    else if(VariableName1.Contains("c_rn")) quantity = "cos#kern[-0.9]{ }#theta_{1}^{r}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{n}";
    else if(VariableName1.Contains("c_rr")) quantity = "cos#kern[-0.9]{ }#theta_{1}^{r}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{r}";
    else if(VariableName1.Contains("c_kj")) quantity = "cos#kern[-0.9]{ }#theta_{1}^{k}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{k*}";
    else if(VariableName1.Contains("c_rq")) quantity = "cos#kern[-0.9]{ }#theta_{1}^{r}#kern[-0.8]{ }cos#kern[-0.9]{ }#theta_{2}^{r*}";

    else if(VariableName1.Contains("c_Prk")) quantity = "cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}";
    else if(VariableName1.Contains("c_Mrk")) quantity = "cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}";
    else if(VariableName1.Contains("c_Pnr")) quantity = "cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}";
    else if(VariableName1.Contains("c_Mnr")) quantity = "cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}";
    else if(VariableName1.Contains("c_Pnk")) quantity = "cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}";
    else if(VariableName1.Contains("c_Mnk")) quantity = "cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}";
    else if(VariableName1.Contains("c_Prj")) quantity = "cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r}";
    else if(VariableName1.Contains("c_Mrj")) quantity = "cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r}";

    else if(VariableName1.Contains("cHel")) quantity = "cos#kern[-0.8]{ }#varphi";
    else if(VariableName1.Contains("cLab")) quantity = "cos#kern[-0.8]{ }#varphi_{lab}";
    else if(VariableName1.Contains("kNorm")) quantity = "sin#theta";
    else if(VariableName1.Contains("rNorm")) quantity = "sin#theta";
    else quantity="";

  }

  return quantity;

}


bool MatrixUnf::isCtype(TString VariableName1){
  if(VariableName1.BeginsWith("c_k") && !(VariableName1.Contains("L") || VariableName1.Contains("M") || VariableName1.Contains("P")) ||  
     VariableName1.BeginsWith("c_r") && !(VariableName1.Contains("L") || VariableName1.Contains("M") || VariableName1.Contains("P")) || 
     VariableName1.BeginsWith("c_n") && !(VariableName1.Contains("L") || VariableName1.Contains("M") || VariableName1.Contains("P")) ) {
    return kTRUE;
  }
  else return kFALSE;
}

bool MatrixUnf::isCPMtype(TString VariableName1){
  if(VariableName1.BeginsWith("c_M") || VariableName1.BeginsWith("c_P") || VariableName1.BeginsWith("ll_rNorm") || VariableName1.BeginsWith("ll_kNorm") ){
    return kTRUE;
  }
  else return kFALSE;
}


bool MatrixUnf::isBPMtype(TString VariableName1){
  if(VariableName1.BeginsWith("b_M") || VariableName1.BeginsWith("b_P")){
    return kTRUE;
  }
  else return kFALSE;
}

bool MatrixUnf::isBtype(TString VariableName1){
  if(VariableName1.BeginsWith("b1") || VariableName1.BeginsWith("b2")){
    return kTRUE;
  }
  else if(VariableName1.Contains("costheta") ){
    return kTRUE;
  }
  else return kFALSE;
}

bool MatrixUnf::isDtype(TString VariableName1){
  if(VariableName1.BeginsWith("ll_cH") ){
    return kTRUE;
  }
  else return kFALSE;
}

bool MatrixUnf::isOtherSymType(TString VariableName1){
  // Amandeep : Allowing both eta and phi here
  if(VariableName1.BeginsWith("llbar_delta_") || VariableName1.BeginsWith("ll_cLab")){
    return kTRUE;
  }
  else return kFALSE;
}

// Amandeep : Adding these for new linear combination variables
bool MatrixUnf::isLinearCombtype(TString VariableName1){
  if( (VariableName1.BeginsWith("c_") && VariableName1.Contains("L")) || 
      (VariableName1.BeginsWith("c_") && VariableName1.Contains("M")) || 
      (VariableName1.BeginsWith("c_") && VariableName1.Contains("P")) || 
      (VariableName1.BeginsWith("c_han")) || (VariableName1.BeginsWith("c_sca")) || (VariableName1.BeginsWith("c_tra"))) {
    return kTRUE;
  }
  else return kFALSE;
}
// End

double MatrixUnf::AtoCfactor(TString VariableName1){
  // Factors to convert inclusive asymmetry to the coefficient
  
  double Afactor    =  1;
  double AfactorC   = -4;
  double AfactorB   =  2;
  double AfactorD   = -2;
  double AfactorBPM =  3;
  double AfactorCPM = -16/TMath::Pi();

  if     ( MatrixUnf::isCtype(VariableName1) )   Afactor = AfactorC;
  else if( MatrixUnf::isCPMtype(VariableName1) ) Afactor = AfactorCPM;
  else if( MatrixUnf::isBPMtype(VariableName1) ) Afactor = AfactorBPM;
  else if( MatrixUnf::isBtype(VariableName1) )   Afactor = AfactorB;
  else if( MatrixUnf::isDtype(VariableName1) )   Afactor = AfactorD;

  // Amandeep : This is an assumption
  // If xsection looks like (1 + [0]x) then it should be AfactorB; otherwise AfactorD
  else if( MatrixUnf::isLinearCombtype(VariableName1)) Afactor = AfactorD;
  // End

  else Afactor=1;
  
return Afactor;
}

/**
 * @brief Calculates a conversion factor to convert between asymmetries and coefficients. 
 * 
 * @param VariableName1 The name of the coefficient that is being converted to
 * @param genhistogram The generator-level distribution being used to extract the coefficient
 * @param factor_in The returned per-asymmetry conversion factors
 * @param numCoefficients The number of coefficients in the distribution
 * @param ithCoefficient Which part of the distribution we are extracting the ith coefficient from
 * @return true 
 * @return false 
 */
bool MatrixUnf::AtoCfactor(TString VariableName1, TH1D* genhistogram, std::vector<double>& factor_in, int numCoefficients, int ithCoefficient) {
  int nbinsrhoi             = genhistogram->GetNbinsX();
  int numBinsPerCoefficient = nbinsrhoi / numCoefficients;
  int offset                = ithCoefficient * numBinsPerCoefficient;

  std::vector<double> factorOther{1., 1., 1.};
  bool factorsExist = kTRUE;

  if (numBinsPerCoefficient != 6) {
    factor_in = factorOther;
    factorsExist = kFALSE;
    return factorsExist;
  }

  std::vector<double> oneOverSMshape_sym_default;
  oneOverSMshape_sym_default.clear();

  for (int i_bin = 0; i_bin < numBinsPerCoefficient; ++i_bin) {
    oneOverSMshape_sym_default.push_back(
        2. * (genhistogram->Integral(offset + 1, offset + numBinsPerCoefficient) / double(numBinsPerCoefficient)) /
             (genhistogram->GetBinContent(offset + i_bin + 1) + genhistogram->GetBinContent(offset + numBinsPerCoefficient - i_bin)));
  }

  // Amandeep : How are these computed ??

  // Factors to convert binned asymmetries to the coefficient
  std::vector<double> factorC{-0.77408770278181780  , -0.47339539118637620, -0.12695788682995793};
  std::vector<double> factorCPM{-0.39965929221482027, -0.24483790876271952, -0.07997510138258578};

  // std::vector<double> factorBPMzeroC {7./9.,13./27.,7./45.};
  // The factors to go from asymmetry to B coefficient depend also on the C coefficient of the distribution 
  // (which affects the shape of the symmetric part of the b_* distribution)
  std::vector<double> factorBPM{
      oneOverSMshape_sym_default[0] * 7. / 27.,
      oneOverSMshape_sym_default[1] * 13./ 27.,
      oneOverSMshape_sym_default[2] * 7. / 27.};  
  std::vector<double> factorB{ 5. / 6.,  0.5,  1. / 6.};
  std::vector<double> factorD{-5. / 6., -0.5, -1. / 6.};

  if (MatrixUnf::isCtype(VariableName1))
    factor_in = factorC;
  else if (MatrixUnf::isCPMtype(VariableName1))
    factor_in = factorCPM;
  else if (MatrixUnf::isBPMtype(VariableName1))
    factor_in = factorBPM;
  else if (MatrixUnf::isBtype(VariableName1))
    factor_in = factorB;
  else if (MatrixUnf::isDtype(VariableName1))
    factor_in = factorD;
  
  // Amandeep : Same assumption as we had in the MatrixUnf::AtoCfactor() method
  // If xsection has a plus then factorB, else factorD
  else if (MatrixUnf::isLinearCombtype(VariableName1))
    factor_in = factorD; 
  // End

  else {
    factor_in = factorOther;
    factorsExist = kFALSE;
  }

  return factorsExist;
}


void MatrixUnf::Chi2tests(TUnfoldDensity &unfold, TUnfoldDensity &unfold_rebinnedB, TH1D* TUnfResult, TH1D* TUnfResult_rebinnedA, TH1D* TUnfResult_rebinnedB){
  TH1D *FoldedBack_tempclone = (TH1D*) unfold.GetFoldedOutput(fVariableName1+"InputFoldedBack_tempclone",0,"ttbarreco","*[UO]", 0);
  TH1D *OrigInput = (TH1D*) unfold.GetInput(fVariableName1+"OrigInput",0,"ttbarreco","*[UO]", 0);

  TH1D* BiasDist    = (TH1D*) unfold.GetBias(fVariableName1+"BiasDist");
  TH1D* LxMinusBias = (TH1D*) unfold.GetLxMinusBias(fVariableName1+"LxMinusBias");
  TH2D* LMatrix     = (TH2D*) unfold.GetL(fVariableName1+"LMatrix");
  TH2D* ProbabilityMatrix = (TH2D*) unfold.GetProbabilityMatrix(fVariableName1+"ProbabilityMatrix");

  for (Int_t cmj = 0; cmj < fNBinsY; cmj++)
    {
      FoldedBack_tempclone->SetBinError(cmj + 1,0);
    }

  cout<<"manual Chi2 (default): "<<FoldedBack_tempclone->Chi2Test(OrigInput,"WWP")<<endl;

  // Manually recalculate Chi2L
  TMatrixD m_L(fNBinsX, fNBinsX-2);
  TMatrixD m_LT(fNBinsX-2, fNBinsX);
  TMatrixD m_delta(1,fNBinsX);
  TMatrixD m_deltaT(fNBinsX,1);

  for (Int_t cmi = 0; cmi < fNBinsX; cmi++)
    {
      m_delta(0,cmi)  = TUnfResult->GetBinContent(cmi + 1) - BiasDist->GetBinContent(cmi + 1);
      m_deltaT(cmi,0) = TUnfResult->GetBinContent(cmi + 1) - BiasDist->GetBinContent(cmi + 1);

    for (Int_t cmj = 0; cmj < fNBinsX-2; cmj++)
      {
        m_L(cmi, cmj)  = LMatrix->GetBinContent(cmi + 1, cmj + 1);
        m_LT(cmj, cmi) = LMatrix->GetBinContent(cmi + 1, cmj + 1);
      }
    }

  TMatrixD LXsquared = m_delta*m_L*m_LT*m_deltaT;

  cout<<"manual Chi2L (default) excinctau2: "<<LXsquared(0,0)<<" "<<LXsquared(0,0)*unfold.GetTau()*unfold.GetTau()<<endl;

  TH1D *FoldedBack_rebinnedB  = (TH1D*) unfold_rebinnedB.GetFoldedOutput(fVariableName1+"InputFoldedBack_rebinnedB",0,"ttbarreco_rebinnedB","*[UO]", 0);
  TH1D *OrigInput_rebinnedB   = (TH1D*) unfold_rebinnedB.GetInput(fVariableName1+"OrigInput_rebinnedB",0,"ttbarreco_rebinnedB","*[UO]", 0);
  TH1D* BiasDist_rebinnedB    = (TH1D*) unfold_rebinnedB.GetBias(fVariableName1+"BiasDist_rebinnedB");
  TH1D* LxMinusBias_rebinnedB = (TH1D*) unfold_rebinnedB.GetLxMinusBias(fVariableName1+"LxMinusBias_rebinnedB");
  TH2D* LMatrix_rebinnedB     = (TH2D*) unfold_rebinnedB.GetL(fVariableName1+"LMatrix_rebinnedB");
  TH2D* ProbabilityMatrix_rebinnedB = (TH2D*) unfold_rebinnedB.GetProbabilityMatrix(fVariableName1+"ProbabilityMatrix_rebinnedB");  

  
  int nbinsrhoi = TUnfResult_rebinnedB->GetNbinsX();
  int nbinsfine = TUnfResult->GetNbinsX();
  int rebinfine = nbinsfine/nbinsrhoi;

  for (Int_t cmj = 0; cmj < nbinsrhoi*2; cmj++)
    {
      FoldedBack_rebinnedB->SetBinError(cmj + 1,0);
    }

  // cout<<"manual Chi2 (rebinnedbeforeunf): "<<FoldedBack_rebinnedB->Chi2Test(OrigInput_rebinnedB,"WWP")<<endl;


  // manually recalculate FoldedBack
  TMatrixD m_Probability_rebinnedB(nbinsrhoi, nbinsrhoi*2);
  TMatrixD m_Unfolded_rebinnedB(1,nbinsrhoi);
  for (Int_t cmi = 0; cmi < nbinsrhoi; cmi++)
    {
      m_Unfolded_rebinnedB(0,cmi) = TUnfResult_rebinnedB->GetBinContent(cmi + 1);
    for (Int_t cmj = 0; cmj < nbinsrhoi*2; cmj++)
      {
        m_Probability_rebinnedB(cmi, cmj) = ProbabilityMatrix_rebinnedB->GetBinContent(cmi + 1, cmj + 1);
      }
    }

  TMatrixD m_FoldedBack_rebinnedB = m_Unfolded_rebinnedB*m_Probability_rebinnedB;

  for (Int_t cmi = 0; cmi < nbinsrhoi*2; cmi++)
    {
      FoldedBack_rebinnedB->SetBinContent(cmi + 1,m_FoldedBack_rebinnedB(0,cmi));
      FoldedBack_rebinnedB->SetBinError(cmi + 1,0);
    }

  cout<<"manual Chi2 (rebinnedbeforeunf): "<<FoldedBack_rebinnedB->Chi2Test(OrigInput_rebinnedB,"WWP")<<endl;

  //manually recalculate Chi2L
  TMatrixD m_L_rebinnedB(nbinsrhoi, nbinsrhoi-2);
  TMatrixD m_LT_rebinnedB(nbinsrhoi-2, nbinsrhoi);
  TMatrixD m_delta_rebinnedB(1,nbinsrhoi);
  TMatrixD m_deltaT_rebinnedB(nbinsrhoi,1);
  for (Int_t cmi = 0; cmi < nbinsrhoi; cmi++)
    {
      m_delta_rebinnedB(0,cmi) = TUnfResult_rebinnedB->GetBinContent(cmi + 1) - BiasDist_rebinnedB->GetBinContent(cmi + 1);
      m_deltaT_rebinnedB(cmi,0) = TUnfResult_rebinnedB->GetBinContent(cmi + 1) - BiasDist_rebinnedB->GetBinContent(cmi + 1);
    for (Int_t cmj = 0; cmj < nbinsrhoi-2; cmj++)
      {
        m_L_rebinnedB(cmi, cmj) = LMatrix_rebinnedB->GetBinContent(cmi + 1, cmj + 1);
        m_LT_rebinnedB(cmj, cmi) = LMatrix_rebinnedB->GetBinContent(cmi + 1, cmj + 1);
      }
    }

  TMatrixD LXsquared_rebinnedB = m_delta_rebinnedB*m_L_rebinnedB*m_LT_rebinnedB*m_deltaT_rebinnedB;

  cout<<"manual Chi2L (rebinnedbeforeunf) excinctau2: "<<LXsquared_rebinnedB(0,0)<<" "<<LXsquared_rebinnedB(0,0)*unfold_rebinnedB.GetTau()*unfold_rebinnedB.GetTau()<<endl;



  // TUnfResult_rebinnedA foldedback Chi2
  TH1D *FoldedBack_rebinnedA = (TH1D*) OrigInput_rebinnedB->Clone(fVariableName1+"InputFoldedBack_rebinnedA");
  TMatrixD m_Unfolded_rebinnedA(1,nbinsrhoi);
  for (Int_t cmi = 0; cmi < nbinsrhoi; cmi++)
    {
      m_Unfolded_rebinnedA(0,cmi) = TUnfResult_rebinnedA->GetBinContent(cmi + 1);
    }

  TMatrixD m_FoldedBack_rebinnedA = m_Unfolded_rebinnedA*m_Probability_rebinnedB;

  for (Int_t cmi = 0; cmi < nbinsrhoi*2; cmi++)
    {
      FoldedBack_rebinnedA->SetBinContent(cmi + 1,m_FoldedBack_rebinnedA(0,cmi));
      FoldedBack_rebinnedA->SetBinError(cmi + 1,0);
    }

  cout<<"manual Chi2 (rebinnedafterunf): "<<FoldedBack_rebinnedA->Chi2Test(OrigInput_rebinnedB,"WWP")<<endl;

  //manually recalculate Chi2L
  TMatrixD m_delta_rebinnedA(1,nbinsrhoi);
  TMatrixD m_deltaT_rebinnedA(nbinsrhoi,1);
  for (Int_t cmi = 0; cmi < nbinsrhoi; cmi++)
    {
      m_delta_rebinnedA(0,cmi)  = TUnfResult_rebinnedA->GetBinContent(cmi + 1) - BiasDist_rebinnedB->GetBinContent(cmi + 1);
      m_deltaT_rebinnedA(cmi,0) = TUnfResult_rebinnedA->GetBinContent(cmi + 1) - BiasDist_rebinnedB->GetBinContent(cmi + 1);
    }

  TMatrixD LXsquared_rebinnedA = m_delta_rebinnedA*m_L_rebinnedB*m_LT_rebinnedB*m_deltaT_rebinnedA;

  cout<<"manual Chi2L (rebinnedafterunf) excinctau2: "<<LXsquared_rebinnedA(0,0)<<" "<<LXsquared_rebinnedA(0,0)*unfold_rebinnedB.GetTau()*unfold_rebinnedB.GetTau()<<endl;

  //cout<<"rrint: "<<TUnfResult_rebinnedA->Integral()/TUnfResult->Integral()<<" "<<TUnfResult_rebinnedA->Integral()/TUnfResult_rebinnedB->Integral()<<" "<<OrigInput->Integral()/OrigInput_rebinnedB->Integral()<<" "<<FoldedBack_rebinnedA->Integral()/FoldedBack_rebinnedB->Integral()<<" "<<FoldedBack_rebinnedB->Integral()/OrigInput_rebinnedB->Integral()<<" "<<FoldedBack_rebinnedA->Integral()/OrigInput->Integral()<<" "<<FoldedBack_tempclone->Integral()/OrigInput->Integral()<<endl;
  //TUnfResult_rebinnedA->Divide(TUnfResult_rebinnedB);
  //TUnfResult_rebinnedA->Print("all");


  TH1D* FoldedBack_tempclone_rebinnedB = (TH1D*) FoldedBack_rebinnedB->Clone();
  // TH1D* FoldedBack_tempclone = (TH1D*) FoldedBack->Clone();
  TH1D* FoldedBack_tempclone_rebinnedA = (TH1D*) FoldedBack_rebinnedA->Clone();
  TH1D* OrigInput_tempclone_rebinnedB  = (TH1D*) OrigInput_rebinnedB->Clone();
  TH1D* OrigInput_tempclone = (TH1D*) OrigInput->Clone();

  if (rebinfine > 1) {
    MatrixUnf::rebinMultidimensionalInput(FoldedBack_tempclone,
                                        detectorBinning_rebinnedB->FindNode("ttbarreco_rebinnedB"),
                                        detectorBinning->FindNode("ttbarreco"),
                                        frebinfine);
    MatrixUnf::rebinMultidimensionalInput(OrigInput_tempclone,
                                        detectorBinning_rebinnedB->FindNode("ttbarreco_rebinnedB"),
                                        detectorBinning->FindNode("ttbarreco"),
                                        frebinfine);
    cout << "manual Chi2 (rebinnedfoldedback): " << FoldedBack_tempclone->Chi2Test(OrigInput_tempclone, "WWP") << endl;
  }

  MatrixUnf::rebinMultidimensionalInput(FoldedBack_tempclone_rebinnedB,
                                      detectorBinning_rebinnedB->FindNode("ttbarreco_rebinnedB"),
                                      detectorBinning->FindNode("ttbarreco"),
                                      2);
  MatrixUnf::rebinMultidimensionalInput(FoldedBack_tempclone,
                                      detectorBinning_rebinnedB->FindNode("ttbarreco_rebinnedB"),
                                      detectorBinning->FindNode("ttbarreco"),
                                      2);
  MatrixUnf::rebinMultidimensionalInput(FoldedBack_tempclone_rebinnedA,
                                      detectorBinning_rebinnedB->FindNode("ttbarreco_rebinnedB"),
                                      detectorBinning->FindNode("ttbarreco"),
                                      2);
  MatrixUnf::rebinMultidimensionalInput(OrigInput_tempclone_rebinnedB,
                                      detectorBinning_rebinnedB->FindNode("ttbarreco_rebinnedB"),
                                      detectorBinning->FindNode("ttbarreco"),
                                      2);
  MatrixUnf::rebinMultidimensionalInput(OrigInput_tempclone,
                                      detectorBinning_rebinnedB->FindNode("ttbarreco_rebinnedB"),
                                      detectorBinning->FindNode("ttbarreco"),
                                      2);

  cout<<"manual Chi2 (rebinned2beforeunf): "<<FoldedBack_tempclone_rebinnedB->Chi2Test(OrigInput_tempclone_rebinnedB,"WWP")<<endl;
  cout<<"manual Chi2 (rebinned2foldedback): "<<FoldedBack_tempclone->Chi2Test(OrigInput_tempclone,"WWP")<<endl;
  cout<<"manual Chi2 (rebinned2afterunf): "<<FoldedBack_tempclone_rebinnedA->Chi2Test(OrigInput_tempclone_rebinnedB,"WWP")<<endl;

  delete FoldedBack_tempclone_rebinnedB;
  delete FoldedBack_tempclone;
  delete FoldedBack_tempclone_rebinnedA;
  delete OrigInput_tempclone_rebinnedB;
  delete OrigInput_tempclone;
}


void MatrixUnf::Chi2testsForPEs(TUnfoldDensity &unfold, TUnfoldDensity &unfold_rebinnedB, TH1D* TUnfResult_rebinnedA, std::vector<double> &chi2unf){

  chi2unf.clear();

  TH1D *FoldedBack_tempclone           = (TH1D*) unfold.GetFoldedOutput(fVariableName1+"InputFoldedBack_tempclone",0,"ttbarreco","*[UO]",0);
  TH1D *FoldedBack_tempclone_rebinnedB = (TH1D*) unfold_rebinnedB.GetFoldedOutput(fVariableName1+"InputFoldedBack_tempclone_rebinnedB",0,"ttbarreco_rebinnedB","*[UO]",0);
  TH1D *OrigInput_tempclone_rebinnedB  = (TH1D*) unfold_rebinnedB.GetInput(fVariableName1+"OrigInput_tempclone_rebinnedB",0,"ttbarreco_rebinnedB","*[UO]",0);

  TH2D* ProbabilityMatrix_tempclone_rebinnedB = (TH2D*) unfold_rebinnedB.GetProbabilityMatrix(fVariableName1+"ProbabilityMatrix_tempclone_rebinnedB");  
 
  int nbinsrhoi = TUnfResult_rebinnedA->GetNbinsX();
  int rebinfine = fNBinsX/nbinsrhoi;

  for (Int_t cmj = 0; cmj < fNBinsY; cmj++)
    {
      FoldedBack_tempclone->SetBinError(cmj + 1,0);
    }

  for (Int_t cmj = 0; cmj < nbinsrhoi*2; cmj++)
    {
      FoldedBack_tempclone_rebinnedB->SetBinError(cmj + 1,0);
    }

  //rebinned probability matrix
  TMatrixD m_Probability_rebinnedB(nbinsrhoi, nbinsrhoi*2);
  for (Int_t cmi = 0; cmi < nbinsrhoi; cmi++)
    {
    for (Int_t cmj = 0; cmj < nbinsrhoi*2; cmj++)
      {
        m_Probability_rebinnedB(cmi, cmj) = ProbabilityMatrix_tempclone_rebinnedB->GetBinContent(cmi + 1, cmj + 1);
      }
    }

  // TUnfResult_rebinnedA foldedback 
  TH1D *FoldedBack_tempclone_rebinnedA = (TH1D*) OrigInput_tempclone_rebinnedB->Clone(fVariableName1+"InputFoldedBack_tempclone_rebinnedA");
  TMatrixD m_Unfolded_rebinnedA(1,nbinsrhoi);
  for (Int_t cmi = 0; cmi < nbinsrhoi; cmi++)
    {
      m_Unfolded_rebinnedA(0,cmi) = TUnfResult_rebinnedA->GetBinContent(cmi + 1);
    }

  TMatrixD m_FoldedBack_rebinnedA = m_Unfolded_rebinnedA*m_Probability_rebinnedB;

  for (Int_t cmi = 0; cmi < nbinsrhoi*2; cmi++)
    {
      FoldedBack_tempclone_rebinnedA->SetBinContent(cmi + 1,m_FoldedBack_rebinnedA(0,cmi));
      FoldedBack_tempclone_rebinnedA->SetBinError(cmi + 1,0);
    }

  if (frebinfine > 1) {
    MatrixUnf::rebinMultidimensionalInput(FoldedBack_tempclone,
                                        detectorBinning_rebinnedB->FindNode("ttbarreco_rebinnedB"),
                                        detectorBinning->FindNode("ttbarreco"),
                                        frebinfine);
  }

  double chi2unfRB = FoldedBack_tempclone->Chi2Test(OrigInput_tempclone_rebinnedB,"WW CHI2");
  double chi2unfRB_rebinnedA = FoldedBack_tempclone_rebinnedA->Chi2Test(OrigInput_tempclone_rebinnedB,"WW CHI2");
  double chi2unfRB_rebinnedB = FoldedBack_tempclone_rebinnedB->Chi2Test(OrigInput_tempclone_rebinnedB,"WW CHI2");

  //now rebin everything to nbinsrhoi bins for a "fair" comparison of the Chi2
// now rebin everything to nbinsrhoi bins for a "fair" comparison of the Chi2
  MatrixUnf::rebinMultidimensionalInput(FoldedBack_tempclone_rebinnedB,
                                      detectorBinning_rebinnedB->FindNode("ttbarreco_rebinnedB"),
                                      detectorBinning->FindNode("ttbarreco"),
                                      2);
  MatrixUnf::rebinMultidimensionalInput(FoldedBack_tempclone,
                                      detectorBinning_rebinnedB->FindNode("ttbarreco_rebinnedB"),
                                      detectorBinning->FindNode("ttbarreco"),
                                      2);
  MatrixUnf::rebinMultidimensionalInput(FoldedBack_tempclone_rebinnedA,
                                      detectorBinning_rebinnedB->FindNode("ttbarreco_rebinnedB"),
                                      detectorBinning->FindNode("ttbarreco"),
                                      2);
  MatrixUnf::rebinMultidimensionalInput(OrigInput_tempclone_rebinnedB,
                                      detectorBinning_rebinnedB->FindNode("ttbarreco_rebinnedB"),
                                      detectorBinning->FindNode("ttbarreco"),
                                      2);

  double chi2unfRB2 = FoldedBack_tempclone->Chi2Test(OrigInput_tempclone_rebinnedB,"WW CHI2");
  double chi2unfRB2_rebinnedA = FoldedBack_tempclone_rebinnedA->Chi2Test(OrigInput_tempclone_rebinnedB,"WW CHI2");
  double chi2unfRB2_rebinnedB = FoldedBack_tempclone_rebinnedB->Chi2Test(OrigInput_tempclone_rebinnedB,"WW CHI2");

  std::vector<double> chi2unf_temp {chi2unfRB, chi2unfRB_rebinnedA, chi2unfRB_rebinnedB, chi2unfRB2, chi2unfRB2_rebinnedA, chi2unfRB2_rebinnedB};

  chi2unf = chi2unf_temp;

  //cout<<"manual Chi2 (rebinned2beforeunf): "<<chi2unfRB_rebinnedB<<endl;
  //cout<<"manual Chi2 (rebinned2foldedback): "<<chi2unfRB<<endl;
  //cout<<"manual Chi2 (rebinned2afterunf): "<<chi2unfRB_rebinnedA<<endl;

  delete FoldedBack_tempclone_rebinnedB;
  delete FoldedBack_tempclone;
  delete FoldedBack_tempclone_rebinnedA;
  delete OrigInput_tempclone_rebinnedB;
  delete ProbabilityMatrix_tempclone_rebinnedB;

}


/**
 * @brief Fills the histograms AfbResult and OtherResult with information regarding the processed unfolded results.
 * 
 * @param TUnfResult The unfolded result that is going to be processed
 * @param Ematrix The covariance matrix returned from the unfolding procedure
 * @param GenVecHist The gen-level truth information that was used to construct the response matrix
 * @param VariableName1 The variable that is being unfolded
 * @param AfbResult The histogram used to store the results from the unfolding. The general structure repeats numCoefficients times.
 *  bin ith_variable * 8 = a_optimal_reco[0] +/- a_optimal_reco[0]
 *  bin ith_variable * 8 + 1 = Afbvec[ith_variable * 3] +/- AfbErrvec[ith_variable * 3] --> asymmetry measured in first and last bin
 *  bin ith_variable * 8 + 2 = Afbvec[ith_variable * 3 + 1] +/- AfbErrvec[ith_variable * 3 + 1] --> asymmetry measured in second and second-to-last bin
 *  bin ith_variable * 8 + 3 = Afbvec[ith_variable * 3 + 2] +/- AfbErrvec[ith_variable * 3 + 2] --> asymmetry measured in third and third-to-last bin
 *  bin ith_variable * 8 + 4 = Percent difference of unfolded result to generator level truth with uncertainty
 *  bin ith_variable * 8 + 5 = Ceof +/- CoefErr --> coefficient converted from inclusive asymmetry
 *  bin ith_variable * 8 + 6 = Afb +/- AfbErr --> inclusive asymmetry
 *  bin ith_variable * 8 + 7 = a_optimal_reco[1] +/- a_optimal_reco[1]
 * @param OtherResult The histogram used to store information about the hyperparameters of the unfolding and the chi^2 values
 *  bin 0 = log10(optimal tau)
 *  bin ith_variable * 2 + 1 = m_Afb(0, 1) for ith variable
 *  bin ith_variable * 2 + 2 = (m_Afb(0, 1) + m_Afb(0, 2) + m_Afb(1, 2)) / 3 for ith variable
 *  bin numCoefficients * 2 + 1 = rho average at optimal tau
 *  bin numCoefficients * 2 + 2 = a_optimal_reco[2]  TODO: This is probably incorrect currently...
 *  bin numCoefficients * 2 + 3 = unfolding chi^2 = chi^2_A + chi^2_L
 *  bin numCoefficients * 2 + 4 = chi2unfRB2 +/- chi2unfRB2_rebinnedA
 *  bin numCoefficients * 2 + 5 = chi2unfRB +/- chi2unfRB_rebinnedA
 * @param ensTest Whether or not this is part of an ensemble test (using pseudoexperiments)
 * @param isgenlevel Whether this unfolded result was from a gen-level input
 * @param numCoefficients The number of variables we will be measuring from this distribution
 */
void MatrixUnf::fillAsymsAndCoefficient(TH1D* TUnfResult,
                                        TH2D* Ematrix,
                                        TH1D* GenVecHist,
                                        TString VariableName1,
                                        TH1D* AfbResult,
                                        TH1D* OtherResult,
                                        bool ensTest,
                                        bool isgenlevel,
                                        int numCoefficients) {
  AfbResult->Reset();
  if (OtherResult)
    OtherResult->Reset();
  int nbinsX = TUnfResult->GetNbinsX();
  int numBinsPerCoefficient = nbinsX / numCoefficients;
  int numAsymmetries = numBinsPerCoefficient / 2;

  Double_t Afb[numCoefficients], AfbErr[numCoefficients], Coef[numCoefficients], CoefErrOpt[numCoefficients], CoefErrTest[numCoefficients];
  std::vector<double> Afbvec, AfbErrvec, Coefvec, a_optimal_reco;
  TMatrixD m_C(nbinsX / 2, nbinsX / 2);

  double binsumerr[numCoefficients];
  double binsum[numCoefficients];
  for (int i = 0; i < numCoefficients; i++) {
    binsumerr[i] = 0.0;
    binsum[i] = 0.0;
  }

  // covariance matrix filling
  TMatrixD m_E(nbinsX, nbinsX);
  if (!isgenlevel) {
    for (int i = 0; i < numCoefficients; i++) {
      for (Int_t cmi = 0; cmi < numBinsPerCoefficient; cmi++) {
        int bin_i = i * numBinsPerCoefficient + cmi;
        for (Int_t cmj = 0; cmj < numBinsPerCoefficient; cmj++) {
          int bin_j = i * numBinsPerCoefficient + cmj;
          m_E(bin_i, bin_j) = Ematrix->GetBinContent(bin_i + 1, bin_j + 1);
          binsumerr[i] += Ematrix->GetBinContent(bin_i + 1, bin_j + 1);
        }
        binsum[i] += TUnfResult->GetBinContent(bin_i + 1);
      }
      binsumerr[i] = sqrt(binsumerr[i]);
      binsum[i] /= GenVecHist->Integral(i * numBinsPerCoefficient + 1, (i + 1) * numBinsPerCoefficient);
      binsumerr[i] /= GenVecHist->Integral(i * numBinsPerCoefficient + 1, (i + 1) * numBinsPerCoefficient);
    }
  } else {
    for (int i = 0; i < numCoefficients; i++) {
      for (Int_t cmi = 0; cmi < numBinsPerCoefficient; ++cmi) {
        int bin_i = i * numBinsPerCoefficient + cmi;
        for (Int_t cmj = 0; cmj < numBinsPerCoefficient; cmj++) {
          int bin_j = i * numBinsPerCoefficient + cmj;
          if (cmi == cmj)
            m_E(bin_i, bin_j) = pow(TUnfResult->GetBinError(bin_i + 1), 2);
          else
            m_E(bin_i, bin_j) = 0;
        }
      }
      binsum[i] = 1.;
      binsumerr[i] = 0.;
    }
  }

  MatrixUnf::GetAfbWithCorrelations(TUnfResult, m_E, Afb, AfbErr, numCoefficients);
  TMatrixD m_AFB(nbinsX / 2, nbinsX / 2);
  MatrixUnf::GetAfbWithCorrelationsBinByBin(TUnfResult, m_E, Afbvec, AfbErrvec, m_AFB, numCoefficients);
  bool isCoefficient = MatrixUnf::GetOptimalCoefficient(Afbvec, m_AFB, GenVecHist, VariableName1, Coef, CoefErrOpt, CoefErrTest, Coefvec, m_C, a_optimal_reco, numCoefficients);

  // if the variable is not a spin density matrix variable, fill the asymmetry of the distribution instead of the calculated spin coefficient
  if (!isCoefficient) {
    for (int i = 0; i < numCoefficients; i++) {
      int binOffset = i * numBinsPerCoefficient;
      int asymOffset = i * numAsymmetries;
      Coef[i] = Afb[i];
      CoefErrTest[i] = AfbErr[i];
      CoefErrOpt[i] = AfbErr[i];
      // to reconstruct the inclusive asymmetry from the binned asymmetries, the appropriate factors are the fractions of events in each bin
      double integral = TUnfResult->Integral(i * numBinsPerCoefficient + 1, (i + 1) * numBinsPerCoefficient);
      for (int j = 0; j < numAsymmetries; ++j) {
        a_optimal_reco[asymOffset + j] = (TUnfResult->GetBinContent(binOffset + j + 1) + TUnfResult->GetBinContent(binOffset + numBinsPerCoefficient - j)) / integral;
      }
    }
  }

  for (int i = 0; i < numCoefficients; i++) {
    int binOffset = i * (numBinsPerCoefficient + 2);
    int AsymOffset = i * numAsymmetries;

    // fill a_optimal_reco[0]
    AfbResult->SetBinContent(binOffset, a_optimal_reco.at(AsymOffset));
    AfbResult->SetBinError(binOffset, a_optimal_reco.at(AsymOffset));

    // fill the asymmetries
    for (int j = 0; j < numAsymmetries; ++j) {
      AfbResult->SetBinContent(binOffset + j + 1, Afbvec[AsymOffset + j]);
      AfbResult->SetBinError(binOffset + j + 1, AfbErrvec[AsymOffset + j]);
      if (j > 0)
        printf(" AFB%s%i: %.6f+-%.6f ", isgenlevel ? "gen" : "", j, Afbvec[AsymOffset + j], AfbErrvec[AsymOffset + j]);
    }

    // percent difference of unfolded result to generator level truth with uncertainty
    AfbResult->SetBinContent(binOffset + numAsymmetries + 1, binsum[i] - 1.);
    AfbResult->SetBinError(binOffset + numAsymmetries + 1, binsumerr[i]);

    // Ceof +/- CoefErr --> coefficient converted from inclusive asymmetry
    AfbResult->SetBinContent(binOffset + numAsymmetries + 2, Coef[i]);
    AfbResult->SetBinError(binOffset + numAsymmetries + 2, CoefErrOpt[i]);

    // Afb +/- AfbErr --> inclusive asymmetry
    AfbResult->SetBinContent(binOffset + numAsymmetries + 3, Afb[i]);
    AfbResult->SetBinError(binOffset + numAsymmetries + 3, AfbErr[i]);

    // a_optimal_reco[1] +/- a_optimal_reco[1]
    AfbResult->SetBinContent(binOffset + numAsymmetries + 4, a_optimal_reco.at(AsymOffset + 1));
    AfbResult->SetBinError(binOffset + numAsymmetries + 4, a_optimal_reco.at(AsymOffset + 1));
    
    printf("MatrixUnf::prepRunUnfolding: %s: AFB%s: %.6f+-%.6f \n",
          VariableName1.Data(),
          isgenlevel ? "gen" : "",
          Afb[i],
          AfbErr[i]);
    printf("MatrixUnf::prepRunUnfolding: %s: AFB%s0: %.6f+-%.6f ",
          VariableName1.Data(),
          isgenlevel ? "gen" : "",
          Afbvec[AsymOffset],
          AfbErrvec[AsymOffset]);
    printf("\n");
    printf("MatrixUnf::prepRunUnfolding: %s: Coef%s: %.6f+-%.6f \n",
          VariableName1.Data(),
          isgenlevel ? "gen" : "",
          Coef[i],
          CoefErrOpt[i]);    

    // TODO: This might be incorrect since we should only be normalizing by part of the diagonal for each submatrix
    m_AFB.NormByDiag(TMatrixDDiag(m_AFB));

    if (OtherResult) {
      // fill rough measure of the average correlation between AFB bins (including sign, unlike rhoavg)
      OtherResult->SetBinContent(i * 2 + 1, m_AFB(AsymOffset, AsymOffset + 1));
      if (nbinsX >= 6) {
        OtherResult->SetBinContent(i * 2 + 2, (m_AFB(AsymOffset, AsymOffset + 1) + m_AFB(AsymOffset, AsymOffset + 2) + m_AFB(AsymOffset + 1, AsymOffset + 2)) / 3.);
        OtherResult->SetBinContent(numCoefficients * 2 + 2, a_optimal_reco.at(2)); // TODO: Needs fixed to grab all a_optimal_recos
      }
    }

    if (!ensTest && !isgenlevel) {
      m_AFB.Print("f=%1.5g ");

      m_E.NormByDiag(TMatrixDDiag(m_E));
      m_E.Print("f=%1.5g ");
    }
  }
}

/**
 * @brief Calculates and stores the bin factor function in the BFFVec parameter. The bin factor function is to ensure that the measured distribution 
 * only varies linearly in changes to the observed spin coefficient. In general, the distribution may vary strongly non-linearly on the spin coefficient.
 * By dividing out the non-linearity of the dependency we can create the regularized distribution to only vary linearly. We know the functional form and are
 * therefore able to do this.
 * 
 * @param GenVecHist_forBias A histogram which stores the GEN information of the distribution we will be unfolding.
 * @param BFFvec The vector that will be used in the regularization to store linearity in the distribution
 * @param binning The TUnfoldBinning object that stores the relevant binning information for the distribution
 */
void MatrixUnf::calculateBFF(TH1D* GenVecHist_forBias, TVectorD* BFFvec, TUnfoldBinning *binning) {
  // grab the spin coefficient name
  TString coefficientName = binning->GetDistributionAxisLabel(0);
  std::cout << "calculating bin factor function for " << coefficientName << std::endl;
  
  // calculate the number of bins along spin coefficient distribution and other dimensions
  int nbinsX  = 0;
  int numBins = 1;

  for (int i = 0; i < binning->GetDistributionDimension(); ++i) {
    const TVectorD* vec_binEdges = binning->GetDistributionBinning(i);
    if (i == 0) { // 0th distribution is the spin correlation variable
      nbinsX = vec_binEdges->GetNoElements() - 1;
    } 
    else { // number of bins of the other dimensions
      numBins *= vec_binEdges->GetNoElements() - 1;
    }
  }

  // copy bin edges because they are const doubles currently which ruins everything
  double xBinEdges[nbinsX + 1];
  for (int i = 0; i < nbinsX + 1; ++i) {
    xBinEdges[i] = binning->GetDistributionBinning(0)->GetMatrixArray()[i];
  }

  std::cout << "Bin factor function = [";
  // compute and store the BFF
  for (int i = 0; i < numBins; i++) {
    int binOffset = i * nbinsX;

    // analytic functions for BinFactorFunctions
    // and calculate oneOverSMdensity_sym analytically
    string symmetrizedSMDensityFunction;
    if      (MatrixUnf::isCtype(coefficientName)  ) symmetrizedSMDensityFunction = "-log(abs(x))";
    else if (MatrixUnf::isCPMtype(coefficientName)) symmetrizedSMDensityFunction = "acos(abs(x))";
    else if (MatrixUnf::isBPMtype(coefficientName)) symmetrizedSMDensityFunction = "(2 - abs(x))"; // note, the factor here depends on whether we're using the distribution to measure a C or B coefficient. Currently set for B.
    else if (MatrixUnf::isBtype(coefficientName) || MatrixUnf::isDtype(coefficientName) || MatrixUnf::isLinearCombtype(coefficientName)) symmetrizedSMDensityFunction = "1"; // FIXME : New variables are also flat 
    else symmetrizedSMDensityFunction = "1";
    // else if( MatrixUnf::isBPMtype(coefficientName) )  SMdensity_sym = new TF1("SMdensity_sym","( (2 - abs(x)) * (1 + sqrt(3) + abs(x)) )/(sqrt(3) + 5/3)", GenVecHist_forBias->GetXaxis()->GetBinLowEdge(1), GenVecHist_forBias->GetXaxis()->GetBinUpEdge(nbinsX) ); //note, the factor here depends on whether we're using the distribution to measure a C or B coefficient. Currently set for C, where the factor makes the function linear in Abs(x), i.e. the curvature is unavoidably non-zero between the middle two bins.
    
    TF1 *SMdensity_sym    = new TF1("SMdensity_sym", symmetrizedSMDensityFunction.c_str(), xBinEdges[0], xBinEdges[nbinsX]);

    // 1/(symmetrised gen-level distribution), oneOverSMdensity_sym, is the same as the analytic functions (above) 
    // for all the spin density matrix variables except for b_P/M*
    // but the values are slightly more precise 
    // (because the function is naturally integrated over the bin width rather than taking the value at the bin centre)

    TVectorD oneOverSMdensity(nbinsX);
    TVectorD oneOverSMdensity_sym(nbinsX);
    
    double gen_sum = 0.;
    
    for (int j = 0; j < nbinsX; ++j) {
      double jthBinWidth    = xBinEdges[j + 1] - xBinEdges[j];
      double genyield       = GenVecHist_forBias->GetBinContent(binOffset + j + 1);
      double gendensity     = GenVecHist_forBias->GetBinContent(binOffset + j + 1) / jthBinWidth;
      double gendensity_sym = 
          (GenVecHist_forBias->GetBinContent(binOffset + j + 1) / jthBinWidth
        + GenVecHist_forBias->GetBinContent(binOffset + nbinsX - j) / (xBinEdges[nbinsX - j] - xBinEdges[nbinsX - j - 1]))
        / 2.;
      oneOverSMdensity(j)     = 1. / gendensity;
      oneOverSMdensity_sym(j) = 1. / gendensity_sym;
      gen_sum += genyield;
    }

    oneOverSMdensity     *= gen_sum / (xBinEdges[nbinsX] - xBinEdges[0]);
    oneOverSMdensity_sym *= gen_sum / (xBinEdges[nbinsX] - xBinEdges[0]);
    
    SMdensity_sym->SetNpx(10000 * nbinsX);
    TH1D* SMdensity_sym_histo = (TH1D*)SMdensity_sym->GetHistogram();
    SMdensity_sym_histo->Rebin(10000);
    SMdensity_sym_histo->Scale(0.0001);

    // Amandeep : Confirm with AJ MatrixUnf::isLinearCombtype() should be added here ?
    if (MatrixUnf::isCtype(coefficientName)   || MatrixUnf::isCPMtype(coefficientName) ||
        MatrixUnf::isBtype(coefficientName)   || MatrixUnf::isDtype(coefficientName)   ||
        MatrixUnf::isBPMtype(coefficientName) || MatrixUnf::isLinearCombtype(coefficientName)) {
          for (int j = 0; j < nbinsX; ++j) {
            if (SMdensity_sym_histo->GetBinContent(j + 1) == 0) {
              std::cout << "Error, bin content = 0 and trying to divide by 0." << std::endl;
            }
            oneOverSMdensity_sym(j) = 1. / SMdensity_sym_histo->GetBinContent(j + 1);
          }
    }

    // Amandeep : Confirm with AJ MatrixUnf::isLinearCombtype() should be added here ?
    if (MatrixUnf::isCtype(coefficientName)   || MatrixUnf::isCPMtype(coefficientName) ||
        MatrixUnf::isBtype(coefficientName)   || MatrixUnf::isDtype(coefficientName)   ||
        MatrixUnf::isBPMtype(coefficientName) || MatrixUnf::isOtherSymType(coefficientName) || 
        MatrixUnf::isLinearCombtype(coefficientName)) {
          for (int j = 0; j < nbinsX; ++j) {
            (*BFFvec)(binOffset + j) = oneOverSMdensity_sym(j);
            std::cout << " " << oneOverSMdensity_sym(j) << ", ";
            // std::cout << ""  << SMdensity_sym_histo->GetBinContent(j + 1) << " | ";
          }
    } 

    else {
      for (int j = 0; j < nbinsX; ++j) {
        // use unsymmetrised SM for the non-spin variables, where the shape is not motivated by the shape of the NP 
        // but gives unbiased regularisation for the case dist(x)->(1+const*x)*dist(x)
        (*BFFvec)(binOffset + j) = oneOverSMdensity(j); 
        std::cout << " " << oneOverSMdensity(j);
      }
    }
  }
  std::cout << " ]" << std::endl;
}

/**
 * @brief TODO: Fix this for multidimensional unfolding
 * 
 * @param GenVecHist_forBias The distribution to be reweighted to test for bias in unfolding
 */
void MatrixUnf::ReweightBiasTest(TH1D* GenVecHist_forBias) {
    double BiasDistIntegral = GenVecHist_forBias->Integral();

    for (int i = 1; i <= GenVecHist_forBias->GetNbinsX(); ++i) {
      double x = ((i - GenVecHist_forBias->GetNbinsX() / 2) - 0.5) / (GenVecHist_forBias->GetNbinsX() / 2);
      double weightbias = 1. + 0.5 * x * x; 
      GenVecHist_forBias->SetBinContent(i, weightbias * GenVecHist_forBias->GetBinContent(i));
      GenVecHist_forBias->SetBinError(i, weightbias * GenVecHist_forBias->GetBinError(i));
    }
    
    GenVecHist_forBias->Scale(BiasDistIntegral / GenVecHist_forBias->Integral());
}


/**
 * @brief Rebins a 1D histogram that is an unraveled 2D+ histogram correctly based on it originally being subdivided
 *        numSubdivisions times.
 * 
 * @param unraveledHistogram The 1D histogram that is a 2D+ distribution that has been unraveled into 1D
 * @param rebinnedBinning The binning object that contains the binning for the rebinned distribution
 * @param originalBinning The binning object that contains the original fine binning
 * @param rebinFactor The amount to rebin by
 */
void MatrixUnf::rebinMultidimensionalInput(TH1* unraveledHistogram,
                                         const TUnfoldBinning* rebinnedBinning,
                                         const TUnfoldBinning* originalBinning,
                                         int rebinFactor) {
  std::cout << "Base rebin factor of " << rebinFactor << std::endl;
  int numBins = 1;
  int totalRebinFactor = 1;

  // totalRebinFactor = rebinFactor^numDimensions
  // where numDimensions is the number of axis the unfolding is taking place on
  for (int i = 0; i < originalBinning->GetDistributionDimension(); ++i) {
    const TVectorD* vec_binEdges = originalBinning->GetDistributionBinning(i);
    numBins *= vec_binEdges->GetNoElements() - 1;
    totalRebinFactor *= rebinFactor;
  }

  int binMap[numBins + 2];
  MatrixUnf::CreateBinMap(originalBinning, rebinFactor, binMap);

  TH1* rebinnedHistogram = rebinnedBinning->CreateHistogram(Form("%s_rebinnedB", unraveledHistogram->GetName()));

  for (int i = 0; i < unraveledHistogram->GetNbinsX() + 1; ++i) {
    int binIndex = binMap[i];
    if (binIndex == -1) {
      continue;
    }
    double currentBinContent = rebinnedHistogram->GetBinContent(binIndex);
    rebinnedHistogram->SetBinContent(binIndex, currentBinContent + unraveledHistogram->GetBinContent(i));

    double currentBinError = pow(rebinnedHistogram->GetBinError(binIndex), 2.0);
    rebinnedHistogram->SetBinError(binIndex, sqrt(currentBinError + pow(unraveledHistogram->GetBinError(i), 2.0)));
  }

  // copy contents over
  std::cout << "Rebinning the unraveled histogram " << unraveledHistogram->GetName() << " with originally "
            << unraveledHistogram->GetNbinsX() << " bins by factor " << totalRebinFactor << std::endl;
  unraveledHistogram->Rebin(totalRebinFactor);
  std::cout << "Finished rebinning the unraveled histogram" << std::endl;
  unraveledHistogram->Reset();
  for (int i = 0; i < unraveledHistogram->GetNbinsX() + 1; ++i) {
    unraveledHistogram->SetBinContent(i, rebinnedHistogram->GetBinContent(i));
    unraveledHistogram->SetBinError(i, rebinnedHistogram->GetBinError(i));
  }
  delete rebinnedHistogram;
}

/**
 * @brief Rebins a 2D histogram that is an unraveled 2D+ histogram correctly based on it originally being subdivided
 *        numSubdivisions times along each axis
 * 
 * @param unraveledMigrationMatrix The unraveled 2D histogram that has 1D+ distributions along each x and y axis
 * @param generatorRebinnedBinning The binning with the rebinning applied along the generator axis
 * @param detectorRebinnedBinning The binning with the rebinning applied along the detector axis
 * @param generatorBinning The original fine binning along the generator axis
 * @param detectorBinning The original fine binning along the detector axis
 * @param rebinFactor The amount to rebin by
 */
void MatrixUnf::rebinMultidimMigrationMatrix(TH2* unraveledMigrationMatrix,
                                             const TUnfoldBinning* generatorRebinnedBinning,
                                             const TUnfoldBinning* detectorRebinnedBinning,
                                             const TUnfoldBinning* generatorBinning,
                                             const TUnfoldBinning* detectorBinning,
                                             int rebinFactor) {
  TH2* rebinnedMigrationMatrix = TUnfoldBinning::CreateHistogramOfMigrations(
      generatorRebinnedBinning, detectorRebinnedBinning, Form("%s_rebinnedB", unraveledMigrationMatrix->GetName()));

  int numGenBins       = 1;
  int numDetectorBins  = 1;
  int totalRebinFactor = 1;
  
  // Amandeep : Ask AJ if this is correct ?
  for (int i = 0; i < generatorBinning->GetDistributionDimension(); ++i) {
    const TVectorD* vec_binEdges = generatorBinning->GetDistributionBinning(i);
    numGenBins       *= vec_binEdges->GetNoElements() - 1;
    totalRebinFactor *= rebinFactor;
  }

  for (int i = 0; i < detectorBinning->GetDistributionDimension(); ++i) {
    const TVectorD* vec_binEdges = detectorBinning->GetDistributionBinning(i);
    numDetectorBins *= vec_binEdges->GetNoElements() - 1;
  }

  int genBinMap[numGenBins + 2], detectorBinMap[numDetectorBins + 2];
  MatrixUnf::CreateBinMap(generatorBinning, rebinFactor, genBinMap);
  MatrixUnf::CreateBinMap(detectorBinning, rebinFactor, detectorBinMap);

  for (int i = 0; i < unraveledMigrationMatrix->GetNbinsX() + 1; ++i) {
    int genBinIndex = genBinMap[i];
    if (genBinIndex == -1) {
      continue;
    }
    for (int j = 0; j < unraveledMigrationMatrix->GetNbinsY() + 1; ++j) {
      int detectorBinIndex = detectorBinMap[j];
      if (detectorBinIndex == -1) {
        continue;
      }

      double currentBinContent = rebinnedMigrationMatrix->GetBinContent(genBinIndex, detectorBinIndex);
      rebinnedMigrationMatrix->SetBinContent(
          genBinIndex, detectorBinIndex, currentBinContent + unraveledMigrationMatrix->GetBinContent(i, j));

      double currentBinError = pow(rebinnedMigrationMatrix->GetBinError(genBinIndex, detectorBinIndex), 2.0);
      rebinnedMigrationMatrix->SetBinError(
          genBinIndex, detectorBinIndex, sqrt(currentBinError + pow(unraveledMigrationMatrix->GetBinError(i, j), 2.0)));
    }
  }

  // copy contents over
  unraveledMigrationMatrix->Rebin2D(totalRebinFactor, totalRebinFactor);
  unraveledMigrationMatrix->Reset();
  for (int i = 0; i < unraveledMigrationMatrix->GetNbinsX() + 1; ++i) {
    for (int j = 0; j < unraveledMigrationMatrix->GetNbinsY() + 1; ++j) {
      unraveledMigrationMatrix->SetBinContent(i, j, rebinnedMigrationMatrix->GetBinContent(i, j));
      unraveledMigrationMatrix->SetBinError(i, j, rebinnedMigrationMatrix->GetBinError(i, j));
    }
  }
  delete rebinnedMigrationMatrix;
}

/**
 * @brief Rebins a 2D histogram that is an unraveled 2D+ histogram along a single axis correctly based on it originally being subdivided
 *        numSubdivisions times along a single axis
 * 
 * @param unraveledMigrationMatrix The unraveled 2D histogram that has 1D+ distributions along each x and y axis
 * @param rebinnedBinning The binning with the rebinning applied
 * @param binning The original fine binning
 * @param rebinFactor The amount to rebin by
 * @param axis The axis in which to rebin along. 0 -> x-axis, 1 -> y-axis
 */
void MatrixUnf::rebinMultidimMigrationMatrix(TH2* unraveledMigrationMatrix,
                                             const TUnfoldBinning* rebinnedBinning,
                                             const TUnfoldBinning* binning,
                                             int rebinFactor,
                                             int axis) {
  std::cout << "Base rebin factor of " << rebinFactor << std::endl;
  if (axis != 0 && axis != 1) {
    std::cout << "ERROR: Unspecified axis to rebin along for N-D migration matrix!!" << std::endl;
    return;
  }
  
  TH2* rebinnedMigrationMatrix = TUnfoldBinning::CreateHistogramOfMigrations(
      rebinnedBinning, rebinnedBinning, Form("%s_rebinned", unraveledMigrationMatrix->GetName()));

  int numBins = 1;
  int totalRebinFactor = 1;
  for (int i = 0; i < binning->GetDistributionDimension(); ++i) {
    const TVectorD* vec_binEdges = binning->GetDistributionBinning(i);
    numBins *= vec_binEdges->GetNoElements() - 1;
    totalRebinFactor *= rebinFactor;
  }

  int binMap[numBins + 2];
  MatrixUnf::CreateBinMap(binning, rebinFactor, binMap);
  for (int i = 0; i <= unraveledMigrationMatrix->GetNbinsX() + 1; ++i) {
    for (int j = 0; j <= unraveledMigrationMatrix->GetNbinsY() + 1; ++j) {
      if (axis == 0) {
        int binIndex = binMap[i];
        if (binIndex == -1) {
          binIndex = 0;
        }
        double currentBinContent = rebinnedMigrationMatrix->GetBinContent(binIndex, j);
        rebinnedMigrationMatrix->SetBinContent(
            binIndex, j, currentBinContent + unraveledMigrationMatrix->GetBinContent(i, j));

        double currentBinError = pow(rebinnedMigrationMatrix->GetBinError(binIndex, j), 2.0);
        rebinnedMigrationMatrix->SetBinError(
            binIndex, j, sqrt(currentBinError + pow(unraveledMigrationMatrix->GetBinError(i, j), 2.0)));
      } else if (axis == 1) {
        int binIndex = binMap[j];
        if (binIndex == -1) {
          binIndex = 0;
        }
        double currentBinContent = rebinnedMigrationMatrix->GetBinContent(i, binIndex);
        rebinnedMigrationMatrix->SetBinContent(
            i, binIndex, currentBinContent + unraveledMigrationMatrix->GetBinContent(i, j));

        double currentBinError = pow(rebinnedMigrationMatrix->GetBinError(i, binIndex), 2.0);
        rebinnedMigrationMatrix->SetBinError(
            i, binIndex, sqrt(currentBinError + pow(unraveledMigrationMatrix->GetBinError(i, j), 2.0)));
      }
    }
  }

  // copy contents over
  if (axis == 0) {
    std::cout << "Rebinning the unraveled histogram " << unraveledMigrationMatrix->GetName() << " with originally "
            << unraveledMigrationMatrix->GetNbinsX() << " bins along x-axis by factor " << totalRebinFactor << std::endl;
    unraveledMigrationMatrix->RebinX(totalRebinFactor);
  } else if (axis == 1) {
    std::cout << "Rebinning the unraveled histogram " << unraveledMigrationMatrix->GetName() << " with originally "
            << unraveledMigrationMatrix->GetNbinsY() << " bins along y-axis by factor " << totalRebinFactor << std::endl;
    unraveledMigrationMatrix->RebinY(totalRebinFactor);
  }
  unraveledMigrationMatrix->Reset();
  for (int i = 0; i <= unraveledMigrationMatrix->GetNbinsX() + 1; ++i) {
    for (int j = 0; j <= unraveledMigrationMatrix->GetNbinsY() + 1; ++j) {
      unraveledMigrationMatrix->SetBinContent(i, j, rebinnedMigrationMatrix->GetBinContent(i, j));
      unraveledMigrationMatrix->SetBinError(i, j, rebinnedMigrationMatrix->GetBinError(i, j));
    }
  }
  delete rebinnedMigrationMatrix;
}
