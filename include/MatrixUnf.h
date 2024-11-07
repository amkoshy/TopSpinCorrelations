/********************************************************************************
File:         MatrixUnf.h

Author:       Andreas Jung 
 
Description:  Class for calculation the unfolding stuff. 
		This class will calculate that for a given variable. 
              

History:      	10.08.09:  first appearance of this file.
		28.04.11:  Make it independent of an environment
*********************************************************************************/
#ifndef MatrixUnf_h
#define MatrixUnf_h 

#include "stdInc.h"
#include "TUnfoldDensity.h"

#define maxUnfBins 35000

//class only to store the fit results (just a completely stupid helper class):  
class  MatrixUnf :  public TNamed { 
  
  Int_t fNBinsX;            //number of X-bins in histos (generated);
  Int_t fNBinsY;            //number of Y-bins in histos (reconstruced);
  Int_t dim;                // dim-dimensional unfolding is taking place

  TUnfoldBinning *generatorBinning;
  TUnfoldBinning *detectorBinning;

  TUnfoldBinning *generatorBinning_withBFF;
  TUnfoldBinning *detectorBinning_rebinnedB;
  TUnfoldBinning *generatorBinning_rebinnedB;

  TH1D *GenVecHist_forBias;
  TH1D *GenVecHist_forBias_rebinnedB;

  Double_t fBinArray[maxUnfBins];    
  Double_t fBin2DArray[maxUnfBins];  
  Double_t fBin2DFitArray[maxUnfBins];  
  Double_t fBin2DGenArray[maxUnfBins];  
  Double_t fBin2DRecArray[maxUnfBins];  
  Double_t fBR; 
 
  TString fVariableName1;
  TString fVariableName2; 
  TString fXAxisTitle; 
  TString fYAxisTitle; 
  TString fZAxisTitle; 
  
  TObjArray *fBkgHistArray;

  //The needed histos for the calculation: 
  TObjArray *fHistUnfArray; 
  TObjArray *fResultUnfHists; 
  TObjArray *fResultUnf2DHists; 
  TObjArray *fResSysHist;
  TObjArray *fAequiResSysHist;
  TObjArray *fTUnfDeltaSys;
  TObjArray *fTUnfDeltaSys_rebinnedA;
  TObjArray *fTUnfDeltaSys_rebinnedB;
  TObjArray *fTUnfCovMatSys;
  TObjArray *fTUnfCovMatSys_rebinnedA;
  TObjArray *fTUnfCovMatSys_rebinnedB;
  TObjArray *fTUnfCovMatSysNorm;
  TObjArray *fTUnfCovMatSysNorm_rebinnedA;
  TObjArray *fTUnfCovMatSysNorm_rebinnedB;

  // some bools
  Bool_t fsetBiasDistUnfold;
  Bool_t fsetSyst;
  Bool_t fEnsTest;
  Bool_t freweight_bias_test = kFALSE;

  Int_t frebinfine;

  //helper functions:

public: 
  //constructors and so: 
  MatrixUnf();
  MatrixUnf(TString VariableName1, TString xAxisTitle, TString VariableName2, TString yAxisTitle,  
			 TString zAxisTitle, TString binFileName, Bool_t InputIs2D, Bool_t setBiasDistUnfold, TString channel, Bool_t setSyst, TH1D* GenVec, Int_t rebinfine);
  void AddObjectsToFile();
  ~MatrixUnf();

  //getter 
  TH1D *GetGenVecHist(){return  ((TH1D*)fHistUnfArray->At(0)); };  //
  TH2D *GetRecResHist(){return  ((TH2D*)fHistUnfArray->At(1)); };  //
  TH1D *GetVisGenVecHist(){return  ((TH1D*)fHistUnfArray->At(2)); };  //
  TH1D *GetFullRecVec(){return  ((TH1D*)fHistUnfArray->At(3)); }; //rec from the input
  TH1D *GetFullRecBgVec(){return  ((TH1D*)fHistUnfArray->At(12)); }; 
  TH1D *GetMissRecVec(){return  ((TH1D*)fHistUnfArray->At(14)); };

  Int_t GetDim(){return dim;};
  TUnfoldBinning *GetGeneratorBinning(){ return generatorBinning;};
  TUnfoldBinning *GetDetectorBinning(){ return detectorBinning;};
  TUnfoldBinning *GetGeneratorBinningWithBFF(){ return generatorBinning_withBFF;};
  TUnfoldBinning *GetRebinnedDetectorBinning(){ return detectorBinning_rebinnedB;};
  TUnfoldBinning *GetRebinnedGeneratorBinning(){ return generatorBinning_rebinnedB;};

  TObjArray *GetBgHists(){return fBkgHistArray;}

  TObjArray *GetInputHists(){ return fHistUnfArray;}
  TObjArray *GetResultHists(){ return fResultUnfHists;}
  TObjArray *GetResSysHists(){ return fResSysHist;}
  TObjArray *GetAequiResSysHists(){ return fAequiResSysHist;}
  TObjArray *GetTUnfDeltaSys(){ return fTUnfDeltaSys;}
  TObjArray *GetTUnfCovMatSys(){ return fTUnfCovMatSys;}

  //results: 
  TH1D *GetAequiGen(){return ((TH1D*)fResultUnfHists->At(0)); };
  TH1D *GetEffHist(){return ((TH1D*)fResultUnfHists->At(1)); }; 
  TH2D *GetResMat(){return ((TH2D*)fResultUnfHists->At(2)); }; 

  TH1D *GetAequiMeas(){return ((TH1D*)fResultUnfHists->At(3)); }; //aequidistant measurement vector
  TH2D *GetAequiResMat(){return ((TH2D*)fResultUnfHists->At(4)); }; 
  TH1D *GetAccepHist(){return ((TH1D*)fResultUnfHists->At(5)); }; 

  TGraph *GetLogTauX(){return  ((TGraph*)fResultUnfHists->At(7)); }; 
  TGraph *GetLogTauXUser(){return  ((TGraph*)fResultUnfHists->At(8)); }; 
  TGraph *GetLCurve(){return  ((TGraph*)fResultUnfHists->At(9)); }; 
  TGraph *GetLCurveUser(){return  ((TGraph*)fResultUnfHists->At(10)); }; 
  TH1D *GetInputFoldBack(){return  ((TH1D*)fResultUnfHists->At(11)); }; 
  TH2D *GetRhoIJ(){return  ((TH2D*)fResultUnfHists->At(12)); }; 
  TH2D *GetErrMatrixData(){return  ((TH2D*)fResultUnfHists->At(13)); }; 
  TH1D *GetTUnfResult(){return  ((TH1D*)fResultUnfHists->At(14)); }; 
  TH1D *GetTUnfResultCorXsec(){return  ((TH1D*)fResultUnfHists->At(15)); }; 
  TH1D *GetTotEffHist(){return  ((TH1D*)fResultUnfHists->At(19)); }; 

  //TSpline *GetRhoIJScanTauHist(){return  ((TSpline*)fResultUnfHists->At(22)); }; 
  TGraph *GetRhoIJScanTauHist(){return  ((TGraph*)fResultUnfHists->At(22)); }; 
  TGraph *GetRhoIJOptimalHist(){return  ((TGraph*)fResultUnfHists->At(23)); }; 
  TH1D *GetTheoryXsec(){return ((TH1D*)fResultUnfHists->At(24)); };

  TH1D *GetTruthMigrationPurity(){return ((TH1D*)fResultUnfHists->At(28)); };
  TH1D *GetStability(){return ((TH1D*)fResultUnfHists->At(29)); };
  TH1D *GetAequiMeasBgSubtracted(){return ((TH1D*)fResultUnfHists->At(30)); };
  TH1D *GetResolutionOffset(){return ((TH1D*)fResultUnfHists->At(31)); };  
  TH1D *GetResolutionMean(){return ((TH1D*)fResultUnfHists->At(32)); };  
  TH1D *GetResolutionStdev(){return ((TH1D*)fResultUnfHists->At(33)); };  

  TH2D *GetEmatrix_rebinnedA(){return  ((TH2D*)fResultUnfHists->At(34)); }; 
  TH1D *GetTUnfResult_rebinnedA(){return  ((TH1D*)fResultUnfHists->At(35)); }; 
  TH1D *GetAfbResult_rebinnedA(){return  ((TH1D*)fResultUnfHists->At(36)); }; 
  TH2D *GetEmatrix_rebinnedB(){return  ((TH2D*)fResultUnfHists->At(37)); }; 
  TH1D *GetTUnfResult_rebinnedB(){return  ((TH1D*)fResultUnfHists->At(38)); }; 
  TH1D *GetAfbResult_rebinnedB(){return  ((TH1D*)fResultUnfHists->At(39)); }; 

  void SetName(TString NameTUnfObjX){ fVariableName1=NameTUnfObjX; };
  Bool_t SetInputDists(TH1D* GenVec, TH1D* VisGenVec, TH1D* RecVec, std::vector<TH1D*> BgVec, TH1D* MissRecVec, TH2D* ResMat, std::vector<TH2D*>vecResMatSys);

  Bool_t prepRunUnfolding(TH1D* measInput, Int_t RegMode, Int_t skipOuterBins, Double_t tauMin, Double_t tauMax, Double_t Lumi, Bool_t subtrBg, Int_t  useAreaConstrain, Bool_t ensTest, Bool_t minrhoavg, std::vector<TH1D*>&pseudoUnfResults, std::vector<TH1D*>&pseudoUnfResultsRebinnedA, std::vector<TH1D*>&pseudoUnfResultsRebinnedB, std::vector<TH1D*>&pseudoUnfAfbResultsRebinnedA, std::vector<TH1D*>&pseudoUnfAfbResultsRebinnedB, std::vector<TH1D*>&pseudoUnfOtherResultsRebinnedA, std::vector<TH1D*>&pseudoUnfOtherResultsRebinnedB);

  //calculate the Results: 
  Bool_t createInputResults(Bool_t ensTest); //create and fill result histos 
  Bool_t createResMat();
  Bool_t createResMatSys();
  Bool_t createEffBbBVec();
  Bool_t createAequiResMat(); 
  Bool_t createAequiResMatSys();
  Bool_t createAequiGenVec();
  Bool_t createAequiMeasVec(TH1D* MeasVec);
  Bool_t createAequiMeasVecBgSubtracted(TH1D* MeasVec);
  TH1D*  calcDiffHisto(TH1D* in1, TH1D* in2);
  Bool_t createTruthMigrationPurityStability();
  Bool_t createResolution();

  void calculateBFF(TH1D* GenVecHist_forBias, TVectorD *BFFvec, TUnfoldBinning *binning);
  void ReweightBiasTest(TH1D* GenVecHist_forBias);

  void Chi2tests(TUnfoldDensity &unfold, TUnfoldDensity &unfold_rebinnedB, TH1D* TUnfResult, TH1D* TUnfResult_rebinnedA, TH1D* TUnfResult_rebinnedB);
  void Chi2testsForPEs(TUnfoldDensity &unfold, TUnfoldDensity &unfold_rebinnedB, TH1D* TUnfResult_rebinnedA, std::vector<double> &chi2unf);
 void rebinMultidimMigrationMatrix(TH2* unraveledMigrationMatrix, const TUnfoldBinning *generatorRebinnedBinning, const TUnfoldBinning *detectorRebinnedBinning, const TUnfoldBinning *generatorBinning, const TUnfoldBinning *detectorBinning, int rebinFactor);
  void rebinMultidimMigrationMatrix(TH2* unraveledMigrationMatrix, const TUnfoldBinning* rebinnedBinning, const TUnfoldBinning* binning, int rebinFactor, int axis);

  void DumpObjArr(); 

  static void CreateBinMap(const TUnfoldBinning *binning, int n_sub_bins, int binMap[]);
  static void rebinMultidimensionalInput(TH1* unraveledHistogram, const TUnfoldBinning *rebinnedBinning, const TUnfoldBinning *originalBinning, int rebinFactor);
  static void GetAfb(TH1D* h, Double_t &afb, Double_t  &afberr);
  static void GetAfbWithCorrelations(TH1* histogram, TMatrixD &covarianceM, Double_t *afb, Double_t *afberr, int numAsymmetries);
  static void GetAfbBinByBin(TH1D* histogram, std::vector<double> &myafb, std::vector<double> &myerr, int numCoefficients);
  static void GetAfbWithCorrelationsBinByBin(TH1* histogram, TMatrixD &covarianceM, std::vector<double> &myafb, std::vector<double> &myerr, TMatrixD &AFBcovarianceM, int numAsymmetries);
  static void GetAfbCorrMAllVars(TH1D* histAllVars, TMatrixD &corrM, std::vector<double> &myafb, std::vector<double> &myerr, TMatrixD &AFBcorrM, bool isCovInput = false);
  static void GetAfbBinsCorrMAllVars(TH1D* histAllVars, TMatrixD &corrM, std::vector<std::vector<double>> &afbvecvec, std::vector<std::vector<double>> &errvecvec, TMatrixD &AFBBinscorrM, TH1D* hist_acombfactor[], TMatrixD &CcorrM, bool isCovInput = false, int numCoefficients = 1);
  static void GetBinSumCorrMAllVars(TH1D* histAllVars, TMatrixD &corrM, std::vector<double> &binsum, std::vector<double> &binsumerr, TMatrixD &BinSumcorrM, bool isCovInput = false);
  static void GetBinSumCorrMSnake(TH1D* histAllVars, TMatrixD &corrM, std::vector<double> &binsum, std::vector<double> &binsumerr, TMatrixD &BinSumcorrM, bool isCovInput = false);
  static Bool_t GetOptimalCoefficient(std::vector<double> &Afbvec, TMatrixD &m_AFB, TH1D* genhistogram, TString VariableName1, Double_t *Coef, Double_t *CoefErrOpt, Double_t *CoefErrTest, std::vector<double> &Coefvec, TMatrixD &m_C, std::vector<double> &a_optimal_reco, int numCoefficients);
  static void fillAsymsAndCoefficient(TH1D* TUnfResult, TH2D* Ematrix, TH1D* GenVecHist, TString VariableName1, TH1D* AfbResult, TH1D* OtherResult, bool ensTest, bool isgenlevel = false, int numCoefficients = 1);
  static void GetAbsoluteCrossSec(TH1D*& uncorrHist, Double_t BR, Double_t Lumi);
  static void GetBinWidthCorrectedCrossSec(TH1D*& uncorrHist);
  // Amandeep : Overloaded method for upto 3D bin width correction
  static void GetBinWidthCorrectedCrossSec(TH1D*& uncorrHist, TString VariableName);
  // End
  static void GetAbsoluteCovarianceMatrix(TH1D* unfHist, TH2D*& covMat, Double_t BR, Double_t Lumi);
  static TVectorD TH1toTVectorD(TH1* deltahist);
  static TMatrixD TH1toTMatrixD(TH1* deltahist);
  static TMatrixD TH2toTMatrixD(TH2* covMhist);
  static TH2D* TMatrixDtoTH2(TMatrixD &covM, TString name);
  //static void GetNormalisedCovarianceMatrix(TH1D* histogram, TMatrixD &covarianceM, TMatrixD &NormcovarianceM);
  static void GetNormalisedCovarianceMatrix(TH1D* histogram, TH2D* covMhist, TH2D*& NormcovMhist);
  static void GetNormalisedCovarianceMatrixAllVars(TH1D* histogram, TH2D* covMhist, TH2D*& NormcovMhist);
  static Double_t AtoCfactor(TString VariableName1);
  static Bool_t AtoCfactor(TString VariableName1, TH1D* genhistogram, std::vector<double> &factor_in, int numCoefficients, int ithCoefficient);
  static bool isCtype(TString VariableName1);
  static bool isCPMtype(TString VariableName1);
  static bool isBPMtype(TString VariableName1);
  static bool isBtype(TString VariableName1);
  static bool isDtype(TString VariableName1);
  static bool isOtherSymType(TString VariableName1);
  static bool isLinearCombtype(TString VariableName1);
  static TString CoefName(TString VariableName1, bool islatex, double lowerBoundary = -999.0, double upperBoundary = -999.0);
  static TString ObsName(TString VariableName1, bool islatex);

//  ClassDef (MatrixUnf, 0)
};


#endif
