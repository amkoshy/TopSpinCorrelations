#ifndef utils_h
#define utils_h

#include <vector>
#include <set>
#include <utility>
#include <fstream>
#include <string>

#include "../../common/include/classes.h"
#include "../../common/include/classesFwd.h"

#include <TMath.h>
#include <TROOT.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TColor.h>
#include <TLegend.h>
#include <TGaxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>

class TString;
class UnfoldingXSecHandler;

namespace styleUtils{

    void setHHStyle(TStyle& HHStyle);
    void setResultLegendStyle(TLegend *leg, const bool result);
}


namespace utils{

    void addAndDelete_or_Assign(TH1*& addToThis, TH1* addThis);
    void drawRatio(TH1* histNumerator, TH1* histDenominator, const TH1* uncband,const Double_t& ratioMin, const Double_t& ratioMax, TString ratioType = "cp",
                   TH1D* hist=NULL, TGraphAsymmErrors* graph = NULL ); //TH1D* hist=NULL - needed for all bins in 2d plotting

    // OZ 15.01.2017 for PlotterForTUnfold
    void drawRatio2d(TH1* histNumerator, TH1* histDenominator, const Double_t& ratioMin, const Double_t& ratioMax, TString ratioType = "cp",
                   TH1D* hist=NULL, TGraphAsymmErrors* graph = NULL ); //TH1D* hist=NULL - needed for all bins in 2d plotting
    void drawRatio2d_Eban(TH1* histNumerator, TH1* histDenominator, const Double_t& ratioMin, const Double_t& ratioMax, TString ratioType = "cp",TGraphAsymmErrors* graph1 = NULL,TGraphAsymmErrors* graph2 = NULL );

    void rebin2d(TH2D*& histOut,TH2D* histIn,TString name,TString xAxisName_,TString yAxisName_,const int rebinX_,const int rebinY_, const int nbinX_, const double* x_binsArr_, const int nbinY_, const double* y_binsArr_);
    TString numToString(double val);
    TString makeBinTitle(TString axisName,double x1,double x2);
    TString makeTitleBins(TString plotNameUnits,std::vector<double>& v_bin,int underflow, int overflow);

    ///Read line of numbers from txt file and write them to the vector
    void readLineToVector(const TString& file, const TString& keyWord,std::vector<double>& outVector);

    ///Printing out whole file
    void cat(const TString& file);


    /// Fill Up and Down systematic variations
    void fillSysUpDown(const std::vector<double>& vDiffUp,const std::vector<double>& vDiffDown,
                                     std::vector<double>& vSystErrUp, std::vector<double>& vSystErrDown);

    ///Squared sum of two vectors
    void addSqrVect(std::vector<double>& a, const std::vector<double>& b);

    ///Sum of two vectors
    std::vector<double> addVect(const std::vector<double>& a,const std::vector<double>& b,
                                const double scale=1, const int precision=-999);

    ///Difference of two vectors
    std::vector<double> diffVect(const std::vector<double>& a,const std::vector<double>& b,
                                 const double scale=1, const int precision=-999);

    ///Divide vectors
    std::vector<double> divideVect(const std::vector<double>& a,const std::vector<double>& b,
                                   const double scale=1, const int precision=-999);

    ///Multiply vectors
    std::vector<double> multiplyVect(const std::vector<double>& a,const std::vector<double>& b,
                                      const double scale=1, const int precision=-999);

    ///Relative difference of two vectors
    std::vector<double> relativeDiff(const std::vector<double>& a,const std::vector<double>& b,
                                     const double scale=1, const int precision=-999);

    ///Set precision for all vector elements
    void setPrecision(std::vector<double>& a,int precision);

    /// Prepare canvas and legend
    TCanvas* setCanvas();

    TLegend* setLegend(const double x1 = 0.53,const double y1 =0.60,const double x2 =0.90,const double y2 =0.85);

    /// Return the path where relevant input data is stored
    const std::string DATA_PATH_DILEPTONIC();

    /// Get vector of pairs with < stepName, objectName from ROOT-file (e.g. histo or tree) >
    std::vector<std::pair<TString, TString> > nameStepPairs(const char* filename,
                                                            const char* objectNameFragment,
                                                            const std::vector<TString>& selectedSteps =std::vector<TString>());



    /// Assign a step name for a given short name of the step, and potentially for a specific category
    TString stepName(const TString& stepShort, const int& category =-1);

    /// Assign a step name for a given short name of the step, and potentially for a set of categories
    TString stepName(const TString& stepShort, const std::vector<int>& v_category);



    /// Get from a TString the selection step of pattern "step*"
    TString extractSelectionStep(const TString& name);

    /// Get from a TString the jet category of pattern "cate*"
    TString extractJetCategory(const TString& name);

    /// Get from a TString the selection step of pattern "step*", combined with 0, 1 or several categories of pattern "cate*"
    TString extractSelectionStepAndJetCategory(const TString& name);

    /// Helper functions only needed by functions defined in this file
    namespace helper{

        /// Helper function to get the fragment containing a searchPattern,
        /// fragments are separated by given token,
        /// Flag allowMultiple to switch whether pattern can exist several times or only once
        TString searchFragmentByToken(const TString& name, const TString& searchPattern, const TString& token, const bool allowMultiple =false);
        uint GetNumberOfNonZeroBins(TH1& inputHisto);
    }



    /// Assign a folder depending on channel and systematic
    TString assignFolder(const char* baseDir, const TString& channel, const TString& systematic);

    /// Access an already existing input folder
    TString accessFolder(const char* baseDir, const TString& channel,
                         const TString& systematic, const bool allowNonexisting =false);

    /// Read NameList file  and fill plotNames vector
    void setPlotNames(const std::string& nameListFile, std::vector<std::vector<TString>>& varNames);

    TString GetNameFromUpDownVariation(TString name_with_up);

    void selectIndicesFromFlags(std::vector<int>&, const std::vector<bool>&);
    std::vector<bool> flagsForPassingCutInPtRange(const std::vector<int>&, const int, const std::vector<LV>&, const float, const float);

    // ------------------------------------
    /* --- JetVetoMap ----
     *
     *  Description:
     *
     *   - this class serves as interface to the ROOT file(s)
     *     provided by the JetMET POG for Jet VetoMaps
     *
     *   - the class loads the SF histogram
     *     and returns the the veto information
     *     based on the provided input (jet kinematics)
     *
     *  Class Members:
     *
     *   - th2_format  [string] = format of input TH2 object (observables used for x-axis and y-axis)
     *
     */
    class JetVetoMap
    {
    public:
      JetVetoMap(const std::string& name, const bool verbose=false);
      virtual ~JetVetoMap();

      void set_name(const std::string& name){ name_ = name; }
      void set_verbose(const bool v){ verbose_ = v; }
      float veto(const LV&, const float& cutValue) const;
      void initMap(const std::string&, const std::string&, const std::string&);

    protected:
      std::string name_;
      bool verbose_;

      class Map {

      public:
      Map() : th2_(nullptr), th2_format_("") {}
        Map(const std::string&, const std::string&, const std::string&, const bool vebose=false);
        ~Map(){}

        void getValues(float&, const LV&) const;

      protected:
        TH2* th2_;
	std::string th2_format_;
      };

      void initMap(Map&, const std::string&, const std::string&, const std::string&);

      Map map_;
    };

    std::vector<bool> flagsForPassingJetVetoMap(const JetVetoMap* vetoMap, const std::vector<LV>& v_p4, const float& cutValue);


    namespace ttmd
    {
      void HCovToFArrays(const TH1* h, const TH2* hCov, float* d, float* dCovDiag, float* dCovUndiag);

      void FArraysToHCov(const int n, const float* d, const float* dCovDiag, const float* dCovUndiag, TH1* h, TH2* hCov);

      template <class T>
      std::vector<T*> CreateVectorPtr(const int n1, const T* ptrInit = NULL)
        {
          std::vector<T*> ret(n1);
          if(ptrInit)
	    for(int j = 0; j < n1; j++)
	      ret[j] = (T*) ptrInit->Clone();
          return ret;
        }

      template <class T>
      std::vector<std::vector<T*> > CreateVectorPtr2D(const int n1, const int n2, const T* ptrInit = NULL)
        {
          std::vector<std::vector<T*> > ret(n1);
          for(int i = 0; i < n1; i++)
	    ret[i] = CreateVectorPtr<T>(n2, ptrInit);
          return ret;
        }

      template <class T>
      std::vector<std::vector<std::vector<T*> > > CreateVectorPtr3D(const int n1, const int n2, const int n3, const T* ptrInit = NULL)
        {
          std::vector<std::vector<std::vector<T*> > > ret(n1);
          for(int i = 0; i < n1; i++)
	    ret[i] = CreateVectorPtr2D<T>(n2, n3, ptrInit);
          return ret;
        }

      template <class T>
      void ClearVectorWithContent2D(std::vector<std::vector<T*> >& vv)
        {
          for(unsigned int i = 0; i < vv.size(); i++)
	    for(int j = 0; j < vv[i].size(); j++)
	      delete vv[i][j];
          vv = std::vector<std::vector<T*> >();
        }

      template <class T>
      void ClearVectorWithContent(std::vector<T*>& vv)
        {
          for(unsigned int i = 0; i < vv.size(); i++)
	    delete vv[i];
          vv = std::vector<T*>();
        }

      double DVectorAppendNWidth(std::vector<double>& vec, const double min, const int n, const double width);

      std::pair<double, double> GetMinMaxFromTH1(const TH1* h);

      std::pair<double, double> GetMinMaxFromTH1DVector(const std::vector<TH1D*>& vH);

      std::pair<double, double> GetMinMaxFromTH1DVector2D(const std::vector<std::vector<TH1D*> >& vH);

      TH2D* MakeProbMatrix(const TH2D* hResp);

      TH2D* MakeCorrMatrix(const TH2D* hCov);

      TH2D* MakeCovMatrixFromCorr(const TH2D* hCor, const TH1D* hVec);

      void DivideSquareWH(const TCanvas* fCanvas, const int n, int& w, int& h);

      void SetErrorsToZero(TH1* h);

      double Chi2(const TH1* hData, const TH2* covmat, const TH1* hGen, int skip = 0);
      double Chi2(const TH1D* hData, const TH2D* covmat, const TH1* hGen, const std::vector<int>& vSkip);

      TH1D* Convolute(const TH1* hGen, const TH2* hResp);

      TH2D* MakeCovTH2(const TH1* hData);

      TGraphAsymmErrors* DrawAsGraph(const std::pair<TH1D*, TH1D*> pairUD, const TH1D* h, const bool flagUnc = true, const double offset = 0.5, const TString& option = "ZP0");
      TGraphAsymmErrors* DrawAsGraph(const TH1D* h, const bool flagUnc = true, const double offset = 0.5, const TString& option = "ZP0");
      TGraphAsymmErrors* DrawBand(const std::pair<TH1D*, TH1D*>& hPairUD, const TH1D* h);

      void AddVarInQuadrature(double& u, double&d, const double& var);
      void AddVarInQuadrature(double& u, double&d, const std::vector<double>& var);
      void AddVarInQuadrature(std::pair<TH1D*, TH1D*> hud, const TH1D* hv);
      void AddVarInQuadrature(std::pair<TH1D*, TH1D*> hud, std::pair<TH1D*, TH1D*> hv);

      void AddVarToEnvelope(double& u, double&d, const double& var);
      void AddVarToEnvelope(double& u, double&d, const std::vector<double>& var);
      void AddVarToEnvelope(std::pair<TH1D*, TH1D*> hud, const TH1D* hv);

      void ScaleRelativeHistoVar(TH1* hvar, const double scale, const TH1* hNom, const bool flagKeepArea = false);

      double GetMeanDeviationFromRef(const TH1* h, const TH1* hRef);
      double GetMeanDeviationFromRef(const TH1* h1, const TH1* h2, const TH1* hRef);

      void AddVarToUncTH1D(std::vector<std::vector<TH1D*> >& vSysUnc, std::vector<std::vector<TH1D*> > vVarUnc);

      void Normalise(TH1* h, TH2* covmat);

      bool IsEqual(const double val1, const double val2, const double eps = 1e-6);

      extern std::vector<TString> gTCanvasExt;
      void SaveCanvas(const TVirtualPad* c, const TString& baseName);

      void ScaleAxisFonts(TAxis* axis, const double scale);

      void ScaleHistoFonts(TH1* h, const double scale);
      // TODO name is confusing
      void ScaleHistoFonts(TGraph* h, const double scale);

      void SetAxisFonts(TAxis* axis, const int font = 63, const int size = 20);

      void SetHistoAxisFonts(TH1* h, const int font = 63);
      // TODO name is confusing
      void SetHistoAxisFonts(TGraph* h, const int font = 63);

      void AdjustDividedCanvasHistograms(TH2D* hr, const int p, const int npads, const int w/*, const int h*/, const double scaleOffset = 1.0);

      void SetOZStyle();

      void SetBWPalette(int ncol = 99);

      void SetNicePalette(int ncol = 99);

      void SetNicePaletteEBW();

      //int GetNiceColor(int i, int n);

      int GetDistinctColor(const int i);
      int GetDistinctColor(const TString& pdf);

      int GetDistinctMarkerStyle(int i);

      int GetDistinctLineStyle(int i);

      // 0 file exists, 1 no file but could be created, 2 dir was created
      int CheckFile(const TString& fileName, const bool flagCreateDir = true);

      bool FileExists(const TString& name, const bool flagThrow = false);

      //void ClearDir(const TString& dirName);

      template <class T1, class T2>
	std::pair<T1*, T2*> CopyPairPtr(const std::pair<T1*, T2*> in)
        {
          std::pair<T1*, T2*> out;
          out.first = new T1(*in.first);
          out.second = new T1(*in.second);
          return out;
        }

      template<class T>
	void AddBranchVerbosely(TTree* tree, const char* name, T* ptr)
        {
          printf("Adding branch %20s ... ", name);
          auto ret = tree->Branch(name, ptr);
          if(ret)
	    printf("success\n");
          else
	    {
	      printf("failed\n");
	      throw std::logic_error(TString::Format("Error while adding branch %s: returned %d\n", name, ret).Data());
	    }
        }

      template<class T>
	void SetBranchAddressVerbosely(TTree* tree, const char* name, T& ptr, TBranch** br)
        {
          //printf("Accessing branch %20s ... ", name);
          tree->SetBranchStatus(name, 1);
          auto ret = tree->SetBranchAddress(name, &ptr, br);
          if(ret == 0)
	    ;//printf("success\n");
          else
	    {
	      printf("Accessing branch %20s ... ", name);printf("failed\n");
	      throw std::logic_error(TString::Format("Error while accessing branch %s: returned %d\n", name, ret).Data());
	    }
        }

      double GetBinContentForValueNoUF(TH2* h, const double x, const double y, int* counterPtr = NULL);

      void ErFgets();

      int ReadNumberOfLines(const TString& fileName);
      int ReadNumberOfColumns(const TString& fileName, const int lineStart = 1);
      void ReadVectorKFactor(std::vector<double>& v, const TString& fileName, const int column, const int lineStart, const bool flagMakeNorm = false);
      void ReadHistoKFactor(TH1D* h, const TString& fileName, const int column, const int lineStart, const bool flagMakeNorm = false);

      double Rhoj(const double mttj, const int njet = 0, const double ptjMin = 0.0);

      TString GetNameExtraJetVar(const TString& varName, const TString& def, const double iso, const double pt, const TString& leadStr = "");

      TString Execute(const TString& cmd);

      TH1D* CloneTH1D(const TH1D* h, const double content = 0.0);

      std::vector<double> MakeSubbins(const std::vector<double>& v, const int nsubbins);

      const TString TexAdoptVarTitle(const TString& str);
      const TString TexAdoptBinLabel(const TString& str);
      const TString GetDigitTimesFormat(const double val, const int nDig);
      const TString TexAdoptVarYAxisTitle(const TString& str);
      const TString TexXsecTitle(const TString& str, TString& strExtraBins);
      const TString TexXsecLabel(const TString& str, const int nj);

      void DoRjHisto(TH1* h, TH2* covmat = NULL);
    }

    std::vector<double> GetVectorFromTStringData(TString vector_tstring);

}

#endif
