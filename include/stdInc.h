#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>
#include <stdlib.h> 

//root Header
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TString.h>
#include <TObject.h>
#include <TCollectionProxyFactory.h>
#include <TStreamerInfo.h>

// inveritng matrices
#include "Riostream.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
#include "TMatrixDLazy.h"

// from root
#include <TObjArray.h>
#include <TDatabasePDG.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h> 
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TGraphAsymmErrors.h"
#include <TLatex.h>
#include <TMath.h> 
#include <TPostScript.h>
#include <TLine.h> 

#include "TAttLine.h"
#include "TAttFill.h"
#include "TAttMarker.h"
#include "TAttText.h"
#include "TClass.h"

#include "TMinuit.h" 
#include "TMatrixD.h"
#include "TVectorF.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TVirtualFitter.h"
#include "TEllipse.h"
#include "TPave.h"
