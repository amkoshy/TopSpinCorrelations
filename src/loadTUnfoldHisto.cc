#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include "../../common/include/CommandLineParameters.h"

#include <TROOT.h>
#include <TUnixSystem.h>
#include <TChain.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TString.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include "makeTUnfoldHisto.h"

using std::cerr;
using std::cout;
using std::endl;
using std::getline;
using std::ifstream;
using std::ios;
using std::string;
using std::vector;

int main(int argc, char** argv) {
  CLParameter<std::string> opt_y("y", "Era", true, 1, 1);
  CLParameter<std::string> opt_s("s", "Systematic", true, 1, 1);
  CLParameter<std::string> opt_c("c", "Channel", true, 1, 1);
  CLParameter<std::string> opt_f("f", "Filename", false, 1, 1);
  CLAnalyser::interpretGlobal(argc, argv);

  // Call the executable like this, with the filename being optional
  // ./install/bin/loadTUnfoldHisto -y <era> -s <systematics> -c <channel> -f <filename>

  
  ifstream input;
  string systematics;
  string era(opt_y[0]); // The era unfolding is running on. Current options: 2017, 2018, 2016preVFP, and 2016postVFP
  string outputsyst(opt_s[0]);
  string channel(opt_c[0]);
  TString validFilenamePattern = opt_f.isSet() ? opt_f[0] : "";

  double dyscale(1.0); 
  double bgscale(1.0); 
  double scalingval(1.0); 
  bool symmetrizeSpinVar;

  
  const int n_sub_bins = 4;       //number of sub-bins filled for each of the configured bins


  if (outputsyst.find("DY_UP") != string::npos) {
    systematics = "Nominal"; 
    dyscale=1.3;
  } else if (outputsyst.find("DY_DOWN") != string::npos) {
	  systematics = "Nominal";  //Just vary the nominal DY and BG systematics by 30%
    dyscale=0.7;
  } else if (outputsyst.find("BG_UP") != string::npos) {
    systematics = "Nominal"; 
    bgscale=1.3;
  } else if (outputsyst.find("BG_DOWN") != string::npos) {
	  systematics = "Nominal";  //Just vary the nominal DY and BG systematics by 30%
    bgscale=0.7;
  } else {
    systematics = outputsyst;
  } 

  input.open("TUnfoldFileLists_" + era  +"/TUnfoldFileList_" + systematics + "_" + channel + ".txt");
  if (!input.is_open()) {
    std::cout << "error opening file TUnfoldFileLists_" + era  +"/TUnfoldFileList_" + systematics + "_" + channel + ".txt\n";
    exit(0);
  }

  vector<string> pathplusfileVec;
  vector<string> outputVec;

  while (!input.eof()) {
    string filename;

    getline(input, filename);

    if (filename == "" || filename[0]=='#') continue;  // skip empty lines

    if(opt_f.isSet() && (filename.find(validFilenamePattern) == string::npos)) continue;

    pathplusfileVec.push_back(filename);
    string filebegins;
    if (channel != "combined") {
	    filebegins = channel + "_";
      std::size_t pos = filename.find(filebegins);
      string justfilename = filename.substr(pos);
	    outputVec.push_back("histosTUnfold_"+justfilename);
    }
  }

  input.clear();  //restart reading input list from the beginning
  input.seekg(0, ios::beg);

  for (unsigned int i = 0; i < pathplusfileVec.size(); i++) {
    string filename = pathplusfileVec[i];
    std::cout << "Full filename: " << filename.c_str() << endl;

    TFile *fin = TFile::Open(filename.c_str());
    TTree* treestep0 = (TTree*)fin->Get("ttBar_treeVariables_step0");
    TTree* treestep8 = (TTree*)fin->Get("ttBar_treeVariables_step8");
    TH1D* wtdEvts    = (TH1D*) fin->Get("weightedEvents");

    // Set to true for blinding
    // symmetrizeSpinVar = true;
    // Amandeep : Changing here for unblinded 2016 results and matrix
    symmetrizeSpinVar    = true;
    // End

    if (filename.find("dy") != string::npos) { 
      scalingval =  dyscale * bgscale;
    } 
    else if (filename.find("run") != string::npos) {
      // Amandeep : Changing here for unblinded 2016 results and matrix
      // symmetrizeSpinVar = true;
      symmetrizeSpinVar    = true;
      // End

      if (outputsyst.find("Nominal") != string::npos) { 
	    scalingval = 1;
      } 
      else { continue; }
    } 
    else if (filename.find("ttbarsignal") != string::npos) {
      scalingval = 1;
    } 
    else {
      // Scale all the backgrounds by 30% up and down for up-down distributions 
      // For nominal case scale by default value of 1
      scalingval = bgscale; 
    }

    TChain* chainstep0 = new TChain("ttBar_treeVariables_step0"); 
    TChain* chainstep8 = new TChain("ttBar_treeVariables_step8"); 
    string outputfilename; 

    if (channel == "combined") {
      if (filename.find("ee/ee_") != string::npos) {
	std::size_t pos   = filename.find("ee/ee_");
	string samplename = filename.substr(pos+6);
	outputfilename = "histosTUnfold_combined_"+samplename; 
	string path    = filename.substr(0,pos);
	//TFile fin_emu((path+"emu/emu_"+samplename).c_str(),"READ");
	string fin_emu = path+"emu/emu_"+samplename;
	
  chainstep0->Add(fin_emu.c_str(), 0);
	chainstep8->Add(fin_emu.c_str(), 0);

	//TFile fin_mumu((path+"mumu/mumu_"+samplename).c_str(),"READ");
	string fin_mumu = path+"mumu/mumu_"+samplename;
	chainstep0->Add(fin_mumu.c_str(), 0);
	chainstep8->Add(fin_mumu.c_str(), 0);
      } else {
	continue;
      } 
    }

    chainstep0->Add(filename.c_str(), 0);
    chainstep8->Add(filename.c_str(), 0);

    //makeTUnfoldHisto a(treestep0, treestep8);
    makeTUnfoldHisto a(chainstep0, chainstep8);
    if (channel != "combined") {
      a.SetOutputFileName(outputVec[i], outputsyst, channel, era);  
    } else {
      a.SetOutputFileName(outputfilename, outputsyst, channel, era);  
    }
    a.BookHisto(n_sub_bins);
    a.Loop(wtdEvts, channel, symmetrizeSpinVar);
    a.Terminate(scalingval);
  }

  return 0;
}
