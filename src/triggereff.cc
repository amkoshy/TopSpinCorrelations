#include "../include/triggereff.h"
#include <string>
#include "../../common/include/CommandLineParameters.h"

////
// Code to calculate trigger efficiencies. To run in root: 
// Data ex: root -l -q triggereff.cc++'("2017_UL", "Data", "All")'
// MC ex: root -l -q triggereff.cc++'("2017_UL", "MC", "ttbarsignalplustau_fromDilepton")'

////

/// Main Function



void triggereff(std::string era_, std::string MCorData_, std::string sample_, std::string sample_number_){


  system("mkdir -p triggereff");

  if (sample_ == "ttbarsignalplustau_fromDilepton") {
    suffix = "_"+MCorData_+"_"+sample_+"_"+sample_number_;
  }
  else {
    suffix = "_"+MCorData_+"_"+sample_;
  }


  TChain *chain = new TChain("writeNTuple/NTuple");

  if(MCorData_ == "MC") {

    isMC = true;

    if (era_ == "2016") {

      if (sample_ == "ttbarsignalplustau_fromDilepton") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/ttbarsignalplustau_fromDilepton_2016_Reza.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/ttbarsignalplustau_fromDilepton_2016_Reza.root"; 
      }
      else if (sample_ == "ttbarbg") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/ttbarbg.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/ttbarbg.root"; 
      }
      else if (sample_ == "dy1050") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/dy1050.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/dy1050.root"; 
      }
      else if (sample_ == "dy50inf") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/dy50inf.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/dy50inf.root"; 
      }
      else if (sample_ == "singleantitop_tw") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/singleantitop_tw.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/singleantitop_tw.root"; 
      }
      else if (sample_ == "singletop_tw") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/singletop_tw.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/singletop_tw.root"; 
      }
      else {
	cout << "Third argument should be ttbarsignalplustau_fromDilepton, ttbarbg, dy1050, dy50inf, singleantitop_tw, or singletop_tw." << endl;
	std::exit(1);
      }
    }

    else if (era_ == "2017") {

      if (sample_ == "ttbarsignalplustau_fromDilepton") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2017/TriggerEff_Sep2020/ttbarsignalplustau_fromDilepton_2017_Reza.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/2017/TriggerEff_Sep2020/ttbarsignalplustau_fromDilepton_2017_Reza.root"; 
      }
      else if (sample_ == "ttbarbg_fromLjets") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2017/TriggerEff_Sep2020/ttbarbg_fromLjets.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/2017/TriggerEff_Sep2020/ttbarbg_fromLjets.root"; 
      }
      else if (sample_ == "ttbarbg_fromHadronic") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2017/TriggerEff_Sep2020/ttbarbg_fromHadronic.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/2017/TriggerEff_Sep2020/ttbarbg_fromHadronic.root"; 
      }
      else if (sample_ == "dy1050") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2017/TriggerEff_Sep2020/dy1050.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/2017/TriggerEff_Sep2020/dy1050.root"; 
      }
      else if (sample_ == "dy50inf") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2017/TriggerEff_Sep2020/dy50inf_amcatnlofxfx.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/2017/TriggerEff_Sep2020/dy50inf_amcatnlofxfx.root"; 
      }
      else if (sample_ == "singleantitop_tw") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2017/TriggerEff_Sep2020/singleantitop_tw.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/2017/TriggerEff_Sep2020/singleantitop_tw.root"; 
      }
      else if (sample_ == "singletop_tw") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2017/TriggerEff_Sep2020/singletop_tw.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/2017/TriggerEff_Sep2020/singletop_tw.root"; 
      }
    }

    else if (era_ == "2018") {

      if (sample_ == "ttbarsignalplustau_fromDilepton") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/ttbarsignalplustau_fromDilepton.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/ttbarsignalplustau_fromDilepton.root"; 
      }
      else if (sample_ == "ttbarbg_fromLjets") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/ttbarbg_fromLjets.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/ttbarbg_fromLjets.root"; 
      }
      else if (sample_ == "ttbarbg_fromHadronic") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/ttbarbg_fromHadronic.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/ttbarbg_fromHadronic.root"; 
      }
      else if (sample_ == "dy1050") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/dy1050.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/dy1050.root"; 
      }
      else if (sample_ == "dy50inf") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/dy50inf.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/dy50inf.root"; 
      }
      else if (sample_ == "singleantitop_tw") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/singleantitop_tw.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/singleantitop_tw.root"; 
      }
      else if (sample_ == "singletop_tw") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/singletop_tw.root");
	filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/singletop_tw.root"; 
      }
      else {
	cout << "Third argument should be ttbarsignalplustau_fromDilepton, ttbarbg, dy1050, dy50inf, singleantitop_tw, or singletop_tw." << endl;
	std::exit(1);
      }
    }

    else if (era_ == "2016preVFP_UL") {

      if (sample_ == "ttbarsignalplustau_fromDilepton") {
	chain->Add(("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016preVFP/fullRun2_UL_November2022/ttbarsignalplustau_fromDilepton_2016ULpreVFP_"+sample_number_+".root").c_str());
	filename = ("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016preVFP/fullRun2_UL_November2022/ttbarsignalplustau_fromDilepton_2016ULpreVFP_"+sample_number_+".root").c_str(); 
      }
    }

    else if (era_ == "2016postVFP_UL") {

      if (sample_ == "ttbarsignalplustau_fromDilepton") {
	chain->Add(("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016postVFP/fullRun2_UL_November2022/ttbarsignalplustau_fromDilepton_2016ULpostVFP_"+sample_number_+".root").c_str());
	filename = ("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016postVFP/fullRun2_UL_November2022/ttbarsignalplustau_fromDilepton_2016ULpostVFP_"+sample_number_+".root").c_str(); 
      }
    }

    else if (era_ == "2017_UL") {

      if (sample_ == "ttbarsignalplustau_fromDilepton") {
	chain->Add(("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2017/fullRun2_UL_November2022/ttbarsignalplustau_fromDilepton_2017UL_"+sample_number_+".root").c_str());
	filename = ("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2017/fullRun2_UL_November2022/ttbarsignalplustau_fromDilepton_2017UL_"+sample_number_+".root").c_str(); 
      }
    }

    else if (era_ == "2018_UL") {

      if (sample_ == "ttbarsignalplustau_fromDilepton") {
	chain->Add(("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2018/fullRun2_UL_November2022/ttbarsignalplustau_fromDilepton_2018UL_"+sample_number_+".root").c_str());
	filename = ("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2018/fullRun2_UL_November2022/ttbarsignalplustau_fromDilepton_2018UL_"+sample_number_+".root").c_str(); 
      }
    }


  }

  else if(MCorData_ == "Data") {

    isMC = false;

    if (era_ == "2016") {

      filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/met_run2016B.root"; 

      if (sample_ == "All") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/met_run2016B.root");
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/met_run2016C.root");
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/met_run2016D.root");
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/met_run2016E.root");
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/met_run2016F1.root");   
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/met_run2016F2.root");   
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/met_run2016G.root");   
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/met_run2016H.root");         
      }
      else if (sample_ == "run2016B") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/met_run2016B.root");
      }
      else if (sample_ == "run2016C") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/met_run2016C.root");
      }
      else if (sample_ == "run2016D") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/met_run2016D.root");
      }
      else if (sample_ == "run2016E") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/met_run2016E.root");
      }
      else if (sample_ == "run2016F") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/met_run2016F1.root");
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/met_run2016F2.root");
      }
      else if (sample_ == "run2016G") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/met_run2016G.root");
      }
      else if (sample_ == "run2016H") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/2016/TriggerEff_Sep2020/met_run2016H.root");
      }
      else {
	cout << "Third argument should be All or run201#*." << endl;
	std::exit(1);
      }
    }

    else if (era_ == "2017") {

      filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2017/TriggerEff_preUL_Dec2021/met_run2017B.root"; 

      if (sample_ == "All") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2017/TriggerEff_preUL_Dec2021/met_run2017B.root");
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2017/TriggerEff_preUL_Dec2021/met_run2017C.root");
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2017/TriggerEff_preUL_Dec2021/met_run2017D.root");
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2017/TriggerEff_preUL_Dec2021/met_run2017E.root");
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2017/TriggerEff_preUL_Dec2021/met_run2017F.root");
      } 
      else if (sample_ == "run2017B") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2017/TriggerEff_preUL_Dec2021/met_run2017B.root");
      }
      else if (sample_ == "run2017C") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2017/TriggerEff_preUL_Dec2021/met_run2017C.root");
      }
      else if (sample_ == "run2017D") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2017/TriggerEff_preUL_Dec2021/met_run2017D.root");
      }
      else if (sample_ == "run2017E") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2017/TriggerEff_preUL_Dec2021/met_run2017E.root");
      }
      else if (sample_ == "run2017F") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2017/TriggerEff_preUL_Dec2021/met_run2017F.root");
      }
      else {
	cout << "Third argument should be All or run201#*." << endl;
	std::exit(1);
      }
    }

    else if (era_ == "2018") {

      filename = "/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/met_run2018A.root"; 

      if (sample_ == "All") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/met_run2018A.root");
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/met_run2018B.root");
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/met_run2018C.root");
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/met_run2018D.root");
      } 
      else if (sample_ == "run2018A") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/met_run2018A.root");
      }
      else if (sample_ == "run2018B") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/met_run2018B.root");
      }
      else if (sample_ == "run2018C") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/met_run2018C.root");
      }
      else if (sample_ == "run2018D") {
	chain->Add("/mnt/hadoop/store/group/local/cmstop/jthiema/ntuples2018/TriggerEff_Sep2020/met_run2018D.root");
      }
      else {
	cout << "Third argument should be All or run201#*." << endl;
	std::exit(1);
      }
    }
 
    else if (era_ == "2016preVFP_UL") {

      filename = "root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016preVFP/fullRun2_UL_April2022/met_run2016B.root"; 

      if (sample_ == "All") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016preVFP/fullRun2_UL_April2022/met_run2016B.root");
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016preVFP/fullRun2_UL_April2022/met_run2016C.root");
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016preVFP/fullRun2_UL_April2022/met_run2016D.root");
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016preVFP/fullRun2_UL_April2022/met_run2016E.root");
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016preVFP/fullRun2_UL_April2022/met_run2016F1.root");
      } 
      else if (sample_ == "run2016_ULB") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016preVFP/fullRun2_UL_April2022/met_run2016B.root");
      }
      else if (sample_ == "run2016_ULC") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016preVFP/fullRun2_UL_April2022/met_run2016C.root");
      }
      else if (sample_ == "run2016_ULD") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016preVFP/fullRun2_UL_April2022/met_run2016D.root");
      }
      else if (sample_ == "run2016_ULE") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016preVFP/fullRun2_UL_April2022/met_run2016E.root");
      }
      else if (sample_ == "run2016_ULF1") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016preVFP/fullRun2_UL_April2022/met_run2016F1.root");
      }
      else {
	cout << "Third argument should be All or run201#*." << endl;
	std::exit(1);
      }
    }

    else if (era_ == "2016postVFP_UL") {

      filename = "root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016postVFP/fullRun2_UL_April2022/met_run2016F2.root"; 

      if (sample_ == "All") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016postVFP/fullRun2_UL_April2022/met_run2016F2.root");
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016postVFP/fullRun2_UL_April2022/met_run2016G.root");
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016postVFP/fullRun2_UL_April2022/met_run2016H.root");
      } 
      else if (sample_ == "run2016_ULF2") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016postVFP/fullRun2_UL_April2022/met_run2016F2.root");
      }
      else if (sample_ == "run2016_ULG") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016postVFP/fullRun2_UL_April2022/met_run2016G.root");
      }
      else if (sample_ == "run2016_ULH") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2016postVFP/fullRun2_UL_April2022/met_run2016H.root");
      }
      else {
	cout << "Third argument should be All or run201#*." << endl;
	std::exit(1);
      }
    }

    else if (era_ == "2017_UL") {

      filename = "root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2017/fullRun2_UL_April2022/met_run2017B.root"; 

      if (sample_ == "All") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2017/fullRun2_UL_April2022/met_run2017B.root");
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2017/fullRun2_UL_April2022/met_run2017C.root");
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2017/fullRun2_UL_April2022/met_run2017D.root");
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2017/fullRun2_UL_April2022/met_run2017E.root");
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2017/fullRun2_UL_April2022/met_run2017F.root");
      } 
      else if (sample_ == "run2017_ULB") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2017/fullRun2_UL_April2022/met_run2017B.root");
      }
      else if (sample_ == "run2017_ULC") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2017/fullRun2_UL_April2022/met_run2017C.root");
      }
      else if (sample_ == "run2017_ULD") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2017/fullRun2_UL_April2022/met_run2017D.root");
      }
      else if (sample_ == "run2017_ULE") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2017/fullRun2_UL_April2022/met_run2017E.root");
      }
      else if (sample_ == "run2017_ULF") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2017/fullRun2_UL_April2022/met_run2017F.root");
      }
      else {
	cout << "Third argument should be All or run201#*." << endl;
	std::exit(1);
      }
    }

    else if (era_ == "2018_UL") {

      filename = "root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2018/fullRun2_UL_April2022/met_run2018A.root"; 

      if (sample_ == "All") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2018/fullRun2_UL_April2022/met_run2018A.root");
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2018/fullRun2_UL_April2022/met_run2018B.root");
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2018/fullRun2_UL_April2022/met_run2018C.root");
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2018/fullRun2_UL_April2022/met_run2018D.root");
      } 
      else if (sample_ == "run2018_ULA") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2018/fullRun2_UL_April2022/met_run2018A.root");
      }
      else if (sample_ == "run2018_ULB") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2018/fullRun2_UL_April2022/met_run2018B.root");
      }
      else if (sample_ == "run2018_ULC") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2018/fullRun2_UL_April2022/met_run2018C.root");
      }
      else if (sample_ == "run2018_ULD") {
	chain->Add("root://xrootd.rcac.purdue.edu//store/group/local/cmstop/jthiema/ntuples2018/fullRun2_UL_April2022/met_run2018D.root");
      }
      else {
	cout << "Third argument should be All or run201#*." << endl;
	std::exit(1);
      }
    }

  }

  else {
    cout << "Second argument should be MC or Data." << endl;
    std::exit(1);
  }


/// To get HLTPaths
  TFile* file = TFile::Open(filename);
  vector<string>* p_HLTPaths = (vector<string>*) file->Get("writeNTuple/HLTPaths");
  SetHLTBitsManager(BitsManager<string>(*p_HLTPaths));
  cout << "print HLTBitsManager_:" << endl;
  HLTBitsManager_.print();
///

  TH1D *weightedEvents = (TH1D*) file->Get("EventsBeforeSelection/weightedEvents");
  TH1D *unweightedEvents = (TH1D*) file->Get("EventsBeforeSelection/unweightedEvents");


  chain->SetBranchAddress("lepPdgId",& lepPdgId, & b_lepPdgId);
  chain->SetBranchAddress("leptons", & v_leptons, & b_leptons);
  chain->SetBranchAddress("jets", & v_jets, & b_jets);
  chain->SetBranchAddress("NPV_all" , & NPV_all, & b_NPV_all ); 
  chain->SetBranchAddress("NPV_good", & NPV_good, & b_NPV_good); 
  chain->SetBranchAddress("met_T1", & met, & b_met); 
  chain->SetBranchAddress("lepID_MuonTight"   , & lepID_MuonTight   , & b_lepID_MuonTight);
  chain->SetBranchAddress("lepID_ElecCutBased_original", &lepID_ElecCutBased, & b_lepID_ElecCutBased); 
  chain->SetBranchAddress("HLTBits", & HLTBits, & b_HLTBits);
  chain->SetBranchAddress("metFilters", &metFilters, &b_metFilters);
  chain->SetBranchAddress("lepPfIso", &lepPfIso, &b_lepPfIso);
  chain->SetBranchAddress("lepSCEta", &lepSCEta, & b_lepSCEta);
  chain->SetBranchAddress("lepID_ElecMVAIso"  , & lepID_ElecMVAIso, & b_lepID_ElecMVAIso);
  chain->SetBranchAddress("lepID_MVATTH"  , & lepID_MVATTH , & b_lepID_MVATTH); 
  chain->SetBranchAddress("jetPFID"   , & jetPFID   , & b_jetPFID);
  chain->SetBranchAddress("jetBTagDeepCSV"   , & jetBTagDeepCSV   , & b_jetBTagDeepCSV);
  chain->SetBranchAddress("jetBTagDeepJet"   , & jetBTagDeepJet   , & b_jetBTagDeepJet);
  chain->SetBranchAddress("lepEnergyCorr"   , & lepEnergyCorr   , & b_lepEnergyCorr);
  chain->SetBranchAddress("timestamp"   , & timestamp   , & b_timestamp);

  if (isMC) chain->SetBranchAddress("weightGenerator"   , & GenWeight   , & b_weightGenerator);

  if (isMC) chain->SetBranchAddress("vertMultiTrue", &vertMultiTrue_, &b_vertMultiTrue);


  string channels[3] = {"ee", "emu", "mumu"}; 
  string selections[4] = {"afterSelections", "dilepTriggers", "crossTriggers", "crossAndDilepTriggers"};
  string systematics[23] = {"nominal", "low2_njets", "high2_njets",  "low3_njets", "high3_njets",  "low4_njets", "high4_njets",  "low5_njets", "high5_njets", "low20_nvtx", "high20_nvtx","low25_nvtx", "high25_nvtx",  "low30_nvtx", "high30_nvtx", "low35_nvtx", "high35_nvtx",  "low40_nvtx", "high40_nvtx", "low45_nvtx", "high45_nvtx", "low50_nvtx", "high50_nvtx"};

  TH1F *h_met[3][4][23];
  TH1F *h_sidereal[3][4][23];
  TH1F *h_nvtx[3][4][23];
  TH1F *h_nvtx_good[3][4][23];
  TH1F *h_njets[3][4][23];
  TH1F *h_lepApt[3][4][23];
  TH1F *h_lepAeta[3][4][23];
  TH1F *h_lepAphi[3][4][23];
  TH1F *h_lepBpt[3][4][23];
  TH1F *h_lepBeta[3][4][23];
  TH1F *h_lepBphi[3][4][23];

  TH1F *h_lepApt_eeIn[3][4][23];
  TH1F *h_lepBpt_eeIn[3][4][23];
  TH1F *h_lepApt_eeSplit[3][4][23];
  TH1F *h_lepBpt_eeSplit[3][4][23];
  TH1F *h_lepApt_eeOut[3][4][23];
  TH1F *h_lepBpt_eeOut[3][4][23];
  TH1F *h_lepApt_emuIn[3][4][23];
  TH1F *h_lepBpt_emuIn[3][4][23];
  TH1F *h_lepApt_emuOut[3][4][23];
  TH1F *h_lepBpt_emuOut[3][4][23];

  TH2F *h_lepABeta_onebin[3][2][23];
  TH2F *h_lepABeta[3][2][23];
  TH2F *h_lepABpt[3][2][23];

  TH2F *h_lepABpt_eeIn[3][2][23];
  TH2F *h_lepABpt_eeSplit[3][2][23];
  TH2F *h_lepABpt_eeOut[3][2][23];
  TH2F *h_lepABpt_emuIn[3][2][23];
  TH2F *h_lepABpt_emuOut[3][2][23];

  Double_t siderealbins[24+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
  Double_t metbins[5+1] = {100,120,140,160,180,200};
  Double_t nvtxbins[12+1] = {0,5,10,15,20,25,30,35,40,45,50,55,60};
  Double_t njetsbins[7+1] = {0,1,2,3,4,5,6,7};
  Double_t lepAptbins[7+1] = {20,40,60,80,100,150,200,500};

  Double_t lepAptbins_eeIn[6+1] = {20,40,65,100,150,200,500};
  Double_t lepBptbins_eeIn[5+1] = {15,40,70,100,150,500};

  Double_t lepAptbins_eeSplit[5+1] = {20,40,60,100,200,500};
  Double_t lepBptbins_eeSplit[3+1] = {15,50,100,500};

  Double_t lepAptbins_eeOut[3+1] = {20,50,100,500};
  Double_t lepBptbins_eeOut[3+1] = {15,50,100,500};

  Double_t lepAptbins_emuIn[5+1] = {15,40,70,100,150,500};
  Double_t lepBptbins_emuIn[5+1] = {15,40,70,100,150,500};

  Double_t lepAptbins_emuOut[4+1] = {15,40,70,100,500};
  Double_t lepBptbins_emuOut[4+1] = {15,40,70,150,500};

  Double_t lepAetabins[10+1] = {-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5};
  Double_t lepAphibins[8+1] = {-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.};
  Double_t lepBptbins[8+1] = {15,30,45,60,80,100,150,200,500};
  Double_t lepBetabins[10+1] = {-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5};
  Double_t lepBphibins[8+1] = {-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.};
  Double_t lepABetabins[2+1] = {0.0,1.2,2.5};
  Double_t lepABetaonebin[1+1] = {0.0,2.5};

  for (unsigned short i(0); i < 3; i++) {
    for (unsigned short j(0); j < 4; j++) {
      for (unsigned short k(0); k < 23; k++) {
	h_sidereal[i][j][k] = new TH1F((channels[i]+"_sidereal_"+systematics[k]+"_"+selections[j]).c_str(),"",24,siderealbins);
	h_met[i][j][k] = new TH1F((channels[i]+"_met_"+systematics[k]+"_"+selections[j]).c_str(),"",5,metbins);
	h_nvtx[i][j][k] = new TH1F((channels[i]+"_nvtx_"+systematics[k]+"_"+selections[j]).c_str(),"",12,nvtxbins);
	h_nvtx_good[i][j][k] = new TH1F((channels[i]+"_nvtx_good_"+systematics[k]+"_"+selections[j]).c_str(),"",12,nvtxbins);
	h_njets[i][j][k] = new TH1F((channels[i]+"_njets_"+systematics[k]+"_"+selections[j]).c_str(),"",7,njetsbins);
	h_lepApt[i][j][k] = new TH1F((channels[i]+"_lepApt_"+systematics[k]+"_"+selections[j]).c_str(),"",7,lepAptbins);
	h_lepAeta[i][j][k] = new TH1F((channels[i]+"_lepAeta_"+systematics[k]+"_"+selections[j]).c_str(),"",10,lepAetabins);
	h_lepAphi[i][j][k] = new TH1F((channels[i]+"_lepAphi_"+systematics[k]+"_"+selections[j]).c_str(),"",8,lepAphibins);
	h_lepBpt[i][j][k] = new TH1F((channels[i]+"_lepBpt_"+systematics[k]+"_"+selections[j]).c_str(),"",8,lepBptbins);
	h_lepBeta[i][j][k] = new TH1F((channels[i]+"_lepBeta_"+systematics[k]+"_"+selections[j]).c_str(),"",10,lepBetabins);
	h_lepBphi[i][j][k] = new TH1F((channels[i]+"_lepBphi_"+systematics[k]+"_"+selections[j]).c_str(),"",8,lepBphibins);

	h_lepApt_eeIn[i][j][k] = new TH1F((channels[i]+"_lepApt_eeIn_"+systematics[k]+"_"+selections[j]).c_str(),"",6,lepAptbins_eeIn);
	h_lepBpt_eeIn[i][j][k] = new TH1F((channels[i]+"_lepBpt_eeIn_"+systematics[k]+"_"+selections[j]).c_str(),"",5,lepBptbins_eeIn);

	h_lepApt_eeOut[i][j][k] = new TH1F((channels[i]+"_lepApt_eeOut_"+systematics[k]+"_"+selections[j]).c_str(),"",3,lepAptbins_eeOut);
	h_lepBpt_eeOut[i][j][k] = new TH1F((channels[i]+"_lepBpt_eeOut_"+systematics[k]+"_"+selections[j]).c_str(),"",3,lepBptbins_eeOut);

	h_lepApt_eeSplit[i][j][k] = new TH1F((channels[i]+"_lepApt_eeSplit_"+systematics[k]+"_"+selections[j]).c_str(),"",5,lepAptbins_eeSplit);
	h_lepBpt_eeSplit[i][j][k] = new TH1F((channels[i]+"_lepBpt_eeSplit_"+systematics[k]+"_"+selections[j]).c_str(),"",3,lepBptbins_eeSplit);

	h_lepApt_emuIn[i][j][k] = new TH1F((channels[i]+"_lepApt_emuIn_"+systematics[k]+"_"+selections[j]).c_str(),"",5,lepAptbins_emuIn);
	h_lepBpt_emuIn[i][j][k] = new TH1F((channels[i]+"_lepBpt_emuIn_"+systematics[k]+"_"+selections[j]).c_str(),"",5,lepBptbins_emuIn);

	h_lepApt_emuOut[i][j][k] = new TH1F((channels[i]+"_lepApt_emuOut_"+systematics[k]+"_"+selections[j]).c_str(),"",4,lepAptbins_emuOut);
	h_lepBpt_emuOut[i][j][k] = new TH1F((channels[i]+"_lepBpt_emuOut_"+systematics[k]+"_"+selections[j]).c_str(),"",4,lepBptbins_emuOut);

	if (j>1){	  
	  h_lepABeta_onebin[i][j-2][k] = new TH2F((channels[i]+"_lepABeta_onebin_"+systematics[k]+"_"+selections[j]).c_str(), "", 1,lepABetaonebin,1,lepABetaonebin);
	  h_lepABeta[i][j-2][k] = new TH2F((channels[i]+"_lepABeta_"+systematics[k]+"_"+selections[j]).c_str(), "", 2,lepABetabins,2,lepABetabins);

	  if (i == 1) h_lepABpt[i][j-2][k] = new TH2F((channels[i]+"_lepABpt_"+systematics[k]+"_"+selections[j]).c_str(), "", 8,lepBptbins,8,lepBptbins);      
	  else h_lepABpt[i][j-2][k] = new TH2F((channels[i]+"_lepABpt_"+systematics[k]+"_"+selections[j]).c_str(), "", 7,lepAptbins,8,lepBptbins);      

	  h_lepABpt_eeIn[i][j-2][k] = new TH2F((channels[i]+"_lepABpt_eeIn_"+systematics[k]+"_"+selections[j]).c_str(), "", 6,lepAptbins_eeIn,5,lepBptbins_eeIn);      
	  h_lepABpt_eeOut[i][j-2][k] = new TH2F((channels[i]+"_lepABpt_eeOut_"+systematics[k]+"_"+selections[j]).c_str(), "", 3,lepAptbins_eeOut,3,lepBptbins_eeOut);      
	  h_lepABpt_eeSplit[i][j-2][k] = new TH2F((channels[i]+"_lepABpt_eeSplit_"+systematics[k]+"_"+selections[j]).c_str(), "", 5,lepAptbins_eeSplit,3,lepBptbins_eeSplit);      

	  h_lepABpt_emuIn[i][j-2][k] = new TH2F((channels[i]+"_lepABpt_emuIn_"+systematics[k]+"_"+selections[j]).c_str(), "", 5,lepAptbins_emuIn,5,lepBptbins_emuIn);      
	  h_lepABpt_emuOut[i][j-2][k] = new TH2F((channels[i]+"_lepABpt_emuOut_"+systematics[k]+"_"+selections[j]).c_str(), "", 4,lepAptbins_emuOut,4,lepBptbins_emuOut);      

	}
      }
    }
  }

  TH1F *h_LepSF = new TH1F(("LepSF"+suffix).c_str(), "Lep. Id and Isol. SF per event", 200, 0.75, 1.25);
  TH1F *h_PUSF = new TH1F(("PUSF"+suffix).c_str(), "PU SF per event", 200, 0.5, 1.5);



  /// for 2016 
  if (era_ == "2016") {

    t0 = 1451606400; // Reference le 1 janv 2016

    SF_electron_ID_File = "electron/94X_2016/SF_cutbasedID/2016LegacyReReco_ElectronTight_Fall17V2.root";
    SF_electron_Reco_File   = "electron/94X_2016/SF_reco/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root"; 
    SF_electron_Reco_File_LowPt   = "electron/94X_2016/SF_reco/EGM2D_BtoH_low_RecoSF_Legacy2016.root";

    SF_muon_ID_File       = "muon/94X_2016/SF_IDnISO/RunBCDEFGH_SF_ID.root";
    SF_muon_ID_Histo      = "NUM_TightID_DEN_genTracks_eta_pt";
    SF_muon_ID_Format     = "pt_vs_eta";
    SF_muon_ID_Histo_Stat = "NUM_TightID_DEN_genTracks_eta_pt_stat";
    SF_muon_ID_Histo_Sys  = "NUM_TightID_DEN_genTracks_eta_pt_syst";

    SF_muon_ISO_File       = "muon/94X_2016/SF_IDnISO/RunBCDEFGH_SF_ISO.root";
    SF_muon_ISO_Histo      = "NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt";
    SF_muon_ISO_Format     = "pt_vs_eta";
    SF_muon_ISO_Histo_Stat = "NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt_stat";
    SF_muon_ISO_Histo_Sys  = "NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt_syst";

    pileup_DATA_File = "pileup/pileupDATA__Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON__max120_n120__MinBias69200.root";

    HLTPathsOR_ee = {    
      "HLT_Ele27_WPTight_Gsf_v*",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
      "HLT_DoubleEle33_CaloIdL_MW_v*",
      "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*",
    };

    HLTPathsOR_emu = {
      "HLT_IsoMu24_v*",   
      "HLT_IsoTkMu24_v*",  
      "HLT_Ele27_WPTight_Gsf_v*",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", 
    };

    HLTPathsOR_mumu = {  
      "HLT_IsoMu24_v*",   
      "HLT_IsoTkMu24_v*", 
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*", 
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
    };

    HLTPathsOR_met = {    
      "HLT_PFMET300_v*",  
      "HLT_MET200_v*",    
      "HLT_PFHT300_PFMET110_v*",    
      "HLT_PFMET170_HBHECleaned_v*",
      "HLT_PFMET120_PFMHT120_IDTight_v*",
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*",
    };
  }
  /// 

  /// for 2017
  else if (era_ == "2017") {

    t0 = 1451606400; // Reference le 1 janv 2016
    //    t0 = 1483228800; // Reference le 1 janv 2017

    SF_electron_ID_File = "electron/94X_2017/SF_cutbasedID/v2/2017_ElectronTight.root";
    SF_electron_Reco_File   = "electron/94X_2017/SF_reco/v2/egammaEffi.txt_EGM2D.root"; 
    SF_electron_Reco_File_LowPt   = "electron/94X_2017/SF_reco/v2/egammaEffi.txt_EGM2D_low.root";

    SF_muon_ID_File       = "muon/94X_2017/SF_IDnISO/RunBCDEF_SF_ID_syst.root";
    SF_muon_ID_Histo      = "NUM_TightID_DEN_genTracks_pt_abseta";
    SF_muon_ID_Format     = "absEta_vs_pt";
    SF_muon_ID_Histo_Stat = "NUM_TightID_DEN_genTracks_pt_abseta_stat";
    SF_muon_ID_Histo_Sys  = "NUM_TightID_DEN_genTracks_pt_abseta_syst";

    SF_muon_ISO_File       = "muon/94X_2017/SF_IDnISO/RunBCDEF_SF_ISO_syst.root";
    SF_muon_ISO_Histo      = "NUM_LooseRelIso_DEN_TightIDandIPCut_pt_abseta";
    SF_muon_ISO_Format     = "absEta_vs_pt";
    SF_muon_ISO_Histo_Stat = "NUM_LooseRelIso_DEN_TightIDandIPCut_pt_abseta_stat";
    SF_muon_ISO_Histo_Sys  = "NUM_LooseRelIso_DEN_TightIDandIPCut_pt_abseta_syst";

    pileup_DATA_File = "pileup/pileupDATA__Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1__max120_n120__MinBias69200.root";

    HLTPathsOR_ee = {    
      "HLT_Ele35_WPTight_Gsf_v*",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
      "HLT_DoubleEle33_CaloIdL_MW_v*",
    };

    HLTPathsOR_emu = {
      "HLT_IsoMu27_v*",
      "HLT_Ele35_WPTight_Gsf_v*",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",
    };

    HLTPathsOR_mumu = {  
      "HLT_IsoMu27_v*",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*",
    };

    HLTPathsOR_met = {             
      "HLT_PFHT500_PFMET100_PFMHT100_IDTight_v*",
      "HLT_PFHT700_PFMET85_PFMHT85_IDTight_v*",
      "HLT_PFHT800_PFMET75_PFMHT75_IDTight_v*",
      "HLT_PFMET120_PFMHT120_IDTight_v*",
      "HLT_PFMET200_HBHECleaned_v*",
      "HLT_PFMET200_HBHE_BeamHaloCleaned_v*",
      "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v*",
      "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v*",
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*",
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v*",
    };

  }
  /// 

  /// for 2018 
  else if (era_ == "2018") {

    t0 = 1451606400; // Reference le 1 janv 2016

    SF_electron_ID_File = "electron/102X_2018/Autumn18/SF_cutbasedID/2018_ElectronTight.root";
    SF_electron_Reco_File   = "electron/102X_2018/Autumn18/SF_reco/egammaEffi.txt_EGM2D_updatedAll.root"; 
    SF_electron_Reco_File_LowPt   = "electron/102X_2018/Autumn18/SF_reco/egammaEffi.txt_EGM2D_updatedAll.root";

    SF_muon_ID_File       = "muon/102X_2018/Autumn18/SF_IDnISO/RunABCD_SF_ID_sys.root";
    SF_muon_ID_Histo      = "NUM_TightID_DEN_TrackerMuons_abseta_pt";
    SF_muon_ID_Format     = "pt_vs_absEta";
    SF_muon_ID_Histo_Stat = "NUM_TightID_DEN_TrackerMuons_abseta_pt_stat";
    SF_muon_ID_Histo_Sys  = "NUM_TightID_DEN_TrackerMuons_abseta_pt_syst";

    SF_muon_ISO_File       = "muon/102X_2018/Autumn18/SF_IDnISO/RunABCD_SF_ISO_sys.root";
    SF_muon_ISO_Histo      = "NUM_LooseRelIso_DEN_TightIDandIPCut_abseta_pt";
    SF_muon_ISO_Format     = "pt_vs_absEta";
    SF_muon_ISO_Histo_Stat = "NUM_LooseRelIso_DEN_TightIDandIPCut_abseta_pt_stat";
    SF_muon_ISO_Histo_Sys  = "NUM_LooseRelIso_DEN_TightIDandIPCut_abseta_pt_syst";

    pileup_DATA_File = "pileup/pileupDATA__Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON__max120_n120__MinBias69200.root";

    HLTPathsOR_ee = {    
      "HLT_Ele32_WPTight_Gsf_v*",  
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
      "HLT_DoubleEle25_CaloIdL_MW_v*",
    };

    HLTPathsOR_emu = {
      "HLT_IsoMu24_v*",          
      "HLT_Ele32_WPTight_Gsf_v*", 
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",
    };

    HLTPathsOR_mumu = {  
      "HLT_IsoMu24_v*",  
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*",
    };

    HLTPathsOR_met = {             
      "HLT_PFHT500_PFMET100_PFMHT100_IDTight_v*",
      "HLT_PFHT700_PFMET85_PFMHT85_IDTight_v*",
      "HLT_PFHT800_PFMET75_PFMHT75_IDTight_v*",
      "HLT_PFMET120_PFMHT120_IDTight_v*",
      "HLT_PFMET200_HBHECleaned_v*",
      "HLT_PFMET200_HBHE_BeamHaloCleaned_v*",
      "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v*",
      "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v*",
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*",
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v*",
    };

  }
  /// 

  /// for 2016preVFP UL
  if (era_ == "2016preVFP_UL") {

    t0 = 1451606400; // Reference le 1 janv 2016

    SF_electron_ID_File = "electron/106x_2016_preVFP/SF_cutbasedID/egammaEffi.txt_Ele_Tight_preVFP_EGM2D.root";
    SF_electron_Reco_File   = "electron/106x_2016_preVFP/SF_reco/egammaEffi_ptAbove20.txt_EGM2D_UL2016preVFP.root";
    SF_electron_Reco_File_LowPt   = "electron/106x_2016_preVFP/SF_reco/egammaEffi_ptBelow20.txt_EGM2D_UL2016preVFP.root"; 

    SF_muon_ID_File       = "muon/106X_2016_preVFP/SF_IDnISO/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root";
    SF_muon_ID_Histo      = "NUM_TightID_DEN_TrackerMuons_abseta_pt";
    SF_muon_ID_Format     = "pt_vs_absEta";
    SF_muon_ID_Histo_Stat = "NUM_TightID_DEN_TrackerMuons_abseta_pt_stat";
    SF_muon_ID_Histo_Sys  = "NUM_TightID_DEN_TrackerMuons_abseta_pt_syst";

    SF_muon_ISO_File       = "muon/106X_2016_preVFP/SF_IDnISO/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root";
    SF_muon_ISO_Histo      = "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt";
    SF_muon_ISO_Format     = "pt_vs_absEta";
    SF_muon_ISO_Histo_Stat = "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_stat";
    SF_muon_ISO_Histo_Sys  = "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_syst";

    pileup_DATA_File = "pileup/pileupDATA__Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON__max120_n120__MinBias69200.root";

    HLTPathsOR_ee = {    
      "HLT_Ele27_WPTight_Gsf_v*",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
      "HLT_DoubleEle33_CaloIdL_MW_v*",
      "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*",
    };

    HLTPathsOR_emu = {
      "HLT_IsoMu24_v*",   
      "HLT_IsoTkMu24_v*",  
      "HLT_Ele27_WPTight_Gsf_v*",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
    };

    HLTPathsOR_mumu = {  
      "HLT_IsoMu24_v*",   
      "HLT_IsoTkMu24_v*", 
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*", 
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
    };

    HLTPathsOR_met = {    
      "HLT_PFMET300_v*",  
      "HLT_MET200_v*",    
      "HLT_PFHT300_PFMET110_v*",    
      "HLT_PFMET170_HBHECleaned_v*",
      "HLT_PFMET120_PFMHT120_IDTight_v*",
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*",
    };
  }
  /// 


  /// for 2016postVFP UL
  if (era_ == "2016postVFP_UL") {

    t0 = 1451606400; // Reference le 1 janv 2016

    SF_electron_ID_File = "electron/106x_2016_postVFP/SF_cutbasedID/egammaEffi.txt_Ele_Tight_postVFP_EGM2D.root";
    SF_electron_Reco_File   = "electron/106x_2016_postVFP/SF_reco/egammaEffi_ptAbove20.txt_EGM2D_UL2016postVFP.root"; 
    SF_electron_Reco_File_LowPt   = "electron/106x_2016_postVFP/SF_reco/egammaEffi_ptBelow20.txt_EGM2D_UL2016postVFP.root";

    SF_muon_ID_File       = "muon/106X_2016_postVFP/SF_IDnISO/Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root";
    SF_muon_ID_Histo      = "NUM_TightID_DEN_TrackerMuons_abseta_pt";
    SF_muon_ID_Format     = "pt_vs_absEta";
    SF_muon_ID_Histo_Stat = "NUM_TightID_DEN_TrackerMuons_abseta_pt_stat";
    SF_muon_ID_Histo_Sys  = "NUM_TightID_DEN_TrackerMuons_abseta_pt_syst";

    SF_muon_ISO_File       = "muon/106X_2016_postVFP/SF_IDnISO/Efficiencies_muon_generalTracks_Z_Run2016_UL_ISO.root";
    SF_muon_ISO_Histo      = "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt";
    SF_muon_ISO_Format     = "pt_vs_absEta";
    SF_muon_ISO_Histo_Stat = "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_stat";
    SF_muon_ISO_Histo_Sys  = "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_syst";

    pileup_DATA_File = "pileup/pileupDATA__Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON__max120_n120__MinBias69200.root";

    HLTPathsOR_ee = {    
      "HLT_Ele27_WPTight_Gsf_v*",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
      "HLT_DoubleEle33_CaloIdL_MW_v*",
      "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*",
    };

    HLTPathsOR_emu = {
      "HLT_IsoMu24_v*",   
      "HLT_IsoTkMu24_v*",  
      "HLT_Ele27_WPTight_Gsf_v*",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", 
    };

    HLTPathsOR_mumu = {  
      "HLT_IsoMu24_v*",   
      "HLT_IsoTkMu24_v*", 
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*", 
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",
    };

    if(sample_ == "run2016H") {
      HLTPathsOR_mumu = {  
	"HLT_IsoMu24_v*",   
	"HLT_IsoTkMu24_v*", 
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
	"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
      };
    }

    HLTPathsOR_met = {    
      "HLT_PFMET300_v*",  
      "HLT_MET200_v*",    
      "HLT_PFHT300_PFMET110_v*",    
      "HLT_PFMET170_HBHECleaned_v*",
      "HLT_PFMET120_PFMHT120_IDTight_v*",
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*",
    };
  }
  /// 

  /// for 2017 UL
  else if (era_ == "2017_UL") {

    t0 = 1451606400; // Reference le 1 janv 2016

    SF_electron_Reco_File   = "electron/106X_2017/SF_reco/egammaEffi_ptAbove20.txt_EGM2D_UL2017.root"; 
    SF_electron_Reco_File_LowPt   = "electron/106X_2017/SF_reco/egammaEffi_ptBelow20.txt_EGM2D_UL2017.root";
    SF_electron_ID_File = "electron/106X_2017/SF_cutbasedID/egammaEffi.txt_EGM2D_Tight_UL17.root";


    SF_muon_ID_File       = "muon/106X_2017/SF_IDnISO/Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root";
    SF_muon_ID_Histo      = "NUM_TightID_DEN_TrackerMuons_abseta_pt";
    SF_muon_ID_Format     = "pt_vs_absEta";
    SF_muon_ID_Histo_Stat = "NUM_TightID_DEN_TrackerMuons_abseta_pt_stat";
    SF_muon_ID_Histo_Sys  = "NUM_TightID_DEN_TrackerMuons_abseta_pt_syst";

    SF_muon_ISO_File       = "muon/106X_2017/SF_IDnISO/Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root";
    SF_muon_ISO_Histo      = "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt";
    SF_muon_ISO_Format     = "pt_vs_absEta";
    SF_muon_ISO_Histo_Stat = "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_stat";
    SF_muon_ISO_Histo_Sys  = "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_syst";

    pileup_DATA_File = "pileup/pileupDATA__Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON__max120_n120__MinBias69200.root";

    HLTPathsOR_ee = {    
      "HLT_Ele35_WPTight_Gsf_v*",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
      "HLT_DoubleEle33_CaloIdL_MW_v*",
    };

    HLTPathsOR_emu = {
      "HLT_IsoMu27_v*",
      "HLT_Ele35_WPTight_Gsf_v*",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",
    };

    HLTPathsOR_mumu = {  
      "HLT_IsoMu27_v*",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*",
    };

    if(sample_ == "run2017B") {
      HLTPathsOR_mumu = {  
	"HLT_IsoMu27_v*",
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
      };
    }

    HLTPathsOR_met = {             
      "HLT_PFHT500_PFMET100_PFMHT100_IDTight_v*",
      "HLT_PFHT700_PFMET85_PFMHT85_IDTight_v*",
      "HLT_PFHT800_PFMET75_PFMHT75_IDTight_v*",
      "HLT_PFMET120_PFMHT120_IDTight_v*",
      "HLT_PFMET200_HBHECleaned_v*",
      "HLT_PFMET200_HBHE_BeamHaloCleaned_v*",
      "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v*",
      "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v*",
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*",
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v*",
    };

  }
  /// 


  /// for 2018 UL
  else if (era_ == "2018_UL") {

    t0 = 1451606400; // Reference le 1 janv 2016

    SF_electron_ID_File = "electron/106X_2018/SF_cutbasedID/egammaEffi.txt_Ele_Tight_EGM2D.root";
    SF_electron_Reco_File   = "electron/106X_2018/SF_reco/egammaEffi_ptAbove20.txt_EGM2D_UL2018.root"; 
    SF_electron_Reco_File_LowPt   = "electron/106X_2018/SF_reco/egammaEffi_ptBelow20.txt_EGM2D_UL2018.root";

    SF_muon_ID_File       = "muon/106X_2018/SF_IDnISO/Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root";
    SF_muon_ID_Histo      = "NUM_TightID_DEN_TrackerMuons_abseta_pt";
    SF_muon_ID_Format     = "pt_vs_absEta";
    SF_muon_ID_Histo_Stat = "NUM_TightID_DEN_TrackerMuons_abseta_pt_stat";
    SF_muon_ID_Histo_Sys  = "NUM_TightID_DEN_TrackerMuons_abseta_pt_syst";

    SF_muon_ISO_File       = "muon/106X_2018/SF_IDnISO/Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root";
    SF_muon_ISO_Histo      = "NUM_LooseRelIso_DEN_TightIDandIPCut_abseta_pt";
    SF_muon_ISO_Format     = "pt_vs_absEta";
    SF_muon_ISO_Histo_Stat = "NUM_LooseRelIso_DEN_TightIDandIPCut_abseta_pt_stat";
    SF_muon_ISO_Histo_Sys  = "NUM_LooseRelIso_DEN_TightIDandIPCut_abseta_pt_syst";

    pileup_DATA_File = "pileup/pileupDATA__Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON__max120_n120__MinBias69200.root";

    HLTPathsOR_ee = {    
      "HLT_Ele32_WPTight_Gsf_v*",  
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
      "HLT_DoubleEle25_CaloIdL_MW_v*",
    };

    HLTPathsOR_emu = {
      "HLT_IsoMu24_v*",          
      "HLT_Ele32_WPTight_Gsf_v*", 
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",
    };

    HLTPathsOR_mumu = {  
      "HLT_IsoMu24_v*",  
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*",
    };

    HLTPathsOR_met = {             
      "HLT_PFHT500_PFMET100_PFMHT100_IDTight_v*",
      "HLT_PFHT700_PFMET85_PFMHT85_IDTight_v*",
      "HLT_PFHT800_PFMET75_PFMHT75_IDTight_v*",
      "HLT_PFMET120_PFMHT120_IDTight_v*",
      "HLT_PFMET200_HBHECleaned_v*",
      "HLT_PFMET200_HBHE_BeamHaloCleaned_v*",
      "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v*",
      "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v*",
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*",
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v*",
    };

  }
  /// 


/// SF implementation 
  if (isMC) {
    prepareLeptonSF(); 
    preparePUSF(filename);
  }
/// 


  float nEventsPassingCrossTriggers_ee = 0.; 
  float nEventsPassingCrossTriggersAndDilepTriggers_ee = 0.;
  float nEventsPassingCrossTriggers_emu = 0.; 
  float nEventsPassingCrossTriggersAndDilepTriggers_emu = 0.;
  float nEventsPassingCrossTriggers_mumu = 0.; 
  float nEventsPassingCrossTriggersAndDilepTriggers_mumu = 0.;


  Long64_t nEvents = chain->GetEntries();
  //  nEvents = 10000; 

  //Event loop                                                                                                                                                           
  for (Long64_t i = 0; i < nEvents; i++) {

    int channelPdgIdProduct = 0;
    if (i % (nEvents/10) == 0) cout << "=== Event " << i/(nEvents/10) * 10 << "%" << endl;
      //  cout << "=== Event " << i  << endl;

    chain->GetEntry(i);
    
    if(v_leptons->size() < 2) continue;
    if(sign(lepPdgId->at(0)) == sign(lepPdgId->at(1))) continue;   

    for (unsigned int k = 0; k < v_leptons->size(); k++) {

      v_leptons->at(k) *= lepEnergyCorr->at(k);
    }
    
    if(v_leptons->at(0).Pt() < v_leptons->at(1).Pt()) {
      leadinglepton = 1;
      subleadinglepton = 0;
      //      cout << "Impossible event with subleading pT greater than leading pT... " << endl;
      //      cout << "Leading pt " << v_leptons->at(0).Pt() << endl;
      //      cout << "Subleading pt " << v_leptons->at(1).Pt() << endl;
      //      continue;
    }
    else {
      leadinglepton = 0;
      subleadinglepton = 1;
    }

    if(v_leptons->at(leadinglepton).Pt() < 25.0 || v_leptons->at(subleadinglepton).Pt() < 15.0) continue;
    channelPdgIdProduct = lepPdgId->at(leadinglepton)*lepPdgId->at(subleadinglepton);
    if(abs(lepPdgId->at(leadinglepton)) == 11) { if(lepID_ElecCutBased->at(leadinglepton) != 4) continue; }
    if(abs(lepPdgId->at(subleadinglepton)) == 11) { if(lepID_ElecCutBased->at(subleadinglepton) != 4) continue; }
    if(abs(lepPdgId->at(leadinglepton)) == 13) { if(lepID_MuonTight->at(leadinglepton) != 1) continue; }
    if(abs(lepPdgId->at(subleadinglepton)) == 13) { if(lepID_MuonTight->at(subleadinglepton) != 1) continue; }
    if(abs(lepPdgId->at(leadinglepton)) == 13) { if(lepPfIso->at(leadinglepton) > 0.15) continue; }
    if(abs(lepPdgId->at(subleadinglepton)) == 13) { if(lepPfIso->at(subleadinglepton) > 0.15) continue; }
    if(abs(lepPdgId->at(leadinglepton)) == 13) { if(fabs(v_leptons->at(leadinglepton).Eta()) > 2.4) continue;}
    if(abs(lepPdgId->at(subleadinglepton)) == 13) { if(fabs(v_leptons->at(subleadinglepton).Eta()) > 2.4) continue;}
    if(abs(lepPdgId->at(leadinglepton)) == 11) { if(fabs(lepSCEta->at(leadinglepton)) > 2.4) continue; }
    if(abs(lepPdgId->at(subleadinglepton)) == 11) { if(fabs(lepSCEta->at(subleadinglepton)) > 2.4) continue; }
    if(abs(lepPdgId->at(leadinglepton)) == 11) { if(fabs(lepSCEta->at(leadinglepton)) > 1.4442 && fabs(lepSCEta->at(leadinglepton)) < 1.5660) continue; }
    if(abs(lepPdgId->at(subleadinglepton)) == 11) { if(fabs(lepSCEta->at(subleadinglepton)) > 1.4442 && fabs(lepSCEta->at(subleadinglepton)) < 1.5660) continue; }
    LeptonA.SetPtEtaPhiM(v_leptons->at(leadinglepton).Pt(), v_leptons->at(leadinglepton).Eta(), v_leptons->at(leadinglepton).Phi(), v_leptons->at(leadinglepton).M());
    LeptonB.SetPtEtaPhiM(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(subleadinglepton).Eta(), v_leptons->at(subleadinglepton).Phi(), v_leptons->at(subleadinglepton).M());
    if((LeptonA + LeptonB).M() < 20) continue;
    if(met->Pt() < 100) continue;
    if(metFilters == 1) continue;



    HLTBitsManager_.set_results(*HLTBits);



    int bjets = 0;    
    int njets = 0;
    for (unsigned int k = 0; k < v_jets->size(); k++) {
      if(fabs(v_jets->at(k).Eta()) > 2.4) continue;
      if(v_jets->at(k).Pt() < 30.0) continue;
      if(jetPFID->at(k) != 3) continue;
      if(deltaR(v_leptons->at(leadinglepton).Eta(),v_leptons->at(leadinglepton).Phi(),v_jets->at(k).Eta(),v_jets->at(k).Phi())  < 0.4) continue;
      if(deltaR(v_leptons->at(subleadinglepton).Eta(),v_leptons->at(subleadinglepton).Phi(),v_jets->at(k).Eta(),v_jets->at(k).Phi())  < 0.4) continue;
      //      cout << "jet " << k  << endl;
      //      cout << "  eta " << v_jets->at(k).Eta()  << endl;
      //      cout << "  pt " << v_jets->at(k).Pt()  << endl;
      njets = njets + 1;

      if (era_ == "2016") {
	if(jetBTagDeepCSV->at(k) > 0.2217) bjets = bjets + 1;
      }
      else if (era_ == "2017") {
	if(jetBTagDeepCSV->at(k) > 0.1522) bjets = bjets + 1;
      }
      else if (era_ == "2018") {
	if(jetBTagDeepCSV->at(k) > 0.1241) bjets = bjets + 1;
      }
      else if (era_ == "2016preVFP_UL") {
	//	if(jetBTagDeepCSV->at(k) > 0.2027) bjets = bjets + 1;
	if(jetBTagDeepJet->at(k) > 0.2598) bjets = bjets + 1;
      }
      else if (era_ == "2016postVFP_UL") {
	//	if(jetBTagDeepCSV->at(k) > 0.1918) bjets = bjets + 1;
	if(jetBTagDeepJet->at(k) > 0.2489) bjets = bjets + 1;
      }
      else if (era_ == "2017_UL") {
	//	if(jetBTagDeepCSV->at(k) > 0.1355) bjets = bjets + 1;
	if(jetBTagDeepJet->at(k) > 0.3040) bjets = bjets + 1;
      }
      else if (era_ == "2018_UL") {
	//	if(jetBTagDeepCSV->at(k) > 0.1208) bjets = bjets + 1;
	if(jetBTagDeepJet->at(k) > 0.2783) bjets = bjets + 1;
      }

    }
    if(bjets < 1) continue;


    if(MCorData_ == "Data") {
      t_sidereal = GetSiderealTime( Omega_UNIX, timestamp,  t0,  phi_UNIX,  phi_longitude, Omega_Sidereal);
      hour_sidereal = (t_sidereal / 3600);
      hour_sidereal = std::fmod(hour_sidereal, 24);

      //      cout<<"timestamp = " << timestamp << endl;
      //      cout<<"hour_sidereal = " << hour_sidereal << endl;
    }




/// SF implementation 
    const double weightLeptonSF = getWeightLeptonSF(0, 1, v_leptons, lepPdgId, lepSCEta);
    // cout<<"weightLeptonSF = " << weightLeptonSF << endl;
    h_LepSF->Fill(weightLeptonSF, 1);

    const double weightPU = getWeightPileup(i); 
    // cout << "weightPU = " << weightPU << endl;
    h_PUSF->Fill(weightPU, 1);


    double weight = 1;

    if (isMC) weight = GenWeight*weightLeptonSF*weightPU; 
    // cout << "weight = " << weight << endl;

///


    if (channelPdgIdProduct == (-11* 11)) { 
      // cout << "Check \"HLTPathsOR_ee\" : " << HLTBitsManager_.passes_OR(HLTPathsOR_ee) << endl; 
      // cout << "Check \"HLTPathsOR_met\" : " << HLTBitsManager_.passes_OR(HLTPathsOR_met) << endl; 
      h_sidereal[0][0][0]->Fill(hour_sidereal, weight);
      h_met[0][0][0]->Fill(met->Pt(), weight);
      h_nvtx[0][0][0]->Fill(NPV_all, weight);
      h_nvtx_good[0][0][0]->Fill(NPV_good, weight);
      h_njets[0][0][0]->Fill(njets, weight);
      h_lepApt[0][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
      h_lepBpt[0][0][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeIn[0][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeIn[0][0][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeSplit[0][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeSplit[0][0][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeSplit[0][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeSplit[0][0][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeOut[0][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeOut[0][0][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
      h_lepAeta[0][0][0]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
      h_lepBeta[0][0][0]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
      h_lepAphi[0][0][0]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
      h_lepBphi[0][0][0]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);

      if(HLTBitsManager_.passes_OR(HLTPathsOR_ee) == 1) {
	h_sidereal[0][1][0]->Fill(hour_sidereal, weight);
        h_met[0][1][0]->Fill(met->Pt(), weight);
        h_nvtx[0][1][0]->Fill(NPV_all, weight);
        h_nvtx_good[0][1][0]->Fill(NPV_good, weight);
        h_njets[0][1][0]->Fill(njets, weight);
        h_lepApt[0][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
        h_lepBpt[0][1][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeIn[0][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeIn[0][1][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeSplit[0][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeSplit[0][1][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeSplit[0][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeSplit[0][1][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeOut[0][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeOut[0][1][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
        h_lepAeta[0][1][0]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
        h_lepBeta[0][1][0]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
        h_lepAphi[0][1][0]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
        h_lepBphi[0][1][0]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
      }
      if(HLTBitsManager_.passes_OR(HLTPathsOR_met) == 1) {
        nEventsPassingCrossTriggers_ee++;
	h_sidereal[0][2][0]->Fill(hour_sidereal, weight);
        h_met[0][2][0]->Fill(met->Pt(), weight);
        h_nvtx[0][2][0]->Fill(NPV_all, weight);
        h_nvtx_good[0][2][0]->Fill(NPV_good, weight);
        h_njets[0][2][0]->Fill(njets, weight);
        h_lepApt[0][2][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
        h_lepBpt[0][2][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeIn[0][2][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeIn[0][2][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeSplit[0][2][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeSplit[0][2][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeSplit[0][2][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeSplit[0][2][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeOut[0][2][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeOut[0][2][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
        h_lepAeta[0][2][0]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
        h_lepBeta[0][2][0]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
        h_lepAphi[0][2][0]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
        h_lepBphi[0][2][0]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
        h_lepABpt[0][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
        h_lepABeta[0][0][0]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
        h_lepABeta_onebin[0][0][0]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeIn[0][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeSplit[0][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeSplit[0][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeOut[0][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);


	for (int k = 2; k < 6; k++) {

	  if (njets < k) { 
	    h_sidereal[0][2][2*k-3]->Fill(hour_sidereal, weight);
	    h_met[0][2][2*k-3]->Fill(met->Pt(), weight);
	    h_nvtx[0][2][2*k-3]->Fill(NPV_all, weight);
	    h_nvtx_good[0][2][2*k-3]->Fill(NPV_good, weight);
	    h_njets[0][2][2*k-3]->Fill(njets, weight);
	    h_lepApt[0][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	    h_lepBpt[0][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeIn[0][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeIn[0][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeSplit[0][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeSplit[0][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeSplit[0][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeSplit[0][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeOut[0][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeOut[0][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepAeta[0][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	    h_lepBeta[0][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	    h_lepAphi[0][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	    h_lepBphi[0][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	    h_lepABpt[0][0][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepABeta[0][0][2*k-3]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    h_lepABeta_onebin[0][0][2*k-3]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeIn[0][0][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeSplit[0][0][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeSplit[0][0][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeOut[0][0][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	  }
	  else {
	    h_sidereal[0][2][2*k-2]->Fill(hour_sidereal, weight);
	    h_met[0][2][2*k-2]->Fill(met->Pt(), weight);
	    h_nvtx[0][2][2*k-2]->Fill(NPV_all, weight);
	    h_nvtx_good[0][2][2*k-2]->Fill(NPV_good, weight);
	    h_njets[0][2][2*k-2]->Fill(njets, weight);
	    h_lepApt[0][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	    h_lepBpt[0][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeIn[0][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeIn[0][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeSplit[0][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeSplit[0][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeSplit[0][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeSplit[0][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeOut[0][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeOut[0][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepAeta[0][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	    h_lepBeta[0][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	    h_lepAphi[0][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	    h_lepBphi[0][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	    h_lepABpt[0][0][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepABeta[0][0][2*k-2]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    h_lepABeta_onebin[0][0][2*k-2]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeIn[0][0][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeSplit[0][0][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeSplit[0][0][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeOut[0][0][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	  }

	}

	for (int k = 2; k < 9; k++) {

	  if (NPV_all < 5*k + 10) {
	    h_sidereal[0][2][2*k+5]->Fill(hour_sidereal, weight);
	    h_met[0][2][2*k+5]->Fill(met->Pt(), weight);
	    h_nvtx[0][2][2*k+5]->Fill(NPV_all, weight);
	    h_nvtx_good[0][2][2*k+5]->Fill(NPV_good, weight);
	    h_njets[0][2][2*k+5]->Fill(njets, weight);
	    h_lepApt[0][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	    h_lepBpt[0][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeIn[0][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeIn[0][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeSplit[0][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeSplit[0][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeSplit[0][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeSplit[0][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeOut[0][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeOut[0][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepAeta[0][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	    h_lepBeta[0][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	    h_lepAphi[0][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	    h_lepBphi[0][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	    h_lepABpt[0][0][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepABeta[0][0][2*k+5]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    h_lepABeta_onebin[0][0][2*k+5]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeIn[0][0][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeSplit[0][0][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeSplit[0][0][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeOut[0][0][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	  }
	  else { 
	    h_sidereal[0][2][2*k+6]->Fill(hour_sidereal, weight);
	    h_met[0][2][2*k+6]->Fill(met->Pt(), weight);
	    h_nvtx[0][2][2*k+6]->Fill(NPV_all, weight);
	    h_nvtx_good[0][2][2*k+6]->Fill(NPV_good, weight);
	    h_njets[0][2][2*k+6]->Fill(njets, weight);
	    h_lepApt[0][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	    h_lepBpt[0][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeIn[0][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeIn[0][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeSplit[0][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeSplit[0][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeSplit[0][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeSplit[0][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeOut[0][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeOut[0][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepAeta[0][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	    h_lepBeta[0][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	    h_lepAphi[0][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	    h_lepBphi[0][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	    h_lepABpt[0][0][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepABeta[0][0][2*k+6]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    h_lepABeta_onebin[0][0][2*k+6]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeIn[0][0][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeSplit[0][0][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeSplit[0][0][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeOut[0][0][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	  }

	}

      }
      
   
      if(HLTBitsManager_.passes_OR(HLTPathsOR_ee) == 1 && HLTBitsManager_.passes_OR(HLTPathsOR_met) == 1) {
        nEventsPassingCrossTriggersAndDilepTriggers_ee++;
	h_sidereal[0][3][0]->Fill(hour_sidereal, weight);
        h_met[0][3][0]->Fill(met->Pt(), weight);
        h_nvtx[0][3][0]->Fill(NPV_all, weight);
        h_nvtx_good[0][3][0]->Fill(NPV_good, weight);
        h_njets[0][3][0]->Fill(njets, weight);
        h_lepApt[0][3][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
        h_lepBpt[0][3][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeIn[0][3][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeIn[0][3][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeSplit[0][3][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeSplit[0][3][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeSplit[0][3][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeSplit[0][3][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeOut[0][3][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeOut[0][3][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
        h_lepAeta[0][3][0]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
        h_lepBeta[0][3][0]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
        h_lepAphi[0][3][0]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
        h_lepBphi[0][3][0]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
        h_lepABpt[0][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
        h_lepABeta[0][1][0]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);        
        h_lepABeta_onebin[0][1][0]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeIn[0][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeSplit[0][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeSplit[0][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeOut[0][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);

	for (int k = 2; k < 6; k++) {

	  if (njets < k) { 
	    h_sidereal[0][3][2*k-3]->Fill(hour_sidereal, weight);
	    h_met[0][3][2*k-3]->Fill(met->Pt(), weight);
	    h_nvtx[0][3][2*k-3]->Fill(NPV_all, weight);
	    h_nvtx_good[0][3][2*k-3]->Fill(NPV_good, weight);
	    h_njets[0][3][2*k-3]->Fill(njets, weight);
	    h_lepApt[0][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	    h_lepBpt[0][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeIn[0][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeIn[0][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeSplit[0][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeSplit[0][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeSplit[0][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeSplit[0][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeOut[0][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeOut[0][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepAeta[0][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	    h_lepBeta[0][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	    h_lepAphi[0][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	    h_lepBphi[0][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	    h_lepABpt[0][1][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepABeta[0][1][2*k-3]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    h_lepABeta_onebin[0][1][2*k-3]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeIn[0][1][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeSplit[0][1][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeSplit[0][1][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeOut[0][1][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	  }
	  else { 
	    h_sidereal[0][3][2*k-2]->Fill(hour_sidereal, weight);
	    h_met[0][3][2*k-2]->Fill(met->Pt(), weight);
	    h_nvtx[0][3][2*k-2]->Fill(NPV_all, weight);
	    h_nvtx_good[0][3][2*k-2]->Fill(NPV_good, weight);
	    h_njets[0][3][2*k-2]->Fill(njets, weight);
	    h_lepApt[0][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	    h_lepBpt[0][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeIn[0][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeIn[0][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeSplit[0][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeSplit[0][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeSplit[0][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeSplit[0][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeOut[0][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeOut[0][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepAeta[0][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	    h_lepBeta[0][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	    h_lepAphi[0][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	    h_lepBphi[0][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	    h_lepABpt[0][1][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepABeta[0][1][2*k-2]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    h_lepABeta_onebin[0][1][2*k-2]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeIn[0][1][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeSplit[0][1][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeSplit[0][1][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeOut[0][1][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	  }

	}

	for (int k = 2; k < 9; k++) {
	  
	  if (NPV_all < 5*k+10) { 
	    h_sidereal[0][3][2*k+5]->Fill(hour_sidereal, weight);
	    h_met[0][3][2*k+5]->Fill(met->Pt(), weight);
	    h_nvtx[0][3][2*k+5]->Fill(NPV_all, weight);
	    h_nvtx_good[0][3][2*k+5]->Fill(NPV_good, weight);
	    h_njets[0][3][2*k+5]->Fill(njets, weight);
	    h_lepApt[0][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	    h_lepBpt[0][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeIn[0][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeIn[0][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeSplit[0][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeSplit[0][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeSplit[0][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeSplit[0][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeOut[0][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeOut[0][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepAeta[0][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	    h_lepBeta[0][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	    h_lepAphi[0][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	    h_lepBphi[0][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	    h_lepABpt[0][1][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepABeta[0][1][2*k+5]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    h_lepABeta_onebin[0][1][2*k+5]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeIn[0][1][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeSplit[0][1][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeSplit[0][1][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeOut[0][1][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	  }
	  else { 
	    h_sidereal[0][3][2*k+6]->Fill(hour_sidereal, weight);
	    h_met[0][3][2*k+6]->Fill(met->Pt(), weight);
	    h_nvtx[0][3][2*k+6]->Fill(NPV_all, weight);
	    h_nvtx_good[0][3][2*k+6]->Fill(NPV_good, weight);
	    h_njets[0][3][2*k+6]->Fill(njets, weight);
	    h_lepApt[0][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	    h_lepBpt[0][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeIn[0][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeIn[0][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeSplit[0][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeSplit[0][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepApt_eeSplit[0][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepBpt_eeSplit[0][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepApt_eeOut[0][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepBpt_eeOut[0][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepAeta[0][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	    h_lepBeta[0][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	    h_lepAphi[0][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	    h_lepBphi[0][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	    h_lepABpt[0][1][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepABeta[0][1][2*k+6]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    h_lepABeta_onebin[0][1][2*k+6]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeIn[0][1][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeSplit[0][1][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) <= 1.479) h_lepABpt_eeSplit[0][1][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && fabs(lepSCEta->at(subleadinglepton)) > 1.479) h_lepABpt_eeOut[0][1][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	  }

	}

      }
    }

    else if (channelPdgIdProduct == (-11* 13)) { 

      // For the emu channel, lep A will be the electron and lep B the muon

      // cout << "Check \"HLTPathsOR_emu\" : " << HLTBitsManager_.passes_OR(HLTPathsOR_emu) << endl;
      // cout << "Check \"HLTPathsOR_met\" : " << HLTBitsManager_.passes_OR(HLTPathsOR_met) << endl; 
      h_sidereal[1][0][0]->Fill(hour_sidereal, weight);
      h_met[1][0][0]->Fill(met->Pt(), weight);
      h_nvtx[1][0][0]->Fill(NPV_all, weight);
      h_nvtx_good[1][0][0]->Fill(NPV_good, weight);
      h_njets[1][0][0]->Fill(njets, weight);
      if(abs(lepPdgId->at(leadinglepton)) == 11 && abs(lepPdgId->at(subleadinglepton)) == 13) {
	h_lepApt[1][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	h_lepAeta[1][0][0]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	h_lepAphi[1][0][0]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	h_lepBpt[1][0][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	h_lepBeta[1][0][0]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	h_lepBphi[1][0][0]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuOut[1][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuOut[1][0][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuIn[1][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuIn[1][0][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
      }
      else if(abs(lepPdgId->at(leadinglepton)) == 13 && abs(lepPdgId->at(subleadinglepton)) == 11) {
	h_lepApt[1][0][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);  
	h_lepAeta[1][0][0]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	h_lepAphi[1][0][0]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	h_lepBpt[1][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	h_lepBeta[1][0][0]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	h_lepBphi[1][0][0]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuOut[1][0][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuOut[1][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuIn[1][0][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuIn[1][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
      }

      if(HLTBitsManager_.passes_OR(HLTPathsOR_emu) == 1) {
	h_sidereal[1][1][0]->Fill(hour_sidereal, weight);
        h_met[1][1][0]->Fill(met->Pt(), weight);
        h_nvtx[1][1][0]->Fill(NPV_all, weight);
        h_nvtx_good[1][1][0]->Fill(NPV_good, weight);
        h_njets[1][1][0]->Fill(njets, weight);
	if(abs(lepPdgId->at(leadinglepton)) == 11 && abs(lepPdgId->at(subleadinglepton)) == 13) {
	  h_lepApt[1][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	  h_lepAeta[1][1][0]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	  h_lepAphi[1][1][0]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	  h_lepBpt[1][1][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	  h_lepBeta[1][1][0]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	  h_lepBphi[1][1][0]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	  if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuOut[1][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuOut[1][1][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuIn[1][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuIn[1][1][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	}
	else if(abs(lepPdgId->at(leadinglepton)) == 13 && abs(lepPdgId->at(subleadinglepton)) == 11) {
	  h_lepApt[1][1][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);  
	  h_lepAeta[1][1][0]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	  h_lepAphi[1][1][0]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	  h_lepBpt[1][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	  h_lepBeta[1][1][0]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	  h_lepBphi[1][1][0]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	  if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuOut[1][1][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuOut[1][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuIn[1][1][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuIn[1][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	}
      }
      if(HLTBitsManager_.passes_OR(HLTPathsOR_met) == 1) {
        nEventsPassingCrossTriggers_emu++;
	h_sidereal[1][2][0]->Fill(hour_sidereal, weight);
        h_met[1][2][0]->Fill(met->Pt(), weight);
        h_nvtx[1][2][0]->Fill(NPV_all, weight);
        h_nvtx_good[1][2][0]->Fill(NPV_good, weight);
        h_njets[1][2][0]->Fill(njets, weight);
	if(abs(lepPdgId->at(leadinglepton)) == 11 && abs(lepPdgId->at(subleadinglepton)) == 13) {
	  h_lepApt[1][2][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	  h_lepAeta[1][2][0]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	  h_lepAphi[1][2][0]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	  h_lepBpt[1][2][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuOut[1][2][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuOut[1][2][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuIn[1][2][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuIn[1][2][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	  h_lepBeta[1][2][0]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	  h_lepBphi[1][2][0]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	  h_lepABpt[1][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	  h_lepABeta[1][0][0]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	  h_lepABeta_onebin[1][0][0]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight); 
	  if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuOut[1][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuIn[1][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	}
	else if(abs(lepPdgId->at(leadinglepton)) == 13 && abs(lepPdgId->at(subleadinglepton)) == 11) {
	  h_lepApt[1][2][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);  
	  h_lepAeta[1][2][0]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	  h_lepAphi[1][2][0]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	  h_lepBpt[1][2][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	  h_lepBeta[1][2][0]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	  h_lepBphi[1][2][0]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	  if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuOut[1][2][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuOut[1][2][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuIn[1][2][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuIn[1][2][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	  h_lepABpt[1][0][0]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	  h_lepABeta[1][0][0]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	  h_lepABeta_onebin[1][0][0]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	  if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuOut[1][0][0]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuIn[1][0][0]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	}


	for (int k = 2; k < 6; k++) {
	  if (njets < k) { 
	    h_sidereal[1][2][2*k-3]->Fill(hour_sidereal, weight);
	    h_met[1][2][2*k-3]->Fill(met->Pt(), weight);
	    h_nvtx[1][2][2*k-3]->Fill(NPV_all, weight);
	    h_nvtx_good[1][2][2*k-3]->Fill(NPV_good, weight);
	    h_njets[1][2][2*k-3]->Fill(njets, weight);
	    if(abs(lepPdgId->at(leadinglepton)) == 11 && abs(lepPdgId->at(subleadinglepton)) == 13) {
	      h_lepApt[1][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	      h_lepAeta[1][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	      h_lepAphi[1][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	      h_lepBpt[1][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepBeta[1][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	      h_lepBphi[1][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuOut[1][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuOut[1][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuIn[1][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuIn[1][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepABpt[1][0][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepABeta[1][0][2*k-3]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	      h_lepABeta_onebin[1][0][2*k-3]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuOut[1][0][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuIn[1][0][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    }
	    else if(abs(lepPdgId->at(leadinglepton)) == 13 && abs(lepPdgId->at(subleadinglepton)) == 11) {
	      h_lepApt[1][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);  
	      h_lepAeta[1][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	      h_lepAphi[1][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	      h_lepBpt[1][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepBeta[1][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	      h_lepBphi[1][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuOut[1][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuOut[1][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuIn[1][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuIn[1][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepABpt[1][0][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepABeta[1][0][2*k-3]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	      h_lepABeta_onebin[1][0][2*k-3]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuOut[1][0][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuIn[1][0][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	    }
	  }
	  else { 
	    h_sidereal[1][2][2*k-2]->Fill(hour_sidereal, weight);
	    h_met[1][2][2*k-2]->Fill(met->Pt(), weight);
	    h_nvtx[1][2][2*k-2]->Fill(NPV_all, weight);
	    h_nvtx_good[1][2][2*k-2]->Fill(NPV_good, weight);
	    h_njets[1][2][2*k-2]->Fill(njets, weight);
	    if(abs(lepPdgId->at(leadinglepton)) == 11 && abs(lepPdgId->at(subleadinglepton)) == 13) {
	      h_lepApt[1][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	      h_lepAeta[1][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	      h_lepAphi[1][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	      h_lepBpt[1][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepBeta[1][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	      h_lepBphi[1][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuOut[1][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuOut[1][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuIn[1][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuIn[1][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepABpt[1][0][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepABeta[1][0][2*k-2]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	      h_lepABeta_onebin[1][0][2*k-2]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuOut[1][0][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuIn[1][0][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    }
	    else if(abs(lepPdgId->at(leadinglepton)) == 13 && abs(lepPdgId->at(subleadinglepton)) == 11) {
	      h_lepApt[1][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);  
	      h_lepAeta[1][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	      h_lepAphi[1][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	      h_lepBpt[1][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepBeta[1][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	      h_lepBphi[1][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuOut[1][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuOut[1][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuIn[1][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuIn[1][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepABpt[1][0][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepABeta[1][0][2*k-2]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	      h_lepABeta_onebin[1][0][2*k-2]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuOut[1][0][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuIn[1][0][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	    }
	  }

	}

	for (int k = 2; k < 9; k++) {

	  if (NPV_all < 5*k+10) { 
	    h_sidereal[1][2][2*k+5]->Fill(hour_sidereal, weight);
	    h_met[1][2][2*k+5]->Fill(met->Pt(), weight);
	    h_nvtx[1][2][2*k+5]->Fill(NPV_all, weight);
	    h_nvtx_good[1][2][2*k+5]->Fill(NPV_good, weight);
	    h_njets[1][2][2*k+5]->Fill(njets, weight);
	    if(abs(lepPdgId->at(leadinglepton)) == 11 && abs(lepPdgId->at(subleadinglepton)) == 13) {
	      h_lepApt[1][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	      h_lepAeta[1][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	      h_lepAphi[1][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	      h_lepBpt[1][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepBeta[1][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	      h_lepBphi[1][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuOut[1][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuOut[1][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuIn[1][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuIn[1][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepABpt[1][0][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepABeta[1][0][2*k+5]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	      h_lepABeta_onebin[1][0][2*k+5]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuOut[1][0][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuIn[1][0][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    }
	    else if(abs(lepPdgId->at(leadinglepton)) == 13 && abs(lepPdgId->at(subleadinglepton)) == 11) {
	      h_lepApt[1][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);  
	      h_lepAeta[1][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	      h_lepAphi[1][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	      h_lepBpt[1][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepBeta[1][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	      h_lepBphi[1][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuOut[1][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuOut[1][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuIn[1][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuIn[1][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepABpt[1][0][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepABeta[1][0][2*k+5]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	      h_lepABeta_onebin[1][0][2*k+5]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuOut[1][0][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuIn[1][0][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	    }
	  }
	  else { 
	    h_sidereal[1][2][2*k+6]->Fill(hour_sidereal, weight);
	    h_met[1][2][2*k+6]->Fill(met->Pt(), weight);
	    h_nvtx[1][2][2*k+6]->Fill(NPV_all, weight);
	    h_nvtx_good[1][2][2*k+6]->Fill(NPV_good, weight);
	    h_njets[1][2][2*k+6]->Fill(njets, weight);
	    if(abs(lepPdgId->at(leadinglepton)) == 11 && abs(lepPdgId->at(subleadinglepton)) == 13) {
	      h_lepApt[1][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	      h_lepAeta[1][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	      h_lepAphi[1][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	      h_lepBpt[1][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepBeta[1][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	      h_lepBphi[1][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuOut[1][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuOut[1][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuIn[1][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuIn[1][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepABpt[1][0][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepABeta[1][0][2*k+6]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	      h_lepABeta_onebin[1][0][2*k+6]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuOut[1][0][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuIn[1][0][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    }
	    else if(abs(lepPdgId->at(leadinglepton)) == 13 && abs(lepPdgId->at(subleadinglepton)) == 11) {
	      h_lepApt[1][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);  
	      h_lepAeta[1][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	      h_lepAphi[1][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	      h_lepBpt[1][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepBeta[1][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	      h_lepBphi[1][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuOut[1][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuOut[1][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuIn[1][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuIn[1][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepABpt[1][0][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepABeta[1][0][2*k+6]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	      h_lepABeta_onebin[1][0][2*k+6]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuOut[1][0][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuIn[1][0][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	    }
	  }
	}
      }
      if(HLTBitsManager_.passes_OR(HLTPathsOR_emu) == 1 && HLTBitsManager_.passes_OR(HLTPathsOR_met) == 1) {
        nEventsPassingCrossTriggersAndDilepTriggers_emu++;
	h_sidereal[1][3][0]->Fill(hour_sidereal, weight);
        h_met[1][3][0]->Fill(met->Pt(), weight);
        h_nvtx[1][3][0]->Fill(NPV_all, weight);
        h_nvtx_good[1][3][0]->Fill(NPV_good, weight);
        h_njets[1][3][0]->Fill(njets, weight);
	if(abs(lepPdgId->at(leadinglepton)) == 11 && abs(lepPdgId->at(subleadinglepton)) == 13) {
	  h_lepApt[1][3][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	  h_lepAeta[1][3][0]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	  h_lepAphi[1][3][0]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	  h_lepBpt[1][3][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	  h_lepBeta[1][3][0]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	  h_lepBphi[1][3][0]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	  if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuOut[1][3][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuOut[1][3][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuIn[1][3][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuIn[1][3][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	  h_lepABpt[1][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	  h_lepABeta[1][1][0]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	  h_lepABeta_onebin[1][1][0]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	  if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuOut[1][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuIn[1][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	}
	else if(abs(lepPdgId->at(leadinglepton)) == 13 && abs(lepPdgId->at(subleadinglepton)) == 11) {
	  h_lepApt[1][3][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);  
	  h_lepAeta[1][3][0]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	  h_lepAphi[1][3][0]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	  h_lepBpt[1][3][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	  h_lepBeta[1][3][0]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	  h_lepBphi[1][3][0]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	  if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuOut[1][3][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuOut[1][3][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuIn[1][3][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuIn[1][3][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	  h_lepABpt[1][1][0]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	  h_lepABeta[1][1][0]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	  h_lepABeta_onebin[1][1][0]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	  if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuOut[1][1][0]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	  if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuIn[1][1][0]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	}


	for (int k = 2; k < 6; k++) {

	  if (njets < k) { 
	    h_sidereal[1][3][2*k-3]->Fill(hour_sidereal, weight);
	    h_met[1][3][2*k-3]->Fill(met->Pt(), weight);
	    h_nvtx[1][3][2*k-3]->Fill(NPV_all, weight);
	    h_nvtx_good[1][3][2*k-3]->Fill(NPV_good, weight);
	    h_njets[1][3][2*k-3]->Fill(njets, weight);
	    if(abs(lepPdgId->at(leadinglepton)) == 11 && abs(lepPdgId->at(subleadinglepton)) == 13) {
	      h_lepApt[1][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	      h_lepAeta[1][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	      h_lepAphi[1][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	      h_lepBpt[1][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepBeta[1][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	      h_lepBphi[1][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuOut[1][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuOut[1][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuIn[1][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuIn[1][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepABpt[1][1][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepABeta[1][1][2*k-3]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	      h_lepABeta_onebin[1][1][2*k-3]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuOut[1][1][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuIn[1][1][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    }
	    else if(abs(lepPdgId->at(leadinglepton)) == 13 && abs(lepPdgId->at(subleadinglepton)) == 11) {
	      h_lepApt[1][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);  
	      h_lepAeta[1][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	      h_lepAphi[1][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	      h_lepBpt[1][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepBeta[1][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	      h_lepBphi[1][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuOut[1][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuOut[1][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuIn[1][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuIn[1][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepABpt[1][1][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepABeta[1][1][2*k-3]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	      h_lepABeta_onebin[1][1][2*k-3]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuOut[1][1][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuIn[1][1][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	    }
	  }
	  else { 
	    h_sidereal[1][3][2*k-2]->Fill(hour_sidereal, weight);
	    h_met[1][3][2*k-2]->Fill(met->Pt(), weight);
	    h_nvtx[1][3][2*k-2]->Fill(NPV_all, weight);
	    h_nvtx_good[1][3][2*k-2]->Fill(NPV_good, weight);
	    h_njets[1][3][2*k-2]->Fill(njets, weight);
	    if(abs(lepPdgId->at(leadinglepton)) == 11 && abs(lepPdgId->at(subleadinglepton)) == 13) {
	      h_lepApt[1][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	      h_lepAeta[1][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	      h_lepAphi[1][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	      h_lepBpt[1][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepBeta[1][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	      h_lepBphi[1][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuOut[1][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuOut[1][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuIn[1][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuIn[1][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepABpt[1][1][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepABeta[1][1][2*k-2]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	      h_lepABeta_onebin[1][1][2*k-2]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuOut[1][1][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuIn[1][1][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    }
	    else if(abs(lepPdgId->at(leadinglepton)) == 13 && abs(lepPdgId->at(subleadinglepton)) == 11) {
	      h_lepApt[1][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);  
	      h_lepAeta[1][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	      h_lepAphi[1][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	      h_lepBpt[1][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepBeta[1][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	      h_lepBphi[1][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuOut[1][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuOut[1][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuIn[1][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuIn[1][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepABpt[1][1][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepABeta[1][1][2*k-2]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	      h_lepABeta_onebin[1][1][2*k-2]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuOut[1][1][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuIn[1][1][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	    }
	  }
	}

	for (int k = 2; k < 9; k++) {

	  if (NPV_all < 5*k+10) { 
	    h_sidereal[1][3][2*k+5]->Fill(hour_sidereal, weight);
	    h_met[1][3][2*k+5]->Fill(met->Pt(), weight);
	    h_nvtx[1][3][2*k+5]->Fill(NPV_all, weight);
	    h_nvtx_good[1][3][2*k+5]->Fill(NPV_good, weight);
	    h_njets[1][3][2*k+5]->Fill(njets, weight);
	    if(abs(lepPdgId->at(leadinglepton)) == 11 && abs(lepPdgId->at(subleadinglepton)) == 13) {
	      h_lepApt[1][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	      h_lepAeta[1][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	      h_lepAphi[1][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	      h_lepBpt[1][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepBeta[1][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	      h_lepBphi[1][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuOut[1][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuOut[1][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuIn[1][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuIn[1][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepABpt[1][1][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepABeta[1][1][2*k+5]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	      h_lepABeta_onebin[1][1][2*k+5]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuOut[1][1][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuIn[1][1][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    }
	    else if(abs(lepPdgId->at(leadinglepton)) == 13 && abs(lepPdgId->at(subleadinglepton)) == 11) {
	      h_lepApt[1][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);  
	      h_lepAeta[1][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	      h_lepAphi[1][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	      h_lepBpt[1][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepBeta[1][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	      h_lepBphi[1][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuOut[1][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuOut[1][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuIn[1][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuIn[1][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepABpt[1][1][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepABeta[1][1][2*k+5]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	      h_lepABeta_onebin[1][1][2*k+5]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuOut[1][1][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuIn[1][1][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	    }
	  }
	  else { 
	    h_sidereal[1][3][2*k+6]->Fill(hour_sidereal, weight);
	    h_met[1][3][2*k+6]->Fill(met->Pt(), weight);
	    h_nvtx[1][3][2*k+6]->Fill(NPV_all, weight);
	    h_nvtx_good[1][3][2*k+6]->Fill(NPV_good, weight);
	    h_njets[1][3][2*k+6]->Fill(njets, weight);
	    if(abs(lepPdgId->at(leadinglepton)) == 11 && abs(lepPdgId->at(subleadinglepton)) == 13) {
	      h_lepApt[1][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	      h_lepAeta[1][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	      h_lepAphi[1][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	      h_lepBpt[1][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepBeta[1][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	      h_lepBphi[1][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuOut[1][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuOut[1][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepApt_emuIn[1][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepBpt_emuIn[1][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepABpt[1][1][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	      h_lepABeta[1][1][2*k+6]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	      h_lepABeta_onebin[1][1][2*k+6]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) > 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuOut[1][1][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(leadinglepton)) <= 1.479 && abs(lepPdgId->at(leadinglepton)) == 11) h_lepABpt_emuIn[1][1][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    }
	    else if(abs(lepPdgId->at(leadinglepton)) == 13 && abs(lepPdgId->at(subleadinglepton)) == 11) {
	      h_lepApt[1][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);  
	      h_lepAeta[1][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	      h_lepAphi[1][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	      h_lepBpt[1][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepBeta[1][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	      h_lepBphi[1][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuOut[1][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuOut[1][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepApt_emuIn[1][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepBpt_emuIn[1][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepABpt[1][1][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	      h_lepABeta[1][1][2*k+6]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	      h_lepABeta_onebin[1][1][2*k+6]->Fill(abs(v_leptons->at(subleadinglepton).Eta()), abs(v_leptons->at(leadinglepton).Eta()), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) > 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuOut[1][1][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	      if(fabs(lepSCEta->at(subleadinglepton)) <= 1.479 && abs(lepPdgId->at(subleadinglepton)) == 11) h_lepABpt_emuIn[1][1][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), v_leptons->at(leadinglepton).Pt(), weight);
	    }
	  }

	}

      }
    }

    else if (channelPdgIdProduct == (-13* 13)) { 
      // cout << "Check \"HLTPathsOR_mumu\" : " << HLTBitsManager_.passes_OR(HLTPathsOR_mumu) << endl; 
      // cout << "Check \"HLTPathsOR_met\" : " << HLTBitsManager_.passes_OR(HLTPathsOR_met) << endl; 
      h_sidereal[2][0][0]->Fill(hour_sidereal, weight);
      h_met[2][0][0]->Fill(met->Pt(), weight);
      h_nvtx[2][0][0]->Fill(NPV_all, weight);
      h_nvtx_good[2][0][0]->Fill(NPV_good, weight);
      h_njets[2][0][0]->Fill(njets, weight);
      h_lepApt[2][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
      h_lepBpt[2][0][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
      h_lepAeta[2][0][0]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
      h_lepBeta[2][0][0]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
      h_lepAphi[2][0][0]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
      h_lepBphi[2][0][0]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);

      if(HLTBitsManager_.passes_OR(HLTPathsOR_mumu) == 1) {
	h_sidereal[2][1][0]->Fill(hour_sidereal, weight);
        h_met[2][1][0]->Fill(met->Pt(), weight);
        h_nvtx[2][1][0]->Fill(NPV_all, weight);
        h_nvtx_good[2][1][0]->Fill(NPV_good, weight);
        h_njets[2][1][0]->Fill(njets, weight);
        h_lepApt[2][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
        h_lepBpt[2][1][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
        h_lepAeta[2][1][0]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
        h_lepBeta[2][1][0]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
        h_lepAphi[2][1][0]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
        h_lepBphi[2][1][0]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
      }
      if(HLTBitsManager_.passes_OR(HLTPathsOR_met) == 1) {
        nEventsPassingCrossTriggers_mumu++;
	h_sidereal[2][2][0]->Fill(hour_sidereal, weight);
        h_met[2][2][0]->Fill(met->Pt(), weight);
        h_nvtx[2][2][0]->Fill(NPV_all, weight);
        h_nvtx_good[2][2][0]->Fill(NPV_good, weight);
        h_njets[2][2][0]->Fill(njets, weight);
        h_lepApt[2][2][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
        h_lepBpt[2][2][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
        h_lepAeta[2][2][0]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
        h_lepBeta[2][2][0]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
        h_lepAphi[2][2][0]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
        h_lepBphi[2][2][0]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
        h_lepABpt[2][0][0]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
        h_lepABeta[2][0][0]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
        h_lepABeta_onebin[2][0][0]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);

	for (int k = 2; k < 6; k++) {

	  if (njets < k) { 
	    h_sidereal[2][2][2*k-3]->Fill(hour_sidereal, weight);
	    h_met[2][2][2*k-3]->Fill(met->Pt(), weight);
	    h_nvtx[2][2][2*k-3]->Fill(NPV_all, weight);
	    h_nvtx_good[2][2][2*k-3]->Fill(NPV_good, weight);
	    h_njets[2][2][2*k-3]->Fill(njets, weight);
	    h_lepApt[2][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	    h_lepBpt[2][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepAeta[2][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	    h_lepBeta[2][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	    h_lepAphi[2][2][2*k-3]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	    h_lepBphi[2][2][2*k-3]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	    h_lepABpt[2][0][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepABeta[2][0][2*k-3]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    h_lepABeta_onebin[2][0][2*k-3]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	  }
	  else { 
	    h_sidereal[2][2][2*k-2]->Fill(hour_sidereal, weight);
	    h_met[2][2][2*k-2]->Fill(met->Pt(), weight);
	    h_nvtx[2][2][2*k-2]->Fill(NPV_all, weight);
	    h_nvtx_good[2][2][2*k-2]->Fill(NPV_good, weight);
	    h_njets[2][2][2*k-2]->Fill(njets, weight);
	    h_lepApt[2][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	    h_lepBpt[2][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepAeta[2][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	    h_lepBeta[2][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	    h_lepAphi[2][2][2*k-2]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	    h_lepBphi[2][2][2*k-2]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	    h_lepABpt[2][0][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepABeta[2][0][2*k-2]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    h_lepABeta_onebin[2][0][2*k-2]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	  }
	}

	for (int k = 2; k < 9; k++) {
	  if (NPV_all < 5*k+10) { 
	    h_sidereal[2][2][2*k+5]->Fill(hour_sidereal, weight);
	    h_met[2][2][2*k+5]->Fill(met->Pt(), weight);
	    h_nvtx[2][2][2*k+5]->Fill(NPV_all, weight);
	    h_nvtx_good[2][2][2*k+5]->Fill(NPV_good, weight);
	    h_njets[2][2][2*k+5]->Fill(njets, weight);
	    h_lepApt[2][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	    h_lepBpt[2][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepAeta[2][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	    h_lepBeta[2][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	    h_lepAphi[2][2][2*k+5]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	    h_lepBphi[2][2][2*k+5]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	    h_lepABpt[2][0][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepABeta[2][0][2*k+5]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    h_lepABeta_onebin[2][0][2*k+5]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	  }
	  else { 
	    h_sidereal[2][2][2*k+6]->Fill(hour_sidereal, weight);
	    h_met[2][2][2*k+6]->Fill(met->Pt(), weight);
	    h_nvtx[2][2][2*k+6]->Fill(NPV_all, weight);
	    h_nvtx_good[2][2][2*k+6]->Fill(NPV_good, weight);
	    h_njets[2][2][2*k+6]->Fill(njets, weight);
	    h_lepApt[2][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	    h_lepBpt[2][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepAeta[2][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	    h_lepBeta[2][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	    h_lepAphi[2][2][2*k+6]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	    h_lepBphi[2][2][2*k+6]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	    h_lepABpt[2][0][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepABeta[2][0][2*k+6]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    h_lepABeta_onebin[2][0][2*k+6]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	  }
	}
      }
      if(HLTBitsManager_.passes_OR(HLTPathsOR_mumu) == 1 && HLTBitsManager_.passes_OR(HLTPathsOR_met) == 1) {
        nEventsPassingCrossTriggersAndDilepTriggers_mumu++;
	h_sidereal[2][3][0]->Fill(hour_sidereal, weight);
        h_met[2][3][0]->Fill(met->Pt(), weight);
        h_nvtx[2][3][0]->Fill(NPV_all, weight);
        h_nvtx_good[2][3][0]->Fill(NPV_good, weight);
        h_njets[2][3][0]->Fill(njets, weight);
        h_lepApt[2][3][0]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
        h_lepBpt[2][3][0]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
        h_lepAeta[2][3][0]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
        h_lepBeta[2][3][0]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
        h_lepAphi[2][3][0]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
        h_lepBphi[2][3][0]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
        h_lepABpt[2][1][0]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
        h_lepABeta[2][1][0]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
        h_lepABeta_onebin[2][1][0]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);

	for (int k = 2; k < 6; k++) {

	  if (njets < k) { 
	    h_sidereal[2][3][2*k-3]->Fill(hour_sidereal, weight);
	    h_met[2][3][2*k-3]->Fill(met->Pt(), weight);
	    h_nvtx[2][3][2*k-3]->Fill(NPV_all, weight);
	    h_nvtx_good[2][3][2*k-3]->Fill(NPV_good, weight);
	    h_njets[2][3][2*k-3]->Fill(njets, weight);
	    h_lepApt[2][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	    h_lepBpt[2][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepAeta[2][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	    h_lepBeta[2][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	    h_lepAphi[2][3][2*k-3]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	    h_lepBphi[2][3][2*k-3]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	    h_lepABpt[2][1][2*k-3]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepABeta[2][1][2*k-3]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    h_lepABeta_onebin[2][1][2*k-3]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	  }
	  else { 
	    h_sidereal[2][3][2*k-2]->Fill(hour_sidereal, weight);
	    h_met[2][3][2*k-2]->Fill(met->Pt(), weight);
	    h_nvtx[2][3][2*k-2]->Fill(NPV_all, weight);
	    h_nvtx_good[2][3][2*k-2]->Fill(NPV_good, weight);
	    h_njets[2][3][2*k-2]->Fill(njets, weight);
	    h_lepApt[2][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	    h_lepBpt[2][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepAeta[2][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	    h_lepBeta[2][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	    h_lepAphi[2][3][2*k-2]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	    h_lepBphi[2][3][2*k-2]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	    h_lepABpt[2][1][2*k-2]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepABeta[2][1][2*k-2]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    h_lepABeta_onebin[2][1][2*k-2]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	  }
	}

	for (int k = 2; k < 9; k++) {

	  if (NPV_all < 5*k+10) { 
	    h_sidereal[2][3][2*k+5]->Fill(hour_sidereal, weight);
	    h_met[2][3][2*k+5]->Fill(met->Pt(), weight);
	    h_nvtx[2][3][2*k+5]->Fill(NPV_all, weight);
	    h_nvtx_good[2][3][2*k+5]->Fill(NPV_good, weight);
	    h_njets[2][3][2*k+5]->Fill(njets, weight);
	    h_lepApt[2][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	    h_lepBpt[2][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepAeta[2][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	    h_lepBeta[2][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	    h_lepAphi[2][3][2*k+5]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	    h_lepBphi[2][3][2*k+5]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	    h_lepABpt[2][1][2*k+5]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepABeta[2][1][2*k+5]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    h_lepABeta_onebin[2][1][2*k+5]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	  }
	  else { 
	    h_sidereal[2][3][2*k+6]->Fill(hour_sidereal, weight);
	    h_met[2][3][2*k+6]->Fill(met->Pt(), weight);
	    h_nvtx[2][3][2*k+6]->Fill(NPV_all, weight);
	    h_nvtx_good[2][3][2*k+6]->Fill(NPV_good, weight);
	    h_njets[2][3][2*k+6]->Fill(njets, weight);
	    h_lepApt[2][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), weight);  
	    h_lepBpt[2][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepAeta[2][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Eta(), weight);
	    h_lepBeta[2][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Eta(), weight);
	    h_lepAphi[2][3][2*k+6]->Fill(v_leptons->at(leadinglepton).Phi(), weight);
	    h_lepBphi[2][3][2*k+6]->Fill(v_leptons->at(subleadinglepton).Phi(), weight);
	    h_lepABpt[2][1][2*k+6]->Fill(v_leptons->at(leadinglepton).Pt(), v_leptons->at(subleadinglepton).Pt(), weight);
	    h_lepABeta[2][1][2*k+6]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	    h_lepABeta_onebin[2][1][2*k+6]->Fill(abs(v_leptons->at(leadinglepton).Eta()), abs(v_leptons->at(subleadinglepton).Eta()), weight);
	  }

	}
      }
    }

    else { continue; }
    
  }
  
  chain->ResetBranchAddresses();

  cout << "nEvents = " << nEvents << endl;

  cout << "nEventsPassingCrossTriggers_ee = " << nEventsPassingCrossTriggers_ee << endl;
  cout << "nEventsPassingCrossTriggersAndDilepTriggers_ee = " << nEventsPassingCrossTriggersAndDilepTriggers_ee << endl;
  cout << "Efficiency_ee = " << nEventsPassingCrossTriggersAndDilepTriggers_ee/nEventsPassingCrossTriggers_ee << endl;

  cout << "nEventsPassingCrossTriggers_emu = " << nEventsPassingCrossTriggers_emu << endl;
  cout << "nEventsPassingCrossTriggersAndDilepTriggers_emu = " << nEventsPassingCrossTriggersAndDilepTriggers_emu << endl;
  cout << "Efficiency_emu = " << nEventsPassingCrossTriggersAndDilepTriggers_emu/nEventsPassingCrossTriggers_emu << endl;

  cout << "nEventsPassingCrossTriggers_mumu = " << nEventsPassingCrossTriggers_mumu << endl;
  cout << "nEventsPassingCrossTriggersAndDilepTriggers_mumu = " << nEventsPassingCrossTriggersAndDilepTriggers_mumu << endl;
  cout << "Efficiency_mumu = " << nEventsPassingCrossTriggersAndDilepTriggers_mumu/nEventsPassingCrossTriggers_mumu << endl;
  

 
  TFile *output_hists = new TFile(("triggereff/output_hists"+suffix+"_"+era_+".root").c_str(),"RECREATE"); 

  weightedEvents->Write();
  unweightedEvents->Write();
  h_LepSF->Write();
  h_PUSF->Write();

  for (unsigned short i(0); i < 3; i++) {
    for (unsigned short j(0); j < 4; j++) {
      for (unsigned short k(0); k < 23; k++) {

	if(MCorData_ == "MC") {
	  for (unsigned short l(1); l < 25; l++) {
	    //              cout << "l: " <<  l << endl;
	    //              cout <<  h_lepABeta_onebin[i][j][k]->GetBinContent(1,1)/4 << endl;
	    h_sidereal[i][j][k]->SetBinContent(l,h_lepABeta_onebin[i][j][k]->GetBinContent(1,1)/24);
	    h_sidereal[i][j][k]->SetBinError(l,0);
	  }
	}

	h_sidereal[i][j][k]->Write();
	h_met[i][j][k]->Write();
	h_nvtx[i][j][k]->Write();
	h_nvtx_good[i][j][k]->Write();
	h_njets[i][j][k]->Write();
	h_lepApt[i][j][k]->Write();
	h_lepAeta[i][j][k]->Write();
	h_lepAphi[i][j][k]->Write();
	h_lepBpt[i][j][k]->Write();
	h_lepBeta[i][j][k]->Write();
	h_lepBphi[i][j][k]->Write();

	h_lepApt_eeIn[i][j][k]->Write();
	h_lepBpt_eeIn[i][j][k]->Write();

	h_lepApt_eeOut[i][j][k]->Write();
	h_lepBpt_eeOut[i][j][k]->Write();

	h_lepApt_eeSplit[i][j][k]->Write();
	h_lepBpt_eeSplit[i][j][k]->Write();

	h_lepApt_emuIn[i][j][k]->Write();
	h_lepBpt_emuIn[i][j][k]->Write();

	h_lepApt_emuOut[i][j][k]->Write();
	h_lepBpt_emuOut[i][j][k]->Write();

	if (j>1){
	  h_lepABeta_onebin[i][j-2][k]->Write();
	  h_lepABeta[i][j-2][k]->Write();
	  h_lepABpt[i][j-2][k]->Write();

	  h_lepABpt_eeIn[i][j-2][k]->Write();
	  h_lepABpt_eeSplit[i][j-2][k]->Write();
	  h_lepABpt_eeOut[i][j-2][k]->Write();
	  h_lepABpt_emuIn[i][j-2][k]->Write();
	  h_lepABpt_emuOut[i][j-2][k]->Write();
       
	}
      }
    }
  }
  
  output_hists->Write();
  output_hists->Close();

}

int main(int argc, char** argv) {

  CLParameter<std::string> opt_y("y", "Era", false, 1, 1);
  CLParameter<std::string> opt_t("t", "MC or Data", false, 1, 1);
  CLParameter<std::string> opt_s("s", "Sample", false, 1, 1);
  CLParameter<std::string> opt_n("n", "Sample Number", false, 0, 1);
  CLAnalyser::interpretGlobal(argc, argv);

  if (opt_n.isSet()) {
    triggereff(opt_y[0], opt_t[0], opt_s[0], opt_n[0]);
  }
  else {
    triggereff(opt_y[0], opt_t[0], opt_s[0], "");
  }
}
