#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;
void MkCondorScriptsNominalVariations(void){
  gSystem->mkdir("condor_jobs");
  gSystem->mkdir("condor_jobs/logs");
  ofstream BigJob("RunNominalVariations.sh");
  
  //string mode="";
  string mode="-m spinCorr";
  string chan[3]={"ee","emu","mumu"}; 
  string sample[9]={"zz", "qcd", "tw.root", "ttbarbg.root", "wtol", "ww", "wz", "ttbarW", "ttgjets"};
  /*string syst[80]={
		   //"JES_UP", "JES_DOWN",
		   //"JER_UP", "JER_DOWN", 
		   //"JESAbsoluteStat_UP", "JESAbsoluteStat_DOWN", 
		   //"JESAbsoluteScale_UP", "JESAbsoluteScale_DOWN", 
		   //"JESAbsoluteFlavMap_UP", "JESAbsoluteFlavMap_DOWN", 
		   //"JESAbsoluteMPFBias_UP", "JESAbsoluteMPFBias_DOWN", 
		   //"JESFragmentation_UP", "JESFragmentation_DOWN", 
		   //"JESSinglePionECAL_UP", "JESSinglePionECAL_DOWN", 
		   //"JESSinglePionHCAL_UP", "JESSinglePionHCAL_DOWN", 
		   //"JESFlavorQCD_UP", "JESFlavorQCD_DOWN", 
		   //"JESTimePtEta_UP", "JESTimePtEta_DOWN", 
		   //"JESRelativeJEREC1_UP", "JESRelativeJEREC1_DOWN", 
		   //"JESRelativeJEREC2_UP", "JESRelativeJEREC2_DOWN", 
		   //"JESRelativeJERHF_UP", "JESRelativeJERHF_DOWN", 
		   //"JESRelativePtBB_UP", "JESRelativePtBB_DOWN", 
		   //"JESRelativePtEC1_UP", "JESRelativePtEC1_DOWN", 
		   //"JESRelativePtEC2_UP", "JESRelativePtEC2_DOWN", 
		   //"JESRelativePtHF_UP", "JESRelativePtHF_DOWN", 
		   //"JESRelativeFSR_UP", "JESRelativeFSR_DOWN", 
		   //"JESRelativeStatFSR_UP", "JESRelativeStatFSR_DOWN", 
		   //"JESRelativeStatEC_UP", "JESRelativeStatEC_DOWN", 
		   "JESRelativeStatHF_UP", "JESRelativeStatHF_DOWN", 
		   "JESPileUpDataMC_UP", "JESPileUpDataMC_DOWN", 
		   "JESPileUpPtRef_UP", "JESPileUpPtRef_DOWN", 
		   "JESPileUpPtEC1_UP", "JESPileUpPtEC1_DOWN", 
		   "JESPileUpPtEC2_UP", "JESPileUpPtEC2_DOWN", 
		   "JESPileUpPtHF_UP", "JESPileUpPtHF_DOWN", 
		   "JESPileUpPtBB_UP", "JESPileUpPtBB_DOWN", 
		   "PU_UP", "PU_DOWN",
		   "TRIG_UP", "TRIG_DOWN", 
		   "TRIG_ETA_UP", "TRIG_ETA_DOWN", 
		   "LEPT_UP", "LEPT_DOWN", 
		   "UNCLUSTERED_UP", "UNCLUSTERED_DOWN", 
		   "KIN_UP", "KIN_DOWN", 
		   "BTAG_UP", "BTAG_DOWN", 
		   "BTAG_LJET_UP", "BTAG_LJET_DOWN", 
		   "BTAG_PT_UP", "BTAG_PT_DOWN", 
		   "BTAG_ETA_UP", "BTAG_ETA_DOWN", 
		   "BTAG_LJET_PT_UP" //, "BTAG_LJET_PT_DOWN", 
		   //"BTAG_LJET_ETA_UP", "BTAG_LJET_ETA_DOWN"
                   };*/

  string syst[11]={
		   "KIN_UP", "KIN_DOWN", 
		   "BTAG_UP", "BTAG_DOWN", 
		   "BTAG_LJET_UP", "BTAG_LJET_DOWN", 
		   "BTAG_PT_UP", "BTAG_PT_DOWN", 
		   "BTAG_ETA_UP", "BTAG_ETA_DOWN", 
		   "BTAG_LJET_PT_UP"}; 

  
  vector<string>sampleVar_arg;
  vector<string>output_arg;
  vector<string>tree_arg;

  for(int j=0; j<11; j++){
    //for(int i=0; i<3; i++){
    for(int i=1; i<3; i++){ //just run over ee for now
      sampleVar_arg.push_back("-f ttbarsignalplustau.root "+mode+" -c "+chan[i]+" -s "+syst[j]);
      sampleVar_arg.push_back("-f ttbarsignalplustau.root "+mode+" -c "+chan[i]+" --bgviatau -s "+syst[j]);
      sampleVar_arg.push_back("-f dy -d 11 "+mode+" -c "+chan[i]+" -s "+syst[j]);
      sampleVar_arg.push_back("-f dy -d 13 "+mode+" -c "+chan[i]+" -s "+syst[j]);
      sampleVar_arg.push_back("-f dy -d 15 "+mode+" -c "+chan[i]+" -s "+syst[j]);
      sampleVar_arg.push_back("-f ttbarZ "+mode+" -c "+chan[i]+" -s "+syst[j]);
    
      for(int k=0; k<6; k++){
	output_arg.push_back("selectionRoot/"+syst[j]+"/"+chan[i]);
	tree_arg.push_back("spinCorrInput/"+syst[j]+"/"+chan[i]);
      }
      
      for(int k=0; k<9; k++){
	sampleVar_arg.push_back("-f "+sample[k]+" "+mode+" -c "+chan[i]+" -s "+syst[j]);
	output_arg.push_back("selectionRoot/"+syst[j]+"/"+chan[i]);
	tree_arg.push_back("spinCorrInput/"+syst[j]+"/"+chan[i]);
      }
    }
  }
  
  BigJob << "#!/bin/bash" <<endl;
  BigJob << "voms-proxy-init -voms cms -valid 100:00:00" <<endl;
  //BigJob << "tar -X excludefromtar.txt -zcvf CMSSW.tgz /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/" <<endl;
  for (int i=0; i<sampleVar_arg.size(); i++){

    char ShellName[600];
    sprintf(ShellName, "condor_jobs/condor_doit_%i.sh",i);
    
    char JDLName[600];
    sprintf(JDLName, "condor_jobs/AnaCondor_%i.jdl",i);
    
    ofstream jdl(JDLName);
    ofstream shel(ShellName);
    BigJob << "condor_submit " << JDLName << endl;
    
    //Make the JDL file:
    
    jdl << "Universe                = vanilla" << endl;
    jdl << "Environment             = CONDORJOBID=$(Process)" << endl;
    jdl << "notification            = Error" << endl;
    jdl << "X509UserProxy = /tmp/x509up_u587091" << endl;
    jdl << "should_transfer_files   = YES" << endl;
    jdl << "Transfer_input_files = condor_jobs/condor_doit_"<<i<<".sh" << endl;
    jdl << "Executable = /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<ShellName << endl;
    jdl << "when_to_transfer_output = ON_EXIT" << endl;
    jdl << "Output = /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/condor_jobs/logs/batch_"<<i<<"_$(cluster)_$(process).stdout" << endl;
    jdl << "Error  = /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/condor_jobs/logs/batch_"<<i<<"_$(cluster)_$(process).stderr" << endl;
    jdl << "Log    = /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/condor_jobs/logs/batch_"<<i<<"_$(cluster)_$(process).condor" << endl;
    jdl << "Queue 1" <<endl; 
    
    //Make the Shell script: 
    shel << "#!/bin/sh" << endl;
    shel << "source /cvmfs/cms.cern.ch/cmsset_default.sh" << endl;
    shel << "cp /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/CMSSW.tgz ." <<endl;
    shel << "tar -zxf CMSSW.tgz" <<endl;
    shel << "rm CMSSW.tgz" <<endl;
    shel << "ls -lhart"<< endl;
    shel << "cd home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src" <<endl;
    shel << "export SCRAM_ARCH=slc6_amd64_gcc530" <<endl;
    shel << "scram b ProjectRename" <<endl;
    shel << "eval `scramv1 runtime -sh`" <<endl;
    shel << "cd TopAnalysis/Configuration/analysis/diLeptonic" <<endl;
    shel << "./install/bin/load_Analysis "<<sampleVar_arg[i].c_str()<<endl;
    //shel << "mkdir -p /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<output_arg[i].c_str()<<endl;
    //shel << "cp "<<output_arg[i].c_str()<<"/*.root /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<output_arg[i].c_str()<<"/." <<endl;
    shel << "gfal-mkdir gsiftp://cms-gridftp.rcac.purdue.edu/store/group/local/cmstop/ajeeta/"<<tree_arg[i].c_str()<<endl; 
    //create a list of files created inside the spinCorrInput folder
    shel << "ls "<<tree_arg[i].c_str()<<"/*.root > filescreated.txt" <<endl;
    shel << "file=filescreated.txt" <<endl;
    shel << "mypwd=$PWD" <<endl;
    shel << "while IFS=\'\' read -r line || [[ -n \"$line\" ]]; do" <<endl;
    shel << "    gfal-copy file:////$mypwd/$line gsiftp://cms-gridftp.rcac.purdue.edu/store/group/local/cmstop/ajeeta/$line " <<endl;
    shel << "    echo \"gfal-copy file:////$mypwd/$line gsiftp://cms-gridftp.rcac.purdue.edu/store/group/local/cmstop/ajeeta/$line\" " <<endl;
    shel << "done < \"$file\" " <<endl;
    shel << "rm -rf CMSSW_8_0_26_patch2" <<endl;
    
  }
}
