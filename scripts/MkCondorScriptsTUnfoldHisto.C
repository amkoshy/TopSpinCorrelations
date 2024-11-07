#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;
void MkCondorScriptsTUnfoldHisto(void){
  gSystem->mkdir("condor_jobs");
  gSystem->mkdir("condor_jobs/logs");
  ofstream BigJob("RunTUnfoldHisto.sh");
  
  string chan[4]={"ee","emu","mumu","combined"}; 
  string syst[64]={
                   "Nominal",
                   "JES_UP", "JES_DOWN",
                   "JER_UP", "JER_DOWN",
                   "BG_UP", "BG_DOWN",
                   "DY_UP", "DY_DOWN",
		   /*"JESAbsoluteStat_UP", "JESAbsoluteStat_DOWN", 
		   "JESAbsoluteScale_UP", "JESAbsoluteScale_DOWN", 
		   "JESAbsoluteFlavMap_UP", "JESAbsoluteFlavMap_DOWN", 
		   "JESAbsoluteMPFBias_UP", "JESAbsoluteMPFBias_DOWN", 
		   "JESFragmentation_UP", "JESFragmentation_DOWN", 
		   "JESSinglePionECAL_UP", "JESSinglePionECAL_DOWN", 
		   "JESSinglePionHCAL_UP", "JESSinglePionHCAL_DOWN", 
		   "JESFlavorQCD_UP", "JESFlavorQCD_DOWN", 
		   "JESTimePtEta_UP", "JESTimePtEta_DOWN", 
		   "JESRelativeJEREC1_UP", "JESRelativeJEREC1_DOWN", 
		   "JESRelativeJEREC2_UP", "JESRelativeJEREC2_DOWN", 
		   "JESRelativeJERHF_UP", "JESRelativeJERHF_DOWN", 
		   "JESRelativePtBB_UP", "JESRelativePtBB_DOWN", 
		   "JESRelativePtEC1_UP", "JESRelativePtEC1_DOWN", 
		   "JESRelativePtEC2_UP", "JESRelativePtEC2_DOWN", 
		   "JESRelativePtHF_UP", "JESRelativePtHF_DOWN", 
		   "JESRelativeFSR_UP", "JESRelativeFSR_DOWN", 
		   "JESRelativeStatFSR_UP", "JESRelativeStatFSR_DOWN", 
		   "JESRelativeStatEC_UP", "JESRelativeStatEC_DOWN", 
		   "JESRelativeStatHF_UP", "JESRelativeStatHF_DOWN", 
		   "JESPileUPDataMC_UP", "JESPileUPDataMC_DOWN", 
		   "JESPileUPPtRef_UP", "JESPileUPPtRef_DOWN", 
		   "JESPileUPPtEC1_UP", "JESPileUPPtEC1_DOWN", 
		   "JESPileUPPtEC2_UP", "JESPileUPPtEC2_DOWN", 
		   "JESPileUPPtHF_UP", "JESPileUPPtHF_DOWN", 
		   "JESPileUPPtBB_UP", "JESPileUPPtBB_DOWN",*/
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
		   "BTAG_LJET_PT_UP", "BTAG_LJET_PT_DOWN", 
		   "BTAG_LJET_ETA_UP", "BTAG_LJET_ETA_DOWN",
  		   "PDF_ALPHAS_UP", "PDF_ALPHAS_DOWN",
		   "MESCALE_UP", "MESCALE_DOWN", 
		   "MEFACSCALE_UP", "MEFACSCALE_DOWN",
		   "MERENSCALE_UP", "MERENSCALE_DOWN",
		   "PSISRSCALE_UP", "PSISRSCALE_DOWN",
		   "PSFSRSCALE_UP", "PSFSRSCALE_DOWN",
		   "UETUNE_UP", "UETUNE_DOWN", 
		   "MASS_UP", "MASS_DOWN",
		   "BFRAG_UP", "BFRAG_DOWN",
		   "BFRAG_CENTRAL", "BFRAG_PETERSON", 
		   "BSEMILEP_UP", "BSEMILEP_DOWN", 
		   "ERDON", "ERDONRETUNE", "GLUONMOVETUNE",
		   "MATCH_UP", "MATCH_DOWN",
		   "POWHEG", "POWHEGHERWIG", "MCATNLO","AMCATNLOFXFX"};
  
  vector<string>sampleVar_arg;
  vector<string>output_arg;

  for(int j=0; j<64; j++){
    for(int i=0; i<4; i++){
      sampleVar_arg.push_back(" "+syst[j]+" "+chan[i]);    
      output_arg.push_back(syst[j]+"/"+chan[i]);
    }
  }
  
  BigJob << "#!/bin/bash" <<endl;
  BigJob << "tar -X excludefromtar.txt -zcvf CMSSWTUnfold.tgz /home/akhatiwa/ttbarSpinCor_fullInstallation/CMSSW_8_0_26_patch2/" <<endl;
  BigJob << "voms-proxy-init -voms cms -valid 100:00:00" <<endl;
  for (int i=0; i<sampleVar_arg.size(); i++){

    char ShellName[600];
    sprintf(ShellName, "condor_jobs/condor_TUnfoldHisto_%i.sh",i);
    
    char JDLName[600];
    sprintf(JDLName, "condor_jobs/AnaCondor_TUnfoldHisto_%i.jdl",i);
    
    ofstream jdl(JDLName);
    ofstream shel(ShellName);
    BigJob << "condor_submit " << JDLName << endl;

    //Make the JDL file:
    
    jdl << "Universe                = vanilla" << endl;
    jdl << "Environment             = CONDORJOBID=$(Process)" << endl;
    jdl << "notification            = Error" << endl;
    jdl << "X509UserProxy = /tmp/x509up_u587091" << endl;
    jdl << "should_transfer_files   = YES" << endl;
    jdl << "Transfer_input_files = condor_jobs/condor_TUnfoldHisto_"<<i<<".sh" << endl;
    jdl << "Executable = /home/akhatiwa/ttbarSpinCor_fullInstallation/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<ShellName << endl;
    jdl << "when_to_transfer_output = ON_EXIT" << endl;
    jdl << "Output = /home/akhatiwa/ttbarSpinCor_fullInstallation/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/condor_jobs/logs/batch_"<<i<<"_$(cluster)_$(process).stdout" << endl;
    jdl << "Error  = /home/akhatiwa/ttbarSpinCor_fullInstallation/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/condor_jobs/logs/batch_"<<i<<"_$(cluster)_$(process).stderr" << endl;
    jdl << "Log    = /home/akhatiwa/ttbarSpinCor_fullInstallation/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/condor_jobs/logs/batch_"<<i<<"_$(cluster)_$(process).condor" << endl;
    jdl << "Queue 1" <<endl; 
    
    //Make the Shell script: 
    shel << "#!/bin/sh" << endl;
    shel << "source /cvmfs/cms.cern.ch/cmsset_default.sh" << endl;
    shel << "cp /home/akhatiwa/ttbarSpinCor_fullInstallation/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/CMSSWTUnfold.tgz ." <<endl;
    shel << "tar -zxf CMSSWTUnfold.tgz" <<endl;
    shel << "rm CMSSWTUnfold.tgz" <<endl;
    shel << "ls -lhart"<< endl;
    shel << "cd /home/akhatiwa/ttbarSpinCor_fullInstallation/CMSSW_8_0_26_patch2/src" <<endl;
    shel << "export SCRAM_ARCH=slc6_amd64_gcc530" <<endl;
    shel << "scram b ProjectRename" <<endl;
    shel << "eval `scramv1 runtime -sh`" <<endl;
    shel << "cd TopAnalysis/Configuration/analysis/diLeptonic" <<endl;
    shel << "./install/bin/loadTUnfoldHisto"<<sampleVar_arg[i].c_str()<<endl;
    //shel << "mkdir -p /home/akhatiwa/ttbarSpinCor_fullInstallation/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/UnfoldingHistos/"<<output_arg[i].c_str()<<endl;
    //shel << "cp "<<output_arg[i].c_str()<<"/*.root /home/akhatiwa/ttbarSpinCor_fullInstallation/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<output_arg[i].c_str()<<"/." <<endl;
    //shel << "rm -rf CMSSW_8_0_26_patch2" <<endl;
    shel << "cp -r binning /home/akhatiwa/ttbarSpinCor_fullInstallation/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/." <<endl;
    shel << "gfal-mkdir -p gsiftp://cms-gridftp.rcac.purdue.edu/store/group/local/cmstop/ajeeta/UnfoldingHistos_DataSymmetrized/"<<output_arg[i].c_str()<<endl;
    shel << "cd UnfoldingHistos" <<endl;
    shel << "ls "<<output_arg[i].c_str()<<"/*.root > filescreated.txt" <<endl;
    shel << "file=filescreated.txt" <<endl;
    shel << "mypwd=$PWD" <<endl;
    shel << "while IFS=\'\' read -r line || [[ -n \"$line\" ]]; do" <<endl;
    shel << "    gfal-copy -f file:////$mypwd/$line gsiftp://cms-gridftp.rcac.purdue.edu/store/group/local/cmstop/ajeeta/UnfoldingHistos_DataSymmetrized/$line " <<endl;
    shel << "    echo \"gfal-copy file:////$mypwd/$line gsiftp://cms-gridftp.rcac.purdue.edu/store/group/local/cmstop/ajeeta/UnfoldingHistos_DataSymmetrized/$line\" " <<endl;
    shel << "done < \"$file\" " <<endl;
    shel << "rm -rf CMSSW_8_0_26_patch2" <<endl;
    
  }
}
