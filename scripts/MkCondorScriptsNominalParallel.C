#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;
void MkCondorScriptsNominalParallel(void){
  gSystem->mkdir("condor_jobs");
  gSystem->mkdir("condor_jobs/logs");
  ofstream BigJob("RunNominalParallel.sh");

  //string mode="";
  string mode="-m spinCorr";
  string chan[3]={"ee","emu","mumu"}; 
  string singlechan[2]={"se","smu"}; 
  string sample[9]={"zz", "qcd", "tw.root", "ttbarbg.root", "wtol", "ww", "wz", "ttbarW", "ttgjets"};
  
  vector<string>sampleVar_arg;
  vector<string>output_arg;
  vector<string>tree_arg;
  
  for(int i=0; i<3; i++){
    sampleVar_arg.push_back("-f ttbarsignalplustau.root "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarsignalplustau.root "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f dy -d 11 "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f dy -d 13 "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f dy -d 15 "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarZ "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f "+chan[i]+"_run2016B "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f "+chan[i]+"_run2016C "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f "+chan[i]+"_run2016D "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f "+chan[i]+"_run2016E "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f "+chan[i]+"_run2016F1 "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f "+chan[i]+"_run2016F2 "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f "+chan[i]+"_run2016G "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f "+chan[i]+"_run2016H "+mode+" -c "+chan[i]);
    
    for(int j=0; j<14; j++){
      output_arg.push_back("selectionRoot/Nominal/"+chan[i]);
      tree_arg.push_back("spinCorrInput/Nominal/"+chan[i]);
    }
    
    for(int j=0; j<2; j++){
      sampleVar_arg.push_back("-f "+singlechan[j]+"_run2016B "+mode+" -c "+chan[i]);
      sampleVar_arg.push_back("-f "+singlechan[j]+"_run2016C "+mode+" -c "+chan[i]);
      sampleVar_arg.push_back("-f "+singlechan[j]+"_run2016D "+mode+" -c "+chan[i]);
      sampleVar_arg.push_back("-f "+singlechan[j]+"_run2016E "+mode+" -c "+chan[i]);
      sampleVar_arg.push_back("-f "+singlechan[j]+"_run2016F1 "+mode+" -c "+chan[i]);
      sampleVar_arg.push_back("-f "+singlechan[j]+"_run2016F2 "+mode+" -c "+chan[i]);
      sampleVar_arg.push_back("-f "+singlechan[j]+"_run2016G "+mode+" -c "+chan[i]);
      sampleVar_arg.push_back("-f "+singlechan[j]+"_run2016H "+mode+" -c "+chan[i]);
      
      for(int k=0; k<8; k++){
        output_arg.push_back("selectionRoot/Nominal/"+chan[i]);
        tree_arg.push_back("spinCorrInput/Nominal/"+chan[i]);
      }
      
    }
    
    for(int j=0; j<9; j++){
      sampleVar_arg.push_back("-f "+sample[j]+" "+mode+" -c "+chan[i]);
      output_arg.push_back("selectionRoot/Nominal/"+chan[i]);
      tree_arg.push_back("spinCorrInput/Nominal/"+chan[i]);  
  }
    
  }
  
  BigJob << "#!/bin/bash" <<endl;
  BigJob << "tar -X excludefromtar.txt -zcvf CMSSW.tgz /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/" <<endl;
  BigJob << "voms-proxy-init -voms cms -valid 100:00:00" <<endl;
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
    shel << "mkdir -p /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<output_arg[i].c_str()<<endl;
    shel << "cp "<<output_arg[i].c_str()<<"/*.root /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<output_arg[i].c_str()<<"/." <<endl;
    shel << "gfal-mkdir gsiftp://cms-gridftp.rcac.purdue.edu/store/group/local/cmstop/ajeeta/"<<tree_arg[i].c_str()<<endl;
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
