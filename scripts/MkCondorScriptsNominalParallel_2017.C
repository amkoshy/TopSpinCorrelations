#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;
void MkCondorScriptsNominalParallel_2017(void){
  //Options for file output
  string output_location_hadoop = "gsiftp://cms-gridftp.rcac.purdue.edu/store/group/local/cmstop/jthiema/";

  const int nChan = 3;
  const int nSingleChan = 2;
  const int nSample = 10;


  gSystem->mkdir("condor_jobs");
  gSystem->mkdir("condor_jobs/logs");
  ofstream BigJob("RunNominalParallel_2017.sh");

  const char* env_cmsswbase = std::getenv("CMSSW_BASE");

  //string mode="";
  string mode="-m spinCorr";
  string chan[nChan]={"ee","emu","mumu"}; 
  string singlechan[nSingleChan]={"se","smu"}; 
  string sample[nSample]={"zz", "qcd", "tw.root", "ttbarbg_fromDilepton.root", "ttbarbg_fromHadronic.root","ttbarbg_fromLjets.root", "wtol", "ww", "wz", "ttbarW"};
  
  vector<string>sampleVar_arg;
  vector<string>output_arg;
  vector<string>tree_arg;
  
  for(int i=0; i<nChan; i++){
    sampleVar_arg.push_back("-f ttbarsignalplustau_fromDilepton.root "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarsignalplustau_fromDilepton.root "+mode+" -c "+chan[i]+" --signalviatau");
    sampleVar_arg.push_back("-f dy -d 11 "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f dy -d 13 "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f dy -d 15 "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarZ "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f "+chan[i]+"_run2017B "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f "+chan[i]+"_run2017C "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f "+chan[i]+"_run2017D "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f "+chan[i]+"_run2017E "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f "+chan[i]+"_run2017F "+mode+" -c "+chan[i]);
       
    for(int j=0; j<11; j++){
      output_arg.push_back("selectionRoot_2017/Nominal/"+chan[i]);
      tree_arg.push_back("spinCorrInput/Nominal/"+chan[i]);
    }
    
    for(int j=0; j<nSingleChan; j++){
      sampleVar_arg.push_back("-f "+singlechan[j]+"_run2017B "+mode+" -c "+chan[i]);
      sampleVar_arg.push_back("-f "+singlechan[j]+"_run2017C "+mode+" -c "+chan[i]);
      sampleVar_arg.push_back("-f "+singlechan[j]+"_run2017D "+mode+" -c "+chan[i]);
      sampleVar_arg.push_back("-f "+singlechan[j]+"_run2017E "+mode+" -c "+chan[i]);
      sampleVar_arg.push_back("-f "+singlechan[j]+"_run2017F "+mode+" -c "+chan[i]);

      for(int k=0; k<5; k++){
        output_arg.push_back("selectionRoot_2017/Nominal/"+chan[i]);
        tree_arg.push_back("spinCorrInput/Nominal/"+chan[i]);
      }
      
    }
    
        for(int j=0; j<nSample; j++){
      sampleVar_arg.push_back("-f "+sample[j]+" "+mode+" -c "+chan[i]);
      output_arg.push_back("selectionRoot_2017/Nominal/"+chan[i]);
      tree_arg.push_back("spinCorrInput/Nominal/"+chan[i]);  
    }
    
  }
  
  BigJob << "#!/bin/bash" <<endl;
  BigJob << "voms-proxy-init -voms cms -valid 100:00:00" <<endl;
  BigJob << "tar -X excludefromtar.txt -zcvf CMSSW.tgz ../../../../../../CMSSW_10_2_22" <<endl;
  BigJob << "export X509_USER_CERT=`voms-proxy-info -path`" <<endl;
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
    jdl << "Requirements            = (Distro == \"RHEL7\")" << endl;
    jdl << "request_disk            = 200000MB" << endl;
    //if proxy is not specified it will use environment variable $X509_USER_CERT
    jdl << "X509UserProxy = /tmp/x509up_u607473" << endl;
    jdl << "should_transfer_files   = YES" << endl;
    jdl << "Transfer_input_files = condor_jobs/condor_doit_"<<i<<".sh" << endl;
    jdl << "Executable = "<<env_cmsswbase<<"/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<ShellName << endl;
    jdl << "when_to_transfer_output = ON_EXIT" << endl;
    jdl << "Output = "<<env_cmsswbase<<"/src/TopAnalysis/Configuration/analysis/diLeptonic/condor_jobs/logs/batch_"<<i<<"_$(cluster)_$(process).stdout" << endl;
    jdl << "Error  = "<<env_cmsswbase<<"/src/TopAnalysis/Configuration/analysis/diLeptonic/condor_jobs/logs/batch_"<<i<<"_$(cluster)_$(process).stderr" << endl;
    jdl << "Log    = "<<env_cmsswbase<<"/src/TopAnalysis/Configuration/analysis/diLeptonic/condor_jobs/logs/batch_"<<i<<"_$(cluster)_$(process).condor" << endl;
    jdl << "Queue 1" <<endl; 
    
    //Make the Shell script:
    shel << "#!/bin/sh" << endl;
    shel << "source /cvmfs/cms.cern.ch/cmsset_default.sh" << endl;
    shel << "cp "<<env_cmsswbase<<"/src/TopAnalysis/Configuration/analysis/diLeptonic/CMSSW.tgz ." <<endl;
    shel << "tar -zxf CMSSW.tgz" <<endl;
    shel << "rm CMSSW.tgz" <<endl;
    shel << "ls -lhart"<< endl;
    shel << "cd CMSSW_10_2_22/src" <<endl;
    shel << "export SCRAM_ARCH=slc7_amd64_gcc700" <<endl;
    shel << "scram b ProjectRename" <<endl;
    shel << "eval `scramv1 runtime -sh`" <<endl;
    shel << "cd TopAnalysis/Configuration/analysis/diLeptonic" <<endl;
    shel << "./install/bin/load_Analysis "<<sampleVar_arg[i].c_str()<<endl;
    shel << "mkdir -p "<<env_cmsswbase<<"/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<output_arg[i].c_str()<<endl;
    shel << "cp "<<output_arg[i].c_str()<<"/*.root "<<env_cmsswbase<<"/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<output_arg[i].c_str()<<"/." <<endl;
    shel << "cp "<<tree_arg[i].c_str()<<"/*.root "<<env_cmsswbase<<"/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<tree_arg[i].c_str()<<"/." <<endl;



    //Save minitree in hadoop
    //    shel << "gfal-mkdir -p "<<output_location_hadoop.c_str()<<tree_arg[i].c_str()<<endl;
    //create a list of files created inside the spinCorrInput folder
    //    shel << "ls "<<tree_arg[i].c_str()<<"/*.root > filescreated.txt" <<endl;
    //    shel << "file=filescreated.txt" <<endl;
    //    shel << "mypwd=$PWD" <<endl;
    //    shel << "while IFS=\'\' read -r line || [[ -n \"$line\" ]]; do" <<endl;
    //    shel << "    gfal-copy -f file:////$mypwd/$line "<<output_location_hadoop.c_str()<<"$line " <<endl;
    //    shel << "    echo \"gfal-copy -f file:////$mypwd/$line "<<output_location_hadoop.c_str()<<"$line\" " <<endl;
    //    shel << "done < \"$file\" " <<endl;
    //    shel << "rm -rf CMSSW_10_2_22" <<endl;    
  }
}
