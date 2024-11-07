#!/usr/bin/env python
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("commandlist")
args = parser.parse_args()

CMSSW_BASE = str(os.getenv("CMSSW_BASE"))


i = 0

if os.path.exists("PlotsSlurmJobs"):
    for job in os.listdir("PlotsSlurmJobs"):
        if i < int(job.strip('SlurmJob_').strip('.sh')): 
            i = int(job.strip('SlurmJob_').strip('.sh'))+1

else: os.makedirs("PlotsSlurmJobs")

if not os.path.exists("nohuplogs"): os.makedirs("nohuplogs")
if not os.path.exists("plotslogs"): os.makedirs("plotslogs")


commandlistfilepath = args.commandlist
commandlistfile     = open(commandlistfilepath, "r")
lines = commandlistfile.readlines()


runfile = open("RunPlotsSlurm_"+os.path.basename(commandlistfilepath).rsplit(".",1)[0]+".sh", "w")

runfile.write("#!/bin/sh")
runfile.write("\n")

for j in range(0,len(lines),10) :

    sh_file = "PlotsSlurmJobs/SlurmJob_" + str(i) + ".sh"
    
    with open(sh_file, "w") as cfg :
        cfg.write("#!/bin/sh")
        cfg.write("\n")
        cfg.write("#SBATCH  -A standby")
        cfg.write("\n")
        cfg.write("#SBATCH --ntasks=1")
        cfg.write("\n")
        cfg.write("#SBATCH --cpus-per-task=1")
        cfg.write("\n")
        cfg.write("#SBATCH --time=3:59:59")
        cfg.write("\n")
        cfg.write("#SBATCH --mem-per-cpu=8000")
        cfg.write("\n")
        cfg.write("cd " + CMSSW_BASE + "/src/")
        cfg.write("\n")
        cfg.write("source /cvmfs/cms.cern.ch/cmsset_default.sh")
        cfg.write("\n")
        cfg.write("export BOOST_ROOT=/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/boost/1.63.0-gnimlf")
        cfg.write("\n")
        cfg.write("export SCRAM_ARCH=slc7_amd64_gcc820")
        cfg.write("\n")
        cfg.write("eval `scramv1 runtime -sh`")
        cfg.write("\n")
        cfg.write("cd TopAnalysis/Configuration/analysis/diLeptonic")
        cfg.write("\n")
        if j+0 < len(lines):
            cfg.write(str(lines[j].strip("\n"))+"\n")
        if j+1 < len(lines):
            cfg.write(str(lines[j+1].strip("\n"))+"\n")
        if j+2 < len(lines):
            cfg.write(str(lines[j+2].strip("\n"))+"\n")
        if j+3 < len(lines):
            cfg.write(str(lines[j+3].strip("\n"))+"\n")
        if j+4 < len(lines):
            cfg.write(str(lines[j+4].strip("\n"))+"\n")
        if j+5 < len(lines):
            cfg.write(str(lines[j+5].strip("\n"))+"\n")
        if j+6 < len(lines):
            cfg.write(str(lines[j+6].strip("\n"))+"\n")
        if j+7 < len(lines):
            cfg.write(str(lines[j+7].strip("\n"))+"\n")
        if j+8 < len(lines):
            cfg.write(str(lines[j+8].strip("\n"))+"\n")
        if j+9 < len(lines):
            cfg.write(str(lines[j+9].strip("\n"))+"\n")

    runfile.write("sbatch PlotsSlurmJobs/SlurmJob_" + str(i) + ".sh")
    runfile.write("\n")

    i   += 1
