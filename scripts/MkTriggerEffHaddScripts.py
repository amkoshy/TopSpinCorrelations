#!/usr/bin/env python

eras = [
"2016preVFP_UL",
"2016postVFP_UL",
"2017_UL",
"2018_UL",
]

nSignalFiles = -1

Mkhadd_triggereff_haddfile = open("scripts/hadd_triggereff.sh","w")
Mkhadd_triggereff_haddfile.write("#!/bin/bash \n")

for era in eras:

    if (era == "2016preVFP_UL"):
        nSignalFiles = 20
    elif (era == "2016postVFP_UL"):
        nSignalFiles = 20
    elif (era == "2017_UL"):
        nSignalFiles = 30
    elif (era == "2018_UL"):
        nSignalFiles = 40


    Mkhadd_triggereff_haddfile.write("hadd -f -j 12 triggereff/output_hists_MC_ttbarsignalplustau_fromDilepton_"+era+".root ")

    for i in range (0,nSignalFiles):

        Mkhadd_triggereff_haddfile.write("triggereff/output_hists_MC_ttbarsignalplustau_fromDilepton_"+str(i)+"_"+era+".root ")

    Mkhadd_triggereff_haddfile.write("\n")


Mkhadd_triggereff_haddfile.close()
