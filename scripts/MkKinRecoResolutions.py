#!/usr/bin/env python

eras = [
"2016preVFP",
"2016postVFP",
"2017",
"2018",
]

MkKinRecoResolutions_scriptfile = open("scripts/KinRecoResolutions.sh","w")
MkKinRecoResolutions_scriptfile.write("#!/bin/bash \n")

for era in eras:

    if (era == "2016preVFP" or era == "2016postVFP" or era == "2017" ):

        for i in range (0,20):

            #MkKinRecoResolutions_scriptfile.write("root -l -q src/KinRecoResolutions.cc++\'(\""+era+"\", \"prompt\", \""+str(i)+"\")\' \n")
            MkKinRecoResolutions_scriptfile.write("./install/bin/KinRecoResolutions -y "+era+" -t prompt -n "+str(i)+" \n")

        for i in range (0,20):

            #MkKinRecoResolutions_scriptfile.write("root -l -q src/KinRecoResolutions.cc++\'(\""+era+"\", \"viatau\", \""+str(i)+"\")\' \n")
            MkKinRecoResolutions_scriptfile.write("./install/bin/KinRecoResolutions -y "+era+" -t viatau -n "+str(i)+" \n")

    elif (era == "2018" ):

        for i in range (0,40):

            #MkKinRecoResolutions_scriptfile.write("root -l -q src/KinRecoResolutions.cc++\'(\""+era+"\", \"prompt\", \""+str(i)+"\")\' \n")
            MkKinRecoResolutions_scriptfile.write("./install/bin/KinRecoResolutions -y "+era+" -t prompt -n "+str(i)+" \n")

        for i in range (0,40):

            #MkKinRecoResolutions_scriptfile.write("root -l -q src/KinRecoResolutions.cc++\'(\""+era+"\", \"viatau\", \""+str(i)+"\")\' \n")
            MkKinRecoResolutions_scriptfile.write("./install/bin/KinRecoResolutions -y "+era+" -t viatau -n "+str(i)+" \n")

MkKinRecoResolutions_scriptfile.close()


MkKinRecoResolutions_haddfile = open("scripts/hadd_KinRecoResolutions.sh","w")
MkKinRecoResolutions_haddfile.write("#!/bin/bash \n")

for era in eras:

    if (era == "2016preVFP" or era == "2016postVFP" or era == "2017" ):

        MkKinRecoResolutions_haddfile.write("hadd -f -j 12 KinRecoResolutions/output_hists_prompt_"+era+".root ")

        for i in range (0,20):

            MkKinRecoResolutions_haddfile.write("KinRecoResolutions/output_hists_prompt_"+era+"_"+str(i)+".root ")

        MkKinRecoResolutions_haddfile.write("\n")


        MkKinRecoResolutions_haddfile.write("hadd -f -j 12 KinRecoResolutions/output_hists_viatau_"+era+".root ")

        for i in range (0,20):

            MkKinRecoResolutions_haddfile.write("KinRecoResolutions/output_hists_viatau_"+era+"_"+str(i)+".root ")

        MkKinRecoResolutions_haddfile.write("\n")


    elif (era == "2018" ):

        MkKinRecoResolutions_haddfile.write("hadd -f -j 12 KinRecoResolutions/output_hists_prompt_"+era+".root ")

        for i in range (0,40):

            MkKinRecoResolutions_haddfile.write("KinRecoResolutions/output_hists_prompt_"+era+"_"+str(i)+".root ")

        MkKinRecoResolutions_haddfile.write("\n")


        MkKinRecoResolutions_haddfile.write("hadd -f -j 12 KinRecoResolutions/output_hists_viatau_"+era+".root ")

        for i in range (0,40):

            MkKinRecoResolutions_haddfile.write("KinRecoResolutions/output_hists_viatau_"+era+"_"+str(i)+".root ")

        MkKinRecoResolutions_haddfile.write("\n")

MkKinRecoResolutions_haddfile.close()
