#!/bin/sh
python macros/getLeptonSFs.py

# 2016preVFP
mkdir -p fig2016preVFP
mkdir -p fig2016preVFP/controlplots
mkdir -p fig2016preVFP/tables
mkdir -p fig2016preVFP/KinRecoEff
python scripts/ttr2_createEventCountTables.py -y 2016preVFP -m cp
cd Plots_2016preVFP
cp fit_status_*.pdf ../fig2016preVFP/.
cd Nominal
cp eventCountTables.tex ../../fig2016preVFP/tables/.
cp --parents ./*/*.pdf ../../fig2016preVFP/controlplots/
cd ../..
python macros/getBtagEfficiencies.py -y 2016preVFP
mkdir -p fig2016preVFP/bTagEfficiencies
cp bTagEfficiencies/2016preVFP/*.pdf fig2016preVFP/bTagEfficiencies/.
mkdir -p kinRecoRefResolutions/2016preVFP
root -b -q 'macros/getKinRecoRefResolutions.C("2016ULpreVFP")'
cp -r kinRecoRefResolutions/2016preVFP fig2016preVFP/kinRecoRefResolutions
cp -r LeptonSFs/2016preVFP fig2016preVFP/LeptonSFs
./install/bin/kinRecoEfficienciesAndSF -y 2016preVFP -c ee
./install/bin/kinRecoEfficienciesAndSF -y 2016preVFP -c emu
./install/bin/kinRecoEfficienciesAndSF -y 2016preVFP -c mumu
cd KinRecoEffSFPlots_2016preVFP/Nominal
cp --parents ./*/*.pdf ../../fig2016preVFP/KinRecoEff/
cd ../..
./install/bin/kinRecoFullLooseQualityPlots -y 2016preVFP --kr_flag 0

# 2016postVFP
mkdir -p fig2016postVFP
mkdir -p fig2016postVFP/controlplots
mkdir -p fig2016postVFP/tables
mkdir -p fig2016postVFP/KinRecoEff
python scripts/ttr2_createEventCountTables.py -y 2016postVFP -m cp
cd Plots_2016postVFP
cp fit_status_*.pdf ../fig2016postVFP/.
cd Nominal
cp eventCountTables.tex ../../fig2016postVFP/tables/.
cp --parents ./*/*.pdf ../../fig2016postVFP/controlplots/
cd ../..
python macros/getBtagEfficiencies.py -y 2016postVFP
mkdir -p fig2016postVFP/bTagEfficiencies
cp bTagEfficiencies/2016postVFP/*.pdf fig2016postVFP/bTagEfficiencies/.
mkdir -p kinRecoRefResolutions/2016postVFP
root -b -q 'macros/getKinRecoRefResolutions.C("2016ULpostVFP")'
cp -r kinRecoRefResolutions/2016postVFP fig2016postVFP/kinRecoRefResolutions
cp -r LeptonSFs/2016postVFP fig2016postVFP/LeptonSFs
./install/bin/kinRecoEfficienciesAndSF -y 2016postVFP -c ee
./install/bin/kinRecoEfficienciesAndSF -y 2016postVFP -c emu
./install/bin/kinRecoEfficienciesAndSF -y 2016postVFP -c mumu
cd KinRecoEffSFPlots_2016postVFP/Nominal
cp --parents ./*/*.pdf ../../fig2016postVFP/KinRecoEff/
cd ../..
./install/bin/kinRecoFullLooseQualityPlots -y 2016postVFP --kr_flag 0

# 2017
mkdir -p fig2017
mkdir -p fig2017/controlplots
mkdir -p fig2017/tables
mkdir -p fig2017/KinRecoEff
python scripts/ttr2_createEventCountTables.py -y 2017 -m cp
cd Plots_2017
cp fit_status_*.pdf ../fig2017/.
cd Nominal
cp eventCountTables.tex ../../fig2017/tables/.
cp --parents ./*/*.pdf ../../fig2017/controlplots/
cd ../..
python macros/getBtagEfficiencies.py -y 2017
mkdir -p fig2017/bTagEfficiencies
cp bTagEfficiencies/2017/*.pdf fig2017/bTagEfficiencies/.
mkdir -p kinRecoRefResolutions/2017
root -b -q 'macros/getKinRecoRefResolutions.C("2017UL")'
cp -r kinRecoRefResolutions/2017 fig2017/kinRecoRefResolutions
cp -r LeptonSFs/2017 fig2017/LeptonSFs
./install/bin/kinRecoEfficienciesAndSF -y 2017 -c ee
./install/bin/kinRecoEfficienciesAndSF -y 2017 -c emu
./install/bin/kinRecoEfficienciesAndSF -y 2017 -c mumu
cd KinRecoEffSFPlots_2017/Nominal
cp --parents ./*/*.pdf ../../fig2017/KinRecoEff/
cd ../..
./install/bin/kinRecoFullLooseQualityPlots -y 2017 --kr_flag 0

# 2018
mkdir -p fig2018
mkdir -p fig2018/controlplots
mkdir -p fig2018/tables
mkdir -p fig2018/KinRecoEff
python scripts/ttr2_createEventCountTables.py -y 2018 -m cp
cd Plots_2018
cp fit_status_*.pdf ../fig2018/.
cd Nominal
cp eventCountTables.tex ../../fig2018/tables/.
cp --parents ./*/*.pdf ../../fig2018/controlplots/
cd ../..
python macros/getBtagEfficiencies.py -y 2018
mkdir -p fig2018/bTagEfficiencies
cp bTagEfficiencies/2018/*.pdf fig2018/bTagEfficiencies/.
mkdir -p kinRecoRefResolutions/2018
root -b -q 'macros/getKinRecoRefResolutions.C("2018UL")'
cp -r kinRecoRefResolutions/2018 fig2018/kinRecoRefResolutions
cp -r LeptonSFs/2018 fig2018/LeptonSFs
./install/bin/kinRecoEfficienciesAndSF -y 2018 -c ee
./install/bin/kinRecoEfficienciesAndSF -y 2018 -c emu
./install/bin/kinRecoEfficienciesAndSF -y 2018 -c mumu
cd KinRecoEffSFPlots_2018/Nominal
cp --parents ./*/*.pdf ../../fig2018/KinRecoEff/
cd ../..
./install/bin/kinRecoFullLooseQualityPlots -y 2018 --kr_flag 0

# 2016
mkdir -p fig2016
mkdir -p fig2016/controlplots
mkdir -p fig2016/tables
python scripts/ttr2_createEventCountTables.py -y 2016 -m cp
cd Plots_2016
cp fit_status_*.pdf ../fig2016/.
cd Nominal
cp eventCountTables.tex ../../fig2016/tables/.
cp --parents ./*/*.pdf ../../fig2016/controlplots/
cd ../..

# fullRun2UL
mkdir -p figfullRun2UL
mkdir -p figfullRun2UL/controlplots
mkdir -p figfullRun2UL/tables
python scripts/ttr2_createEventCountTables.py -y fullRun2UL -m cp
cd Plots_fullRun2UL
cp --parents TOP_PT/combined/HypToppT.pdf ../figfullRun2UL/
cp fit_status_*.pdf ../figfullRun2UL/.
cd Nominal
cp eventCountTables.tex ../../figfullRun2UL/tables/.
cp --parents ./*/*.pdf ../../figfullRun2UL/controlplots/
cd ../..

cd MiniTreeAnalysis
source setup.sh
python plotting/plotSmearing.py
cd ..
cp -r MiniTreeAnalysis/SmearingPlots/ figfullRun2UL/.

cp -r KinRecoResolutions/Plots figfullRun2UL/KinRecoResolutions

# scp to lxplus
scp -r fig* jthieman@lxplus.cern.ch:AN-20-037/.
