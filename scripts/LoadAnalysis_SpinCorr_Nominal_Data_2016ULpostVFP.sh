#!/bin/bash

nohup ./install/bin/load_Analysis -f ee_run2016F2 -m spinCorr cp kinReco kinRecoQualityStudies -c ee &> nohuplogs/SpinCorr_run2016F2_ee.out
nohup ./install/bin/load_Analysis -f ee_run2016G -m spinCorr cp kinReco kinRecoQualityStudies -c ee &> nohuplogs/SpinCorr_run2016G_ee.out
nohup ./install/bin/load_Analysis -f ee_run2016H -m spinCorr cp kinReco kinRecoQualityStudies -c ee &> nohuplogs/SpinCorr_run2016H_ee.out

nohup ./install/bin/load_Analysis -f emu_run2016F2 -m spinCorr cp kinReco kinRecoQualityStudies -c emu &> nohuplogs/SpinCorr_run2016F2_emu.out
nohup ./install/bin/load_Analysis -f emu_run2016G -m spinCorr cp kinReco kinRecoQualityStudies -c emu &> nohuplogs/SpinCorr_run2016G_emu.out
nohup ./install/bin/load_Analysis -f emu_run2016H -m spinCorr cp kinReco kinRecoQualityStudies -c emu &> nohuplogs/SpinCorr_run2016H_emu.out

nohup ./install/bin/load_Analysis -f mumu_run2016F2 -m spinCorr cp kinReco kinRecoQualityStudies -c mumu &> nohuplogs/SpinCorr_run2016F2_mumu.out
nohup ./install/bin/load_Analysis -f mumu_run2016G -m spinCorr cp kinReco kinRecoQualityStudies -c mumu &> nohuplogs/SpinCorr_run2016G_mumu.out
nohup ./install/bin/load_Analysis -f mumu_run2016H -m spinCorr cp kinReco kinRecoQualityStudies -c mumu &> nohuplogs/SpinCorr_run2016H_mumu.out

nohup ./install/bin/load_Analysis -f se_run2016F2 -m spinCorr cp kinReco kinRecoQualityStudies -c ee &> nohuplogs/SpinCorr_se_run2016F2_ee.out
nohup ./install/bin/load_Analysis -f se_run2016F2 -m spinCorr cp kinReco kinRecoQualityStudies -c emu &> nohuplogs/SpinCorr_se_run2016F2_emu.out
nohup ./install/bin/load_Analysis -f se_run2016F2 -m spinCorr cp kinReco kinRecoQualityStudies -c mumu &> nohuplogs/SpinCorr_se_run2016F2_mumu.out
nohup ./install/bin/load_Analysis -f se_run2016G -m spinCorr cp kinReco kinRecoQualityStudies -c ee &> nohuplogs/SpinCorr_se_run2016G_ee.out
nohup ./install/bin/load_Analysis -f se_run2016G -m spinCorr cp kinReco kinRecoQualityStudies -c emu &> nohuplogs/SpinCorr_se_run2016G_emu.out
nohup ./install/bin/load_Analysis -f se_run2016G -m spinCorr cp kinReco kinRecoQualityStudies -c mumu &> nohuplogs/SpinCorr_se_run2016G_mumu.out
nohup ./install/bin/load_Analysis -f se_run2016H -m spinCorr cp kinReco kinRecoQualityStudies -c ee &> nohuplogs/SpinCorr_se_run2016H_ee.out
nohup ./install/bin/load_Analysis -f se_run2016H -m spinCorr cp kinReco kinRecoQualityStudies -c emu &> nohuplogs/SpinCorr_se_run2016H_emu.out
nohup ./install/bin/load_Analysis -f se_run2016H -m spinCorr cp kinReco kinRecoQualityStudies -c mumu &> nohuplogs/SpinCorr_se_run2016H_mumu.out

nohup ./install/bin/load_Analysis -f smu_run2016F2 -m spinCorr cp kinReco kinRecoQualityStudies -c ee &> nohuplogs/SpinCorr_smu_run2016F2_ee.out
nohup ./install/bin/load_Analysis -f smu_run2016F2 -m spinCorr cp kinReco kinRecoQualityStudies -c emu &> nohuplogs/SpinCorr_smu_run2016F2_emu.out
nohup ./install/bin/load_Analysis -f smu_run2016F2 -m spinCorr cp kinReco kinRecoQualityStudies -c mumu &> nohuplogs/SpinCorr_smu_run2016F2_mumu.out
nohup ./install/bin/load_Analysis -f smu_run2016G -m spinCorr cp kinReco kinRecoQualityStudies -c ee &> nohuplogs/SpinCorr_smu_run2016G_ee.out
nohup ./install/bin/load_Analysis -f smu_run2016G -m spinCorr cp kinReco kinRecoQualityStudies -c emu &> nohuplogs/SpinCorr_smu_run2016G_emu.out
nohup ./install/bin/load_Analysis -f smu_run2016G -m spinCorr cp kinReco kinRecoQualityStudies -c mumu &> nohuplogs/SpinCorr_smu_run2016G_mumu.out
nohup ./install/bin/load_Analysis -f smu_run2016H -m spinCorr cp kinReco kinRecoQualityStudies -c ee &> nohuplogs/SpinCorr_smu_run2016H_ee.out
nohup ./install/bin/load_Analysis -f smu_run2016H -m spinCorr cp kinReco kinRecoQualityStudies -c emu &> nohuplogs/SpinCorr_smu_run2016H_emu.out
nohup ./install/bin/load_Analysis -f smu_run2016H -m spinCorr cp kinReco kinRecoQualityStudies -c mumu &> nohuplogs/SpinCorr_smu_run2016H_mumu.out
