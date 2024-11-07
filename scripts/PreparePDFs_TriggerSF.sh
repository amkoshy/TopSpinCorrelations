#!/bin/sh
cd triggereff

mkdir -p fig_2016preVFP_Purdue
cd 2016preVFP_UL
cp *.pdf ../fig_2016preVFP_Purdue/.
cp nominal/*.pdf ../fig_2016preVFP_Purdue/.
cp effData/nominal/*.pdf ../fig_2016preVFP_Purdue/.
cp effMC/nominal/*.pdf ../fig_2016preVFP_Purdue/.
cp low3_njets/h2D_*.pdf ../fig_2016preVFP_Purdue/.
cp high3_njets/h2D_*.pdf ../fig_2016preVFP_Purdue/.
cp low25_nvtx/h2D_*.pdf ../fig_2016preVFP_Purdue/.
cp high25_nvtx/h2D_*.pdf ../fig_2016preVFP_Purdue/.
cp --parents alpha/nominal/*.pdf ../fig_2016preVFP_Purdue/
cd ..

mkdir -p fig_2016postVFP_Purdue
cd 2016postVFP_UL
cp *.pdf ../fig_2016postVFP_Purdue/.
cp nominal/*.pdf ../fig_2016postVFP_Purdue/.
cp effData/nominal/*.pdf ../fig_2016postVFP_Purdue/.
cp effMC/nominal/*.pdf ../fig_2016postVFP_Purdue/.
cp low3_njets/h2D_*.pdf ../fig_2016postVFP_Purdue/.
cp high3_njets/h2D_*.pdf ../fig_2016postVFP_Purdue/.
cp low25_nvtx/h2D_*.pdf ../fig_2016postVFP_Purdue/.
cp high25_nvtx/h2D_*.pdf ../fig_2016postVFP_Purdue/.
cp --parents alpha/nominal/*.pdf ../fig_2016postVFP_Purdue/
cd ..

mkdir -p fig_2017_Purdue
cd 2017_UL
cp *.pdf ../fig_2017_Purdue/.
cp nominal/*.pdf ../fig_2017_Purdue/.
cp effData/nominal/*.pdf ../fig_2017_Purdue/.
cp effMC/nominal/*.pdf ../fig_2017_Purdue/.
cp low3_njets/h2D_*.pdf ../fig_2017_Purdue/.
cp high3_njets/h2D_*.pdf ../fig_2017_Purdue/.
cp low35_nvtx/h2D_*.pdf ../fig_2017_Purdue/.
cp high35_nvtx/h2D_*.pdf ../fig_2017_Purdue/.
cp --parents alpha/nominal/*.pdf ../fig_2017_Purdue/
cd ..

mkdir -p fig_2018_Purdue
cd 2018_UL
cp *.pdf ../fig_2018_Purdue/.
cp nominal/*.pdf ../fig_2018_Purdue/.
cp effData/nominal/*.pdf ../fig_2018_Purdue/.
cp effMC/nominal/*.pdf ../fig_2018_Purdue/.
cp low3_njets/h2D_*.pdf ../fig_2018_Purdue/.
cp high3_njets/h2D_*.pdf ../fig_2018_Purdue/.
cp low35_nvtx/h2D_*.pdf ../fig_2018_Purdue/.
cp high35_nvtx/h2D_*.pdf ../fig_2018_Purdue/.
cp --parents alpha/nominal/*.pdf ../fig_2018_Purdue/
cd ..

cd ..
scp -r triggereff/fig_201*_Purdue jthieman@lxplus.cern.ch:TrigEffAN/chapters/fig/.
