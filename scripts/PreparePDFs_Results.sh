#!/bin/sh

# 2016preVFP
mkdir -p fig2016preVFP
mkdir -p fig2016preVFP/unfolding
cd TUnfoldResults_2016preVFP/Nominal
cp --parents ./*/Plots/*.pdf ../../fig2016preVFP/unfolding/
cd ../..

# 2016postVFP
mkdir -p fig2016postVFP
mkdir -p fig2016postVFP/unfolding
cd TUnfoldResults_2016postVFP/Nominal
cp --parents ./*/Plots/*.pdf ../../fig2016postVFP/unfolding/
cd ../..

# 2017
mkdir -p fig2017
mkdir -p fig2017/unfolding
cd TUnfoldResults_2017/Nominal
cp --parents ./*/Plots/*.pdf ../../fig2017/unfolding/
cd ../..

# 2018
mkdir -p fig2018
mkdir -p fig2018/unfolding
cd TUnfoldResults_2018/Nominal
cp --parents ./*/Plots/*.pdf ../../fig2018/unfolding/
cd ../..

# 2016
mkdir -p fig2016
mkdir -p fig2016/unfolding
cd TUnfoldResults_2016/Nominal
cp --parents ./*/Plots/*.pdf ../../fig2016/unfolding/
cd ../..

# fullRun2UL
mkdir -p figfullRun2UL
mkdir -p figfullRun2UL/unfolding
cd TUnfoldResults_fullRun2/Nominal
cp --parents ./*/Plots/*.pdf ../../figfullRun2UL/unfolding/
cd ../..

root -b -q macros/doResPlots_savefits.C
mkdir -p figfullRun2UL/resolution
cp --parents ./UnfoldingHistosResolutionPlots/*.pdf figfullRun2UL/resolution/

# scp to lxplus
scp -r fig* jthieman@lxplus.cern.ch:AN-20-037/.
