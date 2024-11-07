#!/bin/sh
source /cvmfs/cms.cern.ch/cmsset_default.sh
set SCRAM_ARCH=slc6_amd64_gcc530
export SCRAM_ARCH

if [ -z "$LD_LIBRARY_PATH" ]; then
LD_LIBRARY_PATH=.
export LD_LIBRARY_PATH
fi


SCRIPT=$(readlink -f $0)
SCRIPTPATH=$(dirname "$SCRIPT")
#BIN="@CMAKE_INSTALL_PREFIX@/bin"
BIN="/depot/cms/private/users/akhatiwa/TTBarSpinCor_fullInstallation/CMSSW_8_0_26_patch1/src/TopAnalysis/Configuration/analysis/diLeptonic/install/bin"

ARGS=$1
echo ARG: "$ARG"
echo BIN: "$BIN/load_Analysis"
cd /depot/cms/private/users/akhatiwa/TTBarSpinCor_fullInstallation/CMSSW_8_0_26_patch1/src/TopAnalysis/Configuration/analysis/diLeptonic
pwd
echo "cmsenv"
cmsenv
echo "$BIN/load_Analysis $ARG"
$BIN/load_Analysis $ARG
