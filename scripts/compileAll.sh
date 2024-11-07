#!/bin/sh

########################################################################
############################# STEERING #################################
########################################################################
#
# ************** cmake version **************
# use older cmake
cmake="cmake"
# ... or newer cmake which supports more options
#cmake="cmake3"
#
# ********** extra options for make *********
# to suppress messages make[2]: Entering directory make[2]: Leaving directory
#makeExtraOptions="--no-print-directory"
# (cmake extra options can be specified via command line, call this script with -h)
#
# ****** number of jobs for make/scram ******
Nsplit=8
#
# ******** scram build extra options ********
#scramExtraOptions="USER_CXXFLAGS+=\"-g\" USER_CXXFLAGS+=\"-O0\""
#scramExtraOptions="USER_CXXFLAGS+=\"-DNDEBUG\""
#
########################################################################
########################## END OF STEERING #############################
########################################################################

if [ ! ${CMSSW_BASE} ] ; then
    echo "Environment variable CMSSW_BASE not set, cannot continue!"
    echo "First do the cmsenv in your release"
    exit 6
fi


TopAnalysisDir=${CMSSW_BASE}/src/TopAnalysis
analysisDir=${TopAnalysisDir}/Configuration/analysis


ztopDir=${TopAnalysisDir}/ZTopUtils
commonSourceDir=${analysisDir}/common
mainSourceDir=${analysisDir}/diLeptonic


if [ $# == 1 ] && ( [ "$1" == "-h" ] || [ "$1" == "--help" ] ); then
    echo "Usage (default installation path): $0"
    echo "Usage (installation in path <folder>): $0 <folder>"
    echo "Usage (installation in path <folder> with additional cmake options <options>): $0 <folder> <options>"
    exit 1
elif [ $# -ge 1 ]; then
    installDir=$1
    if [ ! -d "${installDir}" ]; then
        echo "Specified install directory not existing: ${installDir}"
	echo "Please create directory before installation"
	exit 2
    fi
    cd ${installDir}
    commonInstallDir=${PWD}
    mainInstallDir=${PWD}
    cd -
    # all other arguments )if any) will be treated as extra options for cmake
    if [ $# -ge 2 ] ; then
        cmakeExtraOptions=${@:2}
        echo "cmake will be used with extra options: ${cmakeExtraOptions}"
    fi
elif [ $# == 0 ]; then
    commonInstallDir=${commonSourceDir}
    mainInstallDir=${mainSourceDir}
fi

if [ -z ${makeExtraOptions} ]; then
  echo "make will be used with extra options: ${makeExtraOptions}"
fi

echo
echo
echo "Compiling project ZTopUtils using scram"
echo "Source code in folder: ${ztopDir}"
echo "Installation in folder: as defined by scram"
echo "Warning about 'Invalid tool lhapdffull' can be ignored"
echo
cd ${ztopDir}
scram b -r ${scramExtraOptions} -j${Nsplit}
if [ $? -eq 0 ] ; then
    echo
    echo "Compilation successful"
else
    echo
    echo "Compilation NOT successful, stopping..."
    exit 3
fi
cd -
echo
echo


echo
echo
echo "Compiling project ntupleCommon using cmake/make"
echo "Source code in folder: ${commonSourceDir}"
echo "Installation in folder: ${commonInstallDir}"
echo
if [ ! -d "${commonInstallDir}/build_ntupleCommon" ]; then
    mkdir ${commonInstallDir}/build_ntupleCommon
fi
cd ${commonInstallDir}/build_ntupleCommon
${cmake} -D CMAKE_INSTALL_PREFIX=${commonInstallDir}/install ${cmakeExtraOptions} ${commonSourceDir}
make -j${Nsplit} ${makeExtraOptions} install
if [ $? -eq 0 ] ; then
    echo
    echo "Compilation successful"
else
    echo
    echo "Compilation NOT successful, stopping..."
    exit 4
fi
cd -
echo
echo


echo
echo
echo "Compiling project diLeptonic using cmake/make"
echo "Source code in folder: ${mainSourceDir}"
echo "Installation in folder: ${mainInstallDir}"
echo
if [ ! -d "${mainInstallDir}/build_diLeptonic" ]; then
    mkdir ${mainInstallDir}/build_diLeptonic
fi
cd ${mainInstallDir}/build_diLeptonic
${cmake} -D CMAKE_INSTALL_PREFIX=${mainInstallDir}/install ${cmakeExtraOptions} ${mainSourceDir}
make -j${Nsplit} ${makeExtraOptions} install
if [ $? -eq 0 ] ; then
    echo
    echo "Compilation successful"
else
    echo
    echo "Compilation NOT successful, stopping..."
    exit 5
fi
cd -
echo
echo






