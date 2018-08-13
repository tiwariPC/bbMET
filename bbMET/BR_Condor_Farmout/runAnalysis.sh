#!/bin/sh
#### FRAMEWORK SANDBOX SETUP ####
# Load cmssw_setup function
source ./cmssw_setup.sh

# Setup CMSSW Base
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

# Download sandbox
wget --no-check-certificate "http://stash.osgconnect.net/+ptiwari/sandbox-CMSSW_9_4_6_patch1-877018e.tar.bz2"

# Setup framework from sandbox
cmssw_setup sandbox-CMSSW_9_4_6_patch1-877018e.tar.bz2

cd $CMSSW_BASE
cmsenv
cd ../../

python bbMETBranchReader.py -a -F -i "$1" -D . -o BROutput.root --deepcsv

exitcode=$?

if [ ! -e "BROutput.root" ]; then
  echo "Error: The python script failed, could not create the output file."

fi
exit $exitcode
