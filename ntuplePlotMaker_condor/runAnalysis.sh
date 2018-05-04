#### FRAMEWORK SANDBOX SETUP ####
# Load cmssw_setup function
source cmssw_setup.sh

# Setup CMSSW Base
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

# Download sandbox
wget --no-check-certificate "http://stash.osgconnect.net/+spmondal/sandbox-CMSSW_8_0_26_patch1-939efad.tar.bz2"

# Setup framework from sandbox
cmssw_setup sandbox-CMSSW_8_0_26_patch1-939efad.tar.bz2

cd $CMSSW_BASE
cmsenv
cd ../../
python ntuplePlotMaker.py "$1"

exitcode=$?

if [ ! -e "ntupleHistos.root" ]; then
  echo "Error: The python script failed, could not create the output file."

fi
exit $exitcode

