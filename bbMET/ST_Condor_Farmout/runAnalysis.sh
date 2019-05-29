#!/bin/sh
#### FRAMEWORK SANDBOX SETUP ####
# Load cmssw_setup function
source ./cmssw_setup.sh

# Setup CMSSW Base
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

# Download sandbox
#wget --no-check-certificate "http://stash.osgconnect.net/+ptiwari/sandbox-CMSSW_8_0_26_patch1-76efecd.tar.bz2"

# Setup framework from sandbox
cmssw_setup sandbox-CMSSW_8_0_26_patch1-76efecd.tar.bz2

cd $CMSSW_BASE
cmsenv
cd ../../

export X509_USER_PROXY=$1
voms-proxy-info -all
voms-proxy-info -all -file $1

python SkimTree.py -F -i "$2"
until xrdcp -f SkimmedTree.root root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/bbMET/2016_skimmer/bbDM_27052019_miniAODv2_ST/"$3"/SkimmedTree_"$4".root; do
  sleep 60
done

exitcode=$?

if [ ! -e "SkimmedTree.root" ]; then
  echo "Error: The python script failed, could not create the output file."

fi
exit $exitcode
