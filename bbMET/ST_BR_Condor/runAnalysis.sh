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
python SkimTree.py "$1"
until xrdcp -f SkimmedTree.root root://se01.indiacms.res.in//dpm/indiacms.res.in/home/cms/store/user/spmondal/t3store2/bbDM_SkimmedTrees_test/"$2"/SkimmedTree_"$3".root; do
  sleep 60
done
python bbMETBranchReader.py -a -i SkimmedTree.root -D . --csv --met

exitcode=$?

if [ ! -e "Output_SkimmedTree.root" ]; then
  echo "Error: The python script failed, could not create the output file."

fi
exit $exitcode

