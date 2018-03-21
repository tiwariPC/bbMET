# bbMET

# I. Run DMAnalyzer
First step is to generate ntuples.
1. Follow instructions from https://github.com/tiwariPC/DMAnaRun2/tree/onlyAK4_80X_puppi_deepCSV (onlyAK4_80X_puppi_deepCSV branch) to setup DelPanj within CMSSW.
2. Use MultiCrab to submit jobs for signal, backgrounds, and data separately. (Detailed instructions to be updated soon).

# II. Run SkimTree
Clone this repository in a location from where HTCondor jobs can be submitted (```login.uscms.org``` for example). We shall refer to this location as the working directory for the rest of this documentation.
## a. Running Locally
SkimTree can be run using:
```bash
python SkimTree.py path_to_ntuple_file
```
##### Example:
```bash
python SkimTree.py root://se01.indiacms.res.in:1094//dpm/indiacms.res.in/home/cms/store/user/zabai/t3store2/bbDM_bkg/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_LegacyMC_20170328/180202_114154/0001/NCUGlobalTuples_1009.root
```
## b. Running on HTCondor
***Note: The Condor framework is configured to run BranchReader whenever one runs SkimTree, since this helps save a lot of time. So one will get the final BranchReader outputs directly after running this step besides getting the Skimmed Trees.***
#### Making a list of ntuple files
The T2ListMaker can be used to produce lists of ntuples (CRAB job outputs) that SkimTree takes as input. Since the CRAB job outputs are stored in T2, the ListMaker can be used from a T3 location which has read access to the T2 location (via `rfdir`).
1. Copy the T2ListMaker directory from your working directory to the T3 location (if different from your working directory), or run the following from a T3 location:
    ```bash
    wget https://raw.githubusercontent.com/mondalspandan/bbMET/master/T2FileListMaker/ListMaker_T3.py
    ```
2. Find the parent directory in T2 where the ntuples are stored.
    Example:
    ```bash
    rfdir /dpm/indiacms.res.in/home/cms/store/user/username/t3store2
    ```
    can be used to view contents of the directory. Navigate to the parent directory where the ntuples are stored.
3. Run `ListMaker_T3.py` using
    ```bash
    python ListMaker_T3.py crab path_to_parent_T2_directory_of_ntuples filelist_tag
    ```
    Example:
    ```bash
    python ListMaker_T3.py crab /dpm/indiacms.res.in/home/cms/store/user/spmondal/t3store2/bbDM_data 180321_data_ntuples
    ```
4. Repeat the above for each of signal, data and bkg. So in the end you will probably have 3 or more different directories (and corresponding log files).
This produces one directory and one log file. In the example above, the directory and logfile will be named `Filelist_180321_data_ntuples` and `log_180321_data_ntuples.txt` respectively. The filelists are now available in your T3 location. Both the directory and the log files now have to be copied back to the working directory.
5. Navigate to `bbMET/bbMET/ST_BR_Condor/` in the working directory.
6. `mkdir Filelists`
7. Now copy all directories and log files (such as `Filelist_180321_data_ntuples` and `log_180321_data_ntuples.txt`) to workingdir/`bbMET/bbMET/ST_BR_Condor/Filelists` using `cp` (or `scp` if the working directory is not in T3).

#### Running SkimTree + BranchReader Condor Jobs
1. Navigate to workingdir/`bbMET/bbMET/ST_BR_Condor/`.
2. Open `runAnalysis.sh` to edit.
3. Line 19 contains an exemplar path to a T2 location. The SkimTree outputs (Skimmed Trees) will be stored in this location. Edit this path (remember to change the username) and specify where you wish to store the outputs. The directory does not necessarily have to exist, it will be created if non-existent. Make sure you have write access to the specified directory.
4. Initiate your voms-proxy using `voms-proxy-init --voms cms --valid 192:00`.
5. Submit jobs using
    ```bash
    . submitjobs.sh
    ```
6. To monitor the job submission process, use:
    ```bash
    tail -f logsubmit.txt
    ```
7. To monitor status of jobs, use `condor_q username`.

#### Retrieving Outputs
* The SkimTree outputs (Skimmed Trees) are saved to the location that was specified in the `runAnalysis.sh` file. These can be used while running BranchReader (next section) alone in future iterations.
* The Condor jobs that were submitted for SkimTree run BranchReader as well as soon as the Skimming part is over. Hence, it transfers the final BranchReader outputs to the working directory as well. All these will be stored in the `bbMET/bbMET/ST_BR_Condor/data` directory, but one needs to combine several fragments of each sample before one can use them in plotting.

**(SkimTree+BranchReader jobs may take upto 2 days to finish.)**
1. Once all SkimTree+BranchReader Condor jobs are complete, one needs to combine the output .root files for each sample. This can be achieved by using the `hadd` command. If no CMSSW or ROOT instance is sourced by default in your working area, go to a CMSSW release base and run `cmsenv`, otherwise the `hadd` command may not work.
2. Navigate to `bbMET/bbMET/ST_BR_Condor/` and run
    ```bash
    python HADD_multi_v2.py
    ```
    The outputs .root files are stitched on a per sample basis and one .root file per sample is produced inside `bbMET/bbMET/ST_BR_Condor/BROutputs`/Filelist_tag (where Filelist_tag would be `Filelist_180321_data_ntuples` according the previous example.
3. These .root files inside `BROutputs` can be used directly as inputs to the plotting code.

# III. Run BranchReader

# IV. Run Plotting Script
