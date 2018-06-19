# bbMET
### Workflow:

![alt text](http://spmondal.web.cern.ch/spmondal/bbMETFlow.png)


# 1. Run DMAnalyzer

First step is to generate ntuples.

1. Follow instructions from https://github.com/tiwariPC/DMAnaRun2/tree/onlyAK4_80X_puppi_deepCSV (onlyAK4_80X_puppi_deepCSV branch) to setup DelPanj within CMSSW.
2. Use MultiCrab to submit jobs for signal, backgrounds, and data separately. (Detailed instructions to be updated soon).

# 2. Run SkimTree

Clone this repository in a location from where HTCondor jobs can be submitted (```login.uscms.org``` for example). We shall refer to this location as the working directory for this section.

## 2.1. Running Locally

SkimTree can be run using:
```bash
python SkimTree.py path_to_ntuple_file
```
##### Example:
```bash
python SkimTree.py root://se01.indiacms.res.in:1094//dpm/indiacms.res.in/home/cms/store/user/zabai/t3store2/bbDM_bkg/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_LegacyMC_20170328/180202_114154/0001/NCUGlobalTuples_1009.root
```
## 2.2. Running on HTCondor

***SkimTree jobs are numerous is number, hence using a Condor network with large number of cores is recommended. Recommended submit node: uscms.***

### 2.2.1. Making a list of ntuple files

The T2ListMaker can be used to produce lists of ntuples (CRAB job outputs) that SkimTree takes as input. Since the CRAB job outputs are stored in T2, the ListMaker can be used from a T3 location which has read access to the T2 location (via `rfdir`).
1. Copy the T2ListMaker directory from your working directory to the T3 location (if different from your working directory), *or* run the following from a T3 location:
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
This produces one directory and one log file. In the example above, the directory and logfile will be named `Filelist_180321_data_ntuples` and `log_180321_data_ntuples.txt` respectively. The filelists are now available in your T3 location. Both the directory and the log files now have to be copied back to the working directory.
4. Repeat the above for each of signal, data and bkg. So in the end you will probably have 3 or more different directories (and corresponding log files).
5. Navigate to `bbMET/bbMET/ST_Condor/` in the working directory.
6. `mkdir Filelists`
7. Now copy all directories and log files (such as `Filelist_180321_data_ntuples` and `log_180321_data_ntuples.txt`) to workingdir/`bbMET/bbMET/ST_Condor/Filelists` using `cp` (or `scp` if the working directory is not in T3).

### 2.2.2. Running SkimTree Condor Jobs

1. Navigate to workingdir/`bbMET/bbMET/ST_Condor/`.
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

### 2.2.3. Outputs

The SkimTree outputs (Skimmed Trees) are saved to the location that was specified in the `runAnalysis.sh` file. These will be used while running BranchReader (next section).

**(SkimTree jobs for all data and background may take upto 1 day to finish.)**

# 3. Run BranchReader

Clone this repository in a location from where HTCondor jobs can be submitted (```ui.indiacms.res.in``` for example). We shall refer to this location as the working directory for this section.

## 3.1. Running Locally

* BranchReader can be run using:
```bash
python bbMETBranchReader.py -a -i path_to_skimmed_tree -D . --csv
```
* The "farmout" mode can be used to combine multiple input root files in one go. If all the root files are listed in a text file named input.txt, one can use:
```bash
python bbMETBranchReader.py -a -F -i input.txt -D . --csv
```
## 3.2. Running on HTCondor

***The number of BranchReader jobs can be adjusted by combining suitable number of input files in the Farmout mode. Recommended submit node: ui.indiacms.***

### 3.2.1. Making a list of skimmed trees

The T2ListMaker can be used to produce lists of skimmed trees (SkimTree outputs) that SkimTree takes as input. Since the skimmed trees are also stored in T2, the ListMaker can be used from a T3 location which has read access to the T2 location (via `rfdir`). For Farmout mode of the BranchReader, it is necessary to use the "nolog" variant of the ListMaker.

1. Copy the T2ListMaker directory from your working directory to the T3 location (if different from your working directory), *or* run the following from a T3 location:
    ```bash
    wget https://raw.githubusercontent.com/mondalspandan/bbMET/master/T2FileListMaker/ListMaker_T3_nolog.py
    ```
2. Find the parent directory in T2 where the skimmed trees are stored.
    Example:
    ```bash
    rfdir /dpm/indiacms.res.in/home/cms/store/user/username/t3store2
    ```
    can be used to view contents of the directory. Navigate to the parent directory where the skimmed trees are stored.
3. Run `ListMaker_T3_nolog.py` using
    ```bash
    python ListMaker_T3_nolog.py st path_to_parent_T2_directory_of_skimmed_trees filelist_tag
    ```
    Example:
    ```bash
    python ListMaker_T3_nolog.py st /dpm/indiacms.res.in/home/cms/store/user/spmondal/t3store2/bbDM_skimmed_trees 180322_skimmed_trees
    ```
This produces one directory. In the example above, the directory will be named `Filelist_180322_skimmed_trees`, The filelists are now available in your T3 location. Both the *contents of the directory* now have to be copied back to the working directory.
4. Repeat the above for each of signal, data and bkg. So in the end you will probably have 3 or more different directories).
5. Navigate to `bbMET/bbMET/BR_Condor_Farmout/` in the working directory.
6. `mkdir Filelists`
7. Now copy *files inside all the directories* (such as `Filelist_180322_skimmed_trees/*.txt) to workingdir/`bbMET/bbMET/BR_Condor_Farmout/Filelists` using `cp` (or `scp` if the working directory is not in T3). This means there will be .txt files inside the Filelists folder, unlike in case of SkimTree jobs.

### 3.2.2. Running BranchReader Condor Jobs

1. Navigate to workingdir/`bbMET/bbMET/BR_Condor_Farmout/`.
2. This framework automatically combines multiple skimmed_tree root files in one job. The number of root files to be combined in each job can be specified by editing L4 of MultiSubmit.py (`maxfilesperjob=100`).
3. Initiate your voms-proxy using `voms-proxy-init --voms cms --valid 192:00`.
4. Submit jobs using
    ```bash
    . submitjobs.sh
    ```
5. To monitor the job submission process, use:
    ```bash
    tail -f logsubmit.txt
    ```
6. To monitor status of jobs, use `condor_q username`.

### 3.2.3. Retrieving Outputs

1. Once all BranchReader Condor jobs are complete, one needs to combine the output .root files for each sample. This can be achieved by using the `hadd` command. If no CMSSW or ROOT instance is sourced by default in your working area, go to a CMSSW release base and run `cmsenv`, otherwise the `hadd` command may not work.
2. Navigate to `bbMET/bbMET/BR_Condor_Farmout/` and run
    ```bash
    python HADD_multi_Farmout.py
    ```
    The outputs .root files are stitched on a per sample basis and one .root file per sample is produced inside `bbMET/bbMET/BR_Condor_Farmout/hadd_outputs`.
3. These .root files inside `hadd_outputs` can be used directly as inputs to the plotting code and will henceforth be referred to as **BranchReader outputs**.ch

# 4. Generate Signal Efficiency Plots

***This section will require BranchReader outputs of only signal samples. Can be any or all of bbDM-NLO, bbDM-LO or ttDM samples.***

### 4.1. 2D Efficiency Plots

1. Navigate to `SignalEffPlots` directory and then `cd` to either NLO, LO or ttDM depending on which you wish to plot.
2. `mkdir signal`
3. Copy **all** BranchReader outputs of this particular signal to this `signal` directory.
4. Run `python SEffPlot.py` to produce the plots.

The required file name format for files inside the signal directory are as follows:

* NLO: Output_pseudo_NLO_Mchi-x_Mphi-y.root
* LO: Output_pseudo_LO_Mchi-x_Mphi-y.root
* ttDM: Output_pseudoscalar_Mchi-x_Mphi-y_TuneZZZ.root

### 4.2. 1D Plots

1D plotting script is manual and only available for NLO.
1. Edit SignalEffPlots/NLO/1D.py L4-5 to set the values.
2. Run using `python 1D.py`

### 4.3. Cut-by-cut efficiencies

The idea is to see the signal efficiencies after each cut is applied, as a function of M_Med for fixed M_DM.
1. Navigate to `CutbycutEff` directory and `mkdir signal`.
2. Copy all signal files with M_DM = 1 GeV (Output_pseudo_NLO_Mchi-1_Mphi-*.root) to `signal` directory.
3. Run `python CutByCutEff.py`.

# 5. Run Plotting Script

1. Copy all the outputs from BranchReader (in `hadd_outputs` directory) to a directory (or, optionally, segregate the files in separate directories named `data`, `bkg`, and `signal`).
2. Open `bbMETplot/Scripts/bbMET_StackFactory.py` and edit L92 to suit the current working directory. Edit L160, L229, and L272 to the path(s) where the BranchReader outputs are stored.
3. Navigate to `bbMETplot/Scripts/test` and run
```
python ../bbMET_StackFactory.py -d MET -s -m -q
python ../bbMET_StackFactory.py -d SE -e
python ../bbMET_StackFactory.py -d SP -p
```
The boolean flags are explained as follows:

* s: Signal Regions (only MC plots)
* m: Muon Control Regions
* q: QCD Control Regions
* e: Electron Control Regions
* p: Photon Control Regions
* The `-d` flag is used to select the appropriate primary dataset for each region.
4. Similarly for getting systematics only, navigate to `bbMETplot/Scripts/syst` and run:
```
python ../bbMET_StackFactory_syst.py -d MET -s -m -q
python ../bbMET_StackFactory_syst.py -d SE -e
python ../bbMET_StackFactory_syst.py -d SP -p
```
The root files are stored inside syst/date/bbMETROOT directory.
Go to syst directory and run the following file(update the directory named by date):
```
. syst_plot.sh
```

### 5.1. CR Summary plots

* Control region summary plots, separately for muon and electron regions, are automatically plotted by the above plotting script. However, the combined summary plot, although produced by the plotting script, will not be correct.
* To make a combined summary plots, the number of data events and bkgsum events for each region have to be copied to bbMETplot/Scripts/CRSummary.py, L5-L6. It can then be run using ```python CRSummary.py```.

# 6. Make a combined .root histogram file

The idea is to make a combined .root file containing MET and Hadronic Recoil histograms of all regions with all systematics with suitable names.

1. Open bbMETplot/Scripts/CombinedRootMaker.py and edit L36 if histograms for systematics are contained in a separate directory. Otherwise point this to read the same directory as next step.
2. Run:
  ```
  python CombinedRootMaker.py path_to_the_plot_script_output_dir path_to_systematics
  ```
  Example:
  ```
  python CombinedRootMaker.py test/22022018 syst/22022018
  ```

This creates a directory named DataCardRootFiles. The file `AllMETHistos.root` inside this directory contains all histograms from all regions with all systematics. Besides this, a .root file for each region is also created. Depending on the signal model either the `AllMETHistos.root` file or all the other files need to be used as input to the limit setting code.

# 7. Preliminary Limits and Fitting

1. Copy all files from bbMETplot/Scripts/DataCardRootFiles directory to LimitsAndFitting/bbDM/data directory.
2. Run:
   ```
   python RunLimitOnAll.py create
   python RunLimitOnAll.py run
   ```
   This will create .txt files inside the bin/ directory.
   Navigate to the ```plotting``` directory.
3. To make DMSimp plots:
   ```
   python TextToTGraph.py
   python plotLimit.py ps
   python plotLimit.py s
   ```
4. To make 2HDM+a plots with variable MH4:
   ```
   python TextToTGraph_2HDMa.py
   plotLimit.py 2h
   ```
5. To make 2HDM+a plots with fixed MH4 (50 or 100 GeV) but variable tanβ or variable sinθ:
   ```
   python TextToTGraph_fixedMH4.py
   python plotLimit.py tanb50
   python plotLimit.py tanb100
   python plotLimit.py sinp50
   ```
6. To plot overlapping cross section and limit plots for various tanβ and sinθ scans:
   ```
   cd variableTanB
   python TextToGraph_tanB_allMH4.py
   python plotLimit_overlap.py tanb
   python plotLimit_overlap.py sinp
   ```
