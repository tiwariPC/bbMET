mkdir -p Systematics/ROOTFiles
cd Systematics
cp ../19062018/bbMETROOT/* ./ROOTFiles
python plot_systematics.py
