import sys 
import os 

datacard=sys.argv[1]

med=datacard.split("_")[-3]

mchi=datacard.split("_")[-1].split(".")[0]

model=datacard.split("_")[-5]


logfile = "log_"+med+"_"+mchi+".del"
command_  = "combine -M Asymptotic --rAbsAcc 0 --rMax 30 -t -1 "+ datacard + " >> " + logfile

os.system(command_)
expected25_="" 
expected16_="" 
expected50_="" 
expected84_="" 
expected975_=""
observed_=""
for ilongline in open(logfile):
    if "Observed Limit: r < " in ilongline:
        observed_ = ilongline.replace("Observed Limit: r < ","").rstrip()
    if "Expected  2.5%: r < " in ilongline:
        expected25_ = ilongline.replace("Expected  2.5%: r < ","").rstrip()
    if "Expected 16.0%: r < " in ilongline:
        expected16_ = ilongline.replace("Expected 16.0%: r < ","").rstrip()
    if "Expected 50.0%: r < " in ilongline:
        expected50_ = ilongline.replace("Expected 50.0%: r < ","").rstrip()
    if "Expected 84.0%: r < " in ilongline:
        expected84_ = ilongline.replace("Expected 84.0%: r < ","").rstrip()
    if "Expected 97.5%: r < " in ilongline:
        expected975_ = ilongline.replace("Expected 97.5%: r < ","").rstrip()

towrite =  med+" "+mchi+" "+expected25_+" "+expected16_+" "+ expected50_+" "+ expected84_+" "+ expected975_+" "+ observed_+"\n"

print towrite
outfile = 'bin/limits_bbDM2016'+model+'.txt'
fout = open(outfile,'a')
fout.write(towrite)


