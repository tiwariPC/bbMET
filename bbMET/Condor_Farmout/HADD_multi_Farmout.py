import os
folder="BR_Outputs"

flist = ['_'.join(f.split('.')[0].split('_')[:-2]) for f in os.listdir(folder) if f.endswith('.root')]

for fl in list(set(flist)):
    filestr = ""
    filestr=folder+'/'+fl+"_*.root"
    os.system("mkdir -p hadd_outputs/")
    command = "hadd hadd_outputs/"+fl+".root "+filestr
    os.system(command)
#        print command
