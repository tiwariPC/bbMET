import os
folder="data"

listfolders = [f for f in os.listdir(folder) if not os.path.isfile(os.path.join(folder, f))]

for fold1 in listfolders:
    fold2list = [f for f in os.listdir(os.path.join(folder, fold1)) if not os.path.isfile(os.path.join(folder, fold1, f))]
    for fold2 in fold2list:
        filelist = sorted([f for f in os.listdir(os.path.join(folder, fold1, fold2)) if os.path.isfile(os.path.join(folder, fold1, fold2, f))])
        filestr = ""
        for filename in filelist:
            filestr += os.path.join(folder, fold1, fold2)+"/"+filename + " "
        os.system("mkdir -p BROutputs/"+fold1)
        command = "hadd BROutputs/"+fold1+"/Output_"+fold2+".root "+filestr
        os.system(command)
#        print command+"\n"
