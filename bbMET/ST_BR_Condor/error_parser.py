import os, subprocess
folder = "error"
filelist = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]

nerr=0
ntot=0

for ifile in filelist:
    line = subprocess.check_output(['tail', '-2', os.path.join(folder,ifile)])
    ntot+=1
    if line.split('\n')[0].strip() != "TFile::Append:0: RuntimeWarning: Replacing existing TH1: h_total (Potential memory leak).":
        n=2
        while line.strip() == "":
            n+=1
            line = subprocess.check_output(['tail', '-'+str(n), os.path.join(folder,ifile)])
        line = subprocess.check_output(['tail', '-'+str(n+1), os.path.join(folder,ifile)])   
        print ifile+" :"
        print line
        print
        nerr+=1

print "\nParsed "+str(ntot)+" error files. Found "+str(nerr)+" errors."
        
