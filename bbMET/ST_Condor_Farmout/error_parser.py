import os, subprocess
folder = "error"
filelist = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]

nerr=0
ntot=0

for ifile in filelist:
    if ntot%1000==0: print "Parsed "+str(ntot)+" files."
    line = subprocess.check_output(['tail', '-2', os.path.join(folder,ifile)])
    ntot+=1
    #print"****"+ line.split('\n')[0].strip()+"****"
    if not "100%" in  line.split('\n')[0].strip(): # != "TTree::Bronch:0: RuntimeWarning: Using split mode on a class: TLorentzVector with a custom Streamer":
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
        
