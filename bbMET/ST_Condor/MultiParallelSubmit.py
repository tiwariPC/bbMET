import os
folder="Filelists"

logs = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]
listfolders = [f for f in os.listdir(folder) if not os.path.isfile(os.path.join(folder, f))]
#print logs
#print listfolders

os.system("chmod +x runAnalysis.sh")
for outdirs in ['error','log','output']:
    os.system("mkdir -p "+outdirs)

count=1

for foldname in listfolders:
#    os.system("mkdir -p data/"+foldname)
    run=foldname[9:]
    logfile=open(os.path.join(folder,"log_"+run+'.txt'),'r')
    sample={run:foldname}
    for line in logfile:
        sample[line.split()[1]]=line.split()[0]
    
    for listname in os.listdir(os.path.join(folder,foldname)):        
        inpfilename=os.path.join(folder,foldname,listname)
        outfoldname="data/"+foldname+"/"+sample[listname]
        os.system("mkdir -p "+outfoldname)
        
        submittemp=open("submit_parallel_temp.sub","w")
        submitfile=open("submit_parallel.sub","r")
        for line in submitfile:
            submittemp.write(line)                        
        submitfile.close()       
        
        rootfile=open(inpfilename,"r")
        for line in rootfile:
            submittemp.write("arguments = \""+line.strip()+" "+outfoldname[5:]+" $(Process)\"\n")    
            submittemp.write("queue\n")
        
        rootfile.close()
        submittemp.close()
        
        print "\n===============================\nSubmitting jobs set #"+str(count)+" from "+inpfilename+": "+ sample[listname]+"\n===============================\n"
        
        os.system("condor_submit submit_parallel_temp.sub")
        count+=1
    logfile.close()

print "\nDone. Submitted "+str(count-1)+" jobs."
