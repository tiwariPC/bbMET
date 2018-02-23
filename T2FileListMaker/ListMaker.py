import os,sys

errmsg="Usage:\n$ python ListMaker.py crab inpfilename.txt outfileprefix\nor\n$ python ListMaker.py st inpfilename.txt outfileprefix"

if len(sys.argv)==4:
    if sys.argv[1]=="crab":
        isCrab=True
    elif sys.argv[1]=="st":
        isCrab=False
    else:
        print errmsg
        sys.exit()
    inpfilename=sys.argv[2]
    filepref=sys.argv[3]
else:
    print errmsg
    sys.exit()
    
f=open(inpfilename,"r")

os.system("mkdir -p Filelist"+"_"+filepref)
pref="root://se01.indiacms.res.in:1094/"
filecount=1
lineminus1=""
lineminus2=""
fileopen=False
failedfile=False
log=open("log_"+filepref+".txt","w")

for line in f:
    if not line=="\n":
        fname=line.split()[-1]
    else:
        fname=""
        
    if fname.endswith(".root") and not fileopen:
        folder=pref+lineminus2[:-2]+"/"
#        print lineminus2
        if lineminus2.split("/")[-1].strip()=="failed:": failedfile=True
        if not failedfile:
            if isCrab:
                log.write(lineminus2.split("/")[-3]+"_"+lineminus2.split("/")[-1][:-2]+" "+filepref+str(filecount)+".txt\n")          # For making a list of CRAB job outputs
            else:
                log.write(lineminus2.split("/")[-1][:-2]+" "+filepref+str(filecount)+".txt\n")      # For making a list of SkimmedTrees
            out=open("Filelist"+"_"+filepref+"/"+filepref+str(filecount)+".txt","w")        
            out.write(folder+fname+"\n")
            filecount+=1
        fileopen=True
    elif fname.endswith(".root"):
        if not failedfile: out.write(folder+fname+"\n")
    elif fileopen:
        if not failedfile: out.close()
        fileopen=False
        failedfile=False
    
    lineminus2=lineminus1
    lineminus1=line
    
log.close()
f.close()

print "Created Filelist_%s directory and log_%s.txt." %(filepref,filepref)
