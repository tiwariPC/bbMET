from ROOT import TGraph, TFile, TGraphAsymmErrors
from array import array
import os
import matplotlib.pyplot as plt

scanlist=["tanB","sinp"]
for scanname in scanlist:
    os.system("mkdir -p "+scanname)
    scale_dict={50:0.0054626565,100:0.0201360544,350:0.0915618061,400:0.0987681779,500:0.1530311305}
    xsec_dict={}
    xsec_dict={}
    xsecfile=open("scan_run_"+scanname+".txt",'r')
    scanvarlist=[]
    Mlist=[]
    for line in xsecfile:
        if line.startswith("#"): continue
        splt=line.split()
        scanvar=float(splt[1])     #This is tan(B)
        MH4=float(splt[2])
        xsec=float(splt[8])*scale_dict[int(MH4)]
        if not scanvar in scanvarlist: scanvarlist.append(scanvar)
        if not int(MH4) in Mlist: Mlist.append(int(MH4))
        xsec_dict[int(MH4),scanvar]=xsec

    xsecfile.close()
    #print xsec_dict[100,10]
    #for i in scanvarlist: print xsec_dict[50,i]
    
    if scanname=="tanB":
        varname=r'tan$\beta$'
        tit=r"bb+DM 2HDM+a: sin$\theta$ = 0.707"
        ticks=[5*i for i in range(11)]
    elif scanname=="sinp":
        varname=r'sin$\theta$'
        tit=r"bb+DM 2HDM+a: tan$\beta$ = 1"
        ticks=[.1*i for i in range(11)]
    
    plt.clf()
    for M in sorted(Mlist):
        if M==400: continue
        xsecs=[]
        for scanvar in sorted(scanvarlist):
            xsecs.append(xsec_dict[M,scanvar])
        plt.plot(scanvarlist,xsecs,'o-',label=r"$M_{H4}$ = "+str(M)+" GeV")
    
    plt.xlabel(varname)
    plt.ylabel("Cross section (pb)")
    plt.legend()
    plt.title(tit)
    #plt.xscale('log')
    plt.yscale('log')
    #plt.xlim(.8,55)       
    plt.xticks(ticks,ticks)
    plt.grid()
    plt.savefig("XSecs_"+scanname+".png")

    plt.clf()
    for scanvar in sorted(scanvarlist):
        xsecs=[]
        for M in sorted(Mlist):
            if M==400: continue
            xsecs.append(xsec_dict[M,scanvar])
        plt.plot([50,100,350,500],xsecs,'o-',label=str(scanvar))
    plt.xlabel(r'$M_{H4}$')
    plt.ylabel("Cross section (pb)")
    plt.legend(ncol=3,title=varname)
    plt.title(tit)
    plt.xscale('log')
    plt.yscale('log')
    #plt.xlim(.8,55)
    ticks=[50,100,350,500]
    plt.xticks(ticks,ticks)
    plt.grid()
    plt.savefig("XSecs_"+scanname+"_M.png")

    log=""
        
    filename = '../../bin/limits_bbDM20162HDM.txt'
    f = open(filename,"r") 

    for line in f:
        xvals=array('f')
        expm2=array('f')
        expm1=array('f')
        expmed=array('f')
        expp1=array('f')
        expp2=array('f')
        obs=array('f')
        errx=array('f')

        M=int(line.rstrip().split()[0])
        if M==400: continue
        log+="MH4="+str(M)+":\n"
        log+="#"+scanname+" expm2 expm1 exp expp1 expp2 obs\n"
        for ival in scanvarlist:
            xsec_ = float(xsec_dict[M,ival])
            
            ival_=float(ival)
            expm2_=float(line.rstrip().split()[2])/xsec_
            expm1_=float(line.rstrip().split()[3])/xsec_
            expmed_=float(line.rstrip().split()[4])/xsec_
            expp1_=float(line.rstrip().split()[5])/xsec_
            expp2_=float(line.rstrip().split()[6])/xsec_
            obs_=float(line.rstrip().split()[7])/xsec_
            xvals.append(ival_)        
            expm2.append(expmed_-expm2_)
            expm1.append(expmed_-expm1_)
            expmed.append(expmed_)
            expp1.append(expp1_-expmed_)
            expp2.append(expp2_-expmed_)
            obs.append(obs_)
            errx.append(0.0)
            
            log+=str(ival_)+" "+str(expm2_)+" "+str(expm1_)+" "+str(expmed_)+" "+str(expp1_)+" "+str(expp2_)+" "+str(obs_)+"\n"
    #        print ival, float(line.rstrip().split()[4]), xsec_, float(line.rstrip().split()[4])/xsec_
        log+="\n"
        
        g_exp2  = TGraphAsymmErrors(int(len(xvals)), xvals, expmed, errx, errx, expm2, expp2 )   ;  g_exp2.SetName("exp2")
        g_exp1  = TGraphAsymmErrors(int(len(xvals)), xvals, expmed, errx, errx, expm1, expp1 )   ;  g_exp1.SetName("exp1")
        g_expmed = TGraphAsymmErrors(int(len(xvals)), xvals, expmed)   ;  g_expmed.SetName("expmed")
        g_obs    = TGraphAsymmErrors(int(len(xvals)), xvals, obs   )   ;  g_obs.SetName("obs")

        f1 = TFile(scanname+'/'+scanname+'_MH4-'+str(M)+'.root','RECREATE')
        g_exp2.Write()
        g_exp1.Write()
        g_expmed.Write()
        g_obs.Write()

        f1.Write()
        f1.Close()
    f.close()
    logfile=open("alllimits_"+scanname+".txt","w")
    logfile.write(log)
    logfile.close()


