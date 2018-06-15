import os,sys
if len(sys.argv)!=2 and len(sys.argv)!=3:
	print "Usage: python condor_killer.py clusterid_min [clusterid_max]"
	sys.exit()

os.system("condor_q spmondal &> htcq.txt")
htcq=open('htcq.txt','r')
joblist=""
for line in htcq:
	try:
		cid=float(line.split()[0])
		if len(sys.argv)==2:
			willKill=cid>=float(sys.argv[1])
		else:
			willKill=cid>=float(sys.argv[1]) and cid<float(sys.argv[2])
		if willKill: joblist+=line.split()[0]+" "
	except:
		pass
	if len(joblist)>=31900:
		print "Executing condor_rm..."
		os.system("condor_rm  "+joblist)
		joblist=""
htcq.close()
#print len(joblist)
os.system("condor_rm  "+joblist)

