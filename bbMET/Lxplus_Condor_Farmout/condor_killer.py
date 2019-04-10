import os,sys
if len(sys.argv)!=4:
	print "Usage: python condor_killer.py clusterid_min clusterid_max username"
	sys.exit()
os.system("condor_rm -constraint 'ClusterId >="+sys.argv[1]+" && ClusterId >="+sys.argv[2]+"' "+sys.argv[3])
