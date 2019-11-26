
# PURPOSE: returns the RA and Dec coordinates of the center of the cluster

def getrd(locclu) :

	
	# read rac and decc
	racfile = open(locclu+'/rac.txt', 'r')
	rac= racfile.readline()
	rac=rac[0:-1]
	racfile.close()
	deccfile = open(locclu+'/decc.txt', 'r')
	decc= deccfile.readline()
	decc=decc[0:-1]
	deccfile.close()


	return rac, decc
