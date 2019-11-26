
# v2: no more spectrum version needed (in hilx/ there used to be verstr='_v2'

# PURPOSE: bins the in or out spectrum in 33 spectral bins.
# OUTPUT: mfe_CLU/spec/binsrcin.txt OR mfe_CLU/spec/binsrcout.txt

import pdb

import numpy as nn

def bin( clu, srctxt ) :


	sfile = open( 'mfe_'+clu+'/spec/'+srctxt, 'r' )



	intermed=sfile.readline()
	ee=[]
	cc=[]
	while intermed != '' and intermed != '\n' :

		if intermed[0] != '#' :
			ee.append( float(intermed.split()[0]) )
			cc.append( float(intermed.split()[1]) )

		intermed=sfile.readline()

	ee=nn.array(ee)
	cc=nn.array(cc)

	nlines=len(ee)

	eelo= nn.linspace(.3, 6.9, 34)
	eehi= eelo+0.2
	eemid=eelo+0.1
	bincc= nn.zeros(34, dtype=nn.float)


	for jj in range(34) : 


		integ_range = nn.where( (ee >= eelo[jj]) & (ee < eehi[jj]) )
		#pdb.set_trace()
		if   len(integ_range[0])> 1 : bincc[jj] = nn.trapz(ee[integ_range], cc[integ_range])
		elif len(integ_range[0])==0 : bincc[jj] = 0.
		elif len(integ_range[0])==1 : bincc[jj] = 0.2 * cc[integ_range]



	# normalize the spectrum:
	normbc = bincc / nn.sum(bincc)


	filebin = open( 'mfe_'+clu+'/spec/bin'+srctxt, 'w')
	for jj in range(34): filebin.write(str(eemid[jj])+' '+str(normbc[jj])+'\n')
	filebin.close()




