
# PURPOSE: create light-curve, clip data for high count rate times, prompt user to 
# PURPOSE: check results

# v2: cleanOBSID.fits: no chip number filtering. include all chips, filter later.
# temporary! : clobber=yes , and on same line indent the dmcopy commnad one more tab

# v3: * use asol3.py to find asol, flt1, pbk, msk, bpix or evt1 files
#     * i dont think 'temporary' statement needs to be addressed.

# v4: * 2017-03-31: new version of asol: asol5_loc

# v5: use the new location of repro3_evtOBSID.fits in mfe_OBSID/

# OUTPUT: cleanOB.fits
# OUTPUT: lc/gtiOB.fits

import os
import glob
import commands
import shutil
import readline

import asol5_loc

from sherpa.astro.ui import *
from lightcurves import *

def lc(clu, obsids) :

	print
	print 'LC for cluster '+clu+' ---------------------------'
	print

	origclu = 'orig_' + clu
	locclu = 'mfe_' + clu

	if not os.path.exists(locclu+'/lc') : os.mkdir(locclu+'/lc')

	racfile = open(locclu+'/rac.txt', 'r')
	rac= racfile.readline()
	rac=rac[0:-1]
	racfile.close()
	deccfile = open(locclu+'/decc.txt', 'r')
	decc= deccfile.readline()
	decc=decc[0:-1]
	deccfile.close()

	shutil.copy(locclu+'/reg/pt0mfe_wcs.reg', locclu+'/reg/exclude.reg')
	reg_file = open(locclu+'/reg/exclude.reg', 'a')
	reg_file.write('circle('+rac+','+decc+',305)')
	reg_file.close()

	for i in range(len(obsids)) :

		if not os.path.exists(locclu+'/lc/gti'+obsids[i]+'.fits') :

			os.system('punlearn ardlib')
			os.system('acis_set_ardlib ' + asol5_loc.asol(obsids[i],'bpix') )

			os.system('punlearn dmcopy')
			os.system('dmcopy "mfe_'+clu+'/repro3_evt'+obsids[i]+'.fits[exclude sky=region('+locclu+'/reg/exclude.reg)]" '+locclu+'/nosrc'+obsids[i]+'.fits')

			os.system('punlearn dmextract')
			os.system('dmextract infile="'+locclu+'/nosrc'+obsids[i]+'.fits[bin time=::259.28]" outfile='+locclu+'/lc/lc'+obsids[i]+'.fits opt=ltc1 mode=h')

			try    : lc_clean(locclu+'/lc/lc'+obsids[i]+'.fits', outfile=locclu+'/lc/gti'+obsids[i]+'.fits' , mean=None , clip=3.0 , scale=1.2, minfrac=0.0 )
			except : lc_clean(locclu+'/lc/lc'+obsids[i]+'.fits', outfile=locclu+'/lc/gti'+obsids[i]+'.fits' , mean=10.0 , clip=3.0 , scale=1.2, minfrac=0.0 )


			cont = 'n'
			while cont == 'n' or cont == 'N' :
				print
				cont = raw_input ('Continue? (y/n) ')
				print
				print
				if cont =='y' or cont=='Y' : break
				else : cont = 'n'
				if os.path.exists(locclu+'/lc/gti'+obsids[i]+'.fits') : os.remove(locclu+'/lc/gti'+obsids[i]+'.fits')
				readline.add_history('lc_clean("'+locclu+'/lc/lc'+obsids[i]+'.fits" , outfile="'+locclu+'/lc/gti'+obsids[i]+'.fits"  , mean=None , clip=3.0 , scale=1.2, minfrac=0.0 )')
				input('New lc_clean ---> ')


		os.system('punlearn dmcopy')
		os.system('dmcopy "mfe_'+clu+'/repro3_evt'+obsids[i]+'.fits[ccd_id=0,1,2,3,5,6,7][@'+locclu+'/lc/gti'+obsids[i]+'.fits]" '+locclu+'/clean'+obsids[i]+'.fits')


	print
	print







