
# PURPOSE: create ARF file for the purpose of calculating surface brighness
# PURPOSE: radial profile

# v2: * use spN_4.fits instead of _3.fits
#     * use asol3.py to find asol, flt1, pbk, msk, bpix or evt1 files

# v3: * added ccd # to name

# v4: asol6_loc

import glob
import os

import asol6_loc

def arf (clu, ob, i, ccdname) :

	locclu = 'mfe_'+clu
	origclu= 'orig_'+clu


	# MKWARF --------------

	pbk = asol6_loc.asol(ob,'pbk')
	acao = asol6_loc.asol(ob, 'asol')

	os.system('punlearn mkwarf')
	os.system('pset mkwarf infile="'+locclu+'/radp'+ob+'/sp'+str(i)+'_ccd'+ccdname+'.fits[WMAP]"')
	os.system('pset mkwarf outfile='+locclu+'/radp'+ob+'/warf'+str(i)+'_ccd'+ccdname+'.fits')
	os.system('pset mkwarf weightfile='+locclu+'/radp'+ob+'/wfef'+str(i)+'_ccd'+ccdname+'.fits')
	#os.system('pset mkwarf pbkfile='+pbk)
	os.system('pset mkwarf egridspec="0.3:11.0:0.01"')
	os.system('pset mkwarf asolfile='+acao)
	os.system('mkwarf verbose=1 mode=h')
