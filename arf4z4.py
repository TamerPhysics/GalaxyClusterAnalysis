
# PURPOSE: create ARF files for the purpose of fitting Z(r) spectra

# v2: * use spN_4.fits instead of _3.fits
#     * use asol3.py to find asol, flt1, pbk, msk, bpix or evt1 files

# arf4kt9 based on arf4radp2.py

# arf4z2 based on arf4radp2.py

# arf4z4 * based on arf4z3loc
#        * no pbkfile parameter

import glob
import os

import asol6_loc

def arf (clu, ob, isp, ccd) :

	locclu = 'mfe_'+clu


	# MKWARF --------------

	pbk = asol6_loc.asol(ob,'pbk')
	acao = asol6_loc.asol(ob, 'asol')

	os.system('punlearn mkwarf')
	os.system('pset mkwarf infile="'+locclu+'/zofr/sp10simple'+ob+'_'+str(isp)+'_ccd'+ccd+'.fits[WMAP]"')
	os.system('pset mkwarf outfile='+locclu+'/zofr/warf10_'+ob+'_'+str(isp)+'_ccd'+ccd+'.fits')
	os.system('pset mkwarf weightfile='+locclu+'/zofr/wfef10_'+ob+'_'+str(isp)+'_ccd'+ccd+'.fits')
#	os.system('pset mkwarf pbkfile='+pbk)
	os.system('pset mkwarf egridspec="0.3:11.0:0.01"')
	os.system('pset mkwarf asolfile='+acao)
	os.system('mkwarf verbose=1 mode=h clobber=yes')

	os.system('rm '+locclu+'/zofr/wfef10_'+ob+'_'+str(isp)+'_ccd'+ccd+'.fits')
