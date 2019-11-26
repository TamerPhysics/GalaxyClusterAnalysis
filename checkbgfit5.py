
# PURPOSE: Create background spectra and fit them

# v2: * if ARF not present, get one, and set it equal to zero
# v3: * correcting ikt vs ktish0 mistake in above version...
#     * more complex fitting for case of multiple kT-components in one Z-shell
# v4: * write reduced statistic to file
#     * reduce the kT components (for a given Z-shell) to those whose uncertainty ranges 
#       do not overlap
# v5: * perform the last fit, with the normalizations of the different kT components unlinked.
#     * change None to NaN in output files

# ktofr9: * based on zofr5.py
#         * normalizations are only linked for spectra that cover 95-105% of the annulus

# ktofr9_norm: * based on ktofr9, 
#              * only fit norm's while getting kT results from ktofr9.py results
#              * OUTPUT IS NOT WELL WRITTEN, ISSUES FIXED IN ktofr9_norm2

# ktofr9_norm2: * writes to ktofr9_norm2.txt the normalization of each CCD in each OBSID
#               * also writes the covered area a fraction of the complete shell area

# ktofr10: * #++# stands for changes in this version, from ktofr9_norm2
#          * nnec changed, after including the simulations with BGfrac=0.8
#          * exclude CCD's which don't have BG files, while making evtrad files.
#          * ccdname was = 'i576' instead of 'i567', this was fixed
#          * Systematic errors of the hien run not used anymore
#          * THIS CHANGE WAS CANCELLED: delete_window replaced by clear as it works when no windows are open
#          * define bgxs to contain BG scaling factor, relative to blank-sky BG count rate
#          * when linking norms with area almost equal to area of shell, use relative backscal, instead 
#            of equating them

# ktofr11: * find which CCD's are in an OBSID in the very begining, to set presccdglob early on and avoid
#            including data in the analysis which will not be used. For example evtrad file should only
#            contain events from OBSIDS which will actually be analyzed. In previous versions, this caused
#            bgfrac to be wrongly calculated. This is done by introducing whichccd2, and moving some of 
#            of region output files from radreg3 to whichccd2, while making radreg4 instead of 3.
#          * write nH
#          * restrict calculating sclbgct to the CCD's in presccdglob, instead of adding BG
#            counts from all CCD's
#          * freeze nH! (as done in Vikhlinin)
#          * kT for nnec comes from previous shell, unless it's shell 0, then use values in cluname_TZ.txt
#          * initial kT and Z for any shell (except ishell=0) comes from previous shell

# ktofr11_fix: * refit the last n bins in the spectra, where n is taken from fix_ktofr11_list.txt

# ktofr12loc: * in the hien scaling, only instrumental BG components are scaled by hienbgscl.
#               sky BG components are scaled by bgscl
#             * use simple box regions to extract spectra
#             * min of .norm param is 0, not emnorm0/1e3

# checkbgfit : based on the first part of ktofr12loc

# checkbgfit5: remove reference to acchif folder

# OUTPUT: ktofrsp/spfix*
# OUTPUT: ktofrsp/warf12_*
# OUTPUT: ktofrsp/kt12an*
# OUTPUT: ktofr12.txt
# OUTPUT: ktofr12_hien.txt
# OUTPUT: ktofr12_norm3.txt

from sherpa.astro.ui import *
from pychips.all import *

import os
import commands
import pdb
import gc

import numpy as n

import centradec
import loadbg4z8simple2_write
import whichccd2z_loc
import radreg4simple2
import asol6_loc

def ktofr(clu, obsids, zz, nnhhlab, checkfit) :

	print ' -----------> checkbgfit4 for '+ clu

	clean()
	gc.collect()
	set_stat('cstat')

	locclu = 'mfe_'+clu

	(rac, decc) = centradec.getrd(locclu)

	if not os.path.exists(locclu+'/bginfield/') : os.mkdir(locclu+'/bginfield/')
	if not os.path.exists(locclu+'/reg/') : os.mkdir(locclu+'/reg/')
	if not os.path.exists(locclu+'/spec/') : os.mkdir(locclu+'/spec/')
	if not os.path.exists(locclu+'/ktofrsp/') : os.mkdir(locclu+'/ktofrsp/')

	ccdname='i567'
	ccdlist=['0,1,2,3','5','6','7']
	ccdindex = {'i':0, '5':1, '6':2, '7':3}

#	dataplot  = ChipsCurve()
#	modelplot = ChipsHistogram()

#	dataplot.symbol.style='plus'
#	dataplot.symbol.size=3
#	dataplot.line.style='none'
#	modelplot.line.color = "red"
#	modelplot.line.thickness = 3

	# create simplechipsreg regions.
	for ob in obsids :
		os.system('punlearn skyfov')
		os.system('skyfov '+locclu+'/clean'+ob+'.fits '+locclu+'/reg/simplechipsreg'+ob+'.fits clobber=yes')

	# which CCD's are present in each OBSID
	presccdglob=[]
	for iob in range(len(obsids)) : presccdglob.append( whichccd2z_loc.find(clu, obsids[iob]) )

	# Remove OBSID's which do not have any good CCD's:
	emptyobs=[]
	for iob in range(len(obsids)) :
		ob = obsids[iob]
		for ccd in str(presccdglob[iob]) : 
			if not os.path.exists(locclu+'/spec/sp'+ob+'_ccd'+ccd+'_center.wrmf') \
                or not os.path.exists(locclu+'/spec/bgarf'+ob+'_ccd'+ccd+'.fits') \
                or not os.path.exists(locclu+'/spec/bgspsimple'+ob+'_ccd'+ccd+'.fits') \
                or not os.path.exists(locclu+'/spec/bgrmf'+ob+'_ccd'+ccd+'.fits') :
				presccdglob[iob] = presccdglob[iob].replace(ccd, '')
		if presccdglob[iob]=='' : emptyobs.append(iob)
	emptyobs.reverse()
	for iob in emptyobs : 
		del obsids[iob]
		del presccdglob[iob]

	anfile = open(locclu+'/bginfield_all.reg', 'r')
	anline=anfile.readline()
	anline=anfile.readline()
	pbg1=float(anline.split(',')[2][0:-1])*60/0.492
	pbg2=float(anline.split(',')[3][0:-3])*60/0.492
	anfile.close()

	# Create BG annulus for each OBSIDS, which excludes the pt srcs
	# Get presccd
	presccd=[]
	for iob in range(len(obsids)) :
		ob = obsids[iob]
		intermed1 = radreg4simple2.radreg(clu, ob, pbg1, pbg2, locclu+'/bginfield/outer'+ob)
		presccd.append( intermed1 )

	# Make sure that some needed files were created. If not, remove from prescccd
	# any CCD wich is missing PI, ARF or RMF file
	for iob in range(len(obsids)) :
		ob = obsids[iob]
		for ccd in str(presccd[iob]) :
			if not os.path.exists(locclu+'/spec/sp'+ob+'_ccd'+ccd+'_center.wrmf') or ccd not in presccdglob[iob] : presccd[iob] = presccd[iob].replace(ccd, '')
			#if not os.path.exists(locclu+'/spec/sp'+ob+'_ccd'+ccd+'_center.wrmf') : presccd[iob] = presccd[iob].replace(ccd, '')

	srcstr=''
	for iob in range(len(obsids)) :
		ob = obsids[iob]
		for iccd in range(len(presccd[iob])) :
			ccd=presccd[iob][iccd]

			os.system('punlearn dmextract')
			os.system('dmextract "'+locclu+'/clean'+ob+'.fits[ccd_id='+ccdlist[ccdindex[ccd]]+'][sky=region('+locclu+'/bginfield/outer'+ob+'_ccd'+ccd+'_xfov_pt_simple.reg)][bin PI]" '+locclu+'/bginfield/sp'+ob+'_ccd'+ccd+'.fits opt=pha1 wmap="[energy=300:2000][bin det=8]" clobber=yes' )

			# Make ARF given the WMAP from the above spectrum file
			#pbk = asol6_loc.asol(ob,'pbk')
			acao = asol6_loc.asol(ob, 'asol')
			os.system('punlearn mkwarf')
			os.system('pset mkwarf infile="'+locclu+'/bginfield/sp'+ob+'_ccd'+ccd+'.fits[WMAP]"')
			os.system('pset mkwarf outfile='+locclu+'/bginfield/warf'+ob+'_ccd'+ccd+'.fits')
			os.system('pset mkwarf weightfile='+locclu+'/bginfield/wfef'+ob+'_ccd'+ccd+'.fits')
			#os.system('pset mkwarf pbkfile='+pbk)
			os.system('pset mkwarf egridspec="0.3:11.0:0.01"')
			os.system('pset mkwarf asolfile='+acao)
			os.system('mkwarf verbose=1 mode=h clobber=yes')
			os.system('rm '+locclu+'/bginfield/wfef'+ob+'_ccd'+ccd+'.fits')


	# Make sure (AGAIN) that all the needed files were created. If not, remove from prescccd
	# any CCD wich is missing PI, ARF or RMF file
	for iob in range(len(obsids)) :
		ob = obsids[iob]
		for ccd in str(presccd[iob]) :
			if not os.path.exists(locclu+'/bginfield/warf'+ob+'_ccd'+ccd+'.fits') or not os.path.exists(locclu+'/bginfield/sp'+ob+'_ccd'+ccd+'.fits') or ccd not in presccdglob[iob] : 
				presccd[iob]     =     presccd[iob].replace(ccd, '')
			if not os.path.exists(locclu+'/spec/sp'+ob+'_ccd'+ccd+'_center.wrmf') :
				presccdglob[iob] = presccdglob[iob].replace(ccd, '')
				presccd[iob]     =     presccd[iob].replace(ccd, '')

	# load infield BG
	for iob in range(len(obsids)) :
		ob = obsids[iob]
		for iccd in range(len(presccd[iob])) :
			ccd=presccd[iob][iccd]

			load_pha( ob+'.'+ccd , locclu+'/bginfield/sp'+ob+'_ccd'+ccd+'.fits' )
			load_rmf( ob+'.'+ccd , locclu+'/spec/sp'+ob+'_ccd'+ccd+'_center.wrmf'  )
			load_arf( ob+'.'+ccd , locclu+'/bginfield/warf'+ob+'_ccd'+ccd+'.fits' )

			srcstr = srcstr + '"'+ob+'.'+ccd + '",'





	# Fit the BG using blank sky AND infield BG data
	loadbg4z8simple2_write.bg(clu, obsids, presccdglob, presccd, srcstr, zz, nnhhlab, checkfit)




	# Make and save plots
	for iob in range(len(obsids)) :
		ob = obsids[iob]
		for iccd in range(len(presccdglob[iob])) :
			ccd=presccdglob[iob][iccd]

			if ccd in presccd[iob]:
				plot_fit_resid(ob+'.'+ccd)
				if os.path.exists('fitplots/checkbgfit4infield_'+clu+'_'+ob+'_ccd'+ccd+'.eps') : os.remove('fitplots/checkbgfit4infield_'+clu+'_'+ob+'_ccd'+ccd+'.eps')
				print_window( 'fitplots/checkbgfit4infield_'+clu+'_'+ob+'_ccd'+ccd, ['format', 'eps', 'orientation', 'landscape'])
				delete_window('all')

			plot_fit_resid('bg'+ob+'.'+ccd)
			if os.path.exists('fitplots/checkbgfit4blanksky_'+clu+'_'+ob+'_ccd'+ccd+'.eps') : os.remove('fitplots/checkbgfit4blanksky_'+clu+'_'+ob+'_ccd'+ccd+'.eps')
			print_window( 'fitplots/checkbgfit4blanksky_'+clu+'_'+ob+'_ccd'+ccd, ['format', 'eps', 'orientation', 'landscape'])
			delete_window('all')





