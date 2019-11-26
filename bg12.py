
# PURPOSE: create background data, from the blank sky data in the Chandra
# PURPOSE: calibration dataset (CALDB)

# v2: for imlr.py : allows for BG for any chip not just from 0:3,6,7
#                 : actually, just chips 0:3,5,6,7 (4, 8 and 9 are messy!)

# v5: * based on v2
#     * introduced my own blank-sky selection method, instead of acis_lookup_bkgrnd
#        --> period C evt data in VFAINT mode, are assigned blank sky data from
#            period D, since all period C blank skys are in FAINT mode.
#     * filter all blank sky data by status=0, except when using a VFAINT BG for a FAINT evt file.
#       In this case use filter status=00000000x00000000000000000000000 to keep BG events in blanksky file

# v6: * added the region reg/hibgcornerOB.reg, which is at the corner of chip #1, and which
#       has a higher (and variable?) BG rate, to be excluded from further analysis
#     * makes BG even if it exists
#     * ACIS-I and ACIS-S are called in sequence, as opposed to the mess in bg5
#     * no ACIS-I chips when analyzing ACIS-S observation and vice-versa
#     * choose ACIS-S or ACIS-I based on hardsoft2, and not the on-the-fly method used in bg5

# v7: * use asol3.py to find asol, flt1, pbk, msk, bpix or evt1 files

# v8: * minor change in screen output when BG file not found. (made after in v7 after it was run)
#     * more minor screen output changes in v8
#     * Main change: all CCD's are processed, for all OBSIDS: if ACIS-I data present in
#       ACIS-S pointing, it is processed.
#     * changing names of ACIS-I bg files to bg/bgOB_ccdir.fits, so that we can consider
#       the 4 ACIS-I chips, as one chip, chip # "i", keeping in minf that chip i, might not
#       always contain all 4 ACIS-I chips, sometimes, it will be a subset of them.
#     * no longer over-writes BG files
#     * make bg file of a given CCD if it has any counts, as opposed to at least 500 in prev version

# v9: * no more _nxt datasets
#     * the dates at the boundaries beteween periods were skipped because < was used instead of <=, 
#       so I changed the inequalities and the date of period B to include all dates.

# v10: * 2017-03-31: new version of asol: asol5_loc

# v11: * uses hardsoft3
#      * adds time column to bg/bgOB_ccdN.fits because it is required by
#        acis_process_events (which is a bug according to CIAO website)

# OUTPUT: bg/bgOB_ccdir.fits for ACIS-I ... or .... bg/bgOB_ccd[567]r.fits for ACIS-S

# INTERMEDIATE: bg/bgOB_ccdN.fits (soft link)
# INTERMEDIATE: bg/bgOBr_allstat.fits
# INTERMEDIATE: bg/bgOB_BADGAIN.fits

import commands
import os
import glob
import pdb

import asol6_loc
import hardsoft3

def bg(clu, ob) :

	print
	print 'Making BG files for OBSID '+ob+' ---------------------------'
	print

	locclu = 'mfe_' + clu

	acao = asol6_loc.asol(ob, 'asol')
	datamode = commands.getoutput('dmkeypar '+locclu+'/clean'+ob+'.fits datamode echo+')

	sd = commands.getoutput ( 'dmkeypar '+locclu+'/clean'+ob+'.fits date-obs echo+' )
	intdat = int(sd[0:4]+sd[5:7]+sd[8:10])
	if intdat < 19990916 : 
		bgd = 'period a' # period A : BG too uncertain because CTI rapidly increasing!
	elif 19990916 <= intdat and intdat <= 20000128 : 
		bgd = '1999-09-16' # period B
	elif 20000129 <= intdat and intdat <= 20001130 : 
		if   datamode == 'FAINT'  : 
			bgd = '2000-01-29' # period C
		elif datamode == 'VFAINT' : 
			bgd = '2000-12-01' # period C by date, but assign to period D because no VFAINT mode in period C blank-sky files
	elif 20001201 <= intdat and intdat <= 20050831 : 
		bgd = '2000-12-01' # period D
	elif 20050901 <= intdat : 
		bgd = '2005-09-01'  # period E

	hardob = hardsoft3.hardsoft(clu, [ob], 'h')
	softob = hardsoft3.hardsoft(clu, [ob], 's')
	if   len(hardob) == 1 : bgaim = 'i'
	elif len(softob) == 1 : bgaim = 's'
	else :
		print '------------ CANT DETERMINE AIM POINT -----------------'
		print
		bgaim = 'unknown'
	


	chipnum = []
	allchips='0123567'
	for c in allchips :
		if float(commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[ccd_id='+c+']" counts')) > 0 : chipnum.append(c)

	####### Make region for the high-BG corner
	if '1' in chipnum :
		os.system('punlearn dmcoords')
		os.system('dmcoords '+locclu+'/clean'+ob+'.fits option=chip chip_id=1 chipx=900 chipy=512 asol='+acao)
		hibgx = commands.getoutput('pget dmcoords x')
		hibgy = commands.getoutput('pget dmcoords y')
		roll = float( commands.getoutput('dmkeypar '+locclu+'/clean'+ob+'.fits roll_pnt echo+') )

		hibgfile=open(locclu+'/reg/hibgcorner'+ob+'.reg', 'w')
		hibgfile.write('rotbox('+hibgx+','+hibgy+',1050,300,'+str(360.-roll)+')\n')
		hibgfile.close()
	############################################3


	if   commands.getoutput ( 'dmkeypar '+locclu+'/clean'+ob+'.fits cti_corr echo+' ) == '1' : bgcti='_cti'
	elif commands.getoutput ( 'dmkeypar '+locclu+'/clean'+ob+'.fits cti_corr echo+' ) == '0' : bgcti=''

	ciaodir=commands.getoutput('echo $ASCDS_INSTALL')
	bgf = ciaodir+'/CALDB/data/chandra/acis/bkgrnd'
	if bgaim != 'unknown' and bgd != 'period a' :

		for c in chipnum :

			bgfound = glob.glob(bgf + '/acis'+c+bgaim+'D'+bgd+'bkgrnd'+bgcti+'N000*.fits')
			bgfound.sort()

			if   len(bgfound) == 0 : 
				print "---- no BG files found for cluster "+clu+", OBSID "+ob+", chip# "+c+" ------------------ "
				print 'The following was not found: '
				print bgf + '/acis'+c+bgaim+'D'+bgd+'bkgrnd'+bgcti+'N000*.fits'
				print
			elif len(bgfound) == 1 : 
				if os.path.exists(locclu+'/bg/bg'+ob+'_ccd'+c+'.fits') : os.remove(locclu+'/bg/bg'+ob+'_ccd'+c+'.fits')
				os.symlink(bgfound[0], locclu+'/bg/bg'+ob+'_ccd'+c+'.fits')
			elif len(bgfound) > 1 :
				print
				print ' -------------- MANY BG FILES FOUND FOR cluster '+clu+', OBSID '+ob+', CCD '+c
				print ' -------------- CHOOSING THE LATEST: '+bgfound[-1]
				print
				if os.path.exists(locclu+'/bg/bg'+ob+'_ccd'+c+'.fits') : os.remove(locclu+'/bg/bg'+ob+'_ccd'+c+'.fits')
				os.symlink(bgfound[-1], locclu+'/bg/bg'+ob+'_ccd'+c+'.fits')



	chipnums=[]
	for s in ['5', '6', '7'] :
		if s in chipnum : chipnums.append(s)

	chipnumi=[]
	for s in ['0', '1', '2', '3'] :
		if s in chipnum : chipnumi.append(s)

	if len(chipnumi) > 0  : 
		makebgi(ob, locclu, datamode, bgd, acao) 
	if len(chipnums) > 0  : 
		makebgs(ob, locclu, chipnums, datamode, bgd, acao) 

	print
	print
















def makebgi(ob, locclu, datamode, bgd, acao) :

	tobemerged = glob.glob(locclu+'/bg/bg'+ob+'_ccd[0123].fits')

	if len(tobemerged) > 0 :

		# BG files to be merged, corresponding to ACIS 0,1,2,3 blank sky files
		mrgstr=''
		for i in range(len(tobemerged)) :
			mrgstr = mrgstr + ',' + tobemerged[i]
		mrgstr=mrgstr[1:]

		# Merge the BG files into 1
		os.system('punlearn dmmerge')
		os.system('dmmerge "'+mrgstr+'" '+locclu+'/bg/bg'+ob+'_ccdi.fits mode=h')
		expos = commands.getoutput('dmkeypar '+tobemerged[0]+' exposure echo+')
		os.system('dmhedit '+locclu+'/bg/bg'+ob+'_ccdi.fits filelist= op=add key=EXPOSURE value='+expos+' unit=s')

		bggain  = commands.getoutput('dmkeypar '+locclu+'/bg/bg'+ob+'_ccdi.fits gainfile echo+')
		evtgain = commands.getoutput('dmkeypar '+locclu+'/clean'+ob+'.fits gainfile echo+')

		# Match the evt and (merged) bg gain files if they're different
		if bggain != evtgain :

			print
			print '---- evt and BG for OBSID '+ob+' have different gain files.. reprocessing BG file:'
			print

			#os.rename(locclu+'/bg/bg'+ob+'_ccdi.fits', locclu+'/bg/bg'+ob+'_ccdi_BADGAIN.fits')
			os.system('punlearn dmtcalc')
			os.system('dmtcalc '+locclu+'/bg/bg'+ob+'_ccdi.fits '+locclu+'/bg/bg'+ob+'_ccdi_BADGAIN.fits expression="time=#nan;expno=0"')
			os.remove(locclu+'/bg/bg'+ob+'_ccdi.fits')

			os.system('punlearn acis_process_events')
			os.system('pset acis_process_events infile='+locclu+'/bg/bg'+ob+'_ccdi_BADGAIN.fits')
			os.system('pset acis_process_events outfile='+locclu+'/bg/bg'+ob+'_ccdi.fits')
			os.system('pset acis_process_events acaofffile=NONE stop="none" doevtgrade=no')
			os.system('pset acis_process_events apply_cti=yes apply_tgain=no calculate_pi=yes')
			os.system('pset acis_process_events gainfile="$CALDB/data/chandra/acis/det_gain/'+evtgain+'"')
			os.system('pset acis_process_events eventdef="{s:ccd_id,s:node_id,i:expno,s:chip,s:tdet,f:det,f:sky,s:phas,l:pha,l:pha_ro,f:energy,l:pi,s:fltgrade,s:grade,x:status}"')
			os.system('acis_process_events mode=h')
			os.remove(locclu+'/bg/bg'+ob+'_ccdi_BADGAIN.fits')

		# reproject events of bg to match evt file
		print
		print "Reprojecting bg to match evt file ...."
		os.system('punlearn reproject_events ')
		os.system('pset reproject_events infile='+locclu+'/bg/bg'+ob+'_ccdi.fits')
		os.system('pset reproject_events mode=h')
		os.system('pset reproject_events outfile='+locclu+'/bg/bg'+ob+'_ccdir_allstat.fits')
		os.system('pset reproject_events aspect='+acao)
		os.system('pset reproject_events match='+locclu+'/clean'+ob+'.fits')
		os.system('pset reproject_events random=0')
		os.system('reproject_events')

		#os.remove(locclu+'/bg/bg'+ob+'_ccdi.fits')

		os.system('punlearn dmcopy')
		if datamode == 'FAINT' and ( bgd == '2000-12-01' or bgd == '2005-09-01' ) : 
			os.system  ('dmcopy "'+locclu+'/bg/bg'+ob+'_ccdir_allstat.fits[status=00000000x00000000000000000000000]" '+locclu+'/bg/bg'+ob+'_ccdir.fits')
		else : os.system('dmcopy "'+locclu+'/bg/bg'+ob+'_ccdir_allstat.fits[status=0]" '+locclu+'/bg/bg'+ob+'_ccdir.fits')

		#os.remove(locclu+'/bg/bg'+ob+'_ccdir_allstat.fits')



















def makebgs(ob, locclu, chipnum, datamode, bgd, acao) :

	# This paragraph added for the case where the BG file is not found
	# in the _nxt banksky dataset
	absentchip=[]
	for c in chipnum :
		if not os.path.exists(locclu+'/bg/bg'+ob+'_ccd'+c+'.fits') : absentchip.append(c)
	for c in absentchip : chipnum.remove(c)

	# Match the evt and (merged) bg gain files if they're different
	evtgain = commands.getoutput('dmkeypar '+locclu+'/clean'+ob+'.fits gainfile echo+')

	for i in range(len(chipnum)) :

		bggain  = commands.getoutput('dmkeypar '+locclu+'/bg/bg'+ob+'_ccd'+chipnum[i]+'.fits gainfile echo+')

		if bggain != evtgain :		

			print '#### evt and BG for OBSID '+ob+' have different gain files.. reprocessing BG file:'
			print

			#os.rename(locclu+'/bg/bg'+ob+'_ccd'+chipnum[i]+'.fits', locclu+'/bg/bg'+ob+'_ccd'+chipnum[i]+'_BADGAIN.fits')
			os.system('punlearn dmtcalc')
			os.system('dmtcalc '+locclu+'/bg/bg'+ob+'_ccd'+chipnum[i]+'.fits '+locclu+'/bg/bg'+ob+'_ccd'+chipnum[i]+'_BADGAIN.fits expression="time=#nan;expno=0"')
			os.remove(locclu+'/bg/bg'+ob+'_ccd'+chipnum[i]+'.fits')

			os.system('punlearn acis_process_events')
			os.system('pset acis_process_events infile='+locclu+'/bg/bg'+ob+'_ccd'+chipnum[i]+'_BADGAIN.fits')
			os.system('pset acis_process_events outfile='+locclu+'/bg/bg'+ob+'_ccd'+chipnum[i]+'.fits')
			os.system('pset acis_process_events acaofffile=NONE stop="none" doevtgrade=no')
			os.system('pset acis_process_events apply_cti=yes apply_tgain=no calculate_pi=yes')
			os.system('pset acis_process_events gainfile="$CALDB/data/chandra/acis/det_gain/'+evtgain+'"')
			os.system('pset acis_process_events eventdef="{s:ccd_id,s:node_id,i:expno,s:chip,s:tdet,f:det,f:sky,s:phas,l:pha,l:pha_ro,f:energy,l:pi,s:fltgrade,s:grade,x:status}"')
			os.system('acis_process_events mode=h')
			os.remove(locclu+'/bg/bg'+ob+'_ccd'+chipnum[i]+'_BADGAIN.fits')



	# reproject events of bg to match evt file
	print
	print "Reprojecting bg to match evt file ...."
	os.system('punlearn reproject_events ')
	os.system('pset reproject_events mode=h')
	os.system('pset reproject_events aspect='+acao)
	os.system('pset reproject_events random=0')
	os.system('pset reproject_events match='+locclu+'/clean'+ob+'.fits')


	for i in range(len(chipnum)) :

		os.system('pset reproject_events infile='+locclu+'/bg/bg'+ob+'_ccd'+chipnum[i]+'.fits')
		os.system('pset reproject_events outfile='+locclu+'/bg/bg'+ob+'_ccd'+chipnum[i]+'r_allstat.fits')
		os.system('reproject_events')

		os.system('punlearn dmcopy')
		if datamode == 'FAINT' and ( bgd == '2000-12-01' or bgd == '2005-09-01' ) : 
			os.system  ('dmcopy "'+locclu+'/bg/bg'+ob+'_ccd'+chipnum[i]+'r_allstat.fits[status=00000000x00000000000000000000000]" '+locclu+'/bg/bg'+ob+'_ccd'+chipnum[i]+'r.fits')
		else : os.system('dmcopy "'+locclu+'/bg/bg'+ob+'_ccd'+chipnum[i]+'r_allstat.fits[status=0]" '+locclu+'/bg/bg'+ob+'_ccd'+chipnum[i]+'r.fits')

		#os.remove(locclu+'/bg/bg'+ob+'_ccd'+chipnum[i]+'r_allstat.fits')
















