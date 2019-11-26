
# PURPOSE: create exposure maps for all observations of a cluster

# v2: for mfe_*/
#  no soft ccds (ACIS 5, 6, 7)

# v3: for entropy/
#  makes either exposure map for inner or outer region, dependin on the input parameter inorout
#  no more spectrum version

# v4: * use asol3.py to find asol, flt1, pbk, msk, bpix or evt1 files
#     * make sure you dont call ccds[0] if ccds is empty.

# v6: * 2017-03-31: new version of asol: asol5_loc

# v7: *uses hardsoft3
#     * no pbk keywords

# INPUT: spec/srcin.txt or spec/srcout.txt (input to binsrc2.pro)
# INPUT: rac.txt and decc.txt
# INPUT: cleanOB.fits
# INPUT: asol, pbk and mask files (from chandra dataset)

# OUTPUT: evtOB_c.fits (energy=[700:2000], with all chips)
# OUTPUT: instmap/asphistOB_ccdCHIP.fits
# OUTPUT: instmap/expmapOB_{in/out}.fits ---> for ACIS-I
# OUTPUT: instmap/expmapOB_ccdCHIP_{in/out}.fits ---> for ACIS-S
# OUTPUT: mfe_CLU/spec/binsrcin.txt OR mfe_CLU/spec/binsrcout.txt (through binsrc2.pro)

import os
import glob
import commands

import asol6_loc
import centradec
import hardsoft3
import binsrc2

def exp(clu, obsids, inorout) :

	locclu = 'mfe_'+clu
	origclu = 'orig_' + clu

	#os.system( 'echo "binsrc2, \''+clu+'\' , \'src'+inorout+'.txt\' " | idl' )

	binsrc2.bin(clu, 'src'+inorout+'.txt')

	(rac, decc) = centradec.getrd(locclu)

	if not os.path.exists(locclu+'/instmap') : os.mkdir(locclu+'/instmap')

	hardobs=hardsoft3.hardsoft(clu, obsids, 'h')
	softobs=hardsoft3.hardsoft(clu, obsids, 's')

	for i in range(len(obsids)) : 
		ob = obsids[i]

		acao = asol6_loc.asol(ob, 'asol')

		os.system('punlearn ardlib')
		os.system('acis_set_ardlib ' + asol6_loc.asol(ob,'bpix') )

		ccds=''

		# if this is the in region, then check if there are counts in in region
		# but if it's the out region, only check if there are counts in the ccd of interest:
		# the reason is we want to include data that is beyond the out region, which stops
		# at 1Mpc
		if inorout == 'in' :
			if int(commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/'+inorout+'.reg)][ccd_id=0]" counts')) > 0 : ccds = ccds + '0'
			if int(commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/'+inorout+'.reg)][ccd_id=1]" counts')) > 0 : ccds = ccds + '1'
			if int(commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/'+inorout+'.reg)][ccd_id=2]" counts')) > 0 : ccds = ccds + '2'
			if int(commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/'+inorout+'.reg)][ccd_id=3]" counts')) > 0 : ccds = ccds + '3'
			if int(commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/'+inorout+'.reg)][ccd_id=5]" counts')) > 0 : ccds = ccds + '5'
			if int(commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/'+inorout+'.reg)][ccd_id=6]" counts')) > 0 : ccds = ccds + '6'
			if int(commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/'+inorout+'.reg)][ccd_id=7]" counts')) > 0 : ccds = ccds + '7'
		else :
			if int(commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[ccd_id=0]" counts')) > 0 : ccds = ccds + '0'
			if int(commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[ccd_id=1]" counts')) > 0 : ccds = ccds + '1'
			if int(commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[ccd_id=2]" counts')) > 0 : ccds = ccds + '2'
			if int(commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[ccd_id=3]" counts')) > 0 : ccds = ccds + '3'
			if int(commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[ccd_id=5]" counts')) > 0 : ccds = ccds + '5'
			if int(commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[ccd_id=6]" counts')) > 0 : ccds = ccds + '6'
			if int(commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[ccd_id=7]" counts')) > 0 : ccds = ccds + '7'

#		ccdscoma=''
#		for j in range(len(ccds)-1) : ccdscoma=ccdscoma+ccds[j]+','
#		ccdscoma=ccdscoma+ccds[-1]

		if not os.path.exists(locclu+'/evt'+ob+'_c.fits') :
			os.system('punlearn dmcopy')
			os.system('dmcopy "'+locclu+'/clean'+ob+'.fits[energy=700:2000]" '+locclu+'/evt'+ob+'_c.fits')
			#[ccd_id='+ccdscoma+']

		for j in range(len(ccds)) :
			if not os.path.exists( locclu+'/instmap/asphist'+ob+'_ccd'+ccds[j]+'.fits' ) :
				os.system('punlearn asphist')
				os.system('pset asphist infile='+acao)
				os.system('pset asphist outfile='+locclu+'/instmap/asphist'+ob+'_ccd'+ccds[j]+'.fits')
				os.system('pset asphist evtfile="'+locclu+'/clean'+ob+'.fits[ccd_id='+ccds[j]+']"')
				os.system('asphist verbose=1 mode=h')


		# Make instrument map
		for j in range(len(ccds)) :
			os.system('punlearn mkinstmap')
			os.system('pset mkinstmap pixelgrid="1:1024:#1024,1:1024:#1024"')
			os.system('pset mkinstmap maskfile='+asol6_loc.asol(ob,'msk1'))
			#os.system('pset mkinstmap pbkfile='+asol6_loc.asol(ob,'pbk'))
			os.system('pset mkinstmap dafile=CALDB')
			os.system('pset mkinstmap spectrumfile='+locclu+'/spec/binsrc'+inorout+'.txt')
			os.system('pset mkinstmap obsfile='+locclu+'/clean'+ob+'.fits')
			os.system('pset mkinstmap detsubsys=ACIS-'+ccds[j])
			os.system('pset mkinstmap outfile='+locclu+'/instmap/winstmap'+ob+'_ccd'+ccds[j]+'_'+inorout+'.fits')
			os.system('mkinstmap mode=h verbose=1')


		os.system('punlearn dmcoords')
		os.system('dmcoords '+locclu+'/clean'+ob+'.fits asolfile='+acao+' option=cel ra='+rac+' dec='+decc+' celfmt=hms')
		xc = commands.getoutput('pget dmcoords x')
		yc = commands.getoutput('pget dmcoords y')	
		xcf = float(xc)
		ycf = float(yc)


		xs=commands.getoutput('dmlist "'+locclu+'/evt'+ob+'_c.fits[cols x]" data,clean')
		xsv=xs.split('\n')
		xsv.remove(xsv[0])
		for ii in range(len(xsv)) : 
			xsv[ii]=float(xsv[ii])

		ys=commands.getoutput('dmlist "'+locclu+'/evt'+ob+'_c.fits[cols y]" data,clean')
		ysv=ys.split('\n')
		ysv.remove(ysv[0])
		for ii in range(len(ysv)) : 
			ysv[ii]=float(ysv[ii])

		xmin = min(xsv)
		xmax = max(xsv)
		ymin = min(ysv)
		ymax = max(ysv)


		# Make exposure map
		os.system('punlearn mkexpmap')
		os.system('pset mkexpmap normalize=yes')
		os.system('pset mkexpmap xygrid='+str(xmin)+':'+str(xmax)+':4,'+str(ymin)+':'+str(ymax)+':4')
		os.system('pset mkexpmap useavgaspect=no')
		os.system('pset mkexpmap mode=h')

		for j in range(len(ccds)) :
			os.system('mkexpmap instmapfile='+locclu+'/instmap/winstmap'+ob+'_ccd'+ccds[j]+'_'+inorout+'.fits outfile='+locclu+'/instmap/expmap'+ob+'_ccd'+ccds[j]+'_'+inorout+'.fits asphistfile='+locclu+'/instmap/asphist'+ob+'_ccd'+ccds[j]+'.fits verbose=1')


		# if ACIS-I ccd's are preset, merge them.
		if '0' in ccds or '1' in ccds or '2' in ccds or '3' in ccds :

			# Reproject exposure map 
			hmrgstr=''
			for ccd in '0123' : 
				if ccd in ccds : hmrgstr = hmrgstr + ',' + locclu+'/instmap/expmap'+ob+'_ccd'+ccd+'_'+inorout+'.fits'
			hmrgstr = hmrgstr[1:]

			os.system('punlearn reproject_image')
			os.system('pset reproject_image infile="'+hmrgstr+'"')
			os.system('pset reproject_image outfile='+locclu+'/instmap/expmap'+ob+'_ccdi_'+inorout+'.fits')
			os.system('pset reproject_image matchfile='+locclu+'/instmap/expmap'+ob+'_ccd'+ccds[0]+'_'+inorout+'.fits')
			os.system('pset reproject_image method=average')
			os.system('reproject_image mode=h')

			os.system('dmhedit '+locclu+'/instmap/expmap'+ob+'_ccdi_'+inorout+'.fits filelist="" op=add key=BUNIT value="cm**2" mode=h')

		# clean up files not needed
		for j in range(len(ccds)) :
			os.remove(locclu+'/instmap/winstmap'+ob+'_ccd'+ccds[j]+'_'+inorout+'.fits')

	# more cleanup
	os.system('rm '+locclu+'/instmap/expmap*_ccd[0123]_*.fits')





