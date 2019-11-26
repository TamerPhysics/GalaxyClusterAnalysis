
# PURPOSE: load BG data for them to be fit

# v3: * based on loadbgkt9
#     * no need to create BG files, they were created when loadbgkt9
#       was called from ktofr9.. just load them.

# v3simple: * use simple region reg/obsimple*.reg

# v6simple: * add analysis for ccd's that don't have infield BG

# v8simple2: remove ref to acchif folder

# INPUT: * bg/bg'+ob+'_ccd'+ccd+'r.fits
# INPUT: * reg/obOBccdN.reg : region of CCD N in OBSIS OB

# OUTPUT: * no new files made

import commands
import os
import pdb

from sherpa.astro.ui import *
from pychips.all import *

import numpy as n

import rmfarf4bg
import bgmodel_flex11

def bg(clu, obsids, presccdglob) : #, srcstr, firstiob, zz, nnhhlab) :

	locclu = 'mfe_'+clu

	ccdname='i567'
	ccdindex = {'i':0, '5':1, '6':2, '7':3}

	# Read outer temperature from cluname_TZ.txt
	tzfile = open('cluname_TZ.txt','r')
	intermed='init'
	while intermed != '' :
		intermed = tzfile.readline()
		if intermed.split()[0] == clu :
			ktout = float( intermed.split()[3] )
			about = float( intermed.split()[4] )
			break
	tzfile.close()

	# Load BG, set BG model, set source model --------------------------------------
	bgstr=''
	for iob in range(len(obsids)) : 
		ob = obsids[iob]
		for iccd in range(len(presccdglob[iob])) :
			ccd=presccdglob[iob][iccd]

			if os.path.exists(locclu+'/bg/bg'+ob+'_ccd'+ccd+'r.fits') :

				# Make BG spectrum from Blank-sky data using entire CCD (except for ACIS-I where we exclude
				# the region of high-background on chip #1)
				os.system('punlearn dmextract')
				os.system('dmextract "'+locclu+'/bg/bg'+ob+'_ccd'+ccd+'r.fits[sky=region('+locclu+'/reg/obsimple'+ob+'ccd'+ccd+'.reg)][bin PI]" '+locclu+'/spec/bgspsimple'+ob+'_ccd'+ccd+'.fits opt=pha1 wmap="[energy=300:2000][bin det=8]" mode=h')

				rmfarf4bg.rmf(locclu+'/bg/bg'+ob+'_ccd'+ccd+'r.fits', \
                                  locclu+'/spec/bgspsimple'+ob+'_ccd'+ccd+'.fits', \
                                  locclu+'/spec/bgarf'+ob+'_ccd'+ccd+'.fits', \
                                  locclu+'/spec/bgrmf'+ob+'_ccd'+ccd+'.fits', \
                                  ccd)

			if os.path.exists(locclu+'/spec/bgarf'+ob+'_ccd'+ccd+'.fits') and \
                  os.path.exists(locclu+'/spec/bgspsimple'+ob+'_ccd'+ccd+'.fits')  and \
                  os.path.exists(locclu+'/spec/bgrmf'+ob+'_ccd'+ccd+'.fits') :

				# Load BG file we just made, with an RMF, no ARF!
				load_pha( 'bg'+ob+'.'+ccd , locclu+'/spec/bgspsimple'+ob+'_ccd'+ccd+'.fits')
				load_rmf( 'bg'+ob+'.'+ccd , locclu+'/spec/bgrmf'+ob+'_ccd'+ccd+'.fits')

				bgstr = bgstr + '"bg'+ob+'.'+ccd + '",'

				# the dummy dataset only made to extract its response, and use it in the BG full model
				# with the Xray BG components (the non-particle component)
				load_pha( 'dummy'+ob+'.'+ccd , locclu+'/radp'+ob+'/sp0_ccd'+ccd+'.fits' )
				load_rmf( 'dummy'+ob+'.'+ccd , locclu+'/spec/bgrmf'+ob+'_ccd'+ccd+'.fits')
				load_arf( 'dummy'+ob+'.'+ccd , locclu+'/spec/bgarf'+ob+'_ccd'+ccd+'.fits')
				get_data('dummy'+ob+'.'+ccd).exposure = float(get_data('bg'+ob+'.'+ccd).exposure)

				# ========= BG model ===============
				exec 'delete_model("bg'+ob+'.'+ccd+'")'
				rsp0  = get_response( 'dummy'+ob+'.'+ccd )
				bgrsp = get_response( 'bg'+ob+'.'+ccd )
				notice_id('bg'+ob+'.'+ccd)
				#group_bins('bg'+ob+'.'+ccd, 120)
				ignore_id('bg'+ob+'.'+ccd, ':.3,7:')
				exec 'set_full_model( "bg'+ob+'.'+ccd+'", bgrsp( powlaw1d.p1_'+ob+'_'+ccd+' + powlaw1d.p2_'+ob+'_'+ccd+' + exp.e1_'+ob+'_'+ccd+' + gauss1d.g1_'+ob+'_'+ccd+' +  gauss1d.g2_'+ob+'_'+ccd+' + gauss1d.g3_'+ob+'_'+ccd+' + gauss1d.g5_'+ob+'_'+ccd+' ) + rsp0( xsphabs.bgab * (powlaw1d.cxbbg_'+ob+'_'+ccd+') + xsapec.gal1bg_'+ob+'_'+ccd+' )  )' 



				# Set for each OBSID, the same for all shells
				bgmodel_flex11.setbgparams(ob, ccd, False)



	infieldccd=[]
	infieldgalperbs=[]
	infieldcxbperbs=[]
	for iob in range(len(obsids)) : 
		ob = obsids[iob]
		infieldccd.append('')
		infieldgalperbs.append([])
		infieldcxbperbs.append([])
		for iccd in range(len(presccdglob[iob])) :
			ccd=presccdglob[iob][iccd]

			exec 'freeze(p1_'+ob+'_'+ccd+')'
			exec 'freeze(p2_'+ob+'_'+ccd+')'
			exec 'freeze(e1_'+ob+'_'+ccd+')'
			exec 'freeze(g2_'+ob+'_'+ccd+')'
			exec 'freeze(g3_'+ob+'_'+ccd+')'
			exec 'freeze(g5_'+ob+'_'+ccd+')'

			exec 'freeze(gal1bg_'+ob+'_'+ccd+'.norm)'
			exec 'freeze(cxbbg_'+ob+'_'+ccd+'.ampl)'

			exec 'freeze(g1_'+ob+'_'+ccd+'.ampl)'

			bgparfile = open(locclu+'/ktofrsp/bgparams8_'+ob+'_ccd'+ccd+'.txt','r')
			intermed=bgparfile.readline()


			#exec 'cxbbg_'+ob+'_'+ccd+'.ampl = float(bgparfile.readline())' # used to be cxb_
			#exec 'gal1bg_'+ob+'_'+ccd+'.norm= float(bgparfile.readline())' # used to be gal1_

			intermed=bgparfile.readline()
			if intermed[0:2]!='--' :
				#exec 'gal1bg_'+ob+'_'+ccd+'.norm = float(intermed) * get_data("'+ob+'.'+ccd+'").backscal'
				infieldgalperbs[iob].append( float(intermed) )
				infieldccd[iob]= infieldccd[iob]+ccd

			intermed=bgparfile.readline()
			if intermed[0:2]!='--' :
				#exec 'cxbbg_'+ob+'_'+ccd+'.ampl = float(intermed) * get_data("'+ob+'.'+ccd+'").backscal'
				infieldcxbperbs[iob].append( float(intermed) )

			exec 'p1_'+ob+'_'+ccd+'.ampl  = float(bgparfile.readline())'
			exec 'p1_'+ob+'_'+ccd+'.gamma = float(bgparfile.readline())'
			exec 'p2_'+ob+'_'+ccd+'.ampl  = float(bgparfile.readline())'
			exec 'p2_'+ob+'_'+ccd+'.gamma = float(bgparfile.readline())'
			exec 'e1_'+ob+'_'+ccd+'.ampl  = float(bgparfile.readline())'
			exec 'e1_'+ob+'_'+ccd+'.coeff = float(bgparfile.readline())'
			exec 'g1_'+ob+'_'+ccd+'.ampl  = float(bgparfile.readline())'
			exec 'g2_'+ob+'_'+ccd+'.ampl  = float(bgparfile.readline())'
			exec 'g3_'+ob+'_'+ccd+'.ampl  = float(bgparfile.readline())'
			exec 'g5_'+ob+'_'+ccd+'.ampl  = float(bgparfile.readline())'
			intermed = bgparfile.readline() # Used to be where g6.ampl was kept. not used anymore.
			exec 'bgab.nh = float(bgparfile.readline())'

			if ccd not in infieldccd[iob] : 
				exec 'gal1bg_'+ob+'_'+ccd+'.norm  = float(bgparfile.readline())'
				exec 'cxbbg_'+ob+'_'+ccd+'.ampl  = float(bgparfile.readline())'

			bgparfile.close()



	return infieldccd, infieldgalperbs, infieldcxbperbs



