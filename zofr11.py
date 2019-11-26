
# PURPOSE: compute metallicity profile

# v2: * if ARF not present, get one, and set it equal to zero
# v3: * correcting ikt vs ktish0 mistake in above version...
#     * more complex fitting for case of multiple kT-components in one Z-shell
# v4: * write reduced statistic to file
#     * reduce the kT components (for a given Z-shell) to those whose uncertainty ranges 
#       do not overlap
# v5: * perform the last fit, with the normalizations of the different kT components unlinked.
#     * change None to NaN in output files

# ktofr9: * based on zofr5.py
#         * normalizations NOT linked

# zofr6: * based on a mix of ktofr9 and zofr5
# zofr7: * ccdname='i567' instead of 'i576'
#        * find which CCD's are in an OBSID in the very begining, to set presccdglob early on and avoid
#          including data in the analysis which will not be used. This was done in ktofr11. This is done
#          by introducing whichccd2z, and moving some of of region output files from radreg3 to whichccd2z,
#          while using radreg4 instead of 3.
#        * nnecz5 instead of nnecz4
#        * write bgxs2.txt
#        * compute the backscal based on hi-energy counts using the data in ctofr/sclOB.txt
#        * restrict calculating sclbgct to the CCD's in presccdglob, instead of adding BG
#          counts from all CCD's
#        * freeze nH
#        * get temperatures (for nnec calculation, and for initial values of fit) from ktofr11_hien.txt
#        * get Z for the for first shell from ktofr11_hien.txt, unless it's zero, then use clulist_TZ.txt
#        * when linking norms with area almost equal to area of shell, use relative backscal, instead 
#          of equating them
#        * save the temperatures obtain from the hien run, not the non-hien run
#        * don't write vikt to zofr7kt_hien.txt
#        * ct <= 10 for the loop of finding nnec
#        * write bgfrac and nnec

# zofr8: * post referee comment, from 2014 paper
#        * deleted unused section calculating srchienct
#        * always 2 temperature components, no more checking for how many kT bins to determine
#          the number of spectral components for each Z bin (ie remove ktcomp)
#        * no bounds on kT set by the T errors from the kT analysis (ie change the limits on 
#          the .kT spectral parameters)
#        * the second component is called the low T component, and labeled 'lowt'
#        * [FOR NOW] the low-T component has free T value and its initial value is kT/2, 
#          where kT is the value obtained in the kT analysis
#        * in the hien scaling, only instrumental BG components are scaled by hienbgscl.
#          sky BG components are scaled by bgscl
#        * 4 different fits: 
#          1- 2 temperatures, and nH from Dickey & Lockman
#          2- 2 temperatures, and nH from LAB
#          3- 1 temperature, and nH from Dickey & Lockman
#          4- 1 temperature, and nH from LAB
#        * a results file is created for each of the above 4 fits, for metallicity values: zofr8_hien_*.txt
#        * temperature values are all saved in one file, zofr8kt_hien.txt
#        * no spectral files created, using previously created ones.
#        * no region files created either, so changed radreg4 to radreg4_noout

# zofr9: * 1T fits first, for each of nH=DL and nH=LAB
#        * use temperature and norm from 1T fit for high-T component of 2T fit
#        * fix the norm ratio between lo and hi T components, across all obs,ccd
#        * write normalizations to file
#        * plot residuals
#        * switch to radreg4simple, using simple rectangular regions for the CCDs
#        * make new spectra based on the above regions, zofr/spsimple*.fits

# zofr10: * create arf files, instead of using old ones
#         * no special treatment for nH for A2204
#         * gal1 and cxb components are obtained from in-field spec fit, when present
#         * gal1 and cxb components' normalization is scaled, instead of scaling the entire
#           component when the model is set.
#         * gal1 and cxb components are scaled by relative backscal and NOT by backscal * exposure
#           like was done before, and like is done with instrumental BG components. this finding
#           was found using the testing script calibrate_spec_scaling3.py
#         * for certain OBSIDs, scale by bgscl and not by high energy scaling to model the BG. 
#           This is because these observations are missing counts at high energy. Their spectra
#           drop to zero, above a certain energy or PHA.

# OUTPUT: ../acchif3/zofr10*.txt
# OUTPUT: zofr/*

from sherpa.astro.ui import *
from pychips.all import *

import os
import commands
import pdb

import numpy as n

import centradec
import radreg4simple3
import arf4z4
import loadbg4z8simple3
import groupzeros
import globmod
import nnecz5
import whichccd2z_loc


def zofr(clu, obsids, r1mpc, zz, nnhh, nnhhlab) :

	print
	print ' ============ ZOFR10 for ', clu
	print

	clean()
	set_stat('cstat')

	locclu = 'mfe_'+clu

	# get coords of center
	(rac, decc) = centradec.getrd(locclu)

	# Don't use these observations
	badhienobs = ['7686','7688','7689','7690','7692','7693','7694','7696','7701']

	# CCD info
	ccdname='i567'
	ccdlist=['0,1,2,3','5','6','7']

	# Plotting settings
	dataplot  = ChipsCurve()
	modelplot = ChipsHistogram()

	dataplot.symbol.style='plus'
	dataplot.symbol.size=3
	dataplot.line.style='none'
	modelplot.line.color = "red"
	modelplot.line.thickness = 3

	if not os.path.exists(locclu+'/zofr') : os.mkdir(locclu+'/zofr')
	if not os.path.exists(locclu+'/spec') : os.mkdir(locclu+'/spec')


	if not os.path.exists(locclu+'/reg') : os.mkdir(locclu+'/reg')
	for ob in obsids :
		os.system('punlearn skyfov')
		os.system('skyfov '+locclu+'/clean'+ob+'.fits '+locclu+'/reg/simplechipsreg'+ob+'.fits')

	# which CCD's are present in each OBSID.. gets the CCD's from the simplechipsreg files
	# which were produced by skyfov. Although skyfov messes up in a few OBSIDs, it only
	# messes up coords of region. So presence or not of a region should be correct for 
	# all OBSIDs.
	presccdglob=[]
	for iob in range(len(obsids)) : presccdglob.append( whichccd2z_loc.find(clu, obsids[iob]) )

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

	# Remove OBSID's which do not have any good CCD's:
	emptyobs.reverse()
	for iob in emptyobs : 
		del obsids[iob]
		del presccdglob[iob]

	# Reset these global variables
	globmod.globstr1=''
	globmod.globstr2=''

	# Get the radial coord of each photon count
	r=[] # in Mpc
	for ob in obsids : 
		rfile = open(locclu+'/evtrad'+ob+'.txt', 'r')
		intermed = rfile.readline()
		intermed = rfile.readline()
		intermed = rfile.readline() # first data-containing line
		while intermed != '' and intermed != '\n' :
			r.append( 0.492 / r1mpc * float(intermed.split(' ')[2][0:-1])**0.5 )
			intermed = rfile.readline()
		rfile.close()
	r.sort()

	# read kT(r) ##########################################
	ktfile = open( locclu+'/ktofr14_hien.txt', 'r')
	intermed = ktfile.readline()
	ktm0=[]
	ktp0=[]
	ktm=[]
	ktp=[]
	r1=[] # in Mpc
	r2=[] # in Mpc
	ktofr=[]
	zktofr=[]
	vktline=[]
	ktline=-1
	while intermed != '\n' and intermed != '' :

		ktline=ktline+1

		r1.append( float( intermed.split()[0] ) )
		r2.append( float( intermed.split()[1] ) )
		ktofr.append( float(intermed.split()[2]) )
		zktofr.append( float(intermed.split()[9]) )
		if   str(ktofr[-1]).lower() == 'nan' and len(ktofr)==1 : ktofr[-1]=5.
		elif str(ktofr[-1]).lower() == 'nan' and len(ktofr)>1  : ktofr[-1]=ktofr[-2]

		rstat0 = float( intermed.split()[5] )
		ktm.append( float( intermed.split()[3] ) )
		ktp.append( float( intermed.split()[4] ) )

		if rstat0 > 2.0 :
			ktm[-1] = -ktofr[-1]
			ktp[-1] = 3.*ktofr[-1]
		else :
			if str(ktm[-1]).lower() == 'nan' : ktm[-1] =   -ktofr[-1]
			if str(ktp[-1]).lower() == 'nan' : ktp[-1] = 3.*ktofr[-1]

		vktline.append(ktline)
		intermed = ktfile.readline()

	ktfile.close()

	r1 = n.array(r1)
	r2 = n.array(r2)

	# Define rmin and rmax
	rmax = r2[-1]
	rmin = r1[0]

	# The first radius, defined by rmin
	for imin in range(len(r)) :
		if r[imin] >= rmin : break # imin is the index where r >= rmin

	# Read abundance from cluname_TZ.txt
	tzfile = open('cluname_TZ.txt','r')
	intermed='init'
	while intermed != '' :
		intermed = tzfile.readline()
		if intermed.split()[0] == clu :
			abin  = float( intermed.split()[2] )
			about = float( intermed.split()[4] )
			break
	tzfile.close()

	# Bin the radial data
	rbnd=[rmin] # the 1st radius
	ct1bin=100
	i=imin+ct1bin
	while i < len(r) :
		rbnd.append(r[i])
		i=i+ct1bin
	rbnd = n.array(rbnd)

	# Load string containing exculsion of all point sources
	excstr=[]
	for ob in obsids :
		excfile = open(locclu+'/reg/exclpt'+ob+'.reg', 'r')
		excstr.append( excfile.readline() )
		excfile.close()

	################### make radial bins which contain enough counts for a best-fit Z with 20% error ########################
	ccdindex = {'i':0, '5':1, '6':2, '7':3}
	ishell = 0
	rz1=[]
	rz2=[]
	vishell=[]

	# the below sets of vectors will contain the fitting results of fits done in different ways
	# as explains within the loop below

	vkt1=[]
	vnorm1=[]
	vabhe1=[]
	vabphe1=[]
	vabmhe1=[]
	vrstathien1=[]

	vkt2=[]
	vnorm2=[]
	vabhe2=[]
	vabphe2=[]
	vabmhe2=[]
	vrstathien2=[]

	vkt3=[]
	vnorm3=[]
	vabhe3=[]
	vabphe3=[]
	vabmhe3=[]
	vrstathien3=[]

	vkt4=[]
	vnorm4=[]
	vabhe4=[]
	vabphe4=[]
	vabmhe4=[]
	vrstathien4=[]

	ktbins =[0] # will contain the bins in the kT analysis which will be included
	            # in a given Z bin. If not enough counts are in 1 kT bin, the next
	            # one(s) will be added to form a Z bin.
	previnec=0

	while ktbins[-1] < len(ktofr) : 

		# compute the size of this radial bin
		ktbins0=[-1]
		ct = 1
		inec=0
		while ktbins != ktbins0 and ct <= 10 and inec < len(rbnd)-1 :

			print '======== loop ct =', ct

			nnec= []
			for ikt in ktbins :

				# Calculate BG fraction in the annulus
				sclbgct = 0
				srcct = 0
				for iob in range(len(obsids)) :
					ob = obsids[iob]
					annfile = open('antestZ.reg','w')
					annfile.write('annulus('+rac+','+decc+','+str(r1[ikt]*r1mpc)+'",'+str(r2[ikt]*r1mpc)+'")\n')
					annfile.write(excstr[iob])
					annfile.close()
					for ccd in presccdglob[iob] :
						sclbgct = sclbgct + float( commands.getoutput('dmlist "'+locclu+'/bg/bg'+ob+'_ccd'+ccd+'r_en.fits[sky=region(antestZ.reg)]" counts') ) / float( commands.getoutput('dmkeypar '+locclu+'/bg/bg'+ob+'_ccd'+ccd+'r_en.fits exposure echo+') ) * float( commands.getoutput('dmkeypar '+locclu+'/evt'+ob+'_b.fits exposure echo+') )
						srcct   = srcct   + float( commands.getoutput('dmlist "'+locclu+'/evt'+ob+'_b.fits[ccd_id='+ccdlist[ccdindex[ccd]]+'][sky=region(antestZ.reg)]" counts') )
					os.remove('antestZ.reg')
				bgfrac = sclbgct / srcct

				if ishell==0 :
					if zktofr[0] != 0 : abun=zktofr[0]
					elif (r1[ikt]+r2[ikt])/2. <= 0.1 or about < 0. : abun=abin
					else : abun=about
				else : abun = vabhe3[ishell-1]
				nnec.append(nnecz5.nnec(ktofr[ikt], bgfrac, abun))

			inec = min(   int( max(nnec) / ct1bin) +previnec +1  ,  len(rbnd)-1  )
			rr1 = r1[ktbins[0]]
			rr2 = rbnd[inec]

			ktbins0=ktbins
			ktbins=[]
			for ikt in range(len(r1)) :
				if r1[ikt] >= rr1 and r1[ikt] < rr2 : ktbins.append(ikt)

			ct=ct+1

		rz1.append( r1[ktbins[0]]  )
		rz2.append( r2[ktbins[-1]] )


		################# fit this radial bin ####################################

		ikt=ktbins[0]

		# which CCD's have data in this ishell
		presccd=[]
		for iob in range(len(obsids)) :
			ob = obsids[iob]
			# radreg4_noout returns the ccd's that are present in this SHELL!
			intermed1 = radreg4simple3.radreg(clu, ob, rz1[ishell]*r1mpc/.492, rz2[ishell]*r1mpc/.492, locclu+'/zofr/zan11_'+ob+'_'+str(ishell))
			presccd.append( intermed1 )

			for ccd in str(presccd[iob]) :

				# Make the spectra for SHELL x CCD
				os.system('punlearn dmextract')
				os.system('dmextract "'+locclu+'/clean'+ob+'.fits[ccd_id='+ccdlist[ccdindex[ccd]]+'][sky=region('+locclu+'/zofr/zan11_'+ob+'_'+str(ishell)+'_ccd'+ccd+'_xfov_pt_simple.reg)][bin PI]" '+locclu+'/zofr/sp10simple'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits opt=pha1 wmap="[energy=300:2000][bin det=8]" clobber=yes' )

				# Make ARF
				arf4z4.arf(clu, ob, ishell, ccd)

		# Make BG spectra and load them, only one for ishell=0
		# This is here because it uses zofr/sp10simple*.reg, created just
		# above in radreg4simple.radreg()
		if ishell==0 : 

			infieldccd, infieldgalperbs, infieldcxbperbs = loadbg4z8simple3.bg(clu, obsids, presccdglob)

			for iob in range(len(obsids)) :
				ob = obsids[iob]
				for ccd in str(presccdglob[iob]) : 
					exec 'gal1bg_'+ob+'_'+ccd+'_normval = float( gal1bg_'+ob+'_'+ccd+'.norm.val )'
					exec 'cxbbg_'+ob+'_'+ccd+'_amplval = float( cxbbg_'+ob+'_'+ccd+'.ampl.val )'

			bginf_present=False
			iob_in=-1
			for iob in range(len(obsids)) : 
				if len(infieldccd[iob])>0 : 
					bginf_present=True
					iob_in = iob

		firstiob=-1
		for iob in range(len(obsids)) :
			ob = obsids[iob]
			for ccd in str(presccd[iob]) :
				if not os.path.exists(locclu+'/zofr/sp10simple'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits') or not os.path.exists(locclu+'/zofr/warf10_'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits') or not os.path.exists(locclu+'/spec/sp'+ob+'_ccd'+ccd+'_center.wrmf') or ccd not in presccdglob[iob] : presccd[iob] = presccd[iob].replace(ccd, '')
			if len(presccd[iob]) > 0 and firstiob==-1 : firstiob=iob

		# Load spectral data 
		srcstr=''
		datasetct=[]
		datasetname=[]
		for iob in range(len(obsids)) :
			ob = obsids[iob]
			for ccd in presccd[iob] :
				# Load data
				load_pha( ob+'.'+str(ishell)+'.'+ccd , locclu+'/zofr/sp10simple'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits' )
				# Load RMF
				load_rmf( ob+'.'+str(ishell)+'.'+ccd , locclu+'/spec/sp'+ob+'_ccd'+ccd+'_center.wrmf'  )
				# Loading ARF:
				load_arf( ob+'.'+str(ishell)+'.'+ccd , locclu+'/zofr/warf10_'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits' )
				# string containing the dataset names
				srcstr = srcstr +'"'+ob+'.'+str(ishell)+'.'+ccd + '",'




		# If there are datasets to be fit FOR THIS SHELL
		if srcstr[0:-1] != '' :

			print
			print 'start fitting ---'
			print 'initial settings for all fits ---'
			print

###############################################
######### INITIAL SETTINGS FOR ALL FITS (start)
			iobiccd1=[-1,-1] # the first (obs,ccd) found to be within 5% of annulus area. used to set initial normalization
			anarea = n.pi * ( (rz2[ishell]*r1mpc/.492)**2. - (rz1[ishell]*r1mpc/.492)**2. ) # used to set initial normalization
			for iob in range(len(obsids)) : 
				ob=obsids[iob]
				for ccd in presccd[iob] :

					newhienbgscl = calc_data_sum(id=ob+'.'+str(ishell)+'.'+ccd, lo=9.5, hi=12) / calc_data_sum(id='bg'+ob+'.'+ccd, lo=9.5, hi=12)
					bgscl = get_data(ob+'.'+str(ishell)+'.'+ccd).backscal * get_data(ob+'.'+str(ishell)+'.'+ccd).exposure / get_data('bg'+ob+'.'+ccd).backscal / get_data('bg'+ob+'.'+ccd).exposure

					#8 Define a model with 2 temperatures, this is just a string that will be used
					#8 to set the model
					src2 = 'xsapec.em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ktbins[0])+' + xsapec.em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ktbins[0])+'lowt'

					exec 'delete_model("'+ob+'.'+str(ishell)+'.'+ccd+'")'
					rsp   = get_response(  ob+'.'+str(ishell)+'.'+ccd )
					bgrsp = get_response( 'bg'+ob+'.'+ccd )
					groupzeros.grp( ob+'.'+str(ishell)+'.'+ccd, ':.3,7:')

					#notice_id( ob+'.'+str(ishell)+'.'+ccd )
					#group_adapt( ob+'.'+str(ishell)+'.'+ccd,  20)
					#ignore_id( ob+'.'+str(ishell)+'.'+ccd ,':.3,7:' )

					if ob not in badhienobs :
						exec 'set_full_model( "'+ob+'.'+str(ishell)+'.'+ccd+'", rsp(xsphabs.ab * ('+src2+') ) + newhienbgscl * bgrsp( powlaw1d.p1_'+ob+'_'+ccd+' + powlaw1d.p2_'+ob+'_'+ccd+' + exp.e1_'+ob+'_'+ccd+' + gauss1d.g1_'+ob+'_'+ccd+' + gauss1d.g2_'+ob+'_'+ccd+' + gauss1d.g3_'+ob+'_'+ccd+' + gauss1d.g5_'+ob+'_'+ccd+' ) + rsp( xsphabs.ab * powlaw1d.cxbbg_'+ob+'_'+ccd+' + xsapec.gal1bg_'+ob+'_'+ccd+' ) )'
					else : 
						exec 'set_full_model( "'+ob+'.'+str(ishell)+'.'+ccd+'", rsp(xsphabs.ab * ('+src2+') ) + bgscl * bgrsp( powlaw1d.p1_'+ob+'_'+ccd+' + powlaw1d.p2_'+ob+'_'+ccd+' + exp.e1_'+ob+'_'+ccd+' + gauss1d.g1_'+ob+'_'+ccd+' + gauss1d.g2_'+ob+'_'+ccd+' + gauss1d.g3_'+ob+'_'+ccd+' + gauss1d.g5_'+ob+'_'+ccd+' ) + rsp( xsphabs.ab * powlaw1d.cxbbg_'+ob+'_'+ccd+' + xsapec.gal1bg_'+ob+'_'+ccd+' ) )'


					# If we fit the in-field BG, then obtain gal1 and cxb components from it
					if bginf_present : 
						exec 'set_par( gal1bg_'+ob+'_'+ccd+'.norm, val=infieldgalperbs[iob_in][0]*get_data("'+ob+'.'+str(ishell)+'.'+ccd+'").backscal )'
						exec 'set_par(  cxbbg_'+ob+'_'+ccd+'.ampl, val=infieldcxbperbs[iob_in][0]*get_data("'+ob+'.'+str(ishell)+'.'+ccd+'").backscal )'

					# if no in-field BG, then scale down the values of gal1 and cxb, which were already
					# read by loadbg4z8simple. scaling is according to ratio of backscal of source spec
					# to backscal of BG spec
					else :
						exec 'set_par( gal1bg_'+ob+'_'+ccd+'.norm, val=get_data("'+ob+'.'+str(ishell)+'.'+ccd+'").backscal * gal1bg_'+ob+'_'+ccd+'_normval )'
						exec 'set_par(  cxbbg_'+ob+'_'+ccd+'.ampl, val=get_data("'+ob+'.'+str(ishell)+'.'+ccd+'").backscal * cxbbg_'+ob+'_'+ccd+'_amplval )'




				# Set redshift - fixed
				for ccd in presccd[iob] :
					exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ktbins[0])+'.redshift, val='+str(zz)+', frozen=True)'
					exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ktbins[0])+'lowt.redshift, val='+str(zz)+', frozen=True)'

				# Set normalization LINKING for high-T components (for when backscal covers 95-105 % of annulus area)
				for iccd in range(len(presccd[iob])) :
					ccd=presccd[iob][iccd]
					bspix2 = float( commands.getoutput('dmkeypar '+locclu+'/zofr/sp10simple'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits backscal echo+') ) *64.0*1024.**2.

					if iobiccd1==[-1,-1] or bspix2/anarea < 0.95 or bspix2/anarea > 1.05 :
						if iobiccd1==[-1,-1] and bspix2/anarea >= 0.95 and bspix2/anarea <= 1.05 : iobiccd1=[iob,iccd]
					else : 
						relbs = bspix2 / float( commands.getoutput('dmkeypar '+locclu+'/zofr/sp10simple'+obsids[iobiccd1[0]]+'_'+str(ishell)+'_ccd'+presccd[iobiccd1[0]][iobiccd1[1]]+'.fits backscal echo+') ) /(64.0*1024.**2.)
						exec 'link(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ktbins[0])+'.norm, em'+str(ishell)+'_'+obsids[iobiccd1[0]]+'_'+presccd[iobiccd1[0]][iobiccd1[1]]+'_'+str(ktbins[0])+'.norm * relbs )'



######### INITIAL SETTINGS FOR ALL FITS (end)
###############################################


###############################################
############# Fit 3 : 1T, nH=DL

			for iob in range(len(obsids)) : 
				ob=obsids[iob]

				for iccd in range(len(presccd[iob])) :
					ccd=presccd[iob][iccd]

					# Freeze low-T component, set its norm=0
					exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'lowt.norm, val=0.0, min=0., max=2., frozen=True)'
					exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'lowt.kt, frozen=True)'
					exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'lowt.abundanc, frozen=True)'

					# abundance, initial value
					if ishell != 0 and n.isfinite(vabhe3[ishell-1]) : 
						exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, val='+str(vabhe3[ishell-1])+', min=0, max=5, frozen=False)'
					elif zktofr[0] != 0 :
						exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, val='+str(zktofr[0])+', min=0, max=5, frozen=False)'
					elif (r1[ikt]+r2[ikt])/2. <= 0.1 or about < 0. : 
						exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, val='+str(abin)+', min=0, max=5, frozen=False)'
					else : 
						exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, val='+str(about)+', min=0, max=5, frozen=False)'
					if iob != firstiob or ccd != presccd[firstiob][0] : exec 'link( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'.abundanc )'

					# temperature, initial value
					if iob == firstiob and ccd == presccd[firstiob][0] : exec 'set_par(em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ikt)+'.kt, val=ktofr[ikt], min=0., max=40., frozen=False )'
					else : exec 'link(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.kt, em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ikt)+'.kt)'

					# normalization initial values, for the spectra whose norm is NOT linked
					bspix2 = float( commands.getoutput('dmkeypar '+locclu+'/zofr/sp10simple'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits backscal echo+') ) *64.0*1024.**2.
					if iobiccd1==[iob,iccd] or bspix2/anarea < 0.95 or bspix2/anarea > 1.05 : 
						exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.norm, val=0.02, min=2e-4, max=2.)'
						exec 'emnorm0=abs(float(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.norm.val)*calc_data_sum(id="'+ob+'.'+str(ishell)+'.'+ccd+'",lo=0.3,hi=7.0)/calc_model_sum(id="'+ob+'.'+str(ishell)+'.'+ccd+'",lo=0.3,hi=7.0))'
						exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.norm, val=float(emnorm0), min=float(emnorm0/1e3), max=float(emnorm0*1e3) )'

			# Set nH - 3
			if   clu == 'a478'  : set_par(ab.nh, val=nnhh, min=nnhh, max=nnhh*100., frozen=False)
			#elif clu == 'a2204' : set_par(ab.nh, val=nnhh, min=nnhh/100, max=nnhh*100., frozen=False)
			else                : set_par(ab.nh, val=nnhh, min=nnhh/10., max=nnhh*10., frozen=True) # <<< DEFAULT

			print
			print 'BEFORE FIT 3'
			exec 'print em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])
			exec 'print em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'lowt'
			print

			# Fitting 3 - nH(DL), 1T
			exec 'fit('+srcstr[0:-1]+')'

			# Plotting 3 - nH(DL), 1T
			if ishell != 0 : delete_window('all')
			x1=[0.,0.5,0.,0.5]
			y1=[0.5,0.5,0.,0.]
			for iob in range(len(obsids)) : 
				ob=obsids[iob]
				if len(presccd[iob]) > 0 : add_window(8.5,8.5,'inches')
				for iccd in range(len(presccd[iob])) :
					ccd=presccd[iob][iccd]
					add_frame(x1[iccd], y1[iccd], x1[iccd]+0.5, y1[iccd]+0.5)
					add_curve( get_data_plot(ob+'.'+str(ishell)+'.'+ccd).x, get_data_plot(ob+'.'+str(ishell)+'.'+ccd).y, dataplot )
					add_histogram( get_model_plot(ob+'.'+str(ishell)+'.'+ccd).xlo, get_model_plot(ob+'.'+str(ishell)+'.'+ccd).xhi, get_model_plot(ob+'.'+str(ishell)+'.'+ccd).y, modelplot )
					set_plot_title(clu+' OB'+ob+' CCD'+ccd)
					split(2)
					add_curve( get_resid_plot(ob+'.'+str(ishell)+'.'+ccd).x, get_resid_plot(ob+'.'+str(ishell)+'.'+ccd).y, dataplot )
					add_hline(0)
				if os.path.exists('fitplots/zofr10_'+clu+'_'+ob+'_'+str(ishell)+'_hien_1T.eps') : os.remove('fitplots/zofr10_'+clu+'_'+ob+'_'+str(ishell)+'_hien_1T.eps')
				if len(presccd[iob]) > 0 : print_window('fitplots/zofr10_'+clu+'_'+ob+'_'+str(ishell)+'_hien_1T', ['format', 'eps', 'orientation', 'landscape'])

			# HIEN results, best-fit, errors 3 - nH(DL), 1T
			vrstathien3.append( get_fit_results().rstat )
			if get_fit_results().rstat < 3 : 
				exec 'proj('+srcstr+' em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'.abundanc)'
				vabhe3.append( get_proj_results().parvals[0] )
				vabphe3.append( get_proj_results().parmaxes[0] )
				vabmhe3.append( get_proj_results().parmins[0] )
				exec 'vkt3.append(em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'.kt.val)'
				exec 'vnorm3.append(em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'.norm.val)'
			else :
				vabhe3.append( n.nan )
				vkt3.append( n.nan )
				vnorm3.append( n.nan )
				vabphe3.append( n.nan )
				vabmhe3.append( n.nan )

############# Fit 3 : 1T, nH=DL (end)
###############################################





###############################################
############# Fit 1 : 2T, nH=DL



			for iob in range(len(obsids)) : 
				ob=obsids[iob]

				for iccd in range(len(presccd[iob])) :
					ccd=presccd[iob][iccd]

					# high-T abundance, initial value
					if ishell != 0 and n.isfinite(vabhe1[ishell-1]) : 
						exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, val='+str(vabhe1[ishell-1])+', min=0, max=5, frozen=False)'
					elif zktofr[0] != 0 :
						exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, val='+str(zktofr[0])+', min=0, max=5, frozen=False)'
					elif (r1[ikt]+r2[ikt])/2. <= 0.1 or about < 0. : 
						exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, val='+str(abin)+', min=0, max=5, frozen=False)'
					else : 
						exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, val='+str(about)+', min=0, max=5, frozen=False)'
					if iob != firstiob or ccd != presccd[firstiob][0] : exec 'link( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'.abundanc )'

					# high-T temperature, initial value = the fit from the 1T model above, ie vkt3
					if   iob == firstiob and ccd == presccd[firstiob][0] and     n.isfinite(vkt3[ishell]) : exec 'set_par(em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ikt)+'.kt, val=vkt3[ishell], min=0., max=40., frozen=False )'
					elif iob == firstiob and ccd == presccd[firstiob][0] and not n.isfinite(vkt3[ishell]) : exec 'set_par(em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ikt)+'.kt, val=ktofr[ikt], min=0., max=40., frozen=False )'
					else : exec 'link(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.kt, em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ikt)+'.kt)'

					# low-T abundance: always linked to the hi-T component
					exec 'link( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'lowt.abundanc, em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'.abundanc )'


					# low-T temperature linked as kt_hi / 2
					exec 'link(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'lowt.kt, em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ikt)+'.kt / 2)'


					# high-T normalization, start w value from 1T fit
					bspix2 = float( commands.getoutput('dmkeypar '+locclu+'/zofr/sp10simple'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits backscal echo+') ) *64.0*1024.**2.
					if iobiccd1==[iob,iccd] or bspix2/anarea < 0.95 or bspix2/anarea > 1.05 :
						#exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.norm, val=0.02, min=2e-4, max=2.)'
						exec 'emnorm0=float(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.norm.val)' # from 1T fit
						exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.norm, val=float(emnorm0), min=float(emnorm0/1e3), max=float(emnorm0*1e3) )'
					else : 
						relbs = bspix2 / float( commands.getoutput('dmkeypar '+locclu+'/zofr/sp10simple'+obsids[iobiccd1[0]]+'_'+str(ishell)+'_ccd'+presccd[iobiccd1[0]][iobiccd1[1]]+'.fits backscal echo+') ) /(64.0*1024.**2.)

					# low-T normalization
					if   iob==firstiob and ccd == presccd[firstiob][0] and     (iobiccd1==[iob,iccd] or bspix2/anarea < 0.95 or bspix2/anarea > 1.05) :
						exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'lowt.norm, val=float(emnorm0/3.), min=0., max=float(emnorm0*1e3), frozen=False )'
					elif iob==firstiob and ccd == presccd[firstiob][0] and not (iobiccd1==[iob,iccd] or bspix2/anarea < 0.95 or bspix2/anarea > 1.05) :
						exec 'link(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'lowt.norm, em'+str(ishell)+'_'+obsids[iobiccd1[0]]+'_'+presccd[iobiccd1[0]][iobiccd1[1]]+'_'+str(ikt)+'lowt.norm * relbs )'
					else :
						# the ratio between hi and lo components should be the same for ea obs,ccd
						exec 'link(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'lowt.norm, ( em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ikt)+'lowt.norm / em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ikt)+'.norm ) * em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.norm)'


			# Set nH - 1
			if   clu == 'a478'  : set_par(ab.nh, val=nnhh, min=nnhh, max=nnhh*100., frozen=False)
			else                : set_par(ab.nh, val=nnhh, min=nnhh/10., max=nnhh*10., frozen=True) # <<< DEFAULT

			print
			print 'BEFORE FIT 1'
			exec 'print em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])
			exec 'print em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'lowt'
			print

			# Fitting 1 - nH(DL), 2T
			exec 'fit('+srcstr[0:-1]+')'

			# Plotting 1 - nH(DL), 2T
			delete_window('all')
			x1=[0.,0.5,0.,0.5]
			y1=[0.5,0.5,0.,0.]
			for iob in range(len(obsids)) : 
				ob=obsids[iob]
				if len(presccd[iob]) > 0 : add_window(8.5,8.5,'inches')
				for iccd in range(len(presccd[iob])) :
					ccd=presccd[iob][iccd]
					add_frame(x1[iccd], y1[iccd], x1[iccd]+0.5, y1[iccd]+0.5)
					add_curve( get_data_plot(ob+'.'+str(ishell)+'.'+ccd).x, get_data_plot(ob+'.'+str(ishell)+'.'+ccd).y, dataplot )
					add_histogram( get_model_plot(ob+'.'+str(ishell)+'.'+ccd).xlo, get_model_plot(ob+'.'+str(ishell)+'.'+ccd).xhi, get_model_plot(ob+'.'+str(ishell)+'.'+ccd).y, modelplot )
					set_plot_title(clu+' OB'+ob+' CCD'+ccd)
					split(2)
					add_curve( get_resid_plot(ob+'.'+str(ishell)+'.'+ccd).x, get_resid_plot(ob+'.'+str(ishell)+'.'+ccd).y, dataplot )
					add_hline(0)
				if os.path.exists('fitplots/zofr10_'+clu+'_'+ob+'_'+str(ishell)+'_hien_2T.eps') : os.remove('fitplots/zofr10_'+clu+'_'+ob+'_'+str(ishell)+'_hien_2T.eps')
				if len(presccd[iob]) > 0 : print_window('fitplots/zofr10_'+clu+'_'+ob+'_'+str(ishell)+'_hien_2T', ['format', 'eps', 'orientation', 'landscape'])

			# HIEN results, best-fit, errors 1 - nH(DL), 2T
			vrstathien1.append( get_fit_results().rstat )
			if get_fit_results().rstat < 3 : 
				exec 'proj('+srcstr+' em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'.abundanc)'
				vabhe1.append( get_proj_results().parvals[0] )
				vabphe1.append( get_proj_results().parmaxes[0] )
				vabmhe1.append( get_proj_results().parmins[0] )
				exec 'vkt1.append(em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'.kt.val)'
				exec 'vkt1.append(em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'lowt.kt.val)'
				exec 'vnorm1.append(em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'.norm.val)'
				exec 'vnorm1.append(em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'lowt.norm.val)'
			else :
				vabhe1.append( n.nan )
				vkt1.append( n.nan )
				vkt1.append( n.nan )
				vnorm1.append( n.nan )
				vnorm1.append( n.nan )
				vabphe1.append( n.nan )
				vabmhe1.append( n.nan )


############# Fit 1 : 2T, nH=DL (end)
###############################################







###############################################
############# Fit 4 : 1T, nH=LAB

			for iob in range(len(obsids)) : 
				ob=obsids[iob]

				for iccd in range(len(presccd[iob])) :
					ccd=presccd[iob][iccd]

					# Freeze low-T component, set its norm=0
					exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'lowt.norm, val=0.0, min=0., max=2., frozen=True)'
					exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'lowt.kt, frozen=True)'
					exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'lowt.abundanc, frozen=True)'

					# abundance, initial value
					if ishell != 0 and n.isfinite(vabhe4[ishell-1]) : 
						exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, val='+str(vabhe4[ishell-1])+', min=0, max=5, frozen=False)'
					elif zktofr[0] != 0 :
						exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, val='+str(zktofr[0])+', min=0, max=5, frozen=False)'
					elif (r1[ikt]+r2[ikt])/2. <= 0.1 or about < 0. : 
						exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, val='+str(abin)+', min=0, max=5, frozen=False)'
					else : 
						exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, val='+str(about)+', min=0, max=5, frozen=False)'
					if iob != firstiob or ccd != presccd[firstiob][0] : exec 'link( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'.abundanc )'

					# temperature, initial value
					if iob == firstiob and ccd == presccd[firstiob][0] : exec 'set_par(em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ikt)+'.kt, val=ktofr[ikt], min=0., max=40., frozen=False )'
					else : exec 'link(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.kt, em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ikt)+'.kt)'

					# normalization initial values, for the spectra whose norm is NOT linked
					bspix2 = float( commands.getoutput('dmkeypar '+locclu+'/zofr/sp10simple'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits backscal echo+') ) *64.0*1024.**2.
					if iobiccd1==[iob,iccd] or bspix2/anarea < 0.95 or bspix2/anarea > 1.05 : 
						exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.norm, val=0.02, min=2e-4, max=2.)'
						exec 'emnorm0=abs(float(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.norm.val)*calc_data_sum(id="'+ob+'.'+str(ishell)+'.'+ccd+'",lo=0.3,hi=7.0)/calc_model_sum(id="'+ob+'.'+str(ishell)+'.'+ccd+'",lo=0.3,hi=7.0))'
						exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.norm, val=float(emnorm0), min=float(emnorm0/1e3), max=float(emnorm0*1e3) )'

			# Set nH - 4
			if   clu == 'a478'  : set_par(ab.nh, val=nnhhlab, min=nnhhlab, max=nnhhlab*100., frozen=False)
			else                : set_par(ab.nh, val=nnhhlab, min=nnhhlab/10., max=nnhhlab*10., frozen=True) # <<< DEFAULT

			print
			print 'BEFORE FIT 4'
			exec 'print em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])
			exec 'print em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'lowt'
			print

			# Fitting 4 - nH(LAB), 1T
			exec 'fit('+srcstr[0:-1]+')'

			# Plotting 4 - nH(LAB), 1T
			delete_window('all')
			x1=[0.,0.5,0.,0.5]
			y1=[0.5,0.5,0.,0.]
			for iob in range(len(obsids)) : 
				ob=obsids[iob]
				if len(presccd[iob]) > 0 : add_window(8.5,8.5,'inches')
				for iccd in range(len(presccd[iob])) :
					ccd=presccd[iob][iccd]
					add_frame(x1[iccd], y1[iccd], x1[iccd]+0.5, y1[iccd]+0.5)
					add_curve( get_data_plot(ob+'.'+str(ishell)+'.'+ccd).x, get_data_plot(ob+'.'+str(ishell)+'.'+ccd).y, dataplot )
					add_histogram( get_model_plot(ob+'.'+str(ishell)+'.'+ccd).xlo, get_model_plot(ob+'.'+str(ishell)+'.'+ccd).xhi, get_model_plot(ob+'.'+str(ishell)+'.'+ccd).y, modelplot )
					set_plot_title(clu+' OB'+ob+' CCD'+ccd)
					split(2)
					add_curve( get_resid_plot(ob+'.'+str(ishell)+'.'+ccd).x, get_resid_plot(ob+'.'+str(ishell)+'.'+ccd).y, dataplot )
					add_hline(0)
				if os.path.exists('fitplots/zofr10_'+clu+'_'+ob+'_'+str(ishell)+'_hien_1T_lab.eps') : os.remove('fitplots/zofr10_'+clu+'_'+ob+'_'+str(ishell)+'_hien_1T_lab.eps')
				if len(presccd[iob]) > 0 : print_window('fitplots/zofr10_'+clu+'_'+ob+'_'+str(ishell)+'_hien_1T_lab', ['format', 'eps', 'orientation', 'landscape'])

			# HIEN results, best-fit, errors 4 - nH(LAB), 1T
			vrstathien4.append( get_fit_results().rstat )
			if get_fit_results().rstat < 3 : 
				exec 'proj('+srcstr+' em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'.abundanc)'
				vabhe4.append( get_proj_results().parvals[0] )
				vabphe4.append( get_proj_results().parmaxes[0] )
				vabmhe4.append( get_proj_results().parmins[0] )
				exec 'vkt4.append(em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'.kt.val)'
				exec 'vnorm4.append(em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'.norm.val)'
			else :
				vabhe4.append( n.nan )
				vkt4.append( n.nan )
				vnorm4.append( n.nan )
				vabphe4.append( n.nan )
				vabmhe4.append( n.nan )

############# Fit 4 : 1T, nH=LAB (end)
###############################################











###############################################
############# Fit 2 : 2T, nH=LAB



			for iob in range(len(obsids)) : 
				ob=obsids[iob]

				for iccd in range(len(presccd[iob])) :
					ccd=presccd[iob][iccd]

					# high-T abundance, initial value
					if ishell != 0 and n.isfinite(vabhe2[ishell-1]) : 
						exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, val='+str(vabhe2[ishell-1])+', min=0, max=5, frozen=False)'
					elif zktofr[0] != 0 :
						exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, val='+str(zktofr[0])+', min=0, max=5, frozen=False)'
					elif (r1[ikt]+r2[ikt])/2. <= 0.1 or about < 0. : 
						exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, val='+str(abin)+', min=0, max=5, frozen=False)'
					else : 
						exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, val='+str(about)+', min=0, max=5, frozen=False)'
					if iob != firstiob or ccd != presccd[firstiob][0] : exec 'link( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.abundanc, em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'.abundanc )'

					# high-T temperature, initial value = the fit from the 1T model, ie vkt4
					if   iob == firstiob and ccd == presccd[firstiob][0] and     n.isfinite(vkt4[ishell]) : exec 'set_par(em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ikt)+'.kt, val=vkt4[ishell], min=0., max=40., frozen=False )'
					elif iob == firstiob and ccd == presccd[firstiob][0] and not n.isfinite(vkt4[ishell]) : exec 'set_par(em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ikt)+'.kt, val=ktofr[ikt], min=0., max=40., frozen=False )'
					else : exec 'link(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.kt, em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ikt)+'.kt)'

					# low-T abundance: always linked to the hi-T component
					exec 'link( em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'lowt.abundanc, em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'.abundanc )'

					# low-T temperature linked as kt_hi / 2
					exec 'link(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'lowt.kt, em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ikt)+'.kt / 2)'

					# high-T normalization, initial value from the 1T fit
					bspix2 = float( commands.getoutput('dmkeypar '+locclu+'/zofr/sp10simple'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits backscal echo+') ) *64.0*1024.**2.
					if iobiccd1==[iob,iccd] or bspix2/anarea < 0.95 or bspix2/anarea > 1.05 :
						#exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.norm, val=0.02, min=2e-4, max=2.)'
						exec 'emnorm0=float(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.norm.val)' # from the 1T fit
						exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.norm, val=float(emnorm0), min=float(emnorm0/1e3), max=float(emnorm0*1e3) )'
					else : 
						relbs = bspix2 / float( commands.getoutput('dmkeypar '+locclu+'/zofr/sp10simple'+obsids[iobiccd1[0]]+'_'+str(ishell)+'_ccd'+presccd[iobiccd1[0]][iobiccd1[1]]+'.fits backscal echo+') ) /(64.0*1024.**2.)

					# low-T normalization
					if   iob==firstiob and ccd == presccd[firstiob][0] and     (iobiccd1==[iob,iccd] or bspix2/anarea < 0.95 or bspix2/anarea > 1.05) :
						exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'lowt.norm, val=float(emnorm0/3.), min=float(emnorm0/1e3), max=float(emnorm0*1e3), frozen=False )'
					elif iob==firstiob and ccd == presccd[firstiob][0] and not (iobiccd1==[iob,iccd] or bspix2/anarea < 0.95 or bspix2/anarea > 1.05) :
						exec 'link(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'lowt.norm, em'+str(ishell)+'_'+obsids[iobiccd1[0]]+'_'+presccd[iobiccd1[0]][iobiccd1[1]]+'_'+str(ikt)+'lowt.norm * relbs )'
					else :
						# the ratio between hi and lo components should be the same for ea obs,ccd
						exec 'link(em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'lowt.norm, ( em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ikt)+'lowt.norm / em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ikt)+'.norm ) * em'+str(ishell)+'_'+ob+'_'+ccd+'_'+str(ikt)+'.norm)'

			# Set nH - 2
			if   clu == 'a478'  : set_par(ab.nh, val=nnhhlab, min=nnhhlab, max=nnhhlab*100., frozen=False)
			else                : set_par(ab.nh, val=nnhhlab, min=nnhhlab/10., max=nnhhlab*10., frozen=True) # <<< DEFAULT

			print
			print 'BEFORE FIT 2'
			exec 'print em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])
			exec 'print em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'lowt'
			print

			# Fitting 2 - nH(LAB), 2T
			exec 'fit('+srcstr[0:-1]+')'

			# Plotting 2 - nH(LAB), 2T
			delete_window('all')
			x1=[0.,0.5,0.,0.5]
			y1=[0.5,0.5,0.,0.]
			for iob in range(len(obsids)) : 
				ob=obsids[iob]
				if len(presccd[iob]) > 0 : add_window(8.5,8.5,'inches')
				for iccd in range(len(presccd[iob])) :
					ccd=presccd[iob][iccd]
					add_frame(x1[iccd], y1[iccd], x1[iccd]+0.5, y1[iccd]+0.5)
					add_curve( get_data_plot(ob+'.'+str(ishell)+'.'+ccd).x, get_data_plot(ob+'.'+str(ishell)+'.'+ccd).y, dataplot )
					add_histogram( get_model_plot(ob+'.'+str(ishell)+'.'+ccd).xlo, get_model_plot(ob+'.'+str(ishell)+'.'+ccd).xhi, get_model_plot(ob+'.'+str(ishell)+'.'+ccd).y, modelplot )
					set_plot_title(clu+' OB'+ob+' CCD'+ccd)
					split(2)
					add_curve( get_resid_plot(ob+'.'+str(ishell)+'.'+ccd).x, get_resid_plot(ob+'.'+str(ishell)+'.'+ccd).y, dataplot )
					add_hline(0)
				if os.path.exists('fitplots/zofr10_'+clu+'_'+ob+'_'+str(ishell)+'_hien_2T_lab.eps') : os.remove('fitplots/zofr10_'+clu+'_'+ob+'_'+str(ishell)+'_hien_2T_lab.eps')
				if len(presccd[iob]) > 0 : print_window('fitplots/zofr10_'+clu+'_'+ob+'_'+str(ishell)+'_hien_2T_lab', ['format', 'eps', 'orientation', 'landscape'])

			# HIEN results, best-fit, errors 2 - nH(DL), 2T
			vrstathien2.append( get_fit_results().rstat )
			if get_fit_results().rstat < 3 : 
				exec 'proj('+srcstr+' em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'.abundanc)'
				vabhe2.append( get_proj_results().parvals[0] )
				vabphe2.append( get_proj_results().parmaxes[0] )
				vabmhe2.append( get_proj_results().parmins[0] )
				exec 'vkt2.append(em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'.kt.val)'
				exec 'vkt2.append(em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'lowt.kt.val)'
				exec 'vnorm2.append(em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'.norm.val)'
				exec 'vnorm2.append(em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'_'+str(ktbins[0])+'lowt.norm.val)'
			else :
				vabhe2.append( n.nan )
				vkt2.append( n.nan )
				vkt2.append( n.nan )
				vnorm2.append( n.nan )
				vnorm2.append( n.nan )
				vabphe2.append( n.nan )
				vabmhe2.append( n.nan )

############# Fit 2 : 2T, nH=LAB (end)
###############################################

			vishell.append(ishell)


			globmod.globstr1 = globmod.globstr1 + str(rz1[ishell]) +' '+ str(rz2[ishell]) +' '+ str(vabhe1[ishell]) +' '+ str(vabmhe1[ishell]) +' '+ str(vabphe1[ishell]) + ' ' + str(vrstathien1[ishell]) + ' ' + str(bgfrac) + ' ' + str(max(nnec)) +'\n'
			globmod.globstr2 = globmod.globstr2 + str(rz1[ishell]) +' '+ str(rz2[ishell]) +' '+ str(vabhe2[ishell]) +' '+ str(vabmhe2[ishell]) +' '+ str(vabphe2[ishell]) + ' ' + str(vrstathien2[ishell]) + ' ' + str(bgfrac) + ' ' + str(max(nnec)) +'\n'
			globmod.globstr3 = globmod.globstr3 + str(rz1[ishell]) +' '+ str(rz2[ishell]) +' '+ str(vabhe3[ishell]) +' '+ str(vabmhe3[ishell]) +' '+ str(vabphe3[ishell]) + ' ' + str(vrstathien3[ishell]) + ' ' + str(bgfrac) + ' ' + str(max(nnec)) +'\n'
			globmod.globstr4 = globmod.globstr4 + str(rz1[ishell]) +' '+ str(rz2[ishell]) +' '+ str(vabhe4[ishell]) +' '+ str(vabmhe4[ishell]) +' '+ str(vabphe4[ishell]) + ' ' + str(vrstathien4[ishell]) + ' ' + str(bgfrac) + ' ' + str(max(nnec)) +'\n'




		# end of **ishell** loop variable change:
		if len(n.where(rbnd > rz2[ishell])[0]) == 0 : break
		else : previnec = max( int(inec), n.where(rbnd > rz2[ishell])[0][0] )
		ishell=ishell+1
		ktbins=[ktbins[-1]+1]





	# save the results
	abfile=open(locclu+'/zofr10_hien_2T.txt', 'w')
	abfile.write(globmod.globstr1)
	abfile.close()

	abfile=open(locclu+'/zofr10_hien_2T_lab.txt', 'w')
	abfile.write(globmod.globstr2)
	abfile.close()

	abfile=open(locclu+'/zofr10_hien_1T.txt', 'w')
	abfile.write(globmod.globstr3)
	abfile.close()

	abfile=open(locclu+'/zofr10_hien_1T_lab.txt', 'w')
	abfile.write(globmod.globstr4)
	abfile.close()


	ktoutfile=open(locclu+'/zofr10kt_hien_2T.txt', 'w')
	for ikt in range(len(vishell)) :
		#same order as analysis:
		ktoutfile.write(str(vishell[ikt]) +' '+ str(vkt3[ikt]) +' '+ str(vnorm3[ikt]) +' '+ str(vkt1[2*ikt]) +' '+ str(vnorm1[2*ikt]) +' '+ str(vkt1[2*ikt+1]) +' '+ str(vnorm1[2*ikt+1]) +' '+ str(vkt4[ikt]) +' '+ str(vnorm4[ikt]) +' '+ str(vkt2[2*ikt]) +' '+ str(vnorm2[2*ikt]) +' '+ str(vkt2[2*ikt+1]) +' '+ str(vnorm2[2*ikt+1]) +'\n' )

	ktoutfile.close()

	os.system( 'sed -i s/None/Nan/g '+locclu+'/zofr10kt_hien_2T.txt' )

	os.system( 'sed -i s/None/Nan/g '+locclu+'/zofr10_hien_2T.txt' )
	os.system( 'sed -i s/None/Nan/g '+locclu+'/zofr10_hien_2T_lab.txt' )
	os.system( 'sed -i s/None/Nan/g '+locclu+'/zofr10_hien_1T.txt' )
	os.system( 'sed -i s/None/Nan/g '+locclu+'/zofr10_hien_1T_lab.txt' )

	delete_window('all')

