
# PURPOSE: compute temperature profile for a given cluster using all observations

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

# ktofr13loc: * uses loadbg4z4simple, which uses a re-fitting of the BG spectra. the previous
#               BG fits often had offsets from the data

# ktofr14loc: * gal1 and cxb components are obtained from in-field spec fit, when present
#             * gal1 and cxb components' normalization is scaled, instead of scaling the entire
#               component when the model is set.
#             * gal1 and cxb components are scaled by relative backscal and NOT by backscal * exposure
#               like was done before, and like is done with instrumental BG components. this finding
#               was found using the testing script calibrate_spec_scaling3.py
#             * for certain OBSIDs, scale by bgscl and not by high energy scaling to model the BG. 
#               This is because these observations are missing counts at high energy. Their spectra
#               drop to zero, above a certain energy or PHA.

# OUTPUT: ktofrsp/spfix*
# OUTPUT: ktofrsp/warf14_*
# OUTPUT: ktofrsp/kt14an*
# OUTPUT: ktofr14.txt
# OUTPUT: ktofr14_hien.txt
# OUTPUT: ktofr14_norm3.txt
# OUTPUT: fitplots/ktofr14_CLU_OB_ISHELL

from sherpa.astro.ui import *
from pychips.all import *

import os
import commands
import pdb
import gc

import numpy as n

import centradec
import asol6_loc
import radreg4simple3
import rmfarf4bg
import loadbg4z8simple3
import groupzeros
import globmod
import whichccd2z_loc


def ktofr(clu, obsids, r1mpc, zz, nnhhlab) :

	print ' -----------> ktofr14loc for '+ clu

	clean()
	gc.collect()
	set_stat('cstat')

	# local directory for this cluster
	locclu = 'mfe_'+clu

	# coords of center of the cluster
	(rac, decc) = centradec.getrd(locclu)

	# don't use these observations
	badhienobs = ['7686','7688','7689','7690','7692','7693','7694','7696','7701']

	if not os.path.exists(locclu+'/ktofrsp/') : os.mkdir(locclu+'/ktofrsp/')

	# Chandra CCDs to look for data in
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

	# create simplechipsreg regions.
	for ob in obsids :
		os.system('punlearn skyfov')
		os.system('skyfov '+locclu+'/clean'+ob+'.fits '+locclu+'/reg/simplechipsreg'+ob+'.fits clobber=yes')

	# which CCD's are present in each OBSID
	presccdglob=[]
	for iob in range(len(obsids)) : presccdglob.append( whichccd2z_loc.find(clu, obsids[iob]) )

	# Reset these global variables
	globmod.globstr1=''
	globmod.globstr2=''

	for iob in range(len(obsids)) : 

		ob = obsids[iob]

		# Make file with only counts in energy range .3-7keV
		os.system('punlearn dmcopy')
		os.system('dmcopy "'+locclu+'/clean'+ob+'.fits[energy=300:7000]" '+locclu+'/evt'+ob+'_b.fits')

		# get the coordinates of the center in SKY units
		os.system('punlearn dmcoords')
		os.system('dmcoords '+locclu+'/clean'+ob+'.fits asolfile='+asol6_loc.asol(ob, 'asol')+' ra='+rac+' dec='+decc+' option=cel celfmt=hms')
		xs = commands.getoutput('pget dmcoords x')
		ys = commands.getoutput('pget dmcoords y')

		# Make region containing all ptsrcs and hi-BG corner region
		os.system('cp '+locclu+'/reg/pt0mfe_wcs.reg '+locclu+'/reg/pt_hibg'+ob+'.reg')
		if os.path.exists(locclu+'/reg/hibgcorner'+ob+'.reg') : os.system('more '+locclu+'/reg/hibgcorner'+ob+'.reg >> '+locclu+'/reg/pt_hibg'+ob+'.reg')


		# Make file with the radial coordinate of each photon count
		os.system('dmtcalc "'+locclu+'/evt'+ob+'_b.fits[col -time,-ccd_id,-node_id,-expno,-chip,-tdet,-det,-phas,-pha_ro,-energy,-pi,-fltgrade,-grade,-status][exclude sky=region('+locclu+'/reg/pt_hibg'+ob+'.reg)]" "'+locclu+'/evtrad'+ob+'.txt[opt kernel=text/simple]" expr="r2=(('+xs+'-sky[0])^2)+(('+ys+'-sky[1])^2)"')

	# Get the radial coord of EACH photon count!
	# To be used to determine if we have enough counts
	# for a given radial bin to get 10% error on best-fit
	# kT value.
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

	# Define rmin and rmax
	rmax = r[-1]
	rmin = 1.2 / r1mpc 

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
			ktin  = float( intermed.split()[1] )
			ktout = float( intermed.split()[3] )
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

	# Load Hi-energy count in losrc regions
	losrchienct = n.ones( (len(obsids),4) ) * -1
	for iob in range(len(obsids)) :
		ob = obsids[iob]
		sclfile = open( locclu+'/ctofr/scl'+ob+'.txt' , 'r')
		for iccd in range(4) :
			intermed = sclfile.readline()
			if intermed.split()[2].lower() != 'nan' : losrchienct[iob, iccd] = float(intermed.split()[2])
		sclfile.close()

	# Load string containing exculsion of all point sources
	excstr=[]
	for ob in obsids :
		excfile = open(locclu+'/reg/exclpt'+ob+'.reg', 'r')
		excstr.append( excfile.readline() )
		excfile.close()



	################### make radial bins which contain enough counts for a best-fit kT with 10% error ########################
	ccdindex = {'i':0, '5':1, '6':2, '7':3}
	vrstat=[]
	vkt=[]
	vktp=[]
	vktm=[]
	ishell = 0
	rz1=[]
	rz2=[]
	ishell=0 # the index of everything created new, in this loop. only increases, when we get a satisfactory fit.
	i1=0 # index to rbnd indicating inner radius
	i2=1 # index to rbnd indicating outer radius
	while i2 <= len(rbnd)-1 :

		ct = 1
		while ct <= 10 :

			print '======== loop ct =', ct

			# get annulus temperature from cluname_TZ.txt, based on in/out regions
			if ishell==0 :
				if   (rbnd[i1]+rbnd[i2])/2 <= 0.1 and ktin > 0 : tshell = ktin
				elif (rbnd[i1]+rbnd[i2])/2 >  0.1 and ktout> 0 : tshell = ktout
				if   ktin  < 0 and ktout > 0 : tshell = ktout
				elif ktout < 0 and ktin  > 0 : tshell = ktin
				else : tshell=5.
			else : tshell = vkt[ishell-1]

			# Calculate BG fraction in the annulus
			sclbgct = 0
			for iob in range(len(obsids)) :
				ob = obsids[iob]
				annfile = open('antemp.reg','w')
				annfile.write('annulus('+rac+','+decc+','+str(rbnd[i1]*r1mpc)+'",'+str(rbnd[i2]*r1mpc)+'")\n')
				annfile.write(excstr[iob])
				annfile.close()
				for ccd in presccdglob[iob] :
					sclbgct = sclbgct + float( commands.getoutput('dmlist "'+locclu+'/bg/bg'+ob+'_ccd'+ccd+'r_en.fits[sky=region(antemp.reg)]" counts') ) / float( commands.getoutput('dmkeypar '+locclu+'/bg/bg'+ob+'_ccd'+ccd+'r_en.fits exposure echo+') ) * float( commands.getoutput('dmkeypar '+locclu+'/evt'+ob+'_b.fits exposure echo+') )
				os.remove('antemp.reg')
			bgfrac = sclbgct / float(i2-i1) / ct1bin

			# given approx. kT and BGfrac, how many counts do we need for 10% error?
			nnec=( 500 * 10.**(1.976*bgfrac) ) * (tshell/2.)**1.7
			nnec= max(nnec, 153.) # 153counts is the value at BGf=0, kT=1keV. 
                                     # We take this as the min number of counts for making a kT bin

			inec0 = min(   int( nnec / ct1bin) +i1 +1  ,  len(rbnd)-1  )

			# making sure that rhi/rlo >= 1.25:
			i125 = n.where( rbnd >= 1.25*rbnd[i1] )
			if len(i125[0]) == 0 : i125=[[len(rbnd)-1]]
			inec = max(inec0, i125[0][0])

			if inec <= i2 or inec == len(rbnd)-1 : 
				i2=inec
				break
			else : i2=inec

			ct=ct+1

		if ct>=11 : os.system('echo '+str(ishell)+' >> '+locclu+'/shells_nocnvg.txt')

		rz1.append( rbnd[i1] )
		rz2.append( rbnd[i2] )

		################# make spec and fit ####################################

		# Get the CCD's present in this shell
		presccd=[] # the ccd's that are present in this SHELL!
		for iob in range(len(obsids)) :
			ob = obsids[iob]

			# radreg4 returns the ccd's that are present in this SHELL!
			# and makes the regions necessary for them to be used in fits..
			intermed1 = radreg4simple3.radreg(clu, ob, rz1[-1]*r1mpc/.492, rz2[-1]*r1mpc/.492, locclu+'/ktofrsp/kt14an'+ob+'_'+str(ishell))
			presccd.append( intermed1 )

			for ccd in str(presccd[iob]) :

				# Make the spectra for SHELL x CCD
				os.system('punlearn dmextract')
				os.system('dmextract "'+locclu+'/clean'+ob+'.fits[ccd_id='+ccdlist[ccdindex[ccd]]+'][sky=region('+locclu+'/ktofrsp/kt14an'+ob+'_'+str(ishell)+'_ccd'+ccd+'_xfov_pt_simple.reg)][bin PI]" '+locclu+'/ktofrsp/sp14_'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits opt=pha1 wmap="[energy=300:2000][bin det=8]" clobber=yes' )

				# Make ARF given the WMAP from the above spectrum file
				#pbk = asol6_loc.asol(ob,'pbk')
				acao = asol6_loc.asol(ob, 'asol')
				os.system('punlearn mkwarf')
				os.system('pset mkwarf infile="'+locclu+'/ktofrsp/sp14_'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits[WMAP]"')
				os.system('pset mkwarf outfile='+locclu+'/ktofrsp/warf14_'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits')
				os.system('pset mkwarf weightfile='+locclu+'/ktofrsp/wfef14_'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits')
				#os.system('pset mkwarf pbkfile='+pbk)
				os.system('pset mkwarf egridspec="0.3:11.0:0.01"')
				os.system('pset mkwarf asolfile='+acao)
				os.system('mkwarf verbose=1 mode=h clobber=yes')
				os.system('rm '+locclu+'/ktofrsp/wfef14_'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits')

		# Make BG spectra and load them, only once for ishell=0
		# This is here because it uses ktofrsp/ktan14*.reg, created just
		# above in radreg4simple2.radreg()
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

		# Make sure that all the needed files were created. If not, remove from prescccd
		# any CCD wich is missing PI, ARF or RMF file
		firstiob=-1
		for iob in range(len(obsids)) :
			ob = obsids[iob]
			for ccd in str(presccd[iob]) :
				if not os.path.exists(locclu+'/ktofrsp/warf14_'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits') \
           or not os.path.exists(locclu+'/spec/sp'+ob+'_ccd'+ccd+'_center.wrmf') \
           or not os.path.exists(locclu+'/ktofrsp/sp14_'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits') \
           or ccd not in presccdglob[iob] : presccd[iob] = presccd[iob].replace(ccd, '')
			if len(presccd[iob]) > 0 and firstiob==-1 : firstiob=iob


		print
		print '))))))))))))) clu = ', clu
		print '))))))))))))) obsids = ', obsids
		print '))))))))))))) ishell = ', ishell
		print '))))))))))))) presccdglob = ', presccdglob
		print '))))))))))))) presccd = ', presccd
		print

		# Load data we just made
		srcstr=''
		for iob in range(len(obsids)) :
			ob = obsids[iob]

			for ccd in presccd[iob] :

				load_pha( ob+'.'+str(ishell)+'.'+ccd , locclu+'/ktofrsp/sp14_'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits' )
				load_rmf( ob+'.'+str(ishell)+'.'+ccd , locclu+'/spec/sp'+ob+'_ccd'+ccd+'_center.wrmf'  )
				load_arf( ob+'.'+str(ishell)+'.'+ccd , locclu+'/ktofrsp/warf14_'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits' )

				bgscl = get_data(ob+'.'+str(ishell)+'.'+ccd).backscal * get_data(ob+'.'+str(ishell)+'.'+ccd).exposure / get_data('bg'+ob+'.'+ccd).backscal / get_data('bg'+ob+'.'+ccd).exposure
				newhienbgscl = calc_data_sum(id=ob+'.'+str(ishell)+'.'+ccd, lo=9.5, hi=12) / calc_data_sum(id='bg'+ob+'.'+ccd, lo=9.5, hi=12)
				srcstr = srcstr +'"'+ob+'.'+str(ishell)+'.'+ccd + '",'

				# the source model
				rsp   = get_response(  ob+'.'+str(ishell)+'.'+ccd )
				bgrsp = get_response( 'bg'+ob+'.'+ccd )
				groupzeros.grp( ob+'.'+str(ishell)+'.'+ccd, ':.3,7:')

				if ob not in badhienobs :
					exec 'set_full_model( "'+ob+'.'+str(ishell)+'.'+ccd+'", rsp(xsphabs.ab * xsapec.em'+str(ishell)+'_'+ob+'_'+ccd+' ) + newhienbgscl * bgrsp( powlaw1d.p1_'+ob+'_'+ccd+' + powlaw1d.p2_'+ob+'_'+ccd+' + exp.e1_'+ob+'_'+ccd+' + gauss1d.g1_'+ob+'_'+ccd+' + gauss1d.g2_'+ob+'_'+ccd+' + gauss1d.g3_'+ob+'_'+ccd+' + gauss1d.g5_'+ob+'_'+ccd+' ) + rsp( xsphabs.ab * powlaw1d.cxbbg_'+ob+'_'+ccd+' + xsapec.gal1bg_'+ob+'_'+ccd+' ) )'
				else :
					exec 'set_full_model( "'+ob+'.'+str(ishell)+'.'+ccd+'", rsp(xsphabs.ab * xsapec.em'+str(ishell)+'_'+ob+'_'+ccd+' ) + bgscl * bgrsp( powlaw1d.p1_'+ob+'_'+ccd+' + powlaw1d.p2_'+ob+'_'+ccd+' + exp.e1_'+ob+'_'+ccd+' + gauss1d.g1_'+ob+'_'+ccd+' + gauss1d.g2_'+ob+'_'+ccd+' + gauss1d.g3_'+ob+'_'+ccd+' + gauss1d.g5_'+ob+'_'+ccd+' ) + rsp( xsphabs.ab * powlaw1d.cxbbg_'+ob+'_'+ccd+' + xsapec.gal1bg_'+ob+'_'+ccd+' ) )'

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


			# Setting abundance parameter, kT and redshift
			for ccd in presccd[iob] :
				exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'.redshift, val='+str(zz)+')'
				exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'.abundanc, val=0.3, min=0, max=5, frozen=False)' 
				exec 'set_par( em'+str(ishell)+'_'+ob+'_'+ccd+'.kt,       val=5, min=0, max=40, frozen=False)'

				if iob != firstiob or ccd != presccd[firstiob][0] : 
					exec 'link( em'+str(ishell)+'_'+ob+'_'+ccd+'.kt,       em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'.kt )'
					exec 'link( em'+str(ishell)+'_'+ob+'_'+ccd+'.abundanc, em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'.abundanc )'

		if srcstr[0:-1] != '' : # if there is data to fit:

			# Set nH
			if clu != 'a478' : set_par(ab.nh, val=nnhhlab, min=nnhhlab/10., max=nnhhlab*10., frozen=True) 
			else :             set_par(ab.nh, val=nnhhlab*1.1, min=nnhhlab, max=nnhhlab*100., frozen=False) 

			# Setting normalizations, link norms if backscal covers 95-105 % of annulus area
			iobiccd1=[-1,-1]
			anarea = n.pi * ( (rz2[-1]*r1mpc/.492)**2. - (rz1[-1]*r1mpc/.492)**2. )
			for iob in range(len(obsids)) :
				ob = obsids[iob]
				for iccd in range(len(presccd[iob])) :
					ccd=presccd[iob][iccd]
					bspix2 = float( commands.getoutput('dmkeypar '+locclu+'/ktofrsp/sp14_'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits backscal echo+') ) *64.0*1024.**2.

					if iobiccd1==[-1,-1] or bspix2/anarea < 0.95 or bspix2/anarea > 1.05 :
						exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'.norm, val=0.02, min=2e-4, max=2.)'
						exec 'emnorm0=abs(float(em'+str(ishell)+'_'+ob+'_'+ccd+'.norm.val)*calc_data_sum(id="'+ob+'.'+str(ishell)+'.'+ccd+'",lo=0.3,hi=7.0)/calc_model_sum(id="'+ob+'.'+str(ishell)+'.'+ccd+'",lo=0.3,hi=7.0))'
						exec 'set_par(em'+str(ishell)+'_'+ob+'_'+ccd+'.norm, val=float(emnorm0), min=0., max=float(emnorm0*1e3) )'
						if iobiccd1==[-1,-1] and bspix2/anarea >= 0.95 and bspix2/anarea <= 1.05 : iobiccd1=[iob,iccd]
					else : #link to first norm value
						relbs = bspix2 / float( commands.getoutput('dmkeypar '+locclu+'/ktofrsp/sp14_'+obsids[iobiccd1[0]]+'_'+str(ishell)+'_ccd'+presccd[iobiccd1[0]][iobiccd1[1]]+'.fits backscal echo+') ) /(64.0*1024.**2.)
						exec 'link(em'+str(ishell)+'_'+ob+'_'+ccd+'.norm, em'+str(ishell)+'_'+obsids[iobiccd1[0]]+'_'+presccd[iobiccd1[0]][iobiccd1[1]]+'.norm * relbs )'


			# Fitting
			exec 'fit('+srcstr[0:-1]+')'

			# plotting
			if ishell!=0 : delete_window('all')
			x1=[0.,0.5,0.,0.5]
			y1=[0.5,0.5,0.,0.]
			for iob in range(len(obsids)) : 
				ob=obsids[iob]
				if len(presccd[iob]) > 0 : 
					add_window(8.5,8.5,'inches')
				for iccd in range(len(presccd[iob])) :
					ccd=presccd[iob][iccd]
					add_frame(x1[iccd], y1[iccd], x1[iccd]+0.5, y1[iccd]+0.5)
					add_curve( get_data_plot(ob+'.'+str(ishell)+'.'+ccd).x, get_data_plot(ob+'.'+str(ishell)+'.'+ccd).y, dataplot )
					add_histogram( get_model_plot(ob+'.'+str(ishell)+'.'+ccd).xlo, get_model_plot(ob+'.'+str(ishell)+'.'+ccd).xhi, get_model_plot(ob+'.'+str(ishell)+'.'+ccd).y, modelplot )
					set_plot_title(clu+' OB'+ob+' CCD'+ccd)
					split(2)
					add_curve( get_resid_plot(ob+'.'+str(ishell)+'.'+ccd).x, get_resid_plot(ob+'.'+str(ishell)+'.'+ccd).y, dataplot )
					add_hline(0)
				if os.path.exists('fitplots/ktofr14_'+clu+'_'+ob+'_'+str(ishell)+'_hien.eps') : os.remove('fitplots/ktofr14_'+clu+'_'+ob+'_'+str(ishell)+'_hien.eps')
				if len(presccd[iob]) > 0 : print_window( 'fitplots/ktofr14_'+clu+'_'+ob+'_'+str(ishell)+'_hien', ['format', 'eps', 'orientation', 'landscape'])

			normv=-1*n.ones(( len(obsids),4 ))
			if get_fit_results().rstat < 3 : 

				# Error calculation
				exec 'proj('+srcstr+' em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'.kt )'

				#15 if get_proj_results().parmaxes[0]!=None or ishell2==len(r1)-1 :

				vkt.append( get_proj_results().parvals[0] )
				vktp.append( get_proj_results().parmaxes[0] )
				vktm.append( get_proj_results().parmins[0] )
				vrstat.append( get_fit_results().rstat )

				# save normalizations
				for iob in range(len(obsids)) :
					ob = obsids[iob]
					for iccd in range(len(presccd[iob])) :
						ccd=presccd[iob][iccd]
						exec 'normv[iob,ccdindex[ccd]]=em'+str(ishell)+'_'+ob+'_'+ccd+'.norm.val'

				# best-fit Z
				exec 'bestfitz = em'+str(ishell)+'_'+obsids[firstiob]+'_'+presccd[firstiob][0]+'.abundanc.val'
				# best-fit nH
				bestfitnh = ab.nh.val

			else :

				# No error calculation for bad fits
				vkt.append( n.nan )
				vktp.append( n.nan )
				vktm.append( n.nan )
				vrstat.append( get_fit_results().rstat )

				# save normalizations
				for iob in range(len(obsids)) :
					ob = obsids[iob]
					for iccd in range(len(presccd[iob])) :
						ccd=presccd[iob][iccd]
						normv[iob,ccdindex[ccd]]=n.nan

				# best-fit Z
				bestfitz = n.nan
				# best-fit nH
				bestfitnh = n.nan

			globmod.globstr1 = globmod.globstr1 + str(rz1[-1]) +' '+ str(rz2[-1]) +' '+ str(vkt[-1]) +' '+ str(vktm[-1]) +' '+ str(vktp[-1]) +' '+ str(vrstat[-1]) +' '+ str(-1) +' '+ str(-1) +' '+ str(-1) +' '+ str(bestfitz) +' '+ str(bestfitnh) +'\n'



			# save normalization results:

			# rin, rout
			globmod.globstr3 = globmod.globstr3 + str(rz1[-1]) +'\t'+ str(rz2[-1]) +'\t'

			# the normalization for each (ob,ccd)
			for iob in range(len(obsids)) :
				ob = obsids[iob]
				for ccd in presccd[iob] : globmod.globstr3 = globmod.globstr3 +'-1' +' ' #13 used to be: +str(normv[iob,ccdindex[ccd]]) +' '
			globmod.globstr3 = globmod.globstr3 + '\t'

			# the solid angle covered by dataset as a fraction of total
			# annulus solid angle, PI (Rout^2-Rin^2)
			for iob in range(len(obsids)) :
				ob = obsids[iob]
				for ccd in presccd[iob] :
					bspix2 = float( commands.getoutput('dmkeypar '+locclu+'/ktofrsp/sp14_'+ob+'_'+str(ishell)+'_ccd'+ccd+'.fits backscal echo+') ) *64.0*1024.**2.
					globmod.globstr3 = globmod.globstr3 + str(bspix2/anarea) + ' '

			# Reduced statistic
			globmod.globstr3 = globmod.globstr3 +'\t'+ '-1' +'\t' #13 used too be: str(vrstat[-1]) + '\t'

			# Number of net counts, using backscal scaling
			for iob in range(len(obsids)) :
				ob = obsids[iob]
				for ccd in presccd[iob] :
					bgct = get_data(ob+'.'+str(ishell)+'.'+ccd).backscal * get_data(ob+'.'+str(ishell)+'.'+ccd).exposure / get_data('bg'+ob+'.'+ccd).backscal / get_data('bg'+ob+'.'+ccd).exposure * calc_data_sum(id='bg'+ob+'.'+ccd, lo=0.3, hi=7.0) # the number of BG counts to be subtracted
					globmod.globstr3 = globmod.globstr3 + str( max(calc_data_sum( id=ob+'.'+str(ishell)+'.'+ccd, lo=0.3, hi=7.0 ) - bgct,1.) ) + " "
			globmod.globstr3 = globmod.globstr3 +'\t'

			# the --HE-- normalization for each (ob,ccd)
			for iob in range(len(obsids)) :
				ob = obsids[iob]
				for ccd in presccd[iob] : globmod.globstr3 = globmod.globstr3 + str(normv[iob,ccdindex[ccd]]) +' ' #13 normv used to be normhev
			globmod.globstr3 = globmod.globstr3 + '\t'

			# Reduced statistic -- HE --
			globmod.globstr3 = globmod.globstr3 +'\t'+ str(vrstat[-1])  #13 vrstat used to be vrstathien

			globmod.globstr3 = globmod.globstr3 +'\n'



			############### end of fitting if statement ############################


		# For the next ishell loop:
		i1=i2
		i2=i1+1
		ishell=ishell+1

		########################### end of ishell loop ##############################


	kthfile=open(locclu+'/ktofr14_hien.txt', 'w')
	kthfile.write(globmod.globstr1)
	kthfile.close()

	normfile=open(locclu+'/ktofr14_norm3.txt', 'w')
	normfile.write(globmod.globstr3)
	normfile.close()

	os.system( 'sed -i s/None/Nan/g '+locclu+'/ktofr14_hien.txt' )
	os.system( 'sed -i s/None/Nan/g '+locclu+'/ktofr14_norm3.txt' )



