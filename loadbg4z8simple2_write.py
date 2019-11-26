
# PURPOSE: fit the background dataset. Called by checkbgfit.

# v3: * based on loadbgkt9
#     * no need to create BG files, they were created when loadbgkt9
#       was called from ktofr9.. just load them.

# v3simple: * use simple region reg/obsimple*.reg

# v6simple: * add analysis for ccd's that don't have infield BG

# v7simple: * use iobiccd1 instead of firstiob to link kT and Z for em component

# v8simple: * remove bad iobiccd=[-1,-1] line
#           * thaw the cluster emission normalization

# v8simple2_write: remove reference to acchif folder

# INPUT: * bg/bg'+ob+'_ccd'+ccd+'r.fits
# INPUT: * reg/obOBccdN.reg : region of CCD N in OBSIS OB

# OUTPUT: * no new files made

import commands
import os
import pdb

from sherpa.astro.ui import *
from pychips.all import *

import numpy as n

import bgmodel_flex11
import groupzeros

def bg(clu, obsids, presccdglob, presccd, srcstr, zz, nnhhlab, checkfit) :

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

			# Make BG spectrum from Blank-sky data using entire CCD (except for ACIS-I where we exclude
			# the region of high-background on chip #1)
			os.system('punlearn dmextract')
			os.system('dmextract "'+locclu+'/bg/bg'+ob+'_ccd'+ccd+'r.fits[sky=region('+locclu+'/reg/obsimple'+ob+'ccd'+ccd+'.reg)][bin PI]" '+locclu+'/spec/bgspsimple'+ob+'_ccd'+ccd+'.fits opt=pha1 wmap="[energy=300:2000][bin det=8]" mode=h clobber=yes')

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
			get_data('dummy'+ob+'.'+ccd).backscal = float(get_data('bg'+ob+'.'+ccd).backscal)

			# ========= BG model ===============
			exec 'delete_model("bg'+ob+'.'+ccd+'")'
			rsp0  = get_response( 'dummy'+ob+'.'+ccd )
			bgrsp = get_response( 'bg'+ob+'.'+ccd )
			notice_id('bg'+ob+'.'+ccd)
			group_bins('bg'+ob+'.'+ccd, 120)
			ignore_id('bg'+ob+'.'+ccd, ':.3,7:')
			exec 'set_full_model( "bg'+ob+'.'+ccd+'", bgrsp( powlaw1d.p1_'+ob+'_'+ccd+' + powlaw1d.p2_'+ob+'_'+ccd+' + exp.e1_'+ob+'_'+ccd+' + gauss1d.g1_'+ob+'_'+ccd+' +  gauss1d.g2_'+ob+'_'+ccd+' + gauss1d.g3_'+ob+'_'+ccd+' + gauss1d.g5_'+ob+'_'+ccd+' ) + rsp0( xsphabs.bgab * (powlaw1d.cxbbg_'+ob+'_'+ccd+') + xsapec.gal1bg_'+ob+'_'+ccd+' )  )'

			# ========= Source model ===============
			if ccd in presccd[iob] : 

				bgscl = get_data(ob+'.'+ccd).backscal * get_data(ob+'.'+ccd).exposure / get_data('bg'+ob+'.'+ccd).backscal / get_data('bg'+ob+'.'+ccd).exposure
				newhienbgscl = calc_data_sum(id=ob+'.'+ccd, lo=9.5, hi=12) / calc_data_sum(id='bg'+ob+'.'+ccd, lo=9.5, hi=12)
				relbs3 = get_data(ob+'.'+ccd).backscal / get_data('bg'+ob+'.'+ccd).backscal

				print ccd, ob, 'bgscl=', bgscl, '    newhienbgscl=', newhienbgscl, '      relbs3=', relbs3

				exec 'delete_model("'+ob+'.'+ccd+'")'
				rsp   = get_response(  ob+'.'+ccd )
				bgrsp = get_response( 'bg'+ob+'.'+ccd )
				#groupzeros.grp( ob+'.'+ccd, ':.3,7:')
				notice_id(ob+'.'+ccd)
				group_bins(ob+'.'+ccd, 120)
				ignore_id(ob+'.'+ccd, ':.3,7:')
				if ob!='7693' :
					exec 'set_full_model( "'+ob+'.'+ccd+'", rsp(xsphabs.ab * xsapec.em'+ob+'_'+ccd+' ) + newhienbgscl * (bgrsp( powlaw1d.p1_'+ob+'_'+ccd+' + powlaw1d.p2_'+ob+'_'+ccd+' + exp.e1_'+ob+'_'+ccd+' + gauss1d.g1_'+ob+'_'+ccd+' + gauss1d.g2_'+ob+'_'+ccd+' + gauss1d.g3_'+ob+'_'+ccd+' + gauss1d.g5_'+ob+'_'+ccd+' )) + rsp( xsphabs.ab * powlaw1d.cxb_'+ob+'_'+ccd+' + xsapec.gal1_'+ob+'_'+ccd+' ) )'
				else :
					exec 'set_full_model( "'+ob+'.'+ccd+'", rsp(xsphabs.ab * xsapec.em'+ob+'_'+ccd+' ) + bgscl * (bgrsp( powlaw1d.p1_'+ob+'_'+ccd+' + powlaw1d.p2_'+ob+'_'+ccd+' + exp.e1_'+ob+'_'+ccd+' + gauss1d.g1_'+ob+'_'+ccd+' + gauss1d.g2_'+ob+'_'+ccd+' + gauss1d.g3_'+ob+'_'+ccd+' + gauss1d.g5_'+ob+'_'+ccd+' )) + rsp( xsphabs.ab * powlaw1d.cxb_'+ob+'_'+ccd+' + xsapec.gal1_'+ob+'_'+ccd+' ) )'

			# Set for each OBSID, the same for all shells
			if ccd in presccd[iob] : bgmodel_flex11.setbgparams(ob, ccd, True)
			else :                   bgmodel_flex11.setbgparams(ob, ccd, False)

			# fix the BG parameters: make shape of BG fixed to best-fit, by linking norms 
			# of all components wrt to that of g1
			exec 'link(p1_'+ob+'_'+ccd+'.ampl, float(p1_'+ob+'_'+ccd+'.ampl.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
			exec 'link(p2_'+ob+'_'+ccd+'.ampl, float(p2_'+ob+'_'+ccd+'.ampl.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
			exec 'link(e1_'+ob+'_'+ccd+'.ampl, float(e1_'+ob+'_'+ccd+'.ampl.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
			exec 'link(g2_'+ob+'_'+ccd+'.ampl, float(g2_'+ob+'_'+ccd+'.ampl.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
			exec 'link(g3_'+ob+'_'+ccd+'.ampl, float(g3_'+ob+'_'+ccd+'.ampl.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
			exec 'link(g5_'+ob+'_'+ccd+'.ampl, float(g5_'+ob+'_'+ccd+'.ampl.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
			exec 'freeze(bgab.nh)'
			exec 'freeze(p1_'+ob+'_'+ccd+'.gamma)'
			exec 'freeze(p2_'+ob+'_'+ccd+'.gamma)'
			exec 'freeze(e1_'+ob+'_'+ccd+'.coeff)'

	iobiccd1=[-1,-1]
	for iob in range(len(obsids)) :
		ob = obsids[iob]
		for iccd in range(len(presccd[iob])) :
			ccd=presccd[iob][iccd]

			# Setting abundance parameter, kT and redshift
			exec 'set_par( em'+ob+'_'+ccd+'.redshift, val='+str(zz)+')'
			exec 'set_par( em'+ob+'_'+ccd+'.abundanc, val='+str(about)+', min=0, max=5, frozen=True)' 
			exec 'set_par( em'+ob+'_'+ccd+'.kt,       val='+str(ktout)+', min=0, max=40, frozen=True)'

			if iobiccd1==[-1,-1] :
				iobiccd1=[iob,iccd]
				# Set nH
				if clu != 'a478' : set_par(ab.nh, val=nnhhlab, min=nnhhlab/10., max=nnhhlab*10., frozen=True) 
				else :             set_par(ab.nh, val=nnhhlab*1.1, min=nnhhlab, max=nnhhlab*100., frozen=False) 

			else : 
				exec 'link( em'+ob+'_'+ccd+'.kt,       em'+obsids[iobiccd1[0]]+'_'+presccd[iobiccd1[0]][iobiccd1[1]]+'.kt )'
				exec 'link( em'+ob+'_'+ccd+'.abundanc, em'+obsids[iobiccd1[0]]+'_'+presccd[iobiccd1[0]][iobiccd1[1]]+'.abundanc )'


	# Fit each obsids by itself
	for iob in range(len(obsids)) :
		ob = obsids[iob]
		for iccd in range(len(presccdglob[iob])) :
			ccd=presccdglob[iob][iccd]

			print
			print '>>>>>>>>>>>>>>>> 1st fit '+ob+'.'+ccd+' <<<<<<<<<<<<<<<<<<<<<<<<'
			print

			set_method('moncar')
			exec 'fit( "bg'+ob+'.'+ccd+'")'

			exec 'plot_fit_resid( "bg'+ob+'.'+ccd+'")'

			if get_fit_results().rstat > 2. :

				exec 'unlink(p1_'+ob+'_'+ccd+'.ampl)'
				exec 'unlink(p2_'+ob+'_'+ccd+'.ampl)'
				exec 'unlink(e1_'+ob+'_'+ccd+'.ampl)'
				exec 'unlink(g2_'+ob+'_'+ccd+'.ampl)'
				exec 'unlink(g3_'+ob+'_'+ccd+'.ampl)'
				exec 'unlink(g5_'+ob+'_'+ccd+'.ampl)'

				exec 'thaw(p1_'+ob+'_'+ccd+'.ampl)'
				exec 'thaw(p2_'+ob+'_'+ccd+'.ampl)'
				exec 'thaw(e1_'+ob+'_'+ccd+'.ampl)'
				exec 'thaw(g2_'+ob+'_'+ccd+'.ampl)'
				exec 'thaw(g3_'+ob+'_'+ccd+'.ampl)'
				exec 'thaw(g5_'+ob+'_'+ccd+'.ampl)'

				print
				print '>>>>>>>>>>>>>>>> 2nd fit '+ob+'.'+ccd+' <<<<<<<<<<<<<<<<<<<<<<<<'
				print

				set_method('moncar')
				exec 'fit( "bg'+ob+'.'+ccd+'")'

				exec 'plot_fit_resid( "bg'+ob+'.'+ccd+'")'

				exec 'link(p1_'+ob+'_'+ccd+'.ampl, float(p1_'+ob+'_'+ccd+'.ampl.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
				exec 'link(p2_'+ob+'_'+ccd+'.ampl, float(p2_'+ob+'_'+ccd+'.ampl.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
				exec 'link(e1_'+ob+'_'+ccd+'.ampl, float(e1_'+ob+'_'+ccd+'.ampl.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
				exec 'link(g2_'+ob+'_'+ccd+'.ampl, float(g2_'+ob+'_'+ccd+'.ampl.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
				exec 'link(g3_'+ob+'_'+ccd+'.ampl, float(g3_'+ob+'_'+ccd+'.ampl.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
				exec 'link(g5_'+ob+'_'+ccd+'.ampl, float(g5_'+ob+'_'+ccd+'.ampl.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
		

				exec 'link(gal1bg_'+ob+'_'+ccd+'.norm, float(gal1bg_'+ob+'_'+ccd+'.norm.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
				exec 'link(cxbbg_'+ob+'_'+ccd+'.ampl,  float(cxbbg_'+ob+'_'+ccd+'.ampl.val  / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'

			exec 'freeze(g1_'+ob+'_'+ccd+'.ampl)'

#implement this on next version!!!!!!!!!!!!!!!!!!!!!!!!!!!			exec 'bs_rstat_'+ob+'_'+ccd+' = float( get_fit_results().rstat )'


	# Setting normalizations, link norms if backscal covers 95-105 % of annulus area
	for iob in range(len(obsids)) :
		ob = obsids[iob]
		for iccd in range(len(presccd[iob])) :
			ccd=presccd[iob][iccd]

			# Set source emission normalization
			exec 'set_par(em'+ob+'_'+ccd+'.norm, val=0.02, min=2e-4, max=2.)'
			exec 'emnorm0=abs(float(em'+ob+'_'+ccd+'.norm.val)*calc_data_sum(id="'+ob+'.'+ccd+'",lo=0.3,hi=7.0)/calc_model_sum(id="'+ob+'.'+ccd+'",lo=0.3,hi=7.0)/4.)'
			exec 'set_par(em'+ob+'_'+ccd+'.norm, val=float(emnorm0), min=0., max=float(emnorm0*1e3) )'

			if iobiccd1==[iob,iccd] :

				pass

				# thaw gal1 norm parameter!!!!
				#exec 'set_par(gal1_'+ob+'_'+ccd+'.norm, val=float(gal1bg_'+ob+'_'+ccd+'.norm.val), frozen=False)'

				# thaw cxb norm parameter!!!!
				#exec 'set_par(cxb_'+ob+'_'+ccd+'.ampl, val=float(cxbbg_'+ob+'_'+ccd+'.ampl.val), frozen=False)'

				# thaw gal1bg norm parameter!!!!
				#exec 'unlink(gal1bg_'+ob+'_'+ccd+'.norm)'
				#exec 'thaw(gal1bg_'+ob+'_'+ccd+'.norm)'

				# thaw cxbbg norm parameter!!!!
				#exec 'unlink(cxbbg_'+ob+'_'+ccd+'.ampl)'
				#exec 'thaw(cxbbg_'+ob+'_'+ccd+'.ampl)'

				#iobiccd1=[iob,iccd]

			else : 

				#link to first norm value
				relbs2 = get_data(ob+'.'+ccd).backscal / get_data(obsids[iobiccd1[0]]+'.'+presccd[iobiccd1[0]][iobiccd1[1]]).backscal
				#relbsexp= get_data(ob+'.'+ccd).backscal * get_data(ob+'.'+ccd).exposure / get_data(obsids[iobiccd1[0]]+'.'+presccd[iobiccd1[0]][iobiccd1[1]]).backscal / get_data(obsids[iobiccd1[0]]+'.'+presccd[iobiccd1[0]][iobiccd1[1]]).exposure

#8				exec 'link(em'+ob+'_'+ccd+'.norm, em'+obsids[iobiccd1[0]]+'_'+presccd[iobiccd1[0]][iobiccd1[1]]+'.norm )'

				# Link gal1 normalizations
				exec 'link(gal1_'+ob+'_'+ccd+'.norm, gal1_'+obsids[iobiccd1[0]]+'_'+presccd[iobiccd1[0]][iobiccd1[1]]+'.norm * relbs2 )'

				# Link cxb normalizations
				exec 'link(cxb_'+ob+'_'+ccd+'.ampl, cxb_'+obsids[iobiccd1[0]]+'_'+presccd[iobiccd1[0]][iobiccd1[1]]+'.ampl * relbs2 )'

				# Link gal1bg normalizations
				#exec 'link(gal1bg_'+ob+'_'+ccd+'.norm, gal1bg_'+obsids[iobiccd1[0]]+'_'+presccd[iobiccd1[0]][iobiccd1[1]]+'.norm )'

				# Link cxbbg normalizations
				#exec 'link(cxbbg_'+ob+'_'+ccd+'.ampl, cxbbg_'+obsids[iobiccd1[0]]+'_'+presccd[iobiccd1[0]][iobiccd1[1]]+'.ampl )'

	anydataset=False
	for iob in range(len(obsids)) : 
		if len(presccd[iob])>0 : anydataset=True


	if anydataset :
		print
		print '>>>>>>>>>>>>>>>> fit INFIELD <<<<<<<<<<<<<<<<<<<<<<<<'
		print
		set_method('moncar')
		exec 'fit('+srcstr[0:-1]+')'
		if checkfit : pdb.set_trace()

	for iob in range(len(obsids)) : 
		ob = obsids[iob]
		for iccd in range(len(presccdglob[iob])) :
			ccd=presccdglob[iob][iccd]

			exec 'link(p1_'+ob+'_'+ccd+'.ampl, float(p1_'+ob+'_'+ccd+'.ampl.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
			exec 'link(p2_'+ob+'_'+ccd+'.ampl, float(p2_'+ob+'_'+ccd+'.ampl.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
			exec 'link(e1_'+ob+'_'+ccd+'.ampl, float(e1_'+ob+'_'+ccd+'.ampl.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
			exec 'link(g2_'+ob+'_'+ccd+'.ampl, float(g2_'+ob+'_'+ccd+'.ampl.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
			exec 'link(g3_'+ob+'_'+ccd+'.ampl, float(g3_'+ob+'_'+ccd+'.ampl.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
			exec 'link(g5_'+ob+'_'+ccd+'.ampl, float(g5_'+ob+'_'+ccd+'.ampl.val / g1_'+ob+'_'+ccd+'.ampl.val) * g1_'+ob+'_'+ccd+'.ampl)'
			

			#exec 'freeze(g1_'+ob+'_'+ccd+'.ampl)'
			#exec 'freeze(gal1_'+ob+'_'+ccd+'.norm)'

			bgparfile = open(locclu+'/ktofrsp/bgparams8_'+ob+'_ccd'+ccd+'.txt','w')

			#if ccd in presccd[iob] :
			bgparfile.write(str(get_fit_results().rstat)+'\n')
			#else :
			#	exec 'bgparfile.write(str( bs_rstat_'+ob+'_'+ccd+' )+"\\n")'

			if ccd in presccd[iob] :
				exec 'bgparfile.write(str(gal1_'+ob+'_'+ccd+'.norm.val / get_data("'+ob+'.'+ccd+'").backscal )+"\\n")'
				exec 'bgparfile.write(str(cxb_'+ob+'_'+ccd+'.ampl.val / get_data("'+ob+'.'+ccd+'").backscal )+"\\n")'
			else :
				bgparfile.write('--no infield BG--\n')
				bgparfile.write('--no infield BG--\n')

			exec 'bgparfile.write(str(p1_'+ob+'_'+ccd+'.ampl.val)+"\\n")'
			exec 'bgparfile.write(str(p1_'+ob+'_'+ccd+'.gamma.val)+"\\n")'
			exec 'bgparfile.write(str(p2_'+ob+'_'+ccd+'.ampl.val)+"\\n")'
			exec 'bgparfile.write(str(p2_'+ob+'_'+ccd+'.gamma.val)+"\\n")'
			exec 'bgparfile.write(str(e1_'+ob+'_'+ccd+'.ampl.val)+"\\n")'
			exec 'bgparfile.write(str(e1_'+ob+'_'+ccd+'.coeff.val)+"\\n")'
			exec 'bgparfile.write(str(g1_'+ob+'_'+ccd+'.ampl.val)+"\\n")'
			exec 'bgparfile.write(str(g2_'+ob+'_'+ccd+'.ampl.val)+"\\n")'
			exec 'bgparfile.write(str(g3_'+ob+'_'+ccd+'.ampl.val)+"\\n")'
			exec 'bgparfile.write(str(g5_'+ob+'_'+ccd+'.ampl.val)+"\\n")'
			bgparfile.write('--placeholder for g6--\n')
			exec 'bgparfile.write(str(bgab.nh.val)+"\\n")'
			exec 'bgparfile.write(str(gal1bg_'+ob+'_'+ccd+'.norm.val / get_data("bg'+ob+'.'+ccd+'").backscal )+"\\n")'
			exec 'bgparfile.write(str(cxbbg_'+ob+'_'+ccd+'.ampl.val / get_data("bg'+ob+'.'+ccd+'").backscal )+"\\n")'
			bgparfile.close()

	set_method('levmar')






