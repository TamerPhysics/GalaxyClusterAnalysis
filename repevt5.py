
# PURPOSE: Reprocess all observations with the installed CALDB


# v2: * fixed a typo: the option check_vf_pha=yes is only now applied to acis_process_events
#     * new filename repro3_evtOBSID.fits
#     * replaced file.write(comands) by os.system(comands)
#     * no more gunzipping, as CIAO can handle gzipped files
#     * status filtering is different for data in FAINT and VFAINT modes

# v3: * use asol3.py to find asol, flt1, pbk, msk, bpix or evt1 files

# v4: * 2017-03-31: new version of asol: asol5_loc

import os
import glob
import commands

import asol6_loc

import pdb

def repevt(clu,obsid):

	print
	print 'Reprocessing evt file for OBSID '+obsid+' ------------------------'
	print

	if not os.access('mfe_'+clu+'/repro3_evt'+obsid+'.fits', os.F_OK) :

		acao = asol6_loc.asol(obsid, 'asol')

		bpix = asol6_loc.asol(obsid,'bpix')
		evt1 = asol6_loc.asol(obsid,'evt1')
		flt1 = asol6_loc.asol(obsid,'flt1')
		mtl  = asol6_loc.asol(obsid,'mtl')

		readmode = commands.getoutput('dmkeypar '+evt1+' readmode echo+')
		datamode = commands.getoutput('dmkeypar '+evt1+' datamode echo+')

		if readmode == 'TIMED' :
			if datamode == 'VFAINT' or datamode == 'FAINT' : edef = 'stdlev1'
			if datamode == 'GRADED' : edef = 'grdlev1'

		if readmode == 'CONTINUOUS' :
			if 'FAINT' in datamode : edef = 'cclev1'
			if 'GRADED' in datamode : edef = 'ccgrdlev1'

		# ACIS_PROCESS_EVENTS
		os.system('punlearn acis_process_events')
		os.system('pset acis_process_events infile='+evt1)
		os.system('pset acis_process_events outfile=mfe_'+clu+'/repro3_evt1_'+obsid+'.fits')
		os.system('pset acis_process_events badpixfile='+bpix)
		os.system('pset acis_process_events acaofffile='+acao)
		os.system('pset acis_process_events eventdef=")'+edef+'"')
		#os.system('pset acis_process_events rand_pix_size=0.5') # deprecated keyword
		os.system('pset acis_process_events rand_pha=yes')
		os.system('pset acis_process_events apply_tgain=yes')
		os.system('pset acis_process_events apply_cti=yes')
		os.system('pset acis_process_events mtlfile='+mtl)
		os.system('pset acis_process_events mode=h')
		if datamode == 'VFAINT' : os.system('pset acis_process_events check_vf_pha=yes')
		os.system('acis_process_events clobber=yes') # clobber=yes here is permanent

		# Grade and Status filtering
		if   datamode == 'VFAINT' : sstat='0'
		elif datamode == 'FAINT'  : sstat='00000000x00000000000000000000000' 
		#                           this status string chosen so that we stay consistent with blank sky file filtering
		#                           even if for FAINT mode, this status filtering is the same as status=0
		os.system('punlearn dmcopy')
		os.system('dmcopy "mfe_'+clu+'/repro3_evt1_'+obsid+'.fits[EVENTS][grade=0,2,3,4,6,status='+sstat+']" mfe_'+clu+'/repro3_flt_evt1_'+obsid+'.fits clobber=yes') # clobber=yes here is permanent

		# More filtering
		os.system('punlearn dmcopy')
		os.system('dmcopy "mfe_'+clu+'/repro3_flt_evt1_'+obsid+'.fits[EVENTS][@'+flt1+'][cols -pha]" mfe_'+clu+'/repro3_evt'+obsid+'.fits')

	else : print 'The repro3cessed file already exists'

	print '------ End of repevt for cluster ' + clu + ', OBSID ' + obsid + ' ----------------'
	print
	print	

