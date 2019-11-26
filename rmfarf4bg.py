
# v2: * use asol3.py to find asol, flt1, pbk, msk, bpix or evt1 files
# v3: * based on rmf4entsoft2.py
#     * for both ACIS-I AND ACIS-S

# rmfarf4bg: * based on rmf4ent3, with major modifications
#            * added on 1/4/12: check that the identical RMF is not 
#              itself the one we want to make 

# PURPOSE:
# make ARF and RMF to be used with the blank-sky spectra


# INPUT: evtfile (should be cleanOB.fits)
# INPUT: spfile  (should be spec/bgspOB_ccdN.fits)
# INPUT: arffile: desired name for output ARF
# INPUT: rmffile: desired name for output RMF

# OUTPUT: ARF and RMF files

import commands
import os
import pdb

import globrmf

def rmf(evtfile, spfile, arffile, rmffile, ccdstr) : 

	caldbfolder = commands.getoutput('echo $CALDB')

	# ARF
	if os.path.exists('temp_wfef.fits') : os.remove('temp_wfef.fits')
	os.system('punlearn mkwarf')
	os.system('pset mkwarf infile='+spfile+'[WMAP]')
	os.system('pset mkwarf outfile='+arffile)
	os.system('pset mkwarf weightfile=temp_wfef.fits')
	#os.system('pset mkwarf pbkfile=none')
	os.system('pset mkwarf egridspec="0.3:11.0:0.01"')
	os.system('pset mkwarf dafile=none')
	os.system('mkwarf verbose=1 mode=h')

	# RMF -----------------------
	fptemp = float(commands.getoutput('dmkeypar '+evtfile+' fp_temp echo+'))
	tgaincor = commands.getoutput('dmkeypar '+evtfile+' tgaincor echo+')
	cticor  = commands.getoutput('dmkeypar '+evtfile+' cti_corr echo+')
	ctiapp  = commands.getoutput('dmkeypar '+evtfile+' cti_app echo+')


	# using mkacisrmf or...
	if (fptemp < 158. or (fptemp >= 158. and ccdstr in ['5','7']) ) and ( ctiapp != 'NNNNNNNNNN' or cticor == 'TRUE' or cticor == '1' ) and (tgaincor == 'T') :

		gainfile  = commands.getoutput('dmkeypar '+evtfile+' gainfile echo+')
		gainver   = gainfile.split('.fits')[0][-1]
		caldbls   = commands.getoutput('ls $CALDB/data/chandra/acis/p2_resp/'+gainfile[0:15]+'*'+gainver+'.fits')
		caldbresp = caldbls.split('\n')

		if len(caldbresp) != 1 : 
			print
			print
			print '  ############### MORE OR LESS THAN ONE P2RESP FILES FOUND !!!!!!!!!! ###################'
			print
			print

#		if gainfile in globrmf.gain and caldbresp[0] in globrmf.scat :
		rmflinked=False
		for iglob in range(len(globrmf.gain)) :
			if gainfile == globrmf.gain[iglob].split('/')[-1] and caldbresp[0] == globrmf.scat[iglob] and ccdstr == globrmf.name[iglob][-6] and rmffile != globrmf.name[iglob][globrmf.name[iglob].rfind('mfe_'):] : 
				if os.path.exists(rmffile) : os.remove(rmffile)
				os.symlink(globrmf.name[iglob], rmffile)
				rmflinked=True
				break


		if not rmflinked :
			os.system('punlearn mkacisrmf')
			os.system('pset mkacisrmf infile='+caldbresp[0])
			os.system('pset mkacisrmf outfile='+rmffile)
			os.system('pset mkacisrmf energy="0.3:11.0:0.01" channel="1:1024:1"')
			os.system('pset mkacisrmf chantype="PI"')
			os.system('pset mkacisrmf wmap='+spfile+'[WMAP]')
			os.system('pset mkacisrmf gain='+caldbfolder+'/data/chandra/acis/det_gain/'+gainfile)
			os.system('pset mkacisrmf asolfile=none')
			os.system('mkacisrmf mode=h verbose=1')

			globrmf.gain.append(caldbfolder+'/data/chandra/acis/det_gain/'+gainfile)
			globrmf.scat.append(caldbresp[0])
			globrmf.name.append('/caviar/home/research/acchif/'+rmffile)

	# using mkrmf
	else :

		os.system('punlearn mkrmf')
		os.system('pset mkrmf infile=CALDB')
		os.system('pset mkrmf outfile='+rmffile)
		os.system('pset mkrmf axis1="energy=0.3:11.0:0.01"')
		os.system('pset mkrmf axis2="pi=1:1024:1"')
		os.system('pset mkrmf weights=temp_wfef.fits')
		os.system('mkrmf verbose=1 mode=h')

	if os.path.exists('temp_wfef.fits') : os.remove('temp_wfef.fits')
