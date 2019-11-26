
# PURPOSE: merge evt files
# PURPOSE: compute the coordinates of the center of the cluster.

# OUTPUT: repro3_evt_nophasOBS.fits (same as orig_CLU/repro3_evtOBS.fits minus phas column)
# OUTPUT: mrgevt.fits
# OUTPUT: rac.txt and decc.txt

# v3: for ent.py
#     also, changed repro2 to repro3
#     in hicount/ deleted unused code and comments

# v4: * use asol3.py to find asol, flt1, pbk, msk, bpix or evt1 files

# v5: * 2017-03-31 simple bug fixed, where there was 2 consecutive ' characters
#     * 2017-03-31: new version of asol: asol5_loc

# v6: use the new location of repro3_evtOBSID.fits in mfe_OBSID/

import os
import commands
import numpy
import shutil

import asol6_loc
import hardsoft3

def mrg( clu, obsids, longest ) :

	origclu = 'orig_' + clu
	locclu = 'mfe_' + clu

	acao = asol6_loc.asol(longest, 'asol')

	if not os.path.exists(locclu+'/mrgevt.fits')  :

		print
		print 'Merging evt files for cluster '+clu+' ------------------------'
		print

		hardob = hardsoft3.hardsoft(clu, [longest], 'h')
		softob = hardsoft3.hardsoft(clu, [longest], 's')

		# determining which chips to use in computing the centroid
		if len(hardob)==1 : chipstr0 = '[ccd_id=0,1,2,3]' # ACIS-I: straightforward: chips 0,1,2,3
		else : 
			# choose the soft chip with most count to compute centroid
			chipcount=numpy.zeros(3)
			softchip=['5','6','7']
			for ichip in range(3)  : chipcount[ichip] = int(commands.getoutput('dmlist "mfe_'+clu+'/repro3_evt'+longest+'.fits[ccd_id='+softchip[ichip]+']" counts'))
			chipstr0='[ccd_id='+softchip[numpy.argmax(chipcount)]+']'

		# (x1, y1)
		os.system('punlearn dmstat')
		os.system('dmstat "mfe_'+clu+'/repro3_evt'+longest+'.fits'+chipstr0+'[energy=300:7000][cols x]" median=yes')
		x1=commands.getoutput('pget dmstat out_median')
		sx1=commands.getoutput('pget dmstat out_sigma')
		os.system('punlearn dmstat')
		os.system('dmstat "mfe_'+clu+'/repro3_evt'+longest+'.fits'+chipstr0+'[energy=300:7000][cols y]" median=yes')
		y1=commands.getoutput('pget dmstat out_median')
		sy1=commands.getoutput('pget dmstat out_sigma')

		# (x2, y2)
		os.system('punlearn dmstat')
		os.system('dmstat "mfe_'+clu+'/repro3_evt'+longest+'.fits'+chipstr0+'[energy=300:7000][sky=ellipse('+x1+','+y1+','+str(2.*float(sx1))+','+str(2.*float(sy1))+',0)][cols x]" median=yes')
		x2=commands.getoutput('pget dmstat out_median')
		sx2=commands.getoutput('pget dmstat out_sigma')
		os.system('punlearn dmstat')
		os.system('dmstat "mfe_'+clu+'/repro3_evt'+longest+'.fits'+chipstr0+'[energy=300:7000][sky=ellipse('+x1+','+y1+','+str(2.*float(sx1))+','+str(2.*float(sy1))+',0)][cols y]" median=yes')
		y2=commands.getoutput('pget dmstat out_median')
		sy2=commands.getoutput('pget dmstat out_sigma')

		# (x3, y3)
		os.system('punlearn dmstat')
		os.system('dmstat "mfe_'+clu+'/repro3_evt'+longest+'.fits'+chipstr0+'[energy=300:7000][sky=ellipse('+x2+','+y2+','+str(1.*float(sx2))+','+str(1.*float(sy2))+',0)][cols x]" median=yes')
		x3=commands.getoutput('pget dmstat out_median')
		sx3=commands.getoutput('pget dmstat out_sigma')
		os.system('punlearn dmstat')
		os.system('dmstat "mfe_'+clu+'/repro3_evt'+longest+'.fits'+chipstr0+'[energy=300:7000][sky=ellipse('+x2+','+y2+','+str(1.*float(sx2))+','+str(1.*float(sy2))+',0)][cols y]" median=yes')
		y3=commands.getoutput('pget dmstat out_median')
		sy3=commands.getoutput('pget dmstat out_sigma')

		# (x4, y4)
		os.system('punlearn dmstat')
		os.system('dmstat "mfe_'+clu+'/repro3_evt'+longest+'.fits'+chipstr0+'[energy=300:7000][sky=ellipse('+x3+','+y3+','+str(0.75*float(sx3))+','+str(0.75*float(sy3))+',0)][cols x]" median=yes')
		x4=commands.getoutput('pget dmstat out_median')
		sx4=commands.getoutput('pget dmstat out_sigma')
		os.system('punlearn dmstat')
		os.system('dmstat "mfe_'+clu+'/repro3_evt'+longest+'.fits'+chipstr0+'[energy=300:7000][sky=ellipse('+x3+','+y3+','+str(0.75*float(sx3))+','+str(0.75*float(sy3))+',0)][cols y]" median=yes')
		y4=commands.getoutput('pget dmstat out_median')
		sy4=commands.getoutput('pget dmstat out_sigma')

		if len(obsids) > 1 :

			os.system('punlearn merge_obs')

			evtstr = ''
			asolstr= ''
			for j in range(len(obsids)) :
				evtstr = evtstr + ','+locclu+'/repro3_evt_nophas'+obsids[j]+'.fits'
				asolstr= asolstr+ ','+asol6_loc.asol(obsids[j], 'asol')
				os.system('dmcopy "mfe_'+clu+'/repro3_evt'+obsids[j]+'.fits[ccd_id=0,1,2,3,5,6,7][col -phas]" '+locclu+'/repro3_evt_nophas'+obsids[j]+'.fits')

			evtstr = evtstr[1:]
			asolstr= asolstr[1:]

			os.system('pset merge_obs infiles="'+evtstr+'"')
			#os.system('pset merge_all chip="0,1,2,3,5,6,7" \n')

			os.system('pset merge_obs refcoord=mfe_'+clu+'/repro3_evt'+longest+'.fits')
			os.system('pset merge_obs outroot='+locclu+'/mrgevt') 
			os.system('pset merge_obs asolfiles='+asolstr)
			os.system('merge_obs mode=h ')

			os.system('cp '+locclu+'/mrgevt_merged_evt.fits '+locclu+'/mrgevt.fits')

			#mrg_cmd_file.close()

			#os.system('source temp_merge_all.csh')
			#os.remove('temp_merge_all.csh')

		elif len(obsids) == 1 : shutil.copy('mfe_'+clu+'/repro3_evt'+obsids[0]+'.fits',  locclu+'/mrgevt.fits')

		# (xc, yc)
		os.system('punlearn dmstat')
		os.system('dmstat "'+locclu+'/mrgevt.fits[sky=ellipse('+x4+','+y4+','+sx4+','+sy4+',0.)][cols x]" median=yes')
		xc=commands.getoutput('pget dmstat out_median')
		xcf = float(xc)
		os.system('punlearn dmstat')
		os.system('dmstat "'+locclu+'/mrgevt.fits[sky=ellipse('+x4+','+y4+','+sx4+','+sy4+',0.)][cols y]" median=yes')
		yc=commands.getoutput('pget dmstat out_median')
		ycf=float(yc)

		# (rac, decc)
		os.system('punlearn dmcoords')
		os.system('dmcoords '+locclu+'/mrgevt.fits asolfile='+acao+' option=sky x='+xc+' y='+yc+' celfmt=hms')
		rac  = commands.getoutput('pget dmcoords ra')
		decc = commands.getoutput('pget dmcoords dec')		

		# writing rac, dec to file
		racfile = open(locclu+'/rac.txt', 'w')
		racfile.write(rac+'\n')
		racfile.close()
		deccfile = open(locclu+'/decc.txt', 'w')
		deccfile.write(decc+'\n')
		deccfile.close()

	else :
		print
		print '------- Merged file exists already, skipping mrg procedure'

	print
	print

