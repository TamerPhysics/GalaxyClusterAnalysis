
# v3 : * changed pixel size from 0.496" to 0.492", which is correct value.
#      * use asol3.py to find asol, flt1, pbk, msk, bpix or evt1 files

# v5 : change to asol6_loc and hardsoft3

# PURPOSE: call radial-profile-making codes, and merge ACIS-S observations and ACIS-I observations seperately

# OUTPUT: clean_nophasOBSID_en.fits 
# OUTPUT: mrgclean_en.fits

import commands
import glob
import os
import numpy
import shutil

import radial10
import hardsoft3
import centradec
import asol6_loc

def preradial(clu, obsids, r1mpc) :

	locclu = 'mfe_'+clu

	exptime = numpy.zeros(len( obsids ))
	for j in range(len( obsids )) :
		exptime[j] = float ( commands.getoutput('dmkeypar '+locclu+'/clean'+obsids[j]+'.fits livetime echo+') )


	longob = obsids[numpy.argmax(exptime)]

	if not os.path.exists(locclu+'/ctofr') : os.mkdir(locclu+'/ctofr')

	if len(obsids) > 1 : 

		os.system('punlearn obs')

		evtstr = ''
		asolstr= ''
		for j in range(len(obsids)) :
			evtstr = evtstr + ','+locclu+'/clean_nophas'+obsids[j]+'_en.fits'
			asolstr= asolstr+ ','+asol6_loc.asol(obsids[j], 'asol')
			os.system('punlearn dmcopy')
			os.system('dmcopy "'+locclu+'/clean'+obsids[j]+'.fits[energy=300:7000][ccd_id=0,1,2,3,5,6,7][col -phas]" '+locclu+'/clean_nophas'+obsids[j]+'_en.fits')
		evtstr = evtstr[1:]
		asolstr= asolstr[1:]

		os.system('pset merge_obs infiles="'+evtstr+'"')
		os.system('pset merge_obs asolfiles='+asolstr)
		os.system('pset merge_obs refcoord='+locclu+'/clean'+longob+'.fits')
		os.system('pset merge_obs outroot='+locclu+'/mrgclean_en') 
		os.system('merge_obs mode=h ')

		os.system('cp '+locclu+'/mrgclean_en_merged_evt.fits '+locclu+'/mrgclean_en.fits')

		# clean up
		os.system('rm '+locclu+'/clean_nophas[1-9]*_en.fits')

	else : 
		os.system('punlearn dmcopy')
		os.system('dmcopy "'+locclu+'/clean'+obsids[0]+'.fits[ccd_id=0,1,2,3,5,6,7][energy=300:7000][col -phas]" '+locclu+'/mrgclean_en.fits')




	acao=asol6_loc.asol(longob, 'asol')


	# get the coordinates of the center in SKY units
	os.system('punlearn dmcoords')
	os.system('dmcoords '+locclu+'/mrgclean_en.fits asolfile='+acao+' ra=`more '+locclu+'/rac.txt` dec=`more '+locclu+'/decc.txt` option=cel celfmt=hms')
	xs = commands.getoutput('pget dmcoords x')
	ys = commands.getoutput('pget dmcoords y')
	centerccd = int(commands.getoutput('pget dmcoords chip_id'))

	# Make file with the radial coordinate of each photon count
	os.system('punlearn dmtcalc')
	os.system('dmtcalc "'+locclu+'/mrgclean_en.fits[col -time,-ccd_id,-node_id,-expno,-chip,-tdet,-det,-phas,-pha_ro,-energy,-pi,-fltgrade,-grade,-status][exclude sky=region('+locclu+'/reg/pt0mfe_wcs.reg)]" "'+locclu+'/evtrad_mrgclean.txt[opt kernel=text/simple]" expr="r2=(('+xs+'-sky[0])^2)+(('+ys+'-sky[1])^2)" clobber=no')

	r2=[]
	rfile = open(locclu+'/evtrad_mrgclean.txt', 'r')
	intermed = rfile.readline()
	intermed = rfile.readline()
	intermed = rfile.readline() # first data-containing line
	while intermed != '' and intermed != '\n' :
		r2.append( float(intermed.split(' ')[2][0:-1]) )
		intermed = rfile.readline()
	rfile.close()

	rmax = (max(r2))**0.5

	os.system( 'rm '+locclu+'/radbins5.txt'+' '+locclu+'/hienrlo.txt' )
	for ob in obsids : radial10.radial(clu, ob, r1mpc/0.492, rmax, obsids)

