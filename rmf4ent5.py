
# v2: * use asol3.py to find asol, flt1, pbk, msk, bpix or evt1 files
# v3: * based on rmf4entsoft2.py
#     * for both ACIS-I AND ACIS-S

# v5: * 2017-03-31: new version of asol: asol6_loc

# PURPOSE:
# make an RMF evaluated at the center of the cluster without surface brightness weighting

# HOW:
#  get locclu variable (=string containing folder name) and ob string containing OBSID
#  get RA and DEC using centradec function
#  use either mkacisrmf or mkrmf according to keywords in the fits file

# INPUT: rac.txt and decc.txt (through call of centradec function)
# INPUT: cleanOBS.fits
# INPUT: spec/spOB_center.pi

# OUTPUT: spec/spOB_ccdN_center.wrmf


import centradec
import asol6_loc

import commands
import os

def rmf(locclu, ob, ccdstr) : 

	acao = asol6_loc.asol(ob,'asol')
	(rac, decc) = centradec.getrd(locclu)
	caldbfolder = commands.getoutput('echo $CALDB')

	if ccdstr=='i' : ccdlist='0,1,2,3'
	else : ccdlist=ccdstr

	os.system('punlearn dmcoords')
	os.system('dmcoords '+locclu+'/clean'+ob+'.fits '+acao+' ra='+rac+' dec='+decc+' option=cel celfmt=hms mode=h')
	chipx = commands.getoutput('pget dmcoords chipx')
	chipy = commands.getoutput('pget dmcoords chipy')
	centerccd = commands.getoutput('pget dmcoords chip_id')


	if centerccd not in ccdlist :

		os.system('punlearn dmstat')
		os.system('dmstat "'+locclu+'/clean'+ob+'.fits[ccd_id='+ccdlist+'][cols x]" median=yes')
		xx=commands.getoutput('pget dmstat out_median')
		os.system('punlearn dmstat')
		os.system('dmstat "'+locclu+'/clean'+ob+'.fits[ccd_id='+ccdlist+'][cols y]" median=yes')
		yy=commands.getoutput('pget dmstat out_median')

		os.system('punlearn dmcoords')
		os.system('dmcoords '+locclu+'/clean'+ob+'.fits '+acao+' x='+xx+' y='+yy+' option=sky mode=h')
		chipx = commands.getoutput('pget dmcoords chipx')
		chipy = commands.getoutput('pget dmcoords chipy')
		centerccd = commands.getoutput('pget dmcoords chip_id')


	if float(chipx) < 1   : chipx='1'
	if float(chipx) >1024 : chipx='1024'
	if float(chipy) < 1   : chipy='1'
	if float(chipy) >1024 : chipy='1024'


	# RMF -----------------------
	fptemp = float(commands.getoutput('dmkeypar '+locclu+'/clean'+ob+'.fits fp_temp echo+'))
	tgaincor = commands.getoutput('dmkeypar '+locclu+'/clean'+ob+'.fits tgaincor echo+')
	cticor  = commands.getoutput('dmkeypar '+locclu+'/clean'+ob+'.fits cti_corr echo+')
	ctiapp  = commands.getoutput('dmkeypar '+locclu+'/clean'+ob+'.fits cti_app echo+')


	# using mkacisrmf or...
	if (fptemp < 158. or (fptemp >= 158. and ccdstr in ['5','7']) ) and ( ctiapp != 'NNNNNNNNNN' or cticor == 'TRUE' or cticor == '1' ) and (tgaincor == 'T') :

		gainfile  = commands.getoutput('dmkeypar '+locclu+'/clean'+ob+'.fits gainfile echo+')
		gainver   = gainfile.split('.fits')[0][-1]
		caldbls   = commands.getoutput('ls $CALDB/data/chandra/acis/p2_resp/'+gainfile[0:15]+'*'+gainver+'.fits')
		caldbresp = caldbls.split('\n')

		if len(caldbresp) != 1 : 
			print
			print
			print '  ############### MORE OR LESS THAN ONE P2RESP FILES FOUND !!!!!!!!!! ###################'
			print
			print

		os.system('punlearn mkacisrmf')
		os.system('pset mkacisrmf infile='+caldbresp[0])
		os.system('pset mkacisrmf outfile='+locclu+'/spec/sp'+ob+'_ccd'+ccdstr+'_center.wrmf')
		os.system('pset mkacisrmf energy="0.3:11.0:0.01" channel="1:1024:1"')
		os.system('pset mkacisrmf chantype="PI"')
		os.system('pset mkacisrmf wmap=none')
		os.system('pset mkacisrmf chipx='+chipx)
		os.system('pset mkacisrmf chipy='+chipy)
		os.system('pset mkacisrmf ccd_id='+centerccd)
		os.system('pset mkacisrmf gain='+caldbfolder+'/data/chandra/acis/det_gain/'+gainfile)
		os.system('pset mkacisrmf asolfile=none')
		os.system('mkacisrmf mode=h verbose=1')

	# using mkrmf
	else :

		os.system('punlearn acis_fef_lookup')
		os.system('acis_fef_lookup '+locclu+'/clean'+ob+'.fits chipx='+chipx+' chipy='+chipy+' chipid='+centerccd)

		os.system('punlearn mkrmf')
		os.system('pset mkrmf infile='+commands.getoutput('pget acis_fef_lookup outfile') )
		os.system('pset mkrmf outfile='+locclu+'/spec/sp'+ob+'_ccd'+ccdstr+'_center.wrmf')
		os.system('pset mkrmf axis1="energy=0.3:11.0:0.01"')
		os.system('pset mkrmf axis2="pi=1:1024:1"')
		os.system('pset mkrmf weights=none')
		os.system('mkrmf verbose=1 mode=h')


