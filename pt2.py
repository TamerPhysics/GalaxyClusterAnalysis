
# PURPOSE: Detect point sources with the wavdetect tool

# v2: for imlr.py: some name changes + making mrgimg_en.fits files

# OUTPUT: reg/pt0mfe.reg
# OUTPUT: mrgimg_en.fits

import os
import commands

def pt(clu) :     #, obsids) : 

	print
	print 'PTsrc finder for cluster '+clu+' ---------------------------'
	print

	locclu = 'mfe_' + clu

	if not os.path.exists(locclu+'/reg') : os.mkdir(locclu+'/reg')

	xs=commands.getoutput('dmlist "'+locclu+'/mrgevt.fits[cols x]" data,clean')
	xsv=xs.split('\n')
	xsv.remove(xsv[0])
	for ii in range(len(xsv)) : 
		xsv[ii]=float(xsv[ii])

	ys=commands.getoutput('dmlist "'+locclu+'/mrgevt.fits[cols y]" data,clean')
	ysv=ys.split('\n')
	ysv.remove(ysv[0])
	for ii in range(len(ysv)) : 
		ysv[ii]=float(ysv[ii])

	xmin = min(xsv)
	xmax = max(xsv)
	ymin = min(ysv)
	ymax = max(ysv)

	os.system('punlearn dmcopy')
	os.system('dmcopy "'+locclu+'/mrgevt.fits[energy=300:7000][bin x='+str(xmin)+':'+str(xmax)+':2,y='+str(ymin)+':'+str(ymax)+':2]" '+locclu+'/mrgimg_en.fits')

	os.system('rm -f wavzebala/*')

	os.system('punlearn wavdetect')
	os.system('pset wavdetect infile='+locclu+'/mrgimg_en.fits')
	os.system('pset wavdetect outfile=temp_wav_out.fits')
	os.system('pset wavdetect regfile='+locclu+'/reg/pt0mfe.reg')
	os.system('pset wavdetect scellfile=temp_scell.fits')
	os.system('pset wavdetect imagefile=temp_img.fits')
	os.system('pset wavdetect defnbkgfile=temp_bg.fits')
	os.system('pset wavdetect scales="2.0 4.0"')
	os.system('pset wavdetect ellsigma=3')
	os.system('pset wavdetect interdir=wavzebala')
	os.system('wavdetect mode=h')

	# clear useless files:
	if os.path.exists('temp_scell.fits') :    os.remove('temp_scell.fits')
	if os.path.exists('temp_img.fits') :      os.remove('temp_img.fits')
	if os.path.exists('temp_bg.fits') :       os.remove('temp_bg.fits')
	if os.path.exists('temp_wav_out.fits') :  os.remove('temp_wav_out.fits')




