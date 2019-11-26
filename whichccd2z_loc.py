
# PURPOSE: return all CCDs that have enough data in a given observation

# v2: add most of the code from radreg4 which is used to make region file

# whichccd2z : * based on whichccd2
#              * without the making of obOBccdN.reg

# OUTPUT: return value is presccdglob

import numpy

import commands
import os

import centradec
import asol6_loc

def find(clu,ob) :

	locclu = 'mfe_'+clu
	(rac, decc) = centradec.getrd(locclu)

	acao = asol6_loc.asol(ob,'asol')

	ccdname='i567'
	ccdlist=['0,1,2,3','5','6','7']

	# Get the SKY coodinates of the center.
	os.system('punlearn dmcoords')
	os.system('dmcoords '+locclu+'/clean'+ob+'.fits asolfile="'+acao+'" option=cel ra='+rac+' dec='+decc+' celfmt=hms')
	xc = commands.getoutput('pget dmcoords x')
	yc = commands.getoutput('pget dmcoords y')
	xcf = float(xc)
	ycf = float(yc)

	# ra, dec in degrees
	os.system('punlearn dmcoords')
	os.system('dmcoords '+locclu+'/clean'+ob+'.fits asolfile="'+acao+'" option=sky x='+xc+' y='+yc+' celfmt=deg')
	racd = float(commands.getoutput('pget dmcoords ra'))
	decd = float(commands.getoutput('pget dmcoords dec'))

	# Make a list of point sources, to be able to choose which ones are close
	# to a given shell for it to be excluded.
	ptfile = open(locclu+'/reg/pt0mfe_wcs.reg', 'r')
	regline = ptfile.readline()
	regv = []
	ptr=[]  # how far from (rac,decc) in pixels
	ptdr=[] # 'radius' of circle/ellipse in pixels
	while regline != '\n' and regline != '' :

		if regline[0] != '#' :
			
			regv.append(regline)
			
			s1=regline.partition('(')

			s2=s1[2].partition(':')
			hh=float(s2[0])
			s3=s2[2].partition(':')
			mm=float(s3[0])
			s4=s3[2].partition(',')
			ss=float(s4[0])

			ptra = 15.*(hh+mm/60.+ss/3600.)

			s5=s4[2].partition(':')
			dd=float(s5[0])
			if s5[0][0]=='+' : ff=1.
			elif s5[0][0]=='-' : ff=-1.
			s6=s5[2].partition(':')
			am=float(s6[0])
			s7=s6[2].partition(',')
			ac=float(s7[0])

			ptdec= dd + ff * (am/60. + ac/3600.)

			d2 = ( ( (ptra-racd) * numpy.cos(decd/180.*numpy.pi) ) ** 2. + (ptdec-decd) ** 2. ) * (3600./0.492)**2.

			ptr.append( d2**0.5 )

			if s1[0] == 'ellipse' :
				s8=s7[2].partition("',")
				s9=s8[2].partition("',")
				ptdr.append(60./0.492*max(float(s8[0]),float(s9[0])))
			elif s1[0] == 'circle' :
				s8=s7[2].partition("'")
				ptdr.append(60./0.492*float(s8[0]))
			else : 
				print
				print 'BAAAAAAAAAAAAAD REGIOOOOOOOOOOOON'
				print

		regline = ptfile.readline()
	ptfile.close()


	# get CCD regions
	chipx=commands.getoutput('dmlist "'+locclu+'/reg/chipsreg'+ob+'.fits[cols x]" data,clean')
	t1 = chipx.split('\n')
	t1.remove(t1[0])
	chipy=commands.getoutput('dmlist "'+locclu+'/reg/chipsreg'+ob+'.fits[cols y]" data,clean')
	u1 = chipy.split('\n')
	u1.remove(u1[0])
	ccdid=commands.getoutput('dmlist "'+locclu+'/reg/chipsreg'+ob+'.fits[cols ccd_id]" data,clean')
	c1 = ccdid.split('\n')
	c1.remove(c1[0])

	# goodchips[0] contains indices to u1, t1 or c1 that pertain to chip I (ie ACIS-I)
	# goodchips[1] contains indices to u1, t1 or c1 that pertain to chip 5
	# goodchips[2] contains indices to u1, t1 or c1 that pertain to chip 6
	# goodchips[3] contains indices to u1, t1 or c1 that pertain to chip 7
	goodchips=[[],[],[],[]]
	#           I  5  6  7
	for i in range(len(c1)) : 
		if int(c1[i]) in (0,1,2,3) : goodchips[0].append(i)
		if int(c1[i]) == 5 : goodchips[1].append(i)
		if int(c1[i]) == 6 : goodchips[2].append(i)
		if int(c1[i]) == 7 : goodchips[3].append(i)



	presccdglob=''
	for iccd in range(4) :
		if len(goodchips[iccd]) > 0 and os.path.exists(locclu+'/bg/bg'+ob+'_ccd'+ccdname[iccd]+'r_en.fits') : presccdglob=presccdglob+ccdname[iccd]

	return presccdglob









