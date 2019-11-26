
# v2: for ktofr9
# v3: exclude from presccdglob CCD where bg/bgOB_ccdCCDr_en.fits doesn't exist
# v4: * only return presccd, not presccdglob
#     * remove the making of reg/obOBccdN.reg, since it is now made in whichccd2

# v4simple: * use simple regions simplechipsregOB.reg
#           * make reg/obsimple*.reg

# v4simple2: * make the _xfov_ptsrc.reg files even for OBSIDS in badskyfovob

# PURPOSE: this function returns the ccd's that are present for this shell,
# PURPOSE: and makes ccd regions AND shell regions

# INPUT: rlo, rhi in pixels

# OUTPUT: * return value: a string with the CCD's which have counts in the provided
#           range [rlo,rhi]
# OUTPUT: * prefix+'xfov_pt_simple.reg'

import numpy

import commands
import os

import centradec
import asol6_loc

import pdb

def radreg(clu, ob, rlo, rhi, prefix) :

	locclu = 'mfe_'+clu
	(rac, decc) = centradec.getrd(locclu)

	badskyfovob = ['3182', '897', '11708']

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
	chipx=commands.getoutput('dmlist "'+locclu+'/reg/simplechipsreg'+ob+'.fits[cols x]" data,clean')
	t1 = chipx.split('\n')
	t1.remove(t1[0])
	chipy=commands.getoutput('dmlist "'+locclu+'/reg/simplechipsreg'+ob+'.fits[cols y]" data,clean')
	u1 = chipy.split('\n')
	u1.remove(u1[0])
	ccdid=commands.getoutput('dmlist "'+locclu+'/reg/simplechipsreg'+ob+'.fits[cols ccd_id]" data,clean')
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

	# ccdreg will contain the regions (in string format) corresponding to
	# the CCD's I, 5, 6 and 7
	ccdreg=[[],[],[],[]]
	ccdregfileappend=[False,False,False,False]
	for iccd in range(4) :

		if len(goodchips[iccd]) > 0 :
			for igood in range(len(goodchips[iccd])) :

				i = goodchips[iccd][igood]

				polyx = []
				t2=t1[i]
				t3=t2.lstrip()
				t4=t3.partition(' ')
				while ('NaN' not in t4[0]) and (t4[1] != '') :
					polyx.append(float(t4[0]))
					t2=t4[2]
					t3=t2.lstrip()
					t4=t3.partition(' ')

				polyy = []
				u2=u1[i]
				u3=u2.lstrip()
				u4=u3.partition(' ')
				while ('NaN' not in u4[0]) and (u4[1] != '') :
					polyy.append(float(u4[0]))
					u2=u4[2]
					u3=u2.lstrip()
					u4=u3.partition(' ')

				ccdreg[iccd].append('polygon(')
				for j in range(len(polyx)) : 
					ccdreg[iccd][igood] = ccdreg[iccd][igood] + str(polyx[j]) + ',' + str(polyy[j]) + ','
				ccdreg[iccd][igood] = ccdreg[iccd][igood][0:-1]
				ccdreg[iccd][igood] = ccdreg[iccd][igood] + ')'

			if not os.path.exists(locclu+'/reg/obsimple'+ob+'ccd'+ccdname[iccd]+'.reg') and \
                        ob not in badskyfovob : # <--- list of obsids where ccd region was created mannually, because
                                                #      skyfov returned wrong region
				ccdregfile = open(locclu+'/reg/obsimple'+ob+'ccd'+ccdname[iccd]+'.reg', 'w')
				for igood in range(len(goodchips[iccd])) : ccdregfile.write(' +'+ccdreg[iccd][igood])
				ccdregfile.close()
				ccdregfileappend[iccd]=True

	# Hi BG corner region
	if os.path.exists(locclu+'/reg/hibgcorner'+ob+'.reg') :
		hibgfile = open(locclu+'/reg/hibgcorner'+ob+'.reg','r')
		hibgstr = hibgfile.readline()
		hibgfile.close()
		hibgx = float( hibgstr[7:-1].split(',')[0] )
		hibgy = float( hibgstr[7:-1].split(',')[1] )
		hibgr = ( (hibgx-xcf)**2. + (hibgy-ycf)**2. )**0.5


	# make region of intersection of annulus with chip regions MINUS PTSRCs
	presccd=''
	for iccd in range(4) : 
		if len(ccdreg[iccd]) > 0 :

			# write region to file
			if ob not in badskyfovob :
				fovptfile = open(prefix+'_ccd'+str(ccdname[iccd])+'_xfov_pt_simple.reg', 'w')
				for ii in range(len(ccdreg[iccd])) : 
					fovptfile.write('+annulus('+xc+','+yc+','+str(rlo)+','+str(rhi)+') * '+ccdreg[iccd][ii]+'\n')
					for jj in range(len(ptr)) : # add point source regions which lie within this annulus
						if ptr[jj]-ptdr[jj] <= rhi and ptr[jj]+ptdr[jj] >= rlo : fovptfile.write('-'+regv[jj][0:-1]+' ')
					if iccd == 0 and os.path.exists(locclu+'/reg/hibgcorner'+ob+'.reg') and hibgr-1095. <= rhi and hibgr+1095. >= rlo : fovptfile.write( '-'+hibgstr[0:-1] )
					fovptfile.write('\n')
				fovptfile.close()

			elif os.path.exists(locclu+'/reg/obsimple'+ob+'ccd'+ccdname[iccd]+'.reg') :
				fovptfile = open(prefix+'_ccd'+str(ccdname[iccd])+'_xfov_pt_simple.reg', 'w')
				for ii in range(len(ccdreg[iccd])) : 
					obsimplefile = open(locclu+'/reg/obsimple'+ob+'ccd'+ccdname[iccd]+'.reg', 'r')
					obsimpleline = obsimplefile.readline().replace('\n','').replace(' +',' ')
					fovptfile.write('+annulus('+xc+','+yc+','+str(rlo)+','+str(rhi)+') * '+obsimpleline+'\n')
					for jj in range(len(ptr)) : # add point source regions which lie within this annulus
						if ptr[jj]-ptdr[jj] <= rhi and ptr[jj]+ptdr[jj] >= rlo : fovptfile.write('-'+regv[jj][0:-1]+' ')
					if iccd == 0 and os.path.exists(locclu+'/reg/hibgcorner'+ob+'.reg') and hibgr-1095. <= rhi and hibgr+1095. >= rlo : fovptfile.write( '-'+hibgstr[0:-1] )
					fovptfile.write('\n')
					obsimplefile.close()
				fovptfile.close()

			# add a good CCD to presccd (radial shell dependent
			if float(commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[ccd_id='+ccdlist[iccd]+'][sky=region('+prefix+'_ccd'+str(ccdname[iccd])+'_xfov_pt_simple.reg)]" counts')) > 0 and os.path.exists(locclu+'/bg/bg'+ob+'_ccd'+ccdname[iccd]+'r_en.fits') : presccd=presccd+ccdname[iccd]
			else : 
				if ob not in badskyfovob : os.remove(prefix+'_ccd'+str(ccdname[iccd])+'_xfov_pt_simple.reg')

	# presccd is a string containing the CCD's found in this range rlo -> rhi
	return presccd

