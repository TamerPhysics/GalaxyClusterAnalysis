# v4: radial bins are now equally spaced on log scale. See 'HOW' section for details.
#     changed name of output files by adding a '2'
#     statistical error on surface brightness = 0.2 * surface brightness
#                                             = sqrt [ (0.2 * surface brightness)^2 + poisson^2 ]
# v5: take into account the CCD areas in computing solid angle of each shell, by
#     dmextracting with anNxfov_pt3.reg (instead of anN_2.reg)
# v6: * exclude hi-BG corner region
#     * in the making of radial bins, if first bin is empty of counts, don't expand the outer
#       radius, but move the starting inner radius further away from the center
#     * make an ARF for each radial bin! -- WAS NOT APPLIED TO DATA!
#     * add WMAP to spectra to be used for the making of weighted ARF -- WAS NOT APPLIED TO DATA!
# v7: * applied ARF changes made above which were not applied to data in v6
#     * pixel size changed from 0.496" to 0.492"
#     * because the size of the pixel was changed from 0.496" to 0.492", the dimensions
#       in pixel in all the below files will change, and so they are given new names.
#       This time all the files are numbered '4' to easily recognise them : 
#        - radbins_hard.txt                      --> radbins_hard4.txt
#        - radpOB/spN_3.fits                     --> radpOB/spN_4.fits
#        - radpOB/anNxfov2.reg                   --> radpOB/anNxfov4.reg
#        - radpOB/anN_2.reg                      --> radpOB/anN_4.reg
#        - radpOB/anNxfov_pt3.reg                --> radpOB/anNxfov_pt4.reg
#        - radpOB/spN_3.fits                     --> radpOB/spN_4.fits
#        - ctofr/radprofintOB_hard3.txt          --> ctofr/radprofintOB_hard4.txt
#        - ctofr/radprofintOB_hard3pois.txt      --> ctofr/radprofintOB_hard4pois.txt
#        - ctofr/radprofintpixOB_hard3.txt       --> ctofr/radprofintpixOB_hard4.txt
#        - ctofr/radprofintpixOB_hard3pois.txt   --> ctofr/radprofintpixOB_hard4pois.txt
#        - ctofr/radprofintOB_hard3err2.txt      --> ctofr/radprofintOB_hard4err2.txt
#        - ctofr/omegaOB_hard3.txt               --> ctofr/omegaOB_hard4.txt
#     * [energy=700:2000] filter for counts extraction, because this energy better represents
#       surface brightness
#     * use asol3.py to find asol, flt1, pbk, msk, bpix or evt1 files
#     * check that area != 0 before division, to avoid zero division
# v8: * generalize to include any CCD in any OBSID
#     * WMAP in spISHELL.fits from 0.7-2keV, instead of 0.3-2keV
#     * add 4rd parameter: rmax 
# v9: * make losrc*.reg similar to ktofrsp/*_xfov_pt.reg by adding more return characters
# v10: * change to asol6_loc

# PURPOSE: make radial surface brightness profile

# HOW:
#  * radial bins are equally spaced in log space, from rmin=(1.2"/.492)pix to rmax=min(r1mpcpix,1420) :
#    - for r<200kpc, r(i+1)=1.25  *r(i)
#    - UNLESS a radial bin has less than 100counts, in mrgclean_en.fits, then the next bin is added to it.
#  * making the region files corresponding to radial bins: get point source position and distance to center (ptr)
#    from each annulus region, subtract the relevant ptsrc's. this will make filtering evt files faster.

# INPUT: rac.txt and decc.txt (through call of centradec function)
# INPUT: cleanOBS.fits
# INPUT: evtOB_c.fits (energy=[700:2000], with all chips)
# INPUT: instmap/expmapOB_a.fits
# INPUT: reg/pt0_mfe.reg
# INPUT: bg/bgOBr.fits

# OUTPUT: radbins5.txt
# OUTPUT: radpOB/anN_ccdM_xfov.reg 
# OUTPUT: radpOB/anN_ccdM_xfov_pt.reg
# OUTPUT: radpOB/spN_ccdM.fits 
# OUTPUT: radpOB/warfN_ccdM.fits through arf4radp3.py
# OUTPUT: ctofr/srcprofOB_ccdN.txt
# OUTPUT: ctofr/bkgprofOB_ccdN.txt
# OUTPUT: ctofr/omegaOB_ccdM.txt
# OUTPUT: bg/bgOBr_en.fits
# OUTPUT: hienOB.fits: cleanOB.fits filtered in energy=9500:12000
# OUTPUT: hienrlo.txt

import commands
import os
import numpy
import pdb

import asol6_loc
import centradec
import arf4radp4
import xeqy

def radial(clu, ob, r1mpcpix, rmax, obsids) : 

	# obsids param only used to get bg count to compare with counts from merged file.
	# all operations done on ob, which is one element of obsids.

	locclu= 'mfe_'+clu

	if not os.path.exists(locclu+'/radp'+ob) : os.mkdir(locclu+'/radp'+ob)

	# create high-energy photons only event file, to use to estimate BG level
	os.system('punlearn dmcopy')
	os.system('dmcopy "'+locclu+'/clean'+ob+'.fits[energy=9500:12000]" '+locclu+'/hien'+ob+'.fits')

	# center of the cluster
	acao = asol6_loc.asol(ob,'asol')
	(rac, decc) = centradec.getrd(locclu)
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

	# read all point source areas to subtract them later from observed area
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

	# get chip regions
	chipx=commands.getoutput('dmlist "'+locclu+'/reg/chipsreg'+ob+'.fits[cols x]" data,clean')
	t1 = chipx.split('\n')
	t1.remove(t1[0])

	chipy=commands.getoutput('dmlist "'+locclu+'/reg/chipsreg'+ob+'.fits[cols y]" data,clean')
	u1 = chipy.split('\n')
	u1.remove(u1[0])

	# find chips 0,1,2,3
	ccdid=commands.getoutput('dmlist "'+locclu+'/reg/chipsreg'+ob+'.fits[cols ccd_id]" data,clean')
	c1 = ccdid.split('\n')
	c1.remove(c1[0])

	goodchips=[[],[],[],[]]
	#           I  5  6  7
	for i in range(len(c1)) : 
		if int(c1[i]) in (0,1,2,3) : goodchips[0].append(i)
		if int(c1[i]) == 5 : goodchips[1].append(i)
		if int(c1[i]) == 6 : goodchips[2].append(i)
		if int(c1[i]) == 7 : goodchips[3].append(i)

	ccdreg=[[],[],[],[]]
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

	ccdname=['i','5','6','7']
	for iccd in range(4) :
		if len(ccdreg[iccd]) > 0 :
			os.system('punlearn dmcopy')
			os.system('dmcopy "'+locclu+'/bg/bg'+ob+'_ccd'+ccdname[iccd]+'r.fits[energy=300:7000]" '+locclu+'/bg/bg'+ob+'_ccd'+ccdname[iccd]+'r_en.fits')
			if not os.path.exists(locclu+'/bg/bg'+ob+'_ccd'+ccdname[iccd]+'r_en.fits') : ccdreg[iccd]=[]

	exposure = float(commands.getoutput('dmkeypar '+locclu+'/evt'+ob+'_c.fits exposure echo+'))


	# Making radbins5.txt
	vrlo=[]
	vrhi=[]
	if not os.path.exists(locclu+'/radbins5.txt') :

		rmin = 1.2 / .492 # in pixels.  1.2" = angular size of 10kpc at z=1

		rlo=rmin
		rhi=rmin*1.25
		hienrlo1=-1
		hienrlo2=-1
		hienrlo3=-1
		ct1=0
		ct2=0
		ct3=0
		print 'rmax ', rmax
		while rhi <= rmax :

			# make annulus-minus-ptsrc region for the shell between rlo and rhi.
			regfile = open(locclu+'/an_test.reg', 'w')
			regfile.write('annulus('+xc+','+yc+','+str(rlo)+','+str(rhi)+')\n')
			for ii in range(len(ptr)) :
				if ptr[ii]-ptdr[ii] <= rhi and ptr[ii]+ptdr[ii] >= rlo : regfile.write('-'+regv[ii])
			regfile.close()

			ct=float(commands.getoutput('dmlist "'+locclu+'/mrgclean_en.fits[sky=region('+locclu+'/an_test.reg)]" counts'))

			sclbgct=0.
			for oobb in obsids :
				for ccdn in ['i','5','6','7'] :

					if os.path.exists(locclu+'/bg/bg'+oobb+'_ccd'+ccdn+'r_en.fits') : 
						sclbgct = sclbgct + float(commands.getoutput('dmkeypar '+locclu+'/clean'+oobb+'.fits exposure echo+')) / float(commands.getoutput('dmkeypar '+locclu+'/bg/bg'+oobb+'_ccd'+ccdn+'r_en.fits exposure echo+')) * float(commands.getoutput('dmlist "'+locclu+'/bg/bg'+oobb+'_ccd'+ccdn+'r_en.fits[sky=region('+locclu+'/an_test.reg)]" counts'))
						

			shellct = ct - sclbgct

			os.remove(locclu+'/an_test.reg')
			print rlo , rhi , shellct

			# shell has enough counts, OR we've reached rmax
			if shellct > 100 or rhi > 0.99 * rmax :
				vrlo.append(rlo)
				vrhi.append(rhi)

				if rhi > 0.99 * rmax : break

				# radius where we start measuring hi-energy counts
				if ct <= 1.5 * sclbgct : 
					if ct1 == 0 : hienrlo1 = float(rlo) # first choice
					ct1=ct1+int(ct)
				if ct <= 2.0 * sclbgct : 
					if ct2 == 0 : hienrlo2 = float(rlo) # second choice, in case 1st doesnt exist
					ct2=ct2+int(ct)
				if ct <= 3.0 * sclbgct : 
					if ct3 == 0 : hienrlo3 = float(rlo) # 3rd choice
					ct3=ct3+int(ct)

				rlo=float(rhi)
				rhi=rhi*1.25

			else :
				if ct == 0 and len(vrlo) == 0 : rlo=float(rhi)
				rhi=rhi*1.25

			if rhi > rmax : rhi=rmax



		hienrlofile = open(locclu+'/hienrlo.txt','w')
		hienrlofile.write('1.5 '+str(hienrlo1)+' '+str(ct1)+'\n')
		hienrlofile.write('2 '+str(hienrlo2)+' '+str(ct2)+'\n')
		hienrlofile.write('3 '+str(hienrlo3)+' '+str(ct3)+'\n')
		hienrlofile.close()


		binfile = open(locclu+'/radbins5.txt','w')
		for ir in range(len(vrlo)) : 
			binfile.write(str(vrlo[ir])+' '+str(vrhi[ir])+'\n')
		binfile.close()

		xinfile = open(locclu+'/ctofr/xin.txt','w')
		for ir in range(len(vrlo)) : 
			xinfile.write(str(vrlo[ir]/r1mpcpix)+'\n')
		xinfile.close()

		xoutfile = open(locclu+'/ctofr/xout.txt','w')
		for ir in range(len(vrhi)) : 
			xoutfile.write(str(vrhi[ir]/r1mpcpix)+'\n')
		xoutfile.close()


	# the case where radbins5.txt has been created in a previous ob:
	else : 

		# Loading radial bin edges, from radbins5.txt
		binfile=open(locclu+'/radbins5.txt', 'r')
		linestr = binfile.readline()
		while linestr != '' :
			intermed = linestr.split(' ')
			vrlo.append( float(intermed[0]) )
			vrhi.append( float(intermed[1][0:-1]) )
			linestr = binfile.readline()
		binfile.close()


	# Hi BG corner region
	if os.path.exists(locclu+'/reg/hibgcorner'+ob+'.reg') :
		hibgfile = open(locclu+'/reg/hibgcorner'+ob+'.reg','r')
		hibgstr = hibgfile.readline()
		hibgfile.close()
		hibgx = float( hibgstr[7:-1].split(',')[0] )
		hibgy = float( hibgstr[7:-1].split(',')[1] )
		hibgr = ( (hibgx-xcf)**2. + (hibgy-ycf)**2. )**0.5

	hien=[[] ,[] ,[] ,[] ]
	bghien=[[] ,[] ,[] ,[] ]
	hienomega=[[] ,[] ,[] ,[] ]
	ccdlist=['0,1,2,3','5','6','7']
	scl=[-1, -1, -1, -1]
	for iccd in range(4) :

		if len(ccdreg[iccd]) > 0 :

			scl[iccd] = exposure / float(commands.getoutput('dmkeypar "'+locclu+'/bg/bg'+ob+'_ccd'+ccdname[iccd]+'r_en.fits" exposure echo+'))

			radintfile2err = open(locclu+'/ctofr/srcprof'+ob+'_ccd'+ccdname[iccd]+'.txt','w') # poisson and systematic errors added in quadrature
			radintfile2err.write("# R1 R2 SB SBERR\n")
			bgfile = open(locclu+'/ctofr/bkgprof'+ob+'_ccd'+ccdname[iccd]+'.txt','w')
			bsfile = open(locclu+'/ctofr/omega'+ob+'_ccd'+ccdname[iccd]+'.txt','w')

			bgexp = float(commands.getoutput('dmkeypar '+locclu+'/bg/bg'+ob+'_ccd'+ccdname[iccd]+'r_en.fits exposure echo+'))

			for i in range(len(vrlo)) :

				# make region of intersection of annulus with chip regions
				fovfile = open(locclu+'/radp'+ob+'/an'+str(i)+'_ccd'+ccdname[iccd]+'_xfov.reg', 'w')
				for ii in range(len(ccdreg[iccd])) : 
					fovfile.write('+annulus('+xc+','+yc+','+str(vrlo[i])+','+str(vrhi[i])+') * '+ccdreg[iccd][ii]+'\n')
				fovfile.close()

				# make region of intersection of annulus with chip regions MINUS PTSRCs
				fovptfile = open(locclu+'/radp'+ob+'/an'+str(i)+'_ccd'+ccdname[iccd]+'_xfov_pt.reg', 'w')
				for ii in range(len(ccdreg[iccd])) : 
					fovptfile.write('+annulus('+xc+','+yc+','+str(vrlo[i])+','+str(vrhi[i])+') * '+ccdreg[iccd][ii]+'\n')
					for jj in range(len(ptr)) :
						if ptr[jj]-ptdr[jj] <= vrhi[i] and ptr[jj]+ptdr[jj] >= vrlo[i] : fovptfile.write('-'+regv[jj][0:-1]+' ')
					if os.path.exists(locclu+'/reg/hibgcorner'+ob+'.reg') and hibgr-1095. <= vrhi[i] and hibgr+1095. >= vrlo[i] : fovptfile.write( '-'+hibgstr[0:-1] )
					fovptfile.write('\n')
				fovptfile.close()

				# calculate counts/area/solidangle/second
				ct=float(commands.getoutput('dmlist "'+locclu+'/evt'+ob+'_c.fits[sky=region('+locclu+'/radp'+ob+'/an'+str(i)+'_ccd'+ccdname[iccd]+'_xfov_pt.reg)]" counts')) 
				bgct=float(commands.getoutput('dmlist "'+locclu+'/bg/bg'+ob+'_ccd'+ccdname[iccd]+'r_en.fits[energy=700:2000][sky=region('+locclu+'/radp'+ob+'/an'+str(i)+'_ccd'+ccdname[iccd]+'_xfov_pt.reg)]" counts')) # more exact than above (deleted) line, but could be slower if many ptsrcs
				#netct = ct - bgct * scl[iccd]

				# make spectrum of annulus 
				os.system('punlearn dmextract')
				os.system('dmextract "'+locclu+'/clean'+ob+'.fits[ccd_id='+ccdlist[iccd]+'][sky=region('+locclu+'/radp'+ob+'/an'+str(i)+'_ccd'+ccdname[iccd]+'_xfov_pt.reg)][bin PI]" '+locclu+'/radp'+ob+'/sp'+str(i)+'_ccd'+ccdname[iccd]+'.fits opt=pha1 wmap="[energy=700:2000][bin det=8]"')

				# calculate hi-energy SB
				hien[iccd].append( float(commands.getoutput('dmlist "'+locclu+'/hien'+ob+'.fits[ccd_id='+ccdlist[iccd]+'][sky=region('+locclu+'/radp'+ob+'/an'+str(i)+'_ccd'+ccdname[iccd]+'_xfov_pt.reg)]" counts')) )
				bghien[iccd].append( float(commands.getoutput('dmlist "'+locclu+'/bg/bg'+ob+'_ccd'+ccdname[iccd]+'r.fits[energy=9500:12000][sky=region('+locclu+'/radp'+ob+'/an'+str(i)+'_ccd'+ccdname[iccd]+'_xfov_pt.reg)]" counts')) )
				hienomega[iccd].append( float(commands.getoutput('dmkeypar '+locclu+'/radp'+ob+'/sp'+str(i)+'_ccd'+ccdname[iccd]+'.fits backscal echo+')) )

				# make ARF for radial bin
				arf4radp4.arf(clu, ob, i, ccdname[iccd])

				# calculate solid angle from BACKSCAL keyword in spectrum file
				bs = float(commands.getoutput('dmkeypar '+locclu+'/radp'+ob+'/sp'+str(i)+'_ccd'+ccdname[iccd]+'.fits backscal echo+'))

				# calculate effective area
				os.system('punlearn dmstat')
				# this if statement handles the case when the annulus thickness is less than a pixel,
				# which causes dmstat to return area=0
				if vrhi[i]-vrlo[i] > 1.5 : 
					if vrhi[i] <= 0.1*r1mpcpix : os.system('dmstat "'+locclu+'/instmap/expmap'+ob+'_ccd'+ccdname[iccd]+'_in.fits[sky=region('+locclu+'/radp'+ob+'/an'+str(i)+'_ccd'+ccdname[iccd]+'_xfov.reg)]" centroid=no')
					else                       : os.system('dmstat "'+locclu+'/instmap/expmap'+ob+'_ccd'+ccdname[iccd]+'_out.fits[sky=region('+locclu+'/radp'+ob+'/an'+str(i)+'_ccd'+ccdname[iccd]+'_xfov.reg)]" centroid=no')
				else : 
					if vrhi[i] <= 0.1*r1mpcpix : os.system('dmstat "'+locclu+'/instmap/expmap'+ob+'_ccd'+ccdname[iccd]+'_in.fits[sky=circle('+xc+','+yc+','+str(vrhi[i])+')]" centroid=no')
					else                       : os.system('dmstat "'+locclu+'/instmap/expmap'+ob+'_ccd'+ccdname[iccd]+'_out.fits[sky=circle('+xc+','+yc+','+str(vrhi[i])+')]" centroid=no')

				try : area = float(commands.getoutput('pget dmstat out_mean'))
				except ValueError : area = 0
				# the above try-exept is for the case, where the expmap needed to calulate area doesnt exist.
				# this can happen when the there is no intersection between the in/out region and a ccd.

				if bs != 0 and area != 0 and ct != 0 :

					sb = ct / bs / area / exposure
					#sberr = ( ct + scl**2.*bgct )**.5 / bs / area / exposure # poisson error
					sberr = ct**.5 / bs / area / exposure # poisson error
					radintfile2err.write(str(vrlo[i]/r1mpcpix)+' '+str(vrhi[i]/r1mpcpix)+' '+str(sb)+' '+str((sberr**2.+sb**2./25.)**0.5)+'\n') # poisson and 1/5 fraction, added in quadrature

					bgsb = bgct / bs / area / bgexp
					bgsberr = bgct**0.5 / bs / area / bgexp
					bgfile.write(str(vrlo[i]/r1mpcpix) + ' ' + str(vrhi[i]/r1mpcpix)+' '+str(bgsb)+' '+str(bgsberr)+'\n')
					bsfile.write(str(vrlo[i]/r1mpcpix) + ' ' + str(bs) + ' ' + str(bs*64.0*1024.**2./(numpy.pi)/(vrhi[i]**2-vrlo[i]**2)) + '\n' )

			radintfile2err.close()
			bgfile.close()
			bsfile.close()


	# write BG scaling factors from hi-energy counts
	# get the 1000 outermost counts, and save their region
	hientot  =[0.,0.,0.,0.]
	bghientot=[0.,0.,0.,0.]
	hienomegatot=[0.,0.,0.,0.]
	sclfile = open(locclu+'/ctofr/scl'+ob+'.txt', 'w')
	for iccd in range(4) :
		if len(ccdreg[iccd]) > 0 :
			irlo=len(vrlo)-1
			while hientot[iccd] < 1e3 and irlo >= 0 :
				hientot[iccd]     =     hientot[iccd]+  hien[iccd][irlo]
				bghientot[iccd]   =   bghientot[iccd]+bghien[iccd][irlo]
				hienomegatot[iccd]=hienomegatot[iccd]+hienomega[iccd][irlo]
				irlo=irlo-1
			irlo=irlo+1

			if bghientot[iccd] != 0 : 

				# write the BG counts and backscales
				sclfile.write(str(scl[iccd])+' '+str(hientot[iccd]/bghientot[iccd])+' '+str(hientot[iccd])+' '+str(bghientot[iccd])+' '+str(hienomegatot[iccd])+' '+str(irlo)+' '+str(vrlo[irlo])+'\n' )

				# make region out of which we estimate the Hi-energy cts in the src dataset
				losrcregfile = open(locclu+'/reg/losrc'+ob+'_ccd'+ccdname[iccd]+'.reg', 'w')
				for polyg in ccdreg[iccd] :
					losrcregfile.write('+annulus('+xc+','+yc+','+str(vrlo[irlo])+','+str(vrhi[-1])+') * '+polyg+'\n')
					for jj in range(len(ptr)) :
						if ptr[jj]-ptdr[jj] <= vrhi[-1] and ptr[jj]+ptdr[jj] >= vrlo[irlo] : losrcregfile.write('-'+regv[jj][0:-1]+' ')
					if iccd==0 and os.path.exists(locclu+'/reg/hibgcorner'+ob+'.reg') and hibgr-1095. <= vrhi[-1] and hibgr+1095. >= vrlo[irlo] : losrcregfile.write( '-'+hibgstr[0:-1] )
					losrcregfile.write('\n')
				losrcregfile.close()

				# extract spectrum from the above region, to get the BACKSCAL value from the spectrum
				os.system('punlearn dmextract')
				os.system('dmextract "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/losrc'+ob+'_ccd'+ccdname[iccd]+'.reg)][bin PI]" '+locclu+'/spec/losrc'+ob+'_ccd'+ccdname[iccd]+'.fits opt=pha1 wmap="[energy=700:2000][bin det=8]"')

			else : sclfile.write(str(scl[iccd])+' '+'NaN '+str(hientot[iccd])+' '+str(bghientot[iccd])+' '+str(hienomegatot[iccd])+' '+str(irlo)+' '+str(vrlo[irlo])+'\n' )

		# if this ccd is not in this obsid
		else : sclfile.write('-1 -1 -1 -1 -1 -1 -1\n')
	sclfile.close()
















