
# PURPOSE: combinr point source regions with field of view regions, 
# PURPOSE: to create a region with the right observed area

import centradec
import pycrates as cr

import numpy as nn

import os
import pdb

def make(clu, ob, r1mpc, ctsIN0, ctsIN1, ctsIN2, ctsIN3, ctsIN5, ctsIN6, ctsIN7, ctsOUT0, ctsOUT1, ctsOUT2, ctsOUT3, ctsOUT5, ctsOUT6, ctsOUT7) : 

	(rac, decc) = centradec.getrd('mfe_'+clu)


	instr =  'circle('+rac+','+decc+','+str(r1mpc/10)+'")'
	outstr= 'annulus('+rac+','+decc+','+str(r1mpc/10)+'",'+str(r1mpc)+'")'
	bgstr =  'circle('+rac+','+decc+','+str(r1mpc)+'")'


	regfits = cr.read_file('mfe_'+clu+'/reg/chipsreg'+ob+'.fits')

	rx   = cr.copy_colvals(regfits, 'x')
	ry   = cr.copy_colvals(regfits, 'y')

	rccd = cr.copy_colvals(regfits, 'ccd_id')

	# POINT SOURCES WITH A MINUS SIGN FOR EXCLUSION
	filept = open('mfe_'+clu+'/reg/pt0mfe_wcs.reg','r')
	linept=filept.readline()
	linept=filept.readline()
	exclpt=''
	while linept != '' and linept != '\n' :
		exclpt = exclpt + ' -' + linept
		linept=filept.readline()
	filept.close()

	#HI-BG corner region
	hibgbool=False
	if os.path.exists('mfe_'+clu+'/reg/hibgcorner'+ob+'.reg') : 
		hibgbool=True
		filecorner=open('mfe_'+clu+'/reg/hibgcorner'+ob+'.reg', 'r')
		hibgreg=filecorner.readline()
		filecorner.close()

	# combine ptsrcs and hi bg corner
	fileexc = open('mfe_'+clu+'/reg/exclpt'+ob+'.reg', 'w')
	fileexc.write(exclpt)
	if hibgbool : fileexc.write(hibgreg)
	fileexc.close()

	acisiccds='0123'
	# OBSERVATIONS WITH CHIPS 0,1,2,3 -- IN region
	if ctsIN0+ctsIN1+ctsIN2+ctsIN3 >= 40 :

		file1= open( 'mfe_'+clu+'/reg/in2_'+ob+'_ccdi.reg', 'w')

		for ii in range(len(rccd)) :

			if str(rccd[ii]) in acisiccds : 

				sstr = '+polygon('
				for jj in range(rx.shape[1]) :

					if nn.isfinite(rx[ii,jj]) : 
						if jj==0 : sstr = sstr + str(rx[ii,jj]) + ',' + str(ry[ii,jj])
						else     : sstr = sstr + ',' + str(rx[ii,jj]) + ',' + str(ry[ii,jj])

				sstr = sstr + ')'


				file1.write( sstr + ' * ' + instr )
				if hibgbool : file1.write( exclpt + ' -' + hibgreg) 
				else        : file1.write( exclpt )


		file1.close()


	# OBSERVATIONS WITH CHIPS 0,1,2,3 -- OUT/BG regions
	if ctsOUT0+ctsOUT1+ctsOUT2+ctsOUT3 >= 500 :

			file2= open( 'mfe_'+clu+'/reg/out2_'+ob+'_ccdi.reg', 'w')
			file3= open( 'mfe_'+clu+'/reg/bg2_'+ob+'_ccdi.reg', 'w')

			for ii in range(len(rccd)) :

				if str(rccd[ii]) in acisiccds : 

					sstr = '+polygon('
					for jj in range(rx.shape[1]) :

						if nn.isfinite(rx[ii,jj]) : 
							if jj==0 : sstr = sstr + str(rx[ii,jj]) + ',' + str(ry[ii,jj])
							else     : sstr = sstr + ',' + str(rx[ii,jj]) + ',' + str(ry[ii,jj])

					sstr = sstr + ')'


					file2.write( sstr + ' * ' + outstr )
					if hibgbool : file2.write( exclpt + ' -' + hibgreg) 
					else        : file2.write( exclpt )

					file3.write( sstr + ' * ' + bgstr )
					if hibgbool : file2.write( exclpt + ' -' + hibgreg) 
					else        : file2.write( exclpt )


			file2.close()
			file3.close()

	# OBSERVATIONS WITH CHIPS 5, 6 AND 7

	chipv=[5,6,7] # ---> must be integer!!
	regionname=['in2', 'out2', 'bg2']
	regionstr=[instr, outstr, bgstr]
	mincounts=[ 40,    500,   500  ]

	cts = nn.zeros( (len(regionname),len(chipv)) )

	cts[0,0]= ctsIN5
	cts[1,0]= ctsOUT5
	if ctsIN5 >= 40 or ctsOUT5 >= 500 : cts[2,0]= 500+1
	cts[0,1]= ctsIN6
	cts[1,1]= ctsOUT6
	if ctsIN6 >= 40 or ctsOUT6 >= 500 : cts[2,1]= 500+1
	cts[0,2]= ctsIN7
	cts[1,2]= ctsOUT7
	if ctsIN7 >= 40 or ctsOUT7 >= 500 : cts[2,2]= 500+1

	for ichip in range(len( chipv )) :

		chip=chipv[ichip]


		for ireg in range(len(regionname)) :

			if cts[ireg, ichip] >= mincounts[ireg] :


				for ii in range(len(rccd)) :

					if rccd[ii] == chip : 


						file1 = open( 'mfe_'+clu+'/reg/'+regionname[ireg]+'_'+ob+'_ccd'+str(chip)+'.reg', 'w')

						sstr='+polygon('


						for jj in range(rx.shape[1]) :
							if nn.isfinite(rx[ii,jj]) :
								if jj == 0 : sstr = sstr + str(rx[ii,jj]) + ',' + str(ry[ii,jj])
								else : sstr = sstr + ',' + str(rx[ii,jj]) + ',' + str(ry[ii,jj]) 
				 
						sstr = sstr + ')'

						file1.write( sstr + ' * ' + regionstr[ireg] )
						file1.write( exclpt )

						file1.close()






