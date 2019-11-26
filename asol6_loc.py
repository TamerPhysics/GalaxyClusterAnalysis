
# PURPOSE: Return the location of various data files that come with the 
# PURPOSE: observarion dataset.

# v3: * based on v1
#     * make sure that if both asolfile.fits and asolfile.fits.gz, are present
#       only the first is used.
#     * write absolute paths in asol_list.txt

# v6_loc: add MTL file types

import pdb

import glob
import os
import numpy as n

def asol(ob, ftype) :

	# Identifying the ASOL file(s).. (or more generally asol, flt1, pbk, msk, bpix or evt1 files)

	########################################## ASOL
	if ftype == 'asol' :

		asolfits = sorted(glob.glob('orig_*/'+ob+'/primary/*asol*.fits*'))

		if len(asolfits) != 0 : obpath = asolfits[0].partition('/primary/')[0]

		if len(asolfits) == 0 : 
			print 'NO ASOL FILE FOUND FOR OBSID ' + ob
			acao=[]
		elif len(asolfits) == 1 : acao = asolfits[0]
		else :  # --------- Multiple ASOL files OR .fits and .fits.gz both present

			#if asolfits[0]+'.gz' == asolfits[1] : acao = asolfits[0]

			newasolfits=[asolfits[0]]
			for i in range(1,len(asolfits)) :
				if asolfits[i-1]+'.gz' != asolfits[i] : newasolfits.append(asolfits[i])

			if len(newasolfits) == 1 : acao = newasolfits[0]
			else : 

				vernum = n.zeros(len(newasolfits), dtype=n.int)
				for i in range(len(newasolfits)) : 
					vernum[i] = int( newasolfits[i].split('_asol')[0][-1] )
				allequal=True
				for i in range(1,len(vernum)) :
					if vernum[i] != vernum[0] : allequal=False

				if allequal :
					if not os.path.exists('asol_list'+ob+'.txt') :
						asoltxtfile = open('asol_list'+ob+'.txt', 'w')
						for i in range(len(newasolfits)) :
							asoltxtfile.write( newasolfits[i]+'\n' )
						asoltxtfile.close()
					acao = '@'+'asol_list'+ob+'.txt'
				else : 
					maxwh = n.where(vernum == max(vernum))[0]
					if len(maxwh) == 1 : acao = newasolfits[maxwh[0]]
					else : 
						if not os.path.exists('asol_list'+ob+'.txt') :
							asoltxtfile = open('asol_list'+ob+'.txt', 'w')
							for wh in range(len(maxwh)) :
								asoltxtfile.write( newasolfits[wh]+'\n' )
							asoltxtfile.close()
						acao = '@'+'asol_list'+ob+'.txt'

		return acao





	########################################## anything other than ASOL:
	elif ftype == 'bpix' or ftype == 'evt1' or ftype == 'pbk' or ftype == 'flt1' or ftype == 'msk1' or ftype == 'mtl' :

		if ftype == 'bpix' : matches = sorted( glob.glob('orig_*/'+ob+'/primary/*'+ftype+'*') )
		else               : matches = sorted( glob.glob('orig_*/'+ob+'/secondary/*'+ftype+'*') )

		if len(matches) == 0 : 
			print 'NO '+ftype+' FILE FOUND FOR OBSID ' + ob
			acao=[]
		elif len(matches) == 1 : acao = matches[0]
		else :  # --------- Multiple ASOL files OR .fits and .fits.gz both present

			newmatches=[matches[0]]
			for i in range(1,len(matches)) :
				if matches[i-1]+'.gz' != matches[i] : newmatches.append(matches[i])

			if len(newmatches) == 1 : acao = newmatches[0]
			else :                    acao = newmatches[-1]


		return acao




	else : 
		print
		print '============ WRONG USAGE OF ASOL3.PY ========================='
		print
		return ''






