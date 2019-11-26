
# PURPOSE: return necessary counts for 20% error on Z, given (kT, bgf, Z)
# PURPOSE: and interpolating from pre-measured values in the above
# PURPOSE: parameter space

# v2 : changed the weights in the equation to find dnm and nmr
# v3 : changed the weights in the equation to find dnm and nmr again: weights
#      are 1/(normalized distance to point in the nnec array)**2
# v4 : weights = 1 / (dbg dkt dZ)
# v5 : * new bgf=0.8 simulations were made, and their values introduced in 
#        the array nnec
#      * 0.5 --> 0.8

import numpy as n

import pdb


def nnec(kt, bgf, z) :
#if True : 

	kt=float(kt)
	bgf=float(bgf)
	z=float(z)

	nnec = n.array([[[   1082.43004875,   2034.85328726],  \
                      [   3718.49570441,   2430.3070594 ],  \
                      [   8655.72660581,   3010.78026299],  \
                      [  20152.33124375,   6183.61073771]], \
                     [[   1198.12932844,   3727.01166979],  \
                      [   6096.52931254,   3095.94715995],  \
                      [  13634.49424594,   4560.82705006],  \
                      [  30580.58372612,   7825.17006602]], \
                     [[   4773.03517921,  10384.06328573],  \
                      [  18335.72096012,   8786.35135173],  \
                      [  62951.39816004,  15222.35502477],  \
                      [ 169951.87157461,  36317.91752711]], \
                     [[  22193.58610688,  44691.91347428],  \
                      [  98769.91785845,  46842.53369011],  \
                      [ 473686.47409494, 104049.102006  ],  \
                      [1884986.06062666, 274179.42216742]]] )


	ktv  =[1., 2., 5., 10.]
	zv   =[0.3,0.8]
	bgfv =[0.0, 0.1, 0.5, 0.8]

	bgf0= max(0.0, bgf)
	bgf0= min(0.8, bgf0)
	kt0 = max(1.0, kt)
	kt0 = min(10.0, kt0)
	z0  = max(0.3, z)
	z0  = min(0.8, z0)

	# If the input (bgf, kt, z) are on the grid of available (bgf, kt, z)'s,
	# then no need for interpolation, just use the value in the nnec array
	for ibg in range(len(bgfv)) :
		for it in range(len(ktv)) :
			for iz in range(len(zv)) :
				if bgf==bgfv[ibg]  and kt==ktv[it]  and z0==zv[iz] : return nnec[ibg,it,iz]
				if bgf0==bgfv[ibg] and kt0==ktv[it] and z0==zv[iz] :
					nn = nnec[ibg,it,iz]
					if kt  > 10. : nn = nn * (kt/10.)**1.2
					if bgf > 0.8 : nn = nn * (1-0.8)/(1-bgf)
					return nn



	# -------------------------- Interpolation --------------------------------
	for ibg in range(len(bgfv)-1) :
		if bgf0 >= bgfv[ibg] and bgf0 < bgfv[ibg+1] : ibg1 = ibg
	if bgf0 == 0.8 : ibg1 = len(bgfv)-2 # this expression is more general than the one in nnecz4

	for ikt in range(len(ktv)-1) :
		if kt0 >= ktv[ikt] and kt0 < ktv[ikt+1] : it1 = ikt
	if kt0 == 10.0 : it1 = len(ktv)-2

	nmr=0. # numerator   = SUM (weights * N)
	dnm=0. # denominator = SUM (weights)
	dbg = bgfv[ibg1+1] - bgfv[ibg1]
	dt  = ktv[it1+1] - ktv[it1]
	dz  = 0.8 - 0.3
	ibgv = [ibg1, ibg1+1]
	itv  = [it1, it1+1]
	for ibg in ibgv : 
		for it in itv :
			for iz in [0, 1] :

				dnm0 = 1.

				if bgf0 not in bgfv : dnm0 = dnm0 / abs(bgf0-bgfv[ibg])
				elif bgf0 != bgfv[ibg] : dnm0 = 0.

				if kt0 not in ktv : dnm0 = dnm0 / abs(kt0-ktv[it])
				elif kt0 != ktv[it] : dnm0 = 0.

				if z0 not in zv : dnm0 = dnm0 / abs(z0-zv[iz])
				elif z0 != zv[iz] : dnm0 = 0.

				nmr = nmr + dnm0 * nnec[ibg, it, iz]
				dnm = dnm + dnm0

	if bgf <= 0.8 and kt <= 10. : return nmr/dnm
	else :
		nn = nmr/dnm
		if kt  > 10. : nn = nn * (kt/10.)**1.2
		if bgf > 0.8 : nn = nn * (1-0.8)/(1-bgf)
		return nn



























