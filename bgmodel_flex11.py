
# PURPOSE: set the parameters for the BG model

# v3: changes from v2: 
#     - got rid of g2
#     - fwhm of g1=.09 (instead of .18)
#     - min for gal1 and gal2 is zero, instead of 1e-4
# v4: * add setbgparams_nxt which is identical to setbgparams but treats the nxt model components
#     * fix bgab.nh at 0.0209
# v5 : add g2 and p1, which seem to be necessary for ACIS S blanksky
# v6 : * remove p1 and g2 component for ACIS-I observations
#      * start the fit with approximate values obtained from a previous fit
# v7 : * model names have ccd# in them
#      * "if hard" --> "if ccd=='i' or ccd=='6'"

# v10: * remove gaussians 6 and above
#      * set mins and maxes for all params

from sherpa.astro.ui import *

def setbgparams(ob, ccd, src) : 

	if src : exec 'cxb = cxb_'+ob+'_'+ccd
	exec 'cxbbg = cxbbg_'+ob+'_'+ccd
	if src : exec 'gal1 = gal1_'+ob+'_'+ccd
	exec 'gal1bg = gal1bg_'+ob+'_'+ccd
#	exec 'bgab = bgab_'+ob+'_'+ccd
	exec 'p2 = p2_'+ob+'_'+ccd
	exec 'p1 = p1_'+ob+'_'+ccd
	exec 'g1 = g1_'+ob+'_'+ccd
	exec 'g2 = g2_'+ob+'_'+ccd
	exec 'g3 = g3_'+ob+'_'+ccd
	exec 'g5 = g5_'+ob+'_'+ccd
	exec 'e1 = e1_'+ob+'_'+ccd

	if src : unlink(cxb.ampl)
	unlink(cxbbg.ampl)
	if src : unlink(gal1.norm)
	unlink(gal1bg.norm)
	unlink(p1.ampl)
	unlink(p2.ampl)
	unlink(g1.ampl)
	unlink(g2.ampl)
	unlink(g3.ampl)
	unlink(g5.ampl)
	unlink(e1.ampl)

	if src : 
		cxb.gamma=1.4
		set_par(cxb.ampl, val=1e-5, min=0, max=0.1)
		freeze(cxb)
		thaw(cxb.ampl)

	cxbbg.gamma=1.4
	set_par(cxbbg.ampl, val=1e-5, min=0, max=0.1)
	freeze(cxbbg)
	thaw(cxbbg.ampl)

	if src : 
		gal1.kt=0.2
		gal1.abundanc=1.0
		gal1.redshift=0.0
		set_par(gal1.norm, val=0.3, min=0, max=20)
		freeze(gal1)
		thaw(gal1.norm)

	gal1bg.kt=0.2
	gal1bg.abundanc=1.0
	gal1bg.redshift=0.0
	set_par(gal1bg.norm, val=0.3, min=0, max=20)
	freeze(gal1bg)
	thaw(gal1bg.norm)

#	freeze(bgab)
	set_par(bgab.nh, val=0.0209, min=0.004, max=0.2, frozen=True)

	# P1
	if ccd=='i' or ccd=='6' or ccd=='5' :
		set_par(p1.ampl, val=0., min=0., max=100., frozen=True)
		freeze(p1)
	else :
		set_par(p1.gamma, val=2.0, min=1.0, max=7.0)
		set_par(p1.ampl, val=0.1, min=0., max=100.)


	# P2
	set_par(p2.gamma, val=-0.4, min=-10, max=0)
	set_par(p2.ampl, val=0.02, min=0, max=10)
	thaw(p2)

	# G1	
	set_par(g1.pos, val=2.146, frozen=True)
	set_par(g1.fwhm, val=0.09 , frozen=True)
	set_par(g1.ampl, val=0.0612935, min=0., max=10, frozen=False)

	# G2
	g2.pos=1.787
	g2.fwhm=0.127

	# G3
	g3.pos=1.477
	g3.fwhm=0.00528215

	# G5
	g5.pos=2.562
	g5.fwhm=0.299358

	##########################################

	#freeze(g1)
	freeze(g2)
	freeze(g3)
	freeze(g5)

	if ccd=='i' or ccd=='6' : set_par(g2.ampl, min=0., val=0., max=10, frozen=True)
	else                    : set_par(g2.ampl, min=0., max=10, frozen=False)
	set_par(g3.ampl, min=0., max=50, frozen=False)
	set_par(g5.ampl, min=0., max=10, frozen=False)
#	set_par(g6.ampl, min=0., frozen=False)

	set_par(e1.coeff, val=-0.470387, min=-2, max=-0.1)
	freeze(e1)
	thaw(e1.ampl)
	#set_par( e1.ampl, val = 0.90132468447123559 * float(g1.ampl.val), min=0 )
	set_par( e1.ampl, min=0, max=10 )




	# ccd I
	if ccd == 'i' :
		p2.gamma  = -0.192965   
		p2.ampl   = 0.0100117   
		e1.ampl   = 0.0251768   
		g1.ampl   = 0.0464402   
		g3.ampl   = 0.232027    
		g5.ampl   = 0.00207296  
		if src : cxb.ampl  = 1e-5
		cxbbg.ampl= 1e-5
		if src : gal1.norm = 1e-5 
		gal1bg.norm = 1e-5 

	if ccd == '6' :
		p2.ampl   = 0.0105718
		e1.ampl   = 0.0244297
		g3.ampl   = 0.185136
		g5.ampl   = 0.00235066
		if src : cxb.ampl  = 8.4969e-06
		cxbbg.ampl= 8.4969e-06
		if src : gal1.norm = 3.22292e-05
		gal1bg.norm = 3.22292e-05

	# ccd 7
	if ccd == '7' : 
		p1.gamma   = 1.67118
		p1.ampl    = 0.00863298
		p2.gamma   =-1.57918
		p2.ampl    = 0.00113592
		e1.coeff   =-0.387611
		e1.ampl    = 0.0829656
		g1.ampl    = 0.0780609
		g2.ampl    = 0.0473212
		g3.ampl    = 0.156844
		g5.ampl    = 0
		if src : cxb.ampl   = 0
		cxbbg.ampl = 0
		if src : gal1.norm  = 4.21633e-05
		gal1bg.norm  = 4.21633e-05

	if ccd=='5' :
		p2.gamma  = -5.58987
		p2.ampl   = 2.00203e-06
		e1.ampl   = 0.0988879
		g1.ampl   = 0.0923749
		g2.ampl   = 0.0329375
		g3.ampl   = 0.0112692
		g5.ampl   = 0.00717219
		if src : cxb.ampl  = 0
		cxbbg.ampl= 0
		if src : gal1.norm = 2.39515e-05 
		gal1bg.norm = 2.39515e-05 


