
# PURPOSE: Create spectra for inner and outer regions for a given cluster
# PURPOSE: and OBSID.

# v3 for entropy/mfe_*
# v5: * creates the .pi files for the next blank sky period BG files, when available.
#     * uses bkg4ent2 instead of bkg4ent
# v6: * use asol3.py to find asol, flt1, pbk, msk, bpix or evt1 files
#     * arf4ent --> arf4ent2
# v7: * treat acis-i like acis-s chips
#     * make RMF file if BG evt file exists
#     * make sp/warf/bgsp if both region file (from make_chips5) and 
#       BG evt file exist.
# v8: * stop making _nxt spectra for BG from blank-sky of the next period
#     * use make_chips_reg6.pro, as it has a bug fix from previous version
# v9: * 2017-03-31: new version of asol: asol6_loc
#     * remove pbk keyword

# OUTPUT: reg/chipsregOB.fits
# OUTPUT: reg/in.reg reg/out.reg reg/bg.reg
# OUTPUT: through make_chips_reg4.pro: reg/in2_OBSID.reg reg/out2_OBSID.reg reg/bg2_OBSID.reg
# OUTPUT: most files in spec/ 

import glob
import os
import commands
import pdb

import asol6_loc
import centradec
import rmf4ent5
import chipreg

def sp (clu, ob, r1mpc) :

	locclu = 'mfe_'+clu
	origclu = 'orig_'+clu

	os.system('punlearn ardlib')

	acao = asol6_loc.asol(ob, 'asol')

	if not os.path.exists(locclu+'/spec') : os.mkdir(locclu+'/spec')

	(rac, decc) = centradec.getrd(locclu)

	# simple in, out and bg regions (w/o point source exclusion)
	reginfile = open( locclu+'/reg/in.reg', 'w')
	reginfile.write('circle('+rac+','+decc+','+str(r1mpc/10)+'")\n')
	reginfile.close()

	regoutfile = open( locclu+'/reg/out.reg', 'w')
	regoutfile.write('annulus('+rac+','+decc+','+str(r1mpc/10)+'",'+str(r1mpc)+'")\n')
	regoutfile.close()

	regbgfile = open( locclu+'/reg/bg.reg', 'w')
	regbgfile.write('circle('+rac+','+decc+','+str(r1mpc)+'")\n')
	regbgfile.close()

	# make chip regions using skyfov
	os.system('punlearn skyfov')
	os.system('skyfov '+locclu+'/clean'+ob+'.fits '+locclu+'/reg/chipsreg'+ob+'.fits aspect='+acao+' mskfile='+asol6_loc.asol(ob,'msk1'))

	pbk = asol6_loc.asol(ob,'pbk')

	regv = ['in2','out2']

	ctsin0  = commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/in.reg)][ccd_id=0]" counts')
	ctsin1  = commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/in.reg)][ccd_id=1]" counts')
	ctsin2  = commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/in.reg)][ccd_id=2]" counts')
	ctsin3  = commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/in.reg)][ccd_id=3]" counts')
	ctsin5  = commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/in.reg)][ccd_id=5]" counts')
	ctsin6  = commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/in.reg)][ccd_id=6]" counts')
	ctsin7  = commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/in.reg)][ccd_id=7]" counts')

	ctsout0 = commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/out.reg)][ccd_id=0]" counts')
	ctsout1 = commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/out.reg)][ccd_id=1]" counts')
	ctsout2 = commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/out.reg)][ccd_id=2]" counts')
	ctsout3 = commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/out.reg)][ccd_id=3]" counts')
	ctsout5 = commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/out.reg)][ccd_id=5]" counts')
	ctsout6 = commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/out.reg)][ccd_id=6]" counts')
	ctsout7 = commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/out.reg)][ccd_id=7]" counts')

	# Create regions file combining field-of-view regions and point source exclusion regions
	chipreg.make(clu, ob, r1mpc, float(ctsin0), float(ctsin1), float(ctsin2), float(ctsin3), float(ctsin5), float(ctsin6), float(ctsin7), float(ctsout0), float(ctsout1), float(ctsout2), float(ctsout3), float(ctsout5), float(ctsout6), float(ctsout7))


	ccdv = ['i','5','6','7']
	ccdlist=['0,1,2,3','5','6','7']

	for iccd in range(len(ccdv)) :
		# RMF -----------------------
		#if not os.path.exists(locclu+'/spec/sp'+ob+'_ccd'+ccdv[iccd]+'_center.wrmf') and 
		#if float(commands.getoutput('dmlist "'+locclu+'/clean'+ob+'.fits[ccd_id='+ccdlist[iccd]+']" counts')) > 0 
		if os.path.exists(locclu+'/bg/bg'+ob+'_ccd'+ccdv[iccd]+'r.fits') : rmf4ent5.rmf(locclu, ob, ccdv[iccd])

	for ireg in range(len(regv)) :
		for iccd in range(len(ccdv)) :

			# the region file, below, is created if there are enough counts in the region
			# of interest.
			if os.path.exists(locclu+'/reg/'+regv[ireg]+'_'+ob+'_ccd'+ccdv[iccd]+'.reg') and os.path.exists(locclu+'/bg/bg'+ob+'_ccd'+ccdv[iccd]+'r.fits') :

				# DMEXTRACT -----------
				os.system('punlearn dmextract')
				os.system('pset dmextract infile="'+locclu+'/clean'+ob+'.fits[sky=region('+locclu+'/reg/'+regv[ireg]+'_'+ob+'_ccd'+ccdv[iccd]+'.reg)][bin PI]"')
				os.system('pset dmextract outfile='+locclu+'/spec/sp'+ob+'_'+regv[ireg]+'_ccd'+ccdv[iccd]+'.pi')
				os.system('pset dmextract opt=pha1')
				os.system('pset dmextract wmap="[energy=300:2000][bin det=8]"')
				os.system('dmextract verbose=1 mode=h')


				# DMGROUP -----------
				os.system('punlearn dmgroup')
				os.system('pset dmgroup infile='+locclu+'/spec/sp'+ob+'_'+regv[ireg]+'_ccd'+ccdv[iccd]+'.pi')
				os.system('pset dmgroup outfile='+locclu+'/spec/sp'+ob+'_'+regv[ireg]+'_ccd'+ccdv[iccd]+'_grp.pi')
				os.system('pset dmgroup grouptype=ADAPTIVE grouptypeval=30')
				os.system('pset dmgroup xcolumn="channel" ycolumn="counts"')
				os.system('pset dmgroup binspec=""')
				os.system('dmgroup verbose=1 mode=h')

				# MKWARF --------------
				os.system('punlearn mkwarf')
				os.system('pset mkwarf infile="'+locclu+'/spec/sp'+ob+'_'+regv[ireg]+'_ccd'+ccdv[iccd]+'.pi[WMAP]"')
				os.system('pset mkwarf outfile='+locclu+'/spec/sp'+ob+'_'+regv[ireg]+'_ccd'+ccdv[iccd]+'.warf')
				os.system('pset mkwarf weightfile='+locclu+'/spec/sp'+ob+'_'+regv[ireg]+'_ccd'+ccdv[iccd]+'.wfef')
				#os.system('pset mkwarf pbkfile='+pbk)
				os.system('pset mkwarf egridspec="0.3:11.0:0.01"')
				os.system('pset mkwarf asolfile='+acao)
				os.system('mkwarf verbose=1 mode=h')


				# DMEXTRACT BG -----------
				os.system('punlearn dmextract')
				os.system('pset dmextract infile="'+locclu+'/bg/bg'+ob+'_ccd'+ccdv[iccd]+'r.fits[sky=region('+locclu+'/reg/'+regv[ireg]+'_'+ob+'_ccd'+ccdv[iccd]+'.reg)][bin PI]"')
				os.system('pset dmextract outfile='+locclu+'/spec/sp'+ob+'_'+regv[ireg]+'_ccd'+ccdv[iccd]+'_bgbg2.pi')
				os.system('pset dmextract opt=pha1')
				os.system('dmextract verbose=1 mode=h')




