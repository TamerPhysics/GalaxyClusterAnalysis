
# PURPOSE: fit inner and outer regions of a cluster

# v3 for entropy/mfe_*
#  modification: no more the '_v2' at the end of file names.

# v3sherpa42: modified save_data --> save_table, because in 
#             sherpa 4.2, save_data produced fits files instead
#             of text files

# v5: * based on v3sherpa42
#     * changed the guess statement (suspected to have error) to a manual
#       way to estimate normalization, taken from ktofr8.py
#     * make sure when you call out0 that it's defined. in clusters where 
#       there is no out region defined, the variable out0 is never
#       assigned.
#     * fit3_*.ps --> fit5_*.ps

# v6: * same method for both ACIS-S and ACIS-I. for easiness, used same code
#       from fit5.py, keeping only the ACIS-S part, and adding the new chip#
#       "i" to 5,6,7
#     * include spectrum if WARF(ob,ccd,reg) file exists (as opposed to checking
#       existence of RMF file)

# OUTPUT: LOCCLU/bestfitnh.txt
# OUTPUT: cluname_TZ.txt and TZ.txt (results appended to file)
# OUTPUT: cluname_normTZ.txt (results appended to file)
# OUTPUT: chi0TZ.txt (results appended to file)
# OUTPUT: errTZ.txt (results appended to file)
# OUTPUT: LOCCLU/spec/srcin.txt and LOCCLU/spec/srcout.txt and LOCCLU/spec/srctot.txt

import glob
import os
import commands

import sherpa
from sherpa.astro.ui import *
from pychips.all import *

def clufit(clu, z, nh) :

	locclu = 'mfe_' + clu
	emin=0.3
	clean()


	in2soft = glob.glob(locclu+'/spec/sp*_in2_ccd[i567]_grp.pi')
	out2soft = glob.glob(locclu+'/spec/sp*_out2_ccd[i567]_grp.pi')


	################################################################################################################
	#                     SOFT CHIPS 
	################################################################################################################

	# Analysis for observations on CCD's 5, 6 and 7 -- region=in
	softinobs  = []
	nobs=1
	for i in range(len(in2soft)) :

		obs = in2soft[i].rpartition('/')[2].partition('_')[0][2:]
		ccds= in2soft[i].rpartition('/')[2].partition('_ccd')[2][0]

		scnt=commands.getoutput('dmlist "'+in2soft[i]+'[cols counts]" data,clean')
		cntv=scnt.split('\n')
		cntv.remove(cntv[0])
		for ii in range(len(cntv)) : 
			cntv[ii]=float(cntv[ii])

		if sum(cntv) > 100 and os.path.exists(locclu+'/spec/sp'+obs+'_in2_ccd'+ccds+'.warf'):

			softinobs.append(nobs)

			print 'LOADING .... '+in2soft[i]+' ------------'
			load_pha(nobs, in2soft[i] )
			load_rmf(nobs, locclu+'/spec/sp'+obs+'_ccd'+ccds+'_center.wrmf')
			load_arf(nobs, locclu+'/spec/sp'+obs+'_in2_ccd'+ccds+'.warf')
			load_bkg(nobs, locclu+'/spec/sp'+obs+'_in2_ccd'+ccds+'_bgbg2.pi')

			if len(softinobs) == 1 : 
				in0 = nobs
				exec('set_source(nobs, xsphabs.abs1 * xsapec.in' +str(nobs)+')')
				exec('in'+str(nobs)+'.kt = 5')
				exec('in'+str(nobs)+'.abundanc = .3')
				exec('thaw(in'+str(nobs)+'.abundanc)')
				exec('in'+str(nobs)+'.redshift = '+str(z) )
			else :
				exec('set_source(nobs, abs1 * xsapec.in' +str(nobs)+')')
				exec('in'+str(nobs)+'.kt = '+'in'+str(in0)+'.kt')
				exec('in'+str(nobs)+'.abundanc = '+'in'+str(in0)+'.abundanc')
				exec('in'+str(nobs)+'.redshift = '+'in'+str(in0)+'.redshift')
			exec 'emnorm0=abs(float(in'+str(nobs)+'.norm.val)*calc_data_sum(id=nobs,lo=0.3,hi=0.7)/calc_model_sum(id=nobs,lo=0.3,hi=0.7))'
			exec 'set_par(in'+str(nobs)+'.norm, val=float(emnorm0), min=float(emnorm0/1e3), max=float(emnorm0*1e3) )'
			subtract(nobs)

			nobs=nobs+1

	# Analysis for observations on CCD's 5, 6 and 7 -- region=out
	softoutobs  = []
	for i in range(len(out2soft)) :

		obs = out2soft[i].rpartition('/')[2].partition('_')[0][2:]
		ccds= out2soft[i].rpartition('/')[2].partition('_ccd')[2][0]

		scnt=commands.getoutput('dmlist "'+out2soft[i]+'[cols counts]" data,clean')
		cntv=scnt.split('\n')
		cntv.remove(cntv[0])
		for ii in range(len(cntv)) : 
			cntv[ii]=float(cntv[ii])

		if sum(cntv) > 100 and os.path.exists(locclu+'/spec/sp'+obs+'_out2_ccd'+ccds+'.warf'):

			softoutobs.append(nobs)

			print 'LOADING .... '+out2soft[i]+' ------------'
			load_pha(nobs, out2soft[i] )
			load_rmf(nobs, locclu+'/spec/sp'+obs+'_ccd'+ccds+'_center.wrmf')
			load_arf(nobs, locclu+'/spec/sp'+obs+'_out2_ccd'+ccds+'.warf')
			load_bkg(nobs, locclu+'/spec/sp'+obs+'_out2_ccd'+ccds+'_bgbg2.pi')

			if len(softoutobs) == 1 : 
				out0 = nobs
				exec('set_source(nobs, xsphabs.abs1 * xsapec.out' +str(nobs)+')')
				exec('out'+str(nobs)+'.kt = 5')
				exec('out'+str(nobs)+'.abundanc = .3')
				exec('thaw(out'+str(nobs)+'.abundanc)')
				exec('out'+str(nobs)+'.redshift = '+str(z) )
			else :
				exec('set_source(nobs, abs1 * xsapec.out' +str(nobs)+')')
				exec('out'+str(nobs)+'.kt = '+'out'+str(out0)+'.kt')
				exec('out'+str(nobs)+'.abundanc = '+'out'+str(out0)+'.abundanc')
				exec('out'+str(nobs)+'.redshift = '+'out'+str(out0)+'.redshift')
			exec 'emnorm0=abs(float(out'+str(nobs)+'.norm.val)*calc_data_sum(id=nobs,lo=0.3,hi=0.7)/calc_model_sum(id=nobs,lo=0.3,hi=0.7))'
			exec 'set_par(out'+str(nobs)+'.norm, val=float(emnorm0), min=float(emnorm0/1e3), max=float(emnorm0*1e3) )'
			subtract(nobs)

			nobs=nobs+1



	################################################################################################################


	if nobs > 1 :
		abs1.nh = nh
		abs1.nh.min= nh/5.
		abs1.nh.max= nh*5.
#		freeze(abs1.nh)

		ignore('0:'+str(emin)+', 7:')

		fit()

		absfile = open(locclu+'/bestfitnh.txt', 'w')
		absfile.write(str(abs1.nh.val)+'\n')
		absfile.close()

		if get_fit_results().rstat > 3 :
			print ' ------------------ BAD FIT.. TRY MANUAL FIT --------------------'
			cont=''
			while cont != 'xxx' :
				cont = raw_input ('Enter command>> ')
				if cont =='xxx' : break
				else : 
					try : exec cont
					except : pass

		inobs  = softinobs 
		outobs = softoutobs 

		# Write the best fit Tin, Zin, Tout Zout and spectrum normalization
		print
		print 'Writing results to file.. '
		print
		nresfile = open('cluname_TZ.txt', 'a')
		resfile  = open('TZ.txt', 'a')
		normfile = open('cluname_normTZ.txt', 'a')
		if 'out0' in locals() :
			exec 'nresfile.write ( clu + " " + str(in'+str(in0)+'.kt.val) + " " + str(in'+str(in0)+'.abundanc.val) + " " + str(out'+str(out0)+'.kt.val) + " " + str(out'+str(out0)+'.abundanc.val) + "\\n" )'
			exec 'resfile.write ( str(in'+str(in0)+'.kt.val) + " " + str(in'+str(in0)+'.abundanc.val) + " " + str(out'+str(out0)+'.kt.val) + " " + str(out'+str(out0)+'.abundanc.val) + "\\n" )'
			exec 'normfile.write ( clu + " " + str(in'+str(in0)+'.norm.val) + " " + str(out'+str(out0)+'.norm.val) + "\\n" )'
		else : 
			exec 'nresfile.write ( clu + " " + str(in'+str(in0)+'.kt.val) + " " + str(in'+str(in0)+'.abundanc.val) + " -1 -1\\n" )'
			exec 'resfile.write ( str(in'+str(in0)+'.kt.val) + " " + str(in'+str(in0)+'.abundanc.val) + " -1 -1\\n" )'
			exec 'normfile.write ( clu + " " + str(in'+str(in0)+'.norm.val) + " -1\\n" )'

		nresfile.close()
		resfile.close()
		normfile.close()

		# Write the reduced chi squared:
		chifile  = open('chi0TZ.txt', 'a')
		chifile.write(str(get_fit_results().rstat)+'\n')
		chifile.close()


		# Uncertainties in the parameters
		print 'Calculating uncertainties ...'
		print
		instr=''
		for i in range(len(inobs)) : instr = instr + str(inobs[i]) + ','
		outstr=''
		for i in range(len(outobs)) : outstr = outstr + str(outobs[i]) + ','

		class backupproj(object) :
			parmins=['NaN','NaN']
			parmaxes=['NaN','NaN']

		try : 
			exec 'proj('+instr+' in'+str(in0)+'.kt, in'+str(in0)+'.abundanc)'
			resin2 = get_proj_results()
		except : 
			resin2=backupproj()

		try : 
			exec 'proj('+outstr+' out'+str(out0)+'.kt, out'+str(out0)+'.abundanc)'
			resout2 = get_proj_results()
		except : 
			resout2=backupproj()


		uncfile = open('errTZ.txt', 'a')
		# T_in_min T_in_max Z_in_min Z_in_max T_out_min T_out_max Z_out_min Z_out_max 
		uncfile.write( str(resin2.parmins[0]) + ' ' + str(resin2.parmaxes[0]) + ' ' + str(resin2.parmins[1]) + ' ' + str(resin2.parmaxes[1]) + ' ' +   str(resout2.parmins[0]) + ' ' + str(resout2.parmaxes[0]) + ' ' + str(resout2.parmins[1]) + ' ' + str(resout2.parmaxes[1]) + '\n' )
		uncfile.close()



		# Write the spectra to txt file
		print 'Writing source spectrum to file ...'
		print
		load_arrays('srcin', get_source_plot(in0).xlo, get_source_plot(in0).y ) 
		save_table('srcin', locclu+'/spec/srcin.txt[opt kernel=text/raw]')
		if 'out0' in locals() :
			load_arrays('srcout', get_source_plot(out0).xlo, get_source_plot(out0).y ) 
			save_table('srcout', locclu+'/spec/srcout.txt[opt kernel=text/raw]')
		else : os.system('cp '+locclu+'/spec/srcin.txt '+locclu+'/spec/srcout.txt')

		# Save a .PS version of the plots
		add_window(11,8.5,'inches')
		allobs = inobs + outobs
		allobs.sort()
		for obs in allobs : 
			plot_fit_resid(obs)
			if os.path.exists('fitplots/fit5_'+clu+'_'+str(obs)+'.ps') : os.remove('fitplots/fit5_'+clu+'_'+str(obs)+'.ps')
			print_window ('fitplots/fit5_'+clu+'_'+str(obs), ['format', 'ps', 'orientation', 'landscape'])

		# Write the spectrum of entire cluster (in+out regions) to txt file
		if 'out0' in locals() : 

			try : exec 'in'+str(in0)+'.kt = out'+str(out0)+'.kt'
			except NameError : pass

			try : exec 'in'+str(in0)+'.abundanc = out'+str(out0)+'.abundanc'
			except NameError : pass

			exec 'fit'+str(tuple(allobs))

			load_arrays('srctot', get_source_plot(in0).xlo, get_source_plot(in0).y ) 
			save_table('srctot', locclu+'/spec/srctot.txt[opt kernel=text/raw]')

		else : os.system('cp '+locclu+'/spec/srcin.txt '+locclu+'/spec/srctot.txt')



		delete_window()




