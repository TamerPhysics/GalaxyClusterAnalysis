
# This Python script calls sequentially all the functions needed to 
# analyze Chandra X-ray Observatory data of galaxy clusters.
# The end result of this analysis is to produce, for each galaxy 
# cluster:
# - a measurement of temperature for as many radial bins as data
#   allow
# - a measurement of metallicity for as many radial bins as data
#   allow
# This requires join model fitting of all observations in a given
# radial bin, for a given cluster. I also attempt to be efficient
# with the photon counts in the data, using the least possible
# counts for an estimated minimum of 10% error on temperature
# measurements, and 20% error on metallicity measurements.

import os
import glob
import commands
import numpy as nn
import pdb
import multiprocessing as mpr

import repevt5
import mrg6
import pt2 # followed by manual step
import lc5
import bg12
import sp9 
import fit6
import exp7

import preradial5
import globmod
import globrmf
import checkbgfit5
import ktofr15loc
import zofr11


# Hubble constant for calculating R500, angular size only!
h0=72.

# Load list of clusters to be analyzed
clufile = open ( 'clulist.txt' , 'r' )

# alphabetically ordered cluster names, needed for nH file reading
clufileab = open ( 'clulist_alphab.txt' , 'r' )

# file containing values of Hydrogen column density measured in the direction
# of each cluster. 
# nH values obtained from https://heasarc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl
nhfile = open('nh_colden.txt', 'r')

# This array will contain the list of clusters to analyze
clulist = []

# this dictionary will contain the observation IDs of each cluster
obsids={}

# this dictionary variable will contain values of nH from the LAB survey
nhdlab = {}

# Read first lines of the files
cluline = clufile.readline()
clulineab = clufileab.readline()
nhline = nhfile.readline()

while cluline != '' and cluline != '\n' :

	cluname = cluline[0:-1]

	# populate the clulist variable
	clulist.append(cluname)

	# populate the nhdlab variable
	clunameab = clulineab[0:-1]
	nhdlab[clunameab] = float(nhline.split()[8])/100.

	# Create local directory link to data directory
	if not os.path.exists('orig_'+cluname) :
		os.symlink('/home/user/chandradata/'+cluname, 'orig_'+cluname)

	# Find all observations of this cluster
	thiscluobsids = glob.glob('orig_'+cluname+'/[1-9]*')
	for j in range(len(thiscluobsids)) :
		thiscluobsids[j] = thiscluobsids[j][6+len(cluname):len(thiscluobsids[j])]

	# populate obsids
	obsids[cluname] = thiscluobsids
	obsids[cluname].sort()

	# read next line
	cluline   =   clufile.readline()
	clulineab = clufileab.readline()
	nhline = nhfile.readline()	


# close all files used above
clufile.close()
clufileab.close()
nhfile.close()

# sort clulist array
clulist.sort()

# Make local folder for each cluster, which will contain
# all the analysis products
for clu in list(clulist) :
	if not os.path.exists('mfe_'+clu) : os.mkdir('mfe_'+clu)

# Read redshift value of each cluster and nH from Dickey & Lockman survey
znhfile = open('hiflugcs_mod.txt','r') #!
zd = {}
nhd = {}
znhline = znhfile.readline()
while znhline != '' and znhline != '\n' :
	zd[  znhline.split('\t')[1] ] = float(znhline.split('\t')[4]) 
	nhd[ znhline.split('\t')[1] ] = float(znhline.split('\t')[5])/100. 
	# in units 10^22 atoms/cm^2
	znhline = znhfile.readline()
znhfile.close()

# Read angular diameter distance, and compute r1mpcd = the angle
# on the sky corresponding to 1Mpc
dad={} # in Mpc
r1mpcd={} # in arcseconds
with open('da_mpc.txt', 'r') as dafile :
	dalines = dafile.readlines()

for line in dalines : 
	dad[line.split()[0]] = float(line.split()[1])
	r1mpcd[line.split()[0]] = 1/dad[line.split()[0]]*180/nn.pi*3600

# Can choose a subset of clulist to analyze here:	
#clulist=['a2244']

# this directory contains output from wavdetect tool, which is not needed in this analysis
if not os.path.exists('wavzebala') : os.mkdir('wavzebala') 

# This directory will contain many of the output plots further down in the analysis
if not os.path.exists('fitplots') : os.mkdir('fitplots')

longestexp = []
for i in range(len(clulist)) :

	clu = clulist[i]

	# Reprocess all observations with the installed CALDB. This applies
	# the latest Chandra calibration data to the data.
	for j in range(len(obsids[clulist[i]])) :

		ob = obsids[clulist[i]][j]
		repevt5.repevt(clulist[i], obsids[clulist[i]][j])


	# Finding the OBSID with the longest exposure for each cluster
	if len(obsids[clulist[i]]) == 1 : 
		longestexp.append(obsids[clulist[i]][0])
	else :
		exptime = nn.zeros(len( obsids[clulist[i]] ))
		for j in range(len( obsids[clulist[i]] )) :
			exptime[j] = float ( commands.getoutput('dmkeypar mfe_'+clulist[i]+'/repro3_evt'+obsids[clulist[i]][j]+'.fits livetime echo+') )
		longestexp.append(obsids[clulist[i]][nn.argmax(exptime)])


	# Making merged evt lists + 256x256 images of size 1Mpc x 1Mpc to be used for 
	# point source detection
	mrg6.mrg( clulist[i], obsids[clulist[i]], longestexp[i] ) #<>

	# Detect point sources
	pt2.pt(clulist[i]) #<>

# As per Chandra team recommendation:
# remove datasets taken in period A, ie before 1999-09-16
# remove datasets from period B, as well, ie before 2000-09-17
for clu in list(clulist) :

	for ob in list(obsids[clu]) :

		f = glob.glob('orig_'+clu+'/secondary/*evt1*')
		if len(f)==0 : f = glob.glob('mfe_'+clu+'/repro3_evt'+ob+'.fits')

		# Read date of observation		
		sd = commands.getoutput ( 'dmkeypar '+f[0]+' date-obs echo+' )
		intdat = int(sd[0:4]+sd[5:7]+sd[8:10])

		# Periods A and B
		if intdat < 20000128: 
			obsids[clu].remove(ob)

	# remove clusters which have no obsids due to above filtering
	if len(obsids[clu])==0 : clulist.remove(clu)

# make file containing list of OBSIDS
for clu in clulist :
	obfile = open('mfe_'+clu+'/obsids.txt', 'w')
	for ob in obsids[clu] : obfile.write(ob+'\n')
	obfile.close()

##########################################################################
######## These are the only required manual parts in the script ##########
for i in range(len(clulist)) :

	# Check the point detection result, and fix it for any anomalies
	print
	print 'Save the corrected region as mfe_'+clulist[i]+'/reg/pt0mfe_wcs.reg'
	print '----------------------------------------------------------------'
	print
	if not os.path.exists('mfe_'+clulist[i]+'/reg/pt0mfe_wcs.reg') : os.system('ds9 mfe_'+clulist[i]+'/mrgimg_en.fits -regions load mfe_'+clulist[i]+'/reg/pt0mfe.reg -regions format ciao -regions system wcs') #<>

# Filter data for any periods of high count rate, which is usually due
# to solar activity and other local signals. Check the automatic filtering
# and fix it manually at the prompt
for i in range(len(clulist)) : 
	lc5.lc(clulist[i], obsids[clulist[i]]) #<>

print
print '============ Create In-field region for BG extraction ==============='
print '============        mfe_CLU/bginfield_all.reg         ==============='
print
print
##########################################################################
##########################################################################

# Longest exposure for the lightcurve-cleaned data
cleanlongest={}
for clu in clulist :
	cleantime = nn.zeros(len( obsids[clu] ))
	for j in range(len( obsids[clu] )) :
		cleantime[j] = float ( commands.getoutput('dmkeypar mfe_'+clu+'/clean'+obsids[clu][j]+'.fits livetime echo+') )
	cleanlongest[clu] = obsids[clu][nn.argmax(cleantime)] 

# RMF files take a long time to create. In order to save time, keep
# track of all available RMF files, in order to not run the RMF 
# creation commands multiple times.
globrmf.name = []
globrmf.scat = []
globrmf.gain = []
rmfs = glob.glob('mfe_*/spec/bgrmf*.fits')
for rmf in rmfs :
	if not os.path.islink(rmf) :
		globrmf.name.append(rmf)
		globrmf.scat.append(commands.getoutput('dmkeypar '+rmf+' scatfile echo+'))
		globrmf.gain.append(commands.getoutput('dmkeypar '+rmf+' gainfile echo+'))

for i in range(len(clulist)) :

	# Make BG files from blank sky observations
	if not os.path.exists('mfe_'+clulist[i]+'/bg') : os.mkdir('mfe_'+clulist[i]+'/bg')
	for j in range(len(obsids[clulist[i]])) : 
		bg12.bg(clulist[i], obsids[clulist[i]][j]) #<>

	# Make spectra for inner and outer regions, demarkated by radius = 100kpc
	for j in range(len(obsids[clulist[i]])) : 
		sp9.sp(clulist[i], obsids[clulist[i]][j], r1mpcd[clulist[i]]) #<>

	# Fit spectra for inner and outer regions, demarkated by radius = 100kpc
	fit6.clufit(clulist[i], zd[clulist[i]], nhdlab[clulist[i]]) #<>
	os.system('sed -i s/None/NaN/g errTZ.txt') #<>

	# Make exposure maps for the inner and outer regions of the cluster
	os.system('rm mfe_'+clulist[i]+'/evt*_c.fits') 
	exp7.exp(clulist[i], obsids[clulist[i]], 'in') #<>
	exp7.exp(clulist[i], obsids[clulist[i]], 'out') #<>

	# Make radial intensity profiles
	if not os.path.exists('mfe_'+clulist[i]+'/ctofr') :
		preradial5.preradial(clulist[i], obsids[clulist[i]], r1mpcd[clulist[i]]) #<>

	globmod.globstr1 = ''
	globmod.globstr2 = ''
	globmod.globstr3 = ''

	# Fit the background with the possibility of manually checking the fit
	if os.path.exists('mfe_'+clulist[i]+'/bginfield_all.reg') and \
          len(glob.glob('mfe_'+clulist[i]+'/ktofrsp/bgparams8_*.txt'))==0 : 
		checkbgfit5.ktofr(clulist[i], obsids[clulist[i]], zd[clulist[i]], nhd[clulist[i]], False) #<>

	# Make radial temperature profiles
	ktofr15loc.ktofr(clulist[i], obsids[clulist[i]], r1mpcd[clulist[i]], zd[clulist[i]], nhdlab[clulist[i]]) #<>

	# Make radial metallicity profiles
	if os.path.exists('mfe_'+clulist[i]+'/ktofr14_hien.txt') : 
		zofr11.zofr(clulist[i], obsids[clulist[i]], r1mpcd[clulist[i]], zd[clulist[i]], nhd[clulist[i]], nhdlab[clulist[i]]) #<>







