
# PURPOSE: group data such that all adjacent zeros are group in 1 bin

from sherpa.astro.ui import *

import numpy

def grp(dataid, energystr) :

	#dataid = databgid[0]
	#bgid   = databgid[1]

	# group the data so that adjacent zero-count bins are grouped 
	# together.
	sust=False
	cc = get_data(dataid).counts
	grpv = numpy.zeros(len(cc), dtype=numpy.int)
	for ii in range(len(cc)) :
		if sust and cc[ii] == 0 : grpv[ii] = -1
		else : grpv[ii]=1
		if cc[ii] == 0 : sust=True
		else : sust=False
	if get_data(dataid).grouped : ungroup(dataid)
	notice_id(dataid)
	set_grouping(dataid, grpv)
	group(dataid)
	ignore_id(dataid, energystr)

	# group the BG the same was as the data
	#if get_data(bgid).grouped : ungroup(bgid)
	#notice_id(bgid)
	#set_grouping(bgid, grpv)
	#ignore_id(bgid, energystr)
	#group(bgid)
