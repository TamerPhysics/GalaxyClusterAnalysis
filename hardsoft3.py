
# PURPOSE: return whether the aim point of the observation is on 
# PURPOSE: ACIS-I (h) or ACIS-S (s)

# v3: use the new location of repro3_evtOBSID.fits in mfe_OBSID/

import commands

def hardsoft(clu, obsids, hors) :

	hardobs=[]
	softobs=[]

	origclu = 'orig_'+clu

	for i in range(len(obsids)) :
		simz = float(commands.getoutput('dmkeypar mfe_'+clu+'/repro3_evt'+obsids[i]+'.fits sim_z echo+ '))
		if abs( simz - (-225.) ) <= abs( simz - (-190.) ): hardobs.append(obsids[i])
		else : softobs.append(obsids[i])



	if   hors == 'h' : return hardobs
	elif hors == 's' : return softobs
