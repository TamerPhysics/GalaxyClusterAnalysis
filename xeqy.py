
# PURPOSE: check that two floats are within 10^-8 of each other

def xeqy(x,y) :

	if abs(x-y) < 1e-8 : return True
	else               : return False

