import os
import sys
import itertools
import time
import pdb
import DSR

def de(MD):
	
	####################################################################
	###					Construct Scenarios
	####################################################################
	
	DSR.de(MD.parameters, MD.sets, MD.Uncertain, MD.Realizations, MD.make)

	### Generate Set of Outcomes

	
	
	return