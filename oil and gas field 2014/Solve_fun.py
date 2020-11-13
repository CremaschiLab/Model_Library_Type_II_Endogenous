import os
import sys
import itertools
import time
import pdb
import oil_field

def de(MD):
	
	####################################################################
	###					Construct Scenarios
	####################################################################
	
	oil_field.de(MD.parameters, MD.sets, MD.Uncertain, MD.Realizations, MD.make)

	### Generate Set of Outcomes

	
	
	return