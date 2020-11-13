import os
import sys
import itertools
import time
import pdb
import lot_size_MSSP

def de(MD):
	
	####################################################################
	###					Construct Scenarios
	####################################################################
	
	lot_size_MSSP.de(MD.parameters, MD.sets, MD.Uncertain, MD.Realizations, MD.make)

	### Generate Set of Outcomes

	
	
	return