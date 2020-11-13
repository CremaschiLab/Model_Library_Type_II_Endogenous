import os
import sys
import itertools
import time
import pdb
import open_pit_ming

def de(MD):
	
	####################################################################
	###					Construct Scenarios
	####################################################################
	
	open_pit_ming.de(MD.parameters, MD.sets, MD.Uncertain, MD.Realizations, MD.make)

	### Generate Set of Outcomes

	
	
	return