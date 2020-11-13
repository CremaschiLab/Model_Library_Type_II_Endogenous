import os
import sys
import itertools
import time
import pdb
import portfolios

def de(MD):
	
	####################################################################
	###					Construct Scenarios
	####################################################################
	
	portfolios.de(MD.parameters, MD.sets, MD.Uncertain, MD.Realizations, MD.make)

	### Generate Set of Outcomes

	
	
	return