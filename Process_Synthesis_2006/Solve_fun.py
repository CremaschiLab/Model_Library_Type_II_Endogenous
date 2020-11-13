import os
import sys
import itertools
import time
import pdb
import Network_Synthesis

def de(MD):
	
	####################################################################
	###					Construct Scenarios
	####################################################################
	
	Network_Synthesis.de(MD.parameters, MD.sets, MD.Uncertain, MD.Realizations, MD.make)

	### Generate Set of Outcomes

	
	
	return