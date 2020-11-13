import os
import sys
import itertools
import time
import pdb
import gas_model

def de(MD):
	
	####################################################################
	###					Construct Scenarios
	####################################################################
	# List_of_Possible_Realizations = [range(len(MD.Realizations['outcomes'][r])) for r in MD.Realizations['outcomes']]

	# print(List_of_Possible_Realizations)
	
	gas_model.de(MD.parameters, MD.sets, MD.Uncertain, MD.Realizations)

	### Generate Set of Outcomes
	# Outcomes = itertools.product(*List_of_Possible_Realizations)
	# Outcomes = tuple(Outcomes)
	
	
	
	return