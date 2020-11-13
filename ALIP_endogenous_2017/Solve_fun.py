import os
import sys
import itertools
import time
import pdb
import ALIP_EN

def de(MD):
	
	####################################################################
	###					Construct Scenarios
	####################################################################
	List_of_Possible_Realizations = [range(len(MD.Realizations['Qgrec'][r])) for r in MD.Realizations['Qgrec']]
	
	# print(List_of_Possible_Realizations)
	ALIP_EN.de(MD.parameters, MD.sets)

	### Generate Set of Outcomes
	Outcomes = itertools.product(*List_of_Possible_Realizations)
	Outcomes = tuple(Outcomes)
	
	
	
	return