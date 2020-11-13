import os
import sys
from pyomo.environ import *
from pyomo.opt import SolverFactory
import pdb
import numpy
import itertools
import datetime
import time
import Solve_fun
import importlib as ipl

def solve_function(solve_options):
	
	file_name = solve_options[-1]


	### Import Problem Data

	MD = get_data_info(file_name)

		
	### Set output directory
	current_directory = os.path.dirname(os.path.realpath(__file__))
	current_date = time.strftime('%m_%d_%Y', time.gmtime())

	output_directory = current_directory + '/Solutions/' + str(file_name) + '/' +  current_date + '/'

	Solve_fun.de(MD)

	
	return
	
def get_data_info(optcmd):
				
	## If the file is in the problem files folder then the filename is the sub-directory plus the filename
	try:
		fn = 'Cases.' + optcmd
		MD = ipl.import_module(fn) 
			
	### Otherwise assume the file is in the current directory
	except:
		try:	
			MD = ipl.import_module(optcmd)
		
		except:
			print("Error in File Import")
	

	return MD
	
if __name__ == "__main__":
	solve_options = sys.argv
	# print('solve_options',solve_options)
	solve_function(solve_options)
