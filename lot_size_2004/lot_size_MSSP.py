import os
import sys
from pyomo.environ import *
from pyomo.opt import SolverFactory
import pdb
import itertools
import datetime
import time
from random import randint
today = datetime.date.today()


def de(parameters, sets, Uncertain, Realizations, make):

	# generate solution file
	current_directory = os.path.dirname(os.path.realpath(__file__))
	current_date = time.strftime('%m_%d_%Y', time.gmtime())
	output_directory = current_directory + '/Solutions/' + current_date + '/'
	
	
	outcomes = Realizations['outcomes']
	outcomes_list = list(itertools.product(*outcomes))
	S_number = range(0, len(outcomes_list))

	start_time = time.time()

	TT = sets['T']
	I = sets['I']
	I_uncertain = Uncertain['I'] 
	#############################################################################################
	#############################################################################################
	# concrete model
	model = ConcreteModel()
	#############################################################################################
	#############################################################################################
	# sets
	model.I = I
	model.T = TT
	model.S = S_number
	#############################################################################################
	#############################################################################################
	# generated Parameters
	c = {}
	if make['c'] == []:
		for t in TT:
			for s in S_number:
				c[t,s] = 30000
	else:
		c = make['c']
	
	pro = {}
	if make['probability'] == []:
		for s in S_number:
			pro[s] = 1/len(S_number)
	else:
		pro = make['probability']

	Demand = {}
	if make['Demand'] == []:
		for s in model.S:
			for t in TT:
				Demand[t, s] = randint(6000, 8000)
	else:
		Demand = make['Demand']
	#############################################################################################
	#############################################################################################
	M = parameters['Big_M']
	rho = parameters['rho']
	sigma = parameters['sigma']
	# unit production cost for size i. P[i,s]
	p = {}
	for s in S_number:
		p[1,s] = 1
		p[2,s] = outcomes_list[s][0]
		p[3,s] = outcomes_list[s][1]

	save_file = "scenario_specific"
	f = open(os.path.join(save_file),	"w")
	f.write('------------------------------------------production cost ----------------------------------------------' + '\n')
	for i in I_uncertain:
		f.write('----------------------------------' + '\n')
		for s in S_number:
			f.write('p' + ' ' + str((i,s)) + ' ' + str(p[i,s]) + '\n')

			
	f.write('------------------------------------------demand ----------------------------------------------' + '\n')
	for i in I:
		f.write('----------------------------------' + '\n')
		for s in S_number:
			f.write('demand' + ' ' + str((i,s)) + ' ' + str(Demand[i,s]) + '\n')
	#############################################################################################
	#############################################################################################
	D = {}
	for s in S_number:
		for sp in S_number:
			if s<sp:
				D[s,sp] = []
				for i in I:
					if p[i,s] != p[i,sp]:
						D[s,sp].append(i)

	L1 = {}
	for s in S_number:
		for sp in S_number:
			if s<sp:
				if len(D[s,sp]) == 1:
					L1[s,sp] = 1
			
	f.write('------------------------------------------D[s,sp] ----------------------------------------------' + '\n')
	for s in S_number:
		for sp in S_number:
			if s<sp:
				f.write('D[s,sp]' + ' ' + str((s,sp)) + ' ' + str(D[s,sp]) + '\n')	
	#############################################################################################
	#############################################################################################
	# Variable
	model.z = Var(model.I,model.T,model.S, within=Binary)
	model.y = Var(model.I,model.T,model.S, within=NonNegativeIntegers)
	model.x = Var(model.I,model.I,model.T,model.S, within=NonNegativeIntegers)
	model.Z_diff = Var(model.T,model.S,model.S, within=Binary)
		
	#############################################################################################
	# objective function
	#############################################################################################
		
	def objectivemodel(model):
		return sum(pro[s]* sum(( sum(sigma*model.z[i,t,s] + p[i,s]*model.y[i,t,s] for i in model.I) + rho*sum(model.x[i,j,t,s] for i in model.I for j in model.I if j<=i)) for t in model.T) for s in model.S)
	model.objective = Objective(sense=minimize, rule=objectivemodel)

		
	def con_1(model,j,t,s):
		return sum(model.x[i,j,t,s] for i in model.I if i>=j) >= Demand[t,s]
	model.con_1m = Constraint(model.I, model.T, model.S, rule = con_1)
		
	def con_2(model,i,t,s):
		return sum(sum(model.x[i,j,tprim,s] for j in model.I if j<=i) - model.y[i,tprim,s] for tprim in model.T if tprim <=t) <= 0
	model.con_2m = Constraint(model.I, model.T, model.S,  rule = con_2)
	
	def con_3(model,i,t,s):
		return model.y[i,t,s] - M*model.z[i,t,s] <= 0
	model.con_3m = Constraint(model.I, model.T, model.S, rule = con_3)
	
	def con_4(model,t,s):
		return sum(model.y[i,t,s] for i in model.I)<= c[t,s]
	model.con_4m = Constraint(model.T, model.S, rule = con_4)
	
	# initial NACs
	def con_5(model,i,s,sprime):
		return model.y[i,1,s] == model.y[i,1,sprime]
	model.con_5m = Constraint(model.I, model.S, model.S, rule = con_5)
	
	def con_6(model,i,s,sprime):
		return model.z[i,1,s] == model.z[i,1,sprime]
	model.con_6m = Constraint(model.I, model.S, model.S, rule = con_6)
	
	## conditonal endogenous NACs
	def con_indicate1(model,i,t,tao,s,sp):
		if s<sp and (s,sp) in L1 and tao<=t and i in D[s,sp]:
			return model.Z_diff[t,s,sp] <= 1 - model.z[i,tao,s] 
		else:
			return Constraint.Skip
	model.indicate1m = Constraint(model.I, model.T, model.T, model.S, model.S, rule = con_indicate1)
	
	def con_indicate2(model,i,t,s,sp):
		if s < sp and (s,sp) in L1 and i in D[s,sp]:
			return model.Z_diff[t,s,sp] >= 1 - sum(model.z[i,tao,s] for tao in range(1,t+1))
		else:
			return Constraint.Skip
	model.indicate2m = Constraint(model.I, model.T, model.S, model.S, rule = con_indicate2)
	
	def NAC1a(model, i, j, t, s, sp):
		if s < sp and (s,sp) in L1 and j<=i:
			return model.x[i,j,t,s] -  model.x[i,j,t,sp] <= 1 - model.Z_diff[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC1am = Constraint(model.I,model.I, model.T, model.S, model.S, rule = NAC1a)
	
	def NAC1b(model, i, j, t, s, sp):
		if s < sp and (s,sp) in L1 and j<=i:
			return model.x[i,j,t,s] -  model.x[i,j,t,sp] >= - 1 + model.Z_diff[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC1bm = Constraint(model.I,model.I, model.T, model.S, model.S, rule = NAC1b)
	
	def NAC2a(model, i, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1 <= TT[-1]:
			return model.y[i,t+1,s] -  model.y[i,t+1,sp] <= 1 - model.Z_diff[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC2am = Constraint(model.I, model.T, model.S, model.S, rule = NAC2a)
	
	def NAC2b(model, i, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1 <= TT[-1]:
			return model.y[i,t+1,s] -  model.y[i,t+1,sp] >= -1 + model.Z_diff[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC2bm = Constraint(model.I, model.T, model.S, model.S, rule = NAC2b)
	
	def NAC3a(model, i, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1 <= TT[-1]:
			return model.z[i,t+1,s] -  model.z[i,t+1,sp] <= 1 - model.Z_diff[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC3am = Constraint(model.I, model.T, model.S, model.S, rule = NAC3a)
	
	def NAC3b(model, i, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1 <= TT[-1]:
			return model.z[i,t+1,s] -  model.z[i,t+1,sp] >= -1 + model.Z_diff[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC3bm = Constraint(model.I, model.T, model.S, model.S, rule = NAC3b)

	opt = SolverFactory('cplex')
	results= opt.solve(model)
	model.solutions.load_from(results)
	
	end_time = time.time()
	solve_time = end_time - start_time

	### Open save file
	if not os.path.exists(output_directory):
		os.makedirs(output_directory)

	save_file = 'MSSP (size problem)' + ' ' + 'time horizon' + ' ' + str(TT[-1]) + ' ' + '(' + str(today) + ')'
	f = open(os.path.join(output_directory, save_file),	"w")
		
	MST =  'solving Time:' + ' ' + str(solve_time) + " " 

	f.write(MST + '\n')

	f.write('Objective value:' + ' ' + str(model.objective()) + '\n')

	for i in model.I:
		for j in model.I:
			for t in TT:
				for s in model.S:
					if model.x[i,j,t,s].value !=0:
						f.write('x ' +  ' ' + str((i,j,t,s)) + ' ' + str(model.x[i,j,t,s].value) + '\n')
		
	f.write('--------------------------------------------------------------' + '\n')
		
	for i in model.I:
		for t in TT:
			for s in model.S:
				f.write('y ' +  ' ' + str((i,t,s)) + ' ' + str(model.y[i,t,s].value) + '\n')
		
	f.write('--------------------------------------------------------------' + '\n')
		
	for i in model.I:
		for t in TT:
			for s in model.S:
				if model.z[i,t,s].value>0.9:
					f.write('z ' +  ' ' + str((i,t,s)) + ' ' + str(model.z[i,t,s].value) + '\n')
		

	return
	
