import os
import sys
from pyomo.environ import *
from pyomo.opt import SolverFactory
import pdb
import numpy
import itertools
import datetime

today = datetime.date.today()


def de(parameters, sets, Uncertain, Realizations, make):
	
	
	T = sets['T']
	I = sets['I']
	outcomes = Realizations['outcomes']

	outcomes_list = list(itertools.product(*outcomes))
	S_numbers = range(1, len(outcomes_list)+1)

	### define the distinguishble scenario sets for endogenous uncertainties:
	opt = SolverFactory("cplex")

	a0 = {}
	a1 = {}
	g = {}
	for s in S_numbers:
		for i in I:
			a0[i,s] = outcomes_list[s-1][2*i-2]
			a1[i,s] = outcomes_list[s-1][2*i-1]
			g[i,s] = a0[i,s]/a1[i,s]
			
	D = {}
	for i in I:
		for s in S_numbers:
			for sp in S_numbers:
				# if s<sp:
				D[s,sp] = []
				if g[i,s] != g[i,sp]:
					D[s,sp].append(i)
						

	
	#############################################################################################
	# concrete model
	model = ConcreteModel()
	#############################################################################################
	#############################################################################################
	#############################################################################################
	# sets

	model.T = T
	model.I = I
	model.S = S_numbers

	#############################################################################################
	#############################################################################################
	# Parameters
	pp = make['pp']
	c1 = make['c1']
	c_mng = make['c_mng']
	c_proc = make['c_proc']
	P_set = make['P_set'] 
	P = make['P']
	M = make['M']
	
	# check if enter the values, if not generate default values
	if pp == {}:
		for s in S_numbers:
			pp[s] = 1/len(S_numbers)

	for t in T:
		if t not in c_mng:
			c_mng[t] = 1
		if t not in c_proc:
			c_proc[t] = 9
		if t not in c1:
			c1[t] = 40
		if t not in M:
			M[t] = 1
		if t not in P:
			P[t] = 1
		
	if P_set == {}:
		P_set = {1:2, 2:3, 3:[]}
		
	#############################################################################################
	#############################################################################################
	# Variable
	model.x = Var(model.I, model.T, model.S, within = Binary)
	model.y = Var(model.I, model.T, model.S, within = NonNegativeReals, bounds = (0,1))
	model.z = Var(model.I, model.T, model.S, within = NonNegativeReals, bounds = (0,1))
		
	#############################################################################################
	# objective function
	#############################################################################################
		
	def objectivemodel(model):
		return sum(pp[s]*((c1[t]*a1[i,s] - c_proc[t]*a0[i,s]) * model.z[i,t,s] - a0[i,s]*c_mng[t]*model.y[i,t,s]) for i in model.I for t in model.T for s in model.S)
	model.objective = Objective(sense=maximize, rule=objectivemodel)

	def con_21(model,i,t,s):
		return model.z[i,t,s] <= model.y[i,t,s]
	model.con_21m = Constraint(model.I, model.T, model.S, rule = con_21)	
	
	def con_22(model,t,s):
		return sum(a0[i,s]*model.z[i,t,s] for i in model.I)  <= P[t]
	model.con_22m = Constraint(model.T, model.S, rule = con_22)
	
	def con_23(model,t,s):
		return sum(a0[i,s]*model.y[i,t,s] for i in model.I)  <= M[t]
	model.con_23m = Constraint(model.T, model.S, rule = con_23)

	def con_24(model,i,j,t,s):
		return model.x[i,t,s] <= sum(model.y[j,tprime,s] for tprime in range(1,t+1))
	model.con_24m = Constraint(model.I, P_set[i], model.T, model.S, rule = con_24)
	
	def con_25(model,i,t,s):
		return sum(model.y[i,tprime,s] for tprime in range(1,t+1)) <= model.x[i,t,s]
	model.con_25m = Constraint(model.I, model.T, model.S, rule = con_25)
	
	def con_26(model,i,t,s):
		if t+1 <= T[-1]:
			return model.x[i,t,s] <= model.x[i,t+1,s]
		else:
			return Constraint.Skip
	model.con_26m = Constraint(model.I, model.T, model.S, rule = con_26)
	
	#############################################################################################
	# NACs 31-33, 36-38, 41 and 42
	#############################################################################################
	def con_31(model,i,s,sp):
		return model.x[i,1,s] == model.x[i,1,sp] 
	model.con_31m = Constraint(model.I, model.S, model.S, rule = con_31)	
	
	def con_32(model,i,t,s,sp):
		if t>2:
			return model.x[i,t,s] - model.x[i,t,sp] <= sum(model.x[j,t-1,s] for j in D[sp,s])
		else:
			return Constraint.Skip
	model.con_32m = Constraint(model.I, model.T, model.S, model.S, rule = con_32)
	
	def con_33(model,i,t,s,sp):
		if t>2:
			return model.x[i,t,sp] - model.x[i,t,s] <= sum(model.x[j,t-1,s] for j in D[sp,s])
		else:
			return Constraint.Skip
	model.con_33m = Constraint(model.I, model.T, model.S, model.S, rule = con_33)
	
	#
	def con_36(model,i,s,sp):
		return model.y[i,1,s] == model.y[i,1,sp] 
	model.con_36m = Constraint(model.I, model.S, model.S, rule = con_36)	
	
	def con_37(model,i,t,s,sp):
		if t>2:
			return model.y[i,t,s] - model.y[i,t,sp] <= sum(model.x[j,t-1,s] for j in D[sp,s])
		else:
			return Constraint.Skip
	model.con_37m = Constraint(model.I, model.T, model.S, model.S, rule = con_37)
	
	def con_38(model,i,t,s,sp):
		if t>2:
			return model.y[i,t,sp] - model.y[i,t,s] <= sum(model.x[j,t-1,s] for j in D[sp,s])
		else:
			return Constraint.Skip
	model.con_38m = Constraint(model.I, model.T, model.S, model.S, rule = con_38)
	
	#
	def con_41(model,i,t,s,sp):
		if t>2:
			return model.z[i,t,s] - model.z[i,t,sp] <= sum(model.x[j,t,s] for j in D[sp,s])
		else:
			return Constraint.Skip
	model.con_41m = Constraint(model.I, model.T, model.S, model.S, rule = con_41)
	
	def con_42(model,i,t,s,sp):
		if t>2:
			return model.z[i,t,sp] - model.z[i,t,s] <= sum(model.x[j,t,s] for j in D[sp,s])
		else:
			return Constraint.Skip
	model.con_42m = Constraint(model.I, model.T, model.S, model.S, rule = con_42)
	
	
	results= opt.solve(model)
	model.solutions.load_from(results)
	
	save_file1 = "SMILP" + '_' + str(T[-1]) + ' ' + '(' + str(today) + ')'
	save_file = "SMILPOutput" + '_' + str(T[-1]) + ' ' + '(' + str(today) + ')'
	dir_path = os.path.dirname(os.path.realpath(__file__))


	f = open(os.path.join(save_file),	"w")
	results.write(filename=save_file1)
	
	for i in I:
		for t in T:
			for s in S_numbers:
				if model.x[i,t,s].value >0.9: 
					f.write('x' + ' ' + str((i,t,s)) + ' ' + str(model.x[i,t,s].value) + '\n')

	f.write('-------------------------------------------------------------------------------------')

	for i in I:
		for t in T:
			for s in S_numbers:
				if model.y[i,t,s].value >0: 
					f.write('y' + ' ' + str((i,t,s)) + ' ' + str(model.y[i,t,s].value) + '\n')

	return
