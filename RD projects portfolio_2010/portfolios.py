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

	Tend = sets['Tend']
	T = range(1, Tend+1)
	outcomes = Realizations['outcomes']
	outcomes_list = list(itertools.product(*outcomes))
	S_numbers = range(1, len(outcomes_list)+1)

	# sets
	I = sets['I']
	
	### define the endogenous uncertainty values:
	theta = {}
	Z = {}
	Z_intermediate = {}
	for s in S_numbers:
		count = 1
		for i in I:
			theta[i,s] = outcomes_list[s-1][2*count-2]
			Z_intermediate[i,s] =  outcomes_list[s-1][2*count-1]
			Z[i,s] = outcomes_list[s-1][2*count-1]
			count += 1
			
	H = {}
	Y = {}
	for i in I:
		for s in S_numbers:
			for sp in S_numbers:
				# if s<sp:
				H[s,sp] = []
				Y[s,sp] = []
				if Z_intermediate[i,s] != Z_intermediate[i,sp]:
					H[s,sp].append(i)
				if theta[i,s] != theta[i,sp]:
					Y[s,sp].append(i)
				

	opt = SolverFactory("cplex")

	#############################################################################################
	#############################################################################################
	# concrete model
	model = ConcreteModel()

	#############################################################################################
	#############################################################################################
	# Parameters
	delta =  parameters['delta']
	f = parameters['f']
	r = parameters['r']
	theta_theta = parameters['theta_theta'] 
	theta_Z = parameters['theta_Z']
	D =  parameters['D'] 
	M = parameters['Big_M']
	#############################################################################################
	#############################################################################################
	# Generate Parameters
	pp ={}
	if make['pp'] == []:
		for s in S_numbers:
			pp[s] = 1/len(S_numbers)
	else:
		pp = make['pp']
		
	delta_ave = {}
	if make['delta_ave'] == []:
		for i in I:
			for j in I:
				delta_ave[i,j] = max(delta[i], delta[j])
	else:
		delta_ave = make['delta_ave']
		
	delta_max = max(delta.values())
	
	B = {}
	if make['B'] == []:
		for t in T:
			B[t] = 3
	else:
		B = make['B']
		
	#############################################################################################
	#############################################################################################
	# sets
	model.T = T
	model.S = S_numbers
	model.I = I

	#############################################################################################
	#############################################################################################
	# Variable
	model.x = Var(model.I, model.T, model.S, within=NonNegativeReals)
	model.tao = Var(model.I, model.T, model.S, within=NonNegativeReals)
	
	model.y = Var(model.I, model.T, model.S, within=Binary)
	model.h = Var(model.I, model.T, model.S, within=Binary)
	model.alpha = Var(model.I, model.T, model.S, within=Binary)
	
	model.beta = Var(model.I, range(min(delta.values()), Tend + max(delta.values())+1), model.S, within=Binary)
	model.gama = Var(model.I, range(1, Tend+1), model.S, within=Binary)
	model.sai = Var(model.I, model.I, range(1,Tend + max(delta.values())+1), model.S, within=NonNegativeReals)	

	model.obj1 = Var(model.I,  model.S, within=Reals)
	model.obj2 = Var(model.I,  model.S, within=Reals)
	model.M = M
	

	#############################################################################################
	# objective function
	#############################################################################################
		
	def objectivemodel(model):
		return sum(pp[s]*sum(model.obj1[i,s] + model.obj2[i,s] for i in model.I) for s in model.S)
	model.objective = Objective(sense=maximize, rule=objectivemodel)

		### obj constraints
	def objpart1(model,i,s):
		return model.obj1[i,s]== sum(model.beta[i,t+delta[i],s] * Z[i,s]*(1+r)**(-t-delta[i]) for t in model.T if t<= Tend - 1) + \
									model.beta[i, Tend+delta[i], s] * Z[i,s] * ((1+r)**(-Tend-delta_max)/r + sum((1+r)**(-Tend-delta[i]-l) for l in range(0, delta_max - delta[i] )))
	model.objpart1m = Constraint(model.I, model.S, rule = objpart1)
		
	def objpart2(model,i,s):
		return model.obj2[i,s] == sum(sum(model.sai[i,j,t + delta_ave[i,j],s] * 1 * (1+r)**(-t-delta_ave[i,j]) for t in model.T if t <= Tend - 1) + \
										   model.sai[i,j,Tend + delta_ave[i,j],s]* 1 * ((1+r)**(-Tend-delta_max)/r + sum((1+r)**(-Tend-delta_ave[i,j]-l) for l in range(0, delta_max - delta_ave[i,j]))) \
										for j in D[i] if j>i)
	model.objpart12 = Constraint(model.I, model.S, rule = objpart2)
		
	def con_3(model,i,t,s):
		return model.alpha[i,t,s] - model.beta[i,t+delta[i],s] >= 0
	model.con_3m = Constraint(model.I, model.T, model.S, rule = con_3)
		
	def con_4(model,i,t,s):
		return sum(model.x[i,tprime,s] for tprime in model.T if tprime<=t) - max(theta[i,s] + t*f[i], max(B[tp] for tp in range(1,t+1))) * model.alpha[i,t,s] <= 0
	model.con_4m = Constraint(model.I, model.T, model.S, rule = con_4)
	
	def con_5(model,t,s):
		return sum(model.x[i,t,s] for i in model.I) <= B[t]
	model.con_5m = Constraint(model.T, model.S, rule = con_5)	

	def con_6(model,i,t,s):
		return model.x[i,t,s] - B[t]*(model.alpha[i,t,s] - model.beta[i,t+delta[i]-1,s] - model.gama[i,t,s]) <=0
	model.con_6m = Constraint(model.I, model.T, model.S, rule = con_6)
	
	def con_7(model,i,j,t,s):
		if j>i and j in D[i]:
			return model.beta[i,t+delta_ave[i,j],s] + model.beta[j,t+delta_ave[i,j],s] - model.sai[i,j,t + delta_ave[i,j],s] <= 1
		else:
			return Constraint.Skip
	model.con_7m = Constraint(model.I, model.I, model.T, model.S, rule = con_7)	
	
	def con_8(model,i,j,t,s):
		if j>i and j in D[i]:
			return model.beta[i,t+delta_ave[i,j],s] + model.beta[j,t+delta_ave[i,j],s] - 2*model.sai[i,j,t + delta_ave[i,j],s] >= 0
		else:
			return Constraint.Skip
	model.con_8m = Constraint(model.I, model.I, model.T, model.S, rule = con_8)	

	def con_9(model,i,j,t,s):
		if t>= 2:
			return model.tao[i,t,s] - model.tao[i,t-1,s] + model.x[i,t,s] - f[i]*(model.alpha[i,t,s] - model.beta[i,t+delta[i]-1,s] - model.gama[i,t,s]) >= 0
		else:
			return Constraint.Skip
	model.con_9m = Constraint(model.I, model.I, model.T, model.S, rule = con_9)	

	def con_10(model,i,t,s):
		return model.x[i,t,s] - f[i]*(model.alpha[i,t,s] - model.beta[i,t+delta[i]-1,s] - model.gama[i,t,s]) >= 0
	model.con_10m = Constraint(model.I,  model.T, model.S, rule = con_10)	
	
	def con_11(model,i,t,s):
		return model.tao[i,t,s] + theta[i,s]*model.beta[i,t+delta[i],s] <= theta[i,s]
	model.con_11m = Constraint(model.I, model.T, model.S, rule = con_11)

	#############################################################################################
	# define indicator variables
	#############################################################################################
	def con_12(model,i,t,s):
		if t != 1:
			return sum(model.x[i,tprime,s] - f[i]*(model.alpha[i,tprime,s] - model.beta[i,tprime+delta[i]-1,s] - model.gama[i,tprime,s]) for tprime in model.T if tprime<= t) - theta_theta[i]*model.y[i,t,s] >=0
		else:
			return Constraint.Skip
	model.con_12m = Constraint(model.I, model.T, model.S, rule = con_12)
	
	def con_13(model,i,t,s):
		if t != 1:
			return sum(model.x[i,tprime,s] for tprime in model.T if tprime<= t) - ( min(sum(B[tp] for tp in range(1,t+1)), theta[i,s] + (t-1)*f[i]) - theta_theta[i]) * model.y[i,t,s] <= theta_theta[i]
		else:
			return Constraint.Skip
	model.con_13m = Constraint(model.I, model.T, model.S, rule = con_13)	

	
	def con_14(model,i,t,s):
		if t != 1:
			return sum(model.x[i,tprime,s] - f[i]*(model.alpha[i,tprime,s] - model.beta[i,tprime+delta[i]-1,s] - model.gama[i,tprime,s]) for tprime in model.T if tprime<= t) - theta_Z[i]*model.h[i,t,s] >=0
		else:
			return Constraint.Skip
	model.con_14m = Constraint(model.I, model.T, model.S, rule = con_14)
	
	def con_15(model,i,t,s):
		if t != 1:
			return sum(model.x[i,tprime,s] for tprime in model.T if tprime<= t) -( min(sum(B[tp] for tp in range(1,t+1)), theta[i,s] + (t-1)*f[i]) - theta_Z[i]) * model.h[i,t,s] <= theta_Z[i]
		else:
			return Constraint.Skip
	model.con_15m = Constraint(model.I, model.T, model.S, rule = con_15)
	
	#############################################################################################
	# define indicator variables
	#############################################################################################
	def con_16(model,i,t,s):
		if t != 1:
			return  model.beta[i,t+delta[i],s] + model.gama[i,t,s] <=1
		else:
			return Constraint.Skip
	model.con_16m = Constraint(model.I, model.T, model.S, rule = con_16)	
	
	def con_17(model,i,t,s):
		if t+1 <= Tend:
			return  model.alpha[i,t,s] - model.gama[i,t+1,s] >= 0
		else:
			return Constraint.Skip
	model.con_17m = Constraint(model.I, model.T, model.S, rule = con_17)	
	
	#############################################################################################
	# first stage nonanticipativity constriants
	#############################################################################################
	def con_18(model,i,s):
		return  model.x[i,1,s] - sum(pp[sprime]*model.x[i,1,sprime] for sprime in model.S) ==0
	model.con_18m = Constraint(model.I, model.S, rule = con_18)	
	
	def con_19(model,i,t,s,sprime):
		if t != 1:
			return  model.x[i,t,s] - model.x[i,t,sprime] + B[t]*(sum(model.y[j,t,s] + model.y[j,t,sprime] for j in Y[s,sprime]) + sum(model.h[j,t,s] + model.h[j,t,sprime] for j in H[s,sprime])) >=0
		else:
			return Constraint.Skip
	model.con_19m = Constraint(model.I, model.T, model.S, model.S, rule = con_19)
	
	
	results = opt.solve(model)
	model.solutions.load_from(results)
	
	save_file1 = "SMILP" + '_' + str(T[-1]) + ' ' + '(' + str(today) + ')'
	save_file = "SMILPOutput" + '_' + str(T[-1]) + ' ' + '(' + str(today) + ')'
	dir_path = os.path.dirname(os.path.realpath(__file__))
	# TempfileManager.tempdir = dir_path

	f = open(os.path.join(save_file),	"w")
	results.write(filename=save_file1)

	for i in I:
		for t in T:
			for s in S_numbers:
				f.write('x' + ' ' + str((i,t,s)) + ' ' + str(model.x[i,t,s].value) + '\n')

	f.write('-------------------------------------------------------------------------------------'+ '\n')

	for i in I:
		for t in T:
			for s in S_numbers:
				if model.alpha[i,t,s].value >0: 
					f.write('alpha' + ' ' + str((i,t,s)) + ' ' + str(model.alpha[i,t,s].value) + '\n')
	return