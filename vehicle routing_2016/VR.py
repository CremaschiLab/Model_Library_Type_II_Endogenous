import os
import sys
from pyomo.environ import *
from pyomo.opt import SolverFactory
import pdb
import numpy
import itertools
import datetime
import pandas as pd
from collections import defaultdict

today = datetime.date.today()

def de(parameters, sets, Uncertain, Realizations, make):


	# uncertain parameters
	d = parameters['d']
	if d == {}:
		dd= pd.read_csv("d.csv", header= None)  
		d = dict((tuple((a, b)), c) for a,b,c in dd.values)



	### define the distinguishble scenario sets for endogenous uncertainties:
	opt = SolverFactory("cplex")

	#############################################################################################
	#############################################################################################
	# concrete model
	model = ConcreteModel()
	#############################################################################################
	#############################################################################################
	#############################################################################################
	# sets
	S = sets['S']
	J = sets['J']
	K = sets['K']
	I = sets['I']

	model.K = K
	model.J = J
	model.S = S

	f = make['f']
		if f == {}:
		for j in J:
			for jp in J:
				if j<=0 and jp<=0:
					f[j,jp] = 0
				else:
					f[j,jp] = 1
	#############################################################################################
	#############################################################################################
	# Parameters
	R = parameters['R']
	r_set = parameters['r_set']
	Q = parameters['Q']
	
	#############################################################################################
	#############################################################################################

	p = make['p']
	if p == {}:
		for s in S:
			p[s] = 1/len(S)

	Iend = I[-1]
	Kend = K[-1]
	#############################################################################################
	#############################################################################################
	A = {}
	# (5)
	for j in I:
		for r in [1,2,3]:
			for k in range(2*r, Iend+r+1):
				A[j,-r,k]= 1
				

	# (6)
	for r in [1,2,3]:
		for jp in I:
			for k in range(2*r +1, Iend+r+1):
				A[-r,jp,k]= 1
				
				
	# (2)
	for j in J:
		if j<=Iend and j>=1:
						A[0,j,1] = 1

	# (3)					
	for j in J:
		if j == -R:
			A[j,0,Kend] = 1
		if j<=Iend and j>=1:
			A[j,0,Kend] = 1

	# (4)					
	for (j,jp) in r_set:
		A[j,jp,-j+Iend+1] = 1
						
	for j in [1,2,3]:
		for jp in [1,2,3]:
			for k in K:
				if j != jp and k>1:
					A[j,jp,k] = 1

	D = {}
	for s in S:
		for sp in S:
			for j in I:
				if d[j,s] != d[j,sp]:
					try:
						D[s,sp].append(j)
					except KeyError:
						D[s,sp] = [j]
					
	n = {}
	for s in S:
		for sp in S:
			if s != sp:
				n[s,sp] = Iend - len(D[s,sp])
				


	kss = {}
	for s in S:
		for sp in S:
			if s != sp:
				kss[s,sp] = min(2*n[s,sp], n[s,sp] + R) + 1
			if s == sp:
				kss[s,sp] = Kend

	AA = {(j,jp,k) for (j,jp,k) in A}
	A_j_k = {(a,c) for a,b,c in AA}
	A_jp_k = {(b,c) for a,b,c in AA}

	# print(A_j_k)
	# print(A_jp_k)

	JK = {}
	for j in J:
		for k in K:
			if (j,k) in A_j_k:
				if (j,k) in A_jp_k:
					JK[j,k] = 1
				
	# testj = {(a,b) for a,b in JK}
	# print(testj)
					
	model.A_set = Set(initialize=[(j,jp,k) for (j,jp,k) in A])

	model.JK = Set(initialize=[(j,jp) for (j,jp) in JK])
	model.x = Var(model.A_set, model.S, within=Binary)

	#############################################################################################
	#############################################################################################
	# Variable

	model.y = Var(model.K, model.S, within=NonNegativeReals)

		
	#############################################################################################
	# objective function
	#############################################################################################
		
	def objectivemodel(model):
		return sum(p[s]*sum(f[j,jp]*model.x[j,jp,k,s] for j in model.J for jp in model.J for k in model.K if (j,jp,k) in model.A_set) for s in model.S)
	model.objective = Objective(sense=minimize, rule=objectivemodel)


	def con_2(model,s):
		return sum(model.x[0,j,1,s] for j in model.J if (0,j,1) in AA) == 1
	model.con_2m = Constraint(model.S, rule = con_2)
		
	def con_3(model,s):
		return sum(model.x[j,0,Kend,s] for j in model.J if (j,0,Kend) in AA) == 1
	model.con_3m = Constraint(model.S, rule = con_3)

	def con_4(model,j,s):
		if j!= 0:
			return sum(model.x[jp,j,k,s] for jp in model.J for k in model.K if (jp,j,k) in AA) == 1
		else:
			return Constraint.Skip
	model.con_4m = Constraint(model.J, model.S, rule = con_4)

	def con_5(model,j,s):
		if j!= 0:
			return sum(model.x[j,jp,k,s] for jp in model.J for k in model.K if (j,jp,k) in AA) == 1
		else:
			return Constraint.Skip
	model.con_5m = Constraint(model.J, model.S, rule = con_5)


	def con_6ddd(model,j,k,s):
		if j!= 0 and k < Kend:
			return sum(model.x[jp,j,k,s] for jp in model.J if (jp,j,k) in model.A_set) == sum(model.x[j,jp,k+1,s] for jp in model.J if (j,jp,k+1) in model.A_set)
		else:
			return Constraint.Skip
	model.con_6m = Constraint(model.JK, model.S, rule = con_6ddd)

	# def con_7(model,j,s):
		# if j!= 0:
			# return sum(k*sum(model.x[j,jp,k,s] for jp in model.J if (j,jp,k) in AA) for k in K if k>=2 and (j,jp,k) in AA) - sum(k*sum(model.x[jp,j,k,s] for jp in model.J if (jp,j,k) in AA) for k in K if (jp,j,k) in AA) == 1
		# else:
			# return Constraint.Skip
	# model.con_7m = Constraint(model.JK, rule = con_7)	

	def con_8(model,jp,s):
		if jp <0 and jp>-R:
			return sum(k*model.x[j,jp,k,s] for j in model.J for k in model.K if (j,jp,k) in AA) <= sum(k*model.x[j,jp-1,k,s] for j in model.J for k in model.K if (j,jp-1,k) in AA)
		else:
			return Constraint.Skip
	model.con_8m = Constraint(model.J, model.S, rule = con_8)	

	def con_9(model,s):
		return model.y[1,s] == Q - sum(d[j,s]*model.x[0,j,1,s] for j in model.J if (0,j,1) in AA)
	model.con_9m = Constraint(model.S, rule = con_9)	
		
		
	def con_10(model,k,s):
		if k>1:
			return model.y[k,s] <= model.y[k-1,s] - sum(d[jp,s]*model.x[j,jp,k,s] for j in model.J for jp in model.J if (j,jp,k) in AA)
		else:
			return Constraint.Skip
	model.con_10m = Constraint(model.K, model.S, rule = con_10)	

	def con_11(model,k,s):
		if k>1:
			return model.y[k,s] <= Q
		else:
			return Constraint.Skip
	model.con_11m = Constraint(model.K, model.S, rule = con_11)	


	###### check!!!!
	# def con_12(model,k,s):
		# if k>1:
			# return model.y[k,s] >= Q*sum(model.x[j,jp,k,s] for j in model.J for jp in model.J if jp<=0 and (j,jp,k) in AA)
		# else:
			# return Constraint.Skip
	# model.con_12m = Constraint(model.K, model.S, rule = con_12)	

	#############################################################################################
	# NACs
	#############################################################################################
	def NAC1(model,j,s,sp):
		if (0,j,1) in AA:
			return model.x[0,j,1,s] == model.x[0,j,1,sp]
		else:
			return Constraint.Skip
	model.NAC1m = Constraint(model.J, model.S, model.S,  rule = NAC1)	


	def NAC2(model,jpp,jppp,k,s,sp):
		if k>1 and k<= kss[s,sp] and s<sp:
			return model.x[jpp,jppp,k,s] - model.x[jpp,jppp,k,sp] <= sum(model.x[j,jp,kp,s] for kp in model.K if kp<k for j in model.J for jp in D[s,sp] if (j,jp,kp) in model.A_set)
		else:
			return Constraint.Skip
	model.NAC2m = Constraint(model.A_set, model.S, model.S, rule = NAC2)

	def NAC22(model,jpp,jppp,k,s,sp):
		if k>1 and k<= kss[s,sp] and s<sp:
			return -model.x[jpp,jppp,k,s] + model.x[jpp,jppp,k,sp] <= sum(model.x[j,jp,kp,s] for kp in model.K if kp<k for j in model.J for jp in D[s,sp] if (j,jp,kp) in model.A_set)
		else:
			return Constraint.Skip
	model.NAC22m = Constraint(model.A_set, model.S, model.S, rule = NAC22)

	results= opt.solve(model)
	model.solutions.load_from(results)
	save_file1 = "Vehicle routing" 
	results.write(filename=save_file1)

	save_file = "binary"
	f = open(os.path.join(save_file),	"w")


	for j in model.J:
		for jp in model.J:
			for k in model.K:
				for s in model.S:
					if (j,jp,k) in model.A_set and model.x[j,jp,k,s].value >0:
						f.write(str((j,jp,k,s)) + ' ' + str(model.x[j,jp,k,s].value) + '\n')
						
	f.write('--------------------------------------------------------------------------------------------')
						
	f.close
	
	return