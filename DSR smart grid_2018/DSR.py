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

	opt = SolverFactory("cplex")

	#############################################################################################
	#############################################################################################
	# concrete model
	model = ConcreteModel()
	#############################################################################################
	#############################################################################################
	#############################################################################################
	# sets
	Omega_S = sets['Omega_S']
	Omega_DG = sets['Omega_DG']
	Omega_E = sets['Omega_E']
	Omega_Es =sets['Omega_Es']
	Omega_K = sets['Omega_K']
	Omega_L = sets['Omega_L']
	Omega_G = sets['Omega_G']
	Omega_N = sets['Omega_N']
	Omega_O = sets['Omega_O']
	Omega_T = sets['Omega_T']
	#############################################################################################
	#############################################################################################
	# Parameters
	rI = parameters['rI']
	rO = parameters['rO']
	gamaD = parameters['gamaD']
	gama = parameters['gama']
	w = parameters['w']
	cDSR = parameters['cDSR']
	cDG = parameters['cDG']
	K = parameters['K']
	sai = parameters['sai']
	kl = parameters['kl']
	kd = parameters['kd']
	Dmax_load = parameters['Dmax_load']
	d = parameters['d']
	X = parameters['X']
	F = parameters['F']
	delta = parameters['delta']
	
	# uncertain parameters
	f  = parameters['f']
	
	# differentiate sets
	D = parameters['D']

	#############################################################################################

	# Variable
	model.faiI = Var(Omega_S, within=Reals)
	model.faio = Var(Omega_S, within=Reals)
	model.B = Var(Omega_S, Omega_E, Omega_L, Omega_O, within = Binary)
	model.D = Var(Omega_S, Omega_E, Omega_N, within = Binary)
	model.Xi = Var(Omega_S, Omega_E, Omega_T, Omega_N, within=Reals)
	model.F_a = Var(Omega_S, Omega_E, Omega_L, Omega_O, within=Reals)
	model.F = Var(Omega_S, Omega_E, Omega_L, Omega_O, within=Reals)
	model.D_a = Var(Omega_S, Omega_E, Omega_N, within=Reals)
	model.D = Var(Omega_S, Omega_E, Omega_N, within=Reals)
	model.Z = Var(Omega_E, Omega_S, Omega_S, within=Binary)
	model.P_g = Var(Omega_S, Omega_E, Omega_T, Omega_G, within=Reals)
	model.P_l = Var(Omega_S, Omega_E, Omega_T, Omega_L, within=Reals)
	model.theta_u1 = Var(Omega_S, Omega_E, Omega_T, within=Reals)
	model.theta_v1 = Var(Omega_S, Omega_E, Omega_T, within=Reals)
	
	#############################################################################################
	# objective function
	#############################################################################################
		
	def objectivemodel(model):
		return sum(p[s] * (model.faiI[s] + model.faio[s]) for s in Omega_S)
	model.objective = Objective(sense=minimize, rule=objectivemodel)


	def con_2(model,s):
		return model.faiI[s] == sum(rI[e]*(sum(model.B[s,e,l,o]*gama[o] for l in Omega_L for o in Omega_O) + sum(model.D[s,e,n]*gamaD for n in Omega_N)) for e in Omega_E)
	model.con_2m = Constraint(Omega_S, rule = con_2)
		
	def con_3(model,s):
		return model.faio[s] == sum(rO[e]*w[t]*model.Xi[s,e,t,n]*cDSR for e in Omega_E for t in Omega_T for n in Omega_N) + \
								sum(rO[e]*cDG*(K[s,e,g]*sai[t,g] - model.P_g[s,e,t,g]) for e in Omega_E for t in Omega_T for g in Omega_DG)
	model.con_3m = Constraint(Omega_S, rule = con_3)
	
	def con_4(model,s,e,l,o):
		return model.F_a[s,e,l,o] == sum(model.F[s,z,l,o] for z in range(1, e-kl+1))
	model.con_4m = Constraint(Omega_S, Omega_E, Omega_L, Omega_O, rule = con_4)
	
	def con_5(model,s,e,n):
		return model.D_a[s,e,n] == sum(model.D[s,z,n] for z in range(1, e-kd+1))
	model.con_5m = Constraint(Omega_S, Omega_E, Omega_N, rule = con_5)
	
	def con_6(model,s,e,l,o):
		return model.F[s,e,l,o] == Q[o] * model.B[s,e,l,o]
	model.con_6m = Constraint(Omega_S, Omega_E, Omega_L, Omega_O, rule = con_6)

	#############################################################################################
	# Initial NACs
	#############################################################################################
	def con_7(model,s,sp,l,o):
		return model.B[s,1,l,o] == model.B[sp,1,l,o]
	model.con_7m = Constraint(Omega_S, Omega_S, Omega_L, Omega_O, rule = con_7)

	def con_8(model,s,sp,n):
		return model.D[s,1,n] == model.D[sp,1,n]
	model.con_8m = Constraint(Omega_S, Omega_S, Omega_N, rule = con_8)
	#############################################################################################
	# Equation 9
	#############################################################################################
	def con_18(model,s,sp,n):
		if n in D[s,sp] and s<sp:
			return model.Z[e,s,sp] <= 1- model.D_a[s,e,n]
		else:
			return Constraint.Skip
	model.con_18m = Constraint(Omega_S, Omega_S, Omega_N, rule = con_18)
	
	def con_19(model,s,sp,e):
		return model.Z[e,s,sp] >= sum((1 - model.D_a[s,e,n]) for n in D[s,sp]) - (len(D[s,sp]) - 1)
	model.con_19m = Constraint(Omega_S, Omega_S, Omega_E, rule = con_19)

	#############################################################################################
	# Conditional NACs Equation 10
	#############################################################################################
	def con_20(model,s,sp,l,o,e):
		if s != sp:
			return model.B[s,e+1,l,o] - model.B[sp,e+1,l,o] + M*(1 - model.Z[e,s,sp]) >= 0
		else:
			return Constraint.Skip
	model.con_20m = Constraint(Omega_S, Omega_S, Omega_L, Omega_O, Omega_Es, rule = con_20)
	
	def con_21(model,s,sp,n,e):
		if s != sp:
			return model.D[s,e+1,n] - model.D[sp,e+1,n] + M*(1 - model.Z[e,s,sp]) >= 0
		else:
			return Constraint.Skip
	model.con_21m = Constraint(Omega_S, Omega_S, Omega_N, Omega_Es, rule = con_21)
	#############################################################################################
	def con_11(model,s,e,t,n):
		return model.T[s,e,t,n] <= model.D_a[s,e,n] * Dmax_load[n]
	model.con_11m = Constraint(Omega_S, Omega_E, Omega_T, Omega_N, rule = con_11)
	
	def con_12(model,s,e,t,n):
		return model.Xi[s,e,t,n] <= model.D_a[s,e,n] * f[s,n] * d[s,e,t,n]
	model.con_12m = Constraint(Omega_S, Omega_E, Omega_T, Omega_N, rule = con_12)
	
	def con_13(model,s,e,n):
		return sum((model.T[s,e,k,n] - model.Xi[s,e,k,n]) for k in range(1,k+delta))
	model.con_13m = Constraint(Omega_S, Omega_E, Omega_N, rule = con_13)
	
	def con_14(model,s,e,t,g):
		return model.P_g[s,e,t,g] <= K[s,e,g]*sai[t,g]
	model.con_14m = Constraint(Omega_S, Omega_E, Omega_T, Omega_G, rule = con_14)
	
	def con_15(model,s,e,t,l):
		return model.P_l[s,e,t,g] <= (model.theta_u1[s,e,t] - model.theta_v1[s,e,t])/X[l]
	model.con_15m = Constraint(Omega_S, Omega_E, Omega_T, Omega_L, rule = con_15)
	
	def con_16b(model,s,e,t,l):
		return model.P_l[s,e,t,g] >= -F[l] - sum(model.F_a[s,e,l,o] for o in Omega_O)
	model.con_16bm = Constraint(Omega_S, Omega_E, Omega_T, Omega_L, rule = con_16b)

	def con_17(model,s,e,t,l):
		return sum(model.P_g[s,e,t,g] for g in Omega_G) + sum(model.P_l[s,e,t,l]*L[n,t] for l in Omega_L) == model.T[s,e,t,n] - model.Xi[s,e,t,n] + d[s,e,t,n]
	model.con_17m = Constraint(Omega_S, Omega_E, Omega_T, Omega_L, rule = con_17)

	return