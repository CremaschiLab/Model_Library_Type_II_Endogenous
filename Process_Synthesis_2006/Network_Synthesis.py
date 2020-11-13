import os
import sys
from pyomo.environ import *
from pyomo.opt import SolverFactory
import pdb
import numpy
import itertools
import datetime


def de(parameters, sets, Uncertain, Realizations, make):

	opt = SolverFactory("cplex")
	
	outcomes = Realizations['outcomes']
	outcomes_list = list(itertools.product(*outcomes))
	S_number = range(1, len(outcomes_list)+1)

	uu = {}
	for s in S_number:
		count = 1
		for i in Uncertain['I_uncertain']:
			uu[i,s] = outcomes_list[s-1][count-1]
			count += 1

	D = {}
	for s in S_number:
		for sp in S_number:
			if s<sp:
				D[s,sp] = []
				for i in Uncertain['I_uncertain']:
					if uu[i,s] != uu[i,sp]:
						D[s,sp].append(i)
						
	L1 = {}
	for s in S_number:
		for sp in S_number:
			if s<sp:
				if len(D[s,sp]) == 1:
					L1[s,sp] = 1

	vv = {}
	for (s,sp) in L1:
		vv[s,sp] = D[s,sp][0]
			

	# generate scenario specific informations
	save_file = "scenario_specific"
	f = open(os.path.join(save_file),	"w")
	f.write('------------------------------------------uu ----------------------------------------------' + '\n')
	for i in Uncertain['I_uncertain']:
		f.write('----------------------------------' + '\n')
		for s in S_number:
			f.write('uu' + ' ' + str((i,s)) + ' ' + str(uu[i,s]) +  '\n')

	#############################################################################################
	#############################################################################################
	# concrete model
	model = ConcreteModel()
	#############################################################################################
	#############################################################################################
	# sets
	model.I1 = sets['I1']
	model.I2 = Uncertain['I_uncertain']
	model.K126 = sets['K126']
	model.T = sets['T']
	model.K = sets['K']
	model.S = sets['S']
	model.Kstream = sets['Kstream']
	#############################################################################################
	#############################################################################################
	# Parameters
	model.FE = parameters['FE']
	model.d = parameters['d'] 
	model.VE = parameters['VE']
	model.FO = parameters['FO']
	model.FIPP = parameters['FIPP']
	model.UQE = parameters['UQE'] 
	model.LQE = parameters['LQE']
	model.alpha = parameters['alpha']
	model.beta = parameters['beta']
	model.gamma = parameters['gamma']
	model.VO = parameters['VO']
	model.theta3 = parameters['theta3']
	model.T2 = parameters['T2']
	model.M = parameters['Big_M']
	model.Winital = parameters['Winital']
	model.theta = uu
	model.phi = vv
	model.PH = L1
	
	# Generate Parameters
	pro ={}
	if make['pro'] == {}:
		for s in S_number:
			pro[s] = 1/len(S_number)
	else:
		pro = make['pro']
		
	model.pro = pro
	#############################################################################################
	#############################################################################################
	# Variables
	model.Xpurch = Var(model.T, model.S, within =NonNegativeReals)
	model.Xsales = Var(model.T, model.S, within =NonNegativeReals)
	model.Wcap = Var(model.I1, model.T, model.S, within =NonNegativeReals)
	model.WQE = Var(model.I1, model.T, model.S, within =NonNegativeReals)
	#guan jian
	model.Wrate = Var(model.K, model.T, model.S, within =NonNegativeReals)
	model.Winv = Var(model.T, model.S, within =NonNegativeReals)
	model.NPV = Var(model.S)
	### Binary variable
	model.Yexp = Var(model.I1, model.T, model.S, within = Binary)
	model.Yoper = Var(model.I1, model.T, model.S, within = Binary)
	model.Z = Var(model.T, model.S, model.S, within = Binary)
	#############################################################################################
	# objective function
	#############################################################################################
	def objectivemodel(model):
		return sum(model.pro[s] * (model.NPV[s]) for s in model.S)
	model.objective = Objective(sense=maximize, rule=objectivemodel)

	def NPV_value(model,s):
		return model.NPV[s] == -sum(model.FE[i]*model.Yexp[i,t,s] + model.VE[i]*model.WQE[i,t,s] + model.FO[i]*model.Yoper[i,t,s] for i in model.I1 for t in model.T)  - sum(model.VO[k]*model.Wrate[k,t,s] for k in model.K126 for t in model.T) - sum(model.alpha*model.Xpurch[t,s] - model.beta*model.Xsales[t,s] for t in model.T) - sum(model.gamma*model.Winv[t,s] for t in model.T)
	model.NPV_value_constraints = Constraint(model.S, rule=NPV_value)
	# sum(0.001*model.Wrate[k,t,s] for k in model.Kstream for t in model.T)

	def B3(model,s,t):
		return model.Wrate[3,t,s] == model.theta[1,s]* model.Wrate[1,t,s]
	model.B3_constriants = Constraint(model.S, model.T, rule = B3)

	def B5(model,s,t):
		return model.Wrate[4,t,s] == model.theta[2,s]* model.Wrate[2,t,s]
	model.B5_constriants = Constraint(model.S, model.T, rule = B5)

	#######################################################
	def B14(model,s,t):
		return model.Wrate[8,t,s] == model.theta3*model.Wrate[7,t,s]
	model.B14_constriants = Constraint(model.S, model.T, rule = B14)

	def B15(model,s,t):
		return model.Wrate[5,t,s] == model.Wrate[3,t,s] + model.Wrate[4,t,s]
	model.B15_constriants = Constraint(model.S, model.T, rule = B15)

	def B16(model,s,t):
		return model.Wrate[7,t,s] == model.Wrate[5,t,s] + model.Wrate[6,t,s]
	model.B16_constriants = Constraint(model.S, model.T, rule = B16)

	#######################################################
	def B17(model,s,t):
		if t>= 2:
			return model.Winv[t,s] == model.Winv[t-1,s] + model.Wrate[8,t,s] + model.Xpurch[t,s] - model.Xsales[t,s]
		else:
			return model.Winv[t,s] == model.Wrate[8,t,s] + model.Xpurch[t,s] - model.Xsales[t,s]
	model.B17_constriants = Constraint(model.S, model.T, rule = B17)

	def B18(model,s,t):
		return model.Xsales[t,s] == model.d[t]
	model.B18_constriants = Constraint(model.S, model.T, rule = B18)

	#######################################################
	def B19(model,s,t):
		return model.Wrate[3,t,s] <= model.Wcap[1,t,s]
	model.B19_constriants = Constraint(model.S, model.T, rule = B19)

	def B20(model,s,t):
		return model.Wrate[4,t,s] <= model.Wcap[2,t,s]
	model.B20_constriants = Constraint(model.S, model.T, rule = B20)

	def B21(model,s,t):
		return model.Wrate[8,t,s] <= model.Wcap[3,t,s]
	model.B21_constriants = Constraint(model.S, model.T, rule = B21)

	#######################################################
	def B22(model,i,s,t):
		if t>=2:
			return model.Wcap[i,t,s] == model.Wcap[i,t-1,s] + model.WQE[i,t,s]
		else:
			return model.Wcap[i,t,s] == model.Winital[i] + model.WQE[i,t,s]
	model.B22_constriants = Constraint(model.I1, model.S, model.T, rule = B22)

	def B23a(model,i,s,t):
		return model.Yexp[i,t,s] <= model.WQE[i,t,s]
	model.B23a_constriants = Constraint(model.I1, model.S, model.T, rule = B23a)

	def B23b(model,i,s,t):
		return model.Yexp[i,t,s]*10 >= model.WQE[i,t,s]
	model.B23b_constriants = Constraint(model.I1, model.S, model.T, rule = B23b)

	def B24a(model,s,t):
		return 0.5 * model.Yoper[1,t,s] <= model.Wrate[3,t,s]
	model.B24a_constriants = Constraint(model.S, model.T, rule = B24a)

	def B24b(model,s,t):
		return model.Yoper[1,t,s]*10 >= model.Wrate[3,t,s]
	model.B24b_constriants = Constraint(model.S, model.T, rule = B24b)

	def B25a(model,s,t):
		return 0.5 * model.Yoper[2,t,s] <= model.Wrate[4,t,s]
	model.B25a_constriants = Constraint(model.S, model.T, rule = B25a)

	def B25b(model,s,t):
		return model.Yoper[2,t,s]*10 >= model.Wrate[4,t,s]
	model.B25b_constriants = Constraint(model.S, model.T, rule = B25b)

	def B26a(model,s,t):
		return 0.5 * model.Yoper[3,t,s] <= model.Wrate[8,t,s]
	model.B26a_constriants = Constraint(model.S, model.T, rule = B26a)

	def B26b(model,s,t):
		return model.Yoper[3,t,s]*10 >= model.Wrate[8,t,s]
	model.B26b_constriants = Constraint(model.S, model.T, rule = B26b)


	########################################################
	def B27(model,i,s,t):
		return sum(model.Yexp[i,tt,s]  for tt in model.T if tt<=t) >= model.Yoper[i,t,s]
	model.B27_constriants = Constraint(model.I2, model.S, model.T, rule = B27)

	def B28(model,i,s,t):
		return model.Yoper[i,t,s] >= model.Yexp[i,t,s]
	model.B28_constriants = Constraint(model.I1, model.S, model.T, rule = B28)

	####################################################
	########### Initial Non-anticipativity constraints
	#####################################################
	
	def B29(model,i,s,sp):
		return model.Yoper[i,1,s] == model.Yoper[i,1,sp]
	model.B29_constriants = Constraint(model.I1, model.S, model.S, rule = B29)


	def B31(model,i,s,sp):
		return model.Yexp[i,1,s] == model.Yexp[i,1,sp]
	model.B31_constriants = Constraint(model.I1, model.S, model.S, rule = B31)

	def B32(model,i,s,sp):
		return model.WQE[i,1,s] == model.WQE[i,1,sp]
	model.B32_constriants = Constraint(model.I1, model.S, model.S, rule = B32)

	def B33(model,k,s,sp):
		return model.Wrate[k,1,s] == model.Wrate[k,1,sp]
	model.B33_constriants = Constraint(model.K126, model.S, model.S, rule = B33)

	####################################################

	def B39(model,s,sp,t,tao):
		if tao <= t and (s,sp) in model.PH:
			return model.Z[t,s,sp]  <= 1 - model.Yoper[model.phi[s,sp],tao,s]
		else:
			return Constraint.Skip
	model.B39_constriants = Constraint(model.S, model.S, model.T,model.T, rule = B39)

	def B40(model,s,sp,t):
		if (s,sp) in model.PH:
			return model.Z[t,s,sp]  >= 1 - sum(model.Yoper[model.phi[s,sp],tao,s] for tao in model.T if tao<=t)
		else:
			return Constraint.Skip
	model.B40_constriants = Constraint(model.S, model.S, model.T, rule = B40)

	#################################################################
	########### Conditional Non-anticipativity constraints
	#################################################################

	def B41(model,s,sp,t):
		return model.Xpurch[t,s] <= model.Xpurch[t,sp] + model.M*(1 - model.Z[t,s,sp])
	model.B41_constriants = Constraint(model.S, model.S, model.T, rule = B41)

	def B42(model,s,sp,t):
		return model.Xpurch[t,s] >= model.Xpurch[t,sp] - model.M*(1 - model.Z[t,s,sp])
	model.B42_constriants = Constraint(model.S, model.S, model.T, rule = B42)

	def B43(model,s,sp,t):
		return model.Xsales[t,s] <= model.Xsales[t,sp] + model.M*(1 - model.Z[t,s,sp])
	model.B43_constriants = Constraint(model.S, model.S, model.T, rule = B43)

	def B44(model,s,sp,t):
		return model.Xsales[t,s] >= model.Xsales[t,sp] - model.M*(1 - model.Z[t,s,sp])
	model.B44_constriants = Constraint(model.S, model.S, model.T, rule = B44)

	def B45(model,i,s,sp,t):
		return model.Yoper[i,t+1,s] <= model.Yoper[i,t+1,sp] + (1 - model.Z[t,s,sp])
	model.B45_constriants = Constraint(model.I1, model.S, model.S, model.T2, rule = B45)

	def B46(model,i,s,sp,t):
		return model.Yoper[i,t+1,s] >= model.Yoper[i,t+1,sp] - (1 - model.Z[t,s,sp])
	model.B46_constriants = Constraint(model.I1, model.S, model.S, model.T2, rule = B46)

	def B47(model,i,s,sp,t):
		return model.Yexp[i,t+1,s] <= model.Yexp[i,t+1,sp] + (1 - model.Z[t,s,sp])
	model.B47_constriants = Constraint(model.I1, model.S, model.S, model.T2, rule = B47)

	def B48(model,i,s,sp,t):
		return model.Yexp[i,t+1,s] >= model.Yexp[i,t+1,sp] - (1 - model.Z[t,s,sp])
	model.B48_constriants = Constraint(model.I1, model.S, model.S, model.T2, rule = B48)


	def B51(model,i,s,sp,t):
		return model.WQE[i,t+1,s] <= model.WQE[i,t+1,sp] + model.M*(1 - model.Z[t,s,sp])
	model.B51_constriants = Constraint(model.I1, model.S, model.S, model.T2, rule = B51)

	def B52(model,i,s,sp,t):
		return model.WQE[i,t+1,s] >= model.WQE[i,t+1,sp] - model.M*(1 - model.Z[t,s,sp])
	model.B52_constriants = Constraint(model.I1, model.S, model.S, model.T2, rule = B52)

	def B53(model,k,s,sp,t):
		return model.Wrate[k,t+1,s] <= model.Wrate[k,t+1,sp] + model.M*(1 - model.Z[t,s,sp])
	model.B53_constriants = Constraint(model.K126, model.S, model.S, model.T2, rule = B53)

	def B54(model,k,s,sp,t):
		return model.Wrate[k,t+1,s] >= model.Wrate[k,t+1,sp] - model.M*(1 - model.Z[t,s,sp])
	model.B54_constriants = Constraint(model.K126, model.S, model.S, model.T2, rule = B54)

	results= opt.solve(model,keepfiles=True)
	model.solutions.load_from(results)	
	save_file1 = "Network_synthesis"

	save_file_binary = "Network_synthesis_variables" 
	fb = open(os.path.join(save_file_binary),	"w")

	### Generate New File Name
	results.write(filename=save_file1)
	fb.write('------------------------------------ Binary Variables---------------------------------------' + '\n')
	for i in model.I1:
		for t in model.T:
			for s in model.S:
				if model.Yoper[i,t,s].value != 0:
					fb.write('y_oper' + ' ' + str((i,t,s)) + ' ' + str(model.Yoper[i,t,s].value)  + '\n')
			
	for i in model.I1:
		for t in model.T:
			for s in model.S:
				if model.Yexp[i,t,s].value != 0:
					fb.write('y_exp' + ' ' + str((i,t,s)) + ' ' + str(model.Yexp[i,t,s].value)  + '\n')

	for s in model.S:
		for sp in model.S:
			for t in model.T:
				fb.write('Z' + ' ' + str((t,s,sp)) + ' ' + str(model.Z[t,s,sp].value)  + '\n')

	fb.write('------------------------------------ Continuous Variables---------------------------------------' + '\n')
	for t in model.T:
		for s in model.S:
			fb.write('x_purchase' + ' ' + str((t,s)) + ' ' + str(model.Xpurch[t,s].value)  + '\n')
						
	for t in model.T:
		for s in model.S:
			fb.write('x_sale' + ' ' + str((t,s)) + ' ' + str(model.Xsales[t,s].value)  + '\n')

	for i in model.I1:
		for t in model.T:
			for s in model.S:
				fb.write('W_cap' + ' ' + str((i,t,s)) + ' ' + str(model.Wcap[i,t,s].value)  + '\n')
				
	for i in model.I1:
		for t in model.T:
			for s in model.S:
				fb.write('W_QE' + ' ' + str((i,t,s)) + ' ' + str(model.WQE[i,t,s].value)  + '\n')
				
	for i in model.K:
		for t in model.T:
			for s in model.S:
				fb.write('W_rate' + ' ' + str((i,t,s)) + ' ' + str(model.Wrate[i,t,s].value)  + '\n')
				
	for t in model.T:
		for s in model.S:
			fb.write('W_in' + ' ' + str((t,s)) + ' ' + str(model.Winv[t,s].value)  + '\n')
			
	return