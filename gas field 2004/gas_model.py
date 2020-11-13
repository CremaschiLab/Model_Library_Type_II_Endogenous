import os
import sys
from pyomo.environ import *
from pyomo.opt import SolverFactory
import pdb
import numpy
import itertools
import datetime
import time

today = datetime.date.today()


def de(parameters, sets, Uncertain, Realizations):
	
	# generate solution file
	current_directory = os.path.dirname(os.path.realpath(__file__))
	current_date = time.strftime('%m_%d_%Y', time.gmtime())
	output_directory = current_directory + '/Solutions/' + current_date + '/'
	
	
	
	# sets
	WP = sets['WP']
	PP = sets['PP']
	T = sets['T']

	# Uncertain Parameters
	WP_uncertain = Uncertain['WP_uncertain']
	outcomes = Realizations['outcomes']
	outcomes_list = list(itertools.product(*outcomes))
	
	S = range(1, len(outcomes_list)+1)

	deliverbility = parameters['deliverbility']
	size = parameters['size']


	theta_2 = {}
	theta_1 = {}
	for wp in WP:
		for s in S:
			theta_2[wp,s] = deliverbility[wp]
				
	for s in S:
		for wp in WP:
			if wp == 'F':
				theta_1[wp,s] = outcomes_list[s-1][0]
			else:
				theta_1[wp,s] = size[wp]

	pro = parameters['probability']

	D = {}
	for wp in WP_uncertain:
		for s in S:
			for sp in S:
				if s<sp:
					D[s,sp] = []
					if theta_1[wp,s] != theta_1[wp,sp]:
						D[s,sp].append(wp)
	
	L1 = {}
	for s in S:
		for sp in S:
			if s<sp:
				if len(D[s,sp]) == 1:
					L1[s,sp] = 1
				
	# print('D',D)
	# print('L1',L1)
	### define the distinguishble scenario sets for endogenous uncertainties:

	opt = SolverFactory("cplex")
	#############################################################################################
	#############################################################################################
	# concrete model
	model = ConcreteModel()

	#############################################################################################
	#############################################################################################
	# Parameters
	delta = {}
	P = {}
	for t in T:
		delta[t] = 100
		P[t] = 100
		
	shrink = parameters['shrink']
	M = parameters['Big_M']
	
	FCC = {}
	VCC = {}
	FOC = {}
	VOC = {}
	for wp in WP:
		FCC[wp] = 1
		VCC[wp] = 1
		FOC[wp] = 1
		VOC[wp] = 1
	for pp in PP:
		FCC[pp] = 1
		VCC[pp] = 1
		FOC[pp] = 1
		VOC[pp] = 1
	
	FCC1 = {}
	for wp in WP:
		for wpp in WP:
			FCC1[wp,wpp] = 1
		for p in PP:
			FCC1[wp,p] = 1


	###########################################################################################
	model.b = Var(WP, T, S, within=Binary)
	model.b2 = Var(PP, T, S, within=Binary)
	model.b3 = Var(WP, WP, T, S, within=Binary)
	model.b4 = Var(WP, PP, T, S, within=Binary)
	model.Z = Var(T, S, S, within=Binary)
	
	model.q_out = Var(WP, T, S, within=NonNegativeReals)
	model.q2_out = Var(PP, T, S, within=NonNegativeReals)
	model.q3_out = Var(WP, WP, T, S, within=NonNegativeReals)
	model.q4_out = Var(WP, PP, T, S, within=NonNegativeReals)
	
	model.e1 = Var(WP, T, S, within=NonNegativeReals)
	model.e2 = Var(PP, T, S, within=NonNegativeReals)
	
	model.q_cum = Var(WP, T, S, within=NonNegativeReals)
	model.q_prod = Var(WP, T, S, within=NonNegativeReals)
	model.q_deliv = Var(WP, T, S, within=NonNegativeReals)
	
	model.cap1 = Var(WP, T, S, within=NonNegativeReals)
	model.cap2 = Var(PP, T, S, within=NonNegativeReals)
	
	model.q_shr = Var(T, S, within=NonNegativeReals)
	
	model.CC = Var(T, S, within=NonNegativeReals)
	model.QC = Var(T, S, within=NonNegativeReals)
	model.Rev = Var(T, S, within=NonNegativeReals)
	
	model.NPV = Var(S)
	#############################################################################################
	# objective function
	#############################################################################################
		
	def objectivemodel(model):
		return sum(pro[s]*model.NPV[s] for s in S)
	model.objective = Objective(sense=maximize, rule=objectivemodel)


	def con_2(model, wp, t, s):
		return model.q_deliv[wp,t,s]/theta_2[wp,s] + model.q_cum[wp,t,s]/theta_1[wp,s] == 1
	model.con_2m = Constraint(WP, T, S, rule = con_2)
		
	def con_3(model, wp, t, s):
		return model.q_cum[wp,t,s] == sum(model.q_prod[wp,t,s]*delta[tao] for tao in range(1,t+1))
	model.con_3m = Constraint(WP, T, S, rule = con_3)
	
	def con_4(model, wp, t, s):
		return model.q_prod[wp,t,s] <= model.q_deliv[wp,t,s]
	model.con_4m = Constraint(WP, T, S,  rule = con_4)
	
	def con_5(model, wp, t, s):
		return model.q_out[wp,t,s] == model.q_prod[wp,t,s] + sum(model.q3_out[wprime,wp,t,s] for wprime in WP if wprime != wp)
	model.con_5m = Constraint(WP, T, S, rule = con_5)

	def con_6(model, wp, t, s):
		return model.q_out[wp,t,s] == sum(model.q3_out[wprime,wp,t,s] for wprime in WP if wprime != wp) + sum(model.q4_out[wp,pp,t,s] for pp in PP)
	model.con_6m = Constraint(WP, T, S, rule = con_6)

	def con_7(model, pp, t, s):
		return model.q2_out[pp,t,s] == sum(model.q4_out[wp,pp,t,s] for wp in WP)
	model.con_7m = Constraint(PP, T, S, rule = con_7)
	
	def con_8(model, t, s):
		return model.q_shr[t,s] == (1-shrink)*sum(model.q2_out[pp,t,s] for pp in PP)
	model.con_8m = Constraint(T, S, rule = con_8)
	
	#############################################################################################
	def con_9(model,wp, t, s):
		return model.q_out[wp,t,s] <= model.cap1[wp,t,s]
	model.con_9m = Constraint(WP, T, S, rule = con_9)
	
	def con_10(model,pp, t, s):
		return model.q2_out[pp,t,s] <= model.cap2[pp,t,s]
	model.con_10m = Constraint(PP, T, S, rule = con_10)
	
	def con_11(model,wp, t, s):
		if t-1 >= 1:
			return model.cap1[wp,t,s] == model.cap1[wp,t-1,s] + model.e1[wp,t,s]
		else:
			return model.cap1[wp,t,s] == model.e1[wp,t,s]
	model.con_11m = Constraint(WP, T, S, rule = con_11)
	
	def con_12(model,pp, t, s):
		if t-1>=1:
			return model.cap2[pp,t,s] == model.cap2[pp,t-1,s] + model.e2[pp,t,s]
		else:
			return model.cap2[pp,t,s] == model.e2[pp,t,s]
	model.con_12m = Constraint(PP, T, S, rule = con_12)
	
	def con_13(model,wp, t, s):
		return model.e1[wp,t,s] <= M*model.b[wp,t,s]
	model.con_13m = Constraint(WP, T, S, rule = con_13)
	
	def con_14(model,pp, t, s):
		return model.e2[pp,t,s] <= M*model.b2[pp,t,s]
	model.con_14m = Constraint(PP, T, S, rule = con_14)
	
	def con_15(model,wp,wprime,t,s):
		return  model.q3_out[wp,wprime,t,s]<= M*sum(model.b3[wp,wprime,tao,s] for tao in range(1,t+1))
	model.con_15m = Constraint(WP, WP, T, S, rule = con_15)
	
	def con_16(model,wp,pp,t,s):
		return  model.q4_out[wp,pp,t,s]<= M*sum(model.b4[wp,pp,tao,s] for tao in range(1,t+1))
	model.con_16m = Constraint(WP, PP, T, S, rule = con_16)
	#############################################################################################
	def con_17(model,wp,s):
		return  sum(model.b[wp,t,s] for t in T) <= 1
	model.con_17m = Constraint(WP, S, rule = con_17)
	
	def con_18(model,pp,s):
		return  sum(model.b2[pp,t,s] for t in T) <= 1
	model.con_18m = Constraint(PP, S, rule = con_18)
	
	def con_19(model,wp,wprime,s):
		return  sum(model.b3[wp,wprime,t,s] for t in T) <= 1
	model.con_19m = Constraint(WP, WP, S, rule = con_19)
	
	def con_20(model,wp,pp,s):
		return  sum(model.b4[wp,pp,t,s] for t in T) <= 1
	model.con_20m = Constraint(WP, PP, S, rule = con_20)
	
	def con_21(model, wp, t, s):
		return model.b[wp,t,s] == sum(model.b3[wp,wprime,t,s] for wprime in WP if wprime != wp) + sum(model.b4[wp,pp,t,s] for pp in PP)
	model.con_21m = Constraint(WP, T, S, rule = con_21)
	
	def con_22(model, wp,wprime, t, s):
		return model.b3[wp,wprime,t,s] <= sum(model.b[wprime,tao,s] for tao in range(1,t+1))
	model.con_22m = Constraint(WP, WP, T, S, rule = con_22)
	
	def con_23(model, wp, pp, t, s):
		return model.b4[wp,pp,t,s] <= sum(model.b2[pp,tao,s] for tao in range(1,t+1))
	model.con_23m = Constraint(WP, PP, T, S, rule = con_23)
	
	def con_24(model, wp, wprime, t, s):
		return model.b3[wp,wprime,t,s] + model.b3[wprime,wp,t,s] <= 1
	model.con_24m = Constraint(WP, WP, T, S, rule = con_24)
	
	#############################################################################################
	
	def con_25(model, t, s):
		return model.CC[t,s] == sum(FCC[wp]*model.b[wp,t,s] + VCC[wp]*model.e1[wp,t,s] + \
										sum(FCC1[wp,wprime]*model.b3[wp,wprime,t,s] for wprime in WP) + \
										sum(FCC1[wp,pp]*model.b4[wp,pp,t,s] for pp in PP) for wp in WP) + \
								sum(FCC[pp]*model.b2[pp,t,s] + VCC[pp]*model.e2[pp,t,s] for pp in PP)
	model.con_25m = Constraint(T, S, rule = con_25)
	
	def con_26(model, t, s):
		return model.QC[t,s] == sum(FOC[wp]*model.b[wp,t,s] + VOC[wp]*model.q_prod[wp,t,s] for wp in WP) + \
								sum(FOC[pp]*model.b2[pp,t,s] + VOC[pp]*model.q2_out[pp,t,s] for pp in PP)
	model.con_26m = Constraint(T, S, rule = con_26)
	
	def con_27(model, t, s):
		return model.Rev[t,s] == sum(P[tao]*delta[tao]*model.q_shr[tao,s] for tao in T)
	model.con_27m = Constraint(T, S, rule = con_27)
	
	def con_28(model, s):
		return model.NPV[s] == sum(model.Rev[t,s] - model.CC[t,s] - model.QC[t,s] for t in T)
	model.con_28m = Constraint(S, rule = con_28)
	#############################################################################################
	# Initial NACs
	#############################################################################################
	def con_29a(model, wp,s,sprime):
		if s < sprime:
			return model.b[wp,1,s] == model.b[wp,1,sprime]
		else:
			return Constraint.Skip
	model.con_29am = Constraint(WP, S, S, rule = con_29a)
	
	def con_29b(model, wp,s,sprime):
		if s < sprime:
			return model.e1[wp,1,s] == model.e1[wp,1,sprime]
		else:
			return Constraint.Skip
	model.con_29bm = Constraint(WP, S, S, rule = con_29b)
	
	def con_29c(model, pp,s,sprime):
		if s < sprime:
			return model.b2[pp,1,s] == model.b2[pp,1,sprime]
		else:
			return Constraint.Skip
	model.con_29cm = Constraint(PP, S, S, rule = con_29c)
	
	def con_29d(model, pp,s,sprime):
		if s < sprime:
			return model.e2[pp,1,s] == model.e2[pp,1,sprime]
		else:
			return Constraint.Skip
	model.con_29dm = Constraint(PP, S, S, rule = con_29d)
	
	def con_29e(model, wp,wprime,s,sprime):
		if s < sprime:
			return model.b3[wp,wprime,1,s] == model.b3[wp,wprime,1,sprime]
		else:
			return Constraint.Skip
	model.con_29em = Constraint(WP, WP, S, S, rule = con_29e)
	
	def con_29f(model, wp,pp,s,sprime):
		if s < sprime:
			return model.b4[wp,pp,1,s] == model.b4[wp,pp,1,sprime]
		else:
			return Constraint.Skip
	model.con_29fm = Constraint(WP, PP, S, S, rule = con_29f)
	#############################################################################################
	# Equation 30
	#############################################################################################
	def con_30a(model, wp, t, tao,s, sprime):
		if s < sprime and (s,sprime) in L1 and tao<=t and wp in D[s,sprime]:
			return model.Z[t,s,sprime] <= 1 - model.b[wp,tao,s]
		else:
			return Constraint.Skip
	model.con_30am = Constraint(WP, T, T, S, S, rule = con_30a)
	
	def con_30b(model, wp, t, s, sprime):
		if s < sprime and (s,sprime) in L1 and wp in D[s,sprime]:
			return model.Z[t,s,sprime] >= 1 - sum(model.b[wp,tao,s] for tao in range(1,t+1))
		else:
			return Constraint.Skip
	model.con_30bm = Constraint(WP, T, S, S, rule = con_30b)
	#############################################################################################
	# Conditional NACs Equation 31
	#############################################################################################
	# b_wp_t_s
	def NAC1a(model, wp, t, s, sprime):
		if s < sprime and (s,sprime) in L1 and wp not in D[s,sprime] and t+1<=T[-1]:
			return model.b[wp,t+1,s] -  model.b[wp,t+1,sprime]<= 1 - model.Z[t,s,sprime]
		else:
			return Constraint.Skip
	model.NAC1am = Constraint(WP, T, S, S, rule = NAC1a)
	
	def NAC1b(model, wp, t, s, sprime):
		if s < sprime and (s,sprime) in L1 and wp not in [s,sprime] and t+1<=T[-1]:
			return model.b[wp,t+1,s] -  model.b[wp,t+1,sprime] >= -1 + model.Z[t,s,sprime]
		else:
			return Constraint.Skip
	model.NAC1bm = Constraint(WP, T, S, S, rule = NAC1b)
	
	# e_wp_t_s
	def NAC2a(model, wp, t, s, sprime):
		if s < sprime and (s,sprime) in L1 and wp not in D[s,sprime] and t+1<=T[-1]:
			return model.e1[wp,t+1,s] -  model.e1[wp,t+1,sprime]<= 1 - model.Z[t,s,sprime]
		else:
			return Constraint.Skip
	model.NAC2am = Constraint(WP, T, S, S, rule = NAC2a)
	
	def NAC2b(model, wp, t, s, sprime):
		if s < sprime and (s,sprime) in L1 and wp not in D[s,sprime] and t+1<=T[-1]:
			return model.e1[wp,t+1,s] -  model.e1[wp,t+1,sprime] >= -1 + model.Z[t,s,sprime]
		else:
			return Constraint.Skip
	model.NAC2bm = Constraint(WP, T, S, S, rule = NAC2b)
	
	# b_pp_t_s
	def NAC3a(model, pp, t, s, sprime):
		if s < sprime and (s,sprime) in L1 and t+1<=T[-1]:
			return model.b2[pp,t+1,s] -  model.b2[pp,t+1,sprime]<= 1 - model.Z[t,s,sprime]
		else:
			return Constraint.Skip
	model.NAC3am = Constraint(PP, T, S, S, rule = NAC3a)
	
	def NAC3b(model, pp, t, s, sprime):
		if s < sprime and (s,sprime) in L1 and t+1<=T[-1]:
			return model.b2[pp,t+1,s] -  model.b2[pp,t+1,sprime] >= -1 + model.Z[t,s,sprime]
		else:
			return Constraint.Skip
	model.NAC3bm = Constraint(PP, T, S, S, rule = NAC3b)
	
	# e_pp_t_s
	def NAC4a(model, pp, t, s, sprime):
		if s < sprime and (s,sprime) in L1 and t+1<=T[-1]:
			return model.e2[pp,t+1,s] -  model.e2[pp,t+1,sprime]<= 1 - model.Z[t,s,sprime]
		else:
			return Constraint.Skip
	model.NAC4am = Constraint(PP, T, S, S, rule = NAC4a)
	
	def NAC4b(model, pp, t, s, sprime):
		if s < sprime and (s,sprime) in L1 and t+1<=T[-1]:
			return model.e2[pp,t+1,s] -  model.e2[pp,t+1,sprime] >= -1 + model.Z[t,s,sprime]
		else:
			return Constraint.Skip
	model.NAC4bm = Constraint(PP, T, S, S, rule = NAC4b)
	
	# b_wp_wp'_t
	def NAC5a(model, wp, wprime, t, s, sprime):
		if s < sprime and (s,sprime) in L1 and t+1<=T[-1]:
			return model.b3[wp,wprime,t+1,s] -  model.b3[wp,wprime,t+1,sprime]<= 1 - model.Z[t,s,sprime]
		else:
			return Constraint.Skip
	model.NAC5am = Constraint(WP, WP, T, S, S, rule = NAC5a)
	
	def NAC5b(model, wp, wprime, t, s, sprime):
		if s < sprime and (s,sprime) in L1 and t+1<=T[-1]:
			return model.b3[wp,wprime,t+1,s] -  model.b3[wp,wprime,t+1,sprime] >= -1 + model.Z[t,s,sprime]
		else:
			return Constraint.Skip
	model.NAC5bm = Constraint(WP, WP, T, S, S, rule = NAC5b)
	
	# b_wp_pp_t
	def NAC6a(model, wp, pp, t, s, sprime):
		if s < sprime and (s,sprime) in L1 and t+1<=T[-1]:
			return model.b4[wp,pp,t+1,s] -  model.b4[wp,p,t+1,sprime]<= 1 - model.Z[t,s,sprime]
		else:
			return Constraint.Skip
	model.NAC6am = Constraint(WP, PP, T, S, S, rule = NAC6a)
	
	def NAC6b(model, wp, pp, t, s, sprime):
		if s < sprime and (s,sprime) in L1 and t+1<=T[-1]:
			return model.b4[wp,p,t+1,s] -  model.b4[wp,p,t+1,sprime] >= -1 + model.Z[t,s,sprime]
		else:
			return Constraint.Skip
	model.NAC6bm = Constraint(WP, PP, T, S, S, rule = NAC6b)
	
	# q_prod
	def NAC7a(model, wp, t, s, sprime):
		if s < sprime and (s,sprime) in L1 and t+1<=T[-1]:
			return model.q_prod[wp,t,s] -  model.q_prod[wp,t,sprime]<= 1 - model.Z[t,s,sprime]
		else:
			return Constraint.Skip
	model.NAC7am = Constraint(WP, T, S, S, rule = NAC7a)
	
	def NAC7b(model, wp, t, s, sprime):
		if s < sprime and (s,sprime) in L1 and t+1<=T[-1]:
			return model.q_prod[wp,t,s] -  model.q_prod[wp,t,sprime] >= -1 + model.Z[t,s,sprime]
		else:
			return Constraint.Skip
	model.NAC7bm = Constraint(WP, T, S, S, rule = NAC7b)


	results= opt.solve(model)
	model.solutions.load_from(results)
	
	### Open save file
	if not os.path.exists(output_directory):
		os.makedirs(output_directory)

	
	save_file1 ="Oil and gas-field development" + '_' + str(T[-1]) + ' ' + '(' + str(today) + ')'
	save_file = "Oil and gas-field development output" + '_' + str(T[-1]) + ' ' + '(' + str(today) + ')'
	
	dir_path = os.path.dirname(os.path.realpath(__file__))
	# TempfileManager.tempdir = dir_path

	f = open(os.path.join(output_directory, save_file),	"w")
	
	results.write(filename = os.path.join(output_directory, save_file1))
	
	for wp in WP:
		for t in T:
			for s in S:
				if model.b[wp,t,s].value >0.9: 
					f.write('b' + ' ' + str((wp,t,s)) + ' ' + str(model.b[wp,t,s].value) + '\n')


	return
	