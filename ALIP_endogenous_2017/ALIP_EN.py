import os
import sys
from pyomo.environ import *
from pyomo.opt import SolverFactory
import pdb
import numpy
import itertools
import datetime

today = datetime.date.today()


def de(parameters, sets):
	
	T = sets['T']
	time_planning = T[-1]
	TT = range(1,time_planning+1)
	
	outcomes = [[1.176, 0.784],[1.26, 0.84], [1.296, 0.864]]
	outcomes_list = list(itertools.product(*outcomes))

	S_numbers = range(0, len(outcomes_list))


	v = {}
	vv = {}
	### define the endogenous uncertainty values:
	for s in S_numbers:
		v[5,s] = outcomes_list[s][0]
		v[6,s] = outcomes_list[s][1]
		v[7,s] = outcomes_list[s][2]
		
	### define the distinguishble scenario sets for endogenous uncertainties:
	for s in S_numbers:
		for sp in S_numbers:
			if sp>s:
				if v[7,s] == v[7,sp]:
					vv[7,s,sp]=1
				if v[6,s] == v[6,sp]:
					vv[6,s,sp]=1
				if v[5,s] == v[5,sp]:
					vv[5,s,sp]=1



	opt = SolverFactory("cplex")

	#############################################################################################
	#############################################################################################
	# concrete model
	model = ConcreteModel()
	#############################################################################################
	#############################################################################################
	# sets
	model.I = sets['I']


	model.T = TT
	model.S = S_numbers
	model.phi = vv
	model.endo_set = [i for i in model.phi]

	#############################################################################################
	#############################################################################################
	# Parameters
	pp ={}

	for s in S_numbers:
		pp[s] = 1/len(S_numbers)
		

	#############################################################################################
	#############################################################################################
	# Parameters

	model.pro = Param(model.S, initialize = pp)

	model.TubingSize = parameters['TubingSize']
	
	model.Packer = parameters['Packer']
	model.CICHP = parameters['CICHP']
	model.CITHP = parameters['CITHP']
	model.FLP = parameters['FLP']
	model.FTHP = parameters['FTHP']
	model.LF = parameters['LF']
	model.WellDeviation = parameters['WellDeviation']
	model.DLS = parameters['DLS']
	model.PTD = parameters['PTD']
	model.TubingUniformity = parameters['TubingUniformity']
	model.TAC = parameters['TAC']
	model.MCF = parameters['MCF']
	model.WellboreAccess = parameters['WellboreAccess']
	model.WellDepth = parameters['WellDepth']
	model.Temp = parameters['Temp']
	model.Offshore = parameters['Offshore']
	model.WaterCut = parameters['WaterCut']
	model.TDS = parameters['TDS']
	model.BHP = parameters['BHP']
	model.CGA = parameters['CGA']
	model.GVF = parameters['GVF']
	model.CasingDiameter = parameters['CasingDiameter']
	model.EPA = parameters['EPA']
	model.EPEF = parameters['EPEF']
	model.GPA = parameters['GPA']
	model.GPEF = parameters['GPEF']
	model.AromaticContent = parameters['AromaticContent']
	model.pg = parameters['pg']
	model.po = parameters['po']
	model.png = parameters['png']


	model.Qgbar = Param(model.I, model.S,initialize=v, default=0)

		### New Parameters
	model.WI = parameters['WI']
	model.M = parameters['M']
	model.MROR = parameters['MROR']
	model.FT = parameters['FT']
		### End New

	model.Cm = parameters['Cm']
	model.Co = parameters['Co']
	model.Ce = parameters['Ce']

	model.b = parameters['b']
	model.D = parameters['D']
	model.LocTax = parameters['LocTax']
	model.Royalty = parameters['Royalty']
	model.DownholeSeparator = parameters['DownholeSeparator']
		#############################################################################################
		#############################################################################################
		
		# Variables
	model.y = Var(model.I, model.T, model.T, model.S, within=Binary)
		# start
	model.W = Var(model.I, model.T, model.S, within=Binary)
		# end
	model.Z = Var(model.I, model.T, model.S, within=Binary)

		
	model.Qg = Var(model.T, model.S, within=NonNegativeReals)
	model.Qo = Var(model.T, model.S, within=NonNegativeReals)
	model.Qng = Var(model.T, model.S, within=NonNegativeReals)
	model.LFR = Var(model.T, model.S, within=NonNegativeReals)
	model.LFRy = Var(model.I, model.T, model.T, model.T, model.S,within=NonNegativeReals)
	model.LFRy1 = Var(model.I, model.T, model.T, model.T, model.S,within=NonNegativeReals)
	model.yQgest = Var(model.I, model.T, model.T, model.T, model.S,within=NonNegativeReals)
	model.yQgest1 = Var(model.I, model.T, model.T, model.T, model.S,within=NonNegativeReals)
	model.yQngest = Var(model.I, model.T, model.T, model.T, model.S,within=NonNegativeReals)
	model.yQngest1 = Var(model.I, model.T, model.T, model.T, model.S, within=NonNegativeReals)
	model.yQoest = Var(model.I, model.T, model.T, model.T, model.S, within=NonNegativeReals)
	model.yQoest1 = Var(model.I, model.T, model.T, model.T, model.S, within=NonNegativeReals)

	model.Qgest = Var(model.I, model.T, model.T, model.S, within=NonNegativeReals)
	model.Qngest = Var(model.I, model.T, model.T, model.S, within=NonNegativeReals)
	model.Qoest = Var(model.I, model.T, model.T, model.S, within=NonNegativeReals)
	### New Variables
	model.GIr = Var(model.T, model.S, within = NonNegativeReals)
	model.Revr = Var(model.T, model.S, within = NonNegativeReals)
	model.CC = Var(model.S, within = NonNegativeReals)
	model.TIr = Var(model.T, model.S, within = NonNegativeReals)
	model.Depr = Var(model.T, model.S, within = NonNegativeReals)
	model.XTIr = Var(model.T, model.S, within = NonNegativeReals)
	model.XTI1r = Var(model.T, model.S, within = NegativeReals)
	model.X = Var(model.T, model.S, within = Binary)
		### End New
		
		
		
		#############################################################################################
		# objective function
		#############################################################################################
		
	def objectivemodel(model):
		return sum(model.pro[s] * (sum((model.GIr[r,s] - model.XTIr[r,s] * model.FT) for r in model.T) - model.CC[s]*(1 - model.FT)) for s in model.S)
	model.objective = Objective(sense=maximize, rule=objectivemodel)

	   
	### Old Objective model.Qg[r]*180+model.Qo[r]*3000+model.Qng[r]*1500)*1/((1+model.MROR)**r)*(1-model.LocTax)*(1-model.Royalty)- sum(model.Cm[i] * model.y[i,t,p]*1/((1+model.MROR)**r) for i in model.I for t in model.T for p in model.T if t>=r and p<=r) - sum(model.Co[i]*model.y[i,t,p]*1/((1+model.MROR)**r) for i in model.I for t in model.T for p in model.T)for r in model.T)
	#############################################################################################
	#(model.Qg[r]*6+model.Qo[r]*100+model.Qng[r]*50)*1/((1+model.MROR)**r)*(1-model.LocTax)*(1-model.Royalty) for r in model.T
	#############################################################################################
	# - sum(model.Cm[i] * model.y[i,t,p] for i in model.I for t in model.T for p in model.T if t>=r and p<=r)*1/((1+model.MROR)**r) - sum(model.Co[i]*model.y[i,t,p] for i in model.I for t in model.T for p in model.T)*1/((1+model.MROR)**r)
	#initial value

	#((30*Qg(r,s)*Pg + 30*Qo(r,s)*Po + 30*Qng(r,s)*Png)*(1-Royalty)*(1-LocTax)*workinginterest
	#                    -sum((i,t,p)$(ord(t)>=ord(r) and ord(p)<=ord(r)),Cm(i)*w(i,t,p,s))*workinginterest)*1/(1+MROR)**rr(r);
		### New constraints
	def gross_income(model,r,s):
		return model.GIr[r,s]==(model.Revr[r,s]*(1 - model.Royalty)*(1 - model.LocTax)*model.WI - sum(model.Cm[i] * model.y[i,t,p,s] for i in model.I for t in model.T for p in model.T if t>=r and p<=r)*model.WI)*1 / ((1+model.MROR)**r)
	model.gross = Constraint(model.T, model.S, rule = gross_income)
		
	def revenue(model,r,s):
		return model.Revr[r,s] ==(30*model.Qg[r,s]*model.pg + 30*model.po*model.Qo[r,s] + 30*model.png* model.Qng[r,s])
	model.rev = Constraint(model.T, model.S, rule = revenue)
		
	def taxable_income(model,r,s):
		return model.GIr[r,s] - model.Depr[r,s] == model.TIr[r,s]
	model.taxable = Constraint(model.T, model.S, rule = taxable_income)
		
	def taxable_income1(model,r,s):
		return(model.X[r,s] * model.M) >= model.TIr[r,s]
	model.taxable1 = Constraint(model.T, model.S, rule = taxable_income1)
		
	def taxable_income2(model,r,s):
		return((model.X[r,s] - 1) * model.M) <= model.TIr[r,s]
	model.taxable2 = Constraint(model.T, model.S, rule = taxable_income2)
		
	def depreciation_function(model,r,s):
		return model.Depr[r,s] == sum(model.Ce[i]/5 * model.y[i,t,p,s] for i in model.I for t in model.T for p in model.T if t>=r and p+4>=r and p<=r and r<=t)*1/((1+model.MROR)**r)
	model.depreciation = Constraint(model.T, model.S, rule = depreciation_function)
		
	#sum((i,t,p)$(ord(p)+4>=ord(r) and ord(p)<=ord(r) and ord(t)>= ord(r)),Ce(i)/5*w(i,t,p,s))*1/(1+MROR)**rr(r); 
	 
	def capital_cost(model,s):
		return sum(model.Co[i] * model.y[i,t,p,s]/((1 + model.MROR)**p) for i in model.I for t in model.T for p in model.T)*model.WI == model.CC[s]
	model.capital = Constraint(model.S, rule = capital_cost)

	#Capitalcost(s) =e= sum((i,t,p),Co(i)*w(i,t,p,s)*1/(1+MROR)**pp(p))*workinginterest;  
	#

	def linearization_1(model,r,s):
		return(model.XTIr[r,s] + model.XTI1r[r,s]) == model.TIr[r,s]
	model.linear_1 = Constraint(model.T, model.S, rule = linearization_1)
		
	def linearization_2(model,r,s):
		return(model.X[r,s] * model.M) >= model.XTIr[r,s]
	model.linear_2 = Constraint(model.T, model.S, rule = linearization_2)

	def linearization_3(model,r,s):
		return(-(1 - model.X[r,s]) * model.M) <= model.XTI1r[r,s]
	model.linear_3 = Constraint(model.T, model.S, rule = linearization_3)


		### End New

	def initial_value_Qg(model,s):
		return model.Qg[1,s] - 800 == 0
	model.initial_Qg = Constraint(model.S,rule=initial_value_Qg)

	def initial_value_Qo(model,s):
		return model.Qo[1,s] - 20 == 0
	model.initial_Qo = Constraint(model.S,rule=initial_value_Qo)

	def initial_value_Qng(model,s):
		return model.Qng[1,s] - 500 == 0
	model.initial_Qng = Constraint(model.S,rule=initial_value_Qng)

	def initial_value_LFR(model,s):
		return model.LFR[1,s] - 520 == 0
	model.initial_LFR = Constraint(model.S,rule=initial_value_LFR)

	#logic constraints
	def logic_constraints1(model,i,s):
		return sum(model.y[i,t,p,s] for t in model.T for p in model.T) <= 1
	model.logic_constraints11 = Constraint(model.I, model.S, rule=logic_constraints1)

	def logic_constraints2(model,t,s):
		return sum(model.y[i,t,p,s] for i in model.I for p in model.T) <= 1
	model.logic_constraints22 = Constraint(model.T, model.S, rule=logic_constraints2)

	def logic_constraints3(model,p,s):
		return sum(model.y[i,t,p,s] for i in model.I for t in model.T) <= 1
	model.logic_constraints33 = Constraint(model.T, model.S, rule=logic_constraints3)


	def logic_constraints4(model,i,t,p,s):
		if t <= p+3:
			return model.y[i,t,p,s] == 0
		else:
			return Constraint.Skip
	model.logic_constraints44 = Constraint(model.I, model.T, model.T, model.S, rule=logic_constraints4)

	def logic_constraints5(model,i,t,p,j,k,l,s):
		if i!= j and p+1 <= l and l<= t:
			return model.y[i,t,p,s] + model.y[j,k,l,s] <= 1
		else:
			return Constraint.Skip
	model.logic_constraints55 = Constraint(model.I, model.T, model.T,model.I, model.T, model.T, model.S, rule=logic_constraints5)

	def logic_constraints6(model,r,s):
		return model.LFR[r,s] - model.Qo[r,s] - model.Qng[r,s] == 0
	model.logic_constraints66 = Constraint(model.T, model.S, rule=logic_constraints6)
		
		
	def logic_constraints7(model,i,s):
		return sum(model.Z[i,t,s] for t in model.T) == sum(model.W[i,p,s] for p in model.T)
	model.logic_constraints77 = Constraint(model.I,model.S,rule=logic_constraints7)
		

	#############################################################################################
	#############################################################################################
	# decline curve model
	#Qg
	def decline_curve_constraints1(model,r,s):
		if r>=2:
			return model.Qg[r,s] - sum(model.yQgest[i,p,r,t,s] for i in model.I for t in model.T for p in model.T if t>=r and r>=p and p>=2) == 0
		else:
			return Constraint.Skip
	model.decline_curves1 = Constraint(model.T, model.S, rule=decline_curve_constraints1)

	def decline_curve_constraints2(model,i,t,p,r,s):
		if t>=r and r>=p and p>=2:
			return model.yQgest[i,p,r,t,s] - 13070* 1.296 *model.y[i,t,p,s] <= 0
		else:
			return Constraint.Skip
	model.decline_curves2 = Constraint(model.I, model.T, model.T, model.T, model.S, rule=decline_curve_constraints2)

	def decline_curve_constraints3(model,i,t,p,r,s):
		if t>=r and r>=p and p>=2:
			return model.yQgest1[i,p,r,t,s] - 15370* 1.296 *(1 - model.y[i,t,p,s]) <= 0
		else:
			return Constraint.Skip
	model.decline_curves3 = Constraint(model.I, model.T, model.T, model.T, model.S, rule=decline_curve_constraints3)

	def decline_curve_constraints4(model,i,p,r,s):
		if r>=p and p>=2:
			return model.Qgest[i,p,r,s] - model.Qg[p-1,s] * model.Qgbar[i,s] * (1 + 0.6357*0.0507*(r-p+1))**(-1/0.6357) == 0
		else:
			return Constraint.Skip
	model.decline_curves4 = Constraint(model.I, model.T, model.T,  model.S, rule=decline_curve_constraints4)

	def decline_curve_constraints5(model,i,t,p,r,s):
		if t>=r and r>=p and p>=2:
			return model.yQgest[i,p,r,t,s] + model.yQgest1[i,p,r,t,s] - model.Qgest[i,p,r,s] == 0
		else:
			return Constraint.Skip
	model.decline_curves5 = Constraint(model.I, model.T, model.T, model.T, model.S, rule=decline_curve_constraints5)

	############################################ Qo
	def decline_curve_constraints6(model,r,s):
		if r>=2:
			return model.Qo[r,s] - sum(model.yQoest[i,p,r,t,s] for i in model.I for t in model.T for p in model.T if t>=r and r>=p and p>=2) == 0
		else:
			return Constraint.Skip
	model.decline_curves6 = Constraint(model.T,  model.S, rule=decline_curve_constraints6)

	def decline_curve_constraints7(model,i,t,p,r,s):
		if t>=r and r>=p and p>=2:
			return model.yQoest[i,p,r,t,s] - 330* 1.296 *model.y[i,t,p,s] <= 0
		else:
			return Constraint.Skip
	model.decline_curves7 = Constraint(model.I, model.T, model.T, model.T, model.S, rule=decline_curve_constraints7)

	def decline_curve_constraints8(model,i,t,p,r,s):
		if t>=r and r>=p and p>=2:
			return model.yQoest1[i,p,r,t,s] - 390* 1.296 *(1 - model.y[i,t,p,s]) <= 0
		else:
			return Constraint.Skip
	model.decline_curves8 = Constraint(model.I, model.T, model.T, model.T, model.S, rule=decline_curve_constraints8)

	def decline_curve_constraints9(model,i,p,r,s):
		if r>=p and p>=2:
			return model.Qoest[i,p,r,s] - model.Qo[p-1,s] * model.Qgbar[i,s] * (1 + 0.6357*0.0507*(r-p+1))**(-1/0.6357) == 0
		else:
			return Constraint.Skip
	model.decline_curves9 = Constraint(model.I, model.T, model.T, model.S, rule=decline_curve_constraints9)

	def decline_curve_constraints10(model,i,t,p,r,s):
		if t>=r and r>=p and p>=2:
			return model.yQoest[i,p,r,t,s] + model.yQoest1[i,p,r,t,s] - model.Qoest[i,p,r,s] == 0
		else:
			return Constraint.Skip
	model.decline_curves10 = Constraint(model.I, model.T, model.T, model.T, model.S, rule=decline_curve_constraints10)

	############################################## Qng
	def decline_curve_constraints11(model,r,s):
		if r>=2:
			return model.Qng[r,s] - sum(model.yQngest[i,p,r,t,s] for i in model.I for t in model.T for p in model.T if t>=r and r>=p and p>=2) == 0
		else:
			return Constraint.Skip
	model.decline_curves11 = Constraint(model.T, model.S, rule=decline_curve_constraints11)

	def decline_curve_constraints12(model,i,t,p,r,s):
		if t>=r and r>=p and p>=2:
			return model.yQngest[i,p,r,t,s] - 8170* 1.296 *model.y[i,t,p,s] <= 0
		else:
			return Constraint.Skip	
	model.decline_curves12 = Constraint(model.I, model.T, model.T, model.T, model.S, rule=decline_curve_constraints12)

	def decline_curve_constraints13(model,i,t,p,r,s):
		if t>=r and r>=p and p>=2:
			return model.yQngest1[i,p,r,t,s] - 9610* 1.296 *(1 - model.y[i,t,p,s]) <= 0
		else:
			return Constraint.Skip	
	model.decline_curves13 = Constraint(model.I, model.T, model.T, model.T, model.S, rule=decline_curve_constraints13)

	def decline_curve_constraints14(model,i,p,r,s):
		if r>=p and p>=2:
			return model.Qngest[i,p,r,s] - model.Qng[p-1,s] * model.Qgbar[i,s] * (1 + 0.6357*0.0507*(r-p+1))**(-1/0.6357) == 0
		else:
			return Constraint.Skip	
	model.decline_curves14 = Constraint(model.I, model.T, model.T, model.S, rule=decline_curve_constraints14)

	def decline_curve_constraints15(model,i,t,p,r,s):
		if t>=r and r>=p and p>=2:
			return model.yQngest[i,p,r,t,s] + model.yQngest1[i,p,r,t,s] - model.Qngest[i,p,r,s] == 0
		else:
			return Constraint.Skip	
	model.decline_curves15 = Constraint(model.I, model.T, model.T, model.T, model.S, rule=decline_curve_constraints15)
		
		
		# Inequality constraints
	def Inequality_constraints_1(model,r,s):
		if r>=2:
			return model.LFR[r,s] <= 520 * 1.296 * 1.296 * 1.26 * (1 + 0.6357*0.05*(r-1))**(-1/0.6357)
		else:
			return Constraint.Skip	
	model.Inequality_1 = Constraint(model.T, model.S, rule=Inequality_constraints_1)
		
	def Inequality_constraints_2(model,r,s):
		if r>=2:
			return model.Qg[r,s] <= 800 * 1.296 * 1.296 * 1.26 * (1 + 0.6357*0.05*(r-1))**(-1/0.6357)
		else:
			return Constraint.Skip	
	model.Inequality_2 = Constraint(model.T, model.S, rule=Inequality_constraints_2)
		
	def Inequality_constraints_3(model,r,s):
		if r>=2:
			return model.Qo[r,s] <= 20 * 1.296 * 1.296 * 1.26 * (1 + 0.6357*0.05*(r-1))**(-1/0.6357)
		else:
			return Constraint.Skip	
	model.Inequality_3 = Constraint(model.T, model.S, rule=Inequality_constraints_3)
		
	def Inequality_constraints_4(model,r,s):
		if r>=2:
			return model.Qng[r,s] <= 500 * 1.296 * 1.296 * 1.26 * (1 + 0.6357*0.05*(r-1))**(-1/0.6357)
		else:
			return Constraint.Skip	
	model.Inequality_4 = Constraint(model.T, model.S, rule=Inequality_constraints_4)
		
		
		
		# Relationship
	def Net_Constraint_rule(model,i,t,p,s):
		return  model.y[i,t,p,s] >= model.W[i,p,s] + model.Z[i,t,s]-1
	model.Net_Constraint = Constraint(model.I, model.T,  model.T, model.S,rule=Net_Constraint_rule)
		
	def Net_Constraint1_rule(model,i,t,p,s):
		return  model.W[i,p,s] >= model.y[i,t,p,s]
	model.Net_Constraint1 = Constraint(model.I, model.T,  model.T, model.S,rule=Net_Constraint1_rule)

	def Net_Constraint2_rule(model,i,t,p,s):
		return  model.Z[i,t,s] >= model.y[i,t,p,s]
	model.Net_Constraint2 = Constraint(model.I, model.T,  model.T, model.S,rule=Net_Constraint2_rule)


		
		
		# initial NACs
	def initial_NACs(model,i,s):
		return model.W[i,2,0] == model.W[i,2,s]
	model.initial = Constraint(model.I, model.S, rule=initial_NACs)

		# exogenous NACs
	# print(model.phi)
		# method 5,6,7
	def NACs_1(model,i,j,t,p,s,sp):
		if sp > s and p>=2 and t >= p + 1 and j!=i and (j,s,sp) in vv:
			return model.W[i,t,s] - model.W[i,t,sp] <= 1 - model.Z[j,p,s]
		else:
			return Constraint.Skip
	model.N1 = Constraint(model.I, model.I, model.T, model.T, model.S,model.S, rule=NACs_1)
		
	def NACs_2(model,i,j,t,p,s,sp):
		if sp > s and p>=2 and t >= p + 1 and j!=i and (j,s,sp) in vv:
			return model.W[i,t,s] - model.W[i,t,sp] >= model.Z[j,p,s] - 1
		else:
			return Constraint.Skip
	model.N11 = Constraint(model.I, model.I, model.T, model.T,  model.S,model.S, rule=NACs_2)
		


	def physical_constraints(model,t,p,r,s):
		if r >= p and r<=t:
			return model.y[7,t,p,s]*300 - model.LFR[r,s] <=0
		else:
			return Constraint.Skip
	model.physical = Constraint(model.T, model.T, model.T, model.S, rule = physical_constraints)



	yy = {}
	results= opt.solve(model)
	model.solutions.load_from(results)	
	save_file1 = "ALIP" + '_' + str(TT[-1]) + ' ' + '(' + str(today) + ')'
	save_file = "ALIPOutput" + '_' + str(TT[-1]) + ' ' + '(' + str(today) + ')'

	dir_path = os.path.dirname(os.path.realpath(__file__))
	# TempfileManager.tempdir = dir_path

	f = open(os.path.join(save_file),	"w")
		### Generate New File Name
	for s in model.S:
		for i in model.I:
			for p in model.T:
				for t in model.T:
					if model.y[i,t,p,s].value >0.9: 
						f.write('y' + ' ' + str((i,t,p,s)) + ' ' + str(model.y[i,t,p,s].value) + '\n')
	f.write('------------------------------------------------------------------------'+'\n')
	for s in model.S:
		for i in model.I:
			for p in model.T:
				if model.W[i,p,s].value > 0.9:
					f.write('w' + ' ' + str((i,p,s)) + ' ' + str(model.W[i,p,s].value) + '\n')
	f.write('------------------------------------------------------------------------' + '\n')
	for s in model.S:
		for i in model.I:
			for p in model.T:
				if model.Z[i,p,s].value > 0.9:
					f.write('z' + ' ' + str((i,p,s)) + ' ' + str(model.Z[i,p,s].value) + '\n')
	results.write(filename=save_file1)
	f.close
	return
	