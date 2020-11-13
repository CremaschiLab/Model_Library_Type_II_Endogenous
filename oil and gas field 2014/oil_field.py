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
	
	# sets
	S = sets['S']
	T = sets['T'] 
	RF = sets['RF']
	F = sets['F']
	FPSO = sets['FPSO']
	F_rf = sets['F_rf']
	F_fpso = sets['F_fpso']

	# uncertain parameters
	D = parameters['D']
	L1 = parameters['L1']
	
	opt = SolverFactory("cplex")
	#############################################################################################
	#############################################################################################
	# concrete model
	model = ConcreteModel()

	#############################################################################################
	#############################################################################################
	# Parameters
	dis = parameters['dis'] 
	F = parameters['F'] 
	FC_FPSO = parameters['FC_FPSO']
	VC_liq = parameters['VC_liq']
	VC_gas = parameters['VC_gas']
	FC = parameters['FC']
	FC_well = parameters['FC_well']
	REC = parameters['REC']
	delta = parameters['delta']
	OC_gas = parameters['OC_gas']
	f_tax = parameters['f_tax']
	f_po = parameters['f_po']
	U_oil = parameters['U_oil']
	L_oil = parameters['L_oil']
	alpha = parameters['alpha']
	f_CR = parameters['f_CR']
	L = parameters['L']
	alpha_oil = parameters['alpha_oil']
	alpha_gc = parameters['alpha_gc']
	alpha_wc = parameters['alpha_wc']
	a = parameters['a']
	b = parameters['b']
	c = parameters['c']
	d = parameters['d']
	M_wc = parameters['M_wc']
	M = parameters['M']
	U_well_oil = parameters['U_well_oil']
	UN_well = parameters['UN_well']
	N1 = parameters['N1']
	N2 = parameters['N2']
	l1 = parameters['l1']
	l2 = parameters['l2']
	###########################################################################################
	# Variable
	model.ENPV = Var()
	model.TConSh = Var(T, S, within = Reals)
	model.TCAP  = Var(T, S, within = Reals)
	model.TOPER = Var(T, S, within = Reals)
	model.ConSh = Var(RF, T, S, within = Reals)
	model.CAP  = Var(RF, T, S, within = Reals)
	model.OPER = Var(RF, T, S, within = Reals)
	model.CAP1  = Var(RF, T, S, within = Reals)
	model.CAP2 = Var(RF, T, S, within = Reals)
	model.DFPSOC = Var(RF, FPSO,  T, S, within = Reals)
	model.FPSOC = Var(FPSO,  T, S, within = Reals)
	model.QI_liq = Var(FPSO,  T, S, within = Reals)
	model.QE_liq = Var(FPSO,  T, S, within = Reals)
	model.QI_gas  = Var(FPSO,  T, S, within = Reals)
	model.QE_gas = Var(FPSO,  T, S, within = Reals)
	model.ConSh_aftertax = Var(RF, T, S, within = Reals)
	model.CO = Var(RF, T, S, within = Reals)
	model.DConSh_beforetax = Var(RF, I, T, S, within = Reals)
	model.Q_d_well = Var(F, FPSO, T, S, within = Reals)
	model.wc = Var(F, FPSO, T, S, within = Reals)
	model.Q_wc = Var(F, FPSO, T, S, within = Reals)
	model.x_c = Var(F, FPSO, T, S, within = Reals)
	model.x_well = Var(F, FPSO, T, S, within = Reals)
	model.N_well = Var(F, T, S, within = Reals)
	model.DFPSOC_field =  Var(F, FPSO, T, S, within = Reals)
	
	model.w_tot = Var(RF, T, S, within = Reals)
	model.x_tot = Var(RF, T, S, within = Reals)
	model.g_tot = Var(RF, T, S, within = Reals)
	model.Tax = Var(RF, T, S, within = Reals)
	model.ConSh_beforetax  = Var(RF, T, S, within = Reals)
	model.DPO = Var(RF, i, T, S, within = Reals)
	model.PO = Var(RF, T, S, within = Reals)
	model.xc = Var(RF, T, S, within = Reals)
	model.Dxc = Var(RF, i, T, S, within = Reals)
	model.xc_field[f,t,s] = Var(F, T, S, within = Reals)
	model.REV = Var(RF, T, S, within = Reals)
	model.x_f = Var(F, T, S, within = Reals)
	model.CR = Var(RF, T, S, within = Reals)
	model.CO = Var(RF, T, S, within = Reals)
	model.CRF = Var(RF, T, S, within = Reals)
	model.fc = Var(F, T, S, within = Reals)
	model.Q_gc = Var(F, FPSO, T, S, within = Reals)
	model.gc = Var(F, FPSO, T, S, within = Reals)
	model.w = Var(F, FPSO, T, S, within = Reals)
	model.g = Var(F, FPSO, T, S, within = Reals)
	model.x_fpso = Var(FPSO, T, S, within = Reals)
	model.w_fpso = Var(FPSO, T, S, within = Reals)
	model.g_fpso = Var(FPSO, T, S, within = Reals)
	
	model.Q_oil = Var(FPSO, T, S, within = Reals)
	model.Q_liq = Var(FPSO, T, S, within = Reals)
	model.Q_gas = Var(FPSO, T, S, within = Reals)
	model.QI_oil = Var(FPSO, T, S, within = Reals) 
	model.QE_oil = Var(FPSO, T, S, within = Reals)
	model.QI_liq = Var(FPSO, T, S, within = Reals)
	model.QE_liq = Var(FPSO, T, S, within = Reals)
	model.QI_gas = Var(FPSO, T, S, within = Reals)
	model.QE_gas = Var(FPSO, T, S, within = Reals)
	
	
	# Binary Variables
	model.b =  Var(FPSO,  T, S, within = Binary)
	model.b_c = Var(RF, FPSO,  T, S, within = Binary)
	model.b_on = Var(F, FPSO,  S, within = Binary)
	model.b_ex = Var(FPSO,  T, S, within = Binary)
	model.I_well = Var(F, T, S, within = NonNegativeIntegers)
	model.Z = Var(RF, i, T, S, within = Reals)
	model.b_prod = Var(F, T, S, within = Reals)
	
	model.w1 = Var(F, T, S, within = Reals)
	model.w2 = Var(F, T, S, within = Reals)
	model.w3 = Var(F, T, S, within = Reals)
	
	# auxiliary FPSO  OPER delta
	model.ZD_field = Var(FP, F, FPSO, T, S, within = Reals)
	model.ZD1_field = Var(FP, F, FPSO, T, S, within = Reals)
	model.ZD = Var(F, FPSO, T, S, within = Reals)
	
	#############################################################################################
	# objective function
	#############################################################################################
		
	def objectivemodel(model):
		return ENPV
	model.objective = Objective(sense=maximize, rule=objectivemodel)

	def con_2(model):
		return model.ENPV == sum(pro[s]* sum(dis[t] * (model.TConSh[t,s] - model.TCAP[t,s] - model.TOPER[t,s]) for t in T) for s in S)
	model.con_2m = Constraint(rule = con_2)
		
	def con_3(model, t, s):
		return model.TConSh[t,s] == sum(model.ConSh[rf,t,s] for rf in RF)
	model.con_3m = Constraint(T, S, rule = con_3)
	
	def con_4(model, t, s):
		return model.TCAP[t,s] == sum(model.CAP[rf,t,s] for rf in RF)
	model.con_4m = Constraint(T, S,  rule = con_4)
	
	def con_5(model, t, s):
		return model.TOPER[t,s] == sum(model.OPER[rf,t,s] for rf in RF)
	model.con_5m = Constraint(T, S, rule = con_5)
	
	# (ii) cost calculations
	def con_6(model, rf, t, s):
		return model.CAP[rf,t,s] == model.CAP1[rf,t,s] + model.CAP2[rf,t,s]
	model.con_6m = Constraint(RF, T, S, rule = con_6)

	def con_7(model, rf, t, s):
		return model.CAP1[rf,t,s] == sum(FC[f,fpso,t] * model.b_c[f,fpso,t,s] for fpso in FPSO for f in F[rf]) +\
									 sum(FC_well[f,t] * model.I_well[f,t,s] for f in F[rf])
	model.con_7m = Constraint(RF, T, S, rule = con_7)
	
	def con_8(model, rf, t, s):
		return model.CAP2[rf,t,s] == sum(model.DFPSOC[rf,fpso,t,s] for fpso in FPSO)
	model.con_8m = Constraint(RF, T, S, rule = con_8)
	
	def con_9(model,fpso, t, s):
		return model.FPSOC[fpso,t,s] == FC_FPSO[fpso,t]*model.b[fpso,t,s] + VC_liq[fpso,t]*(model.QI_liq[fpso,t,s] + model.QE_liq[fpso,t,s]) + \
										VC_gas[fpso,t]*(model.QI_gas[fpso,t,s] + model.QE_gas[fpso,t,s])
	model.con_9m = Constraint(FPSO, T, S, rule = con_9)
	
	def con_10(model,fpso, t, s):
		return model.FPSOC[fpso,t,s] == sum(model.DFPSOC_field[f,fpso,t,s] for f in F_fpso[fpso])
	model.con_10m = Constraint(FPSO, T, S, rule = con_10)
	
	def con_11(model,rf, fpso, t, s):
		return model.DFPSOC[rf,fpso,t,s] == sum(model.DFPSOC_field[f,fpso,t,s] for f in F_rf[rf])
	model.con_11m = Constraint(RF, FPSO, T, S, rule = con_11)
	
	def con_12(model,f, fpso, s):
		return model.b_on[f,fpso,s] == sum(model.b_c[f,fpso,t,s] for t in T)
	model.con_12m = Constraint(F, FPSO, S, rule = con_12)
	
	def con_13(model,f, fpso,t, s):
		return model.DFPSOC_field[f,fpso,t,s] <= M * model.b_on[f,fpso,s]
	model.con_13m = Constraint(F, FPSO, T, S, rule = con_13)
	
	#############################################################################################
	# Equation 14
	#############################################################################################
	def con_15(model,f, fpso,t, s):
		return sum(model.ZD_field[fp,f,fpso,t,s] * REC[fp,s] for fp in F_fpso[fpso]) == model.ZD[f,fpso,t,s] * REC[f,s]
	model.con_14m = Constraint(F, FPSO, T, S, rule = con_14)
	
	def con_16(model,f,fpso,t,fp,s):
		return model.ZD_field[fp,f,fpso,t,s] + model.ZD1_field[fp,f,fpso,t,s] == model.DFPSOC_field[f,fpso,t,s]
	model.con_16m = Constraint(F, FPSO, T, F_fpso[fpso], S, rule = con_16)

	def con_17(model,f,fpso,t,fp,s):
		return model.ZD_field[fp,f,fpso,t,s] <= U*model.b_on[fp,fpso,s]
	model.con_17m = Constraint(F, FPSO, T, F_fpso[fpso], S, rule = con_17)
	
	def con_18(model,f,fpso,t,fp,s):
		return model.ZD1_field[fp,f,fpso,t,s] <= U*(1-model.b_on[fp,fpso,s])
	model.con_18m = Constraint(F, FPSO, T, F_fpso[fpso], S, rule = con_18)
	
	def con_20(model,f,fpso,t,s):
		return model.ZD[f,fpso,t,s] + model.ZD1[f,fpso,t,s] == model.FPSOC[fpso,t,s]
	model.con_20m = Constraint(F,FPSO,T,S, rule = con_20)
	
	def con_21(model, f,fpso,t,s):
		return model.ZD[f,fpso,t,s] <= U*model.b_on[f,fpso,s]
	model.con_21m = Constraint(F,FPSO,T,S, rule = con_21)
	
	def con_22(model, f,fpso,t,s):
		return model.ZD1[f,fpso,t,s] <= U*(1-model.b_on[f,fpso,s])
	model.con_22m = Constraint(F,FPSO,T,S, rule = con_22)
	
	# operating cost
	def con_24(model, rf,t,s):
		return model.OPER[rf,t,s] == delta[t]*(OC_liq[rf,t]*(model.x_tot[rf,t,s] + model.w_tot[rf,t,s]) + OC_gas[rf,t]*model.g_tot[rf,t,s])
	model.con_24m = Constraint(RF,T,S, rule = con_24)
	
	# (iii) total contractor share calculations
	def con_25(model, rf,t,s):
		return model.ConSh[rf,t,s] == model.ConSh_aftertax[rf,t,s] + model.CO[rf,t,s]
	model.con_25m = Constraint(RF,T,S, rule = con_25)
	
	def con_26(model, rf,t,s):
		return model.ConSh_aftertax[rf,t,s] == model.ConSh_beforetax[rf,t,s] - model.Tax[rf,t,s]
	model.con_26m = Constraint(RF,T,S, rule = con_26)
	
	def con_27(model, rf,t,s):
		return model.Tax[rf,t,s] == f_tax[rf,t]*model.ConSh_beforetax[rf,t,s]
	model.con_27m = Constraint(RF,T,S, rule = con_27)
	
	#############################################################################################
	# Equation 28
	#############################################################################################
	def con_29(model, rf,t,s):
		return model.ConSh_beforetax[rf,t,s] == sum(model.DConSh_beforetax[rf,i,t, s] for i in I)
	model.con_29m = Constraint(RF,T,S, rule = con_29)
	
	def con_30(model, rf,t,s):
		return model.PO[rf,t,s] == sum(model.DPO[rf,i,t, s] for i in I)
	model.con_30m = Constraint(RF,T,S, rule = con_30)
	
	def con_31(model, rf,t,s):
		return model.xc[rf,t,s] == sum(model.Dxc[rf,i,t, s] for i in I)
	model.con_31m = Constraint(RF,T,S, rule = con_31)
	
	def con_32(model, rf,i,t,s):
		return model.DConSh_beforetax[rf,i,t,s] == f_po[rf,i] * model.DPO[rf,i,t, s] 
	model.con_32m = Constraint(RF,I,T,S, rule = con_32)
	
	def con_33(model, rf,i,t,s):
		return model.DConSh_beforetax[rf,i,t,s] <= M * model.Z[rf,i,t, s] 
	model.con_33m = Constraint(RF,I,T,S, rule = con_33)
	
	def con_34(model, rf,i,t,s):
		return model.DPO[rf,i,t, s] <= M * model.Z[rf,i,t,s] 
	model.con_34m = Constraint(RF,I,T,S, rule = con_34)
	
	def con_35a(model, rf,i,t,s):
		return model.Dxc[rf,i,t, s] <= U_oil[rf,i] * model.Z[rf,i,t,s] 
	model.con_35am = Constraint(RF,I,T,S, rule = con_35a)
	
	def con_35b(model, rf,i,t,s):
		return model.Dxc[rf,i,t, s] >= L_oil[rf,i] * model.Z[rf,i,t,s] 
	model.con_35bm = Constraint(RF,I,T,S, rule = con_35b)
	
	def con_36(model, rf,t,s):
		return sum(model.Z[rf,i,t,s] for i in I) == 1
	model.con_36m = Constraint(RF,T,S, rule = con_36)
	
	def con_37(model, rf,t,s):
		return model.xc[rf,t,s] == sum(model.xc_field[f,t,s] for f in F_rf[rf])
	model.con_37m = Constraint(RF,T,S, rule = con_37)
	
	def con_38(model, rf,t,s):
		return model.PO[rf,t,s] == model.REV[rf,t,s] - model.CO[rf,t,s]
	model.con_38m = Constraint(RF,T,S, rule = con_38)
	
	def con_39(model, rf,t,s):
		return model.REV[rf,t,s] == delta[t] * alpha[t] * model.x_tot[rf,t,s]
	model.con_39m = Constraint(RF,T,S, rule = con_39)
	
	def con_40(model, rf,t,s):
		return model.x_tot[rf,t,s] == sum(model.x_f[f,t,s] for f in F_rf[rf])
	model.con_40m = Constraint(RF,T,S, rule = con_40)
	
	#############################################################################################
	# Equation 41
	#############################################################################################
	def con_42(model, rf,t,s):
		return model.CO[rf,t,s] <= model.CR[rf,t,s] + M * (1 - model.b_co[rf,t,s])
	model.con_42m = Constraint(RF,T,S, rule = con_42)
	
	def con_43(model, rf,t,s):
		return model.CO[rf,t,s] >= model.CR[rf,t,s] - M * (1 - model.b_co[rf,t,s])
	model.con_43m = Constraint(RF,T,S, rule = con_43)

	def con_44(model, rf,t,s):
		return model.CO[rf,t,s] <= f_CR[rf,t] * model.REV[rf,t,s] + M * model.b_co[rf,t,s]
	model.con_44m = Constraint(RF,T,S, rule = con_44)

	def con_45(model, rf,t,s):
		return model.CO[rf,t,s] >= f_CR[rf,t] * model.REV[rf,t,s] - M * model.b_co[rf,t,s]
	model.con_45m = Constraint(RF,T,S, rule = con_45)

	def con_46(model, rf,t,s):
		return model.CO[rf,t,s] <= model.CR[rf,t,s]
	model.con_46m = Constraint(RF,T,S, rule = con_46)

	def con_47(model, rf,t,s):
		return model.CO[rf,t,s] <= f_CR[rf,t] * model.REV[rf,t,s]
	model.con_47m = Constraint(RF,T,S, rule = con_47)
	
	def con_48(model, rf,t,s):
		return model.CR[rf,t,s] == model.CAP[rf,t,s] + model.OPER[rf,t,s] + model.CRF[rf,t-1,s]
	model.con_48m = Constraint(RF,T,S, rule = con_48)
	
	def con_49(model, rf,t,s):
		return model.CRF[rf,t,s] == model.CR[rf,t,s] - model.CO[rf,t,s]
	model.con_49m = Constraint(RF,T,S, rule = con_49)

	#############################################################################################
	# Equation 50 and 51
	#############################################################################################
	def con_52(model, rf, i, ip, t, tao, s):
		if ip < i and t<=tao<=T[-1]:
			return model.Z[rf,i,t,s] + model.Z[rf,ip,tao,s] <= 1
		else:
			return Constraint.Skip
	model.con_52m = Constraint(RF, I, I, T, T, S, rule = con_52)

	def con_53(model, rf, i, ip, t, tao, s):
		if ip > i and 1<=tao<=t:
			return model.Z[rf,i,t,s] + model.Z[rf,ip,tao,s] <= 1
		else:
			return Constraint.Skip
	model.con_53m = Constraint(RF, I, I, T, T, S, rule = con_53)

	def con_54(model, rf, i, t, s):
		return sum(model.Contsh_beforetax[rf,tao,s] for tao in range(1,t+1)) <= sum((f_PO[rf,ip] - f_PO[rf,ip-1])*(model.xc[rf,t,s] - L[rf,ip]) for ip in range(1, i+1)) -\
																				f_PO[rf,I[-1]] * sum(model.CO[rf,tao,s]/alpha[tao] for tao in range(1, t+1))
	model.con_54m = Constraint(RF, I, T, S, rule = con_54)

	# (v) reservoir constraints
	def con_55(model, f, fpso, t, s):
		return model.x_well[f,fpso,t,s] <= model.Q_d_well[f,fpso,t,s]
	model.con_55m = Constraint(F, FPSO, T, S, rule = con_55)

	def con_56a(model, f, fpso, s):
		return model.Q_d_well[f,fpso,1,s] == alpha_oil[s] * d[1,f,fpso]
	model.con_56am = Constraint(F, FPSO, S, rule = con_56a)

	def con_56b(model, f, fpso, t, s):
		if t < T[-1]:
			return model.Q_d_well[f,fpso,t+1,s] == alpha_oil[s] * (a[1,f,fpso] * model.fc[f,t,s]**3 + \
																b[1,f,fpso] * model.fc[f,t,s]**2 + \
																c[1,f,fpso] * model.fc[f,t,s] + d[1,f,fpso])
		else:
			return Constraint.Skip
	model.con_56bm = Constraint(F, FPSO, T, S, rule = con_56b)

	def con_57(model, f, fpso, t, s):
		return model.Q_wc[f,fpso,t,s] == alpha_wc[s] * (a[2,f,fpso] * model.fc[f,t,s]**4 + \
														b[2,f,fpso] * model.fc[f,t,s]**3 + \
														c[2,f,fpso] * model.fc[f,t,s]**2 + \
														d[2,f,fpso] * model.fc[f,t,s])
	model.con_57m = Constraint(F, FPSO, T, S, rule = con_57)

	def con_58(model, f, fpso, t, s):
		return model.Q_gc[f,fpso,t,s] == alpha_gc[s] * (a[3,f,fpso] * model.fc[f,t,s]**4 + \
														b[3,f,fpso] * model.fc[f,t,s]**3 + \
														c[3,f,fpso] * model.fc[f,t,s]**2 + \
														d[3,f,fpso] * model.fc[f,t,s])
	model.con_58m = Constraint(F, FPSO, T, S, rule = con_58)
	
	def con_59(model, f, fpso, t, s):
		return model.wc[f,fpso,t,s] <= model.Q_wc[f,fpso,t,s] + M * (1 - sum(model.b_c[f,fpso,tao,s] for tao in range(1,t+1)))
	model.con_59m = Constraint(F, FPSO, T, S, rule = con_59)

	def con_60(model, f, fpso, t, s):
		return model.wc[f,fpso,t,s] >= model.Q_wc[f,fpso,t,s] - M * (1 - sum(model.b_c[f,fpso,tao,s] for tao in range(1,t+1)))
	model.con_60m = Constraint(F, FPSO, T, S, rule = con_60)

	def con_61(model, f, fpso, t, s):
		return model.wc[f,fpso,t,s] <= M_wc[f,fpso,s] * sum(model.b_c[f,fpso,tao,s] for tao in range(1,t+1))
	model.con_61m = Constraint(F, FPSO, T, S, rule = con_61)

	def con_62(model, f, fpso, t, s):
		return model.wc[f,fpso,t,s] >= -M_wc[f,fpso,s] * sum(model.b_c[f,fpso,tao,s] for tao in range(1,t+1))
	model.con_62m = Constraint(F, FPSO, T, S, rule = con_62)

	def con_63(model, f, fpso, t, s):
		return model.gc[f,fpso,t,s] <= model.Q_gc[f,fpso,t,s] + M_gc[f,fpso,s] * (1 - sum(model.b_c[f,fpso,tao,s] for tao in range(1,t+1)))
	model.con_63m = Constraint(F, FPSO, T, S, rule = con_63)

	def con_64(model, f, fpso, t, s):
		return model.gc[f,fpso,t,s] >= model.Q_gc[f,fpso,t,s] - M_gc[f,fpso,s] * (1 - sum(model.b_c[f,fpso,tao,s] for tao in range(1,t+1)))
	model.con_64m = Constraint(F, FPSO, T, S, rule = con_64)

	def con_65(model, f, fpso, t, s):
		return model.gc[f,fpso,t,s] <= M_gc[f,fpso,s] * sum(model.b_c[f,fpso,tao,s] for tao in range(1,t+1))
	model.con_65m = Constraint(F, FPSO, T, S, rule = con_65)

	def con_66(model, f, fpso, t, s):
		return model.gc[f,fpso,t,s] >= -M_gc[f,fpso,s] * sum(model.b_c[f,fpso,tao,s] for tao in range(1,t+1))
	model.con_66m = Constraint(F, FPSO, T, S, rule = con_66)

	def con_67(model, f, fpso, t, s):
		return model.w[f,fpso,t,s] == (model.wc[f,fpso,t,s] - model.wc[f,fpso,t-1,s])/delta[t]
	model.con_67m = Constraint(F, FPSO, T, S, rule = con_67)

	def con_68(model, f, fpso, t, s):
		return model.g[f,fpso,t,s] == (model.gc[f,fpso,t,s] - model.gc[f,fpso,t-1,s])/delta[t]
	model.con_68m = Constraint(F, FPSO, T, S, rule = con_68)

	# (vi) Field-FPSO flow constraints
	def con_69(model, f, t, s):
		return model.x_f[f,t,s] == sum(model.x_con[f,fpso,t,s] for fpso in FPSO)
	model.con_69m = Constraint(F, T, S, rule = con_69)

	# not know N
	def con_70(model, f, fpso, t, s):
		return model.x_c[f,fpso,t,s] == model.N_well[f,t,s] * model.x_well[f,fpso,t,s]
	model.con_70m = Constraint(F, FPSO, T, S, rule = con_70)

	def con_72(model, f, t, s):
		return model.fc[f,t,s] == sum(model.x_f[f,tao,s]*delta[tao] for tao in range(1,t+1))/REC[f,s]
	model.con_72m = Constraint(F, T, S, rule = con_72)

	def con_73(model, f, t, s):
		return sum(model.x_f[f,tao,s]*delta[tao] for tao in range(1,t+1)) <= REC[f,s]
	model.con_73m = Constraint(F, T, S, rule = con_73)

	def con_74(model, fpso, t, s):
		return model.x_fpso[fpso,t,s] == sum(model.x_c[f,fpso,t,s] for f in F)
	model.con_74m = Constraint(FPSO, T, S, rule = con_74)

	def con_75(model, fpso, t, s):
		return model.w_fpso[fpso,t,s] == sum(model.w[f,fpso,t,s] for f in F)
	model.con_75m = Constraint(FPSO, T, S, rule = con_75)

	def con_76(model, fpso, t, s):
		return model.g_fpso[fpso,t,s] == sum(model.g[f,fpso,t,s] for f in F)
	model.con_76m = Constraint(FPSO, T, S, rule = con_76)

	def con_77(model, t, s):
		return model.x_tot[t,s] == sum(model.x_fpso[fpso,t,s] for fpso in FPSO)
	model.con_77m = Constraint(T, S, rule = con_77)

	def con_78(model, t, s):
		return model.w_tot[t,s] == sum(model.w_fpso[fpso,t,s] for fpso in FPSO)
	model.con_78m = Constraint(T, S, rule = con_78)

	def con_79(model, t, s):
		return model.g_tot[t,s] == sum(model.g_fpso[fpso,t,s] for fpso in FPSO)
	model.con_79m = Constraint(T, S, rule = con_79)

	# (vii) FPSO capacity constraints
	def con_80(model, fpso, t, s):
		return model.x_fpso[fpso,t,s] <= model.Q_oil[fpso,t,s]
	model.con_80m = Constraint(FPSO, T, S, rule = con_80)

	def con_81(model, fpso, t, s):
		return model.x_fpso[fpso,t,s] + model.w_fpso[fpso,t,s] <= model.Q_liq[fpso,t,s]
	model.con_81m = Constraint(FPSO, T, S, rule = con_81)

	def con_82(model, fpso, t, s):
		return model.g_fpso[fpso,t,s] <= model.Q_gas[fpso,t,s]
	model.con_82m = Constraint(FPSO, T, S, rule = con_82)
	
	#
	def con_83(model, fpso, t, s):
		return model.Q_oil[fpso,t,s] == model.Q_oil[fpso,t-1,s] + model.QI_oil[fpso,t-l1,s] + model.QE_oil[fpso,t-l2,s]
	model.con_83m = Constraint(FPSO, T, S, rule = con_83)
	
	def con_84(model, fpso, t, s):
		return model.Q_liq[fpso,t,s] == model.Q_liq[fpso,t-1,s] + model.QI_liq[fpso,t-l1,s] + model.QE_liq[fpso,t-l2,s]
	model.con_84m = Constraint(FPSO, T, S, rule = con_84)
	
	def con_85(model, fpso, t, s):
		return model.Q_gas[fpso,t,s] == model.Q_gas[fpso,t-1,s] + model.QI_gas[fpso,t-l1,s] + model.QE_gas[fpso,t-l2,s]
	model.con_85m = Constraint(FPSO, T, S, rule = con_85)
	
	# (viii) logic constraints
	def con_86(model, fpso, s):
		return sum(model.b[fpso,t,s] for t in T) <= 1
	model.con_86m = Constraint(FPSO, S, rule = con_86)
	
	def con_87(model, fpso, s):
		return sum(model.b_ex[fpso,t,s] for t in T) <= 1
	model.con_87m = Constraint(FPSO, S, rule = con_87)

	def con_88(model, f, fpso, s):
		return sum(model.b_c[f, fpso,t,s] for t in T) <= 1
	model.con_88m = Constraint(f, FPSO, S, rule = con_88)

	#
	def con_89(model, f, t, s):
		return sum(model.b_c[f, fpso,t,s] for fpso in FPSO) <= 1
	model.con_89m = Constraint(F, T, S, rule = con_89)

	def con_90(model, f,  s):
		return sum(model.b_c[f, fpso,t,s] for fpso in FPSO for t in T) <= 1
	model.con_90m = Constraint(F, S, rule = con_90)

	#
	def con_91(model, fpso, t, s):
		return model.b_ex[fpso,t,s]  <= sum(model.b[fpso,tao,s] for tao in range(1,t+1))
	model.con_91m = Constraint(FPSO, T, S, rule = con_91)

	def con_92(model, f, fpso, t, s):
		return model.b_c[f,fpso,t,s]  <= sum(model.b[fpso,tao,s] for tao in range(1,t+1))
	model.con_92m = Constraint(F, FPSO, T, S, rule = con_92)

	# (ix) upper bounding constraints
	def con_93(model, f, fpso, t, s):
		return model.x_well[f,fpso,t,s]  <= U_well_oil[f,fpso] * sum(model.b_c[f,fpso,tao,s] for tao in range(1,t+1))
	model.con_93m = Constraint(F, FPSO, T, S, rule = con_93)

	#
	def con_94(model, fpso, t, s):
		return model.QI_oil[fpso,t,s]  <= U_oil[fpso] * model.b[fpso,t,s]
	model.con_94m = Constraint(FPSO, T, S, rule = con_94)

	def con_95(model, fpso, t, s):
		return model.QI_liq[fpso,t,s]  <= U_liq[fpso] * model.b[fpso,t,s]
	model.con_95m = Constraint(FPSO, T, S, rule = con_95)

	def con_96(model, fpso, t, s):
		return model.QI_gas[fpso,t,s]  <= U_gas[fpso] * model.b[fpso,t,s]
	model.con_96m = Constraint(FPSO, T, S, rule = con_96)

	def con_97(model, fpso, t, s):
		return model.QE_oil[fpso,t,s]  <= U_oil[fpso] * model.b_ex[fpso,t,s]
	model.con_97m = Constraint(FPSO, T, S, rule = con_97)

	def con_98(model, fpso, t, s):
		return model.QE_liq[fpso,t,s]  <= U_liq[fpso] * model.b_ex[fpso,t,s]
	model.con_98m = Constraint(FPSO, T, S, rule = con_98)

	def con_99(model, fpso, t, s):
		return model.QE_gas[fpso,t,s]  <= U_gas[fpso] * model.b_ex[fpso,t,s]
	model.con_99m = Constraint(FPSO, T, S, rule = con_99)

	#
	def con_100(model, fpso, t, s):
		return model.QE_oil[fpso,t,s]  <= u * model.Q_oil[fpso,t-1,s]
	model.con_100m = Constraint(FPSO, T, S, rule = con_100)

	def con_101(model, fpso, t, s):
		return model.QE_liq[fpso,t,s]  <= u * model.Q_liq[fpso,t-1,s]
	model.con_101m = Constraint(FPSO, T, S, rule = con_101)

	def con_102(model, fpso, t, s):
		return model.QE_gas[fpso,t,s]  <= u * model.Q_gas[fpso,t-1,s]
	model.con_102m = Constraint(FPSO, T, S, rule = con_102)

	# (x) well drilling limitations
	def con_103(model, f, t, s):
		return model.N_well[f,t,s]  == model.N_well[f,t-1,s] + model.I_well[f,t,s]
	model.con_103m = Constraint(F, T, S, rule = con_103)

	#
	def con_104(model, t, s):
		return sum(model.I_well[f,t,s] for f in F) <= UI_well[t]
	model.con_104m = Constraint(T, S, rule = con_104)

	def con_105(model, f, t, s):
		return model.N_well[f,t,s] <= UN_well[f]
	model.con_105m = Constraint(F, T, S, rule = con_105)

	#############################################################################################
	# Initial NACs
	#############################################################################################
	def con_106(model, fpso, t, s, sp):
		return model.b[fpso,t,s] == model.b[fpso,t,sp]
	model.con_106m = Constraint(FPSO, T1, S, S, rule = con_106)

	def con_107(model, fpso, t, s, sp):
		return model.b_ex[fpso,t,s] == model.b_ex[fpso,t,sp]
	model.con_107m = Constraint(FPSO, T1, S, S, rule = con_107)

	def con_108(model, f,fpso, t, s, sp):
		return model.b_c[f,fpso,t,s] == model.b_c[f,fpso,t,sp]
	model.con_108m = Constraint(F, FPSO, T1, S, S, rule = con_108)
	
	def con_109(model, f, t, s, sp):
		return model.I_well[f,t,s] == model.I_well[f,t,sp]
	model.con_109m = Constraint(F, T1, S, S, rule = con_109)

	def con_110(model, fpso, t, s, sp):
		return model.QI_oil[fpso,t,s] == model.QI_oil[fpso,t,sp]
	model.con_110m = Constraint(FPSO, T1, S, S, rule = con_110)

	def con_111(model, fpso, t, s, sp):
		return model.QI_liq[fpso,t,s] == model.QI_liq[fpso,t,sp]
	model.con_111m = Constraint(FPSO, T1, S, S, rule = con_111)

	def con_112(model, fpso, t, s, sp):
		return model.QI_gas[fpso,t,s] == model.QI_gas[fpso,t,sp]
	model.con_112m = Constraint(FPSO, T1, S, S, rule = con_112)
	
	def con_113(model, fpso, t, s, sp):
		return model.QE_oil[fpso,t,s] == model.QE_oil[fpso,t,sp]
	model.con_113m = Constraint(FPSO, T1, S, S, rule = con_113)

	def con_114(model, fpso, t, s, sp):
		return model.QE_liq[fpso,t,s] == model.QE_liq[fpso,t,sp]
	model.con_114m = Constraint(FPSO, T1, S, S, rule = con_114)

	def con_115(model, fpso, t, s, sp):
		return model.QE_gas[fpso,t,s] == model.QE_gas[fpso,t,sp]
	model.con_115m = Constraint(FPSO, T1, S, S, rule = con_115)
	
	#############################################################################################
	# Equation 116
	#############################################################################################
	
	def con_116a(model, f, t, s):
		return model.N_well[f,t-1,s] - N1 + 1 <= M * (1 - model.w1[f,t,s])
	model.con_116am = Constraint(F, T, S, rule = con_116a)

	def con_116b(model, f, t, s):
		return model.N_well[f,t-1,s] - N1 + 1 >= -M * model.w1[f,t,s] + 0.1
	model.con_116bm = Constraint(F, T, S, rule = con_116b)
	
	#############################################################################################
	# Equation 117
	#############################################################################################
	
	def con_117a(model, f, t, s):
		return sum(model.b_prod[f,tao,s] for tao in range(1,t)) - N2 + 1 <= M * (1 - model.w2[f,t,s])
	model.con_117am = Constraint(F, T, S, rule = con_117a)

	def con_117b(model, f, t, s):
		return sum(model.b_prod[f,tao,s] for tao in range(1,t)) - N2 + 1 >= -M * model.w2[f,t,s] + 0.1
	model.con_117bm = Constraint(F, T, S, rule = con_117b)
	
	#############################################################################################
	# Equation 118
	#############################################################################################

	def con_118a(model, f, t, s):
		return sai - model.x_f[f,t,s] <= M * (1 - model.b_prod[f,t,s])
	model.con_118am = Constraint(F, T, S, rule = con_118a)

	def con_118b(model, f, t, s):
		return sai - model.x_f[f,t,s] >= -M * model.b_prod[f,t,s] + 0.1
	model.con_118bm = Constraint(F, T, S, rule = con_118b)
	
	#############################################################################################
	# Equation 119
	#############################################################################################
	
	def con_119a(model, f, t, s):
		return model.w1[f,t,s] >= model.w3[f,t,s]
	model.con_119am = Constraint(F, T, S, rule = con_119a)
	
	def con_119b(model, f, t, s):
		return model.w2[f,t,s] >= model.w3[f,t,s]
	model.con_119bm = Constraint(F, T, S, rule = con_119b)

	def con_119c(model, f, t, s):
		return model.w3[f,t,s] >= model.w1[f,t,s] + model.w2[f,t,s] - 1
	model.con_119cm = Constraint(F, T, S, rule = con_119c)
	
	#############################################################################################
	# Equation 120
	#############################################################################################
	
	def con_120a(model, t, s, sp):
		if s < sp and (s,sp) in L1 and tao<=t and f in D[s,sp]:
			return model.Z[t,s,sp] <= model.w3[f,t,s]
		else:
			return Constraint.Skip
	model.con_120am = Constraint(T, S, S, rule = con_120a)
	
	def con_120b(model, t, s, sp):
		if s < sp and (s,sp) in L1 and f in D[s,sp]:
			return model.Z[t,s,sp] >= 1 - t + sum(model.w3[f,tao,s] for tao in range(1,t+1))
		else:
			return Constraint.Skip
	model.con_120bm = Constraint(T, S, S, rule = con_120b)
	
	#############################################################################################
	# Conditional NACs 121
	#############################################################################################
	def NAC1a(model, fpso, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return model.b[fpso,t,s] -  model.b[fpso,t,sp] <= 1 - model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC1am = Constraint(FPSO, T, S, S, rule = NAC1a)
	
	def NAC1b(model, fpso, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return  model.b[fpso,t,s] -  model.b[fpso,t,sp] >= -1 + model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC1bm = Constraint(FPSO, T, S, S, rule = NAC1b)

	def NAC2a(model, fpso, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return model.b_ex[fpso,t,s] - model.b_ex[fpso,t,sp] <= 1 - model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC2am = Constraint(FPSO, T, S, S, rule = NAC2a)
	
	def NAC2b(model, fpso, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return  model.b_ex[fpso,t,s] - model.b_ex[fpso,t,sp] >= -1 + model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC2bm = Constraint(FPSO, T, S, S, rule = NAC2b)

	def NAC3a(model, fpso, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return model.b_c[f,fpso,t,s] - model.b_c[f,fpso,t,sp] <= 1 - model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC3am = Constraint(FPSO, T, S, S, rule = NAC3a)
	
	def NAC3b(model, fpso, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return  model.b_c[f,fpso,t,s] - model.b_c[f,fpso,t,sp] >= -1 + model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC3bm = Constraint(FPSO, T, S, S, rule = NAC3b)
	
	def NAC4a(model, f, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return model.I_well[f,t,s] - model.I_well[f,t,sp] <= 1 - model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC4am = Constraint(F, T, S, S, rule = NAC4a)
	
	def NAC4b(model, f, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return  model.I_well[f,t,s] - model.I_well[f,t,sp] >= -1 + model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC4bm = Constraint(F, T, S, S, rule = NAC4b)

	def NAC5a(model, fpso, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return model.QI_oil[fpso,t,s] - model.QI_oil[fpso,t,sp] <= 1 - model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC5am = Constraint(FPSO, T, S, S, rule = NAC5a)
	
	def NAC5b(model, fpso, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return model.QI_oil[fpso,t,s] - model.QI_oil[fpso,t,sp] >= -1 + model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC5bm = Constraint(FPSO, T, S, S, rule = NAC5b)
	
	def NAC6a(model, fpso, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return model.QI_liq[fpso,t,s] - model.QI_liq[fpso,t,sp] <= 1 - model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC6am = Constraint(FPSO, T, S, S, rule = NAC6a)
	
	def NAC6b(model, fpso, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return model.QI_liq[fpso,t,s] - model.QI_liq[fpso,t,sp] >= -1 + model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC6bm = Constraint(FPSO, T, S, S, rule = NAC6b)
	
	def NAC7a(model, fpso, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return model.QI_gas[fpso,t,s] - model.QI_gas[fpso,t,sp] <= 1 - model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC7am = Constraint(FPSO, T, S, S, rule = NAC7a)
	
	def NAC7b(model, fpso, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return model.QI_gas[fpso,t,s] - model.QI_gas[fpso,t,sp] >= -1 + model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC7bm = Constraint(FPSO, T, S, S, rule = NAC7b)

	def NAC8a(model, fpso, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return model.QE_oil[fpso,t,s] - model.QE_oil[fpso,t,sp] <= 1 - model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC8am = Constraint(FPSO, T, S, S, rule = NAC8a)
	
	def NAC8b(model, fpso, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return model.QE_oil[fpso,t,s] - model.QE_oil[fpso,t,sp] >= -1 + model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC8bm = Constraint(FPSO, T, S, S, rule = NAC8b)

	def NAC9a(model, fpso, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return model.QE_liq[fpso,t,s] - model.QE_liq[fpso,t,sp] <= 1 - model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC9am = Constraint(FPSO, T, S, S, rule = NAC9a)
	
	def NAC9b(model, fpso, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return model.QE_liq[fpso,t,s] - model.QE_liq[fpso,t,sp] >= -1 + model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC9bm = Constraint(FPSO, T, S, S, rule = NAC9b)

	def NAC10a(model, fpso, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return model.QE_gas[fpso,t,s] - model.QE_gas[fpso,t,sp] <= 1 - model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC10am = Constraint(FPSO, T, S, S, rule = NAC10a)
	
	def NAC10b(model, fpso, t, s, sp):
		if s < sp and (s,sp) in L1 and t+1<=T[-1]:
			return model.QE_gas[fpso,t,s] - model.QE_gas[fpso,t,sp] >= -1 + model.Z[t,s,sp]
		else:
			return Constraint.Skip
	model.NAC10bm = Constraint(FPSO, T, S, S, rule = NAC10b)

	results= opt.solve(model)
	model.solutions.load_from(results)
	
	save_file1 = "SMILP" + '_' + str(T[-1]) + ' ' + '(' + str(today) + ')'
	save_file = "SMILPOutput" + '_' + str(T[-1]) + ' ' + '(' + str(today) + ')'
	dir_path = os.path.dirname(os.path.realpath(__file__))
	# TempfileManager.tempdir = dir_path

	f = open(os.path.join(save_file),	"w")
	results.write(filename=save_file1)
	
	return
