####################################################################
#### 						Sets
####################################################################
sets = {}
sets['WP'] = ['A','B','C','D','E','F']
sets['PP'] = ['PP']
sets['T'] = range(1,16)
####################################################################
#### 					Parameters
####################################################################
parameters = {}

parameters['deliverbility'] = {'A':130,'B':200,'C':100,'D':100,'E':130,'F':130}
parameters['size'] = {'A':400,'B':400,'C':350,'D':200,'E':290}
parameters['probability'] = {1:0.3, 2:0.4, 3:0.3}
parameters['shrink'] = 0.01
parameters['Big_M'] = 4

####################################################################
###          Uncertainty Information
####################################################################
Realizations = {}

Realizations['outcomes'] = [[100, 300, 470]]

Uncertain = {}
Uncertain['WP_uncertain'] = ['F']


