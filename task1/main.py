# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 22:04:12 2023

@author: barnabaspiri
"""
#%%
from IPython import get_ipython
ipython = get_ipython()

try:
    from IPython import get_ipython
    get_ipython().magic('clear')
    get_ipython().magic('reset -f')
except:
    pass

#%%
# Setting the path to acces modules
import os
import sys

file_dir = os.path.dirname(__file__)

print(f"Working directory: {file_dir}")
sys.path.append(file_dir)

import numpy as np

# Geometric data
L = 3000 # mm
A = 20 # mm
B = 50 # mm

# Material data
E = 2.15e5 # MPa
nu = 0.32 # 1

# Load data
F = 100 # N

# DOF number of the system
DOF = 10

# Beam divided into 4 equal elements
# No. 1.: 0 - 750
# No. 2.: 750 - 1500
# No. 3.: 750 - 2250
# No. 4.: 2250 - 3000

Le = L/4
Ie = A**3*B / 12
Ee = E

ecs = np.matrix([[1,2],[2,3],[3,4],[4,5]]) # element - node table

from stiffness_matrix import Ke, KGe

K1 = Ke(Ie, Ee, Le)
K2 = Ke(Ie, Ee, Le)
K3 = Ke(Ie, Ee, Le)
K4 = Ke(Ie, Ee, Le)

KG1 = KGe(-2*F, Ee, Le)
KG2 = KGe(-2*F, Ee, Le)
KG3 = KGe(-F, Ee, Le)
KG4 = KGe(-F, Ee, Le)

from nodal_disp import eDOF

eDOF1 = eDOF(ecs[0])
eDOF2 = eDOF(ecs[1])
eDOF3 = eDOF(ecs[2])
eDOF4 = eDOF(ecs[3])

from glob_cond import ExtMatrix, SubMatrix

K = ExtMatrix(K1, eDOF1, DOF) + ExtMatrix(K2, eDOF2, DOF) + ExtMatrix(K3, eDOF3, DOF) + ExtMatrix(K4, eDOF4, DOF)

KG = ExtMatrix(KG1, eDOF1, DOF) + ExtMatrix(KG2, eDOF2, DOF) + ExtMatrix(KG3, eDOF3, DOF) + ExtMatrix(KG4, eDOF4, DOF)

freeDOF = np.matrix([3,4,5,6,7,8,9,10])

K_cond = SubMatrix(K, freeDOF)
KG_cond = SubMatrix(KG, freeDOF)

A = np.linalg.inv(KG_cond) @ K_cond

(eigenVALS, eigenVECS) = np.linalg.eig(A)

lambdas = -eigenVALS
lambdas_sort = np.sort(lambdas)

print(lambdas_sort)

forces_sort = F * lambdas_sort

print(forces_sort)

#%% symbolic eigenvals.

#import sympy as sp

#K_cond_sym = sp.Matrix(K_cond)
#KG_cond_sym = sp.Matrix(KG_cond)

#lam = sp.symbols('lambda')

#eq = sp.det(K_cond_sym + lam*KG_cond_sym)

#sol = sp.solve(eq,lam)

#print(sol)
