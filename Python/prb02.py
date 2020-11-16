# -*- coding: utf-8 -*-
"""
Problem 2.
   Xylene, styrene, toluene, and benzene are to be separated with the
   array of distillation columns.  F, D, B, D1, B1, D2, and
   B2 are the molar flow rates in mol/min. Find stream molar flow rates
   and composition

"""


# -----------------------------------------------------------------------------
# Mohammad Rahmani
# Chemical Engineering Department
# Amirkabir University of Technology
# Tehran, Iran
# m.rahmani@aut.ac.ir
#
# Rev 0.1
# Nov 16, 2020


import numpy as np


# Material balances on individual components on the overall separation train

# coefficient matrix
A=np.array([
    [0.07, 0.18, 0.15, 0.24],
    [0.04, 0.24, 0.10, 0.65],
    [0.54, 0.42, 0.54, 0.10],
    [0.35, 0.16, 0.21, 0.01]
    ]);

# the right hand side
f = np.array([
    0.15*70,
    0.25*70,
    0.40*70,
    0.2*70
    ]);

# Solve using numpy.linlang
X=np.linalg.solve(A, f)

# Store results
D1 =X[0]
B1 =X[1]
D2 =X[2]
B2 =X[3]

print('\nSeparation Train ')
print('D1 (mol/min): {:6.4f}'.format(D1))
print('B1 (mol/min): {:6.4f}'.format(B1))
print('D2 (mol/min): {:6.4f}'.format(D2))
print('B2 (mol/min): {:6.4f}'.format(B2))

#Use the balance over column 2
D=D1+B1						#43.75 mol/min

#Solve for Column 2
X_Dx=(0.07*D1+0.18*B1)/D	#0.114 mole fraction
X_Ds=(0.04*D1+0.24*B1)/D	#0.120 mole fraction
X_Dt=(0.54*D1+0.42*B1)/D	#0.492 mole fraction
X_Db=(0.35*D1+0.16*B1)/D	#0.274 mole fraction


print('\nColumn 2 ')
print('D (mol/min): {:6.4f}'.format(D))
print('Xylene,  X_Dx: {:6.4f}'.format(X_Dx))
print('Styrene, X_Ds: {:6.4f}'.format(X_Ds))
print('Toluene, X_Dt: {:6.4f}'.format(X_Dt))
print('Benzene, X_Db: {:6.4f}'.format(X_Db))


#Use the balance over column 3
B=D2+B2						#26.25 mol/min

#Solve for Column 3
X_Bx=(0.15*D2+0.24*B2)/B	#0.2100 mole fraction
X_Bs=(0.10*D2+0.65*B2)/B	#0.4667 mole fraction
X_Bt=(0.54*D2+0.10*B2)/B	#0.2467 mole fraction
X_Bb=(0.21*D2+0.01*B2)/B	#0.0767 mole fraction


# Print results
print('\nColumn 3 ')
print('B (mol/min): {:6.4f}'.format(B))
print('Xylene,  X_Bx: {:6.4f}'.format(X_Bx))
print('Styrene, X_Bs: {:6.4f}'.format(X_Bs))
print('Toluene, X_Bt: {:6.4f}'.format(X_Bt))
print('Benzene, X_Bb: {:6.4f}'.format(X_Bb))