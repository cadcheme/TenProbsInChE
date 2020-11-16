# -*- coding: utf-8 -*-
"""
Problem 1.
Use of the van der Waals equation of state to calculate molar volume and 
compressibility factor for a gas.

The problem has 3 parts: 
a. Calculate the molar volume and compressibility factor for gaseous ammonia 
   at a pressure P = 56 atm and a temperature T = 450 K using the van der Waals
   equation of state.
b. Repeat the calculations for the following reduced pressures: 
   Pr = 1, 2, 4, 10, and 20.
c. How does the compressibility factor vary as a function of Pr.?


P = pressure in atm                                                
V = molar volume in liters/g-mol                                   
T = temperature in K                                               
R = gas constant (R = 0.08206 atm.liter/g-mol.K)                   
Tc = critical temperature (405.5 K for ammonia)                    
Pc = critical pressure (111.3 atm for ammonia)
Pr = reduced pressure
Tr = Reduced temperature
"""


# -----------------------------------------------------------------------------
# Mohammad Rahmani
# Chemical Engineering Department
# Amirkabir University of Technology
# Tehran, Iran
# m.rahmani@aut.ac.ir
#
# Rev 0.2
# Nov 15, 2020



from   scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt


# Problem data and parameters
R = 0.08206 #general gas constant, atm.lit/gmol.K

Tc = 405.5 #ammonia critical temp, K
Pc = 111.3 #ammonia critical pres, atm

# van der Waals EOS parameters
a = ( 27. / 64. ) * ( ( ( R**2 ) * ( Tc**2 ) ) / Pc )
b = ( R * Tc ) / ( 8. * Pc )

# van der Waals EOS
def van_der_waals(V):
    f=(P+a/(V**2))*(V-b)-R*T
    return f


# Part a
#-----------------------------------------------------------------------------
P = 56   # Pres, atm
T = 450  # Temp, K

# Initial guess for Vm
V0 = P/(R*T)

# Solve Cubic van der Waals EOS
Vm = fsolve(van_der_waals,V0)
Z = (P*Vm)/(R*T) 

print("Part a\n") #empty line
print("Molar volume, lit/gmol", Vm[0])
print("Compresibilty factor  ", Z[0])
print("\n") #empty line

# Part b
#-----------------------------------------------------------------------------
Pr = np.array([1., 2., 4., 10., 20.])   
n  =np.size(Pr)
Vm = np.zeros(n) 
Z  = np.zeros(n) 

for i in range(n):
    P=Pr[i]*Pc
    
    # Initial guess for Vm
    V0 = P/(R*T)
    # Solve Cubic van der Waals EOS
    Vm[i] = fsolve(van_der_waals,V0)
    Z[i] = (P*Vm[i])/(R*T) 

np.set_printoptions(precision=4) # for pretty print
print("Part b\n") #empty line
print("Pr:             ", Pr)
print("Vm, (lit/gmol): ", Vm)
print("Z:              ", Z)
print("\n") #empty line
    
# Part c
#-----------------------------------------------------------------------------
plt.plot(Pr,Z,'ro-')
plt.title ('Compressibility factor')
plt.xlabel ('Reduced pressure (Pr)')
plt.ylabel ('Compressibility factor, Z')
