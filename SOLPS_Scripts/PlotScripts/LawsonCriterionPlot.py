# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 12:54:26 2023

@author: Richard Reksoatmodjo

Lawson Criterion plot
"""
import numpy as np
import matplotlib.pyplot as plt

plt.rc('font',size=20)
plt.rc('lines',linewidth=5)

Te=np.array([1,2,5,10,20,50,100,200,500,1000])

D_D=np.array([1.5E-22,5.4E-21,1.8E-19,1.2E-18,5.2E-18,2.1E-17,4.5E-17,8.8E-17,1.8E-16,2.2E-16])

D_T=np.array([5.5e-21,2.6e-19,1.3e-17,1.1e-16,4.2e-16,8.7e-16,8.5e-16,6.3e-16,3.7e-16,2.7e-16])

D_He3=np.array([1e-26,1.4e-23,6.7e-21,2.3e-19,3.8e-18,5.4e-17,1.6e-16,2.4e-16,2.3e-16,1.8e-16])

Ech_D_D = (4.85E3)/2

Ech_D_T = 3.5E3

Ech_D_He3 = 18.3E3

L_D_D=(12E6*Te)/(Ech_D_D*D_D)

L_D_T=(12E6*Te)/(Ech_D_T*D_T)

L_D_He3=(12E6*Te)/(Ech_D_He3*D_He3)

figL, axL = plt.subplots()

axL.loglog(Te, L_D_D, 'g-', label='$D-D$')
axL.loglog(Te, L_D_T, 'b-', label='$D-T$')
axL.loglog(Te, L_D_He3, 'r-', label='$D-He^3$')

ax2=axL.secondary_xaxis('top', functions=(lambda x: 11.6*x, lambda x : x/11.6))
ax2.set_xlabel('Temperature (million Kelvin)')

axL.set_xlabel('Temperature (keV)')
axL.set_ylabel('$n_e\tau_e\;(m^{-3}*s)$')
axL.set_title('Lawson Criterion for select fusion reactions', pad=25)
axL.legend()
