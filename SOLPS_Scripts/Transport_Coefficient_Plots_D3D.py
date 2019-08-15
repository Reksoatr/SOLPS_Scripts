# -*- coding: utf-8 -*-
"""
Created on Thu May 30 13:18:13 2019

@author: rmreksoatmodjo
"""
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from VesselPlotter import SOLPSPlotter
import D3DPreProcess

plt.rc('font',size=25)
plt.rc('lines',linewidth=5,markersize=10)

Shot='d3d'
Attempt=72
JXA=56
SEP = 18
TC_Psin = TC_psin
TC_Flux = TC_pflux

jxa = JXA-1
SEP = SEP-1

D = SOLPSPlotter(Shot,[Attempt],'DN','Export',SEP=21)
Chi_e = SOLPSPlotter(Shot,[Attempt],'KYE','Export',SEP=21)
Chi_i = SOLPSPlotter(Shot,[Attempt],'KYI','Export',SEP=21)
IonFlx = SOLPSPlotter(Shot,[Attempt],'IonFlx','Export',SEP=21)
#RadFlx = SOLPSPlotter(Shot,[Attempt],'RadFlx','Export',SEP=21)
#MolFlx = SOLPSPlotter(Shot,[Attempt],'MolFlx','Export',SEP=21)
NeuDen = SOLPSPlotter(Shot,[Attempt],'NeuDen','Export',SEP=21)
VOL = SOLPSPlotter(Shot,[Attempt],'VOL','Export',SEP=21)
SY = SOLPSPlotter(Shot,[Attempt],'SY','Export',SEP=21)

'''
fig1 = plt.figure(figsize=(14,7))

Q1 = Chi_e['PARAM'].loc[:,jxa,Attempt]/Chi_i['PARAM'].loc[:,jxa,Attempt]

plt.plot(D['RR'].loc[:,jxa,Attempt], Q1)

plt.xlabel('$\psi_N$')
plt.ylabel('$\chi_e$ / $\chi_i$ ($m^2/s$)')
plt.title('Midplane Radial Transport Coefficient $\chi_e$ / $\chi_i$')
Pmin = float(Q1.min())
Pmax = float(Q1.max())
plt.plot([D['RR'].loc[SEP,jxa,Attempt], D['RR'].loc[SEP,jxa,Attempt]],[Pmin, Pmax],color='Black',linewidth=3)


fig2 = plt.figure(figsize=(14,7))

Q2 = D['PARAM'].loc[:,jxa,Attempt]/Chi_e['PARAM'].loc[:,jxa,Attempt]

plt.plot(D['RR'].loc[:,jxa,Attempt], Q2)

plt.xlabel('$\psi_N$')
plt.ylabel('$D$ / $\chi_e$ ($m^2/s$)')
plt.title('Midplane Radial Transport Coefficient $D$ / $\chi_e$')
Pmin = float(Q2.min())
Pmax = float(Q2.max())
plt.plot([D['RR'].loc[SEP,jxa,Attempt], D['RR'].loc[SEP,jxa,Attempt]],[Pmin, Pmax],color='Black',linewidth=3)


fig3 = plt.figure(figsize=(14,7))

Q3 = D['PARAM'].loc[:,jxa,Attempt]/(Chi_e['PARAM'].loc[:,jxa,Attempt] + Chi_i['PARAM'].loc[:,jxa,Attempt])

plt.plot(D['RR'].loc[:,jxa,Attempt], Q3)

plt.xlabel('$\psi_N$')
plt.ylabel('$D$ / ($\chi_e+\chi_i$) ($m^2/s$)')
plt.title('Midplane Radial Transport Coefficient $D$ / ($\chi_e+\chi_i$)')
Pmin = float(Q3.min())
Pmax = float(Q3.max())
plt.plot([D['RR'].loc[SEP,jxa,Attempt], D['RR'].loc[SEP,jxa,Attempt]],[Pmin, Pmax],color='Black',linewidth=3)
'''
fig4 = plt.figure(figsize=(14,7))
Q4 = np.sum(IonFlx['PARAM'].values[:,:,0], axis=1) / np.sum(SY['PARAM'].values[:,:,0], axis=1)
#Q5 = IonFlx['PARAM'].loc[:,jxa,Attempt] / SY['PARAM'].loc[:,jxa,Attempt]

plt.plot(IonFlx['RR'].loc[:,jxa,Attempt], Q4,'x')
plt.plot(TC_Psin,TC_Flux)

plt.legend(['SOLPS','ONETWO'])
plt.xlabel('$\psi_N$')
plt.ylabel('Radial Particle Flux ($m^{-2}s^{-1}$)')
plt.title('Radial Particle Flux')
Pmin = np.amin([Q4])
Pmax = np.amax([Q4])
plt.plot([IonFlx['RR'].loc[SEP,jxa,Attempt], IonFlx['RR'].loc[SEP,jxa,Attempt]],[Pmin, Pmax],color='Black',linewidth=3)


fig5 = plt.figure(figsize=(14,7))

Q6 = np.sum((NeuDen['PARAM'].values[:,:,0]*VOL['PARAM'].values[:,:,0]), axis=1) / np.sum(VOL['PARAM'].values[:,:,0], axis=1)
plt.semilogy(NeuDen['RR'].loc[:,jxa,Attempt], Q6)

plt.xlabel('$\psi_N$')
plt.ylabel('Average Neutral Density ($m^{-3}$)')
plt.title('Volume-Averaged Neutral Density')
Pmin = float(Q6.min())
Pmax = float(Q6.max())
plt.plot([NeuDen['RR'].loc[SEP,jxa,Attempt], NeuDen['RR'].loc[SEP,jxa,Attempt]],[Pmin, Pmax],color='Black',linewidth=3)
