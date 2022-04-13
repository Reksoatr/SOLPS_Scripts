# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 04:39:17 2019

@author: rmreksoatmodjo
"""

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import colors, cm, rc
from VesselPlotter import SOLPSPlotter

plt.rc('font',size=25)
plt.rc('lines',linewidth=5,markersize=10)

Attempt12 = 57
Attempt25 = 129
PARAM = 'Ne'
jxa = 55
sep = 20
LegLoc = 3
LS = 0
Rad = 'psin'

Shot012 = SOLPSPlotter(12, [Attempt12], PARAM, 'Export',PsinOffset=-0.01,LOG10=LS, RRad=Rad)
Shot025 = SOLPSPlotter(25, [Attempt25], PARAM, 'Export',PsinOffset=-0.015,LOG10=LS, RRad=Rad)

fig1 = plt.figure(figsize=(14,7))
cmap = cm.Set1.colors
plt.plot(Shot012['RR'].loc[:,jxa,Attempt12], Shot012['PARAM'].loc[:,jxa,Attempt12], color=cmap[0])
plt.plot(Shot025['RR'].loc[:,jxa,Attempt25],Shot025['PARAM'].loc[:,jxa,Attempt25], color=cmap[1])
plt.title('Radial Profile of ' + Shot012['PARAM'].name)
if PARAM == 'Ne':
    plt.errorbar(Shot012['Rexp'],Shot012['NemidAvg'],yerr=Shot012['ErrNe'],fmt='o',color=cmap[0],linewidth=3,capsize=7)
    plt.errorbar(Shot025['Rexp'],Shot025['NemidAvg'],yerr=Shot025['ErrNe'],fmt='o',color=cmap[1],linewidth=3,capsize=7)
if PARAM == 'Te':
    plt.errorbar(Shot012['Rexp'],Shot012['TemidAvg'],yerr=Shot012['ErrTe'],fmt='o',color=cmap[0],linewidth=3,capsize=7)
    plt.errorbar(Shot025['Rexp'],Shot025['TemidAvg'],yerr=Shot025['ErrTe'],fmt='o',color=cmap[1],linewidth=3,capsize=7)
a = plt.gca()
#a.set_xticklabels(['%.2f' % i for i in a.get_xticks()], fontsize=25)
#a.set_yticklabels(['%.1e' % j for j in a.get_yticks()], fontsize=25)
#plt.xlabel(Shot012['RR'].name,fontsize=25)
if LS == 1:
    plt.ylabel(r'$Log_{10}$ of ' + Shot012['PARAM'].name)
else:
    plt.ylabel(Shot012['PARAM'].name)
plt.legend(['High Current Ip','Low Current Ip'],loc=LegLoc)
Pmin = np.amin([a.get_yticks()[1],np.amin([Shot025['PARAM'].loc[:,jxa,Attempt25],Shot012['PARAM'].loc[:,jxa,Attempt12]])])
Pmax = np.amax([a.get_yticks()[-2],np.amax([Shot025['PARAM'].loc[:,jxa,Attempt25],Shot012['PARAM'].loc[:,jxa,Attempt12]])])
plt.plot([Shot025['RR'].loc[sep,jxa,Attempt25], Shot025['RR'].loc[sep,jxa,Attempt25]],[Pmin, Pmax],color='Black')
plt.grid(linestyle='-.')

'''
Shot012KYE = SOLPSPlotter(12, [Attempt12], 'KYE', 'Export')
Shot025KYE = SOLPSPlotter(25, [Attempt25], 'KYE', 'Export')

fig2 = plt.figure(figsize=(14,10))
plt.plot(Shot012KYE['RR'].loc[:,jxa,Attempt12],Shot012KYE['PARAM'].loc[:,jxa,Attempt12],'b',Shot025KYE['RR'].loc[:,jxa,Attempt25],Shot025KYE['PARAM'].loc[:,jxa,Attempt25],'g')
plt.title('Comparison of Thermal Diffusion Profiles of Discharge 012 vs Discharge 025')
plt.xlabel('Radial Coordinate')
plt.ylabel('Thermal Diffusion Coefficient')
plt.legend(['Shot012','Shot025'])
Pmin = float(Shot025KYE['PARAM'].values.min())
Pmax = float(Shot025KYE['PARAM'].values.max())
plt.plot([Shot025KYE['RR'].loc[sep,jxa,Attempt25], Shot025KYE['RR'].loc[sep,jxa,Attempt25]],[Pmin, Pmax],color='Black')
plt.grid()
'''
