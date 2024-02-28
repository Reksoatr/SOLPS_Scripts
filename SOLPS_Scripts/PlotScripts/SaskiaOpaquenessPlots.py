# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 13:23:09 2024

@author: 18313
"""

import numpy as np
import matplotlib.pyplot as plt
from TOOLS import SET_WDIR
import json
import os

BASEDRT, TOPDRT = SET_WDIR('','')

### Input Fields ###

#List parameters for input shots sequentially in different lists

shots=['d3d','012','025']

attempts=['86','46N','21N']

labels=['DIII-D 175060','C-Mod 1160718012', 'C-Mod 1160718025']

colors=['red','green','blue']

dataDict={}

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

for i, s in enumerate(shots):
    if s == 'd3d':
        file='{}/{}/Attempt{}/efold_data_{}.json'.format(BASEDRT,s,attempts[i],attempts[i])

    else:
        file='{}/cmod/{}home/Attempt{}/efold_data_{}.json'.format(BASEDRT,s,attempts[i],attempts[i])
    
    with open(file) as f:
        dataDict[s]=json.load(f)
         
    efold_adj = np.array(list(dataDict[s][2]['flux expansion adjusted efolding length (mm)'].values()))
    delta_ne = np.array(list(dataDict[s][3]['n_e pedestal width (mm)'].values()))
    opaqueness = np.array(list(dataDict[s][4]['opaqueness (delta_ne/efold_adj)'].values())) 

    X_dict = dict(zip(dataDict['d3d'][5]['coords']['Poloidal Index']['data'],np.array(dataDict[s][5]['data'])[:,1]))
    X = np.array([X_dict[float(k)] for k in dataDict[s][2]['flux expansion adjusted efolding length (mm)'].keys()])

    ax1.plot(X,efold_adj,'v:',color=colors[i],label='{} adjusted e-folding length'.format(labels[i]))
    ax1.plot(X,delta_ne, '^:',color=colors[i],label='{} pedestal width'.format(labels[i]))
    ax1.set_xlabel('Normalized poloidal distance from X-point, clockwise along separatrix')
    ax1.set_ylabel('Length (mm)')
    
    ax2.plot(X,opaqueness,'*:',color=colors[i],label='{} opaqueness'.format(labels[i]))
    ax2.set_xlabel('Normalized poloidal distance from X-point, clockwise along separatrix')
    ax2.set_ylabel('Neutral Opaqueness')

ax1.axvline(x=0.25,color='orange',label='Inner Midplane')
ax1.axvline(x=0.7,color='red',label='Outer Midplane')

ax2.axvline(x=0.25,color='orange',label='Inner Midplane')
ax2.axvline(x=0.7,color='red',label='Outer Midplane')

ax1.legend()    
ax2.legend()