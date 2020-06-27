# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 16:20:03 2020

@author: 18313
"""

import matplotlib.pyplot as plt
import json
from TOOLS import SET_WDIR
#import numpy as np
#from VesselPlotterNew import SOLPSPLOT

BASEDRT, TOPDRT = SET_WDIR('','')

Shot025=[[65,69],[1,1]]
Shot012=[[203],[0]]

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

for N, Attempt in enumerate(Shot025[0]):
    if Shot025[1][N]==1:
        path='gaspuff/'
    else:
        path=''
    file='{}{}cmod/025home/Attempt{}/efold_data_{}.json'.format(BASEDRT,path,Attempt,Attempt)
    with open(file) as fp:
        data=json.load(fp)
    if data[0]['balloon']==1:
        fill_style='full'
    else:
        fill_style='none'
        
    ax1.plot(data[0]['gaslvl'],data[0]['LFS'],fillstyle=fill_style,color='red',marker='>',markersize=15)
    ax1.plot(data[0]['gaslvl'],data[0]['HFS'],fillstyle=fill_style,color='red',marker='<',markersize=15)

    ax2.semilogy(data[0]['gaslvl'],data[0]['LFS_NeuDen'],fillstyle=fill_style,color='red',marker='>',markersize=15)
    ax2.semilogy(data[0]['gaslvl'],data[0]['HFS_NeuDen'],fillstyle=fill_style,color='red',marker='<',markersize=15)

for N, Attempt in enumerate(Shot012[0]):
    if Shot012[1][N]==1:
        path='gaspuff/'
    else:
        path=''
    file='{}{}cmod/012home/Attempt{}/efold_data_{}.json'.format(BASEDRT,path,Attempt,Attempt)
    with open(file) as fp:
        data=json.load(fp)
    if data[0]['balloon']==1:
        fill_style='full'
    else:
        fill_style='none'
        
    ax1.plot(data[0]['gaslvl'],data[0]['LFS'],fillstyle=fill_style,color='blue',marker='>',markersize=15)
    ax1.plot(data[0]['gaslvl'],data[0]['HFS'],fillstyle=fill_style,color='blue',marker='<',markersize=15)

    ax2.semilogy(data[0]['gaslvl'],data[0]['LFS_NeuDen'],fillstyle=fill_style,color='blue',marker='>',markersize=15)
    ax2.semilogy(data[0]['gaslvl'],data[0]['HFS_NeuDen'],fillstyle=fill_style,color='blue',marker='<',markersize=15)

ax1.set_xlabel('D2 gas puff strength (TorrL)')
ax1.set_ylabel('Neutral e-folding length (mm)')

ax2.set_xlabel('D2 gas puff strength (TorrL)')
ax2.set_ylabel(r'Neutral Density at Separatrix ($m^{-3}$)')