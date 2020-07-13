# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 16:20:03 2020

@author: 18313
"""

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import json
from TOOLS import SET_WDIR
import numpy as np
#import numpy as np
#from VesselPlotterNew import SOLPSPLOT

BASEDRT, TOPDRT = SET_WDIR('','')

Shot025=[['11N','12N','15N','16N','13N','14N'],[0,0,0,0,0,0]]
Shot012=[['36N','37N','38N','39N','40N','41N','42N','43N'],[0,0,0,0,0,0,0,0]]

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()
fig4, ax4 = plt.subplots()

NB_012={'x':[],'y':[],'NeuDen':[]}
SB_012={'x':[],'y':[],'NeuDen':[]}
NB_025={'x':[],'y':[],'NeuDen':[]}
SB_025={'x':[],'y':[],'NeuDen':[]}

for N, Attempt in enumerate(Shot025[0]):
    if Shot025[1][N]==1:
        path='gaspuff/'
    else:
        path=''
    file='{}{}cmod/025home/Attempt{}/efold_data_{}.json'.format(BASEDRT,path,Attempt,Attempt)
    with open(file) as fp:
        data=json.load(fp)
    X_PED=np.mean([data[0]['LFS_NePED'],data[0]['HFS_NePED']])
    Y_EFOLD=np.mean([data[0]['LFS'],data[0]['HFS']])
    NEUDEN=np.mean([data[0]['LFS_NeuDen'],data[0]['HFS_NeuDen']])
    
    if data[0]['balloon']==1:
        fill_style='full'
        SB_025['x'].append(X_PED)
        SB_025['y'].append(Y_EFOLD)
        SB_025['NeuDen'].append(NEUDEN)
    else:
        fill_style='none'
        NB_025['x'].append(X_PED)
        NB_025['y'].append(Y_EFOLD)
        NB_025['NeuDen'].append(NEUDEN)
        
    ax1.plot(data[0]['gaslvl'],data[0]['LFS'],fillstyle=fill_style,color='red',marker='>',markersize=15)
    ax1.plot(data[0]['gaslvl'],data[0]['HFS'],fillstyle=fill_style,color='red',marker='<',markersize=15)

    ax2.semilogy(data[0]['gaslvl'],data[0]['LFS_NeuDen'],fillstyle=fill_style,color='red',marker='>',markersize=15)
    ax2.semilogy(data[0]['gaslvl'],data[0]['HFS_NeuDen'],fillstyle=fill_style,color='red',marker='<',markersize=15)

    ax3.plot(data[0]['LFS_NePED'],data[0]['LFS'],fillstyle=fill_style,color='red',marker='>',markersize=15)
    ax3.plot(data[0]['HFS_NePED'],data[0]['HFS'],fillstyle=fill_style,color='red',marker='<',markersize=15)

    ax3.plot([data[0]['LFS_NePED'],data[0]['HFS_NePED']],[data[0]['LFS'],data[0]['HFS']],color='red',linestyle='dotted')

    ax4.plot(data[0]['LFS_NePED'],data[0]['LFS_NeuDen'],fillstyle=fill_style,color='red',marker='>',markersize=15)
    ax4.plot(data[0]['HFS_NePED'],data[0]['HFS_NeuDen'],fillstyle=fill_style,color='red',marker='<',markersize=15)

    ax4.plot([data[0]['LFS_NePED'],data[0]['HFS_NePED']],[data[0]['LFS_NeuDen'],data[0]['HFS_NeuDen']],color='red',linestyle='dotted')


for N, Attempt in enumerate(Shot012[0]):
    if Shot012[1][N]==1:
        path='gaspuff/'
    else:
        path=''
    file='{}{}cmod/012home/Attempt{}/efold_data_{}.json'.format(BASEDRT,path,Attempt,Attempt)
    with open(file) as fp:
        data=json.load(fp)
    X_PED=np.mean([data[0]['LFS_NePED'],data[0]['HFS_NePED']])
    Y_EFOLD=np.mean([data[0]['LFS'],data[0]['HFS']])
    NEUDEN=np.mean([data[0]['LFS_NeuDen'],data[0]['HFS_NeuDen']])
    
    if data[0]['balloon']==1:
        fill_style='full'
        SB_012['x'].append(X_PED)
        SB_012['y'].append(Y_EFOLD)
        SB_012['NeuDen'].append(NEUDEN)
    else:
        fill_style='none'
        NB_012['x'].append(X_PED)
        NB_012['y'].append(Y_EFOLD)
        NB_012['NeuDen'].append(NEUDEN)
        
    ax1.plot(data[0]['gaslvl'],data[0]['LFS'],fillstyle=fill_style,color='blue',marker='>',markersize=15)
    ax1.plot(data[0]['gaslvl'],data[0]['HFS'],fillstyle=fill_style,color='blue',marker='<',markersize=15)

    ax2.semilogy(data[0]['gaslvl'],data[0]['LFS_NeuDen'],fillstyle=fill_style,color='blue',marker='>',markersize=15)
    ax2.semilogy(data[0]['gaslvl'],data[0]['HFS_NeuDen'],fillstyle=fill_style,color='blue',marker='<',markersize=15)

    ax3.plot(data[0]['LFS_NePED'],data[0]['LFS'],fillstyle=fill_style,color='blue',marker='>',markersize=15)
    ax3.plot(data[0]['HFS_NePED'],data[0]['HFS'],fillstyle=fill_style,color='blue',marker='<',markersize=15)

    ax3.plot([data[0]['LFS_NePED'],data[0]['HFS_NePED']],[data[0]['LFS'],data[0]['HFS']],color='blue',linestyle='dotted')

    ax4.semilogy(data[0]['LFS_NePED'],data[0]['LFS_NeuDen'],fillstyle=fill_style,color='blue',marker='>',markersize=15)
    ax4.semilogy(data[0]['HFS_NePED'],data[0]['HFS_NeuDen'],fillstyle=fill_style,color='blue',marker='<',markersize=15)

    ax4.semilogy([data[0]['LFS_NePED'],data[0]['HFS_NePED']],[data[0]['LFS_NeuDen'],data[0]['HFS_NeuDen']],color='blue',linestyle='dotted')


ax1.set_xlabel('D2 gas puff strength (TorrL)')
ax1.set_ylabel('Neutral e-folding length (mm)')

ax2.set_xlabel('D2 gas puff strength (TorrL)')
ax2.set_ylabel(r'Neutral Density at Separatrix ($m^{-3}$)')

ax3.plot(NB_012['x'],NB_012['y'],'b--')
ax3.plot(SB_012['x'],SB_012['y'],'b-')
ax3.plot(NB_025['x'],NB_025['y'],'r--')
ax3.plot(SB_025['x'],SB_025['y'],'r-')

ax3.set_xlabel(r'Pedestal electron density $n_{e,PED}\;(m^{-3})$')
ax3.set_ylabel('Adjusted neutral e-folding length (mm)')

ax4.semilogy(NB_012['x'],NB_012['NeuDen'],'b--')
ax4.semilogy(SB_012['x'],SB_012['NeuDen'],'b-')
ax4.semilogy(NB_025['x'],NB_025['NeuDen'],'r--')
ax4.semilogy(SB_025['x'],SB_025['NeuDen'],'r-')

ax4.set_xlabel(r'Pedestal electron density $n_{e,PED}\;(m^{-3})$')
ax4.set_ylabel(r'Neutral density at separatrix ($m^{-3}$)')

### Create Legend Objects ###

Red = mpatches.Patch(color='red',label='Low Ip')
Blue = mpatches.Patch(color='blue', label='High Ip')
In_HFS = mlines.Line2D([],[], color='black',marker='<',fillstyle='full',linestyle='none',markersize=15,label='HFS')
Out_LFS = mlines.Line2D([],[], color='black',marker='>',fillstyle='full',linestyle='none',markersize=15,label='LFS')
NB = mlines.Line2D([],[], color='black',marker='>',fillstyle='none',linestyle='none',markersize=15,label='Ballooning Off')
SB = mlines.Line2D([],[], color='black',marker='>',fillstyle='full',linestyle='none',markersize=15,label='Ballooning On')

ax3.legend(handles=[Red,Blue,Out_LFS,In_HFS,NB,SB])
ax4.legend(handles=[Red,Blue,Out_LFS,In_HFS,NB,SB])