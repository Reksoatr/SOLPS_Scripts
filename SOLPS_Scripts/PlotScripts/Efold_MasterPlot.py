# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 16:20:03 2020

@author: 18313
"""

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
import json
from TOOLS import SET_WDIR
import numpy as np
#from SOLPS_Plotter import SOLPSPLOT

BASEDRT, TOPDRT = SET_WDIR('','')

Shot025=[['21N','19N'],[0,0]]
Shot012=[['46N','48N','50N'],[0,0,0]]

#fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots(nrows=2,ncols=1,sharex=True)
plt.subplots_adjust(hspace=0.0)
#fig3, ax3 = plt.subplots()
fig4, ax4 = plt.subplots()

file='{}cmod/JerryV1.txt'.format(BASEDRT)
Jerry_LD=np.loadtxt(file,usecols=(0,1,2),skiprows=1)

file='{}cmod/JerryV2.txt'.format(BASEDRT)
Jerry_PedWid=np.loadtxt(file,usecols=(0,1,2),skiprows=1)

PolLegend=[]
JXA=55
JXI=37

plt.rc('axes',titlesize=30,labelsize=30)
plt.rc('lines',linewidth=5,markersize=20,markeredgewidth=2,linestyle='solid')
plt.rc('legend',fontsize=25,title_fontsize=25)
plt.rc('xtick',labelsize=25)
plt.rc('ytick',labelsize=25)
plt.rc('xtick.major',size=10,width=3)
plt.rc('ytick.major',size=10,width=3)

'''
NB_012={'x':[],'y':[],'NeuDen':[]}
SB_012={'x':[],'y':[],'NeuDen':[]}
NB_025={'x':[],'y':[],'NeuDen':[]}
SB_025={'x':[],'y':[],'NeuDen':[]}
'''

for N in range(len(Jerry_LD[:,0])):
    if Jerry_LD[N,0]==0.46:
        Marker='s'
        Color='grey'
    elif Jerry_LD[N,0]==0.76:
        Marker='D'
        Color='grey'
    elif Jerry_LD[N,0]==0.98:
        Marker='o'
        Color='grey'
        
    ax4.plot(Jerry_LD[N,1],1000*Jerry_LD[N,2],marker=Marker,markerfacecolor=Color,markeredgewidth=2,markeredgecolor='k',fillstyle='full')
    
for N in range(len(Jerry_PedWid[:,0])):
    if Jerry_PedWid[N,0]==0.46:
        Marker='s'
    elif Jerry_PedWid[N,0]==0.77:
        Marker='D'
    elif Jerry_PedWid[N,0]==0.98:
        Marker='o'
        
    ax4.plot(Jerry_PedWid[N,1],Jerry_PedWid[N,2],marker=Marker,markerfacecolor='w',markeredgewidth=2,markeredgecolor='k',fillstyle='none')
        
for N, Attempt in enumerate(Shot012[0]):
    if Shot012[1][N]==1:
        path='gaspuff/'
    else:
        path=''
    file='{}{}cmod/012home/Attempt{}/efold_data_{}.json'.format(BASEDRT,path,Attempt,Attempt)
    with open(file) as fp:
        data=json.load(fp)
        
    if data[0]['gaslvl']==8.25:
        Line_style='solid'
    elif data[0]['gaslvl']==77.8:
        Line_style='dashed'
    elif data[0]['gaslvl']==39.5:
        Line_style='dotted'
    
    efold_lengths=np.array(sorted(data[2].items())).astype(float)
    efold_err=np.array(sorted(data[5].items())).astype(float)[:,1]
    idx=efold_lengths[:,0].astype(int)
    neuden=np.array(data[4]['data'])[idx]
    pol_coords=np.array(data[3]['data'])[idx-24]
    
    '''
    X_PED=np.mean([data[0]['LFS_NePED'],data[0]['HFS_NePED']])
    Y_EFOLD=np.mean([data[0]['LFS'],data[0]['HFS']])
    NEUDEN=np.mean([data[0]['LFS_NeuDen'],data[0]['HFS_NeuDen']])
    
    if data[0]['balloon']==1:
        fill_style='none'
        SB_012['x'].append(X_PED)
        SB_012['y'].append(Y_EFOLD)
        SB_012['NeuDen'].append(NEUDEN)
    else:
        fill_style='full'
        NB_012['x'].append(X_PED)
        NB_012['y'].append(Y_EFOLD)
        NB_012['NeuDen'].append(NEUDEN)
    
    
    ax1.plot(data[0]['gaslvl'],data[0]['LFS'],fillstyle=fill_style,color='blue',marker='>',markersize=15)
    ax1.plot(data[0]['gaslvl'],data[0]['HFS'],fillstyle=fill_style,color='blue',marker='<',markersize=15)
    '''
    ax2[0].semilogy(pol_coords[:,0],neuden,linestyle=Line_style,color='blue') #marker='v',markersize=15)
    #ax2[1].plot(pol_coords[:,0],efold_lengths[:,1],linestyle=Line_style,color='blue') #marker='v',markersize=15)
    #ax2[1].fill_between(pol_coords[:,0],efold_lengths[:,1]-efold_err,efold_lengths[:,1]+efold_err,alpha=0.1,edgecolor='b',facecolor='c')

    ax2[1].errorbar(pol_coords[:,0],efold_lengths[:,1],efold_err,linestyle=Line_style,color='blue',capsize=5)
    
    PolLegend.append('{} TorrL'.format(data[0]['gaslvl']))
    '''
    ax3.plot(data[0]['LFS_NePED'],data[0]['LFS'],fillstyle=fill_style,color='blue',marker='>',markersize=15)
    ax3.plot(data[0]['HFS_NePED'],data[0]['HFS'],fillstyle=fill_style,color='blue',marker='<',markersize=15)

    ax3.plot([data[0]['LFS_NePED'],data[0]['HFS_NePED']],[data[0]['LFS'],data[0]['HFS']],color='blue',linestyle='dotted')
    '''
    ax4.plot(data[0]['LFS_NePED'],data[0]['LFS_PedWidth'],marker='^',markerfacecolor='lightblue',markeredgecolor='k')
    ax4.plot(data[0]['LFS_NePED'],data[0]['LFS Gradient_Scale_Length'],marker='^',markerfacecolor='b',markeredgecolor='k')
    ax4.plot(data[0]['LFS_NePED'],data[0]['LFS'],marker='^',markerfacecolor='navy',markeredgecolor='k')
    

for N, Attempt in enumerate(Shot025[0]):
    if Shot025[1][N]==1:
        path='gaspuff/'
    else:
        path=''
    file='{}{}cmod/025home/Attempt{}/efold_data_{}.json'.format(BASEDRT,path,Attempt,Attempt)
    with open(file) as fp:
        data=json.load(fp)
    
    if data[0]['gaslvl']==6.77:
        Line_style='solid'
    elif data[0]['gaslvl']==72.2:
        Line_style='dashed'
    
    efold_lengths=np.array(sorted(data[2].items())).astype(float)
    efold_err=np.array(sorted(data[5].items())).astype(float)[:,1]
    idx=efold_lengths[:,0].astype(int)
    neuden=np.array(data[4]['data'])[idx]
    pol_coords=np.array(data[3]['data'])[idx-24]
    
    '''    
    X_PED=np.mean([data[0]['LFS_NePED'],data[0]['HFS_NePED']])
    Y_EFOLD=np.mean([data[0]['LFS'],data[0]['HFS']])
    NEUDEN=np.mean([data[0]['LFS_NeuDen'],data[0]['HFS_NeuDen']])
    
    if data[0]['balloon']==1:
        fill_style='none'
        SB_025['x'].append(X_PED)
        SB_025['y'].append(Y_EFOLD)
        SB_025['NeuDen'].append(NEUDEN)
    else:
        fill_style='full'
        NB_025['x'].append(X_PED)
        NB_025['y'].append(Y_EFOLD)
        NB_025['NeuDen'].append(NEUDEN)
        
    ax1.plot(data[0]['gaslvl'],data[0]['LFS'],fillstyle=fill_style,color='red',marker='>',markersize=15)
    ax1.plot(data[0]['gaslvl'],data[0]['HFS'],fillstyle=fill_style,color='red',marker='<',markersize=15)
    '''
    
    ax2[0].semilogy(pol_coords[:,0],neuden,linestyle=Line_style,color='r') #marker='v',markersize=15)
    #ax2[1].plot(pol_coords[:,0],efold_lengths[:,1].astype(float),linestyle=Line_style,color='r') #marker='v',markersize=15)
    #ax2[1].fill_between(pol_coords[:,0],efold_lengths[:,1]-efold_err,efold_lengths[:,1]+efold_err,alpha=0.1,edgecolor='r',facecolor='m')
    
    ax2[1].errorbar(pol_coords[:,0],efold_lengths[:,1],efold_err,linestyle=Line_style,color='red',capsize=5)
    
    PolLegend.append('{} TorrL'.format(data[0]['gaslvl']))
    
    #ax3[0].plot(data[0]['LFS_NePED'],data[0]['LFS'],fillstyle=fill_style,color='red',marker='>',markersize=15)
    #ax3[1].plot(data[0]['HFS_NePED'],data[0]['HFS'],fillstyle=fill_style,color='red',marker='>',markersize=15)
    #ax3.plot([data[0]['LFS_NePED'],data[0]['HFS_NePED']],[data[0]['LFS'],data[0]['HFS']],color='red',linestyle='dotted')
    
    ax4.plot(data[0]['LFS_NePED'],data[0]['LFS_PedWidth'],marker='v',markerfacecolor='salmon',markeredgecolor='k')
    ax4.plot(data[0]['LFS_NePED'],data[0]['LFS Gradient_Scale_Length'],marker='v',markerfacecolor='r',markeredgecolor='k')
    ax4.plot(data[0]['LFS_NePED'],data[0]['LFS'],marker='v',markerfacecolor='maroon',markeredgecolor='k')


'''
ax1.set_xlabel('D2 gas puff strength (TorrL)')
ax1.set_ylabel('Neutral e-folding length (mm)')

ax3.plot(NB_012['x'],NB_012['y'],'b--')
ax3.plot(SB_012['x'],SB_012['y'],'b-')
ax3.plot(NB_025['x'],NB_025['y'],'r--')
ax3.plot(SB_025['x'],SB_025['y'],'r-')

ax3.set_xlabel(r'Pedestal electron density $n_{e,PED}\;(m^{-3})$',fontsize=20)
ax3.set_ylabel('Adjusted neutral e-folding length (mm)',fontsize=20)

ax4.plot(NB_012['x'],NB_012['NeuDen'],'b--')
ax4.plot(SB_012['x'],SB_012['NeuDen'],'b-')
ax4.plot(NB_025['x'],NB_025['NeuDen'],'r--')
ax4.plot(SB_025['x'],SB_025['NeuDen'],'r-')
'''

OMP=ax2[0].axvline(data[3]['data'][np.where(np.array(data[3]['coords']['Poloidal Index']['data'])==JXA)[0][0]][0],color='purple',linewidth=3,linestyle='-',label='Outer midplane')
IMP=ax2[0].axvline(data[3]['data'][np.where(np.array(data[3]['coords']['Poloidal Index']['data'])==JXI)[0][0]][0],color='orange',linewidth=3,linestyle='-',label='Inner midplane')
XPT=ax2[0].axvline(data[3]['data'][0][0],color='black',linewidth=2,linestyle='-',label='X-point')
ax2[0].axvline(data[3]['data'][-1][0],color='black',linestyle='-',linewidth=2)
ax2[0].set_ylabel(r'Neutral density ($m^{-3}$)')
ax2[0].legend(PolLegend,loc=(0.33,0.67),ncol=2,title='Low Ip (1.0MA)   High Ip (1.3MA)')

#secax=ax2[0].secondary_xaxis('top',functions=(lambda x: x / np.max(pol_coords[:,0]), lambda x: x * np.max(pol_coords[:,0])))
#secax.set_xlabel('Normalized poloidal distance from x-point along separatrix, clockwise',fontsize=20)

ax2[1].axvline(data[3]['data'][np.where(np.array(data[3]['coords']['Poloidal Index']['data'])==JXA)[0][0]][0],color='purple',linestyle='-',linewidth=3)
ax2[1].axvline(data[3]['data'][np.where(np.array(data[3]['coords']['Poloidal Index']['data'])==JXI)[0][0]][0],color='orange',linestyle='-',linewidth=3)
ax2[1].axvline(data[3]['data'][0][0],color='black',linestyle='-',linewidth=2)
ax2[1].axvline(data[3]['data'][-1][0],color='black',linestyle='-',linewidth=2)
ax2[1].set_xlabel(r'Poloidal distance from x-point along separatrix, clockwise (m)')
ax2[1].set_ylabel(r'Neutral e-folding length (mm)')
ax2[1].legend(handles=[XPT,IMP,OMP],loc=(0.53,0.01))

ax4.set_xlabel(r'Pedestal electron density $n_{e,PED}\;(m^{-3})$',fontsize=30)
ax4.set_ylabel(r'Length (mm)',fontsize=30)
xmax=ax4.get_xlim()[1]
ymax=ax4.get_ylim()[1]
ax4.set_xlim([0,xmax])
ax4.set_ylim([0,ymax])
### Create Legend Objects ###

PedWidth_46=mlines.Line2D([],[], marker='s',markerfacecolor='w',markeredgecolor='k',label='0.46 MA',linestyle='none')
PedWidth_76=mlines.Line2D([],[], marker='D',markerfacecolor='w',markeredgecolor='k',label='0.76 MA',linestyle='none')
PedWidth_98=mlines.Line2D([],[], marker='o',markerfacecolor='w',markeredgecolor='k',label='0.98 MA',linestyle='none')
PedWidth_012=mlines.Line2D([],[], marker='^',markerfacecolor='lightblue',markeredgecolor='k',label='1.3 MA',linestyle='none')
PedWidth_025=mlines.Line2D([],[], marker='v',markerfacecolor='salmon',markeredgecolor='k',label='1.0 MA',linestyle='none')

L_D_46=mlines.Line2D([],[], marker='s',markerfacecolor='grey',markeredgecolor='k',label='0.48 MA',linestyle='none')
L_D_76=mlines.Line2D([],[], marker='D',markerfacecolor='grey',markeredgecolor='k',label='0.76 MA',linestyle='none')
L_D_98=mlines.Line2D([],[], marker='o',markerfacecolor='grey',markeredgecolor='k',label='0.98 MA',linestyle='none')
L_D_012=mlines.Line2D([],[], marker='^',markerfacecolor='b',markeredgecolor='k',label='1.3 MA',linestyle='none')
L_D_025=mlines.Line2D([],[], marker='v',markerfacecolor='r',markeredgecolor='k',label='1.0 MA',linestyle='none')

EFold_Space=mlines.Line2D([],[], marker='o',markerfacecolor='w',markeredgecolor='w',linestyle='none')
EFold_012=mlines.Line2D([],[], marker='^',markerfacecolor='navy',markeredgecolor='k',label='1.3 MA',linestyle='none')
EFold_025=mlines.Line2D([],[], marker='v',markerfacecolor='maroon',markeredgecolor='k',label='1.0 MA',linestyle='none')

legend_cols= [(PedWidth_46,PedWidth_76,PedWidth_98,PedWidth_025,PedWidth_012),(L_D_46,L_D_76,L_D_98,L_D_025,L_D_012), (EFold_Space,EFold_Space,EFold_Space,EFold_025,EFold_012)]
legend_title='         0.46 MA   0.76 MA   0.98 MA   1.0 MA   1.3 MA'
legend_keys=[r'$\Delta n$',r'$<L_D>$',r'$\lambda$']
ax4.legend(legend_cols,legend_keys,title=legend_title,ncol=1,numpoints=1,markerfirst=False,shadow=True,handlelength=24,labelspacing=1,handler_map={tuple: HandlerTuple(ndivide=None)})

'''
Red = mpatches.Patch(color='red',label='Low Ip (1.0 MA)')
Blue = mpatches.Patch(color='blue', label='High Ip (1.3 MA)')
In_HFS = mlines.Line2D([],[], color='black',marker='<',fillstyle='full',linestyle='none',markersize=15,label='HFS')
Out_LFS = mlines.Line2D([],[], color='black',marker='>',fillstyle='full',linestyle='none',markersize=15,label='LFS')
NB = mlines.Line2D([],[], color='black',marker='>',fillstyle='full',linestyle='none',markersize=15,label='Ballooning Off')
SB = mlines.Line2D([],[], color='black',marker='>',fillstyle='none',linestyle='none',markersize=15,label='Ballooning On')
L_D=mlines.Line2D([],[],color='#1f77b4ff',marker='*',linestyle='none',markersize=15,label=r'$<L_D>$, various*')

ax3.legend(handles=[Red,Blue,Out_LFS,In_HFS],fontsize=20)
ax4.legend(handles=[Red,Blue,Out_LFS,In_HFS,NB,SB],fontsize=20)
'''