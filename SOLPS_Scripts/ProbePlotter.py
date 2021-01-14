# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 12:59:44 2019

@author: rmreksoatmodjo
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from VesselPlotterNew import SOLPSPLOT
from TOOLS import SET_WDIR

Shot012=1
Shot025=1

Attempt012='46N'
Attempt025='21N'
GasPuffTxt='Low'

#Good X-Limits for this data seem to be [-0.01,0.023]

plt.rc('axes',titlesize=30,labelsize=30)
plt.rc('lines',linewidth=5,markersize=20,markeredgewidth=2,linestyle='solid')
plt.rc('legend',fontsize=25,title_fontsize=25)
plt.rc('xtick',labelsize=25)
plt.rc('ytick',labelsize=25)
plt.rc('xtick.major',size=10,width=3)
plt.rc('ytick.major',size=10,width=3)

BASEDRT, TOPDRT = SET_WDIR('','')

if Shot012 == 1:
    
    COLOR='navy'
    
    Base012=SOLPSPLOT('12',[Attempt012])
    
    Probe012=loadmat('{}gfileProcessing/cmod_files/1160718012_ProbeData.mat'.format(TOPDRT))
    
    AVG = np.mean(Probe012['Ne_O'],axis=0)
    
    Fig1, Rad012Plot = plt.subplots(nrows=2,ncols=2,sharex=True,sharey='row')
    
    Base012.RadProf('Ne',JXA=96,RADC='rrsep',AX=Rad012Plot[0,1],PlotScheme=['b'],Markers=False,LOG10=2)
    Rad012Plot[0,1].semilogy(Probe012['rho_ON']/1000,Probe012['Ne_O'],'.',markersize=3,color=COLOR)
    Rad012Plot[0,1].set_title('High Density, {} Gas Puff Outer Divertor'.format(GasPuffTxt))
    Rad012Plot[0,1].set_ylabel('')
    Rad012Plot[0,1].set_xlabel('')
    Rad012Plot[0,1].get_legend().remove()
    
    Base012.RadProf('Te',JXA=96,RADC='rrsep',AX=Rad012Plot[1,1],PlotScheme=['b'],Markers=False)
    Rad012Plot[1,1].plot(Probe012['rho_OT']/1000,Probe012['Te_O'],'.',markersize=3,color=COLOR)
    Rad012Plot[1,1].set_title('')
    Rad012Plot[1,1].set_ylabel('')
    Rad012Plot[1,1].set_xlabel('R-Rsep (m)')
    Rad012Plot[1,1].get_legend().remove()
    
    Base012.RadProf('Ne',JXA=1,RADC='rrsep',AX=Rad012Plot[0,0],PlotScheme=['b'],Markers=False,LOG10=2)
    Rad012Plot[0,0].semilogy(Probe012['rho_IN']/1000,Probe012['Ne_I'],'.',markersize=3,color=COLOR)
    Rad012Plot[0,0].set_title('High Density, {} Gas Puff Inner Divertor'.format(GasPuffTxt))
    Rad012Plot[0,0].set_xlabel('')
    Rad012Plot[0,0].set_ylabel('$n_e\;(m^{-3})$')
    #Rad012Plot[0,0].get_legend().remove()
    
    Base012.RadProf('Te',JXA=1,RADC='rrsep',AX=Rad012Plot[1,0],PlotScheme=['b'],Markers=False)
    Rad012Plot[1,0].plot(Probe012['rho_IT']/1000,Probe012['Te_I'],'.',markersize=3,color=COLOR)
    Rad012Plot[1,0].set_title('')
    Rad012Plot[1,0].set_xlabel('R-Rsep (m)')
    Rad012Plot[1,0].set_ylabel('$T_e\;(eV)$')
    Rad012Plot[1,0].get_legend().remove()
    
    plt.subplots_adjust(wspace=0, hspace=0)

if Shot025 == 1:

    COLOR='maroon'    

    Base025=SOLPSPLOT('25',[Attempt025])
    
    Probe025=loadmat('{}gfileProcessing/cmod_files/1160718025_ProbeData.mat'.format(TOPDRT))
    
    Fig2, Rad025Plot = plt.subplots(nrows=2,ncols=2,sharex=True,sharey='row')
    
    Base025.RadProf('Ne',JXA=96,RADC='rrsep',AX=Rad025Plot[0,1],PlotScheme=['r'],Markers=False,LOG10=2)
    Rad025Plot[0,1].semilogy(Probe025['rho_ON']/1000,Probe025['Ne_O'],'.', markersize=3,color=COLOR)
    Rad025Plot[0,1].set_title('Low Density, {} Gas Puff Outer Divertor'.format(GasPuffTxt))
    Rad025Plot[0,1].set_ylabel('')
    Rad025Plot[0,1].set_xlabel('')
    Rad025Plot[0,1].get_legend().remove()
    
    Base025.RadProf('Te',JXA=96,RADC='rrsep',AX=Rad025Plot[1,1],PlotScheme=['r'],Markers=False)
    Rad025Plot[1,1].plot(Probe025['rho_OT']/1000,Probe025['Te_O'],'.',markersize=3,color=COLOR)
    Rad025Plot[1,1].set_title('')
    Rad025Plot[1,1].set_ylabel('')
    Rad025Plot[1,1].set_xlabel('R-Rsep (m)')
    Rad025Plot[1,1].get_legend().remove()
    
    Base025.RadProf('Ne',JXA=1,RADC='rrsep',AX=Rad025Plot[0,0],PlotScheme=['r'],Markers=False,LOG10=2)
    Rad025Plot[0,0].semilogy(Probe025['rho_IN']/1000,Probe025['Ne_I'],'.',markersize=3,color=COLOR)
    Rad025Plot[0,0].set_title('Low Density, {} Gas Puff Inner Divertor'.format(GasPuffTxt))
    Rad025Plot[0,0].set_xlabel('')
    #Rad025Plot[0,0].get_legend().remove()
    
    Base025.RadProf('Te',JXA=1,RADC='rrsep',AX=Rad025Plot[1,0],PlotScheme=['r'],Markers=False)
    Rad025Plot[1,0].plot(Probe025['rho_IT']/1000,Probe025['Te_I'],'.',markersize=3,color=COLOR)
    Rad025Plot[1,0].set_title('')
    Rad025Plot[1,0].set_xlabel('R-Rsep (m)')
    Rad025Plot[1,0].get_legend().remove()
    
    plt.subplots_adjust(wspace=0, hspace=0)

