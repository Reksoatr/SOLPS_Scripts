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

Base012=SOLPSPLOT('12',[97])
#Base025=SOLPSPLOT('025',[154])

BASEDRT, TOPDRT = SET_WDIR('','')

Probe012=loadmat('{}gfileProcessing/cmod_files/1160718012_ProbeData.mat'.format(TOPDRT))
#Probe025=loadmat('{}gfileProcessing/cmod_files/1160718025_ProbeData.mat'.format(TOPDRT))

AVG = np.mean(Probe012['Ne_O'],axis=0)

Fig1, Rad012Plot = plt.subplots(nrows=2,ncols=2,sharex=True,sharey='row')
#Fig2, Rad025Plot = plt.subplots(nrows=2,ncols=2,sharex=True,sharey='row')

Base012.RadProf('Ne',JXA=96,RADC='rrsep',AX=Rad012Plot[0,1],PlotScheme=['r'],Markers=False)
Rad012Plot[0,1].plot(Probe012['rho_ON']/1000,Probe012['Ne_O'],'.',color='orange',markersize=1)
Rad012Plot[0,1].set_title('High Ip Outer Divertor Ne')
Rad012Plot[0,1].set_ylabel('Electron Density $m^{-3}$')
Rad012Plot[0,1].get_legend().remove()
'''
Base025.RadProf('Ne',JXA=96,RADC='rrsep',AX=Rad025Plot[0,1],PlotScheme=['b'],Markers=False)
Rad025Plot[0,1].plot(Probe025['rho_ON']/1000,Probe025['Ne_O'],'.', markersize=1)
Rad025Plot[0,1].set_title('Low Ip Outer Divertor Ne')
Rad025Plot[0,1].set_ylabel('Electron Density $m^{-3}$')
Rad025Plot[0,1].get_legend().remove()
'''
Base012.RadProf('Te',JXA=96,RADC='rrsep',AX=Rad012Plot[1,1],PlotScheme=['r'],Markers=False)
Rad012Plot[1,1].plot(Probe012['rho_OT']/1000,Probe012['Te_O'],'.',color='orange',markersize=1)
Rad012Plot[1,1].set_title('High Ip Outer Divertor Te')
Rad012Plot[1,1].set_ylabel('Electron Temperature $eV$')
Rad012Plot[1,1].set_xlabel('R-Rsep (m)')
Rad012Plot[1,1].get_legend().remove()
'''
Base025.RadProf('Te',JXA=96,RADC='rrsep',AX=Rad025Plot[1,1],PlotScheme=['b'],Markers=False)
Rad025Plot[1,1].plot(Probe025['rho_OT']/1000,Probe025['Te_O'],'.',markersize=1)
Rad025Plot[1,1].set_title('Low Ip Outer Divertor Te')
Rad025Plot[1,1].set_ylabel('Electron Temperature $eV$')
Rad025Plot[1,1].set_xlabel('R-Rsep (m)')
Rad025Plot[1,1].get_legend().remove()
'''
Base012.RadProf('Ne',JXA=1,RADC='rrsep',AX=Rad012Plot[0,0],PlotScheme=['r'],Markers=False)
Rad012Plot[0,0].plot(Probe012['rho_IN']/1000,Probe012['Ne_I'],'.',color='orange',markersize=1)
Rad012Plot[0,0].set_title('High Ip Inner Divertor Ne')
Rad012Plot[0,0].get_legend().remove()
'''
Base025.RadProf('Ne',JXA=1,RADC='rrsep',AX=Rad025Plot[0,0],PlotScheme=['b'],Markers=False)
Rad025Plot[0,0].plot(Probe025['rho_IN']/1000,Probe025['Ne_I'],'.',markersize=1)
Rad025Plot[0,0].set_title('Low Ip Inner Divertor Ne')
Rad025Plot[0,0].get_legend().remove()
'''
Base012.RadProf('Te',JXA=1,RADC='rrsep',AX=Rad012Plot[1,0],PlotScheme=['r'],Markers=False)
Rad012Plot[1,0].plot(Probe012['rho_IT']/1000,Probe012['Te_I'],'.',color='orange',markersize=1)
Rad012Plot[1,0].set_title('High Ip Inner Divertor Te')
Rad012Plot[1,0].set_xlabel('R-Rsep (m)')
Rad012Plot[1,0].get_legend().remove()
'''
Base025.RadProf('Te',JXA=1,RADC='rrsep',AX=Rad025Plot[1,0],PlotScheme=['b'],Markers=False)
Rad025Plot[1,0].plot(Probe025['rho_IT']/1000,Probe025['Te_I'],'.',markersize=1)
Rad025Plot[1,0].set_title('Low Ip Inner Divertor Te')
Rad025Plot[1,0].set_xlabel('R-Rsep (m)')
Rad025Plot[1,0].get_legend().remove()
'''