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

#Base012=SOLPSPLOT('012',[65])
#Base025=SOLPSPLOT('025',[153])

BASEDRT, TOPDRT = SET_WDIR()

Probe012=loadmat('{}gfileProcessing/cmod_files/1160718012_ProbeData.mat'.format(TOPDRT))

Probe025=loadmat('{}gfileProcessing/cmod_files/1160718025_ProbeData.mat'.format(TOPDRT))

Fig, RadPlot = plt.subplots(nrows=2,ncols=2,sharex=True,sharey='row')

RadPlot[0,0].plot(Probe012['rho_ON'],Probe012['Ne_O'])
RadPlot[0,0].plot(Probe025['rho_ON'],Probe025['Ne_O'])
RadPlot[0,0].set_title('Outer Divertor Ne')
RadPlot[0,0].set_ylabel('Electron Density $m^{-3}$')

RadPlot[1,0].plot(Probe012['rho_OT'],Probe012['Te_O'])
RadPlot[1,0].plot(Probe025['rho_OT'],Probe025['Te_O'])
RadPlot[1,0].set_title('Outer Divertor Te')
RadPlot[1,0].set_ylabel('Electron Temperature $eV$')
RadPlot[1,0].set_xlabel('R-Rsep (mm)')

RadPlot[0,1].plot(Probe012['rho_IN'],Probe012['Ne_I'])
RadPlot[0,1].plot(Probe025['rho_IN'],Probe025['Ne_I'])
RadPlot[0,1].set_title('Inner Divertor Ne')

RadPlot[1,1].plot(Probe012['rho_IT'],Probe012['Te_I'])
RadPlot[1,1].plot(Probe025['rho_IT'],Probe025['Te_I'])
RadPlot[1,1].set_title('Inner Divertor Te')
RadPlot[1,1].set_xlabel('R-Rsep (mm)')
