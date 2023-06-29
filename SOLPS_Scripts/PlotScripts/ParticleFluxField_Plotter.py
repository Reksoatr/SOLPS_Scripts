# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 14:00:45 2023

@author: Richard Reksoatmodjo

Particle Flux Field Plotter

Experiments with Matplotlib's built-in streamplot() and quiver() methods with
SOLPS-ITER particle flux simulation data
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from SOLPS_Plotter import SOLPSPLOT

StreamPlot = False
Quiver = True

Shot = '1100305023'
Attempt = '24Rf2.0_split2'
JXA = 55
JXI = 38
SEP = 18
LINTHRESH = 1E17
QThresh = 5E19
SPC = 2
DensityColor = False
QLog = True
Param = 'Ne' # 'IonPol'

AA = SOLPSPLOT(Shot, Attempt, LINTHRESH=LINTHRESH, JXA=JXA, JXI=JXI, SEP=SEP)

X=AA.PolVec.loc[:,:,Attempt,'XXLoc'].values
Y=AA.RadCoords['YYLoc'].loc[:,:,Attempt].values
U=AA.PARAM['IonPol'].loc[:,:,Attempt].values
V=AA.PARAM['IonFlx'].loc[:,:,Attempt].values

if Quiver:
    figQ, axQ = plt.subplots()
       
    AA.Contour(Param,LOG10=2,PHYS=False,POLC='X',RADC='Y',AX=axQ)

    if DensityColor:
        M=np.log10(AA.PARAM['Ne'].loc[:,:,Attempt].values[::SPC,::SPC])
        CMAP='viridis'
    else:
        M=np.zeros(X[::SPC,::SPC].shape)
        CMAP='binary_r'
     
    if QLog:
        SLM=colors.SymLogNorm(QThresh,vmin=np.min([U,V]),vmax=np.max([U,V]))
        UU=SLM(U)
        VV=SLM(V)
        UU=UU-0.5
        VV=VV-0.5
    else:
        UU = U
        VV = V

    axQ.quiver(X[::SPC,::SPC],Y[::SPC,::SPC],
               UU[::SPC,::SPC],VV[::SPC,::SPC],
               M, pivot='mid', cmap=CMAP)

if StreamPlot:
    
    figS, axS = plt.subplots()
       
    AA.Contour(Param,LOG10=2,PHYS=False,POLC='X',RADC='Y',AX=axS)
    
    LB = AA.ContKW['CoreBound'][0]
    RB = AA.ContKW['CoreBound'][1]+1
    
    axS.streamplot(X[SEP:,:],Y[SEP:,:],U[SEP:,:],V[SEP:,:],color='black',linewidth=2)
    axS.streamplot(X[:,LB:RB],Y[:,LB:RB],U[:,LB:RB],V[:,LB:RB], color='black',linewidth=2)
    axS.streamplot(X[:,RB:],Y[:,RB:],U[:,RB:],V[:,RB:], color='black',linewidth=2)
    axS.streamplot(X[:,:LB],Y[:,:LB],U[:,:LB],V[:,:LB], color='black',linewidth=2)
    
    
    
    
    
    
    