# -*- coding: utf-8 -*-
"""
Created on Sun May 12 10:59:21 2019

@author: Richard Reksoatmodjo
"""

import numpy as np
import equilibrium as eq
from scipy.io import loadmat
import matplotlib.pyplot as plt

def PsiNtoR(GFILE,Psin,Z=0,Side=1):
    GF = eq.equilibrium(gfile=GFILE)
    if type(Psin)==int or type(Psin)==float:
        Rset = GF.get_fluxsurface(Psin)
        R = Rset[0][Rset[1]==Z][Side]
        return R
    else:
        R = np.zeros(len(Psin))
        for nn in range(len(Psin)):
            Rset = GF.get_fluxsurface(Psin[nn])
            R[nn]=Rset[0][Rset[1]==Z][Side]
        return R
    
def RhotoPsiN(rho):
    PSIN = rho**2
    return PSIN

if __name__ == '__main__':

    PRINT = 0
    PLOT = 0
    
    Z_mid = 0
    R_OMcore = 2.176
    R_IMcore = 1.2354 
    R_OMSOL = 2.305
    
    ExpData = loadmat('gfileProcessing/d3d_files/175060_data_SOLPS.mat')
    
    TC = loadmat('gfileProcessing/d3d_files/flow_transport.mat')
    
    GFILE = 'gfileProcessing/d3d_files/g175060.02512'
    
    GF = eq.equilibrium(gfile=GFILE)
    
    psin_core = GF.psiN(R_OMcore,Z_mid)
    psin_SOL = GF.psiN(R_OMSOL,Z_mid)
    
    print(psin_core)
    print(psin_SOL)
    
    i = 0
    j = 0
    k = 0
    TC_i = 0
    TC_step = 2
    
    while ExpData['psin_ne'][i] < psin_core:
        i += 1
    
    print('Core ne = ',ExpData['Ne'][i])
        
    while ExpData['psin_te'][j] < psin_core:
        j += 1
        
    print('Core Te = ',ExpData['Te'][j])
     
    while ExpData['psin_ti'][k] < psin_core:
        k += 1
    
    print('Core Ti = ',ExpData['Ti'][k])
    
    TC_psin = RhotoPsiN(TC['rho'])
     
    while TC_psin[TC_i] < psin_core:
        TC_i += 1
    TC_i -= 1
    
    TC_pflux = TC['pflux'][TC_i:]
    
    print('Core Particle Flux = ',TC_pflux[0])
    
    TC_psin = TC_psin[TC_i:]
    
    R_sep = PsiNtoR(GFILE,TC_psin[-1])
    
    TC_RR = PsiNtoR(GFILE,TC_psin)
    
    TC_R = TC_RR-R_sep
    
    if PLOT==1:    
        figl = plt.figure(figsize=(14,10))
        plt.plot(ExpData['psin_ne'][i:],ExpData['Ne'][i:])
        
        fig2 = plt.figure(figsize=(14,10))
        plt.plot(ExpData['psin_te'][j:],ExpData['Te'][j:]) 
        
        fig3 = plt.figure(figsize=(14,10))
        plt.plot(ExpData['psin_ti'][k:],ExpData['Ti'][k:]) 
        
        fig4 = plt.figure(figsize=(14,10))
        plt.plot(TC_R[0::TC_step],TC['D'][TC_i::TC_step]) 
        
        fig5 = plt.figure(figsize=(14,10))
        plt.plot(TC_R[0::TC_step],TC['chie'][TC_i::TC_step]) 
        
        fig6 = plt.figure(figsize=(14,10))
        plt.plot(TC_R[0::TC_step],TC['chii'][TC_i::TC_step])
        
        fig7 = plt.figure(figsize=(14,10))
        plt.plot(TC_R[0::TC_step],TC['pflux'][TC_i::TC_step])
    
    if PRINT==1:
        n = len(TC_R[0::TC_step])
        
        print(' ndata(1, 1, 1)= {0},'.format(n))
        for m in range(n):
            print(' tdata(1, {0}, 1, 1)= {1}, tdata(2, {0}, 1, 1)= {2},'.format(m+1,'%.3f' % float(TC_R[0::TC_step][m]),'%.3f' %TC['D'][TC_i::TC_step][m]))
        
        print(' ndata(1, 3, 1)= {0},'.format(n))
        for m in range(n):
            print(' tdata(1, {0}, 3, 1)= {1}, tdata(2, {0}, 3, 1)= {2},'.format(m+1,'%.3f' %TC_R[0::TC_step][m],'%.3f' %TC['chie'][TC_i::TC_step][m]))
        
        print(' ndata(1, 4, 1)= {0},'.format(n))
        for m in range(n):
            print(' tdata(1, {0}, 4, 1)= {1}, tdata(2, {0}, 4, 1)= {2},'.format(m+1,'%.3f' %TC_R[0::TC_step][m],'%.3f' %TC['chii'][TC_i::TC_step][m]))
     
    
    
    
    
    
    
    
    
    
    
