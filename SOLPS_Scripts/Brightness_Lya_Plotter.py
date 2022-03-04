# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 00:17:59 2021

@author: 18313

Script to plot Lyman-alpha brightness and emissivity profiles for validation
study discharges, and compare to SOLPS profiles 
"""

import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import equilibrium as eq
from TOOLS import SET_WDIR, gaussian_shading
from VesselPlotterNew import SOLPSPLOT

BASE,TOP = SET_WDIR('','')
Plot= True
RADC='psin' #'radial' #
EMISS='tomo' #'tree'#

if EMISS=='tree':
    Remiss_idx=2
    emiss_idx=3
elif EMISS=='tomo':
    Remiss_idx=4
    emiss_idx=5
    
# Load Experimental Data

GFCMOD='{}gfileProcessing/cmod_files/'.format(TOP)

GF308=eq.equilibrium(gfile='{}g1100308004.01000'.format(GFCMOD))
GF108=eq.equilibrium(gfile='{}g1080416025.01000'.format(GFCMOD))
GF305=eq.equilibrium(gfile='{}g1100305023.01000'.format(GFCMOD))

bright308=pkl.load(open('lya_brightness_1100308004.pkl','rb'))
bright108=pkl.load(open('lya_brightness_1080416025.pkl','rb'))
bright305=pkl.load(open('lya_brightness_1100305023.pkl','rb'))

nn308=pkl.load(open('{}Brightness_Profiles/lyman_data_1100308004.pkl'.format(GFCMOD),'rb'))
nn108=pkl.load(open('{}Brightness_Profiles/lyman_data_1080416025.pkl'.format(GFCMOD),'rb'))
nn305=pkl.load(open('{}Brightness_Profiles/lyman_data_1100305023.pkl'.format(GFCMOD),'rb'))

nn308[1]=1e6*nn308[1]
nn308[2]=1e6*nn308[2]
nnlogerr308=0.434*nn308[2]/nn308[1]

nn108[1]=1e6*nn108[1]
nn108[2]=1e6*nn108[2]
nnlogerr108=0.434*nn108[2]/nn108[1]

nn305[1]=1e6*nn305[1]
nn305[2]=1e6*nn305[2]
nnlogerr305=0.434*nn305[2]/nn305[1]


solps308=SOLPSPLOT('1100308004', ['14Rf0.7'], Markers=False, JXA=40, JXI=59, PsinOffset=0.02)
solps108=SOLPSPLOT('1080416025', ['18Rf0.6'], Markers=False, JXA=40, JXI=58, PsinOffset=0.03)
solps305=SOLPSPLOT('1100305023', ['24Rf1.5'], Markers=False, PsinOffset=-0.014)

LYMID_coords=pkl.load(open('{}lya_coords_v3.pkl'.format(GFCMOD),'rb'))
RLYMID=np.sqrt(LYMID_coords['tang']['X']**2 + LYMID_coords['tang']['Y']**2)
RLYMID=np.flip(RLYMID)

LYMID308=np.loadtxt('{}SOLPS_2d_prof/cmod/1100308004home/Attempt14Rf0.7/Output/LyaBrightW14Rf0.7'.format(TOP),
                    skiprows=2).T
LYMID108=np.loadtxt('{}SOLPS_2d_prof/cmod/1080416025home/Attempt18Rf0.6/Output/LyaBrightW18Rf0.6'.format(TOP),
                    skiprows=2).T
LYMID305=np.loadtxt('{}SOLPS_2d_prof/cmod/1100305023home/Attempt24Rf1.5/Output/LyaBrightW24Rf1.5'.format(TOP),
                    skiprows=2).T

if Plot:

    if RADC=='radial':
        Rbright308=bright308[0]
        Rbright108=bright108[0]
        Rbright305=bright305[0]
        
        Remiss308=bright308[Remiss_idx]
        Remiss108=bright108[Remiss_idx]
        Remiss305=bright305[Remiss_idx]
        
        Rnn308=np.array([[j for j,k in zip(*GF308.get_fluxsurface(i)) if j > GF308.axis.r and k==0] 
                         for i in nn308[0]])[:,0]
        Rnn108=np.array([[j for j,k in zip(*GF108.get_fluxsurface(i)) if j > GF108.axis.r and k==0] 
                         for i in nn108[0]])[:,0]
        Rnn305=np.array([[j for j,k in zip(*GF305.get_fluxsurface(i)) if j > GF305.axis.r and k==0] 
                         for i in nn305[0]])[:,0]    
    
        Sep308=GF308.get_fluxsurface(1.0)[0].max()
        Sep108=GF108.get_fluxsurface(1.0)[0].max()
        Sep305=GF305.get_fluxsurface(1.0)[0].max()
        
        RLYMID308=RLYMID
        RLYMID108=RLYMID
        RLYMID305=RLYMID
        
        Rlabel=r'Major Radius [m]'
        
    elif RADC=='psin':
        Rbright308=GF308.psiN(bright308[0],0)[0]
        Rbright108=GF108.psiN(bright108[0],0)[0]
        Rbright305=GF305.psiN(bright305[0],0)[0]
        
        Remiss308=GF308.psiN(bright308[Remiss_idx],0)[0]
        Remiss108=GF108.psiN(bright108[Remiss_idx],0)[0]
        Remiss305=GF305.psiN(bright305[Remiss_idx],0)[0]
        
        Rnn308=nn308[0]
        Rnn108=nn108[0]
        Rnn305=nn305[0]
        
        Sep308=1
        Sep108=1
        Sep305=1
        
        RLYMID308=GF308.psiN(RLYMID,0)[0]
        RLYMID108=GF108.psiN(RLYMID,0)[0]
        RLYMID305=GF305.psiN(RLYMID,0)[0]
        
        Rlabel=r'$\Psi_n$'
    
    fig308,ax308=plt.subplots(5,1,sharex=True)
    fig108,ax108=plt.subplots(5,1,sharex=True)
    fig305,ax305=plt.subplots(5,1,sharex=True)
    
    #L-mode plots
    solps308.RadProf('LyaEmissW',RADC=RADC,AX=ax308[1])
    solps308.RadProf('NeuDen',RADC=RADC,AX=ax308[2],LOG10=1)
    ax308[0].plot(RLYMID308,LYMID308[1],color='#1f77b4ff',linestyle='-',linewidth='3.0')
    ax308[0].plot(Rbright308,bright308[1],'--*',color='#ff7f0eff')
    ax308[0].legend(['SOLPS','LYMID'])
    ax308[1].plot(Remiss308,bright308[emiss_idx],'--*')
    ax308[1].legend(['SOLPS','LYMID'])
    gaussian_shading(ax308[2],Rnn308,np.log10(nn308[1]),nnlogerr308,c='#ff7f0eff')
    ax308[2].legend(['SOLPS','LYMID'])
    solps308.RadProf('Ne',RADC=RADC,AX=ax308[3])
    ax308[3].legend(['SOLPS','LYMID'])
    solps308.RadProf('Te',RADC=RADC,AX=ax308[4])
    ax308[4].legend(['SOLPS','LYMID'])
    
    ax308[0].axvline(Sep308,linestyle=':')
    ax308[1].axvline(Sep308,linestyle=':')
    ax308[2].axvline(Sep308,linestyle=':')
    ax308[3].axvline(Sep308,linestyle=':')
    ax308[4].axvline(Sep308,linestyle=':')
    ax308[0].set_title(r'1100308004')
    ax308[0].set_ylabel(r'Brightness ($W/m^2$)')
    ax308[1].set_title('')
    ax308[1].set_ylabel(r'Emissivity ({}) ($W/m^3$)'.format(EMISS))
    ax308[1].set_xlabel('')
    ax308[2].set_title('')
    ax308[2].set_ylabel(r'$Log_{10}$ of Neutral Density ($10^y\;m^{-3}$)')
    ax308[2].set_xlabel('')
    ax308[3].set_title('')
    ax308[3].set_xlabel('')
    ax308[4].set_title('')
    ax308[4].set_xlabel(Rlabel)
    
    #I-mode plots
    solps108.RadProf('LyaEmissW',RADC=RADC,AX=ax108[1])
    solps108.RadProf('NeuDen',RADC=RADC,AX=ax108[2],LOG10=1)
    ax108[0].plot(RLYMID108,LYMID108[1],color='#1f77b4ff',linestyle='-',linewidth='3.0')
    ax108[0].plot(Rbright108,bright108[1],'--*',color='#ff7f0eff')
    ax108[0].legend(['SOLPS','LYMID'])
    ax108[1].plot(Remiss108,bright108[emiss_idx],'--*')
    ax108[1].legend(['SOLPS','LYMID'])
    gaussian_shading(ax108[2],Rnn108,np.log10(nn108[1]),nnlogerr108,c='#ff7f0eff')
    ax108[2].legend(['SOLPS','LYMID'])
    solps108.RadProf('Ne',RADC=RADC,AX=ax108[3])
    ax108[3].legend(['SOLPS','LYMID'])
    solps108.RadProf('Te',RADC=RADC,AX=ax108[4])
    ax108[4].legend(['SOLPS','LYMID'])
    
    ax108[0].axvline(Sep108,linestyle=':')
    ax108[1].axvline(Sep108,linestyle=':')
    ax108[2].axvline(Sep108,linestyle=':')
    ax108[3].axvline(Sep108,linestyle=':')
    ax108[4].axvline(Sep108,linestyle=':')
    ax108[0].set_title(r'1080416025')
    ax108[0].set_ylabel(r'Brightness ($W/m^2$)')
    ax108[1].set_title('')
    ax108[1].set_ylabel(r'Emissivity ({}) ($W/m^3$)'.format(EMISS))
    ax108[1].set_xlabel('')
    ax108[2].set_title('')
    ax108[2].set_ylabel(r'$Log_{10}$ of Neutral Density ($10^y\;m^{-3}$)')
    ax108[2].set_xlabel('')
    ax108[3].set_title('')
    ax108[3].set_xlabel('')
    ax108[4].set_title('')
    ax108[4].set_xlabel(Rlabel)
    
    #H-mode plots
    solps305.RadProf('LyaEmissW',RADC=RADC,AX=ax305[1])
    solps305.RadProf('NeuDen',RADC=RADC,AX=ax305[2],LOG10=1)                   
    ax305[0].plot(RLYMID305,LYMID305[1],color='#1f77b4ff',linestyle='-',linewidth='3.0')
    ax305[0].plot(Rbright305,bright305[1],'--*',color='#ff7f0eff')
    ax305[0].legend(['SOLPS','LYMID'])
    ax305[1].plot(Remiss305,bright305[emiss_idx],'--*')
    ax305[1].legend(['SOLPS','LYMID'])
    gaussian_shading(ax305[2],Rnn305,np.log10(nn305[1]),nnlogerr305,c='#ff7f0eff')
    ax305[2].legend(['SOLPS','LYMID'])
    solps305.RadProf('Ne',RADC=RADC,AX=ax305[3])
    ax305[3].legend(['SOLPS','LYMID'])
    solps305.RadProf('Te',RADC=RADC,AX=ax305[4])
    ax305[4].legend(['SOLPS','LYMID'])
    
    ax305[0].axvline(Sep305,linestyle=':')
    ax305[1].axvline(Sep305,linestyle=':')
    ax305[2].axvline(Sep305,linestyle=':')
    ax305[3].axvline(Sep305,linestyle=':')
    ax305[4].axvline(Sep305,linestyle=':')
    ax305[0].set_title(r'1100305023')
    ax305[0].set_ylabel(r'Brightness ($W/m^2$)')
    ax305[1].set_title('')
    ax305[1].set_ylabel(r'Emissivity ({}) ($W/m^3$)'.format(EMISS))
    ax305[1].set_xlabel('')
    ax305[2].set_title('')
    ax305[2].set_ylabel(r'$Log_{10}$ of Neutral Density ($10^y\;m^{-3}$)')
    ax305[2].set_xlabel(Rlabel)
    ax305[3].set_title('')
    ax305[3].set_xlabel('')
    ax305[4].set_title('')
    ax305[4].set_xlabel(Rlabel)