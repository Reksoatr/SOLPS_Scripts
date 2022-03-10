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
from B2TransportParser import R2PsiN, PsiN2R

BASE,TOP = SET_WDIR('','')
Plot= True
RADC='psin' #'radial' #
EMISS='tomo' #'tree'#
LCFS_Te=[83,90,93]
ATTEMPTS308=['14Rf0.35','14Rf0.7','14Rf1.05']
ATTEMPTS108=['18Rf0.3','18Rf0.6','18Rf0.9']
ATTEMPTS305=['24Rf1.0','24Rf2.0','24Rf3.0']

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


solps308=SOLPSPLOT('1100308004', ATTEMPTS308, Markers=False, JXA=40, JXI=59, 
                   PlotScheme=['b--','b-','b--'])
solps108=SOLPSPLOT('1080416025', ATTEMPTS108, Markers=False, JXA=40, JXI=58, 
                   PlotScheme=['b--','b-','b--'], PsinOffset=0.02, RadOffset=0.005)
solps305=SOLPSPLOT('1100305023', ATTEMPTS305, Markers=False, PlotScheme=['b--','b-','b--'],
                   PsinOffset=-0.014, RadOffset=-0.007)

LYMID_coords=pkl.load(open('{}lya_coords_v3.pkl'.format(GFCMOD),'rb'))
RLYMID=np.sqrt(LYMID_coords['tang']['X']**2 + LYMID_coords['tang']['Y']**2)
RLYMID=np.flip(RLYMID)

LYMID308 = LYMID108 = LYMID305 = np.zeros((len(RLYMID),3))

for n,i in enumerate(ATTEMPTS308):
    LYMID308[:,n]=np.loadtxt('{}SOLPS_2d_prof/cmod/1100308004home/Attempt{}/Output/LyaBrightW{}'.format(TOP,i,i),
                    skiprows=2)[:,1]

for n,i in enumerate(ATTEMPTS108):    
    LYMID108[:,n]=np.loadtxt('{}SOLPS_2d_prof/cmod/1080416025home/Attempt{}/Output/LyaBrightW{}'.format(TOP,i,i),
                    skiprows=2)[:,1]

for n,i in enumerate(ATTEMPTS305):
    LYMID305[:,n]=np.loadtxt('{}SOLPS_2d_prof/cmod/1100305023home/Attempt{}/Output/LyaBrightW{}'.format(TOP,i,i),
                    skiprows=2)[:,1]

if Plot:

    if RADC=='radial':
        Rbright308=bright308[0]
        Rbright108=bright108[0]
        Rbright305=bright305[0]
        
        Remiss308=bright308[Remiss_idx]
        Remiss108=bright108[Remiss_idx]
        Remiss305=bright305[Remiss_idx]
        
        Rnn308=PsiN2R(GF308,nn308[0])
        Rnn108=PsiN2R(GF108,nn108[0])
        Rnn305=PsiN2R(GF305,nn305[0])
        
        Sep308=PsiN2R(GF308,1.0)
        Sep108=PsiN2R(GF108,1.0)
        Sep305=PsiN2R(GF308,1.0)
        
        RLYMID308=RLYMID
        RLYMID108=RLYMID
        RLYMID305=RLYMID
        
        Rlabel=r'Major Radius [m]'
        
        '''
        Rnn308=np.array([[j for j,k in zip(*GF308.get_fluxsurface(i)) if j > GF308.axis.r and k==0] 
                         for i in nn308[0]])[:,0]
        Rnn108=np.array([[j for j,k in zip(*GF108.get_fluxsurface(i)) if j > GF108.axis.r and k==0] 
                         for i in nn108[0]])[:,0]
        Rnn305=np.array([[j for j,k in zip(*GF305.get_fluxsurface(i)) if j > GF305.axis.r and k==0] 
                         for i in nn305[0]])[:,0]    
        
        Sep308=GF308.get_fluxsurface(1.0)[0].max()
        Sep108=GF108.get_fluxsurface(1.0)[0].max()
        Sep305=GF305.get_fluxsurface(1.0)[0].max()
        '''
        
    elif RADC=='psin':
        Rbright308=R2PsiN(GF308,bright308[0])
        Rbright108=R2PsiN(GF108,bright108[0])
        Rbright305=R2PsiN(GF305,bright305[0])
        
        Remiss308=R2PsiN(GF308,bright308[Remiss_idx])
        Remiss108=R2PsiN(GF108,bright108[Remiss_idx])
        Remiss305=R2PsiN(GF305,bright305[Remiss_idx])
        
        Rnn308=nn308[0]
        Rnn108=nn108[0]
        Rnn305=nn305[0]
        
        Sep308=1
        Sep108=1
        Sep305=1
        
        RLYMID308=R2PsiN(GF308,RLYMID)
        RLYMID108=R2PsiN(GF108,RLYMID)
        RLYMID305=R2PsiN(GF305,RLYMID)
        
        Rlabel=r'$\Psi_n$'
    
    #fig308,ax308=plt.subplots(5,1,sharex=True)
    #fig108,ax108=plt.subplots(5,1,sharex=True)
    fig305,ax305=plt.subplots(6,3,sharex=True,sharey='row')
    
    #L-mode plots
    solps308.RadProf('LyaEmissW',RADC=RADC,AX=ax305[1,0], Publish=['$\pm50\%$','Base',None])
    solps308.RadProf('NeuDen',RADC=RADC,AX=ax305[2,0],LOG10=1, Publish=['$\pm50\%$','Base',None])
    ax305[0,0].plot(RLYMID308,LYMID308[:,0],'b--',linewidth='2.0',label='+/-50%')
    ax305[0,0].plot(RLYMID308,LYMID308[:,1],'b-',linewidth='2.0',label='Base')
    ax305[0,0].plot(RLYMID308,LYMID308[:,2],'b--',linewidth='2.0',label=None)
    ax305[0,0].plot(Rbright308,bright308[1],'--*',color='red',label='LYMID')
    ax305[0,0].legend()
    ax305[1,0].plot(Remiss308,bright308[emiss_idx],'--*',color='red',label='LYMID')
    ax305[1,0].legend()
    gaussian_shading(ax305[2,0],Rnn308,np.log10(nn308[1]),nnlogerr308,c='red')
    ax305[2,0].legend()
    solps308.RadProf('Ne',RADC=RADC,AX=ax305[3,0],PlotScheme=['r--','b-','b--'],Publish=['$\pm50\%$','Base',None])
    ax305[3,0].legend(loc=2)
    solps308.RadProf('Te',RADC=RADC,AX=ax305[4,0],PlotScheme=['r--','b-','b--'],Publish=['$\pm50\%$','Base',None])
    ax305[4,0].axhline(LCFS_Te[0],linestyle=':',color='orange',label='2PM Te')
    ax305[4,0].legend(loc=1)
    solps308.RadProf('DN',RADC=RADC,AX=ax305[5,0],PlotScheme=['g--','g-','g--'],Publish=['$\pm 50\%\;D_n$','$D_n$',None])
    solps308.RadProf('KYE',RADC=RADC,AX=ax305[5,0],PlotScheme=['m--','m-','m--'],Publish=['$\pm 50\%\;\chi_e$','$\chi_e$',None])
    ax305[5,0].legend(loc=1)
    
    ax305[0,0].axvline(Sep308,linestyle=':')
    ax305[1,0].axvline(Sep308,linestyle=':')
    ax305[2,0].axvline(Sep308,linestyle=':')
    ax305[3,0].axvline(Sep308,linestyle=':')
    ax305[4,0].axvline(Sep308,linestyle=':')
    ax305[5,0].axvline(Sep308,linestyle=':')
    
    ax305[0,0].set_title(r'1100308004', pad=10.0)
    ax305[0,0].set_yticklabels([])
    ax305_bright=ax305[0,0].secondary_yaxis('left',functions=(lambda x: x/1000, lambda x: x*1000))
    ax305_bright.set_ylabel('Brightness\n ($kW/m^2$)')
    ax305[1,0].set_title('')
    ax305[1,0].set_yticklabels([])
    ax305_emiss=ax305[1,0].secondary_yaxis('left',functions=(lambda x: x/1000, lambda x: x*1000))
    ax305_emiss.set_ylabel('Emissivity\n ($kW/m^3$)')
    ax305[1,0].set_ylabel('')
    ax305[1,0].set_xlabel('')
    ax305[2,0].set_title('')
    ax305[2,0].set_ylabel('$Log_{10}$ of Neutral\n Density ($10^y\;m^{-3}$)')
    ax305[2,0].set_xlabel('')
    ax305[3,0].set_title('')
    ax305[3,0].set_ylabel('Electron\n Density $n_e\;(m^{-3})$')
    ax305[3,0].set_xlabel('')
    ax305[4,0].set_title('')
    ax305[4,0].set_ylabel('Electron\n Temperature $T_e\;(eV)$')
    ax305[4,0].set_xlabel('')
    ax305[5,0].set_title('')
    ax305[5,0].set_ylabel('Transport\n Coefficient $(m^2/s)$')
    ax305[5,0].set_xlabel(Rlabel)
    
    #I-mode plots
    solps108.RadProf('LyaEmissW',RADC=RADC,AX=ax305[1,1], Publish=['$\pm50\%$','Base',None])
    solps108.RadProf('NeuDen',RADC=RADC,AX=ax305[2,1],LOG10=1, Publish=['$\pm50\%$','Base',None])
    ax305[0,1].plot(RLYMID108,LYMID108[:,0],'b--',linewidth='2.0',label='+/-50%')
    ax305[0,1].plot(RLYMID108,LYMID108[:,1],'b-',linewidth='2.0',label='Base')
    ax305[0,1].plot(RLYMID108,LYMID108[:,2],'b--',linewidth='2.0',label=None)
    ax305[0,1].plot(Rbright108,bright108[1],'--*',color='red',label='LYMID')
    ax305[0,1].legend()
    ax305[1,1].plot(Remiss108,bright108[emiss_idx],'--*',color='red',label='LYMID')
    ax305[1,1].legend()
    gaussian_shading(ax305[2,1],Rnn108,np.log10(nn108[1]),nnlogerr108,c='red')
    ax305[2,1].legend()
    solps108.RadProf('Ne',RADC=RADC,AX=ax305[3,1],PlotScheme=['r--','b-','b--'],Publish=['$\pm50\%$','Base',None])
    ax305[3,1].legend(loc=2)
    solps108.RadProf('Te',RADC=RADC,AX=ax305[4,1],PlotScheme=['r--','b-','b--'], Publish=['$\pm50\%$','Base',None])
    ax305[4,1].axhline(LCFS_Te[1],linestyle=':',color='orange',label='2PM Te')
    ax305[4,1].legend(loc=1)
    solps108.RadProf('DN',RADC=RADC,AX=ax305[5,1],PlotScheme=['g--','g-','g--'],Publish=['$\pm 50\%\;D_n$','$D_n$',None])
    solps108.RadProf('KYE',RADC=RADC,AX=ax305[5,1],PlotScheme=['m--','m-','m--'],Publish=['$\pm 50\%\;\chi_e$','$\chi_e$',None])
    ax305[5,1].legend(loc=1)
    
    ax305[0,1].axvline(Sep108,linestyle=':')
    ax305[1,1].axvline(Sep108,linestyle=':')
    ax305[2,1].axvline(Sep108,linestyle=':')
    ax305[3,1].axvline(Sep108,linestyle=':')
    ax305[4,1].axvline(Sep108,linestyle=':')
    ax305[5,1].axvline(Sep108,linestyle=':')
    
    ax305[0,1].set_title(r'1080416025', pad=10.0)
    ax305[0,1].set_ylabel('')
    ax305[1,1].set_title('')
    ax305[1,1].set_ylabel('')
    ax305[1,1].set_xlabel('')
    ax305[2,1].set_title('')
    ax305[2,1].set_ylabel('')
    ax305[2,1].set_xlabel('')
    ax305[3,1].set_title('')
    ax305[3,1].set_ylabel('')
    ax305[3,1].set_xlabel('')
    ax305[4,1].set_title('')
    ax305[4,1].set_ylabel('')
    ax305[4,1].set_xlabel('')
    ax305[5,1].set_title('')
    ax305[5,1].set_ylabel('')
    ax305[5,1].set_xlabel(Rlabel)
    
    #H-mode plots
    solps305.RadProf('LyaEmissW',RADC=RADC,AX=ax305[1,2], Publish=['$\pm50\%$','Base',None])
    solps305.RadProf('NeuDen',RADC=RADC,AX=ax305[2,2],LOG10=1, Publish=['$\pm50\%$','Base',None])                   
    ax305[0,2].plot(RLYMID305,LYMID305[:,0],'b--',linewidth='2.0',label='+/-50%')
    ax305[0,2].plot(RLYMID305,LYMID305[:,1],'b-',linewidth='2.0',label='Base')
    ax305[0,2].plot(RLYMID305,LYMID305[:,2],'b--',linewidth='2.0',label=None)
    ax305[0,2].plot(Rbright305,bright305[1],'--*',color='red',label='LYMID')
    ax305[0,2].legend()
    ax305[1,2].plot(Remiss305,bright305[emiss_idx],'--*',color='red',label='LYMID')
    ax305[1,2].legend()
    gaussian_shading(ax305[2,2],Rnn305,np.log10(nn305[1]),nnlogerr305,c='red')
    ax305[2,2].legend()
    solps305.RadProf('Ne',RADC=RADC,AX=ax305[3,2],PlotScheme=['r--','b-','b--'], Publish=['$\pm50\%$','Base',None])
    ax305[3,2].legend(loc=2)
    solps305.RadProf('Te',RADC=RADC,AX=ax305[4,2],PlotScheme=['r--','b-','b--'], Publish=['$\pm50\%$','Base',None])
    ax305[4,2].axhline(LCFS_Te[2],linestyle=':',color='orange',label='2PM Te')
    ax305[4,2].legend(loc=1)
    solps305.RadProf('DN',RADC=RADC,AX=ax305[5,2],PlotScheme=['g--','g-','g--'],Publish=['$\pm 50\%\;D_n$','$D_n$',None])
    solps305.RadProf('KYE',RADC=RADC,AX=ax305[5,2],PlotScheme=['m--','m-','m--'],Publish=['$\pm 50\%\;\chi_e$','$\chi_e$',None])
    ax305[5,2].legend(loc=1)
    
    ax305[0,2].axvline(Sep305,linestyle=':')
    ax305[1,2].axvline(Sep305,linestyle=':')
    ax305[2,2].axvline(Sep305,linestyle=':')
    ax305[3,2].axvline(Sep305,linestyle=':')
    ax305[4,2].axvline(Sep305,linestyle=':')
    ax305[5,2].axvline(Sep305,linestyle=':')
    
    ax305[0,2].set_title(r'1100305023', pad=10.0)
    ax305[0,2].set_ylabel('')
    ax305[0,2].secondary_yaxis('right', functions=(lambda x: x/1000, lambda x: x*1000))
    ax305[1,2].set_title('')
    ax305[1,2].set_ylabel('')
    ax305[1,2].secondary_yaxis('right', functions=(lambda x: x/1000, lambda x: x*1000))
    ax305[1,2].set_xlabel('')
    ax305[2,2].set_title('')
    ax305[2,2].set_ylabel('')
    ax305[2,2].secondary_yaxis('right')
    ax305[2,2].set_xlabel('')
    ax305[3,2].set_title('')
    ax305[3,2].set_ylabel('')
    ax305[3,2].secondary_yaxis('right')
    ax305[3,2].set_xlabel('')
    ax305[4,2].set_title('')
    ax305[4,2].set_ylabel('')
    ax305[4,2].secondary_yaxis('right')
    ax305[4,2].set_xlabel('')
    ax305[5,2].set_title('')
    ax305[5,2].set_ylabel('')
    ax305[5,2].secondary_yaxis('right')
    ax305[5,2].set_xlabel(Rlabel)
    
    ax2_308=ax305[0,0].secondary_xaxis('top',functions=(lambda x: PsiN2R(GF108,x), lambda x: x))
    ax2_108=ax305[0,1].secondary_xaxis('top',functions=(lambda x: PsiN2R(GF108,x), lambda x: x))
    ax2_305=ax305[0,2].secondary_xaxis('top',functions=(lambda x: PsiN2R(GF305,x), lambda x: x))
    
    ax2_308.set_xlabel('Major Radius [m]')
    ax2_108.set_xlabel('Major Radius [m]')
    ax2_305.set_xlabel('Major Radius [m]')
    
    fig305.subplots_adjust(wspace=0, hspace=0)
    