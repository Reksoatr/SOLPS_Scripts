# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 20:36:33 2021

@author: 18313
"""
from aurora.solps import solps_case
from SOLPS_Scripts.TOOLS import SET_WDIR
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl

plt.rc('axes',titlesize=25,labelsize=25)
plt.rc('lines',linewidth=5,markersize=20,markeredgewidth=2,linestyle='solid')
plt.rc('legend',fontsize=25,title_fontsize=25)
plt.rc('xtick',labelsize=25)
plt.rc('ytick',labelsize=25)
plt.rc('xtick.major',size=10,width=3)
plt.rc('ytick.major',size=10,width=3)

BASE,TOP = SET_WDIR('','')

fig10,ax10 = plt.subplots(nrows=4, ncols=3, sharex=True, sharey='row')

shot_list = [1100308004,1080416025,1100305023]
labels = ['L-mode', 'I-mode', 'H-mode']

run_list = ['/Attempt14R/Output','/Attempt18R/Output','/Attempt24R/Output']

path_list = ['{}cmod/{}home'.format(BASE,s) for s in shot_list]

gfilepath_list = ['{}gfileProcessing/cmod_files/g{}.01000'.format(TOP,s) for s in shot_list]

label_list = ['{} ({})'.format(*k) for k in zip(labels,shot_list)]

LCFS_Te=[83,90,93]
ne_max=None
Te_max=None

CONTOUR=False
A_max=None
A_min=None

try: SO
except NameError: 
    print('Creating SO Dictionary')
    SO={}

for i,n in enumerate(shot_list):
    
    #Initiate solps_case object and load data
    
    if n not in SO.keys():
        SO[n] = solps_case(path_list[i], gfilepath_list[i], solps_run=run_list[i], form='full')
        
    rhop_fsa, nn_fsa, rhop_LFS, nn_LFS, rhop_HFS, nn_HFS = \
        SO[n].get_radial_prof(SO[n].quants['nn'], dz_mm=1, plot=False, method='cubic')
    
    nerhop_fsa, ne_fsa, nerhop_LFS, ne_LFS, nerhop_HFS, ne_HFS = \
        SO[n].get_radial_prof(SO[n].quants['ne'], dz_mm=1, plot=False, method='cubic')

    Terhop_fsa, Te_fsa, Terhop_LFS, Te_LFS, Terhop_HFS, Te_HFS = \
        SO[n].get_radial_prof(SO[n].quants['Te'], dz_mm=1, plot=False, method='cubic')
    
    dnarhop_fsa, dna_fsa, dnarhop_LFS, dna_LFS, dnarhop_HFS, dna_HFS = \
        SO[n].get_radial_prof(SO[n].b2fplasmf['dna0'][1], dz_mm=1, plot=False, method='cubic')
    
    hcerhop_fsa, hce_fsa, hcerhop_LFS, hce_LFS, hcerhop_HFS, hce_HFS = \
        SO[n].get_radial_prof(SO[n].b2fplasmf['hce0'][0]/1e20, dz_mm=1, plot=False, method='cubic')
    
    #Plot 1D LFS Profiles
    
    ax10[0,i].plot(dnarhop_LFS**2, dna_LFS, linewidth=2,color='green',label=r'$D$')
    ax10[0,i].plot(hcerhop_LFS**2, hce_LFS, linewidth=2,color='red',label=r'$\chi_e=\chi_i$')
    ax10[1,i].plot(nerhop_LFS**2, ne_LFS, label=label_list[i],linewidth=2,color='black')
    ax10[2,i].plot(Terhop_LFS**2, Te_LFS, label=label_list[i],linewidth=2,color='black')
    ax10[3,i].semilogy(rhop_LFS**2, nn_LFS, label=label_list[i],linewidth=2,color='black')
    
    #Load and Plot Experiment TS Data
    
    with open('{}gfileProcessing/cmod_files/{}_CORE.pkl'.format(TOP,shot_list[i]),'rb') as f:
        u = pkl._Unpickler(f)
        u.encoding = 'latin1'
        CTS = u.load()
    
        #Align Data to LCFS Temperature If Specified
    
        if LCFS_Te!=0:
            index=(np.abs(CTS[4]*1000-LCFS_Te[i])).argmin()
            PsinOffset=1-CTS[3][index]
            print('Psin offset for shot {} is {}'.format(shot_list[i],PsinOffset))
    
        CoreNe = ax10[1,i].errorbar(CTS[0]+PsinOffset,CTS[1]*1e20,yerr=CTS[2]*1e20,
                                    linestyle='',capsize=5,marker='.',color='blue')
        CoreTe = ax10[2,i].errorbar(CTS[3]+PsinOffset,CTS[4]*1000,yerr=CTS[5]*1000,
                                    linestyle='',capsize=5,marker='.',color='blue')

    #Plot Separatrix

    ax10[0,i].axvline(1.0,linestyle='dotted')
    ax10[1,i].axvline(1.0,linestyle='dotted')
    ax10[2,i].axvline(1.0,linestyle='dotted')
    ax10[3,i].axvline(1.0,linestyle='dotted')
    
    #Find Max Value of SOLPS Ne and Te Profiles, Set Upper Y-Limit to 1.05 Max Value

    if ne_max: 
        ne_max=np.maximum(ne_max,1.05*np.nanmax(ne_LFS))
    else:
        ne_max=1.05*np.nanmax(ne_LFS)
        
    if Te_max:
        Te_max=np.maximum(Te_max,1.05*np.nanmax(Te_LFS))
    else:
        Te_max=1.05*np.nanmax(Te_LFS)

    ax10[1,i].set_ylim(ymin=0,ymax=ne_max)
    ax10[2,i].set_ylim(ymin=0,ymax=Te_max)
    
    ax10[-1,i].set_xlim(xmin=np.min(rhop_LFS**2),xmax=np.max(rhop_LFS**2))
    ax10[-1,i].set_xlabel(r'$\psi_n$')
    
    #Determine Max and Min of All Contours
    if CONTOUR:
        #PDENA=SO[n].fort46['pdena']*1e6
        PDENA=SO[n].quants['ne']
        if np.any(PDENA==0): PDENA[PDENA==0.0] = np.unique(PDENA)[1]
        
        if A_max:
            A_max=np.maximum(PDENA.max(),A_max)
        else:
            A_max=PDENA.max()
            
        if A_min:
            A_min=np.minimum(PDENA.min(),A_min)
        else:
            A_min=PDENA.min()
            
        print('Contour minimum is now {}'.format(A_min))
        print('Contour maximum is now {}'.format(A_max))    
    
#Plot 2D Contours with Shared Limits
    
if CONTOUR:
    
    figA,axA = plt.subplots(nrows=1,ncols=3,sharex=True,sharey=True)

    figA.suptitle(r'Electron Density ($m^{-3}$)',fontsize=25)
    
    for i,n in enumerate(shot_list):
        #PDENA=SO[n].fort46['pdena']*1e6
        SO[n].plot2d_b2(SO[n].quants['ne'],ax=axA[i],lb=A_min,ub=A_max)
        #SO[n].plot2d_eirene(PDENA,ax=axA[i],lb=A_min,ub=A_max)
        SO[n].geqdsk.plot(only2D=True,ax=axA[i],color='grey')
        axA[i].set_title(label_list[i])
        
        axA[i].plot()
        
    # Remove Redundant Colorbars

    for ii in range(i):
        figA.get_children()[i+2].remove()
        axA[ii+1].set_ylabel('')

ax10[0,1].legend(loc=9,ncol=2)

ax10[0,0].set_ylabel('Transport \n coeff. $(m^2/s)$')
ax10[1,0].set_ylabel(r'$n_e\;(m^{-3})$')
ax10[2,0].set_ylabel(r'$T_e\;(eV)$')
ax10[3,0].set_ylabel(r'$n_n\;(m^{-3})$')

ax10[0,0].set_title(label_list[0])
ax10[0,1].set_title(label_list[1])
ax10[0,2].set_title(label_list[2])

fig10.subplots_adjust(wspace=0,hspace=0)
