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
from scipy.interpolate import interp1d
from TOOLS import SET_WDIR, gaussian_shading
from SOLPS_Plotter import SOLPSPLOT
from B2TransportParser import R2PsiN, PsiN2R

plt.rc('font',size=30)
plt.rc('lines',linewidth=5,markersize=15)

Plot= True
Error_Analysis=True
RADC='psin' #'radial' #
EMISS='tomo' #'tree'#

SHOTS = ['1100308004','1080416025','1100305023'] #,
ATTEMPTS = [['14Rf0.35','14Rf0.7','14Rf1.05'],['18Rf0.3','18Rf0.6','18Rf0.9'],['24Rf1.0','24Rf2.0','24Rf3.0']] # ,
AN=len(SHOTS)

LCFS_Te=[83,90,93]
JXA=[40,40,55]
JXI=[59,58,37]
PsinOffset=[0,0.02,-0.014]
RadOffset=[0,0.005,-0.007]

if EMISS=='tree':
    Remiss_idx=2
    emiss_idx=3
elif EMISS=='tomo':
    Remiss_idx=4
    emiss_idx=5
    
# Load Experimental Data

BASE,TOP = SET_WDIR('','')

GFCMOD='{}gfileProcessing/cmod_files/'.format(TOP)

GFiles=[]

for A in SHOTS:
    GFiles.append(eq.equilibrium(gfile='{}g{}.01000'.format(GFCMOD,A)))

bright = [pkl.load(open('lya_brightness_{}.pkl'.format(A),'rb')) for A in SHOTS]

nn = [pkl.load(open('{}Brightness_Profiles/lyman_data_{}.pkl'.format(GFCMOD,A),'rb')) 
      for A in SHOTS]

nnlogerr=[]

for N in nn:
    N[1]=1e6*N[1]
    N[2]=1e6*N[2]
    nnlogerr.append(0.434*N[2]/N[1])

solps=[SOLPSPLOT(SHOTS[i], ATTEMPTS[i], Markers=False, JXA=JXA[i], JXI=JXI[i], 
                 PlotScheme=['b--','b-','b--'], PsinOffset=PsinOffset[i], RadOffset=RadOffset[i]) 
                 for i in range(AN)]

LYMID_coords=pkl.load(open('{}lya_coords_v3.pkl'.format(GFCMOD),'rb'))
RLYMID0=np.sqrt(LYMID_coords['tang']['X']**2 + LYMID_coords['tang']['Y']**2)
RLYMID0=np.flip(RLYMID0)

LYMID = np.zeros((AN,len(RLYMID0),3))

for m,a in enumerate(SHOTS):
    for n,i in enumerate(ATTEMPTS[m]):
        LYMID[m,:,n]=np.loadtxt('{}SOLPS_2d_prof/cmod/{}home/Attempt{}/Output/LyaBrightW{}'.format(TOP,a,i,i),
                        skiprows=2)[:,1]    

if Plot:
    
    Rbright=[]
    Remiss=[]    
    Rnn=[]
    Sep=[]
    RLYMID=[]

    if RADC=='radial':
        for m in range(AN):
            
            Rbright.append(bright[m][0])
            Remiss.append(bright[m][Remiss_idx])
            Rnn.append(PsiN2R(GFiles[m],nn[m][0]))
            Sep.append(PsiN2R(GFiles[m],nn[m][0]))
            RLYMID.append(RLYMID0)

            '''
            Rnn[m]=np.array([[j for j,k in zip(*GFiles[m].get_fluxsurface(i)) if j > GFiles[m].axis.r and k==0] 
                             for i in nn[m][0]])[:,0]  
            Sep[m]=GFiles[m].get_fluxsurface(1.0)[0].max()
            '''
            
        Rlabel=r'Major Radius [m]'
        
    elif RADC=='psin':
        for m in range(AN):
            
            Rbright.append(R2PsiN(GFiles[m],bright[m][0]))
            Remiss.append(R2PsiN(GFiles[m],bright[m][Remiss_idx]))
            Rnn.append(nn[m][0])
            Sep.append(1)
            RLYMID.append(R2PsiN(GFiles[m],RLYMID0))
        
        Rlabel=r'$\Psi_n$'
    
    figA,axA=plt.subplots(3,AN,sharex=True,sharey='row')
    figB,axB=plt.subplots(3,AN,sharex=True,sharey='row')
    
    if axA.shape==(6,):
        axA=np.array([axA]).T #Make sure axA always has 2 dimensions for indexing sake
    
    for p in range(AN):
        
        solps[p].RadProf('DN',RADC=RADC,AX=axA[0,p],PlotScheme=['g--','g-','g--'],Publish=['$\pm 50\%\;D_n$','$D_n$',None])
        solps[p].RadProf('KYE',RADC=RADC,AX=axA[0,p],PlotScheme=['m--','m-','m--'],Publish=['$\pm 50\%\;\chi_e$','$\chi_e$',None])
        axA[0,p].legend()
        
        solps[p].RadProf('Te',RADC=RADC,AX=axA[1,p],PlotScheme=['r--','b-','b--'],Publish=['$\pm50\%$','Base',None])
        axA[1,p].axhline(LCFS_Te[0],linestyle=':',color='orange',label='2PM Te')
        axA[1,p].legend()
        
        solps[p].RadProf('Ne',RADC=RADC,AX=axA[2,p],PlotScheme=['r--','b-','b--'],Publish=['$\pm50\%$','Base',None])
        axA[2,p].legend()
        
        axA[0,p].axvline(Sep[p],linestyle=':')
        axA[1,p].axvline(Sep[p],linestyle=':')
        axA[2,p].axvline(Sep[p],linestyle=':')
        
        axA[0,p].set_title(SHOTS[p], pad=10.0)
        axA[0,p].set_xlabel('')
        axA[0,p].set_ylabel('')
        axA[1,p].set_title('')
        axA[1,p].set_ylabel('')
        axA[1,p].set_xlabel('')
        axA[2,p].set_title('')
        axA[2,p].set_ylabel('')
        axA[2,p].set_xlabel(Rlabel)
        
        solps[p].RadProf('NeuDen',RADC=RADC,AX=axB[0,p],LOG10=1, Publish=['$\pm50\%$','Base',None])
        gaussian_shading(axB[0,p],Rnn[p],np.log10(nn[p][1]),nnlogerr[p],c='red')
        axB[0,p].legend(['$\pm50\%$','Base','_nolabel','LYMID'],loc=2)
        
        solps[p].RadProf('LyaEmissW',RADC=RADC,AX=axB[1,p], Publish=['$\pm50\%$','Base',None])
        axB[1,p].plot(Remiss[p],bright[p][emiss_idx],'--*',color='red',label='LYMID')
        axB[1,p].legend()
        
        axB[2,p].plot(RLYMID[p],LYMID[p][:,0],'b--',linewidth='2.0',label='+/-50%')
        axB[2,p].plot(RLYMID[p],LYMID[p][:,1],'b-',linewidth='2.0',label='Base')
        axB[2,p].plot(RLYMID[p],LYMID[p][:,2],'b--',linewidth='2.0',label=None)
        axB[2,p].plot(Rbright[p],bright[p][1],'--*',color='red',label='LYMID')
        axB[2,p].legend()
        
        axB[0,p].axvline(Sep[p],linestyle=':')
        axB[1,p].axvline(Sep[p],linestyle=':')
        axB[2,p].axvline(Sep[p],linestyle=':')
        
        axB[0,p].set_title(SHOTS[p], pad=10.0)
        axB[0,p].set_ylabel('')
        axB[0,p].set_xlabel('')
        axB[1,p].set_title('')
        axB[1,p].set_ylabel('')
        axB[1,p].set_xlabel('')
        axB[2,p].set_title('')
        axB[2,p].set_ylabel('')
        axB[2,p].set_xlabel(Rlabel)

        GF=GFiles[p]
        ax2A=axA[0,p].secondary_xaxis('top',functions=(lambda x: PsiN2R(GF,x), lambda x: x))
        ax2A.set_xlabel('Major Radius [m]')
        
        ax2B=axB[0,p].secondary_xaxis('top',functions=(lambda x: PsiN2R(GF,x), lambda x: x))
        ax2B.set_xlabel('Major Radius [m]')
        
        RR=axA[0,p].get_lines()
        axA[0,p].fill_between(RR[0].get_xdata(),RR[0].get_ydata(),RR[2].get_ydata(),alpha=0.5,color='green')
        axA[0,p].fill_between(RR[3].get_xdata(),RR[3].get_ydata(),RR[5].get_ydata(),alpha=0.5,color='magenta')
        
        for i in range(1,3):
            RR=axA[i,p].get_lines()
            axA[i,p].fill_between(RR[0].get_xdata(),RR[0].get_ydata(),RR[2].get_ydata(),alpha=0.5,color='blue')

        for i in range(0,3):
            RR=axB[i,p].get_lines()
            axB[i,p].fill_between(RR[0].get_xdata(),RR[0].get_ydata(),RR[2].get_ydata(),alpha=0.5,color='blue')
           
    axA[0,0].set_ylabel('Transport\n Coefficient $(m^2/s)$') 
    axA[1,0].set_ylabel('Electron\n Temperature $T_e\;(eV)$')
    axA[2,0].set_ylabel('Electron\n Density $n_e\;(m^{-3})$')
    
    axA[0,p].secondary_yaxis('right')
    axA[1,p].secondary_yaxis('right')
    axA[2,p].secondary_yaxis('right')

    axB[0,0].set_ylabel('$Log_{10}$ of Neutral\n Density ($10^y\;m^{-3}$)')
    axB[1,0].set_yticklabels([])
    axB_emiss=axB[1,0].secondary_yaxis('left',functions=(lambda x: x/1000, lambda x: x*1000))
    axB_emiss.set_ylabel('Emissivity\n ($kW/m^3$)')
    axB[2,0].set_yticklabels([])
    axB_bright=axB[2,0].secondary_yaxis('left',functions=(lambda x: x/1000, lambda x: x*1000))
    axB_bright.set_ylabel('Brightness\n ($kW/m^2$)')
       
    axB[0,p].secondary_yaxis('right')
    axB[1,p].secondary_yaxis('right', functions=(lambda x: x/1000, lambda x: x*1000))
    axB[2,p].secondary_yaxis('right', functions=(lambda x: x/1000, lambda x: x*1000))
        
    figA.subplots_adjust(wspace=0, hspace=0)
    figB.subplots_adjust(wspace=0, hspace=0)
    
    if Error_Analysis:
        
        Labels=['Lmode','Imode','Hmode']
        Styles=['r-','g-','b-']
        Plus_Style=['r^-','g^-','b^-']
        Minus_Style=['rv-','gv-','bv-']
        SOLPS_ne=[]
        SOLPS_te=[]
        SOLPS_nn=[]
        SOLPS_emiss=[]
        SOLPS_R=[]
        Interp_Exp_nn=[]
        Interp_Exp_emiss=[]
        Interp_Exp_bright=[]
        nn_res=[]
        nn_res_abs=[]
        emiss_res=[]
        emiss_res_abs=[]
        bright_res=[]
        bright_res_abs=[]
        
        ne_plus=[]
        ne_minus=[]
        te_plus=[]
        te_minus=[]
        nn_plus=[]
        nn_minus=[]
        emiss_plus=[]
        emiss_minus=[]
        bright_plus=[]
        bright_minus=[]
        
        firg_res,ax_res=plt.subplots(nrows=3,sharex=True)
        ax_nn=ax_res[0]
        ax_emiss=ax_res[1]
        ax_bright=ax_res[2]
        
        fig_err,ax_err=plt.subplots(nrows=5,sharex=True)
        ax_ne_err=ax_err[0]
        ax_te_err=ax_err[1]
        ax_nn_err=ax_err[2]
        ax_emiss_err=ax_err[3]
        ax_bright_err=ax_err[4]
        
        ax_ne_err.axhline(0,color='black',linestyle='-')
        ax_ne_err.axhline(0.5,color='black',linestyle=':')
        ax_ne_err.axhline(-0.5,color='black',linestyle=':')
        ax_te_err.axhline(0,color='black',linestyle='-')
        ax_te_err.axhline(0.5,color='black',linestyle=':')
        ax_te_err.axhline(-0.5,color='black',linestyle=':')
        ax_nn_err.axhline(0,color='black',linestyle='-')
        ax_nn_err.axhline(0.5,color='black',linestyle=':')
        ax_nn_err.axhline(-0.5,color='black',linestyle=':')
        ax_emiss_err.axhline(0,color='black',linestyle='-')
        ax_emiss_err.axhline(0.5,color='black',linestyle=':')
        ax_emiss_err.axhline(-0.5,color='black',linestyle=':')
        ax_bright_err.axhline(0,color='black',linestyle='-')
        ax_bright_err.axhline(0.5,color='black',linestyle=':')
        ax_bright_err.axhline(-0.5,color='black',linestyle=':')
        
        for p in range(AN):
            SOLPS_R.append(solps[p].GetRadCoords(RADC)[0].loc[:,JXA[p],:].values[:,0])
            SOLPS_ne.append(solps[p].PARAM['Ne'].loc[:,JXA[p],:].values)
            SOLPS_te.append(solps[p].PARAM['Te'].loc[:,JXA[p],:].values)
            SOLPS_nn.append(solps[p].PARAM['NeuDen'].loc[:,JXA[p],:].values)
            SOLPS_emiss.append(solps[p].PARAM['LyaEmissW'].loc[:,JXA[p],:].values)
            
            Interp_nn=interp1d(Rnn[p],nn[p][1])
            Interp_emiss=interp1d(Remiss[p],bright[p][emiss_idx]) 
            Interp_bright=interp1d(Rbright[p],bright[p][1])
            
            i = np.min([x for x in range(len(SOLPS_R[p])) if SOLPS_R[p][x] > Rnn[p][0]])
            Interp_Exp_nn.append(Interp_nn(SOLPS_R[p][i:]))
            j = np.min([x for x in range(len(SOLPS_R[p])) if SOLPS_R[p][x] > Remiss[p][0]])                     
            Interp_Exp_emiss.append(Interp_emiss(SOLPS_R[p][j:]))
            k = np.min([x for x in range(len(RLYMID[p])) if RLYMID[p][x] > Rbright[p][0]])
            Interp_Exp_bright.append(Interp_bright(RLYMID[p][k:]))
            
            nn_res.append(np.abs(Interp_Exp_nn[p]-SOLPS_nn[p][i:,1]))
            nn_res_abs.append(nn_res[p]/Interp_Exp_nn[p])
            
            emiss_res.append(np.abs(Interp_Exp_emiss[p]-SOLPS_emiss[p][j:,1]))
            emiss_res_abs.append(emiss_res[p]/Interp_Exp_emiss[p])

            bright_res.append(np.abs(Interp_Exp_bright[p]-LYMID[p][k:,1]))
            bright_res_abs.append(bright_res[p]/Interp_Exp_bright[p])            
            
            ne_plus.append((SOLPS_ne[p][:,2]-SOLPS_ne[p][:,1])/SOLPS_ne[p][:,1])
            ne_minus.append((SOLPS_ne[p][:,0]-SOLPS_ne[p][:,1])/SOLPS_ne[p][:,1])
            
            te_plus.append((SOLPS_te[p][:,2]-SOLPS_te[p][:,1])/SOLPS_te[p][:,1])
            te_minus.append((SOLPS_te[p][:,0]-SOLPS_te[p][:,1])/SOLPS_te[p][:,1])
            
            nn_plus.append((SOLPS_nn[p][:,2]-SOLPS_nn[p][:,1])/SOLPS_nn[p][:,1])
            nn_minus.append((SOLPS_nn[p][:,0]-SOLPS_nn[p][:,1])/SOLPS_nn[p][:,1])
            
            emiss_plus.append((SOLPS_emiss[p][:,2]-SOLPS_emiss[p][:,1])/SOLPS_emiss[p][:,1])
            emiss_minus.append((SOLPS_emiss[p][:,0]-SOLPS_emiss[p][:,1])/SOLPS_emiss[p][:,1])
            
            bright_plus.append((LYMID[p][:,2]-LYMID[p][:,1])/LYMID[p][:,1])
            bright_minus.append((LYMID[p][:,0]-LYMID[p][:,1])/LYMID[p][:,1])
            
            ax_nn.plot(SOLPS_R[p][i:],nn_res_abs[p],Styles[p],label=Labels[p])
            ax_emiss.plot(SOLPS_R[p][j:],emiss_res_abs[p],Styles[p],label=Labels[p])
            ax_bright.plot(RLYMID[p][k:],bright_res_abs[p],Styles[p],label=Labels[p])

            ax_ne_err.plot(SOLPS_R[p],ne_plus[p],Plus_Style[p],SOLPS_R[p],ne_minus[p],Minus_Style[p],label=Labels[p])
            ax_te_err.plot(SOLPS_R[p],te_plus[p],Plus_Style[p],SOLPS_R[p],te_minus[p],Minus_Style[p],label=Labels[p])
            ax_nn_err.plot(SOLPS_R[p],nn_plus[p],Plus_Style[p],SOLPS_R[p],nn_minus[p],Minus_Style[p],label=Labels[p])
            ax_emiss_err.plot(SOLPS_R[p],emiss_plus[p],Plus_Style[p],SOLPS_R[p],emiss_minus[p],Minus_Style[p],label=Labels[p])
            ax_bright_err.plot(RLYMID[p],bright_plus[p],Plus_Style[p],RLYMID[p],bright_minus[p],Minus_Style[p],label=Labels[p])
            
            ax_nn.set_title('Normalized Neutral Density Residuals',pad=-10.0)
            ax_nn.set_xlabel('$\Psi_n$')
            ax_nn.set_ylabel('Normalized Residual')
            ax_nn.legend()
            ax_emiss.set_title('Normalized Lyman-alpha Emissivity Residuals')                
            ax_emiss.set_xlabel('$\Psi_n$')
            ax_emiss.set_ylabel('Normalized Residual')
            ax_emiss.legend()
            ax_bright.set_title('Normalized Lyman-alpha Brightness Residuals')
            ax_bright.set_xlabel('$\Psi_n$')
            ax_bright.set_ylabel('Normalized Residual')
            ax_bright.legend()
            
            ax_ne_err.set_title('Electron Density Percent Variation')
            ax_ne_err.set_xlabel('$\Psi_n$')
            ax_ne_err.set_ylabel('Percent Variation')
            ax_ne_err.legend()
            ax_te_err.set_title('Electron Temperature Percent Variation')
            ax_te_err.set_xlabel('$\Psi_n$')
            ax_te_err.set_ylabel('Percent Variation')
            #ax_te_err.legend()
            ax_nn_err.set_title('Neutral Density Percent Variation')
            ax_nn_err.set_xlabel('$\Psi_n$')
            ax_nn_err.set_ylabel('Percent Variation')
            #ax_nn_err.legend()
            ax_emiss_err.set_title('Lyman-alpha Emissivity Percent Variation')                
            ax_emiss_err.set_xlabel('$\Psi_n$')
            ax_emiss_err.set_ylabel('Percent Variation')
            #ax_emiss_err.legend()
            ax_bright_err.set_title('Lyman-alpha Brightness Percent Variation')
            ax_bright_err.set_xlabel('$\Psi_n$')
            ax_bright_err.set_ylabel('Percent Variation')
            #ax_bright_err.legend()
            
            
            
            