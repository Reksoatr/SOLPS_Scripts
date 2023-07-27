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
import matplotlib.ticker as mtick
import equilibrium as eq
from scipy.interpolate import interp1d
from TOOLS import SET_WDIR, R2PsiN, PsiN2R, gaussian_shading, ErrorQuant
from SOLPS_Plotter import SOLPSPLOT
from PARAMDICT import EireneDict

plt.rc('font',size=30)
plt.rc('lines',linewidth=2,markersize=7)

Plot=True
Error_Analysis=True
HydrContrib=True

CoreTS=True
Filter=True
EIR_Data=True

RADC='radial' #'psin' #
EMISS='tomo' #'tree'# 7/7/23 - Now hard-wired for tomo data

SHOTS = ['1100308004','1080416025','1100305023'] #,
ATTEMPTS = [['14Rf0.35_split2','14Rf0.7_split2','14Rf1.05_split2'],
            ['18Rf0.3_split2','18Rf0.6_split2','18Rf0.9_split2'],
            ['24Rf1.0_split2','24Rf2.0_split2','24Rf3.0_split2']]

AN=len(SHOTS)

LCFS_Te=[83,90,93]
JXA=[41,41,55]
JXI=[59,58,38]
LYMID_PolSlc=[41,41,56]
PsinOffset=[0,0.02,-0.014]
RadOffset=[0,0.005,-0.007]
WallShadow=[0.89089,0.89060,0.90130]
LyaBrightErr = 0.05

if EMISS=='tree':
    Remiss_idx=2
    emiss_idx=3
elif EMISS=='tomo':
    Remiss_idx=4
    emiss_idx=5
    
## Load Experimental Data ##

BASE,TOP = SET_WDIR('','')

GFCMOD='{}gfileProcessing/cmod_files/'.format(TOP)                                              #gfile directory under TOP directory

GFiles=[]
bright=[]
lyabrighterr=[]
nn=[]
nnlogerr=[]
TS_Data=[]

for A in SHOTS:
    GFiles.append(eq.equilibrium(gfile='{}g{}.01000'.format(GFCMOD,A)))                         #Load gfiles into list
    
    with open('{}Brightness_Profiles/lya_brightness_{}_v2.pkl'.format(GFCMOD,A),'rb') as file1:
        bright.append(pkl.load(file1))                                                          #Experimental Ly-a brightness and emissivity data

    lyabrighterr.append(LyaBrightErr*bright[-1][1])   #Assume 5% nominal uncertainty in Ly-a Brightness measurements

    with open('{}Brightness_Profiles/lyman_data_{}.pkl'.format(GFCMOD,A),'rb') as file2:
        nn.append(pkl.load(file2))                                                              #Experimental neutral density profile data

    Mask=np.ma.masked_invalid(nn[-1][2]).mask                   #Mask any NaN values
     
    nn[-1][0]=np.ma.array(nn[-1][0], mask=Mask).compressed()         #Convert neutral density data from units of cm^-3 to m^-3
    nn[-1][1]=1e6*np.ma.array(nn[-1][1], mask=Mask).compressed()
    nn[-1][2]=1e6*np.ma.array(nn[-1][2], mask=Mask).compressed()    
    
    nnlogerr.append(0.434*nn[-1][2]/nn[-1][1])    #Process neutral density error for semilogy plotting
    
    

if CoreTS:                                                                      #Load Core TS Data; adapted from TS_Plotter.py
    
    for n, i in enumerate(SHOTS):
        with open('{}{}_CORE.pkl'.format(GFCMOD,i),'rb') as f:
            u = pkl._Unpickler(f)
            u.encoding = 'latin1'
            CTS = u.load()
            #CTS=pkl.load(f)
        
        TS_Data.append({})
        TS_Data[n]['ne']=np.array([CTS[0],CTS[1]*1e20,CTS[2]*1e20])
        TS_Data[n]['te']=np.array([CTS[3],CTS[4]*1e3,CTS[5]*1e3])
        
        TS_Data[n]['ne']=TS_Data[n]['ne'][:,TS_Data[n]['ne'][0,:].argsort()]
        TS_Data[n]['te']=TS_Data[n]['te'][:,TS_Data[n]['te'][0,:].argsort()]            
        
        if LCFS_Te!=0:
            index=(np.abs(TS_Data[n]['te'][1]-LCFS_Te[n])).argmin()
            PsinOffset[n]=1-TS_Data[n]['te'][0][index]
            print('Psin offset for shot {} is {}'.format(i,PsinOffset[n]))
            
        if Filter:
            for v in ['ne','te']:
                for t in range(1,len(TS_Data[n][v][1])):
                    if ((TS_Data[n][v][1][t]-TS_Data[n][v][2][t]) > 2*(TS_Data[n][v][1][t-1]) 
                        or (TS_Data[n][v][1][t]+TS_Data[n][v][2][t]) < 0.5*(TS_Data[n][v][1][t-1])):
                        
                        TS_Data[n][v][0][t]=TS_Data[n][v][1][t]=TS_Data[n][v][2][t]=np.nan
                
                TS_Data[n][v]=np.array([np.ma.masked_invalid(TS_Data[n][v][0]).compressed(),
                                        np.ma.masked_invalid(TS_Data[n][v][1]).compressed(),
                                        np.ma.masked_invalid(TS_Data[n][v][2]).compressed()])

## Load SOLPS simulation data ##

solps=[SOLPSPLOT(SHOTS[i], ATTEMPTS[i], Parameters=['Ne','Te','DN','KYE','KYI','NeuDen','LyaEmissW','IonFlx','IonPol'], 
                 Markers=False, JXA=JXA[i], JXI=JXI[i], TimeRange=[0.95,1.05], PlotScheme=['b--','b-','b--'], 
                 PsinOffset=PsinOffset[i], RadOffset=RadOffset[i],EXP=(not CoreTS)) 
                 for i in range(AN)]    #Load nominal SOLPS B2-matrix data

with open('{}/NeuDenEIR.pkl'.format(GFCMOD),'rb') as file5:
    NeuDenEIR=pkl.load(file5)

with open('{}/LyaEmissEIR.pkl'.format(GFCMOD),'rb') as file4:
    LyaEmissEIR=pkl.load(file4)

with open('{}/Chords/lya_coords_new.pkl'.format(GFCMOD),'rb') as file3:
    LYMID_coords=pkl.load(file3)   #Load in LYMID XYZ coordinates from pickle file
    
RLYMID0=np.flip(np.sqrt(LYMID_coords['tangent']['X']**2 + LYMID_coords['tangent']['Y']**2))    #Calculate R coordinates of tangent points (R^2=X^2+Y^2)
ZLYMID0=np.flip(LYMID_coords['start']['Z'])                                                    #Pull Z coordinates of tangent points from start coordinates     

LYMID = np.zeros((AN,len(RLYMID0),3))
LYMID_FULL = np.zeros((AN,len(RLYMID0),7))

for m,a in enumerate(SHOTS):
    for n,i in enumerate(ATTEMPTS[m]):
        LYMID[m,:,n]=np.loadtxt('{}SOLPS_2d_prof/cmod/{}home/Attempt{}/EirOutput/LyaBrightEIR{}'.format(TOP,a,i,i),
                        skiprows=2)[:,1]*100*100*2.052E-17     # Load Brightness profile simulation data and convert units from PHOTONS/S/CM^2/SR to W/M^2
        
        if n==1:
            LYMID_FULL[m,:,:]=np.loadtxt('{}SOLPS_2d_prof/cmod/{}home/Attempt{}/EirOutput/LyaBrightEIR_FULL{}'.format(TOP,a,i,i),
                            skiprows=2)[:,1:]*100*100*2.052E-17 # Load Full Lya Brightness Contributions Data
            LYMID_FULL[m,:,0]=np.sum(LYMID_FULL[m,:,1:],axis=1) # Calculate SUM OVER CONTRIBUTIONS

Rbright=[]
Remiss=[]    
Rnn=[]
Sep=[]
RLYMID=[]

if RADC=='radial':                                          # Convert all data coordinates to R (meters)
    for m in range(AN):
        
        Rbright.append(bright[m][0])                        # Brightness profile coordinates already given in R
        Remiss.append(bright[m][Remiss_idx])                # Emissivity profile coordinates already given in R
        Rnn.append(PsiN2R(GFiles[m],nn[m][0],Z=ZLYMID0[0])) # Convert neutral density coordinates from PsiN to R with gfile
        Sep.append(PsiN2R(GFiles[m],1.0,Z=ZLYMID0[0]))      # Calculate R coordinate of Separatrix from gfile
        RLYMID.append(RLYMID0)                              # SOLPS brightness profile coordinates already given in R

        '''
        Rnn[m]=np.array([[j for j,k in zip(*GFiles[m].get_fluxsurface(i)) if j > GFiles[m].axis.r and k==0] 
                         for i in nn[m][0]])[:,0]  
        Sep[m]=GFiles[m].get_fluxsurface(1.0)[0].max()
        '''
        
    Rlabel=r'R [m]'
    
elif RADC=='psin':                                                          # Convert all data coordinates to PsiN
    for m in range(AN):
        
        Rbright.append(R2PsiN(GFiles[m],bright[m][0],Z=ZLYMID0))
        Remiss.append(R2PsiN(GFiles[m],bright[m][Remiss_idx],Z=ZLYMID0))
        Rnn.append(nn[m][0])
        Sep.append(1)
        RLYMID.append(R2PsiN(GFiles[m],RLYMID0,Z=ZLYMID0))
    
    Rlabel=r'$\Psi_n$'

### PLOTTING ROUTINE ###

if Plot:
    
    figA,axA=plt.subplots(3,AN,sharex='col',sharey='row')
    figB,axB=plt.subplots(3,AN,sharex='col',sharey='row')
    
    if axA.shape==(6,):
        axA=np.array([axA]).T #Make sure axA always has 2 dimensions for indexing sake
    
    for p in range(AN):
        
        ## Plasma Spread (axA)
        
        solps[p].RadProf('DN',RADC='psin',AX=axA[0,p],PlotScheme=['g--','g-','g--'],Publish=['$\pm 50\%\;D_n$','$D_n$',None])
        solps[p].RadProf('KYE',RADC='psin',AX=axA[0,p],PlotScheme=['m--','m-','m--'],Publish=['$\pm 50\%\;\chi_e$','$\chi_e$',None])
        axA[0,p].legend()
        
        solps[p].RadProf('Te',RADC='psin',AX=axA[1,p],PlotScheme=['b--','b-','b--'],Publish=[None,'Base','$\pm50\%$'])
        #axA[1,p].get_lines()[0].set_color('blue')
        #axA[1,p].axhline(LCFS_Te[0],linestyle=':',color='orange',label='2PM Te')
        if CoreTS:
            axA[1,p].errorbar(TS_Data[p]['te'][0]+PsinOffset[p],TS_Data[p]['te'][1],
                                         yerr=TS_Data[p]['te'][2],fmt='r.',capsize=5,label='TS')
        axA[1,p].legend()
        
        solps[p].RadProf('Ne',RADC='psin',AX=axA[2,p],PlotScheme=['b--','b-','b--'],Publish=[None,'Base','$\pm50\%$'])
        #axA[2,p].get_lines()[0].set_color('blue')
        if CoreTS:
            axA[2,p].errorbar(TS_Data[p]['ne'][0]+PsinOffset[p],TS_Data[p]['ne'][1],
                                         yerr=TS_Data[p]['ne'][2],fmt='r.',capsize=5,label='TS')
        axA[2,p].legend()
        
        axA[0,p].axvline(1.0,linestyle=':',color='k')
        axA[1,p].axvline(1.0,linestyle=':',color='k')
        axA[2,p].axvline(1.0,linestyle=':',color='k')
        
        axA[0,p].set_title(SHOTS[p], pad=10.0)
        axA[0,p].set_xlabel('')
        axA[0,p].set_ylabel('')
        axA[1,p].set_title('')
        axA[1,p].set_ylabel('')
        axA[1,p].set_xlabel('')
        axA[2,p].set_title('')
        axA[2,p].set_ylabel('')
        axA[2,p].set_xlabel('$\Psi_n$')
        
        xA_min=axA[0,p].get_lines()[0].get_xdata().min()
        xA_max=axA[0,p].get_lines()[0].get_xdata().max()
        
        axA[0,p].set_xlim(xA_min,xA_max)
        
        ## Neutral Spread (axB)
        
        # Neutral D Density Plots
        
        if EIR_Data:
            SortKeyND=NeuDenEIR[p]['Xcoords']['R0'].argsort()
            axB[0,p].semilogy(NeuDenEIR[p]['Xcoords']['R0'][SortKeyND],
                          np.power(10,NeuDenEIR[p][ATTEMPTS[p][0]]['0']['Data'][SortKeyND]),
                          color='blue',marker=3,linestyle='None',label='$\pm50\%$')
            axB[0,p].semilogy(NeuDenEIR[p]['Xcoords']['R0'][SortKeyND],
                          np.power(10,NeuDenEIR[p][ATTEMPTS[p][1]]['0']['Data'][SortKeyND]),
                          color='blue',marker='*',linestyle='--',label='Base')
            axB[0,p].semilogy(NeuDenEIR[p]['Xcoords']['R0'][SortKeyND],
                          np.power(10,NeuDenEIR[p][ATTEMPTS[p][2]]['0']['Data'][SortKeyND]),
                          color='blue',marker=2,linestyle='None')
        else:
            solps[p].RadProf('NeuDen',RADC=RADC,AX=axB[0,p],LOG10=2, Publish=['$\pm50\%$','Base',None],PolSlc=LYMID_PolSlc[p])
        axB[0,p].semilogy(Rnn[p],nn[p][1]-nn[p][2],'r--',linewidth=1)
        axB[0,p].semilogy(Rnn[p],nn[p][1],'r.--',linewidth=1,label='LYMID')
        axB[0,p].semilogy(Rnn[p],nn[p][1]+nn[p][2],'r--',linewidth=1)

        # Ly-alpha Emissivity Plots

        if EIR_Data:
            SortKeyLE=LyaEmissEIR[p]['Xcoords']['R0'].argsort()
            axB[1,p].plot(LyaEmissEIR[p]['Xcoords']['R0'][SortKeyLE],
                          np.power(10,LyaEmissEIR[p][ATTEMPTS[p][0]]['0']['Data'][SortKeyLE]),
                          color='blue',marker=3,linestyle='None',label='$\pm50\%$')
            axB[1,p].plot(LyaEmissEIR[p]['Xcoords']['R0'][SortKeyLE],
                          np.power(10,LyaEmissEIR[p][ATTEMPTS[p][1]]['0']['Data'][SortKeyLE]),
                          color='blue',marker='*',linestyle='--',label='$\pm50\%$')
            axB[1,p].plot(LyaEmissEIR[p]['Xcoords']['R0'][SortKeyLE],
                          np.power(10,LyaEmissEIR[p][ATTEMPTS[p][2]]['0']['Data'][SortKeyLE]),
                          color='blue',marker=2,linestyle='None',label='$\pm50\%$')
        else:
            solps[p].RadProf('LyaEmissW',RADC=RADC,AX=axB[1,p], Publish=['$\pm50\%$','Base',None],PolSlc=LYMID_PolSlc[p])
        axB[1,p].plot(Remiss[p],bright[p][5]-bright[p][6],'r--',linewidth=1)
        axB[1,p].plot(Remiss[p],bright[p][5],'r.--',linewidth=1,label='LYMID')
        axB[1,p].plot(Remiss[p],bright[p][5]+bright[p][6],'r--',linewidth=1)
        
        # Ly-alpha Brightness Plots
        
        axB[2,p].plot(RLYMID[p],LYMID[p][:,0],color='blue',marker=3,linestyle='None',markersize='7.0',label='+/-50%')
        axB[2,p].plot(RLYMID[p],LYMID[p][:,1],color='blue',marker='*',linestyle='--',markersize='7.0',label='Base')
        axB[2,p].plot(RLYMID[p],LYMID[p][:,2],color='blue',marker=2,linestyle='None',markersize='7.0',label=None)
        axB[2,p].errorbar(Rbright[p],bright[p][1],lyabrighterr[p],fmt='r.--',elinewidth=1,capsize=5,label='LYMID')
        
        # Separatrix and Limiter Shadow Markers
        
        axB[0,p].axvline(Sep[p],linestyle=':',color='k')
        axB[0,p].axvspan(WallShadow[p],axB[0,p].get_xlim()[1],alpha=0.4,color='grey',zorder=0.1)
        axB[1,p].axvline(Sep[p],linestyle=':',color='k')
        axB[1,p].axvspan(WallShadow[p],axB[1,p].get_xlim()[1],alpha=0.4,color='grey',zorder=0.1)
        axB[2,p].axvline(Sep[p],linestyle=':',color='k')
        axB[2,p].axvspan(WallShadow[p],axB[2,p].get_xlim()[1],alpha=0.4,color='grey',zorder=0.1)
        
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
        
        RGrid=np.linspace(GF.axis.r,GF.R.max(),100)
        PsiNGrid=R2PsiN(GF,RGrid)
        
        def PsiN2Rax(x):
            return np.interp(x, PsiNGrid, RGrid)
        
        def R2PsiNax(x):
            return np.interp(x, RGrid, PsiNGrid)
        
        ax2A=axA[0,p].secondary_xaxis('top',functions=(PsiN2Rax, R2PsiNax))
        ax2A.set_xlabel('R [m]')
        
        if RADC=='psin':
            ax2B=axB[0,p].secondary_xaxis('top',functions=(PsiN2Rax, R2PsiNax))
            ax2B.set_xlabel('R [m]')
        
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
        
    axB[0,1].text(0.6,0.03,'Limiter Shadow',color='0.5',transform=axB[0,1].transAxes)
    axB[1,1].legend()
    
    figA.subplots_adjust(wspace=0, hspace=0)
    figB.subplots_adjust(wspace=0, hspace=0)
    
    if HydrContrib:
        
        ColorCyc=['darkblue','cyan','blue','lightblue','turquoise','teal','purple']
        
        figHE, axHE = plt.subplots(ncols=3,sharey=True,sharex=True)
        figHB, axHB = plt.subplots(ncols=3,sharey=True,sharex=True)
        
        for p in range(AN):
            SortKeyLE=LyaEmissEIR[p]['Xcoords']['R0'].argsort()
            
            for cc in range(7):
                axHE[p].plot(LyaEmissEIR[p]['Xcoords']['R0'][SortKeyLE],
                              np.power(10,LyaEmissEIR[p][ATTEMPTS[p][1]][str(cc)]['Data'][SortKeyLE]),
                              label=EireneDict[str(cc)]['sym'],linestyle='-',linewidth=2,color=ColorCyc[cc])
                axHE[p].axvline(Sep[p],linestyle=':',color='k')
                axHE[p].set_xlabel('R [m]')
                axHE[p].set_title(SHOTS[p])
                
                axHB[p].plot(RLYMID0,LYMID_FULL[p,:,cc],
                              label=EireneDict[str(cc)]['sym'],linestyle='-',linewidth=2,color=ColorCyc[cc])
                axHB[p].axvline(Sep[p],linestyle=':',color='k')
                axHB[p].set_xlabel('R [m]')
                axHB[p].set_title(SHOTS[p])
        
            axHE[p].plot(Remiss[p],bright[p][5]-bright[p][6],'r--',linewidth=1)
            axHE[p].plot(Remiss[p],bright[p][5],'r.--',linewidth=1,label='LYMID')
            axHE[p].plot(Remiss[p],bright[p][5]+bright[p][6],'r--',linewidth=1)
            
            axHB[p].errorbar(Rbright[p],bright[p][1],lyabrighterr[p],fmt='r.--',elinewidth=1,capsize=5,label='LYMID')           
        
        axHE[0].set_yticklabels([])
        axHE2=axHE[0].secondary_yaxis('left',functions=(lambda x: x/1000, lambda x: x*1000))
        axHE2.set_ylabel('Ly-alpha Emissivity ($kW/m^3$)')
        axHB[0].set_yticklabels([])
        axHB2=axHB[0].secondary_yaxis('left',functions=(lambda x: x/1000, lambda x: x*1000))
        axHB2.set_ylabel('Ly-alpha Brightness ($kW/m^2$)')
        
        axHE[1].legend()
        axHB[1].legend()     
                
    
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
        SOLPS_Psin=[]
        
        Interp_ne=[]
        Interp_te=[]
        Interp_nn=[]
        Interp_emiss=[]
        Interp_bright=[]
        
        ne_res=[]
        te_res=[]
        nn_res=[]
        emiss_res=[]
        bright_res=[]
        
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
        
        fig_res1,ax_res1=plt.subplots(nrows=2,sharex=True)
        ax_ne=ax_res1[0]
        ax_te=ax_res1[1]
        
        fig_res2,ax_res2=plt.subplots(nrows=3,sharex=True)
        ax_nn=ax_res2[0]
        ax_emiss=ax_res2[1]
        ax_bright=ax_res2[2]
        
        fig_err1,ax_err1=plt.subplots(nrows=2,sharex=True)
        ax_ne_err=ax_err1[0]
        ax_te_err=ax_err1[1]
        
        fig_err2,ax_err2=plt.subplots(nrows=3,sharex=True)
        ax_nn_err=ax_err2[0]
        ax_emiss_err=ax_err2[1]
        ax_bright_err=ax_err2[2]
        
        ax_ne.axhline(0,color='black',linestyle='--')
        ax_te.axhline(0,color='black',linestyle='--')
        ax_nn.axhline(0,color='black',linestyle='--')
        ax_emiss.axhline(0,color='black',linestyle='--')
        ax_bright.axhline(0,color='black',linestyle='--')
        
        ax_ne_err.axhline(0,color='black',linestyle='--')
        ax_ne_err.axhline(0.5,color='magenta',linestyle='dashdot')
        ax_ne_err.axhline(-0.5,color='magenta',linestyle='dashdot')
        
        ax_te_err.axhline(0,color='black',linestyle='--')
        ax_te_err.axhline(0.5,color='magenta',linestyle='dashdot')
        ax_te_err.axhline(-0.5,color='magenta',linestyle='dashdot')
        
        ax_nn_err.axhline(0,color='black',linestyle='--')
        ax_nn_err.axhline(0.5,color='magenta',linestyle='dashdot')
        ax_nn_err.axhline(-0.5,color='magenta',linestyle='dashdot')
        
        ax_emiss_err.axhline(0,color='black',linestyle='--')
        ax_emiss_err.axhline(0.5,color='magenta',linestyle='dashdot')
        ax_emiss_err.axhline(-0.5,color='magenta',linestyle='dashdot')
        
        ax_bright_err.axhline(0,color='black',linestyle='--')
        ax_bright_err.axhline(0.5,color='magenta',linestyle='dashdot')
        ax_bright_err.axhline(-0.5,color='magenta',linestyle='dashdot')
        
        for p in range(AN):
            
            SOLPS_R.append(solps[p].GetRadCoords('radial')[0].loc[:,JXA[p],:].values[:,0])
            SOLPS_Psin.append(solps[p].GetRadCoords('psin')[0].loc[:,JXA[p],:].values[:,0])
            
            SOLPS_ne.append(solps[p].PARAM['Ne'].loc[:,JXA[p],:].values)
            SOLPS_te.append(solps[p].PARAM['Te'].loc[:,JXA[p],:].values)
            
            if EIR_Data:
                SortKeyND = NeuDenEIR[p]['Xcoords']['R0'].argsort()
                SortKeyLE = LyaEmissEIR[p]['Xcoords']['R0'].argsort()
                
                SOLPS_nn_profile = interp1d(NeuDenEIR[p]['Xcoords']['R0'][SortKeyND],
                              np.power(10,NeuDenEIR[p][ATTEMPTS[p][1]]['0']['Data'][SortKeyND]))
                SOLPS_emiss_profile = interp1d(LyaEmissEIR[p]['Xcoords']['R0'][SortKeyLE],
                              np.power(10,LyaEmissEIR[p][ATTEMPTS[p][1]]['0']['Data'][SortKeyLE]))
                
            else:
                SOLPS_nn.append(solps[p].PARAM['NeuDen'].loc[:,LYMID_PolSlc[p],:].values)
                SOLPS_emiss.append(solps[p].PARAM['LyaEmissW'].loc[:,LYMID_PolSlc[p],:].values)
            
            SOLPS_bright_profile=interp1d(RLYMID0,LYMID[p,:,1])
        
            
            Data_ne_idx = [n for n in range(len(TS_Data[p]['ne'][0])) if 
                            TS_Data[p]['ne'][0][n]+PsinOffset[p] > SOLPS_Psin[p].min() and
                            TS_Data[p]['ne'][0][n]+PsinOffset[p] < SOLPS_Psin[p].max()]
            Interp_ne.append(solps[p].Interp1DRad('Ne', TS_Data[p]['ne'][0][Data_ne_idx]+PsinOffset[p], Attempt=ATTEMPTS[p][1]))
            
            Data_te_idx = [n for n in range(len(TS_Data[p]['te'][0])) if 
                            TS_Data[p]['te'][0][n]+PsinOffset[p] > SOLPS_Psin[p].min() and
                            TS_Data[p]['te'][0][n]+PsinOffset[p] < SOLPS_Psin[p].max()]
            Interp_te.append(solps[p].Interp1DRad('Te', TS_Data[p]['te'][0][Data_te_idx]+PsinOffset[p], Attempt=ATTEMPTS[p][1]))
            
            Data_nn_idx = [n for n in range(len(Rnn[p])) if 
                           Rnn[p][n] > NeuDenEIR[p]['Xcoords']['R0'].min() and 
                           Rnn[p][n] < NeuDenEIR[p]['Xcoords']['R0'].max()]
            Interp_nn.append(SOLPS_nn_profile(Rnn[p][Data_nn_idx]))
            #Interp_nn.append(solps[p].Interp1DRad('NeuDen', Rnn[p][Data_nn_idx], Attempt=ATTEMPTS[p][1], PolSlc=LYMID_PolSlc[p]))
            
            Data_emiss_idx = [n for n in range(len(Remiss[p])) if 
                              Remiss[p][n] > LyaEmissEIR[p]['Xcoords']['R0'].min() and 
                              Remiss[p][n] < LyaEmissEIR[p]['Xcoords']['R0'].max()]
            Interp_emiss.append(SOLPS_emiss_profile(Remiss[p][Data_emiss_idx]))            
            #Interp_emiss.append(solps[p].Interp1DRad('LyaEmissW', Remiss[p][Data_emiss_idx], RADC='radial', Attempt=ATTEMPTS[p][1], PolSlc=LYMID_PolSlc[p]))
            
            Data_bright_idx = [n for n in range(len(Rbright[p])) if 
                               Rbright[p][n] > RLYMID0.min() and 
                               Rbright[p][n] < RLYMID0.max()]
            Interp_bright.append(SOLPS_bright_profile(Rbright[p][Data_bright_idx]))

            ne_res.append(ErrorQuant(TS_Data[p]['ne'][1][Data_ne_idx], Interp_ne[p], TS_Data[p]['ne'][2][Data_ne_idx], name='electron density'))

            te_res.append(ErrorQuant(TS_Data[p]['te'][1][Data_te_idx], Interp_te[p], TS_Data[p]['te'][2][Data_te_idx], name='electron temperature'))            

            nn_res.append(ErrorQuant(nn[p][1][Data_nn_idx], Interp_nn[p], nn[p][2][Data_nn_idx], name='neutral density'))
            
            emiss_res.append(ErrorQuant(bright[p][5][Data_emiss_idx], Interp_emiss[p], bright[p][6][Data_emiss_idx], name='Ly-a emissivity'))

            bright_res.append(ErrorQuant(bright[p][1][Data_bright_idx], Interp_bright[p], lyabrighterr[p][Data_bright_idx], name='Ly-a brightness'))   
            
            
            ne_plus.append((SOLPS_ne[p][:,2]-SOLPS_ne[p][:,1])/SOLPS_ne[p][:,1])
            ne_minus.append((SOLPS_ne[p][:,0]-SOLPS_ne[p][:,1])/SOLPS_ne[p][:,1])
            
            te_plus.append((SOLPS_te[p][:,2]-SOLPS_te[p][:,1])/SOLPS_te[p][:,1])
            te_minus.append((SOLPS_te[p][:,0]-SOLPS_te[p][:,1])/SOLPS_te[p][:,1])
            
            '''
            nn_plus.append((SOLPS_nn[p][:,2]-SOLPS_nn[p][:,1])/SOLPS_nn[p][:,1])
            nn_minus.append((SOLPS_nn[p][:,0]-SOLPS_nn[p][:,1])/SOLPS_nn[p][:,1])
            
            emiss_plus.append((SOLPS_emiss[p][:,2]-SOLPS_emiss[p][:,1])/SOLPS_emiss[p][:,1])
            emiss_minus.append((SOLPS_emiss[p][:,0]-SOLPS_emiss[p][:,1])/SOLPS_emiss[p][:,1])
            '''
            
            nn_plus.append((np.power(10,NeuDenEIR[p][ATTEMPTS[p][2]]['0']['Data'][SortKeyND])-
                           np.power(10,NeuDenEIR[p][ATTEMPTS[p][1]]['0']['Data'][SortKeyND]))/
                           np.power(10,NeuDenEIR[p][ATTEMPTS[p][1]]['0']['Data'][SortKeyND]))
            nn_minus.append((np.power(10,NeuDenEIR[p][ATTEMPTS[p][0]]['0']['Data'][SortKeyND])-
                           np.power(10,NeuDenEIR[p][ATTEMPTS[p][1]]['0']['Data'][SortKeyND]))/
                           np.power(10,NeuDenEIR[p][ATTEMPTS[p][1]]['0']['Data'][SortKeyND]))
            
            emiss_plus.append((np.power(10,LyaEmissEIR[p][ATTEMPTS[p][2]]['0']['Data'][SortKeyLE])-
                           np.power(10,LyaEmissEIR[p][ATTEMPTS[p][1]]['0']['Data'][SortKeyLE]))/
                           np.power(10,LyaEmissEIR[p][ATTEMPTS[p][1]]['0']['Data'][SortKeyLE]))
            emiss_minus.append((np.power(10,LyaEmissEIR[p][ATTEMPTS[p][0]]['0']['Data'][SortKeyLE])-
                           np.power(10,LyaEmissEIR[p][ATTEMPTS[p][1]]['0']['Data'][SortKeyLE]))/
                           np.power(10,LyaEmissEIR[p][ATTEMPTS[p][1]]['0']['Data'][SortKeyLE]))
            
            bright_plus.append((LYMID[p][:,2]-LYMID[p][:,1])/LYMID[p][:,1])
            bright_minus.append((LYMID[p][:,0]-LYMID[p][:,1])/LYMID[p][:,1])
            
            
            ax_ne.plot(TS_Data[p]['ne'][0][Data_ne_idx]+PsinOffset[p],ne_res[p]['norm_res'],Styles[p],label=Labels[p])
            ax_te.plot(TS_Data[p]['te'][0][Data_te_idx]+PsinOffset[p],te_res[p]['norm_res'],Styles[p],label=Labels[p])
            
            ax_nn.plot(Rnn[p][Data_nn_idx],nn_res[p]['norm_res'],Styles[p],label=Labels[p])
            ax_emiss.plot(Remiss[p][Data_emiss_idx],emiss_res[p]['norm_res'],Styles[p],label=Labels[p])
            ax_bright.plot(Rbright[p][Data_bright_idx],bright_res[p]['norm_res'],Styles[p],label=Labels[p])


            ax_ne_err.plot(SOLPS_Psin[p],ne_plus[p],Plus_Style[p],
                           SOLPS_Psin[p],ne_minus[p],Minus_Style[p],label=Labels[p])
            ax_te_err.plot(SOLPS_Psin[p],te_plus[p],Plus_Style[p],
                           SOLPS_Psin[p],te_minus[p],Minus_Style[p],label=Labels[p])
            
            ax_nn_err.plot(NeuDenEIR[p]['Xcoords']['R0'][SortKeyND],nn_plus[p],Plus_Style[p],
                           NeuDenEIR[p]['Xcoords']['R0'][SortKeyND],nn_minus[p],Minus_Style[p],label=Labels[p])
            ax_emiss_err.plot(LyaEmissEIR[p]['Xcoords']['R0'][SortKeyLE],emiss_plus[p],Plus_Style[p],
                              LyaEmissEIR[p]['Xcoords']['R0'][SortKeyLE],emiss_minus[p],Minus_Style[p],label=Labels[p])
            ax_bright_err.plot(RLYMID[p],bright_plus[p],Plus_Style[p],
                               RLYMID[p],bright_minus[p],Minus_Style[p],label=Labels[p])

            
            ax_ne.set_title('Radial profiles of Normalized Residuals ($R_N$)')
            ax_ne.set_xlabel('')
            ax_ne.set_ylabel('$R_N (n_e)$')
            ax_ne.axvline(1.0,color='black',linestyle=':')
            ax_ne.legend()
            ax_te.set_title('')                
            ax_te.set_xlabel('$\Psi_n$')
            ax_te.set_ylabel('$R_N (T_e)$')
            ax_te.axvline(1,color='black',linestyle=':')
            #ax_te.legend()
            
            ax_nn.set_title('Radial profiles of Normalized Residuals ($R_N$)')
            ax_nn.set_xlabel('')
            ax_nn.set_ylabel('$R_N (n_D)$')
            ax_nn.axvline(np.mean(Sep),color='black',linestyle=':')
            #ax_nn.legend()
            ax_emiss.set_title('')                
            ax_emiss.set_xlabel('')
            ax_emiss.set_ylabel('$R_N$(Ly-a Emiss.)')
            ax_emiss.axvline(np.mean(Sep),color='black',linestyle=':')
            #ax_emiss.legend()
            ax_bright.set_title('')
            ax_bright.set_xlabel('R [m]')
            ax_bright.set_ylabel('$R_N$(Ly-a Bright.)')
            ax_bright.axvline(np.mean(Sep),color='black',linestyle=':',label='Avg SEP')
            ax_bright.legend()
            
            
            ax_ne_err.set_title('Percent variation from center baseline case')
            ax_ne_err.set_xlabel('$\Psi_n$')
            ax_ne_err.set_ylabel('$n_e$')
            ax_ne_err.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0,decimals=0))
            ax_ne_err.legend()
            ax_te_err.set_title('')
            ax_te_err.set_xlabel('$\Psi_n$')
            ax_te_err.set_ylabel('$T_e$')
            ax_te_err.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0,decimals=0))
            #ax_te_err.legend()
            
            ax_nn_err.set_title('Percent variation from center baseline case')
            ax_nn_err.set_xlabel('R [m]')
            ax_nn_err.set_ylabel('$n_D$')
            ax_nn_err.yaxis.set_major_locator(mtick.MultipleLocator(0.5))
            ax_nn_err.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0,decimals=0))
            #ax_nn_err.legend()
            ax_emiss_err.set_title('')                
            ax_emiss_err.set_xlabel('R [m]')
            ax_emiss_err.set_ylabel('Ly-a Emissivity')
            ax_emiss_err.yaxis.set_major_locator(mtick.MultipleLocator(0.5))
            ax_emiss_err.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0,decimals=0))
            #ax_emiss_err.legend()
            ax_bright_err.set_title('')
            ax_bright_err.set_xlabel('R [m]')
            ax_bright_err.set_ylabel('Ly-a Brightness')
            ax_bright_err.yaxis.set_major_locator(mtick.MultipleLocator(0.5))
            ax_bright_err.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0,decimals=0))
            ax_bright_err.legend()          
            
            