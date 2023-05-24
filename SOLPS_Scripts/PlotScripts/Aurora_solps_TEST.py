# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 20:44:04 2021

@author: Richard Reksoatmodjo

Aurora SOLPS module test script
"""
from aurora.solps import solps_case
from SOLPS_Scripts.TOOLS import SET_WDIR
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib import collections as mc

SHOT='012'#'025'#
ATTEMPT=['65N']#,'66N','67N']#['48N','47N']#

Legend=['8.25TorrL']#,'39.5TorrL','77.8TorrL']#['6.77TorrL','72.2TorrL']#

Plot_Wall=False
Neutral_Energy=False

BASE,TOP = SET_WDIR('','')
GFILE='{}gfileProcessing/cmod_files/g1160718025.01209_974'.format(TOP) # g1080416025.01100'.format(TOP) #

if isinstance(ATTEMPT,str):
    ATTEMPT=[ATTEMPT]
    
#SO={}
for i in ATTEMPT:
    PATH='{}cmod/{}home/Attempt{}/Output'.format(BASE,SHOT,i)
    SO=solps_case(PATH,GFILE,case_num=i,form='full')

if not Legend:
    Legend=ATTEMPT

if Neutral_Energy:
    if len(SO) > 1:
        histfig, histax = plt.subplots(1,1)
            
    for n,i in enumerate(ATTEMPT):
        F46=SO[i].load_fort46()
        
        atom_energy=F46['edena']/F46['pdena']
        mol_energy=F46['edenm']/F46['pdenm']
        
        keV_per=stats.percentileofscore(atom_energy[~np.isnan(atom_energy)],1000)
        
        atom_energy_fig, atom_energy_ax = plt.subplots(1,3)
        mol_energy_fig, mol_energy_ax = plt.subplots(1,3)
        # Combined total neutral energy? How would that be averaged/calculated?
        
        SO[i].plot2d_eirene(F46['pdena'],ax=atom_energy_ax[0])
        SO[i].plot2d_eirene(atom_energy,ax=atom_energy_ax[1],scale='log')
        atom_energy_ax[2].hist(atom_energy,100,log=True,cumulative=0,density=0,range=(0,1000))
        atom_energy_ax[2].yaxis.grid()
        mean_atomic_energy = np.nanmean(atom_energy)
        atom_energy_90 = np.nanpercentile(atom_energy,90)
        atom_energy_ax[2].axvline(mean_atomic_energy,color='red',linewidth=2)
        #atom_energy_ax[2].axvline(np.nanpercentile(atom_energy,50),color='black',linewidth=1)
        atom_energy_ax[2].axvline(atom_energy_90,color='orange',linewidth=1)
        
        SO[i].plot2d_eirene(F46['pdenm'],ax=mol_energy_ax[0])
        SO[i].plot2d_eirene(mol_energy,ax=mol_energy_ax[1],scale='log')
        mol_energy_ax[2].hist(mol_energy,100,log=True,cumulative=0,density=0,range=(0,np.nanpercentile(mol_energy,99)))    
        mol_energy_ax[2].yaxis.grid()
        mol_energy_ax[2].axvline(np.nanmean(mol_energy),color='red',linewidth=2)
        #mol_energy_ax[2].axvline(np.nanpercentile(mol_energy,50),color='black',linewidth=1)
        mol_energy_ax[2].axvline(np.nanpercentile(mol_energy,90),color='orange',linewidth=1)
            
        atom_energy_ax[0].set_title(r'Neutral D density ($cm^{-3}$)')
        atom_energy_ax[1].set_title(r'Neutral D energy ($eV$)')
        atom_energy_ax[2].set_title(r'Neutral D energy ($eV$), Histogram (1 keV at {:.1f}%)'.format(keV_per))
        atom_energy_ax[2].set_xlabel(r'Energy ($eV$)')
        atom_energy_ax[2].set_ylabel(r'Cell Count Tally')
        atom_energy_ax[2].legend(['Mean','90th Percentile'],loc=8)
        
        mol_energy_ax[0].set_title(r'Neutral D2 density ($cm^{-3}$)')
        mol_energy_ax[1].set_title(r'Neutral D2 energy ($eV$)')
        mol_energy_ax[2].set_title(r'Neutral D2 energy ($eV$), Histogram (99%)')
        mol_energy_ax[2].set_xlabel(r'Energy ($eV$)')
        mol_energy_ax[2].set_ylabel(r'Cell Count Tally')
        mol_energy_ax[2].legend(['Mean','90th Percentile'],loc=8)
        
        atom_energy_fig.suptitle('Shot {} Attempt {} Atomic Neutral Data'.format(SHOT,i))
        mol_energy_fig.suptitle('Shot {} Attempt {} Molecular Neutral Data'.format(SHOT,i))
        
        if len(SO) > 1:
            hh=histax.hist(atom_energy,100,log=0,cumulative=0,density=1,range=(0,1000),histtype='step',label=Legend[n])
            histax.axvline(mean_atomic_energy,linewidth=2,color=hh[2][0].get_edgecolor(),label=r'{}, neutral mean energy ($\langle T_n \rangle =${:.1f}eV)'.format(Legend[n],mean_atomic_energy))
            histax.axvline(atom_energy_90,linewidth=1,color=hh[2][0].get_edgecolor(),linestyle='dashed',label=r'{}, 90th percentile neutral energy ($T_n =${:.1f}eV)'.format(Legend[n],atom_energy_90))
            
    histax.set_title(r'Neutral D energy ($eV$) Histogram')
    histax.set_xlabel(r'Energy ($eV$)')
    histax.set_ylabel(r'Normalized Cell Count Tally (log10)') 
    histax.legend()

if Plot_Wall:
    SO.plot_wall_geometry()
    
    ws=SO.WS
    
    j=0
    ID=[0]
    
    for i in range(1,len(ws)):
        if ws[i-1][-1]==ws[i][0]:
            ID.append(j)
            #print(j)
        else:
            j+=1
            ID.append(j)
            #print(j)



'''
## Old Wall Plotting Routine (now built into aurora.solps as of 3/2/21) ##
'''
if Plot_Wall:
    Wall_Seg=[]
    
    EirOut=SO.load_eirene_output(files=['fort.44'])

    RR=EirOut['fort.44']['wall_geometry'][0::2]
    ZZ=EirOut['fort.44']['wall_geometry'][1::2]
    NLIM=len(RR)//2
    
    for i in range(0,NLIM):
        line=[(RR[2*i],ZZ[2*i]),(RR[2*i+1],ZZ[2*i+1])]
        Wall_Seg.append(line)
    
    
    if SHOT=='012':
        MC=slice(0,80)
        OT=slice(80,88)
        DIV=slice(88,99)
        IT=slice(99,NLIM)
    elif SHOT=='025':
        MC=slice(0,79)
        OT=slice(79,87)
        DIV=slice(87,97)
        IT=slice(97,NLIM)
    else:
        MC=slice(0,80)
        OT=slice(80,89)
        DIV=slice(89,97)
        IT=slice(97,NLIM)
    
    MC_wall=mc.LineCollection(Wall_Seg[MC],colors='b',linewidth=2)
    OT_wall=mc.LineCollection(Wall_Seg[OT],colors='r',linewidth=2)
    DIV_wall=mc.LineCollection(Wall_Seg[DIV],colors='k',linewidth=2)
    IT_wall=mc.LineCollection(Wall_Seg[IT],colors='g',linewidth=2)
    
    wallfig, wallax = plt.subplots(1,2)
    
    wallax[0].add_collection(MC_wall)
    wallax[0].add_collection(OT_wall)
    wallax[0].add_collection(IT_wall)
    wallax[0].add_collection(DIV_wall)
    wallax[0].plot([0.85,0.95],[0,0],'m',linewidth=4)
    wallax[0].set_xlim(RR.min()-0.05,RR.max()+0.05)
    wallax[0].set_ylim(ZZ.min()-0.05,ZZ.max()+0.05)
    wallax[0].set_xlabel('Radial Coordinate (m)')
    wallax[0].set_ylabel('Vertical Coordinate (m)')
    wallax[0].set_aspect('equal')
    
    Wall_Idx=np.linspace(1,NLIM,NLIM)
    RF=(EirOut['fort.44']['wldra(0)'])
    EF=(EirOut['fort.44']['wldpa(0)'])
    
    print('Shot {} Attempt {} Recycled Atomic Flux Data:'.format(ATTEMPT,SHOT))
    print('Main Chamber Recycled Flux = {}'.format(np.sum(RF[MC])))
    print('Outer Target Recycled Flux = {}'.format(np.sum(RF[OT])))
    print('Inner Target Recycled Flux = {}'.format(np.sum(RF[IT])))
    print('Divertor Recycled Flux = {}'.format(np.sum(RF[DIV])))
    print('Total Wall Recycled Flux = {}'.format(np.sum(RF[0:NLIM])))
    
    wallax[1].semilogy(Wall_Idx[MC],RF[MC],'b-')
    #wallax[1].semilogy(Wall_Idx[MC],EirOut['fort.44']['wldra(0)'][MC],'b:')
    
    wallax[1].semilogy(Wall_Idx[OT],RF[OT],'r-')
    #wallax[1].semilogy(Wall_Idx[OT],EirOut['fort.44']['wldra(0)'][OT],'r:')
    
    wallax[1].semilogy(Wall_Idx[DIV],RF[DIV],'k-')
    #wallax[1].semilogy(Wall_Idx[DIV],EirOut['fort.44']['wldra(0)'][DIV],'k:')
    
    wallax[1].semilogy(Wall_Idx[IT],RF[IT],'g-')
    #wallax[1].semilogy(Wall_Idx[IT],EirOut['fort.44']['wldra(0)'][IT],'g:')
    wallax[1].axvline(43,color='m',linewidth=4)
    
    wallax[1].set_xlabel('Wall Segment Index')
    wallax[1].set_ylabel('Particles/Second')
    wallax[1].legend(['Particle Flux reflected from MC Walls',
                      'Particle Flux reflected from Outer Target',
                      'Particle Flux reflected from Divertor Walls',
                      'Particle Flux reflected from Inner Target',
                      'Outer Midplane'])
'''    
''' 
    






