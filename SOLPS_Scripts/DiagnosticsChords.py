# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 13:48:41 2022

@author: 18313

SOLPS Diagnostics Chord Plotter
"""

import numpy as np
import pickle as pkl
import os
import matplotlib.pyplot as plt
import matplotlib.patches as ptchs
from TOOLS import OpenRemoteFile, SSH_config, JumpConnect

def SOLPSDiagnostiChorder(filepath, 
                          device='CMOD', plot=True,
                          RemoteSave=False, RemotePath=None,
                          Output=None, EndKey='tang'):
    
    filename,fmt=os.path.splitext(filepath)
    
    ReadType={'.pkl':'rb','.chords':'r'}
    
    if device == 'CMOD':
        r1=0.46
        r2=0.93
        sep=0.91
    
    if RemotePath:
        RR=SSH_config()
        file=OpenRemoteFile('{}{}'.format(RemotePath,filepath),
                            ReadType[fmt],**RR)
    else:
        file=open(filepath,ReadType[fmt])
    
    if fmt == '.pkl':
        HH=pkl.load(file)           
        C0=HH['start']
        C1=HH[EndKey] #Keyword might be 'end' or 'tang'
        n=len(C0['Z'])
        
        if 'R' and 'phi' in C0.keys(): #Convert R,phi to X,Y coords
            C0['X']=[C0['R'][i]*np.cos(np.radians(C0['phi'][i])) for i in range(n)]
            C0['Y']=[C0['R'][i]*np.sin(np.radians(C0['phi'][i])) for i in range(n)]
            C1['X']=[C1['R'][i]*np.cos(np.radians(C1['phi'][i])) for i in range(n)]
            C1['Y']=[C1['R'][i]*np.sin(np.radians(C1['phi'][i])) for i in range(n)]
            
    elif fmt == '.chords': #.chords format ONLY has X,Y,Z coords
        HH=np.loadtxt(file,skiprows=1)
        HH=HH.T
        n=len(HH[0])
        C0={'X':HH[0],'Y':HH[1],'Z':HH[2]}
        C1={'X':HH[3],'Y':HH[4],'Z':HH[5]}
        
    if Output:
        outpath,savename=os.path.split(Output)
        oname,ofmt=os.path.splitext(savename)
        
        ii=[i+1 for i in range(n)]
        if ofmt == '.chords':
            MAT=np.array([C0['X'],C0['Y'],C0['Z'],C1['X'],C1['Y'],C1['Z'],ii])
            MAT=MAT.T
            np.savetxt(Output,MAT,
                       fmt=['%.3f','%.3f','%.3f','%.3f','%.3f','%.3f','%d'],
                       header="'{}' 100".format(oname))
        elif ofmt == '.pkl':
            MAT={'start':C0,'end':C1}
            pkl.dump(MAT,open(Output,'xb'))
        if RemoteSave:
            RR=SSH_config()
            JJ=JumpConnect(**RR)
            JJ.open_sftp().put(Output,'{}{}'.format(RemotePath,savename))            
    
    if plot==True:
        fig1,ax1=plt.subplots()
        ax1.set_xlim(-(r2+0.1),r2+0.1)
        ax1.set_ylim(-(r2+0.1),r2+0.1)
        ax1.add_patch(ptchs.Circle((0,0),r1,fill=False,edgecolor='k',linewidth=3))
        ax1.add_patch(ptchs.Circle((0,0),r2,fill=False,edgecolor='k',linewidth=3))
        ax1.add_patch(ptchs.Circle((0,0),sep,fill=False,edgecolor='k',linewidth=3,linestyle=':'))  
        
        for i in range(n):
            ax1.plot([C0['X'][i],C1['X'][i]],[C0['Y'][i],C1['Y'][i]])
            
        ax1.set_aspect('equal')
        
    return C0,C1

if __name__=='__main__':
    
    
    A=SOLPSDiagnostiChorder('LYA_MID_WALL_Zero.chords', 
                              device='CMOD', plot=True, RemoteSave=False, 
                              RemotePath=None, EndKey='end',
                              Output=None)
    
'''       
Etendue=np.array([4.8e-9,5.5e-9,5.9e-9,6.3e-9,6.7e-9,6.9e-9,7.3e-9,
                  7.5e-9,7.7e-9,7.8e-9,7.9e-9,8.0e-9,7.9e-9,7.8e-9,
                  7.7e-9,7.5e-9,7.4e-9,7.1e-9,6.7e-9,6.5e-9,6.2e-9,5.8e-9])        
'''        
        