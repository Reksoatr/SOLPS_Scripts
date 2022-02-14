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
from TOOLS import WALL_INTERSECT, OpenRemoteFile, SSH_config, JumpConnect

def SOLPSDiagnostiChorder(filepath, 
                          device='CMOD', plot=True,
                          RemoteSave=False, RemotePath=None,
                          Output=None, EndKey='tang', 
                          Extend2wall=False, Reverse=False):
    
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
            C0['X']=np.array([C0['R'][i]*np.cos(np.radians(C0['phi'][i])) for i in range(n)])
            C0['Y']=np.array([C0['R'][i]*np.sin(np.radians(C0['phi'][i])) for i in range(n)])
            C1['X']=np.array([C1['R'][i]*np.cos(np.radians(C1['phi'][i])) for i in range(n)])
            C1['Y']=np.array([C1['R'][i]*np.sin(np.radians(C1['phi'][i])) for i in range(n)])
            
    elif fmt == '.chords': #.chords format ONLY has X,Y,Z coords
        HH=np.loadtxt(file,skiprows=1)
        HH=HH.T
        n=len(HH[0])
        C0={'X':np.array(HH[0]),'Y':np.array(HH[1]),'Z':np.array(HH[2])}
        C1={'X':np.array(HH[3]),'Y':np.array(HH[4]),'Z':np.array(HH[5])}
    
    if Extend2wall:
        P1,P2=WALL_INTERSECT(C0,C1,r2)
        C1['X']=P1['X']
        C1['Y']=P1['Y']
        
    if Output:
        outpath,savename=os.path.split(Output)
        oname,ofmt=os.path.splitext(savename)
        
        for key in C0.keys():
            if type(C0[key]) is float:
                C0[key]=C0[key]*np.ones(n)
        for key in C1.keys():
            if type(C1[key]) is float:
                C1[key]=C1[key]*np.ones(n)                
        
        ii=[i+1 for i in range(n)]
        if ofmt == '.chords':
            if Reverse:
                MAT=np.array([C0['X'],C0['Y'],C0['Z'],C1['X'],C1['Y'],C1['Z']])
                MAT=np.flip(MAT,axis=1)
                MAT=np.vstack((MAT,ii))
            else:                              
                MAT=np.array([C0['X'],C0['Y'],C0['Z'],C1['X'],C1['Y'],C1['Z'],ii])
            MAT=MAT.T
            np.savetxt(Output,MAT,comments='',
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

RemotePath='/sciclone/data10/rmreksoatmodjo/solps-iter/data/chords/XXX'  
RemotePath='/sciclone/data10/rmreksoatmodjo/gfileProcessing/cmod_files/'    
'''        
        