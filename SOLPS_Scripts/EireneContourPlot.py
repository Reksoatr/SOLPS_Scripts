# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:34:21 2020

@author: Richard
"""

from TOOLS import SET_WDIR
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.path import Path

Shot='d3d'
Device='d3d'
Attempt='70'
MeshID='001'

Slice=1
P0=np.array([1.65,-0.65])
P1=np.array([2,-1.1])
Thresh=0.1

BASEDRT, TOPDRT = SET_WDIR('solps-iter/runs/','')
if 'd3d' in Shot:
    BASEDRT = '{}{}'.format(BASEDRT,Shot)
    DRT = '{}/Attempt{}/EirOutput'.format(BASEDRT, str(Attempt))
else:
    BASEDRT = '{}cmod/0{}home'.format(BASEDRT, Shot[-2:])
    DRT = '{}/Attempt{}/EirOutput'.format(BASEDRT, str(Attempt))

tz=np.loadtxt('{}/TriangVertLoc{}'.format(DRT,Attempt),usecols = (2))

tr=np.loadtxt('{}/TriangRadLoc{}'.format(DRT,Attempt),usecols = (2))

VVFILE = np.loadtxt('{}/vvfile.ogr'.format(BASEDRT))

NeuDen=np.loadtxt('{}/EirAtom{}'.format(DRT,Attempt),usecols = (2))
MolDen=np.loadtxt('{}/EirMol{}'.format(DRT,Attempt),usecols = (2))
NeuEng=np.loadtxt('{}/AtomEnergy{}'.format(DRT,Attempt),usecols = (2))
MolEng=np.loadtxt('{}/MolEnergy{}'.format(DRT,Attempt),usecols = (2))

LogVal=np.ma.log10(NeuEng)

Nodes=np.fromfile('{}/{}.tria.{}.nodes'.format(BASEDRT,Device,MeshID),sep=' ')
NN=int(Nodes[0])
XNodes=Nodes[1:NN+1]/100
YNodes=Nodes[NN+1:]/100

Triangles=np.loadtxt('{}/{}.tria.{}.cells'.format(BASEDRT,Device,MeshID),skiprows=1, usecols=(1,2,3))

TP=tri.Triangulation(XNodes,YNodes,triangles=(Triangles-1))

EirFig, EirAx = plt.subplots()
 
'''
Triangles=np.loadtxt('{}/triangles.dat'.format(BASEDRT))
Triangles=Triangles/100
NT=int(Triangles.shape[0]/4)
for n in range(0,NT):
    EirAx.plot(Triangles[n*4:(n*4)+3,0],Triangles[n*4:(n*4)+3,1],'k-')
'''

#IM=EirAx.tripcolor(TP,LogVal)

#EirAx.plot(tr,tz,'r.')

EirAx.plot(VVFILE[:,0]/1000,VVFILE[:,1]/1000,'k-')

#plt.colorbar(IM,ax=EirAx)

if Slice == 1:
    PP=P1-P0
    Theta=np.arctan(PP[1]/PP[0])
    displace=Thresh*np.array([-np.sin(Theta),np.cos(Theta)])
    P0A=P0+displace
    P0B=P0-displace
    P1A=P1+displace
    P1B=P1-displace
    EirAx.plot([P0[0],P1[0]],[P0[1],P1[1]],'b-')
    EirAx.plot([P0A[0],P1A[0]],[P0A[1],P1A[1]],'g-')
    EirAx.plot([P0B[0],P1B[0]],[P0B[1],P1B[1]],'r-')
    Chord=Path([P0A,P1A,P1B,P0B])
    Band=Chord.contains_points(np.vstack((tr,tz)).T)
    IM=EirAx.tripcolor(TP,Band)
    EirAx.set_aspect('equal')
    EirAx.grid()
    
    