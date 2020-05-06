# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:34:21 2020

@author: Richard
"""

from TOOLS import SET_WDIR
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

Shot='d3d'
Device='d3d'
Attempt='70'
MeshID='001'

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

IM=EirAx.tripcolor(TP,LogVal)

#EirAx.plot(tr,tz,'r.')

EirAx.plot(VVFILE[:,0]/1000,VVFILE[:,1]/1000,'k-')

EirAx.set_aspect('equal')
EirAx.grid()
plt.colorbar(IM,ax=EirAx)