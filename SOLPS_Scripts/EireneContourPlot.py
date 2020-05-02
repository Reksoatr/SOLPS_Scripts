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
Attempt='70'

BASEDRT, TOPDRT = SET_WDIR('solps-iter/runs/','')
if 'd3d' in Shot:
    BASEDRT = '{}{}'.format(BASEDRT,Shot)
    DRT = '{}/Attempt{}/EirOutput'.format(BASEDRT, str(Attempt))
else:
    BASEDRT = '{}cmod/0{}home'.format(BASEDRT, Shot[-2:])
    DRT = '{}/Attempt{}/EirOutput'.format(BASEDRT, str(Attempt))

tz=np.loadtxt('{}/TriangVertLoc{}'.format(DRT,Attempt),usecols = (2))

tr=np.loadtxt('{}/TriangRadLoc{}'.format(DRT,Attempt),usecols = (2))

#Triangles=np.loadtxt('{}Sample_IO_Files/template.g011.tria.ogr'.format(TOPDRT))

VVFILE = np.loadtxt('{}/vvfile.ogr'.format(BASEDRT))

NeuDen=np.loadtxt('{}/EirAtom{}'.format(DRT,Attempt),usecols = (2))
MolDen=np.loadtxt('{}/EirMol{}'.format(DRT,Attempt),usecols = (2))
NeuEng=np.loadtxt('{}/AtomEnergy{}'.format(DRT,Attempt),usecols = (2))
MolEng=np.loadtxt('{}/MolEnergy{}'.format(DRT,Attempt),usecols = (2))

LogVal=np.ma.log10(NeuEng)

#TP=tri.Triangulation(Triangles[:,0]/1000,Triangles[:,1]/1000)

EirFig, EirAx = plt.subplots() 

IM=EirAx.tripcolor(tr,tz,LogVal)

#IM=EirAx.tripcolor(TP,LogVal[0:TP.triangles.shape[0]])

#EirAx.plot(Triangles[:,0]/1000,Triangles[:,1]/1000,'k.')

#EirAx.plot(tr,tz,'r.')

EirAx.plot(VVFILE[:,0]/1000,VVFILE[:,1]/1000,'k-')

EirAx.set_aspect('equal')
EirAx.grid()
plt.colorbar(IM,ax=EirAx)