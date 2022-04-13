# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 23:39:31 2018

@author: Richard
"""

import os
import numpy as np
import xarray as xr
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.widgets import Slider
from TOOLS import SET_WDIR

plt.rc('font',size=9)
plt.rc('lines',linewidth=2,markersize=10)

# Will eventually implement as a function. Input:

Shot = '25'           #Shot Number (12, 25, or d3d)
Attempt = '21N'
Geo = 1
Contour = 0
Device='cmod'

CmodReScale = [1.0,0.92,0.835,0.77]#[1.0,0.9,0.835,0.77]# 
D3DReScale = [10, 9, 8.35, 7.7]

jxi = 37    # Poloidal position of Inner Midplane
jxa = 55    # Poloidal position of Outer Midplane
sep = 18    # Radial position of Separatrix

BalRescale = CmodReScale
ballooning = [0.0,1.0,1.5,2.0]#
#Function would begin here

BASEDRT,TOPDRT=SET_WDIR('','')

dir = '{}{}/0{}home'.format(BASEDRT,Device,Shot)    #Generate path

X = np.linspace(24,71,48)
Y = np.linspace(0,35,36)

BB = xr.DataArray(np.loadtxt('{}/bb0{}'.format(dir,Shot),usecols = (3)).reshape((38,98))[1:37,25:73], coords=[Y,X], dims=['Radial_Location','Poloidal_Location'])     #Load Neutral Density data 

bb_ref = np.mean(BB.values)

BB_norm = bb_ref/BB

ResFactor = np.zeros((len(Y),len(X),len(ballooning)))
ResArea = np.zeros((len(Y),len(X),len(ballooning)))
ResDN = np.zeros((len(Y),len(X),len(ballooning)))
ResKYE = np.zeros((len(Y),len(X),len(ballooning)))

for i in range(len(ballooning)):
    ResFactor[:,:,i] = BalRescale[i]*BB_norm**ballooning[i]

Xx, Yy = np.meshgrid(X,Y)   # Create X and Y mesh grid arrays

RadLoc = xr.DataArray(np.zeros((Y.size,X.size)), coords=[Y,X], dims=['Radial_Location','Poloidal_Location'], name = r'Radial Coordinate $m$')
VertLoc = xr.DataArray(np.zeros((Y.size,X.size)), coords=[Y,X], dims=['Radial_Location','Poloidal_Location'], name = r'Vertical Coordinate $m$')
DN = xr.DataArray(np.zeros((Y.size,X.size)), coords=[Y,X], dims=['Radial_Location','Poloidal_Location'], name = r'Particle Diffusion Coefficient $m^{-2}s^{-1}$')
KYE = xr.DataArray(np.zeros((Y.size,X.size)), coords=[Y,X], dims=['Radial_Location','Poloidal_Location'], name = r'Thermal Diffusion Coefficient $m^{-2}s^{-1}$')
SZ = xr.DataArray(np.zeros((Y.size,X.size)), coords=[Y,X], dims=['Radial_Location','Poloidal_Location'], name = r'Poloidal Cross-Sectional Cell Area $m^2$')
FNAY = xr.DataArray(np.zeros((Y.size,X.size)), coords=[Y,X], dims=['Radial_Location','Poloidal_Location'], name = r'Radial Ion Flux $s^{-1}$')

RadLoc.values[:,:] = np.loadtxt(dir + '/Attempt' + str(Attempt) + '/Output/RadLoc' + str(Attempt),usecols = (3)).reshape((38,98))[1:37,25:73]
VertLoc.values[:,:] = np.loadtxt(dir + '/Attempt' + str(Attempt) + '/Output/VertLoc' + str(Attempt),usecols = (3)).reshape((38,98))[1:37,25:73]
DN.values[:,:] = np.loadtxt(dir + '/Attempt' + str(Attempt) + '/Output/DN' + str(Attempt),usecols = (3)).reshape((76,98))[39:75,25:73]
KYE.values[:,:] = np.loadtxt(dir + '/Attempt' + str(Attempt) + '/Output/KYE' + str(Attempt),usecols = (3)).reshape((38,98))[1:37,25:73]
SZ.values[:,:] = np.loadtxt(dir + '/Attempt' + str(Attempt) + '/Output2/SZ' + str(Attempt),usecols = (3)).reshape((38,98))[1:37,25:73]
FNAY.values[:,:] = np.loadtxt(dir + '/Attempt' + str(Attempt) + '/Output/IonFlx' + str(Attempt),usecols = (3)).reshape((76,98))[39:75,25:73]

for g in range(len(ballooning)):
    ResArea[:,:,g] = ResFactor[:,:,g]/SZ.values

for j in range(len(ballooning)):
    ResDN[:,:,j] = ResFactor[:,:,j]*DN#/SZ.values
    
for f in range(len(ballooning)):
    ResKYE[:,:,f] = ResFactor[:,:,f]*KYE#/SZ.values

if Geo==1:
    Xcord = RadLoc.values[:,:]
    Ycord = VertLoc.values[:,:]
else:
    Xcord = Xx
    Ycord = Yy
    
if Contour == 1:
    fig0 = plt.figure(figsize=(14,10))
    levs = np.linspace(BB.values.min(),BB.values.max(),10)
    plt.contourf(Xcord,Ycord,BB.values[:,:],levs,cmap=cm.viridis)
    plt.plot(Xcord[:,(jxa-24)],Ycord[:,(jxa-24)],color='orange',linewidth=3)
    plt.plot(Xcord[:,(jxi-24)],Ycord[:,(jxi-24)],color='red',linewidth=3)
    plt.plot(Xcord[sep,:],Ycord[sep,:],color='black',linewidth=3)
    plt.title('Toroidal Magnetic B-Field (T)')
    plt.xlabel('Radial Location (m)')
    plt.ylabel('Vertical Location (m)')
    #plt.legend(('Outboard Midplane','Inboard Midplane','Separatrix'))
    plt.colorbar()
    a = plt.gca()
    a.set_aspect(1.0)
    plt.grid()
    
    for k in range(len(ballooning)):
        fig1 = plt.figure(figsize=(14,10))
        levs = np.linspace(ResFactor.min(),ResFactor.max(),10)
        plt.contourf(Xcord,Ycord,ResFactor[:,:,k],levs,cmap=cm.viridis)
        plt.plot(Xcord[:,(jxa-24)],Ycord[:,(jxa-24)],color='orange',linewidth=3)
        plt.plot(Xcord[:,(jxi-24)],Ycord[:,(jxi-24)],color='red',linewidth=3)
        plt.plot(Xcord[sep,:],Ycord[sep,:],color='black',linewidth=3)
        plt.title('Transport Rescale Factor (b=' + str(ballooning[k]) + ', b_r=' + str(BalRescale[k]) + ')')
        plt.xlabel('Radial Location (m)')
        plt.ylabel('Vertical Location (m)')
        #plt.legend(('Outboard Midplane','Inboard Midplane','Separatrix'))
        plt.colorbar()
        a = plt.gca()
        a.set_aspect(1.0)
        plt.grid()
    
    for h in range(len(ballooning)):
        fig2 = plt.figure(figsize=(14,10))
        levs = np.linspace(ResDN.min(),ResDN.max(),10)
        plt.contourf(Xcord,Ycord,ResDN[:,:,h],levs,cmap=cm.viridis)
        plt.plot(Xcord[:,(jxa-24)],Ycord[:,(jxa-24)],color='orange',linewidth=3)
        plt.plot(Xcord[:,(jxi-24)],Ycord[:,(jxi-24)],color='red',linewidth=3)
        plt.plot(Xcord[sep,:],Ycord[sep,:],color='black',linewidth=3)
        plt.title(r'$D$ (b=' + str(ballooning[h]) + ', b_r=' + str(BalRescale[h]) + ')')
        plt.xlabel('Radial Location (m)')
        plt.ylabel('Vertical Location (m)')
        #plt.legend(('Outboard Midplane','Inboard Midplane','Separatrix'))
        plt.colorbar()
        a = plt.gca()
        a.set_aspect(1.0)
        plt.grid()
        
    for m in range(len(ballooning)):
        fig2a = plt.figure(figsize=(14,10))
        levs = np.linspace(ResKYE.min(),ResKYE.max(),10)
        plt.contourf(Xcord,Ycord,ResKYE[:,:,m],levs,cmap=cm.viridis)
        plt.plot(Xcord[:,(jxa-24)],Ycord[:,(jxa-24)],color='orange',linewidth=3)
        plt.plot(Xcord[:,(jxi-24)],Ycord[:,(jxi-24)],color='red',linewidth=3)
        plt.plot(Xcord[sep,:],Ycord[sep,:],color='black',linewidth=3)
        plt.title(r'$\chi_{e,i}$ (b=' + str(ballooning[m]) + ', b_r=' + str(BalRescale[m]) + ')')
        plt.xlabel('Radial Location (m)')
        plt.ylabel('Vertical Location (m)')
        #plt.legend(('Outboard Midplane','Inboard Midplane','Separatrix'))
        plt.colorbar()
        a = plt.gca()
        a.set_aspect(1.0)
        plt.grid()

legendtext = ['Outer Midplane', 'Inner Midplane','Control({},{})'.format(ballooning[0],BalRescale[0]),'Ballooning1({},{})'.format(ballooning[1],BalRescale[1]), 'Ballooning2({},{})'.format(ballooning[2],BalRescale[2]), 'Ballooning3({},{})'.format(ballooning[3],BalRescale[3])]
#for u in range(len(ballooning)):
#    legendtext.append('b=' + str(ballooning[u]) + ', b_r=' + str(BalRescale[u]))

BB_scan2 = np.ones((48,len(ballooning)))

FigBalloon, AxBalloon = plt.subplots(nrows=3,ncols=2)

FigFNAY, AxFNAY = plt.subplots()

#fig3a = plt.figure(figsize=(14,7))
AxBalloon[0,0].axvline(jxa,color='orange',linewidth=2)
AxBalloon[0,0].axvline(jxi,color='red',linewidth=2)
AxBalloon[0,0].plot(Xx[sep,:],ResFactor[sep,:,:],linewidth=2)
AxBalloon[0,0].legend(legendtext)
AxBalloon[0,0].set_title('Transport Coefficient Rescale Factor at Separatrix')
AxBalloon[0,0].set_xlabel('Poloidal Coordinate')
AxBalloon[0,0].set_ylabel(r'$F_{TC}$')
AxBalloon[0,0].set_xlim(Xx[sep,0],Xx[sep,-1])
AxBalloon[0,0].grid()

#fig3b = plt.figure(figsize=(14,10))
AxBalloon[1,0].axvline(jxa,color='orange',linewidth=2)
AxBalloon[1,0].axvline(jxi,color='red',linewidth=2)
AxBalloon[1,0].plot(Xx[sep,:],ResArea[sep,:,:],linewidth=2)
AxBalloon[1,0].legend(legendtext)
AxBalloon[1,0].set_title('Area-Corrected Rescale Factor at Separatrix')
AxBalloon[1,0].set_xlabel('Poloidal Coordinate')
AxBalloon[1,0].set_ylabel(r'$F_{TC}/A\;(m^{-2})$')
#a = AxBalloon[1,0].gca()
AxBalloon[1,0].ticklabel_format(axis='y',style='scientific',scilimits=(0,0))
AxBalloon[1,0].grid()

#fig4 = plt.figure(figsize=(14,10))
AxBalloon[0,1].plot(Y,np.sum(ResDN,axis=1),linewidth=2)
AxBalloon[0,1].legend(legendtext[2:6])
AxBalloon[0,1].set_title('Poloidal Sum of Particle Diffusion Coefficient vs Radial Surface')
AxBalloon[0,1].set_xlabel('Radial Coordinate')
AxBalloon[0,1].set_ylabel(r'Poloidal Sum of Diffusion Coefficient $m^2s^{-1}$')
AxBalloon[0,1].set_xlim(Y[0],Y[-1])
AxBalloon[0,1].grid()

#fig5 = plt.figure(figsize=(14,10))
AxBalloon[1,1].plot(Y,np.sum(ResKYE,axis=1),linewidth=2)
AxBalloon[1,1].legend(legendtext[2:6])
AxBalloon[1,1].set_title('Poloidal Sum of Thermal Diffusion Coefficient vs Radial Surface')
AxBalloon[1,1].set_xlabel('Radial Coordinate')
AxBalloon[1,1].set_ylabel(r'Poloidal Sum of Thermal Coefficient $m^2s^{-1}$')
AxBalloon[1,1].grid()

#fig6 = plt.figure(figsize=(14,10))
AxBalloon[2,1].plot(Y,np.sum(ResArea,axis=1),linewidth=2)
AxBalloon[2,1].legend(legendtext[2:6])
AxBalloon[2,1].plot(sep,np.sum(ResArea[sep,:,0],axis=0),'kX',markersize=15)
AxBalloon[2,1].set_title(r'Poloidal Sum of $F_{TC}/A$ vs Radial Surface')
AxBalloon[2,1].set_xlabel('Radial Coordinate')
AxBalloon[2,1].set_ylabel(r'$\sum_x F_{TC}/A\;(m^{-2})$')
#a = plt.gca()
AxBalloon[2,1].ticklabel_format(axis='y',style='scientific',scilimits=(0,0))
AxBalloon[2,1].grid()

AxBalloon[0,0].relim()
AxBalloon[0,0].autoscale_view(True,True,True)

AxBalloon[0,1].relim()
AxBalloon[0,1].autoscale_view(True,True,True)

FNAY_sum=np.sum(FNAY,axis=1)

AxFNAY.plot(Y,FNAY_sum,linewidth=2)
AxFNAY.axvline(sep,color='black',linewidth=2)
AxFNAY.legend(['Attempt{}'.format(Attempt),'LCFS'])
AxFNAY.set_title('Shot{} Attempt{} Poloidal Sum of Radial Ion Flux vs Radial Surface'.format(Shot,Attempt))
AxFNAY.set_xlabel('Radial Coordinate')
AxFNAY.set_ylabel(r'Poloidal Sum of Radial Ion Flux $s^{-1}$')
AxFNAY.text(0.01,0.05,'Total Ion Flux at LCFS\n={:.3e} particles/sec'.format(FNAY_sum[sep].values),transform=AxFNAY.transAxes,fontsize=13,verticalalignment='bottom')
AxFNAY.grid()

plt.show()

#print(np.sum(ResDN.values,axis=1)[0] - 32.289)
#print(np.sum(ResDN.values,axis=1)[35] - 8.380)
#print("Ballooning Rescale = " + str(BalRescale))
#print("Maximum Rescale Factor = " + str(np.amax(BB_scan2)))