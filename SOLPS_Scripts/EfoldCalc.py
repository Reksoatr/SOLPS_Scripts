# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 18:25:30 2019

@author: rmreksoatmodjo
"""
import os
import numpy as np
import xarray as xr
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import colors, cm
#from VesselPlotter import SOLPSPlotter
from VesselPlotterNew import SOLPSPLOT

plt.rc('font',size=25)
plt.rc('lines',linewidth=5,markersize=15)

Shot = 'd3d'
Attempt = 86
Jxi = 40
Jxa = 56
sep = 21
ND0 = 9 #Attempt 129 -> 14; Attempt 58 -> 9 
NDF = 26 #Attempt 129 -> 30; Attempt 58 -> 33
IF0 = 19
IFF = 35

#IONFLX = SOLPSPlotter(Shot,[Attempt],'IonFlx','Export',RRad='radial')
#NEUDEN = SOLPSPlotter(Shot,[Attempt],'NeuDen','Export',RRad='radial')

SOLPSOBJ = SOLPSPLOT(Shot,[Attempt],Parameter=['IonFlx','NeuDen'])

Rsep = SOLPSOBJ.RadCoords['RadLoc'].loc[sep,:,Attempt]
Vsep = SOLPSOBJ.RadCoords['VertLoc'].loc[sep,:,Attempt]

RRsep = xr.DataArray(np.zeros(SOLPSOBJ.RadCoords['RadLoc'].shape), coords=[SOLPSOBJ.RadCoords['RadLoc'].coords['Radial_Location'].values,SOLPSOBJ.RadCoords['RadLoc'].coords['Poloidal_Location'].values,[Attempt]], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'R-Rsep $m$')
for i in (SOLPSOBJ.RadCoords['RadLoc'].coords['Radial_Location'].values):
    for j in (SOLPSOBJ.RadCoords['RadLoc'].coords['Poloidal_Location'].values):
        RRsep.loc[i,j,Attempt] = np.sign(i-sep)*np.sqrt((SOLPSOBJ.RadCoords['RadLoc'].loc[i,j,Attempt]-Rsep.loc[j])**2 + (SOLPSOBJ.RadCoords['VertLoc'].loc[i,j,Attempt]-Vsep.loc[j])**2)

NeuDen = SOLPSOBJ.PARAM['NeuDen']
NeuDen.values[NeuDen.values==0]=np.nan

IonFlx = SOLPSOBJ.PARAM['IonFlx']
IonFlx.values[IonFlx.values==0]=np.nan
IonFlxPlus = IonFlx.values[IonFlx.values>0]
IonFlxMinus = np.abs(IonFlx.values[IonFlx.values<0])

NDFit = np.ones((SOLPSOBJ.RadCoords['RadLoc'].coords['Poloidal_Location'].size,2))
IFFit = np.ones((SOLPSOBJ.RadCoords['RadLoc'].coords['Poloidal_Location'].size,2))

print('Inner Midplane at Poloidal Grid Cell ' + str(Jxi))
print('Outer Midplane at Poloidal Grid Cell ' + str(Jxa))
print('')

for jxa in range(24,72): #24 to 72 covers entire core region
    NDFit[jxa-24,:] = np.polyfit(RRsep.loc[ND0:NDF,jxa,Attempt],np.log(NeuDen.loc[ND0:NDF,jxa,Attempt]),1)
    
    IFFit[jxa-24,:] = np.polyfit(RRsep.loc[IF0:IFF,jxa,Attempt],np.log(IonFlx.loc[IF0:IFF,jxa,Attempt]),1)
    
    if IFFit[jxa-24,0]==0:
        IFFit[jxa-24,0]=np.nan
    if NDFit[jxa-24,0]==0:
        NDFit[jxa-24,0]=np.nan
    
    if jxa == Jxi:
        print('INNER MIDPLANE')
        print('Poloidal Grid Cell ' + str(jxa))
        print('NeuDen = ' + str(np.exp(NDFit[jxa-24,1])) + '* e ^ [' + str(NDFit[jxa-24,0]) + ' * X]')
        print('e-folding length = ' + str(1000/NDFit[jxa-24,0]) + ' mm')
        print('IonFlx = ' + str(np.exp(IFFit[jxa-24,1])) + '* e ^ [' + str(IFFit[jxa-24,0]) + ' * X]')
        print('e-folding length = ' + str(1000/IFFit[jxa-24,0]) + ' mm')
        print('')
    if jxa == Jxa:
        print('OUTER MIDPLANE')
        print('Poloidal Grid Cell ' + str(jxa))
        print('NeuDen = ' + str(np.exp(NDFit[jxa-24,1])) + '* e ^ [' + str(NDFit[jxa-24,0]) + ' * X]')
        print('e-folding length = ' + str(1000/NDFit[jxa-24,0]) + ' mm')
        print('IonFlx = ' + str(np.exp(IFFit[jxa-24,1])) + '* e ^ [' + str(IFFit[jxa-24,0]) + ' * X]')
        print('e-folding length = ' + str(1000/IFFit[jxa-24,0]) + ' mm')
        print('')


fig0 = plt.figure(figsize=(14,10))
plt.plot(RRsep.coords['Poloidal_Location'].values,1000/NDFit[:,0],'b+-') #,RRsep.coords['Poloidal_Location'].values,1000/IFFit[:,0],'g+-')
Mmin = np.nanmin(1000/NDFit[:,0]) #IFFit[:,0])
Mmax = np.nanmax(1000/NDFit[:,0])
plt.plot([Jxi, Jxi],[Mmin, Mmax],color='Red', linestyle='solid', linewidth=3)
plt.plot([Jxa, Jxa],[Mmin, Mmax],color='Orange', linestyle='solid', linewidth=3)
plt.title('Poloidal profile of Radial e-folding lengths')
plt.xlabel('Poloidal Grid Cell Index X')
plt.ylabel('e-folding length (mm)')
plt.legend(['Neutral Density e-fold length','Inner Midplane','Outer Midplane'],loc=2) #,'Ionization e-fold length'
a = plt.gca()
#a.set_xticklabels(['%.f' % i for i in a.get_xticks()], fontsize='x-large')
#a.set_yticklabels(['%.f' % j for j in a.get_yticks()], fontsize='x-large')
plt.grid()

jxa = Jxa
NeuDenFit1 = np.exp(NDFit[jxa-24,1]) * np.exp(NDFit[jxa-24,0]*RRsep.loc[ND0:NDF,jxa,Attempt])

fig1 = plt.figure(figsize=(14,10))
plt.plot(RRsep.loc[:,jxa,Attempt],NeuDen.loc[:,jxa,Attempt],'x',RRsep.loc[ND0:NDF,jxa,Attempt],NeuDenFit1.values,'-')
NDmin = float(NeuDenFit1.min())
NDmax = float(NeuDenFit1.max())
plt.plot([RRsep.loc[sep,jxa,Attempt], RRsep.loc[sep,jxa,Attempt]],[NDmin, NDmax],color='Black',linewidth=3)
plt.title('Neutral Density Profile at Outer Midplane')
plt.xlabel('R-Rsep (m)')
plt.ylabel(r'Neutral Density ($m^{-3}$)')
plt.legend(['SOLPS Data', 'Exp. Fit'],loc=2)
a = plt.gca()
a.ticklabel_format(axis='y',style='scientific')
#a.set_xticklabels(['%.2f' % i for i in a.get_xticks()], fontsize='x-large')
#a.set_yticklabels(['%.1e' % j for j in a.get_yticks()], fontsize='x-large')
plt.grid()

'''
IonFlxFit1 = np.exp(IFFit[jxa-24,1]) * np.exp(IFFit[jxa-24,0]*RRsep.loc[IF0:IFF,jxa,Attempt])

fig2 = plt.figure(figsize=(14,10))
plt.plot(RRsep.loc[:,jxa,Attempt],IonFlx.loc[:,jxa,Attempt],'x',RRsep.loc[IF0:IFF,jxa,Attempt],IonFlxFit1.values,'-')
IFmin = float(IonFlxFit1.min())
IFmax = float(IonFlxFit1.max())
plt.plot([RRsep.loc[sep,jxa,Attempt], RRsep.loc[sep,jxa,Attempt]],[IFmin, IFmax],color='Black')
plt.title('Discharge 0' + str(Shot) + ' Attempt ' + str(Attempt) + ' Radial Particle Flux/Ionization Profile at Midplane')
plt.xlabel('R-Rsep (m)')
plt.ylabel('Radial Particle Flux/Ionization')
plt.legend(['SOLPS Data', 'Exp. Fit'])
plt.grid()
'''

