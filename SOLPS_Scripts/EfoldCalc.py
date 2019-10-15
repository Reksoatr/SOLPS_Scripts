# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 18:25:30 2019

@author: rmreksoatmodjo
"""
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import colors, cm
from VesselPlotterNew import SOLPSPLOT

plt.rc('font',size=25)
plt.rc('lines',linewidth=5,markersize=15)

Shot = 'gas012'
Attempt = 17
JXI = 40 - 1
JXA = 56 - 1
Crn = 48 - 1
Div1 = 24 - 1
Div2 = 72 - 1
sep = 21
ND0 = np.arange(7,14) #Attempt 129 -> 14; Attempt 58 -> 9 
NDF = np.arange(20,34) #Attempt 129 -> 30; Attempt 58 -> 33
IF0 = 19
IFF = 35

PI = JXA

fitplot = 1

if fitplot == 1:

    SOLPSOBJ = SOLPSPLOT(Shot,[Attempt],Parameter=['IonFlx','NeuDen'])
    
    Rsep = SOLPSOBJ.RadCoords['RadLoc'].loc[sep,:,Attempt]
    Vsep = SOLPSOBJ.RadCoords['VertLoc'].loc[sep,:,Attempt]
    
    RRsep = xr.DataArray(np.zeros(SOLPSOBJ.RadCoords['RadLoc'].shape), coords=[SOLPSOBJ.RadCoords['RadLoc'].coords['Radial_Location'].values,SOLPSOBJ.RadCoords['RadLoc'].coords['Poloidal_Location'].values,[Attempt]], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'R-Rsep $m$')
    for i in (SOLPSOBJ.RadCoords['RadLoc'].coords['Radial_Location'].values):
        for j in (SOLPSOBJ.RadCoords['RadLoc'].coords['Poloidal_Location'].values):
            RRsep.loc[i,j,Attempt] = np.sign(i-sep)*np.sqrt((SOLPSOBJ.RadCoords['RadLoc'].loc[i,j,Attempt]-Rsep.loc[j])**2 + (SOLPSOBJ.RadCoords['VertLoc'].loc[i,j,Attempt]-Vsep.loc[j])**2)
    
    NeuDen = SOLPSOBJ.PARAM['NeuDen']
    NeuDen.values[NeuDen.values==0]=np.nan

'''    
    IonFlx = SOLPSOBJ.PARAM['IonFlx']
    IonFlx.values[IonFlx.values==0]=np.nan
    IonFlxPlus = IonFlx.values[IonFlx.values>0]
    IonFlxMinus = np.abs(IonFlx.values[IonFlx.values<0])
'''
    
    SZ = len(SOLPSOBJ.RadCoords['RadLoc'].coords['Poloidal_Location'])
    
    NDFit = np.ones((SZ,2))
#    IFFit = np.ones((SZ,2))
    NDTrial = {}
    NDResiduals = np.ones((len(NDF),len(ND0)))
    
    print('Inner Midplane at Poloidal Grid Cell ' + str(JXI))
    print('Outer Midplane at Poloidal Grid Cell ' + str(JXA))
    print('')
    
    for Xx in np.arange(24,72): #range(SZ): #24 to 72 covers entire core region
        for n, N in enumerate(NDF):
            for m, M in enumerate(ND0):
                NDTrial[n,m] = np.polyfit(RRsep.loc[M:N,Xx,Attempt],np.log(NeuDen.loc[M:N,Xx,Attempt]),1,full=True)
                #print('Residual for exp fit from {} to {} at Xx={}: {}'.format(M, N, Xx, NDTrial[n,m][1][0]))
                NDResiduals[n, m] = NDTrial[n,m][1][0] 
        Key = np.unravel_index(np.argmin(NDResiduals),NDResiduals.shape)
        NDFit[Xx,:] = NDTrial[Key][0]
        print(Xx, NDTrial[Key][0], NDTrial[Key][1], ND0[Key[1]], NDF[Key[0]])
        
        '''
        #IFFit[Xx-24,:] = np.polyfit(RRsep.loc[IF0:IFF,Xx,Attempt],np.log(IonFlx.loc[IF0:IFF,Xx,Attempt]),1,full=True)   
        if IFFit[Xx-24,0]==0:
            IFFit[Xx-24,0]=np.nan
        '''    
            
        if NDFit[Xx,0]==1:
            NDFit[Xx,0]=1000
        
        if Xx == JXI:
            print('INNER MIDPLANE')
            print('Poloidal Grid Cell ' + str(Xx))
            print('NeuDen = ' + str(np.exp(NDFit[Xx,1])) + '* e ^ [' + str(NDFit[Xx,0]) + ' * X]')
            print('e-folding length = ' + str(1000/NDFit[Xx,0]) + ' mm')
            '''
            print('IonFlx = ' + str(np.exp(IFFit[Xx,1])) + '* e ^ [' + str(IFFit[Xx,0]) + ' * X]')
            print('e-folding length = ' + str(1000/IFFit[Xx,0]) + ' mm')
            print('')
            '''
        if Xx == JXA:
            print('OUTER MIDPLANE')
            print('Poloidal Grid Cell ' + str(Xx))
            print('NeuDen = ' + str(np.exp(NDFit[Xx,1])) + '* e ^ [' + str(NDFit[Xx,0]) + ' * X]')
            print('e-folding length = ' + str(1000/NDFit[Xx,0]) + ' mm')
            ND0A = ND0[Key[1]]
            NDFA = NDF[Key[0]]
            '''
            print('IonFlx = ' + str(np.exp(IFFit[Xx,1])) + '* e ^ [' + str(IFFit[Xx,0]) + ' * X]')
            print('e-folding length = ' + str(1000/IFFit[Xx,0]) + ' mm')
            print('')
            '''
            
    fig0 = plt.figure(figsize=(14,10))
    plt.plot(RRsep.coords['Poloidal_Location'].values,1000/NDFit[:,0],'b+-') #,RRsep.coords['Poloidal_Location'].values,1000/IFFit[:,0],'g+-')
    Mmin = np.nanmin(1000/NDFit[:,0]) #IFFit[:,0])
    Mmax = np.nanmax(1000/NDFit[:,0])
    plt.plot([JXI, JXI],[Mmin, Mmax],color='Red', linestyle='solid', linewidth=3)
    plt.plot([JXA, JXA],[Mmin, Mmax],color='Orange', linestyle='solid', linewidth=3)
    plt.plot([Crn, Crn],[Mmin, Mmax],color='Cyan', linestyle='solid', linewidth=3)
    plt.plot([Div1, Div1],[Mmin, Mmax],color='Green', linestyle='solid', linewidth=3)
    plt.plot([Div2, Div2],[Mmin, Mmax],color='Green', linestyle='solid', linewidth=3)
    plt.title('Poloidal profile of Radial e-folding lengths')
    plt.xlabel('Poloidal Grid Cell Index X')
    plt.ylabel('e-folding length (mm)')
    plt.legend(['Neutral Density e-fold length','Inner Midplane','Outer Midplane', 'Crown', 'Core Boundary']) #,'Ionization e-fold length'
    a = plt.gca()
    #a.set_xticklabels(['%.f' % i for i in a.get_xticks()], fontsize='x-large')
    #a.set_yticklabels(['%.f' % j for j in a.get_yticks()], fontsize='x-large')
    plt.grid()

Xx = PI
NeuDenFit1 = np.exp(NDFit[Xx,1]) * np.exp(NDFit[Xx,0]*RRsep.loc[ND0A:NDFA,Xx,Attempt])

fig1 = plt.figure(figsize=(14,10))
plt.plot(RRsep.loc[:,Xx,Attempt],NeuDen.loc[:,Xx,Attempt],'x',RRsep.loc[ND0A:NDFA,Xx,Attempt],NeuDenFit1.values,'-')
NDmin = float(NeuDenFit1.min())
NDmax = float(NeuDenFit1.max())
plt.plot([RRsep.loc[sep,Xx,Attempt], RRsep.loc[sep,Xx,Attempt]],[NDmin, NDmax],color='Black',linewidth=3)
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
IonFlxFit1 = np.exp(IFFit[Xx-24,1]) * np.exp(IFFit[Xx-24,0]*RRsep.loc[IF0:IFF,Xx,Attempt])

fig2 = plt.figure(figsize=(14,10))
plt.plot(RRsep.loc[:,Xx,Attempt],IonFlx.loc[:,Xx,Attempt],'x',RRsep.loc[IF0:IFF,Xx,Attempt],IonFlxFit1.values,'-')
IFmin = float(IonFlxFit1.min())
IFmax = float(IonFlxFit1.max())
plt.plot([RRsep.loc[sep,Xx,Attempt], RRsep.loc[sep,Xx,Attempt]],[IFmin, IFmax],color='Black')
plt.title('Discharge 0' + str(Shot) + ' Attempt ' + str(Attempt) + ' Radial Particle Flux/Ionization Profile at Midplane')
plt.xlabel('R-Rsep (m)')
plt.ylabel('Radial Particle Flux/Ionization')
plt.legend(['SOLPS Data', 'Exp. Fit'])
plt.grid()
'''

