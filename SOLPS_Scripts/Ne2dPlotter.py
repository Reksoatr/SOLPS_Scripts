# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 12:27:55 2019

@author: rmreksoatmodjo
"""
import numpy as np
import matplotlib.pyplot as plt

Contour = 0

ATT=[1000,1001,1002,1003,1004]

JXA = int(54/2)
SEP = 18
SOL = 27
PED = 12

NeesepmTimefig, AxNesepm = plt.subplots(nrows=1, ncols=1)
NeSOLmTimefig, AxNeSOLm = plt.subplots(nrows=1, ncols=1)
NePedmTimefig, AxNePedm = plt.subplots(nrows=1, ncols=1)

for Attempt in ATT:

    dirT = 'Attempt{}/Output/'.format(Attempt)
    
    TimeS = np.loadtxt('{}TimeStamps'.format(dirT),unpack=1)
    
    NT = 1000
    
    NeTime = np.loadtxt('{}TimeTraces'.format(dirT),usecols=0)
    
    if len(TimeS) != NT:
        NTi = TimeS[0]
        TimeS = np.zeros(NT)
        TimeS[0]=NTi
        for i in range(NT-1):
            TimeS[i+1]=TimeS[i] + 1e-4
    
    ne2d=np.loadtxt('{}NE2D{}'.format(dirT, Attempt)).reshape((NT,38,49))
    
    if Contour == 1: 
        XX, YY = np.meshgrid(range(49),range(38))
    
        fig1 = plt.figure(figsize=(10,10))
        plt.contourf(XX,YY,ne2d[0,:,:])
    
    AxNesepm.plot(TimeS,ne2d[:,SEP,JXA])
    #AxNesepm.plot(TimeS,NeTime)
    #AxNesepm.legend(('ne2d_{}Y_{}X_Attempt{}'.format(SEP,JXA*2,Attempt),'nesepm_{}'.format(Attempt)))
    
    AxNeSOLm.plot(TimeS,ne2d[:,SOL,JXA])
    
    AxNePedm.plot(TimeS,ne2d[:,PED,JXA])

    
AxNesepm.legend(('Attempt1000','Attempt1001','Attempt1002','Attempt1003','Attempt1004'))
AxNesepm.set_title('Time trace of Ne at Separatrix (Psin=1), Outer Midplane')
AxNesepm.set_xlabel('Time (s)')
AxNesepm.set_ylabel('Electron Density ne (m^-3)')

AxNeSOLm.legend(('Attempt1000','Attempt1001','Attempt1002','Attempt1003','Attempt1004'))
AxNeSOLm.set_title('Time trace of Ne in SOL (Psin = 1.06), Outer Midplane')
AxNeSOLm.set_xlabel('Time (s)')
AxNeSOLm.set_ylabel('Electron Density ne (m^-3)')

AxNePedm.legend(('Attempt1000','Attempt1001','Attempt1002','Attempt1003','Attempt1004'))
AxNePedm.set_title('Time trace of Ne at Pedestal Top (Psin=0.95), Outer Midplane')
AxNePedm.set_xlabel('Time (s)')
AxNePedm.set_ylabel('Electron Density ne (m^-3)')




