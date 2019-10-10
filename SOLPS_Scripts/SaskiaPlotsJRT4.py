# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 14:33:38 2019

@author: rmreksoatmodjo

Saskia's plots, JRT 4
"""
from VesselPlotterNew import SOLPSPLOT
import matplotlib.pyplot as plt
import numpy as np

Base=1
Balloon=0
Gas=0
GasTime=0

# BASELINE PLOTS

if Base == 1:

    Base012=SOLPSPLOT('12',[65])
    Base025=SOLPSPLOT('25',[153])
    #Based3d=SOLPSPLOT('d3d',[86])
    
    Fig, RadPlot = plt.subplots(nrows=6,ncols=1,sharex=True)
    
    #Nefig, NeRadPlot = plt.subplots(nrows=1, ncols=1)
    
    Base012.RadProf('Ne',AX=RadPlot[0],PsinOffset=-0.01,Markers=False,PlotScheme=['b-'])
    Base025.RadProf('Ne',AX=RadPlot[0],PsinOffset=-0.015,Markers=False,Publish=['High Ip','Low Ip'],PlotScheme=['r-'])
    #Based3d.RadProf('Ne',AX=NeRadPlot,Publish=['1160718012','1160718025','175060'])
    
    #Tefig, TeRadPlot = plt.subplots(nrows=1, ncols=1)
    
    Base012.RadProf('Te',AX=RadPlot[2],PsinOffset=-0.01,Markers=False,PlotScheme=['b-'])
    Base025.RadProf('Te',AX=RadPlot[2],PsinOffset=-0.015,Markers=False,Publish=['High Ip','Low Ip'],PlotScheme=['r-'])
    #Based3d.RadProf('Te',AX=TeRadPlot,Publish=['1160718012','1160718025','175060'])
    
    #IonFlxfig, IonFlxRadPlot = plt.subplots(nrows=1, ncols=1)
    
    Base012.RadProf('IonFlx',AX=RadPlot[4],Markers=False,PlotScheme=['b-'])
    Base025.RadProf('IonFlx',AX=RadPlot[4],Markers=False,Publish=['High Ip','Low Ip'],PlotScheme=['r-'])
    #Based3d.RadProf('IonFlx',AX=IonFlxRadPlot,Publish=['1160718012','1160718025','175060'])
    
    #DNfig, DNRadPlot = plt.subplots(nrows=1, ncols=1)
    
    Base012.RadProf('DN',AX=RadPlot[1],Markers=False,PlotScheme=['b-'])
    Base025.RadProf('DN',AX=RadPlot[1],Markers=False,Publish=['High Ip','Low Ip'],PlotScheme=['r-'])
    #Based3d.RadProf('DN',AX=DNRadPlot,Publish=['1160718012','1160718025','175060'])
    
    Base012.RadProf('KYE',AX=RadPlot[3],Markers=False,PlotScheme=['b-'])
    Base025.RadProf('KYE',AX=RadPlot[3],Markers=False,Publish=['High Ip','Low Ip'],PlotScheme=['r-'])
    
    #NeuDenfig, NeuDenRadPlot = plt.subplots(nrows=1, ncols=1)
    
    Base012.RadProf('NeuDen',AX=RadPlot[5],LOG10=1,Markers=False,PlotScheme=['b-'])
    Base025.RadProf('NeuDen',AX=RadPlot[5],LOG10=1,Markers=False,Publish=['High Ip','Low Ip'],PlotScheme=['r-'])
    #Based3d.RadProf('NeuDen',AX=NeuDenRadPlot,LOG10=2,Publish=['1160718012','1160718025','175060'])
    
if Balloon == 1:
    
    Bal012=SOLPSPLOT('12',[58,60,62],Publish=['No Ballooning','Weak Ballooning','Strong Ballooning'], PlotScheme=['b-','m-','c-'])
    Bal025=SOLPSPLOT('25',[129,142,144],Publish=['No Ballooning','Weak Ballooning','Strong Ballooning'], PlotScheme=['b.-','m.-','c.-'])
    Bald3d=SOLPSPLOT('d3d',[72,73,75],Publish=['No Ballooning','Weak Ballooning','Strong Ballooning'], PlotScheme=['b--','m--','c--'])
    
    BalNefig, BalNeRadPlot = plt.subplots(nrows=1, ncols=1)
    
    Bal012.RadProf('Ne',AX=BalNeRadPlot,PsinOffset=-0.01,)
    Bal025.RadProf('Ne',AX=BalNeRadPlot,PsinOffset=-0.015,Markers=False)
    Bald3d.RadProf('Ne',AX=BalNeRadPlot,Markers=False)
    
    NeuDenPolfig, NeuDenPolPlot = plt.subplots(nrows=1, ncols=1)
    
    Bal012.PolPlot('NeuDen',AX=NeuDenPolPlot,LOG10=2,Markers=False)
    Bal025.PolPlot('NeuDen',AX=NeuDenPolPlot,LOG10=2,Markers=False)
    Bald3d.PolPlot('NeuDen',AX=NeuDenPolPlot,LOG10=2,Markers=False)
    
    DNPolfig, DNPolPlot = plt.subplots(nrows=1, ncols=1)
    
    Bal012.PolPlot('DN',AX=DNPolPlot,Markers=False)
    Bal025.PolPlot('DN',AX=DNPolPlot,Markers=False)
    Bald3d.PolPlot('DN',AX=DNPolPlot)
    
if Gas == 1:
    Gas012=SOLPSPLOT('gas012',[19,18,17,8],Publish=['77.8 TorrL total D2','39.5 TorrL total D2','8.25 TorrL total D2','No D2 puff'],Markers=False,PlotScheme=['y','r','b','g'])
    Gas025=SOLPSPLOT('gas025',[19,16,1],Publish=['72.2 TorrL total D2','6.77 TorrL total D2','No D2 puff'],Markers=False,PlotScheme=['r','b','g'])
    #Gasd3d=SOLPSPLOT('gasd3d',[30,29,28,27,86],Publish=['300 TorrL total D2','200 TorrL total D2','150 TorrL total D2','50 TorrL total D2','No D2 puff'])
    
    Gas012.RadProf('Ne',PsinOffset=-0.01)
    Gas025.RadProf('Ne',PsinOffset=-0.015)
    #Gasd3d.RadProf('Ne')
    
    Gas012.RadProf('IonFlx')
    Gas025.RadProf('IonFlx')
    #Gasd3d.RadProf('IonFlx')
    
    Gas012.RadProf('NeuDen',LOG10=2)
    Gas025.RadProf('NeuDen',LOG10=2)
    #Gasd3d.RadProf('NeuDen',LOG10=2)

    Gas012.PolPlot('NeuDen',LOG10=2)
    Gas025.PolPlot('NeuDen',LOG10=2)
    #Gasd3d.PolPlot('NeuDen',LOG10=2)

if GasTime == 1:
    
    GasDict = {'gasd3d':[30,29,28,27,86]} #{'gas012':[19,18,17,8],'gas025':[19,16,1],'gasd3d':[30,29,28,27,86]}
    Legend = ['300 TorrL total D2','200 TorrL total D2','150 TorrL total D2','50 TorrL total D2','No D2 puff']
    
    BDRT = r'C:/Users/rmreksoatmodjo/Desktop/My Drive/School stuff/College of William and Mary/Research/SOLPS Stuff/SOLPS_2d_prof/gaspuff'

    XLabel = 'Time (s)'
    NeLabel = r'Electron Density $n_e (m^{-3})$'

    NeTimefig, AxNe = plt.subplots(nrows=1, ncols=1)
    NeITimefig, AxNeI = plt.subplots(nrows=1, ncols=1)
    NeATimefig, AxNeA = plt.subplots(nrows=1, ncols=1)
    
    NeTime = {}
    NeITime = {}
    NeATime = {}
    TimeS = {}

    for key in GasDict:
        for i in GasDict[key]:
            if 'd3d' in key:
                dirT = BDRT + '/d3d/Attempt' + str(i) + '/Output/'
            else:
                dirT = BDRT + '/cmod/0' + key[-2:] + 'home/Attempt' + str(i) + '/Output/'
    
            TimeS[i] = np.loadtxt(dirT+'TimeStamps',unpack=1)  
                
            NeTime[i] = np.loadtxt(dirT+'TimeTraces',usecols=0)
            if len(NeTime[i]) != len(TimeS[i]):
                AxNe.plot(TimeS[i],NeTime[i][::2])
            else:
                AxNe.plot(TimeS[i],NeTime[i])
                
            NeITime[i] = np.loadtxt(dirT+'TimeTraces',usecols=6)
            if len(NeITime[i]) != len(TimeS[i]):
                AxNeI.plot(TimeS[i],NeITime[i][::2])
            else:
                AxNeI.plot(TimeS[i],NeITime[i])
                
            NeATime[i] = np.loadtxt(dirT+'TimeTraces',usecols=3)
            if len(NeATime[i]) != len(TimeS[i]):
                AxNeA.plot(TimeS[i],NeATime[i][::2])
            else:
                AxNeA.plot(TimeS[i],NeATime[i])
                
    AxNe.legend(Legend)
    AxNe.set_title('Outer Midplane Separatrix Electron Density Time Trace')
    AxNe.set_xlabel(XLabel)  
    AxNe.set_ylabel(NeLabel)     
    
    AxNeI.legend(Legend)
    AxNeI.set_title('Inboard Divertor Separatrix Electron Density Time Trace')
    AxNeI.set_xlabel(XLabel)
    AxNeI.set_ylabel(NeLabel)  
    
    AxNeA.legend(Legend)
    AxNeA.set_title('Outboard Divertor Separatrix Electron Density Time Trace')
    AxNeA.set_xlabel(XLabel)
    AxNeA.set_ylabel(NeLabel)    