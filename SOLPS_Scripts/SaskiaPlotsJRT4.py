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

    Base012=SOLPSPLOT('12',['46N'],PlotScheme=['b-'],PsinOffset=-0.005)
    Base025=SOLPSPLOT('25',['21N'],PlotScheme=['r-'],PsinOffset=-0.01,SEP=21)
    #Based3d=SOLPSPLOT('d3d',[86])
    '''
    ContFig, ContPlot = plt.subplots(nrows=1,ncols=2,sharey=True)
    
    Base012.Contour('NeuDen',LOG10=1,JXA=56,Publish=['High Ip'],AX=ContPlot[0])
    ContPlot[0].set_xlabel('')
    ContPlot[0].set_ylabel('')
    ContPlot[0].set_title('High Ip')
    Base025.Contour('NeuDen',LOG10=1,JXA=56,Publish=['Low Ip'],AX=ContPlot[1])
    ContPlot[1].set_xlabel('')
    ContPlot[1].set_ylabel('')
    ContPlot[1].set_title('Low Ip')    
    '''
    Fig, RadPlot = plt.subplots(nrows=5,ncols=1,sharex=True)
    
    #Nefig, NeRadPlot = plt.subplots(nrows=1, ncols=1)
    
    Base012.RadProf('Ne',AX=RadPlot[0],Markers=False)
    Base025.RadProf('Ne',AX=RadPlot[0],Markers=True,Publish=['High Ip (1.3 MA)','Low Ip (1.0 MA)'])
    RadPlot[0].set_xlabel('')
    RadPlot[0].set_ylabel('')
    RadPlot[1].set_xlim(left=0.8,right=1.055)
    #RadPlot[0].get_legend().remove()
    RadPlot[0].set_title(r'(a) Electron Density $n_e\;(m^{-3})$',position=(0,0),ha='left',fontsize=25)
    #Based3d.RadProf('Ne',AX=NeRadPlot,Publish=['1160718012','1160718025','175060'])
    
    #Tefig, TeRadPlot = plt.subplots(nrows=1, ncols=1)
    
    Base012.RadProf('Te',AX=RadPlot[1],Markers=False)
    Base025.RadProf('Te',AX=RadPlot[1],Markers=True,Publish=['High Ip','Low Ip'])
    RadPlot[1].set_xlabel('')
    RadPlot[1].set_ylabel('')
    RadPlot[1].get_legend().remove()
    RadPlot[1].set_title(r'(b) Electron Temperature $T_e\;(eV)$',position=(0,0),ha='left',fontsize=25)
    #Based3d.RadProf('Te',AX=TeRadPlot,Publish=['1160718012','1160718025','175060'])
    
    #IonFlxfig, IonFlxRadPlot = plt.subplots(nrows=1, ncols=1)
    '''
    Base012.RadProf('IonFlx',AX=RadPlot[5],Markers=False)
    Base025.RadProf('IonFlx',AX=RadPlot[5],Markers=False,Publish=['High Ip','Low Ip'])
    RadPlot[5].set_xlabel('')
    RadPlot[5].set_ylabel('')
    RadPlot[5].get_legend().remove()
    RadPlot[5].set_title( r'(b) Radial Particle Flux $s^{-1}$',position=(0,0),ha='left')
    #Based3d.RadProf('IonFlx',AX=IonFlxRadPlot,Publish=['1160718012','1160718025','175060'])
    '''
    #DNfig, DNRadPlot = plt.subplots(nrows=1, ncols=1)
    
    Base012.RadProf('DN',AX=RadPlot[2],Markers=False)
    Base025.RadProf('DN',AX=RadPlot[2],Markers=True,Publish=['High Ip','Low Ip'])
    RadPlot[2].set_xlabel('')
    RadPlot[2].set_ylabel('')
    RadPlot[2].get_legend().remove()
    RadPlot[2].set_title(r'(c) Particle Density Diffusivity $D\;(m^2/s)$',position=(0,0),ha='left',fontsize=25)
    #Based3d.RadProf('DN',AX=DNRadPlot,Publish=['1160718012','1160718025','175060'])
    
    Base012.RadProf('KYE',AX=RadPlot[3],Markers=False)
    Base025.RadProf('KYE',AX=RadPlot[3],Markers=True,Publish=['High Ip','Low Ip'])
    RadPlot[3].set_xlabel('')
    RadPlot[3].set_ylabel('')
    RadPlot[3].get_legend().remove()
    RadPlot[3].set_title(r'(d) Electron Thermal Diffusivity $\chi_e (m^2/s)$',position=(0,0),ha='left',fontsize=25)
    
    #NeuDenfig, NeuDenRadPlot = plt.subplots(nrows=1, ncols=1)
    
    Base012.RadProf('NeuDen',AX=RadPlot[4],LOG10=2,Markers=False)
    Base025.RadProf('NeuDen',AX=RadPlot[4],LOG10=2,Markers=True,Publish=['High Ip','Low Ip'])
    RadPlot[4].set_ylabel('')
    RadPlot[4].get_legend().remove()
    RadPlot[4].set_title(r'(e) Neutral Atom (D) Density $(m^{-3})$',position=(0,0.85),ha='left',fontsize=25)
    
    RadPlot[0].legend(['High Ip (1.3 MA)','Low Ip (1.0 MA)'],fontsize=20,loc=1)
    
    plt.subplots_adjust(wspace=0, hspace=0)
    
    #Based3d.RadProf('NeuDen',AX=NeuDenRadPlot,LOG10=2,Publish=['1160718012','1160718025','175060'])
    
if Balloon == 1:
    
    Bal012=SOLPSPLOT('12',[65,60,62],Publish=['No Ballooning','Weak Ballooning','Strong Ballooning'], PlotScheme=['r-','r--','r:'],PsinOffset=-0.005,Markers=False)
    Bal025=SOLPSPLOT('25',[153,142,144],Publish=['No Ballooning','Weak Ballooning','Strong Ballooning'], PlotScheme=['b-','b--','b:'], TimeRange=[0.9,1.0],PsinOffset=-0.015,Markers=False)
    #Bald3d=SOLPSPLOT('d3d',[72,73,75],Publish=['No Ballooning','Weak Ballooning','Strong Ballooning'], PlotScheme=['b--','m--','c--'])
    
    BalNefig, BalNeRadPlot = plt.subplots(nrows=1, ncols=1)
    
    Bal012.RadProf('Ne',AX=BalNeRadPlot)
    Bal025.RadProf('Ne',AX=BalNeRadPlot)
    BalNeRadPlot.set_ylabel('')
    BalNeRadPlot.set_title(r'Outer Midplane Electron Density $n_e\;(m^{-3})$')
    #Bald3d.RadProf('Ne',AX=BalNeRadPlot,Markers=False)
    
    BalTefig, BalTeRadPlot = plt.subplots(nrows=1, ncols=1)
    
    Bal012.RadProf('Te',AX=BalTeRadPlot)
    Bal025.RadProf('Te',AX=BalTeRadPlot)
    BalTeRadPlot.set_ylabel('')
    BalTeRadPlot.set_title(r'Outer Midplane Electron Temperature $T_e\;(eV)$')    
    
    NeuDenPolfig, NeuDenPolPlot = plt.subplots(nrows=1, ncols=1)
    
    #Bal012.PolPlot('NeuDen',AX=NeuDenPolPlot,LOG10=0)
    Bal025.PolPlot('NeuDen',AX=NeuDenPolPlot,LOG10=0)
    NeuDenPolPlot.set_ylabel('')
    #Bald3d.PolPlot('NeuDen',AX=NeuDenPolPlot,LOG10=2,Markers=False)
    
    DNPolfig, DNPolPlot = plt.subplots(nrows=1, ncols=1)
    
    Bal012.PolPlot('DN',AX=DNPolPlot)
    Bal025.PolPlot('DN',AX=DNPolPlot)
    DNPolPlot.set_ylabel('')
    #Bald3d.PolPlot('DN',AX=DNPolPlot)
    
if Gas == 1:
    Gas012=SOLPSPLOT('gas012',[19,18,17,8],Publish=['77.8 TorrL total D2','39.5 TorrL total D2','8.25 TorrL total D2','No D2 puff'],Markers=False,PlotScheme=['k','y','r','r:'])
    Gas025=SOLPSPLOT('gas025',[19,16,1],Publish=['72.2 TorrL total D2','6.77 TorrL total D2','No D2 puff'],Markers=False,PlotScheme=['c','b','b:'])
    #Gasd3d=SOLPSPLOT('gasd3d',[30,29,28,27,86],Publish=['300 TorrL total D2','200 TorrL total D2','150 TorrL total D2','50 TorrL total D2','No D2 puff'])
    
    #Gas012.RadProf('Ne',PsinOffset=-0.01)
    #Gas025.RadProf('Ne',PsinOffset=-0.015)
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