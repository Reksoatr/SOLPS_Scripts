# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 14:33:38 2019

@author: rmreksoatmodjo

Saskia's plots, JRT 4
"""
from VesselPlotterNew import SOLPSPLOT
import matplotlib.pyplot as plt

Base=0
Balloon=1
Gas=0

# BASELINE PLOTS

if Base == 1:

    Base012=SOLPSPLOT('12',[65])
    Base025=SOLPSPLOT('25',[153])
    Based3d=SOLPSPLOT('d3d',[86])
    
    Nefig, NeRadPlot = plt.subplots(nrows=1, ncols=1)
    
    Base012.RadProf('Ne',AX=NeRadPlot,PsinOffset=-0.01,SEP=1)
    Base025.RadProf('Ne',AX=NeRadPlot,PsinOffset=-0.015,SEP=1)
    Based3d.RadProf('Ne',AX=NeRadPlot,Publish=['1160718012','1160718025','175060'])
    
    Tefig, TeRadPlot = plt.subplots(nrows=1, ncols=1)
    
    Base012.RadProf('Te',AX=TeRadPlot,PsinOffset=-0.01,SEP=1)
    Base025.RadProf('Te',AX=TeRadPlot,PsinOffset=-0.015,SEP=1)
    Based3d.RadProf('Te',AX=TeRadPlot,Publish=['1160718012','1160718025','175060'])
    
    IonFlxfig, IonFlxRadPlot = plt.subplots(nrows=1, ncols=1)
    
    Base012.RadProf('IonFlx',AX=IonFlxRadPlot,SEP=1)
    Base025.RadProf('IonFlx',AX=IonFlxRadPlot,SEP=1)
    Based3d.RadProf('IonFlx',AX=IonFlxRadPlot,Publish=['1160718012','1160718025','175060'])
    
    DNfig, DNRadPlot = plt.subplots(nrows=1, ncols=1)
    
    Base012.RadProf('DN',AX=DNRadPlot,SEP=1)
    Base025.RadProf('DN',AX=DNRadPlot,SEP=1)
    Based3d.RadProf('DN',AX=DNRadPlot,Publish=['1160718012','1160718025','175060'])
    
    NeuDenfig, NeuDenRadPlot = plt.subplots(nrows=1, ncols=1)
    
    Base012.RadProf('NeuDen',AX=NeuDenRadPlot,LOG10=2,SEP=1)
    Base025.RadProf('NeuDen',AX=NeuDenRadPlot,LOG10=2,SEP=1)
    Based3d.RadProf('NeuDen',AX=NeuDenRadPlot,LOG10=2,Publish=['1160718012','1160718025','175060'])
    
if Balloon == 1:
    
    Bal012=SOLPSPLOT('12',[58,60,62],Publish=['No Ballooning','Weak Ballooning','Strong Ballooning'])
    Bal025=SOLPSPLOT('25',[129,142,144],Publish=['No Ballooning','Weak Ballooning','Strong Ballooning'])
    Bald3d=SOLPSPLOT('d3d',[72,73,75],Publish=['No Ballooning','Weak Ballooning','Strong Ballooning'])
    
    BalNefig, BalNeRadPlot = plt.subplots(nrows=1, ncols=1)
    
    Bal012.RadProf('Ne',AX=BalNeRadPlot,PsinOffset=-0.01,SEP=1)
    Bal025.RadProf('Ne',AX=BalNeRadPlot,PsinOffset=-0.015,SEP=1)
    Bald3d.RadProf('Ne',AX=BalNeRadPlot)
    
    NeuDenPolfig, NeuDenPolPlot = plt.subplots(nrows=1, ncols=1)
    
    Bal012.PolPlot('NeuDen',AX=NeuDenPolPlot,LOG10=2)
    Bal025.PolPlt('NeuDen',AX=NeuDenPolPlot,LOG10=2)
    Bald3d.PolPlot('NeuDen',AX=NeuDenPolPlot,LOG10=2)
    
    DNPolfig, DNPolPlot = plt.subplots(nrows=1, ncols=1)
    
    Bal012.PolPlot('DN',AX=DNPolPlot)
    Bal025.PolPlot('DN',AX=DNPolPlot)
    Bald3d.PolPlot('DN',AX=DNPolPlot)
