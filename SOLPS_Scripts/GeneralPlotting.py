# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 13:14:05 2020

@author: Rmreksoatmodjo
"""
from VesselPlotterNew import SOLPSPLOT
import numpy as np
import matplotlib.pyplot as plt

PS=['.','.','.','.','.','-']

Control012=SOLPSPLOT('12',['101', '102', '103', '104', '105'],PsinOffset=-0.005,AVG=True,PlotScheme=PS)
Bal012=SOLPSPLOT('12',['120', '121', '122', '123', '124'],PsinOffset=-0.005,AVG=True,PlotScheme=PS)
Control025=SOLPSPLOT('25',['161', '162', '163', '164', '165'],PsinOffset=-0.0125,AVG=True,PlotScheme=PS,JXA=55)
Bal025=SOLPSPLOT('25',['166', '167', '168', '169', '170'],PsinOffset=-0.0125,AVG=True,PlotScheme=PS,JXA=55)


fig1, ax1= plt.subplots()

Control012.RadProf('Ne',Markers=False,AX=ax1,PlotScheme=[],GRID=True)
Control025.RadProf('Ne',Markers=False,AX=ax1,PlotScheme=[],GRID=True)

fig2, ax2= plt.subplots()

Control012.RadProf('Te',Markers=False,AX=ax2,PlotScheme=[],GRID=True)
Control025.RadProf('Te',Markers=False,AX=ax2,PlotScheme=[],GRID=True)

fig3, ax3= plt.subplots()

Control012.RadProf('NeuDen',LOG10=2,Markers=False,AX=ax3,PlotScheme=PS,GRID=True)
Control025.RadProf('NeuDen',LOG10=2,Markers=False,AX=ax3,PlotScheme=PS,GRID=True)


fig4, ax4= plt.subplots()

Control025.RadProf('NeuDen',LOG10=2,Markers=False,AX=ax4,PlotScheme=PS,GRID=True)
Bal025.RadProf('NeuDen',LOG10=2,Markers=False,AX=ax4,PlotScheme=PS,GRID=True)

fig5, ax5= plt.subplots()

Control025.PolPlot('NeuDen',LOG10=2,Markers=False,AX=ax5,PlotScheme=PS,GRID=True)
Bal025.PolPlot('NeuDen',LOG10=2,Markers=False,AX=ax5,PlotScheme=PS,GRID=True)

fig6, ax6= plt.subplots()

Control025.PolPlot('IonPol',LOG10=0,SURF=24,Markers=False,AX=ax6,PlotScheme=PS,GRID=True)
Bal025.PolPlot('IonPol',LOG10=0,SURF=24,Markers=False,AX=ax6,PlotScheme=PS,GRID=True)

fig7, ax7= plt.subplots()

Control025.PolPlot('IonPol',LOG10=0,SURF=36,Markers=False,AX=ax7,PlotScheme=PS,GRID=True)
Bal025.PolPlot('IonPol',LOG10=0,SURF=36,Markers=False,AX=ax7,PlotScheme=PS,GRID=True)


Control025.Contour('NeuDen',LOG10=2)
Bal025.Contour('NeuDen',LOG10=2)

Control025.Contour('IonPol',LOG10=1)
Bal025.Contour('IonPol',LOG10=1)

NEUDENDIFF = Control025.PARAM['NeuDen'].loc[:,:,'AVG']-Bal025.PARAM['NeuDen'].loc[:,:,'AVG']
Control025.PARAM['NeuDen'].loc[:,:,'AVG']=NEUDENDIFF
Control025.Contour('NeuDen',LOG10=1)

IonPolDIFF = Control025.PARAM['IonPol'].loc[:,:,'AVG']-Bal025.PARAM['IonPol'].loc[:,:,'AVG']
Control025.PARAM['IonPol'].loc[:,:,'AVG']=IonPolDIFF
Control025.Contour('IonPol',LOG10=1)

