# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 12:17:16 2019

@author: rmreksoatmodjo
"""
#import xarray as xr
import numpy as np
import sys
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from SOLPS_Plotter import SOLPSPLOT
from textwrap import wrap

def SOLPSDataSpread(Shot, Attempt):
	
	Shot = str(Shot)
	TimeRange=[1.1,1.3]

	if '1040122027' in Shot:
		PsinOffset = 0
		RADC = 'radial'
		JXA = 59
		JXI = 35
		SEP = 20
		XDIM = 98
		YDIM = 38
		CoreBound = [24,71]
		ROOTSHOT=''
	elif '1080416025' in Shot:
		PsinOffset = 0.03
		RADC = 'psin'
		JXA = 40
		JXI = 58
		SEP = 20
		XDIM = 98
		YDIM = 38
		CoreBound = [24,71]
		ROOTSHOT=''
		'''
		JXA = 28
		JXI = 39
		SEP = 13
		XDIM = 66
		YDIM = 26
		CoreBound = [16,47]
		ROOTSHOT=''
		'''
	elif '1100305023' in Shot:
		PsinOffset = -0.014
		RADC = 'psin'
		JXA = 55
		JXI = 37
		SEP = 20
		XDIM = 98
		YDIM = 38
		CoreBound = [24,71]
		ROOTSHOT=''
	elif '1100308004' in Shot:
		PsinOffset = 0.02
		RADC = 'psin'
		JXA = 40
		JXI = 59
		SEP = 20
		XDIM = 98
		YDIM = 38
		CoreBound = [24,71]
		ROOTSHOT=''
	elif '1120917011' in Shot:
		PsinOffset = 0
		RADC = 'psin'
		JXA = 56
		JXI = 38
		SEP = 20
		XDIM = 98
		YDIM = 38
		CoreBound = [24,71]
		ROOTSHOT=''
		TimeRange=[0.8,1.2]
	elif '12' in Shot:
		PsinOffset = -0.005
		RADC = 'psin'
		JXA = 55
		JXI = 37
		SEP = 20
		XDIM = 98
		YDIM = 38
		CoreBound = [24,71]
		ROOTSHOT='1160718'
	elif '25' in Shot:
		PsinOffset = -0.01 
		RADC = 'psin'
		JXA = 55
		JXI = 37
		SEP = 20
		XDIM = 98
		YDIM = 38
		CoreBound = [24,71]
		ROOTSHOT='1160718'		
	else:
		RADC = 'psin'
		JXA = 56
		JXI = 40
		SEP = 20
		XDIM = 98
		YDIM = 38
		CoreBound = [24,71]
		ROOTSHOT='175060'
		PsinOffset = 0

	SOLPSOBJ = SOLPSPLOT(Shot, Attempt, PsinOffset=PsinOffset, TimeRange=TimeRange, GRID=True, BASEDRT= 'solps-iter/runs/', TOPDRT='',RADC=RADC,JXA=JXA,JXI=JXI,SEP=SEP,XDIM=XDIM,YDIM=YDIM,CoreBound=CoreBound,ROOTSHOT=ROOTSHOT)

	plt.rc('font',size=20)
	plt.rc('lines',linewidth=5,markersize=5)
	plt.rc('figure',titlesize=30)

	COORDC = [1,2,0,1,2,0,1,2]
	COORDR = [0,0,1,1,1,2,2,2]
	
	N = len(Attempt)

	R=int(N < 2)

	Spreadfig = plt.figure(figsize=(48,27))
	gs = GridSpec(4, (3+(2*R)))
	for i, item in enumerate(['IonFlx', 'NeuDen', 'Ne', 'Te', 'Ti', 'DN', 'KYE', 'RadPinch']):

		if item == 'NeuDen':
			LG = 2
		else:
			LG = 0

		SOLPSOBJ.RadProf(Parameter=item, Publish=Attempt, SEP=21, LOG10=LG, AX=Spreadfig.add_subplot(gs[COORDR[i],COORDC[i]]))
	
	NeTime = {}
	TeTime = {}
	TiTime = {}
	TimeS = {}
	NoteText = ''	

	XLabel = 'Time (s)'
	NeLabel = r'Electron Density $n_e (m^{-3})$'
	TeLabel = r'Electron Temperature $T_e (eV)$'
	TiLabel = r'Ion Temperature $T_i (eV)$'

	AxNe = Spreadfig.add_subplot(gs[3,0], title='Outer Midplane Separatrix Electron Density Time Trace', xlabel=XLabel, ylabel=NeLabel)
	AxTe = Spreadfig.add_subplot(gs[3,1], title='Outer Midplane Seperatrix Electron Temperature Time Trace', xlabel=XLabel, ylabel=TeLabel)
	AxTi = Spreadfig.add_subplot(gs[3,2], title='Outer Midplane Separatrix Ion Temperature Time Trace', xlabel=XLabel, ylabel=TiLabel)

	if 'gas' in Shot:
        	BDRT = 'solps-iter/runs/gaspuff/'
	else:
        	BDRT = 'solps-iter/runs/'

	for i in Attempt:

		if 'd3d' in Shot:
			dirT = '{}d3d/Attempt{}/Output/'.format(BDRT,i)
		elif ROOTSHOT == '':
			dirT = '{}cmod/{}home/Attempt{}/Output/'.format(BDRT,Shot,i)   
		else:
			dirT = '{}cmod/0{}home/Attempt{}/Output/'.format(BDRT,Shot[-2:],i)

		TimeS[i] = np.loadtxt(dirT+'TimeStamps',unpack=1)  
    
		NeTime[i] = np.loadtxt(dirT+'TimeTraces',usecols=0)
		if len(NeTime[i]) != len(TimeS[i]):
			AxNe.plot(TimeS[i],NeTime[i][::2])
		else:
			AxNe.plot(TimeS[i],NeTime[i])
    
		TeTime[i] = np.loadtxt(dirT+'TimeTraces',usecols=1)
		if len(TeTime[i]) != len(TimeS[i]):
			AxTe.plot(TimeS[i],TeTime[i][::2])
		else:
			AxTe.plot(TimeS[i],TeTime[i])
    
		TiTime[i] = np.loadtxt(dirT+'TimeTraces',usecols=2)
		if len(TiTime[i]) != len(TimeS[i]):
			AxTi.plot(TimeS[i],TiTime[i][::2])
		else:
			AxTi.plot(TimeS[i],TiTime[i])

		try:
			NoteFile = open(dirT+'Note','r')
		except:
			print('No Note Found')
		else:
			Notes = 'Attempt ' + str(i) + ' Note: ' + NoteFile.read()
			Notes = '\n'.join(wrap(Notes, width=55+45*int(N>1), subsequent_indent='...'))
			NoteText = NoteText + '\n' + Notes + '\n'
			NoteFile.close()
	
	AxNe.legend(Attempt)
	AxTe.legend(Attempt)
	AxTi.legend(Attempt)	
	AxNe.grid(b=1)
	AxTe.grid(b=1)
	AxTi.grid(b=1)

	TitleAx = Spreadfig.add_subplot(gs[0,0])
	TitleAx.axis('off')	

	if N > 1:
		TitleAx.text(0, 1, 'Shot: ' + Shot + '\nAttempts: ' + ', '.join(Attempt) + '\nTotal Iterations: ' + str(2*len(TimeS[i])) + '\nSimulated Time: ' + '%.5f' % TimeS[i][-1] + ' seconds' + '\nTime Stamp: ' + datetime.datetime.today().strftime('%I:%M %p, %m/%d/%Y') + '\n' + NoteText, ha='left', va='top', transform=TitleAx.transAxes, wrap=True)
		plt.tight_layout()
		plt.savefig('plots/Shot' + Shot + '_Attempts' + '_'.join(Attempt) + '.pdf',dpi=200)
	else:
		SOLPSOBJ.Contour(Parameter='NeuDen', Publish=Attempt, SEP=21, LOG10=2, AX=Spreadfig.add_subplot(gs[:,3:]))
		TitleAx.text(0, 1, 'Shot: ' + Shot + '\nAttempt: ' + Attempt[0] + '\nTotal Iterations: ' + str(2*len(TimeS[i])) + '\nSimulated Time: ' + '%.5f' % TimeS[i][-1] + ' seconds' + '\nTime Stamp: ' + datetime.datetime.today().strftime('%I:%M %p, %m/%d/%Y') + '\n' + NoteText, ha='left', va='top', transform=TitleAx.transAxes, wrap=True)
		plt.tight_layout()
		plt.savefig('plots/Shot' + Shot + '_Attempt' + Attempt[0] + '.pdf', dpi=200)

if __name__ == '__main__':
	SOLPSDataSpread(sys.argv[1],sys.argv[2:])	

