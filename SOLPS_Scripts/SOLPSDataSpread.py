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
from HPCPlotter import SOLPSPlotter

def SOLPSDataSpread(Shot, Attempt):

	plt.rc('font',size=20)
	plt.rc('lines',linewidth=5,markersize=5)
	plt.rc('figure',titlesize=30)

	COORDC = [1,2,0,1,2,0,1,2]
	COORDR = [0,0,1,1,1,2,2,2]
	
	N = len(Attempt)

	R=int(N < 2)

	Spreadfig = plt.figure(figsize=(48,27))
	gs = GridSpec(4, (3+(2*R)))
	for i, item in enumerate(['IonFlx', 'NeuDen', 'Ne', 'Te', 'Ti', 'DN', 'KYE', 'KYI']):

		if item == 'NeuDen':
			LG = 2
		else:
			LG = 0

		SOLPSPlotter(str(Shot), Attempt, item,'RadProf', Publish=Attempt, SEP=21, LOG10=LG, ax=Spreadfig.add_subplot(gs[COORDR[i],COORDC[i]]))
	
	NeTime = {}
	TeTime = {}
	TiTime = {}
	TimeS = {}
	NoteText = ''	

	XLabel = 'Time (s)'
	NeLabel = r'Electron Density $n_e (m^{-3})$'
	TeLabel = r'Electron Temperature $T_e (eV)$'
	TiLabel = r'Ion Temperature $T_i (eV)$'

	AxNe = Spreadfig.add_subplot(gs[3,0], title='Electron Density Time Trace', xlabel=XLabel, ylabel=NeLabel)
	AxTe = Spreadfig.add_subplot(gs[3,1], title='Electron Temperature Time Trace', xlabel=XLabel, ylabel=TeLabel)
	AxTi = Spreadfig.add_subplot(gs[3,2], title='Ion Temperature Time Trace', xlabel=XLabel, ylabel=TiLabel)

	for i in Attempt:

		dirT = 'solps-iter/runs/' + Shot + '/Attempt' + str(i) + '/Output/'
		TimeS[i] = np.loadtxt(dirT+'TimeStamps',unpack=1)  
    
		NeTime[i] = np.loadtxt(dirT+'TimeTraces',usecols=0)
		AxNe.plot(TimeS[i],NeTime[i][::2])
    
		TeTime[i] = np.loadtxt(dirT+'TimeTraces',usecols=1)
		AxTe.plot(TimeS[i],TeTime[i][::2])
    
		TiTime[i] = np.loadtxt(dirT+'TimeTraces',usecols=2)
		AxTi.plot(TimeS[i],TeTime[i][::2])
		try:
			NoteFile = open(dirT+'Note','r')
		except:
			print('No Note Found')
		else:
			NoteText = NoteText + '\nAttempt ' + str(i) + ' Note: ' + NoteFile.read()
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
		SOLPSPlotter(str(Shot), Attempt, 'NeuDen', 'Contour', Publish=Attempt, SEP=21, LOG10=2, ax=Spreadfig.add_subplot(gs[:,3:]))
		TitleAx.text(0, 1, 'Shot: ' + Shot + '\nAttempt: ' + Attempt[0] + '\nTotal Iterations: ' + str(2*len(TimeS[i])) + '\nSimulated Time: ' + '%.5f' % TimeS[i][-1] + ' seconds' + '\nTime Stamp: ' + datetime.datetime.today().strftime('%I:%M %p, %m/%d/%Y') + '\n' + NoteText, ha='left', va='top', transform=TitleAx.transAxes, wrap=True)
		plt.tight_layout()
		plt.savefig('plots/Shot' + Shot + '_Attempt' + Attempt[0] + '.pdf', dpi=200)

if __name__ == '__main__':
	SOLPSDataSpread(sys.argv[1],sys.argv[2:])	

