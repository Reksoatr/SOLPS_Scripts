# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 12:17:16 2019

@author: rmreksoatmodjo
"""
#import xarray as xr
#import numpy as np
import sys
import matplotlib.pyplot as plt
from HPCPlotter import SOLPSPlotter

def SOLPSDataSpread(Shot, Attempt):

	plt.rc('font',size=20)
	plt.rc('lines',linewidth=5,markersize=5)
	plt.rc('figure',titlesize=30)

	#print(Attempt)
	#print(type(Attempt))

	if len(Attempt) > 1:
		R = 0
	else:
		R = 1

	Spreadfig, Spreadax = plt.subplots(4,(2+(2*R)),sharex=True,figsize=(48,27))

	for i, item in enumerate(['Ne', 'NeuDen', 'DN', 'IonFlx', 'Te', 'Ti', 'KYE', 'KYI']):
		
		if item == 'NeuDen':
			LG = 2
		else:
			LG = 0
		if i < 2:
			#Spreadax.flatten()[i].set_xlabel('')
			i = i
		elif i < 4: 
			i = i+(2*R)
		elif i < 6:
			i = i+(4*R)
		else:
			i = i+(6*R)

		SOLPSPlotter(str(Shot), Attempt, item,'RadProf', Publish=Attempt, SEP=21, LOG10=LG, ax=plt.subplot(4,(2+(2*R)),i+1))
			
	if len(Attempt) > 1:
		plt.suptitle('Shot: ' + Shot + '\nAttempts: ' + ', '.join(Attempt))
		plt.tight_layout()
		plt.savefig('plots/Shot' + Shot + '_Attempts' + '_'.join(Attempt) + '.png')
	else:
		SOLPSPlotter(str(Shot), Attempt, 'NeuDen', 'Contour', Publish=Attempt, SEP=21, LOG10=2, ax=plt.subplot(122))
		plt.suptitle('Shot: ' + Shot + '\nAttempt: ' + Attempt[0], x=0.55)
		plt.tight_layout()
		plt.savefig('plots/Shot' + Shot + '_Attempt' + Attempt[0] + '.png')

if __name__ == '__main__':
	SOLPSDataSpread(sys.argv[1],sys.argv[2:])	

