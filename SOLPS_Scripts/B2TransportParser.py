# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 14:34:44 2022

@author: 18313
"""

import re
import matplotlib.pyplot as plt

def B2TransportInputfileParser(file='b2.transport.inputfile', plot=False):
    
    Coefficients = {'1':'Particle density-driven diffusivity',
                       '2': 'Particle pressure-driven diffusivity',
                       '3': 'Ion thermal anomalous diffusivity',
                       '4': 'Electron thermal anomalous diffusivity',
                       '5': 'Poloidal-component of the anomalous ”pinch” velocity',
                       '6': 'Radial-component of the anomalous ”pinch” velocity',
                       '7': 'Anomalous viscosity',
                       '8': 'Anomalous radial electrical conductivity',
                       '9': 'Anomalous radial thermo-electric coefficient'}    
    
    with open(file) as f:
        dataList=f.readlines()
    
    Points={}
    ii=1

    while ii<len(dataList)-2: 
        ndata = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", dataList[ii])
        CoeffID = ndata[1]
        PtNo = int(ndata[3])
        XList = []
        YList = []
        for mm in range(PtNo):
            XList.append(float(re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",dataList[mm+ii+1])[4]))
            YList.append(float(re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",dataList[mm+ii+1])[9]))
        Points[CoeffID] = {'X':XList,'Y':YList}
        ii=ii+PtNo+1
        
    if plot:
        dd=len(Points.keys())
        fig1,ax1=plt.subplots(nrows=dd,ncols=1,sharex=True)
        for ii, jj in enumerate(Points.keys()):
            ax1[ii].plot(Points[jj]['X'],Points[jj]['Y'])
            ax1[ii].set_title(r'Radial profile of {}'.format(Coefficients[jj]),y=0.9)
            ax1[ii].set_ylabel(r'{} $[m^2/s]$'.format(Coefficients[jj]))
            
        ax1[-1].set_xlabel(r'$R-R_{sep}$')
        
    return Points

if __name__=='__main__':
    
    Points=B2TransportInputfileParser(file='b2.transport.inputfile.NewV125',plot=True)