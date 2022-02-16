# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 14:34:44 2022

@author: 18313
"""

import re
import matplotlib.pyplot as plt

def InputfileParser(file='b2.transport.inputfile', plot=False):
    
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
            ax1[ii].set_ylabel(r'{} $[m^2/s]$'.format(Coefficients[jj]))
        
        ax1[0].set_title(file)    
        ax1[-1].set_xlabel(r'$R-R_{sep}$')
        
    return Points

def Generate(trans_pts):
    '''
    Function that is used to turn the radial points into a readable
    b2.transport.inputfile

    Parameters
    ----------
    trans_pts : should be 2d point array, x coordinates being r-r_sep and
    y coordinates the diffusivity at that point

    Returns a data frame for use in the b2.transport.inputfile
    -------
    

    '''
    #J = 1
          #print(self._points)
    n = len(trans_pts)
    m = 0
    i = 1
    j = 1
    r = trans_pts
    print(' ndata(1, {0}, {1})= {2},'.format(i,j,n))
    for m in range(n):
        print(' tdata(1, {0}, {1}, {2})= {3}, tdata(2, {0}, {1}, {2})= {4},'.format(m+1,i,j,round(r[m][0],5),round(r[m][1],5)))


if __name__=='__main__':
    
    Points=InputfileParser(file='b2.transport.inputfile.NewV125',plot=True)