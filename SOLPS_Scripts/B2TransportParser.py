# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 14:34:44 2022

@author: 18313
"""

import re

def B2TransportInputfileParser(file='b2.transport.inputfile'):
    
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
        
    return Points

if __name__=='__main__':
    
    Points=B2TransportInputfileParser(file='b2.transport.inputfile.DekeyserV18')