# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 17:34:26 2019

@author: 18313

Jerry's CSV Data Reader!
"""
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.io import loadmat
from VesselPlotterNew import SOLPSPLOT
from TOOLS import SET_WDIR

DEV = 'cmod'
SHOT=12
ExpID = [23,13,12]#[25,24]#
COLORS = ['g','c','b','k']#['r','y','k']#
ATTEMPTS=['43N','41N','39N','37N']#['16N','14N','12N']#
PUBLISH=['77.8 TorrL total D2','39.5 TorrL total D2','8.25 TorrL total D2','Ballooned no D2 puff']#['6.77 TorrL total D2','72.2 TorrL total D2','Ballooned no D2 puff']#
TimeRange=[1.1,1.3]
PsinOffset=-0.01

JJ=0
GAS=1
JND=0

BASEDRT, TOPDRT = SET_WDIR('','')

Wdir = "{}MatData".format(TOPDRT)

if JJ == 1:
    Jexp = np.genfromtxt('{}\exp_data_1040122027.csv'.format(Wdir),dtype=float,delimiter=',',skip_header=3,names=True)

    Jkn1d_plasma = np.genfromtxt('{}\kn1d_out_1040122027_750to830ms_update.csv'.format(Wdir),dtype=float,delimiter=',',skip_header=10,skip_footer=778,names=True)

    Jkn1d_H = np.genfromtxt('{}\kn1d_out_1040122027_750to830ms_update.csv'.format(Wdir),dtype=float,delimiter=',',skip_header=69,skip_footer=56,names=True)

    Jkn1d_H2 = np.genfromtxt('{}\kn1d_out_1040122027_750to830ms_update.csv'.format(Wdir),dtype=float,delimiter=',',skip_header=791,names=True)

    fig1 = plt.figure()
    plt.plot(Jexp['R__R_LCFS_m'],Jexp['N_E_m3'])
    plt.plot(Jkn1d_plasma['R__R_LCFS_m'],Jkn1d_plasma['N_E_m3'],'--')
    plt.legend(('Exp','KN1D'))
    plt.title('Electron Density')
    
    fig2 = plt.figure()
    plt.plot(Jexp['R__R_LCFS_m'],Jexp['T_E_m3'])
    plt.plot(Jkn1d_plasma['R__R_LCFS_m'],Jkn1d_plasma['T_E_eV'],'--')
    plt.legend(('Exp','KN1D'))
    plt.title('Electron Temperature')
    
    fig3 = plt.figure()
    plt.plot(Jexp['R__R_LCFS_m'],Jexp['S_ION_m3s'])
    plt.plot(Jkn1d_H['R__R_LCFS_m'],Jkn1d_H['S_ION_m3s'],'--')
    plt.legend(('Exp','KN1D'))
    plt.title('Ion Source')
    
    fig4 = plt.figure()
    plt.semilogy(Jexp['R__R_LCFS_m'],Jexp['N_D_m3'])
    plt.semilogy(Jkn1d_H['R__R_LCFS_m'],Jkn1d_H['NH_m3'],':')
    plt.semilogy(Jkn1d_H2['R__R_LCFS_m'],Jkn1d_H2['NH2_m3'],'--')
    plt.legend(('Exp_NeuDen','KN1D_H','KN1D_H2'))
    plt.title('Neutral Densities')

if GAS == 1:
    ExpData={}
    PsinAvg={}
    RmidAvg={}
    NemidAvg={}
    ErrNe={}
    TemidAvg={}
    ErrTe={}
        
    for Shot in ExpID:
    #    BASEDRT = '{}cmod/0{}home'.format(BASEDRT, Shot)
        
    #    GFILE = '{}gfileProcessing/cmod_files/g11607180{}.01209_974'.format(TOPDRT, Shot)
    #    GF = eq.equilibrium(gfile=GFILE)
        
        ExpFile = '11607180{}'.format(Shot)
        ExpData[Shot] = loadmat('{}gfileProcessing/cmod_files/{}.mat'.format(TOPDRT, ExpFile))    
        
        Times = ExpData[Shot]['time'].flatten()
        
        ti = 0
        while TimeRange[0] > Times[ti]:
            ti = ti+1           
        tf = 0
        while TimeRange[1] > Times[tf]:
            tf = tf+1
        
        Psin = ExpData[Shot]['psin'][:,ti:tf]
        PsinAvg[Shot] = np.mean(Psin, axis=1)
        
        Rmid = ExpData[Shot]['rmid'][:,ti:tf]
        RmidAvg[Shot] = np.mean(Rmid, axis=1)
        
        Nemid = ExpData[Shot]['ne'][:,ti:tf]
        #Nemid[Nemid == 0] = np.nan
        NemidAvg[Shot] = np.median(Nemid, axis=1)
        ErrNe[Shot] = np.median(ExpData[Shot]['nerr'][:,ti:tf])
        NeThresh = (ErrNe[Shot]*2)/NemidAvg[Shot]
        '''for NT in range(len(NeThresh)):
            if np.abs(NeThresh[NT]) > 0.5:
                NemidAvg[Shot][NT] = np.nan
                ErrNe[Shot][NT] = np.nan
        '''
        Temid = ExpData[Shot]['te'][:,ti:tf]
        #Temid[Temid == 0] = np.nan
        TemidAvg[Shot] = np.median(Temid, axis=1)
        ErrTe[Shot] = np.median(ExpData[Shot]['terr'][:,ti:tf])
        TeThresh = (ErrTe[Shot]*2)/TemidAvg[Shot]
        '''for TT in range(len(TeThresh)):
            if np.abs(TeThresh[TT]) > 0.5:
                TemidAvg[Shot][TT] = np.nan
                ErrTe[Shot][TT] = np.nan
        '''        
        
    Nefig, NeRadPlot = plt.subplots(nrows=1, ncols=1)
    Tefig, TeRadPlot = plt.subplots(nrows=1, ncols=1)
    
    if SHOT == 12:              
        Gas012 = SOLPSPLOT('12',ATTEMPTS,Publish=PUBLISH,PsinOffset=PsinOffset,Markers=False,PlotScheme=COLORS,EXP=False)
        Gas012.RadProf('Ne',AX=NeRadPlot)
        Gas012.RadProf('Te',AX=TeRadPlot)
        
        for n, Shot in enumerate(ExpID):
            NeRadPlot.errorbar(PsinAvg[Shot]+PsinOffset,NemidAvg[Shot],yerr=ErrNe[Shot],fmt='o',color=COLORS[n],markersize=7,linewidth=3,capsize=7)
            TeRadPlot.errorbar(PsinAvg[Shot]+PsinOffset,TemidAvg[Shot],yerr=ErrTe[Shot],fmt='o',color=COLORS[n],markersize=7,linewidth=3,capsize=7)
    elif SHOT == 25:
        Gas025=SOLPSPLOT('25',ATTEMPTS,Publish=PUBLISH,Markers=False,PlotScheme=COLORS,EXP=False)
        
        Gas025.RadProf('Ne',AX=NeRadPlot)
        Gas025.RadProf('Te',AX=TeRadPlot)
        
        for n, Shot in enumerate(ExpID):
            NeRadPlot.errorbar(PsinAvg[Shot]+PsinOffset,NemidAvg[Shot],yerr=ErrNe[Shot],fmt='o',color=COLORS[n][0],markersize=7,linewidth=3,capsize=7)
            TeRadPlot.errorbar(PsinAvg[Shot]+PsinOffset,TemidAvg[Shot],yerr=ErrTe[Shot],fmt='o',color=COLORS[n][0],markersize=7,linewidth=3,capsize=7)
    
    NeRadPlot.set_ylabel(r'Outer Midplane Electron Density $n_e\;(m^{-3})$')
    #NeRadPlot.set_ylabel('')
    
    TeRadPlot.set_ylabel(r'Outer Midplane Electron Temperature $T_e\;(eV)$')
    #TeRadPlot.set_ylabel('')

if JND == 1:
    NeuDenfig, NeuDRadPlot = plt.subplots(nrows=1, ncols=1)
    NeuDRadPlot.semilogy(Jexp['R__R_LCFS_m'],Jexp['N_D_m3'],'--')
    NeuDRadPlot.semilogy(Jkn1d_H['R__R_LCFS_m'],Jkn1d_H['NH_m3'],':')
    NeuDRadPlot.semilogy(Jkn1d_H2['R__R_LCFS_m'],Jkn1d_H2['NH2_m3'],':')
    '''
    Base012.RadProf('NeuDen',LOG10=2,AX=NeuDRadPlot,Markers=False,PlotScheme=['r-','k-'],RADC='rrsep')
    Base025.RadProf('NeuDen',LOG10=2,AX=NeuDRadPlot,Markers=False,PlotScheme=['b-','g-'],RADC='rrsep')
    '''
    NeuDen_2 = SOLPSPLOT('12',[120],RADC='rrsep')
    NeuDen_3 = SOLPSPLOT('25',[171],RADC='rrsep')
    NeuDen_4 = SOLPSPLOT('12',[101],RADC='rrsep')
    NeuDen_5 = SOLPSPLOT('25',[174],RADC='rrsep')
    #NeuDen_6 = SOLPSPLOT('1040122027',[6],ROOTSHOT='',JXA=59,JXI=35,RADC='rrsep')
    #NeuDen_7 = SOLPSPLOT('1040122027',[7],ROOTSHOT='',JXA=59,JXI=35,RADC='rrsep')
    #NeuDen025F = SOLPSPLOT('gas025',[19])
    NeuDen_2.RadProf('NeuDen',LOG10=2,AX=NeuDRadPlot,Markers=False,PlotScheme=['r-'])
    NeuDen_3.RadProf('NeuDen',LOG10=2,AX=NeuDRadPlot,Markers=False,PlotScheme=['b-'])
    NeuDen_4.RadProf('NeuDen',LOG10=2,AX=NeuDRadPlot,Markers=False,PlotScheme=['g-'])    
    NeuDen_5.RadProf('NeuDen',LOG10=2,AX=NeuDRadPlot,Markers=False,PlotScheme=['k-'])  
    #NeuDen_6.RadProf('NeuDen',LOG10=2,AX=NeuDRadPlot,Markers=False)  
    #NeuDen_7.RadProf('NeuDen',LOG10=2,AX=NeuDRadPlot,Markers=False)   
    #NeuDen025F.RadProf('NeuDen',LOG10=2,AX=NeuDRadPlot,RADC='rrsep',Markers=False)
    
    NeuDRadPlot.legend(['Exp D Density', 'KN1D_H', 'KN1D_H2','UnBallooned High Opacity','Ballooned High Opacity','UnBallooned Low Opacity','Ballooned Low Opacity'])
    
    NeuDRadPlot.grid()




