# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 15:34:30 2018

@author: Rmreksoatmodjo
"""
#import os
import numpy as np
import xarray as xr
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import colors, cm
from scipy.io import loadmat
#import geqdsk
import equilibrium as eq
#from D3DPreProcess import PsiNtoR
#from D3DPreProcess import RhotoPsiN

def SOLPSPlotter(Shot, Attempts, Parameter, *PLOTTYPE,
                 TimeRange=[0.90,1.00],  
                 EXP=1,
                 LOG10=0,
                 GRAD=0,
                 ELEV=75,
                 AZIM=270,
                 JXI=37,
                 JXA=55,
                 SEP=20,
                 XDIM=98,
                 YDIM=38,
                 CoreSize=48,
                 CoreBound=[24,71],
                 Publish = [],
                 PsinOffset=0,
                 RadOffset=0,
                 RRad='psin',
                 RadSel=None,
                 PolSel=None, 
                 Device='cmod',
                 Geo=1,
                 LVN=100,
                 DivReg=0,
                 Save=0,
                 COMPARE=0,
                 TC_Flux=[], 
                 TC_Psin=[],
                 ax=None): 
    '''  
    Plot SOLPS data!!!
    
    REQUIRED:
    ** Setup **
    Shot = Shot Number (Either 12 or 25) or 'd3d' for d3d experiment
    Attempts = Simulation Attempt Number(s) [As a List!]
    Parameter = Parameter to be plotted 
        (Specified by a string:  
        'NeuDen' - Neutral Density
        'IonFlx' - Ionization/Radial Particle Flux
        'Ne' - Electron Density
        'Te' - Electron Temperature
        'DN' - Particle Diffusion Coefficient
        'KYE' - Electron Thermal Diffusion Coefficient
        'KYI' - Ion Thermal Diffusion Coefficient  )  
    ** Plotting Flags **
    'Contour' > Plot 2D spatial Contour Plot
    'Surface' > Plot 3D spatial Suface Plot
    'PolPlot' > Plot Poloidal Contours for every Radial grid point (ONLY TAKES 1 ATTEMPT INPUT) 
    'RadPlot' > Plot Radial Contours for every Poloidal grid point (ONLY TAKES 1 ATTEMPT INPUT)
    'RadProf' > Plot Radial Profile of Parameter at Outer Midplane
    'SumPlot' > Plot Poloidally-Summed Parameter values vs Radial location
    'VeslMesh' > Plot SOLPS DG Mesh Grid
    'Export' > No Plot, Save Data to a Variable [PARAM, RR]
    
    OPTIONAL:
    
    TimeRange = [0.90,1.00] > Time range (in sec) over which experimental data is averaged
    JXI = 37 > Poloidal position of Inner Midplane - default is 37
    JXA = 55 > Poloidal position of Outer Midplane - default is 55
    SEP = 18 > Radial position of SEParatrix - default is 18
    AZIM = 270 > Azimuthal viewing angle for Surface Plot
    ELEV = 75 > Elevation of viewing angle
    LOG10 = 0 > Plot Base-10 Logarithm of Parameter Data
    GRAD = 0 > Calculate Gradient of Parameter Data
    RRad = 'psin' > Set radial coordinate convention - Either 'psin', 'radial', or 'Y'
    Geo = 1 > Map Contour to Vessel Geometry 
    DivReg = 0 > Include Divertor Region Data (May cause loss of resolution)
    Save = 0 > Save plot to a .png file
    '''
    #Function begins here
    
    # Parse Experimental Data
        
    if Shot=='d3d':
        GFILE = 'gfileProcessing/d3d_files/g175060.02512'
        GF = eq.equilibrium(gfile=GFILE)
        
        ExpData = loadmat('gfileProcessing/d3d_files/175060_data_SOLPS.mat')
        TC = loadmat('gfileProcessing/d3d_files/flow_transport.mat')
        
        ii = 0
        jj = 0
        kk = 0
        
        Z_mid = 0
        R_OMcore = 2.176
        psin_core = GF.psiN(R_OMcore,Z_mid)
        
        while ExpData['psin_ne'][ii] < psin_core:
            ii += 1
        while ExpData['psin_te'][jj] < psin_core:
            jj += 1 
        while ExpData['psin_ti'][kk] < psin_core:
            kk += 1
        
        PsinNe = ExpData['psin_ne'][ii:] + PsinOffset
        PsinTe = ExpData['psin_te'][jj:] + PsinOffset
        PsinTi = ExpData['psin_ti'][kk:] + PsinOffset
        PsinAvg = [PsinNe,PsinTe,PsinTi]
        
        Ned3d = ExpData['Ne'][ii:]
        Ted3d = ExpData['Te'][jj:]
        Tid3d = ExpData['Ti'][kk:]
        
        JXI=40
        JXA=56
        
    else:
        GFILE = 'gfileProcessing/' + str(Device) + '_files/g11607180' + str(Shot) + '.01209_974'
        GF = eq.equilibrium(gfile=GFILE)
        
        ExpFile = '11607180' + str(Shot)
        ExpData = loadmat(ExpFile)    
        
        ti = 0
        while TimeRange[0] > ExpData['time'][ti]:
            ti = ti+1           
        tf = 0
        while TimeRange[1] > ExpData['time'][tf]:
            tf = tf+1
        
        Psin = np.transpose(ExpData['psin'])[:,ti:tf]
        PsinAvg = np.mean(Psin, axis=1) + PsinOffset
        
        Rmid = np.transpose(ExpData['rmid'])[:,ti:tf]
        RmidAvg = np.mean(Rmid, axis=1) + RadOffset
        
        Nemid = np.transpose(ExpData['ne'])[:,ti:tf]
        Nemid[Nemid == 0] = np.nan
        NemidAvg = np.nanmean(Nemid, axis=1)
        HiErrNe = np.nanmax(Nemid,axis=1) - NemidAvg
        LoErrNe = NemidAvg - np.nanmin(Nemid,axis=1)
        ErrNe = [LoErrNe,HiErrNe]
        NeThresh = (ErrNe[1]-ErrNe[0])/NemidAvg
        for NT in range(len(NeThresh)):
            if np.abs(NeThresh[NT]) > 0.5:
                NemidAvg[NT] = np.nan
                ErrNe[0][NT] = np.nan
                ErrNe[1][NT] = np.nan
        
        Temid = np.transpose(ExpData['te'])[:,ti:tf]
        Temid[Temid == 0] = np.nan
        TemidAvg = np.nanmean(Temid, axis=1)
        HiErrTe = np.nanmax(Temid,axis=1) - TemidAvg
        LoErrTe = TemidAvg - np.nanmin(Temid,axis=1)
        ErrTe = [LoErrTe,HiErrTe]
        TeThresh = (ErrTe[1]-ErrTe[0])/TemidAvg
        for TT in range(len(TeThresh)):
            if np.abs(TeThresh[TT]) > 0.5:
                TemidAvg[TT] = np.nan
                ErrTe[0][TT] = np.nan
                ErrTe[1][TT] = np.nan
    
    N = len(Attempts)
    
    if Publish!=[]:
        plt.rc('font',size=20)
        plt.rc('lines',linewidth=5,markersize=15)
    else:
        plt.rc('font',size=14)
    
    if DivReg==1:
        XGrid=XDIM-2
        XMin=0
        XMax=XGrid-1
    else:
        XGrid=CoreSize
        XMin=CoreBound[0]
        XMax=CoreBound[1]
    YSurf=YDIM-2
    
    X = np.linspace(XMin,XMax,XGrid)
    Y = np.linspace(0,YSurf-1,YSurf)
    Xx, Yy = np.meshgrid(X,Y)   # Create X and Y mesh grid arrays
    
    CMAP = cm.viridis
    
    YYLoc = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Radial Grid Point $N$')
    RadLoc = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Radial Coordinate $m$')
    VertLoc = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Vertical Coordinate $m$')
    RadCor = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Corner Radial Coordinate $m$')
    VertCor = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Corner Vertical Coordinate $m$')
    PsinLoc = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Normalized Psi $\psi_N$')
    
    for n in range(N):
        Attempt = Attempts[n]
        dir = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output'     #Generate path
        dir2 = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output2'     #Generate Mesh path
        
        YYLoc.values[:,:,n] = Yy
        RadLoc.values[:,:,n] = np.loadtxt(dir + '/RadLoc' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
        VertLoc.values[:,:,n] = np.loadtxt(dir + '/VertLoc' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
        
        try:
            RadCor.values[:,:,n] = np.loadtxt(dir2 + '/Rad0Cor' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
            VertCor.values[:,:,n] = np.loadtxt(dir2 + '/Vert0Cor' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
        except:
            print("Warning, Grid Corner Coordinates Not Found for Attempt" + str(Attempt))
   
        for j in range(len(Y)):
            for i in range(len(X)):
                PsinLoc.values[j,i,n] = GF.psiN(RadLoc.loc[Y[j],X[i],Attempt].values,VertLoc.loc[Y[j],X[i],Attempt].values,)
    
    if RRad == 'psin':
        RR = PsinLoc
        Rexp = PsinAvg
        Rstr = '$\psi_N$'
    elif RRad == 'radial':
        RR = RadLoc
        if Shot != 'd3d':
            Rexp = RmidAvg
        Rstr = 'm'
    elif RRad == 'Y':
        RR = YYLoc
        Rexp = None
        Rstr = 'N'
    else:
        print('Invalid Radial Coordinate specified')
    #NeLast10 = xr.DataArray(np.loadtxt(dir + '/ne3da.last10'))
    
    if Parameter == 'VOL':
        VOL = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Cell Volume VOL $(m^{2})$') 
        for n in range(N):
            Attempt = Attempts[n]
            dir = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output'     #Generate path
            VOL.values[:,:,n] = np.loadtxt(dir + '/VOL' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
        PARAM = VOL
    if Parameter == 'SX':
        SX = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Poloidal Contact Area SX $(m^{2})$') 
        for n in range(N):
            Attempt = Attempts[n]
            dir = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output2'     #Generate path
            SX.values[:,:,n] = np.loadtxt(dir + '/SX' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
        PARAM = SX
    if Parameter == 'SY':
        SY = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Radial Contact Area SY $(m^{2})$') 
        for n in range(N):
            Attempt = Attempts[n]
            dir = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output2'     #Generate path
            SY.values[:,:,n] = np.loadtxt(dir + '/SY' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
        PARAM = SY        
    if Parameter == 'SZ':
        SZ = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Poloidal Cross-Sectional Area SZ $(m^{2})$') 
        for n in range(N):
            Attempt = Attempts[n]
            dir = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output2'     #Generate path
            SZ.values[:,:,n] = np.loadtxt(dir + '/SZ' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
        PARAM = SZ
    if Parameter == 'NeuDen':
        NeuDen = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Neutral Density $(m^{-3})$')     #Load Neutral Density data 
        for n in range(N):
            Attempt = Attempts[n]
            dir = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output'     #Generate path
            NeuDen.values[:,:,n] = np.loadtxt(dir + '/NeuDen' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
        PARAM = NeuDen
    if Parameter == 'RadFlx':
        RadFlx = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Radial Atomic Flux $m^{-2}s^{-1}$')     #Load Radial Atomic Flux data 
        for n in range(N):
            Attempt = Attempts[n]
            dir = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output'     #Generate path
            RadFlx.values[:,:,n] = np.loadtxt(dir + '/RadFlx' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
        PARAM = RadFlx
    if Parameter == 'MolFlx':
        MolFlx = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Radial Molecule Flux $m^{-2}s^{-1}$')     #Load Radial Molecule Flux data 
        for n in range(N):
            Attempt = Attempts[n]
            dir = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output'     #Generate path
            MolFlx.values[:,:,n] = np.loadtxt(dir + '/MolFlx' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
        PARAM = MolFlx
    if Parameter == 'IonFlx':
        IonFlx = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Radial Particle Flux $s^{-1}$')     #Load Radial Particle Flux data 
        for n in range(N):
            Attempt = Attempts[n]
            dir = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output'     #Generate path
            IonFlxS = np.loadtxt(dir + '/IonFlx' + str(Attempt),usecols = (3)).reshape((2*YDIM,XDIM))[1+YDIM:2*YDIM-1,XMin+1:XMax+2] 
            #dir2 = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output2'
            #SY = np.loadtxt(dir2 + '/SY' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2] 
            IonFlx.values[:,:,n] = IonFlxS
        PARAM = IonFlx
    if Parameter == 'IonPol':
        IonPol = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Ionization/Poloidal Particle Flux $m^{-2}s^{-1}$')     #Load Ion Flux data 
        for n in range(N):
            Attempt = Attempts[n]
            dir = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output'     #Generate path
            IonPolS = np.loadtxt(dir + '/IonPol' + str(Attempt),usecols = (3)).reshape((2*YDIM,XDIM))[1+YDIM:2*YDIM-1,XMin+1:XMax+2] 
            #dir2 = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output2'
            #SX = np.loadtxt(dir2 + '/SX' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2] 
            IonPol.values[:,:,n] = IonPolS
        PARAM = IonPol
    if Parameter == 'Ne':
        Ne = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Electron Density $n_e (m^{-3})$')      #Load Electron Density data 
        for n in range(N):
            Attempt = Attempts[n]
            dir = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output'     #Generate path
            Ne.values[:,:,n] = np.loadtxt(dir + '/Ne' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
        PARAM = Ne
    if Parameter == 'Te':
        Te = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Electron Temperature $T_e (eV)$')      #Load Electron Temp data 
        for n in range(N):
            Attempt = Attempts[n]
            dir = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output'     #Generate path
            Te.values[:,:,n] = np.loadtxt(dir + '/Te' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
        PARAM = Te
    if Parameter == 'Ti':
        Ti = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Ion Temperature $T_i (eV)$')      #Load Ion Temp data 
        for n in range(N):
            Attempt = Attempts[n]
            dir = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output'     #Generate path
            Ti.values[:,:,n] = np.loadtxt(dir + '/Ti' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
        PARAM = Ti        
    if Parameter == 'DN':
        DN = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Particle Density Diffusivity $D\;(m^2/s)$')     #Load Particle Diffusion Coeff data 
        for n in range(N):
            Attempt = Attempts[n]
            dir = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output'     #Generate path
            DN.values[:,:,n] = np.loadtxt(dir + '/DN' + str(Attempt),usecols = (3)).reshape((2*YDIM,XDIM))[1+YDIM:2*YDIM-1,XMin+1:XMax+2]
        PARAM = DN
    if Parameter == 'KYE':
        KYE = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Electron Thermal Diffusivity $\chi_e (m^2/s)$')      #Load Heat Transfer Coeff data
        for n in range(N):
            Attempt = Attempts[n]
            dir = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output'     #Generate path
            KYE.values[:,:,n] = np.loadtxt(dir + '/KYE' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
        PARAM = KYE
    if Parameter == 'KYI':
        KYI = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Ion Thermal Diffusivity $\chi_i (m^2/s)$')      #Load Heat Transfer Coeff data
        for n in range(N):
            Attempt = Attempts[n]
            dir = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output'     #Generate path
            KYI.values[:,:,n] = np.loadtxt(dir + '/KYI' + str(Attempt),usecols = (3)).reshape((2*YDIM,XDIM))[1+YDIM:2*YDIM-1,XMin+1:XMax+2]
        PARAM = KYI        
    if Parameter == 'VLY':
        VLY = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Radial Pinch $v_y (m^2/s)$')      #Load Heat Transfer Coeff data
        for n in range(N):
            Attempt = Attempts[n]
            dir = 'SOLPS_2D_prof/Shot0' + str(Shot) + '/Attempt' + str(Attempt) + '/Output'     #Generate path
            VLY.values[:,:,n] = np.loadtxt(dir + '/VLY' + str(Attempt),usecols = (3)).reshape((2*YDIM,XDIM))[1+YDIM:2*YDIM-1,XMin+1:XMax+2]
        PARAM = VLY

    if COMPARE == 1:
        CMAP = cm.coolwarm
        for n in np.arange(1,N):
            PARAM.values[:,:,n] = (PARAM.values[:,:,0]-PARAM[:,:,n])
        PARAM.values[:,:,0] = PARAM.values[:,:,0]-PARAM[:,:,0]

    if LOG10 == 1:
        #PARAM.values[PARAM.values<0] = 0
        PARAM.values[PARAM.values>1] = np.log10(PARAM.values[PARAM.values>1])
        PARAM.values[PARAM.values<-1] = -1*np.log10(np.abs(PARAM.values[PARAM.values<-1]))    
        y_exp = np.arange(np.floor(np.nanmin(PARAM.values)), np.ceil(np.nanmax(PARAM.values))+1,2)
        
    if RadSel == 'all':
        RadSel = PARAM.coords['Radial_Location'].values
    if RadSel == None:
        RadSel = [SEP]
    if PolSel == 'all':
        PolSel = PARAM.coords['Poloidal_Location'].values
    if PolSel == None:
        PolSel = [JXI,JXA]
    
    if 'VeslMesh' in PLOTTYPE and N == 1:
        fig = plt.figure(figsize=(14,10))
        plt.plot(RadCor.loc[:,:,Attempts[0]],VertCor.loc[:,:,Attempts[0]])
        plt.plot(np.transpose(RadCor.loc[:,:,Attempts[0]].values),np.transpose(VertCor.loc[:,:,Attempts[0]].values))
        plt.plot(RadLoc.loc[:,JXA,Attempts[0]],VertLoc.loc[:,JXA,Attempts[0]],color='Orange')
        plt.plot(RadLoc.loc[:,JXI,Attempts[0]],VertLoc.loc[:,JXI,Attempts[0]],color='Red')
        plt.plot(RadLoc.loc[SEP,:,Attempts[0]],VertLoc.loc[SEP,:,Attempts[0]],color='Black')
        plt.title('Vessel Mesh Geometry')
        plt.xlabel('Radial Coordinate r (m)')
        plt.ylabel('Vertical Coordinate z (m)')
        plt.gca().set_aspect(1.0)
        plt.grid()
    
    if 'Contour' in PLOTTYPE:
        if LOG10 == 2:
            NPARAM = np.abs(PARAM.values[PARAM.values<0])
            NPARAM[NPARAM<=1] = np.nan
            PARAM.values[np.abs(PARAM.values)<=1] = np.nan
            if NPARAM.size>0:            
                lev_exp = np.arange(-(np.ceil(np.log10(np.nanmax(NPARAM)))+1), np.ceil(np.log10(np.nanmax(PARAM.values)))+1,5)
                levs = np.sign(lev_exp)*np.power(10, np.abs(lev_exp))
                np.set_printoptions(threshold=np.inf)
                print(levs)
            else:
                lev_exp = np.arange(np.floor(np.log10(np.nanmin(PARAM.values)))-1, np.ceil(np.log10(np.nanmax(PARAM.values)))+1)
            levs = np.power(10, lev_exp)
        elif LOG10 == 1:
            levs = np.arange(np.floor(PARAM.values.min()),np.ceil(PARAM.values.max())) #LVN
        else:
            levs = np.linspace(np.floor(PARAM.values.min()),np.ceil(PARAM.values.max()),LVN)
        for n in range(N):
            if ax is None:
                fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(14,10))
            
            if Geo == 1:
                if LOG10 == 2:
                    IM = ax.contourf(RadLoc.values[:,:,n],VertLoc.values[:,:,n],PARAM.values[:,:,n],levs,cmap=CMAP,norm=colors.LogNorm())
                else:
                    IM = ax.contourf(RadLoc.values[:,:,n],VertLoc.values[:,:,n],PARAM.values[:,:,n],levs,cmap=CMAP)
                ax.plot(RadLoc.values[:,(JXA-XMin),n],VertLoc.values[:,(JXA-XMin),n],color='Orange',linewidth=3)
                ax.plot(RadLoc.values[:,(JXI-XMin),n],VertLoc.values[:,(JXI-XMin),n],color='Red',linewidth=3)
                ax.plot(RadLoc.values[SEP,:,n],VertLoc.values[SEP,:,n],color='Black',linewidth=3)
                ax.set_xlabel('Radial Location (m)')
                ax.set_ylabel('Vertical Location (m)')                
            else:
                if LOG10 == 2:
                    IM = ax.contourf(Xx[:,:],Yy[:,:],PARAM.values[:,:,n],levs,norm=colors.LogNorm(),cmap=CMAP)
                else:
                    IM = ax.contour(Xx[:,:],Yy[:,:],PARAM.values[:,:,n],levs,cmap=CMAP)
                ax.plot(Xx[:,(JXA-XMin-1)],Yy[:,(JXA-XMin-1)],color='Orange',linewidth=3)
                ax.plot(Xx[:,(JXI-XMin-1)],Yy[:,(JXI-XMin-1)],color='Red',linewidth=3)
                ax.plot(Xx[SEP,:],Yy[SEP,:],color='Black',linewidth=3)
                ax.set_xlabel('Poloidal Coordinate')
                ax.set_ylabel('Radial Coordinate') 
                
            #plt.legend(('Outboard Midplane','Inboard Midplane','Separatrix'),fontsize='xx-large')
            if Publish!=[]:
                ax.set_title('Attempt ' + Publish[n] + ' ' + PARAM.name)
            else:
                ax.set_title('Discharge 0' + str(Shot) + ' Attempt ' + str(Attempts[n]) + ' ' + PARAM.name)
            plt.colorbar(IM)
            ax.set_aspect('equal')
            ax.grid(b=1)
            #a.set_xticklabels(['%.1f' % i for i in a.get_xticks()], fontsize='x-large')
            #a.set_yticklabels(['%.1f' % j for j in a.get_yticks()], fontsize='x-large')
            #a.tick_params(labelsize=20)
            if Save == 1:
                ImgName = 'Profiles/' + Parameter + str(Attempts[n]) + 'Contour.png'
                ax.savefig(ImgName, bbox_inches='tight')
    
    if 'Surface' in PLOTTYPE:
        if LOG10 == 2:
            PARAM.values[PARAM.values==0] = np.nan
            PARAM.values = np.log10(PARAM.values)
        Zmax = np.nanmax(PARAM.values) 
        for n in range(N):
            fig = plt.figure(figsize=(18,12))  # Size is width x height
            ax1 = fig1.gca(projection='3d')
        
            SURF = ax1.plot_surface(RadLoc.values[:,:,n],VertLoc.values[:,:,n],PARAM.values[:,:,n],cmap=CMAP,vmax=Zmax)
            OMP = ax1.plot(RadLoc.values[:,(JXA-XMin),n],VertLoc.values[:,(JXA-XMin),n],PARAM.loc[:,JXA,Attempts[n]],color='Orange')
            IMP = ax1.plot(RadLoc.values[:,(JXI-XMin),n],VertLoc.values[:,(JXI-XMin),n],PARAM.loc[:,JXI,Attempts[n]],color='Red')
            SEP = ax1.plot(RadLoc.values[SEP,:,n],VertLoc.values[SEP,:,n],PARAM.loc[SEP,:,Attempts[n]],color='Black')
            plt.legend(('Outboard Midplane','Inboard Midplane','SEParatrix'))            
            
            ax1.view_init(ELEV,AZIM)     # Use (30,135) for Ne,Te,DN,KYE,IonFlx. Use (30,315) for NeuDen
            ax1.set_zlabel(PARAM.name)
            ax1.set_zbound(upper=Zmax)
            plt.title('Discharge 0' + str(Shot) + ' Attempt(s) ' + str(Attempts[n]) + ' ' + PARAM.name)
            plt.xlabel('Radial Location (m)')
            plt.ylabel('Poloidal Location (m)')
            plt.colorbar(SURF)
            if Save == 1:
                ImgName = 'Profiles/' + Parameter + 'SurfA' + str(Attempts[n]) + '.png'
                plt.savefig(ImgName, bbox_inches='tight')
    
    #if 'PolPlot' in PLOTTYPE and N > 1: print('Poloidal Plot not compatible with multiple Attempts!')
    if 'PolPlot' in PLOTTYPE:# and N == 1:
        #with plt.xkcd():
            fig = plt.figure(figsize=(14,7))
            PolVal = PARAM.loc[RadSel[0],:,:].values
            PolVal[PolVal==0]=np.nan
            if LOG10 == 2:
                plt.semilogy(X,PolVal)
            else:
                plt.plot(X,PolVal)
            if Publish==[]:
                plt.legend(Attempts)
                plt.title('Discharge 0' + str(Shot) + ' Attempt(s) ' + str(Attempts) + ' Poloidal ' + PARAM.name)
            else:
                plt.legend(Publish,loc=4)
                plt.title('Poloidal ' + PARAM.name + ' along Separatrix')
            Pmin = np.nanmin(PolVal)
            Pmax = np.nanmax(PolVal)
            plt.plot([int(JXI), int(JXI)],[Pmin, Pmax],color='Red')
            plt.plot([int(JXA), int(JXA)],[Pmin, Pmax],color='Orange')
            plt.xlabel('Poloidal Coordinate')
            plt.ylabel(r'$Log_{10}$ of ' + PARAM.name)
            plt.grid()
    
    if 'RadPlot' in PLOTTYPE and N > 1: print('Radial Plot not compatible with multiple Attempts!')
    if 'RadPlot' in PLOTTYPE and N == 1:
            fig = plt.figure(figsize=(14,10))
            if GRAD == 1:
                PARAM.values = np.gradient(PARAM.values,axis=0)
            if LOG10 == 2:
                plt.semilogy(RR.loc[:,PolSel,Attempts[0]].values, PARAM.loc[:,PolSel,Attempts[0]].values)
            else:
                plt.plot(RR.loc[:,PolSel,Attempts[0]].values, PARAM.loc[:,PolSel,Attempts[0]].values)
            if PolSel==[JXI,JXA]:
                plt.legend(['Inner Midplane','Outer Midplane'])
            else:
                plt.legend(PolSel)
            Pmin = float(PARAM.values.min())
            Pmax = float(PARAM.values.max())
            plt.plot([RR.loc[SEP,JXA,Attempt], RR.loc[SEP,JXA,Attempt]],[Pmin, Pmax],color='Black')
            if LOG10 == 1:
                plt.yticks(y_exp)
                labels = ['10^' + str(int(ex)) for ex in y_exp]
                plt.gca().set_yticklabels(labels)
            
            plt.title('Discharge 0' + str(Shot) + ' Attempt(s) ' + str(Attempts) + ' Radial ' + PARAM.name)
            plt.xlabel('Radial Coordinate ' + Rstr)
            plt.ylabel(PARAM.name)
            plt.grid()
    
    if 'RadProf' in PLOTTYPE:
        #with plt.xkcd():
            if ax is None: # if no axis was supplied to the function create our own
                fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(14,7))
                
            if Parameter == 'RadFlx' and EXP == 1 and Shot == 'd3d':
                ax.plot(TC_Psin,TC_Flux,markersize=7,linewidth=3)
            if GRAD == 1:
                PARAM.values = np.gradient(PARAM.values,axis=0)
            if LOG10 == 2:
                ax.semilogy(RR.loc[:,JXA,:], PARAM.loc[:,JXA,:])
            else:
                ax.plot(RR.loc[:,JXA,:], PARAM.loc[:,JXA,:],linewidth=3)
            #ax.plot(RadLoc.loc[:,JXI], PARAM.loc[:,JXI],color='Red')
            if Parameter == 'Ne' and EXP == 1 and Shot != 'd3d':
                ax.errorbar(Rexp,NemidAvg,yerr=ErrNe,fmt='bo',markersize=7,linewidth=3,capsize=7)
            if Parameter == 'Ne' and EXP == 1 and Shot == 'd3d':
                ax.plot(PsinNe,Ned3d,'bo',markersize=7,linewidth=3)
            if Parameter == 'Te' and EXP == 1 and Shot != 'd3d':
                ax.errorbar(Rexp,TemidAvg,yerr=ErrTe,fmt='bo',markersize=7,linewidth=3,capsize=7)
            if Parameter == 'Te' and EXP == 1 and Shot == 'd3d':
                ax.plot(PsinTe,Ted3d,'bo',markersize=7,linewidth=3)
            if Parameter == 'Ti' and EXP == 1 and Shot == 'd3d':
                ax.plot(PsinTi,Tid3d,'bo',markersize=7,linewidth=3)
            ax.set_xlabel('Radial Coordinate ' + Rstr)
            if LOG10 == 1:
                ax.set_ylabel('Log_10 of ' + PARAM.name)
            else:    
                ax.set_ylabel(PARAM.name)
            if Publish==[]:
                ax.legend(Attempts)
                ax.set_title('Discharge 0' + str(Shot) + ' Attempt(s) ' + str(Attempts) + ' Midplane Radial ' + PARAM.name)
            else:
                ax.legend(Publish,loc=2)
                ax.set_title('Midplane Radial ' + PARAM.name)
            Pmin = float(PARAM.loc[:,JXA,:].min())
            Pmax = float(PARAM.loc[:,JXA,:].max())
            ax.plot([RR.loc[SEP,JXA,Attempt], RR.loc[SEP,JXA,Attempt]],[Pmin, Pmax],color='Black',linewidth=3)
            ax.grid(b=1)
#            a = ax.gca()
#            a.set_xticklabels(['%.2f' % i for i in a.get_xticks()], fontsize='x-large')
#            a.set_yticklabels(['%.2e' % j for j in a.get_yticks()], fontsize='x-large')
    
    if 'SumPlot' in PLOTTYPE:
        #with plt.xkcd():
            fig = plt.figure(figsize=(14,7))
            SumParam = np.sum(PARAM.values, axis=1)
            if LOG10 == 2:
                plt.semilogy(RR.loc[:,JXA,:], SumParam)
            else:
                plt.plot(RR.loc[:,JXA,:], SumParam)
            if Publish==[]:
                plt.legend(Attempts)
                plt.title('Discharge 0' + str(Shot) + ' Attempt(s) ' + str(Attempts) + ' Poloidal Summation of ' + PARAM.name)
            else:
                plt.legend(Publish,loc=2)
                plt.title('Poloidal Summation of ' + PARAM.name)    
            plt.plot([RR.loc[SEP,JXA,Attempt], RR.loc[SEP,JXA,Attempt]],[float(SumParam.min()), float(SumParam.max())],color='Black')
            plt.xlabel('Radial Coordinate ' + Rstr)
            if LOG10 == 1 or LOG10 == 2:
                plt.ylabel('Log_10 of ' + PARAM.name)
            else:    
                plt.ylabel(PARAM.name)
            #plt.ylabel(PARAM.name)
            #plt.grid()
    
    plt.show()
    
    if 'Export' in PLOTTYPE:
        Results = {}
        Results['PARAM'] = PARAM
        Results['RR'] = RR
        Results['RadLoc'] = RadLoc
        Results['VertLoc'] = VertLoc
        if Shot != 'd3d':
            Results['NemidAvg'] = NemidAvg
            Results['TemidAvg'] = TemidAvg
            Results['Rexp'] = Rexp
            Results['ErrNe'] = ErrNe
            Results['ErrTe'] = ErrTe
        return Results
    else:
        return Publish
