# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 14:00:00 2019

@author: Rmreksoatmodjo
"""
import os
import numpy as np
import xarray as xr
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import colors, cm
from scipy.io import loadmat
import geqdsk
import equilibrium as eq
from D3DPreProcess import PsiNtoR
from D3DPreProcess import RhotoPsiN
from TOOLS import SET_WDIR

class SOLPSPLOT(object):
        
    '''  
    Plot SOLPS data!!! 2.0!!!
    
    REQUIRED:
    ** Setup **
    Shot = Shot Number (Either 12 or 25) or 'd3d' for d3d experiment
    Attempts = Simulation Attempt Number(s) [As a List!]
    Parameters = Parameters to be plotted 
        (Specified by a string:  
        'NeuDen' - Neutral Density
        'IonFlx' - Ionization/Radial Particle Flux
        'Ne' - Electron Density
        'Te' - Electron Temperature
        'DN' - Particle Diffusion Coefficient
        'KYE' - Electron Thermal Diffusion Coefficient
        'KYI' - Ion Thermal Diffusion Coefficient  )  
    
    OPTIONAL:

    SAVE = False > Save plot to a .png file        
    GRID = False > Turn on plot grid
    SUBTRACT = False > For a series of Attempts, subtract the matrix data of each Attempt from the first declared Attempt, to compare difference 
    GRAD = False > Calculate Gradient of Parameter Data    
    EXP = True > Plot Experimental Data Points
    GEO = True > Map Contour to Vessel Geometry 
    DIVREG = True > Include Divertor Region Data (May cause loss of resolution)        
    LOG10 = 0 > Plot Base-10 Logarithm of Parameter Data (0, 1, or 2)
    PsinOffset = 0 > Psi_n offset of Experimental Data Points
    RadOffset = 0 > Radial (in meters) offset of Experimental Data Points
    JXI = 37 > Poloidal position of Inner Midplane - default is 37
    JXA = 55 > Poloidal position of Outer Midplane - default is 55
    SEP = 20 > Radial position of SEParatrix - default is 18
    AZIM = 270 > Azimuthal viewing angle for Surface Plot
    ELEV = 75 > Elevation of viewing angle
    XDIM = 98 > Dimension of computational grid in the x (poloidal) direction
    YDIM = 38 > Dimension of computational grid in the y (radial) direction
    LVN = 100 > Number of colorbar levels for contour plots
    CoreBound = [25,72] > X-coordinate grid cell numbers that define the [left, right] bounds of Core Region
    TimeRange = [0.90,1.00] > Time range (in sec) over which experimental data is averaged
    Publish = [] > List of strings to use in legends of publication-quality plots; if not [], changes plotting rc.params 
    RADC = 'psin' > Set radial coordinate convention - Either 'psin', 'radial', 'rrsep' or 'Y'    
    BASEDRT = 'SOLPS_2D_prof/' > Local home directory
    RadSlc = None > Radial surface selection for poloidal plots - Can set specific radial index, 'all', or 'None' defaults to SEP
    PolSlc = None > Poloidal grid line selection for radial plots - Can set specific poloidal index, 'all', or 'None' defaults to [JXA, JXI]
    AX = None > Pass the name of a matplotlib axis for the SOLPSPLOT object to plot on; by default SOLPSPLOT plots on a new axis
        
    PLOT TYPES:
        
    SOLPSPLOT.Contour() > Plot 2D spatial Contour Plot
    SOLPSPLOT.Surface() > Plot 3D spatial Suface Plot
    SOLPSPLOT.PolPlot() > Plot Poloidal Contours for every Radial grid point (ONLY TAKES 1 ATTEMPT INPUT) 
    SOLPSPLOT.RadPlot() > Plot Radial Contours for every Poloidal grid point (ONLY TAKES 1 ATTEMPT INPUT)
    SOLPSPLOT.RadProf() > Plot Radial Profile of Parameter at Outer Midplane
    SOLPSPLOT.SumPlot() > Plot Poloidally-Summed Parameter values vs Radial location
    SOLPSPLOT.VeslMesh() > Plot SOLPS DG Mesh Grid
    SOLPSPLOT.Export() > No Plot, Save Data to a Variable [PARAM, RR]
    
    '''
    #Object begins here

    def __init__(self, Shot, Attempts, Parameters=None, **kwargs):
        
        self._reset_object()
        
        self.Shot = Shot
        
        if isinstance(Attempts,list):
            self.Attempts = Attempts
        else:    
            self.Attempts = [Attempts]
        
        if isinstance(Parameters, list):
            self.Parameter = Parameters
        elif Parameters is not None:
            self.Parameter = [Parameters]
        else:
            self.Parameter = ['Ne','Te','Ti','DN','KYE','KYI','NeuDen','IonFlx']
            
        self.DefaultSettings = {'TimeRange' : [1.10,1.30],  
                     'DEV': 'cmod',
                     'EXP' : True,
                     'LOG10' : 0,
                     'GRAD' : False,
                     'ELEV' : 75,
                     'AZIM' : 270,
                     'JXI' : 38,
                     'JXA' : 57,
                     'SEP' : 20,
                     'XDIM' : 98,
                     'YDIM' : 38,
                     'CoreBound' : [25,72],
                     'Publish' : [],
                     'Markers' : True,
                     'PlotScheme' : [],
                     'PsinOffset' : 0,
                     'RadOffset' : 0,
                     'RADC' : 'psin',
                     'POLC' : 'theta',
                     'RadSlc' : None,
                     'PolSlc' : None, 
                     'GEO' : True,
                     'LVN' : 100,
                     'DIVREG' : True,
                     'SAVE' : False,
                     'SUBTRACT' : False,
                     'TC_Flux' : [], 
                     'TC_Psin' : [],
                     'GRID': False,
                     'AX' : None,
                     'BASEDRT': 'solps-iter/runs/',
                     'TOPDRT' : '' }
        
        for key, value in self.DefaultSettings.items():
            if key not in kwargs.keys():
                kwargs[key] = value
                    
        self.KW = kwargs         
        
        self.PARAMDICT = {'Ne': r'Electron Density $n_e\;(m^{-3})$',
                     'Te': r'Electron Temperature $T_e\;(eV)$',
                     'Ti': r'Ion Temperature $T_i\;(eV)$',
                     'DN': r'Particle Density Diffusivity $D\;(m^2/s)$',
                     'KYE': r'Electron Thermal Diffusivity $\chi_e (m^2/s)$',
                     'KYI': r'Ion Thermal Diffusivity $\chi_i (m^2/s)$',
                     'NeuDen': r'Neutral Atom (D) Density $(m^{-3})$',
                     'MolDen': r'Neutral Molecule (D2) Density $(m^{-3})$',
                     'NeuTemp': r'Neutral Atom (D) Temperature (eV)',
                     'MolTemp': r'Neutral Molecule (D2) Temperature (eV)',
                     'IonFlx': r'Radial Particle Flux $s^{-1}$',
                     'MolFlx': r'Radial Molecule Flux $m^{-2}s^{-1}$',
                     'RadFlx': r'Radial Atomic Flux $m^{-2}s^{-1}$',
                     'IonPol': r'Poloidal Particle Flux $m^{-2}s^{-1}$',
                     'VLY': r'Radial Pinch $v_y (m^2/s)$',
                     'SX': r'Poloidal Contact Area SX $(m^{2})$',
                     'SY': r'Radial Contact Area SY $(m^{2})$',
                     'SZ': r'Poloidal Cross-Sectional Area SZ $(m^{2})$',
                     'VOL': r'Cell Volume VOL $(m^{2})$'}
        
        self._LoadSOLPSData()
        
    def _reset_object(self):
        self.Shot=None
        self.Attempts=None
        self.Parameter=[]
        self.PARAM={}
        self.ExpDict={}
        self.RadCoords={}
        
    def _LoadSOLPSData(self, AddNew=None):
        
        Attempts = self.Attempts
        Shot = self.Shot
        Publish = self.KW['Publish']
        DIVREG = self.KW['DIVREG']
        TimeRange = self.KW['TimeRange']
        PsinOffset = self.KW['PsinOffset']        
        RadOffset = self.KW['RadOffset']
        RadSlc = self.KW['RadSlc']
        PolSlc = self.KW['PolSlc']
        SEP = self.KW['SEP']
        XDIM = self.KW['XDIM']
        YDIM = self.KW['YDIM']
        CoreBound = self.KW['CoreBound']
        BASEDRT = self.KW['BASEDRT']
        TOPDRT = self.KW['TOPDRT']
        DEV = self.KW['DEV']
        
        BASEDRT, TOPDRT = SET_WDIR(BASEDRT,TOPDRT)
        
        # Create Experiment Data Dictionary (ExpDict) -> d3d or cmod?
        
        if 'gas' in Shot:
            BASEDRT = '{}gaspuff/'.format(BASEDRT)
        
        if 'd3d' in Shot or DEV=='d3d':
            
            BASEDRT = '{}d3d'.format(BASEDRT)
            
            self.KW['JXI'] = 40
            self.KW['JXA'] = 56
            self.KW['SEP'] = 22
            
            GFILE = '{}gfileProcessing/d3d_files/g175060.02512'.format(TOPDRT)
            GF = eq.equilibrium(gfile=GFILE)
            
            ExpData = loadmat('{}gfileProcessing/d3d_files/175060_data_SOLPS.mat'.format(TOPDRT))
            ONETWO = loadmat('{}gfileProcessing/d3d_files/flow_transport.mat'.format(TOPDRT))
            
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
            
            PsinNe = ExpData['psin_ne'][ii:]
            PsinTe = ExpData['psin_te'][jj:]
            PsinTi = ExpData['psin_ti'][kk:]
            PsinAvg = [PsinNe,PsinTe,PsinTi]
            
            self.ExpDict['Ned3d'] = ExpData['Ne'][ii:]
            self.ExpDict['Ted3d'] = ExpData['Te'][jj:]
            self.ExpDict['Tid3d'] = ExpData['Ti'][kk:]
            
        elif DEV=='cmod':
            
            BASEDRT = '{}cmod/0{}home'.format(BASEDRT, Shot[-2:])
            
            GFILE = '{}gfileProcessing/cmod_files/g11607180{}.01209_974'.format(TOPDRT, Shot[-2:])
            GF = eq.equilibrium(gfile=GFILE)
            
            ExpFile = '11607180{}'.format(Shot[-2:])
            ExpData = loadmat('{}gfileProcessing/cmod_files/{}.mat'.format(TOPDRT, ExpFile))    
            
            ti = 0
            while TimeRange[0] > ExpData['time'].flatten()[ti]:
                ti = ti+1           
            tf = 0
            while TimeRange[1] > ExpData['time'].flatten()[tf]:
                tf = tf+1
            
            Psin = ExpData['psin'][:,ti:tf]
            PsinAvg = np.mean(Psin, axis=1)
            
            Rmid = ExpData['rmid'][:,ti:tf]
            RmidAvg = np.mean(Rmid, axis=1)
            
            Nemid = ExpData['ne'][:,ti:tf]
            Nemid[Nemid == 0] = np.nan
            NemidAvg = np.nanmean(Nemid, axis=1)
            ErrNe = np.nanmean(ExpData['nerr'][:,ti:tf], axis=1)
            NeThresh = (ErrNe*2)/NemidAvg
            for NT in range(len(NeThresh)):
                if np.abs(NeThresh[NT]) > 0.5:
                    NemidAvg[NT] = np.nan
                    ErrNe[NT] = np.nan
            
            Temid = ExpData['te'][:,ti:tf]
            Temid[Temid == 0] = np.nan
            TemidAvg = np.nanmean(Temid, axis=1)
            ErrTe = np.nanmean(ExpData['terr'][:,ti:tf], axis=1)
            TeThresh = (ErrTe*2)/TemidAvg
            for TT in range(len(TeThresh)):
                if np.abs(TeThresh[TT]) > 0.5:
                    TemidAvg[TT] = np.nan
                    ErrTe[TT] = np.nan
                    
            self.ExpDict['NemidAvg'] = NemidAvg
            self.ExpDict['ErrNe'] = ErrNe
            self.ExpDict['TemidAvg'] = TemidAvg
            self.ExpDict['ErrTe'] = ErrTe        
        
        JXA = self.KW['JXA']
        JXI = self.KW['JXI']
        
        #Get dimensions of simulation grid 
        
        N = len(self.Attempts)
        P = len(self.Parameter)
        
        XGrid=XDIM-2
        XMin=1
        XMax=XGrid
        
        YSurf=YDIM-2
        
        X = np.linspace(XMin,XMax,XGrid)
        Y = np.linspace(1,YSurf,YSurf)
        Xx, Yy = np.meshgrid(X,Y)   # Create X and Y mesh grid arrays
        
        RadLoc = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Radial Coordinate $m$')
        VertLoc = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Vertical Coordinate $m$')

        #RadCo = ['YYLoc', 'RRsep', 'PsiNLoc']
        #RadVec = xr.DataArray(np.zeros((YSurf,XGrid,N,3)), coords=[Y,X,Attempts,RadCo], dims=['Radial_Location','Poloidal_Location','Attempt','Radial Metric'], name = 'Radial Coordinate Data')

        YYLoc = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Radial Grid Point $N$')
        PsinLoc = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Normalized Psi $\psi_N$')
        
        #RadCor = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Corner Radial Coordinate $m$')
        #VertCor = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = r'Corner Vertical Coordinate $m$')

        PolLbl = ['XXLoc', 'Theta', 'DJXA', 'Flux Expansion']
        PolVec = xr.DataArray(np.zeros((YSurf,XGrid,N,4)), coords=[Y,X,Attempts,PolLbl], dims=['Radial_Location','Poloidal_Location','Attempt','Poloidal Metric'], name = 'Poloidal Coordinate Data')            

        for n in range(N):
            Attempt = Attempts[n]
            DRT = '{}/Attempt{}/Output'.format(BASEDRT, str(Attempt))     # MAINPATH DIRECTORY STRING
            #DRT2 = 'SOLPS_2D_prof/Shot0' + Shot + '/Attempt' + str(Attempt) + '/Output2'     #Generate Mesh path
                        
            YYLoc.values[:,:,n] = Yy
            RadLoc.values[:,:,n] = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
            VertLoc.values[:,:,n] = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
            
            for j in range(len(Y)):
                for i in range(len(X)):
                    PsinLoc.values[j,i,n] = GF.psiN(RadLoc.loc[Y[j],X[i],Attempt].values,VertLoc.loc[Y[j],X[i],Attempt].values,)            
            
            PolVec.loc[:,:,Attempt,'XXLoc'] = Xx
            
            YVector=np.zeros((len(X),2))
            YVector[:,0] = RadLoc.values[1,:,n] - RadLoc.values[0,:,n]
            YVector[:,1] = VertLoc.values[1,:,n] - VertLoc.values[0,:,n]            
            
            for i in range(len(X)):
                PolVec.loc[:,X[i],Attempt,'Theta'] = np.degrees(np.math.atan2(np.linalg.det([YVector[JXA,:],YVector[i,:]]),np.dot(YVector[JXA,:],YVector[i,:])))
                if PolVec.loc[:,X[i],Attempt,'Theta'].values[0] < 0 and X[i] < JXA:
                    PolVec.loc[:,X[i],Attempt,'Theta'] = PolVec.loc[:,X[i],Attempt,'Theta'] + 360
            #try:
            #    RadCor.values[:,:,n] = np.loadtxt(DRT2 + '/Rad0Cor' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
            #    VertCor.values[:,:,n] = np.loadtxt(DRT2 + '/Vert0Cor' + str(Attempt),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
            #except:
            #    print("Warning, Grid Corner Coordinates Not Found for Attempt" + str(Attempt))
        
        # Condensed Data Loading Command! :
        
        if AddNew is not None:
            if AddNew not in self.Parameter:
                self.Parameter.append(AddNew)
                P = len(self.Parameter)
        
        for p in range(P):
            if self.Parameter[p] not in self.PARAM.keys(): 
                self.PARAM[self.Parameter[p]] = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,self.Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = self.PARAMDICT[self.Parameter[p]])
                for n in range(N):
                    Attempt = self.Attempts[n]
                    DRT = '{}/Attempt{}'.format(BASEDRT, str(Attempt))   #Generate path
                    try:
                        RawData = np.loadtxt('{}/Output/{}{}'.format(DRT, self.Parameter[p], str(Attempt)),usecols = (3))
                    except Exception as err:
                        print(err)
                        try:
                             RawData = np.loadtxt('{}/Output2/{}{}'.format(DRT, self.Parameter[p], str(Attempt)),usecols = (3))
                        except Exception as err:
                            print(err)
                            print('Parameter {} not found for Attempt {}. Creating NAN Array'.format(self.Parameter[p], str(Attempt)))
                            self.PARAM[self.Parameter[p]].values[:,:,n] = np.nan
                            
                    if RawData.size == 3724:
                        self.PARAM[self.Parameter[p]].values[:,:,n] = RawData.reshape((YDIM,XDIM))[1:YDIM-1,XMin+1:XMax+2]
                    elif RawData.size == 7448:
                        self.PARAM[self.Parameter[p]].values[:,:,n] = RawData.reshape((2*YDIM,XDIM))[1+YDIM:2*YDIM-1,XMin+1:XMax+2]
                        
        if RadSlc == 'all':
            RadSlc = self.PARAM.coords['Radial_Location'].values
        if RadSlc == None:
            RadSlc = [SEP]
            
        if PolSlc == 'all':
            PolSlc = self.PARAM.coords['Poloidal_Location'].values
        if PolSlc == None:
            PolSlc = [JXI,JXA]
                
        if Publish != []:
            plt.rc('font',size=30)
            plt.rc('lines',linewidth=5,markersize=15)
        else:
            plt.rc('font',size=14)
        
        #Save all values into self dictionaries
        
        self.KW['BASEDRT'] = BASEDRT
        self.KW['TOPDRT'] = TOPDRT
        self.VVFILE = np.loadtxt('{}/vvfile.ogr'.format(BASEDRT))
        self.Xx = Xx
        self.Yy = Yy
        self.N = N
        self.P = P
        self.RadCoords['XMin'] = XMin
        self.RadCoords['YYLoc'] = YYLoc        
        self.RadCoords['PsinLoc'] = PsinLoc
        self.RadCoords['PsinAvg'] = PsinAvg
        self.RadCoords['RadLoc'] = RadLoc
        self.RadCoords['VertLoc'] = VertLoc
        self.PolVec = PolVec
        if 'd3d' not in Shot:
                self.RadCoords['RmidAvg'] = RmidAvg
                
    def GetRadCoords(self,RADC, Offset=[0,0]):
        #Separate method to handle radial coordinate switching
        
        #ADAPTIVE SEPARATRIX LOCATOR???

        PsinOffset = Offset[0]
        RadOffset = Offset[1]
        
        if RADC == 'Y':
            RR = self.RadCoords['YYLoc']
            Rexp = None
            Rstr = 'N'        
        elif RADC == 'psin':
            RR = self.RadCoords['PsinLoc']
            Rexp = self.RadCoords['PsinAvg']
            Rexp = Rexp + PsinOffset*np.ones(len(Rexp))
            Rstr = '$\psi_N$'
        elif RADC == 'radial':
            RR = self.RadCoords['RadLoc']
            if self.Shot != 'd3d':
                Rexp = self.RadCoords['RmidAvg']
                Rexp = Rexp + RadOffset*np.ones(len(Rexp))
            else:
                print('No experimental radial coordinates available!')
            Rstr = 'm'
        elif RADC == 'rrsep':
            PsinSep=np.abs(self.RadCoords['PsinLoc']-1)
            Sign = np.sign(self.RadCoords['PsinLoc']-1)
            RADSEP=self.RadCoords['RadLoc'].where(PsinSep==PsinSep.min(axis=0)).values.flatten('F')
            RADSEP=RADSEP[~np.isnan(RADSEP)]
            RRsepRad = self.RadCoords['RadLoc'][:,:,0] - RADSEP
            VERTSEP=self.RadCoords['VertLoc'].where(PsinSep==PsinSep.min(axis=0)).values.flatten('F')
            VERTSEP=VERTSEP[~np.isnan(VERTSEP)]
            RRsepVert = self.RadCoords['VertLoc'][:,:,0] - VERTSEP
            RR = np.sqrt(RRsepRad**2 + RRsepVert**2)*Sign
            Rexp = None
            Rstr = '$R-R_{sep}$ (m)'           
        else:
            print('Invalid Radial Coordinate specified')
            
        return RR, Rexp, Rstr
    
    def VeslMesh(self,Parameter=None,**kwargs):
        figVesl = plt.figure(figsize=(14,10))
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
    
    def Contour(self,Parameter=None,**kwargs):
        if Parameter is None:
            self.PltParams = self.Parameter
        elif isinstance(Parameter,str):
            self.PltParams = [Parameter]
        else:
            self.PltParams = Parameter
            
        #Override Statement for Method Keywords to override over KW Keywords:
            
        for key, value in self.KW.items():
            if key not in kwargs.keys():
                kwargs[key] = value
        
        N = self.N
        Xx = self.Xx
        Yy = self.Yy
        Attempts = self.Attempts
        Shot = self.Shot
        VVFILE = self.VVFILE
        PolVec = self.PolVec
        Markers = kwargs['Markers']
        POLC = kwargs['POLC']
        Publish = kwargs['Publish']
        CoreBound = kwargs['CoreBound']
        DIVREG = kwargs['DIVREG']
        JXA = kwargs['JXA']
        JXI = kwargs['JXI']
        SEP = kwargs['SEP']
        Offset = [kwargs['PsinOffset'],kwargs['RadOffset']]             
        ContKW = kwargs
        
        RadLoc = self.RadCoords['RadLoc']
        VertLoc = self.RadCoords['VertLoc']
        XMin = self.RadCoords['XMin']
        
        RR, Rexp, Rstr = self.GetRadCoords(ContKW['RADC'], Offset)
        
        if POLC == 'theta':
            PP = PolVec.loc[:,:,:,'Theta']
            PolXLbl = r'Poloidal Angle $\theta$' 
        elif POLC == 'X':
            PP = PolVec.loc[:,:,:,'XXLoc']
            PolXLbl = r'Poloidal Cell Index $X$'
        elif POLC == 'djxa':
            PP = PolVec.loc[:,:,:,'DJXA']
            PolXLbl = r'Distance from Outer Midplane $m$'
        
        for pn in self.PltParams:
            try:
                PARAM = self.PARAM[pn].copy()
            except:
                try:
                    print('Plot Parameter {} Not Loaded Into Object. Attempting To Load Data...'.format(pn))
                    self._LoadSOLPSData(AddNew=pn)
                    print('Parameter {} Data Load Successful!'.format(pn))
                    PARAM = self.PARAM[pn].copy()
                except:
                    print('Plot Parameter {} Does Not Exist Or Could Not Be Loaded! Skipping Plot...'.format(pn))
                    pass
            
            #print('Beginning Plot Sequence')
            
            if ContKW['SUBTRACT'] is True:
                CMAP = cm.seismic
                for n in np.arange(1,N):
                    PARAM.values[:,:,n] = (PARAM.values[:,:,0]-PARAM[:,:,n])
                PARAM.values[:,:,0] = PARAM.values[:,:,0]-PARAM[:,:,0]
            else:
                CMAP = cm.viridis
            
            # NEED TO FIX -> PREVENT PARAM FROM BEING OVERWRITTEN EVERY TIME
            
            if ContKW['LOG10'] == 1:
                #PARAM.values[PARAM.values<0] = 0
                PARAM.values[PARAM.values>1] = np.log10(PARAM.values[PARAM.values>1])
                PARAM.values[PARAM.values<-1] = -1*np.log10(np.abs(PARAM.values[PARAM.values<-1]))    
                y_exp = np.arange(np.floor(np.nanmin(PARAM.values)), np.ceil(np.nanmax(PARAM.values))+1,2)
                levs = np.arange(np.floor(PARAM.values.min()),np.ceil(PARAM.values.max()))
            elif ContKW['LOG10'] == 2:
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
            else:
                levs = np.linspace(np.floor(PARAM.values.min()),np.ceil(PARAM.values.max()),ContKW['LVN'])
                
            for n in range(N):
                if ContKW['AX'] is None:
                    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(14,10))
                else:
                    ax = ContKW['AX']
                
                if ContKW['GEO'] is True:
                    if ContKW['LOG10'] == 2:
                        if DIVREG is True:
                            IM2 = ax.contourf(RadLoc.loc[:,1:CoreBound[0]-1,Attempts[n]],VertLoc.loc[:,1:CoreBound[0]-1,Attempts[n]],PARAM.loc[:,1:CoreBound[0]-1,Attempts[n]],levs,cmap=CMAP,norm=colors.LogNorm())
                            IM3 = ax.contourf(RadLoc.loc[:,CoreBound[1]+1:,Attempts[n]],VertLoc.loc[:,CoreBound[1]+1:,Attempts[n]],PARAM.loc[:,CoreBound[1]+1:,Attempts[n]],levs,cmap=CMAP,norm=colors.LogNorm())

                        IM1 = ax.contourf(RadLoc.loc[:,CoreBound[0]:CoreBound[1],Attempts[n]],VertLoc.loc[:,CoreBound[0]:CoreBound[1],Attempts[n]],PARAM.loc[:,CoreBound[0]:CoreBound[1],Attempts[n]],levs,cmap=CMAP,norm=colors.LogNorm())
                                            
                    else:
                        if DIVREG is True:
                            IM2 = ax.contourf(RadLoc.loc[:,1:CoreBound[0]-1,Attempts[n]],VertLoc.loc[:,1:CoreBound[0]-1,Attempts[n]],PARAM.loc[:,1:CoreBound[0]-1,Attempts[n]],levs,cmap=CMAP)
                            IM3 = ax.contourf(RadLoc.loc[:,CoreBound[1]+1:,Attempts[n]],VertLoc.loc[:,CoreBound[1]+1:,Attempts[n]],PARAM.loc[:,CoreBound[1]+1:,Attempts[n]],levs,cmap=CMAP)

                        IM1 = ax.contourf(RadLoc.loc[:,CoreBound[0]:CoreBound[1],Attempts[n]],VertLoc.loc[:,CoreBound[0]:CoreBound[1],Attempts[n]],PARAM.loc[:,CoreBound[0]:CoreBound[1],Attempts[n]],levs,cmap=CMAP)                      
                    
                    if Markers is True:                 
                        ax.plot(RadLoc.values[:,(JXA-XMin),n],VertLoc.values[:,(JXA-XMin),n],color='Orange',linewidth=3)
                        ax.plot(RadLoc.values[:,(JXI-XMin),n],VertLoc.values[:,(JXI-XMin),n],color='Red',linewidth=3)
                    
                    ax.plot(RadLoc.values[SEP,:,n],VertLoc.values[SEP,:,n],color='Black',linewidth=3)
                    ax.plot(VVFILE[:,0]/1000,VVFILE[:,1]/1000)
                    ax.set_xlabel('Radial Location (m)')
                    ax.set_ylabel('Vertical Location (m)')
                    ax.set_aspect('equal')
                else:
                    if ContKW['LOG10'] == 2:
                        if DIVREG is True and POLC != 'theta':
                            IM2 = ax.contourf(PP.loc[:,1:CoreBound[0]-1,Attempts[n]],RR.loc[:,1:CoreBound[0]-1,Attempts[n]],PARAM.loc[:,1:CoreBound[0]-1,Attempts[n]],levs,cmap=CMAP,norm=colors.LogNorm())
                            IM3 = ax.contourf(PP.loc[:,CoreBound[1]+1:,Attempts[n]],RR.loc[:,CoreBound[1]+1:,Attempts[n]],PARAM.loc[:,CoreBound[1]+1:,Attempts[n]],levs,cmap=CMAP,norm=colors.LogNorm())
                        
                        IM1 = ax.contourf(PP.loc[:,CoreBound[0]:CoreBound[1],Attempts[n]],RR.loc[:,CoreBound[0]:CoreBound[1],Attempts[n]],PARAM.loc[:,CoreBound[0]:CoreBound[1],Attempts[n]],levs,norm=colors.LogNorm(),cmap=CMAP)
                    else:
                        if DIVREG is True and POLC != 'theta':
                            IM2 = ax.contourf(PP.loc[:,1:CoreBound[0]-1,Attempts[n]],RR.loc[:,1:CoreBound[0]-1,Attempts[n]],PARAM.loc[:,1:CoreBound[0]-1,Attempts[n]],levs,cmap=CMAP)
                            IM3 = ax.contourf(PP.loc[:,CoreBound[1]+1:,Attempts[n]],RR.loc[:,CoreBound[1]+1:,Attempts[n]],PARAM.loc[:,CoreBound[1]+1:,Attempts[n]],levs,cmap=CMAP)

                        IM1 = ax.contour(PP.loc[:,CoreBound[0]:CoreBound[1],Attempts[n]],RR.loc[:,CoreBound[0]:CoreBound[1],Attempts[n]],PARAM.loc[:,CoreBound[0]:CoreBound[1],Attempts[n]],levs,cmap=CMAP)
                    
                    if Markers is True:
                        ax.plot(PP.loc[:,JXA,Attempts[n]],RR.loc[:,JXA,Attempts[n]],color='Orange',linewidth=3)                         #Outer Midplane
                        ax.plot(PP.loc[:,JXI,Attempts[n]],RR.loc[:,JXI,Attempts[n]],color='Red',linewidth=3)                            #Inner Midplane
                    
                    ax.plot(PP.loc[1:SEP,CoreBound[0],Attempts[n]],RR.loc[1:SEP,CoreBound[0],Attempts[n]],color='Black',linewidth=3)    #Inner PFR Boundary
                    ax.plot(PP.loc[1:SEP,CoreBound[1],Attempts[n]],RR.loc[1:SEP,CoreBound[1],Attempts[n]],color='Black',linewidth=3)    #Outer PFR Boundary
                    ax.plot(PP.loc[SEP,:,Attempts[n]],RR.loc[SEP,:,Attempts[n]],color='Black',linewidth=3)                              #Separatrix
                    ax.set_xlabel(PolXLbl)
                    ax.set_ylabel(Rstr) 

                if ContKW['Publish'] != []:
                    ax.set_title('Attempt {} {}'.format(Publish[n], PARAM.name))
                else:
                    ax.set_title('Discharge 0{} Attempt {} {}'.format(Shot, str(Attempts[n]), PARAM.name))
                    
                plt.colorbar(IM1,ax=ax)

                #a.set_xticklabels(['%.1f' % i for i in a.get_xticks()], fontsize='x-large')
                #a.set_yticklabels(['%.1f' % j for j in a.get_yticks()], fontsize='x-large')
                #a.tick_params(labelsize=20)
                
                if ContKW['GRID'] is True:
                    plt.grid()
                
                if ContKW['SAVE'] is True:
                    ImgName = 'Profiles/{}{}Contour.png'.format(Parameter, str(Attempts[n]))
                    plt.savefig(ImgName, bbox_inches='tight')

            self.ContKW = ContKW

    def PolPlot(self,Parameter=None,**kwargs):
        #with plt.xkcd():
            
            for key, value in self.KW.items():
                if key not in kwargs.keys():
                    kwargs[key] = value

            Shot = self.Shot
            Attempts = self.Attempts
            PolVec = self.PolVec
            Publish = kwargs['Publish']
            PlotScheme = kwargs['PlotScheme']
            Markers = kwargs['Markers']
            JXA = kwargs['JXA']
            JXI = kwargs['JXI']
            SEP = kwargs['SEP']
            POLC = kwargs['POLC']
            CoreBound = kwargs['CoreBound']
            CoreBound[1] = CoreBound[1]-1
            PolKW = kwargs
            
            SURF = 17
            
            if Parameter is None:
                self.PltParams = self.Parameter
            elif isinstance(Parameter,str):
                self.PltParams = [Parameter]
            else:
                self.PltParams = Parameter
                
            for pn in self.PltParams:
                try:
                    PARAM = self.PARAM[pn].copy()
                except:
                    try:
                        print('Plot Parameter Not Loaded Into Object. Attempting To Load Data...')
                        self._LoadSOLPSData(AddNew=pn)
                        print('Parameter Data Load Successful!')
                        PARAM = self.PARAM[pn].copy()
                    except:
                        print('Plot Parameter Does Not Exist Or Could Not Be Loaded! Skipping Plot...')
                        pass

            if POLC == 'theta':
                PP = PolVec.loc[:,:,:,'Theta']
                PolXLbl = r'Poloidal Angle $\theta$' 
            elif POLC == 'X':
                PP = PolVec.loc[:,:,:,'XXLoc']
                PolXLbl = r'Poloidal Cell Index $X$'
            elif POLC == 'djxa':
                PP = PolVec.loc[:,:,:,'DJXA']
                PolXLbl = r'Distance from Outer Midplane $m$'

            #print('Beginning Plot Sequence')                

            if PolKW['AX'] is None: # if no axis was supplied to the function create our own
                fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(14,7))
            else:
                ax = PolKW['AX']
            
            if PolKW['GRAD'] is True:
                PARAM.values = np.gradient(PARAM.values,axis=0)
                
            if PolKW['LOG10'] == 1:
                #PARAM.values[PARAM.values<0] = 0
                PARAM.values[PARAM.values>1] = np.log10(PARAM.values[PARAM.values>1])
                PARAM.values[PARAM.values<-1] = -1*np.log10(np.abs(PARAM.values[PARAM.values<-1]))      
            
            if len(PlotScheme) == self.N:
                for n in range(self.N):
                    if PolKW['LOG10'] == 2:
                        ax.semilogy(PP.loc[SURF,CoreBound[0]:CoreBound[1],Attempts[n]], PARAM.loc[SURF,CoreBound[0]:CoreBound[1],Attempts[n]], PlotScheme[n])
                    else:
                        ax.plot(PP.loc[SURF,CoreBound[0]:CoreBound[1],Attempts[n]], PARAM.loc[SURF,CoreBound[0]:CoreBound[1],Attempts[n]], PlotScheme[n], linewidth=3)
            else:            
                if PolKW['LOG10'] == 2:
                    ax.semilogy(PP.loc[SURF,CoreBound[0]:CoreBound[1],:], PARAM.loc[SURF,CoreBound[0]:CoreBound[1],:])
                else:
                    ax.plot(PP.loc[SURF,CoreBound[0]:CoreBound[1],:], PARAM.loc[SURF,CoreBound[0]:CoreBound[1],:],linewidth=3)

            if Publish==[]:
                plt.legend(Attempts)
                plt.title('Discharge 0' + str(Shot) + ' Attempt(s) ' + str(Attempts) + ' Poloidal ' + PARAM.name)
            else:
                plt.legend(Publish)
                plt.title('Poloidal ' + PARAM.name + ' along Separatrix')
            Pmin = float(PARAM.loc[SURF,CoreBound[0]:CoreBound[1],:].min())
            Pmax = float(PARAM.loc[SURF,CoreBound[0]:CoreBound[1],:].max())    
            
            if Markers == True:
                plt.plot([PP.loc[SURF,JXI,Attempts[0]], PP.loc[SURF,JXI,Attempts[0]]],[Pmin, Pmax],'r--')
                plt.plot([PP.loc[SURF,JXA,Attempts[0]], PP.loc[SURF,JXA,Attempts[0]]],[Pmin, Pmax],'--')
            
            plt.xlabel(PolXLbl)
            plt.ylabel(PARAM.name)
            plt.grid()
    
    def RadProf(self,Parameter=None,**kwargs):
        
        for key, value in self.KW.items():
            if key not in kwargs.keys():
                kwargs[key] = value
        
        Shot = self.Shot
        Attempts = self.Attempts
        Publish = kwargs['Publish']
        JXA = kwargs['JXA']
        SEP = kwargs['SEP']
        RADC = kwargs['RADC']
        PolSlc = kwargs['PolSlc']
        Markers = kwargs['Markers']
        PlotScheme = kwargs['PlotScheme']
        Offset = [kwargs['PsinOffset'],kwargs['RadOffset']]         
        RadProfKW = kwargs
        
        RR, Rexp, Rstr = self.GetRadCoords(RADC,Offset)
        
        if Parameter is None:
            self.PltParams = self.Parameter
        elif isinstance(Parameter,str):
            self.PltParams = [Parameter]
        else:
            self.PltParams = Parameter
            
        for pn in self.PltParams:
            try:
                PARAM = self.PARAM[pn].copy()
            except:
                try:
                    print('Plot Parameter Not Loaded Into Object. Attempting To Load Data...')
                    self._LoadSOLPSData(AddNew=pn)
                    print('Parameter Data Load Successful!')
                    PARAM = self.PARAM[pn].copy()
                except:
                    print('Plot Parameter Does Not Exist Or Could Not Be Loaded! Skipping Plot...')
                    pass

            #print('Beginning Plot Sequence')
            
            if RadProfKW['AX'] is None: # if no axis was supplied to the function create our own
                fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(14,7))
            else:
                ax = RadProfKW['AX']
            
            if RadProfKW['GRAD'] is True:
                PARAM.values = np.gradient(PARAM.values,axis=0)
                
            if RadProfKW['LOG10'] == 1:
                #PARAM.values[PARAM.values<0] = 0
                PARAM.values[PARAM.values>1] = np.log10(PARAM.values[PARAM.values>1])
                PARAM.values[PARAM.values<-1] = -1*np.log10(np.abs(PARAM.values[PARAM.values<-1]))    
                y_exp = np.arange(np.floor(np.nanmin(PARAM.values)), np.ceil(np.nanmax(PARAM.values))+1,2)    
            
            
            if len(PlotScheme) == self.N:
                for n in range(self.N):
                    if RadProfKW['LOG10'] == 2:
                        ax.semilogy(RR.loc[:,JXA,Attempts[n]], PARAM.loc[:,JXA,Attempts[n]],PlotScheme[n])
                    else:
                        ax.plot(RR.loc[:,JXA,Attempts[n]], PARAM.loc[:,JXA,Attempts[n]],PlotScheme[n],linewidth=3)
                        
                    if RadProfKW['EXP'] is True and RadProfKW['RADC'] != 'rrsep' and RadProfKW['RADC'] != 'Y':
                        if 'd3d' not in Shot:
                            if pn == 'Ne':
                                NemidAvg = self.ExpDict['NemidAvg']
                                ErrNe = self.ExpDict['ErrNe']
                                ax.errorbar(Rexp,NemidAvg,yerr=ErrNe,fmt='o',markersize=7,linewidth=3,capsize=7,color=PlotScheme[n][0])
                            elif pn == 'Te':
                                TemidAvg = self.ExpDict['TemidAvg']
                                ErrTe = self.ExpDict['ErrTe']
                                ax.errorbar(Rexp,TemidAvg,yerr=ErrTe,fmt='o',markersize=7,linewidth=3,capsize=7,color=PlotScheme[n][0])
                        
                        if 'd3d' in Shot:
                            if pn == 'Ne':
                                PsinNe = Rexp[0]
                                Ned3d = self.ExpDict['Ned3d']
                                ax.plot(PsinNe,Ned3d,'o',markersize=7,linewidth=3,color=PlotScheme[n][0])
                            elif pn == 'Te':
                                PsinTe = Rexp[1]
                                Ted3d = self.ExpDict['Ted3d']
                                ax.plot(PsinTe,Ted3d,'o',markersize=7,linewidth=3,color=PlotScheme[n][0])
                            elif pn == 'Ti':
                                PsinTi = Rexp[2]
                                Tid3d = self.ExpDict['Tid3d']
                                ax.plot(PsinTi,Tid3d,'o',markersize=7,linewidth=3,color=PlotScheme[n][0])
            else:
                if RadProfKW['LOG10'] == 2:
                    ax.semilogy(RR.loc[:,JXA,:], PARAM.loc[:,JXA,:])
                else:
                    ax.plot(RR.loc[:,JXA,:], PARAM.loc[:,JXA,:],linewidth=3)

                if RadProfKW['EXP'] is True and RadProfKW['RADC'] != 'rrsep' and RadProfKW['RADC'] != 'Y':
                    if 'd3d' not in Shot:
                        if pn == 'Ne':
                            NemidAvg = self.ExpDict['NemidAvg']
                            ErrNe = self.ExpDict['ErrNe']
                            ax.errorbar(Rexp,NemidAvg,yerr=ErrNe,fmt='o',markersize=7,linewidth=3,capsize=7)
                        elif pn == 'Te':
                            TemidAvg = self.ExpDict['TemidAvg']
                            ErrTe = self.ExpDict['ErrTe']
                            ax.errorbar(Rexp,TemidAvg,yerr=ErrTe,fmt='o',markersize=7,linewidth=3,capsize=7)
                    
                    if 'd3d' in Shot:
                        if pn == 'Ne':
                            PsinNe = Rexp[0]
                            Ned3d = self.ExpDict['Ned3d']
                            ax.plot(PsinNe,Ned3d,'o',markersize=7,linewidth=3)
                        elif pn == 'Te':
                            PsinTe = Rexp[1]
                            Ted3d = self.ExpDict['Ted3d']
                            ax.plot(PsinTe,Ted3d,'o',markersize=7,linewidth=3)
                        elif pn == 'Ti':
                            PsinTi = Rexp[2]
                            Tid3d = self.ExpDict['Tid3d']
                            ax.plot(PsinTi,Tid3d,'o',markersize=7,linewidth=3)
            
            ax.set_xlabel(Rstr)
            
            if RadProfKW['LOG10'] == 1:
                ax.set_ylabel('Log_10 of {}'.format(PARAM.name))
            else:    
                ax.set_ylabel(PARAM.name)
            if RadProfKW['Publish']==[]:
                ax.legend(Attempts)
                ax.set_title('Discharge 0{} Attempt(s) {} Midplane Radial {}'.format(str(Shot), str(Attempts), PARAM.name))
            else:
                ax.legend(Publish)
                ax.set_title('Radial Midplane {}'.format(PARAM.name))
            Pmin = float(PARAM.loc[:,JXA,:].min())
            Pmax = float(PARAM.loc[:,JXA,:].max())
            if Markers == True:
                ax.plot([RR.loc[SEP,JXA,Attempts[0]], RR.loc[SEP,JXA,Attempts[0]]],[Pmin, Pmax],color='Black',linewidth=3)
            
            if RadProfKW['GRID'] is True:
                ax.grid(b=1)
                
            self.RadProfKW = RadProfKW


### FUNCTIONS BELOW NEED TO BE WORKED ON


    def RadPlot(self,Parameter=None,**kwargs):
        
            for key, value in self.KW.items():
                if key not in kwargs.keys():
                    kwargs[key] = value
                
            RadKW = kwargs    
        
            fig2b = plt.figure(figsize=(14,10))
            if GRAD == 1:
                PARAM.values = np.gradient(PARAM.values,axis=0)
            if LOG10 == 2:
                plt.semilogy(RR.loc[:,PolSlc,Attempts[0]].values, PARAM.loc[:,PolSlc,Attempts[0]].values)
            else:
                plt.plot(RR.loc[:,PolSlc,Attempts[0]].values, PARAM.loc[:,PolSlc,Attempts[0]].values)
            if PolSlc==[JXI,JXA]:
                plt.legend(['Inner Midplane','Outer Midplane'])
            else:
                plt.legend(PolSlc)
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

    def Surface(self,Parameter=None,**kwargs):
        
        for key, value in self.KW.items():
            if key not in kwargs.keys():
                kwargs[key] = value
                
        SurfKW = kwargs
        
        if LOG10 == 2:
            PARAM.values[PARAM.values==0] = np.nan
            PARAM.values = np.log10(PARAM.values)
        Zmax = np.nanmax(PARAM.values) 
        for n in range(N):
            fig1 = plt.figure(figsize=(18,12))  # Size is width x height
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
            if SAVE == 1:
                ImgName = 'Profiles/' + Parameter + 'SurfA' + str(Attempts[n]) + '.png'
                plt.savefig(ImgName, bbox_inches='tight')
    
    
    def SumPolPlot(self,Parameter=None,**kwargs):
        #with plt.xkcd():
            
            for key, value in self.KW.items():
                if key not in kwargs.keys():
                    kwargs[key] = value
                
            SumPolKW = kwargs
            
            fig4 = plt.figure(figsize=(14,7))
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
    
    #plt.show()
    
    def Export(self):
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
