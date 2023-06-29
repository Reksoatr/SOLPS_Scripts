# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 14:00:00 2019

@author: Rmreksoatmodjo

SOLPS_Plotter V1.0 - Renamed from VesselPlotterNew
Last Revised 6/27/2023
"""
import os
import numpy as np
import xarray as xr
import glob
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import colors, cm
from scipy.io import loadmat
from scipy.interpolate import interp1d, griddata, RectBivariateSpline
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
        (Specified by a string - Defaults are:  
        'Ne' - Electron Density
        'Te' - Electron Temperature
        'Ti' - Ion Temperature
        'DN' - Particle Diffusion Coefficient
        'KYE' - Electron Thermal Diffusion Coefficient
        'KYI' - Ion Thermal Diffusion Coefficient  
        'NeuDen' - Neutral Density
        'IonFlx' - Radial Ion Particle Flux
        'IonPol' - Poloidal Ion Particle Flux
        'RadPinch' - Radial Pinch Velocity)
            
    OPTIONAL:

    TimeRange = [1.10,1.30] > Time range (in sec) over which experimental data is averaged    
    DEV = 'cmod' > Tokamak device being investigated
    EXP = True > Plot Experimental Data Points
    LOG10 = 0 > Plot Base-10 Logarithm of Parameter Data (0, 1, or 2)
        0 : No log plot or log processing of data
        1 : Take base-10 log of data, then plot on linear y-scale (ticks will be 0, 1, 2,...)
        2 : Plot data as is, then change y-axis scale to logarithmic (ticks will be 0, 10, 100,...)
    GRAD = False > Calculate Gradient of Parameter Data
    ELEV = 75 > Elevation of viewing angle
    AZIM = 270 > Azimuthal viewing angle for Surface Plot
    JXI = 37 > Poloidal position of Inner Midplane - default is 37
    JXA = 55 > Poloidal position of Outer Midplane - default is 55
    SEP = 20 > Radial position of Separatrix - default is 20    
    XDIM = 98 > Dimension of computational grid in the x (poloidal) direction
    YDIM = 38 > Dimension of computational grid in the y (radial) direction
    CoreBound = [24,71] > X-coordinate grid cell numbers that define the [left, right] bounds of Core Region
    CMAP = cm.viridis > Colormap used for 2D Contour Maps
    Colors = None > Specific colors to color each contour level of a contour map
    Publish = [] > List of strings to use in legends of publication-quality plots; if not [], changes plotting rc.params 
    Markers = True > Plot separatrix, outer midplane, and/or inner midplane marker line
    PlotScheme = [] > List of matplotlib plot style flags for line plots; must either be empty or a list of strings as long as the number of Attempts    
    PsinOffset = 0 > Psi_n offset of Experimental Data Points
    RadOffset = 0 > Radial (in meters) offset of Experimental Data Points
    RADC = 'psin' > Set radial coordinate convention - Either 'psin', 'radial', 'rrsep' or 'Y'    
    POLC = 'dXP' > Set poloidal coordinate convention - Either 'X', 'theta', 'dXP', or 'dXP_norm'
    RadSlc = None > Radial surface selection for poloidal plots - Can set specific radial index, 'all', or 'None' defaults to SEP
    PolSlc = None > Poloidal grid line selection for radial plots - Can set specific poloidal index, 'all', or 'None' defaults to [JXA, JXI]
    SURF = 20 > Same purpose as PolSlc
    PHYS = True > Map Contour to PHYSical Geometry; if False, plots on rectangular grid     
    LVN = 100 > Number of colorbar levels for contour plots or List [] of specific contour levels
    LINTHRESH = 1.0 > +/- Value at which divergent-logarithmic Contour plots should switch to linear scale going from positive to negative values
    DIVREG = True > Include Divertor Region Data (May cause loss of logarithmic resolution)   
    SAVE = False > Save plot to a .png file 
    SUBTRACT = False > For a series of Attempts, subtract the matrix data of each Attempt from the first declared Attempt, to compare difference 
    AVG = False > Take Average of each parameter dataset over all supplied Attempts, append averaged data to end of matrix dataset
    TC_Flux = [] > For d3d cases only; allows specification of ONETWO-calculated Radial Flux values
    TC_Psin = [] > For d3d cases only; allows specification of Psin coordinates of ONETWO-calculated Radial Flux values
    GRID = False > Turn on plot grid
    AX = None > Pass the name of a matplotlib axis for the SOLPSPLOT object to plot on; by default SOLPSPLOT plots on a new axis
    BASEDRT = 'SOLPS_2D_prof/' > Local base directory for all stored SOLPS simulation run data
    TOPDRT = '' > Local home directory, parent of BASEDRT and 'gfileProcessing' directories
    ROOTSHOT = '1160718' > Optional numerical prefix to designate the experimental campaign that a series of shots may share in common 
        
    PLOT TYPES:
        
    SOLPSPLOT.Contour() > Plot 2D spatial Contour Plot
    ! SOLPSPLOT.Surface() > Plot 3D spatial Suface Plot
    SOLPSPLOT.PolPlot() > Plot Poloidal Contour along Separatrix (SEP) 
    ! SOLPSPLOT.RadPlot() > Plot Radial Contours for every Poloidal grid point (ONLY TAKES 1 ATTEMPT INPUT)
    SOLPSPLOT.RadProf() > Plot Radial Profile of Parameter at Outer Midplane (JXA)
    ! SOLPSPLOT.SumPlot() > Plot Poloidally-Summed Parameter values vs Radial location
    SOLPSPLOT.B2Mesh() > Plot SOLPS B2.5 Curvilinear Grid
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
            self.Parameter = ['Ne','Te','Ti','DN','KYE','KYI','NeuDen','IonFlx','IonPol','RadPinch']
            
        self.DefaultSettings = {'TimeRange' : [1.10,1.30],  
                     'DEV': 'cmod',
                     'EXP' : True,
                     'LOG10' : 0,
                     'GRAD' : False,
                     'ELEV' : 75,
                     'AZIM' : 270,
                     'JXI' : 37,
                     'JXA' : 55,
                     'SEP' : 20,
                     'XDIM' : 98,
                     'YDIM' : 38,
                     'CoreBound' : [24,71],
                     'CMAP' : cm.viridis,
                     'Colors' : None,
                     'Publish' : [],
                     'Markers' : True,
                     'PlotScheme' : [],
                     'PsinOffset' : 0,
                     'RadOffset' : 0,
                     'RADC' : 'psin',
                     'POLC' : 'dXP',
                     'RadSlc' : None,
                     'PolSlc' : None,
                     'SURF' : 20,
                     'PHYS' : True,
                     'LVN' : 100,
                     'LINTHRESH' : 1.0,
                     'DIVREG' : True,
                     'SAVE' : False,
                     'SUBTRACT' : False,
                     'AVG' : False,
                     'TC_Flux' : [], 
                     'TC_Psin' : [],
                     'GRID' : False,
                     'MESHGRID' : None,
                     'AX' : None,
                     'BASEDRT' : 'solps-iter/runs/',
                     'TOPDRT' : '',
                     'ROOTSHOT' : ''} #1160718
        
        for key, value in self.DefaultSettings.items():
            if key not in kwargs.keys():
                kwargs[key] = value
                    
        self.KW = kwargs         
        
        self.PARAMDICT = {'Ne': r'Electron Density $n_e\;(m^{-3})$',
                     'Te': r'Electron Temperature $T_e\;(eV)$',
                     'NI' : r'Ion (D+) Density $n_i\;(m^{-3})$',
                     'Ti': r'Ion (D+) Temperature $T_i\;(eV)$',
                     'DN': r'Particle Density Diffusivity $D\;(m^2/s)$',
                     'KYE': r'Electron Thermal Diffusivity $\chi_e (m^2/s)$',
                     'KYI': r'Ion Thermal Diffusivity $\chi_i (m^2/s)$',
                     'NeuDen': r'Neutral Atom (D) Density $(m^{-3})$',
                     'MolDen': r'Neutral Molecule (D2) Density $(m^{-3})$',
                     'NeuTemp': r'Neutral Atom (D) Temperature (eV)',
                     'MolTemp': r'Neutral Molecule (D2) Temperature (eV)',
                     'IonFlx': r'Radial Ion (D+) Flux $s^{-1}$',
                     'IonPol': r'Poloidal Ion (D+) Flux $s^{-1}$',
                     'MolFlx': r'Radial Molecule Flux $m^{-2}s^{-1}$',
                     'RadFlx': r'Radial Atomic Flux $m^{-2}s^{-1}$',
                     'VLY': r'Radial Pinch $v_y (m^2/s)$',
                     'SX': r'Poloidal Contact Area SX $(m^{2})$',
                     'SY': r'Radial Contact Area SY $(m^{2})$',
                     'SZ': r'Poloidal Cross-Sectional Area SZ $(m^{2})$',
                     'VOL': r'Cell Volume VOL $(m^{3})$',
                     'RadPinch': r'Anomalous Radial Pinch Velocity $v_y\;(m/s)$',
                     'AHAL': r'Atomic $H_\alpha$ emissivity $(photons*m^{-3}s^{-1})$',
                     'MHAL': r'Molecular $H_\alpha$ emissivity $(photons*m^{-3}s^{-1})$',
                     'NTI' : r'Test Ion (D2+) Density $n_{ti}\;(m^{-3})$',
                     'TestIonTemp' : r'Test Ion (D2+) Temperature $T_{testion}\;(eV)$',
                     'LyaEmiss' : r'Lyman-alpha Emissivity $(photons*m^{-3}s^{-1}sr^{-1})$',
                     'LyaEmissW' : r'Converted Lyman-alpha Emissivity $(W*m^{-3})$',
                     'LyaEmissW_IE' : r'Ion-Electron Reaction Lyman-alpha Emissivity $(W*m^{-3})$'}
        
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
        AVG = self.KW['AVG']
        JXA = self.KW['JXA']
        JXI = self.KW['JXI']
        EXP = self.KW['EXP']
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
        ROOTSHOT = self.KW['ROOTSHOT']
        
        BASEDRT, TOPDRT = SET_WDIR(BASEDRT,TOPDRT)
        
        #print(BASEDRT)
        #print(TOPDRT)
        
        Attempts = [str(i) for i in Attempts]
        print('Attempts {} Requested...'.format(Attempts))
        # Create Experiment Data Dictionary (ExpDict) -> d3d or cmod?
        
        if 'gas' in Shot:
            BASEDRT = '{}gaspuff/'.format(BASEDRT)
        
        if 'd3d' in Shot or DEV=='d3d':
            
            BASEDRT = '{}d3d'.format(BASEDRT)
            
            JXI = 39
            JXA = 55
            SEP = 17
            
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
            
            if ROOTSHOT == '':
                
                BASEDRT = '{}cmod/{}home'.format(BASEDRT,Shot)
                GFILE = glob.glob('{}gfileProcessing/cmod_files/g{}*'.format(TOPDRT,Shot))
                print(GFILE)
                GF = eq.equilibrium(gfile=GFILE[-1])
                if EXP:
                    ExpFile = Shot
                    ExpData = loadmat('{}gfileProcessing/cmod_files/{}.mat'.format(TOPDRT, ExpFile))
            
            else: 
                
                BASEDRT = '{}cmod/0{}home'.format(BASEDRT, Shot[-2:])
                
                #GFILE = '{}gfileProcessing/cmod_files/g11607180{}.01209_974'.format(TOPDRT, Shot[-2:])
                GFILE = glob.glob('{}gfileProcessing/cmod_files/g{}*{}*'.format(TOPDRT,ROOTSHOT,Shot[-2:]))
                print(GFILE)
                GF = eq.equilibrium(gfile=GFILE[0])
                
                #ExpFile = '11607180{}'.format(Shot[-2:])
                ExpFile = glob.glob('{}gfileProcessing/cmod_files/{}*{}.mat'.format(TOPDRT,ROOTSHOT,Shot[-2:]))
                #print(ExpFile)
                ExpData = loadmat(ExpFile[0])    
                
            if EXP:
                TR=[ii for ii, tt in enumerate(ExpData['time'][0,:]) if tt > TimeRange[0] and tt < TimeRange[1]]
                
                Rmid = ExpData['rmid'][:,TR]
                RmidAvg = np.mean(Rmid, axis=1)
    
                try:
                    Psin = ExpData['psin'][:,TR]
                    PsinAvg = np.mean(Psin, axis=1)
                except:
                    print('Psin coordinates not found. Attempting to approximate experimental psin from gfile')
                    PsinAvg = GF.psiN(RmidAvg,np.zeros(RmidAvg.shape))[0]
                
                #Robust Statistics -> Use MEDIAN, not MEAN, of TS Data -> Need to filter out zeros before taking Median!
                
                Nemid = ExpData['ne'][:,TR]
                NemidAvg=np.zeros(Nemid.shape[0])
                for ii in range(Nemid.shape[0]):
                    if np.sum(Nemid[ii,:]) == 0:
                        continue
                    else:
                        NemidAvg[ii]=np.median(Nemid[ii][np.nonzero(Nemid[ii])])
                        
                ErrNe=np.median(ExpData['nerr'][:,TR], axis=1)
                '''
                Nemid[Nemid == 0] = np.nan
                NemidAvg = np.nanmean(Nemid, axis=1)
                ErrNe = np.nanmean(ExpData['nerr'][:,ti:tf], axis=1)
                '''
                NeThresh = (ErrNe*2)/NemidAvg
    
                for NT in range(len(NeThresh)):
                    if np.abs(NeThresh[NT]) > 5.0:
                        NemidAvg[NT] = np.nan
                        ErrNe[NT] = np.nan
                
                Temid = ExpData['te'][:,TR]
                TemidAvg=np.zeros(Temid.shape[0])
                for ii in range(Temid.shape[0]):
                    if np.sum(Temid[ii,:]) == 0:
                        continue
                    else:
                        TemidAvg[ii]=np.median(Temid[ii][np.nonzero(Temid[ii])])
                ErrTe=np.median(ExpData['terr'][:,TR], axis=1)
                '''
                Temid[Temid == 0] = np.nan
                TemidAvg = np.nanmean(Temid, axis=1)
                ErrTe = np.nanmean(ExpData['terr'][:,ti:tf], axis=1)
                '''
                TeThresh = (ErrTe*2)/TemidAvg
                
                for TT in range(len(TeThresh)):
                    if np.abs(TeThresh[TT]) > 5.0:
                        TemidAvg[TT] = np.nan
                        ErrTe[TT] = np.nan
                
                
                self.ExpDict['NemidAvg'] = NemidAvg
                self.ExpDict['ErrNe'] = ErrNe
                self.ExpDict['TemidAvg'] = TemidAvg
                self.ExpDict['ErrTe'] = ErrTe        
        
        else:
            
            BASEDRT = '{}/{}/{}'.format(BASEDRT,DEV,Shot)
            GFILE = glob.glob('{}/{}/g{}*'.format(TOPDRT,DEV,Shot))
            print(GFILE)
            GF = eq.equilibrium(gfile=GFILE[-1])
            if EXP:
                ExpFile = Shot
                ExpData = loadmat('{}/{}/{}.mat'.format(TOPDRT,DEV,ExpFile))
            
        #Get dimensions of simulation grid 
        
        if AVG is True and 'AVG' not in Attempts:
            Attempts.append('AVG')
            print('Attempts {} Requested...'.format(Attempts))
        
        N = len(Attempts)
        print('{} Attempt(s) Entered'.format(N))
        P = len(self.Parameter)
        
        XGrid=XDIM-2
        XMin=1
        XMax=XGrid
        X_Core=CoreBound[1]-CoreBound[0]+1
        
        YSurf=YDIM-2
        
        X = np.linspace(XMin,XMax,XGrid)
        Y = np.linspace(1,YSurf,YSurf)
        Xx, Yy = np.meshgrid(X,Y)   # Create X and Y mesh grid arrays
        
        RadLoc = xr.DataArray(np.zeros((YSurf,XGrid,N)), 
                              coords=[Y,X,Attempts], 
                              dims=['Radial_Location','Poloidal_Location','Attempt'], 
                              name = r'Radial Coordinate $m$')
        VertLoc = xr.DataArray(np.zeros((YSurf,XGrid,N)), 
                               coords=[Y,X,Attempts], 
                               dims=['Radial_Location','Poloidal_Location','Attempt'], 
                               name = r'Vertical Coordinate $m$')
        
        Core_Corners = xr.DataArray(np.zeros((YSurf+1,X_Core+1,N,2)), 
                                    coords=[np.concatenate(([0],Y)),np.linspace(CoreBound[0],CoreBound[1]+1,X_Core+1),Attempts,['X','Y']], 
                                    dims=['Radial Index','Poloidal Index','Attempt','Point Coordinates'], 
                                    name = r'Core Corner Coordinates $m$')
        Div1_Corners = xr.DataArray(np.zeros((YSurf+1,CoreBound[0]+1,N,2)), 
                                    coords=[np.concatenate(([0],Y)),np.linspace(0,CoreBound[0],CoreBound[0]+1),Attempts,['X','Y']], 
                                    dims=['Radial Index','Poloidal Index','Attempt','Point Coordinates'], 
                                    name = r'Inner Divertor Corner Coordinates $m$')
        Div2_Corners = xr.DataArray(np.zeros((YSurf+1,XMax-CoreBound[1],N,2)), 
                                    coords=[np.concatenate(([0],Y)),np.linspace(CoreBound[1]+1,XMax,XMax-CoreBound[1]),Attempts,['X','Y']], 
                                    dims=['Radial Index','Poloidal Index','Attempt','Point Coordinates'], 
                                    name = r'Outer Divertor Corner Coordinates $m$')

        YYLoc = xr.DataArray(np.zeros((YSurf,XGrid,N)), 
                             coords=[Y,X,Attempts], 
                             dims=['Radial_Location','Poloidal_Location','Attempt'], 
                             name = r'Radial Grid Point $N$')
        PsinLoc = xr.DataArray(np.zeros((YSurf,XGrid,N)), 
                               coords=[Y,X,Attempts], 
                               dims=['Radial_Location','Poloidal_Location','Attempt'], 
                               name = r'Normalized Psi $\psi_N$')
        
        PolLbl = ['XXLoc', 'Theta', 'dXP','dXP_norm']
        PolVec = xr.DataArray(np.zeros((YSurf,XGrid,N,4)), 
                              coords=[Y,X,Attempts,PolLbl], 
                              dims=['Radial_Location','Poloidal_Location','Attempt','Poloidal Metric'], 
                              name = 'Poloidal Coordinate Data')  

        #RadCo = ['YYLoc', 'RRsep', 'PsiNLoc']
        #RadVec = xr.DataArray(np.zeros((YSurf,XGrid,N,3)), 
        #                      coords=[Y,X,Attempts,RadCo], 
        #                      dims=['Radial_Location','Poloidal_Location','Attempt','Radial Metric'], 
        #                      name = 'Radial Coordinate Data')
          

        for n in range(N):
            Attempt = Attempts[n]
            print('Loading data for Attempt {} now...'.format(Attempt))
            if n == N-1 and Attempt == 'AVG':
                YYLoc.values[:,:,n] = Yy
                RadLoc.values[:,:,n] = RadLoc.values[:,:,0:n-1].mean(2)
                VertLoc.values[:,:,n] = VertLoc.values[:,:,0:n-1].mean(2)
                PsinLoc.values[:,:,n] = PsinLoc.values[:,:,0:n-1].mean(2)
                PolVec.values[:,:,n,:] = PolVec.values[:,:,0:n-1,:].mean(2) 
            else:
                DRT = '{}/Attempt{}/Output'.format(BASEDRT, str(Attempt))     # MAINPATH DIRECTORY STRING
                DRT2 = '{}/Attempt{}/Output2'.format(BASEDRT, str(Attempt))    # Mesh data directory path
                            
                YYLoc.values[:,:,n] = Yy
                RadLoc.values[:,:,n] = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
                VertLoc.values[:,:,n] = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
                
                Rad0Cor = np.loadtxt('{}/Rad0Cor{}'.format(DRT2, str(Attempt)),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
                Vert0Cor = np.loadtxt('{}/Vert0Cor{}'.format(DRT2, str(Attempt)),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]

                Rad1Cor = np.loadtxt('{}/Rad1Cor{}'.format(DRT2, str(Attempt)),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
                Vert1Cor = np.loadtxt('{}/Vert1Cor{}'.format(DRT2, str(Attempt)),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]

                Rad2Cor = np.loadtxt('{}/Rad2Cor{}'.format(DRT2, str(Attempt)),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
                Vert2Cor = np.loadtxt('{}/Vert2Cor{}'.format(DRT2, str(Attempt)),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]

                Rad3Cor = np.loadtxt('{}/Rad3Cor{}'.format(DRT2, str(Attempt)),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
                Vert3Cor = np.loadtxt('{}/Vert3Cor{}'.format(DRT2, str(Attempt)),usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
                
                Core_Corners.values[:-1,:-1,n,0] = Rad0Cor[:,CoreBound[0]:CoreBound[1]+1]
                Core_Corners.values[:-1,-1,n,0] = Rad1Cor[:,CoreBound[1]]
                Core_Corners.values[-1,:-1,n,0] = Rad2Cor[-1,CoreBound[0]:CoreBound[1]+1]
                Core_Corners.values[-1,-1,n,0] = Rad3Cor[-1,CoreBound[1]]
                
                Core_Corners.values[:-1,:-1,n,1] = Vert0Cor[:,CoreBound[0]:CoreBound[1]+1]
                Core_Corners.values[:-1,-1,n,1] = Vert1Cor[:,CoreBound[1]]
                Core_Corners.values[-1,:-1,n,1] = Vert2Cor[-1,CoreBound[0]:CoreBound[1]+1]
                Core_Corners.values[-1,-1,n,1] = Vert3Cor[-1,CoreBound[1]]
                
                Div1_Corners.values[:-1,:-1,n,0] = Rad0Cor[:,XMin-1:CoreBound[0]]
                Div1_Corners.values[:-1,-1,n,0] = Rad1Cor[:,CoreBound[0]-1]
                Div1_Corners.values[-1,:-1,n,0] = Rad2Cor[-1,XMin-1:CoreBound[0]]
                Div1_Corners.values[-1,-1,n,0] = Rad3Cor[-1,CoreBound[0]-1]
                
                Div1_Corners.values[:-1,:-1,n,1] = Vert0Cor[:,XMin-1:CoreBound[0]]
                Div1_Corners.values[:-1,-1,n,1] = Vert1Cor[:,CoreBound[0]-1]
                Div1_Corners.values[-1,:-1,n,1] = Vert2Cor[-1,XMin-1:CoreBound[0]]
                Div1_Corners.values[-1,-1,n,1] = Vert3Cor[-1,CoreBound[0]-1]
                
                Div2_Corners.values[:-1,:-1,n,0] = Rad0Cor[:,CoreBound[1]+1:]
                Div2_Corners.values[:-1,-1,n,0] = Rad1Cor[:,-1]
                Div2_Corners.values[-1,:-1,n,0] = Rad2Cor[-1,CoreBound[1]+1:]
                Div2_Corners.values[-1,-1,n,0] = Rad3Cor[-1,-1]
                
                Div2_Corners.values[:-1,:-1,n,1] = Vert0Cor[:,CoreBound[1]+1:]
                Div2_Corners.values[:-1,-1,n,1] = Vert1Cor[:,-1]
                Div2_Corners.values[-1,:-1,n,1] = Vert2Cor[-1,CoreBound[1]+1:]
                Div2_Corners.values[-1,-1,n,1] = Vert3Cor[-1,-1]
                
                for j in range(len(Y)):
                    for i in range(len(X)):
                        PsinLoc.values[j,i,n] = GF.psiN(RadLoc.loc[Y[j],X[i],Attempt].values,VertLoc.loc[Y[j],X[i],Attempt].values,)            
                
                #Populate Poloidal Coordinate Data array
                
                PolVec.loc[:,:,Attempt,'XXLoc'] = Xx
                
                YVector=np.zeros((len(X),2))
                YVector[:,0] = RadLoc.values[1,:,n] - RadLoc.values[0,:,n]
                YVector[:,1] = VertLoc.values[1,:,n] - VertLoc.values[0,:,n]            
            
                for i in range(len(X)):
                    PolVec.loc[:,X[i],Attempt,'Theta'] = np.degrees(np.math.atan2(np.linalg.det([YVector[JXA-1,:],YVector[i,:]]),np.dot(YVector[JXA-1,:],YVector[i,:])))
                    if PolVec.loc[:,X[i],Attempt,'Theta'].values[0] < 0 and X[i] < JXA:
                        PolVec.loc[:,X[i],Attempt,'Theta'] = PolVec.loc[:,X[i],Attempt,'Theta'] + 360          
                
                XP_range=np.array([CoreBound[0]-1,CoreBound[0],CoreBound[1],CoreBound[1]+1])
                X_xp=np.mean(RadLoc.loc[SEP,XP_range,Attempt].values)
                Y_xp=np.mean(VertLoc.loc[SEP,XP_range,Attempt].values)
                Xpoint=np.array([X_xp,Y_xp])
                
                for index in np.arange(CoreBound[0],CoreBound[1]+1):
                    if index == CoreBound[0]:
                        PolVec.loc[:,index,Attempt,'dXP'] = round(np.sqrt((RadLoc.loc[SEP,index,Attempt].values-X_xp)**2 + (VertLoc.loc[SEP,index,Attempt].values-Y_xp)**2),5)
                    else:
                        NL = np.sqrt((RadLoc.loc[SEP,index,Attempt].values-RadLoc.loc[SEP,index-1,Attempt].values)**2 + (VertLoc.loc[SEP,index,Attempt].values-VertLoc.loc[SEP,index-1,Attempt].values)**2)
                        PolVec.loc[:,index,Attempt,'dXP']= PolVec.loc[:,index-1,Attempt,'dXP']+NL
             
                PolVec.loc[:,:,Attempt,'dXP_norm'] = PolVec.loc[:,:,Attempt,'dXP'].values/np.max(PolVec.loc[:,:,Attempt,'dXP'].values)
        
        # Condensed Data Loading Command! :
        
        if AddNew is not None:
            if AddNew not in self.Parameter:
                self.Parameter.append(AddNew)
                P = len(self.Parameter)
        
        for p in self.Parameter:
            if p not in self.PARAM.keys() or 'AVG' in Attempts: 
                self.PARAM[p] = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = self.PARAMDICT[p])
                for n in range(N):
                    Attempt = Attempts[n]
                    if Attempt == 'AVG':
                        self.PARAM[p].values[:,:,n] = self.PARAM[p].values[:,:,:-1].mean(2)
                    else:    
                        DRT = '{}/Attempt{}'.format(BASEDRT, str(Attempt))   #Generate path
                        try:
                            RawData = np.loadtxt('{}/Output/{}{}'.format(DRT, p, str(Attempt)),usecols = (3))
                        except Exception as err:
                            print(err)
                            try:
                                 RawData = np.loadtxt('{}/Output2/{}{}'.format(DRT, p, str(Attempt)),usecols = (3))
                            except Exception as err:
                                print(err)
                                print('Parameter {} not found for Attempt {}. Creating NAN Array'.format(p, str(Attempt)))
                                self.PARAM[p].values[:,:,n] = np.nan
                                RawData=[]
                        
                        if len(RawData) > 0:        
                            if RawData.size == XDIM*YDIM:
                                self.PARAM[p].values[:,:,n] = RawData.reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
                            elif RawData.size == XDIM*YDIM*2:
                                self.PARAM[p].values[:,:,n] = RawData.reshape((2*YDIM,XDIM))[1+YDIM:2*YDIM-1,XMin:XMax+1]
                            
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
        
        
        try:
            WallFile=np.loadtxt('{}/mesh.extra'.format(BASEDRT))
        except:
            print('mesh.extra file not found! Using vvfile.ogr instead')
            WallFile=None
            
        VVFILE = np.loadtxt('{}/vvfile.ogr'.format(BASEDRT))
        
        #Save all values into self dictionaries
        
        self.KW['BASEDRT'] = BASEDRT
        self.KW['TOPDRT'] = TOPDRT
        self.KW['JXA'] = JXA
        self.KW['JXI'] = JXI
        self.KW['SEP'] = SEP
        self.Attempts = Attempts
        self.VVFILE = VVFILE
        self.WallFile = WallFile
        self.GF = GF
        self.Xx = Xx
        self.Yy = Yy
        self.Xpoint = Xpoint
        self.N = N
        self.P = P
        self.RadCoords['XMin'] = XMin
        self.RadCoords['YYLoc'] = YYLoc        
        self.RadCoords['PsinLoc'] = PsinLoc
        self.RadCoords['RadLoc'] = RadLoc
        self.RadCoords['VertLoc'] = VertLoc
        self.RadCoords['Core_Corners'] = Core_Corners
        self.RadCoords['Div1_Corners'] = Div1_Corners
        self.RadCoords['Div2_Corners'] = Div2_Corners
        self.PolVec = PolVec
        if EXP:
            self.RadCoords['PsinAvg'] = PsinAvg
            if 'd3d' not in Shot:
                    self.RadCoords['RmidAvg'] = RmidAvg
        else:
            self.RadCoords['PsinAvg'] = []
            self.RadCoords['RmidAvg'] = []
        
    def Add(self, P1, P2, Name=None):
        if type(P1) == str:
            if P1 not in self.Parameter:
                self._LoadSOLPSData(AddNew=P1)
            P1 = self.PARAM[P1]
            
        if type(P2) == str:
            if P2 not in self.Parameter:
                self._LoadSOLPSData(AddNew=P2)
            P2 = self.PARAM[P2]
            
        P3 = P1 + P2
        
        if Name is not None:
            self.PARAM[Name] = P3     
        return P3
    
    def Subtract(self, P1, P2, Name=None):
        if type(P1) == str:
            if P1 not in self.Parameter:
                self._LoadSOLPSData(AddNew=P1)
            P1 = self.PARAM[P1]
            
        if type(P2) == str:
            if P2 not in self.Parameter:
                self._LoadSOLPSData(AddNew=P2)
            P2 = self.PARAM[P2]
            
        P3 = P1 - P2
        
        if Name is not None:
            self.PARAM[Name] = P3     
        return P3
    
    def Multiply(self, P1, P2, Name=None):
        if type(P1) == str:
            if P1 not in self.Parameter:
                self._LoadSOLPSData(AddNew=P1)
            P1 = self.PARAM[P1]
            
        if type(P2) == str:
            if P2 not in self.Parameter:
                self._LoadSOLPSData(AddNew=P2)
            P2 = self.PARAM[P2]
            
        P3 = P1 * P2
        
        if Name is not None:
            self.PARAM[Name] = P3     
        return P3
    
    def Divide(self, P1, P2, Name=None):
        if type(P1) == str:
            if P1 not in self.Parameter:
                self._LoadSOLPSData(AddNew=P1)
            P1 = self.PARAM[P1]
            
        if type(P2) == str:
            if P2 not in self.Parameter:
                self._LoadSOLPSData(AddNew=P2)
            P2 = self.PARAM[P2]
            
        P3 = P1 / P2
        
        if Name is not None:
            self.PARAM[Name] = P3     
        return P3
            
    def GetRadCoords(self,RADC, Offset=[0,0]):
        #Separate method to handle radial coordinate switching
        
        #ADAPTIVE SEPARATRIX LOCATOR???

        PsinOffset = Offset[0]
        RadOffset = Offset[1]
        
        if RADC == 'Y':
            RR = self.RadCoords['YYLoc']
            Rexp = None
            Rstr = 'Radial Cell Index Y'        
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
                Rexp=[]
                print('No experimental radial coordinates available!')
            Rstr = 'm'
        elif RADC == 'rrsep':
            RR = self.RadCoords['RadLoc'].copy()
            for n in range(self.N):
                PsinSep=np.abs(self.RadCoords['PsinLoc'].loc[:,:,self.Attempts[n]]-1)
                Sign = np.sign(self.RadCoords['PsinLoc'].loc[:,:,self.Attempts[n]]-1)
                RADSEP=self.RadCoords['RadLoc'].loc[:,:,self.Attempts[n]].where(PsinSep==PsinSep.min(axis=0)).values.flatten('F')
                RADSEP=RADSEP[~np.isnan(RADSEP)]
                RRsepRad = self.RadCoords['RadLoc'].loc[:,:,self.Attempts[n]] - RADSEP
                VERTSEP=self.RadCoords['VertLoc'].loc[:,:,self.Attempts[n]].where(PsinSep==PsinSep.min(axis=0)).values.flatten('F')
                VERTSEP=VERTSEP[~np.isnan(VERTSEP)]
                RRsepVert = self.RadCoords['VertLoc'].loc[:,:,self.Attempts[n]] - VERTSEP
                RR.loc[:,:,self.Attempts[n]] = np.sqrt(RRsepRad**2 + RRsepVert**2)*Sign
            Rexp = None
            Rstr = '$R-R_{sep}$ (m)'           
        else:
            print('Invalid Radial Coordinate specified')
            
        return RR, Rexp, Rstr
    
    def B2Mesh(self,Parameter=None,**kwargs):
        #Plot B2 Meshgrid
        
        for key, value in self.KW.items():
            if key not in kwargs.keys():
                kwargs[key] = value
        
        Attempts = self.Attempts
        Shot = self.Shot
        VVFILE = self.VVFILE
        WallFile = self.WallFile
        
        RadLoc = self.RadCoords['RadLoc']
        VertLoc = self.RadCoords['VertLoc']
        Core_Corners = self.RadCoords['Core_Corners']
        Div1_Corners = self.RadCoords['Div1_Corners']
        Div2_Corners = self.RadCoords['Div2_Corners']
        
        Markers = kwargs['Markers']
        Publish = kwargs['Publish']
        DIVREG = kwargs['DIVREG']
        JXA = kwargs['JXA']
        JXI = kwargs['JXI']
        SEP = kwargs['SEP']   
        B2MeshKW = kwargs
        
        if B2MeshKW['AX']:
            ax = B2MeshKW['AX']
        else:
            fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(14,10))
            
        #Assume B2 Meshgrid is same for all Attempts, so just use data from 1st Attempt    
            
        XD1=Div1_Corners.loc[:,:,Attempts[0],'X'].values
        YD1=Div1_Corners.loc[:,:,Attempts[0],'Y'].values
        
        XD2=Div2_Corners.loc[:,:,Attempts[0],'X'].values
        YD2=Div2_Corners.loc[:,:,Attempts[0],'Y'].values
        
        XCC=Core_Corners.loc[:,:,Attempts[0],'X'].values
        YCC=Core_Corners.loc[:,:,Attempts[0],'Y'].values   
            
        PD1 = np.zeros((XD1.shape[0]-1,XD1.shape[1]-1)) 
        PD2 = np.zeros((XD2.shape[0]-1,XD2.shape[1]-1)) 
        PCC = np.zeros((XCC.shape[0]-1,XCC.shape[1]-1)) 
        
        IM1 = ax.pcolormesh(XCC,YCC,PCC,edgecolors='black',cmap='binary')
        
        if DIVREG:
            
            IM2 = ax.pcolormesh(XD1,YD1,PD1,edgecolors='black',cmap='binary')
            
            IM3 = ax.pcolormesh(XD2,YD2,PD2,edgecolors='black',cmap='binary')
        
        if Markers:
            ax.plot(RadLoc.loc[:,JXA,Attempts[0]],VertLoc.loc[:,JXA,Attempts[0]],color='orange',linewidth=2.0)
            ax.plot(RadLoc.loc[:,JXI,Attempts[0]],VertLoc.loc[:,JXI,Attempts[0]],color='gold',linewidth=2.0)
            
            SEP = int((XCC.shape[0]-1)/2)
            ax.plot(XCC[SEP,:],YCC[SEP,:],'r:',XD1[SEP,:],YD1[SEP,:],'r:',XD2[SEP,:],YD2[SEP,:],'r:',linewidth=2.0)
        
        if WallFile is not None:
            for i in range(len(WallFile[:,0])):
                ax.plot((WallFile[i,0],WallFile[i,2]),(WallFile[i,1],WallFile[i,3]),'k-',linewidth=3.0)
        else:
            ax.plot(VVFILE[:,0]/1000,VVFILE[:,1]/1000)

        if Publish != []:
            ax.set_title('B2.5 Meshgrid \n {}'.format(Publish[0]))
        else:
            ax.set_title('B2.5 Meshgrid \n Discharge {} Attempt {}'.format(Shot, str(Attempts[0])))
            
        ax.set_xlabel('Radial Coordinate r (m)')
        ax.set_ylabel('Vertical Coordinate z (m)')
        ax.set_aspect('equal')
        ax.grid()
    
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
        Attempts = self.Attempts
        Shot = self.Shot
        VVFILE = self.VVFILE
        WallFile = self.WallFile
        PolVec = self.PolVec
        Markers = kwargs['Markers']
        POLC = kwargs['POLC']
        CMAP = kwargs['CMAP']
        Colors = kwargs['Colors']
        Publish = kwargs['Publish']
        CoreBound = kwargs['CoreBound']
        DIVREG = kwargs['DIVREG']
        LINTHRESH = kwargs['LINTHRESH']
        AVG = kwargs['AVG']
        JXA = kwargs['JXA']
        JXI = kwargs['JXI']
        SEP = kwargs['SEP']
        MESHGRID = kwargs['MESHGRID']
        Offset = [kwargs['PsinOffset'],kwargs['RadOffset']]             
        ContKW = kwargs
        
        RadLoc = self.RadCoords['RadLoc']
        VertLoc = self.RadCoords['VertLoc']
        Core_Corners = self.RadCoords['Core_Corners']
        Div1_Corners = self.RadCoords['Div1_Corners']
        Div2_Corners = self.RadCoords['Div2_Corners']
        XMin = self.RadCoords['XMin']
        
        RR, Rexp, Rstr = self.GetRadCoords(ContKW['RADC'], Offset)
        
        N = np.arange(0,N)
        
        if POLC == 'theta':
            PP = PolVec.loc[:,:,:,'Theta']
            PolXLbl = r'Poloidal Angle $\theta$' 
        elif POLC == 'X':
            PP = PolVec.loc[:,:,:,'XXLoc']
            PolXLbl = r'Poloidal Cell Index $X$'
        elif POLC == 'dXP':
            PP = PolVec.loc[:,:,:,'dXP']
            PolXLbl = r'Distance from X-Point $m$'
        elif POLC == 'dXP_norm':
            PP = PolVec.loc[:,:,:,'dXP_norm']
            PolXLbl = r'Normalized distance from X-Point'
        
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
            if AVG:
                if 'AVG' not in Attempts:
                    print('Taking Average over Attempts')
                    Attempts.append('AVG')
                    self._LoadSOLPSData()
                    PARAM = self.PARAM[pn].copy()
                    RadLoc = self.RadCoords['RadLoc']
                    VertLoc = self.RadCoords['VertLoc']
                N=[-1]            
            
            if ContKW['SUBTRACT']:
                for n in np.arange(1,N):
                    PARAM.values[:,:,n] = (PARAM.values[:,:,0]-PARAM[:,:,n])
                PARAM.values[:,:,0] = PARAM.values[:,:,0]-PARAM[:,:,0]
                CMAP = cm.RdBu
                
            # NEED TO FIX -> PREVENT PARAM FROM BEING OVERWRITTEN EVERY TIME
            
            NORM=colors.Normalize()
            
            if ContKW['LOG10'] == 1:
                
                if np.any([p<0 for p in PARAM.values]):
                    PARAM.values[np.abs(PARAM.values) <= LINTHRESH] = 0.0
                    PARAM.values[PARAM.values>1] = np.log10(PARAM.values[PARAM.values>1])
                    PARAM.values[PARAM.values<-1] = -1*np.log10(np.abs(PARAM.values[PARAM.values<-1]))
                    CMAP = cm.RdBu
                    NORM=colors.TwoSlopeNorm(0.0)
                else:
                    PARAM.values[PARAM.values==0.0] = np.min(PARAM.values[PARAM.values!=0.0])
                    PARAM.values = np.log10(PARAM.values)
                levs = np.arange(np.floor(PARAM.values.min()),np.ceil(PARAM.values.max()),0.5)
                    
            elif ContKW['LOG10'] == 2:
                
                if np.any([p<0 for p in PARAM.values]):
                    NPARAM = np.abs(PARAM.values[PARAM.values<0])
                    NPARAM[NPARAM<=0.1] = np.nan
                    PARAM.values[np.abs(PARAM.values)<=0.1] = np.nan
                    #lev_exp = np.arange(-1*(np.ceil(np.log10(np.nanmax(NPARAM)))+1), np.ceil(np.log10(np.nanmax(PARAM.values)))+1,0.5)
                    #levs = np.sign(lev_exp)*np.power(10, np.abs(lev_exp))
                    levs=ContKW['LVN']
                    np.set_printoptions(threshold=np.inf)
                    CMAP = cm.RdBu
                    NORM=colors.SymLogNorm(LINTHRESH)
                else:
                    PARAM.values[PARAM.values==0.0] = np.min(PARAM.values[PARAM.values!=0.0])
                    lev_exp = np.arange(np.floor(np.log10(np.nanmin(PARAM.values)))-1, np.ceil(np.log10(np.nanmax(PARAM.values)))+1,0.5)
                    levs = np.power(10, lev_exp)
                    NORM=colors.LogNorm()
                    
                
            elif type(ContKW['LVN'])==list:
                levs=ContKW['LVN']
                
            else:
                
                if np.any([p<0 for p in PARAM.values]):
                    CMAP = cm.RdBu
                    NORM=colors.CenteredNorm(0.0)
                levs = np.linspace(np.floor(PARAM.values.min()),np.ceil(PARAM.values.max()),ContKW['LVN'])
            
            for n in N:
                if ContKW['AX']:
                    ax = ContKW['AX']
                else:
                    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(14,10))
                
                if ContKW['PHYS']:
                    
                    XD1=Div1_Corners.loc[:,:,Attempts[n],'X'].values
                    YD1=Div1_Corners.loc[:,:,Attempts[n],'Y'].values
                    
                    XD2=Div2_Corners.loc[:,:,Attempts[n],'X'].values
                    YD2=Div2_Corners.loc[:,:,Attempts[n],'Y'].values
                    
                    XCC=Core_Corners.loc[:,:,Attempts[n],'X'].values
                    YCC=Core_Corners.loc[:,:,Attempts[n],'Y'].values
                    
                    if DIVREG:
                        
                        NORM.vmin = np.nanmin(PARAM.values)
                        NORM.vmax = np.nanmax(PARAM.values)
                        
                        IM2 = ax.pcolormesh(XD1,YD1,PARAM.loc[:,1.0:CoreBound[0],Attempts[n]],
                                            edgecolors=MESHGRID,cmap=CMAP,norm=NORM)
                        
                        IM3 = ax.pcolormesh(XD2,YD2,PARAM.loc[:,CoreBound[1]+2:,Attempts[n]],
                                            edgecolors=MESHGRID,cmap=CMAP,norm=NORM)
                    else:
                        
                        NORM.vmin = np.nanmin(PARAM.loc[:,CoreBound[0]+1:CoreBound[1]+1,Attempts[n]].values)
                        NORM.vmax = np.nanmax(PARAM.loc[:,CoreBound[0]+1:CoreBound[1]+1,Attempts[n]].values)
                            
                    IM1 = ax.pcolormesh(XCC,YCC,PARAM.loc[:,CoreBound[0]+1:CoreBound[1]+1,Attempts[n]],
                                        edgecolors=MESHGRID,cmap=CMAP,norm=NORM)
                    
                    if Markers:                 
                        ax.plot(RadLoc.values[:,(JXA-XMin),n],VertLoc.values[:,(JXA-XMin),n],color='orange',linewidth=2)
                        ax.plot(RadLoc.values[:,(JXI-XMin),n],VertLoc.values[:,(JXI-XMin),n],color='gold',linewidth=2)
                    
                    if WallFile is not None:
                        for i in range(len(WallFile[:,0])):
                            ax.plot((WallFile[i,0],WallFile[i,2]),(WallFile[i,1],WallFile[i,3]),'k-',linewidth=3.0)
                    else:
                        ax.plot(VVFILE[:,0]/1000,VVFILE[:,1]/1000)
                    
                    SEP = int((XCC.shape[0]-1)/2)
                    
                    ax.plot(XCC[SEP,:],YCC[SEP,:],'r:',XD1[SEP,:],YD1[SEP,:],'r:',XD2[SEP,:],YD2[SEP,:],'r:',linewidth=2)
                    ax.set_xlabel('Radial Location (m)')
                    ax.set_ylabel('Vertical Location (m)')
                    ax.set_aspect('equal')
                
                else:
                    '''
                    if ContKW['LOG10'] == 2:
                        if DIVREG and POLC != 'theta':
                            IM2 = ax.contourf(PP.loc[:,1:CoreBound[0]-1,Attempts[n]],
                                              RR.loc[:,1:CoreBound[0]-1,Attempts[n]],
                                              PARAM.loc[:,1:CoreBound[0]-1,Attempts[n]],
                                              levs,colors=Colors,cmap=CMAP,norm=colors.LogNorm())
                            
                            IM3 = ax.contourf(PP.loc[:,CoreBound[1]+1:,Attempts[n]],
                                              RR.loc[:,CoreBound[1]+1:,Attempts[n]],
                                              PARAM.loc[:,CoreBound[1]+1:,Attempts[n]],
                                              levs,colors=Colors,cmap=CMAP,norm=colors.LogNorm())
                        
                        IM1 = ax.contourf(PP.loc[:,CoreBound[0]:CoreBound[1],Attempts[n]],
                                          RR.loc[:,CoreBound[0]:CoreBound[1],Attempts[n]],
                                          PARAM.loc[:,CoreBound[0]:CoreBound[1],Attempts[n]],
                                          levs,norm=colors.LogNorm(),colors=Colors,cmap=CMAP)
                    else:
                    '''
                    if DIVREG and POLC != 'theta':
                        NORM.vmin = np.nanmin(PARAM.values)
                        NORM.vmax = np.nanmax(PARAM.values)
                        
                        ax.plot(PP.loc[SEP,:,Attempts[n]],
                                RR.loc[SEP,:,Attempts[n]],
                                color='red',linewidth=2,
                                linestyle=':',label='Separatrix')                                                                        #Full Separatrix

                        
                        IM2 = ax.contourf(PP.loc[:,1:CoreBound[0],Attempts[n]],
                                          RR.loc[:,1:CoreBound[0],Attempts[n]],
                                          PARAM.loc[:,1:CoreBound[0],Attempts[n]],
                                          levels=levs,cmap=CMAP,norm=NORM)
                            
                        IM3 = ax.contourf(PP.loc[:,CoreBound[1]+2:,Attempts[n]],
                                          RR.loc[:,CoreBound[1]+2:,Attempts[n]],
                                          PARAM.loc[:,CoreBound[1]+2:,Attempts[n]],
                                          levels=levs,cmap=CMAP,norm=NORM)
                    else:
                        NORM.vmin = np.nanmin(PARAM.loc[:,CoreBound[0]+1:CoreBound[1]+1,Attempts[n]].values)
                        NORM.vmax = np.nanmax(PARAM.loc[:,CoreBound[0]+1:CoreBound[1]+1,Attempts[n]].values)
                        
                        ax.plot(PP.loc[SEP,CoreBound[0]+1:CoreBound[1]+1,Attempts[n]],
                                RR.loc[SEP,CoreBound[0]+1:CoreBound[1]+1,Attempts[n]],
                                color='red',linewidth=2,linestyle=':',label='Separatrix')                                                #Core-only Separatrix


                    IM1 = ax.contourf(PP.loc[:,CoreBound[0]+1:CoreBound[1]+1,Attempts[n]],
                                     RR.loc[:,CoreBound[0]+1:CoreBound[1]+1,Attempts[n]],
                                     PARAM.loc[:,CoreBound[0]+1:CoreBound[1]+1,Attempts[n]],
                                     levels=levs,cmap=CMAP,norm=NORM)
                    
                    if Markers:
                        ax.plot(PP.loc[:,JXA,Attempts[n]],RR.loc[:,JXA,Attempts[n]],color='orange',linewidth=2)                          #Outer Midplane
                        ax.plot(PP.loc[:,JXI,Attempts[n]],RR.loc[:,JXI,Attempts[n]],color='gold',linewidth=2)                            #Inner Midplane
                    
                    ax.plot(PP.loc[1:SEP,CoreBound[0]+1,Attempts[n]],RR.loc[1:SEP,CoreBound[0]+1,Attempts[n]],color='Blue',linewidth=2)    #Inner PFR Boundary
                    ax.plot(PP.loc[1:SEP,CoreBound[1]+1,Attempts[n]],RR.loc[1:SEP,CoreBound[1]+1,Attempts[n]],color='Blue',linewidth=2)    #Outer PFR Boundary
                    
                    ax.set_xlabel(PolXLbl)
                    ax.set_ylabel(Rstr) 

                if Publish != []:
                    ax.set_title('Attempt {} {}'.format(Publish[n], PARAM.name))
                else:
                    ax.set_title('Discharge {} Attempt {}\n{}'.format(Shot, str(Attempts[n]), PARAM.name))
                    
                SM=SM=cm.ScalarMappable(NORM,CMAP)    
                plt.colorbar(SM,ax=ax)

                #a.set_xticklabels(['%.1f' % i for i in a.get_xticks()], fontsize='x-large')
                #a.set_yticklabels(['%.1f' % j for j in a.get_yticks()], fontsize='x-large')
                #a.tick_params(labelsize=20)
                
                if ContKW['GRID']:
                    plt.grid()
                
                if ContKW['SAVE']:
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
            AVG = kwargs['AVG']
            JXA = kwargs['JXA']
            JXI = kwargs['JXI']
            SEP = kwargs['SEP']
            POLC = kwargs['POLC']
            SURF = kwargs['SURF']
            CoreBound = kwargs['CoreBound']
            PolKW = kwargs
            
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

                if PolKW['AVG']:
                    if 'AVG' not in Attempts:
                        print('Taking Average over Attempts')
                        Attempts.append('AVG')
                        self._LoadSOLPSData()
                        PARAM = self.PARAM[pn].copy()
                    N=[-1]  

            if POLC == 'theta':
                PP = PolVec.loc[:,:,:,'Theta']
                PolXLbl = r'Poloidal Angle $\theta$' 
            elif POLC == 'X':
                PP = PolVec.loc[:,:,:,'XXLoc']
                PolXLbl = r'Poloidal Cell Index $X$'
            elif POLC == 'dXP':
                PP = PolVec.loc[:,:,:,'dXP']
                PolXLbl = r'Poloidal distance along separatrix, clockwise from X-point ($m$)'
            elif POLC == 'dXP_norm':
                PP = PolVec.loc[:,:,:,'dXP_norm']
                PolXLbl = r'Normalized poloidal distance along separatrix, clockwise from X-point'

            #print('Beginning Plot Sequence')                

            if PolKW['AX'] is None: # if no axis was supplied to the function create our own
                fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(14,7))
            else:
                ax = PolKW['AX']
            
            if PolKW['GRAD']:
                PARAM.values = np.gradient(PARAM.values,axis=0)
                
            if PolKW['LOG10'] == 1:
                PARAM.values[PARAM.values>1] = np.log10(PARAM.values[PARAM.values>1])
                PARAM.values[PARAM.values<-1] = -1*np.log10(np.abs(PARAM.values[PARAM.values<-1]))      
            
            if len(PlotScheme) == self.N:
                for n in range(self.N):
                    if PolKW['LOG10'] == 2:
                        ax.semilogy(PP.loc[SURF,CoreBound[0]:CoreBound[1],Attempts[n]], PARAM.loc[SURF,CoreBound[0]:CoreBound[1],Attempts[n]], PlotScheme[n], label=Publish)
                    else:
                        ax.plot(PP.loc[SURF,CoreBound[0]:CoreBound[1],Attempts[n]], PARAM.loc[SURF,CoreBound[0]:CoreBound[1],Attempts[n]], PlotScheme[n], linewidth=3, label=Publish)
            else:            
                if PolKW['LOG10'] == 2:
                    ax.semilogy(PP.loc[SURF,CoreBound[0]:CoreBound[1],:], PARAM.loc[SURF,CoreBound[0]:CoreBound[1],:], label=Publish)
                else:
                    ax.plot(PP.loc[SURF,CoreBound[0]:CoreBound[1],:], PARAM.loc[SURF,CoreBound[0]:CoreBound[1],:],linewidth=3, label=Publish)

            if Publish:
                plt.legend(Publish)
                plt.title('Poloidal ' + PARAM.name + ' along Separatrix')
            else:
                plt.legend(Attempts)
                plt.title('Discharge 0' + str(Shot) + ' Attempt(s) ' + str(Attempts) + ' Poloidal ' + PARAM.name)

            if Markers:
                plt.axvline(PP.loc[SURF,JXI,Attempts[0]],color='red',linestyle='dashed')
                plt.axvline(PP.loc[SURF,JXA,Attempts[0]],color='orange',linestyle='dashed')
                plt.axvline(PP.loc[SURF,CoreBound[0],Attempts[0]],color='black',linestyle='dashed')
                plt.axvline(PP.loc[SURF,CoreBound[1],Attempts[0]],color='black',linestyle='dashed')
                  
            plt.xlabel(PolXLbl)
            plt.ylabel(PARAM.name)
            
            if PolKW['GRID']:
                plt.grid()
            
            self.PolKW = PolKW
    
    def RadProf(self,Parameter=None,**kwargs):
        
        for key, value in self.KW.items():
            if key not in kwargs.keys():
                kwargs[key] = value
        
        N = self.N
        Shot = self.Shot
        Attempts = self.Attempts
        Publish = kwargs['Publish']
        AVG = kwargs['AVG']
        JXA = kwargs['JXA']
        SEP = kwargs['SEP']
        RADC = kwargs['RADC']
        PolSlc = kwargs['PolSlc']
        Markers = kwargs['Markers']
        PlotScheme = kwargs['PlotScheme']
        Offset = [kwargs['PsinOffset'],kwargs['RadOffset']]         
        RadProfKW = kwargs
        
        RR, Rexp, Rstr = self.GetRadCoords(RADC,Offset)
        
        if not Publish:
            Publish=Attempts
        
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
                    quit()
            
            if RadProfKW['AVG'] is True:
                if 'AVG' not in Attempts:
                    print('Taking Average over Attempts')
                    Attempts.append('AVG')
                    self._LoadSOLPSData()
                    PARAM = self.PARAM[pn].copy()
                
            #print('Beginning Plot Sequence')
            
            if RadProfKW['AX'] is None: # if no axis was supplied to the function create our own
                fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(14,7))
            else:
                ax = RadProfKW['AX']
                    
            if RadProfKW['LOG10'] == 1:
                PARAM.values[PARAM.values>1] = np.log10(PARAM.values[PARAM.values>1])
                PARAM.values[PARAM.values<-1] = -1*np.log10(np.abs(PARAM.values[PARAM.values<-1]))    
                #y_exp = np.arange(np.floor(np.nanmin(PARAM.values)), np.ceil(np.nanmax(PARAM.values))+1,2)    
            
            if len(PlotScheme) == N:
                for n in range(N):
                    if RadProfKW['GRAD'] is True: #Needs To Be Debugged, think it's a problem with PARAM.values not being copied
                        if RADC == 'rrsep':
                            PARAM.values = np.gradient(PARAM.loc[:,JXA,Attempts[n]].values,RR.loc[:,JXA,Attempts[n]].values,axis=0)
                        else:
                            print('Gradient can not be calculated unless radial coordinate is in meters (RADC=rrsep)')
                    
                    if RadProfKW['LOG10'] == 2:
                        ax.semilogy(RR.loc[:,JXA,Attempts[n]], PARAM.loc[:,JXA,Attempts[n]],PlotScheme[n], label=Publish[n])
                    else:
                        ax.plot(RR.loc[:,JXA,Attempts[n]], PARAM.loc[:,JXA,Attempts[n]],PlotScheme[n],linewidth=3, label=Publish[n])
                        
                if RadProfKW['EXP'] is True and RadProfKW['RADC'] != 'rrsep' and RadProfKW['RADC'] != 'Y':
                    if 'd3d' not in Shot:
                        if pn == 'Ne':
                            NemidAvg = self.ExpDict['NemidAvg']
                            ErrNe = self.ExpDict['ErrNe']
                            ax.errorbar(Rexp,NemidAvg,yerr=ErrNe,fmt='o',markersize=7,linewidth=3,capsize=7,color=PlotScheme[0][0], label='Exp')
                        elif pn == 'Te':
                            TemidAvg = self.ExpDict['TemidAvg']
                            ErrTe = self.ExpDict['ErrTe']
                            ax.errorbar(Rexp,TemidAvg,yerr=ErrTe,fmt='o',markersize=7,linewidth=3,capsize=7,color=PlotScheme[0][0], label='Exp')
                    
                    if 'd3d' in Shot:
                        if pn == 'Ne':
                            PsinNe = Rexp[0]
                            Ned3d = self.ExpDict['Ned3d']
                            ax.plot(PsinNe,Ned3d,'o',markersize=7,linewidth=3,color=PlotScheme[0][0], label='Exp')
                        elif pn == 'Te':
                            PsinTe = Rexp[1]
                            Ted3d = self.ExpDict['Ted3d']
                            ax.plot(PsinTe,Ted3d,'o',markersize=7,linewidth=3,color=PlotScheme[0][0], label='Exp')
                        elif pn == 'Ti':
                            PsinTi = Rexp[2]
                            Tid3d = self.ExpDict['Tid3d']
                            ax.plot(PsinTi,Tid3d,'o',markersize=7,linewidth=3,color=PlotScheme[0][0], label='Exp')
            else:
                if RadProfKW['LOG10'] == 2:
                    ax.semilogy(RR.loc[:,JXA,:], PARAM.loc[:,JXA,:],linewidth=3, label=Publish)
                else:
                    ax.plot(RR.loc[:,JXA,:], PARAM.loc[:,JXA,:],linewidth=3, label=Publish)

                if RadProfKW['EXP'] is True and RadProfKW['RADC'] != 'rrsep' and RadProfKW['RADC'] != 'Y':
                    if 'd3d' not in Shot:
                        if pn == 'Ne':
                            NemidAvg = self.ExpDict['NemidAvg']
                            ErrNe = self.ExpDict['ErrNe']
                            ax.errorbar(Rexp,NemidAvg,yerr=ErrNe,fmt='o',markersize=7,linewidth=3,capsize=7, label='Exp')
                        elif pn == 'Te':
                            TemidAvg = self.ExpDict['TemidAvg']
                            ErrTe = self.ExpDict['ErrTe']
                            ax.errorbar(Rexp,TemidAvg,yerr=ErrTe,fmt='o',markersize=7,linewidth=3,capsize=7, label='Exp')
                    
                    if 'd3d' in Shot:
                        if pn == 'Ne':
                            PsinNe = Rexp[0]
                            Ned3d = self.ExpDict['Ned3d']
                            ax.plot(PsinNe,Ned3d,'o',markersize=7,linewidth=3, label='Exp')
                        elif pn == 'Te':
                            PsinTe = Rexp[1]
                            Ted3d = self.ExpDict['Ted3d']
                            ax.plot(PsinTe,Ted3d,'o',markersize=7,linewidth=3, label='Exp')
                        elif pn == 'Ti':
                            PsinTi = Rexp[2]
                            Tid3d = self.ExpDict['Tid3d']
                            ax.plot(PsinTi,Tid3d,'o',markersize=7,linewidth=3, label='Exp')
            
            if RadProfKW['AX'] is None: #End of current axes reached if not pre-supplied
                
                ax.set_xlabel(Rstr)
                
                if RadProfKW['LOG10'] == 1:
                    ax.set_ylabel('Log_10 of {}'.format(PARAM.name))
                else:    
                    ax.set_ylabel(PARAM.name)
    
                ax.legend([*Publish])
                ax.set_title('Radial Midplane {}'.format(PARAM.name))
                
                if Markers:
                    if RADC == 'psin':
                        ax.axvline(1.0,color='red',linewidth=2,linestyle='dashed')
                    elif RADC == 'rrsep':
                        ax.axvline(0.0,color='red',linewidth=2,linestyle='dashed')    
                    else:
                        ax.axvline(RR.loc[SEP,JXA,Attempts[0]],color='red',linewidth=2,linestyle='dashed')
                
                if RadProfKW['GRID'] is True:
                    ax.grid(b=1)
                
            self.RadProfKW = RadProfKW

        if RadProfKW['AX'] is not None:
            
            ax.set_xlabel(Rstr)
            
            if RadProfKW['LOG10'] == 1:
                ax.set_ylabel('Log_10 of {}'.format(PARAM.name))
            else:    
                ax.set_ylabel(PARAM.name)

            ax.legend([*Publish])
            ax.set_title('Attempts: {}\nRadial Midplane Profile of {}'.format(Attempts,PARAM.name))
            
            if Markers:
                if RADC == 'psin':
                    ax.axvline(1.0,color='red',linewidth=2,linestyle='dashed')
                elif RADC == 'rrsep':
                    ax.axvline(0.0,color='red',linewidth=2,linestyle='dashed')    
                else:
                    ax.axvline(RR.loc[SEP,JXA,Attempts[0]],color='red',linewidth=2,linestyle='dashed')
            
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
