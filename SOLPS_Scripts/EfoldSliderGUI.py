# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 16:15:42 2019

@author: rmreksoatmodjo
"""

from VesselPlotterNew import SOLPSPLOT
from TOOLS import TANH
from TOOLS import EXPFIT
import numpy as np
import xarray as xr
from scipy.optimize import curve_fit
from scipy.stats import binned_statistic
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.path import Path
from matplotlib.widgets import Slider, Button, CheckButtons
import json

### Input Fields ###

Shot = '12'
Attempt = ['51N']
GasLvl = 77.8
Balloon = 1

PS=['.','.','.','.','.','x']

### Setting up Base Variables ###
NeuDen = SOLPSPLOT(Shot,Attempt,['Ne','NeuDen'],EXP=True,AVG=False,PlotScheme='')
JXA = NeuDen.KW['JXA']
JXI = NeuDen.KW['JXI']
SEP = 18
CoreBound = [26,70]
Rmax = 0.01
Rmin = -0.01
Thresh=0.01
Mag=0.165
RadLoc = NeuDen.RadCoords['RadLoc']
VertLoc = NeuDen.RadCoords['VertLoc']
RR = NeuDen.GetRadCoords('rrsep',[0,0])[0]
Psin = NeuDen.GetRadCoords('psin',[0,0])[0]
DRT='{}/Attempt{}/'.format(NeuDen.KW['BASEDRT'],Attempt[0])

XP_range=np.array([CoreBound[0]-1,CoreBound[0],CoreBound[1],CoreBound[1]+1])
X_xp=np.mean(RadLoc.loc[SEP,XP_range,Attempt[-1]].values) #(RadLoc.loc[36,JXA,Attempt[-1]].values+RadLoc.loc[36,JXI,Attempt[-1]].values)/2
Y_xp=np.mean(VertLoc.loc[SEP,XP_range,Attempt[-1]].values) #(VertLoc.loc[36,JXA,Attempt[-1]].values+VertLoc.loc[36,JXI,Attempt[-1]].values)/2

f0 = JXA
p0 = [0,3.5e20,0.005,1e18,1e21]
e0 = [1e15,100]
bins=np.linspace(1,36,36)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
axcolor = 'wheat'
log = 1

RR_SEP_avg = {}
Psin_avg = {}
Ne_SOLPS_med = {}
NeuDen_SOLPS_med = {}

efold={}
efold_adj={}
yparam={}
eparam={}
x0 = {}
xi = {}
fluxpsn = {}
fluxpsnparam={}

if 'AVG' in NeuDen.Attempts:
    Attempt.append('AVG')

# Calculate Poloidal Coordinates dXP, Distance from X-Point along Separatrix!

PolCoords=['Distance from X-Point along Separatrix (m)','Normalized Distance from X-Point']
dXP = xr.DataArray(np.zeros([CoreBound[1]-CoreBound[0]+1,2]),coords=[np.arange(CoreBound[0],CoreBound[1]+1),PolCoords],dims=['Poloidal Index','Poloidal Coordinate Convention'])

for index in np.arange(CoreBound[0],CoreBound[1]+1):
    if index == CoreBound[0]:
        dXP.loc[index,PolCoords[0]] = round(np.sqrt((RadLoc.loc[SEP,index,Attempt[-1]].values-X_xp)**2 + (VertLoc.loc[SEP,index,Attempt[-1]].values-Y_xp)**2),5)
    else:
        NL = np.sqrt((RadLoc.loc[SEP,index,Attempt[-1]].values-RadLoc.loc[SEP,index-1,Attempt[-1]].values)**2 + (VertLoc.loc[SEP,index,Attempt[-1]].values-VertLoc.loc[SEP,index-1,Attempt[-1]].values)**2)
        dXP.loc[index,PolCoords[0]]=dXP.loc[index-1,PolCoords[0]]+NL
             
dXP.loc[:,PolCoords[1]]=dXP.loc[:,PolCoords[0]].values/np.max(dXP.loc[:,PolCoords[0]].values)
        
### Setting up Canvas ###

fig, ax = plt.subplots()
ax.set_frame_on(False)
ax.set_axis_off()
gs=gridspec.GridSpec(9,2,width_ratios=[5,3],height_ratios=[1,1,1,1,1,1,1,1,1],hspace=0.3)

axcontour = fig.add_subplot(gs[0:8,1])
NeuDen.Contour('NeuDen',LOG10=1,AX=axcontour, Markers=False)
axcontour.set_title('Neutral Density Contour')
axcontour.plot(X_xp,Y_xp,'X')
axcontour.margins(x=0)

fluxpsnprofile = fig.add_subplot(gs[0:3,0])
FluxPsnProf, = fluxpsnprofile.plot(np.nan,np.nan,'rX')
Lin_Fit, = fluxpsnprofile.plot(np.nan,np.nan)
FluxPsnV1 = fluxpsnprofile.axvline(np.nan)
FluxPsnV2 = fluxpsnprofile.axvline(np.nan)
fluxpsnprofile.axvline(0.0,color='k')
fluxpsnprofile.axhline(1.0,color='k')
TXT_FLXPSN=fluxpsnprofile.text(0.02,0.95,'',transform=fluxpsnprofile.transAxes,verticalalignment='top', bbox=props)
fluxpsnprofile.set_ylabel(r'$\psi_n$')
fluxpsnprofile.set_title('Shot {}, Attempt {}'.format(Shot,*Attempt))
'''
fluxpsnprofile.plot(RR.loc[:,f0,Attempt[-1]].values,Psin.loc[:,f0,Attempt[-1]].values,'*')
fluxpsnprofile.set_xlim(XLim)
fluxpsnYLim = fluxpsnprofile.get_ylim()
'''

neprofile = fig.add_subplot(gs[3:6,0],sharex=fluxpsnprofile)
NeProf, = neprofile.plot(np.nan,np.nan,'bX')
Tanh_Fit, = neprofile.plot(np.nan,np.nan)
NeV1 = neprofile.axvline(np.nan)
NeV2 = neprofile.axvline(np.nan)
TXT_TANH=neprofile.text(0.02,0.2,'',transform=neprofile.transAxes,verticalalignment='top', bbox=props)
neprofile.set_ylabel(r'Electron Density ($m^{-3}$)')
XLim = neprofile.get_xlim()
'''
NeuDen.RadProf('Ne',AX=neprofile,Markers=False,RADC='rrsep',JXA=f0)
neprofile.set_xlim(XLim)
NeYLim = neprofile.get_ylim()
'''

neudenprofile = fig.add_subplot(gs[6:9,0],sharex=neprofile)
NeuDenProf, = neudenprofile.semilogy(np.nan,np.nan,'gX')
Exp_Fit, = neudenprofile.semilogy(np.nan,np.nan)
NDV1 = neudenprofile.axvline(np.nan)
NDV2 = neudenprofile.axvline(np.nan)
TXT_EFOLD=neudenprofile.text(0.02,0.75,'',transform=neudenprofile.transAxes,verticalalignment='top', bbox=props)
TXT_EADJ=neudenprofile.text(0.02,0.6,'',transform=neudenprofile.transAxes,verticalalignment='top', bbox=props)
neudenprofile.set_xlabel(r'$R-R_{SEP}$ (m)')
neudenprofile.set_ylabel(r'Neutral D Density ($m^{-3}$)')
'''
NeuDen.RadProf('NeuDen',LOG10=log,AX=neudenprofile,Markers=False,RADC='rrsep',JXA=f0)  #,PlotScheme=['x'])
XLim = neudenprofile.get_xlim()
YLim = neudenprofile.get_ylim()
'''

# Initiate Chord Slicing Variables and Data Vectors

Slice, = axcontour.plot(np.nan,np.nan,color='orange',linewidth=1)
ChordXY, =axcontour.plot(np.nan,np.nan,'r.')
PR=RadLoc.loc[:,:,Attempt[-1]].values.flatten()
PZ=VertLoc.loc[:,:,Attempt[-1]].values.flatten()
Psin_ALL=Psin.loc[:,:,Attempt[-1]].values.flatten()
NeuDen_ALL=NeuDen.PARAM['NeuDen'].loc[:,:,Attempt[-1]].values.flatten()
Ne_ALL=NeuDen.PARAM['Ne'].loc[:,:,Attempt[-1]].values.flatten()

### Methods and Functions ###

def update(val):
    
    log=1
    
    if LogOn.get_status()[0] is True:
        log = 1
    elif LogOn.get_status()[0] is False:
        log = 0
        
    PolPos = sslide.val
    
    #RR_SOLPS = RR.loc[:,PolPos,Attempt[-1]].values
    
        #P0 = Grab X,Y Point of Inner Most (Core) Cell at PolPos
    
    P0=np.array([RadLoc.loc[1,PolPos,Attempt[-1]].values,VertLoc.loc[1,PolPos,Attempt[-1]].values])
    
        #P1 = Grab X,Y Point of Separatrix at PolPos
    
    P1=np.array([RadLoc.loc[SEP,PolPos,Attempt[-1]].values,VertLoc.loc[SEP,PolPos,Attempt[-1]].values])
    
        #Form Chord Slice between P0 and P1, extend beyond SOL to P2
    
    PP=P1-P0
    Theta=np.arctan2(PP[1],PP[0])
    displace=Thresh*np.array([-np.sin(Theta),np.cos(Theta)])
    P2=np.array([(Mag*np.cos(Theta)+P0[0]),(Mag*np.sin(Theta)+P0[1])])
    P0A=P0+displace
    P0B=P0-displace
    P2A=P2+displace
    P2B=P2-displace
    
    Bounds=Path([P0A,P2A,P2B,P0B])
    Mask=Bounds.contains_points(np.vstack((PR,PZ)).T)
    ChordX=np.ma.array(PR,mask=~Mask).compressed()
    ChordY=np.ma.array(PZ,mask=~Mask).compressed()
    
        #Mask and Compress all Coordinate, Neutral Density, and Electron Density data not contained by Mask
    
    R_SEP=np.sqrt((P1[0]-P0[0])**2 + (P1[1]-P0[1])**2)
    RR_SEP=np.sqrt((ChordX-P1[0])**2 + (ChordY-P1[1])**2)*np.sign(np.sqrt((ChordX-P0[0])**2 + (ChordY-P0[1])**2)-R_SEP)
    Psin_SOLPS = np.ma.array(Psin_ALL,mask=~Mask).compressed()
    NeuDen_SOLPS = np.ma.array(NeuDen_ALL,mask=~Mask).compressed() 
    Ne_SOLPS = np.ma.array(Ne_ALL,mask=~Mask).compressed()    

        #Use histogram-like binned statistics to take median over chord data
    
    RR = binned_statistic(RR_SEP, RR_SEP,statistic='mean',bins=36)[0]
    RR_SEP_avg[PolPos] = RR[~np.isnan(RR)]
    
    PS = binned_statistic(RR_SEP,Psin_SOLPS,statistic='mean',bins=36)[0] 
    Psin_avg[PolPos] = PS[~np.isnan(PS)]
    
    ND = binned_statistic(RR_SEP,NeuDen_SOLPS,statistic='median',bins=36)[0]
    NeuDen_SOLPS_med[PolPos] = ND[~np.isnan(ND)]
    
    NE = binned_statistic(RR_SEP,Ne_SOLPS,statistic='median',bins=36)[0]
    Ne_SOLPS_med[PolPos] = NE[~np.isnan(NE)]
    
        #Set GUI plot data to new values 
    
    ChordXY.set_data(ChordX,ChordY)
    Slice.set_data(np.array([P0,P2]).transpose())
    FluxPsnProf.set_data(RR_SEP_avg[PolPos],Psin_avg[PolPos])
    NeuDenProf.set_data(RR_SEP_avg[PolPos],NeuDen_SOLPS_med[PolPos])
    NeProf.set_data(RR_SEP_avg[PolPos],Ne_SOLPS_med[PolPos])

    if PolPos in x0.keys(): 
        FluxPsnV1.set_xdata(x0[PolPos])
        FluxPsnV2.set_xdata(xi[PolPos])
        
        NeV1.set_xdata(x0[PolPos])
        NeV2.set_xdata(xi[PolPos])
        
        NDV1.set_xdata(x0[PolPos])
        NDV2.set_xdata(xi[PolPos])
    else:
        FluxPsnV1.set_xdata(np.nan)
        FluxPsnV2.set_xdata(np.nan)
        
        NeV1.set_xdata(np.nan)
        NeV2.set_xdata(np.nan)
        
        NDV1.set_xdata(np.nan)
        NDV2.set_xdata(np.nan)         
        
    if PolPos in yparam.keys():
        Tanh_Fit.set_data(RR_SEP_avg[PolPos],TANH(RR_SEP_avg[PolPos],*yparam[PolPos]))
        TXT_TANH.set_text(r'$n_{{e,PED}}$={:.3e}$m^{{-3}}$, $\Delta n_e$={:.3f}mm'.format(yparam[PolPos][1]+yparam[PolPos][3],2000*yparam[PolPos][2]))
    else:
        Tanh_Fit.set_data(np.nan,np.nan)
        TXT_TANH.set_text('')

    if PolPos in eparam.keys():
        Exp_Fit.set_data(RR_SEP_avg[PolPos], EXPFIT(RR_SEP_avg[PolPos],*eparam[PolPos]))
        TXT_EFOLD.set_text('e-folding length={:.3f}mm'.format(efold[PolPos]))
        TXT_EADJ.set_text('Adjusted e-folding length={:.3f}mm'.format(efold_adj[PolPos]))
        
        Lin_Fit.set_data(RR_SEP_avg[PolPos],(fluxpsnparam[PolPos][0]*RR_SEP_avg[PolPos]+fluxpsnparam[PolPos][1]))
        TXT_FLXPSN.set_text('Flux Expansion={:.3f}mm'.format(fluxpsn[PolPos]))
    else:
        Exp_Fit.set_data(np.nan,np.nan)
        TXT_EFOLD.set_text('')
        TXT_EADJ.set_text('')
        
        Lin_Fit.set_data(np.nan,np.nan)
        TXT_FLXPSN.set_text('')
        '''    
    if Fixed.get_status()[0] is True:
        neudenprofile.set_xlim(XLim)
        neudenprofile.set_ylim(YLim)
        neprofile.set_xlim(XLim)
        neprofile.set_ylim(NeYLim)
        fluxpsnprofile.set_xlim(XLim)
        fluxpsnprofile.set_ylim(fluxpsnYLim)
        '''
    
    neprofile.relim()
    neprofile.autoscale_view(True,True,True)
    neudenprofile.relim()
    neudenprofile.autoscale_view(True,True,True)
    fluxpsnprofile.relim()
    fluxpsnprofile.autoscale_view(True,True,True)
    fig.canvas.draw()

### Basic Buttons and Sliders ###

axslide = fig.add_subplot(gs[8,1], facecolor=axcolor) #plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
sslide = Slider(axslide, 'Poloidal\nSurface', CoreBound[0], CoreBound[1], valinit=f0, valfmt='%0.0f', valstep=1.0)

sslide.on_changed(update)

def arrowclick(event):
    if event.key == 'right' and sslide.val<CoreBound[1]:
        sslide.set_val(sslide.val+1)
    elif event.key == 'left' and sslide.val>CoreBound[0]:
        sslide.set_val(sslide.val-1)
    else:
        pass        

cid = fig.canvas.mpl_connect('key_press_event', arrowclick)

resetax = plt.axes([0.125, 0.025, 0.05, 0.05])
Reset = Button(resetax, 'Reset (JXA)', color=axcolor, hovercolor='0.975')

def reset(event):
    sslide.reset()
    
Reset.on_clicked(reset)

jxiax = plt.axes([0.185, 0.025, 0.05, 0.05])
JXI_Button = Button(jxiax, 'JXI', color=axcolor, hovercolor='0.975')

def jxi_set(event):
    sslide.set_val(JXI)
    
JXI_Button.on_clicked(jxi_set)

logonax = plt.axes([0.245, 0.025, 0.075, 0.05], facecolor=axcolor)
LogOn = CheckButtons(logonax, [r'Log$_{10}(n_D)$'],[True])
LogOn.on_clicked(update)

fixedax = plt.axes([0.330, 0.025, 0.05, 0.05], facecolor=axcolor)
Fixed = CheckButtons(fixedax, ['Fix\nAxes'],[True])
Fixed.on_clicked(update)

### Modified TANH Fitting Function ###

tanhfitax = plt.axes([0.390, 0.025, 0.075, 0.05])
TanhFit = Button(tanhfitax, 'Create\nTanh Fit', color=axcolor, hovercolor='0.975')

def tanhfit(event):
    PolPos=sslide.val
    #RR_SOLPS = RR_SEP_avg[PolPos]
    
    if PolPos not in yparam.keys():
        #RR_SOLPS = RR.loc[:,PolPos,Attempt[-1]].values
        #Ne_SOLPS = NeuDen.PARAM['Ne'].loc[:,PolPos,Attempt[-1]].values
        yfit=curve_fit(TANH, RR_SEP_avg[PolPos], Ne_SOLPS_med[PolPos],p0)
        yparam[PolPos] = yfit[0]
        x0[PolPos] = yparam[PolPos][0]+yparam[PolPos][2]
        xi[PolPos] = yparam[PolPos][0]-yparam[PolPos][2]

    print('Poloidal Slice {:0.0f}: r0={:.3e}m, h={:.3e}m^-3, d={:.3e}m, b={:.3e}m^-3, m={:.3e}m^-4'.format(PolPos,*yparam[PolPos]))
    TXT_TANH.set_text(r'$n_{{e,PED}}$={:.3e}$m^{{-3}}$, $\Delta n_e$={:.3f}mm'.format(yparam[PolPos][1]+yparam[PolPos][3],2000*yparam[PolPos][2]))
    
    Tanh_Fit.set_data(RR_SEP_avg[PolPos],TANH(RR_SEP_avg[PolPos],*yparam[PolPos]))
    
    FluxPsnV1.set_xdata(x0[PolPos])
    FluxPsnV2.set_xdata(xi[PolPos])
        
    NeV1.set_xdata(x0[PolPos])
    NeV2.set_xdata(xi[PolPos])
        
    NDV1.set_xdata(x0[PolPos])
    NDV2.set_xdata(xi[PolPos]) 
        
TanhFit.on_clicked(tanhfit)

### Exponential Fitting Function ###

expfitax = plt.axes([0.475, 0.025, 0.075, 0.05])
ExpFit = Button(expfitax, 'Create\nExp. Fit', color=axcolor, hovercolor='0.975')

def expfit(event):
    PolPos=sslide.val
    #xr = x0[PolPos]
    #xri = xi[PolPos]
    #RR_SOLPS = RR.loc[:,PolPos,Attempt[-1]].values
    RR_i = np.where((RR_SEP_avg[PolPos]>(xi[PolPos])) & (RR_SEP_avg[PolPos]<(x0[PolPos])))[0]
    RR_i = np.arange((RR_i[0]-1),(RR_i[-1]+2),1)
    RR_exp=RR_SEP_avg[PolPos][RR_i]
    if len(RR_exp) < 3:
        RR_i = np.arange((RR_i[0]-1),(RR_i[-1]+2),1)
        RR_exp= RR.loc[RR_i,PolPos,Attempt[-1]].values
        
    if PolPos not in eparam.keys():
        
        exfit=curve_fit(EXPFIT,RR_exp,NeuDen_SOLPS_med[PolPos][RR_i],e0)
        eparam[PolPos] = exfit[0]
        efold[PolPos] = 1000/eparam[PolPos][1]
        
        #Psin_SOLPS = Psin.loc[RR_i,PolPos,Attempt[-1]].values
        fluxpsnparam[PolPos] = np.polyfit(RR_exp,Psin_avg[PolPos][RR_i],1)
        fluxpsn[PolPos] = fluxpsnparam[PolPos][0]/fluxpsnparam[JXA][0]
        efold_adj[PolPos] = fluxpsn[PolPos]*efold[PolPos]
                
    print('Exponential fit from r-r_sep={:.3e}m to r-r_sep={:.3e}m'.format(RR_exp[0],RR_exp[-1]))
    print('A0={:.3e}, lambda={:.3f}'.format(*eparam[PolPos]))
    print('Poloidal Slice {:0.0f}: e-folding length={:.3f}mm'.format(PolPos,efold[PolPos]))
    print('Slope={:.3f}, Flux Expansion={:.3f}, Adjusted e-folding length={:.3f}'.format(fluxpsnparam[PolPos][0],fluxpsn[PolPos],efold_adj[PolPos]))

    Exp_Fit.set_data(RR_exp, EXPFIT(RR_exp,*eparam[PolPos]))
    TXT_EFOLD.set_text('e-folding length={:.3f}mm'.format(efold[PolPos]))
    TXT_EADJ.set_text('Adjusted e-folding length={:.3f}mm'.format(efold_adj[PolPos]))
    Lin_Fit.set_data(RR_exp,(fluxpsnparam[PolPos][0]*RR_exp+fluxpsnparam[PolPos][1]))
    TXT_FLXPSN.set_text('Flux Expansion={:.3f}mm'.format(fluxpsn[PolPos]))

ExpFit.on_clicked(expfit)

### Poloidal Neutral e-Folding Length Plotting Function ###

wholeax = plt.axes([0.575, 0.025, 0.135, 0.05])
WholeFit = Button(wholeax, 'WHOLE POLOIDAL PLOT', color=axcolor, hovercolor='0.975')

bothax = plt.axes([0.720, 0.025, 0.125, 0.05], facecolor=axcolor)
BothPlot = CheckButtons(bothax, ['Plot Raw e-Fold'],[False])

def wholefit(event):
    wholeFig, wholeAx = plt.subplots()
    reset(event)
    tanhfit(event)
    expfit(event)
    for n in range(CoreBound[0],CoreBound[1]+1):
        sslide.set_val(n)
        print('Calculating e-fold length for Poloidal Position {}'.format(n))
        tanhfit(event)
        expfit(event)
  
    x_adj,y_adj = zip(*sorted(efold_adj.items()))
    wholeAx.plot(dXP.loc[:,PolCoords[0]].values,y_adj,'b^')
    wholeAx.set_title('Shot {} Attempt {} neutral e-folding lengths'.format(Shot,Attempt[-1]))
    wholeAx.set_xlabel(PolCoords[0])
    wholeAx.set_ylabel('e_folding length (mm)')
    wholeAx.axvline(dXP.loc[JXA,PolCoords[0]].values,color='black')
    wholeAx.axvline(dXP.loc[JXI,PolCoords[0]].values,color='orange')
    wholeAx.axvline(dXP.loc[CoreBound[0],PolCoords[0]].values,color='red')
    
    if BothPlot.get_status()[0] == True:
        x,y = zip(*sorted(efold.items()))  
        wholeAx.plot(dXP.loc[:,PolCoords[0]].values,y,'r^-')
        wholeAx.legend(['Adjusted e-folding length','Outer Midplane', 'Inner Midplane','X-Point','Raw e-folding length'])
    else:
        wholeAx.legend(['Adjusted e-folding length','Outer Midplane', 'Inner Midplane','X-Point'])
    
    #wholeAx.xaxis.set_ticks(np.arange(20,75,5))
    starty, endy = wholeAx.get_ylim()
    wholeAx.yaxis.set_ticks(np.linspace(0,np.round(endy),11))
    wholeAx.axvline(dXP.loc[CoreBound[1],PolCoords[0]].values,color='red')
    wholeAx.grid()
    
    secax=wholeAx.secondary_xaxis('top',functions=(lambda x: x / np.max(dXP[:,0].values), lambda x: x * np.max(dXP[:,0].values)))
    secax.set_xlabel(PolCoords[1])
    plt.show

exportax = plt.axes([0.855, 0.025, 0.1, 0.05])
ExportButton = Button(exportax, 'Export\nData', color=axcolor, hovercolor='0.975')

def export(event):    
    # Save LFS and HFS Adjusted e-Folding lengths in an external .json file
    if JXA or JXI not in yparam.keys():
        
        print('Processing Inner and Outer Midplane Neutral Data')
        
        reset(event)
        tanhfit(event)
        expfit(event)
        
        jxi_set(event)
        tanhfit(event)
        expfit(event)
    
    print('Formatting data...')
    
    efold_plot={}
    efold_plot['gaslvl']=GasLvl
    efold_plot['balloon']=Balloon
    efold_plot['LFS'] = efold_adj[JXA]
    efold_plot['HFS'] = efold_adj[JXI]
    efold_plot['LFS_NeuDen'] = float(NeuDen.PARAM['NeuDen'].loc[SEP,JXA,Attempt[-1]].values)
    efold_plot['HFS_NeuDen'] = float(NeuDen.PARAM['NeuDen'].loc[SEP,JXI,Attempt[-1]].values)
    efold_plot['LFS_NePED'] = yparam[JXA][1]+yparam[JXA][3]
    efold_plot['HFS_NePED'] = yparam[JXI][1]+yparam[JXI][3]
    efold_plot['LFS_PedWidth'] = 2000*yparam[JXA][2]
    efold_plot['HFS_PedWidth'] = 2000*yparam[JXI][2]
    
    dXP_dict=dXP.to_dict()
    
    export_data=[efold_plot,efold,efold_adj,dXP_dict]
    
    with open('{}efold_data_{}.json'.format(DRT,Attempt[0]),'w') as fp:
        json.dump(export_data,fp,indent=2)
        
    print('Data exported succesfully!')
    
WholeFit.on_clicked(wholefit)
ExportButton.on_clicked(export)

update(JXA)

plt.show()