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
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.path import Path
from matplotlib.widgets import Slider, Button, CheckButtons
import json

### Input Fields ###

Shot = 'gas025'
Attempt = ['69']
GasLvl = 6.77
Balloon = 0

PS=['.','.','.','.','.','x']

### Setting up Base Variables ###
polfig, polax = plt.subplots()

NeuDen = SOLPSPLOT(Shot,Attempt,['Ne','NeuDen'],EXP=True,AVG=False,PlotScheme='')
JXA = NeuDen.KW['JXA']
JXI = NeuDen.KW['JXI']
SEP = 18
CoreBound = NeuDen.KW['CoreBound']
Rmax = 0.01
Rmin = -0.01
Thresh=0.01
Mag=0.2
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
efold={}
efold_adj={}
yparam={}
eparam={}
x0 = {}
xi = {}
fluxpsn = {}
fluxpsnparam={}
log = 2

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
axcolor = 'wheat'

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
fluxpsnprofile.plot(RR.loc[:,f0,Attempt[-1]].values,Psin.loc[:,f0,Attempt[-1]].values,'*')
fluxpsnprofile.axvline(0.0,color='k')
fluxpsnprofile.axhline(1.0,color='k')
#fluxpsnprofile.set_xlim(XLim)
fluxpsnYLim = fluxpsnprofile.get_ylim()
fluxpsnprofile.set_ylabel(r'$\psi_n$')
fluxpsnprofile.set_xlabel('')
fluxpsnprofile.set_title('Shot {}, Attempt {}'.format(Shot,*Attempt))

neprofile = fig.add_subplot(gs[3:6,0],sharex=fluxpsnprofile)
NeProf, = neprofile.plot(np.nan,np.nan,'X')
'''
NeuDen.RadProf('Ne',AX=neprofile,Markers=False,RADC='rrsep',JXA=f0)
neprofile.set_xlim(XLim)
NeYLim = neprofile.get_ylim()
'''
neprofile.set_xlabel('')
neprofile.set_title('')

neudenprofile = fig.add_subplot(gs[6:9,0],sharex=neprofile)
NeuDenProf, = neudenprofile.plot(np.nan,np.nan,'X')
'''
NeuDen.RadProf('NeuDen',LOG10=log,AX=neudenprofile,Markers=False,RADC='rrsep',JXA=f0)  #,PlotScheme=['x'])
XLim = neudenprofile.get_xlim()
YLim = neudenprofile.get_ylim()
'''
neudenprofile.set_title('')

# Implement Chord Slicing from EireneContourPlot here, to get linear cuts of data!

Slice, = axcontour.plot(np.nan,np.nan,color='orange',linewidth=1)
ChordXY, =axcontour.plot(np.nan,np.nan,'r.')
PR=RadLoc.loc[:,:,Attempt[-1]].values.flatten()
PZ=VertLoc.loc[:,:,Attempt[-1]].values.flatten()
NeuDen_ALL=NeuDen.PARAM['NeuDen'].loc[:,:,Attempt[-1]].values.flatten()
Ne_ALL=NeuDen.PARAM['Ne'].loc[:,:,Attempt[-1]].values.flatten()

### Methods and Functions ###

def update(val):
    
    log=2
    
    if LogOn.get_status()[0] is True:
        log = 2
    elif LogOn.get_status()[0] is False:
        log = 0
        
    PolPos = sslide.val
    RR_SOLPS = RR.loc[:,PolPos,Attempt[-1]].values
    fluxpsnprofile.clear()
    
    #P0 = Grab X,Y Point of Inner Most (Core) Cell at PolPos
    
    P0=np.array([RadLoc.loc[1,PolPos,Attempt[-1]].values,VertLoc.loc[1,PolPos,Attempt[-1]].values])
    
    #P1 = Grab X,Y Point of Separatrix at PolPos
    
    P1=np.array([RadLoc.loc[SEP,PolPos,Attempt[-1]].values,VertLoc.loc[SEP,PolPos,Attempt[-1]].values])
    
    #Form Chord Line between P0 and P1, extend beyond SOL to P2
    
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
    
    ChordXY.set_data(ChordX,ChordY)
    
    R_SEP=np.sqrt((P1[0]-P0[0])**2 + (P1[1]-P0[1])**2)
    RR_SEP=np.sqrt((ChordX-P1[0])**2 + (ChordY-P1[1])**2)*np.sign(np.sqrt((ChordX-P0[0])**2 + (ChordY-P0[1])**2)-R_SEP)
    NeuDen_SOLPS = np.ma.array(NeuDen_ALL,mask=~Mask).compressed() 
    Ne_SOLPS = np.ma.array(Ne_ALL,mask=~Mask).compressed()    
    
    Slice.set_data(np.array([P0,P2]).transpose())
    NeuDenProf.set_data(RR_SEP,NeuDen_SOLPS)
    NeProf.set_data(RR_SEP,Ne_SOLPS)
    
    #NeuDen.RadProf('NeuDen',LOG10=log,AX=neudenprofile,Markers=False,RADC='rrsep',JXA=PolPos) #,PlotScheme=['x'])
    #NeuDen.RadProf('Ne',AX=neprofile,Markers=False,RADC='rrsep',JXA=PolPos)
    
    fluxpsnprofile.plot(RR.loc[:,PolPos,Attempt[-1]].values,Psin.loc[:,PolPos,Attempt[-1]].values,'*')
    fluxpsnprofile.axvline(0.0,color='k')
    fluxpsnprofile.axhline(1.0,color='k')
    fluxpsnprofile.set_ylabel(r'$\psi_n$')
    
    if PolPos in x0.keys(): 
        neprofile.axvline(x0[PolPos])
        neudenprofile.axvline(x0[PolPos])
        fluxpsnprofile.axvline(x0[PolPos])
        neprofile.axvline(xi[PolPos])
        neudenprofile.axvline(xi[PolPos])
        fluxpsnprofile.axvline(xi[PolPos])
        
    if PolPos in yparam.keys():
        neprofile.plot(RR_SOLPS,TANH(RR_SOLPS,*yparam[PolPos]))
        neprofile.text(0.02,0.2,r'$n_{{e,PED}}$={:.3e}$m^{{-3}}$, $\Delta n_e$={:.3f}mm'.format(yparam[PolPos][1]+yparam[PolPos][3],2000*yparam[PolPos][2]),transform=neprofile.transAxes,verticalalignment='top', bbox=props)
        
    if PolPos in eparam.keys():
        neudenprofile.plot(RR_SOLPS, EXPFIT(RR_SOLPS,*eparam[PolPos]))
        neudenprofile.text(0.02,0.75,'e-folding length={:.3f}mm'.format(efold[PolPos]),transform=neudenprofile.transAxes,verticalalignment='top', bbox=props)
        neudenprofile.text(0.02,0.6,'Adjusted e-folding length={:.3f}mm'.format(efold_adj[PolPos]),transform=neudenprofile.transAxes,verticalalignment='top', bbox=props)
        
        fluxpsnprofile.plot(RR_SOLPS,(fluxpsnparam[PolPos][0]*RR_SOLPS+fluxpsnparam[PolPos][1]))
        fluxpsnprofile.text(0.02,0.95,'Flux Expansion={:.3f}mm'.format(fluxpsn[PolPos]),transform=fluxpsnprofile.transAxes,verticalalignment='top', bbox=props)
        '''    
    if Fixed.get_status()[0] is True:
        neudenprofile.set_xlim(XLim)
        neudenprofile.set_ylim(YLim)
        neprofile.set_xlim(XLim)
        neprofile.set_ylim(NeYLim)
        fluxpsnprofile.set_xlim(XLim)
        fluxpsnprofile.set_ylim(fluxpsnYLim)
        '''
    fluxpsnprofile.set_title('Shot {}, Attempt {}'.format(Shot,*Attempt))
    fluxpsnprofile.set_xlabel('')
    neprofile.set_title('')
    neprofile.set_xlabel('')
    neudenprofile.set_title('')
    
    neprofile.relim()
    neprofile.autoscale_view(True,True,True)
    neudenprofile.relim()
    neudenprofile.autoscale_view(True,True,True)
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
Reset = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

def reset(event):
    sslide.reset()
    
Reset.on_clicked(reset)

logonax = plt.axes([0.185, 0.025, 0.075, 0.05], facecolor=axcolor)
LogOn = CheckButtons(logonax, [r'Log$_{10}(n_D)$'],[True])
LogOn.on_clicked(update)

fixedax = plt.axes([0.270, 0.025, 0.075, 0.05], facecolor=axcolor)
Fixed = CheckButtons(fixedax, ['Fix Axes'],[True])
Fixed.on_clicked(update)

### Modified TANH Fitting Function ###

tanhfitax = plt.axes([0.355, 0.025, 0.1, 0.05])
TanhFit = Button(tanhfitax, 'Create Tanh Fit', color=axcolor, hovercolor='0.975')

def tanhfit(event):
    PolPos=sslide.val
    RR_SOLPS = RR.loc[:,PolPos,Attempt[-1]].values
    
    if PolPos not in yparam.keys():
        RR_SOLPS = RR.loc[:,PolPos,Attempt[-1]].values
        Ne_SOLPS = NeuDen.PARAM['Ne'].loc[:,PolPos,Attempt[-1]].values
        yfit=curve_fit(TANH, RR_SOLPS, Ne_SOLPS,p0)
        yparam[PolPos] = yfit[0]
        x0[PolPos] = yparam[PolPos][0]+yparam[PolPos][2]
        xi[PolPos] = yparam[PolPos][0]-yparam[PolPos][2]

    print('Poloidal Slice {:0.0f}: r0={:.3f}m, h={:.3e}m^-3, d={:.3f}m, b={:.3e}m^-3, m={:.3e}m^-4'.format(PolPos,*yparam[PolPos]))
    neprofile.text(0.02,0.2,r'$n_{{e,PED}}$={:.3e}$m^{{-3}}$, $\Delta n_e$={:.3f}mm'.format(yparam[PolPos][1]+yparam[PolPos][3],2000*yparam[PolPos][2]),transform=neprofile.transAxes,verticalalignment='top', bbox=props)
    
    neprofile.plot(RR_SOLPS,TANH(RR_SOLPS,*yparam[PolPos]))
    neprofile.axvline(x0[PolPos])
    neudenprofile.axvline(x0[PolPos])
    fluxpsnprofile.axvline(x0[PolPos])
    neprofile.axvline(xi[PolPos])
    neudenprofile.axvline(xi[PolPos])
    fluxpsnprofile.axvline(xi[PolPos])
        
TanhFit.on_clicked(tanhfit)

### Exponential Fitting Function ###

expfitax = plt.axes([0.465, 0.025, 0.1, 0.05])
ExpFit = Button(expfitax, 'Create Exp. Fit', color=axcolor, hovercolor='0.975')

def expfit(event):
    PolPos=sslide.val
    xr = x0[PolPos]
    xri = xi[PolPos]
    RR_SOLPS = RR.loc[:,PolPos,Attempt[-1]].values
    RR_i = np.where((RR_SOLPS>(xri)) & (RR_SOLPS<(xr)))[0]
    RR_SOLPS=RR_SOLPS[RR_i]
    if len(RR_SOLPS) < 3:
        RR_i = np.arange((RR_i[0]-1),(RR_i[-1]+2),1)
        RR_SOLPS= RR.loc[RR_i,PolPos,Attempt[-1]].values
        
    if PolPos not in eparam.keys():
        #P0 = Grab X,Y Point of Inner Most (Core) Cell at PolPos
        P0=np.array([RadLoc.loc[1,PolPos,Attempt[-1]].values,VertLoc.loc[1,PolPos,Attempt[-1]].values])
        #P1 = Grab X,Y Point of Separatrix at PolPos
        P1=np.array([RadLoc.loc[SEP,PolPos,Attempt[-1]].values,VertLoc.loc[SEP,PolPos,Attempt[-1]].values])
        #Form Chord Line between P0 and P1, extend beyond SOL
        Thresh=0.01
        print('Threshhold={}m'.format(Thresh))
        PP=P1-P0
        Theta=np.arctan(PP[1]/PP[0])
        displace=Thresh*np.array([-np.sin(Theta),np.cos(Theta)])
        P0A=P0+displace
        P0B=P0-displace
        P1A=P1+displace
        P1B=P1-displace
        
        Bounds=Path([P0A,P1A,P1B,P0B])
    
        #NeuDen_SOLPS = Grab all neutral density cell values within Threshhold of Chord Line
        
        NeuDen_SOLPS = NeuDen.PARAM['NeuDen'].loc[RR_i,PolPos,Attempt[-1]].values        
        
        #RR_SOLPS = R-Rsep value for each NeuDen_SOLPS value
        
        

        exfit=curve_fit(EXPFIT,RR_SOLPS,NeuDen_SOLPS,e0)
        eparam[PolPos] = exfit[0]
        efold[PolPos] = 1000/eparam[PolPos][1]
        
        Psin_SOLPS = Psin.loc[RR_i,PolPos,Attempt[-1]].values
        fluxpsnparam[PolPos] = np.polyfit(RR_SOLPS,Psin_SOLPS,1)
        fluxpsn[PolPos] = fluxpsnparam[PolPos][0]/fluxpsnparam[JXA][0]
        efold_adj[PolPos] = fluxpsn[PolPos]*efold[PolPos]
                
    print('Exponential fit from r-r_sep={:.3f}m to r-r_sep={:.3f}m'.format(RR_SOLPS[0],RR_SOLPS[-1]))
    print('A0={:.3e}, lambda={:.3f}'.format(*eparam[PolPos]))
    print('Poloidal Slice {:0.0f}: e-folding length={:.3f}mm'.format(PolPos,efold[PolPos]))
    print('Slope={:.3f}, Flux Expansion={:.3f}, Adjusted e-folding length={:.3f}'.format(fluxpsnparam[PolPos][0],fluxpsn[PolPos],efold_adj[PolPos]))

    neudenprofile.plot(RR_SOLPS, EXPFIT(RR_SOLPS,*eparam[PolPos]))
    neudenprofile.text(0.02,0.75,'e-folding length={:.3f}mm'.format(efold[PolPos]),transform=neudenprofile.transAxes,verticalalignment='top', bbox=props)
    neudenprofile.text(0.02,0.6,'Adjusted e-folding length={:.3f}mm'.format(efold_adj[PolPos]),transform=neudenprofile.transAxes,verticalalignment='top', bbox=props)
    fluxpsnprofile.plot(RR_SOLPS,(fluxpsnparam[PolPos][0]*RR_SOLPS+fluxpsnparam[PolPos][1]))
    fluxpsnprofile.text(0.02,0.95,'Flux Expansion={:.3f}mm'.format(fluxpsn[PolPos]),transform=fluxpsnprofile.transAxes,verticalalignment='top', bbox=props)

ExpFit.on_clicked(expfit)

### Poloidal Neutral e-Folding Length Plotting Function ###

wholeax = plt.axes([0.635, 0.025, 0.135, 0.05])
WholeFit = Button(wholeax, 'WHOLE POLOIDAL PLOT', color=axcolor, hovercolor='0.975')

bothax = plt.axes([0.780, 0.025, 0.125, 0.05], facecolor=axcolor)
BothPlot = CheckButtons(bothax, ['Plot Raw e-Fold'],[False])

def wholefit(event):
    wholeAx = polax
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
    
    # Save LFS and HFS Adjusted e-Folding lengths in an external .json file
    
    efold_plot={}
    efold_plot['gaslvl']=GasLvl
    efold_plot['balloon']=Balloon
    efold_plot['LFS'] = efold_adj[JXA]
    efold_plot['HFS'] = efold_adj[JXI]
    efold_plot['LFS_NeuDen'] = float(NeuDen.PARAM['NeuDen'].loc[SEP,JXA,Attempt[-1]].values)
    efold_plot['HFS_NeuDen'] = float(NeuDen.PARAM['NeuDen'].loc[SEP,JXI,Attempt[-1]].values)
    
    dXP_dict=dXP.to_dict()
    
    export_data=[efold_plot,efold,efold_adj,dXP_dict]
    
    with open('{}efold_data_{}.json'.format(DRT,Attempt[0]),'w') as fp:
        json.dump(export_data,fp,indent=2)
    
WholeFit.on_clicked(wholefit)

plt.show()