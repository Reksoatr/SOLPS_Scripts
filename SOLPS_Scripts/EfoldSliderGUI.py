# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 16:15:42 2019

@author: rmreksoatmodjo
"""

from VesselPlotterNew import SOLPSPLOT
from TOOLS import TANH
from TOOLS import EXPFIT
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, Button, CheckButtons

Shot = '12'
Attempt = ['101','102','103','104','105']
PS=['.','.','.','.','.','-']

NeuDen = SOLPSPLOT(Shot,Attempt,['Ne','NeuDen'],EXP=False,AVG=True,PlotScheme=PS)
JXA = NeuDen.KW['JXA']
JXI = NeuDen.KW['JXI']
CoreBound = NeuDen.KW['CoreBound']
Rmax = 0.01
Rmin = -0.01
RadLoc = NeuDen.RadCoords['RadLoc']
VertLoc = NeuDen.RadCoords['VertLoc']
RR = NeuDen.GetRadCoords('rrsep',[0,0])[0]
Psin = NeuDen.GetRadCoords('psin',[0,0])[0]

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

fig, ax = plt.subplots()
ax.set_frame_on(False)
ax.set_axis_off()
gs=gridspec.GridSpec(9,2,width_ratios=[5,3],height_ratios=[1,1,1,1,1,1,1,1,1],hspace=0.3)

axcontour = fig.add_subplot(gs[0:8,1])
NeuDen.Contour('NeuDen',LOG10=1,AX=axcontour, Markers=False)
axcontour.set_title('Neutral Density Contour')
l, = axcontour.plot(RadLoc.loc[:,f0,Attempt[-1]],VertLoc.loc[:,f0,Attempt[-1]],color='Red',linewidth=3)
axcontour.margins(x=0)

neudenprofile = fig.add_subplot(gs[0:3,0]) #plt.axes([0.25, 0.2, 0.4, 0.6], facecolor=axcolor)
NeuDen.RadProf('NeuDen',LOG10=log,AX=neudenprofile,Markers=False,RADC='rrsep',JXA=f0)  #,PlotScheme=['x'])
XLim = neudenprofile.get_xlim()
YLim = neudenprofile.get_ylim()

neprofile = fig.add_subplot(gs[3:6,0])
NeuDen.RadProf('Ne',AX=neprofile,Markers=False,RADC='rrsep',JXA=f0)
neprofile.set_xlim(XLim)
NeYLim = neprofile.get_ylim()
neprofile.set_title('')

fluxpsnprofile = fig.add_subplot(gs[6:9,0])
fluxpsnprofile.plot(RR.loc[:,f0,Attempt[-1]].values,Psin.loc[:,f0,Attempt[-1]].values,'*')
fluxpsnprofile.axvline(0.0,color='k')
fluxpsnprofile.axhline(1.0,color='k')
fluxpsnprofile.set_xlim(XLim)
fluxpsnYLim = fluxpsnprofile.get_ylim()
fluxpsnprofile.set_xlabel(r'$R-R_{sep}$ (m)')
fluxpsnprofile.set_ylabel(r'$\psi_n$')

axslide = fig.add_subplot(gs[8,1], facecolor=axcolor) #plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
sslide = Slider(axslide, 'Poloidal\nSurface', CoreBound[0], CoreBound[1], valinit=f0, valfmt='%0.0f', valstep=1.0)
#axslide.set

def update(val):
    
    log=2
    
    if LogOn.get_status()[0] is True:
        log = 2
    elif LogOn.get_status()[0] is False:
        log = 0
        
    PolPos = sslide.val
    RR_SOLPS = RR.loc[:,PolPos,Attempt[-1]].values
    l.set_xdata(RadLoc.loc[:,PolPos,Attempt[-1]])
    l.set_ydata(VertLoc.loc[:,PolPos,Attempt[-1]])
    neudenprofile.clear()
    neprofile.clear()
    fluxpsnprofile.clear()
    
    NeuDen.RadProf('NeuDen',LOG10=log,AX=neudenprofile,Markers=False,RADC='rrsep',JXA=PolPos) #,PlotScheme=['x'])
    NeuDen.RadProf('Ne',AX=neprofile,Markers=False,RADC='rrsep',JXA=PolPos)
    fluxpsnprofile.plot(RR.loc[:,PolPos,Attempt[-1]].values,Psin.loc[:,PolPos,Attempt[-1]].values,'*')
    fluxpsnprofile.axvline(0.0,color='k')
    fluxpsnprofile.axhline(1.0,color='k')
    fluxpsnprofile.set_xlabel(r'$R-R_{sep}$ (m)')
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
        
    if PolPos in eparam.keys():
        neudenprofile.plot(RR_SOLPS, EXPFIT(RR_SOLPS,*eparam[PolPos]))
        neudenprofile.text(0.01,0.85,'e-folding length={:.3f}mm'.format(efold[PolPos]),transform=neudenprofile.transAxes,verticalalignment='top', bbox=props)
        neudenprofile.text(0.01,0.75,'Adjusted e-folding length={:.3f}mm'.format(efold_adj[PolPos]),transform=neudenprofile.transAxes,verticalalignment='top', bbox=props)
        fluxpsnprofile.plot(RR_SOLPS,(fluxpsnparam[PolPos][0]*RR_SOLPS+fluxpsnparam[PolPos][1]))
        fluxpsnprofile.text(0.01,0.95,'Flux Expansion={:.3f}mm'.format(fluxpsn[PolPos]),transform=fluxpsnprofile.transAxes,verticalalignment='top', bbox=props)
    
    if Fixed.get_status()[0] is True:
        neudenprofile.set_xlim(XLim)
        neudenprofile.set_ylim(YLim)
        neprofile.set_xlim(XLim)
        neprofile.set_ylim(NeYLim)
        fluxpsnprofile.set_xlim(XLim)
        fluxpsnprofile.set_ylim(fluxpsnYLim)

    neprofile.set_title('')
    fig.canvas.draw_idle()


sslide.on_changed(update)

def arrowclick(event):
    if event.key == 'right' and sslide.val<CoreBound[1]:
        sslide.set_val(sslide.val+1)
    elif event.key == 'left' and sslide.val>CoreBound[0]:
        sslide.set_val(sslide.val-1)
    else:
        pass        

cid = fig.canvas.mpl_connect('key_press_event', arrowclick)

resetax = plt.axes([0.125, 0.025, 0.05, 0.06])
Reset = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

logonax = plt.axes([0.185, 0.025, 0.075, 0.06], facecolor=axcolor)
LogOn = CheckButtons(logonax, [r'Log$_{10}$ Scale'],[True])
LogOn.on_clicked(update)

fixedax = plt.axes([0.270, 0.025, 0.075, 0.06], facecolor=axcolor)
Fixed = CheckButtons(fixedax, ['Fix Axes'],[True])
Fixed.on_clicked(update)

def reset(event):
    sslide.reset()
    
Reset.on_clicked(reset)

tanhfitax = plt.axes([0.355, 0.025, 0.1, 0.06])
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
    neprofile.plot(RR_SOLPS,TANH(RR_SOLPS,*yparam[PolPos]))
    neprofile.axvline(x0[PolPos])
    neudenprofile.axvline(x0[PolPos])
    fluxpsnprofile.axvline(x0[PolPos])
    neprofile.axvline(xi[PolPos])
    neudenprofile.axvline(xi[PolPos])
    fluxpsnprofile.axvline(xi[PolPos])
        
TanhFit.on_clicked(tanhfit)

expfitax = plt.axes([0.465, 0.025, 0.1, 0.06])
ExpFit = Button(expfitax, 'Create Exp. Fit', color=axcolor, hovercolor='0.975')

def expfit(event):
    PolPos=sslide.val
    xr = x0[PolPos]
    xri = xi[PolPos]
    RR_SOLPS = RR.loc[:,PolPos,Attempt[-1]].values
    RR_i = np.where((RR_SOLPS>(xri)) & (RR_SOLPS<(xr)))[0]
    RR_SOLPS=RR_SOLPS[RR_i]
        
    if PolPos not in eparam.keys():
        NeuDen_SOLPS = NeuDen.PARAM['NeuDen'].loc[RR_i,PolPos,Attempt[-1]].values
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
    neudenprofile.text(0.01,0.85,'e-folding length={:.3f}mm'.format(efold[PolPos]),transform=neudenprofile.transAxes,verticalalignment='top', bbox=props)
    neudenprofile.text(0.01,0.75,'Adjusted e-folding length={:.3f}mm'.format(efold_adj[PolPos]),transform=neudenprofile.transAxes,verticalalignment='top', bbox=props)
    fluxpsnprofile.plot(RR_SOLPS,(fluxpsnparam[PolPos][0]*RR_SOLPS+fluxpsnparam[PolPos][1]))
    fluxpsnprofile.text(0.01,0.95,'Flux Expansion={:.3f}mm'.format(fluxpsn[PolPos]),transform=fluxpsnprofile.transAxes,verticalalignment='top', bbox=props)

ExpFit.on_clicked(expfit)

wholeax = plt.axes([0.75, 0.025, 0.1, 0.06])
WholeFit = Button(wholeax, 'WHOLE POLOIDAL PLOT', color=axcolor, hovercolor='0.975')

def wholefit(event):
    wholeFig, wholeAx = plt.subplots(1,1)
    reset(event)
    tanhfit(event)
    expfit(event)
    for n in range(CoreBound[0],CoreBound[1]+1):
        sslide.set_val(n)
        print('Calculating e-fold length for Poloidal Position {}'.format(n))
        tanhfit(event)
        expfit(event)
    x,y = zip(*sorted(efold.items()))    
    x_adj,y_adj = zip(*sorted(efold_adj.items()))
    wholeAx.plot(x,y,'r')
    wholeAx.plot(x_adj,y_adj,'b')
    wholeAx.set_title('Shot {} Attempt {} neutral e-folding lengths'.format(Shot,Attempt[-1]))
    wholeAx.set_xlabel('Poloidal Grid Index')
    wholeAx.set_ylabel('e_folding length (mm)')
    wholeAx.axvline(JXA,color='black')
    wholeAx.axvline(JXI,color='orange')    
    wholeAx.legend(['Raw e-folding length','Adjusted e-folding length','Outer Midplane', 'Inner Midplane'])
    wholeAx.xaxis.set_ticks(np.arange(20,75,5))
    starty, endy = wholeAx.get_ylim()
    wholeAx.yaxis.set_ticks(np.arange(0,np.round(endy),5))
    wholeAx.grid()
    plt.show
    
WholeFit.on_clicked(wholefit)

plt.show()