# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 16:15:42 2019

@author: rmreksoatmodjo
"""

from VesselPlotterNew import SOLPSPLOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

Shot = '25'
Attempt = '153'

NeuDen = SOLPSPLOT(Shot,Attempt,'NeuDen')
JXA = NeuDen.KW['JXA']
JXI = NeuDen.KW['JXI']
RadLoc = NeuDen.RadCoords['RadLoc']
VertLoc = NeuDen.RadCoords['VertLoc']

f0 = JXA
axcolor = 'lightgoldenrodyellow'

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)

axcontour = fig.add_subplot(122)
NeuDen.Contour('NeuDen',LOG10=1,AX=axcontour, Markers=False)
axcontour.set_title('Neutral Density')
l, = axcontour.plot(RadLoc.loc[:,f0,Attempt],VertLoc.loc[:,f0,Attempt],color='Red',linewidth=3)
axcontour.margins(x=0)

axprofile = fig.add_subplot(221) #plt.axes([0.25, 0.2, 0.4, 0.6], facecolor=axcolor)
NeuDen.RadProf('NeuDen',LOG10=2,AX=axprofile,Markers=False,JXA=f0)

axslide = fig.add_subplot(223, facecolor=axcolor) #plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
sslide = Slider(axslide, 'Poloidal Surface', 1, 96, valinit=f0, valfmt='%0.0f', valstep=1.0)


def update(val):

    PolPos = sslide.val
    l.set_xdata(RadLoc.loc[:,PolPos,Attempt])
    l.set_ydata(VertLoc.loc[:,PolPos,Attempt])
    fig.canvas.draw_idle()


sslide.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    sslide.reset()
    
button.on_clicked(reset)

'''
rax = plt.axes([0.025, 0.5, 0.15, 0.15], facecolor=axcolor)
radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)


def colorfunc(label):
    l.set_color(label)
    fig.canvas.draw_idle()
radio.on_clicked(colorfunc)
'''

plt.show()