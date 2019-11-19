# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 16:15:42 2019

@author: rmreksoatmodjo
"""

from VesselPlotterNew import SOLPSPLOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

Shot = '12'
Attempt = '65'

NeuDen = SOLPSPLOT(Shot,Attempt,'NeuDen')
JXA = NeuDen.KW['JXA']
JXI = NeuDen.KW['JXI']
RadLoc = NeuDen.RadCoords['RadLoc']
VertLoc = NeuDen.RadCoords['VertLoc']


fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)

NeuDen.Contour('NeuDen',LOG10=1,AX=ax, Markers=False)

#t = np.arange(0.0, 1.0, 0.001)

#a0 = 5
f0 = JXA
#delta_f = 1.0
#s = a0 * np.sin(2 * np.pi * f0 * t)
l, = ax.plot(RadLoc.loc[:,f0,Attempt],VertLoc.loc[:,f0,Attempt],color='Red',linewidth=3)
ax.margins(x=0)

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
#axamp = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

sfreq = Slider(axfreq, 'Poloidal Surface', 24.0, 71.0, valinit=f0, valfmt='%0.0f', valstep=1.0)
#samp = Slider(axamp, 'Amp', 0.1, 10.0, valinit=a0)


def update(val):
    #amp = samp.val
    freq = sfreq.val
    l.set_xdata(RadLoc.loc[:,freq,Attempt])
    l.set_ydata(VertLoc.loc[:,freq,Attempt])
    #ax.plot(RadLoc.loc[:,freq,Attempt],VertLoc.loc[:,freq,Attempt],color='Red',linewidth=3)
    fig.canvas.draw_idle()


sfreq.on_changed(update)
#samp.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    sfreq.reset()
    #samp.reset()
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