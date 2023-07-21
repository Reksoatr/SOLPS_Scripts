# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 14:59:41 2019

@author: rmreksoatmodjo

Collection of general Tools to perform oft-repeated SOLPS data analyis and post-processing tasks
"""
import os
import numpy as np
import paramiko
from scipy import stats

def SET_WDIR(BASEDRT,TOPDRT): #Function to set correct Working Directory Path depending on which machine is in use
    if os.environ['OS'] == 'Windows_NT':
        if os.environ['USERNAME'] == 'rmreksoatmodjo':
            BASEDRT = r"C:/Users/rmreksoatmodjo/Desktop/WMGDrive/College of William and Mary/Research/SOLPS Stuff/SOLPS_2D_prof/"
            TOPDRT = r"C:/Users/rmreksoatmodjo/Desktop/WMGDrive/College of William and Mary/Research/SOLPS Stuff/"
        elif os.environ['USERNAME'] == '18313':
            BASEDRT = r"G:/My Drive/College of William and Mary/Research/SOLPS Stuff/SOLPS_2D_prof/"
            TOPDRT = r"G:/My Drive/College of William and Mary/Research/SOLPS Stuff/"
        elif os.environ['USERNAME'] == 'Richard':
            BASEDRT = r"C:/Users/Richard/WMGDrive/College of William and Mary/Research/SOLPS Stuff/SOLPS_2D_prof/"
            TOPDRT = r"C:/Users/Richard/WMGDrive/College of William and Mary/Research/SOLPS Stuff/"
        elif os.environ['USERNAME'] == 'Yi-Cheng':
            BASEDRT = r"C:/Users/Yi-Cheng/Documents/SOLPS_Data/Simulation_Data"
            TOPDRT = r"C:/Users/Yi-Cheng/Documents/SOLPS_Data/Experimental_Data"
    else:
        BASEDRT=BASEDRT
        TOPDRT=TOPDRT
    
    return BASEDRT, TOPDRT

def TANH(r,r0,h,d,b,m):
    return b+(h/2)*(np.tanh((r0-r)/d)+1) + m*(r0-r-d)*np.heaviside(r0-r-d, 1)

def EXPFIT(x,A,l):  #Removed vertical displacement variable B; seemed to cause 'overfitting'
    return A*np.exp(l*x)

def R2PsiN(GF,R,Z=0):
    '''Uses equilibrium to convert from R to PsiN
        Must provide gfile (GF) loaded as an equilibrium object
        Default assumes Z=0'''
    PP=GF.psiN(R,Z).flatten()
    
    return PP

def PsiN2R(GF,psin,Z=0):
    '''uses equilibrium to convert from PsiN to R
       Must provide gfile (GF) loaded as an equilibrium object
       Default assumes Z=0'''
    Rlfs=[i for i in GF.R if i>GF.axis.r]
    RR=np.interp(psin,GF.psiN(Rlfs,Z).flatten(),Rlfs)
    
    return RR

def WALL_INTERSECT(C0,C1,r):
    '''
    Calculate the intersection between a line or set of lines specified by points
    C0 and C1, and a circle of radius r centered at the origin
    
    C0 : dict of floats or arrays
        Starting coordinates of line(s), in the format C0={'X': ...,'Y': ...}
    C1 : dict of floats or arrays
        Ending coordinates of line(s), in the format C1={'X': ...,'Y': ...}
    r : float
        Radius of circle    
    '''
    dx=C1['X']-C0['X']
    dy=C1['Y']-C0['Y']
    dr=np.sqrt(dx**2+dy**2)
    D=C0['X']*C1['Y']-C1['X']*C0['Y']
    delta=(r**2)*(dr**2)-D**2
    
    X1=(D*dy+np.sign(dy)*dx*np.sqrt(delta))/dr**2
    Y1=(-D*dx+np.abs(dy)*np.sqrt(delta))/dr**2
    X2=(D*dy-np.sign(dy)*dx*np.sqrt(delta))/dr**2
    Y2=(-D*dx-np.abs(dy)*np.sqrt(delta))/dr**2
        
    P1={'X':X1,'Y':Y1}
    P2={'X':X2,'Y':Y2}
    
    return P1,P2  

def gaussian_shading(ax, x, y, y_unc, c='k', min_val=0.0):
    ''' Plot profile with uncertainties displayed as a shading whose color intensity represents a 
    gaussian PDF.
    Adapted from Francesco's method in lyman_single.py
    '''
    norm_val = stats.norm.pdf(0)
    
    num=50  # discrete number of shades    
    for ij in np.arange(num):

        # below mean
        ax.fill_between(x,
                        np.maximum(y - 5*y_unc*(ij-1)/num, min_val),
                        np.maximum(y - 5*y_unc*ij/num, min_val),
                        alpha=stats.norm.pdf(5*ij/num)/norm_val,
                        linewidth=0.0,
                        color=c)

    # start looping from 2 to avoid overshading the same region
    for ij in np.arange(2,num):
        # above mean
        ax.fill_between(x, 
                        y + 5*y_unc*(ij-1.)/num,
                        y + 5*y_unc*ij/num,
                        alpha=stats.norm.pdf(5*ij/num)/norm_val,
                        linewidth=0.0,
                        color=c)

def ErrorQuant(exp_data, model_data, exp_unc=None, name=None):
    '''Calculate a collection of Goodness-of-Fit Error Quantification metrics 
    for a set of experiemntal data and corresponding model prediction data'''
    
    if len(exp_data) == len(model_data):
        N = len(exp_data)
        err = exp_data - model_data
        exp_mean = np.mean(exp_data)
        exp_var = exp_data - exp_mean
        
        MAE = np.sum(np.abs(err))/N
        RMSE = np.sqrt(np.sum(err**2)/N)
        
        RAE = np.sum(np.abs(err))/np.sum(np.abs(exp_var))
        RSE = np.sqrt(np.sum(err**2)/np.sum(exp_var**2))
        
        if exp_unc is not None and len(exp_unc) == len(exp_data):
            norm_res = err/exp_unc
            SER = np.sqrt(np.sum(norm_res**2)/N)
            
            return {'MAE':MAE, 'RMSE':RMSE, 'RAE':RAE, 'RSE':RSE, 'norm_res':norm_res, 'SER':SER, 'Quantity':name}
            
        else:
            
            return {'MAE':MAE, 'RMSE':RMSE, 'RAE':RAE, 'RSE':RSE, 'Quantity':name}
        
    else:
        print("Experimental data array and model data array must be the same length")
        quit()            


def JumpConnect(host, user, ssh_home, jumphost,port=22):
    client = paramiko.SSHClient()
    client.load_system_host_keys(filename='{}known_hosts'.format(ssh_home))
    if jumphost:
        jh_client = paramiko.SSHClient()
        jh_client.load_system_host_keys(filename='{}known_hosts'.format(ssh_home))
        jh_client.connect(jumphost, username=user, key_filename='{}id_rsa'.format(ssh_home))
        sock = jh_client.get_transport().open_channel(
                'direct-tcpip', (host, 22), ('', 0)
                )
        kwargs = dict(
            hostname=host,
            port=port,
            username=user,
            key_filename='{}id_rsa'.format(ssh_home),
            sock=sock,
        )
    else:
        kwargs = dict(
            hostname=host,
            port=port,
            username=user,
            key_filename='{}id_rsa'.format(ssh_home),
        )
    client.connect(**kwargs)
    return client

def OpenRemoteFile(filepath, readtype='r',
                   host='bora.sciclone.wm.edu', 
                   user='rmreksoatmodjo', 
                   ssh_home='C:/cygwin64/home/18313/.ssh/', 
                   jumphost=None):
    if jumphost:
        client=JumpConnect(host, user, ssh_home, jumphost)
    else:
        client = paramiko.SSHClient()
        client.load_system_host_keys(filename='{}known_hosts'.format(ssh_home))
        client.connect(host,username=user,key_filename='{}id_rsa'.format(ssh_home))
    sftp_client=client.open_sftp()
    file=sftp_client.file(filepath, readtype)
    return file

def SSH_config(server):
    if os.name == 'nt':
        if os.environ['USERNAME'] == '18313':
            if server == 'bora':
                Kwargs = dict(
                    host='bora.sciclone.wm.edu',
                    user='rmreksoatmodjo',
                    ssh_home='C:/cygwin64/home/18313/.ssh/',
                    jumphost='stat.wm.edu')
            elif server == 'cmod':
                Kwargs = dict(
                    host='mfews08.psfc.mit.edu',
                    user='reksoatr',
                    ssh_home='C:/cygwin64/home/18313/.ssh/',
                    jumphost=None,
                    port=9224)
    return Kwargs
            
'''
class FORT44(object): #Class of methods used to parse and organize SOLPS fort.44 data
    
    def __init__(self,Shot,Attempt,Parameter,**kwargs):
        
        self._reset_object()
        
        self.Shot=Shot
        
        self.Attempt=Attempt
        
        self.Parameter=Parameter
        
    def _reset_object(self):
        self.Shot=None
        self.Attempts=None
        self.Parameter=None
        
    def _OpenFort44(self):
        
    def _FindParameter(self):
        
    def PlotWalls(self):
'''        
        