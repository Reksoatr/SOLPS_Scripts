# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 13:35:33 2022

@author: james
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import re
from B2TransportParser import InputfileParser, Generate, WriteInputfile, batch_writer, R2PsiN, PsiN2R
from scipy.interpolate import InterpolatedUnivariateSpline
import os
from equilibrium import equilibrium
import csv
#The beginning functions for optimization
def Trainer(x, a=1.05, b=2.5, c=.002,d=2.4,e=1,f=1):
    y = -a*(np.exp(-x**2/c)+1)-b*(x)+d
    return y

def T_Lit(x, a=0, b=1, c=3,d=0,e=.1,f=0):
    y= .5*(a+b*x**c)*(1-np.tanh((x-d)/e))+f
    return y

def DoubleGauss(x, a=1.6, b=0.006, c=0.3,d=0.5,e=0.0007):
    '''
    Double-Gaussian Function
    a = Maximum (base) value of transport coefficient (typically 1.0)
    b = Width parameter of major gaussian b>e
    c = Minimum (trough) value of transport coefficient 0<c<a
    d = Fraction of max coefficient value where minor gaussian begins (c/a)<d<1
    e = Width parameter of minor gaussian e<b
    '''
    y = -(a-d*a)*(np.exp(-x**2/b))-(d*a-c)*(np.exp(-x**2/e))+a
    return y

#joint methods in SOLPS
#compare input D to output D
# look at individual error plot, decide if certain areas need increased weight
#gets points at the same location as true values
def point_finder(x, func, y_only = False):
    func_val = []
    for i in x:
        if y_only == False:
            temp = [i, func(i)]
        elif y_only == True:
            temp = func(i)
        func_val.append((temp))
    return np.array(func_val)


def mean_squared_error(y_true, y_predicted):
     
    # Calculating the loss or cost
    cost = np.sqrt(np.sum((y_true-y_predicted)**2)) / len(y_true)
    return cost


def Loss(exper_shot, sol_run):
    ius = InterpolatedUnivariateSpline(exper_shot[0], exper_shot[1])
    sol_pts = point_finder(sol_run[0],ius, y_only = True)
    loss = mean_squared_error(sol_run[1], sol_pts)
    return loss

def Setup(func, params, steps = 4):
    '''Sets up and runs many runs over the given parameter space, with steps
    determining how many grid points in each direction.'''
#    n = len(params)
    space = []
    for i in params:
        ticks = (i[1] - i[0])/steps
        meep = []
        for j in range(steps+1):
            meep.append(i[0] +j*ticks)
        space.append(meep)
    print(space)
    x = np.linspace(-.14, .08, 25)
    for i_ct, i in enumerate(space[0]):
        for j_ct, j in enumerate(space[1]):
            for k_ct, k in enumerate(space[2]):
#                enter = 'cd Attempt_{}{}{}'.format(i,j,k)
                diff = func(x, a = i, b= j, e=k)
                #os.system('nano b2.transport.inputfile')
                Points0 = InputfileParser(file='b2.transport.inputfile.vi')
                D_Points={'1' : np.array([x,diff])} #This is where the optimization method comes in
                Full_Points={'1':D_Points['1'],'3':Points0['3'],'4':Points0['4']}
                mkdir = 'cp -r base Attempt_{}{}{}'.format(i_ct,j_ct,k_ct)            
                os.system(mkdir)
                WriteInputfile(file='/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_03/Attempt_{}{}{}/b2.transport.inputfile'.format(i_ct,j_ct,k_ct),points=Full_Points)
                path_name = 'cd /sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_03/Attempt_{}{}{}'.format(i_ct,j_ct,k_ct) #finish adding mk0
                batch_writer(path_name, i_ct, j_ct, k_ct)
                os.system('cp batch_use  /sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_03/Attempt_{}{}{}/batch'.format(i_ct,j_ct,k_ct))
                batch_run = 'qsub /sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_03/Attempt_{}{}{}/batch'.format(i_ct,j_ct,k_ct)
                os.system(batch_run)
                os.system('cd ../')
                    
def Loss_Analysis(params, exper_shot, gfilen, run_step = 1, steps = 1):
    '''Post Step Analysis using a comparison of a given experimental shot
    to analyize loss and provided desired run for further optimization.'''
#    n = len(params)
    space = []
    loss_pts = []
    eq = equilibrium(gfile=gfilen)
    for i in params:
        ticks = (i[1] - i[0])/steps
        meep = []
        for j in range(steps+1):
            meep.append(i[0] +j*ticks)
        space.append(meep)                    
    STARTING = .75
    ENDING = 1.03
    exp_data = np.loadtxt(exper_shot, usecols = (0,1))
    exp_new = []
    for R in exp_data:
        if R[0] > STARTING:
            if R[0] < ENDING:
                exp_new.append(R)
    exp_Data = np.array(exp_new)
    exp_data = exp_Data.T
    print(exp_data)
    tick = 0
    for i_ct, i in enumerate(space[0]):
        for j_ct, j in enumerate(space[1]):
            for k_ct, k in enumerate(space[2]):
                enter = f'/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_03/Attempt_{i_ct}{j_ct}{k_ct}_mk{run_step}' 
                os.chdir(enter)
                os.system('pwd')
                os.system('rm *.last10')
                os.system('2d_profiles')
                print('Attempt_{}{}{}'.format(i_ct,j_ct,k_ct))
                try:
                    Attempt = np.loadtxt('ne3da.last10')
                except:
                    continue
                if len(Attempt) != 0:
                    Attempt = Attempt.T
                    R_sep = PsiN2R(eq, 1.0)
                    for R in Attempt[0]:
                        R = R2PsiN(eq,R+R_sep)
                    l = Loss(exp_data, Attempt)
                    loss_pts.append([l,i,j,k,tick])
                tick += 1
    b = np.amin(loss_pts, axis = 0)
    for i in loss_pts:
        if b[0] == i[0]: 
            b_new = i 
            print('Minimum loss is at:')
            print(i)
    params_new = []
    b_star = [b_new[1], b_new[2], b_new[3]]
    params_new.append(b_star)
    temp = int(b_new[4])
    loss_ptsb = np.delete(loss_pts, temp,0)
    b1 = np.amin(loss_ptsb, axis = 0)
    for i in loss_ptsb:
        if b1[0] == i[0]: 
            b_new = i 
    b_star1 = [b_new[1], b_new[2], b_new[3]]
    params_new.append(b_star1)
    params_new = np.array(params_new)
    params_news =  params_new.T
    print(params_news)
    return params_news
#add last10 notes to look at different last10 file, check if last10 files need deleted
#use mv command rm b2mn.prt  
#ls -al
    with open('loss_over_iteration.csv', 'a', encoding='UTF8') as f:
        writer = csv.writer(f)
        writer.writerows([b])
#need to add error/iteration graph,, all parameters need to be optimized
def Loss_Graph(csv):
    y = np.loadtxt(csv, usecols = 0)
    x= range(0, len(y))
    fig, axs = plt.subplots(1,1,dpi = 200)
    axs.plot(x,y)
    axs.set_xlabel('Iterations')
    axs.set_ylabel('Loss from Error')
    

def Further_Steps(func, params, alpha = .2, run_step=2, Post_Analysis = True, exper_shot = None, gfilen = None):
    space = []
    if Post_Analysis == False:
        params = Loss_Analysis(params, exper_shot, gfilen)
    for i in params:
        step = alpha*(i[1]-i[0])
        i[1] = step
        i.append(2*alpha*(i[1]-i[0]))
    space = params
    x = np.linspace(-.14, .08, 25)
    os.system('cp base/b2fstate base/b2fstati')
    for i_ct, i in enumerate(space[0]):
        for j_ct, j in enumerate(space[1]):
            for k_ct, k in enumerate(space[2]):
                diff = func(x, a = i, b= j, e=k)
                Points0 = InputfileParser(file='b2.transport.inputfile.vi')
                D_Points={'1' : np.array([x,diff])} #This is where the optimization method comes in
                Full_Points={'1':D_Points['1'],'3':Points0['3'],'4':Points0['4']}
                mkdir = f'cp -r base Attempt_{i_ct}{j_ct}{k_ct}_mk{run_step}'         
                os.system(mkdir)
                WriteInputfile(file=f'/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_03/Attempt_{i_ct}{j_ct}{k_ct}_mk{run_step}/b2.transport.inputfile',points=Full_Points)
                path_name = f'cd /sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_03/Attempt_{i_ct}{j_ct}{k_ct}_mk{run_step}'
                batch_writer(path_name, i_ct, j_ct, k_ct)
                os.system(f'cp batch_use  /sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_03/Attempt_{i_ct}{j_ct}{k_ct}_mk{run_step}/batch')
                batch_run = f'qsub /sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_03/Attempt_{i_ct}{j_ct}{k_ct}_mk{run_step}/batch'
                os.system(batch_run)

   #check errors if they are going down/flat space for convergence check initial run
def Single_Guess(func, guess, alpha = .2, run_step=2, lib=5, Post_Analysis = False, exper_shot = None, gfilen = None):
    if Post_Analysis ==False:
        x = np.linspace(-.14, .08, 15)
        y = np.linspace(.003, .006, 10)
        for i_ct, i in enumerate(y):
            diff = func(x, a = guess[0], b= i, c=guess[2], d = guess[3], e = guess[4])
            Points0 = InputfileParser(file='b2.transport.inputfile.vi')
            D_Points={'1' : np.array([x,diff])} #This is where the optimization method comes in
            Full_Points={'1':D_Points['1'],'3':Points0['3'],'4':Points0['4']}
            mkdir = f'cp -r base Attempt_mk{run_step}{i_ct}'         
            os.system(mkdir)
            WriteInputfile(file=f'/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_0{lib}/Attempt_mk{run_step}{i_ct}/b2.transport.inputfile',points=Full_Points)
            path_name = f'cd /sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_0{lib}/Attempt_mk{run_step}{i_ct}'
            batch_writer(path_name, run_step, 0, 0)
            os.system(f'cp batch_use  /sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_0{lib}/Attempt_mk{run_step}{i_ct}/batch')
            batch_run = f'qsub /sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_0{lib}/Attempt_mk{run_step}{i_ct}/batch'
            os.system(batch_run)



    
MAST_params = [[1,2],
          [.002,.0075],
          [.0005,.003]]

MAST_params_it = [[1.000e+00, 1.250e+00],
                  [7.500e-03, 6.125e-03],
                  [1.125e-03, 1.750e-03]]
guess_init=[1.6, 0.005, 0.25,0.5,0.0007]
#Initial Case, for optimization algorithm, plus verification plots

# Gradient Descent Function
# Here iterations, learning_rate, stopping_threshold
# are hyperparameters that can be tuned

if __name__ == '__main__':
    Single_Guess(DoubleGauss, guess_init, run_step = 1)
    '''
    initializing = input('Is this before your first run? (y or n)')
    if initializing == 'y':
        Setup(DoubleGauss, MAST_params)
    if initializing == 'n':
        blep = input('What Iteration is this?')
        blep =int(blep)
        data_analysis = input('Is this data analysis after a run? (y or n)')
        if data_analysis == 'y':
            Loss_Analysis(MAST_params, '/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_03/yag.txt', 'g027205.00275_efitpp', run_step = blep)
        elif data_analysis == 'n':

x = np.linspace(-.25,.20)
y = DoubleGauss(x, c=.3)
Points = InputfileParser('b2.transport.inputfile.dblgausstest')
test = Points['1'][0]
test_y = Points['1'][1]
fig, axs = plt.subplots(1,1, dpi = 200)
axs.plot(test, test_y, color = 'b', label = 'Original')
axs.plot(x, y, color = 'g', label = 'Training')
#axs.plot(x, y_Lit, color = 'm', label = 'Literature')
axs.legend()
'''
