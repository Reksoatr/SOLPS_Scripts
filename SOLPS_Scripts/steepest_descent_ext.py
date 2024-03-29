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
def Trainer(x, a=1.5, b=3, c=.004,d=.3):
    y = -a*(np.exp(-(x-.005)**2/c))-b*(x)*np.exp(-x**2/.01)+a+d
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
    y = -(a-d*a)*(np.exp(-(x-.01)**2/b))-(d*a-c)*(np.exp(-(x-.01)**2/e))+a
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


def error(y_true, y_predicted):
     
    # Calculating the loss or cost
    cost = (np.sum(np.abs(y_true-y_predicted)/y_true))/100000#/len(y_true)
    return cost


def Loss(exper_shot, sol_run, plot =False, ice=0):
    ius = InterpolatedUnivariateSpline(exper_shot[0], exper_shot[1])
    exp_pts = point_finder(sol_run[0],ius, y_only = True)
    loss = error(sol_run[1], exp_pts)
    if plot == True:
        plt.figure()
        plt.plot(sol_run[0], np.abs(exp_pts-sol_run[0])/exp_pts )
        plt.savefig(f'error_graph{ice}')
        plt.show()
    return loss


#add last10 notes to look at different last10 file, check if last10 files need deleted
#use mv command rm b2mn.prt  
#ls -al
        
        
def Further_Analysis(params, exper_shot, gfilen, lib = 22, alpha =.3, run_step = 1, steps = 4, learn = .3):
    '''Post Step Analysis using a comparison of a given experimental shot
    to analyize loss and provided desired run for further optimization.'''
#    n = len(params)
    eq = equilibrium(gfile=gfilen)
    space = []
    loss_pts = []
    for i in params:
        l =[]
        step_0 = learn*alpha*i+i
        step_1 = -1*learn*alpha*i+i
        l.append(i)
        l.append(step_0)
        l.append(step_1)
        space.append(l)
    space = np.array(space).T
    print(space)
    STARTING = .62
    ENDING = 1.06
    exp_data = np.loadtxt(exper_shot, usecols = (0,1))
    exp_new = []
    for R in exp_data:
        if R[0] > STARTING:
            if R[0] < ENDING:
                exp_new.append(R)
    exp_Data = np.array(exp_new)
    exp_data = exp_Data.T
    #print(exp_data)
    tick = 0
    for i_ct, i in enumerate(space):
        enter = f'/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_{lib}/Attempt_{i_ct}_mk{run_step}'
        try:
            os.chdir(enter)
        except:
            continue
        os.system('pwd')
        os.system('rm *.last10')
        os.system('2d_profiles')
        print(f'Attempt_{i_ct}')
        try:
            Attempt = np.loadtxt('ne3da.last10')
        except:
            continue
        if len(Attempt) != 0:
            Attempt = Attempt.T
            R_sep = PsiN2R(eq, 1.0)
            for R in Attempt[0]:
                R = R2PsiN(eq,R+R_sep)
            l = Loss(exp_data, Attempt, plot=True,ice=i_ct)
            loss_pts.append([l,i[0],i[1], i[2], i[3], i_ct])
    b = np.amin(loss_pts, axis = 0)
    print('initial guess is:', loss_pts[0])
    print('Difference in loss is:', loss_pts[0][0]-b[0])
    new_step = 10*loss_pts[0][0]-b[0]
    if b[0] == loss_pts[0][0]:
        params_news = [loss_pts[0][1], loss_pts[0][2], loss_pts[0][3], loss_pts[0][4]]
        new_step = learn/2
        print('guess too far')
    elif b[0] == loss_pts[1][0]:
        params_news = [loss_pts[1][1], loss_pts[1][2], loss_pts[1][3], loss_pts[1][4]]
        print('go right')
    elif b[0] == loss_pts[2][0]:
        params_news = [loss_pts[2][1], loss_pts[2][2], loss_pts[2][3], loss_pts[2][4]]
        print('go left')
    f = open(f'error{run_step}.csv', 'w')
    f.writelines(f'{run_step}   {b}')
    f.close()
    new_loss = b[0]
    print(params_news)
    os.chdir(f'/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_{lib}/')
    return params_news, new_loss, new_step

#ls -al
def Loss_Graph(cst):
    y = np.loadtxt(cst, usecols = 0)
    x= range(0, len(y))
    fig, axs = plt.subplots(1,1,dpi = 200)
    axs.plot(x,y)
    axs.set_xlabel('Iterations')
    axs.set_ylabel('Loss from Error')


def Further_Steps(func, params, alpha = .3, run_step=2, lib = 22,Post_Analysis = True, exper_shot = None, gfilen = None, learn = .3):
    space = []
    for i in params:
        l =[]
        step_0 = learn*alpha*i+i
        step_1 = -1*learn*alpha*i+i
        l.append(i)
        l.append(step_0)
        l.append(step_1)
        space.append(l)
    space = np.array(space).T
    print(space)
    x_1 = np.linspace(-.12, -.03, 5)
    x_2 = np.linspace(-.02, .02, 10)
    x = np.append(x_1, x_2)
    #os.system('cp base/b2fstate base/b2fstati')
    for i_ct, i in enumerate(space):
        diff = func(x, a = i[0], b= i[1], c = i[2],d = i[3])
        Points0 = InputfileParser(file='b2.transport.inputfile.vi')
        D_Points={'1' : np.array([x,diff])} #This is where the optimization method comes in
        Full_Points={'1':D_Points['1'],'3':Points0['3'],'4':Points0['4']}
        mkdir = f'cp -r base Attempt_{i_ct}_mk{run_step}'         
        os.system(mkdir)
        WriteInputfile(file=f'/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_{lib}/Attempt_{i_ct}_mk{run_step}/b2.transport.inputfile',points=Full_Points)
    print(space)




    
MAST_params = [[1,2],
          [.002,.0075],
          [.0005,.003]]

MAST_params_it = [[2.125000e+00, 2.375000e+00],
                  [2.759375e-03, 2.959375e-03],
                  [3.000e-04, 3.25000e-04]]

loss_val = .93
guess_init=[1.90831261, 3.81662522, 0.00508883, 0.38166253]
learning_rate=.3#[2.25, 0.002759375, 0.0003]
#Initial Case, for optimization algorithm, plus verification plots

# Gradient Descent Function
# Here iterations, learning_rate, stopping_threshold
# are hyperparameters that can be tuned

if __name__ == '__main__':
    #Single_Guess(DoubleGauss, guess_init, run_step =1)
    blep = input('What Iteration is this?')
    blep =int(blep)
    losm = np.loadtxt('params.txt')
    guess_init = [losm[0], losm[1],losm[2],losm[3]]
    loss_val = losm[4]
    learning_rate = losm[5]
    h=1

    guess_init, loss_val, learning_rate = Further_Analysis(guess_init, '/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_03/yag.txt', 'g027205.00275_efitpp', run_step = blep,alpha=loss_val, learn=learning_rate)
    f = open('params.txt', 'w')
    f.writelines(f'{guess_init[0]} {guess_init[1]} {guess_init[2]} {guess_init[3]} ')
    f.writelines(f'{loss_val} {learning_rate}')
    f.close()
    blep += 1
    Further_Steps(Trainer, guess_init, run_step=blep, alpha = loss_val, learn = learning_rate)        

'''
data = [[1, 1.018],
        [2,	1.0182801375187802-0.029898925908849794],
        [3, .931],
        [4, 0.9054375472903425],
        [5,0.8608029514399158],
        [6,0.],
        [7,0.7912938107323133],
        [8,0.7640928445496504],
        [9, 0.7420399341877436],
        [10, 0.7124840953627268],
        [11, 0.7122870767295693],
        [12, 0.7122870767295693-0.01340096684815062]]
data = np.array(data).T
plt.figure()
plt.plot(data[0], data[1], '-')
plt.title('Error Graph')
plt.x_label('iterations')
plt.y_label('error (normalized)')

x = np.linspace(-.14,.1)
y = Trainer(x,)
Points = InputfileParser('b2.transport.inputfile.dblgausstest')


test = Points['1'][0]
test_y = Points['1'][1]
fig, axs = plt.subplots(1,1, dpi = 200)
axs.plot(test, test_y, color = 'b', label = 'Original')
axs.plot(x, y, color = 'g', label = 'Training')
#axs.plot(x, y_Lit, color = 'm', label = 'Literature')
axs.legend()

'''