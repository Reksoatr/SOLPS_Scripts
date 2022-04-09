# -*- coding: utf-8 -*-
"""
Created on Sat Apr  9 01:39:22 2022

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
    cost = (np.sum(np.abs(y_true-y_predicted)/y_true))/100#/len(y_true)
    return cost


def Loss(exper_shot, sol_run, plot =False, ice=0, lib = 3, run_step =1):
    ius = InterpolatedUnivariateSpline(exper_shot[0], exper_shot[1])
    exp_pts = point_finder(sol_run[0],ius, y_only = True)
    loss = error(exp_pts, sol_run[1])
    if plot == True:
        y = np.abs(exp_pts-sol_run[1])/exp_pts
        plt.figure()
        plt.plot(sol_run[0],  y)
        plt.savefig(f'error_graph{ice}')
        f = open(f'/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_{lib}/Attempt_{run_step}/graph.csv', 'w')
        for i in range(len(sol_run[0])):    
            f.writelines(f'{sol_run[0][i]}   {y[i]}\n')
        f.close()
        f = open(f'/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_{lib}/Attempt_{run_step}/exp_pts', 'w')
        for i in exp_pts:    
            f.writelines(f'{i}\n')
        f.close()
        #plt.show()
    return loss



def Further_Analysis(params, exper_shot, gfilen, lib = 3, alpha =.3, run_step = 1, steps = 4):
    '''Post Step Analysis using a comparison of a given experimental shot
    to analyize loss and provided desired run for further optimization.'''
#    n = len(params)
    eq = equilibrium(gfile=gfilen)
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
    #print(exp_data)
    params_news = []
    enter = f'/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_{lib}/Attempt_{run_step}'
    os.chdir(enter)
    os.system('pwd')
    os.system('rm *.last10')
    os.system('2d_profiles')
    print(f'Attempt_{run_step}')
    try:
        Attempt = np.loadtxt('ne3da.last10')
    except:
        print('Run Failed Unexpectedly')
        return
    if len(Attempt) != 0:
        Attempt = Attempt.T
        print('length of R is: ', len(Attempt[0]))
        print(Attempt[0])
        R_sep = PsiN2R(eq, 1.0)
        print(R_sep)
        new_R = []
        for R in Attempt[0]:
            print('step')
            A = R + R_sep
            new_R.append(A)
        print('New')
        print(A)
        f = open(f'/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_{lib}/Attempt_{run_step}/sol_pts', 'w')
        for i in Attempt:    
            f.writelines(f'{i}\n')
        f.close()


if __name__ == '__main__':
    #Single_Guess(DoubleGauss, guess_init, run_step =1)
    blep = input('What Iteration is this?')
    trip = input('Is This Data Analysis?')
    blep =int(blep)
    losm = np.loadtxt('/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_3/params.txt')
    guess_init = [losm[0], losm[1],losm[2],losm[3], losm[4]]
    loss_val = losm[5]
    if  trip == 'y':
        Further_Analysis(guess_init, '/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_03/yag.txt', '/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_03/g027205.00275_efitpp', run_step = blep,alpha=loss_val)
