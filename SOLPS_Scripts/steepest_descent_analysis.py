# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 18:30:32 2022

@author: Jameson Crouse
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
    cost = np.sum((np.abs(y_true-y_predicted)/y_true))/len(y_true)
    return cost


def Loss(exper_shot, sol_run, plot =False, ice=0, lib = 4, run_step =1):
    ius = InterpolatedUnivariateSpline(exper_shot[0], exper_shot[1])
    exp_pts = point_finder(sol_run[0],ius, y_only = True)
    loss = error(exp_pts, sol_run[1])
    if plot == True:
        y = (exp_pts-sol_run[1])/exp_pts
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
        f = open(f'/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_{lib}/Attempt_{run_step}/sol_pts', 'w')
        for i in sol_run[1]:    
            f.writelines(f'{i}\n')
        f.close()
        #plt.show()
    return loss



def Further_Analysis(params, exper_shot, gfilen, lib = 4, alpha =.3, run_step = 1, steps = 4, dire = 1):
    '''Post Step Analysis using a comparison of a given experimental shot
    to analyize loss and provided desired run for further optimization.'''
#    n = len(params)
    eq = equilibrium(gfile=gfilen)
    STARTING = .75
    ENDING = 1.02
    exp_data = np.loadtxt(exper_shot, usecols = (0,1))
    exp_new = []
    for R in exp_data:
        if R[0] > STARTING-.1:
            if R[0] < ENDING+.04:
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
        R_sep = PsiN2R(eq, 1.0)
        new_R = []
        for R in Attempt[0]:
            A = R + R_sep
            B= R2PsiN(eq,A)
            new_R.append(float(B))
        Attempt[0]=new_R
        Attempt = Attempt.T
        Att_new = []
        for R in Attempt:
            if R[0] > STARTING:
                if R[0] < ENDING:
                    Att_new.append(R)
        attempt = np.array(Att_new)
        Attempt = attempt.T
        for R in Attempt:
            if R[0] > STARTING:
                if R[0] < ENDING:
                    Att_new.append(R)
        attempt = np.array(Att_new)
        Attempt = attempt.T
        l = Loss(exp_data, Attempt, plot=True, ice = run_step, lib = lib, run_step=run_step)
        print(l)
        b = (alpha-l)*10
        while np.abs(b) < .1:
            b = b*10
        if run_step == 1:
            b= .4
            
    print('Difference in loss is:', b)
    if b==0:
        params_news = params
        print('guess too far')
    elif b<0:
        new_dire = -1*dire
        if dire == 1 or 1.0:
            for j in params:
                params_news.append(float(j)+float(j)*b)
        elif dire == -1 or -1.0:
            for j in params:
                params_news.append(float(j)-float(j)*b)
    elif b>0:
        new_dire = dire
        if dire == 1 or 1.0:
            for j in params:
                params_news.append(float(j)+float(j)*b)
        elif dire == -1 or -1.0:
            for j in params:
                params_news.append(float(j)-float(j)*b)
    if run_step == 1:
        new_dire = np.sign(b)
    f = open(f'/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_{lib}/error.csv', 'a')
    f.writelines(f'\n{run_step}   {l}')
    f.close()
    new_loss = l
    print(params_news)
    os.chdir(f'/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_{lib}')
    return params_news, new_loss, new_dire

def Further_Steps(func, params, run_step=1, lib = 4,Post_Analysis = True, exper_shot = None, gfilen = None):
    print(params)
    x_1 = np.linspace(-.12, -.03, 5)
    x_2 = np.linspace(-.02, .02, 10)
    x = np.append(x_1, x_2)
    #os.system('cp base/b2fstate base/b2fstati')
    diff = func(x, a = params[0], b= params[1], c = params[2],d = params[3], e = params[4])
    Points0 = InputfileParser(file='b2.transport.inputfile.vi')
    D_Points={'1' : np.array([x,diff])} #This is where the optimization method comes in
    Full_Points={'1':D_Points['1'],'3':Points0['3'],'4':Points0['4']}
    mkdir = f'cp -r /sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_{lib}/base /sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_{lib}/Attempt_{run_step}'         
    os.system(mkdir)
    WriteInputfile(file=f'/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_{lib}/Attempt_{run_step}/b2.transport.inputfile',points=Full_Points)

if __name__ == '__main__':
    #Single_Guess(DoubleGauss, guess_init, run_step =1)
    blep = input('What Iteration is this?')
    trip = input('Is This Data Analysis?')
    libr = int(input('What Directory?'))
    blep =int(blep)
    losm = np.loadtxt(f'/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_{libr}/params.txt')
    guess_init = [losm[0], losm[1],losm[2],losm[3], losm[4]]
    loss_val = losm[5]
    direction = losm[6]
    if  trip == 'y':
        guess_init, loss_val, direction = Further_Analysis(guess_init, '/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_03/yag.txt', 'g027205.00275_efitpp', run_step = blep,alpha=loss_val,lib=libr, dire =direction)
        f = open('params.txt', 'w')
        f.writelines(f'{guess_init[0]} {guess_init[1]} {guess_init[2]} {guess_init[3]} {guess_init[4]}')
        f.writelines(f' {loss_val} {direction}')
        f.close()
    elif trip == 'n':
        Further_Steps(DoubleGauss, guess_init, run_step=blep, lib = libr)  