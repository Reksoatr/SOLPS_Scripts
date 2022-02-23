# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 13:35:33 2022

@author: james
"""

import numpy as np
import scipy.optimize as spop
import scipy as sp
import matplotlib.pyplot as plt
import re
from B2TransportParser import InputfileParser, Generate, WriteInputfile, replace_line
from scipy.interpolate import InterpolatedUnivariateSpline
import os
#The beginning functions for optimization
def Trainer(x, a=1.05, b=2.5, c=.002,d=2.4,e=1,f=1):
    y = -a*(np.exp(-x**2/c)+1)-b*(x)+d
    return y

def T_Lit(x, a=0, b=1, c=3,d=0,e=.1,f=0):
    y= .5*(a+b*x**c)*(1-np.tanh((x-d)/e))+f
    return y

def DoubleGauss(x, a=1.0, b=0.005, c=0.1,d=0.5,e=0.0001):
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
def point_finder(x, func):
    func_val = []
    for i in x:
        temp = [i, func(i)]

        func_val.append((temp))
    return np.array(func_val)


def mean_squared_error(y_true, y_predicted, ):
     
    # Calculating the loss or cost
    cost = np.sum((y_true-y_predicted)**2) / len(y_true)
    return cost


def Loss(func, exper_shot, sol_run, iterations = 100):
    ius = InterpolatedUnivariateSpline(sol_run[0], sol_run[1])
    sol_pts = point_finder(exper_shot[0],ius)
    loss = mean_squared_error(exper_shot[1], sol_pts)
    return loss

def Opt_BruteForce(func, params, exper_shot, sol_run,attempt_path,initial = True, points = 50, steps = 4):
    '''Algorithm to create multiple runs across parameter space'''
#    n = len(params)
    space = []
    for i in params:
        ticks = (i[1] - i[0])/steps
        meep = []
        for j in range(steps+1):
            meep.append(i[0] +j*ticks)
        space.append(meep)
    for i in space[0]:
        for j in space[1]:
            for k in space[2]:
                x = np.linspace(-.15, .1, 50)
                diff = func(x)
              #  os.system('nano b2.transport.inputfile')
                Points0 = InputfileParser(file='b2.transport.inputfile.V0')
                D_Points={'1' : [[x],[diff]]} #This is where the optimization method comes in
                Full_Points={'1':D_Points['1'],'3':Points0['3'],'4':Points0['4']}
                mkdir = 'cp -r base Attempt_{}{}{}'.format(i,j,k)            
                enter = 'cd Attempt_{}{}{}'.format(i,j,k)
                os.system(mkdir)
                os.system(enter)
                WriteInputfile(file='b2.transport.inputfile',points=Full_Points)
                path_name = '/sciclone/scr20/gjcrouse/SOLPS/runs/OPT_TEST_01'
                Attempt = 'Attempt_{}{}{}'.format(i,j,k)
                #replaces the name and directory lines
                replace_line('batch', 4, Attempt)
                replace_line('batch', 11, path_name)
                batch_run = 'qsub batch_{}{}{}'.format(i,j,k)
                os.system(batch_run)
                
#need to add error/iteration graph

#Initial Case, for optimization algorithm, plus verification plots

# Gradient Descent Function
# Here iterations, learning_rate, stopping_threshold
# are hyperparameters that can be tuned
if __name__ == '__main__':
    x = np.linspace(-.15, .15, 20)
    y = DoubleGauss(x)
    gen_test = point_finder(x,DoubleGauss)
    Generate(gen_test)
'''
y_Lit = T_Lit(x)
Points = B2TransportInputfileParser()
test = Points['1']['X']
test_y = Points['1']['Y']
pointy = point_finder(test,Trainer)
#print(test)
fig, axs = plt.subplots(1,1, dpi = 200)
axs.plot(test, test_y, color = 'b', label = 'Original')
axs.plot(x, y, color = 'g', label = 'Training')
#axs.plot(x, y_Lit, color = 'm', label = 'Literature')
axs.legend()

#Testing
#x=  np.linspace(-.2,.1)
#plt.plot(x,Trainer(x))



def gradient_descent(x, y, iterations = 1000, learning_rate = 0.0001,
                     stopping_threshold = 1e-6):
     
   
    # Estimation of optimal parameters
    for i in range(iterations):
         
        # Making predictions
        y_predicted = (current_weight * x) + current_bias
         
        # Calculationg the current cost
        current_cost = mean_squared_error(y, y_predicted)
 
        # If the change in cost is less than or equal to
        # stopping_threshold we stop the gradient descent
        if previous_cost and abs(previous_cost-current_cost)<=stopping_threshold:
            break
         
        previous_cost = current_cost
 
        costs.append(current_cost)
        weights.append(current_weight)
         
        # Calculating the gradients
        weight_derivative = -(2/n) * sum(x * (y-y_predicted))
        bias_derivative = -(2/n) * sum(y-y_predicted)
         
        # Updating weights and bias
        current_weight = current_weight - (learning_rate * weight_derivative)
        current_bias = current_bias - (learning_rate * bias_derivative)
                 
        # Printing the parameters for each 1000th iteration
        print(f"Iteration {i+1}: Cost {current_cost}, Weight \
        {current_weight}, Bias {current_bias}")'''