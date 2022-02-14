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
import matplotlib.pyplot as plt


#Richard's File Input
def Generate(trans_pts):
    '''
    Function that is used to turn the radial points into a readable
    b2.transport.inputfile

    Parameters
    ----------
    trans_pts : should be 2d point array, x coordinates being r-r_sep and
    y coordinates the diffusivity at that point

    Returns a data frame for use in the b2.transport.inputfile
    -------
    

    '''
    #J = 1
          #print(self._points)
    n = len(trans_pts)
    m = 0
    i = 1
    j = 1
    r = trans_pts
    print(' ndata(1, {0}, {1})= {2},'.format(i,j,n))
    for m in range(n):
        print(' tdata(1, {0}, {1}, {2})= {3}, tdata(2, {0}, {1}, {2})= {4},'.format(m+1,i,j,round(r[m][0],5),round(r[m][1],5)))
                
#Richard's file output
def B2TransportInputfileParser(file='b2.transport.inputfile', plot=False):
    
    Coefficients = {'1':'Particle density-driven diffusivity',
                       '2': 'Particle pressure-driven diffusivity',
                       '3': 'Ion thermal anomalous diffusivity',
                       '4': 'Electron thermal anomalous diffusivity',
                       '5': 'Poloidal-component of the anomalous ”pinch” velocity',
                       '6': 'Radial-component of the anomalous ”pinch” velocity',
                       '7': 'Anomalous viscosity',
                       '8': 'Anomalous radial electrical conductivity',
                       '9': 'Anomalous radial thermo-electric coefficient'}    
    
    with open(file) as f:
        dataList=f.readlines()
    
    Points={}
    ii=1

    while ii<len(dataList)-2: 
        ndata = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", dataList[ii])
        CoeffID = ndata[1]
        PtNo = int(ndata[3])
        XList = []
        YList = []
        for mm in range(PtNo):
            XList.append(float(re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",dataList[mm+ii+1])[4]))
            YList.append(float(re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",dataList[mm+ii+1])[9]))
        Points[CoeffID] = {'X':XList,'Y':YList}
        ii=ii+PtNo+1
        
    if plot:
        dd=len(Points.keys())
        fig1,ax1=plt.subplots(nrows=1,ncols=dd, )
        for ii, jj in enumerate(Points.keys()):
            ax1[ii].plot(Points[jj]['X'],Points[jj]['Y'])
            ax1[ii].set_title(r'Radial profile of {}'.format(Coefficients[jj]))
            ax1[ii].set_xlabel(r'$R-R_{sep}$')
            ax1[ii].set_ylabel(r'{} $[m^2/s]$'.format(Coefficients[jj]))
        
    return Points
#The beginning functions for optimization
def Trainer(x, a=1.05, b=2.5, c=.002,d=2.4,e=1,f=1):
    y = -a*(np.exp(-x**2/c)+1)-b*(x)+d
    return y
def T_Lit(x, a=0, b=1, c=3,d=0,e=.1,f=0):
    y= .5*(a+b*x**c)*(1-np.tanh((x-d)/e))+f
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



#Initial Case, for optimization algorithm, plus verification plots

# Gradient Descent Function
# Here iterations, learning_rate, stopping_threshold
# are hyperparameters that can be tuned

x = np.linspace(-.15, .1, 20)
y = Trainer(x)
gen_test = point_finder(x,Trainer)
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