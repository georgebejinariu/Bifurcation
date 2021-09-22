# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 15:54:16 2020

@author: Asus
"""

# We study the system: du/dt = alpha/(1+v^beta) - u
#                      dv/dt = alpha/(1+u^beta) - v

from functools import partial 
from collections import defaultdict 
import numpy as np # Numerical computing library
import matplotlib.pyplot as plt # Plotting library
import scipy.integrate #Integration library
from mpl_toolkits.mplot3d import axes3d #Used for the 3d bifurcation plot
import matplotlib.patches as mpatches #used to write custom legends


# a Cauchy problem ca be solved on the interval [0,T] with scipy.integrage.odeint(f,y0,time)

#_____________________________________________________________
# - Simulate this system using scipy.integrate.odeint
# - Draw the trajectories using matplotlib.pyplot.plot

# We will look at those set of parameters
scenarios = [{'alpha':1, 'beta':2}, 
             {'alpha':1, 'beta':10}]

# On this timespan
time = np.linspace(0, 20, 1000)

# Here is a list of interesting initial conditions:
initial_conditions = [(.1,1), (2,2),(1,1.3),(2,3),(2,1),(1,2)]

def cellular_switch(y,t,alpha,beta):
    """ Flow of Gardner's bistable cellular switch
    Args:
        y (array): (concentration of u, concentration of v)
        t (float): time (not used, autonomous system)
        alpha (float): maximum rate of repressor synthesis 
        beta (float): degree of cooperative behavior.
    Return: dy/dt
    """
    u, v = y
    return np.array([(alpha/(1+v**beta)) -u,
                     (alpha/(1+u**beta)) - v])
    
# Do the simulations
#Remember that we define f as the partial application of cellular_switch
trajectory = {}
for i,param in enumerate(scenarios):
    for j,ic in enumerate(initial_conditions):
        trajectory[i,j] = scipy.integrate.odeint(partial(cellular_switch, **param),
                                                 y0=ic,
                                                 t=time)


#Draw the traj
fig, ax = plt.subplots(2,2,figsize=(20,10))
for i,param in enumerate(scenarios):
    for j,ic in enumerate(initial_conditions):
        ax[i][0].set(xlabel='Time', ylabel='Concentration of u', title='Trajectory of u, {}'.format(param))
        ax[i][1].set(xlabel='Time', ylabel='Concentration of v', title='Trajectory of v, {}'.format(param))
        l = ax[i][0].plot(time,trajectory[i,j][:,0], label=ic)
        ax[i][1].plot(time,trajectory[i,j][:,1], color=l[0].get_color())
    ax[i][0].legend(title='Initial conditions')


#___________________________________________________________________
#Phase Diagram
    #Isoclines

#Null isiclines are the manifolds on which one component of the flow is null
    #Find the equation of null-isoclines for u and v

uspace = np.linspace(0,2,100)
vspace = np.linspace(0,2,100)

def plot_isocline(ax,uspace,vspace,alpha,beta,color = 'k', style = '--', opacity = .5):
    """ plot the isclones of the symetric cellular switch system"""
    ax.plot(uspace,alpha/(1+uspace**beta),style,color = color,alpha = opacity)
    ax.plot(alpha/(1+vspace**beta),vspace,style,color = color,alpha = opacity)
    ax.set(xlabel = 'u',ylabel = 'v')



fig1, ax1 = plt.subplots(3,3,figsize=(20,10))
    
    
    
    

