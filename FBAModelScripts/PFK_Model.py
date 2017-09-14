# -*- coding: utf-8 -*-
#
#  open systems model of phophofructokinase as a dual substrate, unordered binding
#  catalytic system
#

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def func (y, t, k1, k2, k3, k4, k5, k6, k7, k8, k9, Fatp, Ff6p, outf6p, outatp, outf16bp):
    pfk, f6p, atp, f6ppfk, atppfk, f6ppfkatp, f16bp = y
    dydt = [ -k1*pfk*f6p + k2*f6ppfk - k3*pfk*atp + k4*atppfk + k9*f6ppfkatp,
             -k1*pfk*f6p + k2*f6ppfk - k8*f6p*atppfk + k7*f6ppfkatp + Ff6p - outf6p*f6p,
             -k3*pfk*atp + k4*atppfk - k5*atp*atppfk + k6*f6ppfkatp + Fatp - outatp*atp, 
              k1*pfk*f6p - k2*f6ppfk - k5*atp*f6ppfk + k6*f6ppfkatp,
              k3*pfk*atp - k4*atppfk - k8*f6p*atppfk + k7*f6ppfkatp,
              k5*atp*f6ppfk + k8*f6p*atppfk - (k5 + k6 + k9)*f6ppfkatp,
              k9*f6ppfkatp - outf16bp*f16bp]
    return dydt

def pfk_model(max_time=10):
    """
    Simulate Phosphofructokinase-mediated conversion of fructose-6-phosphate into
       1,6 fructose bisphophate, the third step in glycolysis, and plot the integration results
       
    Parameters
    ----------
    max_time: number
        maximum time for integration interval
        
        
    Returns
    -------
    7 graphs describing the results of integrating the model ode's
    
    
    Examples
    --------
    >>>pfk_model()
    integrates the model ode's and plots results for the default time interval
    (0,10)
    
    >>>pfk_model(50)
    integrates the model ode's and plots results for the time interval
    (0,50)
    
    """
#
# rate constants
#
    k1 = 1.0
    k2 = 0.1
    k3 = 1.0
    k4 = 0.1
    k5 = 1.0
    k6 = 0.1
    k7 = 1.0
    k8 = 0.1
    k9 = 1.0
#
#  input flux
#
    Fatp=50.0
    Ff6p=50.0
#
# output diffusion/ degradation constants
#
    outf6p = 1.0;
    outatp = 1.0;
    outf16bp = 1.0;
#
# 
#

    y0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    
    t = np.linspace(0, max_time, 10*max_time)
    
    sol = odeint(func, y0, t, args=(k1, k2, k3, k4, k5, k6, k7, k8, k9, Fatp, Ff6p, outf6p, outatp, outf16bp))
    
    fig = plt.figure(1)
    ax1 = fig.add_subplot(121)
    plt.plot(t, sol[:, 0])
    plt.xlabel('t')
    plt.ylabel('pfk')
    
    ax2 = fig.add_subplot(122)
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    plt.plot(t, sol[:, 1])
    plt.xlabel('t')
    plt.ylabel('f6p')
    
    fig = plt.figure(2)
    ax1 = fig.add_subplot(121)
    plt.plot(t, sol[:, 2], 'b')
    plt.xlabel('t')
    plt.ylabel('atp')
    
    ax2 = fig.add_subplot(122)
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    plt.plot(t, sol[:, 3], 'b')
    plt.xlabel('t')
    plt.ylabel('f6p-pfk')
    
    fig = plt.figure(3)
    ax1 = fig.add_subplot(121)
    plt.plot(t, sol[:, 4], 'b')
    plt.xlabel('t')
    plt.ylabel('pfk-atp')
    
    ax2 = fig.add_subplot(122)
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    plt.plot(t, sol[:, 5], 'b')
    plt.xlabel('t')
    plt.ylabel('f6p-pfk-atp')
    
    plt.figure(4)
    plt.plot(t, sol[:, 6], 'b')
    plt.xlabel('t')
    plt.ylabel('f1,6bp')
    
    plt.show()


pfk_model()

