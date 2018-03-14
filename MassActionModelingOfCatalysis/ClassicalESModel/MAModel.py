import numpy as np

def func(y,t,k1,kO1,k2,inletS):
    E,S,ES,P = y
    dydt = [-k1*E*S + k2*ES + kO1*ES ,
            inletS - k1*E*S + kO1*ES ,
            k1*E*S - k2*ES - kO1*ES,
            k2*ES-P]
    return dydt
def params2(Eo, So, v, sigma):
    kneg1 = sigma
    Ces = v[1]/kneg1
    k1 = v[0]/(Eo-Ces)/So
    k2 = v[2] / Ces
    return (k1, kneg1, k2)
def func2(y,t,k1,kO1,k2,inletS):
    E,S,ES,P = y
    #k1 = k1+.05*k1*np.sin(3*t)+np.cos(5*t)
    dydt = [-(k1+.05*k1*(np.sin(3*t)+np.cos(5*t)))*E*S + k2*ES + kO1*ES ,
            inletS - (k1+.05*k1*(np.sin(3*t)+np.cos(5*t)))*E*S + kO1*ES ,
            (k1+.05* k1*(np.sin(3*t)+np.cos(5*t)))*E*S - k2*ES - kO1*ES,
            k2*ES-P]
    return dydt

def func3(y,t,k1,kO1,k2,inletS):
    E,S,ES,P = y
    dydt = [-k1*E*S + k2*ES + (kO1+.5*kO1*np.sin(t))*ES ,
            inletS - k1*E*S + (kO1+.5*kO1*np.sin(t))*ES ,
            k1*E*S - k2*ES - (kO1+.5*kO1*np.sin(t))*ES,
            k2*ES-P]
    return dydt