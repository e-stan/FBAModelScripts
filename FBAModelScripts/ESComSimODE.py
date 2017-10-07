import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

import cobra.io

myModel = cobra.io.read_sbml_model('toyModel.xml')
myModel.solver = 'gurobi'

sol = myModel.optimize()

fvaSol = cobra.flux_analysis.flux_variability_analysis(myModel,loopless=True)
print fvaSol
#print sol.fluxes
for x in myModel.reactions:
    print str(x.id) + ' : ' + str(x.reaction) + ' : '  + str(sol[x.id])


def func(y,t,k1,kO1,k2,inletS):
    E,S,ES,P = y
    dydt = [-k1*E*S + k2*ES + kO1*ES ,
            inletS - k1*E*S + kO1*ES ,
            k1*E*S - k2*ES - kO1*ES,
            k2*ES]
    return dydt

def params(Eo,So,v):
    Ces = Eo * So * v[1] / ( v[4] + v[1] * (So  +1 ))
    k1 = v[0] / So / (Eo - Ces)
    kneg1 =v[1] / Ces
    k2 = v[2] / Ces
    return (k1,kneg1,k2)


y0 = [5.0, 1000.0, 0.0 , 0.0]

t = np.linspace(0, 5000, 5000)
sol2 = [sol[x.id] for x in myModel.reactions]
#print sol
k1,k01,k2 = params(y0[0],y0[1],sol)
#print type(sol)
print k1
print k01
print k2

inletS = 5
sol = odeint(func, y0, t, args=(k1,k01,k2,inletS))

fig = plt.figure(1)
ax1 = fig.add_subplot(221)
plt.plot(t, sol[:, 0])
plt.xlabel('t')
plt.ylabel('E')

ax2 = fig.add_subplot(222)
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()
plt.plot(t, sol[:, 1])
plt.xlabel('t')
plt.ylabel('S')

fig = plt.figure(1)
ax1 = fig.add_subplot(223)
plt.plot(t, sol[:, 2])
plt.xlabel('t')
plt.ylabel('ES')

ax2 = fig.add_subplot(224)
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()
plt.plot(t, sol[:, 3])
plt.xlabel('t')
plt.ylabel('P')

finalConc = sol[-1,:]
finalFlux = [k1*finalConc[0]*finalConc[1]]
finalFlux.append(k01*finalConc[2])
finalFlux.append(k2*finalConc[2])

print finalFlux
print np.mean((np.subtract(sol2[:-2],finalFlux))**2)
plt.show()