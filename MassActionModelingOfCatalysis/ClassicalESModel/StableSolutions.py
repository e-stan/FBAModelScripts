import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random
import cobra.io
from MAModel import *

def params(Eo,So,v,sigma):
    k1 = sigma
    Ces = Eo - v[0]/So/k1
    kneg1 =v[1] / Ces
    k2 = v[2] / Ces
    return (k1,kneg1,k2)



def params2(Eo, So, v, sigma):
    kneg1 = sigma
    Ces = v[1]/kneg1
    k1 = v[0]/(Eo-Ces)/So
    k2 = v[2] / Ces
    return (k1, kneg1, k2)
myModel = cobra.io.read_sbml_model('ClassicalESModel.xml')
myModel.solver = 'gurobi'

fig2,ax2 = plt.subplots(2,2)


sol = myModel.optimize()
t = np.linspace(0, 25, 1000)


fvaSol = cobra.flux_analysis.flux_variability_analysis(myModel,loopless=True)
print fvaSol
#print sol.fluxes
for x in myModel.reactions:
    print str(x.id) + ' : ' + str(x.reaction) + ' : '  + str(sol[x.id])

sol2 = [sol[x.id] for x in myModel.reactions]

y0 = [5.0, 1000.0, 0.0, 0.0]

k1,k01,k2 = params2(y0[0],y0[1],sol2,141.4)

inletS = sol2[3]
sol = odeint(func, y0, t, args=(k1, k01, k2, sol2[3]))

ax2[0][0].plot(t, sol[:, 0])
ax2[0][0].set_xlabel('t')
ax2[0][0].set_ylabel('[E]')
#ax1[0][0].set_yticks([x-.5 for x in range(7)])
# ax1[0][0].set_title('             Network Dynamics Sampling\n k1 = ['+ str(param_min) +  ',' + str(k1_max)+'] n = '+ str(nsamples))


ax2[0][1].yaxis.set_label_position("right")
ax2[0][1].yaxis.tick_right()
ax2[0][1].plot(t, sol[:, 1])
ax2[0][1].set_xlabel('t')
ax2[0][1].set_ylabel('[S]')
#ax1[0][1].set_yticks([x+996 for x in range(5)])

ax2[1][0].plot(t, sol[:, 2])
ax2[1][0].set_xlabel('t')
ax2[1][0].set_ylabel('[ES]')
#ax1[1][0].set_yticks([x-.5 for x in range(7)])

ax2[1][1].yaxis.set_label_position("right")
ax2[1][1].yaxis.tick_right()
ax2[1][1].plot(t, sol[:, 3])
ax2[1][1].set_xlabel('t')
ax2[1][1].set_ylabel('[P]')
#ax1[1][1].set_yticks('')

k1,k01,k2 = params(y0[0],y0[1],sol2,12.24)

sol = odeint(func, y0, t, args=(k1, k01, k2, sol2[3]))

ax2[0][0].plot(t, sol[:, 0])
ax2[0][0].set_xlabel('t')
ax2[0][0].set_ylabel('[E]')
# ax1[0][0].set_title('             Network Dynamics Sampling\n k1 = ['+ str(param_min) +  ',' + str(k1_max)+'] n = '+ str(nsamples))


ax2[0][1].yaxis.set_label_position("right")
ax2[0][1].yaxis.tick_right()
ax2[0][1].plot(t, sol[:, 1])
ax2[0][1].set_xlabel('t')
ax2[0][1].set_ylabel('[S]')

ax2[1][0].plot(t, sol[:, 2])
ax2[1][0].set_xlabel('t')
ax2[1][0].set_ylabel('[ES]')

ax2[1][1].yaxis.set_label_position("right")
ax2[1][1].yaxis.tick_right()
ax2[1][1].plot(t, sol[:, 3])
ax2[1][1].set_xlabel('t')
ax2[1][1].set_ylabel('[P]')

fig2.tight_layout()

plt.show()