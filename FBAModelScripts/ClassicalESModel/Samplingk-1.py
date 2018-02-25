import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random
import cobra.io
from MAModel import *

###QUESTIONS TO BE ASKED#####

# In sampling postive bounded real space, do single fixed points emerge independent of k1
# Is there differential senstivity? If so, where?
# What about k-1 (or k2) what happens there?

# Other possibilites k-1 + k2 = 1 or sigma (probability)



myModel = cobra.io.read_sbml_model('ClassicalESModel.xml')
myModel.solver = 'gurobi'

sol = myModel.optimize()

fvaSol = cobra.flux_analysis.flux_variability_analysis(myModel, loopless=True)
print fvaSol
# print sol.fluxes
for x in myModel.reactions:
    print str(x.id) + ' : ' + str(x.reaction) + ' : ' + str(sol[x.id])

def params(Eo, So, v, sigma):
    kneg1 = sigma
    Ces = v[1]/kneg1
    k1 = v[0]/(Eo-Ces)/So
    k2 = v[2] / Ces
    return (k1, kneg1, k2)


sol2 = [sol[x.id] for x in myModel.reactions]
fig1, ax1 = plt.subplots(2, 2)
fig2, ax2 = plt.subplots(2, 2)

t = np.linspace(0, 1000, 1000)
file = open('RandK2Testing.txt', 'w')
meanError = []

kneg1_max = .4*10e3
param_min = .4
nsamples = 1000

p = [random.uniform(param_min, kneg1_max) for _ in range(nsamples)]
# p = [10**x+.001 for x in range(10)]

for z in p:

    sol3 = [x * 1. for x in sol2]

    y0 = [5.0, 1000., 0.0, 0.0]

    file.write(str(z) + '\n')
    k1, k01, k2 = params(y0[0], y0[1], sol3, z)
    # print type(sol)
    file.write(str(k1) + ' ')
    file.write(str(k01) + ' ')
    file.write(str(k2) + '\n')

    inletS = sol3[3]
    sol = odeint(func, y0, t, args=(k1, k01, k2, sol3[3]))

    finalConc = sol[-1, :]
    finalFlux = [k1 * finalConc[0] * finalConc[1]]
    finalFlux.append(k01 * finalConc[2])
    finalFlux.append(k2 * finalConc[2])
    print('k-1 = ' + str(z))
    print finalFlux
    print '\n\n'

    [file.write(str(x) + ' ') for x in finalFlux]
    file.write('\n' + str(np.mean(np.subtract(sol3[:-2], finalFlux) ** 2)) + '\n')
    meanError.append(np.sqrt(np.mean(np.subtract(sol3[:-2], finalFlux) ** 2)))
    if meanError[-1] < .01:
        ax1[0][0].plot(t, sol[:, 0])
        ax1[0][0].set_xlabel('t')
        ax1[0][0].set_ylabel('[E]')
        #ax1[0][0].set_title(
        #    '             Network Dynamics Sampling\n k-1 = [' + str(param_min) + ',' + str(kneg1_max) + '] n = ' + str(nsamples))
        ax1[0][0].set_ylim(-1,7)

        ax1[0][1].yaxis.set_label_position("right")
        ax1[0][1].yaxis.tick_right()
        ax1[0][1].plot(t, sol[:, 1])
        ax1[0][1].set_xlabel('t')
        ax1[0][1].set_ylabel('[S]')

        ax1[1][0].plot(t, sol[:, 2])
        ax1[1][0].set_xlabel('t')
        ax1[1][0].set_ylabel('[ES]')
        ax1[1][0].set_ylim(-1,7)


        ax1[1][1].yaxis.set_label_position("right")
        ax1[1][1].yaxis.tick_right()
        ax1[1][1].plot(t, sol[:, 3])
        ax1[1][1].set_xlabel('t')
        ax1[1][1].set_ylabel('[P]')

        ax2[0][0].scatter(z, sol[-1, 0])
        ax2[0][0].set_xlabel('k-1')
        ax2[0][0].set_ylabel('[E]')
      #  ax2[0][0].set_title('           Phase Plane: Final Concentration vs k-1')

        ax2[0][1].yaxis.set_label_position("right")
        ax2[0][1].yaxis.tick_right()
        ax2[0][1].scatter(z, sol[-1, 1])
        ax2[0][1].set_xlabel('k-1')
        ax2[0][1].set_ylabel('[S]')

        ax2[1][0].scatter(z, sol[-1, 2])
        ax2[1][0].set_xlabel('k-1')
        ax2[1][0].set_ylabel('[ES]')

        ax2[1][1].yaxis.set_label_position("right")
        ax2[1][1].yaxis.tick_right()
        ax2[1][1].scatter(z, sol[-1, 3])
        ax2[1][1].set_xlabel('k-1')
        ax2[1][1].set_ylabel('[P]')
        # """

fig1.tight_layout()
fig2.tight_layout()
pp = PdfPages('ClassicalESModelRandomK-1_'  + str(param_min) +'-' + str(kneg1_max) + '.pdf')
pp.savefig(fig1)
pp.savefig(fig2)
fig = plt.figure(3)
plt.scatter(p, meanError)
plt.title('Root Mean Squared Error')
plt.xlabel('k-1')
plt.ylabel('Error')
fig.tight_layout()
pp.savefig(fig)
plt.show()

pp.close()