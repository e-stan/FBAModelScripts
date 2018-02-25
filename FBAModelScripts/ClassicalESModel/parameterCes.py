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

k1_max = 5-1e-5
para_min = 1e-5
nsamples = 1000

myModel = cobra.io.read_sbml_model('ClassicalESModel.xml')
myModel.solver = 'gurobi'

sol = myModel.optimize()

fvaSol = cobra.flux_analysis.flux_variability_analysis(myModel, loopless=True)
print fvaSol
# print sol.fluxes
for x in myModel.reactions:
    print str(x.id) + ' : ' + str(x.reaction) + ' : ' + str(sol[x.id])




def params(Eo, So, v, sigma):
    Ces = sigma
    kneg1 = v[1] / Ces
    k2 = v[2] / Ces
    k1 = v[0]/(Eo - Ces)/So
    return (k1, kneg1, k2)


sol2 = [sol[x.id] for x in myModel.reactions]
fig1, ax1 = plt.subplots(1, 1)
#fig2, ax2 = plt.subplots(2, 2)

t = np.linspace(0, 1000, 1000)
file = open('RandK1Testing.txt', 'w')
meanError = []

p = [random.uniform(para_min, k1_max) for _ in range(nsamples)]
# p = [10**x+.001 for x in range(10)]
y0 = [5.0, 1000., 0.0, 0.0]
k1,k01,k2 =params(y0[0], y0[1], sol2, para_min)
kone = []
knegone = []
p.sort()
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
    print('Ces = ' + str(z))
    print finalFlux
    print '\n\n'

    [file.write(str(x) + ' ') for x in finalFlux]
    file.write('\n' + str(np.mean(np.subtract(sol3[:-2], finalFlux) ** 2)) + '\n')
    meanError.append(np.sqrt(np.mean(np.subtract(sol3[:-2], finalFlux) ** 2)))
    if meanError[-1] < .01:
        ax1.scatter(1/k1, 1/k01)
        ax1.set_xlabel('1/k1')
        ax1.set_ylabel('1/k-1')
        ax1.set_title(
            '             Rate Constant Relationships as a function of Ces[B]\n Ces[B] = [' + str(para_min) +',' + 'Eo] n = ' + str(nsamples))
        knegone.append((1/k01))
        kone.append(1/k1)

ax1.plot(kone,knegone)

z = np.polyfit(kone, knegone, 1)
#p = np.poly1d(z)
#ax1.plot(kone,p(x),"r--")
# the line equation:
print "y=%.6fx+(%.6f)"%(z[0],z[1])
fig1.tight_layout()
pp = PdfPages('ClassicalESModelk-1K10-' + str(k1_max) + '.pdf')
pp.savefig(fig1)



def params_2(Eo, So, v, a , b):
    k1 = (v[0] + So*v[1]*a)/(So*Eo - So* v[1] * b)
    Ces = Eo - v[0] / So / k1
    kneg1 = v[1] / Ces
    k2 = v[2] / Ces
    return (k1, kneg1, k2)


fig3, ax3 = plt.subplots(2, 2)

t = np.linspace(0, 1000, 1000)
file = open('RandK1Testing.txt', 'w')
meanError = []



sol3 = [x * 1. for x in sol2]

y0 = [5.0, 1000., 0.0, 0.0]

file.write(str(z) + '\n')
print sol3
k1, k01, k2 = params_2(y0[0], y0[1], sol3, z[0],z[1])
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
print('k1 = ' + str(z))
print finalFlux
print '\n\n'

ax3[0][0].plot(t, sol[:, 0])
ax3[0][0].set_xlabel('t')
ax3[0][0].set_ylabel('E')
ax3[0][0].set_title(
    '             Network Dynamics Sampling\n k1 = ' + str(k1))

ax3[0][1].yaxis.set_label_position("right")
ax3[0][1].yaxis.tick_right()
ax3[0][1].plot(t, sol[:, 1])
ax3[0][1].set_xlabel('t')
ax3[0][1].set_ylabel('S')

ax3[1][0].plot(t, sol[:, 2])
ax3[1][0].set_xlabel('t')
ax3[1][0].set_ylabel('ES')

ax3[1][1].yaxis.set_label_position("right")
ax3[1][1].yaxis.tick_right()
ax3[1][1].plot(t, sol[:, 3])
ax3[1][1].set_xlabel('t')
ax3[1][1].set_ylabel('P')

print k1
print k01
print k2

plt.show()

pp.close()