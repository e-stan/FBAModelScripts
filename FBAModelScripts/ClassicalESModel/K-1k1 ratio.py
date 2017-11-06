import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random
import cobra.io

###QUESTIONS TO BE ASKED#####

# In sampling postive bounded real space, do single fixed points emerge independent of k1
# Is there differential senstivity? If so, where?
# What about k-1 (or k2) what happens there?

# Other possibilites k-1 + k2 = 1 or sigma (probability)

k1_max = 2000
para_min = 100
nsamples = 1000

myModel = cobra.io.read_sbml_model('ClassicalESModel.xml')
myModel.solver = 'gurobi'

sol = myModel.optimize()

fvaSol = cobra.flux_analysis.flux_variability_analysis(myModel, loopless=True)
print fvaSol
# print sol.fluxes
for x in myModel.reactions:
    print str(x.id) + ' : ' + str(x.reaction) + ' : ' + str(sol[x.id])


def func(y, t, k1, kO1, k2, inletS):
    E, S, ES, P = y
    dydt = [-k1 * E * S + k2 * ES + kO1 * ES,
            inletS - k1 * E * S + kO1 * ES,
            k1 * E * S - k2 * ES - kO1 * ES,
            k2 * ES]
    return dydt


def params(Eo, So, v, sigma):
    k1 = sigma
    Ces = Eo - v[0] / So / k1
    kneg1 = v[1] / Ces
    k2 = v[2] / Ces
    return (k1, kneg1, k2)


sol2 = [sol[x.id] for x in myModel.reactions]
fig1, ax1 = plt.subplots(2, 1)
#fig2, ax2 = plt.subplots(2, 2)

t = np.linspace(0, 1000, 1000)
file = open('RandK1Testing.txt', 'w')
meanError = []

p = [random.uniform(para_min, k1_max) for _ in range(nsamples)]
# p = [10**x+.001 for x in range(10)]
y0 = [5.0, 1000., 0.0, 0.0]
k1,k01,k2 =params(y0[0], y0[1], sol2, para_min)
kone = [k01]
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
    print('k1 = ' + str(z))
    print finalFlux
    print '\n\n'

    [file.write(str(x) + ' ') for x in finalFlux]
    file.write('\n' + str(np.mean(np.subtract(sol3[:-2], finalFlux) ** 2)) + '\n')
    meanError.append(np.sqrt(np.mean(np.subtract(sol3[:-2], finalFlux) ** 2)))
    if meanError[-1] < .01:
        ax1[0].scatter(k1, k01-kone[-1])
        ax1[0].set_xlabel('k1')
        ax1[0].set_ylabel('k1/k-1')
        ax1[0].set_title(
            '             Rate Constant Relationships\n k1,k-1 = [' + str(para_min) +',' + str(k1_max) + '] n = ' + str(nsamples))
        kone.append(k01)


def params(Eo, So, v, sigma):
    kneg1 = sigma
    Ces = v[1]/kneg1
    k1 = v[0]/(Eo-Ces)/So
    k2 = v[2] / Ces
    return (k1, kneg1, k2)
k1,k01,k2 =params(y0[0], y0[1], sol2, para_min)
kneg1 = [k1]
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

        ax1[1].scatter(k01,k1-kneg1[-1])
        ax1[1].set_xlabel('k-1')
        ax1[1].set_ylabel('k-1/k1')
        kneg1.append(k1)


print "std k1 = " + str(np.std(kneg1))
print "std k-1 = " + str(np.std(kone))

fig1.tight_layout()
#pp = PdfPages('ClassicalESModelk-1/K10-' + str(k1_max) + '.pdf')
#pp.savefig(fig1)

plt.show()

#pp.close()