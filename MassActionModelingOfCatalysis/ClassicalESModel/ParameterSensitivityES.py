import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random
import cobra.io
from MAModel import *


def dCdk1(k1,v1,s0):
    return v1/s0/k1**2
def dCdk_1(k_1,v2):
    return v2/k_1**2
def dCdk2(k2,v3):
    return v3/k2**2

myModel = cobra.io.read_sbml_model('ClassicalESModel.xml')
myModel.solver = 'gurobi'

sol = myModel.optimize()

fvaSol = cobra.flux_analysis.flux_variability_analysis(myModel, loopless=True)
print fvaSol
# print sol.fluxes
for x in myModel.reactions:
    print str(x.id) + ' : ' + str(x.reaction) + ' : ' + str(sol[x.id])


sol2 = [sol[x.id] for x in myModel.reactions]


k1_min = .003
k_1_min = .4
k2_min = 2.6

paramax = 10

k1 = np.linspace(k1_min,paramax,1000)
k_1 = np.linspace(k_1_min,paramax,1000)
k2 = np.linspace(k2_min,paramax,1000)


y0 = [5.0, 1000., 0.0, 0.0]


plt.figure()
plt.plot(k1,[dCdk1(x,sol2[0],y0[1]) for x in k1])
plt.xticks(range(10))
plt.xlabel("k1")
plt.ylabel("Magnitude of the Rate of Change of [ES]")
plt.savefig("k1Sens.png")


plt.figure()
plt.plot(k_1,[dCdk_1(x,sol2[1]) for x in k_1])
plt.xticks(range(10))
plt.xlabel("k-1")
#plt.ylabel("Magnitude of the Rate of Change of [ES]")
plt.savefig("k_1Sens.png")


plt.figure()
plt.plot(k2,[dCdk2(x,sol2[3]) for x in k2])
plt.xticks(range(10))
plt.xlabel("k2")
#plt.ylabel("Magnitude of the Rate of Change of [ES]")
plt.savefig("k2Sens.png")


plt.show()