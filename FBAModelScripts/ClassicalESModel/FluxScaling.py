import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import cobra.io

myModel = cobra.io.read_sbml_model('ClassicalESModel.xml')
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

sol2 = [sol[x.id] for x in myModel.reactions]
fig = plt.figure(1)

t = np.linspace(0, 1000, 1000)
file = open('ICTesting.txt','w')
meanError = []
#p = [1/10.,1,2,3,4,5,6,7,8,9,10,100,1000]
p = [5,50,500,1000]
p = [.000001,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0]
#p = [5+2*x for x in range(500)]
for z in p:

    sol3 = [x * z * 1000000 for x in sol2]

    y0 = [5.0, 1000, 0 , 0.0]


    file.write(str(z)+ '\n')
    k1,k01,k2 = params(y0[0],y0[1],sol3)
    #print type(sol)
    file.write(str(k1)+' ')
    file.write(str(k01)+ ' ')
    file.write(str(k2)+ '\n')

    inletS = sol3[3]
    sol = odeint(func, y0, t, args=(k1,k01,k2,sol3[3]))

    ax1 = fig.add_subplot(221)
    plt.plot(t, sol[:, 0])
    plt.xlabel('t')
    plt.ylabel('E')
    plt.title('           ES Dyanmics V =  V0*k')


    ax2 = fig.add_subplot(222)
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    plt.plot(t, sol[:, 1])
    plt.xlabel('t')
    plt.ylabel('S')

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
    print('So = '+str(z))
    print finalFlux
    print '\n\n'

    [file.write(str(x) + ' ') for x in finalFlux]
    file.write('\n'+str(np.mean(np.subtract(sol3[:-2],finalFlux)**2))+ '\n')
    meanError.append(np.sqrt(np.mean(np.subtract(sol3[:-2],finalFlux)**2)))

#fig.tight_layout()
pp = PdfPages('ClassicalESModelVScaling1-1000000.pdf')
pp.savefig(fig)
fig = plt.figure(2)
plt.plot([x*1000000 for x in p],meanError)
plt.title('Root Mean Squared Error')
plt.xlabel('k')
plt.ylabel('Error')
fig.tight_layout()
pp.savefig(fig)
plt.show()

pp.close()