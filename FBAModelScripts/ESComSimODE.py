import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def func(y,t,k1,kO1,k2,dS,dP,inletE,inletS):
    E,S,ES,P = y
    dydt = [inletE - k1*E*S + k2*ES + kO1*ES ,
            inletS - k1*E*S - dS*S + kO1*ES ,
            k1*E*S - k2*ES - kO1*ES,
            k2*ES - dP *P]
    return dydt


y0 = [0.0, 0.0, 0.0]

t = np.linspace(0, 1000, 100)

k1 =
kO1 =
k2 =
dS =
dP =
inletE =
inletS =
sol = odeint(func, y0, t, args=(k1,kO1,k2,dS,dP,inletE,inletS))

fig = plt.figure(1)
ax1 = fig.add_subplot(121)
plt.plot(t, sol[:, 0])
plt.xlabel('t')
plt.ylabel('E')

ax2 = fig.add_subplot(122)
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()
plt.plot(t, sol[:, 1])
plt.xlabel('t')
plt.ylabel('S')

fig = plt.figure(1)
ax1 = fig.add_subplot(121)
plt.plot(t, sol[:, 3])
plt.xlabel('t')
plt.ylabel('ES')

ax2 = fig.add_subplot(122)
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()
plt.plot(t, sol[:, 4])
plt.xlabel('t')
plt.ylabel('P')

plt.show()