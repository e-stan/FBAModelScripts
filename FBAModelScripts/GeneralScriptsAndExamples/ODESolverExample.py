import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def func(y,t,k1,k2,inletA,inletB):
    A,B,AB = y
    dydt = [inletA - k1*A*B,inletB-k1*A*B,k1*A*B - k2*AB]
    return dydt


y0 = [0.0, 0.0, 0.0]

t = np.linspace(0, 1000, 100)

sol = odeint(func, y0, t, args=(.2,1,5,5))

fig = plt.figure(1)
ax1 = fig.add_subplot(121)
plt.plot(t, sol[:, 0])
plt.xlabel('t')
plt.ylabel('A')

ax2 = fig.add_subplot(122)
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()
plt.plot(t, sol[:, 1])
plt.xlabel('t')
plt.ylabel('B')

fig = plt.figure(2)
ax1 = fig.add_subplot(121)
plt.plot(t, sol[:, 2], 'b')
plt.xlabel('t')
plt.ylabel('AB')

plt.show()