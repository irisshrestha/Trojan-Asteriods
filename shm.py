import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

def shm(t,y):
    dydt = y[1]
    dvdt= -y[0]
    return dydt,dvdt

sol = scipy.integrate.solve_ivp(shm,(0,10),(0,1),t_eval= np.linspace(0,10,100))

plt.plot(sol.t,sol.y[0])
plt.show()

