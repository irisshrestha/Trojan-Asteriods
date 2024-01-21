import numpy as np
import scipy.integrate
from matplotlib import animation
from matplotlib import pyplot as plt



"""
y=(x1,y1,vx1,vy1,x2,y2,vx2,vy2)
red is particle 1
blue is particle 2
"""
def gravity(t,y):
    m1=0.003
    m2=0.3

    dx1dt=y[2]
    dy1dt=y[3]

    dx2dt=y[6]
    dy2dt=y[7]

    dvx1dt= -m2*(y[0]-y[4])/((y[0]-y[4])**2+(y[1]-y[5])**2)**(3/2)
    dvy1dt= -m2*(y[1]-y[5])/((y[0]-y[4])**2+(y[1]-y[5])**2)**(3/2)

    dvx2dt= -m1*(y[4]-y[0])/((y[0]-y[4])**2+(y[1]-y[5])**2)**(3/2)
    dvy2dt= -m1*(y[5]-y[1])/((y[0]-y[4])**2+(y[1]-y[5])**2)**(3/2)

    return dx1dt,dy1dt,dvx1dt,dvy1dt,dx2dt,dy2dt,dvx2dt,dvy2dt

x1, y1, vx1, vy1 = -0.05, 0,0.5,0.5
x2, y2, vx2, vy2 = 0.05 , 0, -0.5 ,-0.5
solution = scipy.integrate.solve_ivp(gravity, (0,20),(x1,y1,vx1,vy1,x2,y2,vx2,vy2), t_eval=np.linspace(0,20,500))

fig = plt.figure(figsize=(10,6))


ax = plt.axes(xlim=(-0.5, 0.5), ylim=(-0.5, 0.5))
point, = ax.plot([solution.y[0][0],solution.y[1][0]],[solution.y[4][0],solution.y[5][0]],'bo')
#p1, = ax.plot(solution.y[0][0],solution.y[1][0], 'ro')
#p2, =  ax.plot(solution.y[4][0],solution.y[5][0], 'bo')
def animate(i):

    #p1.set_ydata(solution.y[1][i])
    #p1.set_xdata(solution.y[0][i])
    #p2.set_ydata(solution.y[5][i])
    #p2.set_xdata(solution.y[4][i])
    x=[solution.y[0][i],solution.y[4][i]]
    y=[solution.y[1][i],solution.y[5][i]]
    point.set_data(x,y)
    return point,

anim = animation.FuncAnimation(fig, animate,
                               frames=200, interval=50)



plt.show()