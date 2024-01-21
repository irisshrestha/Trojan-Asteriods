import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
from matplotlib import animation

"""
All in units of solar system units
"""
#Create class orbits which will contain Jupiter and the asteroids
class Orbits:
     def __init__(self,x,y,vx,vy):
         self.x =x
         self.y =y
         self.vx = vx
         self.vy = vy

     def getx(self):
        return self.x

#Function for gravity to solve 2 body ODE
def gravity(y,t,m): #y is an array consisting of [x1 , y1 , vx1, vy1,x2, y2 , vx2, vy2]
    global G
    x1,y1,x2,y2=y[0],y[1],y[4],y[5]
    vx1,vy1,vx2,vy2 = y[2],y[3],y[6],y[7]
    r=(x2-x1)**2+(y2-y1)**2

    m1=m[0]
    m2=m[1]

    dx1dt, dy1dt = vx1, vy1
    dx2dt, dy2dt = vx2, vy2

    dvx1dt = G*m2*(x2-x1)/r**(3/2)
    dvy1dt = G*m2*(y2-y1)/r**(3/2)

    dvx2dt = G*m1*(x1-x2)/r**(3/2)
    dvy2dt = G*m1*(y1-y2)/r**(3/2)

    return dx1dt, dy1dt, dvx1dt, dvy1dt, \
           dx2dt, dy2dt, dvx2dt, dvy2dt

#Set up constants
G=4*np.pi**2
massSun = 1
massJup = 0.1
jupRadius=5.2


#Create orbit for Jupiter and Sun
#Set initial conditions
mass = [massJup,massSun]
#Initial position with (0,0) for sun

CoM= jupRadius*massJup/(massSun+massJup)
jupx=jupRadius*massSun/(massSun+massJup)
#Initial velocities for circular orbit (Constant distance apart)
jupInitial_v=np.sqrt(G*massSun**2/(jupRadius*(massSun+massJup)))
sunInitial_v=-jupInitial_v*massJup/massSun
#Period for orbits. Should be the same
jupPeriod= 2*np.pi*jupx/jupInitial_v


t = np.linspace(0,100*jupPeriod,100*100)
#Solve ODE with sun's gravity
solution = scipy.integrate.odeint(gravity,[jupx,0,0,jupInitial_v,-CoM,0,0,sunInitial_v],t,args=(mass,))
#Create object jupiter
jupiter = Orbits(solution[:,0],solution[:,1],solution[:,2],solution[:,3])
sun = Orbits(solution[:,4],solution[:,5],solution[:,6],solution[:,7])

fig = plt.figure(figsize=(10,6))
ax = plt.axes(xlim=(-6, 6), ylim=(-6, 6))
plt.axis('square')
plt.title("Two body orbits around their Centre of Mass")
plt.plot(jupiter.x,jupiter.y)
plt.plot(sun.x,sun.y)




pointjup, = ax.plot(jupiter.x[0], jupiter.y[0], 'o')
pointsun, = ax.plot(sun.x[0], sun.y[0], 'o', mfc="yellow", ms="15")
"""
def animatejup(i):
    pointjup.set_ydata(jupiter.y[i])
    pointjup.set_xdata(jupiter.x[i])
    pointsun.set_ydata(sun.y[i])
    pointsun.set_xdata(sun.x[i])
    return pointsun, pointjup,

anim = animation.FuncAnimation(fig, animatejup,
                               frames=2000, interval=0.1)
"""
plt.show()