import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import pyplot as plt
import matplotlib; matplotlib.use('Qt5Agg')

def checkcollision(x1,y1,x2,y2):
    #checks of particles share same x and y co-ords within a limit
    if abs(x1-x2)<0.012 and abs(y1-y2)<0.012:
        return True
    else:
        return False

def velocityaftercollision(vx1,vy1,vx2,vy2):
#calculates new velocities after collision. For same mass elastic collision, they just swap velocities? I think
    xvelocityOfZMFrame= (vx1+vx2)/2
    yvelocityOfZMFrame = (vy1 + vy2)/2

    ZMFvx1= vx1-xvelocityOfZMFrame
    ZFMvx2= vx2 - xvelocityOfZMFrame
    ZFMvy1 = vy1 - yvelocityOfZMFrame
    ZFMvy2 = vy2 - yvelocityOfZMFrame

    newvx1= -ZMFvx1 + xvelocityOfZMFrame
    newvx2 = -ZFMvx2 + xvelocityOfZMFrame
    newvy1 = -ZFMvy1 + yvelocityOfZMFrame
    newvy2 = -ZFMvy2 + yvelocityOfZMFrame


    return [newvx1,newvy1,newvx2,newvy2]


#checks particles are within bounds
def xboundaries(x,vx):
    if x<=0 or x>=1:
        vx=-vx
    else:
        vx=vx

    return vx

def yboundaries(y,vy):
    if y <= 0 or y >= 1:
        vy = -vy
    else:
        vy=vy

    return vy

#set number of particles. Simulation is scalable
numberOfParticles = 20

#generate random initial conditions. Position can be from 0 to 1 in x and y. Velocity from -0.1 to 0.1
position = np.random.uniform(low=0, high=1, size=(numberOfParticles,2))
velocity= np.random.uniform(low=-0.1, high=0.1, size=(numberOfParticles,2))
rint= np.concatenate((position,velocity),axis=1)

#create array to store positions/velocities of each particle at each time dt
rcords= []
for i in range(numberOfParticles):
    value=[[],[],[],[]]
    rcords.append(value)


t, tmax, dt = 0, 1000, 0.03
tcords=[]

#array which is developing with time
r= np.copy(rint)

#simulate what happends after each time dt
while t<tmax:
    t += dt
    tcords.append(t)

    for i in range((len(r))):
        x=r[i][0]
        y=r[i][1]
        vx=r[i][2]
        vy=r[i][3]
        #travelling linearly
        x += vx * dt
        y += vy * dt
        rcords[i][0].append(x)
        rcords[i][1].append(y)

        #check boundaries
        vx=xboundaries(x,vx)
        vy=yboundaries(y,vy)
        rcords[i][2].append(vx)
        rcords[i][3].append(vy)
        r[i][0], r[i][1], r[i][2], r[i][3] = x, y, vx, vy

#check for collisions
    for i in range(len(r)):
        for p in range(len(r)):
            if p<i:

                if checkcollision(r[i][0],r[i][1],r[p][0],r[p][1]):


                    newvelocity=velocityaftercollision(r[i][2],r[i][3],r[p][2],r[p][3])
                    r[i][2] = newvelocity[0]
                    r[i][3] = newvelocity[1]
                    r[p][2] = newvelocity[2]
                    r[p][3] = newvelocity[3]






#animating and plotting
fig = plt.figure()
ax = plt.axes(xlim=(0, 1), ylim=(0, 1))
xplot=[]
yplot=[]
for i in range(len(rint)):

    xplot.append(rint[i][0])
for i in range(len(rint)):
    yplot.append(rint[i][1])

points, = ax.plot(xplot,yplot,'ro')


# animation function.  This is called sequentially
def animate(i):
    xplot = []
    yplot = []
    for l in range(len(rcords)):

        xplot.append(rcords[l][0][i])
    for l in range(len(rcords)):
        yplot.append(rcords[l][1][i])

    points.set_data(xplot, yplot)

    return points,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate,
                               frames=5000, interval=3)

plt.show()