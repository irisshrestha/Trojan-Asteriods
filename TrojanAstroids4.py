import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
from matplotlib import animation

"""
All in units of solar system units
"""
#Set up constants
G = 4*np.pi**2
massSun = 1
massJup = 0.001
massTotal= massSun + massJup
radius = 5.2

#Number of oscillations to plot
numOss = 100

#Perturbations
anglePert = 0    #Angular perturbations in Rads
radiusPert = 0    #Radial perturbation in AU



"""
Create class to store information about all bodies
(x position, y position, x velocity, y velocity)
"""
class Bodies:
     def __init__(self,x,y,vx,vy):
         self.x = x
         self.y = y
         self.vx = vx
         self.vy = vy



"""
Create class to store distances between the different bodies
distance from (origin, sun, jupiter, greek, trojan) 
"""
class Distances:
    def __init__(self,origin, sun, jupiter, greek, trojan):
        self.origin = origin
        self.sun = sun
        self.jupiter = jupiter
        self.greek = greek
        self.trojan = trojan


"""
Function for gravity containing all ODEs to solve
y is an array consisting of all positional and velocity information for all bodies
y= [jupiter.x , jupiter.y , jupiter.vx, jupiter.vy
    sun.x, sun.y , sun.vx, sun.vy,
    greek.x, greek.y, greek.vx, greek.vy
    trojan.x, trojan.y, trojan.vx, trojan.vy ]
"""
def gravity(y,t):
    global G, massSun, massJup
    #Create objects to hold information about each body
    jupiter = Bodies( y[0], y[1], y[2], y[3])
    sun= Bodies( y[4], y[5], y[6], y[7])
    greek= Bodies(  y[8], y[9], y[10], y[11] )
    trojan= Bodies( y[12], y[13], y[14], y[15])

    #Calculate relevant distances
    rsj= (sun.x-jupiter.x)**2 + (sun.y-jupiter.y)**2    #Distance from sun to jupiter

    rsg= (sun.x-greek.x)**2 + (sun.y-greek.y)**2    #Distance from sun to Greek asteriods
    rjg= (jupiter.x-greek.x)**2 + (jupiter.y-greek.y)**2    #Distance from jupiter to Greek asteriods

    rst= (sun.x-trojan.x)**2 + (sun.y-trojan.y)**2      #Distance from sun to Trojan asteriods
    rjt= (jupiter.x-trojan.x)**2 + (jupiter.y-trojan.y)**2   #Distance from jupiter to Trojan asteriods

    #Sun and Jupiter two body orbit ODEs
    dx1dt, dy1dt = jupiter.vx, jupiter.vy
    dx2dt, dy2dt = sun.vx, sun.vy
    dvx1dt = G * massSun * (sun.x - jupiter.x) / rsj**(3/2)
    dvy1dt = G * massSun * (sun.y - jupiter.y) / rsj**(3/2)
    dvx2dt = G * massJup * (jupiter.x - sun.x) / rsj**(3/2)
    dvy2dt = G * massJup * (jupiter.y - sun.y) / rsj**(3/2)

    #Greek orbit ODEs
    dx3dt, dy3dt= greek.vx, greek.vy
    dvx3dt = G * ( (massSun * (sun.x - greek.x) / rsg**(3/2)) + (massJup * (jupiter.x - greek.x) / rjg**(3/2)) )
    dvy3dt = G * ( (massSun * (sun.y - greek.y) / rsg**(3/2)) + (massJup * (jupiter.y - greek.y) / rjg**(3/2)) )

    #Trojan Orbit ODEs
    dx4dt, dy4dt = trojan.vx, trojan.vy
    dvx4dt = G * ((massSun * (sun.x - trojan.x) / rst ** (3 / 2)) + (massJup * (jupiter.x - trojan.x) / rjt ** (3 / 2)))
    dvy4dt = G * ((massSun * (sun.y - trojan.y) / rst ** (3 / 2)) + (massJup * (jupiter.y - trojan.y) / rjt ** (3 / 2)))

    return dx1dt, dy1dt, dvx1dt, dvy1dt, \
           dx2dt, dy2dt, dvx2dt, dvy2dt, \
           dx3dt, dy3dt, dvx3dt, dvy3dt, \
           dx4dt, dy4dt, dvx4dt, dvy4dt



"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Calculate Initial Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

#Angluar displacement of asteriods from Sun
greekInitialAngle = np.pi / 3
trojanInitialAngle = -np.pi / 3

"""
Jupiter and Sun Initial Conditions
"""
sunInitialdis = radius * massJup / (massSun + massJup)     #Distance between Sun and CoM
jupInitialdis = radius * massSun / (massSun + massJup)      #Distance between Jupiter and CoM

jupInitialv = np.sqrt(G * massSun**2 / (radius * (massSun + massJup)))   #Jupiter initial velocity
sunInitialv = -jupInitialv * massJup / massSun                           #Sun Initial Velocity

#Period for orbits. Should be the same
angularVelocity = jupInitialv / jupInitialdis
jupPeriod = 2*np.pi*jupInitialdis / jupInitialv
sunPeriod = 2*np.pi*sunInitialdis / sunInitialv

"""
Asteroid Initial Conditions
"""

"""
Function to calculate radius and angle from CoM to get initial positions for asteroids
Input = (Distance of Asteroid from sun , Distance of sun from CoM, Distance of Jupiter from CoM, Angle of asteroids from Sun)
Output = (Distance of Asteroid from CoM , Angle of Asteroids from CoM)
"""
def fromCom(radiusFromSun, sunComDis, jupComDis, angle):
    comdis = np.sqrt(sunComDis**2 + radiusFromSun**2 - 2 * sunComDis * radiusFromSun * np.cos(angle))
    jupdis = np.sqrt(2 * radius**2 * (1 - np.cos(angle)) )
    comangle = np.sign(angle) * np.arccos((comdis**2 + jupComDis**2 - jupdis**2) / (2 * comdis * jupComDis))

    return comdis, comangle

"""
Greek asteriod initial conditions
"""
greekComInitial = fromCom(radius, sunInitialdis, jupInitialdis, greekInitialAngle)
#(Initial distance of Greek from CoM , Initial angle of Greek from CoM)
greekComInitial = (greekComInitial[0] +radiusPert, greekComInitial[1] + anglePert)  #Apply Perturbations. Only greek is perturbed
greekInitialVelocity = jupInitialv * (greekComInitial[0]) / jupInitialdis           #Initial Greek velocity

#Create object to store all initial conditions
greekInitial = Bodies(np.cos(greekComInitial[1] ) * (greekComInitial[0]),
                     np.sin(greekComInitial[1] ) * (greekComInitial[0]),
                     -np.sin(greekComInitial[1] ) * greekInitialVelocity,
                     np.cos(greekComInitial[1] ) * greekInitialVelocity)


"""
Trojan asteriod initial conditions
"""
trojanComInitial = fromCom(radius, sunInitialdis, jupInitialdis, trojanInitialAngle)
#(Initial distance of Trojan from CoM , Initial angle of Trojan from CoM)
trojanInitialVelocity = jupInitialv * trojanComInitial[0]/ jupInitialdis     #Initial Trojan velocity

#Create object to store all initial conditions
trojanInitial = Bodies(np.cos(trojanComInitial[1] ) * trojanComInitial[0] ,
                       np.sin(trojanComInitial[1] ) * trojanComInitial[0],
                       -np.sin(trojanComInitial[1] ) * trojanInitialVelocity,
                       np.cos(trojanComInitial[1] ) * trojanInitialVelocity )




"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Solving ODEs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
#Times to evaluate ODE based on number of ossicillations
t = np.linspace(0, numOss * jupPeriod, 1000*numOss)

#Solve ODEs using odeint
solution = scipy.integrate.odeint(gravity,[jupInitialdis, 0, 0, jupInitialv,
                                           -sunInitialdis, 0, 0, sunInitialv,
                                           greekInitial.x, greekInitial.y, greekInitial.vx, greekInitial.vy,
                                           trojanInitial.x, trojanInitial.y, trojanInitial.vx, trojanInitial.vy]
                                  ,t)

"""
Create objects to store information for each orbits
Each is an array that holds all x, y, vx, vy information for all times t
for each body
"""
jupiterOrbit = Bodies(solution[:,0],solution[:,1],solution[:,2],solution[:,3])
sunOrbit = Bodies(solution[:,4],solution[:,5],solution[:,6],solution[:,7])
greekOrbit = Bodies(solution[:,8],solution[:,9],solution[:,10],solution[:,11])
trojanOrbit = Bodies(solution[:,12],solution[:,13],solution[:,14],solution[:,15])



"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Measurements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

"""
Create objects to store distance between all bodies
distance from (origin, sun, jupiter, greek, trojan)
"""
disSun = Distances(np.sqrt(sunOrbit.x**2 + sunOrbit.y**2),
                       np.sqrt((sunOrbit.x - sunOrbit.x)**2 + (sunOrbit.y - sunOrbit.y)**2),
                       np.sqrt((sunOrbit.x - jupiterOrbit.x)**2 + (sunOrbit.y - jupiterOrbit.y)**2),
                       np.sqrt((sunOrbit.x - greekOrbit.x)**2 + (sunOrbit.y - greekOrbit.y)**2),
                       np.sqrt((sunOrbit.x - trojanOrbit.x)**2 + (sunOrbit.y - trojanOrbit.y)**2) )

disJupiter = Distances(np.sqrt(jupiterOrbit.x**2 + jupiterOrbit.y**2),
                       np.sqrt((jupiterOrbit.x - sunOrbit.x)**2 + (jupiterOrbit.y - sunOrbit.y)**2),
                       np.sqrt((jupiterOrbit.x - jupiterOrbit.x)**2 + (jupiterOrbit.y - jupiterOrbit.y)**2),
                       np.sqrt((jupiterOrbit.x - greekOrbit.x)**2 + (jupiterOrbit.y - greekOrbit.y)**2),
                       np.sqrt((jupiterOrbit.x - trojanOrbit.x)**2 + (jupiterOrbit.y - trojanOrbit.y)**2) )

disGreek = Distances(np.sqrt(greekOrbit.x**2 + greekOrbit.y**2),
                       np.sqrt((greekOrbit.x - sunOrbit.x)**2 + (greekOrbit.y - sunOrbit.y)**2),
                       np.sqrt((greekOrbit.x - jupiterOrbit.x)**2 + (greekOrbit.y - jupiterOrbit.y)**2),
                       np.sqrt((greekOrbit.x - greekOrbit.x)**2 + (greekOrbit.y - greekOrbit.y)**2),
                       np.sqrt((greekOrbit.x - trojanOrbit.x)**2 + (greekOrbit.y - trojanOrbit.y)**2) )

disTrojan = Distances(np.sqrt(trojanOrbit.x**2 + trojanOrbit.y**2),
                       np.sqrt((trojanOrbit.x - sunOrbit.x)**2 + (trojanOrbit.y - sunOrbit.y)**2),
                       np.sqrt((trojanOrbit.x - jupiterOrbit.x)**2 + (trojanOrbit.y - jupiterOrbit.y)**2),
                       np.sqrt((trojanOrbit.x - greekOrbit.x)**2 + (trojanOrbit.y - greekOrbit.y)**2),
                       np.sqrt((trojanOrbit.x - trojanOrbit.x)**2 + (trojanOrbit.y - trojanOrbit.y)**2) )


"""
Calculate the angle from x axis
Function to calculate angle as need to apply corrections depending on which quadrant it is on due to how arctan function works
"""
def angleCalc(x1,y1):
    #First move origin of coordinates to the sun
    deltax = x1
    deltay = y1
    #Theta calculated seperately in each quadrant due to arctan function
    theta = []
    for i in range(len(jupiterOrbit.x)):
        if deltay[i] > 0 and deltax[i] > 0:  # Q1
            theta.append(np.arctan(abs(deltay[i] / deltax[i])))
        elif deltay[i] > 0 and deltax[i] < 0:  # Q2
            theta.append(np.pi - np.arctan(abs(deltay[i] / deltax[i])))
        elif deltay[i] < 0 and deltax[i] < 0:  # Q3
            theta.append(np.pi + np.arctan(abs(deltay[i] / deltax[i])))
        elif deltay[i] < 0 and deltax[i] > 0:  # Q4
            theta.append(2 * np.pi - np.arctan(abs(deltay[i] / deltax[i])))
        else:
            theta.append(np.arctan(deltay[i] / deltax[i]))
    theta = np.array(theta)
    return theta

#Use above function to find angles for each orbits
jupAngle = angleCalc(jupiterOrbit.x, jupiterOrbit.y)
greekAngle = angleCalc(greekOrbit.x, greekOrbit.y)


"""
Calculate "Ideal" Lagrangian Point
"""
greekComIdeal = fromCom(radius, sunInitialdis, jupInitialdis, np.pi/3)
#(Initial distance of Ideal Point from CoM , Initial angle of Ideal Point from CoM)

#Create object to store infromation about Ideal Orbit
#The points are always ideal angle away from Jupiter's Orbit
greekIdealOrbit = Bodies(np.cos(greekComIdeal[1]+jupAngle)*greekComIdeal[0],
                    np.sin(greekComIdeal[1]+jupAngle)*greekComIdeal[0],
                    0,
                    0)

#Create object to store distance infromation for Ideal Orbit
disGreekIdeal =  Distances(np.sqrt(greekIdealOrbit.x**2 + greekIdealOrbit.y**2),
                       np.sqrt((greekIdealOrbit.x - sunOrbit.x)**2 + (greekIdealOrbit.y - sunOrbit.y)**2),
                       np.sqrt((greekIdealOrbit.x - jupiterOrbit.x)**2 + (greekIdealOrbit.y - jupiterOrbit.y)**2),
                       np.sqrt((greekIdealOrbit.x - greekOrbit.x)**2 + (greekIdealOrbit.y - greekOrbit.y)**2),
                       np.sqrt((greekIdealOrbit.x - trojanOrbit.x)**2 + (greekIdealOrbit.y - trojanOrbit.y)**2) )

#Difference in x and y from ideal orbit
disFromIdealx = greekOrbit.x - greekIdealOrbit.x
disFromIdealy = greekOrbit.y - greekIdealOrbit.y


#Calculate angle differences between Jupiter and asteriods from Origin
greekJupAngle= []     #List to store values

for i in range(len(jupAngle)):
    if abs(greekAngle[i] - jupAngle[i]) > np.pi:
        #Need to account for when a body finishes orbit and angle is reset to 0
        greekJupAngle.append(greekAngle[i] + 2*np.pi - jupAngle[i])
    else:
        greekJupAngle.append(greekAngle[i] - jupAngle[i])

greekJupAngle = np.array(greekJupAngle)



"""
Transform into rotating reference frame
Create mew bodies to store relevant information
"""
newFrameJupiterOrbit = Bodies(jupiterOrbit.x *np.cos(angularVelocity*t) + jupiterOrbit.y *np.sin(angularVelocity*t),
                        -jupiterOrbit.x *np.sin(angularVelocity*t) + jupiterOrbit.y *np.cos(angularVelocity*t),
                        jupiterOrbit.vx *np.cos(angularVelocity*t) - jupiterOrbit.x * angularVelocity * np.sin(angularVelocity*t)
                            + jupiterOrbit.vy*np.sin(angularVelocity*t) + jupiterOrbit.y * angularVelocity * np.cos(angularVelocity*t),
                        -jupiterOrbit.vx *np.sin(angularVelocity*t) - jupiterOrbit.x * angularVelocity *np.cos(angularVelocity*t)
                            + jupiterOrbit.vy*np.cos(angularVelocity*t) - jupiterOrbit.y * angularVelocity * np.sin(angularVelocity*t))

newFrameSunOrbit = Bodies(sunOrbit.x *np.cos(angularVelocity*t) + sunOrbit.y *np.sin(angularVelocity*t),
                        -sunOrbit.x *np.sin(angularVelocity*t) + sunOrbit.y *np.cos(angularVelocity*t),
                        sunOrbit.vx *np.cos(angularVelocity*t) - sunOrbit.x * angularVelocity * np.sin(angularVelocity*t)
                            + sunOrbit.vy*np.sin(angularVelocity*t) + sunOrbit.y * angularVelocity * np.cos(angularVelocity*t),
                        -sunOrbit.vx *np.sin(angularVelocity*t) - sunOrbit.x * angularVelocity *np.cos(angularVelocity*t)
                            + sunOrbit.vy*np.cos(angularVelocity*t) - sunOrbit.y * angularVelocity * np.sin(angularVelocity*t))

newframeGreekOrbit = Bodies(greekOrbit.x *np.cos(angularVelocity*t) + greekOrbit.y *np.sin(angularVelocity*t),
                        -greekOrbit.x *np.sin(angularVelocity*t) + greekOrbit.y *np.cos(angularVelocity*t),
                        greekOrbit.vx *np.cos(angularVelocity*t) - greekOrbit.x * angularVelocity * np.sin(angularVelocity*t)
                            + greekOrbit.vy*np.sin(angularVelocity*t) + greekOrbit.y * angularVelocity * np.cos(angularVelocity*t),
                        -greekOrbit.vx *np.sin(angularVelocity*t) - greekOrbit.x * angularVelocity *np.cos(angularVelocity*t)
                            + greekOrbit.vy*np.cos(angularVelocity*t) - greekOrbit.y * angularVelocity * np.sin(angularVelocity*t))


"""
Calculate Max and Min values"""

maxAngle = np.max(greekJupAngle) / np.pi
minAngle = np.min(greekJupAngle) / np.pi
midAngle = (maxAngle + minAngle) / 2
maxDis = np.max(disGreek.origin)
minDis = np.min(disGreek.origin)
midDis = (maxDis + minDis) / 2







"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PLOTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""


#Figure 1
#Plots of x and y coords for all orbits
fig = plt.figure(figsize=(10,6))
ax = plt.axes(xlim=(-6, 6), ylim=(-6, 6))
plt.axis('square')

plt.plot(jupiterOrbit.x,jupiterOrbit.y)                 #Jupiter's Orbit
#plt.plot(sunOrbit.x,sunOrbit.y)                        #Sun's Orbit
plt.plot(greekOrbit.x,greekOrbit.y, "r-")               #Greek's Orbit
#plt.plot(trojanOrbit.x,trojanOrbit.y,'b--')            #Trojan's Orbit
#plt.plot(greekIdealOrbit.x,greekIdealOrbit.y,'r--')    #Ideal Greek's Orbit



numberofOss= t/jupPeriod    #Number of Orbits
#Figure 2
#Plot of Distance from Jupiter and Sun

fig2, ax2 = plt.subplots(2, 1, constrained_layout=True)
ax2[0].plot(numberofOss, disGreek.jupiter, 'r-', label="Greek")     #Distance of Greek from Jupiter
ax2[0].plot(numberofOss, disTrojan.jupiter, 'b--', label="Trojan")  #Distance of Trojan fro Jupiter

ax2[0].set_title('Distance from Jupiter')
ax2[0].set_xlabel('Number of Orbits')
ax2[0].set_ylabel('Distance (AU)')
ax2[0].legend(loc="upper left")
#ax2[0].set_ylim(5.19,5.21)


ax2[1].plot(numberofOss, disGreek.sun, 'r-', label="Greek")         #Distance of Greek from Sun
ax2[1].plot(numberofOss, disTrojan.sun, 'b--', label="Trojan")      #Distance of Trojan from Sun

ax2[1].set_title('Distance from Sun')
ax2[1].set_xlabel('Number of Orbits')
ax2[1].set_ylabel('Distance (AU)')
ax2[1].legend(loc="upper left")
#ax2[1].set_ylim(5.19,5.21)


#Figure 3
#Plot of Anglular Displacement from Jupiter and Radius from CoM

fig3, ax3 = plt.subplots(2, 1, constrained_layout=True)
ax3[0].plot(numberofOss,  greekJupAngle/np.pi, 'r-', label="Perturbed")         #Greek's Angular Displacement
ax3[0].plot(numberofOss,
            [greekComIdeal[1]/np.pi]*len(numberofOss),                          #Ideal Greek's Angular Displacement
            'b--', label="Ideal")

ax3[0].set_title("Angular displacement between Asteroids and Jupiter from CoM")
ax3[0].set_xlabel("Number of Orbits")
ax3[0].set_ylabel("Angular Displacement / π ")
ax3[0].legend(loc="upper left")
#ax3[0].set_ylim(0.32,0.35)


ax3[1].plot(numberofOss, disGreek.origin, 'r-', label="Perturbed")              #Greek's Radius
ax3[1].plot(numberofOss, disGreekIdeal.origin, 'b--', label="Ideal")            #Ideal Greek's Radius

ax3[1].set_title('Radial Distance from CoM')
ax3[1].set_xlabel('Number of Orbits')
ax3[1].set_ylabel('Distance (AU)')
ax3[1].legend(loc="upper left")
#ax3[1].set_ylim(5.19,5.21)


#Figure 4
#Plot of Perturbations around Lagrangian Point

fig4, ax4 = plt.subplots(2,1,constrained_layout=True)
ax4[0].plot(disFromIdealx[0], disFromIdealy[0], 'rx',ms="10", label="Initial Position")    #Initial Position
ax4[0].plot(disFromIdealx, disFromIdealy, 'b', label="Perturbed Orbit")        #Difference in x and y from ideal orbit
ax4[0].set_title('Position relative to The Lagrangian Point')
ax4[0].set_xlabel('x (AU)')
ax4[0].set_ylabel('y (AU)')
ax4[0].legend(loc="upper left")
ax4[0].axis('square')


ax4[1].plot(newframeGreekOrbit.x[0]-greekIdealOrbit.x[0],
            newframeGreekOrbit.y[0]-greekIdealOrbit.y[0], 'bx',ms="10", label="Initial Position")
ax4[1].plot(newframeGreekOrbit.x-greekIdealOrbit.x[0],
            newframeGreekOrbit.y-greekIdealOrbit.y[0],'r-', label="Perturbed Orbit")    #Difference in x and y from ideal position
ax4[1].plot(greekIdealOrbit.x-greekIdealOrbit.x[0],
            greekIdealOrbit.y-greekIdealOrbit.y[0],'b-', label="Ideal Orbit in normal frame")
ax4[1].set_title('Position relative to The Lagrangian Point in Rotating Frame')
ax4[1].set_xlabel('x (AU)')
ax4[1].set_ylabel('y (AU)')
ax4[1].legend(loc="lower left")
#ax4[1].set_ylim(-0.5,0.5)
#ax4[1].set_xlim(-0.9,0.7)
ax4[1].axis('square')



"""
Animation
"""
#Animation to help visualise orbits

pointJup, = ax.plot(jupiterOrbit.x[0], jupiterOrbit.y[0], 'o',mfc="brown",ms="10")
pointSun, = ax.plot(sunOrbit.x[0], sunOrbit.y[0], 'yo',  ms="20")
pointsGreek, = ax.plot(greekOrbit.x[0], greekOrbit.y[0], 'ro')
pointsTrojan, = ax.plot(trojanOrbit.x[0], trojanOrbit.y[0], 'bo')
idealGreek, = ax.plot(greekIdealOrbit.x[0],greekIdealOrbit.y[0],'o')


def animatejup(i):
    speed=1
    pointJup.set_ydata(jupiterOrbit.y[speed*i])
    pointJup.set_xdata(jupiterOrbit.x[speed*i])

    pointSun.set_ydata(sunOrbit.y[speed*i])
    pointSun.set_xdata(sunOrbit.x[speed*i])

    pointsGreek.set_xdata(greekOrbit.x[speed*i])
    pointsGreek.set_ydata(greekOrbit.y[speed*i])

    pointsTrojan.set_xdata(trojanOrbit.x[speed*i])
    pointsTrojan.set_ydata(trojanOrbit.y[speed*i])

    idealGreek.set_xdata(greekIdealOrbit.x[speed*i])
    idealGreek.set_ydata(greekIdealOrbit.y[speed*i])


    return pointSun, pointJup, pointsGreek, pointsTrojan, idealGreek

anim = animation.FuncAnimation(fig, animatejup,
                               frames=2000, interval=15)

plt.show()