#!/usr/bin/python3
import sys
import math
import numpy as np
import scipy
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from  scipy.optimize import leastsq

if sys.argv[0] == "":
    plt.ion()

class Runway():
    def __init__(self, fx, derx, fy, dery, prange):
        self.fx = fx
        self.derx = derx
        self.fy = fy
        self.dery = dery
        self.prange = prange
        self.genpath()

    def genpath(self):
        x = [ self.fx(p) for p in self.prange ]
        y = [ self.fy(p) for p in self.prange ]
        self.path = (x,y)

    def plot(self, axes):
        axes.plot(*self.path)

    def tanvector(self, p):
        vect = np.array([self.derx(p), self.dery(p)])
        return vect/np.linalg.norm(vect)


class Phyobj():
    def __init__(self, mass):
        self.mass = mass

class Ball(Circle, Phyobj):
    def __init__(self, mass = 1, xy = 1, radius = 0.1, **kw):
        Phyobj.__init__(self, mass)
        Circle.__init__(self, xy, radius, **kw)
    def getpos(self):
        return self.center
    def setpos(self, xy):
        self.center = xy

class System():
    """A system to evolve with a runway and a object
    velocity is dependent on frame, so it should be defined here.
    """
    def __init__(self, obj, runway, init_vel):
        self.obj = obj
        self.runway = runway
        self.vel = init_vel  # this dhould be a np array, representing vector
        self.obj.setpos(np.asarray((self.runway.path[0][0],
                                    self.runway.path[1][0])))
        self.p = self.runway.prange[0]

    def plot(self, axes):
        self.runway.plot(axes)
        axes.add_patch(self.obj)

    def vectdistance(self, p, pos):
        """Return the vector form the runway to a point pos as a tuple
        """
        return (self.runway.fx(p) - pos[0], self.runway.fy(p) - pos[1])

    def step(self, dt):
        estimatepos = self.obj.getpos() + self.vel * dt
        print(estimatepos)
        #Find p which give shortest distance between runway and setimated point
        res = scipy.optimize.leastsq(
                lambda p: self.vectdistance(p[0], estimatepos),
                [self.p])
        print(res)
        self.p = res[0][0]
        bestpos = np.asarray((self.runway.fx(self.p), self.runway.fy(self.p)))
        self.vel = (bestpos - self.obj.getpos()) / dt
        self.obj.setpos(bestpos)
        if hasattr(self.obj, "custom_animate"):
            self.obj.custom_animate()
        return self.obj

fig = plt.figure(figsize=(1,2))
ax = fig.add_subplot(111) 

# This is a roller coaster track
myrunway = Runway(
    lambda p: -p/2 + math.sin(p),
    lambda p: -0.5+ math.cos(p),
    lambda p: (p**2)/9 +math.cos(p),
    lambda p: 2*p/9 - math.sin(p),
    np.arange(-5, 5, 0.1))

cir = Ball(color="#ff9197")

mysystem = System(cir, myrunway, (-0.03, -0.03))
mysystem.plot(ax)
"""
dt = 0.01
def animate(i):
    mysystem.step(dt)
anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=200, interval=20, blit=True)
plt.draw()
"""
