#!/usr/bin/python3

"""
Show a ball rolling down a roller coaster runway
The interface defined can adapt to any runway (maybe 3-D, but need some
adjustment), so I may move it to a library later.

The algorithm is:  there is gravity, and normal reaction, make the ball move,
first assume no constraint of runway, then assumetthere is, find the shortest
point to the estimated point, and move the ball there, update velocity, and go
on ....

First order approximation is junk? I don't care, it is good enough.
with 0.001 step, falling from 3.06, the ball end up in 2.83 o the left end.
7.5% error is good enough.

From height 3.06, 0.01 dt step, the ball end in 3.11

TODO:
==========
0. prettier graphic
1. use abc metaclass for moving objects to ensure they implement the nessary
method.
2. customized motion

lower priority:
3. allow multiple objects
4. off-the-track motion
5. command line argument support
999. energy conservation?
"""
import sys
import math
import numpy as np
import scipy
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from  scipy.optimize import leastsq

def rotate(vec, to_vec):
    if not any(to_vec):
        return vec
    else:
        return np.linalg.norm(vec) * to_vec / np.linalg.norm(to_vec)

def colorfactory():
    while True:
        for c in ["#f16745", "#ffc65d", "#7bc8a4", "#93648d", "#404040"]:
            yield c

colorgen = colorfactory()

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
        axes.plot(*self.path, color="#444444", linewidth=1.5)

    # There is no use of this function yet.
    def tanvector(self, p):
        return np.array([self.derx(p), self.dery(p)])


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
    def __init__(self, obj, runway, init_vel=(0,0), gravity=(0,-9.8)):
        self.obj = obj
        self.runway = runway
        self.vel = np.asarray(init_vel)
        self.veldiff = (np.asarray((0, 0)), 1)  #velocity change over time (1)
        self.obj.setpos(np.asarray((self.runway.path[0][0],
                                    self.runway.path[1][0])))
        self.p = self.runway.prange[0]
        self.gravity = gravity

    def plot(self, axes):
        self.runway.plot(axes)
        axes.add_patch(self.obj)

    def vectdistance(self, p, pos):
        """Return the vector form the runway to a point pos as a tuple
        """
        return (self.runway.fx(p) - pos[0], self.runway.fy(p) - pos[1])

    def accelerate(self, obj, unitvec, dt):
        """increase the velocity
        """
        # obj is not used here, but some day when there is multiple obj
        # it will be used.
        # follwing is not most efficient, but more redable.
        accel = np.dot(self.gravity, unitvec) * unitvec  # project vector
        self.vel += accel * dt
        
    def step(self, dt):
        """evolve the system by time dt
        """
        estimatepos = self.obj.getpos() + self.vel * dt
        #Find p which give shortest distance between runway and setimated point
        res = scipy.optimize.leastsq(
                lambda p: self.vectdistance(p[0], estimatepos),
                [self.p])

        self.p = res[0][0]
        bestpos = np.asarray((self.runway.fx(self.p), self.runway.fy(self.p)))
        midvel = (bestpos - self.obj.getpos()) / dt
        self.obj.setpos(bestpos)
        oldvel = self.vel

        # rotate vel to the guessed current velocity direction
        self.vel = rotate(self.vel, midvel + (midvel - oldvel))

        if any(midvel):  # meaning "if it is not null vector", preventing nan
            unitvec = midvel/np.linalg.norm(midvel)  # get the slope as vaetor
        else:
            tanvec = self.runway.tanvector(self.p)
            unitvec = tanvec/np.linalg.norm(tanvec)
        self.accelerate(self.obj, unitvec, dt)

        if hasattr(self.obj, "custom_animate"):
            self.obj.custom_animate()
        return self.obj

    def multistep(self, dt, n):
        for i in range(n):
            self.step(dt)
            #time.sleep(0.2)
            #plt.draw()
            #input()
        plt.draw()


if sys.argv[0] == "":
    plt.ion()
    import time
fig = plt.figure(figsize=(1,2))
ax = fig.add_subplot(111) 
ax.axis("equal")

# This is a roller coaster track
myrunway = Runway(
    lambda p: -p/2 + math.sin(p),
    lambda p: -0.5+ math.cos(p),
    lambda p: (p**2)/9 +math.cos(p),
    lambda p: 2*p/9 - math.sin(p),
    np.arange(-5, 5, 0.1))

cir = Ball(color=next(colorgen))

mysystem = System(cir, myrunway, (0.0, 0.0), (0.0, -3.0))  #On Mars?
mysystem.plot(ax)
"""
dt = 0.01
def animate(i):
    mysystem.step(dt)
anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=200, interval=20, blit=True)
plt.draw()
"""
