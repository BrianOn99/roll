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
from abc import ABCMeta, abstractmethod
import math
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle

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

def notimplement(s):
    raise NotImplementedError(s)

class Runway():
    def __init__(self, prange, fx, fy,
                 derx=lambda: notimplement("derx"),
                 dery=lambda: notimplement("dery")):
        self.fx = fx
        self.fy = fy
        self.derx = derx
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

class Obj_base(metaclass=ABCMeta):
    @abstractmethod
    def getpos(self):
        raise NotImplementedError
    @abstractmethod
    def setpos(self):
        raise NotImplementedError

class Phyobj():
    def __init__(self, mass, init_vel=(0.0,0.0)):
        super().__init__()
        self.mass = mass
        self.vel = np.asarray(init_vel)  #velocity relative to the graph
        self.veldiff = (np.asarray((0, 0)), 1)

    def trymove(self, dt):
        """try to move only, wont really move
        """
        return self.getpos() + self.vel * dt

    def moveto(self, pos, dt):
        midvel = (pos - self.getpos()) / dt
        self.setpos(pos)
        oldvel = self.vel
        # rotate vel to the guessed current velocity direction
        self.vel = rotate(oldvel, midvel + (midvel - oldvel))
        return midvel

    def accelerate(self, unitvec, gravity, dt):
        """increase the velocity
        """
        # follwing is not most efficient, but more redable.
        accel = np.dot(gravity, unitvec) * unitvec  # project vector
        self.vel += accel * dt
        
class Ball(Circle, Phyobj, Obj_base):
    def __init__(self, mass = 1, init_vel=(0.0,0.0), xy = (0,0), radius = 0.1, **kw):
        Phyobj.__init__(self, mass, init_vel)
        Circle.__init__(self, xy, radius, **kw)
    def getpos(self):
        return self.center
    def setpos(self, xy):
        self.center = xy

class System():
    """A system to evolve with a runway and a object
    velocity is dependent on frame, so it should be defined here.
    The system is at rest relative to the graph
    """
    def __init__(self, obj, runway, gravity=(0,-9.8)):
        self.runway = runway
        self.obj = obj
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

        
    def step(self, dt):
        """evolve the system by time dt
        """
        estimatepos = self.obj.trymove(dt)
        #Find p which give shortest distance between runway and setimated point
        res = scipy.optimize.leastsq(
                lambda p: self.vectdistance(p[0], estimatepos),
                [self.p])

        self.p = res[0][0]
        bestpos = np.asarray((self.runway.fx(self.p), self.runway.fy(self.p)))
        print(bestpos)
        midvel = self.obj.moveto(bestpos, dt)

        if any(midvel):  # meaning "if it is not null vector", preventing nan
            unitvec = midvel/np.linalg.norm(midvel)  # get the slope as vaetor
        else:
            tanvec = self.runway.tanvector(self.p)
            unitvec = tanvec/np.linalg.norm(tanvec)
        self.obj.accelerate(unitvec, self.gravity, dt)

        if hasattr(self.obj, "custom_animate"):
            self.obj.custom_animate()
        return self.obj

    def multistep(self, dt, n):
        for i in range(n):
            self.step(dt)
        plt.draw()  # only used in interactive mode


if sys.argv[0] == "":
    plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111) 
ax.axis("equal")

# This is a roller coaster track
myrunway = Runway(
    np.arange(-5, 5, 0.1),
    fx=lambda p: -p/2 + math.sin(p),
    fy=lambda p: (p**2)/9 +math.cos(p),
    derx=lambda p: -0.5+ math.cos(p),
    dery=lambda p: 2*p/9 - math.sin(p))
    

cir = Ball(color=next(colorgen))

mysystem = System(cir, myrunway, (0.0, -3.0))  #On Mars?
mysystem.plot(ax)

dt = 0.01
def animate(i):
    mysystem.multistep(dt, 5)
    return (mysystem.obj,)

anim = animation.FuncAnimation(fig, animate, frames=60, interval=100,
        blit=True, repeat_delay=1000)
# anim.save("roll.mp4")
plt.show()
