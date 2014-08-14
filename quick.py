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
1. customized motion

lower priority:
2. allow multiple objects
3. off-the-track motion
4. command line argument support
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
        #for c in ["#f16745", "#ffc65d", "#7bc8a4", "#93648d", "#404040"]:
        for c in ["#FF7D7D", "#88FF7D", "#FFFC3D", "#7DFFF8", "#D8BAFF"]:
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
        axes.plot(*self.path, color="#EEEEEE", linewidth=2.5)

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

class Phyobj(Obj_base):
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
        self.vel = rotate(oldvel, midvel)
        return midvel

    def accelerate(self, unitvec, gravity, dt):
        """increase the velocity
        """
        # follwing is not most efficient, but more redable.
        accel = np.dot(gravity, unitvec) * unitvec  # project vector
        self.vel += accel * dt
        
class Merry(Phyobj):
    def __init__(self, mass = 1, init_vel=(0.0,0.0), xy = (0,0), radius = 0.15, **kw):
        Phyobj.__init__(self, mass, init_vel)
        self.center = xy
        self.ballprops = {"orbit_radius" : 0,
                          "angle" : 0,
                          "balls" : [ Circle(xy, radius, color=next(colorgen))
                                      for x in range(3) ]
                         }

    def patches(self):
        return self.ballprops["balls"]

    def show(self, axes):
        for b in self.ballprops["balls"]:
            axes.add_patch(b)

    def anglechanged(self, dt):
        return (np.linalg.norm(self.vel) * 6 * dt) ** 1.3

    def genradius(self):
        return np.linalg.norm(self.vel) * 0.07

    def getpos(self):
        return self.center

    def setpos(self, xy):
        self.center = xy

    def custom_animate(self, dt):
        self.ballprops["orbit_radius"] = self.genradius()
        self.ballprops["angle"] += self.anglechanged(dt)
        for i, b in enumerate(self.ballprops["balls"]):
            angle = self.ballprops["angle"] + (2 * math.pi * i/3)
            offset = (self.ballprops["orbit_radius"] *
                         np.array([math.cos(angle), math.sin(angle)]))
            b.center = self.getpos() + offset

class System():
    """The system is at rest relative to the graph
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
        #axes.add_patch(self.obj)
        self.obj.show(axes)

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
        # print(bestpos)
        midvel = self.obj.moveto(bestpos, dt)

        if any(midvel):  # meaning "if it is not null vector", preventing nan
            unitvec = midvel/np.linalg.norm(midvel)  # get the slope as vaetor
        else:
            tanvec = self.runway.tanvector(self.p)
            unitvec = tanvec/np.linalg.norm(tanvec)
        self.obj.accelerate(unitvec, self.gravity, dt)

        if hasattr(self.obj, "custom_animate"):
            self.obj.custom_animate(dt)
        return self.obj.patches()

    def multistep(self, dt, n):
        for i in range(n):
            self.step(dt)
        plt.draw()  # only used in interactive mode


if sys.argv[0] == "":
    plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111, axisbg="#222222") 
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.axis("equal")

# This is a roller coaster track
myrunway = Runway(
    np.arange(-5, 5, 0.1),
    fx=lambda p: -p/2.8 + math.sin(p),
    fy=lambda p: (p**2)/9 +math.cos(p),
    derx=lambda p: -0.5+ math.cos(p),
    dery=lambda p: 2*p/9 - math.sin(p))
'''
myrunway = Runway(
    np.arange(0, 2*3.14, 0.1),
    fx=lambda p: 2.5*math.sin(p),
    fy=lambda p: 3*math.cos(p))
    '''

    
cir = Merry(init_vel=np.array([0.0,0.0]))

mysystem = System(cir, myrunway, (0.0, -2.5))  #On Mars?
mysystem.plot(ax)

dt = 0.005
def animate(i):
    mysystem.multistep(dt, 8)
    return mysystem.obj.patches()

anim = animation.FuncAnimation(fig, animate, frames=250, interval=40,
        blit=True, repeat_delay=3000)
anim.save("roll.mp4", fps=40, extra_args=['-vcodec', 'libx264'])
#plt.show()
