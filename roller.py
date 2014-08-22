#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import runwaylib
import math
import sys

fig = plt.figure()
ax = fig.add_subplot(111, axisbg="#222222") 
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.axis("equal")

# This is a roller coaster track
myrunway = runwaylib.Runway(
    np.arange(-5, 5, 0.1),
    fx=lambda p: -p/2.8 + math.sin(p),
    fy=lambda p: (p**2)/9 + math.cos(p),
    derx=lambda p: -0.5 + math.cos(p),
    dery=lambda p: 2*p/9 - math.sin(p))

cir = runwaylib.Merry(init_vel=np.array([0.0,0.0]))

mysystem = runwaylib.System(cir, myrunway, (0.0, -2.5))  #On Mars?
mysystem.plot(ax)

dt = 0.005
def animate(i):
    mysystem.multistep(dt, 8)
    return mysystem.obj.patches()

anim = animation.FuncAnimation(fig, animate, frames=250, interval=40,
        blit=True, repeat_delay=3000)
anim.save("roll.mp4", fps=40, extra_args=['-vcodec', 'libx264'])

if sys.argv[0] == "":
    plt.ion()
    # plt.show()    # for debug use
