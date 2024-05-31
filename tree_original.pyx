import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random 
from tree_original import Body, rect, QuadtreeNode, compute_force, update_per_body, find_boundary



cimport cython
from cython.parallel import parallel, prange
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random


# Define classes here...

@cython.boundscheck(False)
@cython.wraparound(False)
def update_per_dt(int frame):
    cdef list x_list = []
    cdef list y_list = []
    cdef float x_max, y_max
    x_max, y_max = find_boundary(body_list)
    cdef QuadtreeNode qt = QuadtreeNode(rect(-x_max, -y_max, 2*x_max, 2*y_max), 0)

    with nogil, parallel():
        for i in prange(len(body_list)):
            qt.insert(body_list[i])

    qt._update_mass_distribution()

    with nogil, parallel():
        for i in prange(len(body_list)):
            body_list[i] = update_per_body(body_list[i], qt)
            x_list.append(body_list[i].x)
            y_list.append(body_list[i].y)
    
def init():
    x_list = []
    y_list = []
    for body in body_list:
        x_list.append(body.x)
        y_list.append(body.y)
    ax.plot(x_list, y_list, "bo", ms=5)

if __name__ == "__main__":
    t = 0.0
    dt = 1e-1
    end_time = 10.0
    n = 100
    body_list = []
    pos = np.linspace(-50, 50, 200)
    vel = np.linspace(-1, 1, 40)
#    body_list = [Body(-0.9,0.9,1),Body(-1,-1,1),Body(1,1,1),Body(1,-1,1),Body(-0.8,0.5,1),Body(-0.5,0.8,1),
#                 Body(-0.4,0.7,1),Body(-1.5,0.5,1),Body(-1.5,0.3,1)]
    for i in range(n):
        x_ini = np.random.choice(pos, 1)[0]
        y_ini = np.random.choice(pos, 1)[0]
        vx_ini = np.random.choice(vel, 1)[0]
        vy_ini = np.random.choice(vel, 1)[0]
        body_list.append(Body(x_ini, y_ini, 1, vx_ini, vy_ini))
    fig, ax = plt.subplots()
    nframe = int( np.ceil( end_time/dt ) )
    anim   = animation.FuncAnimation( fig, func=update_per_dt, init_func=init,
                                  frames=nframe, interval=10, repeat=False )

    plt.show()
