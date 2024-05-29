import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random

class Body:
    def __init__(self, x, y, mass, vx, vy):
        self.x = x
        self.y = y
        self.mass = mass
        self.vx = vx
        self.vy = vy
        self.ax = 0
        self.ay = 0

class rect:
    def __init__(self, x, y, h, w):
        self.x = x
        self.y = y
        self.h = h
        self.w = w

class QuadtreeNode:
    def __init__(self, boundary, level):
        self.max_obj = 1
        self.body = None
        self.area = []
        self.level = level
        self.boundary = boundary
    
    def check_boundary(self, body):
        x_b, y_b = body.x, body.y
        x, y, h, w = self.boundary.x, self.boundary.y, self.boundary.h, self.boundary.w
        if self.level != 0:
            if x_b < x or x_b > x+w or y_b < y or y_b > y+h:
                return True

    def insert(self, body):
        if self.check_boundary(body):
            return False
        if self.body == None:
            self.body = body
            return True
        if len(self.area) <= 0:
            self.subdivide()
            for area_ in self.area:
                if area_.insert(self.body):
                    self.body = "empty"
                    break
        for area_ in self.area:
            if area_.insert(body):
                return True

    def subdivide(self):
        x = self.boundary.x
        y = self.boundary.y
        h = self.boundary.h
        w = self.boundary.w

        self.area.append(QuadtreeNode(rect(x+w/2, y, h/2, w/2), self.level+1))
        self.area.append(QuadtreeNode(rect(x, y, h/2, w/2), self.level+1))
        self.area.append(QuadtreeNode(rect(x, y+h/2, h/2, w/2), self.level+1))
        self.area.append(QuadtreeNode(rect(x+w/2, y+h/2, h/2, w/2), self.level+1))

    def _insert_into_children(self, body):
        for i, (x_min, x_max, y_min, y_max) in enumerate(self._get_quadrants()):
            if x_min <= body.x < x_max and y_min <= body.y < y_max:
                if self.children[i] is None:
                    self.children[i] = QuadtreeNode(x_min, x_max, y_min, y_max)
                self.children[i].insert(body)
                break

    def _get_quadrants(self):
        x_mid = (self.x_min + self.x_max) / 2
        y_mid = (self.y_min + self.y_max) / 2
        return [
            (self.x_min, x_mid, self.y_min, y_mid),
            (x_mid, self.x_max, self.y_min, y_mid),
            (self.x_min, x_mid, y_mid, self.y_max),
            (x_mid, self.x_max, y_mid, self.y_max),
        ]

    def calculate_center_of_mass(self, mass = 0, center_of_mass_x = 0, center_of_mass_y = 0):
        if self.body == "empty":
            for area_ in self.area:
                mass, center_of_mass_x, center_of_mass_y = area_.calculate_center_of_mass(mass, center_of_mass_x, center_of_mass_y)
                
        elif self.body != None:
            mass += self.body.mass
            center_of_mass_x += self.body.mass*self.body.x
            center_of_mass_y += self.body.mass*self.body.y


        return mass, center_of_mass_x, center_of_mass_y

    def _update_mass_distribution(self):
        self.mass = 0
        self.center_of_mass_x = None
        self.center_of_mass_y = None
        if self.body != None:
            self.mass, self.center_of_mass_x, self.center_of_mass_y = self.calculate_center_of_mass()
            if self.mass != 0:
                self.center_of_mass_x /= self.mass
                self.center_of_mass_y /= self.mass
            for area_ in self.area:
                area_._update_mass_distribution()
    
    def plot(self, ax):
        x = self.boundary.x
        y = self.boundary.y
        h = self.boundary.h
        w = self.boundary.w
        ax.plot([x, x], [y+h, y], color="black")
        ax.plot([x, x+w], [y, y], color="black")
        ax.plot([x+w, x+w], [y+h, y], color="black")
        ax.plot([x, x+w], [y+h, y+h], color="black")
        if self.center_of_mass_x != None:
            ax.plot([self.center_of_mass_x],[self.center_of_mass_y], "gx")
        if len(self.area) != 0:
            for area_ in self.area:
                area_.plot(ax)
        if self.body != None and self.body != "empty":
            ax.plot([self.body.x],[self.body.y],"bo", ms=5)

def compute_force(body, node, theta=0.5, G=1e-1):
    if node.body != "empty":
        if node.body is not None and node.body != body:
            dx = node.body.x - body.x
            dy = node.body.y - body.y
            distance = math.sqrt(dx**2 + dy**2)
            if distance > 0:
                force = G * body.mass * node.body.mass / (distance**2)
                body.ax += force * dx / (distance*body.mass)
                body.ay += force * dy / (distance*body.mass)
    else:
        dx = node.center_of_mass_x - body.x
        dy = node.center_of_mass_y - body.y
        distance = math.sqrt(dx**2 + dy**2)
        size = node.boundary.w
        if size / distance < theta:
            force = G * body.mass * node.mass / (distance**2)
            body.ax += force * dx / (distance*body.mass)
            body.ay += force * dy / (distance*body.mass)
        else:
            for area_ in node.area:
                compute_force(body, area_, theta, G)
def update_per_body(body, node):
    ##### update orbit (DKD) ######
    # drift
    body.x = body.x + body.vx*0.5*dt
    body.y = body.y + body.vy*0.5*dt

    # kick
    compute_force(body, node)
    body.vx = body.vx + body.ax*dt
    body.vy = body.vy + body.ay*dt

    # drift
    body.x = body.x + body.vx*0.5*dt
    body.y = body.y + body.vy*0.5*dt
    return body

def find_boundary(body_list):
    x_list = []
    y_list = []
    for body in body_list:
        x_list.append(body.x)
        y_list.append(body.y)
    x_max = max(x_list)
    x_min = min(x_list)
    y_max = max(y_list)
    y_min = min(y_list)
    if abs(x_max) < abs(x_min):
        x_max = abs(x_min)
    if abs(y_max) < abs(y_min):
        y_max = abs(y_min)
    return x_max, y_max

def update_per_dt(frame):
    ax.clear()
    ax.set_xlim([-100,100])
    ax.set_ylim([-100,100])
    x_list = []
    y_list = []
    x_max, y_max = find_boundary(body_list)
    qt = QuadtreeNode(rect(-x_max,-y_max,2*x_max,2*y_max),0)
    for body in body_list:
        qt.insert(body)
    qt._update_mass_distribution()
    for i in range(len(body_list)):
        body_list[i] = update_per_body(body_list[i], qt)
        x_list.append(body_list[i].x)
        y_list.append(body_list[i].y)
    x_max = max(x_list)
    y_max = max(y_list) 
    #print(body_list[1].x)
    ax.plot(x_list, y_list, "bo", ms=5)
    
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
    end_time = 5.0
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
