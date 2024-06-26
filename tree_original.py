import math
import numpy as np
import matplotlib.pyplot as plt

class Body:
    def __init__(self, mass, x, y, vx, vy):
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
            if (x_b < x or x_b > x+w or y_b < y or y_b > y+h):
                return True

    def insert(self, body):
        if self.check_boundary(body):
            return False
        if self.body == None:
            self.body = body
            print(self.body)
            return True
        if len(self.area) <= 0:
            self.subdivide()
            print(len(self.area))
            for area_ in self.area:
                print(self.body)
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

    def _update_mass_distribution(self):
        self.mass = 0
        self.center_of_mass_x = 0
        self.center_of_mass_y = 0
        if self.is_leaf():
            if self.body:
                self.mass = self.body.mass
                self.center_of_mass_x = self.body.x
                self.center_of_mass_y = self.body.y
        else:
            for child in self.children:
                if child:
                    self.mass += child.mass
                    self.center_of_mass_x += child.mass * child.center_of_mass_x
                    self.center_of_mass_y += child.mass * child.center_of_mass_y
            if self.mass > 0:
                self.center_of_mass_x /= self.mass
                self.center_of_mass_y /= self.mass
    
    def plot(self, ax):
        x = self.boundary.x
        y = self.boundary.y
        h = self.boundary.h
        w = self.boundary.w
        ax.plot([x, x], [y+h, y], color="black")
        ax.plot([x, x+w], [y, y], color="black")
        ax.plot([x+w, x+w], [y+h, y], color="black")
        ax.plot([x, x+w], [y+h, y+h], color="black")
        if len(self.area) != 0:
            for area_ in self.area:
                area_.plot(ax)
        if self.body != None and self.body != "empty":
            ax.plot([self.body.x],[self.body.y],"bo", ms=5)

def compute_force(body, node, theta=0.5, G=6.67430e-11):
    if node.is_leaf():
        if node.body is not None and node.body != body:
            dx = node.body.x - body.x
            dy = node.body.y - body.y
            distance = math.sqrt(dx**2 + dy**2)
            if distance > 0:
                force = G * body.mass * node.body.mass / (distance**2)
                body.ax += force * dx / distance
                body.ay += force * dy / distance
    else:
        dx = node.center_of_mass_x - body.x
        dy = node.center_of_mass_y - body.y
        distance = math.sqrt(dx**2 + dy**2)
        size = node.x_max - node.x_min
        if size / distance < theta:
            force = G * body.mass * node.mass / (distance**2)
            body.ax += force * dx / distance
            body.ay += force * dy / distance
        else:
            for child in node.children:
                if child:
                    compute_force(body, child, theta, G)

def update_orbit( body, dt ):
    body.x
    # Drift
    body.x += body.vx * dt/2
    body.y += body.vy * dt/2
    # Kick
    body.ax = compute_force(body)/body.mass # compute_force function is still incorrrect
    body.ay = compute_force(body)/body.mass
    body.vx += body.ax * dt
    body.vy += body.ay * dt
    # Drift
    body.x += body.vx * dt/2
    body.y += body.vy * dt/2 

init_body = np.loadtxt('IC16.txt',comments='#')

if __name__ == "__main__":
    qt = QuadtreeNode(rect(0,0,3,3),0)
    for particle_num in range(16):
        qt.insert(Body(init_body[particle_num,0],init_body[particle_num,1],init_body[particle_num,2],init_body[particle_num,4],init_body[particle_num,5]))
    
    fig, ax = plt.subplots()
    qt.plot(ax)

    plt.show()
