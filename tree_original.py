import math

class Body:
    def __init__(self, x, y, mass):
        self.x = x
        self.y = y
        self.mass = mass
        self.vx = 0
        self.vy = 0
        self.ax = 0
        self.ay = 0

class QuadtreeNode:
    def __init__(self, x_min, x_max, y_min, y_max):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.body = None
        self.mass = 0
        self.center_of_mass_x = 0
        self.center_of_mass_y = 0
        self.children = [None, None, None, None]

    def is_leaf(self):
        return all(child is None for child in self.children)

    def insert(self, body):
        if self.is_leaf():
            if self.body is None:
                self.body = body
            else:
                self.subdivide()
                self._insert_into_children(self.body)
                self._insert_into_children(body)
                self.body = None
        else:
            self._insert_into_children(body)
        self._update_mass_distribution()

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

# Example usage:
bodies = [Body(x, y, mass) for x, y, mass in [(0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)]]
root = QuadtreeNode(-2, 2, -2, 2)
for body in bodies:
    root.insert(body)
for body in bodies:
    compute_force(body, root)
    print(f"Body at ({body.x}, {body.y}) has acceleration ({body.ax}, {body.ay})")
