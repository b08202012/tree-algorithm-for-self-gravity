#include <SFML/Graphics.hpp>
#include <iostream>
#include <vector>
#include <cmath>

const double G = 6.67430e-11; // Gravitational constant
const double TIME_STEP = 0.01; // Time step for the simulation

class Body {
public:
    double x, y, mass, vx, vy, ax, ay;

    Body(double x, double y, double mass) : x(x), y(y), mass(mass), vx(0), vy(0), ax(0), ay(0) {}

    void update() {
        // Update velocity
        vx += ax * TIME_STEP;
        vy += ay * TIME_STEP;
        // Update position
        x += vx * TIME_STEP;
        y += vy * TIME_STEP;
        // Reset acceleration
        ax = 0;
        ay = 0;
    }

    void draw(sf::RenderWindow& window) {
        sf::CircleShape shape(2);
        shape.setPosition(x, y);
        shape.setFillColor(sf::Color::White);
        window.draw(shape);
    }
};

class QuadtreeNode {
public:
    double x_min, x_max, y_min, y_max;
    Body* body;
    double mass, center_of_mass_x, center_of_mass_y;
    QuadtreeNode* children[4];

    QuadtreeNode(double x_min, double x_max, double y_min, double y_max)
        : x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max), body(nullptr), mass(0),
          center_of_mass_x(0), center_of_mass_y(0) {
        for (int i = 0; i < 4; ++i) {
            children[i] = nullptr;
        }
    }

    bool isLeaf() const {
        for (int i = 0; i < 4; ++i) {
            if (children[i] != nullptr) {
                return false;
            }
        }
        return true;
    }

    void subdivide() {
        double x_mid = (x_min + x_max) / 2;
        double y_mid = (y_min + y_max) / 2;
        children[0] = new QuadtreeNode(x_min, x_mid, y_min, y_mid);
        children[1] = new QuadtreeNode(x_mid, x_max, y_min, y_mid);
        children[2] = new QuadtreeNode(x_min, x_mid, y_mid, y_max);
        children[3] = new QuadtreeNode(x_mid, x_max, y_mid, y_max);
    }

    void insert(Body* body) {
        if (isLeaf()) {
            if (this->body == nullptr) {
                this->body = body;
            } else {
                subdivide();
                insertIntoChildren(this->body);
                insertIntoChildren(body);
                this->body = nullptr;
            }
        } else {
            insertIntoChildren(body);
        }
        updateMassDistribution();
    }

    void insertIntoChildren(Body* body) {
        for (int i = 0; i < 4; ++i) {
            if (children[i]->contains(body)) {
                children[i]->insert(body);
                break;
            }
        }
    }

    bool contains(Body* body) const {
        return (body->x >= x_min && body->x < x_max && body->y >= y_min && body->y < y_max);
    }

    void updateMassDistribution() {
        mass = 0;
        center_of_mass_x = 0;
        center_of_mass_y = 0;
        if (isLeaf()) {
            if (body != nullptr) {
                mass = body->mass;
                center_of_mass_x = body->x;
                center_of_mass_y = body->y;
            }
        } else {
            for (int i = 0; i < 4; ++i) {
                if (children[i] != nullptr) {
                    mass += children[i]->mass;
                    center_of_mass_x += children[i]->mass * children[i]->center_of_mass_x;
                    center_of_mass_y += children[i]->mass * children[i]->center_of_mass_y;
                }
            }
            if (mass > 0) {
                center_of_mass_x /= mass;
                center_of_mass_y /= mass;
            }
        }
    }

    ~QuadtreeNode() {
        for (int i = 0; i < 4; ++i) {
            delete children[i];
        }
    }
};

void computeForce(Body* body, QuadtreeNode* node, double theta = 0.5) {
    if (node->isLeaf()) {
        if (node->body != nullptr && node->body != body) {
            double dx = node->body->x - body->x;
            double dy = node->body->y - body->y;
            double distance = std::sqrt(dx * dx + dy * dy);
            if (distance > 0) {
                double force = G * body->mass * node->body->mass / (distance * distance);
                body->ax += force * dx / distance;
                body->ay += force * dy / distance;
            }
        }
    } else {
        double dx = node->center_of_mass_x - body->x;
        double dy = node->center_of_mass_y - body->y;
        double distance = std::sqrt(dx * dx + dy * dy);
        double size = node->x_max - node->x_min;
        if (size / distance < theta) {
            double force = G * body->mass * node->mass / (distance * distance);
            body->ax += force * dx / distance;
            body->ay += force * dy / distance;
        } else {
            for (int i = 0; i < 4; ++i) {
                if (node->children[i] != nullptr) {
                    computeForce(body, node->children[i], theta);
                }
            }
        }
    }
}

int main() {
    // Set up the window
    sf::RenderWindow window(sf::VideoMode(800, 800), "Particle Simulation");

    // Create bodies
    std::vector<Body*> bodies = {
        new Body(400, 400, 1e14), new Body(450, 400, 1e14), new Body(400, 450, 1e14), new Body(450, 450, 1e14)
    };

    // Main loop
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        // Insert bodies into the quadtree
        QuadtreeNode* root = new QuadtreeNode(0, 800, 0, 800);
        for (auto body : bodies) {
            root->insert(body);
        }

        // Compute forces and update bodies
        for (auto body : bodies) {
            computeForce(body, root);
            body->update();
        }

        // Clear the window
        window.clear();

        // Draw bodies
        for (auto body : bodies) {
            body->draw(window);
        }

        // Display the contents of the window
        window.display();

        // Clean up the quadtree
        delete root;
    }

    // Clean up bodies
    for (auto body : bodies) {
        delete body;
    }

    return 0;
}
