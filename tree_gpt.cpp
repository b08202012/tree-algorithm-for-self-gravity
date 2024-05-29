#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <array>

const double G = 6.67430e-11; // Gravitational constant
const double TIME_STEP = 0.01; // Time step for simulation
const int NUM_STEPS = 1000; // Number of simulation steps

class Body {
public:
    double x, y, mass, vx, vy, ax, ay;

    Body(double x, double y, double mass, double vx, double vy, double ax, double ay)
        : x(x), y(y), mass(mass), vx(vx), vy(vy), ax(ax), ay(ay) {}

    void updatePosition(double dt) {
        vx += ax * dt;
        vy += ay * dt;
        x += vx * dt;
        y += vy * dt;
    }

    void resetAcceleration() {
        ax = 0;
        ay = 0;
    }
};

class QuadtreeNode {
public:
    double x_min, x_max, y_min, y_max;
    Body* bp;
    double mass, center_of_mass_x, center_of_mass_y;
    bool subdivided;
    std::array<QuadtreeNode, 4> children;

    QuadtreeNode(double x_min, double x_max, double y_min, double y_max)
        : x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max), body(nullptr), mass(0),
          center_of_mass_x(0), center_of_mass_y(0), subdivided(false) {}

    bool contains(Body* bp) const {
        return (bp->x >= x_min && bp->x < x_max && bp->y >= y_min && bp->y < y_max);
    }

    bool isLeaf() const {
        return !subdivided;
    }

    void subdivide() {
        double x_mid = (x_min + x_max) / 2;
        double y_mid = (y_min + y_max) / 2;
        children[0] = QuadtreeNode(x_min, x_mid, y_min, y_mid);
        children[1] = QuadtreeNode(x_mid, x_max, y_min, y_mid);
        children[2] = QuadtreeNode(x_min, x_mid, y_mid, y_max);
        children[3] = QuadtreeNode(x_mid, x_max, y_mid, y_max);
        subdivided = true;
    }

    void insert(Body* bp) {
        if (!contains(body)) {
            return;
        }

        if (isLeaf()) {
            if (this->bp == nullptr) {
                this->bp = body;
            } else {
                subdivide();
                insertIntoChildren(this->bp);
                insertIntoChildren(body);
                this->body = nullptr;
            }
        } else {
            insertIntoChildren(body);
        }
        updateMassDistribution();
    }

    void insertIntoChildren(Body* body) {
        for (auto& child : *children) {
            if (child.contains(body)) {
                child.insert(body);
                break;
            }
        }
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
            for (const auto& child : *children) {
                mass += child.mass;
                center_of_mass_x += child.mass * child.center_of_mass_x;
                center_of_mass_y += child.mass * child.center_of_mass_y;
            }
            if (mass > 0) {
                center_of_mass_x /= mass;
                center_of_mass_y /= mass;
            }
        }
    }

    ~QuadtreeNode() {
        delete children;
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
            for (auto& child : *node->children) {
                computeForce(body, &child, theta);
            }
        }
    }
}

void simulate(std::vector<Body*>& bodies, double timeStep, int steps, std::ofstream& outFile) {
    for (int step = 0; step < steps; ++step) {
        // Create the root of the quadtree
        double minCoord = -1000; // Adjust these bounds as necessary
        double maxCoord = 1000;
        QuadtreeNode root(minCoord, maxCoord, minCoord, maxCoord);

        // Insert all bodies into the quadtree
        for (auto body : bodies) {
            root.insert(body);
        }

        // Compute forces for each body
        for (auto body : bodies) {
            body->resetAcceleration();
            computeForce(body, &root);
        }

        // Update positions and velocities
        for (auto body : bodies) {
            body->updatePosition(timeStep);
        }

        // Write positions to file
        outFile << std::scientific << std::setprecision(16);
        for (auto body : bodies) {
            outFile << body->mass << " "
                    << body->x << " "
                    << body->y << " "
                    << 0.0 << " " // z position (always 0 in 2D)
                    << body->vx << " "
                    << body->vy << " "
                    << 0.0 << " " // z velocity (always 0 in 2D)
                    << 1.0 << " " // Particle type (constant)
                    << body->ax << " "
                    << body->ay << " "
                    << 0.0 << " " // z acceleration (always 0 in 2D)
                    << (step * timeStep) << std::endl;
        }
    }
}

std::vector<Body*> readParticlesFromFile(const std::string& filename) {
    std::ifstream inFile(filename);
    std::vector<Body*> bodies;
    std::string line;
    while (std::getline(inFile, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        double mass, x, y, z, vx, vy, vz, type, ax, ay, az, time;
        iss >> mass >> x >> y >> z >> vx >> vy >> vz >> type >> ax >> ay >> az >> time;
        // Create Body instance with 2D properties
        bodies.push_back(new Body(x, y, mass, vx, vy, ax, ay));
    }
    return bodies;
}

int main() {
    std::ofstream outFile("particles_2d_simulated.txt");

    // Read bodies from file
    std::vector<Body*> bodies = readParticlesFromFile("particles_2d.txt");

    // Simulate
    simulate(bodies, TIME_STEP, NUM_STEPS, outFile);

    // Cleanup
    for (auto body : bodies) {
        delete body;
    }

    outFile.close();
    return 0;
}
