#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

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
        if (!contains(body)) {
            return;
        }

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

void simulate(std::vector<Body*>& bodies, double timeStep, int steps, std::ofstream& outFile) {
    for (int step = 0; step < steps; ++step) {
        // Calculate the bounds dynamically
        double minCoordX = std::numeric_limits<double>::max();
        double maxCoordX = std::numeric_limits<double>::min();
        double minCoordY = std::numeric_limits<double>::max();
        double maxCoordY = std::numeric_limits<double>::min();

        for (const auto& body : bodies) {
            if (body->x < minCoordX) minCoordX = body->x;
            if (body->x > maxCoordX) maxCoordX = body->x;
            if (body->y < minCoordY) minCoordY = body->y;
            if (body->y > maxCoordY) maxCoordY = body->y;
        }

        double margin = 0.1 * std::max(maxCoordX - minCoordX, maxCoordY - minCoordY);
        minCoordX -= margin;
        maxCoordX += margin;
        minCoordY -= margin;
        maxCoordY += margin;

        // Create the root of the quadtree
        QuadtreeNode* root = new QuadtreeNode(minCoordX, maxCoordX, minCoordY, maxCoordY);

        // Insert all bodies into the quadtree
        for (auto body : bodies) {
            root->insert(body);
        }

        // Compute forces for each body
        for (auto body : bodies) {
            body->resetAcceleration();
            computeForce(body, root);
        }

        // Update positions and velocities
        for (auto body : bodies) {
            body->updatePosition(timeStep);
        }

        // Write positions to file
        outFile << std::scientific << std::setprecision(16);
        for (auto body : bodies) {
            outFile << std::setw(24) << body->mass
                    << std::setw(24) << body->x
                    << std::setw(24) << body->y
                    << std::setw(24) << 0.0 // z position (always 0 in 2D)
                    << std::setw(24) << body->vx
                    << std::setw(24) << body->vy
                    << std::setw(24) << 0.0 // z velocity (always 0 in 2D)
                    << std::setw(24) << 1.0 // Particle type (constant)
                    << std::setw(24) << body->ax
                    << std::setw(24) << body->ay
                    << std::setw(24) << 0.0 // z acceleration (always 0 in 2D)
                    << std::setw(24) << (step * timeStep) << std::endl;
        }

        delete root;
    }
}

std::vector<Body*> readParticlesFromFile(const std::string& filename) {
    std::ifstream inFile(filename);
    if (!inFile) {
        throw std::runtime_error("Unable to open file: " + filename);
    }
    std::vector<Body*> bodies;
    std::string line;
    while (std::getline(inFile, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        double mass, x, y, z, vx, vy, vz, type, ax, ay, az, time;
        iss >> mass >> x >> y >> z >> vx >> vy >> vz >> type >> ax >> ay >> az >> time;
        bodies.push_back(new Body(x, y, mass, vx, vy, ax, ay));
    }
    return bodies;
}

int main() {
    try {
        std::ofstream outFile("output.txt");
        if (!outFile) {
            throw std::runtime_error("Unable to open output file");
        }

        // Read bodies from file
        std::vector<Body*> bodies = readParticlesFromFile("IC_16.txt");

        // Simulate
        simulate(bodies, TIME_STEP, NUM_STEPS, outFile);

        // Cleanup
        for (auto body : bodies) {
            delete body;
        }

        outFile.close();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
