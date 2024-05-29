#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

const double G = 6.67430e-11; // Gravitational constant
const double TIME_STEP = 0.1; // Time step for simulation
const int NUM_STEPS = 100; // Number of simulation steps

class Body {
public:
    double mass, x, y, vx, vy, ax, ay;

    Body(double mass, double x, double y, double vx, double vy, double ax, double ay)
        : mass(mass), x(x), y(y), vx(vx), vy(vy), ax(ax), ay(ay) {}

    void updatePosition(double dt) { //DKD
        x += vx * 0.5 * dt;
        y += vy * 0.5 * dt;
        vx += ax * dt;
        vy += ay * dt;
        x += vx * 0.5 * dt;
        y += vy * 0.5 * dt;
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

    void insert(Body* body) {
        if (!inside_node(body)) {
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

    bool inside_node(Body* body) const {
        return (body->x >= x_min && body->x < x_max && body->y >= y_min && body->y < y_max);
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

    void insertIntoChildren(Body* body) {
        for (int i = 0; i < 4; ++i) {
            if (children[i]->inside_node(body)) {
                children[i]->insert(body);
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
            double accel = G * node->mass / (distance * distance);
            body->ax += accel * dx / distance;
            body->ay += accel * dy / distance;
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
        double minCoordX = -1e100;
        double maxCoordX = 1e100;
        double minCoordY = -1e100;
        double maxCoordY = 1e100;

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

        if(std::remainder(step,10) == 1){
        outFile << std::scientific << std::setprecision(8);
        for (auto body : bodies) {
            outFile << body->mass << "  "
                    << body->x << "  "
                    << body->y << "  "
                    << 0.0 << "  " // z position (always 0 in 2D)
                    << body->vx << "  "
                    << body->vy << "  "
                    << 0.0 << "  " // z velocity (always 0 in 2D)
                    << 1.0 << "  " // Particle type (constant)
                    << body->ax << "  "
                    << body->ay << "  "
                    << 0.0 << "  " // z acceleration (always 0 in 2D)
                    << (step * timeStep) << std::endl;
                    }
        }

        delete root;
    }
}

// Function to read the file and store it in a 2D array
std::vector<Body*> readFile(const std::string& filename) {
    std::ifstream inFile(filename);
    if (!inFile) {
        throw std::runtime_error("Unable to open file: " + filename);
    }
    std::vector<Body*> bodies;
    std::string line;
    while (std::getline(inFile, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        double mass, x, y, z, vx, vy, vz, type, ax, ay, az, time;
        iss >> mass >> x >> y >> z >> vx >> vy >> vz >> type >> ax >> ay >> az >> time;
        bodies.push_back(new Body(mass, x, y, vx, vy, ax, ay));
    }
    return bodies;
}

int main() {
        std::ofstream outFile("output.txt");
        if (!outFile) {
            throw std::runtime_error("Unable to open output file");
        }

        // Read bodies from file
        std::vector<Body*> bodies = readFile("IC16.txt");

        // Simulate
        simulate(bodies, TIME_STEP, NUM_STEPS, outFile);

        // Cleanup
        for (auto body : bodies) {
            delete body;
        }
        outFile.close();
    return 0;
}
