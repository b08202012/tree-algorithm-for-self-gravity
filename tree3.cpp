#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

const double G = 1.0; // Gravitational constant
const double TIME_STEP = 0.1; // Time step for simulation
const int NUM_STEPS = 10000; // Number of simulation steps

class Body {
public:
    double  mass, x, y, z, vx, vy, vz, ax, ay, az;

    Body(double mass, double x, double y, double z, double vx, double vy, double vz, double ax, double ay, double az)
        : mass(mass), x(x), y(y), z(z), vx(vx), vy(vy), vz(vz), ax(ax), ay(ay), az(az) {}

    void updatePosition(double dt) { //DKD
        x += vx * 0.5 * dt;
        y += vy * 0.5 * dt;
        z += vz * 0.5 * dt;
        vx += ax * dt;
        vy += ay * dt;
        vz += az * dt;
        x += vx * 0.5 * dt;
        y += vy * 0.5 * dt;
        z += vz * 0.5 * dt;
    }

    void resetAcceleration() {
        ax = 0;
        ay = 0;
        az = 0;
    }
};

class OctreeNode {
public:
    double x_min, x_max, y_min, y_max,z_min, z_max;
    Body* body;
    double mass, center_of_mass_x, center_of_mass_y, center_of_mass_z;
    OctreeNode* children[8];

    OctreeNode(double x_min, double x_max, double y_min, double y_max,  double z_min, double z_max)
        : x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max), z_min(z_min), z_max(z_max), body(nullptr), mass(0),
          center_of_mass_x(0), center_of_mass_y(0), center_of_mass_z(0) {
        for (int i = 0; i < 8; ++i) {
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
                splitnode();
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
        return (   body->x >= x_min 
                && body->x < x_max 
                && body->y >= y_min 
                && body->y < y_max 
                && body->z >= z_min 
                && body->z < z_max);
    }

    bool isLeaf() const {
        for (int i = 0; i < 8; ++i) {
            if (children[i] != nullptr) {
                return false;
            }
        }
        return true;
    }

    void splitnode() {
        double x_mid = (x_min + x_max) / 2;
        double y_mid = (y_min + y_max) / 2;
        double z_mid = (z_min + z_max) / 2;
        children[0] = new OctreeNode(x_min, x_mid, y_min, y_mid, z_min, z_mid);
        children[1] = new OctreeNode(x_mid, x_max, y_min, y_mid, z_min, z_mid);
        children[2] = new OctreeNode(x_min, x_mid, y_mid, y_max, z_min, z_mid);
        children[3] = new OctreeNode(x_mid, x_max, y_mid, y_max, z_min, z_mid);
        children[4] = new OctreeNode(x_min, x_mid, y_min, y_mid, z_mid, z_max);
        children[5] = new OctreeNode(x_mid, x_max, y_min, y_mid, z_mid, z_max);
        children[6] = new OctreeNode(x_min, x_mid, y_mid, y_max, z_mid, z_max);
        children[7] = new OctreeNode(x_mid, x_max, y_mid, y_max, z_mid, z_max);
    }

    void insertIntoChildren(Body* body) {
        for (int i = 0; i < 8; ++i) {
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
        center_of_mass_z = 0;
        if (isLeaf()) {
            if (body != nullptr) {
                mass = body->mass;
                center_of_mass_x = body->x;
                center_of_mass_y = body->y;
                center_of_mass_z = body->z;
            }
        } else {
            for (int i = 0; i < 8; ++i) {
                if (children[i] != nullptr) {
                    mass += children[i]->mass;
                    center_of_mass_x += children[i]->mass * children[i]->center_of_mass_x;
                    center_of_mass_y += children[i]->mass * children[i]->center_of_mass_y;
                    center_of_mass_z += children[i]->mass * children[i]->center_of_mass_z;
                }
            }
            if (mass > 0) {
                center_of_mass_x /= mass;
                center_of_mass_y /= mass;
                center_of_mass_z /= mass;
            }
        }
    }

    ~OctreeNode() {
        for (int i = 0; i < 8; ++i) {
            delete children[i];
        }
    }
};

void computeForce(Body* body, OctreeNode* node, double theta = 0.5) {
    if (node->isLeaf()) {
        if (node->body != nullptr && node->body != body) {
            double dx = node->body->x - body->x;
            double dy = node->body->y - body->y;
            double dz = node->body->z - body->z;
            double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (distance > 0) {
                double accel = G * node->body->mass / (distance * distance);
                body->ax += accel * dx / distance;
                body->ay += accel * dy / distance;
                body->ay += accel * dz / distance;
            }
        }
    } else {
        double dx = node->center_of_mass_x - body->x;
        double dy = node->center_of_mass_y - body->y;
        double dz = node->center_of_mass_z - body->z;
        double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
        double size = node->x_max - node->x_min;
        if (size / distance < theta) {
            double accel = G * node->mass / (distance * distance);
            body->ax += accel * dx / distance;
            body->ay += accel * dy / distance;
            body->ay += accel * dz / distance;
        } else {
            for (int i = 0; i < 8; ++i) {
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
        double minCoordX = -10;
        double maxCoordX = 10;
        double minCoordY = -10;
        double maxCoordY = 10;
        double minCoordZ = -10;
        double maxCoordZ = 10;

        // Create the root of the quadtree
        OctreeNode* root = new OctreeNode(minCoordX, maxCoordX, minCoordY, maxCoordY, minCoordZ, maxCoordZ);

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

        if(std::remainder(step,1000) == 0){
        outFile << std::scientific << std::setprecision(8);
        for (auto body : bodies) {
            outFile << body->mass << " "
                    << body->x << " "
                    << body->y << " "
                    << body->z << " "
                    << body->vx << " "
                    << body->vy << " "
                    << body->vz << " "
                    << body->ax << " "
                    << body->ay << " "
                    << body->az << " "
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
        std::ofstream outFile("output100.txt");
        if (!outFile) {
            throw std::runtime_error("Unable to open output file");
        }

        // Read bodies from file
        std::vector<Body*> bodies = readFile("IC100.txt");

        // Simulate
        simulate(bodies, TIME_STEP, NUM_STEPS, outFile);

        // Cleanup
        for (auto body : bodies) {
            delete body;
        }
        outFile.close();
    return 0;
}
