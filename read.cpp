#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

// Function to read the file and store it in a 2D array
vector<vector<double>> readFileToArray(const string& filename) {
    std::ifstream file(filename);
    std::string line;


    vector<vector<double>> data;
    while (std::getline(file, line)) {

        if (line[0] == '#') {
            continue;
        }

        std::stringstream ss(line);
        vector<double> row;
        double value;
        while (ss >> value) {
            row.push_back(value);
        }
        data.push_back(row);
    }

    file.close();
    return data;
}

int main() {
    string filename = "IC16.txt"; // Replace with your file name
    vector<vector<double>> data = readFileToArray(filename);

    // Print the 2D array
 double var = data[1][3];
 cout << var << endl;
    

    return 0;
}
