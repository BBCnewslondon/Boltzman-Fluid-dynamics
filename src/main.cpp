#include <iostream>
#include "LBM.h"

int main() {
    std::cout << "Starting LBM Simulation..." << std::endl;
    LBM lbm(100, 50, 0.6); // Example size and tau
    
    // Simulation loop
    for (int t = 0; t < 1000; ++t) {
        lbm.step();
        if (t % 100 == 0) {
            std::cout << "Step " << t << " completed." << std::endl;
        }
    }
    
    return 0;
}
