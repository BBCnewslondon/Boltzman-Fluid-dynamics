#include "LBM.h"
#include <iostream>
#include <cmath>

// D2Q9 Constants
// 0: Center, 1: E, 2: N, 3: W, 4: S, 5: NE, 6: NW, 7: SW, 8: SE
const int LBM::ex[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int LBM::ey[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
const double LBM::w[9] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

LBM::LBM(int width, int height, double tau) : width(width), height(height), tau(tau) {
    f.resize(width, std::vector<std::array<double, 9>>(height));
    f_new.resize(width, std::vector<std::array<double, 9>>(height));
    
    // Initialize with some density and zero velocity
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            double rho = 1.0;
            // Add some perturbation or initial velocity if needed
            if (x == width / 2 && y == height / 2) rho = 1.2; // Initial pulse

            for (int i = 0; i < 9; ++i) {
                f[x][y][i] = equilibrium(i, rho, 0.0, 0.0);
            }
        }
    }
}

void LBM::step() {
    stream();
    collide();
    // Swap f and f_new
    f = f_new;
}

void LBM::stream() {
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            for (int i = 0; i < 9; ++i) {
                // Periodic boundary conditions
                int next_x = (x + ex[i] + width) % width;
                int next_y = (y + ey[i] + height) % height;
                f_new[next_x][next_y][i] = f[x][y][i];
            }
        }
    }
}

void LBM::collide() {
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            double rho = 0.0;
            double ux = 0.0;
            double uy = 0.0;
            
            // Calculate macroscopic quantities
            for (int i = 0; i < 9; ++i) {
                rho += f_new[x][y][i]; // Note: using f_new because we streamed into f_new
                ux += f_new[x][y][i] * ex[i];
                uy += f_new[x][y][i] * ey[i];
            }
            
            if (rho > 0) {
                ux /= rho;
                uy /= rho;
            }

            // BGK collision
            for (int i = 0; i < 9; ++i) {
                double feq = equilibrium(i, rho, ux, uy);
                f_new[x][y][i] = f_new[x][y][i] - (f_new[x][y][i] - feq) / tau;
            }
        }
    }
}

double LBM::equilibrium(int i, double rho, double ux, double uy) {
    double cu = 3.0 * (ex[i] * ux + ey[i] * uy);
    double u2 = ux * ux + uy * uy;
    return rho * w[i] * (1.0 + cu + 0.5 * cu * cu - 1.5 * u2);
}
