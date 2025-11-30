#include "LBM.h"
#include <cmath>

// D2Q9 Constants
// 0:Center, 1:E, 2:N, 3:W, 4:S, 5:NE, 6:NW, 7:SW, 8:SE
const int LBM::ex[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int LBM::ey[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
const double LBM::w[9] = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36};
// Opposite directions for bounce-back (0->0, 1->3, 2->4, etc.)
const int LBM::opposite[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

LBM::LBM(int width, int height, double tau) : width(width), height(height), tau(tau) {
    f.resize(width, std::vector<std::array<double, 9>>(height));
    f_new.resize(width, std::vector<std::array<double, 9>>(height));
    barrier.resize(width, std::vector<bool>(height, false));
    
    // 1. Initialize Barrier (Cylinder at x=width/4, y=height/2)
    int cx = width / 4;
    int cy = height / 2;
    int r = height / 9; // Radius

    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            if ((x - cx)*(x - cx) + (y - cy)*(y - cy) < r*r) {
                barrier[x][y] = true;
            }
        }
    }

    // 2. Initialize Flow (Wind blowing to the right)
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            for (int i = 0; i < 9; ++i) {
                // Initial wind speed u = 0.1
                f[x][y][i] = equilibrium(i, 1.0, 0.1, 0.0); 
            }
        }
    }
}

void LBM::step() {
    stream();
    collide();
    f = f_new;
}

void LBM::stream() {
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            // If this cell is a barrier, skip it (it has no fluid)
            if (barrier[x][y]) continue;

            for (int i = 0; i < 9; ++i) {
                int next_x = x + ex[i];
                int next_y = y + ey[i];

                // Periodic Boundary (Wrap around edges)
                if (next_x < 0) next_x = width - 1;
                if (next_x >= width) next_x = 0;
                if (next_y < 0) next_y = height - 1;
                if (next_y >= height) next_y = 0;

                // BOUNCE-BACK RULE:
                // If the target cell is a barrier, the particle bounces back to where it came from.
                if (barrier[next_x][next_y]) {
                    f_new[x][y][opposite[i]] = f[x][y][i];
                } else {
                    // Standard Stream
                    f_new[next_x][next_y][i] = f[x][y][i];
                }
            }
        }
    }
}

void LBM::collide() {
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            if (barrier[x][y]) continue; // Skip barriers

            double rho = 0.0;
            double ux = 0.0;
            double uy = 0.0;
            
            // Calculate Macroscopic moments
            for (int i = 0; i < 9; ++i) {
                rho += f_new[x][y][i];
                ux += f_new[x][y][i] * ex[i];
                uy += f_new[x][y][i] * ey[i];
            }
            
            if (rho > 0) {
                ux /= rho;
                uy /= rho;
            }

            // BGK Collision
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

// Getters for Visualization
double LBM::get_density(int x, int y) const {
    double rho = 0.0;
    for(int i=0; i<9; i++) rho += f[x][y][i];
    return rho;
}

bool LBM::is_barrier(int x, int y) const {
    return barrier[x][y];
}

double LBM::get_curl(int x, int y) const {
    // Simple finite difference curl: d(uy)/dx - d(ux)/dy
    // We need to compute local velocity first
    auto get_vel = [&](int cx, int cy) -> std::pair<double, double> {
        double r=0, u=0, v=0;
        for(int i=0; i<9; i++) {
            r += f[cx][cy][i];
            u += f[cx][cy][i] * ex[i];
            v += f[cx][cy][i] * ey[i];
        }
        if (r > 0) {
            return {u/r, v/r};
        } else {
            return {0.0, 0.0};
        }
    };

    int x1 = (x < width-1) ? x+1 : x;
    int x0 = (x > 0) ? x-1 : x;
    int y1 = (y < height-1) ? y+1 : y;
    int y0 = (y > 0) ? y-1 : y;

    double uy_x1 = get_vel(x1, y).second;
    double uy_x0 = get_vel(x0, y).second;
    double ux_y1 = get_vel(x, y1).first;
    double ux_y0 = get_vel(x, y0).first;

    return (uy_x1 - uy_x0) - (ux_y1 - ux_y0);
}
