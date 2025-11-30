#ifndef LBM_H
#define LBM_H

#include <vector>
#include <array>

class LBM {
public:
    LBM(int width, int height, double tau);
    void step();

private:
    int width, height;
    double tau;

    static constexpr int Q = 9;
    static const int ex[Q];
    static const int ey[Q];
    static const double w[Q];

    // D2Q9 has 9 directions
    // We need two grids: one for current step, one for next step (streaming)
    // f[x][y][i]
    std::vector<std::vector<std::array<double, 9>>> f;
    std::vector<std::vector<std::array<double, 9>>> f_new;
    
    void stream();
    void collide();
    double equilibrium(int i, double rho, double ux, double uy);
};

#endif
