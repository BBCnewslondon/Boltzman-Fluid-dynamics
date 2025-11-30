#ifndef LBM_H
#define LBM_H

#include <vector>
#include <array>

class LBM {
public:
    LBM(int width, int height, double tau);
    void step();
    
    // Getters for visualization
    double get_density(int x, int y) const;
    double get_curl(int x, int y) const; // To see the vortices
    bool is_barrier(int x, int y) const;

    // Simulation Constants
    int width, height; // Made public for easier looping in main

private:
    double tau;

    static constexpr int Q = 9;
    static const int ex[Q];
    static const int ey[Q];
    static const double w[Q];
    static const int opposite[Q]; // NEW: For bouncing off walls

    // State
    std::vector<std::vector<std::array<double, 9>>> f;
    std::vector<std::vector<std::array<double, 9>>> f_new;
    std::vector<std::vector<bool>> barrier; // NEW: The obstacle grid
    
    void stream();
    void collide();
    double equilibrium(int i, double rho, double ux, double uy);
};

#endif
