#include <iostream>
#include <fstream>
#include <string>
#include "LBM.h"

void write_ppm(const std::string& filename, LBM& lbm) {
    std::ofstream out(filename);
    out << "P3\n" << lbm.width << " " << lbm.height << "\n255\n";

    for (int y = lbm.height - 1; y >= 0; --y) {
        for (int x = 0; x < lbm.width; ++x) {
            if (lbm.is_barrier(x, y)) {
                out << "128 128 128\n"; // Gray Barrier
                continue;
            }

            // Visualize Curl (Vorticity)
            // Red = Clockwise Spin, Blue = Counter-Clockwise
            double curl = lbm.get_curl(x, y);
            double val = (curl * 5.0) + 0.5; // Scale and bias to 0.5
            
            // Clamp
            if (val < 0) val = 0;
            if (val > 1) val = 1;

            int r = static_cast<int>(255.0 * val);
            int b = static_cast<int>(255.0 * (1.0 - val));
            int g = static_cast<int>(255.0 * (0.5 * val + 0.2)); // Slight whiteish tint

            out << r << " " << g << " " << b << "\n";
        }
    }
    out.close();
}

int main() {
    const int W = 400;
    const int H = 100;
    const double TAU = 0.6; // Low viscosity for turbulence

    std::cout << "Initializing Wind Tunnel (D2Q9)...\n";
    LBM lbm(W, H, TAU);

    // Warm up
    std::cout << "Warming up simulation...\n";
    for(int t=0; t<1000; t++) lbm.step();

    // Render loop
    std::cout << "Starting Render...\n";
    for (int t = 0; t < 400; ++t) {
        // Run 20 physics steps per frame for speed
        for(int s=0; s<20; s++) lbm.step(); 

        std::string name = "output/fluid_" + std::to_string(t) + ".ppm";
        write_ppm(name, lbm);
        
        if (t % 10 == 0) std::cout << "Frame " << t << " / 400\n";
    }
    
    return 0;
}
