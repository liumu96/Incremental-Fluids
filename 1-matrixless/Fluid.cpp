#include <algorithm>
#include <math.h>
#include <iostream>

#include "../loadpng/lodepng.h"

#include "FluidQuantity.h"
#include "FluidSolver.h"

using namespace std;

int main()
{

    const int sizeX = 128;
    const int sizeY = 128;

    const double density = 0.1;
    const double timestep = 0.005;

    unsigned char *image = new unsigned char[sizeX * sizeY * 4];

    // fluid solver
    FluidSolver *solver = new FluidSolver(sizeX, sizeY, density);

    double time = 0.0;
    int iterations = 0;

    while (time < 8.0)
    {
        // Use four substeps per iteration
        for (int i = 0; i < 4; i++)
        {
            solver->addInflow(0.45, 0.2, 0.1, 0.01, 1.0, 0.0, 3.0);
            solver->update(timestep);
            time += timestep;
            fflush(stdout);
        }

        // output image
        solver->toImage(image);

        char path[256];
        snprintf(path, sizeX, "Frame%05d.png", iterations++);
        lodepng_encode32_file(path, image, sizeX, sizeY);
    }

    return -1;
}