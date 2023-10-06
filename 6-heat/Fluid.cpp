#include <algorithm>
#include <math.h>
#include <iostream>

#include "../loadpng/lodepng.h"

#include "FluidQuantity.h"
#include "FluidSolver.h"
#include "SolidBody.h"

using namespace std;

int main()
{

    const int sizeX = 128;
    const int sizeY = 128;

    const double densityAir = 0.1;
    const double densitySoot = 0.1; /* You can make this larger to get heavy smoke */
    const double diffusion = 0.01;
    const double timestep = 0.005;

    const bool renderHeat = true; /* Set this to false to disable heat rendering */

    unsigned char *image = new unsigned char[sizeX * 2 * sizeY * 4];

    vector<SolidBody *> bodies;
    bodies.push_back(new SolidBox(0.3, 0.6, 0.1, 0.5, -M_PI * 0.05, 0.0, 0.0, 0.0));

    /* Unfortunately we need a second vector here - the first one allows us to
     * modify bodies with the update method, this one is a const vector the
     * fluid solver is allowed to work with. Unfortunately C++ does not allow
     * us to cast a non-const to a const vector, so we need a separate one.
     */
    vector<const SolidBody *> cBodies;
    for (unsigned i = 0; i < bodies.size(); i++)
        cBodies.push_back(bodies[i]);

    // fluid solver
    FluidSolver *solver = new FluidSolver(sizeX, sizeY, densityAir, densitySoot, diffusion, cBodies);

    double time = 0.0;
    int iterations = 0;

    while (time < 8.0)
    {
        // Use four substeps per iteration
        for (int i = 0; i < 10; i++)
        {
            solver->addInflow(0.35, 0.9, 0.1, 0.05, 1.0, solver->ambientT() + 300.0, 0.0, 0.0);
            solver->update(timestep);
            time += timestep;
            fflush(stdout);
        }

        // output image
        solver->toImage(image, renderHeat);

        char path[256];
        snprintf(path, sizeX, "Frame%05d.png", iterations++);
        lodepng_encode32_file(path, image, (renderHeat ? sizeX * 2 : sizeX), sizeY);

        for (unsigned i = 0; i < bodies.size(); i++)
            bodies[i]->update(timestep);
    }

    return -1;
}