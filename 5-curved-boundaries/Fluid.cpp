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

    const double density = 0.1;
    const double timestep = 0.005;

    unsigned char *image = new unsigned char[sizeX * sizeY * 4];

    vector<SolidBody *> bodies;
    bodies.push_back(new SolidBox(0.5, 0.6, 0.7, 0.1, M_PI * 0.25, 0.0, 0.0, 0.0));

    /* Unfortunately we need a second vector here - the first one allows us to
     * modify bodies with the update method, this one is a const vector the
     * fluid solver is allowed to work with. Unfortunately C++ does not allow
     * us to cast a non-const to a const vector, so we need a separate one.
     */
    vector<const SolidBody *> cBodies;
    for (unsigned i = 0; i < bodies.size(); i++)
        cBodies.push_back(bodies[i]);

    // fluid solver
    FluidSolver *solver = new FluidSolver(sizeX, sizeY, density, cBodies);

    double time = 0.0;
    int iterations = 0;

    while (time < 8.0)
    {
        // Use four substeps per iteration
        for (int i = 0; i < 4; i++)
        {
            solver->addInflow(0.45, 0.2, 0.15, 0.03, 1.0, 0.0, 3.0);
            solver->update(timestep);
            time += timestep;
            fflush(stdout);
        }

        // output image
        solver->toImage(image);

        char path[256];
        snprintf(path, sizeX, "Frame%05d.png", iterations++);
        lodepng_encode32_file(path, image, sizeX, sizeY);

        for (unsigned i = 0; i < bodies.size(); i++)
            bodies[i]->update(timestep);
    }

    return -1;
}