#pragma once

#include <iostream>
#include <algorithm>
#include <math.h>

#include "FluidQuantity.h"

using namespace std;

class FluidSolver
{
private:
    /* Fluid quantities */
    FluidQuantity *_d;
    FluidQuantity *_u;
    FluidQuantity *_v;

    /* Width and height */
    int _w;
    int _h;

    /* Grid cell size and fluid density */
    double _hx;
    double _density;

    /* Arrays for: */
    double *_r; /* Right hand side of pressure solve */
    double *_p; /* Pressure solution */

    void buildRhs()
    {
        double scale = 1.0 / _hx;

        for (int y = 0, idx = 0; y < _h; y++)
            for (int x = 0; x < _w; x++, idx++)
            {
                _r[idx] = -scale * (_u->at(x + 1, y) - _u->at(x, y) + _v->at(x, y + 1) - _v->at(x, y));
            }
    }
    /* Performs the pressure solve using Gauss-Seidel.
     * The solver will run as long as it takes to get the relative error below
     * a threshold, but will never exceed `limit' iterations
     */
    void project(int limit, double timestep)
    {
        double scale = timestep / (_density * _hx * _hx);

        double maxDelta;
        for (int iter = 0; iter < limit; iter++)
        {
            maxDelta = 0.0;
            for (int y = 0, idx = 0; y < _h; y++)
                for (int x = 0; x < _w; x++, idx++)
                {
                    int idx = x + y * _w;

                    double diag = 0.0, offDiag = 0.0;

                    if (x > 0)
                    {
                        diag += scale;
                        offDiag -= scale * _p[idx - 1];
                    }
                    if (y > 0)
                    {
                        diag += scale;
                        offDiag -= scale * _p[idx - _w];
                    }
                    if (x < _w - 1)
                    {
                        diag += scale;
                        offDiag -= scale * _p[idx + 1];
                    }
                    if (y < _h - 1)
                    {
                        diag += scale;
                        offDiag -= scale * _p[idx + _w];
                    }

                    double newP = (_r[idx] - offDiag) / diag;

                    maxDelta = max(maxDelta, fabs(_p[idx] - newP));

                    _p[idx] = newP;
                }

            if (maxDelta < 1e-5)
            {
                printf("Exiting solver after %d iterations, maximum change is %f\n", iter, maxDelta);
                return;
            }
        }
        printf("Exceeded budget of %d iterations, maximum change was %f\n", limit, maxDelta);
    }

    /* Applies the computed pressure to the velocity field */
    void applyPressure(double timestep)
    {
        double scale = timestep / (_density * _hx);

        for (int y = 0, idx = 0; y < _h; y++)
            for (int x = 0; x < _w; x++, idx++)
            {
                _u->at(x, y) -= scale * _p[idx];
                _u->at(x + 1, y) += scale * _p[idx];
                _v->at(x, y) -= scale * _p[idx];
                _v->at(x, y + 1) += scale * _p[idx];
            }

        for (int y = 0; y < _h; y++)
            _u->at(0, y) = _u->at(_w, y) = 0.0;
        for (int x = 0; x < _w; x++)
            _v->at(x, 0) = _v->at(x, _h) = 0.0;
    }

public:
    FluidSolver(int w, int h, double density) : _w(w), _h(h), _density(density)
    {
        _hx = 1.0 / min(w, h);

        _d = new FluidQuantity(_w, _h, 0.5, 0.5, _hx);
        _u = new FluidQuantity(_w + 1, _h, 0.0, 0.5, _hx);
        _v = new FluidQuantity(_w, _h + 1, 0.5, 0.0, _hx);

        _r = new double[_w * _h];
        _p = new double[_w * _h];

        memset(_p, 0, _w * _h * sizeof(double));
    }

    ~FluidSolver()
    {
        delete _d;
        delete _u;
        delete _v;

        delete[] _r;
        delete[] _p;
    }

    void update(double timestep)
    {
        buildRhs();
        project(1200, timestep);
        applyPressure(timestep);

        _d->advect(timestep, *_u, *_v);
        _u->advect(timestep, *_u, *_v);
        _v->advect(timestep, *_u, *_v);

        /* Make effect of advection visible, since it's not an in-place operation */
        _d->flip();
        _u->flip();
        _v->flip();
    }

    /* Set density and x/y velocity in given rectangle to d/u/v, respectively */
    void addInflow(double x, double y, double w, double h, double d, double u, double v)
    {
        _d->addInflow(x, y, x + w, y + h, d);
        _u->addInflow(x, y, x + w, y + h, u);
        _v->addInflow(x, y, x + w, y + h, v);
    }

    /* Returns the maximum allowed timestep. Note that the actual timestep
     * taken should usually be much below this to ensure accurate
     * simulation - just never above.
     */
    double maxTimestep()
    {
        double maxVelocity = 0.0;
        for (int y = 0; y < _h; y++)
        {
            for (int x = 0; x < _w; x++)
            {
                /* Average velocity at grid cell center */
                double u = _u->lerp(x + 0.5, y + 0.5);
                double v = _v->lerp(x + 0.5, y + 0.5);

                double velocity = sqrt(u * u + v * v);
                maxVelocity = max(maxVelocity, velocity);
            }
        }

        /* Fluid should not flow more than two grid cells per iteration */
        double maxTimestep = 2.0 * _hx / maxVelocity;

        /* Clamp to sensible maximum value in case of very small velocities */
        return min(maxTimestep, 1.0);
    }

    /* Convert fluid density to RGBA image */
    void toImage(unsigned char *rgba)
    {
        for (int i = 0; i < _w * _h; i++)
        {
            int shade = (int)((1.0 - _d->src()[i]) * 255.0);
            shade = max(min(shade, 255), 0);

            rgba[i * 4 + 0] = shade;
            rgba[i * 4 + 1] = shade;
            rgba[i * 4 + 2] = shade;
            rgba[i * 4 + 3] = 0xFF;
        }
    }
};
