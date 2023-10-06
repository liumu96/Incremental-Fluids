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
    FluidQuantity *_t; /* Temperature */
    FluidQuantity *_u;
    FluidQuantity *_v;

    /* Densities at staggered grid locations */
    double *_uDensity;
    double *_vDensity;

    /* Width and height */
    int _w;
    int _h;

    /* Grid cell size and fluid density */
    double _hx;
    double _densityAir;  /* Density of air */
    double _densitySoot; /* Density of soot */
    double _diffusion;   /* Diffusion rate of heat */

    /* Arrays for: */
    double *_r;      /* Right hand side of pressure solve */
    double *_p;      /* Pressure solution */
    double *_z;      /* Auxiliary vector */
    double *_s;      /* Search vector */
    double *_precon; /* Preconditioner */

    double *_aDiag;  /* Matrix diagonal */
    double *_aPlusX; /* Matrix off-diagonals */
    double *_aPlusY;

    /* Ambient temperature (here room temperature), in Kelvin */
    double _tAmb;
    /* Gravity */
    double _g;

    /* List of solid bodies to consider in the simulation */
    const vector<const SolidBody *> &_bodies;

    /* We now modify the right hand side to "blend" between solid and fluid
     * velocity based on the cell volume occupied by fluid.
     */
    void buildRhs()
    {
        double scale = 1.0 / _hx;
        const uint8_t *cell = _d->cell();
        const uint8_t *body = _d->body();

        for (int y = 0, idx = 0; y < _h; y++)
        {
            for (int x = 0; x < _w; x++, idx++)
            {
                if (cell[idx] == CELL_FLUID)
                {
                    _r[idx] = -scale *
                              (_u->volume(x + 1, y) * _u->at(x + 1, y) - _u->volume(x, y) * _u->at(x, y) +
                               _v->volume(x, y + 1) * _v->at(x, y + 1) - _v->volume(x, y) * _v->at(x, y));

                    double vol = _d->volume(x, y);

                    if (_bodies.empty())
                        continue;

                    if (x > 0)
                        _r[idx] -= (_u->volume(x, y) - vol) * _bodies[body[idx - 1]]->velocityX(x * _hx, (y + 0.5) * _hx);
                    if (y > 0)
                        _r[idx] -= (_v->volume(x, y) - vol) * _bodies[body[idx - _w]]->velocityY((x + 0.5) * _hx, y * _hx);
                    if (x < _w - 1)
                        _r[idx] += (_u->volume(x + 1, y) - vol) * _bodies[body[idx + 1]]->velocityX((x + 1.0) * _hx, (y + 0.5) * _hx);
                    if (y < _h - 1)
                        _r[idx] += (_v->volume(x, y + 1) - vol) * _bodies[body[idx + _w]]->velocityY((x + 0.5) * _hx, (y + 1.0) * _hx);
                }
                else
                    _r[idx] = 0.0;
            }
        }
    }

    /* Computes densities at the staggered grid locations as a function of
     * temperature and smoke concentration.
     */
    void computeDensities()
    {
        double alpha = (_densitySoot - _densityAir) / _densityAir;

        memset(_uDensity, 0, (_w + 1) * _h * sizeof(double));
        memset(_vDensity, 0, _w * (_h + 1) * sizeof(double));

        for (int y = 0; y < _h; y++)
        {
            for (int x = 0; x < _w; x++)
            {
                double density = _densityAir * _tAmb / _t->at(x, y) * (1.0 + alpha * _d->at(x, y));
                density = max(density, 0.05 * _densityAir); /* Clamp dangerously low densities */

                _uDensity[_u->idx(x, y)] += 0.5 * density;
                _vDensity[_v->idx(x, y)] += 0.5 * density;
                _uDensity[_u->idx(x + 1, y)] += 0.5 * density;
                _vDensity[_v->idx(x, y + 1)] += 0.5 * density;
            }
        }
    }

    /* Instead of constant density per cell, the entries must now be modified
     * to account for variable density at individual grid cells.
     */
    void buildPressureMatrix(double timestep)
    {
        double scale = timestep / (_hx * _hx);
        const uint8_t *cell = _d->cell();

        memset(_aDiag, 0, _w * _h * sizeof(double));
        memset(_aPlusX, 0, _w * _h * sizeof(double));
        memset(_aPlusY, 0, _w * _h * sizeof(double));

        for (int y = 0, idx = 0; y < _h; y++)
        {
            for (int x = 0; x < _w; x++, idx++)
            {
                if (cell[idx] != CELL_FLUID)
                    continue;

                if (x < _w - 1 && cell[idx + 1] == CELL_FLUID)
                {
                    double factor = scale * _u->volume(x + 1, y) / _uDensity[_u->idx(x + 1, y)];
                    _aDiag[idx] += factor;
                    _aDiag[idx + 1] += factor;
                    _aPlusX[idx] = -factor;
                }
                if (y < _h - 1 && cell[idx + _w] == CELL_FLUID)
                {
                    double factor = scale * _v->volume(x, y + 1) / _vDensity[_u->idx(x, y + 1)];
                    _aDiag[idx] += factor;
                    _aDiag[idx + _w] += factor;
                    _aPlusY[idx] = -factor;
                }
            }
        }
    }

    /* Constructs the matrix used to compute heat diffusion.
     * This uses the diffusion constant, which specifies how quickly heat
     * diffuses in the domain. Higher values imply faster diffusion and a
     * more difficult systems of equations (CG will converge slower)
     */
    void buildHeatDiffusionMatrix(double timestep)
    {
        /* Initialize diagonal to one. */
        for (int i = 0; i < _w * _h; i++)
            _aDiag[i] = 1.0;

        memset(_aPlusX, 0, _w * _h * sizeof(double));
        memset(_aPlusY, 0, _w * _h * sizeof(double));

        const uint8_t *cell = _d->cell();
        double scale = _diffusion * timestep * 1.0 / (_hx * _hx);

        for (int y = 0, idx = 0; y < _h; y++)
        {
            for (int x = 0; x < _w; x++, idx++)
            {
                if (cell[idx] != CELL_FLUID)
                    continue;

                if (x < _w - 1 && cell[idx + 1] == CELL_FLUID)
                {
                    _aDiag[idx] += scale;
                    _aDiag[idx + 1] += scale;
                    _aPlusX[idx] = -scale;
                }

                if (y < _h - 1 && cell[idx + _w] == CELL_FLUID)
                {
                    _aDiag[idx] += scale;
                    _aDiag[idx + _w] += scale;
                    _aPlusY[idx] = -scale;
                }
            }
        }
    }

    /* Builds the modified incomplete Cholesky preconditioner */
    void buildPreconditioner()
    {
        const double tau = 0.97;
        const double sigma = 0.25;
        const uint8_t *cell = _d->cell();

        for (int y = 0, idx = 0; y < _h; y++)
        {
            for (int x = 0; x < _w; x++, idx++)
            {
                if (cell[idx] != CELL_FLUID)
                    continue;

                double e = _aDiag[idx];

                if (x > 0 && cell[idx - 1] == CELL_FLUID)
                {
                    double px = _aPlusX[idx - 1] * _precon[idx - 1];
                    double py = _aPlusY[idx - 1] * _precon[idx - 1];
                    e = e - (px * px + tau * px * py);
                }
                if (y > 0 && cell[idx - _w] == CELL_FLUID)
                {
                    double px = _aPlusX[idx - _w] * _precon[idx - _w];
                    double py = _aPlusY[idx - _w] * _precon[idx - _w];
                    e = e - (py * py + tau * px * py);
                }

                if (e < sigma * _aDiag[idx])
                    e = _aDiag[idx];

                _precon[idx] = 1.0 / sqrt(e);
            }
        }
    }

    /* Apply preconditioner to vector `a' and store it in `dst' */
    void applyPreconditioner(double *dst, double *a)
    {
        const uint8_t *cell = _d->cell();

        for (int y = 0, idx = 0; y < _h; y++)
        {
            for (int x = 0; x < _w; x++, idx++)
            {
                if (cell[idx] != CELL_FLUID)
                    continue;

                double t = a[idx];

                if (x > 0 && cell[idx - 1] == CELL_FLUID)
                    t -= _aPlusX[idx - 1] * _precon[idx - 1] * dst[idx - 1];
                if (y > 0 && cell[idx - _w] == CELL_FLUID)
                    t -= _aPlusY[idx - _w] * _precon[idx - _w] * dst[idx - _w];

                dst[idx] = t * _precon[idx];
            }
        }

        for (int y = _h - 1, idx = _w * _h - 1; y >= 0; y--)
        {
            for (int x = _w - 1; x >= 0; x--, idx--)
            {
                if (cell[idx] != CELL_FLUID)
                    continue;

                double t = dst[idx];

                if (x < _w - 1 && cell[idx + 1] == CELL_FLUID)
                    t -= _aPlusX[idx] * _precon[idx] * dst[idx + 1];
                if (y < _h - 1 && cell[idx + _w] == CELL_FLUID)
                    t -= _aPlusY[idx] * _precon[idx] * dst[idx + _w];

                dst[idx] = t * _precon[idx];
            }
        }
    }

    /* Returns the dot product of vectors `a' and `b' */
    double dotProduct(double *a, double *b)
    {
        const uint8_t *cell = _d->cell();

        double result = 0.0;
        for (int i = 0; i < _w * _h; i++)
            if (cell[i] == CELL_FLUID)
                result += a[i] * b[i];
        return result;
    }

    /* Multiplies internal pressure matrix with vector `b' and stores the result in `dst' */
    void matrixVectorProduct(double *dst, double *b)
    {
        for (int y = 0, idx = 0; y < _h; y++)
        {
            for (int x = 0; x < _w; x++, idx++)
            {
                double t = _aDiag[idx] * b[idx];

                if (x > 0)
                    t += _aPlusX[idx - 1] * b[idx - 1];
                if (y > 0)
                    t += _aPlusY[idx - _w] * b[idx - _w];
                if (x < _w - 1)
                    t += _aPlusX[idx] * b[idx + 1];
                if (y < _h - 1)
                    t += _aPlusY[idx] * b[idx + _w];

                dst[idx] = t;
            }
        }
    }

    /* Computes `dst' = `a' + `b'*`s' */
    void scaledAdd(double *dst, double *a, double *b, double s)
    {
        const uint8_t *cell = _d->cell();

        for (int i = 0; i < _w * _h; i++)
            if (cell[i] == CELL_FLUID)
                dst[i] = a[i] + b[i] * s;
    }

    /* Returns maximum absolute value in vector `a' */
    double infinityNorm(double *a)
    {
        const uint8_t *cell = _d->cell();

        double maxA = 0.0;
        for (int i = 0; i < _w * _h; i++)
            if (cell[i] == CELL_FLUID)
                maxA = max(maxA, fabs(a[i]));
        return maxA;
    }

    /* Conjugate gradients solver */
    void project(int limit)
    {
        memset(_p, 0, _w * _h * sizeof(double)); /* Initial guess of zeroes */
        applyPreconditioner(_z, _r);
        memcpy(_s, _z, _w * _h * sizeof(double));

        double maxError = infinityNorm(_r);
        if (maxError < 1e-5)
            return;

        double sigma = dotProduct(_z, _r);

        for (int iter = 0; iter < limit; iter++)
        {
            matrixVectorProduct(_z, _s);
            double alpha = sigma / dotProduct(_z, _s);
            scaledAdd(_p, _p, _s, alpha);
            scaledAdd(_r, _r, _z, -alpha);

            maxError = infinityNorm(_r);
            if (maxError < 1e-5)
            {
                printf("Exiting solver after %d iterations, maximum error is %f\n", iter, maxError);
                return;
            }

            applyPreconditioner(_z, _r);

            double sigmaNew = dotProduct(_z, _r);
            scaledAdd(_s, _z, _s, sigmaNew / sigma);
            sigma = sigmaNew;
        }

        printf("Exceeded budget of %d iterations, maximum error was %f\n", limit, maxError);
    }

    /* Similar to the pressure matrix, we cannot assume constant density per
     * cell here either and must modify the equations accordingly.
     */
    void applyPressure(double timestep)
    {
        double scale = timestep / _hx;
        const uint8_t *cell = _d->cell();

        for (int y = 0, idx = 0; y < _h; y++)
        {
            for (int x = 0; x < _w; x++, idx++)
            {
                if (cell[idx] != CELL_FLUID)
                    continue;

                _u->at(x, y) -= scale * _p[idx] / _uDensity[_u->idx(x, y)];
                _v->at(x, y) -= scale * _p[idx] / _vDensity[_v->idx(x, y)];
                _u->at(x + 1, y) += scale * _p[idx] / _uDensity[_u->idx(x + 1, y)];
                _v->at(x, y + 1) += scale * _p[idx] / _vDensity[_v->idx(x, y + 1)];
            }
        }
    }

    /* Add body force due to density and heat difference */
    void addBuoyancy(double timestep)
    {
        double alpha = (_densitySoot - _densityAir) / _densityAir;

        for (int y = 0; y < _h; y++)
        {
            for (int x = 0; x < _w; x++)
            {
                /* Higher density (more soot) causes downward force, high
                 * temperatures causes upward force - and vice-versa.
                 */
                double buoyancy = timestep * _g * (alpha * _d->at(x, y) - (_t->at(x, y) - _tAmb) / _tAmb);

                _v->at(x, y) += buoyancy * 0.5;
                _v->at(x, y + 1) += buoyancy * 0.5;
            }
        }
    }

    /* Sets all velocity cells bordering solid cells to the solid velocity */
    void setBoundaryCondition()
    {
        const uint8_t *cell = _d->cell();
        const uint8_t *body = _d->body();

        for (int y = 0, idx = 0; y < _h; y++)
        {
            for (int x = 0; x < _w; x++, idx++)
            {
                if (cell[idx] == CELL_SOLID)
                {
                    const SolidBody &b = *_bodies[body[idx]];

                    _u->at(x, y) = b.velocityX(x * _hx, (y + 0.5) * _hx);
                    _v->at(x, y) = b.velocityY((x + 0.5) * _hx, y * _hx);
                    _u->at(x + 1, y) = b.velocityX((x + 1.0) * _hx, (y + 0.5) * _hx);
                    _v->at(x, y + 1) = b.velocityY((x + 0.5) * _hx, (y + 1.0) * _hx);
                }
            }
        }

        for (int y = 0; y < _h; y++)
            _u->at(0, y) = _u->at(_w, y) = 0.0;
        for (int x = 0; x < _w; x++)
            _v->at(x, 0) = _v->at(x, _h) = 0.0;
    }

public:
    FluidSolver(int w, int h, double rhoAir, double rhoSoot, double diffusion,
                const vector<const SolidBody *> &bodies) : _w(w), _h(h),
                                                           _densityAir(rhoAir), _densitySoot(rhoSoot), _diffusion(diffusion),
                                                           _bodies(bodies)
    {

        _tAmb = 294.0;
        _g = 9.81;

        _hx = 1.0 / min(w, h);

        _d = new FluidQuantity(_w, _h, 0.5, 0.5, _hx);
        _t = new FluidQuantity(_w, _h, 0.5, 0.5, _hx);
        _u = new FluidQuantity(_w + 1, _h, 0.0, 0.5, _hx);
        _v = new FluidQuantity(_w, _h + 1, 0.5, 0.0, _hx);

        for (int i = 0; i < _w * _h; i++)
            _t->src()[i] = _tAmb;

        _r = new double[_w * _h];
        _p = new double[_w * _h];
        _z = new double[_w * _h];
        _s = new double[_w * _h];
        _aDiag = new double[_w * _h];
        _aPlusX = new double[_w * _h];
        _aPlusY = new double[_w * _h];
        _precon = new double[_w * _h];

        _uDensity = new double[(_w + 1) * _h];
        _vDensity = new double[_w * (_h + 1)];
    }

    ~FluidSolver()
    {
        delete _d;
        delete _t;
        delete _u;
        delete _v;

        delete[] _r;
        delete[] _p;
        delete[] _z;
        delete[] _s;
        delete[] _aDiag;
        delete[] _aPlusX;
        delete[] _aPlusY;
        delete[] _precon;

        delete[] _uDensity;
        delete[] _vDensity;
    }

    void update(double timestep)
    {
        _d->fillSolidFields(_bodies);
        _t->fillSolidFields(_bodies);
        _u->fillSolidFields(_bodies);
        _v->fillSolidFields(_bodies);

        /* Right-hand side of heat equation is the current heat distribution */
        memcpy(_r, _t->src(), _w * _h * sizeof(double));
        buildHeatDiffusionMatrix(timestep);
        buildPreconditioner();
        project(2000);
        /* The solution of the heat equation is the heat distribution in the
         * next timestep.
         */
        memcpy(_t->src(), _p, _w * _h * sizeof(double));

        _t->extrapolate();

        addBuoyancy(timestep);
        setBoundaryCondition();

        buildRhs();
        computeDensities();
        buildPressureMatrix(timestep);
        buildPreconditioner();
        project(2000);
        applyPressure(timestep);

        _d->extrapolate();
        _u->extrapolate();
        _v->extrapolate();

        setBoundaryCondition();

        _d->advect(timestep, *_u, *_v, _bodies);
        _t->advect(timestep, *_u, *_v, _bodies);
        _u->advect(timestep, *_u, *_v, _bodies);
        _v->advect(timestep, *_u, *_v, _bodies);

        /* Make effect of advection visible, since it's not an in-place operation */
        _d->flip();
        _t->flip();
        _u->flip();
        _v->flip();
    }

    /* Set density and x/y velocity in given rectangle to d/u/v, respectively */
    void addInflow(double x, double y, double w, double h, double d, double t, double u, double v)
    {
        _d->addInflow(x, y, x + w, y + h, d);
        _t->addInflow(x, y, x + w, y + h, t);
        _u->addInflow(x, y, x + w, y + h, u);
        _v->addInflow(x, y, x + w, y + h, v);
    }

    double ambientT()
    {
        return _tAmb;
    }

    /* Convert fluid density to RGBA image */
    void toImage(unsigned char *rgba, bool renderHeat)
    {
        for (int y = 0; y < _h; y++)
        {
            for (int x = 0; x < _w; x++)
            {
                int idxl, idxr;
                if (renderHeat)
                {
                    idxl = 4 * (x + y * _w * 2);
                    idxr = 4 * (x + y * _w * 2 + _w);
                }
                else
                    idxr = 4 * (x + y * _w);

                double volume = _d->volume(x, y);

                double shade = (1.0 - _d->at(x, y)) * volume;
                shade = min(max(shade, 0.0), 1.0);
                rgba[idxr + 0] = (int)(shade * 255.0);
                rgba[idxr + 1] = (int)(shade * 255.0);
                rgba[idxr + 2] = (int)(shade * 255.0);
                rgba[idxr + 3] = 0xFF;

                if (renderHeat)
                {
                    double t = (_t->at(x, y) - _tAmb) / 700.0;
                    t = min(max(t, 0.0), 1.0);

                    double r = 1.0 + volume * (min(t * 4.0, 1.0) - 1.0);
                    double g = 1.0 + volume * (min(t * 2.0, 1.0) - 1.0);
                    double b = 1.0 + volume * (max(min(t * 4.0 - 3.0, 1.0), 0.0) - 1.0);

                    rgba[idxl + 0] = (int)(r * 255.0);
                    rgba[idxl + 1] = (int)(g * 255.0);
                    rgba[idxl + 2] = (int)(b * 255.0);
                    rgba[idxl + 3] = 0xFF;
                }
            }
        }
    }
};
