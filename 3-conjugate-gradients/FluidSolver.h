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
    double *_r;      /* Right hand side of pressure solve */
    double *_p;      /* Pressure solution */
    double *_z;      /* Auxiliary vector */
    double *_s;      /* Search vector */
    double *_precon; /* Preconditioner */

    double *_aDiag;  /* Matrix diagonal */
    double *_aPlusX; /* Matrix off-diagonals */
    double *_aPlusY;

    void buildRhs()
    {
        double scale = 1.0 / _hx;

        for (int y = 0, idx = 0; y < _h; y++)
            for (int x = 0; x < _w; x++, idx++)
            {
                _r[idx] = -scale * (_u->at(x + 1, y) - _u->at(x, y) + _v->at(x, y + 1) - _v->at(x, y));
            }
    }

    /* Builds the pressure matrix. Since the matrix is very sparse and
     * symmetric, it allows for memory friendly storage.
     */
    void buildPressureMatrix(double timestep)
    {
        double scale = timestep / (_density * _hx * _hx);

        memset(_aDiag, 0, _w * _h * sizeof(double));

        for (int y = 0, idx = 0; y < _h; y++)
        {
            for (int x = 0; x < _w; x++, idx++)
            {
                if (x < _w - 1)
                {
                    _aDiag[idx] += scale;
                    _aDiag[idx + 1] += scale;
                    _aPlusX[idx] = -scale;
                }
                else
                    _aPlusX[idx] = 0.0;

                if (y < _h - 1)
                {
                    _aDiag[idx] += scale;
                    _aDiag[idx + _w] += scale;
                    _aPlusY[idx] = -scale;
                }
                else
                    _aPlusY[idx] = 0.0;
            }
        }
    }

    /* Builds the modified incomplete Cholesky preconditioner */
    void buildPreconditioner()
    {
        const double tau = 0.97;
        const double sigma = 0.25;

        for (int y = 0, idx = 0; y < _h; y++)
        {
            for (int x = 0; x < _w; x++, idx++)
            {
                double e = _aDiag[idx];

                if (x > 0)
                {
                    double px = _aPlusX[idx - 1] * _precon[idx - 1];
                    double py = _aPlusY[idx - 1] * _precon[idx - 1];
                    e = e - (px * px + tau * px * py);
                }
                if (y > 0)
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
        for (int y = 0, idx = 0; y < _h; y++)
        {
            for (int x = 0; x < _w; x++, idx++)
            {
                double t = a[idx];

                if (x > 0)
                    t -= _aPlusX[idx - 1] * _precon[idx - 1] * dst[idx - 1];
                if (y > 0)
                    t -= _aPlusY[idx - _w] * _precon[idx - _w] * dst[idx - _w];

                dst[idx] = t * _precon[idx];
            }
        }

        for (int y = _h - 1, idx = _w * _h - 1; y >= 0; y--)
        {
            for (int x = _w - 1; x >= 0; x--, idx--)
            {
                idx = x + y * _w;

                double t = dst[idx];

                if (x < _w - 1)
                    t -= _aPlusX[idx] * _precon[idx] * dst[idx + 1];
                if (y < _h - 1)
                    t -= _aPlusY[idx] * _precon[idx] * dst[idx + _w];

                dst[idx] = t * _precon[idx];
            }
        }
    }

    /* Returns the dot product of vectors `a' and `b' */
    double dotProduct(double *a, double *b)
    {
        double result = 0.0;
        for (int i = 0; i < _w * _h; i++)
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
        for (int i = 0; i < _w * _h; i++)
            dst[i] = a[i] + b[i] * s;
    }

    /* Returns maximum absolute value in vector `a' */
    double infinityNorm(double *a)
    {
        double maxA = 0.0;
        for (int i = 0; i < _w * _h; i++)
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
        _z = new double[_w * _h];
        _s = new double[_w * _h];
        _aDiag = new double[_w * _h];
        _aPlusX = new double[_w * _h];
        _aPlusY = new double[_w * _h];
        _precon = new double[_w * _h];
    }

    ~FluidSolver()
    {
        delete _d;
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
    }

    void update(double timestep)
    {
        buildRhs();
        buildPressureMatrix(timestep);
        buildPreconditioner();
        project(600);
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
