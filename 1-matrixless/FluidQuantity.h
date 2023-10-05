#pragma once

#include <iostream>
#include <algorithm>
#include <math.h>

using namespace std;

class FluidQuantity
{
private:
    // Memory buffers for fluid quantity
    double *_src;
    double *_dst;

    // width and height
    int _w;
    int _h;
    /* X and Y offset from top left grid cell.
     * This is (0.5,0.5) for centered quantities such as density,
     * and (0.0, 0.5) or (0.5, 0.0) for jittered quantities like the velocity.
     */
    double _ox;
    double _oy;
    // Grid cell size
    double _hx;

    /* Linear intERPolate between a and b for x ranging from 0 to 1 */
    double lerp(double a, double b, double x) const
    {
        return a * (1.0 - x) + b * x;
    }

    /* Simple forward Euler method for velocity integration in time */
    void euler(double &x, double &y, double timestep, const FluidQuantity &u, const FluidQuantity &v) const
    {
        double uVel = u.lerp(x, y) / _hx;
        double vVel = v.lerp(x, y) / _hx;

        x -= uVel * timestep;
        y -= vVel * timestep;
    }

public:
    FluidQuantity(int w, int h, double ox, double oy, double hx) : _w(w), _h(h), _ox(ox), _oy(oy), _hx(hx)
    {
        _src = new double[_w * _h];
        _dst = new double[_w * _h];

        memset(_src, 0, _w * _h * sizeof(double));
    }
    ~FluidQuantity()
    {
        delete[] _src;
        delete[] _dst;
    }

    void flip()
    {
        swap(_src, _dst);
    }

    const double *src() const
    {
        return _src;
    }

    double at(int x, int y) const
    {
        return _src[x + y * _w];
    }

    double &at(int x, int y)
    {
        return _src[x + y * _w];
    }

    /* Linear intERPolate on grid at coordinates (x, y).
     * Coordinates will be clamped to lie in simulation domain
     */
    double lerp(double x, double y) const
    {
        x = min(max(x - _ox, 0.0), _w - 1.001);
        y = min(max(y - _oy, 0.0), _h - 1.001);
        int ix = (int)x;
        int iy = (int)y;
        x -= ix;
        y -= iy;

        double x00 = at(ix + 0, iy + 0), x10 = at(ix + 1, iy + 0);
        double x01 = at(ix + 0, iy + 1), x11 = at(ix + 1, iy + 1);

        return lerp(lerp(x00, x10, x), lerp(x01, x11, x), y);
    }

    /* Advect grid in velocity field u, v with given timestep */
    void advect(double timestep, const FluidQuantity &u, const FluidQuantity &v)
    {
        for (int iy = 0, idx = 0; iy < _h; iy++)
            for (int ix = 0; ix < _w; ix++, idx++)
            {
                double x = ix + _ox;
                double y = iy + _oy;

                /* First component: Integrate in time */
                euler(x, y, timestep, u, v);

                /* Second component: Interpolate from grid */
                _dst[idx] = lerp(x, y);
            }
    }

    /* Sets fluid quantity inside the given rect to value `v' */
    void addInflow(double x0, double y0, double x1, double y1, double v)
    {
        int ix0 = (int)(x0 / _hx - _ox);
        int iy0 = (int)(y0 / _hx - _oy);
        int ix1 = (int)(x1 / _hx - _ox);
        int iy1 = (int)(y1 / _hx - _oy);

        for (int y = max(iy0, 0); y < min(iy1, _h); y++)
            for (int x = max(ix0, 0); x < min(ix1, _h); x++)
            {
                if (fabs(_src[x + y * _w]) < fabs(v))
                {
                    _src[x + y * _w] = v;
                }
            }
    }
};
