#pragma once

#include <iostream>
#include <algorithm>
#include <math.h>
#include "Utils.h"

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

    /* Cubic intERPolate using samples a through d for x ranging from 0 to 1.
     * A Catmull-Rom spline is used. Over- and undershoots are clamped to
     * prevent blow-up.
     */
    double cerp(double a, double b, double c, double d, double x) const
    {
        double xsq = x * x;
        double xcu = xsq * x;

        double minV = min(a, min(b, min(c, d)));
        double maxV = max(a, max(b, max(c, d)));

        double t =
            a * (0.0 - 0.5 * x + 1.0 * xsq - 0.5 * xcu) +
            b * (1.0 + 0.0 * x - 2.5 * xsq + 1.5 * xcu) +
            c * (0.0 + 0.5 * x + 2.0 * xsq - 1.5 * xcu) +
            d * (0.0 + 0.0 * x - 0.5 * xsq + 0.5 * xcu);

        return min(max(t, minV), maxV);
    }

    /* Third order Runge-Kutta for velocity integration in time */
    void rungeKutta3(double &x, double &y, double timestep, const FluidQuantity &u, const FluidQuantity &v) const
    {
        double firstU = u.lerp(x, y) / _hx;
        double firstV = v.lerp(x, y) / _hx;

        double midX = x - 0.5 * timestep * firstU;
        double midY = y - 0.5 * timestep * firstV;

        double midU = u.lerp(midX, midY) / _hx;
        double midV = v.lerp(midX, midY) / _hx;

        double lastX = x - 0.75 * timestep * midU;
        double lastY = y - 0.75 * timestep * midV;

        double lastU = u.lerp(lastX, lastY);
        double lastV = v.lerp(lastX, lastY);

        x -= timestep * ((2.0 / 9.0) * firstU + (3.0 / 9.0) * midU + (4.0 / 9.0) * lastU);
        y -= timestep * ((2.0 / 9.0) * firstV + (3.0 / 9.0) * midV + (4.0 / 9.0) * lastV);
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

    /* Cubic intERPolate on grid at coordinates (x, y).
     * Coordinates will be clamped to lie in simulation domain
     */
    double cerp(double x, double y) const
    {
        x = min(max(x - _ox, 0.0), _w - 1.001);
        y = min(max(y - _oy, 0.0), _h - 1.001);
        int ix = (int)x;
        int iy = (int)y;
        x -= ix;
        y -= iy;

        int x0 = max(ix - 1, 0), x1 = ix, x2 = ix + 1, x3 = min(ix + 2, _w - 1);
        int y0 = max(iy - 1, 0), y1 = iy, y2 = iy + 1, y3 = min(iy + 2, _h - 1);

        double q0 = cerp(at(x0, y0), at(x1, y0), at(x2, y0), at(x3, y0), x);
        double q1 = cerp(at(x0, y1), at(x1, y1), at(x2, y1), at(x3, y1), x);
        double q2 = cerp(at(x0, y2), at(x1, y2), at(x2, y2), at(x3, y2), x);
        double q3 = cerp(at(x0, y3), at(x1, y3), at(x2, y3), at(x3, y3), x);

        return cerp(q0, q1, q2, q3, y);
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
                rungeKutta3(x, y, timestep, u, v);

                /* Second component: Interpolate from grid */
                _dst[idx] = cerp(x, y);
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
                double l = length(
                    (2.0 * (x + 0.5) * _hx - (x0 + x1)) / (x1 - x0),
                    (2.0 * (y + 0.5) * _hx - (y0 + y1)) / (y1 - y0));
                double vi = cubicPulse(l) * v;
                if (fabs(_src[x + y * _w]) < fabs(vi))
                    _src[x + y * _w] = vi;
            }
    }
};
