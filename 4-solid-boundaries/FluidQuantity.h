#pragma once

#include <iostream>
#include <algorithm>
#include <math.h>
#include "Utils.h"
#include "SolidBody.h"
#include <stack>

using namespace std;

class FluidQuantity
{
private:
    // Memory buffers for fluid quantity
    double *_src;
    double *_dst;

    /* Normal of distance field at grid points */
    double *_normalX;
    double *_normalY;
    /* Designates cells as fluid or solid cells (CELL_FLUID or CELL_SOLID) */
    uint8_t *_cell;
    /* Specifies the index of the rigid body closes to a grid cell */
    uint8_t *_body;
    /* Auxiliary array used for extrapolation */
    uint8_t *_mask;

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
    FluidQuantity(int w, int h, double ox, double oy, double hx)
        : _w(w), _h(h), _ox(ox), _oy(oy), _hx(hx)
    {
        _src = new double[_w * _h];
        _dst = new double[_w * _h];

        _normalX = new double[_w * _h];
        _normalY = new double[_w * _h];

        _cell = new uint8_t[_w * _h];
        _body = new uint8_t[_w * _h];
        _mask = new uint8_t[_w * _h];

        memset(_cell, 0, _w * _h * sizeof(uint8_t));
        memset(_src, 0, _w * _h * sizeof(double));
    }

    ~FluidQuantity()
    {
        delete[] _src;
        delete[] _dst;

        delete[] _normalX;
        delete[] _normalY;

        delete[] _cell;
        delete[] _body;
        delete[] _mask;
    }

    void flip()
    {
        swap(_src, _dst);
    }

    const double *src() const
    {
        return _src;
    }

    const uint8_t *cell() const
    {
        return _cell;
    }

    const uint8_t *body() const
    {
        return _body;
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

    /* If the point (x, y) is inside a solid, project it back out to the
     * closest point on the surface of the solid.
     */
    void backProject(double &x, double &y, const vector<const SolidBody *> &bodies)
    {
        int rx = min(max((int)(x - _ox), 0), _w - 1);
        int ry = min(max((int)(y - _oy), 0), _h - 1);

        if (_cell[rx + ry * _w] != CELL_FLUID)
        {
            x = (x - _ox) * _hx;
            y = (y - _oy) * _hx;
            bodies[_body[rx + ry * _w]]->closestSurfacePoint(x, y);
            x = x / _hx + _ox;
            y = y / _hx + _oy;
        }
    }

    /* Advect grid in velocity field u, v with given timestep */
    void advect(double timestep, const FluidQuantity &u, const FluidQuantity &v, const vector<const SolidBody *> &bodies)
    {
        for (int iy = 0, idx = 0; iy < _h; iy++)
            for (int ix = 0; ix < _w; ix++, idx++)
            {
                if (_cell[idx] == CELL_FLUID)
                {
                    double x = ix + _ox;
                    double y = iy + _oy;

                    rungeKutta3(x, y, timestep, u, v);

                    /* If integrating back in time leaves us inside a solid
                     * boundary (due to numerical error), make sure we
                     * interpolate from a point inside the fluid.
                     */
                    backProject(x, y, bodies);

                    _dst[idx] = cerp(x, y);
                }
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

    /* Fill all solid related fields - that is, _cell, _body and _normalX/Y */
    void fillSolidFields(const vector<const SolidBody *> &bodies)
    {
        if (bodies.empty())
            return;

        for (int iy = 0, idx = 0; iy < _h; iy++)
        {
            for (int ix = 0; ix < _w; ix++, idx++)
            {
                double x = (ix + _ox) * _hx;
                double y = (iy + _oy) * _hx;

                /* Search closest solid */
                _body[idx] = 0;
                double d = bodies[0]->distance(x, y);
                for (unsigned i = 1; i < bodies.size(); i++)
                {
                    double id = bodies[i]->distance(x, y);
                    if (id < d)
                    {
                        _body[idx] = i;
                        d = id;
                    }
                }

                /* If distance to closest solid is negative, this cell must be
                 * inside it
                 */
                if (d < 0.0)
                    _cell[idx] = CELL_SOLID;
                else
                    _cell[idx] = CELL_FLUID;

                bodies[_body[idx]]->distanceNormal(_normalX[idx], _normalY[idx], x, y);
            }
        }
    }

    /* Prepare auxiliary array for extrapolation.
     * The purpose of extrapolation is to extrapolate fluid quantities into
     * solids, where these quantities would normally be undefined. However, we
     * need these values for stable interpolation and boundary conditions.
     *
     * The way these are extrapolated here is by essentially solving a PDE,
     * such that the gradient of the fluid quantity is 0 along the gradient
     * of the distance field. This is essentially a more robust formulation of
     * "Set quantity inside solid to the value at the closest point on the
     * solid-fluid boundary"
     *
     * This PDE has a particular form which makes it very easy to solve exactly
     * using an upwinding scheme. What this means is that we can solve it from
     * outside-to-inside, with information flowing along the normal from the
     * boundary.
     *
     * Specifically, we can solve for the value inside a fluid cell using
     * extrapolateNormal if the two adjacent grid cells in "upstream" direction
     * (where the normal points to) are either fluid cells or have been solved
     * for already.
     *
     * The mask array keeps track of which cells wait for which values. If an
     * entry is 0, it means both neighbours are available and the cell is ready
     * for the PDE solve. If it is 1, the cell waits for the neighbour in
     * x direction, 2 for y-direction and 3 for both.
     */
    void fillSolidMask()
    {
        for (int y = 1; y < _h - 1; y++)
        {
            for (int x = 1; x < _w - 1; x++)
            {
                int idx = x + y * _w;

                if (_cell[idx] == CELL_FLUID)
                    continue;

                double nx = _normalX[idx];
                double ny = _normalY[idx];

                _mask[idx] = 0;
                if (nx != 0.0 && _cell[idx + sgn(nx)] != CELL_FLUID)
                    _mask[idx] |= 1; /* Neighbour in normal x direction is blocked */
                if (ny != 0.0 && _cell[idx + sgn(ny) * _w] != CELL_FLUID)
                    _mask[idx] |= 2; /* Neighbour in normal y direction is blocked */
            }
        }
    }
    /* Solve for value at index idx using values of neighbours in normal x/y
     * direction. The value is computed such that the directional derivative
     * along distance field normal is 0.
     */
    double extrapolateNormal(int idx)
    {
        double nx = _normalX[idx];
        double ny = _normalY[idx];

        double srcX = _src[idx + sgn(nx)];
        double srcY = _src[idx + sgn(ny) * _w];

        return (fabs(nx) * srcX + fabs(ny) * srcY) / (fabs(nx) + fabs(ny));
    }

    /* Given that a neighbour in upstream direction specified by mask (1=x, 2=y)
     * now has been solved for, update the mask appropriately and, if this cell
     * can now be computed, add it to the queue of ready cells
     */
    void freeNeighbour(int idx, stack<int> &border, int mask)
    {
        _mask[idx] &= ~mask;
        if (_cell[idx] != CELL_FLUID && _mask[idx] == 0)
            border.push(idx);
    }

    void extrapolate()
    {
        fillSolidMask();

        /* Queue of cells which can be computed */
        stack<int> border;
        /* Initialize queue by finding all solid cells with mask=0 (ready for
         * extrapolation)
         */
        for (int y = 1; y < _h - 1; y++)
        {
            for (int x = 1; x < _w - 1; x++)
            {
                int idx = x + y * _w;

                if (_cell[idx] != CELL_FLUID && _mask[idx] == 0)
                    border.push(idx);
            }
        }

        while (!border.empty())
        {
            int idx = border.top();
            border.pop();

            /* Solve for value in cell */
            _src[idx] = extrapolateNormal(idx);

            /* Notify adjacent cells that this cell has been computed and can
             * be used as an upstream value
             */
            if (_normalX[idx - 1] > 0.0)
                freeNeighbour(idx - 1, border, 1);
            if (_normalX[idx + 1] < 0.0)
                freeNeighbour(idx + 1, border, 1);
            if (_normalY[idx - _w] > 0.0)
                freeNeighbour(idx - _w, border, 2);
            if (_normalY[idx + _w] < 0.0)
                freeNeighbour(idx + _w, border, 2);
        }
    }
};
