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
    double *_old; /* Contains old quantities at beginning of iteration */

    /* Distance field induced by solids.
     * Since this is used to compute the cell volumes, the samples are offset
     * by (-0.5, -0.5) from the samples in _src and the grid is one larger in
     * each dimension. This way, each sample of fluid quantity has four samples
     * of the distance function surrounding it - perfect for computing the
     * cell volumes.
     */
    double *_phi;
    /* Fractional cell volume occupied by fluid */
    double *_volume;

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

    /* Adds contribution `value' of sample at (x, y) to grid cell at (ix, iy)
     * using a hat filter.
     */
    void addSample(double *weight, double value, double x, double y, int ix, int iy)
    {
        if (ix < 0 || iy < 0 || ix >= _w || iy >= _h)
            return;

        double k = (1.0 - fabs(ix - x)) * (1.0 - fabs(iy - y));
        weight[ix + iy * _w] += k;
        _src[ix + iy * _w] += k * value;
    }

public:
    FluidQuantity(int w, int h, double ox, double oy, double hx)
        : _w(w), _h(h), _ox(ox), _oy(oy), _hx(hx)
    {
        _src = new double[_w * _h];
        _old = new double[_w * _h];

        _phi = new double[(_w + 1) * (_h + 1)];
        _volume = new double[_w * _h];
        _normalX = new double[_w * _h];
        _normalY = new double[_w * _h];

        _cell = new uint8_t[_w * _h];
        _body = new uint8_t[_w * _h];
        _mask = new uint8_t[_w * _h];

        for (int i = 0; i < _w * _h; i++)
        {
            _cell[i] = CELL_FLUID;
            _volume[i] = 1.0;
        }

        memset(_src, 0, _w * _h * sizeof(double));
    }

    ~FluidQuantity()
    {
        delete[] _src;
        delete[] _old;

        delete[] _phi;
        delete[] _volume;
        delete[] _normalX;
        delete[] _normalY;

        delete[] _cell;
        delete[] _body;
        delete[] _mask;
    }

    const double *src() const
    {
        return _src;
    }

    double *src()
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

    int idx(int x, int y) const
    {
        return x + y * _w;
    }

    double at(int x, int y) const
    {
        return _src[x + y * _w];
    }

    double volume(int x, int y) const
    {
        return _volume[x + y * _w];
    }

    double &at(int x, int y)
    {
        return _src[x + y * _w];
    }

    void copy()
    {
        memcpy(_old, _src, _w * _h * sizeof(double));
    }

    /* Computes the change in quantity during the last update */
    void diff(double alpha)
    {
        for (int i = 0; i < _w * _h; i++)
            _src[i] -= (1.0 - alpha) * _old[i];
    }

    /* Reverses the previous transformation - saves memory */
    void undiff(double alpha)
    {
        for (int i = 0; i < _w * _h; i++)
            _src[i] += (1.0 - alpha) * _old[i];
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

        /* Compute distance field first */
        for (int iy = 0, idx = 0; iy <= _h; iy++)
        {
            for (int ix = 0; ix <= _w; ix++, idx++)
            {
                double x = (ix + _ox - 0.5) * _hx;
                double y = (iy + _oy - 0.5) * _hx;

                _phi[idx] = bodies[0]->distance(x, y);
                for (unsigned i = 1; i < bodies.size(); i++)
                    _phi[idx] = min(_phi[idx], bodies[i]->distance(x, y));
            }
        }

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

                /* Compute cell volume from the four adjacent distance samples */
                int idxp = ix + iy * (_w + 1);
                _volume[idx] = 1.0 - occupancy(
                                         _phi[idxp], _phi[idxp + 1],
                                         _phi[idxp + _w + 1], _phi[idxp + _w + 2]);

                /* Clamp dangerously small cell volumes - could break numerical
                 * solver otherwise
                 */
                if (_volume[idx] < 0.01)
                    _volume[idx] = 0.0;

                bodies[_body[idx]]->distanceNormal(_normalX[idx], _normalY[idx], x, y);

                /* Solid cells are now defined as cells with zero fluid volume */
                if (_volume[idx] == 0.0)
                    _cell[idx] = CELL_SOLID;
                else
                    _cell[idx] = CELL_FLUID;
            }
        }
    }

    /* The extrapolation routine is now augmented to also fill in values for
     * cells that ended up with no particles in them. These are marked with
     * CELL_EMPTY. Empty cells are computed as the average value of all
     * available neighbours, and can therefore be computed as soon as at
     * least one neighbouring cell is available.
     */
    void fillSolidMask()
    {
        /* Make sure border is not touched by extrapolation - will be
         * handled separately.
         */
        for (int x = 0; x < _w; x++)
            _mask[x] = _mask[x + (_h - 1) * _w] = 0xFF;
        for (int y = 0; y < _h; y++)
            _mask[y * _w] = _mask[y * _w + _w - 1] = 0xFF;

        for (int y = 1; y < _h - 1; y++)
        {
            for (int x = 1; x < _w - 1; x++)
            {
                int idx = x + y * _w;

                _mask[idx] = 0;
                if (_cell[idx] == CELL_SOLID)
                {
                    double nx = _normalX[idx];
                    double ny = _normalY[idx];

                    if (nx != 0.0 && _cell[idx + sgn(nx)] != CELL_FLUID)
                        _mask[idx] |= 1;
                    if (ny != 0.0 && _cell[idx + sgn(ny) * _w] != CELL_FLUID)
                        _mask[idx] |= 2;
                }
                else if (_cell[idx] == CELL_EMPTY)
                {
                    /* Empty cells with no available neighbours need to be
                     * processed later.
                     */
                    _mask[idx] =
                        _cell[idx - 1] != CELL_FLUID &&
                        _cell[idx + 1] != CELL_FLUID &&
                        _cell[idx - _w] != CELL_FLUID &&
                        _cell[idx + _w] != CELL_FLUID;
                }
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

    /* Computes the extrapolated value as the average of all available
     * neighbouring cells.
     */
    double extrapolateAverage(int idx)
    {
        double value = 0.0;
        int count = 0;

        if (_cell[idx - 1] == CELL_FLUID)
        {
            value += _src[idx - 1];
            count++;
        }
        if (_cell[idx + 1] == CELL_FLUID)
        {
            value += _src[idx + 1];
            count++;
        }
        if (_cell[idx - _w] == CELL_FLUID)
        {
            value += _src[idx - _w];
            count++;
        }
        if (_cell[idx + _w] == CELL_FLUID)
        {
            value += _src[idx + _w];
            count++;
        }
        return value / count;
    }

    void freeSolidNeighbour(int idx, stack<int> &border, int mask)
    {
        if (_cell[idx] == CELL_SOLID)
        {
            _mask[idx] &= ~mask;
            if (_mask[idx] == 0)
                border.push(idx);
        }
    }

    /* At least one free neighbour cell is enough to add this cell to the queue
     * of ready cells.
     */
    void freeEmptyNeighbour(int idx, stack<int> &border)
    {
        if (_cell[idx] == CELL_EMPTY)
        {
            if (_mask[idx] == 1)
            {
                _mask[idx] = 0;
                border.push(idx);
            }
        }
    }

    /* For empty cells on the border of the simulation domain, we simply copy
     * the values of the adjacent cells.
     */
    void extrapolateEmptyBorders()
    {
        for (int x = 1; x < _w - 1; x++)
        {
            int idxT = x;
            int idxB = x + (_h - 1) * _w;

            if (_cell[idxT] == CELL_EMPTY)
                _src[idxT] = _src[idxT + _w];
            if (_cell[idxB] == CELL_EMPTY)
                _src[idxB] = _src[idxB - _w];
        }

        for (int y = 1; y < _h - 1; y++)
        {
            int idxL = y * _w;
            int idxR = y * _w + _w - 1;

            if (_cell[idxL] == CELL_EMPTY)
                _src[idxL] = _src[idxL + 1];
            if (_cell[idxR] == CELL_EMPTY)
                _src[idxR] = _src[idxR - 1];
        }

        int idxTL = 0;
        int idxTR = _w - 1;
        int idxBL = (_h - 1) * _w;
        int idxBR = _h * _w - 1;

        /* Corner cells average the values of the two adjacent border cells */
        if (_cell[idxTL] == CELL_EMPTY)
            _src[idxTL] = 0.5 * (_src[idxTL + 1] + _src[idxTL + _w]);
        if (_cell[idxTR] == CELL_EMPTY)
            _src[idxTR] = 0.5 * (_src[idxTR - 1] + _src[idxTR + _w]);
        if (_cell[idxBL] == CELL_EMPTY)
            _src[idxBL] = 0.5 * (_src[idxBL + 1] + _src[idxBL - _w]);
        if (_cell[idxBR] == CELL_EMPTY)
            _src[idxBR] = 0.5 * (_src[idxBR - 1] + _src[idxBR - _w]);

        for (int i = 0; i < _w * _h; i++)
            if (_cell[i] == CELL_EMPTY)
                _cell[i] = CELL_FLUID;
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

            if (_cell[idx] == CELL_EMPTY)
            {
                _src[idx] = extrapolateAverage(idx);
                _cell[idx] = CELL_FLUID; /* Mark extrapolated empty cells as fluid */
            }
            else
                _src[idx] = extrapolateNormal(idx);

            if (_normalX[idx - 1] > 0.0)
                freeSolidNeighbour(idx - 1, border, 1);
            if (_normalX[idx + 1] < 0.0)
                freeSolidNeighbour(idx + 1, border, 1);
            if (_normalY[idx - _w] > 0.0)
                freeSolidNeighbour(idx - _w, border, 2);
            if (_normalY[idx + _w] < 0.0)
                freeSolidNeighbour(idx + _w, border, 2);

            /* Notify adjacent empty cells */
            freeEmptyNeighbour(idx - 1, border);
            freeEmptyNeighbour(idx + 1, border);
            freeEmptyNeighbour(idx - _w, border);
            freeEmptyNeighbour(idx + _w, border);
        }
        extrapolateEmptyBorders();
    }

    /* Transfers particle values onto grid using a linear filter.
     *
     * In a first step, particle values and filter weights are accumulated on
     * the grid by looping over all particles and adding the particle contribution
     * to the four closest grid cells.
     *
     * In a second step, the actual grid values are obtained by dividing by the
     * filter weights. Cells with weight zero are cells which do not contain any
     * particles and are subsequently marked as empty for extrapolation.
     */
    void fromParticles(double *weight, int count, double *posX, double *posY, double *property)
    {
        memset(_src, 0, _w * _h * sizeof(double));
        memset(weight, 0, _w * _h * sizeof(double));

        for (int i = 0; i < count; i++)
        {
            double x = posX[i] - _ox;
            double y = posY[i] - _oy;
            x = max(0.5, min(_w - 1.5, x));
            y = max(0.5, min(_h - 1.5, y));

            int ix = (int)x;
            int iy = (int)y;

            addSample(weight, property[i], x, y, ix + 0, iy + 0);
            addSample(weight, property[i], x, y, ix + 1, iy + 0);
            addSample(weight, property[i], x, y, ix + 0, iy + 1);
            addSample(weight, property[i], x, y, ix + 1, iy + 1);
        }

        for (int i = 0; i < _w * _h; i++)
        {
            if (weight[i] != 0.0)
                _src[i] /= weight[i];
            else if (_cell[i] == CELL_FLUID)
                _cell[i] = CELL_EMPTY;
        }
    }
};
