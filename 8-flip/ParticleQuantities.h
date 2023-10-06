#pragma once

#include "FluidQuantity.h"

/* Main class processing fluid particles */
class ParticleQuantities
{
private:
    /* Maximum allowed number of particles per cell */
    static const int _MaxPerCell = 12;
    /* Minimum allowed number of particles per cell */
    static const int _MinPerCell = 3;
    /* Initial number of particles per cell */
    static const int _AvgPerCell = 4;

    /* Number of particles currently active */
    int _particleCount;
    /* Maximum number of particles the simulation can handle */
    int _maxParticles;

    /* The usual culprits */
    int _w;
    int _h;
    double _hx;
    const vector<const SolidBody *> &_bodies;

    /* Filter weights (auxiliary array provided to fluid quantities) */
    double *_weight;
    /* Number of particles per cell */
    int *_counts;

    /* Particle positions */
    double *_posX;
    double *_posY;
    /* Particle 'properties', that is, value for each fluid quantity
     * (velocities, density etc.)
     */
    vector<double *> _properties;
    vector<FluidQuantity *> _quantities;

    /* Helper function returning true if a position is inside a solid body */
    bool pointInBody(double x, double y)
    {
        for (unsigned i = 0; i < _bodies.size(); i++)
            if (_bodies[i]->distance(x * _hx, y * _hx) < 0.0)
                return true;

        return false;
    }

    /* Initializes particle positions on randomly jittered grid locations */
    void initParticles()
    {
        int idx = 0;
        for (int y = 0; y < _h; y++)
        {
            for (int x = 0; x < _w; x++)
            {
                for (int i = 0; i < _AvgPerCell; i++, idx++)
                {
                    _posX[idx] = x + frand();
                    _posY[idx] = y + frand();

                    /* Discard particles landing inside solid bodies */
                    if (pointInBody(_posX[idx], _posY[idx]))
                        idx--;
                }
            }
        }

        _particleCount = idx;
    }

    /* Counts the number of particles per cell */
    void countParticles()
    {
        memset(_counts, 0, _w * _h * sizeof(int));
        for (int i = 0; i < _particleCount; i++)
        {
            int ix = (int)_posX[i];
            int iy = (int)_posY[i];

            if (ix >= 0 && iy >= 0 && ix < _w && iy < _h)
                _counts[ix + iy * _w]++;
        }
    }

    /* Decimates particles in crowded cells */
    void pruneParticles()
    {
        for (int i = 0; i < _particleCount; i++)
        {
            int ix = (int)_posX[i];
            int iy = (int)_posY[i];
            int idx = ix + iy * _w;

            if (ix < 0 || iy < 0 || ix >= _w || iy >= _h)
                continue;

            if (_counts[idx] > _MaxPerCell)
            {
                int j = --_particleCount;
                _posX[i] = _posX[j];
                _posY[i] = _posY[j];
                for (unsigned t = 0; t < _quantities.size(); t++)
                    _properties[t][i] = _properties[t][j];

                _counts[idx]--;
                i--;
            }
        }
    }

    /* Adds new particles in cells with dangerously little particles */
    void seedParticles()
    {
        for (int y = 0, idx = 0; y < _h; y++)
        {
            for (int x = 0; x < _w; x++, idx++)
            {
                for (int i = 0; i < _MinPerCell - _counts[idx]; i++)
                {
                    if (_particleCount == _maxParticles)
                        return;

                    int j = _particleCount;

                    _posX[j] = x + frand();
                    _posY[j] = y + frand();

                    /* Reject particle if it lands inside a solid body */
                    if (pointInBody(_posX[idx], _posY[idx]))
                        continue;

                    /* Get current grid values */
                    for (unsigned t = 0; t < _quantities.size(); t++)
                        _properties[t][j] = _quantities[t]->lerp(_posX[j], _posY[j]);

                    _particleCount++;
                }
            }
        }
    }

    /* Pushes particle back into the fluid if they land inside solid bodies */
    void backProject(double &x, double &y)
    {
        double d = 1e30;
        int closestBody = -1;
        for (unsigned i = 0; i < _bodies.size(); i++)
        {
            double id = _bodies[i]->distance(x * _hx, y * _hx);

            if (id < d)
            {
                d = id;
                closestBody = i;
            }
        }

        if (d < -1.0)
        {
            x *= _hx;
            y *= _hx;
            _bodies[closestBody]->closestSurfacePoint(x, y);
            double nx, ny;
            _bodies[closestBody]->distanceNormal(nx, ny, x, y);
            x -= nx * _hx;
            y -= ny * _hx;
            x /= _hx;
            y /= _hx;
        }
    }

    /* The same Runge Kutta interpolation routine as before - only now forward
     * in time instead of backwards.
     */
    void rungeKutta3(double &x, double &y, double timestep, const FluidQuantity &u, const FluidQuantity &v) const
    {
        double firstU = u.lerp(x, y) / _hx;
        double firstV = v.lerp(x, y) / _hx;

        double midX = x + 0.5 * timestep * firstU;
        double midY = y + 0.5 * timestep * firstV;

        double midU = u.lerp(midX, midY) / _hx;
        double midV = v.lerp(midX, midY) / _hx;

        double lastX = x + 0.75 * timestep * midU;
        double lastY = y + 0.75 * timestep * midV;

        double lastU = u.lerp(lastX, lastY);
        double lastV = v.lerp(lastX, lastY);

        x += timestep * ((2.0 / 9.0) * firstU + (3.0 / 9.0) * midU + (4.0 / 9.0) * lastU);
        y += timestep * ((2.0 / 9.0) * firstV + (3.0 / 9.0) * midV + (4.0 / 9.0) * lastV);
    }

public:
    ParticleQuantities(int w, int h, double hx,
                       const vector<const SolidBody *> &bodies) : _w(w), _h(h), _hx(hx), _bodies(bodies)
    {

        _maxParticles = _w * _h * _MaxPerCell;

        _posX = new double[_maxParticles];
        _posY = new double[_maxParticles];

        _weight = new double[(_w + 1) * (_h + 1)];
        _counts = new int[_w * _h];

        initParticles();
    }

    ~ParticleQuantities()
    {
        delete[] _posX;
        delete[] _posY;

        delete[] _weight;
        delete[] _counts;

        for (size_t i = 0; i < _quantities.size(); ++i)
            delete[] _quantities[i];
    }

    /* Adds a new quantity to be carried by the particles */
    void addQuantity(FluidQuantity *q)
    {
        double *property = new double[_maxParticles];
        memset(property, 0, _maxParticles * sizeof(double));

        _quantities.push_back(q);
        _properties.push_back(property);
    }

    /* Interpolates the change in quantity back onto the particles.
     * Mixes in a little bit of the pure Particle-in-cell update using the
     * parameter alpha.
     */
    void gridToParticles(double alpha)
    {
        for (unsigned t = 0; t < _quantities.size(); t++)
        {
            for (int i = 0; i < _particleCount; i++)
            {
                _properties[t][i] *= 1.0 - alpha;
                _properties[t][i] += _quantities[t]->lerp(_posX[i], _posY[i]);
            }
        }
    }

    /* Interpolates particle quantities onto the grid, extrapolates them and
     * spawns/prunes particles where necessary.
     */
    void particlesToGrid()
    {
        for (unsigned t = 0; t < _quantities.size(); t++)
        {
            _quantities[t]->fromParticles(_weight, _particleCount, _posX, _posY, _properties[t]);
            _quantities[t]->extrapolate();
        }

        countParticles();
        pruneParticles();
        seedParticles();

        printf("Particle count: %d\n", _particleCount);
    }

    /* Advects particle in velocity field and clamps resulting positions to
     * the fluid domain */
    void advect(double timestep, const FluidQuantity &u, const FluidQuantity &v)
    {
        for (int i = 0; i < _particleCount; i++)
        {
            rungeKutta3(_posX[i], _posY[i], timestep, u, v);
            backProject(_posX[i], _posY[i]);

            _posX[i] = max(min(_posX[i], _w - 0.001), 0.0);
            _posY[i] = max(min(_posY[i], _h - 0.001), 0.0);
        }
    }
};
