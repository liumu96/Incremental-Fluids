#pragma once

#include "Utils.h"

class SolidBody
{
protected:
    double _posX; /* Position */
    double _posY;
    double _scaleX; /* Scale */
    double _scaleY;
    double _theta; /* Rotation */

    double _velX; /* Lateral velocity */
    double _velY;
    double _velTheta; /* Angular velocity */

    /* Transforms point (x, y) form the global to the local coordinate system */
    void globalToLocal(double &x, double &y) const
    {
        x -= _posX;
        y -= _posY;
        rotate(x, y, -_theta);
        x /= _scaleX;
        y /= _scaleY;
    }

    /* Transforms point (x, y) form the local to the global coordinate system */
    void localToGlobal(double &x, double &y) const
    {
        x *= _scaleX;
        y *= _scaleY;
        rotate(x, y, _theta);
        x += _posX;
        y += _posY;
    }

    SolidBody(double posX, double posY, double scaleX, double scaleY,
              double theta, double velX, double velY, double velTheta) : _posX(posX), _posY(posY), _scaleX(scaleX), _scaleY(scaleY),
                                                                         _theta(theta), _velX(velX), _velY(velY), _velTheta(velTheta) {}

    virtual ~SolidBody(){};

public:
    /* Returns the signed distance from (x, y) to the nearest point on surface
     * of the solid. The distance is negative if (x, y) is inside the solid
     */
    virtual double distance(double x, double y) const = 0;
    /* Changes (x, y) to lie on the closest point on the surface of the solid */
    virtual void closestSurfacePoint(double &x, double &y) const = 0;
    /* Returns the gradient of the distance function at (x, y) in (nx, ny) */
    virtual void distanceNormal(double &nx, double &ny, double x, double y) const = 0;

    /* Evaluates velocities of the solid at a given point */
    double velocityX(double x, double y) const
    {
        return (_posY - y) * _velTheta + _velX;
    }

    double velocityY(double x, double y) const
    {
        return (x - _posX) * _velTheta + _velY;
    }

    void velocity(double &vx, double &vy, double x, double y) const
    {
        vx = velocityX(x, y);
        vy = velocityY(x, y);
    }

    void update(double timestep)
    {
        /* Simple Euler integration - enough for solid bodies, since they
         * are not influenced by the simulation and velocities are typically
         * static
         */
        _posX += _velX * timestep;
        _posY += _velY * timestep;
        _theta += _velTheta * timestep;
    }
};

/* Represents a box (square) of size 1x1. Can be scaled to the appropriate size */
class SolidBox : public SolidBody
{
public:
    SolidBox(double x, double y, double sx, double sy, double t, double vx, double vy, double vt) : SolidBody(x, y, sx, sy, t, vx, vy, vt) {}

    double distance(double x, double y) const
    {
        x -= _posX;
        y -= _posY;
        rotate(x, y, -_theta);
        double dx = fabs(x) - _scaleX * 0.5;
        double dy = fabs(y) - _scaleY * 0.5;

        if (dx >= 0.0 || dy >= 0.0)
            return length(max(dx, 0.0), max(dy, 0.0));
        else
            return max(dx, dy);
    }

    void closestSurfacePoint(double &x, double &y) const
    {
        x -= _posX;
        y -= _posY;
        rotate(x, y, -_theta);
        double dx = fabs(x) - _scaleX * 0.5;
        double dy = fabs(y) - _scaleY * 0.5;

        if (dx > dy)
            x = nsgn(x) * 0.5 * _scaleX;
        else
            y = nsgn(y) * 0.5 * _scaleY;

        rotate(x, y, _theta);
        x += _posX;
        y += _posY;
    }

    void distanceNormal(double &nx, double &ny, double x, double y) const
    {
        x -= _posX;
        y -= _posY;
        rotate(x, y, -_theta);
        if (fabs(x) - _scaleX * 0.5 > fabs(y) - _scaleY * 0.5)
        {
            nx = nsgn(x);
            ny = 0.0;
        }
        else
        {
            nx = 0.0;
            ny = nsgn(y);
        }
        rotate(nx, ny, _theta);
    }
};

/* Represents a sphere (circle) of diameter 1. Can be scaled to the appropriate size */
class SolidSphere : public SolidBody
{
public:
    SolidSphere(double x, double y, double s, double t, double vx, double vy, double vt) : SolidBody(x, y, s, s, t, vx, vy, vt) {}

    double distance(double x, double y) const
    {
        return length(x - _posX, y - _posY) - _scaleX * 0.5;
    }

    void closestSurfacePoint(double &x, double &y) const
    {
        globalToLocal(x, y);

        double r = length(x, y);
        if (r < 1e-4)
        {
            x = 0.5;
            y = 0.0;
        }
        else
        {
            x /= 2.0 * r;
            y /= 2.0 * r;
        }

        localToGlobal(x, y);
    }

    void distanceNormal(double &nx, double &ny, double x, double y) const
    {
        x -= _posX;
        y -= _posY;
        float r = length(x, y);
        if (r < 1e-4)
        {
            nx = 1.0;
            ny = 0.0;
        }
        else
        {
            nx = x / r;
            ny = y / r;
        }
    }
};
