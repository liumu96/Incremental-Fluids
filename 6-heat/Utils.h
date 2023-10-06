#pragma once

#include <iostream>
#include <math.h>

using namespace std;

/* See http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c */
template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

/* Non-zero sgn */
template <typename T>
int nsgn(T val)
{
    return (val < T(0) ? -1 : 1);
}

/* Length of vector (x, y) */
double length(double x, double y)
{
    return sqrt(x * x + y * y);
}

/* Cubic pulse function.
 * Returns a value in range [0, 1].
 * Return value is 0 for x <= -1 and x >= 1; value is 1 for x=0
 * Smoothly interpolates between 0 and 1 between these three points.
 */
double cubicPulse(double x)
{
    x = min(fabs(x), 1.0);
    return 1.0 - x * x * (3.0 - 2.0 * x);
}

/* Rotates point (x, y) by angle phi */
void rotate(double &x, double &y, double phi)
{
    double tmpX = x, tmpY = y;
    x = cos(phi) * tmpX + sin(phi) * tmpY;
    y = -sin(phi) * tmpX + cos(phi) * tmpY;
}

/* For three corners in a 1x1 square, with `in' being adjacent to `out1' and
 * `out2' and all three parameters being distances to a surface, `in' being
 * inside the surface and `out1' and `out2' outside, returns the area of the
 * square occupied by the surface.
 */
double triangleOccupancy(double out1, double in, double out2)
{
    return 0.5 * in * in / ((out1 - in) * (out2 - in));
}

/* For four corners in a 1x1 square, with all parameters being distances to a
 * surface and `in1' and `in2 inside the surface, returns the are of the square
 * occupied by the surface.
 */
double trapezoidOccupancy(double out1, double out2, double in1, double in2)
{
    return 0.5 * (-in1 / (out1 - in1) - in2 / (out2 - in2));
}

/* Given the distance of four corners in a 1x1 square to a surface, returns the
 * area of the part of the square occupied by the surface computed analytically.
 *
 * The basic workings of this algorithm are quite similar to marching squares
 * (2D marching cubes). First, a mask is computed based on which corners are
 * inside and which are outside.
 * Based on this mask, the function differentiates between one of four cases:
 * a) Only one corner is inside
 *   => Compute using triangle area
 * b) Only one corner is outside
 *   => Invert distance field, compute 1 - triangle area
 * c) Two adjacent corners are inside
 *   => Compute using trapezoid area
 * d) Two opposing corners are inside
 *   => Compute as sum of area of two opposed triangles
 *
 * The two remaining cases, all corners outside/inside, can be computed trivially
 */
double occupancy(double d11, double d12, double d21, double d22)
{
    double ds[] = {d11, d12, d22, d21};

    /* Compute mask */
    uint8_t b = 0;
    for (int i = 3; i >= 0; i--)
        b = (b << 1) | (ds[i] < 0.0 ? 1 : 0);

    switch (b)
    {
    /* All outside */
    case 0x0:
        return 0.0;
    /* One inside */
    case 0x1:
        return triangleOccupancy(d21, d11, d12);
    case 0x2:
        return triangleOccupancy(d11, d12, d22);
    case 0x4:
        return triangleOccupancy(d12, d22, d21);
    case 0x8:
        return triangleOccupancy(d22, d21, d11);
    /* One outside */
    case 0xE:
        return 1.0 - triangleOccupancy(-d21, -d11, -d12);
    case 0xD:
        return 1.0 - triangleOccupancy(-d11, -d12, -d22);
    case 0xB:
        return 1.0 - triangleOccupancy(-d12, -d22, -d21);
    case 0x7:
        return 1.0 - triangleOccupancy(-d22, -d21, -d11);
    /* Two adjacent inside */
    case 0x3:
        return trapezoidOccupancy(d21, d22, d11, d12);
    case 0x6:
        return trapezoidOccupancy(d11, d21, d12, d22);
    case 0x9:
        return trapezoidOccupancy(d12, d22, d11, d21);
    case 0xC:
        return trapezoidOccupancy(d11, d12, d21, d22);
    /* Two opposed inside */
    case 0x5:
        return triangleOccupancy(d11, d12, d22) +
               triangleOccupancy(d22, d21, d11);
    case 0xA:
        return triangleOccupancy(d21, d11, d12) +
               triangleOccupancy(d12, d22, d21);
    /* All inside */
    case 0xF:
        return 1.0;
    }

    return 0.0;
}

/* Enum to differentiate fluid and solid cells */
enum CellType
{
    CELL_FLUID,
    CELL_SOLID
};