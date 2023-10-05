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

/* Enum to differentiate fluid and solid cells */
enum CellType
{
    CELL_FLUID,
    CELL_SOLID
};