#pragma once

#include <iostream>
#include <math.h>

using namespace std;

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