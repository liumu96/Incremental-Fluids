## Better Advection

### Basic Algorithms

Start with an initial divergence-free velocity field $\vec u^0$.

- For time step n = 0,1,2,...
- Determine a good time step ∆t to go from time $t_n$ to time $t_{n+1}$.
  - Set $\vec u^A = advect(\vec u^n, \Delta t, \vec u^n)$
  - Add $\vec u^B = \vec u^A + \Delta t \vec g$
  - Set $\vec u^{n+1} = project(\Delta t, \vec u^B)$

### Advection Algorithm

**Semi-Lagrangian Advection**

$$
\vec x_P = \vec x_G - \Delta t \vec u(\vec x_G)
$$

$$
q_G^{n+1} = interpolate(q^n, \vec x_P)
$$

To improve interpolation, we will replace the linear interpolation with a **cubic Catmull-Rom interpolation spline**.

To improve interpolation, forward Euler is replaced with a **third-order Runge-Kutta method**.

One of the things that should be noted is that the Catmull-Rom spline can lead to **oscillations**. One of the main conditions of the method is that the function which is being interpolated is smooth. When this is not the case, such as at borders of the fluid domain or near sharp details in the fluid flow, the interpolation spline can over- or undershoot and lead to curious stripes in the fluid. For now, this is something we have to deal with, but we will be able to improve this situation once we introduce FLIP.

One of the consequences of this is that we will have to **slightly tweak our inflows** - instead of setting fluid quantities to a specific value inside a rectangle, we place smooth "blobs" instead to ensure the fluid quantities near the inflow stay smooth to avoid oscillation.

### Code Implementation

In `FluidQuantity` class:

**Function `addInflow`**

```cpp
/* Set fluid quantity inside the given rect to the specified value, but use
 * a smooth falloff to avoid oscillations
 */
void addInflow(double x0, double y0, double x1, double y1, double v) {
    int ix0 = (int)(x0/_hx - _ox);
    int iy0 = (int)(y0/_hx - _oy);
    int ix1 = (int)(x1/_hx - _ox);
    int iy1 = (int)(y1/_hx - _oy);

    for (int y = max(iy0, 0); y < min(iy1, _h); y++) {
        for (int x = max(ix0, 0); x < min(ix1, _h); x++) {
            double l = length(
                (2.0*(x + 0.5)*_hx - (x0 + x1))/(x1 - x0),
                (2.0*(y + 0.5)*_hx - (y0 + y1))/(y1 - y0)
            );
            double vi = cubicPulse(l)*v;
            if (fabs(_src[x + y*_w]) < fabs(vi))
                _src[x + y*_w] = vi;
        }
    }
}
```

```cpp
/* Cubic pulse function.
 * Returns a value in range [0, 1].
 * Return value is 0 for x <= -1 and x >= 1; value is 1 for x=0
 * Smoothly interpolates between 0 and 1 between these three points.
 */
double cubicPulse(double x) {
    x = min(fabs(x), 1.0);
    return 1.0 - x*x*(3.0 - 2.0*x);
}
```

**Function `Advect`**

```cpp
/* Advect grid in velocity field u, v with given timestep */
void advect(double timestep, const FluidQuantity &u, const FluidQuantity &v) {
    for (int iy = 0, idx = 0; iy < _h; iy++) {
        for (int ix = 0; ix < _w; ix++, idx++) {
            double x = ix + _ox;
            double y = iy + _oy;

            /* First component: Integrate in time */
            rungeKutta3(x, y, timestep, u, v);

            /* Second component: Interpolate from grid */
            _dst[idx] = cerp(x, y);
        }
    }
}
```

**cubic Catmull-Rom interpolation spline**

```cpp
/* Cubic intERPolate using samples a through d for x ranging from 0 to 1.
 * A Catmull-Rom spline is used. Over- and undershoots are clamped to
 * prevent blow-up.
 */
double cerp(double a, double b, double c, double d, double x) const {
    double xsq = x*x;
    double xcu = xsq*x;

    double minV = min(a, min(b, min(c, d)));
    double maxV = max(a, max(b, max(c, d)));

    double t =
        a*(0.0 - 0.5*x + 1.0*xsq - 0.5*xcu) +
        b*(1.0 + 0.0*x - 2.5*xsq + 1.5*xcu) +
        c*(0.0 + 0.5*x + 2.0*xsq - 1.5*xcu) +
        d*(0.0 + 0.0*x - 0.5*xsq + 0.5*xcu);

    return min(max(t, minV), maxV);
}

/* Cubic intERPolate on grid at coordinates (x, y).
 * Coordinates will be clamped to lie in simulation domain
 */
double cerp(double x, double y) const {
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
```

**third-order Runge-Kutta method**

$$
k_1 = f(q^n) \\
$$

$$
k_2 = f(q^n + \frac{1}{2}\Delta t k_1) \\
$$

$$
k_3 = f(q^n + \frac{3}{4}\Delta t k_2) \\
$$

$$
q^{n+1} = q^n + \frac{2}{9}\Delta t k_1 + \frac{3}{9}\Delta t k_2 + \frac{4}{9}\Delta t k_4
$$

````cpp
/* Third order Runge-Kutta for velocity integration in time */
void rungeKutta3(double &x, double &y, double timestep, const FluidQuantity &u, const FluidQuantity &v) const {
    double firstU = u.lerp(x, y)/_hx;
    double firstV = v.lerp(x, y)/_hx;

    double midX = x - 0.5*timestep*firstU;
    double midY = y - 0.5*timestep*firstV;

    double midU = u.lerp(midX, midY)/_hx;
    double midV = v.lerp(midX, midY)/_hx;

    double lastX = x - 0.75*timestep*midU;
    double lastY = y - 0.75*timestep*midV;

    double lastU = u.lerp(lastX, lastY);
    double lastV = v.lerp(lastX, lastY);

    x -= timestep*((2.0/9.0)*firstU + (3.0/9.0)*midU + (4.0/9.0)*lastU);
    y -= timestep*((2.0/9.0)*firstV + (3.0/9.0)*midV + (4.0/9.0)*lastV);
}
```
## Render Result

<div align="center">
<img src="./result/0001-0400.gif" width="30%"/>
</div>
````
