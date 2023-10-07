## Basic Algorithms

Start with an initial divergence-free velocity field $\vec u^0$.

- For time step n = 0,1,2,...
- Determine a good time step ∆t to go from time $t_n$ to time $t_{n+1}$.
  - Set $\vec u^A = advect(\vec u^n, \Delta t, \vec u^n)$
  - Add $\vec u^B = \vec u^A + \Delta t \vec g$
  - Set $\vec u^{n+1} = project(\Delta t, \vec u^B)$

## Algorithm && Code

### Matrixless Gauss-Seidel Solver

- Add sources
- Force incompressibility
  - Build pressure right-hand side
  - Build pressure matrix
  - Solve for pressure
  - Apply pressure
- Advect: Semi-Lagrangian Advection

  We’ll say that the location in space of the grid point we’re looking at is $\vec x_G$. We want to find the new value of q at that point, which we’ll call $q_G^{n+1}$. We know from our understanding of advection that if a hypothetical particle with old value $q_P^n$ ends up at $\vec x_G$ , when it moves through the velocity field for the time step $\Delta t$, then $q_G^{n+1} = q_P^n$ . So the question is, how do we figure out $q_P^n$ ?

  <div align="center">
    <img src=https://github.com/liumu96/Incremental-Fluids/blob/main/1-matrixless/image.png width=30%/>
  </div>
    To find a fluid value at grid point $\vec x_G$ at the new time step, we need to know where the fluid at was one time step ago, position $\vec x_P$ , following the velocity field.

The first step is figuring out where this imaginary particle would have started from, a position we’ll call $\vec x_P$. The particle moves according to the simple ordinary differential equation

$$
\frac{d\vec x}{dt} = \vec u(\vec x)
$$

and ends up at $\vec x_G$ after time $\Delta t$. If we now run time backwards, we can go in reverse from to the start point of the particle—i.e., finding where a particle would end up under the reverse velocity field $\vec u$ “starting” from . Figure (3.1) illustrates this path. The simplest possible way to estimated $\vec x_P$ is to use one step of“forward” Euler going backwards in time:

$$
\vec x_P = \vec x_G - \Delta t \vec u(\vec x_G)
$$

where we use the velocity $\vec u$ evaluated at the grid point to take a $\Delta t$-step backwards through the flow field. It turns out forward Euler is sometimes adequate, but significantly better results can be obtained using a slightly more sophisticated technique such as a higher-order `Runge-Kutta` method.

We now know the position where the imaginary particle started; next we have to figure out what old value of q it had. Most likely $\vec x_P$ is not on the grid, so we don’t have the exact value, but we can get a good approximation by interpolating from $q^n$ at nearby grid points. Trilinear (bilinear in two dimensions) interpolation is often used, though this comes with a serious penalty which we will fix at the end of the chapter.

Putting this together into a formula, our basic semi-Lagrangian formula, assuming the particle-tracing algorithm has tracked back to location $\vec x_P$(typically with RK2 above), is

$$
q_G^{n+1} = interpolate(q^n, \vec x_P)
$$

Note that the particle I’ve described is purely hypothetical. No particle is actually created in the computer: we simply use Lagrangian particles to conceptually figure out the update formula for the Eulerian advection step. Because we are almost using a Lagrangian approach to do an Eulerian calculation, this is called the semi-Lagrangian method.

Just for completeness, let’s illustrate this in one dimension again, using linear interpolation for the semi-Lagrangian operations. For grid point $x_i$, the particle is traced back to $x_P = x_i - \Delta t u$. Assuming this lies in the interval $[x_j, x_{j+1}]$, and letting $\alpha = (x_P - x_j)/\Delta x$ be the fraction of the interval the point lands in, the linear interpolation is $q_P^n = (1 - \alpha)q_j^n + \alpha q_{j+1}^n$ So our update is

$$
q_i^{n+1} = (1 - \alpha)q_j^n + \alpha q_{j+1}^n
$$

In practice we will need to advect the velocity field, and perhaps additional variables such as smoke density or temperature. Usually the additional variables are stored at the grid cell centers, but the velocity components are stored at the staggered grid locations discussed in the previous chapter. In each case, we will need to use the appropriate averaged velocity, given at the end of the previous chapter, to estimate the particle trajectory.

### Simplification:

- inflows: stay almost **the same** for all of our implementations.
- Adding Sources:
  - only add inflows for this implementation
  - no body forces are applied
  - For inflows, we simply set the fluid quantities inside a specified rectangle to a fixed value.
- Forcing Incompressibility
  - The pressure right hand side is quite standard in this case; it is simply **the negative divergence of the velocity fieldv**, which can be computed easily from the MAC grid.
  - As you can probably guess from the title, we **skip explicitly building the matrix for this implementation**, since it has a very simple structure. Instead, we will simply construct the matrix implicitly in the pressure solve, which keeps in line with most of the "simple fluid" literature out there.
  - For the pressure solve itself, we will use a **Gauss-Seidel solver** [5], which is a very simple iterative solver for linear systems of equations. It doesn't show excellent convergence, but it is simple to implement and requires little memory.
  - Finally, applying the pressure is performed by **subtracting the pressure gradient from the velocity field**, which can be easily performed on the MAC grid.

### FluidQuantity

This is the class representing fluid quantities such as density and velocity on the MAC grid. It saves attributes such as offset from the top left grid cell, grid width and height as well as cell size.

It also contains two memory buffers: A source `(_src)` buffer and a destination (\_dst) buffer.

Most operations on fluid quantities can be done in-place; that is, they write to the same buffer they're reading from (which is always `_src`). However, some operations, such as advection, cannot be done in-place. Instead, they will write to the `_dst` buffer. Once the operation is completed, `flip()` can be called to swap the source and destination buffers, such that the result of the operation is visible to subsequent operations.

```cpp
class FluidQuantity {
private:
    /* Memory buffers for fluid quantity */
    double *_src;
    double *_dst;

	/* Width and height */
    int _w;
    int _h;

    /* X and Y offset from top left grid cell.
    * This is (0.5,0.5) for centered quantities such as density,
    * and (0.0, 0.5) or (0.5, 0.0) for jittered quantities like the velocity.
    */
    double _ox;
    double _oy;

	/* Grid cell size */
    double _hx;
public:
	FluidQuantity(int w, int h, double ox, double oy, double hx)
            : _w(w), _h(h), _ox(ox), _oy(oy), _hx(hx) {
        _src = new double[_w*_h];
        _dst = new double[_w*_h];

        memset(_src, 0, _w*_h*sizeof(double));
    }

    ~FluidQuantity() {
        delete[] _src;
        delete[] _dst;
    }
}
```

**Core Function**

Function `addInflow`

```cpp
/* Sets fluid quantity inside the given rect to value `v' */
void addInflow(double x0, double y0, double x1, double y1, double v) {
    int ix0 = (int)(x0/_hx - _ox);
    int iy0 = (int)(y0/_hx - _oy);
    int ix1 = (int)(x1/_hx - _ox);
    int iy1 = (int)(y1/_hx - _oy);

    for (int y = max(iy0, 0); y < min(iy1, _h); y++)
        for (int x = max(ix0, 0); x < min(ix1, _h); x++)
            if (fabs(_src[x + y*_w]) < fabs(v))
                _src[x + y*_w] = v;
}
```

Function `advect`

$$
\vec x_P = \vec x_G - \Delta t \vec u(\vec x_G) \\
q_G^{n+1} = interpolate(q^n, \vec x_P) \\
q_i^{n+1} = (1 - \alpha)q_j^n + \alpha q_{j+1}^n
$$

```cpp
/* Linear intERPolate on grid at coordinates (x, y).
 * Coordinates will be clamped to lie in simulation domain
 */
double lerp(double x, double y) const {
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
/* Simple forward Euler method for velocity integration in time */
void euler(double &x, double &y, double timestep, const FluidQuantity &u, const FluidQuantity &v) const {
    double uVel = u.lerp(x, y)/_hx;
    double vVel = v.lerp(x, y)/_hx;

    x -= uVel*timestep;
    y -= vVel*timestep;
}
/* Advect grid in velocity field u, v with given timestep */
void advect(double timestep, const FluidQuantity &u, const FluidQuantity &v) {
    for (int iy = 0, idx = 0; iy < _h; iy++) {
        for (int ix = 0; ix < _w; ix++, idx++) {
            double x = ix + _ox;
            double y = iy + _oy;

            /* First component: Integrate in time */
            euler(x, y, timestep, u, v);

            /* Second component: Interpolate from grid */
            _dst[idx] = lerp(x, y);
        }
    }
}
```

### FluidSolver

Fluid solver class. Sets up the fluid quantities, forces incompressibility performs advection and adds inflows.

```cpp
class FluidSolver {
private:
    /* Fluid quantities */
    FluidQuantity *_d;
    FluidQuantity *_u;
    FluidQuantity *_v;

    /* Width and height */
    int _w;
    int _h;

    /* Grid cell size and fluid density */
    double _hx;
    double _density;

    /* Arrays for: */
    double *_r; /* Right hand side of pressure solve */
    double *_p; /* Pressure solution */
public:
    /* Set density and x/y velocity in given rectangle to d/u/v, respectively */
    void addInflow(double x, double y, double w, double h, double d, double u, double v) {
        _d->addInflow(x, y, x + w, y + h, d);
        _u->addInflow(x, y, x + w, y + h, u);
        _v->addInflow(x, y, x + w, y + h, v);
    }
}
```

**Core Function**

```cpp
void update(double timestep) {
    buildRhs();
    project(600, timestep);
    applyPressure(timestep);

    _d->advect(timestep, *_u, *_v);
    _u->advect(timestep, *_u, *_v);
    _v->advect(timestep, *_u, *_v);

    /* Make effect of advection visible, since it's not an in-place operation */
    _d->flip();
    _u->flip();
    _v->flip();
}
```

**BuildRhs Function** : Using the obvious central differences (take a look at the MAC grid again), we approximate the two-dimensional divergence in fluid grid cell $(i,j)$ as

$$
(\nabla \cdot \vec u)_{i,j} \approx \frac{u_{i+1/2,j} - u_{i-1/2,j}}{\Delta x} + \frac{v_{i,j+1/2} - v_{i,j-1/2}}{\Delta x}
$$

$$
(\nabla \cdot \vec u)_{i,j} \approx \frac{u_{i+1/2,j} - u_{i-1/2,j}}{\Delta x} + \frac{v_{i,j+1/2} - v_{i,j-1/2}}{\Delta x}
$$

```cpp
/* Builds the pressure right hand side as the negative divergence */
void buildRhs() {
    double scale = 1.0/_hx;

    for (int y = 0, idx = 0; y < _h; y++) {
        for (int x = 0; x < _w; x++, idx++) {
            _r[idx] = -scale*(_u->at(x + 1, y) - _u->at(x, y) +
                                _v->at(x, y + 1) - _v->at(x, y));
        }
    }
}
```

**Project Function**:

The Pressure Equation:

$$
\frac{\Delta t}{\rho}(\frac{4p_{i,j} - p_{i+1,j} - p_{i,j+1} - p_{i-1,j} - p_{i,j-1}}{\Delta x^2}) = -(\frac{u_{i+1/2,j} - u_{i-1/2,j}}{\Delta x} + \frac{v_{i,j+1/2} - v_{i,j-1/2}}{\Delta x})
$$

numerical approximations to the Poisson problem $-\Delta t/\rho\nabla\cdot\nabla p = -\nabla \cdot \vec u$.
Here we adopt Gause-Seidel to solve this equation.

```cpp
/* Performs the pressure solve using Gauss-Seidel.
 * The solver will run as long as it takes to get the relative error below
 * a threshold, but will never exceed `limit' iterations
 */
void project(int limit, double timestep) {
    double scale = timestep/(_density*_hx*_hx);

    double maxDelta;
    for (int iter = 0; iter < limit; iter++) {
        maxDelta = 0.0;
        for (int y = 0, idx = 0; y < _h; y++) {
            for (int x = 0; x < _w; x++, idx++) {
                int idx = x + y*_w;

                double diag = 0.0, offDiag = 0.0;

                /* Here we build the matrix implicitly as the five-point
                    * stencil. Grid borders are assumed to be solid, i.e.
                    * there is no fluid outside the simulation domain.
                    */
                if (x > 0) {
                    diag    += scale;
                    offDiag -= scale*_p[idx - 1];
                }
                if (y > 0) {
                    diag    += scale;
                    offDiag -= scale*_p[idx - _w];
                }
                if (x < _w - 1) {
                    diag    += scale;
                    offDiag -= scale*_p[idx + 1];
                }
                if (y < _h - 1) {
                    diag    += scale;
                    offDiag -= scale*_p[idx + _w];
                }

                double newP = (_r[idx] - offDiag)/diag;

                maxDelta = max(maxDelta, fabs(_p[idx] - newP));

                _p[idx] = newP;
            }
        }

        if (maxDelta < 1e-5) {
            printf("Exiting solver after %d iterations, maximum change is %f\n", iter, maxDelta);
            return;
        }
    }

    printf("Exceeded budget of %d iterations, maximum change was %f\n", limit, maxDelta);
}
```

**The Discrete Pressure Gradient**
Callback for the MAC grid: the staggering makes accurate central differences robust.

The formulas for the pressure update in two dimensions, using the central difference approximations for $\partial p / \partial x$ and $\partial p / \partial y$ :

$$
\begin{align}
u_{i+1/2,j}^n+1 = u_{i+1/2,j} - \Delta t \frac{1}{\rho}\frac{p_{i+1, j} - p_{i, j}}{\Delta x} \\
v_{i,j+1/2}^n+1 = v_{i,j+1/2} - \Delta t \frac{1}{\rho}\frac{p_{i, j+1} - p_{i, j}}{\Delta x} \\
\end{align}
$$

```cpp
/* Applies the computed pressure to the velocity field */
void applyPressure(double timestep) {
    double scale = timestep/(_density*_hx);

    for (int y = 0, idx = 0; y < _h; y++) {
        for (int x = 0; x < _w; x++, idx++) {
            _u->at(x,     y    ) -= scale*_p[idx];
            _u->at(x + 1, y    ) += scale*_p[idx];
            _v->at(x,     y    ) -= scale*_p[idx];
            _v->at(x,     y + 1) += scale*_p[idx];
        }
    }

    for (int y = 0; y < _h; y++)
        _u->at(0, y) = _u->at(_w, y) = 0.0;
    for (int x = 0; x < _w; x++)
        _v->at(x, 0) = _v->at(x, _h) = 0.0;
}
```

**Advect**

```cpp
_d->advect(timestep, *_u, *_v);
_u->advect(timestep, *_u, *_v);
_v->advect(timestep, *_u, *_v);

/* Make effect of advection visible, since it's not an in-place operation */
_d->flip();
_u->flip();
_v->flip();
```

### FluidSimulation

- Initial Setting:

  - Fluid Domain: 128 \* 128
  - Fluid Density: 0.1 → $\rho$
  - Solver TimeStep: 0.005 → $\Delta t$

  ```cpp
  const int sizeX = 128;
  const int sizeY = 128

  const double density = 0.1
  ```

- Define Solver
  ```cpp
  FluidSolver *solver = nwe FluidSolver(sizeX, sizeY, density);
  ```
- Iteration
  ```cpp
  while(time < 8.0) {
  	/* Use four substeps per iteration */
    for (int i = 0; i < 4; i++) {
        solver->addInflow(0.45, 0.2, 0.1, 0.01, 1.0, 0.0, 3.0);
        solver->update(timestep);
        time += timestep;
        fflush(stdout);
     }
  }
  ```

## Render Result

<div align="center">
<img src="./result/0001-0400.gif" width="30%"/>
</div>
