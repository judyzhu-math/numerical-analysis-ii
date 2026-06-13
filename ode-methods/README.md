# ODE Stability: Euler and Trapezoidal Methods

This folder contains a Python implementation comparing three numerical methods for solving the model ordinary differential equation

$$
y' = \lambda y, \quad y(0)=1, \quad x \in [0,1].
$$

In this experiment, we use

$$
\lambda = -10.
$$

The exact solution is

$$
y(x)=e^{\lambda x}.
$$

Therefore, at \(x=1\),

$$
y(1)=e^{-10}.
$$

## Methods

The following methods are implemented:

1. Explicit Euler Method
2. Implicit Euler Method
3. Trapezoidal Rule

For different values of \(N\), the interval \([0,1]\) is divided into \(N\) equal steps, so the step size is

$$
h = \frac{1}{N}.
$$

The program compares the numerical approximation of \(y(1)\) with the exact value and prints the absolute error.

## Numerical Formulas

### Explicit Euler Method

$$
y_{n+1} = (1+h\lambda)y_n
$$

This method is simple, but it may become unstable when the step size is too large.

### Implicit Euler Method

$$
y_{n+1} = \frac{y_n}{1-h\lambda}
$$

This method is more stable for rapidly decaying problems.

### Trapezoidal Rule

$$
y_{n+1}
=
\frac{1+\frac{h\lambda}{2}}
{1-\frac{h\lambda}{2}}
y_n
$$

The trapezoidal rule usually has better accuracy than the Euler methods.

## How to Run

In the terminal, run:

```bash
python euler_trapezoidal_comparison.py
