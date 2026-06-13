import math


def exact_solution(x, lam=-10):
    """
    Exact solution of the model problem

        y' = lambda * y,   y(0) = 1

    The exact solution is

        y(x) = exp(lambda * x)
    """
    return math.exp(lam * x)


def explicit_euler(lam, N):
    """
    Explicit Euler method.

    Formula:
        y_{n+1} = (1 + h * lambda) * y_n
    """
    h = 1.0 / N
    y = 1.0

    for _ in range(N):
        y = (1 + h * lam) * y

    return y


def implicit_euler(lam, N):
    """
    Implicit Euler method.

    Formula:
        y_{n+1} = y_n / (1 - h * lambda)
    """
    h = 1.0 / N
    y = 1.0

    for _ in range(N):
        y = y / (1 - h * lam)

    return y


def trapezoidal_rule(lam, N):
    """
    Trapezoidal rule.

    Formula:
        y_{n+1} = ((1 + h * lambda / 2) / (1 - h * lambda / 2)) * y_n
    """
    h = 1.0 / N
    y = 1.0

    for _ in range(N):
        y = y * (1 + h * lam / 2) / (1 - h * lam / 2)

    return y


def compute_error(numerical_value, exact_value):
    """
    Absolute error between numerical and exact solution.
    """
    return abs(numerical_value - exact_value)


def print_results(lam, Ns):
    """
    Print a comparison table for different numerical methods.
    """
    y_exact = exact_solution(1.0, lam)

    methods = [
        ("Explicit Euler", explicit_euler),
        ("Implicit Euler", implicit_euler),
        ("Trapezoidal Rule", trapezoidal_rule),
    ]

    print("=" * 70)
    print(f"Model problem: y' = {lam}y,  y(0) = 1,  interval [0, 1]")
    print(f"Exact solution at x = 1: y(1) = exp({lam}) = {y_exact:.10e}")
    print("=" * 70)
    print()

    for name, method in methods:
        print("-" * 70)
        print(f"Method: {name}")
        print("-" * 70)
        print(f"{'N':>6} {'h':>12} {'Numerical y(1)':>20} {'Absolute Error':>20}")

        for N in Ns:
            h = 1.0 / N
            y_num = method(lam, N)
            error = compute_error(y_num, y_exact)

            print(f"{N:6d} {h:12.6f} {y_num:20.10e} {error:20.10e}")

        print()


def main():
    lam = -10
    Ns = [2, 4, 8, 16, 32, 64]

    print_results(lam, Ns)


if __name__ == "__main__":
    main()
