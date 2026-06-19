import math
import matplotlib.pyplot as plt


def exact_solution(x):
    """
    Analytic solution of the stiff IVP:

        y'(x) = 2000(cos(x) - y(x)),   y(0) = 0

    The exact solution is:

        y(x) =
        (4000000 / 4000001) * cos(x)
        + (2000 / 4000001) * sin(x)
        - (4000000 / 4000001) * exp(-2000x)
    """
    return (
        (4000000 / 4000001) * math.cos(x)
        + (2000 / 4000001) * math.sin(x)
        - (4000000 / 4000001) * math.exp(-2000 * x)
    )


def implicit_euler(N):
    """
    Implicit Euler method for

        y'(x) = 2000(cos(x) - y(x))

    Formula:

        y_{n+1} = (y_n + 2000h cos(x_{n+1})) / (1 + 2000h)
    """
    a = 0.0
    b = 1.5
    h = (b - a) / N

    xs = [a]
    ys = [0.0]

    y = 0.0

    for n in range(N):
        x_next = a + (n + 1) * h

        y = (y + 2000 * h * math.cos(x_next)) / (1 + 2000 * h)

        xs.append(x_next)
        ys.append(y)

    return xs, ys


def trapezoidal_rule(N):
    """
    Trapezoidal rule for

        y'(x) = 2000(cos(x) - y(x))

    Formula:

        y_{n+1}
        =
        ((1 - 1000h)y_n + 1000h(cos(x_n) + cos(x_{n+1})))
        / (1 + 1000h)
    """
    a = 0.0
    b = 1.5
    h = (b - a) / N

    xs = [a]
    ys = [0.0]

    y = 0.0

    for n in range(N):
        x_current = a + n * h
        x_next = a + (n + 1) * h

        y = (
            (1 - 1000 * h) * y
            + 1000 * h * (math.cos(x_current) + math.cos(x_next))
        ) / (1 + 1000 * h)

        xs.append(x_next)
        ys.append(y)

    return xs, ys


def exact_values(xs):
    """
    Compute exact solution values for a list of x-values.
    """
    return [exact_solution(x) for x in xs]


def plot_results():
    """
    Plot the exact solution together with the numerical solutions
    from Implicit Euler and the Trapezoidal Rule.
    """
    N = 40

    xs_implicit, ys_implicit = implicit_euler(N)
    xs_trap, ys_trap = trapezoidal_rule(N)

    # Use many points for a smooth exact solution curve
    smooth_xs = [1.5 * i / 1000 for i in range(1001)]
    smooth_exact = exact_values(smooth_xs)

    plt.figure(figsize=(9, 5))

    plt.plot(smooth_xs, smooth_exact, label="Exact Solution")
    plt.plot(xs_implicit, ys_implicit, marker="o", label="Implicit Euler")
    plt.plot(xs_trap, ys_trap, marker="s", label="Trapezoidal Rule")

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("A-stability vs. L-stability for a Stiff IVP")
    plt.legend()
    plt.grid(True)

    plt.savefig("stiff_ivp_a_vs_l_stability.png", dpi=300)
    plt.show()


def print_results():
    """
    Print numerical values for comparison.
    """
    N = 40
    h = 1.5 / N

    xs_implicit, ys_implicit = implicit_euler(N)
    xs_trap, ys_trap = trapezoidal_rule(N)

    print("=" * 80)
    print("Stiff IVP:")
    print("y'(x) = 2000(cos(x) - y(x)),   y(0) = 0,   x in [0, 1.5]")
    print("=" * 80)
    print(f"Number of intervals: N = {N}")
    print(f"Step size: h = {h}")
    print("=" * 80)
    print()

    print(f"{'x':>10} {'Exact':>20} {'Implicit Euler':>20} {'Trapezoidal':>20}")
    print("-" * 80)

    for i in range(len(xs_implicit)):
        x = xs_implicit[i]
        y_exact = exact_solution(x)
        y_imp = ys_implicit[i]
        y_trap = ys_trap[i]

        print(f"{x:10.4f} {y_exact:20.10e} {y_imp:20.10e} {y_trap:20.10e}")


def main():
    print_results()
    plot_results()


if __name__ == "__main__":
    main()
