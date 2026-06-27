import os
import numpy as np
import matplotlib.pyplot as plt


def exact_solution(x):
    """
    Analytic solution:

        u(x) = x sin(3 pi x)
    """
    return x * np.sin(3 * np.pi * x)


def right_hand_side(x):
    """
    Right-hand side of the Poisson problem: -u'' = f

    For u(x) = x sin(3 pi x), we have f(x) = -6 pi cos(3 pi x) + 9 pi^2 x sin(3 pi x)
    """
    return -6 * np.pi * np.cos(3 * np.pi * x) + 9 * np.pi**2 * x * np.sin(3 * np.pi * x)


def build_finite_difference_matrix(N):
    """
    Build the finite difference matrix for

        -u'' = f

    on the interval (0, 1) with homogeneous Dirichlet boundary conditions.

    Parameters:
        N: number such that h = 1 / N

    There are N+1 grid points including the boundary points:
        x_0, x_1, ..., x_N

    The unknowns are the interior values:
        u_1, ..., u_{N-1}

    Therefore, the matrix has size (N-1) x (N-1).
    """
    h = 1 / N
    size = N - 1

    main_diagonal = 2 * np.ones(size)
    off_diagonal = -1 * np.ones(size - 1)

    A = (
        np.diag(main_diagonal)
        + np.diag(off_diagonal, k=1)
        + np.diag(off_diagonal, k=-1)
    )

    A = A / h**2

    return A


def solve_bvp(N):
    h = 1 / N

    # Interior grid points: x_1, ..., x_{N-1}
    x_inner = np.linspace(h, 1 - h, N - 1)

    A = build_finite_difference_matrix(N)
    b = right_hand_side(x_inner)

    numerical_solution = np.linalg.solve(A, b)
    exact_values = exact_solution(x_inner)

    error = compute_discrete_l2_error(exact_values, numerical_solution)

    return x_inner, numerical_solution, exact_values, error


def compute_discrete_l2_error(exact_values, numerical_solution):
    """
    Compute the discrete L2 error:

        ||u - u_h||_2 =
        sqrt( 1/(N-1) * sum_i (u(x_i) - u_i)^2 )
    """
    difference = exact_values - numerical_solution
    error = np.sqrt(np.mean(difference**2))

    return error


def estimate_convergence_order(error_coarse, error_fine):
    
    alpha = np.log(error_coarse / error_fine) / np.log(2)

    return alpha


def run_convergence_experiment():
  
    N_values = [8, 16, 32, 64, 128, 256]

    errors = []
    orders = []

    print("Poisson Boundary Value Problem")
    print("------------------------------")
    print("Problem: -u'' = f in (0,1), u(0)=u(1)=0")
    print("Exact solution: u(x) = x sin(3 pi x)")
    print()

    print(f"{'N':>8} {'h':>12} {'error':>18} {'order':>12}")
    print("-" * 55)

    previous_error = None

    for N in N_values:
        h = 1 / N

        x_inner, numerical_solution, exact_values, error = solve_bvp(N)

        errors.append(error)

        if previous_error is None:
            order = None
            orders.append(None)
            order_text = "-"
        else:
            order = estimate_convergence_order(previous_error, error)
            orders.append(order)
            order_text = f"{order:.6f}"

        print(f"{N:8d} {h:12.6f} {error:18.10e} {order_text:>12}")

        previous_error = error

    save_results_table(N_values, errors, orders)
    plot_convergence(N_values, errors)
    plot_solution_comparison(N_values[-1])

    return N_values, errors, orders


def save_results_table(N_values, errors, orders):
    """
    Save the convergence results as a CSV file.
    """
    os.makedirs("results", exist_ok=True)

    file_path = "results/bvp_convergence_results.csv"

    with open(file_path, "w") as file:
        file.write("N,h,error,estimated_order\n")

        for N, error, order in zip(N_values, errors, orders):
            h = 1 / N

            if order is None:
                order_text = ""
            else:
                order_text = str(order)

            file.write(f"{N},{h},{error},{order_text}\n")

    print()
    print(f"Results saved to {file_path}")


def plot_convergence(N_values, errors):
    """
    Plot the convergence behavior on a log-log scale.
    """
    os.makedirs("plots", exist_ok=True)

    h_values = [1 / N for N in N_values]

    plt.figure(figsize=(7, 5))
    plt.loglog(h_values, errors, marker="o", label="Numerical error")

    # Reference line for second-order convergence
    reference = errors[-1] * (np.array(h_values) / h_values[-1]) ** 2
    plt.loglog(h_values, reference, linestyle="--", label="Reference slope 2")

    plt.xlabel("Mesh width h")
    plt.ylabel("Discrete L2 error")
    plt.title("Convergence of Finite Difference Approximation")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig("plots/bvp_convergence_order.png")
    plt.show()

    print("Convergence plot saved to plots/bvp_convergence_order.png")


def plot_solution_comparison(N):
    """
    Plot the numerical and analytic solutions on the finest grid.
    """
    os.makedirs("plots", exist_ok=True)

    x_inner, numerical_solution, exact_values, error = solve_bvp(N)

    plt.figure(figsize=(7, 5))
    plt.plot(x_inner, exact_values, label="Analytic solution")
    plt.plot(x_inner, numerical_solution, marker="o", linestyle="", label="Numerical solution")

    plt.xlabel("x")
    plt.ylabel("u(x)")
    plt.title(f"Analytic vs Numerical Solution, N={N}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig("plots/bvp_solution_comparison.png")
    plt.show()

    print("Solution comparison plot saved to plots/bvp_solution_comparison.png")


if __name__ == "__main__":
    run_convergence_experiment()
