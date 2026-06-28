import numpy as np


def right_hand_side(x):
    """
    Right-hand side:

        f(x) = -6 pi cos(3 pi x) + 9 pi^2 x sin(3 pi x)
    """
    return -6 * np.pi * np.cos(3 * np.pi * x) + 9 * np.pi**2 * x * np.sin(3 * np.pi * x)


def build_matrix(N):
    """
    Build the finite difference matrix for

        -u'' = f

    with u(0)=0 and u(1)=0.

    N is the number of subintervals, so h = 1/N.
    The unknowns are the interior points x_1, ..., x_{N-1}.
    """
    h = 1 / N
    size = N - 1

    A = np.zeros((size, size))

    for i in range(size):
        A[i, i] = 2 / h**2

        if i > 0:
            A[i, i - 1] = -1 / h**2

        if i < size - 1:
            A[i, i + 1] = -1 / h**2

    return A


def build_rhs(N):
    """
    Build the right-hand side vector.
    """
    h = 1 / N
    x = np.linspace(h, 1 - h, N - 1)

    return right_hand_side(x)


def residual_norm(A, x, b):
    """
    Compute ||Ax - b||_2.
    """
    return np.linalg.norm(A @ x - b)


def damped_jacobi(A, b, omega, tol=1e-10, max_iter=100000):
    """
    Damped Jacobi method:

        x_new = x_old + omega * D^{-1}(b - A x_old)
    """
    n = len(b)
    x = np.zeros(n)

    D = np.diag(A)

    for k in range(1, max_iter + 1):
        r = b - A @ x
        x = x + omega * r / D

        if residual_norm(A, x, b) < tol:
            return k, True

    return max_iter, False


def sor(A, b, omega, tol=1e-10, max_iter=10000):
    """
    SOR method.
    """
    n = len(b)
    x = np.zeros(n)

    for k in range(1, max_iter + 1):
        for i in range(n):
            left_sum = 0
            right_sum = 0

            for j in range(i):
                left_sum += A[i, j] * x[j]

            for j in range(i + 1, n):
                right_sum += A[i, j] * x[j]

            gauss_seidel_value = (b[i] - left_sum - right_sum) / A[i, i]
            x[i] = (1 - omega) * x[i] + omega * gauss_seidel_value

        if residual_norm(A, x, b) < tol:
            return k, True

    return max_iter, False


def run_experiment():
    """
    Run a small experiment for several grid sizes.
    """
    N_values = [8, 16, 32, 64]

    jacobi_omega_values = [0.5, 0.8, 1.0]
    sor_omega_values = [1.0, 1.2, 1.5, 1.8]

    print("Jacobi and SOR Experiment")
    print("=========================")
    print()

    for N in N_values:
        A = build_matrix(N)
        b = build_rhs(N)

        condition_number = np.linalg.cond(A)

        print(f"N = {N}")
        print(f"h = {1 / N}")
        print(f"Condition number: {condition_number:.4e}")
        print()

        print("Damped Jacobi:")
        for omega in jacobi_omega_values:
            iterations, converged = damped_jacobi(A, b, omega)

            print(
                f"  omega = {omega:.1f}, "
                f"iterations = {iterations}, "
                f"converged = {converged}"
            )

        print()

        print("SOR:")
        best_omega = None
        best_iterations = None

        for omega in sor_omega_values:
            iterations, converged = sor(A, b, omega)

            print(
                f"  omega = {omega:.1f}, "
                f"iterations = {iterations}, "
                f"converged = {converged}"
            )

            if converged:
                if best_iterations is None or iterations < best_iterations:
                    best_iterations = iterations
                    best_omega = omega

        print()

        if best_omega is not None:
            print(f"Best SOR omega for N={N}: {best_omega}")
            print(f"Best SOR iterations: {best_iterations}")
        else:
            print(f"No SOR parameter converged for N={N}")

        print()
        print("-" * 50)
        print()


if __name__ == "__main__":
    run_experiment()
