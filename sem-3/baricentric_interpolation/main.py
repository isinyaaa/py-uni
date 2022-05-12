from typing import Callable
import math
import numpy as np


def chebyshev_baricentric_weights(n: int) -> list[float]:
    return np.array([0.5] + [
            (-1) ** j
            for j in range(1, n - 1)
        ] + [0.5 * (-1) ** n])


def chebyshev_indexed_points(a: float, b: float, n: int) -> list[float]:
    return np.array([
        - (b - a) / 2 * math.cos(j * math.pi/n) + (a + b) / 2
        for j in range(0, n)
    ])


def lagrange_dp_j(
        x: list[float], w: list[float],
        j: float, n: int) -> Callable[[int], float]:
    return \
       lambda i: - 2 * (w[j] / w[i])/(x[i] - x[j]) \
       * np.sum([
                    (np.delete(w, i)[k] / w[i]) / (x[i] - np.delete(x, i)[k])
                    - 1 / (x[i] - x[j])
                    for k in range(n - 1)
                ])


def interpolate(
        a: float, b: float,
        alpha: float, beta: float,
        q: Callable[[float], float], f: Callable[[float], float],
        n: int) -> Callable[[float], float]:
    # RHS matrix
    A = np.ndarray((n - 2, n - 2))
    lhs = np.ndarray(n - 2)

    x = chebyshev_indexed_points(a, b, n)
    w = chebyshev_baricentric_weights(n)

    # lagrange terms for current row
    cur_row_ldp = [None] * n

    # those terms won't change throughout
    cur_row_ldp[0] = lagrange_dp_j(x, w, 0, n)
    cur_row_ldp[n - 1] = lagrange_dp_j(x, w, n - 1, n)

    # fillup array core with n - 2 functions
    cur_row_ldp[1: n - 1] = [
        lagrange_dp_j(x, w, j, n)
        for j in range(1, n - 1)
    ]

    # function to fillup LHS
    def lhs_const_func(i):
        return f(x[i]) + \
            cur_row_ldp[0](i) * alpha + \
            cur_row_ldp[n - 1](i) * beta

    for i in range(1, n - 1):
        A[i - 1, : n - 2] = np.array([
                cur_row_ldp[j](i) if i != j else 0
                for j in range(1, n - 1)
            ])

        A[i - 1, i - 1] = - np.sum([
                                cur_row_ldp[j](i) if i != j else 0
                                for j in range(0, n)
                            ]) + q(x[i])
        lhs[i - 1] = lhs_const_func(i)

    solution = np.append(np.append([alpha], np.linalg.solve(A, lhs)), [beta])

    def y_approx(x_val):
        res = 0
        for j in range(n):
            prod = 1
            for k, xk in enumerate(x):
                if k == j:
                    continue
                prod *= (x_val - xk) / (x[j] - xk)

            res += solution[j] * prod

        return res

    return y_approx


def test(
        a: float, b: float,
        alpha: float, beta: float,
        q: Callable[[float], float], f: Callable[[float], float],
        u: Callable[[float], float]) -> list[float]:
    """
    Parameters:
        a (float): x on start of interval
        b (float): x on end of interval
        alpha (float): y on start of interval
        beta (float): y on end of interval
        q (Callable): linear functional term on ODE
        f (Callable): independent functional term on ODE
        u (Callable): exact functional solution for ODE
    """

    print(("                (   a  , alpha ) ~ (   b  ,  beta )\n"
          f"using endpoints ({a:.4f}, {alpha:.4f}) ~ ({b:.4f}, {beta:.4f})\n"))

    test_points = np.array([
                a + (j - 0.5) * (b - a) / 1000
                for j in range(1, 1001)
            ])

    for n in list([
                16 * (2 ** i)
                for i in range(4)
            ]):
        print(f"for n = {n:3}:")
        y = interpolate(a, b, alpha, beta, q, f, n)

        error_list = np.array([
                math.fabs(u(t) - y(t))
                for t in test_points
            ])
        largest_error_index = np.argmax(error_list)
        error_x = test_points[largest_error_index]
        y_string = f"{y(error_x):.5f}" \
            if y(error_x) > 0 else f"({y(error_x):.5f})"
        print((f"\tmax error @\t  x = {error_x:.3f}\n"
               "\t"
               f"|{u(error_x):.5f} - {y_string}| = "
               f"{error_list[largest_error_index]:.5f}"))


def main():
    print(f"{' Example function #1 ':=^79}")
    test(-5, 5, 1/26, 1/26,
         lambda x: 6 * (x ** 2) / (1 + x ** 2) ** 2,
         lambda x: 2 / (1 + x ** 2) ** 3,
         lambda x: 1 / (1 + x ** 2))

    print()
    print(f"{' Example function #2 ':=^79}")
    test(0, 1, 0, 1,
         lambda x: 10000,
         lambda x: 0,
         lambda x: math.sinh(100 * x) / math.sinh(100))


if __name__ == '__main__':
    main()
