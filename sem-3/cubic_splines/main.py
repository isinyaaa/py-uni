from typing import Iterator, Union
import math
import matplotlib.pyplot as plt
import numpy as np


def problem1_gen(num: int) -> Iterator[int]:
    for aux in range(num):
        yield 10 * (2 ** aux)


def gen_points(n: int, s: int = 0) -> list[float]:
    return np.array([i/n for i in range(s, n + 1)])


def test_function(x: float) -> float:
    if x < 0 or x > 1:
        raise ValueError("x must be between 0 and 1")
    return 1/(2 - x)


def problem1_vals(n: int):
    x = gen_points(n)

    # function values
    y = np.array([
        test_function(j)
        for j in x
    ])
    # h_i := x_i - x_{i-1}
    h = np.array([0] + [
        xi - xp
        for xi, xp in zip(x[1:], x[:-1])
    ])

    # RHS values
    d = np.zeros(n)
    for i in range(1, n - 1):
        di = 6/(h[i] + h[i + 1]) * \
            ((y[i + 1] - y[i])/h[i + 1] -
             (y[i] - y[i - 1])/h[i])
        d[i] = di

    return x, y, h, d


def tridiagonal_matrix_gen(n: int, h: list[float], f_val: int = 0):
    # 'lower' diagonal
    a = np.zeros(n)
    # diagonal
    b = np.zeros(n)
    # 'upper' diagonal
    c = np.zeros(n)
    # spare values in a and c shouldn't be a problem
    for i in range(1, n - 1):
        a[i - f_val] += h[i]/(h[i] + h[i + 1])
        b[i - f_val] += 2
        c[i - f_val] += h[i + 1]/(h[i] + h[i + 1])

    return a, b, c


def gauss_method(n: int,
                 a: list[float],
                 b: list[float],
                 c: list[float],
                 d: list[float],
                 f_val: int = 0) -> (list[float], list[float]):

    # resulting main diagonal
    u = np.zeros(n - f_val)
    u[0] = b[0]
    # resulting RHS
    y = np.zeros(n - f_val)
    y[0] = d[0]
    for i in range(1, n - f_val):
        mult = a[i]/u[i - 1]
        u[i] = b[i] - mult * c[i - 1]
        y[i] = d[i] - mult * y[i - 1]

    return u, y


def plot1(n: int,
          tl: list[float],
          y: list[float],
          label: str,
          samples: int = 100):

    t_vals = []
    for t in tl[::-1]:
        t_vals.extend(t)

    plt.plot(t_vals, [y(t) for t in t_vals])

    #plt.plot(t_vals, [test_function(t) for t in t_vals])

    plt.xlabel('x')
    plt.ylabel('y')

    plt.title(label)

    plt.show()

def test_pp(x: float) -> float:
    if x < 0 or x > 1:
        raise ValueError("x must be between 0 and 1")
    return -(2 - 2 * x)/(2 - x) ** 4


def estimate_err(x_vals: list[float], S_vals: list[float]):
    return max([math.fabs(test_function(i/1000) - S_vals(i/1000)) for i in range(1001)])


def problem1(num: int, _plot: bool, samples: int = 100):
    if num == -1:
        num = 5

    for n in problem1_gen(num):
        print("n =", n)
        x_vals, y_vals, h_vals, d_vals = problem1_vals(n)

        a_vals, b_vals, c_vals = tridiagonal_matrix_gen(n, h_vals)

        # natural spline: m[0] = m[-1] = 0
        b_vals[0] = b_vals[-1] = 1

        u, y = gauss_method(n, a_vals, b_vals, c_vals, d_vals)

        # then the solution
        m = np.zeros(n)
        m[-1] = y[-1]/u[-1]
        for i in range(n - 2, -1, -1):
            m[i] = (y[i] - c_vals[i] * m[i + 1])/u[i]

        S_vals = get_S(n, x_vals, y_vals, h_vals, m)
        t_vals = []
        for i, (x, xp) in enumerate(zip(x_vals[1:], x_vals[:-1])):
            t_vals.append(np.linspace(x, xp, samples))

        if _plot:
            plot1(n, t_vals, S_vals, 'natural')

        print("erro máximo para o spline natural:", estimate_err(x_vals, S_vals))

        # 'complete' spline:

        c_vals[0] = 0
        b_vals[0] = 1

        b_vals[-1] = 1
        a_vals[-1] = 0

        d_vals[0] = test_pp(x_vals[0])
        d_vals[-1] = test_pp(x_vals[-1])

        u, y = gauss_method(n, a_vals, b_vals, c_vals, d_vals)

        # then the solution
        m = np.zeros(n)
        m[-1] = y[-1]/u[-1]
        for i in range(n - 2, -1, -1):
            m[i] = (y[i] - c_vals[i] * m[i + 1])/u[i]

        S_vals = get_S(n, x_vals, y_vals, h_vals, m)
        t_vals = []
        for i, (x, xp) in enumerate(zip(x_vals[1:], x_vals[:-1])):
            t_vals.append(np.linspace(x, xp, samples))

        if _plot:
            plot1(n, t_vals, S_vals, 'completo')

        print("erro máximo para o spline completo:", estimate_err(x_vals, S_vals))

        # not a knot

        c_vals[1] = 0
        b_vals[1] = 2 + h_vals[1]/h_vals[2]
        a_vals[1] = 1 - h_vals[1]/h_vals[2]

        c_vals[-2] = 1 - h_vals[-2]/h_vals[-1]
        b_vals[-2] = 2 + h_vals[-2]/h_vals[-1]
        a_vals[-2] = 0

        u, y = gauss_method(n - 1, a_vals, b_vals, c_vals, d_vals, 1)

        # then the solution
        m = np.zeros(n)
        m[-2] = y[-2]/u[-2]
        for i in range(n - 3, 0, -1):
            m[i] = (y[i] - c_vals[i] * m[i + 1])/u[i]

        m[0] = (h_vals[1] + h_vals[2])/h_vals[2] * m[0] -\
            h_vals[1]/h_vals[2] * m[2]
        m[-1] = (h_vals[-2] + h_vals[-1])/h_vals[-2] * m[-1] -\
            h_vals[-1]/h_vals[-2] * m[-2]

        S_vals = get_S(n, x_vals, y_vals, h_vals, m)
        t_vals = []
        for i, (x, xp) in enumerate(zip(x_vals[1:], x_vals[:-1])):
            t_vals.append(np.linspace(x, xp, samples))

        if _plot:
            plot1(n, t_vals, S_vals, 'not a knot')

        print("erro máximo para o spline not a knot:", estimate_err(x_vals, S_vals))


def solve_lsystem(n: int,
                  a: list[float],
                  b: list[float],
                  c: list[float],
                  d: list[float]) -> (list[float], list[float]):
    u, y = gauss_method(n, a, b, c, d)

    z = np.zeros(n)
    z[-1] = y[-1]/u[-1]
    for i in range(n - 2, -1, -1):
        z[i] = (y[i] - c[i] * z[i + 1])/u[i]

    return z


def plot2(n: int,
          t: list[float],
          x: list[float],
          y: list[float],
          samples: int = 200):

    t_vals = []
    for tp, tn in zip(t[:-1], t[1:]):
        t_vals.extend(np.linspace(tp, tn, samples, endpoint=False))

    plt.plot([x(t) for t in t_vals], [y(t) for t in t_vals])

    plt.xlabel('x')
    plt.ylabel('y')

    plt.title('periodico')

    plt.show()


def get_periodic(n: int,
                 f: list[float],
                 h: list[float]):

    a_vals, b_vals, c_vals = tridiagonal_matrix_gen(n, h)

    # fill first and last values
    a_vals[1] = h[1]/(h[1] + h[2])
    b_vals[1] = 2
    c_vals[1] = h[2]/(h[1] + h[2])

    a_vals[-1] = h[-1]/(h[1] + h[-1])
    b_vals[-1] = 2
    c_vals[-1] = h[1]/(h[1] + h[-1])

    # RHS
    d_vals = np.zeros(n)
    for i in range(1, n - 1):
        di = 6/(h[i] + h[i + 1]) * \
            ((f[i + 1] - f[i])/h[i + 1] -
             (f[i] - f[i - 1])/h[i])
        d_vals[i] = di

    d_vals[0] = 6/(h[1] + h[2]) * \
        ((f[2] - f[1])/h[2] -
         (f[1] - f[0])/h[1])

    d_vals[-1] = 6/(h[-1] + h[1]) * \
        ((f[1] - f[-1])/h[1] -
         (f[-1] - f[-2])/h[-1])

    # solution for Ãz~ = d~
    z_t = solve_lsystem(n - 2, a_vals[1:-1], b_vals[1:-1], c_vals[1:-1], d_vals[1:-1])

    u_t = np.array([a_vals[1]] + (n - 4) * [0] + [c_vals[-2]])

    y_t = solve_lsystem(n - 2, a_vals[1:-1], b_vals[1:-1], c_vals[1:-1], u_t)

    #v_t = np.array([c_vals[-1]] + (n - 3) * [0] + [a_vals[-1]])
    xn = (d_vals[-1] - (c_vals[-1] * z_t[0] + a_vals[-1] * z_t[-1]))/(b_vals[-1] - (c_vals[-1] * y_t[0] + a_vals[-1] * y_t[-1]))
    xt_vals = np.subtract(z_t, xn * y_t)

    x_res = np.concatenate(([xn], xt_vals, [xn]))

    return x_res


def get_S(n: int,
          x: list[float],
          y: list[float],
          h: list[float],
          m: list[float],
          is_periodic: bool = False):

    adder = 1 if is_periodic else 0
    def f(xv: float) -> float:
        i = 1
        while xv > x[i]:
            i += 1
        return (x[i] - xv)/h[i] * y[i - 1] + (xv - x[i - 1])/h[i] * y[i] +\
        h[i] ** 2 / 6 *\
        (
            (((x[i] - xv)/h[i]) ** 3 - (x[i] - xv)/h[i]) * m[i - 2 + adder] +
            (((xv - x[i - 1])/h[i]) ** 3 - (xv - x[i - 1])/h[i]) * m[i - 1 + adder]
        )

    return f


def problem2(_plot: bool):
    n = 12

    # x(t) vals
    x_vals = np.array([
        25, 19, 13, 9, 5, 2.2, 1, 3, 8, 13, 18, 25
    ])

    # y(t) vals
    y_vals = np.array([
        5, 7.5, 9.1, 9.4, 9, 7.5, 5, 2.1, 2, 3.5, 4.5, 5
    ])

    # t parameter, equivalent to x_vals
    t_vals = np.zeros(n)
    for i in range(1, n):
        t_vals[i] = t_vals[i - 1] +\
                    math.sqrt((x_vals[i] - x_vals[i - 1]) ** 2 +
                              (y_vals[i] - y_vals[i - 1]) ** 2)

    # h values
    h_vals = np.array([0] + [
        xi - xp
        for xi, xp in zip(t_vals[1:], t_vals[:-1])
    ])

    # -- x(t) solution --#

    x_res = get_periodic(n, x_vals, h_vals)

    Sx_vals = get_S(n, t_vals, x_vals, h_vals, x_res, True)

    # -- y(t) solution -- #

    y_res = get_periodic(n, y_vals, h_vals)

    Sy_vals = get_S(n, t_vals, y_vals, h_vals, y_res, True)

    if _plot:
        plot2(n, t_vals, Sx_vals, Sy_vals)


def main(args: list[Union[int, str]]):
    _plot = False
    n = -1

    match len(args):
        case 2:
            _plot = True
            n = int(args[0])
        case 1:
            n = int(args[0])

    problem1(n, _plot)

    problem2(_plot)


if __name__ == '__main__':
    from sys import argv
    main(argv[1:])
