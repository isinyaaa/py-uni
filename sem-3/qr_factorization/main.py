import os
import math
import numpy as np

__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))


def modified_gram_schmidt(A):
    m, n = A.shape
    Q = np.zeros((m, n))
    R = np.zeros((n, n))
    A = np.transpose(A)
    R[0][0] = np.linalg.norm(A[0], 2)
    Q[:, 0] = A[0]/R[0][0]

    for i in range(1, n):
        Q[:, i] = A[i]
        for j in range(0, i):
            R[j][i] = np.matmul(Q[:, j], Q[:, i])
            Q[:, i] -= (R[j][i] * Q[:, j])
        norm = np.linalg.norm(Q[:, i], 2)
        if norm == 0:
            raise Exception("Matrix is not LI")
        Q[:, i] /= norm
        R[i][i] = norm

    return Q, R


def solve_lsystem(A, b):
    m, n = A.shape
    if m < n:
        raise Exception("System is not overdetermined")
    Q, R = modified_gram_schmidt(A)
    z = np.zeros(n, dtype=float)
    rhs = np.copy(b)
    for i in range(n):
        z[i] = np.inner(Q[:, i], rhs)
        rhs -= z[i] * Q[:, i]

    x = np.zeros(n)
    x[-1] = z[-1] / R[-1, -1]

    for i in range(n - 2, -1, -1):
        acc = sum([x[j] * R[i, j] for j in range(i, n)])
        x[i] = (z[i] - acc) / R[i, i]

    return x


def task1():
    with open(os.path.join(__location__, "Exemplo1.txt")) as f:
        lines = f.readlines()
        m, n = lines[0].split(' ')
        m, n = int(m), int(n)
        A = np.ndarray((m, n))
        b = np.zeros(m)
        for i in range(m):
            *A[i, :], b[i] = lines[i + 1].split(' ')

    print("(A b) =")
    print(np.array([list(A[i, :]) + [b[i]] for i in range(A.shape[0])]))
    sol = solve_lsystem(A, b)
    dist = np.linalg.norm(np.matmul(A, sol) - b)
    print(f"\nsolution:\n{sol}\nnorm of the residue:\n{dist}\n")


def task2():
    with open(os.path.join(__location__, "Exemplo2.txt")) as f:
        lines = f.readlines()
        m = int(lines[0])
        time = np.zeros(m)
        b = np.zeros(m)
        for i in range(m):
            time[i], b[i] = lines[i + 1].split(' ')

    A = np.array([
        [
            ((time[i] - 1950) / 50) ** j
            for j in range(4)
        ] for i in range(11)
    ])

    print("(A b) =")
    print(np.array([list(A[i, :]) + [b[i]] for i in range(A.shape[0])]))
    sol = solve_lsystem(A, b)
    dist = np.linalg.norm(np.matmul(A, sol) - b)
    print(f"\nsolution:\n{sol}\nnorm of the residue:\n{dist}\n")

    yys = []
    for t in time:
        s = (t - 1950) / 50
        val = sol[0] + sol[1] * s + sol[2] * (s ** 2) + sol[3] * (s ** 3)
        yys.append(val)

    print("(ti, yi, yyi)")
    for t, y, yy in zip(time, b, yys):
        print(f"({t}, {y}, {yy})")

    next_val = sol[0]
    for i in range(1, 4):
        next_val += sol[i] * ((6 / 5) ** i)

    print(f"\nApproximation for the population in 2010:\n{next_val}")


def task3():
    with open(os.path.join(__location__, "Exemplo2.txt")) as f:
        lines = f.readlines()
        m = int(lines[0])
        xc = np.zeros(m)
        yc = np.zeros(m)
        for i in range(m):
            xc[i], yc[i] = lines[i + 1].split(' ')

    A = np.array([
        [
            xc[i] ** 2,
            xc[i] * yc[i],
            yc[i] ** 2,
            xc[i],
            yc[i]
        ] for i in range(10)
    ])

    b = np.array([-1] * 10, dtype=float)

    print("(A b) =")
    print(np.array([list(A[i, :]) + [b[i]] for i in range(A.shape[0])]))
    sol = solve_lsystem(A, b)
    dist = np.linalg.norm(np.matmul(A, sol) - b)
    print(f"\nsolution:\n{sol}\nnorm of the residue:\n{dist}\n")

    yys = []
    for x, y in zip(xc, yc):
        c1 = sol[2]
        c2 = x * sol[1] + sol[4]
        c3 = sol[0] * (x ** 2) + sol[3] * x + 1
        delta = math.sqrt(c2 ** 2 - 4 * c1 * c3)
        y0 = (-c2 - delta) / (2 * c1)
        y1 = (-c2 + delta) / (2 * c1)
        yys.append(y0 if abs(y - y0) < abs(y - y1) else y1)

    print("(xi, yi, yyi)")
    for x, y, yy in zip(xc, b, yys):
        print(f"({x}, {y}, {yy})")


def task4():
    with open(os.path.join(__location__, "Exemplo4.txt")) as f:
        lines = f.readlines()
        m = int(lines[0])
        xc = np.zeros(m)
        b = np.zeros(m)
        for i in range(m):
            xc[i], b[i] = lines[i + 1].split(' ')

    A = np.array([
        [
            xc[i] ** 3,
            xc[i] ** 2,
            xc[i],
        ] for i in range(10)
    ])

    print("(A b) =")
    print(np.array([list(A[i, :]) + [b[i]] for i in range(A.shape[0])]))
    sol = solve_lsystem(A, b)
    dist = np.linalg.norm(np.matmul(A, sol) - b)
    print(f"\nsolution:\n{sol}\nnorm of the residue:\n{dist}\n")

    next_val = sol[0] * 110 ** 3 + sol[1] * 110 ** 2 + sol[2] * 110

    print(f"Approximation for the function at 110:\n{next_val}")

    yys = [sol[0] * x ** 3 + sol[1] * x ** 2 + sol[2] * x for x in xc]

    print()

    print("(xi, yi, yyi)")
    for x, y, yy in zip(xc, b, yys):
        print(f"({x}, {y}, {yy})")


def main():
    print(f"{' Task 1 ':=^79}")
    task1()
    print()
    print(f"{' Task 2 ':=^79}")
    task2()
    print()
    print(f"{' Task 3 ':=^79}")
    task3()
    print()
    print(f"{' Task 4 ':=^79}")
    task4()


if __name__ == "__main__":
    main()
