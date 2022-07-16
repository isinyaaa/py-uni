import numpy as np

from constants import AVOGADRO_CONSTANT


def stochastic_simulation_algorithm(X0, st_mat, p_functions, t0=0, tf=50,
        vol=None, event=None, record=False):
    N = len(X0)

    if record:
        ts = [t0]
        Xs = [X0]
        if vol != None:
            c = AVOGADRO_CONSTANT * vol
            ys = [[X0[i]/c for i in range(N)]]

    check_intervention = event != None
    t = t0
    X = X0
    while t < tf:
        if check_intervention:
            check_intervention = event(t, X)

        a = [p(X) for p in p_functions]

        asum = np.sum(a)
        # should those two random numbers be the same?
        rand = np.random.random_sample()
        j = np.amin(np.argwhere(rand < np.cumsum(a / asum)))  # lesser element larger than rand
        tau = np.log(1 / rand) / asum
        X += st_mat[j, :]
        t += tau
        if record:
            ts.append(t)
            Xs.append(X)
            if vol != None:
                ys.append([X[i]/c for i in range(N)])

    if record:
        if vol != None:
            return ts, Xs, ys
        return ts, Xs


def euler_maruyama_algorithm(X0, st_mat, p_functions, tf=50, tstep=250,
        vol=None, record=False):
    tau = tf / tstep

    N = len(X0)

    if record:
        ts = [n * tau for n in range(tstep)]
        Xs = [X0]
        if vol != None:
            c = AVOGADRO_CONSTANT * vol
            ys = [[X0[i]/c for i in range(N)]]

    X = X0
    for _ in range(tstep):
        a = [p(X) for p in p_functions]

        rand = np.random.random_sample()
        X += np.inner(tau * a + np.sqrt(np.absolute(tau * a)) * rand,
                np.transpose(st_mat))

        if record:
            Xs.append(X)
            if vol != None:
                ys.append([X[i]/c for i in range(N)])

    if record:
        if vol != None:
            return ts, Xs, ys
        return ts, Xs
